"""Spectrally matched, uncertainty-aware difference imaging of two real SPHEREx epochs.

Two SPHEREx images of one field sample different wavelengths at most sky pixels
(linear-variable filters, different pointings). This module compares only the
pixels whose per-pixel CWAVE values agree within ``spectral_tolerance_fraction``
of the local CBAND, using each observation's own DataLink calibration.

Processing, all on real archive products:
  1. native background plane per epoch on sigma-clipped blank sky (catalog
     sources and SOURCE-flagged pixels excluded), with fit covariance;
  2. optional cross-convolution with the measured ePSFs projected onto common
     sky axes (A * PSF_B and B * PSF_A share the same effective PSF);
  3. bilinear reprojection of B onto A with squared-weight variance;
  4. per-pixel spectral gate; empirical photometric scale from catalog stars in
     the matched region (applied only when inconsistent with unity);
  5. D = B - s A with propagated variance plus an explicit spectral-mismatch term;
  6. validation: blank-sky null test (noise renormalisation, inflate-only),
     field-star fractional residual floor, catalog-star registration residual,
     and injection-recovery of synthetic point sources on the real noise field;
  7. connected residual components classified against the propagated catalog.

Candidates are residuals that need follow-up, never discoveries. Gates and
policy numbers are analysis choices recorded with every result.
"""
import base64
import json
from pathlib import Path

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time
from scipy import ndimage

from spherex_lab.imaging.psf_processing import extract_prepared_epsf, fit_native_background, project_kernel
from spherex_lab.imaging.source_detection import detect_sources
from spherex_lab.imaging.wcs_geometry import reproject_bilinear

SOURCE_BIT = 1 << 21
TARGET_ID = 'TARGET'
MAS_RAD = np.deg2rad(1/3600000)
POLICY = {
    'spectral_tolerance_fraction': 0.10,   # |lambda_B - lambda_A| <= 0.10 x CBAND_A per pixel
    'spectral_slope_bound': 4.0,           # |dlnF/dlnlambda| bound for the mismatch systematic term
    'background_clip_sigma': 3.0,
    'background_max_iterations': 10,
    'catalog_mask_radius_px': 2.5,
    'aperture_half_width_px': 2,           # 5x5 native-pixel box
    'star_min_snr': 20.0,
    'registration_tolerance_px': 0.3,
    'registration_min_stars': 3,
    'detection_sigma': 5.0,
    'detection_min_pixels': 2,
    'association_radius_px': 1.5,
    'edge_margin_px': 2,
    'min_matched_pixels': 200,
    'min_null_pixels': 200,
    'psf_min_kernel_coverage': 0.95,       # masked convolution needs >=95% of |kernel| weight valid
    'max_noise_scale': 3.0,
    'vet_isolation_radius_px': 4,          # no flagged/invalid pixels this close to a variability candidate
    'vet_blend_radius_px': 2.5,            # no other catalog source inside the aperture
    'vet_min_shape_r2': 0.5,
    'min_noise_model_stars': 30,           # stars needed before the source-level noise model is trusted               # residual must be explained by a scaled copy of the source
    'injection_count': 40,
    'injection_snr': 10.0,
    'injection_seed': 20261115,
}


# ----------------------------------------------------------------------------- inputs

def load_epoch(path, publisher_id, calibration_cache_dir, *, product_identity=None, cache=None):
    """Read one real cutout, its ePSF table, and its own DataLink wavelength maps."""
    from spherex_lab.data.calibration_resolver import DEFAULT_CACHE, wavelength_maps
    from spherex_lab.imaging.native_photometry import prepare_observation
    observation = prepare_observation(path)
    wave, band, calibration_valid, provenance = wavelength_maps(
        path, publisher_id, calibration_cache_dir, product_identity=product_identity, cache=cache or DEFAULT_CACHE)
    with fits.open(path, memmap=False) as hdus:
        zodi = np.array(hdus['ZODI'].data, dtype=float) if 'ZODI' in hdus else None
    header = observation['science']['header']
    return {'id': publisher_id, 'path': str(path), 'science': observation['science'], 'wcs': observation['wcs'],
            'epsf': observation['epsf'], 'wave': np.asarray(wave, float), 'band': np.asarray(band, float),
            'calibration_valid': np.asarray(calibration_valid, bool), 'calibration': provenance, 'zodi': zodi,
            'date': str(header['DATE-OBS']), 'mjd': float(Time(header['DATE-OBS'], scale='utc').mjd),
            'detector': int(header['DETECTOR']), 'psf_fwhm_arcsec': header.get('PSF_FWHM')}


def propagate_catalog(sources, date):
    """Vectorised rectilinear PM propagation; identical formula to reference_catalog.source_position."""
    rows = [s for s in sources if s.get('ra') is not None and s.get('dec') is not None]
    ra = np.deg2rad([float(s['ra']) for s in rows]); dec = np.deg2rad([float(s['dec']) for s in rows])
    pmra = np.array([float(s['pmra']) if s.get('pmra') not in (None, '') else 0. for s in rows])
    pmdec = np.array([float(s['pmdec']) if s.get('pmdec') not in (None, '') else 0. for s in rows])
    epoch = np.array([float(s['ref_epoch']) if s.get('ref_epoch') not in (None, '') else np.nan for s in rows])
    dt = np.where(np.isfinite(epoch), float(Time(date).jyear) - epoch, 0.)
    pmra = np.where(np.isfinite(pmra), pmra, 0.); pmdec = np.where(np.isfinite(pmdec), pmdec, 0.)
    east = np.stack([-np.sin(ra), np.cos(ra), np.zeros_like(ra)], axis=-1)
    north = np.stack([-np.sin(dec)*np.cos(ra), -np.sin(dec)*np.sin(ra), np.cos(dec)], axis=-1)
    direction = np.stack([np.cos(dec)*np.cos(ra), np.cos(dec)*np.sin(ra), np.sin(dec)], axis=-1)
    point = direction + (dt*MAS_RAD)[:, None]*(pmra[:, None]*east + pmdec[:, None]*north)
    point /= np.linalg.norm(point, axis=1, keepdims=True)
    return {'ids': [str(s.get('source_id')) for s in rows],
            'g_mag': np.array([float(s['phot_g_mean_mag']) if s.get('phot_g_mean_mag') not in (None, '') else np.nan for s in rows]),
            'ra': np.rad2deg(np.arctan2(point[:, 1], point[:, 0])) % 360,
            'dec': np.rad2deg(np.arcsin(np.clip(point[:, 2], -1, 1)))}


def _catalog_on_grid(catalog, wcs, shape):
    if catalog is None or not len(catalog['ids']):
        return None
    x, y = wcs.world_to_pixel_values(catalog['ra'], catalog['dec'])
    inside = np.isfinite(x) & np.isfinite(y) & (x > -5) & (y > -5) & (x < shape[1]+4) & (y < shape[0]+4)
    return {'ids': [catalog['ids'][i] for i in np.flatnonzero(inside)], 'g_mag': catalog['g_mag'][inside],
            'x': np.asarray(x)[inside], 'y': np.asarray(y)[inside]}


def _disk_mask(shape, x, y, radius):
    mask = np.zeros(shape, bool)
    if x is None or not len(x):
        return mask
    yy, xx = np.indices(shape)
    r = int(np.ceil(radius))
    for cx, cy in zip(x, y):
        x0, x1 = max(int(np.floor(cx))-r, 0), min(int(np.ceil(cx))+r+1, shape[1])
        y0, y1 = max(int(np.floor(cy))-r, 0), min(int(np.ceil(cy))+r+1, shape[0])
        if x0 >= x1 or y0 >= y1:
            continue
        sub = (xx[y0:y1, x0:x1]-cx)**2 + (yy[y0:y1, x0:x1]-cy)**2 <= radius**2
        mask[y0:y1, x0:x1] |= sub
    return mask


# ----------------------------------------------------------------------------- steps

def fit_background(image, variance, blank, policy=POLICY):
    """Iteratively clipped weighted plane on blank native pixels (wraps fit_native_background)."""
    mask = np.asarray(blank, bool) & np.isfinite(image) & np.isfinite(variance) & (variance > 0)
    scale = np.nan
    for iteration in range(1, policy['background_max_iterations']+1):
        fit = fit_native_background(image, variance, mask)
        z = fit['residual']/np.sqrt(np.where(variance > 0, variance, np.nan))
        values = z[mask]
        center = float(np.median(values)); scale = float(1.4826*np.median(np.abs(values-center)))
        keep = mask & np.isfinite(z) & (np.abs(z-center) <= policy['background_clip_sigma']*max(scale, 1e-12))
        if np.array_equal(keep, mask):
            break
        mask = keep
    fit = fit_native_background(image, variance, mask)
    return {**fit, 'iterations': iteration, 'blank_pixel_count': int(mask.sum()),
            'robust_residual_scale': scale}


def _aperture(image, variance, valid, x, y, half):
    cx, cy = int(np.rint(x)), int(np.rint(y))
    if cx-half < 0 or cy-half < 0 or cx+half >= image.shape[1] or cy+half >= image.shape[0]:
        return None
    window = (slice(cy-half, cy+half+1), slice(cx-half, cx+half+1))
    if not valid[window].all():
        return None
    return float(image[window].sum()), float(variance[window].sum())


def _centroid(image, valid, x, y, half):
    cx, cy = int(np.rint(x)), int(np.rint(y))
    if cx-half < 0 or cy-half < 0 or cx+half >= image.shape[1] or cy+half >= image.shape[0]:
        return None
    window = (slice(cy-half, cy+half+1), slice(cx-half, cx+half+1))
    if not valid[window].all():
        return None
    w = np.clip(image[window], 0, None)
    if w.sum() <= 0:
        return None
    yy, xx = np.mgrid[cy-half:cy+half+1, cx-half:cx+half+1]
    return float((w*xx).sum()/w.sum()), float((w*yy).sum()/w.sum())


def convolve_masked(image, variance, valid, kernel, min_coverage):
    """Normalised convolution tolerant of isolated invalid pixels.

    out = sum(k a m)/sum(k m), var = sum(k^2 v m)/sum(k m)^2, defined where the
    valid pixels carry at least ``min_coverage`` of sum|k|. Exact-support
    convolution (psf_processing.convolve_independent) discards every output
    pixel touching a flagged input, which removes most of a real cutout. The
    renormalisation biases point-source profiles slightly where pixels are
    missing; that effect is part of the empirical field-star floor.
    """
    a = np.asarray(image, float); v = np.asarray(variance, float); k = np.asarray(kernel, float)
    m = np.asarray(valid, bool) & np.isfinite(a) & np.isfinite(v) & (v >= 0)
    k = k/k.sum(); am = np.where(m, a, 0.); vm = np.where(m, v, 0.); mf = m.astype(float)
    numerator = ndimage.convolve(am, k, mode='constant'); weight = ndimage.convolve(mf, k, mode='constant')
    coverage = ndimage.convolve(mf, np.abs(k), mode='constant')/np.abs(k).sum()
    var = ndimage.convolve(vm, k**2, mode='constant')
    good = (coverage >= min_coverage) & (weight > 0.5)
    out = np.divide(numerator, weight, out=np.full(a.shape, np.nan), where=good)
    out_var = np.divide(var, weight**2, out=np.full(a.shape, np.nan), where=good)
    return {'image': out, 'variance': out_var, 'valid': good, 'coverage_fraction': float(good.mean())}


def excess_scatter(values, sigmas):
    """sigma_sys such that values/sqrt(sigmas^2+sigma_sys^2) has robust unit scatter.

    The empirical systematic floor of differential photometry: the scatter of
    constant-source residuals beyond their own propagated statistical errors.
    Zero when the statistical errors already explain the scatter.
    """
    v = np.asarray(values, float); e = np.asarray(sigmas, float)
    good = np.isfinite(v) & np.isfinite(e) & (e > 0)
    v, e = v[good], e[good]
    if v.size < 5:
        return None
    center = np.median(v)
    def spread(extra):
        z = (v-center)/np.sqrt(e**2+extra**2)
        return 1.4826*np.median(np.abs(z-np.median(z)))
    if spread(0.) <= 1:
        return 0.
    lo, hi = 0., float(np.max(np.abs(v-center))+np.max(e))
    for _ in range(80):
        mid = (lo+hi)/2
        lo, hi = (mid, hi) if spread(mid) > 1 else (lo, mid)
    return float(hi)


def _robust(values):
    v = np.asarray(values, float); v = v[np.isfinite(v)]
    if not v.size:
        return None, None
    median = float(np.median(v))
    return median, float(1.4826*np.median(np.abs(v-median)))


def _psf_center(a, b, target_coord):
    shape = a['science']['image'].shape
    if target_coord is not None:
        xa, ya = a['wcs'].world_to_pixel(target_coord); xb, yb = b['wcs'].world_to_pixel(target_coord)
        sb = b['science']['image'].shape
        if 0 <= xa <= shape[1]-1 and 0 <= ya <= shape[0]-1 and 0 <= xb <= sb[1]-1 and 0 <= yb <= sb[0]-1:
            return target_coord
    return a['wcs'].pixel_to_world((shape[1]-1)/2, (shape[0]-1)/2)


def _measure_stars(A, Avar, Bimg, Bvar, valid, cat_a, cat_b, half, min_snr):
    stars = []
    if cat_a is None:
        return stars
    index_b = {ident: k for k, ident in enumerate(cat_b['ids'])}
    for k, ident in enumerate(cat_a['ids']):
        j = index_b.get(ident)
        if j is None or ident == TARGET_ID:
            continue
        pa = _aperture(A, Avar, valid, cat_a['x'][k], cat_a['y'][k], half)
        pb = _aperture(Bimg, Bvar, valid, cat_b['x'][j], cat_b['y'][j], half)
        if pa is None or pb is None or pa[0] <= 0 or pa[1] <= 0:
            continue
        snr = pa[0]/np.sqrt(pa[1])
        if snr < min_snr:
            continue
        stars.append({'source_id': ident, 'g_mag': float(cat_a['g_mag'][k]), 'x_a': float(cat_a['x'][k]), 'y_a': float(cat_a['y'][k]),
                      'x_b': float(cat_b['x'][j]), 'y_b': float(cat_b['y'][j]), 'flux_a': pa[0], 'var_a': pa[1],
                      'flux_b': pb[0], 'var_b': pb[1], 'snr_a': float(snr)})
    return stars


def _inject(D, variance, interior, kernel, positions, snr, half):
    """Add PSF-shaped point sources so their aperture S/N equals ``snr``; returns fluxes."""
    out = D.copy(); fluxes = []
    k = kernel/kernel.sum(); ry, rx = np.array(k.shape)//2
    for x, y in positions:
        cx, cy = int(x), int(y)
        window = (slice(cy-half, cy+half+1), slice(cx-half, cx+half+1))
        noise = float(np.sqrt(np.nansum(variance[window])))
        frac = float(k[ry-half:ry+half+1, rx-half:rx+half+1].sum()) if min(ry, rx) >= half else 1.
        flux = snr*noise/frac
        out[cy-ry:cy+ry+1, cx-rx:cx+rx+1] += flux*k
        fluxes.append(flux)
    return out, fluxes


def _vet_variable(star, flagged, valid, dist_to_edge, A, D, catalog, policy):
    """Independent checks that a source-level flux change is not a known systematic."""
    half = policy['aperture_half_width_px']; x, y = star['x_a'], star['y_a']
    cx, cy = int(np.rint(x)), int(np.rint(y)); r = int(policy['vet_isolation_radius_px'])
    window = (slice(max(cy-r, 0), cy+r+1), slice(max(cx-r, 0), cx+r+1))
    yy, xx = np.mgrid[window]
    near = (xx-x)**2+(yy-y)**2 <= r**2
    checks = {}
    checks['isolated_from_flagged_pixels'] = bool(not flagged[window][near].any())
    checks['away_from_region_edge'] = bool(dist_to_edge[cy, cx] > half+policy['edge_margin_px'])
    others = np.hypot(catalog['x']-x, catalog['y']-y) if catalog is not None else np.array([])
    checks['unblended_in_catalog'] = bool(np.sum(others <= policy['vet_blend_radius_px']) <= 1)
    stamp = (slice(cy-half, cy+half+1), slice(cx-half, cx+half+1))
    a, d = A[stamp].ravel(), D[stamp].ravel()
    good = np.isfinite(a) & np.isfinite(d)
    r2 = None
    if good.sum() >= 9 and np.dot(a[good], a[good]) > 0 and np.dot(d[good], d[good]) > 0:
        c = np.dot(a[good], d[good])/np.dot(a[good], a[good])
        r2 = float(1-np.sum((d[good]-c*a[good])**2)/np.sum(d[good]**2))
    checks['psf_shaped_residual_r2'] = r2
    failed = [k for k, v in checks.items() if v is False] + (['psf_shaped_residual'] if r2 is None or r2 < policy['vet_min_shape_r2'] else [])
    return {'checks': checks, 'failed': failed}


def _injection_recovery(D, variance, interior, blank, kernel, policy):
    rng = np.random.default_rng(policy['injection_seed'])
    half = policy['aperture_half_width_px']; ry, rx = np.array(kernel.shape)//2
    pad = max(ry, rx, half)+1
    # Injection centres avoid catalog sources; apertures must lie in the interior.
    eligible = ndimage.binary_erosion(interior, structure=np.ones((3, 3)), iterations=half) & blank
    border = np.zeros_like(eligible); border[pad:-pad, pad:-pad] = True
    eligible &= border
    candidates = np.argwhere(eligible)
    if not len(candidates):
        return {'status': 'UNAVAILABLE', 'reason': 'no blank interior region large enough for injections'}
    rng.shuffle(candidates)
    chosen = []
    for y, x in candidates:
        if all((x-cx)**2+(y-cy)**2 > (2*pad+2)**2 for cx, cy in chosen):
            chosen.append((int(x), int(y)))
        if len(chosen) >= policy['injection_count']:
            break
    injected, fluxes = _inject(D, variance, interior, kernel, chosen, policy['injection_snr'], half)
    found = detect_sources(injected, variance, interior, sigma_threshold=policy['detection_sigma'],
                           min_pixels=policy['detection_min_pixels'], correlated=True)['sources']
    recovered = []
    for (x, y), flux in zip(chosen, fluxes):
        near = [s for s in found if np.hypot(s['x']-x, s['y']-y) <= policy['association_radius_px']]
        if near:
            best = min(near, key=lambda s: np.hypot(s['x']-x, s['y']-y))
            measured = _aperture(injected-D, variance, np.ones_like(interior), x, y, half)
            recovered.append({'x': x, 'y': y, 'offset_px': float(np.hypot(best['x']-x, best['y']-y)),
                              'flux_ratio': float(measured[0]/flux) if measured else None})
    return {'status': 'COMPLETED', 'seed': policy['injection_seed'], 'injected': len(chosen),
            'recovered': len(recovered), 'completeness': float(len(recovered)/len(chosen)) if chosen else None,
            'aperture_snr': policy['injection_snr'],
            'median_position_offset_px': float(np.median([r['offset_px'] for r in recovered])) if recovered else None,
            'method': 'PSF-shaped sources added to the real difference image at blank interior positions; detection rerun with the production threshold.'}


def _classify(component, sign, stars_by_id, cat_a, cat_b, candidates_other_sign, floor, noise_scale, policy, catalog_available, dist_to_edge):
    x, y = component['x'], component['y']
    cx, cy = int(np.rint(x)), int(np.rint(y))
    result = {'x': x, 'y': y, 'sign': 'brightening' if sign > 0 else 'fading', 'pixel_count': component['pixel_count'],
              'peak_sigma': float(sign*component['peak_marginal_sigma'])}
    if 0 <= cy < dist_to_edge.shape[0] and 0 <= cx < dist_to_edge.shape[1] and dist_to_edge[cy, cx] <= 1.5:
        result.update(classification='EDGE_ARTIFACT_LIKELY', reason='Component touches the boundary of the matched/valid region.')
        return result
    # Association first: a known mover produces a dipole by construction.
    best = None
    for grid, epoch in ((cat_a, 'A'), (cat_b, 'B')):
        if grid is None or not len(grid['x']):
            continue
        d = np.hypot(grid['x']-x, grid['y']-y); k = int(np.argmin(d))
        if d[k] <= policy['association_radius_px'] and (best is None or d[k] < best[0]):
            best = (float(d[k]), grid['ids'][k], float(grid['g_mag'][k]), epoch)
    if best is not None:
        separation, ident, g, epoch = best
        result.update(catalog_source_id=ident, catalog_g_mag=g, catalog_separation_px=separation)
        ia = cat_a['ids'].index(ident) if ident in cat_a['ids'] else None
        ib = cat_b['ids'].index(ident) if ident in cat_b['ids'] else None
        motion = float(np.hypot(cat_b['x'][ib]-cat_a['x'][ia], cat_b['y'][ib]-cat_a['y'][ia])) if ia is not None and ib is not None else 0.
        result['catalog_motion_px'] = motion
        if motion >= 0.3:
            who = 'The requested target' if ident == TARGET_ID else 'Catalog proper motion of this source'
            result.update(classification='KNOWN_PROPER_MOTION_RESIDUAL', reason=f'{who} moves {motion:.2f} px between epochs; the residual is expected motion, not variability.')
            return result
    if any(np.hypot(o['x']-x, o['y']-y) <= 2.5 for o in candidates_other_sign):
        result.update(classification='DIPOLE_PSF_OR_REGISTRATION', reason='Adjacent opposite-sign residual: PSF or sub-pixel registration mismatch, not a flux change.')
        return result
    if np.sign(component['aperture_sum']) != sign:
        result.update(classification='DIPOLE_PSF_OR_REGISTRATION', reason='Net 5x5 aperture residual has the opposite sign to the peak: a shape mismatch, not a point-source flux change.')
        return result
    if best is None and not catalog_available:
        result.update(classification='UNCLASSIFIED_NO_CATALOG', reason='No Gaia catalog available for association; only the requested target is checked.')
        return result
    if best is None:
        result.update(classification='UNASSOCIATED_CANDIDATE',
                      reason='No Gaia DR3 source (G<=19) within association radius: possible transient, solar-system object, source below the catalog limit, or artifact. Follow-up required.')
        return result
    if ident == TARGET_ID:
        result.update(classification='TARGET_RESIDUAL', reason='Residual at the requested target; see the target measurement for its flux change.')
        return result
    star = stars_by_id.get(ident)
    if star is None:
        result.update(classification='CATALOG_SOURCE_UNMEASURED', reason='Associated star lacks a complete valid aperture in both epochs.')
        return result
    change = star['fractional_change']; z = star['change_z']
    result.update(fractional_change=float(change), fractional_change_sigma=float(star['total_sigma']/star['scaled_flux_a']), change_z=float(z))
    if abs(z) >= policy['detection_sigma']:
        result.update(classification='CANDIDATE_VARIABLE', reason='Catalog source change exceeds the detection threshold after the empirical field-star floor.')
    else:
        result.update(classification='CONSISTENT_WITH_SYSTEMATIC_FLOOR', reason='Residual is within the empirical field-star systematic floor.')
    return result


def compare(a, b, *, catalog_sources=None, target=None, policy=POLICY, psf_match=True):
    """Compare epoch B with epoch A on A's grid. Returns arrays plus a JSON-safe report."""
    if a['detector'] != b['detector']:
        raise ValueError('epochs are on different detectors; no response transformation is established')
    if b['mjd'] < a['mjd']:
        a, b = b, a
    sa, sb = a['science'], b['science']; shape = sa['image'].shape
    half = policy['aperture_half_width_px']
    catalog_available = catalog_sources is not None
    sources = list(catalog_sources) if catalog_available else []
    if target is not None and target.get('pm_ra_cosdec_mas_yr') is not None and target.get('epoch') == 'J2000':
        sources.append({'source_id': TARGET_ID, 'ra': target['ra_deg'], 'dec': target['dec_deg'], 'pmra': target['pm_ra_cosdec_mas_yr'],
                        'pmdec': target['pm_dec_mas_yr'], 'ref_epoch': 2000.0, 'phot_g_mean_mag': None})
    cat_epoch_a = propagate_catalog(sources, a['date']) if sources else None
    cat_epoch_b = propagate_catalog(sources, b['date']) if sources else None
    cat_a = _catalog_on_grid(cat_epoch_a, a['wcs'], shape)                       # epoch-A positions, A grid
    cat_b = _catalog_on_grid(cat_epoch_b, a['wcs'], shape)                       # epoch-B positions, A grid
    cat_bn = _catalog_on_grid(cat_epoch_b, b['wcs'], sb['image'].shape)          # epoch-B positions, B grid
    radius = policy['catalog_mask_radius_px']

    def blank(epoch, grid):
        s = epoch['science']
        mask = s['valid'] & (s['flags'] == 0) & epoch['calibration_valid']
        if grid is not None:
            mask &= ~_disk_mask(s['image'].shape, grid['x'], grid['y'], radius)
        return mask

    bg_a = fit_background(sa['image'], sa['variance'], blank(a, cat_a), policy)
    bg_b = fit_background(sb['image'], sb['variance'], blank(b, cat_bn), policy)
    a_valid = sa['valid'] & np.isfinite(bg_a['residual'])
    b_valid = sb['valid'] & np.isfinite(bg_b['residual'])
    target_coord = None; target_info = None
    if target is not None:
        target_coord = SkyCoord(float(target['ra_deg']), float(target['dec_deg']), unit='deg')
    center = _psf_center(a, b, target_coord)
    xa, ya = a['wcs'].world_to_pixel(center); xb, yb = b['wcs'].world_to_pixel(center)
    psf_a = extract_prepared_epsf(a['epsf'], float(xa), float(ya)); psf_b = extract_prepared_epsf(b['epsf'], float(xb), float(yb))
    if psf_match:
        kb_on_a = project_kernel(psf_b['kernel'], b['wcs'], a['wcs'], center)['kernel']
        ka_on_b = project_kernel(psf_a['kernel'], a['wcs'], b['wcs'], center)['kernel']
        A = convolve_masked(bg_a['residual'], bg_a['residual_variance'], a_valid, kb_on_a, policy['psf_min_kernel_coverage'])
        Bn = convolve_masked(bg_b['residual'], bg_b['residual_variance'], b_valid, ka_on_b, policy['psf_min_kernel_coverage'])
        effective_kernel = ndimage.convolve(np.pad(psf_a['kernel'], [(s//2, s//2) for s in kb_on_a.shape]), kb_on_a, mode='constant')
    else:
        A = {'image': np.where(a_valid, bg_a['residual'], np.nan), 'variance': np.where(a_valid, bg_a['residual_variance'], np.nan), 'valid': a_valid}
        Bn = {'image': np.where(b_valid, bg_b['residual'], np.nan), 'variance': np.where(b_valid, bg_b['residual_variance'], np.nan), 'valid': b_valid}
        effective_kernel = psf_a['kernel']
    aligned = reproject_bilinear(np.nan_to_num(Bn['image']), np.nan_to_num(Bn['variance']), Bn['valid'], b['wcs'], a['wcs'], shape)
    # Bilinear resampling smooths B by f(1-f) px^2 per axis (f = fractional source offset).
    # Give A the same second moment with a symmetric, unshifted 3-tap kernel.
    yy, xx = np.indices(shape)
    sx, sy = b['wcs'].world_to_pixel_values(*a['wcs'].pixel_to_world_values(xx, yy))
    fx, fy = sx-np.floor(sx), sy-np.floor(sy)
    vx, vy = float(np.nanmean(fx*(1-fx))), float(np.nanmean(fy*(1-fy)))
    smoothing = np.outer([vy/2, 1-vy, vy/2], [vx/2, 1-vx, vx/2])
    A = convolve_masked(A['image'], A['variance'], A['valid'], smoothing, policy['psf_min_kernel_coverage'])
    wave_b = reproject_bilinear(np.nan_to_num(b['wave']), np.zeros_like(b['wave']), b['calibration_valid'], b['wcs'], a['wcs'], shape)
    delta = wave_b['image'] - a['wave']
    overlap = A['valid'] & aligned['valid'] & a['calibration_valid'] & wave_b['valid']
    spectral = overlap & np.isfinite(delta) & (np.abs(delta) <= policy['spectral_tolerance_fraction']*a['band'])
    valid = spectral & np.isfinite(A['image']) & np.isfinite(aligned['image'])
    interior = ndimage.binary_erosion(valid, structure=np.ones((3, 3)), iterations=policy['edge_margin_px'])
    dist_to_edge = ndimage.distance_transform_edt(valid)
    # Region boundary only (interior holes from flagged pixels are handled by the isolation check).
    dist_to_boundary = ndimage.distance_transform_edt(ndimage.binary_fill_holes(valid))
    Aimg = np.where(valid, A['image'], np.nan); Avar = np.where(valid, A['variance'], np.nan)
    Bimg = np.where(valid, aligned['image'], np.nan); Bvar = np.where(valid, aligned['variance'], np.nan)

    report = {'status': None, 'observation_a': a['id'], 'observation_b': b['id'], 'policy': policy,
              'epochs': {key: {'observation_id': e['id'], 'date': e['date'], 'mjd': e['mjd'], 'detector': e['detector'],
                               'calibration': e['calibration'], 'psf_fwhm_arcsec': e['psf_fwhm_arcsec'],
                               'zodi_median_mjy_sr': float(np.nanmedian(e['zodi'])) if e['zodi'] is not None else None,
                               'science_valid_fraction': float(e['science']['valid'].mean())}
                         for key, e in (('a', a), ('b', b))},
              'separation_days': b['mjd']-a['mjd']}
    rel = np.abs(delta)/a['band']
    report['spectral'] = {'overlap_pixels': int(overlap.sum()), 'matched_pixels': int(valid.sum()), 'interior_pixels': int(interior.sum()),
                          'matched_fraction_of_overlap': float(valid.sum()/max(overlap.sum(), 1)),
                          'median_delta_over_band_matched': float(np.median(rel[valid])) if valid.any() else None,
                          'median_delta_over_band_overlap': float(np.nanmedian(rel[overlap])) if overlap.any() else None,
                          'wavelength_range_um': [float(np.nanmin(a['wave'][valid])), float(np.nanmax(a['wave'][valid]))] if valid.any() else None}
    report['background'] = {k: {'coefficients_mjy_sr': fit['coefficients'].tolist(), 'blank_pixels': fit['blank_pixel_count'],
                                'iterations': fit['iterations'], 'robust_residual_scale': fit['robust_residual_scale']}
                            for k, fit in (('a', bg_a), ('b', bg_b))}
    report['psf'] = {'interpolation_smoothing_variance_px2': [vx, vy], 'matched': bool(psf_match), 'method': 'masked normalised cross-convolution with projected measured ePSFs' if psf_match else 'none',
                     'kernel_coverage_policy': policy['psf_min_kernel_coverage'], 'zone_a': psf_a['provenance']['zone_index'], 'zone_b': psf_b['provenance']['zone_index'],
                     'kernel_shape': list(effective_kernel.shape), 'center_deg': [float(center.ra.deg), float(center.dec.deg)]}
    gates = {}
    if valid.sum() < policy['min_matched_pixels']:
        report.update(status='REJECTED', reason=f"Only {int(valid.sum())} spectrally matched valid pixels (< {policy['min_matched_pixels']}).",
                      gates={'spectral_overlap': 'FAIL'})
        return _result(report, a, Aimg, Bimg, np.full(shape, np.nan), np.full(shape, np.nan), np.full(shape, np.nan), valid, interior, np.zeros(shape, bool), rel)
    gates['spectral_overlap'] = 'PASS'

    # Calibration sources: Gaia stars when a catalog is available; otherwise bright,
    # isolated sources detected in epoch A. Those are assumed mostly constant and
    # stationary: they calibrate registration, scale and noise, never identify residuals.
    cal_a, cal_b = cat_a, cat_b
    calibration_source = 'Gaia DR3 (proper motion propagated to each epoch)'
    source_mask = np.zeros(shape, bool)
    if not catalog_available:
        found = detect_sources(np.nan_to_num(Aimg), np.nan_to_num(Avar, nan=1.), valid, sigma_threshold=policy['detection_sigma'],
                               min_pixels=2, correlated=True)['sources']
        calibration_source = 'image-detected sources in epoch A (Gaia unavailable; assumed mostly constant and stationary)'
        if found:
            xy = np.array([[f['x'], f['y']] for f in found])
            source_mask |= _disk_mask(shape, xy[:, 0], xy[:, 1], radius)
            gaps = np.hypot(xy[:, None, 0]-xy[None, :, 0], xy[:, None, 1]-xy[None, :, 1]); np.fill_diagonal(gaps, np.inf)
            keep = np.array([f['peak_marginal_sigma'] >= 20 for f in found]) & (gaps.min(axis=1) >= 5)
            for grid in (cat_a, cat_b):   # never let the requested target calibrate itself
                if grid is not None and len(grid['x']):
                    keep &= np.min(np.hypot(xy[:, None, 0]-grid['x'][None], xy[:, None, 1]-grid['y'][None]), axis=1) > 3
            if keep.any():
                cal_a = cal_b = {'ids': [f'IMG{k:04d}' for k in range(int(keep.sum()))], 'g_mag': np.full(int(keep.sum()), np.nan),
                                 'x': xy[keep, 0], 'y': xy[keep, 1]}
    for grid in (cat_a, cat_b):
        if grid is not None:
            source_mask |= _disk_mask(shape, grid['x'], grid['y'], radius)
    blank_sky = interior & ~source_mask
    report['calibration_sources'] = calibration_source

    # Difference-image background: the two independent sky fits leave a residual
    # pedestal/gradient; a clipped plane on blank sky of B-A removes it before photometry.
    pedestal = {'applied': False, 'reason': 'Too few blank interior pixels.'}
    if blank_sky.sum() >= policy['min_null_pixels']:
        D0 = Bimg - Aimg; V0 = Bvar + Avar
        fit_d = fit_background(np.nan_to_num(D0), np.where(np.isfinite(V0) & (V0 > 0), V0, 1.), blank_sky, policy)
        plane = fit_d['background']
        plane_var = np.einsum('...i,ij,...j->...', fit_d['design'], fit_d['coefficient_covariance'], fit_d['design'])
        Bimg = np.where(valid, Bimg - plane, np.nan); Bvar = np.where(valid, Bvar + plane_var, np.nan)
        pedestal = {'applied': True, 'blank_pixels': fit_d['blank_pixel_count'], 'median_mjy_sr': float(np.median(plane[valid])),
                    'coefficients_mjy_sr': fit_d['coefficients'].tolist()}
    report['difference_background'] = pedestal

    # Photometric scale from high-S/N catalog stars in the matched region.
    stars = _measure_stars(Aimg, Avar, Bimg, Bvar, valid, cal_a, cal_b, half, policy['star_min_snr'])
    bright = [s for s in stars if s['snr_a'] >= 100]
    basis = bright if len(bright) >= 5 else stars
    median_ratio, mad_ratio = _robust([s['flux_b']/s['flux_a'] for s in basis])
    scale, scale_var = 1., 0.
    photometry = {'stars': len(stars), 'scale_stars': len(basis),
                  'scale_star_selection': 'aperture S/N>=100' if basis is bright else 'aperture S/N>=20',
                  'median_ratio_b_over_a': median_ratio}
    if len(basis) >= 5:
        se = 1.2533*mad_ratio/np.sqrt(len(basis))
        photometry['median_ratio_standard_error'] = float(se)
        if abs(median_ratio-1) > 3*se:
            scale, scale_var = median_ratio, se**2
            photometry['decision'] = 'Scale differs from unity by >3 standard errors; applied with propagated uncertainty.'
        else:
            photometry['decision'] = 'Consistent with unity; no rescaling applied.'
        gates['photometric_consistency'] = 'PASS'
    else:
        photometry['decision'] = 'Fewer than five matched-region stars; absolute MJy/sr calibration assumed, unverified.'
        gates['photometric_consistency'] = 'UNVERIFIED'
    photometry['scale_applied'] = float(scale); photometry['scale_variance'] = float(scale_var)
    report['photometry'] = photometry

    D = Bimg - scale*Aimg
    slope_term = policy['spectral_slope_bound']*np.abs(delta)/a['wave']*np.fmax(np.abs(Aimg), np.abs(Bimg))
    V = Bvar + scale**2*Avar + scale_var*Aimg**2 + slope_term**2
    S = np.divide(D, np.sqrt(V), out=np.full(shape, np.nan), where=valid & (V > 0))

    # Null test on blank interior sky; inflate-only renormalisation of pixel variances.
    null = {'pixels': int(blank_sky.sum())}
    noise_scale = 1.
    if blank_sky.sum() >= policy['min_null_pixels']:
        center_s, std_s = _robust(S[blank_sky])
        tail = int(np.sum(np.abs(S[blank_sky]) >= policy['detection_sigma']))
        noise_scale = max(1., std_s)
        null.update(median_sigma=center_s, robust_std_sigma=std_s, noise_scale_applied=noise_scale,
                    tail_count=tail, gaussian_expected_tail=float(blank_sky.sum()*5.733e-7),
                    tail_count_after_rescale=int(np.sum(np.abs(S[blank_sky]/noise_scale) >= policy['detection_sigma'])))
        gates['noise_model'] = 'PASS' if std_s <= policy['max_noise_scale'] else 'FAIL'
    else:
        null['reason'] = 'Too few blank interior pixels for a null test.'
        gates['noise_model'] = 'UNVERIFIED'
    report['null_test'] = null
    S = S/noise_scale; Vcal = V*noise_scale**2

    # Empirical source-level noise model from (mostly constant) catalog stars:
    # var(dF) = stat^2 + beta^2 + (alpha F)^2; alpha (multiplicative) from S/N>=100 stars,
    # beta (additive, aperture flux units) from fainter stars given alpha.
    for s in stars:
        s['scaled_flux_a'] = scale*s['flux_a']
        s['residual'] = s['flux_b'] - s['scaled_flux_a']
        s['stat_sigma'] = float(noise_scale*np.sqrt(s['var_b'] + scale**2*s['var_a'] + scale_var*s['flux_a']**2))
        s['fractional_change'] = s['residual']/s['scaled_flux_a']
    bright = [s for s in stars if s['snr_a'] >= 100]; faint = [s for s in stars if s['snr_a'] < 100]
    alpha = excess_scatter([s['fractional_change'] for s in bright], [s['stat_sigma']/s['scaled_flux_a'] for s in bright]) if len(bright) >= 5 else None
    beta = excess_scatter([s['residual'] for s in faint], [np.hypot(s['stat_sigma'], (alpha or 0.)*s['scaled_flux_a']) for s in faint]) if len(faint) >= 5 else None
    for s in stars:
        s['total_sigma'] = float(np.sqrt(s['stat_sigma']**2 + (beta or 0.)**2 + ((alpha or 0.)*s['scaled_flux_a'])**2))
        s['change_z'] = float(s['residual']/s['total_sigma'])
    check = {}
    for name, members in (('all', stars), ('snr_20_50', [s for s in stars if s['snr_a'] < 50]),
                          ('snr_50_100', [s for s in stars if 50 <= s['snr_a'] < 100]), ('snr_ge_100', bright)):
        if len(members) >= 5:
            check[name] = {'stars': len(members), 'robust_std_z': _robust([s['change_z'] for s in members])[1]}
    gates['source_noise_model'] = 'PASS' if len(stars) >= policy['min_noise_model_stars'] and check.get('all', {}).get('robust_std_z') is not None and 0.7 <= check['all']['robust_std_z'] <= 1.4 else 'UNVERIFIED'
    report['systematic_floor'] = {
        'stars': len(stars), 'alpha_fractional': alpha, 'beta_aperture_mjy_sr_px': beta,
        'median_fractional_change': _robust([s['fractional_change'] for s in stars])[0] if stars else None,
        'calibration_check': check,
        'stars_beyond_threshold': int(sum(abs(s['change_z']) >= policy['detection_sigma'] for s in stars)),
        'method': 'Catalog stars (aperture S/N>=20) in the matched region, assumed mostly constant: robust excess scatter beyond propagated errors; multiplicative term from S/N>=100 stars, additive term from fainter stars.'}
    floor = alpha

    # Registration residual from catalog-star centroids (overlap region, wavelength-independent).
    offsets = []
    if cal_a is not None:
        reg_valid = A['valid'] & aligned['valid']
        Afull = np.where(reg_valid, A['image'], np.nan); Bfull = np.where(reg_valid, aligned['image'], np.nan)
        vfull = np.where(reg_valid, A['variance'], np.nan)
        index_b = {ident: k for k, ident in enumerate(cal_b['ids'])}
        for k, ident in enumerate(cal_a['ids']):
            j = index_b.get(ident)
            if j is None or ident == TARGET_ID:
                continue
            ap = _aperture(Afull, vfull, reg_valid, cal_a['x'][k], cal_a['y'][k], half)
            if ap is None or ap[0] <= 0 or ap[0]/np.sqrt(ap[1]) < 50:
                continue
            ca = _centroid(Afull, reg_valid, cal_a['x'][k], cal_a['y'][k], half)
            cb = _centroid(Bfull, reg_valid, cal_b['x'][j], cal_b['y'][j], half)
            if ca and cb:
                offsets.append(((cb[0]-cal_b['x'][j])-(ca[0]-cal_a['x'][k]), (cb[1]-cal_b['y'][j])-(ca[1]-cal_a['y'][k])))
    if len(offsets) >= policy['registration_min_stars']:
        o = np.asarray(offsets); med = np.median(o, axis=0)
        report['registration'] = {'stars': len(offsets), 'median_offset_px': med.tolist(),
                                  'median_offset_arcsec': (med*_pixel_scale_arcsec(a['wcs'])).tolist(),
                                  'robust_scatter_px': (1.4826*np.median(np.abs(o-med), axis=0)).tolist()}
        gates['registration'] = 'PASS' if np.hypot(*med) <= policy['registration_tolerance_px'] else 'FAIL'
    else:
        report['registration'] = {'stars': len(offsets), 'reason': 'Fewer than the minimum bright calibration sources for a registration check.'}
        gates['registration'] = 'UNVERIFIED'

    # Candidate residual components.
    positive = detect_sources(np.nan_to_num(D), np.nan_to_num(Vcal, nan=1.), interior, sigma_threshold=policy['detection_sigma'],
                              min_pixels=policy['detection_min_pixels'], correlated=True)['sources']
    negative = detect_sources(np.nan_to_num(-D), np.nan_to_num(Vcal, nan=1.), interior, sigma_threshold=policy['detection_sigma'],
                              min_pixels=policy['detection_min_pixels'], correlated=True)['sources']
    stars_by_id = {s['source_id']: s for s in stars}
    candidates = []
    for sign, group, other in ((1, positive, negative), (-1, negative, positive)):
        for component in group:
            cx, cy = int(np.rint(component['x'])), int(np.rint(component['y']))
            window = (slice(max(cy-half, 0), cy+half+1), slice(max(cx-half, 0), cx+half+1))
            component['aperture_sum'] = float(np.nansum(D[window]))
            c = _classify(component, sign, stars_by_id, cat_a, cat_b, other, floor, noise_scale, policy, catalog_available, dist_to_edge)
            ra, dec = a['wcs'].pixel_to_world_values(c['x'], c['y'])
            c.update(ra_deg=float(ra), dec_deg=float(dec), difference_sum_mjy_sr_px=component['aperture_sum'])
            candidates.append(c)
    # Residuals without a catalog counterpart: is a source present in one epoch or both?
    flagged = ~(sa['valid'] & a['calibration_valid'])
    for c in candidates:
        if c['classification'] != 'UNASSOCIATED_CANDIDATE':
            continue
        fa = _aperture(Aimg, Avar, valid, c['x'], c['y'], half); fb = _aperture(Bimg, Bvar, valid, c['x'], c['y'], half)
        if fa is None or fb is None:
            c.update(classification='UNASSOCIATED_RESIDUAL', reason='Aperture not fully valid in both epochs; presence cannot be tested.')
            continue
        snr_a, snr_b = fa[0]/np.sqrt(fa[1])/noise_scale, fb[0]/np.sqrt(fb[1])/noise_scale
        c.update(aperture_snr_a=float(snr_a), aperture_snr_b=float(snr_b))
        if snr_a >= policy['detection_sigma'] and snr_b >= policy['detection_sigma']:
            residual = fb[0] - scale*fa[0]
            total = float(np.sqrt(noise_scale**2*(fb[1]+scale**2*fa[1]+scale_var*fa[0]**2) + (beta or 0.)**2 + ((alpha or 0.)*scale*fa[0])**2))
            star = {'x_a': c['x'], 'y_a': c['y']}
            vetting = _vet_variable(star, flagged, valid, dist_to_boundary, Aimg, D, cal_a, policy)
            if len(stars) < policy['min_noise_model_stars']:
                vetting['failed'].append('noise_model_calibrated_on_too_few_stars')
            c.update(fractional_change=float(residual/(scale*fa[0])), fractional_change_sigma=float(total/(scale*fa[0])),
                     change_z=float(residual/total), vetting=vetting)
            if abs(residual/total) >= policy['detection_sigma'] and not vetting['failed']:
                c.update(classification='UNCATALOGED_VARIABLE_CANDIDATE', reason='Source present in both epochs, absent from Gaia DR3 (G<=19), with a vetted flux change beyond the calibrated noise model.')
            else:
                c.update(classification='PERSISTENT_SOURCE_RESIDUAL', reason='Source present in both epochs (not in Gaia DR3 G<=19); residual explained by the calibrated noise model or failed vetting.')
        elif max(snr_a, snr_b) >= policy['detection_sigma'] and min(snr_a, snr_b) < 3:
            c.update(classification='TRANSIENT_CANDIDATE', present_in='A' if snr_a > snr_b else 'B',
                     reason='Point source detected in only one epoch with no Gaia DR3 counterpart: transient, solar-system object or unflagged artifact. Follow-up required.')
        else:
            c.update(classification='UNASSOCIATED_RESIDUAL', reason='No significant aperture source in either epoch; consistent with a PSF/noise residual.')
    candidates.sort(key=lambda c: -abs(c['peak_sigma']))
    for i, c in enumerate(candidates):
        c['id'] = i+1
    report['candidates'] = candidates
    variables = []
    flagged = ~(sa['valid'] & a['calibration_valid'])
    for star in stars:
        if abs(star['change_z']) < policy['detection_sigma']:
            continue
        vetting = _vet_variable(star, flagged, valid, dist_to_boundary, Aimg, D, cal_a, policy)
        if len(stars) < policy['min_noise_model_stars']:
            vetting['failed'].append('noise_model_calibrated_on_too_few_stars')
        ra, dec = a['wcs'].pixel_to_world_values(star['x_a'], star['y_a'])
        variables.append({'source_id': star['source_id'], 'g_mag': star['g_mag'], 'ra_deg': float(ra), 'dec_deg': float(dec),
                          'x': star['x_a'], 'y': star['y_a'], 'flux_a_mjy_sr_px': star['scaled_flux_a'], 'flux_b_mjy_sr_px': star['flux_b'],
                          'fractional_change': star['fractional_change'], 'fractional_change_sigma': star['total_sigma']/star['scaled_flux_a'],
                          'change_z': star['change_z'], 'snr_a': star['snr_a'], 'vetting': vetting,
                          'classification': 'SOURCE_LEVEL_VARIABILITY_CANDIDATE' if not vetting['failed'] else 'REJECTED_BY_VETTING',
                          'interpretation': 'Gaia source whose matched-wavelength aperture flux changed beyond the empirically calibrated noise model; requires independent confirmation (blending, flagged neighbours, intrinsic variability).'})
    variables.sort(key=lambda v: (v['classification'] != 'SOURCE_LEVEL_VARIABILITY_CANDIDATE', -abs(v['change_z'])))
    report['variability_candidates'] = variables
    report['variability_summary'] = {'stars_tested': len(stars), 'beyond_threshold': len(variables),
                                     'passed_vetting': sum(v['classification'] == 'SOURCE_LEVEL_VARIABILITY_CANDIDATE' for v in variables)}
    report['candidate_summary'] = {k: sum(c['classification'] == k for c in candidates) for k in sorted({c['classification'] for c in candidates})}

    report['injection_recovery'] = _injection_recovery(np.nan_to_num(D), np.nan_to_num(Vcal, nan=1.), interior, blank_sky, effective_kernel, policy)

    if target_coord is not None:
        pa = _project_target(target, a['date'], a['wcs']); pb = _project_target(target, b['date'], a['wcs'])
        target_info = {'position_epoch_a_px': list(pa), 'position_epoch_b_px': list(pb)}
        ia, ib = _aperture(Aimg, Avar, valid, *pa, half), _aperture(Bimg, Bvar, valid, *pb, half)
        if ia is None or ib is None:
            cx, cy = int(np.rint(pa[0])), int(np.rint(pa[1]))
            window = (slice(max(cy-half, 0), cy+half+1), slice(max(cx-half, 0), cx+half+1))
            inside = 0 <= cx < shape[1] and 0 <= cy < shape[0]
            bits = sorted({b for v in np.unique(sa['flags'][window]) for b in range(31) if int(v) & (1 << b) and b != 21}) if inside else []
            if inside and not sa['valid'][window].all():
                target_info.update(status='TARGET_PIXELS_FLAGGED', flag_bits=bits,
                                   reason=f'Target core pixels carry non-SOURCE quality flag bits {bits} (e.g. saturation/non-linearity); excluded by the mask policy.')
            elif inside and not spectral[window].all():
                target_info.update(status='OUTSIDE_MATCHED_REGION', reason='Target lies where the two epochs sample different wavelengths.')
            else:
                target_info.update(status='OUTSIDE_VALID_REGION', reason='Target aperture is outside the jointly valid comparison region.')
        elif ia[0] <= 0:
            target_info.update(status='NOT_DETECTED_IN_A', flux_a=ia[0])
        else:
            residual = ib[0] - scale*ia[0]
            total = float(np.sqrt(noise_scale**2*(ib[1]+scale**2*ia[1]+scale_var*ia[0]**2) + (beta or 0.)**2 + ((alpha or 0.)*scale*ia[0])**2))
            target_info.update(status='MEASURED', fractional_change=float(residual/(scale*ia[0])), fractional_change_sigma=float(total/(scale*ia[0])),
                               change_z=float(residual/total), snr_a=float(ia[0]/np.sqrt(ia[1])))
    report['target'] = target_info
    report['gates'] = gates
    failed = [k for k, v in gates.items() if v == 'FAIL']
    report['status'] = 'REJECTED' if failed else ('VALIDATED' if all(v == 'PASS' for v in gates.values()) else 'PARTIAL')
    report['reason'] = ('Failed gates: '+', '.join(failed)) if failed else (
        'All validation gates passed.' if report['status'] == 'VALIDATED'
        else 'Comparison computed; unverified gates: '+', '.join(k for k, v in gates.items() if v == 'UNVERIFIED'))
    report['limitations'] = [
        'Only pixels whose own-calibration CWAVE values agree within the policy tolerance are compared; elsewhere the epochs sample different wavelengths.',
        'Bilinear reprojection and PSF cross-convolution correlate neighbouring pixels; per-pixel variances are marginal and the blank-sky null test only inflates them.',
        'The spectral-mismatch term bounds |dlnF/dlnlambda| by the policy value; extreme spectra can exceed it.',
        'Single detector-zone ePSFs are used at the field centre; PSF variation across the cutout contributes to the empirical systematic floor.',
        'Catalog association is positional (Gaia DR3, G<=19, proper motion propagated, no parallax); candidates require independent follow-up.',
    ]
    return _result(report, a, Aimg, Bimg, D, Vcal, S, valid, interior, blank_sky, rel, stars=stars)


def _pixel_scale_arcsec(wcs):
    return float(np.sqrt(abs(np.linalg.det(wcs.pixel_scale_matrix)))*3600)


def _project_target(target, date, wcs):
    pm = target.get('pm_ra_cosdec_mas_yr'), target.get('pm_dec_mas_yr')
    if pm[0] is not None and pm[1] is not None and target.get('epoch') == 'J2000':
        p = propagate_catalog([{'ra': target['ra_deg'], 'dec': target['dec_deg'], 'pmra': pm[0], 'pmdec': pm[1], 'ref_epoch': 2000.0, 'source_id': 'target'}], date)
        ra, dec = float(p['ra'][0]), float(p['dec'][0])
    else:
        ra, dec = float(target['ra_deg']), float(target['dec_deg'])
    x, y = wcs.world_to_pixel_values(ra, dec)
    return float(x), float(y)


def _result(report, a, A, B, D, V, S, valid, interior, blank, rel, stars=()):
    report['measured_stars'] = list(stars)
    return {'report': report, 'wcs': a['wcs'], 'wave': a['wave'],
            'arrays': {'image_a': A, 'image_b_aligned': B, 'difference': D, 'variance': V, 'significance': S,
                       'matched_mask': valid, 'interior_mask': interior, 'null_mask': blank, 'delta_over_band': rel}}


# ----------------------------------------------------------------------------- products

def _b64(array):
    return base64.b64encode(np.ascontiguousarray(array, dtype='<f4').tobytes()).decode('ascii')


def _json_safe(value):
    if isinstance(value, dict):
        return {str(k): _json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(v) for v in value]
    if isinstance(value, (np.floating, float)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, np.bool_):
        return bool(value)
    if isinstance(value, np.ndarray):
        return _json_safe(value.tolist())
    return value


def write_products(result, output_dir):
    """FITS, PNG, JSON report and a compact viewer payload. Returns written file names."""
    out = Path(output_dir); out.mkdir(parents=True, exist_ok=True)
    arrays = result['arrays']; report = _json_safe(result['report']); header = result['wcs'].to_header(relax=True)
    hdus = [fits.PrimaryHDU()]
    hdus[0].header['STATUS'] = report['status']; hdus[0].header['OBS_A'] = report['observation_a'][-40:]; hdus[0].header['OBS_B'] = report['observation_b'][-40:]
    hdus[0].header['COMMENT'] = 'Spectrally matched SPHEREx difference products; residuals are not discoveries.'
    for name, key, unit in (('IMAGE_A', 'image_a', 'MJy/sr'), ('IMAGE_B', 'image_b_aligned', 'MJy/sr'), ('DIFF', 'difference', 'MJy/sr'),
                            ('VARIANCE', 'variance', 'MJy2/sr2'), ('SIGNIF', 'significance', '1'), ('DLAMBDA', 'delta_over_band', '1')):
        h = header.copy(); h['BUNIT'] = unit
        hdus.append(fits.ImageHDU(np.asarray(arrays[key], dtype='>f4'), header=h, name=name))
    for name, key in (('MATCHED', 'matched_mask'), ('INTERIOR', 'interior_mask'), ('NULLSKY', 'null_mask')):
        hdus.append(fits.ImageHDU(np.asarray(arrays[key], dtype='uint8'), header=header.copy(), name=name))
    fits.HDUList(hdus).writeto(out/'comparison.fits', overwrite=True, checksum=True)
    (out/'report.json').write_text(json.dumps(report, indent=2, allow_nan=False), encoding='utf-8')
    shape = arrays['image_a'].shape
    yy, xx = np.indices(shape)
    ra, dec = result['wcs'].pixel_to_world_values(xx, yy)
    viewer = {'shape': list(shape), 'encoding': 'base64 little-endian float32, row-major, NaN = invalid',
              'planes': {k: _b64(np.where(np.isfinite(arrays[k]), arrays[k], np.nan)) for k in ('image_a', 'image_b_aligned', 'difference', 'significance', 'delta_over_band')},
              'masks': {k: _b64(arrays[k].astype(float)) for k in ('matched_mask', 'interior_mask')},
              'ra': _b64(ra), 'dec': _b64(dec), 'wave_um': _b64(result['wave'])}
    (out/'viewer.json').write_text(json.dumps(viewer), encoding='utf-8')
    _render(result, out/'comparison.png')
    rows = ['id,classification,sign,peak_sigma,ra_deg,dec_deg,x,y,pixel_count,catalog_source_id,catalog_g_mag,catalog_separation_px,fractional_change,fractional_change_sigma']
    for c in report.get('candidates', []):
        rows.append(','.join('' if c.get(k) is None else str(c.get(k)) for k in rows[0].split(',')))
    (out/'candidates.csv').write_text('\n'.join(rows)+'\n', encoding='utf-8')
    return ['comparison.fits', 'report.json', 'viewer.json', 'comparison.png', 'candidates.csv']


def _render(result, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    arrays = result['arrays']; report = result['report']
    fig, axes = plt.subplots(1, 4, figsize=(17, 4.6))
    finite = np.concatenate([arrays['image_a'][np.isfinite(arrays['image_a'])], arrays['image_b_aligned'][np.isfinite(arrays['image_b_aligned'])]])
    lo, hi = (np.percentile(finite, [2, 99.5]) if finite.size else (0, 1))
    panels = (('image_a', f"A  {report['epochs']['a']['date'][:10]}", 'gray', lo, hi), ('image_b_aligned', f"B on A grid  {report['epochs']['b']['date'][:10]}", 'gray', lo, hi),
              ('difference', 'B - sA (MJy/sr)', 'RdBu_r', None, None), ('significance', 'Significance (sigma)', 'RdBu_r', -6, 6))
    for ax, (key, title, cmap, vmin, vmax) in zip(axes, panels):
        data = arrays[key]
        if key == 'difference':
            finite_d = data[np.isfinite(data)]; v = np.percentile(np.abs(finite_d), 99) if finite_d.size else 1; vmin, vmax = -v, v
        image = ax.imshow(data, origin='lower', cmap=cmap, vmin=vmin, vmax=vmax, interpolation='nearest')
        ax.contour(arrays['matched_mask'], levels=[.5], colors='#f5a623', linewidths=.6)
        ax.set_title(title, fontsize=10); ax.set_xticks([]); ax.set_yticks([])
        fig.colorbar(image, ax=ax, shrink=.75)
    for c in report.get('candidates', [])[:40]:
        axes[3].add_patch(plt.Circle((c['x'], c['y']), 3, fill=False, color='#ffd400' if c['classification'] in ('CANDIDATE_VARIABLE', 'UNASSOCIATED_CANDIDATE') else '#888888', lw=1))
    fig.suptitle(f"{report['status']}: {report['observation_a'].split('?')[-1]} vs {report['observation_b'].split('?')[-1]}  (orange contour: spectrally matched region; residuals are not discoveries)", fontsize=10)
    fig.tight_layout(); fig.savefig(path, dpi=110); plt.close(fig)
