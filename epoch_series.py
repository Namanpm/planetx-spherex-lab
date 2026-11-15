"""Same-wavelength multi-epoch light curve and astrometric trajectory of one target.

Input: real observations on one detector whose screened CWAVE at the target agree
(``spectral_matching.matched_groups``). Each epoch's own calibration confirms the
target wavelength. Measurements are measured-ePSF fits on native pixels.

Two analyses are kept separate:
  * photometry - target flux per epoch divided by a per-epoch zero point from
    field stars that are themselves wavelength-matched to the reference epoch;
  * astrometry - fitted positions in a common tangent plane, corrected by the
    per-epoch frame offset of image-detected field stars, with an empirical
    floor from those stars' epoch-to-epoch repeatability. The trajectory is
    compared with catalog proper motion/parallax (JPL Horizons observer vectors
    when available). It is a test of known motion, never a discovery claim.
"""
import json
import time
from pathlib import Path

import numpy as np
from astropy.coordinates import SkyCoord

from difference_imaging import _disk_mask, _json_safe, _robust, excess_scatter, fit_background, propagate_catalog, _catalog_on_grid
from native_photometry import fit_stamp, local_pixel_area_sr, sample_epsf
from psf_processing import extract_prepared_epsf
from source_detection import detect_sources

POLICY = {
    'tolerance_fraction': 0.10,
    'stamp_half_width_px': 3,
    'search_radius_px': 2.0,
    'field_star_threshold_sigma': 15.0,
    'field_star_isolation_px': 6.0,
    'max_field_stars': 60,
    'min_field_stars': 8,                  # stars for a zero point with a standard error
    'min_zero_point_stars': 3,             # below 8: zero point error = robust star scatter (not /sqrt n)
    'min_frame_stars_per_epoch': 3,
    'detection_sigma': 5.0,
}


def psf_fit(observation, image, variance, valid, x0, y0, policy=POLICY, search=None):
    """Measured-ePSF position+flux fit on native pixels with a local plane.

    Grid search (0.25 px, then 0.05 px) over sub-pixel position; amplitude and
    plane are linear (native_photometry.fit_stamp). Position covariance from the
    chi-square curvature (delta chi2 = 1). Returns None if the stamp is unusable.
    """
    search = policy['search_radius_px'] if search is None else search
    half = policy['stamp_half_width_px']
    cx, cy = int(np.rint(x0)), int(np.rint(y0))
    if cx-half-2 < 0 or cy-half-2 < 0 or cx+half+2 >= image.shape[1] or cy+half+2 >= image.shape[0]:
        return None
    psf = extract_prepared_epsf(observation['epsf'], float(np.clip(x0, 0, image.shape[1]-1)), float(np.clip(y0, 0, image.shape[0]-1)))
    raw = psf['oversampled_epsf']; ox, oy = psf['provenance']['oversampling_xy']
    window = (slice(cy-half, cy+half+1), slice(cx-half, cx+half+1))
    yy, xx = np.mgrid[window]
    img, var, ok = image[window], variance[window], valid[window]
    if ok.sum() < (2*half+1)**2*0.8:
        return None

    def chi2_at(x, y):
        profile, inside = sample_epsf(raw, (ox, oy), xx-x, yy-y)
        total = np.nansum(profile)
        if not np.isfinite(total) or total <= 0:
            return np.inf, None
        try:
            fit = fit_stamp(img, var, ok & inside, profile/total, xx-x, yy-y)
        except ValueError:
            return np.inf, None
        return fit['chi2'], fit

    def grid(xc, yc, radius, step):
        offsets = np.arange(-radius, radius+step/2, step)
        values = np.full((offsets.size, offsets.size), np.inf)
        for i, dy in enumerate(offsets):
            for j, dx in enumerate(offsets):
                values[i, j] = chi2_at(xc+dx, yc+dy)[0]
        i, j = np.unravel_index(np.argmin(values), values.shape)
        return xc+offsets[j], yc+offsets[i], values, offsets, (i, j)

    x1, y1, *_ = grid(float(x0), float(y0), search, 0.25)
    x2, y2, values, offsets, (i, j) = grid(x1, y1, 0.25, 0.05)
    best, fit = chi2_at(x2, y2)
    if fit is None or not np.isfinite(best):
        return None
    # Quadratic fit of chi2 around the minimum (5x5 grid points) -> Hessian.
    lo_i, hi_i = max(i-2, 0), min(i+3, offsets.size); lo_j, hi_j = max(j-2, 0), min(j+3, offsets.size)
    sub = values[lo_i:hi_i, lo_j:hi_j]; dy, dx = np.meshgrid(offsets[lo_i:hi_i]-offsets[i], offsets[lo_j:hi_j]-offsets[j], indexing='ij')
    good = np.isfinite(sub)
    design = np.column_stack([np.ones(good.sum()), dx[good], dy[good], dx[good]**2, dx[good]*dy[good], dy[good]**2])
    try:
        coef = np.linalg.lstsq(design, sub[good], rcond=None)[0]
        hessian = np.array([[2*coef[3], coef[4]], [coef[4], 2*coef[5]]])
        covariance = 2*np.linalg.inv(hessian)
        if not np.all(np.linalg.eigvalsh(covariance) > 0):
            raise np.linalg.LinAlgError
    except np.linalg.LinAlgError:
        return None
    reduced = fit['chi2']/max(fit['dof'], 1)
    covariance *= max(reduced, 1.)  # inflate-only by reduced chi-square
    area = local_pixel_area_sr(observation['wcs'], x2, y2)
    flux = fit['amplitude']*area*1e6; flux_var = fit['variance']*(area*1e6)**2*max(reduced, 1.)
    flagged = observation['science']['flags'][window]
    return {'x': float(x2), 'y': float(y2), 'position_covariance_px2': covariance.tolist(),
            'flux_jy': float(flux), 'flux_variance_jy2': float(flux_var), 'reduced_chi2': float(reduced),
            'core_flag_values': sorted({int(v) for v in np.unique(flagged[half-1:half+2, half-1:half+2])}),
            'psf_zone': psf['provenance']['zone_index']}


def _tangent_arcsec(ra, dec, ra0, dec0):
    ra, dec, ra0, dec0 = map(np.deg2rad, (np.asarray(ra, float), np.asarray(dec, float), ra0, dec0))
    cosc = np.sin(dec0)*np.sin(dec)+np.cos(dec0)*np.cos(dec)*np.cos(ra-ra0)
    xi = np.cos(dec)*np.sin(ra-ra0)/cosc
    eta = (np.cos(dec0)*np.sin(dec)-np.sin(dec0)*np.cos(dec)*np.cos(ra-ra0))/cosc
    return np.rad2deg(xi)*3600, np.rad2deg(eta)*3600


def _pixel_cov_to_sky(wcs, x, y, cov_px):
    """Linear map of a pixel covariance to (east, north) arcsec^2 via a numeric WCS Jacobian."""
    c = wcs.pixel_to_world(x, y)
    jac = []
    for dx, dy in ((0.5, 0), (0, 0.5)):
        p = wcs.pixel_to_world(x+dx, y+dy); m = wcs.pixel_to_world(x-dx, y-dy)
        e_p, n_p = c.spherical_offsets_to(p); e_m, n_m = c.spherical_offsets_to(m)
        jac.append([(e_p-e_m).arcsec/(2*max(dx, dy)), (n_p-n_m).arcsec/(2*max(dx, dy))])
    j = np.array(jac).T
    return (j @ np.asarray(cov_px) @ j.T).tolist()


def analyze_series(epochs, *, target=None, position=None, catalog_sources=None, policy=POLICY, observer_vectors=None):
    """Light curve + trajectory from loaded epochs (difference_imaging.load_epoch)."""
    epochs = sorted(epochs, key=lambda e: e['mjd'])
    if len({e['detector'] for e in epochs}) != 1:
        raise ValueError('series epochs must share one detector')
    if len(epochs) < 2:
        raise ValueError('at least two epochs required')
    reference = epochs[0]
    prepared = []
    for e in epochs:
        s = e['science']
        cat = _catalog_on_grid(propagate_catalog(catalog_sources, e['date']), e['wcs'], s['image'].shape) if catalog_sources else None
        blank = s['valid'] & (s['flags'] == 0) & e['calibration_valid']
        if cat is not None:
            blank &= ~_disk_mask(s['image'].shape, cat['x'], cat['y'], 2.5)
        bg = fit_background(s['image'], s['variance'], blank)
        prepared.append({'epoch': e, 'image': bg['residual'], 'variance': bg['residual_variance'], 'valid': s['valid'] & np.isfinite(bg['residual'])})

    # Target predicted position per epoch (catalog PM; parallax handled in the trajectory model).
    def predicted(e):
        if target is not None and target.get('pm_ra_cosdec_mas_yr') is not None:
            p = propagate_catalog([{'ra': target['ra_deg'], 'dec': target['dec_deg'], 'pmra': target['pm_ra_cosdec_mas_yr'],
                                    'pmdec': target['pm_dec_mas_yr'], 'ref_epoch': 2000.0, 'source_id': 'target'}], e['date'])
            return float(p['ra'][0]), float(p['dec'][0])
        return float(position[0]), float(position[1])

    # Field stars detected in the reference epoch (catalog-free), isolated and unflagged.
    ref = prepared[0]; rs = reference['science']
    det_valid = ref['valid'] & ((rs['flags'] == 0) | (rs['flags'] == (1 << 21)))
    found = detect_sources(np.nan_to_num(ref['image']), np.nan_to_num(ref['variance'], nan=1.), det_valid,
                           sigma_threshold=policy['field_star_threshold_sigma'], min_pixels=3, correlated=True)['sources']
    tra, tdec = predicted(reference); tx, ty = reference['wcs'].world_to_pixel_values(tra, tdec)
    stars = []
    for star in found:
        others = [o for o in found if o is not star]
        if any(np.hypot(o['x']-star['x'], o['y']-star['y']) < policy['field_star_isolation_px'] for o in others):
            continue
        if np.hypot(star['x']-tx, star['y']-ty) < policy['field_star_isolation_px']:
            continue
        stars.append(star)
        if len(stars) >= policy['max_field_stars']:
            break
    star_sky = [reference['wcs'].pixel_to_world_values(s['x'], s['y']) for s in stars]

    measurements = []
    for p in prepared:
        e = p['epoch']; shape = e['science']['image'].shape
        entry = {'observation_id': e['id'], 'date': e['date'], 'mjd': e['mjd'], 'detector': e['detector'],
                 'calibration_url': e['calibration']['calibration_url']}
        ra, dec = predicted(e); x, y = e['wcs'].world_to_pixel_values(ra, dec)
        fit = psf_fit(e, p['image'], p['variance'], p['valid'], float(x), float(y), policy)
        if fit is None:
            cx, cy = int(np.rint(x)), int(np.rint(y))
            inside = 0 <= cx < shape[1] and 0 <= cy < shape[0]
            bits = sorted({b for v in np.unique(e['science']['flags'][max(cy-1, 0):cy+2, max(cx-1, 0):cx+2]) for b in range(31) if int(v) & (1 << b) and b != 21}) if inside else []
            entry['target'] = {'status': 'UNMEASURED', 'reason': f'Target stamp invalid or flagged (non-SOURCE bits {bits}).' if bits else 'Target stamp outside cutout or unusable.'}
        else:
            ix, iy = int(np.rint(fit['x'])), int(np.rint(fit['y']))
            fit['wavelength_um'] = float(e['wave'][iy, ix]); fit['bandwidth_um'] = float(e['band'][iy, ix])
            fit['ra_deg'], fit['dec_deg'] = map(float, e['wcs'].pixel_to_world_values(fit['x'], fit['y']))
            fit['position_covariance_arcsec2'] = _pixel_cov_to_sky(e['wcs'], fit['x'], fit['y'], fit['position_covariance_px2'])
            fit['status'] = 'MEASURED' if all(v in (0, 1 << 21) for v in fit['core_flag_values']) else 'FLAGGED'
            entry['target'] = fit
        star_fits = []
        for k, (sra, sdec) in enumerate(star_sky):
            sx, sy = e['wcs'].world_to_pixel_values(sra, sdec)
            f = psf_fit(e, p['image'], p['variance'], p['valid'], float(sx), float(sy), policy, search=1.0)
            if f is None or not all(v in (0, 1 << 21) for v in f['core_flag_values']):
                star_fits.append(None); continue
            ix, iy = int(np.rint(f['x'])), int(np.rint(f['y']))
            f['wavelength_um'] = float(e['wave'][iy, ix]); f['bandwidth_um'] = float(e['band'][iy, ix])
            f['ra_deg'], f['dec_deg'] = map(float, e['wcs'].pixel_to_world_values(f['x'], f['y']))
            f['position_covariance_arcsec2'] = _pixel_cov_to_sky(e['wcs'], f['x'], f['y'], f['position_covariance_px2'])
            star_fits.append(f)
        entry['field_star_fits'] = star_fits
        measurements.append(entry)

    report = {'policy': policy, 'epoch_count': len(epochs), 'detector': reference['detector'],
              'observations': [m['observation_id'] for m in measurements], 'field_stars_detected': len(stars)}
    report['photometry'] = _light_curve(measurements, policy)
    report['astrometry'] = _trajectory(measurements, target, policy, observer_vectors)
    report['epochs'] = [{k: v for k, v in m.items() if k != 'field_star_fits'} for m in measurements]
    report['limitations'] = [
        'Only epochs whose target wavelength agrees within the tolerance are photometrically comparable; others are reported, not used.',
        'Per-epoch zero points use field stars that are wavelength-matched to the reference epoch; their median ratio assumes most are constant.',
        'Astrometric floor is the epoch-to-epoch repeatability of image-detected field stars (their own small proper motions are absorbed).',
        'Expected positions use SIMBAD astrometry; light-time, aberration and gravitational deflection are omitted.',
        'A trajectory consistent with catalog motion validates the method on a known mover; it is not a discovery.']
    return _json_safe(report)


def _light_curve(measurements, policy):
    ref = measurements[0]
    result = {'status': 'INSUFFICIENT_EVIDENCE', 'points': []}
    ref_target = ref['target'] if ref['target'].get('status') == 'MEASURED' else None
    zero_points = []
    for m in measurements:
        ratios = []
        for f0, f in zip(ref['field_star_fits'], m['field_star_fits']):
            if f0 is None or f is None or f0['flux_jy'] <= 0:
                continue
            if abs(f['wavelength_um']-f0['wavelength_um']) > policy['tolerance_fraction']*min(f['bandwidth_um'], f0['bandwidth_um']):
                continue
            ratios.append(f['flux_jy']/f0['flux_jy'])
        if m is ref:
            zero_points.append({'observation_id': m['observation_id'], 'stars': len(ratios), 'ratio': 1.0, 'ratio_standard_error': 0.0,
                                'method': 'reference epoch (gauge)'})
        elif len(ratios) >= policy['min_field_stars']:
            median, mad = _robust(ratios)
            zero_points.append({'observation_id': m['observation_id'], 'stars': len(ratios), 'ratio': median,
                                'ratio_standard_error': 1.2533*mad/np.sqrt(len(ratios)), 'robust_scatter': mad,
                                'method': 'median of wavelength-matched field-star ratios'})
        elif len(ratios) >= policy['min_zero_point_stars']:
            median, mad = _robust(ratios)
            zero_points.append({'observation_id': m['observation_id'], 'stars': len(ratios), 'ratio': median,
                                'ratio_standard_error': max(mad, 1.2533*mad/np.sqrt(len(ratios))), 'robust_scatter': mad,
                                'method': 'few wavelength-matched stars: error set to the full robust star scatter (conservative)'})
        else:
            zero_points.append({'observation_id': m['observation_id'], 'stars': len(ratios), 'ratio': None,
                                'reason': 'Too few wavelength-matched field stars for a zero point.'})
    result['zero_points'] = zero_points
    # Multiplicative floor: excess scatter of matched field-star ratios after the zero point.
    residuals, sigmas = [], []
    for m, zp in zip(measurements[1:], zero_points[1:]):
        if zp['ratio'] is None:
            continue
        for f0, f in zip(ref['field_star_fits'], m['field_star_fits']):
            if f0 is None or f is None or f0['flux_jy'] <= 0 or abs(f['wavelength_um']-f0['wavelength_um']) > policy['tolerance_fraction']*min(f['bandwidth_um'], f0['bandwidth_um']):
                continue
            residuals.append(f['flux_jy']/f0['flux_jy']/zp['ratio']-1)
            sigmas.append(np.sqrt(f['flux_variance_jy2']/f['flux_jy']**2 + f0['flux_variance_jy2']/f0['flux_jy']**2))
    floor = excess_scatter(residuals, sigmas) if len(residuals) >= 5 else None
    result['fractional_floor'] = floor
    if ref_target is None:
        result['reason'] = 'Target not measurable in the reference epoch (flagged, saturated or outside the cutout).'
        return result
    tolerance = policy['tolerance_fraction']*ref_target['bandwidth_um']
    for m, zp in zip(measurements, zero_points):
        t = m['target']
        point = {'observation_id': m['observation_id'], 'date': m['date'], 'mjd': m['mjd']}
        if t.get('status') != 'MEASURED':
            point.update(status='EXCLUDED', reason=t.get('reason', 'target flagged'))
        elif abs(t['wavelength_um']-ref_target['wavelength_um']) > tolerance:
            point.update(status='EXCLUDED', reason=f"Target wavelength {t['wavelength_um']:.4f} um differs from reference {ref_target['wavelength_um']:.4f} um by more than {policy['tolerance_fraction']} x CBAND.")
        elif zp['ratio'] is None:
            point.update(status='EXCLUDED', reason=zp['reason'])
        else:
            flux = t['flux_jy']/zp['ratio']
            stat = np.sqrt(t['flux_variance_jy2'])/zp['ratio']
            zp_term = flux*zp['ratio_standard_error']/zp['ratio']
            sigma = float(np.sqrt(stat**2+zp_term**2+((floor or 0.)*flux)**2))
            point.update(status='USED', wavelength_um=t['wavelength_um'], measured_flux_jy=t['flux_jy'], calibrated_flux_jy=float(flux),
                         flux_sigma_jy=sigma, statistical_sigma_jy=float(stat), zero_point=zp['ratio'])
        result['points'].append(point)
    used = [p for p in result['points'] if p['status'] == 'USED']
    result['used_epochs'] = len(used)
    if len(used) >= 2:
        from motion_inference import fit_variability
        fit = fit_variability([p['calibrated_flux_jy'] for p in used], [p['flux_sigma_jy']**2 for p in used])
        from scipy.stats import chi2
        fit['pvalue'] = float(chi2.sf(fit['chi2_constant'], fit['dof']))
        result['constant_flux_test'] = fit
        result['max_abs_deviation_sigma'] = float(max(abs(p['calibrated_flux_jy']-fit['constant_flux'])/p['flux_sigma_jy'] for p in used))
        result['status'] = 'VARIABILITY_CANDIDATE' if fit['pvalue'] < 1e-6 else 'CONSISTENT_WITH_CONSTANT'
        result['interpretation'] = ('Flux changes exceed the empirical noise model at p<1e-6; a variability candidate requiring independent confirmation.'
                                    if result['status'] == 'VARIABILITY_CANDIDATE' else 'No flux change detected beyond the empirical noise model at this wavelength.')
    else:
        result['reason'] = 'Fewer than two usable same-wavelength epochs.'
    return result


def _trajectory(measurements, target, policy, observer_vectors):
    result = {'status': 'INSUFFICIENT_EVIDENCE'}
    measured = [m for m in measurements if m['target'].get('status') == 'MEASURED']
    # Frame offsets and floor from field stars measured in >=2 epochs: alternating robust
    # solution of position_ke = star_k + frame_e (frame of the reference epoch fixed at 0).
    n_epochs = len(measurements); n_stars = len(measurements[0]['field_star_fits'])
    all_fits = [f for m in measurements for f in m['field_star_fits'] if f is not None]
    if not all_fits:
        result['reason'] = 'No field stars measured.'
        return result
    ra0, dec0 = np.median([[f['ra_deg'] for f in all_fits], [f['dec_deg'] for f in all_fits]], axis=1)
    pos = np.full((n_stars, n_epochs, 2), np.nan); err = np.full((n_stars, n_epochs, 2), np.nan)
    for e, m in enumerate(measurements):
        for k, f in enumerate(m['field_star_fits']):
            if f is None:
                continue
            x, y = _tangent_arcsec(f['ra_deg'], f['dec_deg'], ra0, dec0)
            pos[k, e] = [x, y]
            err[k, e] = [np.sqrt(f['position_covariance_arcsec2'][0][0]), np.sqrt(f['position_covariance_arcsec2'][1][1])]
    usable = np.sum(np.isfinite(pos[:, :, 0]), axis=1) >= 2
    pos, err = pos[usable], err[usable]
    per_epoch = np.sum(np.isfinite(pos[:, :, 0]), axis=0)
    result['field_stars'] = int(usable.sum()); result['field_stars_per_epoch'] = per_epoch.tolist()
    if usable.sum() < policy['min_field_stars'] or np.any(per_epoch < policy['min_frame_stars_per_epoch']):
        result['reason'] = f'Field-star frame underconstrained: {int(usable.sum())} stars in >=2 epochs, per-epoch counts {per_epoch.tolist()}.'
        return result
    frame = np.zeros((n_epochs, 2))
    for _ in range(50):
        star_mean = np.nanmedian(pos - frame[None], axis=1)
        new = np.nanmedian(pos - star_mean[:, None], axis=0); new -= new[0]
        if np.allclose(new, frame, atol=1e-6):
            break
        frame = new
    residual = pos - star_mean[:, None] - frame[None]
    good = np.isfinite(residual) & np.isfinite(err)
    floor = excess_scatter(residual[good], err[good])
    result['frame_offsets_arcsec'] = frame.tolist()
    result['astrometric_floor_arcsec'] = floor
    result['frame_method'] = 'Alternating robust medians over field stars measured in >=2 epochs; reference epoch frame fixed at zero.'
    measured = [m for m in measurements if m['target'].get('status') == 'MEASURED']
    if len(measured) < 2:
        result['reason'] = 'Target measured in fewer than two epochs.'
        return result
    index = {m['observation_id']: i for i, m in enumerate(measurements)}
    rows = []
    for m in measured:
        t = m['target']; i = index[m['observation_id']]
        x, y = _tangent_arcsec(t['ra_deg'], t['dec_deg'], ra0, dec0)
        cov = np.asarray(t['position_covariance_arcsec2'])
        rows.append({'observation_id': m['observation_id'], 'date': m['date'], 'mjd': m['mjd'],
                     'east_arcsec': float(x - frame[i, 0]), 'north_arcsec': float(y - frame[i, 1]),
                     'sigma_east_arcsec': float(np.sqrt(cov[0, 0] + (floor or 0.)**2)),
                     'sigma_north_arcsec': float(np.sqrt(cov[1, 1] + (floor or 0.)**2))})
    result['positions'] = rows
    # Expected geometric positions from catalog astrometry (with observer parallax when vectors are supplied).
    if target is not None and target.get('pm_ra_cosdec_mas_yr') is not None:
        from astrometry_model import propagate_direction
        expected = []
        for r in rows:
            vector = (observer_vectors or {}).get(r['observation_id'])
            try:
                p = propagate_direction(target, r['date'], observer_position_au=vector['xyz'] if vector else None,
                                        observer_provenance=vector['provenance'] if vector else None)
            except ValueError:
                p = propagate_direction(target, r['date'])
            ex, ey = _tangent_arcsec(p['ra_deg'], p['dec_deg'], ra0, dec0)
            expected.append((float(ex), float(ey), p['coordinate_kind']))
        e = np.array([[x, y] for x, y, _ in expected]); o = np.array([[r['east_arcsec'], r['north_arcsec']] for r in rows])
        s = np.array([[r['sigma_east_arcsec'], r['sigma_north_arcsec']] for r in rows])
        # Absolute zero point is not tied to the catalog frame; compare motion only (weighted mean offset removed).
        w = 1/s**2; shift = np.sum((o-e)*w, axis=0)/np.sum(w, axis=0)
        residual = o - e - shift
        chi = float(np.sum((residual/s)**2))
        from scipy.stats import chi2
        dof = 2*len(rows)-2
        result['catalog_model'] = {'coordinate_kind': expected[0][2], 'constant_offset_arcsec': shift.tolist(),
                                   'residual_arcsec': residual.tolist(), 'chi2': chi, 'dof': dof,
                                   'pvalue': float(chi2.sf(chi, dof)) if dof > 0 else None,
                                   'predicted_displacement_arcsec': float(np.hypot(*(e[-1]-e[0])))}
        for r, (ex, ey, _) in zip(rows, expected):
            r.update(expected_east_arcsec=ex + float(shift[0]), expected_north_arcsec=ey + float(shift[1]))
    if len(rows) >= 3 and len({round(r['mjd'], 3) for r in rows}) >= 3:
        from motion_inference import fit_motion
        motion = fit_motion([r['mjd'] for r in rows], [r['east_arcsec'] for r in rows], [r['north_arcsec'] for r in rows],
                            [r['sigma_east_arcsec'] for r in rows], [r['sigma_north_arcsec'] for r in rows])
        v = np.array(motion['velocity_arcsec_per_year']); cov = np.diag([motion['covariance_x'][1][1], motion['covariance_y'][1][1]])
        motion['speed_arcsec_per_year'] = float(np.hypot(*v))
        motion['velocity_significance'] = float(np.sqrt(v @ np.linalg.solve(cov, v)))
        result['linear_motion_fit'] = motion
        if target is not None and target.get('pm_ra_cosdec_mas_yr') is not None:
            catalog_v = np.array([target['pm_ra_cosdec_mas_yr'], target['pm_dec_mas_yr']])/1000
            result['linear_motion_fit']['catalog_velocity_arcsec_per_year'] = catalog_v.tolist()
            result['linear_motion_fit']['difference_significance'] = float(np.sqrt((v-catalog_v) @ np.linalg.solve(cov, v-catalog_v)))
        result['status'] = 'MOTION_MEASURED' if motion['velocity_significance'] >= policy['detection_sigma'] else 'NO_SIGNIFICANT_MOTION'
        vectors = [((observer_vectors or {}).get(r['observation_id']) or {}).get('xyz') for r in rows]
        if all(v is not None for v in vectors):
            solution = astrometric_solution(rows, vectors, ra0, dec0)
            if solution['status'] == 'SOLVED' and target is not None and target.get('pm_ra_cosdec_mas_yr') is not None:
                pm = np.array([solution['pm_east_arcsec_per_year'], solution['pm_north_arcsec_per_year']])
                catalog = np.array([target['pm_ra_cosdec_mas_yr'], target['pm_dec_mas_yr']])/1000
                cov = np.asarray(solution['covariance'])[2:4, 2:4]
                solution['catalog_pm_arcsec_per_year'] = catalog.tolist()
                solution['catalog_parallax_arcsec'] = (target.get('parallax_mas') or 0)/1000
                solution['pm_significance'] = float(np.sqrt(pm @ np.linalg.solve(cov, pm)))
                solution['pm_difference_from_catalog_sigma'] = float(np.sqrt((pm-catalog) @ np.linalg.solve(cov, pm-catalog)))
                solution['parallax_difference_from_catalog_sigma'] = float(abs(solution['parallax_arcsec']-solution['catalog_parallax_arcsec'])/solution['parallax_sigma_arcsec'])
            result['astrometric_solution'] = solution
            if solution['status'] == 'SOLVED':
                result['status'] = 'MOTION_MEASURED' if solution.get('pm_significance', 0) >= policy['detection_sigma'] else 'NO_SIGNIFICANT_MOTION'
    else:
        result['status'] = 'TWO_EPOCH_DISPLACEMENT_ONLY'
        result['reason'] = 'A trajectory fit needs three distinct epochs.'
    return result


def astrometric_solution(rows, observers, ra0, dec0):
    """Weighted 5-parameter fit: offset (2), proper motion (2), parallax (1).

    Parallax factors use the observer's barycentric ICRF position (AU), e.g. JPL
    Horizons vectors for SPHEREx:
      d_east  = plx * (X sin a - Y cos a)
      d_north = plx * (X cos a sin d + Y sin a sin d - Z cos d)
    Covariance is inflated by the Birge ratio sqrt(chi2/dof) when it exceeds 1.
    """
    a, d = np.deg2rad(ra0), np.deg2rad(dec0)
    t_ref = float(np.mean([r['mjd'] for r in rows]))
    design, values, sigma = [], [], []
    for r, xyz in zip(rows, observers):
        X, Y, Z = xyz; dt = (r['mjd']-t_ref)/365.25
        pe = X*np.sin(a) - Y*np.cos(a)
        pn = X*np.cos(a)*np.sin(d) + Y*np.sin(a)*np.sin(d) - Z*np.cos(d)
        design.append([1, 0, dt, 0, pe]); values.append(r['east_arcsec']); sigma.append(r['sigma_east_arcsec'])
        design.append([0, 1, 0, dt, pn]); values.append(r['north_arcsec']); sigma.append(r['sigma_north_arcsec'])
    X = np.asarray(design)/np.asarray(sigma)[:, None]; y = np.asarray(values)/np.asarray(sigma)
    dof = len(values)-5
    if dof < 1 or np.linalg.matrix_rank(X) < 5:
        return {'status': 'UNDERCONSTRAINED', 'reason': 'Need at least three epochs with non-degenerate parallax factors.'}
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    covariance = np.linalg.inv(X.T @ X)
    chi2 = float(np.sum((y - X @ beta)**2)); birge = max(1., float(np.sqrt(chi2/dof)))
    covariance *= birge**2; err = np.sqrt(np.diag(covariance))
    return {'status': 'SOLVED', 'reference_mjd': t_ref, 'pm_east_arcsec_per_year': float(beta[2]), 'pm_north_arcsec_per_year': float(beta[3]),
            'pm_sigma_arcsec_per_year': [float(err[2]), float(err[3])], 'parallax_arcsec': float(beta[4]), 'parallax_sigma_arcsec': float(err[4]),
            'chi2': chi2, 'dof': dof, 'birge_ratio_applied': birge, 'covariance': covariance.tolist(),
            'model': 'Linear proper motion plus annual parallax from observer barycentric vectors; tangent-plane offsets.'}


def render(report, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.4))
    lc = report['photometry']; ax = axes[0]
    used = [p for p in lc.get('points', []) if p['status'] == 'USED']
    if used:
        ax.errorbar([p['mjd'] for p in used], [p['calibrated_flux_jy']*1e3 for p in used], yerr=[p['flux_sigma_jy']*1e3 for p in used], fmt='o', color='#2f6fde')
        if 'constant_flux_test' in lc:
            ax.axhline(lc['constant_flux_test']['constant_flux']*1e3, color='#888', lw=1, ls='--')
        ax.set(xlabel='MJD', ylabel='Flux density (mJy)', title=f"Same-wavelength light curve: {lc['status']}")
    else:
        ax.text(.5, .5, lc.get('reason', lc['status']), ha='center', va='center', wrap=True, transform=ax.transAxes); ax.set_axis_off()
    ast = report['astrometry']; ax = axes[1]
    rows = ast.get('positions', [])
    if rows:
        ax.errorbar([r['east_arcsec'] for r in rows], [r['north_arcsec'] for r in rows], xerr=[r['sigma_east_arcsec'] for r in rows],
                    yerr=[r['sigma_north_arcsec'] for r in rows], fmt='o', color='#d9480f', label='measured (frame-corrected)')
        if 'expected_east_arcsec' in rows[0]:
            ax.plot([r['expected_east_arcsec'] for r in rows], [r['expected_north_arcsec'] for r in rows], 's--', mfc='none', color='#333', label='catalog prediction')
        for r in rows:
            ax.annotate(r['date'][:10], (r['east_arcsec'], r['north_arcsec']), fontsize=7, xytext=(4, 4), textcoords='offset points')
        ax.invert_xaxis(); ax.set(xlabel='East offset (arcsec)', ylabel='North offset (arcsec)', title=f"Astrometric trajectory: {ast['status']}"); ax.legend(fontsize=8)
        ax.set_aspect('equal', adjustable='datalim')
    else:
        ax.text(.5, .5, ast.get('reason', ast['status']), ha='center', va='center', wrap=True, transform=ax.transAxes); ax.set_axis_off()
    fig.tight_layout(); fig.savefig(path, dpi=110); plt.close(fig)


def run_series(observation_ids, records, *, target=None, position=None, center, output, catalog=None, catalog_provenance=None,
               use_horizons=True, progress=None):
    """Download/reuse real cutouts for a matched group and write series products."""
    from difference_imaging import load_epoch
    from epoch_lab import CALIBRATION_CACHE, CUTOUT_SIZE_DEG, SCIENCE_CACHE, ARTIFACTS, science_files, software_sha256
    from scientific_pipeline import download_cutout

    def emit(stage, **fields):
        if progress:
            progress({'stage': stage, **fields})
    started = time.perf_counter()
    index = {r['obs_publisher_did']: r for r in records}
    missing = [i for i in observation_ids if i not in index]
    if missing:
        raise ValueError('Observation not in the discovered archive records: ' + ', '.join(missing))
    epochs = []
    unavailable = []
    for ident in observation_ids:
        try:
            emit('download', observation_id=ident)
            path = download_cutout(index[ident], center[0], center[1], SCIENCE_CACHE, size=CUTOUT_SIZE_DEG)
            emit('calibration', observation_id=ident)
            epochs.append(load_epoch(path, ident, CALIBRATION_CACHE, product_identity=index[ident].get('access_url')))
        except Exception as exc:
            # One bad epoch must not discard independent, fully calibrated data.
            # Preserve the failure and continue; never replace its calibration.
            unavailable.append({'observation_id': ident, 'status': 'UNAVAILABLE',
                                'stage': 'download_or_calibration',
                                'error_type': type(exc).__name__, 'reason': str(exc)[:500]})
            emit('observation_unavailable', observation_id=ident,
                 reason=f'{type(exc).__name__}: {str(exc)[:300]}')
    vectors = {}
    if use_horizons and target is not None and target.get('parallax_mas'):
        emit('observer_ephemeris')
        from spacecraft_ephemeris import observer_vector
        for e in epochs:
            try:
                xyz, provenance = observer_vector(e['date'], ARTIFACTS/'observer_cache')
                vectors[e['id']] = {'xyz': xyz, 'provenance': provenance}
            except Exception as exc:  # recorded; model falls back to barycentric direction
                vectors[e['id']] = None
                emit('observer_ephemeris_unavailable', reason=str(exc))
    emit('measurement')
    if len(epochs) >= 2:
        report = analyze_series(epochs, target=target, position=position, catalog_sources=catalog, observer_vectors={k: v for k, v in vectors.items() if v})
    else:
        reason = f'Only {len(epochs)} of {len(observation_ids)} requested observations could be downloaded and calibrated; at least two are required for a time series.'
        report = {'policy': POLICY, 'epoch_count': len(epochs),
                  'detector': int(observation_ids[0][-1]), 'observations': [e['id'] for e in epochs],
                  'field_stars_detected': 0,
                  'photometry': {'status': 'INSUFFICIENT_EVIDENCE', 'points': [], 'reason': reason},
                  'astrometry': {'status': 'INSUFFICIENT_EVIDENCE', 'positions': [], 'reason': reason},
                  'epochs': [], 'limitations': [reason]}
    report['analysis_status'] = 'PARTIAL' if unavailable else 'COMPLETE'
    report['requested_epoch_count'] = len(observation_ids)
    report['requested_observations'] = list(observation_ids)
    report['unavailable_observations'] = unavailable
    if unavailable:
        report['limitations'].append('Unavailable observations are omitted from the fitted time series and listed with their processing errors; no calibration or measurement was inferred for them.')
    report['catalog'] = catalog_provenance
    report['science_files'] = science_files(epochs, index)
    report['observer_ephemeris'] = {k: (v['provenance'] if v else 'UNAVAILABLE') for k, v in vectors.items()} if vectors else 'NOT_USED'
    report['software_sha256'] = software_sha256()
    report['elapsed_seconds'] = time.perf_counter()-started
    out = Path(output); out.mkdir(parents=True, exist_ok=True)
    render(report, out/'series.png')
    rows = ['observation_id,date,mjd,status,wavelength_um,calibrated_flux_jy,flux_sigma_jy']
    for p in report['photometry'].get('points', []):
        rows.append(','.join(str(p.get(k, '')) for k in ('observation_id', 'date', 'mjd', 'status', 'wavelength_um', 'calibrated_flux_jy', 'flux_sigma_jy')))
    for failed in unavailable:
        rows.append(','.join([failed['observation_id'], '', '', 'UNAVAILABLE', '', '', '']))
    (out/'light_curve.csv').write_text('\n'.join(rows)+'\n', encoding='utf-8')
    report['files'] = ['series.json', 'series.png', 'light_curve.csv']
    (out/'series.json').write_text(json.dumps(report, indent=2, allow_nan=False), encoding='utf-8')
    return report
