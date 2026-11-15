"""Synthetic SPHEREx-like epochs with known truth, for tests and method validation only.

Builds dictionaries with the same structure as ``difference_imaging.load_epoch``
(native image/variance/flags, celestial WCS, prepared ePSF table, wavelength maps)
so that the production comparison and time-series code can be checked against an
exact answer: injected transients, flux changes, proper motions and noise levels.
Nothing here is used by the archive pipeline.
"""
import numpy as np
from astropy.time import Time
from astropy.wcs import WCS

PIXEL_ARCSEC = 6.15
SIGMA_PX = 0.6


def make_wcs(ra, dec, shape, *, shift_px=(0., 0.), rotation_deg=0.):
    w = WCS(naxis=2)
    w.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    w.wcs.crval = [ra, dec]
    w.wcs.crpix = [shape[1]/2 + 0.5 + shift_px[0], shape[0]/2 + 0.5 + shift_px[1]]
    scale = PIXEL_ARCSEC/3600
    t = np.deg2rad(rotation_deg)
    w.wcs.cd = scale*np.array([[-np.cos(t), np.sin(t)], [np.sin(t), np.cos(t)]])
    return w


def gaussian_epsf(oversampling=5, size=33, sigma_px=SIGMA_PX):
    yy, xx = np.indices((size, size)) - size//2
    psf = np.exp(-(xx**2+yy**2)/(2*(sigma_px*oversampling)**2))
    return psf/psf.sum()


def prepared_epsf(shape, detector=1, oversampling=5):
    raw = gaussian_epsf(oversampling)
    return {'shape': shape, 'x0': 900, 'y0': 900, 'detector': detector,
            'rows': {'XCENTER': np.array([900.+shape[1]/2]), 'YCENTER': np.array([900.+shape[0]/2]), 'EPSF': raw[None],
                     'NSTAR': np.array([50]), 'CWAVE': np.array([0.85]), 'CBAND': np.array([0.02])},
            'oversampling': (oversampling, oversampling), 'history': (), 'zones': {}}


def render(stars, wcs, shape, sigma_px=SIGMA_PX):
    """Point sources sampled at pixel centres (the ePSF convention), fluxes in image units x pixels."""
    image = np.zeros(shape)
    yy, xx = np.indices(shape)
    for ra, dec, flux in stars:
        x, y = wcs.world_to_pixel_values(ra, dec)
        if not (-5 < x < shape[1]+5 and -5 < y < shape[0]+5):
            continue
        profile = np.exp(-((xx-x)**2+(yy-y)**2)/(2*sigma_px**2))
        image += flux*profile/(2*np.pi*sigma_px**2)
    return image


def catalog(ra, dec, n, seed, half_width_deg=0.09, pm_mas=0.):
    rng = np.random.default_rng(seed)
    rows = []
    for i in range(n):
        rows.append({'source_id': str(1000+i), 'ra': ra + rng.uniform(-half_width_deg, half_width_deg)/np.cos(np.deg2rad(dec)),
                     'dec': dec + rng.uniform(-half_width_deg, half_width_deg), 'pmra': pm_mas, 'pmdec': 0., 'ref_epoch': 2016.0,
                     'phot_g_mean_mag': float(rng.uniform(12, 18)), 'flux': float(10**rng.uniform(1.0, 2.6))})
    return rows


def epoch(stars, *, ra, dec, date, shape=(120, 120), seed=0, shift_px=(0., 0.), rotation_deg=0., wave_offset=0.,
          sky=0.3, read_var=0.01, gain=50., detector=1, invalid_fraction=0.005, ident='ivo://irsa.ipac/spherex_dr1?SYN/D1'):
    """One synthetic epoch. ``stars``: (ra, dec, flux) at this epoch. Returns a load_epoch-like dict."""
    rng = np.random.default_rng(seed)
    w = make_wcs(ra, dec, shape, shift_px=shift_px, rotation_deg=rotation_deg)
    yy, xx = np.indices(shape)
    truth = sky + 0.002*(xx-shape[1]/2) + render(stars, w, shape)
    variance = read_var + np.clip(truth, 0, None)/gain
    image = truth + rng.normal(0, np.sqrt(variance))
    flags = np.zeros(shape, np.int32)
    bad = rng.random(shape) < invalid_fraction
    flags[bad] = 1 << 0
    valid = flags == 0
    # Linear-variable-filter wavelength defined on the sky so equal-orientation epochs match.
    sx, _ = make_wcs(ra, dec, shape).world_to_pixel_values(*w.pixel_to_world_values(xx, yy))
    wave = 0.85 + 2e-4*sx + wave_offset
    band = np.full(shape, 0.02)
    return {'id': ident, 'path': '<synthetic>', 'wcs': w,
            'science': {'image': image, 'variance': variance, 'flags': flags, 'valid': valid, 'wcs': w, 'header': {}},
            'epsf': prepared_epsf(shape, detector), 'wave': wave, 'band': band, 'calibration_valid': np.ones(shape, bool),
            'calibration': {'calibration_url': 'synthetic'}, 'zodi': None, 'date': date, 'mjd': float(Time(date).mjd),
            'detector': detector, 'psf_fwhm_arcsec': 2.355*SIGMA_PX*PIXEL_ARCSEC}


def stars_at(rows, date):
    """(ra, dec, flux) at ``date`` from catalog rows (rectilinear PM, as in the pipeline)."""
    from difference_imaging import propagate_catalog
    p = propagate_catalog(rows, date)
    return [(p['ra'][i], p['dec'][i], rows[i]['flux']) for i in range(len(rows))]
