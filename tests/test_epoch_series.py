"""Same-wavelength time series: ePSF fitting, light curves and the astrometric solution."""
import numpy as np
import pytest

import synthetic_sky as syn
from spherex_lab.lab.difference_imaging import fit_background
from spherex_lab.lab.epoch_series import analyze_series, astrometric_solution, psf_fit

RA, DEC = 150.0, 20.0


def test_psf_fit_recovers_subpixel_position_and_flux():
    rows = [{'source_id': '1', 'ra': RA, 'dec': DEC, 'pmra': 0., 'pmdec': 0., 'ref_epoch': 2016.0, 'phot_g_mean_mag': 12., 'flux': 400.}]
    e = syn.epoch(syn.stars_at(rows, '2025-09-01'), ra=RA, dec=DEC, date='2025-09-01', seed=3, invalid_fraction=0.)
    s = e['science']
    bg = fit_background(s['image'], s['variance'], np.ones_like(s['valid']))
    x0, y0 = e['wcs'].world_to_pixel_values(RA, DEC)
    fit = psf_fit(e, bg['residual'], bg['residual_variance'], s['valid'], float(x0)+0.7, float(y0)-0.6)
    assert fit is not None
    assert fit['x'] == pytest.approx(float(x0), abs=0.05) and fit['y'] == pytest.approx(float(y0), abs=0.05)
    sigma = np.sqrt(np.diag(fit['position_covariance_px2']))
    assert np.all(sigma < 0.05)


def _series(n=4, flux_scale=(1., 1., 1., 1.), pm_mas=0., seed=7):
    rows = syn.catalog(RA, DEC, 120, seed=21)
    target = {'id': 'syn', 'ra_deg': RA, 'dec_deg': DEC, 'pm_ra_cosdec_mas_yr': pm_mas, 'pm_dec_mas_yr': 0., 'epoch': 'J2000', 'parallax_mas': None}
    dates = ['2025-05-01T00:00:00', '2025-08-01T00:00:00', '2025-11-01T00:00:00', '2026-02-01T00:00:00'][:n]
    epochs = []
    for k, date in enumerate(dates):
        from spherex_lab.lab.difference_imaging import propagate_catalog
        tp = propagate_catalog([{'ra': RA, 'dec': DEC, 'pmra': pm_mas, 'pmdec': 0., 'ref_epoch': 2000.0, 'source_id': 't'}], date)
        stars = syn.stars_at(rows, date) + [(float(tp['ra'][0]), float(tp['dec'][0]), 300.*flux_scale[k])]
        epochs.append(syn.epoch(stars, ra=RA, dec=DEC, date=date, seed=seed+k, shift_px=(0.3*k, -0.2*k), invalid_fraction=0.,
                                ident=f'ivo://irsa.ipac/spherex_dr1?S{k}/D1'))
    return epochs, target


def test_constant_target_is_consistent_and_variable_target_is_flagged():
    epochs, target = _series()
    report = analyze_series(epochs, target=target)
    assert report['photometry']['status'] == 'CONSISTENT_WITH_CONSTANT'
    assert report['photometry']['used_epochs'] == 4
    epochs, target = _series(flux_scale=(1., 1., 1.6, 1.))
    report = analyze_series(epochs, target=target)
    assert report['photometry']['status'] == 'VARIABILITY_CANDIDATE'


def test_known_proper_motion_is_measured_from_images():
    epochs, target = _series(pm_mas=9000.)   # 9"/yr ~ 1.1 px over 9 months
    report = analyze_series(epochs, target=target)
    motion = report['astrometry']['linear_motion_fit']
    assert report['astrometry']['status'] == 'MOTION_MEASURED'
    assert motion['velocity_arcsec_per_year'][0] == pytest.approx(9.0, abs=0.6)
    assert abs(motion['velocity_arcsec_per_year'][1]) < 0.6
    assert motion['difference_significance'] < 4


def test_astrometric_solution_recovers_motion_and_parallax():
    rng = np.random.default_rng(4)
    mjd = np.array([60790., 60833., 60952., 60958., 61155., 61200.])
    angles = 2*np.pi*(mjd-60800.)/365.25
    observers = [np.array([np.cos(a), np.sin(a)*0.917, np.sin(a)*0.398]) for a in angles]
    ra0, dec0 = 331.0, -56.8
    a, d = np.deg2rad(ra0), np.deg2rad(dec0)
    true = {'pm': (4.0, -2.5), 'plx': 0.27}
    rows = []
    for t, (X, Y, Z) in zip(mjd, observers):
        dt = (t-mjd.mean())/365.25
        pe = X*np.sin(a) - Y*np.cos(a); pn = X*np.cos(a)*np.sin(d) + Y*np.sin(a)*np.sin(d) - Z*np.cos(d)
        rows.append({'mjd': t, 'east_arcsec': 1.0 + true['pm'][0]*dt + true['plx']*pe + rng.normal(0, 0.02),
                     'north_arcsec': -2.0 + true['pm'][1]*dt + true['plx']*pn + rng.normal(0, 0.02),
                     'sigma_east_arcsec': 0.02, 'sigma_north_arcsec': 0.02})
    sol = astrometric_solution(rows, observers, ra0, dec0)
    assert sol['status'] == 'SOLVED'
    assert sol['pm_east_arcsec_per_year'] == pytest.approx(4.0, abs=0.1)
    assert sol['pm_north_arcsec_per_year'] == pytest.approx(-2.5, abs=0.1)
    assert sol['parallax_arcsec'] == pytest.approx(0.27, abs=0.06)


def test_series_rejects_mixed_detectors():
    epochs, target = _series(n=2)
    epochs[1]['detector'] = 2
    with pytest.raises(ValueError, match='one detector'):
        analyze_series(epochs, target=target)
