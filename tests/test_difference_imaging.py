"""Matched-wavelength difference imaging: exact-truth synthetic tests and a real-data check."""
import json
from pathlib import Path

import numpy as np
import pytest

import synthetic_sky as syn
from spherex_lab.lab.difference_imaging import POLICY, compare, convolve_masked, excess_scatter, propagate_catalog, write_products

RA, DEC = 150.0, 20.0
DATE_A, DATE_B = '2025-08-30T00:00:00', '2025-09-30T00:00:00'


def _pair(rows_a=None, rows_b=None, *, wave_offset=0., seed=1, shift=(0.37, -0.21), detector_b=1, extra_b=()):
    rows = syn.catalog(RA, DEC, 140, seed=11)
    rows_a = rows if rows_a is None else rows_a
    rows_b = rows if rows_b is None else rows_b
    a = syn.epoch(syn.stars_at(rows_a, DATE_A), ra=RA, dec=DEC, date=DATE_A, seed=seed)
    b = syn.epoch(syn.stars_at(rows_b, DATE_B) + list(extra_b), ra=RA, dec=DEC, date=DATE_B, seed=seed+100,
                  shift_px=shift, wave_offset=wave_offset, detector=detector_b,
                  ident='ivo://irsa.ipac/spherex_dr1?SYNB/D%d' % detector_b)
    return a, b, rows


def _blank_position(rows, a, margin=12, clearance=6.0):
    cat = propagate_catalog(rows, DATE_A)
    x, y = a['wcs'].world_to_pixel_values(cat['ra'], cat['dec'])
    for py in range(margin, 120-margin, 3):
        for px in range(margin, 120-margin, 3):
            if np.min(np.hypot(x-px, y-py)) > clearance:
                return px, py
    raise AssertionError('no blank position')


def test_no_change_pair_is_validated_and_quiet():
    a, b, rows = _pair()
    report = compare(a, b, catalog_sources=rows)['report']
    assert report['gates']['spectral_overlap'] == 'PASS'
    assert report['gates']['registration'] == 'PASS'
    assert abs(report['photometry']['scale_applied'] - 1) < 0.02
    assert 0.8 < report['null_test']['robust_std_sigma'] < 1.6
    interesting = {'CANDIDATE_VARIABLE', 'UNCATALOGED_VARIABLE_CANDIDATE', 'TRANSIENT_CANDIDATE', 'UNASSOCIATED_CANDIDATE'}
    assert not [c for c in report['candidates'] if c['classification'] in interesting]
    assert report['variability_summary']['passed_vetting'] == 0


def test_injected_transient_is_found_where_it_was_put():
    a, _, rows = _pair()
    px, py = _blank_position(rows, a)
    ra, dec = a['wcs'].pixel_to_world_values(px, py)
    a, b, rows = _pair(extra_b=[(float(ra), float(dec), 60.0)])
    report = compare(a, b, catalog_sources=rows)['report']
    found = [c for c in report['candidates'] if c['classification'] == 'TRANSIENT_CANDIDATE']
    assert found, report['candidate_summary']
    best = min(found, key=lambda c: np.hypot(c['x']-px, c['y']-py))
    assert np.hypot(best['x']-px, best['y']-py) < 1.0
    assert best['sign'] == 'brightening' and best['present_in'] == 'B'


def test_brightened_catalog_star_is_a_vetted_variability_candidate():
    rows = syn.catalog(RA, DEC, 140, seed=11)
    a = syn.epoch(syn.stars_at(rows, DATE_A), ra=RA, dec=DEC, date=DATE_A, seed=1)
    cat = propagate_catalog(rows, DATE_A)
    x, y = a['wcs'].world_to_pixel_values(cat['ra'], cat['dec'])
    # an isolated, moderately bright star well inside the field
    ok = [i for i in range(len(rows)) if 15 < x[i] < 105 and 15 < y[i] < 105 and np.sort(np.hypot(x-x[i], y-y[i]))[1] > 8 and rows[i]['flux'] > 150]
    assert ok
    changed = [dict(r) for r in rows]; target = ok[0]; changed[target]['flux'] *= 1.5
    a, b, _ = _pair(rows_a=rows, rows_b=changed)
    report = compare(a, b, catalog_sources=rows)['report']
    passed = [v for v in report['variability_candidates'] if v['classification'] == 'SOURCE_LEVEL_VARIABILITY_CANDIDATE']
    assert [v['source_id'] for v in passed] == [rows[target]['source_id']]
    assert passed[0]['fractional_change'] == pytest.approx(0.5, abs=0.12)


def test_different_wavelengths_are_rejected_not_compared():
    a, b, rows = _pair(wave_offset=0.01)   # half a bandwidth everywhere
    report = compare(a, b, catalog_sources=rows)['report']
    assert report['status'] == 'REJECTED'
    assert report['gates']['spectral_overlap'] == 'FAIL'
    assert report['spectral']['matched_pixels'] == 0


def test_different_detectors_raise():
    a, b, rows = _pair(detector_b=2)
    with pytest.raises(ValueError, match='different detectors'):
        compare(a, b, catalog_sources=rows)


def test_known_proper_motion_is_explained_not_flagged():
    rows = syn.catalog(RA, DEC, 140, seed=11)
    a = syn.epoch(syn.stars_at(rows, DATE_A), ra=RA, dec=DEC, date=DATE_A, seed=1)
    cat = propagate_catalog(rows, DATE_A)
    x, y = a['wcs'].world_to_pixel_values(cat['ra'], cat['dec'])
    i = int(np.argmax([r['flux'] if 20 < x[k] < 100 and 20 < y[k] < 100 else 0 for k, r in enumerate(rows)]))
    moving = [dict(r) for r in rows]
    moving[i].update(pmra=60000.0, ref_epoch=2025.66)   # 60"/yr from the epoch-A date -> ~0.8 px in a month
    a = syn.epoch(syn.stars_at(moving, DATE_A), ra=RA, dec=DEC, date=DATE_A, seed=1)
    b = syn.epoch(syn.stars_at(moving, DATE_B), ra=RA, dec=DEC, date=DATE_B, seed=101, shift_px=(0.37, -0.21), ident='ivo://irsa.ipac/spherex_dr1?SYNB/D1')
    report = compare(a, b, catalog_sources=moving)['report']
    near = [c for c in report['candidates'] if c.get('catalog_source_id') == moving[i]['source_id']]
    assert near and all(c['classification'] == 'KNOWN_PROPER_MOTION_RESIDUAL' for c in near)


def test_products_are_written_and_json_is_strict(tmp_path):
    a, b, rows = _pair()
    result = compare(a, b, catalog_sources=rows)
    files = write_products(result, tmp_path)
    assert set(files) <= {p.name for p in tmp_path.iterdir()}
    report = json.loads((tmp_path/'report.json').read_text(encoding='utf-8'))
    viewer = json.loads((tmp_path/'viewer.json').read_text(encoding='utf-8'))
    assert report['status'] == result['report']['status']
    assert viewer['shape'] == [120, 120] and len(viewer['planes']['difference']) > 1000
    from astropy.io import fits
    with fits.open(tmp_path/'comparison.fits', checksum=True) as hdus:
        assert hdus['DIFF'].header['BUNIT'] == 'MJy/sr' and hdus['MATCHED'].data.shape == (120, 120)


def test_convolve_masked_preserves_constants_and_tolerates_isolated_bad_pixels():
    image = np.full((30, 30), 2.0); variance = np.full((30, 30), 0.5); valid = np.ones((30, 30), bool)
    valid[15, 15] = False
    kernel = syn.gaussian_epsf(1, 7, 1.0)
    out = convolve_masked(image, variance, valid, kernel, 0.95)
    assert np.allclose(out['image'][out['valid']], 2.0)
    assert not out['valid'][15, 15]                        # the bad pixel itself carries <95% weight
    assert out['valid'][15, 18] and out['valid'][12, 15]   # 3 px away survives; exact-support convolution would drop it
    assert out['valid'].sum() > 30*30 - 9 - 4*(30*3)       # only the bad pixel's core and the borders are lost
    inner = out['variance'][10, 10]
    assert inner == pytest.approx(0.5*np.sum((kernel/kernel.sum())**2), rel=1e-9)


def test_excess_scatter_recovers_a_known_systematic():
    rng = np.random.default_rng(3)
    stat = rng.uniform(0.01, 0.05, 4000)
    values = rng.normal(0, np.sqrt(stat**2 + 0.03**2))
    assert excess_scatter(values, stat) == pytest.approx(0.03, rel=0.1)
    assert excess_scatter(rng.normal(0, stat), stat) < 0.006


def test_vectorised_propagation_matches_reference_catalog_oracle():
    from spherex_lab.data.reference_catalog import source_position
    rows = syn.catalog(RA, DEC, 25, seed=5, pm_mas=0.)
    rng = np.random.default_rng(9)
    for r in rows:
        r['pmra'], r['pmdec'] = float(rng.normal(0, 500)), float(rng.normal(0, 500))
    p = propagate_catalog(rows, '2026-03-04T07:58:55')
    for k, r in enumerate(rows):
        ref = source_position(r, '2026-03-04T07:58:55')
        assert p['ra'][k] == pytest.approx(ref['ra_deg'], abs=1e-10)
        assert p['dec'][k] == pytest.approx(ref['dec_deg'], abs=1e-10)


def _cached_real_pair():
    demo = Path('artifacts/epoch_lab/demos/barnard-pair/report.json')
    if not demo.is_file():
        return None
    report = json.loads(demo.read_text(encoding='utf-8'))
    files = report.get('science_files', {})
    if len(files) != 2 or not all(Path(v['path']).is_file() for v in files.values()):
        return None
    return report


@pytest.mark.skipif(_cached_real_pair() is None, reason='Barnard pair demo and its cached SPHEREx cutouts/calibrations are not present locally; run: python -m scripts.make_demos barnard-pair')
def test_real_barnard_pair_reproduces_recorded_validation():
    """Re-run the recorded real comparison from cached archive files and check its gates."""
    from spherex_lab.lab.difference_imaging import load_epoch
    from spherex_lab.lab.epoch_lab import CALIBRATION_CACHE, field_catalog, field_center
    recorded = _cached_real_pair()
    epochs = [load_epoch(v['path'], k, CALIBRATION_CACHE, product_identity=v.get('access_url')) for k, v in recorded['science_files'].items()]
    catalog, _ = field_catalog(*field_center('barnards-star'), allow_network=False)
    report = compare(*epochs, catalog_sources=catalog)['report']
    assert report['status'] == recorded['status'] == 'VALIDATED'
    assert report['spectral']['matched_pixels'] == recorded['spectral']['matched_pixels']
    assert np.hypot(*report['registration']['median_offset_px']) < POLICY['registration_tolerance_px']
    assert 0.7 <= report['systematic_floor']['calibration_check']['all']['robust_std_z'] <= 1.4


def test_without_a_catalog_image_detected_sources_calibrate_the_gates():
    a, b, rows = _pair()
    report = compare(a, b, catalog_sources=None)['report']
    assert report['calibration_sources'].startswith('image-detected')
    assert report['gates']['registration'] == 'PASS' and report['gates']['photometric_consistency'] == 'PASS'
    assert report['systematic_floor']['stars'] >= 30
    assert all(c['classification'] != 'UNASSOCIATED_CANDIDATE' for c in report['candidates'])   # never claimed without a catalog
