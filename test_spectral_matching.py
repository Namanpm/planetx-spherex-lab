"""Spectral-epoch screening: footprint model against real cutouts, pairing and grouping rules."""
import json
from pathlib import Path

import numpy as np
import pytest

from spectral_matching import (POLICY, apply_transform, detector_maps, footprint_pixel, footprint_transform, matched_groups,
                               matched_pairs, predicted_matched_fraction, screen_observations, validate_footprint_model)


def _square(ra, dec, half=1.745):   # 2040 px x 6.15" = 3.49 deg detector, square on the tangent plane
    from spectral_matching import field_grid
    gra, gdec = field_grid(ra, dec, half, n=2)
    corners = [(gra[0, 0], gdec[0, 0]), (gra[0, 1], gdec[0, 1]), (gra[1, 1], gdec[1, 1]), (gra[1, 0], gdec[1, 0])]
    return 'POLYGON ICRS ' + ' '.join(f'{a:.8f} {d:.8f}' for a, d in corners)


def test_footprint_model_against_real_cutouts():
    cases = json.loads(Path('tests_data/footprint_ground_truth.json').read_text(encoding='utf-8'))['cases']
    assert len(cases) >= 20
    result = validate_footprint_model(cases)
    assert result['median_error_px'] < 4.0
    assert result['max_error_px'] < POLICY['footprint_error_px']


def test_footprint_corners_centre_and_outside():
    region = _square(150.0, 20.0)
    x, y = footprint_pixel(region, 150.0, 20.0)
    assert x == pytest.approx(1019.5, abs=1.0) and y == pytest.approx(1019.5, abs=1.0)
    with pytest.raises(ValueError, match='outside'):
        footprint_pixel(region, 153.0, 20.0)
    with pytest.raises(ValueError, match='four-corner'):
        footprint_pixel('POLYGON ICRS 1 1 2 1 2 2', 1.5, 1.5)


def test_footprint_handles_ra_wrap():
    region = _square(0.01, 10.0)
    x, y = footprint_pixel(region, 359.995, 10.0)
    assert 0 <= x <= 2039 and 0 <= y <= 2039
    xs, ys = apply_transform(footprint_transform(region), np.array([0.01, 359.995]), np.array([10.0, 10.0]))
    assert np.isfinite(xs).all() and xs[1] < xs[0]


def _maps(gradient=1e-4):
    yy, xx = np.indices((2040, 2040))
    return {1: {'wave': (0.8 + gradient*xx).astype(np.float32), 'band': np.full((2040, 2040), 0.02, np.float32), 'url': 'https://example/D1'},
            2: {'wave': (1.2 + gradient*xx).astype(np.float32), 'band': np.full((2040, 2040), 0.03, np.float32), 'url': 'https://example/D2'}}


def _record(ident, mjd, ra, dec, shift_deg=0.0):
    return {'obs_publisher_did': f'ivo://irsa.ipac/spherex_dr1?{ident}', 'obs_id': ident.split('/')[0], 't_min': mjd,
            's_region': _square(ra+shift_deg, dec), 'access_url': 'https://irsa.ipac.caltech.edu/ibe/data/spherex/dr1/level2/x.fits'}


def test_pairs_need_same_detector_distinct_dates_and_matching_wavelength():
    ra, dec = 150.0, 20.0
    records = [_record('A_1/D1', 60900.0, ra, dec), _record('B_1/D1', 60930.0, ra, dec),           # same pointing -> match
               _record('C_1/D1', 60900.3, ra, dec),                                               # same day -> not a pair
               _record('D_1/D1', 60960.0, ra, dec, shift_deg=0.05),                               # moved on detector -> lambda differs
               _record('E_1/D2', 60990.0, ra, dec)]                                               # other detector
    screened = screen_observations(records, ra, dec, _maps())
    assert all(e['status'] == 'SCREENED' for e in screened)
    pairs = matched_pairs(screened, POLICY, maps=_maps(), center=(ra, dec))
    ids = {(p['observation_a'][-6:], p['observation_b'][-6:]) for p in pairs}
    assert ('A_1/D1', 'B_1/D1') in ids and ('C_1/D1', 'B_1/D1') in ids
    assert not any('D_1/D1' in a or 'D_1/D1' in b for a, b in ids)
    assert not any('E_1/D2' in a or 'E_1/D2' in b for a, b in ids)
    assert not any({a, b} == {'A_1/D1', 'C_1/D1'} for a, b in ids)
    best = next(p for p in pairs if p['observation_a'].endswith('A_1/D1') and p['observation_b'].endswith('B_1/D1'))
    assert best['predicted_matched_fraction'] > 0.95 and best['delta_over_bandwidth'] < 1e-6


def test_unscreened_records_are_reported_not_guessed():
    ra, dec = 150.0, 20.0
    maps = _maps(); maps.pop(2)
    records = [_record('A_1/D2', 60900.0, ra, dec), dict(_record('B_1/D1', 60900.0, ra, dec), t_min=None)]
    screened = screen_observations(records, ra, dec, maps)
    reasons = {e['observation_id'][-6:]: e.get('reason') for e in screened}
    assert reasons['A_1/D2'] == 'detector calibration unavailable'
    assert reasons['B_1/D1'] == 'missing observation time'


def test_groups_chain_same_wavelength_epochs():
    ra, dec = 150.0, 20.0
    records = [_record(f'G{k}_1/D1', 60900.0+30*k, ra, dec) for k in range(4)]
    groups = matched_groups(screen_observations(records, ra, dec, _maps()))
    assert groups and groups[0]['epoch_count'] == 4


def test_predicted_fraction_falls_for_a_rotated_revisit():
    ra, dec = 150.0, 20.0
    a = screen_observations([_record('A_1/D1', 60900.0, ra, dec)], ra, dec, _maps(1e-3))[0]
    rotated = _record('B_1/D1', 61080.0, ra, dec)
    verts = np.asarray(rotated['s_region'].split()[2:], float).reshape(-1, 2)
    rotated['s_region'] = 'POLYGON ICRS ' + ' '.join(f'{x:.8f} {y:.8f}' for x, y in np.roll(verts, 2, axis=0))  # 180 deg
    b = screen_observations([rotated], ra, dec, _maps(1e-3))[0]
    area, overlap = predicted_matched_fraction(a, b, _maps(1e-3), (ra, dec), 0.15, POLICY['tolerance_fraction'])
    assert overlap > 0.9 and area < 0.2


def test_detector_maps_records_failures_and_reuses_associations(tmp_path, monkeypatch):
    import spectral_matching as sm
    calls = []

    def fake_discover(pid, detector, client=None):
        calls.append(pid)
        if detector == 2:
            raise ConnectionError('502 from DataLink')
        return 'https://irsa.ipac.caltech.edu/ibe/data/spherex/dr1/spectral_wcs/D1.fits'

    class Cache:
        def arrays(self, url, detector, cache_dir, client=None):
            return np.ones((4, 4)), np.ones((4, 4)), 'sha'
    monkeypatch.setattr(sm, 'discover', fake_discover)
    records = [{'obs_publisher_did': 'ivo://irsa.ipac/spherex_dr1?A_1/D1'}, {'obs_publisher_did': 'ivo://irsa.ipac/spherex_dr1?B_1/D2'}]
    maps = detector_maps(records, tmp_path, cache=Cache())
    assert 'wave' in maps[1] and 'error' in maps[2] and '502' in maps[2]['error'][0]
    calls.clear()
    again = detector_maps(records, tmp_path, cache=Cache())
    assert again[1]['association'] == 'remembered DataLink result'
    assert calls == ['ivo://irsa.ipac/spherex_dr1?B_1/D2']   # only the failed detector is retried
