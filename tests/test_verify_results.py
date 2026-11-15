import json
import shutil
from pathlib import Path

import pytest

from scripts import verify_results
DEMOS = Path('artifacts/epoch_lab/demos')


def failures(result):
    return [(r['run'], r['check'], r['detail']) for r in result['checks'] if r['status'] == 'FAIL']


def test_recorded_runs_pass_independent_checks():
    result = verify_results.verify(DEMOS)
    assert result['runs'], 'no recorded runs found'
    stale = [f for f in failures(result) if f[1] == 'produced by the current processing code']
    assert not stale, f'recorded runs predate the processing code; regenerate them with python -m scripts.make_demos: {stale}'
    assert not failures(result)


def _copy(tmp_path, name):
    shutil.copytree(DEMOS / name, tmp_path / name)
    return tmp_path / name


def _edit(path, change):
    data = json.loads(path.read_text(encoding='utf-8'))
    change(data)
    path.write_text(json.dumps(data), encoding='utf-8')


@pytest.mark.parametrize('change, check', [
    (lambda r: r['spectral'].update(matched_pixels=r['spectral']['matched_pixels'] + 1), 'compared pixel count'),
    (lambda r: r['measured_stars'][0].update(change_z=r['measured_stars'][0]['change_z'] + 0.5), 'star σ = √(stat² + β² + (αF)²) and z = residual / σ'),
    (lambda r: r.update(status='VALIDATED', gates={**r['gates'], 'noise_model': 'UNVERIFIED'}), 'status follows the gates (UNVERIFIED never counts as PASS)'),
    (lambda r: r.update(separation_days=r['separation_days'] + 1), 'time baseline'),
])
def test_corrupted_comparison_is_caught(tmp_path, change, check):
    directory = _copy(tmp_path, 'barnard-pair')
    _edit(directory / 'report.json', change)
    assert check in {f[1] for f in failures(verify_results.verify(tmp_path))}


@pytest.mark.parametrize('change, check', [
    (lambda r: r['astrometry']['astrometric_solution'].update(pm_significance=20.0), 'proper-motion detection significance'),
    (lambda r: r['astrometry']['astrometric_solution'].update(pm_east_arcsec_per_year=6.0), 'measured motion agrees with SIMBAD (known-mover test)'),
    (lambda r: r['photometry']['constant_flux_test'].update(pvalue=0.9), 'constant-flux χ² and p-value'),
])
def test_corrupted_series_is_caught(tmp_path, change, check):
    directory = _copy(tmp_path, 'eps-indi-b-series')
    _edit(directory / 'series.json', change)
    assert check in {f[1] for f in failures(verify_results.verify(tmp_path))}


def test_wrong_pair_is_caught(tmp_path):
    directory = _copy(tmp_path, 'screening-five-targets')
    path = directory / 'screening_barnards-star.json'
    data = json.loads(path.read_text(encoding='utf-8'))
    other = next(o for o in data['observations'] if o['status'] == 'SCREENED' and o['detector'] != data['matched_pairs'][0]['detector'])
    data['matched_pairs'][0]['observation_b'] = other['observation_id']
    path.write_text(json.dumps(data), encoding='utf-8')
    found = failures(verify_results.verify(tmp_path))
    assert ('screening-five-targets / barnards-star', 'every pair: same detector, ≥1 day apart, |Δλ| ≤ 0.1 × CBAND') in {f[:2] for f in found}


def test_software_hash_ignores_line_endings(tmp_path, monkeypatch):
    from spherex_lab.lab import epoch_lab
    for name in epoch_lab.SOFTWARE_FILES:
        (tmp_path / name).parent.mkdir(parents=True, exist_ok=True)
        (tmp_path / name).write_bytes(b'line one\nline two\n')
    monkeypatch.setattr(epoch_lab, 'ROOT', tmp_path)
    unix = epoch_lab.software_sha256()
    for name in epoch_lab.SOFTWARE_FILES:
        (tmp_path / name).write_bytes(b'line one\r\nline two\r\n')
    assert epoch_lab.software_sha256() == unix
    (tmp_path / epoch_lab.SOFTWARE_FILES[0]).write_bytes(b'line one\r\nline 2\r\n')
    assert epoch_lab.software_sha256() != unix
