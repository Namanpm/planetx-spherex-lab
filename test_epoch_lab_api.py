"""HTTP contract of the Epoch Lab API, the single-page front end, and earlier bug fixes."""
import json
import re
import shutil
import subprocess
import threading
import time
import uuid
from pathlib import Path

import numpy as np
import pytest
import requests
from fastapi.testclient import TestClient

import epoch_lab_api as api
import main

PID_A = 'ivo://irsa.ipac/spherex_dr1?2025W35_2A_0474_1/D1'
PID_B = 'ivo://irsa.ipac/spherex_dr1?2025W40_1A_0441_1/D1'


@pytest.fixture
def client(tmp_path, monkeypatch):
    monkeypatch.setattr(api, 'JOB_ROOT', tmp_path/'jobs')
    monkeypatch.setattr(api, 'DEMO_ROOT', tmp_path/'demos')
    with api._lock:
        api._jobs.clear()
    return TestClient(main.app)


def _wait(client, job_id, timeout=10):
    deadline = time.time() + timeout
    while time.time() < deadline:
        state = client.get(f'/api/lab/jobs/{job_id}').json()
        if state['status'] in ('COMPLETED', 'FAILED', 'INTERRUPTED'):
            return state
        time.sleep(0.05)
    raise AssertionError('job did not finish')


@pytest.mark.parametrize('body', [
    {'target_id': 'rw-aurigae', 'ra': 10.0, 'dec': 5.0},
    {'target_id': 'no-such-star'},
    {'ra': 10.0},
    {'ra': 400.0, 'dec': 0.0},
    {'target_id': 'rw-aurigae', 'unexpected': 1},
])
def test_field_validation(client, body):
    assert client.post('/api/lab/screen', json=body).status_code == 422


@pytest.mark.parametrize('a,b', [(PID_A, PID_A), (PID_A, PID_B[:-1]+'2'), ('not-an-id', PID_B)])
def test_compare_validation(client, a, b):
    response = client.post('/api/lab/compare', json={'target_id': 'barnards-star', 'observation_a': a, 'observation_b': b})
    assert response.status_code == 422


def test_series_validation(client):
    assert client.post('/api/lab/series', json={'target_id': 'rw-aurigae', 'observation_ids': [PID_A]}).status_code == 422
    assert client.post('/api/lab/series', json={'target_id': 'rw-aurigae', 'observation_ids': [PID_A, PID_A]}).status_code == 422


def test_job_lifecycle_progress_and_files(client, monkeypatch):
    def fake_execute(kind, body, output, progress):
        progress({'stage': 'comparison'})
        output.mkdir(parents=True, exist_ok=True)
        (output/'report.json').write_text(json.dumps({'status': 'VALIDATED'}), encoding='utf-8')
        return {'status': 'VALIDATED', 'files': ['report.json']}
    monkeypatch.setattr(api, '_execute', fake_execute)
    sub = client.post('/api/lab/compare', json={'target_id': 'barnards-star', 'observation_a': PID_A, 'observation_b': PID_B})
    assert sub.status_code == 202
    state = _wait(client, sub.json()['id'])
    assert state['status'] == 'COMPLETED' and state['result']['status'] == 'VALIDATED'
    assert client.get(f"/api/lab/jobs/{sub.json()['id']}/files/report.json").json() == {'status': 'VALIDATED'}
    assert client.get(f"/api/lab/jobs/{sub.json()['id']}/files/job.json").status_code == 422
    assert client.get(f"/api/lab/jobs/{sub.json()['id']}/files/missing.fits").status_code == 404
    assert client.get(f"/api/lab/jobs/{sub.json()['id']}/files/..%5Cjob.json").status_code in (404, 422)


def test_failures_are_reported_with_specific_codes(client, monkeypatch):
    def failing(kind, body, output, progress):
        raise requests.ConnectionError('502 Bad Gateway from IRSA')
    monkeypatch.setattr(api, '_execute', failing)
    sub = client.post('/api/lab/screen', json={'target_id': 'rw-aurigae'})
    state = _wait(client, sub.json()['id'])
    assert state['status'] == 'FAILED' and state['error']['code'] == 'ARCHIVE_UNAVAILABLE'
    assert client.get(f"/api/lab/jobs/{sub.json()['id']}/files/report.json").status_code == 409


def test_queue_is_bounded(client, monkeypatch):
    gate = threading.Event()

    def blocking(kind, body, output, progress):
        gate.wait(5)
        return {}
    monkeypatch.setattr(api, '_execute', blocking)
    try:
        codes = [client.post('/api/lab/screen', json={'target_id': 'rw-aurigae'}).status_code for _ in range(api.MAX_QUEUED + 1)]
        assert codes[:api.MAX_QUEUED] == [202]*api.MAX_QUEUED and codes[-1] == 429
    finally:
        gate.set()
        for job_id in list(api._jobs):
            _wait(client, job_id)


def test_restarted_server_reports_interrupted_jobs(client, tmp_path):
    job_id = str(uuid.uuid4()); directory = tmp_path/'jobs'/job_id; directory.mkdir(parents=True)
    (directory/'job.json').write_text(json.dumps({'id': job_id, 'status': 'RUNNING'}), encoding='utf-8')
    state = client.get(f'/api/lab/jobs/{job_id}').json()
    assert state['status'] == 'INTERRUPTED' and state['error']['code'] == 'WORKER_RESTARTED'
    assert client.get('/api/lab/jobs/not-a-uuid').status_code == 422
    assert client.get(f'/api/lab/jobs/{uuid.uuid4()}').status_code == 404


def test_resolve_presets_locally_and_reports_resolver_outage(client, monkeypatch):
    assert client.get('/api/lab/resolve', params={'name': "Barnard's Star"}).json()['target_id'] == 'barnards-star'

    def down(*args, **kwargs):
        raise requests.ConnectionError('timeout')
    monkeypatch.setattr(api.requests, 'get', down)
    response = client.get('/api/lab/resolve', params={'name': 'Vega'})
    assert response.status_code == 502 and response.json()['detail']['code'] == 'RESOLVER_UNAVAILABLE'


def test_resolve_parses_sesame(client, monkeypatch):
    class Reply:
        text = '# Vega\n%J 279.23473479 +38.78368896 = 18 36 56.33\n'
        def raise_for_status(self):
            return None
    monkeypatch.setattr(api.requests, 'get', lambda *a, **k: Reply())
    body = client.get('/api/lab/resolve', params={'name': 'Vega'}).json()
    assert body['ra_deg'] == pytest.approx(279.23473479) and body['dec_deg'] == pytest.approx(38.78368896)


def test_demo_listing_and_file_serving(client, tmp_path):
    demo = tmp_path/'demos'/'sample-run'; demo.mkdir(parents=True)
    (demo/'manifest.json').write_text(json.dumps({'name': 'sample-run', 'status': 'COMPLETED'}), encoding='utf-8')
    (demo/'report.json').write_text('{"ok": true}', encoding='utf-8')
    assert client.get('/api/lab/demos').json() == [{'name': 'sample-run', 'status': 'COMPLETED'}]
    assert client.get('/api/lab/demos/sample-run/files/report.json').json() == {'ok': True}
    assert client.get('/api/lab/demos/sample-run/files/secret.txt').status_code == 404
    assert client.get('/api/lab/demos/..%2F..%2Fx/files/report.json').status_code in (404, 422)
    assert client.get('/api/lab/demos/missing/files/report.json').status_code == 404


def test_index_and_static_assets_are_served(client):
    index = client.get('/')
    assert index.status_code == 200 and 'SPHEREx Epoch Lab' in index.text
    assert client.get('/static/app.js').status_code == 200
    assert client.get('/static/app.css').status_code == 200


def test_frontend_references_existing_elements():
    """Every $('id') used by app.js exists in index.html or is created by app.js itself."""
    html = Path('static/index.html').read_text(encoding='utf-8')
    js = Path('static/app.js').read_text(encoding='utf-8')
    used = set(re.findall(r"\$\('([A-Za-z0-9_-]+)'\)", js))
    created = set(re.findall(r"id: (?:interactive \? )?'([A-Za-z0-9_-]+)'", js))
    missing = sorted(i for i in used if f'id="{i}"' not in html and i not in created and not i.startswith('view-'))
    assert not missing, missing
    for view in ('explore', 'runs', 'method'):
        assert f'id="view-{view}"' in html


@pytest.mark.skipif(shutil.which('node') is None, reason='node is not installed')
def test_frontend_script_parses():
    subprocess.run(['node', '--check', 'static/app.js'], check=True)


# ---------------------------------------------------------------- regression tests for audited bugs

def test_null_observation_time_is_a_value_error_not_a_crash():
    from observation_selection import select_observations
    record = {'obs_publisher_did': PID_A, 't_min': None, 'energy_bandpassname': 'SPHEREx-D1'}
    with pytest.raises(ValueError, match='timestamp'):
        select_observations([record, dict(record, obs_publisher_did=PID_B)])


def test_discovery_rejects_records_without_time(monkeypatch):
    import pair_discovery
    region = 'POLYGON ICRS 9.9 9.9 10.1 9.9 10.1 10.1 9.9 10.1'
    base = {'s_region': region, 'access_url': 'u', 'access_format': 'image/fits', 'energy_bandpassname': 'SPHEREx-D1', 'em_min': 1e-6, 'em_max': 2e-6}
    records = [dict(base, obs_id='good', obs_publisher_did=PID_A, t_min=60900.0), dict(base, obs_id='bad', obs_publisher_did=PID_B, t_min=None)]
    monkeypatch.setattr(pair_discovery, 'search', lambda *a: {'query': {'ra_deg': 10.0, 'dec_deg': 10.0}, 'observations': records})
    result = pair_discovery.discover_pairs(10.0, 10.0)
    assert [o['obs_id'] for o in result['observations']] == ['good']
    assert result['rejected'] == [{'obs_id': 'bad', 'reason': 'missing or invalid observation time'}]


def test_search_cache_is_thread_safe(monkeypatch):
    import spherex
    xml = Path('artifacts/search_fixture.xml').read_bytes()

    class Reply:
        content = xml
        def raise_for_status(self):
            return None
    monkeypatch.setattr(spherex, '_get', lambda *a, **k: Reply())
    monkeypatch.setattr(spherex, 'DISK_CACHE_DIR', None)
    spherex._cache.clear()
    errors = []

    def worker(k):
        try:
            for j in range(20):
                spherex.search(10.0 + (k*20 + j) % 40 * 0.001, 5.0, 0.01)
        except Exception as exc:  # any KeyError/RuntimeError from unsynchronised access
            errors.append(exc)
    threads = [threading.Thread(target=worker, args=(k,)) for k in range(8)]
    [t.start() for t in threads]; [t.join() for t in threads]
    assert not errors and len(spherex._cache) <= 32


def test_busy_comparison_output_maps_to_409(client, monkeypatch):
    import scientific_pipeline

    def busy(**kwargs):
        raise ValueError('output directory is busy')
    monkeypatch.setattr(scientific_pipeline, 'run', busy)
    response = client.post('/api/comparisons', json={'ra': 10.0, 'dec': 5.0})
    assert response.status_code == 409


def test_sia_disk_cache_is_used_and_never_serves_after_errors(tmp_path, monkeypatch):
    import spherex
    xml = Path('artifacts/search_fixture.xml').read_bytes(); calls = []

    class Reply:
        content = xml
        def raise_for_status(self):
            return None
    monkeypatch.setattr(spherex, '_get', lambda *a, **k: calls.append(1) or Reply())
    monkeypatch.setattr(spherex, 'DISK_CACHE_DIR', str(tmp_path))
    spherex._cache.clear()
    first = spherex.search(127.69444, -39.1776, 0.01)
    spherex._cache.clear()
    second = spherex.search(127.69444, -39.1776, 0.01)
    assert first['count'] == second['count'] and len(calls) == 1
    monkeypatch.setattr(spherex, 'DISK_TTL', -1)    # expired entry must be refetched
    spherex._cache.clear()

    def fail(*a, **k):
        raise requests.ConnectionError('502')
    monkeypatch.setattr(spherex, '_get', fail)
    with pytest.raises(requests.ConnectionError):
        spherex.search(127.69444, -39.1776, 0.01)
