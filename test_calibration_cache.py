import hashlib
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pytest
import calibration_resolver as resolver
from test_calibration_resolver import PID, URL, Client, xml, row


def product(tmp_path, monkeypatch):
    path = tmp_path / (hashlib.sha256(URL.encode()).hexdigest() + '.fits')
    path.write_bytes(b'calibration-v1')
    calls = []
    def validate(path, detector):
        calls.append(detector)
        if detector != 1:
            raise ValueError('calibration detector mismatch')
        return np.ones((2, 2)), np.ones((2, 2)) * 2
    monkeypatch.setattr(resolver, 'validate_calibration', validate)
    return path, calls


def test_persistent_association_identity_and_expiry(tmp_path, monkeypatch):
    client = Client(xml(row(URL)))
    calls = []
    original = resolver.discover
    monkeypatch.setattr(resolver, 'discover', lambda *a: (calls.append(a), original(*a))[1])
    cache = resolver.CalibrationCache()
    assert cache.association(PID, 1, {'version': 1}, tmp_path, client) == URL
    # New process-equivalent instance reads disk, without network metadata.
    assert resolver.CalibrationCache().association(PID, 1, {'version': 1}, tmp_path, client) == URL
    cache.association(PID, 1, {'version': 2}, tmp_path, client)
    resolver.CalibrationCache(association_ttl=0).association(PID, 1, {'version': 1}, tmp_path, client)
    assert len(calls) == 3
    with pytest.raises(ValueError, match='detector'):
        cache.association(PID, 2, {'version': 1}, tmp_path, client)


def test_array_hits_readonly_detector_and_file_invalidation(tmp_path, monkeypatch):
    path, calls = product(tmp_path, monkeypatch)
    cache = resolver.CalibrationCache()
    first = cache.arrays(URL, 1, tmp_path)
    assert cache.arrays(URL, 1, tmp_path)[0] is first[0]
    with pytest.raises(ValueError, match='read-only'):
        first[0][0, 0] = 5
    with pytest.raises(ValueError, match='detector'):
        cache.arrays(URL, 2, tmp_path)
    path.write_bytes(b'changed-calibration-v2')
    second = cache.arrays(URL, 1, tmp_path)
    assert second[2] != first[2]
    assert calls == [1, 1]
    assert cache.cache_info()['digests'] == 2
    assert cache.cache_info()['array_entries'] == 1
    monkeypatch.setattr(resolver, 'VALIDATOR_VERSION', 999)
    cache.arrays(URL, 1, tmp_path)
    assert calls == [1, 1, 1]


def test_thread_singleflight_and_memory_budget(tmp_path, monkeypatch):
    _, calls = product(tmp_path, monkeypatch)
    cache = resolver.CalibrationCache(max_bytes=64)
    with ThreadPoolExecutor(max_workers=4) as pool:
        values = list(pool.map(lambda _: cache.arrays(URL, 1, tmp_path), range(16)))
    assert len(calls) == 1
    assert all(value[0] is values[0][0] for value in values)
    assert cache.cache_info()['array_bytes'] == 64
    tiny = resolver.CalibrationCache(max_bytes=1)
    tiny.arrays(URL, 1, tmp_path)
    assert tiny.cache_info()['array_bytes'] == 0
    assert tiny.cache_info()['array_entries'] == 0


def test_corrupt_association_rediscovered(tmp_path):
    cache = resolver.CalibrationCache()
    cache.association(PID, 1, {}, tmp_path, Client(xml(row(URL))))
    path = next((tmp_path / 'associations').glob('*.json'))
    path.write_text('{malformed')
    cache.association(PID, 1, {}, tmp_path, Client(xml(row(URL))))
    assert cache.cache_info()['association_misses'] == 2


def test_array_lru_evicts_to_budget(tmp_path, monkeypatch):
    product(tmp_path, monkeypatch)
    second_url = URL.replace('a.fits', 'b.fits')
    second = tmp_path / (hashlib.sha256(second_url.encode()).hexdigest() + '.fits')
    second.write_bytes(b'second')
    cache = resolver.CalibrationCache(max_bytes=64)
    cache.arrays(URL, 1, tmp_path)
    cache.arrays(second_url, 1, tmp_path)
    assert cache.cache_info()['evictions'] == 1
    assert cache.cache_info()['array_bytes'] == 64
    cache.arrays(URL, 1, tmp_path)
    assert cache.cache_info()['array_misses'] == 3


@pytest.mark.parametrize('body, valid', [(xml(row(URL)), True), (b'malformed', False)])
def test_datalink_transport_error_is_not_replayed_by_resolver(body, valid):
    import requests
    class GatewayClient:
        def __init__(self): self.calls = []
        def get(self, url, **kwargs):
            self.calls.append((url, kwargs))
            response = requests.Response()
            response.status_code = 502
            raise requests.HTTPError(response=response)
    client = GatewayClient()
    with pytest.raises(requests.HTTPError): resolver.discover(PID, 1, client)
    assert len(client.calls) == 1
    assert client.calls[0][0] == resolver.DATALINK
    assert client.calls[0][1]['params'] == {'ID': PID}


def test_datalink_normal_request_not_repeated():
    class CountingClient(Client):
        calls = 0
        def get(self, *a, **kw):
            self.calls += 1
            return super().get(*a, **kw)
    client = CountingClient(xml(row(URL)))
    assert resolver.discover(PID, 1, client) == URL
    assert client.calls == 1


def test_datalink_adapter_stops_after_eleven_attempts(monkeypatch):
    from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
    from threading import Thread
    import requests
    monkeypatch.setattr('archive_http.BoundedRetry.sleep', lambda self, response=None: None)
    calls = []
    class Handler(BaseHTTPRequestHandler):
        def do_GET(self):
            calls.append(self.path)
            self.send_response(502); self.end_headers()
        def log_message(self, *args): pass
    server = ThreadingHTTPServer(('127.0.0.1', 0), Handler)
    thread = Thread(target=server.serve_forever, daemon=True); thread.start()
    monkeypatch.setattr(resolver, 'DATALINK', f'http://127.0.0.1:{server.server_port}/datalink')
    try:
        with resolver.session() as client:
            with pytest.raises(requests.exceptions.RetryError): resolver.discover(PID, 1, client)
        assert len(calls) == 11
        assert all(path.startswith('/datalink?ID=') for path in calls)
    finally:
        server.shutdown(); thread.join(); server.server_close()
