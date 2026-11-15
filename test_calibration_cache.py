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
def test_datalink_transport_fallback_preserves_id(body, valid):
    from urllib.parse import parse_qs, urlparse
    import requests
    class FallbackClient:
        def __init__(self): self.calls = []
        def get(self, url, **kwargs):
            self.calls.append((url, kwargs))
            if len(self.calls) == 1:
                response = requests.Response()
                response.status_code = 502
                raise requests.HTTPError(response=response)
            assert parse_qs(urlparse(url).query)['ID'] == [PID]
            return Client(body).get(url)
    client = FallbackClient()
    if valid:
        assert resolver.discover(PID, 1, client) == URL
    else:
        with pytest.raises(Exception): resolver.discover(PID, 1, client)
    assert len(client.calls) == 2
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


@pytest.mark.parametrize('failures', [2, 3])
def test_three_equivalent_transport_variants_bounded(failures):
    from urllib.parse import parse_qs, urlparse
    import requests
    class GatewayClient:
        calls = 0
        def get(self, url, **kwargs):
            self.calls += 1
            actual_id = kwargs.get('params', {}).get('ID')
            if actual_id is None:
                actual_id = parse_qs(urlparse(url).query)['ID'][0]
            assert actual_id == PID
            if self.calls <= failures:
                response = requests.Response()
                response.status_code = 503
                raise requests.HTTPError(response=response)
            return Client(xml(row(URL))).get(url)
    client = GatewayClient()
    if failures == 2:
        assert resolver.discover(PID, 1, client) == URL
    else:
        with pytest.raises(requests.HTTPError): resolver.discover(PID, 1, client)
    assert client.calls == 3
