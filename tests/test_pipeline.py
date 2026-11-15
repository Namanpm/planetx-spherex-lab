import io
import pytest
import requests
from fastapi.testclient import TestClient
import main
from spherex_lab.data import spherex
client = TestClient(main.app)

def test_invalid_coordinates():
    assert client.get('/api/spherex/search', params={'ra': 360, 'dec': 0}).status_code == 422
    assert client.get('/api/spherex/search', params={'ra': 1, 'dec': -91}).status_code == 422
    assert client.get('/api/spherex/search', params={'ra': 1, 'dec': 0, 'radius': 0}).status_code == 422

def test_timeout_is_not_empty_result(monkeypatch):
    def fail(*args, **kwargs): raise requests.Timeout('timeout')
    monkeypatch.setattr(main, 'search', fail)
    assert client.get('/api/spherex/search', params={'ra': 1, 'dec': 0}).status_code == 504

def test_real_votable_parser_and_cache(monkeypatch):
    data = open('tests/data/search_fixture.xml', 'rb').read()
    class Response:
        content = data
        def raise_for_status(self): pass
    calls = []
    def get(*args, **kwargs):
        calls.append(kwargs)
        return Response()
    monkeypatch.setattr(spherex, '_get', get)
    monkeypatch.setattr(spherex, 'DISK_CACHE_DIR', None)
    spherex._cache.clear()
    a = spherex.search(127.69444, -39.17760, .01)
    b = spherex.search(127.69444, -39.17760, .01)
    assert a['count'] == 414 and a == b and len(calls) == 1
    assert a['observations'][0]['access_url'].startswith('https://irsa.ipac.caltech.edu/')
    assert a['observations'][0]['s_ra'] is not None

def test_archive_detects_overflow_status():
    import io
    from astropy.io.votable import parse
    from spherex_lab.data.spherex import check_query_status
    xml = b'''<?xml version="1.0"?><VOTABLE version="1.3" xmlns="http://www.ivoa.net/xml/VOTable/v1.3"><RESOURCE type="results"><INFO name="QUERY_STATUS" value="OK"/><TABLE><FIELD name="a" datatype="int"/><DATA><TABLEDATA/></DATA></TABLE><INFO name="QUERY_STATUS" value="OVERFLOW">truncated</INFO></RESOURCE></VOTABLE>'''
    with pytest.raises(ValueError, match='OVERFLOW'):
        check_query_status(parse(io.BytesIO(xml)))
