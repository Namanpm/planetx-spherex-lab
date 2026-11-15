import pytest

def test_screening_variable_duplicate_invalid():
    from spherex_lab.data.reference_catalog import screen_rows
    base=dict(source_id='1',ra='10',dec='20',pmra='1',pmdec='2',ref_epoch='2016',bp_rp='1.1',phot_variable_flag='NOT_AVAILABLE',duplicated_source='false',ruwe='1.1',phot_g_mean_mag='14')
    sources,rejected=screen_rows([base,{**base,'source_id':'2','phot_variable_flag':'VARIABLE'},{**base,'source_id':'3','duplicated_source':'true'},{**base,'source_id':'4','bp_rp':''}])
    assert [s['source_id'] for s in sources]==['1']
    assert len(rejected)==3 and sources[0]['stability_status']=='not classified variable by Gaia; constancy unproven'

def test_position_epoch_and_validation():
    from spherex_lab.data.reference_catalog import source_position,query_references
    row=dict(ra=10,dec=20,pmra=1000,pmdec=0,ref_epoch=2016)
    position=source_position(row,'2026-01-01')
    assert position['ra_deg']>10 and abs(position['dec_deg']-20)<1e-4
    with pytest.raises(ValueError):query_references(360,0)

def test_cached_catalog_preserves_rejected_neighbors(tmp_path,monkeypatch):
    from spherex_lab.data import reference_catalog as catalog
    header=','.join(catalog.COLUMNS)+'\n'
    values=['7','10','20','1','2','2016','VARIABLE','false','1.1','1','14','1','1','1','1','1','1']
    body=(header+','.join(values)+'\n').encode()
    class Response:
        def __enter__(self):return self
        def __exit__(self,*args):pass
        def raise_for_status(self):pass
        def iter_content(self,*args):yield body
    class Client:
        def __enter__(self):return self
        def __exit__(self,*args):pass
        def get(self,*args,**kwargs):return Response()
    monkeypatch.setattr(catalog,'session',lambda:Client())
    r=catalog.query_references(10,20,cache_dir=tmp_path)
    assert r['sources']==[] and len(r['all_sources'])==1
    assert r['all_sources'][0]['phot_variable_flag']=='VARIABLE'
    assert r['provenance']['truncated'] is False
    assert 'magnitude_selection' in r['provenance']

def test_gaia_network_failure_has_specific_bounded_retry_error(tmp_path,monkeypatch):
    import requests
    from spherex_lab.data import reference_catalog as catalog
    class Client:
        def __enter__(self): return self
        def __exit__(self,*args): pass
        def get(self,*args,**kwargs):
            assert kwargs['timeout']==(10,60) and kwargs['stream'] is True
            raise requests.exceptions.ConnectTimeout('test timeout')
    monkeypatch.setattr(catalog,'session',lambda:Client())
    with pytest.raises(catalog.GaiaArchiveUnavailable,match='bounded retries'):
        catalog.query_references(10,20,cache_dir=tmp_path)
