import pytest
from pathlib import Path

def test_download_refuses_non_archive(tmp_path):
    from scientific_pipeline import download_cutout
    with pytest.raises(ValueError): download_cutout({'access_url':'https://example.com/x'}, 0,0,tmp_path)

def test_run_validation(tmp_path):
    from scientific_pipeline import run
    with pytest.raises(ValueError): run(ra=361,dec=0,output=tmp_path)

def test_comparison_api_input_and_status():
    from fastapi.testclient import TestClient
    from main import app
    c=TestClient(app)
    assert c.post('/api/comparisons',json={'ra':361,'dec':0}).status_code==422
    assert c.get('/api/comparisons/not-a-uuid').status_code==422
    assert c.get('/api/comparisons/00000000-0000-0000-0000-000000000000').status_code==404

def test_checkpoint_resume_and_changed_sky_invalidate(tmp_path,monkeypatch):
    import scientific_pipeline as pipeline
    import pair_discovery, response_calibration, calibration_resolver
    import numpy as np
    from fits_pipeline import read_cutout
    sample='artifacts/sample.fits';image=read_cutout(sample)['image'];calls=[]
    record={'obs_id':'2025W18_1B_0162_1','obs_publisher_did':'ivo://irsa.ipac/spherex_dr1?2025W18_1B_0162_1/D1'}
    monkeypatch.setattr(pair_discovery,'discover_pairs',lambda *a,**k:{'pairs':[{'a':record,'b':record}]})
    monkeypatch.setattr(response_calibration,'fetch_bundle',lambda *a,**k:None)
    def missing(*a,**k):raise ValueError('fixture has no response')
    monkeypatch.setattr(response_calibration,'pixel_response',missing)
    monkeypatch.setattr(pipeline,'wavelength_maps',lambda *a,**k:(np.ones_like(image),np.ones_like(image),np.ones_like(image,dtype=bool),{'unit':'um'}))
    def download(*a,**k):calls.append(1);return sample
    monkeypatch.setattr(pipeline,'download_cutout',download)
    kwargs={'ra':127.69444,'dec':-39.1776,'output':tmp_path,'max_pairs':1}
    first=pipeline.run(**kwargs)
    assert len(calls)==2 and first['attempts'][0]['status']=='REJECTED'
    pipeline.run(**kwargs)
    assert len(calls)==2
    pipeline.run(**{**kwargs,'ra':127.69445})
    assert len(calls)==4
