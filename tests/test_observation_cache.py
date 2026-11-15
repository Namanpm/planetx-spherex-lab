from pathlib import Path
import numpy as np
import pytest
from astropy.io import fits

SAMPLE=Path('tests/data/sample.fits')

def test_prepare_once_reuse_arrays_and_psf_without_reopening(monkeypatch):
    from spherex_lab.imaging.native_photometry import get_observation,fit_source
    calls=[];original=fits.open
    def counted(*args,**kwargs):calls.append(args[0]);return original(*args,**kwargs)
    monkeypatch.setattr(fits,'open',counted)
    cache={};key='2025W18_1B_0162_1/D1'
    observation=get_observation(SAMPLE,key,cache)
    assert get_observation(SAMPLE,key,cache) is observation
    sky=observation['science']['wcs'].celestial.pixel_to_world(29.1,14.3)
    first=fit_source(observation,sky.ra.deg,sky.dec.deg)
    second=fit_source(observation,sky.ra.deg,sky.dec.deg)
    assert first==second and len(calls)==1 and list(cache)==[key]
    assert len(observation['epsf']['zones'])==1
    assert not observation['science']['image'].flags.writeable

def test_cache_mismatch_cannot_reuse_other_detector():
    from spherex_lab.imaging.native_photometry import get_observation
    cache={}
    with pytest.raises(ValueError,match='identity'):get_observation(SAMPLE,'2025W18_1B_0162_1/D2',cache)
    assert not cache

def test_new_cutout_invalidates_same_observation_key(tmp_path):
    from spherex_lab.imaging.native_photometry import get_observation
    cache={};key='2025W18_1B_0162_1/D1'
    first=get_observation(SAMPLE,key,cache)
    path=tmp_path/'other.fits';path.write_bytes(SAMPLE.read_bytes())
    second=get_observation(path,key,cache)
    assert first is not second and len(cache)==1
    np.testing.assert_array_equal(first['science']['image'],second['science']['image'])

def test_cached_path_and_prepared_photometry_identical():
    from spherex_lab.imaging.native_photometry import get_observation,fit_source
    observation=get_observation(SAMPLE,'2025W18_1B_0162_1/D1',{})
    sky=observation['science']['wcs'].celestial.pixel_to_world(29.1,14.3)
    assert fit_source(SAMPLE,sky.ra.deg,sky.dec.deg)==fit_source(observation,sky.ra.deg,sky.dec.deg)

def test_cached_provenance_does_not_mutate_other_source_results():
    from spherex_lab.imaging.native_photometry import get_observation
    from spherex_lab.imaging.psf_processing import extract_prepared_epsf
    o=get_observation(SAMPLE,'2025W18_1B_0162_1/D1',{})
    first=extract_prepared_epsf(o['epsf'],29.1,14.3)
    first['provenance']['limitations'].append('caller annotation')
    first['provenance']['oversampling_xy'][0]+=1
    first['provenance']['zone_center_detector_xy'][0]=-999
    second=extract_prepared_epsf(o['epsf'],29.1,14.3)
    assert 'caller annotation' not in second['provenance']['limitations']
    assert second['provenance']['oversampling_xy']==[5,5]
    assert second['provenance']['zone_center_detector_xy'][0]!=-999

def test_changed_file_refreshes_values_without_changing_old_entry(tmp_path):
    from spherex_lab.imaging.native_photometry import get_observation
    import os
    path=tmp_path/'sample.fits';path.write_bytes(SAMPLE.read_bytes());cache={};key='2025W18_1B_0162_1/D1'
    original=get_observation(path,key,cache);previous=float(original['science']['image'][14,29]);stamp=path.stat().st_mtime_ns
    with fits.open(path,mode='update',memmap=False) as hdus:hdus['IMAGE'].data[14,29]+=7
    os.utime(path,ns=(stamp+1000000000,stamp+1000000000))
    refreshed=get_observation(path,key,cache)
    assert refreshed is not original
    assert float(original['science']['image'][14,29])==previous
    assert float(refreshed['science']['image'][14,29])==pytest.approx(previous+7,rel=1e-6)

def test_same_pointing_different_detectors_have_distinct_entries(tmp_path):
    from spherex_lab.imaging.native_photometry import get_observation
    path=tmp_path/'detector2.fits'
    with fits.open(SAMPLE,memmap=False) as hdus:
        hdus['IMAGE'].header['DETECTOR']=2
        hdus['EPSF'].header['DETECTOR']='2'
        hdus.writeto(path)
    cache={};base='2025W18_1B_0162_1'
    first=get_observation(SAMPLE,base+'/D1',cache)
    second=get_observation(path,base+'/D2',cache)
    assert len(cache)==2
    assert first['epsf']['detector']==1 and second['epsf']['detector']==2
