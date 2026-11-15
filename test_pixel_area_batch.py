import numpy as np
import pytest
from astropy.wcs import WCS, Sip
import native_photometry as native


def make_wcs(ra=359.999,dec=20.,distorted=False):
    w=WCS(naxis=2)
    w.wcs.crpix=[50,50];w.wcs.crval=[ra,dec]
    w.wcs.cdelt=[-6/3600,6/3600];w.wcs.ctype=['RA---TAN','DEC--TAN']
    if distorted:
        a=np.zeros((3,3));b=a.copy();a[2,0]=.0003;b[0,2]=-.0002
        w.sip=Sip(a,b,None,None,w.wcs.crpix)
        w.wcs.ctype=['RA---TAN-SIP','DEC--TAN-SIP']
    return w


@pytest.mark.parametrize('ra,dec,distorted',[(359.999,20.,False),(12.,89.99,False),(30.,-89.99,True),(42.,-20.,True)])
def test_batch_matches_original_spherical_formula_at_wrap_poles_and_sip(ra,dec,distorted):
    w=make_wcs(ra,dec,distorted)
    x=np.array([0.,.5,49.,98.5,99.]);y=np.array([99.,0.,49.,98.5,.5])
    actual=native.local_pixel_areas_sr(w,x,y)
    expected=np.array([native.local_pixel_area_sr(w,a,b) for a,b in zip(x,y)])
    np.testing.assert_allclose(actual,expected,rtol=1e-13,atol=0)
    assert np.all(actual>0)


def test_prepared_repeated_geometry_has_no_extra_transforms(monkeypatch):
    w=make_wcs(distorted=True);obs={'wcs':w,'science':{'image':np.zeros((100,100))}}
    sky=w.pixel_to_world([20.,49.,80.],[30.,49.,70.]);positions=list(zip(sky.ra.deg,sky.dec.deg))
    calls=[];original=native.local_pixel_areas_sr
    def counted(*args):
        calls.append(len(args[1]));return original(*args)
    monkeypatch.setattr(native,'local_pixel_areas_sr',counted)
    assert native.prepare_pixel_areas(obs,positions)==3
    assert native.prepare_pixel_areas(obs,positions)==0
    for ra,dec in positions:
        x,y=w.world_to_pixel(native.SkyCoord(ra,dec,unit='deg'))
        expected=native.local_pixel_area_sr(w,float(x),float(y))
        assert native._source_pixel_area(obs,x,y)==pytest.approx(expected,rel=1e-13)
    assert calls==[3]
    assert obs['_pixel_area_cache']['hits']==3
    assert obs['_pixel_area_cache']['misses']==0
    assert np.ptp(list(obs['_pixel_area_cache']['values'].values()))>0


def test_invalid_and_outside_positions_do_not_poison_valid_batch():
    w=make_wcs();obs={'wcs':w,'science':{'image':np.zeros((100,100))}}
    good=w.pixel_to_world(49.,49.);outside=w.pixel_to_world(-10.,49.)
    positions=[(np.nan,20.),(12.,91.),(-1.,20.),(360.,20.),(outside.ra.deg,outside.dec.deg),(good.ra.deg,good.dec.deg)]
    assert native.prepare_pixel_areas(obs,positions)==1
    assert len(obs['_pixel_area_cache']['values'])==1


def test_geometry_mutation_invalidates_on_new_batch():
    w=make_wcs();obs={'wcs':w,'science':{'image':np.zeros((100,100))}}
    sky=w.pixel_to_world(49.,49.);positions=[(sky.ra.deg,sky.dec.deg)]
    native.prepare_pixel_areas(obs,positions)
    before=list(obs['_pixel_area_cache']['values'].values())[0]
    w.wcs.cdelt*=2
    assert native.prepare_pixel_areas(obs,positions)==1
    after=list(obs['_pixel_area_cache']['values'].values())[0]
    assert after==pytest.approx(before*4,rel=1e-8)


def test_area_cache_is_bounded_and_scalar_fallback_reuses(monkeypatch):
    monkeypatch.setattr(native,'PIXEL_AREA_CACHE_SIZE',2)
    w=make_wcs();obs={'wcs':w,'science':{'image':np.zeros((100,100))}}
    sky=w.pixel_to_world([10.,30.,50.],[20.,40.,60.])
    native.prepare_pixel_areas(obs,zip(sky.ra.deg,sky.dec.deg))
    assert len(obs['_pixel_area_cache']['values'])==2
    value=native._source_pixel_area(obs,1.,2.)
    assert native._source_pixel_area(obs,1.,2.)==value
    assert len(obs['_pixel_area_cache']['values'])==2
    assert obs['_pixel_area_cache']['misses']==1
    assert obs['_pixel_area_cache']['hits']==1


def test_native_fit_batch_and_scalar_paths_match():
    obs=native.prepare_observation('artifacts/sample.fits')
    sky=obs['wcs'].pixel_to_world(29.1,14.3)
    original=native.fit_source(obs,sky.ra.deg,sky.dec.deg)
    obs.pop('_pixel_area_cache')
    native.prepare_pixel_areas(obs,[(sky.ra.deg,sky.dec.deg)])
    batched=native.fit_source(obs,sky.ra.deg,sky.dec.deg)
    assert batched['central_flag']==original['central_flag']
    np.testing.assert_allclose([batched['flux_jy'],batched['variance_jy2']],[original['flux_jy'],original['variance_jy2']],rtol=1e-13,atol=0)


def test_concurrent_geometry_cache_access():
    from concurrent.futures import ThreadPoolExecutor
    w=make_wcs();obs={'wcs':w,'science':{'image':np.zeros((100,100))}}
    with ThreadPoolExecutor(max_workers=4) as pool:
        values=list(pool.map(lambda _:native._source_pixel_area(obs,49.,49.),range(20)))
    assert values==[values[0]]*20
    assert obs['_pixel_area_cache']['misses']==1
    assert obs['_pixel_area_cache']['hits']==19


def test_nonfinite_pixel_geometry_is_isolated():
    areas=native.local_pixel_areas_sr(make_wcs(),np.array([49.,np.nan]),np.array([49.,49.]))
    assert np.isfinite(areas[0]) and np.isnan(areas[1])
