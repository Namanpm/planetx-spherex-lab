import numpy as np
import pytest

def test_signed_distance_identity_scale_and_disjoint():
    from spherex_lab.imaging.reference_response import signed_laboratory_distance
    a=(np.array([1.,2.,3.]),np.array([-.1,1.,-.1]))
    same=signed_laboratory_distance(a,(a[0],a[1]*2))
    assert same['normalized_l2_distance']==pytest.approx(0,abs=1e-14)
    assert same['negative_sample_counts']==[2,2]
    shifted=signed_laboratory_distance(a,(a[0]+4,a[1]))
    assert shifted['normalized_l2_distance']==pytest.approx(np.sqrt(2))
    assert shifted['acceptance']=='NOT_ASSESSED'

def test_field_response_common_sky_and_mask(tmp_path,monkeypatch):
    from astropy.io import fits
    from astropy.wcs import WCS
    from spherex_lab.imaging import reference_response as module
    w=WCS(naxis=2);w.wcs.crpix=[3,3];w.wcs.crval=[10,20];w.wcs.cdelt=[-.001,.001];w.wcs.ctype=['RA---TAN','DEC--TAN']
    h=w.to_header();h['BUNIT']='MJy/sr';h['DETECTOR']=1
    v=h.copy();v['BUNIT']='MJy2/sr2'
    p=tmp_path/'science.fits';fits.HDUList([fits.PrimaryHDU(),fits.ImageHDU(np.ones((5,5)),h,name='IMAGE'),fits.ImageHDU(np.ones((5,5)),v,name='VARIANCE'),fits.ImageHDU(np.zeros((5,5),dtype='int32'),name='FLAGS')]).writeto(p)
    calls=[]
    def response(bundle,header,shape,x,y,center):
        calls.append((x,y,center));return np.array([1.,2.]),np.array([0.,1.]),{'status':'APPROXIMATE'}
    monkeypatch.setattr(module,'pixel_response',response)
    curve=module.field_response(p,np.full((5,5),1.5),None,10,20)
    assert calls==[(2,2,1.5)] and curve[2]['sampled_pixel_xy']==[2,2]
    with pytest.raises(ValueError):module.field_response(p,np.ones((4,4)),None,10,20)

def test_signed_distance_piecewise_linear_oracle_and_bad_samples():
    from spherex_lab.imaging.reference_response import signed_laboratory_distance
    a=(np.array([-1.,0.,1.]),np.array([0.,1.,0.]))
    shifted=(a[0]+.5,a[1])
    # Triangle autocorrelation at lag .5 is 2/3-.5**2+.5**3/2.
    assert signed_laboratory_distance(a,shifted)['normalized_l2_distance']==pytest.approx(.75)
    assert signed_laboratory_distance(([1,2],[1,1]),([2,3],[1,1]))['normalized_l2_distance']==pytest.approx(np.sqrt(2))
    with pytest.raises(ValueError):signed_laboratory_distance(([1,np.nan],[1,1]),a)
