import numpy as np
import pytest
from spherex_lab.models.comparison_engine import compare_aligned, response_compatibility, write_outputs

def compare(a,b,va=None,vb=None,**kw):
    va=np.ones_like(a) if va is None else va
    vb=np.ones_like(b) if vb is None else vb
    return compare_aligned(a,va,b,vb,np.ones_like(a,dtype=bool),np.ones_like(b,dtype=bool),np.ones_like(a,dtype=bool),synthetic=True,**kw)

def test_identical_and_known_variance():
    a=np.arange(36.).reshape(6,6)
    r=compare(a,a)
    np.testing.assert_array_equal(r['difference'],0)
    np.testing.assert_array_equal(r['variance'],2)
    assert r['report']['candidates']==[]

def test_response_difference_rejected():
    assert response_compatibility([1,2,3],[0,1,0],[1,2,3],[0,1,0])
    assert not response_compatibility([1,2,3],[0,1,0],[1,2,3],[0,0,1])

def test_injected_motion_and_variable():
    a=np.zeros((20,20));b=a.copy();a[5,5]=20;b[5,8]=20;b[12,12]=10
    r=compare(a,b)
    assert len(r['report']['candidates'])==3
    assert r['difference'][5,5]<0 and r['difference'][5,8]>0

def test_mask_and_missing_spectrum():
    a=np.ones((4,4)); valid=np.ones_like(a,dtype=bool);valid[1,1]=False
    r=compare_aligned(a,a,a,a,valid,valid,np.ones_like(valid))
    assert np.isnan(r['difference'][1,1])
    with pytest.raises(ValueError,match='overlap'):
        compare_aligned(a,a,a,a,valid,valid,np.zeros_like(valid))

def test_background_plane():
    y,x=np.indices((10,10));a=np.zeros((10,10));b=3+2*x+.1*y
    r=compare(a,b,blank_mask=np.ones_like(a,dtype=bool))
    np.testing.assert_allclose(r['difference'],0,atol=1e-12)
    assert np.nanmin(r['variance'])>=0

def test_outputs(tmp_path):
    r=compare(np.zeros((4,4)),np.ones((4,4)))
    write_outputs(r,tmp_path)
    assert (tmp_path/'comparison.fits').exists() and (tmp_path/'report.json').exists()

def test_defined_noise_distribution_false_positive_count():
    rng=np.random.default_rng(184)
    counts=[]
    for _ in range(10):
        r=compare(rng.normal(size=(64,64)),rng.normal(size=(64,64)),sigma_threshold=5)
        counts.append(len(r['report']['candidates']))
    # A finite deterministic injection-recovery distribution, not a sky false-positive claim.
    assert sum(counts) <= 3

def test_background_fit_covariance_at_blank_pixel():
    shape=(6,6);a=np.zeros(shape);valid=np.ones(shape,dtype=bool)
    r=compare(a,a,blank_mask=valid)
    y,x=np.indices(shape);z=np.stack([np.ones(shape),x,y],axis=-1).reshape(-1,3)
    hat=np.diag(z@np.linalg.inv(z.T@z)@z.T).reshape(shape)
    np.testing.assert_allclose(r['variance'],2*(1-hat),atol=1e-12)

def test_resampled_background_fit_rejected():
    a=np.zeros((4,4))
    with pytest.raises(ValueError,match='covariance'):
        compare(a,a,blank_mask=np.ones_like(a,dtype=bool),resampled=True)

def test_nonfinite_response_grid_rejected():
    with pytest.raises(ValueError): response_compatibility([1,np.inf],[1,1],[1,np.inf],[1,1])
