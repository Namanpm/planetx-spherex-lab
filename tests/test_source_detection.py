import numpy as np
import pytest
from spherex_lab.imaging.source_detection import detect_sources, crossmatch_catalog


def test_isolated_source_centroid_and_conditional_uncertainty():
    image=np.zeros((10,10));variance=np.ones_like(image);valid=np.ones_like(image,dtype=bool)
    image[4,3]=10;image[4,4]=10
    result=detect_sources(image,variance,valid,sigma_threshold=5.)
    assert len(result['sources'])==1
    source=result['sources'][0]
    assert source['x']==3.5 and source['y']==4. and source['flux_sum']==20
    np.testing.assert_allclose(source['centroid_covariance'],[[.00125,0],[0,0]])
    assert source['flux_variance']==2 and source['pixel_count']==2


def test_detection_masks_bad_pixels_and_does_not_merge_separate_sources():
    image=np.zeros((8,8));image[1,1]=20;image[5,5]=10;image[5,6]=30
    valid=np.ones_like(image,dtype=bool);valid[5,6]=False
    result=detect_sources(image,np.ones_like(image),valid,sigma_threshold=5)
    assert len(result['sources'])==2
    assert sorted(s['flux_sum'] for s in result['sources'])==[10,20]
    assert detect_sources(image,np.ones_like(image),valid,sigma_threshold=5,min_pixels=2)['sources']==[]


def test_correlated_detection_does_not_report_independent_uncertainty():
    a=np.zeros((4,4));a[1,1]=10
    source=detect_sources(a,np.ones_like(a),np.ones_like(a,dtype=bool),sigma_threshold=5,correlated=True)['sources'][0]
    assert source['centroid_covariance'] is None and source['flux_variance'] is None


def test_negative_invalid_and_nonpositive_variance_not_sources():
    a=np.zeros((4,4));a[1,1]=-10;a[1,2]=np.nan;a[2,2]=10
    v=np.ones_like(a);v[2,2]=0
    assert detect_sources(a,v,np.ones_like(a,dtype=bool),sigma_threshold=5)['sources']==[]
    with pytest.raises(ValueError): detect_sources(a,v,np.ones_like(a,dtype=bool),sigma_threshold=0)


def test_catalog_spherical_wrap_and_ambiguous_match():
    result=crossmatch_catalog([[359.9999,0],[10,20]],[[.0001,0],[10.0001,20],[9.9999,20]],max_separation_arcsec=1,catalog_provenance='Synthetic catalog')
    assert len(result['matches'])==1 and result['matches'][0]['source_index']==0
    np.testing.assert_allclose(result['matches'][0]['separation_arcsec'],.72,rtol=1e-9)
    assert result['ambiguous_source_indices']==[1]
    with pytest.raises(ValueError): crossmatch_catalog([[0,0]],[[0,0]],max_separation_arcsec=1,catalog_provenance='')

def test_defined_source_injection_recovery_and_false_positive_distribution():
    rng=np.random.default_rng(84712)
    recovered=0;false_positives=0
    yy,xx=np.indices((32,32))
    for _ in range(12):
        noise=rng.normal(size=(32,32))
        blank=detect_sources(noise,np.ones_like(noise),np.ones_like(noise,dtype=bool),sigma_threshold=5.)
        false_positives+=len(blank['sources'])
        injected=noise+30*np.exp(-((xx-16.2)**2+(yy-15.7)**2)/(2*1.1**2))
        result=detect_sources(injected,np.ones_like(injected),np.ones_like(injected,dtype=bool),sigma_threshold=5.)
        recovered+=int(any(np.hypot(s['x']-16.2,s['y']-15.7)<.3 for s in result['sources']))
    assert recovered==12 and false_positives==0
