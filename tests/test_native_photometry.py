import numpy as np
import pytest
from spherex_lab.imaging.native_photometry import fit_stamp,sample_epsf,local_pixel_area_sr
from astropy.wcs import WCS

def test_joint_native_fit_recovers_source_and_plane():
    y,x=np.mgrid[-3:4,-3:4];p=np.exp(-(x*x+y*y)/2);p/=p.sum();a=17*p+4+.1*x-.2*y
    r=fit_stamp(a,np.ones_like(a)*.04,np.ones_like(a,bool),p,x,y)
    assert r['amplitude']==pytest.approx(17)
    np.testing.assert_allclose(r['background_coefficients'],[4,.1,-.2],atol=1e-12)
    assert r['variance']>0 and r['reduced_chi2']<1e-20

def test_subpixel_sampling_shift_and_bounds():
    y,x=np.mgrid[-10:11,-10:11];raw=np.exp(-(x*x+y*y)/8)
    p,good=sample_epsf(raw,(2,2),np.array([0,.5,8]),np.zeros(3))
    assert p[0]==1 and p[1]<p[0] and not good[2] and np.isnan(p[2])

def test_singular_or_masked_fit_rejected():
    a=np.ones((4,4));y,x=np.indices(a.shape)
    with pytest.raises(ValueError):fit_stamp(a,a,a.astype(bool),a,x,y)
    with pytest.raises(ValueError):fit_stamp(a,a,np.zeros_like(a,bool),a,x,y)

def test_wcs_solid_angle():
    w=WCS(naxis=2);w.wcs.crpix=[10,10];w.wcs.crval=[10,20];w.wcs.cdelt=[-1/3600,1/3600];w.wcs.ctype=['RA---TAN','DEC--TAN']
    assert local_pixel_area_sr(w,9,9)==pytest.approx(np.deg2rad(1/3600)**2,rel=1e-6)

def test_covariance_matches_independent_linear_oracle():
    y,x=np.mgrid[-3:4,-3:4];p=np.exp(-(x*x+y*y)/2);p/=p.sum();v=.1+.01*(x+3)
    valid=np.ones_like(p,bool);valid[2,2]=False
    r=fit_stamp(11*p+4,v,valid,p,x,y)
    design=np.stack([p,np.ones_like(p),x,y],axis=-1)[valid]
    expected=np.linalg.inv(design.T@(design/v[valid,None]))
    np.testing.assert_allclose(r['parameter_covariance'],expected,rtol=1e-11,atol=1e-12)

def test_joint_neighbor_fit_recovers_marginal_target_flux():
    y,x=np.mgrid[-3:4,-3:4];p=np.exp(-(x*x+y*y)/2);p/=p.sum();q=np.exp(-((x-2)**2+y*y)/2);q/=q.sum()
    a=17*p+12*q+3+.1*x
    r=fit_stamp(a,np.full(a.shape,.04),np.ones_like(a,bool),p,x,y,np.array([q]))
    assert r['amplitude']==pytest.approx(17,abs=1e-10)
    assert r['neighbor_amplitudes'][0]==pytest.approx(12,abs=1e-10)
    assert r['variance']>0

def test_native_source_reports_fit_bounds_and_central_quality():
    from spherex_lab.imaging.native_photometry import fit_source
    from spherex_lab.imaging.fits_pipeline import read_cutout
    s=read_cutout('tests/data/sample.fits',allow_source=True)
    x,y=29.1,14.3;position=s['wcs'].celestial.pixel_to_world(x,y)
    fit=fit_source('tests/data/sample.fits',position.ra.deg,position.dec.deg)
    assert fit['fit_bounds']==[26,32,12,17]
    assert fit['central_flag']==int(s['flags'][14,29])
    assert fit['central_valid']==bool(s['valid'][14,29])

def test_neighbor_wing_overlap_normalizes_on_full_source_not_signed_wing():
    from spherex_lab.imaging.native_photometry import neighbor_profile
    raw=np.zeros((9,9));raw[4,4]=1;raw[4,1]=-.02
    yy,xx=np.mgrid[0:3,0:4]
    profile=neighbor_profile(raw,(1,1),xx,yy,6.,1.)
    assert profile.sum()<0
    assert profile[1,3]==pytest.approx(-.02/.98)
    assert np.isfinite(profile).all()

def test_duplicate_and_zero_nuisance_do_not_destroy_identifiable_target():
    y,x=np.mgrid[-3:4,-3:4];p=np.exp(-(x*x+y*y)/2);p/=p.sum();q=np.exp(-((x-2)**2+y*y)/2);q/=q.sum()
    a=17*p+12*q+3+.1*x;v=np.full(a.shape,.04);valid=np.ones_like(a,bool)
    single=fit_stamp(a,v,valid,p,x,y,np.array([q]))
    redundant=fit_stamp(a,v,valid,p,x,y,np.array([q,q,np.zeros_like(q)]))
    assert redundant['amplitude']==pytest.approx(17,abs=1e-10)
    assert redundant['variance']==pytest.approx(single['variance'],rel=1e-12)
    assert redundant['dof']==single['dof']
    assert redundant['nuisance_rank']==4 and redundant['dropped_zero_nuisance_columns']==1

def test_coincident_target_neighbor_is_still_rejected():
    y,x=np.mgrid[-3:4,-3:4];p=np.exp(-(x*x+y*y)/2);p/=p.sum();a=17*p+3
    with pytest.raises(ValueError,match='target amplitude'):
        fit_stamp(a,np.ones_like(a),np.ones_like(a,bool),p,x,y,np.array([p]))
