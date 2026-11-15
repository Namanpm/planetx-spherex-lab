import numpy as np
import pytest

def test_cross_convolution_identical_effective_psf():
    from psf_processing import cross_match
    ka=np.array([[0,0,0],[.2,.6,.2],[0,0,0.]])
    kb=np.array([[0,.1,0],[.1,.6,.1],[0,.1,0.]])
    a=np.zeros((21,21));b=a.copy();a[9:12,9:12]=ka;b[9:12,9:12]=kb
    r=cross_match(a,np.ones_like(a),b,np.ones_like(b),np.ones_like(a,dtype=bool),np.ones_like(b,dtype=bool),ka,kb)
    np.testing.assert_allclose(r['a'][r['valid']],r['b'][r['valid']],atol=1e-15)
    assert np.allclose(r['variance_a'][r['valid']],(kb**2).sum())

def test_invalid_support_mask():
    from psf_processing import convolve_independent
    a=np.ones((7,7));mask=np.ones_like(a,dtype=bool); mask[3,3]=False
    r=convolve_independent(a,a,mask,np.ones((3,3))/9)
    assert not r['valid'][2:5,2:5].any()
    assert np.isnan(r['image'][~r['valid']]).all()

def test_real_epsf_selection():
    from psf_processing import extract_epsf
    r=extract_epsf('artifacts/sample.fits')
    assert r['kernel'].shape==(7,7)
    assert np.isclose(r['kernel'].sum(),1)
    assert r['provenance']['detector']==1

def test_bad_kernel():
    from psf_processing import convolve_independent
    a=np.ones((3,3))
    with pytest.raises(ValueError): convolve_independent(a,a,a.astype(bool),np.zeros((3,3)))

def test_native_background_fit_covariance():
    from psf_processing import fit_native_background
    yy,xx=np.indices((8,8)); image=2+.3*xx-.2*yy
    blank=np.ones_like(image,dtype=bool); blank[4,4]=False
    image[4,4]+=10
    r=fit_native_background(image,np.ones_like(image),blank)
    np.testing.assert_allclose(r['background'],2+.3*xx-.2*yy,atol=1e-12)
    assert np.isclose(r['residual'][4,4],10)
    assert r['residual_variance'][4,4]>1
    assert r['residual_variance'][1,1]<1
    with pytest.raises(ValueError): fit_native_background(image,np.ones_like(image),blank,independent_native=False)

def test_project_kernel_identity_and_rotation():
    from astropy.wcs import WCS
    from astropy.coordinates import SkyCoord
    from psf_processing import project_kernel
    w=WCS(naxis=2);w.wcs.crpix=[10,10];w.wcs.crval=[10,20];w.wcs.cd=np.array([[-.001,0],[0,.001]]);w.wcs.ctype=['RA---TAN','DEC--TAN']
    sky=SkyCoord(10,20,unit='deg')
    k=np.zeros((5,5));k[2,2]=.6;k[2,3]=.3;k[1,2]=.1
    np.testing.assert_allclose(project_kernel(k,w,w,sky,shape_out=(5,5))['kernel'],k,atol=1e-10)
    r=w.deepcopy();r.wcs.cd=w.wcs.cd@np.array([[0,-1],[1,0]])
    got=project_kernel(k,w,r,sky,shape_out=(5,5))['kernel']
    np.testing.assert_allclose(got,np.rot90(k),atol=1e-10)
