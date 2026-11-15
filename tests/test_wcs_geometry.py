import numpy as np
import pytest
from astropy.wcs import WCS


def grid(dx=0., dy=0.):
    w = WCS(naxis=2)
    w.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    w.wcs.crval = [20., 30.]
    w.wcs.crpix = [3. + dx, 3. + dy]
    w.wcs.cdelt = [-.001, .001]
    return w


def test_identity_and_edges():
    from spherex_lab.imaging.wcs_geometry import reproject_bilinear
    a = np.arange(25.).reshape(5, 5)
    r = reproject_bilinear(a, np.ones_like(a)*4, np.ones_like(a,dtype=bool), grid(), grid(), a.shape)
    assert r['valid'].all()
    np.testing.assert_allclose(r['image'], a, atol=1e-7)
    np.testing.assert_allclose(r['variance'], 4, atol=1e-7)


def test_half_pixel_variance_and_mask():
    from spherex_lab.imaging.wcs_geometry import reproject_bilinear
    a = np.arange(25.).reshape(5,5)
    m = np.ones_like(a, dtype=bool)
    r = reproject_bilinear(a, np.ones_like(a)*4, m, grid(), grid(.5,.5), a.shape)
    assert not r['valid'][0,0]
    assert r['image'][2,2] == pytest.approx((a[1,1]+a[1,2]+a[2,1]+a[2,2])/4, abs=1e-6)
    assert r['variance'][2,2] == pytest.approx(1., abs=1e-6)
    m[1,1] = False
    s = reproject_bilinear(a, np.ones_like(a)*4, m, grid(), grid(.5,.5), a.shape)
    assert not s['valid'][2,2]
    assert np.isnan(s['image'][2,2])


def test_integer_shift_and_no_overlap():
    from spherex_lab.imaging.wcs_geometry import reproject_bilinear
    a = np.zeros((5,5)); a[2,2] = 10
    r = reproject_bilinear(a, np.ones_like(a), np.ones_like(a,dtype=bool), grid(), grid(1,0), a.shape)
    assert r['image'][2,3] == pytest.approx(10)
    r = reproject_bilinear(a, np.ones_like(a), np.ones_like(a,dtype=bool), grid(), grid(100,0), a.shape)
    assert not r['valid'].any()


def test_bad_variance_not_used():
    from spherex_lab.imaging.wcs_geometry import reproject_bilinear
    a=np.ones((5,5)); v=np.ones((5,5)); v[2,2]=-1
    r=reproject_bilinear(a,v,np.ones_like(a,dtype=bool),grid(),grid(),a.shape)
    assert not r['valid'][2,2]
