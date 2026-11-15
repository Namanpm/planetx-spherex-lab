import numpy as np
import pytest

def test_correlated_plane_covariance_matches_dense_oracle():
    from spherex_lab.models.covariance_processing import remove_background
    n=16; x,y=np.meshgrid(np.arange(4),np.arange(4)); design=np.column_stack([np.ones(n),x.ravel(),y.ravel()])
    c=np.eye(n)*2+np.ones((n,n))*.1
    values=design@np.array([2.,.1,-.3]); blank=np.ones(n,dtype=bool)
    out=remove_background(values,c,design,blank)
    h=design@np.linalg.inv(design.T@np.linalg.solve(c,design))@design.T@np.linalg.inv(c)
    np.testing.assert_allclose(out['residual'],0,atol=1e-12)
    np.testing.assert_allclose(out['covariance'],(np.eye(n)-h)@c@(np.eye(n)-h).T,atol=1e-12)

def test_covariance_rejects_missing_blank_and_indefinite():
    from spherex_lab.models.covariance_processing import remove_background
    with pytest.raises(ValueError): remove_background(np.zeros(4),-np.eye(4),np.ones((4,1)),np.ones(4,dtype=bool))

def test_covariance_aware_subtraction():
    from spherex_lab.models.covariance_processing import compare_with_covariance
    a=np.zeros((4,4)); y,x=np.indices(a.shape);b=2+x*.2
    c=np.eye(16)+np.ones((16,16))*.1;valid=np.ones_like(a,dtype=bool)
    r=compare_with_covariance(a,b,c,c,valid,valid,blank_mask=valid,synthetic=True)
    np.testing.assert_allclose(r['difference'],0,atol=1e-12)
    np.testing.assert_allclose(r['variance'].ravel(),np.diag(r['full_covariance']))
