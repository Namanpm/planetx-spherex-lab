import numpy as np
import pytest
from spherex_lab.models.spectral_model import response_design, fit_spectrum, predict_spectrum
from spherex_lab.models.motion_inference import fit_motion, fit_variability


def test_response_integrals_and_identifiable_fit():
    wavelength=np.linspace(1.,3.,101)
    basis=np.column_stack([np.ones_like(wavelength),wavelength-2])
    responses=np.array([np.exp(-((wavelength-c)/.25)**2) for c in (1.2,1.6,2.2,2.8)])
    design=response_design(wavelength,responses,basis)
    truth=np.array([3.,.4]); variance=np.full(4,.01)
    fit=fit_spectrum(design,design@truth,variance)
    np.testing.assert_allclose(fit['coefficients'],truth,atol=1e-12)
    expected=np.linalg.inv(design.T@(design/variance[:,None]))
    np.testing.assert_allclose(fit['covariance'],expected)
    prediction=predict_spectrum(fit,design[:1])
    np.testing.assert_allclose(prediction['mean'],design[:1]@truth)
    assert prediction['variance'][0]>0


def test_unidentifiable_rejected_even_with_ridge():
    for penalty in (0,10):
        with pytest.raises(ValueError,match='identifiable'):
            fit_spectrum(np.ones((4,2)),np.ones(4),np.ones(4),ridge=penalty)


def test_heldout_model_prediction():
    x=np.array([[1.,-1.],[1.,0.],[1.,1.],[1.,2.]])
    fit=fit_spectrum(x[:3],x[:3]@np.array([2.,3.]),np.ones(3))
    pred=predict_spectrum(fit,x[3:],observed=[8.],measurement_variance=[1.])
    np.testing.assert_allclose(pred['residual'],0,atol=1e-12)
    assert pred['chi2']<1e-20


def test_invalid_spectral_inputs():
    with pytest.raises(ValueError): response_design([1,np.inf],[[1,1]],[[1],[1]])
    with pytest.raises(ValueError): fit_spectrum(np.eye(2),[1,1],[1,0])
    with pytest.raises(ValueError): fit_spectrum(np.eye(2),[1,1],[1,1],ridge=-1)


def test_motion_truth_and_static_comparison():
    times=59000+np.arange(5)*365.25
    x=1.2+.15*np.arange(5);y=-.7-.2*np.arange(5)
    result=fit_motion(times,x,y,np.full(5,.01),np.full(5,.02),reference_mjd=59000)
    np.testing.assert_allclose(result['velocity_arcsec_per_year'],[.15,-.2],atol=1e-12)
    np.testing.assert_allclose(result['position_arcsec'],[1.2,-.7],atol=1e-12)
    assert result['delta_chi2_static_minus_motion']>100
    static=fit_motion(times,np.ones(5),np.ones(5),np.ones(5),np.ones(5))
    assert abs(static['delta_chi2_static_minus_motion'])<1e-12


def test_motion_needs_distinct_times_and_three_epochs():
    with pytest.raises(ValueError): fit_motion([1,1,1],[1]*3,[1]*3,[1]*3,[1]*3)
    with pytest.raises(ValueError): fit_motion([1,2],[1]*2,[1]*2,[1]*2,[1]*2)


def test_variable_flux_constant_model():
    r=fit_variability([3,3,3],[.25,.25,.25])
    assert r['chi2_constant']==0 and r['constant_flux']==3
    r=fit_variability([1,3,5],[.25,.25,.25])
    assert r['chi2_constant']==32 and r['dof']==2

def test_association_rejects_competing_sources():
    from spherex_lab.models.motion_inference import associate_epochs
    r=associate_epochs([[0,0]],[[.1,.1]],[[.2,0],[-.2,0]],[[.1,.1],[.1,.1]],max_separation_arcsec=1.)
    assert r['matches']==[] and r['ambiguous_a']==[0]
    assert r['unmatched_a']==[0]


def test_association_uncertainty_and_static_position():
    from spherex_lab.models.motion_inference import associate_epochs
    r=associate_epochs([[0,0],[10,0]],[[.1,.1],[.1,.1]],[[.2,0],[10,0]],[[.1,.1],[.1,.1]],max_separation_arcsec=1.)
    assert [(m['index_a'],m['index_b']) for m in r['matches']]==[(0,0),(1,1)]
    np.testing.assert_allclose(r['matches'][0]['position_chi2'],2.)
    assert r['matches'][1]['position_chi2']==0

def test_spectral_fit_stable_without_normal_equation_squaring():
    # Dense rotated design with known singular values and coefficients. Forming
    # X.T X loses the weakest direction even though X is within fit policy.
    rng=np.random.default_rng(9834)
    left,_=np.linalg.qr(rng.normal(size=(8,2)))
    right=np.array([[1.,1.],[-1.,1.]])/np.sqrt(2.)
    singular=np.array([1.,1e-8])
    design=left@np.diag(singular)@right.T
    truth=np.array([1.25,-.75])
    result=fit_spectrum(design,design@truth,np.ones(8))
    np.testing.assert_allclose(result['coefficients'],truth,atol=2e-8,rtol=0)
    expected=(right/singular)@(right/singular).T
    np.testing.assert_allclose(result['covariance'],expected,rtol=1e-8,atol=1.)


def test_repeated_two_epoch_measurements_not_three_distinct_epochs():
    with pytest.raises(ValueError,match='distinct'):
        fit_motion([1.,1.,2.,2.],[0.,0.,1.,1.],[0.]*4,[.1]*4,[.1]*4)
