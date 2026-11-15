import numpy as np
import pytest
from spherex_lab.repeatability.empirical_repeatability import fit_repeatability


def test_training_and_heldout_known_scale():
    ids=np.arange(12).astype(str);a=np.linspace(100,300,12);b=a*1.03
    result=fit_repeatability(ids,a,b,np.ones(12),np.ones(12),training_ids=ids[:8],heldout_ids=ids[8:],outlier_sigma=4)
    np.testing.assert_allclose(result['scale_at_reference_color'],1.03,atol=1e-10)
    assert result['heldout']['source_count']==4
    np.testing.assert_allclose(result['heldout']['residual_log_ratio'],0,atol=1e-10)
    assert result['excess_log_scatter']==0


def test_no_leakage_from_heldout_variable():
    ids=np.arange(16).astype(str);a=np.full(16,100.);b=a.copy();b[-1]=200
    result=fit_repeatability(ids,a,b,np.ones(16),np.ones(16),training_ids=ids[:12],heldout_ids=ids[12:],outlier_sigma=4)
    np.testing.assert_allclose(result['scale_at_reference_color'],1.,atol=1e-12)
    assert max(result['heldout']['standardized_residual'])>10
    assert result['heldout']['coverage_95_fraction']<1


def test_overlap_ids_and_missing_sources_rejected():
    args=(['a','b','c','d'],[100]*4,[100]*4,[1]*4,[1]*4)
    with pytest.raises(ValueError,match='disjoint'):fit_repeatability(*args,training_ids=['a','b','c'],heldout_ids=['c'],outlier_sigma=4)
    with pytest.raises(ValueError):fit_repeatability(*args,training_ids=['a','b','missing'],heldout_ids=['d'],outlier_sigma=4)


def test_color_fit_and_positive_scatter():
    rng=np.random.default_rng(82);n=80;ids=np.arange(n).astype(str);color=np.linspace(0,2,n);a=np.full(n,1000.);b=a*np.exp(.03+.04*(color-1)+rng.normal(0,.02,n))
    r=fit_repeatability(ids,a,b,np.ones(n),np.ones(n),training_ids=ids[:60],heldout_ids=ids[60:],color=color,outlier_sigma=5)
    assert abs(r['coefficients'][1]-.04)<.01
    assert .01<r['excess_log_scatter']<.04
    assert r['heldout']['extrapolated_color_source_count']>0


def test_robust_training_rejection_is_explicit():
    n=30;ids=np.arange(n).astype(str);a=np.full(n,100.);b=a.copy();b[0]=200
    r=fit_repeatability(ids,a,b,np.ones(n),np.ones(n),training_ids=ids[:25],heldout_ids=ids[25:],outlier_sigma=4)
    assert '0' in r['rejected_training_ids']
    assert r['scale_at_reference_color']==pytest.approx(1)

def test_shared_correction_induces_heldout_covariance():
    ids=np.arange(12).astype(str);a=np.full(12,100.);b=np.full(12,101.)
    r=fit_repeatability(ids,a,b,np.ones(12),np.ones(12),training_ids=ids[:8],heldout_ids=ids[8:],outlier_sigma=4)
    covariance=np.array(r['heldout']['predictive_covariance'])
    assert covariance[0,1]>0
    assert r['heldout']['dof']==4
    assert 0<=r['heldout']['pvalue_chi2_approximate']<=1

def test_spatial_residual_diagnostic_detects_controlled_gradient():
    from spherex_lab.repeatability.empirical_repeatability import spatial_repeatability_diagnostic
    y,x=np.indices((8,8));xy=np.column_stack([x.ravel(),y.ravel()])
    gradient=x.ravel().astype(float)
    r=spatial_repeatability_diagnostic(xy,gradient,seed=97,permutations=199)
    assert r['moran_I']>.7
    assert r['pvalue_permutation_positive']<=.01
    assert r['status']=='DIAGNOSTIC_ONLY'
    repeat=spatial_repeatability_diagnostic(xy,gradient,seed=97,permutations=199)
    assert r==repeat


def test_spatial_diagnostic_degenerate_and_duplicate():
    from spherex_lab.repeatability.empirical_repeatability import spatial_repeatability_diagnostic
    xy=np.array([[0,0],[0,1],[1,0],[1,1]])
    r=spatial_repeatability_diagnostic(xy,np.ones(4),seed=5,permutations=19)
    assert r['moran_I'] is None
    with pytest.raises(ValueError):spatial_repeatability_diagnostic(np.zeros((4,2)),np.arange(4),seed=5,permutations=19)
