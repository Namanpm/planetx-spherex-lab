import numpy as np
import pytest
from spherex_lab.repeatability.empirical_comparison import stamps_overlap,fit_pair_model

def test_fit_stamps_are_independent_only_when_disjoint():
    assert stamps_overlap([0,6,0,6],[6,12,0,6])
    assert not stamps_overlap([0,6,0,6],[7,13,0,6])

def test_reference_model_uses_source_level_holdout():
    ids=[str(i) for i in range(40)];rng=np.random.default_rng(21)
    colors=np.linspace(.2,2,40);fa=np.full(40,10.);fb=fa*np.exp(.03+.01*(colors-1))+rng.normal(0,.03,40)
    a={i:{'flux_jy':fa[j],'variance_jy2':.001,'color':colors[j]} for j,i in enumerate(ids)}
    b={i:{'flux_jy':fb[j],'variance_jy2':.001,'color':colors[j]} for j,i in enumerate(ids)}
    fit=fit_pair_model(ids,a,b,ids[:20],ids[20:])
    assert fit['training_source_count']>=12 and fit['heldout']['source_count']==20
    assert set(fit['retained_training_ids']).isdisjoint(fit['heldout']['source_ids'])

def test_actual_acceptance_key_and_extrapolation():
    from spherex_lab.repeatability.empirical_comparison import is_consistent
    assert is_consistent({'heldout':{'pvalue_chi2_approximate':.5,'extrapolated_color_source_count':0}})
    assert not is_consistent({'heldout':{'pvalue_chi2_approximate':.5,'extrapolated_color_source_count':1}})
    assert not is_consistent({'heldout':{'pvalue_chi2_approximate':.001,'extrapolated_color_source_count':0}})

def test_repeatability_endpoint_rejects_invalid_without_network():
    from main import app
    from fastapi.testclient import TestClient
    c=TestClient(app)
    assert c.post('/api/spherex/repeatability',json={'target':'unknown','max_observations':2}).status_code==422
    assert c.post('/api/spherex/repeatability',json={'target':'rw-aurigae','max_observations':99}).status_code==422
