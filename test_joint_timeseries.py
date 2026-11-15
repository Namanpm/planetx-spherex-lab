"""Synthetic joint-scale checks; these do not establish real SPHEREx validity."""
import copy
import numpy as np
import pytest
from joint_timeseries import analyze_time_series,_partition


def synthetic(seed=42, variable=False, epochs=4):
    rng=np.random.default_rng(seed)
    train=[str(i) for i in range(80) if _partition(str(i))<2]
    color={str(i):.6+.8*rng.random() for i in range(80)}
    color[train[0]]=0.;color[train[1]]=2.
    packets={}
    for e in range(epochs):
        packet={'date':f'2025-01-{e+1:02d}T00:00:00','measurements':{},
                'common_scale_domain':{'id':'synthetic-linear-color-domain', 'validated':True,
                    'provenance':'Synthetic generator explicitly uses this log-linear response model; no real instrument claim.'}}
        for i in range(80):
            s=str(i);base=10.+i/10;offset=.03*e+.02*e*(color[s]-1)
            flux=base*np.exp(offset+rng.normal(0,.01))
            packet['measurements'][s]={'flux_jy':flux,'variance_jy2':flux**2*.01**2,
                'color':color[s],'fit_bounds':[i*10,i*10+4,0,4], 'central_valid':True,'central_flag':0}
        offset=.03*e+(.3 if variable and e==2 else 0)
        flux=12*np.exp(offset+rng.normal(0,.01))
        packet['measurements']['interest']={'flux_jy':flux,'variance_jy2':flux**2*.01**2,
             'color':1.,'fit_bounds':[900,904,0,4]}
        packets[f'obs-{e}']=packet
    return packets


def test_simultaneous_calibration_and_independent_null_timeline():
    result=analyze_time_series(synthetic(),'interest')
    assert result['status']=='SUCCESS',result['reason']
    assert result['confirmation_tested_once']
    assert all('interest' not in ids for ids in result['partition_source_ids'].values())
    timeline=result['timelines'][0]
    assert timeline['constant_flux_assessment']['pvalue_chi2_approximate']>.05
    cov=np.asarray(timeline['shared_calibration_covariance_jy2'])
    assert cov[1,2]!=0 and np.allclose(cov,cov.T)
    assert np.linalg.eigvalsh(np.asarray(timeline['flux_covariance_jy2'])).min()>0
    assert result['joint_calibration']['rank_identified']


def test_injected_variability_is_detectable_conditional_on_synthetic_model():
    result=analyze_time_series(synthetic(variable=True),'interest')
    assert result['status']=='SUCCESS',result['reason']
    assert result['timelines'][0]['constant_flux_assessment']['pvalue_chi2_approximate']<1e-10
    assert 'automatic astrophysical detection' in result['timelines'][0]['constant_flux_assessment']['interpretation']


def test_unverified_common_domain_does_not_promote_pairwise_evidence():
    data=synthetic()
    for p in data.values():p.pop('common_scale_domain')
    r=analyze_time_series(data,'interest')
    assert r['status']=='INSUFFICIENT_EVIDENCE'
    assert 'confirmation' in r and not r['timelines']
    assert 'domain' in r['reason']


def test_missing_epoch_retained_and_three_distinct_epochs_required():
    data=synthetic(epochs=4);data['obs-3']['measurements'].pop('interest')
    r=analyze_time_series(data,'interest')
    assert r['status']=='SUCCESS',r['reason']
    assert len(r['timelines'][0]['points'])==3
    assert len(r['timelines'][0]['missing_epochs'])==1
    data['obs-2']['measurements'].pop('interest')
    r=analyze_time_series(data,'interest')
    assert r['status']=='INSUFFICIENT_EVIDENCE' and 'fewer than three' in r['reason']
    assert analyze_time_series(dict(list(data.items())[:2]),'interest')['status']=='INSUFFICIENT_EVIDENCE'


def test_joint_rank_deficiency_rejected_without_regularization():
    data=synthetic()
    for p in data.values():
        for m in p['measurements'].values():m['color']=1.
    r=analyze_time_series(data,'interest')
    assert r['status']=='INSUFFICIENT_EVIDENCE'
    assert 'not identifiable' in r['reason']


def test_blended_reference_stamps_cannot_leak_across_partitions():
    data=synthetic()
    for p in data.values():
        for s,m in p['measurements'].items():
            if s!='interest':m['fit_bounds']=[0,4,0,4]
    r=analyze_time_series(data,'interest')
    assert r['status']=='INSUFFICIENT_EVIDENCE'
    assert 'disjoint' in r['reason']
    assert r['excluded_references']


def test_ambiguous_target_and_color_extrapolation_are_missing():
    data=synthetic()
    data['obs-2']['measurements']['interest']['association_status']='AMBIGUOUS'
    data['obs-3']['measurements']['interest']['color']=99.
    r=analyze_time_series(data,'interest')
    assert r['status']=='INSUFFICIENT_EVIDENCE'
    assert len(r['missing_epochs'])==2


def test_confirmation_does_not_tune_coefficients_or_choose_new_epochs():
    data=synthetic();first=analyze_time_series(data,'interest')
    for p in data.values():
        for s,m in p['measurements'].items():
            if s!='interest' and _partition(s)==3 and p is data['obs-2']:
                m['flux_jy']*=1.3
    second=analyze_time_series(data,'interest')
    assert second['status']=='INSUFFICIENT_EVIDENCE'
    assert 'confirmation rejected' in second['reason']
    assert second['joint_calibration']==first['joint_calibration']
    assert second['observation_ids']==first['observation_ids']


def test_no_false_trend_and_approximate_uncertainty_coverage_repeated_nulls():
    # Fixed-seed controlled experiment; retain genuine calibration rejections.
    successes=[];covered=[]
    for seed in range(10,20):
        result=analyze_time_series(synthetic(seed=seed),'interest')
        if result['status']!='SUCCESS':continue
        t=result['timelines'][0]
        successes.append(t['constant_flux_assessment']['pvalue_chi2_approximate']>.05)
        values=np.array([p['calibrated_flux_jy'] for p in t['points']])
        sigma=np.sqrt(np.diag(t['flux_covariance_jy2']))
        covered.extend(abs(values-12)<=1.96*sigma)
    assert len(successes)>=6
    assert np.mean(successes)>=.8
    assert np.mean(covered)>=.8

def test_common_response_applicability_is_not_transitive():
    data=synthetic()
    data['obs-2']['common_scale_domain']['id']='different-independent-domain'
    r=analyze_time_series(data,'interest')
    assert r['status']=='INSUFFICIENT_EVIDENCE'
    assert 'joint_calibration' in r and not r['timelines']


def test_failed_observation_preserved_as_missing_in_valid_subset():
    data=synthetic()
    data['obs-3']['error']={'type':'RetryError','reason':'502 exhausted'}
    r=analyze_time_series(data,'interest')
    assert r['status']=='SUCCESS',r['reason']
    assert len(r['timelines'][0]['points'])==3
    assert r['timelines'][0]['missing_epochs'][0]['observation_id']=='obs-3'
    assert r['rejected_observations'][0]['observation_id']=='obs-3'
