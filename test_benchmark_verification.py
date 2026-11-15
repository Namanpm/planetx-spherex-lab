"""Regression verifier must catch science changes, not just source counts."""
import copy
import json
from benchmark_verification import compare_reports, peak_memory_bytes


def reports(tmp_path):
    packet={'measurements':{'source':{'flux_jy':1.,'variance_jy2':.01,'flag_values':[0,2097152]}},'rejected':[{'source_id':'other','reason':'quality'}]}
    a=tmp_path/'a.json'; b=tmp_path/'b.json'
    a.write_text(json.dumps(packet));b.write_text(json.dumps(packet))
    report={'observations':[{'publisher_id':'obs','measurement_file':str(a)}], 'confirmation':{'selected_pair':['a','b'],'model':{'coefficients':[.1,.2],'heldout':{'chi2':4.8}}},'pair_evaluations':[],'status':'OK'}
    actual=copy.deepcopy(report);actual['observations'][0]['measurement_file']=str(b)
    return report,actual,b


def test_verifier_checks_flux_uncertainty_flags_model_and_confirmation(tmp_path):
    baseline,actual,path=reports(tmp_path)
    assert compare_reports(baseline,actual)['passed']
    packet=json.loads(path.read_text());packet['measurements']['source']['variance_jy2']*=2;packet['measurements']['source']['flag_values']=[1]
    path.write_text(json.dumps(packet));actual['confirmation']['model']['coefficients'][0]+=.01
    result=compare_reports(baseline,actual)
    assert not result['passed']
    assert any('variance_jy2' in p for p in result['differences'])
    assert any('flag_values' in p for p in result['differences'])
    assert any('coefficients' in p for p in result['differences'])
    assert not result['confirmation_exactly_identical']


def test_process_peak_includes_resident_memory():
    assert peak_memory_bytes()>1024*1024
