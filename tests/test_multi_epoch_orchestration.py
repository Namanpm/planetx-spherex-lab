"""Offline orchestration regressions for input gates and replay evidence."""
import copy
import json
from pathlib import Path
import pytest
from spherex_lab.repeatability import multi_epoch
def _no_network(*args, **kwargs):
    raise AssertionError('Invalid request or explicit replay must not query NASA')


@pytest.mark.parametrize('kwargs', [{'mode':'invalid'}, {'maximum':1}, {'maximum':99}])
def test_invalid_configuration_rejected_before_discovery(tmp_path, monkeypatch, kwargs):
    from spherex_lab.data import pair_discovery
    monkeypatch.setattr(pair_discovery,'discover_pairs',_no_network)
    with pytest.raises(ValueError):
        multi_epoch.run_analysis(output=tmp_path/'run',**kwargs)


def test_readiness_does_not_accept_missing_calibration_validation(tmp_path):
    science=tmp_path/'science.fits';science.write_bytes(b'not used by readiness summary')
    packet={'record':{'obs_publisher_did':'ivo://irsa.ipac/spherex_dr1?x/D1','t_min':60000.},
            'science_path':str(science),'measurements':{'source':{'flux_jy':1.}},
            'calibration':{'calibration_url':'https://irsa.ipac.caltech.edu/ibe/data/spherex/dr1/spectral_wcs/missing.fits'},
            '_readiness_checks':{'wcs_valid':True,'epsf_available':True,'variance_available':True,'quality_flags_available':True}}
    state=multi_epoch.readiness(packet)
    assert state['wavelength_map_valid'] is not True
    assert state['status']!='READY'


def _replay_inputs_present():
    """The replay snapshot needs ~150 MB of locally cached real cutouts/calibrations."""
    report=Path('artifacts/rw_empirical/report.json')
    if not report.is_file():
        return False
    observations=json.loads(report.read_text(encoding='utf-8'))['observations']
    return all(Path(o['measurement_file']).is_file() for o in observations)


requires_replay=pytest.mark.skipif(not _replay_inputs_present(),
    reason='Local RW Aurigae replay snapshot (measurement packets, cached science/calibration FITS) is not '
           'present; regenerate with: python repeatability_products.py --target rw-aurigae --output artifacts/rw_empirical')


def _two_packet_snapshot(tmp_path):
    baseline=Path('artifacts/rw_empirical/report.json')
    data=json.loads(baseline.read_text())
    data['observations']=copy.deepcopy(data['observations'][:2])
    for j,o in enumerate(data['observations']):
        packet=json.loads(Path(o['measurement_file']).read_text())
        packet_path=tmp_path/f'packet-{j}.json'
        packet_path.write_text(json.dumps(packet));o['measurement_file']=str(packet_path)
    report=tmp_path/'report.json';report.write_text(json.dumps(data))
    return report,data


@requires_replay
def test_snapshot_rejects_science_detector_observation_identity_mismatch(tmp_path,monkeypatch):
    from spherex_lab.data import calibration_resolver
    monkeypatch.setattr(calibration_resolver,'discover',_no_network)
    report,data=_two_packet_snapshot(tmp_path)
    a_path=Path(data['observations'][0]['measurement_file']);b_path=Path(data['observations'][1]['measurement_file'])
    a=json.loads(a_path.read_text());b=json.loads(b_path.read_text())
    a['science_path']=b['science_path'];a_path.write_text(json.dumps(a))
    _,_,packets,_=multi_epoch.snapshot(report,maximum=2)
    failed=packets[data['observations'][0]['publisher_id']]
    assert failed.get('error')
    assert not failed['measurements']
    assert 'identity' in failed['error']['reason'] or 'observation' in failed['error']['reason']


@requires_replay
def test_snapshot_missing_calibration_is_explicit_in_readiness(tmp_path,monkeypatch):
    from spherex_lab.data import calibration_resolver
    monkeypatch.setattr(calibration_resolver,'discover',_no_network)
    report,data=_two_packet_snapshot(tmp_path)
    _,_,packets,_=multi_epoch.snapshot(report,maximum=2)
    for o in data['observations']:
        state=multi_epoch.readiness(packets[o['publisher_id']])
        assert state['status']!='READY'
        assert state['wavelength_map_valid'] is not True


def test_replay_rejects_wrong_target_without_discovery(tmp_path,monkeypatch):
    from spherex_lab.data import pair_discovery
    monkeypatch.setattr(pair_discovery,'discover_pairs',_no_network)
    monkeypatch.setattr(multi_epoch,'snapshot',lambda *a,**k:({'target':'different-target'},[],{},'hash'))
    with pytest.raises(ValueError,match='target'):
        multi_epoch.run_analysis(target_id='rw-aurigae',selection_mode='replay',
                                 replay_report=tmp_path/'snapshot.json',output=tmp_path/'run')

@requires_replay
def test_valid_real_snapshot_replay_preserves_selected_pair_offline(tmp_path,monkeypatch):
    from spherex_lab.data import pair_discovery
    from spherex_lab.data import calibration_resolver
    monkeypatch.setattr(pair_discovery,'discover_pairs',_no_network)
    monkeypatch.setattr(calibration_resolver,'discover',_no_network)
    selected=['ivo://irsa.ipac/spherex_dr1?2025W36_1A_0446_1/D1',
              'ivo://irsa.ipac/spherex_dr1?2026W10_1A_0445_3/D1']
    result=multi_epoch.run_analysis(selection_mode='replay',replay_report=Path('artifacts/rw_empirical/report.json'),
          selected_ids=selected,maximum=2,mode='best_pair',output=tmp_path/'run')
    assert result['provenance']['selection_mode']=='replay'
    assert result['provenance']['measurement_replay']
    assert result['status']=='SUCCESS'
    assert all(s['status']=='READY' for s in result['readiness'])
    assert len(result['theoretical_pairs'])==len(result['pair_evaluations'])==1
    assert result['best_pair']['selected_pair']==selected
    assert result['best_pair_confirmation_status']=='REFERENCE_MODEL_CONSISTENT_ON_CONFIRMATION'
    assert result['provenance']['input_sha256']
