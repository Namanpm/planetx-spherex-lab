import pytest
import empirical_comparison as engine


def fixtures(n):
    records=[{'obs_publisher_did':f'ivo://irsa.ipac/spherex_dr1?synthetic_{j}/D1','t_min':60000.+j/2,'obs_id':f'synthetic_{j}','energy_bandpassname':'SPHEREx-D1'} for j in range(n)]
    packets={r['obs_publisher_did']:{'record':r,'measurements':{},'rejected':[],'error':{'type':'Synthetic','reason':'deliberate fixture missing calibration'}} for r in records}
    catalog={'sources':[],'all_sources':[],'provenance':{'sha256':'synthetic-fixture'}}
    return records,packets,catalog


@pytest.mark.parametrize('n,expected',[(2,1),(5,10),(8,28)])
def test_replay_enumerates_all_pairs_without_network(tmp_path,monkeypatch,n,expected):
    records,packets,catalog=fixtures(n)
    def forbidden(*a,**k):raise AssertionError('offline replay attempted I/O')
    for name in ['discover_pairs','query_references','download_cutout','wavelength_maps','fit_source']:monkeypatch.setattr(engine,name,forbidden)
    monkeypatch.setattr(engine,'render_sources',lambda *a,**k:None)
    progress=[]
    report=engine.run_selected(records,ra=12.,dec=0.,output=tmp_path,catalog=catalog,premeasured_packets=packets,progress=progress.append)
    assert report['measurement_replay']
    assert report['theoretical_pair_count']==expected
    assert len(report['pair_evaluations'])==expected
    assert len({tuple(sorted([p['observation_a'],p['observation_b']])) for p in report['pair_evaluations']})==expected
    assert all(p['status']=='REJECTED' and p['reason'] for p in report['pair_evaluations'])
    assert len(report['observations'])==n
    assert progress[-1]['stage']=='complete'


def test_replay_identity_mismatch_rejected(tmp_path):
    records,packets,catalog=fixtures(2)
    packets[records[0]['obs_publisher_did']]['record']={**records[0],'obs_publisher_did':'wrong'}
    with pytest.raises(ValueError,match='identity mismatch'):
        engine.run_selected(records,ra=12.,dec=0.,output=tmp_path,catalog=catalog,premeasured_packets=packets)


def test_cached_failed_packets_are_retried(tmp_path,monkeypatch):
    records,packets,catalog=fixtures(2)
    calls=[]
    def fail(record,*a,**kw):calls.append(record['obs_publisher_did']);raise OSError('simulated upstream failure')
    monkeypatch.setattr(engine,'download_cutout',fail)
    monkeypatch.setattr(engine,'render_sources',lambda *a,**k:None)
    for _ in range(2):engine.run_selected(records,ra=12.,dec=0.,output=tmp_path,catalog=catalog)
    assert len(calls)==4


def test_different_detectors_are_explicit_rejection(tmp_path,monkeypatch):
    records,packets,catalog=fixtures(2)
    records[1]['energy_bandpassname']='SPHEREx-D2'
    monkeypatch.setattr(engine,'render_sources',lambda *a,**k:None)
    report=engine.run_selected(records,ra=12.,dec=0.,output=tmp_path,catalog=catalog,premeasured_packets=packets)
    assert 'different detector' in report['pair_evaluations'][0]['reason']
