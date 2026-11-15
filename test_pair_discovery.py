import pytest
from pair_discovery import contains, group_epochs, rank_pairs

def row(n,t,det='SPHEREx-D1'):
    return dict(obs_id=n,t_min=t,energy_bandpassname=det,em_min=1e-6,em_max=2e-6,access_url='https://irsa.ipac.caltech.edu/'+n,access_format='image/fits',s_region='POLYGON ICRS 359 -1 1 -1 1 1 359 1')

def test_footprint_wrap():
    assert contains(row('a',0)['s_region'],0,0)
    assert not contains(row('a',0)['s_region'],180,0)

def test_epochs_dedup_anchor():
    r=[row('a',0),row('a',0),row('b',.5),row('c',1.1)]
    g=group_epochs(r,1)
    assert [len(x) for x in g]==[2,1]

def test_rank_excludes_same_epoch_and_detector():
    g=group_epochs([row('a',0),row('b',10),row('c',20,'SPHEREx-D2')],1)
    pairs=rank_pairs(g)
    assert len(pairs)==1
    assert pairs[0]['a']['obs_id']=='a' and pairs[0]['b']['obs_id']=='b'
    assert pairs[0]['status']=='PRELIMINARY'

def test_invalid_threshold():
    with pytest.raises(ValueError): group_epochs([],0)


def test_closed_footprint():
    assert contains(row('a',0)['s_region']+' 359 -1',0,0)
