import pytest
from observation_selection import select_observations,candidate_pairs,validate_coordinates


def records(n):
    return [dict(obs_publisher_did=f'ivo://irsa.ipac/spherex_dr1?obs_{i}/D1',energy_bandpassname='SPHEREx-D1',t_min=60000.+i) for i in range(n)]


@pytest.mark.parametrize('n,pairs',[(2,1),(5,10),(8,28)])
def test_exact_selection_and_all_unordered_pairs(n,pairs):
    r=records(n)
    selected=select_observations(r,[x['obs_publisher_did'] for x in reversed(r)],maximum=n)
    assert selected==r
    assert len(candidate_pairs(selected))==pairs


def test_no_silent_substitution_or_duplicate_identity():
    r=records(5)
    with pytest.raises(ValueError,match='absent'):
        select_observations(r,[r[0]['obs_publisher_did'],'ivo://irsa.ipac/spherex_dr1?missing/D1'])
    with pytest.raises(ValueError,match='distinct'):
        select_observations(r,[r[0]['obs_publisher_did']]*2)
    with pytest.raises(ValueError,match='Conflicting'):
        select_observations(r+[dict(r[0],t_min=1.)])


def test_mixed_detectors_not_silently_filtered():
    r=records(2);r[1]['obs_publisher_did']=r[1]['obs_publisher_did'][:-1]+'2';r[1]['energy_bandpassname']='SPHEREx-D2'
    assert len(select_observations(r))==2
    r[1]['energy_bandpassname']='SPHEREx-D1'
    with pytest.raises(ValueError,match='mismatch'):select_observations(r)


@pytest.mark.parametrize('ra,dec',[(360,0),(0,-91),(None,3),(float('nan'),0)])
def test_coordinates(ra,dec):
    with pytest.raises(ValueError):validate_coordinates(ra,dec)
