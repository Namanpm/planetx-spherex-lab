import pytest
from spherex_lab.data.spacecraft_ephemeris import parse_vector

def test_vector_positive_and_mismatch():
    payload={'result':'Target body name: SPHEREx (-163182)\nSolar System Barycenter ICRF\nOutput units : AU-D\nGeometric state vectors have NO corrections\n$$SOE\n2460000.0, date,1,2,3,0,0,0\n$$EOE'}
    assert parse_vector(payload,2460000)==[1,2,3]
    with pytest.raises(ValueError):parse_vector(payload,2460001)
    with pytest.raises(ValueError):parse_vector({'error':'outside range'},2460000)

def test_wrong_observer_rejected():
    with pytest.raises(ValueError):parse_vector({'result':'Target body name: Earth'},2460000)
