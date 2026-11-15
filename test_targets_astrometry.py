import pytest
from targets import targets, barycentric_position

def test_wise_simbad_astrometry():
    wise = next(t for t in targets() if t['id'] == 'wise-0855')
    assert wise['ra_deg'] == pytest.approx(133.81890825)
    assert wise['dec_deg'] == pytest.approx(-7.247048, abs=1e-7)
    assert wise['pm_ra_cosdec_mas_yr'] == -8123.7
    assert wise['parallax_mas'] == 439.0

def test_proper_motion_propagation_reference_epoch():
    target = next(t for t in targets() if t['id'] == 'wise-0855')
    result = barycentric_position('wise-0855', 'J2000')
    assert result['ra_deg'] == pytest.approx(target['ra_deg'], abs=1e-9)
    assert result['dec_deg'] == pytest.approx(target['dec_deg'], abs=1e-9)
    assert 'parallax' in result['limitations']

