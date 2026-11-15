import numpy as np
import pytest
from astrometry_model import propagate_direction


def target():
    return {'ra_deg':0.,'dec_deg':0.,'epoch':'J2000','pm_ra_cosdec_mas_yr':0.,'pm_dec_mas_yr':0.,'parallax_mas':1000.}

def test_geometric_one_parsec_parallax():
    t=target()
    result=propagate_direction(t,'J2000',observer_position_au=[0,1,0],observer_provenance='Synthetic 1 AU transverse observer')
    # One AU transverse displacement at one parsec subtends one arcsecond.
    ra=((result['ra_deg']+180)%360)-180
    np.testing.assert_allclose(ra*3600,-1.,rtol=1e-10)
    assert result['observer_provenance'].startswith('Synthetic')

def test_missing_geometry_is_explicit_and_not_earth_default():
    result=propagate_direction(target(),'J2000')
    assert result['observer_provenance'] is None
    assert result['coordinate_kind']=='barycentric geometric direction'
    with pytest.raises(ValueError): propagate_direction(target(),'J2000',observer_position_au=[0,1,0])

def test_motion_and_unavailable_distance():
    t=target();t['pm_ra_cosdec_mas_yr']=1000.
    result=propagate_direction(t,'J2001')
    np.testing.assert_allclose(result['ra_deg']*3600,1.,rtol=1e-10)
    t['parallax_mas']=None
    with pytest.raises(ValueError,match='distance'): propagate_direction(t,'J2000',observer_position_au=[0,1,0],observer_provenance='synthetic')
