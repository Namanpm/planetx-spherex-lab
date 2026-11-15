import numpy as np
import pytest
from source_association import match_position_candidates,catalog_associations
from reference_catalog import source_position


def candidate(ident,ra):return dict(id=ident,ra_deg=ra,dec_deg=0.,covariance_arcsec2=np.eye(2).tolist())


def test_ra_wrap_unique_and_ambiguous():
    c=candidate('moving',.00001)
    r=match_position_candidates((359.99999,0),np.eye(2),[c])
    assert r['status']=='POSITIONAL_CANDIDATE' and not r['identification_verified']
    assert match_position_candidates((359.99999,0),np.eye(2),[c,candidate('neighbor',.00002)])['status']=='AMBIGUOUS'
    assert match_position_candidates((359.99999,0),np.eye(2),[candidate('far',1.)])['status']=='NO_MATCH'


def test_invalid_covariance_is_not_zero_uncertainty():
    with pytest.raises(ValueError):match_position_candidates((0,0),np.zeros((2,2)),[])


def test_catalog_logical_identity_propagation_and_missing_errors():
    row=dict(source_id='gaia1',ra=20.,dec=30.,pmra=500.,pmdec=-300.,ref_epoch=2016.)
    position=source_position(row,'2026-01-01');measurement=dict(position,flux_jy=1.,variance_jy2=.01)
    packets={'obs':{'date':'2026-01-01','measurements':{'gaia1':measurement},'rejected':[]}}
    a=catalog_associations(packets,{'all_sources':[row]})['obs']['gaia1']
    assert a['status']=='CATALOG_FORCED' and a['position_covariance_arcsec2'] is None
    packets['obs']['measurements']['gaia1']['ra_deg']+=.01
    assert catalog_associations(packets,{'all_sources':[row]})['obs']['gaia1']['status']=='REJECTED'
