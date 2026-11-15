"""Geometric astrometry with explicit caller-supplied observer vectors.

No spacecraft ephemeris is guessed or silently replaced by Earth's orbit.
This computes geometric directions, not aberration/light-time-corrected apparent
coordinates. An ephemeris adapter must provide ICRS barycentric AU positions at
the requested epoch and preserve authoritative source provenance.
"""
import numpy as np
from astropy.time import Time

AU_KM=149597870.7
YEAR_SECONDS=365.25*86400.
MAS_RAD=np.deg2rad(1./3600000.)


def propagate_direction(target, obstime, *, observer_position_au=None, observer_provenance=None):
    """Propagate supplied ICRS astrometry under rectilinear space motion.

    Radial velocity, when absent, contributes zero to this conditional model and
    is explicitly reported missing. Without distance and observer vector the
    result is angular-only barycentric propagation. Uncertainty propagation,
    light travel time, gravitational deflection and aberration are not included.
    """
    if target.get('frame','ICRS')!='ICRS': raise ValueError('ICRS target required')
    ra=float(target['ra_deg']);dec=float(target['dec_deg'])
    pmra=target.get('pm_ra_cosdec_mas_yr');pmdec=target.get('pm_dec_mas_yr')
    if not np.isfinite([ra,dec]).all() or not 0<=ra<360 or not -90<=dec<=90: raise ValueError('invalid coordinates')
    if pmra is None or pmdec is None or not np.isfinite([pmra,pmdec]).all(): raise ValueError('proper motion unavailable')
    dt=float((Time(obstime).tdb-Time(target['epoch']).tdb).to_value('day')/365.25)
    if not np.isfinite(dt): raise ValueError('invalid epoch')
    a,d=np.deg2rad([ra,dec]);direction=np.array([np.cos(d)*np.cos(a),np.cos(d)*np.sin(a),np.sin(d)])
    east=np.array([-np.sin(a),np.cos(a),0.]);north=np.array([-np.sin(d)*np.cos(a),-np.sin(d)*np.sin(a),np.cos(d)])
    tangent=MAS_RAD*(float(pmra)*east+float(pmdec)*north)
    parallax=target.get('parallax_mas');rv=target.get('radial_velocity_km_s')
    limitations=['Geometric direction only: no light time, aberration or gravitational deflection.',
                 'Astrometric covariance is not propagated.']
    if rv is not None and not np.isfinite(rv): raise ValueError('invalid radial velocity')
    distance_available=parallax is not None and np.isfinite(parallax) and parallax>0
    if distance_available:
        distance_au=1/(float(parallax)*MAS_RAD)
        velocity=distance_au*tangent
        if rv is not None: velocity += direction*float(rv)*YEAR_SECONDS/AU_KM
        else: limitations.append('Radial velocity unavailable; perspective acceleration is conditional on zero radial velocity.')
        position=distance_au*direction+dt*velocity
    else:
        position=direction+dt*tangent
        limitations.append('Distance unavailable: angular-only propagation; perspective motion not modeled.')
    if observer_position_au is not None:
        if not distance_available: raise ValueError('observer parallax requires trustworthy distance')
        observer=np.asarray(observer_position_au,dtype=float)
        if observer.shape!=(3,) or not np.isfinite(observer).all() or not observer_provenance:
            raise ValueError('finite ICRS barycentric observer vector and provenance required')
        position=position-observer
        kind='observer-relative geometric direction'
    else:
        if observer_provenance is not None: raise ValueError('observer provenance supplied without vector')
        limitations.append('No observer vector supplied; no observer parallax applied.')
        kind='barycentric geometric direction'
    norm=np.linalg.norm(position)
    if not np.isfinite(norm) or norm<=0: raise ValueError('invalid propagated position')
    unit=position/norm
    return {'ra_deg':float(np.rad2deg(np.arctan2(unit[1],unit[0]))%360),
            'dec_deg':float(np.rad2deg(np.arcsin(np.clip(unit[2],-1,1)))),
            'epoch':Time(obstime).isot,'frame':'ICRS','coordinate_kind':kind,
            'observer_provenance':observer_provenance,'limitations':limitations}
