"""Catalog identities and uncertainty-aware positional candidate associations.

Forced catalog photometry identifies its input catalog row, not an independently
measured astrometric centroid or a moving-object identification.
"""
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.time import Time
from spherex_lab.data.reference_catalog import source_position


def catalog_associations(packets, catalog):
    rows={str(s['source_id']):s for s in catalog['all_sources']}
    associations={}
    for obs_id,packet in packets.items():
        epochs={}
        for ident,measurement in packet.get('measurements',{}).items():
            row=rows.get(str(ident));date=packet.get('date',measurement.get('date'))
            result={'logical_source_id':str(ident),'observation_id':obs_id,'epoch':date,
                    'position_kind':'catalog-propagated forced position; no fitted centroid',
                    'status':'CATALOG_FORCED','association_confidence':'Catalog identity retained; image identity not independently verified',
                    'neighbors_fitted':measurement.get('neighbor_count',0),
                    'blending_status':'JOINT_NEIGHBOR_FIT' if measurement.get('neighbor_count',0) else 'NO_LISTED_NEIGHBOR',
                    'limitations':['No inferred proper-motion detection; omitted parallax and unknown astrometric errors remain limitations.']}
            try:
                if row is None:raise ValueError('Source absent from recorded catalog')
                prediction=source_position(row,date)
                actual=SkyCoord(measurement['ra_deg'],measurement['dec_deg'],unit='deg')
                expected=SkyCoord(prediction['ra_deg'],prediction['dec_deg'],unit='deg')
                separation=float(actual.separation(expected).arcsec)
                result.update(predicted_position=prediction,forced_position_offset_arcsec=separation)
                # Stored positions must agree with their logical catalog identity.
                if separation>1e-5:
                    result.update(status='REJECTED',reason='Recorded forced position does not match propagated catalog identity')
                dt=float(Time(date).jyear)-float(row['ref_epoch'])
                uncertainty=[]
                for position_error,pm_error in [('ra_error','pmra_error'),('dec_error','pmdec_error')]:
                    values=[row.get(position_error),row.get(pm_error)]
                    if any(x in (None,'') for x in values):uncertainty=[];break
                    p,m=map(float,values)
                    if not np.isfinite([p,m]).all() or min(p,m)<0:uncertainty=[];break
                    uncertainty.append((p*p+dt*dt*m*m)/1e6)
                result['position_covariance_arcsec2']=np.diag(uncertainty).tolist() if uncertainty else None
                result['astrometric_uncertainty_policy']='Diagonal propagated catalog approximation; position/PM correlations unavailable' if uncertainty else 'Not supplied by this catalog snapshot; no numerical uncertainty invented'
            except (ValueError,TypeError,KeyError) as exc:
                if result['status']!='REJECTED':result.update(status='UNDERCONSTRAINED',reason=str(exc))
                result['position_covariance_arcsec2']=None
            epochs[str(ident)]=result
        for rejected in packet.get('rejected',[]):
            ident=str(rejected['source_id'])
            epochs.setdefault(ident,{'logical_source_id':ident,'observation_id':obs_id,'status':'REJECTED','reason':rejected['reason']})
        associations[obs_id]=epochs
    return associations


def match_position_candidates(prediction, covariance_arcsec2, candidates, *, chi2_max=9.):
    """Find unique positional candidate given an epoch-specific prediction/covariance.

    Caller must propagate the trajectory to this epoch. Passing this gate only
    yields a candidate, never an astronomical identification or detection.
    """
    covariance=np.asarray(covariance_arcsec2,float)
    if covariance.shape!=(2,2) or not np.isfinite(covariance).all() or not np.allclose(covariance,covariance.T) or np.linalg.eigvalsh(covariance).min()<=0:
        raise ValueError('Positive-definite tangent-plane covariance required')
    if not np.isfinite(chi2_max) or chi2_max<=0:raise ValueError('Invalid association gate')
    if len(candidates)>10000:raise ValueError('Candidate budget exceeded')
    center=SkyCoord(*prediction,unit='deg');accepted=[]
    for c in candidates:
        c_cov=np.asarray(c['covariance_arcsec2'],float)
        if c_cov.shape!=(2,2) or not np.isfinite(c_cov).all() or not np.allclose(c_cov,c_cov.T) or np.linalg.eigvalsh(c_cov).min()<0:
            raise ValueError('Invalid candidate covariance')
        east,north=center.spherical_offsets_to(SkyCoord(c['ra_deg'],c['dec_deg'],unit='deg'))
        residual=np.array([east.arcsec,north.arcsec]);stat=float(residual@np.linalg.solve(covariance+c_cov,residual))
        if stat<=chi2_max:accepted.append({'candidate_id':str(c['id']),'chi2':stat})
    return {'status':'POSITIONAL_CANDIDATE' if len(accepted)==1 else ('AMBIGUOUS' if accepted else 'NO_MATCH'),
            'candidates':accepted,'identification_verified':False,'chi2_max':chi2_max}
