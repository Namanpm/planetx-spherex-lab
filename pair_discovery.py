"""Preliminary pair discovery. Ranking never establishes spectral equivalence."""
import re
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.time import Time
from spherex import search
from targets import targets, barycentric_position


def contains(region, ra, dec):
    """Containment for small convex ICRS spherical polygon footprints only."""
    parts = str(region).split()
    if parts[:2] != ['POLYGON', 'ICRS'] or len(parts)<8 or len(parts)%2:
        raise ValueError('unsupported archive footprint')
    xy=np.asarray(parts[2:],dtype=float).reshape(-1,2)
    if not np.isfinite(xy).all(): raise ValueError('invalid footprint')
    if np.array_equal(xy[0],xy[-1]): xy=xy[:-1]
    if len(xy)<3: raise ValueError('degenerate footprint')
    vertices=SkyCoord(xy[:,0],xy[:,1],unit='deg').cartesian.xyz.value.T
    center=vertices.sum(axis=0); center /= np.linalg.norm(center)
    normals=np.cross(vertices,np.roll(vertices,-1,axis=0))
    orientation=normals@center
    if np.any(np.abs(orientation)<1e-12): raise ValueError('degenerate footprint')
    if np.any((normals@vertices.T)*np.sign(orientation)[:,None] < -1e-12):
        raise ValueError('nonconvex footprint')
    point=SkyCoord(ra,dec,unit='deg').cartesian.xyz.value
    return bool(np.all((normals@point)*np.sign(orientation)>=-1e-12))


def group_epochs(records, span_days=1):
    """Groups have max span <= span_days from first start, never chained gaps."""
    if not np.isfinite(span_days) or span_days<=0: raise ValueError('invalid epoch span')
    unique={}
    for r in records:
        if r.get('t_min') is None or not np.isfinite(r['t_min']): continue
        if not r.get('access_url') or r.get('access_format')!='image/fits': continue
        key=r.get('obs_publisher_did') or (r.get('obs_id'),r.get('energy_bandpassname'))
        unique.setdefault(key,r)
    groups=[]
    for r in sorted(unique.values(),key=lambda x:(x['t_min'],x['obs_id'])):
        if not groups or r['t_min']-groups[-1][0]['t_min']>span_days: groups.append([])
        groups[-1].append(r)
    return groups


def rank_pairs(groups, limit=100):
    if not isinstance(limit,int) or not 1<=limit<=1000: raise ValueError('invalid pair limit')
    candidates=[]
    for i,ga in enumerate(groups):
        for gb in groups[i+1:]:
            for a in ga:
                for b in gb:
                    if a['energy_bandpassname']!=b['energy_bandpassname']: continue
                    if any(x.get(k) is None for x in (a,b) for k in ('em_min','em_max')): continue
                    if min(a['em_max'],b['em_max'])<=max(a['em_min'],b['em_min']): continue
                    candidates.append(dict(a=a,b=b,separation_days=b['t_min']-a['t_min'],status='PRELIMINARY',limitations='Detector-wide wavelength filter only. Calibrated per-pixel responses, overlap, PSF, flags and background must be validated.'))
    candidates.sort(key=lambda p:(-p['separation_days'],p['a']['obs_id'],p['b']['obs_id']))
    return candidates[:limit]


def discover_pairs(ra=None,dec=None,target_id=None,radius=.05,epoch_days=1,limit=100):
    if target_id is not None:
        target=next((t for t in targets() if t['id']==target_id),None)
        if target is None: raise ValueError('unknown target')
        # Survey midpoint for initial search. Returned positions propagated separately.
        position=barycentric_position(target_id,'2025-11-01')
        ra,dec=position['ra_deg'],position['dec_deg']
    result=search(ra,dec,radius,'spherex_dr1')
    accepted=[]; rejected=[]
    for r in result['observations']:
        position = barycentric_position(target_id,Time(r['t_min'],format='mjd')) if target_id and r.get('t_min') is not None else dict(ra_deg=ra,dec_deg=dec)
        try:
            if not contains(r.get('s_region'),position['ra_deg'],position['dec_deg']):
                rejected.append(dict(obs_id=r['obs_id'],reason='target outside footprint')); continue
        except ValueError as exc:
            rejected.append(dict(obs_id=r['obs_id'],reason=str(exc))); continue
        accepted.append({**r,'target_position':position})
    groups=group_epochs(accepted,epoch_days)
    return dict(query=result['query'],observations=accepted,epochs=groups,pairs=rank_pairs(groups,limit),rejected=rejected,astrometry='Fixed ICRS or approximate barycentric proper motion; no observer parallax or radial velocity. Footprint matches are not detections.')
