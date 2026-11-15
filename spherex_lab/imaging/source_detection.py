"""Positive connected excursions and provenance-preserving positional candidates.

This is not PSF photometry or a validated astronomical source catalog. Centroid
uncertainties describe a conditional fixed-pixel flux statistic, not astrometric
precision; threshold selection, blending, PSF and intra-pixel uncertainty remain.
"""
import numpy as np


def detect_sources(image, variance, valid, *, sigma_threshold, min_pixels=1, correlated=False, max_components=10000):
    """Find 8-connected positive excursions on a background-subtracted image.

    sigma_threshold is an explicit analysis policy. Input variance is marginal
    per-pixel measurement variance. Correlated inputs can still be segmented,
    but independent-pixel summed-flux and centroid covariance are withheld.
    Background-estimation error must already be represented in the input model.
    """
    a=np.asarray(image,dtype=float);v=np.asarray(variance,dtype=float);mask=np.asarray(valid,dtype=bool)
    if a.ndim!=2 or v.shape!=a.shape or mask.shape!=a.shape or a.size>1000000:
        raise ValueError('matching 2D arrays with at most one million pixels required')
    if not np.isfinite(sigma_threshold) or sigma_threshold<=0 or not isinstance(min_pixels,int) or min_pixels<1:
        raise ValueError('positive explicit threshold and component size required')
    if not isinstance(max_components,int) or max_components<1: raise ValueError('positive component resource limit required')
    good=mask & np.isfinite(a) & np.isfinite(v) & (v>0)
    sigma=np.divide(a,np.sqrt(np.where(good,v,1)),out=np.zeros_like(a),where=good)
    eligible=good & (sigma>=sigma_threshold)
    seen=np.zeros_like(mask);sources=[];component_count=0
    height,width=a.shape
    for root_y,root_x in np.argwhere(eligible):
        if seen[root_y,root_x]: continue
        component_count+=1
        if component_count>max_components: raise ValueError('connected component resource limit exceeded')
        stack=[(int(root_y),int(root_x))];seen[root_y,root_x]=True;component=[]
        while stack:
            y,x=stack.pop();component.append((y,x))
            for dy in (-1,0,1):
                for dx in (-1,0,1):
                    yy,xx=y+dy,x+dx
                    if 0<=yy<height and 0<=xx<width and eligible[yy,xx] and not seen[yy,xx]:
                        seen[yy,xx]=True;stack.append((yy,xx))
        if len(component)<min_pixels: continue
        yy,xx=np.asarray(component).T;flux=a[yy,xx];vv=v[yy,xx];total=float(np.sum(flux))
        cx=float(np.sum(xx*flux)/total);cy=float(np.sum(yy*flux)/total)
        covariance=None
        if not correlated and len(component)>1:
            # Delta method for fixed component membership and pixel locations.
            jacobian=np.column_stack([xx-cx,yy-cy])/total
            covariance=jacobian.T@(vv[:,None]*jacobian)
        sources.append({'x':cx,'y':cy,'pixel_count':len(component),'flux_sum':total,
                        'flux_variance':None if correlated else float(np.sum(vv)),
                        'peak_marginal_sigma':float(np.max(sigma[yy,xx])),
                        'centroid_covariance':None if covariance is None else covariance.tolist(),
                        'astrometric_uncertainty_available':False,
                        'classification':'unclassified positive excursion'})
    sources.sort(key=lambda s:s['peak_marginal_sigma'],reverse=True)
    return {'sources':sources,'threshold_sigma':float(sigma_threshold),'min_pixels':min_pixels,
            'correlated_input':bool(correlated),'valid_fraction':float(good.mean()) if a.size else 0.,
            'limitations':['Excursions are not confirmed sources; segmentation may blend objects.',
            'Centroid covariance conditions on fixed threshold-selected component membership.',
            'Centroid covariance excludes PSF, intra-pixel position and astrometric systematics.',
            'No independent-input uncertainty interpretation is given for correlated inputs.']}


def crossmatch_catalog(source_ra_dec_deg, catalog_ra_dec_deg, *, max_separation_arcsec, catalog_provenance):
    """Return unambiguous spherical positional candidates, never verified identity.

    Both catalogs use ICRS degrees at a caller-established compatible epoch.
    Motion propagation and uncertainty-aware identity modeling are not implied.
    """
    source=np.asarray(source_ra_dec_deg,dtype=float);catalog=np.asarray(catalog_ra_dec_deg,dtype=float)
    for a in (source,catalog):
        if a.ndim!=2 or a.shape[1]!=2 or not np.isfinite(a).all() or np.any((a[:,0]<0)|(a[:,0]>=360)) or np.any(np.abs(a[:,1])>90):
            raise ValueError('finite ICRS RA/Dec degree arrays required')
    if source.shape[0]*catalog.shape[0]>1000000: raise ValueError('catalog match exceeds pair resource limit')
    if not catalog_provenance or not np.isfinite(max_separation_arcsec) or not 0<max_separation_arcsec<=648000:
        raise ValueError('catalog provenance and positive explicit search radius required')
    a=np.deg2rad(source);b=np.deg2rad(catalog)
    delta_ra=a[:,None,0]-b[None,:,0];delta_dec=a[:,None,1]-b[None,:,1]
    h=np.sin(delta_dec/2)**2+np.cos(a[:,None,1])*np.cos(b[None,:,1])*np.sin(delta_ra/2)**2
    h=np.clip(h,0,1);separation=np.rad2deg(2*np.arctan2(np.sqrt(h),np.sqrt(1-h)))*3600
    eligible=separation<=max_separation_arcsec
    ns=eligible.sum(axis=1);nc=eligible.sum(axis=0)
    accepted=eligible & (ns[:,None]==1) & (nc[None,:]==1)
    matches=[{'source_index':int(i),'catalog_index':int(j),'separation_arcsec':float(separation[i,j]),
              'association':'positional candidate, not confirmed identity'} for i,j in np.argwhere(accepted)]
    return {'matches':matches,'ambiguous_source_indices':np.flatnonzero(ns>1).tolist(),
            'catalog_provenance':catalog_provenance,'max_separation_arcsec':float(max_separation_arcsec),
            'limitations':['Epoch compatibility and catalog astrometry require caller verification.',
            'Crowding, chance coincidences and catalog identity are not modeled.']}
