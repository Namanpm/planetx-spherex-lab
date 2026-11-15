"""Native-grid point-source repeatability measurements, conditional on measured ePSF.

Jointly fits a sampled ePSF and local plane without resampling science pixels.
The WCS solid-angle flux conversion is an approximation, not the official pixel
solid-angle calibration map. Absolute photometric accuracy is not established.
"""
import numpy as np
from astropy.coordinates import SkyCoord
from pathlib import Path
from astropy.io import fits
from psf_processing import extract_epsf,prepare_epsf,extract_prepared_epsf
from fits_pipeline import read_cutout,read_cutout_hdus
from collections import OrderedDict
from threading import RLock
from astropy.wcs import NoConvergence


PIXEL_AREA_CACHE_SIZE=4096
_AREA_CACHE_LOCK=RLock()


def local_pixel_areas_sr(wcs,x,y):
    """Vector form of the existing spherical-offset Jacobian, in steradians.

    Coordinates retain zero-based native pixel centers and half-pixel offsets.
    Invalid geometries return NaN independently; no constant-area approximation.
    """
    x,y=np.broadcast_arrays(np.asarray(x,float),np.asarray(y,float))
    center=wcs.pixel_to_world(x,y)
    left=wcs.pixel_to_world(x-.5,y);right=wcs.pixel_to_world(x+.5,y)
    down=wcs.pixel_to_world(x,y-.5);up=wcs.pixel_to_world(x,y+.5)
    el,nl=center.spherical_offsets_to(left);er,nr=center.spherical_offsets_to(right)
    ed,nd=center.spherical_offsets_to(down);eu,nu=center.spherical_offsets_to(up)
    area=np.abs((er-el).rad*(nu-nd).rad-(eu-ed).rad*(nr-nl).rad)
    return np.where(np.isfinite(area)&(area>0),area,np.nan)


def _area_cache(observation,*,check_geometry=False):
    """Run-owned bounded cache. Prepared WCS objects are immutable during fitting.

    A new preparation batch checks the full WCS header (including SIP), clearing
    old values if callers deliberately replace or edit their WCS between batches.
    The lock also serializes cache creation and LRU access across fitting threads.
    """
    with _AREA_CACHE_LOCK:
        cache=observation.get('_pixel_area_cache')
        if cache is None or check_geometry:
            fingerprint=observation['wcs'].to_header(relax=True).tostring()
            if cache is None or cache['geometry']!=fingerprint:
                cache={'geometry':fingerprint,'values':OrderedDict(),'hits':0,'misses':0,'batch_points':0}
                observation['_pixel_area_cache']=cache
        return cache


def prepare_pixel_areas(observation,positions):
    """Populate position-dependent areas for RA/Dec pairs before source fitting.

    Scalar sky-to-pixel calls deliberately preserve the original fit coordinates.
    Only pixel-to-sky and spherical-offset work is vectorized. Invalid/outside
    catalog positions are skipped, leaving per-source validation unchanged.
    Returns the number of new valid areas retained in the bounded LRU.
    """
    with _AREA_CACHE_LOCK:
        cache=_area_cache(observation,check_geometry=True)
        height,width=observation['science']['image'].shape
        pending=OrderedDict()
        for position in positions:
            try:
                ra,dec=map(float,position)
                if not np.isfinite([ra,dec]).all() or not 0<=ra<360 or not -90<=dec<=90:continue
                x,y=observation['wcs'].world_to_pixel(SkyCoord(ra,dec,unit='deg'))
                x,y=float(x),float(y)
                if not np.isfinite([x,y]).all() or not (0<=x<width and 0<=y<height):continue
                key=(x,y)
                if key not in cache['values']:pending[key]=None
            except (ValueError,TypeError,NoConvergence):
                continue
        # Bound temporary coordinate arrays as well as retained geometry values.
        keys=list(pending)
        added=0
        for start in range(0,len(keys),PIXEL_AREA_CACHE_SIZE):
            group=keys[start:start+PIXEL_AREA_CACHE_SIZE]
            if not group:continue
            xy=np.asarray(group)
            try:
                areas=local_pixel_areas_sr(observation['wcs'],xy[:,0],xy[:,1])
            except (ValueError,TypeError):
                # Unusual WCS domains may fail as a batch. Scalar fitting remains
                # available and still raises the original source-specific error.
                continue
            cache['batch_points']+=len(group)
            for key,area in zip(group,areas):
                if not np.isfinite(area):continue
                cache['values'][key]=float(area);added+=1
                while len(cache['values'])>PIXEL_AREA_CACHE_SIZE:cache['values'].popitem(last=False)
        return added


def _source_pixel_area(observation,x,y):
    with _AREA_CACHE_LOCK:
        cache=_area_cache(observation)
        key=(float(x),float(y))
        if key in cache['values']:
            cache['hits']+=1;cache['values'].move_to_end(key)
            return cache['values'][key]
        cache['misses']+=1
        area=local_pixel_area_sr(observation['wcs'],*key)
        cache['values'][key]=area
        while len(cache['values'])>PIXEL_AREA_CACHE_SIZE:cache['values'].popitem(last=False)
        return area


def get_observation(path,observation_id,cache):
    """Run-owned dictionary keyed by observation/detector ID, guarded by cutout identity.

    A different path, size or modification timestamp replaces the entry. There is
    no global retention; each run owns and releases its prepared observations.
    """
    path=Path(path).resolve();stat=path.stat()
    fingerprint=(str(path),stat.st_size,stat.st_mtime_ns)
    existing=cache.get(observation_id)
    if existing is not None and existing['fingerprint']==fingerprint:return existing
    observation=prepare_observation(path)
    short=f"{observation['science']['header']['OBSID']}/D{observation['science']['header']['DETECTOR']}"
    if observation_id not in (short,'ivo://irsa.ipac/spherex_dr1?'+short,'ivo://irsa.ipac/spherex_qr3?'+short):
        raise ValueError('observation identity does not match FITS detector/OBSID')
    observation['fingerprint']=fingerprint
    cache[observation_id]=observation
    return observation


def prepare_observation(path):
    """Read science/WCS and all ePSF rows in one FITS open, then close the file."""
    with fits.open(path,memmap=False) as hdus:
        science=read_cutout_hdus(hdus,allow_source=True)
        epsf=prepare_epsf(hdus)
    for name in ('image','variance','flags','valid'):science[name].setflags(write=False)
    return {'science':science,'wcs':science['wcs'].celestial,'epsf':epsf}


def sample_epsf(raw,oversampling,x_offsets,y_offsets):
    """Bilinear measured oversampled ePSF at native pixel/source offsets."""
    raw=np.asarray(raw,float);ox,oy=oversampling
    x=np.asarray(x_offsets,float)*ox+(raw.shape[1]-1)/2
    y=np.asarray(y_offsets,float)*oy+(raw.shape[0]-1)/2
    inside=(x>=0)&(x<=raw.shape[1]-1)&(y>=0)&(y<=raw.shape[0]-1)
    ix=np.floor(x).astype(int);iy=np.floor(y).astype(int);fx=x-ix;fy=y-iy
    out=np.zeros(x.shape)
    for dx,dy,weight in ((0,0,(1-fx)*(1-fy)),(1,0,fx*(1-fy)),(0,1,(1-fx)*fy),(1,1,fx*fy)):
        out+=raw[np.clip(iy+dy,0,raw.shape[0]-1),np.clip(ix+dx,0,raw.shape[1]-1)]*weight
    return np.where(inside,out,np.nan),inside


def neighbor_profile(raw,oversampling,xx,yy,source_x,source_y):
    """Normalize over full native source support, retaining signed overlap wings."""
    ox,oy=oversampling;rx=(raw.shape[1]-1)/(2*ox);ry=(raw.shape[0]-1)/(2*oy)
    xmin,xmax=int(np.ceil(source_x-rx)),int(np.floor(source_x+rx))
    ymin,ymax=int(np.ceil(source_y-ry)),int(np.floor(source_y+ry))
    full_y,full_x=np.mgrid[ymin:ymax+1,xmin:xmax+1]
    complete,complete_valid=sample_epsf(raw,oversampling,full_x-source_x,full_y-source_y)
    norm=float(np.nansum(complete))
    if not np.isfinite(norm) or norm<=0:raise ValueError('full neighbor PSF has invalid normalization')
    overlap,inside=sample_epsf(raw,oversampling,xx-source_x,yy-source_y)
    return np.where(inside,overlap,0.)/norm


def fit_stamp(image,variance,valid,profile,x_offsets,y_offsets,neighbor_profiles=None):
    """SVD joint amplitude and plane fit with known independent native variance."""
    image,variance,profile=np.asarray(image,float),np.asarray(variance,float),np.asarray(profile,float)
    mask=np.asarray(valid,bool)&np.isfinite(image)&np.isfinite(variance)&(variance>0)&np.isfinite(profile)
    columns=[profile,np.ones_like(profile),x_offsets,y_offsets]
    if neighbor_profiles is not None:
        neighbor_profiles=np.asarray(neighbor_profiles,float)
        if neighbor_profiles.ndim!=3 or neighbor_profiles.shape[1:]!=profile.shape or len(neighbor_profiles)>12 or not np.isfinite(neighbor_profiles).all(): raise ValueError('invalid or excessive neighbor profiles')
        columns.extend(neighbor_profiles)
    design=np.stack(columns,axis=-1)
    if np.sum(mask)<6:raise ValueError('insufficient valid native fit pixels')
    weighted=design[mask]/np.sqrt(variance[mask,None]);values=image[mask]/np.sqrt(variance[mask])
    target=weighted[:,0];nuisance=weighted[:,1:]
    norms=np.linalg.norm(nuisance,axis=0);active=norms>0
    normalized=nuisance[:,active]/norms[active]
    u,s,vt=np.linalg.svd(normalized,full_matrices=False)
    tolerance=max(normalized.shape)*np.finfo(float).eps*(s[0] if s.size else 0.)
    rank=int(np.sum(s>tolerance));basis=u[:,:rank]
    perpendicular=target-basis@(basis.T@target)
    information=float(perpendicular@perpendicular);target_norm=float(np.linalg.norm(target))
    # Only the target amplitude must be identifiable. Nuisance null directions
    # are retained as an explicitly minimum-norm convention, never as constraints.
    if target_norm==0 or not np.isfinite(information) or np.sqrt(max(information,0.))<=np.sqrt(np.finfo(float).eps)*target_norm:
        raise ValueError('unidentifiable or ill-conditioned target amplitude')
    dof=int(mask.sum()-rank-1)
    if dof<2:raise ValueError('insufficient residual degrees of freedom')
    target_operator=perpendicular/information
    amplitude=float(target_operator@values)
    nuisance_pseudoinverse=(vt[:rank].T/s[:rank])@basis.T
    nuisance_operator=nuisance_pseudoinverse@(np.eye(len(values))-np.outer(target,target_operator))
    operator=np.zeros((weighted.shape[1],len(values)));operator[0]=target_operator
    operator[1:][active]=nuisance_operator/norms[active,None]
    beta=operator@values;cov=operator@operator.T
    residual=values-weighted@beta;chi2=float(residual@residual)
    condition=max(float(target_norm/np.sqrt(information)),float(s[0]/s[rank-1]) if rank else 1.)
    return {'amplitude':amplitude,'variance':float(1/information),'background_coefficients':beta[1:4].tolist(),'neighbor_amplitudes':beta[4:].tolist(),'parameter_covariance':cov.tolist(),'chi2':chi2,'dof':dof,'reduced_chi2':chi2/dof,'condition':condition,'valid_fraction':float(mask.mean()),'fit_pixels':int(mask.sum()),'nuisance_rank':rank,'nuisance_columns':int(nuisance.shape[1]),'dropped_zero_nuisance_columns':int(np.sum(~active)),'nuisance_coefficients':'minimum-norm convention if nuisance rank deficient; only target identifiability required','parameter_covariance_interpretation':'Linear-estimator covariance; nonidentified nuisance combinations are not measured constraints.'}


def local_pixel_area_sr(wcs,x,y):
    """Small-angle WCS Jacobian area around source, including supplied distortion."""
    center=wcs.pixel_to_world(x,y)
    left=wcs.pixel_to_world(x-.5,y);right=wcs.pixel_to_world(x+.5,y)
    down=wcs.pixel_to_world(x,y-.5);up=wcs.pixel_to_world(x,y+.5)
    el,nl=center.spherical_offsets_to(left);er,nr=center.spherical_offsets_to(right)
    ed,nd=center.spherical_offsets_to(down);eu,nu=center.spherical_offsets_to(up)
    area=abs((er-el).rad*(nu-nd).rad-(eu-ed).rad*(nr-nl).rad)
    if not np.isfinite(area) or area<=0:raise ValueError('invalid local WCS pixel solid angle')
    return float(area)


def fit_source(science_path,ra,dec,neighbor_positions=None,*,deblend=False):
    """Fit one source at fixed ICRS position. Neighbor list is RA/Dec degrees.

    Default rejects overlapping catalog neighbors. With deblend=True, fixed
    neighbor amplitudes are jointly fitted using the same local ePSF. No fit-chi2
    cutoff is silently imposed. Caller must
    screen fit quality and independent catalogs before repeatability inference.
    """
    observation=science_path if isinstance(science_path,dict) else prepare_observation(science_path)
    science=observation['science'];wcs=observation['wcs']
    if not np.isfinite([ra,dec]).all() or not 0<=ra<360 or not -90<=dec<=90:raise ValueError('invalid coordinates')
    x,y=wcs.world_to_pixel(SkyCoord(ra,dec,unit='deg'))
    if not np.isfinite([x,y]).all() or not (0<=x<science['image'].shape[1] and 0<=y<science['image'].shape[0]):raise ValueError('source outside cutout')
    psf=extract_prepared_epsf(observation['epsf'],float(x),float(y));raw=psf['oversampled_epsf'];ox,oy=psf['provenance']['oversampling_xy']
    rx=(raw.shape[1]-1)/(2*ox);ry=(raw.shape[0]-1)/(2*oy)
    xlo,xhi=int(np.ceil(x-rx)),int(np.floor(x+rx));ylo,yhi=int(np.ceil(y-ry)),int(np.floor(y+ry))
    if xlo<0 or ylo<0 or xhi>=science['image'].shape[1] or yhi>=science['image'].shape[0]:raise ValueError('ePSF fit support clipped by image edge')
    nuisance_positions=[]
    if neighbor_positions is not None:
        neighbors=np.asarray(neighbor_positions,float)
        if neighbors.size:
            if neighbors.ndim!=2 or neighbors.shape[1]!=2 or not np.isfinite(neighbors).all():raise ValueError('neighbor positions must be finite RA/Dec pairs')
            nx,ny=wcs.world_to_pixel(SkyCoord(neighbors[:,0],neighbors[:,1],unit='deg'))
            overlapping=(np.abs(nx-x)<=2*rx)&(np.abs(ny-y)<=2*ry)
            if overlapping.any() and not deblend:raise ValueError('neighbor can contaminate PSF fit support')
            nuisance_positions=list(zip(nx[overlapping],ny[overlapping]))
            if len(nuisance_positions)>12:raise ValueError('more than 12 contaminating neighbors in native fit support')
    yy,xx=np.mgrid[ylo:yhi+1,xlo:xhi+1];dx=xx-x;dy=yy-y
    profile,inside=sample_epsf(raw,(ox,oy),dx,dy);normalization=float(np.nansum(profile))
    if normalization<=0:raise ValueError('invalid sampled PSF normalization')
    profile/=normalization
    sl=(slice(ylo,yhi+1),slice(xlo,xhi+1));flags=science['flags'][sl]
    nuisance=[]
    for nx,ny in nuisance_positions:
        nuisance.append(neighbor_profile(raw,(ox,oy),xx,yy,nx,ny))
    result=fit_stamp(science['image'][sl],science['variance'][sl],science['valid'][sl]&inside,profile,dx,dy,np.array(nuisance) if nuisance else None)
    area=_source_pixel_area(observation,float(x),float(y));factor=area*1e6
    central_x,central_y=int(np.rint(x)),int(np.rint(y))
    if not (0<=central_x<science['image'].shape[1] and 0<=central_y<science['image'].shape[0]):raise ValueError('source center outside native grid')
    result.update(fit_bounds=[xlo,xhi,ylo,yhi],central_flag=int(science['flags'][central_y,central_x]),central_valid=bool(science['valid'][central_y,central_x]),neighbor_count=len(nuisance),deblend=bool(deblend),ra_deg=float(ra),dec_deg=float(dec),x=float(x),y=float(y),flux_jy=result['amplitude']*factor,variance_jy2=result['variance']*factor**2,
                  pixel_solid_angle_sr=area,amplitude_unit='MJy/sr times normalized finite-stamp pixel sum',flux_conversion='approximate local WCS pixel solid angle; official solid-angle map not used',
                  obs_id=science['header'].get('OBSID'),date=science['header'].get('DATE-OBS'),detector=int(science['header']['DETECTOR']),
                  flag_values=[int(v) for v in np.unique(flags)],bad_pixel_count=int(np.count_nonzero(~science['valid'][sl])),psf_provenance=psf['provenance'],
                  spectral_calibration='not retrieved here; caller associates wavelength_maps using observation publisher ID',
                  limitations=['Native errors treated independent; PSF, astrometric and photometric systematic uncertainty omitted.','Fixed-position ePSF plus local plane; optional catalog neighbors use same local ePSF approximation and truncated support.','Finite-support phase normalization and WCS solid-angle conversion require empirical validation.'])
    return result

