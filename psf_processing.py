"""Measured ePSF extraction and stable cross-convolution, without deconvolution.

Sources: https://caltech-ipac.github.io/irsa-tutorials/spherex-psf/ and
https://caltech-ipac.github.io/irsa-tutorials/spherex-intro/ . ePSF is already
pixel-convolved: a native centered kernel samples it, rather than integrates it.
Real sky matching additionally needs PSF rotation/reprojection onto common axes.
"""
import hashlib
import copy
import numpy as np
from astropy.io import fits
from spectral_calibration import detector_indices


def prepare_epsf(hdus):
    """Own a detached table and geometry; no FITS handles survive preparation."""
    image=hdus['IMAGE'];shape=image.data.shape
    yy,xx=detector_indices(image.header,shape)
    ext=hdus['EPSF'];table=ext.data
    required={'XCENTER','YCENTER','EPSF','NSTAR','CWAVE','CBAND'}
    if table is None or not required<=set(table.names):raise ValueError('unsupported EPSF schema')
    if str(ext.header.get('DETECTOR','')).strip()!=str(image.header['DETECTOR']).strip():raise ValueError('EPSF detector mismatch')
    rows={name:np.array(table[name],dtype=float if name in ('XCENTER','YCENTER') else None,copy=True) for name in required}
    for name in ('XCENTER','YCENTER'):
        if not np.isfinite(rows[name]).all():raise ValueError('invalid EPSF zone coordinates')
    for value in rows.values():value.setflags(write=False)
    return {'shape':shape,'x0':int(xx[0]),'y0':int(yy[0]),'detector':int(image.header['DETECTOR']),
            'rows':rows,'oversampling':(int(ext.header['OVSMPX']),int(ext.header['OVSMPY'])),
            'history':tuple(ext.header.get('HISTORY',[])),'zones':{}}


def extract_prepared_epsf(prepared,x=None,y=None):
    """Reuse detached ePSF rows and memoize each requested detector zone."""
    shape=prepared['shape'];rows=prepared['rows']
    x=(shape[1]-1)/2 if x is None else float(x)
    y=(shape[0]-1)/2 if y is None else float(y)
    if not np.isfinite(x+y) or not (0<=x<=shape[1]-1 and 0<=y<=shape[0]-1):raise ValueError('PSF point outside cutout')
    xorig=prepared['x0']+x;yorig=prepared['y0']+y
    cx=rows['XCENTER'];cy=rows['YCENTER']
    row_index=int(np.argmin((cx-xorig)**2+(cy-yorig)**2))
    if row_index not in prepared['zones']:
        raw=np.array(rows['EPSF'][row_index],dtype=float)
        ox,oy=prepared['oversampling']
        if raw.ndim!=2 or any(n%2!=1 for n in raw.shape) or min(ox,oy)<1 or not np.isfinite(raw).all():raise ValueError('invalid oversampled EPSF')
        if int(rows['NSTAR'][row_index])<=0:raise ValueError('PSF has no calibration stars')
        iy=np.arange(-(raw.shape[0]//2//oy),raw.shape[0]//2//oy+1)*oy+raw.shape[0]//2
        ix=np.arange(-(raw.shape[1]//2//ox),raw.shape[1]//2//ox+1)*ox+raw.shape[1]//2
        kernel=raw[np.ix_(iy,ix)];total=kernel.sum()
        if not np.isfinite(total) or total<=0:raise ValueError('invalid EPSF normalization')
        kernel=kernel/total;raw.setflags(write=False);kernel.setflags(write=False)
        provenance={'detector':prepared['detector'],'zone_index':row_index,'zone_center_detector_xy':[float(cx[row_index]),float(cy[row_index])],
            'oversampling_xy':[ox,oy],'nstar':int(rows['NSTAR'][row_index]),'wave_um':float(rows['CWAVE'][row_index]),'band_um':float(rows['CBAND'][row_index]),
            'raw_sha256':hashlib.sha256(raw.tobytes()).hexdigest(),'sampling':'pixel-centered samples of measured pixel-convolved ePSF',
            'limitations':['Single detector zone; not a spatially varying convolution.','Detector axes must be transformed to common sky axes before real cross-matching.','Sampled finite support does not model unmeasured wings or PSF uncertainty.']}
        prepared['zones'][row_index]={'kernel':kernel,'oversampled_epsf':raw,'provenance':provenance}
    zone=prepared['zones'][row_index]
    provenance={**copy.deepcopy(zone['provenance']),'requested_detector_xy':[xorig,yorig],'source_history':list(prepared['history'])}
    return {'kernel':zone['kernel'],'oversampled_epsf':zone['oversampled_epsf'],'provenance':provenance}


def extract_epsf(path,x=None,y=None):
    """Path-based compatibility entry point; batch callers prepare once instead."""
    with fits.open(path,memmap=False) as hdus:
        prepared=prepare_epsf(hdus)
    return extract_prepared_epsf(prepared,x,y)


def convolve_independent(image,variance,valid,kernel):
    """Linear convolution on independent native pixels; variance uses squared weights.

    Every nonzero kernel-support pixel must be valid. Borders are invalid rather
    than zero-extended. Output diagonal variance omits induced neighboring-pixel
    covariance; applying this to resampled correlated pixels requires full covariance.
    Signed empirical PSF lobes are allowed, but positive total normalization is required.
    """
    a=np.asarray(image,float);v=np.asarray(variance,float);mask=np.asarray(valid,bool);k=np.asarray(kernel,float)
    if a.ndim!=2 or v.shape!=a.shape or mask.shape!=a.shape: raise ValueError('invalid image arrays')
    if k.ndim!=2 or any(n%2!=1 for n in k.shape) or not np.isfinite(k).all() or k.sum()<=0: raise ValueError('invalid kernel')
    k=k/k.sum(); mask=mask & np.isfinite(a) & np.isfinite(v) & (v>=0)
    ry,rx=np.array(k.shape)//2
    aa=np.pad(np.where(mask,a,0),((ry,ry),(rx,rx)));vv=np.pad(np.where(mask,v,0),((ry,ry),(rx,rx)));mm=np.pad(mask,((ry,ry),(rx,rx)))
    out=np.zeros_like(a);var=np.zeros_like(a);good=np.ones_like(mask)
    for iy,ix in np.argwhere(k!=0):
        weight=k[iy,ix]; y0=2*ry-iy;x0=2*rx-ix
        sl=(slice(y0,y0+a.shape[0]),slice(x0,x0+a.shape[1]))
        out+=weight*aa[sl];var+=weight**2*vv[sl];good &= mm[sl]
    out[~good]=np.nan;var[~good]=np.nan
    return {'image':out,'variance':var,'valid':good,'kernel':k,'covariance':'diagonal only; convolution creates output correlations'}


def cross_match(a,va,b,vb,valid_a,valid_b,psf_a,psf_b):
    """Convolve A by measured B PSF and B by measured A PSF on common axes.

    Effective PSF is A*B; no unstable division/deconvolution is used. Supplied
    kernels must already have matched pixel scale, centering and sky orientation.
    """
    ra=convolve_independent(a,va,valid_a,psf_b)
    rb=convolve_independent(b,vb,valid_b,psf_a)
    return {'a':ra['image'],'b':rb['image'],'variance_a':ra['variance'],'variance_b':rb['variance'],'valid':ra['valid'] & rb['valid'],'limitations':['PSFs treated as known; kernel uncertainty omitted.','Convolution generates neighboring-pixel covariance.']}


def fit_native_background(image,variance,blank_mask,*,independent_native=True):
    """Weighted plane on caller-identified blank native pixels, with fit covariance.

    This is not automatic source masking, clipping or a robust-sky estimator. The
    supplied mask must exclude SOURCE and artifacts. Off-diagonal residual
    covariance exists even for native pixels; return plane coefficient covariance
    so downstream model fits can account for it. No resampled input is permitted.
    """
    if not independent_native: raise ValueError('background requires independent native pixels or a full covariance model')
    a=np.asarray(image,float);v=np.asarray(variance,float);blank=np.asarray(blank_mask,bool)
    if a.ndim!=2 or v.shape!=a.shape or blank.shape!=a.shape: raise ValueError('invalid background arrays')
    finite=np.isfinite(a)&np.isfinite(v)&(v>0)
    blank=blank & finite
    y,x=np.indices(a.shape,dtype=float)
    # Normalize coordinates to keep the native-grid plane well conditioned.
    x=(x-(a.shape[1]-1)/2)/max(a.shape[1],1);y=(y-(a.shape[0]-1)/2)/max(a.shape[0],1)
    design=np.stack([np.ones_like(a),x,y],axis=-1);z=design[blank]
    if len(z)<4: raise ValueError('insufficient independent blank pixels')
    weighted=z/np.sqrt(v[blank,None])
    if np.linalg.matrix_rank(weighted)!=3: raise ValueError('blank sky does not constrain a plane')
    covariance=np.linalg.inv(weighted.T@weighted)
    coefficients=covariance@(z.T@(a[blank]/v[blank]))
    background=design@coefficients
    prediction_variance=np.einsum('...i,ij,...j->...',design,covariance,design)
    residual_variance=v+prediction_variance-2*prediction_variance*blank
    residual=np.where(finite,a-background,np.nan)
    residual_variance=np.where(finite,np.maximum(residual_variance,0),np.nan)
    return {'background':background,'residual':residual,'residual_variance':residual_variance,'coefficients':coefficients,'coefficient_covariance':covariance,'blank_mask':blank,'design':design,'limitations':['Marginal variances only; residuals share fitted background covariance.','Mask supplied by caller; no automatic artifact rejection or zodiacal-model correction.']}


def project_kernel(kernel,source_wcs,target_wcs,sky_coord,shape_out=None):
    """Transfer a finite native PSF to common sky axes by geometric sampling.

    Evaluates source detector kernel offsets corresponding to each target pixel
    offset at sky_coord. This is an approximate interpolation of a finite native
    sampled ePSF, not a reconstruction from oversampled data. Kernel uncertainty,
    detector-zone variation and intra-pixel phase are not propagated. Signed
    empirical lobes remain signed. No extrapolation or deconvolution is used.
    """
    k=np.asarray(kernel,dtype=float)
    if k.ndim!=2 or any(n%2!=1 for n in k.shape) or not np.isfinite(k).all() or k.sum()<=0:
        raise ValueError('invalid source kernel')
    if not source_wcs.has_celestial or not target_wcs.has_celestial or source_wcs.pixel_n_dim!=2 or target_wcs.pixel_n_dim!=2:
        raise ValueError('two-dimensional celestial WCS required')
    sx,sy=source_wcs.world_to_pixel(sky_coord);tx,ty=target_wcs.world_to_pixel(sky_coord)
    if not np.isfinite([sx,sy,tx,ty]).all(): raise ValueError('sky position not representable')
    cy,cx=np.array(k.shape)//2
    if shape_out is None:
        # Sample all finite support perimeter points, not only corners, because
        # celestial projection/SIP mappings need not be affine across support.
        yy,xx=np.indices(k.shape,dtype=float)
        perimeter=(yy==0)|(xx==0)|(yy==k.shape[0]-1)|(xx==k.shape[1]-1)
        world=source_wcs.pixel_to_world(sx+xx[perimeter]-cx,sy+yy[perimeter]-cy)
        px,py=target_wcs.world_to_pixel(world)
        if not np.isfinite(px).all() or not np.isfinite(py).all(): raise ValueError('kernel projection undefined')
        extent_x=int(np.ceil(np.max(np.abs(px-tx))));extent_y=int(np.ceil(np.max(np.abs(py-ty))))
        shape_out=(2*extent_y+1,2*extent_x+1)
    if len(shape_out)!=2 or any(not isinstance(n,(int,np.integer)) or n<=0 or n%2!=1 or n>101 for n in shape_out):
        raise ValueError('output kernel must have odd positive dimensions <=101')
    yy,xx=np.indices(shape_out,dtype=float)
    world=target_wcs.pixel_to_world(tx+xx-shape_out[1]//2,ty+yy-shape_out[0]//2)
    px,py=source_wcs.world_to_pixel(world)
    px=px-sx+cx;py=py-sy+cy
    finite=np.isfinite(px)&np.isfinite(py)
    px=np.where(finite,px,-100.);py=np.where(finite,py,-100.)
    px=np.where(np.abs(px-np.rint(px))<1e-7,np.rint(px),px)
    py=np.where(np.abs(py-np.rint(py))<1e-7,np.rint(py),py)
    ix=np.floor(px).astype(int);iy=np.floor(py).astype(int);fx=px-ix;fy=py-iy
    out=np.zeros(shape_out,dtype=float)
    for dx,dy,weight in ((0,0,(1-fx)*(1-fy)),(1,0,fx*(1-fy)),(0,1,(1-fx)*fy),(1,1,fx*fy)):
        x=ix+dx;y=iy+dy
        valid=finite&(x>=0)&(x<k.shape[1])&(y>=0)&(y<k.shape[0])
        out+=np.where(valid,k[np.clip(y,0,k.shape[0]-1),np.clip(x,0,k.shape[1]-1)],0)*weight
    normalization=out.sum()
    if not np.isfinite(normalization) or normalization<=0: raise ValueError('projected kernel has no positive normalization')
    return {'kernel':out/normalization,'normalization_before_scaling':float(normalization),'source_shape':list(k.shape),'target_shape':list(shape_out),'limitations':['Finite native PSF interpolation; oversampled subpixel response is not reconstructed.','Kernel uncertainty and detector-zone variation are omitted.','Output normalization does not validate PSF accuracy or response compatibility.']}
