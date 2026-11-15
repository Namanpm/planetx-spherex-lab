"""Common-sky fiducial response sampling and signed laboratory shape distances.

A shape distance ranks candidates; it never proves response equivalence or
supplies an unmeasured error budget. Laboratory negative samples are preserved.
"""
import numpy as np
from astropy.coordinates import SkyCoord
from spherex_lab.imaging.fits_pipeline import read_cutout
from spherex_lab.imaging.response_calibration import pixel_response


def field_response(path,wave_map,bundle,ra,dec):
    """Sample closest science/calibration-valid native pixel to shared sky point.

    Requested point must lie in the cutout footprint. The actual sampled pixel,
    sky coordinate and separation are returned; no unreported sky substitution.
    """
    if not np.isfinite([ra,dec]).all() or not 0<=ra<360 or not -90<=dec<=90:raise ValueError('invalid response sky coordinate')
    science=read_cutout(path,allow_source=True);wave=np.asarray(wave_map,float)
    if wave.shape!=science['image'].shape:raise ValueError('wavelength grid does not match science grid')
    wcs=science['wcs'].celestial
    px,py=wcs.world_to_pixel_values(float(ra),float(dec))
    height,width=wave.shape
    if not np.isfinite([px,py]).all() or not -.5<=px<width-.5 or not -.5<=py<height-.5:raise ValueError('requested sky coordinate outside cutout')
    good=science['valid'] & np.isfinite(wave) & (wave>0)
    y,x=np.nonzero(good)
    if not len(x):raise ValueError('no valid wavelength/science sample')
    nearest=int(np.argmin((x-px)**2+(y-py)**2));ix,iy=int(x[nearest]),int(y[nearest])
    wavelengths,response,provenance=pixel_response(bundle,science['header'],wave.shape,ix,iy,float(wave[iy,ix]))
    actual_ra,actual_dec=wcs.pixel_to_world_values(ix,iy)
    offset=SkyCoord(float(ra),float(dec),unit='deg').separation(SkyCoord(float(actual_ra),float(actual_dec),unit='deg')).arcsec
    return wavelengths,response,{**provenance,'requested_sky_deg':[float(ra),float(dec)],'sampled_pixel_xy':[ix,iy],'sampled_sky_deg':[float(actual_ra),float(actual_dec)],'sampling_offset_arcsec':float(offset),'sampling':'nearest native science/calibration-valid pixel; not interpolated sky response','compatibility_status':'UNVERIFIED'}


def _samples(curve):
    w,r=np.asarray(curve[0],float),np.asarray(curve[1],float)
    if w.ndim!=1 or r.shape!=w.shape or len(w)<2 or not np.isfinite(w).all() or not np.isfinite(r).all() or not np.all(np.diff(w)>0):raise ValueError('invalid signed laboratory response')
    return w,r


def signed_laboratory_distance(curve_a,curve_b):
    """L2-normalized signed sample shape difference; no clipping or threshold.

    Curves are piecewise linear within their measured support and zero outside.
    Squared/product integrals are evaluated exactly on each linear segment.
    Distances ignore overall amplitude, preserve wavelength location and negative
    lab noise, and must not be interpreted as physical throughput overlap.
    """
    wa,a=_samples(curve_a);wb,b=_samples(curve_b)
    grid=np.unique(np.r_[wa,wb]);lo,hi=grid[:-1],grid[1:];mid=(lo+hi)/2;width=hi-lo
    def edges(w,r):
        within=(mid>=w[0])&(mid<=w[-1])
        return np.where(within,np.interp(lo,w,r),0),np.where(within,np.interp(hi,w,r),0)
    a0,a1=edges(wa,a);b0,b1=edges(wb,b)
    norm_a=float(np.sum(width*(a0*a0+a0*a1+a1*a1)/3));norm_b=float(np.sum(width*(b0*b0+b0*b1+b1*b1)/3))
    if norm_a<=0 or norm_b<=0:raise ValueError('zero response L2 norm')
    product=float(np.sum(width*(2*a0*b0+a0*b1+a1*b0+2*a1*b1)/6))
    cosine=float(np.clip(product/np.sqrt(norm_a*norm_b),-1,1))
    distance=float(np.sqrt(max(0,2-2*cosine)))
    return {'normalized_l2_distance':distance,'signed_normalized_inner_product':cosine,'negative_sample_counts':[int((a<0).sum()),int((b<0).sum())],'semantics':'Signed piecewise-linear laboratory shape distance; amplitude normalized, wavelength location retained. Not physical throughput overlap.','acceptance':'NOT_ASSESSED','systematic_uncertainty':'not estimated by this metric'}
