"""Official SPHEREx fiducial responses, explicitly distinct from exact pixel responses.

Official definition: https://irsa.ipac.caltech.edu/data/SPHEREx/docs/SPHERExSpectralFiles_202607.pdf
The channel middle-pixel curve may be shifted in wavelength for off-center
pixels (§3). This is an approximation with no delivered per-pixel error bound.
It must not, on its own, authorize a calibrated direct subtraction.
"""
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import tarfile
import tempfile
import numpy as np
from spherex_lab.data.calibration_resolver import session

DOCUMENTATION='https://irsa.ipac.caltech.edu/data/SPHEREx/docs/SPHERExSpectralFiles_202607.pdf'
BUNDLE_URL='https://irsa.ipac.caltech.edu/data/SPHEREx/docs/spectral_files_20260812/SPHERExSpectralFiles_202607.tar'
MAX_BUNDLE_BYTES=450*1024*1024


def _curve(wavelength, response, allow_negative=False):
    w,r=np.asarray(wavelength,dtype=float),np.asarray(response,dtype=float)
    if w.ndim!=1 or r.shape!=w.shape or len(w)<2 or not np.isfinite(w).all() or not np.isfinite(r).all() or not np.all(np.diff(w)>0) or (np.any(r<0) and not allow_negative) or np.trapezoid(r,w)<=0:
        raise ValueError('invalid wavelength/response samples')
    return w,r


def parse_response_csv(content):
    """Read documented wavelength [um], normalized transmission CSV samples."""
    text=content.decode('utf-8-sig') if isinstance(content,bytes) else content
    rows=list(csv.reader(io.StringIO(text)))
    rows=[r for r in rows if r and any(s.strip() for s in r)]
    if len(rows)<3 or len(rows[0])!=2:raise ValueError('invalid response CSV')
    header=[''.join(c.lower() for c in s if c.isalnum()) for s in rows[0]]
    if header[0] not in ('wavelengthum','wavelength') or header[1] not in ('normtransmission','normalizedtransmission'):
        raise ValueError('unknown response CSV fields')
    try: data=np.array([[float(v) for v in row] for row in rows[1:]])
    except (ValueError,TypeError) as exc:raise ValueError('invalid response CSV data') from exc
    return _curve(data[:,0],data[:,1],allow_negative=True)


def response_metrics(wavelength_a,response_a,wavelength_b,response_b):
    """Area-normalized overlap and centroid diagnostics; no acceptance threshold.

    Piecewise-linear interpolation onto a common grid with zero outside each
    measured response support. Does not quantify absent calibration systematics.
    """
    wa,a=_curve(wavelength_a,response_a);wb,b=_curve(wavelength_b,response_b)
    grid=np.unique(np.r_[wa,wb])
    left,right=grid[:-1],grid[1:];mid=(left+right)/2
    area_a=np.trapezoid(a,wa);area_b=np.trapezoid(b,wb)
    a0=np.interp(left,wa,a)/area_a;a1=np.interp(right,wa,a)/area_a
    b0=np.interp(left,wb,b)/area_b;b1=np.interp(right,wb,b)/area_b
    inside_a=(mid>=wa[0])&(mid<=wa[-1]);inside_b=(mid>=wb[0])&(mid<=wb[-1])
    a0=np.where(inside_a,a0,0);a1=np.where(inside_a,a1,0)
    b0=np.where(inside_b,b0,0);b1=np.where(inside_b,b1,0)
    d0=a0-b0;d1=a1-b1
    crossing=d0*d1<0
    fraction=np.divide(d0,d0-d1,out=np.zeros_like(d0),where=crossing)
    shared=a0+fraction*(a1-a0)
    noncross=(np.minimum(a0,b0)+np.minimum(a1,b1))/2
    cross=fraction*(np.minimum(a0,b0)+shared)/2+(1-fraction)*(shared+np.minimum(a1,b1))/2
    overlap=np.sum((right-left)*np.where(crossing,cross,noncross))
    def centroid(w,r,area):
        return np.sum(np.diff(w)*((2*w[:-1]+w[1:])*r[:-1]+(w[:-1]+2*w[1:])*r[1:])/6)/area
    return {'normalized_overlap':float(overlap),
            'centroid_difference_um':float(centroid(wb,b,area_b)-centroid(wa,a,area_a)),
            'exact_samples_equal':bool(np.array_equal(wa,wb) and np.array_equal(a,b)),
            'acceptance':'not assessed; overlap is a diagnostic, not an equivalence proof'}



def shifted_fiducial(wavelength,response,*,target_center_um,fiducial_center_um):
    """Apply documented approximate wavelength translation only; no invented stretch."""
    w,r=_curve(wavelength,response,allow_negative=True)
    if not np.isfinite([target_center_um,fiducial_center_um]).all() or min(target_center_um,fiducial_center_um)<=0:raise ValueError('invalid centers')
    return w+target_center_um-fiducial_center_um,r.copy(),{
        'status':'APPROXIMATE','method':'middle-pixel fiducial wavelength translation',
        'documentation':DOCUMENTATION,'response_systematic_uncertainty':'unquantified',
        'observation_specific_response':'not established'}


def fetch_bundle(cache_dir,client=None):
    """Bounded streamed atomic cache. This release is approximately 437 MB."""
    cache=Path(cache_dir);cache.mkdir(parents=True,exist_ok=True);path=cache/'response_bundle.tar'
    if path.exists():
        with tarfile.open(path) as tar:tar.getmembers()
        return path
    own=client is None;client=client or session();temporary=None
    try:
        with client.get(BUNDLE_URL,stream=True,timeout=(10,90),allow_redirects=False) as r:
            r.raise_for_status()
            if r.status_code!=200:raise ValueError('unexpected response bundle status')
            if int(r.headers.get('Content-Length','0'))>MAX_BUNDLE_BYTES:raise ValueError('response bundle too large')
            with tempfile.NamedTemporaryFile(dir=cache,suffix='.part',delete=False) as f:
                temporary=Path(f.name);size=0
                for chunk in r.iter_content(1024*1024):
                    size+=len(chunk)
                    if size>MAX_BUNDLE_BYTES:raise ValueError('response bundle too large')
                    f.write(chunk)
        with tarfile.open(temporary) as tar:tar.getmembers()
        os.replace(temporary,path);temporary=None
        return path
    finally:
        if temporary is not None:temporary.unlink(missing_ok=True)
        if own:client.close()


def load_fiducial(bundle_path,detector,channel,subchannels=17):
    """Load exact named MiddlePixel curve, without extracting arbitrary tar paths.

    Channel numbering must come from this release's geometric channel map.
    Neither CWAVE nor CBAND identifies a geometric channel by itself.
    """
    if detector not in range(1,7) or subchannels not in (17,34,51,68) or not 0<=channel<=subchannels+1:raise ValueError('invalid detector/channel')
    suffix=f'MiddlePixel/SubChannel{subchannels}/Band{detector}/FiducialSpec_Band{detector}_SubCh{channel:02d}_{subchannels}_Cent.csv'
    with tarfile.open(bundle_path) as tar:
        matches=[m for m in tar.getmembers() if m.name.endswith(suffix) and m.isfile()]
        if len(matches)!=1:raise ValueError('missing or ambiguous fiducial response')
        m=matches[0]
        if m.size>4*1024*1024:raise ValueError('response CSV too large')
        content=tar.extractfile(m).read()
    w,r=parse_response_csv(content)
    return w,r,{'source_url':BUNDLE_URL,'documentation':DOCUMENTATION,'member':m.name,
        'member_sha256':hashlib.sha256(content).hexdigest(),'detector':detector,'channel':channel,
        'response_kind':'representative channel middle-pixel average','exact_pixel_response':False,'negative_samples':int(np.count_nonzero(r<0)), 'negative_policy':'preserved laboratory values; no clipping or interpretation as physical negative throughput'}

def pixel_response(bundle_path, science_header, shape, x, y, target_center_um, subchannels=17):
    """Documented shifted middle-pixel approximation at one cutout pixel.

    Detector-map orientation is reversed along X for MWIR detectors 4–6 as
    specified by IRSA. Delivered summary center, not an invented centroid, sets
    the wavelength shift. Edge channels lacking a summary are rejected.
    """
    from astropy.io import fits
    from spherex_lab.imaging.spectral_calibration import detector_indices
    detector=int(science_header['DETECTOR'])
    if detector not in range(1,7) or not isinstance(x,(int,np.integer)) or not isinstance(y,(int,np.integer)) or not 0<=x<shape[1] or not 0<=y<shape[0]:raise ValueError('invalid detector or cutout pixel')
    iy,ix=detector_indices(science_header,shape)
    original_x=int(ix[x]); original_y=int(iy[y])
    lab_x=2039-original_x if detector>=4 else original_x
    with tarfile.open(bundle_path) as tar:
        map_name=f'ChannelMaps/spectral_channels_{subchannels}ch_2040.fits'
        m=tar.getmember(map_name)
        if not m.isfile() or m.size>10*1024*1024:raise ValueError('invalid channel map member')
        with fits.open(io.BytesIO(tar.extractfile(m).read())) as f:
            data=f[1].data
            if data.shape!=(2040,2040) or data.dtype.kind not in 'iu':raise ValueError('invalid channel map geometry')
            channel=int(data[original_y,lab_x])
        summary=f'FidSpecSummary/SubChannel{subchannels}/Band{detector}/FidSpecSummary_Band{detector}_SpecCh{subchannels}.csv'
        m=tar.getmember(summary)
        if not m.isfile() or m.size>1024*1024:raise ValueError('invalid summary member')
        rows=list(csv.DictReader(io.StringIO(tar.extractfile(m).read().decode('utf-8-sig'))))
    matches=[r for r in rows if int(r['SpecCh'])==channel]
    if len(matches)!=1:raise ValueError('no fiducial summary for mapped channel')
    center=float(matches[0]['BandCenter_um'])
    wave,response,provenance=load_fiducial(bundle_path,detector,channel,subchannels)
    wave,response,approximation=shifted_fiducial(wave,response,target_center_um=target_center_um,fiducial_center_um=center)
    provenance.update(approximation)
    provenance.update(detector_x=original_x,detector_y=original_y,laboratory_x=lab_x,fiducial_center_um=center,target_center_um=float(target_center_um))
    return wave,response,provenance
