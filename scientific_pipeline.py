"""Bounded automatic DR1 pair investigation; rejection is not a detection result."""
import argparse
import hashlib
import json
import os
import platform
import tempfile
import time
from pathlib import Path
from urllib.parse import urlparse, urlencode
import requests
import numpy as np
from astropy.io import fits
from calibration_resolver import session, wavelength_maps
from fits_pipeline import read_cutout


def download_cutout(record, ra, dec, cache, size=0.05):
    if not np.isfinite(size) or not 0<size<=0.3: raise ValueError('cutout size must be (0,0.3] degrees')
    url = record['access_url']; parsed = urlparse(url)
    if parsed.scheme != 'https' or parsed.netloc != 'irsa.ipac.caltech.edu' or not parsed.path.startswith('/ibe/data/spherex/dr1/level2/') or parsed.query:
        raise ValueError('unsupported science product URL')
    url += '?' + urlencode({'center':f'{ra},{dec}', 'size':size})
    cache=Path(cache); cache.mkdir(parents=True,exist_ok=True)
    path=cache/(hashlib.sha256(url.encode()).hexdigest()+'.fits')
    if path.exists(): read_cutout(path); return path
    temp=None
    try:
        with session() as client, client.get(url,stream=True,timeout=(10,90),allow_redirects=False) as response:
            response.raise_for_status()
            if response.status_code != 200: raise ValueError('unexpected archive HTTP status')
            with tempfile.NamedTemporaryFile(dir=cache,suffix='.part',delete=False) as stream:
                temp=Path(stream.name); size=0
                for chunk in response.iter_content(1024*1024):
                    size+=len(chunk)
                    if size>80*1024*1024: raise ValueError('cutout exceeds 80 MiB')
                    stream.write(chunk)
        read_cutout(temp); os.replace(temp,path); temp=None
    finally:
        if temp is not None: temp.unlink(missing_ok=True)
    return path


def _run(*, ra, dec, output, max_pairs=3, epoch_days=1, target_id=None):
    if not (np.isfinite(ra) and 0<=ra<360 and np.isfinite(dec) and -90<=dec<=90): raise ValueError('invalid coordinates')
    if not 1<=max_pairs<=10: raise ValueError('max_pairs must be 1..10')
    from pair_discovery import discover_pairs
    started=time.perf_counter(); out=Path(output); out.mkdir(parents=True,exist_ok=True)
    discovery=discover_pairs(ra,dec,target_id=target_id,epoch_days=epoch_days)
    (out/'discovery.json').write_text(json.dumps(discovery,indent=2),encoding='utf-8')
    from response_calibration import fetch_bundle, pixel_response, response_metrics
    bundle=fetch_bundle(Path(__file__).parent/'artifacts/response_calibration')
    report={'status':'PARTIAL','coordinates':[ra,dec],'attempts':[], 'python':platform.python_version(), 'limits':{'max_pairs':max_pairs}, 'limitations':['Fixed sky grid; target propagation must be requested separately.','No full calibrated spectral response or PSF equivalence established; real differences are not authorized.']}
    source_hash=hashlib.sha256(b''.join(Path(__file__).with_name(name).read_bytes() for name in ('scientific_pipeline.py','psf_processing.py','fits_pipeline.py','wcs_geometry.py','response_calibration.py','spacecraft_ephemeris.py','astrometry_model.py'))).hexdigest()
    configuration={'ra':ra,'dec':dec,'target_id':target_id,'epoch_days':epoch_days,'size':0.05,'allow_source':True,'software_hash':source_hash}
    for index,pair in enumerate(discovery['pairs'][:max_pairs]):
        pair_id=hashlib.sha256(json.dumps([pair['a'],pair['b'],configuration,index],sort_keys=True).encode()).hexdigest()
        checkpoint=out/(pair_id+'.json')
        if checkpoint.exists():
            saved=json.loads(checkpoint.read_text(encoding='utf-8'))
            if saved.get('status')!='FAILED' and saved.get('artifact_sha256') and all((out/name).is_file() and hashlib.sha256((out/name).read_bytes()).hexdigest()==digest for name,digest in saved['artifact_sha256'].items()):
                report['attempts'].append(saved)
                continue
        attempt={'status':'FAILED'}
        try:
            items=[]
            for record in (pair['a'],pair['b']):
                path=download_cutout(record,ra,dec,out/'science_cache')
                science=read_cutout(path,allow_source=True)
                observer_astrometry={'status':'fixed sky coordinates'}
                if target_id:
                    try:
                        from spacecraft_ephemeris import observer_vector
                        from astrometry_model import propagate_direction
                        from targets import targets
                        vector,observer_provenance=observer_vector(science['header']['DATE-OBS'],out/'observer_cache')
                        target=next(t for t in targets() if t['id']==target_id)
                        observer_astrometry=propagate_direction(target,science['header']['DATE-OBS'],observer_position_au=vector,observer_provenance=observer_provenance)
                        observer_astrometry['usage']='Post-download geometric coverage diagnostic; initial discovery still uses approximate barycentric astrometry.'
                        tx,ty=science['wcs'].celestial.world_to_pixel_values(observer_astrometry['ra_deg'],observer_astrometry['dec_deg'])
                        observer_astrometry['inside_cutout']=bool(np.isfinite([tx,ty]).all() and -.5<=tx<science['image'].shape[1]-.5 and -.5<=ty<science['image'].shape[0]-.5)
                    except (ValueError,OSError,KeyError,requests.RequestException) as exc:
                        observer_astrometry={'status':'UNAVAILABLE','reason':str(exc)}
                pid=record.get('obs_publisher_did') or f"ivo://irsa.ipac/spherex_dr1?{record['obs_id']}/D{science['header']['DETECTOR']}"
                wave,band,valid,provenance=wavelength_maps(path,pid,out/'calibration_cache')
                from psf_processing import extract_epsf
                try:
                    psfx,psfy=science['wcs'].celestial.world_to_pixel_values(ra,dec)
                    psf=extract_epsf(path,float(psfx),float(psfy))
                except (ValueError,KeyError,OSError) as exc:
                    psf={'provenance':{'status':'UNAVAILABLE','reason':str(exc)}}
                science.update(wave=wave,band=band,calibration_valid=valid,psf=psf)
                with fits.open(path) as hdus:
                    zodi=np.asarray(hdus['ZODI'].data,dtype=float) if 'ZODI' in hdus else None
                blank=(science['flags']==0) & science['valid']
                background={'blank_pixels':int(blank.sum()),'method':'unsubtracted native source-free median diagnostic'}
                if blank.any():
                    background['median_mjy_sr']=float(np.median(science['image'][blank]))
                    background['mad_mjy_sr']=float(np.median(np.abs(science['image'][blank]-np.median(science['image'][blank]))))
                if zodi is not None and np.isfinite(zodi).any(): background['modeled_zodi_median_mjy_sr']=float(np.nanmedian(zodi))
                science['background_diagnostic']=background
                response_info={'status':'UNAVAILABLE'}
                try:
                    px,py=science['wcs'].celestial.world_to_pixel_values(ra,dec)
                    if not np.isfinite([px,py]).all(): raise ValueError('requested sky coordinate has invalid detector projection')
                    cx,cy=int(np.rint(px)),int(np.rint(py))
                    if not (0<=cx<science['image'].shape[1] and 0<=cy<science['image'].shape[0]): raise ValueError('requested sky coordinate outside cutout')
                    rw,rr,rp=pixel_response(bundle,science['header'],science['image'].shape,int(cx),int(cy),float(wave[cy,cx]))
                    np.savetxt(out/f'pair_{index}_{len(items)}_response.csv',np.column_stack([rw,rr]),delimiter=',',header='wavelength_um,measured_lab_response')
                    science['response_curve']=(rw,rr)
                    response_info={**rp,'requested_sky_deg':[ra,dec],'sampling':'nearest native pixel to common sky position; not full response at interpolated sky pixel'}
                except (ValueError,KeyError,OSError) as exc: response_info={'status':'UNAVAILABLE','reason':str(exc)}
                science['response_info']=response_info
                items.append(science)
                attempt.setdefault('observations',[]).append({'obs_id':record['obs_id'],'date':science['header'].get('DATE-OBS'),'observer_relative_astrometry':observer_astrometry,'response':response_info,'psf':psf['provenance'],'background':background,'calibration':provenance,'mask_policy':'zero or SOURCE-only; all other bits excluded; contamination remains possible','quality_valid_fraction':float(science['valid'].mean()),'wavelength_range_um':[float(np.nanmin(wave)),float(np.nanmax(wave))]})
            # Preserve calibrated intermediate products even when comparison cannot be authorized.
            for j,item in enumerate(items):
                header=item['wcs'].celestial.to_header(relax=True)
                hdus=[fits.PrimaryHDU()]
                for name,key,unit in [('IMAGE','image','MJy/sr'),('CWAVE','wave','um'),('CBAND','band','um')]:
                    h=header.copy(); h['BUNIT']=unit; hdus.append(fits.ImageHDU(item[key],header=h,name=name))
                fits.HDUList(hdus).writeto(out/f'pair_{index}_{j}_calibrated.fits',overwrite=True)
            from wcs_geometry import reproject_bilinear
            a,b=items
            aligned=reproject_bilinear(b['image'],b['variance'],b['valid'],b['wcs'].celestial,a['wcs'].celestial,a['image'].shape)
            spectral=reproject_bilinear(b['wave'],np.zeros_like(b['wave']),b['calibration_valid'],b['wcs'].celestial,a['wcs'].celestial,a['image'].shape)
            overlap=aligned['footprint'] & a['calibration_valid'] & spectral['valid']
            attempt['geometric_overlap_fraction']=float(aligned['footprint'].mean())
            attempt['joint_quality_fraction']=float((aligned['valid'] & a['valid']).mean())
            attempt['median_absolute_wavelength_delta_um']=float(np.median(np.abs(spectral['image'][overlap]-a['wave'][overlap]))) if overlap.any() else None
            header=a['wcs'].celestial.to_header(relax=True)
            fits.HDUList([fits.PrimaryHDU(),fits.ImageHDU(aligned['image'],header=header,name='ALIGNED_B'),fits.ImageHDU(aligned['variance'],header=header,name='VARIANCE_B'),fits.ImageHDU(overlap.astype('uint8'),header=header,name='OVERLAP')]).writeto(out/f'pair_{index}_alignment.fits',overwrite=True)
            if 'response_curve' in a and 'response_curve' in b:
                try: attempt['fiducial_response_overlap']=response_metrics(*a['response_curve'],*b['response_curve'])
                except ValueError as exc: attempt['fiducial_response_overlap']={'status':'NOT_PHYSICAL_THROUGHPUT','reason':str(exc),'negative_lab_values_preserved':True}
            if 'kernel' in a['psf'] and 'kernel' in b['psf']:
                try:
                    from psf_processing import project_kernel
                    from astropy.coordinates import SkyCoord
                    projected=project_kernel(b['psf']['kernel'],b['wcs'].celestial,a['wcs'].celestial,SkyCoord(ra,dec,unit='deg'))
                    attempt['psf_projection']={key:value for key,value in projected.items() if key!='kernel'}
                    fits.HDUList([fits.PrimaryHDU(),fits.ImageHDU(a['psf']['kernel'],name='PSF_A'),fits.ImageHDU(projected['kernel'],name='PSF_B_ON_A')]).writeto(out/f'pair_{index}_psf.fits',overwrite=True)
                except ValueError as exc: attempt['psf_projection']={'status':'UNAVAILABLE','reason':str(exc)}
            attempt['gates']={'spatial_overlap':bool(overlap.any()),'joint_quality':bool((aligned['valid'] & a['valid']).any()),'full_response_equivalence':'UNVERIFIED','psf_common_sky_orientation':'UNVERIFIED','background_correction':'NOT_APPLIED'}
            attempt.update(status='REJECTED' ,reason='CWAVE/CBAND do not establish full response equivalence; PSF/background consistency remains unverified. No scientific subtraction authorized.')
        except Exception as exc:
            attempt.update(error_type=type(exc).__name__,reason=str(exc))
        report['attempts'].append(attempt)
        attempt['artifact_sha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in list(out.glob(f'pair_{index}_*.fits'))+list(out.glob(f'pair_{index}_*.csv'))}
        temporary=checkpoint.with_suffix('.tmp')
        temporary.write_text(json.dumps(attempt,indent=2),encoding='utf-8'); os.replace(temporary,checkpoint)
        (out/'report.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    report['elapsed_seconds']=time.perf_counter()-started
    report['output_bytes']=sum(p.stat().st_size for p in out.rglob('*') if p.is_file())
    from importlib.metadata import version
    report['dependencies']={name:version(name) for name in ('astropy','numpy','requests','fastapi','matplotlib')}
    report['configuration_sha256']=hashlib.sha256(json.dumps({'ra':ra,'dec':dec,'max_pairs':max_pairs,'epoch_days':epoch_days,'target_id':target_id},sort_keys=True).encode()).hexdigest()
    report['evaluation_ranking']=sorted(range(len(report['attempts'])),key=lambda i:(report['attempts'][i].get('median_absolute_wavelength_delta_um') is None,report['attempts'][i].get('median_absolute_wavelength_delta_um') or 0,-report['attempts'][i].get('joint_quality_fraction',0)))
    report['search_exhaustive']=False
    report['returned_candidates']=len(discovery['pairs'])
    (out/'report.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    return report


def run(**kwargs):
    """Exclusive output-directory writer; OS lock is released even on process exit."""
    out=Path(kwargs['output']);out.mkdir(parents=True,exist_ok=True)
    with open(out/'.run.lock','a+b') as lock:
        lock.seek(0); lock.write(b'0'); lock.flush();lock.seek(0)
        if os.name=='nt':
            import msvcrt
            try: msvcrt.locking(lock.fileno(),msvcrt.LK_NBLCK,1)
            except OSError as exc: raise ValueError('output directory is busy') from exc
        else:
            import fcntl
            fcntl.flock(lock.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
        try: return _run(**kwargs)
        finally:
            lock.seek(0)
            if os.name=='nt': msvcrt.locking(lock.fileno(),msvcrt.LK_UNLCK,1)
            else: fcntl.flock(lock.fileno(),fcntl.LOCK_UN)


def main():
    p=argparse.ArgumentParser(); group=p.add_mutually_exclusive_group(required=True); group.add_argument('--ra',type=float); group.add_argument('--target'); p.add_argument('--dec',type=float); p.add_argument('--output',default='artifacts/comparison_run'); p.add_argument('--max-pairs',type=int,default=3)
    a=p.parse_args()
    if a.target:
        from targets import barycentric_position
        c=barycentric_position(a.target,'2025-07-01'); a.ra,a.dec=c['ra_deg'],c['dec_deg']
    if a.dec is None: p.error('--dec is required with --ra')
    print(json.dumps(run(ra=a.ra,dec=a.dec,output=a.output,max_pairs=a.max_pairs,target_id=a.target),indent=2))
if __name__=='__main__': main()


def survey(target_ids, output='artifacts/survey', pairs_per_target=1):
    """Sequential multi-target search with explicit total budget of ten pairs."""
    from targets import barycentric_position
    ids=list(target_ids)
    if not ids or len(ids)*pairs_per_target>10 or not 1<=pairs_per_target<=10: raise ValueError('survey exceeds ten-pair budget')
    directory=Path(output);directory.mkdir(parents=True,exist_ok=True);results=[]
    for target in ids:
        center=barycentric_position(target,'2025-07-01')
        try:
            result=run(ra=center['ra_deg'],dec=center['dec_deg'],target_id=target,max_pairs=pairs_per_target,output=directory/target)
            results.append({'target':target,'status':result['status'],'report':str(directory/target/'report.json')})
        except (ValueError,OSError,requests.RequestException) as exc:
            results.append({'target':target,'status':'FAILED','reason':str(exc)})
        (directory/'survey.json').write_text(json.dumps(results,indent=2),encoding='utf-8')
    return results
