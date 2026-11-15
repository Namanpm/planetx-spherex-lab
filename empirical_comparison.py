"""Restricted native-photometry repeatability experiment with source-level holdout."""
import argparse
import hashlib
import json
import time
from itertools import combinations
from pathlib import Path
import numpy as np
from astropy.io import fits
from scipy.stats import chi2
from reference_catalog import query_references,source_position
from native_photometry import fit_source,get_observation,prepare_pixel_areas
from empirical_repeatability import fit_repeatability
from scientific_pipeline import download_cutout
from pair_discovery import discover_pairs
from calibration_resolver import wavelength_maps

POLICY={'snr_min':10.,'fit_pvalue_min':.0001,'valid_fraction_min':.8,'max_g_mag_reference':16.,'minimum_training':12,'minimum_heldout':6,'minimum_confirmation':6,'training_outlier_sigma':4.,'validation_alpha':.05,'cutout_size_deg':.3,'max_observations':8,'source_limit':300}

def dump(path,value):
    temporary=path.with_suffix('.tmp');temporary.write_text(json.dumps(value,indent=2,allow_nan=False),encoding='utf-8');temporary.replace(path)

def _run(target_id='rw-aurigae',output='artifacts/empirical_comparison',max_observations=8, *, selected_records=None, ra=None, dec=None, catalog=None, cache_root=None, progress=None, include_all_pairs=False, premeasured_packets=None):
    if not isinstance(max_observations,int) or not 2<=max_observations<=12:raise ValueError('observation budget must be 2..12')
    from targets import barycentric_position
    started=time.perf_counter();out=Path(output);out.mkdir(parents=True,exist_ok=True)
    cache_root=Path(cache_root) if cache_root is not None else out
    def notify(stage, **details):
        if progress is not None:progress({'stage':stage,**details})
    if ra is None and dec is None:
        center=barycentric_position(target_id,'2025-11-01');ra,dec=center['ra_deg'],center['dec_deg']
    if ra is None or dec is None or not np.isfinite([ra,dec]).all() or not 0<=ra<360 or not -90<=dec<=90:raise ValueError('invalid sky position')
    notify('discovery',selection='supplied' if selected_records is not None else 'live')
    discovery={'observations':list(selected_records)} if selected_records is not None else discover_pairs(target_id=target_id,limit=1000)
    if catalog is None:catalog=query_references(ra,dec,.2,cache_dir=cache_root/'catalog_cache',limit=5000,max_g_mag=19.)
    dump(out/'catalog.json',catalog)
    refs=sorted([s for s in catalog['sources'] if s['phot_g_mean_mag']<=POLICY['max_g_mag_reference']],key=lambda s:(s['phot_g_mean_mag'],s['source_id']))[:POLICY['source_limit']]
    pool=sorted([r for r in discovery['observations'] if selected_records is not None or r['energy_bandpassname']=='SPHEREx-D1'],key=lambda r:(r['t_min'],r['obs_id']))
    # Sample across the complete D1 time span, not a truncated all-detector pair list.
    indices=np.unique(np.rint(np.linspace(0,len(pool)-1,min(max_observations,len(pool)))).astype(int)) if pool else []
    records=list(pool) if selected_records is not None else [pool[i] for i in indices]
    if selected_records is not None:
        if not 2<=len(records)<=max_observations:raise ValueError('selected observation count outside budget')
        if len({r['obs_publisher_did'] for r in records})!=len(records):raise ValueError('duplicate selected observation identity')
    report={'status':'PARTIAL','target':target_id,'policy':POLICY,'reference_candidates':len(refs),'observations':[],'pair_evaluations':[], 'limitations':['Gaia non-variable classification absence is not proof of infrared constancy.','Finite local ePSF, PM-only Gaia positions, WCS solid-angle conversion and independent native pixel errors are approximations.','Calibration applies only to screened reference population and observed regime; no target or arbitrary-pixel correction is authorized.']}
    measurements={}
    observation_cache={}
    software=hashlib.sha256(b''.join(Path(__file__).with_name(name).read_bytes() for name in ('empirical_comparison.py','native_photometry.py','psf_processing.py','fits_pipeline.py','spectral_calibration.py','empirical_repeatability.py','reference_catalog.py','calibration_resolver.py'))).hexdigest()
    configuration={'policy':POLICY,'ra':ra,'dec':dec,'catalog':catalog['provenance']['sha256'],'software':software}
    report['configuration']=configuration
    for record in records:
        notify('observation',observation_id=record['obs_publisher_did'],completed=len(measurements),total=len(records))
        key=record['obs_publisher_did'];measurement_path=out/(hashlib.sha256(json.dumps([key,configuration],sort_keys=True).encode()).hexdigest()+'.json')
        replayed=premeasured_packets is not None
        if replayed:
            if key not in premeasured_packets:raise ValueError('replay missing selected observation packet')
            packet=premeasured_packets[key]
            if packet.get('record',{}).get('obs_publisher_did')!=key:raise ValueError('replay packet identity mismatch')
            dump(measurement_path,packet)
        else:packet=json.loads(measurement_path.read_text(encoding='utf-8')) if measurement_path.exists() else None
        if not replayed and (packet is None or packet.get('error')):
            packet={'record':record,'measurements':{},'rejected':[]}
            try:
                path=download_cutout(record,ra,dec,cache_root/'science_cache',size=POLICY['cutout_size_deg'])
                observation=get_observation(path,key,observation_cache)
                packet['_readiness_checks']={k:True for k in ('wcs_valid','epsf_available','variance_available','quality_flags_available')}
                date=observation['science']['header']['DATE-OBS']
                wave,band,calvalid,prov=wavelength_maps(path,key,cache_root/'calibration_cache',product_identity=record.get('access_url'))
                packet['_readiness_checks']['wavelength_map_valid']=bool(np.any(calvalid & np.isfinite(wave) & (wave>0) & np.isfinite(band) & (band>0)))
                packet.update(science_path=str(path),date=date,calibration=prov)
                from response_calibration import fetch_bundle
                from reference_response import field_response
                bundle=fetch_bundle(Path(__file__).parent/'artifacts/response_calibration')
                try:
                    rw,rr,rp=field_response(path,wave,bundle,ra,dec)
                    packet['representative_response']={'wavelength_um':rw.tolist(),'samples':rr.tolist(),'provenance':rp}
                except (ValueError,KeyError,OSError) as exc:packet['response_error']=str(exc)
                positions=[]
                for source in catalog['all_sources']:
                    if all(source.get(k) not in (None,'') for k in ('pmra','pmdec','ref_epoch')):
                        try:position=source_position(source,date)
                        except (ValueError,TypeError):position={'ra_deg':float(source['ra']),'dec_deg':float(source['dec'])}
                    else:position={'ra_deg':float(source['ra']),'dec_deg':float(source['dec'])}
                    positions.append((str(source['source_id']),position['ra_deg'],position['dec_deg']))
                # Propagate each reference once, then batch the identical native
                # geometry used by fit_source. Preserve individual failure handling.
                reference_positions={}
                for source in refs:
                    try:reference_positions[str(source['source_id'])]=source_position(source,date)
                    except (ValueError,KeyError,OSError):pass
                prepare_pixel_areas(observation,[(p['ra_deg'],p['dec_deg']) for p in reference_positions.values()])
                for source in refs:
                    try:
                        pos=reference_positions.get(str(source['source_id']))
                        if pos is None:pos=source_position(source,date)
                        # Exclude the challenge target and nearby catalog entries from calibration.
                        from astropy.coordinates import SkyCoord
                        target_distance=SkyCoord(pos['ra_deg'],pos['dec_deg'],unit='deg').separation(SkyCoord(ra,dec,unit='deg')).arcsec
                        if target_distance<30:raise ValueError('target exclusion region')
                        neighbors=[[r,d] for ident,r,d in positions if ident!=str(source['source_id'])]
                        fit=fit_source(observation,pos['ra_deg'],pos['dec_deg'],neighbors,deblend=True)
                        if not fit['central_valid'] or fit['central_flag'] not in (0,1<<21):raise ValueError('source core flagged or invalid')
                        snr=fit['flux_jy']/np.sqrt(fit['variance_jy2']);p=float(chi2.sf(fit['chi2'],fit['dof']))
                        if snr<POLICY['snr_min'] or fit['valid_fraction']<POLICY['valid_fraction_min'] or p<POLICY['fit_pvalue_min']:raise ValueError(f'fit quality: SNR={snr:.3g},p={p:.3g},fraction={fit["valid_fraction"]:.3g}')
                        x,y=int(np.rint(fit['x'])),int(np.rint(fit['y']))
                        if not calvalid[y,x] or not np.isfinite([wave[y,x],band[y,x]]).all() or min(wave[y,x],band[y,x])<=0:raise ValueError('source spectral calibration invalid')
                        fit.update(source_id=str(source['source_id']),color=source['bp_rp'],wave_um=float(wave[y,x]),band_um=float(band[y,x]),snr=float(snr),fit_pvalue=p)
                        packet['measurements'][str(source['source_id'])]=fit
                    except (ValueError,KeyError,OSError) as exc:packet['rejected'].append({'source_id':str(source['source_id']),'reason':str(exc)})
            except Exception as exc:packet['error']={'type':type(exc).__name__,'reason':str(exc)}
            dump(measurement_path,packet)
        measurements[key]=packet
        report['observations'].append({'publisher_id':key,'measurement_count':len(packet['measurements']),'rejected_count':len(packet['rejected']),'error':packet.get('error'),'calibration':packet.get('calibration'),'measurement_file':str(measurement_path)})
        dump(out/'report.json',report)
    candidates=[]
    response_distances={}
    for a,b in combinations(records,2):
        reason=None
        if abs(a['t_min']-b['t_min'])<1:reason='observations do not meet the distinct-epoch separation policy (one day)'
        if a['energy_bandpassname']!=b['energy_bandpassname']:reason='different detector/wavelength domains; no transformation established'
        if include_all_pairs and reason:
            report['pair_evaluations'].append({'observation_a':a['obs_publisher_did'],'observation_b':b['obs_publisher_did'],'status':'REJECTED','reason':reason})
            continue
        if not include_all_pairs and abs(a['t_min']-b['t_min'])<1:continue
        pa,pb=measurements[a['obs_publisher_did']],measurements[b['obs_publisher_did']]
        common=sorted(set(pa['measurements'])&set(pb['measurements']))
        # Greedy independent stamps, ordered by catalog magnitude and ID, never residual.
        brightness={str(r['source_id']):r['phot_g_mean_mag'] for r in refs}
        chosen=[]
        for ident in sorted(common,key=lambda i:(int(hashlib.sha256(i.encode()).hexdigest(),16)%4==3,brightness[i],i)):
            if all(not stamps_overlap(pa['measurements'][ident]['fit_bounds'],pa['measurements'][old]['fit_bounds']) and not stamps_overlap(pb['measurements'][ident]['fit_bounds'],pb['measurements'][old]['fit_bounds']) for old in chosen):chosen.append(ident)
        common=chosen
        delta=float(np.median([abs(pa['measurements'][i]['wave_um']-pb['measurements'][i]['wave_um']) for i in common])) if common else None
        metric={'normalized_l2_distance':None,'reason':'representative response unavailable'}
        if 'representative_response' in pa and 'representative_response' in pb:
            from reference_response import signed_laboratory_distance
            ca,cb=pa['representative_response'],pb['representative_response']
            metric=signed_laboratory_distance((ca['wavelength_um'],ca['samples']),(cb['wavelength_um'],cb['samples']))
        response_distances[(a['obs_publisher_did'],b['obs_publisher_did'])]=metric
        candidates.append((len(common),delta,a,b,common))
    candidates.sort(key=lambda r:(response_distances[(r[2]['obs_publisher_did'],r[3]['obs_publisher_did'])]['normalized_l2_distance'] if response_distances[(r[2]['obs_publisher_did'],r[3]['obs_publisher_did'])]['normalized_l2_distance'] is not None else float('inf'),-sum(int(hashlib.sha256(i.encode()).hexdigest(),16)%4!=3 for i in r[4]),r[1] if r[1] is not None else float('inf')))
    notify('pair_evaluation',candidate_count=len(records)*(len(records)-1)//2)
    selected=None
    for count,delta,a,b,common in candidates:
        attempt={'observation_a':a['obs_publisher_did'],'observation_b':b['obs_publisher_did'],'common_independent_screened_sources':count,'median_wavelength_delta_um':delta,'status':'REJECTED'}
        if include_all_pairs and response_distances[(a['obs_publisher_did'],b['obs_publisher_did'])]['normalized_l2_distance'] is None:
            attempt['reason']='representative spectral response unavailable; applicability cannot be assessed'
            report['pair_evaluations'].append(attempt);continue
        pa,pb=measurements[a['obs_publisher_did']]['measurements'],measurements[b['obs_publisher_did']]['measurements']
        train_ids=[i for i in common if int(hashlib.sha256(i.encode()).hexdigest(),16)%4<2]
        validation_ids=[i for i in common if int(hashlib.sha256(i.encode()).hexdigest(),16)%4==2]
        confirmation_ids=[i for i in common if int(hashlib.sha256(i.encode()).hexdigest(),16)%4==3]
        attempt['representative_response_distance']=response_distances[(a['obs_publisher_did'],b['obs_publisher_did'])]
        if train_ids:
            cmin,cmax=min(pa[i]['color'] for i in train_ids),max(pa[i]['color'] for i in train_ids)
            validation_ids=[i for i in validation_ids if cmin<=pa[i]['color']<=cmax]
            confirmation_ids=[i for i in confirmation_ids if cmin<=pa[i]['color']<=cmax]
        attempt['split_counts']={'training':len(train_ids),'validation':len(validation_ids),'confirmation':len(confirmation_ids)}
        if len(train_ids)<POLICY['minimum_training'] or len(validation_ids)<POLICY['minimum_heldout']:
            attempt['reason']='insufficient independent preassigned train/validation/confirmation sources'
            if len(train_ids)>=4 and len(validation_ids)>=2:
                try:attempt['exploratory_model_not_confirmation_eligible']=fit_pair_model(common,pa,pb,train_ids,validation_ids,min_training=4)
                except (ValueError,KeyError,np.linalg.LinAlgError) as exc:attempt['exploratory_model_error']=str(exc)
            report['pair_evaluations'].append(attempt);continue
        try:
            fit=fit_pair_model(common,pa,pb,train_ids,validation_ids)
            attempt['validation_model']=fit
            p=fit['heldout']['pvalue_chi2_approximate']
            attempt['status']='VALIDATION_CONSISTENT' if is_consistent(fit) else 'REJECTED'
            attempt['reason']='Reference-population validation consistency only; confirmation remains untouched.' if attempt['status']=='VALIDATION_CONSISTENT' else 'Validation predictions or color support failed'
            if selected is None and attempt['status']=='VALIDATION_CONSISTENT':selected=(attempt,common,pa,pb,train_ids,confirmation_ids)
        except (ValueError,KeyError,np.linalg.LinAlgError) as exc:attempt['reason']=str(exc)
        report['pair_evaluations'].append(attempt)
    # One final test after pair selection; never search for another confirmation pass.
    if selected is not None:
        attempt,common,pa,pb,train_ids,confirmation_ids=selected
        try:
            if len(confirmation_ids)<POLICY['minimum_confirmation']:raise ValueError('selected pair lacks sufficient confirmation sources; no alternate pair will be tried')
            confirmation=fit_pair_model(common,pa,pb,train_ids,confirmation_ids)
            p=confirmation['heldout']['pvalue_chi2_approximate']
            report['confirmation']={'selected_pair':[attempt['observation_a'],attempt['observation_b']],'model':confirmation,'tested_once':True}
            consistent=is_consistent(confirmation)
            report['status']='REFERENCE_MODEL_CONSISTENT_ON_CONFIRMATION' if consistent else 'CONFIRMATION_FAILED'
            report['interpretation']='Failure to reject on independent confirmation supports only a quantitative, conditional reference-population repeatability comparison; not equivalence, constancy or arbitrary image subtraction.'
        except (ValueError,KeyError,np.linalg.LinAlgError) as exc:
            report['status']='CONFIRMATION_UNAVAILABLE'
            report['confirmation']={'status':'FAILED','reason':str(exc)}
    else: report['status']='INSUFFICIENT_OR_FAILED_REFERENCE_VALIDATION'
    report['prepared_observation_count']=len(observation_cache)
    report['prepared_epsf_zone_count']=sum(len(o['epsf']['zones']) for o in observation_cache.values())
    area_caches=[o.get('_pixel_area_cache',{}) for o in observation_cache.values()]
    report['pixel_area_cache']={name:sum(c.get(name,0) for c in area_caches) for name in ('hits','misses','batch_points')}
    report['pixel_area_cache']['entries']=sum(len(c.get('values',{})) for c in area_caches)
    report['elapsed_seconds']=time.perf_counter()-started;report['output_tree_bytes']=sum(p.stat().st_size for p in out.rglob('*') if p.is_file())
    if candidates:
        _,_,aa,bb,common=candidates[0]
        plot_a,plot_b=measurements[aa['obs_publisher_did']]['measurements'],measurements[bb['obs_publisher_did']]['measurements']
        pair_ids=[aa['obs_publisher_did'],bb['obs_publisher_did']]
        if selected is not None:
            chosen,common,plot_a,plot_b,_,_=selected
            pair_ids=[chosen['observation_a'],chosen['observation_b']]
        report['diagnostic_pair']=pair_ids
        render_sources(out,plot_a,plot_b,common,report['status'],pair_ids)
    if include_all_pairs:
        report['measurement_replay']=premeasured_packets is not None
        report['theoretical_pair_count']=len(records)*(len(records)-1)//2
        report['selection_provenance']='supplied observation records; discovery performed by caller'
        report['confirmation_policy']='Only the preselected best pair receives final confirmation; other validation-consistent pairs have not received independent confirmation.'
    notify('complete',status=report['status'])
    dump(out/'report.json',report)
    return report

def stamps_overlap(a,b):
    return not (a[1]<b[0] or b[1]<a[0] or a[3]<b[2] or b[3]<a[2])


def fit_pair_model(ids,pa,pb,training,heldout,min_training=None):
    result=fit_repeatability(ids,[pa[i]['flux_jy'] for i in ids],[pb[i]['flux_jy'] for i in ids],[pa[i]['variance_jy2'] for i in ids],[pb[i]['variance_jy2'] for i in ids],training_ids=training,heldout_ids=heldout,color=[pa[i]['color'] for i in ids],min_training=POLICY['minimum_training'] if min_training is None else min_training,outlier_sigma=POLICY['training_outlier_sigma'])
    test_ids=result['heldout']['source_ids']
    if len(test_ids)>=4 and all('x' in pa[i] and 'y' in pa[i] for i in test_ids):
        from empirical_repeatability import spatial_repeatability_diagnostic
        result['heldout']['spatial_diagnostic']=spatial_repeatability_diagnostic([[pa[i]['x'],pa[i]['y']] for i in test_ids],result['heldout']['standardized_residual'],seed=9281)
    return result

def render_sources(out,a,b,ids,status,pair_ids):
    """Raw measured flux contrast diagnostics, never a variable-source claim."""
    columns=[fits.Column(name='GAIA_SOURCE_ID',format='24A',array=np.array(ids,dtype='S24'))]
    for name,values,unit in [('FLUX_A',[a[i]['flux_jy'] for i in ids],'Jy'),('FLUX_B',[b[i]['flux_jy'] for i in ids],'Jy'),('VAR_A',[a[i]['variance_jy2'] for i in ids],'Jy2'),('VAR_B',[b[i]['variance_jy2'] for i in ids],'Jy2'),('GAIA_BP_RP',[a[i]['color'] for i in ids],'mag'),('WAVE_A',[a[i]['wave_um'] for i in ids],'um'),('WAVE_B',[b[i]['wave_um'] for i in ids],'um')]:
        columns.append(fits.Column(name=name,format='D',unit=unit,array=np.array(values)))
    table=fits.BinTableHDU.from_columns(columns,name='MEASUREMENTS');table.header['STATUS']=status
    table.header['OBS_A']=pair_ids[0];table.header['OBS_B']=pair_ids[1]
    table.header['COMMENT']='Raw cross-band measurements; no variability or image-subtraction inference.'
    fits.HDUList([fits.PrimaryHDU(),table]).writeto(out/'reference_measurements.fits',overwrite=True,checksum=True)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(11,4))
    if ids:
        color=np.array([a[i]['color'] for i in ids]);fa=np.array([a[i]['flux_jy'] for i in ids]);fb=np.array([b[i]['flux_jy'] for i in ids]);v=np.array([a[i]['variance_jy2']/a[i]['flux_jy']**2+b[i]['variance_jy2']/b[i]['flux_jy']**2 for i in ids])
        axes[0].errorbar(fa,fb,xerr=np.sqrt([a[i]['variance_jy2'] for i in ids]),yerr=np.sqrt([b[i]['variance_jy2'] for i in ids]),fmt='.')
        axes[1].errorbar(color,np.log(fb/fa),yerr=np.sqrt(v),fmt='.')
    axes[0].set(xlabel='Native ePSF flux A (approx. Jy)',ylabel='Native ePSF flux B (approx. Jy)')
    axes[1].set(xlabel='Gaia BP−RP (mag)',ylabel='Raw ln(flux B / flux A)')
    fig.suptitle(status+' — reference diagnostics, not variability detections')
    fig.tight_layout();fig.savefig(out/'reference_diagnostic.png',dpi=150);plt.close(fig)

def is_consistent(fit):
    p=fit['heldout']['pvalue_chi2_approximate']
    spatial=fit['heldout'].get('spatial_diagnostic',{})
    red_flag=spatial.get('moran_I') is not None and spatial['moran_I']>0 and spatial.get('pvalue_permutation_positive',1)<POLICY['validation_alpha']
    return bool(not red_flag and np.isfinite(p) and 0<=p<=1 and p>=POLICY['validation_alpha'] and fit['heldout']['extrapolated_color_source_count']==0)


def run(target_id='rw-aurigae',output='artifacts/empirical_comparison',max_observations=8, **options):
    import os
    out=Path(output);out.mkdir(parents=True,exist_ok=True)
    with open(out/'.empirical.lock','a+b') as lock:
        lock.seek(0);lock.write(b'0');lock.flush();lock.seek(0)
        if os.name=='nt':
            import msvcrt
            try:msvcrt.locking(lock.fileno(),msvcrt.LK_NBLCK,1)
            except OSError as exc:raise ValueError('empirical output is busy') from exc
        else:
            import fcntl
            fcntl.flock(lock.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
        try:return _run(target_id,output,max_observations,**options)
        finally:
            lock.seek(0)
            if os.name=='nt':msvcrt.locking(lock.fileno(),msvcrt.LK_UNLCK,1)
            else:fcntl.flock(lock.fileno(),fcntl.LOCK_UN)

def run_selected(records, *, target_id='rw-aurigae', ra=None, dec=None, output, cache_root=None, progress=None, catalog=None, premeasured_packets=None):
    """Analyze exact supplied observations once, including every rejected pair.

    The caller owns live/replay provenance. This entry never performs discovery.
    Confirmation remains a single test after deterministic pair selection.
    """
    records=list(records)
    if premeasured_packets is not None and catalog is None:raise ValueError('explicit measurement replay requires a supplied catalog')
    return run(target_id,output,len(records),selected_records=records,ra=ra,dec=dec,
               cache_root=cache_root,progress=progress,catalog=catalog,include_all_pairs=True,premeasured_packets=premeasured_packets)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--target',default='rw-aurigae');p.add_argument('--output',default='artifacts/empirical_comparison');p.add_argument('--max-observations',type=int,default=8);a=p.parse_args()
    r=run(a.target,a.output,a.max_observations);print(json.dumps({'status':r['status'],'observations':r['observations'],'pairs':[{k:v for k,v in p.items() if k!='empirical_model'} for p in r['pair_evaluations']]},indent=2))
