"""Multi-observation orchestration with explicit live versus snapshot replay provenance."""
import hashlib
import json
import time
from pathlib import Path
from observation_selection import select_observations, candidate_pairs, validate_coordinates


def snapshot(report_path, selected_ids=None, maximum=8):
    """Load an explicit measured snapshot; never infer replay after network failure."""
    path=Path(report_path)
    report=json.loads(path.read_text(encoding='utf-8'))
    packets={};records=[];input_hashes={}
    raw=[(o,Path(o['measurement_file']),json.loads(Path(o['measurement_file']).read_text(encoding='utf-8'))) for o in report['observations']]
    inspection_ids={r['obs_publisher_did'] for r in select_observations([p['record'] for _,_,p in raw],selected_ids,maximum=maximum)}
    for observation,packet_path,packet in raw:
        input_hashes[str(packet_path)]=hashlib.sha256(packet_path.read_bytes()).hexdigest()
        record=packet['record'];ident=record['obs_publisher_did']
        if ident!=observation['publisher_id']:raise ValueError('Replay packet identity mismatch')
        # Preserve snapshot failures explicitly. A failed epoch is never replaced.
        if not packet.get('error') and not Path(packet.get('science_path','!')).is_file():
            packet={**packet,'error':{'type':'MissingProduct','reason':'Recorded science file missing'}}
            packet['measurements']={}
        if ident in inspection_ids and packet.get('science_path') and Path(packet['science_path']).is_file():
            from native_photometry import get_observation
            try:
                prepared=get_observation(packet['science_path'],ident,{})
                packet['_readiness_checks']={k:True for k in ('wcs_valid','epsf_available','variance_available','quality_flags_available')}
                input_hashes[packet['science_path']]=hashlib.sha256(Path(packet['science_path']).read_bytes()).hexdigest()
                cal=packet.get('calibration',{}) or {}
                cal_url=cal.get('calibration_url')
                if not cal_url:raise ValueError('Replay calibration association missing')
                if cal_url:
                    cal_dir=path.parent/'calibration_cache'
                    cal_path=cal_dir/(hashlib.sha256(cal_url.encode()).hexdigest()+'.fits')
                    if not cal_path.is_file():raise ValueError('Replay calibration product missing; offline replay cannot retrieve it')
                    if cal_path.is_file():
                        from calibration_resolver import DEFAULT_CACHE
                        from spectral_calibration import sample_calibration
                        import numpy as np
                        wave,band,digest=DEFAULT_CACHE.arrays(cal_url,int(ident[-1]),cal_dir)
                        if cal.get('calibration_sha256') and digest!=cal['calibration_sha256']:raise ValueError('Calibration checksum differs from saved measurement provenance')
                        wave,band,valid=sample_calibration(prepared['science']['header'],prepared['science']['image'].shape,wave,band)
                        packet['_readiness_checks']['wavelength_map_valid']=bool(np.any(valid & np.isfinite(wave) & (wave>0) & np.isfinite(band) & (band>0)))
                        input_hashes[str(cal_path)]=digest
                if not packet['_readiness_checks'].get('wavelength_map_valid'):raise ValueError('Replay wavelength map has no valid pixels')
                del prepared
            except (ValueError,KeyError,OSError) as exc:
                packet={**packet,'measurements':{},'error':{'type':type(exc).__name__,'reason':str(exc)}}
        packets[ident]=packet;records.append(record)
    report['_replay_input_sha256']=input_hashes
    return report,records,packets,hashlib.sha256(path.read_bytes()).hexdigest()


def readiness(packet):
    record=packet['record'];measurements=packet.get('measurements',{})
    failure=packet.get('error');cal=packet.get('calibration');checks=packet.get('_readiness_checks',{})
    return {'observation_id':record['obs_publisher_did'],'detector':int(record['obs_publisher_did'][-1]),
            'epoch':packet.get('date'),'mjd':record['t_min'],
            'status':'FAILED' if failure else ('READY' if measurements and all(checks.get(k) for k in ('wcs_valid','epsf_available','variance_available','quality_flags_available','wavelength_map_valid')) else 'PARTIAL'),
            'science_file_available':Path(packet.get('science_path','!')).is_file(),
            'calibration_resolved':bool(cal),'wavelength_map_valid':checks.get('wavelength_map_valid'),
            'wcs_valid':checks.get('wcs_valid'),'epsf_available':checks.get('epsf_available'),
            'variance_available':checks.get('variance_available'),'quality_flags_available':checks.get('quality_flags_available'),
            'validation_evidence':'Replay science FITS inspected directly where checks true; null means no independent inspection recorded',
            'valid_measurement_count':len(measurements),'failure':failure,
            'rejections':packet.get('rejected',[]),'calibration':cal}


def run_analysis(*, target_id='rw-aurigae',ra=None,dec=None,selected_ids=None,mode='best_pair',maximum=8,
                 output,selection_mode='live',replay_report=None,cache_root=None,source_id=None,progress=None):
    from empirical_comparison import run_selected
    from pair_discovery import discover_pairs
    from calibration_resolver import DEFAULT_CACHE
    from benchmark_verification import peak_memory_bytes
    cache_before=DEFAULT_CACHE.cache_info();timings={}
    started=time.perf_counter();out=Path(output);out.mkdir(parents=True,exist_ok=True)
    validate_coordinates(ra,dec)
    from observation_selection import MODES,MAX_OBSERVATIONS,PID
    if mode not in MODES:raise ValueError('Unknown analysis mode')
    if not isinstance(maximum,int) or isinstance(maximum,bool) or not 2<=maximum<=MAX_OBSERVATIONS:raise ValueError('Observation count must be 2..8')
    if selected_ids is not None and (not 2<=len(selected_ids)<=maximum or len(set(selected_ids))!=len(selected_ids) or any(not PID.fullmatch(i) for i in selected_ids)):raise ValueError('Invalid selected observation IDs')
    if selection_mode not in ('live','replay'):raise ValueError('selection_mode must be live or replay')
    def emit(stage,**fields):
        if progress:progress({'stage':stage,**fields})
    emit('observation_selection')
    catalog=None;packets=None;digest=None
    if selection_mode=='replay':
        if replay_report is None:raise ValueError('Explicit replay report required')
        prior,records,packets,digest=snapshot(replay_report,selected_ids,maximum)
        if target_id is not None and prior.get('target') not in (None,target_id):raise ValueError('Replay target does not match snapshot target')
        target_id=prior.get('target',target_id)
        catalog_path=Path(replay_report).parent/'catalog.json'
        catalog=json.loads(catalog_path.read_text(encoding='utf-8'))
        prior['_replay_input_sha256'][str(catalog_path)]=hashlib.sha256(catalog_path.read_bytes()).hexdigest()
        config=prior.get('configuration',{});ra=config.get('ra',ra);dec=config.get('dec',dec)
    else:
        records=discover_pairs(target_id=target_id if ra is None else None,ra=ra,dec=dec,limit=1000)['observations']
    selected=select_observations(records,selected_ids,maximum=maximum,mode=mode)
    timings['selection_seconds']=time.perf_counter()-started
    phase_started=time.perf_counter()
    if packets is not None:
        from source_association import catalog_associations
        initial_associations=catalog_associations(packets,catalog)
        for obs_id,packet in packets.items():
            for ident in list(packet.get('measurements',{})):
                a=initial_associations[obs_id][ident]
                if a['status'] in ('REJECTED','UNDERCONSTRAINED'):
                    packet['measurements'].pop(ident)
                    packet.setdefault('rejected',[]).append({'source_id':ident,'reason':'Association: '+a.get('reason',a['status'])})
    emit('measurement_and_pair_analysis',selected_observations=len(selected))
    empirical=run_selected(selected,target_id=target_id,ra=ra,dec=dec,output=out/'empirical',
                           cache_root=cache_root,catalog=catalog,progress=progress,premeasured_packets=packets)
    used={o['publisher_id']:json.loads(Path(o['measurement_file']).read_text(encoding='utf-8')) for o in empirical['observations']}
    timings['measurement_and_pairs_seconds']=time.perf_counter()-phase_started
    phase_started=time.perf_counter()
    if catalog is None:catalog=json.loads((out/'empirical'/'catalog.json').read_text(encoding='utf-8'))
    from source_association import catalog_associations
    associations=catalog_associations(used,catalog)
    for ident,packet in used.items():
        for source,m in packet.get('measurements',{}).items():
            m['association_status']=associations[ident][source]['status']
    states=[readiness(used[r['obs_publisher_did']]) for r in selected]
    pairs=empirical['pair_evaluations']
    emit('time_series_analysis')
    if mode=='time_series':
        from joint_timeseries import analyze_time_series
        timeline=analyze_time_series(used,source_id=source_id)
    else:timeline={'status':'NOT_REQUESTED'}
    accepted=[p for p in pairs if p['status']=='VALIDATION_CONSISTENT']
    report={'status':'PARTIAL' if any(s['status']!='READY' for s in states) else ('SUCCESS' if accepted else 'INSUFFICIENT_EVIDENCE'),
            'configuration':{'target':target_id,'ra':ra,'dec':dec,'mode':mode,'maximum':maximum,'source_id':source_id},
            'provenance':{'selection_mode':selection_mode,'measurement_replay':packets is not None,'snapshot_sha256':digest,
                          'snapshot_path':str(replay_report) if replay_report else None,
                          'input_sha256':prior.get('_replay_input_sha256',{}) if packets is not None else {}},
            'selected_observation_ids':[r['obs_publisher_did'] for r in selected],
            'source_associations':associations,'readiness':states,'theoretical_pairs':candidate_pairs(selected),'pair_evaluations':pairs,
            'accepted_pair_count':len(accepted),'rejected_pair_count':len(pairs)-len(accepted),
            'best_pair':empirical.get('confirmation'),'time_series':timeline,
            'empirical_report':str(out/'empirical'/'report.json'),
            'elapsed_seconds':time.perf_counter()-started,
            'limitations':['Pair validation consistency does not prove physical response equivalence.',
                            'Snapshot replay reuses recorded measurement outcomes; it does not validate a new processing version.']}
    report['software_sha256']=hashlib.sha256(b''.join(Path(__file__).with_name(n).read_bytes() for n in ('multi_epoch.py','empirical_comparison.py','source_association.py','joint_timeseries.py','calibration_resolver.py'))).hexdigest()
    report['scientific_policy']=empirical['policy']
    report['motion_analysis']={'status':'INSUFFICIENT_EVIDENCE','reason':'Catalog-forced positions are not measured astrometric trajectories. Positional matching returns candidates only.'}
    if report['best_pair'] is None:
        diagnostic=next((p for p in pairs if 'validation_model' in p),pairs[0] if pairs else None)
        report['best_pair']={'status':'REJECTED','diagnostic_candidate':diagnostic,'reason':'No scientifically accepted selected pair'}
    report['pair_acceptance_scope']='Reference validation consistency; final confirmation tested once on selected pair only'
    report['best_pair_confirmation_status']=empirical['status']
    confirmed=empirical.get('confirmation',{});model=confirmed.get('model',{});heldout=model.get('heldout',{})
    if heldout:
        import numpy as np
        selected_pair=confirmed.get('selected_pair')
        selected_evaluation=next((p for p in pairs if [p['observation_a'],p['observation_b']]==selected_pair),{})
        report['selected_pair_summary']={'selected_pair':selected_pair,'split_counts':selected_evaluation.get('split_counts'),
            'rms_log_ratio':heldout['rms_log_ratio'],
            'fractional_residual_rms':float(np.sqrt(np.mean(np.expm1(heldout['residual_log_ratio'])**2))),
            'confirmation_pvalue':heldout['pvalue_chi2_approximate'],'confirmation_coverage_95_fraction':heldout['coverage_95_fraction']}
    if mode=='best_pair' and empirical['status']!='REFERENCE_MODEL_CONSISTENT_ON_CONFIRMATION':report['status']='INSUFFICIENT_EVIDENCE'
    if mode=='time_series' and timeline.get('status')!='SUCCESS':report['status']='INSUFFICIENT_EVIDENCE'
    timings['association_and_time_series_seconds']=time.perf_counter()-phase_started
    report['phase_timings']=timings
    cache_after=DEFAULT_CACHE.cache_info()
    report['calibration_cache']={k:(v-cache_before[k] if k not in ('array_entries','array_bytes','max_bytes') else v) for k,v in cache_after.items()}
    report['peak_process_rss_bytes']=peak_memory_bytes()
    from multi_epoch_products import export
    emit('artifact_generation')
    report['artifacts']=['report.json']+export(report,used,out)
    report['elapsed_seconds']=time.perf_counter()-started
    report['peak_process_rss_bytes']=peak_memory_bytes()
    (out/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False),encoding='utf-8')
    emit('complete',status=report['status'])
    return report


def main():
    import argparse
    parser=argparse.ArgumentParser(description='Selected-observation multi-epoch SPHEREx analysis')
    parser.add_argument('--target',default='rw-aurigae')
    parser.add_argument('--ra',type=float);parser.add_argument('--dec',type=float)
    parser.add_argument('--mode',choices=('best_pair','all_valid_pairs','time_series'),default='best_pair')
    parser.add_argument('--maximum',type=int,default=8)
    parser.add_argument('--selected-id',action='append')
    parser.add_argument('--replay-report',type=Path)
    parser.add_argument('--cache-root',type=Path)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    result=run_analysis(target_id=args.target,ra=args.ra,dec=args.dec,mode=args.mode,
                        maximum=args.maximum,selected_ids=args.selected_id,
                        selection_mode='replay' if args.replay_report else 'live',
                        replay_report=args.replay_report,cache_root=args.cache_root,output=args.output)
    print(json.dumps({'status':result['status'],'observations':len(result['readiness']),
                      'pairs':len(result['pair_evaluations']),'accepted':result['accepted_pair_count'],
                      'seconds':result['elapsed_seconds'],'report':str(args.output/'report.json')},indent=2))


if __name__=='__main__':main()
