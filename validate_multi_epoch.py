"""Reproducible real-data multi-epoch regression and explicitly labeled benchmarks."""
import argparse
import functools
import hashlib
import json
import math
from pathlib import Path
import time

from astropy.io import fits
import empirical_comparison as engine
from benchmark_verification import peak_memory_bytes
from calibration_resolver import DEFAULT_CACHE
from multi_epoch import run_analysis
from observation_selection import select_observations


def read(path):return json.loads(Path(path).read_text(encoding='utf-8'))


def write(path,value):Path(path).write_text(json.dumps(value,indent=2,allow_nan=False),encoding='utf-8')


def pair_key(row):return tuple(sorted((row['observation_a'],row['observation_b'])))


def regression(baseline,actual):
    original={r['publisher_id']:read(r['measurement_file']) for r in baseline['observations']}
    current={r['publisher_id']:read(r['measurement_file']) for r in actual['observations']}
    source_checks={ident:{'measurements_exact':original[ident]['measurements']==current[ident]['measurements'],
                          'rejections_exact':original[ident]['rejected']==current[ident]['rejected'],
                          'source_count':len(current[ident]['measurements'])} for ident in original}
    old_pairs={pair_key(p):p for p in baseline['pair_evaluations']}
    new_pairs={pair_key(p):p for p in actual['pair_evaluations']}
    pair_checks={str(key):old==new_pairs.get(key) for key,old in old_pairs.items()}
    extra=[p for key,p in new_pairs.items() if key not in old_pairs]
    confirmation_exact=baseline.get('confirmation')==actual.get('confirmation')
    selected=actual.get('confirmation',{}).get('selected_pair')
    chosen=next((p for p in actual['pair_evaluations'] if [p['observation_a'],p['observation_b']]==selected),{})
    residuals=actual.get('confirmation',{}).get('model',{}).get('heldout',{}).get('residual_log_ratio',[])
    return {'passed':all(all(v[k] for k in ('measurements_exact','rejections_exact')) for v in source_checks.values()) and all(pair_checks.values()) and confirmation_exact,
            'source_checks':source_checks,'accepted_measurements':sum(v['source_count'] for v in source_checks.values()),
            'prior_pair_records_exact':pair_checks,'additional_pair_records':extra,'confirmation_exact':confirmation_exact,
            'selected_pair':selected,'selected_pair_split_counts':chosen.get('split_counts'),
            'confirmation_fractional_residual_rms':math.sqrt(sum(math.expm1(x)**2 for x in residuals)/len(residuals)) if residuals else None,
            'status':actual['status']}


def artifacts(report,out):
    manifest=read(out/'artifact_manifest.json');checks={}
    for name,entry in manifest.items():
        path=out/name
        valid=path.is_file() and hashlib.sha256(path.read_bytes()).hexdigest()==entry['sha256'] and path.stat().st_size==entry['bytes']
        if path.suffix=='.fits':
            with fits.open(path,checksum=True) as hdus:
                hdus.verify('exception')
                valid=valid and all(h.verify_checksum()==1 and h.verify_datasum()==1 for h in hdus)
        checks[name]=valid
    return {'all_passed':all(checks.values()),'files':checks}


def summary(report,out,elapsed):
    return {'status':report['status'],'wall_seconds_including_outputs':elapsed,
            'peak_process_rss_bytes':peak_memory_bytes(),
            'selected_observation_ids':report['selected_observation_ids'],
            'readiness':[{k:r.get(k) for k in ('observation_id','status','mjd','valid_measurement_count','failure')} for r in report['readiness']],
            'theoretical_pairs':report['theoretical_pairs'],'evaluated_pair_count':len(report['pair_evaluations']),
            'accepted_pair_count':report['accepted_pair_count'],'rejected_pair_count':report['rejected_pair_count'],
            'pair_decisions':[{k:p.get(k) for k in ('observation_a','observation_b','status','reason','split_counts')} for p in report['pair_evaluations']],
            'best_pair':report['best_pair'],'time_series':report['time_series'],
            'phase_timings':report.get('phase_timings'), 'calibration_cache':report.get('calibration_cache'),
            'provenance':report['provenance'],'artifact_integrity':artifacts(report,out)}


def remeasure(records,catalog,baseline,out):
    """Time actual fresh native photometry, restoring every instrumentation hook."""
    stages={};counts={};saved={}
    labels={'download_cutout':'science_access','get_observation':'science_preparation','wavelength_maps':'calibration',
            'prepare_pixel_areas':'pixel_area_preparation','fit_source':'native_photometry','source_position':'position_propagation','fit_repeatability':'empirical_fitting','render_sources':'rendering'}
    for name,label in labels.items():
        original=getattr(engine,name);saved[name]=original
        def make_wrapper(fn,key):
            @functools.wraps(fn)
            def wrapper(*args,**kwargs):
                started=time.perf_counter()
                try:return fn(*args,**kwargs)
                finally:
                    stages[key]=stages.get(key,0.)+time.perf_counter()-started
                    counts[key]=counts.get(key,0)+1
            return wrapper
        setattr(engine,name,make_wrapper(original,label))
    before=DEFAULT_CACHE.cache_info();started=time.perf_counter()
    try:
        actual=engine.run_selected(records,output=out,catalog=catalog,cache_root=BASELINE.parent,
            ra=baseline['configuration']['ra'],dec=baseline['configuration']['dec'])
        elapsed=time.perf_counter()-started
    finally:
        for name,fn in saved.items():setattr(engine,name,fn)
    original={o['publisher_id']:read(o['measurement_file']) for o in baseline['observations']}
    checks={}
    for observation in actual['observations']:
        ident=observation['publisher_id'];packet=read(observation['measurement_file'])
        checks[ident]={'measurements_exact':packet['measurements']==original[ident]['measurements'],'rejections_exact':packet['rejected']==original[ident]['rejected']}
    after=DEFAULT_CACHE.cache_info()
    cache={k:after[k]-before[k] for k in after if k not in ('array_entries','array_bytes','max_bytes')}
    return {'status':actual['status'],'source_equivalence':checks,'wall_seconds':elapsed,'stage_seconds':stages,'stage_calls':counts,
            'other_orchestration_seconds':elapsed-sum(stages.values()),'peak_process_rss_bytes':peak_memory_bytes(),
            'cache_counters':cache,'array_bytes':after['array_bytes'],'observation_count':len(actual['observations']),
            'source_measurement_counts':[o['measurement_count'] for o in actual['observations']],
            'observation_errors':[o['error'] for o in actual['observations']],
            'evaluated_pair_count':len(actual['pair_evaluations']),'selected_pair':actual.get('confirmation',{}).get('selected_pair'),
            'scope':'Fresh source photometry using real cached science/catalog/calibration products and supplied records. Discovery not timed or executed. Calibration metadata HTTP remains possible if association cache misses; counters report actual lookups.'}


BASELINE=Path('artifacts/rw_empirical/report.json')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True,type=Path)
    parser.add_argument('--remeasure',action='store_true',help='Also recompute 2/5/8 observation native photometry, twice per count')
    args=parser.parse_args();out=args.output
    if out.exists():raise ValueError('Use a fresh validation output directory')
    out.mkdir(parents=True)
    baseline=read(BASELINE);catalog=read(BASELINE.parent/'catalog.json')
    packets={o['publisher_id']:read(o['measurement_file']) for o in baseline['observations']}
    records=[p['record'] for p in packets.values()]
    result={'scope':'Real cached RW Aurigae inputs. Regression and replay reuse explicitly recorded measurements, inspect FITS and calibration locally, and never claim fresh photometry speed from replay timings.',
            'baseline_report':str(BASELINE),'baseline_report_sha256':hashlib.sha256(BASELINE.read_bytes()).hexdigest(),
            'peak_memory_policy':'Process peak RSS including imports and native arrays; cumulative high-water mark within this validation process, not isolated per-run peaks.',
            'replay':{},'remeasure':{}}
    def checkpoint():write(out/'validation.json',result)
    actual=engine.run_selected(records,output=out/'legacy_regression',ra=baseline['configuration']['ra'],dec=baseline['configuration']['dec'],catalog=catalog,premeasured_packets=packets)
    result['legacy_regression']=regression(baseline,actual);checkpoint()
    for n in (2,5,8):
        started=time.perf_counter();path=out/f'replay_{n}'
        report=run_analysis(output=path,mode='all_valid_pairs',maximum=n,selection_mode='replay',replay_report=BASELINE)
        result['replay'][str(n)]=summary(report,path,time.perf_counter()-started);checkpoint()
    started=time.perf_counter();path=out/'time_series_5'
    report=run_analysis(output=path,mode='time_series',maximum=5,selection_mode='replay',replay_report=BASELINE)
    result['time_series_5']=summary(report,path,time.perf_counter()-started);checkpoint()
    if args.remeasure:
        for n in (2,5,8):
            selected=select_observations(records,maximum=n)
            for repeat in range(2):
                if repeat==0:DEFAULT_CACHE.clear()
                key=f'{n}_'+('cold_process_arrays' if repeat==0 else 'warm_process_arrays')
                try:result['remeasure'][key]=remeasure(selected,catalog,baseline,out/f'remeasure_{key}')
                except Exception as exc:result['remeasure'][key]={'status':'FAILED','error':{'type':type(exc).__name__,'reason':str(exc)}}
                checkpoint()
    result['verification_passed']=result['legacy_regression']['passed'] and all(r['artifact_integrity']['all_passed'] and r['evaluated_pair_count']==int(n)*(int(n)-1)//2 for n,r in result['replay'].items())
    if args.remeasure:
        result['verification_passed']=result['verification_passed'] and all(r.get('source_equivalence') and all(all(v.values()) for v in r['source_equivalence'].values()) and not any(r['observation_errors']) for r in result['remeasure'].values())
    checkpoint()
    print(json.dumps({'verification_passed':result['verification_passed'],'result':str(out/'validation.json'),
                      'replay_seconds':{n:r['wall_seconds_including_outputs'] for n,r in result['replay'].items()},
                      'remeasure_seconds':{n:r.get('wall_seconds') for n,r in result['remeasure'].items()}},indent=2))
    if not result['verification_passed']:raise SystemExit(2)


if __name__=='__main__':main()
