"""Profile an entire eight-observation workload without reusing measured flux caches."""
import argparse
import cProfile
import functools
import io
import json
from pathlib import Path
import pstats
import time
import numpy as np
from spherex_lab.repeatability import empirical_comparison as pipeline
from spherex_lab.imaging import native_photometry
from spherex_lab.imaging import psf_processing
from spherex_lab.imaging import response_calibration
from spherex_lab.imaging import reference_response
from spherex_lab.repeatability import empirical_repeatability
from spherex_lab.data import calibration_resolver
from spherex_lab.repeatability.benchmark_verification import peak_memory_bytes, compare_reports


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--mode',choices=['wall','cprofile'],default='wall');parser.add_argument('--output',required=True);parser.add_argument('--calibration-cache',type=Path);parser.add_argument('--frozen-observations',action='store_true',help='Replay the eight real archived baseline records; isolate processing from SIA gateway failures');args=parser.parse_args()
    out=Path(args.output)
    if out.exists():raise ValueError('Use a fresh profile directory to avoid reusing measured-flux caches')
    out.mkdir(parents=True)
    baseline=Path('artifacts/rw_empirical')
    reference=json.loads((baseline/'report.json').read_text(encoding='utf-8'))
    timings={}; calls={}
    def hook(module,name,label,redirect=None):
        original=getattr(module,name)
        @functools.wraps(original)
        def wrapper(*a,**kw):
            started=time.perf_counter()
            try:
                if redirect is not None:return redirect(original,a,kw)
                return original(*a,**kw)
            finally:
                timings[label]=timings.get(label,0.)+time.perf_counter()-started
                calls[label]=calls.get(label,0)+1
        setattr(module,name,wrapper)
    if args.frozen_observations:
        records=[json.loads(Path(o['measurement_file']).read_text())['record'] for o in reference['observations']]
        hook(pipeline,'discover_pairs','archive_discovery',lambda fn,a,k:{'observations':records})
    else:hook(pipeline,'discover_pairs','archive_discovery')
    hook(pipeline,'query_references','catalog_access',lambda fn,a,k:fn(*a,**{**k,'cache_dir':baseline/'catalog_cache'}))
    hook(pipeline,'download_cutout','science_cache_access',lambda fn,a,k:fn(*a[:3],baseline/'science_cache',**k))
    hook(pipeline,'wavelength_maps','calibration_resolution',lambda fn,a,k:fn(*a[:2],args.calibration_cache or baseline/'calibration_cache',**k))
    hook(calibration_resolver.CalibrationCache,'association','nested_calibration_association')
    hook(calibration_resolver.CalibrationCache,'arrays','nested_calibration_arrays')
    hook(response_calibration,'fetch_bundle','response_bundle_access')
    hook(reference_response,'field_response','response_curve_sampling')
    hook(pipeline,'source_position','catalog_position_propagation')
    hook(pipeline,'fit_source','native_photometry')
    hook(pipeline,'prepare_pixel_areas','pixel_area_batch_preparation')
    hook(native_photometry,'local_pixel_areas_sr','nested_batched_pixel_area_geometry')
    if hasattr(pipeline,'get_observation'):
        hook(pipeline,'get_observation','observation_preparation')
        hook(native_photometry,'read_cutout_hdus','nested_prepare_FITS_WCS')
        hook(native_photometry,'prepare_epsf','nested_prepare_ePSF_table')
        hook(native_photometry,'extract_prepared_epsf','nested_cached_ePSF_selection')
    hook(pipeline,'fit_repeatability','empirical_model_fitting')
    hook(reference_response,'signed_laboratory_distance','response_pair_ranking')
    hook(empirical_repeatability,'spatial_repeatability_diagnostic','spatial_residual_diagnostics')
    hook(pipeline,'render_sources','output_rendering')
    hook(pipeline,'dump','json_output')
    # Nested diagnostic timings: these are contained in native_photometry, not additive phases.
    hook(native_photometry,'read_cutout','nested_FITS_read_and_WCS')
    hook(native_photometry,'extract_epsf','nested_ePSF_extraction')
    hook(native_photometry,'fit_stamp','nested_numerical_source_fit')
    hook(native_photometry,'local_pixel_area_sr','nested_pixel_area_geometry')
    prof=cProfile.Profile() if args.mode=='cprofile' else None
    started=time.perf_counter()
    if prof:prof.enable()
    result=pipeline.run('rw-aurigae',out/'run',8)
    if prof:prof.disable()
    total=time.perf_counter()-started
    if prof:
        prof.dump_stats(str(out/'full_workflow.prof'))
        buffer=io.StringIO();pstats.Stats(prof,stream=buffer).strip_dirs().sort_stats('cumulative').print_stats(50)
        (out/'cumulative.txt').write_text(buffer.getvalue(),encoding='utf-8')
        buffer=io.StringIO();pstats.Stats(prof,stream=buffer).strip_dirs().sort_stats('tottime').print_stats(50)
        (out/'self_time.txt').write_text(buffer.getvalue(),encoding='utf-8')
    phases={k:{'seconds':v,'percent_total':100*v/total,'calls':calls[k]} for k,v in timings.items() if not k.startswith('nested_')}
    counted=sum(timings[k] for k in phases)
    phases['other_orchestration']={'seconds':total-counted,'percent_total':100*(total-counted)/total,'calls':1}
    nested={}
    for key,value in timings.items():
        if not key.startswith('nested_'):continue
        parent='calibration_resolution' if key.startswith('nested_calibration_') else ('observation_preparation' if key.startswith('nested_prepare_') else ('pixel_area_batch_preparation' if key=='nested_batched_pixel_area_geometry' else 'native_photometry'))
        nested[key]={'seconds':value,'contained_in':parent,'percent_of_parent':100*value/timings[parent],'calls':calls[key]}
    selected=result.get('confirmation',{}).get('selected_pair')
    equivalent=selected==reference['confirmation']['selected_pair'] and [r['measurement_count'] for r in result['observations']]==[r['measurement_count'] for r in reference['observations']]
    if equivalent:equivalent=bool(np.allclose(result['confirmation']['model']['coefficients'],reference['confirmation']['model']['coefficients'],rtol=1e-12,atol=1e-12))
    summary={'mode':args.mode,'total_wall_seconds':total,'phases':phases,'nested_native_diagnostics':nested,'native_fit_calls':calls.get('native_photometry',0),'prepared_observation_count':result.get('prepared_observation_count'),'prepared_epsf_zone_count':result.get('prepared_epsf_zone_count'),'observation_count':len(result['observations']),'pair_count':len(result['pair_evaluations']),'same_measurement_counts_selected_pair_and_model_as_original':equivalent,'status':result['status'],'cache_policy':'Existing raw science/calibration/catalog files reused; all source measurements recomputed in a new output directory. Discovery and DataLink metadata requests still performed.','timing_policy':'Nested times are contained within native photometry. cProfile mode adds instrumentation overhead. Imports precede timed run.'}
    summary['frozen_observations']=args.frozen_observations
    summary['peak_process_rss_bytes']=peak_memory_bytes()
    summary['calibration_cache']=calibration_resolver.DEFAULT_CACHE.cache_info()
    summary['pixel_area_cache']=result.get('pixel_area_cache')
    summary['numerical_verification']=compare_reports(reference,result)
    summary['cache_policy']='All source measurements recomputed; raw science and catalog cached. Calibration cache: '+str(args.calibration_cache or baseline/'calibration_cache')+'; see cache counters for network and array reuse.'
    (out/'profile_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in summary.items() if k!='numerical_verification'},indent=2))
    print('Scientific equivalence:',summary['numerical_verification']['passed'],'differences:',len(summary['numerical_verification']['differences']))
if __name__=='__main__':main()
