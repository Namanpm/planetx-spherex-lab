"""Reproducible Monte Carlo validation of conditional inference, never sky detections."""
import json
from pathlib import Path
import time
import tracemalloc
import numpy as np
from spectral_model import fit_spectrum, predict_spectrum
from motion_inference import fit_motion


def validate(output='artifacts/inference_validation', repeats=200, seed=7821):
    if not 10<=repeats<=1000: raise ValueError('repeats must be 10..1000')
    started=time.perf_counter(); tracemalloc.start(); rng=np.random.default_rng(seed)
    design=np.column_stack([np.ones(8),np.linspace(-1,1,8)]); beta=np.array([2.,.5]); noise=.1
    times=60000+np.arange(6)*180; dt=(times-times.mean())/365.25; truth=np.array([.2,-.1]); sigma=.03
    spectral_errors=[];spectral_cover=[];motion_errors=[];motion_cover=[];heldout=[]
    for _ in range(repeats):
        observed=design@beta+rng.normal(0,noise,8)
        fit=fit_spectrum(design[:6],observed[:6],np.full(6,noise**2))
        spectral_errors.append(fit['coefficients']-beta)
        spectral_cover.append(np.abs(fit['coefficients']-beta)<=1.959963984540054*np.sqrt(np.diag(fit['covariance'])))
        heldout.append(predict_spectrum(fit,design[6:],observed=observed[6:],measurement_variance=np.full(2,noise**2))['chi2'])
        x=1+truth[0]*dt+rng.normal(0,sigma,6);y=-1+truth[1]*dt+rng.normal(0,sigma,6)
        motion=fit_motion(times,x,y,np.full(6,sigma),np.full(6,sigma))
        err=np.array(motion['velocity_arcsec_per_year'])-truth;motion_errors.append(err)
        errors=np.sqrt([motion['covariance_x'][1][1],motion['covariance_y'][1][1]])
        motion_cover.append(np.abs(err)<=1.959963984540054*errors)
    current,peak=tracemalloc.get_traced_memory();tracemalloc.stop()
    report={'status':'SYNTHETIC','seed':seed,'repeats':repeats,
            'spectral_bias':np.mean(spectral_errors,axis=0).tolist(), 'spectral_95pct_interval_coverage':np.mean(spectral_cover,axis=0).tolist(),
            'mean_heldout_chi2_two_measurements':float(np.mean(heldout)),
            'motion_bias_arcsec_per_year':np.mean(motion_errors,axis=0).tolist(), 'motion_95pct_interval_coverage':np.mean(motion_cover,axis=0).tolist(),
            'elapsed_seconds':time.perf_counter()-started,'peak_traced_python_bytes':peak,
            'distribution':'Independent Gaussian noise; correct linear spectral basis, six true-associated motion epochs. Training first six spectral rows, held-out last two.',
            'limitations':['Conditional model validation only; no real response, PSF, association or calibration uncertainty.','Traced Python allocations are not peak process RAM.']}
    out=Path(output);out.mkdir(parents=True,exist_ok=True);(out/'report.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    from covariance_processing import compare_with_covariance
    from comparison_engine import write_outputs
    shape=(8,8); y,x=np.indices(shape); a=np.zeros(shape); b=2+.1*x+.05*y
    covariance=np.eye(a.size)+np.ones((a.size,a.size))*.1; valid=np.ones(shape,dtype=bool)
    diagnostic=compare_with_covariance(a,b,covariance,covariance,valid,valid,blank_mask=valid,synthetic=True)
    write_outputs(diagnostic,out/'covariance_example',provenance={'status':'SYNTHETIC','model':'known plane and supplied correlated covariance'})
    return report
if __name__=='__main__': print(json.dumps(validate(),indent=2))
