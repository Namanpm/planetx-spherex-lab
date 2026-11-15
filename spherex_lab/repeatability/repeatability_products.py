"""Export conditional source contrasts and confirmation diagnostics."""
import json
from pathlib import Path
import numpy as np
from astropy.io import fits

def export(report_path):
    path=Path(report_path);report=json.loads(path.read_text(encoding='utf-8'));out=path.parent
    selected=report.get('confirmation',{}).get('selected_pair',report.get('diagnostic_pair'))
    if not selected:return None
    packets={r['publisher_id']:json.loads(Path(r['measurement_file']).read_text(encoding='utf-8')) for r in report['observations']}
    a,b=[packets[i]['measurements'] for i in selected]
    model=report.get('confirmation',{}).get('model')
    chosen=next((r for r in report['pair_evaluations'] if [r['observation_a'],r['observation_b']]==selected),{})
    if model is None:model=chosen.get('validation_model',chosen.get('exploratory_model_not_confirmation_eligible'))
    if model is None:return None
    ids=model['heldout']['source_ids'];cov=np.asarray(model['heldout']['predictive_covariance'])
    prediction=np.array([model['coefficients'][0]+(model['coefficients'][1]*(a[i]['color']-model['reference_color']) if len(model['coefficients'])>1 else 0) for i in ids])
    residual=np.asarray(model['heldout']['residual_log_ratio'])
    table=fits.BinTableHDU.from_columns([fits.Column(name='SOURCE_ID',format='24A',array=np.array(ids,dtype='S24')),fits.Column(name='PREDICTED_FLUX_B',format='D',unit='Jy',array=np.array([a[i]['flux_jy'] for i in ids])*np.exp(prediction)),fits.Column(name='OBSERVED_FLUX_B',format='D',unit='Jy',array=np.array([b[i]['flux_jy'] for i in ids])),fits.Column(name='RESIDUAL_LN_RATIO',format='D',array=residual),fits.Column(name='SIGMA_LN_RATIO',format='D',array=np.sqrt(np.diag(cov)))],name='HELDOUT_CONTRASTS')
    table.header['STATUS']=report['status'];table.header['OBS_A']=selected[0];table.header['OBS_B']=selected[1]
    table.header['COMMENT']='Conditional reference-population contrast; not image subtraction or variability discovery.'
    fits.HDUList([fits.PrimaryHDU(),table,fits.ImageHDU(cov,name='LOG_RATIO_COVARIANCE')]).writeto(out/'heldout_contrasts.fits',overwrite=True,checksum=True)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(12,4.5))
    training=model['retained_training_ids'];validation=chosen.get('validation_model',{}).get('heldout',{}).get('source_ids',[])
    for group,label,color in [(training,'Training','0.65'),(validation,'Validation','tab:blue'),(ids,'Confirmation','tab:orange')]:
        c=np.array([a[i]['color'] for i in group]);raw=np.array([np.log(b[i]['flux_jy']/a[i]['flux_jy']) for i in group])
        axes[0].scatter(c,raw,label=f'{label} (n={len(group)})',color=color,s=20)
    grid=np.linspace(min(a[i]['color'] for i in training),max(a[i]['color'] for i in training),100)
    line=np.full_like(grid,model['coefficients'][0])
    if len(model['coefficients'])>1:line+=model['coefficients'][1]*(grid-model['reference_color'])
    axes[0].plot(grid,line,color='black',label='Training-only mean model')
    axes[0].set(xlabel='Gaia BP-RP (mag)',ylabel='Raw ln(flux B / flux A)');axes[0].legend(fontsize=8)
    axes[1].axhspan(-1.959963984540054,1.959963984540054,color='0.92',label='Nominal conditional 95% range')
    axes[1].axhline(0,color='black',lw=.8)
    axes[1].errorbar([a[i]['color'] for i in ids],model['heldout']['standardized_residual'],yerr=1,fmt='o',color='tab:orange')
    axes[1].set(xlabel='Gaia BP-RP (mag)',ylabel='Confirmation residual / predicted sigma');axes[1].legend(fontsize=8)
    fig.suptitle('Restricted reference-population comparison: no target variability inference')
    fig.tight_layout();fig.savefig(out/'confirmation_diagnostic.png',dpi=150);plt.close(fig)
    return str(out/'heldout_contrasts.fits')

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('report',nargs='?');p.add_argument('--target');p.add_argument('--output',default='artifacts/reference_experiment');p.add_argument('--max-observations',type=int,default=8);args=p.parse_args()
    if args.target:
        from spherex_lab.repeatability.empirical_comparison import run
        run(args.target,args.output,args.max_observations);args.report=str(Path(args.output)/'report.json')
    if not args.report:p.error('provide report path or --target')
    print(export(args.report))
