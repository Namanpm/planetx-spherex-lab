"""Scientific artifacts from explicit measurement packets and auditable decisions."""
import hashlib
import json
from pathlib import Path
import shutil

import numpy as np
from astropy.io import fits
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def _json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False), encoding='utf-8')


def _table(path, name, rows, fields):
    columns=[]
    for key, fmt, unit in fields:
        default=np.nan if fmt=='D' else ''
        values=[r.get(key,default) if r.get(key) is not None else default for r in rows]
        if fmt.endswith('A'):values=[str(v).encode('ascii','backslashreplace').decode('ascii') for v in values]
        columns.append(fits.Column(name=key.upper(), format=fmt, unit=unit, array=values))
    table=fits.BinTableHDU.from_columns(columns,name=name)
    table.header['COMMENT']='Measured data and explicit missing/rejected states; no interpolated observations.'
    fits.HDUList([fits.PrimaryHDU(),table]).writeto(path,overwrite=True,checksum=True)


def export(report, packets, output):
    """Return flat artifact filenames; manifest records SHA256 of their bytes.

    Only a SUCCESS time-series model can produce calibrated timeline points.
    Missing data remain NaN in FITS and explicit states in JSON/plots.
    """
    out=Path(output);out.mkdir(parents=True,exist_ok=True);names=[]
    def remember(name):names.append(name);return out/name
    states=report.get('readiness',[])
    _json(remember('observation_readiness.json'),states)
    readiness=[]
    for state in states:
        failure=state.get('failure') or {}
        readiness.append({'observation_id':state['observation_id'],'epoch':state.get('epoch'),'mjd':state.get('mjd'),
            'status':state['status'],'detector':str(state.get('detector','')),'measurement_count':float(state.get('valid_measurement_count',0)),
            'reason':failure.get('reason','')})
    _table(remember('observation_readiness.fits'),'READINESS',readiness,
           [('observation_id','128A',None),('epoch','40A',None),('mjd','D','d'),('status','24A',None),
            ('detector','8A',None),('measurement_count','D',None),('reason','1024A',None)])
    all_sources=sorted({str(s) for packet in packets.values() for s in packet.get('measurements',{})} |
                       {str(r['source_id']) for packet in packets.values() for r in packet.get('rejected',[]) if 'source_id' in r})
    rows=[]
    for ident,packet in packets.items():
        accepted=packet.get('measurements',{});rejected={str(r['source_id']):r.get('reason','') for r in packet.get('rejected',[]) if 'source_id' in r}
        failure=packet.get('error') or {}
        for source in all_sources:
            m=accepted.get(source)
            row={'observation_id':ident,'source_id':source,'epoch':packet.get('date'),
                 'mjd':packet.get('record',{}).get('t_min'),'status':'MEASURED' if m else ('REJECTED' if source in rejected or failure else 'MISSING'),
                 'reason':'' if m else rejected.get(source,failure.get('reason','source not measured in this epoch'))}
            if m:
                row.update({k:m.get(k) for k in ('ra_deg','dec_deg','flux_jy','variance_jy2')})
                row['flags']=json.dumps(m.get('flag_values',[]))
                row['association_status']=m.get('association_status','CATALOG_PROPAGATED')
            rows.append(row)
    _table(remember('source_measurements.fits'),'MEASUREMENTS',rows,
           [('observation_id','128A',None),('source_id','40A',None),('epoch','40A',None),('mjd','D','d'),
            ('status','16A',None),('ra_deg','D','deg'),('dec_deg','D','deg'),('flux_jy','D','Jy'),
            ('variance_jy2','D','Jy2'),('flags','256A',None),('association_status','40A',None),('reason','1024A',None)])
    pairs=report.get('pair_evaluations',[])
    _table(remember('pair_compatibility.fits'),'PAIR_DECISIONS',
           [{'observation_a':p['observation_a'],'observation_b':p['observation_b'],'status':p['status'],'reason':p.get('reason','')} for p in pairs],
           [('observation_a','128A',None),('observation_b','128A',None),('status','40A',None),('reason','1024A',None)])
    ids=report.get('selected_observation_ids',list(packets));index={ident:i for i,ident in enumerate(ids)}
    grid=np.full((len(ids),len(ids)),np.nan)
    for p in pairs:
        a,b=index[p['observation_a']],index[p['observation_b']]
        grid[a,b]=grid[b,a]=1. if p['status']=='VALIDATION_CONSISTENT' else 0.
    fig,ax=plt.subplots(figsize=(7,6));cmap=plt.get_cmap('RdYlGn').with_extremes(bad='#dddddd')
    ax.imshow(grid,vmin=0,vmax=1,cmap=cmap)
    for a in range(len(ids)):
        for b in range(len(ids)):
            if a!=b and np.isfinite(grid[a,b]):ax.text(b,a,'V' if grid[a,b] else 'R',ha='center',va='center')
    ax.set(xticks=range(len(ids)),yticks=range(len(ids)),xticklabels=range(1,len(ids)+1),yticklabels=range(1,len(ids)+1),
           xlabel='Selected observation index',ylabel='Selected observation index',title='V: validation consistent (not final confirmation)   R: rejected')
    fig.tight_layout();fig.savefig(remember('pair_compatibility.png'),dpi=140);plt.close(fig)
    fig,ax=plt.subplots(figsize=(11,4))
    for i,state in enumerate(states):
        mjd=state.get('mjd');color={'READY':'#18794e','PARTIAL':'#a65f00','FAILED':'#b42318','REJECTED':'#b42318'}.get(state['status'],'gray')
        if mjd is not None and np.isfinite(mjd):
            ax.scatter(mjd,i,c=color,marker='s',s=55);ax.annotate(state['status'],(mjd,i),xytext=(7,4),textcoords='offset points',fontsize=8)
    ax.set(yticks=range(len(states)),yticklabels=[s['observation_id'].split('?')[-1] for s in states],xlabel='Observation date (MJD days)',title='Selected real observations and readiness; missing epochs are not interpolated')
    fig.tight_layout();fig.savefig(remember('observation_timeline.png'),dpi=140);plt.close(fig)
    time_series=report.get('time_series',{'status':'NOT_REQUESTED'})
    _json(remember('time_series_status.json'),time_series)
    timelines=time_series.get('timelines',[]) if time_series.get('status')=='SUCCESS' else []
    fig,ax=plt.subplots(figsize=(9,4))
    if timelines:
        for timeline in timelines:
            points=timeline['points'];cov=np.asarray(timeline['flux_covariance_jy2'],float)
            if len(points)<3 or len({p['mjd'] for p in points})<3:
                raise ValueError('Accepted timeline requires three distinct measured epochs')
            if cov.shape!=(len(points),len(points)) or not np.isfinite(cov).all() or np.any(np.diag(cov)<0):
                raise ValueError('Invalid accepted timeline covariance')
            ax.errorbar([p['mjd'] for p in points],[p['calibrated_flux_jy'] for p in points],yerr=np.sqrt(np.diag(cov)),fmt='o',linestyle='none',label=timeline['source_id'])
        ax.legend(title='Logical source ID');ax.set(xlabel='Observation date (MJD days)',ylabel='Calibrated flux (Jy)',title='Conditional common-scale measurements; error bars: marginal 1 sigma')
        hdus=[fits.PrimaryHDU()]
        for j,timeline in enumerate(timelines):
            hdu=fits.ImageHDU(np.asarray(timeline['flux_covariance_jy2']),name=f'COV{j}');hdu.header['BUNIT']='Jy2';hdu.header['SOURCEID']=timeline['source_id'];hdus.append(hdu)
        fits.HDUList(hdus).writeto(remember('timeline_covariance.fits'),overwrite=True,checksum=True)
        timeline_rows=[{**p,'total_variance_jy2':float(t['flux_covariance_jy2'][j][j])} for t in timelines for j,p in enumerate(t['points'])]
        _table(remember('timeline_measurements.fits'),'TIMELINE',timeline_rows,[('observation_id','128A',None),('source_id','40A',None),('date','40A',None),('mjd','D','d'),('measured_flux_jy','D','Jy'),('calibrated_flux_jy','D','Jy'),('statistical_variance_jy2','D','Jy2'),('calibration_variance_jy2','D','Jy2'),('total_variance_jy2','D','Jy2')])
    else:
        ax.axis('off');ax.text(.5,.6,time_series.get('status','INSUFFICIENT_EVIDENCE'),ha='center',va='center',transform=ax.transAxes,fontsize=15)
        ax.text(.5,.35,time_series.get('reason','No validated common-scale timeline requested or available.'),ha='center',va='center',transform=ax.transAxes,wrap=True,fontsize=9)
    fig.tight_layout();fig.savefig(remember('light_curve.png'),dpi=140);plt.close(fig)
    for source,name in [('reference_measurements.fits','selected_pair_measurements.fits'),('reference_diagnostic.png','selected_pair_diagnostic.png'),('heldout_contrasts.fits','selected_pair_contrasts.fits'),('confirmation_diagnostic.png','selected_pair_confirmation.png')]:
        path=out/'empirical'/source
        if path.is_file():shutil.copyfile(path,remember(name))
    manifest={name:{'sha256':hashlib.sha256((out/name).read_bytes()).hexdigest(),'bytes':(out/name).stat().st_size} for name in names}
    _json(remember('artifact_manifest.json'),manifest)
    return names
