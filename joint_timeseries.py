"""Conditional joint epoch calibration; never chains pair corrections.

The anchor fixes the gauge, while every training source gets a nuisance log flux.
All epoch intercept/color coefficients are fitted at once. A source of interest is
excluded from every reference partition and from overlapping detector stamps.
"""
import hashlib
import numpy as np
from scipy.linalg import block_diag
from scipy.stats import chi2
from astropy.time import Time
from spectral_model import fit_spectrum

MAX_OBSERVATIONS = 8
MAX_REFERENCES = 300
MIN_TRAINING = 12
MIN_HELDOUT = 6
ALPHA = .05


def _partition(source):
    return int(hashlib.sha256(str(source).encode()).hexdigest(), 16) % 4


def _usable(measurement):
    try:
        flux = float(measurement['flux_jy']); variance = float(measurement['variance_jy2'])
        color = float(measurement['color'])
        return (np.isfinite([flux, variance, color]).all() and flux > 0 and variance > 0
                and np.sqrt(variance)/flux <= .1
                and measurement.get('central_valid', True)
                and measurement.get('central_flag', 0) in (0, 1 << 21)
                and measurement.get('association_status', 'CATALOG_PROPAGATED')
                not in ('AMBIGUOUS', 'REJECTED', 'FAILED'))
    except (KeyError, TypeError, ValueError):
        return False


def _overlap(a, b):
    # Stamp bounds are inclusive pixel centers, consistent with the native fits.
    return not (a[1] < b[0] or b[1] < a[0] or a[3] < b[2] or b[3] < a[2])


def _design(epoch, color, n_epochs):
    row = np.zeros(2*(n_epochs-1))
    if epoch:
        row[2*(epoch-1):2*epoch] = (1., color)
    return row


def _test_group(ids, measurements, colors, n_epochs, coefficients, covariance):
    residuals=[]; designs=[]; noise=[]
    for source in ids:
        rows=measurements[source]
        flux=np.array([m['flux_jy'] for m in rows],float)
        variance=np.array([m['variance_jy2'] for m in rows],float)/flux**2
        design=np.array([_design(e, colors[source], n_epochs) for e in range(1,n_epochs)])
        residuals.extend(np.log(flux[1:]/flux[0])-design@coefficients)
        designs.extend(design)
        noise.append(np.diag(variance[1:])+variance[0]*np.ones((n_epochs-1,n_epochs-1)))
    h=np.asarray(designs); residual=np.asarray(residuals)
    predictive=block_diag(*noise)+h@covariance@h.T
    statistic=float(residual@np.linalg.solve(predictive,residual))
    return {'source_ids':ids,'source_count':len(ids),'chi2':statistic,'dof':len(residual),
            'pvalue_chi2_approximate':float(chi2.sf(statistic,len(residual))),
            'residual_log_contrasts':residual.tolist(),
            'predictive_covariance':predictive.tolist(),
            'coverage_95_fraction':float(np.mean(np.abs(residual)<=1.959963984540054*np.sqrt(np.diag(predictive)))),
            'covariance_policy':'Shared anchor measurement within each source and shared joint calibration across every held-out source/epoch included.'}


def analyze_time_series(packets, source_id=None):
    """Fit a simultaneous joint scale and assess an independent source timeline.

    Input packets use existing empirical measurement fields. An accepted timeline
    additionally requires an explicit common_scale_domain in every packet with
    identical id, validated=True and provenance. This is an upstream applicability
    contract, not a validation invented from pair acceptance or representative
    response similarity. No real-data applicability evidence is assumed here.

    Log errors use the high-S/N delta approximation (fractional sigma <= .1).
    Native statistical errors are assumed independent, except for all derived
    shared-reference/anchor correlations calculated explicitly here. Unknown
    response, PSF and astrometric systematics are not manufactured or set to zero.
    """
    result={'status':'INSUFFICIENT_EVIDENCE','reason':None,'timelines':[],
            'source_id':str(source_id) if source_id is not None else None,
            'limitations':['Conditional reference-source model; a consistency pass does not prove source constancy or absolute calibration.',
                'Native statistical errors assumed independent; log errors use a high-S/N delta approximation.',
                'No unknown response, ePSF, astrometric or photometric systematic covariance is invented.',
                'No moving-object identification or RW Aurigae variability claim is authorized by this reference model.']}
    if not isinstance(packets,dict) or not 2<=len(packets)<=MAX_OBSERVATIONS:
        result['reason']='Require 2..8 explicitly selected observation packets.';return result
    ready=[]; failures=[]
    for ident,packet in packets.items():
        try:
            if packet.get('error') or not packet.get('measurements'):raise ValueError('failed observation or no usable measurements')
            date=packet.get('date') or next(iter(packet['measurements'].values())).get('date')
            mjd=float(Time(date).mjd)
            if not np.isfinite(mjd):raise ValueError('invalid timestamp')
            ready.append((mjd,str(ident),packet,date))
        except (TypeError,ValueError,KeyError) as exc:
            failures.append({'observation_id':str(ident),'reason':str(exc)})
    ready.sort(key=lambda r:(r[0],r[1]));result['rejected_observations']=failures
    result['observation_ids']=[r[1] for r in ready]
    if len(ready)<3 or len(set(r[0] for r in ready))<3:
        result['reason']='At least three distinct usable observation epochs are required.';return result
    n=len(ready); measurements={};source_id=str(source_id) if source_id is not None else None
    common=set.intersection(*(set(r[2]['measurements']) for r in ready))
    excluded=[];chosen=[]
    for source in sorted(common):
        if source==source_id:continue
        rows=[r[2]['measurements'][source] for r in ready]
        if not all(_usable(m) for m in rows):continue
        if not all('fit_bounds' in m for m in rows):
            excluded.append({'source_id':source,'reason':'missing native stamp bounds'});continue
        if not np.allclose([m['color'] for m in rows],rows[0]['color'],rtol=0,atol=1e-12):
            excluded.append({'source_id':source,'reason':'inconsistent catalog color identity'});continue
        overlaps=False
        for e,r in enumerate(ready):
            interest=r[2]['measurements'].get(source_id) if source_id else None
            if interest and ('fit_bounds' not in interest or _overlap(rows[e]['fit_bounds'],interest['fit_bounds'])):overlaps=True
            if any(_overlap(rows[e]['fit_bounds'],measurements[old][e]['fit_bounds']) for old in chosen):overlaps=True
        if overlaps:
            excluded.append({'source_id':source,'reason':'overlapping detector stamps across source partitions or source of interest'});continue
        chosen.append(source);measurements[source]=rows
        if len(chosen)==MAX_REFERENCES:break
    train=[s for s in chosen if _partition(s)<2]
    validation=[s for s in chosen if _partition(s)==2]
    confirmation=[s for s in chosen if _partition(s)==3]
    result['partition_source_ids']={'training':train,'validation':validation,'confirmation':confirmation}
    result['excluded_references']=excluded
    if len(train)<MIN_TRAINING or len(validation)<MIN_HELDOUT or len(confirmation)<MIN_HELDOUT:
        result['reason']='Insufficient globally stamp-disjoint training/validation/confirmation sources.';return result
    # Bound full heldout covariance allocations at <= 512 contrasts per group.
    if max(len(validation),len(confirmation))*(n-1)>512:
        result['reason']='Held-out full covariance exceeds 512-contrast budget.';return result
    reference=float(np.median([measurements[s][0]['color'] for s in train]))
    colors={s:float(measurements[s][0]['color']-reference) for s in chosen}
    lo=min(colors[s] for s in train);hi=max(colors[s] for s in train)
    if any(not lo<=colors[s]<=hi for s in validation+confirmation):
        result['reason']='Independent references require color extrapolation beyond training domain.';return result
    design=[];log_flux=[];variances=[]
    for j,source in enumerate(train):
        for e,m in enumerate(measurements[source]):
            row=np.zeros(len(train)+2*(n-1));row[j]=1.
            row[len(train):]=_design(e,colors[source],n)
            design.append(row);log_flux.append(np.log(m['flux_jy']));variances.append(m['variance_jy2']/m['flux_jy']**2)
    try:
        fit=fit_spectrum(design,log_flux,variances)
        coefficients=fit['coefficients'][len(train):]
        covariance=fit['covariance'][len(train):,len(train):]
        result['joint_calibration']={'model':'Simultaneous source log-flux nuisance parameters plus epoch intercept/catalog-color offsets; first epoch fixed as anchor.',
            'anchor_observation_id':ready[0][1],'reference_color':reference,
            'coefficients':coefficients.tolist(),'coefficient_covariance':covariance.tolist(),
            'training_chi2':fit['chi2'],'training_dof':fit['dof'],'training_source_count':len(train),
            'condition':fit['condition'],'rank_identified':True,'training_outlier_rejection':'none; no heldout-based tuning',
            'color_domain':[lo+reference,hi+reference]}
        validation_result=_test_group(validation,measurements,colors,n,coefficients,covariance)
        result['validation']=validation_result
        if validation_result['pvalue_chi2_approximate']<ALPHA:
            result['reason']='Independent joint-scale validation rejected; confirmation remains untouched.';return result
        # One confirmation after fixed model/epoch selection; no alternate subset search.
        confirmation_result=_test_group(confirmation,measurements,colors,n,coefficients,covariance)
        result['confirmation']=confirmation_result
        result['confirmation_tested_once']=True
        if confirmation_result['pvalue_chi2_approximate']<ALPHA:
            result['reason']='Independent joint-scale confirmation rejected; no alternate subset tried.';return result
        domains=[r[2].get('common_scale_domain',{}) for r in ready]
        if not all(d.get('validated') is True and d.get('id') and d.get('provenance') and d['id']==domains[0].get('id') for d in domains):
            result['reason']='Joint reference diagnostics pass, but no independently validated common response/source applicability domain was supplied.';return result
        result['common_scale_domain']=domains[0]
        if source_id is None:
            result['reason']='Joint reference scale validated conditionally; supply an independent source of interest for a timeline.';return result
        points=[];rows=[];epochs=[];missing=[dict(f, date=packets[f['observation_id']].get('date')) for f in failures]
        for e,(mjd,ident,packet,date) in enumerate(ready):
            m=packet['measurements'].get(source_id)
            if m is None or not _usable(m) or not lo<=m['color']-reference<=hi:
                missing.append({'observation_id':ident,'date':date,'reason':'missing/invalid measurement or source outside training color domain'});continue
            row=_design(e,m['color']-reference,n);correction=float(row@coefficients)
            calibrated=float(m['flux_jy']*np.exp(-correction))
            points.append({'observation_id':ident,'date':date,'mjd':mjd,'source_id':source_id,
                'measured_flux_jy':m['flux_jy'],'calibrated_flux_jy':calibrated,
                'statistical_variance_jy2':float(m['variance_jy2']*np.exp(-2*correction)),
                'calibration_variance_jy2':float(calibrated**2*(row@covariance@row)),
                'quality_flags':m.get('flag_values',[]),'association_status':m.get('association_status','CATALOG_PROPAGATED')})
            rows.append(row);epochs.append(e)
        if len(points)<3 or len(set(p['mjd'] for p in points))<3:
            result['reason']='Source of interest has fewer than three valid distinct epochs.'
            result['missing_epochs']=missing;return result
        h=np.asarray(rows);values=np.array([p['calibrated_flux_jy'] for p in points])
        shared=(values[:,None]*values[None,:])*(h@covariance@h.T)
        full=np.diag([p['statistical_variance_jy2'] for p in points])+shared
        one=np.ones(len(points));inv_values=np.linalg.solve(full,values);inv_one=np.linalg.solve(full,one)
        mean=float(one@inv_values/(one@inv_one));residual=values-mean
        statistic=float(residual@np.linalg.solve(full,residual))
        timeline={'source_id':source_id,'points':points,'missing_epochs':missing,
            'flux_covariance_jy2':full.tolist(),'shared_calibration_covariance_jy2':shared.tolist(),
            'constant_flux_assessment':{'weighted_mean_jy':mean,'chi2':statistic,'dof':len(points)-1,
                'pvalue_chi2_approximate':float(chi2.sf(statistic,len(points)-1)),
                'interpretation':'Conditional consistency diagnostic; small p-value requires follow-up, not an automatic astrophysical detection.'}}
        result.update(status='SUCCESS',reason='Independent joint reference validation and confirmation passed within supplied common-scale applicability domain.',timelines=[timeline])
    except (ValueError,np.linalg.LinAlgError,KeyError) as exc:
        result['reason']='Joint model is not identifiable or covariance is invalid: '+str(exc)
    return result
