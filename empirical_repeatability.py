"""Restricted-regime empirical repeatability with source-disjoint validation.

This models log(F_B/F_A), using the high-S/N delta-method error approximation.
It is not a universal SPHEREx calibration or an assertion that references are
constant. Source screening and matching conditions must be established upstream.
"""
from statistics import NormalDist
import numpy as np
from spectral_model import fit_spectrum


def _profile(x,y,v,scatter_variance):
    fit=fit_spectrum(x,y,v+scatter_variance)
    residual=y-x@fit['coefficients']
    likelihood=.5*np.sum(np.log(v+scatter_variance)+residual**2/(v+scatter_variance))
    return float(likelihood),fit


def _scatter_fit(x,y,v):
    # Search all grid-bracketed minima; include zero-scatter boundary explicitly.
    baseline=_profile(x,y,v,0.)
    upper=max(float(np.var(y)),float(np.max(v)),1e-16)
    while True:
        left=_profile(x,y,v,upper)[0];right=_profile(x,y,v,upper*4)[0]
        if right>left:break
        upper*=4
        if not np.isfinite(upper):raise ValueError('scatter likelihood has no finite bracket')
    grid=np.r_[0.,np.geomspace(max(upper*1e-12,np.finfo(float).tiny),upper*4,96)]
    scores=np.array([_profile(x,y,v,s)[0] for s in grid])
    candidates=[(baseline[0],0.,baseline[1])]
    for i in range(1,len(grid)-1):
        if scores[i]<=scores[i-1] and scores[i]<=scores[i+1]:
            lo,hi=grid[i-1],grid[i+1]
            for _ in range(64):
                p=hi-(hi-lo)*.6180339887498949;q=lo+(hi-lo)*.6180339887498949
                if _profile(x,y,v,p)[0]<_profile(x,y,v,q)[0]:hi=q
                else:lo=p
            s=(lo+hi)/2;score,fit=_profile(x,y,v,s);candidates.append((score,s,fit))
    score,s,fit=min(candidates,key=lambda item:item[0])
    # Profile-likelihood upper limit is diagnostic; asymptotic, not exact at zero boundary.
    threshold=score+NormalDist().inv_cdf(.975)**2/2
    lo=s;hi=max(upper,s+float(np.max(v)))
    while _profile(x,y,v,hi)[0]<threshold:
        hi*=4
        if not np.isfinite(hi):raise ValueError('scatter upper limit is unbounded')
    for _ in range(64):
        middle=(lo+hi)/2
        if _profile(x,y,v,middle)[0]<threshold:lo=middle
        else:hi=middle
    return fit,float(s),float(hi)


def fit_repeatability(source_ids,flux_a,flux_b,variance_a,variance_b,*,training_ids,heldout_ids,outlier_sigma,color=None,min_training=None):
    """Fit conditional multiplicative/color correction and test untouched sources.

    Required outlier_sigma is an explicit analyst policy, not an instrument
    tolerance. Rejection uses only training residuals. Held-out measurements
    never alter coefficients, scatter, source selection, or fit policy.
    Delta-method log errors require sufficiently high S/N; measured fractional
    errors are reported rather than inventing a universal S/N validity cutoff.
    """
    ids=np.asarray(source_ids,dtype=str);a,b,va,vb=[np.asarray(v,dtype=float) for v in (flux_a,flux_b,variance_a,variance_b)]
    if ids.ndim!=1 or len(ids)>10000 or len(set(ids.tolist()))!=len(ids) or any(v.shape!=ids.shape for v in (a,b,va,vb)):
        raise ValueError('unique source IDs and matching vectors required; max 10000 sources')
    if not all(np.isfinite(v).all() for v in (a,b,va,vb)) or any((v<=0).any() for v in (a,b,va,vb)):
        raise ValueError('positive flux and strictly positive variance required')
    if not np.isfinite(outlier_sigma) or outlier_sigma<=0:raise ValueError('positive explicit outlier policy required')
    training_ids=list(map(str,training_ids));heldout_ids=list(map(str,heldout_ids))
    if len(set(training_ids))!=len(training_ids) or len(set(heldout_ids))!=len(heldout_ids):raise ValueError('duplicate split IDs')
    if set(training_ids)&set(heldout_ids):raise ValueError('training and heldout IDs must be disjoint')
    if not set(training_ids+heldout_ids)<=set(ids.tolist()):raise ValueError('split contains missing source IDs')
    if not heldout_ids:raise ValueError('independent heldout sources required')
    if len(heldout_ids)>1024:raise ValueError('heldout full-covariance budget is 1024 sources')
    train=np.isin(ids,training_ids);test=np.isin(ids,heldout_ids)
    reference=0.
    x=np.ones((len(ids),1))
    if color is not None:
        c=np.asarray(color,dtype=float)
        if c.shape!=ids.shape or not np.isfinite(c).all():raise ValueError('finite colors required')
        reference=float(np.median(c[train])) if train.any() else 0.
        x=np.column_stack([np.ones(len(ids)),c-reference])
    needed=x.shape[1]+2 if min_training is None else min_training
    if not isinstance(needed,int) or needed<x.shape[1]+2:raise ValueError('minimum training requires parameters plus two residual degrees of freedom')
    if train.sum()<needed:raise ValueError('insufficient training sources')
    y=np.log(b/a);v=va/a**2+vb/b**2
    retained=train.copy()
    for _ in range(20):
        fit=fit_spectrum(x[retained],y[retained],v[retained])
        residual=y-x@fit['coefficients']
        center=float(np.median(residual[retained]));mad=float(np.median(np.abs(residual[retained]-center)))/NormalDist().inv_cdf(.75)
        keep=retained & (np.abs(residual-center)<=outlier_sigma*np.sqrt(v+mad**2))
        if keep.sum()<needed:raise ValueError('robust screening leaves insufficient references')
        if np.array_equal(keep,retained):break
        retained=keep
    else:raise ValueError('training outlier screening did not converge')
    fit,s2,upper=_scatter_fit(x[retained],y[retained],v[retained])
    prediction=x[test]@fit['coefficients'];residual=y[test]-prediction
    predictive_covariance=x[test]@fit['covariance']@x[test].T+np.diag(v[test]+s2)
    z=residual/np.sqrt(np.diag(predictive_covariance))
    extrapolated=0
    if color is not None:
        extrapolated=int(np.count_nonzero((c[test]<c[retained].min()) | (c[test]>c[retained].max())))
    from scipy.stats import chi2
    heldout_chi2=float(residual@np.linalg.solve(predictive_covariance,residual))
    return {'status':'EMPIRICAL_DIAGNOSTIC','model':'log flux ratio intercept'+(' plus catalog color' if color is not None else ''),
        'coefficients':fit['coefficients'].tolist(),'coefficient_covariance':fit['covariance'].tolist(),
        'reference_color':reference,'scale_at_reference_color':float(np.exp(fit['coefficients'][0])),
        'excess_log_scatter':float(np.sqrt(s2)),'excess_log_scatter_upper_95_profile':float(np.sqrt(upper)),
        'retained_training_ids':ids[retained].tolist(),'rejected_training_ids':ids[train & ~retained].tolist(),
        'training_source_count':int(retained.sum()),'training_rms_log_ratio':float(np.sqrt(np.mean(fit['residual']**2))),'training_chi2':fit['chi2'],'training_dof':fit['dof'],'outlier_sigma_policy':float(outlier_sigma),
        'max_fractional_measurement_error':float(max(np.max(np.sqrt(va)/a),np.max(np.sqrt(vb)/b))),
        'heldout':{'source_ids':ids[test].tolist(),'source_count':int(test.sum()),'residual_log_ratio':residual.tolist(),'standardized_residual':z.tolist(),
            'predictive_covariance':predictive_covariance.tolist(),'chi2':heldout_chi2,'dof':int(test.sum()),'pvalue_chi2_approximate':float(chi2.sf(heldout_chi2,int(test.sum()))),
            'rms_log_ratio':float(np.sqrt(np.mean(residual**2))),'coverage_95_fraction':float(np.mean(np.abs(z)<=NormalDist().inv_cdf(.975))),
            'extrapolated_color_source_count':extrapolated},
        'acceptance':'not automatically assessed; evidence is restricted to screened sources and tested observing regime',
        'limitations':['References must be screened upstream; this function does not establish astrophysical constancy.',
            'Log-flux uncertainty uses high-S/N delta approximation; fractional errors are reported.',
            'Independent measurement noise assumed; common correction covariance is included in heldout predictions.',
            'Intrinsic scatter is conditional on robust training selection and does not include unidentified correlated systematics.',
            'Coefficient covariance and heldout chi-square condition on fitted scatter; scatter-estimation uncertainty is not marginalized.',
            'Profile upper limit is asymptotic and not an exact confidence bound at the zero-scatter boundary.',
            'No target-specific correction is justified outside reference spectral/color and observing-condition support.']}


def spatial_repeatability_diagnostic(xy,residual,seed=0,permutations=499):
    """Nearest-neighbor Moran-style residual diagnostic, never acceptance policy.

    xy must use a common tangent-plane or native-pixel coordinate system.
    All equidistant nearest neighbors are included and the graph is symmetrized.
    Residuals should be externally standardized when measurement errors differ.
    Permutation p-values condition on exchangeability; fitted nuisance models,
    heteroscedasticity and shared instrumental terms can violate that assumption.
    No universal instrument covariance/error bound is inferred from this test.
    """
    positions=np.asarray(xy,dtype=float);r=np.asarray(residual,dtype=float)
    if positions.ndim!=2 or positions.shape[1]!=2 or r.shape!=(len(positions),) or not 4<=len(r)<=1000:
        raise ValueError('spatial diagnostic requires 4..1000 matching 2D positions')
    if not np.isfinite(positions).all() or not np.isfinite(r).all():raise ValueError('finite spatial positions and residuals required')
    if not isinstance(permutations,int) or not 19<=permutations<=10000:raise ValueError('permutation budget must be 19..10000')
    distances=np.sum((positions[:,None,:]-positions[None,:,:])**2,axis=2)
    np.fill_diagonal(distances,np.inf)
    if np.any(distances==0):raise ValueError('duplicated source positions do not define independent spatial samples')
    nearest=distances.min(axis=1)
    # Floating-point tie handling only; not a scientific distance tolerance.
    graph=distances<=nearest[:,None]*(1+16*np.finfo(float).eps)
    graph |= graph.T
    left,right=np.nonzero(graph);n=len(r);weight_sum=len(left)
    centered=r-r.mean();denominator=float(centered@centered)
    result={'method':'symmetric nearest-neighbor Moran I','source_count':n,'directed_edge_count':weight_sum,
        'seed':int(seed),'permutations':permutations,'expected_exchangeable_I':-1/(n-1),
        'limitations':['A residual spatial diagnostic, not a calibrated covariance model or instrument error budget.',
        'Permutation probabilities require exchangeable residuals; shared fitted corrections and unequal errors can violate this.',
        'Source screening, graph choice and repeated inspection affect inferential interpretation; no automatic acceptance threshold.']}
    if denominator==0:
        return {**result,'status':'DEGENERATE_CONSTANT_RESIDUAL','moran_I':None,'pvalue_permutation_positive':None,'pvalue_permutation_two_sided':None}
    def statistic(values):return float(n/weight_sum*np.sum(values[left]*values[right])/denominator)
    observed=statistic(centered);rng=np.random.default_rng(seed);null=np.array([statistic(rng.permutation(centered)) for _ in range(permutations)])
    expectation=-1/(n-1)
    result.update(status='DIAGNOSTIC_ONLY',moran_I=observed,
        pvalue_permutation_positive=float((1+np.count_nonzero(null>=observed))/(permutations+1)),
        pvalue_permutation_two_sided=float((1+np.count_nonzero(np.abs(null-expectation)>=abs(observed-expectation)))/(permutations+1)),
        permutation_null_mean=float(null.mean()),permutation_null_std=float(null.std(ddof=1)))
    return result
