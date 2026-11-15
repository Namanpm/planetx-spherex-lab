"""Identifiable linear spectral inference for caller-validated response curves.

The model is conditional on the supplied basis and calibrated response samples.
It does not reconstruct spectral information absent from the measurements or
establish that an externally supplied response applies to a SPHEREx observation.
"""
import numpy as np


def response_design(wavelength_um, responses, basis):
    """Integrate R_i(lambda) basis_j(lambda) using trapezoidal quadrature.

    No response normalization is silently performed. The caller must supply
    physically consistent responses and a basis on the same wavelength grid.
    """
    w=np.asarray(wavelength_um,dtype=float)
    r=np.asarray(responses,dtype=float)
    b=np.asarray(basis,dtype=float)
    if w.ndim!=1 or w.size<2 or not np.isfinite(w).all() or not (np.diff(w)>0).all():
        raise ValueError('wavelength grid must be finite and increasing')
    if r.ndim!=2 or r.shape[1]!=w.size or b.ndim!=2 or b.shape[0]!=w.size:
        raise ValueError('response and basis grids are inconsistent')
    if not np.isfinite(r).all() or not np.isfinite(b).all() or np.any(r<0) or not np.all(np.sum(r,axis=1)>0):
        raise ValueError('invalid response or basis values')
    return np.trapezoid(r[:,:,None]*b[None,:,:],x=w,axis=1)


def fit_spectrum(design, measurements, variance, *, ridge=0., max_condition=1e12):
    """Weighted least squares with optional zero-centered quadratic penalty.

    Rank/conditioning are checked before regularization. Regularization does not
    turn unidentified data into an identified measurement. Covariance is the
    frequentist conditional covariance of this estimator under independent
    Gaussian measurement errors; it excludes spectral-model/response error and
    regularization bias. The condition limit is a numerical policy, not a
    scientific calibration threshold.
    """
    x=np.asarray(design,dtype=float); y=np.asarray(measurements,dtype=float); v=np.asarray(variance,dtype=float)
    if x.ndim!=2 or x.shape[1]<1 or y.ndim!=1 or v.shape!=y.shape or x.shape[0]!=y.size:
        raise ValueError('inconsistent design, measurement or variance shapes')
    if not np.isfinite(x).all() or not np.isfinite(y).all() or not np.isfinite(v).all() or np.any(v<=0):
        raise ValueError('finite measurements and strictly positive variance required')
    if not np.isfinite(ridge) or ridge<0 or not np.isfinite(max_condition) or max_condition<1:
        raise ValueError('invalid numerical fit policy')
    weighted=x/np.sqrt(v[:,None]); target=y/np.sqrt(v)
    if not np.isfinite(weighted).all() or not np.isfinite(target).all():
        raise ValueError('whitening exceeds finite numerical range')
    u,values,vt=np.linalg.svd(weighted,full_matrices=False)
    tolerance=(values[0]*max(weighted.shape)*np.finfo(float).eps) if values.size else 0.
    rank=int(np.count_nonzero(values>tolerance))
    if rank<x.shape[1]: raise ValueError('spectral model is not identifiable from data')
    condition=float(values[0]/values[-1])
    if condition>max_condition: raise ValueError('spectral design is numerically ill conditioned')
    # SVD avoids squaring the design condition number in normal equations.
    # These filters give the frequentist ridge-estimator covariance, not a
    # Bayesian posterior covariance or an estimate of regularization bias.
    inverse_values=1/values
    filters=inverse_values/(1+ridge*inverse_values**2) if ridge else inverse_values
    coefficients=vt.T@(filters*(u.T@target))
    estimator=vt.T*filters
    covariance=estimator@estimator.T
    if not np.isfinite(coefficients).all() or not np.isfinite(covariance).all():
        raise ValueError('spectral estimator exceeds finite numerical range')
    residual=y-x@coefficients
    return {'coefficients':coefficients,'covariance':covariance,'chi2':float(np.sum(residual**2/v)),
            'dof':int(y.size-x.shape[1]),'rank':rank,'condition':condition,'ridge':float(ridge),
            'residual':residual,'measurement_count':int(y.size),
            'limitations':['Conditional on supplied calibrated response and spectral basis.',
            'Independent Gaussian measurement errors assumed.',
            'Covariance excludes response error, model error and ridge bias.']}


def predict_spectrum(fit, design, *, observed=None, measurement_variance=None):
    """Held-out prediction; optional residual chi2 includes predictive variance.

    Observed test data must not have been used in fit. This function cannot
    establish sample independence; the caller records that provenance.
    """
    x=np.asarray(design,dtype=float); beta=np.asarray(fit['coefficients']); cov=np.asarray(fit['covariance'])
    if x.ndim!=2 or x.shape[1]!=beta.size or not np.isfinite(x).all(): raise ValueError('invalid prediction design')
    mean=x@beta; variance=np.einsum('ij,jk,ik->i',x,cov,x)
    result={'mean':mean,'variance':variance}
    if (observed is None)!=(measurement_variance is None): raise ValueError('observed values require measurement variance')
    if observed is not None:
        y=np.asarray(observed,dtype=float); v=np.asarray(measurement_variance,dtype=float)
        if y.shape!=mean.shape or v.shape!=mean.shape or not np.isfinite(y).all() or not np.isfinite(v).all() or np.any(v<=0):
            raise ValueError('invalid held-out observations')
        residual=y-mean; total_cov=x@cov@x.T+np.diag(v)
        result.update(residual=residual,total_covariance=total_cov,chi2=float(residual@np.linalg.solve(total_cov,residual)))
    return result
