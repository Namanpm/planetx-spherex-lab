"""Small-grid exact linear covariance propagation; explicit 1024-pixel budget."""
import numpy as np


def propagate(operator, variance):
    """L diag(v) L.T for independent input errors, without discarding covariance."""
    l=np.asarray(operator,dtype=float); v=np.asarray(variance,dtype=float)
    if l.ndim!=2 or l.shape[1]!=v.size or max(l.shape)>4096 or l.shape[0]>1024:
        raise ValueError('operator geometry/resource limit')
    if not np.isfinite(l).all() or not np.isfinite(v).all() or (v<0).any():
        raise ValueError('invalid operator or variances')
    return (l*v.ravel())@l.T


def remove_background(values,covariance,design,blank):
    """GLS fit on predeclared source-free pixels; propagate shared fit covariance.

    Assumes supplied covariance includes all known noise and is positive definite.
    Unknown calibration systematics remain excluded; no automatic sky declaration.
    """
    y=np.asarray(values,dtype=float); c=np.asarray(covariance,dtype=float); x=np.asarray(design,dtype=float); blank=np.asarray(blank,dtype=bool)
    n=y.size
    if y.ndim!=1 or n>1024 or c.shape!=(n,n) or x.ndim!=2 or x.shape[0]!=n or blank.shape!=(n,): raise ValueError('invalid geometry/resource limit')
    if not all(np.isfinite(a).all() for a in (y,c,x)) or not np.allclose(c,c.T): raise ValueError('invalid covariance or values')
    try: np.linalg.cholesky(c)
    except np.linalg.LinAlgError as exc: raise ValueError('covariance must be positive definite') from exc
    idx=np.flatnonzero(blank)
    if idx.size<=x.shape[1]: raise ValueError('insufficient blank measurements')
    cb=c[np.ix_(idx,idx)]; xb=x[idx]
    precision_x=np.linalg.solve(cb,xb)
    normal=xb.T@precision_x
    if np.linalg.matrix_rank(normal)<x.shape[1]: raise ValueError('unidentifiable background')
    parameter_covariance=np.linalg.inv(normal)
    fit_operator=parameter_covariance@precision_x.T
    selection=np.eye(n)[idx]
    residual_operator=np.eye(n)-x@fit_operator@selection
    coefficients=fit_operator@y[idx]
    residual_covariance=residual_operator@c@residual_operator.T
    return {'residual':residual_operator@y,'covariance':residual_covariance,'coefficients':coefficients,'parameter_covariance':parameter_covariance,'operator':residual_operator}


def compare_with_covariance(image_a,image_b,covariance_a,covariance_b,valid,spectral_mask,*,blank_mask=None,synthetic=False):
    """Common-grid subtraction using supplied full within-image covariance.

    A/B cross-exposure covariance is assumed zero and must be justified by caller.
    Full covariance must include reprojection/convolution already performed.
    """
    from comparison_engine import compare_aligned
    a=np.asarray(image_a,float);b=np.asarray(image_b,float);n=a.size
    ca=np.asarray(covariance_a,float);cb=np.asarray(covariance_b,float)
    if n>1024 or a.shape!=b.shape or ca.shape!=(n,n) or cb.shape!=(n,n): raise ValueError('invalid covariance image geometry/resource limit')
    mask=np.asarray(valid,bool)&np.asarray(spectral_mask,bool)
    if mask.shape!=a.shape: raise ValueError('invalid mask')
    idx=np.flatnonzero(mask)
    if not len(idx): raise ValueError('no valid comparison pixels')
    if not np.isfinite(ca).all() or not np.isfinite(cb).all() or not np.allclose(ca,ca.T) or not np.allclose(cb,cb.T): raise ValueError('invalid image covariance')
    for c in (ca,cb):
        try: np.linalg.cholesky(c)
        except np.linalg.LinAlgError as exc: raise ValueError('image covariance must be positive definite') from exc
    total=ca+cb;corrected=b.copy();residual_covariance=total.copy()
    background=np.zeros_like(a)
    if blank_mask is not None:
        y,x=np.indices(a.shape);design=np.column_stack([np.ones(n),x.ravel(),y.ravel()])
        fit=remove_background((b-a).ravel()[idx],total[np.ix_(idx,idx)],design[idx],np.asarray(blank_mask,bool).ravel()[idx])
        background=(design@fit['coefficients']).reshape(a.shape);corrected=b-background
        residual_covariance[np.ix_(idx,idx)]=fit['covariance']
    result=compare_aligned(a,np.diag(ca).reshape(a.shape),corrected,np.diag(cb).reshape(a.shape),valid,valid,spectral_mask,synthetic=synthetic,resampled=True)
    result['variance']=np.where(mask,np.diag(residual_covariance).reshape(a.shape),np.nan)
    result['uncertainty']=np.sqrt(np.maximum(result['variance'],0))
    result['significance']=np.divide(result['difference'],result['uncertainty'],out=np.full_like(a,np.nan),where=mask & (result['uncertainty']>0))
    result['background']=background;result['full_covariance']=residual_covariance
    result['report']['candidates']=[]
    result['report']['background_model']='full-covariance GLS plane' if blank_mask is not None else 'none'
    result['report']['limitations'].append('Full supplied within-exposure covariance propagated; cross-exposure covariance and unknown calibration systematics excluded. Candidate extraction suppressed after covariance correction.')
    return result
