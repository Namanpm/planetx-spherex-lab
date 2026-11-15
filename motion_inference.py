"""Conditional motion/variability fits to already-associated measurements.

No source association, sky projection, catalog identity, or discovery inference
is implied by fitting these supplied data. Coordinates must be in one common
small-angle tangent plane, not raw RA/Dec degrees.
"""
import numpy as np


def _fit(design, values, sigma):
    z=design/sigma[:,None]; target=values/sigma
    if not np.isfinite(z).all() or not np.isfinite(target).all():
        raise ValueError('motion whitening exceeds finite numerical range')
    u,s,vt=np.linalg.svd(z,full_matrices=False)
    tolerance=s[0]*max(z.shape)*np.finfo(float).eps
    if np.count_nonzero(s>tolerance)<design.shape[1]:
        raise ValueError('motion fit requires distinct constrained times')
    inverse=1/s
    beta=vt.T@(inverse*(u.T@target))
    estimator=vt.T*inverse
    covariance=estimator@estimator.T
    return beta,covariance,float(np.sum(((values-design@beta)/sigma)**2))



def fit_motion(times_mjd, x_arcsec, y_arcsec, sigma_x_arcsec, sigma_y_arcsec, *, reference_mjd=None):
    """Compare linear motion and static fits under independent Gaussian errors.

    At least three epochs are needed, so the linear model retains residual
    degrees of freedom. Reported delta chi2 has two additional parameters, but
    no classification threshold is imposed. The fit excludes observer parallax,
    acceleration, association uncertainty, and common astrometric systematics.
    """
    arrays=[np.asarray(a,dtype=float) for a in (times_mjd,x_arcsec,y_arcsec,sigma_x_arcsec,sigma_y_arcsec)]
    times,x,y,sx,sy=arrays
    if times.ndim!=1 or times.size<3 or any(a.shape!=times.shape or not np.isfinite(a).all() for a in arrays):
        raise ValueError('at least three finite, matching epoch measurements required')
    if np.unique(times).size<3: raise ValueError('at least three distinct observation epochs required')
    if np.any(sx<=0) or np.any(sy<=0): raise ValueError('position uncertainty must be positive')
    ref=float(np.mean(times) if reference_mjd is None else reference_mjd)
    if not np.isfinite(ref): raise ValueError('invalid reference epoch')
    dt=(times-ref)/365.25; design=np.column_stack([np.ones(times.size),dt])
    bx,cx,chi_x=_fit(design,x,sx); by,cy,chi_y=_fit(design,y,sy)
    _,_,static_x=_fit(design[:,:1],x,sx); _,_,static_y=_fit(design[:,:1],y,sy)
    return {'reference_mjd':ref,'position_arcsec':[float(bx[0]),float(by[0])],
            'velocity_arcsec_per_year':[float(bx[1]),float(by[1])],
            'covariance_x':cx.tolist(),'covariance_y':cy.tolist(),
            'chi2_motion':chi_x+chi_y,'chi2_static':static_x+static_y,
            'delta_chi2_static_minus_motion':static_x+static_y-chi_x-chi_y,
            'additional_parameters':2,'dof_motion':int(2*times.size-4),
            'limitations':['Already-associated common tangent-plane measurements required.',
            'Independent Gaussian position errors; no cross-coordinate covariance.',
            'No parallax, acceleration, catalog confirmation or discovery classification.']}


def fit_variability(flux, variance):
    """Constant-flux hypothesis diagnostic for compatible-band fluxes.

    The caller must establish spectral and photometric compatibility across
    epochs. This fit cannot distinguish variability from those systematics.
    """
    y=np.asarray(flux,dtype=float); v=np.asarray(variance,dtype=float)
    if y.ndim!=1 or y.size<2 or v.shape!=y.shape or not np.isfinite(y).all() or not np.isfinite(v).all() or np.any(v<=0):
        raise ValueError('at least two finite fluxes with positive variance required')
    weights=1/v; constant=float(np.sum(y*weights)/np.sum(weights))
    return {'constant_flux':constant,'variance_constant':float(1/np.sum(weights)),
            'chi2_constant':float(np.sum((y-constant)**2/v)),'dof':int(y.size-1),
            'limitations':['Requires already established band and photometric compatibility.',
            'Independent errors; calibration covariance and intrinsic source model absent.']}


def associate_epochs(positions_a, sigma_a, positions_b, sigma_b, *, max_separation_arcsec):
    """Associate only one-to-one isolated positional matches in a tangent plane.

    The separation radius is a caller-selected motion/search policy, not a
    calibration constant or a posterior identity probability. If either source
    has multiple candidates inside the radius, no association is made. Returned
    chi2 describes the static-position hypothesis under independent diagonal
    errors; high chi2 alone does not establish motion. Source detections and
    uncertainty calibration are supplied externally, not inferred here.
    """
    a,b=np.asarray(positions_a,dtype=float),np.asarray(positions_b,dtype=float)
    sa,sb=np.asarray(sigma_a,dtype=float),np.asarray(sigma_b,dtype=float)
    if a.ndim!=2 or b.ndim!=2 or a.shape[1]!=2 or b.shape[1]!=2 or sa.shape!=a.shape or sb.shape!=b.shape:
        raise ValueError('positions and uncertainties must be N by 2 tangent-plane arrays')
    if not all(np.isfinite(x).all() for x in (a,b,sa,sb)) or np.any(sa<=0) or np.any(sb<=0):
        raise ValueError('finite positions and positive positional uncertainties required')
    if not np.isfinite(max_separation_arcsec) or max_separation_arcsec<=0:
        raise ValueError('positive explicit association separation required')
    # The pairwise matrix is deliberately bounded to avoid hidden quadratic RAM.
    if a.shape[0]*b.shape[0]>1000000: raise ValueError('association exceeds one million candidate pairs')
    delta=a[:,None,:]-b[None,:,:]; distance=np.linalg.norm(delta,axis=2)
    eligible=distance<=max_separation_arcsec
    count_a=eligible.sum(axis=1);count_b=eligible.sum(axis=0)
    accepted=eligible & (count_a[:,None]==1) & (count_b[None,:]==1)
    matches=[]
    for i,j in np.argwhere(accepted):
        variance=sa[i]**2+sb[j]**2
        matches.append({'index_a':int(i),'index_b':int(j),
                        'separation_arcsec':float(distance[i,j]),
                        'position_chi2':float(np.sum(delta[i,j]**2/variance))})
    used_a={m['index_a'] for m in matches};used_b={m['index_b'] for m in matches}
    return {'matches':matches,
            'unmatched_a':[i for i in range(len(a)) if i not in used_a],
            'unmatched_b':[i for i in range(len(b)) if i not in used_b],
            'ambiguous_a':np.flatnonzero(count_a>1).tolist(),
            'ambiguous_b':np.flatnonzero(count_b>1).tolist(),
            'max_separation_arcsec':float(max_separation_arcsec),
            'limitations':['Rejects positional ambiguity rather than selecting a preferred identity.',
            'No flux, color, PSF, blending or known-object catalog confirmation.',
            'Pairwise motion hypotheses require later multi-epoch fits.']}
