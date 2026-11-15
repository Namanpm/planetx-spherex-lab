"""Numerical difference diagnostics; response equality requires full supplied curves.

CWAVE and CBAND alone are never accepted as evidence of equal responses.
"""
from pathlib import Path
import json
import numpy as np
from astropy.io import fits


def response_compatibility(wavelength_a, response_a, wavelength_b, response_b):
    """Exact equality of caller-validated response samples on identical grids.

    No central-wavelength tolerance or unsupported response interpolation is used.
    The caller must validate provenance, response sampling and normalization.
    """
    wa, wb = np.asarray(wavelength_a), np.asarray(wavelength_b)
    a, b = np.asarray(response_a), np.asarray(response_b)
    if not np.isfinite(wa).all() or not np.isfinite(wb).all() or wa.ndim != 1 or wa.size < 2 or not np.all(np.diff(wa) > 0):
        raise ValueError('response wavelengths must be an increasing grid')
    if a.shape != b.shape or a.shape[-1] != wa.size or wa.shape != wb.shape:
        raise ValueError('incompatible response shapes')
    if not np.array_equal(wa, wb):
        return np.zeros(a.shape[:-1], dtype=bool)
    valid = np.isfinite(a).all(axis=-1) & np.isfinite(b).all(axis=-1)
    valid &= (a >= 0).all(axis=-1) & (b >= 0).all(axis=-1)
    valid &= (a.sum(axis=-1) > 0) & (b.sum(axis=-1) > 0)
    return valid & (a == b).all(axis=-1)


def compare_aligned(image_a, variance_a, image_b, variance_b, valid_a, valid_b,
                    spectral_mask, *, blank_mask=None, min_overlap=0.0,
                    psf_status='unverified', synthetic=False, sigma_threshold=5.0, resampled=False):
    """B minus A on a common grid, with independent-input marginal variance.

    min_overlap and sigma_threshold are explicit analysis policies, not mission
    calibration constants. Significance is exploratory and does not include
    unknown calibration covariance or PSF systematics. Background fitting uses
    only caller-identified blank sky; no source field is automatically flattened.
    """
    a, va, b, vb = [np.asarray(x, dtype=float) for x in (image_a, variance_a, image_b, variance_b)]
    if a.ndim != 2 or any(x.shape != a.shape for x in (va,b,vb,valid_a,valid_b,spectral_mask)):
        raise ValueError('aligned arrays and masks must have identical 2D shapes')
    if not 0 <= min_overlap <= 1 or not np.isfinite(sigma_threshold) or sigma_threshold <= 0:
        raise ValueError('invalid overlap or candidate threshold')
    spatial = np.asarray(valid_a,dtype=bool) & np.asarray(valid_b,dtype=bool)
    valid = spatial & np.asarray(spectral_mask,dtype=bool)
    valid &= np.isfinite(a) & np.isfinite(b) & np.isfinite(va) & np.isfinite(vb) & (va >= 0) & (vb >= 0)
    overlap = float(spatial.sum()/a.size)
    if not valid.any() or overlap < min_overlap:
        raise ValueError('insufficient valid spatial/spectral overlap')
    difference = b-a
    variance = va+vb
    background = np.zeros_like(a)
    background_covariance = None
    if blank_mask is not None and resampled:
        raise ValueError('background fitting on resampled pixels requires full covariance')
    if blank_mask is not None:
        blank = np.asarray(blank_mask,dtype=bool)
        if blank.shape != a.shape: raise ValueError('blank mask shape mismatch')
        usable = blank & valid & (variance > 0)
        y,x = np.indices(a.shape)
        design = np.stack([np.ones_like(a),x,y],axis=-1)
        z = design[usable]
        if len(z) < 4: raise ValueError('insufficient independent blank pixels')
        weights = 1/variance[usable]
        normal = z.T @ (weights[:,None]*z)
        if np.linalg.matrix_rank(normal) != 3: raise ValueError('blank mask does not constrain background plane')
        background_covariance = np.linalg.inv(normal)
        coefficients = background_covariance @ (z.T @ (weights*difference[usable]))
        background = design @ coefficients
        difference = difference-background
        prediction_var = np.einsum('...i,ij,...j->...',design,background_covariance,design)
        # Fit and residual share the blank pixels: subtract twice their covariance.
        variance = variance + prediction_var - 2*prediction_var*usable
        variance = np.maximum(variance,0)
    difference = np.where(valid,difference,np.nan)
    variance = np.where(valid,variance,np.nan)
    uncertainty = np.sqrt(variance)
    significance = np.divide(difference,uncertainty,out=np.full_like(a,np.nan),where=valid & (uncertainty>0))
    candidates = []
    for y,x in np.argwhere(valid & (np.abs(significance) >= sigma_threshold)):
        candidates.append({'x':int(x),'y':int(y),'signed_sigma':float(significance[y,x]),'classification':'unclassified residual pixel'})
    candidates.sort(key=lambda p: abs(p['signed_sigma']), reverse=True)
    report = {'status':'SYNTHETIC' if synthetic else 'PARTIAL',
              'spatial_valid_fraction':overlap,'usable_fraction':float(valid.mean()),
              'background_model':'weighted blank-sky plane' if blank_mask is not None else 'none',
              'background_covariance':None if background_covariance is None else background_covariance.tolist(),
              'psf_status':psf_status,'candidate_threshold_sigma':sigma_threshold,
              'candidates':candidates,'limitations':['Independent observation errors assumed.',
              'Significance excludes unknown calibration covariance and PSF systematics.',
              'Neighboring resampled pixels can be correlated; candidates are pixels, not detections.']}
    return {'aligned_a':a,'aligned_b':b,'difference':difference,'variance':variance,
            'uncertainty':uncertainty,'significance':significance,'valid':valid,
            'spatial_mask':spatial,'spectral_mask':np.asarray(spectral_mask,dtype=bool),
            'background':background,'report':report}


def write_outputs(result, output_dir, wcs=None, provenance=None):
    """Write numerical diagnostics with units and celestial WCS where supplied."""
    out=Path(output_dir); out.mkdir(parents=True,exist_ok=True)
    header = wcs.to_header(relax=True) if wcs is not None else fits.Header()
    hdus=[fits.PrimaryHDU()]
    for name in ('aligned_a','aligned_b','difference','variance','uncertainty','significance','valid','spatial_mask','spectral_mask','background'):
        data=result[name]
        h=header.copy()
        if data.dtype.kind=='b': data=data.astype('uint8')
        h['BUNIT']='MJy2/sr2' if name=='variance' else ('1' if name in ('significance','valid','spatial_mask','spectral_mask') else 'MJy/sr')
        hdus.append(fits.ImageHDU(data,header=h,name=name.upper()))
    if 'full_covariance' in result:
        covariance_hdu=fits.ImageHDU(result['full_covariance'],name='COVARIANCE')
        covariance_hdu.header['BUNIT']='MJy2/sr2'
        covariance_hdu.header['COVORD']='row-major flattened image pixels'
        hdus.append(covariance_hdu)
    fits.HDUList(hdus).writeto(out/'comparison.fits',overwrite=True,checksum=True)
    report={**result['report'],'provenance':provenance or {}}
    (out/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False),encoding='utf-8')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,3,figsize=(12,4))
    for ax,key,label in zip(axes,('aligned_a','aligned_b','difference'),('A (MJy/sr)','B (MJy/sr)','B − A (MJy/sr)')):
        fig.colorbar(ax.imshow(result[key],origin='lower'),ax=ax,shrink=.65)
        ax.set_title(label)
    fig.suptitle(report['status']+' — diagnostic; residuals are not discoveries')
    fig.tight_layout(); fig.savefig(out/'diagnostic.png',dpi=120); plt.close(fig)
    return str(out/'report.json')
