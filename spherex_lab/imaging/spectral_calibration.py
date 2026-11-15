"""SPHEREx cutout detector-coordinate and quality helpers.

See https://caltech-ipac.github.io/spherex-archive-documentation/spherex-data-products/
for the cutout alternate-WCS convention and https://irsa.ipac.caltech.edu/
onlinehelp/spherex/spherex/sp.html for the SOURCE bit.
"""
import numpy as np

SOURCE = 1 << 21

def detector_indices(header, shape, detector_shape=(2040, 2040)):
    """Return zero-based original detector (y, x) indices for a cutout grid.

    Fail closed for non-integral or unsupported alternate-WCS mappings.
    """
    if len(shape) != 2 or any(not isinstance(v, (int, np.integer)) or v <= 0 for v in shape):
        raise ValueError('invalid cutout shape')
    try:
        if header['CTYPE1A'].strip() != 'LINEAR' or header['CTYPE2A'].strip() != 'LINEAR':
            raise ValueError('unsupported detector coordinate system')
        if float(header['CDELT1A']) != 1 or float(header['CDELT2A']) != 1 or float(header['CRVAL1A']) != 0 or float(header['CRVAL2A']) != 0:
            raise ValueError('unsupported detector coordinate transform')
        if any(k in header for k in ('CD1_1A','CD1_2A','CD2_1A','CD2_2A')):
            raise ValueError('unsupported detector CD matrix')
        if any(float(header.get(k, v)) != v for k, v in [('PC1_1A',1),('PC1_2A',0),('PC2_1A',0),('PC2_2A',1)]):
            raise ValueError('unsupported detector coordinate matrix')
        x0 = 1 - float(header['CRPIX1A'])
        y0 = 1 - float(header['CRPIX2A'])
    except (KeyError, TypeError) as exc:
        raise ValueError('missing detector coordinate metadata') from exc
    if not np.isfinite(x0) or not np.isfinite(y0) or not x0.is_integer() or not y0.is_integer():
        raise ValueError('non-integral detector offset')
    x = np.arange(shape[1], dtype=np.int64) + int(x0)
    y = np.arange(shape[0], dtype=np.int64) + int(y0)
    if x[0] < 0 or y[0] < 0 or x[-1] >= detector_shape[1] or y[-1] >= detector_shape[0]:
        raise ValueError('cutout extends outside calibration detector grid')
    return y, x


def sample_calibration(header, shape, cwave, cband):
    """Extract already-verified associated detector calibration arrays in microns.

    Association between calibration and science image MUST be independently verified
    by caller; filename/detector equality alone does not establish association.
    """
    cwave, cband = np.asarray(cwave), np.asarray(cband)
    if cwave.ndim != 2 or cwave.shape != cband.shape:
        raise ValueError('incompatible calibration array shapes')
    y, x = detector_indices(header, shape, cwave.shape)
    wave = cwave[np.ix_(y, x)].astype(float)
    band = cband[np.ix_(y, x)].astype(float)
    valid = np.isfinite(wave) & np.isfinite(band) & (wave > 0) & (band > 0)
    return wave, band, valid


def summarize_flags(flags):
    """Report exact bit counts and explicitly identify the SOURCE bit.

    SOURCE alone is nominal in IRSA spectrophotometry, but can indicate
    contamination; this function does not silently change a mask policy.
    """
    a = np.asarray(flags)
    if a.dtype.kind not in 'iu' or np.any(a < 0):
        raise ValueError('flags must be nonnegative integers')
    a = a.astype(np.uint64)
    counts = {str(bit): int(np.count_nonzero(a & np.uint64(1 << bit))) for bit in range(64) if np.any(a & np.uint64(1 << bit))}
    return {'total': int(a.size), 'zero': int(np.count_nonzero(a == 0)), 'source_only': int(np.count_nonzero(a == SOURCE)), 'counts_by_bit': counts}
