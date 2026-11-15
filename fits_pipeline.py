"""Read validated SPHEREx science arrays without assuming spectral compatibility."""
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

def read_cutout(path, allow_source=False):
    with fits.open(path, memmap=False) as hdus:
        return read_cutout_hdus(hdus, allow_source)


def read_cutout_hdus(hdus, allow_source=False):
    """Copy validated science arrays from a caller-owned open FITS product."""
    try:
        image = np.array(hdus['IMAGE'].data, dtype=float)
        variance = np.array(hdus['VARIANCE'].data, dtype=float)
        flags = np.array(hdus['FLAGS'].data)
        header = hdus['IMAGE'].header.copy()
    except (KeyError, TypeError) as exc:
        raise ValueError('Missing required SPHEREx IMAGE, VARIANCE or FLAGS extension') from exc
    if image.ndim != 2 or variance.shape != image.shape or flags.shape != image.shape or flags.dtype.kind not in 'iu':
        raise ValueError('Incompatible image, variance or integer flag arrays')
    if header.get('BUNIT', '').replace(' ', '').lower() != 'mjy/sr':
        raise ValueError('Unsupported science image unit')
    if hdus['VARIANCE'].header.get('BUNIT', '').replace(' ', '').lower() != 'mjy2/sr2':
        raise ValueError('Unsupported variance unit')
    valid = np.isfinite(image) & np.isfinite(variance) & (variance >= 0) & ((flags == 0) | ((flags == (1 << 21)) & allow_source))
    return {'image': image, 'variance': variance, 'flags': flags, 'valid': valid, 'wcs': WCS(header), 'header': header}
