"""Surface-brightness WCS resampling with explicit independent-pixel variance.

Bilinear interpolation is not PSF matching and introduces covariance between
neighboring output pixels. Returned variance is the diagonal only; significance
of sums or fitted sources requires covariance treatment beyond this operator.
"""
import numpy as np


def reproject_bilinear(image, variance, valid, source_wcs, target_wcs, shape_out):
    """Resample to target WCS; require every nonzero-weight input to be valid.

    No extrapolation and no interpolation across invalid neighbors. `footprint`
    describes geometric support; `valid` additionally includes science masks.
    Both WCS arguments must be celestial two-dimensional Astropy WCS objects.
    """
    image = np.asarray(image, dtype=float)
    variance = np.asarray(variance, dtype=float)
    valid = np.asarray(valid, dtype=bool)
    if image.ndim != 2 or image.shape != variance.shape or image.shape != valid.shape:
        raise ValueError('image, variance and valid must share a two-dimensional grid')
    if len(shape_out) != 2 or any(not isinstance(n, (int, np.integer)) or n <= 0 for n in shape_out):
        raise ValueError('invalid output shape')
    if not source_wcs.has_celestial or not target_wcs.has_celestial or source_wcs.pixel_n_dim != 2 or target_wcs.pixel_n_dim != 2:
        raise ValueError('two-dimensional celestial WCS required')
    valid = valid & np.isfinite(image) & np.isfinite(variance) & (variance >= 0)
    yy, xx = np.indices(shape_out, dtype=float)
    world = target_wcs.pixel_to_world(xx, yy)
    x, y = source_wcs.world_to_pixel(world)
    finite = np.isfinite(x) & np.isfinite(y)
    # WCS floating-point roundoff must not discard an exactly coincident edge.
    x = np.where(finite, x, -100.)
    y = np.where(finite, y, -100.)
    x = np.where(np.abs(x-np.rint(x)) < 1e-7, np.rint(x), x)
    y = np.where(np.abs(y-np.rint(y)) < 1e-7, np.rint(y), y)
    x0 = np.floor(x).astype(np.int64)
    y0 = np.floor(y).astype(np.int64)
    fx, fy = x-x0, y-y0
    footprint = finite.copy()
    accepted = finite.copy()
    result = np.zeros(shape_out)
    var_result = np.zeros(shape_out)
    for dx, dy, weight in ((0,0,(1-fx)*(1-fy)),(1,0,fx*(1-fy)),(0,1,(1-fx)*fy),(1,1,fx*fy)):
        ix, iy = x0+dx, y0+dy
        needed = weight > 0
        inside = (ix >= 0) & (ix < image.shape[1]) & (iy >= 0) & (iy < image.shape[0])
        cx, cy = np.clip(ix,0,image.shape[1]-1), np.clip(iy,0,image.shape[0]-1)
        good = inside & valid[cy,cx]
        footprint &= ~needed | inside
        accepted &= ~needed | good
        result += np.where(good & needed, image[cy,cx], 0.) * weight
        var_result += np.where(good & needed, variance[cy,cx], 0.) * weight**2
    result[~accepted] = np.nan
    var_result[~accepted] = np.nan
    return {'image':result, 'variance':var_result, 'valid':accepted, 'footprint':footprint}
