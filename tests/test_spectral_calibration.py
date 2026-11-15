import numpy as np
import pytest
from astropy.io import fits
from spherex_lab.imaging.spectral_calibration import detector_indices, sample_calibration, summarize_flags

def test_real_cutout_offsets_and_flags():
    with fits.open('tests/data/sample.fits') as f:
        y, x = detector_indices(f['IMAGE'].header, f['IMAGE'].data.shape)
        assert (x[0], x[-1], y[0], y[-1]) == (742, 799, 2012, 2039)
        s = summarize_flags(f['FLAGS'].data)
        assert s['zero'] == 54
        assert s['source_only'] > 0

def test_indexing_and_fail_closed():
    h = {'CTYPE1A':'LINEAR','CTYPE2A':'LINEAR','CDELT1A':1,'CDELT2A':1,'CRVAL1A':0,'CRVAL2A':0,'CRPIX1A':0,'CRPIX2A':0}
    a = np.arange(36).reshape(6,6) + 1
    w,b,v = sample_calibration(h,(2,3),a,a+1)
    np.testing.assert_array_equal(w, a[1:3,1:4])
    assert v.all()
    with pytest.raises(ValueError): detector_indices({**h,'PC1_2A':1},(2,3))
    with pytest.raises(ValueError): detector_indices(h,(6,6), (6,6))
