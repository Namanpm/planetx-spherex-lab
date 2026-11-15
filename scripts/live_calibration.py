import json
from pathlib import Path
import numpy as np
from astropy.io import fits
from spherex_lab.data.calibration_resolver import wavelength_maps
w,b,v,p = wavelength_maps('tests/data/sample.fits','ivo://irsa.ipac/spherex_dr1?2025W18_1B_0162_1/D1','artifacts/calibration_cache')
fits.HDUList([fits.PrimaryHDU(), fits.ImageHDU(w,name='CWAVE'), fits.ImageHDU(b,name='CBAND'), fits.ImageHDU(v.astype('uint8'),name='VALID')]).writeto('artifacts/wavelength_maps.fits',overwrite=True)
p['shape'] = list(w.shape); p['valid_pixels'] = int(v.sum()); p['wave_range_um'] = [float(w[v].min()), float(w[v].max())]; p['band_range_um'] = [float(b[v].min()), float(b[v].max())]
Path('artifacts/research/calibration_provenance.json').write_text(json.dumps(p,indent=2),encoding='utf-8')
print(json.dumps(p,indent=2))
