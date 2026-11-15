import io
import numpy as np
import pytest

def test_parse_curve_and_invalid():
 from spherex_lab.imaging.response_calibration import parse_response_csv
 w,r=parse_response_csv(b'Wavelength_um,NormTransmission\n1,0\n2,1\n3,0\n')
 np.testing.assert_array_equal(w,[1,2,3]);np.testing.assert_array_equal(r,[0,1,0])
 with pytest.raises(ValueError):parse_response_csv(b'Wavelength_um,NormTransmission\n2,1\n1,0\n')

def test_overlap_does_not_assert_equivalence():
 from spherex_lab.imaging.response_calibration import response_metrics
 m=response_metrics([1,2,3],[0,1,0],[1,2,3],[0,1,0])
 assert m['normalized_overlap']==pytest.approx(1) and m['exact_samples_equal']
 m=response_metrics([1,2,3],[0,1,0],[3,4,5],[0,1,0])
 assert m['normalized_overlap']==0 and not m['exact_samples_equal']

def test_fiducial_shift_preserves_shape():
 from spherex_lab.imaging.response_calibration import shifted_fiducial
 w,r,p=shifted_fiducial([1,2,3],[0,1,0],target_center_um=4,fiducial_center_um=2)
 np.testing.assert_array_equal(w,[3,4,5]);assert p['status']=='APPROXIMATE'
 assert p['response_systematic_uncertainty']=='unquantified'

def test_measured_negative_response_samples_preserved():
 from spherex_lab.imaging.response_calibration import parse_response_csv,response_metrics
 w,r=parse_response_csv(b'Wavelength_um,NormTransmission\n1,-0.01\n2,1\n3,0\n')
 assert r[0]==-.01
 with pytest.raises(ValueError):response_metrics(w,r,w,r)

def test_pixel_mapping_and_mwir_flip(tmp_path):
 from spherex_lab.imaging.response_calibration import pixel_response
 import tarfile
 from astropy.io import fits
 blob=io.BytesIO(); channel=np.ones((2040,2040),dtype='uint8');channel[1,2038]=2
 fits.HDUList([fits.PrimaryHDU(),fits.ImageHDU(channel)]).writeto(blob)
 path=tmp_path/'bundle.tar'
 def add(t,name,data):
  m=tarfile.TarInfo(name);m.size=len(data);t.addfile(m,io.BytesIO(data))
 with tarfile.open(path,'w') as t:
  add(t,'ChannelMaps/spectral_channels_17ch_2040.fits',blob.getvalue())
  add(t,'FidSpecSummary/SubChannel17/Band6/FidSpecSummary_Band6_SpecCh17.csv',b'SpecCh,BandCenter_um\n2,2\n')
  add(t,'MiddlePixel/SubChannel17/Band6/FiducialSpec_Band6_SubCh02_17_Cent.csv',b'Wavelength_um,NormTransmission\n1,0\n2,1\n3,0\n')
 h={'DETECTOR':6,'CTYPE1A':'LINEAR','CTYPE2A':'LINEAR','CDELT1A':1,'CDELT2A':1,'CRVAL1A':0,'CRVAL2A':0,'CRPIX1A':0,'CRPIX2A':0}
 w,r,p=pixel_response(path,h,(2,2),0,0,4)
 assert p['channel']==2 and p['laboratory_x']==2038
 np.testing.assert_array_equal(w,[3,4,5])
 with pytest.raises(ValueError):pixel_response(path,h,(2,2),5,0,4)

def test_overlap_discontinuous_support_exact_integral():
 from spherex_lab.imaging.response_calibration import response_metrics
 m=response_metrics([1,3],[1,1],[2,4],[1,1])
 assert m['normalized_overlap']==pytest.approx(.5)
 assert m['centroid_difference_um']==pytest.approx(1)
