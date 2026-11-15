import numpy as np
import pytest
from astropy.io import fits
from calibration_resolver import discover, validate_calibration

PID='ivo://irsa.ipac/spherex_dr1?2025W18_1B_0162_1/D1'
class Response:
    def __init__(self, body): self.content=body
    def raise_for_status(self): pass
class Client:
    def __init__(self, body): self.body=body
    def get(self, *a, **k): return Response(self.body)

def xml(rows):
    fields=''.join(f'<FIELD name="{n}" datatype="char" arraysize="*"/>' for n in ['id','semantics','description','access_url','error_message'])
    return ('<VOTABLE version="1.3" xmlns="http://www.ivoa.net/xml/VOTable/v1.3"><RESOURCE><TABLE>'+fields+'<DATA><TABLEDATA>'+rows+'</TABLEDATA></DATA></TABLE></RESOURCE></VOTABLE>').encode()
def row(url):
    return '<TR>'+''.join('<TD>'+s+'</TD>' for s in [PID,'#calibration','spectral_wcs_D1',url,''])+'</TR>'
URL='https://irsa.ipac.caltech.edu/ibe/data/spherex/dr1/spectral_wcs/a.fits'
def test_datalink_positive(): assert discover(PID,1,Client(xml(row(URL)))) == URL
def test_missing_and_ambiguous():
    for r in ('',row(URL)*2):
        with pytest.raises(ValueError): discover(PID,1,Client(xml(r)))
def test_bad_response_and_untrusted_url():
    with pytest.raises(Exception): discover(PID,1,Client(b'not xml'))
    with pytest.raises(ValueError): discover(PID,1,Client(xml(row('https://evil.example/a'))))
def test_wrong_detector():
    with pytest.raises(ValueError): discover(PID,2,Client(xml(row(URL))))
def test_calibration_detector_and_missing(tmp_path):
    p=tmp_path/'x.fits'
    h=fits.ImageHDU(np.zeros((2040,2040),dtype='float32'),name='CWAVE'); h.header['DETECTOR']='2'; h.header['BUNIT']='um'
    h2=h.copy(); h2.name='CBAND'; fits.HDUList([fits.PrimaryHDU(),h,h2]).writeto(p)
    with pytest.raises(ValueError,match='detector'): validate_calibration(p,1)
    with pytest.raises(ValueError,match='missing'): 
        q=tmp_path/'bad.fits'; fits.PrimaryHDU().writeto(q); validate_calibration(q,1)
