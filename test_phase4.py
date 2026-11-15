import numpy as np
import pytest
from astropy.io import fits
from fastapi.testclient import TestClient
from fits_pipeline import read_cutout
from main import app


def test_target_discovery():
    response = TestClient(app).get('/api/targets')
    assert response.status_code == 200
    assert len(response.json()) == 5
    wise = next(t for t in response.json() if t['id'] == 'wise-0855')
    assert wise['ra_deg'] > 0 and wise['pm_ra_cosdec_mas_yr'] < 0

def test_real_cutout_schema_and_conservative_mask():
    data = read_cutout('artifacts/sample.fits')
    assert data['image'].shape == (28, 58)
    assert np.all(data['flags'][data['valid']] == 0)
    assert np.isfinite(data['variance'][data['valid']]).all()

def test_missing_extension_fails(tmp_path):
    path = tmp_path / 'missing.fits'
    fits.PrimaryHDU().writeto(path)
    with pytest.raises(ValueError, match='Missing required'):
        read_cutout(path)

def test_source_only_policy_is_explicit():
    default=read_cutout('artifacts/sample.fits')
    source=read_cutout('artifacts/sample.fits',allow_source=True)
    assert source['valid'].sum() > default['valid'].sum()
    assert np.all(np.isin(source['flags'][source['valid']],[0,1<<21]))

def test_failed_fits_validation_releases_windows_file(tmp_path):
    path=tmp_path/'invalid.fits'
    image=fits.ImageHDU(np.ones((4,4),dtype='float32'),name='IMAGE')
    variance=fits.ImageHDU(np.ones((4,4),dtype='float32'),name='VARIANCE')
    flags=fits.ImageHDU(np.zeros((4,4),dtype='int32'),name='FLAGS')
    fits.HDUList([fits.PrimaryHDU(),image,variance,flags]).writeto(path)
    with pytest.raises(ValueError) as held_exception: read_cutout(path)
    # Keep the exception traceback alive while checking there is no mapped-file lock.
    assert held_exception.value is not None
    path.rename(tmp_path/'released.fits')
