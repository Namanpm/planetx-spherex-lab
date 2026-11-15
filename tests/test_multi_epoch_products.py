import hashlib
import json
import numpy as np
from astropy.io import fits
from spherex_lab.repeatability.multi_epoch_products import export


def data():
    ids=['obsA','obsB']
    report={'selected_observation_ids':ids,'readiness':[{'observation_id':'obsA','epoch':'2025-01-01','mjd':60676.,'status':'READY','detector':1,'valid_measurement_count':1},
            {'observation_id':'obsB','epoch':'2025-01-02','mjd':60677.,'status':'FAILED','detector':1,'valid_measurement_count':0,'failure':{'reason':'no calibration'}}],
            'pair_evaluations':[{'observation_a':'obsA','observation_b':'obsB','status':'REJECTED','reason':'failed calibration'}],
            'time_series':{'status':'INSUFFICIENT_EVIDENCE','reason':'No validated common scale','timelines':[]}}
    packets={'obsA':{'record':{'t_min':60676.},'date':'2025-01-01','measurements':{'s1':{'ra_deg':1.,'dec_deg':2.,'flux_jy':.003,'variance_jy2':1e-8,'flag_values':[0]}}},
             'obsB':{'record':{'t_min':60677.},'date':'2025-01-02','measurements':{},'error':{'reason':'no calibration'}}}
    return report,packets


def test_artifact_integrity_preserves_missing_rejected_and_raw_flux(tmp_path):
    report,packets=data();names=export(report,packets,tmp_path)
    assert all('/' not in name and '\\' not in name for name in names)
    manifest=json.loads((tmp_path/'artifact_manifest.json').read_text())
    for name,m in manifest.items():
        assert hashlib.sha256((tmp_path/name).read_bytes()).hexdigest()==m['sha256']
        assert (tmp_path/name).stat().st_size==m['bytes']
    with fits.open(tmp_path/'source_measurements.fits',checksum=True) as hdus:
        hdus.verify('exception');table=hdus[1].data
        assert len(table)==2
        assert table[0]['FLUX_JY']==.003
        assert table[1]['STATUS']=='REJECTED' and np.isnan(table[1]['FLUX_JY'])
        assert 'no calibration' in table[1]['REASON']
    with fits.open(tmp_path/'observation_readiness.fits',checksum=True) as hdus:
        assert hdus[1].data['STATUS'].tolist()==['READY','FAILED']
    assert not (tmp_path/'timeline_covariance.fits').exists()
    assert (tmp_path/'light_curve.png').read_bytes().startswith(b'\x89PNG')


def test_accepted_timeline_exports_real_points_and_covariance(tmp_path):
    report,packets=data()
    point={'observation_id':'obsA','source_id':'s1','date':'2025-01-01','mjd':60676.,'measured_flux_jy':.003,'calibrated_flux_jy':.002,'statistical_variance_jy2':1e-8,'calibration_variance_jy2':2e-8}
    report['time_series']={'status':'SUCCESS','timelines':[{'source_id':'s1','points':[{**point,'observation_id':f'obs{j}','mjd':60676.+j} for j in range(3)],'flux_covariance_jy2':(np.eye(3)*3e-8).tolist(),'missing_epochs':[]}]}
    names=export(report,packets,tmp_path)
    assert 'timeline_covariance.fits' in names and 'timeline_measurements.fits' in names
    with fits.open(tmp_path/'timeline_covariance.fits',checksum=True) as hdus:np.testing.assert_array_equal(hdus[1].data,np.eye(3)*3e-8)
    with fits.open(tmp_path/'timeline_measurements.fits',checksum=True) as hdus:
        assert hdus[1].data[0]['CALIBRATED_FLUX_JY']==.002
        assert hdus[1].data[0]['TOTAL_VARIANCE_JY2']==3e-8
