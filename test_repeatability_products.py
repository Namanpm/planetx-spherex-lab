import json
import numpy as np
from astropy.io import fits
from repeatability_products import export

def test_export_predictions_and_covariance(tmp_path):
    a={i:{'flux_jy':2.,'color':color} for i,color in [('train',1.),('val',1.5),('confirm',2.)]}
    b={i:{'flux_jy':2*np.exp(.1+.2*(r['color']-1))} for i,r in a.items()}
    for name,data in [('a',a),('b',b)]: (tmp_path/(name+'.json')).write_text(json.dumps({'measurements':data}),encoding='utf-8')
    model={'coefficients':[.1,.2],'reference_color':1.,'retained_training_ids':['train'],'heldout':{'source_ids':['confirm'],'predictive_covariance':[[.04]],'residual_log_ratio':[0.],'standardized_residual':[0.]}}
    report={'status':'FIXTURE','observations':[{'publisher_id':n,'measurement_file':str(tmp_path/(n+'.json'))} for n in ['a','b']],'confirmation':{'selected_pair':['a','b'],'model':model},'pair_evaluations':[{'observation_a':'a','observation_b':'b','validation_model':{'heldout':{'source_ids':['val']}}}]}
    path=tmp_path/'report.json';path.write_text(json.dumps(report),encoding='utf-8');export(path)
    with fits.open(tmp_path/'heldout_contrasts.fits',checksum=True) as f:
        np.testing.assert_allclose(f[1].data['PREDICTED_FLUX_B'],2*np.exp(.3))
        np.testing.assert_allclose(f[2].data,[[.04]])
    assert (tmp_path/'confirmation_diagnostic.png').exists()
