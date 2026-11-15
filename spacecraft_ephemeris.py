"""Authoritative JPL Horizons geometric SPHEREx observer vectors, ICRF AU/TDB."""
import re
import csv
import hashlib
import json
import os
import tempfile
from pathlib import Path
import numpy as np
from astropy.time import Time
from calibration_resolver import session
URL='https://ssd.jpl.nasa.gov/api/horizons.api'
SPACECRAFT='-163182' # JPL Horizons name lookup SPHEREx, saved in artifacts.


def parse_vector(payload, jd_tdb):
    if payload.get('error'): raise ValueError('Horizons: '+str(payload['error']))
    text=payload.get('result','')
    if 'Target body name: SPHEREx' not in text or '(-163182)' not in text or 'Solar System Barycenter' not in text or 'ICRF' not in text: raise ValueError('unexpected ephemeris identity/frame')
    if not re.search(r'Output units\s*:\s*AU-D',text) or 'Geometric state vectors have NO corrections' not in text:
        raise ValueError('unexpected observer units or correction')
    if '$$SOE' not in text or '$$EOE' not in text: raise ValueError('no spacecraft vector coverage')
    lines=[line.strip() for line in text.split('$$SOE')[1].split('$$EOE')[0].splitlines() if line.strip()]
    if len(lines)!=1: raise ValueError('expected one observer epoch')
    row=next(csv.reader(lines)); epoch=float(row[0]); xyz=np.array([float(v) for v in row[2:5]])
    if len(xyz)!=3 or not np.isfinite(xyz).all() or abs(epoch-jd_tdb)>1e-7: raise ValueError('observer epoch/vector mismatch')
    return xyz.tolist()


def observer_vector(obstime,cache_dir):
    jd=float(Time(obstime).tdb.jd)
    parameters={'format':'json','COMMAND':"'-163182'",'CENTER':"'500@0'",'EPHEM_TYPE':"'VECTORS'",'TLIST':str(jd),'REF_PLANE':"'FRAME'",'REF_SYSTEM':"'ICRF'",'OUT_UNITS':"'AU-D'",'VEC_CORR':"'NONE'",'CSV_FORMAT':"'YES'"}
    directory=Path(cache_dir);directory.mkdir(parents=True,exist_ok=True)
    key=hashlib.sha256(json.dumps(parameters,sort_keys=True).encode()).hexdigest();path=directory/(key+'.json')
    if path.exists(): payload=json.loads(path.read_text(encoding='utf-8'))
    else:
        with session() as client:
            response=client.get(URL,params=parameters,timeout=(10,30));response.raise_for_status();payload=response.json()
        parse_vector(payload,jd)
        with tempfile.NamedTemporaryFile(mode='w',encoding='utf-8',dir=directory,delete=False) as f:
            json.dump(payload,f); temporary=f.name
        os.replace(temporary,path)
    xyz=parse_vector(payload,jd)
    return xyz,{'source':URL,'target':'SPHEREx (-163182)','epoch_jd_tdb':jd,'frame':'ICRF','origin':'Solar System Barycenter','unit':'AU','correction':'geometric NONE','response_sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
