"""Bounded authoritative Gaia DR3 reference-candidate discovery.

NOT_AVAILABLE variability is not a claim of constant flux. Science-image masks,
blending/saturation checks and independent repeatability validation remain required.
"""
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import tempfile
from datetime import datetime,timezone
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.time import Time
import astropy.units as u
from archive_http import retry_session as session

ENDPOINT='https://gea.esac.esa.int/tap-server/tap/sync'
MAX_BYTES=8*1024*1024
COLUMNS=('source_id','ra','dec','pmra','pmdec','ref_epoch','phot_variable_flag','duplicated_source','bp_rp','ruwe','phot_g_mean_mag','ra_error','dec_error','pmra_error','pmdec_error','parallax','parallax_error')


def screen_rows(rows):
    sources=[];rejected=[];seen=set()
    for row in rows:
        ident=str(row.get('source_id',''))
        reason=None
        try:
            record={**row,**{key:float(row[key]) for key in ('ra','dec','pmra','pmdec','ref_epoch','bp_rp','ruwe','phot_g_mean_mag')}}
            if not all(np.isfinite(record[key]) for key in ('ra','dec','pmra','pmdec','ref_epoch','bp_rp','ruwe','phot_g_mean_mag')) or not 0<=record['ra']<360 or not -90<=record['dec']<=90:reason='invalid astrometry/color/photometry'
            elif str(row.get('phot_variable_flag','')).upper()=='VARIABLE':reason='Gaia classified variable'
            elif str(row.get('duplicated_source','')).lower() not in ('false','0'):reason='duplicated or unknown duplicate flag'
            elif ident in seen:reason='duplicate source id'
            elif not ident:reason='missing source id'
            else:
                for key in ('ra_error','dec_error','pmra_error','pmdec_error','parallax','parallax_error'):
                    try:value=float(row.get(key,''))
                    except (ValueError,TypeError):value=None
                    record[key]=value if value is not None and np.isfinite(value) else None
        except (ValueError,TypeError,KeyError): reason='missing or malformed astrometry/color/photometry'
        if reason:rejected.append({'source_id':ident,'reason':reason});continue
        seen.add(ident)
        record.update(source_id=ident,frame='ICRS',stability_status='not classified variable by Gaia; constancy unproven')
        sources.append(record)
    return sources,rejected


def query_references(ra,dec,radius=.05,cache_dir='artifacts/reference_catalog',limit=2000,max_g_mag=19.):
    if not np.isfinite([ra,dec,radius]).all() or not 0<=ra<360 or not -90<=dec<=90 or not 0<radius<=.2:raise ValueError('invalid catalog query coordinates/radius')
    if not isinstance(limit,int) or not 1<=limit<=10000:raise ValueError('catalog limit must be 1..10000')
    if not np.isfinite(max_g_mag):raise ValueError('invalid magnitude selection')
    query=f"SELECT TOP {limit+1} {','.join(COLUMNS)} FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{float(ra)},{float(dec)},{float(radius)})) AND phot_g_mean_mag<={float(max_g_mag)} ORDER BY phot_g_mean_mag ASC"
    parameters={'REQUEST':'doQuery','LANG':'ADQL','FORMAT':'csv','QUERY':query}
    directory=Path(cache_dir);directory.mkdir(parents=True,exist_ok=True)
    key=hashlib.sha256(json.dumps(parameters,sort_keys=True).encode()).hexdigest();path=directory/(key+'.json')
    if path.exists():
        cached=json.loads(path.read_text(encoding='utf-8'))
        content=cached['csv'].encode('utf-8');retrieved=cached['retrieved_utc']
    else:
        try:
            with session() as client,client.get(ENDPOINT,params=parameters,timeout=(10,60),stream=True) as response:
                response.raise_for_status();chunks=[];size=0
                for chunk in response.iter_content(65536):
                    size+=len(chunk)
                    if size>MAX_BYTES:raise ValueError('Gaia response exceeds8MiB')
                    chunks.append(chunk)
                content=b''.join(chunks)
        except Exception as exc:
            import requests
            if isinstance(exc, (requests.exceptions.Timeout, requests.exceptions.ConnectionError,
                                requests.exceptions.RetryError, requests.exceptions.HTTPError)):
                raise GaiaArchiveUnavailable('Gaia DR3 catalog service did not respond after bounded retries; retry later, or select saved replay data.') from exc
            raise
        retrieved=datetime.now(timezone.utc).isoformat()
        text=content.decode('utf-8-sig')
        reader=csv.DictReader(io.StringIO(text))
        if reader.fieldnames is None or not set(COLUMNS)<=set(reader.fieldnames):raise ValueError('malformed Gaia TAP CSV or upstream error')
        rows=list(reader)
        if len(rows)>limit:raise ValueError('Gaia catalog limit exceeded; narrow field or raise explicit limit')
        with tempfile.NamedTemporaryFile(mode='w',encoding='utf-8',dir=directory,delete=False) as f:
            json.dump({'csv':text,'retrieved_utc':retrieved},f);temporary=f.name
        os.replace(temporary,path)
    reader=csv.DictReader(io.StringIO(content.decode('utf-8-sig')))
    if reader.fieldnames is None or not set(COLUMNS)<=set(reader.fieldnames):raise ValueError('malformed cached Gaia CSV')
    rows=list(reader)
    sources,rejected=screen_rows(rows)
    all_sources=[]
    for row in rows:
        try:
            ra_value,dec_value=float(row['ra']),float(row['dec'])
            if not np.isfinite([ra_value,dec_value]).all():continue
            item={**row,'ra':ra_value,'dec':dec_value}
            for key in ('pmra','pmdec','ref_epoch','phot_g_mean_mag','bp_rp'):
                try:value=float(row[key])
                except (ValueError,TypeError,KeyError):value=None
                item[key]=value if value is not None and np.isfinite(value) else None
            all_sources.append(item)
        except (ValueError,TypeError,KeyError):continue
    return {'sources':sources,'all_sources':all_sources,'rejected':rejected,'provenance':{'endpoint':ENDPOINT,'catalog':'Gaia DR3 gaia_source','adql':query,'sha256':hashlib.sha256(content).hexdigest(),'retrieved_utc':retrieved,'reference_epoch':'per-source Julian year (usually2016.0)','query_radius_deg':radius,'magnitude_selection':{'phot_g_mean_mag_max':max_g_mag,'status':'Explicit sample selection, not catalog completeness guarantee'},'truncated':False,'variability_policy':'Exclude VARIABLE; unknown variability does not prove constancy.','quality_policy':'Exclude duplicates/missing required astrometry/color/photometry. RUWE preserved without invented cutoff.','limitations':['Catalog sources are candidate references only.','No image blending/saturation or infrared variability screen is established by Gaia alone.','Neighbors fainter than the explicit Gaia G magnitude selection are not represented.']}}


class GaiaArchiveUnavailable(RuntimeError):
    """Finite-retry failure from the external Gaia TAP service."""


def source_position(row,obstime):
    """ICRS PM propagation; no observer parallax/RV, explicit approximation."""
    ra,dec=np.deg2rad([float(row['ra']),float(row['dec'])])
    dt=float(Time(obstime).jyear)-float(row['ref_epoch'])
    east=np.array([-np.sin(ra),np.cos(ra),0.])
    north=np.array([-np.sin(dec)*np.cos(ra),-np.sin(dec)*np.sin(ra),np.cos(dec)])
    direction=np.array([np.cos(dec)*np.cos(ra),np.cos(dec)*np.sin(ra),np.sin(dec)])
    tangent=np.deg2rad(1/3600000)*(float(row['pmra'])*east+float(row['pmdec'])*north)
    point=direction+dt*tangent
    if not np.isfinite(point).all():raise ValueError('invalid source proper motion')
    point/=np.linalg.norm(point)
    return {'ra_deg':float(np.rad2deg(np.arctan2(point[1],point[0]))%360),'dec_deg':float(np.rad2deg(np.arcsin(point[2]))),'limitations':'Rectilinear angular PM-only prediction; no observer parallax, perspective acceleration or radial velocity.'}
