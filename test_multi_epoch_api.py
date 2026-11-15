"""API contract tests use fake analysis jobs; no archive queries."""
import json
import sys
import threading
import time
import types
import uuid
from pathlib import Path

import pytest
from fastapi import FastAPI
from fastapi.testclient import TestClient
import multi_epoch_api as api


@pytest.fixture
def client(tmp_path,monkeypatch):
    monkeypatch.setattr(api,'RUN_ROOT',tmp_path/'runs')
    monkeypatch.setattr(api,'_slots',threading.BoundedSemaphore(1))
    monkeypatch.setattr(api,'_jobs',{})
    app=FastAPI();app.include_router(api.router)
    return TestClient(app)


def await_complete(client, ident):
    end=time.monotonic()+3
    while time.monotonic()<end:
        data=client.get('/api/multi-epoch/'+ident).json()
        if data['status'] in ('COMPLETED','FAILED'):return data
        time.sleep(.01)
    pytest.fail('Fake worker did not terminate')


def test_job_progress_results_and_artifacts(client,monkeypatch):
    def run(**kw):
        kw['progress']({'stage':'measuring','observation_id':'real-id'})
        directory=kw['output'];directory.mkdir(exist_ok=True,parents=True)
        (directory/'report.json').write_text('{"status":"SUCCESS"}')
        return {'status':'SUCCESS','readiness':[{'status':'READY'}],'pair_evaluations':[{'status':'REJECTED','reason':'incompatible'}],'time_series':{'status':'INSUFFICIENT_EVIDENCE'}}
    monkeypatch.setitem(sys.modules,'multi_epoch',types.SimpleNamespace(run_analysis=run))
    r=client.post('/api/multi-epoch',json={'selection_mode':'replay'})
    assert r.status_code==202
    ident=r.json()['id'];done=await_complete(client,ident)
    assert done['status']=='COMPLETED' and done['selection_mode']=='replay'
    assert client.get(f'/api/multi-epoch/{ident}/readiness').json()==[{'status':'READY'}]
    assert client.get(f'/api/multi-epoch/{ident}/pairs').json()[0]['status']=='REJECTED'
    assert client.get(f'/api/multi-epoch/{ident}/time-series').json()['status']=='INSUFFICIENT_EVIDENCE'
    assert client.get(f'/api/multi-epoch/{ident}/artifacts/report.json').status_code==200
    assert client.get(f'/api/multi-epoch/{ident}/artifacts/job.json').status_code==422
    assert client.get(f'/api/multi-epoch/{ident}/artifacts/nope.fits').status_code==404


def test_busy_has_no_queue_and_failure_is_structured(client,monkeypatch):
    entered=threading.Event();release=threading.Event()
    def run(**kw):
        entered.set();assert release.wait(3)
        raise ValueError('missing calibration')
    monkeypatch.setitem(sys.modules,'multi_epoch',types.SimpleNamespace(run_analysis=run))
    first=client.post('/api/multi-epoch',json={});ident=first.json()['id']
    assert entered.wait(2)
    try:
        assert client.post('/api/multi-epoch',json={}).status_code==429
        assert client.get(f'/api/multi-epoch/{ident}/pairs').status_code==409
    finally:release.set()
    final=await_complete(client,ident)
    assert final['status']=='FAILED' and final['error']['code']=='ValueError'


@pytest.mark.parametrize('body',[{'maximum':9},{'maximum':1},{'selected_ids':['bad','also bad']},{'ra':1},{'mode':'guess'},{'cache_root':'C:/secret'},{'replay_report':'secret.json'}])
def test_invalid_inputs_rejected(client,body):
    assert client.post('/api/multi-epoch',json=body).status_code==422


def test_discovery_declares_replay_and_live_failures(client,monkeypatch,tmp_path):
    record={'obs_publisher_did':'ivo://irsa.ipac/spherex_dr1?x/D1'}
    packet=tmp_path/'packet.json';packet.write_text(json.dumps({'record':record}))
    report=tmp_path/'report.json';report.write_text(json.dumps({'observations':[{'measurement_file':str(packet)}]}))
    monkeypatch.setattr(api,'REPLAY_REPORT',report)
    replay=client.get('/api/multi-epoch/discovery?selection_mode=replay').json()
    assert replay['selection_mode']=='replay' and replay['observations']==[record]
    def fail(**kw):raise api.requests.HTTPError('502 upstream gateway')
    monkeypatch.setattr(api,'discover_pairs',fail)
    live=client.get('/api/multi-epoch/discovery')
    assert live.status_code==502 and live.json()['detail']['code']=='ARCHIVE_UNAVAILABLE'
    from reference_catalog import GaiaArchiveUnavailable
    monkeypatch.setattr(api,'discover_pairs',lambda **kw:(_ for _ in ()).throw(GaiaArchiveUnavailable('Gaia timed out after bounded retries')))
    gaia=client.get('/api/multi-epoch/discovery')
    assert gaia.status_code==503 and gaia.json()['detail']['code']=='GAIA_CATALOG_UNAVAILABLE'


def test_artifact_symlink_and_bad_ids_rejected(client,tmp_path):
    ident=str(uuid.uuid4());directory=api.RUN_ROOT/ident;directory.mkdir(parents=True)
    api._jobs[ident]={'status':'COMPLETED','result':{}}
    assert client.get('/api/multi-epoch/not-a-uuid').status_code==422
    assert client.get('/api/multi-epoch/'+str(uuid.uuid4())).status_code==404
    # Windows path separators must never escape the run directory.
    assert client.get(f'/api/multi-epoch/{ident}/artifacts/..%5Csecret.json').status_code==422
    assert client.get('/multi-epoch').status_code==200


def test_frontend_select_options_are_valid_html():
    from html.parser import HTMLParser
    class Options(HTMLParser):
        def __init__(self):
            super().__init__(); self.current=None; self.options={}
        def handle_starttag(self,tag,attrs):
            a=dict(attrs)
            if tag=='select':
                self.current=a['id']; self.options[self.current]=[]
            elif tag=='option':
                assert self.current is not None
                self.options[self.current].append(a.get('value'))
        def handle_endtag(self,tag):
            if tag=='select': self.current=None
    markup=Path('static/multi_epoch.html').read_text(encoding='utf-8')
    assert '</option value=' not in markup
    parsed=Options();parsed.feed(markup)
    assert parsed.options['target']==['rw-aurigae','barnards-star','epsilon-indi-b','wise-0855','r-draconis','']
    assert parsed.options['source']==['live','replay']
    assert parsed.options['mode']==['best_pair','all_valid_pairs','time_series']
    assert parsed.options['filter']==['all','accepted','rejected']
