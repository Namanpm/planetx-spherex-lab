"""Bounded asynchronous multi-epoch API and a small same-origin interface."""
import json
import threading
import uuid
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Annotated, Literal

import requests
from fastapi import APIRouter, HTTPException, Query
from fastapi.responses import FileResponse
from pydantic import BaseModel, ConfigDict, Field, StringConstraints, model_validator
from observation_selection import MAX_OBSERVATIONS, PID, validate_coordinates
from pair_discovery import discover_pairs

router = APIRouter()
_ROOT = Path(__file__).resolve().parent
RUN_ROOT = _ROOT / 'artifacts' / 'multi_epoch_runs'
REPLAY_REPORT = _ROOT / 'artifacts' / 'rw_empirical' / 'report.json'
CACHE_ROOT = _ROOT / 'artifacts' / 'rw_empirical'
_slots = threading.BoundedSemaphore(1)
_executor = ThreadPoolExecutor(max_workers=1, thread_name_prefix='multi-epoch')
_lock = threading.RLock()
_jobs = {}


class AnalysisRequest(BaseModel):
    model_config = ConfigDict(extra='forbid')
    target_id: str | None = Field(default='rw-aurigae', max_length=80)
    ra: float | None = Field(default=None, ge=0, lt=360, allow_inf_nan=False)
    dec: float | None = Field(default=None, ge=-90, le=90, allow_inf_nan=False)
    selected_ids: list[Annotated[str, StringConstraints(max_length=256)]] | None = Field(default=None, min_length=2, max_length=MAX_OBSERVATIONS)
    mode: Literal['best_pair', 'all_valid_pairs', 'time_series'] = 'best_pair'
    maximum: int = Field(default=8, ge=2, le=MAX_OBSERVATIONS)
    selection_mode: Literal['live', 'replay'] = 'live'
    source_id: str | None = Field(default=None, max_length=80)

    @model_validator(mode='after')
    def check_selection(self):
        validate_coordinates(self.ra, self.dec)
        if self.selection_mode == 'replay' and (self.target_id != 'rw-aurigae' or self.ra is not None):
            raise ValueError('The configured replay snapshot is the RW Aurigae field; select that target without coordinate overrides')
        if not self.target_id and self.ra is None:
            raise ValueError('Specify a target or sky coordinates')
        if self.selected_ids is not None:
            if len(set(self.selected_ids)) != len(self.selected_ids) or len(self.selected_ids) > self.maximum:
                raise ValueError('Selected observation IDs must be distinct and within maximum')
            if any(not PID.fullmatch(i) for i in self.selected_ids):
                raise ValueError('Invalid observation publisher ID')
        return self


def _directory(run_id):
    try:
        canonical = str(uuid.UUID(run_id))
    except ValueError as exc:
        raise HTTPException(422, detail={'code':'INVALID_RUN_ID'}) from exc
    return RUN_ROOT / canonical


def _save_state(run_id, state):
    path = _directory(run_id); path.mkdir(parents=True, exist_ok=True)
    temporary = path / 'job.tmp'
    temporary.write_text(json.dumps(state, default=str), encoding='utf-8')
    temporary.replace(path / 'job.json')


def _state(run_id):
    directory = _directory(run_id)
    with _lock:
        if run_id in _jobs:
            return dict(_jobs[run_id])
        path = directory / 'job.json'
        if not path.exists():
            raise HTTPException(404, detail={'code':'RUN_NOT_FOUND'})
        state = json.loads(path.read_text(encoding='utf-8'))
        if state['status'] in ('QUEUED','RUNNING'):
            state.update(status='INTERRUPTED', error={'code':'WORKER_RESTARTED','message':'Worker is no longer present; submit a new run.'})
        return state


def _execute(run_id, request):
    def progress(event):
        with _lock:
            _jobs[run_id].update(status='RUNNING', progress=event)
            _save_state(run_id, _jobs[run_id])
    try:
        from multi_epoch import run_analysis
        progress({'stage':'starting'})
        report = run_analysis(**request.model_dump(), output=_directory(run_id),
                              replay_report=REPLAY_REPORT if request.selection_mode == 'replay' else None,
                              cache_root=CACHE_ROOT, progress=progress)
        with _lock:
            _jobs[run_id].update(status='COMPLETED', result=report, progress={'stage':'completed'})
    except Exception as exc:
        from reference_catalog import GaiaArchiveUnavailable
        with _lock:
            if isinstance(exc, GaiaArchiveUnavailable):
                error={'code':'GAIA_CATALOG_UNAVAILABLE','message':str(exc)}
            else:
                error={'code':type(exc).__name__, 'message':str(exc)[:1000]}
            _jobs[run_id].update(status='FAILED', error=error)
    finally:
        try:
            with _lock:
                _save_state(run_id, _jobs[run_id])
        finally:
            _slots.release()


@router.get('/api/multi-epoch/discovery')
def discovery(target_id: str | None = 'rw-aurigae', ra: float | None = Query(default=None,ge=0,lt=360), dec: float | None = Query(default=None,ge=-90,le=90), selection_mode: Literal['live','replay']='live'):
    try:
        validate_coordinates(ra,dec)
        if selection_mode == 'replay':
            if target_id != 'rw-aurigae' or ra is not None:
                raise ValueError('Configured replay snapshot is the RW Aurigae field')
            report = json.loads(REPLAY_REPORT.read_text(encoding='utf-8'))
            observations = [json.loads(Path(o['measurement_file']).read_text(encoding='utf-8'))['record'] for o in report['observations']]
            return {'selection_mode':'replay','observations':observations,'warning':'Saved real archive records; live coverage has not been queried.'}
        result = discover_pairs(target_id=target_id if ra is None else None, ra=ra,dec=dec,limit=1000)
        return {**result,'selection_mode':'live'}
    except requests.RequestException as exc:
        raise HTTPException(502,detail={'code':'ARCHIVE_UNAVAILABLE','message':'NASA archive discovery failed after bounded retries. Try again later; no cached result was substituted.'}) from exc
    except RuntimeError as exc:
        from reference_catalog import GaiaArchiveUnavailable
        if isinstance(exc, GaiaArchiveUnavailable):
            raise HTTPException(503,detail={'code':'GAIA_CATALOG_UNAVAILABLE','message':str(exc)}) from exc
        raise
    except (ValueError,OSError) as exc:
        raise HTTPException(422,detail={'code':'INVALID_DISCOVERY','message':str(exc)}) from exc


@router.post('/api/multi-epoch', status_code=202)
def submit(body: AnalysisRequest):
    if not _slots.acquire(blocking=False):
        raise HTTPException(429,detail={'code':'WORKER_BUSY','message':'One analysis is already running.'})
    run_id = str(uuid.uuid4())
    try:
        with _lock:
            # Retain at most 64 completed states in RAM; all remain on disk.
            for old in list(_jobs):
                if len(_jobs) < 64: break
                if _jobs[old]['status'] in ('COMPLETED','FAILED'):
                    del _jobs[old]
            _jobs[run_id] = {'id':run_id,'status':'QUEUED','selection_mode':body.selection_mode,'progress':{'stage':'queued'}}
            _save_state(run_id,_jobs[run_id])
        _executor.submit(_execute,run_id,body)
    except Exception:
        _slots.release()
        raise
    return {'id':run_id,'status':'QUEUED','status_url':f'/api/multi-epoch/{run_id}'}


@router.get('/api/multi-epoch/{run_id}')
def status(run_id: str):
    return _state(run_id)


def _result(run_id):
    state = _state(run_id)
    if state['status'] != 'COMPLETED':
        raise HTTPException(409,detail={'code':'RESULT_NOT_READY','status':state['status']})
    return state['result']


@router.get('/api/multi-epoch/{run_id}/results')
def results(run_id: str):
    return _result(run_id)


@router.get('/api/multi-epoch/{run_id}/readiness')
def readiness(run_id: str):
    r = _result(run_id)
    return r.get('readiness',r.get('observations',[]))


@router.get('/api/multi-epoch/{run_id}/pairs')
def pairs(run_id: str):
    r = _result(run_id)
    return r.get('pairs',r.get('pair_evaluations',[]))


@router.get('/api/multi-epoch/{run_id}/time-series')
def time_series(run_id: str):
    r = _result(run_id)
    return r.get('time_series',{'status':'INSUFFICIENT_EVIDENCE','reason':'No validated common scale available'})


@router.get('/api/multi-epoch/{run_id}/artifacts/{name}')
def artifact(run_id: str, name: str):
    _result(run_id)
    directory = _directory(run_id).resolve()
    if name != Path(name).name or '\\' in name or name in ('job.json','job.tmp'):
        raise HTTPException(422,detail={'code':'INVALID_ARTIFACT'})
    path = (directory / name).resolve()
    if path.parent != directory or path.suffix not in ('.fits','.json','.png') or not path.is_file():
        raise HTTPException(404,detail={'code':'ARTIFACT_NOT_FOUND'})
    return FileResponse(path)


@router.get('/multi-epoch', include_in_schema=False)
def interface():
    return FileResponse(Path(__file__).parent/'static'/'multi_epoch.html',media_type='text/html')
