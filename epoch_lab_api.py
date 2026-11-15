"""HTTP API for the SPHEREx Epoch Lab: target resolution, screening, comparisons, series.

Long work runs on one background worker with a bounded queue; job state is
written atomically to disk so status survives page reloads, and jobs that were
running when the server stopped are reported as INTERRUPTED, never as done.
"""
import json
import threading
import uuid
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Annotated

import requests
from fastapi import APIRouter, HTTPException, Query
from fastapi.responses import FileResponse
from pydantic import BaseModel, ConfigDict, Field, StringConstraints, model_validator

from observation_selection import PID
from targets import targets

router = APIRouter(prefix='/api/lab')
ROOT = Path(__file__).resolve().parent
JOB_ROOT = ROOT / 'artifacts' / 'epoch_lab' / 'jobs'
DEMO_ROOT = ROOT / 'artifacts' / 'epoch_lab' / 'demos'
SESAME_MIRRORS = ('https://cds.unistra.fr/cgi-bin/nph-sesame/-oI/SNV?', 'https://vizier.cfa.harvard.edu/viz-bin/nph-sesame/-oI/SNV?')
MAX_QUEUED = 4
SERVED_SUFFIXES = {'.json', '.png', '.fits', '.csv'}
_executor = ThreadPoolExecutor(max_workers=1, thread_name_prefix='epoch-lab')
_lock = threading.RLock()
_jobs = {}
ObservationId = Annotated[str, StringConstraints(max_length=256)]


class FieldRequest(BaseModel):
    model_config = ConfigDict(extra='forbid')
    target_id: str | None = Field(default=None, max_length=80)
    ra: float | None = Field(default=None, ge=0, lt=360, allow_inf_nan=False)
    dec: float | None = Field(default=None, ge=-90, le=90, allow_inf_nan=False)

    @model_validator(mode='after')
    def one_field(self):
        if self.target_id is not None:
            if self.ra is not None or self.dec is not None:
                raise ValueError('Give a preset target or coordinates, not both')
            if self.target_id not in {t['id'] for t in targets()}:
                raise ValueError('Unknown preset target')
        elif self.ra is None or self.dec is None:
            raise ValueError('Give a preset target or both RA and Dec')
        return self


class CompareRequest(FieldRequest):
    observation_a: ObservationId
    observation_b: ObservationId
    psf_match: bool = True

    @model_validator(mode='after')
    def observations(self):
        for ident in (self.observation_a, self.observation_b):
            if not PID.fullmatch(ident):
                raise ValueError('Invalid observation publisher ID')
        if self.observation_a == self.observation_b:
            raise ValueError('Select two different observations')
        if self.observation_a[-1] != self.observation_b[-1]:
            raise ValueError('Observations must be on the same detector')
        return self


class SeriesRequest(FieldRequest):
    observation_ids: list[ObservationId] = Field(min_length=2, max_length=8)

    @model_validator(mode='after')
    def observations(self):
        if len(set(self.observation_ids)) != len(self.observation_ids):
            raise ValueError('Observation IDs must be distinct')
        if any(not PID.fullmatch(i) for i in self.observation_ids):
            raise ValueError('Invalid observation publisher ID')
        if len({i[-1] for i in self.observation_ids}) != 1:
            raise ValueError('Series observations must share one detector')
        return self


def _job_dir(job_id):
    try:
        canonical = str(uuid.UUID(job_id))
    except ValueError as exc:
        raise HTTPException(422, detail={'code': 'INVALID_JOB_ID'}) from exc
    return JOB_ROOT / canonical


def _save(job_id, state):
    directory = JOB_ROOT / job_id; directory.mkdir(parents=True, exist_ok=True)
    temporary = directory / 'job.tmp'
    temporary.write_text(json.dumps(state, default=str), encoding='utf-8')
    temporary.replace(directory / 'job.json')


def _state(job_id):
    directory = _job_dir(job_id)
    with _lock:
        if directory.name in _jobs:
            return dict(_jobs[directory.name])
    path = directory / 'job.json'
    if not path.is_file():
        raise HTTPException(404, detail={'code': 'JOB_NOT_FOUND'})
    state = json.loads(path.read_text(encoding='utf-8'))
    if state.get('status') in ('QUEUED', 'RUNNING'):
        state.update(status='INTERRUPTED', error={'code': 'WORKER_RESTARTED', 'message': 'The server stopped before this job finished; submit it again.'})
    return state


def _field(body):
    from epoch_lab import field_center
    target = next((t for t in targets() if t['id'] == body.target_id), None) if body.target_id else None
    center = field_center(body.target_id, body.ra, body.dec)
    return target, center


def _error(exc):
    from reference_catalog import GaiaArchiveUnavailable
    if isinstance(exc, GaiaArchiveUnavailable):
        return {'code': 'GAIA_CATALOG_UNAVAILABLE', 'message': str(exc)}
    if isinstance(exc, requests.Timeout):
        return {'code': 'ARCHIVE_TIMEOUT', 'message': 'A NASA/IRSA request timed out after bounded retries; try again later.'}
    if isinstance(exc, requests.RequestException):
        return {'code': 'ARCHIVE_UNAVAILABLE', 'message': f'NASA/IRSA request failed after bounded retries: {str(exc)[:300]}'}
    if isinstance(exc, ValueError):
        return {'code': 'ANALYSIS_REJECTED', 'message': str(exc)[:1000]}
    return {'code': type(exc).__name__, 'message': str(exc)[:1000]}


def _run(job_id, kind, body):
    def progress(event):
        with _lock:
            _jobs[job_id].update(status='RUNNING', progress=event)
            _save(job_id, _jobs[job_id])
    try:
        progress({'stage': 'starting'})
        result = _execute(kind, body, JOB_ROOT / job_id, progress)
        with _lock:
            _jobs[job_id].update(status='COMPLETED', result=result, progress={'stage': 'completed'})
    except Exception as exc:  # reported to the client with a specific code
        with _lock:
            _jobs[job_id].update(status='FAILED', error=_error(exc))
    finally:
        with _lock:
            _save(job_id, _jobs[job_id])


def _execute(kind, body, output, progress):
    from epoch_lab import CALIBRATION_CACHE, field_catalog
    from pair_discovery import discover_pairs
    target, center = _field(body)
    if kind == 'screen':
        from spectral_matching import screen
        result = screen(target_id=body.target_id, ra=None if body.target_id else center[0], dec=None if body.target_id else center[1],
                        cache_dir=CALIBRATION_CACHE, progress=progress)
        (output/'screening.json').write_text(json.dumps(result, default=str), encoding='utf-8')
        return result
    progress({'stage': 'archive_discovery'})
    records = discover_pairs(target_id=body.target_id, ra=None if body.target_id else center[0], dec=None if body.target_id else center[1], limit=1000)['observations']
    progress({'stage': 'catalog'})
    catalog, provenance = field_catalog(*center)
    if kind == 'compare':
        from epoch_lab import run_comparison
        return run_comparison(body.observation_a, body.observation_b, records, target=target, center=center, output=output,
                              catalog=catalog, catalog_provenance=provenance, psf_match=body.psf_match, progress=progress)
    from epoch_series import run_series
    return run_series(body.observation_ids, records, target=target, position=center, center=center, output=output,
                      catalog=catalog, catalog_provenance=provenance, progress=progress)


def _submit(kind, body):
    with _lock:
        active = sum(j['status'] in ('QUEUED', 'RUNNING') for j in _jobs.values())
        if active >= MAX_QUEUED:
            raise HTTPException(429, detail={'code': 'QUEUE_FULL', 'message': f'{active} analyses are already queued or running; try again shortly.'})
        job_id = str(uuid.uuid4())
        _jobs[job_id] = {'id': job_id, 'kind': kind, 'status': 'QUEUED', 'request': body.model_dump(), 'progress': {'stage': 'queued'}}
        for old in [k for k, v in _jobs.items() if v['status'] in ('COMPLETED', 'FAILED')][:-32]:
            del _jobs[old]  # bounded memory; all states remain on disk
        _save(job_id, _jobs[job_id])
    _executor.submit(_run, job_id, kind, body)
    return {'id': job_id, 'status': 'QUEUED', 'status_url': f'/api/lab/jobs/{job_id}'}


@router.get('/targets')
def lab_targets():
    return targets()


@router.get('/resolve')
def resolve(name: Annotated[str, Query(min_length=1, max_length=80)]):
    """Resolve an object name with CDS Sesame (SIMBAD/NED/VizieR); presets match locally first."""
    key = name.strip().lower().replace("'", '').replace(' ', '-')
    for t in targets():
        if key in (t['id'], t['name'].lower().replace("'", '').replace(' ', '-')):
            return {'name': t['name'], 'ra_deg': t['ra_deg'], 'dec_deg': t['dec_deg'], 'target_id': t['id'], 'source': 'preset (SIMBAD)'}
    response = None
    for mirror in SESAME_MIRRORS:  # official CDS service, then its CfA mirror
        try:
            response = requests.get(mirror + requests.utils.quote(name.strip()), timeout=(5, 12), allow_redirects=True,
                                    headers={'User-Agent': 'spherex-epoch-lab'})
            response.raise_for_status()
            break
        except requests.RequestException:
            response = None
    if response is None:
        raise HTTPException(502, detail={'code': 'RESOLVER_UNAVAILABLE', 'message': 'CDS Sesame and its mirror did not respond; enter RA/Dec instead.'})
    for line in response.text.splitlines():
        if line.startswith('%J '):
            parts = line.split()
            try:
                ra, dec = float(parts[1]), float(parts[2])
            except (IndexError, ValueError):
                break
            if 0 <= ra < 360 and -90 <= dec <= 90:
                return {'name': name.strip(), 'ra_deg': ra, 'dec_deg': dec, 'target_id': None, 'source': 'CDS Sesame'}
    raise HTTPException(404, detail={'code': 'NAME_NOT_RESOLVED', 'message': f'Sesame could not resolve "{name.strip()}".'})


@router.post('/screen', status_code=202)
def submit_screen(body: FieldRequest):
    return _submit('screen', body)


@router.post('/compare', status_code=202)
def submit_compare(body: CompareRequest):
    return _submit('compare', body)


@router.post('/series', status_code=202)
def submit_series(body: SeriesRequest):
    return _submit('series', body)


@router.get('/jobs/{job_id}')
def job_status(job_id: str):
    return _state(job_id)


def _serve(directory, name):
    if name != Path(name).name or '\\' in name or name.startswith('job.'):
        raise HTTPException(422, detail={'code': 'INVALID_FILE'})
    path = (directory / name).resolve()
    if path.parent != directory.resolve() or path.suffix not in SERVED_SUFFIXES or not path.is_file():
        raise HTTPException(404, detail={'code': 'FILE_NOT_FOUND'})
    return FileResponse(path)


@router.get('/jobs/{job_id}/files/{name}')
def job_file(job_id: str, name: str):
    state = _state(job_id)
    if state['status'] != 'COMPLETED':
        raise HTTPException(409, detail={'code': 'RESULT_NOT_READY', 'status': state['status']})
    return _serve(_job_dir(job_id), name)


@router.get('/demos')
def demos():
    """Recorded real-data runs (see make_demos.py); each lists its provenance."""
    out = []
    for manifest in sorted(DEMO_ROOT.glob('*/manifest.json')):
        try:
            out.append(json.loads(manifest.read_text(encoding='utf-8')))
        except (OSError, ValueError):
            continue
    return out


@router.get('/demos/{name}/files/{file}')
def demo_file(name: str, file: str):
    if name != Path(name).name or not name.replace('-', '').replace('_', '').isalnum():
        raise HTTPException(422, detail={'code': 'INVALID_DEMO'})
    directory = DEMO_ROOT / name
    if not (directory / 'manifest.json').is_file():
        raise HTTPException(404, detail={'code': 'DEMO_NOT_FOUND'})
    return _serve(directory, file)
