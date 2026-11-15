"""FastAPI entry point: thin HTTP routes over the scientific modules."""
import json
import uuid
from pathlib import Path
from threading import BoundedSemaphore

import requests
from fastapi import FastAPI, HTTPException, Query
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field

from pair_discovery import discover_pairs
from spherex import search
from targets import targets

ROOT = Path(__file__).resolve().parent
RUNS = ROOT / 'artifacts' / 'runs'
STATIC = ROOT / 'static'

app = FastAPI(title='SPHEREx Epoch Lab', description='Planet X x SPHEREx: same-sky, same-wavelength temporal analysis of real SPHEREx observations.')


def _upstream_error(exc):
    """Map archive/catalog failures to HTTP; never turn them into empty results."""
    from reference_catalog import GaiaArchiveUnavailable
    if isinstance(exc, requests.Timeout):
        return HTTPException(status_code=504, detail='Archive request timed out')
    if isinstance(exc, GaiaArchiveUnavailable):
        return HTTPException(status_code=503, detail=str(exc))
    return HTTPException(status_code=502, detail=f'Upstream archive request failed: {exc}')


@app.get('/api/spherex/search')
def spherex_search(ra: float = Query(ge=0, lt=360), dec: float = Query(ge=-90, le=90), radius: float = Query(default=0.01, gt=0, le=0.1), collection: str = Query(default='spherex_dr1', pattern='^spherex_(dr1|qr3)$')):
    try:
        return search(ra, dec, radius, collection)
    except requests.Timeout as exc:
        raise HTTPException(status_code=504, detail='IRSA request timed out') from exc
    except (requests.RequestException, ValueError, OSError) as exc:
        raise HTTPException(status_code=502, detail=f'IRSA request failed: {exc}') from exc


@app.get('/api/targets')
def list_targets():
    return targets()


@app.get('/api/targets/{target_id}/observations')
def target_observations(target_id: str, radius: float = Query(default=0.05, gt=0, le=0.1)):
    target = next((t for t in targets() if t['id'] == target_id), None)
    if target is None:
        raise HTTPException(status_code=404, detail='Unknown target')
    if 'ra_deg' not in target:
        raise HTTPException(status_code=422, detail='Catalog position not verified')
    result = spherex_search(target['ra_deg'], target['dec_deg'], radius, 'spherex_dr1')
    return {'target': target, 'warning': 'Search centered on J2000 position; proper motion and parallax not applied.', **result}


class ComparisonRequest(BaseModel):
    ra: float = Field(ge=0, lt=360, allow_inf_nan=False)
    dec: float = Field(ge=-90, le=90, allow_inf_nan=False)
    max_pairs: int = Field(default=2, ge=1, le=10)


@app.get('/api/spherex/pairs')
def pairs(ra: float = Query(ge=0, lt=360), dec: float = Query(ge=-90, le=90)):
    try:
        return discover_pairs(ra, dec, limit=100)
    except (ValueError, requests.RequestException) as exc:
        raise _upstream_error(exc) from exc


@app.post('/api/comparisons')
def submit_comparison(body: ComparisonRequest):
    from scientific_pipeline import run
    run_id = str(uuid.uuid4())
    try:
        report = run(ra=body.ra, dec=body.dec, max_pairs=body.max_pairs, output=RUNS / run_id)
        return {'id': run_id, **report}
    except BlockingIOError as exc:
        raise HTTPException(status_code=409, detail='Output directory is busy') from exc
    except requests.RequestException as exc:
        raise _upstream_error(exc) from exc
    except ValueError as exc:
        status = 409 if 'busy' in str(exc) else 422
        raise HTTPException(status_code=status, detail=str(exc)) from exc


@app.get('/api/comparisons/{run_id}')
def comparison_status(run_id: str):
    try:
        canonical = str(uuid.UUID(run_id))
    except ValueError as exc:
        raise HTTPException(status_code=422, detail='Invalid run ID') from exc
    path = RUNS / canonical / 'report.json'
    if not path.exists():
        raise HTTPException(status_code=404, detail='Run not found')
    return json.loads(path.read_text(encoding='utf-8'))


_repeatability_slots = BoundedSemaphore(1)


class RepeatabilityRequest(BaseModel):
    target: str
    max_observations: int = Field(default=8, ge=2, le=12)


@app.post('/api/spherex/repeatability')
def submit_repeatability(body: RepeatabilityRequest):
    if body.target not in {t['id'] for t in targets()}:
        raise HTTPException(status_code=422, detail='Unknown target')
    if not _repeatability_slots.acquire(blocking=False):
        raise HTTPException(status_code=429, detail='Repeatability worker is busy')
    try:
        from empirical_comparison import run
        from repeatability_products import export
        run_id = str(uuid.uuid4())
        report = run(body.target, RUNS / run_id, body.max_observations)
        export(RUNS / run_id / 'report.json')
        return {'id': run_id, **report}
    except (requests.RequestException, RuntimeError, OSError) as exc:
        raise _upstream_error(exc) from exc
    except ValueError as exc:
        status = 409 if 'busy' in str(exc) else 422
        raise HTTPException(status_code=status, detail=str(exc)) from exc
    finally:
        _repeatability_slots.release()


from multi_epoch_api import router as multi_epoch_router  # noqa: E402
from epoch_lab_api import router as epoch_lab_router  # noqa: E402
app.include_router(multi_epoch_router)
app.include_router(epoch_lab_router)
app.mount('/static', StaticFiles(directory=STATIC), name='static')


@app.get('/', include_in_schema=False)
def index():
    return FileResponse(STATIC / 'index.html', media_type='text/html')
