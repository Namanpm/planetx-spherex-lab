from fastapi import FastAPI, HTTPException, Query
import requests
from spherex import search

app = FastAPI(title='SPHEREx search')

@app.get('/api/spherex/search')
def spherex_search(ra: float = Query(ge=0, lt=360), dec: float = Query(ge=-90, le=90), radius: float = Query(default=0.01, gt=0, le=0.1), collection: str = Query(default='spherex_dr1', pattern='^spherex_(dr1|qr3)$')):
    try:
        return search(ra, dec, radius, collection)
    except requests.Timeout as exc:
        raise HTTPException(status_code=504, detail='IRSA request timed out') from exc
    except (requests.RequestException, ValueError, OSError) as exc:
        raise HTTPException(status_code=502, detail=f'IRSA request failed: {exc}') from exc

from targets import targets

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

from pydantic import BaseModel, Field
from pathlib import Path
import json
import uuid
from pair_discovery import discover_pairs

class ComparisonRequest(BaseModel):
    ra: float = Field(ge=0, lt=360, allow_inf_nan=False)
    dec: float = Field(ge=-90, le=90, allow_inf_nan=False)
    max_pairs: int = Field(default=2, ge=1, le=10)

@app.get('/api/spherex/pairs')
def pairs(ra: float = Query(ge=0,lt=360), dec: float = Query(ge=-90,le=90)):
    try:
        return discover_pairs(ra,dec,limit=100)
    except requests.Timeout as exc:
        raise HTTPException(status_code=504,detail='Archive timeout') from exc
    except (ValueError,requests.RequestException) as exc:
        raise HTTPException(status_code=502,detail=str(exc)) from exc

@app.post('/api/comparisons')
def submit_comparison(body: ComparisonRequest):
    from scientific_pipeline import run
    run_id=str(uuid.uuid4())
    try:
        report=run(ra=body.ra,dec=body.dec,max_pairs=body.max_pairs,output=Path('artifacts/runs')/run_id)
        return {'id':run_id,**report}
    except requests.Timeout as exc:
        raise HTTPException(status_code=504,detail='Archive timeout') from exc
    except (ValueError,requests.RequestException) as exc:
        raise HTTPException(status_code=502,detail=str(exc)) from exc

@app.get('/api/comparisons/{run_id}')
def comparison_status(run_id: str):
    try:
        canonical=str(uuid.UUID(run_id))
    except ValueError as exc:
        raise HTTPException(status_code=422,detail='Invalid run ID') from exc
    path=Path('artifacts/runs')/canonical/'report.json'
    if not path.exists(): raise HTTPException(status_code=404,detail='Run not found')
    return json.loads(path.read_text(encoding='utf-8'))

from threading import BoundedSemaphore
_repeatability_slots=BoundedSemaphore(1)

class RepeatabilityRequest(BaseModel):
    target: str
    max_observations: int = Field(default=8, ge=2, le=12)

@app.post('/api/spherex/repeatability')
def submit_repeatability(body: RepeatabilityRequest):
    if body.target not in {t['id'] for t in targets()}:
        raise HTTPException(status_code=422,detail='Unknown target')
    if not _repeatability_slots.acquire(blocking=False):
        raise HTTPException(status_code=429,detail='Repeatability worker is busy')
    try:
        from empirical_comparison import run
        run_id=str(uuid.uuid4())
        report=run(body.target,Path('artifacts/runs')/run_id,body.max_observations)
        from repeatability_products import export
        export(Path('artifacts/runs')/run_id/'report.json')
        return {'id':run_id,**report}
    except requests.Timeout as exc:
        raise HTTPException(status_code=504,detail='Archive/catalog timeout') from exc
    except requests.RequestException as exc:
        raise HTTPException(status_code=502,detail=str(exc)) from exc
    except ValueError as exc:
        raise HTTPException(status_code=422,detail=str(exc)) from exc
    finally:
        _repeatability_slots.release()

from multi_epoch_api import router as multi_epoch_router
app.include_router(multi_epoch_router)
