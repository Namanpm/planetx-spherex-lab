# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

Flat, single-directory Python project (no package, no pyproject) that queries NASA/IRSA SPHEREx spectral-image archives, downloads real FITS cutouts, and runs conservative, heavily-validated comparison/photometry experiments. It is a research pipeline, not a discovery tool: the README repeatedly stresses that nothing here establishes a detection, and that rejections are not proof of absence. Preserve that conservative framing when editing code or docs.

## Commands

Windows/PowerShell; use the venv interpreter.

```powershell
python -m venv .venv
.\.venv\Scripts\python -m pip install -r requirements.txt

.\.venv\Scripts\python -m pytest -q                       # full suite (offline; uses saved fixtures)
.\.venv\Scripts\python -m pytest test_pair_discovery.py -q            # one file
.\.venv\Scripts\python -m pytest test_pipeline.py::test_name -q       # one test

.\.venv\Scripts\python -m uvicorn main:app --host 127.0.0.1 --port 8000   # API, docs at /docs

.\.venv\Scripts\python spherex.py --ra 127.69444 --dec -39.17760 --radius 0.01   # IRSA search -> JSON
.\.venv\Scripts\python visualize.py --ra 127.69444 --dec -39.17760                # sample FITS + PNG
.\.venv\Scripts\python live_check.py                                              # live target coverage checks
.\.venv\Scripts\python live_calibration.py                                        # DataLink calibration lookup
.\.venv\Scripts\python scientific_pipeline.py --target barnards-star --max-pairs 2 --output artifacts/barnard_comparison
.\.venv\Scripts\python repeatability_products.py --target rw-aurigae --max-observations 8 --output artifacts/rw_empirical
.\.venv\Scripts\python inference_validation.py                                    # synthetic coverage experiment
.\.venv\Scripts\python profile_workflow.py --mode wall --output artifacts/<fresh_dir>
```

No linter or formatter is configured. Tests need no network; live scripts do and can fail from IRSA 502/503s independent of the code. Live scripts download large files (the response bundle is ~437 MB, cached under `artifacts/response_calibration/`).

## Architecture

Two largely independent scientific paths share the same data-access and calibration layers.

**Data access / calibration layer**
- `spherex.py` (IRSA SIA search, XML -> records), `archive_http.py` / `calibration_resolver.py` (HTTP session with retries, DataLink association, validated + atomically written cache, 128 MiB shared CWAVE/CBAND LRU), `fits_pipeline.read_cutout` (validates IMAGE/VARIANCE/FLAG shapes+units; conservative mask excludes all nonzero flags unless `allow_source=True`), `spectral_calibration.py`, `response_calibration.py` (official response bundle), `spacecraft_ephemeris.py` (JPL Horizons), `reference_catalog.py` (Gaia DR3 TAP), `targets.py` (SIMBAD target presets).
- `main.py` is the FastAPI app; it is just thin routes over the modules above. Run outputs go to `artifacts/runs/<uuid>/`. Repeatability uses a single-slot semaphore (429 if busy); search results are cached in-process for 15 min.

**Path A: pixel-level pair comparison** (`scientific_pipeline.py`) — `pair_discovery` -> download cutouts + calibration -> `wcs_geometry.reproject_bilinear` -> `comparison_engine` / `covariance_processing` / `psf_processing` / `spectral_model`. It intentionally **rejects real subtraction** unless response equivalence and PSF/background controls hold; the reports record rejection reasons. Checkpointed, locks its output dir, validates cached artifact hashes.

**Path B: source-level empirical repeatability** (`empirical_comparison.py`, `empirical_repeatability.py`, `native_photometry.py`, `repeatability_products.py`) — Gaia reference stars, native ePSF photometry, response-shape-ranked pair selection, and a color-dependent log-flux-ratio model with train/validation/confirmation splits (confirmation tested once, no fallback pair). `native_photometry` keeps a per-run `observation_cache` of prepared observations and uses batched WCS pixel-area geometry (`prepare_pixel_areas`); optimizations were verified to give exactly identical measurements, so keep numerical behavior bit-identical when refactoring.

Downstream/experimental: `source_detection.py`, `motion_inference.py`, `astrometry_model.py`, `inference_validation.py`, `benchmark_verification.py`. Real automatic motion inference is not integrated.

## Conventions worth knowing

- Cache keys and checkpoints include processing config and software hashes; changing helper code invalidates measured-flux caches by design.
- Direct science-vs-calibration UUID equality is explicitly unverified throughout; do not claim it.
- `artifacts/` holds saved fixtures (e.g. `search_fixture.xml`) that tests depend on alongside many recorded live-run outputs; don't delete or regenerate it casually.
- Tests are flat `test_*.py` files at the repo root, one per module, using synthetic data.
