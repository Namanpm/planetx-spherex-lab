# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

SPHEREx Epoch Lab (Planet X × SPHEREx): a flat, single-directory Python project (no package, no pyproject) that compares real NASA/IRSA SPHEREx observations of the same sky **at the same wavelength** across epochs. SPHEREx uses linear variable filters, so two images of one star are only comparable where their per-pixel CWAVE agree; the code screens, gates and validates on that basis. It is a conservative research instrument: residuals are candidates, never discoveries, and unverifiable gates are reported as UNVERIFIED rather than passed. Preserve that framing in code and docs.

## Commands

Windows/PowerShell; use the venv interpreter.

```powershell
python -m venv .venv
.\.venv\Scripts\python -m pip install -r requirements.txt

.\.venv\Scripts\python -m uvicorn main:app --host 127.0.0.1 --port 8000   # UI at /, API docs at /docs
.\.venv\Scripts\python -m pytest -q                                       # full suite (offline)
.\.venv\Scripts\python -m pytest test_difference_imaging.py -q           # one file
.\.venv\Scripts\python -m pytest test_epoch_lab_api.py::test_job_lifecycle_progress_and_files -q   # one test
.\.venv\Scripts\python make_demos.py                  # regenerate recorded real-data demos (network)
.\.venv\Scripts\python make_demos.py barnard-pair     # one demo
.\.venv\Scripts\python verify_results.py              # independent re-check of recorded runs; rewrites docs/VERIFICATION.md
.\.venv\Scripts\python build_static_demo.py           # server-free copy of recorded runs into site/ (published by .github/workflows/pages.yml)
```

No linter or formatter is configured. Tests need no network; a few real-data tests skip when the local archive caches (`artifacts/science_cache`, `artifacts/calibration_cache`, `artifacts/rw_empirical/*`) are absent. IRSA, Gaia TAP and CDS intermittently fail (502/timeouts); network code retries a bounded number of times and otherwise fails with a specific error. Never substitute data on failure.

## Architecture (see docs/ARCHITECTURE.md, docs/METHODOLOGY.md)

**Epoch Lab (current app)**
- `spectral_matching.py`: footprint-polygon homography (vertex order = detector corners, validated against `tests_data/footprint_ground_truth.json`) + per-detector DataLink CWAVE/CBAND map → wavelength at target for every archive record; same-wavelength pairs (ranked by predicted matched area) and multi-epoch groups.
- `difference_imaging.py`: `load_epoch` (real cutout + own calibration) and `compare`: clipped background planes, masked normalised ePSF cross-convolution, bilinear reprojection with matched interpolation smoothing, per-pixel spectral gate, difference pedestal, photometric scale, empirical noise model (alpha/beta), blank-sky null test (inflate-only), injection–recovery, residual classification + vetting. `write_products` writes FITS/PNG/JSON/viewer payload.
- `epoch_series.py`: ePSF position/flux fits, field-star zero points, same-wavelength light curve, frame solution and PM+parallax astrometric solution (JPL Horizons observer vectors).
- `epoch_lab.py` orchestrates (catalog snapshot reuse, downloads, provenance hashes); `epoch_lab_api.py` is the `/api/lab/*` router (validation, one worker + bounded queue, atomic job state, demo serving); `make_demos.py` regenerates `artifacts/epoch_lab/demos/*`.
- Front end: `static/index.html`, `static/app.js`, `static/app.css` (no build step, no external assets).
- `synthetic_sky.py` builds exact-truth synthetic epochs for tests only.

**Earlier paths (still served, tested)**: `scientific_pipeline.py` (pair investigation that rejects unvalidated subtraction), the Gaia reference-star repeatability experiment (`empirical_comparison.py`, `native_photometry.py`, `empirical_repeatability.py`, `repeatability_products.py`) and its multi-epoch API (`multi_epoch*.py`, page `/multi-epoch`). `native_photometry` optimisations were verified bit-identical; keep numerical behaviour unchanged when refactoring it.

**Shared data layer**: `spherex.py` (SIA search; locked in-process cache, optional disk cache enabled by `epoch_lab`), `pair_discovery.py`, `calibration_resolver.py` (DataLink, validated caches), `fits_pipeline.read_cutout` (all non-SOURCE flags invalid with `allow_source=True`), `psf_processing.py`, `wcs_geometry.py`, `reference_catalog.py` (Gaia DR3), `spacecraft_ephemeris.py` (Horizons), `targets.py`.

## Conventions worth knowing

- Every policy number lives in a module-level `POLICY` dict and is copied into each report; change policies there, not inline.
- Reports carry `software_sha256` of the processing files (LF-normalised, so platform-independent) and SHA-256 of each science FITS; keep `epoch_lab.SOFTWARE_FILES` current when adding processing modules. Editing any of those files (even a comment) makes the recorded demos stale: `test_verify_results.py` fails until you rerun `make_demos.py` (works offline from the local caches) and `verify_results.py`.
- Gaia DR3 query results in `artifacts/catalog_cache/*.json` are reused when their circle covers the field (DR3 is fixed); label that mode, never fabricate a catalog.
- Direct science-vs-calibration UUID equality is explicitly unverified throughout; do not claim it.
- `artifacts/` mostly holds ignored caches; committed items are the fixtures tests read, `artifacts/epoch_lab/demos/`, and the catalog snapshots. Don't delete those.
- `docs/RESEARCH_LOG.md` is the verbatim history of earlier phases and their measured results.
- Tests are flat `test_*.py` files at the repo root.
