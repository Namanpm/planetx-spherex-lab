# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

SPHEREx Epoch Lab (Planet X × SPHEREx): a Python project (one package, `spherex_lab/`, no pyproject; run everything from the repo root) that compares real NASA/IRSA SPHEREx observations of the same sky **at the same wavelength** across epochs. SPHEREx uses linear variable filters, so two images of one star are only comparable where their per-pixel CWAVE agree; the code screens, gates and validates on that basis. It is a conservative research instrument: residuals are candidates, never discoveries, and unverifiable gates are reported as UNVERIFIED rather than passed. Preserve that framing in code and docs.

## Commands

Windows/PowerShell; use the venv interpreter.

```powershell
python -m venv .venv
.\.venv\Scripts\python -m pip install -r requirements.txt

.\.venv\Scripts\python -m uvicorn main:app --host 127.0.0.1 --port 8000   # UI at /, API docs at /docs
.\.venv\Scripts\python -m pytest -q                                       # full suite (offline)
.\.venv\Scripts\python -m pytest tests/test_difference_imaging.py -q           # one file
.\.venv\Scripts\python -m pytest tests/test_epoch_lab_api.py::test_job_lifecycle_progress_and_files -q   # one test
.\.venv\Scripts\python -m scripts.make_demos                  # regenerate recorded real-data demos (network)
.\.venv\Scripts\python -m scripts.make_demos barnard-pair     # one demo
.\.venv\Scripts\python -m scripts.verify_results              # independent re-check of recorded runs; rewrites docs/VERIFICATION.md
.\.venv\Scripts\python -m scripts.build_static_demo           # server-free copy of recorded runs into site/ (published by .github/workflows/pages.yml)
```

No linter or formatter is configured. Tests need no network; a few real-data tests skip when the local archive caches (`artifacts/science_cache`, `artifacts/calibration_cache`, `artifacts/rw_empirical/*`) are absent. IRSA, Gaia TAP and CDS intermittently fail (502/timeouts); network code retries a bounded number of times and otherwise fails with a specific error. Never substitute data on failure.

## Architecture (see docs/ARCHITECTURE.md, docs/METHODOLOGY.md)

**Epoch Lab (current app)**
- `spherex_lab/lab/spectral_matching.py`: footprint-polygon homography (vertex order = detector corners, validated against `tests/data/footprint_ground_truth.json`) + per-detector DataLink CWAVE/CBAND map → wavelength at target for every archive record; same-wavelength pairs (ranked by predicted matched area) and multi-epoch groups.
- `spherex_lab/lab/difference_imaging.py`: `load_epoch` (real cutout + own calibration) and `compare`: clipped background planes, masked normalised ePSF cross-convolution, bilinear reprojection with matched interpolation smoothing, per-pixel spectral gate, difference pedestal, photometric scale, empirical noise model (alpha/beta), blank-sky null test (inflate-only), injection–recovery, residual classification + vetting. `write_products` writes FITS/PNG/JSON/viewer payload.
- `spherex_lab/lab/epoch_series.py`: ePSF position/flux fits, field-star zero points, same-wavelength light curve, frame solution and PM+parallax astrometric solution (JPL Horizons observer vectors).
- `spherex_lab/lab/epoch_lab.py` orchestrates (catalog snapshot reuse, downloads, provenance hashes); `spherex_lab/lab/epoch_lab_api.py` is the `/api/lab/*` router (validation, one worker + bounded queue, atomic job state, demo serving); `scripts/make_demos.py` regenerates `artifacts/epoch_lab/demos/*`.
- Front end: `static/index.html`, `static/app.js`, `static/app.css` (no build step, no external assets).
- `tests/synthetic_sky.py` builds exact-truth synthetic epochs for tests only.

**Earlier paths (still served, tested)**: `spherex_lab/repeatability/scientific_pipeline.py` (pair investigation that rejects unvalidated subtraction), the Gaia reference-star repeatability experiment (`spherex_lab/repeatability/empirical_comparison.py`, `spherex_lab/imaging/native_photometry.py`, `spherex_lab/repeatability/empirical_repeatability.py`, `spherex_lab/repeatability/repeatability_products.py`) and its multi-epoch API (`multi_epoch*.py`, page `/multi-epoch`). `native_photometry` optimisations were verified bit-identical; keep numerical behaviour unchanged when refactoring it.

**Shared data layer**: `spherex_lab/data/spherex.py` (SIA search; locked in-process cache, optional disk cache enabled by `epoch_lab`), `spherex_lab/data/pair_discovery.py`, `spherex_lab/data/calibration_resolver.py` (DataLink, validated caches), `fits_pipeline.read_cutout` (all non-SOURCE flags invalid with `allow_source=True`), `spherex_lab/imaging/psf_processing.py`, `spherex_lab/imaging/wcs_geometry.py`, `spherex_lab/data/reference_catalog.py` (Gaia DR3), `spherex_lab/data/spacecraft_ephemeris.py` (Horizons), `spherex_lab/data/targets.py`.

## Conventions worth knowing

- Every policy number lives in a module-level `POLICY` dict and is copied into each report; change policies there, not inline.
- Reports carry `software_sha256` of the processing files (LF-normalised, so platform-independent) and SHA-256 of each science FITS; keep `epoch_lab.SOFTWARE_FILES` current when adding processing modules. Editing any of those files (even a comment) makes the recorded demos stale: `tests/test_verify_results.py` fails until you rerun `scripts/make_demos.py` (works offline from the local caches) and `scripts/verify_results.py`.
- Gaia DR3 query results in `artifacts/catalog_cache/*.json` are reused when their circle covers the field (DR3 is fixed); label that mode, never fabricate a catalog.
- Direct science-vs-calibration UUID equality is explicitly unverified throughout; do not claim it.
- `artifacts/` mostly holds ignored caches; committed items are `artifacts/epoch_lab/demos/`, the catalog snapshots, `artifacts/rw_empirical/report.json` and `artifacts/research/`; test fixtures live in `tests/data/`. Don't delete those.
- `docs/RESEARCH_LOG.md` is the verbatim history of earlier phases and their measured results.
- Layout: `spherex_lab/{data,imaging,models,lab,repeatability}` (imports are absolute, e.g. `from spherex_lab.lab import epoch_lab`), `scripts/` run as `python -m scripts.<name>`, tests in `tests/` (`pytest.ini` sets `pythonpath = .`), fixtures in `tests/data/`. `spherex_lab.ROOT` is the repo root; code that hashes sibling modules uses `spherex_lab.source_file(name)`.
