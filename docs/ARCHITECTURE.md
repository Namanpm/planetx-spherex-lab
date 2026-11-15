# Architecture

A flat Python package (FastAPI + NumPy/SciPy/Astropy) with a dependency-free
single-page front end. There is no build step: `uvicorn main:app` serves the API,
the UI and the recorded demonstrations.

```
                 ┌─────────────────────── static/ (index.html, app.js, app.css) ─────────────┐
 browser ──HTTP──┤ Explore: field → screening → pair/group → comparison/series results       │
                 └───────────────────────────────┬──────────────────────────────────────────┘
                                                 │ /api/lab/*  (epoch_lab_api.py: validation,
                                                 │  1 worker + bounded queue, atomic job state)
   ┌─────────────────────────────────────────────┴───────────────────────────────────────────┐
   │ epoch_lab.py — orchestration: field centre, Gaia catalog (stored DR3 query or live TAP),  │
   │                cutout download, products, provenance (software SHA-256, file hashes)     │
   ├──────────────────────┬──────────────────────────────┬────────────────────────────────────┤
   │ spectral_matching.py │ difference_imaging.py         │ epoch_series.py                    │
   │ footprint homography │ background, masked PSF        │ ePSF position/flux fits, zero      │
   │ + DataLink CWAVE map │ cross-convolution, reproject, │ points, light curve, frame solve,  │
   │ → λ at target, pairs,│ spectral gate, noise models,  │ PM + parallax solution (Horizons)  │
   │ groups, matched area │ null test, injection, classify│                                    │
   ├──────────────────────┴──────────────────────────────┴────────────────────────────────────┤
   │ existing tested modules reused unchanged in behaviour:                                    │
   │ spherex (SIA search + caches) · pair_discovery · scientific_pipeline.download_cutout      │
   │ calibration_resolver (DataLink, validated caches) · spectral_calibration · fits_pipeline  │
   │ psf_processing · wcs_geometry · source_detection · native_photometry · reference_catalog  │
   │ astrometry_model · spacecraft_ephemeris (JPL Horizons) · motion_inference · targets       │
   └──────────────────────────────────────────────────────────────────────────────────────────┘
```

## Request flow

1. `POST /api/lab/screen {target_id | ra, dec}` → job. `spectral_matching.screen`
   runs SIA discovery (proper motion applied for presets), resolves one DataLink
   calibration map per detector (remembered for 24 h, arrays validated on load),
   and returns per-record wavelengths, ranked pairs and groups.
2. `POST /api/lab/compare {field, observation_a, observation_b}` → job.
   `epoch_lab.run_comparison` downloads both cutouts (validated, atomic cache keyed
   by URL+centre+size), loads each observation's own wavelength maps, calls
   `difference_imaging.compare`, and writes `comparison.fits`, `report.json`,
   `viewer.json` (base64 float32 planes for the canvas viewer), `comparison.png`
   and `candidates.csv`.
3. `POST /api/lab/series {field, observation_ids}` → job. `epoch_series.run_series`
   adds JPL Horizons observer vectors and writes `series.json`, `series.png` and
   `light_curve.csv`.
4. `GET /api/lab/jobs/{id}` returns status and progress; `.../files/{name}` serves
   products (only `.json/.png/.fits/.csv`, no path traversal).
5. `GET /api/lab/demos` lists recorded runs (`make_demos.py`); their files are
   served read-only.

Older routes are unchanged: `/api/spherex/search`, `/api/targets`,
`/api/spherex/pairs`, `/api/comparisons`, `/api/spherex/repeatability`,
`/api/multi-epoch/*` and the reference-star repeatability page at `/multi-epoch`.

## Caches and invalidation

| Cache | Location | Key / invalidation |
|---|---|---|
| SIA responses | process (15 min, 32 entries, locked) + `artifacts/sia_cache` (24 h, Epoch Lab only) | query parameters; an error never falls back to a stale entry |
| Science cutouts | `artifacts/science_cache` | SHA-256 of URL + centre + size; FITS validated on every reuse |
| Calibration maps | `artifacts/calibration_cache` + 128 MiB process LRU | versioned archive URL + file identity; detector checked on hit |
| Screening associations | `artifacts/calibration_cache/screening_associations.json` | 24 h TTL; trusted-URL check |
| Gaia DR3 | `artifacts/catalog_cache/*.json` | ADQL circle must contain the requested field |
| JPL Horizons | `artifacts/observer_cache` | request parameters; response hash recorded |

Every product records the SHA-256 of the processing source files
(`epoch_lab.software_sha256`) and of each science FITS file used.

## Concurrency and failure

* One background worker and at most four queued or running jobs (HTTP 429 beyond that).
  The repeatability and multi-epoch workers keep their own single-slot semaphores.
* Job state is written atomically (`job.tmp` → `job.json`). A job left `RUNNING` by
  a restarted server is reported `INTERRUPTED`, never completed.
* Upstream failures map to specific codes (`ARCHIVE_UNAVAILABLE`, `ARCHIVE_TIMEOUT`,
  `GAIA_CATALOG_UNAVAILABLE`, `ANALYSIS_REJECTED`) and HTTP statuses on the legacy
  routes (502/503/504, 409 for a busy output directory, 422 for invalid input).

## Front end

`static/app.js` is plain ES2020 with no dependencies. Charts are hand-drawn SVG
(coverage, light curve, trajectory) with hover tooltips and table views. The image
viewer paints each plane once into an offscreen canvas and composites
side-by-side, overlay, blink, difference, significance and spectral-match views,
with zoom, pan, a pixel readout (RA/Dec, λ, A, B, D, σ) and residual markers.
Colours: categorical blue/orange validated for colour-vision deficiency on the
dark surface; diverging blue–grey–red for signed images; status colours only for
states.
