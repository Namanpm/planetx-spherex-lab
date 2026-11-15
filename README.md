# SPHEREx data pipeline

A minimal, read-only pipeline for querying NASA/IPAC IRSA spectral image records, retrieving a real FITS cutout, and plotting its sky intensity. It does not yet compare observations or identify moving objects.

## Setup (Windows PowerShell)

```powershell
python -m venv .venv
.\.venv\Scripts\python -m pip install -r requirements.txt
```

## 1. Search IRSA

```powershell
.\.venv\Scripts\python spherex.py --ra 127.69444 --dec -39.17760 --radius 0.01 > search.json
```

Coordinates are ICRS degrees; radius is in degrees, greater than zero and at most 0.1. The documented IRSA example position is used here, not a named challenge target. The result contains the observation ID, image center, start/end MJD, wavelength range in meters, detector, and product URL. Image centers need not equal the searched position: the footprint covers the requested area. `--collection spherex_qr3` searches QR3 instead of the DR1 default. These are wide-survey collections; deep fields are not included.

## 2. API

```powershell
.\.venv\Scripts\python -m uvicorn main:app --host 127.0.0.1 --port 8000
```

Open http://127.0.0.1:8000/docs or request `/api/spherex/search?ra=127.69444&dec=-39.17760&radius=0.01`. Optional `collection` is `spherex_dr1` or `spherex_qr3`. Search results are cached in process for 15 minutes with at most 32 entries; restarting clears the cache. Invalid inputs return 422, upstream timeouts 504, and handled upstream errors 502. The cache is not shared across server processes.

## 3. Download and plot a real FITS cutout

```powershell
.\.venv\Scripts\python visualize.py --ra 127.69444 --dec -39.17760
```

This selects the first returned observation, requests a 0.1-degree cutout, writes `artifacts/sample.fits`, and plots the `IMAGE` extension to `artifacts/sky.png`. These files are overwritten on each run. The plot uses the FITS celestial WCS and surface-brightness units and a 2nd–98th percentile display stretch. Black/colored pixels do not indicate a detection or calibrated source spectrum. The detector wavelength range returned by the archive is not a spectrum at the target position. Flags and variance exist in the FITS but are not applied to this exploratory plot. The first record is not necessarily optimal for the target or wavelength.

Example verified 2026-11-14: DR1 search returned 414 records; first observation `2025W18_1B_0162_1`, image date `2025-04-29T14:14:38.121`, detector D1, archive wavelength range 0.7455–1.1167 micrometers. The downloaded FITS image was 28 × 58 pixels in MJy/sr. This example is not evidence of a moving object.

```powershell
.\.venv\Scripts\python -m pytest -q
```

The test fixture `artifacts/search_fixture.xml` is a saved IRSA response; tests do not require network access. `artifacts/search.json`, `sample.fits`, and `sky.png` are the recorded live-run outputs. IRSA access may change or fail independently of this project.

Sources: [IRSA SPHEREx data access](https://caltech-ipac.github.io/spherex-archive-documentation/spherex-data-access/) and [SPHEREx cutout service](https://irsa.ipac.caltech.edu/data/SPHEREx/docs/cutout_tool.html).

## Phase 4 status: partial

`GET /api/targets` returns five SIMBAD catalog positions and proper motions. Positions are ICRS J2000; proper motions are in mas/year with the RA component multiplied by cos(dec). `GET /api/targets/{id}/observations?radius=0.05` searches the DR1 wide collection at that fixed J2000 position. Proper motion and parallax are not applied, so this is a sky-region search, not a verified current object position. `python live_check.py` repeats the four live coverage checks.

On 2026-11-14 those searches returned 255 Barnard's Star, 362 Epsilon Indi B, 228 RW Aurigae and 451 R Draconis records. These are archive footprint matches, not detections, and their time ranges alone do not establish compatible multi-epoch pairs. These recorded counts did not include WISE.

`fits_pipeline.read_cutout` validates the science, variance, and flag extension shapes and units. Its preliminary valid mask excludes every nonzero flag; it is intentionally conservative and not an implementation of the recommended science-specific flag selection. The existing sample FITS had 54 unflagged finite pixels with nonnegative variance out of 1,624.

**Not yet implemented:** epoch grouping, calibrated wavelength compatibility, reprojection, uncertainty-aware subtraction, comparison endpoint, diagnostics and live end-to-end comparison. Detector-wide wavelength ranges are insufficient for compatibility: SPHEREx linear-variable-filter wavelength varies across the detector. IRSA specifies the separate CWAVE and CBAND calibration layers for science analysis; the cutout WCS-WAVE is intended for visualization. See [IRSA data products](https://caltech-ipac.github.io/spherex-archive-documentation/spherex-data-products/). IMAGE includes zodiacal background, which may change between epochs. No change detection claim can be made from these preliminary searches.

The spectral_calibration module maps a cutout to original detector indices and samples caller-supplied CWAVE/CBAND arrays. It does not discover or authenticate their association with an observation; do not use its output to assert science compatibility without independent verification. The sample contains 1,520 SOURCE-only pixels; the current all-nonzero-flags mask deliberately remains conservative.

## DataLink calibration lookup

Run `.\.venv\Scripts\python live_calibration.py` to query IRSA DataLink for the saved real cutout, cache its associated spectral calibration, and write `artifacts/wavelength_maps.fits` and `artifacts/calibration_provenance.json`. The cache is keyed by URL and validated before reuse. Downloads are streamed, limited to 80 MiB, written atomically, and use connection/read timeouts and HTTP retries. The result checks OBSID, detector, array geometry and units; the DataLink relationship establishes association, while direct UUID equality remains unverified. The maps' VALID extension covers finite positive wavelength and bandwidth only, not science-quality flags.

Live result for 2025W18_1B_0162_1/D1: IRSA returned `spectral_wcs_D1_spx_cal-swcs-v5-2026-191.fits`, SHA-256 `44fb48b236e25a99d42c8986da0064cc58359527e1fc76e182eda52371a1ae23`. The 28 x 58 map has 1,624 valid calibration values; wavelength 0.735262–0.739792 µm, bandwidth 0.015783–0.018853 µm. These are not a source spectrum or a detection. The direct UUID match remains unverified. Source: https://irsa.ipac.caltech.edu/onlinehelp/spherex/spherex/visualization.html

## Automated pair investigation and numerical comparison (current scope)

```powershell
.\.venv\Scripts\python scientific_pipeline.py --target barnards-star --max-pairs 2 --output artifacts/barnard_comparison
.\.venv\Scripts\python scientific_pipeline.py --ra 127.69444 --dec -39.17760 --max-pairs 2 --output artifacts/fixed_comparison
.\.venv\Scripts\python -m pytest -q
```

`pair_discovery.discover_pairs` supports all existing target IDs or fixed coordinates, spherical footprint containment, deduplication, bounded epoch grouping, and preliminary detector/temporal pair ranking. Preset positions use approximate barycentric proper motion at each observation date; parallax and radial velocity are absent. Ranking does not establish spectral compatibility and checks a limited candidate set.

`scientific_pipeline.run` automatically retrieves two observations per candidate, validates associated calibrations, saves calibrated FITS, reprojects B onto A, and records geometric overlap, joint quality fraction, wavelength differences, provenance, dependencies, configuration hash, timings and output bytes. It currently rejects real subtraction: full response equivalence and PSF/background controls are not established. This is a functioning automated investigation path, not a completed validated real comparison engine. The two first Barnard pairs share observation IDs 2025W33_2C_0321_1 and 2026W19_1A_0403_1 on D6 and D3, respectively.

`wcs_geometry.reproject_bilinear` propagates independent-pixel variances using squared weights and requires all contributing pixels valid. It returns geometric footprints separately. Output variances are marginal variances; neighboring resampled pixels are correlated. It is a surface-brightness interpolator, not integrated-flux-conserving mosaicking or PSF matching.

`comparison_engine.compare_aligned` provides numerical subtraction, independent-error variance, spectral/quality masks, optional caller-defined blank-sky plane fitting, and exploratory signed residual pixels. Background fits on resampled inputs are prohibited without covariance. It does not infer motion or classify a residual as a source. `response_compatibility` requires identical caller-validated full response samples; central wavelength and bandwidth alone do not authorize a comparison. `write_outputs` saves FITS/PNG/JSON. Synthetic products are under `artifacts/synthetic_comparison/` and must not be read as real detections.

The optional `read_cutout(..., allow_source=True)` includes SOURCE-only pixels; all other nonzero flags remain excluded. IRSA describes SOURCE-only as potentially nominal for photometry; SOURCE is excluded from blank-sky selection and can still indicate contamination. The previous conservative default remains available. See https://irsa.ipac.caltech.edu/onlinehelp/spherex/spherex/sp.html.

API additions: GET `/api/spherex/pairs?ra=...&dec=...`; POST `/api/comparisons` with JSON `ra`, `dec`, `max_pairs`; GET `/api/comparisons/{id}`. Submission is synchronous and may take minutes. These endpoints expose bounded investigation reports, not scientifically validated discoveries. Original search caching is preserved.

Remaining work: full response calibration ingestion and identifiable response-aware spectral modeling; PSF matching; robust background/photometric consistency for real data; full flag-bit classifications; observer parallax; motion association and catalog crossmatching; covariance-aware source significance; peak-memory and throughput profiling; durable job infrastructure. A pair rejected under this conservative implementation does not prove that no valid SPHEREx comparison exists.


## Extended scientific components and current validation

This implementation remains PARTIAL for scientifically accepted real difference imaging. It now implements the independent components below; none of their synthetic successes authorizes a real SPHEREx subtraction by itself.

- `response_calibration.py`: retrieves the official 436,705,280-byte IRSA response bundle, maps detector pixels including MWIR X reversal to geometric channels, reads measured representative responses and applies the documented middle-pixel wavelength-shift approximation. Negative laboratory samples are preserved. Per-pixel response-shape systematics are not supplied by the release; physical overlap metrics reject signed samples instead of inventing a denoised response.
- `spectral_model.py`: response-integrated linear spectral basis, SVD weighted fits, optional ridge, rank/conditioning rejection, conditional coefficient covariance and held-out predictive covariance. A prior/penalty cannot rescue a rank-deficient measurement. Responses and basis must be scientifically validated by the caller; response uncertainty and model inadequacy are not inferred.
- `psf_processing.py`: extracts the real EPSF zone and centered native sampling; stable cross-convolution on common supplied axes, strict support masking and squared-weight variance. Native independent-pixel background-plane fits and covariance are available. It does not establish a single PSF over a spatially varying field or rotate measured detector PSFs into common sky axes automatically.
- `covariance_processing.py`: small-grid exact linear covariance propagation and GLS background subtraction, including fit/residual covariance. `compare_with_covariance` accepts full image covariance on a common grid; cross-exposure covariance and unknown calibration systematic errors remain excluded. Explicit 1024-output-pixel budget prevents hidden dense quadratic growth.
- `spacecraft_ephemeris.py`: verified JPL Horizons SPHEREx (-163182), SSB ICRF AU vectors at TDB observation time, cached with provenance. `astrometry_model.py` supports rectilinear geometric PM/parallax/RV with supplied observer vectors. These are geometric directions, not aberration/light-time-corrected apparent coordinates, and astrometric covariance is not propagated. Missing catalog RV is explicitly conditional on zero RV.
- `source_detection.py` and `motion_inference.py`: masked positive-excursion components, conditional centroids, ambiguity-rejecting source/catalog associations, three-distinct-epoch weighted motion-vs-static fits and compatible-band constant-flux diagnostics. Conditional centroid covariance is not a calibrated astrometric uncertainty. Catalog matches are positional candidates, not confirmed identities. Real automatic motion inference is not yet integrated.

The real runner now records ePSF provenance, source-free native background/ZODI diagnostics, approximate laboratory-response provenance, Horizons geometry, and response-related rejection reasons. It evaluates more than the original two pairs, writes atomic completed-attempt checkpoints, retries failed attempts, validates cached artifact hashes, and locks its output directory against concurrent writers. Checkpoints include sky center, processing configuration and relevant software hash. Search remains bounded, not exhaustive; no failed operation is converted into absence of sky coverage.

Reproduction:

```powershell
.\.venv\Scripts\python scientific_pipeline.py --target barnards-star --max-pairs 6 --output artifacts/barnard_comparison
.\.venv\Scripts\python inference_validation.py
.\.venv\Scripts\python -m pytest -q
```

The response bundle download is approximately 437 MB and cached at `artifacts/response_calibration/response_bundle.tar`. Live response examples are in `artifacts/response_calibration/live_pixel_response.csv` and `.json`; current live pair reports are in `artifacts/barnard_comparison/report.json`. No full calibrated real difference is asserted. Full-response approximation error, common-sky spatially varying PSF matching, robust joint background/photometric treatment and validated end-to-end motion classification remain open.

Synthetic inference experiment: seed 7821, 200 repetitions with independent Gaussian noise and the correct two-parameter linear spectral basis; first six spectral measurements train and last two are held out. Nominal 95% spectral intervals cover truth in 95.5% of repetitions for each coefficient, mean held-out chi-square is 2.0298 for two measurements. Motion slope interval coverage is 94.5% and 95.0%. See `artifacts/inference_validation/report.json` for measured bias, runtime and peak traced Python allocations (not process RSS). Defined source experiment seed 84712 recovers 12/12 injected Gaussians within 0.3 pixel and yields zero blank excursion components over 12,288 unit-normal pixels. These are conditional synthetic distributions, not real-sky completeness/false-positive rates.

Scientific sources: [IRSA response definitions](https://irsa.ipac.caltech.edu/data/SPHEREx/docs/SPHERExSpectralFiles_202607.pdf), [IRSA spectral calibrations](https://irsa.ipac.caltech.edu/data/SPHEREx/docs/spherex_spectral_calibrations.html), [IRSA ePSF tutorial](https://caltech-ipac.github.io/irsa-tutorials/spherex-psf/), and [JPL Horizons API](https://ssd-api.jpl.nasa.gov/doc/horizons.html).

Final phase evidence: `artifacts/nuclear_phase_report.json`. Full suite: 78 passed, 2 known warnings. Nine real pairs were evaluated across Barnard, RW Aurigae, and R Draconis; zero validated real subtractions. R Draconis completed after retrying archive503 errors and fixing Windows FITS validation locks. Final Barnard evaluation took25.35 seconds with cached science/calibration downloads and recorded101,042,638 output-tree bytes; response bundle storage is separate. These are measured run values, not a throughput or process-RAM benchmark. Observer-relative coordinates are post-download geometric coverage diagnostics, not the initial discovery coordinates. PSF WCS projection is now implemented and tested but remains a finite sampled single-zone approximation. Survey invocation: `python -c "from scientific_pipeline import survey; survey(['rw-aurigae','r-draconis'],pairs_per_target=1)"`.


## Real empirical repeatability experiment

The new workflow performs source-level native ePSF photometry and an empirical color-dependent flux-ratio comparison. It does not attempt an unsupported arbitrary-pixel image subtraction. Run the complete experiment and export with:

```powershell
.\.venv\Scripts\python repeatability_products.py --target rw-aurigae --max-observations 8 --output artifacts/rw_empirical
.\.venv\Scripts\python -m pytest -q
```

For cached experiment outputs only, run `python repeatability_products.py artifacts/rw_empirical/report.json`. POST `/api/spherex/repeatability` accepts target and max_observations; one repeatability worker is allowed per API process, with explicit download/catalog/observation bounds. Existing search routes and cache remain available.

The implementation selects D1 observations from the full discovered observation list, samples the temporal span, and ranks candidate pairs using signed representative laboratory response-shape distance. This preserves negative laboratory samples and does not invent an equivalence threshold. References come from authoritative ESA Gaia DR3 TAP queries, excluding known VARIABLE and duplicated sources; NOT_AVAILABLE does not establish constancy. All catalog neighbors to G=19 are retained for deblending, including sources rejected as references. Fainter/unlisted neighbors remain a limitation.

Native fits use the real local oversampled ePSF at source phase, a simultaneous background plane and neighbor profiles. Target-only identifiability is tested even when nuisance columns are dependent. SOURCE-only pixels are allowed; invalid cores and poor fits are excluded. The flux conversion uses an approximate WCS solid angle, not the official map. These systematics are part of the restricted empirical repeatability experiment, not an absolute calibration claim.

References are selected to have disjoint native pixel stamps in both images. Deterministic source-ID groups reserve training, validation and confirmation. Confirmation references do not influence pair selection, and confirmation is tested once without fallback to another pair. The empirical log-flux-ratio model includes Gaia color, training-only outlier rejection, estimated excess scatter, shared parameter covariance, held-out prediction and spatial residual diagnostics. Numeric policies are recorded in report.json, not presented as instrument calibration constants.

Actual RW-field result: 2025W36_1A_0446_1/D1 (2025-09-02T23:27:24.334) versus 2026W10_1A_0445_3/D1 (2026-03-04T07:58:55.867), using cal-swcs-v5-2026-191. Training17, validation17, confirmation8. Confirmation chi-square4.7977 for8 measurements, approximatep0.779, all8 inside nominal95% predictive intervals. Fractional confirmation residual RMS5.81%; excess logscatter0.06135 with asymptotic profile95% upper0.09730. The scale at BP-RP1.23036 is1.18437 with a fitted color term; it is not a universal photometric gain. Direct UUID equality remains unverified as before.

Eight observations in each of RW and Barnard fields produced51 candidate-pair assessments. Barnard had at most14 common independent screened references and did not meet the fixed three-way confirmation design. The final cached-input runs took78.62s (RW) and82.83s (Barnard). RW output-tree size at run end was117,090,182 bytes including caches and earlier development artifacts; this is not peak RAM or network transfer. Full suite107passed, two known warnings.

Artifacts: artifacts/rw_empirical/report.json, heldout_contrasts.fits, reference_measurements.fits, confirmation_diagnostic.png and native measurement JSONs; artifacts/barnard_empirical/report.json; overall manifest artifacts/empirical_phase_results.json.

This is limited evidence of model consistency for the screened, observed source population. It does not prove stellar constancy, calibrated response equivalence, or applicability to RW Aurigae, brown dwarfs, unusual spectra or arbitrary pixels. The high-S/N log-error and fitted-scatter chi-square calculations are approximate; unknown spatial/systematic covariance and small confirmation sample limit confidence. Training scatter is not a universal instrument error floor.

Sources: https://gea.esac.esa.int/archive/documentation/GDR3/Data_analysis/chap_cu7var/sec_cu7var_intro/ssec_cu7var_dataproducts.html and https://irsa.ipac.caltech.edu/onlinehelp/spherex/spherex/sp.html .


## Prepared observation dictionary optimization

Native photometry now prepares each observation once per processing run and stores it in `observation_cache[publisher_id]`. The full key includes detector/release identity; the same pointing on D1 and D2 cannot collide. Each entry owns detached read-only science/variance/flag/mask arrays, its celestial WCS, detached ePSF rows and detector geometry, plus a per-zone ePSF dictionary. Source fits reuse these objects. The dictionary is local to the run and is released afterwards; it is not an unbounded process-global cache.

`get_observation(path, observation_id, cache)` checks OBSID/detector identity before insertion. A changed cutout path, size or modification timestamp replaces its entry. Failed preparations do not get cached. Returned per-source provenance is copied so editing one fit's metadata cannot alter later fits. `fit_source(path,...)` and `extract_epsf(path,...)` remain supported as compatibility entry points. Persistent measured-flux cache fingerprints now also include PSF/FITS/detector-coordinate helper code.

Measured same eight-image / 26-pair / 2,296-source-attempt workflow: before71.87s; after37.81s and36.84s (median37.33s, observed1.93x end-to-end speedup). Native photometry plus one-time preparation fell from48.04s to11.39s and9.98s. Photometry FITS/WCS preparations fell2,296→8 and ePSF table reads1,565→8;60 selected ePSF zones were memoized across the eight images. All520 accepted source-measurement dictionaries, every rejection list and the final confirmation result were exactly equal to the baseline. No scientific thresholds or fitting model changed. Full suite114passed with two known warnings.

Benchmark evidence is in `artifacts/observation_cache_results.json` and `artifacts/observation_cache_equivalence.json`. The old phase profile remains in `artifacts/profile_full_wall_01/`; optimized profiles are in `artifacts/profile_observation_cache_01/` and `_02/`. These are cached-input timings with live metadata requests; no cold-download or process-RAM comparison is claimed.

Reproduce in a fresh output directory:

```powershell
.\.venv\Scripts\python profile_workflow.py --mode wall --output artifacts/profile_cache_recheck
.\.venv\Scripts\python -m pytest -q
```


## Calibration caching and batched WCS areas

The calibration resolver persists validated DataLink associations for 24 hours, keyed by publisher ID, versioned archive product URL, pipeline version, spectral UUID and processing-history digest. Cache schema changes invalidate associations. TTL expiry rechecks the association; failures are not cached. Calibration products are treated as immutable at their versioned archive URLs. To force a product refresh, remove its cached FITS file; `CalibrationCache.clear()` releases process arrays/counters while preserving disk metadata. Direct science/calibration UUID equality remains explicitly unverified.

Loaded CWAVE/CBAND arrays are read-only and shared across observations in a 128 MiB process LRU, keyed by archive product URL, local file identity/timestamps/size and validator revision. Detector identity is checked on every hit. A lock provides process-local single-flight loading; atomic writes protect shared disk readers. Different processes may independently download the same missing product. Observation-specific cutout mappings are computed separately.

`prepare_pixel_areas(observation, positions)` retains scalar sky-to-pixel coordinates and batches the original five pixel-to-sky transformations and spherical-offset determinant. Each observation owns a locked 4096-entry geometry LRU. Prepared WCS objects are immutable during source fitting; requesting a new preparation batch detects WCS/SIP changes and invalidates old geometry. Invalid positions remain subject to the existing per-source rejection rules. No constant-area approximation, scientific threshold or numerical fitting routine was introduced.

Verified workload: eight RW Aurigae field observations, 2296 source attempts, 520 accepted measurements and 26 candidate pairs. All accepted measurement dictionaries (including flux, propagated variance, flags and fit diagnostics), rejected-source lists, pair evaluations, fitted coefficients and held-out confirmation results are exactly identical to the original baseline. Numeric maximum absolute and relative differences are zero. Pair selection remains 2025W36_1A_0446_1/D1 versus 2026W10_1A_0445_3/D1. This remains a restricted reference-population comparison, not arbitrary-pixel image-subtraction validation.

The complete suite passes: **137 passed**, with the two existing FastAPI/ERFA warnings. New tests cover persistent association hits, version/TTL/schema-record invalidation, changed file/validator identity, detector mismatch, array LRU eviction and concurrency; WCS batching covers SIP distortion, RA wrap, poles, invalid geometry, bounded retention and native-fit equivalence. Gateway retries use at most three equivalent URI encodings of the exact same DataLink publisher ID; malformed responses still fail validation.

### Measured performance

The old total was 36.840 s, including 5.344 s of live SIA discovery. Controlled new runs replay the same eight real archived observation records because SIA returned intermittent gateway502 errors. Therefore the comparable old processing time is **31.496 s**. New raw science/catalog/response inputs remain cached; source measurements are fully recomputed. The cold calibration directory initially did not exist and required eight live DataLink lookups plus one 33,312,960-byte calibration FITS download. The warm run starts a fresh Python process, reuses persistent associations/FITS, and loads/validates arrays once. OS file-cache state is not forcibly cleared.

Successful calibration-cold run: **55.998 s**. Warm runs: **14.048 s and 11.011 s**, corresponding to **2.24–2.86x** faster comparable processing than 31.496 s. This is not a claim that the complete live-discovery workflow always finishes in 11 s. Earlier failed cold/live-discovery attempts are retained and excluded from speed comparisons. Cold wall time is dominated by network transfer and gateway retries.

| Stage (seconds) | Baseline excluding discovery | Calibration-cold | Warm persistent cache |
|---|---:|---:|---:|
| catalog_access | 0.037 | 0.033 | 0.033 |
| json_output | 0.172 | 0.166 | 0.161 |
| science_cache_access | 0.110 | 0.104 | 0.100 |
| observation_preparation | 0.140 | 0.148 | 0.140 |
| calibration_resolution | 14.343 | 44.637 | 0.110 |
| response_bundle_access | 0.242 | 0.234 | 0.233 |
| response_curve_sampling | 0.688 | 0.661 | 0.670 |
| catalog_position_propagation | 2.183 | 1.773 | 1.779 |
| pixel_area_batch_preparation | included in native photometry | 0.918 | 0.915 |
| native_photometry | 9.835 | 4.047 | 3.946 |
| response_pair_ranking | 0.004 | 0.003 | 0.003 |
| empirical_model_fitting | 0.288 | 0.230 | 0.232 |
| spatial_residual_diagnostics | 0.046 | 0.037 | 0.038 |
| output_rendering | 0.574 | 0.422 | 0.389 |
| other_orchestration | 2.833 | 2.585 | 2.261 |
| Total excluding discovery | 31.496 | 55.998 | 11.011 |

Nested timings must not be added to those totals: spherical WCS area geometry fell from 5.436 s (1352 scalar calls) to 0.063 s (eight batches). Full coordinate gathering/batch preparation costs 0.915 s. Numerical source fitting remains unmodified: 0.253 s, versus 0.298 s previously. Calibration association resolution took 19.802 s cold and 0.054 s warm; array access including cold transfer took 24.814 s cold and 0.036 s warm.

Warm counters: eight persistent association hits, zero DataLink lookups/downloads, one array validation/digest, seven shared-array hits; 1581 geometry points prepared and all 1352 requested fit areas hit cache (zero scalar misses). `datalink_requests` counts logical discovery calls; transport adapter retries/URI variants can issue additional HTTP requests.

Peak resident process memory (Windows peak working set, includes NumPy/FITS and imports): **249.38 MiB cold**, **242.19 MiB warm**. The retained calibration arrays occupy 31.75 MiB within the 128 MiB limit. Baseline process RAM was not recorded, so no memory-reduction claim is made.

Evidence: `artifacts/calibration_batch_results.json`, `artifacts/profile_calibration_batch_cold_05/profile_summary.json`, `artifacts/profile_calibration_batch_warm_03/profile_summary.json`; each run includes scientific verification and full source packets. Reproduce from the repository root with fresh output/cache paths:

```powershell
.\.venv\Scripts\python profile_workflow.py --output artifacts/new_calibration_cold --calibration-cache artifacts/new_calibration_cache --frozen-observations
.\.venv\Scripts\python profile_workflow.py --output artifacts/new_calibration_warm --calibration-cache artifacts/new_calibration_cache --frozen-observations
.\.venv\Scripts\python -m pytest -q
```

Omit `--frozen-observations` to include live SIA discovery. Live archive latency and availability remain external dependencies. The benchmark needs the existing `artifacts/rw_empirical` baseline science/catalog files and report; it never substitutes measured-flux caches for real fits.

## Selected observations and multi-epoch analysis

Start the shared backend and minimal interface:

```powershell
.\.venv\Scripts\python -m uvicorn main:app --host 127.0.0.1 --port 8000
```

Open `http://127.0.0.1:8000/multi-epoch`. Select a target or coordinates, choose live discovery or explicitly labeled RW Aurigae replay, select 2–8 observations, then choose Best Pair, All Valid Pairs or Time Series. Cards show actual observation dates and detectors; the worker reports actual stages, and generated outputs are available immediately. Reduced-motion support disables transitions. The same server serves the API and interface; no separate frontend build is needed.

The new engine enumerates every mathematical unordered pair: 1 for 2 observations, 10 for 5, and 28 for 8. Same-day, cross-detector, failed calibration and inadequate-reference comparisons receive explicit rejection records. `VALIDATION_CONSISTENT` is a restricted reference-population consistency result. Independent confirmation is tested once for the deterministically selected pair; it is not claimed for every pair. Original `repeatability_products.py` remains available with its existing time-separated-pair policy.

```powershell
# Explicit replay: reuse historical measured snapshots, inspect cached science FITS,
# verify cached calibration hashes/maps, evaluate pairs without archive requests.
.\.venv\Scripts\python multi_epoch.py --target rw-aurigae --mode all_valid_pairs --maximum 5 --replay-report artifacts/rw_empirical/report.json --output artifacts/my_five_epochs
.\.venv\Scripts\python multi_epoch.py --target rw-aurigae --mode time_series --maximum 5 --replay-report artifacts/rw_empirical/report.json --output artifacts/my_five_timeline

# Live discovery and processing; source failures remain explicit, with bounded retries.
.\.venv\Scripts\python multi_epoch.py --target rw-aurigae --mode best_pair --maximum 5 --cache-root artifacts/rw_empirical --output artifacts/my_live_epochs

# Explicit selection: repeat --selected-id for each archive publisher ID.
# Example ID: ivo://irsa.ipac/spherex_dr1?2025W36_1A_0446_1/D1
.\.venv\Scripts\python -m pytest -q
```

Python entry points are `observation_selection.select_observations`, `empirical_comparison.run_selected`, `multi_epoch.run_analysis`, `source_association.catalog_associations`, `source_association.match_position_candidates`, and `joint_timeseries.analyze_time_series`. The configured selection bound is `observation_selection.MAX_OBSERVATIONS=8`; raising it requires checking corresponding solver/covariance/resource limits before exposing it through the API. Explicit selections preserve detector identity; failed epochs are never replaced automatically.

`run_selected` prepares and measures each observation once, preserving neighbor-aware native ePSF fits, masks, variance propagation and the existing calibration/area caches. Pair evaluation reuses those measurements. Replay reuses a declared historical measurement snapshot; it is not a claim that photometry was recomputed with current software. Reports hash the input report, measurement packets, selected science files, catalog and available calibration products, and record current orchestration software. Historical packets lack an original science-file checksum, so a current replay hash documents the input but cannot independently prove that old measurement bytes correspond to an unchanged original science file. Fresh-photometry regression is a separate check.

The source-association layer preserves logical Gaia IDs and epoch-dependent proper-motion predictions. Reported coordinates are forced catalog positions, not independently measured centroids. Available astrometric uncertainties are propagated with an explicitly diagonal approximation; absent errors remain unknown. Positional matching uses summed tangent-plane covariance and rejects ambiguous candidates, but never calls a positional candidate a moving-object identification. No measured motion trajectory is claimed from forced positions.

Time-series calibration is separate from pair acceptance. It simultaneously fits source log-flux nuisance parameters and epoch intercept/color corrections against a fixed anchor, preserving joint coefficient covariance. Globally stamp-disjoint training, validation and confirmation sources remain independent. At least three distinct epochs, identifiable coefficients, in-domain colors and independent validation/confirmation are required. Shared anchor/reference uncertainty appears in the full covariance. Unknown response/PSF/astrometric systematic errors are not invented. An accepted timeline additionally needs independently justified `common_scale_domain` applicability evidence; representative response similarity or a chain of passing pairs does not provide it. Existing real RW Aurigae packets do not establish that common scale. Synthetic accepted timelines test the stated generator assumptions, not the real instrument.

Outputs include `report.json`, readiness JSON/FITS, raw/rejected/missing source measurement FITS, pair-decision FITS and compatibility PNG, a dated readiness plot, time-series status JSON, and a light-curve plot or explicit insufficient-evidence diagnostic. Genuine accepted synthetic timelines also export calibrated points and covariance FITS. Available selected-pair products are copied to flat artifact paths. `artifact_manifest.json` records artifact byte lengths and SHA256 values. A scientific diagnostic status is distinct from whether the job executed successfully.

API routes:

- `GET /api/multi-epoch/discovery`: target/coordinates and `selection_mode=live|replay`.
- `POST /api/multi-epoch`: selected IDs, mode, target/coordinates, maximum and optional Gaia `source_id`; returns HTTP202 and UUID job ID.
- `GET /api/multi-epoch/{id}`: actual progress and terminal execution status.
- `GET /api/multi-epoch/{id}/results`, `/readiness`, `/pairs`, `/time-series`.
- `GET /api/multi-epoch/{id}/artifacts/{filename}`: guarded flat output files.

One worker is allowed; additional submissions receive HTTP429 instead of entering an unbounded queue. Request bodies cannot supply filesystem paths. API replay is restricted to the configured verified RW snapshot; other targets must use live discovery. The Python/CLI interface can explicitly name a local snapshot and checks its target identity. Live discovery failures are surfaced, never silently replaced by replay.

SIA HTTP retries are finite: four total GET attempts by default, exponential backoff plus jitter, and capped backoff/Retry-After waits of four seconds. Permanent client errors stop immediately. Local HTTP-server tests exercise real 502/503 recovery and exhaustion without contacting NASA. Existing DataLink retrieval retains its separately bounded exact-ID URI variants. Read timeouts and multiple logical archive queries still make an overall live run longer than one request; retry bounds are not a claim of a fixed whole-run deadline.

### Verified multi-epoch results

Run `python validate_multi_epoch.py --output artifacts/<fresh_validation_dir> --remeasure` from the repository root using the venv interpreter. `artifacts/multi_epoch_validation_v1/validation.json` records fresh-photometry timings; `artifacts/multi_epoch_validation_final/validation.json` repeats the complete replay/artifact audit with the final orchestration code. Both report `verification_passed: true`. The original 520 accepted measurements/rejections, all 26 original pair evaluations and selected-pair confirmation are exactly reproduced. Exhaustive eight-image enumeration adds two explicit same-day rejections, making 28 records.

The original September/March result retains **17 training, 17 validation, 8 confirmation** sources and **5.8143208% fractional residual RMS**, with 8/8 inside nominal 95% intervals. The stored `rms_log_ratio=0.06156304693` is a different metric (log-residual RMS); comparing 6.1563% log RMS against 5.8143% fractional RMS is incorrect. `selected_pair_summary` now labels both explicitly.

| Selected images | Pair records | Validation consistent | Rejected | Cold process arrays, fresh fits | Warm arrays, fresh fits |
|---|---:|---:|---:|---:|---:|
| 2 | 1 | 1 | 0 | 7.108 s | 2.865 s |
| 5 | 10 | 6 | 4 | 15.382 s | 6.910 s |
| 8 | 28 | 10 | 18 | 20.230 s | 11.918 s |

These fresh-fit timings recompute source photometry using real cached science/catalog/calibration products and supplied archive records. SIA discovery is excluded. Cold-process runs clear the in-memory array cache; the persistent raw-file caches remain populated. They performed 2/3/2 calibration metadata lookups respectively; all warm runs used zero metadata lookups and zero FITS downloads. They are not fully cold internet-download benchmarks. Cumulative process peak RSS reached **286.49 MiB** across the validation process; isolated per-run RAM was not measured. Major phase timings and cache counters are in `validation.json`.

Snapshot replay with final code took 1.005/1.552/1.903 seconds for 2/5/8 observations, including direct cached FITS/calibration readiness checks, all pair evaluations and generated products. These are snapshot timings, not fresh-photometry times.

The automatic five-observation selection is:

- `2025W36_1A_0095_1/D1`: READY, 55 measurements.
- `2025W36_1A_0606_3/D1`: READY, 101 measurements.
- `2026W10_1A_0445_3/D1`: READY, 70 measurements.
- `2026W10_2B_0137_4/D1`: READY, 63 measurements.
- `2026W11_1A_0645_1/D1`: READY, 74 measurements.

Its six validation-consistent pairs must not be presented as six independently confirmed comparisons. The highest-ranked validation-consistent pair fails its final confirmation (p approximately 0.00903); no fallback pair is tried after that result. Best Pair consequently returns insufficient evidence for this subset. The dedicated five-epoch time series also returns `INSUFFICIENT_EVIDENCE` because globally stamp-disjoint training/validation/confirmation populations are too small. Unknown shared response applicability remains a further limitation. All ten pair IDs, diagnostics and rejection reasons are retained in `replay_5/report.json` and `time_series_5/report.json` under the validation directory.

Browser verification exercised live discovery (226 returned records) and explicit selection of five real replay records through the actual interface, producing a completed asynchronous job with readiness, all ten pair decisions and working artifact links. It found and repaired malformed dropdown options; parsed-HTML regression tests now check every target, source and mode option. The UI test used an explicit five-ID subset different from the automatic evenly spaced five-image selection, so its pair counts by status differ legitimately.
