# Data sources

All science data are real and retrieved at run time from public archives. Nothing
is simulated in the archive pipeline; `synthetic_sky.py` is used only by the tests.

| Source | What is used | How | Notes |
|---|---|---|---|
| NASA/IPAC IRSA SPHEREx DR1 (`spherex_dr1`) | Level-2 spectral images: `IMAGE` (MJy/sr), `VARIANCE`, `FLAGS`, `ZODI`, `EPSF` table, alternate detector WCS | SIA search `https://irsa.ipac.caltech.edu/SIA`, IBE cutout service (`?center=…&size=0.3`) | [Data access](https://caltech-ipac.github.io/spherex-archive-documentation/spherex-data-access/), [cutouts](https://irsa.ipac.caltech.edu/data/SPHEREx/docs/cutout_tool.html), [products](https://caltech-ipac.github.io/spherex-archive-documentation/spherex-data-products/) |
| IRSA DataLink | Spectral calibration products `CWAVE`/`CBAND` (2040×2040 per detector) associated with each observation | `https://irsa.ipac.caltech.edu/datalink/links/spherex?ID=<publisher id>` | Association via DataLink `#calibration`; direct science↔calibration UUID equality is recorded as unverified |
| IRSA SPHEREx spectral response bundle | Fiducial channel responses (existing repeatability path) | `SPHERExSpectralFiles_202607.tar` (~437 MB) | Not needed by the Epoch Lab |
| ESA Gaia DR3 | `gaiadr3.gaia_source` positions, proper motions, G, BP−RP, variability and duplicate flags (G ≤ 19, r = 0.2°) | TAP `https://gea.esac.esa.int/tap-server/tap/sync` | Stored query results are reused when their circle covers the field (DR3 is a fixed release) |
| JPL Horizons | Barycentric ICRF vectors of SPHEREx (`-163182`) at observation times (TDB) | `https://ssd.jpl.nasa.gov/api/horizons.api` | Used for annual-parallax factors |
| CDS Sesame (official and CfA mirror) | Object-name resolution | `nph-sesame/-oI/SNV` | Optional; presets and coordinates do not need it |
| SIMBAD | Preset target positions, proper motions, parallaxes (`targets.py`) | Static values with source links | J2000 positions |

## Preset targets

| ID | Object | Why it is interesting here |
|---|---|---|
| `barnards-star` | Barnard's Star | Large proper motion; field has a stored Gaia snapshot; the star itself is saturated |
| `epsilon-indi-b` | ε Indi B (brown-dwarf binary) | Large proper motion and parallax, unsaturated at 1.1 µm: a known-mover test |
| `rw-aurigae` | RW Aurigae (T Tauri star) | Rotated six-month revisits match only in strips; the target core is flagged |
| `r-draconis` | R Draconis (Mira) | Many same-wavelength groups (18 with ≥3 epochs); very bright in the infrared and expected to saturate; no series has completed because IRSA DataLink returned HTTP 502 |
| `wise-0855` | WISE 0855−0714 | Coldest known brown dwarf; long same-pointing revisits |

## Availability during development

On 2026-11-15 IRSA SIA and DataLink intermittently returned HTTP 502, and the ESA
Gaia TAP and CDS TAP services timed out from this machine. The code retries
idempotent requests a bounded number of times and otherwise fails the job with a
specific error. `make_demos.py` was retried until the archive responded. The
recorded demonstrations use only data actually returned by the archives.

Later on 2026-11-15 the recorded runs were regenerated with the current processing
code from the local caches of those same archive files (science cutouts, DataLink
calibration maps and the stored Gaia DR3 queries). Their scientific values were
identical to the first recording. At that time IRSA SIA returned HTTP 502, the
AAVSO VSX API refused scripted access (HTTP 403) and Gaia TAP timed out, so no new
target could be analysed.
