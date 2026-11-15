# Scientific methodology

This document describes what the SPHEREx Epoch Lab computes, why, and where its
conclusions stop. Every number quoted as a *policy* is an analysis choice stored
in each report (`policy` field), not a mission calibration constant.

## 1. The comparability problem

SPHEREx disperses light with linear variable filters (LVFs). Within one detector
the central wavelength (`CWAVE`) changes continuously across the array, so the
wavelength sampled at a sky position depends on where that position lands on the
detector. Two observations of the same star are therefore *different-colour*
measurements unless the star fell on pixels with the same `CWAVE` both times. The
archive search result only gives a detector-wide range, which cannot decide this.

A pixel-by-pixel difference of two arbitrary SPHEREx images mixes flux change with
colour, and a red star would look "variable". The lab compares only what was
measured at the same wavelength, and rejects the rest instead of correcting it
with an assumed spectrum.

## 2. Screening every archive record (`spectral_matching.py`)

For each DR1 record whose footprint contains the target (via `pair_discovery`,
with proper motion applied for preset targets):

1. The four vertices of the SIA footprint polygon (`s_region`) correspond to
   detector corners (0,0), (2039,0), (2039,2039), (0,2039). A gnomonic
   homography through them maps the target to an estimated detector pixel.
2. The detector's DataLink-associated `CWAVE`/`CBAND` map gives the wavelength and
   bandwidth at that pixel.

The homography ignores optical distortion. Against 26 real downloaded cutouts
(`tests_data/footprint_ground_truth.json`) its error is median 3.1 px and maximum
5.7 px, i.e. a few per cent of a bandwidth on D1. Screening only nominates
candidates; every comparison re-reads the wavelength maps from the downloaded
cutouts' own calibration associations.

**Pairs** require the same detector, observations at least 1 day apart, and
|Δλ| ≤ 0.10 × CBAND at the target. For each pair the predicted fraction of a
0.3° field that is matched is computed by evaluating both footprint models on a
21×21 grid; pairs are ranked by that area. **Groups** chain three or more
distinct epochs at one wavelength (light-curve and trajectory candidates).

## 3. Matched-wavelength difference imaging (`difference_imaging.py`)

Inputs: two real cutouts (0.3°) with `IMAGE`, `VARIANCE`, `FLAGS`, `ZODI`, `EPSF`,
and each observation's own `CWAVE`/`CBAND` sampled onto the cutout grid.

1. **Masks.** Pixels with any non-zero flag other than SOURCE (bit 21) are invalid
   (`fits_pipeline.read_cutout(allow_source=True)`).
2. **Background.** Per epoch, a weighted plane is fitted to blank sky (flags = 0,
   valid calibration, Gaia sources masked within 2.5 px) with iterative 3σ clipping.
   Its coefficient covariance is propagated into the residual variance.
3. **PSF matching.** The measured ePSF of each epoch (zone nearest the target) is
   sampled at native pixel centres, projected onto the other epoch's pixel axes on
   the sky (`psf_processing.project_kernel`), and the images are cross-convolved
   (A ∗ PSF_B, B ∗ PSF_A). The convolution is *masked and normalised*: an output
   pixel is kept when ≥95% of the kernel's absolute weight is valid, with
   variance Σk²v/(Σk)². Exact-support convolution would discard most of a real
   cutout, because isolated flagged pixels are common.
4. **Registration.** B is reprojected onto A's grid by bilinear interpolation with
   squared-weight variance (`wcs_geometry.reproject_bilinear`). Interpolation
   smooths B by f(1−f) px² per axis; A receives a symmetric 3-tap kernel with the
   same second moment so the two effective PSFs match.
5. **Spectral gate.** B's wavelength map is reprojected the same way. Only pixels
   with |λ_B − λ_A| ≤ 0.10 × CBAND_A are compared ("matched region").
6. **Difference background.** A clipped plane fitted to B − A on source-free
   interior pixels removes the residual pedestal left by the two independent sky
   fits.
7. **Photometric scale.** The median B/A ratio of 5×5-pixel apertures on Gaia stars
   with S/N ≥ 100 (S/N ≥ 20 if fewer than five) gives a scale with standard error
   1.2533·MAD/√N. It is applied only if it differs from unity by more than three
   standard errors, and its variance is then propagated.
8. **Difference and uncertainty.** D = B − sA with
   Var(D) = Var_B + s²Var_A + Var(s)·A² + (β_λ · |Δλ|/λ · max(|A|,|B|))², where the
   last term bounds the residual spectral mismatch for |d ln F/d ln λ| ≤ 4.

**Calibration sources.** The scale, registration and noise-model gates use Gaia
DR3 stars (proper motion propagated to each epoch) when a catalog covers the field.
Without one, bright isolated sources detected in epoch A at ≥20σ are used. They are
assumed mostly constant and stationary. They calibrate the gates but never identify
a residual, and all ≥5σ detections are then masked from the null region.

### Validation gates (all computed from the pair's own data)

| Gate | Test | Pass condition |
|---|---|---|
| Spectral overlap | matched valid pixels | ≥ 200 |
| Registration | median Gaia-star centroid offset after resampling (S/N ≥ 50) | ≤ 0.3 px, ≥ 3 stars |
| Photometric scale | high-S/N star ratio | ≥ 5 stars |
| Pixel noise model | robust σ of D/σ on ≥ 200 source-free interior pixels | ≤ 3 (variances inflated by this factor if > 1, never deflated) |
| Source noise model | robust σ of star z-scores under the empirical model | 0.7–1.4 with ≥ 30 stars |

`UNVERIFIED` means the data could not test a gate. It is not a pass, and the
overall status is then `PARTIAL`. Any `FAIL` makes the comparison `REJECTED`.

### Empirical source noise model

Catalog stars in the matched region (mostly constant) calibrate the source-level
error: Var(ΔF) = σ_stat² + β² + (αF)². α (multiplicative) is the robust excess
scatter of S/N ≥ 100 stars beyond their propagated errors; β (additive, aperture
units) is the same for fainter stars given α. The model is checked per S/N bin
(robust σ of z ≈ 1).

### Injection–recovery

Synthetic point sources with the effective PSF, scaled to aperture S/N = 10, are
added to the real difference image at up to 40 interior positions away from
catalog sources (fixed seed). The production detector is rerun and the recovered
fraction is reported as completeness.

### Residual classification

Connected components with |S| ≥ 5σ and ≥ 2 pixels are explained in this order:
edge of the valid region → association with a Gaia source (proper motion
propagated to each epoch) and **known proper motion** (≥ 0.3 px between epochs) →
**dipoles** (an opposite-sign component within 2.5 px, or an aperture sum with the
opposite sign to the peak) → catalog star change measured with the noise model
(**candidate variable** if |z| ≥ 5, otherwise within the systematic floor). A
residual with no Gaia counterpart is measured in both epochs: present in both
epochs means a **persistent non-Gaia source** (tested for change); present in
exactly one means a **transient candidate**; neither means a low-significance
residual.

Source-level variability candidates must also pass **vetting**: no flagged pixel
within 4 px, at least half an aperture plus 2 px from the region boundary, no
second catalog source in the aperture, a residual explained by a scaled copy of
the source (R² ≥ 0.5), and a noise model calibrated on ≥ 30 stars.

## 4. Same-wavelength time series (`epoch_series.py`)

For a group of ≥ 2 epochs on one detector at one wavelength:

* **ePSF fitting.** The target and image-detected field stars (≥ 15σ, isolated) are
  fitted on native pixels with the measured oversampled ePSF plus a local plane,
  by grid search over sub-pixel position (0.25 px, then 0.05 px). The position
  covariance comes from the χ² curvature, inflated by the reduced χ² when > 1.
* **Light curve.** Each epoch's target flux is divided by a zero point: the median
  flux ratio of field stars whose own wavelengths match the reference epoch's.
  With fewer than 8 such stars the full robust star scatter is used as the
  zero-point error. Epochs whose target wavelength differs by more than the
  tolerance are reported and excluded. A χ² test against a constant flux gives
  `VARIABILITY_CANDIDATE` at p < 10⁻⁶.
* **Astrometry.** Field-star positions give per-epoch frame offsets (alternating
  robust medians, reference epoch fixed) and an empirical astrometric floor. The
  target trajectory is compared with the catalog prediction
  (`astrometry_model.propagate_direction` with JPL Horizons SPHEREx vectors) and
  fitted with a five-parameter solution (offset, proper motion, parallax) whose
  covariance is inflated by the Birge ratio √(χ²/dof) when it exceeds 1.

## 5. What the lab will not conclude

* A residual, a candidate or a χ² excess is not a discovery. All candidates need
  independent follow-up.
* Moving-object claims need a measured multi-epoch trajectory. Catalog-forced
  positions are not detections, and the trajectory analysis is presented as a
  test on known movers.
* Different-wavelength epochs are rejected, not transformed with an assumed spectrum.
* Saturated or flagged targets are reported as unmeasurable (with the flag bits).
* External failures (IRSA, Gaia, JPL, CDS) fail the job with that reason. Nothing
  is substituted, except that a stored Gaia DR3 query whose circle covers the
  field is reused, because DR3 is a fixed release. This is labelled in the report.

## 6. Known limitations

* Single ePSF zone per comparison (the zone nearest the target); PSF variation
  across the 0.3° cutout is absorbed by the empirical α/β floor.
* Marginal (diagonal) variances only after convolution and resampling; the null
  test calibrates the pixel noise scale empirically.
* The flux conversion to Jy uses the local WCS pixel solid angle, not the official
  solid-angle map.
* Rotated revisits (e.g. six months apart) are usually matched only in a strip.
* Bright targets (Barnard's Star, RW Aurigae) are saturated or flagged in their
  cores and cannot be measured, although their fields can.
