"""Independently re-check the recorded Epoch Lab results.

Every number shown on the Recorded runs page is recomputed here from the saved
products (FITS planes, per-star tables, covariances, observation lists) with
plain arithmetic, without calling the pipeline functions that produced them.
A check is PASS or FAIL when it can be tested and UNVERIFIED when the input it
needs (for example a science file not present locally) is missing; UNVERIFIED
is never counted as a pass.

    python verify_results.py                 # all recorded runs, rewrites docs/VERIFICATION.md
    python verify_results.py --json out.json # also write machine-readable results
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits
from scipy import stats

ROOT = Path(__file__).resolve().parent
DEMOS = ROOT / 'artifacts' / 'epoch_lab' / 'demos'
REPORT = ROOT / 'docs' / 'VERIFICATION.md'
CALIBRATION_CACHE = ROOT / 'artifacts' / 'calibration_cache'


def sha256(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def close(a, b, rel=1e-6, abs_=1e-12):
    return a is not None and b is not None and abs(a - b) <= max(abs_, rel * max(abs(a), abs(b)))


def short(ident):
    return str(ident).split('?')[-1]


class Checks:
    def __init__(self, run):
        self.run, self.rows = run, []

    def add(self, name, ok, detail):
        """ok: True/False when tested, None when the input needed is missing."""
        status = 'UNVERIFIED' if ok is None else ('PASS' if ok else 'FAIL')
        self.rows.append({'run': self.run, 'check': name, 'status': status, 'detail': detail})


def file_hashes(c, report):
    for ident, info in (report.get('science_files') or {}).items():
        path = ROOT / Path(str(info.get('path', '')).replace('\\', '/'))
        if not path.is_file():
            c.add('science file SHA-256', None, f'{short(ident)}: file not present locally ({path.name})')
            continue
        actual = sha256(path)
        c.add('science file SHA-256', actual == info['sha256'], f'{short(ident)}: recorded {info["sha256"][:12]}…, recomputed {actual[:12]}…')


def check_comparison(name, directory):
    c = Checks(name)
    r = json.loads((directory / 'report.json').read_text(encoding='utf-8'))
    pol, sp = r['policy'], r['spectral']
    tol = pol['spectral_tolerance_fraction']
    with fits.open(directory / 'comparison.fits') as hdul:
        diff, var, sig = (hdul[k].data.astype(float) for k in ('DIFF', 'VARIANCE', 'SIGNIF'))
        rel, matched, blank_sky = hdul['DLAMBDA'].data.astype(float), hdul['MATCHED'].data.astype(bool), hdul['NULLSKY'].data.astype(bool)

    # Wavelength matching: every compared pixel must satisfy the per-pixel gate.
    if matched.any():
        worst = float(np.nanmax(np.abs(rel[matched])))
        c.add('spectral gate holds on every compared pixel', worst <= tol + 1e-6,
              f'{int(matched.sum())} compared pixels, largest |Δλ|/CBAND = {worst:.4f} (limit {tol})')
    else:
        c.add('spectral gate holds on every compared pixel', None, 'no compared pixels')
    c.add('compared pixel count', int(matched.sum()) == sp['matched_pixels'], f'FITS mask {int(matched.sum())}, report {sp["matched_pixels"]}')
    overlap = np.isfinite(rel)
    left_out = overlap & ~matched
    if overlap.any() and matched.any():
        other_colour = left_out & (np.abs(rel) > tol)
        no_data = left_out & ~other_colour & ~np.isfinite(diff)
        c.add('every pixel left out sampled another wavelength or had no valid data', int(other_colour.sum() + no_data.sum()) == int(left_out.sum()),
              f'{int(left_out.sum())} left out: {int(other_colour.sum())} different wavelength, {int(no_data.sum())} flagged or invalid; '
              f'median |Δλ|/CBAND {float(np.median(np.abs(rel[matched]))):.4f} where compared, {float(np.median(np.abs(rel[overlap]))):.3f} over the whole overlap')

    # Time baseline from the epoch MJDs.
    dt = r['epochs']['b']['mjd'] - r['epochs']['a']['mjd']
    c.add('time baseline', close(dt, r['separation_days']), f'{dt:.4f} d from the MJDs, report {r["separation_days"]:.4f} d')

    # Significance plane = difference / sqrt(variance); the stored variance already carries the null-test inflation.
    null = r.get('null_test') or {}
    use = matched & np.isfinite(diff) & np.isfinite(var) & np.isfinite(sig) & (var > 0)
    if use.any():
        err = float(np.max(np.abs(diff[use] / np.sqrt(var[use]) - sig[use])))
        c.add('significance = difference / √variance', err < 1e-3 * max(1.0, float(np.max(np.abs(sig[use])))), f'{int(use.sum())} pixels, max deviation {err:.1e}')
    else:
        c.add('significance = difference / √variance', None, 'no usable pixels')
    blank = blank_sky & np.isfinite(sig)
    if blank.sum() >= 200:
        spread = 1.4826 * float(np.median(np.abs(sig[blank] - np.median(sig[blank]))))
        scale = null.get('noise_scale_applied') or 1.0
        # Inflate-only: after rescaling, blank sky has unit spread, or less if no inflation was needed.
        ok = 0.9 <= spread <= 1.1 if scale > 1 else spread <= 1.1
        c.add('blank-sky null test: significance spread ≈ 1 after rescaling', ok,
              f'{int(blank.sum())} source-free pixels, robust spread {spread:.3f} (raw {null.get("robust_std_sigma", float("nan")):.3f}, uncertainties × {scale:.3f})')
    else:
        c.add('blank-sky null test', None, f'only {int(blank.sum())} source-free pixels (needs ≥200); the report marks the noise gate {r["gates"].get("noise_model")}')

    # Per-star noise model and significances.
    floor, stars = r.get('systematic_floor') or {}, r.get('measured_stars') or []
    alpha, beta = floor.get('alpha_fractional') or 0.0, floor.get('beta_aperture_mjy_sr_px') or 0.0
    vs = r.get('variability_summary') or {}
    if stars:
        total = np.array([np.sqrt(s['stat_sigma'] ** 2 + beta ** 2 + (alpha * s['scaled_flux_a']) ** 2) for s in stars])
        z = np.array([s['residual'] for s in stars]) / total
        recorded = np.array([s['change_z'] for s in stars])
        c.add('star σ = √(stat² + β² + (αF)²) and z = residual / σ', bool(np.allclose(z, recorded, rtol=1e-6, atol=1e-9)),
              f'{len(stars)} stars, α = {alpha:.4f}, β = {beta:.4f}, max |Δz| = {float(np.max(np.abs(z - recorded))):.1e}')
        frac = np.array([s['residual'] / s['scaled_flux_a'] for s in stars])
        c.add('star ΔF/F = residual / scaled flux in A', bool(np.allclose(frac, [s['fractional_change'] for s in stars], rtol=1e-6, atol=1e-12)), f'{len(stars)} stars')
        reported = ((floor.get('calibration_check') or {}).get('all') or {}).get('robust_std_z')
        if len(stars) >= 30 and reported is not None:
            mad = 1.4826 * float(np.median(np.abs(z - np.median(z))))
            c.add('noise model calibrated (robust spread of z ≈ 1)', 0.8 <= mad <= 1.25 and abs(mad - reported) < 0.02, f'recomputed {mad:.3f}, reported {reported:.3f}')
        else:
            c.add('noise model calibrated (robust spread of z ≈ 1)', None, f'only {len(stars)} stars; the source noise model needs ≥30 and the report marks it UNVERIFIED')
        beyond = int(np.sum(np.abs(z) >= pol['detection_sigma']))
        c.add('stars beyond the detection threshold', beyond == vs.get('beyond_threshold'), f'{beyond} with |z| ≥ {pol["detection_sigma"]}, report {vs.get("beyond_threshold")}; {vs.get("passed_vetting", 0)} passed vetting')
    else:
        c.add('per-star noise model', None, 'no measured stars in this report')

    t = r.get('target') or {}
    if t.get('status') == 'MEASURED':
        z_t = t['fractional_change'] / t['fractional_change_sigma']
        c.add('target significance = ΔF/F ÷ σ', close(z_t, t['change_z']), f'{t["fractional_change"]:.4f} / {t["fractional_change_sigma"]:.4f} = {z_t:.3f}σ, report {t["change_z"]:.3f}σ')
    elif t:
        c.add('unmeasurable target is reported with its reason', bool(t.get('reason')), f'{t.get("status")}: {t.get("reason")}')

    gates = r.get('gates') or {}
    values = set(gates.values())
    expected = 'REJECTED' if 'FAIL' in values else ('PARTIAL' if 'UNVERIFIED' in values else 'VALIDATED')
    c.add('status follows the gates (UNVERIFIED never counts as PASS)', r['status'] == expected,
          ', '.join(f'{k} {v}' for k, v in gates.items()) + f' → {r["status"]}')
    for epoch in r['epochs'].values():
        cal = epoch.get('calibration') or {}
        c.add('own calibration recorded for the epoch', bool(cal.get('calibration_url') and cal.get('calibration_sha256') and cal.get('publisher_id') == epoch['observation_id']),
              f'{short(epoch["observation_id"])} → {str(cal.get("calibration_url", "")).rsplit("/", 1)[-1]} (UUID equality {cal.get("uuid_match")})')
    file_hashes(c, r)
    return c.rows


def check_series(name, directory, targets, target_id):
    c = Checks(name)
    r = json.loads((directory / 'series.json').read_text(encoding='utf-8'))
    ph, ast, tol = r['photometry'], r['astrometry'], r['policy']['tolerance_fraction']
    used = [p for p in ph.get('points', []) if p['status'] == 'USED']
    for p in ph.get('points', []):
        if p['status'] != 'USED':
            c.add('excluded epoch has a stated reason', bool(p.get('reason')), f'{short(p["observation_id"])}: {p.get("reason")}')
    ref = next((e['target'] for e in r['epochs'] if (e.get('target') or {}).get('status') == 'MEASURED'), None)
    if ref and used:
        limit = tol * ref['bandwidth_um']
        worst = max(abs(p['wavelength_um'] - ref['wavelength_um']) for p in used)
        c.add('light-curve points share one wavelength', worst <= limit + 1e-9, f'largest |Δλ| {worst * 1e3:.2f} nm, limit {limit * 1e3:.2f} nm (0.1 × CBAND)')
    test = ph.get('constant_flux_test')
    if test and len(used) >= 2:
        f = np.array([p['calibrated_flux_jy'] for p in used])
        s = np.array([p['flux_sigma_jy'] for p in used])
        w = 1 / s ** 2
        mean = float(np.sum(w * f) / np.sum(w))
        chi2 = float(np.sum(((f - mean) / s) ** 2))
        dof = len(f) - 1
        p = float(stats.chi2.sf(chi2, dof))
        ok = close(mean, test['constant_flux']) and close(chi2, test['chi2_constant'], 1e-5) and dof == test['dof'] and close(p, test['pvalue'], 1e-5)
        c.add('constant-flux χ² and p-value', ok, f'weighted mean {mean * 1e3:.3f} mJy, χ² = {chi2:.3f} for {dof} dof, p = {p:.3f} (report p = {test["pvalue"]:.3f})')
    sol = ast.get('astrometric_solution') or {}
    if sol.get('status') != 'SOLVED':
        c.add('astrometric solution', None, f'not solved: {ast.get("reason") or ast.get("status")}')
    else:
        cov = np.array(sol['covariance'])
        pm = np.array([sol['pm_east_arcsec_per_year'], sol['pm_north_arcsec_per_year']])
        # Parameter order in the 5x5 covariance: offsets (0, 1), proper motion (2, 3), parallax (4).
        c.add('quoted errors are the covariance diagonal', bool(np.allclose(np.sqrt(np.diag(cov)[2:4]), sol['pm_sigma_arcsec_per_year'], rtol=1e-6)) and close(float(np.sqrt(cov[4, 4])), sol['parallax_sigma_arcsec']),
              f'σ(μ) = {np.sqrt(cov[2, 2]):.3f}, {np.sqrt(cov[3, 3]):.3f} ″/yr; σ(ϖ) = {np.sqrt(cov[4, 4]):.3f}″')
        block = cov[2:4, 2:4]
        significance = float(np.sqrt(pm @ np.linalg.solve(block, pm)))
        c.add('proper-motion detection significance', close(significance, sol['pm_significance'], 1e-5), f'√(μᵀC⁻¹μ) = {significance:.2f}σ, report {sol["pm_significance"]:.2f}σ')
        target = next((t for t in targets if t['id'] == target_id), None)
        if not target:
            c.add('comparison with an independent catalog', None, f'target {target_id} not found in targets.py')
        else:
            catalog = np.array([target['pm_ra_cosdec_mas_yr'], target['pm_dec_mas_yr']]) / 1000
            c.add('catalog motion is the independent SIMBAD value', bool(np.allclose(catalog, sol['catalog_pm_arcsec_per_year'], atol=1e-9)), f'{target["name"]}: {catalog[0]:.4f}, {catalog[1]:.4f} ″/yr ({target["source"]})')
            d = pm - catalog
            agree = float(np.sqrt(d @ np.linalg.solve(block, d)))
            c.add('measured motion agrees with SIMBAD (known-mover test)', agree < 3 and significance > 5 and close(agree, sol['pm_difference_from_catalog_sigma'], 1e-5),
                  f'measured ({pm[0]:.2f} ± {np.sqrt(cov[2, 2]):.2f}, {pm[1]:.2f} ± {np.sqrt(cov[3, 3]):.2f}) ″/yr vs SIMBAD ({catalog[0]:.2f}, {catalog[1]:.2f}): {agree:.2f}σ apart, detected at {significance:.1f}σ')
            plx = target.get('parallax_mas')
            if plx:
                dz = abs(sol['parallax_arcsec'] - plx / 1000) / sol['parallax_sigma_arcsec']
                c.add('parallax compared with SIMBAD', close(dz, sol['parallax_difference_from_catalog_sigma'], 1e-4),
                      f'{sol["parallax_arcsec"]:.3f} ± {sol["parallax_sigma_arcsec"]:.3f}″ vs {plx / 1000:.3f}″ ({dz:.1f}σ); Birge ratio {sol["birge_ratio_applied"]:.2f}, so the parallax is reported as not reliably constrained')
    missing = r.get('unavailable_observations')
    c.add('requested and unavailable epochs recorded', missing is not None and 'requested_epoch_count' in r,
          f'{r.get("epoch_count")} of {r.get("requested_epoch_count")} requested epochs calibrated; {len(missing or [])} unavailable')
    c.add('every measured epoch has a science-file hash', set(r.get('science_files') or {}) == {e['observation_id'] for e in r['epochs']},
          f'{len(r.get("science_files") or {})} hashes for {len(r["epochs"])} epochs')
    file_hashes(c, r)
    return c.rows


def cwave_check(s, screened):
    """Read each detector's DataLink CWAVE map ourselves and look the wavelength up at the record's detector pixel."""
    maps, hashed = {}, []
    for detector, info in (s.get('detector_calibrations') or {}).items():
        path = CALIBRATION_CACHE / (hashlib.sha256(info['url'].encode()).hexdigest() + '.fits')
        if path.is_file():
            hashed.append(sha256(path) == info.get('sha256'))
            with fits.open(path) as hdul:
                maps[int(detector)] = hdul['CWAVE'].data.astype(float)
    name = 'wavelength at target = DataLink CWAVE map at the detector pixel'
    if len(maps) < len({o['detector'] for o in screened}):
        return name, None, f'calibration maps for detectors {sorted(set(o["detector"] for o in screened) - set(maps))} are not cached locally'
    bad = 0
    for o in screened:
        x, y = o['detector_xy']
        ix, iy = int(np.floor(x)), int(np.floor(y))
        block = maps[o['detector']][max(iy, 0):iy + 2, max(ix, 0):ix + 2]
        block = block[np.isfinite(block)]
        # Interpolated values must lie within the 2 x 2 neighbourhood that brackets the pixel position.
        bad += not block.size or not (block.min() - 1e-5 <= o['wavelength_um'] <= block.max() + 1e-5)
    return name, bad == 0 and all(hashed), f'{len(screened)} records looked up in {len(maps)} maps (SHA-256 {"matches" if all(hashed) else "DIFFERS"}), {bad} disagree'


def check_screening(name, directory, manifest):
    rows = []
    for file in sorted(directory.glob('screening_*.json')):
        tid = file.stem.replace('screening_', '')
        c = Checks(f'{name} / {tid}')
        s = json.loads(file.read_text(encoding='utf-8'))
        pol = s['policy']
        obs = {o['observation_id']: o for o in s['observations']}
        bad = 0
        for p in s['matched_pairs']:
            a, b = obs[p['observation_a']], obs[p['observation_b']]
            band = min(a['bandwidth_um'], b['bandwidth_um'])
            bad += (a['detector'] != b['detector'] or abs(b['mjd'] - a['mjd']) < pol['min_separation_days']
                    or abs(b['wavelength_um'] - a['wavelength_um']) > pol['tolerance_fraction'] * band)
        c.add('every pair: same detector, ≥1 day apart, |Δλ| ≤ 0.1 × CBAND', bad == 0, f'{len(s["matched_pairs"])} pairs re-derived from the observation list, {bad} break the rule')
        screened = [o for o in s['observations'] if o['status'] == 'SCREENED']
        c.add(*cwave_check(s, screened))
        summary = (manifest.get('summary') or {}).get(tid, {})
        groups3 = sum(g['epoch_count'] >= 3 for g in s['matched_groups'])
        c.add('summary counts match the file', summary.get('matched_pairs') == len(s['matched_pairs']) and summary.get('groups_3plus') == groups3,
              f'{len(s["matched_pairs"])} pairs, {groups3} groups of ≥3 epochs')
        rows += c.rows
    return rows


def verify(demo_root=DEMOS):
    sys.path.insert(0, str(ROOT))
    from epoch_lab import software_sha256
    from targets import targets
    target_list, current = targets(), software_sha256()
    results, runs = [], []
    for manifest_path in sorted(Path(demo_root).glob('*/manifest.json')):
        manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
        directory, name = manifest_path.parent, manifest['name']
        runs.append({k: manifest.get(k) for k in ('name', 'kind', 'title', 'status', 'generated_utc', 'software_sha256', 'reproduce', 'catalog_mode')})
        if manifest.get('status') != 'COMPLETED':
            results.append({'run': name, 'check': 'run completed', 'status': 'FAIL', 'detail': manifest.get('error', '')})
            continue
        if manifest['kind'] == 'compare':
            results += check_comparison(name, directory)
        elif manifest['kind'] == 'series':
            results += check_series(name, directory, target_list, (manifest.get('target') or {}).get('id'))
        else:
            results += check_screening(name, directory, manifest)
        results.append({'run': name, 'check': 'produced by the current processing code', 'status': 'PASS' if manifest.get('software_sha256') == current else 'FAIL',
                        'detail': f'run {str(manifest.get("software_sha256"))[:12]}…, current {current[:12]}…'})
    return {'software_sha256': current, 'runs': runs, 'checks': results}


def markdown(result):
    counts = {k: sum(r['status'] == k for r in result['checks']) for k in ('PASS', 'FAIL', 'UNVERIFIED')}
    lines = ['# Verification of the recorded results', '',
             'Generated by `python verify_results.py`. Each recorded run in `artifacts/epoch_lab/demos/` is re-checked from its saved products '
             '(FITS planes, per-star tables, covariances and observation lists) with independent arithmetic, not by calling the pipeline. '
             'UNVERIFIED means the input needed was not available; it is never counted as a pass.', '',
             f'**{counts["PASS"]} passed, {counts["FAIL"]} failed, {counts["UNVERIFIED"]} unverified.** Processing code SHA-256 `{result["software_sha256"][:16]}…`.', '']
    for run in result['runs']:
        rows = [r for r in result['checks'] if r['run'] == run['name'] or r['run'].startswith(run['name'] + ' / ')]
        lines += [f'## {run["title"]}', '',
                  f'`{run["name"]}` · generated {run["generated_utc"]} · catalog: {run.get("catalog_mode") or "n/a"} · reproduce with `{run["reproduce"]}`', '',
                  '| Check | Result | Detail |', '|---|---|---|']
        for r in rows:
            label = r['check'] if r['run'] == run['name'] else f'{r["run"].split(" / ")[1]}: {r["check"]}'
            lines.append(f'| {label} | {r["status"]} | {r["detail"].replace("|", "/")} |')
        lines.append('')
    return '\n'.join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--json', help='also write the results as JSON to this path')
    parser.add_argument('--no-markdown', action='store_true', help='do not rewrite docs/VERIFICATION.md')
    args = parser.parse_args()
    result = verify()
    if not args.no_markdown:
        REPORT.write_text(markdown(result) + '\n', encoding='utf-8')
    if args.json:
        Path(args.json).write_text(json.dumps(result, indent=2), encoding='utf-8')
    for r in result['checks']:
        print(f'{r["status"]:<10} {r["run"]:<38} {r["check"]}')
    failed = sum(r['status'] == 'FAIL' for r in result['checks'])
    print(f'{len(result["checks"])} checks, {failed} failed, {sum(r["status"] == "UNVERIFIED" for r in result["checks"])} unverified')
    sys.exit(1 if failed else 0)


if __name__ == '__main__':
    main()
