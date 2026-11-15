"""Regenerate the recorded real-data demonstrations served by /api/lab/demos.

Every demo re-runs the production pipeline on live NASA/IRSA archive records
(cutouts and calibrations are cached under artifacts/). Nothing is hand-edited:
the manifest summarises the product written by the code, and failures are
recorded as failures.

    python make_demos.py                 # all demos
    python make_demos.py barnard-pair    # one demo
"""
import argparse
import json
import shutil
import sys
import time
import traceback
from datetime import datetime, timezone
from pathlib import Path

from epoch_lab import CALIBRATION_CACHE, field_catalog, field_center, run_comparison, software_sha256
from targets import targets

ROOT = Path(__file__).resolve().parent
DEMO_ROOT = ROOT / 'artifacts' / 'epoch_lab' / 'demos'
PREFIX = 'ivo://irsa.ipac/spherex_dr1?'

DEMOS = [
    {'name': 'barnard-pair', 'kind': 'compare', 'target_id': 'barnards-star',
     'observations': ['2025W35_2A_0474_1/D1', '2025W40_1A_0441_1/D1'],
     'title': "Barnard's Star field: full-field matched pair at 0.84 um (31 days)",
     'why': 'Same pointing and orientation one month apart: every overlapping pixel samples the same wavelength, with a Gaia DR3 catalog for validation.'},
    {'name': 'eps-indi-b-series', 'kind': 'series', 'target_id': 'epsilon-indi-b',
     'observations': ['2025W17_4B_0027_3/D2', '2025W23_2B_0241_2/D2', '2025W40_2B_0304_1/D2', '2025W41_2B_0173_1/D2', '2026W17_2A_0308_3/D2'],
     'title': 'Epsilon Indi B: five epochs at 1.107 um over one year',
     'why': 'A nearby brown dwarf with large known proper motion; tests same-wavelength photometry and astrometric trajectory on a known mover.'},
    {'name': 'eps-indi-b-year-pair', 'kind': 'compare', 'target_id': 'epsilon-indi-b',
     'observations': ['2025W18_1B_0238_1/D3', '2026W18_1A_0457_2/D3'],
     'title': 'Epsilon Indi B: same pointing one year apart at 2.03 um',
     'why': 'A 364-day full-field match; the brown dwarf moves ~4.7 arcsec, so its residual must be recognised as known motion, not change.'},
    {'name': 'rw-aurigae-strip', 'kind': 'compare', 'target_id': 'rw-aurigae',
     'observations': ['2025W36_1A_0606_3/D1', '2026W10_1A_0445_3/D1'],
     'title': 'RW Aurigae field: six months apart, opposite orientation',
     'why': 'Field rotation between visits leaves only a narrow strip where wavelengths agree; shows why whole-image differencing would be wrong.'},
    {'name': 'screening-five-targets', 'kind': 'screen', 'target_ids': ['rw-aurigae', 'barnards-star', 'r-draconis', 'epsilon-indi-b', 'wise-0855'],
     'title': 'Spectral-epoch screening of all preset targets',
     'why': 'Every archive record screened for the wavelength sampled at the target; lists same-wavelength pairs and multi-epoch groups.'},
]


def _summary(kind, report):
    if kind == 'compare':
        return {'status': report['status'], 'reason': report.get('reason'), 'gates': report.get('gates'),
                'matched_pixels': report['spectral']['matched_pixels'], 'separation_days': report['separation_days'],
                'wavelength_range_um': report['spectral'].get('wavelength_range_um'),
                'candidate_summary': report.get('candidate_summary'), 'variability_summary': report.get('variability_summary'),
                'target_status': (report.get('target') or {}).get('status')}
    if kind == 'series':
        sol = report['astrometry'].get('astrometric_solution', {})
        return {'photometry_status': report['photometry']['status'], 'used_epochs': report['photometry'].get('used_epochs'),
                'astrometry_status': report['astrometry']['status'],
                'proper_motion_arcsec_per_year': [sol.get('pm_east_arcsec_per_year'), sol.get('pm_north_arcsec_per_year')] if sol else None,
                'proper_motion_sigma': sol.get('pm_sigma_arcsec_per_year') if sol else None,
                'catalog_proper_motion': sol.get('catalog_pm_arcsec_per_year') if sol else None}
    return {target: {'observations': r['observation_count'], 'screened': r['screened_count'], 'matched_pairs': len(r['matched_pairs']),
                     'groups_3plus': sum(g['epoch_count'] >= 3 for g in r['matched_groups'])} for target, r in report.items()}


def run_demo(demo):
    from pair_discovery import discover_pairs
    out = DEMO_ROOT / demo['name']
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    started = time.perf_counter()
    manifest = {'name': demo['name'], 'kind': demo['kind'], 'title': demo['title'], 'why': demo['why'],
                'generated_utc': datetime.now(timezone.utc).isoformat(timespec='seconds'), 'software_sha256': software_sha256(),
                'reproduce': f'python make_demos.py {demo["name"]}', 'data_source': 'NASA/IPAC IRSA SPHEREx DR1 (spherex_dr1) via SIA, cutout service and DataLink'}
    try:
        if demo['kind'] == 'screen':
            from spectral_matching import screen
            results = {}
            for tid in demo['target_ids']:
                results[tid] = screen(target_id=tid, cache_dir=CALIBRATION_CACHE)
                (out/f'screening_{tid}.json').write_text(json.dumps(results[tid], default=str), encoding='utf-8')
            manifest['files'] = [f'screening_{t}.json' for t in demo['target_ids']]
            manifest['summary'] = _summary('screen', results)
        else:
            tid = demo['target_id']; center = field_center(tid)
            target = next(t for t in targets() if t['id'] == tid)
            records = discover_pairs(target_id=tid, limit=1000)['observations']
            catalog, provenance = field_catalog(*center)
            ids = [PREFIX + o for o in demo['observations']]
            if demo['kind'] == 'compare':
                report = run_comparison(ids[0], ids[1], records, target=target, center=center, output=out, catalog=catalog, catalog_provenance=provenance)
            else:
                from epoch_series import run_series
                report = run_series(ids, records, target=target, position=center, center=center, output=out, catalog=catalog, catalog_provenance=provenance)
            manifest['files'] = report['files']
            manifest['summary'] = _summary(demo['kind'], report)
            manifest['target'] = {'id': tid, 'name': target['name']}
            manifest['observations'] = ids
            manifest['catalog_mode'] = (provenance or {}).get('mode')
        manifest['status'] = 'COMPLETED'
    except Exception as exc:  # recorded, never hidden
        manifest.update(status='FAILED', error=f'{type(exc).__name__}: {exc}', traceback=traceback.format_exc()[-2000:])
    manifest['elapsed_seconds'] = round(time.perf_counter()-started, 2)
    (out/'manifest.json').write_text(json.dumps(manifest, indent=2, default=str), encoding='utf-8')
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('names', nargs='*', help='demo names (default: all)')
    args = parser.parse_args()
    chosen = [d for d in DEMOS if not args.names or d['name'] in args.names]
    if args.names and len(chosen) != len(args.names):
        parser.error('unknown demo name; choose from ' + ', '.join(d['name'] for d in DEMOS))
    failed = 0
    for demo in chosen:
        manifest = run_demo(demo)
        failed += manifest['status'] != 'COMPLETED'
        print(json.dumps({k: manifest.get(k) for k in ('name', 'status', 'elapsed_seconds', 'summary', 'error')}, default=str))
    sys.exit(1 if failed else 0)


if __name__ == '__main__':
    main()
