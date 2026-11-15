"""SPHEREx Epoch Lab orchestration: screening -> matched-pair comparison -> products.

Every step uses real archive records. Network failures propagate as errors;
a missing catalog degrades association (recorded) but never fabricates one.
"""
import hashlib
import json
import time
from pathlib import Path

import numpy as np

from spherex_lab import ROOT  # noqa: E402
ARTIFACTS = ROOT / 'artifacts'
SCIENCE_CACHE = ARTIFACTS / 'science_cache'
CALIBRATION_CACHE = ARTIFACTS / 'calibration_cache'
CATALOG_CACHE = ARTIFACTS / 'catalog_cache'
CUTOUT_SIZE_DEG = 0.3

from spherex_lab.data import spherex  # noqa: E402
spherex.DISK_CACHE_DIR = str(ARTIFACTS / 'sia_cache')
CATALOG_RADIUS_DEG = 0.2
SOFTWARE_FILES = tuple(f'spherex_lab/{name}.py' for name in (
    'lab/epoch_lab', 'lab/difference_imaging', 'lab/spectral_matching', 'lab/epoch_series', 'imaging/psf_processing',
    'imaging/wcs_geometry', 'imaging/fits_pipeline', 'data/calibration_resolver', 'imaging/spectral_calibration',
    'imaging/source_detection'))  # repository-relative


def science_files(epochs, index):
    """Provenance of each science cutout used: cache path, SHA-256 of its bytes, archive URL."""
    return {e['id']: {'path': str(Path(e['path']).relative_to(ROOT)) if Path(e['path']).is_relative_to(ROOT) else e['path'],
                      'sha256': hashlib.sha256(Path(e['path']).read_bytes()).hexdigest(),
                      'access_url': index[e['id']].get('access_url')} for e in epochs}


def software_sha256():
    digest = hashlib.sha256()
    for name in SOFTWARE_FILES:
        path = ROOT / name
        if path.is_file():
            # Hash LF-normalised text so a Windows (CRLF) and a Linux checkout of the same commit agree.
            digest.update(path.read_bytes().replace(b'\r\n', b'\n'))
    return digest.hexdigest()


def field_center(target_id=None, ra=None, dec=None):
    """Search/cutout centre. Preset targets use the proper-motion position at 2025-11-01."""
    if target_id is not None:
        from spherex_lab.data.targets import barycentric_position
        position = barycentric_position(target_id, '2025-11-01')
        return float(position['ra_deg']), float(position['dec_deg'])
    if ra is None or dec is None or not np.isfinite([ra, dec]).all() or not 0 <= ra < 360 or not -90 <= dec <= 90:
        raise ValueError('Supply a preset target or finite RA [0,360) and Dec [-90,90]')
    return float(ra), float(dec)


def _angular_distance_deg(ra1, dec1, ra2, dec2):
    a = np.deg2rad([ra1, dec1, ra2, dec2])
    h = np.sin((a[3]-a[1])/2)**2 + np.cos(a[1])*np.cos(a[3])*np.sin((a[2]-a[0])/2)**2
    return float(np.rad2deg(2*np.arcsin(np.sqrt(min(max(h, 0), 1)))))


def _cached_catalog(ra, dec, needed_radius):
    import re
    for path in sorted(CATALOG_CACHE.glob('*.json')):
        try:
            snapshot = json.loads(path.read_text(encoding='utf-8'))
            match = re.search(r"CIRCLE\('ICRS',([-\d.]+),([-\d.]+),([-\d.]+)\)", snapshot['provenance']['adql'])
        except (OSError, ValueError, KeyError, TypeError):
            continue
        if not match:
            continue
        cra, cdec, cradius = map(float, match.groups())
        if _angular_distance_deg(ra, dec, cra, cdec) + needed_radius <= cradius + 1e-9:
            return snapshot, path.name
    return None, None


def field_catalog(ra, dec, *, radius=CATALOG_RADIUS_DEG, needed_radius=CUTOUT_SIZE_DEG*0.6, allow_network=True):
    """Gaia DR3 sources covering the field. Returns (sources or None, provenance).

    Gaia DR3 is a fixed release, so a stored query result whose circle contains
    the requested field is used first (labelled with its retrieval time). Otherwise
    the live ESA TAP service is queried and the result stored. Failure leaves the
    catalog unavailable; it is never replaced by another field's sources.
    """
    from spherex_lab.data.reference_catalog import GaiaArchiveUnavailable, query_references
    CATALOG_CACHE.mkdir(parents=True, exist_ok=True)
    snapshot, name = _cached_catalog(ra, dec, needed_radius)
    if snapshot is not None:
        return snapshot['all_sources'], {**snapshot['provenance'], 'mode': 'stored Gaia DR3 query result', 'snapshot_file': name}
    errors = []
    if allow_network:
        try:
            catalog = query_references(ra, dec, radius, cache_dir=CATALOG_CACHE/'queries', limit=5000, max_g_mag=19.)
            path = CATALOG_CACHE/(hashlib.sha256(catalog['provenance']['adql'].encode()).hexdigest()+'.json')
            path.write_text(json.dumps({'all_sources': catalog['all_sources'], 'provenance': catalog['provenance']}), encoding='utf-8')
            return catalog['all_sources'], {**catalog['provenance'], 'mode': 'live Gaia DR3 TAP query', 'snapshot_file': path.name}
        except (GaiaArchiveUnavailable, OSError, ValueError) as exc:
            errors.append(f'{type(exc).__name__}: {exc}')
    return None, {'mode': 'UNAVAILABLE', 'errors': errors or ['live query disabled'],
                  'consequence': 'Residuals cannot be associated with catalog sources; classification is limited.'}


def _record_index(records):
    return {r['obs_publisher_did']: r for r in records if r.get('obs_publisher_did')}


def run_comparison(observation_a, observation_b, records, *, target=None, center, output, catalog=None,
                   catalog_provenance=None, psf_match=True, progress=None):
    """Download (or reuse) both real cutouts, compare, and write products."""
    from spherex_lab.lab.difference_imaging import compare, load_epoch, write_products
    from spherex_lab.repeatability.scientific_pipeline import download_cutout

    def emit(stage, **fields):
        if progress:
            progress({'stage': stage, **fields})
    started = time.perf_counter(); timings = {}
    index = _record_index(records)
    missing = [i for i in (observation_a, observation_b) if i not in index]
    if missing:
        raise ValueError('Observation not in the discovered archive records: ' + ', '.join(missing))
    if observation_a == observation_b:
        raise ValueError('Select two different observations')
    epochs = []
    for ident in (observation_a, observation_b):
        emit('download', observation_id=ident)
        t0 = time.perf_counter()
        path = download_cutout(index[ident], center[0], center[1], SCIENCE_CACHE, size=CUTOUT_SIZE_DEG)
        timings[f'download_{len(epochs)}'] = time.perf_counter()-t0
        emit('calibration', observation_id=ident)
        t0 = time.perf_counter()
        epochs.append(load_epoch(path, ident, CALIBRATION_CACHE, product_identity=index[ident].get('access_url')))
        timings[f'calibration_{len(epochs)-1}'] = time.perf_counter()-t0
    emit('comparison')
    t0 = time.perf_counter()
    result = compare(epochs[0], epochs[1], catalog_sources=catalog, target=target, psf_match=psf_match)
    timings['comparison'] = time.perf_counter()-t0
    report = result['report']
    report['catalog'] = catalog_provenance
    report['science_files'] = science_files(epochs, index)
    report['software_sha256'] = software_sha256()
    report['cutout'] = {'center_deg': list(center), 'size_deg': CUTOUT_SIZE_DEG}
    report['timings_seconds'] = timings
    emit('products')
    files = write_products(result, output)
    report['elapsed_seconds'] = time.perf_counter()-started
    report['files'] = files
    from spherex_lab.lab.difference_imaging import _json_safe
    (Path(output)/'report.json').write_text(json.dumps(_json_safe(report), indent=2, allow_nan=False), encoding='utf-8')
    return _json_safe(report)
