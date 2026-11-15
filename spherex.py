import argparse
import hashlib
import io
import json
import os
import tempfile
import threading
import time
from collections import OrderedDict

import numpy as np
from astropy.io.votable import parse
from archive_http import retry_session

ENDPOINT = 'https://irsa.ipac.caltech.edu/SIA'
COLLECTIONS = ('spherex_dr1', 'spherex_qr3')
_cache = OrderedDict()
_cache_lock = threading.Lock()
# Optional on-disk cache of raw SIA VOTables (DR1 records change rarely). Disabled
# unless a directory is set; entries older than DISK_TTL are refetched, and a failed
# request never falls back to a stale entry.
DISK_CACHE_DIR = None
DISK_TTL = 86400
TTL = 900
_session = retry_session(total=3, backoff=0.5)


def _get(url, **kwargs):
    """Indirection kept patchable for transport-free callers and unit tests."""
    return _session.get(url, **kwargs)


def check_query_status(votable):
    infos = [info for resource in votable.resources for info in resource.infos if info.name == 'QUERY_STATUS']
    if not infos or any(info.value != 'OK' for info in infos):
        raise ValueError('IRSA query failed or returned an incomplete result: ' + ', '.join(str(i.value) for i in infos))


def search(ra, dec, radius, collection='spherex_dr1'):
    if not (np.isfinite(ra) and 0 <= ra < 360 and np.isfinite(dec) and -90 <= dec <= 90 and np.isfinite(radius) and 0 < radius <= 0.1):
        raise ValueError('ra must be [0, 360), dec [-90, 90], and radius (0, 0.1] degrees')
    if collection not in COLLECTIONS:
        raise ValueError('unsupported collection')
    key = (ra, dec, radius, collection)
    now = time.monotonic()
    with _cache_lock:
        entry = _cache.get(key)
        if entry is not None and now - entry[0] < TTL:
            _cache.move_to_end(key)
            return entry[1]
    # SIA discovery uses a bounded idempotent-GET policy: initial try + 3
    # retries for connection, timeout, rate-limit, and upstream 5xx failures.
    params = {'COLLECTION': collection, 'POS': f'CIRCLE {ra} {dec} {radius}'}
    content = None; disk = None
    if DISK_CACHE_DIR is not None:
        disk = os.path.join(DISK_CACHE_DIR, hashlib.sha256(json.dumps(params, sort_keys=True).encode()).hexdigest() + '.xml')
        try:
            if time.time() - os.path.getmtime(disk) < DISK_TTL:
                with open(disk, 'rb') as f:
                    content = f.read()
        except OSError:
            content = None
    if content is None:
        response = _get(ENDPOINT, params=params, timeout=(10, 60))
        response.raise_for_status()
        content = response.content
        if disk is not None:
            os.makedirs(DISK_CACHE_DIR, exist_ok=True)
            with tempfile.NamedTemporaryFile(dir=DISK_CACHE_DIR, suffix='.part', delete=False) as f:
                f.write(content); temporary = f.name
            os.replace(temporary, disk)
    votable = parse(io.BytesIO(content))
    check_query_status(votable)
    table = votable.get_first_table().to_table(use_names_over_ids=True)
    fields = ('obs_id', 's_ra', 's_dec', 't_min', 't_max', 'em_min', 'em_max', 'energy_bandpassname', 'access_url', 'access_format', 's_region', 'obs_publisher_did')
    missing = set(fields) - set(table.colnames)
    if missing:
        raise ValueError(f'IRSA response missing fields: {sorted(missing)}')
    def value(v):
        if np.ma.is_masked(v): return None
        if isinstance(v, np.generic): v = v.item()
        if isinstance(v, bytes): v = v.decode('utf-8')
        if isinstance(v, float) and not np.isfinite(v): return None
        return v
    records = [{field: value(row[field]) for field in fields} for row in table]
    result = {'query': {'ra_deg': ra, 'dec_deg': dec, 'radius_deg': radius, 'collection': collection}, 'count': len(records), 'units': {'t_min': 'MJD days', 't_max': 'MJD days', 'em_min': 'm', 'em_max': 'm', 's_ra': 'deg', 's_dec': 'deg'}, 'observations': records}
    with _cache_lock:
        _cache[key] = (now, result)
        _cache.move_to_end(key)
        while len(_cache) > 32: _cache.popitem(last=False)
    return result


def main():
    parser = argparse.ArgumentParser(description='Query IRSA SPHEREx spectral images')
    parser.add_argument('--ra', type=float, required=True)
    parser.add_argument('--dec', type=float, required=True)
    parser.add_argument('--radius', type=float, default=0.01, help='Search radius in degrees (max 0.1)')
    parser.add_argument('--collection', choices=COLLECTIONS, default=COLLECTIONS[0])
    args = parser.parse_args()
    print(json.dumps(search(args.ra, args.dec, args.radius, args.collection), indent=2, allow_nan=False))

if __name__ == '__main__': 
    main()
