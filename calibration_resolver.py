"""Resolve archive-associated SPHEREx spectral calibration, fail closed on ambiguity."""
import hashlib
import io
import json
import os
import re
import tempfile
import time
import threading
from collections import OrderedDict, Counter
from pathlib import Path
from urllib.parse import urlparse

import numpy as np
import requests
from archive_http import retry_session
from astropy.io import fits
from astropy.io.votable import parse
from spectral_calibration import sample_calibration

DATALINK = 'https://irsa.ipac.caltech.edu/datalink/links/spherex'
MAX_BYTES = 80 * 1024 * 1024

def session():
    """Use a finite, slightly larger budget for intermittent DataLink 5xxs."""
    return retry_session(total=10, backoff=0.5)

def discover(publisher_id, detector, client=None):
    if not re.fullmatch(r'ivo://irsa\.ipac/spherex_dr1\?[A-Za-z0-9_]+/D[1-6]', publisher_id) or not publisher_id.endswith('/D' + str(detector)):
        raise ValueError('invalid publisher ID or detector mismatch')
    client = client or session()
    # Publisher IDs are query values. Requests performs the required escaping;
    # retrying equivalent encodings after a gateway error only multiplies load.
    # One call uses the finite policy from session() (at most eleven network attempts).
    r = client.get(DATALINK, params={'ID': publisher_id}, timeout=(10, 60))
    r.raise_for_status()
    v = parse(io.BytesIO(r.content))
    infos = [i for resource in v.resources for i in resource.infos if i.name == 'QUERY_STATUS']
    if any(i.value != 'OK' for i in infos): raise ValueError('DataLink returned an error or incomplete result')
    t = v.get_first_table().to_table(use_names_over_ids=True)
    if not {'id', 'semantics', 'description', 'access_url', 'error_message'} <= set(t.colnames): raise ValueError('malformed DataLink response')
    candidates = []
    for row in t:
        if str(row['id']) != publisher_id or str(row['semantics']) != '#calibration': continue
        if not str(row['description']).startswith('spectral_wcs_'): continue
        if not np.ma.is_masked(row['error_message']) and str(row['error_message']).strip(): raise ValueError('DataLink calibration error')
        url = str(row['access_url'])
        parsed = urlparse(url)
        if parsed.scheme != 'https' or parsed.netloc != 'irsa.ipac.caltech.edu' or not parsed.path.startswith('/ibe/data/spherex/dr1/spectral_wcs/') or parsed.query or parsed.fragment:
            raise ValueError('unsafe calibration URL')
        candidates.append(url)
    if len(candidates) != 1: raise ValueError(f'expected one spectral calibration, found {len(candidates)}')
    return candidates[0]

def validate_calibration(path, detector):
    with fits.open(path, checksum=True, memmap=False) as f:
        for name in ('CWAVE', 'CBAND'):
            if name not in f: raise ValueError('missing ' + name)
            h = f[name]
            if h.data is None or h.data.shape != (2040, 2040): raise ValueError('invalid calibration geometry')
            if str(h.header.get('DETECTOR', '')).strip() != str(detector): raise ValueError('calibration detector mismatch')
            if str(h.header.get('BUNIT', '')).strip().lower() not in ('um', 'micron') : raise ValueError('invalid calibration unit')
        return np.array(f['CWAVE'].data, copy=True), np.array(f['CBAND'].data, copy=True)

def fetch(url, detector, cache_dir, client=None, *, _validated=None):
    cache = Path(cache_dir); cache.mkdir(parents=True, exist_ok=True)
    path = cache / (hashlib.sha256(url.encode()).hexdigest() + '.fits')
    if path.exists():
        validate_calibration(path, detector)
        return path
    client = client or session()
    temp = None
    try:
        with client.get(url, stream=True, timeout=(10, 90), allow_redirects=False) as r:
            r.raise_for_status()
            if r.status_code != 200: raise ValueError('unexpected calibration HTTP status')
            if int(r.headers.get('Content-Length', '0')) > MAX_BYTES: raise ValueError('calibration exceeds size limit')
            with tempfile.NamedTemporaryFile(dir=cache, suffix='.part', delete=False) as f:
                temp = Path(f.name); size = 0
                for chunk in r.iter_content(1024 * 1024):
                    size += len(chunk)
                    if size > MAX_BYTES: raise ValueError('calibration exceeds size limit')
                    f.write(chunk)
        arrays = validate_calibration(temp, detector)
        if _validated is not None:
            _validated.extend(arrays)
        os.replace(temp, path); temp = None
        return path
    finally:
        if temp is not None: temp.unlink(missing_ok=True)


CACHE_SCHEMA = 1
VALIDATOR_VERSION = 1


def _trusted_url(url):
    parsed = urlparse(url)
    return (parsed.scheme == 'https' and parsed.netloc == 'irsa.ipac.caltech.edu'
            and parsed.path.startswith('/ibe/data/spherex/dr1/spectral_wcs/')
            and not parsed.query and not parsed.fragment)


class CalibrationCache:
    """Bounded process array cache and expiring, version-aware disk associations.

    The lock provides single-flight loading within a process. Atomic replacement
    permits independent processes to share disk records without partial reads.
    Array identities include versioned archive URL and local inode/size/timestamps;
    replacing or touching a file invalidates validation and its content digest.
    """
    def __init__(self, max_bytes=128 * 1024 * 1024, association_ttl=86400):
        if max_bytes < 0 or association_ttl < 0:
            raise ValueError('cache limits must be nonnegative')
        self.max_bytes = max_bytes
        self.association_ttl = association_ttl
        self._lock = threading.RLock()
        self._arrays = OrderedDict()
        self._bytes = 0
        self._stats = Counter()

    def clear(self):
        """Release process memory/counters; persistent associations are retained."""
        with self._lock:
            self._arrays.clear()
            self._bytes = 0
            self._stats.clear()

    def cache_info(self):
        with self._lock:
            counters = {name: self._stats[name] for name in (
                'association_hits', 'association_misses', 'datalink_requests', 'downloads',
                'array_hits', 'array_misses', 'validations', 'digests', 'evictions')}
            return dict(counters, array_entries=len(self._arrays), array_bytes=self._bytes,
                        max_bytes=self.max_bytes)

    def association(self, publisher_id, detector, identity, cache_dir, client=None):
        # A hit must never bypass the caller's detector/identity invariant.
        if not publisher_id.endswith('/D' + str(detector)):
            raise ValueError('invalid publisher ID or detector mismatch')
        key = {'schema': CACHE_SCHEMA, 'publisher_id': publisher_id,
               'science_product': identity, 'datalink_url': DATALINK}
        encoded = json.dumps(key, sort_keys=True, separators=(',', ':'))
        directory = Path(cache_dir) / 'associations'
        path = directory / (hashlib.sha256(encoded.encode()).hexdigest() + '.json')
        with self._lock:
            try:
                record = json.loads(path.read_text())
                age = time.time() - record['retrieved_at']
                if (record['key'] == key and 0 <= age < self.association_ttl
                        and isinstance(record['url'], str) and _trusted_url(record['url'])):
                    self._stats['association_hits'] += 1
                    return record['url']
            except (OSError, ValueError, KeyError, TypeError):
                pass
            self._stats['association_misses'] += 1
            self._stats['datalink_requests'] += 1
            url = discover(publisher_id, detector, client)
            directory.mkdir(parents=True, exist_ok=True)
            temp = None
            try:
                with tempfile.NamedTemporaryFile(mode='w', dir=directory, suffix='.part',
                                                 delete=False) as output:
                    temp = Path(output.name)
                    json.dump({'key': key, 'url': url, 'retrieved_at': time.time()}, output)
                    output.flush()
                    os.fsync(output.fileno())
                os.replace(temp, path)
            finally:
                if temp is not None:
                    temp.unlink(missing_ok=True)
            return url

    def arrays(self, url, detector, cache_dir, client=None):
        if not _trusted_url(url):
            raise ValueError('unsafe calibration URL')
        with self._lock:
            path = Path(cache_dir) / (hashlib.sha256(url.encode()).hexdigest() + '.fits')
            downloaded = not path.exists()
            loaded = []
            if downloaded:
                self._stats['downloads'] += 1
                # fetch validates streamed bytes before publication; loaded arrays
                # below remain the authoritative in-memory validated instance.
                path = fetch(url, detector, cache_dir, client, _validated=loaded)
            stat = path.stat()
            key = (url, str(path.resolve()), stat.st_dev, stat.st_ino,
                   stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns, VALIDATOR_VERSION)
            if key in self._arrays:
                wave, band, digest, stored_detector = self._arrays[key]
                if stored_detector != detector:
                    raise ValueError('calibration detector mismatch')
                self._arrays.move_to_end(key)
                self._stats['array_hits'] += 1
                return wave, band, digest
            self._stats['array_misses'] += 1
            wave, band = loaded if loaded else validate_calibration(path, detector)
            self._stats['validations'] += 1
            with path.open('rb') as stream:
                digest = hashlib.file_digest(stream, 'sha256').hexdigest()
            self._stats['digests'] += 1
            wave.flags.writeable = False
            band.flags.writeable = False
            size = wave.nbytes + band.nbytes
            # Discard stale identities of this path even when under budget.
            for old in list(self._arrays):
                if old[1] == key[1]:
                    value = self._arrays.pop(old)
                    self._bytes -= value[0].nbytes + value[1].nbytes
            while self._arrays and self._bytes + size > self.max_bytes:
                _, value = self._arrays.popitem(last=False)
                self._bytes -= value[0].nbytes + value[1].nbytes
                self._stats['evictions'] += 1
            if size <= self.max_bytes:
                self._arrays[key] = (wave, band, digest, detector)
                self._bytes += size
            return wave, band, digest


DEFAULT_CACHE = CalibrationCache()

def wavelength_maps(science_path, publisher_id, cache_dir, client=None, *, cache=DEFAULT_CACHE, product_identity=None):
    with fits.open(science_path, memmap=False) as f:
        h = f['IMAGE'].header.copy(); shape = f['IMAGE'].data.shape
        pipeline_version = str(f[0].header.get('VERSION', ''))
        detector = int(h['DETECTOR'])
        if publisher_id != f"ivo://irsa.ipac/spherex_dr1?{h['OBSID']}/D{detector}":
            raise ValueError('publisher ID does not match science observation')
        history = [str(x) for x in h.get('HISTORY', [])]
        uuids = sorted(set(m.group(1) for line in history if (m := re.search(r'\[CALIB\]\s+spectral_wcs\s+([0-9a-fA-F-]+)', line))))
        if len(uuids) != 1: raise ValueError('missing or ambiguous science spectral calibration UUID')
    # Pipeline version and complete processing history
    # invalidate associations after science reprocessing, independently of cutout.
    identity = {'pipeline_version': pipeline_version, 'archive_product': product_identity,
                'spectral_uuid': uuids[0],
                'history_sha256': hashlib.sha256('\n'.join(history).encode()).hexdigest()}
    url = cache.association(publisher_id, detector, identity, cache_dir, client)
    wave, band, digest = cache.arrays(url, detector, cache_dir, client)
    wave, band, valid = sample_calibration(h, shape, wave, band)
    return wave, band, valid, {'publisher_id': publisher_id, 'datalink_url': DATALINK, 'calibration_url': url, 'calibration_sha256': digest, 'detector': detector, 'science_calibration_uuid': uuids[0], 'uuid_match': 'unverified', 'association': 'IRSA DataLink #calibration', 'unit': 'um'}
