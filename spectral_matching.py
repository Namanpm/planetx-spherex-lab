"""Find epochs in which SPHEREx observed the same sky position at the same wavelength.

SPHEREx disperses light with linear variable filters, so the wavelength sampled
by a pixel depends on where the sky position lands on the detector. Two
observations of one position are spectrally comparable only where their
per-pixel CWAVE values agree. Archive search results carry a detector-wide
wavelength range, which is not sufficient.

Screening uses two archive facts per record:
  * the SIA footprint polygon (``s_region``), whose four vertices correspond to
    detector corners (0,0), (2039,0), (2039,2039), (0,2039) in that order; a
    gnomonic homography through them estimates the detector pixel of a sky point;
  * the detector's DataLink-associated CWAVE/CBAND calibration map.

The homography ignores optical distortion. Its error was measured against real
downloaded cutouts (``validate_footprint_model``, fixture
tests_data/footprint_ground_truth.json): median 3.1 px, maximum 5.7 px over 26
distinct observations. Screening therefore nominates candidates; every comparison
is confirmed on the downloaded cutouts with their own calibration association.
"""
import numpy as np

from calibration_resolver import DEFAULT_CACHE, discover

DETECTOR_SIZE = 2040
# Vertex order validated against downloaded cutouts (see module docstring).
FOOTPRINT_CORNERS = np.array([[0., 0.], [2039., 0.], [2039., 2039.], [0., 2039.]])
ASSOCIATION_TTL = 86400
POLICY = {
    'tolerance_fraction': 0.10,       # |delta lambda| <= 0.10 x CBAND at the target pixel
    'min_separation_days': 1.0,       # distinct epochs only
    'footprint_error_px': 6.0,        # conservative bound from measured 5.7 px maximum
    'max_pairs': 200,
}


def _tangent(ra, dec, ra0, dec0):
    ra, dec, ra0, dec0 = (np.deg2rad(np.asarray(v, dtype=float)) for v in (ra, dec, ra0, dec0))
    cosc = np.sin(dec0)*np.sin(dec) + np.cos(dec0)*np.cos(dec)*np.cos(ra-ra0)
    if np.any(cosc <= 0):
        raise ValueError('point is more than 90 degrees from the footprint centre')
    xi = np.cos(dec)*np.sin(ra-ra0)/cosc
    eta = (np.cos(dec0)*np.sin(dec) - np.sin(dec0)*np.cos(dec)*np.cos(ra-ra0))/cosc
    return xi, eta


def _homography(src, dst):
    rows = []
    for (x, y), (u, v) in zip(src, dst):
        rows.append([x, y, 1, 0, 0, 0, -u*x, -u*y, -u])
        rows.append([0, 0, 0, x, y, 1, -v*x, -v*y, -v])
    _, singular, vt = np.linalg.svd(np.asarray(rows))
    if singular[-2] <= singular[0]*1e-12:
        raise ValueError('degenerate footprint polygon')
    return vt[-1].reshape(3, 3)


def footprint_transform(s_region):
    """Tangent-plane centre and homography for one four-corner detector footprint."""
    parts = str(s_region).split()
    if parts[:2] != ['POLYGON', 'ICRS']:
        raise ValueError('unsupported footprint')
    vertices = np.asarray(parts[2:], dtype=float).reshape(-1, 2)
    if len(vertices) == 5 and np.allclose(vertices[0], vertices[-1]):
        vertices = vertices[:-1]
    if len(vertices) != 4 or not np.isfinite(vertices).all():
        raise ValueError('footprint is not a four-corner detector polygon')
    ra0, dec0 = vertices.mean(axis=0)
    if np.ptp(vertices[:, 0]) > 180:  # RA wrap: average on the unit circle
        ra0 = float(np.rad2deg(np.angle(np.exp(1j*np.deg2rad(vertices[:, 0])).mean())) % 360)
    src = np.column_stack(_tangent(vertices[:, 0], vertices[:, 1], ra0, dec0))
    return float(ra0), float(dec0), _homography(src, FOOTPRINT_CORNERS)


def apply_transform(transform, ra, dec):
    """Vectorised detector (x, y); NaN outside the detector."""
    ra0, dec0, h = transform
    xi, eta = _tangent(ra, dec, ra0, dec0)
    p = h @ np.vstack([np.ravel(xi), np.ravel(eta), np.ones(np.size(xi))])
    x, y = p[0]/p[2], p[1]/p[2]
    inside = np.isfinite(x) & np.isfinite(y) & (x >= -.5) & (x <= DETECTOR_SIZE-.5) & (y >= -.5) & (y <= DETECTOR_SIZE-.5)
    return np.where(inside, x, np.nan).reshape(np.shape(xi)), np.where(inside, y, np.nan).reshape(np.shape(xi))


def footprint_pixel(s_region, ra, dec):
    """Estimated zero-based detector (x, y) of an ICRS position inside a footprint."""
    x, y = apply_transform(footprint_transform(s_region), np.asarray([ra], float), np.asarray([dec], float))
    if not np.isfinite(x[0]) or not np.isfinite(y[0]):
        raise ValueError('position falls outside the detector footprint')
    return float(x[0]), float(y[0])


def field_grid(ra, dec, half_width_deg, n=21):
    """Sky positions on a square tangent-plane grid centred on (ra, dec)."""
    offsets = np.deg2rad(np.linspace(-half_width_deg, half_width_deg, n))
    xi, eta = np.meshgrid(offsets, offsets)
    ra0, dec0 = np.deg2rad(ra), np.deg2rad(dec)
    rho = np.hypot(xi, eta); c = np.arctan(rho)
    with np.errstate(invalid='ignore', divide='ignore'):
        dec_out = np.arcsin(np.cos(c)*np.sin(dec0) + np.where(rho > 0, eta*np.sin(c)*np.cos(dec0)/rho, 0))
        ra_out = ra0 + np.arctan2(xi*np.sin(c), rho*np.cos(dec0)*np.cos(c) - eta*np.sin(dec0)*np.sin(c))
    return np.rad2deg(ra_out) % 360, np.rad2deg(dec_out)


def _sample(calibration, transform, ra, dec):
    x, y = apply_transform(transform, ra, dec)
    good = np.isfinite(x) & np.isfinite(y)
    ix = np.clip(np.rint(np.nan_to_num(x)).astype(int), 0, DETECTOR_SIZE-1); iy = np.clip(np.rint(np.nan_to_num(y)).astype(int), 0, DETECTOR_SIZE-1)
    return np.where(good, calibration['wave'][iy, ix], np.nan), np.where(good, calibration['band'][iy, ix], np.nan)


def predicted_matched_fraction(entry_a, entry_b, maps, center, half_width_deg, tolerance_fraction, n=21):
    """Fraction of a cutout-sized grid where both epochs' estimated CWAVE agree."""
    ra, dec = field_grid(center[0], center[1], half_width_deg, n)
    wa, ba = _sample(maps[entry_a['detector']], footprint_transform(entry_a['s_region']), ra, dec)
    wb, _ = _sample(maps[entry_b['detector']], footprint_transform(entry_b['s_region']), ra, dec)
    both = np.isfinite(wa) & np.isfinite(wb) & (ba > 0)
    matched = both & (np.abs(wb-wa) <= tolerance_fraction*ba)
    return float(matched.sum()/ra.size), float(both.sum()/ra.size)


def detector_maps(records, cache_dir, *, client=None, cache=DEFAULT_CACHE, attempts=3):
    """Resolve one DataLink calibration map per detector present in ``records``.

    A detector whose association cannot be resolved is reported with its error;
    its observations are excluded from screening, never assigned another map.
    """
    by_detector = {}
    for record in records:
        pid = record.get('obs_publisher_did') or ''
        if pid[-3:-1] == '/D' and pid[-1] in '123456':
            by_detector.setdefault(int(pid[-1]), []).append(pid)
    maps = {}
    # Per-detector associations are remembered for ASSOCIATION_TTL seconds (as in
    # calibration_resolver); arrays are validated on every load by cache.arrays.
    import json, os, tempfile, time
    from calibration_resolver import _trusted_url
    record_path = os.path.join(str(cache_dir), 'screening_associations.json')
    try:
        with open(record_path, encoding='utf-8') as f:
            remembered = json.load(f)
    except (OSError, ValueError):
        remembered = {}
    changed = False
    for detector, ids in sorted(by_detector.items()):
        errors = []
        entry = remembered.get(str(detector))
        if entry and 0 <= time.time() - entry.get('retrieved_at', -1) < ASSOCIATION_TTL and _trusted_url(entry.get('url', '')):
            try:
                wave, band, digest = cache.arrays(entry['url'], detector, cache_dir, client)
                maps[detector] = {'wave': wave, 'band': band, 'url': entry['url'], 'sha256': digest,
                                  'association_observation': entry['observation'], 'association': 'remembered DataLink result'}
                continue
            except Exception as exc:
                errors.append(f"remembered association: {type(exc).__name__}: {exc}")
        for pid in ids[:attempts]:
            try:
                url = discover(pid, detector, client)
                wave, band, digest = cache.arrays(url, detector, cache_dir, client)
                maps[detector] = {'wave': wave, 'band': band, 'url': url, 'sha256': digest,
                                  'association_observation': pid, 'association': 'live DataLink'}
                remembered[str(detector)] = {'url': url, 'observation': pid, 'retrieved_at': time.time()}
                changed = True
                break
            except Exception as exc:  # recorded per detector; never silently substituted
                errors.append(f'{pid}: {type(exc).__name__}: {exc}')
        else:
            maps[detector] = {'error': errors}
    if changed:
        os.makedirs(str(cache_dir), exist_ok=True)
        with tempfile.NamedTemporaryFile('w', dir=str(cache_dir), suffix='.part', delete=False, encoding='utf-8') as f:
            json.dump(remembered, f); temporary = f.name
        os.replace(temporary, record_path)
    return maps


def _local_gradient(wave, x, y):
    ix, iy = int(np.clip(round(x), 1, wave.shape[1]-2)), int(np.clip(round(y), 1, wave.shape[0]-2))
    gx = (float(wave[iy, ix+1]) - float(wave[iy, ix-1]))/2
    gy = (float(wave[iy+1, ix]) - float(wave[iy-1, ix]))/2
    return float(np.hypot(gx, gy))


def screen_observations(records, ra, dec, maps, policy=POLICY):
    """Per-record estimated wavelength at the (optionally proper-motion propagated) position."""
    screened = []
    for record in records:
        pid = record.get('obs_publisher_did') or ''
        entry = {'observation_id': pid, 'obs_id': record.get('obs_id'), 'mjd': record.get('t_min'),
                 'detector': int(pid[-1]) if pid[-1:] in tuple('123456') else None,
                 'access_url': record.get('access_url'), 's_region': record.get('s_region')}
        position = record.get('target_position') or {'ra_deg': ra, 'dec_deg': dec}
        entry['position_deg'] = [float(position['ra_deg']), float(position['dec_deg'])]
        try:
            if entry['detector'] is None:
                raise ValueError('publisher ID lacks detector')
            if entry['mjd'] is None or not np.isfinite(entry['mjd']):
                raise ValueError('missing observation time')
            calibration = maps.get(entry['detector'], {})
            if 'wave' not in calibration:
                raise ValueError('detector calibration unavailable')
            x, y = footprint_pixel(record['s_region'], *entry['position_deg'])
            ix, iy = int(round(x)), int(round(y))
            wave, band = float(calibration['wave'][iy, ix]), float(calibration['band'][iy, ix])
            if not (np.isfinite(wave) and np.isfinite(band) and wave > 0 and band > 0):
                raise ValueError('calibration map invalid at estimated pixel')
            gradient = _local_gradient(calibration['wave'], x, y)
            entry.update(status='SCREENED', detector_xy=[x, y], wavelength_um=wave, bandwidth_um=band,
                         wavelength_uncertainty_um=gradient*policy['footprint_error_px'],
                         calibration_url=calibration['url'])
        except (ValueError, KeyError, IndexError) as exc:
            entry.update(status='UNSCREENED', reason=str(exc))
        screened.append(entry)
    screened.sort(key=lambda e: (e['mjd'] if e['mjd'] is not None else np.inf, e['observation_id']))
    return screened


def matched_pairs(screened, policy=POLICY, *, maps=None, center=None, half_width_deg=0.15):
    """Same-detector, distinct-epoch pairs whose screened wavelengths agree at the target.

    With ``maps`` and ``center``, each pair also gets the predicted fraction of a
    cutout-sized field that is spectrally matched; pairs are ranked by that area.
    """
    good = [e for e in screened if e['status'] == 'SCREENED']
    pairs = []
    for i, a in enumerate(good):
        for b in good[i+1:]:
            if a['detector'] != b['detector'] or abs(b['mjd']-a['mjd']) < policy['min_separation_days']:
                continue
            band = min(a['bandwidth_um'], b['bandwidth_um'])
            delta = b['wavelength_um'] - a['wavelength_um']
            if abs(delta) > policy['tolerance_fraction']*band:
                continue
            first, second = (a, b) if a['mjd'] <= b['mjd'] else (b, a)
            pairs.append({'observation_a': first['observation_id'], 'observation_b': second['observation_id'],
                          'detector': a['detector'], 'separation_days': float(second['mjd']-first['mjd']),
                          'wavelength_um': float((a['wavelength_um']+b['wavelength_um'])/2),
                          'delta_wavelength_um': float(delta), 'delta_over_bandwidth': float(abs(delta)/band),
                          'screening_uncertainty_um': float(np.hypot(a['wavelength_uncertainty_um'], b['wavelength_uncertainty_um'])),
                          'status': 'SCREENED_CANDIDATE'})
            if maps is not None and center is not None:
                area, common = predicted_matched_fraction(first, second, maps, center, half_width_deg, policy['tolerance_fraction'])
                pairs[-1].update(predicted_matched_fraction=area, predicted_overlap_fraction=common)
    pairs.sort(key=lambda p: (-p.get('predicted_matched_fraction', 0.), p['delta_over_bandwidth'], -p['separation_days']))
    return pairs[:policy['max_pairs']]


def matched_groups(screened, policy=POLICY):
    """Chains of distinct epochs on one detector at one wavelength (light-curve candidates)."""
    good = sorted((e for e in screened if e['status'] == 'SCREENED'), key=lambda e: (e['detector'], e['wavelength_um']))
    groups = []
    for e in good:
        g = groups[-1] if groups else None
        if g and g['detector'] == e['detector'] and abs(e['wavelength_um']-g['wavelength_um']) <= policy['tolerance_fraction']*min(g['bandwidth_um'], e['bandwidth_um'])/2:
            g['members'].append(e)
        else:
            groups.append({'detector': e['detector'], 'wavelength_um': e['wavelength_um'], 'bandwidth_um': e['bandwidth_um'], 'members': [e]})
    result = []
    for g in groups:
        epochs = []
        for m in sorted(g['members'], key=lambda m: m['mjd']):
            if not epochs or m['mjd'] - epochs[-1]['mjd'] >= policy['min_separation_days']:
                epochs.append(m)
        if len(epochs) >= 2:
            result.append({'detector': g['detector'], 'wavelength_um': float(np.mean([m['wavelength_um'] for m in epochs])),
                           'epoch_count': len(epochs), 'observation_ids': [m['observation_id'] for m in epochs],
                           'mjd': [m['mjd'] for m in epochs]})
    result.sort(key=lambda g: (-g['epoch_count'], g['detector'], g['wavelength_um']))
    return result


def screen(*, target_id=None, ra=None, dec=None, cache_dir, client=None, policy=POLICY, progress=None):
    """Discovery plus spectral screening for one target or fixed position."""
    from pair_discovery import discover_pairs
    def emit(stage, **fields):
        if progress:
            progress({'stage': stage, **fields})
    emit('archive_discovery')
    discovery = discover_pairs(ra=ra, dec=dec, target_id=target_id, limit=1000)
    records = discovery['observations']
    ra, dec = discovery['query']['ra_deg'], discovery['query']['dec_deg']
    emit('calibration_maps', observations=len(records))
    maps = detector_maps(records, cache_dir, client=client)
    emit('screening', detectors=sorted(maps))
    screened = screen_observations(records, ra, dec, maps, policy)
    pairs = matched_pairs(screened, policy, maps=maps, center=(ra, dec))
    return {'query': discovery['query'], 'target_id': target_id, 'policy': policy,
            'observation_count': len(records), 'rejected_records': discovery['rejected'],
            'screened_count': sum(e['status'] == 'SCREENED' for e in screened),
            'detector_calibrations': {str(d): {k: v for k, v in m.items() if k not in ('wave', 'band')} for d, m in maps.items()},
            'observations': screened, 'matched_pairs': pairs, 'matched_groups': matched_groups(screened, policy),
            'method': 'Footprint-homography detector pixel + DataLink CWAVE/CBAND map; screening estimate, confirmed per pair on downloaded cutouts.',
            'limitations': ['Screening ignores optical distortion (measured error: median 3.1 px, max 5.7 px over 26 real cutouts).',
                            'One calibration association per detector is used for screening; each comparison re-resolves its own.',
                            'Archive coverage is bounded by the SIA query radius; absence of a match is not absence of data.']}


def validate_footprint_model(cases):
    """Compare footprint estimates with true detector pixels from downloaded cutouts.

    ``cases`` items: ``{'s_region', 'ra_deg', 'dec_deg', 'detector_xy'}`` where
    ``detector_xy`` was read from a real cutout's alternate WCS.
    """
    errors = []
    for case in cases:
        x, y = footprint_pixel(case['s_region'], case['ra_deg'], case['dec_deg'])
        errors.append(float(np.hypot(x-case['detector_xy'][0], y-case['detector_xy'][1])))
    errors = np.asarray(errors)
    return {'count': int(errors.size), 'median_error_px': float(np.median(errors)),
            'max_error_px': float(errors.max()), 'errors_px': errors.tolist()}
