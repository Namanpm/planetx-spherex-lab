"""Explicit, bounded selection of real archive observation records."""
import re
from itertools import combinations
import numpy as np

MAX_OBSERVATIONS = 8
MODES = ('best_pair', 'all_valid_pairs', 'time_series')
PID = re.compile(r'^ivo://irsa\.ipac/spherex_dr1\?[A-Za-z0-9_]+/D[1-6]$')


def validate_coordinates(ra, dec):
    if ra is None and dec is None:
        return
    if ra is None or dec is None or not np.isfinite([ra, dec]).all() or not 0 <= ra < 360 or not -90 <= dec <= 90:
        raise ValueError('Supply finite RA [0,360) and declination [-90,90] together')


def select_observations(records, selected_ids=None, *, maximum=8, mode='best_pair', upper_bound=MAX_OBSERVATIONS):
    if mode not in MODES:
        raise ValueError('Unknown analysis mode')
    if not isinstance(maximum, int) or isinstance(maximum, bool) or not 2 <= maximum <= upper_bound:
        raise ValueError(f'Observation count must be between 2 and {upper_bound}')
    if len(records) > 10000:
        raise ValueError('Archive record budget exceeded')
    indexed = {}
    for record in records:
        ident = record.get('obs_publisher_did', '')
        if not PID.fullmatch(ident):
            raise ValueError('Invalid observation publisher identity')
        if ident in indexed:
            if record != indexed[ident]:
                raise ValueError('Conflicting archive records for the same observation')
            continue
        if not np.isfinite(record.get('t_min', np.nan)):
            raise ValueError('Observation timestamp is missing or invalid')
        detector = int(ident[-1])
        if record.get('energy_bandpassname') != f'SPHEREx-D{detector}':
            raise ValueError('Observation detector metadata mismatch')
        indexed[ident] = record
    if selected_ids is not None:
        if not 2 <= len(selected_ids) <= maximum or len(set(selected_ids)) != len(selected_ids):
            raise ValueError('Select distinct observation IDs within the configured limit')
        missing = sorted(set(selected_ids) - indexed.keys())
        if missing:
            raise ValueError('Selected observations absent from supplied archive records: ' + ', '.join(missing))
        selected = [indexed[i] for i in selected_ids]
    else:
        ordered = sorted(indexed.values(), key=lambda r:(r['t_min'], r['obs_publisher_did']))
        indices = np.unique(np.rint(np.linspace(0,len(ordered)-1,min(maximum,len(ordered)))).astype(int)) if ordered else []
        selected = [ordered[i] for i in indices]
    if len(selected) < 2:
        raise ValueError('At least two real observations are required')
    return sorted(selected, key=lambda r:(r['t_min'],r['obs_publisher_did']))


def candidate_pairs(records):
    return [{'observation_a':a['obs_publisher_did'], 'observation_b':b['obs_publisher_did']} for a,b in combinations(records,2)]
