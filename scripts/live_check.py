from concurrent.futures import ThreadPoolExecutor
from spherex_lab.data.targets import targets
from spherex_lab.data.spherex import search

def run(t):
    if 'ra_deg' not in t: return (t['id'], 'unverified coordinates')
    try:
        r = search(t['ra_deg'], t['dec_deg'], .05)
        times = [x['t_min'] for x in r['observations'] if x['t_min'] is not None]
        return (t['id'], r['count'], min(times) if times else None, max(times) if times else None)
    except Exception as exc: return (t['id'], type(exc).__name__, str(exc))
with ThreadPoolExecutor(max_workers=4) as pool:
    for x in pool.map(run, targets()): print(x)
