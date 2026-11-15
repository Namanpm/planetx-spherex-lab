import argparse
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import requests
from astropy.io import fits
from astropy.wcs import WCS
from spherex_lab.data.spherex import search


def main():
    p = argparse.ArgumentParser(description='Download and plot one SPHEREx spectral-image cutout')
    p.add_argument('--ra', type=float, required=True)
    p.add_argument('--dec', type=float, required=True)
    p.add_argument('--radius', type=float, default=0.01)
    p.add_argument('--size', type=float, default=0.1, help='Cutout side length in degrees')
    p.add_argument('--collection', default='spherex_dr1', choices=('spherex_dr1', 'spherex_qr3'))
    p.add_argument('--output', default='artifacts')
    args = p.parse_args()
    if not 0 < args.size <= 0.1: p.error('size must be in (0, 0.1] degrees')
    result = search(args.ra, args.dec, args.radius, args.collection)
    if not result['observations']: raise SystemExit('No observations returned')
    row = result['observations'][0]
    response = requests.get(row['access_url'], params={'center': f'{args.ra},{args.dec}', 'size': args.size}, timeout=(10, 90))
    response.raise_for_status()
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    path = out / 'sample.fits'
    path.write_bytes(response.content)
    with fits.open(path) as hdus:
        hdu = hdus['IMAGE']
        data = hdu.data
        valid = np.isfinite(data)
        if not valid.any(): raise ValueError('No finite science image pixels')
        vmin, vmax = np.percentile(data[valid], [2, 98])
        fig = plt.figure(figsize=(8, 7))
        ax = fig.add_subplot(projection=WCS(hdu.header))
        im = ax.imshow(data, origin='lower', cmap='inferno', vmin=vmin, vmax=vmax)
        ax.coords[0].set_axislabel('Right ascension (ICRS)')
        ax.coords[1].set_axislabel('Declination (ICRS)')
        ax.set_title(f"SPHEREx {row['energy_bandpassname']} | {hdu.header.get('DATE-OBS', 'date unknown')}\n{row['obs_id']}")
        fig.colorbar(im, ax=ax, label=f"Surface brightness ({hdu.header.get('BUNIT', 'unit unknown')}); 2–98% display stretch")
        fig.tight_layout()
        fig.savefig(out / 'sky.png', dpi=150)
        print(f"records={result['count']} obs_id={row['obs_id']} date={hdu.header.get('DATE-OBS')} band={row['energy_bandpassname']} wavelength_range_m={row['em_min']},{row['em_max']} shape={data.shape} unit={hdu.header.get('BUNIT')} fits={path} plot={out / 'sky.png'} url={response.url}")

if __name__ == '__main__': main()
