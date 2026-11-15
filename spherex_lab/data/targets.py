"""Catalog positions are SIMBAD ICRS positions at J2000; proper motion is mas/year."""
from astropy.coordinates import SkyCoord
import astropy.units as u

_RAW = {
 'wise-0855': ('WISE 0855−0714', '08 55 16.53798', '-07 14 49.3728', -8123.7, 673.2, 'https://simbad.u-strasbg.fr/simbad/sim-id?Ident=WISEA+J085510.74-071442.5'),
 'barnards-star': ("Barnard's Star", '17 57 48.4984700685', '+04 41 36.113879676', -801.551, 10362.394, "https://simbad.cds.unistra.fr/simbad/sim-basic?Ident=Barnard%27s+star"),
 'epsilon-indi-b': ('Epsilon Indi B', '22 04 10.5951812451', '-56 46 58.239102418', 3981.977, -2466.832, 'https://simbad.u-strasbg.fr/simbad/sim-id?Ident=Epsilon+Indi+B'),
 'rw-aurigae': ('RW Aurigae', '05 07 49.5662483', '+30 24 05.177426', 7.23, -20.33, 'https://simbad.cds.unistra.fr/simbad/sim-id?Ident=RW+Aur'),
 'r-draconis': ('R Draconis', '16 32 40.2265999080', '+66 45 17.925901992', -20.643, 6.149, 'https://simbad.u-strasbg.fr/simbad/sim-id?Ident=R+Dra'),
}

def targets():
    result = []
    for key, (name, ra, dec, pmra, pmdec, source) in _RAW.items():
        c = SkyCoord(ra, dec, unit=(u.hourangle, u.deg), frame='icrs')
        result.append(dict(id=key, name=name, ra_deg=c.ra.deg, dec_deg=c.dec.deg, epoch='J2000', pm_ra_cosdec_mas_yr=pmra, pm_dec_mas_yr=pmdec, source=source))
    parallax = {'wise-0855': (439.0, 2.4), 'barnards-star': (546.9759, 0.0401), 'epsilon-indi-b': (270.6580, 0.6896), 'rw-aurigae': (15.42, 3.69), 'r-draconis': (1.2509, 0.0648)}
    for target in result:
        target['frame'] = 'ICRS'
        target['parallax_mas'], target['parallax_uncertainty_mas'] = parallax[target['id']]
        target['position_uncertainty'] = None
        target['astrometry_note'] = 'SIMBAD position at J2000; parallax and proper motion as listed by SIMBAD; observer parallax not applied.'
    return result


def barycentric_position(target_id, obstime):
    """Approximate ICRS barycentric direction with proper motion; no parallax or radial velocity.

    obstime accepts any Astropy Time value. Output is not an apparent spacecraft position.
    """
    from astropy.time import Time
    target = next((t for t in targets() if t['id'] == target_id), None)
    if target is None:
        raise ValueError('Unknown target')
    coord = SkyCoord(ra=target['ra_deg'] * u.deg, dec=target['dec_deg'] * u.deg,
        pm_ra_cosdec=target['pm_ra_cosdec_mas_yr'] * u.mas / u.yr,
        pm_dec=target['pm_dec_mas_yr'] * u.mas / u.yr, obstime=Time('J2000'), frame='icrs')
    result = coord.apply_space_motion(new_obstime=Time(obstime))
    return {'ra_deg': result.ra.deg, 'dec_deg': result.dec.deg, 'frame': 'ICRS',
            'time': Time(obstime).isot, 'limitations': 'Barycentric proper motion only; parallax and radial velocity not applied.'}
