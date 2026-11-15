"""SPHEREx Epoch Lab processing code.

    data/           archive access: SIA search, DataLink calibration, Gaia DR3, JPL Horizons, preset targets
    imaging/        pixel-level work: cutout reading, ePSF, WCS, source detection, spectral and photometric calibration
    models/         fitted models: astrometry, motion, spectra, covariance, comparison statistics
    lab/            Epoch Lab: same-wavelength screening, difference imaging, time series, /api/lab router
    repeatability/  earlier paths: pair investigation and the Gaia reference-star repeatability experiment
"""
from pathlib import Path

PACKAGE = Path(__file__).resolve().parent
ROOT = PACKAGE.parent  # repository root: artifacts/, static/, docs/


def source_file(name):
    """Path of a processing module from its file name, e.g. 'psf_processing.py' (used for software hashes)."""
    matches = list(PACKAGE.glob(f'*/{name}'))
    if len(matches) != 1:
        raise FileNotFoundError(f'{name}: expected one module in {PACKAGE}, found {len(matches)}')
    return matches[0]
