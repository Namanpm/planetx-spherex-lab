# NASA Space Apps submission: SPHEREx Epoch Lab (Planet X × SPHEREx)

## One-paragraph description

SPHEREx observes every part of the sky many times, but each visit usually sees a star through a different colour, because its linear variable filters change wavelength across the detector. The SPHEREx Epoch Lab finds the visits that sampled the same sky at the same wavelength. It uses the archive's own footprints and calibration maps to screen all records in seconds, then downloads the real images and compares only the pixels where the wavelengths match. Every comparison is validated against its own data (registration, photometric scale, noise null test, injection–recovery). Every residual is explained before anything is called a candidate. On real SPHEREx DR1 data it produces a fully validated Barnard's Star field comparison, measures the proper motion of the brown dwarf ε Indi B from five SPHEREx epochs (11σ, consistent with its catalog value), and shows why naive whole-image differencing would mislead.

## Why it matters

Time-domain science needs like-for-like measurements. SPHEREx's all-sky, 102-colour survey revisits each field, and same-wavelength revisits turn it into a variability and motion survey. The lab finds those revisits and makes them usable without overclaiming.

## Demo script (about 3 minutes)

1. **Problem (20 s).** Open *Method*: SPHEREx's filters mean different visits usually see different colours.
2. **Screening (40 s).** *Explore* → preset *Barnard's Star* → *Find same-wavelength epochs*. Show the coverage chart: wavelength at the target over time, with detectors as bands. Orange points are same-wavelength pairs, and the table ranks them by matched area.
3. **Comparison (60 s).** *Recorded runs* → *Barnard's Star field*. Step through side by side, blink, difference and significance; hover for the pixel readout. Show the validation cards with all gates passing, and the residual table where every residual is explained and none needs follow-up.
4. **Known mover (40 s).** *ε Indi B five epochs*: the trajectory against the catalog prediction, and the proper motion measured from SPHEREx at 11σ. Then *one year apart*: the brown dwarf's dipole is labelled known proper motion, not variability.
5. **Honesty (20 s).** *RW Aurigae*: rotation leaves only a strip at the same wavelength, so the lab compares only that strip. The target is saturated and the lab says so.

## Screenshots

`docs/screenshots/`: home, recorded runs, time series, trajectory, comparison (side by side and significance), validation, and live screening.
