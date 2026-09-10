# Raman peak analysis - D, G and 2D bands

## Peak positions and intensities

| Sample | Band | Position (cm-1) | Shift vs nominal | Intensity (counts) | FWHM (cm-1) | Area |
|---|---|---|---|---|---|---|
| Control (uncoated) | D | 1352.2 | +2.2 | 559.7 | 97.8 | 66,382 |
| Control (uncoated) | G | 1581.5 | +1.5 | 806.9 | 65.0 | 67,121 |
| Control (uncoated) | 2D | 2685.3 | -14.7 | 333.1 | 85.9 | 33,802 |
| Treated (coating A) | D | 1372.7 | +22.7 | 476.7 | n/a | 75,669 |
| Treated (coating A) | G | 1585.7 | +5.7 | 727.9 | n/a | 73,697 |
| Treated (coating A) | 2D | not detected | - | - | - | - |

## Intensity ratios

| Sample | ID/IG | I2D/IG | ID/IG (area) | I2D/IG (area) |
|---|---|---|---|---|
| Control (uncoated) | 0.694 ± 0.017 | 0.413 ± 0.008 | 0.989 | 0.504 |
| Treated (coating A) | 0.655 ± 0.034 | n/a | 1.027 | n/a |

## Treated vs control

| Metric | Treated | Control | Δ (treated − control) | % change |
|---|---|---|---|---|
| ID/IG | 0.655 | 0.694 | -0.039 | not comparable * |
| D position (cm-1) | 1372.7 | 1352.2 | +20.5 | - |
| G position (cm-1) | 1585.7 | 1581.5 | +4.2 | - |

\* ID/IG was measured by Lorentzian fit on one sample and by peak maximum on the other, because the fit was rejected for the other. Those are different measurements, so the difference between them is not a result. Compare the samples using a single estimator instead.

## Spectra analysed

- **Control (uncoated)** - `data/control_01.txt`, 1938 points, 1-3000 cm-1
- **Treated (coating A)** - `data/treated_01.txt`, 7060 points, 3-3000 cm-1
  - note: D band: Lorentzian fit rejected as implausible, height taken from the smoothed maximum instead
  - note: G band: Lorentzian fit rejected as implausible, height taken from the smoothed maximum instead
  - note: 2D band not detected - peak stands only 3.7 sigma above the noise, below the 5 sigma threshold
