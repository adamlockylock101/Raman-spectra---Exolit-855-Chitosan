# Raman peak analysis - D, G and 2D bands

## Peak positions and intensities

| Sample | Band | Position (cm-1) | Shift vs nominal | Intensity (counts) | FWHM (cm-1) | Area |
|---|---|---|---|---|---|---|
| Control (uncoated) | D | 1352.5 | +2.5 | 1,143.5 | 66.3 | 99,257 |
| Control (uncoated) | G | 1586.1 | +6.1 | 1,015.4 | 56.8 | 73,890 |
| Control (uncoated) | 2D | 2701.8 | +1.8 | 279.1 | 95.9 | 21,797 |
| Treated (coating A) | D | 1348.5 | -1.5 | 716.9 | 52.4 | 56,050 |
| Treated (coating A) | G | 1583.1 | +3.1 | 1,010.9 | 48.6 | 66,126 |
| Treated (coating A) | 2D | 2694.7 | -5.3 | 466.6 | 87.7 | 37,667 |

## Intensity ratios

| Sample | ID/IG | I2D/IG | ID/IG (area) | I2D/IG (area) |
|---|---|---|---|---|
| Control (uncoated) | 1.126 ± 0.006 | 0.275 ± 0.003 | 1.343 | 0.295 |
| Treated (coating A) | 0.709 ± 0.005 | 0.462 ± 0.003 | 0.848 | 0.570 |

## Treated vs control

| Metric | Treated | Control | Δ (treated − control) | % change |
|---|---|---|---|---|
| ID/IG | 0.709 | 1.126 | -0.417 | -37.0% |
| I2D/IG | 0.462 | 0.275 | +0.187 | +67.9% |
| D position (cm-1) | 1348.5 | 1352.5 | -4.0 | - |
| G position (cm-1) | 1583.1 | 1586.1 | -2.9 | - |
| 2D position (cm-1) | 2694.7 | 2701.8 | -7.1 | - |

## Spectra analysed

- **Control (uncoated)** - `examples/synthetic/synthetic_control.txt`, 2401 points, 800-3200 cm-1
- **Treated (coating A)** - `examples/synthetic/synthetic_treated.txt`, 2401 points, 800-3200 cm-1
