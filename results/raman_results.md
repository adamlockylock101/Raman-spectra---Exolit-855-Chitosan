# Raman peak analysis - D, G and 2D bands

## Peak positions and intensities

| Sample | Band | Position (cm-1) | Shift vs nominal | Intensity (counts) | FWHM (cm-1) | Area |
|---|---|---|---|---|---|---|
| Ch/Ex | D | 1348.5 | -1.5 | 724.9 | 54.9 | 56,050 |
| Ch/Ex | G | 1583.0 | +3.0 | 997.3 | 45.0 | 66,126 |
| Ch/Ex | 2D | 2694.8 | -5.2 | 469.4 | 88.5 | 37,667 |
| Control | D | 1352.6 | +2.6 | 1,156.6 | 69.0 | 99,257 |
| Control | G | 1585.7 | +5.7 | 1,001.2 | 52.8 | 73,890 |
| Control | 2D | 2701.7 | +1.7 | 281.2 | 98.3 | 21,797 |

## Intensity ratios

| Sample | ID/IG | I2D/IG | ID/IG (area) | I2D/IG (area) |
|---|---|---|---|---|
| Ch/Ex | 0.727 ± 0.007 | 0.471 ± 0.013 | 0.848 | 0.570 |
| Control | 1.155 ± 0.010 | 0.281 ± 0.013 | 1.343 | 0.295 |

## Treated vs control

| Metric | Treated | Control | Δ (treated − control) | % change |
|---|---|---|---|---|
| ID/IG | 0.727 | 1.155 | -0.428 | -37.1% |
| I2D/IG | 0.471 | 0.281 | +0.190 | +67.6% |
| D position (cm-1) | 1348.5 | 1352.6 | -4.0 | - |
| G position (cm-1) | 1583.0 | 1585.7 | -2.7 | - |
| 2D position (cm-1) | 2694.8 | 2701.7 | -6.9 | - |

## Spectra analysed

- **Ch/Ex** - `examples/synthetic/synthetic_chex.txt`, 2401 points, 800-3200 cm-1
- **Control** - `examples/synthetic/synthetic_control.txt`, 2401 points, 800-3200 cm-1
