# Raman analysis — LIG on CH/Ex vs Control

`analyse_raman.py` extracts the D, G and 2D bands from two-column Raman
spectra and reports the ID/IG and I2D/IG intensity ratios as a comparison
table.

## Usage

```bash
pip install numpy scipy
cp /path/to/*.txt data/            # one .txt per sample
python3 analyse_raman.py           # analyses every data/*.txt
```

Outputs `results/raman_results.md` (tables) and `results/raman_results.csv`
(machine-readable), and prints the report to stdout.

Sample names are inferred from the filename — a name containing `chex`,
`ch-ex`, `treated` or `coated` is labelled Ch/Ex, one containing `control`,
`ctrl` or `reference` is labelled Control, which is what drives the
Ch/Ex-vs-Control delta table. Override with `--label` in file order:

```bash
python3 analyse_raman.py a.txt b.txt --label "Ch/Ex" --label "Control"
```

## Input format

Two numeric columns, wavenumber then intensity. The loader auto-detects the
delimiter (tab, comma, semicolon, whitespace), skips headers and instrument
preamble blocks, ignores extra columns, and handles descending wavenumber
order. Bands outside the measured range are reported as "not detected"
rather than guessed at.

## Bands and method

| Band | Search window (cm⁻¹) | Nominal centre |
|---|---|---|
| D | 1250–1450 | 1350 |
| G | 1500–1650 | 1580 |
| 2D | 2550–2850 | 2700 |

Each region gets a local linear baseline anchored on the flat shoulders
either side of it — 1050–1180 and 1750–1900 for the D/G doublet (one
baseline for both, since they sit on a shared background), 2350–2500 and
2900–3050 for the 2D band. Peak position and height come from the maximum
of the baseline-corrected, Savitzky-Golay-smoothed trace, then a Lorentzian
fit over ±60 cm⁻¹ refines centre, height and FWHM; the raw maximum is kept
if the fit fails or wanders outside the window.

Reported ratios are peak *heights*, the usual convention for ID/IG and
I2D/IG in the LIG literature. Integrated-area ratios are given alongside —
they differ whenever the bands have different widths, so quote whichever
your comparison literature uses and say which. The ± figure propagates a
robust noise estimate taken from the baseline shoulders; it is measurement
noise only and does not cover spot-to-spot variation across the sample.

## Reading the numbers

- **ID/IG** — defect density. Lower means fewer sp³/edge defects and larger
  in-plane crystallite size.
- **I2D/IG** — layer count. Roughly >2 for monolayer, ~1 for bilayer,
  <0.5 for multilayer/turbostratic graphene. LIG is typically multilayer and
  turbostratic, so values well below 1 are normal.

Both are ratios from a single spot. For a claim about the treatment, average
several spots per sample and compare the spread — the noise-based ± here is
much smaller than real sample heterogeneity.

## Validation

The extraction was checked against synthetic spectra built from known
Lorentzians on a sloping background: input ratios of 0.750/1.250 (ID/IG) and
0.583/0.333 (I2D/IG) were recovered as 0.746/1.249 and 0.579/0.338, with
positions and FWHM recovered to better than 1 cm⁻¹.
