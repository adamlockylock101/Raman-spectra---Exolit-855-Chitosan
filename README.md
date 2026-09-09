# Measuring how a coating changed a carbon material, from a noisy 1-D signal

Fire a laser at a material and the scattered light comes back as a noisy
one-dimensional signal: intensity against frequency. Three bumps in that signal
— conventionally called the **D**, **G** and **2D** bands — encode how
disordered the material's carbon structure is. The height of the D bump divided
by the height of the G bump (**ID/IG**) is the standard measure of defect
density; the 2D bump over the G bump (**I2D/IG**) indicates how many atomic
layers thick it is.

I have two samples of laser-induced graphene: one treated with a
Ch/Ex coating, one untreated control. **I need one defensible
number saying whether the treatment changed the structure, and an honest error
bar on it.**

That turns out not to be a one-liner. The bumps sit on a large, curved,
sample-dependent background — the material fluoresces under the laser, and the
glow is often bigger than the signal. Cosmic rays strike the detector during the
exposure and leave spikes taller than any real peak. The D and G bands overlap,
so each one's tail inflates the other. "Measure the height of the peak" is
therefore a chain of judgement calls, and each one moves the answer by several
percent.

![Annotated Raman spectra: raw signal with fitted baseline, and the two samples
after background removal scaled so the G band equals 1.0](results/spectrum.png)

The lower panels are scaled so the G band equals exactly 1.0, which makes the
y-axis read directly as the ratio: the control's D band reaches **1.16** and the
treated sample's only **0.73**. Fewer defects after treatment.

## Try it in thirty seconds

```bash
pip install -r requirements.txt
python3 analyse_raman.py examples/synthetic/*.txt --outdir results
```

Two example spectra are committed, so this produces real output on a fresh
clone with nothing else to set up. `make demo`, `make figure` and `make test`
are shorthands for the same commands.

To analyse your own data, drop two-column `.txt` files into `data/` and run
`python3 analyse_raman.py`. The loader auto-detects the delimiter, skips
instrument preamble blocks, ignores extra columns and handles descending
wavenumber order — the two committed examples are deliberately written in
different awkward formats to exercise exactly that.

Output goes to `results/raman_results.md` (tables), `results/raman_results.csv`
(machine-readable) and stdout:

| Metric | Treated | Control | Δ | % change |
|---|---|---|---|---|
| ID/IG | 0.727 | 1.155 | −0.428 | −37.1% |
| I2D/IG | 0.471 | 0.281 | +0.190 | +67.6% |

## How I know the numbers are right

This is the part I would want to see if I were reviewing someone else's
measurement code, so it is the part that got the most attention.

Every number above comes from code that pattern-matches a plausible-looking
answer out of noisy data. Plausible is not the same as correct, and on real
spectra there is nothing to check the answer against. So the pipeline is scored
against spectra I built myself, from Lorentzian peaks whose heights I chose —
where the right answer is known exactly rather than assumed:

```bash
make test      # 58 tests, ~17 seconds
```

The suite measures four things:

| What it checks | Result |
|---|---|
| **Recovery** — are known ratios recovered? | ID/IG within **3%** (worst case measured: 0.9%), I2D/IG within **8%** (measured: 2.0%) |
| **Robustness** — do awkward file formats, sloping backgrounds and cosmic rays change the answer? | Spikes in the baseline shoulder move the result by **&lt;0.0001%** (the anchor is a median, so they are rejected by construction) |
| **Calibration** — does the reported ± actually cover the truth? | Quoted ± is **0.93–1.00×** the real run-to-run scatter |
| **Cross-method** — do two independent implementations agree? | On position yes, within 2.2 cm⁻¹. On height **no** — see below |

Three findings came out of building it, and they changed the code:

**The error bars were wrong by 5×.** The original propagated raw single-point
noise onto a height that came from a smoothed, ~120-point Lorentzian fit,
ignoring the averaging-down the fit performs. Measured over 300 noise
realisations, the quoted ±0.036 described a quantity that actually scattered by
0.0067, and the bar covered the true value 100% of the time where a 1σ bar
should cover ~68%. `curve_fit` had been returning the correct covariance all
along and the code was discarding it (`popt, _ =`). Using it brings coverage to
57–66%. The bars were too *wide*, which is the safe direction to be wrong in —
but a number whose stated precision is fictional is not a measurement.

**The two implementations disagree, and the ground truth says which to
believe.** This repository contains two independent analyses of the same
quantity: a straight-line baseline plus a Lorentzian fit, and an
asymmetric-least-squares (AsLS) baseline plus a smoothed maximum. Agreement
would have been weak evidence that neither was badly wrong. They don't agree,
which turned out to be far more useful. The AsLS baseline rides up underneath
broad bands and subtracts part of the peak along with the background, biasing
**every height 11–25% low**. ID/IG largely survives that (3–5% low) because D and
G share one background and the error cancels in the ratio; I2D/IG does not
(12–13% low) because the 2D band sits in a different spectral region where
nothing cancels. Worse, the AsLS I2D/IG swings from 0.361 to 0.451 — against a true
0.480 — purely on the `lam` stiffness parameter, which is hand-set and was never
calibrated. Conclusion, now enforced by tests: quote ID/IG from either method,
quote I2D/IG only from the Lorentzian fit, and never quote absolute heights from
the AsLS path.

**A sixth bug, found by the suite rather than by reading.** The despiker divides
by the median absolute deviation of the signal. On a flat stretch the MAD is
zero, every z-score becomes NaN, `abs(NaN) > 7` evaluates to False — and
despiking silently switches itself off. No exception, no warning, no spikes
removed. Five earlier bugs in the same file were found by review (see
[CASE_STUDY.md](CASE_STUDY.md)); this one needed a test.

## Reading the numbers

- **ID/IG** — defect density. Lower means fewer sp³/edge defects and larger
  in-plane crystallite size. This is the more trustworthy of the two: it is
  measured on two strong bands sharing a single background.
- **I2D/IG** — layer count. Roughly >2 for monolayer, ~1 for bilayer, <0.5 for
  multilayer/turbostratic graphene. Laser-induced graphene is typically
  multilayer and turbostratic, so values well below 1 are normal. Carries
  3–6× the scatter of ID/IG, because the 2D band is weaker, broader and sits on
  its own background.

Reported ratios are peak **heights**, the usual convention in the
laser-induced-graphene literature. Integrated-area ratios are given alongside;
they differ whenever the bands have different widths, so quote whichever your
comparison literature uses and say which.

## Honest limitations

- **Both ratios come from a single spot on each sample.** The ± here is
  measurement uncertainty on one spectrum. Real spot-to-spot heterogeneity
  across a sample is typically much larger. A claim about the treatment needs
  several spots per sample and a comparison of the spreads — this tool measures
  one spectrum well, it does not establish an effect on its own.
- **The demo data is synthetic.** The numbers in the table above are recovered
  from spectra I generated, not from experimental measurements. That is
  deliberate: it makes the repository self-contained and the ground truth exact.
  Real spectra go in `data/` and are analysed by the same code path.
- **A straight-line baseline is an approximation.** It leaves a small residual
  bias where the true fluorescence background is curved, which is why the I2D/IG
  tolerance (8%) is looser than the ID/IG one (3%). The tolerances are measured,
  not chosen.

## Method

| Band | Search window (cm⁻¹) | Nominal centre |
|---|---|---|
| D | 1250–1450 | 1350 |
| G | 1500–1650 | 1580 |
| 2D | 2550–2850 | 2700 |

Each region gets a local linear baseline anchored on the flat, band-free
shoulders either side of it — 1050–1180 and 1750–1900 for the D/G doublet (one
baseline for both, since they sit on a shared background), 2350–2500 and
2900–3050 for the 2D band. Anchors are the **median** of each shoulder window,
which is what makes cosmic-ray spikes harmless.

Peak position and height come from the maximum of the baseline-corrected,
Savitzky-Golay-smoothed trace; a Lorentzian fit over ±60 cm⁻¹ then refines
centre, height and FWHM. The raw maximum is kept if the fit fails or wanders
outside the window. Uncertainty on each height is the Lorentzian fit covariance,
falling back to shoulder noise if the fit fails. Bands outside the measured
range are reported as "not detected" rather than guessed at.

Sample names are inferred from the filename — one containing `chex`,
`ch-ex` or `coated` is labelled treated; one containing `ctrl`, `reference` or `control`
is labelled control. That is what drives the treated-vs-control table. Override
with `--label` in file order.

## What's in here

| Path | What it is |
|---|---|
| `analyse_raman.py` | **The analysis.** Linear baseline + Lorentzian fit; writes markdown and CSV |
| `Raman_analysis.py` | The alternative AsLS implementation, kept as an independent cross-check |
| `synthetic.py` | Spectrum generator — the single source of ground truth for the demo and the tests |
| `plot_spectrum.py` | Renders the annotated figure above |
| `make_demo_data.py` | Writes the committed example spectra |
| `tests/` | The eval harness: recovery, robustness, calibration, regressions, cross-method |
| `examples/synthetic/` | Two committed spectra plus `ground_truth.json` |
| `data/` | Where real spectra go |
| `CASE_STUDY.md` | How this was built with AI, and what that caught |
