# Measuring how a coating changed a carbon material, from a noisy 1-D signal

Fire a laser at a material and the scattered light comes back as a noisy 1-D
signal: intensity against frequency. Three bumps in it — the **D**, **G** and
**2D** bands — encode how disordered the carbon structure is. D over G
(**ID/IG**) measures defect density; 2D over G (**I2D/IG**) indicates how many
atomic layers thick it is.

I have two samples of laser-induced graphene, one carrying a surface
coating and one untreated control. **I need one defensible
number saying whether the treatment changed the structure, and an honest error
bar on it.**

Measuring the height of a bump is not one line of code. The bumps sit on a
large, curved, sample-dependent background — the material fluoresces under the
laser, often more brightly than it scatters. Cosmic rays hit the detector
mid-exposure and leave spikes taller than any real peak. D and G overlap, so
each one's tail inflates the other.

![Annotated Raman spectra: raw signal with fitted baseline, and the two samples
after background removal scaled so the G band equals 1.0](results/spectrum.png)

The lower panels are scaled so the G band equals 1.0, making the y-axis read
directly as the ratio: the control's D band reaches **1.16**, the treated
sample's **0.73**. Fewer defects after treatment.

## Try it in thirty seconds

```bash
pip install -r requirements.txt
python3 analyse_raman.py            # analyses the two measured spectra in data/
```

Two measured spectra are committed, so this produces real output on a fresh
clone. `make demo`, `make figure` and `make test` are shorthands. Synthetic
spectra with known ground truth live in `examples/synthetic/` and are what the
test suite scores the pipeline against.

## What the measurement shows

| | Ch/Ex | Control |
|---|---|---|
| D band (defects) | 1372 cm⁻¹ | 1352 cm⁻¹ |
| G band (intact rings) | 1586 cm⁻¹ | 1582 cm⁻¹ |
| G band width (FWHM) | too broad to fit | 65 cm⁻¹ |
| **2D band (stacked layers)** | **not detected** (3.7σ) | **present** (10.5σ), I2D/IG = 0.37 |
| Valley between D and G | 0.48 × G height | 0.29 × G height |
| ID/IG | 0.66 | 0.77 |

Read the three bands as three questions about the carbon. **G** asks how much
intact ring structure is present. **D** asks how many defects and edges break
it up. **2D** asks whether the sheets are stacked in register — it only appears
when the material has real layered order.

**The robust result is the 2D band.** Control has an unmistakable one; Ch/Ex
has nothing above noise. Everything else agrees with that
reading: Ch/Ex's bands are so broad they cannot be fitted, and its
valley between D and G never drops — 0.48 × G height against Control's 0.29.
All three say the same thing. **Ch/Ex is substantially more amorphous;
Control retains genuine layered graphitic order.**

**The ID/IG difference is not a result.** It comes out 0.66 against 0.77, which
looks like a clean 15% drop. But Control's own ID/IG is 0.77 by peak
maximum and 0.69 by Lorentzian fit — an 11% swing from the choice of estimator
alone, on one spectrum. The gap between the samples is the same size as the gap
between two defensible ways of measuring either one of them, so this data
cannot support a claim about defect density in either direction.

Two things this cannot distinguish, and neither should be asserted from it:

- Whether the Ch/Ex **material** is more disordered, or whether the coating
  sits on top and contributes its own broad carbon signal. Both produce this
  spectrum.
- Whether either result is typical. This is **one spot on each sample**. Raman
  spot-to-spot variation across a real sample is routinely larger than the
  difference measured here, so several spots per sample are needed before any
  of this generalises.

What is safe to say: *these two spectra differ in a large and consistent way,
and the difference is in structural order rather than in defect count.*

## How I know the numbers are right

Real spectra don't come with answers, so the pipeline is scored against spectra
I built myself — Lorentzian peaks of heights I chose, on a curved background,
with noise and injected cosmic rays — where the right answer is known exactly:

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
57–66%. Too wide is the safe direction to be wrong in, but ±0.036 was still
a claim that wasn't true.

**The two implementations disagree, and the ground truth says which to
believe.** This repository contains two independent analyses of the same
quantity: a straight-line baseline plus a Lorentzian fit, and an
asymmetric-least-squares (AsLS) baseline plus a smoothed maximum. They don't
agree. The AsLS baseline rides up underneath
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

- **One spot per sample.** The ± is measurement uncertainty on a single
  spectrum. Spot-to-spot variation across a real sample is routinely larger.
  This tool measures one spectrum carefully; it does not establish an effect.
- **The quoted ± is statistical only.** It describes noise, which averages down
  over repeated measurements. The straight-line baseline adds a systematic on
  top that does not — measured on synthetic spectra at under 1% of band height
  for D, ~2% for G, and up to 5% for a weak, broad 2D band. The two are
  reported separately because they behave differently, not blended into one
  reassuring number.
- **The Ch/Ex spectrum cannot be peak-fitted.** Its bands are broad enough
  that a two-Lorentzian model rails against its own width bound, so the
  analysis rejects the fit and falls back to the peak maximum, and says so in
  the output notes. Heavily disordered carbon usually needs a 4–5 component
  model; choosing those components is a modelling decision this tool
  deliberately does not make on the user's behalf.
- **Ground truth is synthetic, and only synthetic.** Real spectra have no known
  answer, so the test suite scores the pipeline entirely against generated
  spectra. That validates the extraction maths; it does not validate that the
  band-assignment windows suit any particular material.

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

Sample names are inferred from the filename — one containing `treated` or
`coated` is labelled treated; one containing `control`, `ctrl` or `reference`
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
| `data/` | The two measured spectra analysed above (`chex_01`, `control_01`) |
| `CASE_STUDY.md` | How this was built with AI, and what that caught |
