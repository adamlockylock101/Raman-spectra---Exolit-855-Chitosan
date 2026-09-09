# How this was built, and what that caught

The analysis started as a ~320-line script: asymmetric-least-squares baseline
correction, `scipy.find_peaks`, an interactive Plotly figure. It ran, and its
numbers were plausible — which is the hard case, because a wrong ID/IG looks
exactly like a right one.

## Round 1 — review

A line-by-line review with Claude found five defects. None raised an exception
on the original data:

| # | Defect | Why it was invisible |
|---|---|---|
| 1 | `update_layout(titlefont=…)` — removed in Plotly 5 | Raised `ValueError` only on newer Plotly |
| 2 | `fixer` and `modified_z_score` defined but **never called** | Dead code; spectra were never despiked, so a cosmic ray could be reported as a peak |
| 3 | `fixer` indexed its repair window against `len(y)`, but `spikes` is one shorter after `np.diff` | `IndexError`, but only for a spike in the last few points |
| 4 | Peak-assignment windows too narrow | A band slightly off nominal went unassigned, and the ratio was then **silently not calculated** |
| 5 | `parse_two_column_txt` split on tabs only | A comma-delimited export reported "no data loaded" on a good file |

Two of them — the dead despiking code, and the windows that skipped the ratio —
changed the result without failing. Each defect now has a regression test in
`tests/test_regressions.py`.

## Two things I didn't let it do

**Claim the algorithms.** The baseline is Eilers & Boelens (2005), the despiking
Whitaker & Hayes (2018). Both are credited by name, year and journal, including
that the sparse-matrix formulation follows a community implementation rather
than being derived here.

**Rewrite what it hadn't validated.** The AsLS baseline was the most intricate
part of the file and the easiest to "improve" plausibly. The commit says so:
*"AsLS baseline logic left unchanged pending separate validation."* Two rounds
later it turned out to be biased 11–25% low. Had it been rewritten first, there
would have been no way to tell whether that bias was original or introduced.

## Round 2 — a second implementation

Rather than keep improving one script, I had a second written from the
specification, using a different baseline (straight line anchored on band-free
shoulders) and a different height estimator (Lorentzian fit, not smoothed
maximum). Two implementations that agree are weak evidence neither is badly
wrong; two that disagree point at something.

## Round 3 — ground truth

Real spectra don't come with answers, so I made spectra that do: Lorentzian
peaks of heights I chose, on a curved fluorescence background, with realistic
noise and injected cosmic rays. Both pipelines are scored against the known
ratio. 58 tests, three findings.

**A sixth bug.** `modified_z_score` divides by the median absolute deviation. On
a flat stretch that is zero, every z-score becomes NaN, `abs(NaN) > 7` is False,
and despiking silently switches itself off. Review had not caught it; a test did
immediately.

**The error bars were 5× too wide.** Uncertainty was propagated from raw
single-point noise onto a height taken from a smoothed, ~120-point fit. Over 300
noise realisations the quoted ±0.036 described a quantity that varied by 0.0067,
covering the true value 100% of the time where a 1σ bar should cover ~68%.
`curve_fit` had been returning the right covariance and the code discarded it
(`popt, _ = curve_fit(...)`). Coverage is now 57–66%. Too wide is the safe
direction to be wrong in, but ±0.036 was still a claim that wasn't true.

**The two implementations disagree.** The AsLS baseline rides up underneath
broad bands and subtracts part of the peak along with the background, biasing
every height 11–25% low. ID/IG mostly survives that (3–5%) because D and G share
one background and the error cancels in the ratio; I2D/IG does not (12–13%),
because the 2D band sits in a different region. Its I2D/IG also moves from 0.361
to 0.451 — true value 0.480 — on the hand-set `lam` parameter alone. Neither
reading the code nor comparing against real data would have shown this.

## What I'd carry over

- AI is fastest where I'm slowest: reading unfamiliar code adversarially, and
  writing enough test scaffolding to make a claim checkable.
- It is least reliable where its output looks most plausible. All three findings
  above came from ground truth rather than inspection — two of them in code AI
  wrote and I had already reviewed.
- "It runs and the figure looks right" was true of the original for its entire
  working life, with dead despiking code in it.
- Scope was mine to set, and the two calls I'm most confident about are both
  refusals.

`git log` records the verification for each change.
