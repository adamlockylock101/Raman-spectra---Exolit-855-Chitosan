# How this was built, and what that caught

*A short account of the working method, written because the method transfers
and the graphene does not.*

## The starting point

A working analysis script: ~320 lines, notebook-style, asymmetric-least-squares
baseline correction, `scipy.find_peaks`, an interactive Plotly figure with peak
labels and FWHM annotations. It ran. It produced figures that looked right. Its
numbers were plausible.

That is the dangerous state for measurement code. Nothing was obviously broken,
so nothing prompted a closer look — and on real spectra there is no ground truth
to check the output against. A wrong ID/IG of 0.69 looks exactly like a right
one.

## Round 1 — review, and what "working" was hiding

I put the script through a line-by-line review with Claude. Five defects, none
of which raised an exception on the author's own data:

| # | Defect | Why it was invisible |
|---|---|---|
| 1 | `update_layout(titlefont=…)` — removed in Plotly 5 | Raised `ValueError` only on newer Plotly; on the original machine it worked |
| 2 | `fixer` and `modified_z_score` defined but **never called** | Dead code. Spectra were never despiked; a cosmic ray inside the search window could be reported as a peak |
| 3 | `fixer` indexed its repair window against `len(y)`, but `spikes` is one shorter after `np.diff` | `IndexError`, but only for a spike in the last few points |
| 4 | Peak-assignment windows too narrow | A real band slightly off nominal went unassigned, and the ratio was then **silently not calculated** |
| 5 | `parse_two_column_txt` split on tabs only | A comma-delimited export reported "no data loaded" on a perfectly good file |

Bugs 2 and 4 are the interesting ones. Neither crashes. Both quietly change or
suppress the result. This is the failure mode that review catches and running
the code does not.

Every one of these now has a regression test in `tests/test_regressions.py`.
None of them had one at the time — which is the whole reason they survived.

## What I did not let it do

Two decisions in that round matter more to me than the fixes.

**I did not let it claim the algorithms as ours.** The baseline correction is
Eilers & Boelens (2005); the despiking is Whitaker & Hayes (2018). Both are now
credited by name, year and journal in a header comment, including the detail
that the sparse-matrix formulation follows a widely circulated community
implementation of the paper rather than being derived independently.

**I did not let it touch the code it had not validated.** The AsLS baseline
logic was the most intricate part of the file and the easiest to "improve"
plausibly. The commit says so explicitly: *"AsLS baseline logic left unchanged
pending separate validation."* Fixing what you can verify and stopping at the
boundary of what you can't is the discipline; an AI assistant will happily
refactor straight past that line if you let it.

That call was vindicated two rounds later. The AsLS baseline turned out to be
biased by 11–25% — and if I had let it be rewritten on vibes first, I would never
have known whether the bias was original or introduced.

## Round 2 — a second implementation, as a cross-check

Rather than keep improving one script, I had a second one written from the
specification: a different baseline strategy (straight line anchored on
band-free shoulders) and a different height estimator (Lorentzian fit rather
than smoothed maximum).

The point was not to replace the first. It was to have two independent paths to
the same physical number. Two implementations that agree are weak evidence that
neither is badly wrong; two that disagree tell you there is something to find.

## Round 3 — ground truth, and what it found

Neither implementation could be checked against real spectra, because real
spectra do not come with answers. So the answers had to be manufactured:
synthetic spectra assembled from Lorentzian peaks with heights I chose, plus a
curved fluorescence background, realistic noise and injected cosmic rays. The
correct ID/IG and I2D/IG are then known exactly, and both pipelines can be
scored against them.

58 tests. Three findings, in ascending order of how much they changed my view:

**1. A sixth bug, invisible to review.** `modified_z_score` divides by the median
absolute deviation. On a flat stretch the MAD is zero, every z-score becomes
NaN, `abs(NaN) > 7` is False, and despiking silently switches itself off. Five
rounds of human and AI reading had not caught it. A test did, immediately.

**2. The error bars were fiction.** The uncertainty was computed by propagating
raw single-point noise onto a height that came from a smoothed, ~120-point
Lorentzian fit. Over 300 noise realisations, the quoted ±0.036 described a
quantity that actually varied by 0.0067 — five times too wide, covering the true
value 100% of the time where a 1σ bar should cover about 68%. The correct
covariance was already being computed by `curve_fit` and thrown away on the
line `popt, _ = curve_fit(...)`. Using it brings coverage to 57–66%.

The bars were too wide, not too narrow, so no published number would have been
wrong. But "±0.036" was a claim about precision that was not true, and I would
rather find that in a test than in a reviewer's question.

**3. The two implementations disagree — and the ground truth adjudicates.** The
AsLS baseline rides up underneath broad bands and subtracts part of the peak
along with the background, biasing every height 11–25% low. ID/IG mostly survives
(3–5%) because D and G share one background and the error cancels in the ratio.
I2D/IG does not (12–13%) because the 2D band sits elsewhere with its own
background. And the AsLS I2D/IG moves from 0.361 to 0.451 — true value 0.480 —
on the `lam` stiffness parameter alone, which is hand-set and had never been
calibrated against anything.

That last one is the result I would not have got any other way. It is not
reachable by reading the code, and it is not reachable by comparing to real
data. It needs a known answer and two independent attempts at it.

## What I take from this

- **AI is fastest at the layer where I am slowest**: reading unfamiliar code
  adversarially, and writing the volume of test scaffolding that makes a claim
  checkable. Five silent bugs in an afternoon.
- **AI is least reliable exactly where it is most confident**: numerical code
  that produces plausible output. Every one of the three findings above came
  from ground truth, not from inspection — including two in code that AI had
  itself written and I had reviewed.
- **The useful unit of verification is a known answer, not a passing run.**
  "The script runs and the figure looks right" was true of the original for its
  entire working life, with dead despiking code and a suppressed ratio in it.
- **Scope discipline is a decision the human has to make.** The two moments I am
  most confident about are both refusals: crediting the published algorithms
  rather than absorbing them, and declining to rewrite the AsLS baseline before
  there was any way to tell whether the rewrite was an improvement.

Commit history, with the verification for each change recorded in the message:
`git log`.
