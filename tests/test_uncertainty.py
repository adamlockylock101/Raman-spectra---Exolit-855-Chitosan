"""
Is the reported +/- honest?

A number with an error bar is a claim about how often the truth falls inside
that bar. A well-calibrated 1-sigma bar covers the true value about 68% of the
time. These tests measure that coverage over many noise realisations, which is
the only way to find out.

This suite caught a real defect. The original code estimated uncertainty by
propagating raw single-point noise onto a height that came from a smoothed,
~120-point Lorentzian fit -- ignoring the averaging-down that the fit performs.
The result was an error bar roughly 5x too wide, covering the truth 100% of the
time. `curve_fit` had been returning the correct covariance all along and the
code was discarding it. The tests below hold the fix in place.
"""

import numpy as np
import pytest

import analyse_raman as ar
import synthetic as syn

N_TRIALS = 120
CASES = [
    ("treated", syn.TREATED, "D", "G"),
    ("treated", syn.TREATED, "2D", "G"),
    ("control", syn.CONTROL, "D", "G"),
    ("control", syn.CONTROL, "2D", "G"),
]


def ratio_with_error(bands, seed, num, den):
    x, y = syn.make_spectrum(bands, seed=seed)
    pk = {}
    for n in (num, den):
        base, sigma = ar.linear_baseline(x, y, ar.BANDS[n]["region"])
        pk[n] = ar.extract_peak(x, y - base, n, sigma)
    a, b = pk[num], pk[den]
    r = a.intensity / b.intensity
    dr = abs(r) * np.hypot(a.height_err / a.intensity, b.height_err / b.intensity)
    return r, dr


def trials(bands, num, den):
    out = [ratio_with_error(bands, s, num, den) for s in range(N_TRIALS)]
    return np.array([r for r, _ in out]), np.array([d for _, d in out])


@pytest.mark.parametrize("label,bands,num,den", CASES)
def test_error_bar_matches_the_observed_scatter(label, bands, num, den):
    """The quoted +/- must be the same size as the real run-to-run spread.

    Within a factor of two either way. Wider than that and the number is
    uselessly vague; narrower and it is overconfident, which is worse.
    """
    ratios, errs = trials(bands, num, den)
    observed = ratios.std()
    reported = errs.mean()
    assert 0.5 < reported / observed < 2.0, (
        f"{label} {num}/{den}: reports +/-{reported:.4f} but actually "
        f"scatters by {observed:.4f} ({reported / observed:.1f}x)")


@pytest.mark.parametrize("label,bands", [("treated", syn.TREATED),
                                         ("control", syn.CONTROL)])
def test_baseline_systematic_stays_within_its_documented_bound(label, bands):
    """The reported +/- is statistical only. This bounds what it leaves out.

    A straight-line baseline cannot follow a curved fluorescence background,
    so every height carries a small systematic bias on top of the noise. That
    bias does not average away over repeated measurements, so it cannot be
    folded into a noise error bar -- and it is not a fixed fraction either
    (well under 1% on the strong D band, several percent on a weak, broad 2D
    band), so a single blanket percentage would overstate it for some bands
    and understate it for others.

    It is therefore reported separately and bounded here. If this fails, the
    baseline model has got worse and the README's stated tolerances are stale.
    """
    worst = 0.0
    for name in ("D", "G", "2D"):
        errs = []
        for seed in range(40):
            x, y = syn.make_spectrum(bands, seed=seed)
            base, sigma = ar.linear_baseline(x, y, ar.BANDS[name]["region"])
            peak = (ar.extract_peak(x, y - base, name, sigma) if name == "2D"
                    else ar.fit_dg(x, y - base, sigma)[name])
            errs.append(peak.intensity - bands[name].height)
        bias = abs(float(np.mean(errs))) / bands[name].height
        worst = max(worst, bias)
        assert bias < 0.06, (
            f"{label} {name}: baseline bias {bias:.1%} of band height, "
            f"above the 6% bound")
    assert worst > 0.0


@pytest.mark.parametrize("label,bands,num,den", CASES)
def test_statistical_bar_is_not_overconfident_about_noise(label, bands, num, den):
    """What the +/- does claim: it describes the run-to-run noise scatter.

    Measured against the mean of the repeated runs rather than the true value,
    so this tests the noise model alone -- the systematic is bounded by the
    test above.
    """
    ratios, errs = trials(bands, num, den)
    coverage = (np.abs(ratios - ratios.mean()) <= 2 * errs).mean()
    assert coverage >= 0.80, (
        f"{label} {num}/{den}: 2-sigma bar covered the run-to-run spread only "
        f"{coverage:.0%} of the time")


def test_id_ig_is_more_precise_than_i2d_ig():
    """States the pecking order between the two headline numbers.

    ID/IG is measured on two strong bands sharing one background; I2D/IG uses
    the weak, broad 2D band in a separate region. If a conclusion needs a
    precise number, it should rest on ID/IG. This ordering is a property of
    the physics, so it should hold for any sample.
    """
    for label, bands in (("treated", syn.TREATED), ("control", syn.CONTROL)):
        id_ig, _ = trials(bands, "D", "G")
        i2d_ig, _ = trials(bands, "2D", "G")
        rel_id = id_ig.std() / id_ig.mean()
        rel_2d = i2d_ig.std() / i2d_ig.mean()
        assert rel_id < rel_2d, (
            f"{label}: ID/IG scatter {rel_id:.1%} should be tighter than "
            f"I2D/IG {rel_2d:.1%}")


def test_noisier_spectra_get_bigger_error_bars():
    """A basic sanity property: more noise in, more uncertainty out."""
    def mean_err(noise):
        errs = []
        for seed in range(30):
            x, y = syn.make_spectrum(syn.TREATED, seed=seed, noise=noise)
            base, sigma = ar.linear_baseline(x, y, "DG")
            errs.append(ar.extract_peak(x, y - base, "G", sigma).height_err)
        return float(np.mean(errs))

    quiet, loud = mean_err(5.0), mean_err(40.0)
    assert loud > quiet * 2, (
        f"8x the noise moved the error bar from {quiet:.2f} to {loud:.2f}")
