"""
Two independent implementations, measured against the same known answer.

This repository contains two separate analyses of the same physical quantity:

  * `analyse_raman.py`  -- straight-line baseline anchored on band-free
                           shoulders, then a Lorentzian fit for the height.
  * `Raman_analysis.py` -- asymmetric-least-squares (AsLS) baseline, then
                           the maximum of a Savitzky-Golay smoothed trace.

Agreement between them would be evidence neither is badly wrong. They do not
agree, and the ground truth says which one to believe. That is worth far more
than agreement would have been, and it is the reason both are kept here.

Findings, measured by the tests below:

  * Both methods locate the bands to within 2.2 cm-1 of each other and of
    truth. Peak *position* is robust to the choice of method.
  * Peak *height* is not. The AsLS baseline rides up underneath broad bands
    and subtracts part of the peak along with the background, biasing every
    height 11-25% low.
  * ID/IG partly survives that (3-5% low) because D and G share one background
    and the error largely cancels in the ratio.
  * I2D/IG does not survive it (12-13% low) because the 2D band sits in a
    different spectral region with its own background, so nothing cancels.
  * I2D/IG from the AsLS path also swings from 0.361 to 0.451 -- against a
    true 0.480 -- purely on the `lam` stiffness parameter, which is hand-set
    and was never calibrated against a known answer.

Conclusion: quote ID/IG from either method, quote I2D/IG from the Lorentzian
fit, and do not quote absolute heights from the AsLS path at all.
"""

import numpy as np
import pytest

import analyse_raman as ar
import Raman_analysis as RA
import synthetic as syn

pytestmark = pytest.mark.filterwarnings("ignore::RuntimeWarning")

SAMPLES = [("treated", syn.TREATED), ("control", syn.CONTROL)]


def heights_lorentzian(bands, seed=1):
    """Method A: linear baseline + Lorentzian fit (analyse_raman.py)."""
    x, y = syn.make_spectrum(bands, seed=seed)
    out = {}
    for name, spec in ar.BANDS.items():
        base, sigma = ar.linear_baseline(x, y, spec["region"])
        peak = ar.extract_peak(x, y - base, name, sigma)
        out[name] = (peak.position, peak.intensity)
    return out


def heights_asls(bands, seed=1, **kw):
    """Method B: AsLS baseline + smoothed maximum (Raman_analysis.py)."""
    x, y = syn.make_spectrum(bands, seed=seed)
    _, _, search_area, _, peaks, positions, _, _, _ = RA.full_analysis2(y, x, **kw)
    out = {}
    for i, pos in enumerate(positions):
        h = float(search_area[peaks[i]])
        for name, spec in ar.BANDS.items():
            lo, hi = spec["window"]
            if lo <= pos <= hi:
                out[name] = (float(pos), h)
    return out


@pytest.mark.parametrize("label,bands", SAMPLES)
def test_both_methods_find_all_three_bands(label, bands):
    a, b = heights_lorentzian(bands), heights_asls(bands)
    assert set(a) == {"D", "G", "2D"}, f"Lorentzian method missed {set(a)}"
    assert set(b) == {"D", "G", "2D"}, f"AsLS method missed {set(b)}"


@pytest.mark.parametrize("label,bands", SAMPLES)
def test_methods_agree_on_band_positions(label, bands):
    """Position is the robust quantity: the two methods must agree to 5 cm-1."""
    a, b = heights_lorentzian(bands), heights_asls(bands)
    for name in ("D", "G", "2D"):
        assert a[name][0] == pytest.approx(b[name][0], abs=5.0), (
            f"{label} {name}: Lorentzian {a[name][0]:.1f} vs "
            f"AsLS {b[name][0]:.1f} cm-1")


@pytest.mark.parametrize("label,bands", SAMPLES)
def test_asls_baseline_biases_every_height_low(label, bands):
    """The documented failure mode, pinned down so it cannot drift unnoticed.

    Measured range is 11% to 25% low; the test allows 5% to 30%. If this
    starts failing, the AsLS path has changed behaviour and the guidance in
    this module's docstring needs revisiting.
    """
    b = heights_asls(bands)
    for name, band in bands.items():
        err = (b[name][1] - band.height) / band.height * 100
        assert -30.0 < err < -5.0, (
            f"{label} {name}: AsLS height error {err:+.1f}%, expected -30%..-5%")


@pytest.mark.parametrize("label,bands", SAMPLES)
def test_lorentzian_method_is_the_more_accurate_of_the_two(label, bands):
    """The reason `analyse_raman.py` is the one the README tells you to run."""
    a, b = heights_lorentzian(bands), heights_asls(bands)
    for name, band in bands.items():
        err_a = abs(a[name][1] - band.height) / band.height
        err_b = abs(b[name][1] - band.height) / band.height
        assert err_a < err_b, (
            f"{label} {name}: Lorentzian off by {err_a:.1%}, "
            f"AsLS off by {err_b:.1%}")


@pytest.mark.parametrize("label,bands", SAMPLES)
def test_id_ig_survives_the_baseline_error_but_i2d_ig_does_not(label, bands):
    """The central finding: which ratio you can trust from which method.

    D and G share one background, so a baseline error largely cancels in
    ID/IG. The 2D band has its own background, so it does not cancel in
    I2D/IG. Same spectrum, same method, very different reliability.
    """
    b = heights_asls(bands)
    true = syn.true_ratios(bands)
    id_ig_err = abs(b["D"][1] / b["G"][1] - true["ID/IG"]) / true["ID/IG"]
    i2d_ig_err = abs(b["2D"][1] / b["G"][1] - true["I2D/IG"]) / true["I2D/IG"]

    assert id_ig_err < 0.08, f"{label}: AsLS ID/IG off by {id_ig_err:.1%}"
    assert i2d_ig_err > id_ig_err, (
        f"{label}: expected I2D/IG ({i2d_ig_err:.1%}) to be the less reliable "
        f"of the two, but ID/IG was off by {id_ig_err:.1%}")


def test_i2d_ig_from_asls_is_hostage_to_an_uncalibrated_parameter():
    """`lam` is hand-set. Show how much of the answer it decides.

    The default is 1e7. Sweeping it over the range a practitioner might
    plausibly choose moves I2D/IG from 0.361 to 0.451 against a true 0.480 --
    nearly a fifth of the value of the number itself. An uncalibrated `lam`
    is not a detail, it is the result.
    """
    true = syn.true_ratios(syn.TREATED)["I2D/IG"]
    got = {}
    for lam in (1e6, 1e7, 1e8):
        h = heights_asls(syn.TREATED, l=lam)
        if {"2D", "G"} <= h.keys():
            got[lam] = h["2D"][1] / h["G"][1]

    assert len(got) == 3, "the lam sweep did not produce three usable results"
    spread = (max(got.values()) - min(got.values())) / true
    assert spread > 0.12, (
        f"expected I2D/IG to be strongly lam-dependent, got a spread of "
        f"{spread:.1%} across {got}")
