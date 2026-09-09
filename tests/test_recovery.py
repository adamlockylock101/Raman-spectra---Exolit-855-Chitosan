"""
Does the analysis recover ratios we know to be true?

This is the test that decides whether any number this repository prints can be
believed. Spectra are built from Lorentzians whose heights are chosen by us, so
the correct ID/IG and I2D/IG are known exactly rather than argued for. Anything
the pipeline reports is then scored against that.
"""

import pytest

import analyse_raman as ar
import synthetic as syn

SAMPLES = [("treated", syn.TREATED), ("control", syn.CONTROL)]


def measure(bands, seed=1, **kw):
    """Run the full analysis path and return the peaks it found."""
    x, y = syn.make_spectrum(bands, seed=seed, **kw)
    peaks = {}
    for name, spec in ar.BANDS.items():
        base, sigma = ar.linear_baseline(x, y, spec["region"])
        peaks[name] = ar.extract_peak(x, y - base, name, sigma)
    return peaks


@pytest.mark.parametrize("label,bands", SAMPLES)
def test_band_positions_recovered(label, bands):
    """Peak centres must land within 2 cm-1 of where we put them."""
    peaks = measure(bands)
    for name, band in bands.items():
        assert peaks[name].position == pytest.approx(band.centre, abs=2.0), (
            f"{label} {name}: found {peaks[name].position:.1f}, "
            f"placed at {band.centre:.1f} cm-1")


@pytest.mark.parametrize("label,bands", SAMPLES)
def test_band_widths_recovered(label, bands):
    """Fitted FWHM must be within 10% of the true width."""
    peaks = measure(bands)
    for name, band in bands.items():
        assert peaks[name].fitted, f"{label} {name}: Lorentzian fit did not take"
        assert peaks[name].fwhm == pytest.approx(band.fwhm, rel=0.10), (
            f"{label} {name}: fitted FWHM {peaks[name].fwhm:.1f}, "
            f"true {band.fwhm:.1f} cm-1")


@pytest.mark.parametrize("label,bands", SAMPLES)
def test_id_ig_recovered(label, bands):
    """ID/IG -- the headline number -- within 3% of truth.

    D and G sit on one shared background, so the residual baseline error
    largely cancels in this ratio. It is the more trustworthy of the two.
    """
    peaks = measure(bands)
    got = peaks["D"].intensity / peaks["G"].intensity
    assert got == pytest.approx(syn.true_ratios(bands)["ID/IG"], rel=0.03)


@pytest.mark.parametrize("label,bands", SAMPLES)
def test_i2d_ig_recovered(label, bands):
    """I2D/IG within 8%.

    The looser tolerance is not slack, it is the measured truth: the 2D band
    is weaker and broader than D or G, and it sits in a different spectral
    region, so the straight-line baseline error does not cancel the way it
    does for ID/IG. test_uncertainty.py quantifies the spread.
    """
    peaks = measure(bands)
    got = peaks["2D"].intensity / peaks["G"].intensity
    assert got == pytest.approx(syn.true_ratios(bands)["I2D/IG"], rel=0.08)


def test_treatment_effect_has_the_right_sign_and_size():
    """The comparison the whole project exists to make.

    Treated has lower ID/IG (fewer defects) and higher I2D/IG (fewer layers)
    by construction. The analysis must reproduce both directions, and the
    size of the ID/IG change to within a few percent.
    """
    t = measure(syn.TREATED)
    c = measure(syn.CONTROL)
    id_ig_t = t["D"].intensity / t["G"].intensity
    id_ig_c = c["D"].intensity / c["G"].intensity
    i2d_ig_t = t["2D"].intensity / t["G"].intensity
    i2d_ig_c = c["2D"].intensity / c["G"].intensity

    assert id_ig_t < id_ig_c, "treated should show fewer defects than control"
    assert i2d_ig_t > i2d_ig_c, "treated should show fewer layers than control"

    true_change = (syn.true_ratios(syn.TREATED)["ID/IG"]
                   / syn.true_ratios(syn.CONTROL)["ID/IG"] - 1) * 100
    got_change = (id_ig_t / id_ig_c - 1) * 100
    assert got_change == pytest.approx(true_change, abs=3.0), (
        f"ID/IG change {got_change:+.1f}%, true {true_change:+.1f}%")


def test_area_ratios_differ_from_height_ratios_as_expected():
    """Area and height ratios are not interchangeable, and must not be conflated.

    They agree only when the two bands have equal width. Treated D is wider
    than G, so the area ratio must come out above the height ratio -- if these
    two ever coincide, the area calculation has silently stopped working.
    """
    peaks = measure(syn.TREATED)
    height_ratio = peaks["D"].intensity / peaks["G"].intensity
    area_ratio = peaks["D"].area / peaks["G"].area
    assert syn.TREATED["D"].fwhm > syn.TREATED["G"].fwhm
    assert area_ratio > height_ratio
