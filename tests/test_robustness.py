"""
Does the loader survive the file formats instruments actually produce?

The README claims the loader auto-detects delimiters, skips instrument
preambles, handles descending wavenumbers and marks missing bands as absent
rather than guessing. Each claim gets a test here; a claim without one is
just a sentence.
"""

import os

import numpy as np
import pytest

import analyse_raman as ar
import synthetic as syn

ROWS = [(1000.0 + i, 100.0 + i * 0.5) for i in range(50)]


def write(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text)
    return str(p)


@pytest.mark.parametrize("sep", ["\t", ",", ";", "   ", " "])
def test_every_delimiter_parses_identically(tmp_path, sep):
    """Tab, comma, semicolon and whitespace must all give the same numbers."""
    text = "".join(f"{x}{sep}{y}\n" for x, y in ROWS)
    x, y = ar.load_spectrum(write(tmp_path, "s.txt", text))
    assert np.allclose(x, [r[0] for r in ROWS])
    assert np.allclose(y, [r[1] for r in ROWS])


def test_instrument_preamble_is_skipped(tmp_path):
    """Comment and metadata lines before the data must not become data points."""
    preamble = ("# Raman export\n"
                "# Laser: 532 nm\n"
                "Wavenumber\tIntensity\n"
                "\n")
    text = preamble + "".join(f"{x}\t{y}\n" for x, y in ROWS)
    x, y = ar.load_spectrum(write(tmp_path, "s.txt", text))
    assert len(x) == len(ROWS)
    assert x[0] == ROWS[0][0]


def test_descending_wavenumbers_are_sorted(tmp_path):
    """Some instruments export high-to-low; the result must not depend on it."""
    up = "".join(f"{x}\t{y}\n" for x, y in ROWS)
    down = "".join(f"{x}\t{y}\n" for x, y in reversed(ROWS))
    xu, yu = ar.load_spectrum(write(tmp_path, "up.txt", up))
    xd, yd = ar.load_spectrum(write(tmp_path, "down.txt", down))
    assert np.allclose(xu, xd)
    assert np.allclose(yu, yd)


def test_extra_columns_are_ignored(tmp_path):
    """A third column (error bars, a second detector) must not shift the data."""
    text = "".join(f"{x}\t{y}\t{y * 0.01}\n" for x, y in ROWS)
    x, y = ar.load_spectrum(write(tmp_path, "s.txt", text))
    assert np.allclose(y, [r[1] for r in ROWS])


def test_too_few_rows_raises_rather_than_returning_nonsense(tmp_path):
    """A truncated or wrong file must fail loudly, not return a bad answer."""
    with pytest.raises(ValueError, match="fewer than 10"):
        ar.load_spectrum(write(tmp_path, "s.txt", "1000\t5\n1001\t6\n"))


def test_cosmic_ray_spikes_do_not_move_the_answer():
    """A cosmic ray in the baseline shoulder must not shift the result.

    The shoulder anchor is a *median*, so a spike far taller than the signal
    is rejected by construction. This is why the committed control spectrum
    ships with two injected spikes -- the demo exercises the claim.
    """
    clean = {}
    spiked = {}
    for name, spec in ar.BANDS.items():
        x, y = syn.make_spectrum(syn.CONTROL, seed=2)
        base, sig = ar.linear_baseline(x, y, spec["region"])
        clean[name] = ar.extract_peak(x, y - base, name, sig).intensity

        x, y = syn.make_spectrum(syn.CONTROL, seed=2, spikes=(1100.0, 2200.0))
        base, sig = ar.linear_baseline(x, y, spec["region"])
        spiked[name] = ar.extract_peak(x, y - base, name, sig).intensity

    for name in clean:
        assert spiked[name] == pytest.approx(clean[name], rel=1e-6), (
            f"{name}: spikes moved the height from {clean[name]:.1f} "
            f"to {spiked[name]:.1f}")


def test_band_outside_the_measured_range_is_reported_absent(tmp_path):
    """A 2D band that was never measured must be 'not detected', not invented."""
    x, y = syn.make_spectrum(syn.TREATED, x_max=2000.0, seed=1)
    text = "".join(f"{a}\t{b}\n" for a, b in zip(x, y))
    path = write(tmp_path, "truncated_sample.txt", text)

    sample = ar.analyse(path, label="truncated")
    assert "D" in sample.peaks and "G" in sample.peaks
    assert "2D" not in sample.peaks, "invented a 2D band that was never measured"
    assert any("2D" in n for n in sample.notes), "absence was not reported"


def test_baseline_slope_does_not_change_the_ratio():
    """A sloping fluorescence background must be removed, not measured.

    Same peaks, wildly different background: the ratio must barely move.
    """
    flat = syn.make_spectrum(syn.TREATED, background=False, noise=0.0, seed=1)
    sloped = syn.make_spectrum(syn.TREATED, background=True, noise=0.0, seed=1)

    def id_ig(xy):
        x, y = xy
        h = {}
        for n in ("D", "G"):
            base, sig = ar.linear_baseline(x, y, ar.BANDS[n]["region"])
            h[n] = ar.extract_peak(x, y - base, n, sig).intensity
        return h["D"] / h["G"]

    assert id_ig(sloped) == pytest.approx(id_ig(flat), rel=0.02)


def _sample(label, fitted_by_band, role):
    """A Sample with controllable per-band fit status, for reporting tests."""
    s = ar.Sample(label=label, path=f"{label}.txt", role=role)
    heights = {"D": 700.0, "G": 1000.0, "2D": 400.0}
    for name, fitted in fitted_by_band.items():
        s.peaks[name] = ar.Peak(name, ar.BANDS[name]["nominal"], heights[name],
                                50.0, 1000.0, fitted, 10.0, 5.0)
    return s


def test_comparison_refuses_to_subtract_two_different_estimators():
    """A ratio fitted on one sample and peak-picked on the other is not a delta.

    When the Lorentzian fit is rejected for one sample, its ID/IG comes from
    the peak maximum while the other sample's comes from the fit. Those are
    different measurements; printing a tidy percentage difference between them
    invites exactly the wrong conclusion.
    """
    treated = _sample("Ch/Ex", {"D": False, "G": False}, "treated")
    control = _sample("Control", {"D": True, "G": True}, "control")
    table = ar.comparison([treated, control])

    assert "not comparable" in table
    assert "%" not in table.split("\\*")[0].split("| ID/IG")[1].split("\n")[0]
    assert "different measurements" in table


def test_comparison_reports_a_percentage_when_estimators_match():
    """The guard must not fire when both samples were measured the same way."""
    treated = _sample("Ch/Ex", {"D": True, "G": True}, "treated")
    control = _sample("Control", {"D": True, "G": True}, "control")
    table = ar.comparison([treated, control])

    assert "not comparable" not in table
    assert "%" in table


def test_renaming_a_sample_does_not_disable_the_comparison():
    """The delta table must key off the sample's role, not its display name.

    It used to filter on the literal word "treated" appearing in the label.
    Renaming the samples for publication silently emptied the whole table --
    no error, just a section that quietly stopped being produced.
    """
    treated = _sample("anything at all", {"D": True, "G": True}, "treated")
    control = _sample("something else", {"D": True, "G": True}, "control")
    assert ar.comparison([treated, control]).strip(), (
        "comparison table vanished when the samples were renamed")


@pytest.mark.parametrize("name,role", [
    ("chex_01.txt", "treated"), ("ch-ex_2.txt", "treated"),
    ("treated_01.txt", "treated"), ("coated_a.txt", "treated"),
    ("control_01.txt", "control"), ("ctrl_02.txt", "control"),
    ("reference.txt", "control"), ("sample_xyz.txt", None),
])
def test_role_is_inferred_from_the_filename(name, role):
    assert ar.guess_role(name) == role
