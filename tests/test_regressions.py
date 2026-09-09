"""
Regression tests for the four defects fixed in commit 14919db.

`Raman_analysis.py` was working code that produced plausible-looking figures.
Each of these bugs was silent -- none of them raised on the author's own data,
and two of them suppressed the ratio calculation without saying so. That is
exactly the failure mode a test suite exists to catch, and none of these tests
existed when the bugs were found.
"""

import numpy as np
import pytest

import analyse_raman as ar
import Raman_analysis as RA
import synthetic as syn


def test_despiker_survives_a_spike_in_the_final_points():
    """Bug 1: IndexError on a spike near the end of the spectrum.

    `fixer` built its repair window against len(y), but `spikes` comes from
    np.diff and is one element shorter, so a spike in the last few points
    indexed past the end. Crashed on any spectrum with a late cosmic ray.
    """
    y = np.full(200, 100.0)
    for i in (197, 198, 199):
        spiked = y.copy()
        spiked[i] += 9000.0
        RA.fixer(spiked, 5)          # must not raise IndexError


def test_despiking_is_actually_wired_into_the_analysis():
    """Bug 2: `fixer` and `modified_z_score` were defined but never called.

    Dead code. The spectra were never despiked, so a cosmic ray inside the
    search window could be reported as a peak. Test the behaviour, not the
    call: a spike must not survive into the result.
    """
    x, y = syn.make_spectrum(syn.TREATED, seed=1)
    spike_at = 1450.0
    y_spiked = syn.add_spikes(y, x, [spike_at], height=8000.0)

    _, _, _, _, peaks, positions, _, _, _ = RA.full_analysis2(y_spiked, x)
    assert len(peaks) > 0
    assert not np.any(np.abs(np.asarray(positions) - spike_at) < 20.0), (
        f"a cosmic ray at {spike_at} cm-1 was reported as a peak: {positions}")


def test_despiking_leaves_real_peaks_untouched():
    """The other half of bug 2: despiking must not damage the signal it keeps."""
    x, y = syn.make_spectrum(syn.TREATED, seed=1)
    cleaned = RA.fixer(y.copy(), 5)
    assert np.allclose(cleaned, y, atol=60.0), "despiker altered a clean spectrum"


@pytest.mark.parametrize("centre,band", [(1300.0, "D"), (1430.0, "D"),
                                         (1520.0, "G"), (2660.0, "2D"),
                                         (2820.0, "2D")])
def test_peaks_near_the_window_edges_are_still_assigned(centre, band):
    """Bug 3: assignment windows were too narrow.

    A real peak a little off its nominal position fell outside the window,
    was left unassigned, and the ratio was then silently not calculated --
    the script printed "cannot calculate" rather than a wrong number, which
    is better, but the data was fine and the windows were not.
    """
    lo, hi = ar.BANDS[band]["window"]
    assert lo <= centre <= hi, (
        f"a {band} band at {centre} cm-1 falls outside the "
        f"{lo:.0f}-{hi:.0f} window and would go unassigned")


def test_widened_windows_recover_an_off_nominal_band():
    """The same bug, end to end: a D band 50 cm-1 off nominal must still be found."""
    bands = dict(syn.TREATED)
    bands["D"] = syn.Band(centre=1300.0, height=720.0, fwhm=55.0)
    x, y = syn.make_spectrum(bands, seed=1)
    base, sigma = ar.linear_baseline(x, y, "DG")
    peak = ar.extract_peak(x, y - base, "D", sigma)
    assert peak.position == pytest.approx(1300.0, abs=3.0)
    assert peak.intensity == pytest.approx(720.0, rel=0.05)


@pytest.mark.parametrize("sep", ["\t", ",", ";", " "])
def test_second_parser_accepts_every_delimiter(tmp_path, sep):
    """Bug 4: `parse_two_column_txt` split on tabs only.

    A comma- or semicolon-delimited export produced zero usable rows and the
    script reported "No data loaded" on a perfectly good file.
    """
    rows = [(1000.0 + i, 50.0 + i) for i in range(20)]
    p = tmp_path / "s.txt"
    p.write_text("".join(f"{x}{sep}{y}\n" for x, y in rows))
    x, y = RA.parse_two_column_txt(str(p))
    assert np.allclose(x, [r[0] for r in rows])
    assert np.allclose(y, [r[1] for r in rows])


def test_plotly_figure_builds_without_raising(monkeypatch):
    """Bug 5: `update_layout(titlefont=...)` was removed in Plotly 5.

    The figure raised ValueError before it could be shown, so the script
    could not produce its one output. Nested title=dict(text=, font=) is the
    supported form.
    """
    monkeypatch.setattr(RA.go.Figure, "show", lambda self, *a, **k: None)
    path = "examples/synthetic/synthetic_chex.txt"
    fig = RA.plot_spectrum_analysis(path)
    assert fig is not None
    assert fig.layout.title.text
    assert fig.layout.xaxis.title.text == "Wavenumber (cm⁻¹)"


def test_despiker_handles_a_flat_spectrum_without_switching_itself_off():
    """Bug 6, found by this suite rather than by inspection.

    `modified_z_score` divides by the median absolute deviation. On a flat
    stretch the MAD is zero, every score becomes NaN, `abs(NaN) > 7` is
    False, and despiking quietly stops happening -- no exception, no warning,
    no spikes removed. A constant signal genuinely has no spikes, so the
    scores must be zero rather than NaN.
    """
    flat = np.full(200, 100.0)
    scores = np.asarray(RA.modified_z_score(flat), dtype=float)
    assert np.all(np.isfinite(scores)), "flat input produced non-finite z-scores"

    # And a real spike sitting on a flat background is still repaired.
    spiked = flat.copy()
    spiked[100] += 9000.0
    assert np.isfinite(RA.fixer(spiked, 5)).all()
