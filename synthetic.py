#!/usr/bin/env python3
"""
Synthetic Raman spectra with known ground truth.

Every number this module puts into a spectrum is chosen by the caller, so the
"right answer" is known exactly rather than assumed. `analyse_raman.py` and
`Raman_analysis.py` are then scored against that answer in `tests/`, which is
the only reason to trust either of them on real data.

A spectrum is built as:

    measured = sum(Lorentzian peaks) + fluorescence background + noise + spikes

The background and the spikes exist to be removed again -- a pipeline that
only works on clean data is not a pipeline.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Band:
    """One Lorentzian band: where it sits, how tall it is, how wide."""
    centre: float      # cm^-1
    height: float      # counts, above background
    fwhm: float        # cm^-1

    @property
    def gamma(self) -> float:
        return self.fwhm / 2.0


def lorentzian(x: np.ndarray, band: Band) -> np.ndarray:
    return band.height * band.gamma**2 / ((x - band.centre) ** 2 + band.gamma**2)


def fluorescence(x: np.ndarray, offset: float = 400.0, slope: float = 0.15,
                 bump: float = 900.0, bump_centre: float = 1900.0,
                 bump_width: float = 700.0) -> np.ndarray:
    """A broad, smooth background of the kind a real sample fluoresces at.

    Deliberately not flat and not linear: a sloping bump is what makes
    baseline subtraction a real step rather than a constant offset.
    """
    return (offset
            + slope * (x - x[0])
            + bump * np.exp(-(((x - bump_centre) / bump_width) ** 2)))


def add_spikes(y: np.ndarray, x: np.ndarray, positions, height: float = 6000.0,
               width_pts: int = 1) -> np.ndarray:
    """Inject cosmic-ray spikes: a few adjacent points, enormously too high.

    Real detectors get hit by cosmic rays during an exposure. The result is a
    1-2 point spike far taller than any band, which will be mistaken for a
    peak by anything that just takes a maximum.
    """
    out = y.copy()
    for pos in positions:
        i = int(np.argmin(np.abs(x - pos)))
        out[i:i + width_pts] += height
    return out


def make_spectrum(bands: dict[str, Band], *, x_min: float = 800.0,
                  x_max: float = 3200.0, step: float = 1.0,
                  noise: float = 12.0, background: bool = True,
                  spikes=(), seed: int = 0
                  ) -> tuple[np.ndarray, np.ndarray]:
    """Build one spectrum from named bands. Returns (wavenumber, intensity)."""
    rng = np.random.default_rng(seed)
    x = np.arange(x_min, x_max + step, step)
    y = np.zeros_like(x)
    for band in bands.values():
        y += lorentzian(x, band)
    if background:
        y += fluorescence(x)
    if noise:
        y += rng.normal(0.0, noise, size=x.shape)
    if len(spikes):
        y = add_spikes(y, x, spikes)
    return x, y


# --------------------------------------------------------------------------
# The two reference samples used by the committed demo and by the tests.
#
# Heights are set so the ratios below are exact by construction. LIG is
# multilayer and turbostratic, so I2D/IG well under 1 is the expected regime.
# --------------------------------------------------------------------------

TREATED = {
    "D":  Band(centre=1348.0, height=720.0, fwhm=55.0),
    "G":  Band(centre=1583.0, height=1000.0, fwhm=45.0),
    "2D": Band(centre=2695.0, height=480.0, fwhm=90.0),
}

CONTROL = {
    "D":  Band(centre=1352.0, height=1150.0, fwhm=68.0),
    "G":  Band(centre=1586.0, height=1000.0, fwhm=52.0),
    "2D": Band(centre=2702.0, height=300.0, fwhm=105.0),
}


def true_ratios(bands: dict[str, Band]) -> dict[str, float]:
    """The height ratios that a correct analysis must recover."""
    return {
        "ID/IG": bands["D"].height / bands["G"].height,
        "I2D/IG": bands["2D"].height / bands["G"].height,
    }
