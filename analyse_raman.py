#!/usr/bin/env python3
"""
Raman peak analysis for laser-induced graphene (LIG) spectra.

Extracts the D, G and 2D band positions and intensities from one or more
two-column Raman spectra (wavenumber, intensity) and reports the
ID/IG and I2D/IG intensity ratios.

Usage
-----
    python3 analyse_raman.py                     # analyses every .txt in ./data
    python3 analyse_raman.py a.txt b.txt         # explicit files
    python3 analyse_raman.py --outdir results    # where tables/plots are written

Method
------
1.  Load two numeric columns (delimiter auto-detected), sort by wavenumber.
2.  Subtract a *local linear baseline* per spectral region:
      - D+G region: anchored on flat shoulders either side of the D/G doublet
        (a single baseline for both bands, since they sit on one background).
      - 2D region:  anchored on the shoulders either side of the 2D band.
    Anchors are the median intensity of each shoulder window; if a shoulder
    lies outside the measured range it is clipped to the nearest available data.
3.  Smooth with a Savitzky-Golay filter (position finding only, order 3,
    window ~ 9 points) so a single noise spike cannot be mistaken for a peak.
4.  Peak position/intensity = maximum of the baseline-corrected, smoothed
    trace inside each search window. A Lorentzian is then fitted over the
    peak's local neighbourhood to refine centre, height and FWHM; if the fit
    fails or wanders outside the window, the raw maximum is kept.
5.  Ratios are peak-height ratios (the convention for ID/IG and I2D/IG in the
    LIG literature). Integrated-area ratios are reported alongside them.
6.  The uncertainty on each ratio comes from the Lorentzian fit covariance --
    the fit's own estimate of how well the amplitude is pinned down by the
    data. If the fit fails, it falls back to the (much more conservative)
    point noise measured on the baseline shoulders. Either way it is
    measurement uncertainty on one spot, not sample heterogeneity.
"""

from __future__ import annotations

import argparse
import csv
import glob
import math
import os
import re
import sys
from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter

# np.trapezoid is NumPy >= 2.0; np.trapz is the pre-2.0 spelling, removed in
# 2.x. Bind whichever this NumPy has so the script runs on either.
_trapezoid = getattr(np, "trapezoid", None) or np.trapz

# --------------------------------------------------------------------------
# Band definitions: search window, and the shoulders used to anchor baselines
# --------------------------------------------------------------------------

BANDS = {
    "D":  {"window": (1250.0, 1450.0), "nominal": 1350.0, "region": "DG"},
    "G":  {"window": (1500.0, 1650.0), "nominal": 1580.0, "region": "DG"},
    "2D": {"window": (2550.0, 2850.0), "nominal": 2700.0, "region": "2D"},
}

# Baseline anchor shoulders per region (flat, band-free stretches)
REGIONS = {
    "DG": {"left": (1050.0, 1180.0), "right": (1750.0, 1900.0)},
    "2D": {"left": (2350.0, 2500.0), "right": (2900.0, 3050.0)},
}

# The span the D and G bands are fitted over, together. It reaches from the end
# of the left baseline shoulder to the start of the right one.
DG_SPAN = (1180.0, 1750.0)

# A real Raman band in this system is tens of wavenumbers wide. Anything
# narrower is a noise spike or a cosmic ray, not a band.
MIN_FWHM = 15.0        # cm^-1
MAX_FWHM = 300.0       # cm^-1

# A band must stand this many noise sigma above the baseline to be reported at
# all. Below it, "not detected" is the honest answer -- fitting a Lorentzian to
# noise always succeeds and always returns a number.
DETECT_SNR = 5.0


# --------------------------------------------------------------------------
# Loading
# --------------------------------------------------------------------------

NUM = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")


def load_spectrum(path: str) -> tuple[np.ndarray, np.ndarray]:
    """Read a two-column Raman spectrum, tolerating headers and any delimiter."""
    xs: list[float] = []
    ys: list[float] = []
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line or line[0] in "#%;'\"":
                continue
            nums = NUM.findall(line.replace(",", " ") if line.count(",") and
                               not re.search(r"\d,\d{3}\b", line) else line)
            if len(nums) < 2:
                continue          # header / metadata line
            try:
                x, y = float(nums[0]), float(nums[1])
            except ValueError:
                continue
            if math.isfinite(x) and math.isfinite(y):
                xs.append(x)
                ys.append(y)
    if len(xs) < 10:
        raise ValueError(f"{path}: fewer than 10 usable data rows found")
    x = np.asarray(xs, float)
    y = np.asarray(ys, float)
    order = np.argsort(x)                      # some instruments write descending
    x, y = x[order], y[order]
    keep = np.concatenate(([True], np.diff(x) > 0))   # drop duplicate abscissae
    return x[keep], y[keep]


# --------------------------------------------------------------------------
# Baseline
# --------------------------------------------------------------------------

def _shoulder(x: np.ndarray, y: np.ndarray, lo: float, hi: float,
              fallback_edge: str) -> tuple[float, float, np.ndarray]:
    """Median (x, y) of a shoulder window, clipped to the measured range."""
    m = (x >= lo) & (x <= hi)
    if m.sum() < 3:
        # shoulder outside the measured range: fall back to the nearest 2% of data
        n = max(3, len(x) // 50)
        m = np.zeros_like(x, dtype=bool)
        if fallback_edge == "left":
            m[:n] = True
        else:
            m[-n:] = True
    return float(np.median(x[m])), float(np.median(y[m])), y[m]


def linear_baseline(x: np.ndarray, y: np.ndarray, region: str
                    ) -> tuple[np.ndarray, float]:
    """Straight baseline through the two shoulders of a region; also noise sigma."""
    lo_lo, lo_hi = REGIONS[region]["left"]
    hi_lo, hi_hi = REGIONS[region]["right"]
    x1, y1, s1 = _shoulder(x, y, lo_lo, lo_hi, "left")
    x2, y2, s2 = _shoulder(x, y, hi_lo, hi_hi, "right")
    if abs(x2 - x1) < 1e-9:
        slope = 0.0
    else:
        slope = (y2 - y1) / (x2 - x1)
    base = y1 + slope * (x - x1)
    # noise: scatter of the shoulders about their own medians, robust estimate
    resid = np.concatenate([s1 - np.median(s1), s2 - np.median(s2)])
    sigma = float(1.4826 * np.median(np.abs(resid))) or float(np.std(resid))
    return base, sigma


# --------------------------------------------------------------------------
# Peak extraction
# --------------------------------------------------------------------------

def lorentzian(x, amp, x0, gamma, offset):
    return offset + amp * (gamma ** 2) / ((x - x0) ** 2 + gamma ** 2)


def two_lorentzians(x, a1, x1, g1, a2, x2, g2, offset):
    """D and G as one model, sharing a single residual-background constant."""
    return (offset
            + a1 * g1 ** 2 / ((x - x1) ** 2 + g1 ** 2)
            + a2 * g2 ** 2 / ((x - x2) ** 2 + g2 ** 2))


@dataclass
class Peak:
    name: str
    position: float          # cm^-1
    intensity: float         # baseline-corrected height, counts
    fwhm: float | None       # cm^-1, from the Lorentzian fit
    area: float              # baseline-corrected integrated area over the window
    fitted: bool             # True if the Lorentzian fit was used
    sigma: float             # noise level in this region
    height_err: float        # 1-sigma uncertainty on `intensity`


def smooth(y: np.ndarray) -> np.ndarray:
    win = min(9, len(y) if len(y) % 2 else len(y) - 1)
    if win < 5:
        return y.copy()
    return savgol_filter(y, win, 3)


def observed_max(x: np.ndarray, y_corr: np.ndarray, name: str
                 ) -> tuple[float, float]:
    """Smoothed maximum inside a band window: (position, height)."""
    lo, hi = BANDS[name]["window"]
    m = (x >= lo) & (x <= hi)
    if m.sum() < 5:
        raise ValueError(f"band {name}: window {lo}-{hi} cm-1 not covered by the data")
    xw, ys = x[m], smooth(y_corr[m])
    i = int(np.argmax(ys))
    return float(xw[i]), float(ys[i])


def band_area(x: np.ndarray, y_corr: np.ndarray, name: str) -> float:
    lo, hi = BANDS[name]["window"]
    m = (x >= lo) & (x <= hi)
    return float(_trapezoid(np.clip(y_corr[m], 0, None), x[m]))


def at_bounds(popt, bounds, rtol: float = 1e-3) -> list[int]:
    """Indices of fitted parameters sitting on their own bound.

    A parameter pinned to a bound is not a measurement -- the optimiser wanted
    to go further and was not allowed to. The value it reports is the bound,
    not the data, so a fit with any railed parameter is rejected.
    """
    railed = []
    for i, (v, lo, hi) in enumerate(zip(popt, bounds[0], bounds[1])):
        span = abs(hi - lo) or 1.0
        if abs(v - lo) <= rtol * span or abs(v - hi) <= rtol * span:
            railed.append(i)
    return railed


def plausible(height: float, fwhm: float, obs_height: float) -> bool:
    """Is a fitted band physically believable, given what is in the data?

    Three ways a fit lies. It returns a height above the baseline larger than
    the largest value actually measured in the window (impossible). It returns
    a width narrower than any real band (it found a noise spike). Or it
    collapses to nothing. Any of these means fall back to the raw maximum.
    """
    if not (math.isfinite(height) and math.isfinite(fwhm)):
        return False
    if not (MIN_FWHM <= fwhm <= MAX_FWHM):
        return False
    return 0.0 < height <= 1.25 * obs_height


def fit_dg(x: np.ndarray, y_corr: np.ndarray, sigma: float) -> dict[str, Peak]:
    """Fit the D and G bands simultaneously.

    Fitting them one at a time, each with its own free offset, is degenerate
    whenever they overlap: the offset slides down to absorb the neighbour's
    tail and the amplitude climbs to compensate. On a disordered sample the
    valley between D and G never reaches the baseline, and a single-band fit
    then returns a D/G height ratio that is simply wrong -- in one measured
    spectrum it put the G height 60% above the largest value in the window.
    One model, two bands, one shared constant.
    """
    peaks: dict[str, Peak] = {}
    obs = {n: observed_max(x, y_corr, n) for n in ("D", "G")}
    areas = {n: band_area(x, y_corr, n) for n in ("D", "G")}

    lo, hi = DG_SPAN
    m = (x >= lo) & (x <= hi)
    xw, yw = x[m], y_corr[m]

    fitted: dict[str, tuple] = {}
    if m.sum() >= 12 and all(h > 0 for _, h in obs.values()):
        (dpos, dh), (gpos, gh) = obs["D"], obs["G"]
        p0 = [dh, dpos, 30.0, gh, gpos, 30.0, 0.0]
        big = 3.0 * max(dh, gh)
        bounds = (
            [0.0, *BANDS["D"]["window"][:1], MIN_FWHM / 2,
             0.0, *BANDS["G"]["window"][:1], MIN_FWHM / 2, -big],
            [big, BANDS["D"]["window"][1], MAX_FWHM / 2,
             big, BANDS["G"]["window"][1], MAX_FWHM / 2, big],
        )
        try:
            popt, pcov = curve_fit(two_lorentzians, xw, yw, p0=p0,
                                   bounds=bounds, maxfev=40000)
            if not at_bounds(popt, bounds):
                err = np.sqrt(np.diag(pcov))
                fitted["D"] = (popt[0], popt[1], 2 * popt[2], err[0])
                fitted["G"] = (popt[3], popt[4], 2 * popt[5], err[3])
        except (RuntimeError, ValueError):
            pass

    for name in ("D", "G"):
        pos, obs_h = obs[name]
        height, fwhm, ok, herr = obs_h, None, False, sigma
        if name in fitted:
            amp, centre, width, amp_err = fitted[name]
            if plausible(amp, width, obs_h):
                height, pos, fwhm, ok = float(amp), float(centre), float(width), True
                if math.isfinite(amp_err) and amp_err > 0:
                    herr = float(amp_err)
        peaks[name] = Peak(name, pos, height, fwhm, areas[name], ok, sigma, herr)
    return peaks


def extract_peak(x: np.ndarray, y_corr: np.ndarray, name: str,
                 sigma: float) -> Peak:
    lo, hi = BANDS[name]["window"]
    m = (x >= lo) & (x <= hi)
    if m.sum() < 5:
        raise ValueError(f"band {name}: window {lo}-{hi} cm-1 not covered by the data")
    xw, yw = x[m], y_corr[m]
    ys = smooth(yw)

    pos, height = observed_max(x, y_corr, name)
    obs_h = height
    area = band_area(x, y_corr, name)
    fwhm, fitted = None, False

    # Fallback uncertainty, used only if the Lorentzian fit fails or is
    # rejected: the raw point noise. Deliberately conservative -- the height
    # comes from a smoothed trace, so the true uncertainty is smaller (see
    # tests/test_uncertainty.py, which measures the coverage of both).
    height_err = sigma

    # Refine with a Lorentzian. The neighbourhood scales with the band, since
    # a fixed window narrower than the FWHM cannot constrain a broad band.
    half = [i for i, v in enumerate(ys) if v >= height / 2.0]
    est_fwhm = float(xw[half[-1]] - xw[half[0]]) if len(half) >= 2 else 50.0
    span = float(np.clip(1.5 * est_fwhm, 60.0, 200.0))
    nb = (xw >= pos - span) & (xw <= pos + span)
    if nb.sum() >= 6 and height > 0:
        p0 = [height, pos, 25.0, 0.0]
        bounds = ([0, lo, MIN_FWHM / 2, -abs(height)],
                  [10 * height + 1e-9, hi, MAX_FWHM / 2, abs(height) + 1e-9])
        try:
            popt, pcov = curve_fit(lorentzian, xw[nb], yw[nb], p0=p0,
                                   bounds=bounds, maxfev=20000)
            amp, x0, gamma, _off = popt
            if (lo <= x0 <= hi and plausible(amp, 2 * gamma, obs_h)
                    and not at_bounds(popt, bounds)):
                pos, height, fwhm, fitted = float(x0), float(amp), float(2 * gamma), True
                # The fit's own covariance is the honest uncertainty on the
                # amplitude: it is estimated from the residuals over every
                # point in the fit window, so it already accounts for the
                # averaging-down that propagating single-point noise ignores.
                var = float(pcov[0, 0])
                if math.isfinite(var) and var > 0:
                    height_err = math.sqrt(var)
        except (RuntimeError, ValueError):
            pass

    return Peak(name, pos, height, fwhm, area, fitted, sigma, height_err)


# --------------------------------------------------------------------------
# Per-sample analysis
# --------------------------------------------------------------------------

@dataclass
class Sample:
    label: str
    path: str
    peaks: dict[str, Peak] = field(default_factory=dict)
    xrange: tuple[float, float] = (0.0, 0.0)
    npts: int = 0
    notes: list[str] = field(default_factory=list)

    def ratio(self, num: str, den: str) -> tuple[float, float]:
        """Intensity ratio and its 1-sigma uncertainty from propagated noise."""
        a, b = self.peaks[num], self.peaks[den]
        r = a.intensity / b.intensity
        dr = abs(r) * math.hypot(a.height_err / a.intensity,
                                 b.height_err / b.intensity)
        return r, dr

    def area_ratio(self, num: str, den: str) -> float:
        return self.peaks[num].area / self.peaks[den].area


def guess_label(path: str) -> str:
    stem = os.path.splitext(os.path.basename(path))[0]
    low = stem.lower()
    if "treated" in low or "coated" in low:
        return "Treated (coating A)"
    if "control" in low or "ctrl" in low or "reference" in low:
        return "Control (uncoated)"
    return stem


def analyse(path: str, label: str | None = None) -> Sample:
    x, y = load_spectrum(path)
    s = Sample(label=label or guess_label(path), path=path,
               xrange=(float(x[0]), float(x[-1])), npts=len(x))

    baselines: dict[str, tuple[np.ndarray, float]] = {}
    for region in ("DG", "2D"):
        baselines[region] = linear_baseline(x, y, region)

    covered = []
    for name, spec in BANDS.items():
        lo, hi = spec["window"]
        if x[0] > lo or x[-1] < hi:
            s.notes.append(
                f"{name} window ({lo:.0f}-{hi:.0f}) only partly covered by data "
                f"({x[0]:.0f}-{x[-1]:.0f} cm-1)")
            if x[-1] < lo or x[0] > hi:
                s.notes.append(f"{name} band absent from the measured range - skipped")
                continue
        covered.append(name)

    found: dict[str, Peak] = {}
    if {"D", "G"} <= set(covered):
        base, sigma = baselines["DG"]
        found.update(fit_dg(x, y - base, sigma))
    else:
        for name in covered:
            if BANDS[name]["region"] == "DG":
                base, sigma = baselines["DG"]
                found[name] = extract_peak(x, y - base, name, sigma)
    if "2D" in covered:
        base, sigma = baselines["2D"]
        found["2D"] = extract_peak(x, y - base, "2D", sigma)

    # Detection gate. Fitting a Lorentzian to noise always succeeds and always
    # returns a number; a band that does not stand clear of the noise has to be
    # reported as absent instead, or the ratio is invented rather than measured.
    for name, peak in found.items():
        _, obs_h = observed_max(x, y - baselines[BANDS[name]["region"]][0], name)
        snr = obs_h / peak.sigma if peak.sigma > 0 else float("inf")
        if snr < DETECT_SNR:
            s.notes.append(
                f"{name} band not detected - peak stands only {snr:.1f} sigma above "
                f"the noise, below the {DETECT_SNR:.0f} sigma threshold")
            continue
        if not peak.fitted:
            s.notes.append(
                f"{name} band: Lorentzian fit rejected as implausible, "
                f"height taken from the smoothed maximum instead")
        s.peaks[name] = peak

    for region, (lo_hi) in REGIONS.items():
        for side, (a, b) in lo_hi.items():
            if x[0] > b or x[-1] < a:
                s.notes.append(
                    f"{region} baseline {side} shoulder ({a:.0f}-{b:.0f}) outside "
                    f"the measured range - clipped to the nearest data")
    return s


# --------------------------------------------------------------------------
# Reporting
# --------------------------------------------------------------------------

def fmt(v, nd=1):
    return "n/a" if v is None else f"{v:,.{nd}f}"


def peak_table(samples: list[Sample]) -> str:
    rows = ["| Sample | Band | Position (cm-1) | Shift vs nominal | Intensity (counts) | FWHM (cm-1) | Area |",
            "|---|---|---|---|---|---|---|"]
    for s in samples:
        for name in ("D", "G", "2D"):
            p = s.peaks.get(name)
            if p is None:
                rows.append(f"| {s.label} | {name} | not detected | - | - | - | - |")
                continue
            d = p.position - BANDS[name]["nominal"]
            rows.append(f"| {s.label} | {name} | {p.position:.1f} | {d:+.1f} | "
                        f"{fmt(p.intensity)} | {fmt(p.fwhm)} | {fmt(p.area, 0)} |")
    return "\n".join(rows)


def ratio_table(samples: list[Sample]) -> str:
    rows = ["| Sample | ID/IG | I2D/IG | ID/IG (area) | I2D/IG (area) |",
            "|---|---|---|---|---|"]
    for s in samples:
        def cell(n, d):
            if n not in s.peaks or d not in s.peaks:
                return "n/a"
            r, dr = s.ratio(n, d)
            return f"{r:.3f} ± {dr:.3f}"

        def acell(n, d):
            if n not in s.peaks or d not in s.peaks:
                return "n/a"
            return f"{s.area_ratio(n, d):.3f}"

        rows.append(f"| {s.label} | {cell('D','G')} | {cell('2D','G')} | "
                    f"{acell('D','G')} | {acell('2D','G')} |")
    return "\n".join(rows)


def comparison(samples: list[Sample]) -> str:
    """Treated-vs-control deltas, when exactly one of each is identifiable."""
    treated = [s for s in samples if "treated" in s.label.lower()]
    control = [s for s in samples if "control" in s.label.lower()]
    if len(treated) != 1 or len(control) != 1:
        return ""
    t, c = treated[0], control[0]
    out = ["| Metric | Treated | Control | Δ (treated − control) | % change |",
           "|---|---|---|---|---|"]
    mixed = []
    for label, num, den in (("ID/IG", "D", "G"), ("I2D/IG", "2D", "G")):
        if not all(k in t.peaks and k in c.peaks for k in (num, den)):
            continue
        rt, _ = t.ratio(num, den)
        rc, _ = c.ratio(num, den)
        pct = (rt - rc) / rc * 100 if rc else float("nan")
        # A ratio measured by Lorentzian fit on one sample and by peak maximum
        # on the other is two different measurements; their difference is not
        # a result. Flag it rather than printing a clean-looking percentage.
        if any(t.peaks[k].fitted != c.peaks[k].fitted for k in (num, den)):
            mixed.append(label)
            out.append(f"| {label} | {rt:.3f} | {rc:.3f} | {rt - rc:+.3f} | "
                       f"not comparable * |")
        else:
            out.append(f"| {label} | {rt:.3f} | {rc:.3f} | {rt - rc:+.3f} | {pct:+.1f}% |")
    for name in ("D", "G", "2D"):
        if name in t.peaks and name in c.peaks:
            pt, pc = t.peaks[name].position, c.peaks[name].position
            out.append(f"| {name} position (cm-1) | {pt:.1f} | {pc:.1f} | {pt - pc:+.1f} | - |")
    if mixed:
        out += ["",
                f"\\* {', '.join(mixed)} was measured by Lorentzian fit on one "
                "sample and by peak maximum on the other, because the fit was "
                "rejected for the other. Those are different measurements, so "
                "the difference between them is not a result. Compare the "
                "samples using a single estimator instead."]
    return "\n".join(out)


def write_csv(samples: list[Sample], path: str) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["sample", "file", "band", "position_cm-1", "intensity_counts",
                    "fwhm_cm-1", "area", "lorentzian_fit", "noise_sigma",
                    "height_err", "ID/IG", "I2D/IG"])
        for s in samples:
            rdg = f"{s.ratio('D','G')[0]:.4f}" if {"D", "G"} <= s.peaks.keys() else ""
            r2g = f"{s.ratio('2D','G')[0]:.4f}" if {"2D", "G"} <= s.peaks.keys() else ""
            for name in ("D", "G", "2D"):
                p = s.peaks.get(name)
                if p is None:
                    w.writerow([s.label, s.path, name, "", "", "", "", "", "", "",
                                rdg, r2g])
                    continue
                w.writerow([s.label, s.path, name, f"{p.position:.2f}",
                            f"{p.intensity:.2f}",
                            f"{p.fwhm:.2f}" if p.fwhm else "",
                            f"{p.area:.1f}", p.fitted, f"{p.sigma:.2f}",
                            f"{p.height_err:.2f}", rdg, r2g])


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="*", help="spectra (.txt); default: data/*.txt")
    ap.add_argument("--label", action="append", default=[],
                    help="override sample name, in file order (repeatable)")
    ap.add_argument("--outdir", default="results", help="output directory")
    args = ap.parse_args(argv)

    files = args.files or sorted(glob.glob(os.path.join("data", "*.txt")))
    if not files:
        print("No spectra found. Put the .txt files in ./data or pass them as "
              "arguments.", file=sys.stderr)
        return 2

    samples = []
    for i, f in enumerate(files):
        label = args.label[i] if i < len(args.label) else None
        try:
            samples.append(analyse(f, label))
        except Exception as exc:                      # noqa: BLE001
            print(f"ERROR reading {f}: {exc}", file=sys.stderr)
    if not samples:
        return 1

    os.makedirs(args.outdir, exist_ok=True)
    parts = ["# Raman peak analysis - D, G and 2D bands", "",
             "## Peak positions and intensities", "", peak_table(samples), "",
             "## Intensity ratios", "", ratio_table(samples), ""]
    cmp_tbl = comparison(samples)
    if cmp_tbl:
        parts += ["## Treated vs control", "", cmp_tbl, ""]
    parts += ["## Spectra analysed", ""]
    for s in samples:
        parts.append(f"- **{s.label}** - `{s.path}`, {s.npts} points, "
                     f"{s.xrange[0]:.0f}-{s.xrange[1]:.0f} cm-1")
        for n in s.notes:
            parts.append(f"  - note: {n}")
    report = "\n".join(parts) + "\n"

    md_path = os.path.join(args.outdir, "raman_results.md")
    csv_path = os.path.join(args.outdir, "raman_results.csv")
    with open(md_path, "w") as fh:
        fh.write(report)
    write_csv(samples, csv_path)
    print(report)
    print(f"Written: {md_path}\n         {csv_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
