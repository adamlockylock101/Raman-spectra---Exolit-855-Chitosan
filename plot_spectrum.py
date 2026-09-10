#!/usr/bin/env python3
"""
Render the explanatory figure used in the README.

    python3 plot_spectrum.py [spectra...] [--out results/spectrum.png]

Defaults to the committed synthetic demo spectra. The figure is built to be
readable by someone who has never seen a Raman spectrum: the top panel shows
what the instrument actually produces and what has to be removed, and the
bottom panels are scaled so that the G band equals 1.0 -- which makes the y
axis literally the ratio being reported.
"""

from __future__ import annotations

import argparse
import glob
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec

import analyse_raman as ar

# Palette: categorical slots 1 and 2, validated for colour-vision deficiency
# separation against the chart surface. Text never wears a series colour.
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
INK_MUTED = "#8a8983"
GRID = "#e4e3df"
SERIES = ["#2a78d6", "#eb6834"]
BAND_TINT = "#f0efe9"


def corrected(x, y, region):
    base, sigma = ar.linear_baseline(x, y, region)
    return y - base, base, sigma


def band_heights(path):
    """Run the full analysis so the figure shows exactly what the tables report.

    Going through `analyse` rather than calling the fitter directly means the
    detection threshold and the fit-plausibility gate apply here too: a band
    the analysis refused to report does not quietly reappear in the picture.
    """
    sample = ar.analyse(path)
    return sample.peaks, sample.notes


def style(ax):
    ax.set_facecolor(SURFACE)
    ax.grid(True, color=GRID, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=INK_2, labelsize=9, length=3, color=GRID)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="*")
    ap.add_argument("--out", default=os.path.join("results", "spectrum.png"))
    args = ap.parse_args(argv)

    # Prefer real spectra in data/; fall back to the committed synthetic pair.
    files = args.files or sorted(glob.glob(os.path.join("data", "*.txt"))) or \
        sorted(glob.glob(os.path.join("examples", "synthetic", "*.txt")))
    if not files:
        print("no spectra found")
        return 2

    samples = []
    for path in files:
        x, y = ar.load_spectrum(path)
        peaks, notes = band_heights(path)
        samples.append((ar.guess_label(path), x, y, peaks, notes))
    # Treated first, so the colour assignment is stable regardless of filename order.
    samples.sort(key=lambda s: "control" in s[0].lower())

    fig = plt.figure(figsize=(11, 7.4), facecolor=SURFACE)
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1.0, 1.15],
                  width_ratios=[1.55, 1.0], hspace=0.42, wspace=0.12,
                  left=0.07, right=0.985, top=0.87, bottom=0.09)

    # ---------------------------------------------------------------- top --
    # What the instrument produces, and what gets subtracted.
    ax = fig.add_subplot(gs[0, :])
    style(ax)
    label, x, y, _, _ = samples[0]
    ax.plot(x, y, color=SERIES[0], linewidth=1.4, zorder=3,
            label="measured spectrum")

    base_dg, _ = ar.linear_baseline(x, y, "DG")
    base_2d, _ = ar.linear_baseline(x, y, "2D")
    for region, base, (lo, hi) in (("DG", base_dg, (1050, 1900)),
                                   ("2D", base_2d, (2350, 3050))):
        m = (x >= lo) & (x <= hi)
        ax.plot(x[m], base[m], color=INK_2, linewidth=1.6, linestyle=(0, (5, 3)),
                zorder=4, label="fitted baseline" if region == "DG" else None)

    for name in ("D", "G", "2D"):
        lo, hi = ar.BANDS[name]["window"]
        ax.axvspan(lo, hi, color=BAND_TINT, zorder=1)
        # Anchored in axes fraction on y so the label sits inside the plot and
        # cannot collide with the subtitle above it.
        ax.text((lo + hi) / 2, 0.93, name, ha="center", va="top",
                transform=ax.get_xaxis_transform(), fontsize=10.5,
                color=INK_2, fontweight="bold", zorder=6)

    ax.set_xlim(max(x[0], 300.0), x[-1])
    ax.set_ylim(top=ax.get_ylim()[1] * 1.08)
    ax.set_ylabel("Intensity (counts)", fontsize=10, color=INK_2)
    ax.set_xlabel("Raman shift (cm$^{-1}$)", fontsize=10, color=INK_2)
    ax.set_title("1.  What the instrument measures",
                 fontsize=11.5, color=INK, loc="left", pad=26, fontweight="bold")
    ax.annotate("the broad slope is the sample glowing (fluorescence), not signal —\n"
                "a straight baseline is fitted across each shaded region and subtracted",
                xy=(0.0, 1.02), xycoords="axes fraction", fontsize=9.5,
                color=INK_2, va="bottom")
    ax.legend(loc="upper left", bbox_to_anchor=(0.60, 1.0), frameon=False,
              fontsize=9.5, labelcolor=INK_2, handlelength=1.8)

    # ------------------------------------------------------------- bottom --
    # Same data, background removed, scaled so G = 1. The y axis is the ratio.
    axes = [fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])]
    spans = [(1150.0, 1800.0), (2500.0, 2900.0)]
    regions = ["DG", "2D"]

    for ax_i, (ax, (lo, hi), region) in enumerate(zip(axes, spans, regions)):
        style(ax)
        for i, (label, x, y, peaks, _notes) in enumerate(samples):
            if "G" not in peaks:
                continue
            corr, _, _ = corrected(x, y, region)
            # Markers and the G normalisation both use the peak-maximum
            # estimator, so every number on the figure is the same measurement
            # made the same way -- and lands on the curve that is drawn.
            gbase, _ = ar.linear_baseline(x, y, ar.BANDS["G"]["region"])
            scale = ar.observed_max(x, y - gbase, "G")[1]
            m = (x >= lo) & (x <= hi)
            ax.plot(x[m], corr[m] / scale, color=SERIES[i], linewidth=1.8,
                    zorder=3, label=label)

            for name in ("D", "G", "2D"):
                if name not in peaks:
                    continue
                nbase, _ = ar.linear_baseline(x, y, ar.BANDS[name]["region"])
                mpos, mh = ar.observed_max(x, y - nbase, name)
                p = peaks[name]
                if not (lo <= mpos <= hi):
                    continue
                p = p.__class__(p.name, mpos, mh, p.fwhm, p.area, p.fitted,
                                p.sigma, p.height_err)
                h = mh / scale
                ax.plot([p.position, p.position], [0, h], color=SERIES[i],
                        linewidth=1.0, linestyle=":", alpha=0.75, zorder=2)
                ax.plot([p.position], [h], marker="o", markersize=5,
                        color=SERIES[i], markeredgecolor=SURFACE,
                        markeredgewidth=1.5, zorder=5)
                if name != "G":
                    # Opposite sides as well as opposite vertical offsets, so
                    # two samples' labels on the same band cannot collide.
                    offset = (8, 7) if i == 0 else (-8, -15)
                    ax.annotate(f"{h:.2f}", xy=(p.position, h),
                                xytext=offset, textcoords="offset points",
                                ha="left" if i == 0 else "right",
                                fontsize=9, color=INK_2, fontweight="bold",
                                zorder=6,
                                bbox=dict(facecolor=SURFACE, edgecolor="none",
                                          boxstyle="round,pad=0.15", alpha=0.9))

        missing = [lab for lab, _x, _y, pk, _n in samples
                   if not any(lo <= ar.BANDS[b]["nominal"] <= hi
                              and b in pk for b in ("D", "G", "2D"))]
        if missing and ax_i == 1:
            ax.text(0.5, 0.62, "not detected in\n" + "\n".join(missing),
                    transform=ax.transAxes, ha="center", va="center",
                    fontsize=9.5, color=INK_MUTED, style="italic")

        ax.axhline(1.0, color=INK_MUTED, linewidth=1.0, linestyle=(0, (2, 3)),
                   zorder=2)
        ax.set_xlim(lo, hi)
        ax.set_ylim(-0.06, 1.42)
        ax.set_xlabel("Raman shift (cm$^{-1}$)", fontsize=10, color=INK_2)

        for name in ("D", "G", "2D"):
            centre = ar.BANDS[name]["nominal"]
            if lo <= centre <= hi:
                ax.text(centre, 1.35, name, ha="center", fontsize=11,
                        color=INK, fontweight="bold")

    axes[0].set_ylabel("Intensity relative to the G band", fontsize=10, color=INK_2)
    axes[1].tick_params(labelleft=False)
    axes[1].text(0.99, 1.02, "G = 1.0", transform=axes[1].get_yaxis_transform(),
                 ha="right", va="bottom", fontsize=8.5, color=INK_MUTED,
                 zorder=6)
    axes[0].set_title("2.  Background removed, scaled so the G band = 1.0 — "
                      "the height of each peak is now the ratio itself",
                      fontsize=11.5, color=INK, loc="left", pad=10,
                      fontweight="bold")
    axes[0].legend(loc="upper left", frameon=False, fontsize=9.5,
                   labelcolor=INK_2, handlelength=1.8)

    # The headline comparison. Both samples are quoted from the same estimator
    # -- the baseline-corrected smoothed maximum -- because a Lorentzian fit
    # does not converge on every spectrum, and quoting a fitted ratio for one
    # sample against a peak-maximum ratio for the other compares two different
    # measurements and calls the difference a result.
    lines = []
    for label, x, y, peaks, _notes in samples:
        try:
            corr = {}
            for name in ("D", "G"):
                base, _ = ar.linear_baseline(x, y, ar.BANDS[name]["region"])
                corr[name] = ar.observed_max(x, y - base, name)[1]
            lines.append(f"{label}:  ID/IG = {corr['D'] / corr['G']:.2f}")
        except (ValueError, KeyError):
            pass
    if len(lines) == 2:
        fig.text(0.985, 0.005,
                 "   ·   ".join(lines) + "   (peak-maximum estimator, both samples)",
                 ha="right", va="bottom", fontsize=9, color=INK_MUTED)

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    fig.savefig(args.out, dpi=160, facecolor=SURFACE)
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
