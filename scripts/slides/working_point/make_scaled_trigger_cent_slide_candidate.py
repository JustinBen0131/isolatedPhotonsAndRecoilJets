#!/usr/bin/env python3
"""Build a slide-ready PNG for the scaled-trigger centrality campaign."""

from __future__ import annotations
# Keep purpose-folder helpers runnable when invoked directly.
import sys as _codex_sys
from pathlib import Path as _CodexPath
_CODEX_THIS_FILE = _CodexPath(__file__).resolve()
_CODEX_SCRIPTS_DIR = next((p for p in _CODEX_THIS_FILE.parents if p.name == "scripts"), _CODEX_THIS_FILE.parent)
_CODEX_SCRIPTS_DIR_STR = str(_CODEX_SCRIPTS_DIR)
if _CODEX_SCRIPTS_DIR_STR not in _codex_sys.path:
    _codex_sys.path.append(_CODEX_SCRIPTS_DIR_STR)
del _CODEX_THIS_FILE, _CODEX_SCRIPTS_DIR, _CODEX_SCRIPTS_DIR_STR

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import ROOT  # type: ignore  # noqa: E402
from matplotlib.patches import FancyBboxPatch  # noqa: E402

from make_scaled_trigger_cent_summary_panels import (  # noqa: E402
    DEFAULT_INPUT,
    DEFAULT_OUTDIR,
    TRIGGERS,
    Hist,
    integral,
    ratio_hist,
    read_hist,
)


OUTDIR = DEFAULT_OUTDIR / "slide_candidates"
OUTPNG = OUTDIR / "scaled_trigger_cent_bins_after_slide4_candidate.png"
SUMMARY = DEFAULT_OUTDIR / "scaledTriggerCentStudy_centrality_summary.csv"
CENT_BINS = [("cent0_20", "0-20%"), ("cent20_50", "20-50%"), ("cent50_80", "50-80%")]


def load_summary() -> dict[str, dict[str, str]]:
    if not SUMMARY.exists():
        return {}
    with SUMMARY.open() as f:
        return {row["suffix"]: row for row in csv.DictReader(f)}


def setup_axis(ax: plt.Axes) -> None:
    ax.tick_params(direction="in", top=True, right=True, labelsize=8, length=4)
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)
    ax.grid(True, which="major", color="#d9dee7", linewidth=0.5, alpha=0.75)


def draw_row(
    fig: plt.Figure,
    grid,
    row: int,
    label: str,
    suffix: str,
    hists: dict[str, Hist],
    summary: dict[str, str],
) -> None:
    ax_l = fig.add_subplot(grid[row, 0])
    ax_r = fig.add_subplot(grid[row, 1])

    for key in ("mbd", "p10", "p12"):
        _, _, legend, color = TRIGGERS[key]
        hist = hists[key]
        lw = 1.55 if key == "mbd" else 1.35
        ax_l.step(hist.edges[:-1], hist.values, where="post", lw=lw, color=color, label=legend)

    positives = np.concatenate([hist.values[hist.values > 0] for hist in hists.values()])
    ymin = max(3.0e4, float(np.nanmin(positives)) * 0.5) if positives.size else 1.0
    ymax = max(float(max(np.nanmax(hist.values) for hist in hists.values())) * 1.8, ymin * 10.0)
    ax_l.set_yscale("log")
    ax_l.set_xlim(1.0, 20.0)
    ax_l.set_ylim(ymin, ymax)
    ax_l.set_ylabel("counts", fontsize=9)
    ax_l.text(0.035, 0.91, rf"$\it{{sPHENIX}}$ Internal  Au+Au, $\sqrt{{s_{{NN}}}}=200$ GeV",
              transform=ax_l.transAxes, fontsize=8.5, va="top")
    ax_l.text(0.035, 0.78, f"Centrality {label}", transform=ax_l.transAxes, fontsize=9, va="top")
    if row == 0:
        ax_l.legend(frameon=False, fontsize=8, loc="upper right", handlelength=2.3)
        ax_l.set_title("Max-cluster energy overlay", fontsize=11, pad=6, fontweight="bold")
    if row == 2:
        ax_l.set_xlabel(r"Max cluster energy [GeV], $E_{\mathrm{clus}}>1$ GeV", fontsize=9)
    else:
        ax_l.set_xticklabels([])

    ratios = {"p10": ratio_hist(hists["p10"], hists["mbd"]), "p12": ratio_hist(hists["p12"], hists["mbd"])}
    for key, marker in (("p10", "o"), ("p12", "s")):
        _, _, _, color = TRIGGERS[key]
        ratio = ratios[key]
        mask = np.isfinite(ratio.values)
        ax_r.errorbar(
            ratio.centers[mask],
            ratio.values[mask],
            yerr=ratio.errors[mask],
            xerr=0.5 * ratio.widths[mask],
            fmt=marker,
            ms=2.4,
            lw=0.7,
            capsize=1.1,
            color=color,
            label=f"{'Photon 10' if key == 'p10' else 'Photon 12'} / MBD",
        )

    ax_r.axhline(1.0, color="#b8bec8", lw=0.9, zorder=0)
    ax_r.set_xlim(1.0, 20.0)
    ax_r.set_ylim(0.0, 1.22)
    ax_r.set_ylabel("trigger / MBD", fontsize=9)
    if row == 0:
        ax_r.legend(frameon=False, fontsize=8, loc="upper right", handlelength=1.5)
        ax_r.set_title("Binned turn-on ratio, no fit", fontsize=11, pad=6, fontweight="bold")
    if row == 2:
        ax_r.set_xlabel(r"Max cluster energy [GeV], $E_{\mathrm{clus}}>1$ GeV", fontsize=9)
    else:
        ax_r.set_xticklabels([])

    mbd_tail = integral(hists["mbd"], 15.0, 20.0)
    p10_tail = integral(hists["p10"], 15.0, 20.0)
    p12_tail = integral(hists["p12"], 15.0, 20.0)
    r10 = p10_tail / mbd_tail if mbd_tail > 0 else float("nan")
    r12 = p12_tail / mbd_tail if mbd_tail > 0 else float("nan")
    stat = f"15<Emax<20 GeV: P10/MBD={r10:.3f}, P12/MBD={r12:.3f}"
    ax_r.text(0.035, 0.08, stat, transform=ax_r.transAxes, fontsize=8.5, va="bottom",
              bbox=dict(boxstyle="round,pad=0.22", facecolor="white", edgecolor="#d8dde6", linewidth=0.6))

    for ax in (ax_l, ax_r):
        setup_axis(ax)


def main() -> int:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    ROOT.gROOT.SetBatch(True)
    summary = load_summary()

    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
        "mathtext.fontset": "dejavuserif",
        "axes.linewidth": 0.8,
    })

    fig = plt.figure(figsize=(16, 9), dpi=180)
    fig.patch.set_facecolor("white")

    fig.text(0.045, 0.945, "Centrality-sliced scaled-trigger check", fontsize=24, fontweight="bold", va="top")
    fig.text(
        0.045,
        0.902,
        "Same 620-run Au+Au scaled-trigger sample, now split into 0-20%, 20-50%, and 50-80% centrality bins.",
        fontsize=13,
        va="top",
    )

    box = FancyBboxPatch(
        (0.57, 0.875),
        0.39,
        0.095,
        boxstyle="round,pad=0.012,rounding_size=0.018",
        transform=fig.transFigure,
        facecolor="#edf5f0",
        edgecolor="#edf5f0",
        zorder=0,
    )
    fig.add_artist(box)
    fig.text(0.588, 0.941, "Readout strategy", fontsize=13.5, fontweight="bold", va="top")
    fig.text(
        0.588,
        0.912,
        "Ratios are shown bin-by-bin with no sigmoid fit.\n"
        "Within fixed centrality, the photon-trigger sample is not expected to plateau at unity.",
        fontsize=10.8,
        va="top",
        linespacing=1.12,
    )

    grid = fig.add_gridspec(
        nrows=3,
        ncols=2,
        left=0.055,
        right=0.965,
        bottom=0.075,
        top=0.835,
        width_ratios=[1.05, 1.0],
        hspace=0.16,
        wspace=0.14,
    )

    root_file = ROOT.TFile.Open(str(DEFAULT_INPUT), "READ")
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"could not open ROOT input: {DEFAULT_INPUT}")
    try:
        for row, (suffix, label) in enumerate(CENT_BINS):
            hists = {key: read_hist(root_file, key, suffix) for key in TRIGGERS}
            draw_row(fig, grid, row, label, suffix, hists, summary.get(suffix, {}))
    finally:
        root_file.Close()

    fig.savefig(OUTPNG, facecolor="white")
    plt.close(fig)
    print(f"Wrote {OUTPNG}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
