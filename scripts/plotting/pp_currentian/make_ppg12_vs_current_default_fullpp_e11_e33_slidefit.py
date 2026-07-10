#!/usr/bin/env python3
"""Make a slide-fit E11/E33 data overlay: PPG12 SDCC vs current default full pp.

This plot uses the validated PPG12 Fig. 13 SDCC JSON extraction and the
completed full-pp RecoilJets output that feeds the default pp shower-shape/BDT
path.  It writes an exact-pixel PNG for slide placement.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import ROOT


BASE = Path("dataOutput/ppg12PhotonYield/ppg12_photon_yield_v1_data_20260620")
DEFAULT_PPG12_JSON = (
    BASE
    / "shower_shape_reference_validation/fig13_e11_e33/"
    / "ppg12_sdcc_fig13_e11_to_e33_histograms.json"
)
DEFAULT_CURRENT_ROOT = Path(
    "InputFiles/pp24/ppg12_photon_yield_v1_data_20260620/pp/"
    "RecoilJets_pp_ALL_jetMinPtScan_dphiScan_vz60_isoR40_isSliding_"
    "preselectionNewPPG12_tightNewPPG12_nonTightNewPPG12.root"
)
DEFAULT_OUT = (
    BASE
    / "shower_shape_reference_validation/fig13_e11_e33/"
    / "ppg12_sdcc_vs_current_default_fullpp_e11_e33_data_overlay_slidefit_772x998.png"
)

HIST_DIR = "PPG12_scaledtrigger30"
CURRENT_HISTS = [
    "h_ss_e11e33_inclusive_pT_22_24",
    "h_ss_e11e33_inclusive_pT_24_26",
    "h_ss_e11e33_inclusive_pT_26_28",
]


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ppg12-json", type=Path, default=DEFAULT_PPG12_JSON)
    ap.add_argument("--current-root", type=Path, default=DEFAULT_CURRENT_ROOT)
    ap.add_argument(
        "--current-hist",
        default=None,
        help=(
            "Optional exact current histogram path inside the ROOT file. "
            "When set, this bypasses the legacy h_ss_e11e33 22-28 fallback sum."
        ),
    )
    ap.add_argument("--current-legend", default="Current default pp data")
    ap.add_argument("--current-note", default="current full pp")
    ap.add_argument("--output", type=Path, default=DEFAULT_OUT)
    return ap.parse_args()


def hist_to_arrays(hist) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    raw_entries = float(hist.Integral())
    first = hist.GetXaxis().FindBin(0.000001)
    last = hist.GetXaxis().FindBin(0.999999)
    norm = float(hist.Integral(first, last))
    if norm <= 0:
        raise RuntimeError("Current pp has no entries in 0 <= E11/E33 <= 1")

    xs = []
    ys = []
    errs = []
    axis = hist.GetXaxis()
    for ibin in range(first, last + 1):
        xs.append(0.5 * (axis.GetBinLowEdge(ibin) + axis.GetBinUpEdge(ibin)))
        ys.append(float(hist.GetBinContent(ibin)) / norm)
        errs.append(float(hist.GetBinError(ibin)) / norm)
    return np.asarray(xs), np.asarray(ys), np.asarray(errs), raw_entries


def infer_count_from_normalized_errors(values: np.ndarray, errors: np.ndarray) -> float:
    """Recover the source count when y=n/N and err=sqrt(n)/N."""
    counts = [
        (float(y) / float(err)) ** 2
        for y, err in zip(values, errors)
        if y > 0 and err > 0
    ]
    return float(sum(counts))


def extract_current(root_path: Path, current_hist: str | None) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    root_file = ROOT.TFile.Open(str(root_path))
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"Could not open current ROOT: {root_path}")

    if current_hist:
        hist = root_file.Get(current_hist)
        if not hist:
            root_file.Close()
            raise RuntimeError(f"Missing exact current histogram: {current_hist}")
        clone = hist.Clone("h_current_exact_e11e33")
        clone.SetDirectory(0)
        root_file.Close()
        return hist_to_arrays(clone)

    combined = None
    missing: list[str] = []
    for hist_name in CURRENT_HISTS:
        full_name = f"{HIST_DIR}/{hist_name}"
        hist = root_file.Get(full_name)
        if not hist:
            missing.append(full_name)
            continue
        raw_entries += float(hist.Integral())
        if combined is None:
            combined = hist.Clone("h_current_fullpp_e11e33_raw")
            combined.SetDirectory(0)
        else:
            combined.Add(hist)
    root_file.Close()

    if missing:
        raise RuntimeError(f"Missing current histograms: {missing}")
    if combined is None or combined.Integral() <= 0:
        raise RuntimeError("Current pp combined E11/E33 histogram is empty")

    rebinned = combined.Rebin(4, "h_current_fullpp_e11e33_rebinned_004")
    rebinned.SetDirectory(0)
    return hist_to_arrays(rebinned)


def main() -> None:
    args = parse_args()
    with args.ppg12_json.open() as f:
        ppg12_payload = json.load(f)["data"]

    ref_x = np.asarray(ppg12_payload["centers"], dtype=float)
    ref_y = np.asarray(ppg12_payload["values"], dtype=float)
    ref_err = np.asarray(ppg12_payload["errors"], dtype=float)
    cur_x, cur_y, cur_err, raw_entries = extract_current(args.current_root, args.current_hist)

    if len(cur_x) != len(ref_x) or not np.allclose(cur_x, ref_x, atol=1e-6):
        raise RuntimeError("Current and PPG12 E11/E33 bin centers do not match")

    ratio = np.divide(cur_y, ref_y, out=np.full_like(cur_y, np.nan), where=ref_y > 0)
    max_dev_percent = float(np.nanmax(np.abs(ratio - 1.0)) * 100.0)
    ref_entries = infer_count_from_normalized_errors(ref_y, ref_err)
    ratio_err = np.zeros_like(ratio)
    for i, (r, y, ey, yr, eyr) in enumerate(zip(ratio, cur_y, cur_err, ref_y, ref_err)):
        if not np.isfinite(r) or y <= 0 or yr <= 0:
            ratio_err[i] = np.nan
            continue
        ratio_err[i] = abs(r) * math.sqrt((ey / y) ** 2 + (eyr / yr) ** 2)

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 15,
            "axes.linewidth": 1.2,
            "xtick.major.size": 6,
            "ytick.major.size": 6,
            "xtick.minor.size": 3,
            "ytick.minor.size": 3,
        }
    )
    fig, (ax, rax) = plt.subplots(
        2,
        1,
        figsize=(7.72, 9.98),
        dpi=100,
        sharex=True,
        gridspec_kw={"height_ratios": [3.2, 1.0], "hspace": 0.04},
    )

    ax.errorbar(
        ref_x,
        ref_y,
        yerr=ref_err,
        fmt="o",
        color="black",
        ms=4.2,
        lw=1.0,
        capsize=0,
        label="PPG12 SDCC data",
        zorder=3,
    )
    ax.errorbar(
        cur_x,
        cur_y,
        yerr=cur_err,
        fmt="s",
        color="#1f77b4",
        markerfacecolor="none",
        markeredgewidth=1.3,
        ms=5.0,
        lw=1.0,
        capsize=0,
        label=args.current_legend,
        zorder=4,
    )

    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 0.18)
    ax.set_ylabel("normalized counts", fontsize=17)
    ax.minorticks_on()
    ax.tick_params(which="both", direction="in", top=True, right=True)
    ax.tick_params(axis="both", which="major", labelsize=15)

    ax.text(
        0.050,
        0.92,
        "sPHENIX",
        transform=ax.transAxes,
        fontsize=15,
        fontstyle="italic",
        fontweight="bold",
    )
    ax.text(0.235, 0.92, "Internal", transform=ax.transAxes, fontsize=15)
    ax.text(0.050, 0.845, r"$p$+$p$ $\sqrt{s}=200$ GeV", transform=ax.transAxes, fontsize=11)
    ax.text(0.050, 0.780, r"$|\eta^\gamma|<0.7$", transform=ax.transAxes, fontsize=11)
    ax.text(0.050, 0.715, r"$22<p_T<28$ GeV", transform=ax.transAxes, fontsize=11)
    ax.text(0.050, 0.650, "no NPB cut", transform=ax.transAxes, fontsize=11)
    ax.legend(
        loc="upper right",
        frameon=False,
        fontsize=16.0,
        handlelength=1.4,
        borderaxespad=0.35,
        labelspacing=0.55,
    )
    ax.text(
        0.590,
        0.790,
        (
            f"PPG12 SDCC N={ref_entries:.0f}\n"
            f"{args.current_note} N={raw_entries:.0f}\n"
            rf"max $|R-1|$ = {max_dev_percent:.1f}%"
        ),
        transform=ax.transAxes,
        fontsize=15.0,
        va="top",
        ha="left",
        linespacing=1.55,
    )

    rax.axhline(1.0, color="black", lw=1.0, ls=(0, (4, 4)))
    rax.errorbar(
        cur_x,
        ratio,
        yerr=ratio_err,
        fmt="o",
        color="#1f77b4",
        ms=4.0,
        lw=0.9,
        capsize=0,
    )
    rax.set_xlim(0.0, 1.0)
    rax.set_ylim(0.0, 2.6)
    rax.set_ylabel("Current / PPG12", fontsize=15)
    rax.set_xlabel("e11_to_e33", fontsize=17)
    rax.minorticks_on()
    rax.tick_params(which="both", direction="in", top=True, right=True)
    rax.tick_params(axis="both", which="major", labelsize=15)

    fig.subplots_adjust(left=0.145, right=0.985, top=0.985, bottom=0.115)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=100)
    plt.close(fig)
    print(args.output)


if __name__ == "__main__":
    main()
