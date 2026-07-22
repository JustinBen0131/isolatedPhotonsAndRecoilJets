#!/usr/bin/env python3
"""Overlay PPG12 Fig.20 SDCC BDT data with a registered pp output."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FixedLocator, FuncFormatter, MultipleLocator


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
CAMPAIGN = "the76_ppg12_parity_full_20260701_003024"
SDCC_JSON = (
    REPO
    / "dataOutput/ppg12PhotonYield/ppg12_photon_yield_v1_data_20260620"
    / "shower_shape_reference_validation/fig20_bdt"
    / "ppg12_sdcc_fig20_bdt_tight_fourcurve_projections.json"
)
CURRENT_ROOT = (
    REPO
    / f"dataOutput/ppg12Parity/{CAMPAIGN}/final_roots/pp_data_hierarchical_v1"
    / "RecoilJets_pp_ALL_PERIOD_COMBINED.root"
)
OUTDIR = REPO / f"dataOutput/ppg12Parity/{CAMPAIGN}/data_bdt_fig20"
CURRENT_DIR = "PPG12_scaledtrigger30"
CURRENT_HIST = f"{CURRENT_DIR}/h2d_bdt_eta0_pt0_cut2"
SDCC_HIST = "h2d_bdt_eta0_pt0_cut2"
CURRENT_LABEL = "This analysis"


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "dejavuserif",
            "axes.linewidth": 1.15,
            "xtick.major.size": 6,
            "ytick.major.size": 6,
            "xtick.minor.size": 3,
            "ytick.minor.size": 3,
            "xtick.direction": "in",
            "ytick.direction": "in",
        }
    )


def rootish_tick(value: float, _pos: int) -> str:
    if abs(value) < 1e-12:
        return "0"
    return f"{value:.2f}".rstrip("0").rstrip(".")


def load_sdcc() -> dict[str, np.ndarray | float | str]:
    payload = json.loads(SDCC_JSON.read_text())
    panel = payload["panel"]
    curve = panel["curves"]["data"]
    return {
        "centers": np.asarray(curve["centers"], dtype=float),
        "edges": np.asarray(curve["edges"], dtype=float),
        "values": np.asarray(curve["values"], dtype=float),
        "errors": np.asarray(curve["errors"], dtype=float),
        "raw_integral": float(curve["raw_integral"]),
        "source_selection": payload["source_selection"],
        "source_root": panel["summary"]["source_files"]["data"],
        "chi2": float(panel["chi2"]),
        "ndf": int(panel["ndf"]),
    }


def rebin_counts_to_edges(src_values: np.ndarray, src_edges: np.ndarray, dst_edges: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    out = np.zeros(len(dst_edges) - 1, dtype=float)
    err2 = np.zeros_like(out)
    for i, count in enumerate(src_values):
        lo = float(src_edges[i])
        hi = float(src_edges[i + 1])
        if hi <= dst_edges[0] or lo >= dst_edges[-1]:
            continue
        width = hi - lo
        if width <= 0:
            continue
        for j in range(len(out)):
            ov = max(0.0, min(hi, float(dst_edges[j + 1])) - max(lo, float(dst_edges[j])))
            if ov <= 0:
                continue
            frac = ov / width
            out[j] += count * frac
            err2[j] += count * frac * frac
    return out, np.sqrt(err2)


def load_current(dst_edges: np.ndarray) -> dict[str, np.ndarray | float]:
    import ROOT

    root_file = ROOT.TFile.Open(str(CURRENT_ROOT))
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"Could not open current ROOT: {CURRENT_ROOT}")
    hist = root_file.Get(CURRENT_HIST)
    if not hist:
        root_file.Close()
        raise RuntimeError(f"Missing current histogram: {CURRENT_HIST}")

    source_class = hist.ClassName()
    if hist.InheritsFrom("TH2"):
        hist_for_x = hist.ProjectionX(f"{hist.GetName()}_px_for_fig20_overlay", 1, hist.GetNbinsY(), "e")
    else:
        hist_for_x = hist

    nbins = hist_for_x.GetNbinsX()
    counts = np.asarray([float(hist_for_x.GetBinContent(i)) for i in range(1, nbins + 1)], dtype=float)
    edges = np.asarray(
        [float(hist_for_x.GetXaxis().GetBinLowEdge(i)) for i in range(1, nbins + 1)]
        + [float(hist_for_x.GetXaxis().GetBinUpEdge(nbins))],
        dtype=float,
    )
    underflow = float(hist_for_x.GetBinContent(0))
    overflow = float(hist_for_x.GetBinContent(nbins + 1))
    entries = float(hist.GetEntries())
    root_file.Close()

    raw, raw_err = rebin_counts_to_edges(counts, edges, dst_edges)
    total = float(np.sum(raw))
    values = raw / total if total > 0 else raw
    errors = raw_err / total if total > 0 else raw_err
    centers = 0.5 * (dst_edges[:-1] + dst_edges[1:])
    return {
        "centers": centers,
        "values": values,
        "errors": errors,
        "raw": raw,
        "raw_errors": raw_err,
        "total_0to1": total,
        "entries": entries,
        "underflow": underflow,
        "overflow": overflow,
        "source_class": source_class,
        "source_axis_low": float(edges[0]),
        "source_axis_high": float(edges[-1]),
    }


def ratio_summary(ref_x: np.ndarray, ref_y: np.ndarray, cmp_x: np.ndarray, cmp_y: np.ndarray, stable_threshold: float) -> dict[str, float | int | None]:
    ref_interp = np.interp(cmp_x, ref_x, ref_y, left=np.nan, right=np.nan)
    ratio = np.divide(cmp_y, ref_interp, out=np.full_like(cmp_y, np.nan), where=ref_interp > 0)
    finite = np.isfinite(ratio) & (ref_interp > 0)
    stable = finite & (ref_interp > stable_threshold)
    stable_ratio = ratio[stable]
    if stable_ratio.size:
        idx = int(np.nanargmax(np.abs(stable_ratio - 1.0)))
        stable_x = cmp_x[stable]
        return {
            "ratio_plot_bin_count": int(finite.sum()),
            "stable_ratio_count": int(stable.sum()),
            "stable_threshold_on_reference": float(stable_threshold),
            "stable_mean_ratio": float(np.nanmean(stable_ratio)),
            "stable_mean_abs_ratio_minus_one": float(np.nanmean(np.abs(stable_ratio - 1.0))),
            "stable_max_abs_ratio_minus_one": float(np.nanmax(np.abs(stable_ratio - 1.0))),
            "stable_max_abs_ratio_minus_one_bdt_center": float(stable_x[idx]),
        }
    return {
        "ratio_plot_bin_count": int(finite.sum()),
        "stable_ratio_count": int(stable.sum()),
        "stable_threshold_on_reference": float(stable_threshold),
        "stable_mean_ratio": None,
        "stable_mean_abs_ratio_minus_one": None,
        "stable_max_abs_ratio_minus_one": None,
        "stable_max_abs_ratio_minus_one_bdt_center": None,
    }


def draw_overlay(sdcc: dict[str, np.ndarray | float | str], current: dict[str, np.ndarray | float], out: Path, ratio_ylim: tuple[float, float]) -> dict[str, object]:
    setup_style()
    ref_x = sdcc["centers"]
    ref_y = sdcc["values"]
    ref_err = sdcc["errors"]
    cur_x = current["centers"]
    cur_y = current["values"]
    cur_err = current["errors"]
    if len(cur_x) != len(ref_x) or not np.allclose(cur_x, ref_x, atol=1e-6):
        raise RuntimeError("Current and PPG12 Fig.20 BDT bin centers do not match")

    stable_threshold = 0.00285
    summary = ratio_summary(ref_x, ref_y, cur_x, cur_y, stable_threshold)
    max_dev = summary["stable_max_abs_ratio_minus_one"]
    max_x = summary["stable_max_abs_ratio_minus_one_bdt_center"]
    stats_text = (
        f"PPG12 SDCC N={sdcc['raw_integral']:.0f}\n"
        f"This analysis N={current['total_0to1']:.0f}\n"
        + (rf"max $|R-1|$ = {100.0 * float(max_dev):.1f}%" if max_dev is not None else r"max $|R-1|$ = n/a")
        + (f"\nnear bdt = {float(max_x):.2f}" if max_x is not None else "")
    )

    fig, (ax, rax) = plt.subplots(
        2,
        1,
        figsize=(7.72, 9.98),
        dpi=100,
        sharex=True,
        gridspec_kw={"height_ratios": [3.35, 1.0], "hspace": 0.05},
    )
    ax.errorbar(ref_x, ref_y, yerr=ref_err, fmt="o", color="black", ms=4.2, lw=1.0, label="PPG12 SDCC ROOT data", zorder=3)
    ax.errorbar(
        cur_x,
        cur_y,
        yerr=cur_err,
        fmt="s",
        color="#1f77b4",
        markerfacecolor="none",
        markeredgewidth=1.2,
        ms=4.6,
        lw=1.0,
        label=CURRENT_LABEL,
        zorder=4,
    )

    ax.text(0.05, 0.93, "sPHENIX", transform=ax.transAxes, ha="left", va="top", fontsize=17, fontweight="bold", fontstyle="italic")
    ax.text(0.265, 0.93, "Internal", transform=ax.transAxes, ha="left", va="top", fontsize=17)
    ax.text(0.05, 0.84, r"$p{+}p\ \sqrt{s}=200$ GeV", transform=ax.transAxes, fontsize=15, ha="left")
    ax.text(0.05, 0.77, r"$|\eta^\gamma| < 0.7$", transform=ax.transAxes, fontsize=15, ha="left")
    ax.text(0.05, 0.70, r"$10 < p_T < 14$ GeV, w/ tight cut", transform=ax.transAxes, fontsize=13, ha="left")
    ax.set_ylabel("normalized counts", fontsize=17)
    ax.set_xlim(0.0, 1.0)
    ymax = 1.12 * max(float(np.nanmax(ref_y + ref_err)), float(np.nanmax(cur_y + cur_err)))
    ax.set_ylim(0.0, ymax)
    ax.yaxis.set_major_locator(MultipleLocator(0.05))
    ax.yaxis.set_minor_locator(MultipleLocator(0.01))
    ax.yaxis.set_major_formatter(FuncFormatter(rootish_tick))
    ax.legend(loc="upper right", frameon=False, fontsize=16.0, handlelength=1.5, borderaxespad=0.35, labelspacing=0.55)
    ax.text(0.595, 0.790, stats_text, transform=ax.transAxes, fontsize=15.0, va="top", ha="left", linespacing=1.55)
    ax.tick_params(labelsize=14, top=True, right=True)
    ax.minorticks_on()

    ref_interp = np.interp(cur_x, ref_x, ref_y, left=np.nan, right=np.nan)
    ratio = np.divide(cur_y, ref_interp, out=np.full_like(cur_y, np.nan), where=ref_interp > 0)
    ratio_err = np.divide(cur_err, ref_interp, out=np.full_like(cur_err, np.nan), where=ref_interp > 0)
    plot_mask = np.isfinite(ratio) & (ref_interp > 0)
    rax.errorbar(cur_x[plot_mask], ratio[plot_mask], yerr=ratio_err[plot_mask], fmt="o", color="#1f77b4", ms=4.0, lw=1.0)
    rax.axhline(1.0, color="0.35", ls="--", lw=1.0)
    rax.set_ylabel("This analysis / PPG12", fontsize=15)
    rax.set_xlabel("bdt", fontsize=17)
    rax.set_ylim(*ratio_ylim)
    rax.xaxis.set_major_locator(FixedLocator(np.arange(0.0, 1.0001, 0.2)))
    rax.xaxis.set_minor_locator(MultipleLocator(0.05))
    rax.xaxis.set_major_formatter(FuncFormatter(rootish_tick))
    rax.tick_params(labelsize=14, top=True, right=True)
    rax.minorticks_on()

    fig.subplots_adjust(left=0.115, right=0.985, top=0.985, bottom=0.075, hspace=0.05)
    fig.savefig(out, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)

    return {
        "artifact": str(out),
        "comparison": "PPG12 Fig.20 SDCC ROOT data projection vs current registered full-stat pp TableQA BDT tight TH2 projection",
        "campaign": CAMPAIGN,
        "ppg12_sdcc_json": str(SDCC_JSON),
        "ppg12_sdcc_root": sdcc["source_root"],
        "ppg12_sdcc_hist": SDCC_HIST,
        "ppg12_sdcc_source_selection": sdcc["source_selection"],
        "ppg12_sdcc_raw_integral_0to1": sdcc["raw_integral"],
        "ppg12_fig20_chi2_ndf_fingerprint": [sdcc["chi2"], sdcc["ndf"], sdcc["chi2"] / sdcc["ndf"]],
        "current_root": str(CURRENT_ROOT),
        "current_hist": CURRENT_HIST,
        "current_hist_class": current["source_class"],
        "current_rebinned_integral_0to1": current["total_0to1"],
        "current_entries": current["entries"],
        "current_underflow": current["underflow"],
        "current_overflow": current["overflow"],
        "current_source_axis": [current["source_axis_low"], current["source_axis_high"]],
        "ratio_ylim": list(ratio_ylim),
        **summary,
    }


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ratio-ymin", type=float, default=0.0)
    ap.add_argument("--ratio-ymax", type=float, default=4.5)
    ap.add_argument("--current-root", type=Path, default=CURRENT_ROOT)
    ap.add_argument("--current-label", default=CURRENT_LABEL)
    ap.add_argument("--campaign-tag", default=CAMPAIGN)
    ap.add_argument(
        "--output",
        type=Path,
        default=OUTDIR / "final_hierarchical_v1_bdt_tight_10_14_tableqa_h2d_pt0_cut2_overlay_slidecrop_allbins_ratio0to4p5.png",
    )
    return ap.parse_args()


def main() -> None:
    global CURRENT_ROOT, CURRENT_LABEL, CAMPAIGN
    args = parse_args()
    CURRENT_ROOT = args.current_root.resolve()
    CURRENT_LABEL = args.current_label
    CAMPAIGN = args.campaign_tag
    OUTDIR.mkdir(parents=True, exist_ok=True)
    sdcc = load_sdcc()
    current = load_current(sdcc["edges"])
    manifest = draw_overlay(sdcc, current, args.output, (args.ratio_ymin, args.ratio_ymax))
    args.output.with_suffix(".manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True))
    print(args.output)


if __name__ == "__main__":
    main()
