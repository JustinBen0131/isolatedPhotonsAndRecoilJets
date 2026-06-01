#!/usr/bin/env python3
"""Plot reco-ABCD and truth-labeled purity from RecoilJets SIM ROOT pairs."""

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

import argparse
import csv
import math
import os
import re
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-codex")

import matplotlib.pyplot as plt
import ROOT


ROOT.gROOT.SetBatch(True)

COARSE_CENTS: dict[str, tuple[str, ...]] = {
    "0_20": ("0_10", "10_20"),
    "20_50": ("20_30", "30_40", "40_50"),
    "50_80": ("50_60", "60_80"),
}

CENT_LABELS = {
    "0_20": "0-20%",
    "20_50": "20-50%",
    "50_80": "50-80%",
}

REGION_PREFIX = {
    "A": "h_isIsolated_isTight",
    "B": "h_notIsolated_isTight",
    "C": "h_isIsolated_notTight",
    "D": "h_notIsolated_notTight",
}

TRUTH_SIGNAL_PREFIX = "h_EisoReco_truthSigMatched_tight"
TRUTH_BACKGROUND_PREFIX = "h_Eiso_tight"


@dataclass(frozen=True)
class Sample:
    label: str
    signal_root: Path
    background_root: Path


def parse_sample(text: str) -> Sample:
    pieces = text.split("|")
    if len(pieces) != 3:
        raise argparse.ArgumentTypeError(
            "--sample must be 'label|signal_root|background_root'"
        )
    return Sample(pieces[0], Path(pieces[1]), Path(pieces[2]))


def open_sim(path: Path):
    handle = ROOT.TFile.Open(str(path), "READ")
    if not handle or handle.IsZombie():
        raise RuntimeError(f"Could not open {path}")
    directory = handle.Get("SIM")
    if not directory:
        handle.Close()
        raise RuntimeError(f"Missing SIM directory in {path}")
    return handle, directory


def hist_sum_and_error(hist) -> tuple[float, float]:
    total = 0.0
    err2 = 0.0
    for ibin in range(0, hist.GetNbinsX() + 2):
        total += float(hist.GetBinContent(ibin))
        err2 += float(hist.GetBinError(ibin)) ** 2
    return total, math.sqrt(err2)


def index_histograms(directory, prefixes: tuple[str, ...], cone: str) -> dict[tuple[str, str, float, float], str]:
    patterns = {
        prefix: re.compile(rf"^{re.escape(prefix)}_{re.escape(cone)}_pT_([0-9]+)_([0-9]+)_cent_(.+)$")
        for prefix in prefixes
    }
    index: dict[tuple[str, str, float, float], str] = {}
    keys = directory.GetListOfKeys()
    for ikey in range(keys.GetEntries()):
        name = keys.At(ikey).GetName()
        for prefix, pattern in patterns.items():
            match = pattern.match(name)
            if not match:
                continue
            lo = float(match.group(1))
            hi = float(match.group(2))
            cent = match.group(3)
            index[(prefix, cent, lo, hi)] = name
    return index


def bins_for_cent(
    index: dict[tuple[str, str, float, float], str],
    prefix: str,
    cent: str,
    pt_min: float,
    pt_max: float,
) -> set[tuple[float, float]]:
    direct = {
        (lo, hi)
        for key_prefix, key_cent, lo, hi in index
        if key_prefix == prefix and key_cent == cent and pt_min <= 0.5 * (lo + hi) <= pt_max
    }
    if direct:
        return direct

    bins: set[tuple[float, float]] = set()
    for fine_cent in COARSE_CENTS[cent]:
        bins.update(
            (lo, hi)
            for key_prefix, key_cent, lo, hi in index
            if key_prefix == prefix
            and key_cent == fine_cent
            and pt_min <= 0.5 * (lo + hi) <= pt_max
        )
    return bins


def sum_count_for_bin(
    directory,
    index: dict[tuple[str, str, float, float], str],
    prefix: str,
    cent: str,
    lo: float,
    hi: float,
) -> tuple[float, float]:
    pieces = (cent,) if (prefix, cent, lo, hi) in index else COARSE_CENTS[cent]

    total = 0.0
    err2 = 0.0
    found = False
    for piece in pieces:
        name = index.get((prefix, piece, lo, hi))
        if not name:
            continue
        hist = directory.Get(name)
        if not hist:
            continue
        value, error = hist_sum_and_error(hist)
        total += value
        err2 += error * error
        found = True
    if not found:
        raise KeyError(f"Missing {prefix}, pT {lo:g}-{hi:g}, cent {cent}")
    return total, math.sqrt(err2)


def raw_abcd_purity(
    counts: dict[str, tuple[float, float]],
) -> tuple[float, float, float]:
    a, ea = counts["A"]
    b, eb = counts["B"]
    c, ec = counts["C"]
    d, ed = counts["D"]
    if a <= 0.0 or d <= 0.0:
        return math.nan, math.nan, math.nan

    bkg_est = b * c / d
    raw = 1.0 - bkg_est / a
    value = max(0.0, raw)

    dp_da = (b * c) / (a * a * d)
    dp_db = -c / (a * d)
    dp_dc = -b / (a * d)
    dp_dd = (b * c) / (a * d * d)
    var = (
        dp_da * dp_da * ea * ea
        + dp_db * dp_db * eb * eb
        + dp_dc * dp_dc * ec * ec
        + dp_dd * dp_dd * ed * ed
    )
    return value, math.sqrt(max(var, 0.0)), bkg_est


def ratio_purity(signal: tuple[float, float], background: tuple[float, float]) -> tuple[float, float]:
    sig, esig = signal
    bkg, ebkg = background
    denom = sig + bkg
    if denom <= 0.0:
        return math.nan, math.nan
    value = sig / denom
    err = math.sqrt((bkg * bkg * esig * esig + sig * sig * ebkg * ebkg) / (denom**4))
    return value, err


def collect_rows(
    sample: Sample,
    abcd_cone: str,
    truth_cone: str,
    pt_min: float,
    pt_max: float,
) -> list[dict]:
    sig_handle, sig_dir = open_sim(sample.signal_root)
    bkg_handle, bkg_dir = open_sim(sample.background_root)
    try:
        region_prefixes = tuple(REGION_PREFIX.values())
        truth_prefixes = (TRUTH_SIGNAL_PREFIX, TRUTH_BACKGROUND_PREFIX)
        sig_region_index = index_histograms(sig_dir, region_prefixes, abcd_cone)
        bkg_region_index = index_histograms(bkg_dir, region_prefixes, abcd_cone)
        sig_truth_index = index_histograms(sig_dir, truth_prefixes, truth_cone)
        bkg_truth_index = index_histograms(bkg_dir, truth_prefixes, truth_cone)

        rows: list[dict] = []
        for cent in COARSE_CENTS:
            candidate_bins = set()
            for prefix in REGION_PREFIX.values():
                candidate_bins |= bins_for_cent(sig_region_index, prefix, cent, pt_min, pt_max)
                candidate_bins |= bins_for_cent(bkg_region_index, prefix, cent, pt_min, pt_max)
            truth_bins = bins_for_cent(sig_truth_index, TRUTH_SIGNAL_PREFIX, cent, pt_min, pt_max)
            truth_bins &= bins_for_cent(bkg_truth_index, TRUTH_BACKGROUND_PREFIX, cent, pt_min, pt_max)

            for lo, hi in sorted(candidate_bins):
                combined_counts: dict[str, tuple[float, float]] = {}
                raw_counts: dict[str, float] = {}
                raw_errors: dict[str, float] = {}
                for region, prefix in REGION_PREFIX.items():
                    sig_count = sum_count_for_bin(sig_dir, sig_region_index, prefix, cent, lo, hi)
                    bkg_count = sum_count_for_bin(bkg_dir, bkg_region_index, prefix, cent, lo, hi)
                    total = sig_count[0] + bkg_count[0]
                    error = math.hypot(sig_count[1], bkg_count[1])
                    combined_counts[region] = (total, error)
                    raw_counts[region] = total
                    raw_errors[region] = error
                value, error, bkg_est = raw_abcd_purity(combined_counts)
                rows.append(
                    {
                        "model": sample.label,
                        "quantity": "reco_abcd_raw_purity",
                        "centrality": cent,
                        "pt_low": lo,
                        "pt_high": hi,
                        "pt_mid": 0.5 * (lo + hi),
                        "value": value,
                        "error": error,
                        "A": raw_counts["A"],
                        "A_error": raw_errors["A"],
                        "B": raw_counts["B"],
                        "B_error": raw_errors["B"],
                        "C": raw_counts["C"],
                        "C_error": raw_errors["C"],
                        "D": raw_counts["D"],
                        "D_error": raw_errors["D"],
                        "abcd_background_estimate": bkg_est,
                        "signal_tight": "",
                        "signal_tight_error": "",
                        "background_tight": "",
                        "background_tight_error": "",
                        "signal_root": str(sample.signal_root),
                        "background_root": str(sample.background_root),
                        "definition": "max(0, A - B*C/D)/A using combined signal-plus-inclusive reco ABCD counts",
                    }
                )

            for lo, hi in sorted(truth_bins):
                sig_tight = sum_count_for_bin(sig_dir, sig_truth_index, TRUTH_SIGNAL_PREFIX, cent, lo, hi)
                bkg_tight = sum_count_for_bin(bkg_dir, bkg_truth_index, TRUTH_BACKGROUND_PREFIX, cent, lo, hi)
                value, error = ratio_purity(sig_tight, bkg_tight)
                rows.append(
                    {
                        "model": sample.label,
                        "quantity": "truth_labeled_A_purity",
                        "centrality": cent,
                        "pt_low": lo,
                        "pt_high": hi,
                        "pt_mid": 0.5 * (lo + hi),
                        "value": value,
                        "error": error,
                        "A": "",
                        "A_error": "",
                        "B": "",
                        "B_error": "",
                        "C": "",
                        "C_error": "",
                        "D": "",
                        "D_error": "",
                        "abcd_background_estimate": "",
                        "signal_tight": sig_tight[0],
                        "signal_tight_error": sig_tight[1],
                        "background_tight": bkg_tight[0],
                        "background_tight_error": bkg_tight[1],
                        "signal_root": str(sample.signal_root),
                        "background_root": str(sample.background_root),
                        "definition": "S_A_truthMatched/(S_A_truthMatched+B_A_inclusive) using tight isolated region A",
                    }
                )
        return rows
    finally:
        sig_handle.Close()
        bkg_handle.Close()


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise RuntimeError("No rows to write")
    fieldnames = list(rows[0].keys())
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def read_csv(path: Path) -> list[dict]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def labels_in_order(rows: list[dict]) -> list[str]:
    labels: list[str] = []
    seen: set[str] = set()
    for row in rows:
        label = row["model"]
        if label in seen:
            continue
        labels.append(label)
        seen.add(label)
    return labels


def sphinx_internal(ax) -> None:
    ax.text(
        0.02,
        0.96,
        r"$\bf{\it{sPHENIX}}$ Internal",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=11,
    )


def plot_rows(
    path: Path,
    rows: list[dict],
    labels: list[str],
    quantity_filter: str = "all",
    style: str = "default",
) -> None:
    def color_for(label: str) -> str:
        if label == labels[0] or "Reference" in label:
            return "#222222"
        if "Global" in label or "global" in label:
            return "#CC79A7"
        return "#0072B2"

    def marker_for(label: str) -> str:
        if label == labels[0] or "Reference" in label:
            return "o"
        if "Global" in label or "global" in label:
            return "D"
        return "s"
    quantities = [
        ("reco_abcd_raw_purity", "Reco ABCD estimate", r"$\max(0,A-BC/D)/A$"),
        ("truth_labeled_A_purity", "Truth-labeled region A purity", r"$S_A/(S_A+B_A)$"),
    ]
    if quantity_filter != "all":
        quantities = [entry for entry in quantities if entry[0] == quantity_filter]
        if not quantities:
            raise RuntimeError(f"Unknown quantity filter: {quantity_filter}")

    if style not in {"default", "slide10"}:
        raise RuntimeError(f"Unknown plot style: {style}")

    nrows = len(quantities)
    if style == "slide10" and nrows != 1:
        raise RuntimeError("--style slide10 is currently defined for one-row plots")

    if style == "slide10":
        plt.rcParams.update(
            {
                "font.family": "DejaVu Sans",
                "axes.linewidth": 1.0,
                "xtick.major.size": 6,
                "ytick.major.size": 6,
                "xtick.minor.size": 3,
                "ytick.minor.size": 3,
                "xtick.direction": "in",
                "ytick.direction": "in",
            }
        )

    fig_height = 4.45 if nrows == 1 else 7.9
    # Match the actual reco-efficiency plot raster embedded on JSTG slide 10.
    fig_size = (3586.0 / 190.0, 1192.0 / 190.0) if style == "slide10" else (15.2, fig_height)
    if style == "slide10":
        fig = plt.figure(figsize=fig_size)
        # Pixel-matched axes boxes from the slide-10 source PNG
        # (basev3e_target80_vs_reference_boxcuts_reco_efficiency_pt15to35_1x3.png).
        width_px = 3586.0
        height_px = 1192.0
        bottom = (height_px - 1012.0) / height_px
        height = (1012.0 - 334.0) / height_px
        boxes = [
            (206.0 / width_px, bottom, (1145.0 - 206.0) / width_px, height),
            (1286.0 / width_px, bottom, (2340.0 - 1286.0) / width_px, height),
            (2481.0 / width_px, bottom, (3518.0 - 2481.0) / width_px, height),
        ]
        axes = [[fig.add_axes(box) for box in boxes]]
    else:
        fig, axes_obj = plt.subplots(nrows, 3, figsize=fig_size, sharex=True, sharey=False)
        axes = [axes_obj] if nrows == 1 else axes_obj
    for row_idx, (quantity, ylabel, subtitle) in enumerate(quantities):
        for col_idx, cent in enumerate(COARSE_CENTS):
            ax = axes[row_idx][col_idx]
            ax.set_xlim(14.5, 35.2)
            ax.set_ylim(0.0, 0.9 if style == "slide10" else 1.05)
            if style == "slide10":
                ax.set_xticks([15, 20, 25, 30, 35])
                ax.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8])
                ax.minorticks_on()
                ax.tick_params(axis="both", which="major", direction="in", labelsize=12, top=True, right=True, length=6, width=1.0)
                ax.tick_params(axis="both", which="minor", direction="in", top=True, right=True, length=3, width=0.8, color="#b3b3b3")
            else:
                ax.grid(True, color="#d8d8d8", linewidth=0.8, alpha=0.65)
                ax.tick_params(axis="both", which="major", labelsize=10)
            ax.set_title(
                CENT_LABELS[cent],
                fontsize=20 if style == "slide10" else 13,
                fontweight="bold",
                pad=10 if style == "slide10" else 8,
            )
            if col_idx == 0:
                ax.set_ylabel(
                    "Raw reco ABCD purity" if quantity == "reco_abcd_raw_purity" and style == "slide10" else ylabel,
                    fontsize=18 if style == "slide10" else 12,
                )
            if row_idx == nrows - 1 and (style != "slide10" or col_idx == 1):
                ax.set_xlabel(r"reco photon $p_{T}$ [GeV]", fontsize=16 if style == "slide10" else 12)
            if row_idx == 0 and col_idx == 0 and style != "slide10":
                sphinx_internal(ax)
            if col_idx == 2 and style != "slide10":
                ax.text(
                    0.97,
                    0.94,
                    subtitle,
                    transform=ax.transAxes,
                    ha="right",
                    va="top",
                    fontsize=10,
                    bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78, "pad": 3},
                )

            for label in labels:
                points = [
                    r
                    for r in rows
                    if r["model"] == label
                    and r["quantity"] == quantity
                    and r["centrality"] == cent
                    and math.isfinite(float(r["value"]))
                ]
                if not points:
                    continue
                x = [float(r["pt_mid"]) for r in points]
                y = [float(r["value"]) for r in points]
                xerr = [0.5 * (float(r["pt_high"]) - float(r["pt_low"])) for r in points]
                yerr = [float(r["error"]) for r in points]
                if style == "slide10":
                    marker_style = "o"
                    marker_size = 6.5
                    edge_color = color_for(label)
                    edge_width = 0.0
                    xerr_values = None
                else:
                    marker_style = marker_for(label)
                    marker_size = 7.5
                    edge_color = "white"
                    edge_width = 0.9
                    xerr_values = xerr
                ax.errorbar(
                    x,
                    y,
                    yerr=yerr,
                    xerr=xerr_values,
                    linestyle="none",
                    marker=marker_style,
                    markersize=marker_size,
                    markerfacecolor=color_for(label),
                    markeredgecolor=edge_color,
                    markeredgewidth=edge_width,
                    ecolor=color_for(label),
                    elinewidth=1.35 if style == "slide10" else 1.45,
                    capsize=2.4 if style == "slide10" else 3.0,
                    label=label,
                    alpha=0.96,
                )

    if style == "slide10":
        fig.text(
            0.070,
            0.995,
            r"$\bf{\it{sPHENIX}}$ Internal",
            ha="left",
            va="top",
            fontsize=18,
        )
        fig.text(
            0.070,
            0.940,
            "Embedded Photon12+20 signal & Jet12+20+30 inclusive background",
            ha="left",
            va="top",
            fontsize=14,
            color="#333333",
        )
        fig.text(
            0.070,
            0.900,
            r"Reco ABCD raw purity = $\max(0,A-BC/D)/A$; sliding isolation, R = 0.3",
            ha="left",
            va="top",
            fontsize=14,
            color="#333333",
        )
        handles, legend_labels = axes[0][2].get_legend_handles_labels()
        if handles:
            legend = axes[0][2].legend(
                handles,
                legend_labels,
                loc="lower right",
                frameon=True,
                framealpha=0.92,
                fontsize=13,
                borderpad=0.55,
                handlelength=1.8,
            )
            legend.get_frame().set_linewidth(0)
    else:
        handles, legend_labels = axes[0][2].get_legend_handles_labels()
        if handles:
            fig.legend(
                handles,
                legend_labels,
                loc="upper center",
                bbox_to_anchor=(0.52, 0.988 if nrows == 1 else 0.986),
                ncol=max(1, len(legend_labels)),
                frameon=False,
                fontsize=12,
            )
        fig.text(
            0.52,
            0.018,
            "Embedded Photon12+20 signal and Jet12+20+30 inclusive background; sliding isolation, R = 0.3",
            ha="center",
            fontsize=11,
        )
        top_rect = 0.88 if nrows == 1 else 0.935
        fig.tight_layout(rect=(0.018, 0.060, 0.995, top_rect), w_pad=1.0, h_pad=1.05)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=190 if style == "slide10" else 220)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", action="append", type=parse_sample, default=[])
    parser.add_argument("--input-csv", action="append", type=Path, default=[])
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--stem", default="reference_vs_ptcent7bdt_reco_truth_raw_abcd_purity_pt15to35")
    parser.add_argument("--pt-min", type=float, default=15.0)
    parser.add_argument("--pt-max", type=float, default=35.0)
    parser.add_argument("--abcd-cone", default="isoR30_isSliding")
    parser.add_argument("--truth-cone", default="isoR30")
    parser.add_argument(
        "--plot-quantity",
        choices=("all", "reco_abcd_raw_purity", "truth_labeled_A_purity"),
        default="all",
    )
    parser.add_argument("--style", choices=("default", "slide10"), default="default")
    args = parser.parse_args()

    all_rows: list[dict] = []
    for input_csv in args.input_csv:
        all_rows.extend(read_csv(input_csv))
    for sample in args.sample:
        all_rows.extend(collect_rows(sample, args.abcd_cone, args.truth_cone, args.pt_min, args.pt_max))
    if not all_rows:
        raise RuntimeError("No rows produced; pass at least one --sample or --input-csv")
    labels = labels_in_order(all_rows)

    csv_path = args.outdir / f"{args.stem}_points.csv"
    layout_suffix = "2x3" if args.plot_quantity == "all" else "1x3"
    png_path = args.outdir / f"{args.stem}_{layout_suffix}.png"
    write_csv(csv_path, all_rows)
    plot_rows(png_path, all_rows, labels, args.plot_quantity, args.style)
    print(csv_path)
    print(png_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
