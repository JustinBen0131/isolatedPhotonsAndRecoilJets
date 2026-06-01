#!/usr/bin/env python3
"""Extract and render the pp reco-cluster ET leakage slide.

This is the pp counterpart of the Au+Au reco-cluster ET leakage check.  It
uses Justin/RecoilJets pp inclusive output trees, not PPG12/Shuhang ROOT files
or truth-stitch CSVs.
"""

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
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np


REPO = Path(__file__).resolve().parents[1]
DEFAULT_LOCAL_BASE = (
    REPO
    / "dataOutput/ppPhotonMLPipeline/ppg12_basev3E_currentIAN_recoEtExt50_20260527_0045"
)
DEFAULT_OUTDIR = DEFAULT_LOCAL_BASE / "validation/reco_cluster_et_leakage"
DEFAULT_SLIDE = DEFAULT_LOCAL_BASE / "slide_assets/pp_reco_cluster_et_leakage_slide7_style.png"

DEFAULT_SOURCE_ROOT = Path(
    "/sphenix/tg/tg01/bulk/jbennett/thesisAna/ppPhotonMLPipeline/"
    "ppg12_basev3E_currentIAN_recoEtExt50_20260527_0045"
)
TREE_NAME = "AuAuPhotonIDTrainingTree"
BRANCH = "cluster_Et"
META_PATH = "SIM/h_ppInclusiveJetStitch_ppg12TruthSpectrum_metadata"


@dataclass(frozen=True)
class Sample:
    key: str
    label: str
    xsec_pb: float
    cap_gev: float
    window: str
    color: str
    marker: str


SAMPLES: tuple[Sample, ...] = (
    Sample("run28_jet8", "jet8", 1.3013e7, 15.0, r"$p_T^{truth\,jet}<14$", "#d22c98", "o"),
    Sample("run28_jet12", "jet12", 1.4903e6, 23.0, r"$14\leq p_T^{truth\,jet}<21$", "#2ca02c", "s"),
    Sample("run28_jet20", "jet20", 6.2623e4, 35.0, r"$21\leq p_T^{truth\,jet}<32$", "#0090ff", "^"),
    Sample("run28_jet30", "jet30", 2.5298e3, 45.0, r"$32\leq p_T^{truth\,jet}<42$", "#ff6b00", "v"),
    Sample("run28_jet40", "jet40", 1.3553e2, 100.0, r"$p_T^{truth\,jet}\geq42$", "#cc00cc", "D"),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    extract = sub.add_parser("extract", help="Read pp inclusive ROOT trees and write compact CSV/JSON.")
    extract.add_argument("--source-root", type=Path, default=DEFAULT_SOURCE_ROOT)
    extract.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    extract.add_argument("--bins", default="5:50:1", help="lo:hi:width in GeV")
    extract.add_argument("--max-files-per-sample", type=int, default=0)
    extract.add_argument("--progress-every", type=int, default=100)

    render = sub.add_parser("render", help="Render the slide-ready PNG from compact CSV/JSON.")
    render.add_argument("--csv", type=Path, default=DEFAULT_OUTDIR / "pp_reco_cluster_et_leakage_components.csv")
    render.add_argument("--summary", type=Path, default=DEFAULT_OUTDIR / "pp_reco_cluster_et_leakage_summary.json")
    render.add_argument("--png", type=Path, default=DEFAULT_SLIDE)

    both = sub.add_parser("both", help="Extract and then render.")
    both.add_argument("--source-root", type=Path, default=DEFAULT_SOURCE_ROOT)
    both.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    both.add_argument("--png", type=Path, default=DEFAULT_SLIDE)
    both.add_argument("--bins", default="5:50:1")
    both.add_argument("--max-files-per-sample", type=int, default=0)
    both.add_argument("--progress-every", type=int, default=100)
    return parser.parse_args()


def parse_bins(spec: str) -> np.ndarray:
    lo_s, hi_s, width_s = spec.split(":")
    lo = float(lo_s)
    hi = float(hi_s)
    width = float(width_s)
    if not (hi > lo and width > 0):
        raise ValueError(f"invalid bin spec: {spec}")
    n = int(round((hi - lo) / width))
    edges = lo + width * np.arange(n + 1, dtype=float)
    if not np.isclose(edges[-1], hi):
        raise ValueError(f"bin spec does not end exactly at high edge: {spec}")
    return edges


def sample_files(source_root: Path, sample: Sample, max_files: int) -> list[Path]:
    sample_dir = (
        source_root
        / "isSimInclusive"
        / "preselectionNewPPG12_tightReference_nonTightReference"
        / sample.key
    )
    paths = sorted(sample_dir.glob("*.root"))
    if max_files > 0:
        paths = paths[:max_files]
    if not paths:
        raise FileNotFoundError(f"no ROOT files found for {sample.key} under {sample_dir}")
    return paths


def read_metadata(uproot_file: Any) -> dict[str, float]:
    if META_PATH not in uproot_file:
        raise KeyError(f"missing {META_PATH}")
    hist = uproot_file[META_PATH]
    values = np.asarray(hist.values(), dtype=float)
    labels = list(hist.axis().labels())
    out: dict[str, float] = {}
    for idx, label in enumerate(labels):
        if idx < len(values) and label:
            out[str(label)] = float(values[idx])
    if "events_processed" not in out:
        raise KeyError(f"{META_PATH} has no events_processed bin label")
    return out


def extract(args: argparse.Namespace) -> tuple[Path, Path]:
    import uproot

    edges = parse_bins(args.bins)
    bin_widths = np.diff(edges)
    if not np.allclose(bin_widths, bin_widths[0]):
        raise ValueError("only uniform bins are supported")
    bin_width = float(bin_widths[0])

    rows: list[dict[str, Any]] = []
    summary_samples: list[dict[str, Any]] = []
    for sample in SAMPLES:
        paths = sample_files(args.source_root, sample, args.max_files_per_sample)
        counts = np.zeros(len(edges) - 1, dtype=np.float64)
        entries_seen = 0
        events_processed = 0.0
        metadata_xsecs: list[float] = []
        missing_tree = 0
        missing_branch = 0
        nonfinite = 0
        above_cap = 0
        et_min = math.inf
        et_max = -math.inf

        for idx, path in enumerate(paths, start=1):
            with uproot.open(path) as root_file:
                meta = read_metadata(root_file)
                events_processed += float(meta["events_processed"])
                if "xsec_pb" in meta and math.isfinite(meta["xsec_pb"]):
                    metadata_xsecs.append(float(meta["xsec_pb"]))
                if TREE_NAME not in root_file:
                    missing_tree += 1
                    continue
                tree = root_file[TREE_NAME]
                if BRANCH not in tree.keys():
                    missing_branch += 1
                    continue
                arr = tree[BRANCH].array(library="np")
                vals = np.asarray(arr, dtype=np.float64)
                finite = np.isfinite(vals)
                nonfinite += int((~finite).sum())
                vals = vals[finite]
                entries_seen += int(vals.size)
                if vals.size:
                    et_min = min(et_min, float(np.min(vals)))
                    et_max = max(et_max, float(np.max(vals)))
                    above_cap += int(np.sum(vals > sample.cap_gev))
                    counts += np.histogram(vals, bins=edges)[0].astype(np.float64)
            if args.progress_every and idx % args.progress_every == 0:
                print(f"[extract] {sample.key}: {idx}/{len(paths)} files", flush=True)

        if events_processed <= 0 or not math.isfinite(events_processed):
            raise RuntimeError(f"bad events_processed denominator for {sample.key}: {events_processed}")
        weight = sample.xsec_pb / events_processed / bin_width
        raw_err = np.sqrt(counts)
        weighted = counts * weight
        weighted_err = raw_err * weight
        for ibin, (lo, hi, raw, err, w, we) in enumerate(
            zip(edges[:-1], edges[1:], counts, raw_err, weighted, weighted_err, strict=True),
            start=1,
        ):
            rows.append(
                {
                    "sample": sample.label,
                    "sample_key": sample.key,
                    "bin_index": ibin,
                    "x_low": float(lo),
                    "x_high": float(hi),
                    "x_center": float(0.5 * (lo + hi)),
                    "x_err_low": float(0.5 * (hi - lo)),
                    "x_err_high": float(0.5 * (hi - lo)),
                    "bin_width": bin_width,
                    "raw_entries": float(raw),
                    "raw_error": float(err),
                    "weighted_entries": float(w),
                    "weighted_error": float(we),
                    "xsec_pb_used": sample.xsec_pb,
                    "events_processed": events_processed,
                    "per_entry_weight": weight,
                    "cluster_et_cap_gev_recorded": sample.cap_gev,
                    "truth_jet_window": sample.window,
                    "source_root": str(args.source_root),
                }
            )
        summary_samples.append(
            {
                "sample": sample.label,
                "sample_key": sample.key,
                "file_count": len(paths),
                "tree_entries_seen": entries_seen,
                "events_processed": events_processed,
                "xsec_pb_used": sample.xsec_pb,
                "metadata_xsec_pb_values": sorted(set(round(x, 8) for x in metadata_xsecs)),
                "per_entry_weight": weight,
                "cluster_et_cap_gev_recorded": sample.cap_gev,
                "entries_above_recorded_cap": above_cap,
                "missing_tree_files": missing_tree,
                "missing_branch_files": missing_branch,
                "nonfinite_cluster_et": nonfinite,
                "cluster_et_min": None if et_min == math.inf else et_min,
                "cluster_et_max": None if et_max == -math.inf else et_max,
            }
        )

    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    csv_path = outdir / "pp_reco_cluster_et_leakage_components.csv"
    summary_path = outdir / "pp_reco_cluster_et_leakage_summary.json"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    summary = build_summary(rows, summary_samples, args.source_root, csv_path, summary_path)
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(f"[extract] wrote {csv_path}", flush=True)
    print(f"[extract] wrote {summary_path}", flush=True)
    return csv_path, summary_path


def build_summary(
    rows: list[dict[str, Any]],
    summary_samples: list[dict[str, Any]],
    source_root: Path,
    csv_path: Path,
    summary_path: Path,
) -> dict[str, Any]:
    by_bin: dict[tuple[float, float], dict[str, float]] = {}
    for row in rows:
        key = (float(row["x_low"]), float(row["x_high"]))
        by_bin.setdefault(key, {})[str(row["sample"])] = float(row["weighted_entries"])
    fraction_sum_max_dev = 0.0
    finite_fraction_bins = 0
    ranges = [(5, 9), (9, 15), (15, 23), (23, 35), (35, 45), (45, 50)]
    ranges_out = []
    for lo, hi in ranges:
        weighted_entries = {s.label: 0.0 for s in SAMPLES}
        for row in rows:
            if float(row["x_low"]) >= lo and float(row["x_high"]) <= hi:
                weighted_entries[str(row["sample"])] += float(row["weighted_entries"])
        total = sum(weighted_entries.values())
        fractions = {
            sample: (value / total if total > 0 else 0.0)
            for sample, value in weighted_entries.items()
        }
        ranges_out.append(
            {
                "x_range": f"{lo:g}-{hi:g}",
                "weighted_sum": total,
                "weighted_entries": weighted_entries,
                "fractions": fractions,
            }
        )
    for sample_values in by_bin.values():
        total = sum(sample_values.values())
        if total > 0:
            frac_sum = sum(value / total for value in sample_values.values())
            fraction_sum_max_dev = max(fraction_sum_max_dev, abs(frac_sum - 1.0))
            finite_fraction_bins += 1

    return {
        "schema": "PP_RECO_CLUSTER_ET_LEAKAGE_V1",
        "source_root": str(source_root),
        "tree": TREE_NAME,
        "branch": BRANCH,
        "metadata_histogram": META_PATH,
        "weight_formula": "weighted_entries = raw_candidate_count * xsec_pb / events_processed / bin_width",
        "xsec_source": "sPHENIX wiki-updated pp inclusive cross sections recorded in local RecoilJets constants",
        "output_csv": str(csv_path),
        "output_summary": str(summary_path),
        "samples": summary_samples,
        "binning": {
            "low_gev": float(rows[0]["x_low"]),
            "high_gev": float(rows[-1]["x_high"]),
            "width_gev": float(rows[0]["bin_width"]),
            "nbins": int(max(int(row["bin_index"]) for row in rows)),
        },
        "ranges_of_interest": ranges_out,
        "qa": {
            "sample_count": len(summary_samples),
            "all_samples_present": len(summary_samples) == len(SAMPLES),
            "all_denominators_positive": all(s["events_processed"] > 0 for s in summary_samples),
            "all_weights_finite": all(math.isfinite(s["per_entry_weight"]) for s in summary_samples),
            "fraction_sum_max_deviation": fraction_sum_max_dev,
            "finite_fraction_bins": finite_fraction_bins,
        },
        "plot_note": (
            "This plot answers which pp inclusive source contributes to the reco-cluster "
            "ET population after applying accepted per-event cross-section stitching weights."
        ),
    }


def load_rows(csv_path: Path) -> tuple[list[dict[str, Any]], Any]:
    import pandas as pd

    df = pd.read_csv(csv_path)
    return [], df


def render(csv_path: Path, summary_path: Path, png_path: Path) -> Path:
    import matplotlib.pyplot as plt
    import pandas as pd
    from matplotlib.patches import FancyBboxPatch

    df = pd.read_csv(csv_path)
    summary = json.loads(summary_path.read_text())

    plt.rcParams.update(
        {
            "font.family": "Times New Roman",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "stix",
            "axes.linewidth": 1.2,
            "axes.labelsize": 16,
            "xtick.labelsize": 13,
            "ytick.labelsize": 13,
            "legend.fontsize": 12,
            "figure.dpi": 192,
            "savefig.dpi": 192,
        }
    )

    labels = [s.label for s in SAMPLES]
    colors = {s.label: s.color for s in SAMPLES}
    markers = {s.label: s.marker for s in SAMPLES}
    window = {s.label: s.window for s in SAMPLES}
    legend_window = {
        "jet8": r"$p_T^{truth\,jet}<14$",
        "jet12": r"$14-21$",
        "jet20": r"$21-32$",
        "jet30": r"$32-42$",
        "jet40": r"$\geq42$",
    }
    pivot = (
        df.pivot_table(
            index=["x_low", "x_high", "x_center", "x_err_low", "x_err_high"],
            columns="sample",
            values=["weighted_entries", "weighted_error"],
            aggfunc="sum",
        )
        .sort_index()
        .reset_index()
    )
    pivot.columns = [
        "_".join(str(x) for x in col if str(x)) if isinstance(col, tuple) else str(col)
        for col in pivot.columns
    ]
    for label in labels:
        for base in ("weighted_entries", "weighted_error"):
            col = f"{base}_{label}"
            if col not in pivot:
                pivot[col] = 0.0
    pivot["weighted_entries_sum"] = sum(pivot[f"weighted_entries_{label}"] for label in labels)
    pivot["weighted_error_sum"] = np.sqrt(sum(pivot[f"weighted_error_{label}"] ** 2 for label in labels))
    denom = pivot["weighted_entries_sum"].replace(0, np.nan)
    for label in labels:
        pivot[f"fraction_{label}"] = pivot[f"weighted_entries_{label}"] / denom

    def add_round_box(fig: plt.Figure, xy: tuple[float, float], wh: tuple[float, float], color: str) -> None:
        fig.add_artist(
            FancyBboxPatch(
                xy,
                wh[0],
                wh[1],
                boxstyle="round,pad=0.012,rounding_size=0.018",
                linewidth=0,
                facecolor=color,
                transform=fig.transFigure,
                zorder=0,
            )
        )

    def range_fraction(range_label: str, sample_label: str) -> float:
        for item in summary.get("ranges_of_interest", []):
            if item.get("x_range") == range_label:
                return float(item.get("fractions", {}).get(sample_label, 0.0))
        return float("nan")

    def pct(value: float) -> str:
        return "n/a" if not math.isfinite(value) else f"{100.0 * value:.0f}%"

    x = pivot["x_center"].to_numpy(float)
    xerr = np.vstack([pivot["x_err_low"].to_numpy(float), pivot["x_err_high"].to_numpy(float)])

    fig = plt.figure(figsize=(2560 / 192, 1440 / 192), constrained_layout=False)
    fig.patch.set_facecolor("white")
    fig.text(
        0.045,
        0.955,
        r"pp Reco-cluster $E_T$ leakage check",
        fontsize=34,
        fontweight="bold",
        ha="left",
        va="top",
    )

    add_round_box(fig, (0.045, 0.742), (0.435, 0.134), "#f3f5f8")
    add_round_box(fig, (0.505, 0.742), (0.45, 0.134), "#eef8f1")
    fig.text(0.063, 0.862, "What is plotted", fontsize=17.4, fontweight="bold", va="top")
    fig.text(
        0.063,
        0.829,
        (
            r"Reco photon-cluster $E_T$ in pp inclusive-jet MC."
            "\n1 GeV bins; black markers are the weighted sum."
            "\nBottom panel shows each sample fraction."
            "\nShown over the populated 12-40 GeV range."
        ),
        fontsize=12.9,
        va="top",
        linespacing=1.10,
    )
    fig.text(0.523, 0.862, "Main takeaway", fontsize=17.4, fontweight="bold", va="top")
    fig.text(
        0.523,
        0.829,
        (
            "This is a pp weighting/merge diagnostic."
            f"\njet8 dominates 9-15 GeV ({pct(range_fraction('9-15', 'jet8'))}); "
            f"jet12 carries 15-23 GeV ({pct(range_fraction('15-23', 'jet12'))})."
            "\nThe high-ET tail is sparse and less smooth than Au+Au."
            "\nNext: compare with the Au+Au-style estimator."
        ),
        fontsize=12.7,
        va="top",
        linespacing=1.07,
    )

    gs = fig.add_gridspec(
        nrows=2,
        ncols=1,
        height_ratios=[3.1, 1.25],
        left=0.105,
        right=0.965,
        bottom=0.130,
        top=0.697,
        hspace=0.055,
    )
    ax = fig.add_subplot(gs[0])
    ax_frac = fig.add_subplot(gs[1], sharex=ax)

    for label in labels:
        y = pivot[f"weighted_entries_{label}"].to_numpy(float)
        yerr = pivot[f"weighted_error_{label}"].to_numpy(float)
        mask = y > 0
        ax.errorbar(
            x[mask],
            y[mask],
            yerr=yerr[mask],
            xerr=xerr[:, mask],
            fmt=markers[label],
            color=colors[label],
            ecolor=colors[label],
            elinewidth=1.8,
            capsize=3.2,
            markersize=8.0,
            linestyle="none",
            label=f"{label}: {legend_window[label]}",
            zorder=4,
        )
    sum_y = pivot["weighted_entries_sum"].to_numpy(float)
    sum_err = pivot["weighted_error_sum"].to_numpy(float)
    sum_mask = sum_y > 0
    ax.errorbar(
        x[sum_mask],
        sum_y[sum_mask],
        yerr=sum_err[sum_mask],
        xerr=xerr[:, sum_mask],
        fmt="D",
        color="#111111",
        ecolor="#111111",
        markerfacecolor="white",
        markeredgewidth=1.6,
        elinewidth=1.9,
        capsize=3.2,
        markersize=8.4,
        linestyle="none",
        label="weighted sum",
        zorder=6,
    )

    for label in labels:
        frac = pivot[f"fraction_{label}"].to_numpy(float)
        ax_frac.plot(
            x,
            frac,
            marker=markers[label],
            markersize=8.2,
            linestyle="none",
            color=colors[label],
            label=label,
        )

    ax.set_yscale("log")
    ax.set_xlim(12, 40)
    positive = sum_y[sum_y > 0]
    y_min = max(float(np.min(positive)) * 0.25, 5.0e-2) if positive.size else 5.0e-2
    y_max = float(np.max(positive)) * 2.7 if positive.size else 1.0e6
    ax.set_ylim(y_min, y_max)
    ax.set_ylabel("weighted entries / 1 GeV bin", labelpad=12)
    ax.grid(which="major", color="#d7dce3", linewidth=0.95, alpha=0.78)
    ax.grid(which="minor", color="#edf0f4", linewidth=0.55, alpha=0.55)
    ax.tick_params(which="both", direction="in", top=True, right=True, length=7)
    ax.tick_params(which="minor", length=3.5)
    ax.tick_params(labelbottom=False)
    ax.legend(
        loc="upper right",
        bbox_to_anchor=(0.72, 0.965),
        ncol=2,
        frameon=True,
        facecolor="white",
        edgecolor="white",
        framealpha=0.9,
        fontsize=11.7,
        handlelength=1.35,
        columnspacing=1.0,
        handletextpad=0.45,
        borderpad=0.45,
        labelspacing=0.30,
    )
    ax.text(
        0.985,
        0.965,
        r"$\it{\bf{sPHENIX}}$ Internal",
        transform=ax.transAxes,
        fontsize=16,
        ha="right",
        va="top",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.84, "pad": 2.0},
        zorder=30,
    )
    ax.text(
        0.985,
        0.895,
        r"PYTHIA8 pp inclusive jet, $\sqrt{s}=200$ GeV",
        transform=ax.transAxes,
        fontsize=13.0,
        ha="right",
        va="top",
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.84, "pad": 1.8},
        zorder=30,
    )

    ax_frac.set_ylim(-0.035, 1.05)
    ax_frac.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])
    ax_frac.set_ylabel("fraction of\nweighted sum", labelpad=12)
    ax_frac.set_xlabel(r"reco photon-cluster $E_T$ ($p_T^\gamma$) [GeV]", labelpad=3)
    ax_frac.grid(which="major", color="#d7dce3", linewidth=0.95, alpha=0.78)
    ax_frac.tick_params(which="both", direction="in", top=True, right=True, length=7)

    fig.text(
        0.085,
        0.049,
        r"1 GeV reco-cluster $E_T$ bins, 12-40 GeV display; follow-up: audit pp weighting/merge against the Au+Au-style estimator",
        fontsize=15.5,
        color="#5e6675",
        ha="left",
        va="top",
    )

    png_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png_path, facecolor="white")
    plt.close(fig)
    print(f"[render] wrote {png_path}", flush=True)
    return png_path


def main() -> None:
    args = parse_args()
    if args.command == "extract":
        extract(args)
    elif args.command == "render":
        render(args.csv, args.summary, args.png)
    elif args.command == "both":
        csv_path, summary_path = extract(args)
        render(csv_path, summary_path, args.png)
    else:
        raise SystemExit(f"unknown command {args.command}")


if __name__ == "__main__":
    main()
