#!/usr/bin/env python3
"""Make slide-ready comparisons for raw isolation-cone BDT ablations."""

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
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


CENT_BINS = [("0_20", "0-20%"), ("20_50", "20-50%"), ("50_80", "50-80%")]

PRODUCT_LABELS = {
    "globalEtCent1535_bdt_noIso_ptCent7": "Baseline routed BDT\n32 inputs",
    "globalEtCent1535_bdt_eisoR30_ptCent7": "Baseline + raw\nR=0.3 $E_T^{iso}$",
    "globalEtCent1535_bdt_eisoR40_ptCent7": "Baseline + raw\nR=0.4 $E_T^{iso}$",
    "globalEtCent1535_bdt_eisoR30R40_ptCent7": "Baseline + raw\nR=0.3 and R=0.4 $E_T^{iso}$",
}

COLORS = {
    "globalEtCent1535_bdt_noIso_ptCent7": "#4B5563",
    "globalEtCent1535_bdt_eisoR30_ptCent7": "#0072B2",
    "globalEtCent1535_bdt_eisoR40_ptCent7": "#F97316",
    "globalEtCent1535_bdt_eisoR30R40_ptCent7": "#CC79A7",
    "signal": "#1F77B4",
    "background": "#D62728",
}


def load_json(path: Path) -> dict:
    return json.loads(path.read_text())


def auc_by_cent(metrics: dict, product: str, cent: str) -> float:
    record = metrics["products"][product]["auc_by_centrality"][cent]
    if isinstance(record, dict):
        return float(record["auc"])
    return float(record)


def auc_inclusive(metrics: dict, product: str) -> float:
    return float(metrics["products"][product]["auc_inclusive"])


def copy_hist_product(hist: dict, product: str) -> dict:
    return hist["products"][product]


def load_baseline_hist_csv(path: Path, product: str) -> dict:
    rows_by_cent: dict[str, list[dict[str, str]]] = {}
    with path.open() as handle:
        for row in csv.DictReader(handle):
            if row["product"] == product:
                rows_by_cent.setdefault(row["cent_bin"], []).append(row)
    out = {}
    for cent, rows in rows_by_cent.items():
        rows.sort(key=lambda row: float(row["bin_lo"]))
        out[cent] = {
            "edges": [float(row["bin_lo"]) for row in rows] + [float(rows[-1]["bin_hi"])],
            "signal_density": [float(row["signal_density"]) for row in rows],
            "background_density": [float(row["background_density"]) for row in rows],
            "auc": float(rows[0]["auc"]),
            "signal_entries": int(rows[0]["signal_entries"]),
            "background_entries": int(rows[0]["background_entries"]),
        }
    if not out:
        raise RuntimeError(f"No rows for {product} in {path}")
    return out


def raw_hist_by_cent(hist: dict, product: str) -> dict:
    out = {}
    edges = hist["bin_edges"]
    for cent, _label in CENT_BINS:
        payload = hist["products"][product]["by_centrality"][cent]
        out[cent] = {
            "edges": edges,
            "signal_density": payload["signal"]["density"],
            "background_density": payload["background"]["density"],
            "auc": payload.get("auc"),
            "signal_entries": payload["signal"]["entries"],
            "background_entries": payload["background"]["entries"],
        }
    return out


def step_density(ax, edges, density, *, color, label):
    edges = np.asarray(edges, dtype=float)
    density = np.asarray(density, dtype=float)
    y = np.r_[density, density[-1]]
    ax.step(edges, y, where="post", color=color, lw=2.1, label=label)
    ax.fill_between(edges, y, step="post", color=color, alpha=0.075)


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "axes.linewidth": 1.05,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "axes.grid": True,
            "grid.color": "0.90",
            "grid.linewidth": 0.55,
        }
    )


def make_auc_summary(
    baseline_metrics: dict,
    raw_metrics: dict,
    r40_metrics: dict,
    outdir: Path,
) -> None:
    rows = []
    products = [
        "globalEtCent1535_bdt_noIso_ptCent7",
        "globalEtCent1535_bdt_eisoR30_ptCent7",
        "globalEtCent1535_bdt_eisoR40_ptCent7",
        "globalEtCent1535_bdt_eisoR30R40_ptCent7",
    ]
    for product in products:
        if product == "globalEtCent1535_bdt_noIso_ptCent7":
            metrics = baseline_metrics
        elif product == "globalEtCent1535_bdt_eisoR40_ptCent7":
            metrics = r40_metrics
        else:
            metrics = raw_metrics
        row = {
            "product": product,
            "label": PRODUCT_LABELS[product].replace("\n", " "),
            "auc_inclusive": auc_inclusive(metrics, product),
        }
        for cent, _ in CENT_BINS:
            row[f"auc_{cent}"] = auc_by_cent(metrics, product, cent)
        rows.append(row)

    csv_path = outdir / "eiso_raw_auc_summary.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    setup_style()
    fig, axes = plt.subplots(1, 2, figsize=(15.4, 6.0), gridspec_kw={"width_ratios": [1.0, 1.35]}, dpi=170)
    fig.patch.set_facecolor("white")

    labels = [PRODUCT_LABELS[row["product"]] for row in rows]
    inclusive = [row["auc_inclusive"] for row in rows]
    colors = [COLORS[row["product"]] for row in rows]
    y = np.arange(len(rows))
    ax = axes[0]
    ax.barh(y, inclusive, color=colors, height=0.58)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=12.2)
    ax.invert_yaxis()
    ax.set_xlim(0.82, 0.945)
    ax.set_xlabel("Inclusive AUC", fontsize=13)
    ax.set_title("Full 15-35 GeV sample", fontsize=14, fontweight="bold")
    for yy, val in zip(y, inclusive):
        ax.text(val + 0.002, yy, f"{val:.3f}", va="center", ha="left", fontsize=12.2, fontweight="bold")

    ax = axes[1]
    x = np.arange(len(CENT_BINS))
    width = 0.19
    for i, row in enumerate(rows):
        vals = [row[f"auc_{cent}"] for cent, _ in CENT_BINS]
        offset = (i - 1.5) * width
        ax.bar(x + offset, vals, width=width, color=COLORS[row["product"]], label=PRODUCT_LABELS[row["product"]].replace("\n", " "))
        for xx, val in zip(x + offset, vals):
            ax.text(xx, val + 0.003, f"{val:.3f}", ha="center", va="bottom", fontsize=10.2, rotation=90)
    ax.set_xticks(x)
    ax.set_xticklabels([label for _, label in CENT_BINS], fontsize=12)
    ax.set_ylim(0.80, 0.955)
    ax.set_ylabel("AUC", fontsize=13)
    ax.set_title("Centrality bins", fontsize=14, fontweight="bold")
    ax.legend(frameon=False, fontsize=10.4, loc="upper left", bbox_to_anchor=(0.02, -0.16), ncol=1)

    fig.text(0.045, 0.970, r"$\bf{\it{sPHENIX}}$ Internal", ha="left", va="top", fontsize=16)
    fig.text(0.045, 0.930, "Photon12+20 signal vs Jet12+20+30 background", ha="left", va="top", fontsize=12)
    fig.text(0.52, 0.970, r"Raw $E_T^{iso}$ cone-input ablation: routed BDT AUC", ha="center", va="top", fontsize=19, fontweight="bold")
    fig.tight_layout(rect=[0.035, 0.155, 0.985, 0.895])
    fig.savefig(outdir / "eiso_raw_auc_summary.png")
    plt.close(fig)


def make_score_overlay_table(
    baseline_hist: dict,
    raw_hist: dict,
    baseline_metrics: dict,
    raw_metrics: dict,
    outdir: Path,
    raw_label: str,
    raw_product: str,
) -> None:
    setup_style()
    fig, axes = plt.subplots(2, 3, figsize=(15.8, 8.1), sharex=True, dpi=170)
    row_defs = [
        ("Baseline routed BDT", "globalEtCent1535_bdt_noIso_ptCent7", baseline_hist, baseline_metrics),
        (raw_label, raw_product, raw_hist, raw_metrics),
    ]
    summary = []
    for irow, (row_label, product, hist, metrics) in enumerate(row_defs):
        for icol, (cent, cent_label) in enumerate(CENT_BINS):
            ax = axes[irow, icol]
            payload = hist[cent]
            step_density(ax, payload["edges"], payload["signal_density"], color=COLORS["signal"], label="Signal")
            step_density(ax, payload["edges"], payload["background_density"], color=COLORS["background"], label="Background")
            auc = auc_by_cent(metrics, product, cent)
            ax.text(
                0.045,
                0.875,
                f"AUC = {auc:.3f}",
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=12.5,
                fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.24", fc="white", ec="0.82", alpha=0.93),
            )
            ymax = max(np.nanmax(payload["signal_density"]), np.nanmax(payload["background_density"]))
            ax.set_ylim(0, float(ymax) * 1.25)
            ax.set_xlim(0, 1)
            if irow == 0:
                ax.set_title(f"Centrality {cent_label}", fontsize=14, fontweight="bold")
            if icol == 0:
                ax.set_ylabel("Area-normalized density", fontsize=12.5)
                ax.text(
                    -0.155,
                    0.5,
                    row_label,
                    transform=ax.transAxes,
                    rotation=90,
                    ha="center",
                    va="center",
                    fontsize=13.2,
                    fontweight="bold",
                )
            if irow == 1:
                ax.set_xlabel("Classifier score", fontsize=12.5)
            if irow == 0 and icol == 2:
                ax.legend(frameon=False, fontsize=11.5, loc="upper right")
            summary.append(
                {
                    "row_label": row_label,
                    "product": product,
                    "centrality": cent_label,
                    "auc": auc,
                    "signal_entries": payload["signal_entries"],
                    "background_entries": payload["background_entries"],
                }
            )

    fig.text(0.045, 0.976, r"$\bf{\it{sPHENIX}}$ Internal", ha="left", va="top", fontsize=15.5)
    fig.text(0.045, 0.942, r"$15 < E_T < 35$ GeV", ha="left", va="top", fontsize=11.5)
    fig.text(0.52, 0.976, "Score separation: baseline vs raw isolation inputs", ha="center", va="top", fontsize=18.5, fontweight="bold")
    fig.text(
        0.52,
        0.942,
        r"same routed 8 $E_T$ x 7 centrality BDT family; bottom row adds raw cone isolation inputs",
        ha="center",
        va="top",
        fontsize=12,
        color="#374151",
    )
    fig.tight_layout(rect=[0.055, 0.055, 0.985, 0.895])
    fig.savefig(outdir / f"baseline_vs_{raw_product.replace('globalEtCent1535_bdt_', '')}_score_separation_2x3.png")
    plt.close(fig)

    csv_path = outdir / f"baseline_vs_{raw_product.replace('globalEtCent1535_bdt_', '')}_score_separation_2x3.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary[0]))
        writer.writeheader()
        writer.writerows(summary)


def make_slide35_gain_plot(
    baseline_metrics: dict,
    old_iso_metrics: dict,
    raw_metrics: dict,
    r40_metrics: dict,
    outdir: Path,
    *,
    include_r30_only: bool = True,
    output_stem: str = "eiso_raw_slide35_style_auc_gain_0_20",
) -> None:
    reference_product = "globalEtCent1535_bdt_noIso_ptCent7"
    reference_auc = auc_by_cent(baseline_metrics, reference_product, "0_20")
    rows = [
        {
            "label": "Reference\n8 $p_T$ x 7 centrality",
            "badge": "56 BDTs\n32 inputs",
            "auc": reference_auc,
            "product": reference_product,
            "color": "#64748B",
            "band": "#FFFFFF",
            "kind": "reference",
        },
        {
            "label": "+ raw R=0.3\n$E_T^{iso}$ input",
            "badge": "56 BDTs\n33 inputs",
            "auc": auc_by_cent(raw_metrics, "globalEtCent1535_bdt_eisoR30_ptCent7", "0_20"),
            "product": "globalEtCent1535_bdt_eisoR30_ptCent7",
            "color": "#0072B2",
            "band": "#EFF6FF",
            "kind": "variant",
        },
        {
            "label": "+ raw R=0.4\n$E_T^{iso}$ input",
            "badge": "56 BDTs\n33 inputs",
            "auc": auc_by_cent(r40_metrics, "globalEtCent1535_bdt_eisoR40_ptCent7", "0_20"),
            "product": "globalEtCent1535_bdt_eisoR40_ptCent7",
            "color": "#F97316",
            "band": "#FFF7ED",
            "kind": "variant",
        },
        {
            "label": "+ raw R=0.3 and R=0.4\n$E_T^{iso}$ inputs",
            "badge": "56 BDTs\n34 inputs",
            "auc": auc_by_cent(raw_metrics, "globalEtCent1535_bdt_eisoR30R40_ptCent7", "0_20"),
            "product": "globalEtCent1535_bdt_eisoR30R40_ptCent7",
            "color": "#CC79A7",
            "band": "#FDF2F8",
            "kind": "variant",
        },
    ]
    if not include_r30_only:
        rows = [row for row in rows if row["product"] != "globalEtCent1535_bdt_eisoR30_ptCent7"]

    csv_path = outdir / f"{output_stem}.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["product", "label", "input_count", "auc_0_20", "relative_gain_percent", "absolute_gain"],
        )
        writer.writeheader()
        for row in rows:
            input_count = row["badge"].split("\n")[-1].split()[0]
            writer.writerow(
                {
                    "product": row["product"],
                    "label": row["label"].replace("\n", " "),
                    "input_count": input_count,
                    "auc_0_20": f"{row['auc']:.9f}",
                    "relative_gain_percent": f"{100.0 * (row['auc'] / reference_auc - 1.0):.6f}",
                    "absolute_gain": f"{row['auc'] - reference_auc:.9f}",
                }
            )

    setup_style()
    plt.rcParams.update({"font.family": "DejaVu Sans"})
    fig, ax = plt.subplots(figsize=(16.0, 9.0), dpi=170)
    fig.patch.set_facecolor("white")

    n = len(rows)
    y_positions = np.arange(n)
    xmin = reference_auc - 0.0065
    xmax = max(row["auc"] for row in rows) + 0.0135
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(-0.8, n - 0.35)
    ax.invert_yaxis()
    ax.set_yticks(y_positions)
    ax.set_yticklabels([row["label"] for row in rows], fontsize=13.2, ha="right")
    ax.tick_params(axis="y", length=0, pad=18)
    ax.tick_params(axis="x", labelsize=13, top=True, direction="in")
    ax.grid(axis="x", color="#D1D5DB", linewidth=0.8, alpha=0.9)
    ax.grid(axis="y", visible=False)
    for spine in ("left", "right", "top"):
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_linewidth(1.1)
    ax.set_xlabel("0-20% centrality validation AUC", fontsize=16)

    for y, row in zip(y_positions, rows):
        ax.axhspan(y - 0.43, y + 0.43, color=row["band"], zorder=0)
    ax.axvline(reference_auc, color="#475569", lw=1.7, ls=(0, (4, 4)), zorder=1)

    for y, row in zip(y_positions, rows):
        color = row["color"]
        if row["kind"] == "reference":
            ax.scatter([row["auc"]], [y], s=145, color=color, edgecolor="#334155", linewidth=1.5, zorder=5)
            xtext = row["auc"] + 0.00055
            ax.text(xtext, y - 0.12, f"AUC {row['auc']:.4f}", ha="left", va="center", fontsize=15.5, fontweight="bold", color="#111827")
            ax.text(xtext, y + 0.16, "reference", ha="left", va="center", fontsize=11.8, color="#475569")
        else:
            ax.plot([reference_auc, row["auc"]], [y, y], color=color, lw=7.5, alpha=0.34, solid_capstyle="round", zorder=2)
            ax.scatter([row["auc"]], [y], s=155, color=color, edgecolor="#0F766E" if color == "#009E73" else color, linewidth=1.8, zorder=5)
            xtext = row["auc"] + 0.00055
            gain = 100.0 * (row["auc"] / reference_auc - 1.0)
            ax.text(xtext, y - 0.12, f"AUC {row['auc']:.4f}", ha="left", va="center", fontsize=15.5, fontweight="bold", color="#111827")
            ax.text(xtext, y + 0.17, f"+{gain:.2f}% gain in AUC", ha="left", va="center", fontsize=12.3, color=color)
        ax.text(
            reference_auc - 0.00025,
            y,
            row["badge"],
            ha="right",
            va="center",
            fontsize=10.6,
            color="#64748B",
            bbox=dict(boxstyle="round,pad=0.20", facecolor="white", edgecolor="#CBD5E1", alpha=0.92),
        )

    fig.text(0.265, 0.934, r"AUC gain from raw $E_T^{iso}$ inputs (0-20%)", ha="left", va="top", fontsize=21.5, fontweight="bold")
    fig.text(0.265, 0.884, r"Photon12+20 signal & Jet12+20+30 inclusive background; $15 < p_T < 35$ GeV", ha="left", va="top", fontsize=13.5, color="#4B5563")
    fig.text(0.265, 0.852, r"All rows use 8 $p_T$ x 7 centrality routing; gains are relative to the current 56-BDT reference", ha="left", va="top", fontsize=12.3, color="#4B5563")
    fig.text(0.780, 0.910, "sPHENIX", ha="left", va="top", fontsize=17, fontstyle="italic", fontweight="bold")
    fig.text(0.861, 0.910, " Internal", ha="left", va="top", fontsize=17)
    fig.text(
        0.265,
        0.060,
        "Raw cone-isolation inputs produce the largest 0-20% AUC gain; these are diagnostic inputs, not the ABCD baseline model.",
        ha="left",
        va="bottom",
        fontsize=10.8,
        color="#6B7280",
    )
    fig.subplots_adjust(left=0.305, right=0.955, top=0.805, bottom=0.145)
    fig.savefig(outdir / f"{output_stem}.png")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--baseline-validation", type=Path, required=True)
    parser.add_argument("--raw-validation", type=Path, required=True)
    parser.add_argument("--r40-validation", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    baseline_metrics = load_json(args.baseline_validation / "validation_metrics.json")
    raw_metrics = load_json(args.raw_validation / "validation_metrics.json")
    r40_metrics = load_json(args.r40_validation / "validation_metrics.json")
    raw_hist_json = load_json(args.raw_validation / "score_histograms.json")
    r40_hist_json = load_json(args.r40_validation / "score_histograms.json")

    baseline_hist = load_baseline_hist_csv(
        args.baseline_validation / "coarse_centrality_score_histograms_binned_noIso.csv",
        "globalEtCent1535_bdt_noIso_ptCent7",
    )
    raw_candidates = [
        ("+ raw R=0.3 $E_T^{iso}$", "globalEtCent1535_bdt_eisoR30_ptCent7", raw_metrics, raw_hist_json),
        ("+ raw R=0.4 $E_T^{iso}$", "globalEtCent1535_bdt_eisoR40_ptCent7", r40_metrics, r40_hist_json),
        ("+ raw R=0.3 and R=0.4 $E_T^{iso}$", "globalEtCent1535_bdt_eisoR30R40_ptCent7", raw_metrics, raw_hist_json),
    ]
    best_label, best_product, best_metrics, best_hist_json = max(
        raw_candidates,
        key=lambda item: auc_by_cent(item[2], item[1], "0_20"),
    )
    raw_hist = raw_hist_by_cent(best_hist_json, best_product)

    make_auc_summary(baseline_metrics, raw_metrics, r40_metrics, args.outdir)
    make_score_overlay_table(
        baseline_hist,
        raw_hist,
        baseline_metrics,
        best_metrics,
        args.outdir,
        raw_label=best_label,
        raw_product=best_product,
    )
    make_slide35_gain_plot(baseline_metrics, baseline_metrics, raw_metrics, r40_metrics, args.outdir)
    make_slide35_gain_plot(
        baseline_metrics,
        baseline_metrics,
        raw_metrics,
        r40_metrics,
        args.outdir,
        include_r30_only=False,
        output_stem="eiso_raw_slide35_style_auc_gain_0_20_no_r30_only",
    )
    print(args.outdir / "eiso_raw_auc_summary.png")
    print(args.outdir / f"baseline_vs_{best_product.replace('globalEtCent1535_bdt_', '')}_score_separation_2x3.png")
    print(args.outdir / "eiso_raw_slide35_style_auc_gain_0_20.png")
    print(args.outdir / "eiso_raw_slide35_style_auc_gain_0_20_no_r30_only.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
