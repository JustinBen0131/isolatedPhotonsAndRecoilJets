#!/usr/bin/env python3
"""Make a slide-10-style Au+Au BDT score panel with pp reference overlays."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


REPO = Path(__file__).resolve().parents[1]
AUAU_REPORT = REPO / "dataOutput/auauTightBDTValidation/model_validation_condor_20260509_192942"
PP_CSV = REPO / "dataOutput/ppPhotonMLPipeline/ppg12_fix4_score_overlays/score_separation/pp_baseline_bdt_score_histograms_compact.csv"
OUT_DIR = PP_CSV.parent / "slide10_pp_overlay"
OUT_PNG = OUT_DIR / "slide10_style_auau_with_pp_reference_score_separation_1x3.png"
OUT_UNIT_PNG = OUT_DIR / "slide10_style_auau_pp_unit_normalized_shape_overlay_1x3.png"
OUT_TWO_ROW_PNG = OUT_DIR / "slide10_style_pp_top_auau_bottom_score_separation_2x3.png"
OUT_META = OUT_DIR / "slide10_style_auau_with_pp_reference_score_separation_1x3_metadata.json"

SCORE_KEY = "score_ptCentDepFine"
CENT_BINS = [(0.0, 20.0, "0-20%"), (20.0, 50.0, "20-50%"), (50.0, 80.0, "50-80%")]
SCORE_BINS = np.linspace(0.0, 1.0, 51)
ET_LOW = 15.0
ET_HIGH = 35.0

BLUE = "#1f77b4"
ORANGE = "#ff7f0e"
PP_BLUE = "#4fb3ff"
PP_ORANGE = "#ffbf66"
TEXT = "#202020"


def auc_from_scores(y_true: np.ndarray, score: np.ndarray) -> float:
    """Rank-based AUC without requiring sklearn."""
    sig = score[y_true == 1]
    bkg = score[y_true == 0]
    if sig.size == 0 or bkg.size == 0:
        return float("nan")
    values = np.concatenate([sig, bkg])
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(values.size, dtype=np.float64)
    ranks[order] = np.arange(1, values.size + 1, dtype=np.float64)

    # Average ranks for ties.
    sorted_values = values[order]
    start = 0
    while start < values.size:
        end = start + 1
        while end < values.size and sorted_values[end] == sorted_values[start]:
            end += 1
        if end - start > 1:
            ranks[order[start:end]] = 0.5 * (start + 1 + end)
        start = end

    rank_sum_sig = ranks[: sig.size].sum()
    return float((rank_sum_sig - sig.size * (sig.size + 1) / 2.0) / (sig.size * bkg.size))


def density_hist(scores: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    hist, edges = np.histogram(scores, bins=SCORE_BINS, density=True)
    return hist.astype(float), edges


def load_auau() -> dict[str, dict[str, object]]:
    paths = sorted((AUAU_REPORT / "score_caches").glob("score_cache_*.npz"))
    if not paths:
        raise SystemExit(f"No score caches found under {AUAU_REPORT / 'score_caches'}")

    accum: dict[str, dict[str, list[np.ndarray]]] = {
        label: {"signal": [], "background": []} for _, _, label in CENT_BINS
    }
    for path in paths:
        data = np.load(path, allow_pickle=True)
        score = data[SCORE_KEY].astype(float)
        y = data["is_signal"].astype(int)
        et = data["cluster_Et"].astype(float)
        cent = data["centrality"].astype(float)
        base = np.isfinite(score) & np.isfinite(et) & np.isfinite(cent) & (et > ET_LOW) & (et < ET_HIGH)
        for lo, hi, label in CENT_BINS:
            mask = base & (cent >= lo) & (cent < hi)
            if np.any(mask & (y == 1)):
                accum[label]["signal"].append(score[mask & (y == 1)])
            if np.any(mask & (y == 0)):
                accum[label]["background"].append(score[mask & (y == 0)])

    out: dict[str, dict[str, object]] = {}
    for _, _, label in CENT_BINS:
        sig = np.concatenate(accum[label]["signal"]) if accum[label]["signal"] else np.array([], dtype=float)
        bkg = np.concatenate(accum[label]["background"]) if accum[label]["background"] else np.array([], dtype=float)
        both_y = np.concatenate([np.ones(sig.size, dtype=int), np.zeros(bkg.size, dtype=int)])
        both_s = np.concatenate([sig, bkg])
        sig_h, edges = density_hist(sig)
        bkg_h, _ = density_hist(bkg)
        out[label] = {
            "signal": sig_h,
            "background": bkg_h,
            "edges": edges,
            "auc": auc_from_scores(both_y, both_s),
            "signal_entries": int(sig.size),
            "background_entries": int(bkg.size),
        }
    return out


def load_pp() -> dict[str, object]:
    rows: dict[str, list[dict[str, float]]] = {"signal": [], "background": []}
    with PP_CSV.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row["scope"] != "integrated":
                continue
            rows[row["class"]].append(
                {
                    "bin_low": float(row["bin_low"]),
                    "bin_high": float(row["bin_high"]),
                    "density": float(row["density"]),
                    "entries": int(float(row["entries"])),
                    "auc": float(row["auc"]),
                }
            )
    for vals in rows.values():
        vals.sort(key=lambda item: item["bin_low"])
    edges = np.array([rows["signal"][0]["bin_low"]] + [r["bin_high"] for r in rows["signal"]], dtype=float)
    return {
        "edges": edges,
        "signal": np.array([r["density"] for r in rows["signal"]], dtype=float),
        "background": np.array([r["density"] for r in rows["background"]], dtype=float),
        "signal_entries": rows["signal"][0]["entries"],
        "background_entries": rows["background"][0]["entries"],
        "auc": rows["signal"][0]["auc"],
    }


def step_xy(edges: np.ndarray, values: np.ndarray, scale: float = 1.0, offset: float = 0.0):
    x = np.repeat(edges, 2)[1:-1]
    y = np.repeat(values * scale + offset, 2)
    return x, y


def make_plot() -> dict[str, object]:
    auau = load_auau()
    pp = load_pp()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.titlesize": 16,
            "axes.labelsize": 14,
            "xtick.labelsize": 11,
            "ytick.labelsize": 11,
            "legend.fontsize": 10,
        }
    )

    fig, axes = plt.subplots(1, 3, figsize=(15.15, 5.01), sharey=True)
    fig.subplots_adjust(left=0.069, right=0.987, bottom=0.16, top=0.77, wspace=0.085)

    pp_ridge_top = 2.85
    pp_max = max(float(np.nanmax(pp["signal"])), float(np.nanmax(pp["background"])))
    pp_scale = pp_ridge_top / pp_max
    aucs = {}

    for ax, (_, _, label) in zip(axes, CENT_BINS):
        block = auau[label]
        edges = block["edges"]

        # pp reference ridge: shape-only band, scaled to avoid hiding Au+Au.
        ax.axhspan(0.0, pp_ridge_top, color="#f5f5f5", alpha=0.92, zorder=0)
        x_pp_sig, y_pp_sig = step_xy(pp["edges"], pp["signal"], pp_scale)
        x_pp_bkg, y_pp_bkg = step_xy(pp["edges"], pp["background"], pp_scale)
        ax.fill_between(x_pp_bkg, 0.0, y_pp_bkg, step="pre", color=PP_ORANGE, alpha=0.27, zorder=1)
        ax.fill_between(x_pp_sig, 0.0, y_pp_sig, step="pre", color=PP_BLUE, alpha=0.24, zorder=2)
        ax.step(pp["edges"][:-1], pp["background"] * pp_scale, where="post", color=PP_ORANGE, lw=1.8, ls=(0, (4, 2)), zorder=3)
        ax.step(pp["edges"][:-1], pp["signal"] * pp_scale, where="post", color=PP_BLUE, lw=1.8, ls=(0, (4, 2)), zorder=4)

        ax.step(edges[:-1], block["signal"], where="post", color=BLUE, lw=2.65, zorder=7)
        ax.step(edges[:-1], block["background"], where="post", color=ORANGE, lw=2.65, zorder=6)

        ax.set_xlim(0.0, 1.0)
        ax.set_ylim(0.0, 14.5)
        ax.set_xticks(np.linspace(0.0, 1.0, 6))
        ax.set_yticks(np.arange(0.0, 15.0, 2.0))
        ax.grid(True, color="#dddddd", alpha=0.62, lw=0.8)
        ax.set_axisbelow(True)
        ax.tick_params(direction="in", top=True, right=True, length=5)
        ax.set_xlabel("BDT score")
        ax.set_title(label, fontweight="bold", pad=12)
        aucs[label] = float(block["auc"])

        ax.text(
            0.047,
            0.885,
            f"Au+Au AUC {block['auc']:.3f}",
            transform=ax.transAxes,
            fontsize=11,
            weight="bold",
            bbox=dict(facecolor="white", edgecolor="#cfcfcf", boxstyle="round,pad=0.25", alpha=0.94),
            zorder=10,
        )
        ax.text(
            0.041,
            0.18,
            "pp ref band\nshape scaled",
            transform=ax.transAxes,
            fontsize=8.8,
            color="#555555",
            ha="left",
            va="center",
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.72, pad=2.0),
            zorder=10,
        )

    axes[0].set_ylabel("Area-normalized density")
    fig.text(0.071, 0.895, r"$\it{sPHENIX}$ Internal", fontsize=14, color=TEXT, ha="left", va="center")
    fig.suptitle("BDT score separation by coarse centrality", fontsize=19, fontweight="bold", y=0.968)
    fig.text(
        0.5,
        0.885,
        "Solid: Au+Au validated BDT score.  Shaded/dashed: pp baseline reference, shape-scaled.",
        ha="center",
        va="center",
        fontsize=10.8,
        color="#555555",
    )

    handles = [
        plt.Line2D([0], [0], color=BLUE, lw=3.0, label="Au+Au signal"),
        plt.Line2D([0], [0], color=ORANGE, lw=3.0, label="Au+Au background"),
        plt.Line2D([0], [0], color=PP_BLUE, lw=2.0, ls=(0, (4, 2)), label=f"pp signal ref, AUC {pp['auc']:.3f}"),
        plt.Line2D([0], [0], color=PP_ORANGE, lw=2.0, ls=(0, (4, 2)), label="pp background ref"),
    ]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.595, 0.065), ncol=4, frameon=False, fontsize=11)

    fig.savefig(OUT_PNG, dpi=170)
    plt.close(fig)

    meta = {
        "output_png": str(OUT_PNG),
        "auau_report": str(AUAU_REPORT),
        "auau_score_key": SCORE_KEY,
        "auau_selection": f"{ET_LOW:g} < cluster_Et < {ET_HIGH:g}, centrality 0-80 split into 0-20/20-50/50-80",
        "pp_source_csv": str(PP_CSV),
        "pp_reference": "ppg12_base_v1E_bdt_noIso integrated histogram; shape band scaled to avoid compressing Au+Au y-axis",
        "pp_scale_factor": pp_scale,
        "pp_auc": float(pp["auc"]),
        "pp_signal_entries": int(pp["signal_entries"]),
        "pp_background_entries": int(pp["background_entries"]),
        "auau_auc_by_centrality": aucs,
    }
    OUT_META.write_text(json.dumps(meta, indent=2) + "\n")
    return meta


def unit_shape(values: np.ndarray) -> np.ndarray:
    peak = float(np.nanmax(values)) if values.size else 0.0
    if peak <= 0.0 or not np.isfinite(peak):
        return values
    return values / peak


def make_unit_shape_plot() -> dict[str, object]:
    auau = load_auau()
    pp = load_pp()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.titlesize": 16,
            "axes.labelsize": 14,
            "xtick.labelsize": 11,
            "ytick.labelsize": 11,
            "legend.fontsize": 10,
        }
    )

    fig, axes = plt.subplots(1, 3, figsize=(15.15, 5.01), sharey=True)
    fig.subplots_adjust(left=0.069, right=0.987, bottom=0.16, top=0.77, wspace=0.085)

    pp_sig_unit = unit_shape(pp["signal"])
    pp_bkg_unit = unit_shape(pp["background"])
    aucs = {}

    for ax, (_, _, label) in zip(axes, CENT_BINS):
        block = auau[label]
        edges = block["edges"]
        auau_sig = unit_shape(block["signal"])
        auau_bkg = unit_shape(block["background"])

        ax.fill_between(
            np.repeat(pp["edges"], 2)[1:-1],
            0.0,
            np.repeat(pp_bkg_unit, 2),
            step="pre",
            color=PP_ORANGE,
            alpha=0.22,
            zorder=1,
        )
        ax.fill_between(
            np.repeat(pp["edges"], 2)[1:-1],
            0.0,
            np.repeat(pp_sig_unit, 2),
            step="pre",
            color=PP_BLUE,
            alpha=0.20,
            zorder=2,
        )
        ax.step(pp["edges"][:-1], pp_bkg_unit, where="post", color=PP_ORANGE, lw=2.0, ls=(0, (4, 2)), zorder=3)
        ax.step(pp["edges"][:-1], pp_sig_unit, where="post", color=PP_BLUE, lw=2.0, ls=(0, (4, 2)), zorder=4)
        ax.step(edges[:-1], auau_bkg, where="post", color=ORANGE, lw=2.75, zorder=6)
        ax.step(edges[:-1], auau_sig, where="post", color=BLUE, lw=2.75, zorder=7)

        ax.set_xlim(0.0, 1.0)
        ax.set_ylim(0.0, 1.12)
        ax.set_xticks(np.linspace(0.0, 1.0, 6))
        ax.set_yticks(np.linspace(0.0, 1.0, 6))
        ax.grid(True, color="#dddddd", alpha=0.62, lw=0.8)
        ax.set_axisbelow(True)
        ax.tick_params(direction="in", top=True, right=True, length=5)
        ax.set_xlabel("BDT score")
        ax.set_title(label, fontweight="bold", pad=12)
        aucs[label] = float(block["auc"])

        ax.text(
            0.047,
            0.885,
            f"Au+Au AUC {block['auc']:.3f}",
            transform=ax.transAxes,
            fontsize=11,
            weight="bold",
            bbox=dict(facecolor="white", edgecolor="#cfcfcf", boxstyle="round,pad=0.25", alpha=0.94),
            zorder=10,
        )
        ax.text(
            0.04,
            0.075,
            "each curve peak-normalized to 1",
            transform=ax.transAxes,
            fontsize=8.8,
            color="#555555",
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.72, pad=2.0),
            zorder=10,
        )

    axes[0].set_ylabel("Unit-normalized shape")
    fig.text(0.071, 0.895, r"$\it{sPHENIX}$ Internal", fontsize=14, color=TEXT, ha="left", va="center")
    fig.suptitle("BDT score separation by coarse centrality", fontsize=19, fontweight="bold", y=0.968)
    fig.text(
        0.5,
        0.885,
        "All curves are independently peak-normalized to show shape alignment; this is not a yield or density comparison.",
        ha="center",
        va="center",
        fontsize=10.8,
        color="#555555",
    )

    handles = [
        plt.Line2D([0], [0], color=BLUE, lw=3.0, label="Au+Au signal"),
        plt.Line2D([0], [0], color=ORANGE, lw=3.0, label="Au+Au background"),
        plt.Line2D([0], [0], color=PP_BLUE, lw=2.0, ls=(0, (4, 2)), label=f"pp signal ref, AUC {pp['auc']:.3f}"),
        plt.Line2D([0], [0], color=PP_ORANGE, lw=2.0, ls=(0, (4, 2)), label="pp background ref"),
    ]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.595, 0.065), ncol=4, frameon=False, fontsize=11)
    fig.savefig(OUT_UNIT_PNG, dpi=170)
    plt.close(fig)

    return {
        "output_png": str(OUT_UNIT_PNG),
        "normalization": "each curve divided by its own maximum bin height",
        "auau_auc_by_centrality": aucs,
        "pp_auc": float(pp["auc"]),
    }


def make_two_row_plot() -> dict[str, object]:
    auau = load_auau()
    pp = load_pp()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.titlesize": 13,
            "axes.labelsize": 12,
            "xtick.labelsize": 9.5,
            "ytick.labelsize": 9.5,
            "legend.fontsize": 9.5,
        }
    )

    fig, axes = plt.subplots(2, 3, figsize=(15.15, 7.7), sharex=True)
    fig.subplots_adjust(left=0.068, right=0.987, bottom=0.105, top=0.845, hspace=0.34, wspace=0.09)

    pp_ymax = 32.0
    auau_ymax = 14.5
    pp_col_titles = ["pp reference", "pp reference", "pp reference"]
    auau_aucs = {}

    for col, ax in enumerate(axes[0]):
        ax.step(pp["edges"][:-1], pp["signal"], where="post", color=BLUE, lw=2.15, zorder=3)
        ax.step(pp["edges"][:-1], pp["background"], where="post", color=ORANGE, lw=2.15, zorder=2)
        ax.set_title(pp_col_titles[col], fontweight="bold", pad=8)
        ax.text(
            0.045,
            0.82,
            f"pp AUC {pp['auc']:.3f}",
            transform=ax.transAxes,
            fontsize=10.2,
            weight="bold",
            bbox=dict(facecolor="white", edgecolor="#cfcfcf", boxstyle="round,pad=0.25", alpha=0.94),
        )
        ax.text(
            0.045,
            0.69,
            "integrated 15 < ET < 35 GeV\nno centrality axis",
            transform=ax.transAxes,
            fontsize=8.6,
            color="#555555",
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.76, pad=2),
        )
        ax.set_ylim(0.0, pp_ymax)
        ax.set_yticks([0, 8, 16, 24, 32])

    for ax, (_, _, label) in zip(axes[1], CENT_BINS):
        block = auau[label]
        edges = block["edges"]
        ax.step(edges[:-1], block["signal"], where="post", color=BLUE, lw=2.15, zorder=3)
        ax.step(edges[:-1], block["background"], where="post", color=ORANGE, lw=2.15, zorder=2)
        ax.set_title(label, fontweight="bold", pad=8)
        ax.text(
            0.045,
            0.82,
            f"Au+Au AUC {block['auc']:.3f}",
            transform=ax.transAxes,
            fontsize=10.2,
            weight="bold",
            bbox=dict(facecolor="white", edgecolor="#cfcfcf", boxstyle="round,pad=0.25", alpha=0.94),
        )
        ax.set_ylim(0.0, auau_ymax)
        ax.set_yticks(np.arange(0, 15, 2))
        ax.set_xlabel("BDT score")
        auau_aucs[label] = float(block["auc"])

    for ax in axes.flat:
        ax.set_xlim(0.0, 1.0)
        ax.set_xticks(np.linspace(0.0, 1.0, 6))
        ax.grid(True, color="#dddddd", alpha=0.62, lw=0.8)
        ax.set_axisbelow(True)
        ax.tick_params(direction="in", top=True, right=True, length=4)

    axes[0, 0].set_ylabel("pp area-normalized density")
    axes[1, 0].set_ylabel("Au+Au area-normalized density")
    axes[0, 0].text(0.03, 0.965, r"$\it{sPHENIX}$ Internal", transform=axes[0, 0].transAxes, fontsize=12.5, ha="left", va="top")

    fig.suptitle("BDT score separation: pp reference and Au+Au centrality panels", fontsize=19, fontweight="bold", y=0.968)
    fig.text(
        0.5,
        0.905,
        "Top: pp baseline BDT score separation. Bottom: same Au+Au slide-10 centrality panels. Curves are area-normalized within each class.",
        ha="center",
        va="center",
        fontsize=10.5,
        color="#555555",
    )
    handles = [
        plt.Line2D([0], [0], color=BLUE, lw=2.6, label="Signal"),
        plt.Line2D([0], [0], color=ORANGE, lw=2.6, label="Background"),
    ]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.058), ncol=2, frameon=False, fontsize=11)
    fig.savefig(OUT_TWO_ROW_PNG, dpi=170)
    plt.close(fig)

    return {
        "output_png": str(OUT_TWO_ROW_PNG),
        "normalization": "area-normalized density within each signal/background class",
        "pp_auc": float(pp["auc"]),
        "pp_signal_entries": int(pp["signal_entries"]),
        "pp_background_entries": int(pp["background_entries"]),
        "auau_auc_by_centrality": auau_aucs,
    }


if __name__ == "__main__":
    result = make_plot()
    result["unit_shape_overlay"] = make_unit_shape_plot()
    result["two_row_pp_top_auau_bottom"] = make_two_row_plot()
    OUT_META.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
