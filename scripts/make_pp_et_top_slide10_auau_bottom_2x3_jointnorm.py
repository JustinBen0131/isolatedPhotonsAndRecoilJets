#!/usr/bin/env python3
"""Make the pp/AuAu 2x3 BDT score plot with signal+bkg joint normalization."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402


REPO = Path(__file__).resolve().parents[1]
PP_HIST = (
    REPO
    / "dataOutput/ppPhotonMLPipeline/ppg12_fix4_score_overlays/score_separation/exact_et_bins/"
    / "pp_baseline_bdt_score_histograms_15_20_20_25_25_35.csv"
)
AUAU_HIST = (
    REPO
    / "dataOutput/auauMLDiagnosticRuns/global_etcent_inclusive3_sixpack_20260516_135439/"
    / "validation/basev3e_controls_20260518_1110/basev3e_with_centrality_score_histograms.csv"
)
AUAU_AUC = (
    REPO
    / "dataOutput/auauMLDiagnosticRuns/global_etcent_inclusive3_sixpack_20260516_135439/"
    / "slideReady/basev3e_roc/baseBDT_v3E_withCentrality_roc_auc_override_for_score_split.csv"
)
OUT_DIR = REPO / "dataOutput/ppPhotonMLPipeline/ppg12_fix4_score_overlays/score_separation/slide10_pp_top"
OUT_PNG = OUT_DIR / "pp_et_bins_top_slide10_auau_baseline_bottom_2x3_joint_signal_bkg_norm.png"
OUT_META = OUT_DIR / "pp_et_bins_top_slide10_auau_baseline_bottom_2x3_joint_signal_bkg_norm_metadata.json"

PP_ET_BINS = [("15-20", 15.0, 20.0), ("20-25", 20.0, 25.0), ("25-35", 25.0, 35.0)]
AUAU_CENT_BINS = [("0-20", 0.0, 20.0), ("20-50", 20.0, 50.0), ("50-80", 50.0, 80.0)]
BLUE = "#1f77b4"
ORANGE = "#ff7f0e"


def auc_from_counts(sig: np.ndarray, bkg: np.ndarray) -> float:
    n_sig = float(sig.sum())
    n_bkg = float(bkg.sum())
    if n_sig <= 0.0 or n_bkg <= 0.0:
        return float("nan")
    bkg_below = np.r_[0.0, np.cumsum(bkg)[:-1]]
    return float(((sig * bkg_below).sum() + 0.5 * (sig * bkg).sum()) / (n_sig * n_bkg))


def joint_density(counts: np.ndarray, total: float, widths: np.ndarray) -> np.ndarray:
    denom = total * widths
    return np.divide(counts, denom, out=np.zeros_like(counts), where=denom > 0.0)


def step(ax, edges: np.ndarray, density: np.ndarray, color: str) -> None:
    ax.step(edges, np.r_[density, density[-1]], where="post", lw=2.2, color=color)


def read_pp_hist() -> dict[str, dict[str, object]]:
    lines = PP_HIST.read_text().splitlines()
    header_idx = next(i for i, line in enumerate(lines) if line.startswith("scope,et_low,et_high,"))
    by_bin: dict[str, dict[str, list[dict[str, float]]]] = {
        label: {"signal": [], "background": []} for label, _, _ in PP_ET_BINS
    }
    for row in csv.DictReader(lines[header_idx:]):
        scope = row["scope"]
        if scope not in by_bin:
            continue
        by_bin[scope][row["class"]].append(
            {
                "bin_low": float(row["bin_low"]),
                "bin_high": float(row["bin_high"]),
                "count": float(row["count"]),
                "entries": int(float(row["entries"])),
                "auc": float(row["auc"]),
            }
        )

    out: dict[str, dict[str, object]] = {}
    for label, _, _ in PP_ET_BINS:
        for cls in ("signal", "background"):
            by_bin[label][cls].sort(key=lambda item: item["bin_low"])
        sig_rows = by_bin[label]["signal"]
        bkg_rows = by_bin[label]["background"]
        edges = np.array([sig_rows[0]["bin_low"]] + [item["bin_high"] for item in sig_rows], dtype=float)
        widths = np.diff(edges)
        sig = np.array([item["count"] for item in sig_rows], dtype=float)
        bkg = np.array([item["count"] for item in bkg_rows], dtype=float)
        total = float(sig.sum() + bkg.sum())
        out[label] = {
            "edges": edges,
            "signal_density": joint_density(sig, total, widths),
            "background_density": joint_density(bkg, total, widths),
            "signal_entries": int(sig_rows[0]["entries"]),
            "background_entries": int(bkg_rows[0]["entries"]),
            "signal_integral": float(sig.sum() / total) if total > 0.0 else float("nan"),
            "background_integral": float(bkg.sum() / total) if total > 0.0 else float("nan"),
            "auc": float(sig_rows[0]["auc"]),
        }
    return out


def read_auau_auc_overrides() -> dict[str, float]:
    if not AUAU_AUC.is_file():
        return {}
    df = pd.read_csv(AUAU_AUC)
    needed = {"model", "centrality_bin", "auc_entries_weighted"}
    if not needed.issubset(df.columns):
        return {}
    out = {}
    for row in df[df["model"] == "baseBDT_v3E_withCentrality"].to_dict("records"):
        out[str(row["centrality_bin"]).replace("_", "-")] = float(row["auc_entries_weighted"])
    return out


def read_auau_slide10() -> dict[str, dict[str, object]]:
    df = pd.read_csv(AUAU_HIST)
    base = df[
        (df["product"] == "baseBDT_v3E_withCentrality")
        & (df["et_low"] >= 15.0)
        & (df["et_high"] <= 35.0)
    ].copy()
    if base.empty:
        raise SystemExit(f"No baseBDT_v3E_withCentrality rows found in {AUAU_HIST}")
    auc_override = read_auau_auc_overrides()

    out: dict[str, dict[str, object]] = {}
    for label, lo, hi in AUAU_CENT_BINS:
        part = base[(base["centrality_low"] == lo) & (base["centrality_high"] == hi)].copy()
        if part.empty:
            part = base[(base["centrality_low"] >= lo) & (base["centrality_high"] <= hi)].copy()
        if part.empty:
            raise SystemExit(f"No AuAu rows for centrality {label}")
        grouped = (
            part.groupby(["score_low", "score_high"], as_index=False)[["signal_count", "background_count"]]
            .sum()
            .sort_values("score_low")
        )
        edges = np.r_[grouped["score_low"].to_numpy(float), float(grouped["score_high"].iloc[-1])]
        widths = np.diff(edges)
        sig = grouped["signal_count"].to_numpy(float)
        bkg = grouped["background_count"].to_numpy(float)
        total = float(sig.sum() + bkg.sum())
        out[label] = {
            "edges": edges,
            "signal_density": joint_density(sig, total, widths),
            "background_density": joint_density(bkg, total, widths),
            "signal_entries": int(round(sig.sum())),
            "background_entries": int(round(bkg.sum())),
            "signal_integral": float(sig.sum() / total) if total > 0.0 else float("nan"),
            "background_integral": float(bkg.sum() / total) if total > 0.0 else float("nan"),
            "auc": float(auc_override.get(label, auc_from_counts(sig, bkg))),
            "pooled_auc": auc_from_counts(sig, bkg),
        }
    return out


def main() -> None:
    pp = read_pp_hist()
    auau = read_auau_slide10()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.titlesize": 13,
            "axes.labelsize": 12,
            "xtick.labelsize": 9.5,
            "ytick.labelsize": 9.5,
            "legend.fontsize": 10,
        }
    )
    fig, axes = plt.subplots(2, 3, figsize=(15.15, 7.7), sharex=True)
    fig.subplots_adjust(left=0.07, right=0.987, bottom=0.106, top=0.88, hspace=0.34, wspace=0.09)

    pp_ymax = max(max(v["signal_density"].max(), v["background_density"].max()) for v in pp.values()) * 1.12
    auau_ymax = max(max(v["signal_density"].max(), v["background_density"].max()) for v in auau.values()) * 1.18

    meta = {"pp": {}, "auau_slide10_basev3e_withCentrality": {}}
    for col, (label, _, _) in enumerate(PP_ET_BINS):
        ax = axes[0, col]
        block = pp[label]
        step(ax, block["edges"], block["signal_density"], BLUE)
        step(ax, block["edges"], block["background_density"], ORANGE)
        ax.set_title(f"pp {label} GeV", fontweight="bold", pad=8)
        ax.set_ylim(0.0, pp_ymax)
        ax.text(
            0.045,
            0.82,
            f"AUC {block['auc']:.3f}",
            transform=ax.transAxes,
            fontsize=10.2,
            weight="bold",
            bbox=dict(facecolor="white", edgecolor="#cfcfcf", boxstyle="round,pad=0.25", alpha=0.94),
        )
        ax.text(
            0.045,
            0.68,
            "baseline pp BDT\njoint signal+bkg norm",
            transform=ax.transAxes,
            fontsize=8.5,
            color="#555555",
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.76, pad=2),
        )
        meta["pp"][label] = {
            "auc": float(block["auc"]),
            "signal_entries": int(block["signal_entries"]),
            "background_entries": int(block["background_entries"]),
            "signal_integral": float(block["signal_integral"]),
            "background_integral": float(block["background_integral"]),
        }

    for col, (label, _, _) in enumerate(AUAU_CENT_BINS):
        ax = axes[1, col]
        block = auau[label]
        step(ax, block["edges"], block["signal_density"], BLUE)
        step(ax, block["edges"], block["background_density"], ORANGE)
        ax.set_title(rf"{label}%, $p_T = 15$-$35$ GeV", fontweight="bold", pad=8)
        ax.set_ylim(0.0, auau_ymax)
        ax.set_xlabel("BDT score")
        ax.text(
            0.045,
            0.82,
            f"AUC {block['auc']:.3f}",
            transform=ax.transAxes,
            fontsize=10.2,
            weight="bold",
            bbox=dict(facecolor="white", edgecolor="#cfcfcf", boxstyle="round,pad=0.25", alpha=0.94),
        )
        ax.text(
            0.045,
            0.68,
            "Au+Au slide-10 baseline\njoint signal+bkg norm",
            transform=ax.transAxes,
            fontsize=8.5,
            color="#555555",
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.76, pad=2),
        )
        meta["auau_slide10_basev3e_withCentrality"][label] = {
            "auc_display": float(block["auc"]),
            "auc_pooled_histogram": float(block["pooled_auc"]),
            "signal_entries": int(block["signal_entries"]),
            "background_entries": int(block["background_entries"]),
            "signal_integral": float(block["signal_integral"]),
            "background_integral": float(block["background_integral"]),
        }

    for ax in axes.flat:
        ax.set_xlim(0.0, 1.0)
        ax.set_xticks(np.linspace(0.0, 1.0, 6))
        ax.grid(True, color="#dddddd", alpha=0.62, lw=0.8)
        ax.set_axisbelow(True)
        ax.tick_params(direction="in", top=True, right=True, length=4)

    axes[0, 0].set_ylabel("pp joint-normalized density")
    axes[1, 0].set_ylabel("Au+Au joint-normalized density")
    axes[0, 0].text(
        0.03,
        0.965,
        r"$\it{sPHENIX}$ Internal",
        transform=axes[0, 0].transAxes,
        fontsize=12.5,
        ha="left",
        va="top",
    )

    fig.suptitle("BDT score separation", fontsize=19, fontweight="bold", y=0.965)
    handles = [
        plt.Line2D([0], [0], color=BLUE, lw=2.6, label="Signal"),
        plt.Line2D([0], [0], color=ORANGE, lw=2.6, label="Background"),
    ]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.057), ncol=2, frameon=False, fontsize=11)
    fig.savefig(OUT_PNG, dpi=170)
    plt.close(fig)

    meta.update(
        {
            "output_png": str(OUT_PNG),
            "pp_hist_csv": str(PP_HIST),
            "auau_hist_csv": str(AUAU_HIST),
            "auau_auc_override_csv": str(AUAU_AUC),
            "normalization": (
                "signal and background counts are normalized together in every panel; "
                "signal area plus background area equals one"
            ),
        }
    )
    OUT_META.write_text(json.dumps(meta, indent=2) + "\n")
    print(OUT_PNG)
    print(OUT_META)


if __name__ == "__main__":
    main()
