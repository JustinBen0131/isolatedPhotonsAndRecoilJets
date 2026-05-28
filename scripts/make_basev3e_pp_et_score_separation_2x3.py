#!/usr/bin/env python3
"""Make pp-vs-AuAu base-v3E score separation in matched ET panels."""

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
OUT_DIR = REPO / "dataOutput/ppPhotonMLPipeline/ppg12_fix4_score_overlays/score_separation/basev3e_pp_et_compare"
OUT_PNG = OUT_DIR / "pp_baseline_top_auau_basev3e_cent020_bottom_score_separation_2x3.png"
OUT_META = OUT_DIR / "pp_baseline_top_auau_basev3e_cent020_bottom_score_separation_2x3_metadata.json"

ET_BINS = [("15-20", 15.0, 20.0), ("20-25", 20.0, 25.0), ("25-35", 25.0, 35.0)]
BLUE = "#1f77b4"
ORANGE = "#ff7f0e"


def binned_auc(sig: np.ndarray, bkg: np.ndarray) -> float:
    sig_total = float(sig.sum())
    bkg_total = float(bkg.sum())
    if sig_total <= 0.0 or bkg_total <= 0.0:
        return float("nan")
    # Bins are sorted low score -> high score. AUC is P(signal score > background score).
    bkg_below = np.r_[0.0, np.cumsum(bkg)[:-1]]
    wins = float((sig * bkg_below).sum())
    ties = float((sig * bkg).sum())
    return (wins + 0.5 * ties) / (sig_total * bkg_total)


def step(ax, edges: np.ndarray, density: np.ndarray, color: str, label: str) -> None:
    y = np.r_[density, density[-1]]
    ax.step(edges, y, where="post", lw=2.25, color=color, label=label)


def load_pp() -> dict[str, dict[str, object]]:
    rows: dict[str, dict[str, list[dict[str, float]]]] = {
        label: {"signal": [], "background": []} for label, _, _ in ET_BINS
    }
    lines = PP_HIST.read_text().splitlines()
    try:
        header_idx = next(i for i, line in enumerate(lines) if line.startswith("scope,et_low,et_high,"))
    except StopIteration as exc:
        raise SystemExit(f"Could not find CSV header in {PP_HIST}") from exc
    for row in csv.DictReader(lines[header_idx:]):
        scope = row["scope"]
        if scope not in rows:
            continue
        rows[scope][row["class"]].append(
            {
                "bin_low": float(row["bin_low"]),
                "bin_high": float(row["bin_high"]),
                "density": float(row["density"]),
                "count": float(row["count"]),
                "entries": int(float(row["entries"])),
                "auc": float(row["auc"]),
            }
        )
    out: dict[str, dict[str, object]] = {}
    for label, _, _ in ET_BINS:
        for cls in ("signal", "background"):
            rows[label][cls].sort(key=lambda r: r["bin_low"])
        edges = np.array([rows[label]["signal"][0]["bin_low"]] + [r["bin_high"] for r in rows[label]["signal"]])
        out[label] = {
            "edges": edges,
            "signal_density": np.array([r["density"] for r in rows[label]["signal"]], dtype=float),
            "background_density": np.array([r["density"] for r in rows[label]["background"]], dtype=float),
            "signal_entries": rows[label]["signal"][0]["entries"],
            "background_entries": rows[label]["background"][0]["entries"],
            "auc": rows[label]["signal"][0]["auc"],
        }
    return out


def load_auau() -> dict[str, dict[str, object]]:
    df = pd.read_csv(AUAU_HIST)
    df = df[
        (df["product"] == "baseBDT_v3E_withCentrality")
        & (df["centrality_bin"] == "0-20")
        & (df["et_bin"].isin([label for label, _, _ in ET_BINS]))
    ].copy()
    if df.empty:
        raise SystemExit("No AuAu baseV3E centrality 0-20 rows found")
    out: dict[str, dict[str, object]] = {}
    for label, _, _ in ET_BINS:
        part = df[df["et_bin"] == label].sort_values("score_low")
        if part.empty:
            raise SystemExit(f"No AuAu rows for ET bin {label}")
        edges = np.r_[part["score_low"].to_numpy(float), float(part["score_high"].iloc[-1])]
        widths = part["score_high"].to_numpy(float) - part["score_low"].to_numpy(float)
        sig = part["signal_count"].to_numpy(float)
        bkg = part["background_count"].to_numpy(float)
        sig_total = float(sig.sum())
        bkg_total = float(bkg.sum())
        out[label] = {
            "edges": edges,
            "signal_density": np.divide(sig, sig_total * widths, out=np.zeros_like(sig), where=(sig_total * widths) > 0),
            "background_density": np.divide(bkg, bkg_total * widths, out=np.zeros_like(bkg), where=(bkg_total * widths) > 0),
            "signal_entries": int(sig_total),
            "background_entries": int(bkg_total),
            "auc": binned_auc(sig, bkg),
        }
    return out


def main() -> None:
    pp = load_pp()
    auau = load_auau()
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
    fig.subplots_adjust(left=0.07, right=0.987, bottom=0.108, top=0.845, hspace=0.34, wspace=0.09)

    pp_ymax = max(max(v["signal_density"].max(), v["background_density"].max()) for v in pp.values()) * 1.10
    auau_ymax = max(max(v["signal_density"].max(), v["background_density"].max()) for v in auau.values()) * 1.14
    meta = {"pp": {}, "auau_basev3e_cent020": {}, "et_bins": [b[0] for b in ET_BINS]}

    for col, (label, lo, hi) in enumerate(ET_BINS):
        top = pp[label]
        ax = axes[0, col]
        step(ax, top["edges"], top["signal_density"], BLUE, "Signal")
        step(ax, top["edges"], top["background_density"], ORANGE, "Background")
        ax.set_title(f"{label} GeV", fontweight="bold", pad=8)
        ax.set_ylim(0.0, pp_ymax)
        ax.text(
            0.045,
            0.82,
            f"pp AUC {top['auc']:.3f}",
            transform=ax.transAxes,
            fontsize=10.1,
            weight="bold",
            bbox=dict(facecolor="white", edgecolor="#cfcfcf", boxstyle="round,pad=0.25", alpha=0.94),
        )
        ax.text(
            0.045,
            0.68,
            "baseline pp BDT\nno centrality input",
            transform=ax.transAxes,
            fontsize=8.5,
            color="#555555",
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.76, pad=2),
        )
        meta["pp"][label] = {
            "auc": float(top["auc"]),
            "signal_entries": int(top["signal_entries"]),
            "background_entries": int(top["background_entries"]),
        }

        bottom = auau[label]
        ax = axes[1, col]
        step(ax, bottom["edges"], bottom["signal_density"], BLUE, "Signal")
        step(ax, bottom["edges"], bottom["background_density"], ORANGE, "Background")
        ax.set_title(f"0-20%, {label} GeV", fontweight="bold", pad=8)
        ax.set_ylim(0.0, auau_ymax)
        ax.set_xlabel("BDT score")
        ax.text(
            0.045,
            0.82,
            f"Au+Au AUC {bottom['auc']:.3f}",
            transform=ax.transAxes,
            fontsize=10.1,
            weight="bold",
            bbox=dict(facecolor="white", edgecolor="#cfcfcf", boxstyle="round,pad=0.25", alpha=0.94),
        )
        ax.text(
            0.045,
            0.68,
            "baseV3E + centrality\n0-20% only",
            transform=ax.transAxes,
            fontsize=8.5,
            color="#555555",
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.76, pad=2),
        )
        meta["auau_basev3e_cent020"][label] = {
            "auc_binned": float(bottom["auc"]),
            "signal_entries": int(bottom["signal_entries"]),
            "background_entries": int(bottom["background_entries"]),
        }

    for ax in axes.flat:
        ax.set_xlim(0.0, 1.0)
        ax.set_xticks(np.linspace(0.0, 1.0, 6))
        ax.grid(True, color="#dddddd", alpha=0.62, lw=0.8)
        ax.set_axisbelow(True)
        ax.tick_params(direction="in", top=True, right=True, length=4)

    axes[0, 0].set_ylabel("pp area-normalized density")
    axes[1, 0].set_ylabel("Au+Au area-normalized density")
    axes[0, 0].text(0.03, 0.965, r"$\it{sPHENIX}$ Internal", transform=axes[0, 0].transAxes, fontsize=12.5, ha="left", va="top")

    fig.suptitle("BDT score separation by photon ET", fontsize=19, fontweight="bold", y=0.968)
    fig.text(
        0.5,
        0.905,
        "Top: pp baseline no-centrality BDT. Bottom: Au+Au baseV3E + centrality BDT, 0-20% centrality only. Curves are area-normalized within each class.",
        ha="center",
        va="center",
        fontsize=10.3,
        color="#555555",
    )
    handles = [
        plt.Line2D([0], [0], color=BLUE, lw=2.6, label="Signal"),
        plt.Line2D([0], [0], color=ORANGE, lw=2.6, label="Background"),
    ]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.058), ncol=2, frameon=False, fontsize=11)
    fig.savefig(OUT_PNG, dpi=170)
    plt.close(fig)

    meta.update(
        {
            "output_png": str(OUT_PNG),
            "pp_hist_csv": str(PP_HIST),
            "auau_hist_csv": str(AUAU_HIST),
            "normalization": "signal and background are area-normalized separately inside each row/ET panel",
            "note": "ET binning follows available baseV3E validation table and exact pp cache extraction: 15-20, 20-25, 25-35 GeV.",
        }
    )
    OUT_META.write_text(json.dumps(meta, indent=2) + "\n")
    print(OUT_PNG)
    print(OUT_META)


if __name__ == "__main__":
    main()
