#!/usr/bin/env python3
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

from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(
    "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/auauMLDiagnosticRuns/"
    "global_etcent_inclusive3_sixpack_20260516_135439"
)
OUTDIR = ROOT / "slideReady" / "truth_bdt_performance_variants"
DIRECT_CSV = OUTDIR / "basev3e_cent_reference_vs_bdt_truth_purity_efficiency_direct_2x3.csv"

SUMMARY_CSV = OUTDIR / "basev3e_cent_reference_vs_bdt_truth_performance_slide_summary.csv"
SUMMARY_PNG = OUTDIR / "basev3e_cent_reference_vs_bdt_truth_performance_slide_summary.png"
ARROW_PNG = OUTDIR / "basev3e_cent_reference_vs_bdt_truth_operating_point_arrows.png"

CENT_ORDER = ["0_20", "20_50", "50_80"]
CENT_LABEL = {"0_20": "0-20%", "20_50": "20-50%", "50_80": "50-80%"}

REF = "Reference box cuts"
BDT_FULL = "Base v3E + centrality BDT target-80"
BDT_SHORT = "Base v3E + centrality BDT"

COL_REF = "#6B7280"
COL_BDT = "#009E73"
COL_CENT = {"0_20": "#3B82F6", "20_50": "#14B8A6", "50_80": "#F97316"}


def style_axis(ax, grid_axis: str = "y") -> None:
    ax.tick_params(direction="in", top=True, right=True, labelsize=11.5)
    ax.grid(True, axis=grid_axis, color="#D1D5DB", lw=0.8, alpha=0.75)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_linewidth(1.0)


def sphinx_label(fig, x: float, y: float, size: float = 15.0) -> None:
    fig.text(x, y, "sPHENIX", ha="right", fontsize=size, fontstyle="italic", fontweight="bold")
    fig.text(x + 0.092, y, "Internal", ha="right", fontsize=size)


def weighted_mean(values: np.ndarray, errors: np.ndarray) -> tuple[float, float]:
    values = np.asarray(values, dtype=float)
    errors = np.asarray(errors, dtype=float)
    good = np.isfinite(values) & np.isfinite(errors) & (errors > 0)
    if not np.any(good):
        return float(np.nanmean(values)), float(np.nanstd(values) / np.sqrt(max(len(values), 1)))
    weights = 1.0 / np.square(errors[good])
    mean = float(np.sum(values[good] * weights) / np.sum(weights))
    err = float(np.sqrt(1.0 / np.sum(weights)))
    return mean, err


def build_summary() -> pd.DataFrame:
    if not DIRECT_CSV.exists():
        raise SystemExit(f"Missing input CSV: {DIRECT_CSV}")
    df = pd.read_csv(DIRECT_CSV)
    rows = []
    for cent in CENT_ORDER:
        for selection in [REF, BDT_FULL]:
            g = df[(df["centrality"] == cent) & (df["selection"] == selection)].copy()
            if len(g) != 8:
                raise SystemExit(f"Expected 8 ET rows for {cent} / {selection}, found {len(g)}")
            purity, purity_err = weighted_mean(g["truth_purity"].to_numpy(), g["truth_purity_error"].to_numpy())
            recovery, recovery_err = weighted_mean(
                g["truth_efficiency"].to_numpy(), g["truth_efficiency_error"].to_numpy()
            )
            rows.append(
                {
                    "centrality": cent,
                    "centrality_label": CENT_LABEL[cent],
                    "selection": BDT_SHORT if selection == BDT_FULL else REF,
                    "truth_purity_weighted_mean": purity,
                    "truth_purity_weighted_error": purity_err,
                    "truth_recovery_weighted_mean": recovery,
                    "truth_recovery_weighted_error": recovery_err,
                    "n_et_bins": len(g),
                    "et_range": "15-35 GeV",
                    "summary_method": "inverse-variance weighted average of direct 8-bin points; no threshold retuning",
                }
            )

    out = pd.DataFrame(rows)
    ref = out[out["selection"] == REF].set_index("centrality")
    for idx, row in out.iterrows():
        base = ref.loc[row["centrality"]]
        out.loc[idx, "purity_delta_vs_reference"] = (
            row["truth_purity_weighted_mean"] - base["truth_purity_weighted_mean"]
        )
        out.loc[idx, "recovery_delta_vs_reference"] = (
            row["truth_recovery_weighted_mean"] - base["truth_recovery_weighted_mean"]
        )
    return out


def pct(x: float) -> str:
    return f"{100.0 * x:.1f}%"


def pp(x: float) -> str:
    return f"{100.0 * x:+.1f} pp"


def plot_summary(summary: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(13.7, 5.75), dpi=190, gridspec_kw={"wspace": 0.24})

    x = np.arange(len(CENT_ORDER))
    width = 0.34
    ref = summary[summary["selection"] == REF].set_index("centrality").loc[CENT_ORDER]
    bdt = summary[summary["selection"] == BDT_SHORT].set_index("centrality").loc[CENT_ORDER]

    ax = axes[0]
    ax.bar(
        x - width / 2,
        ref["truth_recovery_weighted_mean"],
        width,
        yerr=ref["truth_recovery_weighted_error"],
        capsize=3,
        color=COL_REF,
        alpha=0.82,
        label=REF,
    )
    ax.bar(
        x + width / 2,
        bdt["truth_recovery_weighted_mean"],
        width,
        yerr=bdt["truth_recovery_weighted_error"],
        capsize=3,
        color=COL_BDT,
        alpha=0.9,
        label=BDT_SHORT,
    )
    ax.set_ylim(0.0, 0.53)
    ax.set_ylabel("Truth-photon recovery", fontsize=13.0)
    ax.set_xticks(x, [CENT_LABEL[c] for c in CENT_ORDER], fontsize=12.5)
    ax.set_title("Truth-photon recovery", fontsize=14.8, fontweight="bold", pad=10)
    style_axis(ax)
    for i, cent in enumerate(CENT_ORDER):
        r = float(ref.loc[cent, "truth_recovery_weighted_mean"])
        b = float(bdt.loc[cent, "truth_recovery_weighted_mean"])
        ax.text(i - width / 2, r + 0.017, pct(r), ha="center", va="bottom", fontsize=11.3, color="#374151")
        ax.text(i + width / 2, b + 0.017, pct(b), ha="center", va="bottom", fontsize=11.3, color="#065F46")
        ax.text(
            i,
            max(r, b) + 0.063,
            pp(b - r),
            ha="center",
            va="bottom",
            fontsize=11.4,
            fontweight="bold",
            color="#065F46",
        )

    ax = axes[1]
    for i, cent in enumerate(CENT_ORDER):
        r = float(ref.loc[cent, "truth_purity_weighted_mean"])
        b = float(bdt.loc[cent, "truth_purity_weighted_mean"])
        re = float(ref.loc[cent, "truth_purity_weighted_error"])
        be = float(bdt.loc[cent, "truth_purity_weighted_error"])
        ax.plot([i - 0.14, i + 0.14], [r, b], color="#9CA3AF", lw=2.0, zorder=1)
        ax.errorbar(i - 0.14, r, yerr=re, fmt="o", ms=8.0, capsize=3, color=COL_REF, zorder=3)
        ax.errorbar(i + 0.14, b, yerr=be, fmt="s", ms=8.0, capsize=3, color=COL_BDT, zorder=4)
        ax.text(i - 0.14, r + 0.0045, pct(r), ha="center", va="bottom", fontsize=10.8, color="#374151")
        ax.text(i + 0.14, b + 0.0045, pct(b), ha="center", va="bottom", fontsize=10.8, color="#065F46")
        ax.text(
            i,
            0.846,
            f"{pp(b - r)}",
            ha="center",
            va="center",
            fontsize=10.8,
            fontweight="bold",
            color="#065F46" if b >= r else "#B45309",
        )
    ax.axhspan(0.875, 0.905, color="#ECFDF5", alpha=0.85, zorder=0)
    ax.text(
        2.58,
        0.903,
        "high-purity band",
        ha="right",
        va="top",
        fontsize=10.8,
        color="#047857",
    )
    ax.set_ylim(0.84, 0.915)
    ax.set_xlim(-0.55, 2.55)
    ax.set_ylabel("Truth purity of selected photons", fontsize=13.0)
    ax.set_xticks(x, [CENT_LABEL[c] for c in CENT_ORDER], fontsize=12.5)
    ax.set_title("Truth purity", fontsize=14.8, fontweight="bold", pad=10)
    style_axis(ax)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.50, 0.928),
        ncol=2,
        frameon=False,
        fontsize=12.2,
        handlelength=1.4,
    )
    fig.text(
        0.075,
        0.955,
        "15 < cluster $E_T$ < 35 GeV; existing selections, no efficiency matching or score retuning",
        fontsize=12.0,
        color="#4B5563",
    )
    sphinx_label(fig, 0.86, 0.955, size=13.6)
    fig.subplots_adjust(left=0.08, right=0.985, bottom=0.14, top=0.79, wspace=0.24)
    fig.savefig(SUMMARY_PNG)
    plt.close(fig)


def plot_arrows(summary: pd.DataFrame) -> None:
    ref = summary[summary["selection"] == REF].set_index("centrality").loc[CENT_ORDER]
    bdt = summary[summary["selection"] == BDT_SHORT].set_index("centrality").loc[CENT_ORDER]

    fig, ax = plt.subplots(figsize=(10.8, 6.4), dpi=190)
    ax.axhspan(0.875, 0.905, color="#ECFDF5", alpha=0.95, zorder=0)
    ax.text(0.455, 0.903, "high-purity selected-photon band", ha="right", va="top", fontsize=11.2, color="#047857")

    for cent in CENT_ORDER:
        x0 = float(ref.loc[cent, "truth_recovery_weighted_mean"])
        y0 = float(ref.loc[cent, "truth_purity_weighted_mean"])
        x1 = float(bdt.loc[cent, "truth_recovery_weighted_mean"])
        y1 = float(bdt.loc[cent, "truth_purity_weighted_mean"])
        color = COL_CENT[cent]
        ax.scatter(x0, y0, s=92, marker="o", color=COL_REF, edgecolor="white", linewidth=1.1, zorder=3)
        ax.scatter(x1, y1, s=110, marker="s", color=color, edgecolor="#111827", linewidth=0.8, zorder=4)
        ax.annotate(
            "",
            xy=(x1, y1),
            xytext=(x0, y0),
            arrowprops=dict(arrowstyle="-|>", color=color, lw=3.0, mutation_scale=14, shrinkA=7, shrinkB=8),
            zorder=2,
        )
        ax.text(
            x1 + 0.006,
            y1 + (0.004 if cent != "0_20" else -0.006),
            f"{CENT_LABEL[cent]}",
            fontsize=12.3,
            fontweight="bold",
            color=color,
            va="center",
        )

    ax.scatter([], [], s=92, marker="o", color=COL_REF, edgecolor="white", linewidth=1.1, label=REF)
    ax.scatter([], [], s=110, marker="s", color=COL_BDT, edgecolor="#111827", linewidth=0.8, label=BDT_SHORT)
    ax.legend(loc="lower right", frameon=False, fontsize=12.1)

    ax.set_xlim(0.20, 0.47)
    ax.set_ylim(0.855, 0.912)
    ax.set_xlabel("Truth-photon recovery", fontsize=13.5)
    ax.set_ylabel("Truth purity of selected photons", fontsize=13.5)
    ax.set_title("BDT moves the selection toward more recovered photons", fontsize=17.0, fontweight="bold", pad=14)
    ax.text(
        0.202,
        0.909,
        "Each arrow: reference box cuts -> base v3E + centrality BDT target-80",
        ha="left",
        va="top",
        fontsize=11.4,
        color="#4B5563",
    )
    style_axis(ax, grid_axis="both")
    sphinx_label(fig, 0.80, 0.942, size=14.2)
    fig.tight_layout(rect=[0.04, 0.05, 0.98, 0.92])
    fig.savefig(ARROW_PNG)
    plt.close(fig)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    summary = build_summary()
    summary.to_csv(SUMMARY_CSV, index=False)
    plot_summary(summary)
    plot_arrows(summary)
    print(f"Wrote {SUMMARY_CSV}")
    print(f"Wrote {SUMMARY_PNG}")
    print(f"Wrote {ARROW_PNG}")


if __name__ == "__main__":
    main()
