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

PURITY_OVERLAY = (
    ROOT
    / "slideReady"
    / "basev3e_reco_efficiency_overlay"
    / "basev3e_target80_vs_reference_boxcuts_raw_sim_purity_recoPt15to35_points.csv"
)
EFFICIENCY_OVERLAY = (
    ROOT
    / "slideReady"
    / "basev3e_reco_efficiency_overlay"
    / "basev3e_target80_vs_reference_boxcuts_reco_efficiency_pt15to35_points.csv"
)
REFERENCE_TRUTH_COUNTS = (
    ROOT
    / "slideReady"
    / "abcd_raw_purity_reference_vs_basev3e"
    / "reference_boxcuts_vs_basev3e_centrality_target80_reco_abcd_raw_purity_pt15to35_points.csv"
)
BDT_TRUTH_COUNTS = (
    ROOT
    / "slideReady"
    / "basev3e_reco_efficiency_overlay"
    / "basev3e_target80_raw_sim_purity_recoPt15to35_points.csv"
)

PNG_OUT = OUTDIR / "basev3e_cent_reference_vs_bdt_truth_purity_efficiency_direct_2x3.png"
CSV_OUT = OUTDIR / "basev3e_cent_reference_vs_bdt_truth_purity_efficiency_direct_2x3.csv"

CENT_ORDER = ["0_20", "20_50", "50_80"]
CENT_LABEL = {"0_20": "0-20%", "20_50": "20-50%", "50_80": "50-80%"}
ET_BINS = [
    ("15-17", 15.0, 17.0, 16.0),
    ("17-19", 17.0, 19.0, 18.0),
    ("19-21", 19.0, 21.0, 20.0),
    ("21-23", 21.0, 23.0, 22.0),
    ("23-25", 23.0, 25.0, 24.0),
    ("25-27", 25.0, 27.0, 26.0),
    ("27-30", 27.0, 30.0, 28.5),
    ("30-35", 30.0, 35.0, 32.5),
]


def style_ax(ax, grid_axis: str = "both") -> None:
    ax.tick_params(direction="in", top=True, right=True, labelsize=10.5)
    ax.grid(True, axis=grid_axis, color="#D1D5DB", lw=0.7, alpha=0.75)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_linewidth(1.0)


def sphinx_label(fig, x=0.865, y=0.957, size=14.5) -> None:
    fig.text(x, y, "sPHENIX", ha="right", fontsize=size, fontstyle="italic", fontweight="bold")
    fig.text(x + 0.095, y, "Internal", ha="right", fontsize=size)


def overlap_fraction(source_low: float, source_high: float, target_low: float, target_high: float) -> float:
    overlap = max(0.0, min(source_high, target_high) - max(source_low, target_low))
    return overlap / max(source_high - source_low, 1e-12)


def ratio_and_error(signal: float, background: float) -> tuple[float, float]:
    denom = signal + background
    if denom <= 0.0:
        return np.nan, np.nan
    value = signal / denom
    # Binomial approximation on weighted truth-labeled selected counts. This is
    # used only for the rebinned reference purity points.
    error = np.sqrt(max(value * (1.0 - value) / denom, 0.0))
    return float(value), float(error)


def load_reference_purity_rebinned() -> pd.DataFrame:
    df = pd.read_csv(REFERENCE_TRUTH_COUNTS)
    df = df[(df["model"] == "Reference box cuts") & (df["quantity"] == "truth_labeled_A_purity")].copy()
    rows = []
    for cent in CENT_ORDER:
        gcent = df[df["centrality"].astype(str) == cent].copy()
        for et_label, et_low, et_high, et_mid in ET_BINS:
            signal = 0.0
            background = 0.0
            contributing_bins = []
            for _, row in gcent.iterrows():
                frac = overlap_fraction(float(row["pt_low"]), float(row["pt_high"]), et_low, et_high)
                if frac <= 0.0:
                    continue
                signal += float(row["signal_tight"]) * frac
                background += float(row["background_tight"]) * frac
                contributing_bins.append(f"{row['pt_low']:g}-{row['pt_high']:g}:{frac:.3g}")
            value, error = ratio_and_error(signal, background)
            rows.append(
                {
                    "selection": "Reference box cuts",
                    "centrality": cent,
                    "centrality_label": CENT_LABEL[cent],
                    "et_bin": et_label,
                    "et_low": et_low,
                    "et_high": et_high,
                    "et_mid": et_mid,
                    "truth_purity": value,
                    "truth_purity_error": error,
                    "purity_provenance": "reference truth-count rows rebinned to BDT 8-bin ET grid by bin-overlap fractions",
                    "purity_source_bins": ";".join(contributing_bins),
                }
            )
    return pd.DataFrame(rows)


def load_bdt_purity() -> pd.DataFrame:
    counts = pd.read_csv(BDT_TRUTH_COUNTS)
    counts = counts[counts["source"] == "recoiljets_final_root"].copy()
    rows = []
    for _, row in counts.iterrows():
        rows.append(
            {
                    "selection": "Base v3E + centrality BDT target-80",
                "centrality": str(row["centrality"]),
                "centrality_label": CENT_LABEL[str(row["centrality"])],
                "et_bin": f"{int(row['pt_low'])}-{int(row['pt_high'])}",
                "et_low": float(row["pt_low"]),
                "et_high": float(row["pt_high"]),
                "et_mid": float(row["pt_mid"]),
                "truth_purity": float(row["value"]),
                "truth_purity_error": float(row["error"]),
                "purity_provenance": "base v3E + centrality target-80 RecoilJets truth-count QA",
                "purity_source_bins": f"{row['pt_low']:g}-{row['pt_high']:g}",
            }
        )
    return pd.DataFrame(rows)


def load_efficiency() -> pd.DataFrame:
    df = pd.read_csv(EFFICIENCY_OVERLAY)
    label_map = {
        "reference_box_cuts": "Reference box cuts",
        "basev3e_centrality_target80": "Base v3E + centrality BDT target-80",
    }
    rows = []
    for _, row in df.iterrows():
        source = str(row["source"])
        if source not in label_map:
            continue
        mid = float(row["pt_mid"])
        et_match = [b for b in ET_BINS if abs(b[3] - mid) < 1e-9]
        if not et_match:
            continue
        et_label, et_low, et_high, et_mid = et_match[0]
        rows.append(
            {
                "selection": label_map[source],
                "centrality": str(row["centrality"]),
                "et_bin": et_label,
                "et_low": et_low,
                "et_high": et_high,
                "et_mid": et_mid,
                "truth_efficiency": float(row["value"]),
                "truth_efficiency_error": float(row["error"]),
                "efficiency_provenance": "existing truth-photon recovery QA; no threshold rematching",
            }
        )
    return pd.DataFrame(rows)


def build_table() -> pd.DataFrame:
    purity = pd.concat([load_reference_purity_rebinned(), load_bdt_purity()], ignore_index=True)
    efficiency = load_efficiency()
    df = purity.merge(
        efficiency,
        on=["selection", "centrality", "et_bin", "et_low", "et_high", "et_mid"],
        how="inner",
    )
    df["centrality_label"] = df["centrality"].map(CENT_LABEL)
    df["purity_minus_reference"] = np.nan
    df["efficiency_minus_reference"] = np.nan

    ref = df[df["selection"] == "Reference box cuts"].set_index(["centrality", "et_bin"])
    for idx, row in df.iterrows():
        key = (row["centrality"], row["et_bin"])
        if key not in ref.index:
            continue
        df.loc[idx, "purity_minus_reference"] = row["truth_purity"] - ref.loc[key, "truth_purity"]
        df.loc[idx, "efficiency_minus_reference"] = row["truth_efficiency"] - ref.loc[key, "truth_efficiency"]

    expected = 2 * len(CENT_ORDER) * len(ET_BINS)
    if len(df) != expected:
        raise SystemExit(f"Expected {expected} aligned rows, found {len(df)}")
    return df.sort_values(["centrality", "et_low", "selection"]).reset_index(drop=True)


def plot(df: pd.DataFrame) -> None:
    colors = {"Reference box cuts": "#6B7280", "Base v3E + centrality BDT target-80": "#009E73"}
    markers = {"Reference box cuts": "o", "Base v3E + centrality BDT target-80": "s"}
    labels = [b[0] for b in ET_BINS]
    x = np.arange(len(labels))

    fig, axes = plt.subplots(
        2,
        3,
        figsize=(17.2, 7.9),
        dpi=180,
        sharex=True,
        gridspec_kw={"height_ratios": [1.15, 1.0], "hspace": 0.12, "wspace": 0.17},
    )

    for col, cent in enumerate(CENT_ORDER):
        gcent = df[df["centrality"] == cent]
        for selection in ("Reference box cuts", "Base v3E + centrality BDT target-80"):
            g = gcent[gcent["selection"] == selection].set_index("et_bin").loc[labels].reset_index()
            axes[0, col].errorbar(
                x,
                g["truth_purity"],
                yerr=g["truth_purity_error"],
                fmt=f"{markers[selection]}-",
                ms=4.9,
                lw=1.8 if selection == "Reference box cuts" else 2.2,
                capsize=2.2,
                color=colors[selection],
                label=selection,
                zorder=3 if selection == "Reference box cuts" else 4,
            )
            axes[1, col].errorbar(
                x,
                g["truth_efficiency"],
                yerr=g["truth_efficiency_error"],
                fmt=f"{markers[selection]}-",
                ms=4.9,
                lw=1.8 if selection == "Reference box cuts" else 2.2,
                capsize=2.2,
                color=colors[selection],
                label=selection,
                zorder=3 if selection == "Reference box cuts" else 4,
            )

        axes[0, col].set_title(f"{CENT_LABEL[cent]} centrality", fontsize=13.4, fontweight="bold")
        axes[0, col].set_ylim(0.82, 0.93)
        axes[1, col].set_ylim(0.08, 0.56)
        axes[1, col].set_xticks(x, labels, rotation=35, ha="right")
        axes[1, col].set_xlabel(r"Cluster $E_T$ bin [GeV]")
        if col == 0:
            axes[0, col].set_ylabel("Truth purity of selected sample")
            axes[1, col].set_ylabel("Truth-photon recovery")
            axes[0, col].legend(frameon=False, fontsize=9.0, loc="lower left")
        style_ax(axes[0, col], "y")
        style_ax(axes[1, col], "y")

    fig.text(0.055, 0.965, "Truth performance: box cuts vs base v3E + centrality BDT", fontsize=17.2, fontweight="bold")
    fig.text(
        0.055,
        0.918,
        "Existing selections shown directly; no efficiency matching, threshold rematching, or score retuning is applied.",
        fontsize=10.8,
        color="#4B5563",
    )
    sphinx_label(fig, x=0.865, y=0.958, size=14.5)
    fig.tight_layout(rect=(0.02, 0.02, 0.99, 0.89))
    fig.savefig(PNG_OUT)
    plt.close(fig)


def main() -> int:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    df = build_table()
    df.to_csv(CSV_OUT, index=False)
    plot(df)
    print(PNG_OUT)
    print(CSV_OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
