#!/usr/bin/env python3
"""Make a 0-20% centrality AUC leaderboard including E22 ablations and the stack."""

from __future__ import annotations

from pathlib import Path
from textwrap import fill

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import pandas as pd
import json


OUT_DIR = Path(
    "dataOutput/auauMLDiagnosticRuns/"
    "global_etcent_inclusive3_sixpack_20260516_135439/"
    "slideReady/bdt_centrality_performance"
)
VALIDATION_DIR = OUT_DIR.parents[1] / "validation"
STACK_SCORE_CSV = (
    OUT_DIR.parent
    / "stack_comparison"
    / "score_separation_3x3_bdt_mlp_stack_pt15to35_test.csv"
)
CENTRALITY_KEY = "0_20"
CENTRALITY_LABEL = "0-20%"


COLORS = {
    "BDT default/comparison": "#2F78A8",
    "BDT+MLP stack (held-out test)": "#0F766E",
    "NN comparison": "#2B946D",
    "diagnostic uses isolation": "#8C5AA1",
    "BDT feature-control": "#7EACC7",
    "E22-ratio ablation": "#D97706",
}


def add_sphenix_label(fig: plt.Figure) -> None:
    fig.text(
        0.055,
        0.865,
        "sPHENIX",
        ha="left",
        va="center",
        fontsize=18,
        fontstyle="italic",
        fontweight="bold",
    )
    fig.text(0.132, 0.865, " Internal", ha="left", va="center", fontsize=18)


def wrap_label(text: str) -> str:
    text = text.replace(" + ", " + ")
    return fill(text, width=42, break_long_words=False)


def load_json(path: Path) -> dict:
    return json.loads(path.read_text())


def auc_by_cent(path: Path, product: str, cent_key: str = CENTRALITY_KEY) -> float:
    data = load_json(path)
    record = data["products"][product]["auc_by_centrality"][cent_key]
    if isinstance(record, dict):
        return float(record["auc"])
    return float(record)


def auc_from_summary(path: Path, product: str, cent_key: str = CENTRALITY_KEY) -> float:
    data = load_json(path)
    for record in data["summary"]:
        if record["model"] == product and record["cent_bin"] == cent_key:
            return float(record["auc"])
    raise KeyError(f"{product} / {cent_key} not found in {path}")


def read_csv_after_header(path: Path, header_prefix: str) -> pd.DataFrame:
    lines = path.read_text().splitlines()
    for idx, line in enumerate(lines):
        if line.startswith(header_prefix):
            return pd.read_csv(path, skiprows=idx)
    raise RuntimeError(f"No CSV header starting with {header_prefix!r} in {path}")


def weighted_auc_from_pt_table(path: Path, product: str, cent_key: str = CENTRALITY_KEY) -> float:
    table = read_csv_after_header(path, "product,")
    table = table[(table["product"] == product) & (table["centrality_bin"] == cent_key)]
    table = table[table["entries"] > 0].copy()
    return float((table["auc"] * table["entries"]).sum() / table["entries"].sum())


def current_stack_auc() -> float:
    stack_df = pd.read_csv(STACK_SCORE_CSV)
    stack_df = stack_df[stack_df["model"] == "BDT+MLP stack, NN combiner"].copy()
    if stack_df.empty:
        raise RuntimeError(f"No BDT+MLP stack rows found in {STACK_SCORE_CSV}")
    stack_df = stack_df[(stack_df["cent_lo"] == 0.0) & (stack_df["cent_hi"] == 20.0)]
    return float(stack_df.iloc[0]["auc"])


def current_small_mlp_auc() -> float:
    stack_df = pd.read_csv(STACK_SCORE_CSV)
    stack_df = stack_df[stack_df["model"] == "Input MLP, single small NN"].copy()
    stack_df = stack_df[(stack_df["cent_lo"] == 0.0) & (stack_df["cent_hi"] == 20.0)]
    return float(stack_df.iloc[0]["auc"])


def build_rows() -> list[tuple[str, float, str, str]]:
    bdt_metrics = VALIDATION_DIR / "bdt_finished_only_20260516_180817" / "validation_metrics.json"
    binned_metrics = (
        VALIDATION_DIR
        / "bdt_binned_sidecars_fullstat_20260517_2152"
        / "validation_metrics.json"
    )
    basev3e_metrics = (
        VALIDATION_DIR / "basev3e_controls_20260518_1110" / "validation_deep_diagnostics.json"
    )
    e22_metrics = (
        VALIDATION_DIR
        / "basev3e_w33_e22ratio_20260518_192630"
        / "validation"
        / "validation_metrics.json"
    )
    mlp_large_summary = (
        VALIDATION_DIR
        / "mlp_noIso_large_official_fullstat_20260517"
        / "coarse_centrality_score_histograms_globalEtCent1535_mlp_noIso_summary.json"
    )
    width_table = OUT_DIR / "width_ablation_validation_auc_table.csv"

    return [
        (
            "Binned diagnostic: 8 E_T x 7 cent + isolation inputs",
            auc_by_cent(binned_metrics, "globalEtCent1535_bdt_iso_ptCent7"),
            "diagnostic uses isolation",
            "centrality-resolved validation metrics",
        ),
        (
            "Binned diagnostic: 8 E_T x 3 cent + isolation inputs",
            auc_by_cent(binned_metrics, "globalEtCent1535_bdt_iso_ptCent3"),
            "diagnostic uses isolation",
            "centrality-resolved validation metrics",
        ),
        (
            "Global BDT diagnostic: isolation inputs",
            auc_by_cent(bdt_metrics, "globalEtCent1535_bdt_iso"),
            "diagnostic uses isolation",
            "centrality-resolved validation metrics",
        ),
        (
            "BDT+MLP stack: 8 E_T x 7 cent BDT + small MLP scores",
            current_stack_auc(),
            "BDT+MLP stack (held-out test)",
            "held-out score comparison",
        ),
        (
            "Binned BDT: 8 E_T x 7 cent",
            auc_by_cent(binned_metrics, "globalEtCent1535_bdt_noIso_ptCent7"),
            "BDT default/comparison",
            "centrality-resolved validation metrics",
        ),
        (
            "Small MLP replay: 32 baseline inputs",
            current_small_mlp_auc(),
            "NN comparison",
            "held-out score comparison",
        ),
        (
            "Large MLP: 32 baseline inputs",
            auc_from_summary(mlp_large_summary, "globalEtCent1535_mlp_noIso"),
            "NN comparison",
            "coarse centrality score-hist summary",
        ),
        (
            "Binned BDT: 8 E_T x 3 cent",
            auc_by_cent(binned_metrics, "globalEtCent1535_bdt_noIso_ptCent3"),
            "BDT default/comparison",
            "centrality-resolved validation metrics",
        ),
        (
            "Global BDT: 32-input AuAu baseline",
            auc_by_cent(bdt_metrics, "globalEtCent1535_bdt_noIso"),
            "BDT default/comparison",
            "centrality-resolved validation metrics",
        ),
        (
            "Global BDT: PPG12 25 + centrality",
            weighted_auc_from_pt_table(width_table, "globalEtCent1535_bdt_ppg12PlusCent"),
            "BDT feature-control",
            "entries-weighted over 15-20, 20-25, 25-35 GeV pT bins",
        ),
        (
            "Global BDT: remove long widths/ratios",
            weighted_auc_from_pt_table(width_table, "globalEtCent1535_bdt_dropElongatedRatios"),
            "BDT feature-control",
            "entries-weighted over 15-20, 20-25, 25-35 GeV pT bins",
        ),
        (
            "E22 ablation: base v3E + 3x3 + cent + both ratios",
            auc_by_cent(e22_metrics, "baseBDT_v3E_withCentrality_w33_E22E37_E22E53"),
            "E22-ratio ablation",
            "centrality-resolved validation metrics",
        ),
        (
            "E22 ablation: base v3E + 3x3 + cent + E22/E37",
            auc_by_cent(e22_metrics, "baseBDT_v3E_withCentrality_w33_E22E37"),
            "E22-ratio ablation",
            "centrality-resolved validation metrics",
        ),
        (
            "E22 ablation: base v3E + 3x3 + cent + E22/E53",
            auc_by_cent(e22_metrics, "baseBDT_v3E_withCentrality_w33_E22E53"),
            "E22-ratio ablation",
            "centrality-resolved validation metrics",
        ),
        (
            "Global BDT: base v3E + 3x3 widths + centrality",
            auc_by_cent(basev3e_metrics, "baseBDT_v3E_withCentrality_w33"),
            "BDT feature-control",
            "centrality-resolved validation diagnostics",
        ),
        (
            "Global BDT: base v3E + 3x3 widths",
            auc_by_cent(basev3e_metrics, "baseBDT_v3E_withOutCentraltiy_w33"),
            "BDT feature-control",
            "centrality-resolved validation diagnostics",
        ),
        (
            "Global BDT: base v3E 11 inputs",
            auc_by_cent(basev3e_metrics, "baseBDT_v3E_withOutCentraltiy"),
            "BDT feature-control",
            "centrality-resolved validation diagnostics",
        ),
        (
            "Global BDT: base v3E + centrality",
            auc_by_cent(basev3e_metrics, "baseBDT_v3E_withCentrality"),
            "BDT feature-control",
            "centrality-resolved validation diagnostics",
        ),
    ]


def draw_leaderboard(df: pd.DataFrame, png_path: Path, title: str, csv_path: Path) -> None:
    plot_df = df.sort_values("auc", ascending=False).reset_index(drop=True)
    plot_df = plot_df.drop(columns=["rank"], errors="ignore")
    plot_df.insert(0, "rank", range(1, len(plot_df) + 1))
    plot_df.to_csv(csv_path, index=False)

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.labelsize": 14,
            "xtick.labelsize": 12,
            "ytick.labelsize": 10.0,
        }
    )

    fig, ax = plt.subplots(figsize=(16, 9), dpi=180)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    y = list(range(len(df)))
    bars = ax.barh(
        y,
        plot_df["auc"],
        height=0.66,
        color=[COLORS[c] for c in plot_df["category"]],
        edgecolor="none",
    )
    ax.invert_yaxis()

    ax.set_yticks(y)
    ax.set_yticklabels([wrap_label(x) for x in plot_df["model_label"]])
    ax.tick_params(axis="y", length=0, pad=12)
    xmin = max(0.0, float(plot_df["auc"].min()) - 0.011)
    xmax = min(1.0, float(plot_df["auc"].max()) + 0.009)
    ax.set_xlim(xmin, xmax)
    ax.set_xticks([0.79, 0.81, 0.83, 0.85, 0.87, 0.89])
    ax.set_xlabel(f"{CENTRALITY_LABEL} centrality validation AUC")
    ax.grid(axis="x", color="#D1D5DB", linewidth=1.0, alpha=0.75)
    ax.set_axisbelow(True)

    for spine in ["top", "right", "left"]:
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_color("#111827")

    for bar, auc in zip(bars, plot_df["auc"]):
        ax.text(
            auc + 0.0032,
            bar.get_y() + bar.get_height() / 2,
            f"{auc:.3f}",
            va="center",
            ha="left",
            fontsize=12.7,
            fontweight="bold",
            color="#111827",
        )

    fig.text(
        0.055,
        0.94,
        title,
        ha="left",
        va="top",
        fontsize=25,
        fontweight="bold",
    )
    fig.text(
        0.055,
        0.905,
        "Signal = embedded Photon12+20; background = embedded inclusive Jet12+20+30; 15 < cluster E_T < 35 GeV.",
        ha="left",
        va="top",
        fontsize=14.5,
        color="#374151",
    )
    add_sphenix_label(fig)

    legend_order = [
        "BDT default/comparison",
        "BDT+MLP stack (held-out test)",
        "NN comparison",
        "diagnostic uses isolation",
        "BDT feature-control",
        "E22-ratio ablation",
    ]
    present_categories = set(plot_df["category"])
    handles = [mpatches.Patch(color=COLORS[k], label=k) for k in legend_order if k in present_categories]
    fig.legend(
        handles=handles,
        loc="lower left",
        bbox_to_anchor=(0.055, 0.028),
        ncol=6,
        frameon=False,
        fontsize=11,
        handlelength=1.0,
        columnspacing=1.55,
    )

    fig.subplots_adjust(left=0.365, right=0.88, top=0.80, bottom=0.13)
    fig.savefig(png_path)
    plt.close(fig)
    print(png_path)
    print(csv_path)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(
        build_rows(),
        columns=["model_label", "auc", "category", "source_note"],
    )
    df.insert(0, "centrality_bin", CENTRALITY_LABEL)

    draw_leaderboard(
        df,
        OUT_DIR / "all_validated_jet1230_outputs_auc_leaderboard_with_e22.png",
        f"Photon12+20 vs Jet12+20+30 AUC leaderboard ({CENTRALITY_LABEL} centrality)",
        OUT_DIR / "all_validated_jet1230_outputs_auc_leaderboard_with_e22.csv",
    )

    no_iso_df = df[df["category"] != "diagnostic uses isolation"].copy()
    draw_leaderboard(
        no_iso_df,
        OUT_DIR / "all_validated_jet1230_outputs_auc_leaderboard_with_e22_no_iso_inputs.png",
        f"AUC leaderboard without isolation-input models ({CENTRALITY_LABEL})",
        OUT_DIR / "all_validated_jet1230_outputs_auc_leaderboard_with_e22_no_iso_inputs.csv",
    )


if __name__ == "__main__":
    main()
