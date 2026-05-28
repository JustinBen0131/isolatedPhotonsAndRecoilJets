#!/usr/bin/env python3
from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch, Rectangle


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
RUN_ROOT = (
    REPO
    / "dataOutput/auauMLDiagnosticRuns/global_etcent_inclusive3_sixpack_20260516_135439"
)
MODEL_DIR = RUN_ROOT / "validation/basev3e_w33_e22ratio_20260518_192630/bdt_models"
AUC_CSV = (
    RUN_ROOT
    / "slideReady/basev3e_w33_e22ratio_diagnostics/basev3e_w33_e22ratio_stepwise_auc_summary.csv"
)
OUTDIR = REPO / "dataOutput/slide_assets/WP_GammaJets_5_20_26/slide20_e22_split_gain"
OUT_PNG = OUTDIR / "slide20_e22_split_gain_ablation_candidate.png"
OUT_CSV = OUTDIR / "slide20_e22_split_gain_ablation_summary.csv"


PRODUCTS = {
    "E22/E37": "baseBDT_v3E_withCentrality_w33_E22E37",
    "E22/E53": "baseBDT_v3E_withCentrality_w33_E22E53",
    "both ratios": "baseBDT_v3E_withCentrality_w33_E22E37_E22E53",
}

RATIO_FEATURES = {
    "e22_over_e37": r"$E_{22}/E_{37}$",
    "e22_over_e53": r"$E_{22}/E_{53}$",
}


def read_json_with_preamble(path: Path) -> dict:
    text = path.read_text()
    return json.loads(text[text.find("{") :])


def split_gain(product: str) -> pd.DataFrame:
    meta = json.loads((MODEL_DIR / f"auau_tight_bdt_{product}_tmva.metadata.json").read_text())
    features = list(meta["features"])
    model = read_json_with_preamble(MODEL_DIR / f"auau_tight_bdt_{product}_tmva.xgb.json")
    gain = defaultdict(float)
    count = defaultdict(int)

    for tree in model["learner"]["gradient_booster"]["model"]["trees"]:
        for left, right, idx, loss in zip(
            tree["left_children"],
            tree["right_children"],
            tree["split_indices"],
            tree["loss_changes"],
        ):
            if int(left) < 0 and int(right) < 0:
                continue
            idx = int(idx)
            if 0 <= idx < len(features):
                feature = features[idx]
                gain[feature] += max(float(loss), 0.0)
                count[feature] += 1

    total = sum(gain.values())
    rows = []
    for feature in features:
        rows.append(
            {
                "feature": feature,
                "split_gain": gain[feature],
                "split_count": count[feature],
                "split_gain_fraction": gain[feature] / total if total else 0.0,
            }
        )
    return pd.DataFrame(rows).sort_values("split_gain_fraction", ascending=False)


def add_box(
    fig: plt.Figure,
    xy: tuple[float, float],
    wh: tuple[float, float],
    face: str,
    *,
    edge: str = "none",
    lw: float = 0.0,
    radius: float = 0.018,
    zorder: int = -1,
) -> None:
    fig.patches.append(
        FancyBboxPatch(
            xy,
            wh[0],
            wh[1],
            boxstyle=f"round,pad=0.010,rounding_size={radius}",
            transform=fig.transFigure,
            facecolor=face,
            edgecolor=edge,
            linewidth=lw,
            zorder=zorder,
        )
    )


def add_text(
    fig: plt.Figure,
    x: float,
    y: float,
    text: str,
    *,
    size: float,
    weight: str = "normal",
    color: str = "#111827",
    ha: str = "left",
    va: str = "top",
    linespacing: float = 1.08,
) -> None:
    fig.text(
        x,
        y,
        text,
        fontsize=size,
        fontweight=weight,
        color=color,
        ha=ha,
        va=va,
        linespacing=linespacing,
    )


def build_summary() -> pd.DataFrame:
    auc = pd.read_csv(AUC_CSV)
    w33_auc = float(auc.loc[auc["short"].eq("w33"), "auc"].iloc[0])
    rows = []

    for label, product in PRODUCTS.items():
        gain_df = split_gain(product)
        auc_row = auc.loc[auc["short"].eq(label)].iloc[0]
        for feature, pretty in RATIO_FEATURES.items():
            feature_row = gain_df.loc[gain_df["feature"].eq(feature)]
            if feature_row.empty:
                fraction = 0.0
                count = 0
            else:
                fraction = float(feature_row["split_gain_fraction"].iloc[0])
                count = int(feature_row["split_count"].iloc[0])
            rows.append(
                {
                    "model": label,
                    "product": product,
                    "feature": feature,
                    "feature_label": pretty,
                    "auc": float(auc_row["auc"]),
                    "delta_auc_vs_w33": float(auc_row["auc"]) - w33_auc,
                    "split_gain_fraction": fraction,
                    "split_gain_percent": 100.0 * fraction,
                    "split_count": count,
                }
            )

    out = pd.DataFrame(rows)
    OUTDIR.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_CSV, index=False)
    return out


def draw() -> None:
    summary = build_summary()

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "dejavuserif",
            "axes.unicode_minus": False,
        }
    )

    fig = plt.figure(figsize=(16, 9), dpi=180, facecolor="white")
    fig.patches.append(Rectangle((0, 0), 1, 1, transform=fig.transFigure, facecolor="white", zorder=-5))

    navy = "#101828"
    gray = "#475467"
    panel_edge = "#D0D5DD"
    soft_gray = "#F2F4F7"
    soft_yellow = "#FFF2B8"
    soft_blue = "#EAF2FB"
    blue = "#1F77B4"
    orange = "#F26A21"
    green = "#189A5B"

    add_text(
        fig,
        0.035,
        0.955,
        "E22 ratios are heavily used, and both give the best validation AUC",
        size=28.5,
        weight="bold",
        color=navy,
    )
    add_text(
        fig,
        0.035,
        0.907,
        "0-20% centrality, 15 < cluster $E_T$ < 35 GeV; each row adds only the listed E22 ratio(s) to base v3E + centrality + 3x3 widths.",
        size=14.6,
        color=gray,
    )

    add_box(fig, (0.035, 0.255), (0.455, 0.585), "white", edge=panel_edge, lw=1.0, radius=0.012)
    add_box(fig, (0.520, 0.255), (0.445, 0.585), "white", edge=panel_edge, lw=1.0, radius=0.012)

    add_text(fig, 0.058, 0.805, "Validation response", size=20.0, weight="bold", color=navy)
    add_text(fig, 0.058, 0.772, r"AUC gain relative to base v3E + centrality + 3x3 widths.", size=12.9, color=gray)

    auc_ax = fig.add_axes([0.090, 0.410, 0.360, 0.300])
    model_order = ["E22/E37", "E22/E53", "both ratios"]
    auc_data = (
        summary[summary["feature"].eq("e22_over_e37")]
        .drop_duplicates("model")
        .set_index("model")
        .loc[model_order]
        .reset_index()
    )
    x = np.arange(len(auc_data))
    heights = 1000.0 * auc_data["delta_auc_vs_w33"].to_numpy()
    colors = [blue, orange, green]
    auc_ax.bar(x, heights, color=colors, width=0.62, edgecolor="none")
    auc_ax.axhline(0, color="#667085", linewidth=1.0)
    auc_ax.set_xticks(x)
    auc_ax.set_xticklabels([r"+$E_{22}/E_{37}$", r"+$E_{22}/E_{53}$", "both"], fontsize=12.5)
    auc_ax.set_ylabel(r"$\Delta$AUC vs 3x3 baseline [$\times 10^{-3}$]", fontsize=12.3)
    auc_ax.set_ylim(0, max(heights) * 1.28)
    auc_ax.grid(axis="y", color="#E4E7EC", linewidth=0.9)
    auc_ax.set_axisbelow(True)
    auc_ax.tick_params(axis="y", labelsize=11.5, direction="in", right=True)
    auc_ax.tick_params(axis="x", direction="in", top=True)
    for spine in auc_ax.spines.values():
        spine.set_color("#344054")
        spine.set_linewidth(0.9)

    for i, (_, row) in enumerate(auc_data.iterrows()):
        auc = float(row["auc"])
        delta = 1000.0 * float(row["delta_auc_vs_w33"])
        auc_ax.text(
            i,
            heights[i] + max(heights) * 0.035,
            f"+{delta:.1f}\nAUC {auc:.4f}",
            ha="center",
            va="bottom",
            fontsize=11.8,
            fontweight="bold",
            linespacing=0.93,
            color=navy,
        )

    add_box(fig, (0.070, 0.300), (0.385, 0.078), soft_gray, radius=0.016)
    add_text(
        fig,
        0.088,
        0.352,
        "Result",
        size=14.0,
        weight="bold",
        color=navy,
    )
    add_text(
        fig,
        0.088,
        0.326,
        "The paired-ratio model gives the largest AUC increase, not just the largest feature list.",
        size=12.2,
        color="#1F2937",
    )

    add_text(fig, 0.545, 0.805, "Where the BDT spends split gain", size=20.0, weight="bold", color=navy)
    add_text(
        fig,
        0.545,
        0.772,
        "Gain sums loss reduction over all tree splits using each added ratio.",
        size=12.9,
        color=gray,
    )

    gain_ax = fig.add_axes([0.565, 0.395, 0.350, 0.330])
    gain_ax.set_xlim(0, 64)
    gain_ax.set_ylim(-0.65, 2.65)
    gain_ax.set_yticks([2, 1, 0])
    gain_ax.set_yticklabels([r"+$E_{22}/E_{37}$", r"+$E_{22}/E_{53}$", "both"], fontsize=12.7)
    gain_ax.set_xlabel("Added-ratio split gain [%]", fontsize=12.4, labelpad=6)
    gain_ax.grid(axis="x", color="#E4E7EC", linewidth=0.9)
    gain_ax.set_axisbelow(True)
    gain_ax.tick_params(axis="x", labelsize=11.5, direction="in", top=True)
    gain_ax.tick_params(axis="y", direction="in", right=True)
    for spine in gain_ax.spines.values():
        spine.set_color("#344054")
        spine.set_linewidth(0.9)

    def feature_value(model: str, feature: str) -> tuple[float, int]:
        row = summary[(summary["model"].eq(model)) & (summary["feature"].eq(feature))].iloc[0]
        return float(row["split_gain_percent"]), int(row["split_count"])

    e37_alone, e37_alone_count = feature_value("E22/E37", "e22_over_e37")
    e53_alone, e53_alone_count = feature_value("E22/E53", "e22_over_e53")
    e37_both, e37_both_count = feature_value("both ratios", "e22_over_e37")
    e53_both, e53_both_count = feature_value("both ratios", "e22_over_e53")

    gain_ax.barh(2, e37_alone, height=0.48, color=blue, edgecolor="none")
    gain_ax.barh(1, e53_alone, height=0.48, color=orange, edgecolor="none")
    gain_ax.barh(0, e37_both, height=0.48, color=blue, edgecolor="none")
    gain_ax.barh(0, e53_both, height=0.48, left=e37_both, color=orange, edgecolor="none")

    gain_ax.text(e37_alone + 1.2, 2, f"{e37_alone:.1f}%\n({e37_alone_count} splits)", va="center", fontsize=11.5, fontweight="bold", linespacing=0.95, color=navy)
    gain_ax.text(e53_alone + 1.2, 1, f"{e53_alone:.1f}%\n({e53_alone_count} splits)", va="center", fontsize=11.5, fontweight="bold", linespacing=0.95, color=navy)
    gain_ax.text(e37_both / 2, 0, f"{e37_both:.1f}%", ha="center", va="center", fontsize=11.0, fontweight="bold", color="white")
    gain_ax.text(e37_both + e53_both / 2, 0, f"{e53_both:.1f}%", ha="center", va="center", fontsize=11.0, fontweight="bold", color="white")
    total_both = e37_both + e53_both
    total_both_count = e37_both_count + e53_both_count
    gain_ax.text(
        min(total_both, 58.8),
        0.40,
        f"{total_both:.1f}% combined\n{total_both_count} total splits",
        va="bottom",
        ha="center",
        fontsize=10.8,
        fontweight="bold",
        color=navy,
        linespacing=0.92,
    )

    fig.patches.append(Rectangle((0.682, 0.318), 0.012, 0.012, transform=fig.transFigure, color=blue, zorder=2))
    add_text(fig, 0.700, 0.330, r"$E_{22}/E_{37}$ gain", size=11.6, color=gray, va="center")
    fig.patches.append(Rectangle((0.790, 0.318), 0.012, 0.012, transform=fig.transFigure, color=orange, zorder=2))
    add_text(fig, 0.808, 0.330, r"$E_{22}/E_{53}$ gain", size=11.6, color=gray, va="center")

    add_box(fig, (0.035, 0.065), (0.445, 0.145), soft_yellow, radius=0.020)
    add_text(fig, 0.055, 0.175, "Interpretation", size=18.6, weight="bold", color=navy)
    add_text(
        fig,
        0.055,
        0.139,
        r"$E_{22}/E_{37}$ dominates when added alone; $E_{22}/E_{53}$ still adds"
        "\ncomplementary ranking power, so the two-ratio model has the best AUC.",
        size=12.9,
        color="#1F2937",
        linespacing=1.04,
    )

    add_box(fig, (0.520, 0.060), (0.445, 0.150), soft_blue, radius=0.020)
    add_text(fig, 0.542, 0.177, "PPG12 Strategy", size=19.4, weight="bold", color=navy)
    add_text(
        fig,
        0.542,
        0.143,
        "PPG12 used split-count importance:\n"
        "how often each input appears in tree cuts.\n"
        "Gain is loss-weighted; the ranking should\n"
        "be similar if the BDT is stable.",
        size=11.8,
        color="#1F2937",
        linespacing=1.04,
    )
    add_text(
        fig,
        0.752,
        0.143,
        "Next cross-check: permutation importance.\n"
        "Shuffle one input at a time on validation\n"
        "and measure the AUC drop.",
        size=11.8,
        color="#1F2937",
        linespacing=1.04,
    )

    fig.savefig(OUT_PNG, dpi=180, facecolor="white")
    plt.close(fig)
    print(OUT_PNG)
    print(OUT_CSV)


if __name__ == "__main__":
    draw()
