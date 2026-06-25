#!/usr/bin/env python3
"""Build a full-slide PNG for the first AuAu BDT purity-closure probe."""

from __future__ import annotations

import json
import math
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-codex")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
import numpy as np
import pandas as pd


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
INPUT = REPO / (
    "dataOutput/auauMLDiagnosticRuns/global_etcent_inclusive3_sixpack_20260516_135439/"
    "slideReady/current_best_target80_recoiljets/"
    "reference_vs_global32_vs_ptcent7bdt_target80_reco_abcd_raw_purity_pt15to35_slide10_style_points.csv"
)
OUTDIR = REPO / (
    "dataOutput/auauTightBDTValidation/THE59_model_comparison_remakes_20260620/"
    "purity_closure_20260622"
)
PNG = OUTDIR / "the8_purity_closure_target80_sim_v1.png"
CSV = OUTDIR / "the8_purity_closure_target80_sim_v1_datapoints.csv"
LAYOUT = OUTDIR / "the8_purity_closure_target80_sim_v1.layout_nodes.json"
MANIFEST = OUTDIR / "the8_purity_closure_target80_sim_v1.manifest.json"
SCRIPT = OUTDIR / "the8_purity_closure_target80_sim_v1_speaker_script.md"

W, H = 2560, 1440
DPI = 160
INK = "#111827"
MUTED = "#4B5563"
GRID = "#DCE3EC"
EDGE = "#CBD5E1"
RED = "#DC2626"
BG = "#FFFFFF"
PANEL_BG = "#F8FAFC"
MODEL_ORDER = ["Reference box cuts", "Global BDT target-80", "8x7 routed BDT target-80"]
MODEL_DISPLAY = {
    "Reference box cuts": "Reference\nbox cuts",
    "Global BDT target-80": "Global BDT\n target-80",
    "8x7 routed BDT target-80": "8x7 routed BDT\n target-80",
}
MODEL_COLORS = {
    "Reference box cuts": "#6B7280",
    "Global BDT target-80": "#0284C7",
    "8x7 routed BDT target-80": "#15803D",
}
CENT_ORDER = ["0_20", "20_50", "50_80"]
CENT_LABELS = {"0_20": "0-20%", "20_50": "20-50%", "50_80": "50-80%"}


def add_node(nodes, *, name, kind, bbox, text=None, role="audience", **extra):
    node = {"name": name, "kind": kind, "bbox": [float(x) for x in bbox], "role": role}
    if text is not None:
        node["text"] = text
    node.update(extra)
    nodes.append(node)


def raw_abcd_from_counts(rows: pd.DataFrame) -> float:
    a = float(rows["A"].sum())
    b = float(rows["B"].sum())
    c = float(rows["C"].sum())
    d = float(rows["D"].sum())
    if a <= 0 or d <= 0:
        return math.nan
    return max(0.0, 1.0 - (b * c / d) / a)


def truth_purity_from_counts(rows: pd.DataFrame) -> float:
    sig = float(rows["signal_tight"].sum())
    bkg = float(rows["background_tight"].sum())
    if sig + bkg <= 0:
        return math.nan
    return sig / (sig + bkg)


def build_datapoints() -> pd.DataFrame:
    df = pd.read_csv(INPUT)
    rows = []
    for model in MODEL_ORDER:
        for cent in CENT_ORDER:
            g = df[(df["model"] == model) & (df["centrality"] == cent)]
            reco = raw_abcd_from_counts(g[g["quantity"] == "reco_abcd_raw_purity"])
            truth = truth_purity_from_counts(g[g["quantity"] == "truth_labeled_A_purity"])
            rows.append(
                {
                    "model": model,
                    "model_display": MODEL_DISPLAY[model].replace("\n", " "),
                    "centrality": cent,
                    "centrality_label": CENT_LABELS[cent],
                    "truth_region_a_purity": truth,
                    "reco_abcd_raw_purity": reco,
                    "closure_residual_abcd_minus_truth": reco - truth,
                    "truth_region_a_purity_percent": 100.0 * truth,
                    "reco_abcd_raw_purity_percent": 100.0 * reco,
                    "closure_residual_abcd_minus_truth_pp": 100.0 * (reco - truth),
                    "definition": "integrated over 15<E_T<35 GeV; ABCD counts integrated before max(0,A-BC/D)/A; truth counts integrated before S_A/(S_A+B_A)",
                }
            )
    return pd.DataFrame(rows)


def style_axis(ax):
    ax.set_facecolor(PANEL_BG)
    ax.grid(axis="y", color=GRID, linewidth=1.5)
    ax.set_axisbelow(True)
    ax.tick_params(axis="both", labelsize=16, colors=INK, length=0)
    for spine in ax.spines.values():
        spine.set_color(EDGE)
        spine.set_linewidth(2.0)


def grouped_bars(ax, data: pd.DataFrame, column: str, *, ylim, ylabel: str, value_format: str, zero_line=False):
    x = np.arange(len(CENT_ORDER))
    width = 0.23
    offsets = np.linspace(-width, width, len(MODEL_ORDER))
    for offset, model in zip(offsets, MODEL_ORDER, strict=True):
        vals = [
            float(data[(data["model"] == model) & (data["centrality"] == cent)][column].iloc[0])
            for cent in CENT_ORDER
        ]
        bars = ax.bar(
            x + offset,
            vals,
            width=width * 0.92,
            color=MODEL_COLORS[model],
            edgecolor="white",
            linewidth=1.6,
            label=MODEL_DISPLAY[model].replace("\n", " "),
            zorder=3,
        )
        for bar, val in zip(bars, vals, strict=True):
            if column == "closure_residual_abcd_minus_truth_pp":
                y = val - 2.2
                va = "top"
            else:
                y = val + 1.2
                va = "bottom"
            ax.text(
                bar.get_x() + bar.get_width() / 2.0,
                y,
                value_format.format(val),
                ha="center",
                va=va,
                fontsize=15.0,
                fontweight="bold",
                color=MODEL_COLORS[model],
                clip_on=False,
            )
    ax.set_xticks(x, [CENT_LABELS[c] for c in CENT_ORDER], fontsize=18, fontweight="bold")
    ax.set_ylim(*ylim)
    ax.set_ylabel(ylabel, fontsize=16.5, fontweight="bold", color=INK, labelpad=12)
    if zero_line:
        ax.axhline(0.0, color="#0F172A", linewidth=2.0, alpha=0.75, zorder=2)
    style_axis(ax)


def draw_header(fig, nodes):
    fig.text(
        0.045,
        0.952,
        "BDT changes the purity estimator more than the truth purity",
        ha="left",
        va="top",
        fontsize=39,
        fontweight="bold",
        color=INK,
        family="Times New Roman",
    )
    add_node(
        nodes,
        name="title",
        kind="text",
        role="title",
        bbox=[115, 58, 1840, 115],
        text="BDT changes the purity estimator more than the truth purity",
        font_px=38 * DPI / 72.0,
        title_anchor=True,
    )

    ax = fig.add_axes([0, 0, 1, 1], zorder=-1)
    ax.axis("off")
    box1 = (0.049, 0.812, 0.515, 0.078)
    box2 = (0.595, 0.812, 0.360, 0.078)
    for x, y, w, h, color in (*[(*box1, RED)], *[(*box2, "#94A3B8")]):
        ax.add_patch(
            FancyBboxPatch(
                (x, y),
                w,
                h,
                boxstyle="round,pad=0.007,rounding_size=0.008",
                facecolor="white",
                edgecolor=color,
                linewidth=2.7,
                transform=fig.transFigure,
            )
        )
    fig.text(
        box1[0] + 0.018,
        box1[1] + box1[3] * 0.62,
        "Existing RecoilJets target-80 SIM products",
        ha="left",
        va="center",
        fontsize=18.5,
        fontweight="bold",
        color=RED,
        family="Times New Roman",
    )
    fig.text(
        box1[0] + 0.018,
        box1[1] + box1[3] * 0.29,
        r"Photon12+20 signal, Jet12+20+30 inclusive, integrated over 15 < $E_T$ < 35 GeV",
        ha="left",
        va="center",
        fontsize=17.5,
        fontweight="bold",
        color="#991B1B",
        family="Times New Roman",
    )
    add_node(
        nodes,
        name="input_contract_box",
        kind="box",
        bbox=[125, 158, 1445, 270],
        role="audience",
    )
    add_node(
        nodes,
        name="input_contract_text",
        kind="text",
        bbox=[170, 185, 1425, 252],
        text="Existing RecoilJets target-80 SIM products; Photon12+20 signal, Jet12+20+30 inclusive, integrated over 15 < E_T < 35 GeV",
        font_px=18.5 * DPI / 72.0,
    )

    fig.text(
        box2[0] + 0.018,
        box2[1] + box2[3] * 0.64,
        r"ABCD = $\max(0, A-BC/D)/A$",
        ha="left",
        va="center",
        fontsize=17.5,
        fontweight="bold",
        color=INK,
        family="Times New Roman",
    )
    fig.text(
        box2[0] + 0.018,
        box2[1] + box2[3] * 0.34,
        r"Truth = $S_A/(S_A+B_A)$; residual = ABCD - truth",
        ha="left",
        va="center",
        fontsize=17.5,
        fontweight="bold",
        color=INK,
        family="Times New Roman",
    )
    add_node(nodes, name="definition_box", kind="box", bbox=[1523, 158, 2445, 270], role="audience")
    add_node(
        nodes,
        name="definition_text",
        kind="text",
        bbox=[1568, 185, 2415, 252],
        text="ABCD = max(0, A-BC/D)/A; Truth = S_A/(S_A+B_A); residual = ABCD - truth",
        font_px=17.5 * DPI / 72.0,
    )


def draw_legend(fig, nodes):
    ax = fig.add_axes([0.090, 0.735, 0.820, 0.050])
    ax.axis("off")
    x = 0.02
    for model in MODEL_ORDER:
        ax.add_patch(Rectangle((x, 0.30), 0.032, 0.32, color=MODEL_COLORS[model], transform=ax.transAxes))
        ax.text(
            x + 0.044,
            0.46,
            MODEL_DISPLAY[model].replace("\n", " "),
            ha="left",
            va="center",
            fontsize=18.0,
            fontweight="bold",
            color=INK,
            transform=ax.transAxes,
            family="Times New Roman",
        )
        x += 0.305
    fig.text(
        0.5,
        0.711,
        "Closure residual uses percentage points; less negative means the reco ABCD estimate is closer to the truth-labeled purity.",
        ha="center",
        va="center",
        fontsize=17.0,
        color=MUTED,
        fontweight="bold",
        family="Times New Roman",
    )
    add_node(
        nodes,
        name="legend_row",
        kind="text",
        bbox=[230, 1000, 2330, 1090],
        text="Reference box cuts, Global BDT target-80, 8x7 routed BDT target-80; less negative residual = better closure",
        font_px=18 * DPI / 72.0,
    )


def write_manifest(data: pd.DataFrame) -> None:
    manifest = {
        "schema": "THE8_PURITY_CLOSURE_TARGET80_SIM_SLIDE_V1",
        "png": str(PNG),
        "csv": str(CSV),
        "input_csv": str(INPUT),
        "selection": "existing RecoilJets SIM target-80 products; 15<E_T<35 GeV; Photon12+20 signal and Jet12+20+30 inclusive background",
        "models": MODEL_ORDER,
        "metrics": {
            "truth_region_a_purity": "S_A/(S_A+B_A), after integrating truth-selected region-A counts in 15<E_T<35",
            "reco_abcd_raw_purity": "max(0,A-BC/D)/A, after integrating reco ABCD counts in 15<E_T<35",
            "closure_residual": "reco_abcd_raw_purity - truth_region_a_purity, percentage points on slide",
        },
        "datapoint_summary": data.to_dict(orient="records"),
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2), encoding="utf-8")


def write_script(data: pd.DataFrame) -> None:
    def get(model, cent, col):
        return float(data[(data["model"] == model) & (data["centrality"] == cent)][col].iloc[0])

    lines = [
        "# Speaker script: purity closure first pass",
        "",
        "This slide is the first physics-level check beyond AUC and score separation. It uses the existing RecoilJets target-80 simulation products, not a new production run.",
        "",
        "The top row is truth-labeled region-A purity. Across the three centrality regions, all selections are already near ninety percent truth purity, so the BDT is not radically changing the true selected-sample composition.",
        "",
        "The middle row is the raw reconstructed ABCD estimate, max of zero and A minus BC over D divided by A. This is where the BDT changes the picture: the reference box cuts give a much lower ABCD purity estimate, especially in central AuAu, while the BDT selections lift the estimate substantially.",
        "",
        "The bottom row is the closure residual, ABCD estimate minus truth purity. Negative means the ABCD estimate is conservative relative to truth. The BDT does not make the residual vanish, but it moves the estimate much closer to truth in 0-20 and 20-50 percent.",
        "",
        "Numerically, in 0-20 percent centrality the closure residual improves from about "
        f"{get('Reference box cuts','0_20','closure_residual_abcd_minus_truth_pp'):.0f} percentage points for the reference cuts to "
        f"{get('Global BDT target-80','0_20','closure_residual_abcd_minus_truth_pp'):.0f} and "
        f"{get('8x7 routed BDT target-80','0_20','closure_residual_abcd_minus_truth_pp'):.0f} percentage points for the BDT choices.",
        "",
        "The interpretation I would carry forward is: score separation is not enough; this first purity-probe says the BDT makes the raw ABCD estimate more compatible with truth, but we still need the routed corrected THE8 products before claiming the full current model family is purity-closed.",
    ]
    SCRIPT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    data = build_datapoints()
    data.to_csv(CSV, index=False)

    plt.rcParams.update(
        {
            "font.family": "Times New Roman",
            "mathtext.fontset": "dejavuserif",
            "figure.facecolor": BG,
            "axes.titleweight": "bold",
            "axes.labelcolor": INK,
            "xtick.color": INK,
            "ytick.color": INK,
        }
    )
    fig = plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI)
    nodes = []
    draw_header(fig, nodes)
    draw_legend(fig, nodes)

    axes = [
        fig.add_axes([0.120, 0.545, 0.815, 0.135]),
        fig.add_axes([0.120, 0.335, 0.815, 0.135]),
        fig.add_axes([0.120, 0.125, 0.815, 0.135]),
    ]
    grouped_bars(
        axes[0],
        data,
        "truth_region_a_purity_percent",
        ylim=(83.5, 94.0),
        ylabel="Truth purity (%)",
        value_format="{:.1f}",
    )
    axes[0].set_xticklabels([])
    fig.text(
        0.120,
        0.690,
        "Truth-labeled selected sample stays near 90%",
        fontsize=24,
        fontweight="bold",
        color=INK,
        ha="left",
        va="bottom",
        family="Times New Roman",
    )
    grouped_bars(
        axes[1],
        data,
        "reco_abcd_raw_purity_percent",
        ylim=(18.0, 73.0),
        ylabel="ABCD purity (%)",
        value_format="{:.0f}",
    )
    axes[1].set_xticklabels([])
    fig.text(
        0.120,
        0.480,
        "Raw ABCD estimate rises strongly under target-80 BDT selections",
        fontsize=24,
        fontweight="bold",
        color=INK,
        ha="left",
        va="bottom",
        family="Times New Roman",
    )
    grouped_bars(
        axes[2],
        data,
        "closure_residual_abcd_minus_truth_pp",
        ylim=(-70.0, 8.0),
        ylabel="ABCD - truth (pp)",
        value_format="{:.0f}",
        zero_line=True,
    )
    fig.text(
        0.120,
        0.270,
        "Closure residual remains conservative, but is less severe for the BDTs",
        fontsize=24,
        fontweight="bold",
        color=INK,
        ha="left",
        va="bottom",
        family="Times New Roman",
    )

    add_node(nodes, name="truth_purity_panel", kind="plot", bbox=[307, 461, 2394, 655], role="audience")
    add_node(nodes, name="abcd_purity_panel", kind="plot", bbox=[307, 763, 2394, 957], role="audience")
    add_node(nodes, name="closure_panel", kind="plot", bbox=[307, 1066, 2394, 1260], role="audience")

    fig.savefig(PNG)
    plt.close(fig)
    write_manifest(data)
    write_script(data)
    LAYOUT.write_text(json.dumps({"nodes": nodes, "canvas": {"width": W, "height": H}}, indent=2), encoding="utf-8")
    print(PNG)
    print(CSV)
    print(SCRIPT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
