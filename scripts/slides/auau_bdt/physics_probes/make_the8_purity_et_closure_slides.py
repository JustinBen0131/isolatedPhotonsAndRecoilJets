#!/usr/bin/env python3
"""Make ET-dependent purity/closure full-slide candidates from existing SIM outputs."""

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
W, H = 2560, 1440
DPI = 160
INK = "#111827"
MUTED = "#475569"
EDGE = "#CBD5E1"
GRID = "#DCE5EE"
PANEL_BG = "#F8FAFC"
RED = "#DC2626"

REFERENCE = "Reference box cuts"
BDT_MODELS = ["Global BDT target-80", "8x7 routed BDT target-80"]
MODEL_TITLES = {
    "Global BDT target-80": "Global BDT target-80",
    "8x7 routed BDT target-80": "8×7 routed BDT target-80",
}
CENT_ORDER = ["0_20", "20_50", "50_80"]
CENT_LABELS = {"0_20": "0-20%", "20_50": "20-50%", "50_80": "50-80%"}
CENT_COLORS = {"0_20": "#B91C1C", "20_50": "#2563EB", "50_80": "#15803D"}
CENT_MARKERS = {"0_20": "o", "20_50": "s", "50_80": "^"}


def add_node(nodes, *, name, kind, bbox, text=None, role="audience", font_px=None, **extra):
    node = {"name": name, "kind": kind, "bbox": [float(x) for x in bbox], "role": role}
    if text is not None:
        node["text"] = text
    if font_px is not None:
        node["font_px"] = float(font_px)
    node.update(extra)
    nodes.append(node)


def raw_abcd_from_counts(a: float, b: float, c: float, d: float) -> float:
    if a <= 0 or d <= 0:
        return math.nan
    return max(0.0, 1.0 - (b * c / d) / a)


def truth_purity_from_counts(sig: float, bkg: float) -> float:
    if sig + bkg <= 0:
        return math.nan
    return sig / (sig + bkg)


def overlap_weight(src_low: float, src_high: float, dst_low: float, dst_high: float) -> float:
    overlap = max(0.0, min(src_high, dst_high) - max(src_low, dst_low))
    width = max(0.0, src_high - src_low)
    if width <= 0:
        return 0.0
    return overlap / width


def reference_values_for_bin(ref_raw: pd.DataFrame, cent: str, lo: float, hi: float) -> tuple[float, float]:
    ref_cent = ref_raw[ref_raw["centrality"] == cent]
    reco_rows = ref_cent[ref_cent["quantity"] == "reco_abcd_raw_purity"]
    truth_rows = ref_cent[ref_cent["quantity"] == "truth_labeled_A_purity"]

    abcd_counts = {"A": 0.0, "B": 0.0, "C": 0.0, "D": 0.0}
    for _, src in reco_rows.iterrows():
        weight = overlap_weight(float(src["pt_low"]), float(src["pt_high"]), lo, hi)
        if weight <= 0:
            continue
        for col in abcd_counts:
            abcd_counts[col] += float(src[col]) * weight
    ref_abcd = raw_abcd_from_counts(abcd_counts["A"], abcd_counts["B"], abcd_counts["C"], abcd_counts["D"])

    sig = 0.0
    bkg = 0.0
    for _, src in truth_rows.iterrows():
        weight = overlap_weight(float(src["pt_low"]), float(src["pt_high"]), lo, hi)
        if weight <= 0:
            continue
        sig += float(src["signal_tight"]) * weight
        bkg += float(src["background_tight"]) * weight
    ref_truth = truth_purity_from_counts(sig, bkg)
    return ref_truth, ref_abcd


def load_points() -> pd.DataFrame:
    df = pd.read_csv(INPUT)
    keep = df["model"].isin([REFERENCE, *BDT_MODELS])
    df = df[keep].copy()
    numeric_cols = [
        "pt_low",
        "pt_high",
        "pt_mid",
        "value",
        "A",
        "B",
        "C",
        "D",
        "signal_tight",
        "background_tight",
    ]
    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df["et_label"] = df.apply(lambda r: f"{int(r['pt_low'])}-{int(r['pt_high'])}", axis=1)
    df["et_mid"] = df["pt_mid"].astype(float)
    pivot = df.pivot_table(
        index=["model", "centrality", "pt_low", "pt_high", "pt_mid", "et_label"],
        columns="quantity",
        values="value",
        aggfunc="first",
    ).reset_index()
    pivot.columns.name = None
    ref_raw = df[df["model"] == REFERENCE].copy()
    rows = []
    for _, row in pivot.iterrows():
        lo = float(row["pt_low"])
        hi = float(row["pt_high"])
        truth = float(row["truth_labeled_A_purity"])
        abcd = float(row["reco_abcd_raw_purity"])
        ref_truth, ref_abcd = reference_values_for_bin(ref_raw, row["centrality"], lo, hi)
        residual = abcd - truth
        ref_residual = ref_abcd - ref_truth
        rows.append(
            {
                "model": row["model"],
                "centrality": row["centrality"],
                "centrality_label": CENT_LABELS[row["centrality"]],
                "pt_low": lo,
                "pt_high": hi,
                "pt_mid": float(row["pt_mid"]),
                "et_label": row["et_label"],
                "truth_region_a_purity": truth,
                "reco_abcd_raw_purity": abcd,
                "reference_truth_region_a_purity": ref_truth,
                "reference_reco_abcd_raw_purity": ref_abcd,
                "truth_purity_gain_pp": 100.0 * (truth - ref_truth),
                "abcd_purity_gain_pp": 100.0 * (abcd - ref_abcd),
                "closure_residual_abcd_minus_truth_pp": 100.0 * residual,
                "reference_closure_residual_abcd_minus_truth_pp": 100.0 * ref_residual,
                "closure_abs_residual_improvement_pp": 100.0 * (abs(ref_residual) - abs(residual)),
            }
        )
    return pd.DataFrame(rows).sort_values(["model", "centrality", "pt_mid"]).reset_index(drop=True)


def style_axis(ax, *, ylabel: str | None = None, xlabel: bool = False):
    ax.set_facecolor(PANEL_BG)
    ax.grid(axis="y", color=GRID, linewidth=1.35)
    ax.set_axisbelow(True)
    ax.tick_params(axis="both", labelsize=15, colors=INK, width=1.4, length=4)
    for spine in ax.spines.values():
        spine.set_color(EDGE)
        spine.set_linewidth(2.0)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=17.5, fontweight="bold", color=INK, labelpad=14)
    if xlabel:
        ax.set_xlabel(r"reco photon $E_T$ bin [GeV]", fontsize=17.5, fontweight="bold", color=INK, labelpad=10)


def draw_header(fig, nodes, *, title: str, metric_line: str, metric_detail: str):
    fig.text(
        0.047,
        0.954,
        title,
        ha="left",
        va="top",
        fontsize=37,
        fontweight="bold",
        color=INK,
        family="Times New Roman",
    )
    add_node(nodes, name="title", kind="text", role="title", bbox=[120, 55, 2300, 118], text=title, font_px=37 * DPI / 72.0, title_anchor=True)
    ax = fig.add_axes([0, 0, 1, 1], zorder=-1)
    ax.axis("off")
    box1 = (0.052, 0.824, 0.492, 0.074)
    box2 = (0.585, 0.824, 0.365, 0.074)
    for x, y, w, h, color in [(*box1, RED), (*box2, "#94A3B8")]:
        ax.add_patch(
            FancyBboxPatch(
                (x, y),
                w,
                h,
                boxstyle="round,pad=0.006,rounding_size=0.007",
                facecolor="white",
                edgecolor=color,
                linewidth=2.6,
                transform=fig.transFigure,
            )
        )
    fig.text(
        box1[0] + 0.017,
        box1[1] + box1[3] * 0.62,
        "Existing RecoilJets target-80 SIM products",
        ha="left",
        va="center",
        fontsize=18.0,
        fontweight="bold",
        color=RED,
        family="Times New Roman",
    )
    fig.text(
        box1[0] + 0.017,
        box1[1] + box1[3] * 0.30,
        r"Photon12+20 signal; Jet12+20+30 inclusive MC; source ROOT $E_T$ bins",
        ha="left",
        va="center",
        fontsize=17.0,
        fontweight="bold",
        color="#991B1B",
        family="Times New Roman",
    )
    fig.text(
        box2[0] + 0.017,
        box2[1] + box2[3] * 0.62,
        metric_line,
        ha="left",
        va="center",
        fontsize=17.0,
        fontweight="bold",
        color=INK,
        family="Times New Roman",
    )
    fig.text(
        box2[0] + 0.017,
        box2[1] + box2[3] * 0.30,
        metric_detail,
        ha="left",
        va="center",
        fontsize=16.0,
        fontweight="bold",
        color=MUTED,
        family="Times New Roman",
    )
    add_node(nodes, name="input_contract", kind="text", bbox=[155, 165, 1395, 262], text="Existing RecoilJets target-80 SIM products; Photon12+20 signal, Jet12+20+30 inclusive", font_px=18 * DPI / 72.0)
    add_node(nodes, name="metric_definition", kind="text", bbox=[1500, 165, 2415, 262], text=metric_line, font_px=17.5 * DPI / 72.0)


def draw_legend(fig, nodes):
    ax = fig.add_axes([0.120, 0.762, 0.760, 0.043])
    ax.axis("off")
    x = 0.02
    for cent in CENT_ORDER:
        ax.plot([x, x + 0.040], [0.52, 0.52], color=CENT_COLORS[cent], lw=4.0, transform=ax.transAxes)
        ax.plot([x + 0.020], [0.52], marker=CENT_MARKERS[cent], color=CENT_COLORS[cent], ms=8, transform=ax.transAxes)
        ax.text(
            x + 0.052,
            0.52,
            f"{CENT_LABELS[cent]} centrality",
            ha="left",
            va="center",
            fontsize=17.5,
            fontweight="bold",
            color=INK,
            family="Times New Roman",
            transform=ax.transAxes,
        )
        x += 0.285
    add_node(nodes, name="centrality_legend", kind="text", bbox=[305, 1080, 2245, 1140], text="0-20%, 20-50%, and 50-80% centrality curves", font_px=17.5 * DPI / 72.0)


def plot_column(ax_top, ax_bottom, data: pd.DataFrame, model: str, *, top_col: str, gain_col: str, top_ylabel: str | None, bottom_ylabel: str | None, top_ylim, gain_ylim, gain_label_fmt="{:+.1f}"):
    gmodel = data[data["model"] == model]
    for cent in CENT_ORDER:
        g = gmodel[gmodel["centrality"] == cent].sort_values("pt_mid")
        color = CENT_COLORS[cent]
        marker = CENT_MARKERS[cent]
        ax_top.plot(
            g["pt_mid"],
            100.0 * g[top_col],
            marker=marker,
            ms=7.5,
            lw=3.0,
            color=color,
            label=CENT_LABELS[cent],
        )
        ax_bottom.plot(
            g["pt_mid"],
            g[gain_col],
            marker=marker,
            ms=7.5,
            lw=3.0,
            color=color,
        )
        for _, row in g.iloc[[0, -1]].iterrows():
            ax_bottom.text(
                row["pt_mid"],
                row[gain_col] + (0.45 if row[gain_col] >= 0 else -0.65),
                gain_label_fmt.format(row[gain_col]),
                ha="center",
                va="bottom" if row[gain_col] >= 0 else "top",
                fontsize=12.5,
                fontweight="bold",
                color=color,
                clip_on=False,
            )
    ax_top.set_title(MODEL_TITLES[model], fontsize=22, fontweight="bold", color=INK, pad=8, family="Times New Roman")
    ax_top.set_ylim(*top_ylim)
    ax_bottom.set_ylim(*gain_ylim)
    ax_bottom.axhline(0.0, color="#0F172A", lw=1.8, alpha=0.7)
    xticks = sorted(gmodel["pt_mid"].unique())
    xlabels = [str(gmodel[gmodel["pt_mid"] == x]["et_label"].iloc[0]) for x in xticks]
    for ax in (ax_top, ax_bottom):
        ax.set_xticks(xticks, xlabels)
        ax.set_xlim(min(xticks) - 0.8, max(xticks) + 0.8)
        style_axis(ax)
    ax_top.set_xticklabels([])
    if top_ylabel:
        ax_top.set_ylabel(top_ylabel, fontsize=17.5, fontweight="bold", color=INK, labelpad=12)
    if bottom_ylabel:
        ax_bottom.set_ylabel(bottom_ylabel, fontsize=17.5, fontweight="bold", color=INK, labelpad=12)
    ax_bottom.set_xlabel(r"reco photon $E_T$ bin [GeV]", fontsize=17.5, fontweight="bold", color=INK, labelpad=10)


def render_truth_slide(data: pd.DataFrame) -> tuple[Path, Path, Path, Path]:
    png = OUTDIR / "the8_truth_purity_et_gain_target80_sim_v1.png"
    layout = OUTDIR / "the8_truth_purity_et_gain_target80_sim_v1.layout_nodes.json"
    manifest = OUTDIR / "the8_truth_purity_et_gain_target80_sim_v1.manifest.json"
    script = OUTDIR / "the8_truth_purity_et_gain_target80_sim_v1_speaker_script.md"
    fig = plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI)
    nodes = []
    draw_header(
        fig,
        nodes,
        title=r"Target-80 truth purity versus $E_T$ by centrality",
        metric_line=r"Truth purity = $S_A/(S_A+B_A)$ in isolated-tight region A",
        metric_detail="Gain = BDT-selected purity minus reference box-cut purity.",
    )
    draw_legend(fig, nodes)
    axes = [
        fig.add_axes([0.092, 0.495, 0.400, 0.226]),
        fig.add_axes([0.550, 0.495, 0.400, 0.226]),
        fig.add_axes([0.092, 0.145, 0.400, 0.252]),
        fig.add_axes([0.550, 0.145, 0.400, 0.252]),
    ]
    plot_column(
        axes[0],
        axes[2],
        data,
        "Global BDT target-80",
        top_col="truth_region_a_purity",
        gain_col="truth_purity_gain_pp",
        top_ylabel="Truth purity (%)",
        bottom_ylabel="Gain vs ref. (pp)",
        top_ylim=(84, 94.5),
        gain_ylim=(-2.5, 7.0),
    )
    plot_column(
        axes[1],
        axes[3],
        data,
        "8x7 routed BDT target-80",
        top_col="truth_region_a_purity",
        gain_col="truth_purity_gain_pp",
        top_ylabel=None,
        bottom_ylabel=None,
        top_ylim=(84, 94.5),
        gain_ylim=(-2.5, 7.0),
    )
    add_node(nodes, name="plot_grid", kind="plot", bbox=[235, 405, 2432, 1288])
    fig.savefig(png)
    plt.close(fig)
    manifest.write_text(json.dumps({"schema": "THE8_TRUTH_PURITY_ET_GAIN_TARGET80_SIM_V1", "input_csv": str(INPUT), "png": str(png), "selection": "existing target-80 RecoilJets SIM products; reference box cuts used as denominator for gains"}, indent=2), encoding="utf-8")
    script.write_text(
        "# Speaker script: ET-dependent truth purity\n\n"
        "This slide turns the purity question into an ET-dependent diagnostic. Each column is one available BDT target-80 product, and each curve is a centrality class.\n\n"
        "The top row shows the truth-labeled purity of the selected isolated-tight region A sample. The bottom row subtracts the reference box-cut purity in the same ET and centrality bin, so positive values mean the BDT-selected sample is more truth-prompt dominated.\n\n"
        "The important read is not just whether the BDT is higher on average. The ET dependence tells us whether the improvement is broad or concentrated at the hard end of the photon window, and whether central AuAu behaves differently from peripheral AuAu.\n",
        encoding="utf-8",
    )
    layout.write_text(json.dumps({"nodes": nodes, "canvas": {"width": W, "height": H}}, indent=2), encoding="utf-8")
    return png, layout, manifest, script


def render_closure_slide(data: pd.DataFrame) -> tuple[Path, Path, Path, Path]:
    png = OUTDIR / "the8_abcd_closure_et_gain_target80_sim_v1.png"
    layout = OUTDIR / "the8_abcd_closure_et_gain_target80_sim_v1.layout_nodes.json"
    manifest = OUTDIR / "the8_abcd_closure_et_gain_target80_sim_v1.manifest.json"
    script = OUTDIR / "the8_abcd_closure_et_gain_target80_sim_v1_speaker_script.md"
    fig = plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI)
    nodes = []
    draw_header(
        fig,
        nodes,
        title=r"ABCD purity closure versus $E_T$ by centrality",
        metric_line=r"Closure gain = |ref residual| - |BDT residual|",
        metric_detail="Positive = BDT ABCD estimate is closer to truth.",
    )
    draw_legend(fig, nodes)
    axes = [
        fig.add_axes([0.092, 0.495, 0.400, 0.226]),
        fig.add_axes([0.550, 0.495, 0.400, 0.226]),
        fig.add_axes([0.092, 0.145, 0.400, 0.252]),
        fig.add_axes([0.550, 0.145, 0.400, 0.252]),
    ]
    plot_column(
        axes[0],
        axes[2],
        data,
        "Global BDT target-80",
        top_col="reco_abcd_raw_purity",
        gain_col="closure_abs_residual_improvement_pp",
        top_ylabel="ABCD purity (%)",
        bottom_ylabel="Closure gain (pp)",
        top_ylim=(20, 86),
        gain_ylim=(-31, 42),
    )
    plot_column(
        axes[1],
        axes[3],
        data,
        "8x7 routed BDT target-80",
        top_col="reco_abcd_raw_purity",
        gain_col="closure_abs_residual_improvement_pp",
        top_ylabel=None,
        bottom_ylabel=None,
        top_ylim=(20, 86),
        gain_ylim=(-31, 42),
    )
    add_node(nodes, name="plot_grid", kind="plot", bbox=[235, 405, 2432, 1288])
    fig.savefig(png)
    plt.close(fig)
    manifest.write_text(json.dumps({"schema": "THE8_ABCD_CLOSURE_ET_GAIN_TARGET80_SIM_V1", "input_csv": str(INPUT), "png": str(png), "selection": "existing target-80 RecoilJets SIM products; closure gain is relative to reference box cuts"}, indent=2), encoding="utf-8")
    script.write_text(
        "# Speaker script: ET-dependent ABCD closure\n\n"
        "This slide checks whether the purity estimate itself behaves better, not only whether the selected sample is cleaner.\n\n"
        "The top row is the raw reconstructed ABCD purity estimate. The bottom row converts closure into one positive number: how many percentage points closer the BDT residual is to truth than the reference-box-cut residual in the same ET and centrality bin.\n\n"
        "A positive bottom curve means the BDT makes the ABCD estimate less conservative relative to the truth-labeled purity. This is a more physics-facing diagnostic than AUC, because it asks whether the photon-ID choice produces a purity estimate that we can actually use downstream.\n",
        encoding="utf-8",
    )
    layout.write_text(json.dumps({"nodes": nodes, "canvas": {"width": W, "height": H}}, indent=2), encoding="utf-8")
    return png, layout, manifest, script


def main() -> int:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "font.family": "Times New Roman",
            "mathtext.fontset": "dejavuserif",
            "figure.facecolor": "white",
            "axes.labelcolor": INK,
            "xtick.color": INK,
            "ytick.color": INK,
        }
    )
    data = load_points()
    csv = OUTDIR / "the8_purity_et_closure_target80_sim_v1_datapoints.csv"
    data.to_csv(csv, index=False)
    outputs = [render_truth_slide(data), render_closure_slide(data)]
    summary = {
        "schema": "THE8_PURITY_ET_CLOSURE_SLIDE_SET_V1",
        "input_csv": str(INPUT),
        "datapoints_csv": str(csv),
        "outputs": [
            {"png": str(png), "layout": str(layout), "manifest": str(manifest), "speaker_script": str(script)}
            for png, layout, manifest, script in outputs
        ],
        "available_model_scope": [REFERENCE, *BDT_MODELS],
        "scope_caveat": "This uses available RecoilJets purity products. It is not yet the corrected C3/C7/ET/ETxC routing-family purity closure because those local ROOT/ABCD products were not found.",
    }
    (OUTDIR / "the8_purity_et_closure_slide_set_v1.manifest.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    for png, layout, _, script in outputs:
        print(png)
        print(layout)
        print(script)
    print(csv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
