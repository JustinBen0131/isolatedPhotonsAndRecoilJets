#!/usr/bin/env python3
"""Generate THE-57 working-point replacement slide candidates.

Inputs are the A1 default AuAu BDT working-point JSON produced by THE-57.
Outputs are full-slide 16:9 PNG candidates that parallel the old working-point
deck slide 10/11 content without mutating Google Slides.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch, Rectangle


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
BASE = REPO / "dataOutput/auauTightBDTValidation/THE57_baseline_validation_20260615"
WP_JSON = BASE / "a1_withcut/bdt_working_points_target80.json"
OUT_DIR = BASE / "slide10_11_replacements"


BLUE = "#2f5e91"
LIGHT_BLUE = "#dcecfb"
MID_BLUE = "#79a8d6"
GREY = "#666666"
DARK = "#1b1b1b"
GRID_GREY = "#d7d7d7"
GOLD = "#d98c00"


def load_wp() -> dict:
    with WP_JSON.open() as f:
        payload = json.load(f)
    product_name = "centAsFeatBase3x3_pt15to35"
    if "products" not in payload or product_name not in payload["products"]:
        raise KeyError(f"{WP_JSON} does not contain products/{product_name}")
    wp = dict(payload["products"][product_name])
    wp["product"] = product_name
    wp["source_file"] = str(WP_JSON)
    wp["top_level_working_point_mode"] = payload.get("working_point_mode")
    return wp


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "DejaVu Serif", "Times"],
            "font.size": 18,
            "axes.labelsize": 18,
            "axes.titlesize": 20,
            "xtick.labelsize": 15,
            "ytick.labelsize": 15,
            "legend.fontsize": 15,
            "savefig.dpi": 160,
        }
    )


def fig_title(fig: plt.Figure, title: str, subtitle: str) -> None:
    fig.text(0.052, 0.94, title, fontsize=33, fontweight="bold", color=DARK)
    fig.text(0.052, 0.895, subtitle, fontsize=20, color=GREY)
    fig.add_artist(Rectangle((0.052, 0.862), 0.895, 0.012, color=BLUE, transform=fig.transFigure))


def add_note_box(
    fig: plt.Figure,
    x: float,
    y: float,
    w: float,
    h: float,
    title: str,
    lines: list[str],
    face: str = "#f7fbff",
    edge: str = MID_BLUE,
    title_color: str = BLUE,
) -> None:
    box = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.012,rounding_size=0.014",
        transform=fig.transFigure,
        facecolor=face,
        edgecolor=edge,
        linewidth=1.4,
    )
    fig.add_artist(box)
    fig.text(x + 0.018, y + h - 0.05, title, fontsize=21, fontweight="bold", color=title_color)
    yy = y + h - 0.095
    for line in lines:
        fig.text(x + 0.024, yy, line, fontsize=15.5, color=DARK)
        yy -= 0.040


def centrality_labels(edges: list[float]) -> list[str]:
    return [f"{int(edges[i])}-{int(edges[i + 1])}%" for i in range(len(edges) - 1)]


def calc_summary(wp: dict) -> dict:
    cuts = np.asarray(wp["grid_thresholds"], dtype=float)
    eff = np.asarray(wp["grid_signal_efficiency"], dtype=float)
    fake = np.asarray(wp["grid_background_fake_rate"], dtype=float)
    cent_edges = np.asarray(wp["cent_edges"], dtype=float)
    pt_edges = np.asarray(wp["pt_edges"], dtype=float)
    cent_centers = 0.5 * (cent_edges[:-1] + cent_edges[1:])
    pt_centers = 0.5 * (pt_edges[:-1] + pt_edges[1:])
    means = cuts.mean(axis=1)
    mins = cuts.min(axis=1)
    maxs = cuts.max(axis=1)
    stds = cuts.std(axis=1)
    fake_mean = fake.mean(axis=1)
    eff_mean = eff.mean(axis=1)
    line_coef = np.polyfit(cent_centers, means, 1)
    line_slope = float(line_coef[0])
    line_intercept = float(line_coef[1])
    line_pred = line_intercept + line_slope * cent_centers
    line_rms = float(np.sqrt(np.mean((means - line_pred) ** 2)))
    return {
        "cuts": cuts,
        "eff": eff,
        "fake": fake,
        "cent_edges": cent_edges,
        "pt_edges": pt_edges,
        "cent_centers": cent_centers,
        "pt_centers": pt_centers,
        "cent_labels": centrality_labels(cent_edges.tolist()),
        "means": means,
        "mins": mins,
        "maxs": maxs,
        "stds": stds,
        "fake_mean": fake_mean,
        "eff_mean": eff_mean,
        "line_intercept": line_intercept,
        "line_slope": line_slope,
        "line_rms": line_rms,
    }


def draw_slide10(wp: dict, summary: dict) -> Path:
    fig = plt.figure(figsize=(16, 9), facecolor="white")
    fig_title(
        fig,
        "Flat summaries of the 80% signal-efficiency BDT threshold",
        r"Default Au+Au baseline: 14 features, low-calo cleaned, 15 < $E_T$ < 35 GeV",
    )

    cuts = summary["cuts"]
    pt_centers = summary["pt_centers"]
    labels = summary["cent_labels"]
    means = summary["means"]
    fake_mean = summary["fake_mean"]
    eff_mean = summary["eff_mean"]
    colors = ["#1f77b4", "#2ca02c", "#9467bd"]

    lefts = [0.065, 0.37, 0.675]
    for i, left in enumerate(lefts):
        ax = fig.add_axes([left, 0.345, 0.255, 0.43])
        ax.plot(pt_centers, cuts[i], "o", color=colors[i], ms=8, label="WP80 cells")
        ax.axhline(means[i], color=GOLD, lw=2.4, label="flat summary")
        ax.fill_between([15, 35], summary["mins"][i], summary["maxs"][i], color=colors[i], alpha=0.10, lw=0)
        ax.set_xlim(14.5, 35.5)
        ax.set_ylim(0.50, 0.64)
        ax.set_title(labels[i], fontsize=21, fontweight="bold", pad=10)
        ax.set_xlabel(r"cluster $E_T$ [GeV]")
        if i == 0:
            ax.set_ylabel(r"BDT threshold $T_{80}$")
        else:
            ax.set_yticklabels([])
        ax.grid(True, color=GRID_GREY, lw=0.8, alpha=0.75)
        ax.tick_params(direction="in", top=True, right=True)
        ax.text(
            0.04,
            0.93,
            f"flat = {means[i]:.4f}\n"
            f"avg eff = {eff_mean[i]:.5f}\n"
            f"avg fake = {fake_mean[i]:.3f}",
            transform=ax.transAxes,
            fontsize=13.5,
            va="top",
            bbox=dict(facecolor="white", edgecolor="#bbbbbb", boxstyle="round,pad=0.25", alpha=0.9),
        )

    add_note_box(
        fig,
        0.062,
        0.095,
        0.875,
        0.17,
        "Working-point readout",
        [
            f"Exact working point is a 3 x 8 centrality-pT grid; max abs cell efficiency error = {wp['fit_quality']['max_abs_cell_efficiency_error']:.6f}.",
            "Flat lines summarize each centrality bin for slide comparison; final application uses the grid lookup, not a single flat number.",
            f"Inclusive A1 performance: signal efficiency = {wp['inclusive']['signal_efficiency']:.6f}, inclusive-MC fake rate = {wp['inclusive']['background_fake_rate']:.3f}.",
        ],
        face="#f9fbfd",
        edge="#c9d7e6",
        title_color=BLUE,
    )

    out = OUT_DIR / "the57_slide10_flat_points_replacement.png"
    fig.savefig(out, facecolor="white")
    plt.close(fig)
    return out


def draw_slide11(wp: dict, summary: dict) -> Path:
    fig = plt.figure(figsize=(16, 9), facecolor="white")
    fig_title(
        fig,
        "The 80% efficiency grid defines the tight BDT selection",
        "A centrality trend is retained for comparison with the previous working-point slide",
    )

    ax = fig.add_axes([0.075, 0.19, 0.60, 0.59])
    cent_centers = summary["cent_centers"]
    means = summary["means"]
    yerr = np.vstack([means - summary["mins"], summary["maxs"] - means])
    ax.errorbar(
        cent_centers,
        means,
        yerr=yerr,
        fmt="o",
        ms=10,
        mfc=BLUE,
        mec="white",
        mew=1.4,
        ecolor=BLUE,
        elinewidth=2.2,
        capsize=6,
        label="centrality summary",
    )
    for i, c in enumerate(cent_centers):
        jitter = np.linspace(-2.6, 2.6, summary["cuts"].shape[1])
        ax.scatter(
            np.full(summary["cuts"].shape[1], c) + jitter,
            summary["cuts"][i],
            s=34,
            color=MID_BLUE,
            alpha=0.55,
            edgecolor="none",
            label="pT cells" if i == 0 else None,
        )

    xfit = np.linspace(0, 80, 200)
    yfit = summary["line_intercept"] + summary["line_slope"] * xfit
    ax.plot(xfit, yfit, color=GOLD, lw=3.0, label="linear summary")
    ax.set_xlim(-2, 82)
    ax.set_ylim(0.50, 0.65)
    ax.set_xlabel("centrality percentile c")
    ax.set_ylabel(r"BDT threshold $T_{80}$")
    ax.grid(True, color=GRID_GREY, lw=0.8, alpha=0.75)
    ax.tick_params(direction="in", top=True, right=True)
    ax.legend(frameon=False, loc="upper left")
    ax.text(
        0.05,
        0.19,
        f"centrality summary:\n"
        f"T80(c) = {summary['line_intercept']:.4f} + {summary['line_slope']:.6f} c\n"
        f"RMS residual = {summary['line_rms']:.4f}",
        transform=ax.transAxes,
        fontsize=16,
        bbox=dict(facecolor="white", edgecolor="#bbbbbb", boxstyle="round,pad=0.35", alpha=0.92),
    )

    fit = wp.get("fit_quality", {})
    plane = fit.get("plane_fit", {})
    plane_line = (
        f"T80(c,E_T) = {plane.get('intercept', math.nan):.4f} "
        f"+ {plane.get('centrality_slope', math.nan):.6f} c"
    )
    plane_pt_line = f"{plane.get('pt_slope', math.nan):+.6f} E_T"
    add_note_box(
        fig,
        0.715,
        0.505,
        0.245,
        0.275,
        "Final working definition",
        [
            r"Tight photon ID:",
            r"BDT score > T80(c,$E_T$)",
            "Exact object: 3 x 8 grid lookup",
            "A1 is the default baseline Au+Au model",
            f"Grid max residual to plane = {plane.get('max_abs_residual', math.nan):.4f}",
        ],
        face="#f6fbff",
        edge=MID_BLUE,
        title_color=BLUE,
    )
    add_note_box(
        fig,
        0.715,
        0.190,
        0.245,
        0.265,
        "Diagnostic plane summary",
        [
            plane_line,
            plane_pt_line,
            "c = centrality percentile",
            r"$E_T$ in GeV",
            f"$E_T$ bins: {len(summary['pt_centers'])}, centrality bins: {len(summary['cent_centers'])}",
        ],
        face="#fffaf2",
        edge="#e7c88f",
        title_color="#a06500",
    )

    out = OUT_DIR / "the57_slide11_wp80_fit_replacement.png"
    fig.savefig(out, facecolor="white")
    plt.close(fig)
    return out


def write_manifest(wp: dict, summary: dict, outputs: list[Path]) -> Path:
    manifest = {
        "inputs": {
            "wp_json": str(WP_JSON),
            "product": wp.get("product"),
            "working_point_mode": wp.get("working_point_mode"),
            "model_dir": wp.get("model_dir"),
        },
        "outputs": [str(p) for p in outputs],
        "summary": {
            "centrality_labels": summary["cent_labels"],
            "pt_edges": summary["pt_edges"].tolist(),
            "centrality_edges": summary["cent_edges"].tolist(),
            "flat_thresholds": [float(x) for x in summary["means"]],
            "threshold_min": [float(x) for x in summary["mins"]],
            "threshold_max": [float(x) for x in summary["maxs"]],
            "centrality_linear_summary": {
                "intercept": summary["line_intercept"],
                "slope": summary["line_slope"],
                "rms_residual": summary["line_rms"],
            },
            "inclusive": wp.get("inclusive", {}),
            "fit_quality": wp.get("fit_quality", {}),
        },
    }
    out = OUT_DIR / "the57_slide10_11_replacements_manifest.json"
    out.write_text(json.dumps(manifest, indent=2) + "\n")
    return out


def main() -> None:
    setup_style()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    wp = load_wp()
    summary = calc_summary(wp)
    outputs = [draw_slide10(wp, summary), draw_slide11(wp, summary)]
    manifest = write_manifest(wp, summary, outputs)
    print("Generated:")
    for out in outputs:
        print(out)
    print(manifest)


if __name__ == "__main__":
    main()
