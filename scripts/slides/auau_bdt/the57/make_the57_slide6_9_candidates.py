#!/usr/bin/env python3
"""Generate THE-57 slide 6-9 candidates from the corrected A1/A2 validation outputs."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np


REPO = Path(__file__).resolve().parents[4]
BASE = REPO / "dataOutput/auauTightBDTValidation/THE57_baseline_validation_20260615"
A1 = BASE / "a1_withcut"
A2 = BASE / "a2_nocut"
OUT = BASE / "slides6_9_candidates"
PRODUCT = "centAsFeatBase3x3_pt15to35"


plt.rcParams.update(
    {
        "font.family": ["Times New Roman", "DejaVu Serif"],
        "axes.titlesize": 17,
        "axes.labelsize": 13.5,
        "xtick.labelsize": 11.5,
        "ytick.labelsize": 11.5,
        "legend.fontsize": 11.2,
    }
)


def read_json(path: Path) -> dict:
    with path.open() as f:
        return json.load(f)


def parse_summary(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    for line in path.read_text().splitlines():
        if "=" in line:
            key, value = line.split("=", 1)
            out[key.strip()] = value.strip()
    return out


def f(summary: dict[str, str], key: str) -> float:
    return float(summary[key])


def cent_label(token: str) -> str:
    return token.replace("_", "-") + "%"


def setup_fig(title: str, bullets: list[str]) -> plt.Figure:
    fig = plt.figure(figsize=(16, 9), dpi=160)
    fig.patch.set_facecolor("white")
    fig.text(0.047, 0.952, "sPHENIX", fontsize=15, fontstyle="italic", fontweight="bold", ha="left", va="top")
    fig.text(0.134, 0.952, "Internal", fontsize=15, ha="left", va="top")
    fig.text(0.047, 0.904, title, fontsize=25, fontweight="bold", ha="left", va="top")
    y = 0.846
    for bullet in bullets:
        fig.text(0.063, y, "\u25b8", fontsize=15, ha="left", va="top", color="#224b7a")
        fig.text(0.082, y, bullet, fontsize=14.8, ha="left", va="top", color="#111111")
        y -= 0.043
    return fig


def density(ax, edges, values, *, color: str, label: str, lw: float = 2.5, ls: str = "-", fill: bool = False, alpha: float = 0.16) -> None:
    edges_arr = np.asarray(edges, dtype=float)
    vals_arr = np.asarray(values, dtype=float)
    if fill:
        ax.stairs(vals_arr, edges_arr, color=color, linewidth=1.3, fill=True, alpha=alpha, label=label)
        ax.stairs(vals_arr, edges_arr, color=color, linewidth=1.7)
    else:
        ax.stairs(vals_arr, edges_arr, color=color, linewidth=lw, linestyle=ls, label=label)


def make_slide6(a1_summary, a2_summary, a1_scores, a2_scores, a1_deep, a2_deep) -> Path:
    fig = setup_fig(
        "THE-58 cleaning preserves the baseline BDT hierarchy",
        [
            "A1 is the default trained with the THE-58 low-calo event cut; A2 is the no-cut diagnostic only.",
            "The score ordering remains stable in all centrality bins: signal high, inclusive MC broad with low-score support.",
            "This checks the event-quality cut did not manufacture the classifier separation.",
        ],
    )
    edges = a1_scores["bin_edges"]
    a1 = a1_scores["products"][PRODUCT]["by_centrality"]
    a2 = a2_scores["products"][PRODUCT]["by_centrality"]
    axes = fig.subplots(1, 3, gridspec_kw={"left": 0.058, "right": 0.98, "bottom": 0.145, "top": 0.655, "wspace": 0.18})
    for ax, token in zip(axes, ["0_20", "20_50", "50_80"]):
        density(ax, edges, a1[token]["background"]["density"], color="#1f77b4", label="A1 inclusive", fill=True, alpha=0.12)
        density(ax, edges, a2[token]["background"]["density"], color="#1f77b4", label="A2 inclusive", lw=2.0, ls="--")
        density(ax, edges, a1[token]["signal"]["density"], color="#e66f00", label="A1 signal", lw=3.0)
        density(ax, edges, a2[token]["signal"]["density"], color="#e66f00", label="A2 signal", lw=2.1, ls="--")
        auc1 = float(a1_deep["products"][PRODUCT]["auc_by_centrality"][token]["auc"])
        auc2 = float(a2_deep["products"][PRODUCT]["auc_by_centrality"][token]["auc"])
        ax.set_title(f"{cent_label(token)}  A1/A2 AUC={auc1:.3f}/{auc2:.3f}", pad=8)
        ax.set_xlim(0, 1)
        ax.set_yscale("log")
        ax.set_ylim(8e-4, 35)
        ax.grid(alpha=0.22)
        ax.set_xlabel("BDT score")
        if ax is axes[0]:
            ax.set_ylabel("Normalized candidates")
        else:
            ax.set_yticklabels([])
    axes[-1].legend(loc="upper right", frameon=True, ncols=1)
    fig.text(
        0.058,
        0.06,
        f"Global AUC: A1={f(a1_summary, PRODUCT + '_auc'):.3f}, A2={f(a2_summary, PRODUCT + '_auc'):.3f}; "
        f"score-gap: A1={f(a1_summary, PRODUCT + '_signal_score_mean') - f(a1_summary, PRODUCT + '_background_score_mean'):.3f}, "
        f"A2={f(a2_summary, PRODUCT + '_signal_score_mean') - f(a2_summary, PRODUCT + '_background_score_mean'):.3f}.",
        fontsize=13.5,
        ha="left",
    )
    out = OUT / "the57_slide6_a1_vs_a2_cut_effect.png"
    fig.savefig(out, facecolor="white")
    plt.close(fig)
    return out


def make_slide7(a1_summary, a1_scores, a1_deep) -> Path:
    fig = setup_fig(
        "Default 14-feature Au+Au BDT has the expected score separation",
        [
            "A1 is the baseline model to hand off: baseV3E + centrality + weta33/wphi33, trained with THE-58 cleaning.",
            "Signal MC concentrates toward high score; inclusive MC retains the low-score background component.",
            "The behavior matches the working-state BDT pattern we wanted to recover.",
        ],
    )
    edges = a1_scores["bin_edges"]
    by_cent = a1_scores["products"][PRODUCT]["by_centrality"]
    auc_by_cent = a1_deep["products"][PRODUCT]["auc_by_centrality"]
    axes = fig.subplots(1, 3, gridspec_kw={"left": 0.058, "right": 0.98, "bottom": 0.145, "top": 0.655, "wspace": 0.18})
    for ax, token in zip(axes, ["0_20", "20_50", "50_80"]):
        density(ax, edges, by_cent[token]["background"]["density"], color="#1f77b4", label="Inclusive MC", fill=True, alpha=0.18)
        density(ax, edges, by_cent[token]["signal"]["density"], color="#e66f00", label="Signal MC", lw=3.0)
        ax.set_title(f"{cent_label(token)}  AUC={auc_by_cent[token]['auc']:.3f}", pad=8)
        ax.set_xlim(0, 1)
        ax.set_yscale("log")
        ax.set_ylim(8e-4, 35)
        ax.grid(alpha=0.22)
        ax.set_xlabel("BDT score")
        if ax is axes[0]:
            ax.set_ylabel("Normalized candidates")
        else:
            ax.set_yticklabels([])
        ax.text(
            0.035,
            0.065,
            f"S={int(by_cent[token]['signal']['entries']):,}\nB={int(by_cent[token]['background']['entries']):,}",
            transform=ax.transAxes,
            fontsize=11.2,
            bbox=dict(facecolor="white", edgecolor="#cccccc", boxstyle="round,pad=0.24", alpha=0.92),
        )
    axes[-1].legend(loc="upper right", frameon=True)
    fig.text(
        0.058,
        0.06,
        f"Validation AUC={f(a1_summary, PRODUCT + '_auc'):.3f}; score means signal/inclusive="
        f"{f(a1_summary, PRODUCT + '_signal_score_mean'):.3f}/{f(a1_summary, PRODUCT + '_background_score_mean'):.3f}.",
        fontsize=13.5,
        ha="left",
    )
    out = OUT / "the57_slide7_default_bdt_score_hierarchy.png"
    fig.savefig(out, facecolor="white")
    plt.close(fig)
    return out


def wp_grids(wp: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    prod = wp["products"][PRODUCT]
    pt_edges = np.asarray(prod["pt_edges"], dtype=float)
    cent_edges = np.asarray(prod["cent_edges"], dtype=float)
    thresholds = np.asarray(prod["grid_thresholds"], dtype=float)
    eff = np.full_like(thresholds, np.nan, dtype=float)
    fake = np.full_like(thresholds, np.nan, dtype=float)
    for cell in prod["cells"]:
        i = list(cent_edges).index(float(cell["centrality_min"]))
        j = list(pt_edges).index(float(cell["pt_min"]))
        eff[i, j] = float(cell["signal_efficiency"])
        fake[i, j] = float(cell["background_fake_rate"])
    return pt_edges, cent_edges, thresholds, eff, fake


def make_heatmap(ax, fig, data, xedges, yedges, *, title, cmap, vmin, vmax, cbar_label, fmt=".2f") -> None:
    mesh = ax.pcolormesh(xedges, yedges, data, cmap=cmap, norm=Normalize(vmin=vmin, vmax=vmax), shading="flat")
    ax.invert_yaxis()
    ax.set_xlabel(r"Cluster $E_T$ [GeV]")
    ax.set_ylabel("Centrality [%]")
    ax.set_title(title, pad=8)
    cbar = fig.colorbar(mesh, ax=ax, fraction=0.046, pad=0.03)
    cbar.set_label(cbar_label)
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            val = data[i, j]
            x = (xedges[j] + xedges[j + 1]) / 2
            y = (yedges[i] + yedges[i + 1]) / 2
            text_color = "white" if (val - vmin) / max(vmax - vmin, 1e-9) < 0.62 else "black"
            ax.text(x, y, format(val, fmt), ha="center", va="center", fontsize=9.5, color=text_color)


def make_slide8(wp80) -> Path:
    fig = setup_fig(
        "WP80 is re-derived from A1 as a centrality-pT grid",
        [
            "The old linear cut is not reused; every cell is set from the A1 signal score distribution.",
            "Thresholds rise toward peripheral bins where the BDT separation is strongest.",
            "This grid is the working-point object to stage with the Blair packet.",
        ],
    )
    prod = wp80["products"][PRODUCT]
    xedges, yedges, thresholds, _, _ = wp_grids(wp80)
    axes = fig.subplots(1, 2, gridspec_kw={"left": 0.065, "right": 0.965, "bottom": 0.145, "top": 0.66, "wspace": 0.23})
    make_heatmap(axes[0], fig, thresholds, xedges, yedges, title="A1 WP80 BDT threshold", cmap="viridis", vmin=0.51, vmax=0.63, cbar_label="BDT cut")
    residual = thresholds - np.asarray(prod["fit_quality"]["plane_fit"]["intercept"], dtype=float)
    residual = thresholds - (
        prod["fit_quality"]["plane_fit"]["intercept"]
        + prod["fit_quality"]["plane_fit"]["centrality_slope"] * ((yedges[:-1] + yedges[1:]) / 2)[:, None]
        + prod["fit_quality"]["plane_fit"]["pt_slope"] * ((xedges[:-1] + xedges[1:]) / 2)[None, :]
    )
    make_heatmap(axes[1], fig, residual, xedges, yedges, title="Residual to simple plane", cmap="coolwarm", vmin=-0.04, vmax=0.04, cbar_label="cut residual", fmt=".3f")
    fig.text(
        0.065,
        0.06,
        f"Mode={prod['mode']}; target={prod['target_signal_efficiency']:.2f}; max |cell eff-0.80|={prod['fit_quality']['max_abs_cell_efficiency_error']:.4g}; min cell S/B={prod['fit_quality']['min_cell_signal_entries']:,}/{prod['fit_quality']['min_cell_background_entries']:,}.",
        fontsize=13.2,
        ha="left",
    )
    out = OUT / "the57_slide8_wp80_threshold_grid.png"
    fig.savefig(out, facecolor="white")
    plt.close(fig)
    return out


def make_slide9(wp80) -> Path:
    fig = setup_fig(
        "WP80 closure: 80% signal efficiency with measured fake rate",
        [
            "The A1 grid lands on 80% signal efficiency in each centrality-pT cell.",
            "The inclusive-MC fake rate is the physics-facing cost of the WP80 photon-ID cut.",
            "This is the corrected replacement for the old centrality-linear threshold slide.",
        ],
    )
    prod = wp80["products"][PRODUCT]
    xedges, yedges, _, eff, fake = wp_grids(wp80)
    axes = fig.subplots(1, 2, gridspec_kw={"left": 0.065, "right": 0.965, "bottom": 0.145, "top": 0.66, "wspace": 0.24})
    make_heatmap(axes[0], fig, eff, xedges, yedges, title="Signal efficiency after WP80", cmap="YlGnBu", vmin=0.7997, vmax=0.8003, cbar_label="signal eff.", fmt=".4f")
    make_heatmap(axes[1], fig, fake, xedges, yedges, title="Inclusive MC passing fraction", cmap="magma", vmin=0.12, vmax=0.36, cbar_label="fake rate", fmt=".2f")
    fig.text(
        0.065,
        0.06,
        f"Inclusive signal efficiency={prod['inclusive']['signal_efficiency']:.6f}; inclusive fake rate={prod['inclusive']['background_fake_rate']:.3f}; mean |cell eff-0.80|={prod['fit_quality']['mean_abs_cell_efficiency_error']:.2e}.",
        fontsize=13.2,
        ha="left",
    )
    out = OUT / "the57_slide9_wp80_efficiency_fake_closure.png"
    fig.savefig(out, facecolor="white")
    plt.close(fig)
    return out


def write_script(path: Path, title: str, lines: list[str]) -> None:
    path.write_text("# " + title + "\n\n" + "\n".join(lines) + "\n")


def write_scripts(slides: list[Path]) -> list[Path]:
    scripts = []
    script_specs = [
        (
            "the57_slide6_a1_vs_a2_cut_effect_script.md",
            "Slide 6 spoken script",
            [
                "This slide checks that the THE-58 event-quality cut is not artificially creating the BDT effect.",
                "The solid curves are the default A1 model trained with the cut, and the dashed curves are the A2 no-cut diagnostic.",
                "In all three centrality bins the signal remains concentrated at high score while inclusive MC keeps the lower-score background support.",
                "That is the working-state behavior we wanted: the cut stabilizes the input sample, but it does not invert or manufacture the classifier hierarchy.",
            ],
        ),
        (
            "the57_slide7_default_bdt_score_hierarchy_script.md",
            "Slide 7 spoken script",
            [
                "This is the default 14-feature Au+Au baseline BDT after the reconstruction-contract fix and THE-58 cleaning.",
                "The signal distribution is high-score dominated, and the inclusive MC is broad with a clear low-score component.",
                "The centrality dependence is sensible: the AUC is lower in central events and improves toward peripheral events.",
                "This is the BDT behavior I would expect from the working-state pipeline.",
            ],
        ),
        (
            "the57_slide8_wp80_threshold_grid_script.md",
            "Slide 8 spoken script",
            [
                "This slide replaces the old linear WP80 line with the A1-derived centrality-pT grid.",
                "Each cell is determined from the A1 signal score distribution at 80 percent signal efficiency.",
                "The right panel shows how far the grid departs from a simple plane; that is a diagnostic, not the working-point formula.",
                "The object to stage with the model is the grid/runtime fragment, not the old centrality-only line.",
            ],
        ),
        (
            "the57_slide9_wp80_efficiency_fake_closure_script.md",
            "Slide 9 spoken script",
            [
                "This is the closure check for the A1 WP80 grid.",
                "The left panel shows the signal efficiency lands on 80 percent cell by cell.",
                "The right panel shows the inclusive MC passing fraction, which is the background cost of the working point.",
                "This is the clean basis for the Blair handoff and for the next photon-side application pass.",
            ],
        ),
    ]
    for filename, title, lines in script_specs:
        script = OUT / filename
        write_script(script, title, lines)
        scripts.append(script)
    return scripts


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    a1_summary = parse_summary(A1 / "validation_summary.txt")
    a2_summary = parse_summary(A2 / "validation_summary.txt")
    a1_scores = read_json(A1 / "score_histograms.json")
    a2_scores = read_json(A2 / "score_histograms.json")
    a1_deep = read_json(A1 / "validation_deep_diagnostics.json")
    a2_deep = read_json(A2 / "validation_deep_diagnostics.json")
    wp80 = read_json(A1 / "bdt_working_points_target80.json")

    slides = [
        make_slide6(a1_summary, a2_summary, a1_scores, a2_scores, a1_deep, a2_deep),
        make_slide7(a1_summary, a1_scores, a1_deep),
        make_slide8(wp80),
        make_slide9(wp80),
    ]
    scripts = write_scripts(slides)
    manifest = {
        "schema": "THE57_SLIDES_6_9_CANDIDATES_V1",
        "status": "READY",
        "product": PRODUCT,
        "default_model": "A1_with_THE58_low_calo_event_cut",
        "diagnostic_model": "A2_without_THE58_cut",
        "source_dirs": {"a1": str(A1), "a2": str(A2)},
        "slides": [str(p) for p in slides],
        "speaker_scripts": [str(p) for p in scripts],
        "notes": [
            "Generated from local pulled A1/A2 validation outputs only.",
            "No Google Slides mutation performed.",
            "A1 is the only default model; A2 is diagnostic only.",
        ],
    }
    manifest_path = OUT / "the57_slides6_9_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print("THE57_SLIDES_6_9_READY")
    print(f"manifest={manifest_path}")
    for p in slides:
        print(f"slide={p}")
    for p in scripts:
        print(f"script={p}")


if __name__ == "__main__":
    main()
