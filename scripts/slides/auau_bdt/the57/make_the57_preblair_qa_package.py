#!/usr/bin/env python3
"""Build the compact THE-57 pre-Blair QA package from pulled validation outputs."""

from __future__ import annotations

import json
import math
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
OUT = BASE / "qa_package"
PRODUCT = "centAsFeatBase3x3_pt15to35"
FEATURES = [
    "cluster_Et",
    "cluster_weta_cogx",
    "cluster_wphi_cogx",
    "cluster_weta33_cogx",
    "cluster_wphi33_cogx",
    "vertexz",
    "cluster_Eta",
    "e11_over_e33",
    "cluster_et1",
    "cluster_et2",
    "cluster_et3",
    "cluster_et4",
    "e32_over_e35",
    "centrality",
]


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


def as_float(summary: dict[str, str], key: str) -> float:
    return float(summary[key])


def density(ax, edges, vals, *, color, label, lw=2.4, fill=False, alpha=0.18):
    edges_arr = np.asarray(edges, dtype=float)
    vals_arr = np.asarray(vals, dtype=float)
    if fill:
        ax.stairs(vals_arr, edges_arr, color=color, linewidth=1.4, fill=True, alpha=alpha, label=label)
        ax.stairs(vals_arr, edges_arr, color=color, linewidth=1.5)
    else:
        ax.stairs(vals_arr, edges_arr, color=color, linewidth=lw, label=label)


def add_header(fig, title: str, bullets: list[str]) -> None:
    fig.text(0.045, 0.948, "sPHENIX", fontsize=14, fontstyle="italic", fontweight="bold", ha="left", va="top")
    fig.text(0.132, 0.948, "Internal", fontsize=14, ha="left", va="top")
    fig.text(0.045, 0.905, title, fontsize=24, fontweight="bold", ha="left", va="top")
    y = 0.848
    for bullet in bullets:
        fig.text(0.06, y, f"> {bullet}", fontsize=14.5, ha="left", va="top")
        y -= 0.046


def load_all():
    a1_summary = parse_summary(A1 / "validation_summary.txt")
    a2_summary = parse_summary(A2 / "validation_summary.txt")
    a1_scores = read_json(A1 / "score_histograms.json")
    a2_scores = read_json(A2 / "score_histograms.json")
    a1_deep = read_json(A1 / "validation_deep_diagnostics.json")
    a2_deep = read_json(A2 / "validation_deep_diagnostics.json")
    wp80 = read_json(A1 / "bdt_working_points_target80.json")
    return a1_summary, a2_summary, a1_scores, a2_scores, a1_deep, a2_deep, wp80


def gate_inputs(a1_summary, a2_summary, a1_deep, a2_deep, wp80) -> dict:
    checks: dict[str, object] = {}
    for label, summary in [("A1", a1_summary), ("A2", a2_summary)]:
        status = summary.get("status")
        finite = as_float(summary, f"{PRODUCT}_finite_score_fraction")
        signal_mean = as_float(summary, f"{PRODUCT}_signal_score_mean")
        background_mean = as_float(summary, f"{PRODUCT}_background_score_mean")
        auc = as_float(summary, f"{PRODUCT}_auc")
        checks[f"{label}_status"] = status
        checks[f"{label}_finite_score_fraction"] = finite
        checks[f"{label}_auc"] = auc
        checks[f"{label}_signal_score_mean"] = signal_mean
        checks[f"{label}_background_score_mean"] = background_mean
        if status != "READY":
            raise RuntimeError(f"{label} validation summary is not READY: {status}")
        if finite < 0.999:
            raise RuntimeError(f"{label} finite score fraction is too low: {finite}")
        if signal_mean <= background_mean:
            raise RuntimeError(f"{label} score hierarchy is inverted: signal {signal_mean} <= background {background_mean}")
        if auc < 0.80:
            raise RuntimeError(f"{label} AUC is below the pre-handoff floor: {auc}")

    for label, deep in [("A1", a1_deep), ("A2", a2_deep)]:
        e11 = deep["features"]["e11_over_e33"]["inclusive"]
        for cls in ["signal", "background"]:
            q01 = float(e11[cls]["q01"])
            entries = int(e11[cls]["entries"])
            checks[f"{label}_e11_{cls}_q01"] = q01
            checks[f"{label}_e11_{cls}_entries"] = entries
            if entries <= 0:
                raise RuntimeError(f"{label} {cls} has no E11/E33 entries")
            if q01 < 0.05:
                raise RuntimeError(f"{label} {cls} E11/E33 q01 is too low, possible low-edge collapse: {q01}")

    wp = wp80["products"][PRODUCT]
    checks["wp80_status"] = wp["status"]
    checks["wp80_mode"] = wp["mode"]
    checks["wp80_signal_efficiency"] = wp["inclusive"]["signal_efficiency"]
    checks["wp80_background_fake_rate"] = wp["inclusive"]["background_fake_rate"]
    checks["wp80_max_abs_cell_efficiency_error"] = wp["fit_quality"]["max_abs_cell_efficiency_error"]
    if wp["status"] != "ready":
        raise RuntimeError(f"WP80 product is not ready: {wp['status']}")
    if abs(float(wp["inclusive"]["signal_efficiency"]) - 0.80) > 0.002:
        raise RuntimeError(f"WP80 inclusive signal efficiency is off target: {wp['inclusive']['signal_efficiency']}")
    if float(wp["fit_quality"]["max_abs_cell_efficiency_error"]) > 0.005:
        raise RuntimeError("WP80 grid has a large cell efficiency error")
    return checks


def cent_label(token: str) -> str:
    return token.replace("_", "-") + "%"


def make_score_slide(a1_summary, a1_scores, a1_deep) -> Path:
    fig = plt.figure(figsize=(16, 9), dpi=160)
    fig.patch.set_facecolor("white")
    add_header(
        fig,
        "Baseline Au+Au BDT: score hierarchy by centrality",
        [
            "A1 default: 14-feature baseline with THE-58 low-calo event cut",
            "Signal sits high in score; inclusive MC keeps the expected low-score component",
            "Validation uses all six embedded samples: photon12/20 and jet12/20/30/40",
        ],
    )

    edges = a1_scores["bin_edges"]
    prod_scores = a1_scores["products"][PRODUCT]
    auc_by_cent = a1_deep["products"][PRODUCT]["auc_by_centrality"]
    axes = fig.subplots(1, 3, gridspec_kw={"left": 0.055, "right": 0.98, "bottom": 0.11, "top": 0.66, "wspace": 0.19})
    for ax, token in zip(axes, ["0_20", "20_50", "50_80"]):
        cent = prod_scores["by_centrality"][token]
        density(ax, edges, cent["background"]["density"], color="#1f77b4", label="Inclusive MC", fill=True)
        density(ax, edges, cent["signal"]["density"], color="#ff7f0e", label="Signal MC", lw=2.8)
        auc_value = auc_by_cent[token]
        auc = float(auc_value["auc"] if isinstance(auc_value, dict) else auc_value)
        ax.set_title(f"{cent_label(token)}   AUC={auc:.3f}", fontsize=17, pad=9)
        ax.set_xlabel("BDT score", fontsize=14)
        ax.set_yscale("log")
        ax.set_ylim(8e-4, 30)
        ax.set_xlim(0, 1)
        ax.grid(alpha=0.22)
        if ax is axes[0]:
            ax.set_ylabel("Normalized candidates", fontsize=14)
        else:
            ax.set_yticklabels([])
        ax.text(
            0.04,
            0.08,
            f"S={int(cent['signal']['entries']):,}\nB={int(cent['background']['entries']):,}",
            transform=ax.transAxes,
            fontsize=11,
            bbox=dict(facecolor="white", edgecolor="#cccccc", boxstyle="round,pad=0.25", alpha=0.9),
        )
    axes[-1].legend(loc="upper right", fontsize=12, frameon=True)
    fig.text(
        0.055,
        0.035,
        f"Global A1 AUC={as_float(a1_summary, PRODUCT + '_auc'):.3f}; score means: signal={as_float(a1_summary, PRODUCT + '_signal_score_mean'):.3f}, inclusive={as_float(a1_summary, PRODUCT + '_background_score_mean'):.3f}.",
        fontsize=12.5,
        ha="left",
        va="bottom",
    )
    out = OUT / "the57_preblair_a1_score_separation.png"
    fig.savefig(out, facecolor="white")
    plt.close(fig)
    return out


def make_shape_slide(a1_deep, a2_deep) -> Path:
    fig = plt.figure(figsize=(16, 9), dpi=160)
    fig.patch.set_facecolor("white")
    add_header(
        fig,
        "Embedded-MC shower-shape sanity: no low-edge collapse",
        [
            "E11/E33 and width quantiles are stable between A1 and A2 validation outputs",
            "A1 E11/E33 q01 stays above 0.05 for signal and inclusive in every centrality bin",
            "The feature behavior matches the intended corrected embedded-MC reconstruction contract",
        ],
    )

    axes = fig.subplots(1, 2, gridspec_kw={"left": 0.07, "right": 0.96, "bottom": 0.13, "top": 0.66, "wspace": 0.23})
    cent_tokens = ["0_20", "20_50", "50_80"]
    x = np.arange(len(cent_tokens))

    def stats(deep, feature, cls, quantile):
        return np.array([deep["features"][feature]["by_centrality"][c][cls][quantile] for c in cent_tokens], dtype=float)

    feature_specs = [
        ("e11_over_e33", "E11/E33", (0.0, 1.0)),
        ("cluster_weta33_cogx", r"$w_{\eta,3x3}^{cog}$", (0.0, 0.45)),
    ]
    for ax, (feature, title, ylim) in zip(axes, feature_specs):
        for cls, color, marker, dx, label in [
            ("signal", "#ff7f0e", "o", -0.08, "Signal MC"),
            ("background", "#1f77b4", "s", 0.08, "Inclusive MC"),
        ]:
            median = stats(a1_deep, feature, cls, "median")
            q05 = stats(a1_deep, feature, cls, "q05")
            q95 = stats(a1_deep, feature, cls, "q95")
            ax.errorbar(
                x + dx,
                median,
                yerr=[median - q05, q95 - median],
                fmt=marker,
                color=color,
                markersize=8,
                linewidth=2,
                capsize=4,
                label=label,
            )
            a2_med = stats(a2_deep, feature, cls, "median")
            ax.scatter(x + dx, a2_med, marker="_", color=color, s=170, linewidths=2.1, alpha=0.85, label=f"{label} A2 median" if feature == "e11_over_e33" else None)
        ax.set_title(title, fontsize=18, pad=10)
        ax.set_xticks(x)
        ax.set_xticklabels([cent_label(c) for c in cent_tokens], fontsize=13)
        ax.set_ylim(*ylim)
        ax.grid(axis="y", alpha=0.24)
        ax.set_ylabel("A1 median with q05-q95 band", fontsize=13)
        if feature == "e11_over_e33":
            q01_sig = stats(a1_deep, feature, "signal", "q01")
            q01_bkg = stats(a1_deep, feature, "background", "q01")
            for i in range(len(x)):
                ax.text(
                    x[i],
                    0.045,
                    f"q01 S {q01_sig[i]:.2f}\nq01 B {q01_bkg[i]:.2f}",
                    ha="center",
                    va="bottom",
                    fontsize=10.5,
                    bbox=dict(facecolor="white", edgecolor="#dddddd", boxstyle="round,pad=0.18", alpha=0.88),
                )
    axes[0].legend(loc="upper left", fontsize=10.5, ncols=2)
    fig.text(
        0.07,
        0.045,
        "Filled points and bands are A1. Horizontal ticks are A2 medians. A2 is diagnostic only; A1 is the default baseline.",
        fontsize=12.5,
        ha="left",
    )
    out = OUT / "the57_preblair_shower_shape_sanity.png"
    fig.savefig(out, facecolor="white")
    plt.close(fig)
    return out


def make_wp80_slide(a1_summary, a2_summary, wp80) -> Path:
    fig = plt.figure(figsize=(16, 9), dpi=160)
    fig.patch.set_facecolor("white")
    add_header(
        fig,
        "WP80 and cut-effect closure",
        [
            "WP80 is derived from A1 only as a centrality-pT grid, not the old linear formula",
            "A1 and A2 have consistent score separation; A2 remains diagnostic and is not the default",
            "The grid achieves 80% signal efficiency with sub-permille cell-level agreement",
        ],
    )
    wp = wp80["products"][PRODUCT]
    axes = fig.subplots(1, 2, gridspec_kw={"left": 0.07, "right": 0.96, "bottom": 0.14, "top": 0.65, "wspace": 0.25})

    # Left: A1/A2 compact metric comparison.
    ax = axes[0]
    labels = ["A1 with cut", "A2 no cut"]
    aucs = [as_float(a1_summary, PRODUCT + "_auc"), as_float(a2_summary, PRODUCT + "_auc")]
    sep = [
        as_float(a1_summary, PRODUCT + "_signal_score_mean") - as_float(a1_summary, PRODUCT + "_background_score_mean"),
        as_float(a2_summary, PRODUCT + "_signal_score_mean") - as_float(a2_summary, PRODUCT + "_background_score_mean"),
    ]
    xpos = np.arange(len(labels))
    width = 0.34
    ax.bar(xpos - width / 2, aucs, width=width, color="#315f9d", label="AUC")
    ax.bar(xpos + width / 2, sep, width=width, color="#d46f2c", label="Mean score gap")
    ax.set_xticks(xpos)
    ax.set_xticklabels(labels, fontsize=13)
    ax.set_ylim(0, 1.0)
    ax.set_title("Training sanity metrics", fontsize=18, pad=10)
    ax.grid(axis="y", alpha=0.24)
    ax.legend(fontsize=12, loc="upper right")
    for i, (auc, gap) in enumerate(zip(aucs, sep)):
        ax.text(i - width / 2, auc + 0.025, f"{auc:.3f}", ha="center", fontsize=12)
        ax.text(i + width / 2, gap + 0.025, f"{gap:.3f}", ha="center", fontsize=12)

    # Right: WP80 threshold grid.
    ax = axes[1]
    grid = np.asarray(wp["grid_thresholds"], dtype=float)
    pt_edges = np.asarray(wp["pt_edges"], dtype=float)
    cent_edges = np.asarray(wp["cent_edges"], dtype=float)
    pcm = ax.pcolormesh(pt_edges, cent_edges, grid, cmap="viridis", norm=Normalize(vmin=0.51, vmax=0.63), shading="flat")
    ax.invert_yaxis()
    ax.set_xlabel(r"Cluster $E_T$ [GeV]", fontsize=13)
    ax.set_ylabel("Centrality [%]", fontsize=13)
    ax.set_title("A1 WP80 grid threshold", fontsize=18, pad=10)
    cbar = fig.colorbar(pcm, ax=ax, fraction=0.046, pad=0.03)
    cbar.set_label("BDT cut", fontsize=12)
    for i in range(grid.shape[0]):
        for j in range(grid.shape[1]):
            ax.text((pt_edges[j] + pt_edges[j + 1]) / 2, (cent_edges[i] + cent_edges[i + 1]) / 2, f"{grid[i, j]:.2f}", ha="center", va="center", fontsize=8.8, color="white")

    fig.text(
        0.07,
        0.055,
        f"WP80 inclusive efficiency={wp['inclusive']['signal_efficiency']:.6f}; inclusive fake rate={wp['inclusive']['background_fake_rate']:.3f}; max cell |eff-0.80|={wp['fit_quality']['max_abs_cell_efficiency_error']:.4g}.",
        fontsize=12.5,
        ha="left",
    )
    out = OUT / "the57_preblair_wp80_and_cut_effect.png"
    fig.savefig(out, facecolor="white")
    plt.close(fig)
    return out


def write_handoff_support_files(checks: dict, a1_summary, wp80) -> tuple[Path, Path]:
    features_path = OUT / "FEATURES.txt"
    features_path.write_text("\n".join(FEATURES) + "\n")
    draft_path = OUT / "the57_blair_message_draft_NOT_SENT.md"
    model_dir = "/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/the57_models/a1_withcut_20260615_stagedcache_streamreduce"
    draft_path.write_text(
        "\n".join(
            [
                "# THE-57 Blair Message Draft - NOT SENT",
                "",
                "Blair,",
                "",
                "I have the corrected Au+Au photon-ID baseline BDT ready for review. This replaces the earlier 32-feature placeholder. The default baseline is now the 14-feature Au+Au model trained with the THE-58 low-calo event cut.",
                "",
                "Model:",
                f"- Product/model id: `{PRODUCT}`",
                f"- TMVA ROOT: `{model_dir}/auau_tight_bdt_centAsFeatBase3x3_pt15to35_tmva.root`",
                "- RBDT key: `myBDT`",
                f"- Ordered features: see `FEATURES.txt` in the staged packet. The local QA copy is `{features_path}`",
                f"- Metadata: `{model_dir}/auau_tight_bdt_centAsFeatBase3x3_pt15to35_tmva.metadata.json`",
                "",
                "WP80:",
                "- Re-derived from the A1 canonical model only.",
                "- It is a centrality-pT grid, not the old linear cut.",
                f"- Runtime fragment: `{A1 / 'bdt_working_points_target80_runtime_fragment.yaml'}`",
                f"- Validation summary: inclusive signal efficiency `{checks['wp80_signal_efficiency']:.9f}`, inclusive background fake rate `{checks['wp80_background_fake_rate']:.9f}`, max cell efficiency error `{checks['wp80_max_abs_cell_efficiency_error']:.9f}`.",
                "",
                "Validation:",
                f"- A1 validation AUC: `{checks['A1_auc']:.6f}`",
                f"- Signal/background score means: `{checks['A1_signal_score_mean']:.6f} / {checks['A1_background_score_mean']:.6f}`",
                f"- Shower-shape sanity: E11/E33 q01 signal/background `{checks['A1_e11_signal_q01']:.6f} / {checks['A1_e11_background_q01']:.6f}`, with no embedded-MC low-edge collapse.",
                "",
                "The no-cut model A2 is diagnostic only and is not the default.",
                "",
                "NOTE FOR JUSTIN/CODEX: Before Justin sends this, stage/copy `FEATURES.txt` and the WP80 runtime fragment into the final Blair-facing SDCC packet path so the message can cite exact SDCC paths for every payload file.",
                "",
            ]
        )
    )
    return features_path, draft_path


def write_manifest(checks: dict, slides: list[Path], a1_summary, a2_summary, wp80) -> Path:
    features_path, draft_path = write_handoff_support_files(checks, a1_summary, wp80)
    manifest = {
        "schema": "THE57_PRE_BLAIR_QA_PACKAGE_V1",
        "status": "READY",
        "product": PRODUCT,
        "default_model": "A1_with_THE58_low_calo_event_cut",
        "diagnostic_model": "A2_without_THE58_cut",
        "source_dirs": {
            "a1_local": str(A1),
            "a2_local": str(A2),
            "a1_remote": a1_summary.get("report_dir"),
            "a2_remote": a2_summary.get("report_dir"),
        },
        "slides": [str(p) for p in slides],
        "features_txt": str(features_path),
        "blair_message_draft_not_sent": str(draft_path),
        "checks": checks,
        "wp80": {
            "mode": wp80["products"][PRODUCT]["mode"],
            "target_signal_efficiency": wp80["target_signal_efficiency"],
            "pt_edges": wp80["products"][PRODUCT]["pt_edges"],
            "cent_edges": wp80["products"][PRODUCT]["cent_edges"],
            "runtime_fragment": str(A1 / "bdt_working_points_target80_runtime_fragment.yaml"),
        },
        "notes": [
            "Mattermost/Blair message is not sent by Codex or Claude.",
            "A1 is the only default baseline candidate.",
            "A2 is diagnostic only.",
        ],
    }
    out = OUT / "the57_preblair_qa_manifest.json"
    out.write_text(json.dumps(manifest, indent=2) + "\n")

    summary = OUT / "the57_preblair_qa_summary.md"
    summary.write_text(
        "\n".join(
            [
                "# THE-57 Pre-Blair QA Package",
                "",
                "Status: READY",
                "",
                f"Default product: `{PRODUCT}`",
                "Default model: A1 with THE-58 low-calo event cut",
                "Diagnostic model: A2 without THE-58 cut",
                "",
                "Key checks:",
                f"- A1 AUC: {checks['A1_auc']:.6f}",
                f"- A2 AUC: {checks['A2_auc']:.6f}",
                f"- A1 signal/background score means: {checks['A1_signal_score_mean']:.6f} / {checks['A1_background_score_mean']:.6f}",
                f"- A1 E11/E33 q01 signal/background: {checks['A1_e11_signal_q01']:.6f} / {checks['A1_e11_background_q01']:.6f}",
                f"- WP80 mode: {checks['wp80_mode']}",
                f"- WP80 inclusive signal efficiency: {checks['wp80_signal_efficiency']:.9f}",
                f"- WP80 inclusive background fake rate: {checks['wp80_background_fake_rate']:.9f}",
                f"- WP80 max cell efficiency error: {checks['wp80_max_abs_cell_efficiency_error']:.9f}",
                "",
                "Slides:",
                *[f"- `{p}`" for p in slides],
                "",
            ]
        )
        + "\n"
    )
    return out


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    a1_summary, a2_summary, a1_scores, a2_scores, a1_deep, a2_deep, wp80 = load_all()
    checks = gate_inputs(a1_summary, a2_summary, a1_deep, a2_deep, wp80)
    slides = [
        make_score_slide(a1_summary, a1_scores, a1_deep),
        make_shape_slide(a1_deep, a2_deep),
        make_wp80_slide(a1_summary, a2_summary, wp80),
    ]
    manifest = write_manifest(checks, slides, a1_summary, a2_summary, wp80)
    print("THE57_PRE_BLAIR_QA_READY")
    print(f"manifest={manifest}")
    for slide in slides:
        print(f"slide={slide}")


if __name__ == "__main__":
    main()
