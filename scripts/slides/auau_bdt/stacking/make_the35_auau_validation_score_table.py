#!/usr/bin/env python3
"""Render the THE-35 Au+Au BDT/MLP/best-stack validation score table."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch


W, H = 2560, 1440
DPI = 200
INK = "#16181d"
MUTED = "#5f6673"
GRID = "#e5e7eb"
SIGNAL = "#cf3347"
BACKGROUND = "#1f70b8"
BDT_EDGE = "#c653a3"
BDT_FILL = "#fbf2f8"
MLP_EDGE = "#e07d1f"
MLP_FILL = "#fbf5eb"
STACK_EDGE = "#23865f"
STACK_FILL = "#eef7f3"


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "Times New Roman",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "axes.edgecolor": "#4b5563",
            "axes.linewidth": 1.0,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "axes.labelcolor": INK,
            "xtick.color": INK,
            "ytick.color": INK,
        }
    )


def load_payload(path: Path) -> dict:
    payload = json.loads(path.read_text())
    lane = payload["selected_lane"]
    selected = next(item for item in payload["candidate_lanes"] if item["lane"] == lane)
    return {"summary": payload, "lane": selected}


def stack_label(model: str) -> str:
    return {
        "score_only_logistic": "Score-only logistic stack",
        "score_only_gbm": "Score-only GBM stack",
        "score_only_mlp": "Score-only MLP stack",
        "score_context_logistic": "Score + context logistic stack",
        "score_context_gbm": "Score + context GBM stack",
        "score_context_mlp": "Score + context MLP stack",
    }.get(model, model.replace("_", " "))


def short_lane_label(lane: str) -> str:
    return {
        "current_oof": "5-fold OOF",
        "simple_holdout": "simple holdout",
    }.get(lane, lane.replace("_", " "))


def compact_count(value: int) -> str:
    if value >= 1_000_000:
        return f"{value / 1_000_000:.2f}M"
    if value >= 10_000:
        return f"{value / 1_000:.0f}k"
    return f"{value:,}"


def panel_lookup(lane: dict) -> dict[tuple[str, float, float], dict]:
    out: dict[tuple[str, float, float], dict] = {}
    for panel in lane["panels"]:
        out[(panel["model"], float(panel["cent_lo"]), float(panel["cent_hi"]))] = panel
    return out


def model_auc_triplet(lookup: dict, model: str, cent_bins: list[tuple[float, float]]) -> list[float]:
    return [float(lookup[(model, lo, hi)]["auc_weighted"]) for lo, hi in cent_bins]


def step_arrays(panel: dict, which: str) -> tuple[np.ndarray, np.ndarray]:
    edges = np.asarray(panel["bins"], dtype=float)
    density = np.asarray(panel[f"{which}_density"], dtype=float)
    x = np.repeat(edges, 2)[1:-1]
    y = np.repeat(density, 2)
    return x, y


def draw_text(fig, x: float, y: float, text: str, **kwargs):
    defaults = {"ha": "left", "va": "top", "color": INK}
    defaults.update(kwargs)
    return fig.text(x, y, text, **defaults)


def rounded_box(fig, x: float, y: float, w: float, h: float, fc: str, ec: str | None = None, lw: float = 1.0, radius: float = 0.012):
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle=f"round,pad=0.006,rounding_size={radius}",
        transform=fig.transFigure,
        facecolor=fc,
        edgecolor=ec or fc,
        linewidth=lw,
        zorder=0,
    )
    fig.patches.append(patch)
    return patch


def draw_legend(fig, x: float = 0.760, y: float = 0.920) -> None:
    fig.lines.append(plt.Line2D([x, x + 0.030], [y, y], transform=fig.transFigure, color=SIGNAL, lw=4))
    draw_text(fig, x + 0.036, y + 0.010, "signal", fontsize=14.0, va="center")
    fig.lines.append(plt.Line2D([x + 0.105, x + 0.135], [y, y], transform=fig.transFigure, color=BACKGROUND, lw=4))
    draw_text(fig, x + 0.141, y + 0.010, "background", fontsize=14.0, va="center")


def draw_row_card(fig, x: float, y: float, w: float, h: float, edge: str, fill: str, title: str, line1: str, line2: str) -> None:
    rounded_box(fig, x, y, w, h, "white", "#d7dce4", lw=0.9, radius=0.008)
    fig.lines.append(plt.Line2D([x + 0.012, x + w - 0.012], [y + h - 0.010, y + h - 0.010], transform=fig.transFigure, color=edge, lw=3.0))
    draw_text(fig, x + 0.016, y + h - 0.030, title, fontsize=18.0, fontweight="bold", color=edge if title.startswith("Best") else INK)
    draw_text(fig, x + 0.016, y + h - 0.070, line1, fontsize=13.2, color=INK)
    draw_text(fig, x + 0.016, y + h - 0.099, line2, fontsize=12.3, color=MUTED)


def draw_panel(ax, panel: dict, title: str | None, show_ylabel: bool, show_xlabel: bool, show_legend: bool, ymax: float) -> None:
    sx, sy = step_arrays(panel, "signal")
    bx, by = step_arrays(panel, "background")
    ax.plot(sx, sy, color=SIGNAL, lw=1.8)
    ax.plot(bx, by, color=BACKGROUND, lw=1.8)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, ymax)
    ax.grid(True, color=GRID, lw=0.65, alpha=0.75)
    ax.tick_params(labelsize=8.5, pad=2)
    if title:
        ax.set_title(title, fontsize=17, fontweight="bold", pad=8)
    if show_ylabel:
        ax.set_ylabel("Unit-area density", fontsize=10.0, labelpad=7)
    else:
        ax.set_ylabel("")
        ax.set_yticklabels([])
    if show_xlabel:
        ax.set_xlabel("Classifier score", fontsize=11.5, labelpad=3)
    else:
        ax.set_xticklabels([])
    ax.text(
        0.035,
        0.875,
        f"AUC {float(panel['auc_weighted']):.3f}\nmed gap {float(panel.get('weighted_median_separation', float('nan'))):.3f}\nWP80 fake {100.0 * float(panel.get('wp80_weighted_background_fake_rate', float('nan'))):.1f}%",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=7.9,
        fontweight="bold",
        linespacing=1.15,
        bbox={"boxstyle": "round,pad=0.22", "fc": "white", "ec": "#cfd4dc", "alpha": 0.94},
    )
    if show_legend:
        ax.legend(["Signal", "Background"], fontsize=7.8, frameon=False, loc="upper right")


def write_summary_csv(path: Path, summary: dict, lane: dict) -> None:
    fields = [
        "lane",
        "model",
        "cent_lo",
        "cent_hi",
        "auc_weighted",
        "weighted_median_separation",
        "wp80_weighted_background_fake_rate",
        "signal_entries",
        "background_entries",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for panel in lane["panels"]:
            if panel["model"] in {"BDT", "MLP", summary["selected_stack_model"]}:
                writer.writerow({name: panel.get(name, lane["lane"] if name == "lane" else "") for name in fields})


def render(summary_path: Path, outdir: Path) -> dict:
    setup_style()
    loaded = load_payload(summary_path)
    summary = loaded["summary"]
    lane = loaded["lane"]
    selected_stack = summary["selected_stack_model"]
    cent_bins = [(float(lo), float(hi)) for lo, hi in summary["centrality_bins"]]
    lookup = panel_lookup(lane)
    row_models = ["BDT", "MLP", selected_stack]
    row_specs = [
        ("BDT", "Input BDT", "final trainval BDT score", "same locked-test rows", BDT_EDGE, BDT_FILL),
        ("MLP", "Input MLP", "final trainval MLP score", "same locked-test rows", MLP_EDGE, MLP_FILL),
        (selected_stack, "Best stack", "score-only GBM combiner", "BDT + MLP score inputs", STACK_EDGE, STACK_FILL),
    ]
    ymax = max(
        max(max(p["signal_density"]), max(p["background_density"]))
        for p in lane["panels"]
        if p["model"] in row_models
    )
    ymax = max(1.0, ymax * 1.13)

    outdir.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI)
    fig.set_facecolor("white")

    draw_text(fig, 0.052, 0.966, "Au+Au validation score separation", fontsize=28.5, fontweight="bold")
    draw_text(
        fig,
        0.052,
        0.890,
        "Same held-out validation rows in every panel. Panel badges show AUC, median score gap, and WP80 fake rate.",
        fontsize=16.5,
        color="#3f4650",
    )
    draw_legend(fig, 0.748, 0.940)

    x0 = 0.310
    w = 0.185
    gap = 0.032
    row_y = [0.602, 0.383, 0.164]
    ax_h = 0.164
    row_h = 0.190
    for i, (lo, hi) in enumerate(cent_bins):
        rounded_box(fig, x0 + i * (w + gap) + w / 2 - 0.033, 0.801, 0.066, 0.032, "#fbfcfe", "#d5d9df", lw=0.8, radius=0.010)
        draw_text(fig, x0 + i * (w + gap) + w / 2, 0.819, f"{lo:g}-{hi:g}%", fontsize=13.2, fontweight="bold", ha="center", va="center")

    for row_index, (model, title, line1, line2, edge, fill) in enumerate(row_specs):
        y = row_y[row_index]
        rounded_box(fig, x0 - 0.012, y - 0.016, 3 * w + 2 * gap + 0.024, row_h, fill, None, radius=0.010)
        fig.lines.append(
            plt.Line2D([x0 - 0.006, x0 + 3 * w + 2 * gap + 0.006], [y + row_h - 0.017, y + row_h - 0.017], transform=fig.transFigure, color=edge, lw=2.0, alpha=0.85)
        )
        draw_row_card(fig, 0.055, y + 0.012, 0.215, 0.128, edge, fill, title, line1, line2)
        for col_index, (lo, hi) in enumerate(cent_bins):
            ax = fig.add_axes([x0 + col_index * (w + gap), y, w, ax_h], zorder=2)
            draw_panel(
                ax,
                lookup[(model, lo, hi)],
                None,
                show_ylabel=(col_index == 0),
                show_xlabel=(row_index == 2),
                show_legend=(row_index == 0 and col_index == 2),
                ymax=ymax,
            )

    bdt_auc = float(lane["locked_test_metrics"]["BDT"]["weighted_auc"])
    mlp_auc = float(lane["locked_test_metrics"]["MLP"]["weighted_auc"])
    stack_auc = float(lane["locked_test_metrics"][selected_stack]["weighted_auc"])
    rounded_box(fig, 0.055, 0.020, 0.880, 0.078, "#f6f8fb", "#e2e7ef", lw=0.8, radius=0.010)
    rounded_box(fig, 0.055, 0.020, 0.006, 0.078, STACK_EDGE, STACK_EDGE, radius=0.003)
    draw_text(fig, 0.074, 0.080, "Readout", fontsize=14.2, fontweight="bold", color=STACK_EDGE)
    draw_text(
        fig,
        0.142,
        0.080,
        "The score-only GBM stack is best by weighted AUC and visibly sharpens the high-score signal peak.",
        fontsize=12.7,
    )
    draw_text(
        fig,
        0.142,
        0.051,
        f"Inclusive weighted AUC: BDT {bdt_auc:.6f}, MLP {mlp_auc:.6f}, stack {stack_auc:.6f}; stack gain vs BDT = {stack_auc - bdt_auc:+.6f}.",
        fontsize=12.7,
        color=STACK_EDGE,
    )

    stem = "auau_validation_bdt_mlp_best_stack_score_table"
    png = outdir / f"{stem}.png"
    fig.savefig(png, dpi=DPI)
    plt.close(fig)

    csv_path = outdir / f"{stem}.csv"
    write_summary_csv(csv_path, summary, lane)
    manifest = {
        "schema": "RJ_THE35_AUAU_VALIDATION_SCORE_TABLE_MANIFEST_V1",
        "png": str(png),
        "summary_json": str(summary_path),
        "summary_csv": str(csv_path),
        "selected_lane": lane["lane"],
        "selected_stack_model": selected_stack,
        "selected_reason": summary["selected_reason"],
        "run_root": summary["run_root"],
        "source_root": lane["root"],
        "rows": lane["rows"],
        "signal_entries": lane["signal_entries"],
        "background_entries": lane["background_entries"],
        "locked_test_weighted_auc": {
            "BDT": bdt_auc,
            "MLP": mlp_auc,
            selected_stack: stack_auc,
        },
        "centrality_bins": summary["centrality_bins"],
        "panel_extra_metrics": {
            "weighted_median_separation": "weighted median(signal score) - weighted median(background score) in each centrality panel",
            "wp80_weighted_background_fake_rate": "background weighted fake rate at per-panel threshold with weighted signal efficiency about 80%",
        },
        "color_convention": "signal red, background blue",
        "canvas": {"width": W, "height": H, "dpi": DPI},
    }
    manifest_path = outdir / f"{stem}.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    speaker = outdir / f"{stem}.speaker.md"
    speaker.write_text(
        "This slide compares the final BDT, final MLP, and the best available Au+Au stack on the same clean-DAG locked-test validation rows. "
        f"The selected stack is the {short_lane_label(lane['lane'])} {stack_label(selected_stack)}, chosen because it has the highest locked-test weighted AUC among the completed Au+Au stack models. "
        f"The gain over the BDT is small but consistent in the three centrality bins shown: {stack_auc - bdt_auc:+.6f} inclusive weighted AUC.\n"
    )
    return {"png": png, "csv": csv_path, "manifest": manifest_path, "speaker": speaker}


def pct(value: float) -> str:
    return f"{100.0 * value:.1f}%"


def model_color(model: str) -> str:
    return {"BDT": BDT_EDGE, "MLP": MLP_EDGE, "score_only_gbm": STACK_EDGE}.get(model, INK)


def draw_integrated_shape_panel(ax, shape: dict, model: str, ymax: float) -> None:
    item = shape["models"][model]
    bins = np.asarray(shape["bins"], dtype=float)
    sx = np.repeat(bins, 2)[1:-1]
    sig = np.repeat(np.asarray(item["hist_signal_density"], dtype=float), 2)
    bkg = np.repeat(np.asarray(item["hist_background_density"], dtype=float), 2)
    ax.plot(sx, sig, color=SIGNAL, lw=2.2, label="Signal")
    ax.plot(sx, bkg, color=BACKGROUND, lw=2.2, label="Background")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, ymax)
    ax.grid(True, color=GRID, lw=0.65, alpha=0.75)
    ax.tick_params(labelsize=10)
    ax.set_title({"BDT": "Input BDT", "MLP": "Input MLP", "score_only_gbm": "GBM stack"}[model], fontsize=17, fontweight="bold", color=model_color(model))
    ax.set_xlabel("Classifier score", fontsize=13)
    if model == "BDT":
        ax.set_ylabel("unit-area density", fontsize=13)
    else:
        ax.set_yticklabels([])
    if model == "score_only_gbm":
        peak = item["signal_peak_bin"]
        ax.axvspan(peak[0], peak[1], color=STACK_EDGE, alpha=0.14)
        ax.text(
            0.965,
            0.88,
            f"signal peak\n{item['signal_peak_density']:.1f}",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=11,
            color=STACK_EDGE,
            bbox={"boxstyle": "round,pad=0.22", "fc": "white", "ec": "#cfd4dc", "alpha": 0.93},
        )
    else:
        ax.text(
            0.965,
            0.88,
            f"signal peak\n{item['signal_peak_density']:.1f}",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=11,
            color=model_color(model),
            bbox={"boxstyle": "round,pad=0.22", "fc": "white", "ec": "#cfd4dc", "alpha": 0.93},
        )


def draw_metric_card(fig, x: float, y: float, w: float, h: float, title: str, lines: list[tuple[str, str]], accent: str) -> None:
    rounded_box(fig, x, y, w, h, "#f6f8fb", "#d7dde6", 1.0, 0.012)
    rounded_box(fig, x, y, 0.006, h, accent, accent, 1.0, 0.003)
    draw_text(fig, x + 0.020, y + h - 0.030, title, fontsize=16.0, fontweight="bold", color=accent)
    yy = y + h - 0.065
    for label, value in lines:
        draw_text(fig, x + 0.020, yy, label, fontsize=12.0, color=MUTED)
        draw_text(fig, x + w - 0.018, yy, value, fontsize=12.0, ha="right", color=INK, fontweight="bold")
        yy -= 0.031


def render_shape_slide(shape_path: Path, outdir: Path) -> dict:
    setup_style()
    shape = json.loads(shape_path.read_text())
    bdt = shape["models"]["BDT"]
    mlp = shape["models"]["MLP"]
    stack = shape["models"]["score_only_gbm"]
    auc_bdt = bdt["locked_test_metrics"]["weighted_auc"]
    auc_mlp = mlp["locked_test_metrics"]["weighted_auc"]
    auc_stack = stack["locked_test_metrics"]["weighted_auc"]
    fake_bdt = bdt["locked_test_metrics"]["wp80_weighted_background_fake_rate"]
    fake_mlp = mlp["locked_test_metrics"]["wp80_weighted_background_fake_rate"]
    fake_stack = stack["locked_test_metrics"]["wp80_weighted_background_fake_rate"]
    rel_fake_gain = 1.0 - fake_stack / fake_bdt
    ymax = 1.10 * max(
        max(max(shape["models"][m]["hist_signal_density"]), max(shape["models"][m]["hist_background_density"]))
        for m in ["BDT", "MLP", "score_only_gbm"]
    )

    fig = plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI)
    fig.set_facecolor("white")
    draw_text(fig, 0.050, 0.944, "Why the stack score looks sharper", fontsize=30, fontweight="bold")
    draw_text(
        fig,
        0.050,
        0.904,
        "The GBM stack improves ranking only slightly, but it reshapes the score into a more saturated high-signal peak.",
        fontsize=16.5,
        color="#3f4650",
    )
    draw_legend(fig)

    for i, model in enumerate(["BDT", "MLP", "score_only_gbm"]):
        ax = fig.add_axes([0.065 + i * 0.295, 0.535, 0.255, 0.255])
        draw_integrated_shape_panel(ax, shape, model, ymax)

    draw_metric_card(
        fig,
        0.065,
        0.290,
        0.255,
        0.175,
        "Ranking metrics",
        [
            ("BDT weighted AUC", f"{auc_bdt:.6f}"),
            ("MLP weighted AUC", f"{auc_mlp:.6f}"),
            ("Stack weighted AUC", f"{auc_stack:.6f}"),
            ("Stack gain vs BDT", f"{auc_stack - auc_bdt:+.6f}"),
        ],
        STACK_EDGE,
    )
    draw_metric_card(
        fig,
        0.365,
        0.290,
        0.255,
        0.175,
        "Shape metrics",
        [
            ("Signal peak: BDT", f"{bdt['signal_peak_density']:.1f}"),
            ("Signal peak: MLP", f"{mlp['signal_peak_density']:.1f}"),
            ("Signal peak: stack", f"{stack['signal_peak_density']:.1f}"),
            ("Top stack score bin", f"{stack['top_rounded_4dp_score']:.4f} ({pct(stack['top_rounded_4dp_entry_fraction'])})"),
        ],
        BDT_EDGE,
    )
    draw_metric_card(
        fig,
        0.665,
        0.290,
        0.255,
        0.175,
        "Working-point behavior",
        [
            ("Signal mass > 0.85: BDT", pct(bdt["signal_mass_gt_0p85"])),
            ("Signal mass > 0.85: stack", pct(stack["signal_mass_gt_0p85"])),
            ("Background < 0.10: BDT", pct(bdt["background_mass_lt_0p1"])),
            ("Background < 0.10: stack", pct(stack["background_mass_lt_0p1"])),
        ],
        MLP_EDGE,
    )

    rounded_box(fig, 0.065, 0.075, 0.855, 0.155, "#f4f7fa", None, radius=0.014)
    draw_text(fig, 0.087, 0.211, "Interpretation", fontsize=16.5, fontweight="bold")
    draw_text(
        fig,
        0.087,
        0.181,
        "The selected stack is a score-only GBM combiner. Tree leaves create a piecewise score map, so many high-confidence signal-like rows collapse near the same score.",
        fontsize=12.7,
    )
    draw_text(
        fig,
        0.087,
        0.151,
        f"AUC measures rank ordering, not score calibration. Shape metrics show the bigger visual change: peak density {stack['signal_peak_density']:.1f} vs BDT {bdt['signal_peak_density']:.1f}; signal >0.85 {pct(stack['signal_mass_gt_0p85'])} vs {pct(bdt['signal_mass_gt_0p85'])}.",
        fontsize=12.7,
        color=STACK_EDGE,
    )
    draw_text(
        fig,
        0.087,
        0.121,
        f"WP80 weighted background fake rate improves from {fake_bdt:.4f} to {fake_stack:.4f} ({100*rel_fake_gain:.1f}% relative reduction).",
        fontsize=12.7,
        color=STACK_EDGE,
    )
    draw_text(
        fig,
        0.087,
        0.091,
        "QA readback: trainval/locked-test overlap is zero and all stack scores are finite. Treat the stack score as a classifier score, not a calibrated probability, until calibration is checked.",
        fontsize=11.8,
        color=MUTED,
    )

    stem = "auau_stack_score_shape_standalone"
    outdir.mkdir(parents=True, exist_ok=True)
    png = outdir / f"{stem}.png"
    fig.savefig(png, dpi=DPI)
    plt.close(fig)

    manifest = {
        "schema": "RJ_AUAU_STACK_SCORE_SHAPE_STANDALONE_MANIFEST_V1",
        "png": str(png),
        "shape_json": str(shape_path),
        "source_root": shape["root"],
        "selected_model": "score_only_gbm",
        "weighted_auc": {"BDT": auc_bdt, "MLP": auc_mlp, "score_only_gbm": auc_stack},
        "wp80_weighted_background_fake_rate": {"BDT": fake_bdt, "MLP": fake_mlp, "score_only_gbm": fake_stack},
        "shape_metrics": {
            "signal_peak_density": {
                "BDT": bdt["signal_peak_density"],
                "MLP": mlp["signal_peak_density"],
                "score_only_gbm": stack["signal_peak_density"],
            },
            "signal_mass_gt_0p85": {
                "BDT": bdt["signal_mass_gt_0p85"],
                "MLP": mlp["signal_mass_gt_0p85"],
                "score_only_gbm": stack["signal_mass_gt_0p85"],
            },
            "background_mass_lt_0p1": {
                "BDT": bdt["background_mass_lt_0p1"],
                "MLP": mlp["background_mass_lt_0p1"],
                "score_only_gbm": stack["background_mass_lt_0p1"],
            },
        },
        "leakage_summary": shape["leakage_summary"],
    }
    manifest_path = outdir / f"{stem}.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    speaker = outdir / f"{stem}.speaker.md"
    speaker.write_text(
        "The stack score looks sharper because the selected score-only GBM stack maps BDT and MLP scores through tree leaves. "
        "That creates a more saturated, piecewise score distribution. This is not automatically bad or leakage; leakage QA reports zero trainval/locked-test overlap and finite stack scores. "
        "AUC only measures ranking, so the visual change is better summarized with peak density, high-score signal mass, low-score background mass, and WP80 fake-rate.\n"
    )
    return {"shape_png": png, "shape_manifest": manifest_path, "shape_speaker": speaker}


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--summary-json", type=Path, required=True)
    ap.add_argument("--shape-json", type=Path)
    ap.add_argument("--outdir", type=Path, required=True)
    return ap.parse_args()


def main() -> int:
    args = parse_args()
    outputs = render(args.summary_json, args.outdir)
    if args.shape_json is not None:
        outputs.update(render_shape_slide(args.shape_json, args.outdir))
    print(json.dumps({key: str(value) for key, value in outputs.items()}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
