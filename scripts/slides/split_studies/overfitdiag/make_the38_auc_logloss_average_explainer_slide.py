#!/usr/bin/env python3
"""Render a THE-38 slide explaining AUC non-additivity and logloss additivity."""

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

import argparse
import csv
import json
import math
import re
import textwrap
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch

from slides.common.slide_defaults import SLIDE_DPI, slide_figsize


DEFAULT_ROOT = Path("dataOutput/auauMLDiagnosticRuns/THE8_jet40_centinput_splitstudy_20260528")
DEFAULT_INTERNAL_METRICS = DEFAULT_ROOT / "slideReady/splitstudy_score_separation/the8_jet40_internal_metrics.json"
DEFAULT_REPORT_ROOT = DEFAULT_ROOT / "remote_reports"
DEFAULT_OUTDIR = DEFAULT_ROOT / "slideReady/the38_auc_logloss_average_explainer"
SPLIT_ORDER = ["90/10", "50/50", "10/90"]
SPLIT_TO_REPORT_GLOB = {
    "90/10": "*train90_test10_quickstat_20260528_1450",
    "50/50": "*train50_test50_quickstat_20260528_1450",
    "10/90": "*train10_test90_quickstat_20260528_1450",
}

COLORS = {
    "ink": "#171717",
    "muted": "#5b6472",
    "line": "#d7dde7",
    "grid": "#e7ebf2",
    "panel": "#f7f9fc",
    "panel_edge": "#cbd6e4",
    "train": "#d97706",
    "holdout": "#2563ad",
    "full": "#2f855a",
    "auc_gap": "#3274b9",
    "logloss_gap": "#c94f43",
    "gold": "#b7791f",
    "green_bg": "#edf6f0",
    "blue_bg": "#edf2fb",
    "yellow_bg": "#fff5d6",
    "warn_bg": "#fff2f0",
    "auc_metric": "#5d4b8a",
    "auc_metric_bg": "#f5f2fb",
    "auc_metric_edge": "#c9bddf",
    "logloss_metric": "#006f7a",
    "logloss_metric_bg": "#eef8f8",
    "logloss_metric_edge": "#a6d2d5",
    "split_metric": "#6f5b2a",
}


@dataclass(frozen=True)
class SplitMetric:
    split: str
    train_fraction: float
    train_rows: int
    holdout_rows: int
    full_rows: int
    train_auc: float
    holdout_auc: float
    full_auc: float
    train_logloss: float
    holdout_logloss: float
    weighted_logloss: float
    auc_gap: float
    logloss_gap: float


def parse_summary_value(summary_text: str, key: str) -> str:
    match = re.search(rf"^{re.escape(key)}=(.+)$", summary_text, flags=re.MULTILINE)
    if not match:
        raise SystemExit(f"Missing {key}= in validation_summary.txt")
    return match.group(1).strip()


def parse_float(value: Any, label: str) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError) as exc:
        raise SystemExit(f"Non-numeric {label}: {value!r}") from exc
    if not math.isfinite(out):
        raise SystemExit(f"Non-finite {label}: {value!r}")
    return out


def validation_summary_for_split(report_root: Path, split: str) -> Path:
    matches = sorted(report_root.glob(SPLIT_TO_REPORT_GLOB[split] + "/validation_summary.txt"))
    if len(matches) != 1:
        raise SystemExit(f"Could not locate one validation_summary for {split}: {matches}")
    return matches[0]


def load_metrics(internal_metrics: Path, report_root: Path) -> list[SplitMetric]:
    payload = json.loads(internal_metrics.read_text())
    out: list[SplitMetric] = []
    for split in SPLIT_ORDER:
        item = payload[split]
        summary_path = validation_summary_for_split(report_root, split)
        summary_text = summary_path.read_text()
        train_rows = int(item["train_rows"])
        holdout_rows = int(item["test_rows"])
        full_rows = int(parse_summary_value(summary_text, "scored_entries"))
        train_logloss = parse_float(item["train_logloss"], f"{split} train_logloss")
        holdout_logloss = parse_float(item["holdout_logloss"], f"{split} holdout_logloss")
        weighted_logloss = (train_rows * train_logloss + holdout_rows * holdout_logloss) / (train_rows + holdout_rows)
        train_fraction = train_rows / (train_rows + holdout_rows)
        out.append(
            SplitMetric(
                split=split,
                train_fraction=train_fraction,
                train_rows=train_rows,
                holdout_rows=holdout_rows,
                full_rows=full_rows,
                train_auc=parse_float(item["train_auc"], f"{split} train_auc"),
                holdout_auc=parse_float(item["auc"], f"{split} holdout_auc"),
                full_auc=parse_float(parse_summary_value(summary_text, "centInput_pt1535_auc"), f"{split} full_auc"),
                train_logloss=train_logloss,
                holdout_logloss=holdout_logloss,
                weighted_logloss=weighted_logloss,
                auc_gap=parse_float(item["train_auc"], f"{split} train_auc") - parse_float(item["auc"], f"{split} holdout_auc"),
                logloss_gap=holdout_logloss - train_logloss,
            )
        )
        if full_rows != 20190290:
            raise SystemExit(f"{split} full scored rows changed unexpectedly: {full_rows}")
    return out


def add_card(fig, xywh, *, face=COLORS["panel"], edge=COLORS["panel_edge"], radius=0.014, lw=1.1):
    x, y, w, h = xywh
    fig.patches.append(
        FancyBboxPatch(
            (x, y),
            w,
            h,
            transform=fig.transFigure,
            boxstyle=f"round,pad=0.010,rounding_size={radius}",
            facecolor=face,
            edgecolor=edge,
            linewidth=lw,
            zorder=-10,
        )
    )


def style_axis(ax, *, ygrid=True):
    if ygrid:
        ax.grid(axis="y", color=COLORS["grid"], linewidth=1.0)
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color("#9aa7b8")
    ax.tick_params(labelsize=12.0, colors=COLORS["muted"])
    ax.xaxis.label.set_color(COLORS["ink"])
    ax.yaxis.label.set_color(COLORS["ink"])


def plot_panel(fig, xywh, title: str, *, face=COLORS["panel"], edge=COLORS["panel_edge"], title_color=COLORS["ink"]):
    add_card(fig, xywh, face=face, edge=edge)
    x, y, w, h = xywh
    fig.text(x + 0.020, y + h - 0.026, title, fontsize=15.9, fontweight="bold", color=title_color, va="top")
    return fig.add_axes([x + 0.052, y + 0.070, w - 0.080, h - 0.135])


def claim_card(fig, xywh, number: str, title: str, body: str, *, face: str, edge: str, accent: str):
    add_card(fig, xywh, face=face, edge=edge, radius=0.011)
    x, y, w, h = xywh
    fig.text(
        x + 0.018,
        y + h - 0.018,
        number,
        fontsize=18.0,
        fontweight="bold",
        color=accent,
        va="top",
        ha="left",
    )
    fig.text(
        x + 0.052,
        y + h - 0.016,
        title,
        fontsize=13.6,
        fontweight="bold",
        color=COLORS["ink"],
        va="top",
        ha="left",
    )
    fig.text(
        x + 0.052,
        y + h - 0.046,
        "\n".join(textwrap.wrap(body, width=37)),
        fontsize=11.35,
        color=COLORS["ink"],
        va="top",
        ha="left",
        linespacing=1.10,
    )


def rule_card(fig, xywh, lead: str, body: str, *, face: str, edge: str, accent: str):
    add_card(fig, xywh, face=face, edge=edge, radius=0.010)
    x, y, w, h = xywh
    fig.text(x + 0.018, y + h - 0.020, lead, fontsize=13.5, fontweight="bold", color=accent, va="top", ha="left")
    fig.text(
        x + 0.140,
        y + h - 0.020,
        "\n".join(textwrap.wrap(body, width=52)),
        fontsize=12.0,
        color=COLORS["ink"],
        va="top",
        ha="left",
        linespacing=1.08,
    )


def same_green_card(fig, xywh):
    add_card(fig, xywh, face=COLORS["green_bg"], edge="#9fceb0", radius=0.010, lw=1.25)
    x, y, w, h = xywh
    ax = fig.add_axes([x + 0.022, y + 0.024, 0.055, h - 0.048])
    ax.plot([0.12, 0.88], [0.5, 0.5], color=COLORS["full"], marker="o", linewidth=3.0, markersize=7.5)
    ax.set_axis_off()
    fig.text(
        x + 0.088,
        y + h - 0.021,
        "Same green curve in both plots = full scored sample",
        fontsize=15.2,
        fontweight="bold",
        color=COLORS["ink"],
        va="top",
    )
    fig.text(
        x + 0.088,
        y + h - 0.052,
        "The evaluation target is the same; the difference is how each metric aggregates train and holdout rows.",
        fontsize=12.3,
        color=COLORS["ink"],
        va="top",
    )
    fig.text(x + w - 0.265, y + h - 0.025, "AUC ranks pairs", fontsize=12.4, fontweight="bold", color=COLORS["auc_metric"], ha="left", va="top")
    fig.text(x + w - 0.148, y + h - 0.025, "|", fontsize=12.4, fontweight="bold", color=COLORS["muted"], ha="center", va="top")
    fig.text(x + w - 0.136, y + h - 0.025, "logloss sums rows", fontsize=12.4, fontweight="bold", color=COLORS["logloss_metric"], ha="left", va="top")


def plot_line(ax, x, values, *, color, label, marker="o", linestyle="-", lw=2.6):
    ax.plot(x, values, color=color, label=label, marker=marker, markersize=7.0, linewidth=lw, linestyle=linestyle)


def add_math_card(fig, xywh, title: str, body_lines: list[tuple[str, str]], face: str, edge: str):
    add_card(fig, xywh, face=face, edge=edge, radius=0.012)
    x, y, w, h = xywh
    fig.text(x + 0.018, y + h - 0.028, title, fontsize=16.0, fontweight="bold", color=COLORS["ink"], va="top")
    yy = y + h - 0.070
    for lead, body in body_lines:
        fig.text(x + 0.018, yy, lead, fontsize=12.4, fontweight="bold", color=COLORS["ink"], va="top")
        fig.text(
            x + 0.105,
            yy,
            "\n".join(textwrap.wrap(body, width=44)),
            fontsize=12.2,
            color=COLORS["ink"],
            va="top",
            linespacing=1.15,
        )
        yy -= 0.060 if len(body) < 75 else 0.076


def render_slide(metrics: list[SplitMetric], args: argparse.Namespace) -> Path:
    out_png = args.outdir / "the38_auc_logloss_average_explainer.png"
    x = np.arange(len(metrics))
    labels = [f"{m.split}\n{m.train_rows / 1e6:.2f}M train rows" for m in metrics]
    train_auc = np.asarray([m.train_auc for m in metrics])
    holdout_auc = np.asarray([m.holdout_auc for m in metrics])
    full_auc = np.asarray([m.full_auc for m in metrics])
    train_logloss = np.asarray([m.train_logloss for m in metrics])
    holdout_logloss = np.asarray([m.holdout_logloss for m in metrics])
    full_logloss = np.asarray([m.weighted_logloss for m in metrics])
    auc_gap = np.asarray([m.auc_gap for m in metrics])
    logloss_gap = np.asarray([m.logloss_gap for m in metrics])

    plt.rcParams.update(
        {
            "font.family": "Times New Roman",
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "axes.facecolor": "white",
            "axes.titleweight": "bold",
        }
    )
    fig = plt.figure(figsize=slide_figsize(SLIDE_DPI), dpi=SLIDE_DPI)
    fig.text(0.050, 0.942, "Same full-sample check, different metric algebra", fontsize=27.5, fontweight="bold", color=COLORS["ink"], va="center")
    fig.text(
        0.050,
        0.902,
        "Depth-4 split study: green is the full scored sample in both panels; only logloss is additive over rows.",
        fontsize=15.0,
        color=COLORS["muted"],
        va="center",
    )
    fig.text(
        0.948,
        0.925,
        "Jet40-inclusive Branch A\ncentInput_pt1535",
        fontsize=12.5,
        color=COLORS["muted"],
        ha="right",
        va="top",
        linespacing=1.15,
    )

    same_green_card(fig, (0.046, 0.790, 0.909, 0.080))

    ax_auc = plot_panel(
        fig,
        (0.046, 0.458, 0.440, 0.292),
        "AUC: pairwise ranking, not additive",
        face=COLORS["auc_metric_bg"],
        edge=COLORS["auc_metric_edge"],
        title_color=COLORS["auc_metric"],
    )
    ax_loss = plot_panel(
        fig,
        (0.535, 0.458, 0.420, 0.292),
        "Logloss: row-wise loss, additive",
        face=COLORS["logloss_metric_bg"],
        edge=COLORS["logloss_metric_edge"],
        title_color=COLORS["logloss_metric"],
    )
    ax_gap = plot_panel(fig, (0.046, 0.138, 0.486, 0.252), "Both gaps grow as training rows shrink")

    plot_line(ax_auc, x, train_auc, color=COLORS["train"], label="Train")
    plot_line(ax_auc, x, holdout_auc, color=COLORS["holdout"], label="Holdout")
    plot_line(ax_auc, x, full_auc, color=COLORS["full"], label="Full scored sample")
    ax_auc.set_xticks(x, labels)
    ax_auc.set_ylabel("AUC  (higher is better)", fontsize=12.2)
    ax_auc.set_ylim(0.8690, 0.8794)
    ax_auc.set_xlim(-0.20, 2.26)
    ax_auc.legend(frameon=False, fontsize=11.0, loc="upper center", bbox_to_anchor=(0.55, 1.15), ncol=3, handlelength=1.7, columnspacing=0.85)
    style_axis(ax_auc)

    plot_line(ax_loss, x, train_logloss, color=COLORS["train"], label="Train")
    plot_line(ax_loss, x, holdout_logloss, color=COLORS["holdout"], label="Holdout")
    plot_line(ax_loss, x, full_logloss, color=COLORS["full"], label="Full scored sample", marker="o", lw=2.4)
    ax_loss.set_xticks(x, labels)
    ax_loss.set_ylabel("Logloss  (lower is better)", fontsize=12.2)
    ax_loss.set_ylim(0.4210, 0.4304)
    ax_loss.set_xlim(-0.20, 2.26)
    ax_loss.legend(frameon=False, fontsize=10.6, loc="upper center", bbox_to_anchor=(0.52, 1.15), ncol=3, handlelength=1.5, columnspacing=0.70)
    style_axis(ax_loss)

    width = 0.33
    ax_gap.bar(x - width / 2, auc_gap, width=width, color=COLORS["auc_gap"], label="AUC gap")
    ax_gap.bar(x + width / 2, logloss_gap, width=width, color=COLORS["logloss_gap"], label="Logloss gap")
    ax_gap.axhline(0, color="#9aa7b8", linewidth=1.0)
    ax_gap.set_xticks(x, [f"{m.split}\n{m.train_rows / 1e6:.2f}M train" for m in metrics])
    ax_gap.set_ylabel("Gap size", fontsize=12.2)
    ax_gap.set_ylim(-0.0016, 0.0078)
    ax_gap.legend(frameon=False, fontsize=10.8, loc="upper center", bbox_to_anchor=(0.58, 1.19), ncol=2, handlelength=1.7, columnspacing=1.00)
    style_axis(ax_gap)
    for xx, yy in zip(x - width / 2, auc_gap):
        y_text = 0.00035 if abs(yy) < 0.00025 else yy + 0.00020
        ax_gap.text(xx, y_text, f"{yy:.4f}", ha="center", va="bottom", fontsize=10.0, color=COLORS["auc_gap"])
    for xx, yy in zip(x + width / 2, logloss_gap):
        if yy < 0:
            ax_gap.text(xx, -0.00062, f"{yy:.4f}", ha="center", va="center", fontsize=10.0, color=COLORS["logloss_gap"])
        else:
            ax_gap.text(xx, yy + 0.00018, f"{yy:.4f}", ha="center", va="bottom", fontsize=10.0, color=COLORS["logloss_gap"])

    add_card(fig, (0.570, 0.138, 0.385, 0.252), face="#fbfcfe", edge=COLORS["panel_edge"], radius=0.012)
    fig.text(0.762, 0.356, "How to read the green curve", fontsize=16.2, fontweight="bold", color=COLORS["ink"], ha="center", va="top")
    readout_rows = [
        ("AUC", COLORS["auc_metric"], "Full AUC is not a weighted average;\ncross-pairs enter the ranking."),
        ("Logloss", COLORS["logloss_metric"], "Full logloss is the weighted train+holdout loss;\nrow losses add."),
        ("Split", COLORS["split_metric"], "Use the gap bars: both gaps open\nas training rows shrink."),
    ]
    for yy, (lead, color, body) in zip([0.304, 0.246, 0.188], readout_rows):
        fig.text(0.640, yy, lead, fontsize=12.9, fontweight="bold", color=color, ha="right", va="top")
        fig.text(
            0.665,
            yy,
            body,
            fontsize=12.0,
            color=COLORS["ink"],
            ha="left",
            va="top",
            linespacing=1.08,
        )

    add_card(fig, (0.046, 0.050, 0.909, 0.054), face="#fff9e8", edge="#d9b45c", radius=0.010)
    fig.text(
        0.065,
        0.078,
        "Takeaway:",
        fontsize=13.7,
        fontweight="bold",
        color=COLORS["ink"],
        va="center",
    )
    fig.text(
        0.143,
        0.078,
        "the widening train-holdout gaps, not the green full-sample curve alone, are what flag the reduced-training overfit risk.",
        fontsize=13.05,
        color=COLORS["ink"],
        va="center",
    )

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=SLIDE_DPI)
    plt.close(fig)
    return out_png


def metric_rows(metrics: list[SplitMetric]) -> list[dict[str, Any]]:
    return [
        {
            "split": m.split,
            "train_fraction": m.train_fraction,
            "train_rows": m.train_rows,
            "holdout_rows": m.holdout_rows,
            "full_scored_rows": m.full_rows,
            "train_auc": m.train_auc,
            "holdout_auc": m.holdout_auc,
            "full_scored_auc": m.full_auc,
            "train_logloss": m.train_logloss,
            "holdout_logloss": m.holdout_logloss,
            "full_scored_logloss": m.weighted_logloss,
            "full_scored_logloss_source": "weighted_train_holdout_identity_from_internal_metrics",
            "weighted_train_holdout_logloss": m.weighted_logloss,
            "auc_gap_train_minus_holdout": m.auc_gap,
            "logloss_gap_holdout_minus_train": m.logloss_gap,
        }
        for m in metrics
    ]


def write_sidecars(metrics: list[SplitMetric], args: argparse.Namespace, out_png: Path) -> tuple[Path, Path, Path, Path]:
    rows = metric_rows(metrics)
    metrics_json = args.outdir / "the38_auc_logloss_average_explainer_metrics.json"
    metrics_csv = args.outdir / "the38_auc_logloss_average_explainer_metrics.csv"
    manifest_json = args.outdir / "the38_auc_logloss_average_explainer_manifest.json"
    script_md = args.outdir / "the38_auc_logloss_average_explainer_script.md"

    metrics_json.write_text(json.dumps({"schema": "THE38_AUC_LOGLOSS_AVERAGE_EXPLAINER_METRICS_V1", "rows": rows}, indent=2, sort_keys=True) + "\n")
    with metrics_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    manifest = {
        "schema": "THE38_AUC_LOGLOSS_AVERAGE_EXPLAINER_MANIFEST_V1",
        "output_png": str(out_png),
        "metrics_json": str(metrics_json),
        "metrics_csv": str(metrics_csv),
        "script_md": str(script_md),
        "internal_metrics": str(args.internal_metrics),
        "report_root": str(args.report_root),
        "render_size_px": [2560, 1440],
        "model": "Jet40-inclusive Branch A centInput_pt1535 depth-4 BDT split study",
        "interpretation_boundary": (
            "AUC full scored sample is an independently computed ranking metric and is not a weighted average. "
            "The plotted green logloss is labeled full scored sample because logloss is row-additive; "
            "the value is derived from the train+holdout weighted-average identity using the available internal diagnostics. "
            "The pulled validation summaries report full scored-sample AUC, but do not contain a separate full-cache logloss field."
        ),
        "full_cache_logloss_direct_field_available": False,
    }
    manifest_json.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")

    script_md.write_text(
        "\n".join(
            [
                "# THE-38 Script - AUC and Logloss Average Check",
                "",
                "This slide is meant to resolve a subtle but important point from the split-study comparison: the green curve means the same thing in both top panels. It is the full scored sample for that split model. What changes is not the sample definition, but the algebra of the metric.",
                "",
                "For AUC, the full-sample value is not a weighted average of the train and holdout AUCs. AUC is a pair-ranking metric: it asks how often a signal row ranks above a background row. When I evaluate the full scored sample, I also create cross-pairs across the train and holdout boundary, so the subset AUCs are not additive.",
                "",
                "For logloss, the situation is different. Logloss is row-wise, so the full-sample loss is just the weighted sum of the row losses. Under matching row and weight definitions, that means full loss equals W train times train loss plus W holdout times holdout loss, divided by the total weight.",
                "",
                "The lower-left panel is the actual overfitting readout. As the training fraction shrinks from ninety percent to ten percent, the train metrics improve because the same model capacity is being fit to fewer rows. But the holdout does not improve, and the train-holdout gaps open.",
                "",
                "So the conclusion is not that 10/90 is better because the training score looks better, and it is not decided just by the green full-sample curve. The widening train-holdout gaps are the warning sign: reduced training statistics make the model look better on its training rows while generalization does not improve.",
                "",
            ]
        ),
        encoding="utf-8",
    )
    return metrics_json, metrics_csv, manifest_json, script_md


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--internal-metrics", type=Path, default=DEFAULT_INTERNAL_METRICS)
    parser.add_argument("--report-root", type=Path, default=DEFAULT_REPORT_ROOT)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    metrics = load_metrics(args.internal_metrics, args.report_root)
    out_png = render_slide(metrics, args)
    metrics_json, metrics_csv, manifest_json, script_md = write_sidecars(metrics, args, out_png)
    print("THE38_AUC_LOGLOSS_AVERAGE_EXPLAINER_READY")
    print(f"png={out_png}")
    print(f"metrics_json={metrics_json}")
    print(f"metrics_csv={metrics_csv}")
    print(f"manifest_json={manifest_json}")
    print(f"script_md={script_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
