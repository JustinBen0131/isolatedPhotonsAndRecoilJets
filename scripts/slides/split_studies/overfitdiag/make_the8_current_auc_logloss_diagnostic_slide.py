#!/usr/bin/env python3
"""Render an enhanced current-data AUC/logloss split diagnostic slide."""

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
import json
import math
import textwrap
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.offsetbox import AnnotationBbox, TextArea, VPacker
from matplotlib.patches import FancyBboxPatch

from slides.common.slide_defaults import SLIDE_DPI, slide_figsize


DEFAULT_OUTDIR = Path(
    "dataOutput/auauTightBDTValidation/THE59_model_comparison_remakes_20260620/"
    "auc_logloss_diagnostic_current_20260622"
)
DEFAULT_METRICS = DEFAULT_OUTDIR / "the8_current_true_full_validation_metrics_20260622.json"

COLORS = {
    "ink": "#121826",
    "muted": "#5f6877",
    "grid": "#e2e8f0",
    "panel": "#f8fafc",
    "panel_edge": "#cbd5e1",
    "train": "#d97706",
    "holdout": "#2563ad",
    "full": "#2f855a",
    "trap": "#8a8f99",
    "auc": "#5d4b8a",
    "auc_bg": "#f6f2fb",
    "auc_edge": "#c9bddf",
    "loss": "#006f7a",
    "loss_bg": "#eef8f8",
    "loss_edge": "#a6d2d5",
    "gap_auc": "#3274b9",
    "gap_loss": "#c94f43",
    "yellow_bg": "#fff8df",
    "yellow_edge": "#d9b45c",
    "green_bg": "#edf7f1",
    "green_edge": "#9fceb0",
    "blue_bg": "#edf3fb",
    "blue_edge": "#adc4e7",
    "red_bg": "#fff2ef",
    "red_edge": "#e5aaa0",
}


@dataclass(frozen=True)
class SplitMetric:
    label: str
    train_rows: int
    holdout_rows: int
    train_auc: float
    holdout_auc: float
    full_auc: float
    train_logloss: float
    holdout_logloss: float
    full_logloss: float

    @property
    def auc_gap(self) -> float:
        return self.train_auc - self.holdout_auc

    @property
    def logloss_gap(self) -> float:
        return self.holdout_logloss - self.train_logloss


def parse_float(value: Any, label: str) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError) as exc:
        raise SystemExit(f"Non-numeric {label}: {value!r}") from exc
    if not math.isfinite(out):
        raise SystemExit(f"Non-finite {label}: {value!r}")
    return out


def load_metrics(path: Path) -> tuple[list[SplitMetric], dict[str, Any]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    rows = payload.get("rows")
    if not isinstance(rows, list) or not rows:
        raise SystemExit(f"{path} must contain non-empty rows[]")
    metrics: list[SplitMetric] = []
    for row in rows:
        label = str(row["split"])
        metrics.append(
            SplitMetric(
                label=label,
                train_rows=int(row["train_rows"]),
                holdout_rows=int(row["holdout_rows"]),
                train_auc=parse_float(row["train_auc"], f"{label} train_auc"),
                holdout_auc=parse_float(row["holdout_auc"], f"{label} holdout_auc"),
                full_auc=parse_float(row["full_auc"], f"{label} full_auc"),
                train_logloss=parse_float(row["train_logloss"], f"{label} train_logloss"),
                holdout_logloss=parse_float(row["holdout_logloss"], f"{label} holdout_logloss"),
                full_logloss=parse_float(row["full_logloss"], f"{label} full_logloss"),
            )
        )
    order = {"90/10": 0, "50/50": 1, "10/90": 2}
    metrics.sort(key=lambda m: order.get(m.label, 999))
    return metrics, payload


def add_card(fig, xywh, *, face=COLORS["panel"], edge=COLORS["panel_edge"], radius=0.014, lw=1.15):
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


def add_panel(
    fig,
    xywh,
    title: str,
    *,
    face: str,
    edge: str,
    color: str,
    title_size: float = 20.5,
    legend_column: bool = False,
):
    add_card(fig, xywh, face=face, edge=edge, radius=0.014, lw=1.25)
    x, y, w, h = xywh
    fig.text(x + 0.020, y + h - 0.030, title, fontsize=title_size, fontweight="bold", color=color, va="top")
    right_margin = 0.200 if legend_column else 0.095
    return fig.add_axes([x + 0.060, y + 0.078, w - right_margin, h - 0.152])


def series_legend(fig, x: float, y: float, items: list[tuple[str, str, str]], *, fontsize: float = 12.7, dy: float = 0.034) -> None:
    for i, (label, color, style) in enumerate(items):
        yy = y - i * dy
        linestyle = "--" if style == "dash" else "-"
        marker = "o" if style != "bar" else "s"
        fig.lines.append(
            Line2D(
                [x, x + 0.020],
                [yy, yy],
                transform=fig.transFigure,
                color=color,
                linewidth=3.0,
                marker=marker,
                markersize=7.0,
                linestyle=linestyle,
                clip_on=False,
                zorder=20,
            )
        )
        fig.text(x + 0.026, yy, label, fontsize=fontsize, color=COLORS["ink"], va="center", ha="left")


def style_axis(ax, *, xlabel: bool = True, tick_size: float = 14.0):
    ax.grid(axis="y", color=COLORS["grid"], linewidth=1.0)
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color("#9aa7b8")
        ax.spines[spine].set_linewidth(1.0)
    ax.tick_params(labelsize=tick_size, colors=COLORS["muted"], pad=5)
    if not xlabel:
        ax.tick_params(labelbottom=False)
    ax.xaxis.label.set_color(COLORS["ink"])
    ax.yaxis.label.set_color(COLORS["ink"])


def fmt_rows(n: int) -> str:
    return f"{n / 1e6:.2f}M"


def metric_rows(metrics: list[SplitMetric]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for metric in metrics:
        rows.append(
            {
                "split": metric.label,
                "train_rows": metric.train_rows,
                "holdout_rows": metric.holdout_rows,
                "full_validation_rows": metric.train_rows + metric.holdout_rows,
                "train_auc": metric.train_auc,
                "holdout_auc": metric.holdout_auc,
                "full_validation_auc": metric.full_auc,
                "train_logloss": metric.train_logloss,
                "holdout_logloss": metric.holdout_logloss,
                "full_validation_logloss": metric.full_logloss,
                "auc_gap_train_minus_holdout": metric.auc_gap,
                "logloss_gap_holdout_minus_train": metric.logloss_gap,
            }
        )
    return rows


def wrap_text(text: str, width: int) -> str:
    return "\n".join(textwrap.wrap(text, width=width, break_long_words=False))


def callout(fig, xywh, heading: str, body: str, *, face: str, edge: str, accent: str, wrap: int = 38):
    add_card(fig, xywh, face=face, edge=edge, radius=0.012, lw=1.15)
    x, y, w, h = xywh
    fig.text(x + 0.020, y + h - 0.026, heading, fontsize=16.4, fontweight="bold", color=accent, va="top")
    fig.text(
        x + 0.020,
        y + h - 0.066,
        wrap_text(body, wrap),
        fontsize=12.6,
        color=COLORS["ink"],
        va="top",
        linespacing=1.12,
    )


def centered_callout(
    fig,
    xywh,
    heading: str,
    body: str,
    *,
    face: str,
    edge: str,
    accent: str,
    wrap: int = 34,
):
    add_card(fig, xywh, face=face, edge=edge, radius=0.012, lw=1.15)
    x, y, w, h = xywh
    heading_area = TextArea(
        heading,
        textprops={
            "fontsize": 16.0,
            "fontweight": "bold",
            "color": accent,
            "ha": "center",
            "multialignment": "center",
        },
    )
    body_area = TextArea(
        wrap_text(body, wrap),
        textprops={
            "fontsize": 12.4,
            "color": COLORS["ink"],
            "ha": "center",
            "multialignment": "center",
            "linespacing": 1.12,
        },
    )
    packed = VPacker(children=[heading_area, body_area], align="center", pad=0, sep=5)
    fig.add_artist(
        AnnotationBbox(
            packed,
            (x + w / 2, y + h / 2),
            xycoords=fig.transFigure,
            frameon=False,
            box_alignment=(0.5, 0.5),
            pad=0,
            zorder=30,
        )
    )


def render_slide(metrics: list[SplitMetric], payload: dict[str, Any], args: argparse.Namespace) -> Path:
    plt.rcParams.update(
        {
            "font.family": "Times New Roman",
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "axes.facecolor": "white",
            "axes.titleweight": "bold",
        }
    )
    out_png = args.outdir / args.output_name
    x = np.arange(len(metrics))
    labels = [m.label for m in metrics]
    train_auc = np.asarray([m.train_auc for m in metrics])
    holdout_auc = np.asarray([m.holdout_auc for m in metrics])
    full_auc = np.asarray([m.full_auc for m in metrics])
    train_loss = np.asarray([m.train_logloss for m in metrics])
    holdout_loss = np.asarray([m.holdout_logloss for m in metrics])
    full_loss = np.asarray([m.full_logloss for m in metrics])
    auc_gap = np.asarray([m.auc_gap for m in metrics])
    loss_gap = np.asarray([m.logloss_gap for m in metrics])

    fig = plt.figure(figsize=slide_figsize(SLIDE_DPI), dpi=SLIDE_DPI)
    fig.text(
        0.046,
        0.938,
        "AUC and logloss split diagnostics",
        fontsize=33.5,
        fontweight="bold",
        color=COLORS["ink"],
        va="center",
    )

    add_card(fig, (0.047, 0.826, 0.906, 0.068), face=COLORS["green_bg"], edge=COLORS["green_edge"], radius=0.012, lw=1.20)
    fig.text(
        0.066,
        0.861,
        "Current corrected Au+Au",
        fontsize=16.0,
        fontweight="bold",
        color=COLORS["ink"],
        va="center",
    )
    fig.text(
        0.315,
        0.861,
        "Full validation = train+holdout rows scored together",
        fontsize=14.2,
        color=COLORS["ink"],
        va="center",
    )
    fig.text(0.680, 0.861, "AUC pair-rank", fontsize=14.6, fontweight="bold", color=COLORS["auc"], va="center")
    fig.text(0.785, 0.861, "|", fontsize=14.6, color=COLORS["muted"], va="center")
    fig.text(0.805, 0.861, "logloss row-loss", fontsize=14.6, fontweight="bold", color=COLORS["loss"], va="center")

    auc_panel = (0.047, 0.480, 0.443, 0.300)
    loss_panel = (0.537, 0.480, 0.416, 0.300)
    gap_panel = (0.047, 0.176, 0.906, 0.224)
    ax_auc = add_panel(
        fig,
        auc_panel,
        "AUC: train, holdout, full validation",
        face=COLORS["auc_bg"],
        edge=COLORS["auc_edge"],
        color=COLORS["auc"],
        title_size=19.0,
        legend_column=True,
    )
    ax_loss = add_panel(
        fig,
        loss_panel,
        "Logloss: train, holdout, full validation",
        face=COLORS["loss_bg"],
        edge=COLORS["loss_edge"],
        color=COLORS["loss"],
        title_size=19.0,
        legend_column=True,
    )
    ax_gap = add_panel(
        fig,
        gap_panel,
        "Overfit diagnostic: the 10/90 lane opens the largest gap",
        face=COLORS["panel"],
        edge=COLORS["panel_edge"],
        color=COLORS["ink"],
        title_size=20.0,
    )

    ax_auc.plot(x, train_auc, color=COLORS["train"], marker="o", linewidth=3.1, markersize=8.2, label="Train")
    ax_auc.plot(x, holdout_auc, color=COLORS["holdout"], marker="o", linewidth=3.1, markersize=8.2, label="Holdout")
    ax_auc.plot(
        x,
        full_auc,
        color=COLORS["full"],
        marker="o",
        linewidth=2.6,
        markersize=7.4,
        linestyle="--",
        label="Full validation",
    )
    ax_auc.set_xticks(x, labels)
    ax_auc.set_ylabel("AUC  (higher is better)", fontsize=16.0)
    ax_auc.set_ylim(0.8580, 0.8682)
    ax_auc.set_xlim(-0.18, len(metrics) - 0.82)
    style_axis(ax_auc, tick_size=15.0)
    series_legend(
        fig,
        auc_panel[0] + auc_panel[2] - 0.135,
        auc_panel[1] + 0.198,
        [
            ("Train", COLORS["train"], "solid"),
            ("Holdout", COLORS["holdout"], "solid"),
            ("Full validation", COLORS["full"], "dash"),
        ],
        fontsize=11.6,
        dy=0.031,
    )
    fig.text(
        auc_panel[0] + auc_panel[2] - 0.134,
        auc_panel[1] + 0.088,
        "green is scored on\nall split rows",
        fontsize=10.8,
        color=COLORS["auc"],
        va="top",
        ha="left",
        linespacing=1.05,
    )

    ax_loss.plot(x, train_loss, color=COLORS["train"], marker="o", linewidth=3.1, markersize=8.2, label="Train")
    ax_loss.plot(x, holdout_loss, color=COLORS["holdout"], marker="o", linewidth=3.1, markersize=8.2, label="Holdout")
    ax_loss.plot(
        x,
        full_loss,
        color=COLORS["full"],
        marker="o",
        linewidth=3.1,
        markersize=8.2,
        linestyle="--",
        label="Full validation",
    )
    ax_loss.set_xticks(x, labels)
    ax_loss.set_ylabel("Logloss  (lower is better)", fontsize=16.0)
    ax_loss.set_ylim(0.4350, 0.4442)
    ax_loss.set_xlim(-0.18, len(metrics) - 0.82)
    style_axis(ax_loss, tick_size=15.0)
    series_legend(
        fig,
        loss_panel[0] + loss_panel[2] - 0.128,
        loss_panel[1] + 0.198,
        [
            ("Train", COLORS["train"], "solid"),
            ("Holdout", COLORS["holdout"], "solid"),
            ("Full validation", COLORS["full"], "dash"),
        ],
        fontsize=11.5,
        dy=0.031,
    )
    fig.text(
        loss_panel[0] + loss_panel[2] - 0.128,
        loss_panel[1] + 0.088,
        "green is scored on\nall split rows",
        fontsize=10.8,
        color=COLORS["loss"],
        va="top",
        ha="left",
        linespacing=1.05,
    )

    width = 0.30
    ax_gap.bar(x - width / 2, auc_gap, width=width, color=COLORS["gap_auc"], label="AUC gap: train - holdout")
    ax_gap.bar(x + width / 2, loss_gap, width=width, color=COLORS["gap_loss"], label="Logloss gap: holdout - train")
    ax_gap.axhline(0, color="#9aa7b8", linewidth=1.1)
    ax_gap.set_xticks(x, [f"{m.label}\ntrain {fmt_rows(m.train_rows)}" for m in metrics])
    ax_gap.set_ylabel("Gap size", fontsize=15.6)
    gap_max = max(float(np.max(np.abs(auc_gap))), float(np.max(np.abs(loss_gap))))
    ax_gap.set_ylim(-0.0030, max(0.0092, gap_max * 1.72))
    fig.text(
        gap_panel[0] + gap_panel[2] - 0.236,
        gap_panel[1] + gap_panel[3] - 0.046,
        "AUC gap",
        fontsize=12.8,
        color=COLORS["gap_auc"],
        fontweight="bold",
        va="center",
        ha="left",
    )
    fig.text(
        gap_panel[0] + gap_panel[2] - 0.145,
        gap_panel[1] + gap_panel[3] - 0.046,
        "Logloss gap",
        fontsize=12.8,
        color=COLORS["gap_loss"],
        fontweight="bold",
        va="center",
        ha="left",
    )
    style_axis(ax_gap, tick_size=14.0)
    for xx, yy in zip(x - width / 2, auc_gap):
        if yy >= 0:
            y_text = yy + 0.00018
            va = "bottom"
        else:
            y_text = min(yy - 0.00020, -0.00100)
            va = "top"
        ax_gap.text(xx, y_text, f"{yy:+.4f}", ha="center", va=va, fontsize=11.8, color=COLORS["gap_auc"], fontweight="bold")
    for xx, yy in zip(x + width / 2, loss_gap):
        if yy >= 0:
            y_text = yy + 0.00018
            va = "bottom"
        else:
            y_text = min(yy - 0.00020, -0.00100)
            va = "top"
        ax_gap.text(xx, y_text, f"{yy:+.4f}", ha="center", va=va, fontsize=11.8, color=COLORS["gap_loss"], fontweight="bold")

    centered_callout(
        fig,
        (0.047, 0.038, 0.284, 0.090),
        "Metric rule",
        "AUC ranks score pairs; full AUC is scored directly.",
        face=COLORS["auc_bg"],
        edge=COLORS["auc_edge"],
        accent=COLORS["auc"],
        wrap=34,
    )
    centered_callout(
        fig,
        (0.358, 0.038, 0.284, 0.090),
        "Additive rule",
        "Logloss is row-wise; full loss is scored on all rows.",
        face=COLORS["loss_bg"],
        edge=COLORS["loss_edge"],
        accent=COLORS["loss"],
        wrap=34,
    )
    largest = max(metrics, key=lambda m: m.auc_gap + m.logloss_gap)
    centered_callout(
        fig,
        (0.669, 0.038, 0.284, 0.090),
        "Current readout",
        f"{largest.label}: AUC {largest.auc_gap:+.4f}, logloss {largest.logloss_gap:+.4f}.",
        face=COLORS["yellow_bg"],
        edge=COLORS["yellow_edge"],
        accent="#8a5a00",
        wrap=34,
    )

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=SLIDE_DPI)
    plt.close(fig)
    return out_png


def write_sidecars(metrics: list[SplitMetric], payload: dict[str, Any], args: argparse.Namespace, out_png: Path) -> tuple[Path, Path]:
    metrics_rows = metric_rows(metrics)
    metrics_out = args.outdir / "the8_current_auc_logloss_diagnostic_metrics.json"
    manifest_out = args.outdir / "the8_current_auc_logloss_diagnostic_manifest.json"
    script_out = args.outdir / "the8_current_auc_logloss_diagnostic_speaker.md"

    metrics_out.write_text(
        json.dumps(
            {
                "schema": "THE8_CURRENT_AUC_LOGLOSS_DIAGNOSTIC_METRICS_V1",
                "source_metrics_json": str(args.metrics_json),
                "source_schema": payload.get("schema"),
                "rows": metrics_rows,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    manifest_out.write_text(
        json.dumps(
            {
                "schema": "THE8_CURRENT_AUC_LOGLOSS_DIAGNOSTIC_MANIFEST_V1",
                "output_png": str(out_png),
                "metrics_json": str(metrics_out),
                "speaker_md": str(script_out),
                "source_metrics_json": str(args.metrics_json),
                "source_campaign_tag": payload.get("campaign_tag"),
                "source_note": payload.get("note"),
                "target_slide_object_id": "g3ebf6faceae_1_25",
                "target_source_image": "the38_auc_logloss_average_explainer.png",
                "render_size_px": [2560, 1440],
                "audience_caveat": (
                    "Green lines are true full-validation metrics from scoring the full train+holdout row "
                    "population on SDCC, not row-weighted averages of train and holdout metrics. "
                    "Train/holdout recomputation closed against registry metrics before use."
                ),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    best = max(metrics, key=lambda m: m.auc_gap + m.logloss_gap)
    script_out.write_text(
        "\n".join(
            [
                "# Speaker Notes - Current AUC/Logloss Diagnostic",
                "",
                "This remake uses the corrected Au+Au split-lane registries rather than the old May split-study outputs.",
                "",
                "The green curve is a real full-validation metric: the full train-plus-holdout row population was scored directly on SDCC for each split model. It is not a row-weighted average of the train and holdout AUCs.",
                "",
                "The same full row scoring is used for logloss. The train and holdout metrics shown on the orange and blue curves still come from the split registries and were reproduced by the SDCC check before using the full metrics.",
                "",
                f"The practical diagnostic is the gap panel. The {best.label} split has the largest warning in the current corrected output: AUC gap {best.auc_gap:+.4f} and logloss gap {best.logloss_gap:+.4f}. The 90/10 baseline is effectively closed by comparison.",
                "",
            ]
        ),
        encoding="utf-8",
    )
    return metrics_out, manifest_out


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metrics-json", type=Path, default=DEFAULT_METRICS)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--output-name", default="the8_current_auc_logloss_diagnostic_enhanced_v10.png")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    metrics, payload = load_metrics(args.metrics_json)
    out_png = render_slide(metrics, payload, args)
    metrics_out, manifest_out = write_sidecars(metrics, payload, args, out_png)
    print("THE8_CURRENT_AUC_LOGLOSS_DIAGNOSTIC_READY")
    print(f"png={out_png}")
    print(f"metrics_json={metrics_out}")
    print(f"manifest_json={manifest_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
