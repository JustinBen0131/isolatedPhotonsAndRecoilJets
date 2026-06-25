#!/usr/bin/env python3
"""Render the THE-38 tree-depth capacity-control full-slide PNG."""

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
import textwrap
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch

from slides.common.slide_defaults import SLIDE_DPI, SLIDE_HEIGHT_PX, SLIDE_WIDTH_PX, slide_figsize


DEFAULT_OUTDIR = Path("dataOutput/auauMLDiagnosticRuns/THE38_tree_depth_capacity/slideReady")
DEFAULT_TITLE = "Tree depth separates underfit from overfit"
DEFAULT_SUBTITLE = (
    "Hold the 90/10 split fixed and change only BDT tree depth; the overfit signal is the train-holdout gap."
)
MODEL_ID = "centInput_pt1535"

COLORS = {
    "ink": "#171717",
    "muted": "#5b6472",
    "line": "#d7dde7",
    "grid": "#e5e9f0",
    "panel": "#f7f9fc",
    "panel_edge": "#cbd6e4",
    "train": "#d97706",
    "holdout": "#2563ad",
    "auc_gap": "#3274b9",
    "logloss_gap": "#c94f43",
    "baseline": "#111827",
    "low": "#edf6f0",
    "mid": "#edf2fb",
    "high": "#fff2d6",
}


@dataclass(frozen=True)
class DepthMetric:
    depth: int
    registry: Path
    status: str
    train_auc: float
    holdout_auc: float
    train_logloss: float
    holdout_logloss: float
    auc_gap: float
    logloss_gap: float
    train_rows: int
    holdout_rows: int
    history_csv: str
    xgboost_max_depth: int
    observed_samples: tuple[str, ...]
    sample_validation_source: str
    sample_check: str


def parse_registry_entry(text: str) -> tuple[int, Path]:
    if "=" in text:
        depth_text, path_text = text.split("=", 1)
    elif ":" in text:
        depth_text, path_text = text.split(":", 1)
    else:
        raise argparse.ArgumentTypeError(f"Registry entry must be DEPTH=PATH, got: {text}")
    try:
        depth = int(depth_text.strip().lstrip("dD"))
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"Invalid depth in registry entry: {text}") from exc
    return depth, Path(path_text).expanduser()


def finite_float(value: Any, label: str, registry: Path) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError) as exc:
        raise SystemExit(f"{registry} has non-numeric {label}: {value!r}") from exc
    if not math.isfinite(out):
        raise SystemExit(f"{registry} has non-finite {label}: {value!r}")
    return out


def observed_samples_from_report(report: dict[str, Any]) -> tuple[tuple[str, ...], str]:
    locations = [
        ("weighting.sample_validation", (report.get("weighting") or {}).get("sample_validation") or {}),
        (
            "ppg12_exact_closure.sample_validation",
            (report.get("ppg12_exact_closure") or {}).get("sample_validation") or {},
        ),
        ("sample_validation", report.get("sample_validation") or {}),
    ]
    for source, sample_validation in locations:
        observed = sample_validation.get("observed_samples") if isinstance(sample_validation, dict) else None
        if observed:
            return tuple(str(item) for item in observed), source
    return (), "not_reported"


def load_registry_metric(
    depth: int,
    registry: Path,
    *,
    model_id: str,
    required_samples: set[str],
    require_ready: bool,
) -> DepthMetric:
    if not registry.is_file():
        raise SystemExit(f"Missing registry for depth {depth}: {registry}")
    payload = json.loads(registry.read_text())
    status = str(payload.get("status", ""))
    if require_ready and status not in {"READY", "READY_WITH_SKIPS"}:
        raise SystemExit(f"{registry} status is {status!r}, expected READY")
    models = payload.get("models") or []
    if not isinstance(models, list) or not models:
        raise SystemExit(f"{registry} has no models")
    model = next((item for item in models if item.get("model_id") == model_id), None)
    if model is None:
        raise SystemExit(f"{registry} does not contain model_id={model_id}")
    report = model.get("report") or {}
    if not isinstance(report, dict):
        raise SystemExit(f"{registry} model report is missing or invalid")
    overfit = report.get("overfit_diagnostics") or {}
    if not isinstance(overfit, dict) or not overfit:
        raise SystemExit(f"{registry} is missing report.overfit_diagnostics")
    xgboost = report.get("xgboost") or {}
    max_depth = int(xgboost.get("max_depth", -1))
    if max_depth != depth:
        raise SystemExit(f"{registry} reports xgboost.max_depth={max_depth}, expected {depth}")

    observed, sample_validation_source = observed_samples_from_report(report)
    sample_check = "not_reported"
    if observed:
        missing_samples = sorted(required_samples - set(observed))
        if missing_samples:
            raise SystemExit(f"{registry} is missing required observed sample(s): {missing_samples}")
        sample_check = "ok"

    train_auc = finite_float(overfit.get("train_auc"), "train_auc", registry)
    holdout_auc = finite_float(overfit.get("holdout_auc"), "holdout_auc", registry)
    train_logloss = finite_float(overfit.get("train_logloss"), "train_logloss", registry)
    holdout_logloss = finite_float(overfit.get("holdout_logloss"), "holdout_logloss", registry)
    auc_gap = finite_float(overfit.get("auc_gap_train_minus_holdout"), "auc_gap_train_minus_holdout", registry)
    logloss_gap = finite_float(overfit.get("logloss_gap_holdout_minus_train"), "logloss_gap_holdout_minus_train", registry)
    return DepthMetric(
        depth=depth,
        registry=registry,
        status=status,
        train_auc=train_auc,
        holdout_auc=holdout_auc,
        train_logloss=train_logloss,
        holdout_logloss=holdout_logloss,
        auc_gap=auc_gap,
        logloss_gap=logloss_gap,
        train_rows=int(overfit.get("train_rows", 0) or 0),
        holdout_rows=int(overfit.get("holdout_rows", 0) or 0),
        history_csv=str(overfit.get("history_csv") or report.get("training_history_csv") or ""),
        xgboost_max_depth=max_depth,
        observed_samples=observed,
        sample_validation_source=sample_validation_source,
        sample_check=sample_check,
    )


def padded_limits(values: np.ndarray, *, frac: float = 0.18, min_pad: float = 0.001) -> tuple[float, float]:
    lo = float(np.nanmin(values))
    hi = float(np.nanmax(values))
    if math.isclose(lo, hi):
        return lo - min_pad, hi + min_pad
    pad = max(min_pad, (hi - lo) * frac)
    return lo - pad, hi + pad


def add_card(fig, xywh, face=COLORS["panel"], edge=COLORS["panel_edge"], radius=0.016):
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
            linewidth=1.2,
            zorder=-10,
        )
    )


def style_axis(ax):
    ax.set_axisbelow(True)
    ax.grid(color=COLORS["grid"], linewidth=1.0)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color("#9aa7b8")
    ax.tick_params(labelsize=14.2, colors=COLORS["muted"])
    ax.xaxis.label.set_color(COLORS["ink"])
    ax.yaxis.label.set_color(COLORS["ink"])


def add_depth_regions(ax, depths: np.ndarray, baseline_depth: int) -> None:
    min_depth = float(np.min(depths)) - 0.45
    max_depth = float(np.max(depths)) + 0.45
    ax.axvspan(min_depth, baseline_depth - 0.5, color=COLORS["low"], zorder=0)
    ax.axvspan(baseline_depth - 0.5, baseline_depth + 0.5, color=COLORS["mid"], zorder=0)
    ax.axvspan(baseline_depth + 0.5, max_depth, color=COLORS["high"], zorder=0)
    ax.axvline(baseline_depth, color=COLORS["baseline"], linewidth=1.5, linestyle=(0, (4, 3)), zorder=1)


def plot_line_panel(ax, depths, train, holdout, ylabel, title, baseline_depth, lower_is_better=False, show_xlabel=True):
    add_depth_regions(ax, depths, baseline_depth)
    ax.plot(depths, train, color=COLORS["train"], marker="o", linewidth=2.4, markersize=6.5, label="Train")
    ax.plot(depths, holdout, color=COLORS["holdout"], marker="o", linewidth=2.4, markersize=6.5, label="Holdout")
    ax.set_title(title, fontsize=17.2, fontweight="bold", pad=10, color=COLORS["ink"])
    ax.set_xlabel("XGBoost max depth" if show_xlabel else "", fontsize=14.3)
    ax.set_ylabel(ylabel, fontsize=14.3)
    ax.set_xticks(depths)
    ax.set_ylim(*padded_limits(np.concatenate([train, holdout]), min_pad=0.0007 if "AUC" in ylabel else 0.0009))
    ax.legend(frameon=False, fontsize=14.0, loc="best")
    style_axis(ax)
    direction = "lower is better" if lower_is_better else "higher is better"
    note_y = 0.075 if lower_is_better else 0.925
    note_va = "bottom" if lower_is_better else "top"
    ax.text(
        0.015,
        note_y,
        direction,
        transform=ax.transAxes,
        ha="left",
        va=note_va,
        fontsize=12.7,
        color=COLORS["muted"],
        bbox=dict(boxstyle="round,pad=0.20", facecolor="white", edgecolor="none", alpha=0.82),
    )


def plot_gap_panel(ax, depths, auc_gap, logloss_gap, baseline_depth):
    add_depth_regions(ax, depths, baseline_depth)
    width = 0.34
    x = depths.astype(float)
    ax.bar(x - width / 2, auc_gap, width=width, color=COLORS["auc_gap"], label="AUC gap")
    ax.bar(x + width / 2, logloss_gap, width=width, color=COLORS["logloss_gap"], label="Logloss gap")
    ax.axhline(0.0, color="#9aa7b8", linewidth=1.0)
    ax.set_title("Train-holdout gaps show overfit risk", fontsize=17.0, fontweight="bold", pad=10, color=COLORS["ink"])
    ax.set_xlabel("", fontsize=14.3)
    ax.set_ylabel("Gap size", fontsize=14.3)
    ax.set_xticks(depths)
    ax.set_ylim(*padded_limits(np.concatenate([auc_gap, logloss_gap]), min_pad=0.0008))
    ax.legend(frameon=False, fontsize=14.0, loc="upper right", bbox_to_anchor=(0.985, 1.36), ncol=2)
    style_axis(ax)
    for xx, yy in zip(x - width / 2, auc_gap):
        offset = 0.00030 if yy >= 0 else -0.00022
        ax.text(
            xx,
            yy + offset,
            f"{yy:.4f}",
            ha="center",
            va="bottom" if yy >= 0 else "top",
            fontsize=10.8,
            color=COLORS["auc_gap"],
        )
    for xx, yy in zip(x + width / 2, logloss_gap):
        offset = 0.00030 if yy >= 0 else -0.00022
        ax.text(
            xx,
            yy + offset,
            f"{yy:.4f}",
            ha="center",
            va="bottom" if yy >= 0 else "top",
            fontsize=10.8,
            color=COLORS["logloss_gap"],
        )


def add_reading_boxes(fig) -> None:
    boxes = [
        (
            0.054,
            0.038,
            0.265,
            0.128,
            COLORS["low"],
            "#9fc8ad",
            "Underfit cue",
            "Train and holdout are both weak, and the gap can still be small.",
        ),
        (
            0.367,
            0.038,
            0.265,
            0.128,
            COLORS["mid"],
            "#9eb5d8",
            "Useful capacity cue",
            "Holdout improves or peaks while the train-holdout gap stays controlled.",
        ),
        (
            0.680,
            0.038,
            0.265,
            0.128,
            COLORS["high"],
            "#d9b45c",
            "Overfit cue",
            "Train keeps improving while holdout stalls or worsens, so the gaps open.",
        ),
    ]
    for x, y, w, h, face, edge, title, body in boxes:
        add_card(fig, (x, y, w, h), face=face, edge=edge, radius=0.010)
        fig.text(x + 0.018, y + h - 0.026, title, fontsize=17.0, fontweight="bold", color=COLORS["ink"], va="top")
        fig.text(
            x + 0.018,
            y + h - 0.068,
            "\n".join(textwrap.wrap(body, width=40)),
            fontsize=13.5,
            color=COLORS["ink"],
            va="top",
            linespacing=1.20,
        )


def best_capacity_sentence(metrics: list[DepthMetric]) -> str:
    best_holdout = max(metrics, key=lambda item: item.holdout_auc)
    smallest_auc_gap = min(metrics, key=lambda item: abs(item.auc_gap))
    return (
        f"Holdout AUC peaks at depth {best_holdout.depth}; "
        f"smallest absolute AUC gap is depth {smallest_auc_gap.depth}."
    )


def render_slide(metrics: list[DepthMetric], args: argparse.Namespace) -> Path:
    out_png = args.outdir / args.output_basename
    depths = np.asarray([m.depth for m in metrics], dtype=int)
    train_auc = np.asarray([m.train_auc for m in metrics], dtype=float)
    holdout_auc = np.asarray([m.holdout_auc for m in metrics], dtype=float)
    train_logloss = np.asarray([m.train_logloss for m in metrics], dtype=float)
    holdout_logloss = np.asarray([m.holdout_logloss for m in metrics], dtype=float)
    auc_gap = np.asarray([m.auc_gap for m in metrics], dtype=float)
    logloss_gap = np.asarray([m.logloss_gap for m in metrics], dtype=float)

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
    fig.text(0.050, 0.925, args.title, fontsize=args.title_font_size, fontweight="bold", color=COLORS["ink"], ha="left", va="center")
    if args.subtitle:
        fig.text(0.050, 0.882, args.subtitle, fontsize=15.2, color=COLORS["muted"], ha="left", va="center")
    right_note = args.right_note if args.right_note is not None else f"Fixed split: 90/10 row holdout\nCurrent baseline: max depth = {args.baseline_depth}"
    if right_note:
        fig.text(0.950, 0.932, right_note, fontsize=12.6, color=COLORS["muted"], ha="right", va="top", linespacing=1.20)

    add_card(fig, (0.045, 0.462, 0.430, 0.373))
    add_card(fig, (0.525, 0.462, 0.430, 0.373))
    add_card(fig, (0.045, 0.195, 0.910, 0.218))

    ax_auc = fig.add_axes([0.086, 0.515, 0.348, 0.265])
    ax_log = fig.add_axes([0.566, 0.515, 0.348, 0.265])
    ax_gap = fig.add_axes([0.103, 0.243, 0.814, 0.135])
    plot_line_panel(
        ax_auc,
        depths,
        train_auc,
        holdout_auc,
        "AUC",
        "AUC: train versus unseen holdout",
        args.baseline_depth,
        show_xlabel=False,
    )
    plot_line_panel(
        ax_log,
        depths,
        train_logloss,
        holdout_logloss,
        "Logloss",
        "Logloss: probability calibration check",
        args.baseline_depth,
        lower_is_better=True,
        show_xlabel=False,
    )
    plot_gap_panel(ax_gap, depths, auc_gap, logloss_gap, args.baseline_depth)

    fig.text(0.050, 0.431, best_capacity_sentence(metrics), fontsize=15.2, color=COLORS["ink"], ha="left", va="center")
    add_reading_boxes(fig)

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=SLIDE_DPI)
    plt.close(fig)
    return out_png


def metric_dict(metric: DepthMetric) -> dict[str, Any]:
    return {
        "depth": metric.depth,
        "registry": str(metric.registry),
        "status": metric.status,
        "train_auc": metric.train_auc,
        "holdout_auc": metric.holdout_auc,
        "train_logloss": metric.train_logloss,
        "holdout_logloss": metric.holdout_logloss,
        "auc_gap_train_minus_holdout": metric.auc_gap,
        "logloss_gap_holdout_minus_train": metric.logloss_gap,
        "train_rows": metric.train_rows,
        "holdout_rows": metric.holdout_rows,
        "training_history_csv": metric.history_csv,
        "xgboost_max_depth": metric.xgboost_max_depth,
        "observed_samples": list(metric.observed_samples),
        "sample_validation_source": metric.sample_validation_source,
        "sample_check": metric.sample_check,
    }


def write_sidecars(metrics: list[DepthMetric], args: argparse.Namespace, out_png: Path) -> tuple[Path, Path, Path, Path]:
    rows = [metric_dict(metric) for metric in metrics]
    metrics_json = args.outdir / "the38_tree_depth_capacity_metrics.json"
    metrics_csv = args.outdir / "the38_tree_depth_capacity_metrics.csv"
    manifest_json = args.outdir / "the38_tree_depth_capacity_manifest.json"
    script_md = args.outdir / "the38_tree_depth_capacity_script.md"

    metrics_json.write_text(json.dumps({"schema": "THE38_TREE_DEPTH_CAPACITY_METRICS_V1", "rows": rows}, indent=2, sort_keys=True) + "\n")
    with metrics_csv.open("w", encoding="utf-8", newline="") as handle:
        fieldnames = [
            "depth",
            "train_auc",
            "holdout_auc",
            "train_logloss",
            "holdout_logloss",
            "auc_gap_train_minus_holdout",
            "logloss_gap_holdout_minus_train",
            "train_rows",
            "holdout_rows",
            "xgboost_max_depth",
            "registry",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key) for key in fieldnames})

    manifest = {
        "schema": "THE38_TREE_DEPTH_CAPACITY_SLIDE_MANIFEST_V1",
        "tag": args.tag,
        "title": args.title,
        "subtitle": args.subtitle,
        "metrics_source": str(args.metrics_json) if args.metrics_json else "model_registry_json",
        "metrics_section": args.metrics_section,
        "model_id": args.model_id,
        "baseline_depth": args.baseline_depth,
        "required_samples": sorted(args.require_sample),
        "output_png": str(out_png),
        "metrics_json": str(metrics_json),
        "metrics_csv": str(metrics_csv),
        "script_md": str(script_md),
        "render_size_px": [SLIDE_WIDTH_PX, SLIDE_HEIGHT_PX],
        "interpretation_boundary": (
            "This slide defines overfitting from train versus holdout diagnostics. "
            "Full score-cache validation, if run later, is an external scored-sample sanity check."
        ),
        "registries": [str(metric.registry) for metric in metrics],
    }
    manifest_json.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")

    best = max(metrics, key=lambda item: item.holdout_auc)
    script_md.write_text(
        "\n".join(
            [
                "# THE-38 Slide Script - Tree Depth Capacity Control",
                "",
                "This slide holds the sample split fixed at ninety percent training and ten percent holdout, and changes only the maximum tree depth of the BDT.",
                "",
                "The left panel shows AUC. The orange curve is the training rows and the blue curve is the holdout rows that the model did not train on. If both curves are weak, that points to underfitting. If the training curve keeps improving but the holdout curve stops improving, that points to overfitting.",
                "",
                "The right panel shows logloss, where lower is better. This is the probability-quality check, so it helps distinguish a model that ranks events well from one that is also assigning sensible probabilities.",
                "",
                "The bottom panel is the direct overfitting diagnostic. The AUC gap is train minus holdout, and the logloss gap is holdout minus train. The current depth-four recipe is marked as the baseline reference.",
                "",
                f"In this sweep, depth {best.depth} has the highest holdout AUC. The decision is not just which training score is largest; the useful capacity point is where holdout performance is strong and the gap is still controlled.",
                "",
                "The full scored-sample validation can still be useful afterward, but it is not what defines the train-holdout overfitting gap on this slide.",
                "",
            ]
        ),
        encoding="utf-8",
    )
    return metrics_json, metrics_csv, manifest_json, script_md


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--registry", action="append", type=parse_registry_entry, default=[], help="Depth registry entry as DEPTH=/path/model_registry.json")
    parser.add_argument("--registry-root", type=Path, default=None, help="Optional root containing d<depth>/model_registry.json folders")
    parser.add_argument("--metrics-json", type=Path, default=None, help="Optional JSON snapshot containing depth rows with train/holdout metrics")
    parser.add_argument("--metrics-section", default="rows", help="Section inside --metrics-json; use 'depth' for THE37 snapshot depth rows")
    parser.add_argument("--depths", default="2,3,4,5,6")
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--output-basename", default="the38_tree_depth_capacity_control_slide.png")
    parser.add_argument("--tag", default="THE38_tree_depth_capacity")
    parser.add_argument("--title", default=DEFAULT_TITLE)
    parser.add_argument("--title-font-size", type=float, default=29.0)
    parser.add_argument("--subtitle", default=DEFAULT_SUBTITLE)
    parser.add_argument("--right-note", default=None)
    parser.add_argument("--baseline-depth", type=int, default=4)
    parser.add_argument("--model-id", default=MODEL_ID)
    parser.add_argument("--require-sample", action="append", default=["run28_embeddedJet40"])
    parser.add_argument("--allow-non-ready", action="store_true", help="Allow registries whose status is not READY; intended only for synthetic smoke tests.")
    return parser.parse_args()


def registry_entries(args: argparse.Namespace) -> list[tuple[int, Path]]:
    entries = list(args.registry)
    if entries:
        return entries
    if args.registry_root is None:
        raise SystemExit("Provide --registry DEPTH=PATH entries or --registry-root")
    out = []
    for item in args.depths.split(","):
        depth = int(item.strip())
        candidates = [
            args.registry_root / f"d{depth}" / "model_registry.json",
            args.registry_root / f"depth{depth}" / "model_registry.json",
        ]
        matches = [path for path in candidates if path.exists()]
        if not matches:
            matches = sorted(args.registry_root.glob(f"*d{depth}*/model_registry.json"))
        if len(matches) != 1:
            raise SystemExit(f"Could not uniquely locate registry for depth {depth} under {args.registry_root}: {matches}")
        out.append((depth, matches[0]))
    return out


def _row_float(row: dict[str, Any], names: tuple[str, ...], label: str, source: Path) -> float:
    for name in names:
        if name in row and row[name] is not None:
            return finite_float(row[name], label, source)
    raise SystemExit(f"{source} row for depth {row.get('depth')} is missing {label}; tried {names}")


def load_metrics_json_metrics(args: argparse.Namespace) -> list[DepthMetric]:
    if args.metrics_json is None:
        raise SystemExit("Internal error: load_metrics_json_metrics called without --metrics-json")
    source = args.metrics_json.expanduser()
    payload = json.loads(source.read_text())
    rows_obj: Any
    if args.metrics_section == ".":
        rows_obj = payload
    else:
        rows_obj = payload.get(args.metrics_section)
    if not isinstance(rows_obj, list):
        raise SystemExit(f"{source} section {args.metrics_section!r} is not a list of metric rows")
    wanted = {int(item.strip()) for item in args.depths.split(",") if item.strip()}
    metrics: list[DepthMetric] = []
    for row_obj in rows_obj:
        if not isinstance(row_obj, dict):
            raise SystemExit(f"{source} contains a non-object metric row: {row_obj!r}")
        if "depth" not in row_obj:
            continue
        depth = int(row_obj["depth"])
        if wanted and depth not in wanted:
            continue
        registry_text = str(row_obj.get("registry") or row_obj.get("source_registry") or source)
        max_depth = int(row_obj.get("xgboost_max_depth") or row_obj.get("max_depth") or depth)
        if max_depth != depth:
            raise SystemExit(f"{source} row depth {depth} reports xgboost max_depth {max_depth}")
        observed = row_obj.get("observed_samples") or []
        if not isinstance(observed, list):
            observed = []
        metrics.append(
            DepthMetric(
                depth=depth,
                registry=Path(registry_text),
                status=str(row_obj.get("status", "READY")),
                train_auc=_row_float(row_obj, ("train_auc",), "train_auc", source),
                holdout_auc=_row_float(row_obj, ("holdout_auc",), "holdout_auc", source),
                train_logloss=_row_float(row_obj, ("train_logloss",), "train_logloss", source),
                holdout_logloss=_row_float(row_obj, ("holdout_logloss",), "holdout_logloss", source),
                auc_gap=_row_float(row_obj, ("auc_gap_train_minus_holdout", "auc_gap"), "auc_gap", source),
                logloss_gap=_row_float(row_obj, ("logloss_gap_holdout_minus_train", "logloss_gap"), "logloss_gap", source),
                train_rows=int(row_obj.get("train_rows", 0) or 0),
                holdout_rows=int(row_obj.get("holdout_rows", 0) or 0),
                history_csv=str(row_obj.get("training_history_csv") or row_obj.get("history_csv") or ""),
                xgboost_max_depth=max_depth,
                observed_samples=tuple(str(item) for item in observed),
                sample_validation_source=str(row_obj.get("sample_validation_source", "snapshot")),
                sample_check=str(row_obj.get("sample_check", "snapshot")),
            )
        )
    missing = sorted(wanted - {metric.depth for metric in metrics})
    if missing:
        raise SystemExit(f"{source} section {args.metrics_section!r} is missing requested depth(s): {missing}")
    return metrics


def main() -> int:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    if args.metrics_json:
        metrics = load_metrics_json_metrics(args)
    else:
        required_samples = set(args.require_sample or [])
        metrics = [
            load_registry_metric(
                depth,
                registry,
                model_id=args.model_id,
                required_samples=required_samples,
                require_ready=not args.allow_non_ready,
            )
            for depth, registry in registry_entries(args)
        ]
    metrics = sorted(metrics, key=lambda item: item.depth)
    if len({metric.depth for metric in metrics}) != len(metrics):
        raise SystemExit("Duplicate depth entries were provided")
    if args.baseline_depth not in {metric.depth for metric in metrics}:
        raise SystemExit(f"Baseline depth {args.baseline_depth} is not present in the sweep")
    out_png = render_slide(metrics, args)
    metrics_json, metrics_csv, manifest_json, script_md = write_sidecars(metrics, args, out_png)
    print("THE38_TREE_DEPTH_CAPACITY_SLIDE_READY")
    print(f"png={out_png}")
    print(f"metrics_json={metrics_json}")
    print(f"metrics_csv={metrics_csv}")
    print(f"manifest_json={manifest_json}")
    print(f"script_md={script_md}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
