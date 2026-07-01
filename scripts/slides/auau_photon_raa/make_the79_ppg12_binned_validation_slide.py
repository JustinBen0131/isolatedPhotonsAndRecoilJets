#!/usr/bin/env python3
"""Build the PPG12-binned AuAu BDT validation comparison slide."""

from __future__ import annotations

import argparse
import csv
import json
import math
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch


THIS_FILE = Path(__file__).resolve()
REPO = next((p for p in THIS_FILE.parents if (p / "AGENTS.md").exists()), THIS_FILE.parents[4])
SCRIPTS = REPO / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.append(str(SCRIPTS))

from slides.common.slide_defaults import SLIDE_DPI, SLIDE_HEIGHT_PX, SLIDE_WIDTH_PX, slide_figsize  # noqa: E402


DEFAULT_OUT_DIR = REPO / "dataOutput/auauBDTModelSelection/THE79_ppg12_binned_validation_20260627"
PPG12_EDGES = [5.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 32.0, 36.0, 40.0]
PT_BINS = list(zip(PPG12_EDGES[:-1], PPG12_EDGES[1:]))
PT_LABELS = [f"{lo:g}-{hi:g}" for lo, hi in PT_BINS]
CENT_RANGE = (0.0, 20.0)

PPG12_SCORE_KEY = "score_base14_perEtCent7"
BASELINE_SCORE_KEY = "score_centAsFeatBase3x3_pt15to35"
EXCLUDED_PT_LABELS: list[str] = []
DEFAULT_TITLE = r"PPG12-binned 5-40 routing vs default 15-35 BDT"
DEFAULT_PPG12_LABEL = "PPG12 5-40 routed"
DEFAULT_BASELINE_LABEL = "15-35 default"
PLOT_LABELS = {"ppg12": DEFAULT_PPG12_LABEL, "baseline": DEFAULT_BASELINE_LABEL}

TIMES_DIR = Path("/System/Library/Fonts/Supplemental")
TIMES_FONTS = [
    TIMES_DIR / "Times New Roman.ttf",
    TIMES_DIR / "Times New Roman Bold.ttf",
    TIMES_DIR / "Times New Roman Italic.ttf",
    TIMES_DIR / "Times New Roman Bold Italic.ttf",
]

INK = "#0F172A"
MUTED = "#475569"
GRID = "#D9E2EF"
EDGE = "#C9D4E5"
PANEL = "#FBFDFF"
SOFT_BLUE = "#EEF5FC"
GREEN = "#15803D"
ORANGE = "#E66100"
RED = "#B91C1C"


@dataclass(frozen=True)
class ModelSpec:
    key: str
    label: str
    report_dir: Path
    score_key: str
    color: str
    product: str


def setup_style() -> None:
    for font in TIMES_FONTS:
        if font.exists():
            font_manager.fontManager.addfont(str(font))
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "custom",
            "mathtext.rm": "Times New Roman",
            "mathtext.it": "Times New Roman:italic",
            "mathtext.bf": "Times New Roman:bold",
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "axes.facecolor": PANEL,
            "axes.edgecolor": EDGE,
            "axes.linewidth": 1.0,
            "xtick.color": INK,
            "ytick.color": INK,
            "axes.unicode_minus": False,
        }
    )


def font_px(points: float) -> float:
    return points * SLIDE_DPI / 72.0


def text_bbox_px(fig: plt.Figure, artist) -> list[float]:
    fig.canvas.draw()
    bbox = artist.get_window_extent(renderer=fig.canvas.get_renderer())
    scale_x = fig.bbox.width / SLIDE_WIDTH_PX
    scale_y = fig.bbox.height / SLIDE_HEIGHT_PX
    height = fig.bbox.height
    return [
        float(bbox.x0 / scale_x),
        float((height - bbox.y1) / scale_y),
        float(bbox.x1 / scale_x),
        float((height - bbox.y0) / scale_y),
    ]


def add_text(
    fig: plt.Figure,
    nodes: list[dict],
    name: str,
    x: float,
    y: float,
    text: str,
    *,
    size: float,
    role: str = "audience",
    weight: str = "normal",
    color: str = INK,
    ha: str = "left",
    va: str = "center",
    linespacing: float = 1.0,
    **flags,
) -> None:
    artist = fig.text(
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
    nodes.append(
        {
            "name": name,
            "kind": "text",
            "role": role,
            "text": text,
            "font_px": font_px(size),
            "bbox": text_bbox_px(fig, artist),
            "title_anchor": role == "title",
            **flags,
        }
    )


def rounded_box(
    fig: plt.Figure,
    nodes: list[dict],
    name: str,
    x: float,
    y: float,
    w: float,
    h: float,
    *,
    face: str,
    edge: str,
    lw: float = 1.3,
    radius: float = 0.006,
    role: str = "card",
) -> None:
    fig.patches.append(
        FancyBboxPatch(
            (x, y),
            w,
            h,
            transform=fig.transFigure,
            boxstyle=f"round,pad=0.006,rounding_size={radius}",
            facecolor=face,
            edgecolor=edge,
            linewidth=lw,
            zorder=3,
        )
    )
    nodes.append(
        {
            "name": name,
            "kind": "rect",
            "role": role,
            "bbox": [
                x * SLIDE_WIDTH_PX,
                (1.0 - y - h) * SLIDE_HEIGHT_PX,
                (x + w) * SLIDE_WIDTH_PX,
                (1.0 - y) * SLIDE_HEIGHT_PX,
            ],
        }
    )


def auc_rank(y_true: np.ndarray, score: np.ndarray) -> float:
    y = np.asarray(y_true, dtype=bool)
    s = np.asarray(score, dtype="float64")
    finite = np.isfinite(s)
    y = y[finite]
    s = s[finite]
    n_sig = int(y.sum())
    n_bkg = int((~y).sum())
    if n_sig == 0 or n_bkg == 0:
        return math.nan
    order = np.argsort(s, kind="mergesort")
    sorted_s = s[order]
    ranks = np.empty(len(s), dtype="float64")
    i = 0
    while i < len(s):
        j = i + 1
        while j < len(s) and sorted_s[j] == sorted_s[i]:
            j += 1
        ranks[order[i:j]] = (i + 1 + j) / 2.0
        i = j
    return float((ranks[y].sum() - n_sig * (n_sig + 1) / 2.0) / (n_sig * n_bkg))


def binary_logloss(y_true: np.ndarray, score: np.ndarray, *, balanced: bool = True) -> float:
    y = np.asarray(y_true, dtype="float64")
    s = np.asarray(score, dtype="float64")
    finite = np.isfinite(y) & np.isfinite(s)
    y = y[finite]
    s = s[finite]
    if len(y) == 0:
        return math.nan
    eps = 1.0e-6
    p = np.clip(s, eps, 1.0 - eps)
    loss = -(y * np.log(p) + (1.0 - y) * np.log(1.0 - p))
    if not balanced:
        return float(np.mean(loss))
    n_sig = float(np.sum(y == 1.0))
    n_bkg = float(np.sum(y == 0.0))
    if n_sig == 0.0 or n_bkg == 0.0:
        return math.nan
    weights = np.where(y == 1.0, 0.5 / n_sig, 0.5 / n_bkg)
    return float(np.sum(weights * loss))


def wp80_background_fake_rate(y_true: np.ndarray, score: np.ndarray) -> tuple[float, float]:
    y = np.asarray(y_true, dtype=bool)
    s = np.asarray(score, dtype="float64")
    finite = np.isfinite(s)
    y = y[finite]
    s = s[finite]
    sig = s[y]
    bkg = s[~y]
    if len(sig) == 0 or len(bkg) == 0:
        return math.nan, math.nan
    threshold = float(np.quantile(sig, 0.20, method="lower"))
    return threshold, float(np.mean(bkg >= threshold))


def score_cache_files(report_dir: Path) -> list[Path]:
    files = sorted((report_dir / "score_caches").glob("score_cache_*.npz"))
    if not files:
        raise FileNotFoundError(f"No score caches found under {report_dir / 'score_caches'}")
    return files


def require_score_key(cache, key: str, path: Path) -> str:
    if key in cache.files:
        return key
    matches = [item for item in cache.files if item == key or item.endswith(key.replace("score_", ""))]
    if len(matches) == 1:
        return matches[0]
    available = ", ".join(sorted(cache.files)[:30])
    raise KeyError(f"Missing score key {key} in {path}; first keys: {available}")


def load_binned_arrays(spec: ModelSpec) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    y_parts: dict[str, list[np.ndarray]] = {label: [] for label in PT_LABELS}
    s_parts: dict[str, list[np.ndarray]] = {label: [] for label in PT_LABELS}
    for path in score_cache_files(spec.report_dir):
        with np.load(path, allow_pickle=False) as cache:
            score_key = require_score_key(cache, spec.score_key, path)
            y = np.asarray(cache["is_signal"], dtype=bool)
            et = np.asarray(cache["cluster_Et"], dtype="float64")
            cent = np.asarray(cache["centrality"], dtype="float64")
            score = np.asarray(cache[score_key], dtype="float64")
        base = (cent >= CENT_RANGE[0]) & (cent < CENT_RANGE[1]) & np.isfinite(et) & np.isfinite(score)
        for (lo, hi), label in zip(PT_BINS, PT_LABELS):
            mask = base & (et >= lo) & (et < hi)
            if mask.any():
                y_parts[label].append(y[mask])
                s_parts[label].append(score[mask])

    out: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for label in PT_LABELS:
        if y_parts[label]:
            out[label] = (np.concatenate(y_parts[label]), np.concatenate(s_parts[label]))
        else:
            out[label] = (np.array([], dtype=bool), np.array([], dtype="float64"))
    return out


def compute_metrics(models: list[ModelSpec]) -> list[dict]:
    rows: list[dict] = []
    for spec in models:
        binned = load_binned_arrays(spec)
        for lo_hi, label in zip(PT_BINS, PT_LABELS):
            y, score = binned[label]
            threshold, fake = wp80_background_fake_rate(y, score)
            rows.append(
                {
                    "model_key": spec.key,
                    "model_label": spec.label,
                    "product": spec.product,
                    "pt_bin": label,
                    "pt_low_GeV": lo_hi[0],
                    "pt_high_GeV": lo_hi[1],
                    "entries": int(len(y)),
                    "signal_entries": int(y.sum()) if len(y) else 0,
                    "background_entries": int((~y).sum()) if len(y) else 0,
                    "auc": auc_rank(y, score),
                    "balanced_score_logloss": binary_logloss(y, score, balanced=True),
                    "wp80_threshold": threshold,
                    "wp80_background_pass_rate": fake,
                    "signal_score_mean": float(np.nanmean(score[y])) if len(y) and y.any() else math.nan,
                    "background_score_mean": float(np.nanmean(score[~y])) if len(y) and (~y).any() else math.nan,
                    "finite_score_fraction": float(np.isfinite(score).mean()) if len(score) else math.nan,
                }
            )
    lookup = {(row["model_key"], row["pt_bin"]): row for row in rows}
    for label in PT_LABELS:
        ppg = lookup[("ppg12", label)]
        base = lookup[("baseline", label)]
        for field in ("auc", "balanced_score_logloss", "wp80_background_pass_rate"):
            delta = math.nan
            if math.isfinite(float(ppg[field])) and math.isfinite(float(base[field])):
                delta = float(ppg[field]) - float(base[field])
            ppg[f"delta_{field}_vs_baseline"] = delta
            base[f"delta_{field}_vs_baseline"] = 0.0 if math.isfinite(float(base[field])) else math.nan
    return rows


def write_csv(path: Path, rows: list[dict]) -> None:
    fields = sorted({key for row in rows for key in row})
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def read_csv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def filter_rows_for_current_bins(rows: list[dict]) -> list[dict]:
    allowed = set(PT_LABELS)
    return [row for row in rows if row.get("pt_bin") in allowed]


def update_plot_labels(rows: list[dict], *, ppg12_label: str | None = None, baseline_label: str | None = None) -> None:
    PLOT_LABELS.update(
        {
            "ppg12": ppg12_label or DEFAULT_PPG12_LABEL,
            "baseline": baseline_label or DEFAULT_BASELINE_LABEL,
        }
    )
    for row in rows:
        key = row.get("model_key")
        if key in PLOT_LABELS:
            row["model_label"] = PLOT_LABELS[key]


def vector(rows: list[dict], model_key: str, field: str) -> np.ndarray:
    lookup = {(row["model_key"], row["pt_bin"]): row for row in rows}
    return np.asarray([float(lookup[(model_key, label)].get(field, math.nan)) for label in PT_LABELS], dtype="float64")


def decorate_axis(ax: plt.Axes, *, title: str, y_label: str | None = None, x_label: str | None = None) -> None:
    ax.set_title(title, loc="left", fontsize=17.5, fontweight="bold", color=INK, pad=8)
    if y_label:
        ax.set_ylabel(y_label, fontsize=14.0, color=INK, labelpad=8)
    if x_label:
        ax.set_xlabel(x_label, fontsize=16.0, fontweight="bold", color=INK, labelpad=8)
    ax.grid(axis="y", color=GRID, linewidth=1.1, zorder=0)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color(EDGE)
    ax.tick_params(labelsize=12.5, colors=INK)


def add_missing_baseline_band(ax: plt.Axes, baseline_auc: np.ndarray) -> None:
    missing = ~np.isfinite(baseline_auc)
    if not missing.any():
        return
    x = np.arange(len(PT_LABELS))
    for idx, is_missing in enumerate(missing):
        if is_missing:
            ax.axvspan(idx - 0.45, idx + 0.45, color="#F8FAFC", alpha=0.85, zorder=0)


def plot_panels(fig: plt.Figure, rows: list[dict]) -> None:
    x = np.arange(len(PT_LABELS))
    auc_ppg = vector(rows, "ppg12", "auc")
    auc_base = vector(rows, "baseline", "auc")
    loss_ppg = vector(rows, "ppg12", "balanced_score_logloss")
    loss_base = vector(rows, "baseline", "balanced_score_logloss")

    ax = fig.add_axes([0.070, 0.435, 0.855, 0.250])
    for idx in range(len(PT_LABELS)):
        if idx % 2 == 0:
            ax.axvspan(idx - 0.5, idx + 0.5, color=SOFT_BLUE, alpha=0.52, zorder=0)
    add_missing_baseline_band(ax, auc_base)
    ax.plot(x, auc_ppg, color=GREEN, marker="o", markersize=6.8, linewidth=2.8, zorder=4)
    finite = np.isfinite(auc_base)
    ax.plot(x[finite], auc_base[finite], color=ORANGE, marker="s", markersize=6.5, linewidth=2.8, zorder=5)
    ax.axhline(0.5, color="#94A3B8", linewidth=1.0, linestyle=(0, (4, 4)), zorder=1)
    finite_auc = np.concatenate([auc_ppg[np.isfinite(auc_ppg)], auc_base[np.isfinite(auc_base)]])
    y_min = max(0.45, float(np.nanmin(finite_auc)) - 0.085)
    y_max = min(0.98, float(np.nanmax(finite_auc)) + 0.055)
    ax.set_ylim(y_min, y_max)
    ax.set_xticks(x, [])
    decorate_axis(ax, title="ROC AUC on unseen true-5-40 validation", y_label="AUC")
    ax.legend(
        handles=[
            Line2D([0], [0], color=GREEN, marker="o", markersize=7, linewidth=2.8, label=PLOT_LABELS["ppg12"]),
            Line2D([0], [0], color=ORANGE, marker="s", markersize=7, linewidth=2.8, label=PLOT_LABELS["baseline"]),
        ],
        frameon=True,
        facecolor="white",
        edgecolor=EDGE,
        framealpha=0.92,
        loc="lower left",
        fontsize=13.2,
        ncol=2,
        handlelength=1.1,
        columnspacing=1.0,
    )

    ax2 = fig.add_axes([0.070, 0.155, 0.855, 0.205])
    for idx in range(len(PT_LABELS)):
        if idx % 2 == 0:
            ax2.axvspan(idx - 0.5, idx + 0.5, color=SOFT_BLUE, alpha=0.52, zorder=0)
    add_missing_baseline_band(ax2, loss_base)
    ax2.plot(x, loss_ppg, color=GREEN, marker="o", markersize=6.8, linewidth=2.8, zorder=4)
    finite = np.isfinite(loss_base)
    ax2.plot(x[finite], loss_base[finite], color=ORANGE, marker="s", markersize=6.5, linewidth=2.8, zorder=5)
    finite_loss = np.concatenate([loss_ppg[np.isfinite(loss_ppg)], loss_base[np.isfinite(loss_base)]])
    y_min = max(0.0, float(np.nanmin(finite_loss)) - 0.10)
    y_max = float(np.nanmax(finite_loss)) + 0.14
    ax2.set_ylim(y_min, y_max)
    ax2.set_xticks(x, PT_LABELS, rotation=35, ha="right", fontsize=11.4, fontweight="bold")
    decorate_axis(ax2, title="Class-balanced score logloss; lower is better", y_label="logloss", x_label="cluster $E_T$ bin (GeV)")


def choose_readout(rows: list[dict]) -> tuple[str, dict, dict]:
    lookup = {(row["model_key"], row["pt_bin"]): row for row in rows}
    for label in reversed(PT_LABELS):
        ppg = lookup[("ppg12", label)]
        base = lookup[("baseline", label)]
        if math.isfinite(float(ppg["auc"])) and math.isfinite(float(base["auc"])):
            return label, ppg, base
    label = PT_LABELS[-1]
    return label, lookup[("ppg12", label)], lookup[("baseline", label)]


def add_header(fig: plt.Figure, nodes: list[dict], rows: list[dict], args: argparse.Namespace) -> None:
    omitted = ""
    if EXCLUDED_PT_LABELS:
        omitted = f"; omit {', '.join(EXCLUDED_PT_LABELS)} GeV"
    sample_note = omitted.lstrip("; ") or "all PPG12 pT bins shown"
    add_text(
        fig,
        nodes,
        "title",
        0.050,
        0.947,
        args.title,
        size=30.0,
        role="title",
        weight="bold",
        va="top",
    )
    rounded_box(fig, nodes, "validation_callout_box", 0.058, 0.755, 0.884, 0.115, face="white", edge="#EF4444", lw=1.6)
    add_text(fig, nodes, "validation_callout_lead", 0.076, 0.838, "UNSEEN TRUE-5-40 VALIDATION", size=15.2, weight="bold", color=RED)
    add_text(
        fig,
        nodes,
        "validation_callout_scope",
        0.076,
        0.807,
        "0-20% centrality; same 14-feature base AuAu BDT input",
        size=14.5,
        weight="bold",
        color=INK,
    )
    add_text(
        fig,
        nodes,
        "validation_callout_body",
        0.076,
        0.778,
        r"Binned model = 98 BDTs (14 $E_T$ bins $\times$ 7 centrality bins); global model = 1 BDT.",
        size=14.5,
        weight="bold",
        color=RED,
    )


def write_manifest(out_dir: Path, rows: list[dict], args: argparse.Namespace) -> None:
    manifest = {
        "description": "PPG12-extended 5-40 AuAu BDT unseen-validation slide.",
        "centrality_range_percent": list(CENT_RANGE),
        "pt_edges_GeV": PPG12_EDGES,
        "target_product": "base14_perEtCent7",
        "target_model_count": 98,
        "target_binning": "14 pT bins x 7 centrality bins",
        "comparison_product": args.baseline_product,
        "comparison_label": PLOT_LABELS["baseline"],
        "comparison_note": args.comparison_note,
        "feature_statement": "Same 14 features: baseV3E + centrality + weta33/wphi33.",
        "training_weight_statement": "PPG12-extended campaign used ppg12-exact ET/eta training weights and no event/cross-section training weights.",
        "sample_statement": "Embedded photon signal samples Photon12/Photon20 and embedded jet background samples Jet12/Jet20/Jet30/Jet40.",
        "plotted_pt_bins": PT_LABELS,
        "excluded_pt_bins": EXCLUDED_PT_LABELS,
        "ppg12_report": str(args.ppg12_report) if args.ppg12_report else None,
        "baseline_report": str(args.baseline_report) if args.baseline_report else None,
        "metrics_csv_input": str(args.metrics_csv) if args.metrics_csv else None,
        "metrics": {
            "auc": "Rank-based ROC AUC from score caches within each pT bin and 0-20% centrality.",
            "balanced_score_logloss": "Binary logloss after clipping scores to [1e-6, 1 - 1e-6], with signal and background each contributing half the bin weight.",
            "wp80_background_pass_rate": "Background fraction passing the threshold chosen for 80% signal efficiency in the same bin.",
        },
    }
    (out_dir / "the79_ppg12_binned_validation_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    lookup = {(row["model_key"], row["pt_bin"]): row for row in rows}
    common = [label for label in PT_LABELS if math.isfinite(float(lookup[("baseline", label)]["auc"]))]
    lines = [
        "# Speaker Notes",
        "",
        "This is an apples-to-apples validation-source comparison, not a data-running result.",
        "Both lines use the same true-5-40 embedded validation source and the same 14 input features.",
        "The green curve is the PPG12-extended routed product base14_perEtCent7: 98 BDTs from 14 pT bins times 7 centrality bins.",
        f"The orange curve is {PLOT_LABELS['baseline']}; {args.comparison_note or 'see slide note for its training binning'}.",
    ]
    if common:
        last = common[-1]
        ppg = lookup[("ppg12", last)]
        base = lookup[("baseline", last)]
        lines.append(
            f"In the highest common finite bin, {last} GeV, AUC changes from {float(base['auc']):.3f} to {float(ppg['auc']):.3f} and balanced logloss changes from {float(base['balanced_score_logloss']):.2f} to {float(ppg['balanced_score_logloss']):.2f}."
        )
    (out_dir / "the79_ppg12_binned_validation_speaker_script.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_slide(args: argparse.Namespace) -> dict[str, Path]:
    global PT_BINS, PT_LABELS, EXCLUDED_PT_LABELS
    if args.exclude_pt_bin:
        excluded = set(args.exclude_pt_bin)
        kept = [(pt_bin, label) for pt_bin, label in zip(PT_BINS, PT_LABELS) if label not in excluded]
        if not kept:
            raise ValueError(f"All pT bins were excluded: {sorted(excluded)}")
        PT_BINS = [item[0] for item in kept]
        PT_LABELS = [item[1] for item in kept]
        EXCLUDED_PT_LABELS = [label for label in args.exclude_pt_bin if label in excluded]
    setup_style()
    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    if args.metrics_csv:
        rows = read_csv(args.metrics_csv)
    else:
        if not args.ppg12_report or not args.baseline_report:
            raise ValueError("--ppg12-report and --baseline-report are required unless --metrics-csv is supplied")
        models = [
            ModelSpec("ppg12", args.ppg12_label, args.ppg12_report, args.ppg12_score_key, GREEN, "base14_perEtCent7"),
            ModelSpec("baseline", args.baseline_label, args.baseline_report, args.baseline_score_key, ORANGE, args.baseline_product),
        ]
        rows = compute_metrics(models)
    rows = filter_rows_for_current_bins(rows)
    update_plot_labels(rows, ppg12_label=args.ppg12_label, baseline_label=args.baseline_label)
    csv_path = out_dir / "the79_ppg12_binned_validation_metrics.csv"
    write_csv(csv_path, rows)

    fig = plt.figure(figsize=slide_figsize(), dpi=SLIDE_DPI)
    nodes: list[dict] = []
    add_header(fig, nodes, rows, args)
    plot_panels(fig, rows)
    png_path = out_dir / "the79_ppg12_binned_validation_slide.png"
    layout_path = out_dir / "the79_ppg12_binned_validation_layout_nodes.json"
    audit_path = out_dir / "the79_ppg12_binned_validation_slide_audit.json"
    fig.savefig(png_path, dpi=SLIDE_DPI)
    plt.close(fig)
    layout_path.write_text(json.dumps({"nodes": nodes}, indent=2) + "\n", encoding="utf-8")
    write_manifest(out_dir, rows, args)

    subprocess.run(
        [
            sys.executable,
            str(SCRIPTS / "slides/common/post_render_slide_audit.py"),
            "--png",
            str(png_path),
            "--layout-nodes",
            str(layout_path),
            "--output",
            str(audit_path),
            "--require-pass",
        ],
        check=True,
    )
    return {
        "png": png_path,
        "csv": csv_path,
        "layout": layout_path,
        "audit": audit_path,
        "manifest": out_dir / "the79_ppg12_binned_validation_manifest.json",
        "speaker_script": out_dir / "the79_ppg12_binned_validation_speaker_script.md",
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ppg12-report", type=Path, help="Validation report directory for PPG12-extended binned14 model.")
    parser.add_argument("--baseline-report", type=Path, help="Validation report directory for default 15-35 baseline.")
    parser.add_argument("--metrics-csv", type=Path, help="Precomputed per-bin metrics CSV produced from the same score-cache schema.")
    parser.add_argument("--exclude-pt-bin", action="append", default=[], help="Drop a pT-bin label from the plotted view, for example 5-8.")
    parser.add_argument("--title", default=DEFAULT_TITLE)
    parser.add_argument("--ppg12-label", default=DEFAULT_PPG12_LABEL)
    parser.add_argument("--baseline-label", default=DEFAULT_BASELINE_LABEL)
    parser.add_argument("--ppg12-score-key", default=PPG12_SCORE_KEY)
    parser.add_argument("--baseline-score-key", default=BASELINE_SCORE_KEY)
    parser.add_argument("--baseline-product", default="centAsFeatBase3x3_pt15to35")
    parser.add_argument("--comparison-note", default="")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--json", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    paths = build_slide(args)
    if args.json:
        print(json.dumps({key: str(value) for key, value in paths.items()}, sort_keys=True))
    else:
        print(f"slide_png={paths['png']}")
        print(f"metrics_csv={paths['csv']}")
        print(f"audit_json={paths['audit']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
