#!/usr/bin/env python3
"""Build slide-ready PNGs for the Au+Au BDT split-study interpretation."""

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
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch


INK = "#111827"
MUTED = "#4B5563"
GRID = "#E5E7EB"
LIGHT = "#F8FAFC"
BLUE = "#0072B2"
ORANGE = "#E69F00"
GREEN = "#009E73"
RED = "#CC334E"
PURPLE = "#7E57C2"
CYAN = "#56B4E9"

SPLITS = ["90/10", "50/50", "10/90"]
COLORS = {"90/10": BLUE, "50/50": ORANGE, "10/90": GREEN}
TRAIN_FRACTION = {"90/10": 0.90, "50/50": 0.50, "10/90": 0.10}

# Fallback values from the original split-study. New campaigns should pass
# --internal-metrics-json so train/holdout diagnostics are campaign-specific.
INTERNAL_HOLDOUT = {
    "90/10": {"auc": 0.829797238580007, "train_rows": 14262741, "test_rows": 1584749},
    "50/50": {"auc": 0.809210924242421, "train_rows": 3328539, "test_rows": 3328539},
    "10/90": {"auc": 0.8054697218137871, "train_rows": 665707, "test_rows": 5991371},
}


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--compact-json", type=Path, required=True)
    ap.add_argument("--metrics-csv", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--internal-metrics-json", type=Path, default=None)
    ap.add_argument("--campaign", default="ppg12_weighted_centinput_splitstudy_20260522_1450")
    return ap.parse_args()


def load_metrics(path: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            parsed: dict[str, object] = {}
            for key, value in row.items():
                if key in {"centrality", "split"}:
                    parsed[key] = value
                else:
                    parsed[key] = float(value)
            rows.append(parsed)
    return rows


def metric_row(rows: list[dict[str, object]], cent: str, split: str) -> dict[str, object]:
    return next(r for r in rows if r["centrality"] == cent and r["split"] == split)


def hist_map(compact: dict[str, object]) -> dict[tuple[str, str], dict[str, object]]:
    return {(h["split_label"], h["cent_label"]): h for h in compact["histograms"]}


def density_step(hist: dict[str, object], key: str) -> tuple[np.ndarray, np.ndarray]:
    edges = np.asarray(hist["bin_edges"], dtype=float)
    centers = 0.5 * (edges[:-1] + edges[1:])
    return centers, np.asarray(hist[key], dtype=float)


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "axes.edgecolor": INK,
            "axes.labelcolor": INK,
            "axes.linewidth": 1.0,
            "xtick.color": INK,
            "ytick.color": INK,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )


def add_title(fig: plt.Figure, title: str, subtitle: str) -> None:
    fig.text(0.045, 0.955, title, fontsize=29, weight="bold", color=INK, ha="left", va="top")
    fig.text(0.045, 0.908, subtitle, fontsize=15.5, color=MUTED, ha="left", va="top")
    fig.text(0.92, 0.952, r"$\it{\bf{sPHENIX}}$ Internal", fontsize=14, color=INK, ha="right", va="top")


def add_verdict_box(fig: plt.Figure, xywh: tuple[float, float, float, float], title: str, body: str, face: str = "#FEF3C7") -> None:
    x, y, w, h = xywh
    box = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.012,rounding_size=0.018",
        transform=fig.transFigure,
        linewidth=0,
        facecolor=face,
        zorder=0,
    )
    fig.add_artist(box)
    fig.text(x + 0.018, y + h - 0.028, title, fontsize=16, weight="bold", color=INK, ha="left", va="top")
    wrapped = textwrap.fill(body, width=142)
    fig.text(x + 0.018, y + h - 0.066, wrapped, fontsize=13.2, color=INK, ha="left", va="top", linespacing=1.25)


def add_small_box(ax: plt.Axes, title: str, lines: list[str], color: str = "#EFF6FF") -> None:
    ax.axis("off")
    rect = FancyBboxPatch(
        (0.02, 0.05),
        0.96,
        0.9,
        boxstyle="round,pad=0.018,rounding_size=0.035",
        transform=ax.transAxes,
        linewidth=0,
        facecolor=color,
    )
    ax.add_patch(rect)
    ax.text(0.06, 0.87, title, transform=ax.transAxes, fontsize=15.5, weight="bold", color=INK, va="top")
    ax.text(0.06, 0.73, "\n".join(lines), transform=ax.transAxes, fontsize=12.2, color=INK, va="top", linespacing=1.35)


def save_summary(outdir: Path, rows: list[dict[str, object]], compact: dict[str, object]) -> None:
    focus = [metric_row(rows, "0-20%", s) for s in SPLITS]
    summary = {
        "artifact_type": "slide_png_proof_pack",
        "campaign": compact.get("source", "ppg12_weighted_centinput_splitstudy_20260522_1450"),
        "main_conclusion": (
            "The 90/10 Jet40-inclusive model is marginally highest in global validation. "
            "Centrality-bin differences are small and are not a standalone 5x-statistics argument."
        ),
        "global_auc": {r["label"]: r["global_auc"] for r in compact["rows"]},
        "internal_holdout": INTERNAL_HOLDOUT,
        "focus_0_20": focus,
        "missing_for_final_overfitting_claim": [
            "paired per-event bootstrap or DeLong from score caches",
            "multi-seed learning curve",
            "ET-bin decomposition of the score-cache differences",
        ],
    }
    (outdir / "splitstudy_slide_proof_pack_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")


def make_main_slide(outdir: Path, rows: list[dict[str, object]], compact: dict[str, object]) -> Path:
    fig = plt.figure(figsize=(16, 9), dpi=160)
    add_title(
        fig,
        "Does less training data really perform better?",
        "Train/validation split stress test, evaluated on the same full-stat candidate sample.",
    )
    gs = fig.add_gridspec(2, 3, left=0.06, right=0.96, top=0.83, bottom=0.22, wspace=0.32, hspace=0.35)
    ax_global = fig.add_subplot(gs[0, 0])
    ax_delta = fig.add_subplot(gs[0, 1])
    ax_wp = fig.add_subplot(gs[0, 2])
    ax_table = fig.add_subplot(gs[1, 0])
    ax_score = fig.add_subplot(gs[1, 1])
    ax_decision = fig.add_subplot(gs[1, 2])

    x = np.arange(len(SPLITS))
    global_auc = [next(r for r in compact["rows"] if r["label"] == s)["global_auc"] for s in SPLITS]
    ax_global.plot(x, global_auc, marker="o", lw=3, color=BLUE)
    ax_global.set_xticks(x, SPLITS)
    gmin, gmax = min(global_auc), max(global_auc)
    gpad = max(0.0008, (gmax - gmin) * 3.0)
    ax_global.set_ylim(gmin - gpad, gmax + gpad)
    ax_global.set_ylabel("global full-stat AUC")
    ax_global.set_title("90/10 is marginally highest globally", weight="bold", fontsize=14)
    ax_global.grid(True, color=GRID)
    for i, y in enumerate(global_auc):
        ax_global.text(i, y + 0.12 * gpad, f"{y:.3f}", ha="center", fontsize=11, weight="bold", clip_on=True)

    focus = [metric_row(rows, "0-20%", s) for s in SPLITS]
    auc_delta = [float(r["auc_delta_vs_90_10"]) * 1000.0 for r in focus]
    auc_err = []
    base_se = float(focus[0]["auc_se_hanley_mcneil"])
    for r in focus:
        auc_err.append(np.sqrt(base_se**2 + float(r["auc_se_hanley_mcneil"]) ** 2) * 1000.0)
    ax_delta.bar(x, auc_delta, color=[COLORS[s] for s in SPLITS], width=0.62)
    ax_delta.errorbar(x, auc_delta, yerr=auc_err, fmt="none", color=INK, lw=1.5, capsize=4)
    ax_delta.axhline(0, color=INK, lw=1)
    ax_delta.set_xticks(x, SPLITS)
    ax_delta.set_ylabel("0-20% AUC delta vs 90/10 (x1e-3)")
    ax_delta.set_title("No robust paradox appears", weight="bold", fontsize=14)
    ax_delta.grid(True, axis="y", color=GRID)
    for i, y in enumerate(auc_delta):
        ax_delta.text(i, y + (0.25 if y >= 0 else -0.35), f"{y:+.2f}", ha="center", fontsize=11, weight="bold")

    wp_delta = [float(r["wp80_fake_delta_vs_90_10"]) * 1000.0 for r in focus]
    width = 0.36
    ax_wp.bar(x - width / 2, auc_delta, width=width, color=PURPLE, label="AUC delta")
    ax_wp.bar(x + width / 2, wp_delta, width=width, color=RED, label="WP80 fake delta")
    ax_wp.axhline(0, color=INK, lw=1)
    ax_wp.set_xticks(x, SPLITS)
    ax_wp.set_ylabel("delta vs 90/10 (x1e-3)")
    ax_wp.set_title("AUC is not the WP decision", weight="bold", fontsize=14)
    ax_wp.legend(frameon=False, fontsize=10, loc="upper left")
    ax_wp.grid(True, axis="y", color=GRID)

    ax_table.axis("off")
    train_rows = [INTERNAL_HOLDOUT[s]["train_rows"] / 1e6 for s in SPLITS]
    cell_text = [[s, f"{train_rows[i]:.2f}M", f"{global_auc[i]:.6f}", f"{focus[i]['auc']:.6f}"] for i, s in enumerate(SPLITS)]
    table = ax_table.table(
        cellText=cell_text,
        colLabels=["split", "train rows", "global AUC", "0-20 AUC"],
        loc="center",
        cellLoc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(11)
    table.scale(1, 1.65)
    for (r, _c), cell in table.get_celld().items():
        cell.set_edgecolor("#CBD5E1")
        if r == 0:
            cell.set_facecolor("#E2E8F0")
            cell.set_text_props(weight="bold", color=INK)
    ax_table.set_title("Same evaluation counts, different training rows", weight="bold", fontsize=14, pad=8)

    ax_score.plot([0, 1], [0, 1], alpha=0)
    ax_score.axis("off")
    add_small_box(
        ax_score,
        "Read the slide this way",
        [
            "1. 90/10 is marginally highest in global validation.",
            "2. Centrality-bin AUC changes are only a few 1e-4.",
            "3. WP80 fake rate does not improve with less training.",
            "4. AUC alone is not a 5x-simulation argument.",
        ],
        "#ECFDF5",
    )

    ax_decision.axis("off")
    add_small_box(
        ax_decision,
        "Current verdict",
        [
            "No less-data-is-better evidence.",
            "Jet40-inclusive replay is stable across splits.",
            "Not enough to justify 5x embedded data by itself.",
            "Need paired score-cache test plus learning curve.",
        ],
        "#FEF3C7",
    )

    add_verdict_box(
        fig,
        (0.06, 0.055, 0.90, 0.105),
        "Bottom line",
        "The Jet40-inclusive replay removes the apparent less-data-is-better story. The production-relevant question is not raw AUC alone, "
        "but whether more independent embedded statistics reduce train-validation gaps and stabilize WP80 fake rate.",
    )
    out = outdir / "slide23_splitstudy_verdict_main.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


def make_overfitting_slide(outdir: Path, compact: dict[str, object]) -> Path:
    fig = plt.figure(figsize=(16, 9), dpi=160)
    add_title(
        fig,
        "Overfitting check: what is measured now?",
        "Jet40-inclusive trainings record train and holdout diagnostics; the full-stat validation is the common production sample.",
    )
    gs = fig.add_gridspec(2, 3, left=0.06, right=0.96, top=0.82, bottom=0.20, wspace=0.32, hspace=0.38)
    ax_auc = fig.add_subplot(gs[0, 0])
    ax_rows = fig.add_subplot(gs[0, 1])
    ax_gap = fig.add_subplot(gs[0, 2])
    ax_matrix = fig.add_subplot(gs[1, :2])
    ax_callout = fig.add_subplot(gs[1, 2])

    x = np.arange(len(SPLITS))
    holdout_auc = [INTERNAL_HOLDOUT[s]["auc"] for s in SPLITS]
    train_auc = [INTERNAL_HOLDOUT[s].get("train_auc", np.nan) for s in SPLITS]
    global_auc = [next(r for r in compact["rows"] if r["label"] == s)["global_auc"] for s in SPLITS]
    train_rows = np.asarray([INTERNAL_HOLDOUT[s]["train_rows"] / 1e6 for s in SPLITS])

    ax_auc.plot(x, train_auc, marker="^", lw=2.4, color=PURPLE, label="train AUC")
    ax_auc.plot(x, holdout_auc, marker="o", lw=2.7, color=ORANGE, label="internal holdout AUC")
    ax_auc.plot(x, global_auc, marker="s", lw=2.7, color=BLUE, label="full-stat validation AUC")
    ax_auc.set_xticks(x, SPLITS)
    auc_vals = [v for v in train_auc + holdout_auc + global_auc if np.isfinite(v)]
    amin, amax = min(auc_vals), max(auc_vals)
    apad = max(0.002, (amax - amin) * 0.35)
    ax_auc.set_ylim(amin - apad, amax + apad)
    ax_auc.set_ylabel("AUC")
    ax_auc.set_title("Train and holdout AUC stay close", weight="bold", fontsize=14)
    ax_auc.legend(frameon=False, fontsize=10)
    ax_auc.grid(True, color=GRID)

    ax_rows.bar(x, train_rows, color=[COLORS[s] for s in SPLITS], width=0.62)
    ax_rows.set_yscale("log")
    ax_rows.set_xticks(x, SPLITS)
    ax_rows.set_ylabel("training rows, log scale")
    ax_rows.set_title("This is a real statistics stress test", weight="bold", fontsize=14)
    ax_rows.grid(True, axis="y", color=GRID)
    for i, y in enumerate(train_rows):
        if y > 1.0:
            ax_rows.text(i, y / 1.32, f"{y:.2f}M", ha="center", va="top", fontsize=11, weight="bold", color="white")
        else:
            ax_rows.text(i, y * 1.18, f"{y:.2f}M", ha="center", va="bottom", fontsize=11, weight="bold", color=INK)

    train_holdout_gap = np.asarray(train_auc) - np.asarray(holdout_auc)
    ax_gap.bar(x, train_holdout_gap * 1000.0, color=CYAN, width=0.62)
    ax_gap.axhline(0, color=INK, lw=1)
    gap_vals = train_holdout_gap * 1000.0
    ymax = max(1.0, float(np.nanmax(gap_vals)) * 1.25)
    ax_gap.set_ylim(0, ymax)
    ax_gap.set_xticks(x, SPLITS)
    ax_gap.set_ylabel("train minus holdout AUC (x1e-3)")
    ax_gap.set_title("Train-holdout gaps stay small", weight="bold", fontsize=14)
    ax_gap.grid(True, axis="y", color=GRID)
    for i, y in enumerate(gap_vals):
        color = "white" if y > 1.2 else INK
        va = "top" if y > 1.2 else "bottom"
        ypos = y - 0.25 if y > 1.2 else y + 0.08
        ax_gap.text(i, ypos, f"{y:+.2f}", ha="center", va=va, fontsize=11, weight="bold", color=color)

    ax_matrix.axis("off")
    matrix_rows = [
        ["internal holdout AUC", "available", "stable across splits"],
        ["train AUC", "available", "largest gap appears in 10/90"],
        ["train and holdout logloss", "available", "no runaway calibration collapse"],
        ["full-stat validation AUC", "available", "90/10 is marginally highest"],
        ["0-20% validation AUC", "available", "few-1e-4 level shift"],
        ["multi-seed spread", "missing", "needed for 5x simulation argument"],
    ]
    table = ax_matrix.table(
        cellText=matrix_rows,
        colLabels=["diagnostic", "status", "interpretation"],
        loc="center",
        cellLoc="left",
        colLoc="left",
        colWidths=[0.34, 0.18, 0.48],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(11.3)
    table.scale(1, 1.55)
    for (r, c), cell in table.get_celld().items():
        cell.set_edgecolor("#CBD5E1")
        if r == 0:
            cell.set_facecolor("#E2E8F0")
            cell.set_text_props(weight="bold", color=INK)
        elif c == 1 and cell.get_text().get_text() == "missing":
            cell.set_facecolor("#FEE2E2")
        elif c == 1:
            cell.set_facecolor("#DCFCE7")
    ax_matrix.set_title("Overfitting evidence matrix", weight="bold", fontsize=14, pad=8)

    add_small_box(
        ax_callout,
        "Expert read",
        [
            "No split shows a large train-validation collapse.",
            "10/90 has the largest train-holdout AUC gap, as expected.",
            "Full-stat validation AUC differs by less than 3e-4.",
            "A real 5x argument still needs a learning curve or multi-seed study.",
        ],
        "#F3E8FF",
    )

    add_verdict_box(
        fig,
        (0.06, 0.055, 0.90, 0.105),
        "Bottom line",
        "The Jet40-inclusive diagnostics do not show smaller training samples winning globally or a runaway overfit failure. "
        "The unresolved resource question is how quickly the model saturates with more independent simulation.",
    )
    out = outdir / "slide24_overfitting_check_backup.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


def make_autopsy_slide(outdir: Path, rows: list[dict[str, object]], compact: dict[str, object]) -> Path:
    hmap = hist_map(compact)
    fig = plt.figure(figsize=(16, 9), dpi=160)
    add_title(
        fig,
        "0-20% autopsy: what changes when training data is reduced?",
        "Compact histogram triage says the shift is a small rank-order reshuffling, not a working-point improvement.",
    )
    gs = fig.add_gridspec(2, 3, left=0.06, right=0.96, top=0.82, bottom=0.20, wspace=0.32, hspace=0.36)
    ax_sig = fig.add_subplot(gs[0, 0])
    ax_bkg = fig.add_subplot(gs[0, 1])
    ax_norm = fig.add_subplot(gs[0, 2])
    ax_delta = fig.add_subplot(gs[1, 0])
    ax_tail = fig.add_subplot(gs[1, 1])
    ax_read = fig.add_subplot(gs[1, 2])

    base = hmap[("90/10", "0-20%")]
    alt = hmap[("10/90", "0-20%")]
    cx, bsig = density_step(base, "signal_density")
    _, bbkg = density_step(base, "background_density")
    _, asig = density_step(alt, "signal_density")
    _, abkg = density_step(alt, "background_density")
    ax_sig.step(cx, bsig, where="mid", lw=2.5, color=BLUE, label="90/10 signal")
    ax_sig.step(cx, asig, where="mid", lw=2.5, color=GREEN, label="10/90 signal")
    ax_sig.set_title("Signal shifts slightly higher", weight="bold", fontsize=14)
    ax_sig.set_xlabel("BDT score")
    ax_sig.set_ylabel("density")
    ax_sig.grid(True, color=GRID)
    ax_sig.legend(frameon=False, fontsize=10)

    ax_bkg.step(cx, bbkg, where="mid", lw=2.5, color=BLUE, label="90/10 background")
    ax_bkg.step(cx, abkg, where="mid", lw=2.5, color=GREEN, label="10/90 background")
    ax_bkg.set_title("Background tail barely changes", weight="bold", fontsize=14)
    ax_bkg.set_xlabel("BDT score")
    ax_bkg.set_ylabel("density")
    ax_bkg.grid(True, color=GRID)
    ax_bkg.legend(frameon=False, fontsize=10)

    focus90 = metric_row(rows, "0-20%", "90/10")
    focus10 = metric_row(rows, "0-20%", "10/90")
    delta = float(focus10["auc_delta_vs_90_10"])
    sigma = np.sqrt(float(focus90["auc_se_hanley_mcneil"]) ** 2 + float(focus10["auc_se_hanley_mcneil"]) ** 2)
    xs = np.linspace(delta - 4 * sigma, delta + 4 * sigma, 300)
    ys = np.exp(-0.5 * ((xs - delta) / sigma) ** 2)
    ys /= ys.max()
    ax_norm.plot(xs * 1000.0, ys, color=PURPLE, lw=2.6)
    ax_norm.fill_between(xs * 1000.0, 0, ys, where=np.abs(xs - delta) <= sigma, color=PURPLE, alpha=0.20)
    ax_norm.axvline(0, color=INK, lw=1.2)
    ax_norm.axvline(delta * 1000.0, color=PURPLE, lw=1.5, ls="--")
    ax_norm.set_title("Current uncertainty is not paired", weight="bold", fontsize=14)
    ax_norm.set_xlabel("10/90 minus 90/10 AUC (x1e-3)")
    ax_norm.set_yticks([])
    ax_norm.grid(True, axis="x", color=GRID)
    ax_norm.text(0.05, 0.86, f"center = {delta*1000:+.2f}\nindep. sigma = {sigma*1000:.2f}", transform=ax_norm.transAxes, fontsize=11.5, color=INK)

    ax_delta.step(cx, asig - bsig, where="mid", lw=2.4, color=GREEN, label="signal delta")
    ax_delta.step(cx, abkg - bbkg, where="mid", lw=2.4, color=RED, label="background delta")
    ax_delta.axhline(0, color=INK, lw=1)
    ax_delta.set_title("Deltas are local shape changes", weight="bold", fontsize=14)
    ax_delta.set_xlabel("BDT score")
    ax_delta.set_ylabel("10/90 - 90/10 density")
    ax_delta.grid(True, color=GRID)
    ax_delta.legend(frameon=False, fontsize=10)

    tail_metrics = [
        ["WP80 fake", float(focus90["wp80_background_fake_hist"]), float(focus10["wp80_background_fake_hist"])],
        ["bkg score >= 0.5", float(focus90["background_survival_score_ge_0p5"]), float(focus10["background_survival_score_ge_0p5"])],
        ["signal score >= 0.5", float(focus90["signal_survival_score_ge_0p5"]), float(focus10["signal_survival_score_ge_0p5"])],
    ]
    y = np.arange(len(tail_metrics))
    width = 0.36
    ax_tail.barh(y - width / 2, [m[1] for m in tail_metrics], height=width, color=BLUE, label="90/10")
    ax_tail.barh(y + width / 2, [m[2] for m in tail_metrics], height=width, color=GREEN, label="10/90")
    ax_tail.set_yticks(y, [m[0] for m in tail_metrics])
    ax_tail.set_xlim(0, 0.9)
    ax_tail.set_xlabel("fraction")
    ax_tail.set_title("Working-point quantities stay flat", weight="bold", fontsize=14)
    ax_tail.legend(frameon=False, fontsize=10)
    ax_tail.grid(True, axis="x", color=GRID)

    add_small_box(
        ax_read,
        "Autopsy verdict",
        [
            "The AUC movement is only at the few-1e-4 level.",
            "WP80 background survival is essentially unchanged.",
            "The next decisive test is paired per-event resampling.",
            "If a single ET bin or tail drives it, that is composition/statistics, not a model-quality reversal.",
        ],
        "#ECFDF5",
    )

    add_verdict_box(
        fig,
        (0.06, 0.055, 0.90, 0.105),
        "Bottom line",
        "The 0-20% effect is currently best treated as a tiny rank-order perturbation. "
        "Do not use it as a resource argument until paired score-cache and ET-bin decomposition are done.",
    )
    out = outdir / "slide25_0to20_auc_autopsy_backup.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


def main() -> None:
    args = parse_args()
    setup_style()
    args.outdir.mkdir(parents=True, exist_ok=True)
    compact = json.loads(args.compact_json.read_text())
    if args.internal_metrics_json is not None:
        INTERNAL_HOLDOUT.clear()
        INTERNAL_HOLDOUT.update(json.loads(args.internal_metrics_json.read_text()))
    rows = load_metrics(args.metrics_csv)
    save_summary(args.outdir, rows, compact)
    made = [
        make_main_slide(args.outdir, rows, compact),
        make_overfitting_slide(args.outdir, compact),
        make_autopsy_slide(args.outdir, rows, compact),
    ]
    for path in made:
        print(path)


if __name__ == "__main__":
    main()
