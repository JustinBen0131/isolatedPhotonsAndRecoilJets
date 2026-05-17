#!/usr/bin/env python3
"""Build an expert-facing AuAu MLP / BDT+MLP-stack validation plot pack.

The plots are meant to answer the questions ML reviewers usually ask first:
did training converge, is the score useful in the target phase space, is the
working point stable, is calibration sane, and does the score correlate with
physics-control observables such as isolation.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import train_auau_stacked_bdt_mlp_sweep as stack  # noqa: E402
from make_auau_bdt_mlp_stack_roc_overlay import route_score  # noqa: E402


DEFAULT_OUTDIR = Path("dataOutput/auauMLPDiagnosticPlots/expert_validation_pack_20260515")
DEFAULT_MLP_CACHE = Path(
    "dataOutput/auauMLPDiagnosticPlots/stack_score_separation_fineEt_cent7_20260514/local_aligned_mlp_score_caches.list"
)
DEFAULT_BDT_CACHE = Path(
    "dataOutput/auauMLPDiagnosticPlots/stack_score_separation_fineEt_cent7_20260514/local_aligned_bdt_score_caches.list"
)
DEFAULT_STACK_ARTIFACT = Path(
    "dataOutput/auauBDTMLPStackPromotion/"
    "bdt_mlp_stack_nn_wp80_uncapped_diagnostic_20260514_1915_nnstack_wp80_uncapped/"
    "artifacts/ptFine15to35_cent7_full_nn.json"
)
DEFAULT_MLP_HISTORY = Path(
    "dataOutput/auauMLPDiagnosticPlots/aligned_stack_score_and_training_20260513/"
    "training_curves/mlp_training_history_combined.csv"
)
DEFAULT_MLP_METRICS = Path(
    "dataOutput/auauTightMLPValidation/"
    "mlp_model_validation_condor_deep_primary_ratios_nostat_fullval_20260512_145449/"
    "validation_metrics.json"
)
DEFAULT_STACK_RANK = Path(
    "dataOutput/auauBDTMLPStackFullStat/"
    "stacked_bdt_mlp_targeted_fullstat_ptFine15to35_cent7_full_20260513_123934/"
    "stacked_sweep_rank_table.csv"
)
DEFAULT_WP_CELLS = Path(
    "dataOutput/auauBDTMLPStackPromotion/"
    "bdt_mlp_stack_nn_wp80_uncapped_diagnostic_20260514_1915_nnstack_wp80_uncapped/"
    "stack_working_points_target80_cells.csv"
)
DEFAULT_WP_DIAG_DIR = Path(
    "dataOutput/auauBDTMLPStackPromotion/"
    "bdt_mlp_stack_nn_wp80_uncapped_diagnostic_20260514_1915_nnstack_wp80_uncapped/diagnostics"
)
DERIVED_STACK_FEATURES = {
    "cluster_weta_over_wphi",
    "cluster_weta33_over_wphi33",
    "mlp_logit",
    "bdt_is_finite",
    "mlp_is_finite",
    "log_cluster_Et",
    "centrality_scaled",
    "bdt_x_mlp_logit",
    "bdt_x_logEt",
    "mlp_logit_x_logEt",
    "bdt_x_cent",
    "mlp_logit_x_cent",
    "logEt_x_cent",
}


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def fval(row: dict[str, str], key: str, default: float = math.nan) -> float:
    try:
        value = float(row.get(key, ""))
    except Exception:
        return default
    return value if math.isfinite(value) else default


def finite_corr(a: np.ndarray, b: np.ndarray) -> float:
    mask = np.isfinite(a) & np.isfinite(b)
    if int(mask.sum()) < 3:
        return math.nan
    return float(np.corrcoef(a[mask], b[mask])[0, 1])


def auc_score(y_true: np.ndarray, score: np.ndarray) -> float:
    y = np.asarray(y_true, dtype="int32")
    s = np.asarray(score, dtype="float64")
    mask = np.isfinite(s) & np.isin(y, [0, 1])
    y = y[mask]
    s = s[mask]
    n_pos = int(np.sum(y == 1))
    n_neg = int(np.sum(y == 0))
    if n_pos == 0 or n_neg == 0:
        return math.nan
    order = np.argsort(s, kind="mergesort")
    ranks = np.empty(len(s), dtype="float64")
    ranks[order] = np.arange(1, len(s) + 1, dtype="float64")
    _, inverse, counts = np.unique(s[order], return_inverse=True, return_counts=True)
    if np.any(counts > 1):
        starts = np.cumsum(np.r_[0, counts[:-1]]) + 1.0
        avg = starts + (counts - 1.0) / 2.0
        ranks[order] = avg[inverse]
    pos_ranks = float(np.sum(ranks[y == 1]))
    return (pos_ranks - n_pos * (n_pos + 1) / 2.0) / float(n_pos * n_neg)


def roc_curve(y_true: np.ndarray, score: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    y = np.asarray(y_true, dtype="int32")
    s = np.asarray(score, dtype="float64")
    mask = np.isfinite(s) & np.isin(y, [0, 1])
    y = y[mask]
    s = s[mask]
    n_pos = int(np.sum(y == 1))
    n_neg = int(np.sum(y == 0))
    if n_pos == 0 or n_neg == 0:
        return np.array([0.0, 1.0]), np.array([0.0, 1.0]), math.nan
    order = np.argsort(-s, kind="mergesort")
    ys = y[order]
    tp = np.cumsum(ys == 1, dtype="float64")
    fp = np.cumsum(ys == 0, dtype="float64")
    tpr = np.concatenate(([0.0], tp / n_pos, [1.0]))
    fpr = np.concatenate(([0.0], fp / n_neg, [1.0]))
    return fpr, tpr, float(np.trapezoid(tpr, fpr))


def threshold_for_eff(y_true: np.ndarray, score: np.ndarray, target: float = 0.80) -> float:
    y = np.asarray(y_true, dtype="int32")
    s = np.asarray(score, dtype="float64")
    sig = s[(y == 1) & np.isfinite(s)]
    if len(sig) == 0:
        return math.nan
    return float(np.quantile(sig, max(0.0, min(1.0, 1.0 - target))))


def calibration_bins(y_true: np.ndarray, score: np.ndarray, bins: np.ndarray) -> tuple[list[dict], float, float]:
    y = np.asarray(y_true, dtype="float64")
    s = np.asarray(score, dtype="float64")
    rows = []
    total = 0
    ece = 0.0
    brier_terms = []
    for lo, hi in zip(bins[:-1], bins[1:], strict=True):
        mask = np.isfinite(s) & (s >= lo) & (s < hi if hi < bins[-1] else s <= hi)
        n = int(mask.sum())
        if n == 0:
            rows.append({"lo": lo, "hi": hi, "entries": 0, "confidence": math.nan, "signal_fraction": math.nan})
            continue
        conf = float(np.mean(s[mask]))
        frac = float(np.mean(y[mask]))
        rows.append({"lo": float(lo), "hi": float(hi), "entries": n, "confidence": conf, "signal_fraction": frac})
        total += n
        ece += n * abs(conf - frac)
        brier_terms.append(np.mean((s[mask] - y[mask]) ** 2) * n)
    return rows, (ece / total if total else math.nan), (float(np.sum(brier_terms) / total) if total else math.nan)


def setup_style():
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "axes.linewidth": 1.1,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )


def label_sphenix(fig, y=0.985):
    fig.text(0.025, y, "sPHENIX", ha="left", va="top", fontsize=15, fontstyle="italic", fontweight="bold")
    fig.text(0.132, y, "Internal", ha="left", va="top", fontsize=15)


def plot_learning_curves(args: argparse.Namespace, artifact: dict, outdir: Path) -> Path:
    import matplotlib.pyplot as plt

    mlp_rows = read_csv_rows(args.mlp_history)
    clean_rows = [r for r in mlp_rows if "High-pT balanced clean MLP" in r.get("label", "")]
    if not clean_rows:
        clean_rows = mlp_rows
    stack_histories = []
    for route in artifact.get("routes", []):
        hist = route.get("model", {}).get("training", {}).get("history", [])
        if hist:
            stack_histories.append(hist)

    fig, axes = plt.subplots(2, 1, figsize=(10.8, 8.4), sharex=False)
    ax = axes[0]
    if clean_rows:
        epochs = np.array([fval(r, "epoch") for r in clean_rows])
        train = np.array([fval(r, "train_loss") for r in clean_rows])
        val = np.array([fval(r, "validation_loss") for r in clean_rows])
        ax.plot(epochs, train, color="#1f77b4", lw=2.0, marker="o", ms=2.8, label="Training loss")
        ax.plot(epochs, val, color="#d95f02", lw=2.0, marker="s", ms=2.8, label="Validation loss")
        best = int(epochs[np.nanargmin(val)])
        ax.axvline(best, color="0.35", lw=1.4, ls=":")
        ax.text(0.98, 0.86, f"best val epoch {best}", transform=ax.transAxes, ha="right", fontsize=12, color="0.25")
    ax.set_title("Clean MLP input model", fontsize=15, fontweight="bold")
    ax.set_ylabel("Loss", fontsize=13)
    ax.grid(True, color="0.90")
    ax.legend(frameon=False, fontsize=12, loc="best")

    ax = axes[1]
    if stack_histories:
        max_epoch = max(max(int(h["epoch"]) for h in hist) for hist in stack_histories)
        epochs = np.arange(1, max_epoch + 1)
        train_matrix = np.full((len(stack_histories), max_epoch), np.nan)
        val_matrix = np.full_like(train_matrix, np.nan)
        auc_matrix = np.full_like(train_matrix, np.nan)
        for i, hist in enumerate(stack_histories):
            for h in hist:
                j = int(h["epoch"]) - 1
                train_matrix[i, j] = float(h["train_loss"])
                val_matrix[i, j] = float(h["val_loss"])
                auc_matrix[i, j] = float(h.get("val_auc", math.nan))
        train_med = np.nanmedian(train_matrix, axis=0)
        val_med = np.nanmedian(val_matrix, axis=0)
        train_lo, train_hi = np.nanpercentile(train_matrix, [16, 84], axis=0)
        val_lo, val_hi = np.nanpercentile(val_matrix, [16, 84], axis=0)
        ax.fill_between(epochs, train_lo, train_hi, color="#1f77b4", alpha=0.16, lw=0)
        ax.fill_between(epochs, val_lo, val_hi, color="#d95f02", alpha=0.16, lw=0)
        ax.plot(epochs, train_med, color="#1f77b4", lw=2.2, label="Median route training loss")
        ax.plot(epochs, val_med, color="#d95f02", lw=2.2, label="Median route validation loss")
        best_epochs = [route["model"]["training"]["best_epoch"] for route in artifact.get("routes", []) if route["model"].get("training")]
        ax.axvline(float(np.median(best_epochs)), color="0.35", lw=1.4, ls=":")
        ax.text(
            0.02,
            0.88,
            f"42 routed NN stackers; median best epoch {int(np.median(best_epochs))}",
            transform=ax.transAxes,
            ha="left",
            fontsize=11.5,
            color="0.25",
        )
        ax2 = ax.twinx()
        ax2.plot(epochs, np.nanmedian(auc_matrix, axis=0), color="#6a3d9a", lw=1.6, ls="--", label="Median route val AUC")
        ax2.set_ylabel("Validation AUC", fontsize=12, color="#6a3d9a")
        ax2.tick_params(axis="y", labelcolor="#6a3d9a")
    ax.set_title("BDT+MLP NN stacker routed submodels", fontsize=15, fontweight="bold")
    ax.set_xlabel("Epoch", fontsize=13)
    ax.set_ylabel("Loss", fontsize=13)
    ax.grid(True, color="0.90")
    ax.legend(frameon=False, fontsize=11, loc="lower left")

    label_sphenix(fig, y=0.965)
    fig.text(0.54, 0.955, "Au+Au MLP / NN-stack training diagnostics", ha="center", va="top", fontsize=20, fontweight="bold")
    fig.text(0.54, 0.895, "Validation loss flattening, route spread shown as 16-84% band", ha="center", va="top", fontsize=12, color="0.35")
    fig.tight_layout(rect=[0.02, 0.03, 1.0, 0.88], h_pad=2.2)
    path = outdir / "01_training_validation_loss_mlp_and_nnstack.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    return path


def split_description(split: str) -> str:
    return {
        "test": "held-out test split",
        "val": "held-out validation split",
        "all": "all uncapped validation rows",
    }.get(split, split)


def plot_roc_and_score(frame: dict, stack_score: np.ndarray, split_mask: np.ndarray, outdir: Path, split: str) -> list[Path]:
    import matplotlib.pyplot as plt

    y = frame["is_signal"].astype("int32")
    et = np.asarray(frame["cluster_Et"], dtype="float64")
    cent = np.asarray(frame["centrality"], dtype="float64")
    base = split_mask & (et >= 15.0) & (et < 35.0) & np.isfinite(cent)
    cent_bins = [(0, 20), (20, 50), (50, 80)]
    scores = {
        "BDT": np.asarray(frame["bdt_score"], dtype="float64"),
        "MLP": np.asarray(frame["mlp_score"], dtype="float64"),
        "BDT+MLP NN stack": stack_score,
    }
    colors = {"BDT": "#1f77b4", "MLP": "#d95f02", "BDT+MLP NN stack": "#6a3d9a"}
    paths: list[Path] = []

    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.5), sharex=True, sharey=True)
    roc_rows = []
    for ax, (lo, hi) in zip(axes, cent_bins, strict=True):
        mask = base & (cent >= lo) & (cent < hi)
        for name, score in scores.items():
            fpr, tpr, auc = roc_curve(y[mask], score[mask])
            ax.plot(fpr, tpr, color=colors[name], lw=2.25, label=f"{name} AUC={auc:.3f}")
            roc_rows.append({"model": name, "cent_lo": lo, "cent_hi": hi, "auc": auc, "entries": int(mask.sum())})
        ax.plot([0, 1], [0, 1], color="0.70", ls="--", lw=1.0)
        ax.set_title(f"{lo}-{hi}% centrality", fontsize=13, fontweight="bold")
        ax.set_xlabel("False positive rate / background efficiency", fontsize=11)
        ax.grid(True, color="0.90")
        ax.legend(fontsize=8.6, loc="lower right", frameon=True, framealpha=0.92, edgecolor="0.75")
    axes[0].set_ylabel("True positive rate / signal efficiency", fontsize=11)
    label_sphenix(fig, y=0.955)
    fig.text(0.52, 0.945, "ROC curves by centrality", ha="center", va="top", fontsize=18, fontweight="bold")
    fig.text(0.52, 0.875, rf"$15 < E_T < 35$ GeV; {split_description(split)}; exact NN-stack artifact", ha="center", fontsize=11, color="0.35")
    fig.tight_layout(rect=[0.035, 0.06, 1.0, 0.78], w_pad=1.8)
    path = outdir / "02_roc_by_centrality_bdt_mlp_nnstack.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    paths.append(path)

    fig, axes = plt.subplots(3, 3, figsize=(14.5, 8.9), sharex=True)
    bins = np.linspace(0, 1, 30)
    for row, (name, score) in enumerate(scores.items()):
        for col, (lo, hi) in enumerate(cent_bins):
            ax = axes[row, col]
            mask = base & (cent >= lo) & (cent < hi) & np.isfinite(score)
            sig = score[mask & (y == 1)]
            bkg = score[mask & (y == 0)]
            auc = auc_score(y[mask], score[mask])
            ax.hist(sig, bins=bins, histtype="step", density=True, color="#1f77b4", lw=1.7, label="Signal")
            ax.hist(bkg, bins=bins, histtype="step", density=True, color="#b33b4b", lw=1.7, label="Background")
            ax.text(0.05, 0.84, f"AUC = {auc:.3f}", transform=ax.transAxes, fontsize=11)
            if row == 0:
                ax.set_title(f"{lo}-{hi}% centrality", fontsize=12.5, fontweight="bold")
            if col == 0:
                ax.set_ylabel(f"{name}\nDensity", fontsize=10.5, fontweight="bold")
            if row == 2:
                ax.set_xlabel("Classifier score", fontsize=11)
            ax.grid(True, color="0.93", lw=0.6)
            if row == 0 and col == 2:
                ax.legend(frameon=False, fontsize=9)
    label_sphenix(fig, y=0.965)
    fig.text(0.54, 0.955, r"Signal/background score separation, $15 < E_T < 35$ GeV", ha="center", va="top", fontsize=18, fontweight="bold")
    fig.text(0.54, 0.895, f"{split_description(split)}; exact score definitions as ROC overlay", ha="center", fontsize=11, color="0.35")
    fig.tight_layout(rect=[0.04, 0.05, 1.0, 0.84], h_pad=1.3, w_pad=1.1)
    path = outdir / "03_score_separation_3x3_bdt_mlp_nnstack.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    paths.append(path)

    with (outdir / "roc_score_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["model", "cent_lo", "cent_hi", "auc", "entries"])
        writer.writeheader()
        writer.writerows(roc_rows)
    return paths


def plot_calibration_and_stability(frame: dict, stack_score: np.ndarray, split_mask: np.ndarray, outdir: Path, split: str) -> list[Path]:
    import matplotlib.pyplot as plt

    y = frame["is_signal"].astype("int32")
    et = np.asarray(frame["cluster_Et"], dtype="float64")
    cent = np.asarray(frame["centrality"], dtype="float64")
    eiso = np.asarray(frame["reco_eiso"], dtype="float64")
    base = split_mask & (et >= 15.0) & (et < 35.0)
    scores = {
        "MLP": np.asarray(frame["mlp_score"], dtype="float64"),
        "BDT+MLP NN stack": stack_score,
    }
    colors = {"MLP": "#d95f02", "BDT+MLP NN stack": "#6a3d9a"}
    paths = []

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.5), sharex=True, sharey=True)
    bins = np.linspace(0, 1, 11)
    calib_rows = []
    for ax, (name, score) in zip(axes, scores.items(), strict=True):
        rows, ece, brier = calibration_bins(y[base], score[base], bins)
        x = [r["confidence"] for r in rows if r["entries"] > 0]
        yy = [r["signal_fraction"] for r in rows if r["entries"] > 0]
        sizes = [min(1350, max(24, 4.0 * math.sqrt(r["entries"]))) for r in rows if r["entries"] > 0]
        ax.scatter(x, yy, s=sizes, color=colors[name], alpha=0.85, edgecolor="white", lw=0.7)
        ax.plot([0, 1], [0, 1], color="0.55", ls="--", lw=1.2)
        ax.set_title(f"{name}\nECE={ece:.3f}, Brier={brier:.3f}", fontsize=13, fontweight="bold")
        ax.set_xlabel("Mean raw classifier score", fontsize=11)
        ax.grid(True, color="0.90")
        for r in rows:
            calib_rows.append({"model": name, "ece": ece, "brier": brier, **r})
    axes[0].set_ylabel("Observed signal fraction", fontsize=11)
    label_sphenix(fig, y=0.955)
    fig.text(0.54, 0.945, "Reliability / calibration check", ha="center", va="top", fontsize=18, fontweight="bold")
    fig.text(
        0.54,
        0.875,
        rf"$15 < E_T < 35$ GeV; {split_description(split)}; raw-score calibration diagnostic",
        ha="center",
        fontsize=11,
        color="0.35",
    )
    fig.tight_layout(rect=[0.04, 0.06, 1.0, 0.78], w_pad=2.0)
    path = outdir / "04_reliability_calibration_mlp_vs_nnstack.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    paths.append(path)

    with (outdir / "calibration_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(calib_rows[0].keys()))
        writer.writeheader()
        writer.writerows(calib_rows)

    et_bins = np.array([15, 18, 20, 22.5, 25, 30, 35], dtype="float64")
    cent_bins = np.array([0, 10, 20, 30, 40, 50, 60, 80], dtype="float64")
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.3), sharex="col")
    stability_rows = []
    for col, (name, score) in enumerate(scores.items()):
        auc_et = []
        fake_et = []
        xcent = []
        auc_cent = []
        for lo, hi in zip(et_bins[:-1], et_bins[1:], strict=True):
            mask = base & (et >= lo) & (et < hi)
            thr = threshold_for_eff(y[mask], score[mask], 0.80)
            fake = float(np.mean(score[mask & (y == 0)] >= thr)) if math.isfinite(thr) and np.any(mask & (y == 0)) else math.nan
            auc = auc_score(y[mask], score[mask])
            auc_et.append(auc)
            fake_et.append(fake)
            stability_rows.append({"model": name, "axis": "et", "lo": lo, "hi": hi, "auc": auc, "wp80_fake": fake})
        for lo, hi in zip(cent_bins[:-1], cent_bins[1:], strict=True):
            mask = base & (cent >= lo) & (cent < hi)
            auc = auc_score(y[mask], score[mask])
            xcent.append(0.5 * (lo + hi))
            auc_cent.append(auc)
            stability_rows.append({"model": name, "axis": "centrality", "lo": lo, "hi": hi, "auc": auc, "wp80_fake": math.nan})
        centers = 0.5 * (et_bins[:-1] + et_bins[1:])
        axes[0, col].plot(centers, auc_et, marker="o", color=colors[name], lw=2.1)
        axes[1, col].plot(centers, fake_et, marker="s", color=colors[name], lw=2.1)
        axes[0, col].set_title(name, fontsize=13.5, fontweight="bold")
        axes[0, col].set_ylabel("AUC", fontsize=11)
        axes[1, col].set_ylabel("WP80 background efficiency", fontsize=11)
        axes[1, col].set_xlabel(r"$E_T$ bin center [GeV]", fontsize=11)
        for ax in axes[:, col]:
            ax.grid(True, color="0.90")
    label_sphenix(fig, y=0.965)
    fig.text(0.54, 0.955, r"Performance stability versus $E_T$", ha="center", va="top", fontsize=18, fontweight="bold")
    fig.text(0.54, 0.895, f"{split_description(split)}; WP80 thresholds are re-derived within each plotted bin", ha="center", fontsize=11, color="0.35")
    fig.tight_layout(rect=[0.04, 0.06, 1.0, 0.84], h_pad=1.5)
    path = outdir / "05_auc_and_wp80_fake_vs_et.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    paths.append(path)

    with (outdir / "stability_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(stability_rows[0].keys()))
        writer.writeheader()
        writer.writerows(stability_rows)

    eiso_bins = np.quantile(eiso[base & np.isfinite(eiso)], np.linspace(0, 1, 9))
    eiso_bins = np.unique(eiso_bins)
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.4), sharey=True)
    for ax, (name, score) in zip(axes, scores.items(), strict=True):
        xs, sig_mean, bkg_mean = [], [], []
        for lo, hi in zip(eiso_bins[:-1], eiso_bins[1:], strict=True):
            mask = base & np.isfinite(eiso) & (eiso >= lo) & (eiso <= hi)
            xs.append(float(np.nanmedian(eiso[mask])))
            sig_mean.append(float(np.nanmean(score[mask & (y == 1)])))
            bkg_mean.append(float(np.nanmean(score[mask & (y == 0)])))
        ax.plot(xs, sig_mean, color="#1f77b4", marker="o", lw=2.0, label="Signal")
        ax.plot(xs, bkg_mean, color="#b33b4b", marker="s", lw=2.0, label="Background")
        corr = finite_corr(score[base], eiso[base])
        ax.set_title(f"{name}\nscore-isolation corr. {corr:.3f}", fontsize=13, fontweight="bold")
        ax.set_xlabel("Isolation quantile-bin median", fontsize=11)
        ax.grid(True, color="0.90")
        ax.legend(frameon=False, fontsize=10)
    axes[0].set_ylabel("Mean classifier score", fontsize=11)
    label_sphenix(fig, y=0.955)
    fig.text(0.54, 0.945, "Score response versus isolation", ha="center", va="top", fontsize=18, fontweight="bold")
    fig.text(0.54, 0.875, f"{split_description(split)}; isolation is not used as an input here", ha="center", fontsize=11, color="0.35")
    fig.tight_layout(rect=[0.04, 0.06, 1.0, 0.78], w_pad=2.0)
    path = outdir / "06_score_vs_isolation_mlp_vs_nnstack.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    paths.append(path)
    return paths


def plot_rank_and_wp(args: argparse.Namespace, outdir: Path) -> list[Path]:
    import matplotlib.pyplot as plt

    paths = []
    if args.stack_rank.exists():
        rows = [r for r in read_csv_rows(args.stack_rank) if r.get("split") == "test"]
        rows = sorted(rows, key=lambda r: fval(r, "wp80_fake"))
        names = [r["model"].replace("ptFine15to35_cent7_full_", "NN stack ") for r in rows[:6]]
        auc = [fval(r, "auc") for r in rows[:6]]
        fake = [fval(r, "wp80_fake") for r in rows[:6]]
        fig, axes = plt.subplots(1, 2, figsize=(12.8, 4.6))
        y = np.arange(len(names))
        axes[0].barh(y, auc, color="#6a3d9a")
        axes[0].set_yticks(y, names)
        axes[0].invert_yaxis()
        axes[0].set_xlabel("Held-out AUC", fontsize=11)
        axes[0].grid(True, axis="x", color="0.90")
        axes[1].barh(y, fake, color="#c44e52")
        axes[1].set_yticks(y, [])
        axes[1].invert_yaxis()
        axes[1].set_xlabel("WP80 background efficiency", fontsize=11)
        axes[1].grid(True, axis="x", color="0.90")
        label_sphenix(fig, y=0.955)
        fig.text(0.55, 0.945, "BDT+MLP stacker held-out ranking", ha="center", va="top", fontsize=18, fontweight="bold")
        fig.text(0.55, 0.875, "Lower WP80 background efficiency is better; stackers use BDT score as an input", ha="center", fontsize=11, color="0.35")
        fig.tight_layout(rect=[0.06, 0.08, 1.0, 0.78])
        path = outdir / "07_stacker_rank_auc_wp80.png"
        fig.savefig(path, dpi=220)
        plt.close(fig)
        paths.append(path)

    for name in [
        "stack_wp80_threshold_grid.png",
        "stack_wp80_signal_efficiency_grid.png",
        "stack_wp80_fake_rate_grid.png",
        "stack_wp80_quadratic_fit_by_fine_centrality.png",
    ]:
        src = args.wp_diag_dir / name
        if src.exists():
            dst = outdir / f"08_{name}"
            shutil.copy2(src, dst)
            paths.append(dst)
    return paths


def heatmap_text_color(value: float, vmin: float, vmax: float, cmap_name: str = "viridis") -> str:
    if not math.isfinite(value) or vmax <= vmin:
        return "black"
    import matplotlib.pyplot as plt

    rgba = plt.get_cmap(cmap_name)((value - vmin) / (vmax - vmin))
    luminance = 0.2126 * rgba[0] + 0.7152 * rgba[1] + 0.0722 * rgba[2]
    return "black" if luminance > 0.58 else "white"


def plot_pt_cent_heatmaps(frame: dict, stack_score: np.ndarray, split_mask: np.ndarray, outdir: Path, split: str) -> list[Path]:
    import matplotlib.pyplot as plt

    y = frame["is_signal"].astype("int32")
    et = np.asarray(frame["cluster_Et"], dtype="float64")
    cent = np.asarray(frame["centrality"], dtype="float64")
    base = split_mask & (et >= 15.0) & (et < 35.0) & np.isfinite(cent)
    et_bins = np.array([15, 18, 20, 22.5, 25, 30, 35], dtype="float64")
    cent_bins = np.array([0, 10, 20, 30, 40, 50, 60, 80], dtype="float64")
    models = {
        "MLP": np.asarray(frame["mlp_score"], dtype="float64"),
        "BDT+MLP NN stack": stack_score,
    }
    rows = []
    auc_maps = {}
    fake_maps = {}
    entry_maps = {}
    for name, score in models.items():
        auc_grid = np.full((len(cent_bins) - 1, len(et_bins) - 1), np.nan)
        fake_grid = np.full_like(auc_grid, np.nan)
        entry_grid = np.zeros_like(auc_grid)
        for iy, (clo, chi) in enumerate(zip(cent_bins[:-1], cent_bins[1:], strict=True)):
            for ix, (elo, ehi) in enumerate(zip(et_bins[:-1], et_bins[1:], strict=True)):
                mask = base & (cent >= clo) & (cent < chi) & (et >= elo) & (et < ehi)
                auc = auc_score(y[mask], score[mask])
                thr = threshold_for_eff(y[mask], score[mask], 0.80)
                fake = float(np.mean(score[mask & (y == 0)] >= thr)) if math.isfinite(thr) and np.any(mask & (y == 0)) else math.nan
                auc_grid[iy, ix] = auc
                fake_grid[iy, ix] = fake
                entry_grid[iy, ix] = int(mask.sum())
                rows.append(
                    {
                        "model": name,
                        "et_lo": elo,
                        "et_hi": ehi,
                        "cent_lo": clo,
                        "cent_hi": chi,
                        "entries": int(mask.sum()),
                        "signal_entries": int(np.sum(mask & (y == 1))),
                        "background_entries": int(np.sum(mask & (y == 0))),
                        "auc": auc,
                        "wp80_fake": fake,
                    }
                )
        auc_maps[name] = auc_grid
        fake_maps[name] = fake_grid
        entry_maps[name] = entry_grid

    with (outdir / "pt_cent_heatmap_summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    paths = []
    for metric_name, maps, cmap, vmin, vmax, label, subtitle in [
        (
            "auc",
            auc_maps,
            "viridis",
            0.60,
            0.92,
            "AUC",
            f"{split_description(split)}; AUC computed independently in each centrality x $E_T$ cell",
        ),
        (
            "wp80_fake",
            fake_maps,
            "magma_r",
            0.10,
            0.65,
            "WP80 background efficiency",
            f"{split_description(split)}; signal-efficiency threshold re-derived in each cell",
        ),
    ]:
        fig = plt.figure(figsize=(14.2, 5.7))
        grid_spec = fig.add_gridspec(
            1,
            3,
            width_ratios=[1.0, 1.0, 0.045],
            left=0.075,
            right=0.92,
            bottom=0.16,
            top=0.76,
            wspace=0.16,
        )
        axes = [fig.add_subplot(grid_spec[0, 0]), fig.add_subplot(grid_spec[0, 1])]
        color_axis = fig.add_subplot(grid_spec[0, 2])
        for ax, (name, grid) in zip(axes, maps.items(), strict=True):
            im = ax.imshow(grid, origin="lower", aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax)
            ax.set_title(name, fontsize=13.0, fontweight="bold", pad=8)
            ax.set_xticks(np.arange(len(et_bins) - 1), [f"{lo:g}-{hi:g}" for lo, hi in zip(et_bins[:-1], et_bins[1:], strict=True)], rotation=35, ha="right")
            ax.set_yticks(np.arange(len(cent_bins) - 1), [f"{lo:g}-{hi:g}" for lo, hi in zip(cent_bins[:-1], cent_bins[1:], strict=True)])
            ax.set_xlabel(r"$E_T$ bin [GeV]", fontsize=11)
            ax.tick_params(axis="both", labelsize=10)
            for iy in range(grid.shape[0]):
                for ix in range(grid.shape[1]):
                    value = grid[iy, ix]
                    if math.isfinite(value):
                        ax.text(ix, iy, f"{value:.2f}", ha="center", va="center", fontsize=8.5, color=heatmap_text_color(value, vmin, vmax, cmap))
        axes[0].set_ylabel("Centrality bin [%]", fontsize=11)
        axes[1].tick_params(labelleft=False)
        cbar = fig.colorbar(im, cax=color_axis)
        cbar.set_label(label, fontsize=11)
        label_sphenix(fig, y=0.965)
        fig.text(0.56, 0.955, f"{label} by centrality and $E_T$", ha="center", va="top", fontsize=18, fontweight="bold")
        fig.text(0.56, 0.895, subtitle, ha="center", fontsize=11, color="0.35")
        path = outdir / f"07_{metric_name}_heatmap_pt_cent_mlp_vs_nnstack.png"
        fig.savefig(path, dpi=220)
        plt.close(fig)
        paths.append(path)
    return paths


def write_manifest(outdir: Path, args: argparse.Namespace, plots: list[Path], cache_info: dict, summary: dict) -> None:
    payload = {
        "schema": "AUAU_MLP_EXPERT_VALIDATION_PACK_V1",
        "outdir": str(outdir),
        "plots": [str(p) for p in plots],
        "inputs": {
            "mlp_cache": str(args.mlp_cache),
            "bdt_cache": str(args.bdt_cache),
            "stack_artifact": str(args.stack_artifact),
            "mlp_history": str(args.mlp_history),
            "mlp_metrics": str(args.mlp_metrics),
            "stack_rank": str(args.stack_rank),
            "wp_cells": str(args.wp_cells),
            "wp_diag_dir": str(args.wp_diag_dir),
        },
        "cache_info": cache_info,
        "summary": summary,
        "caveat": "BDT+MLP stacker uses BDT score as a runtime input; treat as diagnostic/closure-pending.",
    }
    (outdir / "expert_validation_manifest.json").write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def stack_raw_feature_names(artifact: dict) -> set[str]:
    names: set[str] = set()
    if "model" in artifact:
        names.update(artifact["model"].get("feature_names", []))
    for route in artifact.get("routes", []):
        names.update(route.get("model", {}).get("feature_names", []))
    names.difference_update(DERIVED_STACK_FEATURES)
    names.difference_update({"bdt_score", "mlp_score"})
    names.update(["is_signal", "cluster_Et", "centrality", "reco_eiso"])
    return names


def load_stream_scored_frame(args: argparse.Namespace, artifact: dict) -> tuple[dict, dict]:
    mlp_paths = stack.read_manifest(args.mlp_cache)
    bdt_paths = stack.read_manifest(args.bdt_cache)
    if len(mlp_paths) != len(bdt_paths):
        raise SystemExit(f"Cache manifest length mismatch: MLP {len(mlp_paths)} vs BDT {len(bdt_paths)}")
    if args.max_shards > 0:
        mlp_paths = mlp_paths[: args.max_shards]
        bdt_paths = bdt_paths[: args.max_shards]

    wanted = stack_raw_feature_names(artifact)
    pieces: dict[str, list[np.ndarray]] = {key: [] for key in ["is_signal", "cluster_Et", "centrality", "reco_eiso", "mlp_score", "bdt_score", "stack_score"]}
    mlp_key = args.mlp_score
    bdt_key = args.bdt_score
    missing_seen: set[str] = set()
    rows = 0
    for idx, (mlp_path, bdt_path) in enumerate(zip(mlp_paths, bdt_paths, strict=True), 1):
        if idx == 1 or idx % 10 == 0 or idx == len(mlp_paths):
            print(f"[expertValidation] stream-scoring shard {idx}/{len(mlp_paths)}", flush=True)
        with np.load(mlp_path, allow_pickle=True) as mlp_cache, np.load(bdt_path, allow_pickle=True) as bdt_cache:
            if idx == 1:
                mlp_key = stack.auto_score_key(mlp_cache, args.mlp_score, "MLP cache")
                bdt_key = stack.auto_score_key(bdt_cache, args.bdt_score, "BDT cache")
                print(f"[expertValidation] mlp_score_key={mlp_key}", flush=True)
                print(f"[expertValidation] bdt_score_key={bdt_key}", flush=True)
            stack.compare_alignment(mlp_cache, bdt_cache, f"shard {idx}")
            frame = {}
            for key in wanted:
                if key in mlp_cache.files:
                    frame[key] = np.asarray(mlp_cache[key])
                elif key in bdt_cache.files:
                    frame[key] = np.asarray(bdt_cache[key])
                elif key not in DERIVED_STACK_FEATURES:
                    missing_seen.add(key)
            if missing_seen:
                raise SystemExit(f"Missing required raw stack features: {sorted(missing_seen)[:12]}")
            frame["mlp_score"] = np.asarray(mlp_cache[mlp_key], dtype="float64")
            frame["bdt_score"] = np.asarray(bdt_cache[bdt_key], dtype="float64")
            frame["is_signal"] = np.asarray(frame["is_signal"], dtype="int32")
            stack.add_derived_features(frame)
            score = route_score(artifact, frame)
            for key in ["is_signal", "cluster_Et", "centrality", "reco_eiso", "mlp_score", "bdt_score"]:
                pieces[key].append(np.asarray(frame[key]))
            pieces["stack_score"].append(np.asarray(score, dtype="float64"))
            rows += len(frame["is_signal"])
    frame_out = {key: np.concatenate(parts) for key, parts in pieces.items()}
    if args.max_rows > 0 and len(frame_out["is_signal"]) > args.max_rows:
        rng = np.random.default_rng(args.random_seed)
        chosen = rng.choice(len(frame_out["is_signal"]), size=args.max_rows, replace=False)
        chosen.sort()
        frame_out = {key: value[chosen] for key, value in frame_out.items()}
    return frame_out, {
        "mlp_cache": str(args.mlp_cache),
        "bdt_cache": str(args.bdt_cache),
        "mlp_score_key": mlp_key,
        "bdt_score_key": bdt_key,
        "shards": len(mlp_paths),
        "rows": rows,
        "stream_scoring": True,
        "columns": sorted(frame_out.keys()),
    }


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    ap.add_argument("--mlp-cache", type=Path, default=DEFAULT_MLP_CACHE)
    ap.add_argument("--bdt-cache", type=Path, default=DEFAULT_BDT_CACHE)
    ap.add_argument("--stack-artifact", type=Path, default=DEFAULT_STACK_ARTIFACT)
    ap.add_argument("--mlp-history", type=Path, default=DEFAULT_MLP_HISTORY)
    ap.add_argument("--mlp-metrics", type=Path, default=DEFAULT_MLP_METRICS)
    ap.add_argument("--stack-rank", type=Path, default=DEFAULT_STACK_RANK)
    ap.add_argument("--wp-cells", type=Path, default=DEFAULT_WP_CELLS)
    ap.add_argument("--wp-diag-dir", type=Path, default=DEFAULT_WP_DIAG_DIR)
    ap.add_argument("--split", choices=["test", "val", "all"], default="test")
    ap.add_argument("--random-seed", type=int, default=24681357)
    ap.add_argument("--max-shards", type=int, default=0)
    ap.add_argument("--max-rows", type=int, default=0)
    ap.add_argument("--stream-scoring", action="store_true", help="Score stack shard-by-shard and keep only minimal arrays.")
    return ap.parse_args()


def main() -> None:
    import matplotlib

    matplotlib.use("Agg")
    setup_style()
    args = parse_args()
    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    artifact = json.loads(args.stack_artifact.read_text())
    load_args = argparse.Namespace(
        mlp_cache=args.mlp_cache,
        bdt_cache=args.bdt_cache,
        mlp_score=stack.DEFAULT_MLP_SCORE,
        bdt_score=stack.DEFAULT_BDT_SCORE,
        max_shards=args.max_shards,
        max_rows=args.max_rows,
        random_seed=args.random_seed,
    )
    if args.stream_scoring:
        frame, cache_info = load_stream_scored_frame(load_args, artifact)
        stack_score = np.asarray(frame.pop("stack_score"), dtype="float64")
    else:
        frame, cache_info = stack.load_aligned_caches(load_args)
        stack_score = route_score(artifact, frame)
    y = frame["is_signal"].astype("int32")
    masks = stack.stratified_split(y, args.random_seed, 0.60, 0.20)
    split_mask = np.ones(len(y), dtype=bool) if args.split == "all" else masks[args.split]
    plots: list[Path] = []
    plots.append(plot_learning_curves(args, artifact, outdir))
    plots.extend(plot_roc_and_score(frame, stack_score, split_mask, outdir, args.split))
    plots.extend(plot_calibration_and_stability(frame, stack_score, split_mask, outdir, args.split))
    plots.extend(plot_pt_cent_heatmaps(frame, stack_score, split_mask, outdir, args.split))
    plots.extend(plot_rank_and_wp(args, outdir))
    summary = {
        "rows_loaded": int(len(y)),
        "split": args.split,
        "split_rows": int(split_mask.sum()),
        "mlp_auc_15_35_split": auc_score(
            y[split_mask & (frame["cluster_Et"] >= 15.0) & (frame["cluster_Et"] < 35.0)],
            frame["mlp_score"][split_mask & (frame["cluster_Et"] >= 15.0) & (frame["cluster_Et"] < 35.0)],
        ),
        "stack_auc_15_35_split": auc_score(
            y[split_mask & (frame["cluster_Et"] >= 15.0) & (frame["cluster_Et"] < 35.0)],
            stack_score[split_mask & (frame["cluster_Et"] >= 15.0) & (frame["cluster_Et"] < 35.0)],
        ),
    }
    write_manifest(outdir, args, plots, cache_info, summary)
    print(f"[expertValidation] outdir={outdir}", flush=True)
    for path in plots:
        print(f"[expertValidation] plot={path}", flush=True)
    print(f"[expertValidation] manifest={outdir / 'expert_validation_manifest.json'}", flush=True)


if __name__ == "__main__":
    main()
