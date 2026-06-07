#!/usr/bin/env python3
"""Render slide-ready PNGs for the fresh pp/AuAu BDT+MLP OOF stack campaign."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch


W, H = 2560, 1440
DPI = 200
TITLE = 36
SUBTITLE = 21
BODY = 19
SMALL = 16
COLORS = {
    "ink": "#171717",
    "muted": "#5f6368",
    "blue": "#2f6f9f",
    "green": "#2e7d5b",
    "purple": "#6b4ea3",
    "gold": "#d6a21d",
    "gray": "#eef1f4",
    "line": "#d6d9de",
    "warn": "#f8e7a3",
}


def setup():
    plt.rcParams.update(
        {
            "font.family": "Times New Roman",
            "axes.edgecolor": COLORS["line"],
            "axes.labelcolor": COLORS["ink"],
            "xtick.color": COLORS["ink"],
            "ytick.color": COLORS["ink"],
            "axes.titlecolor": COLORS["ink"],
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )


def fig_ax():
    fig = plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    return fig, ax


def add_text(ax, x, y, text, size=BODY, weight="normal", color=None, ha="left", va="top", **kwargs):
    return ax.text(x, y, text, transform=ax.transAxes, fontsize=size, fontweight=weight, color=color or COLORS["ink"], ha=ha, va=va, **kwargs)


def box(ax, x, y, w, h, fc, ec=None, lw=1.2, radius=0.012):
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle=f"round,pad=0.008,rounding_size={radius}",
        transform=ax.transAxes,
        linewidth=lw,
        edgecolor=ec or fc,
        facecolor=fc,
    )
    ax.add_patch(patch)
    return patch


def write_json(path: Path, payload) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def save(fig, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=DPI)
    plt.close(fig)


def read_json(path: Path) -> dict:
    return json.loads(path.read_text())


def load_domain(root: Path) -> dict:
    manifest = read_json(root / "campaign_manifest.json")
    metrics = pd.read_csv(root / "model_metrics.csv")
    strat = pd.read_csv(root / "stratified_metrics.csv")
    corr = read_json(root / "score_correlations.json")
    hist = pd.read_csv(root / "overlay_histograms.csv")
    return {"root": root, "manifest": manifest, "metrics": metrics, "strat": strat, "corr": corr, "hist": hist}


def locked_metrics(domain: dict) -> pd.DataFrame:
    df = domain["metrics"]
    return df[df["region"] == "locked_test"].copy()


def best_stack_name(domain: dict) -> str:
    df = locked_metrics(domain)
    stack = df[df["model"].str.contains("score_", regex=False)].copy()
    if stack.empty:
        return "MLP" if float(df[df["model"] == "MLP"]["weighted_auc"].iloc[0]) > float(df[df["model"] == "BDT"]["weighted_auc"].iloc[0]) else "BDT"
    stack = stack.sort_values("weighted_auc", ascending=False)
    return str(stack.iloc[0]["model"])


def metric_value(domain: dict, model: str, column: str = "weighted_auc") -> float:
    df = locked_metrics(domain)
    row = df[df["model"] == model]
    if row.empty:
        return math.nan
    return float(row.iloc[0][column])


def model_label(name: str) -> str:
    return {
        "BDT": "BDT",
        "MLP": "MLP",
        "score_only_logistic": "Logistic stack",
        "score_only_gbm": "GBM stack",
        "score_only_mlp": "MLP stack",
        "score_context_logistic": "Logistic + context",
        "score_context_gbm": "GBM + context",
        "score_context_mlp": "MLP + context",
    }.get(name, name.replace("_", " "))


def draw_metric_bars(fig, domain: dict, rect, title: str, emphasize: str | None = None):
    ax = fig.add_axes(rect)
    df = locked_metrics(domain)
    order = ["BDT", "MLP", "score_only_logistic", "score_only_gbm", "score_only_mlp", "score_context_logistic", "score_context_gbm", "score_context_mlp"]
    df = df[df["model"].isin(order)].copy()
    df["order"] = df["model"].map({name: i for i, name in enumerate(order)})
    df = df.sort_values("order")
    vals = df["weighted_auc"].to_numpy(dtype=float)
    colors = [COLORS["blue"] if name == "BDT" else COLORS["green"] if name == "MLP" else COLORS["purple"] for name in df["model"]]
    if emphasize:
        colors = [COLORS["gold"] if name == emphasize else c for name, c in zip(df["model"], colors)]
    ax.bar(np.arange(len(df)), vals, color=colors, width=0.72)
    ax.set_title(title, fontsize=SUBTITLE, pad=14)
    ax.set_ylim(max(0.0, np.nanmin(vals) - 0.04), min(1.0, np.nanmax(vals) + 0.025))
    ax.set_ylabel("Weighted AUC", fontsize=SMALL)
    ax.set_xticks(np.arange(len(df)))
    ax.set_xticklabels([model_label(x) for x in df["model"]], rotation=26, ha="right", fontsize=14)
    ax.tick_params(axis="y", labelsize=SMALL)
    ax.grid(axis="y", color=COLORS["line"], linewidth=0.8)
    for i, v in enumerate(vals):
        ax.text(i, v + 0.002, f"{v:.3f}", ha="center", va="bottom", fontsize=16)
    return ax


def draw_contract_slide(pp, auau, outdir: Path):
    fig, ax = fig_ax()
    add_text(ax, 0.045, 0.94, "Fresh OOF Stack Contract", TITLE, "bold")
    add_text(ax, 0.045, 0.875, "One event-level split controls BDT, MLP, and stack training for pp and Au+Au.", SUBTITLE, color=COLORS["muted"])
    box(ax, 0.055, 0.66, 0.43, 0.20, COLORS["gray"])
    add_text(ax, 0.075, 0.82, "Training classes", 27, "bold")
    add_text(ax, 0.075, 0.77, "Signal = is_signal == 1\nBackground = is_signal == 0\nsource_sample records provenance only", BODY, linespacing=1.28)
    box(ax, 0.515, 0.66, 0.43, 0.20, "#eee8f7")
    add_text(ax, 0.535, 0.82, "Overlay diagnostic", 27, "bold")
    add_text(ax, 0.535, 0.77, "Signal MC = Photon+Jet source, is_signal == 1\nInclusive MC = inclusive-jet source\nwith no truth filter", BODY, linespacing=1.22)
    for x, domain, label in [(0.055, pp, "pp Current-IAN"), (0.515, auau, "Au+Au embedded")]:
        qa = domain["manifest"]["partition_qa"]
        box(ax, x, 0.35, 0.43, 0.27, "white", COLORS["line"])
        add_text(ax, x + 0.02, 0.585, label, 29, "bold")
        train = qa["partitions"]["trainval"]["label_counts"]
        test = qa["partitions"]["locked_test"]["label_counts"]
        add_text(
            ax,
            x + 0.02,
            0.535,
            f"Event key: source_sample/run/evt\nLocked test: {qa['partitions']['locked_test']['events']:,} events\nOOF folds: {qa['folds']} folds over train region\nTest S/B labels: {test['is_signal_1']:,} / {test['is_signal_0']:,}\nTrain S/B labels: {train['is_signal_1']:,} / {train['is_signal_0']:,}",
            SMALL,
            linespacing=1.22,
        )
    box(ax, 0.055, 0.12, 0.89, 0.15, COLORS["warn"], "#d8bf5f")
    add_text(ax, 0.075, 0.235, "Leakage control", 27, "bold")
    add_text(ax, 0.075, 0.19, "Stackers train on fold-held-out BDT/MLP scores.\nLocked-test rows are scored only by final base models trained on non-test events.", BODY, linespacing=1.25)
    path = outdir / "fresh_oof_stack_slide01_split_contract.png"
    save(fig, path)
    write_json(path.with_suffix(".json"), metadata_payload("split_contract", [pp, auau], path))
    (path.with_suffix(".speaker.md")).write_text("Training labels are truth labels. Source samples are provenance for QA and the separate Signal-MC versus Inclusive-MC overlay.\n")
    return path


def draw_pp_slide(pp, outdir: Path):
    best = best_stack_name(pp)
    fig, ax = fig_ax()
    add_text(ax, 0.045, 0.94, "pp Baseline Stack Check", TITLE, "bold")
    add_text(ax, 0.045, 0.875, "PPG12-equivalent Current-IAN sample; locked 20% event test.", SUBTITLE, color=COLORS["muted"])
    draw_metric_bars(fig, pp, [0.08, 0.25, 0.58, 0.50], "Locked-test model comparison", best)
    bdt = metric_value(pp, "BDT")
    mlp = metric_value(pp, "MLP")
    bst = metric_value(pp, best)
    box(ax, 0.69, 0.51, 0.27, 0.30, "#eef6f2", "#b8d6c7")
    add_text(ax, 0.715, 0.765, "Primary read", 28, "bold")
    add_text(ax, 0.715, 0.705, f"BDT AUC = {bdt:.3f}\nMLP AUC = {mlp:.3f}\nBest stack = {model_label(best)}\nStack gain vs BDT = {bst - bdt:+.3f}", BODY, linespacing=1.22)
    box(ax, 0.69, 0.24, 0.27, 0.20, COLORS["gray"], "#ccd2d8")
    add_text(ax, 0.715, 0.395, "Sample wording", 26, "bold")
    add_text(ax, 0.715, 0.35, "Photon+Jet plus inclusive-jet\nCurrent-IAN contract;\nJet40 is overlay-only.", SMALL, linespacing=1.18)
    path = outdir / "fresh_oof_stack_slide02_pp_comparison.png"
    save(fig, path)
    write_json(path.with_suffix(".json"), metadata_payload("pp_comparison", [pp], path))
    (path.with_suffix(".speaker.md")).write_text(f"Compare BDT, MLP, and stackers on the locked pp test set. Best stack gain versus BDT is {bst - bdt:+.3f} weighted AUC.\n")
    return path


def draw_auau_slide(auau, outdir: Path):
    best = best_stack_name(auau)
    fig, ax = fig_ax()
    add_text(ax, 0.045, 0.94, "Au+Au Stack Check", TITLE, "bold")
    add_text(ax, 0.045, 0.875, "Embedded Photon+Jet 12+20 with embedded inclusive Jet 12+20+30+40; locked 20% event test.", SUBTITLE, color=COLORS["muted"])
    draw_metric_bars(fig, auau, [0.08, 0.25, 0.58, 0.50], "0-80% locked-test model comparison", best)
    bdt = metric_value(auau, "BDT")
    mlp = metric_value(auau, "MLP")
    bst = metric_value(auau, best)
    box(ax, 0.69, 0.52, 0.27, 0.29, "#f1eef8", "#cabde1")
    add_text(ax, 0.715, 0.765, "Primary read", 28, "bold")
    add_text(ax, 0.715, 0.705, f"BDT AUC = {bdt:.3f}\nMLP AUC = {mlp:.3f}\nBest stack = {model_label(best)}\nStack gain vs BDT = {bst - bdt:+.3f}", BODY, linespacing=1.22)
    ax2 = fig.add_axes([0.69, 0.23, 0.27, 0.21])
    cent = auau["strat"]
    cent = cent[(cent["stratification"] == "centrality") & (cent["model"].isin(["BDT", best]))].copy()
    for model, color in [("BDT", COLORS["blue"]), (best, COLORS["gold"])]:
        rows = cent[cent["model"] == model].sort_values("bin_lo")
        if rows.empty:
            continue
        centers = 0.5 * (rows["bin_lo"].to_numpy(float) + rows["bin_hi"].to_numpy(float))
        ax2.plot(centers, rows["weighted_auc"].to_numpy(float), marker="o", linewidth=2.2, label=model_label(model), color=color)
    ax2.set_title("Centrality trend", fontsize=SMALL)
    ax2.set_xlabel("Centrality (%)", fontsize=19)
    ax2.set_ylabel("AUC", fontsize=19)
    ax2.tick_params(labelsize=18)
    ax2.grid(color=COLORS["line"])
    ax2.legend(fontsize=17, frameon=False)
    path = outdir / "fresh_oof_stack_slide03_auau_comparison.png"
    save(fig, path)
    write_json(path.with_suffix(".json"), metadata_payload("auau_comparison", [auau], path))
    (path.with_suffix(".speaker.md")).write_text(f"Compare the Au+Au BDT, MLP, and stackers on the locked 0-80% test set. Best stack gain versus BDT is {bst - bdt:+.3f} weighted AUC.\n")
    return path


def draw_et_slide(pp, auau, outdir: Path):
    pp_best = best_stack_name(pp)
    au_best = best_stack_name(auau)
    fig, ax = fig_ax()
    add_text(ax, 0.045, 0.94, "E_T Dependence", TITLE, "bold")
    add_text(ax, 0.045, 0.875, "AUC change relative to the base BDT, evaluated on each locked test slice.", SUBTITLE, color=COLORS["muted"])
    for i, (domain, best, label, color) in enumerate([(pp, pp_best, "pp", COLORS["green"]), (auau, au_best, "Au+Au", COLORS["purple"])]):
        pax = fig.add_axes([0.08 + i * 0.48, 0.25, 0.38, 0.58])
        strat = domain["strat"]
        et = strat[(strat["stratification"] == "Et") & (strat["model"].isin(["BDT", best]))].copy()
        bdt = et[et["model"] == "BDT"].set_index(["bin_lo", "bin_hi"])["weighted_auc"]
        bst = et[et["model"] == best].set_index(["bin_lo", "bin_hi"])["weighted_auc"]
        keys = sorted(set(bdt.index).intersection(set(bst.index)))
        centers = np.asarray([0.5 * (k[0] + k[1]) for k in keys], dtype=float)
        deltas = np.asarray([bst.loc[k] - bdt.loc[k] for k in keys], dtype=float)
        pax.axhline(0.0, color=COLORS["line"], linewidth=1.2)
        pax.plot(centers, deltas, marker="o", linewidth=2.5, color=color)
        pax.set_title(f"{label}: {model_label(best)} - BDT", fontsize=SUBTITLE)
        pax.set_xlabel(r"$E_T$ bin center (GeV)", fontsize=SMALL)
        pax.set_ylabel("AUC change", fontsize=SMALL)
        pax.tick_params(labelsize=SMALL)
        pax.grid(color=COLORS["line"])
    path = outdir / "fresh_oof_stack_slide04_et_dependence.png"
    save(fig, path)
    write_json(path.with_suffix(".json"), metadata_payload("et_dependence", [pp, auau], path))
    (path.with_suffix(".speaker.md")).write_text("This slide asks whether the stack helps in a narrow E_T region rather than only in the inclusive average.\n")
    return path


def draw_heatmap_slide(auau, outdir: Path):
    best = best_stack_name(auau)
    fig, ax = fig_ax()
    add_text(ax, 0.045, 0.94, "Au+Au E_T x Centrality", TITLE, "bold")
    add_text(ax, 0.045, 0.875, f"Cell value is weighted-AUC change for {model_label(best)} relative to BDT.", SUBTITLE, color=COLORS["muted"])
    strat = auau["strat"]
    grid = strat[(strat["stratification"] == "Et_x_centrality") & (strat["model"].isin(["BDT", best]))].copy()
    bdt = grid[grid["model"] == "BDT"].set_index(["et_lo", "et_hi", "cent_lo", "cent_hi"])["weighted_auc"]
    bst = grid[grid["model"] == best].set_index(["et_lo", "et_hi", "cent_lo", "cent_hi"])["weighted_auc"]
    keys = sorted(set(bdt.index).intersection(set(bst.index)))
    et_bins = sorted({(k[0], k[1]) for k in keys})
    cent_bins = sorted({(k[2], k[3]) for k in keys})
    mat = np.full((len(cent_bins), len(et_bins)), np.nan)
    for j, cb in enumerate(cent_bins):
        for i, eb in enumerate(et_bins):
            key = (eb[0], eb[1], cb[0], cb[1])
            if key in bdt.index and key in bst.index:
                mat[j, i] = float(bst.loc[key] - bdt.loc[key])
    hax = fig.add_axes([0.10, 0.22, 0.80, 0.62])
    vmax = max(0.01, float(np.nanmax(np.abs(mat))) if np.isfinite(mat).any() else 0.01)
    im = hax.imshow(mat, cmap="PRGn", vmin=-vmax, vmax=vmax, aspect="auto")
    hax.set_xticks(np.arange(len(et_bins)))
    hax.set_xticklabels([f"{lo:g}-{hi:g}" for lo, hi in et_bins], fontsize=SMALL)
    hax.set_yticks(np.arange(len(cent_bins)))
    hax.set_yticklabels([f"{lo:g}-{hi:g}%" for lo, hi in cent_bins], fontsize=SMALL)
    hax.set_xlabel(r"$E_T$ bin (GeV)", fontsize=BODY)
    hax.set_ylabel("Centrality", fontsize=BODY)
    for j in range(mat.shape[0]):
        for i in range(mat.shape[1]):
            if np.isfinite(mat[j, i]):
                hax.text(i, j, f"{mat[j, i]:+.3f}", ha="center", va="center", fontsize=19, color=COLORS["ink"])
    cax = fig.add_axes([0.92, 0.22, 0.018, 0.62])
    cb = fig.colorbar(im, cax=cax)
    cb.ax.tick_params(labelsize=SMALL)
    cb.set_label("AUC change", fontsize=SMALL)
    path = outdir / "fresh_oof_stack_slide05_auau_et_centrality_heatmap.png"
    save(fig, path)
    write_json(path.with_suffix(".json"), metadata_payload("auau_et_centrality_heatmap", [auau], path))
    (path.with_suffix(".speaker.md")).write_text("This heatmap localizes any stack gain by both E_T and centrality. Small values mean the stack is not adding a qualitatively new handle.\n")
    return path


def draw_correlation_slide(pp, auau, outdir: Path):
    fig, ax = fig_ax()
    add_text(ax, 0.045, 0.94, "BDT vs MLP Independence", TITLE, "bold")
    add_text(ax, 0.045, 0.875, "Score correlations test whether the stack has independent information to combine.", SUBTITLE, color=COLORS["muted"])
    rows = []
    for domain, label in [(pp, "pp"), (auau, "Au+Au")]:
        corr = domain["corr"]["pairs"].get("BDT__vs__MLP", {})
        rows.append((label, corr))
    for i, (label, corr) in enumerate(rows):
        x = 0.075 + i * 0.47
        box(ax, x, 0.54, 0.40, 0.25, "#f7f8fa", COLORS["line"])
        add_text(ax, x + 0.025, 0.745, label, 32, "bold")
        lines = []
        for key, shown in [
            ("all_locked_test", "All locked-test"),
            ("truth_signal_locked_test", "Truth signal"),
            ("truth_background_locked_test", "Truth background"),
            ("overlay_inclusive_mc_locked_test", "Inclusive MC overlay"),
        ]:
            item = corr.get(key, {})
            val = item.get("pearson")
            rows_n = item.get("rows", 0)
            lines.append(f"{shown}: r = {val:.3f}  (n={rows_n:,})" if val is not None and math.isfinite(float(val)) else f"{shown}: n/a")
        add_text(ax, x + 0.025, 0.69, "\n".join(lines), BODY, linespacing=1.28)
    box(ax, 0.075, 0.22, 0.87, 0.15, COLORS["warn"], "#d8bf5f")
    add_text(ax, 0.10, 0.325, "Interpretation", 29, "bold")
    add_text(ax, 0.10, 0.279, "High BDT-MLP correlation means a stack gain should be small; a visible gain needs a low-correlation slice or a different failure mode.", BODY, linespacing=1.22)
    path = outdir / "fresh_oof_stack_slide06_score_correlation.png"
    save(fig, path)
    write_json(path.with_suffix(".json"), metadata_payload("score_correlation", [pp, auau], path))
    (path.with_suffix(".speaker.md")).write_text("Use the correlation numbers to explain whether the stack is combining meaningfully different model responses.\n")
    return path


def metadata_payload(kind: str, domains: list[dict], path: Path) -> dict:
    return {
        "schema": "RJ_FRESH_OOF_STACK_SLIDE_METADATA_V1",
        "slide_kind": kind,
        "png": str(path),
        "training_class_definition": "Signal = is_signal == 1; background = is_signal == 0",
        "overlay_class_definition": "Signal MC = Photon+Jet source with is_signal == 1; Inclusive MC = inclusive-jet source with no truth-background filter",
        "source_artifacts": [str(d["root"]) for d in domains],
        "domain_manifests": [d["manifest"] for d in domains],
        "canvas": {"width": W, "height": H},
    }


def combined_provenance_payload(pp: dict, auau: dict, paths: list[Path]) -> dict:
    def domain_payload(domain: dict) -> dict:
        manifest = domain["manifest"]
        root = domain["root"]
        return {
            "domain": manifest["domain"],
            "root": str(root),
            "manifest": str(root / "campaign_manifest.json"),
            "metrics_csv": str(root / "model_metrics.csv"),
            "metrics_json": str(root / "model_metrics.json"),
            "stratified_metrics_csv": str(root / "stratified_metrics.csv"),
            "overlay_histograms_csv": str(root / "overlay_histograms.csv"),
            "score_correlations_json": str(root / "score_correlations.json"),
            "partition_qa_json": str(root / "partition_qa.json"),
            "leakage_qa_json": str(root / "leakage_qa.json"),
            "feature_contract": manifest.get("feature_contract", {}),
            "sample_definitions": manifest.get("sample_definitions", {}),
            "training_class_definition": manifest.get("training_class_definition"),
            "overlay_class_definition": manifest.get("overlay_class_definition"),
            "partition_seed": manifest.get("partition_seed"),
            "locked_test_fraction": manifest.get("locked_test_fraction"),
            "oof_fold_count": manifest.get("oof_fold_count"),
            "metric_convention": manifest.get("metric_convention", {}),
        }

    return {
        "schema": "RJ_FRESH_OOF_STACK_COMBINED_PROVENANCE_V1",
        "training_class_definition": "Signal = is_signal == 1; background = is_signal == 0",
        "overlay_class_definition": "Signal MC = Photon+Jet source with is_signal == 1; Inclusive MC = inclusive-jet source with no truth-background filter",
        "domains": {
            "pp": domain_payload(pp),
            "auau": domain_payload(auau),
        },
        "slide_pngs": [str(p) for p in paths],
        "slide_metadata_jsons": [str(p.with_suffix(".json")) for p in paths],
        "canvas": {"width": W, "height": H},
    }


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pp-dir", type=Path, required=True)
    ap.add_argument("--auau-dir", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    return ap.parse_args()


def main() -> int:
    setup()
    args = parse_args()
    pp = load_domain(args.pp_dir)
    auau = load_domain(args.auau_dir)
    args.outdir.mkdir(parents=True, exist_ok=True)
    paths = [
        draw_contract_slide(pp, auau, args.outdir),
        draw_pp_slide(pp, args.outdir),
        draw_auau_slide(auau, args.outdir),
        draw_et_slide(pp, auau, args.outdir),
        draw_heatmap_slide(auau, args.outdir),
        draw_correlation_slide(pp, auau, args.outdir),
    ]
    provenance_path = args.outdir / "fresh_oof_stack_combined_provenance.json"
    write_json(provenance_path, combined_provenance_payload(pp, auau, paths))
    write_json(
        args.outdir / "fresh_oof_stack_slide_bundle.json",
        {
            "schema": "RJ_FRESH_OOF_STACK_SLIDE_BUNDLE_V1",
            "slides": [str(p) for p in paths],
            "combined_provenance_json": str(provenance_path),
        },
    )
    print(json.dumps({"status": "READY", "outdir": str(args.outdir), "slides": [str(p) for p in paths], "combined_provenance_json": str(provenance_path)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
