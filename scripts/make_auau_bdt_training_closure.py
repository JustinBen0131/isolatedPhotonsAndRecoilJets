#!/usr/bin/env python3
"""Training-vs-validation BDT score-separation closure plots."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import train_auau_photon_bdt as trainer  # noqa: E402
import validate_auau_tight_bdt_on_sim as validator  # noqa: E402


def parse_edges(text: str) -> list[float]:
    edges = [float(item.strip()) for item in text.split(",") if item.strip()]
    if len(edges) < 2:
        raise SystemExit(f"Need at least two centrality edges, got: {text}")
    if any(hi <= lo for lo, hi in zip(edges[:-1], edges[1:])):
        raise SystemExit(f"Centrality edges must increase: {text}")
    return edges


def auc_rank(y, score) -> float:
    y = np.asarray(y, dtype=np.int32)
    score = np.asarray(score, dtype=np.float64)
    mask = np.isin(y, [0, 1]) & np.isfinite(score)
    y = y[mask]
    score = score[mask]
    n_pos = int((y == 1).sum())
    n_neg = int((y == 0).sum())
    if n_pos == 0 or n_neg == 0:
        return math.nan
    order = np.argsort(score, kind="mergesort")
    ranks = np.empty_like(order, dtype=np.float64)
    ranks[order] = np.arange(1, len(score) + 1, dtype=np.float64)
    # Average ranks for ties.
    sorted_score = score[order]
    start = 0
    while start < len(score):
        stop = start + 1
        while stop < len(score) and sorted_score[stop] == sorted_score[start]:
            stop += 1
        if stop - start > 1:
            ranks[order[start:stop]] = 0.5 * (start + 1 + stop)
        start = stop
    pos_rank_sum = float(ranks[y == 1].sum())
    return (pos_rank_sum - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)


def find_model(registry: dict, model_id: str | None, product: str | None) -> dict:
    models = registry.get("models", [])
    if model_id:
        matches = [m for m in models if m.get("model_id") == model_id]
    elif product:
        matches = [m for m in models if m.get("product") == product]
    else:
        raise SystemExit("Pass --model-id or --product")
    if not matches:
        key = model_id or product
        raise SystemExit(f"No model found in registry for {key!r}")
    ready = [m for m in matches if (m.get("report") or {}).get("status", "trained") == "trained"]
    return ready[0] if ready else matches[0]


def default_training_args(model: dict, registry: dict) -> SimpleNamespace:
    report = model.get("report") or {}
    defaults = registry.get("defaults") or {}
    bkg = report.get("background_subsampling") or {}
    return SimpleNamespace(
        min_rows_per_class=10,
        test_size=0.10,
        random_seed=13,
        use_event_weight=True,
        weight_branch="event_weight",
        et_reweight=True,
        eta_reweight=True,
        flatten_bins=20,
        max_flatten_factor=8.0,
        max_total_weight_factor=50.0,
        majority_cap_ratio=float(report.get("majority_cap_ratio", defaults.get("majority_cap_ratio", 4.0)) or 4.0),
        majority_cap_et_bins=12,
        majority_cap_cent_bins=8,
        background_subsample_fraction=float(bkg.get("fraction", 1.0) if bkg.get("enabled") else 1.0),
        background_subsample_et_threshold=float(bkg.get("et_threshold", 15.0)),
        background_subsample_bins=int(bkg.get("bins", 20) or 20),
        background_subsample_seed=42,
        background_subsample_flatten=bool(bkg.get("enabled", False) and bkg.get("flatten", False)),
    )


def load_training_scores(registry_path: Path, train_cache: Path, model_id: str | None, product: str | None):
    import xgboost as xgb
    from sklearn.model_selection import train_test_split

    registry = json.loads(registry_path.read_text())
    model = find_model(registry, model_id, product)
    report = model.get("report") or {}
    features = list(model.get("features") or report.get("features") or [])
    label = report.get("label_branch") or "is_signal"
    if not features:
        raise SystemExit(f"Model {model.get('model_id')} has no feature list")
    xgb_path = Path(report.get("output_xgb_json") or Path(report.get("output_tmva", "")).with_suffix(".xgb.json"))
    if not xgb_path.is_file():
        raise SystemExit(f"Missing XGBoost JSON model: {xgb_path}")

    frame = trainer.load_frame_cache(train_cache)
    frame = trainer.add_derived_features(frame)
    frame = trainer.frame_for_spec(frame, model)
    cols = features + [label]
    frame = frame.loc[trainer.finite_mask(frame, cols)].copy()
    frame[label] = frame[label].astype(int)
    frame = frame[frame[label].isin([0, 1])].copy()
    closure_args = default_training_args(model, registry)
    frame, _ = trainer.ppg12_background_subsample(frame, label, closure_args)
    frame, _ = trainer.adaptive_majority_cap(frame, label, closure_args, report)

    y = frame[label].to_numpy(dtype="int32")
    weights, _ = trainer.compute_weights(frame, label, closure_args)
    stratify = y if min(int((y == 0).sum()), int((y == 1).sum())) >= 2 else None
    indices = np.arange(len(frame), dtype=np.int64)
    train_idx, test_idx = train_test_split(
        indices,
        test_size=closure_args.test_size,
        random_state=closure_args.random_seed,
        stratify=stratify,
    )

    booster = xgb.Booster()
    booster.load_model(str(xgb_path))
    x = frame[features].to_numpy(dtype="float32")
    score = booster.predict(xgb.DMatrix(x))
    score = np.asarray(score, dtype=np.float64).reshape(-1)
    out = {
        "is_signal": y,
        "score": score,
        "cluster_Et": frame["cluster_Et"].to_numpy(dtype="float64"),
        "centrality": frame["centrality"].to_numpy(dtype="float64"),
        "split": np.full(len(frame), "train", dtype=object),
    }
    out["split"][test_idx] = "internal_test"
    metadata = {
        "registry": str(registry_path),
        "train_cache": str(train_cache),
        "model_id": model.get("model_id"),
        "product": model.get("product"),
        "features": features,
        "xgb_model": str(xgb_path),
        "n_training_closure_rows": int(len(frame)),
        "n_train_split_rows": int(len(train_idx)),
        "n_internal_test_rows": int(len(test_idx)),
    }
    return out, metadata


def load_validation_scores(manifest: Path, product: str, max_rows: int = 0):
    frame, scores, counts = validator.load_score_caches(manifest)
    key = product
    if key not in scores and product.startswith("score_"):
        key = product[6:]
    if key not in scores:
        raise SystemExit(f"Validation caches do not contain score_{product}")
    n = len(frame["is_signal"])
    if max_rows > 0 and n > max_rows:
        rng = np.random.default_rng(24681357)
        idx = np.sort(rng.choice(n, size=max_rows, replace=False))
    else:
        idx = np.arange(n)
    return {
        "is_signal": np.asarray(frame["is_signal"], dtype="int32")[idx],
        "score": np.asarray(scores[key], dtype="float64")[idx],
        "cluster_Et": np.asarray(frame["cluster_Et"], dtype="float64")[idx],
        "centrality": np.asarray(frame["centrality"], dtype="float64")[idx],
        "split": np.full(len(idx), "external_validation", dtype=object),
    }, counts


def density_hist(ax, signal, background, bins: int):
    edges = np.linspace(0.0, 1.0, bins + 1)
    ax.hist(signal, bins=edges, density=True, histtype="step", lw=1.8, color="#1f77b4", label="Signal")
    ax.hist(background, bins=edges, density=True, histtype="step", lw=1.8, color="#b23a48", label="Background")


def plot_and_write(training, validation, metadata, args):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "axes.linewidth": 1.0,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )

    args.outdir.mkdir(parents=True, exist_ok=True)
    cent_edges = parse_edges(args.cent_bins)
    datasets = [
        ("Training split", training, "train"),
        ("Internal test split", training, "internal_test"),
        ("Independent validation", validation, "external_validation"),
    ]
    fig, axes = plt.subplots(len(datasets), len(cent_edges) - 1, figsize=(14.2, 8.4), sharex=True)
    summary = []
    for irow, (row_label, data, split_name) in enumerate(datasets):
        y = np.asarray(data["is_signal"], dtype=np.int32)
        score = np.asarray(data["score"], dtype=np.float64)
        et = np.asarray(data["cluster_Et"], dtype=np.float64)
        cent = np.asarray(data["centrality"], dtype=np.float64)
        split = np.asarray(data["split"], dtype=object)
        base = (
            (split == split_name)
            & np.isfinite(score)
            & np.isfinite(et)
            & np.isfinite(cent)
            & (et >= args.pt_min)
            & (et < args.pt_max)
        )
        for icol, (clo, chi) in enumerate(zip(cent_edges[:-1], cent_edges[1:])):
            ax = axes[irow, icol]
            mask = base & (cent >= clo) & (cent < chi)
            sig = score[mask & (y == 1)]
            bkg = score[mask & (y == 0)]
            auc = auc_rank(y[mask], score[mask])
            density_hist(ax, sig, bkg, args.bins)
            ax.text(0.04, 0.88, f"AUC = {auc:.3f}", transform=ax.transAxes, ha="left", va="top", fontsize=11.5)
            ax.text(0.04, 0.78, f"S {sig.size:,}  B {bkg.size:,}", transform=ax.transAxes, ha="left", va="top", fontsize=8.6, color="0.35")
            if irow == 0:
                ax.set_title(f"{clo:g}-{chi:g}% centrality", fontsize=13, fontweight="bold")
            if icol == 0:
                ax.set_ylabel("Density", fontsize=11)
                ax.text(-0.19, 0.5, row_label, transform=ax.transAxes, rotation=90, ha="center", va="center", fontsize=12.5, fontweight="bold")
            if irow == len(datasets) - 1:
                ax.set_xlabel("BDT score", fontsize=11)
            if irow == 0 and icol == len(cent_edges) - 2:
                ax.legend(loc="upper right", fontsize=9, frameon=False)
            ax.grid(True, color="0.90", lw=0.5)
            ax.set_xlim(0, 1)
            summary.append(
                {
                    "split": split_name,
                    "centrality": f"{clo:g}-{chi:g}",
                    "pt_min": args.pt_min,
                    "pt_max": args.pt_max,
                    "entries": int(mask.sum()),
                    "signal_entries": int(sig.size),
                    "background_entries": int(bkg.size),
                    "auc": auc,
                    "signal_score_mean": float(np.mean(sig)) if sig.size else math.nan,
                    "background_score_mean": float(np.mean(bkg)) if bkg.size else math.nan,
                }
            )

    fig.text(0.035, 0.915, "sPHENIX", ha="left", va="top", fontsize=14, fontweight="bold", fontstyle="italic")
    fig.text(0.122, 0.915, "Internal", ha="left", va="top", fontsize=14)
    title = f"BDT training closure: signal/background score separation, {args.pt_min:g} < $E_T$ < {args.pt_max:g} GeV"
    fig.text(0.52, 0.965, title, ha="center", va="top", fontsize=16.5, fontweight="bold")
    subtitle = f"{metadata['model_id']}: training split vs internal test vs independent embedded validation"
    fig.text(0.52, 0.932, subtitle, ha="center", va="top", fontsize=10.5, color="0.35")
    fig.tight_layout(rect=[0.075, 0.045, 0.995, 0.89], h_pad=1.0, w_pad=1.6)

    tag = metadata["model_id"].replace("/", "_")
    pt_tag = f"{args.pt_min:g}to{args.pt_max:g}".replace(".", "p")
    png = args.outdir / f"bdt_training_closure_{tag}_pt{pt_tag}_cent3.png"
    fig.savefig(png, dpi=220)
    plt.close(fig)

    csv_path = args.outdir / f"bdt_training_closure_{tag}_summary.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary[0]))
        writer.writeheader()
        writer.writerows(summary)
    meta = {**metadata, "plot": str(png), "summary_csv": str(csv_path), "selection": f"{args.pt_min:g} <= cluster_Et < {args.pt_max:g}", "centrality_edges": cent_edges}
    meta_path = args.outdir / f"bdt_training_closure_{tag}_metadata.json"
    meta_path.write_text(json.dumps(meta, indent=2, sort_keys=True) + "\n")
    print(f"[bdtClosure] plot={png}")
    print(f"[bdtClosure] summary={csv_path}")
    print(f"[bdtClosure] metadata={meta_path}")
    for row in summary:
        print(
            f"[bdtClosure] {row['split']} cent={row['centrality']} "
            f"auc={row['auc']:.5f} S={row['signal_entries']} B={row['background_entries']}"
        )


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--registry", type=Path, required=True)
    ap.add_argument("--train-cache", type=Path, required=True)
    ap.add_argument("--validation-cache", type=Path, required=True, help="score_caches.list from the independent validation report")
    ap.add_argument("--model-id", default="")
    ap.add_argument("--product", default="")
    ap.add_argument("--validation-product", default="", help="Validation score product name; defaults to --product or --model-id")
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--pt-min", type=float, default=15.0)
    ap.add_argument("--pt-max", type=float, default=35.0)
    ap.add_argument("--cent-bins", default="0,20,50,80")
    ap.add_argument("--bins", type=int, default=34)
    ap.add_argument("--validation-max-rows", type=int, default=0)
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    model_key = args.model_id or args.product
    validation_product = args.validation_product or args.product or args.model_id
    if not model_key:
        raise SystemExit("Pass --model-id or --product")
    if not validation_product:
        raise SystemExit("Pass --validation-product, --product, or --model-id")
    training, metadata = load_training_scores(args.registry, args.train_cache, args.model_id or None, args.product or None)
    validation, validation_counts = load_validation_scores(args.validation_cache, validation_product, args.validation_max_rows)
    metadata["validation_cache"] = str(args.validation_cache)
    metadata["validation_product"] = validation_product
    metadata["validation_counts"] = validation_counts
    plot_and_write(training, validation, metadata, args)


if __name__ == "__main__":
    main()
