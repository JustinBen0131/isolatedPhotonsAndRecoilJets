#!/usr/bin/env python3
"""Make centrality-binned ROC overlays for the AuAu BDT+MLP stack.

This script is intentionally provenance-heavy because the slide-facing plot must
distinguish an exact exported stack artifact from earlier local diagnostic
stack fits.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import train_auau_stacked_bdt_mlp_sweep as stack  # noqa: E402


DEFAULT_PT_RANGE = (15.0, 35.0)
DEFAULT_CENT_BINS = "0,20,50,80"


def parse_edges(text: str) -> list[float]:
    vals = [float(tok) for tok in text.replace(";", ",").split(",") if tok.strip()]
    if len(vals) < 2 or any(vals[i] >= vals[i + 1] for i in range(len(vals) - 1)):
        raise SystemExit(f"Bad bin edges: {text}")
    return vals


def write_json(path: Path, payload) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(json_ready(payload), indent=2, sort_keys=True) + "\n")


def json_ready(value):
    if isinstance(value, dict):
        return {str(k): json_ready(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_ready(v) for v in value]
    if isinstance(value, np.ndarray):
        return json_ready(value.tolist())
    if isinstance(value, (np.floating, float)):
        return None if not math.isfinite(float(value)) else float(value)
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    return value


def sigmoid(x: np.ndarray | float) -> np.ndarray:
    x = np.asarray(x, dtype="float64")
    out = np.empty_like(x)
    pos = x >= 0.0
    out[pos] = 1.0 / (1.0 + np.exp(-x[pos]))
    expx = np.exp(x[~pos])
    out[~pos] = expx / (1.0 + expx)
    return out


def roc_curve_auc(y_true: np.ndarray, score: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    y = np.asarray(y_true, dtype="int32")
    s = np.asarray(score, dtype="float64")
    mask = np.isfinite(s) & np.isin(y, [0, 1])
    y = y[mask]
    s = s[mask]
    n_sig = int(np.sum(y == 1))
    n_bkg = int(np.sum(y == 0))
    if n_sig == 0 or n_bkg == 0:
        return np.array([0.0, 1.0]), np.array([0.0, 1.0]), math.nan
    order = np.argsort(-s, kind="mergesort")
    y_sorted = y[order]
    tp = np.cumsum(y_sorted == 1, dtype="float64")
    fp = np.cumsum(y_sorted == 0, dtype="float64")
    tpr = np.concatenate(([0.0], tp / n_sig, [1.0]))
    fpr = np.concatenate(([0.0], fp / n_bkg, [1.0]))
    auc = float(np.trapezoid(tpr, fpr))
    return fpr, tpr, auc


def tree_raw(tree: dict, x: np.ndarray) -> np.ndarray:
    children_left = np.asarray(tree["children_left"], dtype="int32")
    children_right = np.asarray(tree["children_right"], dtype="int32")
    feature = np.asarray(tree["feature"], dtype="int32")
    threshold = np.asarray(tree["threshold"], dtype="float64")
    value = np.asarray(tree["value"], dtype="float64")
    node = np.zeros(x.shape[0], dtype="int32")
    active = np.ones(x.shape[0], dtype=bool)
    row_index = np.arange(x.shape[0])
    guard = 0
    while active.any():
        if guard > 100000:
            raise RuntimeError("Tree traversal guard tripped")
        guard += 1
        active_rows = row_index[active]
        nodes = node[active]
        is_leaf = (children_left[nodes] < 0) & (children_right[nodes] < 0)
        if is_leaf.all():
            break
        branch_rows = active_rows[~is_leaf]
        branch_nodes = nodes[~is_leaf]
        fidx = feature[branch_nodes]
        xv = np.zeros(len(branch_rows), dtype="float64")
        good = (fidx >= 0) & (fidx < x.shape[1])
        xv[good] = x[branch_rows[good], fidx[good]]
        go_left = xv <= threshold[branch_nodes]
        node[branch_rows] = np.where(go_left, children_left[branch_nodes], children_right[branch_nodes])
        active[active_rows[is_leaf]] = False
    return value[node]


def score_submodel(model: dict, frame: dict[str, np.ndarray]) -> np.ndarray:
    names = list(model["feature_names"])
    x = np.column_stack([np.asarray(frame[name], dtype="float64") for name in names])
    impute = np.asarray(model.get("impute", np.zeros(len(names))), dtype="float64")
    x = np.where(np.isfinite(x), x, impute)
    kind = model.get("kind", "")
    if kind == "logistic":
        mean = np.asarray(model["mean"], dtype="float64")
        scale = np.asarray(model["scale"], dtype="float64")
        coef = np.asarray(model["coef"], dtype="float64")
        scale = np.where(np.abs(scale) > 1.0e-12, scale, 1.0)
        z = float(model["intercept"]) + ((x - mean) / scale) @ coef
        return sigmoid(z)
    if kind == "gradient_boosting_classifier":
        raw = np.full(x.shape[0], float(model["initial_raw_score"]), dtype="float64")
        learning_rate = float(model["learning_rate"])
        for tree in model["trees"]:
            raw += learning_rate * tree_raw(tree, x)
        return sigmoid(raw)
    if kind == "mlp":
        mean = np.asarray(model["mean"], dtype="float64")
        scale = np.asarray(model["scale"], dtype="float64")
        scale = np.where(np.abs(scale) > 1.0e-12, scale, 1.0)
        a = np.nan_to_num((x - mean) / scale, nan=0.0, posinf=0.0, neginf=0.0)
        for layer in model["layers"][:-1]:
            a = np.maximum(a @ np.asarray(layer["weight"], dtype="float64") + np.asarray(layer["bias"], dtype="float64"), 0.0)
        raw = a @ np.asarray(model["layers"][-1]["weight"], dtype="float64") + np.asarray(model["layers"][-1]["bias"], dtype="float64")
        return sigmoid(raw[:, 0])
    raise SystemExit(f"Unsupported stack submodel kind: {kind}")


def route_score(artifact: dict, frame: dict[str, np.ndarray]) -> np.ndarray:
    stack.add_interactions(frame)
    n = len(frame["is_signal"])
    score = np.full(n, np.nan, dtype="float64")
    if artifact.get("routes"):
        et = np.asarray(frame["cluster_Et"], dtype="float64")
        cent = np.asarray(frame["centrality"], dtype="float64")
        for route in artifact["routes"]:
            mask = np.isfinite(et) & (et >= float(route["pt_lo"])) & (et < float(route["pt_hi"]))
            cent_lo = route.get("cent_lo")
            cent_hi = route.get("cent_hi")
            if cent_lo is not None and cent_hi is not None:
                mask &= np.isfinite(cent) & (cent >= float(cent_lo)) & (cent < float(cent_hi))
            if not mask.any():
                continue
            subframe = {key: value[mask] for key, value in frame.items()}
            score[mask] = score_submodel(route["model"], subframe)
        return score
    if "model" not in artifact:
        raise SystemExit("Stack artifact has neither routes nor a global model")
    return score_submodel(artifact["model"], frame)


def make_plot(frame, stack_score, split_mask, args, cache_info, artifact):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    y = np.asarray(frame["is_signal"], dtype="int32")
    et = np.asarray(frame["cluster_Et"], dtype="float64")
    cent = np.asarray(frame["centrality"], dtype="float64")
    mlp_score = np.asarray(frame["mlp_score"], dtype="float64")
    bdt_score = np.asarray(frame["bdt_score"], dtype="float64")
    cent_edges = parse_edges(args.cent_bins)
    pt_lo, pt_hi = args.pt_min, args.pt_max
    base = split_mask & np.isfinite(et) & np.isfinite(cent) & (et >= pt_lo) & (et < pt_hi)

    colors = {
        "BDT": "#1f77b4",
        "MLP": "#d95f02",
        "BDT+MLP stack": "#6a3d9a",
    }
    scores = {
        "BDT": bdt_score,
        "MLP": mlp_score,
        "BDT+MLP stack": stack_score,
    }

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
    fig, axes = plt.subplots(1, len(cent_edges) - 1, figsize=(14.4, 4.45), sharex=True, sharey=True)
    if len(cent_edges) == 2:
        axes = [axes]
    rows: list[dict] = []
    for ax, clo, chi in zip(axes, cent_edges[:-1], cent_edges[1:], strict=True):
        cmask = base & (cent >= clo) & (cent < chi)
        n_sig = int(np.sum(cmask & (y == 1)))
        n_bkg = int(np.sum(cmask & (y == 0)))
        for label, values in scores.items():
            fpr, tpr, auc = roc_curve_auc(y[cmask], values[cmask])
            ax.plot(fpr, tpr, lw=2.25, color=colors[label], label=f"{label}  AUC={auc:.3f}")
            rows.append(
                {
                    "model": label,
                    "cent_lo": clo,
                    "cent_hi": chi,
                    "split": args.split,
                    "entries": int(np.sum(cmask)),
                    "signal_entries": n_sig,
                    "background_entries": n_bkg,
                    "finite_fraction": float(np.mean(np.isfinite(values[cmask]))) if np.sum(cmask) else math.nan,
                    "auc": auc,
                }
            )
        ax.plot([0, 1], [0, 1], color="0.70", lw=1.0, ls="--", zorder=0)
        ax.set_title(f"{clo:g}-{chi:g}% centrality", fontsize=12.5, fontweight="bold")
        ax.set_xlabel("False positive rate / background efficiency", fontsize=10.5)
        ax.grid(True, color="0.88", lw=0.6)
        ax.text(
            0.045,
            0.055,
            f"N={int(np.sum(cmask)):,}\nS={n_sig:,}, B={n_bkg:,}",
            transform=ax.transAxes,
            ha="left",
            va="bottom",
            fontsize=8.2,
            color="0.25",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.78, "pad": 1.8},
        )
        ax.legend(loc="lower right", fontsize=8.6, frameon=True, framealpha=0.92, edgecolor="0.75")
    axes[0].set_ylabel("True positive rate / signal efficiency", fontsize=10.5)
    for ax in axes:
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)

    fig.text(0.045, 0.965, r"$\bf{\it{sPHENIX}}$ Internal", ha="left", va="top", fontsize=14)
    fig.text(0.045, 0.918, "Au+Au embedded validation", ha="left", va="top", fontsize=10.2)
    fig.text(
        0.50,
        0.965,
        "ROC curves by centrality",
        ha="center",
        va="top",
        fontsize=16,
        fontweight="bold",
    )
    split_label = {"val": "held-out validation split", "test": "held-out test split", "all": "all scored rows"}[args.split]
    fig.text(
        0.50,
        0.918,
        f"{pt_lo:g} < E_T < {pt_hi:g} GeV; {split_label}; exact BDT+MLP stack artifact",
        ha="center",
        va="top",
        fontsize=10.8,
        color="0.35",
    )
    fig.text(
        0.50,
        0.020,
        "BDT+MLP stack uses the exported JSON artifact evaluated on aligned enriched validation caches.",
        ha="center",
        va="bottom",
        fontsize=8.3,
        color="0.35",
    )
    fig.tight_layout(rect=[0.035, 0.065, 0.995, 0.86], w_pad=1.8)

    pt_tag = f"{args.pt_min:g}to{args.pt_max:g}".replace(".", "p")
    png = outdir / f"roc_overlay_bdt_mlp_exact_stack_by_centrality_pt{pt_tag}_{args.split}.png"
    fig.savefig(png, dpi=220)
    plt.close(fig)

    csv_path = outdir / "roc_overlay_auc_summary.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    metadata = {
        "schema": "RECOILJETS_AUAU_BDT_MLP_STACK_ROC_OVERLAY_V1",
        "selection": f"{pt_lo:g} < cluster_Et < {pt_hi:g} GeV",
        "centrality_bins": cent_edges,
        "split": args.split,
        "split_description": split_label,
        "artifact": str(args.artifact),
        "artifact_name": artifact.get("name", ""),
        "mlp_cache": str(args.mlp_cache),
        "bdt_cache": str(args.bdt_cache),
        "mlp_score": cache_info["mlp_score_key"],
        "bdt_score": cache_info["bdt_score_key"],
        "rows_loaded": cache_info["rows"],
        "plot": str(png),
        "summary_csv": str(csv_path),
        "note": "The BDT+MLP curve is the exact exported stack JSON artifact, not the older local logistic side diagnostic.",
    }
    write_json(outdir / "roc_overlay_metadata.json", metadata)
    return png, csv_path, outdir / "roc_overlay_metadata.json"


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mlp-cache", type=Path, required=True)
    ap.add_argument("--bdt-cache", type=Path, required=True)
    ap.add_argument("--artifact", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--mlp-score", default=stack.DEFAULT_MLP_SCORE)
    ap.add_argument("--bdt-score", default=stack.DEFAULT_BDT_SCORE)
    ap.add_argument("--pt-min", type=float, default=DEFAULT_PT_RANGE[0])
    ap.add_argument("--pt-max", type=float, default=DEFAULT_PT_RANGE[1])
    ap.add_argument("--cent-bins", default=DEFAULT_CENT_BINS)
    ap.add_argument("--split", choices=["val", "test", "all"], default="test")
    ap.add_argument("--train-fraction", type=float, default=0.60)
    ap.add_argument("--val-fraction", type=float, default=0.20)
    ap.add_argument("--random-seed", type=int, default=24681357)
    ap.add_argument("--max-shards", type=int, default=0)
    ap.add_argument("--max-rows", type=int, default=0)
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    if args.pt_min >= args.pt_max:
        raise SystemExit("--pt-min must be smaller than --pt-max")
    artifact = json.loads(args.artifact.read_text())
    load_args = argparse.Namespace(
        mlp_cache=args.mlp_cache,
        bdt_cache=args.bdt_cache,
        mlp_score=args.mlp_score,
        bdt_score=args.bdt_score,
        max_shards=args.max_shards,
        max_rows=args.max_rows,
        random_seed=args.random_seed,
    )
    frame, cache_info = stack.load_aligned_caches(load_args)
    y = np.asarray(frame["is_signal"], dtype="int32")
    masks = stack.stratified_split(y, args.random_seed, args.train_fraction, args.val_fraction)
    split_mask = np.ones(len(y), dtype=bool) if args.split == "all" else masks[args.split]
    stack_score = route_score(artifact, frame)
    png, csv_path, metadata = make_plot(frame, stack_score, split_mask, args, cache_info, artifact)
    print(f"[stackROC] plot={png}", flush=True)
    print(f"[stackROC] summary={csv_path}", flush=True)
    print(f"[stackROC] metadata={metadata}", flush=True)


if __name__ == "__main__":
    main()
