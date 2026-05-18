#!/usr/bin/env python3
"""Score pp photon-ID training tables with BDT/MLP registries.

This is the pp table-level counterpart of the AuAu validation cache step.  It
keeps the output cache contract expected by the BDT+MLP stacker while avoiding
the AuAu runtime centrality selection assumptions.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path
from typing import Iterable

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from train_auau_photon_bdt import add_derived_features, expand_input_paths, expand_required_columns, load_frame  # noqa: E402
from train_auau_photon_mlp import (  # noqa: E402
    MODEL_SPECS,
    auc_score,
    load_mlp_artifact,
    predict_mlp_array,
    threshold_for_signal_efficiency,
)
from train_auau_stacked_bdt_mlp_sweep import FULL_FEATURES, ISOLATION_CONTEXT_FEATURES  # noqa: E402


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", nargs="+", type=Path, required=True, help="ROOT files or @manifest.list files.")
    ap.add_argument("--tree", default="AuAuPhotonIDTrainingTree")
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--kind", choices=["bdt", "mlp", "both"], default="both")
    ap.add_argument("--bdt-registry", type=Path, default=None)
    ap.add_argument("--mlp-registry", type=Path, default=None)
    ap.add_argument("--label-branch", default="is_signal")
    ap.add_argument("--weight-branch", default="event_weight")
    ap.add_argument("--pt-range", default="6:35")
    ap.add_argument("--centrality-range", default="-1:0")
    ap.add_argument("--pt-bins", default="6,10,15,20,25,35")
    ap.add_argument("--target-signal-efficiency", type=float, default=0.80)
    ap.add_argument("--max-rows", type=int, default=0, help="Optional deterministic row cap for smoke tests.")
    ap.add_argument(
        "--max-load-rows-per-class",
        type=int,
        default=0,
        help="Optional deterministic per-label row cap applied while reading ROOT inputs.",
    )
    ap.add_argument(
        "--max-load-rows",
        type=int,
        default=0,
        help="Optional deterministic total row cap applied while reading ROOT inputs.",
    )
    ap.add_argument("--random-seed", type=int, default=42)
    ap.add_argument(
        "--skip-missing-tree",
        action="store_true",
        help="Skip input ROOT files that do not contain the requested training tree instead of failing immediately.",
    )
    return ap.parse_args()


def parse_range(text: str) -> tuple[float, float]:
    lo_s, hi_s = str(text).split(":", 1)
    lo = float(lo_s)
    hi = float(hi_s)
    if hi <= lo:
        raise SystemExit(f"Bad range: {text}")
    return lo, hi


def parse_edges(text: str) -> list[float]:
    edges = [float(tok) for tok in str(text).replace(";", ",").split(",") if tok.strip()]
    if len(edges) < 2 or any(edges[i] >= edges[i + 1] for i in range(len(edges) - 1)):
        raise SystemExit(f"Bad bin edges: {text}")
    return edges


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


def load_registry(path: Path | None, label: str) -> dict:
    if path is None:
        raise SystemExit(f"--{label}-registry is required for this validation kind")
    if not path.is_file():
        raise SystemExit(f"Missing {label} registry: {path}")
    return json.loads(path.read_text())


def registry_models(registry: dict) -> list[dict]:
    models = registry.get("models", [])
    if not isinstance(models, list) or not models:
        raise SystemExit("Registry has no models")
    return models


def model_features(model: dict) -> list[str]:
    features = model.get("features")
    if not isinstance(features, list) or not features:
        raise SystemExit(f"Model {model.get('model_id') or model.get('product')} has no features")
    return [str(item) for item in features]


def required_columns_for(args, bdt_registry: dict | None, mlp_registry: dict | None) -> list[str]:
    columns: set[str] = {
        args.label_branch,
        "cluster_Et",
        "cluster_Eta",
        "centrality",
        "vertexz",
        "reco_eiso",
    }
    columns.update(expand_required_columns(FULL_FEATURES + ISOLATION_CONTEXT_FEATURES))
    if bdt_registry:
        for model in registry_models(bdt_registry):
            columns.update(expand_required_columns(model_features(model)))
    if mlp_registry:
        for model in registry_models(mlp_registry):
            columns.update(expand_required_columns(model_features(model)))
    return sorted(columns)


def select_rows(frame, args):
    label = np.asarray(frame[args.label_branch], dtype="int32")
    et = np.asarray(frame["cluster_Et"], dtype="float64")
    cent = np.asarray(frame["centrality"], dtype="float64")
    pt_lo, pt_hi = parse_range(args.pt_range)
    cent_lo, cent_hi = parse_range(args.centrality_range)
    mask = (
        np.isin(label, [0, 1])
        & np.isfinite(et)
        & (et >= pt_lo)
        & (et < pt_hi)
        & np.isfinite(cent)
        & (cent >= cent_lo)
        & (cent < cent_hi)
    )
    selected = frame.loc[mask].copy()
    if args.max_rows > 0 and len(selected) > args.max_rows:
        rng = np.random.default_rng(args.random_seed)
        keep = np.sort(rng.choice(np.arange(len(selected)), size=args.max_rows, replace=False))
        selected = selected.iloc[keep].copy()
    return selected


def feature_matrix(frame, features: Iterable[str]) -> np.ndarray:
    names = list(features)
    missing = [name for name in names if name not in frame.columns]
    if missing:
        raise SystemExit("Frame is missing model features: " + ", ".join(missing))
    return np.column_stack([frame[name].to_numpy(dtype="float32") for name in names])


def score_bdt_models(frame, registry: dict) -> dict[str, np.ndarray]:
    try:
        import xgboost as xgb
    except ImportError as exc:
        raise SystemExit("BDT validation needs xgboost in the active Python environment") from exc

    scores: dict[str, np.ndarray] = {}
    for model in registry_models(registry):
        product = str(model.get("product") or model.get("model_id"))
        xgb_path = model.get("output_xgb_json") or model.get("report", {}).get("output_xgb_json")
        if not xgb_path:
            raise SystemExit(f"BDT model {product} has no output_xgb_json in registry")
        path = Path(xgb_path)
        if not path.is_file():
            raise SystemExit(f"Missing BDT XGBoost artifact for {product}: {path}")
        booster = xgb.Booster()
        booster.load_model(str(path))
        x = feature_matrix(frame, model_features(model))
        scores[f"score_{product}"] = booster.predict(xgb.DMatrix(x)).astype("float32")
    return scores


def score_mlp_models(frame, registry: dict) -> dict[str, np.ndarray]:
    scores: dict[str, np.ndarray] = {}
    for model in registry_models(registry):
        product = str(model.get("product"))
        artifact_path = model.get("output_json")
        if not artifact_path and product in MODEL_SPECS:
            registry_path = Path(model.get("metadata", ".")).parent
            artifact_path = str(registry_path / MODEL_SPECS[product]["filename"])
        if not artifact_path:
            raise SystemExit(f"MLP model {product} has no output_json in registry")
        path = Path(artifact_path)
        if not path.is_file():
            raise SystemExit(f"Missing MLP artifact for {product}: {path}")
        artifact = load_mlp_artifact(path)
        x = feature_matrix(frame, artifact["features"])
        scores[f"score_{product}"] = predict_mlp_array(artifact, x).astype("float32")
    return scores


def metric_rows(frame, scores: dict[str, np.ndarray], args) -> tuple[list[dict], dict]:
    labels = frame[args.label_branch].to_numpy(dtype="int32")
    et = frame["cluster_Et"].to_numpy(dtype="float64")
    weights = frame[args.weight_branch].to_numpy(dtype="float64") if args.weight_branch in frame.columns else None
    edges = parse_edges(args.pt_bins)
    rows: list[dict] = []
    summary: dict[str, dict] = {}
    for score_name, score in sorted(scores.items()):
        wp = threshold_for_signal_efficiency(labels, score, args.target_signal_efficiency)
        item = {
            "model": score_name,
            "bin": "inclusive",
            "entries": int(len(labels)),
            "signal": int((labels == 1).sum()),
            "background": int((labels == 0).sum()),
            "auc": auc_score(labels, score, weights),
            "wp80_threshold": None if wp is None else wp["threshold"],
            "wp80_fake_rate": None if wp is None else wp["background_fake_rate"],
        }
        rows.append(item)
        summary[score_name] = dict(item)
        for lo, hi in zip(edges[:-1], edges[1:], strict=True):
            mask = np.isfinite(et) & (et >= lo) & (et < hi)
            if not mask.any():
                continue
            wp_bin = threshold_for_signal_efficiency(labels[mask], score[mask], args.target_signal_efficiency)
            rows.append(
                {
                    "model": score_name,
                    "bin": f"{lo:g}-{hi:g}",
                    "entries": int(mask.sum()),
                    "signal": int((labels[mask] == 1).sum()),
                    "background": int((labels[mask] == 0).sum()),
                    "auc": auc_score(labels[mask], score[mask], None if weights is None else weights[mask]),
                    "wp80_threshold": None if wp_bin is None else wp_bin["threshold"],
                    "wp80_fake_rate": None if wp_bin is None else wp_bin["background_fake_rate"],
                }
            )
    return rows, summary


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def write_score_cache(frame, scores: dict[str, np.ndarray], args, outdir: Path) -> Path:
    payload = {}
    cache_columns = [
        args.label_branch,
        "cluster_Et",
        "centrality",
        "reco_eiso",
        "cluster_Eta",
        "vertexz",
        "run",
        "evt",
    ]
    for name in FULL_FEATURES + ISOLATION_CONTEXT_FEATURES:
        if name not in cache_columns:
            cache_columns.append(name)
    for name in cache_columns:
        if name in frame.columns:
            payload[name] = frame[name].to_numpy()
    if args.label_branch != "is_signal":
        payload["is_signal"] = frame[args.label_branch].to_numpy(dtype="int32")
    payload.update(scores)
    payload["counts_json"] = np.array(
        json.dumps(
            {
                "rows": int(len(frame)),
                "signal": int((frame[args.label_branch].to_numpy(dtype="int32") == 1).sum()),
                "background": int((frame[args.label_branch].to_numpy(dtype="int32") == 0).sum()),
            },
            sort_keys=True,
        )
    )
    path = outdir / "score_caches" / "pp_table_score_cache_full.npz"
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, **payload)
    manifest = outdir / "score_caches.list"
    manifest.write_text(str(path) + "\n")
    return path


def write_plots(frame, scores: dict[str, np.ndarray], args, outdir: Path) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return

    labels = frame[args.label_branch].to_numpy(dtype="int32")
    fig, ax = plt.subplots(figsize=(7.0, 5.0))
    for score_name, score in sorted(scores.items()):
        try:
            from sklearn.metrics import roc_curve

            fpr, tpr, _ = roc_curve(labels, score)
        except Exception:
            continue
        ax.plot(fpr, tpr, lw=2.0, label=f"{score_name} AUC={auc_score(labels, score):.3f}")
    ax.plot([0, 1], [0, 1], color="0.6", lw=1.0, ls="--")
    ax.set_xlabel("Background efficiency")
    ax.set_ylabel("Signal efficiency")
    ax.set_title("pp photon-ID ROC on held table rows")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "roc_summary.png", dpi=180)
    plt.close(fig)

    for score_name, score in sorted(scores.items()):
        fig, ax = plt.subplots(figsize=(7.0, 5.0))
        ax.hist(score[labels == 0], bins=60, range=(0, 1), histtype="step", lw=1.8, density=True, label="background")
        ax.hist(score[labels == 1], bins=60, range=(0, 1), histtype="step", lw=1.8, density=True, label="signal")
        ax.set_xlabel("Classifier score")
        ax.set_ylabel("Unit-normalized entries")
        ax.set_title(score_name)
        ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(outdir / f"{score_name}_score_separation.png", dpi=180)
        plt.close(fig)


def main() -> int:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    bdt_registry = load_registry(args.bdt_registry, "bdt") if args.kind in ("bdt", "both") else None
    mlp_registry = load_registry(args.mlp_registry, "mlp") if args.kind in ("mlp", "both") else None
    paths = expand_input_paths(args.input)
    required = required_columns_for(args, bdt_registry, mlp_registry)
    frame, optional_seen = load_frame(
        paths,
        args.tree,
        required,
        [args.weight_branch, "run", "evt"],
        args.label_branch,
        None,
        skip_missing_tree=args.skip_missing_tree,
        label_branch=args.label_branch,
        max_load_rows_per_class=args.max_load_rows_per_class,
        max_load_rows=args.max_load_rows,
        load_sample_seed=args.random_seed,
    )
    frame = add_derived_features(frame)
    frame = select_rows(frame, args)
    if len(frame) == 0:
        raise SystemExit("No rows selected for validation")
    scores: dict[str, np.ndarray] = {}
    if bdt_registry is not None:
        scores.update(score_bdt_models(frame, bdt_registry))
    if mlp_registry is not None:
        scores.update(score_mlp_models(frame, mlp_registry))
    rows, summary = metric_rows(frame, scores, args)
    cache_path = write_score_cache(frame, scores, args, args.outdir)
    write_csv(args.outdir / "pp_table_validation_metrics.csv", rows)
    write_json(
        args.outdir / "pp_table_validation_summary.json",
        {
            "schema": "RJ_PP_PHOTON_ML_TABLE_VALIDATION_V1",
            "kind": args.kind,
            "tree": args.tree,
            "input_files": [str(path) for path in paths],
            "rows": int(len(frame)),
            "optional_branches_seen": optional_seen,
            "pt_range": args.pt_range,
            "centrality_range": args.centrality_range,
            "score_cache": str(cache_path),
            "models": summary,
        },
    )
    write_plots(frame, scores, args, args.outdir)
    print(f"[validatePPPhotonML] wrote {args.outdir} rows={len(frame)} scores={len(scores)} cache={cache_path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
