#!/usr/bin/env python3
"""Train a diagnostic BDT+MLP stacked photon-ID calibrator from score caches.

This is a validation-side ceiling test, not a production runtime mode.  It uses
already scored, row-compatible validation caches to test whether the best BDT
and an MLP contain complementary information.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import pickle
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from train_auau_photon_mlp import (  # noqa: E402
    auc_score,
    calibration_report,
    parse_bin_spec,
    threshold_for_signal_efficiency,
    write_json,
)


DEFAULT_MLP_SCORE = "score_centInputBase3x3WidthRatiosMLP_pt1535"
DEFAULT_BDT_SCORE = "score_ptFine_cent7"
DEFAULT_PT_BINS = "15,18,20,22.5,25,30,35"
REQUIRED_ALIGNMENT_KEYS = ("is_signal", "cluster_Et", "centrality", "reco_eiso")
EPS = 1.0e-6


@dataclass
class LinearModel:
    feature_names: list[str]
    mean: np.ndarray
    scale: np.ndarray
    coef: np.ndarray
    intercept: float


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mlp-cache", type=Path, required=True)
    ap.add_argument("--bdt-cache", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--mlp-score", default=DEFAULT_MLP_SCORE)
    ap.add_argument("--bdt-score", default=DEFAULT_BDT_SCORE)
    ap.add_argument("--pt-bins", default=DEFAULT_PT_BINS)
    ap.add_argument("--max-shards", type=int, default=0)
    ap.add_argument("--max-rows", type=int, default=0)
    ap.add_argument("--random-seed", type=int, default=24681357)
    ap.add_argument("--train-fraction", type=float, default=0.60)
    ap.add_argument("--val-fraction", type=float, default=0.20)
    ap.add_argument("--models", default="logistic,gbm")
    ap.add_argument("--gbm-max-iter", type=int, default=120)
    ap.add_argument("--gbm-max-leaf-nodes", type=int, default=8)
    ap.add_argument("--gbm-learning-rate", type=float, default=0.045)
    ap.add_argument("--l2", type=float, default=2.0e-3)
    ap.add_argument("--max-linear-steps", type=int, default=2400)
    ap.add_argument("--target-signal-efficiency", type=float, default=0.80)
    ap.add_argument("--bdt-inclusive-anchor", type=float, default=0.886)
    ap.add_argument("--primary-mlp-inclusive-anchor", type=float, default=0.871)
    ap.add_argument("--note", default="")
    return ap.parse_args()


def read_manifest(path: Path) -> list[Path]:
    if not path.is_file():
        raise SystemExit(f"Missing cache manifest: {path}")
    paths = [Path(line.strip()) for line in path.read_text().splitlines() if line.strip()]
    if not paths:
        raise SystemExit(f"Empty cache manifest: {path}")
    return paths


def load_score_cache(path: Path):
    if not path.is_file():
        raise SystemExit(f"Missing score cache: {path}")
    return np.load(path, allow_pickle=True)


def auto_score_key(cache, requested: str, label: str) -> str:
    if requested in cache.files:
        return requested
    score_keys = [k for k in cache.files if k.startswith("score_")]
    if len(score_keys) == 1:
        print(f"[stackBDTMLP] {label}: requested {requested} not found; using only score key {score_keys[0]}", flush=True)
        return score_keys[0]
    raise SystemExit(
        f"{label}: requested score key {requested} not found. Available score keys: "
        + ", ".join(score_keys)
    )


def assert_aligned(mlp_cache, bdt_cache, shard_label: str):
    for key in REQUIRED_ALIGNMENT_KEYS:
        if key not in mlp_cache.files or key not in bdt_cache.files:
            raise SystemExit(f"{shard_label}: missing alignment key {key}")
        a = mlp_cache[key]
        b = bdt_cache[key]
        if len(a) != len(b):
            raise SystemExit(f"{shard_label}: row count mismatch for {key}: {len(a)} vs {len(b)}")
        if np.issubdtype(a.dtype, np.number):
            if not np.allclose(a.astype("float64"), b.astype("float64"), rtol=0.0, atol=1.0e-6, equal_nan=True):
                diff = float(np.nanmax(np.abs(a.astype("float64") - b.astype("float64"))))
                raise SystemExit(f"{shard_label}: row alignment mismatch for {key}; max abs diff {diff}")
        elif not np.array_equal(a, b):
            raise SystemExit(f"{shard_label}: row alignment mismatch for {key}")


def load_aligned_caches(args: argparse.Namespace):
    mlp_paths = read_manifest(args.mlp_cache)
    bdt_paths = read_manifest(args.bdt_cache)
    if len(mlp_paths) != len(bdt_paths):
        raise SystemExit(f"Cache manifest length mismatch: MLP {len(mlp_paths)} vs BDT {len(bdt_paths)}")
    if args.max_shards > 0:
        mlp_paths = mlp_paths[: args.max_shards]
        bdt_paths = bdt_paths[: args.max_shards]

    pieces = {key: [] for key in REQUIRED_ALIGNMENT_KEYS}
    pieces["mlp_score"] = []
    pieces["bdt_score"] = []
    mlp_key = args.mlp_score
    bdt_key = args.bdt_score
    for idx, (mlp_path, bdt_path) in enumerate(zip(mlp_paths, bdt_paths, strict=True), 1):
        if idx == 1 or idx % 20 == 0 or idx == len(mlp_paths):
            print(f"[stackBDTMLP] reading shard {idx}/{len(mlp_paths)}", flush=True)
        mlp_cache = load_score_cache(mlp_path)
        bdt_cache = load_score_cache(bdt_path)
        if idx == 1:
            mlp_key = auto_score_key(mlp_cache, args.mlp_score, "MLP cache")
            bdt_key = auto_score_key(bdt_cache, args.bdt_score, "BDT cache")
            print(f"[stackBDTMLP] mlp_score_key={mlp_key}", flush=True)
            print(f"[stackBDTMLP] bdt_score_key={bdt_key}", flush=True)
        assert_aligned(mlp_cache, bdt_cache, f"shard {idx}")
        for key in REQUIRED_ALIGNMENT_KEYS:
            pieces[key].append(np.asarray(mlp_cache[key]))
        pieces["mlp_score"].append(np.asarray(mlp_cache[mlp_key], dtype="float64"))
        pieces["bdt_score"].append(np.asarray(bdt_cache[bdt_key], dtype="float64"))

    frame = {key: np.concatenate(parts) for key, parts in pieces.items()}
    if args.max_rows > 0 and len(frame["is_signal"]) > args.max_rows:
        rng = np.random.default_rng(args.random_seed)
        idx = rng.choice(len(frame["is_signal"]), size=args.max_rows, replace=False)
        idx.sort()
        frame = {key: value[idx] for key, value in frame.items()}
    frame["is_signal"] = frame["is_signal"].astype("int32")
    return frame, {"mlp_score_key": mlp_key, "bdt_score_key": bdt_key, "shards": len(mlp_paths)}


def sigmoid(x):
    x = np.asarray(x, dtype="float64")
    out = np.empty_like(x)
    pos = x >= 0
    out[pos] = 1.0 / (1.0 + np.exp(-x[pos]))
    expx = np.exp(x[~pos])
    out[~pos] = expx / (1.0 + expx)
    return out


def logit_prob(p):
    p = np.asarray(p, dtype="float64")
    p = np.clip(p, EPS, 1.0 - EPS)
    return np.log(p / (1.0 - p))


def impute_numeric(values):
    values = np.asarray(values, dtype="float64")
    finite = np.isfinite(values)
    fill = float(np.nanmedian(values[finite])) if finite.any() else 0.0
    clean = np.where(finite, values, fill)
    return clean, finite.astype("float64")


def make_feature_matrix(frame, kind: str):
    et = np.asarray(frame["cluster_Et"], dtype="float64")
    cent = np.asarray(frame["centrality"], dtype="float64")
    bdt, bdt_finite = impute_numeric(frame["bdt_score"])
    mlp, mlp_finite = impute_numeric(frame["mlp_score"])
    mlp_logit = logit_prob(mlp)
    log_et = np.log(np.clip(et, 1.0e-3, None))
    cent_scaled = np.clip(cent, 0.0, 100.0) / 100.0

    columns: list[tuple[str, np.ndarray]]
    if kind == "bdt_only":
        columns = [("bdt_score", bdt), ("bdt_is_finite", bdt_finite)]
    elif kind == "mlp_only":
        columns = [("mlp_score", mlp), ("mlp_logit", mlp_logit), ("mlp_is_finite", mlp_finite)]
    elif kind == "stack_scores":
        columns = [
            ("bdt_score", bdt),
            ("bdt_is_finite", bdt_finite),
            ("mlp_score", mlp),
            ("mlp_logit", mlp_logit),
            ("mlp_is_finite", mlp_finite),
        ]
    elif kind == "stack_context":
        columns = [
            ("bdt_score", bdt),
            ("bdt_is_finite", bdt_finite),
            ("mlp_score", mlp),
            ("mlp_logit", mlp_logit),
            ("mlp_is_finite", mlp_finite),
            ("log_cluster_Et", log_et),
            ("centrality_scaled", cent_scaled),
        ]
    elif kind == "stack_interactions":
        columns = [
            ("bdt_score", bdt),
            ("bdt_is_finite", bdt_finite),
            ("mlp_score", mlp),
            ("mlp_logit", mlp_logit),
            ("mlp_is_finite", mlp_finite),
            ("log_cluster_Et", log_et),
            ("centrality_scaled", cent_scaled),
            ("bdt_x_mlp_logit", bdt * mlp_logit),
            ("bdt_x_logEt", bdt * log_et),
            ("mlp_logit_x_logEt", mlp_logit * log_et),
            ("bdt_x_cent", bdt * cent_scaled),
            ("mlp_logit_x_cent", mlp_logit * cent_scaled),
        ]
    else:
        raise ValueError(f"Unknown feature kind: {kind}")
    names = [name for name, _ in columns]
    matrix = np.column_stack([value for _, value in columns]).astype("float64")
    return names, matrix


def stratified_split(y, seed: int, train_fraction: float, val_fraction: float):
    if not (0.0 < train_fraction < 1.0) or not (0.0 <= val_fraction < 1.0):
        raise SystemExit("Train/validation fractions are invalid")
    if train_fraction + val_fraction >= 1.0:
        raise SystemExit("Train fraction + validation fraction must leave a positive test split")
    rng = np.random.default_rng(seed)
    y = np.asarray(y, dtype="int32")
    split = np.full(len(y), "test", dtype=object)
    for cls in (0, 1):
        idx = np.flatnonzero(y == cls)
        rng.shuffle(idx)
        n_train = int(round(train_fraction * len(idx)))
        n_val = int(round(val_fraction * len(idx)))
        split[idx[:n_train]] = "train"
        split[idx[n_train : n_train + n_val]] = "val"
    return {
        "train": split == "train",
        "val": split == "val",
        "test": split == "test",
        "all": np.ones(len(y), dtype=bool),
    }


def class_balanced_weights(y):
    y = np.asarray(y, dtype="int32")
    w = np.ones(len(y), dtype="float64")
    for cls in (0, 1):
        mask = y == cls
        if mask.any():
            w[mask] = 0.5 * len(y) / float(mask.sum())
    return w


def standardize_train_apply(x_train, x):
    mean = np.nanmean(x_train, axis=0)
    scale = np.nanstd(x_train, axis=0)
    scale = np.where(np.isfinite(scale) & (scale > 1.0e-9), scale, 1.0)
    xz = (x - mean) / scale
    return np.nan_to_num(xz, nan=0.0, posinf=0.0, neginf=0.0), mean, scale


def fit_logistic_numpy(feature_names, x, y, train_mask, l2: float, max_steps: int, seed: int) -> LinearModel:
    x_train = x[train_mask]
    y_train = y[train_mask].astype("float64")
    x_train_z, mean, scale = standardize_train_apply(x_train, x_train)
    weights = class_balanced_weights(y_train.astype("int32"))
    rng = np.random.default_rng(seed)
    coef = rng.normal(0.0, 0.01, size=x_train_z.shape[1])
    intercept = 0.0
    m_coef = np.zeros_like(coef)
    v_coef = np.zeros_like(coef)
    m_int = 0.0
    v_int = 0.0
    lr = 0.035
    beta1 = 0.9
    beta2 = 0.999
    for step in range(1, max_steps + 1):
        prob = sigmoid(x_train_z @ coef + intercept)
        err = (prob - y_train) * weights / max(1.0, float(weights.mean()))
        grad_coef = (x_train_z.T @ err) / len(y_train) + l2 * coef
        grad_int = float(err.mean())
        m_coef = beta1 * m_coef + (1.0 - beta1) * grad_coef
        v_coef = beta2 * v_coef + (1.0 - beta2) * (grad_coef * grad_coef)
        m_int = beta1 * m_int + (1.0 - beta1) * grad_int
        v_int = beta2 * v_int + (1.0 - beta2) * (grad_int * grad_int)
        coef -= lr * (m_coef / (1.0 - beta1**step)) / (np.sqrt(v_coef / (1.0 - beta2**step)) + 1.0e-8)
        intercept -= lr * (m_int / (1.0 - beta1**step)) / (math.sqrt(v_int / (1.0 - beta2**step)) + 1.0e-8)
    return LinearModel(list(feature_names), mean, scale, coef, float(intercept))


def predict_linear(model: LinearModel, x):
    xz = np.nan_to_num((x - model.mean) / model.scale, nan=0.0, posinf=0.0, neginf=0.0)
    return sigmoid(xz @ model.coef + model.intercept)


def fit_logistic_sklearn_or_numpy(name, feature_names, x, y, train_mask, args):
    try:
        from sklearn.linear_model import LogisticRegression
        from sklearn.pipeline import make_pipeline
        from sklearn.preprocessing import StandardScaler

        clf = make_pipeline(
            StandardScaler(),
            LogisticRegression(
                C=max(1.0e-6, 1.0 / max(args.l2, 1.0e-9)),
                class_weight="balanced",
                max_iter=2000,
                solver="lbfgs",
                random_state=args.random_seed,
            ),
        )
        clf.fit(x[train_mask], y[train_mask])
        return {
            "name": name,
            "kind": "logistic_sklearn",
            "feature_names": list(feature_names),
            "model": clf,
            "predict": lambda xx: clf.predict_proba(xx)[:, 1],
        }
    except Exception as exc:
        print(f"[stackBDTMLP] sklearn logistic unavailable for {name}; using numpy optimizer: {exc}", flush=True)
        model = fit_logistic_numpy(feature_names, x, y, train_mask, args.l2, args.max_linear_steps, args.random_seed)
        return {
            "name": name,
            "kind": "logistic_numpy",
            "feature_names": list(feature_names),
            "model": model,
            "predict": lambda xx: predict_linear(model, xx),
        }


def fit_gbm(name, feature_names, x, y, train_mask, args):
    from sklearn.ensemble import HistGradientBoostingClassifier

    clf = HistGradientBoostingClassifier(
        loss="log_loss",
        learning_rate=args.gbm_learning_rate,
        max_iter=args.gbm_max_iter,
        max_leaf_nodes=args.gbm_max_leaf_nodes,
        l2_regularization=0.10,
        early_stopping=True,
        validation_fraction=0.20,
        random_state=args.random_seed,
    )
    clf.fit(x[train_mask], y[train_mask], sample_weight=class_balanced_weights(y[train_mask]))
    return {
        "name": name,
        "kind": "hist_gradient_boosting",
        "feature_names": list(feature_names),
        "model": clf,
        "predict": lambda xx: clf.predict_proba(xx)[:, 1],
    }


def corrcoef_or_nan(x, y):
    x = np.asarray(x, dtype="float64")
    y = np.asarray(y, dtype="float64")
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 3:
        return math.nan
    if np.std(x[mask]) <= 0.0 or np.std(y[mask]) <= 0.0:
        return math.nan
    return float(np.corrcoef(x[mask], y[mask])[0, 1])


def split_metrics(frame, score, mask, target_eff, pt_bins):
    y = frame["is_signal"][mask].astype("int32")
    score_m = np.asarray(score[mask], dtype="float64")
    finite = np.isfinite(score_m) & np.isin(y, [0, 1])
    wp = threshold_for_signal_efficiency(y[finite], score_m[finite], target_eff)
    metrics = {
        "entries": int(mask.sum()),
        "signal_entries": int(np.sum(y == 1)),
        "background_entries": int(np.sum(y == 0)),
        "finite_fraction": float(np.mean(np.isfinite(score_m))) if len(score_m) else math.nan,
        "auc": auc_score(y, score_m),
        "wp80_fake": None if wp is None else wp["background_fake_rate"],
        "wp80_threshold": None if wp is None else wp["threshold"],
        "ece": calibration_report(y[finite], score_m[finite]).get("ece", math.nan) if finite.any() else math.nan,
        "score_vs_eiso_corr": corrcoef_or_nan(score_m, frame["reco_eiso"][mask]),
        "pt_bins": [],
    }
    et = np.asarray(frame["cluster_Et"][mask], dtype="float64")
    for lo, hi in pt_bins:
        bmask = np.isfinite(et) & (et >= lo) & (et < hi)
        if not bmask.any():
            continue
        wp_bin = threshold_for_signal_efficiency(y[bmask], score_m[bmask], target_eff)
        metrics["pt_bins"].append(
            {
                "lo": float(lo),
                "hi": float(hi),
                "entries": int(bmask.sum()),
                "auc": auc_score(y[bmask], score_m[bmask]),
                "wp80_fake": None if wp_bin is None else wp_bin["background_fake_rate"],
            }
        )
    return metrics


def flatten_rank_row(model_name, model_kind, split_name, feature_names, metrics, anchors):
    row = {
        "model": model_name,
        "kind": model_kind,
        "split": split_name,
        "features": "+".join(feature_names),
        "entries": metrics["entries"],
        "auc": metrics["auc"],
        "wp80_fake": metrics["wp80_fake"],
        "ece": metrics["ece"],
        "finite_fraction": metrics["finite_fraction"],
        "score_vs_eiso_corr": metrics["score_vs_eiso_corr"],
        "delta_auc_vs_broad_bdt": metrics["auc"] - anchors["bdt_inclusive"] if math.isfinite(metrics["auc"]) else math.nan,
        "delta_auc_vs_primary_mlp": metrics["auc"] - anchors["primary_mlp_inclusive"] if math.isfinite(metrics["auc"]) else math.nan,
    }
    for item in metrics["pt_bins"]:
        label = f"pt_{item['lo']:g}_{item['hi']:g}".replace(".", "p")
        row[f"{label}_auc"] = item["auc"]
        row[f"{label}_wp80_fake"] = item["wp80_fake"]
    return row


def write_csv(path: Path, rows: list[dict]):
    if not rows:
        return
    fields = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def json_safe_model(model_info):
    model = model_info["model"]
    if isinstance(model, LinearModel):
        return {
            "kind": model_info["kind"],
            "feature_names": model.feature_names,
            "mean": model.mean.tolist(),
            "scale": model.scale.tolist(),
            "coef": model.coef.tolist(),
            "intercept": model.intercept,
        }
    return {
        "kind": model_info["kind"],
        "feature_names": model_info["feature_names"],
        "note": "Model object stored in stacked_calibrator_models.pkl for diagnostic reproduction.",
    }


def main():
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    pt_bins = parse_bin_spec(args.pt_bins)
    if not pt_bins:
        raise SystemExit(f"Invalid --pt-bins: {args.pt_bins}")
    requested_models = {x.strip() for x in args.models.split(",") if x.strip()}

    print("[stackBDTMLP] RECOILJETS_AUAU_STACKED_BDT_MLP_CALIBRATOR_V1", flush=True)
    print(f"[stackBDTMLP] outdir={args.outdir}", flush=True)
    frame, cache_info = load_aligned_caches(args)
    y = frame["is_signal"].astype("int32")
    masks = stratified_split(y, args.random_seed, args.train_fraction, args.val_fraction)
    print(
        f"[stackBDTMLP] rows={len(y)} train={masks['train'].sum()} "
        f"val={masks['val'].sum()} test={masks['test'].sum()}",
        flush=True,
    )

    feature_kinds = ["bdt_only", "mlp_only", "stack_scores", "stack_context", "stack_interactions"]
    models = []
    for kind in feature_kinds:
        feature_names, x = make_feature_matrix(frame, kind)
        if "logistic" in requested_models:
            models.append((fit_logistic_sklearn_or_numpy(f"{kind}_logistic", feature_names, x, y, masks["train"], args), x))
        if "gbm" in requested_models and kind.startswith("stack"):
            try:
                models.append((fit_gbm(f"{kind}_tiny_gbm", feature_names, x, y, masks["train"], args), x))
            except Exception as exc:
                print(f"[stackBDTMLP] skipping GBM for {kind}: {exc}", flush=True)

    anchors = {"bdt_inclusive": args.bdt_inclusive_anchor, "primary_mlp_inclusive": args.primary_mlp_inclusive_anchor}
    all_metrics = {
        "schema": "RJ_AUAU_STACKED_BDT_MLP_CALIBRATOR_V1",
        "note": args.note,
        "cache_info": cache_info,
        "anchors": anchors,
        "splits": {name: int(mask.sum()) for name, mask in masks.items()},
        "models": {},
    }
    rows = []
    fitted_for_pickle = {}
    for model_info, x in models:
        model_name = model_info["name"]
        score = np.asarray(model_info["predict"](x), dtype="float64")
        all_metrics["models"][model_name] = {
            "kind": model_info["kind"],
            "features": model_info["feature_names"],
            "artifact": json_safe_model(model_info),
            "metrics": {},
        }
        fitted_for_pickle[model_name] = {
            "kind": model_info["kind"],
            "feature_names": model_info["feature_names"],
            "model": model_info["model"],
        }
        for split_name in ("val", "test", "all"):
            metrics = split_metrics(frame, score, masks[split_name], args.target_signal_efficiency, pt_bins)
            all_metrics["models"][model_name]["metrics"][split_name] = metrics
            rows.append(flatten_rank_row(model_name, model_info["kind"], split_name, model_info["feature_names"], metrics, anchors))
        test_auc = all_metrics["models"][model_name]["metrics"]["test"]["auc"]
        test_fake = all_metrics["models"][model_name]["metrics"]["test"]["wp80_fake"]
        print(f"[stackBDTMLP] model={model_name} kind={model_info['kind']} test_auc={test_auc:.5f} wp80_fake={test_fake}", flush=True)

    rank_path = args.outdir / "stacked_calibrator_rank_table.csv"
    write_csv(rank_path, rows)
    write_json(args.outdir / "stacked_calibrator_metrics.json", all_metrics)
    with (args.outdir / "stacked_calibrator_models.pkl").open("wb") as handle:
        pickle.dump(fitted_for_pickle, handle)

    test_rows = [row for row in rows if row["split"] == "test"]
    best = max(test_rows, key=lambda row: (-999.0 if not math.isfinite(float(row["auc"])) else float(row["auc"])))
    summary = {
        "best_test_model": best["model"],
        "best_test_auc": best["auc"],
        "best_test_wp80_fake": best["wp80_fake"],
        "rank_table": str(rank_path),
        "metrics_json": str(args.outdir / "stacked_calibrator_metrics.json"),
        "model_pickle": str(args.outdir / "stacked_calibrator_models.pkl"),
    }
    write_json(args.outdir / "stacked_calibrator_summary.json", summary)
    print(f"[stackBDTMLP] DONE_STACKED_BDT_MLP_OUTDIR={args.outdir}", flush=True)
    print(f"[stackBDTMLP] DONE_STACKED_BDT_MLP_RANK_TABLE={rank_path}", flush=True)
    print(f"[stackBDTMLP] BEST_TEST_MODEL={best['model']} BEST_TEST_AUC={best['auc']} BEST_TEST_WP80_FAKE={best['wp80_fake']}", flush=True)


if __name__ == "__main__":
    main()
