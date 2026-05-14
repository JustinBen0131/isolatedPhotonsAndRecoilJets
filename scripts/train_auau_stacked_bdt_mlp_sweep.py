#!/usr/bin/env python3
"""Train a full-feature diagnostic AuAu BDT+MLP stacker sweep.

The stackers use the runtime BDT score and MLP score as inputs.  They are
therefore diagnostic MC ceiling tests unless promoted through an explicit
physics review.  The serious variants keep the full photon-ID feature family:
BDT base inputs, standard/3x3 width moments, width ratios, and extended
runtime-safe shower/energy-ratio features when present in the score caches.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import pickle
import shutil
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from train_auau_photon_mlp import (  # noqa: E402
    auc_score,
    calibration_report,
    parse_bin_spec,
    threshold_for_signal_efficiency,
)


DEFAULT_MLP_SCORE = "score_centInputBase3x3WidthRatiosMLP_pt1535"
DEFAULT_BDT_SCORE = "score_ptFine_cent7"
DEFAULT_REPORT_PT_BINS = "5,10,15,18,20,22.5,25,30,35"
DEFAULT_REPORT_CENT_BINS = "0,10,20,30,40,50,60,80"
ALIGNMENT_KEYS = ("is_signal", "cluster_Et", "centrality", "reco_eiso")
EPS = 1.0e-6

BDT_BASE_FEATURES = [
    "cluster_Et",
    "cluster_weta_cogx",
    "cluster_wphi_cogx",
    "vertexz",
    "cluster_Eta",
    "e11_over_e33",
    "cluster_et1",
    "cluster_et2",
    "cluster_et3",
    "cluster_et4",
    "e32_over_e35",
]
WIDTH_3X3_FEATURES = ["cluster_weta33_cogx", "cluster_wphi33_cogx"]
WIDTH_RATIO_FEATURES = ["cluster_weta_over_wphi", "cluster_weta33_over_wphi33"]
EXTENDED_SHOWER_FEATURES = [
    "cluster_weta35_cogx",
    "cluster_wphi53_cogx",
    "cluster_w32",
    "cluster_w52",
    "cluster_w72",
    "e11_over_e22",
    "e11_over_e13",
    "e11_over_e15",
    "e11_over_e17",
    "e11_over_e31",
    "e11_over_e51",
    "e11_over_e71",
    "e22_over_e33",
    "e22_over_e35",
    "e22_over_e37",
    "e22_over_e53",
]
FULL_FEATURES = []
for _feature in BDT_BASE_FEATURES + WIDTH_3X3_FEATURES + WIDTH_RATIO_FEATURES + EXTENDED_SHOWER_FEATURES:
    if _feature not in FULL_FEATURES:
        FULL_FEATURES.append(_feature)


@dataclass(frozen=True)
class RouteSpec:
    label: str
    pt_lo: float
    pt_hi: float
    cent_lo: float | None = None
    cent_hi: float | None = None


@dataclass(frozen=True)
class VariantSpec:
    name: str
    pt_lo: float
    pt_hi: float
    context: str
    include_centrality: bool
    routes: tuple[RouteSpec, ...] = ()


@dataclass
class LinearArtifact:
    feature_names: list[str]
    impute: np.ndarray
    mean: np.ndarray
    scale: np.ndarray
    coef: np.ndarray
    intercept: float


@dataclass
class FittedModel:
    name: str
    algorithm: str
    feature_names: list[str]
    artifact: dict
    model_object: object
    history: list[dict] = field(default_factory=list)

    def predict(self, x: np.ndarray) -> np.ndarray:
        if self.algorithm == "logistic":
            model = self.model_object
            x_clean = np.where(np.isfinite(x), x, model.impute)
            xz = np.nan_to_num((x_clean - model.mean) / model.scale, nan=0.0, posinf=0.0, neginf=0.0)
            return sigmoid(xz @ model.coef + model.intercept)
        if self.algorithm == "nn":
            return nn_predict_from_artifact(self.artifact, x)
        x_clean = np.where(np.isfinite(x), x, np.asarray(self.artifact["impute"], dtype="float64"))
        return self.model_object.predict_proba(x_clean)[:, 1]


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mlp-cache", type=Path, default=None)
    ap.add_argument("--bdt-cache", type=Path, default=None)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--mlp-score", default=DEFAULT_MLP_SCORE)
    ap.add_argument("--bdt-score", default=DEFAULT_BDT_SCORE)
    ap.add_argument("--algorithms", default="logistic,gbm")
    ap.add_argument("--variants", default="", help="Comma-separated variant-name filter. Empty means use the selected sweep preset.")
    ap.add_argument("--sweep", choices=["full", "compact"], default="full")
    ap.add_argument("--include-controls", action="store_true")
    ap.add_argument("--allow-missing-full-features", action="store_true")
    ap.add_argument("--preflight-only", action="store_true")
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--max-shards", type=int, default=0)
    ap.add_argument("--max-rows", type=int, default=0)
    ap.add_argument("--train-fraction", type=float, default=0.60)
    ap.add_argument("--val-fraction", type=float, default=0.20)
    ap.add_argument("--random-seed", type=int, default=24681357)
    ap.add_argument("--target-signal-efficiency", type=float, default=0.80)
    ap.add_argument("--report-pt-bins", default=DEFAULT_REPORT_PT_BINS)
    ap.add_argument("--report-cent-bins", default=DEFAULT_REPORT_CENT_BINS)
    ap.add_argument("--top-n", type=int, default=4)
    ap.add_argument("--l2", type=float, default=2.0e-3)
    ap.add_argument("--linear-backend", choices=["numpy", "sklearn"], default="numpy")
    ap.add_argument("--max-linear-steps", type=int, default=1800)
    ap.add_argument("--gbm-estimators", type=int, default=90)
    ap.add_argument("--gbm-learning-rate", type=float, default=0.045)
    ap.add_argument("--gbm-max-depth", type=int, default=3)
    ap.add_argument("--gbm-max-leaf-nodes", type=int, default=8)
    ap.add_argument("--nn-hidden", default="64,32")
    ap.add_argument("--nn-epochs", type=int, default=180)
    ap.add_argument("--nn-patience", type=int, default=24)
    ap.add_argument("--nn-batch-size", type=int, default=8192)
    ap.add_argument("--nn-learning-rate", type=float, default=1.5e-3)
    ap.add_argument("--nn-l2", type=float, default=1.0e-3)
    ap.add_argument("--bdt-inclusive-anchor", type=float, default=0.886)
    ap.add_argument("--primary-mlp-inclusive-anchor", type=float, default=0.871)
    ap.add_argument("--bdt-20-25-anchor", type=float, default=0.778)
    ap.add_argument("--bdt-25-35-anchor", type=float, default=0.834)
    ap.add_argument("--note", default="")
    return ap.parse_args()


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


def safe_label(value: float) -> str:
    return f"{value:g}".replace(".", "p")


def range_label(lo: float, hi: float) -> str:
    return f"{safe_label(lo)}to{safe_label(hi)}"


def read_manifest(path: Path) -> list[Path]:
    if path is None or not path.is_file():
        raise SystemExit(f"Missing score-cache manifest: {path}")
    paths = [Path(x.strip()) for x in path.read_text().splitlines() if x.strip()]
    if not paths:
        raise SystemExit(f"Empty score-cache manifest: {path}")
    return paths


def auto_score_key(cache, requested: str, label: str) -> str:
    if requested in cache.files:
        return requested
    score_keys = [key for key in cache.files if key.startswith("score_")]
    if len(score_keys) == 1:
        print(f"[stackSweep] {label}: requested {requested} not found; using {score_keys[0]}", flush=True)
        return score_keys[0]
    raise SystemExit(f"{label}: missing requested score key {requested}; available scores: {', '.join(score_keys)}")


def compare_alignment(mlp_cache, bdt_cache, shard_label: str) -> None:
    for key in ALIGNMENT_KEYS:
        if key not in mlp_cache.files or key not in bdt_cache.files:
            raise SystemExit(f"{shard_label}: missing alignment key {key}")
        a = np.asarray(mlp_cache[key])
        b = np.asarray(bdt_cache[key])
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

    pieces: dict[str, list[np.ndarray]] = {"mlp_score": [], "bdt_score": []}
    mlp_key = args.mlp_score
    bdt_key = args.bdt_score
    copied_columns: list[str] = []
    for idx, (mlp_path, bdt_path) in enumerate(zip(mlp_paths, bdt_paths, strict=True), 1):
        if idx == 1 or idx % 20 == 0 or idx == len(mlp_paths):
            print(f"[stackSweep] reading shard {idx}/{len(mlp_paths)}", flush=True)
        with np.load(mlp_path, allow_pickle=True) as mlp_cache, np.load(bdt_path, allow_pickle=True) as bdt_cache:
            if idx == 1:
                mlp_key = auto_score_key(mlp_cache, args.mlp_score, "MLP cache")
                bdt_key = auto_score_key(bdt_cache, args.bdt_score, "BDT cache")
                copied_columns = sorted(
                    key
                    for key in set(mlp_cache.files) | set(bdt_cache.files)
                    if not key.startswith("score_") and key != "counts_json"
                )
                for key in copied_columns:
                    pieces.setdefault(key, [])
                print(f"[stackSweep] mlp_score_key={mlp_key}", flush=True)
                print(f"[stackSweep] bdt_score_key={bdt_key}", flush=True)
            compare_alignment(mlp_cache, bdt_cache, f"shard {idx}")
            for key in copied_columns:
                if key in mlp_cache.files:
                    pieces[key].append(np.asarray(mlp_cache[key]))
                elif key in bdt_cache.files:
                    pieces[key].append(np.asarray(bdt_cache[key]))
                elif key in ALIGNMENT_KEYS:
                    raise SystemExit(f"shard {idx}: missing required column {key}")
            pieces["mlp_score"].append(np.asarray(mlp_cache[mlp_key], dtype="float64"))
            pieces["bdt_score"].append(np.asarray(bdt_cache[bdt_key], dtype="float64"))

    frame = {key: np.concatenate(parts) for key, parts in pieces.items() if parts}
    if args.max_rows > 0 and len(frame["is_signal"]) > args.max_rows:
        rng = np.random.default_rng(args.random_seed)
        chosen = rng.choice(len(frame["is_signal"]), size=args.max_rows, replace=False)
        chosen.sort()
        frame = {key: value[chosen] for key, value in frame.items()}
    frame["is_signal"] = frame["is_signal"].astype("int32")
    add_derived_features(frame)
    return frame, {
        "mlp_cache": str(args.mlp_cache),
        "bdt_cache": str(args.bdt_cache),
        "mlp_score_key": mlp_key,
        "bdt_score_key": bdt_key,
        "shards": len(mlp_paths),
        "rows": int(len(frame["is_signal"])),
        "columns": sorted(frame.keys()),
    }


def add_derived_features(frame) -> None:
    def ratio(num: str, den: str):
        a = np.asarray(frame[num], dtype="float64")
        b = np.asarray(frame[den], dtype="float64")
        out = np.full(len(a), np.nan, dtype="float64")
        mask = np.isfinite(a) & np.isfinite(b) & (np.abs(b) > 1.0e-12)
        out[mask] = a[mask] / b[mask]
        return out

    if "cluster_weta_over_wphi" not in frame and {"cluster_weta_cogx", "cluster_wphi_cogx"}.issubset(frame):
        frame["cluster_weta_over_wphi"] = ratio("cluster_weta_cogx", "cluster_wphi_cogx")
    if "cluster_weta33_over_wphi33" not in frame and {"cluster_weta33_cogx", "cluster_wphi33_cogx"}.issubset(frame):
        frame["cluster_weta33_over_wphi33"] = ratio("cluster_weta33_cogx", "cluster_wphi33_cogx")
    if "mlp_logit" not in frame and "mlp_score" in frame:
        frame["mlp_logit"] = logit_prob(frame["mlp_score"])
    if "bdt_is_finite" not in frame and "bdt_score" in frame:
        frame["bdt_is_finite"] = np.isfinite(frame["bdt_score"]).astype("float64")
    if "mlp_is_finite" not in frame and "mlp_score" in frame:
        frame["mlp_is_finite"] = np.isfinite(frame["mlp_score"]).astype("float64")
    if "log_cluster_Et" not in frame and "cluster_Et" in frame:
        frame["log_cluster_Et"] = np.log(np.clip(np.asarray(frame["cluster_Et"], dtype="float64"), 1.0e-3, None))
    if "centrality_scaled" not in frame and "centrality" in frame:
        frame["centrality_scaled"] = np.clip(np.asarray(frame["centrality"], dtype="float64"), 0.0, 100.0) / 100.0


def add_interactions(frame) -> None:
    add_derived_features(frame)
    pairs = {
        "bdt_x_mlp_logit": ("bdt_score", "mlp_logit"),
        "bdt_x_logEt": ("bdt_score", "log_cluster_Et"),
        "mlp_logit_x_logEt": ("mlp_logit", "log_cluster_Et"),
        "bdt_x_cent": ("bdt_score", "centrality_scaled"),
        "mlp_logit_x_cent": ("mlp_logit", "centrality_scaled"),
        "logEt_x_cent": ("log_cluster_Et", "centrality_scaled"),
    }
    for out, (a, b) in pairs.items():
        if out not in frame and a in frame and b in frame:
            frame[out] = np.asarray(frame[a], dtype="float64") * np.asarray(frame[b], dtype="float64")


def sigmoid(x):
    x = np.asarray(x, dtype="float64")
    out = np.empty_like(x)
    pos = x >= 0
    out[pos] = 1.0 / (1.0 + np.exp(-x[pos]))
    expx = np.exp(x[~pos])
    out[~pos] = expx / (1.0 + expx)
    return out


def parse_hidden_layers(text: str) -> list[int]:
    layers = [int(tok) for tok in text.replace("x", ",").split(",") if tok.strip()]
    if not layers or any(layer <= 0 for layer in layers):
        raise SystemExit(f"Invalid --nn-hidden specification: {text}")
    return layers


def logit_prob(p):
    p = np.asarray(p, dtype="float64")
    p = np.clip(p, EPS, 1.0 - EPS)
    return np.log(p / (1.0 - p))


def class_balanced_weights(y):
    y = np.asarray(y, dtype="int32")
    weights = np.ones(len(y), dtype="float64")
    for cls in (0, 1):
        mask = y == cls
        if mask.any():
            weights[mask] = 0.5 * len(y) / float(mask.sum())
    return weights


def train_apply_impute(x_train, x):
    impute = np.nanmedian(np.where(np.isfinite(x_train), x_train, np.nan), axis=0)
    impute = np.where(np.isfinite(impute), impute, 0.0)
    return np.where(np.isfinite(x), x, impute), impute


def standardize_train_apply(x_train, x):
    mean = np.nanmean(x_train, axis=0)
    scale = np.nanstd(x_train, axis=0)
    mean = np.where(np.isfinite(mean), mean, 0.0)
    scale = np.where(np.isfinite(scale) & (scale > 1.0e-9), scale, 1.0)
    return (x - mean) / scale, mean, scale


def fit_logistic_numpy(name, feature_names, x, y, train_mask, args, seed, val_mask=None) -> FittedModel:
    x_clean, impute = train_apply_impute(x[train_mask], x)
    x_train_z, mean, scale = standardize_train_apply(x_clean[train_mask], x_clean[train_mask])
    x_all_z = np.nan_to_num((x_clean - mean) / scale, nan=0.0, posinf=0.0, neginf=0.0)
    y_train = y[train_mask].astype("float64")
    weights = class_balanced_weights(y[train_mask])
    val_mask = np.zeros(len(y), dtype=bool) if val_mask is None else val_mask
    y_val = y[val_mask].astype("float64")
    val_weights = class_balanced_weights(y[val_mask]) if val_mask.any() else np.ones(0, dtype="float64")
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
    history = []
    record_every = max(1, int(args.max_linear_steps) // 120)
    for step in range(1, args.max_linear_steps + 1):
        prob = sigmoid(x_train_z @ coef + intercept)
        err = (prob - y_train) * weights / max(1.0, float(weights.mean()))
        grad_coef = (x_train_z.T @ err) / len(y_train) + args.l2 * coef
        grad_int = float(err.mean())
        m_coef = beta1 * m_coef + (1.0 - beta1) * grad_coef
        v_coef = beta2 * v_coef + (1.0 - beta2) * (grad_coef * grad_coef)
        m_int = beta1 * m_int + (1.0 - beta1) * grad_int
        v_int = beta2 * v_int + (1.0 - beta2) * (grad_int * grad_int)
        coef -= lr * (m_coef / (1.0 - beta1**step)) / (np.sqrt(v_coef / (1.0 - beta2**step)) + 1.0e-8)
        intercept -= lr * (m_int / (1.0 - beta1**step)) / (math.sqrt(v_int / (1.0 - beta2**step)) + 1.0e-8)
        if step == 1 or step == args.max_linear_steps or step % record_every == 0:
            train_prob = sigmoid(x_all_z[train_mask] @ coef + intercept)
            l2_value = 0.5 * float(args.l2) * float(np.sum(coef * coef))
            row = {
                "step": step,
                "train_loss": binary_cross_entropy(y_train, train_prob, weights, l2_value),
                "train_auc": auc_score(y[train_mask], train_prob),
            }
            if val_mask.any():
                val_prob = sigmoid(x_all_z[val_mask] @ coef + intercept)
                row["val_loss"] = binary_cross_entropy(y_val, val_prob, val_weights, l2_value)
                row["val_auc"] = auc_score(y[val_mask], val_prob)
            else:
                row["val_loss"] = math.nan
                row["val_auc"] = math.nan
            history.append(row)
    artifact = LinearArtifact(list(feature_names), impute, mean, scale, coef, float(intercept))
    artifact_json = linear_artifact_json(artifact)
    artifact_json["training"] = {
        "backend": "numpy_adam",
        "steps_run": int(args.max_linear_steps),
        "learning_rate": lr,
        "l2": float(args.l2),
        "history": history,
    }
    return FittedModel(name, "logistic", list(feature_names), artifact_json, artifact, history)


def fit_logistic(name, feature_names, x, y, train_mask, args, seed, val_mask=None) -> FittedModel:
    if args.linear_backend == "numpy":
        return fit_logistic_numpy(name, feature_names, x, y, train_mask, args, seed, val_mask=val_mask)
    try:
        from sklearn.linear_model import LogisticRegression
        from sklearn.preprocessing import StandardScaler

        x_clean, impute = train_apply_impute(x[train_mask], x)
        scaler = StandardScaler()
        xz = scaler.fit_transform(x_clean[train_mask])
        clf = LogisticRegression(
            C=max(1.0e-6, 1.0 / max(args.l2, 1.0e-9)),
            class_weight="balanced",
            max_iter=2000,
            solver="lbfgs",
            random_state=seed,
        )
        clf.fit(xz, y[train_mask])
        artifact = LinearArtifact(
            list(feature_names),
            impute,
            np.asarray(scaler.mean_, dtype="float64"),
            np.asarray(scaler.scale_, dtype="float64"),
            np.asarray(clf.coef_[0], dtype="float64"),
            float(clf.intercept_[0]),
        )
        artifact_json = linear_artifact_json(artifact)
        x_all_z = np.nan_to_num((x_clean - artifact.mean) / artifact.scale, nan=0.0, posinf=0.0, neginf=0.0)
        train_prob = sigmoid(x_all_z[train_mask] @ artifact.coef + artifact.intercept)
        train_weights = class_balanced_weights(y[train_mask])
        history = [
            {
                "step": int(getattr(clf, "n_iter_", [0])[0]),
                "train_loss": binary_cross_entropy(y[train_mask], train_prob, train_weights, 0.5 * float(args.l2) * float(np.sum(artifact.coef * artifact.coef))),
                "train_auc": auc_score(y[train_mask], train_prob),
            }
        ]
        val_mask = np.zeros(len(y), dtype=bool) if val_mask is None else val_mask
        if val_mask.any():
            val_prob = sigmoid(x_all_z[val_mask] @ artifact.coef + artifact.intercept)
            val_weights = class_balanced_weights(y[val_mask])
            history[0]["val_loss"] = binary_cross_entropy(y[val_mask], val_prob, val_weights)
            history[0]["val_auc"] = auc_score(y[val_mask], val_prob)
        else:
            history[0]["val_loss"] = math.nan
            history[0]["val_auc"] = math.nan
        artifact_json["training"] = {
            "backend": "sklearn_lbfgs",
            "steps_run": int(getattr(clf, "n_iter_", [0])[0]),
            "l2": float(args.l2),
            "history": history,
        }
        return FittedModel(name, "logistic", list(feature_names), artifact_json, artifact, history)
    except Exception as exc:
        print(f"[stackSweep] sklearn logistic unavailable for {name}; using numpy optimizer: {exc}", flush=True)
        return fit_logistic_numpy(name, feature_names, x, y, train_mask, args, seed, val_mask=val_mask)


def fit_gbm(name, feature_names, x, y, train_mask, args, seed, val_mask=None) -> FittedModel:
    from sklearn.ensemble import GradientBoostingClassifier

    x_clean, impute = train_apply_impute(x[train_mask], x)
    clf = GradientBoostingClassifier(
        n_estimators=args.gbm_estimators,
        learning_rate=args.gbm_learning_rate,
        max_depth=args.gbm_max_depth,
        max_leaf_nodes=args.gbm_max_leaf_nodes,
        subsample=0.92,
        random_state=seed,
    )
    clf.fit(x_clean[train_mask], y[train_mask], sample_weight=class_balanced_weights(y[train_mask]))
    artifact = gbm_artifact_json(clf, feature_names, impute)
    train_weights = class_balanced_weights(y[train_mask])
    val_mask = np.zeros(len(y), dtype=bool) if val_mask is None else val_mask
    val_weights = class_balanced_weights(y[val_mask]) if val_mask.any() else np.ones(0, dtype="float64")
    history = []
    for iteration, train_prob in enumerate(clf.staged_predict_proba(x_clean[train_mask]), 1):
        train_score = np.asarray(train_prob[:, 1], dtype="float64")
        row = {
            "iteration": iteration,
            "train_loss": binary_cross_entropy(y[train_mask], train_score, train_weights),
            "train_auc": auc_score(y[train_mask], train_score),
        }
        if val_mask.any():
            # scikit-learn staged_predict_proba returns a generator, so build the
            # matching validation sequence lazily below after fitting.
            row["val_loss"] = math.nan
            row["val_auc"] = math.nan
        else:
            row["val_loss"] = math.nan
            row["val_auc"] = math.nan
        history.append(row)
    if val_mask.any():
        for row, val_prob in zip(history, clf.staged_predict_proba(x_clean[val_mask]), strict=False):
            val_score = np.asarray(val_prob[:, 1], dtype="float64")
            row["val_loss"] = binary_cross_entropy(y[val_mask], val_score, val_weights)
            row["val_auc"] = auc_score(y[val_mask], val_score)
    artifact["training"] = {
        "backend": "sklearn_gradient_boosting_classifier",
        "iterations_run": int(len(history)),
        "n_estimators": int(args.gbm_estimators),
        "learning_rate": float(args.gbm_learning_rate),
        "max_depth": int(args.gbm_max_depth),
        "max_leaf_nodes": int(args.gbm_max_leaf_nodes),
        "history": history,
    }
    return FittedModel(name, "tiny_gbm", list(feature_names), artifact, clf, history)


def binary_cross_entropy(y_true, prob, weights, l2_value: float = 0.0) -> float:
    y_true = np.asarray(y_true, dtype="float64")
    prob = np.clip(np.asarray(prob, dtype="float64"), EPS, 1.0 - EPS)
    weights = np.asarray(weights, dtype="float64")
    denom = max(float(np.sum(weights)), 1.0)
    loss = -np.sum(weights * (y_true * np.log(prob) + (1.0 - y_true) * np.log(1.0 - prob))) / denom
    return float(loss + l2_value)


def init_nn_params(input_dim: int, hidden_layers: list[int], seed: int) -> list[dict[str, np.ndarray]]:
    rng = np.random.default_rng(seed)
    dims = [input_dim] + list(hidden_layers) + [1]
    params = []
    for fan_in, fan_out in zip(dims[:-1], dims[1:]):
        scale = math.sqrt(2.0 / max(1, fan_in)) if fan_out != 1 else math.sqrt(1.0 / max(1, fan_in))
        params.append(
            {
                "weight": rng.normal(0.0, scale, size=(fan_in, fan_out)).astype("float64"),
                "bias": np.zeros(fan_out, dtype="float64"),
            }
        )
    return params


def nn_forward(xz: np.ndarray, params: list[dict[str, np.ndarray]]):
    activations = [xz]
    preacts = []
    a = xz
    for layer in params[:-1]:
        z = a @ layer["weight"] + layer["bias"]
        preacts.append(z)
        a = np.maximum(z, 0.0)
        activations.append(a)
    z = a @ params[-1]["weight"] + params[-1]["bias"]
    preacts.append(z)
    prob = sigmoid(z[:, 0])
    return prob, activations, preacts


def nn_predict_from_artifact(model: dict, x: np.ndarray) -> np.ndarray:
    impute = np.asarray(model["impute"], dtype="float64")
    mean = np.asarray(model["mean"], dtype="float64")
    scale = np.asarray(model["scale"], dtype="float64")
    scale = np.where(np.abs(scale) > 1.0e-12, scale, 1.0)
    x_clean = np.where(np.isfinite(x), x, impute)
    a = np.nan_to_num((x_clean - mean) / scale, nan=0.0, posinf=0.0, neginf=0.0)
    for layer in model["layers"][:-1]:
        a = np.maximum(a @ np.asarray(layer["weight"], dtype="float64") + np.asarray(layer["bias"], dtype="float64"), 0.0)
    out = a @ np.asarray(model["layers"][-1]["weight"], dtype="float64") + np.asarray(model["layers"][-1]["bias"], dtype="float64")
    return sigmoid(out[:, 0])


def fit_nn(name, feature_names, x, y, train_mask, val_mask, args, seed) -> FittedModel:
    hidden_layers = parse_hidden_layers(args.nn_hidden)
    x_clean, impute = train_apply_impute(x[train_mask], x)
    x_train_z, mean, scale = standardize_train_apply(x_clean[train_mask], x_clean[train_mask])
    x_all_z = np.nan_to_num((x_clean - mean) / scale, nan=0.0, posinf=0.0, neginf=0.0)
    y_train = y[train_mask].astype("float64")
    w_train = class_balanced_weights(y[train_mask])
    if val_mask is None or val_mask.sum() < 20 or len(np.unique(y[val_mask])) < 2:
        rng_split = np.random.default_rng(seed + 7919)
        train_idx_all = np.flatnonzero(train_mask)
        val_take = []
        keep_take = []
        for cls in (0, 1):
            cls_idx = train_idx_all[y[train_idx_all] == cls]
            rng_split.shuffle(cls_idx)
            n_val = max(1, int(round(0.15 * len(cls_idx)))) if len(cls_idx) >= 10 else 0
            val_take.extend(cls_idx[:n_val].tolist())
            keep_take.extend(cls_idx[n_val:].tolist())
        if val_take and keep_take:
            train_mask = np.zeros(len(y), dtype=bool)
            train_mask[np.asarray(keep_take, dtype=int)] = True
            val_mask = np.zeros(len(y), dtype=bool)
            val_mask[np.asarray(val_take, dtype=int)] = True
            y_train = y[train_mask].astype("float64")
            w_train = class_balanced_weights(y[train_mask])
    val_mask = np.zeros(len(y), dtype=bool) if val_mask is None else val_mask
    y_val = y[val_mask].astype("float64")
    w_val = class_balanced_weights(y[val_mask]) if val_mask.any() else np.ones(0, dtype="float64")
    train_indices = np.flatnonzero(train_mask)
    if train_indices.size < 20:
        raise ValueError(f"{name}: insufficient neural-stack train rows")

    params = init_nn_params(x.shape[1], hidden_layers, seed)
    m = [{"weight": np.zeros_like(layer["weight"]), "bias": np.zeros_like(layer["bias"])} for layer in params]
    v = [{"weight": np.zeros_like(layer["weight"]), "bias": np.zeros_like(layer["bias"])} for layer in params]
    best_params = [{k: arr.copy() for k, arr in layer.items()} for layer in params]
    best_val = math.inf
    best_epoch = 0
    history = []
    rng = np.random.default_rng(seed + 17)
    step = 0
    batch_size = max(128, int(args.nn_batch_size))
    lr = float(args.nn_learning_rate)
    beta1 = 0.9
    beta2 = 0.999
    for epoch in range(1, max(1, int(args.nn_epochs)) + 1):
        rng.shuffle(train_indices)
        for start in range(0, len(train_indices), batch_size):
            batch = train_indices[start : start + batch_size]
            xb = x_all_z[batch]
            yb = y[batch].astype("float64")
            wb = class_balanced_weights(y[batch])
            prob, activations, preacts = nn_forward(xb, params)
            delta = ((prob - yb) * wb / max(float(np.sum(wb)), 1.0))[:, None]
            grad_w = []
            grad_b = []
            for layer_idx in reversed(range(len(params))):
                a_prev = activations[layer_idx]
                gw = a_prev.T @ delta
                if args.nn_l2 > 0.0:
                    gw = gw + float(args.nn_l2) * params[layer_idx]["weight"]
                gb = delta.sum(axis=0)
                grad_w.append(gw)
                grad_b.append(gb)
                if layer_idx > 0:
                    delta = delta @ params[layer_idx]["weight"].T
                    delta = delta * (preacts[layer_idx - 1] > 0.0)
            grad_w.reverse()
            grad_b.reverse()
            step += 1
            for layer_idx in range(len(params)):
                for key, grad in (("weight", grad_w[layer_idx]), ("bias", grad_b[layer_idx])):
                    m[layer_idx][key] = beta1 * m[layer_idx][key] + (1.0 - beta1) * grad
                    v[layer_idx][key] = beta2 * v[layer_idx][key] + (1.0 - beta2) * (grad * grad)
                    m_hat = m[layer_idx][key] / (1.0 - beta1**step)
                    v_hat = v[layer_idx][key] / (1.0 - beta2**step)
                    params[layer_idx][key] -= lr * m_hat / (np.sqrt(v_hat) + 1.0e-8)
        train_prob = nn_forward(x_all_z[train_mask], params)[0]
        l2_value = 0.5 * float(args.nn_l2) * sum(float(np.sum(layer["weight"] ** 2)) for layer in params)
        train_loss = binary_cross_entropy(y_train, train_prob, w_train, l2_value)
        if val_mask.any():
            val_prob = nn_forward(x_all_z[val_mask], params)[0]
            val_loss = binary_cross_entropy(y_val, val_prob, w_val, l2_value)
            val_auc = auc_score(y[val_mask], val_prob)
        else:
            val_loss = train_loss
            val_auc = math.nan
        history.append({"epoch": epoch, "train_loss": train_loss, "val_loss": val_loss, "val_auc": val_auc})
        if val_loss < best_val - 1.0e-5:
            best_val = val_loss
            best_epoch = epoch
            best_params = [{k: arr.copy() for k, arr in layer.items()} for layer in params]
        elif epoch - best_epoch >= max(1, int(args.nn_patience)):
            break

    artifact = {
        "kind": "mlp",
        "feature_names": list(feature_names),
        "impute": np.asarray(impute, dtype="float64").tolist(),
        "mean": np.asarray(mean, dtype="float64").tolist(),
        "scale": np.asarray(scale, dtype="float64").tolist(),
        "activation": "relu",
        "output": "sigmoid",
        "hidden_layers": hidden_layers,
        "layers": [
            {"weight": layer["weight"].astype("float64").tolist(), "bias": layer["bias"].astype("float64").tolist()}
            for layer in best_params
        ],
        "training": {
            "epochs_run": len(history),
            "best_epoch": best_epoch,
            "best_val_loss": best_val,
            "last_train_loss": history[-1]["train_loss"] if history else None,
            "last_val_loss": history[-1]["val_loss"] if history else None,
            "last_val_auc": history[-1]["val_auc"] if history else None,
            "history": history,
            "learning_rate": float(args.nn_learning_rate),
            "l2": float(args.nn_l2),
            "batch_size": int(args.nn_batch_size),
        },
    }
    return FittedModel(name, "nn", list(feature_names), artifact, artifact)


def linear_artifact_json(model: LinearArtifact) -> dict:
    return {
        "kind": "logistic",
        "feature_names": model.feature_names,
        "impute": model.impute.tolist(),
        "mean": model.mean.tolist(),
        "scale": model.scale.tolist(),
        "coef": model.coef.tolist(),
        "intercept": model.intercept,
    }


def export_tree(tree) -> dict:
    t = tree.tree_
    values = t.value.reshape((t.node_count, -1))[:, 0]
    return {
        "children_left": t.children_left.astype("int32").tolist(),
        "children_right": t.children_right.astype("int32").tolist(),
        "feature": t.feature.astype("int32").tolist(),
        "threshold": t.threshold.astype("float64").tolist(),
        "value": values.astype("float64").tolist(),
    }


def gbm_artifact_json(clf, feature_names, impute) -> dict:
    try:
        init_raw = float(clf._raw_predict_init(np.zeros((1, len(feature_names)), dtype="float64"))[0, 0])
    except Exception:
        prior = float(getattr(clf.init_, "class_prior_", [0.5, 0.5])[1])
        init_raw = float(logit_prob(np.array([prior]))[0])
    trees = [export_tree(est[0]) for est in clf.estimators_]
    return {
        "kind": "gradient_boosting_classifier",
        "feature_names": list(feature_names),
        "impute": np.asarray(impute, dtype="float64").tolist(),
        "initial_raw_score": init_raw,
        "learning_rate": float(clf.learning_rate),
        "trees": trees,
    }


def route_mask(frame, route: RouteSpec) -> np.ndarray:
    et = np.asarray(frame["cluster_Et"], dtype="float64")
    mask = np.isfinite(et) & (et >= route.pt_lo) & (et < route.pt_hi)
    if route.cent_lo is not None and route.cent_hi is not None:
        cent = np.asarray(frame["centrality"], dtype="float64")
        mask &= np.isfinite(cent) & (cent >= route.cent_lo) & (cent < route.cent_hi)
    return mask


def pt_range_mask(frame, lo: float, hi: float) -> np.ndarray:
    et = np.asarray(frame["cluster_Et"], dtype="float64")
    return np.isfinite(et) & (et >= lo) & (et < hi)


def bins_from_edges(edges: list[float]) -> list[tuple[float, float]]:
    return [(float(edges[i]), float(edges[i + 1])) for i in range(len(edges) - 1)]


def route_specs_from_bins(pt_bins, cent_bins=None) -> tuple[RouteSpec, ...]:
    routes = []
    for plo, phi in pt_bins:
        if cent_bins is None:
            routes.append(RouteSpec(f"pt_{safe_label(plo)}_{safe_label(phi)}", plo, phi))
        else:
            for clo, chi in cent_bins:
                routes.append(
                    RouteSpec(
                        f"pt_{safe_label(plo)}_{safe_label(phi)}_cent_{safe_label(clo)}_{safe_label(chi)}",
                        plo,
                        phi,
                        clo,
                        chi,
                    )
                )
    return tuple(routes)


def route_specs_from_cent_bins(pt_lo: float, pt_hi: float, cent_bins) -> tuple[RouteSpec, ...]:
    routes = []
    for clo, chi in cent_bins:
        routes.append(
            RouteSpec(
                f"cent_{safe_label(clo)}_{safe_label(chi)}",
                float(pt_lo),
                float(pt_hi),
                float(clo),
                float(chi),
            )
        )
    return tuple(routes)


def route_payload(route: RouteSpec) -> dict:
    return {
        "label": route.label,
        "pt_lo": route.pt_lo,
        "pt_hi": route.pt_hi,
        "cent_lo": route.cent_lo,
        "cent_hi": route.cent_hi,
    }


def variant_payload(variant: VariantSpec) -> dict:
    return {
        "name": variant.name,
        "pt_lo": variant.pt_lo,
        "pt_hi": variant.pt_hi,
        "context": variant.context,
        "include_centrality": variant.include_centrality,
        "routes": [route_payload(route) for route in variant.routes],
    }


def build_variants(args) -> list[VariantSpec]:
    ranges = [(5.0, 35.0), (15.0, 35.0)]
    cent3 = [(0.0, 20.0), (20.0, 50.0), (50.0, 80.0)]
    cent7 = [(0.0, 10.0), (10.0, 20.0), (20.0, 30.0), (30.0, 40.0), (40.0, 50.0), (50.0, 60.0), (60.0, 80.0)]
    variants: list[VariantSpec] = []
    for lo, hi in ranges:
        tag = range_label(lo, hi)
        coarse_edges = [5.0, 15.0, 25.0, 35.0] if lo < 10.0 else [15.0, 20.0, 25.0, 35.0]
        fine_edges = [5.0, 10.0, 15.0, 18.0, 20.0, 22.5, 25.0, 30.0, 35.0] if lo < 10.0 else [15.0, 18.0, 20.0, 22.5, 25.0, 30.0, 35.0]
        coarse_bins = bins_from_edges(coarse_edges)
        fine_bins = bins_from_edges(fine_edges)
        variants.extend(
            [
                VariantSpec(f"global{tag}_EtCent_full", lo, hi, "et_cent", True),
                VariantSpec(f"global{tag}_EtOnly_full", lo, hi, "et_only", False),
                VariantSpec(f"global{tag}_CentOnly_full", lo, hi, "cent_only", True),
                VariantSpec(f"ptCoarse{tag}_centInput_full", lo, hi, "cent_only", True, route_specs_from_bins(coarse_bins)),
                VariantSpec(f"ptFine{tag}_centInput_full", lo, hi, "cent_only", True, route_specs_from_bins(fine_bins)),
                VariantSpec(f"cent3{tag}_EtInput_full", lo, hi, "et_only", False, route_specs_from_cent_bins(lo, hi, cent3)),
                VariantSpec(f"cent7{tag}_EtInput_full", lo, hi, "et_only", False, route_specs_from_cent_bins(lo, hi, cent7)),
                VariantSpec(f"ptCoarse{tag}_cent3_full", lo, hi, "score_only", False, route_specs_from_bins(coarse_bins, cent3)),
                VariantSpec(f"ptFine{tag}_cent3_full", lo, hi, "score_only", False, route_specs_from_bins(fine_bins, cent3)),
                VariantSpec(f"ptCoarse{tag}_cent7_full", lo, hi, "score_only", False, route_specs_from_bins(coarse_bins, cent7)),
                VariantSpec(f"ptFine{tag}_cent7_full", lo, hi, "score_only", False, route_specs_from_bins(fine_bins, cent7)),
            ]
        )
    if args.sweep == "compact":
        keep = {"global15to35_EtCent_full", "global15to35_EtOnly_full", "ptFine15to35_centInput_full", "ptFine15to35_cent7_full"}
        variants = [variant for variant in variants if variant.name in keep]
    if args.include_controls:
        variants.extend(
            [
                VariantSpec("control15to35_scoreOnly", 15.0, 35.0, "score_only", False),
                VariantSpec("control15to35_bdtOnly", 15.0, 35.0, "bdt_only", False),
                VariantSpec("control15to35_mlpOnly", 15.0, 35.0, "mlp_only", False),
            ]
        )
    requested = [item.strip() for item in getattr(args, "variants", "").split(",") if item.strip()]
    if requested:
        known = {variant.name for variant in variants}
        missing = [name for name in requested if name not in known]
        if missing:
            raise SystemExit(f"Unknown --variants entry/entries: {missing}. Known variants include: {', '.join(sorted(known))}")
        keep = set(requested)
        variants = [variant for variant in variants if variant.name in keep]
    return variants


def feature_names_for_variant(frame, variant: VariantSpec, args) -> tuple[list[str], list[str]]:
    add_interactions(frame)
    names = ["bdt_score", "bdt_is_finite", "mlp_score", "mlp_logit", "mlp_is_finite"]
    if variant.context == "bdt_only":
        names = ["bdt_score", "bdt_is_finite"]
    elif variant.context == "mlp_only":
        names = ["mlp_score", "mlp_logit", "mlp_is_finite"]
    else:
        for feature in FULL_FEATURES:
            if feature == "centrality":
                continue
            if feature not in names:
                names.append(feature)
        if variant.include_centrality and "centrality" not in names:
            names.append("centrality")
        if variant.context in ("et_cent", "et_only") and "log_cluster_Et" not in names:
            names.append("log_cluster_Et")
        if variant.context in ("et_cent", "cent_only") and variant.include_centrality and "centrality_scaled" not in names:
            names.append("centrality_scaled")
        if variant.context == "et_cent":
            names += ["bdt_x_mlp_logit", "bdt_x_logEt", "mlp_logit_x_logEt", "bdt_x_cent", "mlp_logit_x_cent", "logEt_x_cent"]
        elif variant.context == "et_only":
            names += ["bdt_x_mlp_logit", "bdt_x_logEt", "mlp_logit_x_logEt"]
        elif variant.context == "cent_only":
            names += ["bdt_x_mlp_logit", "bdt_x_cent", "mlp_logit_x_cent"]
    deduped = []
    for name in names:
        if name not in deduped:
            deduped.append(name)
    missing = [name for name in deduped if name not in frame]
    if missing and not args.allow_missing_full_features:
        raise SystemExit(
            f"Variant {variant.name} needs missing full-feature column(s): {missing}. "
            "Regenerate MLP/BDT validation score caches with the updated validators, "
            "or pass --allow-missing-full-features for a clearly limited diagnostic."
        )
    return [name for name in deduped if name in frame], missing


def feature_matrix(frame, names: list[str]) -> np.ndarray:
    return np.column_stack([np.asarray(frame[name], dtype="float64") for name in names]).astype("float64")


def stratified_split(y, seed: int, train_fraction: float, val_fraction: float):
    if train_fraction + val_fraction >= 1.0:
        raise SystemExit("train-fraction + val-fraction must leave a positive test split")
    rng = np.random.default_rng(seed)
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


def fit_one(name, algorithm, feature_names, x, y, train_mask, args, seed, val_mask=None) -> FittedModel | None:
    if train_mask.sum() < 20 or len(np.unique(y[train_mask])) < 2:
        print(f"[stackSweep][WARN] skipping {name}: insufficient train rows/classes", flush=True)
        return None
    try:
        if algorithm == "logistic":
            return fit_logistic(name, feature_names, x, y, train_mask, args, seed, val_mask=val_mask)
        if algorithm == "gbm":
            return fit_gbm(name, feature_names, x, y, train_mask, args, seed, val_mask=val_mask)
        if algorithm in ("nn", "mlp"):
            return fit_nn(name, feature_names, x, y, train_mask, val_mask, args, seed)
    except Exception as exc:
        print(f"[stackSweep][WARN] skipping {name}: {exc}", flush=True)
        return None
    raise ValueError(f"Unknown algorithm: {algorithm}")


def score_variant(frame, variant: VariantSpec, algorithm: str, masks, args, seed: int):
    y = np.asarray(frame["is_signal"], dtype="int32")
    feature_names, missing_features = feature_names_for_variant(frame, variant, args)
    x = feature_matrix(frame, feature_names)
    eval_mask = pt_range_mask(frame, variant.pt_lo, variant.pt_hi)
    score = np.full(len(y), np.nan, dtype="float64")
    route_artifacts = []
    fitted_for_pickle = {}
    histories = []
    model_name = f"{variant.name}_{algorithm}"
    if variant.routes:
        for route_idx, route in enumerate(variant.routes):
            mask = eval_mask & route_mask(frame, route)
            train_mask = masks["train"] & mask
            val_mask = masks["val"] & mask
            fitted = fit_one(f"{model_name}_{route.label}", algorithm, feature_names, x, y, train_mask, args, seed + route_idx + 1, val_mask=val_mask)
            if fitted is None:
                continue
            pred_mask = mask
            score[pred_mask] = fitted.predict(x[pred_mask])
            route_artifacts.append(
                {
                    "label": route.label,
                    "pt_lo": route.pt_lo,
                    "pt_hi": route.pt_hi,
                    "cent_lo": route.cent_lo,
                    "cent_hi": route.cent_hi,
                    "model": fitted.artifact,
                }
            )
            fitted_for_pickle[route.label] = {"algorithm": fitted.algorithm, "feature_names": fitted.feature_names, "model": fitted.model_object}
            for row in fitted.history:
                out_row = {"route": route.label, "submodel": f"{model_name}_{route.label}", "algorithm": fitted.algorithm}
                out_row.update(row)
                histories.append(out_row)
        if not route_artifacts:
            return None
        artifact = {
            "schema": "RJ_AUAU_STACKED_BDT_MLP_ARTIFACT_V1",
            "name": model_name,
            "variant": variant_payload(variant),
            "diagnostic_only": True,
            "uses_runtime_bdt_score": True,
            "algorithm": algorithm,
            "feature_names": feature_names,
            "missing_feature_columns": missing_features,
            "routes": route_artifacts,
        }
    else:
        train_mask = masks["train"] & eval_mask
        val_mask = masks["val"] & eval_mask
        fitted = fit_one(model_name, algorithm, feature_names, x, y, train_mask, args, seed, val_mask=val_mask)
        if fitted is None:
            return None
        score[eval_mask] = fitted.predict(x[eval_mask])
        artifact = {
            "schema": "RJ_AUAU_STACKED_BDT_MLP_ARTIFACT_V1",
            "name": model_name,
            "variant": variant_payload(variant),
            "diagnostic_only": True,
            "uses_runtime_bdt_score": True,
            "algorithm": algorithm,
            "feature_names": feature_names,
            "missing_feature_columns": missing_features,
            "model": fitted.artifact,
        }
        fitted_for_pickle["global"] = {"algorithm": fitted.algorithm, "feature_names": fitted.feature_names, "model": fitted.model_object}
        for row in fitted.history:
            out_row = {"route": "global", "submodel": model_name, "algorithm": fitted.algorithm}
            out_row.update(row)
            histories.append(out_row)
    return {
        "name": model_name,
        "variant": variant,
        "algorithm": algorithm,
        "score": score,
        "eval_mask": eval_mask,
        "artifact": artifact,
        "pickle": fitted_for_pickle,
        "feature_names": feature_names,
        "missing_features": missing_features,
        "histories": histories,
    }


def corrcoef_or_nan(x, y):
    x = np.asarray(x, dtype="float64")
    y = np.asarray(y, dtype="float64")
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 3 or np.std(x[mask]) <= 0.0 or np.std(y[mask]) <= 0.0:
        return math.nan
    return float(np.corrcoef(x[mask], y[mask])[0, 1])


def split_metrics(frame, score, mask, target_eff, pt_bins, cent_bins):
    y = np.asarray(frame["is_signal"][mask], dtype="int32")
    score_m = np.asarray(score[mask], dtype="float64")
    finite = np.isfinite(score_m) & np.isin(y, [0, 1])
    wp = threshold_for_signal_efficiency(y[finite], score_m[finite], target_eff) if finite.any() else None
    metrics = {
        "entries": int(mask.sum()),
        "signal_entries": int(np.sum(y == 1)),
        "background_entries": int(np.sum(y == 0)),
        "finite_fraction": float(np.mean(np.isfinite(score_m))) if len(score_m) else math.nan,
        "auc": auc_score(y, score_m),
        "wp80_fake": None if wp is None else wp["background_fake_rate"],
        "wp80_threshold": None if wp is None else wp["threshold"],
        "ece": calibration_report(y[finite], score_m[finite]).get("ece", math.nan) if finite.any() else math.nan,
        "score_vs_eiso_corr": corrcoef_or_nan(score_m, frame["reco_eiso"][mask]) if "reco_eiso" in frame else math.nan,
        "pt_bins": [],
        "centrality_bins": [],
    }
    et = np.asarray(frame["cluster_Et"][mask], dtype="float64")
    cent = np.asarray(frame["centrality"][mask], dtype="float64")
    for lo, hi in pt_bins:
        bmask = np.isfinite(et) & (et >= lo) & (et < hi)
        if not bmask.any():
            continue
        wp_bin = threshold_for_signal_efficiency(y[bmask], score_m[bmask], target_eff)
        metrics["pt_bins"].append({"lo": lo, "hi": hi, "entries": int(bmask.sum()), "auc": auc_score(y[bmask], score_m[bmask]), "wp80_fake": None if wp_bin is None else wp_bin["background_fake_rate"]})
    for lo, hi in cent_bins:
        bmask = np.isfinite(cent) & (cent >= lo) & (cent < hi)
        if not bmask.any():
            continue
        wp_bin = threshold_for_signal_efficiency(y[bmask], score_m[bmask], target_eff)
        metrics["centrality_bins"].append({"lo": lo, "hi": hi, "entries": int(bmask.sum()), "auc": auc_score(y[bmask], score_m[bmask]), "wp80_fake": None if wp_bin is None else wp_bin["background_fake_rate"]})
    return metrics


def flatten_row(model_name, variant: VariantSpec, algorithm, split_name, feature_names, metrics, highpt, anchors):
    row = {
        "model": model_name,
        "variant": variant.name,
        "algorithm": algorithm,
        "split": split_name,
        "pt_range": f"{variant.pt_lo:g}:{variant.pt_hi:g}",
        "context": variant.context,
        "routes": len(variant.routes),
        "features": "+".join(feature_names),
        "entries": metrics["entries"],
        "auc": metrics["auc"],
        "wp80_fake": metrics["wp80_fake"],
        "ece": metrics["ece"],
        "finite_fraction": metrics["finite_fraction"],
        "score_vs_eiso_corr": metrics["score_vs_eiso_corr"],
        "highpt_auc_20_35": highpt["auc"],
        "highpt_wp80_fake_20_35": highpt["wp80_fake"],
        "delta_auc_vs_broad_bdt": metrics["auc"] - anchors["bdt_inclusive"] if math.isfinite(metrics["auc"]) else math.nan,
        "delta_auc_vs_primary_mlp": metrics["auc"] - anchors["primary_mlp_inclusive"] if math.isfinite(metrics["auc"]) else math.nan,
    }
    for item in metrics["pt_bins"]:
        label = f"pt_{safe_label(item['lo'])}_{safe_label(item['hi'])}"
        row[f"{label}_auc"] = item["auc"]
        row[f"{label}_wp80_fake"] = item["wp80_fake"]
    for item in metrics["centrality_bins"]:
        label = f"cent_{safe_label(item['lo'])}_{safe_label(item['hi'])}"
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
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def rank_key(row: dict):
    fake = row.get("highpt_wp80_fake_20_35")
    high_auc = row.get("highpt_auc_20_35")
    auc = row.get("auc")
    ece = row.get("ece")
    fake = 999.0 if fake in (None, "") or not math.isfinite(float(fake)) else float(fake)
    high_auc = -999.0 if high_auc in (None, "") or not math.isfinite(float(high_auc)) else float(high_auc)
    auc = -999.0 if auc in (None, "") or not math.isfinite(float(auc)) else float(auc)
    ece = 999.0 if ece in (None, "") or not math.isfinite(float(ece)) else float(ece)
    return (fake, -high_auc, -auc, ece)


def make_synthetic_cache(outdir: Path, n: int, seed: int) -> tuple[Path, Path]:
    rng = np.random.default_rng(seed)
    y = rng.integers(0, 2, size=n).astype("int32")
    et = rng.uniform(5.0, 35.0, size=n)
    cent = rng.uniform(0.0, 80.0, size=n)
    frame = {
        "is_signal": y,
        "cluster_Et": et.astype("float32"),
        "centrality": cent.astype("float32"),
        "cluster_Eta": rng.normal(0.0, 0.35, size=n).astype("float32"),
        "vertexz": rng.normal(0.0, 4.0, size=n).astype("float32"),
        "reco_eiso": (rng.normal(0.0, 4.0, size=n) + 2.2 * (1 - y)).astype("float32"),
    }
    for idx, feature in enumerate(FULL_FEATURES):
        if feature in frame or feature in WIDTH_RATIO_FEATURES:
            continue
        signal_shift = 0.16 * y - 0.05 * (et - 20.0) / 15.0 + 0.03 * (cent - 40.0) / 40.0
        frame[feature] = (rng.normal(0.0, 1.0, size=n) + signal_shift + 0.02 * idx).astype("float32")
    add_derived_features(frame)
    linear = 1.3 * y + 0.5 * (20.0 - et) / 15.0 - 0.15 * cent / 80.0 + rng.normal(0.0, 0.9, size=n)
    mlp_score = sigmoid(linear - 0.6)
    bdt_score = sigmoid(0.9 * y - 0.45 * (et - 20.0) / 15.0 + 0.25 * frame["cluster_weta_cogx"] + rng.normal(0.0, 0.8, size=n))
    outdir.mkdir(parents=True, exist_ok=True)
    mlp_npz = outdir / "synthetic_mlp_cache_000.npz"
    bdt_npz = outdir / "synthetic_bdt_cache_000.npz"
    common = {k: np.asarray(v) for k, v in frame.items()}
    np.savez_compressed(mlp_npz, **common, **{DEFAULT_MLP_SCORE: mlp_score.astype("float32")})
    np.savez_compressed(bdt_npz, **common, **{DEFAULT_BDT_SCORE: bdt_score.astype("float32")})
    mlp_manifest = outdir / "synthetic_mlp_score_caches.list"
    bdt_manifest = outdir / "synthetic_bdt_score_caches.list"
    mlp_manifest.write_text(str(mlp_npz) + "\n")
    bdt_manifest.write_text(str(bdt_npz) + "\n")
    return mlp_manifest, bdt_manifest


def preflight(frame, cache_info, args) -> dict:
    missing_full = [feature for feature in FULL_FEATURES if feature not in frame]
    payload = {
        "schema": "RJ_AUAU_STACKED_BDT_MLP_PREFLIGHT_V1",
        "cache_info": cache_info,
        "rows": int(len(frame["is_signal"])),
        "available_full_features": [feature for feature in FULL_FEATURES if feature in frame],
        "missing_full_features": missing_full,
        "allow_missing_full_features": bool(args.allow_missing_full_features),
        "status": "OK" if (not missing_full or args.allow_missing_full_features) else "MISSING_FULL_FEATURES",
    }
    write_json(args.outdir / "stacked_sweep_preflight.json", payload)
    if missing_full and not args.allow_missing_full_features:
        raise SystemExit(
            "Full-feature stacker preflight failed; missing columns: "
            + ", ".join(missing_full)
            + ". Regenerate BDT/MLP score caches with the updated validators or use --allow-missing-full-features for a limited diagnostic."
        )
    return payload


def main():
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    if args.self_test:
        synthetic_dir = args.outdir / "synthetic_inputs"
        if synthetic_dir.exists():
            shutil.rmtree(synthetic_dir)
        args.mlp_cache, args.bdt_cache = make_synthetic_cache(synthetic_dir, 6000, args.random_seed)
        args.max_rows = 0
        args.allow_missing_full_features = False
        if not args.variants:
            args.sweep = "compact"
        args.include_controls = True

    print("[stackSweep] RECOILJETS_AUAU_STACKED_BDT_MLP_FULL_FEATURE_SWEEP_V1", flush=True)
    print(f"[stackSweep] outdir={args.outdir}", flush=True)
    frame, cache_info = load_aligned_caches(args)
    preflight_payload = preflight(frame, cache_info, args)
    if args.preflight_only:
        print(f"[stackSweep] DONE_PREFLIGHT={args.outdir / 'stacked_sweep_preflight.json'}", flush=True)
        return

    algorithms = [item.strip() for item in args.algorithms.split(",") if item.strip()]
    pt_bins = parse_bin_spec(args.report_pt_bins)
    cent_bins = parse_bin_spec(args.report_cent_bins)
    y = np.asarray(frame["is_signal"], dtype="int32")
    masks = stratified_split(y, args.random_seed, args.train_fraction, args.val_fraction)
    variants = build_variants(args)
    anchors = {
        "bdt_inclusive": args.bdt_inclusive_anchor,
        "primary_mlp_inclusive": args.primary_mlp_inclusive_anchor,
        "bdt_20_25": args.bdt_20_25_anchor,
        "bdt_25_35": args.bdt_25_35_anchor,
    }

    artifact_dir = args.outdir / "artifacts"
    pickle_dir = args.outdir / "pickles"
    history_dir = args.outdir / "training_history"
    rows: list[dict] = []
    all_history_rows: list[dict] = []
    metrics_payload = {
        "schema": "RJ_AUAU_STACKED_BDT_MLP_FULL_FEATURE_SWEEP_V1",
        "note": args.note,
        "diagnostic_only": True,
        "uses_runtime_bdt_score": True,
        "preflight": preflight_payload,
        "anchors": anchors,
        "splits": {key: int(mask.sum()) for key, mask in masks.items()},
        "models": {},
    }
    for v_idx, variant in enumerate(variants, 1):
        print(f"[stackSweep] variant {v_idx}/{len(variants)}: {variant.name}", flush=True)
        for a_idx, algorithm in enumerate(algorithms, 1):
            try:
                result = score_variant(frame, variant, algorithm, masks, args, args.random_seed + 1000 * v_idx + 17 * a_idx)
            except Exception as exc:
                print(f"[stackSweep][WARN] failed {variant.name}/{algorithm}: {exc}", flush=True)
                continue
            if result is None:
                continue
            model_name = result["name"]
            artifact_path = artifact_dir / f"{model_name}.json"
            pickle_path = pickle_dir / f"{model_name}.pkl"
            history_path = history_dir / f"{model_name}.training_history.csv"
            write_json(artifact_path, result["artifact"])
            history_rows = []
            for row in result["histories"]:
                hist_row = {"model": model_name, "variant": variant.name}
                hist_row.update(row)
                history_rows.append(hist_row)
            write_csv(history_path, history_rows)
            all_history_rows.extend(history_rows)
            pickle_path.parent.mkdir(parents=True, exist_ok=True)
            with pickle_path.open("wb") as handle:
                pickle.dump(result["pickle"], handle)
            metrics_payload["models"][model_name] = {
                "variant": variant_payload(variant),
                "algorithm": algorithm,
                "artifact": str(artifact_path),
                "pickle": str(pickle_path),
                "training_history_csv": str(history_path),
                "features": result["feature_names"],
                "missing_features": result["missing_features"],
                "metrics": {},
            }
            for split_name in ("val", "test", "all"):
                split_mask = masks[split_name] & result["eval_mask"]
                highpt_mask = split_mask & pt_range_mask(frame, max(20.0, variant.pt_lo), min(35.0, variant.pt_hi))
                metrics = split_metrics(frame, result["score"], split_mask, args.target_signal_efficiency, pt_bins, cent_bins)
                highpt = split_metrics(frame, result["score"], highpt_mask, args.target_signal_efficiency, [(20.0, 25.0), (25.0, 35.0)], cent_bins)
                metrics_payload["models"][model_name]["metrics"][split_name] = {"inclusive": metrics, "highpt_20_35": highpt}
                rows.append(flatten_row(model_name, variant, algorithm, split_name, result["feature_names"], metrics, highpt, anchors))
            test = metrics_payload["models"][model_name]["metrics"]["test"]["inclusive"]
            high = metrics_payload["models"][model_name]["metrics"]["test"]["highpt_20_35"]
            print(
                f"[stackSweep] model={model_name} test_auc={test['auc']:.5f} "
                f"test_wp80_fake={test['wp80_fake']} highpt_auc={high['auc']:.5f} highpt_fake={high['wp80_fake']}",
                flush=True,
            )

    rank_path = args.outdir / "stacked_sweep_rank_table.csv"
    write_csv(rank_path, rows)
    write_csv(args.outdir / "stacked_sweep_training_history.csv", all_history_rows)
    write_json(args.outdir / "stacked_sweep_metrics.json", metrics_payload)
    test_rows = [row for row in rows if row["split"] == "test"]
    test_rows.sort(key=rank_key)
    top = test_rows[: max(1, args.top_n)]
    top_payload = {
        "schema": "RJ_AUAU_STACKED_BDT_MLP_TOP_VARIANTS_V1",
        "rank_table": str(rank_path),
        "top_n": args.top_n,
        "variants": top,
        "artifacts": [
            str((artifact_dir / f"{row['model']}.json"))
            for row in top
        ],
    }
    write_json(args.outdir / "stacked_sweep_top4.json", top_payload)
    print(f"[stackSweep] DONE_STACKED_SWEEP_OUTDIR={args.outdir}", flush=True)
    print(f"[stackSweep] DONE_STACKED_SWEEP_RANK_TABLE={rank_path}", flush=True)
    if top:
        best = top[0]
        print(
            f"[stackSweep] BEST_STACKED_SWEEP_MODEL={best['model']} "
            f"TEST_AUC={best['auc']} HIGHPT_FAKE={best['highpt_wp80_fake_20_35']}",
            flush=True,
        )


if __name__ == "__main__":
    main()
