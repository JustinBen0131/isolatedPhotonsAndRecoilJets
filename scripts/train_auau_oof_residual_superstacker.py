#!/usr/bin/env python3
"""Train an out-of-fold NN super-stacker for AuAu photon ID.

This is a diagnostic ceiling test.  It can either train the historical residual
NN on BDT/MLP score inputs, or a stricter tri-score superlearner: fold-train
logistic/GBM/NN base learners on the full photon feature family, evaluate each
base learner on held-out folds, then train a heavier final MLP on those honest
out-of-fold scores.  The locked test split is evaluated with base learners
trained only on the train+validation split.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import train_auau_stacked_bdt_mlp_sweep as stack  # noqa: E402
from train_auau_photon_mlp import (  # noqa: E402
    auc_score,
    calibration_report,
    threshold_for_signal_efficiency,
)

EPS = 1.0e-6
DEFAULT_PT_BINS = "15,18,20,22.5,25,30,35"
DEFAULT_CENT_BINS = "0,20,50,80"


@dataclass
class ResidualNN:
    feature_names: list[str]
    impute: np.ndarray
    mean: np.ndarray
    scale: np.ndarray
    layers: list[dict[str, np.ndarray]]
    base_score_name: str
    correction_scale: float
    history: list[dict]


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mlp-cache", type=Path)
    ap.add_argument("--bdt-cache", type=Path)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--mlp-score", default=stack.DEFAULT_MLP_SCORE)
    ap.add_argument("--bdt-score", default=stack.DEFAULT_BDT_SCORE)
    ap.add_argument("--pt-min", type=float, default=15.0)
    ap.add_argument("--pt-max", type=float, default=35.0)
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument("--lower-algorithms", default="logistic,gbm,nn")
    ap.add_argument("--base-score", default="nn", help="Lower stacker used as residual base, e.g. nn or lower_nn.")
    ap.add_argument("--lower-feature-mode", choices=["score_context", "full_features", "full_features_plus_scores"], default="score_context")
    ap.add_argument("--super-feature-mode", choices=["scores_context", "scores_plus_full_features", "scores_only"], default="scores_context")
    ap.add_argument("--final-mode", choices=["residual", "direct"], default="residual")
    ap.add_argument("--include-isolation-context", action="store_true")
    ap.add_argument("--model-name", default="")
    ap.add_argument("--hidden", default="16")
    ap.add_argument("--epochs", type=int, default=180)
    ap.add_argument("--patience", type=int, default=28)
    ap.add_argument("--batch-size", type=int, default=32768)
    ap.add_argument("--learning-rate", type=float, default=1.2e-3)
    ap.add_argument("--l2", type=float, default=2.0e-3)
    ap.add_argument("--residual-l2", type=float, default=3.0e-2)
    ap.add_argument("--correction-scale", type=float, default=0.55)
    ap.add_argument("--max-shards", type=int, default=0)
    ap.add_argument("--max-rows", type=int, default=0)
    ap.add_argument("--require-full-stat", action="store_true")
    ap.add_argument("--expected-shards", type=int, default=80)
    ap.add_argument("--train-fraction", type=float, default=0.60)
    ap.add_argument("--val-fraction", type=float, default=0.20)
    ap.add_argument("--random-seed", type=int, default=24681357)
    ap.add_argument("--target-signal-efficiency", type=float, default=0.80)
    ap.add_argument("--report-pt-bins", default=DEFAULT_PT_BINS)
    ap.add_argument("--report-cent-bins", default=DEFAULT_CENT_BINS)
    ap.add_argument("--correlation-max-rows", type=int, default=600000)
    ap.add_argument("--correlation-top-n", type=int, default=80)
    ap.add_argument("--skip-correlation-diagnostics", action="store_true")
    ap.add_argument("--self-test", action="store_true")
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


def parse_edges(text: str) -> list[float]:
    vals = [float(tok) for tok in text.replace(";", ",").split(",") if tok.strip()]
    if len(vals) < 2 or any(vals[i] >= vals[i + 1] for i in range(len(vals) - 1)):
        raise SystemExit(f"Bad bin edges: {text}")
    return vals


def edge_pairs(text: str) -> list[tuple[float, float]]:
    vals = parse_edges(text)
    return list(zip(vals[:-1], vals[1:], strict=True))


def sigmoid(x):
    return stack.sigmoid(x)


def logit_prob(x):
    return stack.logit_prob(x)


def finite_auc(y, score) -> float:
    return auc_score(np.asarray(y, dtype="int32"), np.asarray(score, dtype="float64"))


def corrcoef_or_nan(x, y):
    x = np.asarray(x, dtype="float64")
    y = np.asarray(y, dtype="float64")
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 3 or np.std(x[mask]) <= 0.0 or np.std(y[mask]) <= 0.0:
        return math.nan
    return float(np.corrcoef(x[mask], y[mask])[0, 1])


def make_self_test_frame(seed: int = 12345, n: int = 9000):
    rng = np.random.default_rng(seed)
    et = rng.uniform(15.0, 35.0, size=n)
    cent = rng.uniform(0.0, 80.0, size=n)
    latent = rng.normal(0.0, 1.0, size=n) + 0.35 * (et - 22.0) / 8.0 - 0.20 * cent / 80.0
    y = (latent + rng.normal(0.0, 0.85, size=n) > 0.0).astype("int32")
    bdt = sigmoid(0.7 * latent + rng.normal(0.0, 0.75, size=n))
    mlp = sigmoid(0.55 * latent + 0.15 * np.sin(et / 3.0) + rng.normal(0.0, 0.85, size=n))
    frame = {
        "is_signal": y,
        "cluster_Et": et,
        "centrality": cent,
        "reco_eiso": rng.normal(0.0, 1.0, size=n) - 0.25 * latent,
        "bdt_score": bdt,
        "mlp_score": mlp,
    }
    stack.add_interactions(frame)
    return frame, {"rows": n, "shards": 0, "mlp_score_key": "synthetic_mlp", "bdt_score_key": "synthetic_bdt"}


def load_frame(args):
    if args.self_test:
        return make_self_test_frame(args.random_seed)
    if args.mlp_cache is None or args.bdt_cache is None:
        raise SystemExit("--mlp-cache and --bdt-cache are required unless --self-test is used")
    load_args = argparse.Namespace(
        mlp_cache=args.mlp_cache,
        bdt_cache=args.bdt_cache,
        logreg_cache=None,
        mlp_score=args.mlp_score,
        bdt_score=args.bdt_score,
        logreg_score=stack.DEFAULT_LOGREG_SCORE,
        max_shards=args.max_shards,
        max_rows=args.max_rows,
        require_full_stat=args.require_full_stat,
        expected_shards=args.expected_shards,
        random_seed=args.random_seed,
    )
    return stack.load_aligned_caches(load_args)


def selected_mask(frame, args):
    et = np.asarray(frame["cluster_Et"], dtype="float64")
    cent = np.asarray(frame["centrality"], dtype="float64")
    return np.isfinite(et) & np.isfinite(cent) & (et >= args.pt_min) & (et < args.pt_max) & (cent >= 0.0) & (cent < 80.0)


def dedup_existing(frame, names: list[str]) -> list[str]:
    out = []
    for name in names:
        if name in frame and name not in out:
            out.append(name)
    return out


def score_context_names(frame) -> list[str]:
    stack.add_interactions(frame)
    return dedup_existing(
        frame,
        [
            "bdt_score",
            "mlp_score",
            "mlp_logit",
            "log_cluster_Et",
            "centrality_scaled",
            "bdt_x_mlp_logit",
            "bdt_x_logEt",
            "mlp_logit_x_logEt",
            "bdt_x_cent",
            "mlp_logit_x_cent",
            "logEt_x_cent",
        ],
    )


def full_feature_names(frame) -> list[str]:
    stack.add_interactions(frame)
    names = []
    for feature in stack.FULL_FEATURES:
        if feature in frame and feature not in names:
            names.append(feature)
    for feature in ("log_cluster_Et", "centrality", "centrality_scaled"):
        if feature in frame and feature not in names:
            names.append(feature)
    return names


def full_feature_names_for_args(frame, args=None) -> list[str]:
    names = full_feature_names(frame)
    if args is not None and getattr(args, "include_isolation_context", False):
        for feature in stack.ISOLATION_CONTEXT_FEATURES:
            if feature in frame and feature not in names:
                names.append(feature)
    return names


def lower_feature_names(frame, mode: str, args=None) -> list[str]:
    if mode == "score_context":
        return score_context_names(frame)
    if mode == "full_features":
        return full_feature_names_for_args(frame, args)
    if mode == "full_features_plus_scores":
        return dedup_existing(frame, score_context_names(frame) + full_feature_names_for_args(frame, args))
    raise ValueError(f"Unknown lower feature mode: {mode}")


def historical_stack_feature_names(frame) -> list[str]:
    names = [
        "bdt_score",
        "mlp_score",
        "mlp_logit",
        "log_cluster_Et",
        "centrality_scaled",
        "bdt_x_mlp_logit",
        "bdt_x_logEt",
        "mlp_logit_x_logEt",
        "bdt_x_cent",
        "mlp_logit_x_cent",
        "logEt_x_cent",
    ]
    return dedup_existing(frame, names)


def matrix_from_frame(frame, names: list[str]) -> np.ndarray:
    return np.column_stack([np.asarray(frame[name], dtype="float64") for name in names]).astype("float64")


def stable_unit_interval(seed: int, *parts) -> float:
    h = hashlib.blake2b(digest_size=8)
    h.update(str(seed).encode())
    for part in parts:
        h.update(b"|")
        h.update(str(part).encode())
    return int.from_bytes(h.digest(), "big") / float(1 << 64)


def stable_split(frame, y, seed: int, train_fraction: float, val_fraction: float):
    if train_fraction + val_fraction >= 1.0:
        raise SystemExit("train-fraction + val-fraction must leave a positive test split")
    if "run" in frame and "evt" in frame:
        run = np.asarray(frame["run"])
        evt = np.asarray(frame["evt"])
        u = np.asarray([stable_unit_interval(seed, int(r), int(e)) for r, e in zip(run, evt, strict=False)])
        test_fraction = 1.0 - train_fraction - val_fraction
        test = u < test_fraction
        val = (u >= test_fraction) & (u < test_fraction + val_fraction)
        train = ~(test | val)
        return {"train": train, "val": val, "test": test, "all": np.ones(len(y), dtype=bool), "mode": "event_hash_run_evt"}
    split = stack.stratified_split(y, seed, train_fraction, val_fraction)
    split["mode"] = "row_stratified"
    return split


def make_folds(frame, y, mask, n_folds: int, seed: int) -> list[np.ndarray]:
    if n_folds < 2:
        raise SystemExit("--folds must be >= 2")
    if "run" in frame and "evt" in frame:
        run = np.asarray(frame["run"])
        evt = np.asarray(frame["evt"])
        fold_id = np.full(len(y), -1, dtype="int32")
        idx = np.flatnonzero(mask)
        for row in idx:
            fold_id[row] = int(stable_unit_interval(seed, int(run[row]), int(evt[row])) * n_folds) % n_folds
        folds = [(fold_id == fold) for fold in range(n_folds)]
        if any(fold.sum() == 0 for fold in folds):
            raise SystemExit("At least one event-hash fold is empty")
        return folds
    rng = np.random.default_rng(seed)
    fold_id = np.full(len(y), -1, dtype="int32")
    for cls in (0, 1):
        idx = np.flatnonzero(mask & (y == cls))
        rng.shuffle(idx)
        for i, row in enumerate(idx):
            fold_id[row] = i % n_folds
    folds = [(fold_id == fold) for fold in range(n_folds)]
    if any(fold.sum() == 0 for fold in folds):
        raise SystemExit("At least one fold is empty")
    return folds


def lower_args(args) -> argparse.Namespace:
    return argparse.Namespace(
        l2=2.0e-3,
        linear_backend="numpy",
        max_linear_steps=1100,
        gbm_estimators=80,
        gbm_learning_rate=0.045,
        gbm_max_depth=3,
        gbm_max_leaf_nodes=8,
        nn_hidden="48,24",
        nn_epochs=130,
        nn_patience=18,
        nn_batch_size=max(4096, args.batch_size),
        nn_learning_rate=1.4e-3,
        nn_l2=1.5e-3,
    )


def train_lower_models_oof(frame, y, trainval_mask, test_mask, args):
    names = lower_feature_names(frame, args.lower_feature_mode, args)
    if not names:
        raise SystemExit(f"No lower-model features available for mode {args.lower_feature_mode}")
    x = matrix_from_frame(frame, names)
    algorithms = [tok.strip() for tok in args.lower_algorithms.split(",") if tok.strip()]
    if not algorithms:
        raise SystemExit("--lower-algorithms is empty")
    folds = make_folds(frame, y, trainval_mask, args.folds, args.random_seed + 1001)
    largs = lower_args(args)
    oof_scores = {alg: np.full(len(y), np.nan, dtype="float64") for alg in algorithms}
    final_scores = {alg: np.full(len(y), np.nan, dtype="float64") for alg in algorithms}
    final_artifacts = {}
    fold_rows = []
    for alg in algorithms:
        print(f"[superStack] lower algorithm={alg} folds={args.folds}", flush=True)
        for ifold, fold_mask in enumerate(folds):
            train_mask = trainval_mask & ~fold_mask
            val_mask = fold_mask
            fitted = stack.fit_one(
                f"lower_{alg}_fold{ifold}",
                alg,
                names,
                x,
                y,
                train_mask,
                largs,
                args.random_seed + 2000 + ifold * 17 + len(alg),
                val_mask=val_mask,
            )
            if fitted is None:
                raise SystemExit(f"Failed lower {alg} fold {ifold}")
            oof_scores[alg][fold_mask] = fitted.predict(x[fold_mask])
            fold_rows.append(
                {
                    "algorithm": alg,
                    "fold": ifold,
                    "train_entries": int(train_mask.sum()),
                    "oof_entries": int(fold_mask.sum()),
                    "oof_auc": finite_auc(y[fold_mask], oof_scores[alg][fold_mask]),
                }
            )
        final = stack.fit_one(
            f"lower_{alg}_final_trainval",
            alg,
            names,
            x,
            y,
            trainval_mask,
            largs,
            args.random_seed + 3000 + len(alg),
            val_mask=test_mask,
        )
        if final is None:
            raise SystemExit(f"Failed lower {alg} final model")
        final_scores[alg][test_mask] = final.predict(x[test_mask])
        final_scores[alg][trainval_mask] = oof_scores[alg][trainval_mask]
        final_artifacts[alg] = final.artifact
    return names, oof_scores, final_scores, final_artifacts, fold_rows


def make_super_features(frame, lower_scores: dict[str, np.ndarray], algorithms: list[str], mode: str, args=None) -> tuple[list[str], np.ndarray]:
    cols = []
    names = []
    for alg in algorithms:
        score_name = f"lower_{alg}"
        score = np.asarray(lower_scores[alg], dtype="float64")
        names.append(score_name)
        cols.append(score)
        names.append(f"{score_name}_logit")
        cols.append(logit_prob(score))
    if len(algorithms) >= 2:
        score_stack = np.column_stack([np.asarray(lower_scores[alg], dtype="float64") for alg in algorithms])
        names.append("lower_score_mean")
        cols.append(np.nanmean(np.where(np.isfinite(score_stack), score_stack, np.nan), axis=1))
        names.append("lower_score_spread")
        cols.append(np.nanstd(np.where(np.isfinite(score_stack), score_stack, np.nan), axis=1))
    if "log_cluster_Et" in frame:
        log_et = np.asarray(frame["log_cluster_Et"], dtype="float64")
        for alg in algorithms:
            names.append(f"lower_{alg}_x_logEt")
            cols.append(np.asarray(lower_scores[alg], dtype="float64") * log_et)
    if "centrality_scaled" in frame:
        cent = np.asarray(frame["centrality_scaled"], dtype="float64")
        for alg in algorithms:
            names.append(f"lower_{alg}_x_cent")
            cols.append(np.asarray(lower_scores[alg], dtype="float64") * cent)
    if mode == "scores_context":
        extra = score_context_names(frame)
    elif mode == "scores_plus_full_features":
        extra = full_feature_names_for_args(frame, args)
    elif mode == "scores_only":
        extra = []
    else:
        raise ValueError(f"Unknown super feature mode: {mode}")
    for name in extra:
        if name not in names:
            names.append(name)
            cols.append(np.asarray(frame[name], dtype="float64"))
    if not cols:
        raise SystemExit(f"No final-superlearner features available for mode {mode}")
    return names, np.column_stack(cols).astype("float64")


def impute_standardize(x, train_mask):
    x_train = x[train_mask]
    impute = np.nanmedian(np.where(np.isfinite(x_train), x_train, np.nan), axis=0)
    impute = np.where(np.isfinite(impute), impute, 0.0)
    clean = np.where(np.isfinite(x), x, impute)
    mean = np.nanmean(clean[train_mask], axis=0)
    scale = np.nanstd(clean[train_mask], axis=0)
    mean = np.where(np.isfinite(mean), mean, 0.0)
    scale = np.where(np.isfinite(scale) & (scale > 1.0e-9), scale, 1.0)
    return np.nan_to_num((clean - mean) / scale, nan=0.0, posinf=0.0, neginf=0.0), impute, mean, scale


def parse_hidden(text: str) -> list[int]:
    vals = [int(tok) for tok in text.replace("x", ",").split(",") if tok.strip()]
    if not vals or any(v <= 0 for v in vals):
        raise SystemExit(f"Bad hidden spec: {text}")
    return vals


def init_layers(dim: int, hidden: list[int], seed: int):
    rng = np.random.default_rng(seed)
    dims = [dim] + hidden + [1]
    layers = []
    for fan_in, fan_out in zip(dims[:-1], dims[1:], strict=True):
        scale = math.sqrt(2.0 / max(1, fan_in)) if fan_out != 1 else math.sqrt(1.0 / max(1, fan_in))
        layers.append({"weight": rng.normal(0.0, scale, size=(fan_in, fan_out)), "bias": np.zeros(fan_out)})
    return layers


def nn_raw(xz, layers):
    a = xz
    activations = [a]
    masks = []
    for layer in layers[:-1]:
        z = a @ layer["weight"] + layer["bias"]
        mask = z > 0.0
        masks.append(mask)
        a = np.where(mask, z, 0.0)
        activations.append(a)
    raw = (a @ layers[-1]["weight"] + layers[-1]["bias"])[:, 0]
    return raw, activations, masks


def residual_score(base_score, delta_raw, correction_scale):
    return sigmoid(logit_prob(base_score) + correction_scale * delta_raw)


def weighted_bce(y, p, w):
    p = np.clip(p, EPS, 1.0 - EPS)
    return float(-np.sum(w * (y * np.log(p) + (1.0 - y) * np.log(1.0 - p))) / max(float(np.sum(w)), 1.0))


def train_residual_nn(x, y, base_score, train_mask, val_mask, args, feature_names):
    xz, impute, mean, scale = impute_standardize(x, train_mask)
    hidden = parse_hidden(args.hidden)
    layers = init_layers(xz.shape[1], hidden, args.random_seed + 5000)
    m = [{"weight": np.zeros_like(l["weight"]), "bias": np.zeros_like(l["bias"])} for l in layers]
    v = [{"weight": np.zeros_like(l["weight"]), "bias": np.zeros_like(l["bias"])} for l in layers]
    beta1, beta2 = 0.9, 0.999
    train_idx = np.flatnonzero(train_mask)
    val_idx = np.flatnonzero(val_mask)
    w_all = stack.class_balanced_weights(y)
    rng = np.random.default_rng(args.random_seed + 6000)
    best = None
    best_epoch = 0
    bad = 0
    history = []
    step = 0
    for epoch in range(1, args.epochs + 1):
        rng.shuffle(train_idx)
        for start in range(0, len(train_idx), args.batch_size):
            idx = train_idx[start : start + args.batch_size]
            if len(idx) == 0:
                continue
            step += 1
            raw, acts, relu_masks = nn_raw(xz[idx], layers)
            pred = residual_score(base_score[idx], raw, args.correction_scale)
            err = (pred - y[idx].astype("float64")) * w_all[idx] / max(float(np.mean(w_all[idx])), 1.0)
            # derivative wrt residual raw; correction scale limits correction capacity.
            d_raw = err * args.correction_scale / max(1, len(idx))
            d_raw += (args.residual_l2 / max(1, len(idx))) * raw
            grads = [None for _ in layers]
            dz = d_raw[:, None]
            for ilayer in reversed(range(len(layers))):
                a_prev = acts[ilayer]
                grad_w = a_prev.T @ dz + args.l2 * layers[ilayer]["weight"]
                grad_b = np.sum(dz, axis=0)
                grads[ilayer] = {"weight": grad_w, "bias": grad_b}
                if ilayer > 0:
                    dz = (dz @ layers[ilayer]["weight"].T) * relu_masks[ilayer - 1]
            for ilayer, grad in enumerate(grads):
                for key in ("weight", "bias"):
                    m[ilayer][key] = beta1 * m[ilayer][key] + (1.0 - beta1) * grad[key]
                    v[ilayer][key] = beta2 * v[ilayer][key] + (1.0 - beta2) * (grad[key] * grad[key])
                    mh = m[ilayer][key] / (1.0 - beta1**step)
                    vh = v[ilayer][key] / (1.0 - beta2**step)
                    layers[ilayer][key] -= args.learning_rate * mh / (np.sqrt(vh) + 1.0e-8)
        train_raw, _, _ = nn_raw(xz[train_mask], layers)
        val_raw, _, _ = nn_raw(xz[val_mask], layers)
        train_pred = residual_score(base_score[train_mask], train_raw, args.correction_scale)
        val_pred = residual_score(base_score[val_mask], val_raw, args.correction_scale)
        row = {
            "epoch": epoch,
            "train_loss": weighted_bce(y[train_mask], train_pred, w_all[train_mask]),
            "val_loss": weighted_bce(y[val_mask], val_pred, w_all[val_mask]),
            "train_auc": finite_auc(y[train_mask], train_pred),
            "val_auc": finite_auc(y[val_mask], val_pred),
            "mean_abs_train_delta_raw": float(np.mean(np.abs(train_raw))),
            "mean_abs_val_delta_raw": float(np.mean(np.abs(val_raw))),
        }
        history.append(row)
        if epoch == 1 or epoch % 10 == 0:
            print(
                f"[superStack] epoch={epoch} train_loss={row['train_loss']:.5f} "
                f"val_loss={row['val_loss']:.5f} val_auc={row['val_auc']:.5f}",
                flush=True,
            )
        score = row["val_loss"]
        if best is None or score < best[0]:
            best = (score, json_ready(layers))
            best_epoch = epoch
            bad = 0
        else:
            bad += 1
            if bad >= args.patience:
                print(f"[superStack] early_stop epoch={epoch} best_epoch={best_epoch}", flush=True)
                break
    if best is not None:
        layers = [
            {"weight": np.asarray(layer["weight"], dtype="float64"), "bias": np.asarray(layer["bias"], dtype="float64")}
            for layer in best[1]
        ]
    return ResidualNN(feature_names, impute, mean, scale, layers, args.base_score, args.correction_scale, history)


def predict_residual_nn(model: ResidualNN, x, base_score):
    clean = np.where(np.isfinite(x), x, model.impute)
    xz = np.nan_to_num((clean - model.mean) / model.scale, nan=0.0, posinf=0.0, neginf=0.0)
    raw, _, _ = nn_raw(xz, model.layers)
    return residual_score(base_score, raw, model.correction_scale), raw


def final_nn_args(args) -> argparse.Namespace:
    return argparse.Namespace(
        nn_hidden=args.hidden,
        nn_epochs=args.epochs,
        nn_patience=args.patience,
        nn_batch_size=args.batch_size,
        nn_learning_rate=args.learning_rate,
        nn_l2=args.l2,
        linear_backend="numpy",
        max_linear_steps=1200,
        l2=args.l2,
        gbm_estimators=80,
        gbm_learning_rate=0.045,
        gbm_max_depth=3,
        gbm_max_leaf_nodes=8,
    )


def split_metrics(frame, y, score, mask, args) -> dict:
    score = np.asarray(score, dtype="float64")
    yy = y[mask]
    ss = score[mask]
    finite = np.isfinite(ss) & np.isin(yy, [0, 1])
    wp = threshold_for_signal_efficiency(yy[finite], ss[finite], args.target_signal_efficiency) if finite.any() else None
    out = {
        "entries": int(mask.sum()),
        "signal_entries": int(np.sum(yy == 1)),
        "background_entries": int(np.sum(yy == 0)),
        "finite_fraction": float(np.mean(np.isfinite(ss))) if ss.size else math.nan,
        "auc": finite_auc(yy, ss),
        "wp80_fake": None if wp is None else wp["background_fake_rate"],
        "wp80_threshold": None if wp is None else wp["threshold"],
        "ece": calibration_report(yy[finite], ss[finite]).get("ece", math.nan) if finite.any() else math.nan,
        "score_vs_eiso_corr": corrcoef_or_nan(ss, frame["reco_eiso"][mask]) if "reco_eiso" in frame else math.nan,
        "pt_bins": [],
        "centrality_bins": [],
    }
    et = np.asarray(frame["cluster_Et"][mask], dtype="float64")
    cent = np.asarray(frame["centrality"][mask], dtype="float64")
    for lo, hi in edge_pairs(args.report_pt_bins):
        bmask = np.isfinite(et) & (et >= lo) & (et < hi)
        if bmask.any():
            wpb = threshold_for_signal_efficiency(yy[bmask], ss[bmask], args.target_signal_efficiency)
            out["pt_bins"].append({"lo": lo, "hi": hi, "entries": int(bmask.sum()), "auc": finite_auc(yy[bmask], ss[bmask]), "wp80_fake": None if wpb is None else wpb["background_fake_rate"]})
    for lo, hi in edge_pairs(args.report_cent_bins):
        bmask = np.isfinite(cent) & (cent >= lo) & (cent < hi)
        if bmask.any():
            wpb = threshold_for_signal_efficiency(yy[bmask], ss[bmask], args.target_signal_efficiency)
            out["centrality_bins"].append({"lo": lo, "hi": hi, "entries": int(bmask.sum()), "auc": finite_auc(yy[bmask], ss[bmask]), "wp80_fake": None if wpb is None else wpb["background_fake_rate"]})
    return out


def highpt_mask(frame):
    et = np.asarray(frame["cluster_Et"], dtype="float64")
    return np.isfinite(et) & (et >= 20.0) & (et < 35.0)


def write_rank_table(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_rows(path: Path, rows: list[dict]) -> None:
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
        writer.writerows(rows)


def sample_mask_for_correlations(mask: np.ndarray, max_rows: int, seed: int) -> np.ndarray:
    idx = np.flatnonzero(mask)
    if max_rows > 0 and idx.size > max_rows:
        rng = np.random.default_rng(seed)
        idx = rng.choice(idx, size=max_rows, replace=False)
        idx.sort()
    out = np.zeros(len(mask), dtype=bool)
    out[idx] = True
    return out


def clean_matrix_for_correlation(x: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    clean = np.asarray(x, dtype="float64").copy()
    clean[~np.isfinite(clean)] = np.nan
    impute = np.nanmedian(clean, axis=0)
    impute = np.where(np.isfinite(impute), impute, 0.0)
    inds = np.where(~np.isfinite(clean))
    if inds[0].size:
        clean[inds] = impute[inds[1]]
    mean = np.mean(clean, axis=0)
    scale = np.std(clean, axis=0)
    scale = np.where(np.isfinite(scale) & (scale > 1.0e-12), scale, 1.0)
    z = (clean - mean) / scale
    z = np.nan_to_num(z, nan=0.0, posinf=0.0, neginf=0.0)
    return z, mean, scale


def correlation_vector(a, b) -> float:
    a = np.asarray(a, dtype="float64")
    b = np.asarray(b, dtype="float64")
    mask = np.isfinite(a) & np.isfinite(b)
    if mask.sum() < 3:
        return math.nan
    aa = a[mask]
    bb = b[mask]
    if np.std(aa) <= 1.0e-12 or np.std(bb) <= 1.0e-12:
        return math.nan
    return float(np.corrcoef(aa, bb)[0, 1])


def write_correlation_heatmap(path: Path, names: list[str], corr: np.ndarray) -> bool:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:
        print(f"[superStack][WARN] matplotlib unavailable; skipping correlation heatmap: {exc}", flush=True)
        return False
    n = len(names)
    fig_w = max(9.0, min(22.0, 0.32 * n + 4.0))
    fig_h = max(7.5, min(22.0, 0.32 * n + 3.0))
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=170)
    im = ax.imshow(corr, vmin=-1.0, vmax=1.0, cmap="coolwarm", interpolation="nearest")
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(names, rotation=90, fontsize=6)
    ax.set_yticklabels(names, fontsize=6)
    ax.set_title("OOF superstacker input correlations", fontsize=12)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Pearson r")
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path)
    plt.close(fig)
    return True


def write_correlation_diagnostics(
    outdir: Path,
    prefix: str,
    frame,
    y: np.ndarray,
    super_names: list[str],
    super_x: np.ndarray,
    score_map: dict[str, np.ndarray],
    final_name: str,
    selected_mask_values: np.ndarray,
    args,
) -> dict:
    sample_mask = sample_mask_for_correlations(selected_mask_values, args.correlation_max_rows, args.random_seed + 424242)
    idx = np.flatnonzero(sample_mask)
    if idx.size < 10:
        return {"status": "SKIPPED", "reason": "too_few_rows", "rows": int(idx.size)}

    target_columns = {
        "truth_is_signal": np.asarray(y, dtype="float64"),
        "cluster_Et": np.asarray(frame["cluster_Et"], dtype="float64") if "cluster_Et" in frame else None,
        "centrality": np.asarray(frame["centrality"], dtype="float64") if "centrality" in frame else None,
        "reco_eiso": np.asarray(frame["reco_eiso"], dtype="float64") if "reco_eiso" in frame else None,
    }
    for name, score in score_map.items():
        target_columns[f"score_{name}"] = np.asarray(score, dtype="float64")
    target_columns = {name: col for name, col in target_columns.items() if col is not None}

    matrix_names = list(super_names)
    matrix_values = [np.asarray(super_x[:, i], dtype="float64") for i in range(super_x.shape[1])]
    for name, col in target_columns.items():
        if name not in matrix_names:
            matrix_names.append(name)
            matrix_values.append(col)

    matrix = np.column_stack([col[idx] for col in matrix_values])
    z, _, _ = clean_matrix_for_correlation(matrix)
    corr = (z.T @ z) / max(1, z.shape[0] - 1)
    corr = np.clip(corr, -1.0, 1.0)

    matrix_path = outdir / f"{prefix}_input_correlation_matrix.csv"
    matrix_path.parent.mkdir(parents=True, exist_ok=True)
    with matrix_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["feature"] + matrix_names)
        for name, row in zip(matrix_names, corr, strict=True):
            writer.writerow([name] + [f"{value:.8g}" for value in row])

    target_rows = []
    for i, feature in enumerate(super_names):
        values = np.asarray(super_x[:, i], dtype="float64")[idx]
        row = {"feature": feature}
        for target_name, target_values in target_columns.items():
            row[f"pearson_{target_name}"] = correlation_vector(values, np.asarray(target_values, dtype="float64")[idx])
        target_rows.append(row)
    target_path = outdir / f"{prefix}_input_target_correlations.csv"
    write_rows(target_path, target_rows)

    pair_rows = []
    n_features = len(super_names)
    for i in range(n_features):
        for j in range(i + 1, n_features):
            pair_rows.append(
                {
                    "feature_a": super_names[i],
                    "feature_b": super_names[j],
                    "pearson": float(corr[i, j]),
                    "abs_pearson": abs(float(corr[i, j])),
                }
            )
    pair_rows.sort(key=lambda row: row["abs_pearson"], reverse=True)
    top_path = outdir / f"{prefix}_input_pairwise_top_correlations.csv"
    write_rows(top_path, pair_rows[: max(1, int(args.correlation_top_n))])

    heatmap_path = outdir / f"{prefix}_input_correlation_heatmap.png"
    heatmap_written = write_correlation_heatmap(heatmap_path, matrix_names, corr)
    payload = {
        "schema": "RJ_AUAU_OOF_SUPERSTACKER_CORRELATION_DIAGNOSTICS_V1",
        "rows_available": int(selected_mask_values.sum()),
        "rows_sampled": int(idx.size),
        "super_feature_count": int(len(super_names)),
        "matrix_column_count": int(len(matrix_names)),
        "final_model": final_name,
        "matrix_csv": str(matrix_path),
        "target_correlations_csv": str(target_path),
        "top_pairwise_correlations_csv": str(top_path),
        "heatmap_png": str(heatmap_path) if heatmap_written else "",
    }
    write_json(outdir / f"{prefix}_correlation_diagnostics.json", payload)
    return payload


def export_residual_artifact(model: ResidualNN, lower_artifacts: dict, metadata: dict) -> dict:
    return {
        "schema": "RJ_AUAU_STACKED_BDT_MLP_RESIDUAL_SUPERSTACKER_V1",
        "name": "oofResidualSuperNN_pt1535",
        "diagnostic_only": True,
        "uses_runtime_bdt_score": True,
        "uses_runtime_mlp_score": True,
        "feature_names": model.feature_names,
        "impute": model.impute,
        "mean": model.mean,
        "scale": model.scale,
        "base_score_name": model.base_score_name,
        "correction_scale": model.correction_scale,
        "layers": model.layers,
        "lower_stackers": lower_artifacts,
        "training": {
            "kind": "out_of_fold_residual_nn",
            "history": model.history,
        },
        "metadata": metadata,
    }


def export_direct_artifact(model, lower_artifacts: dict, metadata: dict, name: str) -> dict:
    return {
        "schema": "RJ_AUAU_OOF_TRISCORE_HEAVYNN_SUPERLEARNER_V1",
        "name": name,
        "diagnostic_only": True,
        "uses_out_of_fold_lower_scores": True,
        "final_model": model.artifact,
        "lower_models": lower_artifacts,
        "training": {
            "kind": "out_of_fold_direct_heavy_nn",
            "history": model.history,
        },
        "metadata": metadata,
    }


def overfit_report(metrics: dict, fold_rows: list[dict], final_name: str) -> dict:
    test = metrics[final_name]["test"]
    train = metrics[final_name]["train_oof"]
    val = metrics[final_name]["val_oof"]
    fold_aucs = [row["oof_auc"] for row in fold_rows if math.isfinite(row["oof_auc"])]
    train_test_gap = train["auc"] - test["auc"] if math.isfinite(train["auc"]) and math.isfinite(test["auc"]) else math.nan
    val_test_gap = val["auc"] - test["auc"] if math.isfinite(val["auc"]) and math.isfinite(test["auc"]) else math.nan
    return {
        "train_oof_minus_test_auc": train_test_gap,
        "val_oof_minus_test_auc": val_test_gap,
        "fold_oof_auc_mean": float(np.mean(fold_aucs)) if fold_aucs else math.nan,
        "fold_oof_auc_std": float(np.std(fold_aucs)) if fold_aucs else math.nan,
        "flags": {
            "large_train_test_auc_gap": bool(math.isfinite(train_test_gap) and train_test_gap > 0.025),
            "large_val_test_auc_gap": bool(math.isfinite(val_test_gap) and val_test_gap > 0.020),
            "large_fold_auc_spread": bool(fold_aucs and float(np.std(fold_aucs)) > 0.020),
            "low_test_finite_fraction": bool(math.isfinite(test["finite_fraction"]) and test["finite_fraction"] < 0.999),
            "high_abs_isolation_correlation": bool(math.isfinite(test["score_vs_eiso_corr"]) and abs(test["score_vs_eiso_corr"]) > 0.75),
        },
    }


def main() -> None:
    args = parse_args()
    if args.pt_min >= args.pt_max:
        raise SystemExit("--pt-min must be less than --pt-max")
    args.outdir.mkdir(parents=True, exist_ok=True)
    print("[superStack] RECOILJETS_AUAU_OOF_RESIDUAL_SUPERSTACKER_V1", flush=True)
    frame, cache_info = load_frame(args)
    stack.add_interactions(frame)
    y = np.asarray(frame["is_signal"], dtype="int32")
    select = selected_mask(frame, args)
    split = stable_split(frame, y, args.random_seed, args.train_fraction, args.val_fraction)
    train_mask = select & split["train"]
    val_mask = select & split["val"]
    test_mask = select & split["test"]
    trainval_mask = train_mask | val_mask
    print(
        f"[superStack] rows={len(y)} selected={int(select.sum())} "
        f"train={int(train_mask.sum())} val={int(val_mask.sum())} test={int(test_mask.sum())}",
        flush=True,
    )
    lower_names, lower_oof, lower_scores, lower_artifacts, fold_rows = train_lower_models_oof(frame, y, trainval_mask, test_mask, args)
    algorithms = [tok.strip() for tok in args.lower_algorithms.split(",") if tok.strip()]
    if args.base_score.startswith("lower_"):
        args.base_score = args.base_score[len("lower_") :]
    if args.base_score not in lower_scores:
        raise SystemExit(f"--base-score must be one of {sorted(lower_scores)}")
    super_names, super_x = make_super_features(frame, lower_scores, algorithms, args.super_feature_mode, args)
    model_stem = args.model_name or ("oofTriScoreHeavyNN_pt1535" if args.final_mode == "direct" else "oofResidualSuperNN_pt1535")
    final_name = "direct_supernn" if args.final_mode == "direct" else "residual_supernn"
    final_history = []
    final_artifact = None
    if args.final_mode == "direct":
        fitted_final = stack.fit_one(
            model_stem,
            "nn",
            super_names,
            super_x,
            y,
            train_mask,
            final_nn_args(args),
            args.random_seed + 9001,
            val_mask=val_mask,
        )
        if fitted_final is None:
            raise SystemExit("Failed final direct heavy NN superlearner")
        final_score_all = fitted_final.predict(super_x)
        final_history = fitted_final.history
    else:
        residual_model = train_residual_nn(super_x, y, lower_scores[args.base_score], train_mask, val_mask, args, super_names)
        final_score_all, delta_raw = predict_residual_nn(residual_model, super_x, lower_scores[args.base_score])
        final_history = residual_model.history

    score_map = {
        "bdt_input": np.asarray(frame["bdt_score"], dtype="float64"),
        "mlp_input": np.asarray(frame["mlp_score"], dtype="float64"),
    }
    for alg in algorithms:
        score_map[f"lower_{alg}"] = lower_scores[alg]
    score_map[final_name] = final_score_all

    masks = {
        "train_oof": train_mask,
        "val_oof": val_mask,
        "trainval_oof": trainval_mask,
        "test": test_mask,
        "all_selected": select,
    }
    metrics = {}
    rank_rows = []
    highpt = highpt_mask(frame)
    for name, score in score_map.items():
        metrics[name] = {}
        for split_name, mask in masks.items():
            m = split_metrics(frame, y, score, mask, args)
            metrics[name][split_name] = m
            hp = split_metrics(frame, y, score, mask & highpt, args)
            rank_rows.append(
                {
                    "model": name,
                    "split": split_name,
                    "entries": m["entries"],
                    "auc": m["auc"],
                    "wp80_fake": m["wp80_fake"],
                    "ece": m["ece"],
                    "finite_fraction": m["finite_fraction"],
                    "score_vs_eiso_corr": m["score_vs_eiso_corr"],
                    "highpt_auc_20_35": hp["auc"],
                    "highpt_wp80_fake_20_35": hp["wp80_fake"],
                }
            )

    metadata = {
        "cache_info": cache_info,
        "selection": f"{args.pt_min:g} <= cluster_Et < {args.pt_max:g}, 0 <= centrality < 80",
        "folds": args.folds,
        "lower_feature_mode": args.lower_feature_mode,
        "super_feature_mode": args.super_feature_mode,
        "final_mode": args.final_mode,
        "lower_feature_names": lower_names,
        "super_feature_names": super_names,
        "lower_algorithms": algorithms,
        "locked_test_split": {"train_fraction": args.train_fraction, "val_fraction": args.val_fraction, "random_seed": args.random_seed},
        "split_mode": split.get("mode", "unknown"),
        "include_isolation_context": bool(args.include_isolation_context),
    }
    overfit = overfit_report(metrics, fold_rows, final_name)
    if args.final_mode == "direct":
        artifact = export_direct_artifact(fitted_final, lower_artifacts, metadata, model_stem)
        out_prefix = "oof_direct_supernn"
    else:
        artifact = export_residual_artifact(residual_model, lower_artifacts, metadata)
        out_prefix = "oof_residual_supernn"
    correlation_payload = {"status": "DISABLED"}
    if not args.skip_correlation_diagnostics:
        print(f"[superStack] writing correlation diagnostics max_rows={args.correlation_max_rows}", flush=True)
        correlation_payload = write_correlation_diagnostics(
            args.outdir,
            out_prefix,
            frame,
            y,
            super_names,
            super_x,
            score_map,
            final_name,
            select,
            args,
        )
        metadata["correlation_diagnostics"] = correlation_payload
    artifact_path = args.outdir / f"{out_prefix}_artifact.json"
    write_json(artifact_path, artifact)
    write_json(
        args.outdir / f"{out_prefix}_metrics.json",
        {"metrics": metrics, "overfit_report": overfit, "fold_rows": fold_rows, "metadata": metadata, "correlation_diagnostics": correlation_payload},
    )
    write_rank_table(args.outdir / f"{out_prefix}_rank_table.csv", rank_rows)
    write_rank_table(args.outdir / f"{out_prefix}_fold_table.csv", fold_rows)
    with (args.outdir / f"{out_prefix}_training_history.csv").open("w", newline="") as handle:
        fields = list(final_history[0]) if final_history else ["epoch"]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(final_history)
    write_json(
        args.outdir / f"{out_prefix}_manifest.json",
        {
            "artifact": str(artifact_path),
            "rank_table": str(args.outdir / f"{out_prefix}_rank_table.csv"),
            "metrics": str(args.outdir / f"{out_prefix}_metrics.json"),
            "final_model": final_name,
        },
    )
    # Also write stable aliases so older status scripts keep working.
    if out_prefix != "oof_residual_supernn":
        write_json(
            args.outdir / "oof_residual_supernn_metrics.json",
            {"metrics": metrics, "overfit_report": overfit, "fold_rows": fold_rows, "metadata": metadata, "correlation_diagnostics": correlation_payload},
        )
        write_rank_table(args.outdir / "oof_residual_supernn_rank_table.csv", rank_rows)
        write_rank_table(args.outdir / "oof_residual_supernn_fold_table.csv", fold_rows)
        with (args.outdir / "oof_residual_supernn_training_history.csv").open("w", newline="") as handle:
            fields = list(final_history[0]) if final_history else ["epoch"]
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(final_history)
        write_json(
            args.outdir / "oof_residual_supernn_manifest.json",
            {
                "artifact": str(artifact_path),
                "rank_table": str(args.outdir / f"{out_prefix}_rank_table.csv"),
                "metrics": str(args.outdir / f"{out_prefix}_metrics.json"),
                "final_model": final_name,
            },
        )
    test = metrics[final_name]["test"]
    final_test_hp = next(
        (row for row in rank_rows if row["model"] == final_name and row["split"] == "test"),
        {},
    )
    print(
        "[superStack] DONE "
        f"test_auc={test['auc']:.6f} test_wp80_fake={test['wp80_fake']} "
        f"highpt_auc={final_test_hp.get('highpt_auc_20_35', math.nan)} "
        f"overfit_flags={overfit['flags']}",
        flush=True,
    )
    print(f"[superStack] DONE_SUPERSTACK_OUTDIR={args.outdir}", flush=True)
    print(f"[superStack] DONE_SUPERSTACK_RANK_TABLE={args.outdir / f'{out_prefix}_rank_table.csv'}", flush=True)


if __name__ == "__main__":
    main()
