#!/usr/bin/env python3
"""Train an out-of-fold residual NN super-stacker for AuAu photon ID.

This is a diagnostic ceiling test.  It uses BDT/MLP score inputs and lower
stacker outputs, so it is not a clean standalone photon-ID candidate.  The
important guardrail is that the residual NN is trained on out-of-fold lower
stacker predictions and evaluated on a locked test split.
"""

from __future__ import annotations

import argparse
import csv
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
    ap.add_argument("--train-fraction", type=float, default=0.60)
    ap.add_argument("--val-fraction", type=float, default=0.20)
    ap.add_argument("--random-seed", type=int, default=24681357)
    ap.add_argument("--target-signal-efficiency", type=float, default=0.80)
    ap.add_argument("--report-pt-bins", default=DEFAULT_PT_BINS)
    ap.add_argument("--report-cent-bins", default=DEFAULT_CENT_BINS)
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
        mlp_score=args.mlp_score,
        bdt_score=args.bdt_score,
        max_shards=args.max_shards,
        max_rows=args.max_rows,
        random_seed=args.random_seed,
    )
    return stack.load_aligned_caches(load_args)


def selected_mask(frame, args):
    et = np.asarray(frame["cluster_Et"], dtype="float64")
    cent = np.asarray(frame["centrality"], dtype="float64")
    return np.isfinite(et) & np.isfinite(cent) & (et >= args.pt_min) & (et < args.pt_max) & (cent >= 0.0) & (cent < 80.0)


def stack_feature_names(frame) -> list[str]:
    stack.add_interactions(frame)
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
    return [name for name in names if name in frame]


def matrix_from_frame(frame, names: list[str]) -> np.ndarray:
    return np.column_stack([np.asarray(frame[name], dtype="float64") for name in names]).astype("float64")


def make_folds(y, mask, n_folds: int, seed: int) -> list[np.ndarray]:
    if n_folds < 2:
        raise SystemExit("--folds must be >= 2")
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
    names = stack_feature_names(frame)
    x = matrix_from_frame(frame, names)
    algorithms = [tok.strip() for tok in args.lower_algorithms.split(",") if tok.strip()]
    if not algorithms:
        raise SystemExit("--lower-algorithms is empty")
    folds = make_folds(y, trainval_mask, args.folds, args.random_seed + 1001)
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


def make_super_features(frame, lower_scores: dict[str, np.ndarray], algorithms: list[str]) -> tuple[list[str], np.ndarray]:
    cols = []
    names = []
    for name in ["bdt_score", "mlp_score", "mlp_logit", "log_cluster_Et", "centrality_scaled"]:
        if name in frame:
            names.append(name)
            cols.append(np.asarray(frame[name], dtype="float64"))
    for alg in algorithms:
        score_name = f"lower_{alg}"
        score = np.asarray(lower_scores[alg], dtype="float64")
        names.append(score_name)
        cols.append(score)
        names.append(f"{score_name}_logit")
        cols.append(logit_prob(score))
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


def overfit_report(metrics: dict, fold_rows: list[dict]) -> dict:
    test = metrics["residual_supernn"]["test"]
    train = metrics["residual_supernn"]["train_oof"]
    val = metrics["residual_supernn"]["val_oof"]
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
    split = stack.stratified_split(y, args.random_seed, args.train_fraction, args.val_fraction)
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
    super_names, super_x = make_super_features(frame, lower_scores, algorithms)
    residual_model = train_residual_nn(super_x, y, lower_scores[args.base_score], train_mask, val_mask, args, super_names)
    residual_score_all, delta_raw = predict_residual_nn(residual_model, super_x, lower_scores[args.base_score])

    score_map = {
        "bdt_input": np.asarray(frame["bdt_score"], dtype="float64"),
        "mlp_input": np.asarray(frame["mlp_score"], dtype="float64"),
    }
    for alg in algorithms:
        score_map[f"lower_{alg}"] = lower_scores[alg]
    score_map["residual_supernn"] = residual_score_all

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

    overfit = overfit_report(metrics, fold_rows)
    artifact = export_residual_artifact(
        residual_model,
        lower_artifacts,
        {
            "cache_info": cache_info,
            "selection": f"{args.pt_min:g} <= cluster_Et < {args.pt_max:g}, 0 <= centrality < 80",
            "folds": args.folds,
            "lower_feature_names": lower_names,
            "lower_algorithms": algorithms,
            "locked_test_split": {"train_fraction": args.train_fraction, "val_fraction": args.val_fraction, "random_seed": args.random_seed},
        },
    )
    artifact_path = args.outdir / "oof_residual_supernn_artifact.json"
    write_json(artifact_path, artifact)
    write_json(args.outdir / "oof_residual_supernn_metrics.json", {"metrics": metrics, "overfit_report": overfit, "fold_rows": fold_rows})
    write_rank_table(args.outdir / "oof_residual_supernn_rank_table.csv", rank_rows)
    write_rank_table(args.outdir / "oof_residual_supernn_fold_table.csv", fold_rows)
    with (args.outdir / "oof_residual_supernn_training_history.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(residual_model.history[0]))
        writer.writeheader()
        writer.writerows(residual_model.history)
    write_json(args.outdir / "oof_residual_supernn_manifest.json", {"artifact": str(artifact_path), "rank_table": str(args.outdir / "oof_residual_supernn_rank_table.csv"), "metrics": str(args.outdir / "oof_residual_supernn_metrics.json")})
    test = metrics["residual_supernn"]["test"]
    print(
        "[superStack] DONE "
        f"test_auc={test['auc']:.6f} test_wp80_fake={test['wp80_fake']} "
        f"highpt_auc={rank_rows[-2]['highpt_auc_20_35'] if rank_rows else math.nan} "
        f"overfit_flags={overfit['flags']}",
        flush=True,
    )
    print(f"[superStack] DONE_SUPERSTACK_OUTDIR={args.outdir}", flush=True)
    print(f"[superStack] DONE_SUPERSTACK_RANK_TABLE={args.outdir / 'oof_residual_supernn_rank_table.csv'}", flush=True)


if __name__ == "__main__":
    main()
