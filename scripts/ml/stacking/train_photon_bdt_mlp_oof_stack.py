#!/usr/bin/env python3
"""Fresh pp/AuAu BDT+MLP stack campaign driver.

This trains the base BDT and base MLP from the same rows and event-level split,
then trains stackers from honest base scores produced on events excluded from
the corresponding base-model training region.  The supervised training class is
always ``is_signal == 1`` versus ``is_signal == 0``.  Source sample names are
provenance used for QA and the separate PPG12-style Signal-MC-vs-Inclusive-MC
overlay diagnostics.
"""

from __future__ import annotations

import argparse
import csv
import gc
import hashlib
import json
import math
import os
import pickle
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[3]
TRAINING_DIR = REPO_ROOT / "scripts" / "ml" / "training"
STACKING_DIR = REPO_ROOT / "scripts" / "ml" / "stacking"
for _path in (str(TRAINING_DIR), str(STACKING_DIR)):
    if _path not in sys.path:
        sys.path.insert(0, _path)

from train_auau_photon_bdt import (  # noqa: E402
    PPG12_EXACT_WEIGHT_COLUMN,
    PPG12_TIGHT_FEATURES,
    add_derived_features,
    global_sixpack_noiso_features,
    load_or_build_frame,
    prepare_ppg12_exact_global_weights,
    summarize_ppg12_exact_precomputed_weights,
)
from train_auau_photon_mlp import auc_score, threshold_for_signal_efficiency  # noqa: E402
from train_auau_stacked_bdt_mlp_sweep import fit_one, nn_predict_from_artifact  # noqa: E402


SCHEMA = "RJ_FRESH_PP_AUAU_BDT_MLP_STACK_V2"
TRAINING_CLASS_DEFINITION = "supervised class label: signal = is_signal == 1; background = is_signal == 0"
OVERLAY_CLASS_DEFINITION = (
    "Signal MC = Photon+Jet source_sample with is_signal == 1; "
    "Inclusive MC = inclusive-jet source_sample with no truth-background filter"
)

PP_FEATURES_BASE_V3E = list(PPG12_TIGHT_FEATURES)
AUAU_SIGNAL_SOURCES = ["run28_embeddedPhoton12", "run28_embeddedPhoton20"]
AUAU_INCLUSIVE_SOURCES = [
    "run28_embeddedJet12",
    "run28_embeddedJet20",
    "run28_embeddedJet30",
    "run28_embeddedJet40",
]
PP_SIGNAL_SOURCES = ["run28_photonjet5", "run28_photonjet10", "run28_photonjet20"]
PP_TRAIN_INCLUSIVE_SOURCES = ["run28_jet8", "run28_jet12", "run28_jet20", "run28_jet30"]
PP_OVERLAY_INCLUSIVE_SOURCES = ["run28_jet8", "run28_jet12", "run28_jet20", "run28_jet30", "run28_jet40"]
MATRIX_DTYPE = np.float32
SCORE_DTYPE = np.float32


@dataclass(frozen=True)
class EventKeyIndex:
    """Compact event-key representation for row-level partitioning.

    The previous implementation stored one combined source/run/evt Python
    string per row.  For full pp samples that can cost many GB and it stays
    live for the whole training campaign.  This stores one integer event id per
    row and one stable string per unique event.
    """

    row_ids: np.ndarray
    key_strings: tuple[str, ...]

    @property
    def n_events(self) -> int:
        return len(self.key_strings)


@dataclass(frozen=True)
class FeatureContract:
    preset: str
    features: list[str]
    description: str


def json_ready(value):
    if isinstance(value, dict):
        return {str(k): json_ready(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_ready(v) for v in value]
    if isinstance(value, np.ndarray):
        return json_ready(value.tolist())
    if isinstance(value, (np.floating, float)):
        f = float(value)
        return f if math.isfinite(f) else None
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, Path):
        return str(value)
    return value


def write_json(path: Path, payload) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(json_ready(payload), indent=2, sort_keys=True) + "\n")


def parse_csv_list(text: str | None) -> list[str]:
    if not text:
        return []
    return [item.strip() for item in text.split(",") if item.strip()]


def parse_range(text: str | None) -> tuple[float, float] | None:
    if not text:
        return None
    if ":" in text:
        lo, hi = text.split(":", 1)
    elif "," in text:
        lo, hi = text.split(",", 1)
    else:
        raise SystemExit(f"Range must be lo:hi: {text}")
    return float(lo), float(hi)


def parse_edges(text: str) -> list[float]:
    edges = [float(item) for item in text.replace(":", ",").split(",") if item.strip()]
    if len(edges) < 2:
        raise SystemExit(f"Need at least two bin edges: {text}")
    if any(edges[i + 1] <= edges[i] for i in range(len(edges) - 1)):
        raise SystemExit(f"Bin edges must be strictly increasing: {text}")
    return edges


def read_input_paths(inputs: list[str]) -> list[Path]:
    paths: list[Path] = []
    for item in inputs:
        if item.startswith("@"):
            manifest = Path(item[1:])
            if not manifest.is_file():
                raise SystemExit(f"Manifest does not exist: {manifest}")
            for line in manifest.read_text().splitlines():
                line = line.strip()
                if line and not line.startswith("#"):
                    paths.append(Path(line))
        else:
            paths.append(Path(item))
    deduped: list[Path] = []
    seen = set()
    for path in paths:
        key = str(path)
        if key not in seen:
            seen.add(key)
            deduped.append(path)
    if not deduped:
        raise SystemExit("No input ROOT files were provided")
    return deduped


def feature_contract(args: argparse.Namespace) -> FeatureContract:
    if args.feature_preset == "pp_basev3e_noiso":
        return FeatureContract(
            preset=args.feature_preset,
            features=list(PP_FEATURES_BASE_V3E),
            description="PPG12-equivalent pp base-v3E noIso 11-feature set.",
        )
    if args.feature_preset == "auau_global_noiso":
        return FeatureContract(
            preset=args.feature_preset,
            features=list(global_sixpack_noiso_features()),
            description="Current successful AuAu noIso global feature family with centrality.",
        )
    if args.features:
        features = parse_csv_list(args.features)
        return FeatureContract(preset="explicit", features=features, description="Explicit CLI feature list.")
    raise SystemExit("Set --feature-preset or --features")


def stable_unit(seed: int, key: str) -> float:
    digest = hashlib.blake2b(f"{seed}|{key}".encode("utf-8"), digest_size=8).digest()
    return int.from_bytes(digest, "big") / float(1 << 64)


def event_keys(frame) -> EventKeyIndex:
    missing = [col for col in ("source_sample", "run", "evt", "is_signal") if col not in frame.columns]
    if missing:
        raise SystemExit(f"Required split/class column(s) missing: {missing}")
    import pandas as pd

    source_series = frame["source_sample"]
    if isinstance(source_series.dtype, pd.CategoricalDtype):
        source_codes = source_series.cat.codes.to_numpy(dtype="int32", copy=False)
        source_labels = source_series.cat.categories.astype(str)
    else:
        source_codes, source_labels = pd.factorize(source_series.astype(str), sort=True)
    run = pd.to_numeric(frame["run"], errors="raise").to_numpy(dtype="int64", copy=False)
    evt = pd.to_numeric(frame["evt"], errors="raise").to_numpy(dtype="int64", copy=False)
    records = np.empty(
        len(frame),
        dtype=[("source", np.int32), ("run", np.int64), ("evt", np.int64)],
    )
    records["source"] = np.asarray(source_codes, dtype="int32")
    records["run"] = run
    records["evt"] = evt
    unique_records, inverse = np.unique(records, return_inverse=True)
    source_names = [str(value) for value in source_labels.tolist()]
    key_strings = tuple(
        f"{source_names[int(row['source'])]}/{int(row['run'])}/{int(row['evt'])}"
        for row in unique_records
    )
    return EventKeyIndex(row_ids=inverse.astype("int64", copy=False), key_strings=key_strings)


def event_count(keys: EventKeyIndex, mask: np.ndarray) -> int:
    return int(len(np.unique(keys.row_ids[mask])))


def event_key_digest(keys: EventKeyIndex, mask: np.ndarray) -> dict[str, object]:
    unique_ids = np.unique(keys.row_ids[mask])
    unique = [keys.key_strings[int(idx)] for idx in unique_ids]
    digest = hashlib.blake2b("\n".join(unique).encode("utf-8"), digest_size=16).hexdigest()
    return {"events": int(len(unique)), "rows": int(mask.sum()), "digest_blake2b16": digest}


def event_key_overlap_count(keys: EventKeyIndex, left: np.ndarray, right: np.ndarray) -> int:
    left_ids = np.unique(keys.row_ids[left])
    right_ids = np.unique(keys.row_ids[right])
    return int(len(np.intersect1d(left_ids, right_ids, assume_unique=True)))


def source_counts(frame, mask: np.ndarray) -> dict[str, int]:
    if mask.sum() == 0:
        return {}
    vals, counts = np.unique(frame.loc[mask, "source_sample"].astype(str).to_numpy(), return_counts=True)
    return {str(v): int(c) for v, c in zip(vals, counts)}


def label_counts(y: np.ndarray, mask: np.ndarray) -> dict[str, int]:
    return {
        "is_signal_0": int(((y == 0) & mask).sum()),
        "is_signal_1": int(((y == 1) & mask).sum()),
    }


def assign_partitions(frame, seed: int, test_fraction: float, folds: int) -> tuple[np.ndarray, np.ndarray, dict]:
    if not (0.0 < test_fraction < 1.0):
        raise SystemExit("--locked-test-fraction must be in (0, 1)")
    if folds < 2:
        raise SystemExit("--folds must be at least 2")
    keys = event_keys(frame)
    key_to_test = np.zeros(keys.n_events, dtype=bool)
    key_to_fold = np.full(keys.n_events, -1, dtype="int16")
    for idx, key_s in enumerate(keys.key_strings):
        is_test = stable_unit(seed, key_s) < test_fraction
        key_to_test[idx] = bool(is_test)
        key_to_fold[idx] = -1 if is_test else int(min(folds - 1, math.floor(stable_unit(seed + 1000003, key_s) * folds)))
    is_test_arr = key_to_test[keys.row_ids]
    fold_id = key_to_fold[keys.row_ids]
    y = frame["is_signal"].to_numpy(dtype="int32")
    trainval = ~is_test_arr
    qa = {
        "schema": "RJ_EVENT_HASH_PARTITION_QA_V1",
        "event_key": "source_sample/run/evt",
        "seed": int(seed),
        "locked_test_fraction": float(test_fraction),
        "folds": int(folds),
        "n_rows": int(len(frame)),
        "n_events": int(keys.n_events),
        "partitions": {
            "trainval": {
                "rows": int(trainval.sum()),
                "events": event_count(keys, trainval),
                "event_digest": event_key_digest(keys, trainval),
                "label_counts": label_counts(y, trainval),
                "source_counts": source_counts(frame, trainval),
            },
            "locked_test": {
                "rows": int(is_test_arr.sum()),
                "events": event_count(keys, is_test_arr),
                "event_digest": event_key_digest(keys, is_test_arr),
                "label_counts": label_counts(y, is_test_arr),
                "source_counts": source_counts(frame, is_test_arr),
            },
        },
        "folds_detail": [],
        "leakage_check": {"event_key_overlap_trainval_locked_test": 0, "status": "pass"},
    }
    overlap_count = event_key_overlap_count(keys, trainval, is_test_arr)
    qa["leakage_check"]["event_key_overlap_trainval_locked_test"] = int(overlap_count)
    if overlap_count:
        qa["leakage_check"]["status"] = "fail"
        raise SystemExit("Event-key leakage between trainval and locked test")
    for fold in range(folds):
        mask = trainval & (fold_id == fold)
        item = {
            "fold": int(fold),
            "rows": int(mask.sum()),
            "events": event_count(keys, mask),
            "event_digest": event_key_digest(keys, mask),
            "label_counts": label_counts(y, mask),
            "source_counts": source_counts(frame, mask),
        }
        qa["folds_detail"].append(item)
        if item["label_counts"]["is_signal_0"] == 0 or item["label_counts"]["is_signal_1"] == 0:
            raise SystemExit(f"Fold {fold} is missing one training class: {item['label_counts']}")
    for name, mask in (("trainval", trainval), ("locked_test", is_test_arr)):
        counts = label_counts(y, mask)
        if counts["is_signal_0"] == 0 or counts["is_signal_1"] == 0:
            raise SystemExit(f"Partition {name} is missing one training class: {counts}")
    return is_test_arr, fold_id, qa


def np_matrix(frame, features: list[str]) -> np.ndarray:
    return np.column_stack([np.asarray(frame[name], dtype=MATRIX_DTYPE) for name in features]).astype(MATRIX_DTYPE, copy=False)


def compact_campaign_frame(frame):
    """Drop loader-only object columns and downcast numeric columns in-place."""

    import pandas as pd

    drop_cols = [
        col
        for col in ("input_file", "global_event_key", "input_file_index", "input_tree_entry")
        if col in frame.columns
    ]
    if drop_cols:
        frame = frame.drop(columns=drop_cols)
    if "source_sample" in frame.columns:
        frame["source_sample"] = frame["source_sample"].astype("category")
    for col in list(frame.columns):
        if col == "source_sample":
            continue
        dtype = frame[col].dtype
        if pd.api.types.is_float_dtype(dtype):
            frame[col] = pd.to_numeric(frame[col], downcast="float")
        elif pd.api.types.is_integer_dtype(dtype) and col not in ("run", "evt"):
            frame[col] = pd.to_numeric(frame[col], downcast="integer")
    return frame


def fit_imputer(x_train: np.ndarray) -> np.ndarray:
    with np.errstate(all="ignore"):
        med = np.nanmedian(np.where(np.isfinite(x_train), x_train, np.nan), axis=0)
    med = np.where(np.isfinite(med), med, 0.0)
    return med.astype(MATRIX_DTYPE, copy=False)


def apply_imputer(x: np.ndarray, impute: np.ndarray) -> np.ndarray:
    return np.where(np.isfinite(x), x, impute).astype(MATRIX_DTYPE, copy=False)


def predict_in_chunks(predict, x: np.ndarray, chunk_rows: int) -> np.ndarray:
    if chunk_rows <= 0 or len(x) <= chunk_rows:
        return np.asarray(predict(x), dtype=SCORE_DTYPE)
    out = np.empty(len(x), dtype=SCORE_DTYPE)
    for start in range(0, len(x), chunk_rows):
        stop = min(start + chunk_rows, len(x))
        out[start:stop] = np.asarray(predict(x[start:stop]), dtype=SCORE_DTYPE)
    return out


def sigmoid(values: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-np.clip(values, -40.0, 40.0)))


def fit_numpy_logistic(x_train: np.ndarray, y_train: np.ndarray, weights: np.ndarray, seed: int):
    rng = np.random.default_rng(seed)
    mean = x_train.mean(axis=0)
    scale = x_train.std(axis=0)
    scale = np.where(scale > 1.0e-12, scale, 1.0)
    xz = (x_train - mean) / scale
    coef = rng.normal(0.0, 0.01, size=xz.shape[1])
    intercept = 0.0
    w = np.asarray(weights, dtype="float64")
    w = w / max(float(w.mean()), 1.0e-12)
    lr = 0.08
    l2 = 1.0e-3
    for _ in range(260):
        pred = sigmoid(xz @ coef + intercept)
        err = (pred - y_train) * w
        grad = (xz.T @ err) / max(1, len(y_train)) + l2 * coef
        bias_grad = float(err.mean())
        coef -= lr * grad
        intercept -= lr * bias_grad
    return {"mean": mean, "scale": scale, "coef": coef, "intercept": float(intercept)}


def fit_bdt(name: str, x: np.ndarray, y: np.ndarray, train_mask: np.ndarray, weights: np.ndarray, args, outdir: Path, seed: int):
    if train_mask.sum() < 20 or len(np.unique(y[train_mask])) < 2:
        raise SystemExit(f"Cannot train BDT {name}: insufficient rows/classes")
    model_dir = outdir / "base_models"
    model_dir.mkdir(parents=True, exist_ok=True)
    impute = fit_imputer(x[train_mask])
    x_train = apply_imputer(x[train_mask], impute)
    y_train = y[train_mask]
    w_train = weights[train_mask]
    backend = "xgboost"
    history: dict[str, object] = {}
    try:
        import xgboost as xgb

        clf = xgb.XGBClassifier(
            objective="binary:logistic",
            eval_metric="auc",
            n_estimators=args.bdt_estimators,
            max_depth=args.bdt_max_depth,
            learning_rate=args.bdt_learning_rate,
            subsample=args.bdt_subsample,
            colsample_bytree=args.bdt_colsample_bytree,
            tree_method=args.bdt_tree_method,
            reg_alpha=args.bdt_reg_alpha,
            reg_lambda=args.bdt_reg_lambda,
            grow_policy=args.bdt_grow_policy,
            max_bin=args.bdt_max_bin,
            n_jobs=args.n_jobs,
            random_state=seed,
        )
        clf.fit(x_train, y_train, sample_weight=w_train, verbose=False)
        model_path = model_dir / f"{name}.xgb.json"
        clf.save_model(str(model_path))
    except Exception as exc:  # noqa: BLE001
        history["xgboost_fallback_reason"] = str(exc)
        try:
            backend = "sklearn_gradient_boosting"
            from sklearn.ensemble import GradientBoostingClassifier

            clf = GradientBoostingClassifier(
                n_estimators=max(20, min(args.bdt_estimators, 160)),
                learning_rate=args.bdt_learning_rate,
                max_depth=min(args.bdt_max_depth, 4),
                subsample=min(1.0, args.bdt_subsample),
                random_state=seed,
            )
            clf.fit(x_train, y_train, sample_weight=w_train)
            model_path = model_dir / f"{name}.sklearn.pkl"
            with model_path.open("wb") as handle:
                pickle.dump(clf, handle)
        except Exception as sklearn_exc:  # noqa: BLE001
            if not args.self_test and not args.allow_bdt_linear_fallback:
                raise SystemExit(
                    "Neither xgboost nor sklearn is available for base-BDT training. "
                    "Install/use the SDCC ML env or pass --allow-bdt-linear-fallback only for a diagnostic smoke."
                ) from sklearn_exc
            backend = "numpy_logistic_fallback_for_self_test"
            history["sklearn_fallback_reason"] = str(sklearn_exc)
            clf = fit_numpy_logistic(x_train, y_train, w_train, seed)
            model_path = model_dir / f"{name}.numpy_logistic.json"
            write_json(model_path, {key: value.tolist() if isinstance(value, np.ndarray) else value for key, value in clf.items()})
    artifact = {
        "name": name,
        "kind": "bdt",
        "backend": backend,
        "model_path": str(model_path),
        "feature_impute": impute.tolist(),
        "n_train_rows": int(train_mask.sum()),
        "train_label_counts": label_counts(y, train_mask),
        **history,
    }
    write_json(model_dir / f"{name}.metadata.json", artifact)

    def predict(x_pred: np.ndarray) -> np.ndarray:
        if backend == "numpy_logistic_fallback_for_self_test":
            xp = apply_imputer(x_pred, impute)
            xz = (xp - clf["mean"]) / clf["scale"]
            return sigmoid(xz @ clf["coef"] + clf["intercept"]).astype("float64")
        return np.asarray(clf.predict_proba(apply_imputer(x_pred, impute))[:, 1], dtype="float64")

    return artifact, predict


def load_existing_bdt(name: str, outdir: Path, expected_train_mask: np.ndarray | None = None, y: np.ndarray | None = None):
    model_dir = outdir / "base_models"
    metadata_path = model_dir / f"{name}.metadata.json"
    if not metadata_path.is_file():
        return None
    artifact = json.loads(metadata_path.read_text())
    model_path = Path(str(artifact.get("model_path", "")))
    if not model_path.is_file():
        return None
    if expected_train_mask is not None:
        expected_rows = int(expected_train_mask.sum())
        if int(artifact.get("n_train_rows", -1)) != expected_rows:
            raise SystemExit(f"Existing BDT {name} has n_train_rows={artifact.get('n_train_rows')} but expected {expected_rows}")
        if y is not None and artifact.get("train_label_counts") != label_counts(y, expected_train_mask):
            raise SystemExit(f"Existing BDT {name} label counts do not match the current partition")
    backend = str(artifact.get("backend", ""))
    impute = np.asarray(artifact.get("feature_impute", []), dtype=MATRIX_DTYPE)
    if backend == "xgboost" or model_path.suffixes[-2:] == [".xgb", ".json"]:
        import xgboost as xgb

        clf = xgb.XGBClassifier()
        clf.load_model(str(model_path))

        def predict(x_pred: np.ndarray) -> np.ndarray:
            return np.asarray(clf.predict_proba(apply_imputer(x_pred, impute))[:, 1], dtype="float64")

        return artifact, predict
    if backend == "sklearn_gradient_boosting" or model_path.suffix == ".pkl":
        with model_path.open("rb") as handle:
            clf = pickle.load(handle)

        def predict(x_pred: np.ndarray) -> np.ndarray:
            return np.asarray(clf.predict_proba(apply_imputer(x_pred, impute))[:, 1], dtype="float64")

        return artifact, predict
    if backend == "numpy_logistic_fallback_for_self_test":
        payload = json.loads(model_path.read_text())
        coef = np.asarray(payload["coef"], dtype="float64")
        mean = np.asarray(payload["mean"], dtype="float64")
        scale = np.asarray(payload["scale"], dtype="float64")
        intercept = float(payload["intercept"])

        def predict(x_pred: np.ndarray) -> np.ndarray:
            xp = apply_imputer(x_pred, impute)
            xz = (xp - mean) / scale
            return sigmoid(xz @ coef + intercept).astype("float64")

        return artifact, predict
    raise SystemExit(f"Unsupported existing BDT backend for {name}: {backend}")


def stack_args_from(base_args: argparse.Namespace) -> argparse.Namespace:
    return argparse.Namespace(
        l2=base_args.stack_l2,
        linear_backend=base_args.stack_linear_backend,
        max_linear_steps=base_args.stack_max_linear_steps,
        gbm_estimators=base_args.stack_gbm_estimators,
        gbm_learning_rate=base_args.stack_gbm_learning_rate,
        gbm_max_depth=base_args.stack_gbm_max_depth,
        gbm_max_leaf_nodes=base_args.stack_gbm_max_leaf_nodes,
        nn_hidden=base_args.stack_mlp_hidden,
        nn_epochs=base_args.stack_mlp_epochs,
        nn_patience=base_args.stack_mlp_patience,
        nn_batch_size=base_args.stack_mlp_batch_size,
        nn_learning_rate=base_args.stack_mlp_learning_rate,
        nn_l2=base_args.stack_mlp_l2,
    )


def train_mlp_model(name: str, features: list[str], x: np.ndarray, y: np.ndarray, train_mask: np.ndarray, val_mask: np.ndarray, weights: np.ndarray, args, outdir: Path, seed: int):
    fitted = fit_one(
        name,
        "nn",
        features,
        x,
        y,
        train_mask,
        stack_args_from(args),
        seed,
        val_mask=val_mask,
        sample_weights=weights,
    )
    if fitted is None:
        raise SystemExit(f"Cannot train MLP {name}: fit_one returned None")
    model_dir = outdir / "base_models"
    model_dir.mkdir(parents=True, exist_ok=True)
    artifact_path = model_dir / f"{name}.artifact.json"
    pickle_path = model_dir / f"{name}.pkl"
    write_json(artifact_path, fitted.artifact)
    with pickle_path.open("wb") as handle:
        pickle.dump({"algorithm": fitted.algorithm, "feature_names": fitted.feature_names, "model": fitted.model_object}, handle)
    metadata = {
        "name": name,
        "kind": "mlp",
        "algorithm": fitted.algorithm,
        "artifact_path": str(artifact_path),
        "pickle_path": str(pickle_path),
        "n_train_rows": int(train_mask.sum()),
        "n_val_rows": int(val_mask.sum()),
        "train_label_counts": label_counts(y, train_mask),
        "val_label_counts": label_counts(y, val_mask),
        "history": fitted.history,
    }
    write_json(model_dir / f"{name}.metadata.json", metadata)
    return metadata, fitted.predict


def load_existing_mlp(name: str, features: list[str], outdir: Path, expected_train_mask: np.ndarray | None = None, y: np.ndarray | None = None):
    model_dir = outdir / "base_models"
    metadata_path = model_dir / f"{name}.metadata.json"
    if not metadata_path.is_file():
        return None
    metadata = json.loads(metadata_path.read_text())
    artifact_path = Path(str(metadata.get("artifact_path", model_dir / f"{name}.artifact.json")))
    if not artifact_path.is_file():
        return None
    artifact = json.loads(artifact_path.read_text())
    if list(artifact.get("feature_names", [])) != list(features):
        raise SystemExit(f"Existing MLP {name} feature contract does not match the current run")
    if expected_train_mask is not None:
        expected_rows = int(expected_train_mask.sum())
        if int(metadata.get("n_train_rows", -1)) != expected_rows:
            raise SystemExit(f"Existing MLP {name} has n_train_rows={metadata.get('n_train_rows')} but expected {expected_rows}")
        if y is not None and metadata.get("train_label_counts") != label_counts(y, expected_train_mask):
            raise SystemExit(f"Existing MLP {name} label counts do not match the current partition")

    def predict(x_pred: np.ndarray) -> np.ndarray:
        return np.asarray(nn_predict_from_artifact(artifact, x_pred), dtype="float64")

    return metadata, predict


def fit_or_load_bdt(name: str, x: np.ndarray, y: np.ndarray, train_mask: np.ndarray, weights: np.ndarray, args, outdir: Path, seed: int):
    if args.reuse_existing_base_models:
        loaded = load_existing_bdt(name, outdir, train_mask, y)
        if loaded is not None:
            return loaded
    return fit_bdt(name, x, y, train_mask, weights, args, outdir, seed)


def train_or_load_mlp_model(name: str, features: list[str], x: np.ndarray, y: np.ndarray, train_mask: np.ndarray, val_mask: np.ndarray, weights: np.ndarray, args, outdir: Path, seed: int):
    if args.reuse_existing_base_models:
        loaded = load_existing_mlp(name, features, outdir, train_mask, y)
        if loaded is not None:
            return loaded
    return train_mlp_model(name, features, x, y, train_mask, val_mask, weights, args, outdir, seed)


def weighted_threshold_report(y: np.ndarray, score: np.ndarray, weights: np.ndarray, target: float) -> dict:
    base = threshold_for_signal_efficiency(y, score, target) or {}
    mask = np.isfinite(score) & np.isin(y, [0, 1]) & np.isfinite(weights) & (weights > 0)
    sig = mask & (y == 1)
    bkg = mask & (y == 0)
    if sig.sum() == 0 or bkg.sum() == 0:
        return base
    order = np.argsort(score[sig])
    s_scores = score[sig][order]
    s_weights = weights[sig][order]
    cdf = np.cumsum(s_weights) / float(s_weights.sum())
    idx = int(np.searchsorted(cdf, max(0.0, min(1.0, 1.0 - target)), side="left"))
    idx = min(max(idx, 0), len(s_scores) - 1)
    threshold = float(s_scores[idx])
    base.update(
        {
            "weighted_threshold": threshold,
            "weighted_signal_efficiency": float(weights[sig & (score > threshold)].sum() / weights[sig].sum()),
            "weighted_background_fake_rate": float(weights[bkg & (score > threshold)].sum() / weights[bkg].sum()),
        }
    )
    return base


def metric_row(domain: str, region: str, model: str, y: np.ndarray, score: np.ndarray, weights: np.ndarray, mask: np.ndarray, target_eff: float) -> dict:
    valid = mask & np.isfinite(score) & np.isin(y, [0, 1])
    counts = label_counts(y, valid)
    row = {
        "domain": domain,
        "region": region,
        "model": model,
        "rows": int(valid.sum()),
        "is_signal_1": counts["is_signal_1"],
        "is_signal_0": counts["is_signal_0"],
        "weighted_auc": auc_score(y[valid], score[valid], weights[valid]) if valid.any() else math.nan,
        "unweighted_auc": auc_score(y[valid], score[valid], None) if valid.any() else math.nan,
    }
    wp = weighted_threshold_report(y[valid], score[valid], weights[valid], target_eff) if valid.any() else {}
    row.update({f"wp{int(target_eff*100):02d}_{k}": v for k, v in wp.items()})
    return row


def source_mask(frame, samples: list[str]) -> np.ndarray:
    if not samples:
        return np.ones(len(frame), dtype=bool)
    src = [str(value).lower() for value in frame["source_sample"].to_numpy()]
    sample_l = [sample.lower() for sample in samples]
    return np.asarray([any(value == sample or sample in value for sample in sample_l) for value in src], dtype=bool)


def histogram_rows(domain: str, frame, model_scores: dict[str, np.ndarray], weights: np.ndarray, test_mask: np.ndarray, signal_sources: list[str], inclusive_sources: list[str], bins: int = 60) -> list[dict]:
    rows: list[dict] = []
    y = frame["is_signal"].to_numpy(dtype="int32")
    signal_mask = test_mask & source_mask(frame, signal_sources) & (y == 1)
    inclusive_mask = test_mask & source_mask(frame, inclusive_sources)
    edges = np.linspace(0.0, 1.0, bins + 1)
    width = float(edges[1] - edges[0])
    for model, score in model_scores.items():
        for label, mask in (("Signal MC (truth prompt)", signal_mask), ("Inclusive MC (source inclusive jet)", inclusive_mask)):
            valid = mask & np.isfinite(score) & np.isfinite(weights) & (weights > 0)
            hist, _ = np.histogram(score[valid], bins=edges, weights=weights[valid])
            norm = float(hist.sum())
            density = hist / (norm * width) if norm > 0 else np.zeros_like(hist, dtype="float64")
            for i in range(len(hist)):
                rows.append(
                    {
                        "domain": domain,
                        "model": model,
                        "overlay_class": label,
                        "bin_low": float(edges[i]),
                        "bin_high": float(edges[i + 1]),
                        "weighted_count": float(hist[i]),
                        "unit_area_density": float(density[i]),
                        "source_rows": int(valid.sum()),
                    }
                )
    return rows


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: json_ready(row.get(key, "")) for key in fields})


def stratified_metric_rows(domain: str, frame, model_scores: dict[str, np.ndarray], weights: np.ndarray, test_mask: np.ndarray, et_edges: list[float], cent_edges: list[float] | None, target_eff: float) -> list[dict]:
    rows: list[dict] = []
    y = frame["is_signal"].to_numpy(dtype="int32")
    et = frame["cluster_Et"].to_numpy(dtype="float64")
    for lo, hi in zip(et_edges[:-1], et_edges[1:]):
        mask = test_mask & np.isfinite(et) & (et >= lo) & (et < hi)
        for model, score in model_scores.items():
            row = metric_row(domain, f"Et_{lo:g}_{hi:g}", model, y, score, weights, mask, target_eff)
            row.update({"stratification": "Et", "bin_lo": lo, "bin_hi": hi})
            rows.append(row)
    if cent_edges is not None and "centrality" in frame.columns:
        cent = frame["centrality"].to_numpy(dtype="float64")
        for lo, hi in zip(cent_edges[:-1], cent_edges[1:]):
            mask = test_mask & np.isfinite(cent) & (cent >= lo) & (cent < hi)
            for model, score in model_scores.items():
                row = metric_row(domain, f"cent_{lo:g}_{hi:g}", model, y, score, weights, mask, target_eff)
                row.update({"stratification": "centrality", "bin_lo": lo, "bin_hi": hi})
                rows.append(row)
        for elo, ehi in zip(et_edges[:-1], et_edges[1:]):
            for clo, chi in zip(cent_edges[:-1], cent_edges[1:]):
                mask = test_mask & np.isfinite(et) & (et >= elo) & (et < ehi) & np.isfinite(cent) & (cent >= clo) & (cent < chi)
                for model, score in model_scores.items():
                    row = metric_row(domain, f"Et_{elo:g}_{ehi:g}_cent_{clo:g}_{chi:g}", model, y, score, weights, mask, target_eff)
                    row.update(
                        {
                            "stratification": "Et_x_centrality",
                            "et_lo": elo,
                            "et_hi": ehi,
                            "cent_lo": clo,
                            "cent_hi": chi,
                        }
                    )
                    rows.append(row)
    return rows


def score_correlations(frame, scores: dict[str, np.ndarray], test_mask: np.ndarray, signal_sources: list[str], inclusive_sources: list[str]) -> dict:
    y = frame["is_signal"].to_numpy(dtype="int32")
    masks = {
        "all_locked_test": test_mask,
        "truth_signal_locked_test": test_mask & (y == 1),
        "truth_background_locked_test": test_mask & (y == 0),
        "overlay_signal_mc_locked_test": test_mask & source_mask(frame, signal_sources) & (y == 1),
        "overlay_inclusive_mc_locked_test": test_mask & source_mask(frame, inclusive_sources),
    }
    out: dict[str, object] = {"schema": "RJ_SCORE_CORRELATION_QA_V1", "pairs": {}}
    names = list(scores)
    for i, left in enumerate(names):
        for right in names[i + 1 :]:
            pair_key = f"{left}__vs__{right}"
            out["pairs"][pair_key] = {}
            for region, mask in masks.items():
                valid = mask & np.isfinite(scores[left]) & np.isfinite(scores[right])
                if valid.sum() < 3:
                    corr = math.nan
                else:
                    corr = float(np.corrcoef(scores[left][valid], scores[right][valid])[0, 1])
                out["pairs"][pair_key][region] = {"rows": int(valid.sum()), "pearson": corr}
    return out


def self_test_frame(args: argparse.Namespace, features: list[str]):
    import pandas as pd

    rng = np.random.default_rng(args.random_seed)
    n = args.self_test_rows
    if args.domain == "pp":
        signal_sources = PP_SIGNAL_SOURCES
        inclusive_sources = PP_TRAIN_INCLUSIVE_SOURCES
        cent = np.full(n, -1.0)
        et_lo, et_hi = 5.0, 35.0
    else:
        signal_sources = AUAU_SIGNAL_SOURCES
        inclusive_sources = AUAU_INCLUSIVE_SOURCES
        cent = rng.uniform(0.0, 80.0, size=n)
        et_lo, et_hi = 15.0, 35.0
    y = rng.binomial(1, 0.52, size=n).astype("int32")
    source = np.empty(n, dtype=object)
    for idx in range(n):
        source[idx] = rng.choice(signal_sources if y[idx] else inclusive_sources)
    et = rng.uniform(et_lo, et_hi, size=n) + y * rng.normal(1.2, 0.4, size=n)
    data = {
        "source_sample": source,
        "run": rng.integers(1000, 9999, size=n),
        "evt": np.arange(n) // 3,
        "is_signal": y,
        "cluster_Et": et,
        "cluster_Eta": rng.uniform(-0.7, 0.7, size=n),
        "centrality": cent,
    }
    for j, feature in enumerate(features):
        if feature in data:
            continue
        data[feature] = rng.normal(0.0, 1.0, size=n) + y * (0.15 + 0.02 * (j % 5))
    frame = pd.DataFrame(data)
    return add_derived_features(frame)


def load_frame_for_campaign(args: argparse.Namespace, features: list[str]):
    if args.self_test:
        return self_test_frame(args, features), [], False, []
    paths = read_input_paths(args.input)
    # The existing BDT loader derives source_sample from the ROOT path when the
    # tree does not carry that provenance branch. The split check below still
    # requires source_sample/run/evt/is_signal after loading.
    required = sorted(set(features + ["run", "evt", args.label_branch, "cluster_Et", "cluster_Eta"]))
    if args.domain == "auau" or "centrality" in features:
        required.append("centrality")
    optional = ["event_weight", "weight"]
    frame, optional_seen, from_cache = load_or_build_frame(
        paths,
        args.tree,
        required,
        optional,
        args.label_branch,
        None,
        args.cache_file,
        cache_only=False,
        skip_missing_tree=args.skip_missing_tree,
        max_load_rows_per_class=args.max_load_rows_per_class,
        max_load_rows=args.max_load_rows,
        load_sample_seed=args.load_sample_seed,
    )
    frame = add_derived_features(frame)
    return frame, optional_seen, from_cache, [str(path) for path in paths]


def load_overlay_frame(args: argparse.Namespace, features: list[str]):
    if not args.overlay_input:
        return None, []
    overlay_args = argparse.Namespace(**vars(args))
    overlay_args.input = args.overlay_input
    overlay_args.cache_file = args.overlay_cache_file
    frame, _, _, paths = load_frame_for_campaign(overlay_args, features)
    return filter_rows(frame, args), paths


def filter_rows(frame, args: argparse.Namespace):
    if "centrality" not in frame.columns:
        frame = frame.copy()
        frame["centrality"] = -1.0
    y = frame[args.label_branch].to_numpy(dtype="int32")
    mask = np.isin(y, [0, 1])
    pt_range = parse_range(args.pt_range)
    if pt_range is not None:
        et = frame["cluster_Et"].to_numpy(dtype="float64")
        mask &= np.isfinite(et) & (et >= pt_range[0]) & (et < pt_range[1])
    cent_range = parse_range(args.centrality_range)
    if cent_range is not None:
        cent = frame["centrality"].to_numpy(dtype="float64")
        mask &= np.isfinite(cent) & (cent >= cent_range[0]) & (cent < cent_range[1])
    filtered = frame.loc[mask].reset_index(drop=True)
    if len(filtered) == 0:
        raise SystemExit("Selection left zero rows")
    return filtered


def prepare_weights(frame, args: argparse.Namespace) -> tuple[np.ndarray, dict]:
    if args.weight_mode == "none":
        return np.ones(len(frame), dtype="float64"), {"weight_mode": "none"}
    proxy = argparse.Namespace(
        outdir=args.outdir,
        ppg12_exact_closure_dir=args.ppg12_exact_closure_dir,
        ppg12_exact_expected_samples=args.ppg12_exact_expected_samples,
    )
    if PPG12_EXACT_WEIGHT_COLUMN in frame.columns and np.isfinite(frame[PPG12_EXACT_WEIGHT_COLUMN].to_numpy(dtype="float64")).all():
        report = summarize_ppg12_exact_precomputed_weights(frame, args.label_branch, proxy)
        weights = frame[PPG12_EXACT_WEIGHT_COLUMN].to_numpy(dtype="float64")
        return weights, report
    weighted_frame, report = prepare_ppg12_exact_global_weights(frame, args.label_branch, proxy)
    frame[PPG12_EXACT_WEIGHT_COLUMN] = weighted_frame[PPG12_EXACT_WEIGHT_COLUMN].to_numpy(dtype="float64")
    return frame[PPG12_EXACT_WEIGHT_COLUMN].to_numpy(dtype="float64"), report


def train_and_score(args: argparse.Namespace) -> dict:
    start = time.time()
    args.outdir.mkdir(parents=True, exist_ok=True)
    contract = feature_contract(args)
    frame, optional_seen, from_cache, input_paths = load_frame_for_campaign(args, contract.features)
    frame = filter_rows(frame, args)
    frame = compact_campaign_frame(frame)
    missing_features = [feature for feature in contract.features if feature not in frame.columns]
    if missing_features:
        raise SystemExit(f"Selected feature(s) missing after derived-feature build: {missing_features}")
    y = frame[args.label_branch].to_numpy(dtype="int32")
    weights, weight_report = prepare_weights(frame, args)
    test_mask, fold_id, partition_qa = assign_partitions(frame, args.random_seed, args.locked_test_fraction, args.folds)
    trainval_mask = ~test_mask
    write_json(args.outdir / "partition_qa.json", partition_qa)
    feature_metadata = {
        "preset": contract.preset,
        "description": contract.description,
        "features": contract.features,
        "n_features": len(contract.features),
    }
    write_json(args.outdir / "feature_contract.json", feature_metadata)
    if args.preflight_only:
        manifest = {
            "schema": SCHEMA,
            "status": "PREFLIGHT_ONLY",
            "domain": args.domain,
            "input_paths": input_paths,
            "feature_contract": feature_metadata,
            "partition_qa": partition_qa,
            "weight_report": weight_report,
            "training_class_definition": TRAINING_CLASS_DEFINITION,
            "overlay_class_definition": OVERLAY_CLASS_DEFINITION,
            "stack_training_mode": args.stack_training_mode,
            "oof_fold_count": int(args.folds),
            "stack_validation_fold": int(args.stack_validation_fold),
            "simple_stack_fold": int(args.simple_stack_fold if args.simple_stack_fold is not None else args.stack_validation_fold),
        }
        write_json(args.outdir / "campaign_manifest.json", manifest)
        return manifest

    x = np_matrix(frame, contract.features)
    keys = event_keys(frame)
    oof_bdt = np.full(len(frame), np.nan, dtype=SCORE_DTYPE)
    oof_mlp = np.full(len(frame), np.nan, dtype=SCORE_DTYPE)
    stack_val_fold = int(args.stack_validation_fold)
    simple_stack_fold = int(args.simple_stack_fold if args.simple_stack_fold is not None else args.stack_validation_fold)
    if simple_stack_fold < 0 or simple_stack_fold >= args.folds:
        raise SystemExit("--simple-stack-fold must be in [0, folds)")

    base_artifacts: dict[str, object] = {
        "stack_training_mode": args.stack_training_mode,
        "folds": [],
        "simple_holdout": None,
    }
    if args.stack_training_mode == "oof5":
        for fold in range(args.folds):
            pred_mask = trainval_mask & (fold_id == fold)
            fit_mask = trainval_mask & (fold_id != fold)
            bdt_artifact, bdt_predict = fit_or_load_bdt(f"fold{fold}_bdt", x, y, fit_mask, weights, args, args.outdir, args.random_seed + fold + 1)
            oof_bdt[pred_mask] = predict_in_chunks(bdt_predict, x[pred_mask], args.predict_chunk_rows)
            del bdt_predict
            gc.collect()
            mlp_artifact, mlp_predict = train_or_load_mlp_model(
                f"fold{fold}_mlp",
                contract.features,
                x,
                y,
                fit_mask,
                pred_mask,
                weights,
                args,
                args.outdir,
                args.random_seed + 100 + fold,
            )
            oof_mlp[pred_mask] = predict_in_chunks(mlp_predict, x[pred_mask], args.predict_chunk_rows)
            del mlp_predict
            gc.collect()
            base_artifacts["folds"].append({"fold": fold, "heldout_rows": int(pred_mask.sum()), "bdt": bdt_artifact, "mlp": mlp_artifact})
        stack_val_mask = trainval_mask & (fold_id == stack_val_fold)
        stack_train_mask = trainval_mask & (fold_id != stack_val_fold)
        stack_score_mask = trainval_mask
        stack_eval_extra_mask = stack_val_mask
        base_score_metric_region = "oof_trainval"
        stack_training_contract = {
            "mode": "oof5",
            "description": "Stackers train on five-fold out-of-fold BDT/MLP base scores over the trainval region.",
            "base_score_region": "all non-test trainval folds, each scored by the fold model that excluded that fold from base training",
            "stack_train_region": f"OOF folds excluding validation fold {stack_val_fold}",
            "stack_validation_region": f"OOF fold {stack_val_fold}",
            "fold_count": int(args.folds),
            "stack_validation_fold": stack_val_fold,
            "event_key_overlap_checks": [
                {
                    "fold": int(fold),
                    "base_train_vs_scored_fold": event_key_overlap_count(keys, trainval_mask & (fold_id != fold), trainval_mask & (fold_id == fold)),
                }
                for fold in range(args.folds)
            ],
        }
    elif args.stack_training_mode == "simple_holdout":
        pred_mask = trainval_mask & (fold_id == simple_stack_fold)
        fit_mask = trainval_mask & (fold_id != simple_stack_fold)
        bdt_artifact, bdt_predict = fit_or_load_bdt(
            f"simple_holdout_fold{simple_stack_fold}_bdt",
            x,
            y,
            fit_mask,
            weights,
            args,
            args.outdir,
            args.random_seed + 501,
        )
        oof_bdt[pred_mask] = predict_in_chunks(bdt_predict, x[pred_mask], args.predict_chunk_rows)
        del bdt_predict
        gc.collect()
        mlp_artifact, mlp_predict = train_or_load_mlp_model(
            f"simple_holdout_fold{simple_stack_fold}_mlp",
            contract.features,
            x,
            y,
            fit_mask,
            pred_mask,
            weights,
            args,
            args.outdir,
            args.random_seed + 601,
        )
        oof_mlp[pred_mask] = predict_in_chunks(mlp_predict, x[pred_mask], args.predict_chunk_rows)
        del mlp_predict
        gc.collect()
        base_artifacts["simple_holdout"] = {
            "heldout_fold": simple_stack_fold,
            "base_train_rows": int(fit_mask.sum()),
            "stack_training_score_rows": int(pred_mask.sum()),
            "bdt": bdt_artifact,
            "mlp": mlp_artifact,
        }
        stack_val_mask = np.zeros(len(frame), dtype=bool)
        stack_train_mask = pred_mask
        stack_score_mask = pred_mask
        stack_eval_extra_mask = pred_mask
        base_score_metric_region = f"simple_holdout_stack_training_fold{simple_stack_fold}"
        stack_training_contract = {
            "mode": "simple_holdout",
            "description": "Stackers train on one held-out trainval block scored by BDT/MLP base models trained on the other trainval blocks.",
            "base_score_region": f"trainval fold {simple_stack_fold} only",
            "base_train_region": f"trainval folds excluding fold {simple_stack_fold}",
            "stack_train_region": f"scored trainval fold {simple_stack_fold}",
            "stack_validation_region": "stacker MLP uses an internal split inside the simple stack-training block when needed for early stopping",
            "fold_count": int(args.folds),
            "simple_stack_fold": simple_stack_fold,
            "event_key_overlap_checks": [
                {
                    "simple_base_train_vs_stack_training_block": event_key_overlap_count(keys, fit_mask, pred_mask),
                    "simple_stack_training_block_vs_locked_test": event_key_overlap_count(keys, pred_mask, test_mask),
                }
            ],
        }
    else:
        raise SystemExit(f"Unsupported --stack-training-mode: {args.stack_training_mode}")

    final_bdt_artifact, final_bdt_predict = fit_or_load_bdt("final_trainval_bdt", x, y, trainval_mask, weights, args, args.outdir, args.random_seed + 999)
    final_mlp_artifact, final_mlp_predict = train_or_load_mlp_model(
        "final_trainval_mlp",
        contract.features,
        x,
        y,
        trainval_mask,
        trainval_mask & (fold_id == stack_val_fold),
        weights,
        args,
        args.outdir,
        args.random_seed + 1999,
    )
    final_bdt_score = np.full(len(frame), np.nan, dtype=SCORE_DTYPE)
    final_mlp_score = np.full(len(frame), np.nan, dtype=SCORE_DTYPE)
    final_bdt_score[test_mask] = predict_in_chunks(final_bdt_predict, x[test_mask], args.predict_chunk_rows)
    final_mlp_score[test_mask] = predict_in_chunks(final_mlp_predict, x[test_mask], args.predict_chunk_rows)
    base_artifacts["final"] = {"bdt": final_bdt_artifact, "mlp": final_mlp_artifact}
    write_json(args.outdir / "base_model_artifacts.json", base_artifacts)

    stack_bdt_input = np.where(test_mask, final_bdt_score, oof_bdt)
    stack_mlp_input = np.where(test_mask, final_mlp_score, oof_mlp)
    stack_columns = {
        "bdt_score": stack_bdt_input,
        "mlp_score": stack_mlp_input,
        "cluster_Et": frame["cluster_Et"].to_numpy(dtype=MATRIX_DTYPE),
    }
    if args.domain == "auau":
        stack_columns["centrality"] = frame["centrality"].to_numpy(dtype=MATRIX_DTYPE)
    stack_specs = {
        "score_only": ["bdt_score", "mlp_score"],
        "score_context": ["bdt_score", "mlp_score", "cluster_Et"] + (["centrality"] if args.domain == "auau" else []),
    }
    stack_artifacts = []
    stack_scores: dict[str, np.ndarray] = {}
    stack_predictors: dict[str, tuple[list[str], object]] = {}
    algorithms = parse_csv_list(args.stack_algorithms)
    sargs = stack_args_from(args)
    for stack_name, stack_features in stack_specs.items():
        sx = np.column_stack([stack_columns[name] for name in stack_features]).astype(MATRIX_DTYPE, copy=False)
        for alg_index, algorithm in enumerate(algorithms):
            model_name = f"{stack_name}_{algorithm}"
            fitted = fit_one(
                model_name,
                "nn" if algorithm == "mlp" else algorithm,
                stack_features,
                sx,
                y,
                stack_train_mask,
                sargs,
                args.random_seed + 3000 + alg_index,
                val_mask=stack_val_mask,
                sample_weights=weights,
            )
            if fitted is None:
                continue
            score = np.full(len(frame), np.nan, dtype=SCORE_DTYPE)
            eval_mask = test_mask | stack_eval_extra_mask | stack_val_mask
            score[eval_mask] = predict_in_chunks(fitted.predict, sx[eval_mask], args.predict_chunk_rows)
            stack_scores[model_name] = score
            stack_predictors[model_name] = (stack_features, fitted.predict)
            artifact_dir = args.outdir / "stack_models"
            artifact_dir.mkdir(parents=True, exist_ok=True)
            artifact_path = artifact_dir / f"{model_name}.artifact.json"
            pickle_path = artifact_dir / f"{model_name}.pkl"
            write_json(artifact_path, fitted.artifact)
            with pickle_path.open("wb") as handle:
                pickle.dump({"algorithm": fitted.algorithm, "feature_names": fitted.feature_names, "model": fitted.model_object}, handle)
            stack_artifacts.append(
                {
                    "model": model_name,
                    "algorithm": algorithm,
                    "fit_algorithm": fitted.algorithm,
                    "stack_input_set": stack_name,
                    "features": stack_features,
                    "stack_training_mode": args.stack_training_mode,
                    "train_region": stack_training_contract["stack_train_region"],
                    "validation_region": stack_training_contract["stack_validation_region"],
                    "artifact_path": str(artifact_path),
                    "pickle_path": str(pickle_path),
                    "history": fitted.history,
                }
            )
    write_json(args.outdir / "stack_model_artifacts.json", {"schema": "RJ_STACK_MODEL_ARTIFACTS_V1", "models": stack_artifacts})

    model_scores = {
        "BDT": final_bdt_score,
        "MLP": final_mlp_score,
        **{name: score for name, score in stack_scores.items()},
    }
    oof_scores = {"BDT": oof_bdt, "MLP": oof_mlp}
    metric_rows = []
    for name, score in model_scores.items():
        metric_rows.append(metric_row(args.domain, "locked_test", name, y, score, weights, test_mask, args.target_signal_efficiency))
    for name, score in oof_scores.items():
        metric_rows.append(metric_row(args.domain, base_score_metric_region, name, y, score, weights, stack_score_mask, args.target_signal_efficiency))
    if stack_val_mask.any():
        for name, score in stack_scores.items():
            metric_rows.append(metric_row(args.domain, f"oof_stack_validation_fold{stack_val_fold}", name, y, score, weights, stack_val_mask, args.target_signal_efficiency))
    elif stack_eval_extra_mask.any():
        for name, score in stack_scores.items():
            metric_rows.append(
                metric_row(
                    args.domain,
                    f"simple_holdout_stack_training_fold{simple_stack_fold}_fit_diagnostic",
                    name,
                    y,
                    score,
                    weights,
                    stack_eval_extra_mask,
                    args.target_signal_efficiency,
                )
            )
    write_csv(args.outdir / "model_metrics.csv", metric_rows)
    write_json(args.outdir / "model_metrics.json", {"schema": "RJ_MODEL_METRICS_V1", "rows": metric_rows})

    et_edges = parse_edges(args.report_et_bins)
    cent_edges = parse_edges(args.report_cent_bins) if args.domain == "auau" and args.report_cent_bins else None
    strat_rows = stratified_metric_rows(args.domain, frame, model_scores, weights, test_mask, et_edges, cent_edges, args.target_signal_efficiency)
    write_csv(args.outdir / "stratified_metrics.csv", strat_rows)

    signal_sources = parse_csv_list(args.overlay_signal_sources)
    inclusive_sources = parse_csv_list(args.overlay_inclusive_sources)
    overlay_frame, overlay_paths = load_overlay_frame(args, contract.features)
    if overlay_frame is not None:
        overlay_frame = compact_campaign_frame(overlay_frame)
        overlay_weights = np.ones(len(overlay_frame), dtype="float64")
        ox = np_matrix(overlay_frame, contract.features)
        overlay_bdt = predict_in_chunks(final_bdt_predict, ox, args.predict_chunk_rows)
        overlay_mlp = predict_in_chunks(final_mlp_predict, ox, args.predict_chunk_rows)
        overlay_stack_columns = {
            "bdt_score": overlay_bdt,
            "mlp_score": overlay_mlp,
            "cluster_Et": overlay_frame["cluster_Et"].to_numpy(dtype=MATRIX_DTYPE),
        }
        if args.domain == "auau":
            overlay_stack_columns["centrality"] = overlay_frame["centrality"].to_numpy(dtype=MATRIX_DTYPE)
        overlay_scores = {"BDT": overlay_bdt, "MLP": overlay_mlp}
        for name, (stack_features, predict) in stack_predictors.items():
            sx_overlay = np.column_stack([overlay_stack_columns[col] for col in stack_features]).astype(MATRIX_DTYPE, copy=False)
            overlay_scores[name] = predict_in_chunks(predict, sx_overlay, args.predict_chunk_rows)
        overlay_eval_mask = np.ones(len(overlay_frame), dtype=bool)
        hist_rows = histogram_rows(
            args.domain,
            overlay_frame,
            overlay_scores,
            overlay_weights,
            overlay_eval_mask,
            signal_sources,
            inclusive_sources,
            bins=args.overlay_bins,
        )
        overlay_source = {
            "histogram_source": "separate_overlay_input",
            "overlay_input_paths": overlay_paths,
            "overlay_weights": "unit weights for shape-only unit-area density diagnostics",
        }
    else:
        hist_rows = histogram_rows(args.domain, frame, model_scores, weights, test_mask, signal_sources, inclusive_sources, bins=args.overlay_bins)
        overlay_source = {
            "histogram_source": "locked_test_partition",
            "overlay_input_paths": [],
            "overlay_weights": "training/evaluation weights",
        }
    write_csv(args.outdir / "overlay_histograms.csv", hist_rows)
    correlations = score_correlations(frame, {"BDT": final_bdt_score, "MLP": final_mlp_score, **stack_scores}, test_mask, signal_sources, inclusive_sources)
    write_json(args.outdir / "score_correlations.json", correlations)

    table_payload = {
        "y": y[test_mask],
        "weight": weights[test_mask],
        "source_sample": frame.loc[test_mask, "source_sample"].astype(str).to_numpy(),
        "cluster_Et": frame.loc[test_mask, "cluster_Et"].to_numpy(dtype="float64"),
        "centrality": frame.loc[test_mask, "centrality"].to_numpy(dtype="float64"),
        "BDT": final_bdt_score[test_mask],
        "MLP": final_mlp_score[test_mask],
    }
    for name, score in stack_scores.items():
        table_payload[name] = score[test_mask]
    np.savez_compressed(args.outdir / "locked_test_score_table.npz", **table_payload)

    stack_missing_by_model = {
        name: int((test_mask & ~np.isfinite(score)).sum())
        for name, score in stack_scores.items()
    }
    leakage_qa = {
        "schema": "RJ_STACK_LEAKAGE_QA_V2",
        "stack_training_mode": args.stack_training_mode,
        "stack_training_contract": stack_training_contract,
        "base_score_contract": "Stack-training base scores are produced by base models whose gradient-training mask excludes the scored event key.",
        "locked_test_contract": "Locked-test base scores are produced by final base models trained only on non-test events.",
        "training_class_definition": TRAINING_CLASS_DEFINITION,
        "source_sample_used_as_supervised_label": False,
        "event_key_overlap_trainval_locked_test": event_key_overlap_count(keys, trainval_mask, test_mask),
        "stack_score_missing_bdt_rows": int(np.isnan(oof_bdt[stack_score_mask]).sum()),
        "stack_score_missing_mlp_rows": int(np.isnan(oof_mlp[stack_score_mask]).sum()),
        "locked_test_missing_bdt_rows": int(np.isnan(final_bdt_score[test_mask]).sum()),
        "locked_test_missing_mlp_rows": int(np.isnan(final_mlp_score[test_mask]).sum()),
        "locked_test_missing_stack_rows_by_model": stack_missing_by_model,
        "score_columns_finite": bool(
            np.isfinite(oof_bdt[stack_score_mask]).all()
            and np.isfinite(oof_mlp[stack_score_mask]).all()
            and np.isfinite(final_bdt_score[test_mask]).all()
            and np.isfinite(final_mlp_score[test_mask]).all()
            and all(value == 0 for value in stack_missing_by_model.values())
        ),
        "split_counts": {
            "trainval": partition_qa["partitions"]["trainval"],
            "locked_test": partition_qa["partitions"]["locked_test"],
            "stack_score_region": {
                "rows": int(stack_score_mask.sum()),
                "events": event_count(keys, stack_score_mask),
                "event_digest": event_key_digest(keys, stack_score_mask),
                "label_counts": label_counts(y, stack_score_mask),
                "source_counts": source_counts(frame, stack_score_mask),
            },
        },
        "status": "pass",
    }
    missing_keys = (
        "stack_score_missing_bdt_rows",
        "stack_score_missing_mlp_rows",
        "locked_test_missing_bdt_rows",
        "locked_test_missing_mlp_rows",
        "event_key_overlap_trainval_locked_test",
    )
    if any(leakage_qa[key] for key in missing_keys) or not leakage_qa["score_columns_finite"]:
        leakage_qa["status"] = "fail"
        raise SystemExit(f"Score leakage/missing QA failed: {leakage_qa}")
    write_json(args.outdir / "leakage_qa.json", leakage_qa)

    manifest = {
        "schema": SCHEMA,
        "status": "READY",
        "domain": args.domain,
        "created_unix": time.time(),
        "elapsed_seconds": time.time() - start,
        "input_paths": input_paths,
        "from_cache": bool(from_cache),
        "optional_columns_seen": optional_seen,
        "sample_definitions": {
            "training_signal_sources": signal_sources,
            "training_inclusive_sources": parse_csv_list(args.training_inclusive_sources),
            "overlay_signal_sources": signal_sources,
            "overlay_inclusive_sources": inclusive_sources,
            **overlay_source,
        },
        "training_class_definition": TRAINING_CLASS_DEFINITION,
        "overlay_class_definition": OVERLAY_CLASS_DEFINITION,
        "feature_contract": feature_metadata,
        "partition_seed": int(args.random_seed),
        "locked_test_fraction": float(args.locked_test_fraction),
        "stack_training_mode": args.stack_training_mode,
        "stack_training_contract": stack_training_contract,
        "oof_fold_count": int(args.folds),
        "stack_validation_fold": stack_val_fold,
        "simple_stack_fold": simple_stack_fold if args.stack_training_mode == "simple_holdout" else None,
        "metric_convention": {
            "primary_auc": "weighted_auc from train_auau_photon_mlp.auc_score with PPG12-exact weights when requested",
            "also_recorded": ["unweighted_auc", f"WP{int(args.target_signal_efficiency*100)} signal-efficiency threshold/fake-rate"],
        },
        "class_definition_checks": {
            "training_uses_source_sample_as_label": False,
            "source_sample_role": "provenance and overlay/sample QA only",
            "source_sample_source": "ROOT branch when present; otherwise deterministic path-derived sample name from train_auau_photon_bdt.infer_source_sample",
        },
        "artifacts": {
            "partition_qa": str(args.outdir / "partition_qa.json"),
            "feature_contract": str(args.outdir / "feature_contract.json"),
            "base_model_artifacts": str(args.outdir / "base_model_artifacts.json"),
            "stack_model_artifacts": str(args.outdir / "stack_model_artifacts.json"),
            "model_metrics_csv": str(args.outdir / "model_metrics.csv"),
            "model_metrics_json": str(args.outdir / "model_metrics.json"),
            "stratified_metrics_csv": str(args.outdir / "stratified_metrics.csv"),
            "overlay_histograms_csv": str(args.outdir / "overlay_histograms.csv"),
            "score_correlations_json": str(args.outdir / "score_correlations.json"),
            "locked_test_score_table": str(args.outdir / "locked_test_score_table.npz"),
            "leakage_qa": str(args.outdir / "leakage_qa.json"),
        },
        "partition_qa": partition_qa,
        "weight_report": weight_report,
        "leakage_qa": leakage_qa,
    }
    write_json(args.outdir / "campaign_manifest.json", manifest)
    return manifest


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--domain", choices=["pp", "auau"], required=True)
    ap.add_argument("--input", nargs="*", default=[])
    ap.add_argument("--overlay-input", nargs="*", default=[], help="Optional separate source-defined overlay input; never used for training.")
    ap.add_argument("--tree", default="DecayPhotonInfo")
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--feature-preset", choices=["pp_basev3e_noiso", "auau_global_noiso"], default=None)
    ap.add_argument("--features", default="")
    ap.add_argument("--label-branch", default="is_signal")
    ap.add_argument("--pt-range", default="")
    ap.add_argument("--centrality-range", default="")
    ap.add_argument("--cache-file", type=Path, default=None)
    ap.add_argument("--overlay-cache-file", type=Path, default=None)
    ap.add_argument("--skip-missing-tree", action="store_true")
    ap.add_argument("--max-load-rows-per-class", type=int, default=0)
    ap.add_argument("--max-load-rows", type=int, default=0)
    ap.add_argument("--load-sample-seed", type=int, default=42)
    ap.add_argument("--weight-mode", choices=["ppg12-exact", "none"], default="ppg12-exact")
    ap.add_argument("--ppg12-exact-expected-samples", default="")
    ap.add_argument("--ppg12-exact-closure-dir", type=Path, default=None)
    ap.add_argument("--training-inclusive-sources", default="")
    ap.add_argument("--overlay-signal-sources", default="")
    ap.add_argument("--overlay-inclusive-sources", default="")
    ap.add_argument("--locked-test-fraction", type=float, default=0.20)
    ap.add_argument("--folds", type=int, default=5)
    ap.add_argument(
        "--stack-training-mode",
        choices=["oof5", "simple_holdout"],
        default="oof5",
        help="How honest BDT/MLP base scores are produced for stacker training.",
    )
    ap.add_argument("--stack-validation-fold", type=int, default=0)
    ap.add_argument(
        "--simple-stack-fold",
        type=int,
        default=None,
        help="Trainval fold held out for simple-holdout stack training. Defaults to --stack-validation-fold.",
    )
    ap.add_argument("--random-seed", type=int, default=260604)
    ap.add_argument("--target-signal-efficiency", type=float, default=0.80)
    ap.add_argument("--report-et-bins", default="15,17,19,21,23,25,27,30,35")
    ap.add_argument("--report-cent-bins", default="0,20,40,60,80")
    ap.add_argument("--overlay-bins", type=int, default=60)
    ap.add_argument("--n-jobs", type=int, default=4)
    ap.add_argument(
        "--predict-chunk-rows",
        type=int,
        default=750000,
        help="Rows per prediction chunk for memory-bounded scoring. Set <=0 to score each region at once.",
    )
    ap.add_argument(
        "--reuse-existing-base-models",
        action="store_true",
        help="Load completed base BDT/MLP artifacts from --outdir/base_models when their partition contract matches, then train only missing base models.",
    )
    ap.add_argument("--bdt-estimators", type=int, default=750)
    ap.add_argument("--bdt-max-depth", type=int, default=5)
    ap.add_argument("--bdt-learning-rate", type=float, default=0.1)
    ap.add_argument("--bdt-subsample", type=float, default=0.5)
    ap.add_argument("--bdt-colsample-bytree", type=float, default=0.6)
    ap.add_argument("--bdt-tree-method", default="hist")
    ap.add_argument("--bdt-reg-alpha", type=float, default=5.0)
    ap.add_argument("--bdt-reg-lambda", type=float, default=0.3)
    ap.add_argument("--bdt-grow-policy", default="lossguide")
    ap.add_argument("--bdt-max-bin", type=int, default=256)
    ap.add_argument("--allow-bdt-linear-fallback", action="store_true", help="Diagnostic smoke fallback only if xgboost/sklearn are unavailable.")
    ap.add_argument("--stack-algorithms", default="logistic,gbm,mlp")
    ap.add_argument("--stack-l2", type=float, default=2.0e-3)
    ap.add_argument("--stack-linear-backend", choices=["numpy", "sklearn"], default="numpy")
    ap.add_argument("--stack-max-linear-steps", type=int, default=1800)
    ap.add_argument("--stack-gbm-estimators", type=int, default=90)
    ap.add_argument("--stack-gbm-learning-rate", type=float, default=0.045)
    ap.add_argument("--stack-gbm-max-depth", type=int, default=3)
    ap.add_argument("--stack-gbm-max-leaf-nodes", type=int, default=8)
    ap.add_argument("--stack-mlp-hidden", default="256,128,64")
    ap.add_argument("--stack-mlp-epochs", type=int, default=180)
    ap.add_argument("--stack-mlp-patience", type=int, default=24)
    ap.add_argument("--stack-mlp-batch-size", type=int, default=8192)
    ap.add_argument("--stack-mlp-learning-rate", type=float, default=1.5e-3)
    ap.add_argument("--stack-mlp-l2", type=float, default=1.0e-3)
    ap.add_argument("--preflight-only", action="store_true")
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--self-test-rows", type=int, default=5000)
    args = ap.parse_args()
    if args.stack_validation_fold < 0 or args.stack_validation_fold >= args.folds:
        raise SystemExit("--stack-validation-fold must be in [0, folds)")
    if args.simple_stack_fold is not None and (args.simple_stack_fold < 0 or args.simple_stack_fold >= args.folds):
        raise SystemExit("--simple-stack-fold must be in [0, folds)")
    if args.domain == "pp" and not args.overlay_signal_sources:
        args.overlay_signal_sources = ",".join(PP_SIGNAL_SOURCES)
    if args.domain == "pp" and not args.training_inclusive_sources:
        args.training_inclusive_sources = ",".join(PP_TRAIN_INCLUSIVE_SOURCES)
    if args.domain == "pp" and not args.overlay_inclusive_sources:
        args.overlay_inclusive_sources = ",".join(PP_OVERLAY_INCLUSIVE_SOURCES)
    if args.domain == "auau" and not args.overlay_signal_sources:
        args.overlay_signal_sources = ",".join(AUAU_SIGNAL_SOURCES)
    if args.domain == "auau" and not args.training_inclusive_sources:
        args.training_inclusive_sources = ",".join(AUAU_INCLUSIVE_SOURCES)
    if args.domain == "auau" and not args.overlay_inclusive_sources:
        args.overlay_inclusive_sources = ",".join(AUAU_INCLUSIVE_SOURCES)
    if not args.ppg12_exact_expected_samples:
        args.ppg12_exact_expected_samples = ",".join(
            parse_csv_list(args.overlay_signal_sources) + parse_csv_list(args.training_inclusive_sources)
        )
    if args.ppg12_exact_closure_dir is None:
        args.ppg12_exact_closure_dir = args.outdir / "ppg12_exact_reweight_closure"
    return args


def main() -> int:
    args = parse_args()
    manifest = train_and_score(args)
    print(json.dumps({"status": manifest.get("status"), "domain": args.domain, "outdir": str(args.outdir)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
