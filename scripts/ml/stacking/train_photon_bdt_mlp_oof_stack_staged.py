#!/usr/bin/env python3
"""Stage-wise Fresh pp/AuAu BDT+MLP stack campaign runner.

This is the Condor-DAG friendly execution shape for THE-35.  It preserves the
science contract from ``train_photon_bdt_mlp_oof_stack.py`` while splitting the
work into bounded stages:

``build-matrix -> train-base -> score-base -> train-stack -> reduce-domain``.
"""

from __future__ import annotations

import argparse
import gc
import hashlib
import json
import math
import pickle
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[3]
STACKING_DIR = REPO_ROOT / "scripts" / "ml" / "stacking"
if str(STACKING_DIR) not in sys.path:
    sys.path.insert(0, str(STACKING_DIR))

import train_photon_bdt_mlp_oof_stack as mono  # noqa: E402
from train_auau_stacked_bdt_mlp_sweep import FittedModel, fit_one  # noqa: E402


STAGED_SCHEMA = "RJ_FRESH_PP_AUAU_BDT_MLP_STACK_STAGED_V1"
MATRIX_SCHEMA = "RJ_FRESH_STACK_MATRIX_V1"
SCORE_SCHEMA = "RJ_FRESH_STACK_SCORE_SHARD_V1"
MATRIX_DTYPE = np.float32
SCORE_DTYPE = np.float32


@dataclass(frozen=True)
class MatrixStore:
    root: Path
    manifest: dict
    x: np.ndarray
    y: np.ndarray
    weights: np.ndarray
    locked_test_mask: np.ndarray
    fold_id: np.ndarray
    row_event_id: np.ndarray
    event_hash: np.ndarray
    source_code: np.ndarray
    source_labels: list[str]
    cluster_et: np.ndarray
    centrality: np.ndarray
    run: np.ndarray
    evt: np.ndarray

    @property
    def trainval_mask(self) -> np.ndarray:
        return ~self.locked_test_mask


def write_json(path: Path, payload) -> None:
    mono.write_json(path, payload)


def read_json(path: Path) -> dict:
    return json.loads(path.read_text())


def parse_csv_list(text: str | None) -> list[str]:
    return mono.parse_csv_list(text)


def matrix_dir_for(args: argparse.Namespace) -> Path:
    if getattr(args, "matrix_dir", None):
        return Path(args.matrix_dir)
    return args.outdir / "matrix"


def load_matrix_dir(matrix_dir: Path, mmap_mode: str | None = "r") -> MatrixStore:
    manifest_path = matrix_dir / "matrix_manifest.json"
    if not manifest_path.is_file():
        raise SystemExit(f"Missing matrix manifest: {manifest_path}")
    manifest = read_json(manifest_path)
    source_payload = read_json(matrix_dir / "source_labels.json")
    return MatrixStore(
        root=matrix_dir,
        manifest=manifest,
        x=np.load(matrix_dir / "x.npy", mmap_mode=mmap_mode),
        y=np.load(matrix_dir / "y.npy", mmap_mode=mmap_mode),
        weights=np.load(matrix_dir / "weights.npy", mmap_mode=mmap_mode),
        locked_test_mask=np.load(matrix_dir / "locked_test_mask.npy", mmap_mode=mmap_mode),
        fold_id=np.load(matrix_dir / "fold_id.npy", mmap_mode=mmap_mode),
        row_event_id=np.load(matrix_dir / "row_event_id.npy", mmap_mode=mmap_mode),
        event_hash=np.load(matrix_dir / "event_hash.npy", mmap_mode=mmap_mode),
        source_code=np.load(matrix_dir / "source_code.npy", mmap_mode=mmap_mode),
        source_labels=[str(item) for item in source_payload["labels"]],
        cluster_et=np.load(matrix_dir / "cluster_Et.npy", mmap_mode=mmap_mode),
        centrality=np.load(matrix_dir / "centrality.npy", mmap_mode=mmap_mode),
        run=np.load(matrix_dir / "run.npy", mmap_mode=mmap_mode),
        evt=np.load(matrix_dir / "evt.npy", mmap_mode=mmap_mode),
    )


def load_matrix(args: argparse.Namespace, mmap_mode: str | None = "r") -> MatrixStore:
    return load_matrix_dir(matrix_dir_for(args), mmap_mode=mmap_mode)


def sync_static_artifacts(args: argparse.Namespace, matrix: MatrixStore) -> None:
    feature_path = args.outdir / "feature_contract.json"
    partition_path = args.outdir / "partition_qa.json"
    if not feature_path.is_file():
        write_json(feature_path, matrix.manifest["feature_contract"])
    if not partition_path.is_file():
        src = matrix.root.parent / "partition_qa.json"
        if src.is_file():
            partition_path.write_text(src.read_text())


def matrix_to_frame(matrix: MatrixStore):
    import pandas as pd

    return pd.DataFrame(
        {
            "source_sample": pd.Categorical.from_codes(
                np.asarray(matrix.source_code, dtype="int32"),
                categories=matrix.source_labels,
            ),
            "run": np.asarray(matrix.run, dtype="int64"),
            "evt": np.asarray(matrix.evt, dtype="int64"),
            "is_signal": np.asarray(matrix.y, dtype="int32"),
            "cluster_Et": np.asarray(matrix.cluster_et, dtype="float32"),
            "centrality": np.asarray(matrix.centrality, dtype="float32"),
        }
    )


def source_codes(frame) -> tuple[np.ndarray, list[str]]:
    import pandas as pd

    series = frame["source_sample"]
    if isinstance(series.dtype, pd.CategoricalDtype):
        return series.cat.codes.to_numpy(dtype="int16", copy=False), [str(x) for x in series.cat.categories.tolist()]
    codes, labels = pd.factorize(series.astype(str), sort=True)
    return codes.astype("int16", copy=False), [str(x) for x in labels.tolist()]


def event_hashes(keys: mono.EventKeyIndex) -> np.ndarray:
    return np.asarray(
        [hashlib.blake2b(key.encode("utf-8"), digest_size=16).hexdigest() for key in keys.key_strings],
        dtype="U32",
    )


def event_digest_from_ids(matrix: MatrixStore, mask: np.ndarray) -> dict[str, object]:
    unique_ids = np.unique(np.asarray(matrix.row_event_id)[mask])
    hashes = np.asarray(matrix.event_hash)[unique_ids]
    digest = hashlib.blake2b("\n".join(str(x) for x in hashes).encode("utf-8"), digest_size=16).hexdigest()
    return {"events": int(len(unique_ids)), "rows": int(np.asarray(mask, dtype=bool).sum()), "digest_blake2b16": digest}


def event_overlap_from_ids(matrix: MatrixStore, left: np.ndarray, right: np.ndarray) -> int:
    left_ids = np.unique(np.asarray(matrix.row_event_id)[left])
    right_ids = np.unique(np.asarray(matrix.row_event_id)[right])
    return int(len(np.intersect1d(left_ids, right_ids, assume_unique=False)))


def model_spec(args: argparse.Namespace, matrix: MatrixStore) -> tuple[str, np.ndarray, np.ndarray, int]:
    trainval = matrix.trainval_mask
    fold_id = np.asarray(matrix.fold_id)
    stack_val_fold = int(args.stack_validation_fold)
    simple_fold = int(args.simple_stack_fold if args.simple_stack_fold is not None else args.stack_validation_fold)
    if args.model_role == "fold":
        if args.stack_training_mode != "oof5":
            raise SystemExit("--model-role fold is only valid with --stack-training-mode oof5")
        fold = int(args.fold)
        name = f"fold{fold}_{args.model_kind}"
        train_mask = trainval & (fold_id != fold)
        val_mask = trainval & (fold_id == fold)
        seed = int(args.random_seed + (fold + 1 if args.model_kind == "bdt" else 100 + fold))
        return name, train_mask, val_mask, seed
    if args.model_role == "simple":
        if args.stack_training_mode != "simple_holdout":
            raise SystemExit("--model-role simple is only valid with --stack-training-mode simple_holdout")
        name = f"simple_holdout_fold{simple_fold}_{args.model_kind}"
        train_mask = trainval & (fold_id != simple_fold)
        val_mask = trainval & (fold_id == simple_fold)
        seed = int(args.random_seed + (501 if args.model_kind == "bdt" else 601))
        return name, train_mask, val_mask, seed
    if args.model_role == "final":
        name = f"final_trainval_{args.model_kind}"
        train_mask = trainval
        val_mask = trainval & (fold_id == stack_val_fold)
        seed = int(args.random_seed + (999 if args.model_kind == "bdt" else 1999))
        return name, train_mask, val_mask, seed
    raise SystemExit(f"Unsupported model role: {args.model_role}")


def score_region_mask(args: argparse.Namespace, matrix: MatrixStore) -> np.ndarray:
    trainval = matrix.trainval_mask
    fold_id = np.asarray(matrix.fold_id)
    simple_fold = int(args.simple_stack_fold if args.simple_stack_fold is not None else args.stack_validation_fold)
    if args.score_region == "locked_test":
        return np.asarray(matrix.locked_test_mask, dtype=bool)
    if args.score_region == "fold":
        return trainval & (fold_id == int(args.fold))
    if args.score_region == "simple_holdout":
        return trainval & (fold_id == simple_fold)
    raise SystemExit(f"Unsupported score region: {args.score_region}")


def score_path(outdir: Path, model_name: str, region: str) -> Path:
    return outdir / "scores" / f"{model_name}__{region}.npz"


def load_score(outdir: Path, model_name: str, region: str, n_rows: int) -> np.ndarray:
    path = score_path(outdir, model_name, region)
    if not path.is_file():
        raise SystemExit(f"Missing score shard: {path}")
    with np.load(path) as payload:
        row_index = payload["row_index"].astype("int64", copy=False)
        score = payload["score"].astype(SCORE_DTYPE, copy=False)
    out = np.full(n_rows, np.nan, dtype=SCORE_DTYPE)
    out[row_index] = score
    return out


def load_base_scores(args: argparse.Namespace, matrix: MatrixStore) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    n_rows = len(matrix.y)
    oof_bdt = np.full(n_rows, np.nan, dtype=SCORE_DTYPE)
    oof_mlp = np.full(n_rows, np.nan, dtype=SCORE_DTYPE)
    if args.stack_training_mode == "oof5":
        for fold in range(int(args.folds)):
            for kind, target in (("bdt", oof_bdt), ("mlp", oof_mlp)):
                name = f"fold{fold}_{kind}"
                score = load_score(args.outdir, name, "fold", n_rows)
                mask = np.isfinite(score)
                target[mask] = score[mask]
    else:
        simple_fold = int(args.simple_stack_fold if args.simple_stack_fold is not None else args.stack_validation_fold)
        for kind, target in (("bdt", oof_bdt), ("mlp", oof_mlp)):
            name = f"simple_holdout_fold{simple_fold}_{kind}"
            score = load_score(args.outdir, name, "simple_holdout", n_rows)
            mask = np.isfinite(score)
            target[mask] = score[mask]
    final_bdt = load_score(args.outdir, "final_trainval_bdt", "locked_test", n_rows)
    final_mlp = load_score(args.outdir, "final_trainval_mlp", "locked_test", n_rows)
    return oof_bdt, oof_mlp, final_bdt, final_mlp


def load_stack_predictors(outdir: Path) -> dict[str, tuple[list[str], FittedModel]]:
    artifact_path = outdir / "stack_model_artifacts.json"
    if not artifact_path.is_file():
        raise SystemExit(f"Missing stack artifact index: {artifact_path}")
    payload = read_json(artifact_path)
    predictors: dict[str, tuple[list[str], FittedModel]] = {}
    for item in payload.get("models", []):
        pickle_path = Path(item["pickle_path"])
        with pickle_path.open("rb") as handle:
            pickled = pickle.load(handle)
        artifact = read_json(Path(item["artifact_path"]))
        model = FittedModel(
            name=str(item["model"]),
            algorithm=str(pickled["algorithm"]),
            feature_names=list(pickled["feature_names"]),
            artifact=artifact,
            model_object=pickled["model"],
            history=list(item.get("history", [])),
        )
        predictors[str(item["model"])] = (list(item["features"]), model)
    return predictors


def stage_build_matrix(args: argparse.Namespace) -> dict:
    start = time.time()
    args.outdir.mkdir(parents=True, exist_ok=True)
    matrix_dir = args.outdir / "matrix"
    matrix_dir.mkdir(parents=True, exist_ok=True)
    contract = mono.feature_contract(args)
    frame, optional_seen, from_cache, input_paths = mono.load_frame_for_campaign(args, contract.features)
    frame = mono.filter_rows(frame, args)
    frame = mono.compact_campaign_frame(frame)
    missing_features = [feature for feature in contract.features if feature not in frame.columns]
    if missing_features:
        raise SystemExit(f"Selected feature(s) missing after derived-feature build: {missing_features}")
    y = frame[args.label_branch].to_numpy(dtype="int32")
    weights, weight_report = mono.prepare_weights(frame, args)
    locked_test_mask, fold_id, partition_qa = mono.assign_partitions(frame, args.random_seed, args.locked_test_fraction, args.folds)
    keys = mono.event_keys(frame)
    src_codes, src_labels = source_codes(frame)
    feature_metadata = {
        "preset": contract.preset,
        "description": contract.description,
        "features": contract.features,
        "n_features": len(contract.features),
    }
    x = mono.np_matrix(frame, contract.features)
    np.save(matrix_dir / "x.npy", np.asarray(x, dtype=MATRIX_DTYPE))
    np.save(matrix_dir / "y.npy", np.asarray(y, dtype="int8"))
    np.save(matrix_dir / "weights.npy", np.asarray(weights, dtype="float64"))
    np.save(matrix_dir / "locked_test_mask.npy", np.asarray(locked_test_mask, dtype=bool))
    np.save(matrix_dir / "fold_id.npy", np.asarray(fold_id, dtype="int16"))
    np.save(matrix_dir / "row_event_id.npy", np.asarray(keys.row_ids, dtype="int64"))
    np.save(matrix_dir / "event_hash.npy", event_hashes(keys))
    np.save(matrix_dir / "source_code.npy", np.asarray(src_codes, dtype="int16"))
    np.save(matrix_dir / "cluster_Et.npy", frame["cluster_Et"].to_numpy(dtype=MATRIX_DTYPE))
    np.save(matrix_dir / "centrality.npy", frame["centrality"].to_numpy(dtype=MATRIX_DTYPE))
    np.save(matrix_dir / "run.npy", frame["run"].to_numpy(dtype="int64"))
    np.save(matrix_dir / "evt.npy", frame["evt"].to_numpy(dtype="int64"))
    write_json(matrix_dir / "source_labels.json", {"labels": src_labels})
    write_json(args.outdir / "partition_qa.json", partition_qa)
    write_json(args.outdir / "feature_contract.json", feature_metadata)
    manifest = {
        "schema": MATRIX_SCHEMA,
        "status": "READY",
        "domain": args.domain,
        "created_unix": time.time(),
        "elapsed_seconds": time.time() - start,
        "input_paths": input_paths,
        "from_cache": bool(from_cache),
        "optional_columns_seen": optional_seen,
        "n_rows": int(len(frame)),
        "n_events": int(keys.n_events),
        "x_shape": list(x.shape),
        "x_dtype": str(x.dtype),
        "feature_contract": feature_metadata,
        "partition_seed": int(args.random_seed),
        "locked_test_fraction": float(args.locked_test_fraction),
        "folds": int(args.folds),
        "weight_report": weight_report,
        "training_class_definition": mono.TRAINING_CLASS_DEFINITION,
        "overlay_class_definition": mono.OVERLAY_CLASS_DEFINITION,
        "artifacts": {
            "x": str(matrix_dir / "x.npy"),
            "y": str(matrix_dir / "y.npy"),
            "weights": str(matrix_dir / "weights.npy"),
            "locked_test_mask": str(matrix_dir / "locked_test_mask.npy"),
            "fold_id": str(matrix_dir / "fold_id.npy"),
            "row_event_id": str(matrix_dir / "row_event_id.npy"),
            "event_hash": str(matrix_dir / "event_hash.npy"),
        },
    }
    write_json(matrix_dir / "matrix_manifest.json", manifest)
    print(json.dumps({"stage": "build-matrix", "status": "READY", "domain": args.domain, "rows": len(frame)}, sort_keys=True))
    return manifest


def stage_train_base(args: argparse.Namespace) -> dict:
    matrix = load_matrix(args, mmap_mode="r")
    contract = matrix.manifest["feature_contract"]
    name, train_mask, val_mask, seed = model_spec(args, matrix)
    args.self_test = bool(args.self_test or args.allow_bdt_linear_fallback)
    if args.model_kind == "bdt":
        artifact, _ = mono.fit_bdt(name, matrix.x, matrix.y.astype("int32"), train_mask, matrix.weights, args, args.outdir, seed)
    else:
        artifact, _ = mono.train_mlp_model(
            name,
            list(contract["features"]),
            matrix.x,
            matrix.y.astype("int32"),
            train_mask,
            val_mask,
            matrix.weights,
            args,
            args.outdir,
            seed,
        )
    done = {
        "schema": STAGED_SCHEMA,
        "stage": "train-base",
        "status": "READY",
        "domain": args.domain,
        "model_name": name,
        "model_kind": args.model_kind,
        "model_role": args.model_role,
        "fold": args.fold,
        "train_rows": int(train_mask.sum()),
        "val_rows": int(val_mask.sum()),
        "artifact": artifact,
    }
    write_json(args.outdir / "stage_manifests" / f"train_base__{name}.json", done)
    print(json.dumps({"stage": "train-base", "status": "READY", "model": name}, sort_keys=True))
    return done


def load_predictor(args: argparse.Namespace, matrix: MatrixStore, name: str, kind: str) -> Callable[[np.ndarray], np.ndarray]:
    features = list(matrix.manifest["feature_contract"]["features"])
    y = matrix.y.astype("int32")
    if kind == "bdt":
        loaded = mono.load_existing_bdt(name, args.outdir, None, y)
    else:
        loaded = mono.load_existing_mlp(name, features, args.outdir, None, y)
    if loaded is None:
        raise SystemExit(f"Missing trained base model for scoring: {name}")
    return loaded[1]


def stage_score_base(args: argparse.Namespace) -> dict:
    matrix = load_matrix(args, mmap_mode="r")
    if args.model_name:
        name = args.model_name
        kind = args.model_kind
    else:
        name, _, _, _ = model_spec(args, matrix)
        kind = args.model_kind
    mask = score_region_mask(args, matrix)
    predictor = load_predictor(args, matrix, name, kind)
    row_index = np.flatnonzero(mask).astype("int64")
    score = mono.predict_in_chunks(predictor, matrix.x[row_index], args.predict_chunk_rows).astype(SCORE_DTYPE, copy=False)
    out = score_path(args.outdir, name, args.score_region)
    out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        out,
        schema=SCORE_SCHEMA,
        model_name=name,
        model_kind=kind,
        score_region=args.score_region,
        row_index=row_index,
        score=score,
    )
    done = {
        "schema": STAGED_SCHEMA,
        "stage": "score-base",
        "status": "READY",
        "domain": args.domain,
        "model_name": name,
        "model_kind": kind,
        "score_region": args.score_region,
        "rows": int(len(row_index)),
        "score_path": str(out),
    }
    write_json(args.outdir / "stage_manifests" / f"score_base__{name}__{args.score_region}.json", done)
    print(json.dumps({"stage": "score-base", "status": "READY", "model": name, "rows": len(row_index)}, sort_keys=True))
    return done


def stack_masks(args: argparse.Namespace, matrix: MatrixStore) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, dict]:
    trainval = matrix.trainval_mask
    fold_id = np.asarray(matrix.fold_id)
    stack_val_fold = int(args.stack_validation_fold)
    simple_fold = int(args.simple_stack_fold if args.simple_stack_fold is not None else args.stack_validation_fold)
    if args.stack_training_mode == "oof5":
        stack_val_mask = trainval & (fold_id == stack_val_fold)
        stack_train_mask = trainval & (fold_id != stack_val_fold)
        stack_score_mask = trainval
        stack_eval_extra_mask = stack_val_mask
        contract = {
            "mode": "oof5",
            "description": "Stackers train on five-fold out-of-fold BDT/MLP base scores over the trainval region.",
            "base_score_region": "all non-test trainval folds, each scored by the fold model that excluded that fold from base training",
            "stack_train_region": f"OOF folds excluding validation fold {stack_val_fold}",
            "stack_validation_region": f"OOF fold {stack_val_fold}",
            "fold_count": int(args.folds),
            "stack_validation_fold": stack_val_fold,
        }
    elif args.stack_training_mode == "simple_holdout":
        pred_mask = trainval & (fold_id == simple_fold)
        stack_val_mask = np.zeros(len(matrix.y), dtype=bool)
        stack_train_mask = pred_mask
        stack_score_mask = pred_mask
        stack_eval_extra_mask = pred_mask
        contract = {
            "mode": "simple_holdout",
            "description": "Stackers train on one held-out trainval block scored by BDT/MLP base models trained on the other trainval blocks.",
            "base_score_region": f"trainval fold {simple_fold} only",
            "base_train_region": f"trainval folds excluding fold {simple_fold}",
            "stack_train_region": f"scored trainval fold {simple_fold}",
            "stack_validation_region": "stacker MLP uses an internal split inside the simple stack-training block when needed for early stopping",
            "fold_count": int(args.folds),
            "simple_stack_fold": simple_fold,
        }
    else:
        raise SystemExit(f"Unsupported stack mode: {args.stack_training_mode}")
    return stack_train_mask, stack_val_mask, stack_score_mask, stack_eval_extra_mask, contract


def stack_input_columns(args: argparse.Namespace, matrix: MatrixStore, oof_bdt, oof_mlp, final_bdt, final_mlp) -> tuple[dict[str, np.ndarray], dict[str, list[str]]]:
    test_mask = np.asarray(matrix.locked_test_mask, dtype=bool)
    stack_columns = {
        "bdt_score": np.where(test_mask, final_bdt, oof_bdt).astype(SCORE_DTYPE, copy=False),
        "mlp_score": np.where(test_mask, final_mlp, oof_mlp).astype(SCORE_DTYPE, copy=False),
        "cluster_Et": np.asarray(matrix.cluster_et, dtype=MATRIX_DTYPE),
    }
    if args.domain == "auau":
        stack_columns["centrality"] = np.asarray(matrix.centrality, dtype=MATRIX_DTYPE)
    specs = {
        "score_only": ["bdt_score", "mlp_score"],
        "score_context": ["bdt_score", "mlp_score", "cluster_Et"] + (["centrality"] if args.domain == "auau" else []),
    }
    return stack_columns, specs


def stage_train_stack(args: argparse.Namespace) -> dict:
    matrix = load_matrix(args, mmap_mode="r")
    oof_bdt, oof_mlp, final_bdt, final_mlp = load_base_scores(args, matrix)
    stack_train_mask, stack_val_mask, _, stack_eval_extra_mask, stack_contract = stack_masks(args, matrix)
    stack_columns, specs = stack_input_columns(args, matrix, oof_bdt, oof_mlp, final_bdt, final_mlp)
    y = matrix.y.astype("int32")
    algorithms = parse_csv_list(args.stack_algorithms)
    sargs = mono.stack_args_from(args)
    stack_artifacts = []
    for stack_name, stack_features in specs.items():
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
                sample_weights=matrix.weights,
            )
            if fitted is None:
                continue
            eval_mask = np.asarray(matrix.locked_test_mask, dtype=bool) | stack_eval_extra_mask | stack_val_mask
            row_index = np.flatnonzero(eval_mask).astype("int64")
            score = mono.predict_in_chunks(fitted.predict, sx[row_index], args.predict_chunk_rows).astype(SCORE_DTYPE, copy=False)
            score_out = score_path(args.outdir, model_name, "stack_eval")
            score_out.parent.mkdir(parents=True, exist_ok=True)
            np.savez_compressed(score_out, schema=SCORE_SCHEMA, model_name=model_name, score_region="stack_eval", row_index=row_index, score=score)
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
                    "train_region": stack_contract["stack_train_region"],
                    "validation_region": stack_contract["stack_validation_region"],
                    "artifact_path": str(artifact_path),
                    "pickle_path": str(pickle_path),
                    "score_path": str(score_out),
                    "history": fitted.history,
                }
            )
        del sx
        gc.collect()
    write_json(args.outdir / "stack_model_artifacts.json", {"schema": "RJ_STACK_MODEL_ARTIFACTS_V1", "models": stack_artifacts})
    done = {"schema": STAGED_SCHEMA, "stage": "train-stack", "status": "READY", "domain": args.domain, "models": stack_artifacts}
    write_json(args.outdir / "stage_manifests" / "train_stack.json", done)
    print(json.dumps({"stage": "train-stack", "status": "READY", "models": len(stack_artifacts)}, sort_keys=True))
    return done


def collect_base_artifacts(args: argparse.Namespace) -> dict:
    model_dir = args.outdir / "base_models"
    if args.stack_training_mode == "oof5":
        folds = []
        for fold in range(int(args.folds)):
            folds.append(
                {
                    "fold": fold,
                    "bdt": read_json(model_dir / f"fold{fold}_bdt.metadata.json"),
                    "mlp": read_json(model_dir / f"fold{fold}_mlp.metadata.json"),
                }
            )
        simple = None
    else:
        simple_fold = int(args.simple_stack_fold if args.simple_stack_fold is not None else args.stack_validation_fold)
        folds = []
        simple = {
            "heldout_fold": simple_fold,
            "bdt": read_json(model_dir / f"simple_holdout_fold{simple_fold}_bdt.metadata.json"),
            "mlp": read_json(model_dir / f"simple_holdout_fold{simple_fold}_mlp.metadata.json"),
        }
    final = {
        "bdt": read_json(model_dir / "final_trainval_bdt.metadata.json"),
        "mlp": read_json(model_dir / "final_trainval_mlp.metadata.json"),
    }
    return {"stack_training_mode": args.stack_training_mode, "folds": folds, "simple_holdout": simple, "final": final}


def overlay_histograms(args: argparse.Namespace, matrix: MatrixStore, final_bdt_predict, final_mlp_predict, stack_predictors) -> tuple[list[dict], dict]:
    contract = matrix.manifest["feature_contract"]
    signal_sources = parse_csv_list(args.overlay_signal_sources)
    inclusive_sources = parse_csv_list(args.overlay_inclusive_sources)
    overlay_frame, overlay_paths = mono.load_overlay_frame(args, list(contract["features"]))
    if overlay_frame is None:
        frame = matrix_to_frame(matrix)
        model_scores = {
            "BDT": load_score(args.outdir, "final_trainval_bdt", "locked_test", len(matrix.y)),
            "MLP": load_score(args.outdir, "final_trainval_mlp", "locked_test", len(matrix.y)),
        }
        for name in stack_predictors:
            model_scores[name] = load_score(args.outdir, name, "stack_eval", len(matrix.y))
        rows = mono.histogram_rows(args.domain, frame, model_scores, matrix.weights, matrix.locked_test_mask, signal_sources, inclusive_sources, bins=args.overlay_bins)
        return rows, {"histogram_source": "locked_test_partition", "overlay_input_paths": [], "overlay_weights": "training/evaluation weights"}
    overlay_frame = mono.compact_campaign_frame(overlay_frame)
    ox = mono.np_matrix(overlay_frame, list(contract["features"]))
    overlay_bdt = mono.predict_in_chunks(final_bdt_predict, ox, args.predict_chunk_rows)
    overlay_mlp = mono.predict_in_chunks(final_mlp_predict, ox, args.predict_chunk_rows)
    overlay_stack_columns = {
        "bdt_score": overlay_bdt,
        "mlp_score": overlay_mlp,
        "cluster_Et": overlay_frame["cluster_Et"].to_numpy(dtype=MATRIX_DTYPE),
    }
    if args.domain == "auau":
        overlay_stack_columns["centrality"] = overlay_frame["centrality"].to_numpy(dtype=MATRIX_DTYPE)
    overlay_scores = {"BDT": overlay_bdt, "MLP": overlay_mlp}
    for name, (stack_features, fitted) in stack_predictors.items():
        sx_overlay = np.column_stack([overlay_stack_columns[col] for col in stack_features]).astype(MATRIX_DTYPE, copy=False)
        overlay_scores[name] = mono.predict_in_chunks(fitted.predict, sx_overlay, args.predict_chunk_rows)
    rows = mono.histogram_rows(
        args.domain,
        overlay_frame,
        overlay_scores,
        np.ones(len(overlay_frame), dtype="float64"),
        np.ones(len(overlay_frame), dtype=bool),
        signal_sources,
        inclusive_sources,
        bins=args.overlay_bins,
    )
    return rows, {
        "histogram_source": "separate_overlay_input",
        "overlay_input_paths": overlay_paths,
        "overlay_weights": "unit weights for shape-only unit-area density diagnostics",
    }


def stage_reduce_domain(args: argparse.Namespace) -> dict:
    start = time.time()
    matrix = load_matrix(args, mmap_mode="r")
    sync_static_artifacts(args, matrix)
    frame = matrix_to_frame(matrix)
    y = matrix.y.astype("int32")
    weights = matrix.weights
    test_mask = np.asarray(matrix.locked_test_mask, dtype=bool)
    oof_bdt, oof_mlp, final_bdt, final_mlp = load_base_scores(args, matrix)
    _, stack_val_mask, stack_score_mask, stack_eval_extra_mask, stack_contract = stack_masks(args, matrix)
    stack_predictors = load_stack_predictors(args.outdir)
    stack_scores = {
        name: load_score(args.outdir, name, "stack_eval", len(y))
        for name in stack_predictors
    }
    model_scores = {"BDT": final_bdt, "MLP": final_mlp, **stack_scores}
    oof_scores = {"BDT": oof_bdt, "MLP": oof_mlp}
    base_region = "oof_trainval" if args.stack_training_mode == "oof5" else f"simple_holdout_stack_training_fold{int(args.simple_stack_fold if args.simple_stack_fold is not None else args.stack_validation_fold)}"
    metric_rows = []
    for name, score in model_scores.items():
        metric_rows.append(mono.metric_row(args.domain, "locked_test", name, y, score, weights, test_mask, args.target_signal_efficiency))
    for name, score in oof_scores.items():
        metric_rows.append(mono.metric_row(args.domain, base_region, name, y, score, weights, stack_score_mask, args.target_signal_efficiency))
    if stack_val_mask.any():
        for name, score in stack_scores.items():
            metric_rows.append(mono.metric_row(args.domain, f"oof_stack_validation_fold{args.stack_validation_fold}", name, y, score, weights, stack_val_mask, args.target_signal_efficiency))
    elif stack_eval_extra_mask.any():
        simple_fold = int(args.simple_stack_fold if args.simple_stack_fold is not None else args.stack_validation_fold)
        for name, score in stack_scores.items():
            metric_rows.append(mono.metric_row(args.domain, f"simple_holdout_stack_training_fold{simple_fold}_fit_diagnostic", name, y, score, weights, stack_eval_extra_mask, args.target_signal_efficiency))
    mono.write_csv(args.outdir / "model_metrics.csv", metric_rows)
    write_json(args.outdir / "model_metrics.json", {"schema": "RJ_MODEL_METRICS_V1", "rows": metric_rows})
    et_edges = mono.parse_edges(args.report_et_bins)
    cent_edges = mono.parse_edges(args.report_cent_bins) if args.domain == "auau" and args.report_cent_bins else None
    strat_rows = mono.stratified_metric_rows(args.domain, frame, model_scores, weights, test_mask, et_edges, cent_edges, args.target_signal_efficiency)
    mono.write_csv(args.outdir / "stratified_metrics.csv", strat_rows)
    final_bdt_predict = load_predictor(args, matrix, "final_trainval_bdt", "bdt")
    final_mlp_predict = load_predictor(args, matrix, "final_trainval_mlp", "mlp")
    hist_rows, overlay_source = overlay_histograms(args, matrix, final_bdt_predict, final_mlp_predict, stack_predictors)
    mono.write_csv(args.outdir / "overlay_histograms.csv", hist_rows)
    signal_sources = parse_csv_list(args.overlay_signal_sources)
    inclusive_sources = parse_csv_list(args.overlay_inclusive_sources)
    correlations = mono.score_correlations(frame, model_scores, test_mask, signal_sources, inclusive_sources)
    write_json(args.outdir / "score_correlations.json", correlations)
    table_payload = {
        "y": y[test_mask],
        "weight": weights[test_mask],
        "source_sample": frame.loc[test_mask, "source_sample"].astype(str).to_numpy(),
        "cluster_Et": matrix.cluster_et[test_mask].astype("float64"),
        "centrality": matrix.centrality[test_mask].astype("float64"),
        "BDT": final_bdt[test_mask],
        "MLP": final_mlp[test_mask],
    }
    for name, score in stack_scores.items():
        table_payload[name] = score[test_mask]
    np.savez_compressed(args.outdir / "locked_test_score_table.npz", **table_payload)
    stack_missing_by_model = {name: int((test_mask & ~np.isfinite(score)).sum()) for name, score in stack_scores.items()}
    partition_qa = read_json(args.outdir / "partition_qa.json")
    leakage_qa = {
        "schema": "RJ_STACK_LEAKAGE_QA_V2",
        "stack_training_mode": args.stack_training_mode,
        "stack_training_contract": stack_contract,
        "base_score_contract": "Stack-training base scores are produced by base models whose gradient-training mask excludes the scored event key.",
        "locked_test_contract": "Locked-test base scores are produced by final base models trained only on non-test events.",
        "training_class_definition": mono.TRAINING_CLASS_DEFINITION,
        "source_sample_used_as_supervised_label": False,
        "event_key_overlap_trainval_locked_test": event_overlap_from_ids(matrix, matrix.trainval_mask, test_mask),
        "stack_score_missing_bdt_rows": int(np.isnan(oof_bdt[stack_score_mask]).sum()),
        "stack_score_missing_mlp_rows": int(np.isnan(oof_mlp[stack_score_mask]).sum()),
        "locked_test_missing_bdt_rows": int(np.isnan(final_bdt[test_mask]).sum()),
        "locked_test_missing_mlp_rows": int(np.isnan(final_mlp[test_mask]).sum()),
        "locked_test_missing_stack_rows_by_model": stack_missing_by_model,
        "score_columns_finite": bool(
            np.isfinite(oof_bdt[stack_score_mask]).all()
            and np.isfinite(oof_mlp[stack_score_mask]).all()
            and np.isfinite(final_bdt[test_mask]).all()
            and np.isfinite(final_mlp[test_mask]).all()
            and all(value == 0 for value in stack_missing_by_model.values())
        ),
        "split_counts": {
            "trainval": partition_qa["partitions"]["trainval"],
            "locked_test": partition_qa["partitions"]["locked_test"],
            "stack_score_region": {
                "rows": int(stack_score_mask.sum()),
                "events": int(len(np.unique(matrix.row_event_id[stack_score_mask]))),
                "event_digest": event_digest_from_ids(matrix, stack_score_mask),
                "label_counts": mono.label_counts(y, stack_score_mask),
                "source_counts": mono.source_counts(frame, stack_score_mask),
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
    base_artifacts = collect_base_artifacts(args)
    write_json(args.outdir / "base_model_artifacts.json", base_artifacts)
    feature_metadata = read_json(args.outdir / "feature_contract.json")
    manifest = {
        "schema": mono.SCHEMA,
        "staged_schema": STAGED_SCHEMA,
        "status": "READY",
        "domain": args.domain,
        "created_unix": time.time(),
        "elapsed_seconds": time.time() - start,
        "input_paths": matrix.manifest.get("input_paths", []),
        "from_cache": bool(matrix.manifest.get("from_cache", False)),
        "optional_columns_seen": matrix.manifest.get("optional_columns_seen", []),
        "sample_definitions": {
            "training_signal_sources": signal_sources,
            "training_inclusive_sources": parse_csv_list(args.training_inclusive_sources),
            "overlay_signal_sources": signal_sources,
            "overlay_inclusive_sources": inclusive_sources,
            **overlay_source,
        },
        "training_class_definition": mono.TRAINING_CLASS_DEFINITION,
        "overlay_class_definition": mono.OVERLAY_CLASS_DEFINITION,
        "feature_contract": feature_metadata,
        "partition_seed": int(args.random_seed),
        "locked_test_fraction": float(args.locked_test_fraction),
        "stack_training_mode": args.stack_training_mode,
        "stack_training_contract": stack_contract,
        "oof_fold_count": int(args.folds),
        "stack_validation_fold": int(args.stack_validation_fold),
        "simple_stack_fold": int(args.simple_stack_fold if args.stack_training_mode == "simple_holdout" and args.simple_stack_fold is not None else args.stack_validation_fold if args.stack_training_mode == "simple_holdout" else -1),
        "metric_convention": {
            "primary_auc": "weighted_auc from train_auau_photon_mlp.auc_score with PPG12-exact weights when requested",
            "also_recorded": ["unweighted_auc", f"WP{int(args.target_signal_efficiency*100)} signal-efficiency threshold/fake-rate"],
        },
        "class_definition_checks": {
            "training_uses_source_sample_as_label": False,
            "source_sample_role": "provenance and overlay/sample QA only",
            "source_sample_source": "ROOT branch when present; otherwise deterministic path-derived sample name from train_auau_photon_bdt.infer_source_sample",
        },
        "execution_model": {
            "type": "condor_dag_staged_matrix_base_score_stack_reduce",
            "matrix_manifest": str(args.outdir / "matrix" / "matrix_manifest.json"),
            "single_process_monolith": False,
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
        "weight_report": matrix.manifest.get("weight_report", {}),
        "leakage_qa": leakage_qa,
    }
    write_json(args.outdir / "campaign_manifest.json", manifest)
    done = {"schema": STAGED_SCHEMA, "stage": "reduce-domain", "status": "READY", "domain": args.domain, "manifest": manifest}
    write_json(args.outdir / "stage_manifests" / "reduce_domain.json", done)
    print(json.dumps({"stage": "reduce-domain", "status": "READY", "domain": args.domain}, sort_keys=True))
    return manifest


def add_common_args(ap: argparse.ArgumentParser) -> None:
    ap.add_argument("--stage", required=True, choices=["build-matrix", "train-base", "score-base", "train-stack", "reduce-domain"])
    ap.add_argument("--domain", choices=["pp", "auau"], required=True)
    ap.add_argument("--input", nargs="*", default=[])
    ap.add_argument("--overlay-input", nargs="*", default=[])
    ap.add_argument("--tree", default="DecayPhotonInfo")
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--matrix-dir", type=Path, default=None, help="Shared matrix directory. Defaults to --outdir/matrix.")
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
    ap.add_argument("--stack-training-mode", choices=["oof5", "simple_holdout"], default="oof5")
    ap.add_argument("--stack-validation-fold", type=int, default=0)
    ap.add_argument("--simple-stack-fold", type=int, default=None)
    ap.add_argument("--random-seed", type=int, default=260604)
    ap.add_argument("--target-signal-efficiency", type=float, default=0.80)
    ap.add_argument("--report-et-bins", default="15,17,19,21,23,25,27,30,35")
    ap.add_argument("--report-cent-bins", default="0,20,40,60,80")
    ap.add_argument("--overlay-bins", type=int, default=60)
    ap.add_argument("--n-jobs", type=int, default=2)
    ap.add_argument("--predict-chunk-rows", type=int, default=100000)
    ap.add_argument("--model-kind", choices=["bdt", "mlp"], default="bdt")
    ap.add_argument("--model-role", choices=["fold", "simple", "final"], default="fold")
    ap.add_argument("--model-name", default="")
    ap.add_argument("--score-region", choices=["fold", "simple_holdout", "locked_test"], default="fold")
    ap.add_argument("--fold", type=int, default=0)
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
    ap.add_argument("--allow-bdt-linear-fallback", action="store_true")
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
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--self-test-rows", type=int, default=5000)


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    add_common_args(ap)
    args = ap.parse_args()
    if args.stack_validation_fold < 0 or args.stack_validation_fold >= args.folds:
        raise SystemExit("--stack-validation-fold must be in [0, folds)")
    if args.simple_stack_fold is not None and (args.simple_stack_fold < 0 or args.simple_stack_fold >= args.folds):
        raise SystemExit("--simple-stack-fold must be in [0, folds)")
    if args.domain == "pp" and not args.overlay_signal_sources:
        args.overlay_signal_sources = ",".join(mono.PP_SIGNAL_SOURCES)
    if args.domain == "pp" and not args.training_inclusive_sources:
        args.training_inclusive_sources = ",".join(mono.PP_TRAIN_INCLUSIVE_SOURCES)
    if args.domain == "pp" and not args.overlay_inclusive_sources:
        args.overlay_inclusive_sources = ",".join(mono.PP_OVERLAY_INCLUSIVE_SOURCES)
    if args.domain == "auau" and not args.overlay_signal_sources:
        args.overlay_signal_sources = ",".join(mono.AUAU_SIGNAL_SOURCES)
    if args.domain == "auau" and not args.training_inclusive_sources:
        args.training_inclusive_sources = ",".join(mono.AUAU_INCLUSIVE_SOURCES)
    if args.domain == "auau" and not args.overlay_inclusive_sources:
        args.overlay_inclusive_sources = ",".join(mono.AUAU_INCLUSIVE_SOURCES)
    if not args.ppg12_exact_expected_samples:
        args.ppg12_exact_expected_samples = ",".join(
            parse_csv_list(args.overlay_signal_sources) + parse_csv_list(args.training_inclusive_sources)
        )
    if args.ppg12_exact_closure_dir is None:
        args.ppg12_exact_closure_dir = args.outdir / "ppg12_exact_reweight_closure"
    return args


def main() -> int:
    args = parse_args()
    dispatch = {
        "build-matrix": stage_build_matrix,
        "train-base": stage_train_base,
        "score-base": stage_score_base,
        "train-stack": stage_train_stack,
        "reduce-domain": stage_reduce_domain,
    }
    dispatch[args.stage](args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
