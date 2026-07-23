#!/usr/bin/env python3
"""Reconstruct and score an exact THE-134 view-specific event-group holdout."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve()
_CONTRACTS = _HERE.parents[1] / "contracts"
sys.path.insert(0, str(_CONTRACTS))
from the134_h70_contract import (  # noqa: E402
    ALL_SHOWER_VIEWS,
    CONTROL_VIEW_BY_SYSTEM,
    FEATURES_BY_SYSTEM,
    MODEL_DOMAIN_GEV,
    VIEW_NAME,
    extraction_authority_binding,
    expected_model_split,
    full_extraction_authority_checks,
    model_artifact_receipt_checks,
    model_extraction_authority_checks,
    sha256_file,
    shower_semantic_sha256,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--system", choices=sorted(FEATURES_BY_SYSTEM), required=True)
    parser.add_argument("--view", choices=ALL_SHOWER_VIEWS, default=VIEW_NAME)
    parser.add_argument("--control-view", choices=ALL_SHOWER_VIEWS, required=True)
    parser.add_argument("--weighted-matrix", type=Path, required=True)
    parser.add_argument("--audited-matrix", type=Path, required=True)
    parser.add_argument("--extraction-audit", type=Path, required=True)
    parser.add_argument("--model-xgb", type=Path, required=True)
    parser.add_argument("--model-metadata", type=Path, required=True)
    parser.add_argument("--model-receipt", type=Path, required=True)
    parser.add_argument("--holdout-out", type=Path, required=True)
    parser.add_argument("--certificate-out", type=Path, required=True)
    return parser.parse_args()


def stable_seed(*items: object) -> int:
    text = "|".join(str(item) for item in items)
    return int(hashlib.sha256(text.encode("utf-8")).hexdigest()[:8], 16)


def event50_mask(keys: np.ndarray, *, seed: int, model_id: str, test_fraction: float) -> tuple[np.ndarray, dict]:
    keys = np.asarray(keys).astype(str)
    unique = np.unique(keys)
    if len(unique) < 4:
        raise SystemExit(f"event50 needs at least four unique events; got {len(unique)}")
    hashes = np.asarray(
        [stable_seed("event50", seed, model_id, key) for key in unique], dtype=np.uint64
    )
    order = np.argsort(hashes, kind="mergesort")
    n_test = int(round(len(unique) * test_fraction))
    n_test = min(max(n_test, 1), len(unique) - 1)
    test_keys = set(unique[order[:n_test]].tolist())
    mask = np.asarray([key in test_keys for key in keys], dtype=bool)
    overlap = len(set(keys[~mask].tolist()).intersection(keys[mask].tolist()))
    return mask, {
        "mode": "event50",
        "random_seed": seed,
        "model_id": model_id,
        "test_fraction_requested": test_fraction,
        "unique_events": int(len(unique)),
        "train_events": int(len(unique) - n_test),
        "test_events": int(n_test),
        "train_rows": int(np.sum(~mask)),
        "test_rows": int(np.sum(mask)),
        "event_overlap": int(overlap),
    }


def aligned_view_matrix(
    audited,
    aligned: np.ndarray,
    features: tuple[str, ...],
    *,
    view_name: str | None,
) -> np.ndarray:
    columns = [
        np.asarray(
            audited[feature if view_name is None else f"view_{view_name}__{feature}"][aligned],
            dtype=np.float32,
        )
        for feature in features
    ]
    return np.ascontiguousarray(np.column_stack(columns), dtype=np.float32)


def main() -> int:
    args = parse_args()
    if args.control_view != CONTROL_VIEW_BY_SYSTEM[args.system]:
        raise SystemExit(
            f"control view for {args.system} is frozen to "
            f"{CONTROL_VIEW_BY_SYSTEM[args.system]}, not {args.control_view}"
        )
    for path in (
        args.weighted_matrix,
        args.audited_matrix,
        args.extraction_audit,
        args.model_xgb,
        args.model_metadata,
        args.model_receipt,
    ):
        if not path.is_file():
            raise SystemExit(f"missing required input: {path}")
    metadata = json.loads(args.model_metadata.read_text())
    model_receipt = json.loads(args.model_receipt.read_text())
    extraction = json.loads(args.extraction_audit.read_text())
    view_contract = metadata.get("the134_view_contract", {})
    expected_split = expected_model_split(args.system, args.view)
    if expected_split["mode"] != "event50":
        raise SystemExit(
            f"{args.system}/{args.view} is an immutable legacy "
            f"{expected_split['mode']} split control; this event-group holdout tool must "
            "not reconstruct or rebaseline its accepted holdout"
        )
    checks = {
        "features": metadata.get("features") == list(FEATURES_BY_SYSTEM[args.system]),
        "pt_range": metadata.get("pt_range") == list(MODEL_DOMAIN_GEV),
        "weight_mode": metadata.get("weight_mode") == "ppg12-exact",
        "split.mode": metadata.get("split", {}).get("mode")
        == expected_split["mode"],
        "split.seed": metadata.get("split", {}).get("random_seed")
        == expected_split["random_seed"],
        "view": view_contract.get("shower_definition") == args.view,
        "semantic": view_contract.get("shower_semantic_sha256")
        == shower_semantic_sha256(args.view),
        "extraction_status": extraction.get("status") == "PASS",
        "extraction_system": extraction.get("system") == args.system,
        "extraction_view": extraction.get("shower_definition") == args.view,
        "extraction_semantic": extraction.get("shower_semantic_sha256")
        == shower_semantic_sha256(args.view),
        "audited_matrix_sha256": extraction.get("matrix_sha256")
        == sha256_file(args.audited_matrix),
        "extraction_audit_sha256": view_contract.get("extraction_audit_sha256")
        == sha256_file(args.extraction_audit),
        "full_extraction_authority": all(
            full_extraction_authority_checks(
                extraction,
                args.system,
                args.view,
                matrix_sha256=sha256_file(args.audited_matrix),
            ).values()
        ),
        "model_extraction_authority": all(
            model_extraction_authority_checks(metadata).values()
        ),
        "extraction_authority_binding": view_contract.get(
            "extraction_authority"
        )
        == extraction_authority_binding(extraction),
        "model_artifact_receipt": all(
            model_artifact_receipt_checks(
                model_receipt,
                args.system,
                args.view,
                model_xgb_sha256=sha256_file(args.model_xgb),
                model_metadata_sha256=sha256_file(args.model_metadata),
                weighted_matrix_sha256=sha256_file(args.weighted_matrix),
            ).values()
        ),
    }
    if not all(checks.values()):
        raise SystemExit(f"model metadata contract failed: {checks}")

    data = np.load(args.weighted_matrix, allow_pickle=True)
    required = set(FEATURES_BY_SYSTEM[args.system]) | {
        "__ppg12_exact_training_weight",
        "is_signal",
        "global_event_key",
        "source_sample",
        "cluster_Et",
        "cluster_Eta",
        "candidate_id_hi",
        "candidate_id_lo",
    }
    missing = sorted(required - set(data.files))
    if missing:
        raise SystemExit(f"weighted matrix missing fields: {missing}")
    audited = np.load(args.audited_matrix, allow_pickle=True)
    audited_required = set(FEATURES_BY_SYSTEM[args.system]) | {
        "candidate_id_hi",
        "candidate_id_lo",
    }
    audited_required |= {
        f"view_{args.control_view}__{feature}"
        for feature in FEATURES_BY_SYSTEM[args.system]
    }
    audited_missing = sorted(audited_required - set(audited.files))
    if audited_missing:
        raise SystemExit(f"audited matrix missing fields: {audited_missing}")
    audited_ids = list(
        zip(
            np.asarray(audited["candidate_id_hi"], dtype=np.uint64).tolist(),
            np.asarray(audited["candidate_id_lo"], dtype=np.uint64).tolist(),
        )
    )
    if len(audited_ids) != len(set(audited_ids)):
        raise SystemExit("audited matrix contains duplicate candidate identities")
    audited_index = {identity: index for index, identity in enumerate(audited_ids)}
    weighted_ids = list(
        zip(
            np.asarray(data["candidate_id_hi"], dtype=np.uint64).tolist(),
            np.asarray(data["candidate_id_lo"], dtype=np.uint64).tolist(),
        )
    )
    if len(weighted_ids) != len(set(weighted_ids)):
        raise SystemExit("weighted matrix contains duplicate candidate identities")
    missing_identities = [identity for identity in weighted_ids if identity not in audited_index]
    if missing_identities:
        raise SystemExit(
            f"weighted matrix has {len(missing_identities)} candidates absent from audited matrix"
        )
    aligned = np.asarray([audited_index[identity] for identity in weighted_ids], dtype=np.int64)
    audited_selected = aligned_view_matrix(
        audited,
        aligned,
        FEATURES_BY_SYSTEM[args.system],
        view_name=None,
    )
    model_id = str(metadata.get("model_id"))
    seed = int(metadata.get("split", {}).get("random_seed"))
    test_fraction = float(metadata.get("split", {}).get("test_fraction_requested"))
    test, split = event50_mask(
        data["global_event_key"], seed=seed, model_id=model_id, test_fraction=test_fraction
    )
    matrix = np.ascontiguousarray(
        np.column_stack(
            [np.asarray(data[feature][test], dtype=np.float32) for feature in FEATURES_BY_SYSTEM[args.system]]
        ),
        dtype=np.float32,
    )
    if not np.array_equal(matrix, audited_selected[test]):
        raise SystemExit(f"weighted/audited {args.view} feature alignment mismatch")
    control_matrix_all = aligned_view_matrix(
        audited,
        aligned,
        FEATURES_BY_SYSTEM[args.system],
        view_name=args.control_view,
    )
    control_view = args.control_view
    control_matrix = control_matrix_all[test]
    try:
        import xgboost as xgb
    except ImportError as exc:
        raise SystemExit("materialize_the134_h70_holdout.py requires xgboost") from exc
    booster = xgb.Booster()
    booster.load_model(str(args.model_xgb))
    scores = np.asarray(booster.predict(xgb.DMatrix(matrix)), dtype=np.float32)
    labels = np.asarray(data["is_signal"][test], dtype=np.int8)
    weights = np.asarray(data["__ppg12_exact_training_weight"][test], dtype=np.float64)
    if set(labels.tolist()) != {0, 1}:
        raise SystemExit(f"holdout does not contain both classes: {sorted(set(labels.tolist()))}")
    if (
        not np.isfinite(matrix).all()
        or not np.isfinite(control_matrix).all()
        or not np.isfinite(scores).all()
    ):
        raise SystemExit("holdout contains nonfinite model inputs/scores")
    if not np.all(np.isfinite(weights) & (weights > 0.0)):
        raise SystemExit("holdout contains nonfinite/non-positive weights")
    if split["event_overlap"] != 0:
        raise SystemExit(f"event overlap is nonzero: {split['event_overlap']}")

    output = {
        "features": np.asarray(FEATURES_BY_SYSTEM[args.system], dtype=object),
        "x": matrix,
        "x_control": control_matrix,
        "control_features": np.asarray(FEATURES_BY_SYSTEM[args.system], dtype=object),
        "control_shower_definition": np.asarray(control_view),
        "control_shower_semantic_sha256": np.asarray(
            shower_semantic_sha256(control_view)
        ),
        "is_signal": labels,
        "training_weight": weights,
        "score_xgboost": scores,
    }
    for name in (
        "cluster_Et",
        "cluster_Eta",
        "centrality",
        "source_sample",
        "input_file_index",
        "input_tree_entry",
        "run",
        "evt",
        "global_event_key",
        "candidate_id_hi",
        "candidate_id_lo",
        "training_label",
        "nominal_is_signal",
    ):
        if name in data.files:
            output[name] = np.asarray(data[name][test])
    args.holdout_out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.holdout_out, **output)
    reported_split = metadata.get("split", {})
    split_count_checks = {
        name: reported_split.get(name) == split[name]
        for name in ("train_events", "test_events", "train_rows", "test_rows")
    }
    certificate = {
        "schema": "THE134_FACTORIAL_VIEW_EXACT_HOLDOUT_CERTIFICATE_V1",
        "status": "PASS" if all(split_count_checks.values()) else "FAIL",
        "system": args.system,
        "shower_definition": args.view,
        "shower_semantic_sha256": shower_semantic_sha256(args.view),
        "feature_order": list(FEATURES_BY_SYSTEM[args.system]),
        "model_origin": metadata.get("model_origin"),
        "reuse_pinned_hashes": metadata.get("reuse_pinned_hashes"),
        "control_shower_definition": control_view,
        "control_shower_semantic_sha256": shower_semantic_sha256(control_view),
        "split": split,
        "reported_split_count_checks": split_count_checks,
        "extraction_authority": extraction_authority_binding(extraction),
        "class_counts": {str(label): int(np.sum(labels == label)) for label in (0, 1)},
        "provenance": {
            "weighted_matrix": str(args.weighted_matrix),
            "weighted_matrix_sha256": sha256_file(args.weighted_matrix),
            "audited_matrix": str(args.audited_matrix),
            "audited_matrix_sha256": sha256_file(args.audited_matrix),
            "extraction_audit": str(args.extraction_audit),
            "extraction_audit_sha256": sha256_file(args.extraction_audit),
            "model_xgb": str(args.model_xgb),
            "model_xgb_sha256": sha256_file(args.model_xgb),
            "model_metadata": str(args.model_metadata),
            "model_metadata_sha256": sha256_file(args.model_metadata),
            "model_receipt": str(args.model_receipt),
            "model_receipt_sha256": sha256_file(args.model_receipt),
            "holdout": str(args.holdout_out),
            "holdout_sha256": sha256_file(args.holdout_out),
        },
    }
    args.certificate_out.parent.mkdir(parents=True, exist_ok=True)
    args.certificate_out.write_text(json.dumps(certificate, indent=2, sort_keys=True) + "\n")
    print(args.holdout_out)
    print(args.certificate_out)
    return 0 if certificate["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
