#!/usr/bin/env python3
"""Score the full weighted THE-134 population used for Au+Au view WPs.

The Au+Au baseline working points are not holdout quantiles.  They are derived
from the full combined embedded Photon12+20 signal population after the frozen
global class/ET/eta PPG12-exact weight has been computed.  This tool produces a
hash-pinned score sample and a source/label/weight certificate; it does not fit
or train a model.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve()
_CONTRACTS = _HERE.parents[1] / "contracts"
sys.path.insert(0, str(_CONTRACTS))
from the134_h70_contract import (  # noqa: E402
    ALL_SHOWER_VIEWS,
    FEATURES_BY_SYSTEM,
    MODEL_DOMAIN_GEV,
    SIGNAL_SOURCES,
    VIEW_NAME,
    expected_model_origin,
    expected_model_split,
    expected_sources,
    model_artifact_receipt_checks,
    model_extraction_authority_checks,
    sha256_file,
    shower_semantic_sha256,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--system", choices=("auau",), default="auau")
    parser.add_argument("--view", choices=ALL_SHOWER_VIEWS, default=VIEW_NAME)
    parser.add_argument("--weighted-matrix", type=Path, required=True)
    parser.add_argument("--model-xgb", type=Path, required=True)
    parser.add_argument("--model-metadata", type=Path, required=True)
    parser.add_argument("--model-receipt", type=Path, required=True)
    parser.add_argument("--score-sample-out", type=Path, required=True)
    parser.add_argument("--certificate-out", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    for path in (
        args.weighted_matrix,
        args.model_xgb,
        args.model_metadata,
        args.model_receipt,
    ):
        if not path.is_file():
            raise SystemExit(f"missing required input: {path}")
    metadata = json.loads(args.model_metadata.read_text())
    model_receipt = json.loads(args.model_receipt.read_text())
    view = metadata.get("the134_view_contract", {})
    expected_split = expected_model_split(args.system, args.view)
    split = metadata.get("split", {})
    metadata_checks = {
        "features": metadata.get("features") == list(FEATURES_BY_SYSTEM[args.system]),
        "pt_range": metadata.get("pt_range") == list(MODEL_DOMAIN_GEV),
        "weight_mode": metadata.get("weight_mode") == "ppg12-exact",
        "split_mode": metadata.get("split", {}).get("mode")
        == expected_split["mode"],
        "split_fraction": metadata.get("split", {}).get("test_fraction_requested")
        == expected_split["test_fraction_requested"],
        "split_seed": metadata.get("split", {}).get("random_seed")
        == expected_split["random_seed"],
        "split_positive_rows": int(split.get("train_rows", 0)) > 0
        and int(split.get("test_rows", 0)) > 0,
        "split_positive_events": (
            int(split.get("train_events", 0)) > 0
            and int(split.get("test_events", 0)) > 0
            if expected_split["mode"] == "event50"
            else True
        ),
        "accepted_legacy_split_boundary": (
            metadata.get("accepted_split_contract") == expected_split
            if expected_split["mode"] == "row"
            else True
        ),
        "shower_definition": view.get("shower_definition") == args.view,
        "shower_semantic_sha256": view.get("shower_semantic_sha256")
        == shower_semantic_sha256(args.view),
        "old_low_calo_veto_disabled": not bool(
            metadata.get("event_quality_filter", {}).get("enabled", False)
        ),
        "model_origin": metadata.get("model_origin")
        == expected_model_origin(args.system, args.view),
        "full_extraction_authority": all(
            model_extraction_authority_checks(metadata).values()
        ),
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
    weighting = metadata.get("weighting", {})
    metadata_checks.update(
        {
            "event_weight_unused": weighting.get("event_weight_used") is False,
            "vertex_weight_unused": weighting.get("vertex_reweight") is False,
            "centrality_weight_unused": weighting.get("centrality_event_weight") is False,
            "cross_section_weight_unused": weighting.get(
                "cross_section_weight_used_for_training"
            )
            is False,
            "weights_computed_before_binning": weighting.get(
                "weights_computed_before_binning"
            )
            is True,
        }
    )
    if not all(metadata_checks.values()):
        raise SystemExit(
            "model metadata violates full-WP population contract: "
            + json.dumps(metadata_checks, sort_keys=True)
        )

    data = np.load(args.weighted_matrix, allow_pickle=True)
    required = set(FEATURES_BY_SYSTEM[args.system]) | {
        "__ppg12_exact_training_weight",
        "is_signal",
        "source_sample",
        "cluster_Et",
        "cluster_Eta",
        "centrality",
        "candidate_id_hi",
        "candidate_id_lo",
        "minimum_bias_classifier_decision",
        "weight_application_count",
        "label_authority",
    }
    missing = sorted(required - set(data.files))
    if missing:
        raise SystemExit(f"weighted matrix missing fields: {missing}")
    matrix = np.ascontiguousarray(
        np.column_stack(
            [
                np.asarray(data[feature], dtype=np.float32)
                for feature in FEATURES_BY_SYSTEM[args.system]
            ]
        ),
        dtype=np.float32,
    )
    labels = np.asarray(data["is_signal"], dtype=np.int8)
    weights = np.asarray(data["__ppg12_exact_training_weight"], dtype=np.float64)
    sources = np.asarray([str(value) for value in data["source_sample"].tolist()])
    et = np.asarray(data["cluster_Et"], dtype=np.float64)
    eta = np.asarray(data["cluster_Eta"], dtype=np.float64)
    centrality = np.asarray(data["centrality"], dtype=np.float64)
    observed_sources = set(sources.tolist())
    expected = set(expected_sources(args.system))
    source_label = {
        **{source: 1 for source in SIGNAL_SOURCES[args.system]},
        **{source: 0 for source in expected - set(SIGNAL_SOURCES[args.system])},
    }
    label_mismatch = np.asarray(
        [source_label.get(source, -999) != int(label) for source, label in zip(sources, labels)],
        dtype=bool,
    )
    candidate_ids = list(
        zip(
            np.asarray(data["candidate_id_hi"], dtype=np.uint64).tolist(),
            np.asarray(data["candidate_id_lo"], dtype=np.uint64).tolist(),
        )
    )
    population_checks = {
        "source_complete": observed_sources == expected,
        "signal_source_complete": set(sources[labels == 1].tolist())
        == set(SIGNAL_SOURCES[args.system]),
        "source_label_closure": not bool(np.any(label_mismatch)),
        "both_classes": set(labels.tolist()) == {0, 1},
        "finite_model_inputs": bool(np.isfinite(matrix).all()),
        "finite_positive_training_weights": bool(
            np.all(np.isfinite(weights) & (weights > 0.0))
        ),
        "model_domain": bool(
            np.all(
                np.isfinite(et)
                & (et >= MODEL_DOMAIN_GEV[0])
                & (et < MODEL_DOMAIN_GEV[1])
                & np.isfinite(eta)
                & (np.abs(eta) < 0.7)
                & np.isfinite(centrality)
                & (centrality >= 0.0)
                & (centrality < 80.0)
            )
        ),
        "minimum_bias_classifier_pass": bool(
            np.all(np.asarray(data["minimum_bias_classifier_decision"]) == 2)
        ),
        "weight_applied_exactly_once": bool(
            np.all(np.asarray(data["weight_application_count"]) == 1)
        ),
        "label_authority": bool(
            np.all(np.asarray(data["label_authority"]).astype(str) == "PPG12_SOURCE_ROLE")
        ),
        "unique_candidate_identity": len(candidate_ids) == len(set(candidate_ids)),
    }
    if not all(population_checks.values()):
        raise SystemExit(
            "full weighted WP population failed: "
            + json.dumps(population_checks, sort_keys=True)
        )
    try:
        import xgboost as xgb
    except ImportError as exc:
        raise SystemExit("materialize_the134_h70_wp_sample.py requires xgboost") from exc
    booster = xgb.Booster()
    booster.load_model(str(args.model_xgb))
    scores = np.asarray(booster.predict(xgb.DMatrix(matrix)), dtype=np.float32)
    if not np.isfinite(scores).all():
        raise SystemExit("full weighted WP population produced nonfinite scores")

    output = {
        "features": np.asarray(FEATURES_BY_SYSTEM[args.system], dtype=object),
        "x": matrix,
        "is_signal": labels,
        "training_weight": weights,
        "score_xgboost": scores,
    }
    for name in data.files:
        if name != "__columns__" and name not in output:
            output[name] = np.asarray(data[name])
    args.score_sample_out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.score_sample_out, **output)
    certificate = {
        "schema": "THE134_FACTORIAL_VIEW_FULL_WEIGHTED_WP_SAMPLE_CERTIFICATE_V1",
        "status": "PASS",
        "system": args.system,
        "population": "full combined eligible Photon12+20 and Jet12+20+30+40",
        "shower_definition": args.view,
        "shower_semantic_sha256": shower_semantic_sha256(args.view),
        "feature_order": list(FEATURES_BY_SYSTEM[args.system]),
        "model_origin": metadata.get("model_origin"),
        "reuse_pinned_hashes": metadata.get("reuse_pinned_hashes"),
        "extraction_authority": view.get("extraction_authority"),
        "metadata_checks": metadata_checks,
        "population_checks": population_checks,
        "rows": int(len(labels)),
        "signal_rows": int(np.sum(labels == 1)),
        "background_rows": int(np.sum(labels == 0)),
        "observed_sources": sorted(observed_sources),
        "provenance": {
            "weighted_matrix": str(args.weighted_matrix),
            "weighted_matrix_sha256": sha256_file(args.weighted_matrix),
            "model_xgb": str(args.model_xgb),
            "model_xgb_sha256": sha256_file(args.model_xgb),
            "model_metadata": str(args.model_metadata),
            "model_metadata_sha256": sha256_file(args.model_metadata),
            "model_receipt": str(args.model_receipt),
            "model_receipt_sha256": sha256_file(args.model_receipt),
            "score_sample": str(args.score_sample_out),
            "score_sample_sha256": sha256_file(args.score_sample_out),
        },
    }
    args.certificate_out.parent.mkdir(parents=True, exist_ok=True)
    args.certificate_out.write_text(json.dumps(certificate, indent=2, sort_keys=True) + "\n")
    print(args.score_sample_out)
    print(args.certificate_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
