#!/usr/bin/env python3
"""Fail-closed reuse audit for an existing THE-134 factorial-view model.

This audit is mandatory for p+p H70 because THE-116 already owns the frozen
15--35 GeV model.  A new training run is forbidden when the model, feature,
source, trainer, pipeline, configuration, and shower-semantic identities close.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

_HERE = Path(__file__).resolve()
_CONTRACTS = _HERE.parents[1] / "contracts"
sys.path.insert(0, str(_CONTRACTS))
from the134_h70_contract import (  # noqa: E402
    ALL_SHOWER_VIEWS,
    FEATURES_BY_SYSTEM,
    MODEL_DOMAIN_GEV,
    expected_model_split,
    expected_reuse_pinned_hashes,
    is_sha256,
    sha256_file,
    shower_semantic_sha256,
)


PIN_KEYS = (
    "source_sha256",
    "trainer_sha256",
    "pipeline_sha256",
    "config_sha256",
    "working_points_sha256",
)
AUTHORIZED_REUSE = {("pp", "H70"): "REUSED_THE116", ("auau", "H0"): "REUSED_THE111"}
FROZEN_MODEL_HASHES = {
    ("pp", "H70"): {
        "model_xgb_sha256": "6e1ccb0d2c76e10dbdb4ce0dd5aaa7e56bae81eb5b20f3a77ab6bcf22c749727",
        "model_tmva_sha256": "228d4cb73f7dc945a613c5a604add71a372b7540c2dc8c630b533d215bb17b30",
    },
    ("auau", "H0"): {
        "model_xgb_sha256": "8d1e07d5b2ef7b692442661ce0e0e8b05c242e374a21686be0eb485f7d4bfb12",
        "model_tmva_sha256": "d50c69ec98558cb80730ab45fe6801d4accbcf2221c482af91e8899cede1c925",
        "model_metadata_sha256": "93d7dd4fbbb22d4601fa6f57f058852aad40e1583a9c076a49e27e3ef111b219",
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--system", choices=sorted(FEATURES_BY_SYSTEM), required=True)
    parser.add_argument("--view", choices=ALL_SHOWER_VIEWS, required=True)
    parser.add_argument("--model-xgb", type=Path, required=True)
    parser.add_argument("--model-tmva", type=Path, required=True)
    parser.add_argument("--model-metadata", type=Path, required=True)
    parser.add_argument("--authority-certificate", type=Path, required=True)
    parser.add_argument("--json-out", type=Path, required=True)
    return parser.parse_args()


def reuse_split_checks(metadata: dict, system: str, view_name: str) -> dict[str, bool]:
    expected = expected_model_split(system, view_name)
    split = metadata.get("split", {})
    positive_rows = int(split.get("train_rows", 0)) > 0 and int(
        split.get("test_rows", 0)
    ) > 0
    return {
        "split_mode": split.get("mode") == expected["mode"],
        "split_fraction": split.get("test_fraction_requested")
        == expected["test_fraction_requested"],
        "split_seed": split.get("random_seed") == expected["random_seed"],
        "split_population": (
            positive_rows
            and (
                int(split.get("train_events", 0)) > 0
                and int(split.get("test_events", 0)) > 0
                if expected["mode"] == "event50"
                else True
            )
        ),
    }


def main() -> int:
    args = parse_args()
    if (args.system, args.view) not in AUTHORIZED_REUSE:
        raise SystemExit(
            f"no frozen reuse authority for {args.system}/{args.view}; training-required "
            "views cannot be silently reused"
        )
    paths = (
        args.model_xgb,
        args.model_tmva,
        args.model_metadata,
        args.authority_certificate,
    )
    for path in paths:
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"missing/empty reuse input: {path}")
    metadata = json.loads(args.model_metadata.read_text())
    authority = json.loads(args.authority_certificate.read_text())
    reuse_key = (args.system, args.view)
    expected_origin = AUTHORIZED_REUSE[reuse_key]
    expected_pins = expected_reuse_pinned_hashes(args.system, args.view)
    pinned = authority.get("pinned_hashes", {})
    pinned_artifacts = authority.get("pinned_artifacts", {})
    artifact_pin_checks = {}
    if reuse_key == ("pp", "H70"):
        for key in PIN_KEYS:
            artifact_path = Path(str(pinned_artifacts.get(key, "")))
            artifact_pin_checks[key] = (
                artifact_path.is_file()
                and artifact_path.stat().st_size > 0
                and sha256_file(artifact_path) == pinned.get(key)
            )
        authority_identity_checks = {
            "status": authority.get("status") in {"PASS", "READY"},
            "system": authority.get("system") == args.system,
            "shower_definition": authority.get("shower_definition") == args.view,
            "shower_semantic_sha256": authority.get("shower_semantic_sha256")
            == shower_semantic_sha256(args.view),
            "feature_order": authority.get("feature_order")
            == list(FEATURES_BY_SYSTEM[args.system]),
            "model_origin": authority.get("model_origin") == expected_origin,
            "frozen_pinned_hashes": pinned == expected_pins,
            "frozen_pin_artifacts": all(artifact_pin_checks.values()),
        }
    else:
        variant = authority.get("inputs", {}).get("ppg12_variant", {}).get(
            "sha256", {}
        )
        artifact_pin_checks = {
            "authority_certificate_sha256": sha256_file(args.authority_certificate)
            == expected_pins["authority_certificate_sha256"],
            "model_xgb_sha256": sha256_file(args.model_xgb)
            == expected_pins["model_xgb_sha256"],
            "model_tmva_sha256": sha256_file(args.model_tmva)
            == expected_pins["model_tmva_sha256"],
            "model_metadata_sha256": sha256_file(args.model_metadata)
            == expected_pins["model_metadata_sha256"],
            "accepted_holdout_sha256": variant.get("holdout")
            == expected_pins["accepted_holdout_sha256"],
            "working_points_sha256": is_sha256(
                expected_pins["working_points_sha256"]
            ),
        }
        authority_identity_checks = {
            "schema": authority.get("schema")
            == "THE107_AUAU_BDT_LABEL_CONTRACT_VALIDATION_V1",
            "status": authority.get("status")
            == "VALIDATED_SIMULATION_ONLY_NO_PROMOTION",
            "campaign": authority.get("campaign") == "THE-111",
            "strict_production_contract": authority.get(
                "frozen_configuration", {}
            ).get("strict_production_contract")
            is True,
            "frozen_configuration_passed": authority.get(
                "frozen_configuration", {}
            ).get("passed")
            is True,
            "variant_xgb": variant.get("xgboost")
            == expected_pins["model_xgb_sha256"],
            "variant_tmva": variant.get("tmva")
            == expected_pins["model_tmva_sha256"],
            "variant_metadata": variant.get("metadata")
            == expected_pins["model_metadata_sha256"],
            "frozen_pin_artifacts": all(artifact_pin_checks.values()),
        }
        pinned = expected_pins
        pinned_artifacts = {
            "authority_certificate": str(args.authority_certificate),
            "model_xgb": str(args.model_xgb),
            "model_tmva": str(args.model_tmva),
            "model_metadata": str(args.model_metadata),
        }
    frozen_models = FROZEN_MODEL_HASHES[reuse_key]
    authority_checks = {
        **authority_identity_checks,
        "model_xgb_sha256": authority.get("model_xgb_sha256")
        == sha256_file(args.model_xgb)
        if reuse_key == ("pp", "H70")
        else sha256_file(args.model_xgb)
        == frozen_models["model_xgb_sha256"],
        "model_tmva_sha256": authority.get("model_tmva_sha256")
        == sha256_file(args.model_tmva)
        if reuse_key == ("pp", "H70")
        else sha256_file(args.model_tmva)
        == frozen_models["model_tmva_sha256"],
        "model_metadata_sha256": authority.get("model_metadata_sha256")
        == sha256_file(args.model_metadata)
        if reuse_key == ("pp", "H70")
        else sha256_file(args.model_metadata)
        == frozen_models["model_metadata_sha256"],
        "frozen_model_xgb_sha256": frozen_models["model_xgb_sha256"]
        == sha256_file(args.model_xgb),
        "frozen_model_tmva_sha256": frozen_models["model_tmva_sha256"]
        == sha256_file(args.model_tmva),
        "all_frozen_pins_exact": pinned == expected_pins
        and all(is_sha256(value) for value in pinned.values()),
    }
    split_contract = expected_model_split(args.system, args.view)
    metadata_checks = {
        "features": metadata.get("features") == list(FEATURES_BY_SYSTEM[args.system]),
        "pt_range": metadata.get("pt_range") == list(MODEL_DOMAIN_GEV),
        "weight_mode": metadata.get("weight_mode") == "ppg12-exact",
        **reuse_split_checks(metadata, args.system, args.view),
    }
    status = "PASS" if all(authority_checks.values()) and all(metadata_checks.values()) else "FAIL"
    payload = {
        "schema": "THE134_FACTORIAL_VIEW_MODEL_REUSE_AUDIT_V1",
        "status": status,
        "system": args.system,
        "shower_definition": args.view,
        "shower_semantic_sha256": shower_semantic_sha256(args.view),
        "feature_order": list(FEATURES_BY_SYSTEM[args.system]),
        "model_origin": expected_origin,
        "authority_checks": authority_checks,
        "metadata_checks": metadata_checks,
        "accepted_split_contract": split_contract,
        "pinned_hashes": pinned,
        "pinned_artifacts": pinned_artifacts,
        "artifact_pin_checks": artifact_pin_checks,
        "frozen_pin_authority": (
            "embedded accepted predecessor artifact and authority hashes"
        ),
        "artifacts": {
            "xgboost": str(args.model_xgb),
            "xgboost_sha256": sha256_file(args.model_xgb),
            "tmva": str(args.model_tmva),
            "tmva_sha256": sha256_file(args.model_tmva),
            "metadata": str(args.model_metadata),
            "metadata_sha256": sha256_file(args.model_metadata),
            "authority_certificate": str(args.authority_certificate),
            "authority_certificate_sha256": sha256_file(args.authority_certificate),
        },
    }
    args.json_out.parent.mkdir(parents=True, exist_ok=True)
    args.json_out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(args.json_out)
    return 0 if status == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
