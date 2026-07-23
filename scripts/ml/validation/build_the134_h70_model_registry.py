#!/usr/bin/env python3
"""Build one immutable paired p+p/Au+Au THE-134 view registry."""

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
    VIEW_NAME,
    canonical_json_sha256,
    expected_model_origin,
    extraction_authority_binding,
    full_extraction_authority_checks,
    is_sha256,
    sha256_file,
    shower_semantic_sha256,
    valid_reuse_pinned_hashes,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--view", choices=ALL_SHOWER_VIEWS, default=VIEW_NAME)
    for system in ("pp", "auau"):
        parser.add_argument(f"--{system}-extraction", type=Path, required=True)
        parser.add_argument(f"--{system}-validation", type=Path, required=True)
        parser.add_argument(f"--{system}-working-points", type=Path, required=True)
        parser.add_argument(
            f"--{system}-model-origin",
            choices=("TRAINED_THE134", "REUSED_THE116", "REUSED_THE111"),
            required=True,
        )
    parser.add_argument("--json-out", type=Path, required=True)
    parser.add_argument("--public-commit", required=True)
    parser.add_argument("--code-sha256", required=True)
    return parser.parse_args()


def read_json(path: Path) -> dict:
    if not path.is_file():
        raise SystemExit(f"missing JSON input: {path}")
    try:
        return json.loads(path.read_text())
    except json.JSONDecodeError as exc:
        raise SystemExit(f"malformed JSON {path}: {exc}") from exc


def require_sha256(text: str, label: str) -> None:
    if len(text) != 64 or any(char not in "0123456789abcdef" for char in text.lower()):
        raise SystemExit(f"{label} is not a 64-character SHA-256: {text!r}")


def verify_artifact_binding(
    provenance: dict,
    path_key: str,
    sha_key: str,
    label: str,
    failures: list[str],
) -> None:
    """Rehash a registered artifact at the final registry authority boundary."""

    expected_sha = provenance.get(sha_key)
    path = Path(str(provenance.get(path_key, "")))
    if not is_sha256(expected_sha):
        failures.append(f"{label} recorded SHA-256 is malformed")
        return
    if not path.is_file() or path.stat().st_size <= 0:
        failures.append(f"{label} artifact is missing or empty: {path}")
        return
    if sha256_file(path) != expected_sha:
        failures.append(f"{label} artifact SHA-256 changed after certification")


def required_origin(system: str, view_name: str) -> str:
    return expected_model_origin(system, view_name)


def validate_lane(
    system: str,
    view_name: str,
    model_origin: str,
    extraction_path: Path,
    validation_path: Path,
    wp_path: Path,
) -> dict:
    extraction = read_json(extraction_path)
    validation = read_json(validation_path)
    wp = read_json(wp_path)
    expected_semantic = shower_semantic_sha256(view_name)
    failures = []
    if model_origin != required_origin(system, view_name):
        failures.append(
            f"model_origin={model_origin!r} expected={required_origin(system, view_name)!r}"
        )
    expected_schemas = {
        "extraction": "THE134_FACTORIAL_VIEW_TRAINING_MATRIX_AUDIT_V1",
        "validation": "THE134_FACTORIAL_VIEW_MODEL_VALIDATION_V1",
        "working_points": "THE134_FACTORIAL_VIEW_WEIGHTED_WORKING_POINTS_V1",
    }
    for label, payload in (("extraction", extraction), ("validation", validation), ("working_points", wp)):
        if payload.get("schema") != expected_schemas[label]:
            failures.append(f"{label}.schema={payload.get('schema')!r}")
        if payload.get("status") != "PASS":
            failures.append(f"{label}.status={payload.get('status')!r}")
        if payload.get("system") != system:
            failures.append(f"{label}.system={payload.get('system')!r}")
        if payload.get("shower_definition") != view_name:
            failures.append(f"{label}.shower_definition={payload.get('shower_definition')!r}")
        if payload.get("shower_semantic_sha256") != expected_semantic:
            failures.append(f"{label}.shower_semantic_sha256 mismatch")
    if extraction.get("feature_order") != list(FEATURES_BY_SYSTEM[system]):
        failures.append("extraction.feature_order mismatch")
    extraction_authority_checks = full_extraction_authority_checks(
        extraction, system, view_name
    )
    if not all(extraction_authority_checks.values()):
        failures.append(
            "extraction full authority failed: "
            + json.dumps(extraction_authority_checks, sort_keys=True)
        )
    if not validation.get("critical_gates") or not all(
        validation.get("critical_gates", {}).values()
    ):
        failures.append("validation critical gates are incomplete or not all true")
    if not wp.get("gates") or not all(wp.get("gates", {}).values()):
        failures.append("working-point gates are incomplete or not all true")
    if validation.get("feature_order") != list(FEATURES_BY_SYSTEM[system]):
        failures.append("validation.feature_order mismatch")
    if wp.get("feature_order") != list(FEATURES_BY_SYSTEM[system]):
        failures.append("working_points.feature_order mismatch")
    if validation.get("model_domain_gev") != list(MODEL_DOMAIN_GEV):
        failures.append("validation.model_domain_gev mismatch")
    if wp.get("model_domain_gev") != list(MODEL_DOMAIN_GEV):
        failures.append("working_points.model_domain_gev mismatch")
    if validation.get("model_origin") != model_origin:
        failures.append("validation.model_origin mismatch")
    if wp.get("model_origin") != model_origin:
        failures.append("working_points.model_origin mismatch")
    reuse_pins = validation.get("reuse_pinned_hashes") or {}
    if (wp.get("reuse_pinned_hashes") or {}) != reuse_pins:
        failures.append("validation/working_points reuse pins mismatch")
    if model_origin.startswith("REUSED_"):
        if not valid_reuse_pinned_hashes(reuse_pins, system, view_name):
            failures.append("validation reused-model pins incomplete")
    elif reuse_pins:
        failures.append("trained model unexpectedly carries reuse pins")
    validation_provenance = validation.get("provenance", {})
    wp_provenance = wp.get("provenance", {})
    for path_key, sha_key, label in (
        ("model_xgb", "model_xgb_sha256", "model XGBoost"),
        ("model_tmva", "model_tmva_sha256", "model TMVA"),
        ("model_metadata", "model_metadata_sha256", "model metadata"),
        ("model_receipt", "model_receipt_sha256", "model completion receipt"),
        ("holdout", "holdout_sha256", "validation holdout"),
        (
            "holdout_certificate",
            "holdout_certificate_sha256",
            "holdout certificate",
        ),
    ):
        verify_artifact_binding(
            validation_provenance, path_key, sha_key, label, failures
        )
    for path_key, sha_key, label in (
        ("score_sample", "score_sample_sha256", "working-point score sample"),
        (
            "population_certificate",
            "population_certificate_sha256",
            "working-point population certificate",
        ),
    ):
        verify_artifact_binding(wp_provenance, path_key, sha_key, label, failures)
    if system == "pp" and wp_provenance.get("score_sample_sha256") != validation_provenance.get(
        "holdout_sha256"
    ):
        failures.append("pp validation/WP holdout SHA mismatch")
    if wp_provenance.get("model_metadata_sha256") != validation_provenance.get(
        "model_metadata_sha256"
    ):
        failures.append("validation/WP model metadata SHA mismatch")
    if wp_provenance.get("model_xgb_sha256") != validation_provenance.get(
        "model_xgb_sha256"
    ):
        failures.append("validation/WP XGBoost model SHA mismatch")
    if (
        not is_sha256(validation_provenance.get("model_receipt_sha256", ""))
        or wp_provenance.get("model_receipt_sha256")
        != validation_provenance.get("model_receipt_sha256")
    ):
        failures.append("validation/WP model completion receipt SHA mismatch")
    expected_extraction_binding = extraction_authority_binding(extraction)
    if validation.get("extraction_authority") != expected_extraction_binding:
        failures.append("validation/full extraction authority mismatch")
    if wp.get("extraction_authority") != expected_extraction_binding:
        failures.append("working-points/full extraction authority mismatch")
    if validation_provenance.get("extraction_audit_sha256") != sha256_file(
        extraction_path
    ):
        failures.append("validation/extraction audit SHA mismatch")
    if system == "pp" and wp_provenance.get(
        "population_certificate_sha256"
    ) != validation_provenance.get("holdout_certificate_sha256"):
        failures.append("pp validation/WP holdout certificate SHA mismatch")
    if failures:
        raise SystemExit(f"{system} registry lane failed: {failures}")
    return {
        "system": system,
        "status": "READY",
        "model_origin": model_origin,
        "reuse_pinned_hashes": reuse_pins or None,
        "shower_definition": view_name,
        "shower_semantic_sha256": expected_semantic,
        "feature_order": list(FEATURES_BY_SYSTEM[system]),
        "model_domain_gev": list(MODEL_DOMAIN_GEV),
        "model": {
            "xgboost": validation_provenance["model_xgb"],
            "xgboost_sha256": validation_provenance["model_xgb_sha256"],
            "tmva": validation_provenance["model_tmva"],
            "tmva_sha256": validation_provenance["model_tmva_sha256"],
            "metadata": validation_provenance["model_metadata"],
            "metadata_sha256": validation_provenance["model_metadata_sha256"],
        },
        "exact_working_points": wp.get("exact_thresholds"),
        "working_point_axis": wp.get("axis"),
        "working_point_bin_edges": wp.get("bin_edges"),
        "threshold_authority": wp.get("threshold_authority"),
        "runtime_surfaces": wp.get("runtime_surfaces"),
        "certificates": {
            "extraction": str(extraction_path),
            "extraction_sha256": sha256_file(extraction_path),
            "validation": str(validation_path),
            "validation_sha256": sha256_file(validation_path),
            "working_points": str(wp_path),
            "working_points_sha256": sha256_file(wp_path),
            "model_receipt": validation_provenance.get("model_receipt"),
            "model_receipt_sha256": validation_provenance.get(
                "model_receipt_sha256"
            ),
            "holdout": validation_provenance.get("holdout"),
            "holdout_sha256": validation_provenance.get("holdout_sha256"),
            "holdout_certificate": validation_provenance.get("holdout_certificate"),
            "holdout_certificate_sha256": validation_provenance.get(
                "holdout_certificate_sha256"
            ),
            "working_point_score_sample": wp_provenance.get("score_sample"),
            "working_point_score_sample_sha256": wp_provenance.get("score_sample_sha256"),
            "working_point_population_certificate": wp_provenance.get(
                "population_certificate"
            ),
            "working_point_population_certificate_sha256": wp_provenance.get(
                "population_certificate_sha256"
            ),
        },
    }


def main() -> int:
    args = parse_args()
    require_sha256(args.code_sha256, "code_sha256")
    if len(args.public_commit) != 40 or any(
        character not in "0123456789abcdef" for character in args.public_commit.lower()
    ):
        raise SystemExit("public_commit must be an exact 40-character Git commit")
    lanes = []
    for system in ("pp", "auau"):
        lanes.append(
            validate_lane(
                system,
                args.view,
                getattr(args, f"{system}_model_origin"),
                getattr(args, f"{system}_extraction"),
                getattr(args, f"{system}_validation"),
                getattr(args, f"{system}_working_points"),
            )
        )
    payload = {
        "schema": "THE134_FACTORIAL_VIEW_PAIRED_MODEL_REGISTRY_V1",
        "status": "READY",
        "promotion_status": "CANDIDATE_CURRENT_NOT_CANONICAL",
        "model_count": 2,
        "systems": ["pp", "auau"],
        "shower_definition": args.view,
        "shower_semantic_sha256": shower_semantic_sha256(args.view),
        "public_commit": args.public_commit,
        "code_sha256": args.code_sha256,
        "models": lanes,
        "boundaries": [
            "registry readiness does not promote either model to CANONICAL",
            "view-specific direct/writer/cache replay certification remains required",
            "THE-121/THE-122 remain closed until final P5C certification",
        ],
    }
    payload["registry_semantic_sha256"] = canonical_json_sha256(payload)
    args.json_out.parent.mkdir(parents=True, exist_ok=True)
    args.json_out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(args.json_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
