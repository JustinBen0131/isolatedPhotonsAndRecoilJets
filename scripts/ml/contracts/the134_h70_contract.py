#!/usr/bin/env python3
"""Frozen detector-neutral identities for THE-134 factorial-view model pairs.

This module intentionally contains no ROOT, pandas, XGBoost, or sPHENIX
imports.  Extraction, training, validation, and working-point tools import the
same constants so a feature-order or shower-definition drift fails before a
model is trained.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Iterable


TREE_NAME = "RJPhotonTrainingViewV1"
TRAINING_SCHEMA_NAME = "RJ_PHOTON_TRAINING_VIEW_V1"
TRAINING_SCHEMA_VERSION = "1"
VIEW_NAME = "H70"
MODEL_DOMAIN_GEV = (15.0, 35.0)
ETA_ABS_MAX = 0.7
READER_ET_EDGES = (15.0, 17.0, 19.0, 21.0, 23.0, 25.0, 28.0, 32.0, 35.0)
AUAU_CENTRALITY_EDGES = (0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0)
TARGET_EFFICIENCIES = (0.90, 0.80, 0.70)
ALL_SHOWER_VIEWS = ("H70", "H0", "G70", "G0", "O70", "O0", "R70")
CONTROL_VIEW_BY_SYSTEM = {"pp": "H70", "auau": "H0"}
CONTROL_MODEL_ORIGIN_BY_SYSTEM = {"pp": "REUSED_THE116", "auau": "REUSED_THE111"}
CONTROL_TMVA_SHA256_BY_SYSTEM = {
    "pp": "228d4cb73f7dc945a613c5a604add71a372b7540c2dc8c630b533d215bb17b30",
    "auau": "d50c69ec98558cb80730ab45fe6801d4accbcf2221c482af91e8899cede1c925",
}
REUSED_MODEL_ORIGIN = {("pp", "H70"): "REUSED_THE116", ("auau", "H0"): "REUSED_THE111"}
FROZEN_REUSE_PINNED_HASHES = {
    ("pp", "H70"): {
        "source_sha256": "e9852be2112c3d91d9ab248731ee9d199e03213b6dd5c93e6f48bfdaaa206258",
        "trainer_sha256": "06de85c527e1f7d41a48b66046051bad70f238acc4af13bd5520a263a514b2c3",
        "pipeline_sha256": "7395397e70fc834f7534fc44a993c9d55829bdf465de5da8e2fb9c58fe851233",
        "config_sha256": "7a908a954b382e46a17d3fdcd7b8188dc65a68e320c1b5e80557f0d5419087e8",
        "working_points_sha256": "264cf8aaf078fc1c9b6340075772f9ea8809ae78d9a91cb5836122bdf8cd352c",
    },
    ("auau", "H0"): {
        "authority_certificate_sha256": "bccf77b766e9ee652b1ca8d088dd2027612d4d768738d4b5e8857f2abc1bc4b3",
        "model_xgb_sha256": "8d1e07d5b2ef7b692442661ce0e0e8b05c242e374a21686be0eb485f7d4bfb12",
        "model_tmva_sha256": "d50c69ec98558cb80730ab45fe6801d4accbcf2221c482af91e8899cede1c925",
        "model_metadata_sha256": "93d7dd4fbbb22d4601fa6f57f058852aad40e1583a9c076a49e27e3ef111b219",
        "accepted_holdout_sha256": "4663dd99e9dc84812c086b175eb2bc8054194075743512fd2cd63c5e44ea45e4",
        "working_points_sha256": "5cead86e42ced5222b2b82630df8517c66434593f126201616884e3ac925463d",
    },
}
LEGACY_ACCEPTED_REUSE_SPLITS = {
    ("pp", "H70"): {
        "mode": "event50",
        "test_fraction_requested": 0.5,
        "random_seed": 42,
        "boundary": "THE116_ACCEPTED_EVENT_GROUP_HOLDOUT",
    },
    ("auau", "H0"): {
        "mode": "row",
        "test_fraction_requested": 0.1,
        "random_seed": 13,
        "boundary": "THE111_ACCEPTED_LEGACY_CANDIDATE_ROW_HOLDOUT_CONTROL_ONLY",
    },
}

# These are scientific acceptance maxima, not tunable debugging defaults.
MAX_RUNTIME_SCORE_ABS_DIFFERENCE = 2.0e-7
MAX_OVERALL_AUC_REGRESSION = 0.005
MAX_BINNED_AUC_REGRESSION = 0.02
MAX_TRAIN_HOLDOUT_AUC_GAP = 0.12
MAX_WP_EXACT_EFFICIENCY_ERROR = 0.02
MAX_SURFACE_RMS = 0.015
MAX_SURFACE_ABS_RESIDUAL = 0.025
MAX_SURFACE_EFFICIENCY_ERROR = 0.02

PP_FEATURES = (
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
)

AUAU_FEATURES = (
    "cluster_Et",
    "cluster_weta_cogx",
    "cluster_wphi_cogx",
    "cluster_weta33_cogx",
    "cluster_wphi33_cogx",
    "vertexz",
    "cluster_Eta",
    "e11_over_e33",
    "cluster_et1",
    "cluster_et2",
    "cluster_et3",
    "cluster_et4",
    "e32_over_e35",
    "centrality",
)

SYSTEM_CODES = {"pp": 1, "auau": 2}
FEATURES_BY_SYSTEM = {"pp": PP_FEATURES, "auau": AUAU_FEATURES}
SIGNAL_SOURCES = {
    "pp": ("run28_photonjet5", "run28_photonjet10", "run28_photonjet20"),
    "auau": ("run28_embeddedPhoton12", "run28_embeddedPhoton20"),
}
BACKGROUND_SOURCES = {
    "pp": ("run28_jet8", "run28_jet12", "run28_jet20", "run28_jet30"),
    "auau": (
        "run28_embeddedJet12",
        "run28_embeddedJet20",
        "run28_embeddedJet30",
        "run28_embeddedJet40",
    ),
}
DIAGNOSTIC_ONLY_SOURCES = {"pp": ("run28_jet40",), "auau": ()}

CORE_BRANCHES = (
    "training_view_id_hi",
    "training_view_id_lo",
    "source_occurrence_id_hi",
    "source_occurrence_id_lo",
    "event_id_hi",
    "event_id_lo",
    "candidate_id_hi",
    "candidate_id_lo",
    "definition_id_hi",
    "definition_id_lo",
    "definition_name",
    "shower_semantic_sha256",
    "feature_contract_sha256",
    "ordered_features",
    "feature_count",
    "finite_feature_state",
    "system_code",
    "source_lane",
    "source_dataset",
    "source_sample",
    "source_period",
    "source_si_di_role",
    "source_ownership_state",
    "source_manifest_sha256",
    "input_uri_sha256",
    "input_file_sha256",
    "source_run",
    "source_segment",
    "run",
    "event_sequence",
    "encounter_ordinal",
    "cluster_index",
    "cluster_map_key",
    "cluster_Et",
    "cluster_Eta",
    "cluster_Phi",
    "vertexz",
    "centrality",
    "model_domain_state",
    "below15_retention_state",
    "nominal_training_eligible",
    "working_point_state",
    "tag_state",
    "training_label",
    "is_signal",
    "label_authority",
    "source_role",
    "source_sample_code",
    "ppg12_source_role_label",
    "minimum_bias_classifier_decision",
    "event_weight",
    "weight_slice",
    "weight_cross_section",
    "weight_vertex",
    "weight_si_di",
    "weight_period",
    "weight_exposure",
    "weight_final",
    "truth_match_found",
    "truth_photon_class",
    "truth_is_prompt",
    "truth_iso_et",
    "truth_iso_pass",
    "cluster_truth_track_id",
    "cluster_truth_pid",
    "cluster_truth_barcode",
    "truth_energy_contribution",
    "npb_label",
    "is_npb",
    "ppg12_analysis_window_pass",
    "ppg12_response_window_pass",
    "ppg12_truth_window_pass_r04",
    "ppg12_sample_bin",
    "ppg12_xsec_pb",
    "ppg12_xsec_weight",
    "ppg12_window_low",
    "ppg12_window_high",
    "max_truth_jet_pt_r04",
    "weight_application_count",
)


def feature_contract_text(system: str) -> str:
    """Mirror ``RJPhotonTrainingViewV1::featureContractText`` exactly."""

    require_system(system)
    return (
        f"RJ_PHOTON_TRAINING_FEATURE_ORDER_V1|system={system}|"
        + ",".join(FEATURES_BY_SYSTEM[system])
    )


def feature_contract_sha256(system: str) -> str:
    return hashlib.sha256(feature_contract_text(system).encode("utf-8")).hexdigest()


def training_schema_text() -> str:
    """Mirror ``RJPhotonTrainingViewV1::schemaSha256`` input exactly."""

    return (
        "RJ_PHOTON_TRAINING_VIEW_V1|tree=RJPhotonTrainingViewV1|"
        "identity=source,event,candidate,definition|"
        "source_stable_inputs=lane,dataset,sample,period,run,segment,input_uri,input_file,manifest|"
        "event_stable_inputs=lane,sample,run,segment,event_sequence|"
        "candidate_stable_inputs=event,encounter_ordinal,cluster_map_key|"
        "features=ordered+contract|"
        "labels=truth+source+npb|weights=component-ledger|domain=15to35"
    )


def training_schema_sha256() -> str:
    return hashlib.sha256(training_schema_text().encode("utf-8")).hexdigest()


def is_sha256(value: object) -> bool:
    text = str(value).lower()
    return len(text) == 64 and all(character in "0123456789abcdef" for character in text)


def identity128_from_text(canonical_input: str) -> tuple[int, int]:
    digest = hashlib.sha256(canonical_input.encode("utf-8")).digest()
    hi = int.from_bytes(digest[:8], byteorder="big", signed=False)
    lo = int.from_bytes(digest[8:16], byteorder="big", signed=False)
    return (hi, lo if hi or lo else 1)


def identity128_hex(identity: tuple[int, int]) -> str:
    return f"{identity[0]:016x}{identity[1]:016x}"


def shower_semantic_text(view_name: str = VIEW_NAME) -> str:
    """Mirror ``RJShowerFactorialV1::semanticText`` exactly."""

    definitions = {
        "H70": (0, 0, 1, "0.070000"),
        "H0": (0, 0, 1, "0.000000"),
        "G70": (0, 0, 0, "0.070000"),
        "G0": (0, 0, 0, "0.000000"),
        "O70": (0, 1, 1, "0.070000"),
        "O0": (0, 1, 1, "0.000000"),
        "R70": (1, 1, 1, "0.070000"),
    }
    if view_name not in definitions:
        raise ValueError(f"unknown shower definition: {view_name}")
    energy, sums, moments, floor = definitions[view_name]
    return (
        f"RJ_SHOWER_DEFINITION_FACTORIAL_V1|{view_name}"
        f"|energy={energy}|sums={sums}|moments={moments}|floor_gev={floor}"
        "|grid=7x7|tower_quality=TowerInfo_get_isGood_only"
        "|center_excluded_from_cogx_numerator=1"
        "|numeric=float32_PhotonClusterBuilder_row_major"
    )


def shower_semantic_sha256(view_name: str = VIEW_NAME) -> str:
    return hashlib.sha256(shower_semantic_text(view_name).encode("utf-8")).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_json_sha256(payload: object) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def source_input_records_sha256(records: Iterable[dict]) -> str:
    """Hash the exact normalized input set for one frozen source population.

    A source token and a nonzero row count are not evidence of a complete
    source population.  Full extraction authority therefore binds the exact
    path/provenance tuple set independently for every source.
    """

    fields = (
        "path",
        "system",
        "source_sample",
        "input_uri_sha256",
        "input_file_sha256",
        "source_manifest_sha256",
        "config_sha256",
        "code_sha256",
        "training_view_root_sha256",
    )
    normalized = [
        {field: str(record.get(field, "")) for field in fields}
        for record in records
    ]
    normalized.sort(key=lambda record: tuple(record[field] for field in fields))
    return canonical_json_sha256(normalized)


def extraction_authority_binding(audit: dict) -> dict[str, object]:
    """Return the immutable full-extraction lineage copied across stages."""

    closure = audit.get("source_population_closure", {})
    return {
        "schema": "THE134_FULL_EXTRACTION_AUTHORITY_BINDING_V1",
        "scope": audit.get("scope"),
        "full_training_authority": audit.get("full_training_authority"),
        "source_provenance_json_sha256": audit.get(
            "source_provenance_json_sha256"
        ),
        "source_population_closure_sha256": closure.get("semantic_sha256"),
    }


def full_extraction_authority_checks(
    audit: dict,
    system: str,
    view_name: str,
    *,
    matrix_sha256: str | None = None,
) -> dict[str, bool]:
    """Fail-closed checks for an authoritative full training extraction."""

    require_system(system)
    closure = audit.get("source_population_closure", {})
    source_rows = closure.get("sources", {})
    required = sorted(expected_sources(system))
    per_source = bool(source_rows) and set(source_rows) == set(required)
    if per_source:
        per_source = all(
            isinstance(source_rows[source], dict)
            and source_rows[source].get("status") == "PASS"
            and all(source_rows[source].get("checks", {}).values())
            and int(source_rows[source].get("observed_input_count", 0)) > 0
            and int(source_rows[source].get("observed_occurrence_count", 0)) > 0
            and is_sha256(source_rows[source].get("full_source_manifest_sha256"))
            and is_sha256(source_rows[source].get("input_records_sha256"))
            for source in required
        )
    closure_hash = closure.get("semantic_sha256")
    closure_hash_valid = (
        is_sha256(closure_hash)
        and closure_hash == canonical_json_sha256(source_rows)
    )
    checks = {
        "schema": audit.get("schema")
        == "THE134_FACTORIAL_VIEW_TRAINING_MATRIX_AUDIT_V1",
        "status": audit.get("status") == "PASS",
        "system": audit.get("system") == system,
        "scope_full": audit.get("scope") == "full",
        "full_training_authority": audit.get("full_training_authority") == 1,
        "shower_definition": audit.get("shower_definition") == view_name,
        "shower_semantic_sha256": audit.get("shower_semantic_sha256")
        == shower_semantic_sha256(view_name),
        "feature_order": audit.get("feature_order")
        == list(FEATURES_BY_SYSTEM[system]),
        "required_sources": audit.get("required_sources") == required,
        "observed_manifest_sources": audit.get("observed_manifest_sources")
        == required,
        "source_provenance_hash": is_sha256(
            audit.get("source_provenance_json_sha256")
        ),
        "source_population_closure": closure.get("status") == "PASS",
        "source_population_closure_hash": closure_hash_valid,
        "per_source_exact_closure": per_source,
    }
    if matrix_sha256 is not None:
        checks["matrix_sha256"] = (
            is_sha256(matrix_sha256)
            and audit.get("matrix_sha256") == matrix_sha256
        )
    return checks


def model_extraction_authority_checks(metadata: dict) -> dict[str, bool]:
    """Require a model receipt to carry an authoritative extraction binding."""

    binding = metadata.get("the134_view_contract", {}).get(
        "extraction_authority", {}
    )
    return {
        "schema": binding.get("schema")
        == "THE134_FULL_EXTRACTION_AUTHORITY_BINDING_V1",
        "scope_full": binding.get("scope") == "full",
        "full_training_authority": binding.get("full_training_authority") == 1,
        "source_provenance_hash": is_sha256(
            binding.get("source_provenance_json_sha256")
        ),
        "source_population_closure_hash": is_sha256(
            binding.get("source_population_closure_sha256")
        ),
    }


def model_artifact_receipt_checks(
    receipt: dict,
    system: str,
    view_name: str,
    *,
    model_xgb_sha256: str,
    model_metadata_sha256: str,
    weighted_matrix_sha256: str | None = None,
) -> dict[str, bool]:
    """Bind downstream scoring to one completed training/reuse receipt."""

    artifacts = receipt.get("artifacts", {})
    expected_schema = (
        "THE134_FACTORIAL_VIEW_MODEL_REUSE_COMPLETE_V1"
        if expected_model_origin(system, view_name).startswith("REUSED_")
        else "THE134_FACTORIAL_VIEW_MODEL_TRAINING_COMPLETE_V1"
    )
    checks = {
        "schema": receipt.get("schema") == expected_schema,
        "status": receipt.get("status") == "PASS",
        "system": receipt.get("system") == system,
        "shower_definition": receipt.get("shower_definition") == view_name,
        "shower_semantic_sha256": receipt.get("shower_semantic_sha256")
        == shower_semantic_sha256(view_name),
        "model_origin": receipt.get("model_origin")
        == expected_model_origin(system, view_name),
        "model_xgb_sha256": artifacts.get("xgboost_sha256")
        == model_xgb_sha256,
        "model_metadata_sha256": artifacts.get("metadata_sha256")
        == model_metadata_sha256,
        "extraction_authority": receipt.get("extraction_authority", {}).get(
            "full_training_authority"
        )
        == 1,
    }
    if weighted_matrix_sha256 is not None:
        checks["weighted_matrix_sha256"] = artifacts.get(
            "weighted_matrix_sha256"
        ) == weighted_matrix_sha256
    return checks


def expected_sources(system: str, *, include_diagnostic: bool = False) -> tuple[str, ...]:
    require_system(system)
    sources = SIGNAL_SOURCES[system] + BACKGROUND_SOURCES[system]
    if include_diagnostic:
        sources += DIAGNOSTIC_ONLY_SOURCES[system]
    return sources


def require_system(system: str) -> None:
    if system not in SYSTEM_CODES:
        raise ValueError(f"unknown system {system!r}; expected one of {sorted(SYSTEM_CODES)}")


def expected_model_origin(system: str, view_name: str) -> str:
    require_system(system)
    if view_name not in ALL_SHOWER_VIEWS:
        raise ValueError(f"unknown shower definition {view_name!r}")
    return REUSED_MODEL_ORIGIN.get((system, view_name), "TRAINED_THE134")


def expected_model_split(system: str, view_name: str) -> dict[str, object]:
    """Return the immutable accepted split contract for one system/view model.

    THE-111 H0 remains a legacy control with a candidate-row holdout. New
    Au+Au view training, including H70, remains event-grouped 90/10.
    """

    origin = expected_model_origin(system, view_name)
    if (system, view_name) in LEGACY_ACCEPTED_REUSE_SPLITS:
        return dict(LEGACY_ACCEPTED_REUSE_SPLITS[(system, view_name)])
    return {
        "mode": "event50",
        "test_fraction_requested": 0.5 if system == "pp" else 0.1,
        "random_seed": 42 if system == "pp" else 13,
        "boundary": f"{origin}_EVENT_GROUP_HOLDOUT",
    }


def expected_reuse_pinned_hashes(system: str, view_name: str) -> dict[str, str]:
    expected_model_origin(system, view_name)
    return dict(FROZEN_REUSE_PINNED_HASHES.get((system, view_name), {}))


def valid_reuse_pinned_hashes(
    pins: object, system: str, view_name: str
) -> bool:
    expected = expected_reuse_pinned_hashes(system, view_name)
    return (
        bool(expected)
        and isinstance(pins, dict)
        and pins == expected
        and all(is_sha256(value) for value in pins.values())
    )


def source_from_path(path: Path, system: str, *, include_diagnostic: bool = True) -> str:
    """Resolve exactly one frozen source token from a path."""

    candidates = [
        source
        for source in expected_sources(system, include_diagnostic=include_diagnostic)
        if source in path.parts or f"/{source}/" in str(path)
    ]
    if len(candidates) != 1:
        raise ValueError(f"cannot resolve exactly one {system} source from {path}: {candidates}")
    return candidates[0]


def expected_label(system: str, source: str) -> int | None:
    require_system(system)
    if source in SIGNAL_SOURCES[system]:
        return 1
    if source in BACKGROUND_SOURCES[system]:
        return 0
    if source in DIAGNOSTIC_ONLY_SOURCES[system]:
        return None
    raise ValueError(f"source {source!r} is outside the frozen {system} authority")


def expected_source_role_and_code(system: str, source: str) -> tuple[int, int]:
    """Return the frozen signal/background role and generator-bin code."""

    label = expected_label(system, source)
    if label not in (0, 1):
        raise ValueError(f"diagnostic source {source!r} has no training role")
    prefix = source.rstrip("0123456789")
    suffix = source[len(prefix) :]
    if not suffix:
        raise ValueError(f"source {source!r} has no numeric sample code")
    return (1 if label == 1 else 2, int(suffix))


def stable_row_digest(parts: Iterable[object]) -> str:
    return hashlib.sha256("|".join(str(part) for part in parts).encode("utf-8")).hexdigest()
