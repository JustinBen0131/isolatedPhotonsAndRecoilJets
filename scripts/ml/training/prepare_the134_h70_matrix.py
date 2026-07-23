#!/usr/bin/env python3
"""Validate THE-134 multi-view extraction and materialize one view matrix.

The input artifact is deliberately separate from the legacy
``AuAuPhotonIDTrainingTree``.  Each accepted legacy training candidate must
have exactly the seven registered shower views in ``RJPhotonTrainingViewV1``.
Only the selected view is projected into the training cache; labels and source roles
are never reconstructed from shower variables.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable

import numpy as np

_HERE = Path(__file__).resolve()
_CONTRACTS = _HERE.parents[1] / "contracts"
sys.path.insert(0, str(_CONTRACTS))
from the134_h70_contract import (  # noqa: E402
    ALL_SHOWER_VIEWS,
    CORE_BRANCHES,
    ETA_ABS_MAX,
    FEATURES_BY_SYSTEM,
    MODEL_DOMAIN_GEV,
    SYSTEM_CODES,
    TRAINING_SCHEMA_NAME,
    TRAINING_SCHEMA_VERSION,
    TREE_NAME,
    VIEW_NAME,
    canonical_json_sha256,
    expected_label,
    expected_source_role_and_code,
    expected_sources,
    feature_contract_sha256,
    identity128_from_text,
    identity128_hex,
    is_sha256,
    sha256_file,
    shower_semantic_text,
    shower_semantic_sha256,
    source_input_records_sha256,
    source_from_path,
    training_schema_sha256,
)


ROOT_METADATA_KEYS = (
    "rj_photon_training_schema",
    "rj_photon_training_schema_version",
    "schema_sha256",
    "pp_feature_contract_sha256",
    "auau_feature_contract_sha256",
    "source_manifest_sha256",
    "config_sha256",
    "code_sha256",
    "rj_photon_training_complete",
    "rj_photon_training_entries",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--system", choices=sorted(SYSTEM_CODES), required=True)
    parser.add_argument("--view", choices=ALL_SHOWER_VIEWS, default=VIEW_NAME)
    parser.add_argument(
        "--input",
        nargs="+",
        type=Path,
        required=True,
        help="ROOT paths or @manifest files; ordering is part of event identity",
    )
    parser.add_argument("--tree", default=TREE_NAME)
    parser.add_argument("--matrix-out", type=Path, required=True)
    parser.add_argument("--audit-out", type=Path, required=True)
    parser.add_argument(
        "--source-provenance-json",
        type=Path,
        help=(
            "External source-manifest closure. Required for full scope; schema "
            "THE134_SOURCE_PROVENANCE_V1 with one exact-path record per input."
        ),
    )
    parser.add_argument(
        "--scope",
        choices=("full", "smoke"),
        default="full",
        help="full requires the exact frozen source set; smoke permits a subset",
    )
    parser.add_argument(
        "--skip-input-hashes",
        action="store_true",
        help="Allowed only for a local smoke; full certification always hashes inputs",
    )
    return parser.parse_args()


def expand_paths(items: Iterable[Path]) -> list[Path]:
    paths: list[Path] = []
    for item in items:
        text = str(item)
        if text.startswith("@"):
            manifest = Path(text[1:])
            if not manifest.is_file():
                raise SystemExit(f"missing manifest: {manifest}")
            for raw in manifest.read_text().splitlines():
                raw = raw.strip()
                if raw and not raw.startswith("#"):
                    paths.append(Path(raw))
        else:
            paths.append(item)
    if not paths:
        raise SystemExit("input expansion produced no ROOT paths")
    duplicates = [str(path) for path, count in Counter(map(str, paths)).items() if count > 1]
    if duplicates:
        raise SystemExit(f"duplicate input paths are forbidden: {duplicates[:10]}")
    return paths


def _text(value: object) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value)


def _identity_pairs(arrays: dict[str, np.ndarray], prefix: str) -> list[tuple[int, int]]:
    return list(
        zip(
            np.asarray(arrays[f"{prefix}_hi"], dtype=np.uint64).tolist(),
            np.asarray(arrays[f"{prefix}_lo"], dtype=np.uint64).tolist(),
        )
    )


def load_source_provenance(path: Path | None) -> tuple[dict[str, dict], dict[str, dict]]:
    if path is None:
        return {}, {}
    if not path.is_file():
        raise SystemExit(f"missing source provenance JSON: {path}")
    payload = json.loads(path.read_text())
    if payload.get("schema") != "THE134_SOURCE_PROVENANCE_V1":
        raise SystemExit("source provenance schema must be THE134_SOURCE_PROVENANCE_V1")
    records = payload.get("inputs")
    if not isinstance(records, list):
        raise SystemExit("source provenance inputs must be a list")
    indexed: dict[str, dict] = {}
    for record in records:
        if not isinstance(record, dict) or not record.get("path"):
            raise SystemExit("source provenance contains a malformed input record")
        key = str(Path(str(record["path"])))
        if key in indexed:
            raise SystemExit(f"duplicate source provenance path: {key}")
        indexed[key] = record
    coverage = payload.get("source_coverage_authority", {})
    if not isinstance(coverage, dict):
        raise SystemExit("source_coverage_authority must be an object")
    return indexed, coverage


def validate_artifact_metadata(
    metadata: dict[str, str | int],
    arrays: dict[str, np.ndarray],
    *,
    external: dict | None,
) -> list[str]:
    failures: list[str] = []
    expected = {
        "rj_photon_training_schema": TRAINING_SCHEMA_NAME,
        "rj_photon_training_schema_version": TRAINING_SCHEMA_VERSION,
        "schema_sha256": training_schema_sha256(),
        "pp_feature_contract_sha256": feature_contract_sha256("pp"),
        "auau_feature_contract_sha256": feature_contract_sha256("auau"),
        "rj_photon_training_complete": "1",
    }
    for name, value in expected.items():
        if str(metadata.get(name, "")) != value:
            failures.append(f"ROOT metadata {name} mismatch")
    for name in ("source_manifest_sha256", "config_sha256", "code_sha256"):
        if not is_sha256(metadata.get(name, "")):
            failures.append(f"ROOT metadata {name} is not SHA-256")
    try:
        declared_entries = int(str(metadata.get("rj_photon_training_entries", "")))
    except ValueError:
        declared_entries = -1
        failures.append("ROOT metadata rj_photon_training_entries is not an integer")
    observed_entries = len(arrays.get("definition_name", []))
    if declared_entries != observed_entries:
        failures.append(
            "ROOT metadata/tree entry-count mismatch "
            f"declared={declared_entries} observed={observed_entries}"
        )
    if int(metadata.get("tree_num_entries", -1)) != observed_entries:
        failures.append("uproot tree entry count disagrees with materialized arrays")
    if external is not None:
        for name in (
            "input_uri_sha256",
            "input_file_sha256",
            "source_manifest_sha256",
            "config_sha256",
            "code_sha256",
        ):
            expected_value = str(external.get(name, ""))
            if not is_sha256(expected_value):
                failures.append(f"external source provenance {name} is not SHA-256")
                continue
            if name in metadata and str(metadata[name]) != expected_value:
                failures.append(f"external/ROOT metadata {name} mismatch")
            if name in arrays:
                row_values = {_text(item) for item in arrays[name]}
                if row_values != {expected_value}:
                    failures.append(f"external/row {name} mismatch")
    return failures


def validate_file_arrays(
    arrays: dict[str, np.ndarray],
    *,
    system: str,
    source: str,
    view_name: str = VIEW_NAME,
) -> tuple[np.ndarray, list[list[float]], dict]:
    """Validate one file and return the selected-view mask plus ordered vectors."""

    n_rows = len(arrays["definition_name"])
    failures: list[str] = []
    feature_count = len(FEATURES_BY_SYSTEM[system])
    expected_code = SYSTEM_CODES[system]
    definitions = np.asarray([_text(item) for item in arrays["definition_name"]])
    semantics = np.asarray([_text(item) for item in arrays["shower_semantic_sha256"]])
    feature_contracts = np.asarray([_text(item) for item in arrays["feature_contract_sha256"]])
    vectors = [list(map(float, row)) for row in arrays["ordered_features"]]
    candidate_ids = _identity_pairs(arrays, "candidate_id")
    training_view_ids = _identity_pairs(arrays, "training_view_id")
    source_ids = _identity_pairs(arrays, "source_occurrence_id")
    event_ids = _identity_pairs(arrays, "event_id")
    definition_ids = _identity_pairs(arrays, "definition_id")

    if any(hi == 0 and lo == 0 for hi, lo in candidate_ids):
        failures.append("null candidate identity")
    if any(hi == 0 and lo == 0 for hi, lo in training_view_ids):
        failures.append("null training-view identity")
    if any(hi == 0 and lo == 0 for hi, lo in source_ids):
        failures.append("null source-occurrence identity")
    if len(set(source_ids)) != 1:
        failures.append(f"source-occurrence identity count per file={len(set(source_ids))}")
    if any(hi == 0 and lo == 0 for hi, lo in event_ids):
        failures.append("null event identity")
    if any(hi == 0 and lo == 0 for hi, lo in definition_ids):
        failures.append("null definition identity")
    duplicate_views = n_rows - len(set(training_view_ids))
    if duplicate_views:
        failures.append(f"duplicate training-view identities={duplicate_views}")

    by_candidate: dict[tuple[int, int], list[str]] = defaultdict(list)
    for candidate_id, definition in zip(candidate_ids, definitions.tolist()):
        by_candidate[candidate_id].append(definition)
    bad_factorials = [
        (candidate_id, observed)
        for candidate_id, observed in by_candidate.items()
        if len(observed) != len(ALL_SHOWER_VIEWS)
        or set(observed) != set(ALL_SHOWER_VIEWS)
        or len(set(observed)) != len(observed)
    ]
    if bad_factorials:
        failures.append(f"candidate factorial closure failures={len(bad_factorials)}")

    encounter_ordinals = np.asarray(arrays["encounter_ordinal"], dtype=np.int64)
    cluster_indices = np.asarray(arrays["cluster_index"], dtype=np.int64)
    if np.any(cluster_indices != encounter_ordinals):
        failures.append(
            "cluster_index/encounter_ordinal label binding mismatch rows="
            f"{int(np.sum(cluster_indices != encounter_ordinals))}"
        )

    # Every definition row for one detector candidate inherits one immutable
    # label/provenance ledger.  Definition-specific features may differ; the
    # source role, truth label, stable candidate inputs, and weights may not.
    candidate_invariant_branches = (
        "source_occurrence_id_hi",
        "source_occurrence_id_lo",
        "event_id_hi",
        "event_id_lo",
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
        "below15_retention_state",
        "training_label",
        "is_signal",
        "label_authority",
        "source_role",
        "source_sample_code",
        "ppg12_source_role_label",
        "minimum_bias_classifier_decision",
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
        "event_weight",
        "weight_slice",
        "weight_cross_section",
        "weight_vertex",
        "weight_si_di",
        "weight_period",
        "weight_exposure",
        "weight_final",
        "weight_application_count",
    )

    def invariant_token(value: object) -> object:
        if isinstance(value, np.generic):
            value = value.item()
        if isinstance(value, (bytes, str)):
            return _text(value)
        if isinstance(value, float) and math.isnan(value):
            return ("float", "nan")
        return value

    candidate_indices: dict[tuple[int, int], list[int]] = defaultdict(list)
    for index, candidate_id in enumerate(candidate_ids):
        candidate_indices[candidate_id].append(index)
    for branch in candidate_invariant_branches:
        values = np.asarray(arrays[branch])
        drifted = 0
        for indices in candidate_indices.values():
            reference = invariant_token(values[indices[0]])
            if any(invariant_token(values[index]) != reference for index in indices[1:]):
                drifted += 1
        if drifted:
            failures.append(
                f"candidate-invariant {branch} drift candidates={drifted}"
            )

    expected_semantics = np.asarray(
        [shower_semantic_sha256(definition) for definition in definitions], dtype=object
    )
    if np.any(semantics != expected_semantics):
        failures.append(
            "shower semantic hash mismatch rows="
            f"{int(np.sum(semantics != expected_semantics))}"
        )
    expected_definition_ids = [
        identity128_from_text(
            f"shower-definition|{definition}|"
            f"{shower_semantic_text(definition)}"
        )
        for definition in definitions.tolist()
    ]
    if definition_ids != expected_definition_ids:
        failures.append(
            "definition identity/name/semantic mismatch rows="
            f"{sum(left != right for left, right in zip(definition_ids, expected_definition_ids))}"
        )
    expected_training_view_ids = [
        identity128_from_text(
            f"{identity128_hex(candidate_id)}|training-view|{identity128_hex(definition_id)}"
        )
        for candidate_id, definition_id in zip(candidate_ids, definition_ids)
    ]
    if training_view_ids != expected_training_view_ids:
        failures.append(
            "training-view identity/candidate/definition mismatch rows="
            f"{sum(left != right for left, right in zip(training_view_ids, expected_training_view_ids))}"
        )
    candidate_event_pairs: dict[tuple[int, int], set[tuple[int, int]]] = defaultdict(set)
    for candidate_id, event_id in zip(candidate_ids, event_ids):
        candidate_event_pairs[candidate_id].add(event_id)
    if any(len(events) != 1 for events in candidate_event_pairs.values()):
        failures.append("candidate identity maps to multiple event identities")
    event_sequence_pairs: dict[int, set[tuple[int, int]]] = defaultdict(set)
    for event_sequence, event_id in zip(
        np.asarray(arrays["event_sequence"], dtype=np.int64).tolist(), event_ids
    ):
        event_sequence_pairs[int(event_sequence)].add(event_id)
    if any(len(events) != 1 for events in event_sequence_pairs.values()):
        failures.append("event_sequence maps to multiple event identities")
    expected_feature_contract = feature_contract_sha256(system)
    if np.any(feature_contracts != expected_feature_contract):
        failures.append(
            "feature-order contract hash mismatch rows="
            f"{int(np.sum(feature_contracts != expected_feature_contract))}"
        )

    observed_feature_counts = np.asarray(arrays["feature_count"], dtype=np.int64)
    if np.any(observed_feature_counts != feature_count):
        failures.append(
            f"feature_count mismatch rows={int(np.sum(observed_feature_counts != feature_count))}"
        )
    vector_length_bad = sum(len(row) != feature_count for row in vectors)
    if vector_length_bad:
        failures.append(f"ordered feature length mismatch rows={vector_length_bad}")
    system_codes = np.asarray(arrays["system_code"], dtype=np.int64)
    if np.any(system_codes != expected_code):
        failures.append(f"wrong system code rows={int(np.sum(system_codes != expected_code))}")
    source_samples = np.asarray([_text(item) for item in arrays["source_sample"]])
    if np.any(source_samples != source):
        failures.append(f"source_sample/path mismatch rows={int(np.sum(source_samples != source))}")
    for branch in (
        "source_lane",
        "source_dataset",
        "source_period",
        "source_si_di_role",
        "source_ownership_state",
    ):
        values = np.asarray([_text(item) for item in arrays[branch]])
        if np.any(values == ""):
            failures.append(f"empty {branch} rows={int(np.sum(values == ''))}")
        if len(set(values.tolist())) != 1:
            failures.append(f"{branch} is not constant within source occurrence")
    for branch in ("source_manifest_sha256", "input_uri_sha256", "input_file_sha256"):
        values = np.asarray([_text(item) for item in arrays[branch]])
        malformed = np.asarray([not is_sha256(value) for value in values], dtype=bool)
        if np.any(malformed):
            failures.append(f"malformed {branch} rows={int(np.sum(malformed))}")
    source_runs = np.asarray(arrays["source_run"], dtype=np.int64)
    source_segments = np.asarray(arrays["source_segment"], dtype=np.int64)
    if len(set(source_runs.tolist())) != 1 or len(set(source_segments.tolist())) != 1:
        failures.append("source run/segment are not constant within source occurrence")
    if n_rows:
        expected_source_id = identity128_from_text(
            "|".join(
                [
                    _text(arrays["source_lane"][0]),
                    _text(arrays["source_dataset"][0]),
                    _text(arrays["source_sample"][0]),
                    _text(arrays["source_period"][0]),
                    str(int(source_runs[0])),
                    str(int(source_segments[0])),
                    _text(arrays["input_uri_sha256"][0]),
                    _text(arrays["input_file_sha256"][0]),
                    _text(arrays["source_manifest_sha256"][0]),
                ]
            )
        )
        if any(identity != expected_source_id for identity in source_ids):
            failures.append("source-occurrence identity/provenance mismatch")
    expected_event_ids = [
        identity128_from_text(
            "|".join(
                [
                    _text(lane),
                    _text(sample),
                    str(int(run)),
                    str(int(segment)),
                    str(int(event_sequence)),
                ]
            )
        )
        for lane, sample, run, segment, event_sequence in zip(
            arrays["source_lane"],
            arrays["source_sample"],
            arrays["run"],
            arrays["source_segment"],
            arrays["event_sequence"],
        )
    ]
    if event_ids != expected_event_ids:
        failures.append(
            "event identity/provenance mismatch rows="
            f"{sum(left != right for left, right in zip(event_ids, expected_event_ids))}"
        )
    expected_candidate_ids = [
        identity128_from_text(
            f"{identity128_hex(event_id)}|candidate|{int(encounter_ordinal)}|"
            f"{int(cluster_map_key)}"
        )
        for event_id, encounter_ordinal, cluster_map_key in zip(
            event_ids, arrays["encounter_ordinal"], arrays["cluster_map_key"]
        )
    ]
    if candidate_ids != expected_candidate_ids:
        failures.append(
            "candidate identity/event/ordinal/cluster-map-key mismatch rows="
            f"{sum(left != right for left, right in zip(candidate_ids, expected_candidate_ids))}"
        )
    label_authority = np.asarray([_text(item) for item in arrays["label_authority"]])
    if np.any(label_authority != "PPG12_SOURCE_ROLE"):
        failures.append(
            "training label authority mismatch rows="
            f"{int(np.sum(label_authority != 'PPG12_SOURCE_ROLE'))}"
        )
    weight_application_count = np.asarray(
        arrays["weight_application_count"], dtype=np.int64
    )
    if np.any(weight_application_count != 1):
        failures.append(
            "weight application count must be exactly one rows="
            f"{int(np.sum(weight_application_count != 1))}"
        )
    weight_names = (
        "weight_slice",
        "weight_cross_section",
        "weight_vertex",
        "weight_si_di",
        "weight_period",
        "weight_exposure",
        "weight_final",
        "event_weight",
    )
    weights = {
        name: np.asarray(arrays[name], dtype=np.float64) for name in weight_names
    }
    finite_positive_weights = np.logical_and.reduce(
        [np.isfinite(values) & (values > 0.0) for values in weights.values()]
    )
    if np.any(~finite_positive_weights):
        failures.append(
            "weight ledger contains nonfinite/nonpositive rows="
            f"{int(np.sum(~finite_positive_weights))}"
        )

    def close(left: np.ndarray, right: np.ndarray) -> np.ndarray:
        scale = np.maximum.reduce(
            [np.ones(len(left)), np.abs(left), np.abs(right)]
        )
        return np.abs(left - right) <= 1.0e-12 * scale

    if system == "pp":
        slice_alias = close(weights["weight_slice"], weights["weight_cross_section"])
        period_alias = close(weights["weight_period"], weights["weight_exposure"])
        if np.any(~slice_alias):
            failures.append(
                "p+p slice/cross-section provenance alias mismatch rows="
                f"{int(np.sum(~slice_alias))}"
            )
        if np.any(~period_alias):
            failures.append(
                "p+p period/exposure provenance alias mismatch rows="
                f"{int(np.sum(~period_alias))}"
            )
        component_product = (
            weights["weight_slice"]
            * weights["weight_vertex"]
            * weights["weight_si_di"]
            * weights["weight_period"]
        )
    else:
        component_product = (
            weights["weight_slice"]
            * weights["weight_cross_section"]
            * weights["weight_vertex"]
            * weights["weight_si_di"]
            * weights["weight_period"]
            * weights["weight_exposure"]
        )
    component_closure = close(component_product, weights["weight_final"])
    if np.any(~component_closure):
        failures.append(
            f"{system} weight-component product/final mismatch rows="
            f"{int(np.sum(~component_closure))}"
        )
    event_final_closure = close(weights["event_weight"], weights["weight_final"])
    if np.any(~event_final_closure):
        failures.append(
            "event_weight/weight_final mismatch rows="
            f"{int(np.sum(~event_final_closure))}"
        )
    minimum_bias = np.asarray(
        arrays["minimum_bias_classifier_decision"], dtype=np.int64
    )
    if system == "auau" and np.any(minimum_bias != 2):
        failures.append(
            "AuAu embedded MinimumBias classifier decision must equal 2 rows="
            f"{int(np.sum(minimum_bias != 2))}"
        )

    labels = np.asarray(arrays["training_label"], dtype=np.int64)
    is_signal = np.asarray(arrays["is_signal"], dtype=np.int64)
    eligible = np.asarray(arrays["nominal_training_eligible"], dtype=np.int64)
    model_domain_state = np.asarray(arrays["model_domain_state"], dtype=np.int64)
    scalar_et = np.asarray(arrays["cluster_Et"], dtype=np.float64)
    scalar_eta = np.asarray(arrays["cluster_Eta"], dtype=np.float64)
    scalar_domain = (
        np.isfinite(scalar_et)
        & np.isfinite(scalar_eta)
        & (scalar_et >= MODEL_DOMAIN_GEV[0])
        & (scalar_et < MODEL_DOMAIN_GEV[1])
        & (np.abs(scalar_eta) < ETA_ABS_MAX)
    )
    if system == "auau":
        scalar_centrality = np.asarray(arrays["centrality"], dtype=np.float64)
        scalar_domain &= (
            np.isfinite(scalar_centrality)
            & (scalar_centrality >= 0.0)
            & (scalar_centrality < 80.0)
        )
    expected_eligible = (
        np.isin(labels, [0, 1]) & (model_domain_state == 0) & scalar_domain
    ).astype(np.int64)
    if np.any(eligible != expected_eligible):
        failures.append(
            "nominal_training_eligible disagrees with training_label rows="
            f"{int(np.sum(eligible != expected_eligible))}"
        )
    label_expected = expected_label(system, source)
    if label_expected is None:
        failures.append(f"diagnostic-only source cannot enter training: {source}")
    elif system == "pp":
        if np.any(labels != label_expected):
            failures.append(f"source-role label mismatch rows={int(np.sum(labels != label_expected))}")
        if np.any(is_signal != labels):
            failures.append(f"training_label/is_signal mismatch rows={int(np.sum(is_signal != labels))}")
    else:
        ppg12_label = np.asarray(arrays["ppg12_source_role_label"], dtype=np.int64)
        if np.any(~np.isin(labels, [-1, 0, 1])):
            failures.append("AuAu training_label contains values outside {-1,0,1}")
        if np.any(labels != ppg12_label):
            failures.append(
                "AuAu training_label/ppg12_source_role_label mismatch rows="
                f"{int(np.sum(labels != ppg12_label))}"
            )
        allowed = np.isin(labels, [-1, label_expected])
        if np.any(~allowed):
            failures.append(
                f"AuAu source-role label mismatch rows={int(np.sum(~allowed))}"
            )
    if label_expected is not None:
        expected_role, expected_code = expected_source_role_and_code(system, source)
        source_roles = np.asarray(arrays["source_role"], dtype=np.int64)
        source_codes = np.asarray(arrays["source_sample_code"], dtype=np.int64)
        if np.any(source_roles != expected_role):
            failures.append(
                "source_role/source authority mismatch rows="
                f"{int(np.sum(source_roles != expected_role))}"
            )
        if np.any(source_codes != expected_code):
            failures.append(
                "source_sample_code/source authority mismatch rows="
                f"{int(np.sum(source_codes != expected_code))}"
            )

    selected_view = definitions == view_name
    expected_semantic = shower_semantic_sha256(view_name)
    if int(np.sum(selected_view)) != len(by_candidate):
        failures.append(
            f"{view_name} population mismatch rows={int(np.sum(selected_view))} "
            f"candidates={len(by_candidate)}"
        )
    if np.any(semantics[selected_view] != expected_semantic):
        failures.append(
            f"{view_name} semantic hash mismatch rows="
            f"{int(np.sum(semantics[selected_view] != expected_semantic))}"
        )
    finite_flags = np.asarray(arrays["finite_feature_state"], dtype=np.int64)
    selected = selected_view & (eligible == 1)
    selected_view_vectors = [vectors[index] for index in np.flatnonzero(selected_view)]
    finite_rows = np.asarray(
        [
            len(row) == feature_count and all(math.isfinite(value) for value in row)
            for row in selected_view_vectors
        ],
        dtype=bool,
    )
    if np.any(finite_flags[selected_view] != finite_rows.astype(np.int64)):
        failures.append(
            f"finite_feature_state disagrees with {view_name} ordered vector rows="
            f"{int(np.sum(finite_flags[selected_view] != finite_rows.astype(np.int64)))}"
        )

    selected_view_indices = np.flatnonzero(selected_view)
    expected_domain_state = np.asarray(
        [
            3 if not finite else (0 if scalar_domain[index] else 1)
            for index, finite in zip(selected_view_indices, finite_rows)
        ],
        dtype=np.int64,
    )
    if np.any(model_domain_state[selected_view] != expected_domain_state):
        failures.append(
            "model_domain_state disagrees with finite 15-35 domain rows="
            f"{int(np.sum(model_domain_state[selected_view] != expected_domain_state))}"
        )
    below15 = selected_view & np.isfinite(scalar_et) & (scalar_et < MODEL_DOMAIN_GEV[0])
    below15_retention = np.asarray(arrays["below15_retention_state"], dtype=np.int64)
    expected_below15_retention = (
        np.isfinite(scalar_et) & (scalar_et < MODEL_DOMAIN_GEV[0])
    ).astype(np.int64)
    if np.any(below15_retention != expected_below15_retention):
        failures.append(
            "below15_retention_state mismatch rows="
            f"{int(np.sum(below15_retention != expected_below15_retention))}"
        )
    working_point_state = np.asarray(arrays["working_point_state"], dtype=np.int64)
    tag_state = np.asarray(arrays["tag_state"], dtype=np.int64)
    if np.any(working_point_state != -1):
        failures.append(
            "training-view working_point_state must be null rows="
            f"{int(np.sum(working_point_state != -1))}"
        )
    if np.any(tag_state != -1):
        failures.append(
            f"training-view tag_state must be null rows={int(np.sum(tag_state != -1))}"
        )
    if np.any(
        below15
        & ((~np.isin(model_domain_state, [1, 3])) | (eligible != 0))
    ):
        failures.append(
            "below-15 row escaped diagnostic-only/null-training state rows="
            f"{int(np.sum(below15 & ((model_domain_state != 1) | (eligible != 0))))}"
        )

    if selected_view_vectors:
        matrix = np.asarray(selected_view_vectors, dtype=np.float64)
        scalar_checks = {
            "cluster_Et": (
                0,
                np.asarray(arrays["cluster_Et"], dtype=np.float64)[selected_view],
            ),
            "cluster_Eta": (
                FEATURES_BY_SYSTEM[system].index("cluster_Eta"),
                np.asarray(arrays["cluster_Eta"], dtype=np.float64)[selected_view],
            ),
            "vertexz": (
                FEATURES_BY_SYSTEM[system].index("vertexz"),
                np.asarray(arrays["vertexz"], dtype=np.float64)[selected_view],
            ),
        }
        if system == "auau":
            scalar_checks["centrality"] = (
                FEATURES_BY_SYSTEM[system].index("centrality"),
                np.asarray(arrays["centrality"], dtype=np.float64)[selected_view],
            )
        for name, (column, scalar) in scalar_checks.items():
            left = np.asarray(matrix[:, column], dtype=np.float32)
            right = np.asarray(scalar, dtype=np.float32)
            mismatch = ~((left == right) | (np.isnan(left) & np.isnan(right)))
            if np.any(mismatch):
                failures.append(f"ordered vector/scalar {name} mismatch rows={int(np.sum(mismatch))}")

    selected_vectors = [vectors[index] for index in np.flatnonzero(selected)]
    return selected, selected_vectors, {
        "status": "PASS" if not failures else "FAIL",
        "rows": n_rows,
        "candidates": len(by_candidate),
        "selected_view": view_name,
        "selected_view_rows": int(np.sum(selected_view)),
        "selected_view_training_rows": int(np.sum(selected)),
        "selected_view_input_invalid_rows": int(np.sum(~finite_rows)),
        "discarded_selected_view_label_minus_one_rows": int(
            np.sum(selected_view & (labels < 0))
        ),
        "below15_selected_view_rows": int(np.sum(below15)),
        "below15_selected_view_eligible_rows": int(np.sum(below15 & (eligible != 0))),
        "out_of_domain_selected_view_rows": int(np.sum(selected_view & ~scalar_domain)),
        "out_of_domain_selected_view_eligible_rows": int(
            np.sum(selected_view & ~scalar_domain & (eligible != 0))
        ),
        "definition_counts": dict(sorted(Counter(definitions.tolist()).items())),
        "source": source,
        "source_occurrence_ids": sorted(
            identity128_hex(identity) for identity in set(source_ids)
        ),
        "expected_label": label_expected,
        "failures": failures,
    }


def read_tree(
    path: Path, tree_name: str
) -> tuple[dict[str, np.ndarray], list[str], dict[str, str | int]]:
    try:
        import awkward as ak
        import uproot
    except ImportError as exc:
        raise SystemExit("prepare_the134_h70_matrix.py requires uproot and awkward") from exc

    with uproot.open(path) as root_file:
        if tree_name not in root_file:
            raise ValueError(f"missing tree {tree_name}")
        tree = root_file[tree_name]
        branch_names = sorted(str(name) for name in tree.keys())
        missing = sorted(set(CORE_BRANCHES) - set(branch_names))
        if missing:
            raise ValueError(f"missing required branches: {missing}")
        awkward_arrays = tree.arrays(list(CORE_BRANCHES), library="ak")
        arrays: dict[str, np.ndarray] = {}
        for name in CORE_BRANCHES:
            values = awkward_arrays[name]
            if name == "ordered_features":
                arrays[name] = np.asarray(ak.to_list(values), dtype=object)
            elif name in {
                "definition_name",
                "shower_semantic_sha256",
                "feature_contract_sha256",
                "source_sample",
                "source_lane",
                "source_dataset",
                "source_period",
                "source_si_di_role",
                "source_ownership_state",
                "source_manifest_sha256",
                "input_uri_sha256",
                "input_file_sha256",
                "label_authority",
            }:
                arrays[name] = np.asarray(ak.to_list(values), dtype=object)
            else:
                arrays[name] = ak.to_numpy(values)
        metadata: dict[str, str | int] = {"tree_num_entries": int(tree.num_entries)}
        for name in ROOT_METADATA_KEYS:
            if name not in root_file:
                metadata[name] = ""
                continue
            metadata[name] = str(root_file[name].member("fTitle"))
    return arrays, branch_names, metadata


def build_source_population_closure(
    *,
    system: str,
    paths: list[Path],
    path_sources: list[str],
    external_provenance: dict[str, dict],
    source_coverage_authority: dict[str, dict],
    file_reports: list[dict],
    selected_rows_by_source: dict[str, int],
) -> dict:
    """Prove exact full-manifest/input/occurrence closure per source."""

    required = sorted(expected_sources(system))
    sources: dict[str, dict] = {}
    for source in required:
        records = [
            external_provenance.get(str(path), {})
            for path, observed_source in zip(paths, path_sources)
            if observed_source == source
        ]
        reports = [
            report
            for report in file_reports
            if report.get("source") == source and report.get("status") == "PASS"
        ]
        occurrence_ids = {
            str(identity)
            for report in reports
            for identity in report.get("source_occurrence_ids", [])
        }
        authority = source_coverage_authority.get(source, {})
        expected_state = (
            "SOURCE_COMPLETE_ZERO_IN_DOMAIN"
            if int(selected_rows_by_source.get(source, 0)) == 0
            else "SOURCE_COMPLETE"
        )
        full_manifest = str(authority.get("full_source_manifest_sha256", ""))
        observed_manifests = {
            str(record.get("source_manifest_sha256", "")) for record in records
        }
        observed_records_sha = source_input_records_sha256(records)
        try:
            expected_input_count = int(authority.get("expected_input_count", 0))
            expected_occurrence_count = int(
                authority.get("expected_occurrence_count", 0)
            )
        except (TypeError, ValueError):
            expected_input_count = 0
            expected_occurrence_count = 0
        checks = {
            "authority_state": authority.get("state") == expected_state,
            "full_source_manifest_sha256": is_sha256(full_manifest),
            "all_inputs_bound_to_full_manifest": observed_manifests
            == {full_manifest},
            "expected_input_count_positive": expected_input_count > 0,
            "input_count_exact": len(records) == expected_input_count,
            "expected_occurrence_count_positive": expected_occurrence_count > 0,
            "occurrence_count_exact": len(occurrence_ids)
            == expected_occurrence_count,
            "input_records_sha256": is_sha256(
                authority.get("input_records_sha256", "")
            ),
            "input_record_set_exact": authority.get("input_records_sha256")
            == observed_records_sha,
            "record_system_exact": all(
                record.get("system") == system for record in records
            ),
            "record_source_exact": all(
                record.get("source_sample") == source for record in records
            ),
            "training_view_root_hashes": bool(records)
            and all(
                is_sha256(record.get("training_view_root_sha256", ""))
                for record in records
            ),
        }
        sources[source] = {
            "status": "PASS" if all(checks.values()) else "FAIL",
            "state": authority.get("state"),
            "full_source_manifest_sha256": full_manifest,
            "input_records_sha256": observed_records_sha,
            "expected_input_count": expected_input_count,
            "observed_input_count": len(records),
            "expected_occurrence_count": expected_occurrence_count,
            "observed_occurrence_count": len(occurrence_ids),
            "selected_view_training_rows": int(
                selected_rows_by_source.get(source, 0)
            ),
            "checks": checks,
        }
    status = (
        "PASS"
        if set(source_coverage_authority) == set(required)
        and all(record["status"] == "PASS" for record in sources.values())
        else "FAIL"
    )
    return {
        "schema": "THE134_EXACT_SOURCE_POPULATION_CLOSURE_V1",
        "status": status,
        "sources": sources,
        "semantic_sha256": canonical_json_sha256(sources),
    }


def main() -> int:
    args = parse_args()
    if args.scope == "full" and args.skip_input_hashes:
        raise SystemExit("--skip-input-hashes is forbidden for --scope full")
    if args.scope == "full" and args.source_provenance_json is None:
        raise SystemExit("--source-provenance-json is required for --scope full")
    external_provenance, source_coverage_authority = load_source_provenance(
        args.source_provenance_json
    )
    paths = expand_paths(args.input)
    if external_provenance:
        expected_path_set = {str(path) for path in paths}
        observed_path_set = set(external_provenance)
        if expected_path_set != observed_path_set:
            raise SystemExit(
                "source provenance/input path mismatch: "
                f"missing={sorted(expected_path_set - observed_path_set)} "
                f"extra={sorted(observed_path_set - expected_path_set)}"
            )
    path_sources: list[str] = []
    source_failures: list[str] = []
    for path in paths:
        try:
            path_sources.append(source_from_path(path, args.system, include_diagnostic=True))
        except ValueError as exc:
            source_failures.append(str(exc))
            path_sources.append("UNRESOLVED")
    observed_manifest_sources = set(path_sources) - {"UNRESOLVED"}
    required_sources = set(expected_sources(args.system))
    if args.scope == "full" and observed_manifest_sources != required_sources:
        source_failures.append(
            "full manifest source mismatch: "
            f"missing={sorted(required_sources - observed_manifest_sources)} "
            f"extra={sorted(observed_manifest_sources - required_sources)}"
        )

    pieces: dict[str, list[np.ndarray]] = defaultdict(list)
    file_reports: list[dict] = []
    input_hashes: dict[str, str | None] = {}
    for input_index, (path, source) in enumerate(zip(paths, path_sources)):
        report = {"path": str(path), "input_file_index": input_index, "source": source}
        try:
            if not path.is_file():
                raise ValueError("input ROOT does not exist")
            input_hashes[str(path)] = None if args.skip_input_hashes else sha256_file(path)
            if source == "UNRESOLVED":
                raise ValueError("source unresolved from path")
            external_record = external_provenance.get(str(path))
            if external_record is not None:
                if external_record.get("system") != args.system:
                    raise ValueError("external provenance system mismatch")
                if external_record.get("source_sample") != source:
                    raise ValueError("external provenance source mismatch")
                if args.scope == "full" and external_record.get(
                    "training_view_root_sha256"
                ) != input_hashes[str(path)]:
                    raise ValueError(
                        "external provenance/output ROOT SHA-256 mismatch"
                    )
            arrays, branches, metadata = read_tree(path, args.tree)
            metadata_failures = validate_artifact_metadata(
                metadata,
                arrays,
                external=external_provenance.get(str(path)) if external_provenance else None,
            )
            selected, vectors, validation = validate_file_arrays(
                arrays, system=args.system, source=source, view_name=args.view
            )
            validation["failures"] = metadata_failures + validation["failures"]
            validation["status"] = "PASS" if not validation["failures"] else "FAIL"
            report.update(validation)
            report["branches"] = branches
            report["root_metadata"] = metadata
            if validation["status"] != "PASS":
                file_reports.append(report)
                continue
            if not vectors:
                # A source shard may legitimately contain accepted candidates
                # that all carry AuAu cross-role label -1.  Preserve it as a
                # valid-empty file; aggregate source coverage is checked after
                # every shard has been inspected.
                file_reports.append(report)
                continue
            matrix = np.asarray(vectors, dtype=np.float32)
            for column_index, feature in enumerate(FEATURES_BY_SYSTEM[args.system]):
                pieces[feature].append(matrix[:, column_index])
            candidate_ids = _identity_pairs(arrays, "candidate_id")
            selected_indices = np.flatnonzero(selected)
            definitions = np.asarray([_text(item) for item in arrays["definition_name"]])
            for retained_view in ALL_SHOWER_VIEWS:
                index_by_candidate = {
                    candidate_id: index
                    for index, (candidate_id, definition) in enumerate(
                        zip(candidate_ids, definitions)
                    )
                    if definition == retained_view
                }
                retained_indices = np.asarray(
                    [index_by_candidate[candidate_ids[index]] for index in selected_indices],
                    dtype=np.int64,
                )
                retained_matrix = np.asarray(
                    [arrays["ordered_features"][index] for index in retained_indices],
                    dtype=np.float32,
                )
                for column_index, feature in enumerate(FEATURES_BY_SYSTEM[args.system]):
                    pieces[f"view_{retained_view}__{feature}"].append(
                        retained_matrix[:, column_index]
                    )
            for name in (
                "run",
                "event_sequence",
                "encounter_ordinal",
                "cluster_index",
                "cluster_map_key",
                "training_label",
                "is_signal",
                "source_role",
                "source_sample_code",
                "ppg12_source_role_label",
                "candidate_id_hi",
                "candidate_id_lo",
                "event_id_hi",
                "event_id_lo",
                "finite_feature_state",
                "model_domain_state",
                "nominal_training_eligible",
                "minimum_bias_classifier_decision",
                "weight_application_count",
                "event_weight",
                "weight_slice",
                "weight_cross_section",
                "weight_vertex",
                "weight_si_di",
                "weight_period",
                "weight_exposure",
                "weight_final",
            ):
                pieces[name].append(np.asarray(arrays[name])[selected])
            pieces["label_authority"].append(
                np.asarray([_text(item) for item in arrays["label_authority"]])[selected]
            )
            for name in (
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
            ):
                pieces[name].append(np.asarray(arrays[name])[selected])
            n_selected = int(np.sum(selected))
            pieces["source_sample"].append(np.full(n_selected, source, dtype="U40"))
            pieces["input_file_index"].append(np.full(n_selected, input_index, dtype=np.int64))
            pieces["input_tree_entry"].append(np.flatnonzero(selected).astype(np.int64))
            pieces["definition_name"].append(np.full(n_selected, args.view, dtype="U4"))
            pieces["shower_semantic_sha256"].append(
                np.full(n_selected, shower_semantic_sha256(args.view), dtype="U64")
            )
        except Exception as exc:  # preserve the exact first-bad file
            report.update({"status": "FAIL", "failures": [f"{type(exc).__name__}: {exc}"]})
        file_reports.append(report)

    failures = list(source_failures)
    failures.extend(
        f"{report['path']}: {failure}"
        for report in file_reports
        if report.get("status") != "PASS"
        for failure in report.get("failures", ["unknown validation failure"])
    )
    selected_rows_by_source = {
        source: int(
            sum(
                int(report.get("selected_view_training_rows", 0))
                for report in file_reports
                if report.get("status") == "PASS" and report.get("source") == source
            )
        )
        for source in sorted(required_sources)
    }
    if args.scope == "full":
        source_population_closure = build_source_population_closure(
            system=args.system,
            paths=paths,
            path_sources=path_sources,
            external_provenance=external_provenance,
            source_coverage_authority=source_coverage_authority,
            file_reports=file_reports,
            selected_rows_by_source=selected_rows_by_source,
        )
        if source_population_closure["status"] != "PASS":
            failures.append(
                "full source population closure failed; exact per-source manifest, "
                "input-count, occurrence-count, and input-record-set authority is required"
            )
    else:
        source_population_closure = {
            "schema": "THE134_EXACT_SOURCE_POPULATION_CLOSURE_V1",
            "status": "NOT_APPLICABLE_SMOKE",
            "sources": {},
            "semantic_sha256": canonical_json_sha256({}),
        }
    if not pieces:
        failures.append(f"no valid {args.view} rows were materialized")

    payload: dict[str, np.ndarray] = {}
    if pieces:
        payload = {name: np.concatenate(values) for name, values in pieces.items()}
        payload["evt"] = np.asarray(payload["event_sequence"], dtype=np.int64)
        payload["global_event_key"] = np.asarray(
            [
                identity128_hex((int(hi), int(lo)))
                for hi, lo in zip(
                    payload["event_id_hi"],
                    payload["event_id_lo"],
                )
            ],
            dtype="U32",
        )
        et = np.asarray(payload["cluster_Et"], dtype=np.float64)
        eta = np.asarray(payload["cluster_Eta"], dtype=np.float64)
        domain = (
            np.isfinite(et)
            & np.isfinite(eta)
            & (et >= MODEL_DOMAIN_GEV[0])
            & (et < MODEL_DOMAIN_GEV[1])
            & (np.abs(eta) < ETA_ABS_MAX)
        )
        for name in list(payload):
            payload[name] = payload[name][domain]
        # The generic trainer consumes ``is_signal``.  Preserve the legacy
        # nominal label separately and expose the already-audited frozen
        # training label without recomputing or relabeling it downstream.
        payload["nominal_is_signal"] = np.asarray(payload["is_signal"], dtype=np.int8)
        payload["is_signal"] = np.asarray(payload["training_label"], dtype=np.int8)
        candidate_pairs = list(
            zip(payload["candidate_id_hi"].tolist(), payload["candidate_id_lo"].tolist())
        )
        duplicate_candidates = len(candidate_pairs) - len(set(candidate_pairs))
        if duplicate_candidates:
            failures.append(
                f"duplicate {args.view} candidate identities after domain filter="
                f"{duplicate_candidates}"
            )
        observed_labels = set(np.asarray(payload["training_label"], dtype=np.int64).tolist())
        if observed_labels != {0, 1}:
            failures.append(f"matrix does not contain exactly both classes: {sorted(observed_labels)}")
        payload["__columns__"] = np.asarray(sorted(payload), dtype=object)

    status = "PASS" if not failures else "FAIL"
    audit = {
        "schema": "THE134_FACTORIAL_VIEW_TRAINING_MATRIX_AUDIT_V1",
        "status": status,
        "system": args.system,
        "scope": args.scope,
        "full_training_authority": int(
            status == "PASS"
            and args.scope == "full"
            and source_population_closure.get("status") == "PASS"
        ),
        "tree": args.tree,
        "shower_definition": args.view,
        "shower_semantic_sha256": shower_semantic_sha256(args.view),
        "feature_order": list(FEATURES_BY_SYSTEM[args.system]),
        "model_domain_gev": list(MODEL_DOMAIN_GEV),
        "eta_abs_max": ETA_ABS_MAX,
        "required_sources": sorted(required_sources),
        "observed_manifest_sources": sorted(observed_manifest_sources),
        "selected_view_training_rows_by_required_source": selected_rows_by_source,
        "source_coverage_authority": source_coverage_authority,
        "source_population_closure": source_population_closure,
        "input_count": len(paths),
        "input_hashes": input_hashes,
        "source_provenance_json": (
            str(args.source_provenance_json) if args.source_provenance_json else None
        ),
        "source_provenance_json_sha256": (
            sha256_file(args.source_provenance_json)
            if args.source_provenance_json is not None
            else None
        ),
        "file_reports": file_reports,
        "selected_rows": int(len(payload.get("cluster_Et", []))),
        "selected_class_counts": {
            str(label): int(np.sum(payload.get("training_label", np.asarray([])) == label))
            for label in (0, 1)
        },
        "failures": failures,
    }
    args.audit_out.parent.mkdir(parents=True, exist_ok=True)
    if status == "PASS":
        args.matrix_out.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(args.matrix_out, **payload)
        audit["matrix"] = str(args.matrix_out)
        audit["matrix_sha256"] = sha256_file(args.matrix_out)
    args.audit_out.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n")
    print(args.audit_out)
    if status != "PASS":
        print(json.dumps({"status": status, "failures": failures[:20]}, indent=2), file=sys.stderr)
        return 2
    print(args.matrix_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
