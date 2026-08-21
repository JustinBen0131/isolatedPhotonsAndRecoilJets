#!/usr/bin/env python3
"""Validate one populated p+p THE-134 photon-training sidecar.

This validator implements the established artifact-specific
``RJ_ARTIFACT_HEALTH_PROFILE_V1 / photon_training_multiview_v1`` contract.  The
sidecar is a highly compressed normalized table, so its authority is structural
rather than a generic 50 kB file-size heuristic.  A small-size signal is
reported as a truncation diagnostic, while ROOT health and semantic structure
remain authoritative.  The existing 50 kB rule for analysis/writer ROOT outputs
is deliberately outside this validator and is not changed here.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import uproot


REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.ml.contracts.the134_h70_contract import (  # noqa: E402
    ALL_SHOWER_VIEWS,
    CORE_BRANCHES,
    TRAINING_SCHEMA_NAME,
    TRAINING_SCHEMA_VERSION,
    feature_contract_sha256,
    identity128_from_text,
    identity128_hex,
    is_sha256,
    shower_semantic_sha256,
    shower_semantic_text,
    training_schema_sha256,
)


PROFILE_SCHEMA = "RJ_ARTIFACT_HEALTH_PROFILE_V1"
PROFILE_NAME = "photon_training_multiview_v1"
CERTIFICATE_SCHEMA = "THE134_ARTIFACT_HEALTH_CERTIFICATE_V1"
TREE_NAME = "RJPhotonTrainingViewV1"
GROSS_TRUNCATION_FLOOR_BYTES = 4096
PP_SYSTEM_CODE = 1
PP_FEATURE_COUNT = 11
VALIDATED_DOMAIN = 0
DIAGNOSTIC_EXTRAPOLATION = 1

METADATA_KEYS = (
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
EXPECTED_ROOT_KEYS = frozenset((*METADATA_KEYS, TREE_NAME))
EXPECTED_BRANCHES = frozenset(CORE_BRANCHES)
TEXT_BRANCHES = frozenset(
    {
        "definition_name",
        "shower_semantic_sha256",
        "feature_contract_sha256",
        "source_lane",
        "source_dataset",
        "source_sample",
        "source_period",
        "source_si_di_role",
        "source_ownership_state",
        "source_manifest_sha256",
        "input_uri_sha256",
        "input_file_sha256",
        "label_authority",
    }
)
IDENTITY_PREFIXES = (
    "training_view_id",
    "source_occurrence_id",
    "event_id",
    "candidate_id",
    "definition_id",
)
WEIGHT_BRANCHES = (
    "weight_slice",
    "weight_cross_section",
    "weight_vertex",
    "weight_si_di",
    "weight_period",
    "weight_exposure",
    "weight_final",
    "event_weight",
)


@dataclass(frozen=True)
class ExpectedContract:
    config_sha256: str
    code_sha256: str
    source_manifest_sha256: str
    source_lane: str
    source_dataset: str
    source_sample: str
    source_period: str
    source_si_di_role: str
    source_ownership_state: str
    source_run: int
    source_segment: int
    input_uri_sha256: str
    input_file_sha256: str


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _text(value: Any) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8")
    if isinstance(value, np.generic):
        value = value.item()
    return str(value)


def _identity_pairs(
    arrays: Mapping[str, Sequence[Any]], prefix: str
) -> list[tuple[int, int]]:
    return [
        (int(hi), int(lo))
        for hi, lo in zip(arrays[f"{prefix}_hi"], arrays[f"{prefix}_lo"])
    ]


def _relative_close(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    scale = np.maximum.reduce(
        [np.ones(len(left), dtype=np.float64), np.abs(left), np.abs(right)]
    )
    return np.abs(left - right) <= 1.0e-12 * scale


def _constant_text(
    arrays: Mapping[str, Sequence[Any]],
    branch: str,
    expected: str,
    failures: list[str],
) -> None:
    values = np.asarray([_text(value) for value in arrays[branch]], dtype=object)
    if np.any(values != expected):
        failures.append(
            f"{branch}_mismatch_rows={int(np.sum(values != expected))}"
        )


def _constant_int(
    arrays: Mapping[str, Sequence[Any]],
    branch: str,
    expected: int,
    failures: list[str],
) -> None:
    values = np.asarray(arrays[branch], dtype=np.int64)
    if np.any(values != expected):
        failures.append(
            f"{branch}_mismatch_rows={int(np.sum(values != expected))}"
        )


def validate_payload(
    *,
    metadata: Mapping[str, str],
    branch_names: Sequence[str],
    arrays: Mapping[str, Sequence[Any]],
    expected: ExpectedContract,
) -> dict[str, Any]:
    """Validate an already decoded sidecar payload.

    Keeping this logic independent of ROOT I/O makes every invariant
    mutation-testable while file health remains a separate mandatory gate.
    """

    failures: list[str] = []
    metadata_keys = set(metadata)
    if metadata_keys != set(METADATA_KEYS):
        failures.append(
            "metadata_inventory_mismatch:"
            f"missing={sorted(set(METADATA_KEYS) - metadata_keys)}:"
            f"extra={sorted(metadata_keys - set(METADATA_KEYS))}"
        )
    observed_branches = set(branch_names)
    if observed_branches != EXPECTED_BRANCHES:
        failures.append(
            "branch_inventory_mismatch:"
            f"missing={sorted(EXPECTED_BRANCHES - observed_branches)}:"
            f"extra={sorted(observed_branches - EXPECTED_BRANCHES)}"
        )

    expected_metadata = {
        "rj_photon_training_schema": TRAINING_SCHEMA_NAME,
        "rj_photon_training_schema_version": TRAINING_SCHEMA_VERSION,
        "schema_sha256": training_schema_sha256(),
        "pp_feature_contract_sha256": feature_contract_sha256("pp"),
        "auau_feature_contract_sha256": feature_contract_sha256("auau"),
        "source_manifest_sha256": expected.source_manifest_sha256,
        "config_sha256": expected.config_sha256,
        "code_sha256": expected.code_sha256,
        "rj_photon_training_complete": "1",
    }
    for name, value in expected_metadata.items():
        if metadata.get(name) != value:
            failures.append(
                f"metadata_mismatch:{name}:{metadata.get(name)!r}!={value!r}"
            )

    missing_arrays = sorted(EXPECTED_BRANCHES - set(arrays))
    if missing_arrays:
        failures.append(f"decoded_arrays_missing={missing_arrays}")
        return {
            "status": "FAIL",
            "profile_schema": PROFILE_SCHEMA,
            "profile": PROFILE_NAME,
            "failures": failures,
        }

    n_rows = len(arrays["training_view_id_hi"])
    if n_rows <= 0:
        failures.append("populated_profile_has_zero_entries")
    for branch in EXPECTED_BRANCHES:
        if len(arrays[branch]) != n_rows:
            failures.append(
                f"branch_length_mismatch:{branch}:{len(arrays[branch])}!={n_rows}"
            )
    try:
        recorded_entries = int(metadata.get("rj_photon_training_entries", ""))
    except ValueError:
        recorded_entries = -1
        failures.append("invalid_completion_entry_count")
    if recorded_entries != n_rows:
        failures.append(f"completion_entry_count_mismatch:{recorded_entries}!={n_rows}")

    for value in (
        expected.config_sha256,
        expected.code_sha256,
        expected.source_manifest_sha256,
        expected.input_uri_sha256,
        expected.input_file_sha256,
    ):
        if not is_sha256(value):
            failures.append(f"malformed_expected_sha256:{value!r}")

    identities = {
        prefix: _identity_pairs(arrays, prefix) for prefix in IDENTITY_PREFIXES
    }
    for prefix, values in identities.items():
        null_count = sum(value == (0, 0) for value in values)
        if null_count:
            failures.append(f"null_{prefix}_rows={null_count}")
    duplicate_views = n_rows - len(set(identities["training_view_id"]))
    if duplicate_views:
        failures.append(f"duplicate_training_view_ids={duplicate_views}")

    source_ids = identities["source_occurrence_id"]
    if len(set(source_ids)) != 1:
        failures.append(f"source_occurrence_count={len(set(source_ids))}")
    source_fields = {
        "source_lane": expected.source_lane,
        "source_dataset": expected.source_dataset,
        "source_sample": expected.source_sample,
        "source_period": expected.source_period,
        "source_si_di_role": expected.source_si_di_role,
        "source_ownership_state": expected.source_ownership_state,
        "source_manifest_sha256": expected.source_manifest_sha256,
        "input_uri_sha256": expected.input_uri_sha256,
        "input_file_sha256": expected.input_file_sha256,
    }
    for branch, value in source_fields.items():
        _constant_text(arrays, branch, value, failures)
    _constant_int(arrays, "source_run", expected.source_run, failures)
    _constant_int(arrays, "source_segment", expected.source_segment, failures)
    _constant_int(arrays, "run", expected.source_run, failures)
    _constant_int(arrays, "system_code", PP_SYSTEM_CODE, failures)

    expected_source_id = identity128_from_text(
        "|".join(
            (
                expected.source_lane,
                expected.source_dataset,
                expected.source_sample,
                expected.source_period,
                str(expected.source_run),
                str(expected.source_segment),
                expected.input_uri_sha256,
                expected.input_file_sha256,
                expected.source_manifest_sha256,
            )
        )
    )
    source_identity_mismatch = sum(value != expected_source_id for value in source_ids)
    if source_identity_mismatch:
        failures.append(
            f"source_occurrence_identity_mismatch_rows={source_identity_mismatch}"
        )

    event_ids = identities["event_id"]
    expected_event_ids = [
        identity128_from_text(
            "|".join(
                (
                    _text(lane),
                    _text(sample),
                    str(int(run)),
                    str(int(segment)),
                    str(int(sequence)),
                )
            )
        )
        for lane, sample, run, segment, sequence in zip(
            arrays["source_lane"],
            arrays["source_sample"],
            arrays["run"],
            arrays["source_segment"],
            arrays["event_sequence"],
        )
    ]
    event_identity_mismatch = sum(
        observed != wanted
        for observed, wanted in zip(event_ids, expected_event_ids)
    )
    if event_identity_mismatch:
        failures.append(f"event_identity_mismatch_rows={event_identity_mismatch}")

    candidate_ids = identities["candidate_id"]
    expected_candidate_ids = [
        identity128_from_text(
            f"{identity128_hex(event_id)}|candidate|{int(ordinal)}|{int(map_key)}"
        )
        for event_id, ordinal, map_key in zip(
            event_ids, arrays["encounter_ordinal"], arrays["cluster_map_key"]
        )
    ]
    candidate_identity_mismatch = sum(
        observed != wanted
        for observed, wanted in zip(candidate_ids, expected_candidate_ids)
    )
    if candidate_identity_mismatch:
        failures.append(
            f"candidate_identity_mismatch_rows={candidate_identity_mismatch}"
        )
    event_by_candidate: dict[tuple[int, int], set[tuple[int, int]]] = defaultdict(set)
    for candidate_id, event_id in zip(candidate_ids, event_ids):
        event_by_candidate[candidate_id].add(event_id)
    if any(len(events) != 1 for events in event_by_candidate.values()):
        failures.append("candidate_to_event_foreign_key_defect")

    definitions = [_text(value) for value in arrays["definition_name"]]
    semantics = [_text(value) for value in arrays["shower_semantic_sha256"]]
    definition_ids = identities["definition_id"]
    expected_definition_ids: list[tuple[int, int]] = []
    bad_definition_names = 0
    bad_semantics = 0
    for name, semantic in zip(definitions, semantics):
        if name not in ALL_SHOWER_VIEWS:
            bad_definition_names += 1
            expected_definition_ids.append((0, 0))
            continue
        expected_semantic = shower_semantic_sha256(name)
        if semantic != expected_semantic:
            bad_semantics += 1
        expected_definition_ids.append(
            identity128_from_text(
                f"shower-definition|{name}|{shower_semantic_text(name)}"
            )
        )
    if bad_definition_names:
        failures.append(f"unknown_definition_rows={bad_definition_names}")
    if bad_semantics:
        failures.append(f"shower_semantic_mismatch_rows={bad_semantics}")
    definition_identity_mismatch = sum(
        observed != wanted
        for observed, wanted in zip(definition_ids, expected_definition_ids)
    )
    if definition_identity_mismatch:
        failures.append(
            f"definition_identity_mismatch_rows={definition_identity_mismatch}"
        )

    expected_training_ids = [
        identity128_from_text(
            f"{identity128_hex(candidate_id)}|training-view|"
            f"{identity128_hex(definition_id)}"
        )
        for candidate_id, definition_id in zip(candidate_ids, definition_ids)
    ]
    training_identity_mismatch = sum(
        observed != wanted
        for observed, wanted in zip(
            identities["training_view_id"], expected_training_ids
        )
    )
    if training_identity_mismatch:
        failures.append(
            f"training_view_identity_mismatch_rows={training_identity_mismatch}"
        )

    definitions_by_candidate: dict[tuple[int, int], list[str]] = defaultdict(list)
    for candidate_id, definition in zip(candidate_ids, definitions):
        definitions_by_candidate[candidate_id].append(definition)
    expected_views = set(ALL_SHOWER_VIEWS)
    factorial_failures = sum(
        len(observed) != len(ALL_SHOWER_VIEWS)
        or len(set(observed)) != len(observed)
        or set(observed) != expected_views
        for observed in definitions_by_candidate.values()
    )
    if factorial_failures:
        failures.append(f"seven_view_candidate_failures={factorial_failures}")

    feature_contracts = np.asarray(
        [_text(value) for value in arrays["feature_contract_sha256"]], dtype=object
    )
    expected_feature_contract = feature_contract_sha256("pp")
    if np.any(feature_contracts != expected_feature_contract):
        failures.append(
            "pp_feature_contract_mismatch_rows="
            f"{int(np.sum(feature_contracts != expected_feature_contract))}"
        )
    feature_counts = np.asarray(arrays["feature_count"], dtype=np.int64)
    if np.any(feature_counts != PP_FEATURE_COUNT):
        failures.append(
            f"pp_feature_count_mismatch_rows={int(np.sum(feature_counts != PP_FEATURE_COUNT))}"
        )
    feature_vectors = [np.asarray(row, dtype=np.float64) for row in arrays["ordered_features"]]
    bad_vector_size = sum(len(row) != PP_FEATURE_COUNT for row in feature_vectors)
    bad_vector_finite = sum(
        len(row) != PP_FEATURE_COUNT or not np.all(np.isfinite(row))
        for row in feature_vectors
    )
    if bad_vector_size:
        failures.append(f"pp_ordered_feature_length_mismatch_rows={bad_vector_size}")
    if bad_vector_finite:
        failures.append(f"pp_nonfinite_feature_rows={bad_vector_finite}")
    finite_states = np.asarray(arrays["finite_feature_state"], dtype=np.int64)
    if np.any(finite_states != 1):
        failures.append(
            f"finite_feature_state_mismatch_rows={int(np.sum(finite_states != 1))}"
        )

    ordinals = np.asarray(arrays["encounter_ordinal"], dtype=np.int64)
    cluster_indices = np.asarray(arrays["cluster_index"], dtype=np.int64)
    if np.any(ordinals != cluster_indices):
        failures.append(
            "cluster_index_ordinal_mismatch_rows="
            f"{int(np.sum(ordinals != cluster_indices))}"
        )
    sample_digits = ""
    for character in reversed(expected.source_sample):
        if not character.isdigit():
            break
        sample_digits = character + sample_digits
    if sample_digits:
        _constant_int(
            arrays, "source_sample_code", int(sample_digits), failures
        )
    _constant_text(arrays, "label_authority", "PPG12_SOURCE_ROLE", failures)
    labels = np.asarray(arrays["training_label"], dtype=np.int64)
    signal_labels = np.asarray(arrays["is_signal"], dtype=np.int64)
    if np.any(~np.isin(labels, (0, 1))):
        failures.append(
            f"nonbinary_training_label_rows={int(np.sum(~np.isin(labels, (0, 1))))}"
        )
    if np.any(~np.isin(signal_labels, (0, 1))):
        failures.append(
            f"nonbinary_signal_label_rows={int(np.sum(~np.isin(signal_labels, (0, 1))))}"
        )

    application_counts = np.asarray(
        arrays["weight_application_count"], dtype=np.int64
    )
    if np.any(application_counts != 1):
        failures.append(
            "weight_application_count_mismatch_rows="
            f"{int(np.sum(application_counts != 1))}"
        )
    weights = {
        branch: np.asarray(arrays[branch], dtype=np.float64)
        for branch in WEIGHT_BRANCHES
    }
    invalid_weight = np.logical_or.reduce(
        [~np.isfinite(values) | (values <= 0.0) for values in weights.values()]
    )
    if np.any(invalid_weight):
        failures.append(
            f"nonfinite_or_nonpositive_weight_rows={int(np.sum(invalid_weight))}"
        )
    slice_alias = _relative_close(
        weights["weight_slice"], weights["weight_cross_section"]
    )
    period_alias = _relative_close(
        weights["weight_period"], weights["weight_exposure"]
    )
    if np.any(~slice_alias):
        failures.append(
            f"pp_slice_cross_section_alias_mismatch_rows={int(np.sum(~slice_alias))}"
        )
    if np.any(~period_alias):
        failures.append(
            f"pp_period_exposure_alias_mismatch_rows={int(np.sum(~period_alias))}"
        )
    component_product = (
        weights["weight_slice"]
        * weights["weight_vertex"]
        * weights["weight_si_di"]
        * weights["weight_period"]
    )
    final_closure = _relative_close(component_product, weights["weight_final"])
    event_closure = _relative_close(
        weights["event_weight"], weights["weight_final"]
    )
    if np.any(~final_closure):
        failures.append(
            f"pp_weight_component_closure_mismatch_rows={int(np.sum(~final_closure))}"
        )
    if np.any(~event_closure):
        failures.append(
            f"event_final_weight_mismatch_rows={int(np.sum(~event_closure))}"
        )

    cluster_et = np.asarray(arrays["cluster_Et"], dtype=np.float64)
    domain = np.asarray(arrays["model_domain_state"], dtype=np.int64)
    retained_below15 = np.asarray(
        arrays["below15_retention_state"], dtype=np.int64
    )
    eligible = np.asarray(
        arrays["nominal_training_eligible"], dtype=np.int64
    )
    wp_state = np.asarray(arrays["working_point_state"], dtype=np.int64)
    tag_state = np.asarray(arrays["tag_state"], dtype=np.int64)
    below15 = cluster_et < 15.0
    in_domain = (cluster_et >= 15.0) & (cluster_et < 35.0)
    outside_domain = ~in_domain
    if np.any(
        below15
        & (
            (domain != DIAGNOSTIC_EXTRAPOLATION)
            | (retained_below15 != 1)
            | (eligible != 0)
            | (wp_state != -1)
            | (tag_state != -1)
        )
    ):
        failures.append("below15_diagnostic_null_wp_safety_violation")
    if np.any(
        in_domain
        & (
            (domain != VALIDATED_DOMAIN)
            | (retained_below15 != 0)
            | (eligible != 1)
            | (wp_state != -1)
            | (tag_state != -1)
        )
    ):
        failures.append("validated_domain_training_state_violation")
    if np.any(
        outside_domain
        & ~below15
        & (
            (domain != DIAGNOSTIC_EXTRAPOLATION)
            | (retained_below15 != 0)
            | (eligible != 0)
            | (wp_state != -1)
            | (tag_state != -1)
        )
    ):
        failures.append("above_domain_diagnostic_null_wp_safety_violation")

    return {
        "status": "PASS" if not failures else "FAIL",
        "profile_schema": PROFILE_SCHEMA,
        "profile": PROFILE_NAME,
        "entries": n_rows,
        "candidates": len(definitions_by_candidate),
        "views_per_candidate": len(ALL_SHOWER_VIEWS),
        "feature_count": PP_FEATURE_COUNT,
        "failures": failures,
    }


def _pyroot_health(path: Path) -> tuple[bool, str]:
    try:
        import ROOT  # type: ignore
    except Exception as exc:  # pragma: no cover - runtime-specific
        return False, f"pyroot_unavailable:{exc}"
    root_file = ROOT.TFile.Open(str(path), "READ")
    if not root_file:
        return False, "pyroot_open_returned_null"
    try:
        if root_file.IsZombie():
            return False, "root_is_zombie"
        if root_file.TestBit(ROOT.TFile.kRecovered):
            return False, "root_is_recovered"
        return True, "healthy"
    finally:
        root_file.Close()


def _read_payload(
    path: Path,
) -> tuple[dict[str, str], list[str], dict[str, np.ndarray], dict[str, str]]:
    import awkward as ak

    with uproot.open(path) as root_file:
        classnames: dict[str, str] = {}
        for raw_name, class_name in root_file.classnames().items():
            name = str(raw_name).split(";")[0]
            if name in classnames:
                raise ValueError(f"duplicate ROOT key cycles for {name}")
            classnames[name] = str(class_name)
        if set(classnames) != EXPECTED_ROOT_KEYS:
            raise ValueError(
                "ROOT key inventory mismatch:"
                f"missing={sorted(EXPECTED_ROOT_KEYS - set(classnames))}:"
                f"extra={sorted(set(classnames) - EXPECTED_ROOT_KEYS)}"
            )
        if not classnames[TREE_NAME].startswith("TTree"):
            raise ValueError(
                f"{TREE_NAME} class is {classnames[TREE_NAME]!r}, expected TTree"
            )
        for name in METADATA_KEYS:
            if classnames[name] != "TNamed":
                raise ValueError(
                    f"metadata key {name} class is {classnames[name]!r}, expected TNamed"
                )
        metadata = {
            name: _text(root_file[name].member("fTitle")) for name in METADATA_KEYS
        }
        tree = root_file[TREE_NAME]
        branch_names = sorted(str(name) for name in tree.keys())
        awkward_arrays = tree.arrays(list(CORE_BRANCHES), library="ak")
        arrays: dict[str, np.ndarray] = {}
        for name in CORE_BRANCHES:
            values = awkward_arrays[name]
            if name == "ordered_features":
                arrays[name] = np.asarray(ak.to_list(values), dtype=object)
            elif name in TEXT_BRANCHES:
                arrays[name] = np.asarray(ak.to_list(values), dtype=object)
            else:
                arrays[name] = ak.to_numpy(values)
        return metadata, branch_names, arrays, classnames


def validate_file(path: Path, expected: ExpectedContract) -> dict[str, Any]:
    result: dict[str, Any] = {
        "schema": CERTIFICATE_SCHEMA,
        "profile_schema": PROFILE_SCHEMA,
        "profile": PROFILE_NAME,
        "status": "FAIL",
        "path": str(path),
        "gross_truncation_floor_bytes": GROSS_TRUNCATION_FLOOR_BYTES,
        "minimum_bytes_mode": "diagnostic_only",
        "analysis_writer_50000_byte_gate_unchanged": True,
        "diagnostics": [],
        "failures": [],
    }
    if not path.is_file():
        result["failures"].append("artifact_missing")
        return result
    size = path.stat().st_size
    result["bytes"] = size
    result["sha256"] = sha256_file(path)
    if size < GROSS_TRUNCATION_FLOOR_BYTES:
        result["diagnostics"].append(
            f"gross_truncation:{size}<{GROSS_TRUNCATION_FLOOR_BYTES}"
        )

    healthy, health_detail = _pyroot_health(path)
    result["pyroot_health"] = health_detail
    if not healthy:
        result["failures"].append(health_detail)
        return result
    try:
        metadata, branch_names, arrays, classnames = _read_payload(path)
        payload = validate_payload(
            metadata=metadata,
            branch_names=branch_names,
            arrays=arrays,
            expected=expected,
        )
        result["root_key_classes"] = dict(sorted(classnames.items()))
        result["metadata"] = dict(sorted(metadata.items()))
        result["branch_count"] = len(branch_names)
        result["entries"] = payload.get("entries", 0)
        result["candidates"] = payload.get("candidates", 0)
        result["failures"].extend(payload["failures"])
    except Exception as exc:
        result["failures"].append(f"root_payload_error:{exc}")
    result["status"] = "PASS" if not result["failures"] else "FAIL"
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--expected-config-sha256", required=True)
    parser.add_argument("--expected-code-sha256", required=True)
    parser.add_argument("--expected-source-manifest-sha256", required=True)
    parser.add_argument("--expected-source-lane", required=True)
    parser.add_argument("--expected-source-dataset", required=True)
    parser.add_argument("--expected-source-sample", required=True)
    parser.add_argument("--expected-source-period", required=True)
    parser.add_argument("--expected-source-si-di-role", required=True)
    parser.add_argument("--expected-source-ownership-state", required=True)
    parser.add_argument("--expected-source-run", type=int, required=True)
    parser.add_argument("--expected-source-segment", type=int, required=True)
    parser.add_argument("--expected-input-uri-sha256", required=True)
    parser.add_argument("--expected-input-file-sha256", required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    expected = ExpectedContract(
        config_sha256=args.expected_config_sha256,
        code_sha256=args.expected_code_sha256,
        source_manifest_sha256=args.expected_source_manifest_sha256,
        source_lane=args.expected_source_lane,
        source_dataset=args.expected_source_dataset,
        source_sample=args.expected_source_sample,
        source_period=args.expected_source_period,
        source_si_di_role=args.expected_source_si_di_role,
        source_ownership_state=args.expected_source_ownership_state,
        source_run=args.expected_source_run,
        source_segment=args.expected_source_segment,
        input_uri_sha256=args.expected_input_uri_sha256,
        input_file_sha256=args.expected_input_file_sha256,
    )
    report = validate_file(args.input, expected)
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output_json.with_name(
        args.output_json.name + f".tmp.{os.getpid()}"
    )
    temporary.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, args.output_json)
    print(json.dumps(report, sort_keys=True))
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
