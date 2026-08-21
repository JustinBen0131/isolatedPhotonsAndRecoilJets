#!/usr/bin/env python3
"""Mutation tests for the THE-134 training-sidecar health profile."""

from __future__ import annotations

import copy
import sys
import unittest
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
VALIDATOR_DIR = HERE.parent
if str(VALIDATOR_DIR) not in sys.path:
    sys.path.insert(0, str(VALIDATOR_DIR))

import validate_the134_training_sidecar as validator


def _repeat(value: object, count: int) -> np.ndarray:
    return np.asarray([copy.deepcopy(value) for _ in range(count)], dtype=object)


def valid_fixture() -> tuple[
    dict[str, str],
    list[str],
    dict[str, np.ndarray],
    validator.ExpectedContract,
]:
    contract = validator.ExpectedContract(
        config_sha256="1" * 64,
        code_sha256="2" * 64,
        source_manifest_sha256="3" * 64,
        source_lane="pp_sim",
        source_dataset="run28",
        source_sample="run28_photonjet20",
        source_period="0mrad",
        source_si_di_role="SI",
        source_ownership_state="OWNED",
        source_run=54280,
        source_segment=7,
        input_uri_sha256="4" * 64,
        input_file_sha256="5" * 64,
    )
    views = list(validator.ALL_SHOWER_VIEWS)
    n_rows = len(views)
    source_id = validator.identity128_from_text(
        "|".join(
            (
                contract.source_lane,
                contract.source_dataset,
                contract.source_sample,
                contract.source_period,
                str(contract.source_run),
                str(contract.source_segment),
                contract.input_uri_sha256,
                contract.input_file_sha256,
                contract.source_manifest_sha256,
            )
        )
    )
    event_sequence = 9
    event_id = validator.identity128_from_text(
        "|".join(
            (
                contract.source_lane,
                contract.source_sample,
                str(contract.source_run),
                str(contract.source_segment),
                str(event_sequence),
            )
        )
    )
    encounter_ordinal = 2
    cluster_map_key = 91
    candidate_id = validator.identity128_from_text(
        f"{validator.identity128_hex(event_id)}|candidate|"
        f"{encounter_ordinal}|{cluster_map_key}"
    )
    definition_ids = [
        validator.identity128_from_text(
            f"shower-definition|{name}|{validator.shower_semantic_text(name)}"
        )
        for name in views
    ]
    training_ids = [
        validator.identity128_from_text(
            f"{validator.identity128_hex(candidate_id)}|training-view|"
            f"{validator.identity128_hex(definition_id)}"
        )
        for definition_id in definition_ids
    ]

    arrays: dict[str, np.ndarray] = {
        name: np.zeros(n_rows, dtype=np.int64)
        for name in validator.EXPECTED_BRANCHES
        if name not in validator.TEXT_BRANCHES and name != "ordered_features"
    }
    for name in validator.TEXT_BRANCHES:
        arrays[name] = _repeat("", n_rows)
    arrays["ordered_features"] = _repeat([float(index) for index in range(11)], n_rows)

    identity_values = {
        "training_view_id": training_ids,
        "source_occurrence_id": [source_id] * n_rows,
        "event_id": [event_id] * n_rows,
        "candidate_id": [candidate_id] * n_rows,
        "definition_id": definition_ids,
    }
    for prefix, values in identity_values.items():
        arrays[f"{prefix}_hi"] = np.asarray(
            [value[0] for value in values], dtype=np.uint64
        )
        arrays[f"{prefix}_lo"] = np.asarray(
            [value[1] for value in values], dtype=np.uint64
        )

    arrays["definition_name"] = _repeat("", n_rows)
    arrays["definition_name"][:] = views
    arrays["shower_semantic_sha256"] = np.asarray(
        [validator.shower_semantic_sha256(view) for view in views], dtype=object
    )
    arrays["feature_contract_sha256"] = _repeat(
        validator.feature_contract_sha256("pp"), n_rows
    )
    arrays["feature_count"][:] = 11
    arrays["finite_feature_state"][:] = 1
    arrays["system_code"][:] = 1
    for branch, value in {
        "source_lane": contract.source_lane,
        "source_dataset": contract.source_dataset,
        "source_sample": contract.source_sample,
        "source_period": contract.source_period,
        "source_si_di_role": contract.source_si_di_role,
        "source_ownership_state": contract.source_ownership_state,
        "source_manifest_sha256": contract.source_manifest_sha256,
        "input_uri_sha256": contract.input_uri_sha256,
        "input_file_sha256": contract.input_file_sha256,
    }.items():
        arrays[branch] = _repeat(value, n_rows)
    arrays["source_run"][:] = contract.source_run
    arrays["source_segment"][:] = contract.source_segment
    arrays["run"][:] = contract.source_run
    arrays["event_sequence"][:] = event_sequence
    arrays["encounter_ordinal"][:] = encounter_ordinal
    arrays["cluster_index"][:] = encounter_ordinal
    arrays["cluster_map_key"][:] = cluster_map_key
    arrays["cluster_Et"] = np.full(n_rows, 20.0)
    arrays["cluster_Eta"] = np.full(n_rows, 0.1)
    arrays["cluster_Phi"] = np.full(n_rows, 1.2)
    arrays["vertexz"] = np.full(n_rows, -2.0)
    arrays["centrality"] = np.full(n_rows, -1.0)
    arrays["model_domain_state"][:] = validator.VALIDATED_DOMAIN
    arrays["below15_retention_state"][:] = 0
    arrays["nominal_training_eligible"][:] = 1
    arrays["working_point_state"][:] = -1
    arrays["tag_state"][:] = -1
    arrays["training_label"][:] = 1
    arrays["is_signal"][:] = 1
    arrays["label_authority"] = _repeat("PPG12_SOURCE_ROLE", n_rows)
    arrays["source_role"][:] = 1
    arrays["source_sample_code"][:] = 20
    arrays["ppg12_source_role_label"][:] = 1
    arrays["npb_label"][:] = -1
    arrays["is_npb"][:] = -1
    arrays["minimum_bias_classifier_decision"][:] = -1
    for branch in validator.WEIGHT_BRANCHES:
        arrays[branch] = np.ones(n_rows, dtype=np.float64)
    arrays["weight_application_count"][:] = 1

    metadata = {
        "rj_photon_training_schema": validator.TRAINING_SCHEMA_NAME,
        "rj_photon_training_schema_version": validator.TRAINING_SCHEMA_VERSION,
        "schema_sha256": validator.training_schema_sha256(),
        "pp_feature_contract_sha256": validator.feature_contract_sha256("pp"),
        "auau_feature_contract_sha256": validator.feature_contract_sha256("auau"),
        "source_manifest_sha256": contract.source_manifest_sha256,
        "config_sha256": contract.config_sha256,
        "code_sha256": contract.code_sha256,
        "rj_photon_training_complete": "1",
        "rj_photon_training_entries": str(n_rows),
    }
    return metadata, sorted(validator.EXPECTED_BRANCHES), arrays, contract


class TrainingSidecarHealthTests(unittest.TestCase):
    def validate(
        self,
        metadata: dict[str, str],
        branches: list[str],
        arrays: dict[str, np.ndarray],
        contract: validator.ExpectedContract,
    ) -> dict[str, object]:
        return validator.validate_payload(
            metadata=metadata,
            branch_names=branches,
            arrays=arrays,
            expected=contract,
        )

    def test_valid_populated_pp_sidecar_passes(self) -> None:
        metadata, branches, arrays, contract = valid_fixture()
        report = self.validate(metadata, branches, arrays, contract)
        self.assertEqual(report["status"], "PASS", report["failures"])
        self.assertEqual(report["profile_schema"], validator.PROFILE_SCHEMA)
        self.assertEqual(report["profile"], validator.PROFILE_NAME)
        self.assertEqual(report["entries"], 7)
        self.assertEqual(report["candidates"], 1)

    def test_missing_completion_is_rejected(self) -> None:
        metadata, branches, arrays, contract = valid_fixture()
        del metadata["rj_photon_training_complete"]
        report = self.validate(metadata, branches, arrays, contract)
        self.assertEqual(report["status"], "FAIL")
        self.assertTrue(
            any("rj_photon_training_complete" in item for item in report["failures"])
        )

    def test_wrong_frozen_hash_is_rejected(self) -> None:
        metadata, branches, arrays, contract = valid_fixture()
        metadata["code_sha256"] = "f" * 64
        report = self.validate(metadata, branches, arrays, contract)
        self.assertEqual(report["status"], "FAIL")
        self.assertTrue(
            any("metadata_mismatch:code_sha256" in item for item in report["failures"])
        )

    def test_duplicate_view_is_rejected(self) -> None:
        metadata, branches, arrays, contract = valid_fixture()
        arrays["definition_name"][1] = arrays["definition_name"][0]
        report = self.validate(metadata, branches, arrays, contract)
        self.assertEqual(report["status"], "FAIL")
        self.assertTrue(
            any("seven_view_candidate_failures" in item for item in report["failures"])
        )

    def test_orphan_identity_is_rejected(self) -> None:
        metadata, branches, arrays, contract = valid_fixture()
        arrays["event_id_lo"][3] = np.uint64(int(arrays["event_id_lo"][3]) ^ 1)
        report = self.validate(metadata, branches, arrays, contract)
        self.assertEqual(report["status"], "FAIL")
        self.assertTrue(
            any(
                "event_identity_mismatch" in item
                or "candidate_to_event_foreign_key_defect" in item
                for item in report["failures"]
            )
        )

    def test_below15_nonnull_tag_is_rejected(self) -> None:
        metadata, branches, arrays, contract = valid_fixture()
        arrays["cluster_Et"][:] = 12.0
        arrays["model_domain_state"][:] = validator.DIAGNOSTIC_EXTRAPOLATION
        arrays["below15_retention_state"][:] = 1
        arrays["nominal_training_eligible"][:] = 0
        arrays["tag_state"][:] = 1
        report = self.validate(metadata, branches, arrays, contract)
        self.assertEqual(report["status"], "FAIL")
        self.assertIn(
            "below15_diagnostic_null_wp_safety_violation", report["failures"]
        )

    def test_weight_double_application_is_rejected(self) -> None:
        metadata, branches, arrays, contract = valid_fixture()
        arrays["weight_application_count"][0] = 2
        report = self.validate(metadata, branches, arrays, contract)
        self.assertEqual(report["status"], "FAIL")
        self.assertTrue(
            any("weight_application_count" in item for item in report["failures"])
        )


if __name__ == "__main__":
    unittest.main()
