#!/usr/bin/env python3
"""Mutation tests for DirectHistogramReplayObligationV1."""

from __future__ import annotations

import copy
import json
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import uproot


HERE = Path(__file__).resolve().parent
VALIDATOR_DIR = HERE.parent
if str(VALIDATOR_DIR) not in sys.path:
    sys.path.insert(0, str(VALIDATOR_DIR))

import validate_direct_histogram_replay_obligations as validator


def inventory_fixture() -> dict[str, object]:
    entries = [
        {
            "root_key": "qa/photon_et",
            "object_class": "TH1D",
            "directory": "qa",
            "name": "photon_et",
            "title": "Photon E_{T}",
            "axes": [
                {
                    "name": "xaxis",
                    "title": "E_{T} [GeV]",
                    "edges": [0.0, 15.0, 35.0, 50.0],
                }
            ],
        }
    ]
    return {
        "schema": validator.INVENTORY_SCHEMA,
        "schema_version": 1,
        "artifact_id": "pp_data_direct",
        "system": "pp",
        "lane": "data",
        "root_sha256": "1" * 64,
        "histogram_count": len(entries),
        "entries": entries,
        "inventory_sha256": validator.sha256_json(entries),
    }


def obligation_fixture() -> dict[str, object]:
    inventory = inventory_fixture()
    surface = copy.deepcopy(inventory["entries"][0])
    obligation = {
        **surface,
        "source_inputs": [
            {
                "tree": "RJPhotonCandidateV1",
                "branches": ["candidate_id_hi", "candidate_id_lo", "et"],
            },
            {
                "tree": "RJWeightComponentV1",
                "branches": ["candidate_id_hi", "candidate_id_lo", "weight_final"],
            },
        ],
        "selection_rule": "15 <= et && et < 35",
        "fill_rule": "photon_et.Fill(et, weight_final)",
        "ordering_policy": "stable candidate_id ascending",
        "tie_policy": "candidate_id breaks equal-et ties",
        "weight_contract": {
            "components": ["weight_final"],
            "application_count": 1,
            "rule": "use the stored final weight exactly once",
        },
        "flow_policy": {
            "underflow": "ROOT native underflow",
            "overflow": "ROOT native overflow",
        },
        "response_classification": "NOT_APPLICABLE",
        "equality_contract": {
            "status": "REPLAY_EXACT",
            "numerical_mode": "EXACT",
            "compared_fields": sorted(validator.REQUIRED_COMPARE_FIELDS),
        },
        "analysis_consumers": ["THE-124/photon_qa"],
        "ian_consumers": ["IAN/photon_et"],
        "recipe_identity_sha256": "",
    }
    obligation["recipe_identity_sha256"] = validator.obligation_recipe_identity(
        obligation
    )
    return {
        "schema": validator.REGISTRY_SCHEMA,
        "schema_version": 1,
        "direct_reference": {
            "artifact_id": inventory["artifact_id"],
            "system": inventory["system"],
            "lane": inventory["lane"],
            "root_sha256": inventory["root_sha256"],
            "inventory_sha256": inventory["inventory_sha256"],
            "bundle_sha256": "2" * 64,
            "code_sha256": "3" * 64,
            "schema_sha256": "4" * 64,
            "config_sha256": "5" * 64,
            "model_registry_sha256": "6" * 64,
            "source_manifest_sha256": "7" * 64,
        },
        "obligations": [obligation],
    }


class DirectHistogramReplayObligationTests(unittest.TestCase):
    def test_complete_exact_registry_passes(self) -> None:
        report = validator.validate_registry(obligation_fixture(), inventory_fixture())
        self.assertEqual(report["status"], "PASS", report["failures"])
        self.assertEqual(report["coverage_fraction"], 1.0)
        self.assertEqual(report["unregistered_count"], 0)

    def test_missing_histogram_recipe_is_rejected(self) -> None:
        registry = obligation_fixture()
        registry["obligations"] = []
        report = validator.validate_registry(registry, inventory_fixture())
        self.assertEqual(report["status"], "FAIL")
        self.assertTrue(
            any(
                failure.startswith("unregistered_direct_histograms=")
                for failure in report["failures"]
            )
        )

    def test_dst_required_escape_hatch_is_rejected(self) -> None:
        registry = obligation_fixture()
        registry["obligations"][0]["response_classification"] = "DST_REQUIRED"
        registry["obligations"][0][
            "recipe_identity_sha256"
        ] = validator.obligation_recipe_identity(registry["obligations"][0])
        report = validator.validate_registry(registry, inventory_fixture())
        self.assertEqual(report["status"], "FAIL")
        self.assertTrue(
            any(
                failure.startswith("forbidden_state:")
                for failure in report["failures"]
            )
        )

    def test_weight_double_application_is_rejected(self) -> None:
        registry = obligation_fixture()
        registry["obligations"][0]["weight_contract"]["application_count"] = 2
        registry["obligations"][0][
            "recipe_identity_sha256"
        ] = validator.obligation_recipe_identity(registry["obligations"][0])
        report = validator.validate_registry(registry, inventory_fixture())
        self.assertEqual(report["status"], "FAIL")
        self.assertTrue(
            any(
                "weight_application_count_not_one" in failure
                for failure in report["failures"]
            )
        )

    def test_partial_equality_contract_is_rejected(self) -> None:
        registry = obligation_fixture()
        registry["obligations"][0]["equality_contract"]["compared_fields"].remove(
            "sumw2"
        )
        registry["obligations"][0][
            "recipe_identity_sha256"
        ] = validator.obligation_recipe_identity(registry["obligations"][0])
        report = validator.validate_registry(registry, inventory_fixture())
        self.assertEqual(report["status"], "FAIL")
        self.assertTrue(
            any(
                "incomplete_compared_fields" in failure
                for failure in report["failures"]
            )
        )

    def test_recipe_identity_drift_is_rejected(self) -> None:
        registry = obligation_fixture()
        registry["obligations"][0]["selection_rule"] = "et > 15"
        report = validator.validate_registry(registry, inventory_fixture())
        self.assertEqual(report["status"], "FAIL")
        self.assertTrue(
            any(
                "recipe_identity_sha256_mismatch" in failure
                for failure in report["failures"]
            )
        )

    def test_root_semantic_comparator_accepts_identical_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            direct = Path(directory) / "direct.root"
            replay = Path(directory) / "replay.root"
            values = np.asarray([2.0, 3.0, 5.0])
            edges = np.asarray([0.0, 15.0, 35.0, 50.0])
            for path in (direct, replay):
                with uproot.recreate(path) as root_file:
                    root_file["qa/photon_et"] = (values, edges, "Photon E_{T}")
            registry = obligation_fixture()
            inventory = validator.build_inventory(
                direct,
                artifact_id="pp_data_direct",
                system="pp",
                lane="data",
            )
            registry["direct_reference"]["root_sha256"] = inventory["root_sha256"]
            registry["direct_reference"]["inventory_sha256"] = inventory[
                "inventory_sha256"
            ]
            registry["obligations"][0].update(inventory["entries"][0])
            registry["obligations"][0][
                "recipe_identity_sha256"
            ] = validator.obligation_recipe_identity(registry["obligations"][0])
            validation = validator.validate_registry(registry, inventory)
            self.assertEqual(validation["status"], "PASS", validation["failures"])
            comparison = validator.compare_outputs(registry, direct, replay)
            self.assertEqual(comparison["status"], "PASS", comparison["failures"])

    def test_root_semantic_comparator_rejects_bin_mutation(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            direct = Path(directory) / "direct.root"
            replay = Path(directory) / "replay.root"
            edges = np.asarray([0.0, 15.0, 35.0, 50.0])
            with uproot.recreate(direct) as root_file:
                root_file["qa/photon_et"] = (
                    np.asarray([2.0, 3.0, 5.0]),
                    edges,
                    "Photon E_{T}",
                )
            with uproot.recreate(replay) as root_file:
                root_file["qa/photon_et"] = (
                    np.asarray([2.0, 4.0, 5.0]),
                    edges,
                    "Photon E_{T}",
                )
            registry = obligation_fixture()
            inventory = validator.build_inventory(
                direct,
                artifact_id="pp_data_direct",
                system="pp",
                lane="data",
            )
            registry["direct_reference"]["root_sha256"] = inventory["root_sha256"]
            registry["direct_reference"]["inventory_sha256"] = inventory[
                "inventory_sha256"
            ]
            registry["obligations"][0].update(inventory["entries"][0])
            registry["obligations"][0][
                "recipe_identity_sha256"
            ] = validator.obligation_recipe_identity(registry["obligations"][0])
            comparison = validator.compare_outputs(registry, direct, replay)
            self.assertEqual(comparison["status"], "FAIL")
            self.assertIn(
                "semantic_mismatch:qa/photon_et", comparison["failures"]
            )


if __name__ == "__main__":
    unittest.main()
