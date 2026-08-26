from __future__ import annotations

import hashlib
import json
from pathlib import Path
import unittest

import numpy as np

from photonjet.training.features import AUAU_H70_FEATURES, PP_H70_FEATURES
from photonjet.training.identity import candidate_identity, event_identity
from photonjet.training.split import assign_group_split
from photonjet.training.weights import inverse_density_weights, normalize_mean_one


ROOT = Path(__file__).resolve().parents[1]


class TrainingContractsTest(unittest.TestCase):
    def test_feature_order_and_cardinality(self) -> None:
        self.assertEqual(len(PP_H70_FEATURES), 11)
        self.assertEqual(len(AUAU_H70_FEATURES), 14)
        self.assertEqual(PP_H70_FEATURES[0], "cluster_Et")
        self.assertEqual(AUAU_H70_FEATURES[-1], "centrality")

    def test_source_qualified_identities(self) -> None:
        event = event_identity("sample", 2, 42, 7)
        candidate = candidate_identity("sample", 2, 42, 7, 9)
        self.assertEqual(len(event), 64)
        self.assertEqual(len(candidate), 64)
        self.assertNotEqual(event, candidate)
        self.assertNotEqual(event, event_identity("sample", 3, 42, 7))

    def test_split_is_group_stable(self) -> None:
        groups = [f"event-{index}" for index in range(100)]
        first = assign_group_split(groups, seed=131)
        second = assign_group_split(groups, seed=131)
        self.assertEqual(first, second)
        self.assertEqual(len(first), len(groups))
        self.assertTrue({"train", "validation", "holdout"}.issubset(set(first)))

    def test_weight_normalization(self) -> None:
        normalized = normalize_mean_one([1.0, 2.0, 3.0])
        self.assertAlmostEqual(float(normalized.mean()), 1.0)
        weights = inverse_density_weights(
            [0.1, 0.2, 0.8, 0.9],
            [0, 0, 1, 1],
            bins=[0.0, 0.5, 1.0],
        )
        self.assertAlmostEqual(float(weights.mean()), 1.0)
        self.assertTrue(np.all(weights > 0))

    def test_packaged_models_are_explicit_hash_only_references(self) -> None:
        for system in ("pp", "auau"):
            manifest = json.loads(
                (ROOT / f"models/manifests/{system}_h70.json").read_text(encoding="utf-8")
            )
            artifact = ROOT / manifest["artifact_filename"]
            observed = hashlib.sha256(artifact.read_bytes()).hexdigest()
            self.assertEqual(manifest["schema"], "PhotonJetModelReferenceV1")
            self.assertEqual(manifest["status"], "HASH_BOUND_REFERENCE_ONLY")
            self.assertEqual(manifest["usage"], "not_consumed_by_public_v1")
            self.assertEqual(manifest["inference_validation"], "NOT_RUN")
            self.assertEqual(manifest["artifact_sha256"], observed)


if __name__ == "__main__":
    unittest.main()
