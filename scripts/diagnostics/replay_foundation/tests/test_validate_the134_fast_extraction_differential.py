#!/usr/bin/env python3
"""Focused tests for the THE-134 fast-extraction differential validator."""

from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np


HERE = Path(__file__).resolve()
VALIDATOR_PATH = HERE.parents[1] / "validate_the134_fast_extraction_differential.py"
SPEC = importlib.util.spec_from_file_location(
    "validate_the134_fast_extraction_differential", VALIDATOR_PATH
)
assert SPEC is not None and SPEC.loader is not None
validator = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(validator)


class FastExtractionDifferentialTest(unittest.TestCase):
    def test_nested_exact_comparison_accepts_nan_and_rejects_drift(self) -> None:
        left = np.asarray(
            [[1.0, float("nan")], [2.0, -3.0]], dtype=object
        )
        right = np.asarray(
            [[1.0, float("nan")], [2.0, -3.0]], dtype=object
        )
        self.assertTrue(validator.exact_value_equal(left, right))
        right[1, 1] = -3.000000000000001
        self.assertFalse(validator.exact_value_equal(left, right))

    def test_resource_gate_requires_speed_rss_and_scratch_gain(self) -> None:
        payload = {
            "schema": validator.RESOURCE_SCHEMA,
            "rows": {
                system: {
                    "baseline": {
                        "elapsed_seconds": 300.0,
                        "max_rss_kb": 1_200_000,
                        "analysis_output_bytes": 500_000_000,
                        "sidecar_bytes": 20_000,
                    },
                    "candidate": {
                        "elapsed_seconds": 150.0,
                        "max_rss_kb": 1_100_000,
                        "analysis_output_bytes": 25_000,
                        "sidecar_bytes": 20_000,
                    },
                }
                for system in validator.SYSTEMS
            },
        }
        report = validator.validate_resources(payload, minimum_speedup=1.5)
        self.assertEqual(report["status"], "PASS")

        payload["rows"]["auau"]["candidate"]["elapsed_seconds"] = 250.0
        report = validator.validate_resources(payload, minimum_speedup=1.5)
        self.assertEqual(report["status"], "FAIL")
        self.assertTrue(
            any(item.startswith("auau:speedup=") for item in report["failures"])
        )

    def test_resource_contract_rejects_missing_system(self) -> None:
        payload = {
            "schema": validator.RESOURCE_SCHEMA,
            "rows": {"pp": {}},
        }
        with self.assertRaisesRegex(
            validator.ValidationError, "system inventory"
        ):
            validator.validate_resources(payload, minimum_speedup=1.5)

    def test_sidecar_comparison_allows_only_pinned_code_metadata(self) -> None:
        baseline_code = "1" * 64
        candidate_code = "2" * 64
        arrays = {
            "event_id_hi": np.asarray([1, 2], dtype=np.uint64),
            "ordered_features": np.asarray(
                [[1.0, float("nan")], [2.0, 3.0]], dtype=object
            ),
            "weight_final": np.asarray([0.5, 1.0], dtype=np.float64),
        }
        metadata = {
            "tree_num_entries": 2,
            "schema_sha256": "3" * 64,
            "code_sha256": baseline_code,
        }
        calls = iter(
            [
                (arrays, sorted(arrays), metadata),
                (
                    {name: value.copy() for name, value in arrays.items()},
                    sorted(arrays),
                    {**metadata, "code_sha256": candidate_code},
                ),
            ]
        )
        preparer = SimpleNamespace(
            read_tree=lambda _path, _tree: next(calls)
        )
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            baseline = root / "baseline.root"
            candidate = root / "candidate.root"
            baseline.write_bytes(b"baseline")
            candidate.write_bytes(b"candidate")
            with patch.object(
                validator,
                "root_inventory",
                return_value={
                    validator.TREE_NAME: "TTree",
                    "code_sha256": "TNamed",
                },
            ):
                report = validator.compare_sidecars(
                    baseline,
                    candidate,
                    baseline_code_sha256=baseline_code,
                    candidate_code_sha256=candidate_code,
                    preparer=preparer,
                )
        self.assertEqual(report["status"], "PASS")

        mutated = {
            name: value.copy() for name, value in arrays.items()
        }
        mutated["weight_final"][1] = 1.000000000000001
        calls = iter(
            [
                (arrays, sorted(arrays), metadata),
                (
                    mutated,
                    sorted(arrays),
                    {**metadata, "code_sha256": candidate_code},
                ),
            ]
        )
        preparer = SimpleNamespace(
            read_tree=lambda _path, _tree: next(calls)
        )
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            baseline = root / "baseline.root"
            candidate = root / "candidate.root"
            baseline.write_bytes(b"baseline")
            candidate.write_bytes(b"candidate")
            with patch.object(
                validator,
                "root_inventory",
                return_value={validator.TREE_NAME: "TTree"},
            ):
                report = validator.compare_sidecars(
                    baseline,
                    candidate,
                    baseline_code_sha256=baseline_code,
                    candidate_code_sha256=candidate_code,
                    preparer=preparer,
                )
        self.assertEqual(report["status"], "FAIL")
        self.assertIn("branch:weight_final", report["failures"])


if __name__ == "__main__":
    unittest.main()
