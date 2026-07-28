#!/usr/bin/env python3
"""Focused unit tests for the THE-134 sidecar differential comparator."""

from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path

import numpy as np
import uproot


HERE = Path(__file__).resolve()
VALIDATOR_PATH = HERE.parents[1] / "validate_the134_sidecar_differential.py"
SPEC = importlib.util.spec_from_file_location(
    "validate_the134_sidecar_differential", VALIDATOR_PATH
)
assert SPEC is not None and SPEC.loader is not None
validator = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(validator)


class SidecarDifferentialComparatorTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)

    def tearDown(self) -> None:
        self.temporary.cleanup()

    @staticmethod
    def write_histogram(path: Path, values: np.ndarray) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        with uproot.recreate(path) as output:
            output["qa/hExact"] = (
                values.astype(np.float64),
                np.asarray([0.0, 1.0, 2.0, 3.0], dtype=np.float64),
            )

    def test_exact_histogram_neutrality_and_mutation_rejection(self) -> None:
        direct = self.root / "direct.root"
        writer = self.root / "writer.root"
        values = np.asarray([1.0, 4.0, 9.0], dtype=np.float64)
        self.write_histogram(direct, values)
        self.write_histogram(writer, values)
        report = validator.compare_histograms(direct, writer)
        self.assertEqual(report["status"], "PASS")
        self.assertEqual(report["histograms"], 1)

        self.write_histogram(writer, np.asarray([1.0, 4.0, 10.0]))
        report = validator.compare_histograms(direct, writer)
        self.assertEqual(report["status"], "FAIL")
        self.assertTrue(any("values" in item for item in report["failures"]))

    def test_source_identity_is_deterministic_and_row_specific(self) -> None:
        tag = "the134_fixture"
        pp = validator.source_sha256(tag, dict(validator.ROWS[0]))
        auau = validator.source_sha256(tag, dict(validator.ROWS[1]))
        self.assertRegex(pp, r"^[0-9a-f]{64}$")
        self.assertRegex(auau, r"^[0-9a-f]{64}$")
        self.assertNotEqual(pp, auau)
        self.assertEqual(
            pp, validator.source_sha256(tag, dict(validator.ROWS[0]))
        )

    def test_preflight_requires_frozen_direct_and_writer_profiles(self) -> None:
        receipt = self.root / "preflight_receipt.txt"
        exact_keys = ",".join(
            f"{arm}:{row['lane']}:{row['sample']}"
            for row in validator.ROWS
            for arm in ("direct", "writer")
        )
        fields = {
            "tag": "the134_fixture",
            "code_sha256": "1" * 64,
            "schema_sha": "2" * 64,
            "semantic_sha": "3" * 64,
            "only_keys": exact_keys,
            "extra_pp_template": "RJ_PP_PHOTONID_EXTRACT_ONLY=1",
            "extra_auau_template": "RJ_AUAU_BDT_EXTRACT_ONLY=1",
            "writer_extra_common_template": (
                "RJ_THE134_MULTIVIEW_SIDECAR_ONLY_V1=1"
            ),
            "writer_extra_pp_template": "",
            "writer_extra_auau_template": "",
        }
        receipt.write_text(
            "".join(f"{key}={value}\n" for key, value in fields.items()),
            encoding="utf-8",
        )
        parsed = validator.parse_preflight(receipt)
        self.assertEqual(parsed["fields"]["tag"], "the134_fixture")

        fields["extra_pp_template"] = ""
        receipt.write_text(
            "".join(f"{key}={value}\n" for key, value in fields.items()),
            encoding="utf-8",
        )
        with self.assertRaises(validator.ValidationError):
            validator.parse_preflight(receipt)

    def test_direct_sidecar_and_inventory_drift_rejected(self) -> None:
        row = dict(validator.ROWS[0])
        direct_base = (
            self.root / "direct" / row["lane"] / row["sample"]
        )
        writer_base = (
            self.root / "writer" / row["lane"] / row["sample"]
        )
        self.write_histogram(direct_base / "nested" / "analysis.root", np.ones(3))
        self.write_histogram(writer_base / "nested" / "analysis.root", np.ones(3))
        self.write_histogram(writer_base / validator.SIDECAR_NAME, np.ones(3))
        discovered = validator.discover_row_files(self.root, row)
        self.assertEqual(len(discovered["pairs"]), 1)

        self.write_histogram(direct_base / validator.SIDECAR_NAME, np.ones(3))
        with self.assertRaises(validator.ValidationError):
            validator.discover_row_files(self.root, row)


if __name__ == "__main__":
    unittest.main()
