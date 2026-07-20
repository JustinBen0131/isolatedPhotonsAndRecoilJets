#!/usr/bin/env python3
"""Schema regression tests for the THE-112 canary validator."""

from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).resolve().parents[1] / "validate_the112_sideband_canary.py"
SPEC = importlib.util.spec_from_file_location("validate_the112_sideband_canary", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
validator = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = validator
SPEC.loader.exec_module(validator)


def complete_keys(sample_kind: str) -> list[str]:
    categories = (
        ("all",)
        if sample_kind == "data"
        else ("all", "truthPrompt", "truthSignal", "truthBackground")
    )
    return [
        f"ALL/{validator.expected_surface_basename(category, centrality, view)}"
        for category in categories
        for centrality in validator.CENTRALITY_TOKENS
        for view in validator.ISO_VIEWS
    ]


class CanaryValidatorSchemaTest(unittest.TestCase):
    def test_simulation_inventory_requires_truth_background(self) -> None:
        keys = complete_keys("inclusive")
        failures: list[str] = []
        report = validator.validate_v2_inventory(keys, "inclusive", failures)
        self.assertFalse(failures)
        self.assertIn("truthBackground", report["categories"])

        missing_background = [key for key in keys if "truthBackground" not in key]
        failures = []
        validator.validate_v2_inventory(missing_background, "inclusive", failures)
        self.assertTrue(any("truthBackground" in failure for failure in failures))

    def test_data_inventory_remains_audience_neutral_all_only(self) -> None:
        failures: list[str] = []
        report = validator.validate_v2_inventory(complete_keys("data"), "data", failures)
        self.assertFalse(failures)
        self.assertEqual(report["categories"], ["all"])


if __name__ == "__main__":
    unittest.main()
