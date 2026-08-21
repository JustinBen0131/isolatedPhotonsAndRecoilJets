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
    @staticmethod
    def config_payload(cone_r: str, internal_views: str | None = None) -> dict[str, str]:
        views = internal_views or validator.EXPECTED_INTERNAL_ISO_VIEWS
        return {
            "analysis_config_yaml": (
                "isSlidingIso: true\n"
                "isSlidingAndFixed: false\n"
                "fixedGeV: 0.0\n"
                f"coneR: {cone_r}\n"
                f"internal_iso_cone_views: '{views}'\n"
                "isolation_wp:\n"
                "  truthIsoGeV: 4.0\n"
            )
        }

    def test_stamped_nominal_scalar_cone_requires_exact_internal_views(self) -> None:
        failures: list[str] = []
        report = validator.validate_embedded_config(
            self.config_payload("0.4"), "stamped", failures
        )
        self.assertFalse(failures)
        self.assertEqual(report["coneR"], 0.4)

        failures = []
        validator.validate_embedded_config(
            self.config_payload("0.4", "isoR40_isSliding:0.40:true:0.0"),
            "missing-r30",
            failures,
        )
        self.assertTrue(any("internal_iso_cone_views" in failure for failure in failures))
        self.assertTrue(any("coneR" in failure for failure in failures))

    def test_source_ordered_cone_list_remains_valid(self) -> None:
        failures: list[str] = []
        validator.validate_embedded_config(
            self.config_payload("[0.40, 0.30]"), "source", failures
        )
        self.assertFalse(failures)

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
