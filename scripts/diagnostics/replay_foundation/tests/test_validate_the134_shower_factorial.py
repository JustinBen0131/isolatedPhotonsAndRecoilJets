#!/usr/bin/env python3
"""Pure mutation tests for the THE-134 independent shower oracle."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[1]))

import validate_the134_shower_factorial as validator  # noqa: E402


def raw_cell(eta: int, phi: int, energy: float) -> dict[str, object]:
    return {
        "tower_eta_index": eta,
        "tower_phi_index": phi,
        "rawcluster_owned": 1,
        "rawcluster_value_present": 1,
        "rawcluster_map_value": energy,
        "calibrated_energy": energy,
        "is_good": 1,
        "is_zero": int(energy == 0.0),
        "is_negative": int(energy < 0.0),
        "is_nonfinite": int(not np.isfinite(energy)),
        "grid_membership_bitmask": 3,
    }


class RawClusterNativeShapeReplayTest(unittest.TestCase):
    def setUp(self) -> None:
        self.cells = [
            raw_cell(9, 20, 0.05),
            raw_cell(10, 20, 2.0),
            raw_cell(10, 21, 1.0),
            raw_cell(11, 20, 0.5),
            raw_cell(11, 21, 0.25),
        ]

    def test_exact_float32_rawclusterv1_native_shape(self) -> None:
        rebuilt = validator.native_rawcluster_shape_from_cells(0.070, self.cells)
        self.assertTrue(rebuilt["valid"])
        self.assertEqual(rebuilt["center_eta_index"], 10)
        self.assertEqual(rebuilt["center_phi_index"], 20)
        expected = np.asarray([1.0, 0.6, 1.0 / 3.0, 1.0 / 15.0], dtype=np.float32)
        observed = np.asarray(
            [rebuilt[f"native_et{index}"] for index in range(1, 5)],
            dtype=np.float32,
        )
        np.testing.assert_array_equal(observed, expected)

        zero_floor = validator.native_rawcluster_shape_from_cells(0.0, self.cells)
        self.assertTrue(zero_floor["valid"])
        self.assertNotEqual(zero_floor["native_et1"], rebuilt["native_et1"])

    def test_coordinated_native_and_ordered_feature_drift_is_rejected(self) -> None:
        rebuilt = validator.view_from_cells("H70", self.cells)
        stored = {
            f"native_et{index}": float(rebuilt[f"native_et{index}"])
            for index in range(1, 5)
        }
        expected_features = np.asarray(
            [
                20.0,
                rebuilt["weta_cogx"],
                rebuilt["wphi_cogx"],
                0.0,
                0.1,
                rebuilt["e11_over_e33"],
                *(rebuilt[f"native_et{index}"] for index in range(1, 5)),
                rebuilt["e32_over_e35"],
            ],
            dtype=np.float32,
        )
        self.assertFalse(validator.native_et_replay_mismatch(stored, rebuilt, 0.0))

        # This paired mutation was invisible when the validator trusted both
        # the stored native value and the ordered vector built from it.
        stored["native_et1"] = float(np.float32(stored["native_et1"] + 0.25))
        coordinated_features = expected_features.copy()
        coordinated_features[6] = np.float32(stored["native_et1"])
        self.assertTrue(validator.native_et_replay_mismatch(stored, rebuilt, 0.0))
        self.assertFalse(
            np.allclose(
                coordinated_features,
                expected_features,
                rtol=0.0,
                atol=0.0,
                equal_nan=True,
            )
        )

    def test_outside_grid_owned_provenance_is_valid_but_unowned_zero_is_not(self) -> None:
        outside_owned = raw_cell(18, 27, 0.2)
        outside_owned["grid_membership_bitmask"] = 0
        self.assertTrue(validator.valid_grid_or_owned_provenance(outside_owned))

        unowned = dict(outside_owned)
        unowned["rawcluster_owned"] = 0
        unowned["rawcluster_value_present"] = 0
        self.assertFalse(validator.valid_grid_or_owned_provenance(unowned))

        malformed = dict(outside_owned)
        malformed["grid_membership_bitmask"] = 4
        self.assertFalse(validator.valid_grid_or_owned_provenance(malformed))


if __name__ == "__main__":
    unittest.main()
