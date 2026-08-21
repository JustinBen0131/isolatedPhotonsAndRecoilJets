#!/usr/bin/env python3
"""Focused regression tests for the THE-107 Au+Au label construction."""

from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd


TRAINER = Path(__file__).resolve().parents[1] / "train_auau_photon_bdt.py"
SPEC = importlib.util.spec_from_file_location("the107_trainer", TRAINER)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def base_frame() -> pd.DataFrame:
    rows = []
    cases = [
        ("run28_embeddedPhoton12", 1, 12, 1, 1, 1, 1),
        ("run28_embeddedPhoton20", 1, 20, 2, 0, 0, 1),
        ("run28_embeddedPhoton12", 1, 12, 3, 0, 0, -1),
        ("run28_embeddedPhoton20", 1, 20, 0, 0, 0, -1),
        ("run28_embeddedJet12", 2, 12, 1, 1, 1, -1),
        ("run28_embeddedJet20", 2, 20, 2, 0, 0, -1),
        ("run28_embeddedJet30", 2, 30, 3, 0, 0, 0),
        ("run28_embeddedJet40", 2, 40, 0, 0, 0, 0),
    ]
    for source, role, code, photon_class, iso_pass, nominal, ppg12 in cases:
        rows.append(
            {
                "source_sample": source,
                "source_role": role,
                "source_sample_code": code,
                "truth_match_found": int(photon_class != 0),
                "truth_photon_class": photon_class,
                "truth_is_prompt": int(photon_class in (1, 2)),
                "truth_iso_et": 2.0 if iso_pass else 7.0,
                "truth_iso_pass": iso_pass,
                "is_signal": nominal,
                "ppg12_source_role_label": ppg12,
            }
        )
    return pd.DataFrame(rows)


class LabelContractTest(unittest.TestCase):
    def test_legacy_contract_does_not_require_the107_fields(self):
        frame = pd.DataFrame({"is_signal": [0, 1, 0]})
        result, report = MODULE.apply_training_label_contract(
            frame, "extracted-is-signal"
        )
        self.assertIs(result, frame)
        self.assertEqual(report["rows_kept"], 3)
        self.assertEqual(report["signal_rows_kept"], 1)

    def test_nominal_preserves_all_rows_and_labels(self):
        frame = base_frame()
        result, report = MODULE.apply_training_label_contract(
            frame, "nominal-isolated-prompt"
        )
        self.assertEqual(len(result), len(frame))
        self.assertEqual(result["is_signal"].tolist(), frame["is_signal"].tolist())
        self.assertEqual(report["rows_discarded"], 0)
        self.assertEqual(report["nominal_label_closure_mismatches"], 0)

    def test_ppg12_source_role_ignores_truth_isolation(self):
        frame = base_frame()
        frame[MODULE.PPG12_EXACT_WEIGHT_COLUMN] = 3.0
        result, report = MODULE.apply_training_label_contract(frame, "ppg12-source-role")
        self.assertEqual(len(result), 4)
        self.assertEqual(result["is_signal"].tolist(), [1, 1, 0, 0])
        self.assertEqual(report["rows_discarded"], 4)
        self.assertTrue(report["discarded_stale_precomputed_weight"])
        self.assertNotIn(MODULE.PPG12_EXACT_WEIGHT_COLUMN, result.columns)
        # The second retained photon-source row is prompt but non-isolated.
        self.assertEqual(result.iloc[1]["truth_iso_pass"], 0)
        self.assertEqual(result.iloc[1]["is_signal"], 1)

    def test_extracted_contract_mismatch_is_fatal(self):
        frame = base_frame()
        frame.loc[0, "ppg12_source_role_label"] = 0
        with self.assertRaises(SystemExit):
            MODULE.apply_training_label_contract(frame, "ppg12-source-role")

    def test_adding_row_indices_preserves_sklearn_split(self):
        try:
            from sklearn.model_selection import train_test_split
        except ImportError:
            self.skipTest("scikit-learn is not installed in the lightweight local test environment")
        x = np.arange(80, dtype="float32").reshape(40, 2)
        y = np.asarray([0, 1] * 20, dtype="int32")
        w = np.linspace(1.0, 2.0, len(y))
        old = train_test_split(x, y, w, test_size=0.2, random_state=13, stratify=y)
        new = train_test_split(
            x,
            y,
            w,
            np.arange(len(y)),
            test_size=0.2,
            random_state=13,
            stratify=y,
        )
        for old_array, new_array in zip(old, new[:6]):
            np.testing.assert_array_equal(old_array, new_array)

    def test_exact_holdout_cache_retains_label_audit_fields(self):
        frame = base_frame()
        frame["centrality"] = np.arange(len(frame), dtype="float64")
        frame["cluster_Et"] = 20.0
        frame["cluster_Eta"] = 0.1
        frame["input_file_index"] = 0
        frame["input_tree_entry"] = np.arange(len(frame), dtype="int64")
        frame["nominal_is_signal"] = frame["is_signal"]
        indices = np.asarray([1, 6], dtype="int64")
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "holdout.npz"
            MODULE.write_holdout_score_cache(
                path,
                frame,
                indices,
                np.asarray([[1.0], [2.0]], dtype="float32"),
                frame.iloc[indices]["is_signal"].to_numpy(),
                np.ones(2),
                np.asarray([0.8, 0.2]),
                ["cluster_Et"],
            )
            with np.load(path, allow_pickle=True) as payload:
                self.assertEqual(payload["truth_iso_pass"].tolist(), [0, 0])
                self.assertEqual(payload["ppg12_source_role_label"].tolist(), [1, 0])
                self.assertEqual(payload["features"].tolist(), ["cluster_Et"])


if __name__ == "__main__":
    unittest.main()
