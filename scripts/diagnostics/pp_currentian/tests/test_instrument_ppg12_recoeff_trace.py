#!/usr/bin/env python3
"""Contract tests for the source-locked RecoEff trace transform."""

from __future__ import annotations

import importlib.util
import subprocess
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
MODULE_PATH = Path(__file__).resolve().parents[1] / "instrument_ppg12_recoeff_trace.py"
SPEC = importlib.util.spec_from_file_location("recoeff_trace_transform", MODULE_PATH)
assert SPEC and SPEC.loader
TRACE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = TRACE
SPEC.loader.exec_module(TRACE)
REVISION = "29f8223bd9b36dffab07961b597afa94185bbdf1"


def canonical_source() -> str:
    return subprocess.check_output(
        [
            "git",
            "-c",
            f"safe.directory={ROOT / 'ppg12codeGit'}",
            "-C",
            str(ROOT / "ppg12codeGit"),
            "show",
            f"{REVISION}:efficiencytool/RecoEffCalculator_TTreeReader.C",
        ],
        text=True,
    )


class TestRecoEffTraceTransform(unittest.TestCase):
    def test_exact_29f_transform_has_all_candidate_and_weight_evidence(self) -> None:
        source = canonical_source()
        transformed, operations = TRACE.instrument(source)
        self.assertEqual(
            [item["label"] for item in operations],
            [
                "event_identity_reader",
                "trace_stream_setup",
                "capture_raw_isolation",
                "truth_vertex_factor_scope",
                "capture_truth_vertex_factor",
                "trace_model_scope",
                "reuse_scoped_model",
                "reuse_scoped_score",
                "candidate_trace_emit",
                "response_trace_emit",
            ],
        )
        for field in (
            "tree_entry,chain_file_index,local_tree_entry",
            "sample_weight,mix_weight,",
            "lumi_weight,cross_weight,vertex_weight,truth_vertex_weight,",
            "trigger_weight,event_weight,weight",
            "base_E_score,base_v3E_score,selected_score",
            "analysis_window_pass,signal_fill_A,signal_fill_B,signal_fill_C",
        ):
            self.assertIn(field, transformed)
        self.assertIn("oracle_trace\n                    << ientry", transformed)
        self.assertIn("if (tight && iso)\n            {", transformed)
        self.assertEqual(
            source.count('h_tight_iso_cluster_signal[etabin]->Fill'),
            transformed.count('h_tight_iso_cluster_signal[etabin]->Fill'),
        )

    def test_marker_drift_fails_closed(self) -> None:
        drifted = canonical_source().replace(
            '    TTreeReaderValue<int> runnumber(reader, "runnumber");',
            '    TTreeReaderValue<int> runnumber_changed(reader, "runnumber");',
            1,
        )
        with self.assertRaises(TRACE.TransformFailure):
            TRACE.instrument(drifted)


if __name__ == "__main__":
    unittest.main()
