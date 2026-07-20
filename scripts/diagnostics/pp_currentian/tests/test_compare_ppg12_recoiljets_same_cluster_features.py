#!/usr/bin/env python3
"""Pure contract tests for the same-candidate PPG12/RecoilJets oracle."""

from __future__ import annotations

import importlib.util
import csv
import math
import sys
import tempfile
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).resolve().parents[1] / "compare_ppg12_recoiljets_same_cluster_features.py"
SPEC = importlib.util.spec_from_file_location("same_cluster_oracle", MODULE_PATH)
assert SPEC and SPEC.loader
ORACLE = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = ORACLE
SPEC.loader.exec_module(ORACLE)


def valid_row() -> dict[str, float]:
    return {
        "cluster_prob": 0.5,
        "e11_over_e33": 0.9,
        "e32_over_e35": 0.9,
        "cluster_weta_cogx": 0.1,
        "cluster_wphi_cogx": 0.1,
        "cluster_et1": 0.8,
        "cluster_et2": 0.5,
        "cluster_et3": 0.5,
        "cluster_et4": 0.5,
        "npb_score": 0.8,
    }


class TestPpg12SelectionContract(unittest.TestCase):
    @staticmethod
    def identity_row(
        *,
        eta: float,
        phi: float,
        cluster_index: int,
        truth_track_id: int = 7,
    ) -> dict[str, float | int]:
        return {
            "segment": 0,
            "eventnumber": 1,
            "truth_track_id": truth_track_id,
            "cluster_Eta": eta,
            "cluster_Phi": phi,
            "cluster_index": cluster_index,
        }

    def test_model_route_uses_half_open_et_bins(self) -> None:
        cases = {
            7.999999: "base_E",
            8.0: "base_v3E",
            15.0: "base_v3E",
            34.999999: "base_v3E",
            35.0: "base_E",
        }
        for et, expected in cases.items():
            with self.subTest(et=et):
                self.assertEqual(ORACLE.selected_ppg12_model(et), expected)

    def test_slimtree_segment_follows_event_sequence_resets(self) -> None:
        segment = 0
        previous = None
        observed = []
        for eventnumber in (1, 1, 2, 4, 1, 1, 2, 3, 1):
            segment = ORACLE.next_ppg12_segment(segment, previous, eventnumber)
            observed.append(segment)
            previous = eventnumber
        self.assertEqual(observed, [0, 0, 0, 0, 1, 1, 1, 1, 2])

    def test_recoil_identity_uses_processed_ordinal_not_local_eventnumber(self) -> None:
        self.assertEqual(ORACLE.recoil_source_identity(1, 17, 1000), (0, 17))
        self.assertEqual(ORACLE.recoil_source_identity(1000, 998, 1000), (0, 998))
        self.assertEqual(ORACLE.recoil_source_identity(1001, 1, 1000), (1, 1))
        self.assertEqual(ORACLE.recoil_source_identity(4001, 1, 1000), (4, 1))
        with self.assertRaises(RuntimeError):
            ORACLE.recoil_source_identity(0, 1, 1000)

    def test_symmetric_scan_includes_events_absent_from_recoil(self) -> None:
        event = (4, 999)
        self.assertTrue(ORACLE.ppg12_event_in_scan(event, None))
        self.assertFalse(ORACLE.ppg12_event_in_scan(event, {(0, 1)}))

        ppg_only = self.identity_row(
            eta=0.1, phi=0.1, cluster_index=8, truth_track_id=91
        )
        ppg_only["segment"] = event[0]
        ppg_only["eventnumber"] = event[1]
        matches, rj_only, ppg12_only = ORACLE.match_rows(
            [], {event: [ppg_only]}
        )
        self.assertEqual(matches, [])
        self.assertEqual(rj_only, [])
        self.assertEqual(ppg12_only, [ppg_only])

    def test_selected_score_follows_unsmeared_route(self) -> None:
        scores = {"base_E": 0.2, "base_v3E": 0.8}
        self.assertEqual(ORACLE.selected_score(scores, 8.0), ("base_v3E", 0.8))
        self.assertEqual(ORACLE.selected_score(scores, 35.0), ("base_E", 0.2))

    def test_executable_trace_replaces_shadow_only_after_identity_checks(self) -> None:
        ppg = {
            "tree_entry": 4,
            "cluster_index": 2,
            "eventnumber": 9,
            "is_signal": 1,
            "cluster_Et": 20.0,
            "cluster_weta_cogx": 0.1,
            "cluster_wphi_cogx": 0.2,
            "vertexz": 1.0,
            "cluster_Eta": 0.3,
            "e11_over_e33": 0.9,
            "cluster_et1": 0.8,
            "cluster_et2": 0.2,
            "cluster_et3": 0.1,
            "cluster_et4": 0.05,
            "e32_over_e35": 0.95,
            "bdt_base_E": 0.4,
            "bdt_base_v3E": 0.8,
        }
        trace = {
            "tree_entry": "4", "chain_file_index": "0", "local_tree_entry": "4",
            "runnumber": "28", "eventnumber": "9", "cluster_index": "2",
            "truth_track_id": "7", "cluster_Et": "20", "raw_eiso": "0.1",
            "corrected_eiso": "0.22", "iso_threshold": "1.23",
            "noniso_threshold": "2.03", "sample_weight": "2", "mix_weight": "1",
            "lumi_weight": "1", "cross_weight": "2", "vertex_weight": "1",
            "truth_vertex_weight": "1", "trigger_weight": "1", "event_weight": "2",
            "weight": "2", "selected_model": "base_v3E", "base_E_score": "0.4",
            "base_v3E_score": "0.8", "selected_score": "0.8",
            "cluster_weta_cogx": "0.1", "cluster_wphi_cogx": "0.2",
            "vertexz": "1", "cluster_Eta": "0.3", "e11_over_e33": "0.9",
            "cluster_et1": "0.8", "cluster_et2": "0.2", "cluster_et3": "0.1",
            "cluster_et4": "0.05", "e32_over_e35": "0.95", "common_pass": "1",
            "tight": "1", "nontight": "0", "is_iso": "1", "is_noniso": "0",
            "logical_abcd_region": "1", "is_signal": "1", "truth_particle_index": "0",
            "truth_class": "1", "truth_pt": "20", "analysis_window_pass": "1",
            "signal_fill_A": "1", "signal_fill_B": "0", "signal_fill_C": "0",
            "signal_fill_D": "0", "fill_multiplicity": "1",
        }
        response = {
            "tree_entry": "4", "chain_file_index": "0", "local_tree_entry": "4",
            "runnumber": "28", "eventnumber": "9", "cluster_index": "2",
            "response_Et": "19.8", "response_window_pass": "1",
        }
        with tempfile.TemporaryDirectory() as tmp:
            trace_path = Path(tmp) / "trace.csv"
            response_path = Path(tmp) / "response.csv"
            for path, value in ((trace_path, trace), (response_path, response)):
                with path.open("w", newline="") as stream:
                    writer = csv.DictWriter(stream, fieldnames=list(value))
                    writer.writeheader()
                    writer.writerow(value)
            loaded = ORACLE.read_executable_trace(trace_path, response_path)
        ORACLE.bind_executable_trace({(0, 9): [ppg]}, loaded)
        self.assertEqual(ppg["tag_evidence_source"], ORACLE.PRESERVED_EXECUTABLE_EVIDENCE)
        self.assertEqual(ppg["selected_model"], "base_v3E")
        self.assertEqual(ppg["response_Et"], "19.8")
        self.assertEqual(ppg["weight_final"], 2.0)

    def test_infer_stored_model_route_requires_unique_score_match(self) -> None:
        self.assertEqual(ORACLE.infer_stored_model_route(0.2, 0.2, 0.8), "base_E")
        self.assertEqual(ORACLE.infer_stored_model_route(0.8, 0.2, 0.8), "base_v3E")
        self.assertEqual(ORACLE.infer_stored_model_route(0.2, 0.2, 0.2), "ambiguous")
        self.assertEqual(ORACLE.infer_stored_model_route(0.5, 0.2, 0.8), "unmatched")

    def test_report_bins_are_strict_open_intervals(self) -> None:
        self.assertEqual(ORACLE.find_bin(10.0), -1)
        self.assertEqual(ORACLE.find_bin(10.000001), 0)
        self.assertEqual(ORACLE.find_bin(12.0), -1)
        self.assertEqual(ORACLE.find_bin(11.999999), 0)

    def test_kinematic_boundaries_match_executable(self) -> None:
        self.assertTrue(ORACLE.passes_ppg12_kinematics(5.0, 0.0, 60.0))
        self.assertTrue(ORACLE.passes_ppg12_kinematics(5.0, 0.0, -60.0))
        self.assertFalse(ORACLE.passes_ppg12_kinematics(math.nextafter(5.0, 0.0), 0.0, 0.0))
        self.assertFalse(ORACLE.passes_ppg12_kinematics(5.0, -0.7, 0.0))
        self.assertFalse(ORACLE.passes_ppg12_kinematics(5.0, 0.7, 0.0))
        self.assertFalse(ORACLE.passes_ppg12_kinematics(5.0, 0.0, math.nextafter(60.0, math.inf)))

    def test_cpp_ratio_preserves_executable_zero_denominator_semantics(self) -> None:
        self.assertTrue(math.isinf(ORACLE.cpp_float_ratio(1.0, 0.0)))
        self.assertGreater(ORACLE.cpp_float_ratio(1.0, 0.0), 0.0)
        self.assertTrue(math.isnan(ORACLE.cpp_float_ratio(0.0, 0.0)))

    def test_tight_threshold_is_strict(self) -> None:
        row = valid_row()
        et = 20.0
        threshold, _, _ = ORACLE.ppg12_thresholds(et)
        tag_at_boundary, _ = ORACLE.classify(row, threshold, et)
        tag_above_boundary, _ = ORACLE.classify(row, math.nextafter(threshold, math.inf), et)
        self.assertNotEqual(tag_at_boundary, 1)
        self.assertEqual(tag_above_boundary, 1)

    def test_nontight_score_window_is_strict(self) -> None:
        row = valid_row()
        et = 20.0
        _, low, high = ORACLE.ppg12_thresholds(et)
        self.assertEqual(ORACLE.classify(row, low, et)[0], 3)
        self.assertEqual(ORACLE.classify(row, high, et)[0], 3)
        self.assertEqual(ORACLE.classify(row, 0.5 * (low + high), et)[0], 2)

    def test_nominal_mc_isolation_and_abcd_mapping(self) -> None:
        et = 20.0
        isolated = ORACLE.ppg12_iso_assignment(0.0, et)
        nonisolated = ORACLE.ppg12_iso_assignment(3.0, et)
        self.assertEqual(isolated["corrected_eiso"], 0.1)
        self.assertEqual(isolated["is_iso"], 1)
        self.assertEqual(nonisolated["is_noniso"], 1)
        self.assertEqual(ORACLE.abcd_region(1, 1, 0), "A")
        self.assertEqual(ORACLE.abcd_region(1, 0, 1), "B")
        self.assertEqual(ORACLE.abcd_region(2, 1, 0), "C")
        self.assertEqual(ORACLE.abcd_region(2, 0, 1), "D")
        self.assertEqual(ORACLE.abcd_region(3, 1, 0), "none")

    def test_signal_fill_flags_preserve_a_only_window_asymmetry(self) -> None:
        self.assertEqual(
            ORACLE.ppg12_signal_fill_flags("A", 0),
            {"A": 0, "B": 0, "C": 0, "D": 0},
        )
        self.assertEqual(
            ORACLE.ppg12_signal_fill_flags("A", 1),
            {"A": 1, "B": 0, "C": 0, "D": 0},
        )
        self.assertEqual(
            ORACLE.ppg12_signal_fill_flags("C", 0),
            {"A": 0, "B": 0, "C": 1, "D": 0},
        )

    def test_candidate_report_keeps_population_mismatch_rows(self) -> None:
        rj = {
            "segment": 0,
            "eventnumber": 1,
            "truth_track_id": 7,
            "cluster_Et": 20.0,
            "rj_selected_model": "base_v3E",
            "rj_recomputed_tag": 1,
            "rj_recomputed_common_pass": 1,
            "rj_recomputed_abcd_region": "A",
            "rj_stored_tag": 2,
            "rj_stored_common_pass": 1,
            "rj_stored_abcd_region": "C",
        }
        ppg = {
            "segment": 0,
            "eventnumber": 1,
            "truth_track_id": 7,
            "cluster_Et": 20.0,
            "selected_model": "base_v3E",
            "ppg12_recomputed_tag": 1,
            "ppg12_common_pass": 1,
            "ppg12_abcd_region": "A",
        }
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "candidates.csv"
            ORACLE.write_candidate_report(
                [(rj, ppg, 0.0)], [rj], [ppg], output,
                lane_id="photon:photon5:0mrad:si",
                runtime_contract_sha256="a" * 64,
            )
            with output.open(newline="") as handle:
                rows = list(csv.DictReader(handle))
        self.assertEqual([row["match_status"] for row in rows], ["matched", "rj_only", "ppg12_only"])
        self.assertTrue(all(row["lane_id"] == "photon:photon5:0mrad:si" for row in rows))
        self.assertTrue(all(row["runtime_contract_sha256"] == "a" * 64 for row in rows))
        self.assertEqual(rows[0]["stored_common_agree"], "1")
        self.assertEqual(rows[0]["stored_tag_agree"], "0")
        self.assertEqual(rows[0]["rj_stored_abcd_region"], "C")
        self.assertEqual(rows[0]["stored_abcd_agree"], "0")
        self.assertEqual(
            rows[0]["ppg12_tag_evidence_source"],
            ORACLE.UNAVAILABLE_EXECUTABLE_EVIDENCE,
        )

    def test_same_identity_with_large_dr_is_matched_and_flagged(self) -> None:
        rj = self.identity_row(eta=0.10, phi=0.10, cluster_index=1)
        ppg = self.identity_row(eta=0.20, phi=0.20, cluster_index=2)
        matches, rj_only, ppg12_only = ORACLE.match_rows(
            [rj], {(0, 1): [ppg]}
        )
        self.assertEqual(len(matches), 1)
        self.assertGreater(matches[0][2], ORACLE.GEOMETRY_DIAGNOSTIC_DR)
        self.assertEqual(rj_only, [])
        self.assertEqual(ppg12_only, [])

        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp) / "candidates.csv"
            ORACLE.write_candidate_report(matches, rj_only, ppg12_only, output)
            with output.open(newline="") as handle:
                row = next(csv.DictReader(handle))
        self.assertEqual(row["identity_match"], "1")
        self.assertEqual(row["match_dr_exceeds_002"], "1")

    def test_same_identity_multiplicity_surplus_remains_population_mismatch(self) -> None:
        rj_rows = [
            self.identity_row(eta=0.10, phi=0.10, cluster_index=1),
            self.identity_row(eta=0.11, phi=0.11, cluster_index=2),
        ]
        ppg_rows = [self.identity_row(eta=0.20, phi=0.20, cluster_index=3)]
        matches, rj_only, ppg12_only = ORACLE.match_rows(
            rj_rows, {(0, 1): ppg_rows}
        )
        self.assertEqual(len(matches), 1)
        self.assertEqual(len(rj_only), 1)
        self.assertEqual(ppg12_only, [])

    def test_population_rejection_reports_first_failed_gate(self) -> None:
        base = {
            "event_found": 1,
            "vertex_pass": 1,
            "truth_signal": 1,
            "truth_cluster_count": 1,
            "unmasked_cluster_count": 1,
            "et_pass_count": 1,
            "eta_pass_count": 1,
            "accepted_count": 1,
        }
        cases = [
            ({}, "ppg12_event_missing"),
            ({**base, "vertex_pass": 0}, "ppg12_vertex_rejected"),
            ({**base, "truth_signal": 0}, "ppg12_truth_signal_rejected"),
            ({**base, "truth_cluster_count": 0}, "ppg12_truth_matched_cluster_missing"),
            ({**base, "unmasked_cluster_count": 0}, "ppg12_tower_mask_rejected"),
            ({**base, "et_pass_count": 0}, "ppg12_reco_et_rejected"),
            ({**base, "eta_pass_count": 0}, "ppg12_reco_eta_rejected"),
            ({**base, "accepted_count": 0}, "ppg12_population_unknown"),
            (base, "recoiljets_multiplicity_surplus"),
        ]
        for audit, expected in cases:
            with self.subTest(expected=expected):
                self.assertEqual(ORACLE.first_population_rejection(audit), expected)

    def test_population_annotation_is_written_to_candidate_report(self) -> None:
        rj = self.identity_row(eta=0.1, phi=0.1, cluster_index=1)
        ppg = self.identity_row(eta=0.1, phi=0.1, cluster_index=2)
        audit = {
            (0, 1, 7): {
                "event_found": 1,
                "vertex_pass": 1,
                "ppg12_vertexz": 0.5,
                "truth_signal": 1,
                "truth_cluster_count": 1,
                "tower_masked_count": 1,
                "unmasked_cluster_count": 0,
                "et_pass_count": 0,
                "eta_pass_count": 0,
                "accepted_count": 0,
            }
        }
        ORACLE.annotate_population_mismatches([], [rj], [ppg], audit)
        self.assertEqual(rj["population_mismatch_reason"], "ppg12_tower_mask_rejected")
        self.assertEqual(ppg["population_mismatch_reason"], "recoiljets_multiplicity_deficit")


if __name__ == "__main__":
    unittest.main()
