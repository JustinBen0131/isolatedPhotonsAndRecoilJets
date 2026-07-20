#!/usr/bin/env python3
"""Stage-order tests for paired-oracle first-divergence evidence."""

from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).resolve().parents[1] / "audit_ppg12_recoiljets_paired_oracle.py"
SPEC = importlib.util.spec_from_file_location("paired_oracle_audit", MODULE_PATH)
assert SPEC and SPEC.loader
AUDIT = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = AUDIT
SPEC.loader.exec_module(AUDIT)


def passing_row() -> dict[str, str]:
    row = {
        "match_status": "matched",
        "candidate_identity": "seg0:evt1:trk7:rj2:ppg3",
        "identity_match": "1",
        "model_route_agree": "1",
        "rj_selected_model": "base_v3E",
        "rj_inferred_stored_model": "base_v3E",
        "ppg12_selected_model": "base_v3E",
        "ppg12_tag_evidence_source": AUDIT.PRESERVED_EXECUTABLE_EVIDENCE,
        "ppg12_isolation_abcd_evidence_source": AUDIT.PRESERVED_EXECUTABLE_EVIDENCE,
        "ppg12_truth_response_fill_evidence_source": AUDIT.PRESERVED_EXECUTABLE_EVIDENCE,
        "ppg12_weight_evidence_source": AUDIT.PRESERVED_EXECUTABLE_EVIDENCE,
        "base_E_score_delta": "0",
        "base_v3E_score_delta": "0",
        "selected_score_delta": "0",
        "stored_minus_routed_score": "0",
        "rj_stored_bdt_score": "0.8",
        "ppg12_selected_bdt_score": "0.8",
        "common_agree": "1",
        "tag_agree": "1",
        "stored_common_agree": "1",
        "stored_tag_agree": "1",
        "signal_status_agree": "1",
        "rj_is_iso": "1",
        "ppg12_is_iso": "1",
        "rj_is_noniso": "0",
        "ppg12_is_noniso": "0",
        "rj_raw_eiso": "0.2",
        "ppg12_raw_eiso": "0.2",
        "rj_corrected_eiso": "0.34",
        "ppg12_corrected_eiso": "0.34",
        "rj_iso_threshold": "1.23",
        "ppg12_iso_threshold": "1.23",
        "rj_noniso_threshold": "2.03",
        "ppg12_noniso_threshold": "2.03",
        "abcd_agree": "1",
        "stored_abcd_agree": "1",
        "truth_class_agree": "1",
        "rj_truth_class": "1",
        "truth_class": "1",
        "rj_logical_abcd_region": "1",
        "ppg12_logical_abcd_region": "1",
        "rj_analysis_window_pass": "1",
        "ppg12_analysis_window_pass": "1",
        "rj_response_Et": "20",
        "ppg12_response_Et": "20",
        "truth_pt": "20",
        "rj_response_window_pass": "1",
        "ppg12_response_window_pass": "1",
        "ppg12_is_signal": "1",
        "rj_signal_fill_A": "1",
        "rj_signal_fill_B": "0",
        "rj_signal_fill_C": "0",
        "rj_signal_fill_D": "0",
        "ppg12_signal_fill_A": "1",
        "ppg12_signal_fill_B": "0",
        "ppg12_signal_fill_C": "0",
        "ppg12_signal_fill_D": "0",
        "rj_signal_fill_multiplicity": "1",
        "ppg12_fill_multiplicity": "1",
        "rj_weight_lane_code": "1",
        "rj_weight_component_code": "1",
        "rj_weight_slice": "2",
        "ppg12_weight_sample": "2",
        "rj_weight_mix": "1",
        "ppg12_weight_mix": "1",
        "rj_weight_period": "1",
        "ppg12_weight_lumi": "1",
        "ppg12_weight_cross": "2",
        "rj_weight_vertex": "1",
        "ppg12_weight_vertex": "1",
        "ppg12_weight_truth_vertex": "1",
        "ppg12_weight_trigger": "1",
        "ppg12_weight_event": "2",
        "rj_weight_final": "2",
        "ppg12_weight_final": "2",
        "rj_weight_product_delta": "0",
        "rj_event_weight_delta": "0",
    }
    for feature in AUDIT.ORACLE_FEATURE_NAMES:
        row[f"{feature}_rj"] = "1"
        row[f"{feature}_ppg12"] = "1"
        row[f"{feature}_delta"] = "0"
    return row


class TestFirstDivergenceAudit(unittest.TestCase):
    def test_period_labels_map_to_preserved_config_suffixes(self) -> None:
        self.assertEqual(
            AUDIT.PERIOD_CONFIG_VAR_SUFFIX,
            {"0mrad": "0rad", "1p5mrad": "1p5mrad"},
        )

    def test_yaml_cpp_header_tree_drift_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            include_root = root / "include"
            tree = include_root / "yaml-cpp"
            tree.mkdir(parents=True)
            (tree / "yaml.h").write_text("// sealed yaml umbrella\n")
            (tree / "node.h").write_text("// sealed yaml node\n")
            rows = []
            tree_digest = hashlib.sha256()
            for path in sorted(tree.iterdir()):
                relative = path.name
                digest = hashlib.sha256(path.read_bytes()).hexdigest()
                rows.append({"relative_path": relative, "sha256": digest})
                tree_digest.update(relative.encode())
                tree_digest.update(b"\0")
                tree_digest.update(bytes.fromhex(digest))
            receipt = root / "yaml_cpp_header_tree_receipt.json"
            receipt.write_text(
                json.dumps(
                    {
                        "schema_version": 1,
                        "role": "ppg_recoeff_yaml_cpp_header_tree",
                        "include_root": str(include_root),
                        "staged_tree": str(tree),
                        "tree_sha256": tree_digest.hexdigest(),
                        "files": rows,
                    }
                )
            )
            self.assertEqual(
                AUDIT.validate_yaml_cpp_header_receipt(receipt),
                include_root.resolve(),
            )
            (tree / "yaml.h").write_text("// drifted yaml umbrella\n")
            with self.assertRaisesRegex(AUDIT.AuditFailure, "header drifted"):
                AUDIT.validate_yaml_cpp_header_receipt(receipt)

    def analyze(
        self, row: dict[str, str], *, expected_lane_id: str | None = None
    ) -> dict:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "candidates.csv"
            with path.open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=list(row))
                writer.writeheader()
                writer.writerow(row)
            return AUDIT.analyze_candidate_report(
                path, expected_lane_id=expected_lane_id
            )

    def test_weight_component_code_is_lane_aware_for_si_and_di(self) -> None:
        si_lane = "photon:photon5:0mrad:si"
        si = passing_row()
        si["lane_id"] = si_lane
        self.assertEqual(
            self.analyze(si, expected_lane_id=si_lane)["status"],
            "PASS",
        )
        di_lane = "photon:photon5:0mrad:di"
        di = passing_row()
        di["lane_id"] = di_lane
        di["rj_weight_component_code"] = "2"
        self.assertEqual(
            self.analyze(di, expected_lane_id=di_lane)["status"],
            "PASS",
        )
        wrong = passing_row()
        wrong["lane_id"] = di_lane
        report = self.analyze(wrong, expected_lane_id=di_lane)
        self.assertEqual(report["status"], "FAIL")
        self.assertEqual(report["first_divergence"]["field"], "rj_weight_component_code")

    def test_passing_row_reaches_every_stage(self) -> None:
        report = self.analyze(passing_row())
        self.assertEqual(report["status"], "PASS")
        self.assertIsNone(report["first_divergence"])
        self.assertEqual(
            list(report["stage_counts"]), list(AUDIT.FIRST_DIVERGENCE_STAGE_ORDER)
        )
        self.assertTrue(all(value["evaluated"] == 1 for value in report["stage_counts"].values()))

    def test_each_mutation_fails_at_its_stage(self) -> None:
        cases = (
            ("population", "match_status", "rj_only"),
            ("features", "cluster_Et_score_input_delta", "0.01"),
            ("scores", "selected_score_delta", "0.01"),
            ("route", "model_route_agree", "0"),
            ("tags", "tag_agree", "0"),
            ("isolation_abcd", "rj_is_iso", "0"),
            ("truth_response_fills", "rj_signal_fill_A", "0"),
            ("weights", "rj_event_weight_delta", "0.01"),
        )
        for expected_stage, field, value in cases:
            with self.subTest(stage=expected_stage):
                row = passing_row()
                row[field] = value
                report = self.analyze(row)
                self.assertEqual(report["status"], "FAIL")
                self.assertEqual(report["first_divergence"]["stage"], expected_stage)

    def test_stage_order_reports_earliest_not_largest(self) -> None:
        row = passing_row()
        row["cluster_Et_score_input_delta"] = "0.01"
        row["rj_event_weight_delta"] = "100"
        report = self.analyze(row)
        self.assertEqual(report["first_divergence"]["stage"], "features")

    def test_python_shadow_cannot_certify_executable_stages(self) -> None:
        cases = (
            ("tags", "ppg12_tag_evidence_source"),
            ("isolation_abcd", "ppg12_isolation_abcd_evidence_source"),
            ("truth_response_fills", "ppg12_truth_response_fill_evidence_source"),
            ("weights", "ppg12_weight_evidence_source"),
        )
        for expected_stage, field in cases:
            with self.subTest(stage=expected_stage):
                row = passing_row()
                row[field] = "python_shadow_not_executable"
                report = self.analyze(row)
                self.assertEqual(report["status"], "FAIL")
                self.assertEqual(report["first_divergence"]["stage"], expected_stage)
                self.assertEqual(report["first_divergence"]["field"], field)

    def test_asymmetric_executable_abcd_evidence_fails(self) -> None:
        row = passing_row()
        row["ppg12_signal_fill_C"] = "1"
        report = self.analyze(row)
        self.assertEqual(report["status"], "FAIL")
        self.assertEqual(report["first_divergence"]["stage"], "truth_response_fills")
        self.assertEqual(report["first_divergence"]["field"], "rj_signal_fill_C")

    def test_asymmetric_executable_weight_evidence_fails(self) -> None:
        row = passing_row()
        row["ppg12_weight_final"] = "2.25"
        report = self.analyze(row)
        self.assertEqual(report["status"], "FAIL")
        self.assertEqual(report["first_divergence"]["stage"], "weights")
        self.assertEqual(report["first_divergence"]["field"], "rj_weight_final")

    def test_raw_reconstruction_dependency_mutation_fails(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            paths = {}
            for role in AUDIT.RAW_RECONSTRUCTION_FILE_ROLES:
                path = root / role
                path.write_text(f"frozen raw dependency: {role}\n")
                paths[role] = str(path)
            lane = {
                "lane_id": "photon:photon5:1p5mrad:si",
                "sample": "Photon5",
                "period": "1p5mrad",
                "interaction": "SI",
                "rows": 5,
            }
            rng = AUDIT.historical_rng_contract()
            runtime = {
                "profile": AUDIT.EXPECTED_RUNTIME_PROFILE,
                "offline_main": AUDIT.EXPECTED_OFFLINE_MAIN,
                "actual_offline_main": AUDIT.EXPECTED_OFFLINE_MAIN,
            }
            source_pairs = {
                "schema_version": 1,
                "sample": "Photon5",
                "interaction": "SI",
                "row_count": 5,
            }
            receipt = AUDIT.raw_reconstruction_dependency_receipt(
                paths, lane, rng, runtime, source_pairs
            )
            AUDIT.validate_raw_reconstruction_dependency_receipt(
                receipt, paths, lane, rng, runtime, source_pairs
            )
            Path(paths["ppg_macro"]).write_text("mutated reconstruction macro\n")
            with self.assertRaisesRegex(
                AUDIT.AuditFailure, "raw reconstruction dependency receipt differs"
            ):
                AUDIT.validate_raw_reconstruction_dependency_receipt(
                    receipt, paths, lane, rng, runtime, source_pairs
                )


if __name__ == "__main__":
    unittest.main()
