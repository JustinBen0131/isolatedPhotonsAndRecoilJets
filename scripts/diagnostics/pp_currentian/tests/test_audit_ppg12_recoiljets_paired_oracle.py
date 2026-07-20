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


def source_import_fixture(root: Path) -> tuple[dict, dict[str, Path]]:
    current_library = root / "runtime" / "lib" / "libCaloAna24.so"
    provenance = root / "runtime" / "provenance"
    provenance.mkdir(parents=True)
    current_library.parent.mkdir(parents=True, exist_ok=True)
    current_library.write_bytes(b"exact sealed Attempt-12 binary\n")
    origin_receipt_path = provenance / "ppg_source_build_receipt.json"
    origin_manifest_path = provenance / "ppg_source_runtime_manifest.json"
    origin_receipt = {
        "schema_version": 1,
        "ppg12_source": {
            "revision": AUDIT.EXPECTED_PPG_SOURCE_REVISION,
            "working_tree_ignored": True,
            "rebuilt_against_common_runtime": True,
        },
        "staged_rewrites": {
            "archived_ppg12_binary_reused": False,
            "ppg12_source_locked_rebuild": True,
        },
    }
    origin_receipt_path.write_text(json.dumps(origin_receipt, sort_keys=True) + "\n")
    origin_manifest = {
        "schema_version": 1,
        "runtime_profile": AUDIT.EXPECTED_RUNTIME_PROFILE,
        "offline_main": AUDIT.EXPECTED_OFFLINE_MAIN,
        "isolated_build": True,
        "estimator_revision": AUDIT.EXPECTED_ESTIMATOR_REVISION,
        # Immutable origin documents preserve their old absolute paths.  The
        # copied roles and hashes, not those historical paths, are authoritative.
        "build_receipt": "/retired/origin/build_receipt.json",
        "build_receipt_sha256": AUDIT.sha256(origin_receipt_path),
        "files": [{
            "role": "libCaloAna24.so",
            "path": "/retired/origin/runtime/lib/libCaloAna24.so",
            "sha256": AUDIT.sha256(current_library),
        }],
    }
    origin_manifest_path.write_text(json.dumps(origin_manifest, sort_keys=True) + "\n")
    by_role = {
        "libCaloAna24.so": current_library,
        "ppg_source_runtime_manifest": origin_manifest_path,
        "ppg_source_build_receipt": origin_receipt_path,
    }
    receipt = {
        "ppg12_source": {
            "revision": AUDIT.EXPECTED_PPG_SOURCE_REVISION,
            "working_tree_ignored": True,
            "rebuilt_against_common_runtime": True,
            "binary_mode": AUDIT.PPG_BINARY_IMPORT,
            "rebuilt_in_this_runtime": False,
            "source_runtime_import": {
                "runtime_manifest": {
                    "path": str(origin_manifest_path),
                    "sha256": AUDIT.sha256(origin_manifest_path),
                },
                "build_receipt": {
                    "path": str(origin_receipt_path),
                    "sha256": AUDIT.sha256(origin_receipt_path),
                },
                "library_sha256": AUDIT.sha256(current_library),
                "immutable_provenance_documents": True,
            },
        },
        "staged_rewrites": {
            "archived_ppg12_binary_reused": False,
            "ppg12_source_locked_rebuild": False,
            "ppg12_source_locked_binary_import": True,
        },
    }
    return receipt, by_role


def roounfold_runtime_fixture(
    root: Path, *, relocate_pcm: bool = False
) -> tuple[dict, dict[str, Path], str, str, str, Path]:
    libdir = root / "runtime" / "lib"
    include_root = root / "runtime" / "estimator" / "include"
    libdir.mkdir(parents=True)
    include_root.mkdir(parents=True)
    library = libdir / "libRooUnfold.so"
    pcm_root = root / "external" if relocate_pcm else libdir
    pcm_root.mkdir(parents=True, exist_ok=True)
    pcm = pcm_root / "RooUnfoldDict_rdict.pcm"
    library.write_bytes(b"synthetic historical RooUnfold library\n")
    pcm.write_bytes(b"matching synthetic RooUnfold dictionary\n")
    rows = []
    tree_digest = hashlib.sha256()
    for name in AUDIT.EXPECTED_ROOUNFOLD_HEADERS:
        header = include_root / name
        header.write_text(f"// sealed historical {name}\n")
        digest = AUDIT.sha256(header)
        rows.append({"relative_path": name, "sha256": digest})
        tree_digest.update(name.encode())
        tree_digest.update(b"\0")
        tree_digest.update(bytes.fromhex(digest))
    receipt = root / "runtime" / "estimator" / "roounfold_header_tree_receipt.json"
    receipt.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "role": "ppg_recoeff_roounfold_header_tree",
                "include_root": str(include_root),
                "tree_sha256": tree_digest.hexdigest(),
                "files": rows,
            }
        )
        + "\n"
    )
    assets = [library, pcm, receipt] + [
        include_root / name for name in AUDIT.EXPECTED_ROOUNFOLD_HEADERS
    ]
    estimator = {
        "runtime_assets": [
            {"path": str(path), "sha256": AUDIT.sha256(path)} for path in assets
        ]
    }
    by_role = {
        "ppg_recoeff_roounfold": library,
        "ppg_recoeff_roounfold_pcm": pcm,
        "ppg_recoeff_roounfold_header_tree_receipt": receipt,
        "ppg_recoeff_roounfold_response_header": include_root / "RooUnfoldResponse.h",
        "ppg_recoeff_roounfold_bayes_header": include_root / "RooUnfoldBayes.h",
    }
    return (
        estimator,
        by_role,
        AUDIT.sha256(library),
        AUDIT.sha256(pcm),
        tree_digest.hexdigest(),
        include_root,
    )


class TestFirstDivergenceAudit(unittest.TestCase):
    def test_source_locked_binary_import_provenance_passes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            receipt, by_role = source_import_fixture(Path(tmp))
            self.assertEqual(
                AUDIT.validate_ppg_binary_mode(receipt), AUDIT.PPG_BINARY_IMPORT
            )
            AUDIT.validate_source_locked_ppg_import(
                receipt, by_role, AUDIT.EXPECTED_OFFLINE_MAIN
            )

    def test_source_locked_binary_import_rejects_library_drift(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            receipt, by_role = source_import_fixture(Path(tmp))
            by_role["libCaloAna24.so"].write_bytes(b"different binary\n")
            with self.assertRaisesRegex(
                AUDIT.AuditFailure, "differs from the sealed origin binary"
            ):
                AUDIT.validate_source_locked_ppg_import(
                    receipt, by_role, AUDIT.EXPECTED_OFFLINE_MAIN
                )

    def test_source_locked_binary_import_rejects_ambiguous_mode_flags(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            receipt, _ = source_import_fixture(Path(tmp))
            receipt["staged_rewrites"]["ppg12_source_locked_rebuild"] = True
            with self.assertRaisesRegex(AUDIT.AuditFailure, "ambiguous"):
                AUDIT.validate_ppg_binary_mode(receipt)

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

    def test_roounfold_runtime_requires_hashed_colocated_pcm(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            estimator, by_role, lib_hash, pcm_hash, tree_hash, _ = \
                roounfold_runtime_fixture(Path(tmp))
            AUDIT.validate_roounfold_runtime(
                estimator, by_role, lib_hash, pcm_hash, tree_hash
            )

            by_role["ppg_recoeff_roounfold_pcm"].write_bytes(b"drifted dictionary\n")
            with self.assertRaisesRegex(AUDIT.AuditFailure, "PCM.*hash mismatch"):
                AUDIT.validate_roounfold_runtime(
                    estimator, by_role, lib_hash, pcm_hash, tree_hash
                )

    def test_roounfold_runtime_rejects_non_colocated_pcm(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            estimator, by_role, lib_hash, pcm_hash, tree_hash, _ = roounfold_runtime_fixture(
                Path(tmp), relocate_pcm=True
            )
            with self.assertRaisesRegex(AUDIT.AuditFailure, "not co-located"):
                AUDIT.validate_roounfold_runtime(
                    estimator, by_role, lib_hash, pcm_hash, tree_hash
                )

    def test_roounfold_runtime_rejects_header_drift(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            estimator, by_role, lib_hash, pcm_hash, tree_hash, include_root = \
                roounfold_runtime_fixture(Path(tmp))
            (include_root / "RooUnfold.h").write_text("// modern fallback header\n")
            with self.assertRaisesRegex(AUDIT.AuditFailure, "header drifted"):
                AUDIT.validate_roounfold_runtime(
                    estimator, by_role, lib_hash, pcm_hash, tree_hash
                )

    def test_roounfold_runtime_rejects_recomputed_substitute_pcm(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            estimator, by_role, lib_hash, pcm_hash, tree_hash, _ = \
                roounfold_runtime_fixture(Path(tmp))
            pcm = by_role["ppg_recoeff_roounfold_pcm"]
            pcm.write_bytes(b"internally consistent but substituted PCM\n")
            for row in estimator["runtime_assets"]:
                if Path(row["path"]).resolve() == pcm.resolve():
                    row["sha256"] = AUDIT.sha256(pcm)
            with self.assertRaisesRegex(AUDIT.AuditFailure, "PCM.*hash mismatch"):
                AUDIT.validate_roounfold_runtime(
                    estimator, by_role, lib_hash, pcm_hash, tree_hash
                )

    def test_roounfold_runtime_rejects_recomputed_substitute_header_tree(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            estimator, by_role, lib_hash, pcm_hash, tree_hash, include_root = \
                roounfold_runtime_fixture(Path(tmp))
            receipt = by_role["ppg_recoeff_roounfold_header_tree_receipt"]
            (include_root / "RooUnfold.h").write_text(
                "// internally consistent but substituted RooUnfold.h\n"
            )
            data = json.loads(receipt.read_text())
            digest = hashlib.sha256()
            for row in data["files"]:
                header = include_root / row["relative_path"]
                row["sha256"] = AUDIT.sha256(header)
                digest.update(row["relative_path"].encode())
                digest.update(b"\0")
                digest.update(bytes.fromhex(row["sha256"]))
            data["tree_sha256"] = digest.hexdigest()
            receipt.write_text(json.dumps(data) + "\n")
            for row in estimator["runtime_assets"]:
                path = Path(row["path"])
                if path.resolve() in {
                    receipt.resolve(),
                    (include_root / "RooUnfold.h").resolve(),
                }:
                    row["sha256"] = AUDIT.sha256(path)
            with self.assertRaisesRegex(AUDIT.AuditFailure, "pinned digest"):
                AUDIT.validate_roounfold_runtime(
                    estimator, by_role, lib_hash, pcm_hash, tree_hash
                )

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
