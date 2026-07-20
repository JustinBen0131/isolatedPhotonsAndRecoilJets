#!/usr/bin/env python3
"""Pure aggregation/contract tests for the 12-lane photon audit runner."""

from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).resolve().parents[1] / "run_ppg12_photon_oracle_canary_audit.py"
SPEC = importlib.util.spec_from_file_location("photon_oracle_audit", MODULE_PATH)
assert SPEC and SPEC.loader
AUDIT = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = AUDIT
SPEC.loader.exec_module(AUDIT)

ASSEMBLER_PATH = MODULE_PATH.with_name("assemble_ppg12_stitched_purity_manifest.py")
ASSEMBLER_SPEC = importlib.util.spec_from_file_location(
    "stitched_purity_manifest_assembler", ASSEMBLER_PATH
)
assert ASSEMBLER_SPEC and ASSEMBLER_SPEC.loader
ASSEMBLER = importlib.util.module_from_spec(ASSEMBLER_SPEC)
sys.modules[ASSEMBLER_SPEC.name] = ASSEMBLER
ASSEMBLER_SPEC.loader.exec_module(ASSEMBLER)

CONTRACT_PATH = (
    Path(__file__).resolve().parents[4]
    / "agent_context/analysis_contracts/ppg12_stitched_purity_closure.yaml"
)
CONTRACT = json.loads(CONTRACT_PATH.read_text())


def matched_row(lane_id: str | None = None) -> dict[str, str]:
    row = {
        "match_status": "matched",
        "candidate_identity": "seg0:evt1:trk7:rj2:ppg2",
        "identity_match": "1",
        "segment": "0",
        "eventnumber": "1",
        "truth_track_id": "7",
        "rj_cluster_index": "2",
        "ppg12_cluster_index": "2",
        "rj_selected_model": "base_v3E",
        "rj_inferred_stored_model": "base_v3E",
        "ppg12_selected_model": "base_v3E",
        "model_route_agree": "1",
        "ppg12_tag_evidence_source": "preserved_ppg12_executable",
        "ppg12_isolation_abcd_evidence_source": "preserved_ppg12_executable",
        "ppg12_truth_response_fill_evidence_source": "preserved_ppg12_executable",
        "ppg12_weight_evidence_source": "preserved_ppg12_executable",
        "common_agree": "1",
        "stored_common_agree": "1",
        "tag_agree": "1",
        "stored_tag_agree": "1",
        "signal_status_agree": "1",
        "rj_recomputed_common_pass": "1",
        "rj_stored_common_pass": "1",
        "ppg12_common_pass": "1",
        "rj_recomputed_tag": "1",
        "rj_stored_tag": "1",
        "ppg12_recomputed_tag": "1",
        "rj_is_iso": "1",
        "ppg12_is_iso": "1",
        "rj_is_noniso": "0",
        "ppg12_is_noniso": "0",
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
        "rj_weight_component_code": (
            "2" if lane_id is not None and lane_id.endswith(":di") else "1"
        ),
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
        "rj_abcd_region": "A",
        "rj_stored_abcd_region": "A",
        "ppg12_abcd_region": "A",
        "rj_event_weight": "1.0",
        "rj_xsec_pb": "2.0",
        "rj_xsec_weight": "0.5",
        "rj_window_low": "0.0",
        "rj_window_high": "14.0",
    }
    for name in AUDIT.FEATURE_NAMES:
        row[f"{name}_rj"] = "1.0"
        row[f"{name}_ppg12"] = "1.0"
        row[f"{name}_delta"] = "0"
    for _, left, right in AUDIT.SCORE_PAIRS:
        row[left] = "0.8"
        row[right] = "0.8"
    row.update(
        {
            "base_E_score_delta": "0",
            "base_v3E_score_delta": "0",
            "selected_score_delta": "0",
            "stored_minus_routed_score": "0",
            "rj_stored_bdt_score": "0.8",
            "ppg12_selected_bdt_score": "0.8",
        }
    )
    for _, left, right in AUDIT.ISOLATION_PAIRS:
        row[left] = "0.1"
        row[right] = "0.1"
    return row


def file_item(path: Path) -> dict[str, str]:
    return {"path": str(path.resolve()), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def make_bundle(base: Path, lane_id: str, index: int) -> tuple[dict[str, Path], dict[str, str]]:
    lane_dir = base / f"lane_{index:02d}"
    lane_dir.mkdir()
    _, sample_key, period, interaction = lane_id.split(":")
    sample = "Photon" + sample_key.removeprefix("photon")
    trace = lane_dir / "trace.csv"
    response = lane_dir / "response.csv"
    candidate = lane_dir / "candidate.csv"
    aggregate = lane_dir / "aggregate.json"
    contract = lane_dir / "contract.json"
    trace.write_text(f"trace,{lane_id}\n")
    response.write_text(f"response,{lane_id}\n")
    contract.write_text(json.dumps({
        "schema_version": 3,
        "lane": {"lane_id": lane_id, "sample": sample, "period": period,
                 "interaction": interaction.upper(), "rows": 5},
        "paths": {
            "candidate_csv": str(candidate.resolve()),
            "ppg12_executable_aggregate": str(aggregate.resolve()),
            "ppg12_executable_trace": str(trace.resolve()),
            "ppg12_executable_response_trace": str(response.resolve()),
        },
    }, sort_keys=True) + "\n")
    contract_hash = hashlib.sha256(contract.read_bytes()).hexdigest()
    row = matched_row(lane_id)
    row.update(lane_id=lane_id, runtime_contract_sha256=contract_hash)
    with candidate.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(row))
        writer.writeheader(); writer.writerow(row)
    aggregate.write_text(json.dumps({
        "schema_version": 1,
        "evidence_source": "preserved_ppg12_executable_aggregate",
        "status": "PASS", "mode": "full",
        "lane_identity": {"lane_id": lane_id, "runtime_contract_sha256": contract_hash},
        "provenance": {
            "runtime_contract": file_item(contract), "candidate_csv": file_item(candidate),
            "trace_csv": file_item(trace), "response_trace_csv": file_item(response),
        },
    }, sort_keys=True) + "\n")
    bundle = {
        "runtime_contract": contract, "candidate_csv": candidate,
        "executable_aggregate": aggregate, "trace_csv": trace,
        "response_trace_csv": response,
    }
    return bundle, row


class TestCandidateAggregation(unittest.TestCase):
    def test_exact_row_passes_same_event_scope(self) -> None:
        summary = AUDIT.summarize_candidate_rows([matched_row()])
        self.assertIsNone(summary["first_divergence"])
        self.assertTrue(all(value == 0 for value in summary["stage_failure_counts"].values()))

    def test_population_precedes_feature_and_score(self) -> None:
        matched = matched_row()
        matched["cluster_et1_rj"] = "2.0"
        matched["rj_stored_bdt_score"] = "0.2"
        unmatched = dict(matched)
        unmatched["match_status"] = "rj_only"
        summary = AUDIT.summarize_candidate_rows([matched, unmatched])
        self.assertEqual(summary["first_divergence"], "candidate_population")

    def test_feature_precedes_score(self) -> None:
        row = matched_row()
        row["e11_over_e33_rj"] = "1.1"
        row["rj_stored_bdt_score"] = "0.2"
        summary = AUDIT.summarize_candidate_rows([row])
        self.assertEqual(summary["first_divergence"], "features")

    def test_stored_route_tag_and_abcd_are_not_hidden_by_recomputation(self) -> None:
        row = matched_row()
        row["rj_stored_bdt_score"] = "0.2"
        row["rj_inferred_stored_model"] = "base_E"
        row["rj_stored_tag"] = "2"
        row["stored_tag_agree"] = "0"
        row["rj_stored_abcd_region"] = "C"
        row["stored_abcd_agree"] = "0"
        summary = AUDIT.summarize_candidate_rows([row])
        self.assertEqual(summary["first_divergence"], "scores")
        self.assertEqual(summary["stage_failure_counts"]["model_route"], 1)
        self.assertEqual(summary["stage_failure_counts"]["tags"], 1)
        self.assertEqual(summary["stage_failure_counts"]["isolation_abcd"], 1)


class TestCandidateParityArtifact(unittest.TestCase):
    def lane_rows(self) -> dict[str, list[dict[str, str]]]:
        return {
            lane_id: [matched_row(lane_id)]
            for lane_id in AUDIT.expected_lane_ids()
        }

    def test_exact_12_lane_artifact_matches_gate_schema_and_aliases(self) -> None:
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        bundles = {}
        lane_rows = {}
        for index, lane_id in enumerate(sorted(AUDIT.expected_lane_ids())):
            bundle, row = make_bundle(Path(temp.name), lane_id, index)
            bundles[lane_id] = bundle
            lane_rows[lane_id] = [row]
        payload = AUDIT.build_candidate_parity(
            lane_rows,
            created_utc="2026-07-20T00:00:00+00:00",
            trace_bundle_by_lane=bundles,
        )
        self.assertEqual(payload["schema"], AUDIT.CANDIDATE_PARITY_SCHEMA)
        self.assertEqual(
            payload["name_aliases"]["features"],
            {"cluster_Et_score_input": "cluster_Et"},
        )
        self.assertEqual(
            payload["name_aliases"]["scores"], {"base_v3E": "baseV3E"}
        )
        parity = payload["candidate_parity"]
        self.assertEqual(parity["tolerances"], CONTRACT["candidate_parity_tolerances"])
        self.assertEqual(len(parity["lane_coverage"]), 12)
        self.assertEqual(
            [row["lane_id"] for row in parity["lane_coverage"]],
            sorted(AUDIT.expected_lane_ids()),
        )
        self.assertEqual(
            [row["name"] for row in parity["features"]], CONTRACT["required_features"]
        )
        self.assertEqual(
            [row["name"] for row in parity["scores"]],
            CONTRACT["required_score_models"],
        )

        path = Path(temp.name) / "candidate_parity.json"
        path.write_text(json.dumps(payload))
        normalized, _ = ASSEMBLER._validate_candidate_parity(path, CONTRACT)
        self.assertEqual(len(normalized["lane_coverage"]), 12)

    def test_candidate_only_markers_are_inadmissible(self) -> None:
        with self.assertRaisesRegex(RuntimeError, "paths alone are inadmissible"):
            AUDIT.build_candidate_parity(
                self.lane_rows(),
                candidate_row_paths={
                    lane: Path("/tmp/hand-authored.csv") for lane in AUDIT.expected_lane_ids()
                },
            )

    def test_duplicate_trace_assigned_to_two_lanes_fails(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            bundles, rows = {}, {}
            for index, lane_id in enumerate(sorted(AUDIT.expected_lane_ids())):
                bundle, row = make_bundle(Path(tmp), lane_id, index)
                bundles[lane_id], rows[lane_id] = bundle, [row]
            lanes = sorted(bundles)
            bundles[lanes[1]] = dict(bundles[lanes[1]])
            bundles[lanes[1]]["candidate_csv"] = bundles[lanes[0]]["candidate_csv"]
            with self.assertRaises(RuntimeError):
                AUDIT.build_candidate_parity(rows, trace_bundle_by_lane=bundles)

    def test_marker_csv_without_trace_provenance_fails(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            bundles, rows = {}, {}
            for index, lane_id in enumerate(sorted(AUDIT.expected_lane_ids())):
                bundle, row = make_bundle(Path(tmp), lane_id, index)
                bundles[lane_id], rows[lane_id] = bundle, [row]
            lane = sorted(bundles)[0]
            aggregate = json.loads(bundles[lane]["executable_aggregate"].read_text())
            aggregate["provenance"].pop("trace_csv")
            bundles[lane]["executable_aggregate"].write_text(json.dumps(aggregate))
            with self.assertRaisesRegex(RuntimeError, "trace_csv provenance is missing"):
                AUDIT.build_candidate_parity(rows, trace_bundle_by_lane=bundles)

    def test_global_totals_are_exact_lane_sums(self) -> None:
        parity = AUDIT.build_candidate_parity(self.lane_rows())["candidate_parity"]
        lanes = parity["lane_coverage"]
        for name in (
            "reference_count",
            "candidate_count",
            "unmatched_reference",
            "unmatched_candidate",
        ):
            self.assertEqual(
                parity["population"][name],
                sum(row["population"][name] for row in lanes),
            )
        self.assertEqual(
            sum(row["comparisons"] for row in parity["features"]),
            sum(row["features"]["comparisons"] for row in lanes),
        )
        self.assertEqual(
            sum(row["comparisons"] for row in parity["scores"]),
            sum(row["scores"]["comparisons"] for row in lanes),
        )
        for name in ("comparisons", "mismatches"):
            self.assertEqual(
                parity["model_route"][name],
                sum(row["model_route"][name] for row in lanes),
            )
        for tag in CONTRACT["required_tag_checks"]:
            for name in ("comparisons", "mismatches"):
                self.assertEqual(
                    parity["tags"][tag][name],
                    sum(row["tags"][tag][name] for row in lanes),
                )

    def test_missing_lane_fails_closed(self) -> None:
        rows = self.lane_rows()
        rows.pop(AUDIT.expected_lane_ids()[-1])
        with self.assertRaisesRegex(RuntimeError, "not the exact 12 photon lanes"):
            AUDIT.build_candidate_parity(rows)

    def test_di_component_code_one_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            bundles, rows = {}, {}
            for index, lane_id in enumerate(sorted(AUDIT.expected_lane_ids())):
                bundle, row = make_bundle(Path(tmp), lane_id, index)
                bundles[lane_id], rows[lane_id] = bundle, [row]

            lane_id = "photon:photon5:0mrad:di"
            bad_row = dict(rows[lane_id][0])
            bad_row["rj_weight_component_code"] = "1"
            candidate_path = bundles[lane_id]["candidate_csv"]
            with candidate_path.open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=list(bad_row))
                writer.writeheader()
                writer.writerow(bad_row)
            aggregate_path = bundles[lane_id]["executable_aggregate"]
            aggregate = json.loads(aggregate_path.read_text())
            aggregate["provenance"]["candidate_csv"] = file_item(candidate_path)
            aggregate_path.write_text(json.dumps(aggregate, sort_keys=True) + "\n")
            rows[lane_id] = [bad_row]

            with self.assertRaisesRegex(
                RuntimeError, "candidate trace is not an executable parity PASS"
            ):
                AUDIT.build_candidate_parity(
                    rows,
                    created_utc="2026-07-20T00:00:00+00:00",
                    trace_bundle_by_lane=bundles,
                )


class TestPlanContract(unittest.TestCase):
    def test_exact_plan_and_receipts_are_accepted(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            source_manifest = base / "source.json"; source_manifest.write_text("{}\n")
            driver = base / "driver.sh"; driver.write_text("#!/bin/sh\n")
            common = {}
            for role in (
                "setup_script", "apply_bdt", "apply_config", "base_e_model",
                "base_v3e_model", "npb_model", "tower_mask",
                "recoil_runtime_manifest", "recoil_config",
            ):
                path = base / role; path.write_text(role + "\n"); common[role] = file_item(path)
            lanes, receipts = [], []
            for index, lane_id in enumerate(AUDIT.expected_lane_ids()):
                bundle, _ = make_bundle(base, lane_id, index)
                lane_dir = bundle["runtime_contract"].parent
                (lane_dir / "RUN_STATE").write_text("PASS\n")
                _, sample_key, period, interaction = lane_id.split(":")
                sources = {}
                for role in ("ppg_macro", "g4_full_list", "truthjet_full_list"):
                    path = lane_dir / role; path.write_text(f"{lane_id} {role}\n")
                    sources[role] = file_item(path)
                lanes.append({
                    "lane_id": lane_id,
                    "sample": "Photon" + sample_key.removeprefix("photon"),
                    "period": period, "interaction": interaction.upper(), "rows": 5,
                    "execution": "source_locked_paired_executable",
                    "output_base": str(lane_dir), "sources": sources,
                    "expected_evidence": {
                        "runtime_contract": str(bundle["runtime_contract"]),
                        "candidate_csv": str(bundle["candidate_csv"]),
                        "executable_aggregate": str(bundle["executable_aggregate"]),
                        "run_state": str(lane_dir / "RUN_STATE"),
                    },
                })
                receipts.append({
                    "lane_id": lane_id, "output_base": str(lane_dir), "status": "PASS",
                    "runtime_contract": str(bundle["runtime_contract"]),
                    "runtime_contract_sha256": AUDIT.sha256_file(bundle["runtime_contract"]),
                    "candidate_csv": str(bundle["candidate_csv"]),
                    "candidate_csv_sha256": AUDIT.sha256_file(bundle["candidate_csv"]),
                    "executable_aggregate": str(bundle["executable_aggregate"]),
                    "executable_aggregate_sha256": AUDIT.sha256_file(bundle["executable_aggregate"]),
                })
            auth = {
                "schema": AUDIT.PLAN_SCHEMA, "campaign_tag": "unit",
                "output_root": str(base), "source_manifest": file_item(source_manifest),
                "paired_driver": file_item(driver), "common": common, "lanes": lanes,
            }
            doc = dict(auth)
            doc["submission_token"] = "RUN_PPG12_PAIRED_" + AUDIT.canonical_payload_sha256(auth)
            plan = base / "plan.json"; receipt = base / "receipts.tsv"
            plan.write_text(json.dumps(doc))
            with receipt.open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=receipts[0], delimiter="\t")
                writer.writeheader()
                writer.writerows(receipts)
            parsed = AUDIT.validate_plan_and_receipts(plan, receipt)
            self.assertEqual(len(parsed), 12)


if __name__ == "__main__":
    unittest.main()
