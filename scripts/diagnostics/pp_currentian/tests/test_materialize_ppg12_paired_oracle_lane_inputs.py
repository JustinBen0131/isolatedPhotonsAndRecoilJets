#!/usr/bin/env python3
"""Synthetic tests for paired-oracle closure-lane materialization."""

from __future__ import annotations

import copy
import csv
import hashlib
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[4]
MODULE_PATH = (
    REPO
    / "scripts/diagnostics/pp_currentian/materialize_ppg12_paired_oracle_lane_inputs.py"
)
SPEC = importlib.util.spec_from_file_location("paired_lane_materializer", MODULE_PATH)
assert SPEC and SPEC.loader
MATERIALIZER = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = MATERIALIZER
SPEC.loader.exec_module(MATERIALIZER)

CLOSURE_CONTRACT = (
    REPO / "agent_context/analysis_contracts/ppg12_stitched_purity_closure.yaml"
)
CONTRACT = json.loads(CLOSURE_CONTRACT.read_text())


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


class FakeReader:
    def __init__(self, histograms: dict[str, dict[str, object]]):
        self.histograms = copy.deepcopy(histograms)
        self.closed = False

    def paths_for_basename(self, basename: str) -> list[str]:
        return sorted(
            path for path in self.histograms if path.rsplit("/", 1)[-1] == basename
        )

    def cycle_count(self, path: str) -> int:
        return 1 if path in self.histograms else 0

    def histogram(self, path: str) -> dict[str, object]:
        return copy.deepcopy(self.histograms[path])

    def close(self) -> None:
        self.closed = True


class PairedLaneMaterializerTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.run = self.root / "paired_run"
        self.run.mkdir()
        (self.run / "RUN_STATE").write_text("PASS\n")
        self.comparison = self.run / "comparison"
        self.comparison.mkdir()

        self.lane_id = "photon:photon5:0mrad:si"
        self.source_list = self.run / "photon5_g4_full.list"
        self.event_set = self.run / "photon5_g4_first5.list"
        source_rows = [f"/events/g4_segment_{index}.root" for index in range(5)]
        self.source_list.write_text("\n".join(source_rows) + "\n")
        self.event_set.write_text("\n".join(source_rows) + "\n")
        self.config = self.run / "analysis_config.yaml"
        self.config.write_text("canonical: lane-config\n")
        self.ppg_period_config = self.run / "ppg_recoeff_0rad_config.yaml"
        self.ppg_period_config.write_text(
            "canonical: preserved-ppg12-period-estimator\nperiod: 0mrad\n"
        )
        self.ppg_vertex_weight = self.run / "truth_vertex_reweight_0mrad.root"
        self.ppg_vertex_weight.write_bytes(b"sealed-0mrad-truth-vertex-weights")
        self.yaml_header_receipt = self.run / "yaml_cpp_header_tree_receipt.json"
        write_json(
            self.yaml_header_receipt,
            {"schema_version": 1, "role": "ppg_recoeff_yaml_cpp_header_tree"},
        )

        self.reference_root = self.run / "ppg12_efficiency.root"
        self.candidate_root = self.run / "recoiljets.root"
        self.reference_root.write_bytes(b"reference-root" + b"r" * 60_000)
        self.candidate_root.write_bytes(b"candidate-root" + b"c" * 60_000)
        self.trace = self.run / "ppg12_trace.csv"
        self.response_trace = self.run / "ppg12_response_trace.csv"
        self.trace.write_text("entry,cluster\n0,0\n")
        self.response_trace.write_text("entry,cluster\n0,0\n")

        self.runtime = self.run / "runtime"
        self.runtime.mkdir()
        self.receipt = self.runtime / "build_receipt.json"
        write_json(self.receipt, {"status": "PASS"})
        roles = sorted(
            {
                role
                for group in CONTRACT["lane_runtime_contract"]
                ["required_roles_by_lane_field"].values()
                for role in group
                if role != "lane_config"
            }
        )
        runtime_files = []
        for role in roles:
            if role == "ppg_recoeff_period_config":
                path = self.ppg_period_config
            elif role == "ppg_recoeff_truth_vertex_reweight":
                path = self.ppg_vertex_weight
            elif role == "ppg_recoeff_yaml_cpp_header_tree_receipt":
                path = self.yaml_header_receipt
            else:
                path = self.runtime / role.replace("/", "_")
                path.write_text(f"source-locked bytes for {role}\n")
            runtime_files.append(
                {"role": role, "path": str(path), "sha256": sha256(path)}
            )
        if "ppg_recoeff_period_config" not in roles:
            runtime_files.append(
                {
                    "role": "ppg_recoeff_period_config",
                    "path": str(self.ppg_period_config),
                    "sha256": sha256(self.ppg_period_config),
                }
            )
        self.runtime_manifest = self.runtime / "runtime_manifest.json"
        runtime_payload = {
            "schema_version": 1,
            **CONTRACT["lane_runtime_contract"]["required_manifest_values"],
            "build_receipt": str(self.receipt),
            "build_receipt_sha256": sha256(self.receipt),
            "selected_period": "0mrad",
            "files": runtime_files,
        }
        write_json(self.runtime_manifest, runtime_payload)
        self.paired_runtime_manifest = self.runtime / "paired_runtime_manifest.json"
        write_json(self.paired_runtime_manifest, runtime_payload)

        self.candidate_csv = self.comparison / "paired_oracle_candidates.csv"
        self.aggregate_json = self.comparison / "executable_aggregate.json"
        self.contract_path = self.run / "paired_oracle_contract.json"
        self.contract_payload = {
            "schema_version": 3,
            "lane": {
                "lane_id": self.lane_id,
                "sample": "Photon5",
                "period": "0mrad",
                "interaction": "SI",
                "rows": 5,
            },
            "paths": {
                "g4_full_list": str(self.source_list),
                "g4_slice": str(self.event_set),
                "recoil_runtime_manifest": str(self.runtime_manifest),
                "paired_runtime_manifest": str(self.paired_runtime_manifest),
                "recoil_config": str(self.config),
                "ppg_recoeff_period_config": str(self.ppg_period_config),
                "ppg_recoeff_truth_vertex_reweight": str(self.ppg_vertex_weight),
                "ppg_recoeff_yaml_cpp_header_tree_receipt": str(
                    self.yaml_header_receipt
                ),
                "recoil_root": str(self.candidate_root),
                "candidate_csv": str(self.candidate_csv),
                "ppg_recoeff_baseline_eff_root": str(self.reference_root),
                "ppg12_executable_trace": str(self.trace),
                "ppg12_executable_response_trace": str(self.response_trace),
                "ppg12_executable_aggregate": str(self.aggregate_json),
            },
        }
        write_json(self.contract_path, self.contract_payload)
        self.row = self._candidate_row()
        self._write_candidate_rows([self.row])
        self.aggregate_payload = self._aggregate_payload()
        write_json(self.aggregate_json, self.aggregate_payload)
        self.output = self.root / "materialized"

        self.edges = [float(item) for item in MATERIALIZER.AGGREGATE_HELPER.RECO_EDGES]
        self.histograms_by_root = {
            self.reference_root.resolve(): self._histograms(prefix=""),
            self.candidate_root.resolve(): self._histograms(prefix="SIM"),
        }

    def tearDown(self) -> None:
        self.temp.cleanup()

    def _candidate_row(self) -> dict[str, str]:
        contract_sha = sha256(self.contract_path)
        row = {
            "lane_id": self.lane_id,
            "runtime_contract_sha256": contract_sha,
            "match_status": "matched",
            "candidate_identity": "seg0:evt1:trk7:rj0:ppg0",
            "rj_cluster_Et": "11.0",
            "rj_weight_final": "2.0",
            "rj_signal_fill_multiplicity": "1",
            "ppg12_response_Et": "11.0",
            "ppg12_weight_final": "2.0",
            "ppg12_fill_multiplicity": "1",
            "ppg12_tag_evidence_source": "preserved_ppg12_executable",
            "ppg12_isolation_abcd_evidence_source": "preserved_ppg12_executable",
            "ppg12_truth_response_fill_evidence_source": "preserved_ppg12_executable",
            "ppg12_weight_evidence_source": "preserved_ppg12_executable",
        }
        for prefix in ("rj", "ppg12"):
            for region in MATERIALIZER.REGIONS:
                row[f"{prefix}_signal_fill_{region}"] = "1" if region == "A" else "0"
        return row

    def _write_candidate_rows(self, rows: list[dict[str, str]]) -> None:
        with self.candidate_csv.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)

    def _cells(self) -> dict[str, list[dict[str, float]]]:
        cells = {
            region: [
                {"content": 0.0, "sumw2": 0.0}
                for _ in range(len(MATERIALIZER.AGGREGATE_HELPER.RECO_EDGES) + 1)
            ]
            for region in MATERIALIZER.REGIONS
        }
        cells["A"][1] = {"content": 2.0, "sumw2": 4.0}
        return cells

    def _aggregate_payload(self) -> dict:
        links = {
            "runtime_contract": self.contract_path,
            "candidate_csv": self.candidate_csv,
            "runtime_manifest": self.paired_runtime_manifest,
            "baseline_root": self.reference_root,
            "trace_csv": self.trace,
            "response_trace_csv": self.response_trace,
        }
        return {
            "schema_version": 1,
            "evidence_source": "preserved_ppg12_executable_aggregate",
            "status": "PASS",
            "mode": "full",
            "certified_scopes": sorted(MATERIALIZER.REQUIRED_CERTIFIED_SCOPES),
            "uncertified_scopes": ["response_candidate_identity"],
            "root_equivalence": {
                "pass": True,
                "efficiency_root_exact": True,
                "response_root_exact": True,
            },
            "aggregate_comparison": {"pass": True, "cells": [], "leakage": []},
            "oracle_cells": self._cells(),
            "recoil_cells": self._cells(),
            "lane_identity": {
                "lane_id": self.lane_id,
                "runtime_contract_sha256": sha256(self.contract_path),
            },
            "provenance": {
                role: {"path": str(path), "sha256": sha256(path)}
                for role, path in links.items()
            },
        }

    def _histograms(self, *, prefix: str) -> dict[str, dict[str, object]]:
        histograms: dict[str, dict[str, object]] = {}
        for region in MATERIALIZER.REGIONS:
            basename = MATERIALIZER.LANE_EXTRACTOR.REGION_OBJECTS[f"{region}_signal"]
            path = f"{prefix}/{basename}" if prefix else basename
            values = [0.0] * (len(self.edges) - 1)
            sumw2 = [0.0] * (len(self.edges) - 1)
            entries = 0.0
            if region == "A":
                values[0] = 2.0
                sumw2[0] = 4.0
                entries = 1.0
            histograms[path] = {
                "bin_edges": self.edges,
                "sumw": values,
                "sumw2": sumw2,
                "entries": entries,
                "flow_sumw": [0.0, 0.0],
                "flow_sumw2": [0.0, 0.0],
                "classname": "TH1D",
            }
        return histograms

    def reader_factory(self, path: Path) -> FakeReader:
        return FakeReader(self.histograms_by_root[path.resolve()])

    def _refresh_candidate_provenance(self) -> None:
        self.aggregate_payload["provenance"]["candidate_csv"] = {
            "path": str(self.candidate_csv),
            "sha256": sha256(self.candidate_csv),
        }
        write_json(self.aggregate_json, self.aggregate_payload)

    def _rebind_contract_and_rows(self) -> None:
        contract_sha = sha256(self.contract_path)
        self.row["runtime_contract_sha256"] = contract_sha
        self._write_candidate_rows([self.row])
        self.aggregate_payload["lane_identity"]["runtime_contract_sha256"] = contract_sha
        self.aggregate_payload["provenance"]["runtime_contract"] = {
            "path": str(self.contract_path),
            "sha256": contract_sha,
        }
        self.aggregate_payload["provenance"]["candidate_csv"] = {
            "path": str(self.candidate_csv),
            "sha256": sha256(self.candidate_csv),
        }
        write_json(self.aggregate_json, self.aggregate_payload)

    def materialize(self) -> dict:
        return MATERIALIZER.materialize(
            self.contract_path,
            self.output,
            CLOSURE_CONTRACT,
            reader_factory=self.reader_factory,
        )

    def test_materializes_both_sides_and_canonical_extracts(self) -> None:
        result = self.materialize()
        self.assertEqual(result["status"], "PASS")
        for side, object_prefix in (("reference", ""), ("candidate", "SIM")):
            lane_input = json.loads(
                (self.output / side / "lane_input.json").read_text()
            )
            fills = json.loads(
                (self.output / side / "fill_evidence.json").read_text()
            )
            lane = json.loads((self.output / side / "lane.json").read_text())
            self.assertEqual(lane_input["schema"], "ppg12-stitched-purity-lane-input/v1")
            self.assertEqual(lane_input["object_prefix"], object_prefix)
            self.assertEqual(fills["observables"]["A_signal"]["fills"][0], 1)
            self.assertEqual(sum(fills["observables"]["B_signal"]["fills"]), 0)
            self.assertEqual(lane["observables"]["A_signal"]["fills"][0], 1.0)
            self.assertEqual(lane["observables"]["A_signal"]["sumw"][0], 2.0)
        reference_input = json.loads(
            (self.output / "reference" / "lane_input.json").read_text()
        )
        candidate_input = json.loads(
            (self.output / "candidate" / "lane_input.json").read_text()
        )
        expected_common_config = {
            "path": str(self.config.resolve()),
            "sha256": sha256(self.config),
        }
        self.assertEqual(reference_input["evidence"]["config"], expected_common_config)
        self.assertEqual(candidate_input["evidence"]["config"], expected_common_config)
        self.assertEqual(
            [row["role"] for row in reference_input["candidate_parity_evidence"]],
            ["candidate_rows"],
        )
        self.assertEqual(
            [row["role"] for row in candidate_input["candidate_parity_evidence"]],
            ["candidate_rows"],
        )
        reference_lane = json.loads(
            (self.output / "reference" / "lane.json").read_text()
        )
        candidate_lane = json.loads(
            (self.output / "candidate" / "lane.json").read_text()
        )
        self.assertEqual(
            reference_lane["config_sha256"], candidate_lane["config_sha256"]
        )
        manifest = json.loads(
            (self.output / "materialization_manifest.json").read_text()
        )
        self.assertEqual(set(manifest["sides"]), {"reference", "candidate"})
        self.assertEqual(
            manifest["ppg_recoeff_period_config"],
            {
                "path": str(self.ppg_period_config.resolve()),
                "sha256": sha256(self.ppg_period_config),
            },
        )
        self.assertEqual(
            manifest["ppg_recoeff_truth_vertex_reweight"],
            {
                "path": str(self.ppg_vertex_weight.resolve()),
                "sha256": sha256(self.ppg_vertex_weight),
            },
        )
        self.assertEqual(
            manifest["paired_runtime_manifest"],
            {
                "path": str(self.paired_runtime_manifest.resolve()),
                "sha256": sha256(self.paired_runtime_manifest),
            },
        )

    def test_trace_moment_mismatch_fails_without_output(self) -> None:
        self.row["rj_signal_fill_A"] = "0"
        self.row["rj_signal_fill_multiplicity"] = "0"
        self._write_candidate_rows([self.row])
        self._refresh_candidate_provenance()
        with self.assertRaisesRegex(MATERIALIZER.MaterializationError, "does not reproduce"):
            self.materialize()
        self.assertFalse(self.output.exists())

    def test_missing_fill_flag_cannot_fall_back_to_root_entries(self) -> None:
        self.row["ppg12_signal_fill_A"] = ""
        self._write_candidate_rows([self.row])
        self._refresh_candidate_provenance()
        with self.assertRaisesRegex(MATERIALIZER.MaterializationError, "not an integer"):
            self.materialize()
        self.assertFalse(self.output.exists())

    def test_reference_requires_preserved_executable_provenance(self) -> None:
        self.row["ppg12_weight_evidence_source"] = "python_shadow"
        self._write_candidate_rows([self.row])
        self._refresh_candidate_provenance()
        with self.assertRaisesRegex(
            MATERIALIZER.MaterializationError, "preserved-executable"
        ):
            self.materialize()
        self.assertFalse(self.output.exists())

    def test_stale_aggregate_candidate_link_fails(self) -> None:
        self.candidate_csv.write_text(self.candidate_csv.read_text() + "\n")
        with self.assertRaisesRegex(MATERIALIZER.MaterializationError, "hash mismatch"):
            self.materialize()
        self.assertFalse(self.output.exists())

    def test_canonical_extractor_failure_removes_partial_output(self) -> None:
        candidate_histograms = self.histograms_by_root[self.candidate_root.resolve()]
        candidate_histograms["SIM/h_tight_iso_cluster_signal_0"]["entries"] = 2.0
        with self.assertRaisesRegex(
            MATERIALIZER.MaterializationError, "canonical lane extraction rejected candidate"
        ):
            self.materialize()
        self.assertFalse(self.output.exists())

    def test_nonphoton_lane_is_explicitly_rejected(self) -> None:
        self.contract_payload["lane"].update(
            {
                "lane_id": "inclusive:jet8:0mrad:si",
                "sample": "Jet8",
            }
        )
        write_json(self.contract_path, self.contract_payload)
        with self.assertRaisesRegex(MATERIALIZER.MaterializationError, "photon-only"):
            self.materialize()
        self.assertFalse(self.output.exists())

    def test_missing_period_specific_ppg_config_is_rejected(self) -> None:
        del self.contract_payload["paths"]["ppg_recoeff_period_config"]
        write_json(self.contract_path, self.contract_payload)
        with self.assertRaisesRegex(
            MATERIALIZER.MaterializationError, "ppg_recoeff_period_config"
        ):
            self.materialize()
        self.assertFalse(self.output.exists())

    def test_missing_selected_runtime_manifest_is_rejected(self) -> None:
        del self.contract_payload["paths"]["paired_runtime_manifest"]
        write_json(self.contract_path, self.contract_payload)
        with self.assertRaisesRegex(
            MATERIALIZER.MaterializationError, "paired_runtime_manifest"
        ):
            self.materialize()
        self.assertFalse(self.output.exists())

    def test_selected_period_config_hash_drift_is_rejected(self) -> None:
        self.ppg_period_config.write_text("mutated after manifest sealing\n")
        with self.assertRaisesRegex(MATERIALIZER.MaterializationError, "hash mismatch"):
            self.materialize()
        self.assertFalse(self.output.exists())

    def test_selected_truth_vertex_weight_hash_drift_is_rejected(self) -> None:
        self.ppg_vertex_weight.write_bytes(b"mutated after manifest sealing")
        with self.assertRaisesRegex(MATERIALIZER.MaterializationError, "hash mismatch"):
            self.materialize()
        self.assertFalse(self.output.exists())

    def test_cross_period_truth_vertex_weight_swap_is_rejected(self) -> None:
        other = self.run / "truth_vertex_reweight_1p5mrad.root"
        other.write_bytes(b"sealed-1p5mrad-truth-vertex-weights")
        self.contract_payload["paths"]["ppg_recoeff_truth_vertex_reweight"] = str(
            other
        )
        write_json(self.contract_path, self.contract_payload)
        self._rebind_contract_and_rows()
        with self.assertRaisesRegex(
            MATERIALIZER.MaterializationError,
            "path differs from the paired runtime contract",
        ):
            self.materialize()
        self.assertFalse(self.output.exists())

    def test_cross_period_selected_config_swap_is_rejected(self) -> None:
        other = self.run / "ppg_recoeff_1p5mrad_config.yaml"
        other.write_text("canonical: preserved-ppg12-period-estimator\nperiod: 1p5mrad\n")
        self.contract_payload["paths"]["ppg_recoeff_period_config"] = str(other)
        write_json(self.contract_path, self.contract_payload)
        self._rebind_contract_and_rows()
        with self.assertRaisesRegex(
            MATERIALIZER.MaterializationError,
            "path differs from the paired runtime contract",
        ):
            self.materialize()
        self.assertFalse(self.output.exists())

    def test_selected_manifest_rejects_unselected_period_roles(self) -> None:
        payload = json.loads(self.paired_runtime_manifest.read_text())
        payload["files"].append(
            {
                "role": "ppg_recoeff_period_config_1p5mrad",
                "path": str(self.ppg_period_config),
                "sha256": sha256(self.ppg_period_config),
            }
        )
        write_json(self.paired_runtime_manifest, payload)
        self.aggregate_payload["provenance"]["runtime_manifest"] = {
            "path": str(self.paired_runtime_manifest),
            "sha256": sha256(self.paired_runtime_manifest),
        }
        write_json(self.aggregate_json, self.aggregate_payload)
        with self.assertRaisesRegex(
            MATERIALIZER.MaterializationError, "retains unselected period roles"
        ):
            self.materialize()
        self.assertFalse(self.output.exists())

    def test_refuses_to_overwrite_existing_output(self) -> None:
        self.output.mkdir()
        marker = self.output / "owned.txt"
        marker.write_text("do not replace\n")
        with self.assertRaisesRegex(MATERIALIZER.MaterializationError, "refusing to overwrite"):
            self.materialize()
        self.assertEqual(marker.read_text(), "do not replace\n")


if __name__ == "__main__":
    unittest.main()
