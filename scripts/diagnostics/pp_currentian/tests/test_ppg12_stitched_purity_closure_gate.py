#!/usr/bin/env python3
"""Synthetic regression tests for the PPG12 stitched-purity closure gate."""

from __future__ import annotations

import copy
import contextlib
import hashlib
import importlib.util
import io
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock


REPO = Path(__file__).resolve().parents[4]
MODULE_PATH = REPO / "scripts/diagnostics/pp_currentian/ppg12_stitched_purity_closure_gate.py"
CONTRACT_PATH = REPO / "agent_context/analysis_contracts/ppg12_stitched_purity_closure.yaml"
SPEC = importlib.util.spec_from_file_location("ppg12_stitched_purity_closure_gate", MODULE_PATH)
assert SPEC and SPEC.loader
GATE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(GATE)
CONTRACT = json.loads(CONTRACT_PATH.read_text())


class SyntheticCanonicalAssembler:
    """Isolate comparison-policy tests from canonical evidence construction."""

    def __init__(
        self,
        output_artifacts: list[dict[str, str]] | None = None,
        artifacts_by_merge: dict[str, list[dict[str, str]]] | None = None,
    ) -> None:
        self.output_artifacts = output_artifacts or []
        self.artifacts_by_merge = artifacts_by_merge or {}

    def verify_manifest_evidence(self, *_args, **_kwargs):
        return None

    def verify_merge_audit_evidence(self, *args, **_kwargs):
        merge_path = str(Path(args[0]).resolve()) if args else ""
        return {}, self.artifacts_by_merge.get(merge_path, self.output_artifacts)


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def histogram(value: float) -> dict[str, object]:
    return {
        "bin_edges": [10.0, 12.0, 14.0],
        "sumw": [value, value * 0.8],
        "sumw2": [value * 0.2, value * 0.16],
        "fills": [value * 10.0, value * 8.0],
    }


def candidate_parity() -> dict[str, object]:
    lane_coverage = []
    for sample in CONTRACT["families"]["photon"]["samples"]:
        for period in CONTRACT["periods"]:
            for interaction in CONTRACT["interactions"]:
                lane_coverage.append(
                    {
                        "lane_id": f"photon:{sample}:{period}:{interaction}",
                        "population": {
                            "reference_count": 10,
                            "candidate_count": 10,
                            "unmatched_reference": 0,
                            "unmatched_candidate": 0,
                        },
                        "features": {
                            "comparisons": 110,
                            "violations": 0,
                            "max_normalized_delta": 0.0,
                        },
                        "scores": {
                            "comparisons": 20,
                            "mismatches": 0,
                            "max_abs_delta": 0.0,
                        },
                        "model_route": {"comparisons": 10, "mismatches": 0},
                        "tags": {
                            name: {"comparisons": 10, "mismatches": 0}
                            for name in CONTRACT["required_tag_checks"]
                        },
                    }
                )
    return {
        "tolerances": dict(CONTRACT["candidate_parity_tolerances"]),
        "lane_coverage": lane_coverage,
        "population": {
            "reference_count": 120,
            "candidate_count": 120,
            "unmatched_reference": 0,
            "unmatched_candidate": 0,
        },
        "features": [
            {
                "name": name,
                "comparisons": 120,
                "violations": 0,
                "max_abs_delta": 0.0,
                "max_normalized_delta": 0.0,
            }
            for name in CONTRACT["required_features"]
        ],
        "scores": [
            {"name": name, "comparisons": 120, "mismatches": 0, "max_abs_delta": 0.0}
            for name in CONTRACT["required_score_models"]
        ],
        "model_route": {"comparisons": 120, "mismatches": 0},
        "tags": {
            name: {"comparisons": 120, "mismatches": 0}
            for name in CONTRACT["required_tag_checks"]
        },
    }


def manifest(role: str) -> dict[str, object]:
    provenance = {
        field: ("reference-implementation" if field == "implementation_sha256" and role == "reference" else
                "candidate-implementation" if field == "implementation_sha256" else f"shared-{field}")
        for field in CONTRACT["frozen_provenance_fields"]
    }
    lanes: list[dict[str, object]] = []
    for family, family_spec in CONTRACT["families"].items():
        for sample in family_spec["samples"]:
            for period in CONTRACT["periods"]:
                for interaction in CONTRACT["interactions"]:
                    lane_id = f"{family}:{sample}:{period}:{interaction}"
                    event_rows = [f"{lane_id}:source-{index}" for index in range(5)]
                    group_sha = GATE._payload_sha256(
                        {
                            "lane_id": lane_id,
                            "group_index": 0,
                            "event_rows": event_rows,
                        }
                    )
                    groups = [
                        {
                            "group_index": 0,
                            "group_id": f"group-00000-{group_sha[:16]}",
                            "event_rows": event_rows,
                            "group_sha256": group_sha,
                        }
                    ]
                    observables = {
                        name: histogram(10.0 + len(lanes) + index)
                        for index, name in enumerate(family_spec["required_observables"])
                    }
                    lane: dict[str, object] = {
                        "lane_id": lane_id,
                        "family": family,
                        "sample": sample,
                        "period": period,
                        "interaction": interaction,
                        "group_count": 1,
                        "group_size": 5,
                        "event_set_row_count": 5,
                        "groups": groups,
                        "external_scale": 1.0,
                        "observables": observables,
                    }
                    lane.update({field: f"{lane_id}-{field}" for field in CONTRACT["paired_lane_fields"]})
                    lane["group_set_sha256"] = GATE._payload_sha256(groups)
                    lanes.append(lane)
    purity_core = {
        "bin_edges": [10.0, 12.0, 14.0],
        "truth": {"value": [0.7, 0.72], "error": [0.01, 0.01]},
        "raw": {"value": [0.6, 0.62], "error": [0.02, 0.02]},
        "corrected": {"value": [0.68, 0.70], "error": [0.025, 0.025]},
    }
    purity_sha = GATE._payload_sha256(purity_core)
    purity = {
        **purity_core,
        "fixed_seed_repetition": {
            "first_output_sha256": purity_sha,
            "repeated_output_sha256": purity_sha,
        },
    }
    payload: dict[str, object] = {
        "schema": "ppg12-stitched-purity-manifest/v1",
        "role": role,
        "abcd_population": "unsuffixed",
        "external_scale": 1.0,
        "random_seed": 42,
        "toy_count": 20000,
        "provenance": provenance,
        "lanes": lanes,
        "purity": purity,
    }
    if role != "reference":
        payload["candidate_parity"] = candidate_parity()
    return payload


def merge_audit(candidate_path: Path) -> dict[str, object]:
    inclusive = [f"inclusive-input-{index}" for index in range(20)]
    photon = [f"photon-input-{index}" for index in range(12)]
    return {
        "schema": "ppg12-stitched-purity-merge-audit/v1",
        "candidate_manifest": {
            "path": str(candidate_path.resolve()),
            "sha256": file_sha256(candidate_path),
        },
        "audits": [
            {
                "family": "inclusive",
                "status": "PASS",
                "failures": [],
                "max_content_delta": 0.0,
                "max_sumw2_delta": 0.0,
                "inputs_fixed_order": inclusive,
                "inputs_fixed_order_count": len(inclusive),
            },
            {
                "family": "photon",
                "status": "PASS",
                "failures": [],
                "max_content_delta": 0.0,
                "max_sumw2_delta": 0.0,
                "inputs_fixed_order": photon,
                "inputs_fixed_order_count": len(photon),
            },
        ]
    }


class ClosureGateTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.reference = manifest("reference")
        self.candidate = manifest("candidate")
        self.reference_path = self.root / "reference.json"
        self.candidate_path = self.root / "candidate.json"
        self.merge_path = self.root / "merge.json"
        self.outdir = self.root / "out"
        self._flush()

    def tearDown(self) -> None:
        self.temp.cleanup()

    def _flush(self) -> None:
        write_json(self.reference_path, self.reference)
        write_json(self.candidate_path, self.candidate)
        write_json(self.merge_path, merge_audit(self.candidate_path))

    def _run(
        self, command: str = "diagnose", *, enforce_canonical: bool = False
    ) -> tuple[int, dict[str, object]]:
        argv = [
            "--contract", str(CONTRACT_PATH), command,
            "--reference-manifest", str(self.reference_path),
            "--candidate-manifest", str(self.candidate_path),
            "--outdir", str(self.outdir),
        ]
        if command == "admit":
            argv.extend(["--merge-audit", str(self.merge_path)])
        canonical_context = (
            contextlib.nullcontext()
            if enforce_canonical
            else mock.patch.object(
                GATE,
                "_canonical_assembler",
                return_value=SyntheticCanonicalAssembler(),
            )
        )
        with canonical_context, contextlib.redirect_stdout(io.StringIO()):
            code = GATE.main(argv)
        return code, json.loads((self.outdir / "gate_report.json").read_text())

    def _failure_codes(self, report: dict[str, object]) -> set[str]:
        return {row["code"] for row in report["failures"]}

    def _lane(self, payload: dict[str, object], lane_id: str) -> dict[str, object]:
        return next(row for row in payload["lanes"] if row["lane_id"] == lane_id)

    def test_exact_pair_passes_and_is_deterministic(self) -> None:
        code, report = self._run()
        first = (self.outdir / "gate_report.json").read_bytes()
        self.assertEqual(code, 0)
        self.assertEqual(report["status"], "PASS")
        code, _ = self._run()
        self.assertEqual(code, 0)
        self.assertEqual(first, (self.outdir / "gate_report.json").read_bytes())

    def test_pass_shaped_manifests_are_rejected_without_canonical_evidence(self) -> None:
        code, report = self._run(enforce_canonical=True)
        self.assertNotEqual(code, 0)
        self.assertIn("manifest_invalid", self._failure_codes(report))

    def test_admit_emits_manifest_only_on_pass(self) -> None:
        code, report = self._run("admit")
        admission = self.outdir / "admission_manifest.json"
        self.assertEqual(code, 0)
        self.assertEqual(report["status"], "PASS")
        payload = json.loads(admission.read_text())
        self.assertEqual(payload["schema"], "ppg12-stitched-purity-admission/v1")
        self.assertIn("gate_report_payload_sha256", payload)
        self.candidate["external_scale"] = 2.3398
        self._flush()
        code, report = self._run("admit")
        self.assertNotEqual(code, 0)
        self.assertEqual(report["status"], "FAIL")
        self.assertFalse(admission.exists())

    def test_unsuffixed_population_is_required(self) -> None:
        self.candidate["abcd_population"] = "classed"
        self._flush()
        _, report = self._run()
        self.assertIn("candidate_contract_failed", self._failure_codes(report))

    def test_lane_level_jet8_scale_is_forbidden(self) -> None:
        lane = self._lane(self.candidate, "inclusive:jet8:0mrad:si")
        lane["external_scale"] = 2.3398
        self._flush()
        _, report = self._run()
        self.assertIn("candidate_contract_failed", self._failure_codes(report))

    def test_every_lane_requires_complete_five_file_groups(self) -> None:
        for payload in (self.reference, self.candidate):
            lane = self._lane(payload, "inclusive:jet8:0mrad:si")
            lane["group_size"] = 4
        self._flush()
        _, report = self._run()
        self.assertIn("canary_underpopulated", self._failure_codes(report))

    def test_additional_complete_groups_are_allowed(self) -> None:
        lane_id = "inclusive:jet8:0mrad:si"
        for payload in (self.reference, self.candidate):
            lane = self._lane(payload, lane_id)
            event_rows = [f"{lane_id}:source-{index}" for index in range(10)]
            groups = []
            for group_index in range(2):
                rows = event_rows[group_index * 5 : (group_index + 1) * 5]
                group_sha = GATE._payload_sha256(
                    {
                        "lane_id": lane_id,
                        "group_index": group_index,
                        "event_rows": rows,
                    }
                )
                groups.append(
                    {
                        "group_index": group_index,
                        "group_id": f"group-{group_index:05d}-{group_sha[:16]}",
                        "event_rows": rows,
                        "group_sha256": group_sha,
                    }
                )
            lane.update(
                {
                    "group_count": 2,
                    "event_set_row_count": 10,
                    "groups": groups,
                    "group_set_sha256": GATE._payload_sha256(groups),
                }
            )
        self._flush()
        code, report = self._run()
        self.assertEqual(code, 0)
        self.assertEqual(report["status"], "PASS")

    def test_abcd_content_mutation_fails(self) -> None:
        lane = self._lane(self.candidate, "inclusive:jet8:0mrad:si")
        lane["observables"]["A"]["sumw"][0] *= 1.01
        self._flush()
        _, report = self._run()
        self.assertIn("inclusive_abcd_failed", self._failure_codes(report))

    def test_sumw2_only_mutation_fails(self) -> None:
        lane = self._lane(self.candidate, "inclusive:jet8:0mrad:si")
        lane["observables"]["A"]["sumw2"][0] *= 1.01
        self._flush()
        _, report = self._run()
        self.assertIn("inclusive_abcd_failed", self._failure_codes(report))

    def test_photon_leakage_mutation_fails(self) -> None:
        lane = self._lane(self.candidate, "photon:photon5:0mrad:si")
        lane["observables"]["C_signal"]["sumw"][0] *= 1.01
        self._flush()
        _, report = self._run()
        self.assertIn("photon_leakage_failed", self._failure_codes(report))

    def test_truth_component_mutation_fails(self) -> None:
        lane = self._lane(self.candidate, "inclusive:jet12:1p5mrad:di")
        lane["observables"]["A_notmatch"]["sumw"][1] *= 1.01
        self._flush()
        _, report = self._run()
        self.assertIn("inclusive_abcd_failed", self._failure_codes(report))

    def test_feature_and_model_route_mutations_fail(self) -> None:
        self.candidate["candidate_parity"]["features"][0]["violations"] = 1
        self.candidate["candidate_parity"]["model_route"]["mismatches"] = 1
        self._flush()
        _, report = self._run()
        self.assertIn("candidate_contract_failed", self._failure_codes(report))

    def test_candidate_parity_requires_exact_lane_coverage_and_tolerances(self) -> None:
        self.candidate["candidate_parity"]["lane_coverage"].pop()
        self.candidate["candidate_parity"]["tolerances"]["score_abs"] = 1e-5
        self._flush()
        _, report = self._run()
        self.assertIn("candidate_contract_failed", self._failure_codes(report))

    def test_candidate_parity_global_totals_must_equal_lane_sum(self) -> None:
        self.candidate["candidate_parity"]["features"][0]["comparisons"] += 1
        self._flush()
        _, report = self._run()
        self.assertIn("candidate_contract_failed", self._failure_codes(report))

    def test_fixed_seed_repetition_must_bind_exact_purity_payload(self) -> None:
        self.candidate["purity"]["fixed_seed_repetition"]["repeated_output_sha256"] = "0" * 64
        self._flush()
        _, report = self._run()
        self.assertIn("purity_solver_failed", self._failure_codes(report))

    def test_missing_duplicate_and_si_di_lane_mismatch_fail(self) -> None:
        self.candidate["lanes"].pop()
        duplicate = copy.deepcopy(self.candidate["lanes"][0])
        self.candidate["lanes"].append(duplicate)
        self._flush()
        _, report = self._run()
        self.assertIn("source_contract_failed", self._failure_codes(report))

    def test_binning_mismatch_fails(self) -> None:
        lane = self._lane(self.candidate, "photon:photon10:1p5mrad:si")
        lane["observables"]["A_signal"]["bin_edges"] = [10.0, 11.5, 14.0]
        self._flush()
        _, report = self._run()
        self.assertIn("binning_mismatch", self._failure_codes(report))

    def test_missing_root_family_fails(self) -> None:
        lane = self._lane(self.candidate, "photon:photon20:0mrad:di")
        del lane["observables"]["D_signal"]
        self._flush()
        _, report = self._run()
        self.assertIn("root_family_missing", self._failure_codes(report))

    def test_underpopulated_stitched_bin_is_retryable_failure(self) -> None:
        for payload in (self.reference, self.candidate):
            for lane in payload["lanes"]:
                if lane["family"] == "inclusive":
                    for field in ("sumw", "sumw2", "fills"):
                        lane["observables"]["A"][field][0] = 0.0
        self._flush()
        _, report = self._run()
        self.assertIn("canary_underpopulated", self._failure_codes(report))

    def test_zero_photon_leakage_final_bin_is_underpopulated(self) -> None:
        for payload in (self.reference, self.candidate):
            for lane in payload["lanes"]:
                if lane["family"] == "photon":
                    for field in ("sumw", "sumw2", "fills"):
                        lane["observables"]["D_signal"][field][0] = 0.0
        self._flush()
        _, report = self._run()
        self.assertIn("canary_underpopulated", self._failure_codes(report))

    def test_purity_mutation_fails(self) -> None:
        self.candidate["purity"]["corrected"]["value"][0] *= 1.001
        self._flush()
        _, report = self._run()
        self.assertIn("stitched_purity_failed", self._failure_codes(report))

    def test_purity_uncertainty_mutation_fails(self) -> None:
        self.candidate["purity"]["corrected"]["error"][0] *= 1.01
        purity_core = {
            key: self.candidate["purity"][key]
            for key in ("bin_edges", "truth", "raw", "corrected")
        }
        purity_sha = GATE._payload_sha256(purity_core)
        self.candidate["purity"]["fixed_seed_repetition"] = {
            "first_output_sha256": purity_sha,
            "repeated_output_sha256": purity_sha,
        }
        self._flush()
        _, report = self._run()
        self.assertIn("stitched_purity_failed", self._failure_codes(report))

    def test_duplicate_or_omitted_merge_input_blocks_admission(self) -> None:
        payload = merge_audit(self.candidate_path)
        payload["audits"][0]["inputs_fixed_order"][-1] = payload["audits"][0]["inputs_fixed_order"][0]
        write_json(self.merge_path, payload)
        _, report = self._run("admit")
        self.assertIn("merge_arithmetic_failed", self._failure_codes(report))
        self.assertFalse((self.outdir / "admission_manifest.json").exists())

    def test_stale_merge_audit_candidate_binding_blocks_admission(self) -> None:
        stale_candidate = self.root / "stale_candidate.json"
        write_json(stale_candidate, self.candidate)
        payload = merge_audit(stale_candidate)
        write_json(self.merge_path, payload)
        _, report = self._run("admit")
        self.assertIn("merge_arithmetic_failed", self._failure_codes(report))
        self.assertFalse((self.outdir / "admission_manifest.json").exists())

    def test_present_the94_summary_is_intentionally_not_admissible(self) -> None:
        legacy = (
            REPO
            / "InputFiles/the94_ppg12_inclusive_archivedreco_final_20260719_1615EDT"
            / "final_merged_roots/final_root_arithmetic_audit_summary.json"
        )
        if not legacy.exists():
            self.skipTest("local THE-94 final audit is not present")
        self.merge_path = legacy
        _, report = self._run("admit")
        self.assertIn("merge_arithmetic_failed", self._failure_codes(report))
        self.assertFalse((self.outdir / "admission_manifest.json").exists())

    def test_verify_production_passes_and_stale_admission_blocks(self) -> None:
        code, _ = self._run("admit")
        self.assertEqual(code, 0)
        admission_path = self.outdir / "admission_manifest.json"

        production_candidate = copy.deepcopy(self.candidate)
        production_candidate["role"] = "production"
        production_reference = copy.deepcopy(self.reference)
        # Broad production contains more events than its admission canary.
        # This paired event-set change is allowed; source-list and physics
        # contract hashes remain frozen.
        for ref_lane, cand_lane in zip(production_reference["lanes"], production_candidate["lanes"]):
            ref_lane["event_set_sha256"] = "full-production-event-set"
            cand_lane["event_set_sha256"] = "full-production-event-set"
        production_reference_path = self.root / "production_reference.json"
        production_candidate_path = self.root / "production_candidate.json"
        candidate_purity_path = self.root / "production_candidate_purity.json"
        candidate_purity_core = {
            key: copy.deepcopy(production_candidate["purity"][key])
            for key in ("bin_edges", "truth", "raw", "corrected")
        }
        candidate_purity_sha = GATE._payload_sha256(candidate_purity_core)
        write_json(
            candidate_purity_path,
            {
                "schema": "ppg12-stitched-purity-purity/v1",
                "random_seed": 42,
                "toy_count": 20000,
                "purity": candidate_purity_core,
                "fixed_seed_repetition": {
                    "first_output_sha256": candidate_purity_sha,
                    "repeated_output_sha256": candidate_purity_sha,
                },
            },
        )
        production_candidate["assembly"] = {
            "purity": {
                "path": str(candidate_purity_path),
                "sha256": file_sha256(candidate_purity_path),
            }
        }
        write_json(production_reference_path, production_reference)
        write_json(production_candidate_path, production_candidate)
        production_merge_path = self.root / "production_merge.json"
        write_json(production_merge_path, merge_audit(production_candidate_path))
        inclusive_root = self.root / "inclusive.root"
        photon_root = self.root / "photon.root"
        inclusive_root.write_bytes(b"synthetic inclusive ROOT placeholder")
        photon_root.write_bytes(b"synthetic photon ROOT placeholder")
        source_paths = {
            "historical_purity": self.root / "historical_purity.root",
            "historical_abcd": self.root / "historical_abcd.root",
            "historical_leakage": self.root / "historical_purity.root",
            "candidate_purity": candidate_purity_path,
            "candidate_inclusive": inclusive_root,
            "candidate_photon": photon_root,
        }
        for role in ("historical_purity", "historical_abcd"):
            source_paths[role].write_bytes(role.encode())
        source_links = [
            {"role": role, "path": str(path.resolve()), "sha256": file_sha256(path)}
            for role, path in sorted(source_paths.items())
        ]
        historical = {
            "schema": "ppg12-stitched-purity-historical-comparison/v1",
            "final_purity_series": "corrected",
            "chi2_ndf": 1.0,
            "max_abs_pull": 2.0,
            "weighted_mean_ratio": 1.0,
            "weighted_mean_ratio_error": 0.01,
            "abcd_coherent_trend_max_sigma": 1.0,
            "leakage_coherent_trend_max_sigma": 1.0,
            "final_purity": {"point_count": 2},
            "abcd": {name: {} for name in ("A", "B", "C", "D")},
            "leakage": {name: {} for name in ("cB", "cC", "cD")},
            "input": {"mode": "direct_root_objects"},
            "source_links": source_links,
            "source_set_sha256": GATE._payload_sha256(source_links),
        }
        historical_path = self.root / "historical.json"
        write_json(historical_path, historical)
        candidate_artifacts = [
            {"family": "inclusive", "path": str(inclusive_root.resolve()), "sha256": file_sha256(inclusive_root)},
            {"family": "photon", "path": str(photon_root.resolve()), "sha256": file_sha256(photon_root)},
        ]
        candidate_artifacts.sort(key=lambda row: row["family"])
        historical_metric_fields = (
            "chi2_ndf",
            "max_abs_pull",
            "weighted_mean_ratio",
            "weighted_mean_ratio_error",
            "abcd_coherent_trend_max_sigma",
            "leakage_coherent_trend_max_sigma",
        )
        historical_metrics = {
            key: historical[key] for key in historical_metric_fields
        }
        wrapper = {
            "schema": "ppg12-stitched-purity-production/v1",
            "admission_sha256": file_sha256(admission_path),
            "reference_manifest": {"path": str(production_reference_path), "sha256": file_sha256(production_reference_path)},
            "candidate_manifest": {"path": str(production_candidate_path), "sha256": file_sha256(production_candidate_path)},
            "merge_audit": {
                "path": str(production_merge_path.resolve()),
                "sha256": file_sha256(production_merge_path),
            },
            "candidate_artifacts": candidate_artifacts,
            "historical_archive_comparison": historical_metrics,
            "assembly": {
                "contract_sha256": GATE._payload_sha256(CONTRACT),
                "historical_comparison": {
                    "path": str(historical_path),
                    "sha256": file_sha256(historical_path),
                },
                "candidate_artifact_set_sha256": GATE._payload_sha256(candidate_artifacts),
            },
        }
        wrapper_path = self.root / "production_wrapper.json"
        write_json(wrapper_path, wrapper)
        production_out = self.root / "production_out"
        def run_production_gate() -> int:
            with mock.patch.object(
                GATE,
                "_canonical_assembler",
                return_value=SyntheticCanonicalAssembler(
                    artifacts_by_merge={
                        str(self.merge_path.resolve()): [],
                        str(production_merge_path.resolve()): candidate_artifacts,
                    }
                ),
            ), contextlib.redirect_stdout(io.StringIO()):
                return GATE.main([
                    "--contract", str(CONTRACT_PATH), "verify-production",
                    "--admission-manifest", str(admission_path),
                    "--production-manifest", str(wrapper_path),
                    "--merge-audit", str(production_merge_path),
                    "--outdir", str(production_out),
                ])

        code = run_production_gate()
        self.assertEqual(
            code,
            0,
            (production_out / "gate_report.json").read_text()
            if (production_out / "gate_report.json").exists()
            else "missing production gate report",
        )
        production_report = json.loads((production_out / "production_gate_report.json").read_text())
        self.assertEqual(production_report["status"], "PASS")
        self.assertEqual({row["family"] for row in production_report["candidate_artifacts"]}, {"inclusive", "photon"})

        # A self-consistent PASS-shaped admission remains invalid when its
        # embedded gate receipt was not produced by replaying the pair gate.
        pristine_admission = admission_path.read_bytes()
        forged_admission = json.loads(pristine_admission)
        forged_admission["gate_report_payload_sha256"] = "a" * 64
        write_json(admission_path, forged_admission)
        wrapper["admission_sha256"] = file_sha256(admission_path)
        write_json(wrapper_path, wrapper)
        code = run_production_gate()
        failed = json.loads((production_out / "gate_report.json").read_text())
        self.assertNotEqual(code, 0)
        self.assertIn("hash_drift", self._failure_codes(failed))
        admission_path.write_bytes(pristine_admission)
        wrapper["admission_sha256"] = file_sha256(admission_path)
        write_json(wrapper_path, wrapper)

        photon_root.write_bytes(b"mutated photon candidate")
        code = run_production_gate()
        failed = json.loads((production_out / "gate_report.json").read_text())
        self.assertNotEqual(code, 0)
        self.assertIn("hash_drift", self._failure_codes(failed))
        photon_root.write_bytes(b"synthetic photon ROOT placeholder")

        wrapper["historical_archive_comparison"]["leakage_coherent_trend_max_sigma"] = 2.1
        write_json(wrapper_path, wrapper)
        code = run_production_gate()
        failed = json.loads((production_out / "gate_report.json").read_text())
        self.assertNotEqual(code, 0)
        self.assertIn("promotion_blocked", self._failure_codes(failed))
        self.assertIn("hash_drift", self._failure_codes(failed))
        wrapper["historical_archive_comparison"] = historical_metrics
        write_json(wrapper_path, wrapper)

        production_candidate["provenance"]["implementation_sha256"] = "changed-after-admission"
        write_json(production_candidate_path, production_candidate)
        write_json(production_merge_path, merge_audit(production_candidate_path))
        wrapper["candidate_manifest"]["sha256"] = file_sha256(production_candidate_path)
        write_json(wrapper_path, wrapper)
        code = run_production_gate()
        failed = json.loads((production_out / "gate_report.json").read_text())
        self.assertNotEqual(code, 0)
        self.assertIn("hash_drift", self._failure_codes(failed))
        self.assertFalse((production_out / "production_gate_report.json").exists())


if __name__ == "__main__":
    unittest.main()
