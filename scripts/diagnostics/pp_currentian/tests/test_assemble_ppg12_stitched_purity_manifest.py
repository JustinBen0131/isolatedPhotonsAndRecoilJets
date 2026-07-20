#!/usr/bin/env python3
"""Regression tests for deterministic stitched-purity manifest assembly."""

from __future__ import annotations

import contextlib
import copy
import csv
import hashlib
import importlib.util
import io
import json
import subprocess
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[4]
ASSEMBLER_PATH = (
    REPO
    / "scripts/diagnostics/pp_currentian/assemble_ppg12_stitched_purity_manifest.py"
)
GATE_PATH = (
    REPO / "scripts/diagnostics/pp_currentian/ppg12_stitched_purity_closure_gate.py"
)
RUNNER_PATH = (
    REPO / "scripts/diagnostics/pp_currentian/run_ppg12_photon_oracle_canary_audit.py"
)
EXTRACTOR_PATH = (
    REPO / "scripts/diagnostics/pp_currentian/extract_ppg12_stitched_purity_lane.py"
)
CONTRACT_PATH = (
    REPO / "agent_context/analysis_contracts/ppg12_stitched_purity_closure.yaml"
)


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


ASSEMBLER = load_module("assemble_ppg12_stitched_purity_manifest", ASSEMBLER_PATH)
GATE = load_module("ppg12_stitched_purity_closure_gate_for_assembler", GATE_PATH)
RUNNER = load_module("ppg12_photon_oracle_runner_for_assembler", RUNNER_PATH)
EXTRACTOR = load_module("ppg12_lane_extractor_for_assembler", EXTRACTOR_PATH)
CONTRACT = json.loads(CONTRACT_PATH.read_text())


def digest(label: str) -> str:
    return hashlib.sha256(label.encode()).hexdigest()


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def passing_candidate_row(lane_id: str | None = None) -> dict[str, str]:
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
    }
    for name in RUNNER.FEATURE_NAMES:
        row[f"{name}_rj"] = "1"
        row[f"{name}_ppg12"] = "1"
        row[f"{name}_delta"] = "0"
    for _, left, right in RUNNER.SCORE_PAIRS:
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
    return row


class ManifestAssemblerTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name).resolve()
        self.lane_dir = self.root / "lanes"
        self.artifact_dir = self.root / "lane_roots"
        self.lane_paths: list[Path] = []
        self.lane_payloads: dict[str, dict[str, object]] = {}
        self.original_lane_reproducer = ASSEMBLER._reproduce_lane_extract
        self.original_purity_reproducer = ASSEMBLER._reproduce_purity_evidence
        self.original_admission_replayer = ASSEMBLER._execute_canonical_admission_replay
        self.original_gate_canonical_assembler = GATE._canonical_assembler
        self._make_lanes()
        self._install_lane_replay_fixture()
        self._install_purity_replay_fixture()
        GATE._canonical_assembler = lambda _contract: ASSEMBLER
        ASSEMBLER._execute_canonical_admission_replay = self._replay_admission_fixture
        self.provenance_path = self.root / "provenance.json"
        self.purity_path = self.root / "purity.json"
        self.parity_path = self.root / "candidate_parity.json"
        self.index_path = self.root / "lane_index.json"
        self.coverage_path = self.root / "global_minimality.json"
        self.reference_path = self.root / "reference.json"
        self.candidate_path = self.root / "candidate.json"
        self.production_candidate_path = self.root / "production_candidate.json"
        self.merge_path = self.root / "merge_audit.json"
        self._make_supporting_payloads()
        self.assertEqual(self._make_index(), 0, self.last_output)
        self.assertEqual(self._make_coverage(), 0, self.last_output)

    def tearDown(self) -> None:
        ASSEMBLER._reproduce_lane_extract = self.original_lane_reproducer
        ASSEMBLER._reproduce_purity_evidence = self.original_purity_reproducer
        ASSEMBLER._execute_canonical_admission_replay = self.original_admission_replayer
        GATE._canonical_assembler = self.original_gate_canonical_assembler
        self.temp.cleanup()

    @staticmethod
    def _replay_admission_fixture(
        _gate_path: Path,
        contract_path: Path,
        reference_path: str,
        candidate_path: str,
        merge_path: str,
        outdir: Path,
    ) -> tuple[int, str]:
        with contextlib.redirect_stdout(io.StringIO()):
            code = GATE.main(
                [
                    "--contract",
                    str(contract_path),
                    "admit",
                    "--reference-manifest",
                    reference_path,
                    "--candidate-manifest",
                    candidate_path,
                    "--merge-audit",
                    merge_path,
                    "--outdir",
                    str(outdir),
                ]
            )
        return code, "fixture canonical admission replay"

    def _install_lane_replay_fixture(self) -> None:
        replayed = {
            payload["extraction"]["input_sidecar"]["path"]: copy.deepcopy(payload)
            for payload in self.lane_payloads.values()
        }

        def reproduce(_extractor: Path, link: dict[str, str], _contract: dict[str, object]):
            return copy.deepcopy(replayed[link["path"]])

        ASSEMBLER._reproduce_lane_extract = reproduce

    def _install_purity_replay_fixture(self) -> None:
        scientific = {
            "bin_edges": [10.0, 12.0, 14.0],
            "truth": {"value": [0.5, 0.5], "error": [0.01, 0.01]},
            "raw": {"value": [0.8, 0.8], "error": [0.02, 0.02]},
            "corrected": {"value": [0.7, 0.7], "error": [0.02, 0.02]},
        }
        purity_sha = ASSEMBLER._payload_sha256(scientific)
        diagnostics: list[dict[str, object]] = []
        diagnostics_sha = ASSEMBLER._payload_sha256(diagnostics)
        algorithm = {
            "name": "PPG12 CalculatePhotonYield truth/raw/leakage-corrected estimator",
            "random_stream": "one TRandom3(42) stream per full repeated evaluation",
            "toys_per_bin": 20000,
            "abcd_toys": "effective-Poisson throws using sumw^2/sumw2",
            "leakage_toys": "Gaussian cB/cC/cD throws from photon signal TH1::Divide",
            "root_solver": "executable PPG12 myfunc evaluated with ROOT TF1::GetX over [-0.5*A,2*A]",
            "toy_histogram": "1000 bins over [-1,2]",
            "fit_window": "mean-RMS to mean+1.5*RMS",
            "fit_model": "ROOT Gaussian with REMQN options",
            "truth": "inclusive A_signal divided by unsuffixed inclusive A via TGraphAsymmErrors",
            "source": {
                "path": str((REPO / CONTRACT["canonical_tools"]["ppg12_estimator_source"]).resolve()),
                "sha256": file_sha256(
                    REPO / CONTRACT["canonical_tools"]["ppg12_estimator_source"]
                ),
            },
            "producer": {
                "path": str((REPO / CONTRACT["canonical_tools"]["purity_producer"]).resolve()),
                "sha256": file_sha256(
                    REPO / CONTRACT["canonical_tools"]["purity_producer"]
                ),
            },
        }

        def reproduce(
            contract_path: Path,
            contract: dict[str, object],
            _lanes: list[dict[str, object]],
            lane_links: list[dict[str, str]],
        ) -> dict[str, object]:
            normalized_links = sorted(lane_links, key=lambda row: row["lane_id"])
            return {
                "schema": "ppg12-stitched-purity-purity/v1",
                "random_seed": 42,
                "toy_count": 20000,
                "purity": copy.deepcopy(scientific),
                "fixed_seed_repetition": {
                    "first_output_sha256": purity_sha,
                    "repeated_output_sha256": purity_sha,
                },
                "run_diagnostics": {
                    "first": diagnostics,
                    "repeated": diagnostics,
                    "diagnostics_sha256": diagnostics_sha,
                    "repeated_diagnostics_sha256": diagnostics_sha,
                },
                "algorithm": copy.deepcopy(algorithm),
                "inputs": {
                    "contract": {
                        "path": str(contract_path.resolve()),
                        "sha256": ASSEMBLER._payload_sha256(contract),
                    },
                    "lane_index": None,
                    "lanes": normalized_links,
                    "lane_set_sha256": ASSEMBLER._payload_sha256(normalized_links),
                },
            }

        ASSEMBLER._reproduce_purity_evidence = reproduce

    def _run(self, argv: list[str]) -> int:
        captured = io.StringIO()
        with contextlib.redirect_stdout(captured):
            code = ASSEMBLER.main(["--contract", str(CONTRACT_PATH), *argv])
        self.last_output = captured.getvalue()
        return code

    def _make_lanes(self) -> None:
        specs: list[tuple[str, str, str, str, str, Path, int]] = []
        index = 0
        for family, family_spec in CONTRACT["families"].items():
            for sample in family_spec["samples"]:
                for period in CONTRACT["periods"]:
                    for interaction in CONTRACT["interactions"]:
                        lane_id = f"{family}:{sample}:{period}:{interaction}"
                        root_path = self.artifact_dir / f"{lane_id.replace(':', '_')}.root"
                        specs.append(
                            (family, sample, period, interaction, lane_id, root_path, index)
                        )
                        index += 1
        self.assertEqual(index, 32)
        self.artifact_dir.mkdir(parents=True, exist_ok=True)
        root_specs = self.root / "root_specs.json"
        write_json(
            root_specs,
            [{"path": str(row[5]), "weight": row[6] + 1} for row in specs],
        )
        subprocess.run(
            [
                str(REPO.parent / "analysis/env/bin/python3"),
                "-c",
                (
                    "import hashlib,json,ROOT,sys; from pathlib import Path; "
                    "specs=json.loads(Path(sys.argv[1]).read_text()); "
                    "[(lambda f,w: ("
                    "(lambda h: (h.Sumw2(),h.Fill(0.5,w),h.Fill(1.5,2*w),h.Write()))"
                    "(ROOT.TH1D('merge_hist','',2,0,2)),"
                    "ROOT.TNamed('analysis_config_yaml','unit').Write(),"
                    "ROOT.TObjString(''.join(hashlib.sha256(str(i).encode()).hexdigest() "
                    "for i in range(3000))).Write('terminal_padding'),f.Close()))"
                    "(ROOT.TFile(s['path'],'RECREATE'),float(s['weight'])) for s in specs]"
                ),
                str(root_specs),
            ],
            check=True,
            capture_output=True,
            text=True,
        )

        runtime_dir = self.root / "runtime"
        runtime_dir.mkdir(exist_ok=True)
        build_receipt = runtime_dir / "build_receipt.json"
        write_json(build_receipt, {"status": "source-locked-test-runtime"})
        runtime_roles = sorted(
            {
                role
                for roles in CONTRACT["lane_runtime_contract"]
                ["required_roles_by_lane_field"].values()
                for role in roles
                if role != "lane_config"
            }
        )
        runtime_files = []
        for role in runtime_roles:
            runtime_path = runtime_dir / role.replace("/", "_")
            runtime_path.write_text(f"executed bytes for {role}\n")
            runtime_files.append(
                {
                    "role": role,
                    "path": str(runtime_path),
                    "sha256": file_sha256(runtime_path),
                }
            )
        runtime_manifest = runtime_dir / "runtime_manifest.json"
        write_json(
            runtime_manifest,
            {
                "schema_version": 1,
                **CONTRACT["lane_runtime_contract"]["required_manifest_values"],
                "build_receipt": str(build_receipt),
                "build_receipt_sha256": file_sha256(build_receipt),
                "files": runtime_files,
            },
        )
        runtime_manifest_link = {
            "path": str(runtime_manifest),
            "sha256": file_sha256(runtime_manifest),
        }

        self.candidate_paths: dict[str, Path] = {}
        self.candidate_rows: dict[str, list[dict[str, str]]] = {}
        self.trace_bundles: dict[str, dict[str, Path]] = {}
        for family, sample, period, interaction, lane_id, root_path, index in specs:
            lane_evidence_dir = self.root / "lane_evidence" / lane_id.replace(":", "_")
            lane_evidence_dir.mkdir(parents=True, exist_ok=True)
            evidence_links: dict[str, dict[str, str]] = {}
            for evidence_name in sorted(ASSEMBLER.LANE_EVIDENCE_KEYS):
                evidence_path = lane_evidence_dir / f"{evidence_name}.txt"
                if evidence_name in {"source_list", "event_set"}:
                    evidence_path.write_text(
                        "\n".join(f"{lane_id}:source-{row}" for row in range(5)) + "\n"
                    )
                else:
                    evidence_path.write_text(f"{lane_id}:{evidence_name}\n")
                evidence_links[evidence_name] = {
                    "path": str(evidence_path),
                    "sha256": file_sha256(evidence_path),
                }

            candidate_links: list[dict[str, str]] = []
            if family == "photon":
                candidate_path = lane_evidence_dir / "candidate_rows.csv"
                aggregate_path = lane_evidence_dir / "executable_aggregate.json"
                trace_path = lane_evidence_dir / "executable_trace.csv"
                response_trace_path = lane_evidence_dir / "response_trace.csv"
                paired_contract_path = lane_evidence_dir / "paired_runtime_contract.json"
                trace_path.write_text(f"lane_id,trace_value\n{lane_id},{index}\n")
                response_trace_path.write_text(
                    f"lane_id,response_value\n{lane_id},{index}\n"
                )
                write_json(
                    paired_contract_path,
                    {
                        "lane": {"lane_id": lane_id},
                        "paths": {
                            "candidate_csv": str(candidate_path),
                            "ppg12_executable_aggregate": str(aggregate_path),
                            "ppg12_executable_trace": str(trace_path),
                            "ppg12_executable_response_trace": str(response_trace_path),
                        },
                    },
                )
                paired_contract_sha = file_sha256(paired_contract_path)
                row = passing_candidate_row(lane_id)
                row["lane_id"] = lane_id
                row["runtime_contract_sha256"] = paired_contract_sha
                with candidate_path.open("w", newline="") as stream:
                    writer = csv.DictWriter(stream, fieldnames=list(row))
                    writer.writeheader()
                    writer.writerow(row)
                write_json(
                    aggregate_path,
                    {
                        "evidence_source": "preserved_ppg12_executable_aggregate",
                        "status": "PASS",
                        "mode": "full",
                        "lane_identity": {
                            "lane_id": lane_id,
                            "runtime_contract_sha256": paired_contract_sha,
                        },
                        "provenance": {
                            role: {
                                "path": str(evidence_path.resolve()),
                                "sha256": file_sha256(evidence_path),
                            }
                            for role, evidence_path in (
                                ("runtime_contract", paired_contract_path),
                                ("candidate_csv", candidate_path),
                                ("trace_csv", trace_path),
                                ("response_trace_csv", response_trace_path),
                            )
                        },
                    },
                )
                candidate_links = [
                    {
                        "role": "candidate_rows",
                        "path": str(candidate_path),
                        "sha256": file_sha256(candidate_path),
                    }
                ]
                self.candidate_paths[lane_id] = candidate_path
                self.candidate_rows[lane_id] = [row]
                self.trace_bundles[lane_id] = {
                    "runtime_contract": paired_contract_path,
                    "candidate_csv": candidate_path,
                    "executable_aggregate": aggregate_path,
                    "trace_csv": trace_path,
                    "response_trace_csv": response_trace_path,
                }

            observables: dict[str, object] = {}
            family_spec = CONTRACT["families"][family]
            for offset, observable in enumerate(family_spec["required_observables"]):
                base = float(100 + index + offset)
                value = (
                    0.5 * base
                    if observable.endswith("_signal")
                    else 0.4 * base
                    if observable.endswith("_notmatch")
                    else base
                )
                observables[observable] = {
                    "bin_edges": [10.0, 12.0, 14.0],
                    "sumw": [value, value * 0.8],
                    "sumw2": [value * 0.2, value * 0.16],
                    "fills": [value * 10.0, value * 8.0],
                    "object_path": f"SIM/{EXTRACTOR.REGION_OBJECTS[observable]}",
                    "entries": value * 18.0,
                    "flow": {
                        "sumw": [0.0, 0.0],
                        "sumw2": [0.0, 0.0],
                        "fills": [0.0, 0.0],
                    },
                }
            fill_path = lane_evidence_dir / "fills.json"
            write_json(
                fill_path,
                {
                    "schema": "ppg12-stitched-purity-fill-evidence/v1",
                    "lane_id": lane_id,
                    "root_sha256": file_sha256(root_path),
                    "observables": {
                        name: {"object": name} for name in observables
                    },
                },
            )
            sidecar_path = lane_evidence_dir / "input_sidecar.json"
            write_json(
                sidecar_path,
                {"schema": "ppg12-stitched-purity-lane-input/v1", "lane_id": lane_id},
            )
            root_link = {"path": str(root_path), "sha256": file_sha256(root_path)}
            event_rows = [f"{lane_id}:source-{row}" for row in range(5)]
            groups = EXTRACTOR._canonical_groups(lane_id, event_rows, 5)
            group_set_sha256 = EXTRACTOR._payload_sha256(groups)
            runtime_evidence, runtime_source_sets, runtime_hashes = (
                EXTRACTOR._validate_runtime_manifest(
                    runtime_manifest_link,
                    lane_evidence_dir,
                    evidence_links["config"],
                    CONTRACT,
                )
            )
            extraction = {
                "schema": "ppg12-stitched-purity-lane-extraction/v1",
                "extractor": {
                    "path": str(EXTRACTOR_PATH.resolve()),
                    "sha256": file_sha256(EXTRACTOR_PATH),
                },
                "input_sidecar": {
                    "path": str(sidecar_path),
                    "sha256": file_sha256(sidecar_path),
                },
                "root": root_link,
                "fill_evidence": {
                    "path": str(fill_path),
                    "sha256": file_sha256(fill_path),
                },
                "evidence": evidence_links,
                "runtime_evidence": runtime_evidence,
                "runtime_source_sets": runtime_source_sets,
                "candidate_parity_evidence": candidate_links,
                "groups": groups,
                "group_set_sha256": group_set_sha256,
                "group_input_sidecars": [
                    {
                        "group_index": 0,
                        "group_id": groups[0]["group_id"],
                        "path": str(sidecar_path.resolve()),
                        "sha256": file_sha256(sidecar_path),
                    }
                ],
                "object_prefix": "SIM",
                "fill_count_semantics": "exact external per-bin counts checked against TH1 entries",
            }
            payload: dict[str, object] = {
                "schema": "ppg12-stitched-purity-lane/v1",
                "lane_id": lane_id,
                "family": family,
                "sample": sample,
                "period": period,
                "interaction": interaction,
                "group_count": 1,
                "group_size": 5,
                "group_index_start": 0,
                "event_set_row_count": 5,
                "groups": groups,
                "group_set_sha256": group_set_sha256,
                "abcd_population": "unsuffixed",
                "external_scale": 1.0,
                "observables": observables,
                "merge_input": {"input_id": lane_id, **root_link},
                "extraction": extraction,
                "candidate_parity_evidence": candidate_links,
            }
            for evidence_name, field in ASSEMBLER.LANE_DIRECT_EVIDENCE_TO_FIELD.items():
                payload[field] = evidence_links[evidence_name]["sha256"]
            payload.update(runtime_hashes)
            path = self.lane_dir / f"lane_{index:02d}.json"
            write_json(path, payload)
            self.lane_paths.append(path)
            self.lane_payloads[lane_id] = payload

    def _make_supporting_payloads(self) -> None:
        provenance: dict[str, str] = {}
        source_sets: dict[str, object] = {}
        provenance_dir = self.root / "provenance_sources"
        provenance_dir.mkdir(exist_ok=True)
        for field in CONTRACT["frozen_provenance_fields"]:
            source = provenance_dir / f"{field}.txt"
            source.write_text(f"canonical source for {field}\n")
            files = [{"path": str(source), "sha256": file_sha256(source)}]
            set_sha = ASSEMBLER._payload_sha256(files)
            provenance[field] = set_sha
            source_sets[field] = {"files": files, "set_sha256": set_sha}
        write_json(
            self.provenance_path,
            {
                "schema": "ppg12-stitched-purity-provenance/v1",
                "provenance": provenance,
                "source_sets": source_sets,
            },
        )
        lane_links = [
            {
                "lane_id": lane_id,
                "path": str(path.resolve()),
                "sha256": file_sha256(path),
            }
            for lane_id, path in sorted(zip(self.lane_payloads, self.lane_paths))
        ]
        reproduced_purity = ASSEMBLER._reproduce_purity_evidence(
            CONTRACT_PATH,
            CONTRACT,
            sorted(self.lane_payloads.values(), key=lambda row: row["lane_id"]),
            lane_links,
        )
        write_json(
            self.purity_path,
            reproduced_purity,
        )
        parity_payload = RUNNER.build_candidate_parity(
            self.candidate_rows,
            created_utc="2026-07-20T00:00:00+00:00",
            trace_bundle_by_lane=self.trace_bundles,
        )
        write_json(self.parity_path, parity_payload)

    def _make_index(self) -> int:
        argv = ["lane-index"]
        for path in self.lane_paths:
            argv.extend(["--lane-json", str(path)])
        argv.extend(["--output", str(self.index_path)])
        return self._run(argv)

    def _make_coverage(self) -> int:
        return self._run(
            [
                "coverage-evidence",
                "--lane-index",
                str(self.index_path),
                "--output",
                str(self.coverage_path),
            ]
        )

    def _make_manifest(self, role: str, output: Path) -> int:
        argv = [
            "manifest",
            "--role",
            role,
            "--lane-index",
            str(self.index_path),
            "--provenance-json",
            str(self.provenance_path),
            "--purity-json",
            str(self.purity_path),
            "--global-minimality-json",
            str(self.coverage_path),
            "--output",
            str(output),
        ]
        if role != "reference":
            argv.extend(["--candidate-parity-json", str(self.parity_path)])
        return self._run(argv)

    def _rewrite_lane(self, lane_id: str) -> None:
        position = list(self.lane_payloads).index(lane_id)
        write_json(self.lane_paths[position], self.lane_payloads[lane_id])

    def _make_family_summaries(self) -> tuple[Path, Path]:
        paths: list[Path] = []
        self.family_outputs: dict[str, Path] = {}
        for family in ("inclusive", "photon"):
            inputs = [
                Path(lane["merge_input"]["path"])
                for lane in sorted(
                    (
                        lane
                        for lane in self.lane_payloads.values()
                        if lane["family"] == family
                    ),
                    key=lambda lane: lane["lane_id"],
                )
            ]
            output_root = self.root / f"final_{family}.root"
            subprocess.run(
                [
                    str(REPO.parent / "analysis/env/bin/python3"),
                    "-c",
                    (
                        "import ROOT,sys; paths=sys.argv[2:]; "
                        "files=[ROOT.TFile.Open(p,'READ') for p in paths]; "
                        "out=ROOT.TFile(sys.argv[1],'RECREATE'); "
                        "h=files[0].Get('merge_hist').Clone('merge_hist'); h.SetDirectory(0); "
                        "[h.Add(f.Get('merge_hist')) for f in files[1:]]; "
                        "out.cd(); h.Write(); ROOT.TNamed('analysis_config_yaml','unit').Write(); "
                        "out.Close(); [f.Close() for f in files]"
                    ),
                    str(output_root),
                    *[str(path) for path in inputs],
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            path = self.root / f"{family}_summary.json"
            subprocess.run(
                [
                    str(REPO.parent / "analysis/env/bin/python3"),
                    str(REPO / CONTRACT["canonical_tools"]["merge_auditor"]),
                    "--family",
                    family,
                    "--output",
                    str(output_root),
                    *[item for input_path in inputs for item in ("--input", str(input_path))],
                    "--required-token",
                    "merge_hist",
                    "--json",
                    str(path),
                ],
                check=True,
                capture_output=True,
                text=True,
            )
            self.family_outputs[family] = output_root
            paths.append(path)
        return paths[0], paths[1]

    def test_manifest_assembly_is_exact_and_deterministic(self) -> None:
        self.assertEqual(self._make_manifest("reference", self.reference_path), 0)
        self.assertEqual(
            self._make_manifest("candidate", self.candidate_path), 0, self.last_output
        )
        first = self.candidate_path.read_bytes()
        self.assertEqual(self._make_manifest("candidate", self.candidate_path), 0)
        self.assertEqual(first, self.candidate_path.read_bytes())
        payload = json.loads(first)
        self.assertEqual(payload["schema"], "ppg12-stitched-purity-manifest/v1")
        self.assertEqual(len(payload["lanes"]), 32)
        self.assertEqual(payload["abcd_population"], "unsuffixed")
        self.assertEqual(payload["random_seed"], 42)
        self.assertEqual(payload["toy_count"], 20000)

    def test_global_minimality_is_replay_derived_and_rejects_stale_or_forged_artifacts(self) -> None:
        coverage = json.loads(self.coverage_path.read_text())
        self.assertEqual(coverage["baseline"]["group_count"], 32)
        self.assertEqual(coverage["added_groups"], [])
        self.assertEqual(coverage["final_missing_cells"], [])

        coverage["baseline"]["missing_cells"] = ["inclusive:A:bin-1"]
        coverage.pop("verification_payload_sha256")
        coverage["verification_payload_sha256"] = ASSEMBLER._payload_sha256(coverage)
        write_json(self.coverage_path, coverage)
        self.assertNotEqual(self._make_manifest("candidate", self.candidate_path), 0)

        self.assertEqual(self._make_coverage(), 0, self.last_output)
        self.assertEqual(self._make_manifest("candidate", self.candidate_path), 0)
        self.coverage_path.write_text(self.coverage_path.read_text() + "\n")
        with self.assertRaisesRegex(ASSEMBLER.AssemblyError, "hash mismatch"):
            ASSEMBLER.verify_manifest_evidence(
                self.candidate_path, CONTRACT, {"candidate"}
            )

    def test_global_prefix_uses_only_cell_reducing_groups_and_stops_when_complete(self) -> None:
        links = ASSEMBLER._load_lane_links(
            index_path=self.index_path, lane_paths=[], contract=CONTRACT
        )
        lanes, lane_links, bin_edges = ASSEMBLER._load_and_validate_lanes(
            links, CONTRACT, require_merge_input=False
        )
        groups = ASSEMBLER._canonical_group_extracts(lanes, CONTRACT)
        original_group = copy.deepcopy(
            next(row for row in groups if row["family"] == "inclusive")
        )
        for row in groups:
            if row["family"] == "inclusive":
                row["observables"]["A"]["sumw"][0] = 0.0
                row["observables"]["A"]["fills"][0] = 0.0
        extra = original_group
        extra["group_index"] = 1
        extra["group_id"] = "group-00001-useful"
        extra["group_sha256"] = digest("useful-extra")
        groups.append(extra)
        original_builder = ASSEMBLER._canonical_group_extracts
        ASSEMBLER._canonical_group_extracts = lambda *_args, **_kwargs: copy.deepcopy(groups)
        try:
            evidence = ASSEMBLER._build_global_minimality_evidence(
                lanes, lane_links, bin_edges, CONTRACT
            )
            self.assertEqual(len(evidence["added_groups"]), 1)
            self.assertEqual(
                evidence["added_groups"][0]["newly_populated"],
                ["inclusive:A:bin-1"],
            )

            useless = copy.deepcopy(extra)
            useless["lane_id"] = "inclusive:jet8:0mrad:di"
            useless["group_id"] = "group-00001-arbitrary"
            useless["group_sha256"] = digest("arbitrary-extra")
            groups.append(useless)
            with self.assertRaisesRegex(ASSEMBLER.AssemblyError, "arbitrary extras"):
                ASSEMBLER._build_global_minimality_evidence(
                    lanes, lane_links, bin_edges, CONTRACT
                )
        finally:
            ASSEMBLER._canonical_group_extracts = original_builder

    def test_index_hash_drift_fails_and_removes_stale_output(self) -> None:
        self.assertEqual(self._make_manifest("candidate", self.candidate_path), 0)
        self.lane_paths[0].write_text(self.lane_paths[0].read_text() + "\n")
        self.assertNotEqual(self._make_manifest("candidate", self.candidate_path), 0)
        self.assertFalse(self.candidate_path.exists())

    def test_missing_or_duplicate_lane_fails(self) -> None:
        original = self.lane_paths[-1]
        self.lane_paths[-1] = self.lane_paths[0]
        self.assertNotEqual(self._make_index(), 0)
        self.assertFalse(self.index_path.exists())
        self.lane_paths[-1] = original

    def test_classed_abcd_and_jet8_scale_fail(self) -> None:
        lane_id = "inclusive:jet8:0mrad:si"
        self.lane_payloads[lane_id]["abcd_population"] = "class-suffixed"
        self.lane_payloads[lane_id]["external_scale"] = 2.3398
        self._rewrite_lane(lane_id)
        self.assertNotEqual(self._make_index(), 0)

    def test_lane_requires_canonical_extractor_and_group_evidence(self) -> None:
        lane_id = "inclusive:jet8:0mrad:si"
        del self.lane_payloads[lane_id]["extraction"]["extractor"]
        self._rewrite_lane(lane_id)
        self.assertNotEqual(self._make_index(), 0)

        self.setUp_cleanup_rebuild()
        lane_id = "inclusive:jet8:0mrad:si"
        self.lane_payloads[lane_id]["groups"][0]["group_sha256"] = "0" * 64
        self._rewrite_lane(lane_id)
        self.assertNotEqual(self._make_index(), 0)

    def test_event_set_cannot_be_an_independent_hash_valid_list(self) -> None:
        lane_id = "inclusive:jet8:0mrad:si"
        lane = self.lane_payloads[lane_id]
        event_path = Path(lane["extraction"]["evidence"]["event_set"]["path"])
        rows = [f"{lane_id}:source-{row}" for row in (1, 0, 2, 3, 4)]
        event_path.write_text("\n".join(rows) + "\n")
        lane["event_set_sha256"] = file_sha256(event_path)
        lane["extraction"]["evidence"]["event_set"]["sha256"] = file_sha256(
            event_path
        )
        groups = EXTRACTOR._canonical_groups(lane_id, rows, 5)
        lane["groups"] = groups
        lane["group_set_sha256"] = EXTRACTOR._payload_sha256(groups)
        lane["extraction"]["groups"] = groups
        lane["extraction"]["group_set_sha256"] = lane["group_set_sha256"]
        lane["extraction"]["group_input_sidecars"][0]["group_id"] = groups[0][
            "group_id"
        ]
        self._rewrite_lane(lane_id)
        self.assertNotEqual(self._make_index(), 0)
        self.assertIn("production-list slice", self.last_output)

    def test_runtime_contract_hashes_cannot_be_replaced_by_caller_files(self) -> None:
        lane_id = "inclusive:jet8:0mrad:si"
        lane = self.lane_payloads[lane_id]
        model = next(
            row
            for row in lane["extraction"]["runtime_source_sets"]["model_set_sha256"][
                "files"
            ]
            if row["role"] == "ppg_apply_model_base_v3E"
        )
        Path(model["path"]).write_text("mutated executed model\n")
        lane["model_set_sha256"] = digest("caller supplied replacement")
        self._rewrite_lane(lane_id)
        self.assertNotEqual(self._make_index(), 0)
        self.assertIn("runtime evidence is invalid", self.last_output)

    def test_candidate_rejects_reused_physical_root_path_or_hash(self) -> None:
        first_id = "inclusive:jet8:0mrad:si"
        second_id = "inclusive:jet8:0mrad:di"
        first = self.lane_payloads[first_id]
        second = self.lane_payloads[second_id]
        first_link = copy.deepcopy(first["merge_input"])
        second["merge_input"] = {"input_id": second_id, **{
            "path": first_link["path"], "sha256": first_link["sha256"]
        }}
        second["extraction"]["root"] = {
            "path": first_link["path"], "sha256": first_link["sha256"]
        }
        fill_path = Path(second["extraction"]["fill_evidence"]["path"])
        fill = json.loads(fill_path.read_text())
        fill["root_sha256"] = first_link["sha256"]
        write_json(fill_path, fill)
        second["extraction"]["fill_evidence"]["sha256"] = file_sha256(fill_path)
        self._rewrite_lane(second_id)
        self._install_lane_replay_fixture()
        self.assertEqual(self._make_index(), 0, self.last_output)
        self.assertNotEqual(
            self._make_manifest("candidate", self.candidate_path), 0
        )
        self.assertIn("unique resolved ROOT paths", self.last_output)

        self.setUp_cleanup_rebuild()
        first = self.lane_payloads[first_id]
        second = self.lane_payloads[second_id]
        first_path = Path(first["merge_input"]["path"])
        second_path = Path(second["merge_input"]["path"])
        second_path.write_bytes(first_path.read_bytes())
        duplicate_sha = file_sha256(second_path)
        second["merge_input"]["sha256"] = duplicate_sha
        second["extraction"]["root"]["sha256"] = duplicate_sha
        fill_path = Path(second["extraction"]["fill_evidence"]["path"])
        fill = json.loads(fill_path.read_text())
        fill["root_sha256"] = duplicate_sha
        write_json(fill_path, fill)
        second["extraction"]["fill_evidence"]["sha256"] = file_sha256(fill_path)
        self._rewrite_lane(second_id)
        self._install_lane_replay_fixture()
        self.assertEqual(self._make_index(), 0, self.last_output)
        self.assertNotEqual(
            self._make_manifest("candidate", self.candidate_path), 0
        )
        self.assertIn("unique ROOT content hashes", self.last_output)

    def test_missing_observable_and_bin_mismatch_fail(self) -> None:
        lane_id = "photon:photon5:0mrad:si"
        del self.lane_payloads[lane_id]["observables"]["D_signal"]
        self._rewrite_lane(lane_id)
        self.assertNotEqual(self._make_index(), 0)

        self.setUp_cleanup_rebuild()
        lane_id = "inclusive:jet12:0mrad:di"
        self.lane_payloads[lane_id]["observables"]["A"]["bin_edges"] = [10.0, 11.0, 14.0]
        self._rewrite_lane(lane_id)
        self.assertNotEqual(self._make_index(), 0)

    def setUp_cleanup_rebuild(self) -> None:
        # Rebuild only lane fixtures inside the current TemporaryDirectory.
        for path in self.lane_paths:
            path.unlink(missing_ok=True)
        for path in self.artifact_dir.glob("*.root"):
            path.unlink()
        self.lane_paths = []
        self.lane_payloads = {}
        self._make_lanes()
        self._install_lane_replay_fixture()

    def test_bad_provenance_hash_and_seed_fail(self) -> None:
        provenance = json.loads(self.provenance_path.read_text())
        provenance["provenance"][CONTRACT["frozen_provenance_fields"][0]] = "not-a-hash"
        write_json(self.provenance_path, provenance)
        self.assertNotEqual(self._make_manifest("candidate", self.candidate_path), 0)

    def test_linked_provenance_purity_and_trace_bytes_are_revalidated(self) -> None:
        source = next((self.root / "provenance_sources").glob("*.txt"))
        source.write_text(source.read_text() + "mutated\n")
        self.assertNotEqual(self._make_manifest("candidate", self.candidate_path), 0)

        self._make_supporting_payloads()
        purity = json.loads(self.purity_path.read_text())
        purity["algorithm"]["producer"]["sha256"] = "0" * 64
        write_json(self.purity_path, purity)
        self.assertNotEqual(self._make_manifest("candidate", self.candidate_path), 0)

        self._make_supporting_payloads()
        trace = next(iter(self.candidate_paths.values()))
        trace.write_text(trace.read_text() + "\n")
        self.assertNotEqual(self._make_manifest("candidate", self.candidate_path), 0)

        self.setUp_cleanup_rebuild()
        self._make_supporting_payloads()
        purity = json.loads(self.purity_path.read_text())
        purity["toy_count"] = 1000
        write_json(self.purity_path, purity)
        self.assertNotEqual(self._make_manifest("candidate", self.candidate_path), 0)

    def test_self_hashed_purity_mutation_cannot_impersonate_executable_replay(self) -> None:
        purity = json.loads(self.purity_path.read_text())
        purity["purity"]["corrected"]["value"][0] += 0.123
        purity_sha = ASSEMBLER._payload_sha256(purity["purity"])
        purity["fixed_seed_repetition"] = {
            "first_output_sha256": purity_sha,
            "repeated_output_sha256": purity_sha,
        }
        write_json(self.purity_path, purity)
        self.assertNotEqual(self._make_manifest("candidate", self.candidate_path), 0)

    def test_candidate_parity_requires_all_features_and_both_models(self) -> None:
        parity = json.loads(self.parity_path.read_text())
        parity["candidate_parity"]["features"].pop()
        parity["candidate_parity"]["scores"].pop()
        write_json(self.parity_path, parity)
        self.assertNotEqual(self._make_manifest("candidate", self.candidate_path), 0)

    def test_candidate_parity_requires_exact_lane_coverage_and_tolerances(self) -> None:
        parity = json.loads(self.parity_path.read_text())
        parity["candidate_parity"]["lane_coverage"].pop()
        parity["candidate_parity"]["tolerances"]["feature_relative"] = 1e-4
        write_json(self.parity_path, parity)
        self.assertNotEqual(self._make_manifest("candidate", self.candidate_path), 0)

    def test_purity_fixed_seed_repetition_must_match_payload(self) -> None:
        purity = json.loads(self.purity_path.read_text())
        purity["fixed_seed_repetition"]["repeated_output_sha256"] = "0" * 64
        write_json(self.purity_path, purity)
        self.assertNotEqual(self._make_manifest("candidate", self.candidate_path), 0)

    def test_merge_audit_binds_exact_32_inputs_and_detects_file_drift(self) -> None:
        self.assertEqual(self._make_manifest("candidate", self.candidate_path), 0)
        inclusive_summary, photon_summary = self._make_family_summaries()
        argv = [
            "merge-audit",
            "--candidate-manifest",
            str(self.candidate_path),
            "--inclusive-summary",
            str(inclusive_summary),
            "--photon-summary",
            str(photon_summary),
            "--output",
            str(self.merge_path),
        ]
        self.assertEqual(self._run(argv), 0, self.last_output)
        payload = json.loads(self.merge_path.read_text())
        self.assertEqual(
            [row["inputs_fixed_order_count"] for row in payload["audits"]], [20, 12]
        )
        first = self.merge_path.read_bytes()
        self.assertEqual(self._run(argv), 0)
        self.assertEqual(first, self.merge_path.read_bytes())

        first_root = Path(self.lane_payloads["inclusive:jet8:0mrad:si"]["merge_input"]["path"])
        first_root.write_bytes(b"mutated ROOT")
        self.assertNotEqual(self._run(argv), 0)
        self.assertFalse(self.merge_path.exists())

    def test_forged_family_summary_or_final_root_drift_fails(self) -> None:
        self.assertEqual(self._make_manifest("candidate", self.candidate_path), 0)
        inclusive_summary, photon_summary = self._make_family_summaries()
        forged = json.loads(inclusive_summary.read_text())
        forged["output_bytes"] += 1
        unhashed = dict(forged)
        unhashed.pop("verification_payload_sha256")
        forged["verification_payload_sha256"] = ASSEMBLER._payload_sha256(unhashed)
        write_json(inclusive_summary, forged)
        argv = [
            "merge-audit",
            "--candidate-manifest",
            str(self.candidate_path),
            "--inclusive-summary",
            str(inclusive_summary),
            "--photon-summary",
            str(photon_summary),
            "--output",
            str(self.merge_path),
        ]
        self.assertNotEqual(self._run(argv), 0)

        self._make_family_summaries()
        self.family_outputs["photon"].write_bytes(b"forged final ROOT")
        self.assertNotEqual(self._run(argv), 0)

    def test_end_to_end_admission_and_production_wrapper(self) -> None:
        self.assertEqual(self._make_manifest("reference", self.reference_path), 0)
        self.assertEqual(
            self._make_manifest("candidate", self.candidate_path), 0, self.last_output
        )
        inclusive_summary, photon_summary = self._make_family_summaries()
        self.assertEqual(
            self._run(
                [
                    "merge-audit",
                    "--candidate-manifest",
                    str(self.candidate_path),
                    "--inclusive-summary",
                    str(inclusive_summary),
                    "--photon-summary",
                    str(photon_summary),
                    "--output",
                    str(self.merge_path),
                ]
            ),
            0,
        )
        gate_out = self.root / "gate"
        with contextlib.redirect_stdout(io.StringIO()):
            code = GATE.main(
                [
                    "--contract",
                    str(CONTRACT_PATH),
                    "admit",
                    "--reference-manifest",
                    str(self.reference_path),
                    "--candidate-manifest",
                    str(self.candidate_path),
                    "--merge-audit",
                    str(self.merge_path),
                    "--outdir",
                    str(gate_out),
                ]
            )
        self.assertEqual(code, 0)
        admission_path = gate_out / "admission_manifest.json"
        self.assertTrue(admission_path.exists())

        self.assertEqual(
            self._make_manifest("production", self.production_candidate_path), 0
        )
        production_merge_path = self.root / "production_merge_audit.json"
        self.assertEqual(
            self._run(
                [
                    "merge-audit",
                    "--candidate-manifest",
                    str(self.production_candidate_path),
                    "--inclusive-summary",
                    str(inclusive_summary),
                    "--photon-summary",
                    str(photon_summary),
                    "--output",
                    str(production_merge_path),
                ]
            ),
            0,
        )
        inclusive_root = self.family_outputs["inclusive"]
        photon_root = self.family_outputs["photon"]
        source_paths = {
            "historical_purity": self.root / "historical_purity.root",
            "historical_abcd": self.root / "historical_abcd.root",
            "historical_leakage": self.root / "historical_purity.root",
            "candidate_purity": self.purity_path,
            "candidate_inclusive": inclusive_root,
            "candidate_photon": photon_root,
        }
        for role in ("historical_purity", "historical_abcd"):
            source_paths[role].write_bytes(role.encode())
        source_links = [
            {
                "role": role,
                "path": str(path.resolve()),
                "sha256": file_sha256(path),
            }
            for role, path in sorted(source_paths.items())
        ]
        historical_path = self.root / "historical.json"
        write_json(
            historical_path,
            {
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
                "source_set_sha256": ASSEMBLER._payload_sha256(source_links),
            },
        )
        wrapper_path = self.root / "production_wrapper.json"
        argv = [
            "production-wrapper",
            "--admission-manifest",
            str(admission_path),
            "--reference-manifest",
            str(self.reference_path),
            "--candidate-manifest",
            str(self.production_candidate_path),
            "--merge-audit",
            str(production_merge_path),
            "--inclusive-artifact",
            str(inclusive_root),
            "--photon-artifact",
            str(photon_root),
            "--historical-comparison",
            str(historical_path),
            "--output",
            str(wrapper_path),
        ]
        self.assertEqual(self._run(argv), 0, self.last_output)
        wrapper = json.loads(wrapper_path.read_text())
        self.assertEqual(wrapper["admission_sha256"], file_sha256(admission_path))
        self.assertEqual(
            {row["family"] for row in wrapper["candidate_artifacts"]},
            {"inclusive", "photon"},
        )

        # A caller cannot fabricate authority by copying a PASS-shaped
        # admission, mutating frozen physics, and reserializing the JSON.
        forged_admission_path = self.root / "forged_admission.json"
        forged_admission = json.loads(admission_path.read_text())
        forged_lane = sorted(forged_admission["candidate_frozen"]["lanes"])[0]
        forged_admission["candidate_frozen"]["lanes"][forged_lane][
            "model_set_sha256"
        ] = digest("self-rehashed forged model")
        write_json(forged_admission_path, forged_admission)
        forged_wrapper_path = self.root / "forged_production_wrapper.json"
        forged_argv = list(argv)
        forged_argv[forged_argv.index(str(admission_path))] = str(
            forged_admission_path
        )
        forged_argv[forged_argv.index(str(wrapper_path))] = str(forged_wrapper_path)
        self.assertNotEqual(self._run(forged_argv), 0)
        self.assertIn("exact output of the canonical pair gate", self.last_output)
        self.assertFalse(forged_wrapper_path.exists())

        production_out = self.root / "production_gate"
        with contextlib.redirect_stdout(io.StringIO()):
            code = GATE.main(
                [
                    "--contract",
                    str(CONTRACT_PATH),
                    "verify-production",
                    "--admission-manifest",
                    str(admission_path),
                    "--production-manifest",
                    str(wrapper_path),
                    "--merge-audit",
                    str(production_merge_path),
                    "--outdir",
                    str(production_out),
                ]
            )
        self.assertEqual(code, 0)
        self.assertTrue((production_out / "production_gate_report.json").exists())

        # Stale admission/physics hashes fail before a wrapper can be emitted.
        production = json.loads(self.production_candidate_path.read_text())
        production["provenance"]["model_set_sha256"] = digest("changed model")
        write_json(self.production_candidate_path, production)
        self.assertNotEqual(self._run(argv), 0)
        self.assertFalse(wrapper_path.exists())


if __name__ == "__main__":
    unittest.main()
