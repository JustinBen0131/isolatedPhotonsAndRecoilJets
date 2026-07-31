#!/usr/bin/env python3
"""Adversarial tests for the non-submitting THE-134 full controller."""

from __future__ import annotations

import argparse
import ast
import copy
import hashlib
import importlib.util
import json
import os
import subprocess
import sys
import tempfile
import time
import unittest
from pathlib import Path
from unittest import mock


HERE = Path(__file__).resolve().parent
CONTROLLER_PATH = (
    HERE.parent / "materialize_the134_full_multiview_extraction.py"
)
SPEC = importlib.util.spec_from_file_location(
    "the134_full_materializer", CONTROLLER_PATH
)
assert SPEC is not None and SPEC.loader is not None
controller = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(controller)
resolver = controller.resolver
projector = controller.projector
AMENDMENT_TEST_PATH = (
    HERE / "test_build_the134_capacity_count_amendment.py"
)
AMENDMENT_TEST_SPEC = importlib.util.spec_from_file_location(
    "the134_capacity_amendment_fixture_builder", AMENDMENT_TEST_PATH
)
assert (
    AMENDMENT_TEST_SPEC is not None
    and AMENDMENT_TEST_SPEC.loader is not None
)
amendment_test_module = importlib.util.module_from_spec(AMENDMENT_TEST_SPEC)
AMENDMENT_TEST_SPEC.loader.exec_module(amendment_test_module)
BUDGET_BUILDER_PATH = (
    HERE.parent / "build_the134_controller_dry_materialization_budget.py"
)
BUDGET_BUILDER_SPEC = importlib.util.spec_from_file_location(
    "the134_controller_budget_fixture_builder", BUDGET_BUILDER_PATH
)
assert (
    BUDGET_BUILDER_SPEC is not None
    and BUDGET_BUILDER_SPEC.loader is not None
)
budget_builder = importlib.util.module_from_spec(BUDGET_BUILDER_SPEC)
BUDGET_BUILDER_SPEC.loader.exec_module(budget_builder)


def write_json(path: Path, payload: object) -> None:
    path.write_bytes(controller.canonical_json_bytes(payload))


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class FullControllerFixture:
    def __init__(self, root: Path):
        self.root = root
        self.root.mkdir(parents=True)
        self.inputs = root / "inputs"
        self.inputs.mkdir()
        self.release_lib = root / "release_lib"
        self.release_lib64 = root / "release_lib64"
        self.release_lib.mkdir()
        self.release_lib64.mkdir()
        self.artifact_root = root / "bundle_artifacts"
        self.artifact_root.mkdir()
        self.remote_root = root / "remote"
        self.remote_root.mkdir()
        self.tag = "the134_full_controller_fixture_v1"
        self.campaign = {
            "tag": self.tag,
            "output_root": str(self.remote_root / "output" / self.tag),
            "evidence_root": str(self.remote_root / "evidence" / self.tag),
            "submit_root": str(self.remote_root / "submit" / self.tag),
        }
        self.inventory_path = root / "bundle_inventory.json"
        self.bundle_parent = root / "immutable_bundles"
        self.bundle_parent.mkdir()
        self.source_path = root / "sources.json"
        self.partition_path = root / "the134_full_extraction_partition.jsonl"
        self.plan_path = root / "the134_full_extraction_plan.json"
        self.receipt_path = root / "preflight_receipt.json"
        self.budget_path = root / "controller_budget.json"
        self.controller_derivation_path = (
            root / "controller_budget_derivation.json"
        )
        self.storage_manifest_path = root / "storage_manifest.json"
        self.quota_snapshot_path = root / "quota_snapshot.json"
        self.storage_certificate_path = root / "storage_certificate.json"
        self.amendment_path = root / "capacity_amendment.json"
        self.readback_path = root / "capacity_readback.json"
        write_json(self.amendment_path, {"fixture": "capacity amendment"})
        write_json(self.readback_path, {"fixture": "capacity readback"})
        self.bundle, self.by_role, self.materialization = self.make_bundle()
        (
            self.sources,
            self.source_rows,
            self.chunks,
            self.row_partition_data,
        ) = self.make_sources_and_partition()
        write_json(self.source_path, self.sources)
        self.write_partition()
        self.plan = self.make_plan()
        self.write_plan_and_receipt()
        self.write_exact_budget_and_derivation()
        self.rebuild_storage_evidence()
        self.baseline_plan = copy.deepcopy(self.plan)
        self.baseline_budget = json.loads(self.budget_path.read_text())
        self.baseline_controller_derivation = json.loads(
            self.controller_derivation_path.read_text()
        )
        self.baseline_materialization = copy.deepcopy(self.materialization)
        self.baseline_sources = copy.deepcopy(self.sources)

    def make_bundle(self) -> tuple[dict, dict[str, dict], dict]:
        materializer = resolver.bundle_materializer
        builder = materializer.builder
        roles = sorted(builder.REQUIRED_ARTIFACT_ROLES)
        declared_hashes = {
            "code_sha256": "d" * 64,
            "replay_schema_sha256": "e" * 64,
            "training_schema_sha256": "f" * 64,
            "semantic_sha256": "1" * 64,
        }
        records = []
        for index, role in enumerate(roles):
            path = self.artifact_root / f"{role}.artifact"
            if role == "code_manifest":
                path.write_text(
                    json.dumps(declared_hashes, sort_keys=True) + "\n",
                    encoding="utf-8",
                )
            else:
                path.write_text(
                    f"{role}|{index}|fixture\n", encoding="utf-8"
                )
            if role in builder.EXECUTABLE_ROLES:
                path.chmod(0o755)
            records.append(
                {
                    "role": role,
                    "path": str(path),
                    "sha256": sha(path),
                    "size_bytes": path.stat().st_size,
                }
            )
        source_by_role = {record["role"]: record for record in records}
        dependencies = []
        for role in sorted(builder.REQUIRED_DEPENDENCY_PROVIDER_ROLES):
            if role.startswith("pp_"):
                consumers = ["pp_executor"]
            elif role.startswith("auau_"):
                consumers = ["auau_executor"]
            else:
                consumers = ["auau_executor", "pp_executor"]
            dependencies.append(
                {
                    "provider_role": role,
                    "consumer_roles": consumers,
                    "kind": "runtime_or_contract",
                    "required": True,
                    "sha256": source_by_role[role]["sha256"],
                }
            )
        spec = {
            "schema": builder.SPEC_SCHEMA,
            "public_commit": "a" * 40,
            "bundle_parent": str(self.bundle_parent),
            "runtime": {
                "release": builder.RELEASE,
                "offline_main": "/cvmfs/sphenix/ana.560",
                "calo_reco_soname": builder.CALO_RECO_SONAME,
                "request_memory_mb": builder.REQUEST_MEMORY_MB,
                "release_core_lib_dir": str(self.release_lib),
                "release_core_lib64_dir": str(self.release_lib64),
            },
            "declared_hashes": declared_hashes,
            "artifacts": records,
            "dependencies": dependencies,
            "hash_bindings": [
                {
                    "name": name,
                    "artifact_role": "code_manifest",
                    "mode": "json_pointer",
                    "pointer": f"/{name}",
                }
                for name in builder.DECLARED_HASH_NAMES
            ],
        }
        inventory = builder.build_receipt(spec)
        self.inventory_path.write_bytes(builder.canonical_json_bytes(inventory))
        summary = materializer.materialize(
            self.inventory_path, sha(self.inventory_path)
        )
        self.bundle_path = Path(summary["resolver_bundle_receipt"])
        self.materialization_path = Path(
            summary["materialization_receipt"]
        )
        payload = json.loads(self.bundle_path.read_text(encoding="utf-8"))
        materialization = json.loads(
            self.materialization_path.read_text(encoding="utf-8")
        )
        by_role = {
            record["role"]: record for record in payload["artifacts"]
        }
        return payload, by_role, materialization

    def make_sources_and_partition(
        self,
    ) -> tuple[dict, dict[str, dict], list[dict], dict[str, dict]]:
        source_rows: dict[str, dict] = {}
        chunks: list[dict] = []
        row_partition_data: dict[str, dict] = {}
        inventory = resolver.inventory_rows()
        global_index = 0
        for row_index, row in enumerate(inventory):
            tuple_count = 10_000 if row_index < 12 else 9_998
            tuple_record_sha = controller.canonical_sha256(
                {"row_id": row["row_id"], "tuple_count": tuple_count}
            )
            full_source_sha = controller.canonical_sha256(
                {"source": row["row_id"], "count": tuple_count}
            )
            first_sha = controller.canonical_sha256(
                {"source": row["row_id"], "tuple": 0}
            )
            last_sha = controller.canonical_sha256(
                {"source": row["row_id"], "tuple": tuple_count - 1}
            )
            list_dir = self.inputs / row["sample"]
            input_lists = [
                {
                    "role": role,
                    "path": str(list_dir / resolver.LIST_FILENAMES[role]),
                    "sha256": controller.canonical_sha256(
                        {"row": row["row_id"], "role": role}
                    ),
                    "line_count": tuple_count,
                    "executable_count": tuple_count,
                }
                for role in resolver.LIST_ROLES
            ]
            source_rows[row["row_id"]] = {
                "row_id": row["row_id"],
                "system": row["system"],
                "sample": row["sample"],
                "input_lists": input_lists,
                "tuple_count": tuple_count,
                "tuple_records_sha256": tuple_record_sha,
                "full_source_manifest_sha256": full_source_sha,
                "first_tuple_sha256": first_sha,
                "last_tuple_sha256": last_sha,
            }
            row_chunks = []
            for chunk_index, start in enumerate(
                range(0, tuple_count, controller.EXPECTED_GROUP_SIZE)
            ):
                stop = min(
                    start + controller.EXPECTED_GROUP_SIZE, tuple_count
                )
                tuple_hashes = [
                    controller.canonical_sha256(
                        {"row": row["row_id"], "tuple": index}
                    )
                    for index in range(start, stop)
                ]
                base = {
                    "schema": resolver.CHUNK_SCHEMA,
                    "row_id": row["row_id"],
                    "full_source_manifest_sha256": full_source_sha,
                    "group_size": controller.EXPECTED_GROUP_SIZE,
                    "chunk_index": chunk_index,
                    "segment": chunk_index + 1,
                    "tuple_index_start": start,
                    "tuple_index_end_exclusive": stop,
                    "tuple_count": stop - start,
                    "physical_lines": list(range(start + 1, stop + 1)),
                    "tuple_input_sha256s": tuple_hashes,
                    "tuple_records_sha256": controller.canonical_sha256(
                        {
                            "row": row["row_id"],
                            "start": start,
                            "stop": stop,
                            "tuples": tuple_hashes,
                        }
                    ),
                }
                base["chunk_fingerprint_sha256"] = (
                    controller.canonical_sha256(base)
                )
                record = {**base, "global_chunk_index": global_index}
                record["execution_chunk_sha256"] = (
                    controller.canonical_sha256(record)
                )
                row_chunks.append(record)
                chunks.append(record)
                global_index += 1
            source_chunks = [
                {
                    key: value
                    for key, value in record.items()
                    if key
                    not in {"global_chunk_index", "execution_chunk_sha256"}
                }
                for record in row_chunks
            ]
            row_partition_data[row["row_id"]] = {
                "job_count": len(row_chunks),
                "chunk_records_sha256": controller.canonical_sha256(
                    source_chunks
                ),
                "row_partition_sha256": controller.canonical_sha256(
                    {
                        "row_id": row["row_id"],
                        "tuple_count": tuple_count,
                        "chunk_fingerprints": [
                            record["chunk_fingerprint_sha256"]
                            for record in source_chunks
                        ],
                    }
                ),
            }
        self.assert_fixture_counts(source_rows, chunks)
        sources = {
            "schema": resolver.SOURCE_SCHEMA,
            "status": "PASS",
            "authority": {
                "scope": resolver.SOURCE_AUTHORITY_SCOPE,
                "row_count": resolver.SOURCE_AUTHORITY_ROW_COUNT,
                "pp_period": "0mrad",
                "pp_si_di_role": "SI",
                "auau_period": resolver.SOURCE_AUTHORITY_AUAU_PERIOD,
                "auau_si_di_role": (
                    resolver.SOURCE_AUTHORITY_AUAU_SI_DI_ROLE
                ),
                "source_ownership_state": (
                    resolver.SOURCE_AUTHORITY_OWNERSHIP_STATE
                ),
                "diagnostic_sources_excluded": list(
                    resolver.SOURCE_AUTHORITY_DIAGNOSTIC_EXCLUSIONS
                ),
                "scientific_completion_granted": False,
            },
            "rows": [source_rows[row["row_id"]] for row in inventory],
        }
        return sources, source_rows, chunks, row_partition_data

    @staticmethod
    def assert_fixture_counts(
        source_rows: dict[str, dict], chunks: list[dict]
    ) -> None:
        assert len(source_rows) == controller.EXPECTED_ROW_COUNT
        assert (
            sum(row["tuple_count"] for row in source_rows.values())
            == controller.EXPECTED_SOURCE_TUPLES
        )
        assert len(chunks) == controller.EXPECTED_JOB_COUNT

    def write_partition(self) -> None:
        with self.partition_path.open("wb") as stream:
            for record in self.chunks:
                stream.write(controller.canonical_json_bytes(record))

    def normalized_bundle(self) -> dict:
        return {
            "public_commit": self.bundle["public_commit"],
            "bundle_identity_sha256": self.bundle[
                "bundle_identity_sha256"
            ],
            "semantic_fingerprint_sha256": self.bundle[
                "semantic_fingerprint_sha256"
            ],
            "code_sha256": self.bundle["code_sha256"],
            "replay_schema_sha256": self.bundle["replay_schema_sha256"],
            "training_schema_sha256": self.bundle[
                "training_schema_sha256"
            ],
            "semantic_sha256": self.bundle["semantic_sha256"],
            "runtime": self.bundle["runtime"],
            "artifact_by_role": self.by_role,
            "artifacts": list(self.by_role.values()),
        }

    def make_row(self, inventory_row: dict) -> dict:
        row_id = inventory_row["row_id"]
        source = self.source_rows[row_id]
        partition = self.row_partition_data[row_id]
        system = inventory_row["system"]
        sidecar = (
            f"{self.campaign['output_root']}/{row_id}/training_views/"
            "$(Cluster).$(Process).root"
        )
        row = {
            "schema": resolver.ROW_SCHEMA,
            **inventory_row,
            "source_period": "0mrad" if system == "pp" else "AUAU_RUN24",
            "source_si_di_role": "SI" if system == "pp" else "EMBEDDED",
            "training_period_si_contract_sha256": (
                controller.canonical_sha256(
                    resolver.training_period_si_contract("0mrad")
                )
                if system == "pp"
                else None
            ),
            "requested_scope": "full",
            "full_training_authority": 0,
            "authority_state": resolver.PREFLIGHT_AUTHORITY_STATE,
            "input_contract": {
                "group_size": 7,
                "event_limit_per_job": 0,
                "source_tuple_count": source["tuple_count"],
                "expected_chunk_count": partition["job_count"],
                "expected_job_count": partition["job_count"],
                "expected_output_pair_count": partition["job_count"],
                "expected_analysis_output_count": partition["job_count"],
                "expected_sidecar_output_count": partition["job_count"],
                "expected_source_occurrence_count": partition["job_count"],
                "source_occurrences_per_output_pair": 1,
                "row_partition_sha256": partition[
                    "row_partition_sha256"
                ],
                "chunk_records_sha256": partition[
                    "chunk_records_sha256"
                ],
                "partition_artifact_name": self.partition_path.name,
                "partition_artifact_sha256": sha(self.partition_path),
                "full_source_manifest_sha256": source[
                    "full_source_manifest_sha256"
                ],
                "tuple_records_sha256": source["tuple_records_sha256"],
                "first_tuple_sha256": source["first_tuple_sha256"],
                "last_tuple_sha256": source["last_tuple_sha256"],
                "input_lists": source["input_lists"],
            },
            "science_contract": {
                "nominal_model_domain_gev": [15.0, 35.0],
                "loose_capture_et_min_gev": 5.0,
                "extraction_cone_r": 0.4,
                "legacy_training_tree_max_entries": 0,
                "shower_views": list(resolver.SHOWER_VIEWS),
            },
            "bundle_contract": {
                "public_commit": self.bundle["public_commit"],
                "code_sha256": self.bundle["code_sha256"],
                "replay_schema_sha256": self.bundle[
                    "replay_schema_sha256"
                ],
                "training_schema_sha256": self.bundle[
                    "training_schema_sha256"
                ],
                "semantic_sha256": self.bundle["semantic_sha256"],
                "library": self.by_role[f"{system}_library"],
                "model": self.by_role[f"{system}_model"],
                "config": self.by_role[f"{system}_config"],
                "submitter": self.by_role["submitter"],
                "executor": self.by_role[f"{system}_executor"],
            },
            "execution_contract": {
                "schema": resolver.EXECUTION_SCHEMA,
                "state": "PREFLIGHT_ONLY_NO_CONDOR_MUTATION",
                "existing_submitter_argv": [
                    self.by_role["submitter"]["path"],
                    inventory_row["dataset"],
                    "condorDoAll",
                    "groupSize",
                    "7",
                    f"SAMPLE={inventory_row['sample']}",
                ],
                "materialization_requires": "RJ_DAG_DRYRUN=1",
                "submission_requires": (
                    "separate explicit duplicate-guarded campaign authorization"
                ),
                "required_execution_inputs": [
                    "RJ_CODEX_CHAT_NAME",
                    "RJ_CODEX_THREAD_ID",
                ],
                "partition_artifact_name": self.partition_path.name,
                "partition_artifact_sha256": sha(self.partition_path),
                "row_partition_sha256": partition[
                    "row_partition_sha256"
                ],
                "analysis_output_namespace": (
                    f"{self.campaign['output_root']}/{row_id}"
                ),
                "multiview_sidecar_template": sidecar,
                "submit_namespace": (
                    f"{self.campaign['submit_root']}/{row_id}"
                ),
                "evidence_namespace": (
                    f"{self.campaign['evidence_root']}/{row_id}"
                ),
                "artifact_profile": copy.deepcopy(
                    resolver.SIDECAR_ONLY_ARTIFACT_PROFILE
                ),
                "materialization_environment": {},
                "worker_environment": {},
                "forbidden_ambient_environment": [
                    "RJ_PPG12_CROSSING_PERIOD",
                    "RJ_PPG12_PHOTON_YIELD_MIX_WEIGHT",
                    "RJ_PP_VERTEX_REWEIGHT_FILE",
                    "RJ_PPG12_PHOTON_YIELD_TRUTH_VERTEX",
                    "RJ_PPG12_PHOTON_YIELD_BUILDER_TRUTH_VERTEX",
                    "RJ_PPG12_PHOTON_YIELD_RECO_TRUTH_VERTEX",
                    "RJ_THE134_TRUTH_LABEL_DIAGNOSTIC_V1",
                    "RJ_THE134_TRUTH_LABEL_COUNTER_DIAGNOSTIC_V1",
                    "RJ_THE134_TRUTH_LABEL_PREWEIGHT_DIAGNOSTIC_V1",
                ],
            },
        }
        worker = controller.expected_worker_environment(
            row,
            campaign=self.campaign,
            training_period="0mrad",
            bundle=self.normalized_bundle(),
        )
        row["execution_contract"]["worker_environment"] = worker
        row["execution_contract"]["materialization_environment"] = (
            controller.expected_materialization_environment(
                row,
                campaign=self.campaign,
                bundle=self.normalized_bundle(),
                by_role=self.by_role,
                worker_environment=worker,
            )
        )
        row["row_fingerprint_sha256"] = controller.canonical_sha256(row)
        return row

    def make_plan(self) -> dict:
        rows = [self.make_row(row) for row in resolver.inventory_rows()]
        training = resolver.training_period_si_contract("0mrad")
        duplicate = {
            "schema": resolver.DUPLICATE_SCHEMA,
            "purpose": "THE134_SOURCE_COMPLETE_FACTORIAL_VIEW_TRAINING_EXTRACTION",
            "requested_scope": "full",
            "full_training_authority": 0,
            "authority_state": resolver.PREFLIGHT_AUTHORITY_STATE,
            "training_period_si_contract_sha256": (
                controller.canonical_sha256(training)
            ),
            "science_contract": {
                "model_domain_gev": [15.0, 35.0],
                "capture_et_min_gev": 5.0,
                "extraction_cone_r": 0.4,
                "shower_views": list(resolver.SHOWER_VIEWS),
            },
            "bundle_semantic_fingerprint_sha256": self.bundle[
                "semantic_fingerprint_sha256"
            ],
            "sources": [
                {
                    "row_id": row["row_id"],
                    "system": row["system"],
                    "sample": row["sample"],
                    "tuple_count": row["input_contract"][
                        "source_tuple_count"
                    ],
                    "tuple_records_sha256": row["input_contract"][
                        "tuple_records_sha256"
                    ],
                    "full_source_manifest_sha256": row["input_contract"][
                        "full_source_manifest_sha256"
                    ],
                }
                for row in rows
            ],
        }
        partition = {
            "schema": resolver.PARTITION_SCHEMA,
            "group_size": 7,
            "source_tuple_count": 129_998,
            "expected_chunk_count": 18_577,
            "expected_job_count": 18_577,
            "expected_output_pair_count": 18_577,
            "expected_analysis_output_count": 18_577,
            "expected_sidecar_output_count": 18_577,
            "expected_physical_root_artifact_count": 37_154,
            "expected_retained_analysis_output_count": 0,
            "expected_durable_root_artifact_count": 18_577,
            "expected_source_occurrence_count": 18_577,
            "source_occurrences_per_output_pair": 1,
            "row_partition_records_sha256": controller.canonical_sha256(
                [row["input_contract"] for row in rows]
            ),
            "row_partition_sha256s": [
                row["input_contract"]["row_partition_sha256"] for row in rows
            ],
            "partition_artifact": {
                "name": self.partition_path.name,
                "sha256": sha(self.partition_path),
                "record_count": 18_577,
            },
            "basis": (
                "ordered_disjoint_seven_tuple_partition_one_execution_per_chunk"
            ),
            "capacity_canary_required_before_submission": True,
            "capacity_authority_earned": False,
        }
        plan = {
            "schema": resolver.PLAN_SCHEMA,
            "status": "PREFLIGHT_PASS",
            "execution_state": "PREFLIGHT_ONLY_NO_CONDOR_MUTATION",
            "submission_performed": False,
            "campaign": self.campaign,
            "training_period_si_contract": training,
            "training_period_si_contract_sha256": (
                controller.canonical_sha256(training)
            ),
            "authority": {
                "requested_scope": "full",
                "full_training_authority": 0,
                "authority_state": resolver.PREFLIGHT_AUTHORITY_STATE,
            },
            "artifact_profile": copy.deepcopy(
                resolver.SIDECAR_ONLY_ARTIFACT_PROFILE
            ),
            "closure_witness_boundary": resolver.closure_witness_boundary(),
            "source_family_closure": {
                "row_count": 13,
                "pp_signal": 3,
                "pp_background": 4,
                "auau_signal": 2,
                "auau_background": 4,
                "pp_jet40_excluded_as_diagnostic_only": True,
            },
            "execution_partition": partition,
            "input_manifests": {
                "materialization": {
                    "path": str(self.materialization_path),
                    "sha256": sha(self.materialization_path),
                    "bundle_identity_sha256": self.bundle[
                        "bundle_identity_sha256"
                    ],
                    "digest_named_bundle_path": self.materialization[
                        "digest_named_bundle_path"
                    ],
                    "readback": "PASS_SYMLINK_FREE_READONLY_CONTENT_EXACT",
                },
                "bundle": {
                    "path": str(self.bundle_path),
                    "sha256": sha(self.bundle_path),
                    "semantic_fingerprint_sha256": self.bundle[
                        "semantic_fingerprint_sha256"
                    ],
                },
                "sources": {
                    "path": str(self.source_path),
                    "sha256": sha(self.source_path),
                },
            },
            "duplicate_contract": duplicate,
            "duplicate_fingerprint_sha256": (
                controller.canonical_sha256(duplicate)
            ),
            "execution_fingerprint_sha256": "",
            "rows": rows,
        }
        plan["execution_fingerprint_sha256"] = self.execution_fingerprint(plan)
        return plan

    def execution_fingerprint(self, plan: dict) -> str:
        return controller.canonical_sha256(
            {
                "schema": resolver.EXECUTION_SCHEMA,
                "tag": self.campaign["tag"],
                "output_root": self.campaign["output_root"],
                "evidence_root": self.campaign["evidence_root"],
                "submit_root": self.campaign["submit_root"],
                "materialization_receipt_sha256": sha(
                    self.materialization_path
                ),
                "bundle_manifest_sha256": sha(self.bundle_path),
                "source_manifest_sha256": sha(self.source_path),
                "duplicate_fingerprint_sha256": plan[
                    "duplicate_fingerprint_sha256"
                ],
                "partition_artifact_sha256": sha(self.partition_path),
                "execution_partition_sha256": (
                    controller.canonical_sha256(
                        plan["execution_partition"]
                    )
                ),
                "row_fingerprints": [
                    row["row_fingerprint_sha256"] for row in plan["rows"]
                ],
            }
        )

    def write_plan_and_receipt(self) -> None:
        write_json(self.plan_path, self.plan)
        receipt = {
            "schema": resolver.RECEIPT_SCHEMA,
            "status": "PASS",
            "submission_performed": False,
            "row_count": 13,
            "requested_scope": "full",
            "full_training_authority": 0,
            "authority_state": resolver.PREFLIGHT_AUTHORITY_STATE,
            "artifact_profile_sha256": controller.canonical_sha256(
                self.plan["artifact_profile"]
            ),
            "training_period_si_contract_sha256": self.plan[
                "training_period_si_contract_sha256"
            ],
            "closure_witness_boundary_sha256": (
                controller.canonical_sha256(
                    self.plan["closure_witness_boundary"]
                )
            ),
            "execution_partition_sha256": controller.canonical_sha256(
                self.plan["execution_partition"]
            ),
            "bundle_manifest_sha256": sha(self.bundle_path),
            "materialization_receipt_sha256": sha(
                self.materialization_path
            ),
            "source_manifest_sha256": sha(self.source_path),
            "duplicate_fingerprint_sha256": self.plan[
                "duplicate_fingerprint_sha256"
            ],
            "execution_fingerprint_sha256": self.plan[
                "execution_fingerprint_sha256"
            ],
            "artifacts": {
                "plan": {
                    "name": self.plan_path.name,
                    "sha256": sha(self.plan_path),
                },
                "rows": {
                    "name": "the134_full_extraction_rows.jsonl",
                    "sha256": "2" * 64,
                },
                "duplicate_fingerprint": {
                    "name": "the134_full_extraction_duplicate_fingerprint.sha256",
                    "sha256": "3" * 64,
                },
                "partition": {
                    "name": self.partition_path.name,
                    "sha256": sha(self.partition_path),
                    "record_count": 18_577,
                },
            },
        }
        write_json(self.receipt_path, receipt)

    def write_budget(
        self,
        *,
        fixed_bytes: int,
        bytes_per_job: int,
        fixed_inodes: int = 10,
        inodes_per_job: int = 1,
    ) -> None:
        payload = {
            "schema": projector.CONTROLLER_BUDGET_SCHEMA,
            "status": "PASS",
            "submission_performed": False,
            "expected_job_count": controller.EXPECTED_JOB_COUNT,
            "storage_budget": {
                "fixed_bytes": fixed_bytes,
                "fixed_inodes": fixed_inodes,
                "bytes_per_job": bytes_per_job,
                "inodes_per_job": inodes_per_job,
            },
            "full_training_authority": 0,
            "full_extraction_authority": False,
        }
        write_json(self.budget_path, payload)

    def write_exact_budget_and_derivation(self) -> None:
        plan, plan_artifact = budget_builder.load_pinned_json(
            self.plan_path,
            sha(self.plan_path),
            "fixture extraction plan",
        )
        receipt, receipt_artifact = budget_builder.load_pinned_json(
            self.receipt_path,
            sha(self.receipt_path),
            "fixture preflight receipt",
        )
        context = budget_builder.validate_input_chain(
            plan=plan,
            plan_path=self.plan_path,
            plan_artifact=plan_artifact,
            receipt=receipt,
            receipt_artifact=receipt_artifact,
        )
        budget, derivation = budget_builder.derive_budget(
            context,
            budget_output=self.budget_path,
            future_storage_certificate=self.storage_certificate_path,
        )
        write_json(self.controller_derivation_path, derivation)
        write_json(self.budget_path, budget)

    @staticmethod
    def repin_amendment_record(record: dict, path: Path) -> dict:
        updated = copy.deepcopy(record)
        updated.update(
            {
                "path": str(path),
                "resolved_path": str(path.resolve(strict=True)),
                "sha256": sha(path),
                "size_bytes": path.stat().st_size,
            }
        )
        return updated

    def write_valid_amendment_chain(self) -> None:
        fixture_root = self.root / "amendment_fixture"
        fixture_root.mkdir(exist_ok=True)
        helper = amendment_test_module.CapacityCountAmendmentTests(
            methodName="runTest"
        )
        helper.root = fixture_root
        payload = helper.make_valid_amendment_payload()
        corrected = payload["corrected_preflight"]
        immutable = payload["immutable_authority"]
        corrected["plan"] = self.repin_amendment_record(
            corrected["plan"], self.plan_path
        )
        corrected["preflight_receipt"] = self.repin_amendment_record(
            corrected["preflight_receipt"], self.receipt_path
        )
        corrected["partition"] = self.repin_amendment_record(
            corrected["partition"], self.partition_path
        )
        corrected["partition"]["record_count"] = (
            controller.EXPECTED_JOB_COUNT
        )
        corrected["bundle_manifest_sha256"] = sha(self.bundle_path)
        corrected["materialization_receipt_sha256"] = sha(
            self.materialization_path
        )
        payload["legacy_preflight"]["bundle_manifest_sha256"] = (
            corrected["bundle_manifest_sha256"]
        )
        payload["legacy_preflight"][
            "materialization_receipt_sha256"
        ] = corrected["materialization_receipt_sha256"]
        corrected["execution_partition_sha256"] = (
            controller.canonical_sha256(
                self.plan["execution_partition"]
            )
        )
        immutable["bundle_manifest"] = self.repin_amendment_record(
            immutable["bundle_manifest"], self.bundle_path
        )
        immutable["materialization_receipt"] = (
            self.repin_amendment_record(
                immutable["materialization_receipt"],
                self.materialization_path,
            )
        )
        capacity_binding = payload["capacity_to_partition_binding"]
        capacity_binding["corrected_execution_partition_sha256"] = (
            corrected["execution_partition_sha256"]
        )
        capacity_binding["corrected_partition_artifact_sha256"] = (
            corrected["partition"]["sha256"]
        )
        for capacity_row in capacity_binding["selected_rows"]:
            capacity_row["remote_wall_clock_seconds"] = int(
                capacity_row["remote_wall_clock_seconds"]
            )
        payload["amendment_semantic_sha256"] = (
            projector.amendment_tool.semantic_sha256(
                {
                    key: value
                    for key, value in payload.items()
                    if key != "amendment_semantic_sha256"
                }
            )
        )
        validated = projector.amendment_tool.validate_amendment_payload(
            payload
        )
        self.amendment_path.write_bytes(
            projector.amendment_tool.canonical_json_bytes(validated)
        )
        amendment_sha = sha(self.amendment_path)
        readback = {
            "schema": projector.amendment_tool.READBACK_SCHEMA,
            "status": "PASS",
            "authority_state": projector.amendment_tool.AUTHORITY_STATE,
            "amendment": str(self.amendment_path),
            "amendment_sha256": amendment_sha,
            "amendment_semantic_sha256": validated[
                "amendment_semantic_sha256"
            ],
            "all_evidence_rehashed": True,
            "byte_exact_rebuild": True,
            "submission_performed": False,
            "full_training_authority": 0,
            "full_extraction_authority": False,
        }
        self.readback_path.write_bytes(
            projector.amendment_tool.canonical_json_bytes(readback)
        )

    def rebuild_storage_evidence(self) -> None:
        self.write_valid_amendment_chain()
        budget_ref = {
            "path": str(self.budget_path),
            "sha256": sha(self.budget_path),
        }
        normalized_budget, blockers = projector.validate_controller_budget(
            budget_ref
        )
        assert normalized_budget is not None and not blockers
        derivation_ref = {
            "path": str(self.controller_derivation_path),
            "sha256": sha(self.controller_derivation_path),
        }
        bindings = {
            "plan": {
                "path": str(self.plan_path),
                "sha256": sha(self.plan_path),
                "size_bytes": self.plan_path.stat().st_size,
            },
            "preflight_receipt": {
                "path": str(self.receipt_path),
                "sha256": sha(self.receipt_path),
                "size_bytes": self.receipt_path.stat().st_size,
            },
            "capacity_count_amendment": {
                "path": str(self.amendment_path),
                "sha256": sha(self.amendment_path),
                "size_bytes": self.amendment_path.stat().st_size,
            },
            "capacity_count_amendment_readback": {
                "path": str(self.readback_path),
                "sha256": sha(self.readback_path),
                "size_bytes": self.readback_path.stat().st_size,
            },
            "bundle_manifest": {
                "path": str(self.bundle_path),
                "sha256": sha(self.bundle_path),
                "size_bytes": self.bundle_path.stat().st_size,
            },
            "materialization_receipt": {
                "path": str(self.materialization_path),
                "sha256": sha(self.materialization_path),
                "size_bytes": self.materialization_path.stat().st_size,
            },
            "execution_partition_sha256": controller.canonical_sha256(
                self.plan["execution_partition"]
            ),
            "partition_artifact_sha256": sha(self.partition_path),
        }
        rows = [
            {
                "row_id": row["row_id"],
                "system": row["system"],
                "source_tuple_count": row["input_contract"][
                    "source_tuple_count"
                ],
                "expected_job_count": row["input_contract"][
                    "expected_job_count"
                ],
            }
            for row in self.plan["rows"]
        ]
        amendment = json.loads(
            self.amendment_path.read_text(encoding="utf-8")
        )
        witnesses = projector.derive_capacity_witnesses(amendment)
        envelopes = []
        for row in rows:
            witness_id = (
                "pp_background_jet8"
                if row["system"] == "pp"
                else "auau_background_jet12"
            )
            classes = []
            for name, contract in projector.ARTIFACT_CLASS_CONTRACT.items():
                if contract["requires_witness"]:
                    witness = witnesses[witness_id][name]
                    ceiling = projector.ceil_ratio(
                        witness["size_bytes"], 5, 4
                    )
                    witness_ids = [witness_id]
                else:
                    ceiling = 64_000
                    witness_ids = []
                classes.append(
                    {
                        "artifact_class": name,
                        "storage_domain": contract["domain"],
                        "count_per_job": contract["count_per_job"],
                        "bytes_per_artifact_ceiling": ceiling,
                        "basis": contract["basis"],
                        "witness_row_ids": witness_ids,
                    }
                )
            envelopes.append(
                {"row_id": row["row_id"], "artifact_classes": classes}
            )
        spec = {
            "controller_dry_materialization_budget": budget_ref,
            "controller_dry_materialization_derivation": derivation_ref,
            "retry_reserve": {"numerator": 1, "denominator": 10},
            "row_envelopes": envelopes,
        }
        manifest = projector.build_storage_manifest_from_validated(
            spec,
            plan=self.plan,
            rows=rows,
            witnesses=witnesses,
            bindings=bindings,
        )
        self.assert_pass_manifest(manifest)
        write_json(self.storage_manifest_path, manifest)

        now = int(time.time())
        quota_stdout = (
            "Disk quotas for usr fixture:\n"
            "Filesystem kbytes quota limit grace files quota limit grace\n"
            f"{self.remote_root} 1K 1000T 1000T - "
            "10 10000000 10000000 -\n"
        )
        def fake_lfs(argv, **_kwargs):
            return subprocess.CompletedProcess(argv, 0, quota_stdout, "")

        snapshot = projector.build_quota_snapshot(
            self.plan_path,
            "fixture",
            max_age_seconds=900,
            runner=fake_lfs,
            now_seconds=now,
            lfs_binary="/usr/bin/lfs",
        )
        assert snapshot["status"] == projector.PASS_STATUS
        write_json(self.quota_snapshot_path, snapshot)
        certificate = projector.build_projection(
            self.storage_manifest_path,
            self.quota_snapshot_path,
            now_seconds=now,
        )
        write_json(self.storage_certificate_path, certificate)

    @staticmethod
    def assert_pass_manifest(manifest: dict) -> None:
        assert manifest["status"] == projector.PASS_STATUS
        assert manifest["blocker_codes"] == []

    def restore(self) -> None:
        self.materialization_path.chmod(0o644)
        write_json(self.materialization_path, self.baseline_materialization)
        self.materialization_path.chmod(0o444)
        write_json(self.source_path, self.baseline_sources)
        self.plan = copy.deepcopy(self.baseline_plan)
        self.write_plan_and_receipt()
        write_json(self.budget_path, self.baseline_budget)
        write_json(
            self.controller_derivation_path,
            self.baseline_controller_derivation,
        )
        self.rebuild_storage_evidence()

    def args(
        self,
        staging_root: Path,
        *,
        action: str = "preflight",
        profile: str = controller.ARTIFACT_PROFILE,
    ) -> argparse.Namespace:
        return argparse.Namespace(
            action=action,
            plan=self.plan_path,
            plan_sha256=sha(self.plan_path),
            preflight_receipt=self.receipt_path,
            preflight_receipt_sha256=sha(self.receipt_path),
            storage_certificate=self.storage_certificate_path,
            storage_certificate_sha256=sha(self.storage_certificate_path),
            controller_budget=self.budget_path,
            controller_budget_sha256=sha(self.budget_path),
            artifact_profile=profile,
            staging_root=staging_root,
            overwrite=False,
        )

    def command(
        self,
        staging_root: Path,
        *,
        action: str = "preflight",
        profile: str = controller.ARTIFACT_PROFILE,
        overwrite: bool = False,
    ) -> list[str]:
        command = [
            sys.executable,
            str(CONTROLLER_PATH),
            action,
            "--plan",
            str(self.plan_path),
            "--plan-sha256",
            sha(self.plan_path),
            "--preflight-receipt",
            str(self.receipt_path),
            "--preflight-receipt-sha256",
            sha(self.receipt_path),
            "--storage-certificate",
            str(self.storage_certificate_path),
            "--storage-certificate-sha256",
            sha(self.storage_certificate_path),
            "--controller-budget",
            str(self.budget_path),
            "--controller-budget-sha256",
            sha(self.budget_path),
            "--artifact-profile",
            profile,
            "--staging-root",
            str(staging_root),
        ]
        if overwrite:
            command.append("--overwrite")
        return command


class FullControllerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.temporary = tempfile.TemporaryDirectory()
        cls.root = Path(cls.temporary.name).resolve()
        cls.fixture = FullControllerFixture(cls.root / "fixture")

    @classmethod
    def tearDownClass(cls) -> None:
        for path in sorted(
            cls.root.rglob("*"),
            key=lambda value: len(value.parts),
            reverse=True,
        ):
            try:
                if path.is_dir():
                    path.chmod(0o755)
                elif path.exists() and not path.is_symlink():
                    path.chmod(0o644)
            except FileNotFoundError:
                pass
        cls.temporary.cleanup()

    def setUp(self) -> None:
        self.fixture.restore()
        self.staging_root = self.root / f"stage_{self._testMethodName}"
        if self.staging_root.exists():
            shutil = __import__("shutil")
            shutil.rmtree(self.staging_root)

    def validate(self) -> dict:
        return controller.validate_plan_and_evidence(
            self.fixture.args(self.staging_root)
        )

    def repin_mutated_plan(self, *, rebuild_storage: bool = False) -> None:
        self.fixture.write_plan_and_receipt()
        if rebuild_storage:
            self.fixture.rebuild_storage_evidence()

    def test_preflight_and_local_materialization_are_exact_and_inert(self) -> None:
        preflight = subprocess.run(
            self.fixture.command(self.staging_root, action="preflight"),
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        self.assertEqual(preflight.returncode, 0, msg=preflight.stderr)
        preflight_payload = json.loads(preflight.stdout)
        self.assertEqual(preflight_payload["status"], "PASS_PREFLIGHT_NO_WRITE")
        self.assertFalse(preflight_payload["staging_root_created"])
        self.assertFalse(self.staging_root.exists())

        result = subprocess.run(
            self.fixture.command(self.staging_root, action="materialize"),
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        self.assertEqual(result.returncode, 0, msg=result.stderr)
        payload = json.loads(result.stdout)
        self.assertEqual(payload["row_count"], 13)
        self.assertEqual(payload["source_tuple_occurrence_count"], 129_998)
        self.assertEqual(payload["job_count"], 18_577)
        self.assertEqual(payload["physical_root_artifact_count"], 37_154)
        self.assertFalse(payload["submission_performed"])
        self.assertEqual(
            sorted(path.name for path in self.staging_root.iterdir()),
            sorted(controller.OUTPUT_FILENAMES),
        )
        manifest = json.loads(
            (self.staging_root / "materialization_manifest.json").read_text()
        )
        self.assertEqual(
            manifest["artifact_profile"], controller.ARTIFACT_PROFILE
        )
        self.assertFalse(manifest["authority"]["submission_authority"])
        self.assertFalse(manifest["authority"]["physics_output_authority"])
        self.assertFalse(manifest["authority"]["replay_cache_authority"])
        self.assertFalse(Path(self.fixture.campaign["output_root"]).exists())
        self.assertFalse(Path(self.fixture.campaign["submit_root"]).exists())

    def test_one_byte_controller_budget_boundary(self) -> None:
        context = self.validate()
        artifacts = controller.build_staged_artifacts(context)
        required = sum(len(data) for data in artifacts.values())
        exact_budget = {
            "projected_bytes": required,
            "projected_inodes": len(artifacts),
        }
        self.assertEqual(
            controller.enforce_controller_budget(
                artifacts, exact_budget
            )[0],
            required,
        )

        with self.assertRaisesRegex(
            controller.ControllerError, "exceeds certified bytes"
        ):
            controller.enforce_controller_budget(
                artifacts,
                {
                    "projected_bytes": required - 1,
                    "projected_inodes": len(artifacts),
                },
            )

    def test_existing_output_or_submit_namespace_is_rejected(self) -> None:
        for field in ("output_root", "submit_root"):
            with self.subTest(field=field):
                namespace = Path(self.fixture.campaign[field])
                namespace.mkdir(parents=True)
                try:
                    with self.assertRaisesRegex(
                        controller.ControllerError,
                        "namespace already exists",
                    ):
                        self.validate()
                finally:
                    namespace.rmdir()
                    namespace.parent.rmdir()

    def test_missing_controller_budget_certificate_is_rejected(self) -> None:
        args = self.fixture.args(self.staging_root)
        args.controller_budget = self.root / "missing_controller_budget.json"
        args.controller_budget_sha256 = "0" * 64
        with self.assertRaisesRegex(
            controller.ControllerError, "controller storage budget is missing"
        ):
            controller.validate_plan_and_evidence(args)

    def test_duplicate_logical_row_is_rejected(self) -> None:
        self.fixture.plan["rows"][1] = copy.deepcopy(
            self.fixture.plan["rows"][0]
        )
        self.repin_mutated_plan()
        with self.assertRaisesRegex(
            controller.ControllerError, "duplicated"
        ):
            self.validate()

    def test_wrong_global_or_row_count_is_rejected(self) -> None:
        self.fixture.plan["execution_partition"][
            "source_tuple_count"
        ] = 129_997
        self.repin_mutated_plan()
        with self.assertRaisesRegex(
            controller.ControllerError, "partition differs"
        ):
            self.validate()

        self.fixture.restore()
        self.fixture.plan["rows"][0]["input_contract"][
            "source_tuple_count"
        ] -= 1
        self.repin_mutated_plan()
        with self.assertRaisesRegex(
            controller.ControllerError, "tuple count differs"
        ):
            self.validate()

    def test_wrong_training_or_submission_authority_is_rejected(self) -> None:
        self.fixture.plan["authority"]["full_training_authority"] = 1
        self.repin_mutated_plan()
        with self.assertRaisesRegex(
            controller.ControllerError, "full_training_authority"
        ):
            self.validate()

    def test_source_authority_builder_contract_drift_is_rejected(self) -> None:
        sources = json.loads(self.fixture.source_path.read_text())
        sources["authority"]["source_ownership_state"] = "ambiguous"
        self.fixture.source_path.write_text(
            json.dumps(sources, indent=2, sort_keys=True) + "\n"
        )
        self.fixture.plan["input_manifests"]["sources"]["sha256"] = sha(
            self.fixture.source_path
        )
        self.fixture.plan["execution_fingerprint_sha256"] = (
            self.fixture.execution_fingerprint(self.fixture.plan)
        )
        self.repin_mutated_plan()
        with self.assertRaisesRegex(
            controller.ControllerError,
            "source-authority V2 field inventory differs",
        ):
            self.validate()

    def test_replay_cache_and_the121_the122_authority_are_rejected(self) -> None:
        self.fixture.plan["closure_witness_boundary"][
            "replay_cache_authority"
        ] = True
        self.repin_mutated_plan()
        with self.assertRaisesRegex(
            controller.ControllerError, "forbidden authority"
        ):
            self.validate()

        self.fixture.restore()
        self.fixture.plan["closure_witness_boundary"][
            "the121_authority"
        ] = True
        self.repin_mutated_plan()
        with self.assertRaisesRegex(
            controller.ControllerError, "forbidden authority"
        ):
            self.validate()

    def test_profile_hash_bundle_and_materialization_mismatch_fail_closed(self) -> None:
        with self.assertRaisesRegex(
            controller.ControllerError, "artifact profile"
        ):
            controller.validate_plan_and_evidence(
                self.fixture.args(
                    self.staging_root, profile="THE121_REPLAY_CACHE_V1"
                )
            )

        self.fixture.plan["rows"][0]["bundle_contract"][
            "code_sha256"
        ] = "9" * 64
        row = self.fixture.plan["rows"][0]
        del row["row_fingerprint_sha256"]
        row["row_fingerprint_sha256"] = controller.canonical_sha256(row)
        self.fixture.plan["execution_fingerprint_sha256"] = (
            self.fixture.execution_fingerprint(self.fixture.plan)
        )
        self.repin_mutated_plan()
        with self.assertRaisesRegex(controller.ControllerError, "bundle"):
            self.validate()

        self.fixture.restore()
        materialization = json.loads(
            self.fixture.materialization_path.read_text()
        )
        materialization["bundle_identity_sha256"] = "8" * 64
        self.fixture.materialization_path.chmod(0o644)
        write_json(self.fixture.materialization_path, materialization)
        self.fixture.materialization_path.chmod(0o444)
        self.fixture.plan["input_manifests"]["materialization"][
            "sha256"
        ] = sha(self.fixture.materialization_path)
        self.fixture.plan["execution_fingerprint_sha256"] = (
            self.fixture.execution_fingerprint(self.fixture.plan)
        )
        self.repin_mutated_plan()
        with self.assertRaisesRegex(
            controller.ControllerError, "materialization readback failed"
        ):
            self.validate()

    def test_symlink_to_identical_materialized_bytes_is_rejected(self) -> None:
        artifact = self.fixture.by_role["pp_config"]
        bundle_path = Path(artifact["path"])
        original_bytes = bundle_path.read_bytes()
        original_mode = bundle_path.stat().st_mode & 0o777
        replacement = self.root / "identical_pp_config.artifact"
        replacement.write_bytes(original_bytes)
        replacement.chmod(original_mode)
        role_parent = bundle_path.parent
        role_parent.chmod(0o755)
        bundle_path.unlink()
        bundle_path.symlink_to(replacement)
        role_parent.chmod(0o555)
        try:
            with self.assertRaisesRegex(
                controller.ControllerError,
                "mutable artifact path drift|symlink-free regular file",
            ):
                self.validate()
        finally:
            role_parent.chmod(0o755)
            bundle_path.unlink()
            bundle_path.write_bytes(original_bytes)
            bundle_path.chmod(original_mode)
            role_parent.chmod(0o555)
            replacement.unlink()

    def test_no_submission_primitive_and_overwrite_is_forbidden(self) -> None:
        tree = ast.parse(CONTROLLER_PATH.read_text(encoding="utf-8"))
        imported_roots = {
            alias.name.split(".")[0]
            for node in ast.walk(tree)
            if isinstance(node, (ast.Import, ast.ImportFrom))
            for alias in node.names
        }
        self.assertNotIn("subprocess", imported_roots)
        self.assertNotIn("socket", imported_roots)
        self.assertNotIn("urllib", imported_roots)
        calls = [
            node
            for node in ast.walk(tree)
            if isinstance(node, ast.Call)
            and isinstance(node.func, ast.Attribute)
            and node.func.attr in {"system", "popen", "spawn", "execv"}
        ]
        self.assertEqual(calls, [])

        result = subprocess.run(
            self.fixture.command(
                self.staging_root, action="materialize", overwrite=True
            ),
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        self.assertEqual(result.returncode, 2)
        self.assertIn("--overwrite is explicitly forbidden", result.stderr)
        self.assertFalse(self.staging_root.exists())

    def test_staging_dotdot_symlink_parent_and_campaign_overlap_are_rejected(
        self,
    ) -> None:
        dotdot = self.root / "fixture" / ".." / "unsafe_stage"
        with self.assertRaisesRegex(controller.ControllerError, "unsafe"):
            controller.validate_plan_and_evidence(
                self.fixture.args(dotdot)
            )

        real_parent = self.root / "real_stage_parent"
        real_parent.mkdir()
        symlink_parent = self.root / "stage_parent_alias"
        symlink_parent.symlink_to(real_parent, target_is_directory=True)
        try:
            with self.assertRaisesRegex(
                controller.ControllerError, "symlink"
            ):
                controller.validate_plan_and_evidence(
                    self.fixture.args(symlink_parent / "stage")
                )
        finally:
            symlink_parent.unlink()
            real_parent.rmdir()

        overlapping = Path(self.fixture.campaign["evidence_root"]) / "stage"
        overlapping.parent.mkdir(parents=True)
        try:
            with self.assertRaisesRegex(
                controller.ControllerError, "overlaps declared campaign"
            ):
                controller.validate_plan_and_evidence(
                    self.fixture.args(overlapping)
                )
        finally:
            overlapping.parent.rmdir()
            overlapping.parent.parent.rmdir()

    def test_exact_sdcc_scratch_alias_resolves_to_canonical_gpfs(self) -> None:
        alias_root = self.root / "sphenix_u"
        canonical_root = self.root / "gpfs_user"
        (alias_root / "alice").mkdir(parents=True)
        (canonical_root / "alice" / "evidence").mkdir(parents=True)
        (alias_root / "alice" / "scratch").symlink_to(
            canonical_root / "alice",
            target_is_directory=True,
        )
        requested = alias_root / "alice" / "scratch" / "evidence" / "stage"
        with (
            mock.patch.object(controller, "SDCC_USER_ALIAS_ROOT", alias_root),
            mock.patch.object(
                controller,
                "SDCC_CANONICAL_USER_ROOT",
                canonical_root,
            ),
        ):
            self.assertEqual(
                controller.validate_staging_root(requested),
                canonical_root / "alice" / "evidence" / "stage",
            )

    def test_sdcc_scratch_alias_rejects_wrong_user_target(self) -> None:
        alias_root = self.root / "wrong_user_alias"
        canonical_root = self.root / "wrong_user_gpfs"
        (alias_root / "alice").mkdir(parents=True)
        (canonical_root / "alice").mkdir(parents=True)
        (canonical_root / "bob").mkdir(parents=True)
        (alias_root / "alice" / "scratch").symlink_to(
            canonical_root / "bob",
            target_is_directory=True,
        )
        with (
            mock.patch.object(controller, "SDCC_USER_ALIAS_ROOT", alias_root),
            mock.patch.object(
                controller,
                "SDCC_CANONICAL_USER_ROOT",
                canonical_root,
            ),
            self.assertRaisesRegex(
                controller.ControllerError,
                "does not target the exact user GPFS root",
            ),
        ):
            controller.validate_staging_root(
                alias_root / "alice" / "scratch" / "stage"
            )

    def test_sdcc_scratch_alias_rejects_nested_escape(self) -> None:
        alias_root = self.root / "nested_alias"
        canonical_root = self.root / "nested_gpfs_user"
        external_root = self.root / "nested_external"
        (alias_root / "alice").mkdir(parents=True)
        (canonical_root / "alice").mkdir(parents=True)
        external_root.mkdir()
        (alias_root / "alice" / "scratch").symlink_to(
            canonical_root / "alice",
            target_is_directory=True,
        )
        (canonical_root / "alice" / "evidence").symlink_to(
            external_root,
            target_is_directory=True,
        )
        with (
            mock.patch.object(controller, "SDCC_USER_ALIAS_ROOT", alias_root),
            mock.patch.object(
                controller,
                "SDCC_CANONICAL_USER_ROOT",
                canonical_root,
            ),
            self.assertRaisesRegex(
                controller.ControllerError,
                "resolution is unstable or escapes GPFS",
            ),
        ):
            controller.validate_staging_root(
                alias_root / "alice" / "scratch" / "evidence" / "stage"
            )

    def test_cli_help_lists_only_non_submitting_actions(self) -> None:
        result = subprocess.run(
            [sys.executable, str(CONTROLLER_PATH), "--help"],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        self.assertEqual(result.returncode, 0, msg=result.stderr)
        self.assertIn("{preflight,materialize}", result.stdout)
        self.assertNotIn("{preflight,materialize,submit}", result.stdout)


if __name__ == "__main__":
    unittest.main()
