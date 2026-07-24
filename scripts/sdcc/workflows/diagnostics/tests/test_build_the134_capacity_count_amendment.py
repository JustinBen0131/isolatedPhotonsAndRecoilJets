#!/usr/bin/env python3
"""Adversarial tests for the THE-134 capacity-count amendment."""

from __future__ import annotations

import copy
import contextlib
import hashlib
import importlib.util
import io
import json
import os
import stat
import tempfile
import unittest
from pathlib import Path
from unittest import mock


HERE = Path(__file__).resolve().parent
MODULE_PATH = HERE.parent / "build_the134_capacity_count_amendment.py"
SPEC = importlib.util.spec_from_file_location(
    "the134_capacity_count_amendment", MODULE_PATH
)
assert SPEC is not None and SPEC.loader is not None
amendment = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(amendment)


class CapacityCountAmendmentTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name).resolve()

    def tearDown(self) -> None:
        self.temporary.cleanup()

    @staticmethod
    def file_sha256(path: Path) -> str:
        return hashlib.sha256(path.read_bytes()).hexdigest()

    @staticmethod
    def artifact(path: Path) -> dict[str, str]:
        return {
            "path": str(path),
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        }

    def test_submitted_args_hash_reconstructs_condor_classad_bytes(
        self,
    ) -> None:
        args_path = self.root / "row.args"
        args_path.write_text(
            "sample /tmp/chunk.list dataset $(Cluster) 0 1 NONE "
            "/tmp/output_$(Process)_$(ClusterId)_$(ProcId).root\n",
            encoding="utf-8",
        )
        expected = hashlib.sha256(
            b"0 sample /tmp/chunk.list dataset 5616584 0 1 NONE "
            b"/tmp/output_0_5616584_0.root\n"
        ).hexdigest()
        self.assertEqual(
            amendment.submitted_args_sha256_from_template(
                "pp_background_jet8",
                args_path,
                "5616584.0",
            ),
            expected,
        )
        self.assertNotEqual(self.file_sha256(args_path), expected)

    def test_submitted_args_hash_rejects_malformed_identity_and_rows(
        self,
    ) -> None:
        args_path = self.root / "row.args"
        args_path.write_text(
            "sample /tmp/chunk.list dataset $(Cluster) 0 1 NONE "
            "/tmp/output_$(Process).root\n",
            encoding="utf-8",
        )
        for cluster_proc in (
            "0.0",
            "5616584",
            "5616584.-1",
            "5616584.1",
            "5616584.0.0",
            " 5616584.0",
        ):
            with self.subTest(cluster_proc=cluster_proc):
                with self.assertRaises(amendment.AmendmentError):
                    amendment.submitted_args_sha256_from_template(
                        "pp_background_jet8",
                        args_path,
                        cluster_proc,
                    )

        for contents in (
            "",
            "first\nsecond\n",
            "first\n\n",
        ):
            with self.subTest(contents=contents):
                args_path.write_text(contents, encoding="utf-8")
                with self.assertRaisesRegex(
                    amendment.AmendmentError,
                    "exactly one argument row",
                ):
                    amendment.submitted_args_sha256_from_template(
                        "pp_background_jet8",
                        args_path,
                        "5616584.0",
                    )

        for contents in (
            "sample chunk dataset 5616584 0 1 NONE output\n",
            "sample chunk dataset $(Cluster) 1 1 NONE output\n",
            "sample chunk dataset $(Cluster) 0 2 NONE output\n",
            "sample chunk dataset $(Cluster) 0 1 NOT_NONE output\n",
            "sample chunk dataset $(Cluster) 0 1 NONE\n",
            "sample chunk dataset $(Cluster) 0 1 NONE output extra\n",
        ):
            with self.subTest(contract_contents=contents):
                args_path.write_text(contents, encoding="utf-8")
                with self.assertRaisesRegex(
                    amendment.AmendmentError,
                    "one-proc group-7 capacity contract",
                ):
                    amendment.submitted_args_sha256_from_template(
                        "pp_background_jet8",
                        args_path,
                        "5616584.0",
                    )

    def test_submitted_args_hash_changes_on_template_or_cluster_drift(
        self,
    ) -> None:
        args_path = self.root / "row.args"
        args_path.write_text(
            "sample /tmp/chunk.list dataset $(Cluster) 0 1 NONE "
            "/tmp/output_$(Process).root\n",
            encoding="utf-8",
        )
        baseline = amendment.submitted_args_sha256_from_template(
            "pp_background_jet8",
            args_path,
            "5616584.0",
        )
        changed_cluster = amendment.submitted_args_sha256_from_template(
            "pp_background_jet8",
            args_path,
            "5616585.0",
        )
        args_path.write_text(
            "sample /tmp/other.list dataset $(Cluster) 0 1 NONE "
            "/tmp/output_$(Process).root\n",
            encoding="utf-8",
        )
        changed_template = amendment.submitted_args_sha256_from_template(
            "pp_background_jet8",
            args_path,
            "5616584.0",
        )
        self.assertNotEqual(baseline, changed_cluster)
        self.assertNotEqual(baseline, changed_template)

    def test_immutable_authority_preserves_lexical_materialization_alias(
        self,
    ) -> None:
        physical_root = self.root / "physical"
        physical_root.mkdir()
        stable_alias = self.root / "stable_alias"
        stable_alias.symlink_to(physical_root, target_is_directory=True)

        bundle_path = physical_root / "resolver_bundle_receipt.json"
        bundle_path.write_text("{}\n", encoding="utf-8")
        materialization_path = (
            physical_root
            / f"the134_bundle_sha256_{'a' * 64}"
            / "metadata"
            / "materialization_receipt.json"
        )
        materialization_path.parent.mkdir(parents=True)
        materialization_path.write_text("{}\n", encoding="utf-8")

        bundle_alias = stable_alias / bundle_path.name
        materialization_alias = (
            stable_alias
            / materialization_path.relative_to(physical_root)
        )
        bundle_artifact = {
            **self.artifact(bundle_alias),
            "resolved_path": str(bundle_path.resolve(strict=True)),
            "size_bytes": bundle_path.stat().st_size,
        }
        materialization_artifact = {
            **self.artifact(materialization_alias),
            "resolved_path": str(materialization_path.resolve(strict=True)),
            "size_bytes": materialization_path.stat().st_size,
        }
        bundle = {
            "artifact_by_role": {},
            "bundle_identity_sha256": "a" * 64,
            "semantic_fingerprint_sha256": "b" * 64,
            "public_commit": "c" * 40,
            "code_sha256": "d" * 64,
            "replay_schema_sha256": "e" * 64,
            "training_schema_sha256": "f" * 64,
            "semantic_sha256": "1" * 64,
            "runtime": {},
        }
        materialization = {
            "digest_named_bundle_path": str(
                materialization_alias.parent.parent
            ),
            "artifact_count": 0,
            "total_artifact_bytes": 0,
            "readback": "PASS_SYMLINK_FREE_READONLY_CONTENT_EXACT",
        }
        with (
            mock.patch.object(
                amendment,
                "load_json_artifact",
                side_effect=[
                    ({}, bundle_artifact),
                    ({}, materialization_artifact),
                ],
            ),
            mock.patch.object(
                amendment.resolver,
                "validate_bundle",
                return_value=bundle,
            ),
            mock.patch.object(
                amendment.resolver,
                "validate_materialization_binding",
                return_value=materialization,
            ) as binding,
        ):
            observed = amendment.validate_immutable_authority(
                {
                    "bundle_manifest": {},
                    "materialization_receipt": {},
                }
            )

        self.assertEqual(observed["status"], "PASS")
        self.assertEqual(
            binding.call_args.kwargs["materialization_path"],
            materialization_alias,
        )
        self.assertEqual(
            binding.call_args.kwargs["bundle_path"],
            bundle_path.resolve(strict=True),
        )

    def make_source_manifests(
        self,
    ) -> tuple[dict[str, object], dict[str, object]]:
        rows_v2 = []
        for index, row_id in enumerate(amendment.expected_row_ids()):
            tuple_count = (
                9_998
                if row_id == "auau_background_jet12"
                else 10_000
            )
            input_lists = []
            for role in amendment.resolver.LIST_ROLES:
                input_lists.append(
                    {
                        "role": role,
                        "path": f"/source/{row_id}/{role}.list",
                        "resolved_path": f"/resolved/{row_id}/{role}.list",
                        "sha256": f"{index + 1:064x}"[-64:],
                        "size_bytes": 100 + index,
                        "line_count": tuple_count,
                        "executable_count": tuple_count,
                    }
                )
            rows_v2.append(
                {
                    "row_id": row_id,
                    "tuple_count": tuple_count,
                    "tuple_records_sha256": f"{100 + index:064x}",
                    "full_source_manifest_sha256": f"{200 + index:064x}",
                    "first_tuple_sha256": f"{300 + index:064x}",
                    "last_tuple_sha256": f"{400 + index:064x}",
                    "input_lists": input_lists,
                }
            )
        global_check = {
            "total": amendment.EXPECTED_SOURCE_TUPLES,
            "unique": amendment.EXPECTED_SOURCE_TUPLES,
            "duplicate_count": 0,
            "tuple_identity_records_sha256": "a" * 64,
        }
        common = {
            "status": "PASS",
            "authority_state": "MEASURED_CANDIDATE_NOT_SCIENTIFIC_AUTHORITY",
            "authority": {"pp_period": "0mrad"},
            "sim_list_root": "/source",
            "resolved_sim_list_root": "/resolved",
            "row_count": 13,
            "global_five_tuple_check": global_check,
            "source_rows_semantic_sha256": "b" * 64,
            "manifest_semantic_sha256": "c" * 64,
        }
        v2 = {
            "schema": amendment.source_builder.SOURCE_SCHEMA,
            **copy.deepcopy(common),
            "rows": rows_v2,
        }
        rows_v1 = []
        for row in rows_v2:
            rows_v1.append(
                {
                    **copy.deepcopy(row),
                    "expected_input_count": row["tuple_count"],
                    "expected_occurrence_count": row["tuple_count"],
                }
            )
        v1 = {
            "schema": amendment.source_builder.SOURCE_SCHEMA_V1,
            **copy.deepcopy(common),
            "rows": rows_v1,
        }
        return v1, v2

    def make_valid_amendment_payload(self) -> dict[str, object]:
        def pinned(
            name: str,
            *,
            schema: str | None = None,
            extra: dict[str, object] | None = None,
        ) -> dict[str, object]:
            record: dict[str, object] = {
                "label": name,
                "path": str(self.root / name),
                "resolved_path": str(self.root / name),
                "sha256": hashlib.sha256(name.encode("utf-8")).hexdigest(),
                "size_bytes": 100,
            }
            if schema is not None:
                record["schema"] = schema
            if extra is not None:
                record.update(extra)
            return record

        inventory = amendment.resolver.inventory_rows()
        row_ids = amendment.expected_row_ids()
        legacy_rows = []
        corrected_rows = []
        corrected_semantics = []
        crosswalk_rows = []
        for index, descriptor in enumerate(inventory):
            tuple_count = (
                9_998
                if descriptor["row_id"] == "auau_background_jet12"
                else 10_000
            )
            chunk_count = 1_429
            tail_count = 2 if tuple_count == 9_998 else 4
            legacy_rows.append(
                {
                    "row_id": descriptor["row_id"],
                    "source_tuple_count": tuple_count,
                    "expected_chunk_count": chunk_count,
                    "legacy_expected_occurrence_count": tuple_count,
                }
            )
            corrected_rows.append(
                {
                    "row_id": descriptor["row_id"],
                    "source_tuple_count": tuple_count,
                    "expected_chunk_count": chunk_count,
                    "expected_job_count": chunk_count,
                    "expected_output_pair_count": chunk_count,
                    "expected_source_occurrence_count": chunk_count,
                    "tail_tuple_count": tail_count,
                    "row_partition_sha256": f"{1000 + index:064x}",
                    "chunk_records_sha256": f"{1100 + index:064x}",
                    "first_execution_chunk_sha256": f"{1200 + index:064x}",
                    "row_fingerprint_sha256": f"{1300 + index:064x}",
                }
            )
            system = descriptor["system"]
            sample = descriptor["sample"]
            corrected_semantics.append(
                {
                    "row_id": descriptor["row_id"],
                    "system": system,
                    "lane": descriptor["lane"],
                    "dataset": descriptor["dataset"],
                    "sample": sample,
                    "source_role": descriptor["source_role"],
                    "minimum_bias_gate": descriptor[
                        "minimum_bias_gate"
                    ],
                    "period": "0mrad" if system == "pp" else "AUAU_RUN24",
                    "run": 28,
                    "si_di_role": "SI" if system == "pp" else "EMBEDDED",
                    "ownership_state": "source_role_frozen",
                    "code_sha256": "a" * 64,
                    "replay_schema_sha256": "b" * 64,
                    "training_schema_sha256": "c" * 64,
                    "semantic_sha256": "d" * 64,
                    "library_sha256": hashlib.sha256(
                        f"{system}_library".encode("utf-8")
                    ).hexdigest(),
                    "model_sha256": hashlib.sha256(
                        f"{system}_model".encode("utf-8")
                    ).hexdigest(),
                    "base_config_sha256": "b" * 64,
                }
            )
            crosswalk_rows.append(
                {
                    "row_id": descriptor["row_id"],
                    "tuple_count": tuple_count,
                    "tuple_records_sha256": f"{2000 + index:064x}",
                    "full_source_manifest_sha256": f"{2100 + index:064x}",
                    "first_tuple_sha256": f"{2200 + index:064x}",
                    "last_tuple_sha256": f"{2300 + index:064x}",
                    "input_lists": [
                        {
                            "role": role,
                            "path": str(
                                self.root
                                / descriptor["row_id"]
                                / f"{role}.list"
                            ),
                            "resolved_path": str(
                                self.root
                                / descriptor["row_id"]
                                / f"{role}.list"
                            ),
                            "sha256": f"{2400 + index:064x}",
                            "size_bytes": 1,
                            "line_count": tuple_count,
                            "executable_count": tuple_count,
                        }
                        for role in amendment.resolver.LIST_ROLES
                    ],
                    "legacy_tuple_aliases_preserved_as_history_only": True,
                    "corrected_execution_aliases_absent": True,
                }
            )

        immutable_bundle = pinned(
            "bundle.json", schema=amendment.resolver.BUNDLE_SCHEMA
        )
        immutable_materialization = pinned(
            "materialization.json",
            schema=amendment.resolver.MATERIALIZATION_SCHEMA,
        )
        immutable_artifacts = {
            role: {
                "role": role,
                "path": str(self.root / "bundle" / role / role),
                "resolved_path": str(self.root / "bundle" / role / role),
                "sha256": hashlib.sha256(role.encode("utf-8")).hexdigest(),
                "size_bytes": 1,
            }
            for role in amendment.resolver.REQUIRED_ARTIFACT_ROLES
        }
        legacy_preflight = {
            "status": "PASS",
            "plan": pinned(
                "legacy-plan.json", schema=amendment.LEGACY_PLAN_SCHEMA
            ),
            "preflight_receipt": pinned(
                "legacy-receipt.json",
                schema=amendment.LEGACY_RECEIPT_SCHEMA,
            ),
            "source_manifest": pinned(
                "legacy-source.json",
                schema=amendment.source_builder.SOURCE_SCHEMA_V1,
            ),
            "submission_performed": False,
            "full_training_authority": 0,
            "execution_partition_sha256": "1" * 64,
            "bundle_manifest_sha256": immutable_bundle["sha256"],
            "materialization_receipt_sha256": immutable_materialization[
                "sha256"
            ],
            "source_manifest_sha256": "2" * 64,
            "source_tuple_count": amendment.EXPECTED_SOURCE_TUPLES,
            "expected_job_count": amendment.EXPECTED_CHUNKS,
            "expected_output_count": amendment.EXPECTED_CHUNKS,
            "legacy_expected_occurrence_count": (
                amendment.EXPECTED_SOURCE_TUPLES
            ),
            "rows": legacy_rows,
        }
        corrected_preflight = {
            "status": "PASS",
            "plan": pinned(
                "corrected-plan.json", schema=amendment.CORRECTED_PLAN_SCHEMA
            ),
            "preflight_receipt": pinned(
                "corrected-receipt.json",
                schema=amendment.CORRECTED_RECEIPT_SCHEMA,
            ),
            "source_manifest": pinned(
                "corrected-source.json",
                schema=amendment.CORRECTED_SOURCE_SCHEMA,
            ),
            "rows": pinned("corrected-rows.jsonl"),
            "duplicate_fingerprint": pinned("duplicate.sha256"),
            "partition": pinned(
                "partition.jsonl",
                extra={"record_count": amendment.EXPECTED_CHUNKS},
            ),
            "submission_performed": False,
            "full_training_authority": 0,
            "bundle_manifest_sha256": immutable_bundle["sha256"],
            "materialization_receipt_sha256": immutable_materialization[
                "sha256"
            ],
            "source_manifest_sha256": "3" * 64,
            "duplicate_fingerprint_sha256": "4" * 64,
            "execution_partition_sha256": "5" * 64,
            "execution_fingerprint_sha256": "6" * 64,
            "source_tuple_count": amendment.EXPECTED_SOURCE_TUPLES,
            "expected_chunk_count": amendment.EXPECTED_CHUNKS,
            "expected_job_count": amendment.EXPECTED_CHUNKS,
            "expected_output_pair_count": amendment.EXPECTED_CHUNKS,
            "expected_analysis_output_count": amendment.EXPECTED_CHUNKS,
            "expected_sidecar_output_count": amendment.EXPECTED_CHUNKS,
            "expected_physical_root_artifact_count": (
                amendment.EXPECTED_ROOT_ARTIFACTS
            ),
            "expected_source_occurrence_count": amendment.EXPECTED_CHUNKS,
            "source_occurrences_per_output_pair": 1,
            "rows_proof": corrected_rows,
            "capacity_source_semantics": corrected_semantics,
            "bundle_contract_semantic_sha256_by_system": {
                "pp": "7" * 64,
                "auau": "8" * 64,
            },
            "science_contract_semantic_sha256": "9" * 64,
        }
        source_crosswalk = {
            "status": "PASS",
            "row_count": 13,
            "source_tuple_count": amendment.EXPECTED_SOURCE_TUPLES,
            "global_unique_tuple_count": amendment.EXPECTED_SOURCE_TUPLES,
            "source_rows_semantic_sha256": "a" * 64,
            "manifest_semantic_sha256": "b" * 64,
            "global_tuple_identity_records_sha256": "c" * 64,
            "rows": crosswalk_rows,
        }
        capacity_evidence = {
            "resource_certificate": pinned(
                "resource.json", schema=amendment.CAPACITY_SCHEMA
            ),
            "root_join_certificate": pinned(
                "root-join.json", schema=amendment.ROOT_JOIN_SCHEMA
            ),
            "pp_audit": pinned(
                "pp-audit.json", schema=amendment.CAPACITY_AUDIT_SCHEMA
            ),
            "auau_audit": pinned(
                "auau-audit.json", schema=amendment.CAPACITY_AUDIT_SCHEMA
            ),
            "validation_authority": {
                "controller": pinned(
                    "controller.sh",
                    extra={"validation_commit": "d" * 40},
                ),
                "runtime_authority": pinned("runtime.json"),
                "submission_journal": pinned("journal.tsv"),
                "submission_manifest": pinned("manifest.tsv"),
                "submission_receipt": pinned("receipt.tsv"),
            },
        }
        capacity_rows = []
        for row_id in amendment.SELECTED_CAPACITY_ROWS:
            semantics = next(
                row
                for row in corrected_semantics
                if row["row_id"] == row_id
            )
            staged_sha = (
                "d" * 64 if row_id.startswith("pp_") else "e" * 64
            )
            source_manifest_sha = (
                "f" * 64 if row_id.startswith("pp_") else "0" * 63 + "1"
            )
            source_contract = {
                "dataset": semantics["dataset"],
                "input_file_sha256": staged_sha,
                "input_uri_hash": staged_sha,
                "lane": semantics["lane"],
                "ownership_state": "source_role_frozen",
                "period": semantics["period"],
                "run": 28,
                "sample": semantics["sample"],
                "segment": 1,
                "si_di_role": semantics["si_di_role"],
                "source_manifest_sha256": source_manifest_sha,
            }
            source_identity = amendment.canonical_source_identity_sha256(
                source_contract
            )
            source_occurrence = source_identity[:32]
            if source_occurrence == "0" * 32:
                source_occurrence = "0" * 31 + "1"
            full_source_manifest = next(
                row["full_source_manifest_sha256"]
                for row in crosswalk_rows
                if row["row_id"] == row_id
            )
            capacity_rows.append(
                {
                    "row_id": row_id,
                    "system": semantics["system"],
                    "population_state": (
                        amendment.EXPECTED_POPULATION_STATE[row_id]
                    ),
                    "cluster_proc": (
                        "5616584.0"
                        if row_id.startswith("pp_")
                        else "5616585.0"
                    ),
                    "exit_code": 0,
                    "num_job_starts": 1,
                    "num_holds": 0,
                    "request_memory_mb": amendment.REQUEST_MEMORY_MB,
                    "memory_usage_mb": amendment.REQUEST_MEMORY_MB,
                    "resident_set_size_kb": (
                        amendment.REQUEST_MEMORY_MB * 1024
                    ),
                    "remote_wall_clock_seconds": 1.0,
                    "capacity_chunk_index": 1,
                    "corrected_partition_chunk_index": 0,
                    "segment": 1,
                    "tuple_count": amendment.GROUP_SIZE,
                    "staged_chunk_list": str(self.root / f"{row_id}.list"),
                    "staged_chunk_sha256": staged_sha,
                    "corrected_chunk_fingerprint_sha256": "1" * 64,
                    "corrected_execution_chunk_sha256": "2" * 64,
                    "full_source_manifest_sha256": full_source_manifest,
                    "source_contract": source_contract,
                    "source_identity_canonical_sha256": source_identity,
                    "source_occurrence_id_hex": source_occurrence,
                    "base_config_sha256": semantics[
                        "base_config_sha256"
                    ],
                    "capacity_runtime_config": pinned(
                        f"{row_id}-runtime-config.yaml"
                    ),
                    "capacity_materialized_config": pinned(
                        f"{row_id}-materialized-config.yaml"
                    ),
                    "capacity_fanout_contract": pinned(
                        f"{row_id}-fanout.txt"
                    ),
                    "runtime_config_chain_semantic_sha256": "4" * 64,
                    "snapshot_loader_receipt": pinned(
                        f"{row_id}-snapshot-loader.json",
                        schema=amendment.SNAPSHOT_LOADER_SCHEMA,
                    ),
                    "snapshot_manifest": pinned(
                        f"{row_id}-snapshot-manifest.json",
                        schema=amendment.SNAPSHOT_MANIFEST_SCHEMA,
                    ),
                    "capacity_submission_manifest_row_sha256": "5" * 64,
                    "analysis_output": {
                        "path": str(self.root / f"{row_id}-analysis.root"),
                        "sha256": "6" * 64,
                        "size_bytes": 1,
                    },
                    "sidecar_output": {
                        "path": str(self.root / f"{row_id}-sidecar.root"),
                        "sha256": "7" * 64,
                        "size_bytes": 1,
                    },
                    "source_provenance": pinned(f"{row_id}-provenance.json"),
                    "staged_chunk_equals_corrected_partition_chunk0": True,
                }
            )

        payload: dict[str, object] = {
            "schema": amendment.AMENDMENT_SCHEMA,
            "status": "PASS",
            "authority_state": amendment.AUTHORITY_STATE,
            "generated_by": amendment.GENERATED_BY,
            "immutable_authority": {
                "status": "PASS",
                "bundle_manifest": immutable_bundle,
                "materialization_receipt": immutable_materialization,
                "bundle_identity_sha256": "8" * 64,
                "bundle_semantic_fingerprint_sha256": "9" * 64,
                "public_commit": "a" * 40,
                "code_sha256": "a" * 64,
                "replay_schema_sha256": "b" * 64,
                "training_schema_sha256": "c" * 64,
                "semantic_sha256": "d" * 64,
                "runtime": {
                    "release": "ana.560",
                    "offline_main": "/cvmfs/offline",
                    "calo_reco_soname": "libcalo_reco.so.0",
                    "request_memory_mb": amendment.REQUEST_MEMORY_MB,
                    "release_core_lib_dir": "/cvmfs/release/lib",
                    "release_core_lib64_dir": "/cvmfs/release/lib64",
                },
                "bundle_artifacts": immutable_artifacts,
                "materialized_bundle_path": str(self.root / "bundle"),
                "artifact_count": 1,
                "total_artifact_bytes": 1,
                "readback": "PASS_SYMLINK_FREE_READONLY_CONTENT_EXACT",
            },
            "legacy_preflight": legacy_preflight,
            "corrected_preflight": corrected_preflight,
            "source_semantic_crosswalk": source_crosswalk,
            "capacity_evidence": capacity_evidence,
            "count_correction": {
                "basis": (
                    "ordered_disjoint_group_of_seven_partition_one_execution_"
                    "and_source_occurrence_per_output_pair"
                ),
                "group_size": amendment.GROUP_SIZE,
                "source_tuple_count": amendment.EXPECTED_SOURCE_TUPLES,
                "legacy_expected_occurrence_count": (
                    amendment.EXPECTED_SOURCE_TUPLES
                ),
                "corrected_chunk_count": amendment.EXPECTED_CHUNKS,
                "corrected_job_count": amendment.EXPECTED_CHUNKS,
                "corrected_output_pair_count": amendment.EXPECTED_CHUNKS,
                "corrected_analysis_output_count": amendment.EXPECTED_CHUNKS,
                "corrected_sidecar_output_count": amendment.EXPECTED_CHUNKS,
                "physical_root_artifact_count": (
                    amendment.EXPECTED_ROOT_ARTIFACTS
                ),
                "corrected_source_occurrence_count": (
                    amendment.EXPECTED_CHUNKS
                ),
                "source_occurrences_per_output_pair": 1,
                "legacy_occurrence_overstatement": (
                    amendment.LEGACY_OVERSTATEMENT
                ),
                "science_semantics_changed": False,
                "tolerances_changed": False,
            },
            "capacity_to_partition_binding": {
                "status": "PASS",
                "reuse_scope": "group_size_7_resource_runtime_only",
                "corrected_execution_partition_sha256": "1" * 64,
                "corrected_partition_artifact_sha256": corrected_preflight[
                    "partition"
                ]["sha256"],
                "selected_rows": capacity_rows,
                "capacity_authority_earned": True,
                "full_training_authority": 0,
                "full_extraction_authority": False,
                "broad_production_authority": False,
            },
            "checks": {key: True for key in amendment.CHECK_KEYS},
            "submission_performed": False,
            "jobs_rerun": 0,
            "science_rerun_required": False,
            "capacity_reuse_scope": "group_size_7_resource_runtime_only",
            "full_training_authority": 0,
            "full_extraction_authority": False,
            "broad_production_authority": False,
            "canonical_promotion": False,
        }
        payload["capacity_to_partition_binding"][
            "corrected_execution_partition_sha256"
        ] = corrected_preflight["execution_partition_sha256"]
        payload["amendment_semantic_sha256"] = amendment.semantic_sha256(
            payload
        )
        return payload

    def test_source_v1_to_v2_crosswalk_is_exact(self) -> None:
        legacy, corrected = self.make_source_manifests()
        crosswalk = amendment.compare_source_manifests(legacy, corrected)
        self.assertEqual(crosswalk["source_tuple_count"], 129_998)
        self.assertEqual(crosswalk["row_count"], 13)
        self.assertTrue(
            all(
                row["legacy_tuple_aliases_preserved_as_history_only"]
                for row in crosswalk["rows"]
            )
        )

    def test_source_v2_alias_leak_is_rejected(self) -> None:
        legacy, corrected = self.make_source_manifests()
        corrected["rows"][0]["expected_occurrence_count"] = 10_000
        with self.assertRaisesRegex(
            amendment.AmendmentError, "leaks V1 count aliases"
        ):
            amendment.compare_source_manifests(legacy, corrected)

    def test_source_semantic_drift_is_rejected(self) -> None:
        legacy, corrected = self.make_source_manifests()
        corrected["rows"][3]["tuple_records_sha256"] = "f" * 64
        with self.assertRaisesRegex(
            amendment.AmendmentError, "source payload differs"
        ):
            amendment.compare_source_manifests(legacy, corrected)

    def test_first_chunk_materialization_is_byte_exact(self) -> None:
        row_root = self.root / "row"
        row_root.mkdir()
        input_lists = []
        tuples = []
        for role_index, role in enumerate(amendment.resolver.LIST_ROLES):
            path = row_root / amendment.resolver.LIST_FILENAMES[role]
            lines = [
                f"/{role}/input_{index}_{role_index}.root"
                for index in range(amendment.GROUP_SIZE)
            ]
            path.write_text("\n".join(lines) + "\n", encoding="utf-8")
            input_lists.append({"role": role, "path": str(path)})
        for index in range(amendment.GROUP_SIZE):
            tuples.append(
                {
                    "tuple_index": index,
                    "physical_line": index + 1,
                    "inputs": {
                        role: f"/{role}/input_{index}_{role_index}.root"
                        for role_index, role in enumerate(
                            amendment.resolver.LIST_ROLES
                        )
                    },
                }
            )
        partition, chunks = amendment.resolver.build_source_partition(
            "pp_background_jet8",
            tuples,
            full_source_manifest_sha256="d" * 64,
        )
        source = {
            "row_id": "pp_background_jet8",
            "input_lists": input_lists,
            "partition_contract": partition,
            "_chunks": chunks,
        }
        observed, chunk = amendment.materialize_first_chunk_bytes(source)
        expected = b"".join(
            (
                "\t".join(
                    f"/{role}/input_{index}_{role_index}.root"
                    for role_index, role in enumerate(
                        amendment.resolver.LIST_ROLES
                    )
                )
                + "\n"
            ).encode("utf-8")
            for index in range(amendment.GROUP_SIZE)
        )
        self.assertEqual(observed, expected)
        self.assertEqual(chunk["chunk_index"], 0)
        self.assertEqual(chunk["segment"], 1)
        raw_first_tuple_sha256 = hashlib.sha256(
            expected.partition(b"\n")[0] + b"\n"
        ).hexdigest()
        semantic_first_tuple_sha256 = (
            amendment.resolver.canonical_sha256(tuples[0])
        )
        self.assertEqual(
            amendment.first_tuple_tsv_sha256(observed),
            raw_first_tuple_sha256,
        )
        self.assertNotEqual(
            raw_first_tuple_sha256,
            semantic_first_tuple_sha256,
        )

    def test_source_identity_matches_controller_contract(self) -> None:
        contract = {
            "lane": "pp_inclusive_sim",
            "dataset": "isSimInclusive",
            "sample": "run28_jet8",
            "period": "0mrad",
            "run": 28,
            "segment": 1,
            "input_uri_hash": "1" * 64,
            "input_file_sha256": "1" * 64,
            "source_manifest_sha256": "2" * 64,
        }
        canonical = "|".join(
            str(contract[field])
            for field in (
                "lane",
                "dataset",
                "sample",
                "period",
                "run",
                "segment",
                "input_uri_hash",
                "input_file_sha256",
                "source_manifest_sha256",
            )
        )
        self.assertEqual(
            amendment.canonical_source_identity_sha256(contract),
            hashlib.sha256(canonical.encode("utf-8")).hexdigest(),
        )

    def test_plan_authority_inventory_rejects_additive_grants(self) -> None:
        payload = {
            "submission_performed": False,
            "authority": {
                "requested_scope": "full",
                "full_training_authority": 0,
                "authority_state": "PREFLIGHT_RESOLVED_NOT_EARNED",
            },
        }
        amendment.require_preflight_authority("plan", payload)
        for extra_field, value in (
            ("full_extraction_authority", True),
            ("broad_production_authority", True),
            ("canonical_promotion", True),
        ):
            with self.subTest(extra_field=extra_field):
                mutated = copy.deepcopy(payload)
                mutated["authority"][extra_field] = value
                with self.assertRaisesRegex(
                    amendment.AmendmentError, "field inventory differs"
                ):
                    amendment.require_preflight_authority("plan", mutated)

    def test_capacity_audit_non_training_authority_is_exact(self) -> None:
        audit_path = self.root / "pp_capacity_audit.json"
        audit_path.write_text("{}\n", encoding="utf-8")
        expected_matrix = (
            self.root / "pp_h70_source_complete_smoke_matrix.npz"
        )
        payload = {
            "scope": "capacity",
            "authority_state": "CAPACITY_NON_TRAINING_ONLY",
            "row_id": "pp_background_jet8",
            "system": "pp",
            "class_balance_state": (
                "NOT_APPLICABLE_CAPACITY_NON_TRAINING"
            ),
            "failures": [],
            "source_occurrence_count": 1,
            "full_training_authority": 0,
            "matrix_materialized": False,
            "training_matrix_authority_earned": False,
            "matrix_out": str(expected_matrix),
            "checks": {
                key: True for key in amendment.CAPACITY_AUDIT_CHECK_KEYS
            },
        }
        artifact = {"path": str(audit_path)}
        amendment.validate_capacity_non_training_audit(
            "pp_background_jet8", payload, artifact
        )
        mutations = {
            "scope": ("scope", "training"),
            "authority": ("authority_state", "TRAINING"),
            "occurrences": ("source_occurrence_count", 2),
            "training": ("full_training_authority", 1),
            "materialized": ("matrix_materialized", True),
            "matrix_authority": (
                "training_matrix_authority_earned",
                True,
            ),
            "matrix_path": ("matrix_out", str(self.root / "other.npz")),
            "class_balance": ("class_balance_state", "BALANCED"),
            "failures": ("failures", ["failure"]),
        }
        for label, (field, value) in mutations.items():
            with self.subTest(label=label):
                mutated = copy.deepcopy(payload)
                mutated[field] = value
                with self.assertRaises(amendment.AmendmentError):
                    amendment.validate_capacity_non_training_audit(
                        "pp_background_jet8", mutated, artifact
                    )
        mutated = copy.deepcopy(payload)
        mutated["checks"]["matrix_absent"] = False
        with self.assertRaises(amendment.AmendmentError):
            amendment.validate_capacity_non_training_audit(
                "pp_background_jet8", mutated, artifact
            )
        mutated = copy.deepcopy(payload)
        mutated["checks"]["unexpected"] = True
        with self.assertRaisesRegex(
            amendment.AmendmentError, "field inventory differs"
        ):
            amendment.validate_capacity_non_training_audit(
                "pp_background_jet8", mutated, artifact
            )

    def test_capacity_resource_memory_ceiling_is_fail_closed(self) -> None:
        row = {
            "row_id": "pp_background_jet8",
            "cluster_proc": "5616584.0",
            "job_status": 4,
            "exit_code": 0,
            "num_holds": 0,
            "num_job_starts": 1,
            "request_memory_mb": amendment.REQUEST_MEMORY_MB,
            "memory_usage_mb": amendment.REQUEST_MEMORY_MB,
            "memory_usage_resolution": "EVALUATED_FROM_RESIDENT_SET_SIZE_KB",
            "resident_set_size_kb": amendment.REQUEST_MEMORY_MB * 1024,
            "remote_wall_clock_seconds": 1.0,
        }
        self.assertEqual(
            amendment.validate_capacity_resource_row(
                "pp_background_jet8", row
            ),
            row,
        )
        for field, value in (
            ("memory_usage_mb", amendment.REQUEST_MEMORY_MB + 1),
            (
                "resident_set_size_kb",
                amendment.REQUEST_MEMORY_MB * 1024 + 1,
            ),
        ):
            mutated = copy.deepcopy(row)
            mutated[field] = value
            with self.assertRaisesRegex(
                amendment.AmendmentError, field
            ):
                amendment.validate_capacity_resource_row(
                    "pp_background_jet8", mutated
                )

    def test_capacity_source_binding_rejects_plan_drift(self) -> None:
        semantics = {
            "lane": "pp_inclusive_sim",
            "dataset": "isSimInclusive",
            "sample": "run28_jet8",
            "period": "0mrad",
            "run": 28,
            "si_di_role": "SI",
            "ownership_state": "source_role_frozen",
        }
        source_contract = {
            **semantics,
            "segment": 1,
            "input_uri_hash": "1" * 64,
            "input_file_sha256": "1" * 64,
            "source_manifest_sha256": "2" * 64,
        }
        source_execution = {
            "args_file": str(self.root / "args.txt"),
            "chunk_index": 1,
            "run": 28,
            "staged_chunk_list": str(self.root / "chunk.list"),
            "submitted_args_sha256": "3" * 64,
        }
        amendment.validate_capacity_source_binding(
            "pp_background_jet8",
            source_contract,
            source_execution,
            semantics,
        )
        for field in (
            "lane",
            "dataset",
            "sample",
            "period",
            "run",
            "si_di_role",
            "ownership_state",
        ):
            mutated = copy.deepcopy(source_contract)
            mutated[field] = 29 if field == "run" else "drift"
            with self.assertRaisesRegex(
                amendment.AmendmentError, field
            ):
                amendment.validate_capacity_source_binding(
                    "pp_background_jet8",
                    mutated,
                    source_execution,
                    semantics,
                )
        mutated_execution = copy.deepcopy(source_execution)
        mutated_execution["run"] = 29
        with self.assertRaisesRegex(amendment.AmendmentError, "execution run"):
            amendment.validate_capacity_source_binding(
                "pp_background_jet8",
                source_contract,
                mutated_execution,
                semantics,
            )

    def test_source_provenance_requires_exact_file_and_embedded_record(
        self,
    ) -> None:
        sidecar = self.root / "sidecar.root"
        sidecar.write_bytes(b"ROOT")
        provenance_path = self.root / "source_provenance.json"
        staged_sha = "1" * 64
        source_contract = {"source_manifest_sha256": "2" * 64}
        semantics = {
            "system": "pp",
            "sample": "run28_jet8",
            "base_config_sha256": "9" * 64,
            "code_sha256": "4" * 64,
        }
        submission_manifest_row = {
            "resolved_config_sha256": "3" * 64,
            "code_sha256": "4" * 64,
        }
        record = {
            "path": str(sidecar),
            "row_id": "pp_background_jet8",
            "system": "pp",
            "source_sample": "run28_jet8",
            "input_uri_sha256": staged_sha,
            "input_file_sha256": staged_sha,
            "source_manifest_sha256": "2" * 64,
            "config_sha256": "3" * 64,
            "code_sha256": "4" * 64,
        }
        provenance_path.write_text(
            json.dumps(
                {
                    "schema": "THE134_SOURCE_PROVENANCE_V1",
                    "inputs": [record],
                },
                sort_keys=True,
            )
            + "\n",
            encoding="utf-8",
        )
        wrapper = {
            "path": str(provenance_path),
            "record": copy.deepcopy(record),
            "sha256": self.file_sha256(provenance_path),
        }
        amendment.validate_source_provenance_artifact(
            "pp_background_jet8",
            wrapper,
            sidecar_path=sidecar,
            staged_sha256=staged_sha,
            source_contract=source_contract,
            corrected_semantics=semantics,
            submission_manifest_row=submission_manifest_row,
        )
        mutated_wrapper = copy.deepcopy(wrapper)
        mutated_wrapper["record"]["config_sha256"] = "5" * 64
        with self.assertRaisesRegex(
            amendment.AmendmentError, "content differs"
        ):
            amendment.validate_source_provenance_artifact(
                "pp_background_jet8",
                mutated_wrapper,
                sidecar_path=sidecar,
                staged_sha256=staged_sha,
                source_contract=source_contract,
                corrected_semantics=semantics,
                submission_manifest_row=submission_manifest_row,
            )
        payload = json.loads(provenance_path.read_text(encoding="utf-8"))
        payload["inputs"].append(copy.deepcopy(record))
        provenance_path.write_text(
            json.dumps(payload, sort_keys=True) + "\n", encoding="utf-8"
        )
        mutated_wrapper = {
            **wrapper,
            "sha256": self.file_sha256(provenance_path),
        }
        with self.assertRaisesRegex(
            amendment.AmendmentError, "exactly one input"
        ):
            amendment.validate_source_provenance_artifact(
                "pp_background_jet8",
                mutated_wrapper,
                sidecar_path=sidecar,
                staged_sha256=staged_sha,
                source_contract=source_contract,
                corrected_semantics=semantics,
                submission_manifest_row=submission_manifest_row,
            )
        for label, serialized, expected_error in (
            (
                "duplicate_key",
                (
                    '{"schema":"WRONG","schema":'
                    '"THE134_SOURCE_PROVENANCE_V1","inputs":'
                    f"{json.dumps([record], sort_keys=True)}}}\n"
                ),
                "duplicates key",
            ),
            (
                "nonfinite",
                (
                    '{"schema":"THE134_SOURCE_PROVENANCE_V1",'
                    f'"inputs":{json.dumps([record], sort_keys=True)},'
                    '"unexpected":NaN}\n'
                ),
                "non-finite constant",
            ),
        ):
            with self.subTest(label=label):
                provenance_path.write_text(serialized, encoding="utf-8")
                malformed_wrapper = {
                    **wrapper,
                    "sha256": self.file_sha256(provenance_path),
                }
                with self.assertRaisesRegex(
                    amendment.AmendmentError, expected_error
                ):
                    amendment.validate_source_provenance_artifact(
                        "pp_background_jet8",
                        malformed_wrapper,
                        sidecar_path=sidecar,
                        staged_sha256=staged_sha,
                        source_contract=source_contract,
                        corrected_semantics=semantics,
                        submission_manifest_row=submission_manifest_row,
                    )

    def test_submission_manifest_runtime_config_is_independent_and_pinned(
        self,
    ) -> None:
        runtime_config = self.root / "resolved.yaml"
        runtime_config.write_text("runtime: exact\n", encoding="utf-8")
        release_lib = self.root / "release" / "lib"
        release_lib64 = self.root / "release" / "lib64"
        release_lib.mkdir(parents=True)
        release_lib64.mkdir(parents=True)
        runtime_fields = (
            "library",
            "model",
            "photon_cluster_builder_header",
            "calo_reco_library",
            "release_calo_io",
            "release_clusteriso",
            "release_jetbase",
        )
        runtime_paths: dict[str, Path] = {}
        for field in runtime_fields:
            path = self.root / field
            path.write_bytes(field.encode("utf-8"))
            runtime_paths[field] = path
        first_tuple_inputs = {
            role: f"/{role}/input.root"
            for role in amendment.resolver.LIST_ROLES
        }
        raw_first_tuple_sha256 = hashlib.sha256(
            (
                "\t".join(
                    first_tuple_inputs[role]
                    for role in amendment.resolver.LIST_ROLES
                )
                + "\n"
            ).encode("utf-8")
        ).hexdigest()
        semantic_first_tuple_sha256 = amendment.resolver.canonical_sha256(
            {
                "tuple_index": 0,
                "physical_line": 1,
                "inputs": first_tuple_inputs,
            }
        )
        self.assertNotEqual(
            raw_first_tuple_sha256,
            semantic_first_tuple_sha256,
        )
        row = {
            field: "value"
            for field in amendment.SUBMISSION_MANIFEST_FIELDS
        }
        row.update(
            {
                "row_id": "pp_background_jet8",
                "system": "pp",
                "lane": "pp_inclusive_sim",
                "dataset": "isSimInclusive",
                "sample": "run28_jet8",
                "source_role": "background",
                "minimum_bias_gate": "not_applicable",
                "input_files": "7",
                "input_jobs": "1",
                "source_manifest_sha256": "2" * 64,
                "first_input_tuple_sha256": raw_first_tuple_sha256,
                "resolved_config": str(runtime_config),
                "resolved_config_sha256": self.file_sha256(runtime_config),
                "code_sha256": "4" * 64,
                "replay_schema_sha256": "7" * 64,
                "training_schema_sha256": "8" * 64,
                "semantic_sha256": "9" * 64,
                "nominal_et_min_gev": "15.0",
                "nominal_et_max_gev_exclusive": "35.0",
                "loose_capture_et_min_gev": "5.0",
                "legacy_training_tree_max_entries": "0",
                "event_limit_per_job": "0",
                "executable_input_tuple_count": "10000",
                "full_training_authority": "0",
                "release_core_lib_dir": str(release_lib),
                "release_core_lib64_dir": str(release_lib64),
            }
        )
        for field, path in runtime_paths.items():
            row[field] = str(path)
            row[f"{field}_sha256"] = self.file_sha256(path)
        source_contract = {"source_manifest_sha256": "2" * 64}
        semantics = {
            "system": "pp",
            "lane": "pp_inclusive_sim",
            "dataset": "isSimInclusive",
            "sample": "run28_jet8",
            "source_role": "background",
            "minimum_bias_gate": "not_applicable",
            "code_sha256": "4" * 64,
            "base_config_sha256": "9" * 64,
        }
        role_by_field = {
            "library": "pp_library",
            "model": "pp_model",
            "photon_cluster_builder_header":
                "photon_cluster_builder_header",
            "calo_reco_library": "calo_reco_library",
            "release_calo_io": "release_calo_io",
            "release_clusteriso": "release_clusteriso",
            "release_jetbase": "release_jetbase",
        }
        immutable_artifacts = {}
        for field, role in role_by_field.items():
            authority_path = runtime_paths[field]
            if field not in {"library", "photon_cluster_builder_header"}:
                authority_path = self.root / "bundle-copy" / field
                authority_path.parent.mkdir(exist_ok=True)
                authority_path.write_bytes(runtime_paths[field].read_bytes())
            immutable_artifacts[role] = {
                "role": role,
                "path": str(authority_path),
                "resolved_path": str(authority_path),
                "sha256": self.file_sha256(authority_path),
                "size_bytes": authority_path.stat().st_size,
            }
        immutable = {
            "code_sha256": "4" * 64,
            "replay_schema_sha256": "7" * 64,
            "training_schema_sha256": "8" * 64,
            "semantic_sha256": "9" * 64,
            "runtime": {
                "release_core_lib_dir": str(release_lib),
                "release_core_lib64_dir": str(release_lib64),
            },
            "bundle_artifacts": immutable_artifacts,
        }
        runtime_authority = {
            "calo_reco": {
                "library": {
                    "path": str(runtime_paths["calo_reco_library"]),
                },
            },
            "release_dirs": {
                "lib": str(release_lib),
                "lib64": str(release_lib64),
            },
            "providers": {
                family: {
                    "path": str(runtime_paths[field]),
                }
                for family, field in {
                    "libcalo_io.so": "release_calo_io",
                    "libclusteriso.so": "release_clusteriso",
                    "libjetbase.so": "release_jetbase",
                }.items()
            },
        }
        binding = amendment.validate_submission_manifest_record(
            "pp_background_jet8",
            row,
            source_contract=source_contract,
            corrected_semantics=semantics,
            corrected_source={
                "tuple_count": 10_000,
                "first_tuple_sha256": semantic_first_tuple_sha256,
            },
            expected_first_input_tuple_sha256=raw_first_tuple_sha256,
            immutable=immutable,
            runtime_authority=runtime_authority,
        )
        self.assertEqual(
            binding["runtime_config"]["sha256"],
            self.file_sha256(runtime_config),
        )
        self.assertNotEqual(
            binding["runtime_config"]["sha256"],
            semantics["base_config_sha256"],
        )
        with self.assertRaisesRegex(
            amendment.AmendmentError,
            "first-tuple authority differs",
        ):
            amendment.validate_submission_manifest_record(
                "pp_background_jet8",
                row,
                source_contract=source_contract,
                corrected_semantics=semantics,
                corrected_source={
                    "tuple_count": 10_000,
                    "first_tuple_sha256": semantic_first_tuple_sha256,
                },
                expected_first_input_tuple_sha256=(
                    semantic_first_tuple_sha256
                ),
                immutable=immutable,
                runtime_authority=runtime_authority,
            )
        copied_library = self.root / "copied-analysis-library.so"
        copied_library.write_bytes(runtime_paths["library"].read_bytes())
        row_with_copied_library = copy.deepcopy(row)
        row_with_copied_library["library"] = str(copied_library)
        row_with_copied_library["library_sha256"] = self.file_sha256(
            copied_library
        )
        with self.assertRaisesRegex(
            amendment.AmendmentError, "different files"
        ):
            amendment.validate_submission_manifest_record(
                "pp_background_jet8",
                row_with_copied_library,
                source_contract=source_contract,
                corrected_semantics=semantics,
                corrected_source={
                    "tuple_count": 10_000,
                    "first_tuple_sha256": semantic_first_tuple_sha256,
                },
                expected_first_input_tuple_sha256=raw_first_tuple_sha256,
                immutable=immutable,
                runtime_authority=runtime_authority,
            )
        runtime_config.write_text("runtime: drift\n", encoding="utf-8")
        with self.assertRaisesRegex(
            amendment.AmendmentError, "hash drift"
        ):
            amendment.validate_submission_manifest_record(
                "pp_background_jet8",
                row,
                source_contract=source_contract,
                corrected_semantics=semantics,
                corrected_source={
                    "tuple_count": 10_000,
                    "first_tuple_sha256": semantic_first_tuple_sha256,
                },
                expected_first_input_tuple_sha256=raw_first_tuple_sha256,
                immutable=immutable,
                runtime_authority=runtime_authority,
            )

    def test_amendment_payload_is_deterministic_and_fail_closed(self) -> None:
        payload = self.make_valid_amendment_payload()
        first = amendment.validate_amendment_payload(payload)
        second = amendment.validate_amendment_payload(copy.deepcopy(payload))
        self.assertEqual(
            amendment.canonical_json_bytes(first),
            amendment.canonical_json_bytes(second),
        )
        for field, replacement in (
            ("submission_performed", True),
            ("jobs_rerun", 1),
            ("science_rerun_required", True),
            ("full_training_authority", 1),
            ("full_extraction_authority", True),
            ("broad_production_authority", True),
            ("canonical_promotion", True),
        ):
            mutated = copy.deepcopy(payload)
            mutated[field] = replacement
            mutated["amendment_semantic_sha256"] = amendment.semantic_sha256(
                {
                    key: value
                    for key, value in mutated.items()
                    if key != "amendment_semantic_sha256"
                }
            )
            with self.assertRaises(amendment.AmendmentError):
                amendment.validate_amendment_payload(mutated)

    def test_incomplete_or_additive_nested_evidence_is_rejected(self) -> None:
        payload = self.make_valid_amendment_payload()
        incomplete = copy.deepcopy(payload)
        incomplete["legacy_preflight"] = {"status": "PASS"}
        incomplete["amendment_semantic_sha256"] = amendment.semantic_sha256(
            {
                key: value
                for key, value in incomplete.items()
                if key != "amendment_semantic_sha256"
            }
        )
        with self.assertRaisesRegex(
            amendment.AmendmentError, "field inventory differs"
        ):
            amendment.validate_amendment_payload(incomplete)

        additive = copy.deepcopy(payload)
        additive["corrected_preflight"]["unexpected_authority"] = True
        additive["amendment_semantic_sha256"] = amendment.semantic_sha256(
            {
                key: value
                for key, value in additive.items()
                if key != "amendment_semantic_sha256"
            }
        )
        with self.assertRaisesRegex(
            amendment.AmendmentError, "field inventory differs"
        ):
            amendment.validate_amendment_payload(additive)

    def test_unknown_field_and_self_hash_drift_are_rejected(self) -> None:
        payload = self.make_valid_amendment_payload()
        mutated = copy.deepcopy(payload)
        mutated["unexpected"] = True
        with self.assertRaisesRegex(
            amendment.AmendmentError, "field inventory differs"
        ):
            amendment.validate_amendment_payload(mutated)
        mutated = copy.deepcopy(payload)
        mutated["count_correction"]["corrected_job_count"] += 1
        with self.assertRaises(amendment.AmendmentError):
            amendment.validate_amendment_payload(mutated)

    def test_write_once_allows_only_byte_identical_readback(self) -> None:
        payload = self.make_valid_amendment_payload()
        output = self.root / "amendment.json"
        self.assertEqual(
            amendment.write_once_readonly_json(output, payload),
            "CREATED_READONLY",
        )
        self.assertEqual(stat.S_IMODE(output.stat().st_mode), 0o444)
        self.assertEqual(
            amendment.write_once_readonly_json(output, payload),
            "EXISTING_BYTE_IDENTICAL",
        )
        os.chmod(output, 0o644)
        output.write_text("{}\n", encoding="utf-8")
        with self.assertRaisesRegex(
            amendment.AmendmentError, "refusing overwrite"
        ):
            amendment.write_once_readonly_json(output, payload)

    def test_write_once_rejects_existing_external_hardlink(self) -> None:
        payload = self.make_valid_amendment_payload()
        external = self.root / "external.json"
        external.write_bytes(amendment.canonical_json_bytes(payload))
        os.chmod(external, 0o600)
        output = self.root / "amendment-hardlink.json"
        os.link(external, output)
        with self.assertRaisesRegex(
            amendment.AmendmentError, "exactly one hard link"
        ):
            amendment.write_once_readonly_json(output, payload)
        self.assertEqual(stat.S_IMODE(external.stat().st_mode), 0o600)

    def test_concurrent_creator_is_never_overwritten(self) -> None:
        payload = self.make_valid_amendment_payload()
        output = self.root / "race.json"
        competing_bytes = b'{"owner":"other"}\n'

        def create_competing_destination(
            _source: object,
            destination: object,
            *,
            follow_symlinks: bool,
        ) -> None:
            self.assertFalse(follow_symlinks)
            Path(destination).write_bytes(competing_bytes)
            raise FileExistsError(str(destination))

        with (
            mock.patch.object(
                amendment.os,
                "link",
                side_effect=create_competing_destination,
            ),
            self.assertRaisesRegex(
                amendment.AmendmentError, "concurrent amendment differs"
            ),
        ):
            amendment.write_once_readonly_json(output, payload)
        self.assertEqual(output.read_bytes(), competing_bytes)

    def test_concurrent_external_hardlink_is_never_chmodded(self) -> None:
        payload = self.make_valid_amendment_payload()
        output = self.root / "race-hardlink.json"
        external = self.root / "race-external.json"
        external.write_bytes(amendment.canonical_json_bytes(payload))
        os.chmod(external, 0o600)
        real_link = os.link

        def create_competing_hardlink(
            _source: object,
            destination: object,
            *,
            follow_symlinks: bool,
        ) -> None:
            self.assertFalse(follow_symlinks)
            real_link(external, destination)
            raise FileExistsError(str(destination))

        with (
            mock.patch.object(
                amendment.os,
                "link",
                side_effect=create_competing_hardlink,
            ),
            self.assertRaisesRegex(
                amendment.AmendmentError, "exactly one hard link"
            ),
        ):
            amendment.write_once_readonly_json(output, payload)
        self.assertEqual(stat.S_IMODE(external.stat().st_mode), 0o600)

    def test_verify_requires_raw_canonical_bytes_and_is_idempotent(
        self,
    ) -> None:
        payload = self.make_valid_amendment_payload()
        spec_path = self.root / "spec.json"
        spec_path.write_text("{}\n", encoding="utf-8")
        amendment_path = self.root / "amendment-readback.json"
        amendment_path.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        pretty_sha = self.file_sha256(amendment_path)
        common_args = [
            "verify",
            "--spec",
            str(spec_path),
            "--expected-spec-sha256",
            "1" * 64,
            "--amendment",
            str(amendment_path),
            "--expected-amendment-sha256",
            pretty_sha,
        ]
        with (
            mock.patch.object(amendment, "load_spec", return_value={}),
            mock.patch.object(
                amendment, "build_amendment", return_value=payload
            ),
            mock.patch.object(
                amendment, "load_amendment", return_value=payload
            ),
            contextlib.redirect_stdout(io.StringIO()),
            contextlib.redirect_stderr(io.StringIO()),
        ):
            self.assertEqual(amendment.main(common_args), 2)

        amendment_path.write_bytes(amendment.canonical_json_bytes(payload))
        canonical_sha = self.file_sha256(amendment_path)
        output_json = self.root / "verify-readback.json"
        canonical_args = [
            *common_args[:-1],
            canonical_sha,
            "--output-json",
            str(output_json),
            "--allowed-output-root",
            str(self.root),
        ]
        with (
            mock.patch.object(amendment, "load_spec", return_value={}),
            mock.patch.object(
                amendment, "build_amendment", return_value=payload
            ),
            mock.patch.object(
                amendment, "load_amendment", return_value=payload
            ),
            contextlib.redirect_stdout(io.StringIO()),
            contextlib.redirect_stderr(io.StringIO()),
        ):
            self.assertEqual(amendment.main(canonical_args), 0)
            self.assertEqual(amendment.main(canonical_args), 0)
        self.assertEqual(stat.S_IMODE(output_json.stat().st_mode), 0o444)

    def test_symlink_inputs_and_outputs_are_rejected(self) -> None:
        target = self.root / "target.json"
        target.write_text("{}\n", encoding="utf-8")
        link = self.root / "link.json"
        link.symlink_to(target)
        with self.assertRaisesRegex(amendment.AmendmentError, "symlink"):
            amendment.normalize_artifact_ref(
                "linked evidence",
                {"path": str(link), "sha256": self.file_sha256(target)},
            )
        with self.assertRaisesRegex(amendment.AmendmentError, "symlink"):
            amendment.write_once_readonly_json(
                link, self.make_valid_amendment_payload()
            )

    def test_declared_site_alias_root_is_allowed_but_escape_is_rejected(
        self,
    ) -> None:
        physical_root = self.root / "physical"
        evidence_root = physical_root / "evidence"
        evidence_root.mkdir(parents=True)
        alias_root = self.root / "site_alias"
        alias_root.symlink_to(physical_root, target_is_directory=True)
        output = alias_root / "evidence" / "amendment.json"
        self.assertEqual(
            amendment.write_once_readonly_json(
                output,
                self.make_valid_amendment_payload(),
                allowed_output_root=alias_root / "evidence",
            ),
            "CREATED_READONLY",
        )
        self.assertTrue((evidence_root / "amendment.json").is_file())

        outside = self.root / "outside"
        outside.mkdir()
        escape = alias_root / "evidence" / "escape"
        escape.symlink_to(outside, target_is_directory=True)
        with self.assertRaisesRegex(
            amendment.AmendmentError, "resolves outside"
        ):
            amendment.write_once_readonly_json(
                escape / "bad.json",
                self.make_valid_amendment_payload(),
                allowed_output_root=alias_root / "evidence",
            )

    def test_resolved_file_identity_accepts_site_alias_spelling(self) -> None:
        physical_root = self.root / "physical"
        physical_root.mkdir()
        artifact = physical_root / "evidence.json"
        artifact.write_text("{}\n", encoding="utf-8")
        alias_root = self.root / "site_alias"
        alias_root.symlink_to(physical_root, target_is_directory=True)
        left, right = amendment.require_same_resolved_file(
            "site alias",
            alias_root / "evidence.json",
            artifact,
        )
        self.assertEqual(left.resolve(), right.resolve())
        other = self.root / "other.json"
        other.write_text("{}\n", encoding="utf-8")
        with self.assertRaisesRegex(
            amendment.AmendmentError, "different files"
        ):
            amendment.require_same_resolved_file(
                "site alias",
                alias_root / "evidence.json",
                other,
            )

    def test_builder_output_never_grants_execution_authority(self) -> None:
        validated_fixture = self.make_valid_amendment_payload()
        legacy_plan = {
            "schema": amendment.LEGACY_PLAN_SCHEMA,
            "rows": [
                {
                    "row_id": row_id,
                    "system": "pp" if row_id.startswith("pp_") else "auau",
                    "lane": "lane",
                    "dataset": "dataset",
                    "sample": "sample",
                    "source_role": "signal",
                    "minimum_bias_gate": "not_applicable",
                    "photon_id_row_match": "newPPG12",
                    "source_period": "0mrad",
                    "source_si_di_role": "SI",
                    "training_period_si_contract_sha256": "1" * 64,
                    "science_contract": {"model_domain": [15, 35]},
                    "bundle_contract": {"system": row_id.split("_", 1)[0]},
                }
                for row_id in amendment.expected_row_ids()
            ],
        }
        corrected_plan = copy.deepcopy(legacy_plan)
        corrected_plan["schema"] = amendment.CORRECTED_PLAN_SCHEMA
        files: dict[str, Path] = {}
        payloads = {
            "legacy_plan": legacy_plan,
            "legacy_receipt": {
                "schema": amendment.LEGACY_RECEIPT_SCHEMA
            },
            "legacy_source": {
                "schema": amendment.source_builder.SOURCE_SCHEMA_V1
            },
            "corrected_plan": corrected_plan,
            "corrected_receipt": {
                "schema": amendment.CORRECTED_RECEIPT_SCHEMA
            },
            "corrected_source": {
                "schema": amendment.CORRECTED_SOURCE_SCHEMA
            },
            "capacity": {"schema": amendment.CAPACITY_SCHEMA},
            "root": {"schema": amendment.ROOT_JOIN_SCHEMA},
            "pp_audit": {"schema": amendment.CAPACITY_AUDIT_SCHEMA},
            "auau_audit": {"schema": amendment.CAPACITY_AUDIT_SCHEMA},
            "bundle": {"schema": amendment.resolver.BUNDLE_SCHEMA},
            "materialization": {
                "schema": amendment.resolver.MATERIALIZATION_SCHEMA
            },
        }
        for name, payload in payloads.items():
            path = self.root / f"{name}.json"
            path.write_text(json.dumps(payload) + "\n", encoding="utf-8")
            files[name] = path
        rows_path = self.root / "rows.jsonl"
        rows_path.write_text("{}\n", encoding="utf-8")
        duplicate_path = self.root / "duplicate.sha256"
        duplicate_path.write_text("d" * 64 + "\n", encoding="utf-8")
        partition_path = self.root / "partition.jsonl"
        partition_path.write_text("{}\n", encoding="utf-8")
        spec = {
            "schema": amendment.SPEC_SCHEMA,
            "immutable_authority": {
                "bundle_manifest": self.artifact(files["bundle"]),
                "materialization_receipt": self.artifact(
                    files["materialization"]
                ),
            },
            "legacy": {
                "plan": self.artifact(files["legacy_plan"]),
                "preflight_receipt": self.artifact(files["legacy_receipt"]),
                "source_manifest": self.artifact(files["legacy_source"]),
            },
            "corrected": {
                "plan": self.artifact(files["corrected_plan"]),
                "preflight_receipt": self.artifact(
                    files["corrected_receipt"]
                ),
                "source_manifest": self.artifact(files["corrected_source"]),
                "rows": self.artifact(rows_path),
                "duplicate_fingerprint": self.artifact(duplicate_path),
                "partition": self.artifact(partition_path),
            },
            "capacity": {
                "resource_certificate": self.artifact(files["capacity"]),
                "root_join_certificate": self.artifact(files["root"]),
                "pp_audit": self.artifact(files["pp_audit"]),
                "auau_audit": self.artifact(files["auau_audit"]),
            },
        }
        legacy_validated = validated_fixture["legacy_preflight"]
        corrected_validated = validated_fixture["corrected_preflight"]
        with (
            mock.patch.object(
                amendment,
                "validate_immutable_authority",
                return_value=validated_fixture["immutable_authority"],
            ),
            mock.patch.object(
                amendment,
                "load_source_manifest_artifact",
                side_effect=[
                    (
                        payloads["legacy_source"],
                        {
                            **self.artifact(files["legacy_source"]),
                            "resolved_path": str(files["legacy_source"]),
                            "size_bytes": files["legacy_source"].stat().st_size,
                        },
                    ),
                    (
                        payloads["corrected_source"],
                        {
                            **self.artifact(files["corrected_source"]),
                            "resolved_path": str(files["corrected_source"]),
                            "size_bytes": files["corrected_source"].stat().st_size,
                        },
                    ),
                ],
            ),
            mock.patch.object(
                amendment,
                "validate_legacy_preflight",
                return_value=legacy_validated,
            ),
            mock.patch.object(
                amendment,
                "validate_preflight_immutable_binding",
                return_value=None,
            ),
            mock.patch.object(
                amendment,
                "validate_corrected_bundle_authority",
                return_value=None,
            ),
            mock.patch.object(
                amendment.resolver,
                "validate_source_manifest",
                return_value=[],
            ),
            mock.patch.object(
                amendment,
                "validate_corrected_preflight",
                return_value=corrected_validated,
            ),
            mock.patch.object(
                amendment,
                "compare_source_manifests",
                return_value=validated_fixture[
                    "source_semantic_crosswalk"
                ],
            ),
            mock.patch.object(
                amendment,
                "validate_capacity_evidence",
                return_value=(
                    validated_fixture["capacity_evidence"],
                    validated_fixture["capacity_to_partition_binding"][
                        "selected_rows"
                    ],
                ),
            ),
        ):
            first = amendment.build_amendment(spec)
        self.assertFalse(first["submission_performed"])
        self.assertEqual(first["jobs_rerun"], 0)
        self.assertFalse(first["full_extraction_authority"])
        self.assertFalse(first["broad_production_authority"])
        self.assertFalse(first["canonical_promotion"])

    def test_source_has_no_remote_mutation_or_audio_surface(self) -> None:
        source = MODULE_PATH.read_text(encoding="utf-8")
        for forbidden in (
            "subprocess.",
            "os.system(",
            "os.popen(",
            "os.spawn",
            "os.exec",
            "socket.",
            "condor_submit",
            "condor_rm",
            "condor_release",
            "osascript",
            "\a",
            "/usr/bin/say",
            "afplay",
            "NSSound",
        ):
            self.assertNotIn(forbidden, source)


if __name__ == "__main__":
    unittest.main()
