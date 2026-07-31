#!/usr/bin/env python3
"""Focused adversarial tests for the current-head capacity binding."""

from __future__ import annotations

import copy
import hashlib
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock


HERE = Path(__file__).resolve().parent
MODULE_PATH = HERE.parent / "build_the134_capacity_partition_binding.py"
PROJECTOR_PATH = HERE.parent / "project_the134_preextraction_storage_quota.py"


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


binding = load_module("the134_capacity_partition_binding_tested", MODULE_PATH)
projector = load_module("the134_capacity_projector_binding_tested", PROJECTOR_PATH)


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


class CapacityPartitionBindingTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name).resolve()

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def artifact(self, name: str, sha_char: str) -> dict[str, str]:
        return {
            "path": str(self.root / name),
            "sha256": sha_char * 64,
        }

    def spec(self) -> dict[str, object]:
        return {
            "schema": binding.SPEC_SCHEMA,
            "immutable_authority": {
                "bundle_manifest": self.artifact("bundle.json", "1"),
                "materialization_receipt": self.artifact(
                    "materialization.json", "2"
                ),
            },
            "preflight": {
                "plan": self.artifact("plan.json", "3"),
                "preflight_receipt": self.artifact("receipt.json", "4"),
                "source_manifest": self.artifact("source.json", "5"),
                "rows": self.artifact("rows.jsonl", "6"),
                "duplicate_fingerprint": self.artifact(
                    "duplicate.txt", "7"
                ),
                "partition": self.artifact("partition.jsonl", "8"),
            },
            "capacity": {
                "resource_certificate": self.artifact("capacity.json", "9"),
                "root_join_certificate": self.artifact("root_join.json", "a"),
                "pp_audit": self.artifact("pp_audit.json", "b"),
                "auau_audit": self.artifact("auau_audit.json", "c"),
            },
        }

    def normalized_inputs(self):
        spec = self.spec()
        immutable = {
            "status": "PASS",
            "bundle_manifest": {
                **spec["immutable_authority"]["bundle_manifest"],
                "schema": binding.resolver.BUNDLE_SCHEMA,
            },
            "materialization_receipt": {
                **spec["immutable_authority"]["materialization_receipt"],
                "schema": binding.resolver.MATERIALIZATION_SCHEMA,
            },
        }
        preflight = {
            "status": "PASS",
            **{
                field: {**value}
                for field, value in spec["preflight"].items()
            },
            "execution_partition_sha256": "d" * 64,
        }
        preflight["partition"]["schema"] = "jsonl"
        capacity_evidence = {
            "resource_certificate": {
                **spec["capacity"]["resource_certificate"],
                "schema": binding.evidence.CAPACITY_SCHEMA,
            },
            "root_join_certificate": {
                **spec["capacity"]["root_join_certificate"],
                "schema": binding.evidence.ROOT_JOIN_SCHEMA,
            },
            "pp_audit": {
                **spec["capacity"]["pp_audit"],
                "schema": binding.evidence.CAPACITY_AUDIT_SCHEMA,
            },
            "auau_audit": {
                **spec["capacity"]["auau_audit"],
                "schema": binding.evidence.CAPACITY_AUDIT_SCHEMA,
            },
        }
        selected_rows = [
            self.capacity_row("pp_background_jet8", "pp"),
            self.capacity_row("auau_background_jet12", "auau"),
        ]
        return spec, immutable, preflight, capacity_evidence, selected_rows

    @staticmethod
    def capacity_row(row_id: str, system: str) -> dict[str, object]:
        return {
            "row_id": row_id,
            "system": system,
            "population_state": (
                "VALID_EMPTY" if system == "pp" else "POPULATED"
            ),
            "remote_wall_clock_seconds": 100,
            "memory_usage_mb": 500,
            "analysis_output": {
                "path": f"/capacity/{row_id}.root",
                "sha256": "e" * 64,
                "size_bytes": 60_000,
            },
            "sidecar_output": {
                "path": f"/capacity/{row_id}_sidecar.root",
                "sha256": "f" * 64,
                "size_bytes": 70_000,
            },
        }

    def assemble(self) -> dict[str, object]:
        spec, immutable, preflight, capacity_evidence, selected_rows = (
            self.normalized_inputs()
        )
        with (
            mock.patch.object(
                binding.evidence,
                "validate_immutable_authority",
                return_value=immutable,
            ),
            mock.patch.object(
                binding,
                "_load_current_preflight",
                return_value=(
                    preflight,
                    [{"row_id": "source"}],
                    dict(binding.resolver.SIDECAR_ONLY_ARTIFACT_PROFILE),
                ),
            ),
            mock.patch.object(
                binding,
                "_load_capacity",
                return_value=(capacity_evidence, selected_rows),
            ),
        ):
            return binding.build_binding(spec)

    def test_current_binding_has_exact_counts_and_no_legacy_authority(self) -> None:
        payload = self.assemble()
        self.assertEqual(payload["schema"], binding.BINDING_SCHEMA)
        self.assertEqual(payload["status"], "PASS")
        self.assertEqual(
            payload["authority_state"], binding.AUTHORITY_STATE
        )
        self.assertNotIn("legacy_preflight", payload)
        self.assertNotIn("source_semantic_crosswalk", payload)
        self.assertFalse(payload["checks"]["legacy_evidence_consumed"])
        self.assertEqual(
            payload["count_contract"],
            {
                "basis": (
                    "ordered_disjoint_group_of_seven_partition_one_execution_"
                    "and_source_occurrence_per_output_pair"
                ),
                "row_count": 13,
                "source_tuple_count": 129_998,
                "group_size": 7,
                "expected_chunk_count": 18_577,
                "expected_job_count": 18_577,
                "expected_output_pair_count": 18_577,
                "expected_analysis_output_count": 18_577,
                "expected_sidecar_output_count": 18_577,
                "expected_physical_root_artifact_count": 37_154,
                "expected_source_occurrence_count": 18_577,
                "source_occurrences_per_output_pair": 1,
                "request_memory_mb": 8_000,
            },
        )
        self.assertFalse(payload["submission_performed"])
        self.assertEqual(payload["full_training_authority"], 0)
        self.assertFalse(payload["full_extraction_authority"])
        self.assertFalse(payload["broad_production_authority"])
        self.assertFalse(payload["canonical_promotion"])

    def test_capacity_validator_binds_resource_and_chunks_to_same_current_plan(
        self,
    ) -> None:
        current = {"plan": {"sha256": "1" * 64}}
        immutable = {"bundle_manifest": {"sha256": "2" * 64}}
        artifacts = [
            ({}, {"path": "/resource", "sha256": "3" * 64}),
            (
                {
                    **{
                        key: None
                        for key in binding.evidence.ROOT_JOIN_KEYS
                    },
                    "artifact_profile": dict(
                        binding.resolver.SIDECAR_ONLY_ARTIFACT_PROFILE
                    ),
                },
                {"path": "/root", "sha256": "4" * 64},
            ),
            ({}, {"path": "/pp", "sha256": "5" * 64}),
            ({}, {"path": "/auau", "sha256": "6" * 64}),
        ]
        with (
            mock.patch.object(
                binding.evidence,
                "load_json_artifact",
                side_effect=artifacts,
            ),
            mock.patch.object(
                binding.evidence,
                "validate_capacity_evidence",
                return_value=({}, []),
            ) as validate,
        ):
            binding._load_capacity(
                self.spec()["capacity"],
                source_records=[],
                current=current,
                immutable=immutable,
                artifact_profile=dict(
                    binding.resolver.SIDECAR_ONLY_ARTIFACT_PROFILE
                ),
            )
        args = validate.call_args.args
        self.assertNotIn("artifact_profile", args[1])
        self.assertIs(args[5], current)
        self.assertIs(args[6], current)
        self.assertNotIn("legacy", repr(validate.call_args))

    def test_capacity_validator_rejects_artifact_profile_drift(self) -> None:
        current = {"plan": {"sha256": "1" * 64}}
        immutable = {"bundle_manifest": {"sha256": "2" * 64}}
        drifted_profile = dict(
            binding.resolver.SIDECAR_ONLY_ARTIFACT_PROFILE
        )
        drifted_profile["replay_serialization"] = "ENABLED"
        root_join = {
            **{key: None for key in binding.evidence.ROOT_JOIN_KEYS},
            "artifact_profile": drifted_profile,
        }
        artifacts = [
            ({}, {"path": "/resource", "sha256": "3" * 64}),
            (root_join, {"path": "/root", "sha256": "4" * 64}),
            ({}, {"path": "/pp", "sha256": "5" * 64}),
            ({}, {"path": "/auau", "sha256": "6" * 64}),
        ]
        with (
            mock.patch.object(
                binding.evidence,
                "load_json_artifact",
                side_effect=artifacts,
            ),
            self.assertRaisesRegex(
                binding.BindingError,
                "artifact profile differs from current plan",
            ),
        ):
            binding._load_capacity(
                self.spec()["capacity"],
                source_records=[],
                current=current,
                immutable=immutable,
                artifact_profile=dict(
                    binding.resolver.SIDECAR_ONLY_ARTIFACT_PROFILE
                ),
            )

    def test_reconstruct_spec_contains_only_current_pinned_inputs(self) -> None:
        payload = self.assemble()
        reconstructed = binding.reconstruct_spec(payload)
        self.assertEqual(reconstructed, self.spec())
        self.assertNotIn("legacy", reconstructed)

    def test_spec_rejects_duplicate_keys_and_nonfinite_numbers(self) -> None:
        cases = {
            "duplicate_schema": (
                '{"schema":"WRONG","schema":"'
                + binding.SPEC_SCHEMA
                + '","immutable_authority":{},'
                '"preflight":{},"capacity":{}}\n'
            ),
            "nonfinite": (
                '{"schema":"'
                + binding.SPEC_SCHEMA
                + '","immutable_authority":{},'
                '"preflight":{},"capacity":{"value":NaN}}\n'
            ),
        }
        for name, text in cases.items():
            path = self.root / f"{name}.json"
            path.write_text(text, encoding="utf-8")
            with self.subTest(name=name), self.assertRaisesRegex(
                binding.BindingError, "ambiguous JSON"
            ):
                binding.load_spec(path, sha256_bytes(path.read_bytes()))

    def test_validation_rebuilds_all_pinned_inputs_and_rejects_tampering(
        self,
    ) -> None:
        payload = self.assemble()
        with mock.patch.object(
            binding, "_assemble_binding", return_value=copy.deepcopy(payload)
        ) as rebuild:
            observed = binding.validate_binding_payload(payload)
        self.assertEqual(observed, payload)
        rebuild.assert_called_once_with(binding.reconstruct_spec(payload))

        for mutation in ("count", "authority", "hidden_legacy"):
            changed = copy.deepcopy(payload)
            if mutation == "count":
                changed["count_contract"]["expected_job_count"] += 1
            elif mutation == "authority":
                changed["full_extraction_authority"] = True
            else:
                changed["legacy_preflight"] = {"status": "PASS"}
            with self.subTest(mutation=mutation):
                with self.assertRaises(binding.BindingError):
                    binding.validate_binding_payload(changed)

    def write_json(self, name: str, payload: dict[str, object]) -> dict[str, str]:
        path = self.root / name
        path.write_text(
            json.dumps(payload, sort_keys=True) + "\n", encoding="utf-8"
        )
        return {"path": str(path), "sha256": sha256_bytes(path.read_bytes())}

    def test_storage_projector_accepts_current_binding_and_readback(self) -> None:
        selected_rows = [
            self.capacity_row("pp_background_jet8", "pp"),
            self.capacity_row("auau_background_jet12", "auau"),
        ]
        current_binding = {
            "schema": binding.BINDING_SCHEMA,
            "binding_semantic_sha256": "d" * 64,
            "capacity_binding": {"selected_rows": selected_rows},
        }
        binding_ref = self.write_json("binding.json", current_binding)
        readback = {
            "schema": binding.READBACK_SCHEMA,
            "status": "PASS",
            "authority_state": binding.AUTHORITY_STATE,
            "binding": binding_ref["path"],
            "binding_sha256": binding_ref["sha256"],
            "binding_semantic_sha256": "d" * 64,
            "all_evidence_rehashed": True,
            "byte_exact_rebuild": True,
            "legacy_evidence_consumed": False,
            "submission_performed": False,
            "full_training_authority": 0,
            "full_extraction_authority": False,
        }
        readback_ref = self.write_json("readback.json", readback)
        with mock.patch.object(
            projector.capacity_binding_tool,
            "validate_binding_payload",
            return_value=current_binding,
        ) as validator:
            (
                observed,
                observed_artifact,
                observed_readback,
                witnesses,
            ) = projector.load_revalidated_amendment_chain(
                binding_ref, readback_ref
            )
        self.assertEqual(observed, current_binding)
        self.assertEqual(observed_artifact["sha256"], binding_ref["sha256"])
        self.assertEqual(observed_readback["sha256"], readback_ref["sha256"])
        self.assertEqual(
            set(witnesses),
            {"pp_background_jet8", "auau_background_jet12"},
        )
        self.assertEqual(
            witnesses["pp_background_jet8"]["request_memory_mb"], 8_000
        )
        validator.assert_called_once_with(current_binding)

    def test_storage_projector_rejects_readback_claiming_legacy_consumption(
        self,
    ) -> None:
        current_binding = {
            "schema": binding.BINDING_SCHEMA,
            "binding_semantic_sha256": "d" * 64,
            "capacity_binding": {"selected_rows": []},
        }
        binding_ref = self.write_json("binding.json", current_binding)
        readback = {
            "schema": binding.READBACK_SCHEMA,
            "status": "PASS",
            "authority_state": binding.AUTHORITY_STATE,
            "binding": binding_ref["path"],
            "binding_sha256": binding_ref["sha256"],
            "binding_semantic_sha256": "d" * 64,
            "all_evidence_rehashed": True,
            "byte_exact_rebuild": True,
            "legacy_evidence_consumed": True,
            "submission_performed": False,
            "full_training_authority": 0,
            "full_extraction_authority": False,
        }
        readback_ref = self.write_json("readback.json", readback)
        with mock.patch.object(
            projector.capacity_binding_tool,
            "validate_binding_payload",
            return_value=current_binding,
        ):
            with self.assertRaisesRegex(
                projector.ProjectionError,
                "capacity partition binding readback differs",
            ):
                projector.load_revalidated_amendment_chain(
                    binding_ref, readback_ref
                )

    def test_storage_projector_validates_current_binding_plan_and_counts(
        self,
    ) -> None:
        row_ids = projector.expected_row_ids()
        rows = []
        for index, row_id in enumerate(row_ids):
            source_count = 10_000 if index < 12 else 9_998
            jobs = (source_count + projector.EXPECTED_GROUP_SIZE - 1) // (
                projector.EXPECTED_GROUP_SIZE
            )
            rows.append(
                {
                    "row_id": row_id,
                    "system": "pp" if row_id.startswith("pp_") else "auau",
                    "input_contract": {
                        "source_tuple_count": source_count,
                        "group_size": projector.EXPECTED_GROUP_SIZE,
                        "expected_chunk_count": jobs,
                        "expected_job_count": jobs,
                        "expected_output_pair_count": jobs,
                        "expected_analysis_output_count": jobs,
                        "expected_sidecar_output_count": jobs,
                        "expected_source_occurrence_count": jobs,
                        "source_occurrences_per_output_pair": 1,
                    },
                    "execution_contract": {
                        "worker_environment": {
                            "RJ_THE134_EPHEMERAL_ANALYSIS_OUTPUT": "1"
                        },
                        "materialization_environment": {
                            "RJ_REQUEST_MEMORY": "8000MB"
                        }
                    },
                }
            )
        plan = {
            "schema": projector.resolver.PLAN_SCHEMA,
            "status": "PREFLIGHT_PASS",
            "submission_performed": False,
            "rows": rows,
            "execution_partition": {
                "schema": projector.resolver.PARTITION_SCHEMA,
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
            },
        }
        plan_artifact = self.artifact("plan.json", "1")
        receipt_artifact = self.artifact("receipt.json", "2")
        binding_artifact = self.artifact("binding.json", "3")
        readback_artifact = self.artifact("readback.json", "4")
        bundle_artifact = self.artifact("bundle.json", "5")
        materialization_artifact = self.artifact("materialization.json", "6")
        current = {
            "plan": plan_artifact,
            "preflight_receipt": receipt_artifact,
            "partition": {
                "path": "/partition.jsonl",
                "sha256": "8" * 64,
            },
            "execution_partition_sha256": "7" * 64,
        }
        current_binding = {
            "schema": binding.BINDING_SCHEMA,
            "preflight": current,
            "immutable_authority": {
                "bundle_manifest": bundle_artifact,
                "materialization_receipt": materialization_artifact,
            },
            "count_contract": {
                "row_count": 13,
                "group_size": 7,
                "source_tuple_count": 129_998,
                "expected_chunk_count": 18_577,
                "expected_job_count": 18_577,
                "expected_output_pair_count": 18_577,
                "expected_analysis_output_count": 18_577,
                "expected_sidecar_output_count": 18_577,
                "expected_physical_root_artifact_count": 37_154,
                "expected_source_occurrence_count": 18_577,
                "source_occurrences_per_output_pair": 1,
                "request_memory_mb": 8_000,
            },
            "capacity_binding": {
                "execution_partition_sha256": "7" * 64,
                "partition_artifact_sha256": "8" * 64,
            },
        }
        spec = {
            "bindings": {
                "plan": plan_artifact,
                "preflight_receipt": receipt_artifact,
                "capacity_count_amendment": binding_artifact,
                "capacity_count_amendment_readback": readback_artifact,
                "bundle_manifest": bundle_artifact,
                "materialization_receipt": materialization_artifact,
            }
        }
        artifact_loads = [
            (plan, plan_artifact),
            ({"schema": projector.resolver.RECEIPT_SCHEMA}, receipt_artifact),
            ({"schema": "bundle"}, bundle_artifact),
            ({"schema": "materialization"}, materialization_artifact),
        ]
        with (
            mock.patch.object(
                projector, "load_artifact", side_effect=artifact_loads
            ),
            mock.patch.object(
                projector, "validate_preflight_receipt", return_value=None
            ),
            mock.patch.object(
                projector,
                "load_revalidated_amendment_chain",
                return_value=(
                    current_binding,
                    binding_artifact,
                    readback_artifact,
                    {
                        "pp_background_jet8": {},
                        "auau_background_jet12": {},
                    },
                ),
            ),
        ):
            _plan, normalized_rows, _witnesses, bindings_out, blockers = (
                projector.validate_chain(spec)
            )
        self.assertEqual(len(normalized_rows), 13)
        self.assertEqual(blockers, [])
        self.assertEqual(
            bindings_out["execution_partition_sha256"], "7" * 64
        )
        self.assertEqual(bindings_out["partition_artifact_sha256"], "8" * 64)


if __name__ == "__main__":
    unittest.main()
