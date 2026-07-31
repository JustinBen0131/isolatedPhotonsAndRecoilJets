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
            "analysis_health": {
                "path": f"/capacity/{row_id}.root",
                "artifact_state": "EPHEMERAL_VALIDATED_NOT_RETAINED",
                "health_receipt_path": f"/capacity/{row_id}.stdout",
                "health_receipt_sha256": "e" * 64,
                "key_inventory_sha256": "d" * 64,
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
            mock.patch.object(
                binding.evidence,
                "load_json_artifact",
                return_value=({}, {}),
            ),
            mock.patch.object(
                binding,
                "validate_capacity_plan_binding",
                return_value={"mode": "EXACT_PLAN"},
            ),
            mock.patch.object(
                binding,
                "validate_capacity_preflight_binding",
                return_value={"mode": "EXACT_PREFLIGHT"},
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

    def test_capacity_validator_does_not_use_legacy_durable_root_contract(
        self,
    ) -> None:
        current = {"plan": {"sha256": "1" * 64}}
        immutable = {"bundle_manifest": {"sha256": "2" * 64}}
        artifacts = [
            (
                {
                    **{
                        key: None
                        for key in binding.evidence.CAPACITY_KEYS
                    },
                    "schema": binding.evidence.CAPACITY_SCHEMA,
                    "status": "PASS",
                },
                {"path": "/resource", "sha256": "3" * 64},
            ),
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
            ) as legacy_validate,
            self.assertRaises(binding.BindingError),
        ):
            binding._load_capacity(
                self.spec()["capacity"],
                source_records=[],
                current=current,
                immutable=immutable,
                artifact_profile=dict(
                    binding.resolver.SIDECAR_ONLY_ARTIFACT_PROFILE
                ),
                plan_binding={"mode": "EXACT_PLAN"},
                preflight_binding={"mode": "EXACT_PREFLIGHT"},
            )
        legacy_validate.assert_not_called()

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
                plan_binding={"mode": "EXACT_PLAN"},
                preflight_binding={"mode": "EXACT_PREFLIGHT"},
            )

    def test_validation_only_commit_is_distinct_and_explicit(self) -> None:
        observed = binding.validate_science_validation_commit_binding(
            {"public_commit": "1" * 40},
            {"controller": {"validation_commit": "2" * 40}},
        )
        self.assertEqual(
            observed,
            {
                "science_commit": "1" * 40,
                "validation_commit": "2" * 40,
                "authority_mode": "FROZEN_SCIENCE_WITH_VALIDATION_ONLY_COMMIT",
            },
        )

    def test_validation_commit_binding_rejects_malformed_commits(self) -> None:
        cases = (
            (
                {"public_commit": "1" * 39},
                {"controller": {"validation_commit": "2" * 40}},
                "science public_commit",
            ),
            (
                {"public_commit": "1" * 40},
                {"controller": {"validation_commit": "not-a-commit"}},
                "validation_commit",
            ),
        )
        for immutable, authority, message in cases:
            with self.subTest(message=message), self.assertRaisesRegex(
                binding.BindingError, message
            ):
                binding.validate_science_validation_commit_binding(
                    immutable, authority
                )

    def test_capacity_plan_reuses_only_fresh_namespace_equivalent_plan(
        self,
    ) -> None:
        old_tag = "the134_old_attempt"
        new_tag = "the134_new_attempt"

        def plan(tag: str, threshold: float = 5.0) -> dict[str, object]:
            return {
                "schema": "THE134_FULL_EXTRACTION_PLAN_V3",
                "campaign": {
                    "tag": tag,
                    "output_root": f"/output/{tag}",
                    "evidence_root": f"/evidence/{tag}",
                    "submit_root": f"/submit/{tag}",
                },
                "execution_fingerprint_sha256": "1" * 64,
                "rows": [
                    {
                        "row_id": "pp_background_jet8",
                        "row_fingerprint_sha256": "2" * 64,
                        "execution_contract": {
                            "submit_namespace": f"/submit/{tag}/row",
                            "worker_environment": {
                                "RJ_PROFILE_LABEL": f"{tag}_row",
                                "RJ_REPLAY_PHOTON_CAPTURE_ET_MIN": threshold,
                            },
                        },
                    }
                ],
            }

        old_path = self.root / "old_plan.json"
        new_path = self.root / "new_plan.json"
        old_path.write_text(json.dumps(plan(old_tag)), encoding="utf-8")
        new_path.write_text(json.dumps(plan(new_tag)), encoding="utf-8")
        resource = {
            "full_plan": str(old_path),
            "full_plan_sha256": sha256_bytes(old_path.read_bytes()),
        }
        current = {
            "plan": {
                "path": str(new_path),
                "sha256": sha256_bytes(new_path.read_bytes()),
            }
        }
        observed = binding.validate_capacity_plan_binding(resource, current)
        self.assertEqual(
            observed["mode"],
            "FRESH_NAMESPACE_AND_PINNED_PROVIDER_EQUIVALENT_PLAN",
        )
        self.assertEqual(
            observed["namespace_normalized_sha256"],
            binding.semantic_sha256(
                binding._capacity_plan_namespace_normal_form(plan(new_tag))
            ),
        )

        new_path.write_text(
            json.dumps(plan(new_tag, threshold=6.0)), encoding="utf-8"
        )
        current["plan"]["sha256"] = sha256_bytes(new_path.read_bytes())
        with self.assertRaisesRegex(
            binding.BindingError, "differs beyond campaign namespace"
        ):
            binding.validate_capacity_plan_binding(resource, current)

    def test_capacity_plan_allows_only_verified_operational_bundle_refresh(
        self,
    ) -> None:
        old_tag = "the134_old_bundle"
        new_tag = "the134_new_bundle"

        def plan(
            tag: str,
            bundle_id: str,
            code_sha: str,
            submitter_sha: str,
            submitter_size: int,
        ) -> dict[str, object]:
            root = (
                f"/immutable_bundles/the134_bundle_sha256_{bundle_id}"
            )
            return {
                "schema": "THE134_FULL_EXTRACTION_PLAN_V3",
                "campaign": {
                    "tag": tag,
                    "output_root": f"/output/{tag}",
                    "evidence_root": f"/evidence/{tag}",
                    "submit_root": f"/submit/{tag}",
                },
                "duplicate_fingerprint_sha256": "1" * 64,
                "execution_fingerprint_sha256": "2" * 64,
                "input_manifests": {
                    "bundle": {
                        "path": f"{root}/metadata/resolver_bundle_receipt.json",
                        "sha256": bundle_id,
                        "semantic_fingerprint_sha256": bundle_id,
                    },
                    "materialization": {
                        "path": f"{root}/metadata/materialization_receipt.json",
                        "sha256": bundle_id,
                        "bundle_identity_sha256": bundle_id,
                        "digest_named_bundle_path": root,
                    },
                },
                "duplicate_contract": {
                    "bundle_semantic_fingerprint_sha256": bundle_id,
                },
                "rows": [
                    {
                        "row_id": "pp_background_jet8",
                        "row_fingerprint_sha256": "3" * 64,
                        "bundle_contract": {
                            "code_sha256": code_sha,
                            "library": {
                                "path": f"{root}/artifacts/pp_library/lib.so",
                                "sha256": "a" * 64,
                            },
                            "submitter": {
                                "path": (
                                    f"{root}/artifacts/submitter/submit.sh"
                                ),
                                "sha256": submitter_sha,
                                "size_bytes": submitter_size,
                            },
                        },
                        "execution_contract": {
                            "submit_namespace": f"/submit/{tag}/row",
                            "worker_environment": {
                                "RJ_PROFILE_LABEL": f"{tag}_row",
                                "RJ_REPLAY_CODE_SHA256": code_sha,
                                "RJ_REPLAY_PHOTON_CAPTURE_ET_MIN": 5.0,
                            },
                        },
                    }
                ],
            }

        old_plan = plan(
            old_tag, "4" * 64, "5" * 64, "6" * 64, 100
        )
        new_plan = plan(
            new_tag, "7" * 64, "8" * 64, "9" * 64, 200
        )
        authorities = {
            old_tag: {
                "bundle_roots": [
                    f"/immutable_bundles/the134_bundle_sha256_{'4' * 64}"
                ],
                "bundle_identity_sha256": "4" * 64,
                "bundle_file_sha256": "4" * 64,
                "bundle_semantic_fingerprint_sha256": "4" * 64,
                "materialization_file_sha256": "4" * 64,
                "code_sha256": "5" * 64,
                "submitter_sha256": "6" * 64,
                "submitter_size_bytes": 100,
                "science_bundle_fingerprint_sha256": "a" * 64,
            },
            new_tag: {
                "bundle_roots": [
                    f"/immutable_bundles/the134_bundle_sha256_{'7' * 64}"
                ],
                "bundle_identity_sha256": "7" * 64,
                "bundle_file_sha256": "7" * 64,
                "bundle_semantic_fingerprint_sha256": "7" * 64,
                "materialization_file_sha256": "7" * 64,
                "code_sha256": "8" * 64,
                "submitter_sha256": "9" * 64,
                "submitter_size_bytes": 200,
                "science_bundle_fingerprint_sha256": "a" * 64,
            },
        }

        def authority(payload):
            return authorities[payload["campaign"]["tag"]]

        old_path = self.root / "old_bundle_plan.json"
        new_path = self.root / "new_bundle_plan.json"
        old_path.write_text(json.dumps(old_plan), encoding="utf-8")
        new_path.write_text(json.dumps(new_plan), encoding="utf-8")
        resource = {
            "full_plan": str(old_path),
            "full_plan_sha256": sha256_bytes(old_path.read_bytes()),
        }
        current = {
            "plan": {
                "path": str(new_path),
                "sha256": sha256_bytes(new_path.read_bytes()),
            }
        }
        with mock.patch.object(
            binding,
            "_validated_plan_bundle_operational_authority",
            side_effect=authority,
        ):
            observed = binding.validate_capacity_plan_binding(
                resource, current
            )
        self.assertEqual(
            observed["mode"],
            "FRESH_NAMESPACE_AND_PINNED_PROVIDER_EQUIVALENT_PLAN",
        )

        authorities[new_tag][
            "science_bundle_fingerprint_sha256"
        ] = "b" * 64
        with mock.patch.object(
            binding,
            "_validated_plan_bundle_operational_authority",
            side_effect=authority,
        ), self.assertRaisesRegex(
            binding.BindingError, "differs beyond campaign namespace"
        ):
            binding.validate_capacity_plan_binding(resource, current)

    def test_bundle_validation_preserves_materialization_alias_path(
        self,
    ) -> None:
        actual_root = self.root / "gpfs"
        alias_root = self.root / "sphenix_u"
        bundle_name = f"the134_bundle_sha256_{'a' * 64}"
        metadata = actual_root / bundle_name / "metadata"
        metadata.mkdir(parents=True)
        alias_root.symlink_to(actual_root, target_is_directory=True)
        bundle_path = metadata / "resolver_bundle_receipt.json"
        materialization_path = metadata / "materialization_receipt.json"
        bundle_path.write_text("{}\n", encoding="utf-8")
        materialization_path.write_text("{}\n", encoding="utf-8")
        alias_bundle_path = (
            alias_root / bundle_name / "metadata" / bundle_path.name
        )
        alias_materialization_path = (
            alias_root / bundle_name / "metadata" / materialization_path.name
        )
        plan = {
            "input_manifests": {
                "bundle": {
                    "path": str(alias_bundle_path),
                    "sha256": sha256_bytes(bundle_path.read_bytes()),
                },
                "materialization": {
                    "path": str(alias_materialization_path),
                    "sha256": sha256_bytes(materialization_path.read_bytes()),
                },
            }
        }
        validated_bundle = {
            "artifact_by_role": {
                "submitter": {"sha256": "b" * 64, "size_bytes": 10},
                "pp_library": {"sha256": "c" * 64, "size_bytes": 20},
            },
            "bundle_identity_sha256": "a" * 64,
            "semantic_fingerprint_sha256": "d" * 64,
            "code_sha256": "e" * 64,
            "public_commit": "f" * 40,
            "replay_schema_sha256": "1" * 64,
            "training_schema_sha256": "2" * 64,
            "semantic_sha256": "3" * 64,
            "runtime": {
                "release": "ana.561",
                "offline_main": "new",
                "calo_reco_soname": "libcalo_reco.so",
                "request_memory_mb": 8000,
            },
        }

        def validate_materialization(payload, **kwargs):
            self.assertEqual(
                kwargs["materialization_path"], alias_materialization_path
            )
            self.assertEqual(kwargs["bundle_path"], bundle_path)
            return {"digest_named_bundle_path": str(alias_root / bundle_name)}

        with mock.patch.object(
            binding.resolver,
            "validate_bundle",
            return_value=validated_bundle,
        ), mock.patch.object(
            binding.resolver,
            "validate_materialization_binding",
            side_effect=validate_materialization,
        ):
            observed = binding._validated_plan_bundle_operational_authority(
                plan
            )
        self.assertEqual(observed["bundle_identity_sha256"], "a" * 64)
        self.assertEqual(
            observed["materialization_file_sha256"],
            sha256_bytes(materialization_path.read_bytes()),
        )

    def test_capacity_plan_accepts_only_byte_identical_pinned_provider_routing(
        self,
    ) -> None:
        old_tag = "the134_old_provider_attempt"
        new_tag = "the134_new_provider_attempt"
        release_lib = self.root / "release" / "lib"
        release_lib64 = self.root / "release" / "lib64"
        release_lib.mkdir(parents=True)
        release_lib64.mkdir(parents=True)
        bundle_root = (
            self.root
            / "immutable_bundles"
            / f"the134_bundle_sha256_{'a' * 64}"
            / "artifacts"
        )
        provider_rows = (
            ("libcalo_io.so", "release_calo_io", "CALO_IO"),
            ("libclusteriso.so", "release_clusteriso", "CLUSTERISO"),
            ("libjetbase.so", "release_jetbase", "JETBASE"),
        )
        old_paths = {}
        new_paths = {}
        hashes = {}
        for family, role, _suffix in provider_rows:
            release_path = release_lib / family
            release_path.write_text(f"{family}|provider\n", encoding="utf-8")
            bundle_path = bundle_root / role / family
            bundle_path.parent.mkdir(parents=True)
            bundle_path.write_bytes(release_path.read_bytes())
            old_paths[family] = str(bundle_path)
            new_paths[family] = str(release_path)
            hashes[family] = binding.file_sha256(release_path)

        def plan(
            tag: str,
            paths: dict[str, str],
            *,
            mirrored_submitter_bindings: bool = False,
        ) -> dict[str, object]:
            worker_environment = {
                "RJ_PROFILE_LABEL": f"{tag}_row",
                "RJ_THE134_MULTIVIEW_TRAINING_V1": "1",
                "RJ_THE134_MULTIVIEW_SIDECAR_ONLY_V1": "1",
                "RJ_THE134_EPHEMERAL_ANALYSIS_OUTPUT": "1",
            }
            environment = {
                "RJ_RELEASE_CORE_LIB_DIR": str(release_lib),
                "RJ_RELEASE_CORE_LIB64_DIR": str(release_lib64),
            }
            if mirrored_submitter_bindings:
                environment.update(worker_environment)
            for family, _role, suffix in provider_rows:
                environment[f"RJ_PINNED_RELEASE_{suffix}_PATH"] = paths[family]
                environment[f"RJ_PINNED_RELEASE_{suffix}_SHA256"] = hashes[family]
            return {
                "schema": "THE134_FULL_EXTRACTION_PLAN_V3",
                "campaign": {
                    "tag": tag,
                    "output_root": f"/output/{tag}",
                    "evidence_root": f"/evidence/{tag}",
                    "submit_root": f"/submit/{tag}",
                },
                "execution_fingerprint_sha256": "1" * 64,
                "rows": [
                    {
                        "row_id": "pp_background_jet8",
                        "row_fingerprint_sha256": "2" * 64,
                        "execution_contract": {
                            "submit_namespace": f"/submit/{tag}/row",
                            "materialization_environment": environment,
                            "worker_environment": worker_environment,
                        },
                    }
                ],
            }

        old_path = self.root / "old_provider_plan.json"
        new_path = self.root / "new_provider_plan.json"
        old_path.write_text(json.dumps(plan(old_tag, old_paths)), encoding="utf-8")
        new_path.write_text(
            json.dumps(
                plan(
                    new_tag,
                    new_paths,
                    mirrored_submitter_bindings=True,
                )
            ),
            encoding="utf-8",
        )
        resource = {
            "full_plan": str(old_path),
            "full_plan_sha256": binding.file_sha256(old_path),
        }
        current = {
            "plan": {
                "path": str(new_path),
                "sha256": binding.file_sha256(new_path),
            }
        }
        observed = binding.validate_capacity_plan_binding(resource, current)
        self.assertEqual(
            observed["mode"],
            "FRESH_NAMESPACE_AND_PINNED_PROVIDER_EQUIVALENT_PLAN",
        )

        mismatched = plan(
            new_tag,
            new_paths,
            mirrored_submitter_bindings=True,
        )
        mismatched["rows"][0]["execution_contract"][
            "materialization_environment"
        ]["RJ_THE134_MULTIVIEW_TRAINING_V1"] = "0"
        new_path.write_text(json.dumps(mismatched), encoding="utf-8")
        current["plan"]["sha256"] = binding.file_sha256(new_path)
        with self.assertRaisesRegex(
            binding.BindingError, "submitter mirror differs"
        ):
            binding.validate_capacity_plan_binding(resource, current)
        new_path.write_text(
            json.dumps(
                plan(
                    new_tag,
                    new_paths,
                    mirrored_submitter_bindings=True,
                )
            ),
            encoding="utf-8",
        )
        current["plan"]["sha256"] = binding.file_sha256(new_path)

        bundle_authorities = {
            old_tag: {
                "bundle_roots": [str(bundle_root.parent)],
                "bundle_identity_sha256": "a" * 64,
                "bundle_file_sha256": "b" * 64,
                "bundle_semantic_fingerprint_sha256": "c" * 64,
                "materialization_file_sha256": "d" * 64,
                "code_sha256": "e" * 64,
                "submitter_sha256": "f" * 64,
                "submitter_size_bytes": 100,
                "science_bundle_fingerprint_sha256": "1" * 64,
            },
            new_tag: {
                "bundle_roots": [str(self.root / "new_bundle")],
                "bundle_identity_sha256": "2" * 64,
                "bundle_file_sha256": "3" * 64,
                "bundle_semantic_fingerprint_sha256": "4" * 64,
                "materialization_file_sha256": "5" * 64,
                "code_sha256": "6" * 64,
                "submitter_sha256": "7" * 64,
                "submitter_size_bytes": 200,
                "science_bundle_fingerprint_sha256": "1" * 64,
            },
        }

        def authority(payload):
            return bundle_authorities[payload["campaign"]["tag"]]

        with mock.patch.object(
            binding,
            "_validated_plan_bundle_operational_authority",
            side_effect=authority,
        ):
            observed = binding.validate_capacity_plan_binding(
                resource, current
            )
        self.assertEqual(
            observed["mode"],
            "FRESH_NAMESPACE_AND_PINNED_PROVIDER_EQUIVALENT_PLAN",
        )

        (release_lib / "libclusteriso.so").write_text(
            "hash drift\n", encoding="utf-8"
        )
        with self.assertRaisesRegex(
            binding.BindingError, "pinned provider hash differs"
        ):
            binding.validate_capacity_plan_binding(resource, current)

    def test_capacity_preflight_reuses_only_derived_identity_changes(
        self,
    ) -> None:
        old_plan = self.root / "old" / "plan.json"
        old_plan.parent.mkdir()
        old_receipt = old_plan.parent / "preflight_receipt.json"
        new_receipt = self.root / "new_receipt.json"

        def receipt(plan_sha: str, rows_sha: str, partition_sha: str):
            return {
                "schema": "THE134_FULL_EXTRACTION_PREFLIGHT_V3",
                "execution_fingerprint_sha256": "1" * 64,
                "duplicate_fingerprint_sha256": "2" * 64,
                "bundle_manifest_sha256": "3" * 64,
                "materialization_receipt_sha256": "4" * 64,
                "artifacts": {
                    "plan": {"sha256": plan_sha},
                    "rows": {"sha256": rows_sha},
                    "partition": {"sha256": partition_sha},
                    "duplicate_fingerprint": {"sha256": "5" * 64},
                },
            }

        old_receipt.write_text(
            json.dumps(receipt("2" * 64, "3" * 64, "4" * 64)),
            encoding="utf-8",
        )
        new_receipt.write_text(
            json.dumps(receipt("5" * 64, "6" * 64, "4" * 64)),
            encoding="utf-8",
        )
        resource = {
            "full_plan": str(old_plan),
            "preflight_receipt_sha256": sha256_bytes(
                old_receipt.read_bytes()
            ),
        }
        current = {
            "preflight_receipt": {
                "path": str(new_receipt),
                "sha256": sha256_bytes(new_receipt.read_bytes()),
            }
        }
        current_payload = json.loads(new_receipt.read_text(encoding="utf-8"))
        current_payload["duplicate_fingerprint_sha256"] = "6" * 64
        current_payload["bundle_manifest_sha256"] = "7" * 64
        current_payload["materialization_receipt_sha256"] = "8" * 64
        current_payload["artifacts"]["duplicate_fingerprint"]["sha256"] = (
            "9" * 64
        )
        new_receipt.write_text(
            json.dumps(current_payload), encoding="utf-8"
        )
        current["preflight_receipt"]["sha256"] = sha256_bytes(
            new_receipt.read_bytes()
        )
        observed = binding.validate_capacity_preflight_binding(
            resource, current
        )
        self.assertEqual(
            observed["mode"],
            "FRESH_NAMESPACE_ONLY_EQUIVALENT_PREFLIGHT",
        )

        new_receipt.write_text(
            json.dumps(receipt("5" * 64, "6" * 64, "7" * 64)),
            encoding="utf-8",
        )
        current["preflight_receipt"]["sha256"] = sha256_bytes(
            new_receipt.read_bytes()
        )
        with self.assertRaisesRegex(
            binding.BindingError, "namespace-derived identities"
        ):
            binding.validate_capacity_preflight_binding(resource, current)

    def test_resource_bundle_binding_follows_validated_plan_mode(self) -> None:
        resource = {
            "bundle_manifest_sha256": "1" * 64,
            "materialization_receipt_sha256": "2" * 64,
        }
        immutable = {
            "bundle_manifest": {"sha256": "3" * 64},
            "materialization_receipt": {"sha256": "4" * 64},
        }
        self.assertFalse(
            binding.capacity_resource_bundle_binding_valid(
                resource, immutable, {"mode": "EXACT_PLAN"}
            )
        )
        self.assertTrue(
            binding.capacity_resource_bundle_binding_valid(
                resource,
                immutable,
                {
                    "mode": (
                        "FRESH_NAMESPACE_AND_PINNED_PROVIDER_EQUIVALENT_PLAN"
                    )
                },
            )
        )
        self.assertFalse(
            binding.capacity_resource_bundle_binding_valid(
                resource, immutable, {"mode": "UNKNOWN"}
            )
        )
        resource["bundle_manifest_sha256"] = "3" * 64
        resource["materialization_receipt_sha256"] = "4" * 64
        self.assertTrue(
            binding.capacity_resource_bundle_binding_valid(
                resource, immutable, {"mode": "EXACT_PLAN"}
            )
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
