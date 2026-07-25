#!/usr/bin/env python3
"""Adversarial tests for THE-134 pre-extraction storage admission."""

from __future__ import annotations

import copy
import hashlib
import importlib.util
import json
import subprocess
import tempfile
import unittest
from decimal import Decimal
from pathlib import Path
from unittest import mock


HERE = Path(__file__).resolve().parent
MODULE_PATH = HERE.parent / "project_the134_preextraction_storage_quota.py"
SPEC = importlib.util.spec_from_file_location(
    "the134_preextraction_storage", MODULE_PATH
)
assert SPEC is not None and SPEC.loader is not None
projector = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(projector)


def seal(payload: dict, field: str) -> dict:
    payload[field] = projector.semantic_sha256(payload)
    return payload


class PreExtractionStorageTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)
        self.amendment_validator_patch = mock.patch.object(
            projector.amendment_tool,
            "validate_amendment_payload",
            side_effect=lambda payload: dict(payload),
        )
        self.amendment_validator = self.amendment_validator_patch.start()
        bulk_parent = self.root / "bulk"
        submit_parent = self.root / "submit"
        evidence_parent = self.root / "evidence"
        for parent in (bulk_parent, submit_parent, evidence_parent):
            parent.mkdir()
        self.rows = []
        row_ids = projector.expected_row_ids()
        for index, row_id in enumerate(row_ids):
            self.rows.append(
                {
                    "row_id": row_id,
                    "system": "pp" if row_id.startswith("pp_") else "auau",
                    "source_tuple_count": 10_000 if index < 12 else 9_998,
                    "expected_job_count": 1_429 if index < 12 else 1_429,
                }
            )
        # The real distribution has twelve 10,000 rows and one 9,998 row,
        # but per-row chunk counts sum to 18,577 (not 12*1429+1429).
        remaining = projector.EXPECTED_JOBS
        for row in self.rows:
            jobs = (
                row["source_tuple_count"] + projector.EXPECTED_GROUP_SIZE - 1
            ) // projector.EXPECTED_GROUP_SIZE
            row["expected_job_count"] = jobs
            remaining -= jobs
        self.assertEqual(remaining, 0)
        self.witnesses = {
            "pp_background_jet8": self.witness(
                "pp_background_jet8", "pp", 470_880_550, 11_575
            ),
            "auau_background_jet12": self.witness(
                "auau_background_jet12", "auau", 563_791_091, 25_701
            ),
        }
        self.plan = {
            "campaign": {
                "output_root": str(bulk_parent / "the134"),
                "submit_root": str(submit_parent / "the134"),
                "evidence_root": str(evidence_parent / "the134"),
            }
        }
        amendment_binding, readback_binding = self.write_capacity_chain()
        self.bindings = {
            "plan": self.write_binding_artifact("plan"),
            "preflight_receipt": self.write_binding_artifact("preflight"),
            "capacity_count_amendment": amendment_binding,
            "capacity_count_amendment_readback": readback_binding,
            "bundle_manifest": self.write_binding_artifact("bundle"),
            "materialization_receipt": self.write_binding_artifact(
                "materialization"
            ),
            "execution_partition_sha256": "7" * 64,
            "partition_artifact_sha256": "8" * 64,
        }
        self.controller_budget = self.write_controller_budget()

    def tearDown(self) -> None:
        self.amendment_validator_patch.stop()
        self.temporary.cleanup()

    def write_capacity_chain(self) -> tuple[dict, dict]:
        amendment_path = self.root / "amendment.json"
        selected_rows = []
        for row_id in projector.amendment_tool.SELECTED_CAPACITY_ROWS:
            witness = self.witnesses[row_id]
            selected_rows.append(
                {
                    "row_id": row_id,
                    "system": witness["system"],
                    "remote_wall_clock_seconds": witness[
                        "wall_time_seconds"
                    ],
                    "memory_usage_mb": witness["max_rss_mb"],
                    "analysis_output": {
                        "path": witness["analysis_root"]["path"],
                        "sha256": witness["analysis_root"]["sha256"],
                        "size_bytes": witness["analysis_root"]["size_bytes"],
                    },
                    "sidecar_output": {
                        "path": witness["training_sidecar_root"]["path"],
                        "sha256": witness["training_sidecar_root"]["sha256"],
                        "size_bytes": witness[
                            "training_sidecar_root"
                        ]["size_bytes"],
                    },
                }
            )
        amendment = {
            "schema": projector.amendment_tool.AMENDMENT_SCHEMA,
            "status": "PASS",
            "amendment_semantic_sha256": "c" * 64,
            "capacity_to_partition_binding": {
                "selected_rows": selected_rows,
            },
        }
        amendment_path.write_bytes(projector.canonical_json_bytes(amendment))
        amendment_binding = {
            "path": str(amendment_path),
            "sha256": projector.file_sha256(amendment_path),
            "size_bytes": amendment_path.stat().st_size,
        }

        readback_path = self.root / "readback.json"
        readback = {
            "schema": projector.amendment_tool.READBACK_SCHEMA,
            "status": "PASS",
            "authority_state": projector.amendment_tool.AUTHORITY_STATE,
            "amendment": str(amendment_path),
            "amendment_sha256": amendment_binding["sha256"],
            "amendment_semantic_sha256": amendment[
                "amendment_semantic_sha256"
            ],
            "all_evidence_rehashed": True,
            "byte_exact_rebuild": True,
            "submission_performed": False,
            "full_training_authority": 0,
            "full_extraction_authority": False,
        }
        readback_path.write_bytes(projector.canonical_json_bytes(readback))
        readback_binding = {
            "path": str(readback_path),
            "sha256": projector.file_sha256(readback_path),
            "size_bytes": readback_path.stat().st_size,
        }
        return amendment_binding, readback_binding

    def write_binding_artifact(self, label: str) -> dict:
        path = self.root / f"{label}.json"
        path.write_text(
            json.dumps({"fixture": label}, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        return {
            "path": str(path),
            "sha256": projector.file_sha256(path),
            "size_bytes": path.stat().st_size,
        }

    @staticmethod
    def witness(
        row_id: str, system: str, analysis_size: int, sidecar_size: int
    ) -> dict:
        return {
            "witness_id": row_id,
            "row_id": row_id,
            "system": system,
            "tuple_count": 7,
            "request_memory_mb": 8_000,
            "wall_time_seconds": 100,
            "max_rss_mb": 2_000,
            "storage_ceiling_authority": False,
            "analysis_root": {
                "path": f"/capacity/{row_id}/analysis.root",
                "sha256": "a" * 64,
                "size_bytes": analysis_size,
                "health_profile": projector.ANALYSIS_HEALTH_PROFILE,
                "health_status": "PASS",
            },
            "training_sidecar_root": {
                "path": f"/capacity/{row_id}/training.root",
                "sha256": "b" * 64,
                "size_bytes": sidecar_size,
                "health_profile": projector.SIDECAR_HEALTH_PROFILE,
                "health_status": "PASS",
            },
        }

    def write_controller_budget(self) -> dict:
        path = self.root / "controller_budget.json"
        payload = {
            "schema": projector.CONTROLLER_BUDGET_SCHEMA,
            "status": "PASS",
            "submission_performed": False,
            "expected_job_count": projector.EXPECTED_JOBS,
            "storage_budget": {
                "fixed_bytes": 1_000,
                "fixed_inodes": 10,
                "bytes_per_job": 2_000,
                "inodes_per_job": 4,
            },
            "full_training_authority": 0,
            "full_extraction_authority": False,
        }
        path.write_text(json.dumps(payload, sort_keys=True) + "\n")
        return {"path": str(path), "sha256": projector.file_sha256(path)}

    def envelopes(self) -> list[dict]:
        records = []
        for row in self.rows:
            witness_id = (
                "pp_background_jet8"
                if row["system"] == "pp"
                else "auau_background_jet12"
            )
            witness = self.witnesses[witness_id]
            classes = []
            for name, contract in projector.ARTIFACT_CLASS_CONTRACT.items():
                if contract["requires_witness"]:
                    ceiling = projector.ceil_ratio(
                        witness[name]["size_bytes"], 5, 4
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
            records.append({"row_id": row["row_id"], "artifact_classes": classes})
        return records

    def measurement_spec(self) -> dict:
        return {
            "schema": projector.MEASUREMENT_SPEC_SCHEMA,
            "authority": projector.authority_surface(),
            **projector.AUTHORITY_FIELDS,
            "bindings": {},
            "controller_dry_materialization_budget": self.controller_budget,
            "retry_reserve": {"numerator": 1, "denominator": 10},
            "row_envelopes": self.envelopes(),
        }

    def count_plan(self) -> dict:
        rows = []
        for row in self.rows:
            jobs = row["expected_job_count"]
            rows.append(
                {
                    "row_id": row["row_id"],
                    "system": row["system"],
                    "input_contract": {
                        "source_tuple_count": row["source_tuple_count"],
                        "expected_job_count": jobs,
                        "group_size": projector.EXPECTED_GROUP_SIZE,
                        "expected_chunk_count": jobs,
                        "expected_output_pair_count": jobs,
                        "expected_analysis_output_count": jobs,
                        "expected_sidecar_output_count": jobs,
                        "expected_source_occurrence_count": jobs,
                        "source_occurrences_per_output_pair": 1,
                    },
                    "execution_contract": {
                        "materialization_environment": {
                            "RJ_REQUEST_MEMORY": (
                                f"{projector.EXPECTED_REQUEST_MEMORY_MB}MB"
                            )
                        }
                    },
                }
            )
        return {
            "schema": projector.resolver.PLAN_SCHEMA,
            "status": "PREFLIGHT_PASS",
            "submission_performed": False,
            "rows": rows,
            "execution_partition": {
                "schema": projector.resolver.PARTITION_SCHEMA,
                "group_size": projector.EXPECTED_GROUP_SIZE,
                "source_tuple_count": projector.EXPECTED_SOURCE_TUPLES,
                "expected_chunk_count": projector.EXPECTED_JOBS,
                "expected_job_count": projector.EXPECTED_JOBS,
                "expected_output_pair_count": projector.EXPECTED_OUTPUT_PAIRS,
                "expected_analysis_output_count": (
                    projector.EXPECTED_ANALYSIS_OUTPUTS
                ),
                "expected_sidecar_output_count": (
                    projector.EXPECTED_SIDECAR_OUTPUTS
                ),
                "expected_physical_root_artifact_count": (
                    projector.EXPECTED_ROOT_ARTIFACTS
                ),
                "expected_source_occurrence_count": (
                    projector.EXPECTED_SOURCE_OCCURRENCES
                ),
                "source_occurrences_per_output_pair": 1,
            },
        }

    def passed_manifest(self) -> dict:
        return projector.build_storage_manifest_from_validated(
            self.measurement_spec(),
            plan=self.plan,
            rows=self.rows,
            witnesses=self.witnesses,
            bindings=self.bindings,
        )

    def quota_snapshot(
        self,
        manifest: dict,
        *,
        used_bytes: int = 100,
        quota_bytes: int = 20_000_000_000_000,
        used_inodes: int = 100,
        quota_inodes: int = 10_000_000,
        same_quota_domain: bool = True,
    ) -> dict:
        plan_sha = manifest["artifact_measurement"]["bindings"]["plan"]["sha256"]
        plan_path = manifest["artifact_measurement"]["bindings"]["plan"]["path"]
        roots = manifest["artifact_measurement"]["storage_domain_roots"]
        domains = []
        for index, domain in enumerate(projector.STORAGE_DOMAINS):
            filesystem = (
                "/remote"
                if same_quota_domain
                else f"/remote-{index}"
            )
            quota_domain_id = projector.semantic_sha256(
                {"principal": "fixture", "filesystem": filesystem}
            )
            intended_root = roots[domain]
            quota_probe_path = str(
                projector.nearest_existing_parent(Path(intended_root))
            )
            stdout = (
                f"{filesystem} {self.byte_token(used_bytes)} "
                f"{self.byte_token(quota_bytes)} "
                f"{self.byte_token(quota_bytes)} - "
                f"{used_inodes} {quota_inodes} {quota_inodes} -\n"
            )
            stderr = ""
            parsed = projector.parse_lfs_quota_output(stdout, stderr, 0)
            self.assertTrue(parsed["usage_authoritative"])
            domains.append(
                {
                    "storage_domain": domain,
                    "intended_root": intended_root,
                    "namespace_fresh": True,
                    "quota_probe_path": quota_probe_path,
                    "quota_domain_id": quota_domain_id,
                    "filesystem": filesystem,
                    "query_argv": [
                        "lfs",
                        "quota",
                        "-h",
                        "-u",
                        "fixture",
                        quota_probe_path,
                    ],
                    "query_exit_code": 0,
                    "stdout": stdout,
                    "stderr": stderr,
                    "stdout_sha256": hashlib.sha256(
                        stdout.encode("utf-8")
                    ).hexdigest(),
                    "stderr_sha256": hashlib.sha256(
                        stderr.encode("utf-8")
                    ).hexdigest(),
                    "used_bytes": parsed["used_bytes"],
                    "quota_bytes": parsed["quota_bytes"],
                    "used_inodes": parsed["used_inodes"],
                    "quota_inodes": parsed["quota_inodes"],
                    "usage_authoritative": parsed["usage_authoritative"],
                    "blocker_codes": parsed["blocker_codes"],
                }
            )
        payload = projector.apply_authority(
            {
                "schema": projector.QUOTA_SNAPSHOT_SCHEMA,
                "status": projector.PASS_STATUS,
                "gate": projector.GATE,
                "plan": {"path": plan_path, "sha256": plan_sha},
                "principal": "fixture",
                "observed_at": "2026-07-24T00:00:00+00:00",
                "observed_unix_seconds": 1000,
                "max_age_seconds": 900,
                "storage_domains": domains,
                "blocker_codes": [],
            }
        )
        return seal(payload, "snapshot_semantic_sha256")

    @staticmethod
    def byte_token(value: int) -> str:
        return f"{Decimal(value) / Decimal(1024):f}K"

    def set_quota_observation(
        self,
        record: dict,
        *,
        used_bytes: int,
        quota_bytes: int,
        used_inodes: int,
        quota_inodes: int,
        stderr: str = "",
        returncode: int = 0,
    ) -> None:
        stdout = (
            f"{record['filesystem']} {self.byte_token(used_bytes)} "
            f"{self.byte_token(quota_bytes)} "
            f"{self.byte_token(quota_bytes)} - "
            f"{used_inodes} {quota_inodes} {quota_inodes} -\n"
        )
        parsed = projector.parse_lfs_quota_output(
            stdout, stderr, returncode
        )
        record.update(
            {
                "query_exit_code": returncode,
                "stdout": stdout,
                "stderr": stderr,
                "stdout_sha256": hashlib.sha256(
                    stdout.encode("utf-8")
                ).hexdigest(),
                "stderr_sha256": hashlib.sha256(
                    stderr.encode("utf-8")
                ).hexdigest(),
                "used_bytes": parsed.get("used_bytes"),
                "quota_bytes": parsed.get("quota_bytes"),
                "used_inodes": parsed.get("used_inodes"),
                "quota_inodes": parsed.get("quota_inodes"),
                "usage_authoritative": parsed["usage_authoritative"],
                "blocker_codes": parsed["blocker_codes"],
            }
        )

    def write_projection_inputs(
        self, manifest: dict, snapshot: dict
    ) -> tuple[Path, Path]:
        manifest_path = self.root / "storage_manifest.json"
        snapshot_path = self.root / "quota_snapshot.json"
        manifest_path.write_bytes(projector.canonical_json_bytes(manifest))
        snapshot_path.write_bytes(projector.canonical_json_bytes(snapshot))
        return manifest_path, snapshot_path

    def test_exact_measurement_and_shared_quota_projection_pass(self) -> None:
        manifest = self.passed_manifest()
        self.assertEqual(manifest["status"], projector.PASS_STATUS)
        self.assertEqual(manifest["blocker_codes"], [])
        self.assertEqual(
            manifest["artifact_measurement"]["count_contract"],
            {
                "row_count": 13,
                "source_tuple_count": 129_998,
                "group_size": 7,
                "expected_job_count": 18_577,
                "expected_output_pair_count": 18_577,
                "expected_analysis_output_count": 18_577,
                "expected_sidecar_output_count": 18_577,
                "expected_physical_root_artifact_count": 37_154,
                "expected_source_occurrence_count": 18_577,
                "request_memory_mb": 8_000,
            },
        )
        self.assertEqual(manifest["authority"], projector.EXPECTED_AUTHORITY)
        self.assertFalse(manifest["submission_performed"])
        self.assertFalse(manifest["full_extraction_authority"])

        snapshot = self.quota_snapshot(manifest)
        manifest_path, snapshot_path = self.write_projection_inputs(
            manifest, snapshot
        )
        certificate = projector.build_projection(
            manifest_path, snapshot_path, now_seconds=1001
        )
        self.assertEqual(certificate["status"], projector.PASS_STATUS)
        self.assertEqual(len(certificate["quota_domain_projections"]), 1)
        self.assertEqual(
            certificate["quota_domain_projections"][0]["storage_domains"],
            sorted(projector.STORAGE_DOMAINS),
        )
        expected_increment = sum(
            record["projected_increment_bytes"]
            for record in manifest["artifact_measurement"][
                "storage_domain_totals"
            ].values()
        )
        self.assertEqual(
            certificate["quota_domain_projections"][0][
                "projected_increment_bytes"
            ],
            expected_increment,
        )

    def test_resealed_zero_domain_totals_are_rederived_and_rejected(self) -> None:
        manifest = self.passed_manifest()
        for record in manifest["artifact_measurement"][
            "storage_domain_totals"
        ].values():
            record["projected_increment_bytes"] = 0
            record["projected_increment_inodes"] = 0
        del manifest["artifact_measurement"]["measurement_semantic_sha256"]
        seal(
            manifest["artifact_measurement"],
            "measurement_semantic_sha256",
        )
        del manifest["manifest_semantic_sha256"]
        seal(manifest, "manifest_semantic_sha256")
        with self.assertRaisesRegex(
            projector.ProjectionError, "storage-domain totals differ"
        ):
            projector.validate_storage_manifest(manifest)

    def test_resealed_row_and_controller_arithmetic_are_rejected(self) -> None:
        for mutation, expected_error in (
            ("row", "projected artifact arithmetic differs"),
            ("controller", "controller arithmetic differs"),
        ):
            with self.subTest(mutation=mutation):
                manifest = self.passed_manifest()
                measurement = manifest["artifact_measurement"]
                if mutation == "row":
                    measurement["row_envelopes"][0]["artifact_classes"][0][
                        "projected_bytes"
                    ] += 1
                else:
                    measurement["controller_dry_materialization_budget"][
                        "projected_bytes"
                    ] += 1
                del measurement["measurement_semantic_sha256"]
                seal(measurement, "measurement_semantic_sha256")
                del manifest["manifest_semantic_sha256"]
                seal(manifest, "manifest_semantic_sha256")
                with self.assertRaisesRegex(
                    projector.ProjectionError, expected_error
                ):
                    projector.validate_storage_manifest(manifest)

    def test_coordinated_resealed_witness_underprojection_is_rejected(
        self,
    ) -> None:
        forged_witnesses = copy.deepcopy(self.witnesses)
        for witness in forged_witnesses.values():
            witness["analysis_root"]["size_bytes"] = 1
            witness["training_sidecar_root"]["size_bytes"] = 1
        forged_spec = self.measurement_spec()
        for row in forged_spec["row_envelopes"]:
            for artifact_class in row["artifact_classes"]:
                if artifact_class["artifact_class"] in {
                    "analysis_root",
                    "training_sidecar_root",
                }:
                    artifact_class["bytes_per_artifact_ceiling"] = 2

        forged = projector.build_storage_manifest_from_validated(
            forged_spec,
            plan=self.plan,
            rows=self.rows,
            witnesses=forged_witnesses,
            bindings=self.bindings,
        )
        self.assertEqual(forged["status"], projector.PASS_STATUS)
        self.assertEqual(
            forged["artifact_measurement"]["storage_domain_totals"][
                "bulk_science"
            ]["projected_increment_bytes"],
            81_739,
        )
        with self.assertRaisesRegex(
            projector.ProjectionError,
            "capacity witnesses differ from the bound amendment",
        ):
            projector.validate_storage_manifest(forged)

    def test_bound_amendment_loader_and_readback_are_exact(self) -> None:
        amendment, amendment_artifact, readback_artifact, witnesses = (
            projector.load_revalidated_amendment_chain(
                {
                    "path": self.bindings["capacity_count_amendment"]["path"],
                    "sha256": self.bindings[
                        "capacity_count_amendment"
                    ]["sha256"],
                },
                {
                    "path": self.bindings[
                        "capacity_count_amendment_readback"
                    ]["path"],
                    "sha256": self.bindings[
                        "capacity_count_amendment_readback"
                    ]["sha256"],
                },
            )
        )
        self.assertEqual(
            amendment["schema"], projector.amendment_tool.AMENDMENT_SCHEMA
        )
        self.assertEqual(
            amendment_artifact,
            self.bindings["capacity_count_amendment"],
        )
        self.assertEqual(
            readback_artifact,
            self.bindings["capacity_count_amendment_readback"],
        )
        self.assertEqual(witnesses, self.witnesses)
        self.assertGreaterEqual(self.amendment_validator.call_count, 1)

        readback_path = Path(
            self.bindings["capacity_count_amendment_readback"]["path"]
        )
        readback = json.loads(readback_path.read_text(encoding="utf-8"))
        readback["amendment"] = str(self.root / "different-amendment.json")
        readback_path.write_bytes(projector.canonical_json_bytes(readback))
        with self.assertRaisesRegex(
            projector.ProjectionError,
            "amendment readback differs",
        ):
            projector.load_revalidated_amendment_chain(
                {
                    "path": self.bindings["capacity_count_amendment"]["path"],
                    "sha256": self.bindings[
                        "capacity_count_amendment"
                    ]["sha256"],
                },
                {
                    "path": str(readback_path),
                    "sha256": projector.file_sha256(readback_path),
                },
            )

    def test_path_bound_measurement_evidence_is_rehashed(self) -> None:
        manifest = self.passed_manifest()
        plan_path = Path(
            manifest["artifact_measurement"]["bindings"]["plan"]["path"]
        )
        plan_path.write_text('{"fixture":"changed"}\n', encoding="utf-8")
        with self.assertRaisesRegex(projector.ProjectionError, "SHA-256 differs"):
            projector.validate_storage_manifest(manifest)

    def test_snapshot_root_query_and_raw_output_are_bound(self) -> None:
        manifest = self.passed_manifest()

        unrelated_parent = self.root / "unrelated"
        unrelated_parent.mkdir()
        mismatched = self.quota_snapshot(manifest)
        record = mismatched["storage_domains"][0]
        record["intended_root"] = str(unrelated_parent / "the134")
        record["quota_probe_path"] = str(unrelated_parent)
        record["query_argv"][-1] = str(unrelated_parent)
        del mismatched["snapshot_semantic_sha256"]
        seal(mismatched, "snapshot_semantic_sha256")
        manifest_path, snapshot_path = self.write_projection_inputs(
            manifest, mismatched
        )
        certificate = projector.build_projection(
            manifest_path, snapshot_path, now_seconds=1001
        )
        self.assertEqual(certificate["status"], projector.BLOCKED_STATUS)
        self.assertIn("PLAN_BINDING_DRIFT", certificate["blocker_codes"])

        bad_argv = self.quota_snapshot(manifest)
        bad_argv["storage_domains"][0]["query_argv"][0] = "not-lfs"
        del bad_argv["snapshot_semantic_sha256"]
        seal(bad_argv, "snapshot_semantic_sha256")
        with self.assertRaisesRegex(
            projector.ProjectionError, "query_argv is not canonical"
        ):
            projector.validate_quota_snapshot(bad_argv)

        bad_probe = self.quota_snapshot(manifest)
        bad_probe["storage_domains"][0]["quota_probe_path"] = "/"
        bad_probe["storage_domains"][0]["query_argv"][-1] = "/"
        del bad_probe["snapshot_semantic_sha256"]
        seal(bad_probe, "snapshot_semantic_sha256")
        with self.assertRaisesRegex(
            projector.ProjectionError, "nearest existing parent"
        ):
            projector.validate_quota_snapshot(bad_probe)

        bad_raw = self.quota_snapshot(manifest)
        bad_raw["storage_domains"][0]["stdout"] += "tamper\n"
        del bad_raw["snapshot_semantic_sha256"]
        seal(bad_raw, "snapshot_semantic_sha256")
        with self.assertRaisesRegex(
            projector.ProjectionError, "stdout_sha256 differs"
        ):
            projector.validate_quota_snapshot(bad_raw)

    def test_snapshot_builder_preserves_reparseable_raw_quota_output(
        self,
    ) -> None:
        plan = self.count_plan()
        plan["campaign"] = dict(self.plan["campaign"])
        plan_path = self.root / "count_plan.json"
        plan_path.write_text(
            json.dumps(plan, sort_keys=True) + "\n", encoding="utf-8"
        )
        stdout = (
            "/fixture 1G 20T 20T - "
            "100 20000000 20000000 -\n"
        )

        def runner(argv: list[str], **_kwargs: object) -> subprocess.CompletedProcess:
            return subprocess.CompletedProcess(argv, 0, stdout, "")

        snapshot = projector.build_quota_snapshot(
            plan_path,
            "fixture",
            runner=runner,
            now_seconds=1000,
            lfs_binary="/usr/bin/lfs",
        )
        self.assertEqual(snapshot["status"], projector.PASS_STATUS)
        self.assertTrue(
            all(
                record["stdout"] == stdout
                and record["stderr"] == ""
                for record in snapshot["storage_domains"]
            )
        )
        projector.validate_quota_snapshot(snapshot)

    def test_shared_quota_observations_use_conservative_extrema(self) -> None:
        manifest = self.passed_manifest()
        increment = sum(
            record["projected_increment_bytes"]
            for record in manifest["artifact_measurement"][
                "storage_domain_totals"
            ].values()
        )
        quota_bytes = max(10 * increment, 10_000)
        quota_inodes = 20_000_000
        snapshot = self.quota_snapshot(
            manifest,
            used_bytes=0,
            quota_bytes=quota_bytes,
            used_inodes=0,
            quota_inodes=quota_inodes,
        )
        observations = snapshot["storage_domains"]
        self.set_quota_observation(
            observations[0],
            used_bytes=100,
            quota_bytes=quota_bytes,
            used_inodes=10,
            quota_inodes=quota_inodes,
        )
        self.set_quota_observation(
            observations[1],
            used_bytes=quota_bytes * 4 // 5,
            quota_bytes=quota_bytes,
            used_inodes=100,
            quota_inodes=quota_inodes - 1,
        )
        del snapshot["snapshot_semantic_sha256"]
        seal(snapshot, "snapshot_semantic_sha256")
        manifest_path, snapshot_path = self.write_projection_inputs(
            manifest, snapshot
        )
        certificate = projector.build_projection(
            manifest_path, snapshot_path, now_seconds=1001
        )
        self.assertEqual(certificate["status"], projector.FAIL_STATUS)
        group = certificate["quota_domain_projections"][0]
        self.assertEqual(group["quota_observation_count"], 4)
        self.assertEqual(group["used_bytes"], quota_bytes * 4 // 5)
        self.assertEqual(group["quota_inodes"], quota_inodes - 1)
        self.assertIn("BYTE_HEADROOM_LT_20PCT", certificate["blocker_codes"])

    def test_missing_row_and_artifact_class_block(self) -> None:
        missing_row = self.measurement_spec()
        missing_row["row_envelopes"].pop()
        result = projector.build_storage_manifest_from_validated(
            missing_row,
            plan=self.plan,
            rows=self.rows,
            witnesses=self.witnesses,
            bindings=self.bindings,
        )
        self.assertEqual(result["status"], projector.BLOCKED_STATUS)
        self.assertIn("ROW_ENVELOPE_MISSING", result["blocker_codes"])

        missing_class = self.measurement_spec()
        missing_class["row_envelopes"][0]["artifact_classes"].pop()
        result = projector.build_storage_manifest_from_validated(
            missing_class,
            plan=self.plan,
            rows=self.rows,
            witnesses=self.witnesses,
            bindings=self.bindings,
        )
        self.assertIn("ARTIFACT_CLASS_MISSING", result["blocker_codes"])

    def test_capacity_witness_is_not_silently_used_as_ceiling(self) -> None:
        spec = self.measurement_spec()
        first = spec["row_envelopes"][0]["artifact_classes"][0]
        witness = self.witnesses["pp_background_jet8"]["analysis_root"]
        first["bytes_per_artifact_ceiling"] = witness["size_bytes"]
        result = projector.build_storage_manifest_from_validated(
            spec,
            plan=self.plan,
            rows=self.rows,
            witnesses=self.witnesses,
            bindings=self.bindings,
        )
        self.assertEqual(result["status"], projector.BLOCKED_STATUS)
        self.assertIn(
            "CAPACITY_WITNESS_NOT_STORAGE_CEILING",
            result["blocker_codes"],
        )

    def test_controller_budget_is_mandatory(self) -> None:
        spec = self.measurement_spec()
        spec["controller_dry_materialization_budget"] = None
        result = projector.build_storage_manifest_from_validated(
            spec,
            plan=self.plan,
            rows=self.rows,
            witnesses=self.witnesses,
            bindings=self.bindings,
        )
        self.assertEqual(result["status"], projector.BLOCKED_STATUS)
        self.assertIn("CONTROLLER_BUDGET_MISSING", result["blocker_codes"])

    def test_lustre_authority_signals_are_independently_fail_closed(self) -> None:
        bracketed = (
            "Disk quotas for usr fixture:\n"
            "Filesystem kbytes quota limit grace files quota limit grace\n"
            "/sphenix/tg/tg01 [6.881G] 5T 5T - 6792320 10240000 10240000 -\n"
        )
        parsed = projector.parse_lfs_quota_output(
            bracketed, "", 0
        )
        self.assertFalse(parsed["usage_authoritative"])
        self.assertIn("QUOTA_USAGE_UNAUTHORITATIVE", parsed["blocker_codes"])

        clean = bracketed.replace("[6.881G]", "6.881G")
        parsed = projector.parse_lfs_quota_output(clean, "", 5)
        self.assertFalse(parsed["usage_authoritative"])

        parsed = projector.parse_lfs_quota_output(
            clean, "Some data are inaccurate", 0
        )
        self.assertFalse(parsed["usage_authoritative"])

        parsed = projector.parse_lfs_quota_output(clean, "", 0)
        self.assertTrue(parsed["usage_authoritative"])
        self.assertEqual(parsed["quota_bytes"], 5 * 1024**4)
        self.assertEqual(parsed["quota_inodes"], 10_240_000)

    def test_per_row_chunk_redistribution_is_count_drift(self) -> None:
        plan = self.count_plan()
        _rows, blockers = projector.validate_plan_counts(plan)
        self.assertEqual(blockers, [])

        first = plan["rows"][0]["input_contract"]
        second = plan["rows"][1]["input_contract"]
        for field in (
            "expected_job_count",
            "expected_chunk_count",
            "expected_output_pair_count",
            "expected_analysis_output_count",
            "expected_sidecar_output_count",
            "expected_source_occurrence_count",
        ):
            first[field] += 1
            second[field] -= 1
        _rows, blockers = projector.validate_plan_counts(plan)
        self.assertIn("COUNT_CONTRACT_DRIFT", blockers)

    def test_byte_and_inode_capacity_failures_are_distinct(self) -> None:
        manifest = self.passed_manifest()
        snapshot = self.quota_snapshot(
            manifest,
            used_bytes=900,
            quota_bytes=1_000,
            used_inodes=9_000_000,
            quota_inodes=10_000_000,
        )
        manifest_path, snapshot_path = self.write_projection_inputs(
            manifest, snapshot
        )
        certificate = projector.build_projection(
            manifest_path, snapshot_path, now_seconds=1001
        )
        self.assertEqual(certificate["status"], projector.FAIL_STATUS)
        self.assertIn("BYTE_HEADROOM_LT_20PCT", certificate["blocker_codes"])
        self.assertIn("INODE_HEADROOM_LT_20PCT", certificate["blocker_codes"])
        self.assertEqual(
            projector.exit_for_status(certificate["status"]),
            projector.EXIT_CAPACITY,
        )

    def test_exact_twenty_percent_headroom_passes_one_unit_below_fails(self) -> None:
        manifest = self.passed_manifest()
        totals = manifest["artifact_measurement"]["storage_domain_totals"]
        projected_bytes = sum(
            record["projected_increment_bytes"] for record in totals.values()
        )
        projected_inodes = sum(
            record["projected_increment_inodes"] for record in totals.values()
        )
        quota_bytes = max(5 * projected_bytes, 5)
        quota_inodes = max(5 * projected_inodes, 5)
        used_bytes_at_boundary = 4 * quota_bytes // 5 - projected_bytes
        used_inodes_at_boundary = 4 * quota_inodes // 5 - projected_inodes
        self.assertGreaterEqual(used_bytes_at_boundary, 0)
        self.assertGreaterEqual(used_inodes_at_boundary, 0)

        passing_snapshot = self.quota_snapshot(
            manifest,
            used_bytes=used_bytes_at_boundary,
            quota_bytes=quota_bytes,
            used_inodes=used_inodes_at_boundary,
            quota_inodes=quota_inodes,
        )
        manifest_path, snapshot_path = self.write_projection_inputs(
            manifest, passing_snapshot
        )
        passing = projector.build_projection(
            manifest_path, snapshot_path, now_seconds=1001
        )
        self.assertEqual(passing["status"], projector.PASS_STATUS)

        failing_snapshot = self.quota_snapshot(
            manifest,
            used_bytes=used_bytes_at_boundary + 1,
            quota_bytes=quota_bytes,
            used_inodes=used_inodes_at_boundary + 1,
            quota_inodes=quota_inodes,
        )
        manifest_path, snapshot_path = self.write_projection_inputs(
            manifest, failing_snapshot
        )
        failing = projector.build_projection(
            manifest_path, snapshot_path, now_seconds=1001
        )
        self.assertEqual(failing["status"], projector.FAIL_STATUS)
        self.assertIn("BYTE_HEADROOM_LT_20PCT", failing["blocker_codes"])
        self.assertIn("INODE_HEADROOM_LT_20PCT", failing["blocker_codes"])

    def test_stale_snapshot_and_namespace_reuse_block(self) -> None:
        manifest = self.passed_manifest()
        snapshot = self.quota_snapshot(manifest)
        reused = snapshot["storage_domains"][0]
        Path(reused["intended_root"]).mkdir()
        reused["namespace_fresh"] = False
        reused["quota_probe_path"] = reused["intended_root"]
        reused["query_argv"][-1] = reused["intended_root"]
        reused["blocker_codes"] = ["NAMESPACE_NOT_FRESH"]
        snapshot["blocker_codes"] = ["NAMESPACE_NOT_FRESH"]
        snapshot["status"] = projector.BLOCKED_STATUS
        del snapshot["snapshot_semantic_sha256"]
        seal(snapshot, "snapshot_semantic_sha256")
        manifest_path, snapshot_path = self.write_projection_inputs(
            manifest, snapshot
        )
        certificate = projector.build_projection(
            manifest_path, snapshot_path, now_seconds=2000
        )
        self.assertEqual(certificate["status"], projector.BLOCKED_STATUS)
        self.assertIn("QUOTA_SNAPSHOT_STALE", certificate["blocker_codes"])
        self.assertIn("NAMESPACE_NOT_FRESH", certificate["blocker_codes"])
        self.assertEqual(
            projector.exit_for_status(certificate["status"]),
            projector.EXIT_BLOCKED,
        )

    def test_semantic_tamper_and_authority_escalation_are_malformed(self) -> None:
        manifest = self.passed_manifest()
        tampered = copy.deepcopy(manifest)
        tampered["artifact_measurement"]["storage_domain_totals"][
            "bulk_science"
        ]["projected_increment_bytes"] += 1
        with self.assertRaisesRegex(projector.ProjectionError, "semantic"):
            projector.validate_storage_manifest(tampered)

        escalated = copy.deepcopy(manifest)
        escalated["broad_production_authority"] = True
        del escalated["manifest_semantic_sha256"]
        seal(escalated, "manifest_semantic_sha256")
        with self.assertRaisesRegex(projector.ProjectionError, "authority"):
            projector.validate_storage_manifest(escalated)

    def test_p5a_s_input_is_detected(self) -> None:
        spec = self.measurement_spec()
        spec["bindings"]["science_freeze_certificate"] = {
            "gate": "P5A-S"
        }
        self.assertTrue(projector.contains_p5a_s_input(spec))
        for value in (
            "/evidence/P5A-S/pass.json",
            "gate=P5A_S",
            "P5A S certificate",
            "prefix/P5A-S:PASS",
        ):
            with self.subTest(value=value):
                self.assertTrue(
                    projector.contains_p5a_s_input({"evidence": value})
                )
        self.assertTrue(
            projector.contains_p5a_s_input({"p5a_s_receipt": "present"})
        )
        self.assertFalse(
            projector.contains_p5a_s_input({"evidence": "P5A-control"})
        )

    def test_json_output_is_exclusive_and_cannot_create_campaign_roots(
        self,
    ) -> None:
        output_root = Path(self.plan["campaign"]["output_root"])
        forbidden_output = output_root / "snapshot.json"
        with self.assertRaisesRegex(
            projector.ProjectionError, "outside intended campaign root"
        ):
            projector.write_json(
                forbidden_output,
                {"status": "fixture"},
                forbidden_roots=[str(output_root)],
            )
        self.assertFalse(output_root.exists())

        canonical_root = self.root / "canonical_remote"
        canonical_root.mkdir()
        alias_root = self.root / "remote_alias"
        alias_root.symlink_to(canonical_root, target_is_directory=True)
        canonical_forbidden = canonical_root / "campaign"
        aliased_output = alias_root / "campaign" / "snapshot.json"
        with self.assertRaisesRegex(
            projector.ProjectionError, "outside intended campaign root"
        ):
            projector.write_json(
                aliased_output,
                {"status": "fixture"},
                forbidden_roots=[str(canonical_forbidden)],
            )
        self.assertFalse(canonical_forbidden.exists())

        missing_parent = self.root / "missing" / "result.json"
        with self.assertRaisesRegex(
            projector.ProjectionError, "parent must be"
        ):
            projector.write_json(missing_parent, {"status": "fixture"})
        self.assertFalse(missing_parent.parent.exists())

        output = self.root / "result.json"
        projector.write_json(output, {"status": "fixture"})
        self.assertEqual(
            output.read_bytes(),
            projector.canonical_json_bytes({"status": "fixture"}),
        )
        with self.assertRaisesRegex(
            projector.ProjectionError, "must be new"
        ):
            projector.write_json(output, {"status": "replacement"})


if __name__ == "__main__":
    unittest.main()
