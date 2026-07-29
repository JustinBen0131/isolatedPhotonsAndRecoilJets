#!/usr/bin/env python3
"""Focused fail-closed tests for the guarded THE-134 execution controller."""

from __future__ import annotations

import argparse
import copy
import hashlib
import importlib.util
import json
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest import mock


HERE = Path(__file__).resolve().parent
CONTROLLER_PATH = HERE.parent / "the134_full_extraction_controller.py"
FIXTURE_PATH = HERE / "test_materialize_the134_full_multiview_extraction.py"
PREPARER_PATH = (
    CONTROLLER_PATH.parents[3]
    / "ml"
    / "training"
    / "prepare_the134_h70_matrix.py"
)


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


controller = load_module("the134_full_extraction_controller_tested", CONTROLLER_PATH)
fixture_module = load_module("the134_full_extraction_fixture", FIXTURE_PATH)
preparer = load_module("prepare_the134_h70_matrix_for_controller_test", PREPARER_PATH)


def sha256_bytes(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class FullExtractionControllerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.temporary = tempfile.TemporaryDirectory()
        cls.root = Path(cls.temporary.name).resolve()
        cls.fixture = fixture_module.FullControllerFixture(cls.root / "fixture")
        cls.validation_root = cls.root / "validation_must_remain_absent"
        materializer_args = cls.fixture.args(cls.validation_root)
        cls.base_context = controller.materializer.validate_plan_and_evidence(
            materializer_args
        )
        cls.dry_root = cls.root / "dry_materialization"
        controller.materializer.materialize_local(
            cls.dry_root,
            controller.materializer.build_staged_artifacts(cls.base_context),
        )
        cls.args = argparse.Namespace(
            plan=cls.fixture.plan_path,
            plan_sha256=file_sha256(cls.fixture.plan_path),
            preflight_receipt=cls.fixture.receipt_path,
            preflight_receipt_sha256=file_sha256(cls.fixture.receipt_path),
            storage_certificate=cls.fixture.storage_certificate_path,
            storage_certificate_sha256=file_sha256(
                cls.fixture.storage_certificate_path
            ),
            controller_budget=cls.fixture.budget_path,
            controller_budget_sha256=file_sha256(cls.fixture.budget_path),
            artifact_profile=controller.ARTIFACT_PROFILE,
            validation_root=cls.validation_root,
            dry_staging_root=cls.dry_root,
            dry_manifest_sha256=file_sha256(
                cls.dry_root / "materialization_manifest.json"
            ),
        )
        cls.context = controller.materializer_context(cls.args)
        cls.dry_stage = controller.validate_dry_stage(
            cls.context,
            cls.dry_root,
            cls.args.dry_manifest_sha256,
        )
        cls.authorization_payload = cls.make_authorization_payload()
        cls.authorization_path = cls.root / "submission_authorization.json"
        cls.authorization_path.write_bytes(
            controller.canonical_json_bytes(cls.authorization_payload)
        )
        cls.authorization_artifact = {
            "path": str(cls.authorization_path),
            "sha256": file_sha256(cls.authorization_path),
            "size_bytes": cls.authorization_path.stat().st_size,
        }
        cls.authorization = controller.validate_authorization(
            cls.authorization_payload,
            artifact=cls.authorization_artifact,
            context=cls.context,
            dry_stage=cls.dry_stage,
            now_seconds=1_800_000_000,
        )
        cls.execution_payload = controller.build_execution_manifest(
            cls.context,
            cls.dry_stage,
            cls.authorization,
        )
        cls.execution = controller.validate_execution_manifest(
            cls.execution_payload
        )
        cls.execution_artifact = {
            "path": str(cls.root / "execution_manifest.json"),
            "sha256": "b" * 64,
            "size_bytes": len(controller.canonical_json_bytes(cls.execution_payload)),
        }

    @classmethod
    def tearDownClass(cls) -> None:
        for path in sorted(
            cls.root.rglob("*"),
            key=lambda item: len(item.parts),
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

    @classmethod
    def make_authorization_payload(cls) -> dict:
        return {
            "schema": controller.AUTHORIZATION_SCHEMA,
            "status": "PASS_EXACT_SUBMISSION_AUTHORIZED",
            "campaign": copy.deepcopy(cls.context["campaign"]),
            "bindings": {
                "plan_sha256": cls.context["plan_artifact"]["sha256"],
                "dry_materialization_manifest_sha256": cls.dry_stage["manifest"][
                    "sha256"
                ],
                "storage_certificate_sha256": cls.context[
                    "storage_certificate_artifact"
                ]["sha256"],
                "duplicate_fingerprint_sha256": cls.context[
                    "duplicate_fingerprint_sha256"
                ],
                "execution_fingerprint_sha256": cls.context[
                    "execution_fingerprint_sha256"
                ],
            },
            "limits": {
                "row_count": controller.EXPECTED_ROW_COUNT,
                "job_count": controller.EXPECTED_JOB_COUNT,
                "group_size": controller.EXPECTED_GROUP_SIZE,
                "request_memory_mb": controller.EXPECTED_REQUEST_MEMORY_MB,
                "output_pair_count": controller.EXPECTED_OUTPUT_PAIRS,
            },
            "duplicate_guard": {
                "status": "PASS_NO_ACTIVE_DUPLICATE",
                "active_matching_jobs": 0,
                "fingerprint_sha256": cls.context[
                    "duplicate_fingerprint_sha256"
                ],
            },
            "capacity_guard": {
                "status": "PASS",
                "group_size": controller.EXPECTED_GROUP_SIZE,
                "request_memory_mb": controller.EXPECTED_REQUEST_MEMORY_MB,
            },
            "quota_guard": {
                "status": "PASS",
                "authoritative": True,
                "certificate_sha256": cls.context[
                    "storage_certificate_artifact"
                ]["sha256"],
            },
            "approval": {
                "explicit": True,
                "scope": "THE134_FULL_SOURCE_COMPLETE_GROUP7_EXTRACTION",
                "principal": "Justin",
                "evidence_sha256": "c" * 64,
            },
            "provenance": {
                "codex_chat_name": "THE-114-THE-134",
                "codex_thread_id": "019f80b5-dc56-7330-9ee7-56ef417547dc",
            },
            "expires_at_unix": 2_000_000_000,
        }

    @classmethod
    def make_submission_receipt(cls) -> dict:
        receipt = controller.submission_receipt_base(
            cls.execution_artifact,
            cls.execution,
        )
        receipt["status"] = controller.SUBMITTED_STATUS
        receipt["submission_performed"] = True
        receipt["attempt_lock"] = {
            "path": (
                f"{cls.execution['campaign']['evidence_root']}/"
                "the134_full_extraction_submission_attempt.lock"
            ),
            "sha256": "9" * 64,
            "size_bytes": 1000,
        }
        receipt["first_bad"] = None
        receipt["counts"] = {
            "expected_row_count": controller.EXPECTED_ROW_COUNT,
            "expected_job_count": controller.EXPECTED_JOB_COUNT,
            "submitted_row_count": controller.EXPECTED_ROW_COUNT,
            "submitted_job_count": controller.EXPECTED_JOB_COUNT,
        }
        provenance = cls.execution["authorization"]["provenance"]
        for row_index, row in enumerate(cls.execution["_rows"]):
            receipt["rows"].append(
                {
                    "row_id": row["row_id"],
                    "returncode": 0,
                    "started_at_unix": 1_900_000_000 + row_index,
                    "finished_at_unix": 1_900_000_001 + row_index,
                    "argv": list(row["submitter_argv"]),
                    "row_fingerprint_sha256": row["row_fingerprint_sha256"],
                    "inherited_environment": {},
                    "environment_sha256": controller.canonical_sha256(
                        controller.sealed_submit_environment(
                            row,
                            provenance,
                            ambient={},
                        )
                    ),
                    "stdout": {
                        "path": (
                            f"{cls.execution['campaign']['evidence_root']}/"
                            f"submission/{row['row_id']}.stdout.txt"
                        ),
                        "sha256": sha256_bytes(f"stdout:{row['row_id']}"),
                        "size_bytes": 80,
                    },
                    "stderr": {
                        "path": (
                            f"{cls.execution['campaign']['evidence_root']}/"
                            f"submission/{row['row_id']}.stderr.txt"
                        ),
                        "sha256": sha256_bytes(f"stderr:{row['row_id']}"),
                        "size_bytes": 0,
                    },
                    "expected_job_count": row["expected_job_count"],
                    "cluster_id": 7_000_000 + row_index,
                    "submitted_job_count": row["expected_job_count"],
                }
            )
        return receipt

    @classmethod
    def make_terminal_snapshot(
        cls,
        submission: dict,
        submission_artifact: dict,
    ) -> dict:
        jobs = []
        submission_by_row = {
            row["row_id"]: row for row in submission["rows"]
        }
        for row in cls.execution["_rows"]:
            cluster = submission_by_row[row["row_id"]]["cluster_id"]
            for row_job_index in range(row["expected_job_count"]):
                global_index = row["global_job_start"] + row_job_index
                staged_chunk_sha256 = sha256_bytes(
                    f"staged-chunk:{global_index}"
                )
                expected_chunk = cls.execution["_partition_chunks"][global_index]
                sidecar_path = (
                    row["training_sidecar_template"]
                    .replace("$(Cluster)", str(cluster))
                    .replace("$(Process)", str(row_job_index))
                )
                jobs.append(
                    {
                        "row_id": row["row_id"],
                        "global_job_index": global_index,
                        "row_job_index": row_job_index,
                        "cluster_id": cluster,
                        "proc_id": row_job_index,
                        "job_status": 4,
                        "exit_code": 0,
                        "num_job_starts": 1,
                        "staged_chunk_list": {
                            "path": (
                                f"{row['submit_namespace']}/chunks/"
                                f"{cluster}.{row_job_index}.list"
                            ),
                            "size_bytes": 100 + global_index,
                            "sha256": staged_chunk_sha256,
                            "status": "PASS",
                            "readable": True,
                            "execution_chunk_sha256": expected_chunk[
                                "execution_chunk_sha256"
                            ],
                            "tuple_input_sha256s": list(
                                expected_chunk["tuple_input_sha256s"]
                            ),
                        },
                        "analysis_root": {
                            "path": (
                                f"{row['analysis_output_namespace']}/analysis/"
                                f"{cluster}.{row_job_index}.root"
                            ),
                            "size_bytes": 50_000 + global_index,
                            "sha256": sha256_bytes(f"analysis:{global_index}"),
                            "status": "PASS",
                            "readable": True,
                            "zombie": False,
                            "recovered": False,
                        },
                        "training_sidecar_root": {
                            "path": sidecar_path,
                            "size_bytes": 1_000 + global_index,
                            "sha256": sha256_bytes(f"sidecar:{global_index}"),
                            "status": "PASS",
                            "readable": True,
                            "zombie": False,
                            "recovered": False,
                        },
                        "sidecar_health": {
                            "schema": "THE134_FULL_EXTRACTION_SIDECAR_HEALTH_V1",
                            "status": "PASS",
                            "profile": "photon_training_multiview_v1",
                            "certificate_sha256": sha256_bytes(
                                f"certificate:{global_index}"
                            ),
                            "tree_name": "RJPhotonTrainingViewV1",
                            "tree_entries": len(controller.SHOWER_VIEWS) * 2,
                            "source_occurrence_id_hex": f"{global_index + 1:032x}",
                            "input_uri_sha256": staged_chunk_sha256,
                            "input_file_sha256": staged_chunk_sha256,
                            "source_manifest_sha256": row[
                                "full_source_manifest_sha256"
                            ],
                            "config_sha256": row["config_sha256"],
                            "code_sha256": row["code_sha256"],
                            "selected_training_rows_by_view": {
                                view: 2 for view in controller.SHOWER_VIEWS
                            },
                        },
                    }
                )
        return {
            "schema": controller.TERMINAL_SNAPSHOT_SCHEMA,
            "status": "COMPLETE",
            "campaign": copy.deepcopy(cls.execution["campaign"]),
            "execution_manifest_sha256": cls.execution_artifact["sha256"],
            "submission_receipt_sha256": submission_artifact["sha256"],
            "captured_at_unix": 1_900_100_000,
            "jobs": jobs,
        }

    def test_preflight_rebinds_exact_18577_job_dry_stage(self) -> None:
        self.assertEqual(len(self.context["rows"]), controller.EXPECTED_ROW_COUNT)
        self.assertEqual(len(self.context["chunks"]), controller.EXPECTED_JOB_COUNT)
        self.assertEqual(
            self.dry_stage["manifest"]["sha256"],
            self.args.dry_manifest_sha256,
        )
        self.assertEqual(
            self.execution["counts"]["request_memory_mb"],
            controller.EXPECTED_REQUEST_MEMORY_MB,
        )
        self.assertFalse(self.execution["authority"]["submission_performed"])
        self.assertTrue(self.execution["authority"]["submission_authority"])
        self.assertFalse(self.execution["authority"]["full_extraction_authority"])

    def test_dry_stage_byte_mutation_fails_closed(self) -> None:
        mutated_root = self.root / f"mutated_{self._testMethodName}"
        shutil.copytree(self.dry_root, mutated_root)
        jobs_path = mutated_root / "the134_full_multiview_jobs.jsonl"
        jobs_path.write_bytes(jobs_path.read_bytes() + b" ")
        with self.assertRaisesRegex(
            controller.ControllerError,
            "dry staged artifact bytes differ",
        ):
            controller.validate_dry_stage(
                self.context,
                mutated_root,
                self.args.dry_manifest_sha256,
            )

    def test_authorization_mutations_fail_closed(self) -> None:
        mutations = {
            "active duplicate": lambda payload: payload["duplicate_guard"].update(
                {"active_matching_jobs": 1}
            ),
            "non-authoritative quota": lambda payload: payload[
                "quota_guard"
            ].update({"authoritative": False}),
            "implicit approval": lambda payload: payload["approval"].update(
                {"explicit": False}
            ),
            "unknown provenance": lambda payload: payload["provenance"].update(
                {"codex_thread_id": "unknown"}
            ),
            "wrong execution fingerprint": lambda payload: payload[
                "bindings"
            ].update({"execution_fingerprint_sha256": "f" * 64}),
            "extra authority field": lambda payload: payload.update(
                {"full_training_authority": 1}
            ),
        }
        for label, mutate in mutations.items():
            with self.subTest(label=label):
                payload = copy.deepcopy(self.authorization_payload)
                mutate(payload)
                with self.assertRaises(controller.ControllerError):
                    controller.validate_authorization(
                        payload,
                        artifact=self.authorization_artifact,
                        context=self.context,
                        dry_stage=self.dry_stage,
                        now_seconds=1_800_000_000,
                    )
        with self.assertRaisesRegex(controller.ControllerError, "expired"):
            controller.validate_authorization(
                self.authorization_payload,
                artifact=self.authorization_artifact,
                context=self.context,
                dry_stage=self.dry_stage,
                now_seconds=2_000_000_000,
            )

    def test_submit_report_and_sealed_environment_are_exact(self) -> None:
        self.assertEqual(
            controller.parse_submit_report(
                "1429 job(s) submitted to cluster 7654321.\n",
                1429,
            ),
            (7_654_321, 1429),
        )
        for output in (
            "1428 job(s) submitted to cluster 7654321.",
            "no cluster was submitted",
            (
                "1429 job(s) submitted to cluster 7654321.\n"
                "1429 job(s) submitted to cluster 7654322."
            ),
        ):
            with self.subTest(output=output):
                with self.assertRaises(controller.ControllerError):
                    controller.parse_submit_report(output, 1429)

        row = self.execution["_rows"][0]
        provenance = self.execution["authorization"]["provenance"]
        environment = controller.sealed_submit_environment(
            row,
            provenance,
            ambient={
                "PATH": "/usr/bin",
                "HOME": "/tmp/home",
                "RJ_UNSEALED": "bad",
                "CODEX_BAD": "bad",
                "BASH_ENV": "/tmp/hostile.sh",
                "LD_PRELOAD": "/tmp/hostile.so",
            },
        )
        self.assertEqual(environment["PATH"], "/usr/bin")
        self.assertEqual(environment["HOME"], "/tmp/home")
        self.assertNotIn("RJ_UNSEALED", environment)
        self.assertNotIn("CODEX_BAD", environment)
        self.assertNotIn("BASH_ENV", environment)
        self.assertNotIn("LD_PRELOAD", environment)
        self.assertEqual(environment["RJ_DAG_DRYRUN"], "0")
        self.assertEqual(
            environment["RJ_CODEX_THREAD_ID"],
            provenance["codex_thread_id"],
        )

    def test_fresh_evidence_tree_rejects_symlinked_parent(self) -> None:
        real_parent = self.root / "real_evidence_parent"
        real_parent.mkdir()
        symlink_parent = self.root / "symlinked_evidence_parent"
        symlink_parent.symlink_to(real_parent, target_is_directory=True)
        with self.assertRaisesRegex(
            controller.ControllerError,
            "nearest existing parent must be a real directory",
        ):
            controller.safe_fresh_local_tree_output(
                symlink_parent / "qa" / "campaign",
                "campaign evidence root",
            )

    def test_staged_chunk_bytes_bind_exact_partition_membership(self) -> None:
        submit_root = self.root / "staged_chunk_submit_root"
        chunk_root = submit_root / "chunks"
        chunk_root.mkdir(parents=True)
        roles = tuple(controller.materializer.resolver.LIST_ROLES)
        inputs = [
            {
                role: f"/frozen/{role}/tuple{tuple_index}.root"
                for role in roles
            }
            for tuple_index in range(2)
        ]
        data = "".join(
            "\t".join(record[role] for role in roles) + "\n"
            for record in inputs
        ).encode("utf-8")
        chunk_path = chunk_root / "7000000.0.list"
        chunk_path.write_bytes(data)
        tuple_hashes = [
            controller.materializer.resolver.canonical_sha256(record)
            for record in inputs
        ]
        execution_chunk_sha256 = "8" * 64
        expected_chunk = {
            "tuple_count": 2,
            "tuple_input_sha256s": tuple_hashes,
            "execution_chunk_sha256": execution_chunk_sha256,
        }
        artifact = {
            "path": str(chunk_path),
            "size_bytes": len(data),
            "sha256": file_sha256(chunk_path),
            "status": "PASS",
            "readable": True,
            "execution_chunk_sha256": execution_chunk_sha256,
            "tuple_input_sha256s": tuple_hashes,
        }
        validated = controller.validate_staged_chunk_list(
            artifact,
            expected_chunk=expected_chunk,
            row={"submit_namespace": str(submit_root)},
            label="staged chunk",
        )
        self.assertEqual(validated["sha256"], file_sha256(chunk_path))

        wrong_partition = copy.deepcopy(artifact)
        wrong_partition["execution_chunk_sha256"] = "7" * 64
        with self.assertRaisesRegex(
            controller.ControllerError,
            "frozen partition identity differs",
        ):
            controller.validate_staged_chunk_list(
                wrong_partition,
                expected_chunk=expected_chunk,
                row={"submit_namespace": str(submit_root)},
                label="staged chunk",
            )

        chunk_path.write_bytes(data.replace(b"tuple0", b"tuple9", 1))
        changed = copy.deepcopy(artifact)
        changed["sha256"] = file_sha256(chunk_path)
        with self.assertRaisesRegex(
            controller.ControllerError,
            "bytes differ from frozen tuple inputs",
        ):
            controller.validate_staged_chunk_list(
                changed,
                expected_chunk=expected_chunk,
                row={"submit_namespace": str(submit_root)},
                label="staged chunk",
            )

    def test_submit_reopens_every_bound_prerequisite(self) -> None:
        original = self.fixture.plan_path.read_bytes()
        self.fixture.plan_path.write_bytes(original + b" ")
        runner_called = False

        def forbidden_runner(*_args, **_kwargs):
            nonlocal runner_called
            runner_called = True
            raise AssertionError("runner must not be reached")

        try:
            with self.assertRaisesRegex(
                controller.ControllerError,
                "execution binding plan SHA-256 changed",
            ):
                controller.execute_submission(
                    self.execution,
                    self.execution_artifact,
                    Path(self.execution["campaign"]["evidence_root"])
                    / "binding_failure",
                    runner=forbidden_runner,
                    now_seconds=1_900_000_000,
                )
        finally:
            self.fixture.plan_path.write_bytes(original)
        self.assertFalse(runner_called)

    def test_submit_rejects_self_consistent_row_not_in_sealed_plan(self) -> None:
        candidate_payload = copy.deepcopy(self.execution_payload)
        row = candidate_payload["rows"][0]
        replacement = Path("/bin/echo")
        self.assertTrue(replacement.is_file())
        row["submitter"]["path"] = str(replacement)
        row["submitter"]["sha256"] = controller.file_sha256(replacement)
        row["submitter_argv"][0] = str(replacement)
        candidate = controller.validate_execution_manifest(candidate_payload)
        runner_called = False

        def forbidden_runner(*_args, **_kwargs):
            nonlocal runner_called
            runner_called = True
            raise AssertionError("runner must not be reached")

        with self.assertRaisesRegex(
            controller.ControllerError,
            "differs from reconstructed materializer output",
        ):
            controller.execute_submission(
                candidate,
                self.execution_artifact,
                Path(candidate["campaign"]["evidence_root"])
                / "forged_row_failure",
                runner=forbidden_runner,
                now_seconds=1_900_000_000,
            )
        self.assertFalse(runner_called)

    def test_submit_requires_switch_and_fake_runner_is_exact_once(self) -> None:
        execution_path = self.root / "execution_manifest_cli.json"
        execution_path.write_bytes(
            controller.canonical_json_bytes(self.execution_payload)
        )
        refused_receipt_dir = Path(
            self.execution["campaign"]["evidence_root"]
        ) / "refused_submission"
        self.assertEqual(
            controller.main(
                [
                    "submit",
                    "--execution-manifest",
                    str(execution_path),
                    "--execution-manifest-sha256",
                    file_sha256(execution_path),
                    "--receipt-dir",
                    str(refused_receipt_dir),
                ]
            ),
            2,
        )
        self.assertFalse(refused_receipt_dir.exists())

        evidence_root = Path(self.execution["campaign"]["evidence_root"])
        calls: list[list[str]] = []

        def fake_runner(argv, **kwargs):
            expected = self.execution["_rows"][len(calls)]
            self.assertEqual(argv, expected["submitter_argv"])
            self.assertFalse(kwargs["check"])
            self.assertEqual(kwargs["env"]["RJ_DAG_DRYRUN"], "0")
            calls.append(list(argv))
            return subprocess.CompletedProcess(
                argv,
                0,
                (
                    f"{expected['expected_job_count']} job(s) submitted to "
                    f"cluster {8_000_000 + len(calls) - 1}.\n"
                ),
                "",
            )

        receipt = controller.execute_submission(
            self.execution,
            self.execution_artifact,
            evidence_root / "fake_submission",
            runner=fake_runner,
            now_seconds=1_900_000_000,
        )
        self.assertEqual(len(calls), controller.EXPECTED_ROW_COUNT)
        self.assertEqual(receipt["status"], controller.SUBMITTED_STATUS)
        self.assertEqual(
            receipt["counts"]["submitted_job_count"],
            controller.EXPECTED_JOB_COUNT,
        )
        controller.validate_submission_receipt(
            receipt,
            self.execution,
            self.execution_artifact,
        )
        attempt_lock = (
            evidence_root
            / "the134_full_extraction_submission_attempt.lock"
        )
        original_attempt_lock = attempt_lock.read_bytes()
        attempt_lock.write_bytes(original_attempt_lock + b" ")
        with self.assertRaisesRegex(
            controller.ControllerError,
            "attempt-lock artifact differs",
        ):
            controller.validate_submission_receipt(
                receipt,
                self.execution,
                self.execution_artifact,
            )
        attempt_lock.write_bytes(original_attempt_lock)

        replay_called = False

        def replay_runner(*_args, **_kwargs):
            nonlocal replay_called
            replay_called = True
            raise AssertionError("persistent campaign lock must reject replay")

        with self.assertRaisesRegex(
            controller.ControllerError,
            "submission-attempt lock already exists",
        ):
            controller.execute_submission(
                self.execution,
                self.execution_artifact,
                evidence_root / "forbidden_replay",
                runner=replay_runner,
                now_seconds=1_900_000_000,
            )
        self.assertFalse(replay_called)

        # Isolated test reset only. Production code intentionally has no lock
        # removal path. Remove the test-created evidence parent as well so the
        # frozen nearest-existing-parent quota witness remains unchanged for
        # this independent nonzero-return fixture.
        shutil.rmtree(evidence_root.parent)
        failed_calls = 0

        def nonzero_after_submit_runner(argv, **_kwargs):
            nonlocal failed_calls
            failed_calls += 1
            expected = self.execution["_rows"][0]
            return subprocess.CompletedProcess(
                argv,
                9,
                (
                    f"{expected['expected_job_count']} job(s) submitted to "
                    "cluster 9000000.\n"
                ),
                "",
            )

        with self.assertRaisesRegex(
            controller.ControllerError,
            "exited 9",
        ):
            controller.execute_submission(
                self.execution,
                self.execution_artifact,
                evidence_root / "fake_partial_failure",
                runner=nonzero_after_submit_runner,
                now_seconds=1_900_000_000,
            )
        self.assertEqual(failed_calls, 1)
        partial = json.loads(
            (
                evidence_root
                / "fake_partial_failure"
                / "submission_receipt.json"
            ).read_text(encoding="utf-8")
        )
        self.assertEqual(partial["status"], "PARTIAL_FAILED_HARD_STOP")
        self.assertEqual(len(partial["rows"]), 1)
        self.assertEqual(partial["rows"][0]["cluster_id"], 9_000_000)
        self.assertEqual(
            partial["rows"][0]["submitted_job_count"],
            self.execution["_rows"][0]["expected_job_count"],
        )
        self.assertTrue(
            partial["first_bad"]["cluster_identity_preserved"]
        )

    def test_submission_receipt_accounting_is_fail_closed(self) -> None:
        payload = self.make_submission_receipt()
        normalized = controller.validate_submission_receipt(
            payload,
            self.execution,
            self.execution_artifact,
            verify_attempt_lock=False,
        )
        self.assertEqual(len(normalized["_rows"]), controller.EXPECTED_ROW_COUNT)

        mutations = {
            "wrong aggregate count": lambda candidate: candidate["counts"].update(
                {"submitted_job_count": controller.EXPECTED_JOB_COUNT - 1}
            ),
            "duplicate cluster": lambda candidate: candidate["rows"][1].update(
                {"cluster_id": candidate["rows"][0]["cluster_id"]}
            ),
            "wrong row environment": lambda candidate: candidate["rows"][0].update(
                {"environment_sha256": "f" * 64}
            ),
            "unexpected authority": lambda candidate: candidate.update(
                {"science_freeze_authority": True}
            ),
        }
        for label, mutate in mutations.items():
            with self.subTest(label=label):
                candidate = copy.deepcopy(payload)
                mutate(candidate)
                with self.assertRaises(controller.ControllerError):
                    controller.validate_submission_receipt(
                        candidate,
                        self.execution,
                        self.execution_artifact,
                        verify_attempt_lock=False,
                    )

    def test_terminal_and_exact_input_once_logical_aggregate(self) -> None:
        submission_payload = self.make_submission_receipt()
        submission = controller.validate_submission_receipt(
            submission_payload,
            self.execution,
            self.execution_artifact,
            verify_attempt_lock=False,
        )
        submission_artifact = {
            "path": str(self.root / "submission_receipt.json"),
            "sha256": "d" * 64,
            "size_bytes": len(controller.canonical_json_bytes(submission_payload)),
        }
        snapshot_payload = self.make_terminal_snapshot(
            submission_payload,
            submission_artifact,
        )
        snapshot_artifact = {
            "path": str(self.root / "terminal_snapshot.json"),
            "sha256": "e" * 64,
            "size_bytes": len(controller.canonical_json_bytes(snapshot_payload)),
        }

        def validate_snapshot(candidate):
            def synthetic_staged_chunk_readback(
                artifact,
                *,
                expected_chunk,
                row,
                label,
            ):
                record = dict(artifact)
                self.assertTrue(
                    controller.path_below(
                        record["path"], row["submit_namespace"]
                    ),
                    label,
                )
                self.assertEqual(
                    record["execution_chunk_sha256"],
                    expected_chunk["execution_chunk_sha256"],
                )
                self.assertEqual(
                    record["tuple_input_sha256s"],
                    expected_chunk["tuple_input_sha256s"],
                )
                return record

            with mock.patch.object(
                controller,
                "validate_staged_chunk_list",
                side_effect=synthetic_staged_chunk_readback,
            ):
                return controller.validate_terminal_snapshot(
                    candidate,
                    execution=self.execution,
                    execution_artifact=self.execution_artifact,
                    submission=submission,
                    submission_artifact=submission_artifact,
                )

        first_selected = snapshot_payload["jobs"][0]["sidecar_health"][
            "selected_training_rows_by_view"
        ]
        first_selected[controller.SHOWER_VIEWS[0]] += 1
        with self.assertRaisesRegex(
            controller.ControllerError,
            "selected training population differs",
        ):
            validate_snapshot(snapshot_payload)
        first_selected[controller.SHOWER_VIEWS[0]] -= 1
        duplicate_input = snapshot_payload["jobs"][1]["sidecar_health"]
        original_input_uri = duplicate_input["input_uri_sha256"]
        original_input_file = duplicate_input["input_file_sha256"]
        duplicate_input["input_file_sha256"] = "a" * 64
        with self.assertRaisesRegex(
            controller.ControllerError,
            "staged-chunk input identity differs",
        ):
            validate_snapshot(snapshot_payload)
        duplicate_input["input_file_sha256"] = original_input_file
        duplicate_input["input_uri_sha256"] = snapshot_payload["jobs"][0][
            "sidecar_health"
        ]["input_uri_sha256"]
        duplicate_input["input_file_sha256"] = snapshot_payload["jobs"][0][
            "sidecar_health"
        ]["input_file_sha256"]
        with self.assertRaisesRegex(
            controller.ControllerError,
            "sidecar input identity differs from the frozen staged chunk",
        ):
            validate_snapshot(snapshot_payload)
        duplicate_staged_chunk = snapshot_payload["jobs"][1]["staged_chunk_list"]
        original_staged_chunk_sha256 = duplicate_staged_chunk["sha256"]
        duplicate_staged_chunk["sha256"] = snapshot_payload["jobs"][0][
            "staged_chunk_list"
        ]["sha256"]
        with self.assertRaisesRegex(
            controller.ControllerError,
            "duplicates a staged-chunk input identity",
        ):
            validate_snapshot(snapshot_payload)
        duplicate_staged_chunk["sha256"] = original_staged_chunk_sha256
        duplicate_input["input_uri_sha256"] = original_input_uri
        duplicate_input["input_file_sha256"] = original_input_file

        first_health = snapshot_payload["jobs"][0]["sidecar_health"]
        first_health["tree_entries"] = 0
        first_health["selected_training_rows_by_view"] = {
            view: 0 for view in controller.SHOWER_VIEWS
        }
        _, jobs = validate_snapshot(snapshot_payload)
        self.assertEqual(len(jobs), controller.EXPECTED_JOB_COUNT)
        terminal_payload = controller.terminal_receipt(
            execution_artifact=self.execution_artifact,
            submission_artifact=submission_artifact,
            snapshot_artifact=snapshot_artifact,
            jobs=jobs,
        )
        terminal_artifact = {
            "path": str(self.root / "terminal_receipt.json"),
            "sha256": "f" * 64,
            "size_bytes": len(controller.canonical_json_bytes(terminal_payload)),
        }
        controller.validate_terminal_receipt(
            terminal_payload,
            execution_artifact=self.execution_artifact,
            submission_artifact=submission_artifact,
            snapshot_artifact=snapshot_artifact,
        )
        artifacts = controller.build_logical_aggregate(
            execution=self.execution,
            jobs=jobs,
            execution_artifact=self.execution_artifact,
            submission_artifact=submission_artifact,
            snapshot_artifact=snapshot_artifact,
            terminal_artifact=terminal_artifact,
        )
        self.assertEqual(
            set(artifacts),
            {
                "pp_sidecars.list",
                "auau_sidecars.list",
                "pp_source_provenance.json",
                "auau_source_provenance.json",
                "logical_aggregate_receipt.json",
            },
        )
        pp_paths = artifacts["pp_sidecars.list"].decode("utf-8").splitlines()
        auau_paths = artifacts["auau_sidecars.list"].decode("utf-8").splitlines()
        self.assertEqual(
            len(pp_paths) + len(auau_paths),
            controller.EXPECTED_JOB_COUNT,
        )
        self.assertEqual(
            len(set(pp_paths + auau_paths)),
            controller.EXPECTED_JOB_COUNT,
        )
        aggregate = json.loads(artifacts["logical_aggregate_receipt.json"])
        self.assertEqual(aggregate["status"], controller.AGGREGATE_STATUS)
        self.assertFalse(
            aggregate["merge_semantics"]["physical_hadd_performed"]
        )
        for system in ("pp", "auau"):
            provenance = json.loads(
                artifacts[f"{system}_source_provenance.json"]
            )
            self.assertEqual(
                set(provenance["source_coverage_authority"]),
                set(controller.h70_contract.expected_sources(system)),
            )
            self.assertEqual(provenance["schema"], "THE134_SOURCE_PROVENANCE_V1")
            provenance_path = self.root / f"{system}_source_provenance.json"
            provenance_path.write_bytes(
                artifacts[f"{system}_source_provenance.json"]
            )
            indexed, authority = preparer.load_source_provenance(
                provenance_path
            )
            paths = [
                Path(path)
                for path in artifacts[f"{system}_sidecars.list"]
                .decode("utf-8")
                .splitlines()
            ]
            sources, failures = preparer.resolve_path_sources(
                paths,
                system,
                indexed,
            )
            self.assertEqual(failures, [])
            self.assertEqual(len(sources), len(paths))
            for source, source_authority in authority.items():
                records = [
                    record
                    for record in provenance["inputs"]
                    if record["source_sample"] == source
                ]
                self.assertEqual(
                    source_authority["input_records_sha256"],
                    preparer.source_input_records_sha256(records),
                )
                self.assertTrue(
                    all(
                        preparer.IDENTITY128_HEX_RE.fullmatch(
                            record["source_occurrence_id_hex"]
                        )
                        is not None
                        for record in records
                    )
                )


if __name__ == "__main__":
    unittest.main()
