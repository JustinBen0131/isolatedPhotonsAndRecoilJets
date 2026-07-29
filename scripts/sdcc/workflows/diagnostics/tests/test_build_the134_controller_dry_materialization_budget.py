#!/usr/bin/env python3
"""Focused tests for the exact THE-134 controller-budget builder."""

from __future__ import annotations

import copy
import hashlib
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


HERE = Path(__file__).resolve().parent
BUILDER_PATH = (
    HERE.parent / "build_the134_controller_dry_materialization_budget.py"
)
BUILDER_SPEC = importlib.util.spec_from_file_location(
    "the134_controller_budget_builder", BUILDER_PATH
)
assert BUILDER_SPEC is not None and BUILDER_SPEC.loader is not None
builder = importlib.util.module_from_spec(BUILDER_SPEC)
BUILDER_SPEC.loader.exec_module(builder)

FIXTURE_PATH = HERE / "test_materialize_the134_full_multiview_extraction.py"
FIXTURE_SPEC = importlib.util.spec_from_file_location(
    "the134_materializer_fixture", FIXTURE_PATH
)
assert FIXTURE_SPEC is not None and FIXTURE_SPEC.loader is not None
fixture_module = importlib.util.module_from_spec(FIXTURE_SPEC)
FIXTURE_SPEC.loader.exec_module(fixture_module)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class ControllerBudgetBuilderTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.temp = tempfile.TemporaryDirectory()
        cls.root = Path(cls.temp.name)
        original_rebuild = (
            fixture_module.FullControllerFixture.rebuild_storage_evidence
        )
        try:
            # This helper needs the resolver/bundle/source/partition fixture,
            # not the downstream storage certificate.  Keeping that unrelated
            # construction out also proves the budget builder has no circular
            # dependency on a pre-existing budget-derived certificate.
            fixture_module.FullControllerFixture.rebuild_storage_evidence = (
                lambda _self: None
            )
            cls.fixture = fixture_module.FullControllerFixture(
                cls.root / "fixture"
            )
        finally:
            fixture_module.FullControllerFixture.rebuild_storage_evidence = (
                original_rebuild
            )
        plan, plan_artifact = builder.load_pinned_json(
            cls.fixture.plan_path,
            sha(cls.fixture.plan_path),
            "fixture extraction plan",
        )
        receipt, receipt_artifact = builder.load_pinned_json(
            cls.fixture.receipt_path,
            sha(cls.fixture.receipt_path),
            "fixture preflight receipt",
        )
        cls.validation_context = builder.validate_input_chain(
            plan=plan,
            plan_path=cls.fixture.plan_path,
            plan_artifact=plan_artifact,
            receipt=receipt,
            receipt_artifact=receipt_artifact,
        )
        cls.validation_bindings = {
            "plan": cls.validation_context["plan_artifact"],
            "preflight_receipt": cls.validation_context[
                "preflight_receipt_artifact"
            ],
            "bundle_manifest": cls.validation_context["bundle_artifact"],
            "materialization_receipt": cls.validation_context[
                "materialization_artifact"
            ],
            "execution_partition_sha256": (
                builder.materializer.canonical_sha256(
                    plan["execution_partition"]
                )
            ),
            "partition_artifact_sha256": cls.validation_context[
                "partition_artifact"
            ]["sha256"],
        }

    @classmethod
    def tearDownClass(cls) -> None:
        cls.temp.cleanup()

    def invoke(self, name: str, *, plan: Path | None = None) -> subprocess.CompletedProcess:
        root = self.root / name
        root.mkdir()
        plan_path = plan or self.fixture.plan_path
        return subprocess.run(
            [
                sys.executable,
                str(BUILDER_PATH),
                "--plan",
                str(plan_path),
                "--plan-sha256",
                sha(plan_path),
                "--preflight-receipt",
                str(self.fixture.receipt_path),
                "--preflight-receipt-sha256",
                sha(self.fixture.receipt_path),
                "--future-storage-certificate",
                str(root / "future_storage_certificate.json"),
                "--output",
                str(root / "controller_budget.json"),
                "--derivation-output",
                str(root / "controller_budget_derivation.json"),
            ],
            text=True,
            capture_output=True,
            check=False,
        )

    def test_exact_serializers_produce_consumable_auditable_budget(self) -> None:
        result = self.invoke("pass")
        self.assertEqual(result.returncode, 0, result.stderr)
        summary = json.loads(result.stdout)
        root = self.root / "pass"
        budget_path = root / "controller_budget.json"
        derivation_path = root / "controller_budget_derivation.json"
        budget = json.loads(budget_path.read_text())
        derivation = json.loads(derivation_path.read_text())
        self.assertEqual(budget["schema"], builder.BUDGET_SCHEMA)
        self.assertEqual(set(budget), {
            "schema",
            "status",
            "submission_performed",
            "expected_job_count",
            "storage_budget",
            "full_training_authority",
            "full_extraction_authority",
        })
        normalized, blockers = builder.projector.validate_controller_budget(
            {"path": str(budget_path), "sha256": sha(budget_path)}
        )
        self.assertEqual(blockers, [])
        self.assertIsNotNone(normalized)
        assert normalized is not None
        self.assertEqual(
            normalized["projected_bytes"],
            summary["projected_bytes"],
        )
        self.assertEqual(
            normalized["projected_inodes"],
            summary["projected_inodes"],
        )
        self.assertEqual(
            derivation["budget_receipt"]["sha256"],
            sha(budget_path),
        )
        expected_semantic = derivation.pop("derivation_semantic_sha256")
        self.assertEqual(
            expected_semantic,
            builder.projector.semantic_sha256(derivation),
        )
        measurements = derivation["measurements"]
        self.assertEqual(
            measurements["job_record_count"],
            builder.materializer.EXPECTED_JOB_COUNT,
        )
        self.assertEqual(
            measurements["row_record_count"],
            builder.materializer.EXPECTED_ROW_COUNT,
        )
        self.assertLessEqual(
            measurements["exact_envelope_total_bytes"],
            derivation["budget_derivation"]["projected_bytes"],
        )
        self.assertFalse(summary["submission_performed"])
        derivation_artifact, derivation_blockers = (
            builder.projector.validate_controller_derivation(
                {
                    "path": str(derivation_path),
                    "sha256": sha(derivation_path),
                },
                controller=normalized,
                plan=self.validation_context["plan"],
                bindings=self.validation_bindings,
            )
        )
        self.assertEqual(derivation_blockers, [])
        self.assertEqual(
            derivation_artifact,
            builder.artifact_record(derivation_path),
        )

    def test_resealed_fabricated_low_measurement_is_rejected(self) -> None:
        result = self.invoke("fabricated_low_source")
        self.assertEqual(result.returncode, 0, result.stderr)
        source_root = self.root / "fabricated_low_source"
        derivation = json.loads(
            (source_root / "controller_budget_derivation.json").read_text()
        )

        target_root = self.root / "fabricated_low"
        target_root.mkdir()
        budget_path = target_root / "controller_budget.json"
        budget = {
            "schema": builder.BUDGET_SCHEMA,
            "status": "PASS",
            "submission_performed": False,
            "expected_job_count": builder.materializer.EXPECTED_JOB_COUNT,
            "storage_budget": {
                "fixed_bytes": 3,
                "fixed_inodes": 4,
                "bytes_per_job": 1,
                "inodes_per_job": 1,
            },
            "full_training_authority": 0,
            "full_extraction_authority": False,
        }
        budget_path.write_bytes(builder.canonical_json_bytes(budget))
        projected_bytes = 3 + builder.materializer.EXPECTED_JOB_COUNT
        projected_inodes = 4 + builder.materializer.EXPECTED_JOB_COUNT
        derivation["measurements"].update(
            {
                "row_manifest_bytes": 1,
                "job_manifest_exact_bytes": (
                    builder.materializer.EXPECTED_JOB_COUNT
                ),
                "job_record_min_bytes": 1,
                "job_record_max_bytes": 1,
                "job_record_ceiling_bytes": (
                    builder.materializer.EXPECTED_JOB_COUNT
                ),
                "submit_description_bytes": 1,
                "manifest_envelope_bytes": 1,
                "exact_envelope_total_bytes": projected_bytes,
                "exact_envelope_total_inodes": 4,
                "fixed_point_iterations": 1,
            }
        )
        for key, value in (
            ("fixed_bytes", 3),
            ("fixed_inodes", 4),
            ("bytes_per_job", 1),
            ("inodes_per_job", 1),
        ):
            derivation["budget_derivation"][key]["value"] = value
        derivation["budget_derivation"]["projected_bytes"] = projected_bytes
        derivation["budget_derivation"]["projected_inodes"] = projected_inodes
        derivation["future_bindings"]["budget_output_path"] = str(budget_path)
        derivation["budget_receipt"] = builder.artifact_record(budget_path)
        derivation.pop("derivation_semantic_sha256", None)
        derivation["derivation_semantic_sha256"] = (
            builder.projector.semantic_sha256(derivation)
        )
        derivation_path = target_root / "controller_budget_derivation.json"
        derivation_path.write_bytes(builder.canonical_json_bytes(derivation))

        normalized, blockers = builder.projector.validate_controller_budget(
            {"path": str(budget_path), "sha256": sha(budget_path)}
        )
        self.assertEqual(blockers, [])
        self.assertIsNotNone(normalized)
        with self.assertRaisesRegex(
            builder.projector.ProjectionError,
            "exact serializer replay",
        ):
            builder.projector.validate_controller_derivation(
                {
                    "path": str(derivation_path),
                    "sha256": sha(derivation_path),
                },
                controller=normalized,
                plan=self.validation_context["plan"],
                bindings=self.validation_bindings,
            )

    def test_partition_hash_mutation_is_rejected(self) -> None:
        original = self.fixture.partition_path.read_bytes()
        try:
            self.fixture.partition_path.write_bytes(original + b"\n")
            result = self.invoke("partition_mutation")
        finally:
            self.fixture.partition_path.write_bytes(original)
        self.assertEqual(result.returncode, 2)
        self.assertIn("partition", result.stderr.lower())
        self.assertFalse(
            (self.root / "partition_mutation" / "controller_budget.json").exists()
        )

    def test_plan_count_mutation_is_rejected(self) -> None:
        root = self.root / "mutated_plan_input"
        root.mkdir()
        plan = copy.deepcopy(self.fixture.plan)
        plan["rows"][0]["input_contract"]["expected_job_count"] += 1
        mutated = root / "plan.json"
        mutated.write_bytes(builder.canonical_json_bytes(plan))
        result = self.invoke("plan_count_mutation", plan=mutated)
        self.assertEqual(result.returncode, 2)
        self.assertIn("plan", result.stderr.lower())
        self.assertFalse(
            (self.root / "plan_count_mutation" / "controller_budget.json").exists()
        )

    def test_existing_output_is_never_overwritten(self) -> None:
        root = self.root / "existing_output"
        root.mkdir()
        existing = root / "controller_budget.json"
        existing.write_text("preserve\n")
        result = subprocess.run(
            [
                sys.executable,
                str(BUILDER_PATH),
                "--plan",
                str(self.fixture.plan_path),
                "--plan-sha256",
                sha(self.fixture.plan_path),
                "--preflight-receipt",
                str(self.fixture.receipt_path),
                "--preflight-receipt-sha256",
                sha(self.fixture.receipt_path),
                "--future-storage-certificate",
                str(root / "future_storage_certificate.json"),
                "--output",
                str(existing),
                "--derivation-output",
                str(root / "derivation.json"),
            ],
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(result.returncode, 2)
        self.assertEqual(existing.read_text(), "preserve\n")
        self.assertFalse((root / "derivation.json").exists())


if __name__ == "__main__":
    unittest.main()
