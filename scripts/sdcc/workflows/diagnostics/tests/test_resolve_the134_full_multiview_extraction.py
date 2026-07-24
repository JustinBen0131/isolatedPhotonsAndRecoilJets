#!/usr/bin/env python3
"""Focused contract tests for the non-submitting THE-134 full resolver."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


HERE = Path(__file__).resolve()
CONTROLLER = HERE.parents[1] / "resolve_the134_full_multiview_extraction.py"
SPEC = importlib.util.spec_from_file_location("the134_full_resolver", CONTROLLER)
assert SPEC is not None and SPEC.loader is not None
resolver = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(resolver)


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_sha256(payload: object) -> str:
    return hashlib.sha256(
        json.dumps(
            payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True
        ).encode()
    ).hexdigest()


def values_for_key(payload: object, expected_key: str) -> list[object]:
    values: list[object] = []
    if isinstance(payload, dict):
        for key, value in payload.items():
            if key == expected_key:
                values.append(value)
            values.extend(values_for_key(value, expected_key))
    elif isinstance(payload, list):
        for value in payload:
            values.extend(values_for_key(value, expected_key))
    return values


class ResolverFixture:
    def __init__(self, root: Path):
        self.root = root
        self.root.mkdir(parents=True)
        self.artifacts = root / "artifacts"
        self.lists = root / "lists"
        self.artifacts.mkdir()
        self.lists.mkdir()
        self.release_lib = self.artifacts / "release_lib"
        self.release_lib64 = self.artifacts / "release_lib64"
        self.release_lib.mkdir()
        self.release_lib64.mkdir()
        self.bundle = root / "bundle.json"
        self.sources = root / "sources.json"
        self.write_bundle()
        self.write_sources()

    def write_bundle(self) -> None:
        records = []
        for index, role in enumerate(resolver.REQUIRED_ARTIFACT_ROLES, start=1):
            path = self.artifacts / f"{role}.artifact"
            path.write_text(f"{role}|fixture|{index}\n", encoding="utf-8")
            if role in {"submitter", "pp_executor", "auau_executor"}:
                path.chmod(0o755)
            records.append(
                {
                    "role": role,
                    "path": str(path),
                    "sha256": sha256_file(path),
                    "size_bytes": path.stat().st_size,
                }
            )
        payload = {
            "schema": resolver.BUNDLE_SCHEMA,
            "status": "PASS",
            "public_commit": "a" * 40,
            "code_sha256": "b" * 64,
            "replay_schema_sha256": "c" * 64,
            "training_schema_sha256": "d" * 64,
            "semantic_sha256": "e" * 64,
            "runtime": {
                "release": "ana.560",
                "offline_main": "/cvmfs/sphenix/example/ana.560",
                "calo_reco_soname": "libcalo_reco.so.0",
                "request_memory_mb": 8000,
                "release_core_lib_dir": str(self.release_lib),
                "release_core_lib64_dir": str(self.release_lib64),
            },
            "artifacts": records,
        }
        self.bundle.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )

    def write_sources(self) -> None:
        rows = []
        for row in resolver.inventory_rows():
            row_dir = self.lists / row["sample"]
            row_dir.mkdir()
            input_lists = []
            lines_by_role = {}
            for role in resolver.LIST_ROLES:
                lines = [
                    f"/inputs/{row['sample']}/{role}/file_000.root",
                    "# aligned comment",
                    f"/inputs/{row['sample']}/{role}/file_001.root",
                ]
                path = row_dir / resolver.LIST_FILENAMES[role]
                path.write_text("\n".join(lines) + "\n", encoding="utf-8")
                lines_by_role[role] = lines
                input_lists.append(
                    {
                        "role": role,
                        "path": str(path),
                        "sha256": sha256_file(path),
                        "line_count": 3,
                        "executable_count": 2,
                    }
                )
            tuples = []
            for physical_index in (0, 2):
                tuples.append(
                    {
                        "tuple_index": len(tuples),
                        "physical_line": physical_index + 1,
                        "inputs": {
                            role: lines_by_role[role][physical_index]
                            for role in resolver.LIST_ROLES
                        },
                    }
                )
            tuple_records_sha = canonical_sha256(tuples)
            source_semantic_payload = {
                "row_id": row["row_id"],
                "system": row["system"],
                "sample": row["sample"],
                "lists": [
                    {
                        "role": record["role"],
                        "sha256": record["sha256"],
                        "line_count": record["line_count"],
                        "executable_count": record["executable_count"],
                    }
                    for record in input_lists
                ],
                "tuple_count": 2,
                "tuple_records_sha256": tuple_records_sha,
            }
            rows.append(
                {
                    "row_id": row["row_id"],
                    "system": row["system"],
                    "sample": row["sample"],
                    "input_lists": input_lists,
                    "expected_input_count": 2,
                    "expected_occurrence_count": 2,
                    "tuple_records_sha256": tuple_records_sha,
                    "full_source_manifest_sha256": canonical_sha256(
                        source_semantic_payload
                    ),
                }
            )
        payload = {
            "schema": resolver.SOURCE_SCHEMA,
            "status": "PASS",
            "rows": rows,
        }
        self.sources.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )

    def command(self, out_dir: Path, *, period: str = "0mrad") -> list[str]:
        return [
            sys.executable,
            str(CONTROLLER),
            "preflight",
            "--bundle-manifest",
            str(self.bundle),
            "--bundle-sha256",
            sha256_file(self.bundle),
            "--source-manifest",
            str(self.sources),
            "--source-manifest-sha256",
            sha256_file(self.sources),
            "--tag",
            "the134_full_extraction_fixture_v1",
            "--pp-period",
            period,
            "--output-root",
            "/remote/output/the134_full_extraction_fixture_v1",
            "--evidence-root",
            "/remote/evidence/the134_full_extraction_fixture_v1",
            "--submit-root",
            "/remote/submit/the134_full_extraction_fixture_v1",
            "--out-dir",
            str(out_dir),
        ]


class TestFullExtractionResolver(unittest.TestCase):
    def run_command(
        self, arguments: list[str], *, expected_returncode: int = 0
    ) -> subprocess.CompletedProcess[str]:
        result = subprocess.run(
            arguments,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        self.assertEqual(
            result.returncode,
            expected_returncode,
            msg=f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}",
        )
        return result

    def test_inventory_is_exact_source_family(self) -> None:
        result = self.run_command(
            [sys.executable, str(CONTROLLER), "inventory", "--format", "json"]
        )
        payload = json.loads(result.stdout)
        self.assertEqual(payload["row_count"], 13)
        rows = payload["rows"]
        self.assertEqual(len({row["row_id"] for row in rows}), 13)
        self.assertEqual(len({row["sample"] for row in rows}), 13)
        samples = {row["sample"] for row in rows}
        self.assertEqual(
            samples,
            {
                "run28_photonjet5",
                "run28_photonjet10",
                "run28_photonjet20",
                "run28_jet8",
                "run28_jet12",
                "run28_jet20",
                "run28_jet30",
                "run28_embeddedPhoton12",
                "run28_embeddedPhoton20",
                "run28_embeddedJet12",
                "run28_embeddedJet20",
                "run28_embeddedJet30",
                "run28_embeddedJet40",
            },
        )
        self.assertNotIn("run28_jet40", samples)

    def test_preflight_is_deterministic_typed_and_non_submitting(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fixture = ResolverFixture(root / "fixture")
            first = root / "resolved_a"
            second = root / "resolved_b"
            first_result = self.run_command(fixture.command(first))
            self.run_command(fixture.command(second))

            command_result = json.loads(first_result.stdout)
            plan = json.loads((first / "the134_full_extraction_plan.json").read_text())
            receipt = json.loads((first / "preflight_receipt.json").read_text())
            rows = [
                json.loads(line)
                for line in (
                    first / "the134_full_extraction_rows.jsonl"
                ).read_text().splitlines()
            ]
            self.assertEqual(plan["status"], "PREFLIGHT_PASS")
            self.assertFalse(plan["submission_performed"])
            self.assertEqual(
                plan["execution_state"], "PREFLIGHT_ONLY_NO_CONDOR_MUTATION"
            )
            self.assertEqual(plan["authority"]["requested_scope"], "full")
            self.assertEqual(plan["authority"]["full_training_authority"], 0)
            self.assertEqual(
                plan["authority"]["authority_state"],
                "PREFLIGHT_RESOLVED_NOT_EARNED",
            )
            self.assertEqual(receipt["requested_scope"], "full")
            self.assertEqual(receipt["full_training_authority"], 0)
            self.assertEqual(
                receipt["authority_state"], "PREFLIGHT_RESOLVED_NOT_EARNED"
            )
            self.assertEqual(command_result["requested_scope"], "full")
            self.assertEqual(command_result["full_training_authority"], 0)
            self.assertEqual(
                command_result["authority_state"],
                "PREFLIGHT_RESOLVED_NOT_EARNED",
            )
            self.assertEqual(
                set(values_for_key(plan, "full_training_authority")), {0}
            )
            self.assertEqual(
                set(values_for_key(plan, "authority_state")),
                {"PREFLIGHT_RESOLVED_NOT_EARNED"},
            )
            self.assertFalse(receipt["submission_performed"])
            self.assertEqual(len(rows), 13)
            self.assertTrue(
                all(
                    row["requested_scope"] == "full"
                    and row["full_training_authority"] == 0
                    and row["authority_state"]
                    == "PREFLIGHT_RESOLVED_NOT_EARNED"
                    for row in rows
                )
            )
            self.assertTrue(
                all(
                    row["input_contract"]["expected_job_count"] == 2
                    and row["input_contract"]["tuple_count"] == 2
                    for row in rows
                )
            )
            self.assertTrue(
                all(
                    row["execution_contract"]["existing_submitter_argv"][2]
                    == "condorDoAll"
                    for row in rows
                )
            )
            self.assertTrue(
                all(
                    "$(Cluster).$(Process).root"
                    in row["execution_contract"]["multiview_sidecar_template"]
                    for row in rows
                )
            )
            pp_rows = [row for row in rows if row["system"] == "pp"]
            auau_rows = [row for row in rows if row["system"] == "auau"]
            training_contract = plan["training_period_si_contract"]
            self.assertEqual(
                training_contract,
                {
                    "schema": "THE134_PP_TRAINING_PERIOD_SI_CONTRACT_V1",
                    "system": "pp",
                    "period": "0mrad",
                    "si_di_role": "SI",
                    "row_scope": "ALL_PP_TRAINING_ROWS",
                },
            )
            training_contract_sha256 = canonical_sha256(training_contract)
            self.assertEqual(
                plan["training_period_si_contract_sha256"],
                training_contract_sha256,
            )
            self.assertEqual(
                receipt["training_period_si_contract_sha256"],
                training_contract_sha256,
            )
            self.assertTrue(
                all(
                    row["source_period"] == "0mrad"
                    and row["source_si_di_role"] == "SI"
                    and row["training_period_si_contract_sha256"]
                    == training_contract_sha256
                    and row["execution_contract"]["worker_environment"][
                        "RJ_PPG12_PHOTON_YIELD"
                    ]
                    == "1"
                    and row["execution_contract"]["worker_environment"][
                        "RJ_PPG12_PERIOD"
                    ]
                    == "0mrad"
                    for row in pp_rows
                )
            )
            self.assertTrue(
                all(
                    row["source_period"] == "AUAU_RUN24"
                    and row["source_si_di_role"] == "EMBEDDED"
                    and row["training_period_si_contract_sha256"] is None
                    for row in auau_rows
                )
            )
            self.assertEqual(
                plan["closure_witness_boundary"],
                {
                    "schema": (
                        "THE134_FULL_TRAINING_CLOSURE_WITNESS_BOUNDARY_V1"
                    ),
                    "preflight_contains_closure_witness": False,
                    "sole_authority_earning_tool": (
                        "prepare_the134_h70_matrix.py"
                    ),
                    "required_argv_suffix": ["--scope", "full"],
                    "required_audit_schema": (
                        "THE134_FACTORIAL_VIEW_TRAINING_MATRIX_AUDIT_V1"
                    ),
                    "required_audit_status": "PASS",
                    "required_source_population_closure_status": "PASS",
                },
            )
            self.assertEqual(
                (
                    first
                    / "the134_full_extraction_duplicate_fingerprint.sha256"
                ).read_bytes(),
                (
                    second
                    / "the134_full_extraction_duplicate_fingerprint.sha256"
                ).read_bytes(),
            )
            self.assertEqual(
                (
                    first / "the134_full_extraction_rows.jsonl"
                ).read_bytes(),
                (
                    second / "the134_full_extraction_rows.jsonl"
                ).read_bytes(),
            )
            self.assertEqual(
                (
                    first / "the134_full_extraction_plan.json"
                ).read_bytes(),
                (
                    second / "the134_full_extraction_plan.json"
                ).read_bytes(),
            )

    def test_preflight_authority_mutations_are_rejected(self) -> None:
        valid = {
            "requested_scope": "full",
            "full_training_authority": 0,
            "authority_state": "PREFLIGHT_RESOLVED_NOT_EARNED",
        }
        resolver.validate_preflight_authority_payload(valid, label="fixture")
        mutations = (
            (
                "full_training_authority",
                1,
                "must remain full_training_authority=0 during preflight",
            ),
            (
                "authority_state",
                "FULL_TRAINING_AUTHORITY_EARNED",
                "must remain authority_state=PREFLIGHT_RESOLVED_NOT_EARNED",
            ),
            (
                "requested_scope",
                "smoke",
                "must remain requested_scope=full",
            ),
        )
        for key, value, message in mutations:
            with self.subTest(key=key, value=value):
                mutated = {"nested": dict(valid)}
                mutated["nested"][key] = value
                with self.assertRaisesRegex(resolver.ControllerError, message):
                    resolver.validate_preflight_authority_payload(
                        mutated, label="fixture"
                    )

    def test_explicit_pp_period_si_contract_controls_every_pp_row(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fixture = ResolverFixture(root / "fixture")
            out_dir = root / "resolved"
            self.run_command(fixture.command(out_dir, period="1p5mrad"))
            plan = json.loads(
                (out_dir / "the134_full_extraction_plan.json").read_text()
            )
            contract = plan["training_period_si_contract"]
            contract_sha256 = canonical_sha256(contract)
            self.assertEqual(contract["period"], "1p5mrad")
            self.assertEqual(contract["si_di_role"], "SI")
            pp_rows = [
                row for row in plan["rows"] if row["system"] == "pp"
            ]
            self.assertEqual(len(pp_rows), 7)
            self.assertTrue(
                all(
                    row["source_period"] == "1p5mrad"
                    and row["source_si_di_role"] == "SI"
                    and row["training_period_si_contract_sha256"]
                    == contract_sha256
                    and row["execution_contract"]["worker_environment"][
                        "RJ_PPG12_PERIOD"
                    ]
                    == "1p5mrad"
                    and row["execution_contract"]["worker_environment"][
                        "RJ_REPLAY_SI_DI_ROLE"
                    ]
                    == "SI"
                    for row in pp_rows
                )
            )

    def test_bundle_hash_drift_is_rejected_without_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fixture = ResolverFixture(root / "fixture")
            out_dir = root / "resolved"
            command = fixture.command(out_dir)
            command[command.index("--bundle-sha256") + 1] = "0" * 64
            result = self.run_command(command, expected_returncode=2)
            self.assertIn("immutable bundle manifest hash drift", result.stderr)
            self.assertFalse(out_dir.exists())

    def test_bundle_artifact_mutation_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fixture = ResolverFixture(root / "fixture")
            payload = json.loads(fixture.bundle.read_text())
            mutated = Path(payload["artifacts"][0]["path"])
            mutated.write_text(mutated.read_text() + "drift\n")
            out_dir = root / "resolved"
            result = self.run_command(
                fixture.command(out_dir), expected_returncode=2
            )
            self.assertIn("bundle artifact hash drift", result.stderr)
            self.assertFalse(out_dir.exists())

    def test_source_row_omission_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fixture = ResolverFixture(root / "fixture")
            payload = json.loads(fixture.sources.read_text())
            payload["rows"].pop()
            fixture.sources.write_text(
                json.dumps(payload, indent=2, sort_keys=True) + "\n"
            )
            out_dir = root / "resolved"
            result = self.run_command(
                fixture.command(out_dir), expected_returncode=2
            )
            self.assertIn("source-authority row closure differs", result.stderr)
            self.assertFalse(out_dir.exists())

    def test_source_list_mutation_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fixture = ResolverFixture(root / "fixture")
            payload = json.loads(fixture.sources.read_text())
            mutated = Path(payload["rows"][0]["input_lists"][0]["path"])
            mutated.write_text(mutated.read_text() + "/unexpected/root/file.root\n")
            out_dir = root / "resolved"
            result = self.run_command(
                fixture.command(out_dir), expected_returncode=2
            )
            self.assertIn("list hash drift", result.stderr)
            self.assertFalse(out_dir.exists())

    def test_partial_five_file_tuple_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fixture = ResolverFixture(root / "fixture")
            payload = json.loads(fixture.sources.read_text())
            record = payload["rows"][0]["input_lists"][0]
            path = Path(record["path"])
            lines = path.read_text().splitlines()
            lines[2] = "# one list becomes inactive"
            path.write_text("\n".join(lines) + "\n")
            record["sha256"] = sha256_file(path)
            record["executable_count"] = 1
            fixture.sources.write_text(
                json.dumps(payload, indent=2, sort_keys=True) + "\n"
            )
            out_dir = root / "resolved"
            result = self.run_command(
                fixture.command(out_dir), expected_returncode=2
            )
            self.assertIn("partial five-file tuple", result.stderr)
            self.assertFalse(out_dir.exists())

    def test_cross_source_duplicate_tuple_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fixture = ResolverFixture(root / "fixture")
            payload = json.loads(fixture.sources.read_text())
            first = payload["rows"][0]
            second = payload["rows"][1]
            for first_record, second_record in zip(
                first["input_lists"], second["input_lists"]
            ):
                first_path = Path(first_record["path"])
                second_path = Path(second_record["path"])
                second_path.write_bytes(first_path.read_bytes())
                second_record["sha256"] = sha256_file(second_path)
                second_record["line_count"] = first_record["line_count"]
                second_record["executable_count"] = first_record[
                    "executable_count"
                ]
            lines_by_role = {
                record["role"]: Path(record["path"]).read_text().splitlines()
                for record in second["input_lists"]
            }
            tuples = []
            for physical_index in (0, 2):
                tuples.append(
                    {
                        "tuple_index": len(tuples),
                        "physical_line": physical_index + 1,
                        "inputs": {
                            role: lines_by_role[role][physical_index]
                            for role in resolver.LIST_ROLES
                        },
                    }
                )
            second["tuple_records_sha256"] = canonical_sha256(tuples)
            second_semantic = {
                "row_id": second["row_id"],
                "system": second["system"],
                "sample": second["sample"],
                "lists": [
                    {
                        "role": record["role"],
                        "sha256": record["sha256"],
                        "line_count": record["line_count"],
                        "executable_count": record["executable_count"],
                    }
                    for record in second["input_lists"]
                ],
                "tuple_count": 2,
                "tuple_records_sha256": second["tuple_records_sha256"],
            }
            second["full_source_manifest_sha256"] = canonical_sha256(
                second_semantic
            )
            fixture.sources.write_text(
                json.dumps(payload, indent=2, sort_keys=True) + "\n"
            )
            out_dir = root / "resolved"
            result = self.run_command(
                fixture.command(out_dir), expected_returncode=2
            )
            self.assertIn("duplicate five-file source tuples", result.stderr)
            self.assertFalse(out_dir.exists())

    def test_cli_has_no_submit_action(self) -> None:
        result = self.run_command(
            [sys.executable, str(CONTROLLER), "submit"],
            expected_returncode=2,
        )
        self.assertIn("invalid choice", result.stderr)
        source = CONTROLLER.read_text(encoding="utf-8")
        self.assertNotIn("import subprocess", source)
        self.assertNotIn("condor_submit", source)
        self.assertNotIn("os.system(", source)

    def test_preflight_refuses_to_overwrite_existing_evidence(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fixture = ResolverFixture(root / "fixture")
            out_dir = root / "resolved"
            self.run_command(fixture.command(out_dir))
            result = self.run_command(
                fixture.command(out_dir), expected_returncode=2
            )
            self.assertIn("preflight output directory already exists", result.stderr)


if __name__ == "__main__":
    unittest.main()
