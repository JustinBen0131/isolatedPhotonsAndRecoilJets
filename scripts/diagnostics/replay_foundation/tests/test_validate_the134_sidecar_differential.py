#!/usr/bin/env python3
"""Focused unit tests for the THE-134 sidecar differential comparator."""

from __future__ import annotations

import copy
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import uproot


HERE = Path(__file__).resolve()
VALIDATOR_PATH = HERE.parents[1] / "validate_the134_sidecar_differential.py"
SPEC = importlib.util.spec_from_file_location(
    "validate_the134_sidecar_differential", VALIDATOR_PATH
)
assert SPEC is not None and SPEC.loader is not None
validator = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(validator)


class SidecarDifferentialComparatorTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)

    def tearDown(self) -> None:
        self.temporary.cleanup()

    @staticmethod
    def write_histogram(path: Path, values: np.ndarray) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        with uproot.recreate(path) as output:
            output["qa/hExact"] = (
                values.astype(np.float64),
                np.asarray([0.0, 1.0, 2.0, 3.0], dtype=np.float64),
            )

    @staticmethod
    def write_preflight(
        path: Path,
        rows: tuple[dict[str, str], ...],
        *,
        code_sha256: str = "1" * 64,
        schema_sha256: str = "2" * 64,
        semantic_sha256: str = "3" * 64,
        tag: str = "the134_fixture",
        base: Path | None = None,
    ) -> dict[str, str]:
        resolved_base = (
            base.resolve()
            if base is not None
            else (path.parent / f"{tag}-output").resolve()
        )
        exact_keys = ",".join(
            f"{arm}:{row['lane']}:{row['sample']}"
            for row in rows
            for arm in ("direct", "writer")
        )
        fields = {
            "tag": tag,
            "base": str(resolved_base),
            "code_sha256": code_sha256,
            "schema_sha": schema_sha256,
            "semantic_sha": semantic_sha256,
            "only_keys": exact_keys,
            "source_sha_override": "",
            "extra_pp_template": (
                "RJ_PP_PHOTONID_EXTRACT_ONLY=1"
                if any(row["system"] == "pp" for row in rows)
                else ""
            ),
            "extra_auau_template": (
                "RJ_AUAU_BDT_EXTRACT_ONLY=1"
                if any(row["system"] == "auau" for row in rows)
                else ""
            ),
            "writer_extra_common_template": (
                "RJ_THE134_MULTIVIEW_SIDECAR_ONLY_V1=1"
            ),
            "writer_extra_pp_template": "",
            "writer_extra_auau_template": "",
        }
        path.write_text(
            "".join(f"{key}={value}\n" for key, value in fields.items()),
            encoding="utf-8",
        )
        return fields

    def write_terminal_gate(
        self,
        path: Path,
        *,
        preflight: Path,
        output_root: Path,
        tag: str,
    ) -> dict[str, object]:
        initial_queue = path.parent / "initial_queue.tsv"
        initial_queue.write_text(
            (
                f"1001 0 2 --output {output_root.resolve()}/"
                "direct/pp_inclusive_sim/run28_jet8\n"
                f"1002 0 2 --output {output_root.resolve()}/"
                "writer/pp_inclusive_sim/run28_jet8\n"
            ),
            encoding="utf-8",
        )
        payload = {
            "schema": validator.TERMINAL_GATE_SCHEMA,
            "tag": tag,
            "output_root": str(output_root.resolve()),
            "preflight_receipt": str(preflight.resolve()),
            "preflight_receipt_sha256": validator.sha256_file(preflight),
            "initial_queue_tsv": str(initial_queue.resolve()),
            "initial_queue_tsv_sha256": validator.sha256_file(initial_queue),
            "row_count": 2,
            "rows": [
                {
                    "cluster_id": 1001,
                    "proc_id": 0,
                    "cluster_proc": "1001.0",
                    "role": "direct",
                    "lane": "pp_inclusive_sim",
                    "sample": "run28_jet8",
                    "job_status": 4,
                    "exit_code": 0,
                },
                {
                    "cluster_id": 1002,
                    "proc_id": 0,
                    "cluster_proc": "1002.0",
                    "role": "writer",
                    "lane": "pp_inclusive_sim",
                    "sample": "run28_jet8",
                    "job_status": 4,
                    "exit_code": 0,
                },
            ],
        }
        path.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        return payload

    def make_component_payload(
        self,
        *,
        system: str,
        receipt: Path,
        ownership_rows: tuple[dict[str, str], ...],
        label: str,
    ) -> dict[str, object]:
        parsed = validator.parse_preflight(receipt, ownership_rows)
        output_root = Path(parsed["fields"]["base"])
        expected_row = (
            validator.PP_ROWS[0]
            if system == "pp"
            else validator.AUAU_ROWS[0]
        )
        direct_analysis = (
            output_root
            / "direct"
            / expected_row["lane"]
            / expected_row["sample"]
            / f"{label}-analysis.root"
        )
        writer_analysis = (
            output_root
            / "writer"
            / expected_row["lane"]
            / expected_row["sample"]
            / f"{label}-analysis.root"
        )
        sidecar = writer_analysis.parent / validator.SIDECAR_NAME
        direct_analysis.parent.mkdir(parents=True, exist_ok=True)
        writer_analysis.parent.mkdir(parents=True, exist_ok=True)
        direct_analysis.write_bytes(
            b"A" * (validator.ANALYSIS_MIN_BYTES + 1)
        )
        writer_analysis.write_bytes(
            b"A" * (validator.ANALYSIS_MIN_BYTES + 1)
        )
        sidecar.write_bytes(
            b"S" * (validator.SIDECAR_GROSS_TRUNCATION_BYTES + 1)
        )
        payload = {
            "schema": validator.COMPONENT_CERTIFICATE_SCHEMA,
            "status": "PASS",
            "artifact_profile": validator.PROFILE,
            "component_system": system,
            "output_root": str(output_root),
            "preflight_receipt": str(receipt),
            "preflight_receipt_sha256": validator.sha256_file(receipt),
            "receipt_ownership_systems": sorted(
                {row["system"] for row in ownership_rows}
            ),
            "provenance": validator.preflight_provenance(parsed),
            "full_training_authority": 0,
            "full_extraction_authority": False,
            "broad_production_authority": False,
            "rows": [
                {
                    "system": system,
                    "lane": expected_row["lane"],
                    "sample": expected_row["sample"],
                    "status": "PASS",
                    "analysis_health": [
                        {
                            "path": str(direct_analysis),
                            "bytes": direct_analysis.stat().st_size,
                            "sha256": validator.sha256_file(direct_analysis),
                            "readable": True,
                            "zombie": False,
                            "recovered": False,
                        },
                        {
                            "path": str(writer_analysis),
                            "bytes": writer_analysis.stat().st_size,
                            "sha256": validator.sha256_file(writer_analysis),
                            "readable": True,
                            "zombie": False,
                            "recovered": False,
                        }
                    ],
                    "sidecar_health": {
                        "path": str(sidecar),
                        "bytes": sidecar.stat().st_size,
                        "sha256": validator.sha256_file(sidecar),
                        "readable": True,
                        "zombie": False,
                        "recovered": False,
                    },
                    "failures": [],
                }
            ],
            "failures": [],
        }
        if system == "pp":
            terminal_path = receipt.parent / f"{label}-terminal-gate.json"
            self.write_terminal_gate(
                terminal_path,
                preflight=receipt,
                output_root=output_root,
                tag=parsed["fields"]["tag"],
            )
            payload["terminal_gate"] = validator.parse_terminal_gate_receipt(
                terminal_path,
                preflight=parsed,
                output_root=output_root,
            )
        return payload

    def test_exact_histogram_neutrality_and_mutation_rejection(self) -> None:
        direct = self.root / "direct.root"
        writer = self.root / "writer.root"
        values = np.asarray([1.0, 4.0, 9.0], dtype=np.float64)
        self.write_histogram(direct, values)
        self.write_histogram(writer, values)
        report = validator.compare_histograms(direct, writer)
        self.assertEqual(report["status"], "PASS")
        self.assertEqual(report["histograms"], 1)

        self.write_histogram(writer, np.asarray([1.0, 4.0, 10.0]))
        report = validator.compare_histograms(direct, writer)
        self.assertEqual(report["status"], "FAIL")
        self.assertTrue(any("values" in item for item in report["failures"]))

    def test_source_identity_is_deterministic_and_row_specific(self) -> None:
        tag = "the134_fixture"
        pp = validator.source_sha256(tag, dict(validator.ROWS[0]))
        auau = validator.source_sha256(tag, dict(validator.ROWS[1]))
        self.assertRegex(pp, r"^[0-9a-f]{64}$")
        self.assertRegex(auau, r"^[0-9a-f]{64}$")
        self.assertNotEqual(pp, auau)
        self.assertEqual(
            pp, validator.source_sha256(tag, dict(validator.ROWS[0]))
        )

    def test_preflight_requires_frozen_direct_and_writer_profiles(self) -> None:
        receipt = self.root / "preflight_receipt.txt"
        exact_keys = ",".join(
            f"{arm}:{row['lane']}:{row['sample']}"
            for row in validator.ROWS
            for arm in ("direct", "writer")
        )
        fields = {
            "tag": "the134_fixture",
            "base": str((self.root / "full-output").resolve()),
            "code_sha256": "1" * 64,
            "schema_sha": "2" * 64,
            "semantic_sha": "3" * 64,
            "only_keys": exact_keys,
            "extra_pp_template": "RJ_PP_PHOTONID_EXTRACT_ONLY=1",
            "extra_auau_template": "RJ_AUAU_BDT_EXTRACT_ONLY=1",
            "writer_extra_common_template": (
                "RJ_THE134_MULTIVIEW_SIDECAR_ONLY_V1=1"
            ),
            "writer_extra_pp_template": "",
            "writer_extra_auau_template": "",
        }
        receipt.write_text(
            "".join(f"{key}={value}\n" for key, value in fields.items()),
            encoding="utf-8",
        )
        parsed = validator.parse_preflight(receipt)
        self.assertEqual(parsed["fields"]["tag"], "the134_fixture")

        fields["extra_pp_template"] = ""
        receipt.write_text(
            "".join(f"{key}={value}\n" for key, value in fields.items()),
            encoding="utf-8",
        )
        with self.assertRaises(validator.ValidationError):
            validator.parse_preflight(receipt)

    def test_pp_component_preflight_rejects_missing_or_duplicate_keys(self) -> None:
        receipt = self.root / "pp_preflight_receipt.txt"
        fields = self.write_preflight(receipt, validator.PP_ROWS)
        parsed = validator.parse_preflight(receipt, validator.PP_ROWS)
        self.assertEqual(
            set(parsed["fields"]["only_keys"].split(",")),
            {
                "direct:pp_inclusive_sim:run28_jet8",
                "writer:pp_inclusive_sim:run28_jet8",
            },
        )

        fields["only_keys"] += ",direct:pp_inclusive_sim:run28_jet8"
        receipt.write_text(
            "".join(f"{key}={value}\n" for key, value in fields.items()),
            encoding="utf-8",
        )
        with self.assertRaises(validator.ValidationError):
            validator.parse_preflight(receipt, validator.PP_ROWS)

        fields["only_keys"] = "direct:pp_inclusive_sim:run28_jet8"
        receipt.write_text(
            "".join(f"{key}={value}\n" for key, value in fields.items()),
            encoding="utf-8",
        )
        with self.assertRaises(validator.ValidationError):
            validator.parse_preflight(receipt, validator.PP_ROWS)

    def test_preflight_rejects_duplicate_fields_and_pinned_basenames(self) -> None:
        receipt = self.root / "duplicate_preflight_receipt.txt"
        fields = self.write_preflight(receipt, validator.PP_ROWS)
        receipt.write_text(
            receipt.read_text(encoding="utf-8")
            + f"only_keys={fields['only_keys']}\n",
            encoding="utf-8",
        )
        with self.assertRaises(validator.ValidationError):
            validator.parse_preflight(receipt, validator.PP_ROWS)

        self.write_preflight(receipt, validator.PP_ROWS)
        receipt.write_text(
            receipt.read_text(encoding="utf-8")
            + f"{'a' * 64} /one/frozen-model.xml\n"
            + f"{'b' * 64} /two/frozen-model.xml\n",
            encoding="utf-8",
        )
        with self.assertRaises(validator.ValidationError):
            validator.parse_preflight(receipt, validator.PP_ROWS)

    def test_single_receipt_and_pp_component_certification_compatibility(
        self,
    ) -> None:
        full_receipt = self.root / "full_preflight_receipt.txt"
        full_fields = self.write_preflight(full_receipt, validator.ROWS)
        full_rows = [
            {"system": "pp", "status": "PASS", "failures": []},
            {"system": "auau", "status": "PASS", "failures": []},
        ]
        full_args = SimpleNamespace(
            output_root=Path(full_fields["base"]),
            preflight_receipt=full_receipt,
            component_system=None,
            pp_component_certificate=None,
            auau_output_root=None,
            auau_preflight_receipt=None,
        )
        with (
            patch.object(validator, "load_module", return_value=object()),
            patch.object(validator, "validate_rows", return_value=full_rows),
        ):
            full = validator.validate(full_args)
        self.assertEqual(full["schema"], validator.CERTIFICATE_SCHEMA)
        self.assertEqual(full["status"], "PASS")
        self.assertNotIn("component_system", full)
        self.assertEqual(
            full["receipt_ownership_systems"], ["auau", "pp"]
        )

        pp_receipt = self.root / "pp_preflight_receipt.txt"
        pp_fields = self.write_preflight(
            pp_receipt,
            validator.PP_ROWS,
            code_sha256="4" * 64,
            tag="the134_pp_replacement",
        )
        pp_terminal = self.root / "pp_terminal_gate.json"
        self.write_terminal_gate(
            pp_terminal,
            preflight=pp_receipt,
            output_root=Path(pp_fields["base"]),
            tag=pp_fields["tag"],
        )
        pp_args = SimpleNamespace(
            output_root=Path(pp_fields["base"]),
            preflight_receipt=pp_receipt,
            component_system="pp",
            terminal_gate_receipt=pp_terminal,
            pp_component_certificate=None,
            auau_output_root=None,
            auau_preflight_receipt=None,
        )
        with (
            patch.object(validator, "load_module", return_value=object()),
            patch.object(
                validator,
                "validate_rows",
                return_value=[
                    {"system": "pp", "status": "PASS", "failures": []}
                ],
            ),
        ):
            component = validator.validate(pp_args)
        self.assertEqual(
            component["schema"], validator.COMPONENT_CERTIFICATE_SCHEMA
        )
        self.assertEqual(component["component_system"], "pp")
        self.assertEqual(component["receipt_ownership_systems"], ["pp"])
        self.assertEqual(component["provenance"]["code_sha256"], "4" * 64)

    def test_aggregate_binds_distinct_honest_receipts_and_rejects_drift(
        self,
    ) -> None:
        pp_receipt = self.root / "pp_preflight_receipt.txt"
        auau_receipt = self.root / "preserved_full_preflight_receipt.txt"
        self.write_preflight(
            pp_receipt,
            validator.PP_ROWS,
            code_sha256="4" * 64,
            tag="the134_pp_replacement",
        )
        self.write_preflight(
            auau_receipt,
            validator.ROWS,
            code_sha256="5" * 64,
            tag="the134_preserved_four_row",
        )
        pp_payload = self.make_component_payload(
            system="pp",
            receipt=pp_receipt,
            ownership_rows=validator.PP_ROWS,
            label="pp",
        )
        auau_payload = self.make_component_payload(
            system="auau",
            receipt=auau_receipt,
            ownership_rows=validator.ROWS,
            label="auau",
        )
        pp_certificate = self.root / "pp_component_certificate.json"
        pp_certificate.write_text(
            json.dumps(pp_payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        with patch.object(
            validator,
            "certify_single_receipt",
            return_value=copy.deepcopy(pp_payload),
        ):
            pp_binding = validator.load_pp_component_certificate(pp_certificate)
        for field, value in (
            ("lane", "pp_photon_sim"),
            ("sample", "run28_photonjet20"),
        ):
            wrong_row = copy.deepcopy(pp_payload)
            wrong_row["rows"][0][field] = value
            pp_certificate.write_text(
                json.dumps(wrong_row, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
            with self.assertRaises(validator.ValidationError):
                validator.load_pp_component_certificate(pp_certificate)

        wrong_root = copy.deepcopy(pp_payload)
        unrelated_root = self.root / "unrelated-output"
        unrelated_root.mkdir()
        wrong_root["output_root"] = str(unrelated_root)
        pp_certificate.write_text(
            json.dumps(wrong_root, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        with (
            patch.object(
                validator,
                "certify_single_receipt",
                return_value=copy.deepcopy(pp_payload),
            ),
            self.assertRaises(validator.ValidationError),
        ):
            validator.load_pp_component_certificate(pp_certificate)

        outside_artifact = self.root / "outside-analysis.root"
        outside_artifact.write_bytes(
            b"X" * (validator.ANALYSIS_MIN_BYTES + 1)
        )
        outside_binding = copy.deepcopy(pp_payload)
        outside_health = outside_binding["rows"][0]["analysis_health"][0]
        outside_health["path"] = str(outside_artifact)
        outside_health["bytes"] = outside_artifact.stat().st_size
        outside_health["sha256"] = validator.sha256_file(outside_artifact)
        pp_certificate.write_text(
            json.dumps(outside_binding, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        with (
            patch.object(
                validator,
                "certify_single_receipt",
                return_value=copy.deepcopy(pp_payload),
            ),
            self.assertRaises(validator.ValidationError),
        ):
            validator.load_pp_component_certificate(pp_certificate)

        missing_sidecar = copy.deepcopy(pp_payload)
        missing_sidecar["rows"][0].pop("sidecar_health")
        pp_certificate.write_text(
            json.dumps(missing_sidecar, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        with (
            patch.object(
                validator,
                "certify_single_receipt",
                return_value=copy.deepcopy(pp_payload),
            ),
            self.assertRaises(validator.ValidationError),
        ):
            validator.load_pp_component_certificate(pp_certificate)

        with patch.object(
            validator,
            "certify_single_receipt",
            return_value=copy.deepcopy(pp_payload),
        ):
            aggregate = validator.bind_aggregate_certificate(
                pp_binding, auau_payload
            )
        self.assertEqual(
            aggregate["schema"], validator.AGGREGATE_CERTIFICATE_SCHEMA
        )
        self.assertEqual(aggregate["status"], "PASS")
        self.assertTrue(aggregate["independent_code_provenance"])
        self.assertEqual(
            aggregate["component_bindings"]["pp"]["provenance"][
                "code_sha256"
            ],
            "4" * 64,
        )
        self.assertEqual(
            aggregate["component_bindings"]["auau"]["provenance"][
                "code_sha256"
            ],
            "5" * 64,
        )

        dishonest_pp = copy.deepcopy(pp_payload)
        dishonest_pp["provenance"]["code_sha256"] = "6" * 64
        pp_certificate.write_text(
            json.dumps(dishonest_pp, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        with self.assertRaises(validator.ValidationError):
            validator.load_pp_component_certificate(pp_certificate)

        dishonest_auau = copy.deepcopy(auau_payload)
        dishonest_auau["provenance"]["code_sha256"] = "7" * 64
        with (
            patch.object(
                validator,
                "certify_single_receipt",
                return_value=copy.deepcopy(pp_payload),
            ),
            self.assertRaises(validator.ValidationError),
        ):
            validator.bind_aggregate_certificate(pp_binding, dishonest_auau)

        duplicate_auau = copy.deepcopy(auau_payload)
        duplicate_auau["rows"].append(copy.deepcopy(duplicate_auau["rows"][0]))
        with (
            patch.object(
                validator,
                "certify_single_receipt",
                return_value=copy.deepcopy(pp_payload),
            ),
            self.assertRaises(validator.ValidationError),
        ):
            validator.bind_aggregate_certificate(pp_binding, duplicate_auau)

        same_code_receipt = self.root / "same_code_full_preflight_receipt.txt"
        self.write_preflight(
            same_code_receipt,
            validator.ROWS,
            code_sha256="4" * 64,
            tag="the134_preserved_same_code_mutation",
        )
        same_code_auau = self.make_component_payload(
            system="auau",
            receipt=same_code_receipt,
            ownership_rows=validator.ROWS,
            label="same-code-auau",
        )
        with (
            patch.object(
                validator,
                "certify_single_receipt",
                return_value=copy.deepcopy(pp_payload),
            ),
            self.assertRaises(validator.ValidationError),
        ):
            validator.bind_aggregate_certificate(pp_binding, same_code_auau)

        for identity_name, schema_sha, semantic_sha in (
            ("schema", "8" * 64, "3" * 64),
            ("semantic", "2" * 64, "9" * 64),
        ):
            mismatch_receipt = (
                self.root / f"{identity_name}_mismatch_preflight.txt"
            )
            self.write_preflight(
                mismatch_receipt,
                validator.ROWS,
                code_sha256="5" * 64,
                schema_sha256=schema_sha,
                semantic_sha256=semantic_sha,
                tag=f"the134_preserved_{identity_name}_mismatch",
            )
            mismatch_auau = self.make_component_payload(
                system="auau",
                receipt=mismatch_receipt,
                ownership_rows=validator.ROWS,
                label=f"{identity_name}-mismatch-auau",
            )
            with (
                patch.object(
                    validator,
                    "certify_single_receipt",
                    return_value=copy.deepcopy(pp_payload),
                ),
                self.assertRaises(validator.ValidationError),
            ):
                validator.bind_aggregate_certificate(
                    pp_binding, mismatch_auau
                )

    def test_direct_sidecar_and_inventory_drift_rejected(self) -> None:
        row = dict(validator.ROWS[0])
        direct_base = (
            self.root / "direct" / row["lane"] / row["sample"]
        )
        writer_base = (
            self.root / "writer" / row["lane"] / row["sample"]
        )
        self.write_histogram(direct_base / "nested" / "analysis.root", np.ones(3))
        self.write_histogram(writer_base / "nested" / "analysis.root", np.ones(3))
        self.write_histogram(writer_base / validator.SIDECAR_NAME, np.ones(3))
        discovered = validator.discover_row_files(self.root, row)
        self.assertEqual(len(discovered["pairs"]), 1)

        self.write_histogram(direct_base / validator.SIDECAR_NAME, np.ones(3))
        with self.assertRaises(validator.ValidationError):
            validator.discover_row_files(self.root, row)

    def test_terminal_gate_rejects_identity_role_state_and_queue_drift(
        self,
    ) -> None:
        receipt = self.root / "pp_preflight_receipt.txt"
        fields = self.write_preflight(
            receipt,
            validator.PP_ROWS,
            tag="the134_pp_terminal_fixture",
        )
        terminal_path = self.root / "terminal_gate.json"
        valid = self.write_terminal_gate(
            terminal_path,
            preflight=receipt,
            output_root=Path(fields["base"]),
            tag=fields["tag"],
        )
        preflight = validator.parse_preflight(receipt, validator.PP_ROWS)
        parsed = validator.parse_terminal_gate_receipt(
            terminal_path,
            preflight=preflight,
            output_root=Path(fields["base"]),
        )
        self.assertEqual(parsed["roles"], ["direct", "writer"])

        exact_root = Path(fields["base"]).resolve()
        queue_path = self.root / "initial_queue.tsv"
        queue_mutations = (
            (
                f"1001 0 2 --output {exact_root}/direct/pp_inclusive_sim/"
                "run28_jet8_evil\n"
                f"1002 0 2 --output {exact_root}/writer/pp_inclusive_sim/"
                "run28_jet8_evil\n"
            ),
            (
                f"1001 0 2 --note={exact_root}/direct/pp_inclusive_sim/"
                "run28_jet8 unrelated-destination\n"
                f"1002 0 2 --note={exact_root}/writer/pp_inclusive_sim/"
                "run28_jet8 unrelated-destination\n"
            ),
            (
                f"1001 0 2 {exact_root}/writer/pp_inclusive_sim/"
                f"run28_jet8 {exact_root}/direct/pp_inclusive_sim/"
                "run28_jet8\n"
                f"1002 0 2 --output {exact_root}/writer/pp_inclusive_sim/"
                "run28_jet8\n"
            ),
        )
        for queue_text in queue_mutations:
            queue_path.write_text(queue_text, encoding="utf-8")
            queue_mutation = copy.deepcopy(valid)
            queue_mutation["initial_queue_tsv_sha256"] = validator.sha256_file(
                queue_path
            )
            terminal_path.write_text(
                json.dumps(queue_mutation, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
            with self.assertRaises(validator.ValidationError):
                validator.parse_terminal_gate_receipt(
                    terminal_path,
                    preflight=preflight,
                    output_root=Path(fields["base"]),
                )

        valid = self.write_terminal_gate(
            terminal_path,
            preflight=receipt,
            output_root=Path(fields["base"]),
            tag=fields["tag"],
        )
        mutations = []
        duplicate_identity = copy.deepcopy(valid)
        duplicate_identity["rows"][1]["cluster_id"] = 1001
        duplicate_identity["rows"][1]["cluster_proc"] = "1001.0"
        mutations.append(duplicate_identity)
        duplicate_role = copy.deepcopy(valid)
        duplicate_role["rows"][1]["role"] = "direct"
        mutations.append(duplicate_role)
        held_row = copy.deepcopy(valid)
        held_row["rows"][0]["job_status"] = 5
        mutations.append(held_row)
        wrong_exit = copy.deepcopy(valid)
        wrong_exit["rows"][0]["exit_code"] = 8
        mutations.append(wrong_exit)
        wrong_count = copy.deepcopy(valid)
        wrong_count["row_count"] = 1
        mutations.append(wrong_count)
        for mutation in mutations:
            terminal_path.write_text(
                json.dumps(mutation, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
            with self.assertRaises(validator.ValidationError):
                validator.parse_terminal_gate_receipt(
                    terminal_path,
                    preflight=preflight,
                    output_root=Path(fields["base"]),
                )

        self.write_terminal_gate(
            terminal_path,
            preflight=receipt,
            output_root=Path(fields["base"]),
            tag=fields["tag"],
        )
        initial_queue = self.root / "initial_queue.tsv"
        initial_queue.write_text(
            initial_queue.read_text(encoding="utf-8") + "1003 0 2 extra\n",
            encoding="utf-8",
        )
        with self.assertRaises(validator.ValidationError):
            validator.parse_terminal_gate_receipt(
                terminal_path,
                preflight=preflight,
                output_root=Path(fields["base"]),
            )

        mismatched_queue = self.root / "initial_queue.tsv"
        mismatched_queue.write_text(
            (
                f"9001 0 2 --output {Path(fields['base']).resolve()}/"
                "direct/pp_inclusive_sim/run28_jet8\n"
                f"9002 0 2 --output {Path(fields['base']).resolve()}/"
                "writer/pp_inclusive_sim/run28_jet8\n"
            ),
            encoding="utf-8",
        )
        mismatched_receipt = copy.deepcopy(valid)
        mismatched_receipt["initial_queue_tsv_sha256"] = validator.sha256_file(
            mismatched_queue
        )
        terminal_path.write_text(
            json.dumps(mismatched_receipt, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        with self.assertRaisesRegex(
            validator.ValidationError,
            "disagree with the bound initial queue",
        ):
            validator.parse_terminal_gate_receipt(
                terminal_path,
                preflight=preflight,
                output_root=Path(fields["base"]),
            )

    def test_pp_component_requires_fresh_full_artifact_revalidation(
        self,
    ) -> None:
        receipt = self.root / "pp_preflight_receipt.txt"
        self.write_preflight(
            receipt,
            validator.PP_ROWS,
            code_sha256="4" * 64,
            tag="the134_pp_revalidation_fixture",
        )
        payload = self.make_component_payload(
            system="pp",
            receipt=receipt,
            ownership_rows=validator.PP_ROWS,
            label="pp-revalidation",
        )
        certificate = self.root / "pp_component_certificate.json"
        certificate.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )

        # The historical false-PASS fixture is deliberately not a ROOT file.
        # Fresh certification must open it and fail instead of trusting the
        # certificate's stale readable/zombie/recovered booleans.
        with (
            patch.object(validator, "load_module", return_value=object()),
            self.assertRaises(validator.ValidationError),
        ):
            validator.load_pp_component_certificate(certificate)

        direct_base = (
            Path(payload["output_root"])
            / "direct"
            / validator.PP_ROWS[0]["lane"]
            / validator.PP_ROWS[0]["sample"]
        )
        self.write_histogram(
            direct_base / validator.SIDECAR_NAME,
            np.ones(3),
        )
        with (
            patch.object(validator, "load_module", return_value=object()),
            self.assertRaisesRegex(
                validator.ValidationError,
                "direct arm unexpectedly contains a sidecar",
            ),
        ):
            validator.load_pp_component_certificate(certificate)

        (direct_base / validator.SIDECAR_NAME).unlink()
        unexpected = (
            Path(payload["output_root"])
            / "direct"
            / "pp_photon_sim"
            / "unexpected"
            / "extra.root"
        )
        self.write_histogram(unexpected, np.ones(3))
        with (
            patch.object(validator, "load_module", return_value=object()),
            patch.object(
                validator,
                "validate_rows",
                return_value=copy.deepcopy(payload["rows"]),
            ),
            self.assertRaisesRegex(
                validator.ValidationError,
                "exact ROOT namespace differs",
            ),
        ):
            validator.load_pp_component_certificate(certificate)

        # Even with the expensive validator mocked, aggregate-time loading
        # must compare its regenerated payload to the stored certificate.
        regenerated = copy.deepcopy(payload)
        regenerated["rows"][0]["status"] = "FAIL"
        with (
            patch.object(
                validator,
                "certify_single_receipt",
                return_value=regenerated,
            ),
            self.assertRaisesRegex(
                validator.ValidationError,
                "fresh full revalidation",
            ),
        ):
            validator.load_pp_component_certificate(certificate)


if __name__ == "__main__":
    unittest.main()
