#!/usr/bin/env python3
"""Tests for the source-complete THE-134 replay-certificate orchestrator."""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
VALIDATION = ROOT / "scripts" / "ml" / "validation"
sys.path.insert(0, str(VALIDATION))

from build_the134_source_complete_replay_certificates import (  # noqa: E402
    INDEX_NAME,
    INPUT_MANIFEST_SCHEMA,
    OUTPUT_INDEX_SCHEMA,
    PROMOTION_STATUS,
    REQUIRED_SOURCE_VIEW_KEYS,
    SOURCE_VALIDATION_RECEIPT_SCHEMA,
    certificate_name,
)
from build_the134_science_freeze_certificate import (  # noqa: E402
    validate_source_witnesses,
)
from the134_aggregate_contract import (  # noqa: E402
    PAIRED_REGISTRY_SCHEMA,
    REPLAY_CERTIFICATE_GATES,
    REPLAY_CERTIFICATE_SCHEMA,
    REQUIRED_SOURCES_BY_SYSTEM,
    REQUIRED_SYSTEM_VIEW_KEYS,
    REQUIRED_VIEW_KEYS,
    SOURCE_WITNESS_GATES,
    expected_view_semantic,
)


SCRIPT = VALIDATION / "build_the134_source_complete_replay_certificates.py"
PUBLIC_COMMIT = "1" * 40
CODE_SHA256 = "2" * 64


def write_json(path: Path, payload: dict) -> None:
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class SourceCompleteFixture:
    def __init__(self, root: Path):
        self.root = root
        self.registries: list[dict] = []
        self.source_records: list[dict] = []
        self.stable_artifacts: dict[tuple[str, str], dict[str, Path]] = {}

        for view in REQUIRED_VIEW_KEYS:
            path = root / f"registry_{view}.json"
            write_json(
                path,
                {
                    "schema": PAIRED_REGISTRY_SCHEMA,
                    "status": "READY",
                    "promotion_status": PROMOTION_STATUS,
                    "shower_definition": view,
                    "shower_semantic_sha256": expected_view_semantic(view),
                    "public_commit": PUBLIC_COMMIT,
                    "code_sha256": CODE_SHA256,
                },
            )
            self.registries.append(
                {"view": view, "path": str(path), "sha256": sha256_file(path)}
            )

        for system, sources in REQUIRED_SOURCES_BY_SYSTEM.items():
            for source in sources:
                paths = {}
                for artifact_name in ("source_manifest", "direct", "writer"):
                    path = root / f"{system}_{source}_{artifact_name}.bin"
                    path.write_bytes(
                        f"{system}|{source}|{artifact_name}\n".encode("utf-8")
                    )
                    paths[artifact_name] = path
                self.stable_artifacts[(system, source)] = paths

        for system, view, source in REQUIRED_SOURCE_VIEW_KEYS:
            stable = self.stable_artifacts[(system, source)]
            cache = root / f"{system}_{view}_{source}_cache_receipt.json"
            write_json(
                cache,
                {
                    "system": system,
                    "view": view,
                    "source": source,
                    "status": "PASS",
                },
            )
            artifact_records = {
                name: {"path": str(path), "sha256": sha256_file(path)}
                for name, path in stable.items()
            }
            artifact_records["cache_receipt"] = {
                "path": str(cache),
                "sha256": sha256_file(cache),
            }
            receipt = root / f"{system}_{view}_{source}_validation.json"
            write_json(
                receipt,
                {
                    "schema": SOURCE_VALIDATION_RECEIPT_SCHEMA,
                    "status": "PASS",
                    "promotion_status": PROMOTION_STATUS,
                    "system": system,
                    "shower_definition": view,
                    "shower_semantic_sha256": expected_view_semantic(view),
                    "source": source,
                    "public_commit": PUBLIC_COMMIT,
                    "code_sha256": CODE_SHA256,
                    "source_manifest_sha256": artifact_records[
                        "source_manifest"
                    ]["sha256"],
                    "direct_artifact_sha256": artifact_records["direct"][
                        "sha256"
                    ],
                    "writer_artifact_sha256": artifact_records["writer"][
                        "sha256"
                    ],
                    "cache_receipt_sha256": artifact_records["cache_receipt"][
                        "sha256"
                    ],
                    "gates": {gate: True for gate in SOURCE_WITNESS_GATES},
                },
            )
            self.source_records.append(
                {
                    "system": system,
                    "view": view,
                    "source": source,
                    "validation_receipt": {
                        "path": str(receipt),
                        "sha256": sha256_file(receipt),
                    },
                    "artifacts": artifact_records,
                }
            )

        self.manifest = root / "input_manifest.json"
        self.write_manifest()

    def write_manifest(self) -> None:
        write_json(
            self.manifest,
            {
                "schema": INPUT_MANIFEST_SCHEMA,
                "public_commit": PUBLIC_COMMIT,
                "code_sha256": CODE_SHA256,
                "paired_registries": self.registries,
                "source_validations": self.source_records,
            },
        )

    def command(self, output: Path) -> list[str]:
        return [
            sys.executable,
            str(SCRIPT),
            "--input-manifest",
            str(self.manifest),
            "--output-directory",
            str(output),
        ]

    def first_source_record(self) -> dict:
        return self.source_records[0]

    def update_receipt_hash(self, record: dict) -> None:
        receipt = Path(record["validation_receipt"]["path"])
        record["validation_receipt"]["sha256"] = sha256_file(receipt)
        self.write_manifest()


class BuildSourceCompleteReplayCertificatesTest(unittest.TestCase):
    def test_builds_exact_deterministic_fourteen_certificate_bundle(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fixture = SourceCompleteFixture(root)
            output_a = root / "certificates_a"
            output_b = root / "certificates_b"

            result_a = subprocess.run(
                fixture.command(output_a), text=True, capture_output=True
            )
            result_b = subprocess.run(
                fixture.command(output_b), text=True, capture_output=True
            )
            self.assertEqual(result_a.returncode, 0, result_a.stderr)
            self.assertEqual(result_b.returncode, 0, result_b.stderr)

            expected_names = {
                certificate_name(system, view)
                for system, view in REQUIRED_SYSTEM_VIEW_KEYS
            }
            self.assertEqual(
                {path.name for path in output_a.iterdir()},
                expected_names | {INDEX_NAME},
            )
            index = json.loads((output_a / INDEX_NAME).read_text())
            self.assertEqual(index["schema"], OUTPUT_INDEX_SCHEMA)
            self.assertEqual(index["status"], "PASS")
            self.assertEqual(index["certificate_count"], 14)
            self.assertEqual(
                index["input_manifest_sha256"], sha256_file(fixture.manifest)
            )
            self.assertEqual(
                [
                    (record["system"], record["view"])
                    for record in index["certificates"]
                ],
                list(REQUIRED_SYSTEM_VIEW_KEYS),
            )

            for record in index["certificates"]:
                path = Path(record["path"])
                self.assertEqual(record["sha256"], sha256_file(path))
                certificate = json.loads(path.read_text())
                self.assertEqual(
                    certificate["schema"], REPLAY_CERTIFICATE_SCHEMA
                )
                self.assertEqual(certificate["status"], "PASS")
                self.assertEqual(
                    set(certificate["gates"]), set(REPLAY_CERTIFICATE_GATES)
                )
                self.assertTrue(all(certificate["gates"].values()))
                self.assertEqual(
                    len(certificate["source_witnesses"]),
                    len(REQUIRED_SOURCES_BY_SYSTEM[record["system"]]),
                )
                self.assertTrue(
                    all(
                        witness["system"] == record["system"]
                        and witness["shower_definition"] == record["view"]
                        for witness in certificate["source_witnesses"]
                    )
                )
                self.assertEqual(
                    len(
                        validate_source_witnesses(
                            record["system"], record["view"], certificate
                        )
                    ),
                    len(REQUIRED_SOURCES_BY_SYSTEM[record["system"]]),
                )
                self.assertEqual(
                    path.read_bytes(),
                    (
                        output_b
                        / certificate_name(record["system"], record["view"])
                    ).read_bytes(),
                )

            idempotent = subprocess.run(
                fixture.command(output_a), text=True, capture_output=True
            )
            self.assertEqual(idempotent.returncode, 0, idempotent.stderr)

    def test_missing_and_duplicate_source_inventory_fail_before_output(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fixture = SourceCompleteFixture(root)
            missing_record = fixture.source_records.pop()
            fixture.write_manifest()
            missing_output = root / "missing_output"
            result = subprocess.run(
                fixture.command(missing_output), text=True, capture_output=True
            )
            self.assertEqual(result.returncode, 2)
            self.assertIn("inventory mismatch", result.stderr)
            self.assertFalse(missing_output.exists())

            fixture.source_records.append(missing_record)
            fixture.source_records.append(dict(fixture.source_records[0]))
            fixture.write_manifest()
            duplicate_output = root / "duplicate_output"
            result = subprocess.run(
                fixture.command(duplicate_output), text=True, capture_output=True
            )
            self.assertEqual(result.returncode, 2)
            self.assertIn("duplicate keys", result.stderr)
            self.assertFalse(duplicate_output.exists())

    def test_incomplete_or_false_receipt_gate_inventory_fails_closed(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fixture = SourceCompleteFixture(root)
            record = fixture.first_source_record()
            receipt_path = Path(record["validation_receipt"]["path"])
            receipt = json.loads(receipt_path.read_text())
            removed_gate = SOURCE_WITNESS_GATES[0]
            receipt["gates"].pop(removed_gate)
            write_json(receipt_path, receipt)
            fixture.update_receipt_hash(record)

            missing_gate_output = root / "missing_gate_output"
            result = subprocess.run(
                fixture.command(missing_gate_output),
                text=True,
                capture_output=True,
            )
            self.assertEqual(result.returncode, 2)
            self.assertIn("gate inventory mismatch", result.stderr)
            self.assertFalse(missing_gate_output.exists())

            receipt["gates"][removed_gate] = False
            write_json(receipt_path, receipt)
            fixture.update_receipt_hash(record)
            false_gate_output = root / "false_gate_output"
            result = subprocess.run(
                fixture.command(false_gate_output),
                text=True,
                capture_output=True,
            )
            self.assertEqual(result.returncode, 2)
            self.assertIn("gates are not all true", result.stderr)
            self.assertFalse(false_gate_output.exists())

    def test_artifact_drift_and_receipt_binding_mismatch_fail_closed(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fixture = SourceCompleteFixture(root)
            record = fixture.first_source_record()
            direct_path = Path(record["artifacts"]["direct"]["path"])
            direct_path.write_bytes(direct_path.read_bytes() + b"drift\n")

            drift_output = root / "drift_output"
            result = subprocess.run(
                fixture.command(drift_output), text=True, capture_output=True
            )
            self.assertEqual(result.returncode, 2)
            self.assertIn("SHA-256 mismatch", result.stderr)
            self.assertFalse(drift_output.exists())

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fixture = SourceCompleteFixture(root)
            record = fixture.first_source_record()
            receipt_path = Path(record["validation_receipt"]["path"])
            receipt = json.loads(receipt_path.read_text())
            receipt["direct_artifact_sha256"] = "f" * 64
            write_json(receipt_path, receipt)
            fixture.update_receipt_hash(record)

            mismatch_output = root / "binding_mismatch_output"
            result = subprocess.run(
                fixture.command(mismatch_output),
                text=True,
                capture_output=True,
            )
            self.assertEqual(result.returncode, 2)
            self.assertIn(
                "registration disagrees with validation receipt", result.stderr
            )
            self.assertFalse(mismatch_output.exists())


if __name__ == "__main__":
    unittest.main()
