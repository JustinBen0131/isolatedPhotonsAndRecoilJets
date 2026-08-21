#!/usr/bin/env python3
"""Tests for exact THE-134 science-freeze manifest assembly."""

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

from the134_aggregate_contract import (  # noqa: E402
    PAIRED_REGISTRY_SCHEMA,
    REPLAY_CERTIFICATE_SCHEMA,
    REQUIRED_SYSTEM_VIEW_KEYS,
    REQUIRED_VIEW_KEYS,
    SCIENCE_MANIFEST_SCHEMA,
)


SCRIPT = VALIDATION / "assemble_the134_science_freeze_manifest.py"
PUBLIC_COMMIT = "1" * 40
CODE_SHA256 = "2" * 64


def write_json(path: Path, payload: dict) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class ManifestFixture:
    def __init__(self, root: Path):
        self.root = root
        self.registries: dict[str, Path] = {}
        self.replays: dict[tuple[str, str], Path] = {}
        for view in REQUIRED_VIEW_KEYS:
            path = root / f"registry_{view}.json"
            write_json(
                path,
                {
                    "schema": PAIRED_REGISTRY_SCHEMA,
                    "status": "READY",
                    "promotion_status": "CANDIDATE_CURRENT_NOT_CANONICAL",
                    "shower_definition": view,
                    "public_commit": PUBLIC_COMMIT,
                    "code_sha256": CODE_SHA256,
                },
            )
            self.registries[view] = path
        for system, view in REQUIRED_SYSTEM_VIEW_KEYS:
            path = root / f"replay_{system}_{view}.json"
            write_json(
                path,
                {
                    "schema": REPLAY_CERTIFICATE_SCHEMA,
                    "status": "PASS",
                    "promotion_status": "CANDIDATE_CURRENT_NOT_CANONICAL",
                    "system": system,
                    "shower_definition": view,
                    "public_commit": PUBLIC_COMMIT,
                    "code_sha256": CODE_SHA256,
                },
            )
            self.replays[(system, view)] = path

    def command(self, output: Path) -> list[str]:
        command = [
            sys.executable,
            str(SCRIPT),
            "--public-commit",
            PUBLIC_COMMIT,
            "--code-sha256",
            CODE_SHA256,
            "--json-out",
            str(output),
        ]
        for view in reversed(REQUIRED_VIEW_KEYS):
            command.extend(
                ["--registry", f"{view}={self.registries[view]}"]
            )
        for system, view in reversed(REQUIRED_SYSTEM_VIEW_KEYS):
            command.extend(
                [
                    "--replay-certificate",
                    f"{system}:{view}={self.replays[(system, view)]}",
                ]
            )
        return command


class AssembleThe134ScienceFreezeManifestTest(unittest.TestCase):
    def test_complete_manifest_is_deterministic_and_ordered(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fixture = ManifestFixture(root)
            output_a = root / "manifest_a.json"
            output_b = root / "manifest_b.json"
            result_a = subprocess.run(
                fixture.command(output_a), text=True, capture_output=True
            )
            result_b = subprocess.run(
                fixture.command(output_b), text=True, capture_output=True
            )
            self.assertEqual(result_a.returncode, 0, result_a.stderr)
            self.assertEqual(result_b.returncode, 0, result_b.stderr)
            self.assertEqual(output_a.read_bytes(), output_b.read_bytes())
            payload = json.loads(output_a.read_text())
            self.assertEqual(payload["schema"], SCIENCE_MANIFEST_SCHEMA)
            self.assertEqual(
                [item["view"] for item in payload["paired_registries"]],
                list(REQUIRED_VIEW_KEYS),
            )
            self.assertEqual(
                [
                    (item["system"], item["view"])
                    for item in payload["replay_certificates"]
                ],
                list(REQUIRED_SYSTEM_VIEW_KEYS),
            )
            for record in (
                payload["paired_registries"] + payload["replay_certificates"]
            ):
                self.assertEqual(
                    record["sha256"], sha256_file(Path(record["path"]))
                )

    def test_missing_or_duplicate_binding_fails_closed(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fixture = ManifestFixture(root)
            command = fixture.command(root / "missing.json")
            registry_index = command.index("--registry")
            del command[registry_index : registry_index + 2]
            result = subprocess.run(command, text=True, capture_output=True)
            self.assertEqual(result.returncode, 2)
            self.assertIn("inventory mismatch", result.stderr)

            command = fixture.command(root / "duplicate.json")
            command.extend(
                [
                    "--registry",
                    f"{REQUIRED_VIEW_KEYS[0]}={fixture.registries[REQUIRED_VIEW_KEYS[0]]}",
                ]
            )
            result = subprocess.run(command, text=True, capture_output=True)
            self.assertEqual(result.returncode, 2)
            self.assertIn("duplicate key", result.stderr)

    def test_identity_drift_fails_before_manifest_write(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fixture = ManifestFixture(root)
            system, view = REQUIRED_SYSTEM_VIEW_KEYS[0]
            path = fixture.replays[(system, view)]
            payload = json.loads(path.read_text())
            payload["shower_definition"] = "H0" if view != "H0" else "H70"
            write_json(path, payload)
            output = root / "identity_drift.json"
            result = subprocess.run(
                fixture.command(output), text=True, capture_output=True
            )
            self.assertEqual(result.returncode, 2)
            self.assertIn("top-level identity mismatch", result.stderr)
            self.assertFalse(output.exists())

    def test_malformed_hash_and_binding_syntax_fail_closed(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            fixture = ManifestFixture(root)
            command = fixture.command(root / "bad_hash.json")
            command[command.index(CODE_SHA256)] = "not-a-sha"
            result = subprocess.run(command, text=True, capture_output=True)
            self.assertEqual(result.returncode, 2)
            self.assertIn("code_sha256", result.stderr)

            command = fixture.command(root / "bad_binding.json")
            binding_index = command.index("--replay-certificate") + 1
            command[binding_index] = "missing-equals"
            result = subprocess.run(command, text=True, capture_output=True)
            self.assertEqual(result.returncode, 2)
            self.assertIn("KEY=PATH syntax", result.stderr)


if __name__ == "__main__":
    unittest.main()
