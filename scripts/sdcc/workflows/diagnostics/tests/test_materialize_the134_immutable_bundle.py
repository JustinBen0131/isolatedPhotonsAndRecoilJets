#!/usr/bin/env python3
"""Adversarial tests for THE-134 digest-bundle materialization."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import stat
import tempfile
import unittest
from pathlib import Path


HERE = Path(__file__).resolve().parent
MATERIALIZER_PATH = (
    HERE.parent / "materialize_the134_immutable_bundle.py"
)
SPEC = importlib.util.spec_from_file_location(
    "the134_bundle_materializer", MATERIALIZER_PATH
)
assert SPEC is not None and SPEC.loader is not None
materializer = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(materializer)
builder = materializer.builder


class BundleMaterializationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name).resolve()
        self.artifact_root = self.root / "source_artifacts"
        self.artifact_root.mkdir()
        self.release_lib = self.root / "release" / "lib"
        self.release_lib64 = self.root / "release" / "lib64"
        self.release_lib.mkdir(parents=True)
        self.release_lib64.mkdir(parents=True)
        self.bundle_parent = self.root / "bundles"
        self.bundle_parent.mkdir()
        self.declared_hashes = {
            "code_sha256": "1" * 64,
            "replay_schema_sha256": "2" * 64,
            "training_schema_sha256": "3" * 64,
            "semantic_sha256": "4" * 64,
        }
        self.artifact_paths: dict[str, Path] = {}
        for role in builder.REQUIRED_ARTIFACT_ROLES:
            path = self.artifact_root / f"{role}.artifact"
            if role == "code_manifest":
                path.write_text(
                    json.dumps(self.declared_hashes, sort_keys=True) + "\n",
                    encoding="utf-8",
                )
            else:
                path.write_text(
                    f"immutable artifact role={role}\n",
                    encoding="utf-8",
                )
            if role in builder.EXECUTABLE_ROLES:
                path.chmod(0o755)
            self.artifact_paths[role] = path
        self.inventory_receipt = self.make_inventory_receipt()
        self.inventory_path = self.root / "inventory_receipt.json"
        self.inventory_path.write_bytes(
            builder.canonical_json_bytes(self.inventory_receipt)
        )
        self.inventory_sha = self.sha256(self.inventory_path)

    def tearDown(self) -> None:
        for path in sorted(
            self.root.rglob("*"),
            key=lambda value: len(value.parts),
            reverse=True,
        ):
            try:
                if path.is_dir():
                    path.chmod(0o755)
                elif path.exists():
                    path.chmod(0o644)
            except FileNotFoundError:
                pass
        self.temporary.cleanup()

    @staticmethod
    def sha256(path: Path) -> str:
        return hashlib.sha256(path.read_bytes()).hexdigest()

    def make_inventory_receipt(self) -> dict:
        artifacts = []
        for role, path in sorted(self.artifact_paths.items()):
            artifacts.append(
                {
                    "role": role,
                    "path": str(path),
                    "sha256": self.sha256(path),
                    "size_bytes": path.stat().st_size,
                }
            )
        dependencies = []
        for role in sorted(builder.REQUIRED_DEPENDENCY_PROVIDER_ROLES):
            if role.startswith("pp_"):
                consumers = ["pp_executor"]
            elif role.startswith("auau_"):
                consumers = ["auau_executor"]
            else:
                consumers = ["auau_executor", "pp_executor"]
            dependencies.append(
                {
                    "provider_role": role,
                    "consumer_roles": consumers,
                    "kind": "runtime_or_contract",
                    "required": True,
                    "sha256": self.sha256(self.artifact_paths[role]),
                }
            )
        spec = {
            "schema": builder.SPEC_SCHEMA,
            "public_commit": "a" * 40,
            "bundle_parent": str(self.bundle_parent),
            "runtime": {
                "release": builder.RELEASE,
                "offline_main": "/cvmfs/example/ana.560",
                "calo_reco_soname": builder.CALO_RECO_SONAME,
                "request_memory_mb": builder.REQUEST_MEMORY_MB,
                "release_core_lib_dir": str(self.release_lib),
                "release_core_lib64_dir": str(self.release_lib64),
            },
            "declared_hashes": self.declared_hashes,
            "artifacts": artifacts,
            "dependencies": dependencies,
            "hash_bindings": [
                {
                    "name": name,
                    "artifact_role": "code_manifest",
                    "mode": "json_pointer",
                    "pointer": f"/{name}",
                }
                for name in builder.DECLARED_HASH_NAMES
            ],
        }
        return builder.build_receipt(spec)

    def materialize_once(self) -> dict:
        return materializer.materialize(
            self.inventory_path, self.inventory_sha
        )

    def test_materializes_readonly_symlink_free_bundle(self) -> None:
        summary = self.materialize_once()
        target = Path(summary["digest_named_bundle_path"])
        receipt_path = Path(summary["materialization_receipt"])
        self.assertTrue(target.is_dir())
        self.assertFalse(target.is_symlink())
        self.assertEqual(
            stat.S_IMODE(target.stat().st_mode),
            0o555,
        )
        readback = materializer.verify(
            receipt_path, summary["materialization_receipt_sha256"]
        )
        self.assertEqual(readback["status"], "PASS")
        self.assertEqual(
            readback["artifacts_rehashed"],
            len(builder.REQUIRED_ARTIFACT_ROLES),
        )
        resolver = builder.load_receipt(
            Path(summary["resolver_bundle_receipt"]),
            summary["resolver_bundle_receipt_sha256"],
        )
        for artifact in resolver["artifacts"]:
            path = Path(artifact["path"])
            self.assertTrue(path.is_relative_to(target))
            self.assertFalse(path.is_symlink())
            expected_mode = (
                0o555
                if artifact["role"] in builder.EXECUTABLE_ROLES
                else 0o444
            )
            self.assertEqual(
                stat.S_IMODE(path.stat().st_mode), expected_mode
            )

    def test_refuses_existing_digest_bundle(self) -> None:
        self.materialize_once()
        with self.assertRaisesRegex(
            materializer.MaterializationError,
            "refusing to overwrite",
        ):
            self.materialize_once()

    def test_readback_detects_copied_artifact_drift(self) -> None:
        summary = self.materialize_once()
        receipt_path = Path(summary["materialization_receipt"])
        payload = json.loads(receipt_path.read_text(encoding="utf-8"))
        artifact = Path(payload["artifact_inventory"][0]["bundle_path"])
        artifact.parent.chmod(0o755)
        artifact.chmod(0o644)
        artifact.write_bytes(artifact.read_bytes() + b"drift")
        with self.assertRaisesRegex(
            materializer.MaterializationError,
            "hash mismatch|content differs",
        ):
            materializer.verify(
                receipt_path,
                summary["materialization_receipt_sha256"],
            )

    def test_materialization_receipt_cannot_claim_science_authority(
        self,
    ) -> None:
        summary = self.materialize_once()
        path = Path(summary["materialization_receipt"])
        payload = json.loads(path.read_text(encoding="utf-8"))
        payload["authority_state"] = "SCIENTIFIC_AUTHORITY"
        forged = self.root / "forged.json"
        forged.write_bytes(materializer.canonical_json_bytes(payload))
        with self.assertRaisesRegex(
            materializer.MaterializationError,
            "schema/status/authority",
        ):
            materializer.verify(forged, self.sha256(forged))


if __name__ == "__main__":
    unittest.main()
