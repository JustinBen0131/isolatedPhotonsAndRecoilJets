#!/usr/bin/env python3
"""Adversarial tests for the immutable THE-134 bundle receipt."""

from __future__ import annotations

import copy
import hashlib
import importlib.util
import json
import os
import tempfile
import unittest
from pathlib import Path


HERE = Path(__file__).resolve().parent
MODULE_PATH = HERE.parent / "build_the134_immutable_bundle_receipt.py"
SPEC = importlib.util.spec_from_file_location("the134_bundle_builder", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
builder = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(builder)


class ImmutableBundleReceiptTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name).resolve()
        self.artifact_root = self.root / "artifacts"
        self.artifact_root.mkdir()
        self.lib = self.root / "release" / "lib"
        self.lib64 = self.root / "release" / "lib64"
        self.lib.mkdir(parents=True)
        self.lib64.mkdir(parents=True)

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
                path.write_text(f"immutable artifact role={role}\n", encoding="utf-8")
            if role in builder.EXECUTABLE_ROLES:
                os.chmod(path, 0o755)
            self.artifact_paths[role] = path
        self.spec = self.make_spec()

    def tearDown(self) -> None:
        self.temporary.cleanup()

    @staticmethod
    def file_sha256(path: Path) -> str:
        return hashlib.sha256(path.read_bytes()).hexdigest()

    def make_spec(self) -> dict[str, object]:
        artifacts = [
            {
                "role": role,
                "path": str(path),
                "sha256": self.file_sha256(path),
                "size_bytes": path.stat().st_size,
            }
            for role, path in sorted(self.artifact_paths.items())
        ]
        dependencies = []
        for provider in sorted(builder.REQUIRED_DEPENDENCY_PROVIDER_ROLES):
            if provider.startswith("pp_"):
                consumers = ["pp_executor"]
            elif provider.startswith("auau_"):
                consumers = ["auau_executor"]
            else:
                consumers = ["pp_executor", "auau_executor"]
            dependencies.append(
                {
                    "provider_role": provider,
                    "consumer_roles": consumers,
                    "kind": "runtime_or_contract",
                    "required": True,
                    "sha256": self.file_sha256(self.artifact_paths[provider]),
                }
            )
        bindings = [
            {
                "name": name,
                "artifact_role": "code_manifest",
                "mode": "json_pointer",
                "pointer": f"/{name}",
            }
            for name in builder.DECLARED_HASH_NAMES
        ]
        return {
            "schema": builder.SPEC_SCHEMA,
            "public_commit": "a" * 40,
            "bundle_parent": "/sphenix/u/test/immutable-bundles",
            "runtime": {
                "release": "ana.560",
                "offline_main": "/cvmfs/sphenix.example/release/ana.560",
                "calo_reco_soname": "libcalo_reco.so.0",
                "request_memory_mb": 3000,
                "release_core_lib_dir": str(self.lib),
                "release_core_lib64_dir": str(self.lib64),
            },
            "declared_hashes": dict(self.declared_hashes),
            "artifacts": artifacts,
            "dependencies": dependencies,
            "hash_bindings": bindings,
        }

    def test_receipt_is_deterministic_and_digest_named(self) -> None:
        first = builder.build_receipt(self.spec)
        second = builder.build_receipt(self.spec)
        self.assertEqual(
            builder.canonical_json_bytes(first),
            builder.canonical_json_bytes(second),
        )
        identity = first["bundle_identity_sha256"]
        self.assertEqual(first["authority_state"], builder.AUTHORITY_STATE)
        self.assertEqual(
            first["bundle_name"], f"{builder.BUNDLE_NAME_PREFIX}{identity}"
        )
        self.assertTrue(first["digest_named_bundle_path"].endswith(identity))
        self.assertFalse(first["digest_named_bundle_materialized"])
        self.assertEqual(len(first["artifacts"]), 22)
        self.assertEqual(
            len(first["dependencies"]),
            len(builder.REQUIRED_DEPENDENCY_PROVIDER_ROLES),
        )

    def test_pinned_pre_submit_readback_passes(self) -> None:
        receipt = builder.build_receipt(self.spec)
        path = self.root / "bundle_receipt.json"
        builder.atomic_write_json(path, receipt)
        expected = self.file_sha256(path)
        observed = builder.load_receipt(path, expected)
        self.assertEqual(observed, receipt)

    def test_mutable_artifact_path_drift_is_rejected(self) -> None:
        receipt = builder.build_receipt(self.spec)
        path = self.root / "bundle_receipt.json"
        builder.atomic_write_json(path, receipt)
        expected = self.file_sha256(path)
        self.artifact_paths["pp_library"].write_text(
            "mutated after receipt\n", encoding="utf-8"
        )
        with self.assertRaisesRegex(builder.BundleError, "hash mismatch"):
            builder.load_receipt(path, expected)

    def test_declared_artifact_hash_mismatch_is_rejected(self) -> None:
        mutated = copy.deepcopy(self.spec)
        for artifact in mutated["artifacts"]:
            if artifact["role"] == "auau_library":
                artifact["sha256"] = "f" * 64
                break
        with self.assertRaisesRegex(builder.BundleError, "hash mismatch"):
            builder.build_receipt(mutated)

    def test_incomplete_dependency_inventory_is_rejected(self) -> None:
        mutated = copy.deepcopy(self.spec)
        mutated["dependencies"].pop()
        with self.assertRaisesRegex(builder.BundleError, "closure differs"):
            builder.build_receipt(mutated)

    def test_declared_semantic_hash_binding_mismatch_is_rejected(self) -> None:
        mutated = copy.deepcopy(self.spec)
        mutated["declared_hashes"]["semantic_sha256"] = "9" * 64
        with self.assertRaisesRegex(builder.BundleError, "binding mismatch"):
            builder.build_receipt(mutated)

    def test_artifact_role_alias_is_rejected(self) -> None:
        mutated = copy.deepcopy(self.spec)
        pp_model = next(
            artifact
            for artifact in mutated["artifacts"]
            if artifact["role"] == "pp_model"
        )
        auau_model = next(
            artifact
            for artifact in mutated["artifacts"]
            if artifact["role"] == "auau_model"
        )
        auau_model["path"] = pp_model["path"]
        auau_model["sha256"] = pp_model["sha256"]
        auau_model["size_bytes"] = pp_model["size_bytes"]
        with self.assertRaisesRegex(builder.BundleError, "alias"):
            builder.build_receipt(mutated)


if __name__ == "__main__":
    unittest.main()
