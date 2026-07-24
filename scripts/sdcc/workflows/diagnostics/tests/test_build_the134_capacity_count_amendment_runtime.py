#!/usr/bin/env python3
"""Focused runtime-authority tests for the THE-134 count amendment."""

from __future__ import annotations

import copy
import hashlib
import importlib.util
import json
import os
import shutil
import tempfile
import unittest
from pathlib import Path
from typing import Any


HERE = Path(__file__).resolve().parent
MODULE_PATH = HERE.parent / "build_the134_capacity_count_amendment.py"
SPEC = importlib.util.spec_from_file_location(
    "the134_capacity_count_amendment_runtime", MODULE_PATH
)
assert SPEC is not None and SPEC.loader is not None
amendment = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(amendment)


class CapacityCountAmendmentRuntimeTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name).resolve()

    def tearDown(self) -> None:
        self.temporary.cleanup()

    @staticmethod
    def sha256(path: Path) -> str:
        return hashlib.sha256(path.read_bytes()).hexdigest()

    def artifact(self, path: Path) -> dict[str, str]:
        return {"path": str(path), "sha256": self.sha256(path)}

    def bundle_artifact(self, role: str, path: Path) -> dict[str, Any]:
        return {
            "role": role,
            "path": str(path),
            "resolved_path": str(path.resolve(strict=True)),
            "sha256": self.sha256(path),
            "size_bytes": path.stat().st_size,
        }

    @staticmethod
    def config_text(
        values: dict[str, str],
        transition_only: dict[str, str],
        *,
        residual: str = "stable_contract: frozen",
    ) -> str:
        lines = [f"{key}: {value}" for key, value in values.items()]
        lines.extend(
            f"{key}: {value}" for key, value in transition_only.items()
        )
        lines.append(residual)
        return "\n".join(lines) + "\n"

    def make_config_chain(
        self, system: str
    ) -> dict[str, Any]:
        chain_root = self.root / f"{system}_config_chain"
        chain_root.mkdir()
        expected_transition = amendment.EXPECTED_BASE_TO_RESOLVED_VALUES[
            system
        ]
        expected_materialization = (
            amendment.EXPECTED_RESOLVED_MATERIALIZATION_VALUES[system]
        )

        base_materialization = dict(
            expected_materialization["resolved"]
        )
        base_materialization.update(
            {
                key: value
                for key, value in expected_transition["base"].items()
                if key in amendment.MATERIALIZED_CONFIG_KEYS
            }
        )
        resolved_materialization = dict(
            expected_materialization["resolved"]
        )
        materialized_values = dict(
            expected_materialization["materialized"]
        )

        transition_only = {
            key: value
            for key, value in expected_transition["resolved"].items()
            if key not in amendment.MATERIALIZED_CONFIG_KEYS
        }
        base_transition_only = {
            key: value
            for key, value in expected_transition["base"].items()
            if key not in amendment.MATERIALIZED_CONFIG_KEYS
        }

        base_path = chain_root / "base.yaml"
        resolved_path = chain_root / "resolved.yaml"
        materialized_path = chain_root / "materialized.yaml"
        base_path.write_text(
            self.config_text(
                base_materialization,
                base_transition_only,
            ),
            encoding="utf-8",
        )
        resolved_path.write_text(
            self.config_text(
                resolved_materialization,
                transition_only,
            ),
            encoding="utf-8",
        )
        materialized_path.write_text(
            self.config_text(
                materialized_values,
                transition_only,
            ),
            encoding="utf-8",
        )

        output_namespace = chain_root / "outputs"
        output_namespace.mkdir()
        output_name = f"{system}_analysis.root"
        tags = expected_materialization["materialized"]
        fanout_path = chain_root / "fanout.txt"
        fanout_path.write_text(
            "|".join(
                (
                    str(output_namespace / output_name),
                    output_name,
                    tags["preselection"],
                    tags["tight"],
                    tags["nonTight"],
                )
            )
            + "\n",
            encoding="utf-8",
        )
        return {
            "row_id": f"{system}_background_jet8",
            "system": system,
            "base_path": base_path,
            "resolved_path": resolved_path,
            "materialized_path": materialized_path,
            "fanout_path": fanout_path,
            "output_namespace": output_namespace,
        }

    def validate_config_chain(self, fixture: dict[str, Any]) -> str:
        return amendment.validate_runtime_config_chain(
            fixture["row_id"],
            fixture["system"],
            base_config=self.artifact(fixture["base_path"]),
            resolved_config=self.artifact(fixture["resolved_path"]),
            materialized_config=self.artifact(
                fixture["materialized_path"]
            ),
            fanout_contract=self.artifact(fixture["fanout_path"]),
            analysis_output_namespace=str(fixture["output_namespace"]),
        )

    def make_runtime_authority(
        self, name: str
    ) -> tuple[dict[str, str], dict[str, Any], dict[str, Any]]:
        fixture_root = self.root / name
        bundle_root = fixture_root / "bundle_origins"
        authority_root = fixture_root / "authority_copies"
        offline_main = fixture_root / "release" / "offline"
        release_lib = offline_main / "lib"
        release_lib64 = offline_main / "lib64"
        for directory in (
            bundle_root,
            authority_root,
            release_lib,
            release_lib64,
        ):
            directory.mkdir(parents=True, exist_ok=True)

        role_contents = {
            "calo_reco_build_receipt": b"calo-build-receipt\n",
            "calo_reco_library": b"calo-reco-library\n",
            "calo_reco_source_manifest": b"calo-source-manifest\n",
            "release_calo_io": b"release-calo-io\n",
            "release_clusteriso": b"release-clusteriso\n",
            "release_jetbase": b"release-jetbase\n",
        }
        bundle_artifacts: dict[str, dict[str, Any]] = {}
        bundle_paths: dict[str, Path] = {}
        for role, content in role_contents.items():
            path = bundle_root / role
            path.write_bytes(content)
            bundle_paths[role] = path
            bundle_artifacts[role] = self.bundle_artifact(role, path)

        calo_authority_paths = {
            "build_receipt": authority_root / "calo_build_receipt.json",
            "library": authority_root / "libcalo_reco.so.1",
            "source_manifest": authority_root / "calo_source_manifest.json",
        }
        for field, role in (
            ("build_receipt", "calo_reco_build_receipt"),
            ("library", "calo_reco_library"),
            ("source_manifest", "calo_reco_source_manifest"),
        ):
            shutil.copyfile(bundle_paths[role], calo_authority_paths[field])

        provider_paths = {
            "libcalo_io.so": release_lib / "libcalo_io.so",
            "libclusteriso.so": release_lib / "libclusteriso.so",
            "libjetbase.so": release_lib64 / "libjetbase.so",
        }
        for family, role in amendment.RUNTIME_PROVIDER_ROLE_BY_FAMILY.items():
            shutil.copyfile(bundle_paths[role], provider_paths[family])

        payload = {
            "schema": amendment.RUNTIME_AUTHORITY_SCHEMA,
            "status": "PASS",
            "calo_reco": {
                "build_receipt": str(
                    calo_authority_paths["build_receipt"]
                ),
                "build_receipt_sha256": self.sha256(
                    calo_authority_paths["build_receipt"]
                ),
                "library": str(calo_authority_paths["library"]),
                "library_sha256": self.sha256(
                    calo_authority_paths["library"]
                ),
                "soname": "libcalo_reco.so.1",
                "source_manifest": str(
                    calo_authority_paths["source_manifest"]
                ),
                "source_manifest_sha256": self.sha256(
                    calo_authority_paths["source_manifest"]
                ),
            },
            "release": {
                "name": "ana.541",
                "offline_main": str(offline_main),
                "lib": str(release_lib),
                "lib64": str(release_lib64),
                "providers": {
                    family: {
                        "path": str(path),
                        "sha256": self.sha256(path),
                    }
                    for family, path in provider_paths.items()
                },
            },
            "pp_sim_weight_contract": copy.deepcopy(
                amendment.EXPECTED_PP_SIM_WEIGHT_CONTRACT
            ),
        }
        authority_path = authority_root / "runtime_authority.json"
        authority_path.write_text(
            json.dumps(payload, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        immutable_authority_path = (
            bundle_root / "runtime_authority_manifest"
        )
        shutil.copyfile(authority_path, immutable_authority_path)
        bundle_artifacts["runtime_authority_manifest"] = (
            self.bundle_artifact(
                "runtime_authority_manifest",
                immutable_authority_path,
            )
        )
        immutable = {
            "runtime": {
                "release": "ana.541",
                "offline_main": str(offline_main),
                "release_core_lib_dir": str(release_lib),
                "release_core_lib64_dir": str(release_lib64),
                "calo_reco_soname": "libcalo_reco.so.1",
            },
            "bundle_artifacts": bundle_artifacts,
        }
        return self.artifact(authority_path), immutable, payload

    def rewrite_runtime_authority(
        self,
        artifact: dict[str, str],
        immutable: dict[str, Any],
        payload: dict[str, Any],
    ) -> None:
        path = Path(artifact["path"])
        path.write_text(
            json.dumps(payload, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        digest = self.sha256(path)
        artifact["sha256"] = digest
        runtime_role = immutable["bundle_artifacts"][
            "runtime_authority_manifest"
        ]
        runtime_role["sha256"] = digest
        runtime_role["size_bytes"] = path.stat().st_size

    def test_runtime_config_chain_accepts_exact_pp_and_auau(self) -> None:
        for system in ("pp", "auau"):
            with self.subTest(system=system):
                fixture = self.make_config_chain(system)
                observed = self.validate_config_chain(fixture)
                self.assertRegex(observed, r"^[0-9a-f]{64}$")

    def test_runtime_config_chain_rejects_unauthorized_residual_drift(
        self,
    ) -> None:
        for system in ("pp", "auau"):
            with self.subTest(system=system):
                fixture = self.make_config_chain(system)
                resolved = fixture["resolved_path"]
                resolved.write_text(
                    resolved.read_text(encoding="utf-8").replace(
                        "stable_contract: frozen",
                        "stable_contract: drifted",
                    ),
                    encoding="utf-8",
                )
                with self.assertRaisesRegex(
                    amendment.AmendmentError,
                    "scalarization differs",
                ):
                    self.validate_config_chain(fixture)

    def test_runtime_config_chain_rejects_materialized_value_drift(
        self,
    ) -> None:
        mutations = {
            "pp": ("fixedGeV: 2.0", "fixedGeV: 3.0"),
            "auau": ("vz_cut_cm: 10", "vz_cut_cm: 11"),
        }
        for system, (old, new) in mutations.items():
            with self.subTest(system=system):
                fixture = self.make_config_chain(system)
                materialized = fixture["materialized_path"]
                materialized.write_text(
                    materialized.read_text(encoding="utf-8").replace(
                        old, new
                    ),
                    encoding="utf-8",
                )
                with self.assertRaisesRegex(
                    amendment.AmendmentError,
                    "scalarization differs",
                ):
                    self.validate_config_chain(fixture)

    def test_runtime_config_chain_rejects_fanout_mismatch(self) -> None:
        for system in ("pp", "auau"):
            with self.subTest(system=system):
                fixture = self.make_config_chain(system)
                fanout = fixture["fanout_path"]
                fanout.write_text(
                    fanout.read_text(encoding="utf-8").replace(
                        "newPPG12", "wrongSelection", 1
                    ),
                    encoding="utf-8",
                )
                with self.assertRaisesRegex(
                    amendment.AmendmentError,
                    "fanout selection identity differs",
                ):
                    self.validate_config_chain(fixture)

    def test_runtime_authority_accepts_legitimate_byte_copies(self) -> None:
        artifact, immutable, payload = self.make_runtime_authority("valid")
        observed = amendment.validate_runtime_authority_artifact(
            artifact, immutable
        )
        self.assertEqual(observed["payload"], payload)
        for field, role in (
            ("build_receipt", "calo_reco_build_receipt"),
            ("library", "calo_reco_library"),
            ("source_manifest", "calo_reco_source_manifest"),
        ):
            self.assertFalse(
                os.path.samefile(
                    observed["calo_reco"][field]["path"],
                    immutable["bundle_artifacts"][role]["path"],
                )
            )
            self.assertEqual(
                observed["calo_reco"][field]["sha256"],
                immutable["bundle_artifacts"][role]["sha256"],
            )
        for family, role in amendment.RUNTIME_PROVIDER_ROLE_BY_FAMILY.items():
            self.assertFalse(
                os.path.samefile(
                    observed["providers"][family]["path"],
                    immutable["bundle_artifacts"][role]["path"],
                )
            )
            self.assertEqual(
                observed["providers"][family]["sha256"],
                immutable["bundle_artifacts"][role]["sha256"],
            )

    def test_runtime_authority_rejects_contract_mutations(self) -> None:
        def mutate_schema(
            payload: dict[str, Any],
            _artifact: dict[str, str],
            _immutable: dict[str, Any],
        ) -> None:
            payload["schema"] = "WRONG_SCHEMA"

        def mutate_key(
            payload: dict[str, Any],
            _artifact: dict[str, str],
            _immutable: dict[str, Any],
        ) -> None:
            payload["unexpected"] = True

        def mutate_hash(
            payload: dict[str, Any],
            _artifact: dict[str, str],
            _immutable: dict[str, Any],
        ) -> None:
            payload["calo_reco"]["library_sha256"] = "0" * 64

        def mutate_provider(
            payload: dict[str, Any],
            artifact: dict[str, str],
            _immutable: dict[str, Any],
        ) -> None:
            source = Path(
                payload["release"]["providers"]["libcalo_io.so"]["path"]
            )
            outsider = Path(artifact["path"]).parent / "outside_calo_io.so"
            shutil.copyfile(source, outsider)
            payload["release"]["providers"]["libcalo_io.so"]["path"] = str(
                outsider
            )

        def mutate_release(
            payload: dict[str, Any],
            _artifact: dict[str, str],
            _immutable: dict[str, Any],
        ) -> None:
            payload["release"]["name"] = "ana.999"

        def mutate_weight(
            payload: dict[str, Any],
            _artifact: dict[str, str],
            _immutable: dict[str, Any],
        ) -> None:
            payload["pp_sim_weight_contract"]["period_lumi_weight"] = False

        mutations = {
            "schema": mutate_schema,
            "key": mutate_key,
            "hash": mutate_hash,
            "provider": mutate_provider,
            "release": mutate_release,
            "weight": mutate_weight,
        }
        for label, mutator in mutations.items():
            with self.subTest(label=label):
                artifact, immutable, payload = self.make_runtime_authority(
                    f"mutation_{label}"
                )
                mutator(payload, artifact, immutable)
                self.rewrite_runtime_authority(
                    artifact, immutable, payload
                )
                with self.assertRaises(amendment.AmendmentError):
                    amendment.validate_runtime_authority_artifact(
                        artifact, immutable
                    )

    def test_read_tsv_accepts_strict_rectangular_rows(self) -> None:
        path = self.root / "valid.tsv"
        path.write_text("first\tsecond\none\ttwo\n", encoding="utf-8")
        header, rows = amendment.read_tsv(path)
        self.assertEqual(header, ["first", "second"])
        self.assertEqual(rows, [{"first": "one", "second": "two"}])

    def test_read_tsv_rejects_duplicate_overflow_and_empty_cells(
        self,
    ) -> None:
        cases = {
            "duplicate": "first\tfirst\none\ttwo\n",
            "overflow": "first\tsecond\none\ttwo\tthree\n",
            "empty": "first\tsecond\none\t\n",
        }
        for label, text in cases.items():
            with self.subTest(label=label):
                path = self.root / f"{label}.tsv"
                path.write_text(text, encoding="utf-8")
                with self.assertRaises(amendment.AmendmentError):
                    amendment.read_tsv(path)

    def test_strict_json_rejects_duplicate_keys_at_every_depth(self) -> None:
        for text in (
            '{"duplicate":1,"duplicate":2}',
            '{"outer":{"duplicate":1,"duplicate":2}}',
        ):
            with self.subTest(text=text):
                with self.assertRaisesRegex(
                    amendment.AmendmentError,
                    "duplicates key: duplicate",
                ):
                    amendment.strict_json_loads(text)

        path = self.root / "duplicate.json"
        path.write_text('{"schema":"A","schema":"B"}\n', encoding="utf-8")
        with self.assertRaisesRegex(
            amendment.AmendmentError,
            "duplicates key: schema",
        ):
            amendment.load_json_artifact(
                "duplicate JSON artifact",
                self.artifact(path),
            )

    def test_json_contract_rejects_nonfinite_numbers(self) -> None:
        for constant in ("NaN", "Infinity", "-Infinity"):
            with self.subTest(constant=constant):
                with self.assertRaisesRegex(
                    amendment.AmendmentError, "non-finite constant"
                ):
                    amendment.strict_json_loads(
                        f'{{"value":{constant}}}'
                    )
        for value in (float("nan"), float("inf"), float("-inf")):
            with self.subTest(value=value):
                with self.assertRaises(ValueError):
                    amendment.canonical_json_bytes({"value": value})
                with self.assertRaises(ValueError):
                    amendment.semantic_sha256({"value": value})


if __name__ == "__main__":
    unittest.main()
