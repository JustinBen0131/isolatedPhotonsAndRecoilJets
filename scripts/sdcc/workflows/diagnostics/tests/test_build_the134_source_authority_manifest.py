#!/usr/bin/env python3
"""Adversarial tests for the deterministic THE-134 source manifest."""

from __future__ import annotations

import copy
import hashlib
import importlib.util
import tempfile
import unittest
from pathlib import Path


HERE = Path(__file__).resolve().parent
MODULE_PATH = HERE.parent / "build_the134_source_authority_manifest.py"
SPEC = importlib.util.spec_from_file_location("the134_source_builder", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
builder = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(builder)


class SourceAuthorityManifestTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name).resolve()
        self.sim_root = self.root / "sim_lists"
        self.sim_root.mkdir()
        for row in builder.inventory_rows():
            sample_root = self.sim_root / row["sample"]
            sample_root.mkdir()
            for role in builder.LIST_ROLES:
                path = sample_root / builder.LIST_FILENAMES[role]
                path.write_text(
                    "# frozen list\n"
                    f"/input/{row['sample']}/{role}.000.root\n"
                    "\n",
                    encoding="utf-8",
                )

    def tearDown(self) -> None:
        self.temporary.cleanup()

    def test_exact_thirteen_row_manifest_is_deterministic(self) -> None:
        first = builder.build_manifest(self.sim_root, "0mrad")
        second = builder.build_manifest(self.sim_root, "0mrad")
        self.assertEqual(
            builder.canonical_json_bytes(first),
            builder.canonical_json_bytes(second),
        )
        self.assertEqual(first["schema"], builder.SOURCE_SCHEMA)
        self.assertEqual(first["authority_state"], builder.AUTHORITY_STATE)
        self.assertEqual(first["row_count"], 13)
        self.assertEqual(len(first["rows"]), 13)
        self.assertEqual(
            first["global_five_tuple_check"],
            {
                "total": 13,
                "unique": 13,
                "duplicate_count": 0,
                "tuple_identity_records_sha256": first[
                    "global_five_tuple_check"
                ]["tuple_identity_records_sha256"],
            },
        )
        self.assertEqual(
            [row["row_id"] for row in first["rows"]],
            [row["row_id"] for row in builder.inventory_rows()],
        )
        pp_rows = [row for row in first["rows"] if row["system"] == "pp"]
        auau_rows = [row for row in first["rows"] if row["system"] == "auau"]
        self.assertTrue(
            all(
                row["source_period"] == "0mrad"
                and row["source_si_di_role"] == "SI"
                for row in pp_rows
            )
        )
        self.assertTrue(
            all(
                row["source_period"] == "AUAU_RUN24"
                and row["source_si_di_role"] == "EMBEDDED"
                for row in auau_rows
            )
        )

    def test_pinned_readback_passes_before_drift(self) -> None:
        payload = builder.build_manifest(self.sim_root, "0mrad")
        manifest = self.root / "source_authority.json"
        builder.atomic_write_json(manifest, payload)
        expected = hashlib.sha256(manifest.read_bytes()).hexdigest()
        observed = builder.load_manifest(manifest, expected)
        self.assertEqual(observed, payload)

    def test_mutable_source_path_drift_is_rejected(self) -> None:
        payload = builder.build_manifest(self.sim_root, "0mrad")
        manifest = self.root / "source_authority.json"
        builder.atomic_write_json(manifest, payload)
        expected = hashlib.sha256(manifest.read_bytes()).hexdigest()
        target = (
            self.sim_root
            / "run28_photonjet20"
            / builder.LIST_FILENAMES["calo_cluster"]
        )
        target.write_text(
            target.read_text(encoding="utf-8")
            + "/input/run28_photonjet20/calo_cluster.001.root\n",
            encoding="utf-8",
        )
        with self.assertRaisesRegex(builder.ManifestError, "hash drift"):
            builder.load_manifest(manifest, expected)

    def test_duplicate_row_is_rejected(self) -> None:
        payload = builder.build_manifest(self.sim_root, "0mrad")
        mutated = copy.deepcopy(payload)
        mutated["rows"][1] = copy.deepcopy(mutated["rows"][0])
        with self.assertRaisesRegex(builder.ManifestError, "duplicated"):
            builder.validate_manifest_payload(mutated, rehash=False)

    def test_partial_five_file_tuple_is_rejected(self) -> None:
        target = (
            self.sim_root
            / "run28_jet8"
            / builder.LIST_FILENAMES["g4hits"]
        )
        target.write_text("# frozen list\n\n\n", encoding="utf-8")
        with self.assertRaises(builder.ManifestError):
            builder.build_manifest(self.sim_root, "0mrad")

    def test_cross_source_duplicate_tuple_is_rejected(self) -> None:
        first_sample = self.sim_root / "run28_photonjet5"
        second_sample = self.sim_root / "run28_photonjet10"
        for role in builder.LIST_ROLES:
            first = first_sample / builder.LIST_FILENAMES[role]
            second = second_sample / builder.LIST_FILENAMES[role]
            second.write_bytes(first.read_bytes())
        with self.assertRaisesRegex(
            builder.ManifestError, "duplicate five-file source tuples"
        ):
            builder.build_manifest(self.sim_root, "0mrad")


if __name__ == "__main__":
    unittest.main()
