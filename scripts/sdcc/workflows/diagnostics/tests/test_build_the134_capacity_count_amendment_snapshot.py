#!/usr/bin/env python3
"""Adversarial frozen-snapshot tests for the THE-134 count amendment."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import shutil
import tempfile
import unittest
from pathlib import Path
from typing import Callable


HERE = Path(__file__).resolve().parent
MODULE_PATH = HERE.parent / "build_the134_capacity_count_amendment.py"
SPEC = importlib.util.spec_from_file_location(
    "the134_capacity_count_amendment_snapshot_tests", MODULE_PATH
)
assert SPEC is not None and SPEC.loader is not None
amendment = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(amendment)


class FrozenSnapshotReceiptTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name).resolve()
        self.fixture_index = 0
        self.release = self.root / "release"
        (self.release / "lib").mkdir(parents=True)
        (self.release / "lib64").mkdir()
        self.release_alias = self.root / "release_alias"
        self.release_alias.symlink_to(self.release, target_is_directory=True)
        self.providers = {
            "libcalo_io.so": self._write(
                self.release / "lib64/libcalo_io.so", b"calo-io\n"
            ),
            "libclusteriso.so": self._write(
                self.release / "lib/libclusteriso.so", b"cluster-iso\n"
            ),
            "libjetbase.so": self._write(
                self.release / "lib/libjetbase.so", b"jet-base\n"
            ),
        }
        self.external_helper = self._write(
            self.root / "external/libhelper.so", b"external-helper\n"
        )
        self.source = {
            "header": self._write(
                self.root / "authority/PhotonClusterBuilder.h",
                b"builder-header\n",
            ),
            "calo": self._write(
                self.root / "authority/libcalo_reco.so",
                b"calo-reco\n",
            ),
            "pp": self._write(
                self.root / "authority/libRecoilJets.so",
                b"pp-analysis\n",
            ),
            "auau": self._write(
                self.root / "authority/libRecoilJetsAuAu.so",
                b"auau-analysis\n",
            ),
        }

    def tearDown(self) -> None:
        for path in sorted(
            self.root.rglob("*"), key=lambda item: len(item.parts), reverse=True
        ):
            if not path.is_symlink():
                try:
                    path.chmod(0o755 if path.is_dir() else 0o644)
                except (FileNotFoundError, OSError):
                    pass
        self.root.chmod(0o755)
        self.temporary.cleanup()

    @staticmethod
    def _sha(path: Path) -> str:
        return hashlib.sha256(path.read_bytes()).hexdigest()

    @staticmethod
    def _write(path: Path, content: bytes) -> Path:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(content)
        return path

    @staticmethod
    def _json_bytes(payload: object) -> bytes:
        return (
            json.dumps(
                payload,
                sort_keys=True,
                separators=(",", ":"),
                ensure_ascii=True,
            )
            + "\n"
        ).encode("utf-8")

    @staticmethod
    def _unseal(snapshot: Path) -> None:
        snapshot.chmod(0o755)
        for path in snapshot.rglob("*"):
            if not path.is_symlink():
                path.chmod(0o755 if path.is_dir() else 0o644)

    @staticmethod
    def _seal(snapshot: Path) -> None:
        for path in snapshot.rglob("*"):
            if not path.is_symlink():
                path.chmod(0o555 if path.is_dir() else 0o444)
        snapshot.chmod(0o555)

    def _manifest_entries(self, snapshot: Path) -> list[dict[str, object]]:
        entries: list[dict[str, object]] = []
        for path in sorted(
            snapshot.rglob("*"),
            key=lambda item: item.relative_to(snapshot).as_posix(),
        ):
            if path == snapshot / "snapshot_manifest.json":
                continue
            relative = path.relative_to(snapshot).as_posix()
            if path.is_symlink():
                entries.append(
                    {
                        "path": relative,
                        "sha256": self._sha(path.resolve(strict=True)),
                        "symlink_target": os.readlink(path),
                        "type": "symlink",
                    }
                )
            elif path.is_file():
                entries.append(
                    {
                        "path": relative,
                        "sha256": self._sha(path),
                        "size": path.stat().st_size,
                        "type": "file",
                    }
                )
            elif path.is_dir():
                entries.append({"path": relative, "type": "directory"})
            else:
                self.fail(f"unsupported test fixture path: {path}")
        return entries

    def _rewrite_manifest(self, fixture: dict[str, object]) -> None:
        snapshot = fixture["snapshot"]
        assert isinstance(snapshot, Path)
        manifest_path = snapshot / "snapshot_manifest.json"
        payload = {
            "schema": amendment.SNAPSHOT_MANIFEST_SCHEMA,
            "status": "PASS",
            "root": str(snapshot),
            "entries": self._manifest_entries(snapshot),
        }
        manifest_path.write_bytes(self._json_bytes(payload))
        receipt_row = fixture["receipt_row"]
        assert isinstance(receipt_row, dict)
        receipt_row["snapshot_manifest"] = str(manifest_path)
        receipt_row["snapshot_manifest_sha256"] = self._sha(manifest_path)

    def _fixture(
        self,
        system: str,
        *,
        extra_ldd_lines: list[str] | None = None,
        extra_snapshot_files: dict[str, bytes] | None = None,
        marker_only_ldd: bool = False,
        explicit_load_only: bool = False,
    ) -> dict[str, object]:
        self.fixture_index += 1
        snapshot = self.root / f"{system}_snapshot_{self.fixture_index}"
        snapshot_lib = snapshot / "lib"
        snapshot_lib.mkdir(parents=True)
        header = shutil.copy2(
            self.source["header"], snapshot / "PhotonClusterBuilder.h"
        )
        calo = shutil.copy2(
            self.source["calo"], snapshot_lib / "libcalo_reco.so"
        )
        soname = snapshot_lib / "libcalo_reco.so.1"
        soname.symlink_to("libcalo_reco.so")
        analysis_name = (
            "libRecoilJets.so"
            if system == "pp"
            else "libRecoilJetsAuAu.so"
        )
        analysis = shutil.copy2(
            self.source[system], snapshot_lib / analysis_name
        )
        for relative, content in (extra_snapshot_files or {}).items():
            self._write(snapshot / relative, content)

        aliased_providers = {
            "libcalo_io.so": self.release_alias / "lib64/libcalo_io.so",
            "libclusteriso.so": self.release_alias
            / "lib/libclusteriso.so",
            "libjetbase.so": self.release_alias / "lib/libjetbase.so",
        }
        expected_providers = {
            "libcalo_reco.so": calo,
            **aliased_providers,
        }
        target_paths = [
            analysis,
            calo,
            self.providers["libcalo_io.so"],
            self.providers["libclusteriso.so"],
            self.providers["libjetbase.so"],
        ]
        if marker_only_ldd:
            ldd_lines = [f"@@TARGET {path}" for path in target_paths]
        else:
            ldd_lines = [
                f"@@TARGET {target_paths[0]}",
                *(
                    ["linux-vdso.so.1 (0x0001)"]
                    if explicit_load_only
                    else [
                        f"libcalo_reco.so => {calo} (0x0001)",
                        (
                            "libcalo_io.so => "
                            f"{self.providers['libcalo_io.so']} (0x0002)"
                        ),
                        (
                            "libclusteriso.so => "
                            f"{self.providers['libclusteriso.so']} (0x0003)"
                        ),
                        (
                            "libjetbase.so => "
                            f"{self.providers['libjetbase.so']} (0x0004)"
                        ),
                    ]
                ),
                *(extra_ldd_lines or []),
            ]
            for index, path in enumerate(target_paths[1:], start=1):
                ldd_lines.extend(
                    (
                        f"@@TARGET {path}",
                        f"linux-vdso.so.1 (0x000{index + 4})",
                    )
                )
        ldd_path = self._write(
            snapshot / "snapshot_loader_receipt.ldd.txt",
            ("\n".join(ldd_lines) + "\n").encode("utf-8"),
        )

        providers: dict[str, dict[str, object]] = {}
        for family, path in expected_providers.items():
            provider = {
                "observed_resolutions": [str(path)],
                "realpath": str(path.resolve(strict=True)),
                "sha256": self._sha(path),
            }
            if family != "libcalo_reco.so":
                provider["declared_path"] = str(
                    self.providers[family]
                )
            providers[family] = provider
        loader_payload = {
            "schema": amendment.SNAPSHOT_LOADER_SCHEMA,
            "status": "PASS",
            "mode": system,
            "snapshot_lib": str(snapshot_lib),
            # The receipt uses the site alias while immutable authority uses
            # the real release roots; same-file identity must accept this.
            "release_roots": [
                str(self.release_alias / "lib64"),
                str(self.release_alias / "lib"),
            ],
            "providers": providers,
            "targets_inspected": 5,
            "pinned_calo_reco_release_companions": True,
            "ldd_report": {
                "path": "snapshot_loader_receipt.ldd.txt",
                "sha256": self._sha(ldd_path),
            },
        }
        loader_path = self._write(
            snapshot / "snapshot_loader_receipt.json",
            self._json_bytes(loader_payload),
        )

        immutable_artifacts = {
            "photon_cluster_builder_header": {
                "sha256": self._sha(self.source["header"])
            },
            "calo_reco_library": {
                "sha256": self._sha(self.source["calo"])
            },
            "release_calo_io": {
                "sha256": self._sha(self.providers["libcalo_io.so"])
            },
            "release_clusteriso": {
                "sha256": self._sha(self.providers["libclusteriso.so"])
            },
            "release_jetbase": {
                "sha256": self._sha(self.providers["libjetbase.so"])
            },
        }
        fixture: dict[str, object] = {
            "snapshot": snapshot,
            "receipt_row": {
                "snapshot_dir": str(snapshot),
                "snapshot_loader_receipt": str(loader_path),
                "snapshot_loader_receipt_sha256": self._sha(loader_path),
                "snapshot_manifest": str(
                    snapshot / "snapshot_manifest.json"
                ),
                "snapshot_manifest_sha256": "0" * 64,
                "snapshot_builder_header": str(header),
                "snapshot_builder_header_sha256": self._sha(header),
                "snapshot_calo_reco_library": str(calo),
                "snapshot_calo_reco_library_sha256": self._sha(calo),
                "snapshot_analysis_library": str(analysis),
                "snapshot_analysis_library_sha256": self._sha(analysis),
            },
            "manifest_row": {
                "photon_cluster_builder_header_sha256": self._sha(
                    self.source["header"]
                ),
                "calo_reco_library_sha256": self._sha(self.source["calo"]),
                "release_calo_io_sha256": self._sha(
                    self.providers["libcalo_io.so"]
                ),
                "release_clusteriso_sha256": self._sha(
                    self.providers["libclusteriso.so"]
                ),
                "release_jetbase_sha256": self._sha(
                    self.providers["libjetbase.so"]
                ),
                "library_sha256": self._sha(self.source[system]),
            },
            "immutable": {
                "runtime": {
                    "release_core_lib64_dir": str(self.release / "lib64"),
                    "release_core_lib_dir": str(self.release / "lib"),
                    "calo_reco_soname": "libcalo_reco.so.1",
                },
                "bundle_artifacts": immutable_artifacts,
            },
            "runtime_authority": {
                "providers": {
                    family: {
                        "path": str(path),
                        "sha256": self._sha(path),
                        "size_bytes": path.stat().st_size,
                    }
                    for family, path in aliased_providers.items()
                }
            },
        }
        self._rewrite_manifest(fixture)
        self._seal(snapshot)
        return fixture

    def _validate(self, fixture: dict[str, object], system: str) -> object:
        return amendment.validate_frozen_snapshot_receipts(
            f"{system}_fixture",
            system,
            fixture["receipt_row"],
            fixture["manifest_row"],
            fixture["immutable"],
            fixture["runtime_authority"],
        )

    def _mutate(
        self,
        fixture: dict[str, object],
        callback: Callable[[Path], None],
        *,
        rewrite_manifest: bool,
    ) -> None:
        snapshot = fixture["snapshot"]
        assert isinstance(snapshot, Path)
        self._unseal(snapshot)
        callback(snapshot)
        if rewrite_manifest:
            self._rewrite_manifest(fixture)
        self._seal(snapshot)

    def test_pp_and_auau_accept_sealed_copies_and_site_aliases(self) -> None:
        for system in ("pp", "auau"):
            with self.subTest(system=system):
                fixture = self._fixture(system)
                loader, manifest = self._validate(fixture, system)
                self.assertEqual(
                    loader["schema"], amendment.SNAPSHOT_LOADER_SCHEMA
                )
                self.assertEqual(
                    manifest["schema"], amendment.SNAPSHOT_MANIFEST_SCHEMA
                )

    def test_ldd_semantic_failures_are_rejected(self) -> None:
        cases = {
            "not_found": ["libmissing.so => not found"],
            "mutable_private": [
                "libbad.so => "
                "/sphenix/u/test/thesisAnalysis/install/lib/libbad.so "
                "(0x0005)"
            ],
            "staged_escape": [
                f"libhelper.so => {self.external_helper} (0x0006)"
            ],
            "swapped_provider": [
                "libcalo_io.so => "
                f"{self.providers['libclusteriso.so']} (0x0007)"
            ],
        }
        for name, lines in cases.items():
            with self.subTest(case=name):
                fixture = self._fixture(
                    "pp",
                    extra_ldd_lines=lines,
                    extra_snapshot_files=(
                        {"lib/libhelper.so": b"staged-helper\n"}
                        if name == "staged_escape"
                        else None
                    ),
                )
                with self.assertRaises(amendment.AmendmentError):
                    self._validate(fixture, "pp")

    def test_marker_only_ldd_report_is_rejected(self) -> None:
        fixture = self._fixture("pp", marker_only_ldd=True)
        with self.assertRaisesRegex(
            amendment.AmendmentError, "empty target section"
        ):
            self._validate(fixture, "pp")

    def test_explicit_load_companions_need_not_be_dt_needed(self) -> None:
        fixture = self._fixture("pp", explicit_load_only=True)
        self._validate(fixture, "pp")

    def test_noncanonical_ldd_report_alias_is_rejected(self) -> None:
        fixture = self._fixture("pp")

        def mutate(snapshot: Path) -> None:
            canonical = snapshot / "snapshot_loader_receipt.ldd.txt"
            alternate = snapshot / "ldd.actual.txt"
            canonical.rename(alternate)
            canonical.symlink_to(alternate.name)

        self._mutate(fixture, mutate, rewrite_manifest=True)
        with self.assertRaises(amendment.AmendmentError):
            self._validate(fixture, "pp")

    def test_noncanonical_manifest_hardlink_reference_is_rejected(self) -> None:
        fixture = self._fixture("pp")
        snapshot = fixture["snapshot"]
        assert isinstance(snapshot, Path)
        alias = snapshot / "manifest.alias.json"
        self._unseal(snapshot)
        os.link(snapshot / "snapshot_manifest.json", alias)
        receipt_row = fixture["receipt_row"]
        assert isinstance(receipt_row, dict)
        receipt_row["snapshot_manifest"] = str(alias)
        self._seal(snapshot)
        with self.assertRaises(amendment.AmendmentError):
            self._validate(fixture, "pp")

    def test_snapshot_dotdot_spelling_is_rejected(self) -> None:
        fixture = self._fixture("pp")
        alias_component = self.root / "alias_component"
        alias_component.mkdir()
        receipt_row = fixture["receipt_row"]
        snapshot = fixture["snapshot"]
        assert isinstance(receipt_row, dict)
        assert isinstance(snapshot, Path)
        receipt_row["snapshot_dir"] = str(
            alias_component / ".." / snapshot.name
        )
        with self.assertRaisesRegex(
            amendment.AmendmentError, "absolute real directory"
        ):
            self._validate(fixture, "pp")

    def test_manifest_inventory_detects_hardlink_alias(self) -> None:
        fixture = self._fixture("pp")

        def mutate(snapshot: Path) -> None:
            os.link(
                snapshot / "snapshot_loader_receipt.json",
                snapshot / "loader.hardlink.json",
            )

        self._mutate(fixture, mutate, rewrite_manifest=False)
        with self.assertRaises(amendment.AmendmentError):
            self._validate(fixture, "pp")

    def test_external_hardlink_to_snapshot_file_is_rejected(self) -> None:
        fixture = self._fixture("pp")
        snapshot = fixture["snapshot"]
        assert isinstance(snapshot, Path)
        self._unseal(snapshot)
        os.link(
            snapshot / "PhotonClusterBuilder.h",
            self.root / "external_builder_hardlink.h",
        )
        self._seal(snapshot)
        with self.assertRaisesRegex(
            amendment.AmendmentError, "linked snapshot file"
        ):
            self._validate(fixture, "pp")

    def test_forbidden_snapshot_local_provider_is_rejected(self) -> None:
        fixture = self._fixture(
            "pp",
            extra_snapshot_files={"lib/libcalo_io.so": b"local-provider\n"},
        )
        with self.assertRaises(amendment.AmendmentError):
            self._validate(fixture, "pp")

    def test_dangling_symlink_is_controlled_amendment_error(self) -> None:
        fixture = self._fixture("pp")

        def mutate(snapshot: Path) -> None:
            (snapshot / "dangling").symlink_to("missing-target")

        self._mutate(fixture, mutate, rewrite_manifest=False)
        with self.assertRaises(amendment.AmendmentError):
            self._validate(fixture, "pp")


if __name__ == "__main__":
    unittest.main()
