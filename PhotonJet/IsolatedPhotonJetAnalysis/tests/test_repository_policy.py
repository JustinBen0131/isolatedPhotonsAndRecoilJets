from __future__ import annotations

import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import zipfile

import uproot


ROOT = Path(__file__).resolve().parents[1]


class RepositoryPolicyTest(unittest.TestCase):
    def test_public_boundary_and_plot_gate(self) -> None:
        command = [sys.executable, str(ROOT / "tools/validate_repository.py")]
        projection_receipt = ROOT / "PROJECTION_RECEIPT.json"
        if projection_receipt.is_file():
            command.extend(["--projection-receipt", str(projection_receipt)])
        completed = subprocess.run(
            command,
            cwd=ROOT,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stdout + completed.stderr)
        payload = json.loads(completed.stdout)
        self.assertEqual(payload["status"], "PASS")
        self.assertEqual(payload["findings"], [])

    def test_binary_metadata_gate_rejects_machine_paths(self) -> None:
        spec = importlib.util.spec_from_file_location(
            "photonjet_repository_validator", ROOT / "tools/validate_repository.py"
        )
        self.assertIsNotNone(spec)
        self.assertIsNotNone(spec.loader)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        with tempfile.TemporaryDirectory() as temporary:
            target = Path(temporary) / "payload.root"
            target.write_bytes(
                b"root\x00/private" + b"/" + b"tmp/example/output.root\x00/gpfs" + b"02/data/file.root"
            )
            findings = module._binary_sensitive_findings(target)
        self.assertIn("ephemeral path", findings)
        self.assertIn("site filesystem path", findings)

    def test_binary_metadata_gate_covers_account_path_families(self) -> None:
        spec = importlib.util.spec_from_file_location(
            "photonjet_repository_validator_paths", ROOT / "tools/validate_repository.py"
        )
        self.assertIsNotNone(spec)
        self.assertIsNotNone(spec.loader)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        probes = (
            b"/" + b"Users/alice/repository/file.root",
            b"/" + b"home/alice/repository/file.root",
            b"/" + b"sphenix/tg/alice/file.root",
        )
        with tempfile.TemporaryDirectory() as temporary:
            for index, probe in enumerate(probes):
                target = Path(temporary) / f"probe-{index}.root"
                target.write_bytes(b"root\x00" + probe)
                self.assertIn("absolute account path", module._binary_sensitive_findings(target))

    def test_task_identifier_gate_covers_separator_variants(self) -> None:
        spec = importlib.util.spec_from_file_location(
            "photonjet_repository_validator_tasks", ROOT / "tools/validate_repository.py"
        )
        self.assertIsNotNone(spec)
        self.assertIsNotNone(spec.loader)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        pattern = module.SENSITIVE_PATTERNS["task-style identifier"]
        prefix = "T" + "HE"
        for probe in (
            prefix + "-5",
            prefix + "49",
            prefix + "-259",
            prefix + "_259",
            prefix + " 259",
            prefix + "/259",
            prefix + "259",
            prefix.lower() + "-259",
            prefix.lower() + " 259",
            prefix.lower() + "/259",
            prefix.lower() + "263_workstream",
        ):
            self.assertIsNotNone(pattern.search(probe), probe)
        for public_token in (
            "SHA-256", "SHA256", "GPL-3.0", "GPL3", "ROOT-6", "UTF8", "BDT70", "CHI2_NDF"
        ):
            self.assertIsNone(pattern.search(public_token), public_token)

    def test_governance_and_provenance_vocabulary_is_gated(self) -> None:
        spec = importlib.util.spec_from_file_location(
            "photonjet_repository_validator_vocabulary", ROOT / "tools/validate_repository.py"
        )
        self.assertIsNotNone(spec)
        self.assertIsNotNone(spec.loader)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        governance = module.SENSITIVE_PATTERNS["workflow governance vocabulary"]
        provenance = module.SENSITIVE_PATTERNS["assistant or vendor provenance"]
        for probe in (
            "agent" + "-control infrastructure",
            "agent" + "_context/session",
            "control" + " plane",
            "agent " + "workstream",
        ):
            self.assertIsNotNone(governance.search(probe), probe)
        for probe in (
            "co" + "dex",
            "cl" + "aude",
            "chat" + "gpt",
            "open" + "ai",
            "g" + "pt-5",
            "a" + "i assistant",
            "large" + " language " + "model",
            "generated" + " with gpt",
        ):
            self.assertIsNotNone(provenance.search(probe), probe)

    def test_binary_gate_decodes_root_utf16_and_compressed_metadata(self) -> None:
        spec = importlib.util.spec_from_file_location(
            "photonjet_repository_validator_containers", ROOT / "tools/validate_repository.py"
        )
        self.assertIsNotNone(spec)
        self.assertIsNotNone(spec.loader)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        task = "T" + "HE-5"
        governance = "agent" + "_context/session"
        provenance = "co" + "dex"
        with tempfile.TemporaryDirectory() as temporary:
            work = Path(temporary)
            root_path = work / "metadata.root"
            with uproot.recreate(root_path) as root:
                root["metadata"] = uproot.writing.identify.to_TObjString(task)
            self.assertIn("task-style identifier", module._binary_sensitive_findings(root_path))

            encoded = work / "encoded.bin"
            encoded.write_bytes(governance.encode("utf-16le") + provenance.encode("utf-16be"))
            encoded_findings = module._binary_sensitive_findings(encoded)
            self.assertIn("workflow governance vocabulary", encoded_findings)
            self.assertIn("assistant or vendor provenance", encoded_findings)
            self.assertTrue(any(item.startswith("unsupported binary format") for item in encoded_findings))

            gzip_path = work / "metadata.txt.gz"
            with gzip.open(gzip_path, "wb") as stream:
                stream.write(governance.encode("utf-8"))
            self.assertIn("workflow governance vocabulary", module._binary_sensitive_findings(gzip_path))

            zip_path = work / "metadata.zip"
            with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
                archive.writestr("metadata.txt", provenance)
            self.assertIn("assistant or vendor provenance", module._binary_sensitive_findings(zip_path))

    def test_private_identifiers_and_nonportable_paths_are_gated(self) -> None:
        spec = importlib.util.spec_from_file_location(
            "photonjet_repository_validator_private", ROOT / "tools/validate_repository.py"
        )
        self.assertIsNotNone(spec)
        self.assertIsNotNone(spec.loader)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        probes = {
            "private credential material": (
                "-----BEGIN " + "PRIVATE " + "KEY-----",
                "gh" + "p_" + "A" * 24,
                "https" + "://alice:" + "example-password" + "@example.org/path",
            ),
            "private contact identifier": ("alice" + "@" + "example.org",),
            "scheduler or session identifier": ("cluster" + "-12345", "session" + " 9876"),
            "private host or account identifier": (
                "host" + "=" + "sub" + "mit-1.example.org",
            ),
            "nonportable URI or platform path": (
                "C:" + chr(92) + "data" + chr(92) + "file.root",
                "root" + "://private.example.org/file.root",
            ),
        }
        for label, values in probes.items():
            pattern = module.SENSITIVE_PATTERNS[label]
            for value in values:
                self.assertIsNotNone(pattern.search(value), value)

    def test_manifest_paths_are_safe_regular_and_tracked(self) -> None:
        spec = importlib.util.spec_from_file_location(
            "photonjet_repository_validator_manifest", ROOT / "tools/validate_repository.py"
        )
        self.assertIsNotNone(spec)
        self.assertIsNotNone(spec.loader)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        original_repo = module.REPO
        with tempfile.TemporaryDirectory() as temporary:
            base = Path(temporary)
            repository = base / "repository"
            repository.mkdir()
            (repository / "inside.bin").write_bytes(b"inside")
            outside = base / "outside.bin"
            outside.write_bytes(b"outside")
            (repository / "link.bin").symlink_to(outside)
            module.REPO = repository
            try:
                findings: list[str] = []
                self.assertIsNotNone(
                    module._safe_tracked_file(
                        "inside.bin", {"inside.bin"}, "test artifact", findings
                    )
                )
                self.assertIsNone(
                    module._safe_tracked_file(
                        str(outside), {"inside.bin"}, "test artifact", findings
                    )
                )
                self.assertIsNone(
                    module._safe_tracked_file(
                        "../outside.bin", {"inside.bin"}, "test artifact", findings
                    )
                )
                self.assertIsNone(
                    module._safe_tracked_file(
                        "inside.bin", set(), "test artifact", findings
                    )
                )
                self.assertIsNone(
                    module._safe_tracked_file(
                        "link.bin", {"link.bin"}, "test artifact", findings
                    )
                )
                self.assertGreaterEqual(len(findings), 4)
            finally:
                module.REPO = original_repo

    def test_release_matrix_claims_are_mechanically_bound(self) -> None:
        spec = importlib.util.spec_from_file_location(
            "photonjet_repository_validator_release", ROOT / "tools/validate_repository.py"
        )
        self.assertIsNotNone(spec)
        self.assertIsNotNone(spec.loader)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        findings: list[str] = []
        _, branch_total = module._validate_branch_contract(findings)
        self.assertEqual(branch_total, 314)
        module._validate_release_matrix(branch_total, findings)
        self.assertEqual(findings, [])

        original_repo = module.REPO
        with tempfile.TemporaryDirectory() as temporary:
            repository = Path(temporary)
            matrix = json.loads((ROOT / "RELEASE_MATRIX.json").read_text(encoding="utf-8"))
            matrix["reference_artifacts"]["pp_simulation_photonjet"]["sha256"] = "0" * 64
            (repository / "RELEASE_MATRIX.json").write_text(
                json.dumps(matrix), encoding="utf-8"
            )
            module.REPO = repository
            try:
                tampered: list[str] = []
                module._validate_release_matrix(314, tampered)
                self.assertTrue(any("reference hash differs" in item for item in tampered))
            finally:
                module.REPO = original_repo

    def test_projection_receipt_binds_each_target_surface(self) -> None:
        spec = importlib.util.spec_from_file_location(
            "photonjet_repository_validator_projection", ROOT / "tools/validate_repository.py"
        )
        self.assertIsNotNone(spec)
        self.assertIsNotNone(spec.loader)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        internal_extras = {
            "docs/internal_integration.md",
            "manifests/projections.json",
            "tools/project_release.py",
        }

        def write_receipt(repository: Path, target: str, paths: set[str]) -> Path:
            rows = []
            for relative in sorted(paths):
                output = repository / relative
                output.parent.mkdir(parents=True, exist_ok=True)
                output.write_text("synthetic projection row\n", encoding="utf-8")
                os.chmod(output, 0o644)
                rows.append(
                    {
                        "path": relative,
                        "mode": "100644",
                        "sha256": hashlib.sha256(output.read_bytes()).hexdigest(),
                        "size_bytes": output.stat().st_size,
                    }
                )
            content_digest = hashlib.sha256(
                json.dumps(rows, sort_keys=True, separators=(",", ":")).encode("utf-8")
            ).hexdigest()
            manifest = repository / "manifests/projections.json"
            manifest_digest = (
                hashlib.sha256(manifest.read_bytes()).hexdigest()
                if manifest.is_file()
                else "0" * 64
            )
            payload = {
                "schema": "PhotonJetProjectionReceiptV2",
                "target": target,
                "projection_path": {
                    "analysis": "PhotonJet/IsolatedPhotonJetAnalysis",
                    "internal": "canonical/photonjet",
                }[target],
                "source_commit": "0" * 40,
                "source_tree_clean": True,
                "projection_manifest_sha256": manifest_digest,
                "projection_content_sha256": content_digest,
                "file_count": len(rows),
                "files": rows,
            }
            receipt = repository / "PROJECTION_RECEIPT.json"
            receipt.write_text(json.dumps(payload), encoding="utf-8")
            return receipt

        original_repo = module.REPO
        try:
            for forbidden in internal_extras:
                with tempfile.TemporaryDirectory() as temporary:
                    repository = Path(temporary)
                    module.REPO = repository
                    receipt = write_receipt(repository, "analysis", {forbidden})
                    findings: list[str] = []
                    module._validate_projection_receipt(receipt, findings)
                    self.assertTrue(
                        any(item.endswith(forbidden) and "forbidden file" in item for item in findings),
                        findings,
                    )
            for missing in internal_extras:
                with tempfile.TemporaryDirectory() as temporary:
                    repository = Path(temporary)
                    module.REPO = repository
                    receipt = write_receipt(repository, "internal", internal_extras - {missing})
                    findings = []
                    module._validate_projection_receipt(receipt, findings)
                    self.assertIn(
                        f"required internal projection file is absent: {missing}", findings
                    )
        finally:
            module.REPO = original_repo

    def test_git_gate_scans_every_reachable_commit(self) -> None:
        spec = importlib.util.spec_from_file_location(
            "photonjet_repository_validator_history", ROOT / "tools/validate_repository.py"
        )
        self.assertIsNotNone(spec)
        self.assertIsNotNone(spec.loader)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        original_repo = module.REPO
        with tempfile.TemporaryDirectory() as temporary:
            repository = Path(temporary)
            subprocess.run(["git", "init", "-q", "-b", "main", str(repository)], check=True)
            subprocess.run(
                ["git", "-C", str(repository), "config", "user.name", "JustinBen0131"],
                check=True,
            )
            subprocess.run(
                [
                    "git", "-C", str(repository), "config", "user.email",
                    "146116947+JustinBen0131" + "@" + "users.noreply.github.com",
                ],
                check=True,
            )
            tracked = repository / "tracked.txt"
            tracked.write_text("Internal " + "T" + "HE-5" + " tree marker\n", encoding="utf-8")
            subprocess.run(["git", "-C", str(repository), "add", "tracked.txt"], check=True)
            subprocess.run(
                ["git", "-C", str(repository), "commit", "-q", "-m", "Neutral first commit"],
                check=True,
            )
            tracked.write_text("second\n", encoding="utf-8")
            subprocess.run(["git", "-C", str(repository), "add", "tracked.txt"], check=True)
            subprocess.run(
                ["git", "-C", str(repository), "commit", "-q", "-m", "Neutral tip"],
                check=True,
            )
            module.REPO = repository
            try:
                findings: list[str] = []
                module._validate_git_metadata(findings)
                self.assertTrue(any("reachable Git content" in item for item in findings), findings)
            finally:
                module.REPO = original_repo


if __name__ == "__main__":
    unittest.main()
