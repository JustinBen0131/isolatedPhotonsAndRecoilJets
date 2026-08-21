#!/usr/bin/env python3
"""Regression tests for the bounded shared-SDCC permission contract."""

from __future__ import annotations

import importlib.util
import re
import tempfile
import unittest
from pathlib import Path


HERE = Path(__file__).resolve().parent
GUARD_PATH = HERE.parent / "sdcc_scientific_permission_guard.py"
SPEC = importlib.util.spec_from_file_location("sdcc_permission_guard", GUARD_PATH)
assert SPEC is not None and SPEC.loader is not None
guard = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(guard)


class ScientificPermissionGuardTests(unittest.TestCase):
    def test_private_ancestor_and_file_fail_then_shared_modes_pass(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary) / "science"
            root.mkdir(mode=0o700)
            root.chmod(0o700)
            evidence = root / "evidence"
            evidence.mkdir(mode=0o700)
            evidence.chmod(0o700)
            receipt = evidence / "receipt.json"
            receipt.write_text("{}\n", encoding="utf-8")
            receipt.chmod(0o600)
            failed = guard.inspect_paths(root, [receipt])
            self.assertEqual(failed["status"], "FAIL")
            self.assertGreaterEqual(len(failed["failures"]), 2)
            root.chmod(0o2755)
            evidence.chmod(0o2755)
            receipt.chmod(0o644)
            passed = guard.inspect_paths(root, [receipt])
            self.assertEqual(passed["status"], "PASS")
            self.assertFalse(passed["tree_walk_performed"])

    def test_path_limit_prevents_high_cardinality_preflight(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            with self.assertRaisesRegex(
                guard.PermissionContractError,
                "bounded limit",
            ):
                guard.inspect_paths(
                    root,
                    [root / f"item_{index}" for index in range(guard.MAX_PATHS + 1)],
                    allow_missing_leaf=True,
                )

    def test_sdcc_runtime_sources_do_not_reintroduce_private_modes(self) -> None:
        sdcc_root = HERE.parents[2]
        forbidden = re.compile(
            r"umask\s+0?77|"
            r"chmod\s+0?(?:400|500|600|700)\b|"
            r"\.chmod\(0o(?:400|500|600|700)\)|"
            r"mode=0o700|"
            r"os\.open\([^\n]*0o600"
        )
        findings: list[str] = []
        for path in sorted(sdcc_root.rglob("*")):
            if not path.is_file() or "tests" in path.parts:
                continue
            if path.suffix not in {".py", ".sh"}:
                continue
            for line_number, line in enumerate(
                path.read_text(encoding="utf-8", errors="replace").splitlines(),
                start=1,
            ):
                if forbidden.search(line):
                    findings.append(f"{path.relative_to(sdcc_root)}:{line_number}:{line}")
        self.assertEqual(findings, [])


if __name__ == "__main__":
    unittest.main()
