#!/usr/bin/env python3
"""Focused regression tests for bounded RecoilJets login-node preparation."""

from __future__ import annotations

import os
import stat
import subprocess
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[5]
SUBMITTER = REPO / "scripts/sdcc/runtime/condor/RecoilJets_Condor_submit.sh"


def shell_function(name: str, next_marker: str) -> str:
    source = SUBMITTER.read_text(encoding="utf-8")
    start = source.index(f"{name}() {{")
    end = source.index(next_marker, start)
    return source[start:end]


class LoginNodeSafetyTests(unittest.TestCase):
    def run_bash(
        self, body: str, *, environment: dict[str, str] | None = None
    ) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            ["bash", "-c", body],
            env={**os.environ, **(environment or {})},
            text=True,
            capture_output=True,
            check=False,
        )

    def test_group_materialization_is_byte_exact_and_single_process(self) -> None:
        function = shell_function("make_sim_groups", "# Dry-run job count for isSim")
        self.assertNotIn("split -l", function)
        self.assertNotIn('mv "$raw"', function)
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source.list"
            records = [f"record-{index}\n".encode() for index in range(9)]
            records.append(b"record-9")
            source.write_bytes(b"".join(records))
            body = f"""
set -euo pipefail
err() {{ printf 'ERROR: %s\\n' "$*" >&2; }}
say() {{ printf '%s\\n' "$*" >&2; }}
sim_init() {{ :; }}
SIM_CLEAN_LIST={source!s}
SIM_STAGE_DIR={root!s}
SIM_JOB_PREFIX=unit
{function}
make_sim_groups 3 0
"""
            result = self.run_bash(body)
            self.assertEqual(result.returncode, 0, result.stderr)
            paths = [Path(line) for line in result.stdout.splitlines()]
            self.assertEqual(
                [path.name for path in paths],
                [
                    "unit_grp001.list",
                    "unit_grp002.list",
                    "unit_grp003.list",
                    "unit_grp004.list",
                ],
            )
            self.assertEqual(
                [path.read_bytes() for path in paths],
                [b"".join(records[0:3]), b"".join(records[3:6]), b"".join(records[6:9]), records[9]],
            )
            for path in paths:
                self.assertEqual(stat.S_IMODE(path.stat().st_mode), 0o644)

    def test_group_hard_cap_fails_without_partial_files(self) -> None:
        function = shell_function("make_sim_groups", "# Dry-run job count for isSim")
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            source = root / "source.list"
            source.write_text("".join(f"row-{index}\n" for index in range(7)))
            body = f"""
set -uo pipefail
err() {{ printf 'ERROR: %s\\n' "$*" >&2; }}
say() {{ printf '%s\\n' "$*" >&2; }}
sim_init() {{ :; }}
SIM_CLEAN_LIST={source!s}
SIM_STAGE_DIR={root!s}
SIM_JOB_PREFIX=unit
RJ_LOGIN_NODE_MAX_GROUP_FILES_PER_ROW=2
{function}
make_sim_groups 2 0
"""
            result = self.run_bash(body)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("more than 2 group files", result.stderr)
            self.assertEqual(list(root.glob("unit_grp*.list")), [])

    def test_path_validation_default_is_bounded_and_hard_cap_rejects(self) -> None:
        function = shell_function(
            "validate_sim_clean_list_paths", "# Historical replay controls"
        )
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            existing = root / "input.root"
            existing.touch()
            source = root / "source.list"
            row = "\t".join([str(existing)] * 5) + "\n"
            source.write_text(row * 40, encoding="utf-8")
            common = f"""
err() {{ printf 'ERROR: %s\\n' "$*" >&2; }}
say() {{ printf '%s\\n' "$*" >&2; }}
sim_requires_global_lane() {{ return 1; }}
{function}
"""
            bounded = self.run_bash(
                common + f'validate_sim_clean_list_paths {source!s} 0\n'
            )
            self.assertEqual(bounded.returncode, 0, bounded.stderr)
            self.assertIn("validation capped at 32 line(s)", bounded.stderr)
            rejected = self.run_bash(
                common + f'validate_sim_clean_list_paths {source!s} 0\n',
                environment={"RJ_VALIDATE_SIM_INPUT_MAX_LINES": "257"},
            )
            self.assertEqual(rejected.returncode, 24)
            self.assertIn("exceeds login-node hard cap 256", rejected.stderr)

    def test_late_materialization_is_bounded(self) -> None:
        function = shell_function(
            "condor_late_materialization_block", "condor_worker_failure_hold_block"
        )
        body = f"""
set -euo pipefail
err() {{ printf 'ERROR: %s\\n' "$*" >&2; }}
{function}
condor_late_materialization_block
"""
        passed = self.run_bash(
            body,
            environment={
                "RJ_CONDOR_MAX_MATERIALIZE": "20",
                "RJ_CONDOR_MAX_IDLE": "5",
            },
        )
        self.assertEqual(passed.returncode, 0, passed.stderr)
        self.assertEqual(
            passed.stdout,
            "max_materialize = 20\nmax_idle = 5\n",
        )
        rejected = self.run_bash(
            body,
            environment={
                "RJ_CONDOR_MAX_MATERIALIZE": "257",
                "RJ_CONDOR_MAX_IDLE": "5",
            },
        )
        self.assertNotEqual(rejected.returncode, 0)
        self.assertIn("must be in [1,256]", rejected.stderr)


if __name__ == "__main__":
    unittest.main()
