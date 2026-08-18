"""Regression guard: batch jobs must never discard their own error output.

An earlier version of Fun4All_recoilJets_unified_impl.C redirected BOTH
std::cout and std::cerr to /dev/null whenever RJ_VERBOSITY resolved to 0, and
the macro sets 0 automatically for any job with _CONDOR_SCRATCH_DIR or
_CONDOR_JOB_AD in its environment. Every [FATAL] raised on the farm therefore
went to /dev/null: jobs died with empty .err files and no recoverable reason.

That defect cost roughly three weeks of debugging failures that were invisible
rather than hard, and permanently destroyed the root cause of one cluster,
which could not be diagnosed afterwards because the message was never written
anywhere at all.

These tests fail loudly if any of the three conditions that produced it come
back:

  1. the macro redirects std::cerr,
  2. a worker stops forcing an explicit RJ_VERBOSITY, or
  3. a worker sets RJ_VERBOSITY=0.

Suppressing routine stdout volume in batch is legitimate. Discarding error
output is not, at any verbosity level.
"""

from __future__ import annotations

import re
import unittest
from pathlib import Path


REPO = next(p for p in Path(__file__).resolve().parents if (p / "AGENTS.md").is_file())
MACRO = REPO / "macros/Fun4All_recoilJets_unified_impl.C"
WORKERS = (
    REPO / "scripts/sdcc/runtime/condor/run_the121_the122_sim_production_row.sh",
    REPO / "scripts/sdcc/runtime/condor/run_the236_schema10_data_row.sh",
)
# Any C/C++ source that could reintroduce the same suppression.
SOURCE_GLOBS = ("macros/**/*.C", "src/**/*.cc", "src/**/*.h", "src_AuAu/**/*.cc", "src_AuAu/**/*.h")

CERR_REDIRECT = re.compile(r"cerr\s*\.\s*rdbuf\s*\(")


class NoStderrSuppressionTests(unittest.TestCase):
    def test_macro_never_redirects_stderr(self) -> None:
        self.assertTrue(MACRO.is_file(), f"missing macro: {MACRO}")
        text = MACRO.read_text(encoding="utf-8", errors="replace")
        hits = [
            f"{MACRO.name}:{i}"
            for i, line in enumerate(text.splitlines(), 1)
            if CERR_REDIRECT.search(line)
        ]
        self.assertEqual(
            hits, [],
            "std::cerr is redirected in the unified macro. Batch FATALs would be "
            f"discarded again. Offending lines: {hits}",
        )

    def test_macro_still_suppresses_stdout_volume(self) -> None:
        # The legitimate half of the behaviour must survive, otherwise a future
        # edit could 'fix' this by deleting the whole mechanism and flooding the
        # farm logs instead.
        text = MACRO.read_text(encoding="utf-8", errors="replace")
        self.assertIn(
            "cout.rdbuf", text,
            "the stdout volume control was removed; only the stderr redirect "
            "was ever the defect",
        )

    def test_no_source_file_redirects_stderr(self) -> None:
        offenders: list[str] = []
        for pattern in SOURCE_GLOBS:
            for path in REPO.glob(pattern):
                if not path.is_file():
                    continue
                body = path.read_text(encoding="utf-8", errors="replace")
                for i, line in enumerate(body.splitlines(), 1):
                    if CERR_REDIRECT.search(line):
                        offenders.append(f"{path.relative_to(REPO)}:{i}")
        self.assertEqual(
            offenders, [],
            f"std::cerr redirection reintroduced in: {offenders}",
        )

    def test_workers_force_explicit_verbosity_and_never_zero(self) -> None:
        for worker in WORKERS:
            if not worker.is_file():
                continue
            body = worker.read_text(encoding="utf-8", errors="replace")
            with self.subTest(worker=worker.name):
                self.assertRegex(
                    body, r"RJ_VERBOSITY=[1-9]",
                    f"{worker.name} must export a non-zero RJ_VERBOSITY so the "
                    "macro's Condor auto-quiet path cannot engage",
                )
                zeros = [
                    f"{worker.name}:{i}"
                    for i, line in enumerate(body.splitlines(), 1)
                    if re.search(r"RJ_VERBOSITY\s*=\s*0\b", line)
                ]
                self.assertEqual(
                    zeros, [],
                    f"{worker.name} sets RJ_VERBOSITY=0, which re-enables the "
                    f"batch silence path: {zeros}",
                )


if __name__ == "__main__":
    unittest.main()
