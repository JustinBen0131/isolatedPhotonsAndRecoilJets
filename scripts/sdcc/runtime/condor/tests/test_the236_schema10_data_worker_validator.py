#!/usr/bin/env python3
"""Fast fail-closed checks for the sparse-data publication validator."""

from pathlib import Path
import subprocess
import tempfile
import unittest


WORKER = Path(__file__).resolve().parents[1] / "run_the236_schema10_data_row.sh"


class SparseDataWorkerValidatorTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.source = WORKER.read_text(encoding="utf-8")

    def test_shell_syntax(self) -> None:
        result = subprocess.run(["bash", "-n", str(WORKER)], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_validator_uses_certified_pyroot_reader(self) -> None:
        self.assertNotIn("std::stoll", self.source)
        self.assertNotIn("root -b -q -l -e '", self.source)
        self.assertIn("import ROOT", self.source)
        self.assertIn('value.InheritsFrom("TNamed")', self.source)
        self.assertIn('value.InheritsFrom("TTree")', self.source)
        self.assertIn("parse_nonnegative", self.source)
        self.assertIn("invalid integer key=", self.source)
        self.assertIn('re.fullmatch(r"[0-9]+", raw)', self.source)

    def test_every_required_count_uses_checked_parser(self) -> None:
        for key in (
            "processed_events",
            "retained_events",
            "omitted_events",
            "retained_primary_events",
            "retained_extension_events",
            "retained_primary_candidates",
            "retained_extension_candidates",
        ):
            self.assertIn(f'parse_nonnegative("{key}")', self.source)
        self.assertIn('parse_nonnegative(f"{tree}_entries")', self.source)

    def test_zero_row_sparse_tables_remain_supported(self) -> None:
        self.assertIn("if not value:\n        return 0", self.source)
        self.assertIn("if expected != actual:", self.source)

    def test_root_validator_receives_only_root_path_after_setup(self) -> None:
        setup = self.source.index("source /opt/sphenix/core/bin/sphenix_setup.sh")
        validator = self.source.index('python3 - "$candidate_path" <<\'PY\'', setup)
        self.assertLess(setup, validator)

    def test_external_event_count_is_checked_before_publication(self) -> None:
        self.assertNotIn("RJ_EXPECT_PROCESSED", self.source)
        check = self.source.index("if processed_events!=expected_processed:")
        publish = self.source.index('/bin/mv -- "$candidate_path" "$output_path"')
        self.assertLess(check, publish)

    def _measurement_program(self) -> str:
        anchor = 'python3 - "$time_raw" "$measurement_candidate"'
        start = self.source.index(anchor)
        heredoc = self.source.index("<<'PY'\n", start) + len("<<'PY'\n")
        end = self.source.index("\nPY\n", heredoc)
        return self.source[heredoc:end]

    def _run_measurement(self, expected: int, processed: int) -> subprocess.CompletedProcess[str]:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            time_raw = root / "resource.time"
            candidate = root / "candidate.root.part"
            receipt = root / "measurement.json.part"
            time_raw.write_text("Maximum resident set size (kbytes): 250000\n", encoding="utf-8")
            candidate.write_bytes(b"production-shaped-root-candidate")
            summary = f"RJ_DATA_ROOT_CONTRACT_V1 processed={processed} retained=0 primary=0 extension=0"
            argv = [
                "python3", "-", str(time_raw), str(receipt), "pp", "pp_typical",
                "1", "0", "3072", str(candidate), str(root / "final.root"),
                "a" * 64, "b" * 64, str(expected), summary,
            ]
            result = subprocess.run(argv, input=self._measurement_program(), capture_output=True, text=True)
            result.receipt_exists = receipt.exists()  # type: ignore[attr-defined]
            result.receipt_text = receipt.read_text() if receipt.exists() else ""  # type: ignore[attr-defined]
            return result

    def test_production_shaped_measurement_accepts_exact_count(self) -> None:
        result = self._run_measurement(expected=1000, processed=1000)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(result.receipt_exists)  # type: ignore[attr-defined]
        self.assertIn('"processed_events":1000', result.receipt_text)  # type: ignore[attr-defined]

    def test_production_shaped_measurement_rejects_count_mismatch(self) -> None:
        result = self._run_measurement(expected=1000, processed=999)
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(result.receipt_exists)  # type: ignore[attr-defined]
        self.assertIn("processed-event count differs expected=1000 actual=999", result.stderr)


if __name__ == "__main__":
    unittest.main()
