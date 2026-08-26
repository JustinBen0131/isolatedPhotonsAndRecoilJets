#!/usr/bin/env python3
"""Fast fail-closed checks for the sparse-data publication validator."""

from __future__ import annotations

from pathlib import Path
import json
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
        check = self.source.index("shortfall=expected_processed-processed_events")
        publish = self.source.index('/bin/mv -- "$candidate_path" "$output_path"')
        self.assertLess(check, publish)
        self.assertIn('processed_count_contract="EXACT_EOF_V1"', self.source)
        self.assertIn('processed_count_contract="VERIFIED_PAIRED_DST_EOF_MINUS_ONE_V1"', self.source)
        self.assertIn('processed_count_contract="VERIFIED_PAIRED_DST_EOF_MINUS_TWO_V1"', self.source)
        self.assertIn('processed_count_contract="EXACT_BOUNDED_EVENT_RANGE_V1"', self.source)

    def _measurement_program(self) -> str:
        anchor = 'python3 - "$time_raw" "$measurement_candidate"'
        start = self.source.index(anchor)
        heredoc = self.source.index("<<'PY'\n", start) + len("<<'PY'\n")
        end = self.source.index("\nPY\n", heredoc)
        return self.source[heredoc:end]

    def _run_measurement(
        self,
        expected: int,
        processed: int,
        event_limit: int = 0,
        *,
        system: str = "pp",
        event_offset: int = 0,
        source_total_events: int | None = None,
    ) -> subprocess.CompletedProcess[str]:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            time_raw = root / "resource.time"
            candidate = root / "candidate.root.part"
            receipt = root / "measurement.json.part"
            gate_candidate = root / "gate.root.part"
            gate_output = root / "gate.root"
            gate_receipt_candidate = root / "gate.json.part"
            gate_receipt = root / "gate.json"
            time_raw.write_text("Maximum resident set size (kbytes): 250000\n", encoding="utf-8")
            candidate.write_bytes(b"production-shaped-root-candidate")
            gate_args = ["", "", "", ""]
            if system == "auau":
                gate_candidate.write_bytes(b"production-shaped-gate-candidate")
                gate_receipt_candidate.write_text(
                    json.dumps(
                        {
                            "schema": "AuAuInlineEventGateValidationReceiptV1",
                            "status": "PASS",
                        }
                    ),
                    encoding="utf-8",
                )
                gate_args = [
                    str(gate_candidate),
                    str(gate_output),
                    str(gate_receipt_candidate),
                    str(gate_receipt),
                ]
            summary = f"RJ_DATA_ROOT_CONTRACT_V1 processed={processed} retained=0 primary=0 extension=0"
            total = expected if source_total_events is None else source_total_events
            argv = [
                "python3", "-", str(time_raw), str(receipt), system, f"{system}_typical",
                "1", "0", "3072", str(candidate), str(root / "final.root"),
                "a" * 64, "b" * 64, str(expected), summary, str(event_limit),
                *gate_args, str(event_offset), str(total),
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

    def test_production_shaped_measurement_accepts_one_event_eof_shortfall(self) -> None:
        result = self._run_measurement(expected=1000, processed=999)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(result.receipt_exists)  # type: ignore[attr-defined]
        self.assertIn('"processed_event_shortfall":1', result.receipt_text)  # type: ignore[attr-defined]
        self.assertIn(
            '"processed_count_contract":"VERIFIED_PAIRED_DST_EOF_MINUS_ONE_V1"',
            result.receipt_text,  # type: ignore[attr-defined]
        )

    def test_production_shaped_measurement_accepts_two_event_eof_shortfall(self) -> None:
        result = self._run_measurement(expected=1000, processed=998)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(result.receipt_exists)  # type: ignore[attr-defined]
        self.assertIn(
            '"processed_count_contract":"VERIFIED_PAIRED_DST_EOF_MINUS_TWO_V1"',
            result.receipt_text,  # type: ignore[attr-defined]
        )

    def test_production_shaped_measurement_rejects_three_event_shortfall(self) -> None:
        result = self._run_measurement(expected=1000, processed=997)
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(result.receipt_exists)  # type: ignore[attr-defined]
        self.assertIn("processed-event count differs expected=1000 actual=997", result.stderr)

    def test_bounded_measurement_stays_exact(self) -> None:
        result = self._run_measurement(expected=1000, processed=999, event_limit=1000)
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(result.receipt_exists)  # type: ignore[attr-defined]
        self.assertIn("bounded processed-event count differs expected=1000 actual=999", result.stderr)

    def test_terminal_auau_range_accepts_proven_one_event_eof_shortfall(self) -> None:
        result = self._run_measurement(
            expected=35_000,
            processed=34_999,
            event_limit=35_000,
            system="auau",
            event_offset=65_000,
            source_total_events=100_000,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn(
            '"processed_count_contract":"VERIFIED_PAIRED_DST_TERMINAL_RANGE_EOF_MINUS_ONE_V1"',
            result.receipt_text,  # type: ignore[attr-defined]
        )

    def test_nonterminal_auau_range_rejects_eof_shortfall(self) -> None:
        result = self._run_measurement(
            expected=35_000,
            processed=34_999,
            event_limit=35_000,
            system="auau",
            event_offset=0,
            source_total_events=100_000,
        )
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(
            "bounded processed-event count differs expected=35000 actual=34999",
            result.stderr,
        )

    def test_production_shaped_measurement_rejects_overcount(self) -> None:
        result = self._run_measurement(expected=1000, processed=1001)
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(result.receipt_exists)  # type: ignore[attr-defined]
        self.assertIn("processed-event count differs expected=1000 actual=1001", result.stderr)

    def test_auau_companion_is_mandatory_and_base_is_commit_marker(self) -> None:
        self.assertIn('contract.get("schema")!="AuAuInlineEventGateProductionContractV1"', self.source)
        self.assertIn('contract.get("required_for_sparse_auau_data") is not True', self.source)
        self.assertIn('contract.get("photon10_scaled_bit")!=22', self.source)
        self.assertIn('export RJ_AUAU_EVENT_GATE_OUTPUT_CANDIDATE="$gate_candidate_path"', self.source)
        self.assertIn('export RJ_AUAU_EVENT_GATE_ROW_ID="$row_id"', self.source)
        self.assertIn('export RJ_AUAU_EVENT_GATE_SOURCE_PAIRS="$gate_source_pairs"', self.source)
        self.assertIn('mandatory AuAuEventGateV1 candidate is absent', self.source)
        gate_publish = self.source.index('/bin/mv -- "$gate_candidate_path" "$gate_output_path"')
        receipt_publish = self.source.index('/bin/mv -- "$measurement_candidate" "$measurement_path"')
        base_publish = self.source.index('/bin/mv -- "$candidate_path" "$output_path"')
        self.assertLess(gate_publish, receipt_publish)
        self.assertLess(receipt_publish, base_publish)

    def test_auau_dst_ttree_sharding_is_enforced_before_runtime(self) -> None:
        self.assertIn('"max_events_per_job":35000', self.source)
        self.assertIn('"workload_profile":"AUAU_DATA_DST_TO_TTREE"', self.source)
        self.assertIn('[[ "$event_limit" -gt 0 && "$event_limit" -le 35000 ]]', self.source)
        self.assertIn('[[ "$gate_source_pairs" -eq 1 ]]', self.source)
        self.assertIn('export RJ_EVENT_OFFSET="$event_offset"', self.source)



if __name__ == "__main__":
    unittest.main()
