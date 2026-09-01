"""Focused pure tests for the THE-248 schema-10 stitching reducer."""

from __future__ import annotations

import contextlib
import io
import json
import math
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch


sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import reduce_the248_schema10_stitching as reducer  # noqa: E402
from reduce_the248_schema10_stitching import (  # noqa: E402
    BIN_COUNT,
    Event,
    TruthJet,
    TruthPhoton,
    merge_histograms,
    ownership_contains,
    reduce_records,
    sample_contract,
)


def with_input(payload: dict, path: Path, size: int = 10) -> dict:
    payload["inputs"] = [{"path": str(path.resolve()), "size_bytes": size}]
    payload["input_count"] = 1
    return payload


class TestSchema10StitchingReducer(unittest.TestCase):
    def test_all_fourteen_sample_contracts_and_inclusive_upper_edges(self) -> None:
        expected = {
            "pp_photon5": (0.0, 14.0),
            "pp_photon10": (14.0, 22.0),
            "pp_photon20": (22.0, math.inf),
            "pp_jet8": (9.0, 14.0),
            "pp_jet12": (14.0, 21.0),
            "pp_jet20": (21.0, 32.0),
            "pp_jet30": (32.0, 42.0),
            "pp_jet40": (42.0, 100.0),
            "auau_photon12": (12.0, 21.0),
            "auau_photon20": (21.0, math.inf),
            "auau_jet12": (12.0, 21.0),
            "auau_jet20": (21.0, 31.0),
            "auau_jet30": (31.0, 41.0),
            "auau_jet40": (41.0, math.inf),
        }
        self.assertEqual(set(expected), set(reducer.SAMPLE_CONTRACTS))
        for sample_id, (low, high) in expected.items():
            contract = sample_contract(sample_id)
            self.assertEqual((low, high), (contract.ownership_low_gev, contract.ownership_high_gev))
            self.assertTrue(ownership_contains(low, contract))
            if math.isfinite(high):
                self.assertTrue(ownership_contains(high, contract))

    def test_photon_maxima_weights_missing_and_fixed_binning(self) -> None:
        events = [Event("a", 2.0), Event("b", math.nan), Event("c", 0.5)]
        photons = [
            TruthPhoton("a", 5.0, 1),
            TruthPhoton("a", 5.5, 1),
            TruthPhoton("a", 99.0, 0),
            TruthPhoton("b", 14.0, 1),
            TruthPhoton("b", math.inf, 1),
            TruthPhoton("c", -1.0, 1),
        ]
        payload = reduce_records("pp_photon5", events, photons)
        histogram = payload["histogram"]
        self.assertEqual("THE248Schema10StitchingHistogramV1", payload["schema"])
        self.assertEqual("PASS", payload["status"])
        self.assertEqual(BIN_COUNT, len(histogram["raw_counts"]))
        self.assertEqual(401, len(payload["binning"]["edges_gev"]))
        self.assertEqual(1, histogram["raw_counts"][22])  # 5.5 / 0.25
        self.assertEqual(2.0, histogram["sumw"][22])
        self.assertEqual(4.0, histogram["sumw2"][22])
        self.assertEqual(1, histogram["raw_counts"][56])  # upper ownership edge 14 GeV
        self.assertEqual(1.0, histogram["sumw"][56])
        self.assertEqual(1, histogram["missing"]["raw_count"])
        self.assertEqual(3, payload["event_totals"]["generated_events"])
        self.assertEqual(2, payload["event_totals"]["events_with_observable"])
        self.assertEqual(1, payload["event_totals"]["events_with_weight_fallback"])
        self.assertEqual(2, payload["ownership_audit"]["events_in_window"])
        self.assertEqual(0, payload["ownership_audit"]["entries_dropped"])

    def test_jet_filter_max_and_ownership_audit_do_not_drop(self) -> None:
        events = [Event("outside", 1.5), Event("edge", 2.0), Event("overflow", 3.0)]
        jets = [
            TruthJet("outside", "antikt", 0.4, 8.0),
            TruthJet("outside", "kt", 0.4, 50.0),
            TruthJet("outside", "antikt", 0.3, 60.0),
            TruthJet("edge", "antikt", 0.4, 14.0),
            TruthJet("edge", "antikt", 0.4, math.nan),
            TruthJet("overflow", "antikt", 0.4, 100.0),
        ]
        payload = reduce_records("pp_jet8", events, jets)
        self.assertEqual(1, payload["histogram"]["raw_counts"][32])
        self.assertEqual(1, payload["histogram"]["raw_counts"][56])
        self.assertEqual(1, payload["histogram"]["overflow"]["raw_count"])
        self.assertEqual(1, payload["ownership_audit"]["events_in_window"])
        self.assertEqual(2, payload["ownership_audit"]["events_outside_window"])
        self.assertEqual(2, payload["event_totals"]["events_in_histogram_range"])
        self.assertEqual(1, payload["event_totals"]["events_overflow"])
        self.assertEqual(0, payload["ownership_audit"]["entries_dropped"])

    def test_pp_jet40_upper_edge_is_owned_but_histogram_overflow(self) -> None:
        payload = reduce_records(
            "pp_jet40",
            [Event("event", 1.0)],
            [TruthJet("event", "antikt", 0.4, 100.0)],
        )
        self.assertEqual(1, payload["ownership_audit"]["events_in_window"])
        self.assertEqual(1, payload["ownership_audit"]["events_on_upper_edge"])
        self.assertEqual(1, payload["histogram"]["overflow"]["raw_count"])

    def test_merge_is_additive_and_rejects_duplicate_inputs(self) -> None:
        first = with_input(
            reduce_records("auau_photon12", [Event("a", 2.0)], [TruthPhoton("a", 12.0, 1)]),
            Path("/tmp/the248-a.root"),
        )
        second = with_input(
            reduce_records("auau_photon12", [Event("b", 3.0)], [TruthPhoton("b", 20.0, 1)]),
            Path("/tmp/the248-b.root"),
        )
        merged = merge_histograms([first, second])
        self.assertEqual(2, merged["input_count"])
        self.assertEqual(2, merged["event_totals"]["generated_events"])
        self.assertEqual(2, sum(merged["histogram"]["raw_counts"]))
        self.assertEqual(5.0, sum(merged["histogram"]["sumw"]))
        with self.assertRaisesRegex(ValueError, "duplicate input path"):
            merge_histograms([first, first])
        other = reduce_records("auau_photon20", [Event("c", 1.0)], [TruthPhoton("c", 22.0, 1)])
        with self.assertRaisesRegex(ValueError, "incompatible stitching histogram field sample_id"):
            merge_histograms([first, other])

    def test_duplicate_events_or_orphan_truth_records_fail_closed(self) -> None:
        with self.assertRaisesRegex(ValueError, "duplicate event identity"):
            reduce_records("pp_photon10", [Event("e", 1.0), Event("e", 1.0)], [])
        with self.assertRaisesRegex(ValueError, "unknown event"):
            reduce_records("pp_photon10", [], [TruthPhoton("orphan", 15.0, 1)])

    def test_duplicate_cli_input_paths_fail_before_root_import(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "input.root"
            path.touch()
            with self.assertRaisesRegex(ValueError, "duplicate input path"):
                reducer.reduce_root_files("pp_photon5", [path, path])

    def test_missing_required_tree_fails_closed(self) -> None:
        class EmptyDirectory:
            def Get(self, _name):
                return None

        with self.assertRaisesRegex(RuntimeError, "missing required tree ReplayFoundationV1/RJEventV1"):
            reducer._required_tree(EmptyDirectory(), "RJEventV1")

    def test_stdout_mode_emits_only_canonical_json(self) -> None:
        payload = reduce_records("pp_photon20", [Event("event", 1.0)], [TruthPhoton("event", 25.0, 1)])
        with tempfile.TemporaryDirectory() as directory:
            input_path = Path(directory) / "input.root"
            input_path.touch()
            stdout = io.StringIO()
            with patch.object(reducer, "reduce_root_files", return_value=payload):
                with contextlib.redirect_stdout(stdout):
                    self.assertEqual(
                        0,
                        reducer.main(
                            ["--sample-id", "pp_photon20", "--input", str(input_path), "--output", "-"]
                        ),
                    )
            rendered = stdout.getvalue()
            self.assertEqual(payload, json.loads(rendered))
            self.assertEqual(reducer.canonical_json(payload), rendered)


if __name__ == "__main__":
    unittest.main()
