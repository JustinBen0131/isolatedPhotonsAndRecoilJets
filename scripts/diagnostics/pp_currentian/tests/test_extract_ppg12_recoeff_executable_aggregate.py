#!/usr/bin/env python3
"""Pure fail-closed tests for preserved-RecoEff aggregate evidence."""

from __future__ import annotations

import importlib.util
import csv
import hashlib
import sys
import tempfile
import unittest
from pathlib import Path


MODULE_PATH = Path(__file__).resolve().parents[1] / "extract_ppg12_recoeff_executable_aggregate.py"
SPEC = importlib.util.spec_from_file_location("recoeff_aggregate", MODULE_PATH)
assert SPEC and SPEC.loader
AGG = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = AGG
SPEC.loader.exec_module(AGG)


def row(identity: str, region: str, et: float, weight: float) -> dict[str, str]:
    output = {
        "candidate_identity": identity,
        "rj_cluster_Et": str(et),
        "rj_weight_final": str(weight),
        "rj_signal_fill_multiplicity": "1",
    }
    for name in AGG.REGION_HISTS:
        output[f"rj_signal_fill_{name}"] = str(int(name == region))
    return output


class TestRecoEffAggregate(unittest.TestCase):
    def test_content_sumw2_and_leakage_are_deterministic(self) -> None:
        rows = [
            row("a", "A", 10.0, 2.0),
            row("b", "A", 10.5, 3.0),
            row("c", "C", 10.5, 1.0),
            row("d", "D", 36.0, 4.0),
        ]
        cells = AGG.aggregate_recoil_rows(rows)
        self.assertEqual(cells["A"][1], {"content": 5.0, "sumw2": 13.0})
        self.assertEqual(cells["C"][1], {"content": 1.0, "sumw2": 1.0})
        self.assertEqual(cells["D"][-1], {"content": 4.0, "sumw2": 16.0})
        self.assertEqual(AGG.leakage_ratios(cells)["C_over_A"][1], 0.2)
        self.assertTrue(AGG.compare_cells(cells, cells)["pass"])

    def test_one_content_or_sumw2_mutation_fails(self) -> None:
        oracle = AGG.aggregate_recoil_rows([row("a", "A", 12.5, 2.0)])
        mutated_content = AGG.aggregate_recoil_rows([row("a", "A", 12.5, 2.001)])
        self.assertFalse(AGG.compare_cells(oracle, mutated_content)["pass"])
        mutated_sumw2 = {
            region: [dict(cell) for cell in values]
            for region, values in oracle.items()
        }
        mutated_sumw2["A"][2]["sumw2"] += 0.1
        self.assertFalse(AGG.compare_cells(oracle, mutated_sumw2)["pass"])

    def test_duplicate_or_multi_region_candidate_fails_closed(self) -> None:
        with self.assertRaises(AGG.ExtractFailure):
            AGG.aggregate_recoil_rows(
                [row("duplicate", "A", 20.0, 1.0), row("duplicate", "A", 20.0, 1.0)]
            )
        invalid = row("multi", "A", 20.0, 1.0)
        invalid["rj_signal_fill_B"] = "1"
        invalid["rj_signal_fill_multiplicity"] = "2"
        with self.assertRaises(AGG.ExtractFailure):
            AGG.aggregate_recoil_rows([invalid])

    def test_candidate_csv_is_cryptographically_bound_to_lane_and_contract(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            contract = base / "contract.json"
            contract.write_text('{"schema_version":3}\n')
            digest = hashlib.sha256(contract.read_bytes()).hexdigest()
            candidate = row("a", "A", 12.0, 1.0)
            candidate.update(
                lane_id="photon:photon5:0mrad:si",
                runtime_contract_sha256=digest,
            )
            csv_path = base / "candidates.csv"
            with csv_path.open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=list(candidate))
                writer.writeheader()
                writer.writerow(candidate)
            rows = AGG.validate_candidate_identity_binding(
                csv_path,
                lane_id="photon:photon5:0mrad:si",
                runtime_contract_sha256=digest,
            )
            self.assertEqual(len(rows), 1)
            with self.assertRaises(AGG.ExtractFailure):
                AGG.validate_candidate_identity_binding(
                    csv_path,
                    lane_id="photon:photon10:0mrad:si",
                    runtime_contract_sha256=digest,
                )
            with self.assertRaises(AGG.ExtractFailure):
                AGG.validate_candidate_identity_binding(
                    csv_path,
                    lane_id="photon:photon5:0mrad:si",
                    runtime_contract_sha256="0" * 64,
                )


if __name__ == "__main__":
    unittest.main()
