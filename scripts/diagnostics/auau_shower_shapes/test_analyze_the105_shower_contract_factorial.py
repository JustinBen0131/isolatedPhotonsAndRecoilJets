#!/usr/bin/env python3
"""Synthetic end-to-end test for the THE-105 factorial analyzer."""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
import uproot


HERE = Path(__file__).resolve().parent
ANALYZER = HERE / "analyze_the105_shower_contract_factorial.py"
VARIANTS = ("historical", "towerinfo70", "canonical")
POPULATIONS = ("data", "signal", "inclusive")
SHOWER_BRANCHES = (
    "cluster_weta_cogx",
    "cluster_wphi_cogx",
    "e11_over_e33",
    "e32_over_e35",
    "cluster_et1",
    "cluster_et2",
    "cluster_et3",
    "cluster_et4",
)


def fixture_arrays(variant: str, population: str, rows: int = 72) -> dict[str, np.ndarray]:
    index = np.arange(rows)
    centrality = np.tile(np.array([5.0, 15.0, 25.0, 40.0, 55.0, 72.0]), rows // 6)
    population_shift = {"data": 0.00, "signal": -0.02, "inclusive": 0.03}[population]
    route_shift = 0.0
    if population != "data" and variant == "historical":
        route_shift = 0.05
    floor_shift = -0.025 if variant == "canonical" else 0.0
    base = 0.20 + 0.002 * (index % 20) + population_shift + route_shift + floor_shift

    # The candidate identity is deliberately invariant across variants.
    arrays: dict[str, np.ndarray] = {
        "run": np.full(rows, 50001, dtype=np.int32),
        "evt": index.astype(np.int64),
        "event_count": index.astype(np.int64),
        "is_sim": np.full(rows, population != "data", dtype=np.int32),
        "is_sim_embedded": np.full(rows, population != "data", dtype=np.int32),
        "sample_code": np.full(rows, {"data": 0, "signal": 12, "inclusive": 30}[population], dtype=np.int32),
        "cluster_Et": (16.0 + (index % 18)).astype(np.float64),
        "cluster_Eta": (-0.6 + 1.2 * (index % 24) / 23.0).astype(np.float64),
        "cluster_Phi": (-3.0 + 6.0 * (index % 36) / 35.0).astype(np.float64),
        "centrality": centrality.astype(np.float64),
        "event_weight": np.where(population == "data", 1.0, 0.7 + 0.1 * (index % 4)).astype(np.float64),
        "event_calo_total_energy": (100.0 + index).astype(np.float64),
        "reco_eiso": (0.1 * (index % 8)).astype(np.float64),
        "preselection_pass": ((index % 4) != 0).astype(np.int32),
        "baseline_bdt_tight": ((index % 5) >= 2).astype(np.int32),
        "active_tight_tag": ((index % 5) >= 2).astype(np.int32),
        "npb_pass": ((index % 4) != 0).astype(np.int32),
        "auau_tight_bdt_score": (0.45 + 0.01 * (index % 30)).astype(np.float64),
        "baseline_wp80_threshold": np.full(rows, 0.58, dtype=np.float64),
        "cluster_truth_barcode": np.where(population == "data", -1, 100000 + index).astype(np.int64),
    }
    arrays["cluster_weta_cogx"] = base.copy()
    arrays["cluster_wphi_cogx"] = (base * 1.08).copy()
    arrays["e11_over_e33"] = (0.72 + base * 0.35).copy()
    arrays["e32_over_e35"] = (0.78 + base * 0.28).copy()
    arrays["cluster_et1"] = (0.50 + base * 0.20).copy()
    arrays["cluster_et2"] = (0.30 + base * 0.15).copy()
    arrays["cluster_et3"] = (0.15 + base * 0.10).copy()
    arrays["cluster_et4"] = (0.04 + base * 0.04).copy()

    # Exercise explicit boundary and invalid-value accounting.
    arrays["cluster_weta_cogx"][0] = 0.0
    arrays["cluster_wphi_cogx"][1] = np.nan
    arrays["e11_over_e33"][2] = 1.3
    arrays["e32_over_e35"][3] = -0.1
    return arrays


def write_fixture(path: Path, variant: str, population: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with uproot.recreate(path) as root_file:
        root_file["analysis_config_yaml"] = (
            "analysis_config:\n"
            f"  cemc_shower_shape_diagnostic_variant: {variant}\n"
        )
        root_file["AuAuPhotonCandidateSkim"] = fixture_arrays(variant, population)


class FactorialAnalyzerTest(unittest.TestCase):
    def test_synthetic_factorial(self) -> None:
        with tempfile.TemporaryDirectory(prefix="the105-factorial-") as temp:
            root = Path(temp)
            inputs = root / "inputs"
            outputs = root / "outputs"
            for variant in VARIANTS:
                for population in POPULATIONS:
                    write_fixture(inputs / variant / population / "fixture.root", variant, population)

            completed = subprocess.run(
                [
                    sys.executable,
                    str(ANALYZER),
                    "--input-root",
                    str(inputs),
                    "--output-dir",
                    str(outputs),
                    "--strict",
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(completed.returncode, 0, msg=completed.stdout + completed.stderr)
            summary = json.loads((outputs / "validation_summary.json").read_text())
            self.assertEqual(summary["status"], "pass")
            self.assertEqual(len(summary["plots"]), 25)

            boundary = pd.read_csv(outputs / "boundary_and_range_fractions.csv")
            row = boundary[
                (boundary.population == "data")
                & (boundary.variant == "historical")
                & (boundary.centrality == "cent0_20")
                & (boundary.stage == "before")
                & (boundary.variable == "cluster_weta_cogx")
            ].iloc[0]
            self.assertGreater(row.exact_zero_fraction, 0.0)

            deltas = pd.read_csv(outputs / "matched_feature_deltas.csv")
            source = deltas[
                (deltas.population == "signal")
                & (deltas.comparison == "combined_historical_to_canonical")
                & (deltas.variable == "cluster_weta_cogx")
            ]
            self.assertTrue(np.all(source.mean_target_minus_source < 0.0))


if __name__ == "__main__":
    unittest.main()
