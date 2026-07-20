#!/usr/bin/env python3
"""Tests for the canonical stitched-purity numerical evidence producers."""

from __future__ import annotations

import copy
import hashlib
import importlib.util
import json
import tempfile
import unittest
from argparse import Namespace
from array import array
from pathlib import Path


REPO = Path(__file__).resolve().parents[4]
PRODUCER_PATH = (
    REPO
    / "scripts/diagnostics/pp_currentian/produce_ppg12_stitched_purity_evidence.py"
)
CONTRACT_PATH = (
    REPO / "agent_context/analysis_contracts/ppg12_stitched_purity_closure.yaml"
)


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


PRODUCER = load_module("produce_ppg12_stitched_purity_evidence", PRODUCER_PATH)
CONTRACT = json.loads(CONTRACT_PATH.read_text())
try:
    import ROOT  # noqa: F401
except ModuleNotFoundError:
    HAVE_ROOT = False
else:
    HAVE_ROOT = True


def file_sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical_purity() -> dict[str, object]:
    return {
        "bin_edges": [10.0, 12.0, 14.0],
        "truth": {"value": [0.70, 0.72], "error": [0.01, 0.01]},
        "raw": {"value": [0.60, 0.62], "error": [0.02, 0.02]},
        "corrected": {"value": [0.68, 0.70], "error": [0.025, 0.025]},
    }


class RepeatProducerTest(unittest.TestCase):
    def test_zero_region_a_still_consumes_leakage_draws_in_ppg12_order(self) -> None:
        class RecordingRng:
            def __init__(self) -> None:
                self.calls: list[tuple[str, float, float | None]] = []
                self.poisson_values = iter((0.0, 3.0, 8.0, 10.0))

            def PoissonD(self, mean: float) -> float:
                self.calls.append(("PoissonD", mean, None))
                return next(self.poisson_values)

            def Gaus(self, mean: float, sigma: float) -> float:
                self.calls.append(("Gaus", mean, sigma))
                return mean

        rng = RecordingRng()
        toy, leakage_toy = PRODUCER._draw_toy_state(
            rng,
            (1.0, 1.0, 1.0, 1.0),
            (0.001, 3.0, 8.0, 10.0),
            (0.1, 0.2, 0.3),
            (0.01, 0.02, 0.03),
        )
        self.assertEqual(toy[0], 0.0)
        self.assertEqual(leakage_toy, (0.1, 0.2, 0.3))
        self.assertEqual(
            [name for name, _, _ in rng.calls],
            ["PoissonD"] * 4 + ["Gaus"] * 3,
        )

    def test_repeated_payload_binds_identical_estimator_and_diagnostics(self) -> None:
        purity = canonical_purity()
        diagnostics = [{"bin": 1, "accepted": 20000.0}]

        def evaluator(inputs, contract, label):
            del inputs, contract, label
            return copy.deepcopy(purity), copy.deepcopy(diagnostics)

        payload = PRODUCER._build_repeated_purity_payload(
            {}, CONTRACT, evaluator=evaluator
        )
        expected = PRODUCER._payload_sha256(purity)
        self.assertEqual(
            payload["fixed_seed_repetition"]["first_output_sha256"], expected
        )
        self.assertEqual(
            payload["fixed_seed_repetition"]["repeated_output_sha256"], expected
        )
        self.assertEqual(
            payload["run_diagnostics"]["diagnostics_sha256"],
            payload["run_diagnostics"]["repeated_diagnostics_sha256"],
        )

    def test_changed_second_estimator_output_fails(self) -> None:
        calls = 0

        def evaluator(inputs, contract, label):
            nonlocal calls
            del inputs, contract, label
            calls += 1
            purity = canonical_purity()
            if calls == 2:
                purity["corrected"]["value"][0] += 1.0e-12
            return purity, [{"bin": 1, "accepted": 20000.0}]

        with self.assertRaisesRegex(PRODUCER.EvidenceError, "not byte-deterministic"):
            PRODUCER._build_repeated_purity_payload({}, CONTRACT, evaluator=evaluator)

    def test_changed_second_diagnostics_fails(self) -> None:
        calls = 0

        def evaluator(inputs, contract, label):
            nonlocal calls
            del inputs, contract, label
            calls += 1
            return canonical_purity(), [{"bin": 1, "accepted": 20000.0 + calls}]

        with self.assertRaisesRegex(PRODUCER.EvidenceError, "diagnostics"):
            PRODUCER._build_repeated_purity_payload({}, CONTRACT, evaluator=evaluator)


class HistoricalMetricTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.input_path = self.root / "candidate_purity.json"
        self.input_path.write_text("{}\n")
        self.source_paths: dict[str, Path] = {}
        for role in PRODUCER.REQUIRED_HISTORICAL_SOURCE_ROLES:
            path = self.root / f"{role}.dat"
            path.write_bytes(role.encode())
            self.source_paths[role] = path

    def tearDown(self) -> None:
        self.temp.cleanup()

    @staticmethod
    def _series(values: list[float]) -> dict[str, object]:
        return {
            "reference": {"value": values, "error": [0.10, 0.10]},
            "candidate": {"value": list(values), "error": [0.10, 0.10]},
        }

    def _source(self) -> dict[str, object]:
        return {
            "schema": "ppg12-stitched-purity-historical-input/v1",
            "final_purity_series": "corrected",
            "bin_edges": [10.0, 12.0, 14.0],
            "final_purity": self._series([0.70, 0.80]),
            "abcd": {
                region: self._series([10.0 + index, 11.0 + index])
                for index, region in enumerate(("A", "B", "C", "D"))
            },
            "leakage": {
                region: self._series([0.05 + index * 0.01, 0.06 + index * 0.01])
                for index, region in enumerate(("cB", "cC", "cD"))
            },
            "source_links": [
                {
                    "role": role,
                    "path": str(path),
                    "sha256": file_sha256(path),
                }
                for role, path in sorted(self.source_paths.items())
            ],
        }

    def test_unity_inputs_emit_explicit_zero_trends(self) -> None:
        payload = PRODUCER._historical_metrics_payload(
            self._source(), self.input_path, CONTRACT
        )
        self.assertEqual(payload["chi2_ndf"], 0.0)
        self.assertEqual(payload["max_abs_pull"], 0.0)
        self.assertEqual(payload["weighted_mean_ratio"], 1.0)
        self.assertEqual(payload["abcd_coherent_trend_max_sigma"], 0.0)
        self.assertEqual(payload["leakage_coherent_trend_max_sigma"], 0.0)
        self.assertEqual(set(payload["abcd"]), {"A", "B", "C", "D"})
        self.assertEqual(set(payload["leakage"]), {"cB", "cC", "cD"})

    def test_abcd_and_leakage_trends_are_independent(self) -> None:
        source = self._source()
        source["abcd"]["B"]["candidate"]["value"] = [12.1, 13.2]
        payload = PRODUCER._historical_metrics_payload(
            source, self.input_path, CONTRACT
        )
        self.assertGreater(payload["abcd_coherent_trend_max_sigma"], 0.0)
        self.assertEqual(payload["leakage_coherent_trend_max_sigma"], 0.0)

        source = self._source()
        source["leakage"]["cC"]["candidate"]["value"] = [0.08, 0.09]
        payload = PRODUCER._historical_metrics_payload(
            source, self.input_path, CONTRACT
        )
        self.assertEqual(payload["abcd_coherent_trend_max_sigma"], 0.0)
        self.assertGreater(payload["leakage_coherent_trend_max_sigma"], 0.0)

    def test_missing_source_role_and_hash_drift_fail(self) -> None:
        source = self._source()
        source["source_links"].pop()
        with self.assertRaisesRegex(PRODUCER.EvidenceError, "missing roles"):
            PRODUCER._historical_metrics_payload(source, self.input_path, CONTRACT)

        source = self._source()
        source["source_links"][0]["sha256"] = "0" * 64
        with self.assertRaisesRegex(PRODUCER.EvidenceError, "hash mismatch"):
            PRODUCER._historical_metrics_payload(source, self.input_path, CONTRACT)


@unittest.skipUnless(HAVE_ROOT, "PyROOT is unavailable")
class RootEstimatorIntegrationTest(unittest.TestCase):
    def test_zero_a_throw_preserves_the_following_root_rng_state(self) -> None:
        root = PRODUCER._load_root()
        helper_rng = root.TRandom3(42)
        oracle_rng = root.TRandom3(42)
        values = (1.0, 1.0, 1.0, 1.0)
        effective = (0.001, 3.0, 8.0, 10.0)
        leakage = (0.1, 0.2, 0.3)
        errors = (0.01, 0.02, 0.03)

        toy, _ = PRODUCER._draw_toy_state(
            helper_rng, values, effective, leakage, errors
        )
        self.assertEqual(toy[0], 0.0)

        for mean in effective:
            oracle_rng.PoissonD(mean)
        for mean, sigma in zip(leakage, errors):
            oracle_rng.Gaus(mean, sigma)

        helper_next = tuple(helper_rng.PoissonD(mean) for mean in effective)
        oracle_next = tuple(oracle_rng.PoissonD(mean) for mean in effective)
        self.assertEqual(helper_next, oracle_next)

    def test_exact_tf1_solver_and_full_fixed_seed_repeat(self) -> None:
        root = PRODUCER._load_root()
        solver = PRODUCER._new_root_solver(root, "unit_exact_solver", 100.0)
        self.assertAlmostEqual(
            PRODUCER._tf1_root(
                solver, (100.0, 20.0, 30.0, 10.0), (0.0, 0.0, 0.0)
            ),
            40.0,
            places=7,
        )
        inputs = {
            "bin_edges": [10.0, 12.0],
            "inclusive": {
                "A": ([10000.0], [10000.0]),
                "B": ([2000.0], [2000.0]),
                "C": ([3000.0], [3000.0]),
                "D": ([1000.0], [1000.0]),
                "A_signal": ([7000.0], [7000.0]),
                "A_notmatch": ([2000.0], [2000.0]),
            },
            "photon": {
                "A_signal": ([10000.0], [10000.0]),
                "B_signal": ([500.0], [500.0]),
                "C_signal": ([800.0], [800.0]),
                "D_signal": ([100.0], [100.0]),
            },
        }
        payload = PRODUCER._build_repeated_purity_payload(inputs, CONTRACT)
        self.assertEqual(
            payload["fixed_seed_repetition"]["first_output_sha256"],
            payload["fixed_seed_repetition"]["repeated_output_sha256"],
        )
        # Executable PPG12 truth is A_signal/A_unsuffixed = 0.70, not
        # A_signal/(A_signal+A_notmatch) = 0.777...
        self.assertAlmostEqual(payload["purity"]["truth"]["value"][0], 0.70, places=12)

    def test_direct_root_historical_extraction(self) -> None:
        root = PRODUCER._load_root()
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            purity_path = directory / "purity.json"
            purity = {
                "bin_edges": [10.0, 12.0],
                "truth": {"value": [0.7], "error": [0.01]},
                "raw": {"value": [0.6], "error": [0.02]},
                "corrected": {"value": [0.68], "error": [0.025]},
            }
            purity_sha = PRODUCER._payload_sha256(purity)
            purity_path.write_text(
                json.dumps(
                    {
                        "schema": "ppg12-stitched-purity-purity/v1",
                        "random_seed": 42,
                        "toy_count": 20000,
                        "purity": purity,
                        "fixed_seed_repetition": {
                            "first_output_sha256": purity_sha,
                            "repeated_output_sha256": purity_sha,
                        },
                    }
                )
            )

            final_root = directory / "historical_final.root"
            handle = root.TFile(str(final_root), "RECREATE")
            graph = root.TGraphErrors(1)
            graph.SetName("gpurity_leak")
            graph.SetPoint(0, 11.0, 0.68)
            graph.SetPointError(0, 1.0, 0.025)
            graph.Write()
            for name, value in (("B", 0.05), ("C", 0.08), ("D", 0.01)):
                hist = root.TH1D(f"h_leak_{name}", "", 1, array("d", [10.0, 12.0]))
                hist.SetBinContent(1, value)
                hist.SetBinError(1, 0.005)
                hist.Write()
            handle.Close()

            abcd_root = directory / "historical_abcd.root"
            handle = root.TFile(str(abcd_root), "RECREATE")
            abcd_values = {"A": 100.0, "B": 20.0, "C": 30.0, "D": 10.0}
            abcd_names = {
                "A": "h_tight_iso_cluster_0",
                "B": "h_tight_noniso_cluster_0",
                "C": "h_nontight_iso_cluster_0",
                "D": "h_nontight_noniso_cluster_0",
            }
            for region, name in abcd_names.items():
                hist = root.TH1D(name, "", 1, array("d", [10.0, 12.0]))
                hist.SetBinContent(1, abcd_values[region])
                hist.SetBinError(1, abcd_values[region] ** 0.5)
                hist.Write()
            handle.Close()

            inclusive_root = directory / "candidate_inclusive.root"
            handle = root.TFile(str(inclusive_root), "RECREATE")
            sim = handle.mkdir("SIM")
            sim.cd()
            for region, name in abcd_names.items():
                hist = root.TH1D(name, "", 1, array("d", [10.0, 12.0]))
                hist.SetBinContent(1, abcd_values[region])
                hist.SetBinError(1, abcd_values[region] ** 0.5)
                hist.Write()
            handle.Close()

            photon_root = directory / "candidate_photon.root"
            handle = root.TFile(str(photon_root), "RECREATE")
            sim = handle.mkdir("SIM")
            sim.cd()
            photon_values = {
                "h_tight_iso_cluster_signal_0": 100.0,
                "h_tight_noniso_cluster_signal_0": 5.0,
                "h_nontight_iso_cluster_signal_0": 8.0,
                "h_nontight_noniso_cluster_signal_0": 1.0,
            }
            for name, value in photon_values.items():
                hist = root.TH1D(name, "", 1, array("d", [10.0, 12.0]))
                hist.SetBinContent(1, value)
                hist.SetBinError(1, value ** 0.5)
                hist.Write()
            handle.Close()

            args = Namespace(
                candidate_purity=str(purity_path),
                historical_final_root=str(final_root),
                historical_abcd_root=str(abcd_root),
                candidate_inclusive_root=str(inclusive_root),
                candidate_photon_root=str(photon_root),
            )
            source, source_path = PRODUCER._historical_input_from_root(args)
            payload = PRODUCER._historical_metrics_payload(
                source, source_path, CONTRACT
            )
            self.assertEqual(payload["chi2_ndf"], 0.0)
            self.assertEqual(payload["abcd_coherent_trend_max_sigma"], 0.0)
            self.assertEqual(payload["leakage_coherent_trend_max_sigma"], 0.0)
            self.assertEqual(
                {row["role"] for row in payload["source_links"]},
                PRODUCER.REQUIRED_HISTORICAL_SOURCE_ROLES,
            )


if __name__ == "__main__":
    unittest.main()
