#!/usr/bin/env python3
"""Synthetic regression tests for the THE-112 sideband ranker."""

from __future__ import annotations

import importlib.util
import contextlib
import io
import json
import math
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np


MODULE_PATH = Path(__file__).resolve().parents[1] / "rank_the112_sideband_scan.py"
SPEC = importlib.util.spec_from_file_location("rank_the112_sideband_scan", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
ranker = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = ranker
SPEC.loader.exec_module(ranker)


X_EDGES = np.asarray([-0.4, -0.3, -0.2, -0.1, 0.0, 0.1], dtype=float)
Y_V2_EDGES = np.asarray([-1.0, 0.0, 1.0], dtype=float)
Y_LEGACY_EDGES = np.asarray([0.0, 4.0, 8.0], dtype=float)
Z_EDGES = np.asarray([15.0, 35.0], dtype=float)


def make_surface(
    category: str,
    counts: dict[str, float],
    *,
    coordinate_kind: str = "eiso_minus_cut",
    extra_wide: dict[str, float] | None = None,
) -> object:
    y_edges = Y_V2_EDGES if coordinate_kind == "eiso_minus_cut" else Y_LEGACY_EDGES
    shape = (X_EDGES.size + 1, y_edges.size + 1, Z_EDGES.size + 1)
    values = np.zeros(shape, dtype=float)

    # Flow bins are intentional: A uses a regular tight/isolated bin, B uses
    # score and isolation overflow, C uses isolation underflow, and D uses
    # isolation overflow.  A correct region sum must retain all four.
    values[5, 1, 1] = counts["A"]
    values[6, 3, 1] = counts["B"]
    values[3, 0, 1] = counts["C"]
    values[3, 3, 1] = counts["D"]
    if extra_wide:
        values[2, 0, 1] = extra_wide.get("C", 0.0)
        values[2, 3, 1] = extra_wide.get("D", 0.0)
    identity = ranker.SurfaceIdentity("cent_0_20", "isoR40_isSliding", category)
    surface = ranker.Surface(
        identity=identity,
        source_key=f"synthetic/{category}",
        values=values,
        variances=values.copy(),
        x_edges=X_EDGES,
        y_edges=y_edges,
        z_edges=Z_EDGES,
        coordinate_kind=coordinate_kind,
    )
    surface.validate()
    return surface


def mapping(surface: object) -> dict[object, object]:
    return {surface.identity: surface}


class SidebandRankerTest(unittest.TestCase):
    def setUp(self) -> None:
        self.window = ranker.ScanWindow(gap=0.10, lower_width=0.20)
        self.background = make_surface(
            "truthBackground", {"A": 100, "B": 20, "C": 50, "D": 10}
        )
        self.signal = make_surface("truthSignal", {"A": 100, "B": 10, "C": 20, "D": 2})
        # This is exactly the 50%-purity pseudo-data mixture of the two inputs.
        self.data = make_surface("all", {"A": 200, "B": 30, "C": 70, "D": 12})

    def test_flow_bins_and_sumw2_are_retained(self) -> None:
        counts = ranker.region_counts(self.background, self.window, [0], 4.0, 0.0)
        self.assertEqual([counts.by_name(name).sumw for name in "ABCD"], [100, 20, 50, 10])
        self.assertEqual([counts.by_name(name).sumw2 for name in "ABCD"], [100, 20, 50, 10])
        self.assertEqual(counts.B.neff, 20.0)
        self.assertEqual(counts.C.neff, 50.0)

    def test_legacy_raw_eiso_and_v2_delta_coordinates_agree(self) -> None:
        legacy = make_surface(
            "truthBackground",
            {"A": 100, "B": 20, "C": 50, "D": 10},
            coordinate_kind="raw_eiso",
        )
        modern = ranker.region_counts(self.background, self.window, [0], 4.0, 0.0)
        historical = ranker.region_counts(legacy, self.window, [0], 4.0, 0.0)
        self.assertEqual(modern, historical)

    def test_background_factorization_and_fixed_point_injection_close(self) -> None:
        data_counts = ranker.region_counts(self.data, self.window, [0], 4.0, 0.0)
        background_counts = ranker.region_counts(
            self.background, self.window, [0], 4.0, 0.0
        )
        signal_counts = ranker.region_counts(self.signal, self.window, [0], 4.0, 0.0)
        report = ranker.evaluate_cell(data_counts, background_counts, signal_counts, [0.5])
        self.assertAlmostEqual(report["background_predicted_A"], 100.0)
        self.assertAlmostEqual(report["background_factorization_relative_residual"], 0.0)
        injection = report["injected_purity_closure"][0]
        self.assertTrue(injection["converged"])
        self.assertTrue(injection["nonnegative_corrected_remainders"])
        self.assertAlmostEqual(injection["estimated_signal_A"], 100.0, places=6)
        self.assertAlmostEqual(injection["estimated_purity"], 0.5, places=8)
        self.assertAlmostEqual(injection["signal_yield_relative_residual"], 0.0, places=8)

    def test_window_uses_disjoint_truth_background_not_inclusive_all(self) -> None:
        inclusive_all = make_surface(
            "all", {"A": 999, "B": 1, "C": 1, "D": 999}
        )
        inclusive = {
            inclusive_all.identity: inclusive_all,
            self.background.identity: self.background,
        }
        report = ranker.evaluate_window(
            self.window,
            mapping(self.data),
            mapping(self.signal),
            inclusive,
            [0.5],
            4.0,
            0.0,
            ranker.RankingGates(
                min_sideband_neff=1.0,
                max_abs_background_factorization_residual=2.0,
                max_abs_injected_signal_yield_residual=2.0,
                max_prompt_sideband_over_tight=2.0,
            ),
        )
        cell = next(item for item in report["cells"] if item["pt_scope"] == "pt_15_35")
        self.assertEqual(cell["background_regions"]["A"]["sumw"], 100.0)
        self.assertAlmostEqual(cell["background_factorization_relative_residual"], 0.0)

    def test_deterministic_ranking_prefers_closing_narrow_window(self) -> None:
        background = make_surface(
            "truthBackground",
            {"A": 100, "B": 20, "C": 50, "D": 10},
            extra_wide={"C": 50, "D": 1},
        )
        signal = make_surface(
            "truthSignal",
            {"A": 100, "B": 10, "C": 20, "D": 2},
            extra_wide={"C": 10, "D": 1},
        )
        data = make_surface(
            "all",
            {"A": 200, "B": 30, "C": 70, "D": 12},
            extra_wide={"C": 60, "D": 2},
        )
        windows = [
            ranker.ScanWindow(gap=0.10, lower_width=0.30),
            ranker.ScanWindow(gap=0.10, lower_width=0.20),
        ]
        ranked = ranker.rank_windows(
            windows,
            mapping(data),
            mapping(signal),
            mapping(background),
            [0.5],
            4.0,
            0.0,
            ranker.RankingGates(
                min_sideband_neff=1.0,
                max_abs_background_factorization_residual=2.0,
                max_abs_injected_signal_yield_residual=2.0,
                max_prompt_sideband_over_tight=2.0,
            ),
        )
        self.assertEqual(ranked[0]["window"]["key"], "minus0.20_to_minus0.10")
        self.assertLess(
            ranked[0]["summary"]["worst_abs_background_factorization_residual"],
            ranked[1]["summary"]["worst_abs_background_factorization_residual"],
        )

    def test_gates_reject_nonfactorizing_or_statistically_empty_window(self) -> None:
        bad_background = make_surface(
            "truthBackground", {"A": 100, "B": 20, "C": 50, "D": 2}
        )
        report = ranker.evaluate_window(
            self.window,
            mapping(self.data),
            mapping(self.signal),
            mapping(bad_background),
            [0.5],
            4.0,
            0.0,
            ranker.RankingGates(
                min_sideband_neff=5.0,
                max_abs_background_factorization_residual=0.5,
                max_abs_injected_signal_yield_residual=0.5,
                max_prompt_sideband_over_tight=1.0,
            ),
        )
        self.assertFalse(report["passes_gates"])
        self.assertTrue(
            any("Neff" in failure or "factorization" in failure for failure in report["gate_failures"])
        )

    def test_surface_name_parsing_and_schema_aliases(self) -> None:
        modern_key = (
            "ALL/"
            + ranker.V2_TOKEN
            + "truthSignal_isoR30_isSliding_cent_0_20"
        )
        legacy_key = "ALL/" + ranker.LEGACY_TOKEN + "all_cent_20_50"
        modern = ranker._parse_identity(modern_key, "v2")
        legacy = ranker._parse_identity(legacy_key, "legacy")
        self.assertEqual(modern.view, "isoR30_isSliding")
        self.assertEqual(modern.category, "truthSignal")
        background_key = (
            "ALL/"
            + ranker.V2_TOKEN
            + "truthBackground_isoR40_isSliding_cent_20_50"
        )
        background = ranker._parse_identity(background_key, "v2")
        self.assertEqual(background.category, "truthBackground")
        self.assertEqual(legacy.view, "legacy_mixed_internal_views")
        self.assertEqual(ranker._schema_alias("the100"), "legacy")
        self.assertEqual(ranker._schema_alias("the112"), "v2")

    def test_schema_specific_gap_grid_excludes_half_bin_legacy_benchmark(self) -> None:
        self.assertEqual(
            ranker._parse_predeclared_subset(
                "0.02,0.04", ranker.LEGACY_EDGE_ALIGNED_GAPS, "legacy gaps"
            ),
            (0.02, 0.04),
        )
        with self.assertRaisesRegex(Exception, "outside predeclared grid"):
            ranker._parse_predeclared_subset(
                "0.03", ranker.LEGACY_EDGE_ALIGNED_GAPS, "legacy gaps"
            )

    def test_non_edge_window_is_rejected_instead_of_interpolated(self) -> None:
        non_edge = ranker.ScanWindow(gap=0.05, lower_width=0.20)
        with self.assertRaisesRegex(ValueError, "not a histogram edge"):
            ranker.region_counts(self.background, non_edge, [0], 4.0, 0.0)

    def test_json_cleanup_removes_nonstandard_nan(self) -> None:
        cleaned = ranker._finite_json_value({"x": math.nan, "y": [math.inf, 1.0]})
        self.assertEqual(cleaned, {"x": None, "y": [None, 1.0]})

    @unittest.skipIf(ranker.uproot is None, "uproot unavailable")
    def test_root_io_v2_can_nominate_but_legacy_is_advisory_only(self) -> None:
        def write_root(path: Path, schema: str, sample: str) -> None:
            if schema == "v2":
                x_edges = np.asarray([-0.20, -0.16, -0.02, 0.0, 0.10])
                y_edges = np.asarray([-1.0, 0.0, 1.0])
                views = ranker.V2_VIEWS
            else:
                x_edges = np.linspace(-0.8, 0.5, 66)
                y_edges = np.asarray([0.0, 4.0, 8.0])
                views = ("legacy",)
            with ranker.uproot.recreate(path) as root_file:
                if schema == "v2" and sample in ("signal", "inclusive"):
                    categories = (
                        "all",
                        "truthPrompt",
                        "truthSignal",
                        "truthBackground",
                    )
                elif sample == "signal":
                    categories = ("all", "truthSignal")
                else:
                    categories = ("all",)
                for centrality in ranker.CENTRALITY_TOKENS:
                    for view in views:
                        for category in categories:
                            counts = {"A": 200, "B": 30, "C": 70, "D": 12}
                            if sample == "inclusive" and category in ("all", "truthBackground"):
                                counts = {"A": 100, "B": 20, "C": 50, "D": 10}
                            if category == "truthSignal":
                                counts = {"A": 100, "B": 10, "C": 20, "D": 2}
                            array = np.zeros(
                                (x_edges.size - 1, y_edges.size - 1, 1), dtype=float
                            )
                            x_centers = 0.5 * (x_edges[:-1] + x_edges[1:])
                            y_centers = 0.5 * (y_edges[:-1] + y_edges[1:])
                            gap = 0.02 if schema == "v2" else 0.04
                            band = np.flatnonzero((x_centers > -0.16) & (x_centers < -gap))
                            tight = np.flatnonzero(x_centers > 0.0)
                            cut = 0.0 if schema == "v2" else 4.0
                            isolated = np.flatnonzero(y_centers < cut)
                            nonisolated = np.flatnonzero(y_centers > cut)
                            array[tight[0], isolated[-1], 0] = counts["A"]
                            array[tight[0], nonisolated[0], 0] = counts["B"]
                            array[band[0], isolated[-1], 0] = counts["C"]
                            array[band[0], nonisolated[0], 0] = counts["D"]
                            if schema == "v2":
                                name = (
                                    f"ALL/{ranker.V2_TOKEN}{category}_{view}_{centrality}"
                                )
                            else:
                                name = (
                                    f"ALL/{ranker.LEGACY_TOKEN}{category}_isoR30_{centrality}"
                                )
                            root_file[name] = (array, x_edges, y_edges, Z_EDGES)

        with tempfile.TemporaryDirectory(prefix="the112-ranker-test-") as temp:
            root = Path(temp)
            for schema in ("v2", "legacy"):
                inputs: dict[str, Path] = {}
                for sample in ("data", "signal", "inclusive"):
                    inputs[sample] = root / f"{schema}_{sample}.root"
                    write_root(inputs[sample], schema, sample)
                output_json = root / f"{schema}.json"
                output_csv = root / f"{schema}.csv"
                args = [
                    "--schema",
                    schema,
                    "--data",
                    str(inputs["data"]),
                    "--signal",
                    str(inputs["signal"]),
                    "--inclusive",
                    str(inputs["inclusive"]),
                    "--output-json",
                    str(output_json),
                    "--output-csv",
                    str(output_csv),
                    "--gaps",
                    "0.02" if schema == "v2" else "0.04",
                    "--lower-widths",
                    "0.16",
                    "--injected-purities",
                    "0.5",
                    "--min-neff",
                    "1",
                ]
                with contextlib.redirect_stdout(io.StringIO()):
                    self.assertEqual(ranker.main(args), 0)
                payload = json.loads(output_json.read_text())
                self.assertTrue(output_csv.is_file())
                if schema == "v2":
                    self.assertIsNotNone(payload["best_candidate"])
                    self.assertEqual(payload["canonical_eligible_candidate_count"], 1)
                    self.assertEqual(
                        payload["scan_contract"]["view_roles"]["isoR40_isSliding"],
                        "canonical",
                    )
                else:
                    self.assertIsNone(payload["best_candidate"])
                    self.assertIsNotNone(payload["advisory_best_score_band_prior"])
                    self.assertEqual(payload["canonical_eligible_candidate_count"], 0)


if __name__ == "__main__":
    unittest.main()
