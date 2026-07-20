#!/usr/bin/env python3
"""Rank bounded Au+Au photon-ID sidebands without looking at recoil results.

The scan uses candidate-level three-dimensional histograms whose axes are the
BDT score relative to the tight WP80 threshold, reconstructed isolation, and
photon transverse momentum.  It supports both:

* the historical THE-100 surface, which stored raw isolation and mixed the
  internal isolation views in one histogram; and
* the corrected THE-112 V2 surface, which stores ``Eiso - EisoCut`` for the
  canonical centrality-dependent sliding-isolation definition, separately for
  the nominal R=0.4 cone and the R=0.3 robustness view.

THE-100 output is useful only as a historical prior because its views cannot
be unmixed after filling.  Nominal selection must be based on V2 output.

The ranker is deliberately blind to reconstructed recoil jets and unfolded
``xJgamma``.  It reads only the named candidate-level TH3 families.  Candidate
windows are predeclared so a result cannot drive an unbounded scan.  For every
centrality, isolation view, and photon-pT bin it records:

* flow-aware ABCD sums and Sumw2;
* source-role background-only factorization from embedded inclusive simulation;
* prompt leakage fractions from the ``truthSignal`` surface;
* fixed-point closure after injecting known prompt purities into background;
* data and simulation effective entries in the non-tight regions; and
* non-negative corrected-region checks.

The output is a detailed JSON packet and a one-row-per-window CSV.  Passing a
scan gate identifies a canary shortlist; it does not by itself promote a
sideband or authorize a broad production.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

try:
    import uproot
except ModuleNotFoundError:  # Pure unit tests do not require ROOT I/O.
    uproot = None  # type: ignore[assignment]


SCHEMA_VERSION = "THE112_AUAU_SIDEBAND_RANKING_V1"
V2_TOKEN = "h3_auauSidebandScanV2_scoreMinusT80_vs_EisoMinusCut_vs_pT_"
LEGACY_TOKEN = "h3_auauDualView_scoreMinusT80_vs_Eiso_vs_pT_"
CENTRALITY_TOKENS = ("cent_0_20", "cent_20_50", "cent_50_80")
V2_VIEWS = (
    "isoR40_isSliding",
    "isoR30_isSliding",
)
V2_VIEW_ROLES = {
    "isoR40_isSliding": "canonical",
    "isoR30_isSliding": "robustness",
}

# This grid was frozen before looking at THE-112 data.  V2 can represent the
# historical (-0.20,-0.03) benchmark exactly; THE-100 cannot, so its advisory
# scan uses only the edge-aligned subset.  The grid also contains the
# score-only (-0.28,-0.06) prior and nearby systematic candidates.
PREDECLARED_GAPS = (0.02, 0.03, 0.04, 0.05, 0.06, 0.08)
LEGACY_EDGE_ALIGNED_GAPS = (0.02, 0.04, 0.06, 0.08)
PREDECLARED_LOWER_WIDTHS = (0.16, 0.20, 0.24, 0.28, 0.32, 0.40)
PREDECLARED_INJECTED_PURITIES = (0.30, 0.50, 0.70, 0.90)
AXIS_TOLERANCE = 1.0e-10


@dataclass(frozen=True, order=True)
class SurfaceIdentity:
    centrality: str
    view: str
    category: str


@dataclass
class Surface:
    """A flow-inclusive TH3 payload independent of the ROOT implementation."""

    identity: SurfaceIdentity
    source_key: str
    values: np.ndarray
    variances: np.ndarray
    x_edges: np.ndarray
    y_edges: np.ndarray
    z_edges: np.ndarray
    coordinate_kind: str

    def validate(self) -> None:
        expected = (
            self.x_edges.size + 1,
            self.y_edges.size + 1,
            self.z_edges.size + 1,
        )
        if self.values.shape != expected:
            raise ValueError(
                f"{self.source_key}: flow-inclusive shape {self.values.shape} != {expected}"
            )
        if self.variances.shape != expected:
            raise ValueError(
                f"{self.source_key}: variance shape {self.variances.shape} != {expected}"
            )
        for axis_name, edges in (
            ("score", self.x_edges),
            ("isolation", self.y_edges),
            ("photon-pT", self.z_edges),
        ):
            if edges.ndim != 1 or edges.size < 2 or not np.all(np.isfinite(edges)):
                raise ValueError(f"{self.source_key}: invalid {axis_name} edges")
            if not np.all(np.diff(edges) > 0.0):
                raise ValueError(f"{self.source_key}: non-increasing {axis_name} edges")
        if not np.all(np.isfinite(self.values)) or np.any(self.values < -1.0e-12):
            raise ValueError(f"{self.source_key}: non-finite or negative bin contents")
        if not np.all(np.isfinite(self.variances)) or np.any(self.variances < -1.0e-12):
            raise ValueError(f"{self.source_key}: non-finite or negative Sumw2")
        if self.coordinate_kind not in ("raw_eiso", "eiso_minus_cut"):
            raise ValueError(f"{self.source_key}: unknown isolation coordinate")


@dataclass(frozen=True)
class WeightedCount:
    sumw: float
    sumw2: float

    @property
    def neff(self) -> float:
        if self.sumw2 <= 0.0:
            return 0.0
        return self.sumw * self.sumw / self.sumw2


@dataclass(frozen=True)
class RegionCounts:
    A: WeightedCount
    B: WeightedCount
    C: WeightedCount
    D: WeightedCount

    def by_name(self, name: str) -> WeightedCount:
        return getattr(self, name)


@dataclass(frozen=True)
class ScanWindow:
    gap: float
    lower_width: float

    @property
    def lower_offset(self) -> float:
        return -self.lower_width

    @property
    def upper_offset(self) -> float:
        return -self.gap

    @property
    def key(self) -> str:
        return f"minus{self.lower_width:.2f}_to_minus{self.gap:.2f}"


@dataclass
class RankingGates:
    min_sideband_neff: float = 25.0
    max_abs_background_factorization_residual: float = 0.50
    max_abs_injected_signal_yield_residual: float = 0.50
    max_prompt_sideband_over_tight: float = 1.00


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _clean_keys(root_file: Any) -> list[str]:
    return sorted(str(key).split(";", 1)[0] for key in root_file.keys(recursive=True))


def _trigger_matches(key: str, trigger: str) -> bool:
    parts = key.split("/")
    return len(parts) > 1 and trigger in parts[:-1]


def _parse_identity(key: str, schema: str) -> SurfaceIdentity | None:
    token = V2_TOKEN if schema == "v2" else LEGACY_TOKEN
    if token not in key:
        return None
    tail = key.split(token, 1)[1]
    category = next(
        (
            item
            for item in ("truthBackground", "truthSignal", "truthPrompt", "all")
            if tail.startswith(item)
        ),
        None,
    )
    centrality = next((item for item in CENTRALITY_TOKENS if item in key), None)
    if category is None or centrality is None:
        return None
    if schema == "v2":
        view = next((item for item in V2_VIEWS if item in key), None)
        if view is None:
            return None
    else:
        view = "legacy_mixed_internal_views"
    return SurfaceIdentity(centrality=centrality, view=view, category=category)


def _detect_schema(keys_by_input: Mapping[str, Sequence[str]], requested: str) -> str:
    presence = {
        "v2": {label: any(V2_TOKEN in key for key in keys) for label, keys in keys_by_input.items()},
        "legacy": {
            label: any(LEGACY_TOKEN in key for key in keys) for label, keys in keys_by_input.items()
        },
    }
    if requested != "auto":
        missing = [label for label, present in presence[requested].items() if not present]
        if missing:
            raise ValueError(f"schema {requested} surfaces missing from: {', '.join(missing)}")
        return requested
    for candidate in ("v2", "legacy"):
        if all(presence[candidate].values()):
            return candidate
    raise ValueError(f"no common surface schema across inputs: {presence}")


def _surface_from_uproot(
    histogram: Any,
    identity: SurfaceIdentity,
    key: str,
    coordinate_kind: str,
    allow_poisson_variance_fallback: bool,
) -> Surface:
    values = np.asarray(histogram.values(flow=True), dtype=np.float64)
    variance_payload = histogram.variances(flow=True)
    if variance_payload is None:
        if not allow_poisson_variance_fallback:
            raise ValueError(
                f"{key}: no TH3 Sumw2; rerun with proper Sumw2 rather than inferring weighted errors"
            )
        variances = np.clip(values, 0.0, None)
    else:
        variances = np.asarray(variance_payload, dtype=np.float64)
    surface = Surface(
        identity=identity,
        source_key=key,
        values=values,
        variances=variances,
        x_edges=np.asarray(histogram.axis(0).edges(), dtype=np.float64),
        y_edges=np.asarray(histogram.axis(1).edges(), dtype=np.float64),
        z_edges=np.asarray(histogram.axis(2).edges(), dtype=np.float64),
        coordinate_kind=coordinate_kind,
    )
    surface.validate()
    return surface


def load_surfaces(
    path: Path,
    schema: str,
    trigger: str,
    allow_poisson_variance_fallback: bool,
) -> dict[SurfaceIdentity, Surface]:
    if uproot is None:
        raise RuntimeError(
            "uproot is unavailable; use /Users/patsfan753/Desktop/analysis/env/bin/python3"
        )
    if not path.is_file():
        raise FileNotFoundError(path)
    surfaces: dict[SurfaceIdentity, Surface] = {}
    with uproot.open(path) as root_file:
        keys = _clean_keys(root_file)
        trigger_keys = [key for key in keys if _trigger_matches(key, trigger)]
        if not trigger_keys:
            top_levels = sorted({key.split("/", 1)[0] for key in keys if "/" in key})
            raise ValueError(
                f"{path}: trigger directory {trigger!r} absent; available top levels={top_levels}"
            )
        for key in trigger_keys:
            identity = _parse_identity(key, schema)
            if identity is None:
                continue
            if identity in surfaces:
                raise ValueError(
                    f"{path}: duplicate {schema} surface for {identity}: "
                    f"{surfaces[identity].source_key} and {key}"
                )
            surfaces[identity] = _surface_from_uproot(
                root_file[key],
                identity,
                key,
                "eiso_minus_cut" if schema == "v2" else "raw_eiso",
                allow_poisson_variance_fallback,
            )
    if not surfaces:
        raise ValueError(f"{path}: no {schema} surfaces found under trigger {trigger}")
    return surfaces


def validate_required_surface_coverage(
    loaded: Mapping[str, Mapping[SurfaceIdentity, Surface]], schema: str
) -> None:
    views = V2_VIEWS if schema == "v2" else ("legacy_mixed_internal_views",)
    requirements = {
        "data": "all",
        # THE-100 did not store the disjoint source-role category and remains
        # advisory only. Nominal-capable V2 packets must carry it explicitly.
        "inclusive": "truthBackground" if schema == "v2" else "all",
        "signal": "truthSignal",
    }
    failures: list[str] = []
    for sample, category in requirements.items():
        present = set(loaded[sample])
        expected = {
            SurfaceIdentity(centrality, view, category)
            for centrality in CENTRALITY_TOKENS
            for view in views
        }
        missing = sorted(expected - present)
        if missing:
            failures.append(f"{sample}/{category} missing {missing}")
    if failures:
        raise ValueError("required surface coverage incomplete: " + "; ".join(failures))


def _require_boundary(edges: np.ndarray, boundary: float, label: str) -> None:
    if not np.any(np.isclose(edges, boundary, rtol=0.0, atol=AXIS_TOLERANCE)):
        raise ValueError(
            f"{label}={boundary:g} is not a histogram edge; interpolation is intentionally forbidden"
        )


def _interval_mask(edges: np.ndarray, low: float, high: float) -> np.ndarray:
    """Return complete flow/regular bins contained in ``[low, high]``.

    Boundaries must be exact histogram edges.  In particular, the legacy
    0.02-wide score axis is scanned only at edge-aligned gaps.  The historical
    -0.20 to -0.03 band is retained as an external benchmark, not fabricated
    from a half-bin approximation.
    """

    if math.isfinite(low):
        _require_boundary(edges, low, "lower boundary")
    if math.isfinite(high):
        _require_boundary(edges, high, "upper boundary")
    lower = np.concatenate(([-np.inf], edges[:-1], [edges[-1]]))
    upper = np.concatenate(([edges[0]], edges[1:], [np.inf]))
    return (lower >= low - AXIS_TOLERANCE) & (upper <= high + AXIS_TOLERANCE)


def _sum_region(
    surface: Surface,
    x_interval: tuple[float, float],
    y_interval: tuple[float, float],
    z_regular_bins: Sequence[int],
) -> WeightedCount:
    x_mask = _interval_mask(surface.x_edges, *x_interval)
    y_mask = _interval_mask(surface.y_edges, *y_interval)
    z_indices = np.asarray([index + 1 for index in z_regular_bins], dtype=np.int64)
    if z_indices.size == 0:
        return WeightedCount(0.0, 0.0)
    selected_values = surface.values[np.ix_(x_mask, y_mask, z_indices)]
    selected_variances = surface.variances[np.ix_(x_mask, y_mask, z_indices)]
    return WeightedCount(float(np.sum(selected_values)), float(np.sum(selected_variances)))


def region_counts(
    surface: Surface,
    window: ScanWindow,
    z_regular_bins: Sequence[int],
    legacy_iso_cut: float,
    isolation_gap: float,
) -> RegionCounts:
    if surface.coordinate_kind == "eiso_minus_cut":
        iso_boundary = 0.0
    else:
        iso_boundary = legacy_iso_cut
    noniso_boundary = iso_boundary + isolation_gap
    tight = (0.0, math.inf)
    nontight = (window.lower_offset, window.upper_offset)
    isolated = (-math.inf, iso_boundary)
    nonisolated = (noniso_boundary, math.inf)
    return RegionCounts(
        A=_sum_region(surface, tight, isolated, z_regular_bins),
        B=_sum_region(surface, tight, nonisolated, z_regular_bins),
        C=_sum_region(surface, nontight, isolated, z_regular_bins),
        D=_sum_region(surface, nontight, nonisolated, z_regular_bins),
    )


def _safe_ratio(numerator: float, denominator: float) -> float:
    return numerator / denominator if denominator != 0.0 else math.nan


def _abcd_prediction(background: RegionCounts) -> tuple[float, float, float]:
    """Return predicted A, relative residual, and relative stat uncertainty."""

    b, c, d = background.B.sumw, background.C.sumw, background.D.sumw
    if background.A.sumw <= 0.0 or b <= 0.0 or c <= 0.0 or d <= 0.0:
        return math.nan, math.nan, math.nan
    prediction = b * c / d
    residual = (prediction - background.A.sumw) / background.A.sumw
    relative_variance = (
        background.B.sumw2 / (b * b)
        + background.C.sumw2 / (c * c)
        + background.D.sumw2 / (d * d)
    )
    return prediction, residual, math.sqrt(max(relative_variance, 0.0))


def solve_leakage_fixed_point(
    observed: RegionCounts,
    f_b: float,
    f_c: float,
    f_d: float,
    max_iterations: int = 500,
) -> tuple[float, bool, int]:
    """Solve the same damped leakage equation used by RecoilJets."""

    a, b, c, d = (
        observed.A.sumw,
        observed.B.sumw,
        observed.C.sumw,
        observed.D.sumw,
    )
    if a <= 0.0 or not all(math.isfinite(item) for item in (a, b, c, d, f_b, f_c, f_d)):
        return math.nan, False, 0
    estimate = min(max(a - b * c / d, 0.0), a) if d > 0.0 else a
    damping = 0.25
    tolerance = max(1.0e-11, 1.0e-11 * a)
    for iteration in range(1, max_iterations + 1):
        if f_d > 0.0:
            estimate = min(estimate, max(0.0, d / f_d * 0.999))
        denominator = d - f_d * estimate
        if denominator <= 0.0 or not math.isfinite(denominator):
            return estimate, False, iteration
        fixed = a - (b - f_b * estimate) * (c - f_c * estimate) / denominator
        if not math.isfinite(fixed):
            return estimate, False, iteration
        updated = min(max((1.0 - damping) * estimate + damping * fixed, 0.0), a)
        delta = abs(updated - estimate)
        estimate = updated
        if delta < tolerance:
            return estimate, True, iteration
        if iteration > 10 and delta > 0.5 * a and damping > 0.05:
            damping *= 0.5
    return estimate, False, max_iterations


def _add_counts(left: RegionCounts, right: RegionCounts, right_scale: float) -> RegionCounts:
    def add(name: str) -> WeightedCount:
        l_value = left.by_name(name)
        r_value = right.by_name(name)
        return WeightedCount(
            l_value.sumw + right_scale * r_value.sumw,
            l_value.sumw2 + right_scale * right_scale * r_value.sumw2,
        )

    return RegionCounts(**{name: add(name) for name in "ABCD"})


def evaluate_cell(
    data: RegionCounts,
    background: RegionCounts,
    truth_signal: RegionCounts,
    injected_purities: Sequence[float],
) -> dict[str, Any]:
    prediction, factorization_residual, prediction_stat = _abcd_prediction(background)
    signal_a = truth_signal.A.sumw
    leakage = {
        name: _safe_ratio(truth_signal.by_name(name).sumw, signal_a)
        for name in "BCD"
    }
    injections: list[dict[str, Any]] = []
    for purity in injected_purities:
        if not 0.0 < purity < 1.0 or signal_a <= 0.0 or background.A.sumw <= 0.0:
            injections.append(
                {
                    "injected_purity": purity,
                    "status": "invalid_inputs",
                    "converged": False,
                }
            )
            continue
        signal_scale = purity * background.A.sumw / ((1.0 - purity) * signal_a)
        pseudo = _add_counts(background, truth_signal, signal_scale)
        estimate, converged, iterations = solve_leakage_fixed_point(
            pseudo, leakage["B"], leakage["C"], leakage["D"]
        )
        true_signal_a = signal_scale * signal_a
        corrected_remainders = {
            name: pseudo.by_name(name).sumw - leakage[name] * estimate
            for name in "BCD"
        }
        nonnegative = all(value >= -1.0e-8 * max(pseudo.A.sumw, 1.0) for value in corrected_remainders.values())
        injections.append(
            {
                "injected_purity": purity,
                "signal_scale": signal_scale,
                "true_signal_A": true_signal_a,
                "estimated_signal_A": estimate,
                "estimated_purity": _safe_ratio(estimate, pseudo.A.sumw),
                "purity_residual": _safe_ratio(estimate, pseudo.A.sumw) - purity,
                "signal_yield_relative_residual": _safe_ratio(
                    estimate - true_signal_a, true_signal_a
                ),
                "converged": converged,
                "iterations": iterations,
                "corrected_remainders": corrected_remainders,
                "nonnegative_corrected_remainders": nonnegative,
            }
        )

    data_raw_purity = math.nan
    if data.A.sumw > 0.0 and data.D.sumw > 0.0:
        data_raw_purity = 1.0 - data.B.sumw * data.C.sumw / (data.A.sumw * data.D.sumw)
    return {
        "data_regions": {name: asdict(data.by_name(name)) | {"neff": data.by_name(name).neff} for name in "ABCD"},
        "background_regions": {
            name: asdict(background.by_name(name)) | {"neff": background.by_name(name).neff}
            for name in "ABCD"
        },
        "truth_signal_regions": {
            name: asdict(truth_signal.by_name(name)) | {"neff": truth_signal.by_name(name).neff}
            for name in "ABCD"
        },
        "data_raw_abcd_purity": data_raw_purity,
        "background_predicted_A": prediction,
        "background_factorization_relative_residual": factorization_residual,
        "background_abcd_relative_stat_uncertainty": prediction_stat,
        "background_factorization_double_ratio": _safe_ratio(
            background.A.sumw * background.D.sumw,
            background.B.sumw * background.C.sumw,
        ),
        "truth_signal_leakage_fractions_relative_to_A": leakage,
        "prompt_sideband_over_tight": _safe_ratio(truth_signal.C.sumw, signal_a),
        "injected_purity_closure": injections,
    }


def _assert_compatible(reference: Surface, candidate: Surface, context: str) -> None:
    for axis_name in ("x_edges", "y_edges", "z_edges"):
        left = getattr(reference, axis_name)
        right = getattr(candidate, axis_name)
        if not np.array_equal(left, right):
            raise ValueError(f"{context}: mismatched {axis_name}")
    if reference.coordinate_kind != candidate.coordinate_kind:
        raise ValueError(f"{context}: mismatched isolation coordinate")


def _category_map(
    surfaces: Mapping[SurfaceIdentity, Surface], category: str
) -> dict[tuple[str, str], Surface]:
    return {
        (identity.centrality, identity.view): surface
        for identity, surface in surfaces.items()
        if identity.category == category
    }


def evaluate_window(
    window: ScanWindow,
    data_surfaces: Mapping[SurfaceIdentity, Surface],
    signal_surfaces: Mapping[SurfaceIdentity, Surface],
    inclusive_surfaces: Mapping[SurfaceIdentity, Surface],
    injected_purities: Sequence[float],
    legacy_iso_cut: float,
    isolation_gap: float,
    gates: RankingGates,
) -> dict[str, Any]:
    data_map = _category_map(data_surfaces, "all")
    # This is the exact disjoint THE-111 PPG12 source-role background:
    # non-prompt candidates from embedded inclusive-jet sources. The
    # inclusive `all` surface also contains prompt cross-role candidates and
    # is diagnostic only, not a valid background closure template.
    background_map = _category_map(inclusive_surfaces, "truthBackground")
    if not background_map:
        # Historical THE-100 has only the mixed inclusive `all` category and
        # cannot nominate the nominal sideband. V2 coverage validation above
        # forbids this fallback for any nominal-capable packet.
        background_map = _category_map(inclusive_surfaces, "all")
    signal_map = _category_map(signal_surfaces, "truthSignal")
    common = sorted(set(data_map) & set(background_map) & set(signal_map))
    if not common:
        raise ValueError("no common data/background/truthSignal surfaces")
    missing = {
        "data": sorted((set(background_map) & set(signal_map)) - set(data_map)),
        "inclusive": sorted((set(data_map) & set(signal_map)) - set(background_map)),
        "truthSignal": sorted((set(data_map) & set(background_map)) - set(signal_map)),
    }
    if any(missing.values()):
        raise ValueError(f"surface identity mismatch: {missing}")

    cell_reports: list[dict[str, Any]] = []
    for centrality, view in common:
        data_surface = data_map[(centrality, view)]
        background_surface = background_map[(centrality, view)]
        signal_surface = signal_map[(centrality, view)]
        _assert_compatible(data_surface, background_surface, f"{centrality}/{view}/data-background")
        _assert_compatible(data_surface, signal_surface, f"{centrality}/{view}/data-signal")
        n_pt = data_surface.z_edges.size - 1
        scopes = [(f"pt_{data_surface.z_edges[i]:g}_{data_surface.z_edges[i + 1]:g}", [i]) for i in range(n_pt)]
        scopes.append(("pt_all", list(range(n_pt))))
        for pt_scope, z_bins in scopes:
            data_counts = region_counts(
                data_surface, window, z_bins, legacy_iso_cut, isolation_gap
            )
            background_counts = region_counts(
                background_surface, window, z_bins, legacy_iso_cut, isolation_gap
            )
            signal_counts = region_counts(
                signal_surface, window, z_bins, legacy_iso_cut, isolation_gap
            )
            report = evaluate_cell(
                data_counts, background_counts, signal_counts, injected_purities
            )
            report.update(
                {
                    "centrality": centrality,
                    "view": view,
                    "view_role": V2_VIEW_ROLES.get(view, "historical_mixed"),
                    "pt_scope": pt_scope,
                    "is_aggregate_pt_scope": pt_scope == "pt_all",
                }
            )
            cell_reports.append(report)

    # The gate is applied to each actual photon-pT analysis bin. Aggregate-pT
    # rows remain useful summaries but do not dilute a sparse or non-closing bin.
    gate_cells = [item for item in cell_reports if not item["is_aggregate_pt_scope"]]
    canonical_gate_cells = [
        item for item in gate_cells if item["view_role"] in ("canonical", "historical_mixed")
    ]
    failures: list[str] = []
    factorization = [
        abs(item["background_factorization_relative_residual"])
        for item in gate_cells
        if math.isfinite(item["background_factorization_relative_residual"])
    ]
    injected = [
        abs(injection.get("signal_yield_relative_residual", math.nan))
        for item in gate_cells
        for injection in item["injected_purity_closure"]
        if math.isfinite(injection.get("signal_yield_relative_residual", math.nan))
    ]
    min_neff_values = [
        item[f"{population}_regions"][region]["neff"]
        for item in gate_cells
        for population in ("data", "background")
        for region in ("C", "D")
    ]
    leakages = [
        item["prompt_sideband_over_tight"]
        for item in gate_cells
        if math.isfinite(item["prompt_sideband_over_tight"])
    ]
    canonical_factorization = [
        abs(item["background_factorization_relative_residual"])
        for item in canonical_gate_cells
        if math.isfinite(item["background_factorization_relative_residual"])
    ]
    canonical_injected = [
        abs(injection.get("signal_yield_relative_residual", math.nan))
        for item in canonical_gate_cells
        for injection in item["injected_purity_closure"]
        if math.isfinite(injection.get("signal_yield_relative_residual", math.nan))
    ]
    canonical_leakages = [
        item["prompt_sideband_over_tight"]
        for item in canonical_gate_cells
        if math.isfinite(item["prompt_sideband_over_tight"])
    ]
    canonical_neff_values = [
        item[f"{population}_regions"][region]["neff"]
        for item in canonical_gate_cells
        for population in ("data", "background")
        for region in ("C", "D")
    ]
    all_fixed_point_ok = all(
        injection.get("converged", False)
        and injection.get("nonnegative_corrected_remainders", False)
        for item in gate_cells
        for injection in item["injected_purity_closure"]
    )
    expected_injections = len(gate_cells) * len(injected_purities)
    if len(factorization) != len(gate_cells):
        failures.append("undefined background factorization in one or more pT cells")
    if len(injected) != expected_injections or not all_fixed_point_ok:
        failures.append("fixed-point leakage closure failed or produced a negative remainder")
    min_neff = min(min_neff_values, default=0.0)
    worst_factorization = max(factorization, default=math.inf)
    worst_injected = max(injected, default=math.inf)
    max_leakage = max(leakages, default=math.inf)
    if min_neff < gates.min_sideband_neff:
        failures.append(
            f"minimum C/D Neff {min_neff:.6g} < {gates.min_sideband_neff:.6g}"
        )
    if worst_factorization > gates.max_abs_background_factorization_residual:
        failures.append(
            "worst background factorization residual "
            f"{worst_factorization:.6g} > {gates.max_abs_background_factorization_residual:.6g}"
        )
    if worst_injected > gates.max_abs_injected_signal_yield_residual:
        failures.append(
            f"worst injected-purity yield residual {worst_injected:.6g} > "
            f"{gates.max_abs_injected_signal_yield_residual:.6g}"
        )
    if max_leakage > gates.max_prompt_sideband_over_tight:
        failures.append(
            f"maximum prompt C/A leakage {max_leakage:.6g} > "
            f"{gates.max_prompt_sideband_over_tight:.6g}"
        )
    return {
        "window": asdict(window)
        | {
            "lower_offset": window.lower_offset,
            "upper_offset": window.upper_offset,
            "key": window.key,
        },
        "passes_gates": not failures,
        "gate_failures": failures,
        "summary": {
            "n_gated_cells": len(gate_cells),
            "n_canonical_view_gated_cells": len(canonical_gate_cells),
            "n_aggregate_cells": len(cell_reports) - len(gate_cells),
            "minimum_data_or_background_C_D_neff": min_neff,
            "worst_abs_background_factorization_residual": worst_factorization,
            "median_abs_background_factorization_residual": float(np.median(factorization)) if factorization else math.inf,
            "worst_abs_injected_signal_yield_residual": worst_injected,
            "median_abs_injected_signal_yield_residual": float(np.median(injected)) if injected else math.inf,
            "max_prompt_sideband_over_tight": max_leakage,
            "canonical_view_minimum_data_or_background_C_D_neff": min(
                canonical_neff_values, default=0.0
            ),
            "canonical_view_worst_abs_background_factorization_residual": max(
                canonical_factorization, default=math.inf
            ),
            "canonical_view_worst_abs_injected_signal_yield_residual": max(
                canonical_injected, default=math.inf
            ),
            "canonical_view_max_prompt_sideband_over_tight": max(
                canonical_leakages, default=math.inf
            ),
            "all_fixed_points_converged_with_nonnegative_remainders": all_fixed_point_ok,
        },
        "cells": cell_reports,
    }


def _ranking_key(candidate: Mapping[str, Any]) -> tuple[Any, ...]:
    summary = candidate["summary"]
    window = candidate["window"]
    return (
        0 if candidate["passes_gates"] else 1,
        summary["canonical_view_worst_abs_injected_signal_yield_residual"],
        summary["canonical_view_worst_abs_background_factorization_residual"],
        summary["canonical_view_max_prompt_sideband_over_tight"],
        -summary["canonical_view_minimum_data_or_background_C_D_neff"],
        summary["worst_abs_injected_signal_yield_residual"],
        summary["worst_abs_background_factorization_residual"],
        summary["max_prompt_sideband_over_tight"],
        -summary["minimum_data_or_background_C_D_neff"],
        window["gap"],
        window["lower_width"],
    )


def rank_windows(
    windows: Sequence[ScanWindow],
    data_surfaces: Mapping[SurfaceIdentity, Surface],
    signal_surfaces: Mapping[SurfaceIdentity, Surface],
    inclusive_surfaces: Mapping[SurfaceIdentity, Surface],
    injected_purities: Sequence[float],
    legacy_iso_cut: float,
    isolation_gap: float,
    gates: RankingGates,
) -> list[dict[str, Any]]:
    candidates = [
        evaluate_window(
            window,
            data_surfaces,
            signal_surfaces,
            inclusive_surfaces,
            injected_purities,
            legacy_iso_cut,
            isolation_gap,
            gates,
        )
        for window in windows
    ]
    candidates.sort(key=_ranking_key)
    for index, candidate in enumerate(candidates, start=1):
        candidate["rank"] = index
    return candidates


def _parse_predeclared_subset(raw: str, allowed: Sequence[float], label: str) -> tuple[float, ...]:
    requested = tuple(float(item.strip()) for item in raw.split(",") if item.strip())
    if not requested:
        raise argparse.ArgumentTypeError(f"{label} cannot be empty")
    invalid = [
        item
        for item in requested
        if not any(math.isclose(item, allowed_item, rel_tol=0.0, abs_tol=1.0e-12) for allowed_item in allowed)
    ]
    if invalid:
        raise argparse.ArgumentTypeError(
            f"{label} values {invalid} are outside predeclared grid {list(allowed)}"
        )
    return tuple(sorted(set(requested)))


def _finite_json_value(value: Any) -> Any:
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {key: _finite_json_value(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_finite_json_value(item) for item in value]
    return value


def _write_csv(path: Path, candidates: Sequence[Mapping[str, Any]]) -> None:
    fields = (
        "rank",
        "passes_gates",
        "lower_offset",
        "upper_offset",
        "lower_width",
        "gap",
        "n_gated_cells",
        "minimum_data_or_background_C_D_neff",
        "canonical_view_minimum_data_or_background_C_D_neff",
        "worst_abs_background_factorization_residual",
        "canonical_view_worst_abs_background_factorization_residual",
        "median_abs_background_factorization_residual",
        "worst_abs_injected_signal_yield_residual",
        "canonical_view_worst_abs_injected_signal_yield_residual",
        "median_abs_injected_signal_yield_residual",
        "max_prompt_sideband_over_tight",
        "canonical_view_max_prompt_sideband_over_tight",
        "all_fixed_points_converged_with_nonnegative_remainders",
        "gate_failures",
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for candidate in candidates:
            row = {
                "rank": candidate["rank"],
                "passes_gates": candidate["passes_gates"],
                "lower_offset": candidate["window"]["lower_offset"],
                "upper_offset": candidate["window"]["upper_offset"],
                "lower_width": candidate["window"]["lower_width"],
                "gap": candidate["window"]["gap"],
                **candidate["summary"],
                "gate_failures": " | ".join(candidate["gate_failures"]),
            }
            writer.writerow({field: row.get(field, "") for field in fields})


def _schema_alias(value: str) -> str:
    aliases = {"the100": "legacy", "the112": "v2"}
    normalized = aliases.get(value, value)
    if normalized not in ("auto", "legacy", "v2"):
        raise argparse.ArgumentTypeError("schema must be auto, legacy/the100, or v2/the112")
    return normalized


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--signal", type=Path, required=True)
    parser.add_argument("--inclusive", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-csv", type=Path, required=True)
    parser.add_argument(
        "--schema",
        "--input-format",
        dest="schema",
        type=_schema_alias,
        default="auto",
        help="auto, legacy/the100, or v2/the112",
    )
    parser.add_argument(
        "--trigger",
        default=None,
        help=(
            "Shared trigger directory for all inputs (backward-compatible). "
            "Defaults to ALL unless lane-specific triggers are supplied."
        ),
    )
    parser.add_argument(
        "--data-trigger",
        default=None,
        help="Data trigger directory; overrides --trigger for the data input.",
    )
    parser.add_argument(
        "--simulation-trigger",
        default=None,
        help="Simulation trigger directory; overrides --trigger for signal and inclusive inputs.",
    )
    parser.add_argument("--legacy-iso-cut", type=float, default=4.0)
    parser.add_argument("--isolation-gap", type=float, default=0.0)
    parser.add_argument(
        "--gaps",
        default=None,
        help=(
            "Comma-separated subset of the frozen schema-specific gap grid. "
            "Defaults to exact legacy edges or the full V2 grid."
        ),
    )
    parser.add_argument(
        "--lower-widths",
        default=",".join(f"{item:g}" for item in PREDECLARED_LOWER_WIDTHS),
        help="Comma-separated subset of the frozen lower-width grid.",
    )
    parser.add_argument(
        "--injected-purities",
        default=",".join(f"{item:g}" for item in PREDECLARED_INJECTED_PURITIES),
        help="Comma-separated subset of the frozen injected-purity grid.",
    )
    parser.add_argument("--min-neff", type=float, default=25.0)
    parser.add_argument("--max-factorization-residual", type=float, default=0.50)
    parser.add_argument("--max-injected-purity-residual", type=float, default=0.50)
    parser.add_argument("--max-prompt-sideband-over-tight", type=float, default=1.00)
    parser.add_argument(
        "--allow-poisson-variance-fallback",
        action="store_true",
        help="Diagnostic-only fallback for unweighted synthetic/legacy files lacking Sumw2.",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.isolation_gap < 0.0:
        raise ValueError("isolation gap must be non-negative")
    if min(args.min_neff, args.max_factorization_residual, args.max_injected_purity_residual) < 0.0:
        raise ValueError("ranking gates must be non-negative")
    widths = _parse_predeclared_subset(
        args.lower_widths, PREDECLARED_LOWER_WIDTHS, "lower widths"
    )
    injected_purities = _parse_predeclared_subset(
        args.injected_purities,
        PREDECLARED_INJECTED_PURITIES,
        "injected purities",
    )
    paths = {"data": args.data, "signal": args.signal, "inclusive": args.inclusive}
    if uproot is None:
        raise RuntimeError(
            "uproot is unavailable; use /Users/patsfan753/Desktop/analysis/env/bin/python3"
        )
    roots = {label: uproot.open(path) for label, path in paths.items()}
    try:
        keys_by_input = {label: _clean_keys(root) for label, root in roots.items()}
    finally:
        for root in roots.values():
            root.close()
    schema = _detect_schema(keys_by_input, args.schema)
    allowed_gaps = LEGACY_EDGE_ALIGNED_GAPS if schema == "legacy" else PREDECLARED_GAPS
    gaps = (
        tuple(allowed_gaps)
        if args.gaps is None
        else _parse_predeclared_subset(args.gaps, allowed_gaps, f"{schema} gaps")
    )
    shared_trigger = args.trigger or "ALL"
    triggers = {
        "data": args.data_trigger or shared_trigger,
        "signal": args.simulation_trigger or shared_trigger,
        "inclusive": args.simulation_trigger or shared_trigger,
    }
    loaded = {
        label: load_surfaces(
            path,
            schema,
            triggers[label],
            args.allow_poisson_variance_fallback,
        )
        for label, path in paths.items()
    }
    validate_required_surface_coverage(loaded, schema)
    windows = [
        ScanWindow(gap=gap, lower_width=width)
        for gap in gaps
        for width in widths
        if width > gap
    ]
    gates = RankingGates(
        min_sideband_neff=args.min_neff,
        max_abs_background_factorization_residual=args.max_factorization_residual,
        max_abs_injected_signal_yield_residual=args.max_injected_purity_residual,
        max_prompt_sideband_over_tight=args.max_prompt_sideband_over_tight,
    )
    candidates = rank_windows(
        windows,
        loaded["data"],
        loaded["signal"],
        loaded["inclusive"],
        injected_purities,
        args.legacy_iso_cut,
        args.isolation_gap,
        gates,
    )
    for candidate in candidates:
        candidate["eligible_for_nominal_selection"] = (
            schema == "v2" and candidate["passes_gates"]
        )
    interpretation = (
        "historical_prior_only: THE-100 mixed fixed/sliding and cone views in one raw-Eiso "
        "surface; do not promote a nominal sideband from this packet alone"
        if schema == "legacy"
        else "corrected_canary_surface: candidate-level selection evidence only; no recoil or xJ read"
    )
    advisory_best = candidates[0]["window"] if candidates else None
    nominal_best = next(
        (
            candidate["window"]
            for candidate in candidates
            if candidate["eligible_for_nominal_selection"]
        ),
        None,
    )
    report = {
        "schema": SCHEMA_VERSION,
        "status": "historical_prior_ranked" if schema == "legacy" else "ranked",
        "input_schema": schema,
        "interpretation": interpretation,
        "blinding_contract": {
            "objects_read": "candidate-level score/isolation/photon-pT TH3 surfaces only",
            "real_data_recoil_histograms_read": False,
            "real_data_xjgamma_histograms_read": False,
            "unfolded_results_read": False,
            "selection_frozen_before_recoil_result_inspection": True,
        },
        "inputs": {
            label: {"path": str(path), "sha256": _file_sha256(path)}
            for label, path in paths.items()
        },
        "scan_contract": {
            "trigger": (
                shared_trigger if len(set(triggers.values())) == 1 else None
            ),
            "triggers": triggers,
            "gaps": list(gaps),
            "lower_widths": list(widths),
            "injected_purities": list(injected_purities),
            "legacy_iso_cut": args.legacy_iso_cut if schema == "legacy" else None,
            "legacy_score_boundary_policy": (
                "exact 0.02-axis-edge integration only; the historical -0.20 to -0.03 "
                "band is not approximated from half bins"
                if schema == "legacy"
                else None
            ),
            "external_legacy_benchmark_not_scanned": (
                {
                    "lower_offset": -0.20,
                    "upper_offset": -0.03,
                    "reason": "upper edge is not representable on the THE-100 0.02 score grid",
                }
                if schema == "legacy"
                else None
            ),
            "view_roles": (
                V2_VIEW_ROLES
                if schema == "v2"
                else {"legacy_mixed_internal_views": "historical_prior_only"}
            ),
            "isolation_gap": args.isolation_gap,
            "gates": asdict(gates),
            "ranking_order": [
                "passes gates",
                "R=0.4 sliding canonical-view injected-purity closure",
                "R=0.4 sliding canonical-view background factorization",
                "R=0.4 sliding canonical-view prompt leakage and statistics",
                "worst R=0.3 robustness-view behavior",
                "gap and lower width deterministic tie break",
            ],
        },
        "candidate_count": len(candidates),
        "passing_candidate_count": sum(bool(item["passes_gates"]) for item in candidates),
        "canonical_eligible_candidate_count": sum(
            bool(item["eligible_for_nominal_selection"]) for item in candidates
        ),
        "best_candidate": nominal_best,
        "advisory_best_score_band_prior": advisory_best,
        "candidates": candidates,
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(
        json.dumps(_finite_json_value(report), indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    _write_csv(args.output_csv, candidates)
    print(
        json.dumps(
            {
                "status": report["status"],
                "input_schema": schema,
                "candidate_count": len(candidates),
                "passing_candidate_count": report["passing_candidate_count"],
                "best_candidate": report["best_candidate"],
                "advisory_best_score_band_prior": report[
                    "advisory_best_score_band_prior"
                ],
                "output_json": str(args.output_json),
                "output_csv": str(args.output_csv),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        end="",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
