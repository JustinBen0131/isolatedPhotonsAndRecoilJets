#!/usr/bin/env python3
"""Pure numerical correction and unfolding helpers for the THE-243 H70 overlay."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Mapping, Sequence

import numpy as np


PT_EDGES = np.asarray([15.0, 20.0, 25.0, 35.0])
CATEGORIES = ("A", "B", "C", "D")
DATA_SCHEMA = "THE243Schema10InclusiveXJGammaThresholdFanoutV1"
RESPONSE_SCHEMA = "THE243Schema10InclusiveResponseThresholdFanoutV1"


@dataclass(frozen=True)
class ABCDSolution:
    signal_a: float
    purity: float
    background_a_over_c: float
    residual: float
    method: str


@dataclass(frozen=True)
class ResponseBundle:
    system: str
    xj_edges: np.ndarray
    photon_response: np.ndarray
    photon_misses: np.ndarray
    photon_truth: np.ndarray
    photon_reco: np.ndarray
    photon_reco_sumw2: np.ndarray
    photon_matched_reco: np.ndarray
    photon_boundary_fakes_reco: np.ndarray
    combinatoric_normalization_reco: np.ndarray
    combinatoric_normalization_reco_sumw2: np.ndarray
    xj_response: np.ndarray
    xj_misses: np.ndarray
    xj_truth: np.ndarray
    xj_reco: np.ndarray
    xj_fakes_reco: np.ndarray
    xj_detector_fakes_reco: np.ndarray
    xj_detector_fakes_reco_sumw2: np.ndarray
    xj_boundary_fakes_reco: np.ndarray
    combinatoric_reco: np.ndarray
    combinatoric_reco_sumw2: np.ndarray
    leakage: np.ndarray
    leakage_sumw2: np.ndarray
    sample_bindings: Mapping[str, str]


def solve_leakage_abcd(counts: Sequence[float], leakage: Sequence[float]) -> ABCDSolution:
    """Solve factorized ABCD backgrounds with signal leakage relative to A."""
    n = np.asarray(counts, dtype=float)
    l = np.asarray(leakage, dtype=float)
    if n.shape != (4,) or l.shape != (4,) or np.any(~np.isfinite(n)) or np.any(n < 0):
        raise ValueError("ABCD counts/leakage must be four finite nonnegative values")
    if l[0] <= 0:
        raise ValueError("ABCD signal leakage requires positive A normalization")
    l = l / l[0]
    a = l[3] - l[1] * l[2]
    b = -n[0] * l[3] - n[3] + n[1] * l[2] + n[2] * l[1]
    c = n[0] * n[3] - n[1] * n[2]
    upper = n[0]
    for value, ratio in zip(n[1:], l[1:]):
        if ratio > 0:
            upper = min(upper, value / ratio)
    naive = float(np.clip(n[0] - n[1] * n[2] / n[3], 0.0, upper)) if n[3] > 0 else upper
    roots: list[float] = []
    if abs(a) < 1.0e-14:
        if abs(b) > 1.0e-14:
            roots = [-c / b]
    else:
        discriminant = b * b - 4.0 * a * c
        if discriminant >= 0:
            root = math.sqrt(discriminant)
            roots = [(-b - root) / (2.0 * a), (-b + root) / (2.0 * a)]
    physical = [float(value) for value in roots if -1.0e-9 <= value <= upper + 1.0e-9]
    method = "analytic_physical_root"
    if physical:
        signal = min(physical, key=lambda value: abs(value - naive))
    else:
        # Match the established ROOT correction contract exactly: when the
        # leakage-aware fixed-point equation has no physical solution, retain
        # the clipped raw ABCD estimate.  A residual-minimizing boundary point
        # is not an ABCD solution and previously injected a spurious zero
        # 25--35 GeV Au+Au photon bin into the joint unfolding.
        signal = naive
        method = "raw_abcd_fallback_no_physical_leakage_root"
    backgrounds = n - l * signal
    bkg_c = backgrounds[2]
    alpha = backgrounds[0] / bkg_c if bkg_c > 0 else 0.0
    residual = abs(backgrounds[0] * backgrounds[3] - backgrounds[1] * backgrounds[2])
    scale = max(1.0, abs(backgrounds[0] * backgrounds[3]), abs(backgrounds[1] * backgrounds[2]))
    return ABCDSolution(
        signal_a=max(0.0, signal),
        purity=max(0.0, min(1.0, signal / n[0])) if n[0] > 0 else 0.0,
        background_a_over_c=max(0.0, alpha),
        residual=float(residual / scale),
        method=method,
    )


def _edge_indices(source_edges: Sequence[float], target_edges: Sequence[float]) -> list[int]:
    source = np.asarray(source_edges, dtype=float)
    result: list[int] = []
    for target in target_edges:
        matches = np.flatnonzero(np.isclose(source, target, rtol=0.0, atol=1.0e-9))
        if len(matches) != 1:
            raise ValueError(f"target edge {target} is not on source grid")
        result.append(int(matches[0]))
    if any(right <= left for left, right in zip(result[:-1], result[1:])):
        raise ValueError("target edges are not strictly increasing")
    return result


def rebin_1d(values: Sequence[float], source_edges: Sequence[float], target_edges: Sequence[float]) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    indices = _edge_indices(source_edges, target_edges)
    if array.shape != (len(source_edges) - 1,):
        raise ValueError("1D accumulator shape differs")
    return np.asarray([array[left:right].sum() for left, right in zip(indices[:-1], indices[1:])])


def rebin_2d(values: Sequence[Sequence[float]], x_edges: Sequence[float], y_edges: Sequence[float],
             target_x: Sequence[float], target_y: Sequence[float]) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    ix = _edge_indices(x_edges, target_x)
    iy = _edge_indices(y_edges, target_y)
    if array.shape != (len(x_edges) - 1, len(y_edges) - 1):
        raise ValueError("2D accumulator shape differs")
    return np.asarray([[array[x0:x1, y0:y1].sum() for y0, y1 in zip(iy[:-1], iy[1:])]
                       for x0, x1 in zip(ix[:-1], ix[1:])])


def _signal_samples(system: str) -> tuple[str, ...]:
    """Return every disjoint production-owned photon source for a system.

    Source ownership is already enforced during THE-121/122 production using
    the generator-level producer photon (pp: Photon5 <14, Photon10 [14,22),
    Photon20 >=22; AuAu: Photon12 [12,20), Photon20 >=20).  The response must
    therefore add all owned samples.  Selecting another sample from reco or
    truth pT here would discard migrations across those source boundaries.
    """
    if system == "pp":
        return ("pp_photon5", "pp_photon10", "pp_photon20")
    if system == "auau":
        return ("auau_photon12", "auau_photon20")
    raise ValueError(system)


def build_response_bundle(samples: Mapping[str, Mapping[str, Any]], system: str,
                          threshold: int, xj_edges: Sequence[float]) -> ResponseBundle:
    """Add production-exclusive photon samples into one weighted response.

    pp ``RJEventV1.event_weight`` already contains the producer-side photon
    slice factor ``sigma(slice) / sigma(Photon20)``.  Multiplying those rows by
    ``sigma(slice) / Ngen`` again would square the relative slice scaling.
    Therefore pp samples receive the common ``sigma(Photon20) / Ngen`` scale;
    Au+Au *photon-signal* rows carry unit event weights and retain
    ``sigma(slice) / Ngen``.  Embedded inclusive-jet samples are not response
    sources here and their ``sigma_eff/Npass`` weights must never enter this
    photon-signal path.
    """
    tag = f"pt{threshold}"
    required = _signal_samples(system)
    missing = sorted(set(required).difference(samples))
    if missing:
        raise ValueError(f"missing signal response samples: {missing}")
    for sample_id in required:
        if samples[sample_id].get("schema") != RESPONSE_SCHEMA:
            raise ValueError(
                f"{sample_id} is not an inclusive-pair response product; "
                "leading-jet response products cannot enter additive combinatoric subtraction"
            )
    first = samples[required[0]]["reductions"][tag]
    source_xj = first["xj_edges"]
    nx = len(xj_edges) - 1
    npt = len(PT_EDGES) - 1
    photon_response = np.zeros((npt, npt))
    photon_misses = np.zeros(npt)
    photon_truth = np.zeros(npt)
    photon_reco = np.zeros(npt)
    photon_reco_sumw2 = np.zeros(npt)
    photon_raw_fakes = np.zeros(npt)
    photon_boundary_fakes = np.zeros(npt)
    xj_response = np.zeros((npt * nx, npt * nx))
    xj_misses = np.zeros(npt * nx)
    xj_truth = np.zeros((npt, nx))
    xj_reco = np.zeros((npt, nx))
    xj_fakes = np.zeros((npt, nx))
    xj_detector_fakes = np.zeros((npt, nx))
    xj_detector_fakes_sumw2 = np.zeros((npt, nx))
    xj_boundary_fakes = np.zeros((npt, nx))
    combinatoric = np.zeros((npt, nx))
    combinatoric_sumw2 = np.zeros((npt, nx))
    leakage = np.zeros((4, npt))
    leakage_sumw2 = np.zeros((4, npt))

    pp_reference_cross_section = None
    if system == "pp":
        reference = samples["pp_photon20"]
        reference_contract = reference.get("sample_contract")
        if not isinstance(reference_contract, Mapping):
            raise ValueError("pp Photon20 response sample contract missing")
        pp_reference_cross_section = float(reference_contract.get("cross_section_pb", math.nan))
        if not math.isfinite(pp_reference_cross_section) or pp_reference_cross_section <= 0:
            raise ValueError("pp Photon20 reference cross section is invalid")

    def payload(sample_id: str) -> tuple[Mapping[str, Any], float]:
        item = samples[sample_id]
        generated_events = int(item.get("generated_events", 0))
        contract = item.get("sample_contract")
        if generated_events <= 0 or not isinstance(contract, Mapping):
            raise ValueError(f"{sample_id} response normalization contract missing")
        if (
            contract.get("system") != system
            or contract.get("source_class") != "photon_signal"
        ):
            raise ValueError(
                f"{sample_id} is not a {system} photon-signal response source; "
                "embedded inclusive sigma_eff/Npass weights cannot enter this correction path"
            )
        cross_section = float(contract.get("cross_section_pb", math.nan))
        recorded_scale = float(item.get("cross_section_weight_pb_per_event", math.nan))
        expected_scale = cross_section / generated_events
        if (
            not math.isfinite(cross_section)
            or cross_section <= 0
            or not math.isfinite(recorded_scale)
            or not math.isclose(recorded_scale, expected_scale, rel_tol=1.0e-12, abs_tol=0.0)
        ):
            raise ValueError(f"{sample_id} response cross-section normalization differs")
        scale = (
            float(pp_reference_cross_section) / generated_events
            if system == "pp"
            else recorded_scale
        )
        return item["reductions"][tag], scale

    # Every sample has already passed its mutually exclusive generator-level
    # production ownership gate.  Sum all samples for every reco/truth bin so
    # migrations across an ownership edge remain represented exactly once.
    for sample_id in required:
        reduced, weight = payload(sample_id)
        truth_source_edges = reduced["truth_ptgamma_edges"]
        reco_source_edges = reduced["reco_ptgamma_edges"]
        truth_indices = _edge_indices(truth_source_edges, PT_EDGES)
        reco_indices = _edge_indices(reco_source_edges, PT_EDGES)
        truth_inside_bins = set(range(truth_indices[0], truth_indices[-1]))
        reco_inside_bins = set(range(reco_indices[0], reco_indices[-1]))
        truth_outside_bins = [
            index for index in range(len(truth_source_edges) - 1)
            if index not in truth_inside_bins
        ]
        reco_outside_bins = [
            index for index in range(len(reco_source_edges) - 1)
            if index not in reco_inside_bins
        ]
        response = np.asarray(reduced["photon"]["response_truth_x_reco"], dtype=float) * weight
        for truth_pt, (ti0, ti1) in enumerate(zip(truth_indices[:-1], truth_indices[1:])):
            photon_response[truth_pt] += [
                response[ti0:ti1, r0:r1].sum()
                for r0, r1 in zip(reco_indices[:-1], reco_indices[1:])
            ]
            photon_misses[truth_pt] += (
                np.asarray(reduced["photon"]["misses"], dtype=float)[ti0:ti1].sum() * weight
                + response[np.ix_(range(ti0, ti1), reco_outside_bins)].sum()
            )
            photon_truth[truth_pt] += (
                np.asarray(reduced["photon"]["truth"], dtype=float)[ti0:ti1].sum() * weight
            )
        for reco_pt, (r0, r1) in enumerate(zip(reco_indices[:-1], reco_indices[1:])):
            photon_boundary_fakes[reco_pt] += response[
                np.ix_(truth_outside_bins, range(r0, r1))
            ].sum()
        truth_x_indices = _edge_indices(source_xj, xj_edges)
        source_nx = len(source_xj) - 1
        matrix = np.asarray(reduced["xj"]["response_truth_global_x_reco_global"], dtype=float) * weight
        misses = np.asarray(reduced["xj"]["misses"], dtype=float) * weight
        truth = np.asarray(reduced["xj"]["truth"], dtype=float) * weight
        for truth_pt, (ti0, ti1) in enumerate(zip(truth_indices[:-1], truth_indices[1:])):
            for tx, (x0, x1) in enumerate(zip(truth_x_indices[:-1], truth_x_indices[1:])):
                target_truth = truth_pt * nx + tx
                xj_misses[target_truth] += misses[ti0:ti1, x0:x1].sum()
                xj_truth[truth_pt, tx] += truth[ti0:ti1, x0:x1].sum()
                truth_rows = [
                    p * source_nx + x
                    for p in range(ti0, ti1)
                    for x in range(x0, x1)
                ]
                reco_outside_columns = [
                    p * source_nx + x
                    for p in reco_outside_bins
                    for x in range(source_nx)
                ]
                xj_misses[target_truth] += matrix[
                    np.ix_(truth_rows, reco_outside_columns)
                ].sum()
                for reco_pt, (r0, r1) in enumerate(zip(reco_indices[:-1], reco_indices[1:])):
                    for rx, (y0, y1) in enumerate(zip(truth_x_indices[:-1], truth_x_indices[1:])):
                        target_reco = reco_pt * nx + rx
                        rows = [p * source_nx + x for p in range(ti0, ti1) for x in range(x0, x1)]
                        cols = [p * source_nx + x for p in range(r0, r1) for x in range(y0, y1)]
                        xj_response[target_truth, target_reco] += matrix[np.ix_(rows, cols)].sum()
        truth_outside_rows = [
            p * source_nx + x
            for p in truth_outside_bins
            for x in range(source_nx)
        ]
        for reco_pt, (r0, r1) in enumerate(zip(reco_indices[:-1], reco_indices[1:])):
            for rx, (y0, y1) in enumerate(zip(truth_x_indices[:-1], truth_x_indices[1:])):
                reco_columns = [
                    p * source_nx + x
                    for p in range(r0, r1)
                    for x in range(y0, y1)
                ]
                xj_boundary_fakes[reco_pt, rx] += matrix[
                    np.ix_(truth_outside_rows, reco_columns)
                ].sum()

        photon_reco += rebin_1d(
            reduced["photon"]["reco"], reco_source_edges, PT_EDGES
        ) * weight
        photon_reco_sumw2 += rebin_1d(
            reduced["photon"]["reco_sumw2"], reco_source_edges, PT_EDGES
        ) * weight * weight
        photon_raw_fakes += rebin_1d(
            reduced["photon"]["fakes"], reco_source_edges, PT_EDGES
        ) * weight
        rebinned_reco = rebin_2d(reduced["xj"]["reco"], reduced["reco_ptgamma_edges"], source_xj,
                                 PT_EDGES, xj_edges)
        rebinned_comb = rebin_2d(reduced["xj"]["combinatoric_reco"], reduced["reco_ptgamma_edges"], source_xj,
                                 PT_EDGES, xj_edges)
        rebinned_comb_sumw2 = rebin_2d(
            reduced["xj"]["combinatoric_reco_sumw2"],
            reduced["reco_ptgamma_edges"],
            source_xj,
            PT_EDGES,
            xj_edges,
        )
        rebinned_fakes = rebin_2d(reduced["xj"]["fakes"], reduced["reco_ptgamma_edges"], source_xj,
                                  PT_EDGES, xj_edges)
        rebinned_detector_fakes = rebin_2d(
            reduced["xj"]["detector_fakes_reco"],
            reduced["reco_ptgamma_edges"],
            source_xj,
            PT_EDGES,
            xj_edges,
        )
        rebinned_detector_fakes_sumw2 = rebin_2d(
            reduced["xj"]["detector_fakes_reco_sumw2"],
            reduced["reco_ptgamma_edges"],
            source_xj,
            PT_EDGES,
            xj_edges,
        )
        xj_reco += rebinned_reco * weight
        combinatoric += rebinned_comb * weight
        combinatoric_sumw2 += rebinned_comb_sumw2 * weight * weight
        xj_fakes += rebinned_fakes * weight
        xj_detector_fakes += rebinned_detector_fakes * weight
        xj_detector_fakes_sumw2 += rebinned_detector_fakes_sumw2 * weight * weight
        for category_index, category in enumerate(CATEGORIES):
            leakage[category_index] += rebin_1d(
                reduced["signal_leakage"]["truth_matched_abcd"][category],
                reco_source_edges,
                PT_EDGES,
            ) * weight
            leakage_sumw2[category_index] += rebin_1d(
                reduced["signal_leakage"]["truth_matched_abcd_sumw2"][category],
                reco_source_edges,
                PT_EDGES,
            ) * weight * weight

    photon_reco_partition = (
        photon_response.sum(axis=0) + photon_boundary_fakes + photon_raw_fakes
    )
    photon_truth_partition = photon_response.sum(axis=1) + photon_misses
    xj_reco_partition = (
        xj_response.sum(axis=0).reshape(npt, nx) + xj_fakes + xj_boundary_fakes
    )
    xj_truth_partition = xj_response.sum(axis=1).reshape(npt, nx) + xj_misses.reshape(npt, nx)
    for label, observed, partitioned in (
        ("photon reco", photon_reco, photon_reco_partition),
        ("photon truth", photon_truth, photon_truth_partition),
        ("xJ reco", xj_reco, xj_reco_partition),
        ("xJ truth", xj_truth, xj_truth_partition),
    ):
        if not np.allclose(observed, partitioned, rtol=1.0e-9, atol=1.0e-9):
            delta = float(np.max(np.abs(observed - partitioned)))
            raise ValueError(f"{system} weighted {label} response partition differs (max |delta|={delta:g})")

    return ResponseBundle(
        system=system, xj_edges=np.asarray(xj_edges, dtype=float),
        photon_response=photon_response, photon_misses=photon_misses,
        photon_truth=photon_truth, photon_reco=photon_reco,
        photon_reco_sumw2=photon_reco_sumw2,
        # The unmatched-recoil template is defined for every truth-matched
        # Region-A photon in the reconstructed pT bin.  Do not normalize it
        # with the response-matrix column sum: that excludes otherwise valid
        # matched photons whose truth pT migrates across the 15/35 GeV
        # analysis boundaries and would over-subtract the measured recoil.
        photon_matched_reco=leakage[0].copy(),
        photon_boundary_fakes_reco=photon_boundary_fakes,
        # The source-side RooUnfold contract conditions the SIM reconstructed
        # photon spectrum on a matched truth-signal photon.  The combinatoric
        # template is defined on that same matched-photon population, so its
        # photon-yield denominator must be the event-leading Region-A matched
        # population as well.  Using every reconstructed SIM photon dilutes
        # K(xJ) with fake-photon tags that ABCD has already removed from data.
        combinatoric_normalization_reco=leakage[0].copy(),
        combinatoric_normalization_reco_sumw2=leakage_sumw2[0].copy(),
        xj_response=xj_response, xj_misses=xj_misses,
        xj_truth=xj_truth, xj_reco=xj_reco, xj_fakes_reco=xj_fakes,
        xj_detector_fakes_reco=xj_detector_fakes,
        xj_detector_fakes_reco_sumw2=xj_detector_fakes_sumw2,
        xj_boundary_fakes_reco=xj_boundary_fakes,
        combinatoric_reco=combinatoric,
        combinatoric_reco_sumw2=combinatoric_sumw2,
        leakage=leakage,
        leakage_sumw2=leakage_sumw2,
        sample_bindings={
            sample_id: "production-exclusive generator-level ownership; weighted additive contribution"
            for sample_id in required
        },
    )


def aggregate_data(payload: Mapping[str, Any], threshold: int, xj_edges: Sequence[float]) -> dict[str, np.ndarray]:
    if payload.get("schema") != DATA_SCHEMA:
        raise ValueError(
            "data input is not an inclusive-pair reduction; leading-jet spectra "
            "are incompatible with additive combinatoric subtraction"
        )
    reduced = payload["reductions"][f"pt{threshold}"]
    source_pt = reduced["ptgamma_fine_edges"]
    source_xj = reduced["xj_fine_edges"]
    result: dict[str, np.ndarray] = {}
    for field in ("abcd_event_leading_counts_ptgamma", "abcd_event_leading_sumw2_ptgamma"):
        result[field] = np.asarray([
            rebin_1d(reduced[field][category], source_pt, PT_EDGES) for category in CATEGORIES
        ])
    for field in ("inclusive_recoil_spectra_ptgamma_xj", "inclusive_recoil_sumw2_ptgamma_xj"):
        result[field] = np.asarray([
            rebin_2d(reduced[field][category], source_pt, source_xj, PT_EDGES, xj_edges)
            for category in ("A", "C")
        ])
    return result


def correct_measured(
    data: Mapping[str, np.ndarray],
    response: ResponseBundle,
    subtract_unmatched_recoil: bool,
    *,
    combinatoric_reco: np.ndarray | None = None,
    combinatoric_normalization_reco: np.ndarray | None = None,
    purity_strategy: str = "per_pt",
) -> dict[str, Any]:
    counts = data["abcd_event_leading_counts_ptgamma"]
    spectra = data["inclusive_recoil_spectra_ptgamma_xj"]
    spectra_sumw2 = data["inclusive_recoil_sumw2_ptgamma_xj"]
    counts_sumw2 = data["abcd_event_leading_sumw2_ptgamma"]
    corrected = np.zeros_like(spectra[0])
    variance = np.zeros_like(spectra_sumw2[0])
    photons = np.zeros(len(PT_EDGES) - 1)
    solutions: list[dict[str, Any]] = []
    comb_subtracted = np.zeros_like(corrected)
    comb_variance = np.zeros_like(corrected)
    comb_template = (
        response.combinatoric_reco
        if combinatoric_reco is None
        else np.asarray(combinatoric_reco, dtype=float)
    )
    comb_normalization = (
        response.combinatoric_normalization_reco
        if combinatoric_normalization_reco is None
        else np.asarray(combinatoric_normalization_reco, dtype=float)
    )
    if comb_template.shape != corrected.shape:
        raise ValueError("combinatoric template shape differs from measured recoil shape")
    if comb_normalization.shape != photons.shape:
        raise ValueError("combinatoric photon-normalization shape differs")
    if purity_strategy not in {"per_pt", "integrated_15_35_transfer"}:
        raise ValueError(f"unknown purity strategy: {purity_strategy}")

    integrated_solution: ABCDSolution | None = None
    if purity_strategy == "integrated_15_35_transfer":
        # The reported observable is integrated over 15--35 GeV.  When the
        # upper-pT sidebands are too sparse to support three independent ABCD
        # solves, determine one leakage-aware background transfer from the
        # exact integrated A/B/C/D population.  Apply that transfer inside
        # each reconstructed photon-pT bin so the response still sees the
        # measured pT spectrum and its boundary migrations.  This is not a
        # pooled recoil shape: Region A and C spectra remain bin-local.
        integrated_solution = solve_leakage_abcd(
            counts.sum(axis=1),
            response.leakage.sum(axis=1),
        )
        if integrated_solution.residual >= 0.02:
            raise ValueError(
                "integrated 15--35 GeV ABCD transfer does not factorize: "
                f"residual={integrated_solution.residual}"
            )
    for pt in range(len(PT_EDGES) - 1):
        local_solution = solve_leakage_abcd(counts[:, pt], response.leakage[:, pt])
        solution = integrated_solution or local_solution
        transfer = solution.background_a_over_c
        leakage_a = response.leakage[0, pt]
        leakage_c = response.leakage[2, pt]
        if leakage_a <= 0:
            raise ValueError(f"missing Region-A prompt leakage normalization in photon-pT bin {pt}")
        relative_c_leakage = leakage_c / leakage_a
        leakage_denominator = 1.0 - transfer * relative_c_leakage
        if not math.isfinite(leakage_denominator) or leakage_denominator <= 1.0e-6:
            raise ValueError(
                f"nonphysical Region-C leakage denominator in photon-pT bin {pt}: "
                f"1-k*fC={leakage_denominator}"
            )
        if purity_strategy == "integrated_15_35_transfer":
            photons[pt] = (
                counts[0, pt] - transfer * counts[2, pt]
            ) / leakage_denominator
            if not math.isfinite(photons[pt]) or photons[pt] <= 0:
                raise ValueError(
                    "integrated 15--35 GeV ABCD transfer gives a nonpositive "
                    f"signal in photon-pT bin {pt}: {photons[pt]}"
                )
        else:
            photons[pt] = local_solution.signal_a
        numerator = spectra[0, pt] - transfer * spectra[1, pt]
        corrected[pt] = numerator / leakage_denominator
        variance[pt] = (
            spectra_sumw2[0, pt]
            + transfer ** 2 * spectra_sumw2[1, pt]
        ) / (leakage_denominator ** 2)
        if subtract_unmatched_recoil and comb_normalization[pt] > 0:
            # ABCD has already removed fake-photon recoil.  Scale only the
            # response component with a matched prompt photon and an unmatched
            # reconstructed recoil, thereby avoiding double subtraction.
            comb_subtracted[pt] = (
                comb_template[pt]
                * photons[pt] / comb_normalization[pt]
            )
            corrected[pt] -= comb_subtracted[pt]
            # Propagate finite embedded-template and photon-normalization
            # statistics.  Inclusive recoil bins are jet-pair yields and are
            # not a partition of the photon denominator, so no Bernoulli
            # subset covariance is invented here.
            template = comb_template[pt]
            template_var = response.combinatoric_reco_sumw2[pt]
            denominator = comb_normalization[pt]
            denominator_var = response.combinatoric_normalization_reco_sumw2[pt]
            scale = photons[pt] / denominator
            derivative_denominator = -photons[pt] * template / (denominator * denominator)
            comb_variance[pt] = np.clip(
                scale * scale * template_var
                + derivative_denominator * derivative_denominator * denominator_var,
                0.0,
                None,
            )
            variance[pt] += comb_variance[pt]
        solutions.append({
            "pt_low": float(PT_EDGES[pt]), "pt_high": float(PT_EDGES[pt + 1]),
            "signal_a": float(photons[pt]),
            "purity": float(photons[pt] / counts[0, pt]) if counts[0, pt] > 0 else 0.0,
            "background_a_over_c": transfer,
            "relative_region_c_signal_leakage": float(relative_c_leakage),
            "region_c_signal_leakage_restoration": float(1.0 / leakage_denominator),
            "factorization_residual": solution.residual,
            "local_factorization_residual": local_solution.residual,
            "method": (
                "integrated_15_35_leakage_transfer"
                if integrated_solution is not None else local_solution.method
            ),
            "combinatoric_sim_yield_per_photon": (
                float(comb_template[pt].sum() / comb_normalization[pt])
                if subtract_unmatched_recoil and comb_normalization[pt] > 0 else 0.0
            ),
            "combinatoric_scaled_recoil": float(comb_subtracted[pt].sum()),
            "purity_corrected_recoil_before_combinatoric": float(
                (corrected[pt] + comb_subtracted[pt]).sum()
            ),
        })
    return {
        "corrected": corrected,
        "variance": variance,
        "photons": photons,
        "solutions": solutions,
        "combinatoric_subtracted": comb_subtracted,
        "combinatoric_variance": comb_variance,
        "negative_input_sum": float(np.abs(corrected[corrected < 0]).sum()),
        "purity_strategy": purity_strategy,
    }


def correct_integrated_measured(data: Mapping[str, np.ndarray], response: ResponseBundle) -> dict[str, Any]:
    """Purity-correct the integrated 15--35 GeV recoil spectrum.

    The headline observable is integrated over the declared photon range.  A
    single leakage-aware ABCD solve avoids imposing an unstable three-bin
    purity decomposition on the sparse high-pT AuAu sidebands.  Region C is a
    per-tag fake-photon recoil template, so its normalization is the inferred
    Region-A background divided by the *observed* Region-C tag count.  The
    detector unmatched-recoil component is not subtracted here as an absolute
    MC rate; it is applied later as a reconstructed-bin matching fraction.
    """
    counts = np.asarray(data["abcd_event_leading_counts_ptgamma"], dtype=float).sum(axis=1)
    spectra = np.asarray(data["inclusive_recoil_spectra_ptgamma_xj"], dtype=float).sum(axis=1)
    spectra_sumw2 = np.asarray(data["inclusive_recoil_sumw2_ptgamma_xj"], dtype=float).sum(axis=1)
    leakage = np.asarray(response.leakage, dtype=float).sum(axis=1)
    solution = solve_leakage_abcd(counts, leakage)
    background_a = max(0.0, float(counts[0] - solution.signal_a))
    region_c_tags = float(counts[2])
    alpha = background_a / region_c_tags if region_c_tags > 0 else 0.0
    corrected = spectra[0] - alpha * spectra[1]
    variance = spectra_sumw2[0] + alpha * alpha * spectra_sumw2[1]
    negative = float(np.abs(corrected[corrected < 0]).sum())
    return {
        "corrected": corrected,
        "variance": variance,
        "photons": float(solution.signal_a),
        "solutions": [{
            "pt_low": float(PT_EDGES[0]),
            "pt_high": float(PT_EDGES[-1]),
            "signal_a": solution.signal_a,
            "purity": solution.purity,
            "background_a_over_observed_c": alpha,
            "factorization_residual": solution.residual,
            "method": solution.method,
        }],
        "region_a_recoil": spectra[0],
        "region_c_recoil": spectra[1],
        "region_c_scale": alpha,
        "negative_input_sum": negative,
    }


def iterative_bayes(measured: Sequence[float], response_truth_x_reco: np.ndarray,
                    misses: Sequence[float], iterations: int,
                    prior: Sequence[float] | None = None, *,
                    fakes: Sequence[float] | None = None) -> tuple[np.ndarray, np.ndarray]:
    """D'Agostini unfolding with RooUnfold-compatible measured fakes.

    RooUnfoldBayes represents reconstructed entries without a reported truth
    partner as one additional cause whose detector distribution is ``fakes``.
    The fake cause participates in the Bayesian normalization and is dropped
    from the returned truth spectrum, while remaining in the refolded result.
    """
    measured_array = np.asarray(measured, dtype=float)
    matrix = np.asarray(response_truth_x_reco, dtype=float)
    miss = np.asarray(misses, dtype=float)
    if matrix.ndim != 2 or matrix.shape[1] != measured_array.size or matrix.shape[0] != miss.size:
        raise ValueError("unfolding dimensions differ")
    if iterations < 1:
        raise ValueError("at least one Bayesian iteration is required")
    totals = matrix.sum(axis=1) + miss
    conditional = np.divide(matrix, totals[:, None], out=np.zeros_like(matrix), where=totals[:, None] > 0)
    efficiency = conditional.sum(axis=1)
    if prior is None:
        state_counts = totals.copy()
    else:
        state_counts = np.asarray(prior, dtype=float).copy()
        if state_counts.shape != totals.shape:
            raise ValueError("unfolding prior dimensions differ")
    n_truth = matrix.shape[0]
    fake_array = np.zeros(matrix.shape[1], dtype=float) if fakes is None else np.asarray(fakes, dtype=float)
    if fake_array.shape != (matrix.shape[1],) or np.any(~np.isfinite(fake_array)) or np.any(fake_array < 0):
        raise ValueError("unfolding fakes must match reconstructed bins and be finite nonnegative")
    fake_total = float(fake_array.sum())
    if fake_total > 0:
        conditional = np.vstack((conditional, fake_array / fake_total))
        efficiency = np.append(efficiency, 1.0)
        state_counts = np.append(state_counts, fake_total)
    state = state_counts
    if state.sum() <= 0:
        raise ValueError("unfolding prior is empty")
    state /= state.sum()
    unfolded = np.zeros(conditional.shape[0])
    # Mirror RooUnfoldBayes::unfold exactly: signed, background-subtracted
    # measured bins enter nbarC = sum_j M_ij * nEst_j without clipping.  A
    # non-positive folded prior bin receives zero inverse weight on a later
    # iteration, just as RooUnfold's ``Uj > 0 ? 1/Uj : 0`` guard does.
    used = measured_array
    for _ in range(iterations):
        denominator = conditional.T @ state
        posterior = np.divide(conditional * state[:, None], denominator[None, :],
                              out=np.zeros_like(conditional), where=denominator[None, :] > 0)
        unfolded = np.divide(posterior @ used, efficiency,
                             out=np.zeros(conditional.shape[0]), where=efficiency > 0)
        unfolded_total = float(unfolded.sum())
        if not math.isfinite(unfolded_total) or abs(unfolded_total) <= 1.0e-15:
            raise ValueError("Bayesian unfolding produced an empty signed estimate")
        state = unfolded / unfolded_total
    refolded = conditional.T @ unfolded
    return unfolded[:n_truth], refolded


def chi2_ndf(observed: Sequence[float], expected: Sequence[float], variance: Sequence[float]) -> float:
    obs = np.asarray(observed, dtype=float)
    exp = np.asarray(expected, dtype=float)
    var = np.asarray(variance, dtype=float)
    mask = np.isfinite(obs) & np.isfinite(exp) & (var > 0)
    return float(np.sum((obs[mask] - exp[mask]) ** 2 / var[mask]) / max(1, int(mask.sum())))


def draw_joint_data_toy(data: Mapping[str, np.ndarray], rng: np.random.Generator) -> dict[str, np.ndarray]:
    """Fluctuate inclusive-pair sufficient statistics without a false partition.

    Multiple recoil jets may accompany one photon, so the xJ bins do not
    partition the event-leading ABCD counts.  Until an event-bootstrap
    covariance product is supplied, fluctuate the recorded count and pair
    sumw2 terms independently.  This is conservative and explicit; it never
    forces the inclusive jet yield below the photon count.
    """
    counts = np.asarray(data["abcd_event_leading_counts_ptgamma"], dtype=float)
    counts_var = np.asarray(data["abcd_event_leading_sumw2_ptgamma"], dtype=float)
    spectra = np.asarray(data["inclusive_recoil_spectra_ptgamma_xj"], dtype=float)
    spectra_var = np.asarray(data["inclusive_recoil_sumw2_ptgamma_xj"], dtype=float)
    toy_counts = np.clip(
        rng.normal(counts, np.sqrt(np.clip(counts_var, 0.0, None))),
        0.0,
        None,
    )
    toy_spectra = np.clip(
        rng.normal(spectra, np.sqrt(np.clip(spectra_var, 0.0, None))),
        0.0,
        None,
    )
    return {
        **data,
        "abcd_event_leading_counts_ptgamma": toy_counts,
        "inclusive_recoil_spectra_ptgamma_xj": toy_spectra,
    }


def draw_combinatoric_template_toy(
    response: ResponseBundle,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    """Fluctuate inclusive combinatoric pairs and photon normalization.

    The template is a per-photon inclusive jet yield and can exceed one when
    integrated over xJ.  It is therefore never treated as a Bernoulli
    partition of the reconstructed-photon population.
    """
    template = np.asarray(response.combinatoric_reco, dtype=float)
    template_var = np.asarray(response.combinatoric_reco_sumw2, dtype=float)
    denominator = np.asarray(response.combinatoric_normalization_reco, dtype=float)
    denominator_var = np.asarray(response.combinatoric_normalization_reco_sumw2, dtype=float)
    toy_template = np.clip(
        rng.normal(template, np.sqrt(np.clip(template_var, 0.0, None))),
        0.0,
        None,
    )
    toy_denominator = np.clip(
        rng.normal(denominator, np.sqrt(np.clip(denominator_var, 0.0, None))),
        0.0,
        None,
    )
    return toy_template, toy_denominator


def run_chain(data_payload: Mapping[str, Any], samples: Mapping[str, Mapping[str, Any]],
              system: str, threshold: int, xj_edges: Sequence[float], iterations: int,
              *, toys: int = 400, seed: int = 243,
              purity_strategy: str = "per_pt",
              subtract_combinatoric: bool = True) -> dict[str, Any]:
    response = build_response_bundle(samples, system, threshold, xj_edges)
    data = aggregate_data(data_payload, threshold, xj_edges)
    # The optional response-consistent combinatoric correction is an additive,
    # photon-yield-scaled template in each reconstructed photon-pT bin; it is
    # not a global efficiency factor.  THE-254 pp explicitly disables it so
    # the presentation interface does not inherit the distinct AuAu K term.
    correction = correct_measured(
        data,
        response,
        subtract_unmatched_recoil=subtract_combinatoric,
        purity_strategy=purity_strategy,
    )
    measured = correction["corrected"].reshape(-1)
    variance = correction["variance"].reshape(-1)
    unfolded, refolded = iterative_bayes(
        measured,
        response.xj_response,
        response.xj_misses,
        iterations,
        response.xj_truth.reshape(-1),
        fakes=(
            response.xj_detector_fakes_reco
            + response.xj_boundary_fakes_reco
        ).reshape(-1),
    )
    unfolded_photons, refolded_photons = iterative_bayes(
        correction["photons"], response.photon_response, response.photon_misses,
        iterations, response.photon_truth,
        fakes=response.photon_boundary_fakes_reco,
    )
    nx = len(xj_edges) - 1
    projected = unfolded.reshape(len(PT_EDGES) - 1, nx).sum(axis=0)
    photon_total = unfolded_photons.sum()
    widths = np.diff(np.asarray(xj_edges, dtype=float))
    density = np.divide(projected, photon_total * widths, out=np.zeros(nx), where=(photon_total * widths) > 0)

    rng = np.random.default_rng(seed)
    toy_density: list[np.ndarray] = []
    toy_measured: list[np.ndarray] = []
    central_purity_methods = tuple(row["method"] for row in correction["solutions"])
    for _ in range(toys):
        toy_data = draw_joint_data_toy(data, rng)
        if subtract_combinatoric:
            toy_comb, toy_comb_normalization = draw_combinatoric_template_toy(response, rng)
        else:
            toy_comb = response.combinatoric_reco
            toy_comb_normalization = response.combinatoric_normalization_reco
        try:
            # Gaussian sufficient-statistic fluctuations can leave the
            # physical leakage-aware ABCD parameter space.  Reject those
            # draws and count them against the explicit toy-success gate;
            # never convert an invalid denominator into a physics estimate.
            toy = correct_measured(
                toy_data,
                response,
                subtract_unmatched_recoil=subtract_combinatoric,
                combinatoric_reco=toy_comb,
                combinatoric_normalization_reco=toy_comb_normalization,
                purity_strategy=purity_strategy,
            )
            # A fallback or alternate ABCD branch is a different estimator,
            # not a statistical fluctuation of the selected central method.
            # Mixing branches creates pathological non-Gaussian tails that
            # previously dominated the displayed pp uncertainty.  Count such
            # toys against the success gate instead of hiding them inside an
            # RMS error bar.
            if tuple(row["method"] for row in toy["solutions"]) != central_purity_methods:
                continue
            toy_unfolded, _ = iterative_bayes(
                toy["corrected"].reshape(-1),
                response.xj_response,
                response.xj_misses,
                iterations,
                response.xj_truth.reshape(-1),
                fakes=(
                    response.xj_detector_fakes_reco
                    + response.xj_boundary_fakes_reco
                ).reshape(-1),
            )
            toy_photons, _ = iterative_bayes(toy["photons"], response.photon_response,
                                             response.photon_misses, iterations, response.photon_truth,
                                             fakes=response.photon_boundary_fakes_reco)
        except ValueError:
            continue
        total = toy_photons.sum()
        if total > 0:
            toy_measured.append(toy["corrected"].reshape(-1))
            toy_density.append(toy_unfolded.reshape(len(PT_EDGES) - 1, nx).sum(axis=0) / (total * widths))
    toy_array = np.asarray(toy_density)
    if len(toy_array) > 1:
        lower = np.quantile(toy_array, 0.16, axis=0)
        upper = np.quantile(toy_array, 0.84, axis=0)
        density_error = 0.5 * (upper - lower)
        # Build a robust PSD covariance for the same central 68% estimator.
        # Per-bin winsorization prevents a handful of signed-unfolding tail
        # toys from defining the covariance, while the final rescaling makes
        # its diagonal exactly equal to the reported percentile errors.
        winsorized = np.clip(toy_array, lower, upper)
        winsorized_covariance = np.atleast_2d(
            np.cov(winsorized, rowvar=False, ddof=1)
        )
        winsorized_scale = np.sqrt(
            np.clip(np.diag(winsorized_covariance), 0.0, None)
        )
        winsorized_correlation = np.divide(
            winsorized_covariance,
            winsorized_scale[:, None] * winsorized_scale[None, :],
            out=np.eye(nx),
            where=(winsorized_scale[:, None] * winsorized_scale[None, :]) > 0,
        )
        density_covariance = (
            winsorized_correlation * density_error[:, None] * density_error[None, :]
        )
    else:
        lower = np.full(nx, np.nan)
        upper = np.full(nx, np.nan)
        density_error = np.full(nx, np.nan)
        density_covariance = np.full((nx, nx), np.nan)
    covariance_scale = np.sqrt(np.clip(np.diag(density_covariance), 0.0, None))
    density_correlation = np.divide(
        density_covariance,
        covariance_scale[:, None] * covariance_scale[None, :],
        out=np.zeros_like(density_covariance),
        where=(covariance_scale[:, None] * covariance_scale[None, :]) > 0,
    )
    measured_toy_array = np.asarray(toy_measured)
    measured_variance = (
        measured_toy_array.var(axis=0, ddof=1)
        if len(measured_toy_array) > 1
        else np.clip(variance, 1.0e-12, None)
    )

    matched_mc_reco = response.xj_response.sum(axis=0)
    residual_reco_fakes = (
        response.xj_detector_fakes_reco + response.xj_boundary_fakes_reco
    ).reshape(-1)
    closure_mc_reco = matched_mc_reco + residual_reco_fakes
    mc_unfolded, mc_refolded = iterative_bayes(
        closure_mc_reco,
        response.xj_response,
        response.xj_misses,
        iterations,
        response.xj_truth.reshape(-1),
        fakes=residual_reco_fakes,
    )
    closure_truth = response.xj_truth.reshape(-1)
    closure_variance = np.clip(closure_truth + mc_unfolded, 1.0e-12, None)
    xj_totals = response.xj_response.sum(axis=1) + response.xj_misses
    xj_efficiency = np.divide(
        response.xj_response.sum(axis=1),
        xj_totals,
        out=np.zeros_like(xj_totals),
        where=xj_totals > 0,
    )
    photon_totals = response.photon_response.sum(axis=1) + response.photon_misses
    photon_efficiency = np.divide(
        response.photon_response.sum(axis=1),
        photon_totals,
        out=np.zeros_like(photon_totals),
        where=photon_totals > 0,
    )
    xj_supported = closure_truth > 0
    xj_edges_array = np.asarray(xj_edges, dtype=float)
    xj_centers = 0.5 * (xj_edges_array[:-1] + xj_edges_array[1:])
    reported_xj = (xj_edges_array[:-1] >= 0.30) & (xj_centers <= 1.80)
    reported_xj_global = np.tile(reported_xj, len(PT_EDGES) - 1)
    photon_supported = response.photon_truth > 0
    return {
        "system": system, "threshold_gev": threshold, "iterations": iterations,
        "purity_strategy": purity_strategy,
        "xj_edges": list(map(float, xj_edges)), "ptgamma_edges": PT_EDGES.tolist(),
        "density": density.tolist(), "density_error": density_error.tolist(),
        "density_error_low": (density - lower).tolist(),
        "density_error_high": (upper - density).tolist(),
        "density_error_method": "central_68_percentile_branch_matched_toys",
        "density_stat_covariance": density_covariance.tolist(),
        "density_stat_correlation": density_correlation.tolist(),
        "unfolded_counts": projected.tolist(), "unfolded_photons": float(photon_total),
        "purity": correction["solutions"],
        "combinatoric_subtracted": correction["combinatoric_subtracted"].tolist(),
        "diagnostics": {
            "region_a_recoil": data["inclusive_recoil_spectra_ptgamma_xj"][0].sum(axis=0).tolist(),
            "region_c_recoil": data["inclusive_recoil_spectra_ptgamma_xj"][1].sum(axis=0).tolist(),
            "purity_corrected_recoil": (
                correction["corrected"] + correction["combinatoric_subtracted"]
            ).sum(axis=0).tolist(),
            "scaled_combinatoric_recoil": correction["combinatoric_subtracted"].sum(axis=0).tolist(),
            "post_combinatoric_recoil": correction["corrected"].sum(axis=0).tolist(),
            "matched_response_reco": response.xj_response.sum(axis=0).reshape(
                len(PT_EDGES) - 1, nx
            ).sum(axis=0).tolist(),
            "detector_fake_reco": response.xj_detector_fakes_reco.sum(axis=0).tolist(),
            "boundary_migration_fake_reco": response.xj_boundary_fakes_reco.sum(axis=0).tolist(),
            "refolded_reco": refolded.reshape(len(PT_EDGES) - 1, nx).sum(axis=0).tolist(),
        },
        "negative_input_sum": correction["negative_input_sum"],
        "positive_input_sum": float(np.clip(measured, 0.0, None).sum()),
        "negative_input_fraction": float(
            correction["negative_input_sum"]
            / max(correction["negative_input_sum"] + np.clip(measured, 0.0, None).sum(), 1.0e-12)
        ),
        # Empty reconstructed bins carry no defined statistical variance.
        # chi2_ndf deliberately masks them; assigning an artificial 1e-12
        # variance turns harmless response leakage into an O(1e10) failure.
        "refold_chi2_ndf": chi2_ndf(
            measured,
            refolded,
            measured_variance,
        ),
        "photon_refold_chi2_ndf": chi2_ndf(correction["photons"], refolded_photons,
                                             np.clip(correction["photons"], 1.0, None)),
        "mc_closure_chi2_ndf": chi2_ndf(closure_truth, mc_unfolded, closure_variance),
        "mc_refold_chi2_ndf": chi2_ndf(closure_mc_reco, mc_refolded, closure_mc_reco),
        "toy_successes": len(toy_density), "toy_requested": toys,
        "response_observability": {
            "xj_supported_truth_bins": int(xj_supported.sum()),
            "xj_zero_efficiency_supported_bins": int(np.sum(xj_supported & (xj_efficiency <= 0))),
            "xj_zero_efficiency_reported_bins": int(np.sum(
                xj_supported & reported_xj_global & (xj_efficiency <= 0)
            )),
            "xj_min_supported_efficiency": float(np.min(xj_efficiency[xj_supported])) if np.any(xj_supported) else 0.0,
            "photon_supported_truth_bins": int(photon_supported.sum()),
            "photon_zero_efficiency_supported_bins": int(np.sum(photon_supported & (photon_efficiency <= 0))),
            "photon_min_supported_efficiency": float(np.min(photon_efficiency[photon_supported])) if np.any(photon_supported) else 0.0,
        },
        "response_internal_closure": {
            "xj_reco_max_abs_residual": float(np.max(np.abs(
                response.xj_reco
                - (
                    response.xj_response.sum(axis=0).reshape(len(PT_EDGES) - 1, nx)
                    + response.xj_fakes_reco
                    + response.xj_boundary_fakes_reco
                )
            ))),
            "xj_truth_max_abs_residual": float(np.max(np.abs(
                response.xj_truth
                - (
                    response.xj_response.sum(axis=1).reshape(len(PT_EDGES) - 1, nx)
                    + response.xj_misses.reshape(len(PT_EDGES) - 1, nx)
                )
            ))),
        },
        "sample_bindings": dict(response.sample_bindings),
        "method": {
            "purity": "leakage-aware event-leading ABCD solved independently in three photon-pT bins",
            "reco_fakes": (
                "ABCD removes fake-photon recoil before the distinct matched-photon combinatoric template"
            ),
            "combinatoric": (
                "simulated matched-photon inclusive recoil jets without a nearby truth jet are photon-yield scaled per reconstructed pTgamma bin and subtracted additively before unfolding"
                if subtract_combinatoric else
                "disabled for the pp direct-histogram presentation interface; pp does not acquire an AuAu K-subtraction term"
            ),
            "unfolding": (
                "iterative Bayesian joint pTgamma-xJ response with explicit misses and a RooUnfold-compatible "
                "residual detector/order plus reporting-window migration fake cause; fake-photon and unrelated-recoil components are "
                "excluded here because they are removed by ABCD and Au+Au combinatoric subtraction"
            ),
            "negative_bins": (
                "preserved as signed measured inputs exactly as RooUnfoldBayes; only a non-positive "
                "folded-prior denominator receives RooUnfold's zero inverse-weight guard"
            ),
            "data_statistical_covariance": (
                "event-leading ABCD counts and inclusive recoil-pair sumw2 are fluctuated independently; "
                "the compact product does not yet carry event-bootstrap cross-bin covariance"
            ),
        },
    }


def run_integrated_chain(
    data_payload: Mapping[str, Any],
    samples: Mapping[str, Mapping[str, Any]],
    system: str,
    threshold: int,
    xj_edges: Sequence[float],
    iterations: int,
    *,
    toys: int = 400,
    seed: int = 243,
) -> dict[str, Any]:
    """Deprecated global-matching implementation retained only for audit history.

    Multiplying the measured recoil by a global matched/(matched+unmatched)
    fraction is not the accepted Au+Au combinatoric correction.  It suppresses
    the entire data spectrum and can reproduce the exact failure mode seen in
    the rejected overlay.  Keep this symbol fail-closed so old callers cannot
    silently reintroduce that bias; ``run_chain`` is the only accepted entry.
    """
    raise RuntimeError(
        "run_integrated_chain is retired: use run_chain with additive, "
        "photon-yield-scaled Au+Au combinatoric subtraction"
    )
