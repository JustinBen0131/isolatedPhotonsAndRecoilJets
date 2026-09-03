#!/usr/bin/env python3
"""Canonical pp/schema10 normalization with explicit, non-interchangeable channels.

Generator-slice stitching is a generator-spectrum closure diagnostic. It uses
raw event counts and the accepted per-sample cross section. The stored
``RJEventV1.event_weight`` is instead the downstream analysis weight; for pp it
already contains the generator-slice factor relative to a family reference as
well as non-slice weights. Keeping both channels explicitly named prevents a
downstream-weighted spectrum from being mistaken for a stitching diagnostic.

Au+Au embedded inclusive-jet samples are deliberately rejected here.  Their
accepted source stitch is ``sigma_eff / Npass`` in the matching half-open
truth-jet ownership window and is provided only by
``auau_embedded_inclusive_schema10_weighting``.
"""

from __future__ import annotations

import copy
import math
from typing import Any, Mapping


SCHEMA = "CanonicalSchema10StitchingWeightingV3"
THE291_CONTRACT_ROLE = "pp_only_rejects_auau_inclusive"
GENERATOR_STITCHING_CHANNEL = "generator_stitching"
ANALYSIS_WEIGHTED_CHANNEL = "downstream_analysis_weighted"
PP_REFERENCE_CROSS_SECTION_PB = {
    "photon_signal": 130.4461,
    "inclusive_background": 7.3113,
}


class Schema10WeightingError(ValueError):
    """Raised when a sample cannot satisfy the canonical weighting contract."""


def _finite(value: Any, label: str, *, positive: bool = False) -> float:
    if isinstance(value, bool):
        raise Schema10WeightingError(f"{label} is not numeric")
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise Schema10WeightingError(f"{label} is not numeric") from exc
    if not math.isfinite(result) or (positive and result <= 0.0):
        raise Schema10WeightingError(f"{label} is invalid")
    return result


def normalization_contract() -> dict[str, Any]:
    return {
        "schema": SCHEMA,
        "stitching_plot_channel": GENERATOR_STITCHING_CHANNEL,
        "channels": {
            GENERATOR_STITCHING_CHANNEL: {
                "purpose": "generator-slice spectrum closure and stitching plots",
                "source_moment": "raw_counts",
                "density_formula": (
                    "raw_counts * sample_cross_section_pb / generated_events / bin_width_gev"
                ),
                "variance_formula": (
                    "raw_counts * (sample_cross_section_pb / generated_events / bin_width_gev)^2"
                ),
            },
            ANALYSIS_WEIGHTED_CHANNEL: {
                "purpose": "downstream analysis-weighted distributions; not a stitching plot channel",
                "source_moments": ["sumw", "sumw2"],
                "pp_density_formula": (
                    "sumw * family_reference_cross_section_pb / generated_events / bin_width_gev"
                ),
                "auau_density_formula": (
                    "not_applicable_to_embedded_inclusive; consume the canonical "
                    "producer_then_sigma_eff_over_Npass_then_centrality artifact"
                ),
                "pp_reference_cross_section_pb": dict(PP_REFERENCE_CROSS_SECTION_PB),
                "pp_stored_event_weight_includes": [
                    "generator_slice_cross_section_relative_to_family_reference",
                    "truth_vertex_weight",
                    "interaction_mix_weight",
                    "period_luminosity_weight",
                ],
            },
        },
        "prohibitions": [
            "do not use sumw or sumw2 for a generator-stitching plot",
            "do not use raw_counts for a downstream analysis-weighted distribution",
            "do not multiply pp sumw by the per-sample cross section",
            "do not normalize Au+Au embedded inclusive jets with cross section / generated events",
            "do not emit ambiguous weighted_density fields",
        ],
    }


def normalize_stitching_sample(
    merged: Mapping[str, Any],
    catalog_contract: Mapping[str, Any],
    *,
    bin_width_gev: float,
) -> dict[str, Any]:
    """Return a copy with explicit generator and downstream density channels."""

    sample_id = str(merged.get("sample_id", ""))
    system = str(merged.get("system", ""))
    source_class = str(merged.get("source_class", ""))
    if not sample_id or system not in {"pp", "auau"}:
        raise Schema10WeightingError("sample identity or system is invalid")
    if catalog_contract.get("sample_id") not in (None, sample_id):
        raise Schema10WeightingError(f"{sample_id} catalog sample id differs")
    for field in ("system", "source_class"):
        if catalog_contract.get(field) != merged.get(field):
            raise Schema10WeightingError(f"{sample_id} catalog {field} differs")

    if system == "auau" and source_class == "inclusive_background":
        raise Schema10WeightingError(
            f"{sample_id} is an Au+Au embedded inclusive sample; cross section / generated "
            "events is forbidden. Consume the hash-bound sigma_eff/Npass source-stitch "
            "artifact through auau_embedded_inclusive_schema10_weighting."
        )
    if system == "auau" and source_class != "photon_signal":
        raise Schema10WeightingError(
            f"{sample_id} has unsupported Au+Au source class {source_class!r}"
        )

    generated_events = catalog_contract.get("generated_events")
    if type(generated_events) is not int or generated_events <= 0:
        raise Schema10WeightingError("catalog.generated_events is invalid")
    if merged.get("events_total") != generated_events:
        raise Schema10WeightingError(
            f"{sample_id} generated-event denominator differs: "
            f"{merged.get('events_total')} != {generated_events}"
        )
    bin_width = _finite(bin_width_gev, "bin_width_gev", positive=True)
    cross_section = _finite(
        catalog_contract.get("cross_section_pb"), "catalog.cross_section_pb", positive=True
    )
    catalog_weight = _finite(
        catalog_contract.get("cross_section_weight_pb_per_event"),
        "catalog.cross_section_weight_pb_per_event",
        positive=True,
    )
    expected_catalog_weight = cross_section / generated_events
    if not math.isclose(catalog_weight, expected_catalog_weight, rel_tol=2e-12, abs_tol=1e-18):
        raise Schema10WeightingError(f"{sample_id} cross-section weight differs")

    raw_counts = merged.get("raw_counts")
    sumw = merged.get("sumw")
    sumw2 = merged.get("sumw2")
    if (
        not isinstance(raw_counts, list)
        or not isinstance(sumw, list)
        or not isinstance(sumw2, list)
        or len(raw_counts) != len(sumw)
        or len(sumw) != len(sumw2)
    ):
        raise Schema10WeightingError(f"{sample_id} stored moments differ")
    if any(type(value) is not int or value < 0 for value in raw_counts):
        raise Schema10WeightingError(f"{sample_id} raw counts are invalid")

    if system == "pp":
        try:
            analysis_reference_cross_section = PP_REFERENCE_CROSS_SECTION_PB[source_class]
        except KeyError as exc:
            raise Schema10WeightingError(
                f"{sample_id} has unsupported pp source class {source_class!r}"
            ) from exc
        if sample_id == "pp_photon20" and not math.isclose(
            cross_section, analysis_reference_cross_section, rel_tol=0.0, abs_tol=1e-9
        ):
            raise Schema10WeightingError("pp Photon20 reference cross section differs")
        analysis_mode = "pp_stored_analysis_weight_times_common_family_reference"
        stored_weight_includes_slice = True
    else:
        analysis_reference_cross_section = cross_section
        analysis_mode = "auau_stored_analysis_weight_times_sample_cross_section"
        stored_weight_includes_slice = False

    generator_density_scale = catalog_weight / bin_width
    analysis_integrated_scale = analysis_reference_cross_section / generated_events
    analysis_density_scale = analysis_integrated_scale / bin_width
    generator_density = [count * generator_density_scale for count in raw_counts]
    generator_variance = [count * generator_density_scale**2 for count in raw_counts]
    analysis_density = [
        _finite(value, f"{sample_id}.sumw[{index}]") * analysis_density_scale
        for index, value in enumerate(sumw)
    ]
    analysis_variance = [
        _finite(value, f"{sample_id}.sumw2[{index}]") * analysis_density_scale**2
        for index, value in enumerate(sumw2)
    ]
    if any(value < 0.0 for value in analysis_variance):
        raise Schema10WeightingError(f"{sample_id} has negative second moments")

    output = copy.deepcopy(dict(merged))
    for legacy_field in (
        "weighted_density_pb_per_gev",
        "weighted_density_sumw2_pb2_per_gev2",
        "normalization_mode",
        "normalization_reference_cross_section_pb",
        "normalization_integrated_scale_pb_per_stored_weight",
        "normalization_density_scale_pb_per_gev_per_stored_weight",
    ):
        output.pop(legacy_field, None)
    output.update(
        {
            "cross_section_pb": cross_section,
            "generated_events": generated_events,
            "cross_section_weight_pb_per_event": catalog_weight,
            "production_campaign_tag": catalog_contract["production_campaign_tag"],
            "normalization_schema": SCHEMA,
            "generator_stitching_channel": {
                "source_moment": "raw_counts",
                "sample_cross_section_pb": cross_section,
                "generated_events": generated_events,
                "integrated_scale_pb_per_raw_event": catalog_weight,
                "density_scale_pb_per_gev_per_raw_event": generator_density_scale,
            },
            "generator_stitching_density_pb_per_gev": generator_density,
            "generator_stitching_density_sumw2_pb2_per_gev2": generator_variance,
            "analysis_weighted_channel": {
                "source_moments": ["sumw", "sumw2"],
                "mode": analysis_mode,
                "stored_event_weight_includes_generator_slice": stored_weight_includes_slice,
                "reference_cross_section_pb": analysis_reference_cross_section,
                "integrated_scale_pb_per_stored_weight": analysis_integrated_scale,
                "density_scale_pb_per_gev_per_stored_weight": analysis_density_scale,
            },
            "analysis_weighted_density_pb_per_gev": analysis_density,
            "analysis_weighted_density_sumw2_pb2_per_gev2": analysis_variance,
        }
    )
    return output
