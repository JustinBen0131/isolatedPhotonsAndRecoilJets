#!/usr/bin/env python3
"""Fail-closed authority gate for nominal ThesisAnalysis schema10 plots."""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from scripts.data_prep.recoiljets.pp_schema10_weighting import (
    SCHEMA as WEIGHTING_SCHEMA,
    normalize_stitching_sample,
    normalization_contract,
)


REPO = Path(__file__).resolve().parents[2]
CORRECTED_ASSEMBLY_SCHEMA = "THE248Schema10StitchingAssemblyV4"
EXPECTED_PP_SAMPLES = {
    "pp_photon5",
    "pp_photon10",
    "pp_photon20",
    "pp_jet8",
    "pp_jet12",
    "pp_jet20",
    "pp_jet30",
    "pp_jet40",
}


class CanonicalBaselineError(RuntimeError):
    """Raised when nominal plot provenance is absent, stale, or ambiguous."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _load_object(path: Path) -> dict[str, Any]:
    if path.is_symlink() or not path.is_file():
        raise CanonicalBaselineError(f"required canonical artifact is not a regular file: {path}")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise CanonicalBaselineError(f"canonical artifact is not valid JSON: {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise CanonicalBaselineError(f"canonical artifact must contain a JSON object: {path}")
    return value


def _resolve_repo_path(value: Any, label: str) -> Path:
    if not isinstance(value, str) or not value:
        raise CanonicalBaselineError(f"{label}.path is missing")
    candidate = Path(value).expanduser()
    if not candidate.is_absolute():
        candidate = REPO / candidate
    resolved = candidate.resolve()
    try:
        resolved.relative_to(REPO.resolve())
    except ValueError as exc:
        raise CanonicalBaselineError(f"{label}.path escapes the ThesisAnalysis repository") from exc
    return resolved


@dataclass(frozen=True)
class BaselineEvidence:
    authority_path: Path
    authority_sha256: str
    authority: dict[str, Any]
    paths: dict[str, Path]
    payloads: dict[str, dict[str, Any]]

    @property
    def interaction_scope(self) -> str:
        return str(self.authority["interaction_scope"])


def validate_authority(authority_path: Path) -> BaselineEvidence:
    authority_path = authority_path.resolve()
    authority = _load_object(authority_path)
    if authority.get("schema") != "CanonicalSchema10PlotBaselineV1" or authority.get("status") != "PASS":
        raise CanonicalBaselineError("canonical schema10 plot authority schema/status differs")
    if authority.get("collision_system") != "pp":
        raise CanonicalBaselineError("canonical schema10 plot authority collision system differs")
    if authority.get("interaction_scope") != "SINGLE_INTERACTION_ONLY":
        raise CanonicalBaselineError("canonical pp schema10 interaction scope is not SI-only")
    if authority.get("weighting_contract", {}).get("schema") != WEIGHTING_SCHEMA:
        raise CanonicalBaselineError("canonical schema10 weighting authority differs")
    if authority.get("stitching_acceptance_contract") != {
        "criterion": "same-bin density overlap between neighboring generator slices on both sides of every ownership boundary",
        "fit_residual_is_acceptance_gate": False,
        "maximum_fractional_deviation": 0.15,
        "visual_flatness_is_acceptance_gate": False,
    }:
        raise CanonicalBaselineError("canonical schema10 stitching acceptance contract differs")

    bindings = authority.get("bindings")
    if not isinstance(bindings, Mapping):
        raise CanonicalBaselineError("canonical schema10 authority bindings are missing")
    required = {
        "accepted_catalog",
        "raw_stitching_assembly",
        "raw_stitching_assembly_receipt",
        "interaction_scope_decision",
        "source_scope_observation",
    }
    if set(bindings) != required:
        raise CanonicalBaselineError(
            f"canonical schema10 authority binding inventory differs: {sorted(set(bindings) ^ required)}"
        )

    paths: dict[str, Path] = {}
    payloads: dict[str, dict[str, Any]] = {}
    for name in sorted(required):
        row = bindings[name]
        if not isinstance(row, Mapping):
            raise CanonicalBaselineError(f"bindings.{name} is not an object")
        path = _resolve_repo_path(row.get("path"), f"bindings.{name}")
        expected_hash = row.get("sha256")
        actual_hash = sha256_file(path)
        if expected_hash != actual_hash:
            raise CanonicalBaselineError(
                f"canonical schema10 binding drifted for {name}: expected={expected_hash}, actual={actual_hash}"
            )
        paths[name] = path
        payloads[name] = _load_object(path)

    expected = authority.get("expected_contract", {})
    catalog = payloads["accepted_catalog"]
    raw = payloads["raw_stitching_assembly"]
    raw_receipt = payloads["raw_stitching_assembly_receipt"]
    decision = payloads["interaction_scope_decision"]
    observation = payloads["source_scope_observation"]
    if catalog.get("schema") != expected.get("accepted_catalog_schema") or catalog.get("status") != "PASS":
        raise CanonicalBaselineError("accepted schema10 catalog contract differs")
    if raw.get("schema") != expected.get("raw_assembly_schema") or raw.get("status") != "PASS":
        raise CanonicalBaselineError("raw schema10 stitching assembly contract differs")
    for field in (
        "input_root_count",
        "unique_input_root_count",
        "bulk_product_count",
        "foreground_product_count",
        "sample_count",
    ):
        if raw.get(field) != expected.get(field):
            raise CanonicalBaselineError(f"raw schema10 stitching {field} differs")
    if raw.get("reader_contract", {}).get("source") != expected.get("reader_source"):
        raise CanonicalBaselineError("raw schema10 reader source differs")
    if raw.get("reader_contract", {}).get("dst_reads") != 0:
        raise CanonicalBaselineError("raw schema10 reader unexpectedly used DSTs")
    if raw.get("reader_contract", {}).get("ttree_regeneration") is not False:
        raise CanonicalBaselineError("raw schema10 reader unexpectedly regenerated TTrees")
    if raw.get("accepted_catalog_sha256") != sha256_file(paths["accepted_catalog"]):
        raise CanonicalBaselineError("raw assembly is not bound to the accepted schema10 catalog")
    if raw_receipt.get("assembled_sha256") != sha256_file(paths["raw_stitching_assembly"]):
        raise CanonicalBaselineError("raw assembly receipt does not bind the raw schema10 assembly")
    if set(catalog.get("samples", {})) != set(raw.get("samples", {})):
        raise CanonicalBaselineError("catalog/raw schema10 sample inventories differ")
    if set(expected.get("pp_samples", [])) != EXPECTED_PP_SAMPLES:
        raise CanonicalBaselineError("authority pp sample inventory differs")
    if not EXPECTED_PP_SAMPLES.issubset(raw.get("samples", {})):
        raise CanonicalBaselineError("raw assembly is missing canonical pp schema10 samples")

    if (
        decision.get("schema") != "THE121THE122InteractionScopeDecisionV1"
        or decision.get("status") != "ACCEPTED_APPEND_ONLY_SCOPE_DECISION"
        or decision.get("decision", {}).get("current_pp_production_scope") != "SINGLE_INTERACTION_ONLY"
        or decision.get("current_output_contract", {}).get("interaction_label_required") != "SI"
    ):
        raise CanonicalBaselineError("accepted pp interaction-scope decision differs")
    if (
        observation.get("schema") != "THE121THE122FullSISourceObservationV1"
        or observation.get("status") != "PASS_EXACT_SINGLE_INTERACTION_SOURCE_BOUNDARY"
        or observation.get("scope_decision", {}).get("sha256") != sha256_file(paths["interaction_scope_decision"])
        or observation.get("coverage", {}).get("pp_sample_count") != len(EXPECTED_PP_SAMPLES)
    ):
        raise CanonicalBaselineError("pp schema10 SI source observation differs")
    pp_source_rows = [row for row in observation.get("samples", []) if row.get("system") == "pp"]
    if len(pp_source_rows) != len(EXPECTED_PP_SAMPLES) or any(
        row.get("interaction") != "si" for row in pp_source_rows
    ):
        raise CanonicalBaselineError("pp schema10 source observation is not exactly SI")

    return BaselineEvidence(
        authority_path=authority_path,
        authority_sha256=sha256_file(authority_path),
        authority=authority,
        paths=paths,
        payloads=payloads,
    )


def _arrays_close(actual: Any, expected: Any, label: str) -> None:
    try:
        left = np.asarray(actual, dtype=float)
        right = np.asarray(expected, dtype=float)
    except (TypeError, ValueError) as exc:
        raise CanonicalBaselineError(f"{label} is not numeric") from exc
    if left.shape != right.shape or not np.all(np.isfinite(left)) or not np.all(np.isfinite(right)):
        raise CanonicalBaselineError(f"{label} shape or finiteness differs")
    if not np.allclose(left, right, rtol=2.0e-12, atol=1.0e-18):
        raise CanonicalBaselineError(f"{label} differs from canonical dual-channel normalization")


def validate_corrected_assembly(
    assembly_path: Path,
    *,
    authority_path: Path,
    embedded_inclusive_assembly_path: Path,
    embedded_inclusive_receipt_path: Path,
) -> tuple[dict[str, Any], dict[str, Any]]:
    from scripts.data_prep.recoiljets.auau_embedded_inclusive_schema10_weighting import (
        load_source_stitch_artifact,
    )

    evidence = validate_authority(authority_path)
    embedded_source = load_source_stitch_artifact(
        embedded_inclusive_assembly_path,
        embedded_inclusive_receipt_path,
    )
    embedded_assembly = _load_object(embedded_inclusive_assembly_path)
    assembly_path = assembly_path.resolve()
    assembly = _load_object(assembly_path)
    if assembly.get("schema") != CORRECTED_ASSEMBLY_SCHEMA or assembly.get("status") != "PASS":
        raise CanonicalBaselineError("corrected schema10 assembly schema/status differs")
    expected = evidence.authority["expected_contract"]
    for field in (
        "input_root_count",
        "unique_input_root_count",
        "bulk_product_count",
        "foreground_product_count",
        "sample_count",
    ):
        if assembly.get(field) != expected.get(field):
            raise CanonicalBaselineError(f"corrected schema10 assembly {field} differs")
    if assembly.get("accepted_catalog_sha256") != sha256_file(evidence.paths["accepted_catalog"]):
        raise CanonicalBaselineError("corrected assembly catalog binding differs")
    normalization = assembly.get("normalization")
    if (
        not isinstance(normalization, Mapping)
        or normalization.get("schema") != "THE248Schema10SplitStitchingNormalizationV4"
        or normalization.get("pp_and_noninclusive_contract") != normalization_contract()
        or normalization.get("auau_embedded_inclusive", {}).get("source_artifact", {}).get(
            "assembly_sha256"
        )
        != embedded_source.assembly_sha256
        or normalization.get("auau_embedded_inclusive", {}).get(
            "centrality_reweighting_applied"
        )
        is not False
    ):
        raise CanonicalBaselineError("corrected split-normalization contract differs")
    if assembly.get("source_raw_assembly", {}).get("sha256") != sha256_file(
        evidence.paths["raw_stitching_assembly"]
    ):
        raise CanonicalBaselineError("corrected assembly does not bind the canonical raw assembly")
    raw = evidence.payloads["raw_stitching_assembly"]
    if assembly.get("products") != raw.get("products"):
        raise CanonicalBaselineError("corrected assembly product/hash inventory differs from bound input")
    if set(assembly.get("samples", {})) != set(raw.get("samples", {})):
        raise CanonicalBaselineError("corrected/raw schema10 sample inventories differ")

    catalog_samples = evidence.payloads["accepted_catalog"]["samples"]
    for sample_id, actual in assembly["samples"].items():
        source = raw["samples"][sample_id]
        for field in ("bin_edges_gev", "raw_counts", "sumw", "sumw2"):
            if actual.get(field) != source.get(field):
                raise CanonicalBaselineError(f"{sample_id}.{field} differs from bound raw assembly")
        catalog_sample = catalog_samples[sample_id]
        embedded_inclusive = (
            catalog_sample.get("system") == "auau"
            and catalog_sample.get("source_class") == "inclusive_background"
        )
        if embedded_inclusive:
            embedded_source.sample(sample_id)
            expected_embedded = embedded_assembly["samples"].get(sample_id)
            if not isinstance(expected_embedded, Mapping):
                raise CanonicalBaselineError(
                    f"canonical embedded-inclusive sample missing: {sample_id}"
                )
            for field in (
                "generator_stitching_density_pb_per_gev",
                "generator_stitching_density_sumw2_pb2_per_gev2",
            ):
                _arrays_close(actual.get(field), expected_embedded.get(field), f"{sample_id}.{field}")
            for field in (
                "normalization_schema",
                "normalization_denominator_events",
                "normalization_denominator_kind",
                "ownership_effective_cross_section_pb",
                "ownership_window",
                "stitching_weight_pb_per_owned_event",
                "generator_stitching_channel",
                "same_observable_denominator_proof",
            ):
                if actual.get(field) != expected_embedded.get(field):
                    raise CanonicalBaselineError(f"{sample_id}.{field} differs")
            if actual.get("analysis_weight_state") != (
                "CANONICAL_SOURCE_STITCH_ONLY__CENTRALITY_NOT_APPLIED"
            ) or actual.get("nominal_downstream_analysis_ready") is not False:
                raise CanonicalBaselineError(
                    f"{sample_id} must remain explicitly pre-centrality in the stitching assembly"
                )
            for forbidden in (
                "cross_section_weight_pb_per_event",
                "analysis_weighted_channel",
                "analysis_weighted_density_pb_per_gev",
                "analysis_weighted_density_sumw2_pb2_per_gev2",
            ):
                if forbidden in actual:
                    raise CanonicalBaselineError(
                        f"{sample_id} exposes forbidden legacy normalization field {forbidden}"
                    )
            continue
        expected_sample = normalize_stitching_sample(
            source,
            catalog_sample,
            bin_width_gev=float(assembly["binning"]["width_gev"]),
        )
        if actual.get("ownership_window") != source.get("ownership_window"):
            raise CanonicalBaselineError(f"{sample_id}.ownership_window differs from bound raw assembly")
        for field in (
            "generator_stitching_density_pb_per_gev",
            "generator_stitching_density_sumw2_pb2_per_gev2",
            "analysis_weighted_density_pb_per_gev",
            "analysis_weighted_density_sumw2_pb2_per_gev2",
        ):
            _arrays_close(actual.get(field), expected_sample.get(field), f"{sample_id}.{field}")
        for field in (
            "normalization_schema",
            "generator_stitching_channel",
            "analysis_weighted_channel",
        ):
            if actual.get(field) != expected_sample.get(field):
                raise CanonicalBaselineError(f"{sample_id}.{field} differs")
        for legacy_field in (
            "weighted_density_pb_per_gev",
            "weighted_density_sumw2_pb2_per_gev2",
        ):
            if legacy_field in actual:
                raise CanonicalBaselineError(f"{sample_id} exposes ambiguous {legacy_field}")

    receipt = {
        "schema": "CanonicalSchema10PlotBaselineGuardReceiptV1",
        "status": "PASS",
        "baseline_id": evidence.authority["baseline_id"],
        "authority_path": evidence.authority_path.as_posix(),
        "authority_sha256": evidence.authority_sha256,
        "corrected_assembly_path": assembly_path.as_posix(),
        "corrected_assembly_sha256": sha256_file(assembly_path),
        "accepted_catalog_path": evidence.paths["accepted_catalog"].as_posix(),
        "accepted_catalog_sha256": sha256_file(evidence.paths["accepted_catalog"]),
        "raw_assembly_sha256": sha256_file(evidence.paths["raw_stitching_assembly"]),
        "interaction_scope": evidence.interaction_scope,
        "pp_sample_count": len(EXPECTED_PP_SAMPLES),
        "input_root_count": assembly["input_root_count"],
        "weighting_schema": WEIGHTING_SCHEMA,
        "auau_embedded_inclusive_source_stitch": embedded_source.binding_payload(),
        "stitching_plot_channel": "generator_stitching",
        "stitching_acceptance_contract": evidence.authority[
            "stitching_acceptance_contract"
        ],
        "si_plus_di_claim_authorized": False,
    }
    return assembly, receipt


def assert_nominal_label(label: str, receipt: Mapping[str, Any]) -> None:
    normalized = " ".join(label.lower().replace("_", " ").split())
    if "schema10" not in normalized or "si" not in normalized:
        raise CanonicalBaselineError("nominal label must explicitly say canonical schema10 SI")
    if "si+di" in normalized or "si + di" in normalized or "double interaction" in normalized:
        raise CanonicalBaselineError("SI+DI label is not authorized by the accepted schema10 baseline")
    if receipt.get("interaction_scope") != "SINGLE_INTERACTION_ONLY":
        raise CanonicalBaselineError("nominal label receipt interaction scope differs")


def validate_auau_embedded_inclusive_assembly(
    assembly_path: Path,
    authority_path: Path,
    validation_receipt_path: Path,
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Route nominal embedded-inclusive plots through their dedicated authority."""
    from scripts.data_prep.recoiljets.auau_embedded_inclusive_schema10_weighting import (
        EmbeddedInclusiveWeightError,
        load_source_stitch_artifact,
    )

    try:
        source = load_source_stitch_artifact(assembly_path, validation_receipt_path)
        assembly = _load_object(assembly_path)
        receipt = _load_object(validation_receipt_path)
    except EmbeddedInclusiveWeightError as error:
        raise CanonicalBaselineError(str(error)) from error
    if source.authority_path != authority_path.resolve():
        raise CanonicalBaselineError("embedded-inclusive authority path differs")
    if (
        receipt.get("status") != "pass"
        or receipt.get("centrality_reweighting_applied") is not False
    ):
        raise CanonicalBaselineError("embedded-inclusive canonical guard differs")
    return assembly, receipt
