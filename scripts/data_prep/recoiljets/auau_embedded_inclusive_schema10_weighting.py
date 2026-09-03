#!/usr/bin/env python3
"""Canonical Au+Au schema10 embedded-inclusive analysis weights.

This module is the only downstream entry point for embedded inclusive-jet
weights.  It consumes the accepted source-stitch artifact, applies its
``sigma_eff / Npass`` factor to the producer event weight, and then applies a
READY centrality contract.  Plotters and correction code consume the emitted
payload and receipt; they must never reconstruct either factor.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import math
import os
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

try:
    from scripts.data_prep.recoiljets.auau_centrality_weight_contract import (
        APPLICATION_ORDER,
        CANONICAL_CONTRACT_STATUS,
        CentralityWeightContract,
        EventWeightState,
        load_contract_receipt,
    )
except ModuleNotFoundError:  # Support imports with ``scripts`` itself on sys.path.
    from .auau_centrality_weight_contract import (  # type: ignore[no-redef]
        APPLICATION_ORDER,
        CANONICAL_CONTRACT_STATUS,
        CentralityWeightContract,
        EventWeightState,
        load_contract_receipt,
    )


SOURCE_ASSEMBLY_SCHEMA = "CanonicalSchema10EmbeddedInclusiveStitchingAssemblyV1"
SOURCE_VALIDATION_SCHEMA = "CanonicalSchema10EmbeddedInclusiveStitchingValidationV1"
SOURCE_WEIGHTING_SCHEMA = "CanonicalSchema10EmbeddedInclusiveStitchingWeightingV1"
ANALYSIS_PAYLOAD_SCHEMA = "CanonicalSchema10AuAuEmbeddedInclusiveAnalysisWeightPayloadV1"
ANALYSIS_RECEIPT_SCHEMA = "CanonicalSchema10AuAuEmbeddedInclusiveAnalysisWeightReceiptV1"
ANALYSIS_STATUS = "COMPLETE_PRODUCER_STITCH_CENTRALITY"
ANALYSIS_SCOPE = "auau_embedded_inclusive"
SIMULATION_FAMILY = "inclusive_background"
DOWNSTREAM_RECEIPT_SCHEMA = (
    "CanonicalSchema10AuAuEmbeddedInclusiveDownstreamArtifactReceiptV1"
)
DOWNSTREAM_STATUS = "PASS_HASH_BOUND_ANALYSIS_WEIGHT_AND_OUTPUTS"
SAMPLE_IDS = ("auau_jet12", "auau_jet20", "auau_jet30", "auau_jet40")
SAMPLE_FAMILY = "inclusivejet"


class EmbeddedInclusiveWeightError(ValueError):
    """An embedded-inclusive weight artifact or application is invalid."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _sha256_payload(payload: Mapping[str, Any]) -> str:
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _load_object(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise EmbeddedInclusiveWeightError(f"invalid JSON artifact: {path}") from error
    if not isinstance(value, dict):
        raise EmbeddedInclusiveWeightError(f"JSON object required: {path}")
    return value


def _finite(value: Any, label: str, *, positive: bool = False) -> float:
    if isinstance(value, bool):
        raise EmbeddedInclusiveWeightError(f"{label} is not numeric")
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise EmbeddedInclusiveWeightError(f"{label} is not numeric") from error
    if not math.isfinite(result) or (positive and result <= 0.0):
        raise EmbeddedInclusiveWeightError(f"{label} is invalid")
    return result


def _is_sha256(value: Any) -> bool:
    if not isinstance(value, str) or len(value) != 64 or value != value.lower():
        return False
    try:
        int(value, 16)
    except ValueError:
        return False
    return True


@dataclass(frozen=True)
class SourceSampleContract:
    sample_id: str
    ownership_effective_cross_section_pb: float
    normalization_denominator_events: int
    ownership_low_gev: float
    ownership_high_gev: float | None
    stitching_weight_pb_per_owned_event: float

    def owns(self, maximum_truth_jet_pt_gev: float) -> bool:
        value = _finite(maximum_truth_jet_pt_gev, "maximum_truth_jet_pt_gev")
        return value >= self.ownership_low_gev and (
            self.ownership_high_gev is None or value < self.ownership_high_gev
        )

    def payload(self) -> dict[str, Any]:
        return {
            "sample_id": self.sample_id,
            "ownership_effective_cross_section_pb": self.ownership_effective_cross_section_pb,
            "normalization_denominator_events": self.normalization_denominator_events,
            "ownership_window_gev": [self.ownership_low_gev, self.ownership_high_gev],
            "stitching_weight_pb_per_owned_event": self.stitching_weight_pb_per_owned_event,
        }


@dataclass(frozen=True)
class SourceStitchArtifact:
    assembly_path: Path
    assembly_sha256: str
    validation_receipt_path: Path
    validation_receipt_sha256: str
    authority_path: Path
    authority_sha256: str
    samples: tuple[SourceSampleContract, ...]

    def sample(self, sample_id: str) -> SourceSampleContract:
        for sample in self.samples:
            if sample.sample_id == sample_id:
                return sample
        raise EmbeddedInclusiveWeightError(
            f"unsupported embedded inclusive sample {sample_id!r}; expected {list(SAMPLE_IDS)}"
        )

    def binding_payload(self) -> dict[str, Any]:
        return {
            "assembly_path": str(self.assembly_path.resolve()),
            "assembly_sha256": self.assembly_sha256,
            "validation_receipt_path": str(self.validation_receipt_path.resolve()),
            "validation_receipt_sha256": self.validation_receipt_sha256,
            "authority_path": str(self.authority_path.resolve()),
            "authority_sha256": self.authority_sha256,
            "sample_ids": list(SAMPLE_IDS),
            "samples": {sample.sample_id: sample.payload() for sample in self.samples},
        }


def load_source_stitch_artifact(
    assembly_path: Path,
    validation_receipt_path: Path,
    *,
    verify_authority_file: bool = True,
) -> SourceStitchArtifact:
    """Load and fully validate the accepted ``sigma_eff/Npass`` authority chain."""
    assembly_path = Path(assembly_path).resolve()
    validation_receipt_path = Path(validation_receipt_path).resolve()
    assembly = _load_object(assembly_path)
    receipt = _load_object(validation_receipt_path)
    if assembly.get("schema") != SOURCE_ASSEMBLY_SCHEMA or assembly.get("status") != "PASS":
        raise EmbeddedInclusiveWeightError("source-stitch assembly is not accepted")
    if (
        receipt.get("schema") != SOURCE_VALIDATION_SCHEMA
        or receipt.get("status") != "pass"
        or receipt.get("output_schema") != SOURCE_ASSEMBLY_SCHEMA
    ):
        raise EmbeddedInclusiveWeightError("source-stitch validation receipt is not accepted")
    assembly_sha = sha256_file(assembly_path)
    if receipt.get("assembly_sha256") != assembly_sha:
        raise EmbeddedInclusiveWeightError("source-stitch assembly SHA-256 differs from receipt")
    if Path(str(receipt.get("assembly_path"))).resolve() != assembly_path:
        raise EmbeddedInclusiveWeightError("source-stitch assembly path differs from receipt")

    authority_path = Path(str(assembly.get("authority_path"))).resolve()
    authority_sha = assembly.get("authority_sha256")
    if (
        not _is_sha256(authority_sha)
        or receipt.get("authority_sha256") != authority_sha
        or Path(str(receipt.get("authority_path"))).resolve() != authority_path
    ):
        raise EmbeddedInclusiveWeightError("source-stitch authority binding differs")
    if verify_authority_file and (
        not authority_path.is_file() or sha256_file(authority_path) != authority_sha
    ):
        raise EmbeddedInclusiveWeightError("source-stitch authority file is missing or stale")

    normalization = assembly.get("normalization")
    if not isinstance(normalization, Mapping) or (
        normalization.get("schema") != SOURCE_WEIGHTING_SCHEMA
        or normalization.get("centrality_reweighting_applied") is not False
        or normalization.get("denominator_kind")
        != "current_schema10_events_in_same_truth_jet_ownership_window"
    ):
        raise EmbeddedInclusiveWeightError("source-stitch normalization contract differs")
    if assembly.get("sample_family") != "inclusive_background":
        raise EmbeddedInclusiveWeightError("source-stitch sample family differs")
    if tuple(assembly.get("sample_ids", ())) != SAMPLE_IDS or tuple(receipt.get("sample_ids", ())) != SAMPLE_IDS:
        raise EmbeddedInclusiveWeightError("source-stitch sample inventory differs")

    rows = assembly.get("samples")
    denominators = receipt.get("normalization_denominators")
    if not isinstance(rows, Mapping) or not isinstance(denominators, Mapping):
        raise EmbeddedInclusiveWeightError("source-stitch sample contracts are missing")
    samples: list[SourceSampleContract] = []
    previous_high: float | None = None
    for index, sample_id in enumerate(SAMPLE_IDS):
        row = rows.get(sample_id)
        if not isinstance(row, Mapping):
            raise EmbeddedInclusiveWeightError(f"source-stitch sample is missing: {sample_id}")
        window = row.get("ownership_window")
        proof = row.get("same_observable_denominator_proof")
        channel = row.get("generator_stitching_channel")
        if not all(isinstance(value, Mapping) for value in (window, proof, channel)):
            raise EmbeddedInclusiveWeightError(f"{sample_id} source-stitch proof is missing")
        assert isinstance(window, Mapping) and isinstance(proof, Mapping) and isinstance(channel, Mapping)
        if (
            row.get("sample_id") != sample_id
            or row.get("system") != "auau"
            or row.get("source_class") != "inclusive_background"
            or row.get("normalization_schema") != SOURCE_WEIGHTING_SCHEMA
            or window.get("lower_edge_inclusive") is not True
            or window.get("upper_edge_inclusive") is not False
            or channel.get("centrality_reweighting_applied") is not False
        ):
            raise EmbeddedInclusiveWeightError(f"{sample_id} source-stitch identity differs")
        denominator = row.get("normalization_denominator_events")
        if type(denominator) is not int or denominator <= 0 or denominators.get(sample_id) != denominator:
            raise EmbeddedInclusiveWeightError(f"{sample_id} Npass binding differs")
        if any(proof.get(field) != denominator for field in (
            "event_level_ownership_audit_half_open_count", "histogram_owned_raw_count"
        )) or any(proof.get(field) != 0 for field in (
            "source_missing_raw_count", "source_underflow_raw_count", "source_overflow_raw_count"
        )) or proof.get("source_sumw_equals_raw_counts") is not True or proof.get("source_sumw2_equals_raw_counts") is not True:
            raise EmbeddedInclusiveWeightError(f"{sample_id} same-observable denominator proof differs")
        xsec = _finite(row.get("ownership_effective_cross_section_pb"), f"{sample_id}.sigma_eff", positive=True)
        weight = _finite(row.get("stitching_weight_pb_per_owned_event"), f"{sample_id}.stitch_weight", positive=True)
        expected = xsec / denominator
        if not math.isclose(weight, expected, rel_tol=2e-14, abs_tol=0.0):
            raise EmbeddedInclusiveWeightError(f"{sample_id} stitch weight is not sigma_eff/Npass")
        generated = row.get("generated_events")
        if type(generated) is int and generated > 0 and math.isclose(
            weight, xsec / generated, rel_tol=2e-14, abs_tol=0.0
        ):
            raise EmbeddedInclusiveWeightError(f"{sample_id} incorrectly uses sigma_eff/Ngen")
        edges = row.get("bin_edges_gev")
        raw_counts = row.get("raw_counts")
        density = row.get("generator_stitching_density_pb_per_gev")
        variance = row.get("generator_stitching_density_sumw2_pb2_per_gev2")
        if (
            not isinstance(edges, list)
            or not isinstance(raw_counts, list)
            or not isinstance(density, list)
            or not isinstance(variance, list)
            or len(edges) != len(raw_counts) + 1
            or len(raw_counts) != len(density)
            or len(density) != len(variance)
            or not raw_counts
        ):
            raise EmbeddedInclusiveWeightError(f"{sample_id} source-stitch histogram shape differs")
        widths = [
            _finite(right, f"{sample_id}.edge") - _finite(left, f"{sample_id}.edge")
            for left, right in zip(edges[:-1], edges[1:])
        ]
        if any(width <= 0.0 for width in widths) or any(
            not math.isclose(width, widths[0], rel_tol=0.0, abs_tol=1e-12)
            for width in widths
        ):
            raise EmbeddedInclusiveWeightError(f"{sample_id} source-stitch bin widths differ")
        density_scale = weight / widths[0]
        if not math.isclose(
            _finite(channel.get("density_scale_pb_per_gev_per_raw_event"), f"{sample_id}.density_scale"),
            density_scale,
            rel_tol=2e-14,
            abs_tol=0.0,
        ):
            raise EmbeddedInclusiveWeightError(f"{sample_id} density scale differs from sigma_eff/Npass")
        for bin_index, (count, value, variance_value) in enumerate(
            zip(raw_counts, density, variance)
        ):
            if type(count) is not int or count < 0:
                raise EmbeddedInclusiveWeightError(f"{sample_id} raw count {bin_index} is invalid")
            if not math.isclose(
                _finite(value, f"{sample_id}.density[{bin_index}]"),
                count * density_scale,
                rel_tol=2e-12,
                abs_tol=1e-18,
            ) or not math.isclose(
                _finite(variance_value, f"{sample_id}.variance[{bin_index}]"),
                count * density_scale * density_scale,
                rel_tol=2e-12,
                abs_tol=1e-18,
            ):
                raise EmbeddedInclusiveWeightError(
                    f"{sample_id} source-stitch density moments differ in bin {bin_index}"
                )
        low = _finite(window.get("low_gev"), f"{sample_id}.ownership_low")
        high_raw = window.get("high_gev")
        high = None if high_raw is None else _finite(high_raw, f"{sample_id}.ownership_high")
        if (index == len(SAMPLE_IDS) - 1) != (high is None):
            raise EmbeddedInclusiveWeightError(f"{sample_id} final-window semantics differ")
        if previous_high is not None and not math.isclose(low, previous_high, rel_tol=0.0, abs_tol=1e-12):
            raise EmbeddedInclusiveWeightError(f"{sample_id} ownership windows are not contiguous")
        if high is not None and high <= low:
            raise EmbeddedInclusiveWeightError(f"{sample_id} ownership window is invalid")
        previous_high = high
        samples.append(SourceSampleContract(sample_id, xsec, denominator, low, high, weight))
    return SourceStitchArtifact(
        assembly_path,
        assembly_sha,
        validation_receipt_path,
        sha256_file(validation_receipt_path),
        authority_path,
        authority_sha,
        tuple(samples),
    )


@dataclass(frozen=True)
class CompleteAnalysisWeight:
    value: float
    components: tuple[str, ...]
    source_stitch_weight: float
    centrality_weight: float


@dataclass(frozen=True)
class EmbeddedInclusiveAnalysisWeightProvider:
    source: SourceStitchArtifact
    centrality: CentralityWeightContract
    centrality_receipt_path: Path
    centrality_receipt_sha256: str

    @classmethod
    def load(
        cls,
        assembly_path: Path,
        source_validation_receipt_path: Path,
        centrality_receipt_path: Path,
        *,
        expected_centrality_dependency_fingerprint: str,
        verify_dependency_files: bool = False,
    ) -> "EmbeddedInclusiveAnalysisWeightProvider":
        source = load_source_stitch_artifact(assembly_path, source_validation_receipt_path)
        centrality_path = Path(centrality_receipt_path).resolve()
        centrality = load_contract_receipt(
            centrality_path,
            expected_dependency_fingerprint=expected_centrality_dependency_fingerprint,
            require_ready=True,
            verify_dependency_files=verify_dependency_files,
        )
        if centrality.status != CANONICAL_CONTRACT_STATUS:
            raise EmbeddedInclusiveWeightError("centrality contract is not READY")
        centrality.family_map(SAMPLE_FAMILY)
        return cls(source, centrality, centrality_path, sha256_file(centrality_path))

    @property
    def dependency_fingerprint(self) -> str:
        return _sha256_payload(
            {
                "source_stitch": self.source.binding_payload(),
                "centrality_receipt_path": str(self.centrality_receipt_path),
                "centrality_receipt_sha256": self.centrality_receipt_sha256,
                "centrality_dependency_fingerprint": self.centrality.dependency_fingerprint,
                "centrality_contract_fingerprint": self.centrality.contract_fingerprint,
                "application_order": list(APPLICATION_ORDER),
            }
        )

    def apply(
        self,
        *,
        sample_id: str,
        producer_event_weight: float,
        maximum_truth_jet_pt_gev: float,
        centrality_percent: float,
    ) -> CompleteAnalysisWeight:
        producer = _finite(producer_event_weight, "producer_event_weight")
        if producer < 0.0:
            raise EmbeddedInclusiveWeightError("producer_event_weight is negative")
        sample = self.source.sample(sample_id)
        if not sample.owns(maximum_truth_jet_pt_gev):
            raise EmbeddedInclusiveWeightError(
                f"{sample_id} event is outside its half-open truth-jet ownership window"
            )
        stitched = producer * sample.stitching_weight_pb_per_owned_event
        state = EventWeightState.after_canonical_stitch(stitched)
        complete = self.centrality.apply_once(
            state,
            centrality_percent,
            SAMPLE_FAMILY,
            expected_dependency_fingerprint=self.centrality.dependency_fingerprint,
        )
        centrality_weight = complete.value / stitched if stitched > 0.0 else self.centrality.weight(
            centrality_percent, SAMPLE_FAMILY
        )
        return CompleteAnalysisWeight(
            complete.value,
            complete.components,
            sample.stitching_weight_pb_per_owned_event,
            centrality_weight,
        )

    def provenance_payload(self) -> dict[str, Any]:
        return {
            "analysis_scope": ANALYSIS_SCOPE,
            "sample_family": SIMULATION_FAMILY,
            "sample_ids": list(SAMPLE_IDS),
            "state": ANALYSIS_STATUS,
            "application_order": list(APPLICATION_ORDER),
            "dependency_fingerprint": self.dependency_fingerprint,
            "source_stitch": self.source.binding_payload(),
            "centrality_contract": {
                "path": str(self.centrality_receipt_path),
                "sha256": self.centrality_receipt_sha256,
                "dependency_fingerprint": self.centrality.dependency_fingerprint,
                "contract_fingerprint": self.centrality.contract_fingerprint,
                "sample_family": SAMPLE_FAMILY,
            },
        }


def build_reco_cluster_source_fraction_payload(
    rows: Iterable[Mapping[str, Any]],
    provider: EmbeddedInclusiveAnalysisWeightProvider,
    *,
    reco_cluster_pt_edges_gev: Sequence[float],
    input_bindings: Sequence[Mapping[str, Any]],
) -> dict[str, Any]:
    """Build plot-consumable 1D source moments from complete event weights."""
    edges = tuple(_finite(value, "reco_cluster_pt_edge") for value in reco_cluster_pt_edges_gev)
    if len(edges) < 2 or any(right <= left for left, right in zip(edges[:-1], edges[1:])):
        raise EmbeddedInclusiveWeightError("reco-cluster pT edges must be strictly increasing")
    width = edges[1] - edges[0]
    if any(not math.isclose(right - left, width, rel_tol=0.0, abs_tol=1e-12) for left, right in zip(edges[:-1], edges[1:])):
        raise EmbeddedInclusiveWeightError("reco-cluster pT bins must have one constant width")
    if not math.isclose(width, 1.0, rel_tol=0.0, abs_tol=1e-12):
        raise EmbeddedInclusiveWeightError("plot-source payload requires exact 1 GeV reco-cluster pT bins")
    validated_inputs: list[dict[str, Any]] = []
    if not input_bindings:
        raise EmbeddedInclusiveWeightError("plot-source payload requires hash-bound input artifacts")
    seen_paths: set[str] = set()
    for index, binding in enumerate(input_bindings):
        path = Path(str(binding.get("path", ""))).resolve()
        expected_sha = binding.get("sha256")
        size = binding.get("size_bytes")
        if (
            not path.is_file()
            or str(path) in seen_paths
            or not _is_sha256(expected_sha)
            or sha256_file(path) != expected_sha
            or type(size) is not int
            or size != path.stat().st_size
        ):
            raise EmbeddedInclusiveWeightError(f"input binding {index} is missing, duplicate, or stale")
        seen_paths.add(str(path))
        validated_inputs.append({"path": str(path), "sha256": expected_sha, "size_bytes": size})
    nbins = len(edges) - 1
    sums = {sample_id: [[0.0, 0.0, 0] for _ in range(nbins)] for sample_id in SAMPLE_IDS}
    accepted = 0
    for row in rows:
        sample_id = str(row.get("sample_id", ""))
        reco_pt = _finite(row.get("reco_cluster_pt_gev"), "reco_cluster_pt_gev")
        if reco_pt < edges[0] or reco_pt >= edges[-1]:
            continue
        result = provider.apply(
            sample_id=sample_id,
            producer_event_weight=row.get("producer_event_weight"),
            maximum_truth_jet_pt_gev=row.get("maximum_truth_jet_pt_gev"),
            centrality_percent=row.get("centrality_percent"),
        )
        index = min(int((reco_pt - edges[0]) / width), nbins - 1)
        sums[sample_id][index][0] += result.value
        sums[sample_id][index][1] += result.value * result.value
        sums[sample_id][index][2] += 1
        accepted += 1
    if accepted == 0:
        raise EmbeddedInclusiveWeightError(
            "plot-source payload has no accepted rows in the requested reco-cluster pT support"
        )
    total_sumw = [sum(sums[sample_id][i][0] for sample_id in SAMPLE_IDS) for i in range(nbins)]
    total_sumw2 = [sum(sums[sample_id][i][1] for sample_id in SAMPLE_IDS) for i in range(nbins)]
    sources: dict[str, Any] = {}
    for sample_id in SAMPLE_IDS:
        sample_sumw = [row[0] for row in sums[sample_id]]
        sources[sample_id] = {
            "sumw": sample_sumw,
            "sumw2": [row[1] for row in sums[sample_id]],
            "raw_count": [row[2] for row in sums[sample_id]],
            "fraction": [
                value / total if total > 0.0 else None
                for value, total in zip(sample_sumw, total_sumw)
            ],
        }
    return {
        "schema": ANALYSIS_PAYLOAD_SCHEMA,
        "status": "PASS",
        "analysis_scope": ANALYSIS_SCOPE,
        "sample_family": SIMULATION_FAMILY,
        "sample_ids": list(SAMPLE_IDS),
        "weight_state": ANALYSIS_STATUS,
        "collision_system": "Au+Au embedded simulation",
        "observable": "reconstructed cluster pT source fraction",
        "reco_cluster_pt_edges_gev": list(edges),
        "bin_width_gev": 1.0,
        "application_order": list(APPLICATION_ORDER),
        "accepted_row_count": accepted,
        "inputs": validated_inputs,
        "sumw": total_sumw,
        "sumw2": total_sumw2,
        "sources": sources,
        "weight_provenance": provider.provenance_payload(),
    }


def _write_exclusive(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o640)
    with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
        json.dump(payload, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")


def _validate_file_binding(
    binding: Mapping[str, Any],
    label: str,
    *,
    verify_file: bool,
) -> dict[str, Any]:
    if set(binding) != {"path", "sha256", "size_bytes"}:
        raise EmbeddedInclusiveWeightError(f"{label} file-binding fields differ")
    raw_path = binding.get("path")
    expected_sha = binding.get("sha256")
    expected_size = binding.get("size_bytes")
    if not isinstance(raw_path, str) or not raw_path:
        raise EmbeddedInclusiveWeightError(f"{label} path is invalid")
    path = Path(raw_path)
    if not path.is_absolute() or str(path.resolve()) != raw_path:
        raise EmbeddedInclusiveWeightError(f"{label} path is not canonical absolute")
    if not _is_sha256(expected_sha) or type(expected_size) is not int or expected_size <= 0:
        raise EmbeddedInclusiveWeightError(f"{label} hash or size is invalid")
    if verify_file and (
        not path.is_file()
        or path.stat().st_size != expected_size
        or sha256_file(path) != expected_sha
    ):
        raise EmbeddedInclusiveWeightError(f"{label} is missing or stale")
    return {"path": raw_path, "sha256": expected_sha, "size_bytes": expected_size}


def _validate_downstream_audit(
    audit_payload: Mapping[str, Any],
    *,
    downstream_receipt_path: Path,
    analysis_weight_receipt_path: Path,
    analysis_weight_receipt: Mapping[str, Any],
    expected_dependency_fingerprint: str,
    verify_generator_file: bool,
) -> dict[str, Any]:
    provenance = audit_payload.get("analysis_weight_provenance")
    if not isinstance(provenance, Mapping):
        raise EmbeddedInclusiveWeightError(
            "downstream audit lacks exact PlotLabelContract analysis-weight provenance"
        )
    expected_provenance = {
        "path": str(analysis_weight_receipt_path.resolve()),
        "expected_dependency_fingerprint": expected_dependency_fingerprint,
        "receipt": dict(analysis_weight_receipt),
    }
    if dict(provenance) != expected_provenance:
        raise EmbeddedInclusiveWeightError(
            "downstream audit analysis-weight provenance differs from the exact receipt"
        )
    if (
        audit_payload.get("status") != "PASS"
        or audit_payload.get("collision_system") != "Au+Au"
        or audit_payload.get("sample_kind") not in {"simulation", "mixed"}
        or not isinstance(audit_payload.get("plot_kind"), str)
        or audit_payload.get("plot_kind") in {"diagnostic", "bdt_qa"}
        or audit_payload.get("simulation_scope") != "full accepted simulation"
        or audit_payload.get("simulation_family") != SIMULATION_FAMILY
        or audit_payload.get("analysis_weight_receipt_sha256")
        != sha256_file(analysis_weight_receipt_path)
        or audit_payload.get("analysis_weight_dependency_fingerprint")
        != expected_dependency_fingerprint
        or Path(str(audit_payload.get("downstream_artifact_receipt_path"))).resolve()
        != downstream_receipt_path.resolve()
    ):
        raise EmbeddedInclusiveWeightError("downstream audit nominal provenance differs")
    generator_binding = audit_payload.get("generator")
    if not isinstance(generator_binding, Mapping):
        raise EmbeddedInclusiveWeightError("downstream audit lacks generator binding")
    return _validate_file_binding(
        generator_binding,
        "downstream generator",
        verify_file=verify_generator_file,
    )


def validate_analysis_weight_payload(
    payload: Mapping[str, Any],
    provider: EmbeddedInclusiveAnalysisWeightProvider,
    *,
    verify_input_files: bool = True,
) -> None:
    """Deeply validate every numerical and provenance invariant in a payload."""
    expected_keys = {
        "schema",
        "status",
        "analysis_scope",
        "sample_family",
        "sample_ids",
        "weight_state",
        "collision_system",
        "observable",
        "reco_cluster_pt_edges_gev",
        "bin_width_gev",
        "application_order",
        "accepted_row_count",
        "inputs",
        "sumw",
        "sumw2",
        "sources",
        "weight_provenance",
    }
    if set(payload) != expected_keys:
        raise EmbeddedInclusiveWeightError("analysis-weight payload fields differ")
    if (
        payload.get("schema") != ANALYSIS_PAYLOAD_SCHEMA
        or payload.get("status") != "PASS"
        or payload.get("analysis_scope") != ANALYSIS_SCOPE
        or payload.get("sample_family") != SIMULATION_FAMILY
        or tuple(payload.get("sample_ids", ())) != SAMPLE_IDS
        or payload.get("weight_state") != ANALYSIS_STATUS
        or payload.get("collision_system") != "Au+Au embedded simulation"
        or payload.get("observable") != "reconstructed cluster pT source fraction"
        or payload.get("application_order") != list(APPLICATION_ORDER)
    ):
        raise EmbeddedInclusiveWeightError("analysis-weight payload identity differs")

    edges_raw = payload.get("reco_cluster_pt_edges_gev")
    if not isinstance(edges_raw, list) or len(edges_raw) < 2:
        raise EmbeddedInclusiveWeightError("analysis-weight payload bin edges are invalid")
    edges = [_finite(value, "reco_cluster_pt_edge") for value in edges_raw]
    widths = [right - left for left, right in zip(edges[:-1], edges[1:])]
    if (
        any(width <= 0.0 for width in widths)
        or any(not math.isclose(width, 1.0, rel_tol=0.0, abs_tol=1e-12) for width in widths)
        or not math.isclose(
            _finite(payload.get("bin_width_gev"), "bin_width_gev"),
            1.0,
            rel_tol=0.0,
            abs_tol=1e-12,
        )
    ):
        raise EmbeddedInclusiveWeightError("analysis-weight payload requires exact 1 GeV bins")
    nbins = len(edges) - 1

    inputs = payload.get("inputs")
    if not isinstance(inputs, list) or not inputs:
        raise EmbeddedInclusiveWeightError("analysis-weight payload lacks input bindings")
    canonical_inputs: list[dict[str, Any]] = []
    for index, binding in enumerate(inputs):
        if not isinstance(binding, Mapping):
            raise EmbeddedInclusiveWeightError(f"input binding {index} is not an object")
        canonical_inputs.append(
            _validate_file_binding(
                binding,
                f"input binding {index}",
                verify_file=verify_input_files,
            )
        )
    if len({item["path"] for item in canonical_inputs}) != len(canonical_inputs):
        raise EmbeddedInclusiveWeightError("analysis-weight input bindings contain duplicates")

    accepted = payload.get("accepted_row_count")
    if type(accepted) is not int or accepted <= 0:
        raise EmbeddedInclusiveWeightError("analysis-weight accepted_row_count is invalid")
    totals: dict[str, list[float]] = {}
    for moment in ("sumw", "sumw2"):
        raw_values = payload.get(moment)
        if not isinstance(raw_values, list) or len(raw_values) != nbins:
            raise EmbeddedInclusiveWeightError(f"analysis-weight {moment} shape differs")
        values = [_finite(value, f"{moment} value") for value in raw_values]
        if any(value < 0.0 for value in values):
            raise EmbeddedInclusiveWeightError(f"analysis-weight {moment} is negative")
        totals[moment] = values

    sources = payload.get("sources")
    if not isinstance(sources, Mapping) or set(sources) != set(SAMPLE_IDS):
        raise EmbeddedInclusiveWeightError("analysis-weight source inventory differs")
    recomputed_sumw = [0.0] * nbins
    recomputed_sumw2 = [0.0] * nbins
    recomputed_count = 0
    parsed_sources: dict[str, tuple[list[float], list[float], list[int], list[Any]]] = {}
    for sample_id in SAMPLE_IDS:
        source = sources.get(sample_id)
        if not isinstance(source, Mapping) or set(source) != {
            "sumw", "sumw2", "raw_count", "fraction"
        }:
            raise EmbeddedInclusiveWeightError(f"{sample_id} source fields differ")
        sample_moments: dict[str, list[float]] = {}
        for moment in ("sumw", "sumw2"):
            raw_values = source.get(moment)
            if not isinstance(raw_values, list) or len(raw_values) != nbins:
                raise EmbeddedInclusiveWeightError(f"{sample_id}.{moment} shape differs")
            values = [_finite(value, f"{sample_id}.{moment}") for value in raw_values]
            if any(value < 0.0 for value in values):
                raise EmbeddedInclusiveWeightError(f"{sample_id}.{moment} is negative")
            sample_moments[moment] = values
        raw_counts = source.get("raw_count")
        fractions = source.get("fraction")
        if (
            not isinstance(raw_counts, list)
            or len(raw_counts) != nbins
            or any(type(value) is not int or value < 0 for value in raw_counts)
            or not isinstance(fractions, list)
            or len(fractions) != nbins
        ):
            raise EmbeddedInclusiveWeightError(f"{sample_id} count/fraction shape differs")
        recomputed_count += sum(raw_counts)
        recomputed_sumw = [
            total + value for total, value in zip(recomputed_sumw, sample_moments["sumw"])
        ]
        recomputed_sumw2 = [
            total + value for total, value in zip(recomputed_sumw2, sample_moments["sumw2"])
        ]
        parsed_sources[sample_id] = (
            sample_moments["sumw"],
            sample_moments["sumw2"],
            list(raw_counts),
            list(fractions),
        )
    if recomputed_count != accepted:
        raise EmbeddedInclusiveWeightError("analysis-weight accepted_row_count differs from sources")
    for label, actual, expected in (
        ("sumw", totals["sumw"], recomputed_sumw),
        ("sumw2", totals["sumw2"], recomputed_sumw2),
    ):
        if any(
            not math.isclose(left, right, rel_tol=2e-13, abs_tol=1e-18)
            for left, right in zip(actual, expected)
        ):
            raise EmbeddedInclusiveWeightError(f"analysis-weight {label} totals differ")
    for sample_id, (sample_sumw, _, _, fractions) in parsed_sources.items():
        for index, (value, total, fraction) in enumerate(
            zip(sample_sumw, totals["sumw"], fractions)
        ):
            if total == 0.0:
                if fraction is not None:
                    raise EmbeddedInclusiveWeightError(
                        f"{sample_id}.fraction[{index}] must be null for an empty bin"
                    )
            else:
                parsed_fraction = _finite(fraction, f"{sample_id}.fraction[{index}]")
                if not math.isclose(
                    parsed_fraction, value / total, rel_tol=2e-13, abs_tol=1e-18
                ):
                    raise EmbeddedInclusiveWeightError(
                        f"{sample_id}.fraction[{index}] differs from source/total sumw"
                    )
    if payload.get("weight_provenance") != provider.provenance_payload():
        raise EmbeddedInclusiveWeightError("analysis-weight payload provenance differs")


def write_analysis_weight_artifacts(
    payload_path: Path,
    receipt_path: Path,
    payload: Mapping[str, Any],
    provider: EmbeddedInclusiveAnalysisWeightProvider,
) -> dict[str, Any]:
    """Write one payload and one hash-bound receipt without overwriting either."""
    validate_analysis_weight_payload(payload, provider, verify_input_files=True)
    payload_path = Path(payload_path).resolve()
    receipt_path = Path(receipt_path).resolve()
    _write_exclusive(payload_path, payload)
    receipt = {
        "schema": ANALYSIS_RECEIPT_SCHEMA,
        "status": "pass",
        "analysis_scope": ANALYSIS_SCOPE,
        "sample_family": SIMULATION_FAMILY,
        "sample_ids": list(SAMPLE_IDS),
        "weight_state": ANALYSIS_STATUS,
        "application_order": list(APPLICATION_ORDER),
        "dependency_fingerprint": provider.dependency_fingerprint,
        "payload_path": str(payload_path),
        "payload_sha256": sha256_file(payload_path),
        "payload_size_bytes": payload_path.stat().st_size,
        "payload_schema": ANALYSIS_PAYLOAD_SCHEMA,
        "source_stitch": provider.source.binding_payload(),
        "centrality_contract": provider.provenance_payload()["centrality_contract"],
    }
    try:
        _write_exclusive(receipt_path, receipt)
    except Exception:
        payload_path.unlink(missing_ok=True)
        raise
    return receipt


def load_analysis_weight_receipt(
    path: Path,
    *,
    expected_dependency_fingerprint: str,
    verify_files: bool = True,
) -> dict[str, Any]:
    """Validate a composite receipt intended for plot/correction consumption."""
    receipt_path = Path(path).resolve()
    receipt = _load_object(receipt_path)
    expected_keys = {
        "schema",
        "status",
        "analysis_scope",
        "sample_family",
        "sample_ids",
        "weight_state",
        "application_order",
        "dependency_fingerprint",
        "payload_path",
        "payload_sha256",
        "payload_size_bytes",
        "payload_schema",
        "source_stitch",
        "centrality_contract",
    }
    if (
        set(receipt) != expected_keys
        or receipt.get("schema") != ANALYSIS_RECEIPT_SCHEMA
        or receipt.get("status") != "pass"
        or receipt.get("analysis_scope") != ANALYSIS_SCOPE
        or receipt.get("sample_family") != SIMULATION_FAMILY
        or tuple(receipt.get("sample_ids", ())) != SAMPLE_IDS
        or receipt.get("weight_state") != ANALYSIS_STATUS
        or receipt.get("application_order") != list(APPLICATION_ORDER)
        or receipt.get("payload_schema") != ANALYSIS_PAYLOAD_SCHEMA
    ):
        raise EmbeddedInclusiveWeightError("analysis-weight receipt contract differs")
    dependency = receipt.get("dependency_fingerprint")
    if not _is_sha256(expected_dependency_fingerprint) or dependency != expected_dependency_fingerprint:
        raise EmbeddedInclusiveWeightError("analysis-weight dependency fingerprint differs")
    payload_path = Path(str(receipt.get("payload_path"))).resolve()
    if verify_files:
        if (
            not payload_path.is_file()
            or sha256_file(payload_path) != receipt.get("payload_sha256")
            or payload_path.stat().st_size != receipt.get("payload_size_bytes")
        ):
            raise EmbeddedInclusiveWeightError("analysis-weight payload is missing or stale")
    if not _is_sha256(receipt.get("payload_sha256")) or type(
        receipt.get("payload_size_bytes")
    ) is not int:
        raise EmbeddedInclusiveWeightError("analysis-weight payload binding is invalid")
    for role in ("source_stitch", "centrality_contract"):
        if not isinstance(receipt.get(role), Mapping):
            raise EmbeddedInclusiveWeightError(f"analysis-weight receipt lacks {role}")
    source_binding = receipt["source_stitch"]
    centrality_binding = receipt["centrality_contract"]
    assert isinstance(source_binding, Mapping) and isinstance(centrality_binding, Mapping)
    source = load_source_stitch_artifact(
        Path(str(source_binding.get("assembly_path"))),
        Path(str(source_binding.get("validation_receipt_path"))),
        verify_authority_file=verify_files,
    )
    if source.binding_payload() != dict(source_binding):
        raise EmbeddedInclusiveWeightError("analysis-weight source-stitch binding differs")
    centrality_path = Path(str(centrality_binding.get("path"))).resolve()
    if (
        not centrality_path.is_file()
        or sha256_file(centrality_path) != centrality_binding.get("sha256")
    ):
        raise EmbeddedInclusiveWeightError("analysis-weight centrality receipt is missing or stale")
    centrality = load_contract_receipt(
        centrality_path,
        expected_dependency_fingerprint=str(
            centrality_binding.get("dependency_fingerprint", "")
        ),
        require_ready=True,
        verify_dependency_files=verify_files,
    )
    if (
        centrality.contract_fingerprint
        != centrality_binding.get("contract_fingerprint")
        or centrality_binding.get("sample_family") != SAMPLE_FAMILY
    ):
        raise EmbeddedInclusiveWeightError("analysis-weight centrality contract binding differs")
    recomputed_dependency = _sha256_payload({
        "source_stitch": source.binding_payload(),
        "centrality_receipt_path": str(centrality_path),
        "centrality_receipt_sha256": sha256_file(centrality_path),
        "centrality_dependency_fingerprint": centrality.dependency_fingerprint,
        "centrality_contract_fingerprint": centrality.contract_fingerprint,
        "application_order": list(APPLICATION_ORDER),
    })
    if dependency != recomputed_dependency:
        raise EmbeddedInclusiveWeightError("analysis-weight dependency fingerprint is internally inconsistent")
    provider = EmbeddedInclusiveAnalysisWeightProvider(
        source=source,
        centrality=centrality,
        centrality_receipt_path=centrality_path,
        centrality_receipt_sha256=sha256_file(centrality_path),
    )
    if provider.dependency_fingerprint != dependency:
        raise EmbeddedInclusiveWeightError("analysis-weight provider reconstruction differs")
    if receipt.get("source_stitch") != provider.source.binding_payload() or receipt.get(
        "centrality_contract"
    ) != provider.provenance_payload()["centrality_contract"]:
        raise EmbeddedInclusiveWeightError("analysis-weight receipt provenance differs")
    if verify_files:
        validate_analysis_weight_payload(
            _load_object(payload_path), provider, verify_input_files=True
        )
    return receipt


def _bound_file(path: Path) -> dict[str, Any]:
    resolved = Path(path).resolve()
    if not resolved.is_file() or resolved.stat().st_size <= 0:
        raise EmbeddedInclusiveWeightError(f"downstream artifact is missing or empty: {resolved}")
    return {
        "path": str(resolved),
        "sha256": sha256_file(resolved),
        "size_bytes": resolved.stat().st_size,
    }


def _bound_artifact(path: Path, role: str) -> dict[str, Any]:
    if not isinstance(role, str) or not role or not role.replace("_", "").isalnum():
        raise EmbeddedInclusiveWeightError(f"downstream artifact role is invalid: {role!r}")
    return {"role": role, **_bound_file(path)}


def write_downstream_artifact_receipt(
    receipt_path: Path,
    *,
    analysis_weight_receipt_path: Path,
    expected_dependency_fingerprint: str,
    audit_path: Path,
    artifacts: Mapping[str, Path],
    input_artifacts: Sequence[Path],
) -> dict[str, Any]:
    """Bind nominal downstream bytes to the complete analysis-weight chain."""
    destination = Path(receipt_path).resolve()
    if not artifacts or not input_artifacts:
        raise EmbeddedInclusiveWeightError(
            "downstream receipt requires output artifacts and exact input artifacts"
        )
    analysis_path = Path(analysis_weight_receipt_path).resolve()
    analysis_receipt = load_analysis_weight_receipt(
        analysis_path,
        expected_dependency_fingerprint=expected_dependency_fingerprint,
        verify_files=True,
    )
    audit = Path(audit_path).resolve()
    audit_payload = _load_object(audit)
    validated_generator = _validate_downstream_audit(
        audit_payload,
        downstream_receipt_path=destination,
        analysis_weight_receipt_path=analysis_path,
        analysis_weight_receipt=analysis_receipt,
        expected_dependency_fingerprint=expected_dependency_fingerprint,
        verify_generator_file=True,
    )
    artifact_bindings = [_bound_artifact(path, role) for role, path in artifacts.items()]
    if len({item["role"] for item in artifact_bindings}) != len(artifact_bindings) or len(
        {item["path"] for item in artifact_bindings}
    ) != len(artifact_bindings):
        raise EmbeddedInclusiveWeightError("downstream artifact roles or paths are duplicated")
    artifact_bindings.sort(key=lambda item: item["role"])
    input_bindings = [_bound_file(path) for path in input_artifacts]
    input_bindings.sort(key=lambda item: item["path"])
    if len({item["path"] for item in input_bindings}) != len(input_bindings):
        raise EmbeddedInclusiveWeightError("downstream input artifact paths are duplicated")
    payload = {
        "schema": DOWNSTREAM_RECEIPT_SCHEMA,
        "status": "pass",
        "downstream_state": DOWNSTREAM_STATUS,
        "analysis_scope": ANALYSIS_SCOPE,
        "sample_family": SIMULATION_FAMILY,
        "sample_ids": list(SAMPLE_IDS),
        "weight_state": ANALYSIS_STATUS,
        "application_order": list(APPLICATION_ORDER),
        "dependency_fingerprint": expected_dependency_fingerprint,
        "analysis_weight_receipt": {
            "path": str(analysis_path),
            "sha256": sha256_file(analysis_path),
            "size_bytes": analysis_path.stat().st_size,
            "payload_sha256": analysis_receipt["payload_sha256"],
        },
        "audit": {
            "path": str(audit),
            "sha256": sha256_file(audit),
            "size_bytes": audit.stat().st_size,
            "schema": audit_payload.get("schema"),
        },
        "generator": validated_generator,
        "inputs": input_bindings,
        "artifacts": artifact_bindings,
    }
    _write_exclusive(destination, payload)
    return payload


def load_downstream_artifact_receipt(
    path: Path,
    *,
    expected_dependency_fingerprint: str,
    expected_artifacts: Mapping[str, Path] | None = None,
    expected_inputs: Sequence[Path] | None = None,
    expected_audit_path: Path | None = None,
    verify_files: bool = True,
) -> dict[str, Any]:
    """Validate a nominal output receipt and every transitive file binding."""
    receipt_path = Path(path).resolve()
    receipt = _load_object(receipt_path)
    expected_keys = {
        "schema",
        "status",
        "downstream_state",
        "analysis_scope",
        "sample_family",
        "sample_ids",
        "weight_state",
        "application_order",
        "dependency_fingerprint",
        "analysis_weight_receipt",
        "audit",
        "generator",
        "inputs",
        "artifacts",
    }
    if (
        set(receipt) != expected_keys
        or receipt.get("schema") != DOWNSTREAM_RECEIPT_SCHEMA
        or receipt.get("status") != "pass"
        or receipt.get("downstream_state") != DOWNSTREAM_STATUS
        or receipt.get("analysis_scope") != ANALYSIS_SCOPE
        or receipt.get("sample_family") != SIMULATION_FAMILY
        or tuple(receipt.get("sample_ids", ())) != SAMPLE_IDS
        or receipt.get("weight_state") != ANALYSIS_STATUS
        or receipt.get("application_order") != list(APPLICATION_ORDER)
        or receipt.get("dependency_fingerprint") != expected_dependency_fingerprint
        or not _is_sha256(expected_dependency_fingerprint)
    ):
        raise EmbeddedInclusiveWeightError("downstream artifact receipt contract differs")

    analysis_binding = receipt.get("analysis_weight_receipt")
    if not isinstance(analysis_binding, Mapping) or set(analysis_binding) != {
        "path", "sha256", "size_bytes", "payload_sha256"
    }:
        raise EmbeddedInclusiveWeightError("downstream analysis-weight binding differs")
    analysis_file_binding = _validate_file_binding(
        {key: analysis_binding.get(key) for key in ("path", "sha256", "size_bytes")},
        "downstream analysis-weight receipt",
        verify_file=verify_files,
    )
    analysis_receipt = load_analysis_weight_receipt(
        Path(analysis_file_binding["path"]),
        expected_dependency_fingerprint=expected_dependency_fingerprint,
        verify_files=verify_files,
    )
    if analysis_binding.get("payload_sha256") != analysis_receipt.get("payload_sha256"):
        raise EmbeddedInclusiveWeightError("downstream payload binding differs")

    audit_binding = receipt.get("audit")
    if not isinstance(audit_binding, Mapping) or set(audit_binding) != {
        "path", "sha256", "size_bytes", "schema"
    }:
        raise EmbeddedInclusiveWeightError("downstream audit binding differs")
    validated_audit_binding = _validate_file_binding(
        {key: audit_binding.get(key) for key in ("path", "sha256", "size_bytes")},
        "downstream audit",
        verify_file=verify_files,
    )
    audit_path = Path(validated_audit_binding["path"])
    if expected_audit_path is not None and audit_path != Path(expected_audit_path).resolve():
        raise EmbeddedInclusiveWeightError("downstream audit path differs from consumer request")
    if verify_files:
        audit_payload = _load_object(audit_path)
        if audit_payload.get("schema") != audit_binding.get("schema"):
            raise EmbeddedInclusiveWeightError("downstream audit contents differ")
        validated_audit_generator = _validate_downstream_audit(
            audit_payload,
            downstream_receipt_path=receipt_path,
            analysis_weight_receipt_path=Path(analysis_file_binding["path"]),
            analysis_weight_receipt=analysis_receipt,
            expected_dependency_fingerprint=expected_dependency_fingerprint,
            verify_generator_file=True,
        )
    else:
        validated_audit_generator = None

    generator_binding = receipt.get("generator")
    if not isinstance(generator_binding, Mapping):
        raise EmbeddedInclusiveWeightError("downstream generator binding differs")
    validated_generator = _validate_file_binding(
        generator_binding, "downstream generator", verify_file=verify_files
    )
    if verify_files and validated_audit_generator != validated_generator:
        raise EmbeddedInclusiveWeightError("downstream audit generator binding differs")

    inputs = receipt.get("inputs")
    if not isinstance(inputs, list) or not inputs:
        raise EmbeddedInclusiveWeightError("downstream input artifact inventory is empty")
    observed_inputs: list[Path] = []
    for index, binding in enumerate(inputs):
        if not isinstance(binding, Mapping):
            raise EmbeddedInclusiveWeightError(
                f"downstream input artifact binding {index} is not an object"
            )
        validated = _validate_file_binding(
            binding,
            f"downstream input artifact {index}",
            verify_file=verify_files,
        )
        input_path = Path(validated["path"])
        if input_path in observed_inputs:
            raise EmbeddedInclusiveWeightError("downstream input artifact paths are duplicated")
        observed_inputs.append(input_path)
    if observed_inputs != sorted(observed_inputs):
        raise EmbeddedInclusiveWeightError("downstream input artifacts are not canonically sorted")
    if expected_inputs is not None and observed_inputs != sorted(
        Path(value).resolve() for value in expected_inputs
    ):
        raise EmbeddedInclusiveWeightError("downstream input inventory differs from consumer request")

    artifacts = receipt.get("artifacts")
    if not isinstance(artifacts, list) or not artifacts:
        raise EmbeddedInclusiveWeightError("downstream artifact inventory is empty")
    observed: dict[str, Path] = {}
    for index, binding in enumerate(artifacts):
        if not isinstance(binding, Mapping) or set(binding) != {
            "role", "path", "sha256", "size_bytes"
        }:
            raise EmbeddedInclusiveWeightError(f"downstream artifact binding {index} differs")
        role = binding.get("role")
        if not isinstance(role, str) or role in observed:
            raise EmbeddedInclusiveWeightError("downstream artifact roles are invalid or duplicated")
        validated = _validate_file_binding(
            {key: binding.get(key) for key in ("path", "sha256", "size_bytes")},
            f"downstream artifact {role}",
            verify_file=verify_files,
        )
        artifact_path = Path(validated["path"])
        if artifact_path in observed.values():
            raise EmbeddedInclusiveWeightError("downstream artifact paths are duplicated")
        observed[role] = artifact_path
    if expected_artifacts is not None:
        expected = {role: Path(value).resolve() for role, value in expected_artifacts.items()}
        if observed != expected:
            raise EmbeddedInclusiveWeightError("downstream artifact inventory differs from consumer request")
    return receipt


__all__ = [
    "ANALYSIS_SCOPE",
    "ANALYSIS_PAYLOAD_SCHEMA",
    "ANALYSIS_RECEIPT_SCHEMA",
    "ANALYSIS_STATUS",
    "DOWNSTREAM_RECEIPT_SCHEMA",
    "DOWNSTREAM_STATUS",
    "EmbeddedInclusiveAnalysisWeightProvider",
    "EmbeddedInclusiveWeightError",
    "SAMPLE_IDS",
    "SIMULATION_FAMILY",
    "SourceStitchArtifact",
    "build_reco_cluster_source_fraction_payload",
    "load_analysis_weight_receipt",
    "load_downstream_artifact_receipt",
    "load_source_stitch_artifact",
    "write_analysis_weight_artifacts",
    "write_downstream_artifact_receipt",
    "validate_analysis_weight_payload",
]
