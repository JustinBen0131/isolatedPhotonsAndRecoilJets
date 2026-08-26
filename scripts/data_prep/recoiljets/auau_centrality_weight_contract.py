#!/usr/bin/env python3
"""Current-artifact diagnostics and authority-gated Au+Au SIM centrality weights.

The numerical maps in this module are deliberately *not* constants.  A
diagnostic contract can be resolved from current Au+Au DATA, photon+jet, and
inclusive-jet artifact pointers, with every pointer and ROOT input SHA-256
bound.  Current histogram agreement is not nominal authority: only
``build_canonical_contract`` can emit an applicable contract, and it requires
the complete hash-bound authority set from the source audit.

The canonical application order is::

    existing event weight -> canonical source stitch -> centrality reweight

Plotting code should consume an artifact stamped with the resulting contract;
it must not derive or apply these weights locally.
"""

from __future__ import annotations

from bisect import bisect_right
from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, Mapping, Optional, Sequence, Tuple


REPO = Path(__file__).resolve().parents[3]
CURRENT_ARTIFACT_DIR = REPO / "dataOutput/current_recoiljets_artifacts/current"

CONTRACT_SCHEMA = "AuAuCentralityWeightContractV2"
DIAGNOSTIC_CONTRACT_STATUS = "PASS_DIAGNOSTIC_CURRENT_HISTOGRAM_BINDING__NOT_NOMINAL"
CANONICAL_CONTRACT_STATUS = "READY_FOR_CANONICAL_AUAU_SIM_APPLICATION"
# Backward-compatible name for callers that only inspect the current-histogram
# builder.  That builder is intentionally diagnostic after the source audit.
CONTRACT_STATUS = DIAGNOSTIC_CONTRACT_STATUS
FORMULA = "w_cent(c_i) = p_DATA(c_i) / p_SIM_family(c_i)"
SUPPORT_PERCENT = (0.0, 80.0)
BIN_WIDTH_PERCENT = 5.0
BIN_EDGES_PERCENT = tuple(float(value) for value in range(0, 81, 5))

DATA_SAMPLE_KEY = "auau_data_merged"
DATA_HISTOGRAM_KEY = "MBD_NS_geq_2_vtx_lt_150/h_centrality"
FAMILY_SPECS: Mapping[str, Tuple[str, str]] = {
    "photonjet": ("auau_sim_photonjet_merged", "SIM/h_centrality"),
    "inclusivejet": ("auau_sim_inclusivejet_merged", "SIM/h_centrality"),
}

EXISTING_EVENT_WEIGHT_COMPONENT = "existing_event_weight"
CANONICAL_STITCH_COMPONENT = "canonical_source_stitch"
CENTRALITY_COMPONENT_PREFIX = "auau_centrality_reweight:"
APPLICATION_ORDER = (
    EXISTING_EVENT_WEIGHT_COMPONENT,
    CANONICAL_STITCH_COMPONENT,
    "auau_centrality_reweight",
)
REQUIRED_AUTHORITY_ROLES = (
    "data_target",
    "sim_stitch",
    "centrality_calibration",
    "tree_selection",
    "derivation_code",
)

HistogramLoader = Callable[[Path, str], Tuple[Sequence[float], Sequence[float]]]


class CentralityWeightContractError(ValueError):
    """Base class for fail-closed contract errors."""


class CurrentArtifactResolutionError(CentralityWeightContractError):
    """A required current artifact pointer or ROOT dependency is invalid."""


class ContractBuildError(CentralityWeightContractError):
    """Current histograms cannot form a valid nominal centrality contract."""


class UnsupportedSampleFamilyError(CentralityWeightContractError):
    """The requested sample is not an embedded Au+Au SIM family."""


class UnsupportedCentralityError(CentralityWeightContractError):
    """The requested centrality is outside the frozen half-open support."""


class WeightApplicationOrderError(CentralityWeightContractError):
    """Centrality reweighting was requested at the wrong weight stage."""


class DoubleApplicationError(CentralityWeightContractError):
    """A centrality component is already present in the weight provenance."""


class ContractNotReadyError(CentralityWeightContractError):
    """A diagnostic or incomplete contract was used for nominal application."""


class AuthorityValidationError(CentralityWeightContractError):
    """Canonical authority metadata is missing, malformed, or unbound."""


class StaleDependencyError(CentralityWeightContractError):
    """The consumer and contract dependency fingerprints do not agree."""


def sha256_file(path: Path) -> str:
    """Return the SHA-256 digest of one exact file."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _sha256_payload(payload: Mapping[str, Any]) -> str:
    encoded = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=True,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _float_hex(values: Iterable[float]) -> Tuple[str, ...]:
    return tuple(float(value).hex() for value in values)


def _is_sha256(value: Any) -> bool:
    if not isinstance(value, str) or len(value) != 64:
        return False
    try:
        int(value, 16)
    except ValueError:
        return False
    return value == value.lower()


def _canonical_json(payload: Mapping[str, Any]) -> str:
    try:
        return json.dumps(
            payload,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=True,
            allow_nan=False,
        )
    except (TypeError, ValueError) as error:
        raise AuthorityValidationError("authority metadata is not canonical JSON") from error


@dataclass(frozen=True)
class RootDependency:
    path: str
    sha256: str
    size_bytes: int

    def payload(self) -> Dict[str, Any]:
        return {
            "path": self.path,
            "sha256": self.sha256,
            "size_bytes": self.size_bytes,
        }


@dataclass(frozen=True)
class ArtifactDependency:
    sample_key: str
    histogram_key: str
    pointer_path: str
    pointer_sha256: str
    current_entry_id: str
    campaign_tag: str
    canonical_status: str
    roots: Tuple[RootDependency, ...]

    def payload(self) -> Dict[str, Any]:
        return {
            "sample_key": self.sample_key,
            "histogram_key": self.histogram_key,
            "pointer_path": self.pointer_path,
            "pointer_sha256": self.pointer_sha256,
            "current_entry_id": self.current_entry_id,
            "campaign_tag": self.campaign_tag,
            "canonical_status": self.canonical_status,
            "roots": [root.payload() for root in self.roots],
        }


@dataclass(frozen=True)
class AuthorityBinding:
    """One explicit canonical authority with an exact SHA-256 identity."""

    role: str
    authority_id: str
    sha256: str
    metadata_json: str

    def payload(self) -> Dict[str, Any]:
        payload = json.loads(self.metadata_json)
        payload.update(
            {
                "role": self.role,
                "authority_id": self.authority_id,
                "sha256": self.sha256,
            }
        )
        return payload


@dataclass(frozen=True)
class DistributionInput:
    """Exact counts or probabilities supplied to the canonical builder."""

    role: str
    kind: str
    values: Tuple[float, ...]

    def payload(self) -> Dict[str, Any]:
        return {
            "role": self.role,
            "kind": self.kind,
            "values": list(self.values),
        }

    def fingerprint_payload(self) -> Dict[str, Any]:
        return {
            "role": self.role,
            "kind": self.kind,
            "values_hex": _float_hex(self.values),
        }


@dataclass(frozen=True)
class FamilyWeightMap:
    family: str
    sample_key: str
    simulation_probability: Tuple[float, ...]
    weights: Tuple[float, ...]

    def payload(self) -> Dict[str, Any]:
        return {
            "family": self.family,
            "sample_key": self.sample_key,
            "simulation_probability": list(self.simulation_probability),
            "weights": list(self.weights),
        }


@dataclass(frozen=True)
class EventWeightState:
    """A numeric event weight plus its ordered, machine-checkable components."""

    value: float
    components: Tuple[str, ...]

    def __post_init__(self) -> None:
        if not math.isfinite(float(self.value)) or float(self.value) < 0.0:
            raise CentralityWeightContractError("event weight must be finite and nonnegative")
        if not self.components or any(not isinstance(item, str) or not item for item in self.components):
            raise CentralityWeightContractError("weight provenance components must be nonempty strings")
        if len(set(self.components)) != len(self.components):
            raise CentralityWeightContractError("weight provenance contains a duplicate component")

    @classmethod
    def after_canonical_stitch(
        cls,
        value: float,
        prior_components: Sequence[str] = (EXISTING_EVENT_WEIGHT_COMPONENT,),
    ) -> "EventWeightState":
        """Construct the only state accepted by centrality ``apply_once``.

        Callers may name earlier event-level components, but the canonical
        source stitch must be appended here and must remain immediately before
        centrality weighting.
        """
        prior = tuple(prior_components)
        if CANONICAL_STITCH_COMPONENT in prior:
            raise WeightApplicationOrderError("canonical source stitch is already present")
        if any(item.startswith(CENTRALITY_COMPONENT_PREFIX) for item in prior):
            raise DoubleApplicationError("centrality reweight component is already present")
        return cls(float(value), prior + (CANONICAL_STITCH_COMPONENT,))


@dataclass(frozen=True)
class CentralityWeightContract:
    schema: str
    status: str
    dependency_fingerprint: str
    contract_fingerprint: str
    dependencies: Tuple[ArtifactDependency, ...]
    data_probability: Tuple[float, ...]
    family_maps: Tuple[FamilyWeightMap, ...]
    authorities: Tuple[AuthorityBinding, ...] = ()
    distribution_inputs: Tuple[DistributionInput, ...] = ()

    def family_map(self, sample_family: str) -> FamilyWeightMap:
        for family_map in self.family_maps:
            if family_map.family == sample_family:
                return family_map
        raise UnsupportedSampleFamilyError(
            "unsupported Au+Au embedded SIM family {!r}; expected one of {}".format(
                sample_family,
                sorted(item.family for item in self.family_maps),
            )
        )

    def assert_ready_for_application(self, expected_dependency_fingerprint: str) -> None:
        if self.status != CANONICAL_CONTRACT_STATUS:
            raise ContractNotReadyError(
                "contract status {!r} is not ready for canonical Au+Au SIM application".format(
                    self.status
                )
            )
        roles = {authority.role for authority in self.authorities}
        missing = sorted(set(REQUIRED_AUTHORITY_ROLES) - roles)
        if missing:
            raise AuthorityValidationError(
                "canonical contract is missing authority roles: {}".format(missing)
            )
        if not _is_sha256(expected_dependency_fingerprint):
            raise StaleDependencyError(
                "application requires an exact expected dependency SHA-256 fingerprint"
            )
        if expected_dependency_fingerprint != self.dependency_fingerprint:
            raise StaleDependencyError(
                "contract dependency fingerprint is stale: expected {}, contract {}".format(
                    expected_dependency_fingerprint,
                    self.dependency_fingerprint,
                )
            )

    def bin_index(self, centrality_percent: float) -> int:
        try:
            centrality = float(centrality_percent)
        except (TypeError, ValueError) as error:
            raise UnsupportedCentralityError(
                "centrality must be a finite number in [0,80) percent"
            ) from error
        if not math.isfinite(centrality) or not (
            SUPPORT_PERCENT[0] <= centrality < SUPPORT_PERCENT[1]
        ):
            raise UnsupportedCentralityError(
                "centrality {!r} is outside the half-open [0,80) percent support".format(
                    centrality_percent
                )
            )
        index = bisect_right(BIN_EDGES_PERCENT, centrality) - 1
        if not 0 <= index < len(BIN_EDGES_PERCENT) - 1:
            raise UnsupportedCentralityError(
                "centrality {!r} did not resolve to a supported bin".format(centrality_percent)
            )
        return index

    def weight(self, centrality_percent: float, sample_family: str) -> float:
        family_map = self.family_map(sample_family)
        value = family_map.weights[self.bin_index(centrality_percent)]
        if not math.isfinite(value) or value <= 0.0:
            raise ContractBuildError(
                "resolved contract contains an invalid weight for family={} centrality={}".format(
                    sample_family,
                    centrality_percent,
                )
            )
        return value

    def apply_once(
        self,
        state: EventWeightState,
        centrality_percent: float,
        sample_family: str,
        *,
        expected_dependency_fingerprint: str,
    ) -> EventWeightState:
        """Append the family map exactly once after canonical source stitching."""
        self.assert_ready_for_application(expected_dependency_fingerprint)
        if any(item.startswith(CENTRALITY_COMPONENT_PREFIX) for item in state.components):
            raise DoubleApplicationError("centrality reweight component is already present")
        if state.components[-1] != CANONICAL_STITCH_COMPONENT:
            raise WeightApplicationOrderError(
                "centrality reweighting must immediately follow canonical source stitching"
            )
        factor = self.weight(centrality_percent, sample_family)
        result = float(state.value) * factor
        if not math.isfinite(result) or result < 0.0:
            raise CentralityWeightContractError("complete event weight is invalid")
        component = "{}{}:{}".format(
            CENTRALITY_COMPONENT_PREFIX,
            sample_family,
            self.contract_fingerprint,
        )
        return EventWeightState(result, state.components + (component,))

    def to_payload(self) -> Dict[str, Any]:
        """Return a deterministic, serializable current-resolution receipt."""
        return {
            "schema": self.schema,
            "status": self.status,
            "formula": FORMULA,
            "scope": "embedded Au+Au simulation only",
            "numerical_lifetime": (
                "canonical authority dependency fingerprint"
                if self.status == CANONICAL_CONTRACT_STATUS
                else "current artifact dependency fingerprint"
            ),
            "support_percent": list(SUPPORT_PERCENT),
            "bin_width_percent": BIN_WIDTH_PERCENT,
            "bin_edges_percent": list(BIN_EDGES_PERCENT),
            "lower_edge_inclusive": True,
            "upper_edge_exclusive": True,
            "application_order": list(APPLICATION_ORDER),
            "dependency_fingerprint": self.dependency_fingerprint,
            "contract_fingerprint": self.contract_fingerprint,
            "dependencies": [dependency.payload() for dependency in self.dependencies],
            "authorities": [
                authority.payload()
                for authority in sorted(self.authorities, key=lambda item: item.role)
            ],
            "distribution_inputs": [
                distribution.payload()
                for distribution in sorted(self.distribution_inputs, key=lambda item: item.role)
            ],
            "data_probability": list(self.data_probability),
            "family_maps": {
                family_map.family: family_map.payload()
                for family_map in sorted(self.family_maps, key=lambda item: item.family)
            },
        }

    def write_receipt(self, path: Path) -> None:
        """Write this resolved contract atomically when a consumer materializes it."""
        destination = Path(path)
        destination.parent.mkdir(parents=True, exist_ok=True)
        temporary = destination.with_suffix(destination.suffix + ".tmp")
        temporary.write_text(
            json.dumps(self.to_payload(), indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        temporary.replace(destination)


@dataclass(frozen=True)
class CanonicalContractResolution:
    """Result of resolving the one current canonical numerical map.

    ``action`` is deliberately machine-readable so an upstream Tree-to-hist
    stage can distinguish an exact reuse from an input-driven rebuild without
    reconstructing freshness logic itself.
    """

    contract: CentralityWeightContract
    action: str


def _load_pointer(pointer_path: Path, sample_key: str, histogram_key: str) -> ArtifactDependency:
    if not pointer_path.is_file():
        raise CurrentArtifactResolutionError(
            "missing current artifact pointer: {}".format(pointer_path)
        )
    try:
        payload = json.loads(pointer_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise CurrentArtifactResolutionError(
            "invalid current artifact pointer: {}".format(pointer_path)
        ) from error
    if not isinstance(payload, dict):
        raise CurrentArtifactResolutionError("current artifact pointer is not a JSON object")
    if payload.get("sample_key") != sample_key:
        raise CurrentArtifactResolutionError(
            "pointer sample_key mismatch: expected {}, found {!r}".format(
                sample_key,
                payload.get("sample_key"),
            )
        )
    if payload.get("canonical_status") != "canonical":
        raise CurrentArtifactResolutionError(
            "current artifact {} is not canonical: {!r}".format(
                sample_key,
                payload.get("canonical_status"),
            )
        )
    current_entry_id = payload.get("current_entry_id")
    if not isinstance(current_entry_id, str) or not current_entry_id:
        raise CurrentArtifactResolutionError(
            "current artifact {} has no current_entry_id".format(sample_key)
        )
    root_paths = payload.get("root_paths")
    if not isinstance(root_paths, list) or not root_paths:
        raise CurrentArtifactResolutionError(
            "current artifact {} has no ROOT paths".format(sample_key)
        )

    roots = []
    for raw_path in root_paths:
        if not isinstance(raw_path, str) or not raw_path:
            raise CurrentArtifactResolutionError(
                "current artifact {} contains an invalid ROOT path".format(sample_key)
            )
        root_path = Path(raw_path).expanduser()
        if not root_path.is_absolute():
            root_path = (REPO / root_path).resolve()
        if not root_path.is_file():
            raise CurrentArtifactResolutionError(
                "current artifact ROOT does not exist: {}".format(root_path)
            )
        roots.append(
            RootDependency(
                path=str(root_path),
                sha256=sha256_file(root_path),
                size_bytes=root_path.stat().st_size,
            )
        )
    return ArtifactDependency(
        sample_key=sample_key,
        histogram_key=histogram_key,
        pointer_path=str(pointer_path.resolve()),
        pointer_sha256=sha256_file(pointer_path),
        current_entry_id=current_entry_id,
        campaign_tag=str(payload.get("campaign_tag") or ""),
        canonical_status="canonical",
        roots=tuple(roots),
    )


def resolve_current_inputs(
    current_artifact_dir: Path = CURRENT_ARTIFACT_DIR,
) -> Tuple[ArtifactDependency, ...]:
    """Resolve and hash-bind the three required current artifact pointers."""
    current_dir = Path(current_artifact_dir)
    specs = [(DATA_SAMPLE_KEY, DATA_HISTOGRAM_KEY)] + [
        FAMILY_SPECS[family] for family in sorted(FAMILY_SPECS)
    ]
    return tuple(
        _load_pointer(current_dir / sample_key / "current.json", sample_key, histogram_key)
        for sample_key, histogram_key in specs
    )


def _default_histogram_loader(path: Path, histogram_key: str) -> Tuple[Sequence[float], Sequence[float]]:
    try:
        import uproot  # type: ignore
    except ImportError as error:
        raise ContractBuildError(
            "uproot is required to resolve current centrality histograms; use "
            "/Users/patsfan753/Desktop/analysis/env/bin/python3"
        ) from error
    try:
        with uproot.open(path) as root_file:
            histogram = root_file[histogram_key]
            values = tuple(float(value) for value in histogram.values(flow=False))
            edges = tuple(float(value) for value in histogram.axis().edges(flow=False))
    except Exception as error:
        raise ContractBuildError(
            "cannot read centrality histogram {}:{}".format(path, histogram_key)
        ) from error
    return values, edges


def _same_edges(left: Sequence[float], right: Sequence[float]) -> bool:
    return len(left) == len(right) and all(
        math.isclose(float(a), float(b), rel_tol=0.0, abs_tol=1.0e-9)
        for a, b in zip(left, right)
    )


def _combined_histogram(
    dependency: ArtifactDependency,
    histogram_loader: HistogramLoader,
) -> Tuple[Tuple[float, ...], Tuple[float, ...]]:
    combined: Optional[list] = None
    reference_edges: Optional[Tuple[float, ...]] = None
    for root in dependency.roots:
        values_raw, edges_raw = histogram_loader(Path(root.path), dependency.histogram_key)
        values = tuple(float(value) for value in values_raw)
        edges = tuple(float(value) for value in edges_raw)
        if len(edges) != len(values) + 1:
            raise ContractBuildError(
                "histogram edge/value mismatch for {}:{}".format(root.path, dependency.histogram_key)
            )
        if reference_edges is None:
            reference_edges = edges
            combined = [0.0 for _ in values]
        elif not _same_edges(reference_edges, edges):
            raise ContractBuildError(
                "ROOT inputs for {} have incompatible centrality binning".format(
                    dependency.sample_key
                )
            )
        assert combined is not None
        if len(combined) != len(values):
            raise ContractBuildError(
                "ROOT inputs for {} have incompatible centrality bin counts".format(
                    dependency.sample_key
                )
            )
        for index, value in enumerate(values):
            if not math.isfinite(value) or value < 0.0:
                raise ContractBuildError(
                    "{} has a nonfinite or negative centrality bin".format(dependency.sample_key)
                )
            combined[index] += value
    if combined is None or reference_edges is None:
        raise ContractBuildError("no histograms resolved for {}".format(dependency.sample_key))
    return tuple(combined), reference_edges


def _rebin_to_contract(
    values: Sequence[float],
    edges: Sequence[float],
    sample_key: str,
) -> Tuple[float, ...]:
    """Sum any exact subdivision of the frozen 5% bins into 0--80%."""
    if len(edges) != len(values) + 1:
        raise ContractBuildError("{} histogram edge/value mismatch".format(sample_key))
    if any(not math.isfinite(float(edge)) for edge in edges):
        raise ContractBuildError("{} histogram has nonfinite edges".format(sample_key))
    if any(float(right) <= float(left) for left, right in zip(edges[:-1], edges[1:])):
        raise ContractBuildError("{} histogram edges are not strictly increasing".format(sample_key))

    rebinned = []
    tolerance = 1.0e-9
    source_bins = [
        (float(left), float(right), float(value))
        for left, right, value in zip(edges[:-1], edges[1:], values)
    ]
    for target_low, target_high in zip(BIN_EDGES_PERCENT[:-1], BIN_EDGES_PERCENT[1:]):
        selected = [
            row
            for row in source_bins
            if row[0] >= target_low - tolerance and row[1] <= target_high + tolerance
        ]
        cursor = target_low
        total = 0.0
        for low, high, value in selected:
            if not math.isclose(low, cursor, rel_tol=0.0, abs_tol=tolerance):
                raise ContractBuildError(
                    "{} does not exactly cover [{}, {})".format(
                        sample_key, target_low, target_high
                    )
                )
            if high > target_high + tolerance:
                raise ContractBuildError(
                    "{} has a source bin crossing a canonical 5% edge".format(sample_key)
                )
            total += value
            cursor = high
        if not math.isclose(cursor, target_high, rel_tol=0.0, abs_tol=tolerance):
            raise ContractBuildError(
                "{} does not exactly cover [{}, {})".format(sample_key, target_low, target_high)
            )
        if not math.isfinite(total) or total <= 0.0:
            raise ContractBuildError(
                "{} has zero or invalid support in [{}, {})".format(
                    sample_key, target_low, target_high
                )
            )
        rebinned.append(total)
    return tuple(rebinned)


def _normalize(counts: Sequence[float], sample_key: str) -> Tuple[float, ...]:
    total = math.fsum(float(value) for value in counts)
    if not math.isfinite(total) or total <= 0.0:
        raise ContractBuildError("{} has a nonpositive 0--80% integral".format(sample_key))
    probabilities = tuple(float(value) / total for value in counts)
    if any(not math.isfinite(value) or value <= 0.0 for value in probabilities):
        raise ContractBuildError("{} has invalid normalized support".format(sample_key))
    return probabilities


def _validated_probability(values: Sequence[float], role: str) -> Tuple[float, ...]:
    probability = tuple(float(value) for value in values)
    if len(probability) != len(BIN_EDGES_PERCENT) - 1:
        raise ContractBuildError(
            "{} probability must have exactly {} bins".format(
                role,
                len(BIN_EDGES_PERCENT) - 1,
            )
        )
    if any(not math.isfinite(value) or value <= 0.0 for value in probability):
        raise ContractBuildError("{} probability has nonpositive or invalid support".format(role))
    if not math.isclose(math.fsum(probability), 1.0, rel_tol=0.0, abs_tol=1.0e-12):
        raise ContractBuildError("{} probability is not normalized to one".format(role))
    return probability


def _make_distribution_input(
    role: str,
    *,
    counts: Optional[Sequence[float]] = None,
    probability: Optional[Sequence[float]] = None,
) -> Tuple[DistributionInput, Tuple[float, ...]]:
    if (counts is None) == (probability is None):
        raise ContractBuildError(
            "{} must supply exactly one of counts or probability".format(role)
        )
    if counts is not None:
        exact = tuple(float(value) for value in counts)
        if len(exact) != len(BIN_EDGES_PERCENT) - 1:
            raise ContractBuildError(
                "{} counts must have exactly {} bins".format(
                    role,
                    len(BIN_EDGES_PERCENT) - 1,
                )
            )
        if any(not math.isfinite(value) or value <= 0.0 for value in exact):
            raise ContractBuildError("{} counts have nonpositive or invalid support".format(role))
        return DistributionInput(role, "counts", exact), _normalize(exact, role)
    assert probability is not None
    exact_probability = _validated_probability(probability, role)
    return (
        DistributionInput(role, "probability", exact_probability),
        exact_probability,
    )


def _authority_bindings(
    authorities: Mapping[str, Mapping[str, Any]],
) -> Tuple[AuthorityBinding, ...]:
    if not isinstance(authorities, Mapping):
        raise AuthorityValidationError("canonical authorities must be a role-keyed mapping")
    missing = sorted(set(REQUIRED_AUTHORITY_ROLES) - set(authorities))
    if missing:
        raise AuthorityValidationError(
            "canonical authorities are missing required roles: {}".format(missing)
        )
    bindings = []
    for role, raw_payload in sorted(authorities.items()):
        if not isinstance(role, str) or not role:
            raise AuthorityValidationError("authority role must be a nonempty string")
        if not isinstance(raw_payload, Mapping):
            raise AuthorityValidationError(
                "authority {} metadata must be a mapping".format(role)
            )
        payload = dict(raw_payload)
        authority_id = payload.pop("authority_id", None)
        sha256 = payload.pop("sha256", None)
        if not isinstance(authority_id, str) or not authority_id:
            raise AuthorityValidationError(
                "authority {} requires a nonempty authority_id".format(role)
            )
        if not _is_sha256(sha256):
            raise AuthorityValidationError(
                "authority {} requires a lowercase SHA-256 digest".format(role)
            )
        bindings.append(
            AuthorityBinding(
                role=role,
                authority_id=authority_id,
                sha256=sha256,
                metadata_json=_canonical_json(payload),
            )
        )
    return tuple(bindings)


def _policy_payload() -> Dict[str, Any]:
    return {
        "schema": CONTRACT_SCHEMA,
        "formula": FORMULA,
        "support_percent": list(SUPPORT_PERCENT),
        "bin_width_percent": BIN_WIDTH_PERCENT,
        "bin_edges_percent": list(BIN_EDGES_PERCENT),
        "lower_edge_inclusive": True,
        "upper_edge_exclusive": True,
        "application_order": list(APPLICATION_ORDER),
        "families": {
            family: {"sample_key": spec[0], "histogram_key": spec[1]}
            for family, spec in sorted(FAMILY_SPECS.items())
        },
    }


def _dependency_payload(dependencies: Sequence[ArtifactDependency]) -> Dict[str, Any]:
    payload = _policy_payload()
    payload["mode"] = "diagnostic_current_histograms"
    payload["dependencies"] = [dependency.payload() for dependency in dependencies]
    return payload


def _canonical_dependency_payload(
    authorities: Sequence[AuthorityBinding],
    distribution_inputs: Sequence[DistributionInput],
) -> Dict[str, Any]:
    payload = _policy_payload()
    payload["mode"] = "canonical_authority_bound"
    payload["authorities"] = [
        authority.payload() for authority in sorted(authorities, key=lambda item: item.role)
    ]
    payload["distribution_inputs"] = [
        distribution.fingerprint_payload()
        for distribution in sorted(distribution_inputs, key=lambda item: item.role)
    ]
    return payload


def _contract_fingerprint(
    status: str,
    dependency_fingerprint: str,
    data_probability: Sequence[float],
    family_maps: Sequence[FamilyWeightMap],
) -> str:
    numerical_payload = {
        "schema": CONTRACT_SCHEMA,
        "status": status,
        "dependency_fingerprint": dependency_fingerprint,
        "data_probability_hex": _float_hex(data_probability),
        "family_maps": {
            item.family: {
                "simulation_probability_hex": _float_hex(item.simulation_probability),
                "weights_hex": _float_hex(item.weights),
            }
            for item in sorted(family_maps, key=lambda value: value.family)
        },
    }
    return _sha256_payload(numerical_payload)


def build_current_contract(
    current_artifact_dir: Path = CURRENT_ARTIFACT_DIR,
    histogram_loader: Optional[HistogramLoader] = None,
) -> CentralityWeightContract:
    """Resolve current histograms into a diagnostic, non-applicable contract.

    No numerical map is cached or embedded in source.  Callers that materialize
    downstream products should persist :meth:`CentralityWeightContract.to_payload`
    beside them.  Source-histogram agreement alone is not nominal authority, so
    :meth:`CentralityWeightContract.apply_once` always refuses this status.
    """
    dependencies = resolve_current_inputs(current_artifact_dir)
    loader = histogram_loader or _default_histogram_loader

    probabilities: Dict[str, Tuple[float, ...]] = {}
    for dependency in dependencies:
        values, edges = _combined_histogram(dependency, loader)
        counts = _rebin_to_contract(values, edges, dependency.sample_key)
        probabilities[dependency.sample_key] = _normalize(counts, dependency.sample_key)

    data_probability = probabilities[DATA_SAMPLE_KEY]
    family_maps = []
    for family, (sample_key, _histogram_key) in sorted(FAMILY_SPECS.items()):
        simulation_probability = probabilities[sample_key]
        weights = tuple(
            data_value / simulation_value
            for data_value, simulation_value in zip(data_probability, simulation_probability)
        )
        if len(weights) != len(BIN_EDGES_PERCENT) - 1 or any(
            not math.isfinite(value) or value <= 0.0 for value in weights
        ):
            raise ContractBuildError("invalid derived weight map for {}".format(family))
        family_maps.append(
            FamilyWeightMap(
                family=family,
                sample_key=sample_key,
                simulation_probability=simulation_probability,
                weights=weights,
            )
        )

    dependency_fingerprint = _sha256_payload(_dependency_payload(dependencies))
    contract_fingerprint = _contract_fingerprint(
        DIAGNOSTIC_CONTRACT_STATUS,
        dependency_fingerprint,
        data_probability,
        family_maps,
    )
    return CentralityWeightContract(
        schema=CONTRACT_SCHEMA,
        status=DIAGNOSTIC_CONTRACT_STATUS,
        dependency_fingerprint=dependency_fingerprint,
        contract_fingerprint=contract_fingerprint,
        dependencies=dependencies,
        data_probability=data_probability,
        family_maps=tuple(family_maps),
    )


def build_canonical_contract(
    *,
    authorities: Mapping[str, Mapping[str, Any]],
    data_counts: Optional[Sequence[float]] = None,
    data_probability: Optional[Sequence[float]] = None,
    family_counts: Optional[Mapping[str, Sequence[float]]] = None,
    family_probabilities: Optional[Mapping[str, Sequence[float]]] = None,
) -> CentralityWeightContract:
    """Build a READY contract from explicit hash-bound scientific authorities.

    Exactly one DATA representation and one family representation must be
    supplied.  Counts are normalized internally; probabilities must already be
    positive and normalized.  The exact supplied representation is included in
    the dependency fingerprint, so additional data or a calibration/selection/
    code authority change invalidates stale receipts even when a normalized
    shape happens to remain numerically identical.
    """
    authority_bindings = _authority_bindings(authorities)
    data_input, resolved_data_probability = _make_distribution_input(
        "data_target",
        counts=data_counts,
        probability=data_probability,
    )
    if (family_counts is None) == (family_probabilities is None):
        raise ContractBuildError(
            "supply exactly one of family_counts or family_probabilities"
        )
    supplied_families = family_counts if family_counts is not None else family_probabilities
    assert supplied_families is not None
    expected_families = set(FAMILY_SPECS)
    if set(supplied_families) != expected_families:
        raise ContractBuildError(
            "family distributions must contain exactly {}; found {}".format(
                sorted(expected_families),
                sorted(supplied_families),
            )
        )

    distribution_inputs = [data_input]
    family_maps = []
    for family, (sample_key, _histogram_key) in sorted(FAMILY_SPECS.items()):
        arguments: Dict[str, Optional[Sequence[float]]] = {
            "counts": None,
            "probability": None,
        }
        if family_counts is not None:
            arguments["counts"] = supplied_families[family]
        else:
            arguments["probability"] = supplied_families[family]
        distribution, simulation_probability = _make_distribution_input(
            family,
            counts=arguments["counts"],
            probability=arguments["probability"],
        )
        weights = tuple(
            data_value / simulation_value
            for data_value, simulation_value in zip(
                resolved_data_probability,
                simulation_probability,
            )
        )
        if any(not math.isfinite(value) or value <= 0.0 for value in weights):
            raise ContractBuildError("invalid canonical weight map for {}".format(family))
        distribution_inputs.append(distribution)
        family_maps.append(
            FamilyWeightMap(
                family=family,
                sample_key=sample_key,
                simulation_probability=simulation_probability,
                weights=weights,
            )
        )

    dependency_fingerprint = _sha256_payload(
        _canonical_dependency_payload(authority_bindings, distribution_inputs)
    )
    contract_fingerprint = _contract_fingerprint(
        CANONICAL_CONTRACT_STATUS,
        dependency_fingerprint,
        resolved_data_probability,
        family_maps,
    )
    return CentralityWeightContract(
        schema=CONTRACT_SCHEMA,
        status=CANONICAL_CONTRACT_STATUS,
        dependency_fingerprint=dependency_fingerprint,
        contract_fingerprint=contract_fingerprint,
        dependencies=(),
        data_probability=resolved_data_probability,
        family_maps=tuple(family_maps),
        authorities=authority_bindings,
        distribution_inputs=tuple(distribution_inputs),
    )


def ensure_canonical_contract(
    receipt_path: Path,
    *,
    authorities: Mapping[str, Mapping[str, Any]],
    data_counts: Optional[Sequence[float]] = None,
    data_probability: Optional[Sequence[float]] = None,
    family_counts: Optional[Mapping[str, Sequence[float]]] = None,
    family_probabilities: Optional[Mapping[str, Sequence[float]]] = None,
) -> CanonicalContractResolution:
    """Reuse the exact current receipt or rebuild it from changed authorities.

    The caller supplies the current Tree-to-hist sufficient statistics and the
    five exact authority bindings.  This function owns all freshness logic:
    it first constructs the contract that *should* exist, accepts an existing
    receipt only when its dependency and numerical fingerprints match exactly,
    and otherwise atomically replaces it.  Therefore new DATA/SIM counts or a
    calibration, selection, stitch, or derivation-code change cannot silently
    reuse a stale numerical map.

    This function does not reinterpret centrality.  If the selected Trees do
    not carry centrality under the declared calibration, the upstream producer
    must regenerate those centrality-bearing Trees before calling this entry
    point.
    """

    candidate = build_canonical_contract(
        authorities=authorities,
        data_counts=data_counts,
        data_probability=data_probability,
        family_counts=family_counts,
        family_probabilities=family_probabilities,
    )
    destination = Path(receipt_path)
    if destination.is_file():
        try:
            existing = load_contract_receipt(
                destination,
                expected_dependency_fingerprint=candidate.dependency_fingerprint,
                require_ready=True,
                verify_dependency_files=False,
            )
        except CentralityWeightContractError:
            action = "REBUILT_CHANGED_DEPENDENCY"
        else:
            if existing.contract_fingerprint != candidate.contract_fingerprint:
                # This should already be impossible when the full dependency
                # fingerprint agrees, but keep the numerical identity check as
                # an independent fail-closed invariant.
                action = "REBUILT_NUMERICAL_MISMATCH"
            else:
                return CanonicalContractResolution(existing, "REUSED_EXACT")
    else:
        action = "CREATED"

    candidate.write_receipt(destination)
    # Read back through the independent validator.  A successful atomic write
    # is not enough evidence if the serialized receipt cannot reproduce the
    # exact dependency and numerical fingerprints.
    resolved = load_contract_receipt(
        destination,
        expected_dependency_fingerprint=candidate.dependency_fingerprint,
        require_ready=True,
        verify_dependency_files=False,
    )
    return CanonicalContractResolution(resolved, action)


def _parse_artifact_dependencies(payload: Any) -> Tuple[ArtifactDependency, ...]:
    if not isinstance(payload, list):
        raise CentralityWeightContractError("receipt dependencies must be a list")
    dependencies = []
    for raw in payload:
        if not isinstance(raw, Mapping):
            raise CentralityWeightContractError("receipt dependency must be an object")
        roots_raw = raw.get("roots")
        if not isinstance(roots_raw, list) or not roots_raw:
            raise CentralityWeightContractError("receipt dependency has no ROOT bindings")
        roots = []
        for root_raw in roots_raw:
            if not isinstance(root_raw, Mapping):
                raise CentralityWeightContractError("receipt ROOT binding must be an object")
            path = root_raw.get("path")
            sha256 = root_raw.get("sha256")
            size_bytes = root_raw.get("size_bytes")
            if not isinstance(path, str) or not path or not _is_sha256(sha256):
                raise CentralityWeightContractError("receipt ROOT binding is malformed")
            if isinstance(size_bytes, bool) or not isinstance(size_bytes, int) or size_bytes < 0:
                raise CentralityWeightContractError("receipt ROOT size is invalid")
            roots.append(RootDependency(path, sha256, size_bytes))
        required_strings = (
            "sample_key",
            "histogram_key",
            "pointer_path",
            "current_entry_id",
            "canonical_status",
        )
        if any(not isinstance(raw.get(key), str) or not raw.get(key) for key in required_strings):
            raise CentralityWeightContractError("receipt artifact dependency is malformed")
        if not _is_sha256(raw.get("pointer_sha256")):
            raise CentralityWeightContractError("receipt pointer SHA-256 is malformed")
        dependencies.append(
            ArtifactDependency(
                sample_key=raw["sample_key"],
                histogram_key=raw["histogram_key"],
                pointer_path=raw["pointer_path"],
                pointer_sha256=raw["pointer_sha256"],
                current_entry_id=raw["current_entry_id"],
                campaign_tag=str(raw.get("campaign_tag") or ""),
                canonical_status=raw["canonical_status"],
                roots=tuple(roots),
            )
        )
    return tuple(dependencies)


def _parse_authorities(payload: Any) -> Tuple[AuthorityBinding, ...]:
    if not isinstance(payload, list):
        raise AuthorityValidationError("receipt authorities must be a list")
    role_mapping: Dict[str, Mapping[str, Any]] = {}
    for raw in payload:
        if not isinstance(raw, Mapping):
            raise AuthorityValidationError("receipt authority must be an object")
        row = dict(raw)
        role = row.pop("role", None)
        if not isinstance(role, str) or not role or role in role_mapping:
            raise AuthorityValidationError("receipt authority role is missing or duplicated")
        role_mapping[role] = row
    return _authority_bindings(role_mapping)


def _parse_distribution_inputs(payload: Any) -> Tuple[DistributionInput, ...]:
    if not isinstance(payload, list):
        raise ContractBuildError("receipt distribution_inputs must be a list")
    distributions = []
    seen = set()
    for raw in payload:
        if not isinstance(raw, Mapping):
            raise ContractBuildError("receipt distribution input must be an object")
        role = raw.get("role")
        kind = raw.get("kind")
        values = raw.get("values")
        if not isinstance(role, str) or not role or role in seen:
            raise ContractBuildError("receipt distribution role is missing or duplicated")
        if kind not in {"counts", "probability"} or not isinstance(values, list):
            raise ContractBuildError("receipt distribution input is malformed")
        if kind == "counts":
            distribution, _probability = _make_distribution_input(role, counts=values)
        else:
            distribution, _probability = _make_distribution_input(role, probability=values)
        distributions.append(distribution)
        seen.add(role)
    return tuple(distributions)


def _parse_family_maps(
    payload: Any,
    data_probability: Sequence[float],
) -> Tuple[FamilyWeightMap, ...]:
    if not isinstance(payload, Mapping) or set(payload) != set(FAMILY_SPECS):
        raise ContractBuildError("receipt family_maps do not match the canonical families")
    family_maps = []
    for family, (sample_key, _histogram_key) in sorted(FAMILY_SPECS.items()):
        raw = payload[family]
        if not isinstance(raw, Mapping):
            raise ContractBuildError("receipt family map must be an object")
        if raw.get("family") != family or raw.get("sample_key") != sample_key:
            raise ContractBuildError("receipt family map identity is inconsistent")
        simulation_probability = _validated_probability(
            raw.get("simulation_probability", ()),
            family,
        )
        weights = tuple(float(value) for value in raw.get("weights", ()))
        if len(weights) != len(BIN_EDGES_PERCENT) - 1:
            raise ContractBuildError("receipt family map has the wrong weight bin count")
        expected_weights = tuple(
            data_value / simulation_value
            for data_value, simulation_value in zip(data_probability, simulation_probability)
        )
        if any(
            not math.isfinite(value)
            or value <= 0.0
            or not math.isclose(value, expected, rel_tol=1.0e-14, abs_tol=0.0)
            for value, expected in zip(weights, expected_weights)
        ):
            raise ContractBuildError("receipt weights do not equal p_DATA/p_SIM for {}".format(family))
        family_maps.append(
            FamilyWeightMap(family, sample_key, simulation_probability, weights)
        )
    return tuple(family_maps)


def _verify_dependency_files(dependencies: Sequence[ArtifactDependency]) -> None:
    for dependency in dependencies:
        pointer = Path(dependency.pointer_path)
        if not pointer.is_file() or sha256_file(pointer) != dependency.pointer_sha256:
            raise StaleDependencyError(
                "current pointer changed since receipt creation: {}".format(pointer)
            )
        for root in dependency.roots:
            path = Path(root.path)
            if (
                not path.is_file()
                or path.stat().st_size != root.size_bytes
                or sha256_file(path) != root.sha256
            ):
                raise StaleDependencyError(
                    "ROOT dependency changed since receipt creation: {}".format(path)
                )


def validate_contract_receipt(
    payload: Mapping[str, Any],
    *,
    expected_dependency_fingerprint: Optional[str] = None,
    require_ready: bool = True,
    verify_dependency_files: bool = False,
) -> CentralityWeightContract:
    """Validate and reconstruct a serialized V2 contract receipt.

    ``expected_dependency_fingerprint`` is the consumer's calibration/input
    binding.  READY consumers should always provide it; a mismatch is stale and
    fails closed.  Diagnostic inspection may set ``require_ready=False``.
    """
    if not isinstance(payload, Mapping) or payload.get("schema") != CONTRACT_SCHEMA:
        raise CentralityWeightContractError("unsupported centrality contract receipt schema")
    status = payload.get("status")
    if status not in {DIAGNOSTIC_CONTRACT_STATUS, CANONICAL_CONTRACT_STATUS}:
        raise CentralityWeightContractError("unsupported centrality contract receipt status")
    dependency_fingerprint = payload.get("dependency_fingerprint")
    contract_fingerprint = payload.get("contract_fingerprint")
    if not _is_sha256(dependency_fingerprint) or not _is_sha256(contract_fingerprint):
        raise CentralityWeightContractError("receipt fingerprints are malformed")

    data_probability = _validated_probability(
        payload.get("data_probability", ()),
        "data_target",
    )
    family_maps = _parse_family_maps(payload.get("family_maps"), data_probability)
    dependencies = _parse_artifact_dependencies(payload.get("dependencies", []))

    if status == DIAGNOSTIC_CONTRACT_STATUS:
        if not dependencies:
            raise CurrentArtifactResolutionError("diagnostic receipt has no current dependencies")
        authorities = ()
        distribution_inputs = ()
        expected_internal_dependency = _sha256_payload(_dependency_payload(dependencies))
    else:
        if dependencies:
            raise AuthorityValidationError("canonical receipt cannot substitute current pointers for authorities")
        authorities = _parse_authorities(payload.get("authorities", []))
        distribution_inputs = _parse_distribution_inputs(
            payload.get("distribution_inputs", [])
        )
        distribution_by_role = {item.role: item for item in distribution_inputs}
        expected_roles = {"data_target"} | set(FAMILY_SPECS)
        if set(distribution_by_role) != expected_roles:
            raise ContractBuildError(
                "canonical receipt distribution roles must be exactly {}".format(
                    sorted(expected_roles)
                )
            )
        resolved_probabilities = {}
        for role, distribution in distribution_by_role.items():
            if distribution.kind == "counts":
                resolved_probabilities[role] = _normalize(distribution.values, role)
            else:
                resolved_probabilities[role] = _validated_probability(
                    distribution.values,
                    role,
                )
        if any(
            not math.isclose(left, right, rel_tol=0.0, abs_tol=1.0e-15)
            for left, right in zip(resolved_probabilities["data_target"], data_probability)
        ):
            raise ContractBuildError("canonical DATA distribution does not reproduce receipt probability")
        for family_map in family_maps:
            if any(
                not math.isclose(left, right, rel_tol=0.0, abs_tol=1.0e-15)
                for left, right in zip(
                    resolved_probabilities[family_map.family],
                    family_map.simulation_probability,
                )
            ):
                raise ContractBuildError(
                    "canonical {} distribution does not reproduce receipt probability".format(
                        family_map.family
                    )
                )
        expected_internal_dependency = _sha256_payload(
            _canonical_dependency_payload(authorities, distribution_inputs)
        )

    if expected_internal_dependency != dependency_fingerprint:
        raise StaleDependencyError("receipt dependency fingerprint does not match its bindings")
    expected_contract = _contract_fingerprint(
        status,
        dependency_fingerprint,
        data_probability,
        family_maps,
    )
    if expected_contract != contract_fingerprint:
        raise CentralityWeightContractError(
            "receipt contract fingerprint does not match its numerical maps"
        )
    if expected_dependency_fingerprint is not None:
        if not _is_sha256(expected_dependency_fingerprint):
            raise StaleDependencyError("expected dependency fingerprint is not a SHA-256 digest")
        if expected_dependency_fingerprint != dependency_fingerprint:
            raise StaleDependencyError(
                "receipt is stale for the expected dependency fingerprint"
            )
    if verify_dependency_files:
        _verify_dependency_files(dependencies)

    contract = CentralityWeightContract(
        schema=CONTRACT_SCHEMA,
        status=status,
        dependency_fingerprint=dependency_fingerprint,
        contract_fingerprint=contract_fingerprint,
        dependencies=dependencies,
        data_probability=data_probability,
        family_maps=family_maps,
        authorities=authorities,
        distribution_inputs=distribution_inputs,
    )
    if require_ready:
        if expected_dependency_fingerprint is None:
            raise StaleDependencyError(
                "READY receipt validation requires an expected dependency fingerprint"
            )
        contract.assert_ready_for_application(expected_dependency_fingerprint)
    return contract


def load_contract_receipt(
    path: Path,
    *,
    expected_dependency_fingerprint: Optional[str] = None,
    require_ready: bool = True,
    verify_dependency_files: bool = False,
) -> CentralityWeightContract:
    """Load a JSON receipt and apply :func:`validate_contract_receipt`."""
    receipt_path = Path(path)
    if not receipt_path.is_file():
        raise CentralityWeightContractError(
            "centrality contract receipt does not exist: {}".format(receipt_path)
        )
    try:
        payload = json.loads(receipt_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise CentralityWeightContractError(
            "centrality contract receipt is not valid JSON: {}".format(receipt_path)
        ) from error
    return validate_contract_receipt(
        payload,
        expected_dependency_fingerprint=expected_dependency_fingerprint,
        require_ready=require_ready,
        verify_dependency_files=verify_dependency_files,
    )


__all__ = [
    "APPLICATION_ORDER",
    "BIN_EDGES_PERCENT",
    "BIN_WIDTH_PERCENT",
    "CANONICAL_CONTRACT_STATUS",
    "CANONICAL_STITCH_COMPONENT",
    "CENTRALITY_COMPONENT_PREFIX",
    "CONTRACT_SCHEMA",
    "CONTRACT_STATUS",
    "DIAGNOSTIC_CONTRACT_STATUS",
    "AuthorityBinding",
    "AuthorityValidationError",
    "CanonicalContractResolution",
    "CentralityWeightContract",
    "CentralityWeightContractError",
    "ContractBuildError",
    "ContractNotReadyError",
    "CurrentArtifactResolutionError",
    "DoubleApplicationError",
    "EventWeightState",
    "FAMILY_SPECS",
    "REQUIRED_AUTHORITY_ROLES",
    "StaleDependencyError",
    "UnsupportedCentralityError",
    "UnsupportedSampleFamilyError",
    "WeightApplicationOrderError",
    "build_canonical_contract",
    "build_current_contract",
    "ensure_canonical_contract",
    "load_contract_receipt",
    "resolve_current_inputs",
    "sha256_file",
    "validate_contract_receipt",
]
