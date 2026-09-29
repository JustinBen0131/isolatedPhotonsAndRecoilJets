"""Pure AuAu range, centrality-provenance, and finite short-job planning contracts.

No filesystem, scheduler, subprocess, network, clock, or execution backend is
used here. Callers supply verified manifests, hashes, terminal receipts, and a
fresh owner-wide queue snapshot. These checks do not verify remote bytes, grant
submission authority, replace WalltimeShardingAdmissionV1, or assert scientific
acceptance. Coverage applies only to the explicitly supplied range universe.

Historical products without reconstructible per-run calibration hashes remain
UNVERIFIED. A directory name or an old queue total never proves validity.
Superseded products must be excluded explicitly by the caller: overlapping or
duplicate products are rejected, never silently selected or counted twice.
Centrality-contaminated products need exact-range replacement because centrality
affects selections; this module offers no scalar centrality patch.

Augmentation integration (existing qa_augmentation.py / THE-357 preparation):
* Reuse positive_charge/make_columns/build/verify for calibrated finite-positive
  PMT totals, keyed by source occurrence and event identity (plus source index
  in v2). Preserve base hashes, finite values, and event/candidate populations.
* Native MbdOut get_t0/get_time(arm), cluster mean time, and missing calibrated
  calorimeter sums need original paired-DST donors with the pinned runtime.
  Prove physical-event and rebuilt-candidate identity; preserve raw times and
  establish sample-to-ns/run-offset conventions before subtracting times.
* NPB needs the applicable named model/hash, all 25 ordered inputs, score,
  validity, pass state, and distinct reference/active preselection state. The
  retained 14-feature tight-ID vector is insufficient; no zero imputation.
* Inspect existing gl1_trigger_interface, schema14_payload, exporter and
  readback coverage for GL1/truth vertices; absence must not be presumed.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import json
import math
import re
from types import MappingProxyType
from typing import Any, Mapping, Sequence


NOMINAL_SET = "sam_fit_20260913_v1"
AB_ONLY_SET = "sam_datapct_20260918_v1"
# FrozenCalibrationSetsV1, 2026-09-18: the nominal original fit set.
NOMINAL_MANIFEST_SHA256 = "4a11e0e10ae9b82062325c3275679935da2259f4c91f796950a6b5ddd2e1501c"
MAX_EVENTS_PER_JOB = 35_000
RATE_FLOOR = 0.51
MAX_OCCUPANCY = 30_000
LOW_WATER = 20_000
EMERGENCY_WATER = 5_000


class ContractError(ValueError):
    """Input cannot support an unambiguous, bounded proposal."""


def _integer(value: Any, name: str, minimum: int = 0) -> None:
    if type(value) is not int or value < minimum:
        raise ContractError(f"{name} must be an integer >= {minimum}")


def _number(value: Any, name: str, minimum: float = 0) -> None:
    if (type(value) not in (int, float) or not math.isfinite(value)
            or value < minimum):
        raise ContractError(f"{name} must be finite and >= {minimum}")


def _text(value: Any, name: str) -> None:
    if not isinstance(value, str) or not value.strip():
        raise ContractError(f"{name} must be nonempty text")


def _sha(value: Any, name: str) -> None:
    if not isinstance(value, str) or not re.fullmatch(r"[0-9a-f]{64}", value):
        raise ContractError(f"{name} must be a lowercase SHA256")


def _path(value: Any, name: str) -> None:
    _text(value, name)
    if not value.startswith("/") or any(x in (".", "..", "") for x in value[1:].split("/")):
        raise ContractError(f"{name} must be a literal absolute file path")


@dataclass(frozen=True)
class SourceRange:
    row_id: str
    run: int
    segment: int
    jets: str
    jetcalo: str
    source_total_events: int
    event_offset: int
    event_count: int
    sample_id: str = "auau"

    def __post_init__(self) -> None:
        _text(self.row_id, "row_id")
        _integer(self.run, "run", 1)
        _integer(self.segment, "segment")
        _integer(self.source_total_events, "source_total_events", 1)
        _integer(self.event_offset, "event_offset")
        _integer(self.event_count, "event_count", 1)
        _path(self.jets, "jets")
        _path(self.jetcalo, "jetcalo")
        if self.jets == self.jetcalo:
            raise ContractError("paired sources must be distinct")
        if self.sample_id != "auau":
            raise ContractError("this contract is for AuAu DATA only")
        if self.end > self.source_total_events:
            raise ContractError("range exceeds exact source_total_events")

    @property
    def end(self) -> int:
        return self.event_offset + self.event_count

    @property
    def pair(self) -> tuple[str, str]:
        return self.jets, self.jetcalo

    @property
    def identity(self) -> tuple:
        # row_id is an alias; changing it cannot disguise duplicate exposure.
        return (self.run, self.segment, *self.pair, self.source_total_events,
                self.event_offset, self.event_count)


def range_from_manifest(row: Mapping[str, Any]) -> SourceRange:
    """Adapt one existing ranges.jsonl atom without deriving paired filenames."""
    try:
        return SourceRange(row_id=row["row_id"], run=row["run"], segment=row["segment"],
                           jets=row["source_paths"]["jets"],
                           jetcalo=row["source_paths"]["jetcalo"],
                           source_total_events=row["source_total_events"],
                           event_offset=row["event_offset"], event_count=row["event_count"],
                           sample_id=row["sample_id"])
    except (KeyError, TypeError) as exc:
        raise ContractError("incomplete paired-source atom") from exc


def validate_ranges(rows: Sequence[SourceRange]) -> tuple[SourceRange, ...]:
    """Reject aliases, duplicate rows, conflicting pair metadata, and overlap.

    Gaps are permitted: a remainder universe need not contain an entire source.
    Returned rows have deterministic order; no source range is changed.
    """
    by_id: set[str] = set()
    source_metadata: dict[str, tuple] = {}
    pairs: dict[tuple, list[SourceRange]] = {}
    for row in rows:
        if not isinstance(row, SourceRange):
            raise ContractError("ranges must contain SourceRange objects")
        if row.row_id in by_id:
            raise ContractError(f"duplicate row_id: {row.row_id}")
        by_id.add(row.row_id)
        metadata = (row.run, row.segment, row.source_total_events, *row.pair)
        for path in row.pair:
            if path in source_metadata and source_metadata[path] != metadata:
                raise ContractError("source reused with conflicting pair/run/segment/total")
            source_metadata[path] = metadata
        pairs.setdefault(row.pair, []).append(row)
    for group in pairs.values():
        group.sort(key=lambda r: (r.event_offset, r.end))
        for left, right in zip(group, group[1:]):
            if right.event_offset < left.end:
                raise ContractError(f"duplicate/overlapping ranges: {left.row_id}, {right.row_id}")
    return tuple(sorted(rows, key=lambda r: (r.run, r.segment, r.pair, r.event_offset)))


@dataclass(frozen=True)
class Payload:
    path: str
    sha256: str

    def __post_init__(self) -> None:
        _path(self.path, "payload path")
        _sha(self.sha256, "payload sha256")


@dataclass(frozen=True)
class CentralityTriplet:
    run: int
    divs: Payload
    scales: Payload
    vertexscales: Payload

    def __post_init__(self) -> None:
        _integer(self.run, "calibration run", 1)
        for payload in (self.divs, self.scales, self.vertexscales):
            if not isinstance(payload, Payload):
                raise ContractError("centrality requires all three payloads and hashes")
        if len({p.path for p in (self.divs, self.scales, self.vertexscales)}) != 3:
            raise ContractError("centrality triplet payload paths must be distinct")

    @property
    def hashes(self) -> tuple[str, str, str]:
        return self.divs.sha256, self.scales.sha256, self.vertexscales.sha256


@dataclass(frozen=True)
class FallbackApproval:
    run: int
    evidence_sha256: str
    own_triplet_absent: bool
    target_run: int = 68144

    def __post_init__(self) -> None:
        _integer(self.run, "fallback source run", 1)
        _sha(self.evidence_sha256, "fallback legitimacy evidence")
        if self.own_triplet_absent is not True or self.target_run != 68144 or self.run == 68144:
            raise ContractError("fallback requires proven absent own triplet and target 68144")


@dataclass(frozen=True)
class CalibrationManifest:
    set_id: str
    manifest_sha256: str
    triplets: tuple[CentralityTriplet, ...]
    fallbacks: tuple[FallbackApproval, ...] = ()
    _by_run: Mapping[int, CentralityTriplet] = field(init=False, repr=False, compare=False)
    _fallback_by_run: Mapping[int, FallbackApproval] = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        if self.set_id != NOMINAL_SET or self.manifest_sha256 != NOMINAL_MANIFEST_SHA256:
            raise ContractError("production requires the pinned nominal fit manifest; data-percentile is A/B only")
        if type(self.triplets) is not tuple or not self.triplets:
            raise ContractError("provide immutable exact-run triplets from the verified manifest")
        if type(self.fallbacks) is not tuple:
            raise ContractError("fallback approvals must be immutable")
        if any(not isinstance(t, CentralityTriplet) for t in self.triplets):
            raise ContractError("invalid calibration triplet")
        runs = [t.run for t in self.triplets]
        if len(runs) != len(set(runs)):
            raise ContractError("duplicate calibration run")
        if any(not isinstance(a, FallbackApproval) for a in self.fallbacks):
            raise ContractError("invalid fallback approval")
        allowed = [a.run for a in self.fallbacks]
        if len(allowed) != len(set(allowed)) or set(allowed) & set(runs):
            raise ContractError("duplicate fallback or fallback despite available own triplet")
        if allowed and 68144 not in runs:
            raise ContractError("fallback target triplet is missing")
        object.__setattr__(self, "_by_run", MappingProxyType({t.run: t for t in self.triplets}))
        object.__setattr__(self, "_fallback_by_run", MappingProxyType({a.run: a for a in self.fallbacks}))

    def expected(self, run: int) -> CentralityTriplet:
        own = self._by_run.get(run)
        if own is not None:
            return own
        if run in self._fallback_by_run:
            return self._by_run[68144]
        raise ContractError(f"no complete own centrality triplet or approved fallback for run {run}")

    def fallback_evidence(self, run: int) -> str | None:
        approval = self._fallback_by_run.get(run)
        return approval.evidence_sha256 if approval is not None else None


def manifest_from_frozen_bytes(raw_manifest: bytes, required_runs: Sequence[int], *,
                               fallbacks: tuple[FallbackApproval, ...] = ()) -> CalibrationManifest:
    """Bind selected triplets to the actual pinned manifest bytes, not labels.

    This verifies manifest membership/hashes, not present remote payload bytes
    or scientific fallback approval. Callers still verify both before admission.
    Partial own-run triplets fail; components are never mixed with fallback.
    """
    if not isinstance(raw_manifest, bytes) or not 0 < len(raw_manifest) <= 2*1024**2:
        raise ContractError("bounded frozen manifest bytes required")
    digest = hashlib.sha256(raw_manifest).hexdigest()
    if digest != NOMINAL_MANIFEST_SHA256:
        raise ContractError("frozen manifest bytes differ from pinned nominal hash")

    def unique_object(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ContractError("duplicate frozen manifest key")
            result[key] = value
        return result

    data = json.loads(raw_manifest, object_pairs_hook=unique_object)
    if (not isinstance(data, dict) or data.get("schema") != "FrozenCentralityCalibrationV1"
            or type(data.get("copy_mismatches")) is not int or data["copy_mismatches"] != 0
            or not isinstance(data.get("files"), dict)
            or type(data.get("file_count")) is not int or data["file_count"] != len(data["files"])):
        raise ContractError("invalid frozen calibration manifest structure")
    base = data.get("frozen_base")
    _path(base, "frozen calibration base")
    files = data["files"]
    for value in files.values():
        _sha(value, "frozen manifest payload hash")
    runs = tuple(required_runs)
    if not runs or len(runs) != len(set(runs)):
        raise ContractError("nonempty unique required runs needed")
    for run in runs:
        _integer(run, "required run", 1)
    if type(fallbacks) is not tuple or any(not isinstance(x,FallbackApproval) for x in fallbacks):
        raise ContractError("explicit immutable fallback approvals required")
    if any(a.run not in runs for a in fallbacks):
        raise ContractError("fallback outside requested run scope")
    selected = {}
    for run in set(runs) | ({68144} if fallbacks else set()):
        names = (f"divs/cdb_centrality_{run}.root", f"scales/cdb_centrality_scale_{run}.root",
                 f"vertexscales/cdb_centrality_vertex_scale_{run}.root")
        present = [name in files for name in names]
        if any(present) and not all(present):
            raise ContractError(f"partial own centrality triplet for run {run}")
        if all(present):
            selected[run] = CentralityTriplet(run, *(Payload(base+"/"+name,files[name]) for name in names))
    manifest = CalibrationManifest(NOMINAL_SET,digest,tuple(selected[r] for r in sorted(selected)),fallbacks)
    for run in runs:
        manifest.expected(run)  # Missing own triplet never silently becomes a fallback.
    return manifest


@dataclass(frozen=True)
class CentralityBinding:
    requested_run: int
    set_id: str
    manifest_sha256: str
    triplet: CentralityTriplet
    fallback: bool
    evidence_sha256: str

    def __post_init__(self) -> None:
        _integer(self.requested_run, "requested_run", 1)
        _text(self.set_id, "binding set_id")
        _sha(self.manifest_sha256, "binding manifest_sha256")
        _sha(self.evidence_sha256, "binding evidence_sha256")
        if not isinstance(self.triplet, CentralityTriplet) or type(self.fallback) is not bool:
            raise ContractError("binding must retain exact triplet and boolean fallback")


def binding_matches(binding: CentralityBinding, run: int, manifest: CalibrationManifest) -> bool:
    """Compare byte hashes, not directory names; old paths may hold identical bytes."""
    expected = manifest.expected(run)
    return (binding.requested_run == run and binding.set_id == manifest.set_id
            and binding.manifest_sha256 == manifest.manifest_sha256
            and binding.triplet.run == expected.run and binding.triplet.hashes == expected.hashes
            and binding.fallback == (expected.run != run))


@dataclass(frozen=True)
class ProductEvidence:
    """Caller-verified terminal evidence, not a counts-only acceptance shortcut.

    For upstream rejections, source_accounting_sha256 must identify the audit
    binding these exact product bytes and source range to a disjoint, complete
    retained/rejected identity partition (the native finalizer checks it).
    This pure planner cannot authenticate an arbitrary caller-supplied digest.
    """
    source_range: SourceRange
    terminal_success: bool
    processed_events: int
    replay_events: int
    product_sha256: str
    receipt_sha256: str
    binding: CentralityBinding | None
    upstream_rejected_events: int = 0
    source_accounting_sha256: str | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.source_range, SourceRange) or type(self.terminal_success) is not bool:
            raise ContractError("product requires exact range and boolean terminal success")
        _integer(self.processed_events, "processed_events")
        _integer(self.replay_events, "replay_events")
        _integer(self.upstream_rejected_events, "upstream_rejected_events")
        if self.source_accounting_sha256 is not None:
            _sha(self.source_accounting_sha256, "source accounting evidence")
        if self.upstream_rejected_events and self.source_accounting_sha256 is None:
            raise ContractError("upstream rejections require exact source-accounting evidence")
        _sha(self.product_sha256, "product_sha256")
        _sha(self.receipt_sha256, "receipt_sha256")
        if self.binding is not None and not isinstance(self.binding, CentralityBinding):
            raise ContractError("invalid product binding")


@dataclass(frozen=True)
class RangeDecision:
    source_range: SourceRange
    state: str
    reason: str


@dataclass(frozen=True)
class CoverageReport:
    decisions: tuple[RangeDecision, ...]
    coverage_scope: str = "DECLARED_RANGE_UNIVERSE_ONLY"
    scientific_acceptance: bool = False

    @property
    def complete(self) -> bool:
        return bool(self.decisions) and all(d.state == "KEEP_VALID" for d in self.decisions)


def classify_coverage(universe: Sequence[SourceRange], products: Sequence[ProductEvidence],
                      manifest: CalibrationManifest, *,
                      allow_unbound_runs: bool = False) -> CoverageReport:
    if type(allow_unbound_runs) is not bool:
        raise ContractError("allow_unbound_runs must be boolean")
    rows = validate_ranges(universe)
    validate_ranges(tuple(p.source_range for p in products))
    expected = {r.identity: r for r in rows}
    by_range = {}
    for product in products:
        if product.source_range.identity not in expected:
            raise ContractError("product range is outside or not exactly equal to a declared atom")
        by_range[product.source_range.identity] = product
    decisions = []
    for row in rows:
        try:
            manifest.expected(row.run)
        except ContractError:
            if not allow_unbound_runs:
                raise
            # Diagnostic partition only: these rows never enter a job plan.
            decisions.append(RangeDecision(row, "BLOCKED_CALIBRATION",
                "no complete own triplet or explicit approved fallback"))
            continue
        product = by_range.get(row.identity)
        if product is None:
            state, reason = "NEEDS_PRODUCTION", "no terminal product evidence supplied"
        elif (not product.terminal_success or product.processed_events != row.event_count
              or product.replay_events + product.upstream_rejected_events != row.event_count):
            state, reason = "UNVERIFIED", "terminal success and exact processed/retained/rejected counts required"
        elif product.binding is None:
            state, reason = "UNVERIFIED", "per-run triplet binding evidence missing"
        elif not binding_matches(product.binding, row.run, manifest):
            state, reason = "REPLACE_CENTRALITY", "off-nominal binding; replace this exact selected range"
        else:
            state, reason = "KEEP_VALID", "exact range, terminal counts, and nominal triplet hashes match"
        decisions.append(RangeDecision(row, state, reason))
    return CoverageReport(tuple(decisions))


@dataclass(frozen=True)
class PlannedJob:
    source_range: SourceRange
    calibration: CentralityTriplet
    range_proof_sha256: str
    rate_evidence_sha256: str
    skip_rate_evidence_sha256: str
    process_seconds_per_event: float
    skip_seconds_per_event: float
    safety_factor: float
    projected_wall_seconds: float
    target_wall_seconds: float
    fallback_approval_sha256: str | None = None
    calibration_set_id: str = NOMINAL_SET
    calibration_manifest_sha256: str = NOMINAL_MANIFEST_SHA256
    hard_planned_ceiling_seconds: float = 12 * 3600
    execution_authorized: bool = False

    def __post_init__(self) -> None:
        if not isinstance(self.source_range, SourceRange) or not isinstance(self.calibration, CentralityTriplet):
            raise ContractError("planned job requires one source range and one calibration triplet")
        for value in (self.range_proof_sha256, self.rate_evidence_sha256, self.skip_rate_evidence_sha256):
            _sha(value, "planned evidence hash")
        if self.calibration_set_id != NOMINAL_SET or self.calibration_manifest_sha256 != NOMINAL_MANIFEST_SHA256:
            raise ContractError("planned calibration must remain pinned to nominal")
        if self.calibration.run == self.source_range.run:
            if self.fallback_approval_sha256 is not None:
                raise ContractError("own-run calibration cannot carry fallback approval")
        else:
            _sha(self.fallback_approval_sha256, "planned fallback approval")
            if self.calibration.run != 68144:
                raise ContractError("only approved 68144 fallback is supported")
        _number(self.process_seconds_per_event, "planned process rate", RATE_FLOOR)
        _number(self.skip_seconds_per_event, "planned skip rate")
        _number(self.safety_factor, "planned safety_factor", 1.5)
        _number(self.target_wall_seconds, "planned target", 1)
        _number(self.projected_wall_seconds, "planned wall time", 0)
        expected = (self.source_range.event_count * self.process_seconds_per_event
                    + self.source_range.event_offset * self.skip_seconds_per_event) * self.safety_factor
        if (self.source_range.event_count > MAX_EVENTS_PER_JOB or self.target_wall_seconds > 8 * 3600
                or self.projected_wall_seconds != expected or expected > self.target_wall_seconds
                or self.hard_planned_ceiling_seconds != 12 * 3600
                or (self.source_range.event_offset and self.skip_seconds_per_event <= 0)
                or self.execution_authorized is not False):
            raise ContractError("planned job violates the finite short-job envelope")

    @property
    def deterministic_skip_events(self) -> int:
        return self.source_range.event_offset


def centrality_environment(job: PlannedJob, manifest: CalibrationManifest, *,
                           inherited: Mapping[str, str] | None = None) -> dict[str, str]:
    """Bind one job to the native local-file API, never Sam's mutable directory.

    The caller must verify the frozen manifest and payload bytes. Native
    resolve() rechecks all three payload hashes before registering CentralityReco.
    This adapter is not submission or scientific acceptance. Every old centrality
    variable is discarded so an inherited CDB/base/fallback override cannot win.
    """
    if not isinstance(job, PlannedJob) or not isinstance(manifest, CalibrationManifest):
        raise ContractError("typed planned job and pinned calibration manifest required")
    run = job.source_range.run
    expected = manifest.expected(run)
    if (job.calibration != expected or job.calibration_set_id != manifest.set_id
            or job.calibration_manifest_sha256 != manifest.manifest_sha256
            or job.fallback_approval_sha256 != manifest.fallback_evidence(run)):
        raise ContractError("planned centrality differs from frozen per-run binding")
    return centrality_environment_for_run(run,manifest,inherited=inherited)


def centrality_environment_for_run(run: int, manifest: CalibrationManifest, *,
                                  inherited: Mapping[str,str] | None = None) -> dict[str,str]:
    """Shared frozen per-run binding for planners and the hash-sealed batch worker."""
    _integer(run,"requested run",1)
    if not isinstance(manifest,CalibrationManifest):
        raise ContractError("verified frozen manifest required")
    expected=manifest.expected(run)
    specs = (("DIVS", expected.divs, f"cdb_centrality_{expected.run}.root"),
             ("SCALE", expected.scales, f"cdb_centrality_scale_{expected.run}.root"),
             ("VERTEX_SCALE", expected.vertexscales,
              f"cdb_centrality_vertex_scale_{expected.run}.root"))
    for _, payload, basename in specs:
        if payload.path.rsplit("/", 1)[-1] != basename:
            raise ContractError("payload basename does not identify the native calibration run")
    environment = {}
    for key, value in (inherited or {}).items():
        if not isinstance(key, str) or not isinstance(value, str):
            raise ContractError("inherited environment must contain literal string keys/values")
        if not key.startswith("RJ_AUAU_CENTRALITY_"):
            environment[key] = value
    prefix = "RJ_AUAU_CENTRALITY_"
    environment.update({prefix + "SOURCE": "local", prefix + "SET_ID": manifest.set_id,
        prefix + "MANIFEST_SHA256": manifest.manifest_sha256,
        prefix + "PAYLOAD_RUN": str(expected.run)})
    for key, payload, _ in specs:
        environment[prefix + key] = payload.path
        environment[prefix + key + "_SHA256"] = payload.sha256
    if manifest.fallback_evidence(run) is not None:
        environment[prefix + "FALLBACK_EVIDENCE_SHA256"] = manifest.fallback_evidence(run)
    return environment


def plan_short_jobs(rows: Sequence[SourceRange], manifest: CalibrationManifest, *,
                    range_proof_sha256: str, process_seconds_per_event: float,
                    rate_evidence_sha256: str, skip_seconds_per_event: float,
                    skip_rate_evidence_sha256: str, safety_factor: float = 1.5,
                    target_wall_seconds: float = 8 * 3600) -> tuple[PlannedJob, ...]:
    """One existing atomic paired range per job; never group, split, or add kills.

    Offsets model independent deterministic skipping from source entry zero.
    The caller must supply measured skip cost and its evidence, even though the
    existing admission helper currently projects only processed-event cost.
    An over-target atom is rejected for a separately proven partition/access
    strategy; reducing its event count here could not remove its skip overhead.
    """
    ordered = validate_ranges(rows)
    for value, name in ((range_proof_sha256, "range proof"),
                        (rate_evidence_sha256, "process rate evidence"),
                        (skip_rate_evidence_sha256, "skip rate evidence")):
        _sha(value, name)
    _number(process_seconds_per_event, "process_seconds_per_event")
    _number(skip_seconds_per_event, "skip_seconds_per_event")
    _number(safety_factor, "safety_factor", 1.5)
    _number(target_wall_seconds, "target_wall_seconds", 1)
    if target_wall_seconds > 8 * 3600:
        raise ContractError("short-job target must not exceed eight hours")
    if any(r.event_offset for r in ordered) and skip_seconds_per_event <= 0:
        raise ContractError("nonzero offsets require a measured positive deterministic skip cost")
    rate = max(RATE_FLOOR, process_seconds_per_event)
    jobs = []
    for row in ordered:
        if row.event_count > MAX_EVENTS_PER_JOB:
            raise ContractError("atomic range exceeds the 35000-event AuAu cap")
        projected = (row.event_count * rate + row.event_offset * skip_seconds_per_event) * safety_factor
        if not math.isfinite(projected) or projected > target_wall_seconds or projected > 12 * 3600:
            raise ContractError(f"{row.row_id} including deterministic skip exceeds short-job target")
        jobs.append(PlannedJob(row, manifest.expected(row.run), range_proof_sha256,
                               rate_evidence_sha256, skip_rate_evidence_sha256,
                               rate, skip_seconds_per_event, safety_factor,
                               projected, target_wall_seconds, manifest.fallback_evidence(row.run)))
    return tuple(jobs)


@dataclass(frozen=True)
class ContinuationPlan:
    coverage: CoverageReport
    jobs: tuple[PlannedJob, ...]
    execution_authorized: bool = False


def plan_continuation(universe: Sequence[SourceRange], products: Sequence[ProductEvidence],
                      manifest: CalibrationManifest, **planning_parameters: Any) -> ContinuationPlan:
    """Connect exact coverage to finite work without rerunning accepted atoms.

    Unknown counts/provenance are UNVERIFIED, not permission to rerun. Missing
    calibration is isolated as BLOCKED_CALIBRATION, never an implicit fallback.
    KEEP_VALID concerns the source/centrality gate only: its products may still
    need the separately planned augmentation and canonical conversion.
    Replacement namespace/authority, native qualification, and admission remain
    required before an executor can consume the proposed jobs.
    """
    coverage = classify_coverage(universe, products, manifest, allow_unbound_runs=True)
    needed = tuple(d.source_range for d in coverage.decisions
                   if d.state in ("NEEDS_PRODUCTION", "REPLACE_CENTRALITY"))
    jobs = plan_short_jobs(needed, manifest, **planning_parameters)
    return ContinuationPlan(coverage, jobs)


@dataclass(frozen=True)
class QueueSnapshot:
    """Owner-wide queue membership, including jobs not currently runnable.

    All counters are required, even when zero. ``other_queued`` includes
    removing/completed jobs still present in the scheduler and any additional
    reported states. These disjoint counters must account for every observed
    job across the owner's submit hosts; held/suspended jobs are not free slots.
    """
    observed_at: float
    evidence_sha256: str
    idle: int
    running: int
    transferring: int
    held: int
    suspended: int
    other_queued: int
    inflight_ranges: tuple[SourceRange, ...]
    owner_wide_complete: bool

    def __post_init__(self) -> None:
        _number(self.observed_at, "observed_at")
        _sha(self.evidence_sha256, "queue evidence")
        for name in ("idle", "running", "transferring", "held", "suspended", "other_queued"):
            _integer(getattr(self, name), name)
        if self.owner_wide_complete is not True:
            raise ContractError("all owner jobs in every queue state must be counted")
        if type(self.inflight_ranges) is not tuple:
            raise ContractError("inflight_ranges must be immutable")
        validate_ranges(self.inflight_ranges)
        if len(self.inflight_ranges) > self.occupancy:
            raise ContractError("inflight ranges exceed total queue occupancy")

    @property
    def occupancy(self) -> int:
        return (self.idle + self.running + self.transferring + self.held
                + self.suspended + self.other_queued)


@dataclass(frozen=True)
class ReplenishmentProposal:
    jobs: tuple[PlannedJob, ...]
    state: str
    occupancy_before: int
    projected_occupancy: int
    snapshot_sha256: str
    execution_authorized: bool = False
    reserved_other_lane_slots: int = 0


def propose_replenishment(jobs: Sequence[PlannedJob], snapshot: QueueSnapshot, *,
                           previously_submitted: Sequence[SourceRange], now: float,
                           max_new_jobs: int = MAX_OCCUPANCY,
                           max_snapshot_age_seconds: float = 300,
                           reserved_other_lane_slots: int = 0) -> ReplenishmentProposal:
    """Return finite rows only; never re-propose any previously submitted range.

    previously_submitted is the complete exact submission ledger for this plan,
    including drained/held/failed rows, not merely the currently visible queue.
    Replacement/retry authority and namespaces are a separate caller obligation.
    Snapshot completeness, freshness and ledger evidence must be established by
    the caller; this function performs no scheduler observation or submission.
    Other-lane reservations are additional not-yet-queued slots, not a second
    count of pp/SIM jobs already included in owner-wide occupancy. They constrain
    admission only, not Condor priority or guaranteed execution throughput.
    """
    _number(now, "now")
    _number(max_snapshot_age_seconds, "max_snapshot_age_seconds", 1)
    _integer(max_new_jobs, "max_new_jobs", 1)
    _integer(reserved_other_lane_slots, "reserved_other_lane_slots")
    if reserved_other_lane_slots > MAX_OCCUPANCY:
        raise ContractError("other-lane reservation exceeds owner-wide cap")
    if max_new_jobs > MAX_OCCUPANCY or max_snapshot_age_seconds > 300:
        raise ContractError("proposal count/freshness cannot widen the bounded envelope")
    if not 0 <= now - snapshot.observed_at <= max_snapshot_age_seconds:
        raise ContractError("queue snapshot is stale or from the future")
    validate_ranges(tuple(j.source_range for j in jobs))
    submitted = validate_ranges(previously_submitted)
    submitted_ids = {r.identity for r in submitted}
    if any(r.identity not in submitted_ids for r in snapshot.inflight_ranges):
        raise ContractError("inflight range absent from submission ledger")
    # Deduplicate exact identities for overlap validation across plan/history;
    # aliases do not grant a second submission, and partial overlaps fail closed.
    combined = {r.identity: r for r in submitted}
    for job in jobs:
        combined.setdefault(job.source_range.identity, job.source_range)
    validate_ranges(tuple(combined.values()))
    pending = tuple(j for j in jobs if j.source_range.identity not in submitted_ids)
    count = snapshot.occupancy
    if not pending:
        state, selected = "NO_UNSUBMITTED_ROWS", ()
    elif count >= MAX_OCCUPANCY:
        state, selected = "AT_OR_ABOVE_CAP", ()
    elif count >= LOW_WATER:
        state, selected = "ABOVE_LOW_WATER", ()
    elif count + reserved_other_lane_slots >= MAX_OCCUPANCY:
        state, selected = "CAPACITY_RESERVED_FOR_OTHER_LANES", ()
    else:
        state = "EMERGENCY_REFILL_PROPOSAL" if count <= EMERGENCY_WATER else "REFILL_PROPOSAL"
        selected = pending[:min(MAX_OCCUPANCY - count - reserved_other_lane_slots, max_new_jobs)]
    return ReplenishmentProposal(selected, state, count, count + len(selected),
                                 snapshot.evidence_sha256,
                                 reserved_other_lane_slots=reserved_other_lane_slots)
