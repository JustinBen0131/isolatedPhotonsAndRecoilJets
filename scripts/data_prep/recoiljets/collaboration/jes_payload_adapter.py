"""Opt-in adapter for one historical JetCalib source, never current-runtime parity.

Bound source SHA256: a327a8dcafad5c7602943201705983374102a1a41fd6d9d50841a37103e04dae.
This implements that source's z+eta-enabled branch for R=0.2/0.3/0.4.
Runtime identity, disabled EM-fraction correction and externally approved input
applicability are mandatory. No payload, coefficients or scientific approval are
supplied by this module. Injection is restricted to contracts marked TEST_ONLY.

Usage: adapter = HistoricalJetCalibAdapter(payload_bytes, contract, evaluator)
Then pass adapter.payload_for(eta, vertex_z) to physics_repairs.repair_jet_pt.
Persist adapter.route(eta, vertex_z) with each result for z/eta provenance. Use
evaluate_with_provenance for a standalone diagnostic. Calibration alone does not
rebuild downstream products or certify capture support.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
import math
from pathlib import Path
import re
import struct
from typing import Any, Callable, Mapping, Protocol

try:
    from .physics_repairs import JESPayload, RepairError
except ImportError:
    from physics_repairs import JESPayload, RepairError

SNAPSHOT_SHA256 = "a327a8dcafad5c7602943201705983374102a1a41fd6d9d50841a37103e04dae"
HISTORICAL_BRANCH = "HISTORICAL_JETCALIB_Z_AND_ETA_EMFRAC_DISABLED"


def float32(value: float) -> float:
    """C++ float assignment, rejecting nonfinite/unrepresentable input."""
    try:
        result = struct.unpack("=f", struct.pack("=f", float(value)))[0]
    except (TypeError, ValueError, OverflowError, struct.error) as exc:
        raise RepairError("value is not representable as finite float32") from exc
    if isinstance(value, bool) or not math.isfinite(result):
        raise RepairError("value is not representable as finite float32")
    return result


def _radius(value: float) -> float:
    result = float32(value)
    if result not in tuple(float32(item) for item in (.2, .3, .4)):
        raise RepairError("historical adapter supports only R=0.2/0.3/0.4")
    return result


@dataclass(frozen=True)
class Route:
    function_name: str
    z_bin: int
    eta_bin: int
    uses_eta_suffixed_function: bool
    radius_float32: float
    eta_float32: float
    vertex_z_float32: float
    missing_vertex_sentinel: bool
    historical_z_default_fallback: bool
    source_sha256: str = SNAPSHOT_SHA256


def historical_route(radius: float, eta: float, vertex_z: float) -> Route:
    """Route only; this grants no payload/application authority.

    Float32 arithmetic deliberately preserves the R=0.4 quirk: the stored
    float radius is greater than double literal 0.4, so InitRun loads Z-only
    functions into every eta slot. Do not 'fix' it in a historical adapter.
    """
    radius, eta, z = _radius(radius), float32(eta), float32(vertex_z)
    fallback = False
    if -60. <= z < -30.:
        z_bin = 1
    elif -30. <= z < 30.:
        z_bin = 0
    elif 30. <= z < 60.:
        z_bin = 2
    elif -900. <= z < -60. or 60. <= z < 900.:
        z_bin = 3
    elif z < -900.:
        z_bin = 4
    else:
        z_bin, fallback = 0, True
    eta_bin = 0
    if z_bin < 3:
        edges = ((-1.2, 1.2), (-.95, 1.25), (-1.25, .95))[z_bin]
        # Literals in C++ are double, then stored to float; high-low is float.
        low, high = float32(edges[0] + radius), float32(edges[1] - radius)
        span = float32(high - low)
        thresholds = (float32(low + span / 4.), float32(low + span / 2.),
                      float32(low + 3. * span / 4.))
        eta_bin = sum(eta >= threshold for threshold in thresholds)
    radius_code = int(float32(radius * 10))
    name = f"JES_Calib_Func_R0{radius_code}_Z{z_bin}"
    use_eta = z_bin < 3 and radius <= .4
    if use_eta:
        name += f"_Eta{eta_bin}"
    return Route(name, z_bin, eta_bin, use_eta, radius, eta, z, z_bin == 4, fallback)


def required_function_names(radius: float) -> tuple[str, ...]:
    """Exact InitRun inventory, including functions outside a selected z window."""
    r = _radius(radius)
    base = f"JES_Calib_Func_R0{int(float32(r * 10))}"
    return tuple(f"{base}_Z{z}" + (f"_Eta{eta}" if z < 3 and r <= .4 else "")
                 for z in range(5) for eta in (range(4) if z < 3 and r <= .4 else (0,)))


def _bounds(bounds: tuple[float, float], name: str) -> tuple[float, float]:
    if len(bounds) != 2:
        raise RepairError(f"{name} requires an explicit two-sided domain")
    low, high = (float(value) for value in bounds)
    if not math.isfinite(low) or not math.isfinite(high) or low > high:
        raise RepairError(f"invalid {name} domain")
    return low, high


def _within(value: float, bounds: tuple[float, float], name: str) -> None:
    low, high = _bounds(bounds, name)
    if not low <= float(value) <= high:
        raise RepairError(f"{name} outside approved applicability domain")


@dataclass(frozen=True)
class ApplicabilityContract:
    """Exact external decision; strings/flags are references, not proof by themselves.

    The caller must verify the cited approval and runtime evidence before loading
    this contract. This module checks its explicit bindings and never fabricates
    an approval from payload existence, matching labels or equal jet branches.
    Only HISTORICAL_REPLAY or TEST_ONLY is supported; neither claims production
    acceptance, a current JetCalib implementation, or fullseg/retower equivalence.
    """
    contract_id: str
    approval_status: str
    approved_by: str
    approval_evidence: str
    applicability_evidence: str
    runtime_evidence: str
    runtime_source_sha256: str
    em_fraction_state: str
    purpose: str
    payload_sha256: str
    cdb_globaltag: str
    system: str
    simulation: bool
    radius: float
    input_identity: str
    subtraction_identity: str
    timestamp_range: tuple[int, int]
    raw_pt_range: tuple[float, float]
    eta_range: tuple[float, float]
    vertex_z_range: tuple[float, float]

    def validate(self) -> None:
        if self.approval_status != "APPROVED":
            raise RepairError("external approved applicability contract is required")
        for name in ("contract_id", "approved_by", "approval_evidence", "applicability_evidence",
                     "runtime_evidence", "cdb_globaltag", "input_identity"):
            if not isinstance(getattr(self, name), str) or not getattr(self, name).strip():
                raise RepairError(f"explicit contract field required: {name}")
        if self.runtime_source_sha256 != SNAPSHOT_SHA256:
            raise RepairError("runtime is unknown or differs from the bound historical snapshot")
        if self.em_fraction_state != "DISABLED":
            raise RepairError("EM-fraction correction must be explicitly proved DISABLED")
        if self.purpose not in ("HISTORICAL_REPLAY", "TEST_ONLY"):
            raise RepairError("only historical replay or test purpose is supported")
        if not re.fullmatch(r"[0-9a-f]{64}", self.payload_sha256):
            raise RepairError("actual payload SHA256 is required")
        _radius(self.radius)
        if self.system not in ("pp", "auau") or type(self.simulation) is not bool:
            raise RepairError("explicit collision system and DATA/MC applicability required")
        expected = "none" if self.system == "pp" else "SUB1"
        if self.subtraction_identity != expected:
            raise RepairError("adapter supports pp unsubtracted or nominal AuAu SUB1 only")
        for name in ("timestamp_range", "raw_pt_range", "eta_range", "vertex_z_range"):
            _bounds(getattr(self, name), name)
        if self.raw_pt_range[0] < 0:
            raise RepairError("negative pre-JES pT is outside this adapter")
        if any(type(value) is not int for value in self.timestamp_range):
            raise RepairError("timestamp applicability requires integer bounds")

    @property
    def sha256(self) -> str:
        return hashlib.sha256(json.dumps(asdict(self), sort_keys=True, separators=(",", ":"),
                                         allow_nan=False).encode()).hexdigest()


class Evaluator(Protocol):
    payload_sha256: str
    kind: str

    def domain(self, function_name: str) -> tuple[float, float]: ...
    def evaluate(self, function_name: str, raw_pt_float32: float) -> float: ...


@dataclass(frozen=True)
class InjectedTestEvaluator:
    """Toy evaluator; rejected by any non-TEST_ONLY applicability contract."""
    payload_sha256: str
    domains: Mapping[str, tuple[float, float]]
    callback: Callable[[str, float], float]
    kind: str = "INJECTED_TEST"

    def domain(self, function_name: str) -> tuple[float, float]:
        try:
            return self.domains[function_name]
        except KeyError as exc:
            raise RepairError(f"missing TF1: {function_name}") from exc

    def evaluate(self, function_name: str, raw_pt_float32: float) -> float:
        return self.callback(function_name, raw_pt_float32)


class LocalRootTF1Evaluator:
    """Read a caller-supplied local ROOT file; never resolve CDB or fetch remotely.

    Use the repository ROOT environment wrapper. Keep this object alive while
    using adapters; close explicitly when finished. Payload hashes are checked
    before and after opening; exact named objects must be TF1 instances.
    """
    kind = "LOCAL_ROOT_TF1"

    def __init__(self, path: str | Path, expected_sha256: str):
        if "://" in str(path):
            raise RepairError("only an existing local payload file is permitted")
        self.path = Path(path).resolve(strict=True)
        if not self.path.is_file():
            raise RepairError("payload must be a local regular file")
        self.payload_bytes = self.path.read_bytes()
        self.payload_sha256 = hashlib.sha256(self.payload_bytes).hexdigest()
        if self.payload_sha256 != expected_sha256:
            raise RepairError("local ROOT payload hash mismatch")
        import ROOT  # Deliberately lazy; routing and injected tests need no ROOT.
        self._file = ROOT.TFile.Open(str(self.path), "READ")
        if not self._file or self._file.IsZombie():
            raise RepairError("local ROOT payload cannot be opened")
        if hashlib.sha256(self.path.read_bytes()).hexdigest() != expected_sha256:
            self.close()
            raise RepairError("payload changed while opening ROOT file")

    def _function(self, name: str) -> Any:
        if not self._file or not self._file.IsOpen():
            raise RepairError("local ROOT payload is closed")
        function = self._file.Get(name)
        if not function or not function.InheritsFrom("TF1"):
            raise RepairError(f"missing/non-TF1 payload object: {name}")
        return function

    def domain(self, function_name: str) -> tuple[float, float]:
        function = self._function(function_name)
        return float(function.GetXmin()), float(function.GetXmax())

    def evaluate(self, function_name: str, raw_pt_float32: float) -> float:
        return float(self._function(function_name).Eval(raw_pt_float32))

    def close(self) -> None:
        if self._file:
            self._file.Close()


@dataclass(frozen=True)
class Evaluation:
    corrected_pt: float
    raw_pt_float32: float
    route: Route
    payload_sha256: str
    applicability_contract_sha256: str
    purpose: str
    status: str = "HISTORICAL_CALIBRATED_PT_ONLY"


class HistoricalJetCalibAdapter:
    def __init__(self, payload_bytes: bytes, contract: ApplicabilityContract, evaluator: Evaluator):
        contract.validate()
        actual = hashlib.sha256(payload_bytes).hexdigest()
        if not payload_bytes or actual != contract.payload_sha256 or actual != evaluator.payload_sha256:
            raise RepairError("actual payload bytes, contract and evaluator hashes must agree")
        if isinstance(evaluator, InjectedTestEvaluator):
            if evaluator.kind != "INJECTED_TEST" or contract.purpose != "TEST_ONLY":
                raise RepairError("injected evaluator requires TEST_ONLY applicability and kind")
        elif isinstance(evaluator, LocalRootTF1Evaluator):
            if evaluator.kind != "LOCAL_ROOT_TF1":
                raise RepairError("unrecognized local ROOT evaluator kind")
        else:
            raise RepairError("unrecognized payload evaluator")
        for name in required_function_names(contract.radius):
            _bounds(evaluator.domain(name), f"TF1 {name}")
        self.payload_bytes, self.contract, self.evaluator = payload_bytes, contract, evaluator

    def route(self, eta: float, vertex_z: float) -> Route:
        for value, domain, name in ((eta, self.contract.eta_range, "eta"),
                                    (vertex_z, self.contract.vertex_z_range, "vertex_z")):
            _within(value, domain, name)
            _within(float32(value), domain, f"float32 {name}")
        return historical_route(self.contract.radius, eta, vertex_z)

    def evaluate_with_provenance(self, raw_pt: float, eta: float, vertex_z: float) -> Evaluation:
        route = self.route(eta, vertex_z)
        pt = float32(raw_pt)
        _within(raw_pt, self.contract.raw_pt_range, "raw_pt")
        _within(pt, self.contract.raw_pt_range, "float32 raw_pt")
        _within(pt, self.evaluator.domain(route.function_name), "TF1 evaluation pT")
        try:
            corrected = float32(self.evaluator.evaluate(route.function_name, pt))
        except Exception as exc:
            raise RepairError(f"historical TF1 evaluation failed: {exc}") from exc
        if corrected < 0:
            raise RepairError("historical TF1 returned negative corrected pT")
        return Evaluation(corrected, pt, route, self.contract.payload_sha256,
                          self.contract.sha256, self.contract.purpose)

    def payload_for(self, eta: float, vertex_z: float) -> JESPayload:
        """Existing repair API, bound to one exact eta/z point and function.

        Capture ``route(eta, vertex_z)`` alongside the repair receipt. Supplying
        a different float32 eta/z point to this callable is rejected even when
        it would choose the same TF1; build a fresh binding for the next jet.
        """
        route = self.route(eta, vertex_z)

        def evaluate(pt: float, actual_eta: float, actual_z: float) -> float:
            if self.route(actual_eta, actual_z) != route:
                raise RepairError("JESPayload eta/z binding changed")
            return self.evaluate_with_provenance(pt, actual_eta, actual_z).corrected_pt

        c = self.contract
        branch = f"{HISTORICAL_BRANCH}:{c.purpose}:contract={c.sha256}"
        return JESPayload(self.payload_bytes, c.payload_sha256, route.function_name,
                          branch, c.cdb_globaltag, c.system, c.simulation, c.radius,
                          c.input_identity, c.subtraction_identity, c.timestamp_range,
                          c.raw_pt_range, c.eta_range, c.vertex_z_range, evaluate)
