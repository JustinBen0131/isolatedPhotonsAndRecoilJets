"""Pure, opt-in repair primitives for retained schema14 rows.

This module performs no I/O, changes no native producer, and applies no photon
class/isolation or jet matching prescription. A caller must supply the selected
populations separately and persist the returned unknown states and rebuild list.
Input rows use the native RJ*V1 branch names, with one explicitly named source
scope per call (compact IDs from separate ROOT parts must not be mixed).

Photon authority: PPG12 CURRENT_REFERENCE.json -> commit 57cd99850e888d31dc52559fc2377f94662bc739,
CaloAna24.cc:1358 and RecoEffCalculator.C:1983,2035-2054. The dominant primary
comes from the cluster evaluator, not a photon-only search or nearest deltaR.
JES authority: Fun4All_recoilJets_unified_impl.C:5714-5722,5830-5836,5925-5940;
RJDominantTruthWitnessV1.h; RJReplayFoundationV1.h jetCalibRawOrdinal contract.
No calibration coefficients or inferred no-calibration defaults live here.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import math
import operator
from typing import Any, Callable, Collection, Mapping, Sequence

Row = Mapping[str, Any]
Identity = tuple[int, int]
TruthKey = tuple[Identity, int, int]


class RepairError(ValueError):
    """An identity, payload, applicability, or provenance gate failed."""


def _integer(value: Any, name: str) -> int:
    if isinstance(value, bool):
        raise RepairError(f"{name} must be an integer, not bool")
    try:
        return operator.index(value)
    except TypeError as exc:
        raise RepairError(f"{name} must be an integer") from exc


def _id(row: Row, prefix: str) -> Identity:
    try:
        key = tuple(_integer(row[f"{prefix}_{part}"], prefix) for part in ("hi", "lo"))
    except KeyError as exc:
        raise RepairError(f"missing {prefix} identity") from exc
    if min(key) < 0 or max(key) >= 2**64 or key == (0, 0):
        raise RepairError(f"invalid {prefix} identity: {key}")
    return key  # type: ignore[return-value]


def _number(value: Any, name: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise RepairError(f"{name} must be finite") from exc
    if isinstance(value, bool) or not math.isfinite(result):
        raise RepairError(f"{name} must be finite")
    return result


def _unique(rows: Sequence[Row], prefix: str) -> dict[Identity, Row]:
    indexed: dict[Identity, Row] = {}
    for row in rows:
        identity = _id(row, prefix)
        _id(row, "event_id")
        if identity in indexed:
            raise RepairError(f"duplicate {prefix}: {identity}")
        indexed[identity] = row
    return indexed


def _selection(selected: Collection[Identity] | None, index: Mapping[Identity, Row],
               name: str) -> set[Identity]:
    result = set(index) if selected is None else set(selected)
    if result - index.keys():
        raise RepairError(f"{name} contains IDs absent from the retained inventory")
    return result


PHOTON_DONOR_INPUTS = (
    "source scope and event/cluster identities",
    "CaloRawClusterEval maximum-energy primary state and evaluator mode",
    "dominant GEANT track ID, embedding ID, PID and energy contribution",
    "truth-photon event/embedding/GEANT-track inventory, including rejected selections",
    "complete reconstructed-cluster capture receipt for truth-miss assertions",
    "complete truth-photon inventory receipt for response denominator acceptance",
    "explicit event inventory including events with no retained candidates or truth photons",
)


@dataclass(frozen=True)
class PhotonLinkRepair:
    source_scope: str
    associations: tuple[dict[str, Any], ...]
    links: tuple[dict[str, Any], ...]
    required_donor_inputs: tuple[str, ...]
    association_complete: bool
    response_complete: bool


def rebuild_photon_links(
    candidates: Sequence[Row], truth_photons: Sequence[Row], *,
    source_scope: str, simulation: bool,
    selected_candidate_ids: Collection[Identity] | None = None,
    selected_truth_ids: Collection[Identity] | None = None,
    complete_reco_events: Collection[Identity] = (),
    complete_truth_events: Collection[Identity] = (),
    event_ids: Collection[Identity] | None = None,
) -> PhotonLinkRepair:
    """Join energy-dominant identities and rebuild selection-relative links.

    ``None`` selections mean every supplied row; prompt class, isolation,
    embedding population and acceptance must be selected by the caller. Every
    cluster may independently match the same truth photon. Truth misses require
    an explicit complete reco inventory and no selected UNKNOWN cluster in that
    event. Response completeness also needs an explicit event inventory and
    complete truth/reco inventories for every event, including zero-object events.
    Missing truth donors, unavailable/invalid evaluators
    and incomplete reco capture remain UNKNOWN, never ordinary fakes or misses.

    Returned association IDs and link IDs are tuples; the writer must allocate
    fresh deterministic link IDs when serializing. Existing jet links are out of
    scope. A VALID_PRIMARY witness certifies the native evaluator's maximum; this
    function cannot reconstruct that maximum from a photon-only truth list.
    """
    if not isinstance(source_scope, str) or not source_scope.strip():
        raise RepairError("one explicit source_scope is required")
    if type(simulation) is not bool:
        raise RepairError("simulation must be explicit bool")
    reco = _unique(candidates, "candidate_id")
    selected_reco = _selection(selected_candidate_ids, reco, "reco selection")
    if not simulation:
        associations = tuple(dict(candidate_id=key, event_id=_id(row, "event_id"),
                                  state="NOT_APPLICABLE", reason="DATA", truth_id=None)
                             for key, row in reco.items())
        return PhotonLinkRepair(source_scope, associations, (), (), True, True)

    truth = _unique(truth_photons, "truth_photon_id")
    row_events = {_id(row, "event_id") for row in (*candidates, *truth_photons)}
    if event_ids is not None and row_events - set(event_ids):
        raise RepairError("object event absent from the explicit event inventory")
    selected_truth = _selection(selected_truth_ids, truth, "truth selection")
    truth_by_key: dict[TruthKey, Identity] = {}
    for key, row in truth.items():
        try:
            native = _integer(row["native_track_id"], "native_track_id")
            embed = _integer(row["embedding_id"], "embedding_id")
            valid = _integer(row["g4_photon_valid"], "g4_photon_valid")
        except KeyError as exc:
            raise RepairError(f"truth inventory missing donor field: {exc.args[0]}") from exc
        if native <= 0 or valid != 1:
            raise RepairError("truth inventory needs valid G4 photon and GEANT track ID")
        join = (_id(row, "event_id"), embed, native)
        if join in truth_by_key:
            raise RepairError(f"ambiguous event/embedding/GEANT-track truth identity: {join}")
        truth_by_key[join] = key

    associations: list[dict[str, Any]] = []
    links: list[dict[str, Any]] = []
    matched: set[Identity] = set()
    unknown_events: set[Identity] = set()
    missing: set[str] = set()
    for key, row in reco.items():
        event = _id(row, "event_id")
        association: dict[str, Any] = dict(candidate_id=key, event_id=event,
                                           state="UNKNOWN", reason="", truth_id=None)
        try:
            state = _integer(row["dominant_truth_state"], "dominant_truth_state")
            mode = _integer(row["dominant_truth_evaluator_mode"], "dominant_truth_evaluator_mode")
            if state not in (1, 2) or mode not in (1, 2):
                association["reason"] = "EVALUATOR_UNAVAILABLE_OR_INVALID_PRIMARY"
                missing.add(PHOTON_DONOR_INPUTS[1])
            elif state == 1:
                association.update(state="NO_PRIMARY", reason="EVALUATOR_PROVED_NO_PRIMARY")
            else:
                track = _integer(row["dominant_truth_track_id"], "dominant_truth_track_id")
                embed = _integer(row["dominant_truth_embedding_id"], "dominant_truth_embedding_id")
                pid = _integer(row["dominant_truth_pid"], "dominant_truth_pid")
                energy = _number(row["dominant_truth_energy_contribution"], "energy contribution")
                if track <= 0 or pid == 0 or energy < 0:
                    raise RepairError("invalid dominant-primary witness")
                association.update(embedding_id=embed, native_track_id=track,
                                   energy_contribution=energy)
                if pid != 22:
                    association.update(state="NON_PHOTON_PRIMARY", reason="DOMINANT_PRIMARY_NOT_PHOTON")
                else:
                    truth_id = truth_by_key.get((event, embed, track))
                    if truth_id is None:
                        association["reason"] = "DOMINANT_PHOTON_ABSENT_FROM_TRUTH_INVENTORY"
                        missing.add(PHOTON_DONOR_INPUTS[3])
                    else:
                        association.update(state="MATCH", reason="DOMINANT_PRIMARY_IDENTITY",
                                           truth_id=truth_id)
        except (KeyError, RepairError) as exc:
            association.update(state="UNKNOWN", reason=f"MISSING_OR_INVALID_WITNESS: {exc}")
            missing.add(PHOTON_DONOR_INPUTS[2])
        associations.append(association)
        if key not in selected_reco:
            continue
        link = dict(association)
        if association["state"] == "UNKNOWN":
            link["link_class"] = "UNKNOWN"
            unknown_events.add(event)
        elif association["state"] == "MATCH" and association["truth_id"] in selected_truth:
            link["link_class"] = "MATCH"
            matched.add(association["truth_id"])
        else:
            link["link_class"] = "RECO_FAKE"
            if association["state"] == "MATCH":
                link["reason"] = "DOMINANT_PHOTON_OUTSIDE_TRUTH_SELECTION"
                link["association_truth_id"] = association["truth_id"]
                link["truth_id"] = None
        links.append(link)

    complete_events = set(complete_reco_events)
    for key, row in truth.items():
        if key not in selected_truth or key in matched:
            continue
        event = _id(row, "event_id")
        known = event in complete_events and event not in unknown_events
        reason = "NO_SELECTED_RECO_ASSOCIATION" if known else "RECO_ASSOCIATION_OR_CAPTURE_INCOMPLETE"
        links.append(dict(candidate_id=None, event_id=event, truth_id=key,
                          link_class="TRUTH_MISS" if known else "UNKNOWN", reason=reason))
        if event not in complete_events:
            missing.add(PHOTON_DONOR_INPUTS[4])
    required_events = row_events if event_ids is None else set(event_ids)
    reco_complete = required_events <= complete_events
    truth_complete = required_events <= set(complete_truth_events)
    if not reco_complete:
        missing.add(PHOTON_DONOR_INPUTS[4])
    if not truth_complete:
        missing.add(PHOTON_DONOR_INPUTS[5])
    if event_ids is None:
        missing.add(PHOTON_DONOR_INPUTS[6])
    return PhotonLinkRepair(source_scope, tuple(associations), tuple(links), tuple(sorted(missing)),
                            all(a["state"] != "UNKNOWN" for a in associations),
                            event_ids is not None and reco_complete and truth_complete and
                            all(link["link_class"] != "UNKNOWN" for link in links))


# These consumers contain copied/derived values or selections. Updating only
# RJJetV1.corrected_pt never completes a JES repair. The jet matching algorithm
# remains unchanged; pT-dependent ordering/selection and its products need replay.
JES_REBUILD_REQUIRED = (
    "RJJetV1 corrected_pt consistency; preserve native deterministic_order container ordinal",
    "RJPhotonJetPairV1 xjgamma; preserve native jet_rank ordinal; replay any downstream pT-dependent selection",
    "RJRecoTruthLinkV1 jet selection/order-dependent bookkeeping using unchanged jet matching",
    "PhotonJetTree jet_pt, xjgamma and selected pair rows",
    "PhotonJetEventTree jet/pair vectors, leading jets, counts and event acceptance",
    "reco-jet threshold/fiducial/recoil selection and efficiency/fake/miss denominators",
    "jet/photon-jet response matrices, JES/JER diagnostics, alpha and unfolding inputs",
    "all jet-dependent histograms, closure results, cached summaries and publication receipts",
)


@dataclass(frozen=True)
class JESPayload:
    """Caller-supplied, byte-verified calibration plus exact applicability.

    The adapter must bind ``evaluate_corrected_pt(raw_pt, eta, vertex_z)`` to
    this payload and the stated function/branch, including real JetCalib eta/z
    routing. Its return value is corrected pT, never a multiplicative factor.
    The pure primitive verifies bytes and metadata; validating the ROOT adapter
    against real JetCalib remains an integration gate.
    """
    payload_bytes: bytes
    sha256: str
    function_name: str
    branch: str
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
    evaluate_corrected_pt: Callable[[float, float, float], float]


@dataclass(frozen=True)
class JetPTRepair:
    jet: dict[str, Any]
    status: str
    raw_pt_semantics: str
    corrected_pt_semantics: str
    rebuild_required: tuple[str, ...]
    capture_support: str
    provenance: dict[str, Any]


def _in_range(value: float, bounds: tuple[float, float], name: str) -> None:
    if len(bounds) != 2:
        raise RepairError(f"{name} applicability needs two finite bounds")
    low, high = (_number(item, name) for item in bounds)
    if low > high or not low <= value <= high:
        raise RepairError(f"{name} outside calibration applicability: {value}, {bounds}")


def repair_jet_pt(
    jet: Row, *, object_kind: str, context: Row | None = None,
    pre_jes_provenance: Row | None = None, payload: JESPayload | None = None,
) -> JetPTRepair:
    """Evaluate calibrated pT exactly once from explicitly proved pre-JES pT.

    Context requires source_scope, system, simulation, cdb_globaltag, timestamp,
    eta and vertex_z. Provenance must bind source_scope, event_id, jet_id,
    raw_pt, raw_pt_stage, raw_node, corrected_node and evidence_ref. It cannot be
    inferred from equality of old raw/corrected branches or a tree name.
    Truth jets are returned unchanged. No mass/four-vector/pair/selection repair
    is implied by CALIBRATED_PT_ONLY; the complete rebuild list must be discharged
    before export/acceptance. Uncaptured objects are never recovered here.
    """
    if object_kind == "truth":
        return JetPTRepair(dict(jet), "UNMODIFIED_TRUTH", "TRUTH_UNCHANGED", "TRUTH_UNCHANGED",
                           (), "TRUTH_UNCHANGED", {})
    if object_kind != "reco":
        raise RepairError("object_kind must be explicitly reco or truth")
    if context is None or pre_jes_provenance is None or payload is None:
        raise RepairError("explicit context, pre-JES provenance and calibration payload are required")
    try:
        event, identity = _id(jet, "event_id"), _id(jet, "jet_id")
        system = context["system"]
        if system not in ("pp", "auau") or type(context["simulation"]) is not bool:
            raise RepairError("explicit pp/auau system and simulation bool are required")
        scope = context["source_scope"]
        if not isinstance(scope, str) or not scope.strip():
            raise RepairError("source_scope is required")
        raw = _number(jet["raw_pt"], "raw_pt")
        eta = _number(context["eta"], "eta")
        vertex = _number(context["vertex_z"], "vertex_z")
        timestamp = _integer(context["timestamp"], "timestamp")
        radius = _number(jet["radius"], "radius")
        if raw < 0 or radius <= 0 or eta != _number(jet["eta"], "jet eta"):
            raise RepairError("invalid jet input or context eta mismatch")
        stage = "PRE_JES" if system == "pp" else "POST_UE_PRE_JES"
        subtraction = "none" if system == "pp" else "SUB1"
        if jet["subtraction_identity"] != subtraction:
            raise RepairError("unsupported nominal subtraction view for collision system")
        for field, expected in (("source_scope", scope), ("event_id", event),
                                ("jet_id", identity), ("raw_pt_stage", stage)):
            if pre_jes_provenance[field] != expected:
                raise RepairError(f"pre-JES provenance mismatch: {field}")
        if _number(pre_jes_provenance["raw_pt"], "provenance raw_pt") != raw:
            raise RepairError("pre-JES provenance raw_pt mismatch")
        for name in ("raw_node", "corrected_node", "evidence_ref"):
            if not isinstance(pre_jes_provenance[name], str) or not pre_jes_provenance[name].strip():
                raise RepairError(f"pre-JES provenance missing {name}")
        if pre_jes_provenance["raw_node"] == pre_jes_provenance["corrected_node"]:
            raise RepairError("pre-JES and corrected node identities must be distinct")
        if not payload.payload_bytes or hashlib.sha256(payload.payload_bytes).hexdigest() != payload.sha256:
            raise RepairError("calibration payload bytes/hash mismatch")
        for name in ("function_name", "branch", "cdb_globaltag"):
            if not isinstance(getattr(payload, name), str) or not getattr(payload, name).strip():
                raise RepairError(f"calibration missing {name}")
        for name, expected in (("system", system), ("simulation", context["simulation"]),
                               ("radius", radius), ("input_identity", jet["input_identity"]),
                               ("subtraction_identity", subtraction),
                               ("cdb_globaltag", context["cdb_globaltag"])):
            if getattr(payload, name) != expected:
                raise RepairError(f"calibration applicability mismatch: {name}")
        if type(payload.simulation) is not bool:
            raise RepairError("payload simulation applicability must be explicit bool")
        if not isinstance(jet["input_identity"], str) or not jet["input_identity"].strip():
            raise RepairError("jet input identity is required")
        _in_range(timestamp, payload.timestamp_range, "timestamp")
        _in_range(raw, payload.raw_pt_range, "raw_pt")
        _in_range(eta, payload.eta_range, "eta")
        _in_range(vertex, payload.vertex_z_range, "vertex_z")
    except KeyError as exc:
        raise RepairError(f"missing required JES field: {exc.args[0]}") from exc
    try:
        corrected = _number(payload.evaluate_corrected_pt(raw, eta, vertex), "calibrated pT")
    except Exception as exc:
        raise RepairError(f"calibration evaluation failed: {exc}") from exc
    if corrected < 0:
        raise RepairError("calibrated pT must be nonnegative")
    result = dict(jet)
    result["corrected_pt"] = corrected
    return JetPTRepair(
        result, "CALIBRATED_PT_ONLY", stage,
        "POST_JES" if system == "pp" else "POST_UE_POST_JES", JES_REBUILD_REQUIRED,
        "RETAINED_ROWS_ONLY: capture-threshold losses and unrecorded jets are not recovered; "
        "full-population acceptance requires capture support or donor reconstruction",
        dict(source_scope=scope, payload_sha256=payload.sha256, function_name=payload.function_name,
             branch=payload.branch, cdb_globaltag=payload.cdb_globaltag, timestamp=timestamp,
             pre_jes_evidence_ref=pre_jes_provenance["evidence_ref"],
             evaluation_semantics="TF1_EVAL_RETURNS_CORRECTED_PT"),
    )


@dataclass(frozen=True)
class JetPairRepair:
    pairs: tuple[dict[str, Any], ...]
    source_scope: str
    status: str = "RETAINED_PAIR_XJ_ONLY"
    scientific_acceptance: bool = False
    rebuild_required: tuple[str, ...] = JES_REBUILD_REQUIRED


def repair_retained_jet_pairs(pairs: Sequence[Row], photons: Sequence[Row],
                              jets: Sequence[JetPTRepair], *, source_scope: str) -> JetPairRepair:
    """Propagate evaluated JES pT into retained pairs without changing selection.

    Source: pp/AuAu writeReplayFoundationEvent stores xjgamma=corrected_pt/pt,
    and jet_rank=deterministic_order (container ordinal, NOT a sorted pT rank).
    Join by part/event/object identity, never entry position. Zero/nonpositive
    photon ET yields NaN as in the native producer. No new pairs, jets, mass
    correction, response selection or physics acceptance is inferred here.
    Callers must bind the input rows to this source scope before invoking it.
    """
    if not isinstance(source_scope, str) or not source_scope.strip():
        raise RepairError("explicit source_scope required")
    photon_index = _unique(photons, "candidate_id")
    _unique(pairs, "pair_id")
    jet_index = {}
    for repair in jets:
        if not isinstance(repair, JetPTRepair) or repair.status not in (
                "CALIBRATED_PT_ONLY", "JES_INVALID_NATIVE_SIGN_PRESERVED"):
            raise RepairError("pairs require classified reconstructed-jet JES repairs")
        if repair.provenance.get("source_scope") != source_scope:
            raise RepairError("jet repair source scope mismatch")
        row = repair.jet
        key = _id(row, "jet_id")
        _id(row, "event_id")
        if key in jet_index:
            raise RepairError("duplicate repaired jet identity")
        _number(row["corrected_pt"], "corrected jet pT")
        if row["corrected_pt"] < 0 or _integer(row["deterministic_order"], "jet ordinal") < 0:
            raise RepairError("invalid corrected jet pT or ordinal")
        jet_index[key] = repair
    result, seen_links = [], set()
    for pair in pairs:
        event = _id(pair, "event_id")
        candidate, jet = _id(pair, "candidate_id"), _id(pair, "jet_id")
        if candidate not in photon_index or jet not in jet_index:
            raise RepairError("pair references a missing photon or repaired jet")
        photon, jet_repair = photon_index[candidate], jet_index[jet]
        corrected = jet_repair.jet
        if _id(photon, "event_id") != event or _id(corrected, "event_id") != event:
            raise RepairError("pair crosses event identity")
        if (candidate, jet) in seen_links:
            raise RepairError("duplicate photon-jet pair relation")
        seen_links.add((candidate, jet))
        if _integer(pair["jet_rank"], "pair jet_rank") != corrected["deterministic_order"]:
            raise RepairError("pair jet_rank contradicts retained container ordinal")
        et = _number(photon["cluster_et"], "photon cluster_et")
        updated = dict(pair)
        valid_jes = jet_repair.status == "CALIBRATED_PT_ONLY"
        updated["xjgamma"] = corrected["corrected_pt"] / et if valid_jes and et > 0 else math.nan
        if valid_jes and et > 0 and not math.isfinite(updated["xjgamma"]):
            raise RepairError("repaired xjgamma is not finite")
        result.append(updated)
    return JetPairRepair(tuple(result), source_scope)
