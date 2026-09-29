"""Decode additive native captures into the keyed augmentation contract.

Local ROOT bytes and joins are verified. Caller-supplied DST, reconstruction,
calibration and model bindings still require independent runtime acceptance.
No model is run, no selection is changed, and no donor is joined by row order.
"""
from dataclasses import asdict, dataclass, replace
import json
import math
from pathlib import Path
import struct

import uproot

import augmentation_contract as c
from canonical_augmentation_io import BASE, _read
from canonical_repairs_io import _identity, file_sha256

VERSION = "NativeAugmentationReaderV1"
ADDITIVE_CAPTURE_POLICY = 'SAME_INPUTS_AND_CALIBRATION_ADDITIVE_CAPTURE_V1'
FLAGS = ("reference_preselection_state", "active_preselection_state", "preselection_bitmask")
OPTIONAL_FLAGS = ("ppg12_tower_mask_state", "ppg12_source_eligible_state")
MBD_FIELDS = tuple("mbd_native_" + name for name in (
    "capture_version", "available", "is_valid", "t0_ns", "south_time_ns", "north_time_ns",
    "t0_finite", "south_time_finite", "north_time_finite", "south_npmt", "north_npmt"))
CALO_FIELDS = (*c.CALO_FIELDS, "event_calo_capture_version", "event_calo_available_mask",
               "event_calo_valid_mask", "event_calo_require_isgood")
TIMING_FIELDS = ("native_timing_capture_version", "native_mean_time_samples",
    "native_time_energy_denominator", "native_time_energy_numerator", "native_time_contributing_towers",
    "native_time_finite", "native_timing_input_cluster_node", "native_timing_calibrated_tower_node")
NPB_FIELDS = ("capture_version", "evaluation_state", "domain_mode", "in_producer_domain",
    "model_available", "input_order_is_npb25", "inputs_captured", "inputs_finite", "min_et", "max_et",
    "max_abs_eta", "raw_score", "source_cluster_et", "source_eta", "source_phi", "binding_policy",
    "source_identity_available", "source_native_cluster_key", "source_encounter_ordinal",
    "kinematic_input_bits_available", "input_et_bits", "input_eta_bits", "source_cluster_node",
    "model_file", "feature_names", "ordered_inputs", "ordered_input_finite")


@dataclass(frozen=True)
class DecodedAugmentation:
    report: c.OverlayReport
    receipt: dict


def _require(row, fields):
    absent = set(fields) - row.keys()
    if absent:
        raise c.ContractError(f"partial native capture missing fields: {sorted(absent)}")


def _flag(value, name):
    if type(value) is not int or value not in (0, 1):
        raise c.ContractError(f"invalid native boolean: {name}")
    return bool(value)


def _capture_present(row, version, fields):
    if version not in row or row[version] == 0:
        return False
    if type(row[version]) is not int or row[version] != 1:
        raise c.ContractError(f"unsupported native capture version: {version}")
    _require(row, fields)
    return True


def _finite_flag(value, flag, name):
    if math.isfinite(value) != _flag(flag, name):
        raise c.ContractError(f"native finite flag disagrees with value: {name}")


def _measurement(value, valid, reason):
    return c.Measurement("VALID", value) if valid else c.Measurement(
        "INVALID", value if math.isfinite(value) else None, reason)


def _float_bits(value):
    return struct.pack("!f", value).hex()


def decode_mbd(row, definition):
    if not _capture_present(row, "mbd_native_capture_version", MBD_FIELDS):
        return c.Missing("native MbdOut capture absent")
    if not _flag(row["mbd_native_available"], "mbd availability"):
        return c.Missing("native MbdOut node unavailable")
    valid = _flag(row["mbd_native_is_valid"], "MbdOut isValid")
    values = [row["mbd_native_"+name+"_ns"] for name in ("t0", "south_time", "north_time")]
    for name, value in zip(("t0", "south_time", "north_time"), values):
        _finite_flag(value, row["mbd_native_"+name+"_finite"], name)
    hits = [row["mbd_native_"+arm+"_npmt"] for arm in ("south", "north")]
    for value in hits:
        c._int(value, "native MBD arm hits")
    if valid and (not math.isfinite(values[0]) or values[0] <= -9999):
        raise c.ContractError("native MbdOut isValid contradicts its t0")
    measures = [_measurement(values[0], valid, "native MbdOut invalid")]
    measures.extend(_measurement(value, count > 0 and math.isfinite(value) and value > -9999,
        "empty MBD arm or invalid native arm time") for count, value in zip(hits, values[1:]))
    return c.MbdCapture(definition, valid, *measures, *hits)


def decode_calo(row, definition):
    if not _capture_present(row, "event_calo_capture_version", CALO_FIELDS):
        return c.Missing("native calorimeter capture absent")
    if _flag(row["event_calo_require_isgood"], "calorimeter isGood policy") != definition.require_isgood:
        raise c.ContractError("native calorimeter isGood policy differs from bound definition")
    available, valid = row["event_calo_available_mask"], row["event_calo_valid_mask"]
    c._int(available, "calorimeter availability", 0, 7)
    c._int(valid, "calorimeter validity", 0, 7)
    if valid & ~available:
        raise c.ContractError("calorimeter validity contradicts availability")
    values = []
    for index, name in enumerate(c.CALO_FIELDS):
        value = row[name]
        present = bool(available & (1 << index)) if index < 3 else available == 7
        complete = bool(valid & (1 << index)) if index < 3 else valid == 7
        if not present:
            if not math.isnan(value):
                raise c.ContractError("absent calorimeter container has a scalar value")
            values.append(c.Measurement("MISSING", reason="calibrated calorimeter container unavailable"))
            continue
        if index < 3 and complete and not math.isfinite(value):
            raise c.ContractError("valid calorimeter layer has a nonfinite value")
        if math.isfinite(value):
            try:
                stored = struct.unpack("!f", struct.pack("!f", value))[0]
            except (OverflowError, struct.error) as exc:
                raise c.ContractError("native calorimeter scalar exceeds binary32 storage") from exc
            if stored != value:
                raise c.ContractError("native calorimeter scalar is not the specified binary32 storage value")
        values.append(_measurement(value, complete and math.isfinite(value),
                                   "incomplete or nonfinite calibrated calorimeter sum"))
    return c.CaloCapture(definition, available, valid, tuple(values))


def decode_timing(row, definition):
    if not _capture_present(row, "native_timing_capture_version", TIMING_FIELDS):
        return c.Missing("native cluster-time capture absent")
    if (definition.definition != c.NATIVE_TIMING or
            row["native_timing_input_cluster_node"] != definition.cluster_node or
            row["native_timing_calibrated_tower_node"] != definition.tower_node):
        raise c.ContractError("native timing definition/node binding differs")
    raw, denom, numerator = (row[name] for name in ("native_mean_time_samples",
        "native_time_energy_denominator", "native_time_energy_numerator"))
    count = row["native_time_contributing_towers"]
    c._int(count, "native timing tower count")
    _finite_flag(raw, row["native_time_finite"], "cluster time")
    usable = math.isfinite(denom) and denom > 0 and count > 0
    if usable and math.isfinite(numerator):
        # The producer accumulates float numerator/denominator and stores the
        # float quotient. Compare that exact convention, not ideal doubles.
        quotient = numerator / denom
        try:
            expected = struct.unpack("!f", struct.pack("!f", quotient))[0]
        except OverflowError:
            expected = math.copysign(math.inf, quotient)
        if raw != expected:
            raise c.ContractError("native cluster-time quotient contradicts numerator/denominator")
    valid = usable and math.isfinite(numerator) and math.isfinite(raw) and raw != -999
    return c.TimingCapture(definition,
        _measurement(raw, valid, "invalid native mean time or empty tower denominator"),
        _measurement(denom, math.isfinite(denom), "nonfinite native timing denominator"), count)


def decode_npb(row, retained, definition, *, model_file, evidence):
    prefix = "auau_npb_capture_" if definition.score_field == "auau_npb_score" else "npb_capture_"
    fields = tuple(prefix+name for name in NPB_FIELDS)
    if not _capture_present(row, prefix+"capture_version", fields):
        return c.Missing("native NPB capture absent")
    r = {name: row[prefix+name] for name in NPB_FIELDS}
    state = r["evaluation_state"]
    if type(state) is not int or state not in range(1, 6):
        raise c.ContractError("invalid captured NPB evaluation state")
    if not model_file:
        return c.Missing("native NPB model path has no explicit runtime dependency binding")
    if r["model_file"] != model_file:
        raise c.ContractError("native NPB model path differs from bound runtime dependency")
    if r["feature_names"] != list(c.NPB_FEATURES) or r["input_order_is_npb25"] != 1:
        raise c.ContractError("native NPB input order differs from fixed 25-input definition")
    mode = {"NAMED_FLOAT32_V1": 1, "PRIMARY_UNGATED_V1": 2}.get(definition.domain_policy)
    if mode is None or r["domain_mode"] != mode:
        raise c.ContractError("native NPB domain mode differs from bound definition")
    if mode == 1:
        bounds = definition.native_bounds
        expected = [b.value for b in bounds] if bounds else [definition.et_min, definition.et_max, definition.abs_eta_max]
        for name, value in zip(("min_et", "max_et", "max_abs_eta"), expected):
            if _float_bits(r[name]) != _float_bits(value):
                raise c.ContractError(f"native NPB domain bound differs: {name}")
    cut_name = "auau_npb_configured_cut" if definition.score_field == "auau_npb_score" else "npb_configured_cut"
    _require(row, (cut_name,))
    if definition.configured_cut is None:
        if math.isfinite(row[cut_name]):
            raise c.ContractError("native NPB cut is finite but definition left it unbound")
    elif row[cut_name] != definition.configured_cut:
        raise c.ContractError("native NPB configured cut differs")
    captured = _flag(r["inputs_captured"], "NPB inputs captured")
    inputs_finite = _flag(r["inputs_finite"], "NPB inputs finite")
    available = _flag(r["model_available"], "NPB model available")
    inside = _flag(r["in_producer_domain"], "NPB native domain")
    if (not _flag(r["source_identity_available"], "NPB source identity") or
            not _flag(r["kinematic_input_bits_available"], "NPB accessor bits")):
        return c.Missing("NPB scoring-object identity or actual kinematic accessor bits unavailable")
    if r["binding_policy"] not in (1, 2):
        raise c.ContractError("unrecognized native NPB scoring-object binding")
    c._int(r["input_et_bits"], "NPB ET bits", 0, 2**32-1)
    c._int(r["input_eta_bits"], "NPB eta bits", 0, 2**32-1)
    if not all(math.isfinite(c._float32_bits_value(f'{r[name]:08x}'))
               for name in ("input_et_bits", "input_eta_bits")):
        return c.Missing("nonfinite native NPB kinematic inputs; raw capture retained, valid scoring binding unavailable")
    witness = c.ScoringObjectWitness(
        c.CandidateKey(retained.key.event, r["source_native_cluster_key"]), r["source_cluster_node"],
        c.CandidateAnchor(r["source_cluster_et"], r["source_eta"], r["source_phi"]),
        f'{r["input_et_bits"]:08x}', f'{r["input_eta_bits"]:08x}',
        "SAME_OBJECT_V1" if r["binding_policy"] == 1 else "NATIVE_FIRST_KINEMATIC_MATCH_V1",
        evidence, None if r["binding_policy"] == 1 else r["source_encounter_ordinal"])
    witness.validate_binding(retained.key, retained.anchor, definition)
    if inside != definition.contains(retained.key.event.source.sample, witness.input_anchor):
        raise c.ContractError("captured NPB domain decision contradicts bound native policy")
    values, finite_flags = r["ordered_inputs"], r["ordered_input_finite"]
    if captured:
        if len(values) != 25 or len(finite_flags) != 25:
            raise c.ContractError("captured NPB needs exactly 25 values and finite flags")
        for name, value, flag in zip(c.NPB_FEATURES, values, finite_flags):
            _finite_flag(value, flag, name)
        if inputs_finite != all(math.isfinite(v) for v in values):
            raise c.ContractError("NPB all-inputs-finite flag contradicts captured vector")
        inputs = tuple(c.FeatureInput(name, _measurement(value, math.isfinite(value),
            "nonfinite native NPB input"), "NATIVE_PRODUCER", evidence)
            for name, value in zip(c.NPB_FEATURES, values))
    else:
        if inputs_finite:
            raise c.ContractError("uncaptured NPB inputs cannot be all finite")
        inputs = tuple(c.FeatureInput(name, c.Measurement("MISSING", reason="native input not captured"),
            "UNAVAILABLE", None) for name in c.NPB_FEATURES)
    score = r["raw_score"]
    if ((state == 1 and available) or (state == 2 and (not available or inside)) or
            (state == 3 and (not available or not inside or inputs_finite))):
        raise c.ContractError("NPB evaluation state contradicts native invocation flags")
    if state in (4, 5) and (not available or not inside or (state == 4) != math.isfinite(score)):
        raise c.ContractError("NPB evaluated state contradicts native invocation/score")
    applicability = "APPLICABLE" if inside else "NOT_APPLICABLE"
    if state == 4 and captured and inputs_finite and 0 <= score <= 1:
        evaluation, measured = "EVALUATED", c.Measurement("VALID", score)
    elif state in (4, 5):
        evaluation = "EVALUATED_INVALID"
        measured = c.Measurement("INVALID", score if math.isfinite(score) else None,
                                 "native Compute ran with invalid inputs or invalid output")
    elif not inside:
        evaluation, measured = "NOT_APPLICABLE", c.Measurement("NOT_APPLICABLE",
            score if math.isfinite(score) else None, "outside the bound native NPB domain")
    else:
        evaluation = "NOT_EVALUATED" if inputs_finite else "INPUTS_UNAVAILABLE"
        measured = c.Measurement("NOT_EVALUATED", score if math.isfinite(score) else None,
                                 "native Compute did not run")
    result = c.NPBCapture(definition, inputs, applicability, evaluation, measured, witness)
    result.validate_candidate(retained.key, retained.anchor)
    return result


def _inventory(path, binding, max_rows):
    from canonical_donor_native_view import open_native_view
    with open_native_view(path) as root:
        metadata = root.get(BASE+"config_sha256")
        config_hash = metadata.member("fTitle") if metadata is not None else None
        sources = _read(root[BASE+"RJSourceOccurrenceV1"], max_rows,
            ["source_occurrence_id_hi", "source_occurrence_id_lo", "run", "segment", "source_manifest_sha256"],
            ("input_uri_hash",))
        event_tree, candidate_tree = root[BASE+"RJEventV1"], root[BASE+"RJPhotonCandidateV1"]
        events = _read(event_tree, max_rows, ["event_id_hi", "event_id_lo", "source_occurrence_id_hi",
            "source_occurrence_id_lo", "run", "physical_event_sequence", "physical_event_sequence_valid"],
            (*c.EVENT_FIELDS, *MBD_FIELDS, *CALO_FIELDS))
        optional = (*c.CANDIDATE_FIELDS, *OPTIONAL_FLAGS, *TIMING_FIELDS,
            "npb_configured_cut", "auau_npb_configured_cut",
            *(prefix+field for prefix in ("npb_capture_", "auau_npb_capture_") for field in NPB_FIELDS))
        candidates = _read(candidate_tree, max_rows, ["candidate_id_hi", "candidate_id_lo", "event_id_hi",
            "event_id_lo", "native_cluster_key", "cluster_et", "eta", "phi", *FLAGS], optional)
    by_source = {_identity(r, "source_occurrence_id"): r for r in sources}
    # This decoder's contract is one exact source occurrence per part. It must
    # not conflate repeated physical event numbers from different occurrences.
    if len(sources) != 1 or len(by_source) != 1:
        raise c.ContractError("native augmentation needs exactly one bound source occurrence")
    source = sources[0]
    if (source["run"], source["segment"], source["source_manifest_sha256"]) != (
            binding.run, binding.segment, binding.source_manifest_sha256):
        raise c.ContractError("native augmentation source/run/segment/manifest differs")
    if binding.input_tuple_sha256 is not None and source.get('input_uri_hash') != binding.input_tuple_sha256:
        raise c.ContractError('native augmentation frozen input tuple differs')
    ids, by_event = {}, {}
    for row in events:
        eid, physical = _identity(row, "event_id"), row["physical_event_sequence"]
        c._int(physical, "physical event", 0, 2**63-1)
        if (row["physical_event_sequence_valid"] != 1 or eid in ids or physical in by_event or
                _identity(row, "source_occurrence_id") not in by_source or row["run"] != binding.run):
            raise c.ContractError("ambiguous/unavailable native source/physical-event identity")
        ids[eid], by_event[physical] = physical, row
    candidate_ids, by_candidate = set(), {}
    for row in candidates:
        eid, cid = _identity(row, "event_id"), _identity(row, "candidate_id")
        if eid not in ids or cid in candidate_ids:
            raise c.ContractError("duplicate candidate or candidate event unavailable")
        key = (ids[eid], row["native_cluster_key"])
        if key in by_candidate:
            raise c.ContractError("duplicate native cluster key")
        candidate_ids.add(cid)
        by_candidate[key] = row
    return by_event, by_candidate, config_hash


def decode_native_augmentation(source, binding, *, timing, mbd, npb, model_file=None,
                               donor=None, donor_binding=None, retained_semantics=None, max_rows=2_000_000,
                               calo=None, donor_migration_policy=None, donor_event_superset=False):
    """Decode one part or an explicitly selected, exact-source donor subset.

    ``model_file`` binds a captured runtime path to the supplied NPB definition;
    it does not verify model bytes. ``retained_semantics`` is required for any
    finite overlap and must come from independently bound source provenance,
    never a branch-name guess. Native capture absence remains explicit missing.
    A changed producer needs an explicit additive migration policy; exact input
    files and calibration basis still cannot change. An explicitly allowed
    donor event superset is intersected with retained physical events. Extra
    candidate objects inside that intersection are always rejected.
    """
    if type(max_rows) is not int or max_rows <= 0:
        raise c.ContractError("positive bounded max_rows required")
    if type(donor_event_superset) is not bool:
        raise c.ContractError('donor event superset must be explicit boolean')
    if donor_migration_policy not in (None, ADDITIVE_CAPTURE_POLICY):
        raise c.ContractError('unknown donor migration policy')
    source = Path(source).resolve()
    donor = Path(donor).resolve() if donor is not None else source
    source_hash = file_sha256(source)
    if source_hash != binding.retained_file_sha256:
        raise c.ContractError("native augmentation retained-file bytes differ")
    donor_hash = source_hash if donor == source else file_sha256(donor)
    if donor != source:
        if donor_binding is None or donor_binding.retained_file_sha256 != donor_hash:
            raise c.ContractError('separate native donor needs an explicit source binding and exact donor bytes')
        expected = replace(binding, retained_file_sha256=donor_hash)
        if donor_migration_policy == ADDITIVE_CAPTURE_POLICY:
            expected = replace(expected, source_manifest_sha256=donor_binding.source_manifest_sha256,
                               reconstruction_sha256=donor_binding.reconstruction_sha256)
        if donor_binding != expected:
            raise c.ContractError('native donor changes source/dependency binding or calibration basis')
    elif donor_binding is not None and donor_binding != binding:
        raise c.ContractError("same-file donor binding differs")
    events, candidates, config_hash = _inventory(source, binding, max_rows)
    donor_events, donor_candidates, donor_config = ((events, candidates, config_hash) if donor == source
        else _inventory(donor, donor_binding, max_rows))
    excluded_donor_events = len(donor_events.keys() - events.keys())
    excluded_donor_candidates = 0
    if donor_event_superset:
        selected_candidates = {key:row for key,row in donor_candidates.items() if key[0] in events}
        excluded_donor_candidates = len(donor_candidates) - len(selected_candidates)
        donor_candidates = selected_candidates
        donor_events = {key:row for key,row in donor_events.items() if key in events}
    captured = (any(r.get("mbd_native_capture_version", 0) or
                   (calo is not None and r.get("event_calo_capture_version", 0)) for r in donor_events.values()) or
        any(r.get(name, 0) for r in donor_candidates.values() for name in (
            "native_timing_capture_version", "npb_capture_capture_version", "auau_npb_capture_capture_version")))
    definitions = {"timing": timing, "mbd": mbd, "npb": npb}
    if calo is not None:
        c._typed(calo, c.CaloDefinition, "calorimeter definition")
        definitions["calo"] = calo
    if captured and any(donor_config != definition.config_sha256 for definition in definitions.values()):
        raise c.ContractError("native capture config metadata is absent or contradicts definitions")
    if not donor_events.keys() <= events.keys() or not donor_candidates.keys() <= candidates.keys():
        raise c.ContractError("foreign native donor event/cluster population")
    semantics = dict(retained_semantics or {})
    if semantics.keys() - set(c.EVENT_FIELDS+c.CALO_FIELDS+c.CANDIDATE_FIELDS):
        raise c.ContractError("unknown retained semantic field")
    def existing(row, fields):
        return tuple(c.ExistingField(name, row[name] if math.isfinite(row[name]) else None,
            semantics.get(name)) for name in fields if name in row)
    def anchor(row):
        return c.CandidateAnchor(row["cluster_et"], row["eta"], row["phi"])
    def flags(row, other_names=None):
        names = tuple(name for name in OPTIONAL_FLAGS if name in row) if other_names is None else other_names
        _require(row, (*FLAGS, *names))
        return c.SelectionFlags(*(row[name] for name in FLAGS),
            tuple((name, row[name]) for name in names))
    retained_events, retained_candidates, event_donors, candidate_donors = [], [], [], []
    for physical, row in events.items():
        key = c.EventKey(binding, physical)
        fields = c.EVENT_FIELDS + (c.CALO_FIELDS if calo is not None else ())
        retained_events.append(c.RetainedEvent(key, existing(row, fields)))
        if physical in donor_events:
            captured_calo = (decode_calo(donor_events[physical], calo) if calo is not None
                             else c.Missing("calorimeter donor not requested"))
            event_donors.append(c.EventDonor(key, donor_hash, decode_mbd(donor_events[physical], mbd), captured_calo))
    for (physical, cluster), row in candidates.items():
        key = c.CandidateKey(c.EventKey(binding, physical), cluster)
        retained = c.RetainedCandidate(key, anchor(row), flags(row), existing(row, c.CANDIDATE_FIELDS))
        retained_candidates.append(retained)
        d = donor_candidates.get((physical, cluster))
        if d is not None:
            donor_flags = flags(d, tuple(name for name, _ in retained.selection.other_flags))
            if anchor(d) != retained.anchor or donor_flags != retained.selection:
                raise c.ContractError("native donor candidate anchor or selection flags differ")
            candidate_donors.append(c.CandidateDonor(key, anchor(d), donor_flags, donor_hash,
                decode_timing(d, timing), decode_npb(d, retained, npb,
                    model_file=model_file, evidence=donor_hash)))
    report = c.reconcile(binding, tuple(retained_events), tuple(retained_candidates), tuple(event_donors),
                         tuple(candidate_donors), timing=timing, mbd=mbd, npb=npb, calo=calo,
                         capture_source=donor_binding if donor_migration_policy else None)
    if file_sha256(source) != source_hash or file_sha256(donor) != donor_hash:
        raise c.ContractError("native augmentation input changed during decoding")
    receipt = dict(schema=VERSION, source_path=str(source), source_sha256=source_hash,
        donor_path=str(donor), donor_sha256=donor_hash, events=len(events), candidates=len(candidates),
        status="LOCAL_CAPTURE_DECODE_NOT_CERTIFIED", donor_coverage_complete=report.complete,
        scientific_acceptance=False, binding_authority="CALLER_DECLARED_RUNTIME_DEPENDENCIES_REQUIRE_ACCEPTANCE",
        npb_model_file=model_file, definitions={name: value.semantic_sha256 for name, value in definitions.items()})
    receipt.update(donor_migration_policy=donor_migration_policy,
        donor_binding=json.loads(json.dumps(asdict(donor_binding))) if donor_binding else None,
        donor_event_superset=donor_event_superset, excluded_donor_events=excluded_donor_events,
        excluded_donor_candidates=excluded_donor_candidates)
    return DecodedAugmentation(report, receipt)
