"""Materialize validated augmentation inside a canonical ROOT part.

No row-order joins, sentinel-to-valid inference, native mutation, model execution
or selection changes. Scalar values/states and all 25 named NPB inputs are ROOT
branches; repeated definitions/evidence are dictionary entries in one manifest.
The native tables remain recoverable; these tables are keyed within the part.
"""
from dataclasses import asdict
import json
import math
from pathlib import Path

import awkward as ak
import numpy as np
import uproot

import augmentation_contract as contract
from qa_augmentation import make_columns as recover_charge_columns
from canonical_repairs_io import (RepairTables, _column_digest, _columns, _identity,
    _put_id, file_sha256, write_scalar_tables, verify_scalar_tables)

VERSION = "CanonicalAugmentationTablesV3"
MANIFEST = "augmentation_manifest"
BASE = "ReplayFoundationV1/"
STATES = {"MISSING": 0, "VALID": 1, "INVALID": 2, "NOT_EVALUATED": 3, "NOT_APPLICABLE": 4}
KEY_SCHEMA = {f"{key}_{part}": "uint64" for key in ("part_id", "source_occurrence_id", "event_id") for part in ("hi", "lo")}
EVENT_SCHEMA = dict(KEY_SCHEMA, physical_event_sequence="int64", donor_evidence_index="int32")
EVENT_SCHEMA.update(mbd_native_is_valid="int32", mbd_npmt_south="int32", mbd_npmt_north="int32")
EVENT_SCHEMA.update(mbd_total_charge="float64", mbd_total_charge_state="int32", mbd_total_charge_origin="int32")
EVENT_SCHEMA.update(event_calo_capture_version="int32", event_calo_available_mask="uint32",
                    event_calo_valid_mask="uint32", event_calo_require_isgood="int32")
for _name in contract.EVENT_FIELDS + contract.CALO_FIELDS:
    EVENT_SCHEMA[_name], EVENT_SCHEMA[_name + "_state"] = "float64", "int32"
    EVENT_SCHEMA[_name + "_reason_index"] = "int32"
CANDIDATE_SCHEMA = dict(KEY_SCHEMA, candidate_id_hi="uint64", candidate_id_lo="uint64",
    physical_event_sequence="int64", native_cluster_key="uint32", donor_evidence_index="int32",
    cluster_time_raw="float64", cluster_time_raw_state="int32", cluster_time_ns="float64",
    cluster_time_raw_reason_index="int32", cluster_time_denominator_energy="float64",
    cluster_time_denominator_energy_state="int32", cluster_time_denominator_energy_reason_index="int32",
    cluster_time_contributing_towers="int32",
    npb_score="float64", npb_score_state="int32", npb_evaluation_state="int32",
    npb_score_reason_index="int32", npb_applicability_state="int32")
CANDIDATE_SCHEMA.update(npb_configured_cut="float64", npb_pass="int32", npb_pass_state="int32",
                        npb_pass_reason_index="int32")
CANDIDATE_SCHEMA.update(npb_scoring_binding_policy="int32", npb_scoring_node_index="int32",
    npb_scoring_native_cluster_key="uint32", npb_scoring_cluster_pt="float64",
    npb_scoring_eta="float64", npb_scoring_phi="float64", npb_scoring_input_et_bits="uint32",
    npb_scoring_input_eta_bits="uint32", npb_scoring_evidence_index="int32",
    npb_scoring_encounter_ordinal="int64")
for _name in contract.NPB_FEATURES:
    CANDIDATE_SCHEMA["npb_input_" + _name] = "float64"
    CANDIDATE_SCHEMA["npb_input_" + _name + "_state"] = "int32"
    CANDIDATE_SCHEMA["npb_input_" + _name + "_reason_index"] = "int32"
    CANDIDATE_SCHEMA["npb_input_" + _name + "_origin_index"] = "int32"
    CANDIDATE_SCHEMA["npb_input_" + _name + "_evidence_index"] = "int32"
SCHEMAS = {"eventAugmentations": EVENT_SCHEMA, "photonAugmentations": CANDIDATE_SCHEMA}
EVALUATION = {"MISSING": 0, "EVALUATED": 1, "NOT_EVALUATED": 2, "INPUTS_UNAVAILABLE": 3, "NOT_APPLICABLE": 4,
              "EVALUATED_INVALID": 5}
APPLICABILITY = {"UNVERIFIED": 0, "APPLICABLE": 1, "NOT_APPLICABLE": 2}
SCORING_BINDING = {"UNAVAILABLE": 0, "SAME_OBJECT_V1": 1, "NATIVE_FIRST_KINEMATIC_MATCH_V1": 2}


def _read(tree, max_rows, required, optional=()):
    if tree.num_entries > max_rows:
        raise contract.ContractError("one-part augmentation row bound exceeded")
    if set(required) - set(tree.keys()):
        raise contract.ContractError("required native identity/anchor fields absent")
    fields = list(dict.fromkeys([*required, *(n for n in optional if n in tree)]))
    return ak.to_list(tree.arrays(fields, library="ak"))


def _check_existing(native, retained, fields):
    reported = {field.name: field for field in retained.existing}
    for name in fields:
        if name not in native:
            if name in reported:
                raise contract.ContractError(f"reported retained branch is absent: {name}")
            continue
        value = float(native[name])
        if math.isfinite(value) and (name not in reported or reported[name].value != value):
            raise contract.ContractError(f"finite native value omitted or altered in report: {name}")
        if not math.isfinite(value) and name in reported and reported[name].value is not None:
            raise contract.ContractError(f"invented finite native value in report: {name}")


def _charge_rows(tree):
    """Reuse the independently tested finite-positive calibrated-PMT recovery."""
    keys = ("event_id_hi", "event_id_lo")
    fields = [*keys, *(n for n in ("mbd_total_charge", "mbd_pmt_available", "mbd_pmt_id",
        "mbd_pmt_charge", "mbd_pmt_charge_valid") if n in tree)]
    result = {}
    for chunk in tree.iterate(fields, step_size=4096, library="ak"):
        if "mbd_total_charge" not in chunk.fields:
            chunk = ak.with_field(chunk, np.full(len(chunk), np.nan), "mbd_total_charge")
        columns, _ = recover_charge_columns(chunk, keys)
        for i in range(len(chunk)):
            key = (int(columns[keys[0]][i]), int(columns[keys[1]][i]))
            if key in result:
                raise contract.ContractError("duplicate native event in charge recovery")
            result[key] = (float(columns["mbd_total_charge"][i]), int(columns["mbd_total_charge_state"][i]))
    return result


def prepare_augmentation_tables(source, report, *, timing, mbd, npb, max_rows=2_000_000,
                                capture_receipt=None, calo=None):
    """Revalidate a report and bind its exact retained population to native ROOT.

    Input DST/calibration/model hashes are retained claims, not independently
    verified remote bytes. Publication/physics admission must supply that proof;
    this adapter certifies local identity/value readback only.
    """
    if type(max_rows) is not int or max_rows <= 0:
        raise contract.ContractError("positive row bound required")
    source = Path(source).resolve()
    digest = file_sha256(source)
    if digest != report.source.retained_file_sha256:
        raise contract.ContractError("augmentation retained-file bytes differ")
    checked = contract.reconcile(report.source,
        tuple(r.retained for r in report.events), tuple(r.retained for r in report.candidates),
        tuple(r.donor for r in report.events if r.donor is not None),
        tuple(r.donor for r in report.candidates if r.donor is not None), timing=timing, mbd=mbd, npb=npb, calo=calo,
        capture_source=report.capture_source)
    if checked != report:
        raise contract.ContractError("overlay report was changed after reconciliation")
    event_ids = ["event_id_hi", "event_id_lo"]
    source_ids = ["source_occurrence_id_hi", "source_occurrence_id_lo"]
    candidate_ids_fields = ["candidate_id_hi", "candidate_id_lo"]
    other_flags = {name for r in report.candidates for name, _ in r.retained.selection.other_flags}
    event_fields = contract.EVENT_FIELDS + (contract.CALO_FIELDS if calo is not None else ())
    with uproot.open(source) as root:
        events = _read(root[BASE + "RJEventV1"], max_rows,
            event_ids + source_ids + ["run", "physical_event_sequence", "physical_event_sequence_valid"], event_fields)
        candidates = _read(root[BASE + "RJPhotonCandidateV1"], max_rows,
            event_ids + candidate_ids_fields + ["native_cluster_key", "cluster_et", "eta", "phi",
                "reference_preselection_state", "active_preselection_state", "preselection_bitmask"],
            (*contract.CANDIDATE_FIELDS, *sorted(other_flags)))
        sources = _read(root[BASE + "RJSourceOccurrenceV1"], max_rows,
            source_ids + ["run", "segment", "source_manifest_sha256"])
        charges = _charge_rows(root[BASE + "RJEventV1"])
    by_source = {_identity(s, "source_occurrence_id"): s for s in sources}
    if len(by_source) != len(sources):
        raise contract.ContractError("duplicate native source identity")
    by_event, by_physical = {}, {}
    for row in events:
        event, occurrence = _identity(row, "event_id"), _identity(row, "source_occurrence_id")
        physical = row["physical_event_sequence"]
        if row["physical_event_sequence_valid"] != 1 or type(physical) is not int or physical < 0:
            raise contract.ContractError("native physical-event identity unavailable")
        if event in by_event or physical in by_physical or occurrence not in by_source:
            raise contract.ContractError("ambiguous native physical/event/source identity")
        native_source = by_source[occurrence]
        if (row["run"] != report.source.run or native_source["run"] != report.source.run
                or native_source["segment"] != report.source.segment
                or native_source["source_manifest_sha256"] != report.source.source_manifest_sha256):
            raise contract.ContractError("native run/segment/manifest binding differs")
        by_event[event], by_physical[physical] = row, row
    if set(by_physical) != {r.retained.key.physical_event for r in report.events}:
        raise contract.ContractError("report event population differs from native inventory")
    by_candidate = {}
    candidate_ids = set()
    for row in candidates:
        event = by_event.get(_identity(row, "event_id"))
        if event is None:
            raise contract.ContractError("native candidate has no event")
        key = (event["physical_event_sequence"], row["native_cluster_key"])
        identity = _identity(row, "candidate_id")
        if key in by_candidate or identity in candidate_ids:
            raise contract.ContractError("ambiguous native cluster key or candidate identity")
        by_candidate[key] = row
        candidate_ids.add(identity)
    if set(by_candidate) != {(r.retained.key.event.physical_event, r.retained.key.native_cluster_key) for r in report.candidates}:
        raise contract.ContractError("report candidate population differs from native inventory")
    evidence = sorted({r.donor.capture_sha256 for r in report.events + report.candidates if r.donor})
    if report.capture_source is not None and (capture_receipt is None or
            capture_receipt.get('donor_binding') != json.loads(json.dumps(asdict(report.capture_source)))):
        raise contract.ContractError('migrated capture requires its exact native donor binding receipt')
    if capture_receipt is not None:
        expected_definitions = {"timing": timing.semantic_sha256, "mbd": mbd.semantic_sha256,
                                "npb": npb.semantic_sha256}
        if calo is not None:
            expected_definitions["calo"] = calo.semantic_sha256
        if (capture_receipt.get("schema") != "NativeAugmentationReaderV1" or
                capture_receipt.get("source_sha256") != digest or
                capture_receipt.get("source_path") != str(source) or
                capture_receipt.get("definitions") != expected_definitions or
                capture_receipt.get("events") != len(report.events) or
                capture_receipt.get("candidates") != len(report.candidates) or
                capture_receipt.get("status") != "LOCAL_CAPTURE_DECODE_NOT_CERTIFIED" or
                capture_receipt.get("scientific_acceptance") is not False or
                evidence != [capture_receipt.get("donor_sha256")]):
            raise contract.ContractError("native capture receipt contradicts augmentation report")
        # The donor is a live local input until final publication. Do not bind
        # a recorded hash to a path whose bytes have since changed.
        if file_sha256(capture_receipt["donor_path"]) != capture_receipt["donor_sha256"]:
            raise contract.ContractError("native donor bytes changed after decoding")
    part_id = (int(digest[:16], 16), int(digest[16:32], 16))

    def key_fields(physical, donor):
        event = by_physical[physical]
        out = dict(physical_event_sequence=physical,
                   donor_evidence_index=evidence.index(donor.capture_sha256) if donor else -1)
        for name, identity in (("part_id", part_id), ("event_id", _identity(event, "event_id")),
                               ("source_occurrence_id", _identity(event, "source_occurrence_id"))):
            _put_id(out, name, identity)
        return out

    reason_dictionary, origin_dictionary, feature_evidence = [], [], []
    scoring_nodes, scoring_evidence = [], []

    def dictionary_index(values, value):
        if value is None or value == "":
            return -1
        if value not in values:
            values.append(value)
        return values.index(value)

    def measurement(out, name, value):
        out[name] = math.nan if value.value is None else value.value
        out[name + "_state"] = STATES[value.state]
        out[name + "_reason_index"] = dictionary_index(reason_dictionary, value.reason)

    output_events, output_candidates, reasons = [], [], {}
    for row in report.events:
        physical = row.retained.key.physical_event
        _check_existing(by_physical[physical], row.retained, event_fields)
        out = key_fields(physical, row.donor)
        charge, origin = charges[_identity(by_physical[physical], "event_id")]
        out.update(mbd_total_charge=charge, mbd_total_charge_state=1 if origin else 0,
                   mbd_total_charge_origin=origin)
        for field in row.fields:
            measurement(out, field.name, field.measurement)
        capture = row.donor.mbd if row.donor and isinstance(row.donor.mbd, contract.MbdCapture) else None
        out.update(mbd_native_is_valid=int(capture.native_is_valid) if capture else -1,
                   mbd_npmt_south=capture.south_npmt if capture else -1,
                   mbd_npmt_north=capture.north_npmt if capture else -1)
        captured_calo = row.donor.calo if calo is not None and row.donor and isinstance(row.donor.calo, contract.CaloCapture) else None
        out.update(event_calo_capture_version=1 if captured_calo else 0,
                   event_calo_available_mask=captured_calo.available_mask if captured_calo else 0,
                   event_calo_valid_mask=captured_calo.valid_mask if captured_calo else 0,
                   event_calo_require_isgood=int(calo.require_isgood) if calo is not None else -1)
        if calo is None:
            for name in contract.CALO_FIELDS:
                measurement(out, name, contract.Measurement("MISSING", reason="calorimeter augmentation not requested"))
        output_events.append(out)
    for row in report.candidates:
        physical, cluster = row.retained.key.event.physical_event, row.retained.key.native_cluster_key
        native = by_candidate[(physical, cluster)]
        retained = row.retained
        if contract.CandidateAnchor(native["cluster_et"], native["eta"], native["phi"]) != retained.anchor:
            raise contract.ContractError("native candidate kinematic anchor differs")
        flags = retained.selection
        for name, value in (("reference_preselection_state", flags.reference_preselection_state),
                            ("active_preselection_state", flags.active_preselection_state),
                            ("preselection_bitmask", flags.preselection_bitmask), *flags.other_flags):
            if name not in native or native[name] != value:
                raise contract.ContractError("native immutable selection flag differs")
        _check_existing(native, retained, contract.CANDIDATE_FIELDS)
        out = key_fields(physical, row.donor)
        out["native_cluster_key"] = cluster
        _put_id(out, "candidate_id", _identity(native, "candidate_id"))
        for field in row.fields:
            name = "npb_score" if field.name == npb.score_field else field.name
            measurement(out, name, field.measurement)
        out["cluster_time_ns"] = out["cluster_time_raw"] * timing.sample_to_ns if out["cluster_time_raw_state"] == STATES["VALID"] else math.nan
        timing_capture = row.donor.timing if row.donor and isinstance(row.donor.timing, contract.TimingCapture) else None
        measurement(out, "cluster_time_denominator_energy", timing_capture.denominator_energy if timing_capture
                    else contract.Measurement("MISSING", reason="no timing denominator donor"))
        out["cluster_time_contributing_towers"] = (timing_capture.contributing_towers
            if timing_capture and timing_capture.contributing_towers is not None else -1)
        capture = row.donor.npb if row.donor and isinstance(row.donor.npb, contract.NPBCapture) else None
        out["npb_evaluation_state"] = EVALUATION[capture.evaluation_state if capture else "MISSING"]
        out["npb_applicability_state"] = APPLICABILITY[capture.applicability_state if capture else "UNVERIFIED"]
        passed = capture.pass_measurement if capture else contract.Measurement("MISSING", reason="no NPB donor")
        out["npb_configured_cut"] = npb.configured_cut if npb.configured_cut is not None else math.nan
        out["npb_pass"] = int(passed.value) if passed.state == "VALID" else -1
        out["npb_pass_state"] = STATES[passed.state]
        out["npb_pass_reason_index"] = dictionary_index(reason_dictionary, passed.reason)
        witness = capture.scoring_object if capture else None
        # Source and physical-event keys are shared with the retained photon;
        # reconcile() has already rejected any cross-source/event witness.
        out.update(
            npb_scoring_binding_policy=SCORING_BINDING[witness.binding_policy if witness else "UNAVAILABLE"],
            npb_scoring_node_index=dictionary_index(scoring_nodes, witness.cluster_node if witness else None),
            npb_scoring_native_cluster_key=witness.key.native_cluster_key if witness else 0,
            npb_scoring_cluster_pt=witness.anchor.cluster_et if witness else math.nan,
            npb_scoring_eta=witness.anchor.eta if witness else math.nan,
            npb_scoring_phi=witness.anchor.phi if witness else math.nan,
            npb_scoring_input_et_bits=int(witness.input_et_bits, 16) if witness else 0,
            npb_scoring_input_eta_bits=int(witness.input_eta_bits, 16) if witness else 0,
            npb_scoring_evidence_index=dictionary_index(scoring_evidence, witness.evidence_sha256 if witness else None),
            npb_scoring_encounter_ordinal=(witness.matched_encounter_ordinal
                if witness and witness.matched_encounter_ordinal is not None else -1))
        for i, name in enumerate(contract.NPB_FEATURES):
            value = capture.inputs[i].measurement if capture else contract.Measurement("MISSING", reason="no feature donor")
            measurement(out, "npb_input_" + name, value)
            item = capture.inputs[i] if capture else None
            out["npb_input_" + name + "_origin_index"] = dictionary_index(origin_dictionary, item.origin if item else None)
            out["npb_input_" + name + "_evidence_index"] = dictionary_index(feature_evidence, item.evidence_sha256 if item else None)
        output_candidates.append(out)
    # Unresolved reasons are dictionary/counter metadata, not a repeated source
    # object/config blob per photon. Explicit states remain per row.
    for _, name, state, reason in report.unresolved:
        token = f"{name}|{state}|{reason}"
        reasons[token] = reasons.get(token, 0) + 1
    tables = {"eventAugmentations": _columns(output_events, EVENT_SCHEMA),
              "photonAugmentations": _columns(output_candidates, CANDIDATE_SCHEMA)}
    if file_sha256(source) != digest:
        raise contract.ContractError("source changed while preparing augmentation")
    charge_counts = {str(origin): sum(value[1] == origin for value in charges.values()) for origin in (0, 1, 2)}
    missing_pass = sum(r["npb_pass_state"] not in (STATES["VALID"], STATES["NOT_APPLICABLE"]) for r in output_candidates)
    receipt = dict(schema=VERSION, source_path=str(source), source_sha256=digest,
        status="LOCAL_AUGMENTATION_NOT_CERTIFIED",
        coverage_scope=["timing", "mbd_event_timing", "npb", "mbd_total_charge"] + (["calorimeter_sums"] if calo is not None else []),
        coverage_complete=report.complete and charge_counts["0"] == 0 and missing_pass == 0,
        donor_coverage_complete=report.complete, mbd_total_charge_counts_by_origin=charge_counts,
        npb_pass_unresolved_candidates=missing_pass, npb_pass_predicate="VALID_EVALUATED_SCORE_STRICTLY_GREATER_THAN_BOUND_CUT",
        scientific_acceptance=False, input_hash_authority="CALLER_CLAIMS_REQUIRE_BYTE_VERIFICATION",
        native_capture_receipt=capture_receipt,
        capture_source_binding=asdict(report.capture_source) if report.capture_source else None,
        source_binding=asdict(report.source), definitions=dict(timing=asdict(timing), mbd=asdict(mbd), npb=asdict(npb),
                                                            calo=asdict(calo) if calo is not None else None),
        states=STATES, evaluation_states=EVALUATION, applicability_states=APPLICABILITY,
        mbd_total_charge_origins={"0": "UNAVAILABLE", "1": "NATIVE", "2": "RECOVERED_EXACT_FINITE_POSITIVE_PMT_SUM"},
        donor_evidence_sha256=evidence,
        feature_origin_dictionary=origin_dictionary, feature_evidence_sha256=feature_evidence,
        scoring_binding_policies=SCORING_BINDING, scoring_node_dictionary=scoring_nodes,
        scoring_evidence_sha256=scoring_evidence,
        scoring_witness_authority="CAPTURE_CLAIMS_REQUIRE_NATIVE_IDENTITY_AND_FIRST_MATCH_VERIFICATION",
        reason_dictionary=reason_dictionary, unresolved_reason_counts=reasons,
        tables={name: dict(rows=len(next(iter(columns.values()))), schema=SCHEMAS[name], sha256=_column_digest(columns)) for name, columns in tables.items()})
    # Normalize tuple-bearing dataclasses to the actual JSON manifest value.
    return RepairTables(tables, json.loads(json.dumps(receipt, sort_keys=True, allow_nan=False)))


def write_augmentation_tables(directory, prepared, **kwargs):
    return write_scalar_tables(directory, prepared, manifest_name=MANIFEST, **kwargs)


def verify_augmentation_tables(path, prepared):
    return verify_scalar_tables(path, prepared, manifest_name=MANIFEST)
