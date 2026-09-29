"""Source-bound repairs for retained photons and reconstructed jets.

The original mixed photon/jet link table is preserved with its per-row policy
provenance. Its photon rows are not trusted as the repaired association: both
legacy and new native links are rebuilt from dominant-primary witnesses. The
two new tables have a fixed scalar schema readable by ROOT and uproot. Unknown
links have class -1 (also supported by native photon policy 1). Optional evaluated
JES inputs replace only jet corrected_pt and dependent retained-pair xjgamma.

This adapter is intentionally bounded to one native part. It does not infer
capture completeness or promote a response denominator from a storage test.
"""
from __future__ import annotations

from array import array
from collections import Counter
from dataclasses import dataclass
import hashlib
import json
import math
import operator
from pathlib import Path
import time

import awkward as ak
import numpy as np
import uproot

from canonical_storage import file_sha256
from physics_repairs import rebuild_photon_links, RepairError, JetPTRepair, repair_retained_jet_pairs

VERSION = "CanonicalPhotonRepairTablesV1"
BASE = "ReplayFoundationV1/"
ASSOCIATION_STATES = {"NOT_APPLICABLE": -1, "UNKNOWN": 0, "NO_PRIMARY": 1,
                      "MATCH": 2, "NON_PHOTON_PRIMARY": 3}
LINK_CLASSES = {"UNKNOWN": -1, "MATCH": 0, "RECO_FAKE": 1, "TRUTH_MISS": 2}


def native_jet_donor_repairs(source, donor, *, source_binding, donor_binding,
                            capture_binding, max_rows=2_000_000):
    """Reuse calibrated jets from the SAME augmentation capture, never rerun JES.

    The caller must validate the donor's native producer request/terminal receipt.
    Source bindings pin the physical DST files and common non-JES calibration
    basis; producer identities may differ for the additive capture migration.
    Only corrected_pt is taken from the donor. All other old values/identities
    stay unchanged, and the existing writer recomputes retained pair xjgamma.
    A changed raw jet, missing jet, or extra jet in a retained event is an error.
    This is not a claim that old retained populations constitute a full release.
    """
    import augmentation_contract as c
    if (not isinstance(source_binding, c.SourceBinding)
            or not isinstance(donor_binding, c.SourceBinding)):
        raise RepairError('explicit typed source/donor bindings required')
    for name in ('run', 'segment', 'sample', 'inputs', 'calibration_sha256'):
        if getattr(source_binding, name) != getattr(donor_binding, name):
            raise RepairError('jet donor input/calibration binding differs: '+name)
    source, donor = Path(source).resolve(), Path(donor).resolve()
    if source == donor:
        raise RepairError('retained repair requires a distinct augmentation donor')
    source_sha, donor_sha = file_sha256(source), file_sha256(donor)
    if (source_sha != source_binding.retained_file_sha256
            or donor_sha != donor_binding.retained_file_sha256):
        raise RepairError('jet donor source bytes differ from bindings')
    if (not isinstance(capture_binding, dict)
            or capture_binding.get('native_sha256') != donor_sha
            or capture_binding.get('producer_library_sha256') != donor_binding.reconstruction_sha256
            or capture_binding.get('physics_acceptance') is not False):
        raise RepairError('validated donor capture binding required')
    if type(max_rows) is not int or max_rows <= 0:
        raise RepairError('positive finite jet donor row bound required')

    # Also bind the stored occurrences, not just the externally supplied labels.
    occurrence_fields = ('lane', 'dataset', 'sample', 'period', 'run', 'segment',
                         'si_di_role', 'input_uri_hash')
    inventories = []
    for path, binding in ((source, source_binding), (donor, donor_binding)):
        with uproot.open(path, array_cache=None) as root:
            occurrences = _rows(root[BASE+'RJSourceOccurrenceV1'],
                                list(occurrence_fields)+['source_manifest_sha256'], [], max_rows)
            if len(occurrences) != 1:
                raise RepairError('one exact source occurrence per donor part required')
            occurrence = occurrences[0]
            if binding.input_tuple_sha256 is not None and occurrence['input_uri_hash'] != binding.input_tuple_sha256:
                raise RepairError('jet donor frozen input tuple differs')
            if (occurrence['run'], occurrence['segment'], occurrence.pop('source_manifest_sha256')) != (
                    binding.run, binding.segment, binding.source_manifest_sha256):
                raise RepairError('stored jet source binding differs')
            events = _rows(root[BASE+'RJEventV1'], _ids('event_id') +
                ['physical_event_sequence', 'physical_event_sequence_valid'], [], max_rows)
            index = {}
            physical = set()
            for row in events:
                eid, sequence = _identity(row, 'event_id'), row['physical_event_sequence']
                if (row['physical_event_sequence_valid'] != 1 or type(sequence) is not int
                        or sequence < 0 or eid in index or sequence in physical):
                    raise RepairError('unique valid physical event identities required for JES donor')
                index[eid] = sequence
                physical.add(sequence)
            jets = _rows(root[BASE+'RJJetV1'], list(root[BASE+'RJJetV1'].keys()), [], max_rows)
            inventory = {}
            for row in jets:
                eid = _identity(row, 'event_id')
                if eid not in index:
                    raise RepairError('jet has no physical event')
                key = (index[eid], *(row[k] for k in ('algorithm', 'radius', 'input_identity',
                       'subtraction_identity', 'native_raw_jet_key')))
                if key in inventory:
                    raise RepairError('duplicate native raw-jet identity')
                inventory[key] = row
            inventories.append((occurrences[0], physical, inventory))
    old, new = inventories
    if old[0] != new[0] or not old[1] <= new[1]:
        raise RepairError('jet donor source occurrence/event coverage differs')
    selected = {k: v for k, v in new[2].items() if k[0] in old[1]}
    if selected.keys() != old[2].keys():
        raise RepairError('retained-event native jet population differs')
    provenance = dict(source_scope='sha256:'+source_sha,
        evaluation_semantics='EXACT_SAVED_NATIVE_JETCALIB_OUTPUT',
        donor_sha256=donor_sha,
        producer_library_sha256=capture_binding['producer_library_sha256'],
        capture_request_sha256=capture_binding['request_file_sha256'],
        capture_receipt_sha256=capture_binding['receipt_file_sha256'],
        calibration_basis_sha256=source_binding.calibration_sha256,
        correction_application='REPLACE_FROM_NATIVE_DONOR_NEVER_MULTIPLY_CORRECTED_PT')
    repairs = []
    for key, row in old[2].items():
        other = selected[key]
        for name in ('raw_pt', 'deterministic_order', 'native_jet_key'):
            if row[name] != other[name]:
                raise RepairError('jet donor changed raw jet: '+name)
        for name in ('eta', 'phi'):
            difference = row[name]-other[name]
            if name == 'phi':
                difference = math.remainder(difference, 2*math.pi)
            if not math.isfinite(difference) or abs(difference) > 1e-5:
                raise RepairError('jet donor direction differs beyond native float roundoff')
        pt = other['corrected_pt']
        if not math.isfinite(pt) or pt < 0:
            raise RepairError('invalid native calibrated jet pT')
        repairs.append(JetPTRepair(dict(row, corrected_pt=pt), 'CALIBRATED_PT_ONLY',
            'PRE_JES', 'POST_JES', (), 'RETAINED_ROWS_ONLY', dict(provenance)))
    if file_sha256(source) != source_sha or file_sha256(donor) != donor_sha:
        raise RepairError('jet donor/source changed during join')
    return repairs
RECO_CENSUS_FIELDS = tuple("reco_photon_" + name for name in (
    "capture_version", "capture_domain", "capture_state", "inspected",
    "eligible", "written", "unclassifiable"))
TRUTH_CENSUS_FIELDS = tuple("truth_photon_" + name for name in (
    "capture_state", "native_embedded_primary_pid22_count", "serialized_raw_count",
    "rejected_count", "rejected_invalid_kinematics_count",
    "rejected_nonfinite_isolation_count", "duplicate_track_count",
    "isolation_input_incomplete_count", "analysis_signal_count")) + ("truth_denominator_complete",)


def native_truth_complete_events(events, truths):
    """Cross-check saved native photon census against keyed, finite truth rows.

    Selection-neutral inventory evidence, not physics acceptance. Legacy input
    stays unknown. Require zero rejected objects for this stricter complete
    retained-inventory claim, including invalid-kinematics rejections which the
    upstream denominator flag alone permits.
    """
    grouped = {}
    for row in truths:
        grouped.setdefault(_identity(row, "event_id"), []).append(row)
    observed, complete = set(), set()
    for event in events:
        eid = _identity(event, "event_id")
        if eid in observed:
            raise RepairError("duplicate event in native truth census")
        observed.add(eid)
        present = set(TRUTH_CENSUS_FIELDS) & event.keys()
        if not present:
            continue
        if present != set(TRUTH_CENSUS_FIELDS):
            raise RepairError("partial native truth census")
        raw = [event[k] for k in TRUTH_CENSUS_FIELDS]
        try:
            if any(isinstance(v, bool) for v in raw): raise TypeError("boolean")
            state, native, written, rejected, badkin, badiso, duplicate, incomplete, signal, flag = map(operator.index, raw)
        except TypeError as exc:
            raise RepairError("noninteger native truth census") from exc
        if (state not in (-2, -1, 0, 1) or flag not in (0, 1)
                or min(native, written, rejected, badkin, badiso, duplicate, incomplete, signal) < 0
                or max(native, written, rejected, badkin, badiso, duplicate, incomplete, signal) > 2**31-1):
            raise RepairError("invalid native truth census state/counts")
        if state != 1:
            if any(raw[1:]):
                raise RepairError("uncaptured native truth census claims counts/completeness")
            continue
        if native != written + rejected or rejected != badkin + badiso + duplicate or signal > written:
            raise RepairError("native truth census accounting does not close")
        expected_flag = int(incomplete == 0 and duplicate == 0 and badiso == 0)
        if flag != expected_flag:
            raise RepairError("native truth denominator flag contradicts census")
        rows = grouped.get(eid, [])
        if len(rows) != written:
            raise RepairError("native truth census written count differs from retained rows")
        keys, signals = set(), 0
        for row in rows:
            try:
                raw_key = (row['embedding_id'], row['native_track_id'])
                sig = row['analysis_signal_r03']
                if any(isinstance(v, bool) for v in raw_key + (sig,)): raise TypeError("boolean")
                key = tuple(map(operator.index, raw_key)); sig = operator.index(sig)
                pt, eta, phi = (float(row[k]) for k in ('pt', 'eta', 'phi'))
            except (KeyError, TypeError, ValueError) as exc:
                raise RepairError("missing native truth row witnesses") from exc
            if key in keys or sig not in (0, 1) or pt <= 0 or not all(map(math.isfinite, (pt, eta, phi))):
                raise RepairError("native truth key/kinematics/signal contradiction")
            keys.add(key); signals += sig
        if signals != signal:
            raise RepairError("native truth signal count differs from retained rows")
        if flag == 1 and rejected == 0:
            complete.add(eid)
    if grouped.keys() - observed:
        raise RepairError("truth absent from native census event inventory")
    return complete


def native_reco_complete_events(events, candidates, *, domain):
    """Verify persisted domain-1 census against every retained candidate.

    This establishes capture only within finite 5<=pt<40, |eta|<0.7 in the
    reconstructed photon container, not raw clusters or a wider population.
    Missing/legacy evidence stays unknown; malformed positive evidence fails.
    It does not certify truth completeness or physics selections.
    """
    if type(domain) is not int or domain != 1:
        raise RepairError("only explicit native reco capture domain 1 is supported")
    grouped = {}
    for row in candidates:
        grouped.setdefault(_identity(row, "event_id"), []).append(row)
    complete, observed = set(), set()
    for event in events:
        eid = _identity(event, "event_id")
        if eid in observed:
            raise RepairError("duplicate event in native reco census")
        observed.add(eid)
        present = set(RECO_CENSUS_FIELDS) & event.keys()
        if not present:
            continue
        if present != set(RECO_CENSUS_FIELDS):
            raise RepairError("partial native reco census")
        values = [event[k] for k in RECO_CENSUS_FIELDS]
        if any(isinstance(v, bool) for v in values):
            raise RepairError("boolean native reco census")
        try:
            version, saved_domain, state, inspected, eligible, written, invalid = map(operator.index, values)
        except TypeError as exc:
            raise RepairError("noninteger native reco census") from exc
        if version == 0:
            if (saved_domain, state, inspected, eligible, written, invalid) != (0, 0, -1, -1, -1, -1):
                raise RepairError("contradictory unrecorded native reco census")
            continue
        if (version != 1 or saved_domain != domain or state not in (1, 2)
                or min(inspected, eligible, written, invalid) < 0
                or inspected > 2**31-1 or not written <= eligible <= inspected
                or invalid > inspected-eligible):
            raise RepairError("invalid native reco census state/counts/domain")
        rows = grouped.get(eid, [])
        keys = set()
        if len(rows) != written:
            raise RepairError("native reco census written count differs from retained rows")
        for row in rows:
            try:
                native_key = row["native_cluster_key"]
                if isinstance(native_key, bool): raise TypeError("boolean key")
                native_key = operator.index(native_key)
                pt, eta, phi = (float(row[k]) for k in ("cluster_et", "eta", "phi"))
            except (KeyError, TypeError, ValueError) as exc:
                raise RepairError("missing native reco census row witnesses") from exc
            if (not 0 <= native_key < 2**32 or native_key in keys
                    or not all(map(math.isfinite, (pt, eta, phi)))
                    or not 5 <= pt < 40 or not abs(eta) < .7):
                raise RepairError("native reco census key/domain contradiction")
            keys.add(native_key)
        if state == 1:
            if eligible != written or invalid != 0:
                raise RepairError("complete native reco census has missing/invalid objects")
            complete.add(eid)
    if grouped.keys() - observed:
        raise RepairError("candidate absent from native reco census event inventory")
    return complete
ID_SCHEMA = {f"{name}_{part}": "uint64" for name in (
    "part_id", "source_occurrence_id", "event_id", "candidate_id", "truth_photon_id"
) for part in ("hi", "lo")}
ASSOCIATION_SCHEMA = dict(ID_SCHEMA, association_state="int32",
    dominant_energy_contribution="float64", diagnostic_delta_r="float64")
LINK_SCHEMA = dict(ID_SCHEMA, link_id_hi="uint64", link_id_lo="uint64",
                   link_class="int32", association_state="int32")
SCHEMAS = {"photonAssociations": ASSOCIATION_SCHEMA, "photonResponseLinks": LINK_SCHEMA}


def _identity(row, prefix):
    try:
        raw = [row[f"{prefix}_{p}"] for p in ("hi", "lo")]
        if any(isinstance(x, bool) for x in raw):
            raise TypeError("boolean identity")
        values = tuple(operator.index(x) for x in raw)
    except (KeyError, TypeError) as exc:
        raise RepairError(f"invalid {prefix}") from exc
    if any(x < 0 or x >= 2**64 for x in values) or values == (0, 0):
        raise RepairError(f"invalid {prefix}")
    return values


def _put_id(row, name, value):
    row[name + "_hi"], row[name + "_lo"] = value or (0, 0)


def _hash_id(value):
    digest = hashlib.sha256(value.encode()).digest()
    return int.from_bytes(digest[:8], "big"), int.from_bytes(digest[8:16], "big")


def _rows(tree, required, optional, max_rows):
    if tree.num_entries > max_rows:
        raise RepairError(f"one-part repair row bound exceeded: {tree.object_path}")
    missing = set(required) - set(tree.keys())
    if missing:
        raise RepairError(f"required identity fields missing: {sorted(missing)}")
    return ak.to_list(tree.arrays(list(required) + [n for n in optional if n in tree], library="ak"))


def _native_link_provenance(root, max_rows):
    """Inventory policy markers only; a marker is not association/response proof."""
    path = BASE + "RJRecoTruthLinkV1"
    result = {"table_present": path in root, "policy_branch_present": False,
        "photon_policy_counts": {"LEGACY_UNSPECIFIED": 0, "ENERGY_PRIMARY_IDENTITY_V1": 0},
        "jet_rows": 0, "other_rows": 0,
        "qualification": "NOT_CERTIFIED",
        "use": "PRESERVED_ONLY_REBUILD_PHOTONS_FROM_DOMINANT_PRIMARY_WITNESSES"}
    if path not in root:
        return result
    tree = root[path]
    result["policy_branch_present"] = "photon_association_policy" in tree
    for row in _rows(tree, ["reco_type", "truth_type", "link_class"],
                     ["photon_association_policy"], max_rows):
        try:
            values = [row["reco_type"], row["truth_type"], row["link_class"],
                      row.get("photon_association_policy", 0)]
            if any(isinstance(value, bool) for value in values):
                raise TypeError("boolean marker")
            rt, tt, lc, policy = map(operator.index, values)
        except TypeError as exc:
            raise RepairError("native link types/class/policy must be integers") from exc
        if policy not in (0, 1) or rt not in (0, 1, 2) or tt not in (0, 1, 2):
            raise RepairError("unsupported native photon association policy/type")
        if policy == 1:
            permitted = {0: {(1, 1)}, 1: {(1, 0)}, 2: {(0, 1)},
                         -1: {(1, 0), (0, 1)}}
            if (rt, tt) not in permitted.get(lc, set()):
                raise RepairError("invalid energy-primary photon link class/endpoints")
        if rt == 1 or tt == 1:
            name = "ENERGY_PRIMARY_IDENTITY_V1" if policy else "LEGACY_UNSPECIFIED"
            result["photon_policy_counts"][name] += 1
        elif rt == 2 or tt == 2:
            result["jet_rows"] += 1
        else:
            result["other_rows"] += 1
    return result


def _ids(*prefixes):
    return [f"{prefix}_{part}" for prefix in prefixes for part in ("hi", "lo")]


@dataclass(frozen=True)
class PreparedJetRepair:
    replacements: dict
    receipt: dict


def prepare_jet_replacements(source, repairs, *, max_rows=2_000_000):
    """Bind evaluated JES repairs to every retained jet and derive pair xJ.

    The evaluator supplies JetPTRepair objects; this adapter never chooses a
    payload or infers pre-JES state from unity branches. All source jet fields
    except corrected_pt must be present and unchanged in each evaluated row.
    Output arrays follow native entry order, irrespective of evaluator order.
    """
    if type(max_rows) is not int or max_rows <= 0:
        raise RepairError('positive integer row bound required')
    if not isinstance(repairs, (list, tuple)) or len(repairs) > max_rows:
        raise RepairError('bounded explicit jet repair sequence required')
    source = Path(source).resolve()
    source_hash = file_sha256(source)
    scope = 'sha256:' + source_hash
    with uproot.open(source, array_cache=None) as root:
        tree = root[BASE + 'RJJetV1']
        jets = _rows(tree, tree.keys(), (), max_rows)
        photons = _rows(root[BASE + 'RJPhotonCandidateV1'],
                        _ids('event_id', 'candidate_id') + ['cluster_et'], (), max_rows)
        pairs = _rows(root[BASE + 'RJPhotonJetPairV1'],
                      _ids('event_id', 'pair_id', 'candidate_id', 'jet_id') +
                      ['jet_rank', 'xjgamma'], (), max_rows)
    indexed = {}
    provenance = set()
    for repair in repairs:
        if not isinstance(repair, JetPTRepair) or repair.status not in (
                'CALIBRATED_PT_ONLY', 'JES_INVALID_NATIVE_SIGN_PRESERVED'):
            raise RepairError('only classified reconstructed jet repairs are accepted')
        if repair.provenance.get('source_scope') != scope:
            raise RepairError('jet repair is not bound to exact source bytes')
        key = _identity(repair.jet, 'jet_id')
        if key in indexed:
            raise RepairError('duplicate repaired jet identity')
        indexed[key] = repair
        provenance.add(json.dumps(dict(repair.provenance,
            raw_pt_semantics=repair.raw_pt_semantics,
            corrected_pt_semantics=repair.corrected_pt_semantics), sort_keys=True, allow_nan=False))
    ordered, seen = [], set()
    for row in jets:
        key = _identity(row, 'jet_id')
        if key in seen or key not in indexed:
            raise RepairError('missing or duplicate native jet identity')
        seen.add(key)
        repair = indexed[key]
        if repair.status == 'JES_INVALID_NATIVE_SIGN_PRESERVED' and repair.jet.get('corrected_pt') != row['corrected_pt']:
            raise RepairError('invalid-JES jet must preserve original corrected_pt')
        for field, expected in row.items():
            if field == 'corrected_pt':
                continue
            actual = repair.jet.get(field)
            both_nan = (isinstance(expected, float) and isinstance(actual, float)
                        and math.isnan(expected) and math.isnan(actual))
            if not both_nan and actual != expected:
                raise RepairError(f'jet repair altered source field: {field}')
        ordered.append(repair)
    if set(indexed) != seen:
        raise RepairError('repair contains foreign jet identities')
    repaired_pairs = repair_retained_jet_pairs(pairs, photons, ordered, source_scope=scope)
    invalid_ids = [_identity(r.jet, 'jet_id') for r in ordered
                   if r.status == 'JES_INVALID_NATIVE_SIGN_PRESERVED']
    invalid_id_set = set(invalid_ids)
    invalid_pair_ids = [_identity(row, 'pair_id') for row in pairs
                        if _identity(row, 'jet_id') in invalid_id_set]
    replacements = {BASE + 'RJJetV1': {'corrected_pt': np.asarray(
        [r.jet['corrected_pt'] for r in ordered], dtype='float64')},
        BASE + 'RJPhotonJetPairV1': {'xjgamma': np.asarray(
        [r['xjgamma'] for r in repaired_pairs.pairs], dtype='float64')}}
    if file_sha256(source) != source_hash:
        raise RepairError('source changed while preparing JES replacements')
    return PreparedJetRepair(replacements, dict(
        schema='CanonicalRetainedJetRepairV1', source_sha256=source_hash,
        status=('PARTIAL_RETAINED_JES_INVALID_JETS_MARKED'
                if invalid_ids else 'RETAINED_JET_PT_AND_PAIR_XJ_REPAIRED'),
        jets=len(jets), pairs=len(pairs), invalid_jes_jets=len(invalid_ids),
        invalid_jes_jet_ids=[list(x) for x in invalid_ids],
        invalid_jes_pair_ids=[list(x) for x in invalid_pair_ids],
        invalid_jes_policy='PRESERVE_ORIGINAL_JET_CORRECTED_PT_SET_PAIR_XJGAMMA_NAN',
        provenance=[json.loads(value) for value in sorted(provenance)],
        raw_pt_unchanged=True, truth_jets_unchanged=True, populations_unchanged=True,
        scientific_acceptance=False, capture_acceptance=False,
        scope='Retained objects only; missing threshold-crossing objects require capture evidence.'))


def _columns(rows, schema):
    return {name: np.asarray([row[name] for row in rows], dtype=kind)
            for name, kind in schema.items()}


def _column_digest(columns):
    h = hashlib.sha256()
    for name, values in sorted(columns.items()):
        h.update(name.encode() + b"\0" + values.dtype.str.encode() + b"\0")
        h.update(values.astype(values.dtype.newbyteorder("<"), copy=False).tobytes())
    return h.hexdigest()


@dataclass(frozen=True)
class RepairTables:
    tables: dict
    receipt: dict


def prepare_photon_tables(source, *, simulation, selection_name="all_retained_objects",
                          selected_candidate_ids=None, selected_truth_ids=None,
                          complete_reco_events=(), complete_truth_events=(),
                          max_rows=2_000_000, native_reco_capture_domain=None,
                          native_truth_capture=False):
    """Prepare exact keyed associations; no modification or event-order joins.

    Completeness collections are caller claims, not certification. Even when
    they are structurally complete, the receipt remains NOT_CERTIFIED. Selection
    IDs and completeness claims are fingerprinted in the response view identity.
    The default is an inventory diagnostic, NOT the physics signal definition.
    """
    if type(simulation) is not bool or not selection_name.strip():
        raise RepairError("explicit simulation flag and selection name required")
    if type(max_rows) is not int or max_rows <= 0:
        raise RepairError("positive bounded max_rows required")
    if type(native_truth_capture) is not bool:
        raise RepairError("native truth capture must be explicit boolean")
    source = Path(source).resolve()
    digest = file_sha256(source)
    part = (int(digest[:16], 16), int(digest[16:32], 16))
    with uproot.open(source) as root:
        events = _rows(root[BASE + "RJEventV1"], _ids("event_id", "source_occurrence_id"),
                       (list(RECO_CENSUS_FIELDS) if native_reco_capture_domain is not None else []) +
                       (list(TRUTH_CENSUS_FIELDS) if native_truth_capture else []), max_rows)
        sources = _rows(root[BASE + "RJSourceOccurrenceV1"], _ids("source_occurrence_id"), [], max_rows)
        source_ids = {_identity(r, "source_occurrence_id") for r in sources}
        if len(source_ids) != len(sources):
            raise RepairError("duplicate source occurrence")
        event_index = {}
        for event in events:
            eid, sid = _identity(event, "event_id"), _identity(event, "source_occurrence_id")
            if eid in event_index or sid not in source_ids:
                raise RepairError("duplicate event or missing source occurrence")
            event_index[eid] = sid
        candidates = _rows(root[BASE + "RJPhotonCandidateV1"], _ids("candidate_id", "event_id"), [
            "dominant_truth_state", "dominant_truth_evaluator_mode", "dominant_truth_track_id",
            "dominant_truth_embedding_id", "dominant_truth_pid", "dominant_truth_energy_contribution",
            "eta", "phi", "cluster_et", "native_cluster_key"], max_rows)
        truths = _rows(root[BASE + "RJTruthPhotonV1"], _ids("truth_photon_id", "event_id"), [
            "native_track_id", "embedding_id", "g4_photon_valid", "eta", "phi",
            "pt", "analysis_signal_r03"], max_rows) if simulation else []
        native_links = _native_link_provenance(root, max_rows)
    if file_sha256(source) != digest:
        raise RepairError("source changed during photon repair preparation")
    if any(_identity(r, "event_id") not in event_index for r in candidates + truths):
        raise RepairError("object event absent from explicit event inventory")
    reco_capture, truth_capture = set(complete_reco_events), set(complete_truth_events)
    if native_truth_capture:
        if not simulation or truth_capture:
            raise RepairError("native truth census requires simulation and no manual truth completeness override")
        truth_capture = native_truth_complete_events(events, truths)
    if native_reco_capture_domain is not None:
        if not simulation or reco_capture:
            raise RepairError("native reco census requires simulation and no manual reco completeness override")
        reco_capture = native_reco_complete_events(events, candidates, domain=native_reco_capture_domain)
    if (reco_capture | truth_capture) - event_index.keys():
        raise RepairError("capture claim references an absent event")
    selection = {"name": selection_name, "candidates": None if selected_candidate_ids is None else sorted(set(selected_candidate_ids)),
        "truths": None if selected_truth_ids is None else sorted(set(selected_truth_ids)),
        "complete_reco_events": sorted(reco_capture), "complete_truth_events": sorted(truth_capture)}
    if native_reco_capture_domain is not None:
        selection["native_reco_capture_domain"] = native_reco_capture_domain
    if native_truth_capture:
        selection["native_truth_capture"] = True
    selection_digest = hashlib.sha256(json.dumps(selection, sort_keys=True).encode()).hexdigest()
    result = rebuild_photon_links(candidates, truths, source_scope=digest, simulation=simulation,
        selected_candidate_ids=selection["candidates"], selected_truth_ids=selection["truths"],
        complete_reco_events=reco_capture, complete_truth_events=truth_capture, event_ids=set(event_index))
    reco_index = {_identity(r, "candidate_id"): r for r in candidates}
    truth_index = {_identity(r, "truth_photon_id"): r for r in truths}
    associations, links = [], []
    states = {}

    def identity_fields(row):
        value = {}
        for name, identity in (("part_id", part), ("source_occurrence_id", event_index[row["event_id"]]),
            ("event_id", row["event_id"]), ("candidate_id", row["candidate_id"]), ("truth_photon_id", row["truth_id"])):
            _put_id(value, name, identity)
        return value

    for row in result.associations:
        value = identity_fields(row)
        state = ASSOCIATION_STATES[row["state"]]
        states[row["candidate_id"]] = state
        dr = math.nan
        if row["truth_id"] is not None:
            reco, truth = reco_index[row["candidate_id"]], truth_index[row["truth_id"]]
            if all(k in reco and k in truth and math.isfinite(float(reco[k])) and math.isfinite(float(truth[k])) for k in ("eta", "phi")):
                dr = math.hypot(float(reco["eta"]) - float(truth["eta"]),
                    math.remainder(float(reco["phi"]) - float(truth["phi"]), 2 * math.pi))
        value.update(association_state=state, diagnostic_delta_r=dr,
                     dominant_energy_contribution=row.get("energy_contribution", math.nan))
        associations.append(value)
    seen_links = set()
    for row in result.links:
        value = identity_fields(row)
        token = json.dumps([VERSION, digest, selection_digest, row["event_id"], row["candidate_id"], row["truth_id"]])
        link_id = _hash_id(token)
        if link_id in seen_links:
            raise RepairError("duplicate repaired link identity")
        seen_links.add(link_id)
        _put_id(value, "link_id", link_id)
        value.update(link_class=LINK_CLASSES[row["link_class"]], association_state=states.get(row["candidate_id"], 0))
        links.append(value)
    tables = {"photonAssociations": _columns(associations, ASSOCIATION_SCHEMA),
              "photonResponseLinks": _columns(links, LINK_SCHEMA)}
    receipt = {"schema": VERSION, "source_sha256": digest, "source_path": str(source),
        "simulation": simulation, "selection_name": selection_name, "selection_sha256": selection_digest,
        "status": "NOT_CERTIFIED", "association_complete": result.association_complete,
        "response_structurally_complete": result.response_complete,
        "capture_authority": "CALLER_CLAIM_ONLY_REQUIRES_INDEPENDENT_ACCEPTANCE",
        "native_truth_capture": {"enabled": native_truth_capture,
            "verified_complete_events": len(truth_capture) if native_truth_capture else None,
            "scope": "NATIVE_EMBEDDED_PRIMARY_PID22_CENSUS_NO_REJECTED_ROWS_NOT_PHYSICS_ACCEPTANCE"},
        "native_reco_capture": {"domain": native_reco_capture_domain,
            "verified_complete_events": len(reco_capture) if native_reco_capture_domain is not None else None,
            "scope": "RECO_CONTAINER_DOMAIN_ONLY_NOT_TRUTH_OR_PHYSICS_ACCEPTANCE"},
        "association_method": "CLUSTER_MAXIMUM_ENERGY_PRIMARY_IDENTITY",
        "diagnostic_delta_r_semantics": "DIAGNOSTIC_ONLY_NOT_A_MATCH_CUT",
        "native_photon_links": native_links,
        "legacy_photon_links": ("PRESERVED_HISTORICAL_NOT_CANONICAL"
            if native_links["photon_policy_counts"]["LEGACY_UNSPECIFIED"] else "ABSENT"),
        "jet_links": "UNCHANGED", "association_states": ASSOCIATION_STATES, "link_classes": LINK_CLASSES,
        "required_donor_inputs": list(result.required_donor_inputs),
        "association_counts": dict(Counter(r["state"] for r in result.associations)),
        "link_counts": dict(Counter(r["link_class"] for r in result.links)),
        "tables": {name: {"rows": len(next(iter(columns.values()))), "schema": SCHEMAS[name],
                           "sha256": _column_digest(columns)} for name, columns in tables.items()}}
    return RepairTables(tables, receipt)


def write_scalar_tables(directory, prepared, *, manifest_name, progress=None, max_seconds=120):
    """Write only new repair tables to an already-open staged ROOT TDirectory.

    Caller owns atomic publication/cleanup. Existing keys fail before any write;
    use verify_root_tables on the closed stage before publishing.
    """
    import ROOT
    if not math.isfinite(max_seconds) or max_seconds <= 0:
        raise RepairError("finite positive write deadline required")
    started, processed = time.monotonic(), 0
    total = sum(r["rows"] for r in prepared.receipt["tables"].values())

    def update(phase):
        elapsed = time.monotonic() - started
        if elapsed > max_seconds:
            raise RepairError("repair table write deadline exceeded")
        if progress:
            progress(dict(phase=phase, processed=processed, total=total,
                          remaining=total-processed, elapsed_seconds=elapsed))

    update("write_repairs")
    schemas = {name: record["schema"] for name, record in prepared.receipt["tables"].items()}
    if set(schemas) != set(prepared.tables):
        raise RepairError("prepared table inventory mismatch")
    if any(directory.Get(name) for name in (*schemas, manifest_name)):
        raise RepairError("refusing to replace existing repair objects")
    if file_sha256(prepared.receipt["source_path"]) != prepared.receipt["source_sha256"]:
        raise RepairError("source changed before repair write")
    buffers_by_type = {"uint64": ("Q", "l"), "uint32": ("I", "i"),
                       "int64": ("q", "L"), "int32": ("i", "I"), "float64": ("d", "D")}
    for name, schema in schemas.items():
        directory.cd()
        tree = ROOT.TTree(name, prepared.receipt["schema"])
        buffers = {}
        for field, kind in schema.items():
            code, leaf = buffers_by_type[kind]
            buffers[field] = array(code, [0])
            tree.Branch(field, buffers[field], f"{field}/{leaf}")
        columns = prepared.tables[name]
        if _column_digest(columns) != prepared.receipt["tables"][name]["sha256"]:
            raise RepairError("prepared repair columns changed")
        for i in range(prepared.receipt["tables"][name]["rows"]):
            for field in schema:
                buffers[field][0] = columns[field][i].item()
            tree.Fill()
            processed += 1
            if processed % 1024 == 0:
                update("write_repairs")
        tree.Write()
        tree.SetDirectory(0)
    directory.cd()
    ROOT.TObjString(json.dumps(prepared.receipt, sort_keys=True)).Write(manifest_name)
    update("repair_tables_written")
    return prepared.receipt


def verify_scalar_tables(path, prepared, *, manifest_name, prefix=""):
    """Full independent scalar readback, including NaN and uint64 bit patterns."""
    with uproot.open(path) as root:
        if json.loads(str(root[prefix + manifest_name])) != prepared.receipt:
            raise RepairError("repair manifest mismatch")
        for name, record in prepared.receipt["tables"].items():
            schema = record["schema"]
            tree = root[prefix + name]
            if tree.classname != "TTree" or set(tree.keys()) != set(schema):
                raise RepairError("repair physical schema mismatch")
            arrays = tree.arrays(library="np")
            if any(arrays[n].dtype != np.dtype(kind) for n, kind in schema.items()):
                raise RepairError("repair physical type mismatch")
            if _column_digest(arrays) != prepared.receipt["tables"][name]["sha256"]:
                raise RepairError("repair value readback mismatch")
    return {"status": "PASS_REPAIR_TABLE_READBACK", "science_acceptance": False}


def write_root_tables(directory, prepared, *, progress=None, max_seconds=120):
    return write_scalar_tables(directory, prepared, manifest_name="photon_repair_manifest",
                               progress=progress, max_seconds=max_seconds)


def verify_root_tables(path, prepared, *, prefix=""):
    return verify_scalar_tables(path, prepared, manifest_name="photon_repair_manifest", prefix=prefix)
