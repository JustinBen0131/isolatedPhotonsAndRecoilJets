#!/usr/bin/env python3
"""Read-only audit and additive, keyed QA recovery from completed schema14/v2.

No row-order friend joins, reconstruction, calibration changes, or base edits.
Only a finite-positive calibrated PMT charge sum is recoverable here. Timing,
calorimeter energies, and missing scores remain explicit recovery requirements.
"""
from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
import sqlite3
import tempfile

import awkward as ak
import numpy as np
import uproot

from centrality_replay import digest

VERSION = "Schema14QAAugmentationV1"
EVENT_NATIVE = "ReplayFoundationV1/RJEventV1"
PHOTON_NATIVE = "ReplayFoundationV1/RJPhotonCandidateV1"
IDENTITY = ("source_occurrence_id_hi", "source_occurrence_id_lo", "event_id_hi", "event_id_lo")
CALO = ("emcal_total_energy", "ihcal_total_energy", "ohcal_total_energy", "total_calo_energy")
COPIED = ("run", "physical_event_sequence", "physical_event_sequence_valid")


def layout(root):
    native = EVENT_NATIVE in root
    tree = EVENT_NATIVE if native else "events"
    keys = IDENTITY + (() if native else ("source_file_index",))
    missing = set(keys) - set(root[tree].keys())
    if missing:
        raise ValueError(f"Missing identity columns: {missing}")
    return tree, keys


def key_rows(arrays, keys):
    # Decimal strings preserve uint64 values in sqlite without signed overflow.
    return [": ".join(str(int(x)) for x in row)
            for row in zip(*(ak.to_numpy(arrays[k]) for k in keys))]


def copied_rows(arrays, fields):
    if not fields:
        return ["[]"] * len(arrays)
    return [json.dumps([int(x) for x in row])
            for row in zip(*(ak.to_numpy(arrays[k]) for k in fields))]


def positive_charge(arrays):
    """Exact AuAu RecoilJets sum: finite q>0, no centrality/time threshold.

    Require a complete, identity-valid calibrated PMT array before recovering a
    previously missing total. Do not substitute centrality_pmt_charge: embedding
    can bind that to a different node.
    """
    required = ("mbd_pmt_available", "mbd_pmt_id", "mbd_pmt_charge", "mbd_pmt_charge_valid")
    if any(k not in arrays.fields for k in required):
        return np.full(len(arrays), np.nan), np.zeros(len(arrays), dtype=bool)
    q, ids, valid = arrays.mbd_pmt_charge, arrays.mbd_pmt_id, arrays.mbd_pmt_charge_valid
    if not (ak.all(ak.num(q) == ak.num(ids)) and ak.all(ak.num(q) == ak.num(valid))):
        raise ValueError("PMT lengths disagree")
    complete = ak.to_numpy((arrays.mbd_pmt_available == 1) & (ak.num(q) == 128)
                          & ak.all(ids == ak.local_index(ids), axis=1)
                          & ak.all((valid == 1) & np.isfinite(q), axis=1))
    summed = ak.to_numpy(ak.sum(ak.where(np.isfinite(q) & (q > 0), q, 0), axis=1))
    return np.where(complete, summed, np.nan), complete


def make_columns(a, keys):
    derived, complete = positive_charge(a)
    old = ak.to_numpy(a.mbd_total_charge)
    finite = np.isfinite(old)
    comparable = finite & complete
    if np.any(old[comparable] != derived[comparable]):
        raise ValueError("Saved total disagrees with exact finite-positive PMT sum")
    state = np.where(finite, 1, np.where(complete, 2, 0)).astype("int32")
    columns = {k: ak.to_numpy(a[k]) for k in keys}
    columns.update(mbd_total_charge=np.where(finite, old, derived),
                   mbd_total_charge_state=state)
    for name in COPIED:
        if name in a.fields:
            columns[name] = ak.to_numpy(a[name])
    return columns, comparable


def build(source, output, step_size=10000):
    source, output = Path(source).resolve(), Path(output).resolve()
    before = digest(source)
    output.mkdir(parents=True, exist_ok=False)
    friend_path = output / "qa_friend.root"
    receipt = {"schema": VERSION, "source": str(source), "source_sha256": before,
               "status": "INCOMPLETE", "states": {"0": "UNAVAILABLE", "1": "NATIVE",
               "2": "RECOVERED_EXACT_PMT_SUM"}, "base_modified": False}
    stats, runs, by_status = Counter(), Counter(), Counter()
    with uproot.open(source) as root, uproot.recreate(friend_path) as dest, \
            tempfile.TemporaryDirectory() as tmp, sqlite3.connect(str(Path(tmp)/"keys.sqlite")) as db:
        tree, keys = layout(root)
        receipt.update(event_tree=tree, keys=list(keys), root_uuid=str(root.file.uuid))
        metadata = {}
        if "source_files" in root:
            metadata["source_files"] = json.loads(str(root["source_files"]))
        elif "ReplayFoundationV1/RJSourceOccurrenceV1" in root:
            metadata["source_occurrences"] = ak.to_list(root["ReplayFoundationV1/RJSourceOccurrenceV1"].arrays())
        for name in ("analysis_config_yaml",):
            if name in root:
                metadata[name] = str(root[name])
        receipt["source_identity"] = metadata
        schema = {k: "uint64" for k in keys}
        if "source_file_index" in schema:
            schema["source_file_index"] = "int32"
        for k in ("run", "physical_event_sequence", "physical_event_sequence_valid"):
            if k in root[tree].keys():
                schema[k] = "int64"
        schema.update(mbd_total_charge="float64", mbd_total_charge_state="int32")
        friend = dest.mktree("eventQA", schema)
        db.execute("CREATE TABLE seen (identity TEXT PRIMARY KEY)")
        wanted = list(keys) + [k for k in root[tree].keys() if k in {
            "run", "physical_event_sequence", "physical_event_sequence_valid",
            "terminal_status", "candidate_count", "mbd_total_charge", "mbd_pmt_available",
            "mbd_pmt_id", "mbd_pmt_charge", "mbd_pmt_charge_valid", *CALO}]
        for a in root[tree].iterate(wanted, step_size=step_size, library="ak"):
            columns, comparable = make_columns(a, keys)
            try:
                db.executemany("INSERT INTO seen VALUES (?)", [(x,) for x in key_rows(a, keys)])
            except sqlite3.IntegrityError as exc:
                raise ValueError("Duplicate source-scoped event identity") from exc
            friend.extend(columns)
            stats.update({"events": len(a), "existing_totals_verified": int(comparable.sum()),
                "totals_recovered": int(np.sum(columns["mbd_total_charge_state"] == 2)),
                "totals_unavailable": int(np.sum(columns["mbd_total_charge_state"] == 0))})
            runs.update(int(x) for x in a.run)
            for field in CALO:
                if field in a.fields:
                    stats[f"{field}_missing"] += int(ak.sum(~np.isfinite(a[field])))
            if "terminal_status" in a.fields:
                for status, missing in zip(ak.to_list(a.terminal_status),
                                            ak.to_list(~np.isfinite(a.mbd_total_charge))):
                    by_status[f"status_{status}_events"] += 1
                    by_status[f"status_{status}_missing_total"] += int(missing)
            if "candidate_count" in a.fields:
                stats["missing_total_events_with_candidates"] += int(ak.sum(
                    ~np.isfinite(a.mbd_total_charge) & (a.candidate_count > 0)))
        photon_tree = PHOTON_NATIVE if PHOTON_NATIVE in root else "photons"
        if photon_tree in root:
            for field in ("npb_score", "auau_npb_score"):
                if field in root[photon_tree].keys():
                    counts = Counter()
                    for a in root[photon_tree].iterate([field], library="np", step_size=step_size):
                        counts["rows"] += len(a[field]); counts["finite"] += int(np.isfinite(a[field]).sum())
                    receipt[field] = dict(counts)
        receipt["timing_branches"] = {
            k: [b for b in root[k].keys() if "time" in b.lower()]
            for k in (tree, photon_tree) if k in root}
    if digest(source) != before:
        raise ValueError("Source changed during augmentation")
    receipt.update(counts=dict(stats), by_terminal_status=dict(by_status),
                   runs={str(k):v for k,v in sorted(runs.items())},
                   friend_sha256=digest(friend_path), status="READY_FOR_READBACK")
    (output/"receipt.json").write_text(json.dumps(receipt, indent=2)+"\n")
    verify(source, output, step_size)
    receipt["status"] = "PASS_EXACT_CHARGE_AUGMENTATION"
    (output/"receipt.json").write_text(json.dumps(receipt, indent=2)+"\n")
    return receipt


def verify(source, output, step_size=10000):
    """Join by source-scoped key even when friend rows are arbitrarily reordered."""
    source, output = Path(source).resolve(), Path(output)
    r = json.loads((output/"receipt.json").read_text())
    if r["source"] != str(source) or digest(source) != r["source_sha256"]:
        raise ValueError("Source identity or content changed")
    friend_path = output/"qa_friend.root"
    if digest(friend_path) != r["friend_sha256"]:
        raise ValueError("Friend content changed")
    keys = r["keys"]
    with tempfile.TemporaryDirectory() as tmp, sqlite3.connect(str(Path(tmp)/"join.sqlite")) as db, \
            uproot.open(source) as src, uproot.open(friend_path) as friend:
        copied = [k for k in COPIED if k in src[r["event_tree"]].keys()]
        db.execute("CREATE TABLE rows (identity TEXT PRIMARY KEY, val REAL, state INTEGER, copied TEXT)")
        for a in friend["eventQA"].iterate(library="ak", step_size=step_size):
            try:
                db.executemany("INSERT INTO rows VALUES (?,?,?,?)", [
                    (k, float(v) if np.isfinite(v) else None, int(s), metadata)
                    for k,v,s,metadata in zip(key_rows(a,keys), ak.to_numpy(a.mbd_total_charge),
                                             ak.to_numpy(a.mbd_total_charge_state), copied_rows(a,copied))])
            except sqlite3.IntegrityError as exc:
                raise ValueError("Duplicate friend identity") from exc
        seen = 0
        wanted = keys + [k for k in src[r["event_tree"]].keys() if k in {
            "mbd_total_charge", "mbd_pmt_available", "mbd_pmt_id", "mbd_pmt_charge", "mbd_pmt_charge_valid", *COPIED}]
        for a in src[r["event_tree"]].iterate(wanted, library="ak", step_size=step_size):
            columns, _ = make_columns(a, keys)
            for k,v,s,metadata in zip(key_rows(a,keys), columns["mbd_total_charge"],
                                     columns["mbd_total_charge_state"], copied_rows(a,copied)):
                actual = db.execute("SELECT val,state,copied FROM rows WHERE identity=?",(k,)).fetchone()
                expected = (float(v) if np.isfinite(v) else None, int(s), metadata)
                if actual != expected:
                    raise ValueError(f"Keyed value mismatch: {k}")
                db.execute("DELETE FROM rows WHERE identity=?",(k,)); seen += 1
        if db.execute("SELECT COUNT(*) FROM rows").fetchone()[0] or seen != r["counts"]["events"]:
            raise ValueError("Unmatched friend/source rows")
    return {"status":"PASS", "rows":seen, "join":"source-bound identity, not entry order"}


def run_index(source, output):
    """Non-copying run index for an immutable mixed-run v2 part."""
    source, output = Path(source).resolve(), Path(output)
    if output.exists(): raise FileExistsError(output)
    before = digest(source)
    with uproot.open(source) as f:
        tree, keys = layout(f)
        a=f[tree].arrays(["run"],library="np")
        runs={}
        for run in sorted(set(map(int,a["run"]))):
            entries=np.flatnonzero(a["run"]==run)
            breaks=np.flatnonzero(np.diff(entries)!=1)+1
            ranges=[[int(x[0]),int(x[-1])+1] for x in np.split(entries,breaks)]
            runs[str(run)]={"events":len(entries),"entry_ranges":ranges}
    if digest(source)!=before:raise ValueError("Source changed")
    result={"schema":"Schema14RunIndexV1","source":str(source),"source_sha256":before,
        "event_tree":tree,"runs":runs,"selection":"Use event identities to select photons/jets; ranges apply only to event tree.",
        "run_completeness":"Not asserted: requires full production input/receipt reconciliation."}
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,indent=2)+"\n")
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("action",choices=("build","verify","run-index"));p.add_argument("--source",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True);a=p.parse_args()
    result={"build":build,"verify":verify,"run-index":run_index}[a.action](a.source,a.output)
    print(json.dumps({k:v for k,v in result.items() if k in ("status","counts","runs","rows")},indent=2))


if __name__=="__main__":main()
