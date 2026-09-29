"""Offline AuAu DATA continuation/refill manifest; never execution authority.

Call plan_continuation first to obtain jobs, bind explicit output paths with
job_sha256, then build_manifest with the same inputs. QueueSnapshot counters
must be disjoint per host (all owner states, including other campaigns).
The caller verifies source/calibration bytes, receipts, ledger completeness,
output namespace availability, and later obtains executor admission separately.
No files, clocks, scheduler APIs, or submission backends are accessed here.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
import re
from typing import Any, Iterable, Mapping, Sequence

try:
    from . import production_contract as c
except ImportError:  # existing standalone worker import convention
    import production_contract as c


PLANNED_HOSTS = {
    "auau_data": "sphnxuser04.sdcc.bnl.gov",
    "pp_data": "sphnxuser06.sdcc.bnl.gov",
    "pp_sim": "sphnxuser03.sdcc.bnl.gov",
    "auau_sim": "sphnxuser05.sdcc.bnl.gov",
}


def canonical_bytes(value: Any) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":"),
                       allow_nan=False) + "\n").encode("utf-8")


def _digest(value: Any) -> str:
    return hashlib.sha256(canonical_bytes(value)).hexdigest()


def job_sha256(job: c.PlannedJob) -> str:
    """Bind exact range, triplet, timing evidence, and planning envelope."""
    if not isinstance(job, c.PlannedJob):
        raise c.ContractError("typed planned job required")
    return _digest(asdict(job))


@dataclass(frozen=True)
class OutputIdentity:
    source_range: c.SourceRange
    job_sha256: str
    schema14_path: str
    collaborator_v2_path: str
    replacement: bool = False

    def __post_init__(self) -> None:
        if not isinstance(self.source_range, c.SourceRange):
            raise c.ContractError("output requires exact source range")
        c._sha(self.job_sha256, "output job binding")
        c._path(self.schema14_path, "schema14 output")
        c._path(self.collaborator_v2_path, "collaborator v2 output")
        if self.schema14_path == self.collaborator_v2_path:
            raise c.ContractError("output paths must be distinct")
        if type(self.replacement) is not bool:
            raise c.ContractError("replacement must be explicit boolean")


@dataclass(frozen=True)
class SourceCountCorrection:
    """An explicitly measured shorter source, never an inferred EOF repair.

    Normalize total metadata on every sibling, but shorten only a range that
    included the nonexistent tail. Original ordinals and global coordinates
    remain immutable. Evidence authentication belongs to the packet builder.
    """
    run: int
    segment: int
    jets: str
    jetcalo: str
    source_file_ordinal: int
    original_total_events: int
    measured_total_events: int
    evidence_sha256: str

    def __post_init__(self) -> None:
        c._integer(self.run, 'correction run', 1)
        c._integer(self.segment, 'correction segment')
        c._integer(self.source_file_ordinal, 'correction source ordinal', 1)
        c._integer(self.original_total_events, 'original source total', 1)
        c._integer(self.measured_total_events, 'measured source total', 1)
        c._path(self.jets, 'correction jets')
        c._path(self.jetcalo, 'correction jetcalo')
        c._sha(self.evidence_sha256, 'source count measurement')
        if self.jets == self.jetcalo or self.measured_total_events >= self.original_total_events:
            raise c.ContractError('correction requires a distinct pair and measured shorter source')


def normalize_source_record(record: Mapping[str, Any],
                            corrections: Sequence[SourceCountCorrection] = ()) -> dict[str, Any]:
    """Apply only an exact, evidence-bound source-total correction to a copy."""
    row = c.range_from_manifest(record)
    result = dict(record)
    matched = []
    for correction in corrections:
        if not isinstance(correction, SourceCountCorrection):
            raise c.ContractError('typed source count correction required')
        if row.jets == correction.jets or row.jetcalo == correction.jetcalo:
            matched.append(correction)
    if len(matched) > 1:
        raise c.ContractError('duplicate source count correction')
    if not matched:
        return result
    fix = matched[0]
    if ((row.run, row.segment, row.jets, row.jetcalo, row.source_total_events,
         record.get('source_file_ordinal')) !=
        (fix.run, fix.segment, fix.jets, fix.jetcalo, fix.original_total_events,
         fix.source_file_ordinal)):
        raise c.ContractError('source count correction identity mismatch')
    if row.event_offset >= fix.measured_total_events:
        raise c.ContractError('correction would remove an entire range')
    result['source_total_events'] = fix.measured_total_events
    result['event_count'] = min(row.end, fix.measured_total_events) - row.event_offset
    c.range_from_manifest(result)
    return result


def build_manifest(
    universe: Sequence[c.SourceRange], products: Sequence[c.ProductEvidence],
    calibration: c.CalibrationManifest, *,
    planning_parameters: Mapping[str, Any], outputs: Sequence[OutputIdentity],
    host_snapshots: Mapping[str, c.QueueSnapshot],
    previously_submitted: Sequence[c.SourceRange], submission_ledger_sha256: str,
    now: float, reserved_pp_sim_slots: int,
    max_new_jobs: int = c.MAX_OCCUPANCY,
    max_snapshot_age_seconds: float = 300,
) -> dict[str, Any]:
    """Return a canonical hash-bound PLAN_ONLY envelope for a finite universe.

    All four FQDN hosts are required, even when empty. The oldest observation
    bounds aggregate freshness; every host is checked for future timestamps.
    The ledger includes drained/failed/held attempts, not only visible jobs.
    Output paths are mandatory for every planned job, including submitted jobs;
    KEEP_VALID and UNVERIFIED rows may not acquire new output identities here.
    """
    c._sha(submission_ledger_sha256, "submission ledger evidence")
    c._number(now, "now")
    c._number(max_snapshot_age_seconds, "snapshot age", 1)
    if max_snapshot_age_seconds > 300:
        raise c.ContractError("snapshot age cannot exceed 300 seconds")
    if set(host_snapshots) != set(PLANNED_HOSTS.values()):
        raise c.ContractError("exactly the four owner submit hosts are required")
    for snapshot in host_snapshots.values():
        if not isinstance(snapshot, c.QueueSnapshot):
            raise c.ContractError("typed per-host queue snapshot required")
        if not 0 <= now - snapshot.observed_at <= max_snapshot_age_seconds:
            raise c.ContractError("host queue snapshot is stale or from the future")
    rows = c.validate_ranges(universe)
    submitted = c.validate_ranges(previously_submitted)
    plan = c.plan_continuation(rows, products, calibration, **planning_parameters)
    host_evidence = {host: asdict(host_snapshots[host]) for host in sorted(host_snapshots)}
    aggregate = c.QueueSnapshot(
        observed_at=min(s.observed_at for s in host_snapshots.values()),
        evidence_sha256=_digest(host_evidence),
        **{name: sum(getattr(s, name) for s in host_snapshots.values())
           for name in ("idle", "running", "transferring", "held", "suspended", "other_queued")},
        inflight_ranges=tuple(r for host in sorted(host_snapshots)
                              for r in host_snapshots[host].inflight_ranges),
        owner_wide_complete=True,
    )
    refill = c.propose_replenishment(
        plan.jobs, aggregate, previously_submitted=submitted, now=now,
        max_new_jobs=max_new_jobs, max_snapshot_age_seconds=max_snapshot_age_seconds,
        reserved_other_lane_slots=reserved_pp_sim_slots,
    )
    if any(not isinstance(o, OutputIdentity) for o in outputs):
        raise c.ContractError("typed output identities required")
    c.validate_ranges(tuple(o.source_range for o in outputs))
    by_range = {o.source_range.identity: o for o in outputs}
    if set(by_range) != {j.source_range.identity for j in plan.jobs}:
        raise c.ContractError("outputs must exactly cover planned jobs")
    paths = [p for o in outputs for p in (o.schema14_path, o.collaborator_v2_path)]
    if len(paths) != len(set(paths)):
        raise c.ContractError("duplicate output path")
    states = {d.source_range.identity: d.state for d in plan.coverage.decisions}
    selected = {j.source_range.identity for j in refill.jobs}
    jobs = []
    for job in plan.jobs:
        output = by_range[job.source_range.identity]
        if output.job_sha256 != job_sha256(job):
            raise c.ContractError("stale or mixed-basis output job binding")
        if output.replacement != (states[job.source_range.identity] == "REPLACE_CENTRALITY"):
            raise c.ContractError("replacement output identity must match coverage decision")
        jobs.append({"job": asdict(job), "output": asdict(output),
                     "centrality_environment": c.centrality_environment(job, calibration),
                     "planned_host": PLANNED_HOSTS["auau_data"],
                     "selected_for_refill": job.source_range.identity in selected})
    payload = {
        "schema": "ContinuationManifestV1", "state": "PLAN_ONLY",
        "execution_authorized": False, "scientific_acceptance": False,
        "independent_pool_priority": False, "planned_hosts": dict(PLANNED_HOSTS),
        "lane": "auau_data", "coverage": asdict(plan.coverage), "jobs": jobs,
        "replenishment": asdict(refill), "host_snapshots": host_evidence,
        "calibration": {"set_id": calibration.set_id,
            "manifest_sha256": calibration.manifest_sha256,
            "triplets": [asdict(t) for t in sorted(calibration.triplets, key=lambda t: t.run)],
            "fallbacks": [asdict(f) for f in sorted(calibration.fallbacks, key=lambda f: f.run)]},
        "products": [asdict(p) for p in sorted(products, key=lambda p: p.source_range.identity)],
        "submission_ledger": {"sha256": submission_ledger_sha256,
                              "ranges": [asdict(r) for r in submitted]},
        "planning_parameters": dict(planning_parameters),
        "refill_parameters": {"now": now, "max_new_jobs": max_new_jobs,
            "max_snapshot_age_seconds": max_snapshot_age_seconds,
            "reserved_pp_sim_slots": reserved_pp_sim_slots},
        "limits": {"total_cap": c.MAX_OCCUPANCY, "low_water": c.LOW_WATER,
                   "emergency_water": c.EMERGENCY_WATER},
    }
    return {"payload": payload, "sha256": _digest(payload)}


def verify_manifest(envelope: Mapping[str, Any]) -> None:
    """Detect serialization corruption; this is not provenance authentication."""
    if set(envelope) != {"payload", "sha256"} or _digest(envelope["payload"]) != envelope["sha256"]:
        raise c.ContractError("manifest hash mismatch")
    payload = envelope["payload"]
    if (payload.get("schema") != "ContinuationManifestV1" or payload.get("state") != "PLAN_ONLY"
            or payload.get("execution_authorized") is not False):
        raise c.ContractError("manifest is not a PLAN_ONLY proposal")


def bind_capture_rows(envelope: Mapping[str, Any], source_lines: Iterable[bytes], *,
                      expected_source_sha256: str,
                      source_count_corrections: Sequence[SourceCountCorrection] = ()) -> dict[str, Any]:
    """Stream an original ranges.jsonl into exact selected capture bindings.

    Keep the original source ordinal/global coordinate, never enumerate the
    selected subset. Hash every original byte while retaining only selected
    rows. The caller supplies the approved original manifest hash; this does
    not authenticate that approval or admit a packet, payload or submission.
    The existing DATA worker remains the executor; no shell is generated here.
    """
    verify_manifest(envelope)
    c._sha(expected_source_sha256, 'original source manifest hash')
    selected = {}
    for entry in envelope['payload']['jobs']:
        if type(entry['selected_for_refill']) is not bool:
            raise c.ContractError('explicit selected-for-refill boolean required')
        if not entry['selected_for_refill']:
            continue
        row = c.SourceRange(**entry['job']['source_range'])
        if row.identity in selected:
            raise c.ContractError('duplicate selected source identity')
        selected[row.identity] = entry
    if len(selected) > c.MAX_OCCUPANCY:
        raise c.ContractError('capture binding exceeds finite refill cap')
    digest = hashlib.sha256()
    found, ordinals, pairs = {}, {}, {}
    for line in source_lines:
        if not isinstance(line, bytes):
            raise c.ContractError('source records must be original byte lines')
        digest.update(line)
        if not line.strip():
            continue
        original_record = json.loads(line)
        record = normalize_source_record(original_record, source_count_corrections)
        row = c.range_from_manifest(record)
        if row.identity not in selected:
            continue
        if row.identity in found:
            raise c.ContractError('duplicate original source range/alias')
        ordinal = record.get('source_file_ordinal')
        global_begin = record.get('source_global_entry_begin')
        c._integer(ordinal, 'original source ordinal', 1)
        c._integer(global_begin, 'original global entry')
        base = global_begin - row.event_offset
        if base < 0:
            raise c.ContractError('global source coordinate precedes event offset')
        if ordinal in ordinals and ordinals[ordinal] != row.pair:
            raise c.ContractError('original ordinal aliases different source pairs')
        if row.pair in pairs and pairs[row.pair] != (ordinal, base):
            raise c.ContractError('inconsistent coordinates within one source pair')
        ordinals[ordinal] = row.pair
        pairs[row.pair] = (ordinal, base)
        entry = selected[row.identity]
        environment = dict(entry['centrality_environment'])
        environment.update(RJ_REPLAY_RUN=str(row.run), RJ_REPLAY_SEGMENT=str(row.segment),
            RJ_EVENT_OFFSET=str(row.event_offset), RJ_REQUIRE_ORIGINAL_SOURCE_CURSOR='1',
            RJ_REPLAY_SOURCE_FILE_ORDINAL=str(ordinal),
            RJ_REPLAY_SOURCE_ENTRY_BEGIN=str(row.event_offset),
            RJ_REPLAY_SOURCE_GLOBAL_ENTRY_BEGIN=str(global_begin))
        found[row.identity] = dict(planned_row_id=entry['job']['source_range']['row_id'],
            original_row=original_record, effective_row=record,
            source_coordinate_environment=environment,
            output=dict(entry['output']), planned_host=entry['planned_host'],
            planned_job_sha256=_digest(entry['job']))
    if digest.hexdigest() != expected_source_sha256:
        raise c.ContractError('original source manifest bytes changed')
    if set(found) != set(selected):
        raise c.ContractError('selected source range missing from original manifest')
    payload = dict(schema='ContinuationCaptureRowsV1', state='PLAN_ONLY',
        execution_authorized=False, scientific_acceptance=False,
        continuation_manifest_sha256=envelope['sha256'],
        original_ranges_sha256=expected_source_sha256,
        source_count_corrections=[asdict(fix) for fix in source_count_corrections],
        rows=[found[key] for key in sorted(found)],
        note='Packet/runtime hashes and source semantic provenance remain the existing worker sealer responsibility; this is not a submission descriptor.')
    return dict(payload=payload, sha256=_digest(payload))


def _validate_pp_packet_ranges(rows: Sequence[Mapping[str, Any]], sample_id: str) -> None:
    """Validate pp geometry without routing it through AuAu centrality planning."""
    seen_ids, metadata, intervals = set(), {}, {}
    for row in rows:
        try:
            rid = row['row_id']
            c._text(rid, 'row_id')
            if rid in seen_ids or row['sample_id'] != sample_id:
                raise c.ContractError('duplicate row or pp sample identity differs')
            seen_ids.add(rid)
            for key, minimum in (('run',1), ('segment',0), ('source_total_events',1),
                                 ('event_offset',0), ('event_count',1),
                                 ('source_file_ordinal',1), ('source_global_entry_begin',0)):
                c._integer(row[key], key, minimum)
            if row['event_count'] > 100000:
                raise c.ContractError('pp event count exceeds row ceiling')
            start, end = row['event_offset'], row['event_offset'] + row['event_count']
            if end > row['source_total_events']:
                raise c.ContractError('range exceeds exact source_total_events')
            sources = row['source_paths']
            if set(sources) != {'jets', 'jetcalo'}:
                raise c.ContractError('exact pp paired sources required')
            pair = tuple(sources[key] for key in ('jets','jetcalo'))
            for path in pair: c._path(path, 'pp paired source')
            if pair[0] == pair[1]:
                raise c.ContractError('paired sources must be distinct')
            identity = (row['run'], row['segment'], row['source_total_events'],
                        row['source_file_ordinal'], row['source_global_entry_begin']-start, pair)
            for path in pair:
                if path in metadata and metadata[path] != identity:
                    raise c.ContractError('source reused with conflicting identity or coordinates')
                metadata[path] = identity
            intervals.setdefault(pair, []).append((start,end))
        except (KeyError,TypeError) as exc:
            raise c.ContractError('incomplete pp paired-source atom') from exc
    for group in intervals.values():
        ordered = sorted(group)
        if any(right[0] < left[1] for left,right in zip(ordered,ordered[1:])):
            raise c.ContractError('duplicate/overlapping pp ranges')


def compile_data_worker_packet(bound_rows: Mapping[str, Any], *, runtime_manifest: Mapping[str, Any],
                               packet_root: str, sample_id: str,
                               augmentation_templates: Mapping[str, Mapping[str, str]],
                               receipt_paths: Mapping[str, str], request_memory_mb: int = 4096,
                               max_packet_rows: int = 1000, data_lane: str = 'auau_data',
                               submit_host: str | None = None) -> dict[str, Any]:
    """Compile selected paired DATA rows into the EXISTING worker's 20 arguments.

    No scheduler or transfer is called. The returned immutable byte payloads
    are staged by the existing executor only after its separate admission.
    Small packet shards bound per-job manifest hashing; no million-row manifest
    is copied/read by every worker. Output is the canonical v2 path; native
    capture stays in execute-node scratch under the existing DATA worker.
    """
    if (set(bound_rows) != {'payload','sha256'}
            or _digest(bound_rows['payload']) != bound_rows['sha256']
            or bound_rows['payload'].get('schema') != 'ContinuationCaptureRowsV1'
            or bound_rows['payload'].get('execution_authorized') is not False):
        raise c.ContractError('hash-bound capture rows required')
    c._path(packet_root, 'packet root')
    c._integer(request_memory_mb, 'request memory MB', 1)
    c._integer(max_packet_rows, 'packet row limit', 1)
    if request_memory_mb > 4096 or max_packet_rows > 1000:
        raise c.ContractError('packet resource/shard bound exceeded')
    selected = bound_rows['payload']['rows']
    if not selected or len(selected) > max_packet_rows:
        raise c.ContractError('nonempty bounded packet shard required')
    manifest = json.loads(canonical_bytes(runtime_manifest))
    if (manifest.get('schema') != 'THE327Schema14DATAProductionManifestV1'
            or manifest.get('status') != 'PASS_FROZEN_UNSUBMITTED'
            or manifest.get('submission_ready') is not True
            or manifest.get('source_schema_version') != 14
            or set(manifest.get('samples', {})) != {sample_id}
            or not manifest.get('runtime', {}).get('canonical_assembly')):
        raise c.ContractError('sealed canonical DATA runtime manifest required')
    sample = manifest['samples'][sample_id]
    if data_lane not in ('auau_data','pp_data'):
        raise c.ContractError('explicit supported DATA lane required')
    planned_host = PLANNED_HOSTS[data_lane] if submit_host is None else submit_host
    if not isinstance(planned_host, str) or not re.fullmatch(r'sphnxuser0[1-8]\.sdcc\.bnl\.gov', planned_host):
        raise c.ContractError('approved exact DATA submit host required')
    system, dataset = ('auau','isAuAu') if data_lane == 'auau_data' else ('pp','isPP')
    if (sample.get('system'), sample.get('dataset'), sample.get('si_di_role'), sample.get('source_roles')) != (
            system, dataset, 'DATA', ['jets','jetcalo']):
        raise c.ContractError('DATA lane and paired-source sample differ')
    if data_lane == 'pp_data' and (sample.get('lane') != data_lane or
            sample.get('period') not in ('0mrad','1p5mrad') or
            sample.get('normalization') != {'kind':'unweighted_DATA_raw_capture','luminosity_applied':False}):
        raise c.ContractError('exact pp period and unweighted DATA normalization required')
    rows = [entry['effective_row'] for entry in selected]
    if data_lane == 'auau_data':
        c.validate_ranges(tuple(c.range_from_manifest(row) for row in rows))
    else:
        _validate_pp_packet_ranges(rows, sample_id)
    row_ids = {row['row_id'] for row in rows}
    if set(augmentation_templates) != row_ids or set(receipt_paths) != row_ids:
        raise c.ContractError('one template and receipt path per exact row required')
    records, plans, offsets = bytearray(), bytearray(), []
    destinations = set()
    for entry, source in zip(selected, rows):
        row = dict(source)
        row['sample_id'] = sample_id
        rid = source['row_id']
        template = dict(augmentation_templates[rid])
        if set(template) not in ({'path','sha256'}, {'path','sha256','locator'}):
            raise c.ContractError('exact augmentation template pin required')
        c._path(template['path'], 'augmentation template')
        c._sha(template['sha256'], 'augmentation template hash')
        if 'locator' in template:
            match = re.fullmatch(r'([0-9]+):([1-9][0-9]*):([0-9a-f]{64})', str(template['locator']))
            if not match or int(match.group(2)) > 65536 or match.group(3) != template['sha256']:
                raise c.ContractError('exact bounded augmentation template locator required')
        row['canonical_augmentation_template'] = template
        output = entry['output']['collaborator_v2_path']
        receipt = receipt_paths[rid]
        for name, path in [('output',output),('receipt',receipt)]:
            c._path(path, name)
            if (not path.startswith('/sphenix/tg/tg01/bulk/') or path in destinations
                    or path.startswith(packet_root+'/')):
                raise c.ContractError('distinct TG publication outside packet required')
            destinations.add(path)
        if not output.endswith('.root') or not receipt.endswith('.json'):
            raise c.ContractError('ROOT output and JSON receipt required')
        if entry['planned_host'] != planned_host:
            raise c.ContractError('DATA lane submit host differs')
        record = ('\t'.join(source['source_paths'][role] for role in ('jets','jetcalo'))+'\n').encode()
        plan = canonical_bytes(row)
        offsets.append((len(records),len(record),hashlib.sha256(record).hexdigest(),
                        len(plans),len(plan),hashlib.sha256(plan).hexdigest()))
        records.extend(record); plans.extend(plan)
    records, plans = bytes(records), bytes(plans)
    record_sha, plan_sha = hashlib.sha256(records).hexdigest(), hashlib.sha256(plans).hexdigest()
    sample.update(tuple_records_path=packet_root+'/records.tsv',tuple_records_sha256=record_sha,
                  tuple_plan_path=packet_root+'/plan.jsonl',tuple_plan_sha256=plan_sha)
    manifest['continuation_source_binding_sha256'] = bound_rows['sha256']
    raw_manifest = canonical_bytes(manifest)
    manifest_sha = hashlib.sha256(raw_manifest).hexdigest()
    arguments = []
    for entry, row, loc in zip(selected, rows, offsets):
        off, length, chunk_sha, poff, plen, psha = loc
        args = [packet_root+'/RUNTIME_MANIFEST.json',manifest_sha,row['row_id'],sample_id,
            sample['tuple_records_path'],record_sha,str(off),str(length),chunk_sha,
            str(row['event_count']),entry['output']['collaborator_v2_path'],receipt_paths[row['row_id']],
            str(request_memory_mb),str(row['source_file_ordinal']-1),sample['tuple_plan_path'],plan_sha,
            str(row['source_file_ordinal']),str(row['event_offset']),str(row['source_global_entry_begin']),
            f'{poff}:{plen}:{psha}']
        if len(args) != 20 or any(any(ch.isspace() for ch in arg) for arg in args):
            raise c.ContractError('worker arguments require exact whitespace-free tokens')
        arguments.append(args)
    return dict(schema='CanonicalDataWorkerPacketV1',status='LOCAL_COMPILED_NOT_ADMITTED',
        execution_authorized=False, rows=len(rows), source_units_per_row=1,
        planned_host=planned_host, worker_arguments=arguments,
        payloads={'records.tsv':records,'plan.jsonl':plans,'RUNTIME_MANIFEST.json':raw_manifest},
        payload_sha256={'records.tsv':record_sha,'plan.jsonl':plan_sha,'RUNTIME_MANIFEST.json':manifest_sha},
        native_output='EXECUTE_NODE_SCRATCH_ONLY', publication='CANONICAL_V2_ROOT_AND_COMMIT_RECEIPT')
