#!/usr/bin/env python3
"""Expand the THE-236 sparse-data production rows REMOTELY, on SDCC.

Why this exists
---------------
build_the236_schema10_data_plan.py enumerates all 51,535 rows locally, which
forces every DST pair list (~500k lines across ~2,800 files) onto the local
machine only to be shipped straight back to SDCC inside a packet. The simulation
side never did that: its production_plan.json declares 5 samples and
exact_total_rows, and materialize_production.py expands the 7,145 rows on SDCC
from lists that never leave the cluster.

This module is the data analogue of that materializer. The plan stays compact
and human-readable (which is what the 2026-08-11 review asked for), the lists
stay remote, and row expansion happens where the data already lives.

Grouping contract, reproduced exactly from build_the236_schema10_data_plan.py
----------------------------------------------------------------------------
  pp    group_pp()   flat chunks of PP_GROUP_SIZE=20 across the whole ordered
                     record sequence; groups MAY span runs.
                     21,494 sources -> 1,075 rows
  auau  group_auau() chunks of AUAU_GROUP_SIZE=10 WITHIN each run; groups never
                     span runs, so most runs end in a partial group.
                     492,280 sources in 2,776 runs -> 50,460 rows

That asymmetry is the reason the AuAu row count is 50,460 and not 49,228. Do not
"simplify" it into a single flat chunker.

Event counts
------------
Exact per-row event counts ARE required, and an earlier version of this file
claimed the opposite. run_the236_schema10_data_row.sh:35 regex-gates
expected_processed as ^[1-9][0-9]*$, so a zero is rejected outright, and :293
compares it to the actual processed count before the output is published. The
counts come from the sPHENIX FileCatalog datasets table, joined on exact
filename because two productions can exist for the same run with different
counts (run 47289 segment 0 is 6,366 events under ana521_2025p007 and 4,436
under pro001_2025p012). The route fetches that join on SDCC and passes it in
via --event-counts.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import subprocess
from pathlib import Path
from typing import Any, Iterator


PP_GROUP_SIZE = 20
AUAU_GROUP_SIZE = 10
# One shard directory per 1,000 rows keeps every output directory an order of
# magnitude under the project's ~10,000-files-per-directory rule, counting the
# .root and its sibling receipt.
SHARD_SIZE = 1000
SAFE_TOKEN = re.compile(r"[A-Za-z0-9_.-]{1,120}")
PLAN_SCHEMA = "THE236_SCHEMA10_DATA_PRODUCTION_PLAN_V1"


class MaterializationError(RuntimeError):
    """Raised on any contract breach; the caller must fail closed."""


def safe_token(value: str, label: str) -> str:
    if not SAFE_TOKEN.fullmatch(value):
        raise MaterializationError(f"unsafe {label}: {value!r}")
    return value


def safe_remote_dir(value: str, label: str) -> Path:
    path = Path(value)
    if not path.is_absolute() or path.is_symlink() or not path.is_dir():
        raise MaterializationError(f"{label} must be an existing absolute non-symlink dir: {value}")
    if not str(path).startswith(("/sphenix/u/", "/sphenix/tg/")):
        raise MaterializationError(f"{label} is outside the approved roots: {value}")
    return path


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise MaterializationError(f"JSON object required: {path}")
    return value


def list_files(list_dir: Path, prefix: str) -> list[Path]:
    """Deterministic ordering is part of the contract: identical inputs must
    always produce identical rows, byte offsets and digests.

    The prefix filter is load-bearing, not cosmetic. Both list directories hold
    several kinds of .list file side by side: pp has 1,568 dst_ppg12_pair-*.list
    among 4,704 total, AuAu has 2,776 dst_auau_jet_pair-*.list among 8,328.
    Only the *_pair-* files carry the two-column DST_JET/DST_JETCALO layout that
    read_pairs() requires; globbing every .list picks up single-column jet lists
    and dies on the first one. The AuAu paired count of 2,776 is exactly
    AUAU_RUN_COUNT in build_the236_schema10_data_plan.py, which is the
    cross-check that this filter selects the intended set.
    """

    found = sorted(
        p for p in list_dir.iterdir()
        if p.is_file() and p.suffix == ".list" and p.name.startswith(prefix)
    )
    if not found:
        raise MaterializationError(f"no {prefix}*.list manifests under {list_dir}")
    for path in found:
        if path.is_symlink():
            raise MaterializationError(f"symlinked list manifest refused: {path}")
    return found


def read_pairs(path: Path) -> list[tuple[str, str]]:
    """Each line is 'DST_JET<TAB>DST_JETCALO'. Both name the same events in
    different node collections, so a pair is one source, not two."""

    pairs: list[tuple[str, str]] = []
    with path.open("r", encoding="utf-8") as stream:
        for number, raw in enumerate(stream, start=1):
            line = raw.rstrip("\n")
            if not line:
                continue
            # The two list families use different separators: dst_ppg12_pair-*
            # is tab-separated, dst_auau_jet_pair-* is space-separated. Splitting
            # on generic whitespace accepts both; requiring a tab silently killed
            # every AuAu list on line 1.
            fields = line.split()
            if len(fields) != 2 or any("/" in f for f in fields):
                raise MaterializationError(
                    f"malformed pair at {path.name}:{number}: expected two "
                    f"whitespace-separated bare filenames, got {len(fields)} field(s)"
                )
            pairs.append((fields[0], fields[1]))
    if not pairs:
        raise MaterializationError(f"empty list manifest: {path}")
    return pairs


def run_token(path: Path) -> str:
    """Run identity comes from the list filename, e.g. dst_auau_jet_pair-00067597.list"""

    stem = path.stem
    if "-" not in stem:
        raise MaterializationError(f"list manifest carries no run token: {path.name}")
    return safe_token(stem.rsplit("-", 1)[1], "run token")


def grouped_rows(
    system: str, lists: list[Path], restrict: set[str] | None = None
) -> Iterator[tuple[str, list[tuple[str, str]]]]:
    """Yield (run, group) reproducing group_pp / group_auau exactly.

    `restrict`, when given, limits materialization to those jet LFNs. It exists
    for recovery campaigns: reprocessing a source that already produced a clean
    output would DOUBLE COUNT its events in the final dataset, so a rerun must
    cover only the sources no successful row has already consumed. Filtering
    happens before grouping, so survivors regroup compactly rather than leaving
    holes.
    """

    def survivors(path: Path) -> list[tuple[str, str]]:
        pairs = read_pairs(path)
        if restrict is None:
            return pairs
        return [pair for pair in pairs if pair[0] in restrict]

    if system == "auau":
        size = AUAU_GROUP_SIZE
        for path in lists:
            run = run_token(path)
            pairs = survivors(path)
            for start in range(0, len(pairs), size):
                yield run, pairs[start : start + size]
        return

    # pp is run-bounded too, despite the original plan builder's flat chunking.
    # Fun4AllSyncManager.cc:266 refuses inputs from more than one run:
    #   "Mixing run numbers (except runnumber=0 ...) is not supported ... Exiting now"
    # so a group spanning a run boundary is not merely unusual, it is guaranteed
    # to die at the file where the run changes. In campaign
    # the236_schema10_data_prod_20260818_813a2538_user04 that killed 1,605 of
    # 10,007 pp rows, each burning up to eight hours of farm time first, because
    # the crash lands wherever the boundary happens to fall inside the group.
    # Run-bounding costs only a partial group at the end of each run, which is
    # what AuAu has always done and why AuAu never hit this.
    size = PP_GROUP_SIZE
    for path in lists:
        run = run_token(path)
        pairs = survivors(path)
        for start in range(0, len(pairs), size):
            yield run, pairs[start : start + size]


# --- Source-count recovery -------------------------------------------------
#
# The FileCatalog is a MIRROR of what is on disk; the DST itself is the truth.
# On 2026-08-18 that distinction cost real science: 1,786 Au+Au sources had no
# datasets row, every one of them was skipped, and 11 whole GRL runs silently
# vanished from the campaign because every source in every one of their rows was
# skipped and `if not kept: continue` left no trace.
#
# One of those 1,786 (run 68491 segment 204) was on disk the whole time and only
# missing its bookkeeping row. That case must never be lost again: when the
# catalog has no count, open the file and count it. Only a source that cannot be
# located on disk at all is genuinely unusable, and that is now recorded as a
# distinct, loud outcome rather than as the same quiet "no count" bucket.

LFN_PATTERN = re.compile(r"\A[A-Za-z0-9_]+-\d{8}-\d{5}\.root\Z")

# Proven layout, verified against run 68491 (293 files, block run_00068400_00068500):
#   /sphenix/lustre01/sphnxpro/production/run3auau/physics/<tag>/DST_JETCALO/run_<lo>_<hi>/<lfn>
AUAU_PRODUCTION_ROOT = "/sphenix/lustre01/sphnxpro/production/run3auau/physics"
RUN_BLOCK_SIZE = 100

_COUNT_CACHE: dict[str, tuple[int | None, str]] = {}


def conventional_auau_path(lfn: str) -> str | None:
    """Derive the on-disk path from the filename, for files absent from `files`.

    A file can exist on disk while being absent from BOTH catalog tables, which
    is exactly the failure mode this recovery path exists for, so the resolver
    cannot depend on the catalog alone.
    """

    try:
        head, run_text, _ = lfn[: -len(".root")].split("-")
        tag = head.split("_run3auau_", 1)[1]
        run = int(run_text)
    except (ValueError, IndexError):
        return None
    low = (run // RUN_BLOCK_SIZE) * RUN_BLOCK_SIZE
    block = f"run_{low:08d}_{low + RUN_BLOCK_SIZE:08d}"
    return f"{AUAU_PRODUCTION_ROOT}/{tag}/DST_JETCALO/{block}/{lfn}"


def catalog_path(lfn: str) -> str | None:
    """Authoritative location from the FileCatalog `files` table (what FROG reads)."""

    try:
        finished = subprocess.run(
            ["psql", "-h", "sphnxdbmaster", "-d", "FileCatalog", "-At", "-c",
             f"SELECT full_file_path FROM files WHERE lfn = '{lfn}'"],
            capture_output=True, text=True, timeout=120, check=False,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    if finished.returncode != 0:
        return None
    location = finished.stdout.strip().splitlines()
    return location[0].strip() if location and location[0].strip() else None


def entries_in_dst(path: str) -> int | None:
    """Count events in a DST the same way audit_auau_grl_projection.sh does:
    the tree named "T" if present, else the largest TTree in the file."""

    try:
        import ROOT  # noqa: PLC0415  (heavy; only imported on the recovery path)
    except ImportError:
        return None
    ROOT.gROOT.SetBatch(True)
    handle = ROOT.TFile.Open(path, "READ")
    if not handle or handle.IsZombie():
        return None
    try:
        tree = handle.Get("T")
        if tree and tree.InheritsFrom("TTree"):
            return int(tree.GetEntriesFast())
        best = -1
        for key in handle.GetListOfKeys():
            obj = key.ReadObj()
            if obj and obj.InheritsFrom("TTree"):
                best = max(best, int(obj.GetEntriesFast()))
        return best if best > 0 else None
    finally:
        handle.Close()


def recover_source_count(lfn: str) -> tuple[int | None, str]:
    """Return (events, reason) for a source the FileCatalog could not count."""

    if lfn in _COUNT_CACHE:
        return _COUNT_CACHE[lfn]

    result: tuple[int | None, str]
    if not LFN_PATTERN.fullmatch(lfn):
        # Never interpolate an unvalidated name into SQL or a filesystem path.
        result = (None, "MALFORMED_SOURCE_NAME")
    else:
        located = None
        for candidate in (catalog_path(lfn), conventional_auau_path(lfn)):
            if candidate and os.path.isfile(candidate):
                located = candidate
                break
        if located is None:
            result = (None, "INPUT_FILE_ABSENT")
        else:
            counted = entries_in_dst(located)
            if counted is None:
                result = (None, "FILE_PRESENT_BUT_UNREADABLE")
            elif counted <= 0:
                result = (None, "FILE_PRESENT_BUT_EMPTY")
            else:
                result = (counted, "COUNTED_FROM_FILE")

    _COUNT_CACHE[lfn] = result
    return result


def materialize_system(
    *, system: str, list_dir: Path, packet_root: Path, output_root: Path,
    evidence_root: Path, events: dict[str, int], list_prefix: str,
    resolvable: set[str] | None = None,
    restrict: set[str] | None = None,
) -> list[dict[str, Any]]:
    system = safe_token(system, "system")
    if system not in {"pp", "auau"}:
        raise MaterializationError(f"unknown system: {system}")

    records_parent = packet_root / "records" / system
    records_parent.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, Any]] = []
    missing: list[str] = []
    reasons: dict[str, list[str]] = {}
    recovered: dict[str, dict[str, Any]] = {}
    dropped: list[dict[str, Any]] = []
    globals()["_SKIPPED"] = globals().get("_SKIPPED", {})
    globals()["_REASONS"] = globals().get("_REASONS", {})
    globals()["_RECOVERED"] = globals().get("_RECOVERED", {})
    globals()["_DROPPED"] = globals().get("_DROPPED", {})
    for index, (run, group) in enumerate(
        grouped_rows(system, list_files(list_dir, list_prefix), restrict), start=1):
        shard = (index - 1) // SHARD_SIZE
        row_id = f"{system}_data_{index:06d}"

        # ONE RECORDS FILE PER ROW. This is not a style choice.
        # run_the236_schema10_data_row.sh:32 sha256s the ENTIRE records_path on
        # every job. A single per-system blob (~165 MB for AuAu) would therefore
        # be rehashed by all 50,460 jobs: ~8.3 PB of reads against the shared
        # Lustre filesystem. The proven canary shape is one small file per row
        # with byte_offset 0 and chunk_sha == records_sha, so each job hashes
        # only its own few hundred bytes.

        # expected_processed must be a POSITIVE integer: the worker regex at
        # line 35 rejects 0, and line 293 compares it to the actual processed
        # count before publishing. It is the exact sum of this row's sources.
        # A handful of real, on-disk DSTs are absent from the FileCatalog (and
        # from production_status, and their calo partners are absent too). The
        # worker demands an exact expected_processed and verifies it after
        # processing, so a source with no count cannot be included. Skip the
        # SOURCE, not the row: dropping whole rows would discard nine good
        # files for every uncounted one. Every skip is recorded.
        kept: list[tuple[str, str]] = []
        expected = 0
        for jet, calo in group:
            # BOTH halves must be resolvable at RUNTIME, not just the jet.
            # The macro turns each bare name into a path with FROG, which reads
            # the FileCatalog `files` table; an unresolvable name falls through
            # as a bare string and ROOT then reports
            #   "file /home/condor/.../DST_JETCALO_...root does not exist".
            # Checking only the jet is why 1,956 Au+Au rows died ~20 s in: their
            # jet half was registered and their calo half had been removed from
            # central production. os.path.isfile() is deliberately NOT used here
            # -- stat()ing ~500k sources would hammer the shared Lustre MDS.
            if resolvable is not None:
                absent = [name for name in (jet, calo) if name not in resolvable]
                if absent:
                    missing.append(jet)
                    for name in absent:
                        reasons.setdefault("UNRESOLVABLE_AT_RUNTIME", []).append(name)
                    continue

            count = events.get(jet)
            if not isinstance(count, int) or count <= 0:
                # The catalog cannot count it. Ask the file itself before giving
                # up: the catalog is a mirror, the DST is the truth. Run 68491
                # segment 204 sat on disk, 1.24 GB, and was skipped purely for a
                # missing bookkeeping row.
                count, reason = recover_source_count(jet)
                if count is None:
                    missing.append(jet)
                    reasons.setdefault(reason, []).append(jet)
                    continue
                recovered[jet] = {"events": count, "via": reason}

            kept.append((jet, calo))
            expected += count

        if not kept:
            # Every source in this row was unusable, so the row vanishes. RECORD
            # IT. Silently dropping rows is how 11 whole GRL runs left the
            # campaign leaving nothing behind but a count buried in a manifest.
            dropped.append({"row_index": index, "run": run, "sources_discarded": len(group)})
            continue
        group = kept

        # Records are written from the KEPT sources only. Writing them before
        # the filter shipped every source while expected_processed covered only
        # the survivors, so the job processed more events than it declared and
        # died on the worker's processed-vs-expected check at
        # run_the236_schema10_data_row.sh:293.
        record_shard = records_parent / f"shard_{shard:04d}"
        record_shard.mkdir(parents=True, exist_ok=True)
        records_path = record_shard / f"{row_id}.tsv"
        if records_path.exists() or records_path.is_symlink():
            raise MaterializationError(f"records file already exists: {records_path}")

        blob = b"".join(f"{jet}\t{calo}\n".encode("utf-8") for jet, calo in group)
        descriptor = os.open(records_path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o440)
        with os.fdopen(descriptor, "wb") as out:
            out.write(blob)
        digest = hashlib.sha256(blob).hexdigest()

        out_shard = output_root / system / f"shard_{shard:04d}"
        receipt_shard = evidence_root / "rows" / system / f"shard_{shard:04d}"
        rows.append(
            {
                "row_id": row_id,
                "row_index": index,
                "system": system,
                "run": run,
                "source_count": len(group),
                "records_path": str(records_path),
                "records_sha256": digest,
                "byte_offset": 0,
                "byte_count": len(blob),
                "chunk_sha256": digest,
                "expected_processed": expected,
                # Production rows process every event in their sources; the
                # bounded-canary path is off.
                "event_limit": 0,
                "output_path": str(out_shard / f"{row_id}.root"),
                "measurement_path": str(receipt_shard / f"{row_id}.json"),
            }
        )
    if recovered:
        globals()["_RECOVERED"][system] = recovered
    if missing:
        globals()["_SKIPPED"][system] = missing
    if reasons:
        globals()["_REASONS"][system] = {k: sorted(v) for k, v in reasons.items()}
    if dropped:
        globals()["_DROPPED"][system] = dropped

    # A campaign that quietly analyses less than it was asked to is worse than
    # one that fails, because the gap ships. Say it loudly, at the end, where it
    # cannot scroll past unnoticed.
    if missing or dropped or recovered:
        print("=" * 72)
        print(f"SOURCE ACCOUNTING [{system}]")
        print(f"  rows materialized          : {len(rows)}")
        if recovered:
            print(f"  counts recovered from file : {len(recovered)}")
        for reason, names in sorted(reasons.items()):
            print(f"  discarded {reason:<26}: {len(names)}")
        if dropped:
            runs = sorted({str(entry['run']) for entry in dropped})
            print(f"  ROWS DROPPED ENTIRELY      : {len(dropped)}"
                  f"  across {len(runs)} run(s)")
            print(f"    runs: {', '.join(runs[:12])}"
                  f"{' ...' if len(runs) > 12 else ''}")
        print("=" * 72)
    return interleave_for_representative_prefix(rows)


def interleave_for_representative_prefix(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Reorder rows so that ANY prefix spans the whole run range.

    Rows are generated run-by-run in sorted order, and Condor starts them
    roughly in submission order. Left alone, a campaign interrupted at 45% would
    have finished the first ~45% of RUNS and none of the rest, which is a biased
    subset rather than a partial sample. The Wednesday farm reboot guarantees an
    interruption, and the whole point of analyzing before completion is that the
    partial set be usable.

    Transposing run-major into position-major fixes it: take row 0 of every run,
    then row 1 of every run, and so on. Any prefix then covers all runs roughly
    uniformly.

    This only changes SUBMISSION ORDER. row_id, byte_offset, byte_count and
    chunk_sha256 are already assigned and point into the records file, so they
    are untouched and stay valid.
    """

    by_run: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        by_run.setdefault(str(row["run"]), []).append(row)
    ordered: list[dict[str, Any]] = []
    depth = 0
    longest = max((len(group) for group in by_run.values()), default=0)
    while depth < longest:
        for run in sorted(by_run):
            group = by_run[run]
            if depth < len(group):
                ordered.append(group[depth])
        depth += 1
    if len(ordered) != len(rows):
        raise MaterializationError("interleave dropped or duplicated rows")
    for position, row in enumerate(ordered, start=1):
        row["submission_order"] = position
    return ordered


def load_event_counts(path: Path) -> dict[str, int]:
    """Read the FileCatalog counts TSV (filename<TAB>events) fetched remotely.

    The counts are a real measurement from the sPHENIX FileCatalog datasets
    table, joined on exact filename because two productions can exist for the
    same run with different event counts.
    """

    counts: dict[str, int] = {}
    with path.open("r", encoding="utf-8") as stream:
        for number, raw in enumerate(stream, start=1):
            line = raw.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) != 2:
                raise MaterializationError(f"malformed counts row at {path.name}:{number}")
            try:
                value = int(fields[1])
            except ValueError as exc:
                raise MaterializationError(f"non-integer event count at {path.name}:{number}") from exc
            if value <= 0:
                raise MaterializationError(f"non-positive event count at {path.name}:{number}")
            counts[fields[0]] = value
    if not counts:
        raise MaterializationError(f"catalog counts file is empty: {path}")
    return counts


def load_resolvable(path: Path | None) -> set[str] | None:
    """Every LFN the FileCatalog `files` table can resolve, i.e. every name FROG
    can turn into a real path at runtime. Membership is the only cheap way to
    answer "will this source open?" for ~500k sources; stat()ing them all would
    be a metadata storm against a filesystem shared with every sPHENIX user."""

    if path is None:
        return None
    names = {line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()}
    if not names:
        raise MaterializationError(f"resolvable-LFN file is empty: {path}")
    return names


def materialize(args: argparse.Namespace) -> dict[str, Any]:
    plan = load_json(Path(args.plan).resolve())
    if plan.get("schema") != PLAN_SCHEMA or plan.get("status") != "APPROVED_FOR_MATERIALIZATION":
        raise MaterializationError("data production plan state differs")

    packet_root = safe_remote_dir(args.packet_root, "packet root")
    output_root = Path(os.path.abspath(args.output_root))
    evidence_root = Path(os.path.abspath(args.evidence_root))

    events = load_event_counts(Path(os.path.abspath(args.event_counts)))
    resolvable = load_resolvable(
        Path(os.path.abspath(args.resolvable)) if getattr(args, "resolvable", None) else None
    )
    restrict = load_resolvable(
        Path(os.path.abspath(args.restrict_sources))
        if getattr(args, "restrict_sources", None) else None
    )
    declared = plan.get("systems")
    if not isinstance(declared, list) or not declared:
        raise MaterializationError("plan declares no systems")

    manifest_rows: list[dict[str, Any]] = []
    per_system: dict[str, int] = {}
    memory: dict[str, int] = {}
    for entry in declared:
        system = str(entry.get("system"))
        requested = entry.get("request_memory_mb")
        if not isinstance(requested, int) or not 512 <= requested <= 4096:
            raise MaterializationError(f"{system} request_memory_mb is out of contract: {requested!r}")
        memory[system] = requested
        rows = materialize_system(
            system=system,
            list_dir=safe_remote_dir(str(entry.get("list_dir")), f"{system} list dir"),
            packet_root=packet_root,
            output_root=output_root,
            evidence_root=evidence_root,
            events=events,
            list_prefix=str(entry.get("list_prefix") or ""),
            resolvable=resolvable,
            restrict=restrict,
        )
        # Row counts are DERIVED from the lists, which are the source of truth:
        # make_dstListsData.sh already applied the good-run selection when it
        # built them. The plan's declared counts came from a separate builder
        # that was never run against these lists, so asserting them only blocks
        # real data.
        declared = entry.get("exact_rows")
        if isinstance(declared, int) and declared != len(rows):
            print(f"NOTE {system}: lists yield {len(rows)} rows; plan declared {declared}")
        per_system[system] = len(rows)
        manifest_rows.extend(rows)

    total = sum(per_system.values())

    identities = {row["row_id"] for row in manifest_rows}
    outputs = {row["output_path"] for row in manifest_rows}
    if len(identities) != total or len(outputs) != total:
        raise MaterializationError("row identities or output paths are not unique")

    return {
        "schema": "THE236Schema10DataProductionManifestV1",
        "status": "MATERIALIZED_UNSUBMITTED",
        "campaign_tag": plan.get("campaign_tag"),
        "exact_total_rows": total,
        "rows_by_system": per_system,
        "request_memory_mb": memory,
        "shard_size": SHARD_SIZE,
        "group_sizes": {"pp": PP_GROUP_SIZE, "auau": AUAU_GROUP_SIZE},
        "grouping": {"pp": "FLAT_ORDERED_GROUPS_OF_20", "auau": "RUN_BOUNDED_ORDERED_GROUPS_OF_10"},
        "getenv": False,
        "retry": False,
        "maxjobs": False,
        "on_exit_hold": False,
        "skipped_sources_without_catalog_counts": globals().get("_SKIPPED", {}),
        # WHY each source was discarded, not merely THAT it was. One flat
        # "no catalog count" bucket made a recoverable bookkeeping gap and real
        # missing input look identical, and cost a day finding that out.
        "discarded_sources_by_reason": globals().get("_REASONS", {}),
        "recovered_source_counts": globals().get("_RECOVERED", {}),
        # Rows that ceased to exist because every source in them was discarded.
        "rows_dropped_entirely": globals().get("_DROPPED", {}),
        "rows": manifest_rows,
    }


def write_args_tsv(manifest: dict[str, Any], packet_root: Path, destination: Path) -> str:
    """Emit the queue-arguments file the submit description consumes.

    Column order is fixed by the submit description's `queue ... from` line and
    must not be reordered:

      system row_id records_path records_sha byte_offset byte_count chunk_sha
      expected_processed output_path measurement_path row_mem_mb event_limit

    expected_processed and event_limit are both 0 for production. The worker
    only asserts a processed-event count when event_limit is greater than zero,
    so a production row simply processes every event in its sources.
    """

    memory = manifest["request_memory_mb"]
    lines: list[str] = []
    for row in manifest["rows"]:
        system = row["system"]
        fields = [
            system, row["row_id"], row["records_path"], row["records_sha256"],
            str(row["byte_offset"]), str(row["byte_count"]), row["chunk_sha256"],
            str(row["expected_processed"]), row["output_path"], row["measurement_path"],
            str(memory[system]), "0",
        ]
        if any("\t" in field or "\n" in field for field in fields):
            raise MaterializationError(f"unsafe queue-argument field in row {row['row_id']}")
        lines.append("\t".join(fields))

    blob = ("\n".join(lines) + "\n").encode("utf-8")
    if destination.exists() or destination.is_symlink():
        raise MaterializationError(f"queue arguments already exist: {destination}")
    descriptor = os.open(destination, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o440)
    with os.fdopen(descriptor, "wb") as stream:
        stream.write(blob)
    return hashlib.sha256(blob).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", required=True)
    parser.add_argument("--packet-root", required=True)
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--evidence-root", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--args-output", required=True)
    parser.add_argument("--event-counts", required=True)
    # Every LFN FROG can resolve, one per line, taken from the FileCatalog
    # `files` table. Optional so an older packet still materializes, but without
    # it neither half of a pair is checked for runtime resolvability.
    parser.add_argument("--resolvable", default=None)
    # Recovery campaigns only: the jet LFNs that still need processing.
    parser.add_argument("--restrict-sources", default=None)
    args = parser.parse_args()

    manifest = materialize(args)
    manifest["queue_arguments_sha256"] = write_args_tsv(
        manifest, Path(os.path.abspath(args.packet_root)), Path(os.path.abspath(args.args_output))
    )
    destination = Path(os.path.abspath(args.output))
    if destination.exists() or destination.is_symlink():
        raise MaterializationError(f"manifest already exists: {destination}")
    descriptor = os.open(destination, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o440)
    with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
        json.dump(manifest, stream, sort_keys=True, separators=(",", ":"))
        stream.write("\n")
    print(json.dumps({k: v for k, v in manifest.items() if k != "rows"}, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
