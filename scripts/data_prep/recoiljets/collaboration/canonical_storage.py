"""Experimental lossless, streaming ROOT column store with logical TTree views.

Each input branch is mapped explicitly; identical whole columns are stored once.
This is a storage adapter, not a reconstruction or selection stage. Input file
SHA256 and original ROOT typenames are retained, while the reader reconstructs
the original Awkward value types and bits. The physical ROOT layout intentionally
changes (including std::vector to TTree leaf arrays); ROOT file byte identity is
not promised. No full native copy is hidden beside a full convenience export.

Use compact(source, destination) then CanonicalReader(destination).arrays(tree).
All objects in the input are retained or conversion fails. A v2-only input is
explicitly incomplete as a native package. Even passing the minimum native
presence check is not certification of the upstream scientific schema.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import tempfile
import time
from typing import Callable

import awkward as ak
import numpy as np
import uproot

CONTRACT = "CanonicalROOTStoragePrototypeV1"
MANIFEST = "canonical_storage_manifest"
STORE = "canonical_columns"
OBJECTS = "source_objects"
NATIVE_TREES = ("RJShowerCellV1", "RJIsolationConstituentV1", "RJModelEvaluationV1")
UE_BRANCHES = tuple("jet_background_ue_" + name for name in ("emcal", "ihcal", "ohcal"))
ROOT_NATIVE_CONTRACT = "FixedROOTNativeStorageFoundationV1"
ROOT_NATIVE_MANIFEST = "root_native_storage_manifest"
PMT_CANONICAL_PAIRS = (("centrality_pmt_id", "mbd_pmt_id"),
                   ("centrality_pmt_charge", "mbd_pmt_charge"))
# Fixed native_analysis_names_v1 table names, matching canonical_views.TABLES.
# Branch names and physical C++ types remain untouched in this layout.
ROOT_NATIVE_TABLES = {
    "sources": "ReplayFoundationV1/RJSourceOccurrenceV1",
    "events": "ReplayFoundationV1/RJEventV1",
    "photons": "ReplayFoundationV1/RJPhotonCandidateV1",
    "jets": "ReplayFoundationV1/RJJetV1",
    "photonJets": "ReplayFoundationV1/RJPhotonJetPairV1",
    "truthPhotons": "ReplayFoundationV1/RJTruthPhotonV1",
    "truthJets": "ReplayFoundationV1/RJTruthJetV1",
    "recoTruthLinks": "ReplayFoundationV1/RJRecoTruthLinkV1",
    "weights": "ReplayFoundationV1/RJWeightComponentV1",
    "showerCells": "ReplayFoundationV1/RJShowerCellV1",
    "showerViews": "ReplayFoundationV1/RJShowerFeatureViewV1",
    "isolationConstituents": "ReplayFoundationV1/RJIsolationConstituentV1",
    "isolationWitnesses": "ReplayFoundationV1/RJIsolationWitnessV1",
    "models": "ReplayFoundationV1/RJModelEvaluationV1",
    "jetConstituents": "ReplayFoundationV1/RJJetConstituentV1",
    "truthPhotonMisses": "ReplayFoundationV1/RJTruthPhotonMissOccurrenceV1",
    "triggerScalers": "RJTriggerScalersV1",
    "triggerRunInfo": "RJTriggerRunInfoV1",
}


def file_sha256(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def _inventory(root):
    if len(root.keys(recursive=True)) != len(root.keys(recursive=True, cycle=False)):
        raise ValueError("multiple ROOT key cycles need explicit archival handling")
    result = root.classnames(recursive=True, cycle=False)
    if any("RNTuple" in name for name in result.values()):
        raise ValueError("RNTuple input is unsupported; no objects were dropped")
    return result


def _type(array):
    return str(ak.type(array).content)


def _bytes(array):
    data = np.asarray(array)
    return data.astype(data.dtype.newbyteorder("<"), copy=False).tobytes()


def _parts(array):
    """Chunk-independent length and value streams, including NaN payload bits."""
    parts = []

    def walk(layout):
        if isinstance(layout, ak.contents.NumpyArray):
            parts.append(_bytes(layout.data))
        elif isinstance(layout, ak.contents.ListOffsetArray):
            offsets = np.asarray(layout.offsets, dtype="int64")
            parts.append(_bytes(np.diff(offsets)))
            walk(layout.content[offsets[0]:offsets[-1]])
        elif isinstance(layout, ak.contents.RegularArray):
            walk(layout.content[:len(layout) * layout.size])
        else:
            raise ValueError(f"unsupported logical layout: {type(layout).__name__}")

    walk(ak.to_layout(ak.to_packed(array)))
    return parts


def bitwise_equal(a, b):
    return len(a) == len(b) and _type(a) == _type(b) and _parts(a) == _parts(b)


class _ColumnHash:
    def __init__(self, empty):
        self.logical_type = _type(empty)
        self.hashers = [hashlib.sha256() for _ in _parts(empty)]
        self.rows = 0

    def update(self, array):
        parts = _parts(array)
        if _type(array) != self.logical_type or len(parts) != len(self.hashers):
            raise ValueError("column layout changed between chunks")
        for hasher, part in zip(self.hashers, parts):
            hasher.update(part)
        self.rows += len(array)

    def hexdigest(self):
        message = [self.logical_type, self.rows, [h.hexdigest() for h in self.hashers]]
        return hashlib.sha256(json.dumps(message, separators=(",", ":")).encode()).hexdigest()


def _writable_type(empty):
    typ = ak.type(empty).content
    text = str(typ)
    if text == "var * string":
        return "var * uint8", "string_list_bytes_v1"
    node = typ
    if isinstance(node, ak.types.ListType) and text != "string":
        node = node.content
    while isinstance(node, ak.types.RegularType):
        node = node.content
    if text == "string" or isinstance(node, ak.types.NumpyType):
        if isinstance(node, ak.types.NumpyType) and node.primitive not in {
            "bool", "int8", "uint8", "int16", "uint16", "int32", "uint32",
            "int64", "uint64", "float32", "float64",
        }:
            raise ValueError(f"unsupported TTree primitive: {node.primitive}")
        return text, "identity"
    raise ValueError(f"unsupported TTree type: {text}; no branch may be omitted")


def _encode_strings(array):
    """Pack var * string using raw UTF-8 bytes, preserving even embedded NULs."""
    outer = ak.to_layout(ak.to_packed(array))
    inner = outer.content
    offsets, string_offsets = np.asarray(outer.offsets), np.asarray(inner.offsets)
    content = np.asarray(inner.content.data, dtype="uint8")
    payload, row_offsets = bytearray(), [0]
    for start, stop in zip(offsets[:-1], offsets[1:]):
        payload.extend(np.asarray([stop - start], dtype="<u8").tobytes())
        for index in range(start, stop):
            first, last = string_offsets[index:index + 2]
            payload.extend(np.asarray([last - first], dtype="<u8").tobytes())
            payload.extend(content[first:last].tobytes())
        row_offsets.append(len(payload))
    return ak.Array(ak.contents.ListOffsetArray(
        ak.index.Index64(np.asarray(row_offsets, dtype="int64")),
        ak.contents.NumpyArray(np.frombuffer(payload, dtype="uint8"))))


def _decode_strings(array):
    layout = ak.to_layout(ak.to_packed(array))
    offsets, raw = np.asarray(layout.offsets), np.asarray(layout.content.data)
    outer, inner, content = [0], [0], bytearray()
    for start, stop in zip(offsets[:-1], offsets[1:]):
        row, position = raw[start:stop], 0

        def length():
            nonlocal position
            if position + 8 > len(row):
                raise ValueError("truncated string-list encoding")
            value = int(np.frombuffer(row[position:position + 8], dtype="<u8")[0])
            position += 8
            return value

        count = length()
        if count > (len(row) - position) // 8:
            raise ValueError("invalid string-list count")
        for _ in range(count):
            size = length()
            if position + size > len(row):
                raise ValueError("truncated string payload")
            content.extend(row[position:position + size].tobytes())
            position += size
            inner.append(len(content))
        if position != len(row):
            raise ValueError("trailing string-list payload")
        outer.append(len(inner) - 1)
    chars = ak.contents.NumpyArray(np.frombuffer(content, dtype="uint8"),
                                  parameters={"__array__": "char"})
    strings = ak.contents.ListOffsetArray(ak.index.Index64(np.asarray(inner)), chars,
                                         parameters={"__array__": "string"})
    return ak.Array(ak.contents.ListOffsetArray(ak.index.Index64(np.asarray(outer)), strings))


def _encode(array, record):
    if record["encoding"] == "string_list_bytes_v1":
        return _encode_strings(array)
    if record["encoding"] == "float32_exact":
        return ak.values_astype(array, np.float32)
    return array


def _decode(array, record):
    if record["encoding"] == "string_list_bytes_v1":
        return _decode_strings(array)
    if record["encoding"] == "float32_exact":
        return ak.values_astype(array, np.float64)
    if record["encoding"] != "identity":
        raise ValueError("unknown column encoding")
    return array


def _native_coverage(trees):
    basenames, ambiguous = {}, set()
    for name in sorted(trees):
        base = name.rsplit("/", 1)[-1]
        if base in basenames:
            ambiguous.add(base)
        else:
            basenames[base] = name
    for base in ambiguous:
        basenames.pop(base, None)
    missing = [name for name in NATIVE_TREES if name not in basenames]
    missing.extend("ambiguous_native_tree_basename:" + name for name in sorted(ambiguous)
                   if name in (*NATIVE_TREES, "RJEventV1"))
    event = trees.get(basenames.get("RJEventV1"), {})
    missing.extend("RJEventV1/" + branch for branch in UE_BRANCHES
                   if branch not in event.get("branches", {}))
    return {"all_input_objects_retained": True, "unknown_branches_dropped": 0,
            "missing_minimum_native_payload": missing,
            "minimum_native_presence": "FAIL" if missing else "PASS",
            "complete_native_package": "INCOMPLETE" if missing else "NOT_CERTIFIED",
            "note": "Presence only; no upstream schema or scientific certification."}


def _raw_object_hash(root, name):
    chunk, _ = root.key(name).get_uncompressed_chunk_cursor()
    return hashlib.sha256(chunk.raw_data.tobytes()).hexdigest()


def _chunk_entries(tree, step_size, chunk_bytes):
    if tree.num_entries == 0 or chunk_bytes is None:
        return step_size
    return max(1, min(step_size, tree.num_entries_for(chunk_bytes)))


def compact(source, destination, *, step_size=65536, chunk_bytes=16 * 1024**2,
            narrow_float64=True, dedup_scope="tree",
            require_native=False, max_output_bytes=1024**3,
            progress: Callable[[dict], None] | None = None, progress_path=None):
    """Write a new compact file, verify every value, then publish without overwrite.

    Memory scales with one tree chunk, capped by step_size and an uncompressed
    basket-byte estimate. chunk_bytes is a target, not a hard heap bound; a
    single oversized source row is not split. Array caches are disabled.
    Default deduplication stays within each tree for direct ROOT usability.
    cross_tree is experimental and requires the logical reader for joins.
    Progress callbacks run in-process once per chunk.
    Exact narrowing is decided from the whole column, never a sampled prefix.
    Destination must not exist. Temporary failure outputs are removed.
    """
    source, destination = Path(source).resolve(), Path(destination).resolve()
    if step_size < 1 or max_output_bytes < 1:
        raise ValueError("step_size and max_output_bytes must be positive")
    if chunk_bytes is not None and chunk_bytes < 1:
        raise ValueError("chunk_bytes must be positive or None")
    if dedup_scope not in ("tree", "cross_tree"):
        raise ValueError("dedup_scope must be tree or cross_tree")
    if source == destination or destination.exists():
        raise ValueError("destination aliases input or already exists")
    if not destination.parent.is_dir():
        raise ValueError("destination parent must already exist")
    if progress_path is not None:
        progress_path = Path(progress_path).resolve()
        if progress_path in (source, destination) or not progress_path.parent.is_dir():
            raise ValueError("progress path aliases an input/output or lacks its parent")
    source_hash = file_sha256(source)
    manifest = {"contract": CONTRACT, "source_sha256": source_hash,
                "source_bytes": source.stat().st_size, "source_name": source.name,
                "trees": {}, "objects": {}, "columns": {}, "step_size": step_size,
                "chunk_bytes_target": chunk_bytes, "dedup_scope": dedup_scope,
                "narrow_float64": bool(narrow_float64),
                "identity": "source_sha256 + original tree path + original entry index",
                "hash_contract": "logical type, count, and little-endian length/value streams"}
    canonical = {}
    started, completed, total_rows = time.monotonic(), {}, 0

    def tick(phase, name, processed, total):
        if progress or progress_path:
            completed[(phase, name)] = processed
            overall = sum(completed.values())
            elapsed = max(time.monotonic() - started, 1.e-9)
            rate = overall / elapsed
            record = {"phase": phase, "tree": name, "processed": processed,
                      "total": total, "remaining": total - processed,
                      "overall_processed": overall, "overall_total": 3 * total_rows,
                      "overall_remaining": 3 * total_rows - overall,
                      "rate_rows_per_second": rate,
                      "eta_seconds": (3 * total_rows - overall) / rate if rate else None,
                      "source_sha256": source_hash, "source": str(source),
                      "branches": list(manifest["trees"][name]["branches"]),
                      "updated_utc": datetime.now(timezone.utc).isoformat(),
                      "latest_error": None}
            if progress_path:
                fd, temp_name = tempfile.mkstemp(prefix=".progress-", dir=progress_path.parent)
                try:
                    with os.fdopen(fd, "w") as stream:
                        json.dump(record, stream, sort_keys=True)
                    os.replace(temp_name, progress_path)
                finally:
                    Path(temp_name).unlink(missing_ok=True)
            if progress:
                progress(record)

    with uproot.open(source, array_cache=None) as src:
        inventory = _inventory(src)
        manifest["inventory"] = inventory
        total_rows = sum(src[name].num_entries for name, classname in inventory.items()
                         if classname == "TTree")
        for name, classname in inventory.items():
            if classname != "TTree":
                if classname not in ("TDirectory", "TDirectoryFile"):
                    manifest["objects"][name] = {"classname": classname,
                        "sha256": _raw_object_hash(src, name)}
                continue
            tree = src[name]
            if not tree.keys():
                raise ValueError(f"branchless TTree unsupported: {name}")
            info = {"entries": tree.num_entries, "title": tree.title, "branches": {},
                    "chunk_entries": _chunk_entries(tree, step_size, chunk_bytes)}
            manifest["trees"][name] = info
            hashers, eligible = {}, {}
            for branch in tree.keys():
                empty = tree[branch].array(entry_stop=0, library="ak")
                physical_type, encoding = _writable_type(empty)
                hashers[branch] = _ColumnHash(empty)
                eligible[branch] = (narrow_float64 and encoding == "identity" and
                                    "float64" in _type(empty) and tree.num_entries > 0)
                info["branches"][branch] = {"typename": tree[branch].typename,
                    "logical_type": _type(empty), "physical_type": physical_type,
                    "encoding": encoding}
            tick("scan", name, 0, tree.num_entries)
            for start in range(0, tree.num_entries, info["chunk_entries"]):
                stop = min(start + info["chunk_entries"], tree.num_entries)
                arrays = tree.arrays(entry_start=start, entry_stop=stop, library="ak", how=dict)
                for branch, array in arrays.items():
                    hashers[branch].update(array)
                    if eligible[branch]:
                        with np.errstate(over="ignore", invalid="ignore"):
                            back = ak.values_astype(ak.values_astype(array, np.float32), np.float64)
                        eligible[branch] = bitwise_equal(array, back)
                tick("scan", name, stop, tree.num_entries)
            for branch, record in info["branches"].items():
                record["sha256"] = hashers[branch].hexdigest()
                key = (name if dedup_scope == "tree" else None, record["typename"],
                       record["logical_type"], tree.num_entries, record["sha256"])
                if key not in canonical:
                    column = f"c{len(canonical):06d}"
                    canonical[key] = column
                    if eligible[branch]:
                        record["encoding"] = "float32_exact"
                        record["physical_type"] = record["physical_type"].replace("float64", "float32")
                    manifest["columns"][column] = {**record, "tree": name, "branch": branch,
                        "physical_tree": f"{STORE}/{name}", "physical_branch": branch,
                        "entries": tree.num_entries}
                column = canonical[key]
                record.update({key: manifest["columns"][column][key]
                               for key in ("encoding", "physical_type")})
                record["column"] = column
        manifest["coverage"] = _native_coverage(manifest["trees"])
        if require_native and manifest["coverage"]["minimum_native_presence"] != "PASS":
            raise ValueError("minimum native payload absent: " + ", ".join(
                manifest["coverage"]["missing_minimum_native_payload"]))

        fd, temporary_name = tempfile.mkstemp(prefix=".canonical-", suffix=".root", dir=destination.parent)
        os.close(fd)
        temporary = Path(temporary_name)
        try:
            with uproot.recreate(temporary, compression=uproot.ZSTD(5)) as out:
                for name, classname in inventory.items():
                    if classname in ("TDirectory", "TDirectoryFile"):
                        out.mkdir(f"{OBJECTS}/{name}")
                out.copy_from(src, filter_classname=lambda cls: cls not in (
                    "TTree", "TDirectory", "TDirectoryFile"),
                    rename=lambda name: f"{OBJECTS}/{name}", require_matches=False)
                for name, info in manifest["trees"].items():
                    columns = {column: record for column, record in manifest["columns"].items()
                               if record["tree"] == name}
                    if not columns:
                        tick("write", name, info["entries"], info["entries"])
                        continue
                    physical = next(iter(columns.values()))["physical_tree"]
                    schema = {value["physical_branch"]: value["physical_type"]
                              for value in columns.values()}
                    counters = {value["physical_branch"]: "__canonical_count_" + key
                                for key, value in columns.items()}
                    if set(schema) & set(counters.values()):
                        raise ValueError("source branch conflicts with generated count name")
                    writer = out.mktree(physical, schema, title=info["title"],
                                       counter_name=counters.__getitem__)
                    tick("write", name, 0, info["entries"])
                    for start in range(0, info["entries"], info["chunk_entries"]):
                        stop = min(start + info["chunk_entries"], info["entries"])
                        arrays = src[name].arrays([r["branch"] for r in columns.values()],
                            entry_start=start, entry_stop=stop, library="ak", how=dict)
                        writer.extend({record["physical_branch"]: _encode(arrays[record["branch"]], record)
                                       for record in columns.values()})
                        if temporary.stat().st_size > max_output_bytes:
                            raise ValueError("compact output exceeded byte ceiling")
                        tick("write", name, stop, info["entries"])
                manifest["logical_branches"] = sum(len(t["branches"]) for t in manifest["trees"].values())
                manifest["stored_columns"] = len(manifest["columns"])
                manifest["deduplicated_columns"] = manifest["logical_branches"] - manifest["stored_columns"]
                manifest["narrowed_columns"] = sum(r["encoding"] == "float32_exact"
                                                   for r in manifest["columns"].values())
                out[MANIFEST] = json.dumps(manifest, sort_keys=True)
            if temporary.stat().st_size > max_output_bytes:
                raise ValueError("compact output exceeded byte ceiling")
            receipt = validate(source, temporary, step_size=step_size, chunk_bytes=chunk_bytes,
                               progress=lambda item: tick("verify", item["tree"],
                                                          item["processed"], item["total"]))
            if file_sha256(source) != source_hash:
                raise ValueError("source changed during compact conversion")
            # An atomic link cannot overwrite a destination created concurrently.
            os.link(temporary, destination)
            return {**receipt, "output": str(destination),
                    "output_sha256": file_sha256(destination),
                    "source_bytes": source.stat().st_size, "output_bytes": destination.stat().st_size,
                    "dedup_scope": dedup_scope, "step_size": step_size,
                    "chunk_bytes_target": chunk_bytes,
                    "size_ratio": destination.stat().st_size / source.stat().st_size,
                    "size_scope": "this exact input fixture only",
                    **{key: manifest[key] for key in ("logical_branches", "stored_columns",
                                                     "deduplicated_columns", "narrowed_columns")}}
        finally:
            temporary.unlink(missing_ok=True)


class CanonicalReader:
    """Logical tree/branch names remain addressable without materializing a clone."""
    def __init__(self, path):
        self.root = uproot.open(path, array_cache=None)
        try:
            self.manifest = json.loads(str(self.root[MANIFEST]))
            if self.manifest["contract"] != CONTRACT:
                raise ValueError("unknown canonical storage contract")
        except Exception:
            self.root.close()
            raise

    def close(self):
        self.root.close()

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()

    def tree_names(self):
        return list(self.manifest["trees"])

    def typenames(self, tree):
        return {name: record["typename"] for name, record in
                self.manifest["trees"][tree]["branches"].items()}

    def num_entries(self, tree):
        return self.manifest["trees"][tree]["entries"]

    def arrays(self, tree, expressions=None, *, entry_start=0, entry_stop=None):
        info = self.manifest["trees"][tree]
        names = (list(info["branches"]) if expressions is None else
                 [expressions] if isinstance(expressions, str) else list(expressions))
        stop = info["entries"] if entry_stop is None else entry_stop
        if not 0 <= entry_start <= stop <= info["entries"]:
            raise ValueError("entry range is outside logical tree")
        grouped, decoded = {}, {}
        for name in names:
            column = info["branches"][name]["column"]
            record = self.manifest["columns"][column]
            grouped.setdefault(record["physical_tree"], set()).add(column)
        for physical, columns in grouped.items():
            physical_names = {self.manifest["columns"][column]["physical_branch"]: column
                              for column in columns}
            arrays = self.root[physical].arrays(sorted(physical_names), entry_start=entry_start,
                                               entry_stop=stop, library="ak", how=dict)
            decoded.update({physical_names[name]: _decode(array,
                            self.manifest["columns"][physical_names[name]])
                            for name, array in arrays.items()})
        return {name: decoded[info["branches"][name]["column"]] for name in names}

    def iterate(self, tree, expressions=None, *, step_size=8192):
        if step_size < 1:
            raise ValueError("step_size must be positive")
        for start in range(0, self.num_entries(tree), step_size):
            yield self.arrays(tree, expressions, entry_start=start,
                              entry_stop=min(start + step_size, self.num_entries(tree)))

    def object(self, name):
        if name not in self.manifest["objects"]:
            raise KeyError(name)
        return self.root[f"{OBJECTS}/{name}"]


def validate(source, compact_path, *, step_size=65536, chunk_bytes=16 * 1024**2, progress=None):
    """Compare every input row/branch bitwise; recompute chunk-independent hashes."""
    if step_size < 1:
        raise ValueError("step_size must be positive")
    if chunk_bytes is not None and chunk_bytes < 1:
        raise ValueError("chunk_bytes must be positive or None")
    branches = rows = 0
    with uproot.open(source, array_cache=None) as src, CanonicalReader(compact_path) as reader:
        manifest = reader.manifest
        if file_sha256(source) != manifest["source_sha256"] or _inventory(src) != manifest["inventory"]:
            raise ValueError("source identity or object inventory changed")
        expected_trees = {name for name, cls in manifest["inventory"].items() if cls == "TTree"}
        expected_objects = {name for name, cls in manifest["inventory"].items()
                            if cls not in ("TTree", "TDirectory", "TDirectoryFile")}
        if set(manifest["trees"]) != expected_trees or set(manifest["objects"]) != expected_objects:
            raise ValueError("logical inventory omits or adds source objects")
        if manifest["coverage"] != _native_coverage(manifest["trees"]):
            raise ValueError("coverage record differs from logical inventory")
        for column in manifest["columns"].values():
            if reader.root[column["physical_tree"]].num_entries != column["entries"]:
                raise ValueError("stored column row count mismatch")
        for name, cls in manifest["inventory"].items():
            if cls in ("TDirectory", "TDirectoryFile") and f"{OBJECTS}/{name}" not in reader.root:
                raise ValueError(f"source directory missing: {name}")
        for name, info in manifest["trees"].items():
            tree = src[name]
            if (tree.num_entries != info["entries"] or tree.typenames() != reader.typenames(name)
                    or tree.title != info["title"]):
                raise ValueError(f"logical tree schema mismatch: {name}")
            empty = reader.arrays(name, entry_stop=0)
            hashers = {}
            for branch in tree.keys():
                if not bitwise_equal(tree[branch].array(entry_stop=0), empty[branch]):
                    raise ValueError(f"empty/type mismatch: {name}/{branch}")
                hashers[branch] = _ColumnHash(empty[branch])
            chunk_entries = _chunk_entries(tree, step_size, chunk_bytes)
            for start in range(0, tree.num_entries, chunk_entries):
                stop = min(start + chunk_entries, tree.num_entries)
                before = tree.arrays(entry_start=start, entry_stop=stop, library="ak", how=dict)
                after = reader.arrays(name, entry_start=start, entry_stop=stop)
                for branch in before:
                    if not bitwise_equal(before[branch], after[branch]):
                        raise ValueError(f"bitwise mismatch: {name}/{branch} entries {start}:{stop}")
                    hashers[branch].update(after[branch])
                if progress:
                    progress({"phase": "verify", "tree": name, "processed": stop,
                              "total": tree.num_entries, "remaining": tree.num_entries - stop,
                              "source": str(Path(source).resolve()),
                              "source_sha256": manifest["source_sha256"],
                              "branches": list(before),
                              "updated_utc": datetime.now(timezone.utc).isoformat()})
            for branch, hasher in hashers.items():
                if hasher.hexdigest() != info["branches"][branch]["sha256"]:
                    raise ValueError(f"logical hash mismatch: {name}/{branch}")
            if progress and tree.num_entries == 0:
                progress({"phase": "verify", "tree": name, "processed": 0,
                          "total": 0, "remaining": 0})
            rows += tree.num_entries
            branches += len(tree.keys())
        for name, record in manifest["objects"].items():
            expected = _raw_object_hash(src, name)
            if (record["sha256"] != expected or
                    _raw_object_hash(reader.root, f"{OBJECTS}/{name}") != expected):
                raise ValueError(f"object payload mismatch: {name}")
        return {"status": "PASS", "source_sha256": manifest["source_sha256"],
                "trees": len(manifest["trees"]), "branches": branches, "rows": rows,
                "objects": len(manifest["objects"]), "coverage": manifest["coverage"]}


class _NativeProgress:
    def __init__(self, source, total, seconds, callback, path):
        if not np.isfinite(seconds) or seconds <= 0:
            raise ValueError("max_seconds must be finite and positive")
        self.source, self.total = str(source), total
        self.started, self.seconds = time.monotonic(), seconds
        self.callback, self.path, self.completed = callback, path, {}
        self.latest = {}

    def remaining_seconds(self):
        remaining = self.seconds - (time.monotonic() - self.started)
        if remaining <= 0:
            raise TimeoutError("ROOT-native conversion deadline exceeded")
        return remaining

    def emit(self, phase, tree, processed, total, *, status="RUNNING", error=None):
        if status == "RUNNING":
            self.remaining_seconds()
        if phase in ("copy", "verify"):
            self.completed[(phase, tree)] = processed
        done = sum(self.completed.values())
        elapsed = max(time.monotonic() - self.started, 1.e-9)
        rate = done / elapsed
        record = {"contract": ROOT_NATIVE_CONTRACT, "phase": phase, "tree": tree,
            "processed": processed, "total": total, "remaining": total - processed,
            "overall_processed": done, "overall_total": 2 * self.total,
            "overall_remaining": 2 * self.total - done, "source": self.source,
            "status": status, "latest_error": error,
            "rate_rows_per_second": rate,
            "eta_seconds": (2 * self.total - done) / rate if rate else None,
            "updated_utc": datetime.now(timezone.utc).isoformat()}
        self.latest = record
        if self.path is not None:
            fd, name = tempfile.mkstemp(prefix=".native-progress-", dir=self.path.parent)
            try:
                with os.fdopen(fd, "w") as stream:
                    json.dump(record, stream, sort_keys=True)
                os.replace(name, self.path)
            finally:
                Path(name).unlink(missing_ok=True)
        if self.callback:
            self.callback(record)


def _root_native_copier(ROOT):
    # Long copies stay in C++, without Python per-entry serialization. The
    # callback runs in this process at <=4096 rows / ~1 second boundaries.
    if not hasattr(ROOT, "CanonicalNativeStorageCopyV2"):
        if not ROOT.gInterpreter.Declare(r'''
            #include <TTree.h>
            #include <chrono>
            #include <functional>
            #include <stdexcept>
            bool CanonicalNativeStorageCopyV2(TTree *src, TTree *dst,
                double seconds, const std::function<bool(Long64_t)> &progress) {
              using clock = std::chrono::steady_clock;
              const auto started = clock::now();
              auto previous = started;
              const auto total = src->GetEntries();
              for (Long64_t entry = 0; entry < total; ++entry) {
                if (entry % 256 == 0) {
                  const auto now = clock::now();
                  if (std::chrono::duration<double>(now-started).count() >= seconds)
                    throw std::runtime_error("ROOT-native copy deadline exceeded");
                  if (entry % 4096 == 0 ||
                      std::chrono::duration<double>(now-previous).count() >= 1.0) {
                    if (!progress(entry)) return false;
                    previous = now;
                  }
                }
                if (src->GetEntry(entry) < 0)
                  throw std::runtime_error("ROOT-native entry copy failed");
                if (dst->Fill() < 0)
                  throw std::runtime_error("ROOT-native entry fill failed");
              }
              return progress(total);
            }
        '''):
            raise RuntimeError("cannot compile bounded ROOT-native copier")
    return ROOT.CanonicalNativeStorageCopyV2


def _native_schema(tree, omit=()):
    return {"typenames": {k: v for k, v in tree.typenames().items() if k not in omit}, "title": tree.title,
            "branch_titles": {name: tree[name].title for name in tree.keys() if name not in omit},
            "branch_classes": {name: tree[name].classname for name in tree.keys() if name not in omit},
            "aliases": tree.aliases}


def _pmt_projection(tree):
    """Fixed schema: keep the MBD arrays and their existing availability flags."""
    result = {}
    for old, canonical in PMT_CANONICAL_PAIRS:
        if old not in tree.keys():
            continue
        if (canonical not in tree.keys() or 'centrality_pmt_available' not in tree.keys()
                or tree[old].typename != tree[canonical].typename):
            raise ValueError('single PMT array format requires canonical arrays and availability')
        result[old] = canonical
    return result


def _pmt_schema(tree, omit, projection):
    return _native_schema(tree, tuple(omit) + tuple(projection))


def _jet_mass_omissions(name, tree, drop_jet_mass):
    # Exact authorized exception, never a wildcard or content-based suppression.
    if drop_jet_mass and name == "ReplayFoundationV1/RJJetV1" and "mass" in tree.keys():
        if tree.aliases:
            raise ValueError("jet mass omission requires an alias-free native jet tree")
        return ("mass",)
    return ()


def _has_character_leaf(tree):
    return any(leaf.classname == "TLeafC" for name in tree.keys()
               for leaf in tree[name].member("fLeaves"))


def _verify_zstd_baskets(branch, check_deadline=None, expected_tag=b"ZS"):
    """Inspect compression headers, not just misleading cloned fCompress.

    No second basket decompression is needed: read each 9-byte ROOT compression
    header and skip its compressed payload. ROOT may leave tiny baskets raw.
    """
    for index in range(branch.num_baskets):
        if check_deadline:
            check_deadline()
        try:
            key = branch.basket_key(index)
        except ValueError:
            if not branch.basket(index).is_embedded:
                raise
            continue  # embedded baskets are uncompressed
        if key.data_compressed_bytes == key.data_uncompressed_bytes:
            continue
        position = key.data_cursor.index
        end = position + key.data_compressed_bytes
        while position < end:
            if end - position < 9:
                raise ValueError(f"truncated ROOT compression header: {branch.name}")
            header = branch.file.source.chunk(position, position + 9).raw_data.tobytes()
            if header[:2] != expected_tag:
                raise ValueError(f"ROOT-native basket compression policy differs: {branch.name}")
            size = int.from_bytes(header[3:6], "little")
            if size < 1 or position + 9 + size > end:
                raise ValueError(f"invalid ROOT compression block length: {branch.name}")
            position += 9 + size


def _jet_replacement_metadata(replacements, root):
    """Only the coupled reconstructed-JES pT and pair-xJ repair may replace values.

    Numerical/identity provenance belongs to the repair adapter. This storage
    layer independently binds every replacement value, never a manifest-only
    exemption from source readback. Truth, IDs, row populations and raw pT
    cannot be changed through this interface.
    """
    if replacements is None:
        return {}
    expected = {'ReplayFoundationV1/RJJetV1': {'corrected_pt'},
                'ReplayFoundationV1/RJPhotonJetPairV1': {'xjgamma'}}
    if not isinstance(replacements, dict) or set(replacements) != set(expected):
        raise ValueError('JES replacements require both reconstructed jets and pairs')
    result = {}
    for name, fields in expected.items():
        if not isinstance(replacements[name], dict) or set(replacements[name]) != fields:
            raise ValueError('only corrected_pt and xjgamma may be replaced')
        tree = root[name]
        result[name] = {}
        for branch, values in replacements[name].items():
            if (not isinstance(values, np.ndarray) or values.dtype != np.dtype('float64')
                    or values.shape != (tree.num_entries,)
                    or tree[branch].typename != 'double'):
                raise ValueError('replacement dtype/count differs from native scalar branch')
            if branch == 'corrected_pt' and (not np.all(np.isfinite(values)) or np.any(values < 0)):
                raise ValueError('corrected jet pT must be finite and nonnegative')
            if branch == 'xjgamma' and np.any(np.isinf(values)):
                raise ValueError('pair xJ cannot be infinite')
            result[name][branch] = dict(entries=len(values),
                sha256=hashlib.sha256(_bytes(values)).hexdigest(), typename='double')
    return result


def write_root_native(source, destination, *, require_native=False,
                      step_size=65536, chunk_bytes=16 * 1024**2,
                      max_output_bytes=1024**3, max_seconds=300,
                      progress=None, progress_path=None, fast_clone_max_entries=32768,
                      recompress_trees=False, drop_jet_mass=False, jet_replacements=None,
                      recompression_codec="zstd5", compact_pmt=True):
    """Write fixed ROOT-native trees; optionally omit the explicitly excluded jet mass.

    Known trees get fixed friendly table names; unknown trees retain full paths.
    Non-tree objects are copied once under source_objects. This is a storage
    foundation, never a complete-v2 or scientific release claim. No content-
    dependent narrowing, selection, or byte-vector codecs occur. Exactly
    duplicated centrality PMT id/charge vectors are omitted by a fixed schema;
    direct branch access uses the retained mbd_pmt_* names. All other physical
    branches, availability flags and numerical values stay unchanged.
    Standard ROOT TTreeReader/TChain can read native branches directly.
    TChain compatibility proves concatenated iteration only. Native identifiers
    remain scoped to each input part: cross-part joins must include source
    SHA256, native identity domain, and source occurrence where available.
    Never join different parts on their native object/event IDs alone.

    Default: retain input branch compression and fast-clone bounded short trees
    (<=32768 entries and <=32MiB uncompressed). Already-compressed native inputs
    do not justify unconditional decode/recompression. Opt-in recompress_trees
    explicitly rewrites native tree branches using ZSTD level 5, or the
    explicit lzma4 profile. Both retain physical branch types and bitwise values.
    Legacy TLeafC character trees preserve their original baskets because
    ROOT's entry copier does not preserve their physical string encoding.
    Such fallback trees must fit the bounded fast-copy envelope or fail closed.
    Fast-cloning preserves old baskets and does not apply the output file's
    compression setting. Record inherited compression honestly; a ZSTD output
    file setting alone does not prove that copied baskets use ZSTD.
    Rewritten trees use a native C++ loop with progress and deadline checks.
    Runtime checks occur between local I/O calls; they cannot interrupt a
    kernel-stalled file read. The caller's process deadline remains advisable.
    """
    import ROOT
    ROOT.gROOT.SetBatch(True)
    if type(recompress_trees) is not bool:
        raise ValueError("recompress_trees must be a boolean")
    if recompression_codec not in ("zstd5", "lzma4"):
        raise ValueError("unknown lossless recompression codec")
    if recompression_codec != "zstd5" and not recompress_trees:
        raise ValueError("nondefault codec requires explicit tree recompression")
    compression_settings = 204 if recompression_codec == "lzma4" else 505
    compression_mode = "LZMA4_RECOMPRESS" if recompression_codec == "lzma4" else "ZSTD5_RECOMPRESS"
    if type(drop_jet_mass) is not bool:
        raise ValueError("drop_jet_mass must be a boolean")
    if type(compact_pmt) is not bool:
        raise ValueError("compact_pmt must be a boolean")
    source, destination = Path(source).resolve(), Path(destination).resolve()
    if source == destination or destination.exists():
        raise ValueError("destination aliases input or already exists")
    if not destination.parent.is_dir():
        raise ValueError("destination parent must already exist")
    if step_size < 1 or max_output_bytes < 1 or fast_clone_max_entries < 0:
        raise ValueError("invalid ROOT-native resource bound")
    if chunk_bytes is not None and chunk_bytes < 1:
        raise ValueError("chunk_bytes must be positive or None")
    if progress_path is not None:
        progress_path = Path(progress_path).resolve()
        if progress_path in (source, destination) or not progress_path.parent.is_dir():
            raise ValueError("invalid ROOT-native progress path")
    observer = _NativeProgress(source, 0, max_seconds, progress, progress_path)
    temporary = None
    try:
        observer.emit("inventory", "", 0, 0)
        source_hash = file_sha256(source)
        manifest = {"contract": ROOT_NATIVE_CONTRACT, "source_sha256": source_hash,
            "source_bytes": source.stat().st_size, "source_name": source.name,
            "table_mapping_version": "native_analysis_names_v1",
            "trees": {}, "objects": {}, "identity": "source_sha256 + native tree + entry",
            "native_join_scope": "source_sha256 + native identity domain + native IDs; include source occurrence where available",
            "scientific_release": False, "physical_narrowing": False,
            "content_dependent_branch_suppression": False,
            "drop_jet_mass": drop_jet_mass,
            "declared_branch_omissions": {},
            "pmt_storage_policy": 'SINGLE_MBD_ARRAY_V1' if compact_pmt else 'LEGACY_DUPLICATED',
            "pmt_projection": {},
            "tree_compression": ({"mode": compression_mode, "branch_settings": compression_settings}
                if recompress_trees else {"mode": "INHERITED_BASKETS", "file_settings": 505})}
        reverse = {native: friendly for friendly, native in ROOT_NATIVE_TABLES.items()}
        with uproot.open(source, array_cache=None) as src:
            manifest['declared_jet_replacements'] = _jet_replacement_metadata(jet_replacements, src)
            manifest["inventory"] = _inventory(src)
            targets = {ROOT_NATIVE_MANIFEST: "manifest"}
            for name, cls in manifest["inventory"].items():
                observer.remaining_seconds()
                if cls in ("TDirectory", "TDirectoryFile"):
                    continue
                target = reverse.get(name, name) if cls == "TTree" else f"{OBJECTS}/{name}"
                if target == OBJECTS or target.startswith(OBJECTS + "/") and cls == "TTree":
                    raise ValueError(f"tree uses reserved metadata namespace: {name}")
                if target in targets:
                    raise ValueError(f"ROOT-native target collision: {name} and {targets[target]} -> {target}")
                targets[target] = name
                if cls == "TTree":
                    tree = src[name]
                    omitted = _jet_mass_omissions(name, tree, drop_jet_mass)
                    if omitted:
                        manifest["declared_branch_omissions"][name] = list(omitted)
                    projection = _pmt_projection(tree) if compact_pmt else {}
                    if projection:
                        manifest['pmt_projection'][name] = projection
                    schema = _pmt_schema(tree, omitted, projection)
                    manifest["trees"][name] = {**schema, "output_tree": target,
                        "entries": tree.num_entries,
                        "branches": {branch: {"typename": typename}
                                     for branch, typename in schema["typenames"].items()}}
                else:
                    manifest["objects"][name] = {"classname": cls,
                        "output_object": target, "sha256": _raw_object_hash(src, name)}
            for target in targets:
                parts = target.split("/")
                if any("/".join(parts[:i]) in targets for i in range(1, len(parts))):
                    raise ValueError(f"ROOT-native target directory collision: {target}")
            manifest["coverage"] = _native_coverage(manifest["trees"])
            manifest["preserved_character_trees"] = sorted(
                name for name in manifest["trees"] if _has_character_leaf(src[name]))
            if require_native and manifest["coverage"]["minimum_native_presence"] != "PASS":
                raise ValueError("minimum native payload absent")
            schema = {record["output_tree"]: {key: record[key] for key in
                      ("typenames", "branch_titles", "branch_classes", "aliases")}
                      for record in manifest["trees"].values()}
            manifest["physical_schema_sha256"] = hashlib.sha256(json.dumps(
                schema, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
        observer.total = sum(record["entries"] for record in manifest["trees"].values())
        copier = _root_native_copier(ROOT)
        fd, temp_name = tempfile.mkstemp(prefix=".root-native-", suffix=".root", dir=destination.parent)
        os.close(fd)
        temporary = Path(temp_name)
        src = ROOT.TFile.Open(str(source), "READ")
        dst = ROOT.TFile.Open(str(temporary), "RECREATE")
        modes = {}
        try:
            if not src or src.IsZombie() or not dst or dst.IsZombie():
                raise ValueError("cannot open native source/output ROOT file")
            dst.SetCompressionSettings(compression_settings)
            for name, record in manifest["trees"].items():
                observer.emit("copy", name, 0, record["entries"])
                tree = src.Get(name)
                tree.SetBranchStatus("*", 1)
                for branch in manifest["declared_branch_omissions"].get(name, ()):
                    tree.SetBranchStatus(branch, 0)
                for branch in manifest['declared_jet_replacements'].get(name, {}):
                    tree.SetBranchStatus(branch, 0)
                for branch in manifest['pmt_projection'].get(name, {}):
                    tree.SetBranchStatus(branch, 0)
                target = record["output_tree"]
                directory = dst
                for part in target.split("/")[:-1]:
                    directory = directory.GetDirectory(part) or directory.mkdir(part)
                directory.cd()
                character_tree = name in manifest["preserved_character_trees"]
                short_tree = (tree.GetEntries() <= fast_clone_max_entries
                              and tree.GetTotBytes() <= 32 * 1024**2)
                if character_tree and not short_tree:
                    raise ValueError(f"character tree exceeds bounded basket-preserving copy envelope: {name}")
                fast = character_tree or (not recompress_trees and short_tree)
                clone = tree.CloneTree(-1 if fast else 0, "fast" if fast else "")
                if not clone:
                    raise ValueError(f"native tree clone failed: {name}")
                clone.SetName(target.rsplit("/", 1)[-1])
                clone.SetAutoSave(0)
                if recompress_trees and not character_tree:
                    # CloneTree(0) also inherits branch compression settings.
                    # Set every level, including split/nested child branches,
                    # before filling; changing TFile alone is insufficient.
                    def set_compression(branches):
                        for branch in branches:
                            branch.SetCompressionSettings(compression_settings)
                            set_compression(branch.GetListOfBranches())
                    set_compression(clone.GetListOfBranches())
                if not fast:
                    callback_errors = []
                    def copied(processed, native=name, total=record["entries"]):
                        # Never unwind a Python exception through cppyy's C++
                        # std::function trampoline: that can abort ROOT.
                        try:
                            if temporary.stat().st_size > max_output_bytes:
                                raise ValueError("ROOT-native output exceeded byte ceiling")
                            observer.emit("copy", native, int(processed), total)
                            return True
                        except BaseException as error:
                            callback_errors.append(error)
                            return False
                    completed = copier(tree, clone, observer.remaining_seconds(), copied)
                    if callback_errors:
                        raise callback_errors[0]
                    if not completed:
                        raise RuntimeError("ROOT-native copy canceled without a recorded error")
                if clone.GetEntries() != record["entries"]:
                    raise ValueError(f"ROOT-native copied row count differs: {name}")
                # Fill only the new scalar branches, not the cloned tree. The
                # retained population/order and all other native buffers stay
                # untouched. No second full jet/pair table is stored.
                from array import array as scalar_buffer
                for branch in manifest['declared_jet_replacements'].get(name, {}):
                    values = jet_replacements[name][branch]
                    buffer = scalar_buffer('d', [0.])
                    new_branch = clone.Branch(branch, buffer, record['branch_titles'][branch])
                    new_branch.SetCompressionSettings(compression_settings)
                    for index, value in enumerate(values):
                        if index % 1024 == 0:
                            observer.emit('repair', name, index, record['entries'])
                            if temporary.stat().st_size > max_output_bytes:
                                raise ValueError('ROOT-native output exceeded byte ceiling')
                        buffer[0] = float(value)
                        if new_branch.Fill() < 0:
                            raise ValueError('JES scalar branch fill failed')
                    new_branch.ResetAddress()
                if clone.Write(clone.GetName(), ROOT.TObject.kOverwrite) <= 0:
                    raise ValueError(f"ROOT-native tree write failed: {name}")
                modes[name] = ("fast_clone_character_preservation" if character_tree else
                               "fast_clone" if fast else "bounded_cpp_entry_copy")
                observer.emit("copy", name, record["entries"], record["entries"])
                if temporary.stat().st_size > max_output_bytes:
                    raise ValueError("ROOT-native output exceeded byte ceiling")
        finally:
            if dst:
                dst.Close()
            if src:
                src.Close()
        manifest["copy_modes"] = modes
        with uproot.open(source, array_cache=None) as src, uproot.update(temporary) as out:
            for name, cls in manifest["inventory"].items():
                if cls in ("TDirectory", "TDirectoryFile"):
                    out.mkdir(f"{OBJECTS}/{name}")
            out.copy_from(src, filter_classname=lambda cls: cls not in (
                "TTree", "TDirectory", "TDirectoryFile"),
                rename=lambda name: f"{OBJECTS}/{name}", require_matches=False)
            out[ROOT_NATIVE_MANIFEST] = json.dumps(manifest, sort_keys=True)
        if temporary.stat().st_size > max_output_bytes:
            raise ValueError("ROOT-native output exceeded byte ceiling")
        receipt = validate_root_native(source, temporary, step_size=step_size,
            chunk_bytes=chunk_bytes, progress=lambda value: observer.emit(
                "verify", value["tree"], value["processed"], value["total"]),
            check_deadline=observer.remaining_seconds, expected_drop_jet_mass=drop_jet_mass,
            expected_jet_replacements=jet_replacements)
        if file_sha256(source) != source_hash:
            raise ValueError("source changed during ROOT-native conversion")
        observer.remaining_seconds()
        output_hash = file_sha256(temporary)
        observer.remaining_seconds()
        observer.emit("complete", "", 0, 0, status="PASS")
        os.link(temporary, destination)
        return {**receipt, "output": str(destination), "output_sha256": output_hash,
            "source_bytes": source.stat().st_size, "output_bytes": destination.stat().st_size,
            "size_ratio": destination.stat().st_size / source.stat().st_size,
            "elapsed_seconds": time.monotonic() - observer.started,
            "size_scope": "this exact input fixture only", "copy_modes": modes}
    except Exception as error:
        last = observer.latest
        observer.emit("failed", last.get("tree", ""), last.get("processed", 0),
                      last.get("total", 0), status="FAIL", error=str(error))
        raise
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def validate_root_native(source, output, *, step_size=65536, chunk_bytes=16 * 1024**2,
                         progress=None, check_deadline=None,
                         extra_trees=(), extra_objects=None, expected_drop_jet_mass=False,
                         expected_jet_replacements=None):
    """Independent full source-versus-output type, value, order and object audit.

    Staging composition may explicitly allow exact extra TTree paths and an
    exact {metadata_path: ROOT_classname} mapping. Their physics/value validator
    and manifest binding belong to the augmentation owner. This function still
    verifies every original value and records the precise allowed additions;
    it never silently permits other output objects or overrides source trees.
    The caller must explicitly authorize jet-mass omission; all other branches
    remain bitwise checked. A manifest cannot authorize its own exceptions.
    """
    if step_size < 1 or chunk_bytes is not None and chunk_bytes < 1:
        raise ValueError("invalid validation chunk bound")
    if type(expected_drop_jet_mass) is not bool:
        raise ValueError("expected_drop_jet_mass must be a boolean")
    checked_rows = checked_branches = 0
    tree_digests = {}
    extra_trees = tuple(extra_trees)
    extra_objects = dict(extra_objects or {})
    if len(extra_trees) != len(set(extra_trees)) or set(extra_trees) & set(extra_objects):
        raise ValueError("duplicate or ambiguous augmentation paths")
    with uproot.open(source, array_cache=None) as src, uproot.open(output, array_cache=None) as out:
        manifest = json.loads(str(out[ROOT_NATIVE_MANIFEST]))
        if manifest["contract"] != ROOT_NATIVE_CONTRACT:
            raise ValueError("wrong ROOT-native storage contract")
        if manifest.get("drop_jet_mass", False) is not expected_drop_jet_mass:
            raise ValueError("jet mass omission policy differs from caller authorization")
        pmt_policy = manifest.get('pmt_storage_policy', 'LEGACY_DUPLICATED')
        if pmt_policy not in ('SINGLE_MBD_ARRAY_V1', 'LEGACY_DUPLICATED'):
            raise ValueError('unknown PMT storage policy')
        pmt_projection = {name: _pmt_projection(src[name]) for name, cls in _inventory(src).items()
                          if cls == 'TTree'} if pmt_policy == 'SINGLE_MBD_ARRAY_V1' else {}
        pmt_projection = {name: value for name, value in pmt_projection.items() if value}
        if manifest.get('pmt_projection', {}) != pmt_projection:
            raise ValueError('PMT projection differs from fixed single-array schema')
        replacements = _jet_replacement_metadata(expected_jet_replacements, src)
        if manifest.get('declared_jet_replacements', {}) != replacements:
            raise ValueError('JES replacements differ from caller values/authorization')
        compression = manifest.get("tree_compression", {"mode": "LEGACY_UNSPECIFIED"})
        if compression not in ({"mode": "LEGACY_UNSPECIFIED"},
                               {"mode": "ZSTD5_RECOMPRESS", "branch_settings": 505},
                               {"mode": "LZMA4_RECOMPRESS", "branch_settings": 204},
                               {"mode": "INHERITED_BASKETS", "file_settings": 505}):
            raise ValueError("unknown ROOT-native tree compression policy")
        source_inventory, output_inventory = _inventory(src), _inventory(out)
        if source_inventory != manifest["inventory"] or file_sha256(source) != manifest["source_sha256"]:
            raise ValueError("ROOT-native source identity/inventory differs")
        trees = {name for name, cls in source_inventory.items() if cls == "TTree"}
        omissions = {name: list(_jet_mass_omissions(name, src[name], expected_drop_jet_mass))
                     for name in trees}
        omissions = {name: branches for name, branches in omissions.items() if branches}
        if manifest.get("declared_branch_omissions", {}) != omissions:
            raise ValueError("declared branch omissions differ from exact jet mass policy")
        character_trees = sorted(name for name in trees if _has_character_leaf(src[name]))
        if "preserved_character_trees" in manifest and manifest["preserved_character_trees"] != character_trees:
            raise ValueError("ROOT-native character preservation policy differs from source")
        objects = {name for name, cls in source_inventory.items()
                   if cls not in ("TTree", "TDirectory", "TDirectoryFile")}
        if trees != set(manifest["trees"]) or objects != set(manifest["objects"]):
            raise ValueError("ROOT-native manifest omits or adds source objects")
        reverse = {native: friendly for friendly, native in ROOT_NATIVE_TABLES.items()}
        observed_schema, observed_coverage = {}, {}
        for name, record in manifest["trees"].items():
            if record["output_tree"] != reverse.get(name, name):
                raise ValueError("ROOT-native fixed table mapping differs")
            schema = _pmt_schema(src[name], omissions.get(name, ()), pmt_projection.get(name, {}))
            observed_schema[record["output_tree"]] = {key: schema[key] for key in (
                "typenames", "branch_titles", "branch_classes", "aliases")}
            observed_coverage[name] = {"branches": schema["typenames"]}
        fingerprint = hashlib.sha256(json.dumps(observed_schema, sort_keys=True,
                                    separators=(",", ":")).encode()).hexdigest()
        if fingerprint != manifest["physical_schema_sha256"]:
            raise ValueError("ROOT-native physical schema fingerprint differs from source")
        if any(record["output_object"] != f"{OBJECTS}/{name}" or
               record["classname"] != source_inventory[name]
               for name, record in manifest["objects"].items()):
            raise ValueError("ROOT-native fixed metadata mapping differs")
        expected = {ROOT_NATIVE_MANIFEST: "TObjString"}
        expected.update({record["output_tree"]: "TTree" for record in manifest["trees"].values()})
        expected.update({record["output_object"]: record["classname"]
                         for record in manifest["objects"].values()})
        if set(expected) & (set(extra_trees) | set(extra_objects)):
            raise ValueError("augmentation cannot replace original objects")
        expected.update({name: "TTree" for name in extra_trees})
        expected.update(extra_objects)
        actual = {name: cls for name, cls in output_inventory.items()
                  if cls not in ("TDirectory", "TDirectoryFile")}
        if actual != expected:
            raise ValueError("ROOT-native physical inventory differs or contains extra stores")
        for name, cls in source_inventory.items():
            if cls in ("TDirectory", "TDirectoryFile") and f"{OBJECTS}/{name}" not in out:
                raise ValueError(f"ROOT-native source directory missing: {name}")
        for name, record in manifest["trees"].items():
            if check_deadline:
                check_deadline()
            a, b = src[name], out[record["output_tree"]]
            projection = pmt_projection.get(name, {})
            schema = _pmt_schema(a, omissions.get(name, ()), projection)
            if (_native_schema(b) != schema or a.num_entries != b.num_entries or
                    a.num_entries != record["entries"] or
                    any(record[key] != value for key, value in schema.items())):
                raise ValueError(f"ROOT-native physical schema/count differs: {name}")
            hashers = {}
            original_schema = _native_schema(a, omissions.get(name, ()))
            for branch in original_schema["typenames"]:
                physical_branch = projection.get(branch, branch)
                if (compression["mode"] in ("ZSTD5_RECOMPRESS", "LZMA4_RECOMPRESS") and name not in character_trees
                        and b[physical_branch].member("fCompress") != compression["branch_settings"]):
                    raise ValueError(f"ROOT-native branch compression policy differs: {name}/{branch}")
                if compression["mode"] in ("ZSTD5_RECOMPRESS", "LZMA4_RECOMPRESS") and name not in character_trees:
                    _verify_zstd_baskets(b[physical_branch], check_deadline,
                        b"XZ" if compression["mode"] == "LZMA4_RECOMPRESS" else b"ZS")
                empty = a[branch].array(entry_stop=0)
                if not bitwise_equal(empty, b[physical_branch].array(entry_stop=0)):
                    raise ValueError(f"ROOT-native empty/type mismatch: {name}/{branch}")
                hashers[branch] = _ColumnHash(empty)
            chunk = _chunk_entries(a, step_size, chunk_bytes)
            for start in range(0, a.num_entries, chunk):
                if check_deadline:
                    check_deadline()
                stop = min(start + chunk, a.num_entries)
                before = a.arrays(list(original_schema["typenames"]), entry_start=start, entry_stop=stop, library="ak", how=dict)
                after = b.arrays(entry_start=start, entry_stop=stop, library="ak", how=dict)
                if projection:
                    from centrality_replay import restore_pmt_columns
                    after = restore_pmt_columns(after)
                for branch, array in before.items():
                    expected_array = (ak.Array(expected_jet_replacements[name][branch][start:stop])
                        if branch in replacements.get(name, {}) else array)
                    if not bitwise_equal(expected_array, after[branch]):
                        raise ValueError(f"ROOT-native bitwise mismatch: {name}/{branch} {start}:{stop}")
                    hashers[branch].update(expected_array)
                if progress:
                    progress({"tree": name, "processed": stop, "total": a.num_entries})
            if progress and a.num_entries == 0:
                progress({"tree": name, "processed": 0, "total": 0})
            hashes = {branch: value.hexdigest() for branch, value in hashers.items()}
            tree_digests[name] = hashlib.sha256(json.dumps(hashes, sort_keys=True).encode()).hexdigest()
            checked_rows += a.num_entries
            checked_branches += len(schema["typenames"])
        for name, record in manifest["objects"].items():
            if check_deadline:
                check_deadline()
            expected_hash = _raw_object_hash(src, name)
            if record["sha256"] != expected_hash or _raw_object_hash(out, record["output_object"]) != expected_hash:
                raise ValueError(f"ROOT-native object payload mismatch: {name}")
        coverage = _native_coverage(observed_coverage)
        if coverage != manifest["coverage"]:
            raise ValueError("ROOT-native coverage record differs")
        return {"status": "PASS", "contract": ROOT_NATIVE_CONTRACT,
            "declared_branch_omissions": omissions,
            "declared_jet_replacements": replacements,
            "all_retained_values_unchanged": not bool(replacements),
            "all_input_branches_retained": not bool(omissions or pmt_projection),
            "all_input_values_retained": not bool(omissions or replacements),
            "pmt_projection": pmt_projection,
            "tree_compression": compression,
            "preserved_character_trees": character_trees,
            "source_sha256": manifest["source_sha256"],
            "physical_schema_sha256": manifest["physical_schema_sha256"],
            "logical_content_sha256": hashlib.sha256(json.dumps(tree_digests, sort_keys=True).encode()).hexdigest(),
            "trees": len(trees), "branches": checked_branches, "rows": checked_rows,
            "objects": len(objects), "coverage": coverage, "scientific_release": False,
            "authorized_extra_trees": sorted(extra_trees),
            "authorized_extra_objects": extra_objects,
            "augmentation_values_verified": False if extra_trees or extra_objects else None}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("destination", type=Path)
    parser.add_argument("--step-size", type=int, default=65536)
    parser.add_argument("--chunk-bytes", type=int, default=16 * 1024**2)
    parser.add_argument("--dedup-scope", choices=("tree", "cross_tree"), default="tree")
    parser.add_argument("--no-narrow", action="store_true")
    parser.add_argument("--require-native", action="store_true")
    parser.add_argument("--progress-json", type=Path)
    parser.add_argument("--layout", choices=("experimental", "root-native"), default="experimental")
    parser.add_argument("--max-seconds", type=float, default=300)
    args = parser.parse_args()
    common = dict(step_size=args.step_size, chunk_bytes=args.chunk_bytes,
                  require_native=args.require_native, progress_path=args.progress_json,
                  progress=lambda value: print(json.dumps(value), flush=True))
    if args.layout == "root-native":
        receipt = write_root_native(args.source, args.destination,
                                    max_seconds=args.max_seconds, **common)
    else:
        receipt = compact(args.source, args.destination, dedup_scope=args.dedup_scope,
                          narrow_float64=not args.no_narrow, **common)
    print(json.dumps(receipt, sort_keys=True))


if __name__ == "__main__":
    main()
