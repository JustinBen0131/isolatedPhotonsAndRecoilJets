"""Local, read-only admission checks for a complete canonical ROOT package.

This module does not confer scientific acceptance, submit jobs or retire inputs.
Presence is intentionally distinct from value validation and dataset coverage.
Unknown source bytes are None, never zero. A converter must account for newly
encountered branches too, not just a historical schema's hard-coded branch list.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
from pathlib import Path
from typing import Iterable, Mapping


PRODUCTS = (
    "pp_data", "pp_photonjet_sim", "pp_inclusivejet_sim", "pp_di_sim",
    "auau_data", "auau_photonjet_sim", "auau_inclusivejet_sim",
)


@dataclass(frozen=True)
class Requirement:
    name: str
    products: tuple[str, ...]
    obligation: str
    route: str


SIM = tuple(p for p in PRODUCTS if p.endswith("sim"))
AUAU = tuple(p for p in PRODUCTS if p.startswith("auau"))
REQUIREMENTS = (
    Requirement("source_coverage", PRODUCTS,
        "Exact accepted source-occurrence and half-open event-range coverage; no duplicates or unexplained pp source omissions", "inventory"),
    Requirement("retained_information", PRODUCTS,
        "Every required source branch/object recoverable inside this package, including shower cells, isolation constituents, all views, UE, models and diagnostics", "conversion"),
    Requirement("jes", PRODUCTS,
        "Reco raw pT before JES and corrected pT after applicable pinned JES; AuAu nominal both after UE; preserve truth and rebuild dependencies", "tree_repair_or_targeted_dst"),
    Requirement("photon_truth_association", SIM,
        "Cluster max-energy primary identity association without photon deltaR veto; embedding-aware, explicit unknown, fake and miss; selection separate", "tree_repair_or_targeted_dst"),
    Requirement("truth_vertex", SIM,
        "Hard and embedded-MB truth vertex, identity and capture validity for the full accepted truth denominator", "coverage_then_targeted_dst"),
    Requirement("centrality", AUAU,
        "Frozen fit-basis triplet hashes and exact-run/explicit legitimate fallback; dependent quantities consistent; updated data-percentile variant separate", "provenance_or_range_replacement"),
    Requirement("cluster_timing", PRODUCTS,
        "Defined cluster time, native units/conversion and validity; do not substitute native MBD event time or retain duplicate scalar aliases", "targeted_dst"),
    Requirement("mbd_charge", PRODUCTS,
        "Calibrated full PMT charge/channel information and finite-positive total with availability; embedding node distinctions retained", "retained_tree_first"),
    Requirement("mbd_event_timing", PRODUCTS,
        "Native MbdOut t0 and arm timing with units and validity, separately from retained per-PMT time", "targeted_dst"),
    Requirement("npb", PRODUCTS,
        "All 25 named ordered inputs, applicable pinned model score/hash/validity and pass decision; existing reference selection unchanged", "exact_feature_recovery_or_targeted_dst"),
    Requirement("event_objects", PRODUCTS,
        "Blair/PPG12 calorimeter sums, max cluster, event decisions and capture states; preserve early-rejected events", "coverage_then_targeted_dst"),
    Requirement("gl1_exposure", ("pp_data", "auau_data"),
        "Native legacy/live/scaled event words with explicit alias semantics; all 64 raw/live/scaled scaler channels, BCO/event/run identity, trigger run map and prescales; retain preselection exposure", "coverage_then_targeted_dst"),
    Requirement("weights_models", PRODUCTS,
        "Frozen score/feature definitions and normalization components with validity, source ownership and exactly-once weighting; no placeholder promotion", "qualification"),
    Requirement("usability_storage", PRODUCTS,
        "ROOT and Uproot examples, bounded reading, exact byte receipts, self-contained provenance and measured lossless migration", "qualification"),
)


def requirement_matrix() -> dict:
    """Every applicable obligation is initially unresolved, not assumed valid."""
    return {p: {r.name: {"status": "UNVERIFIED", "route": r.route,
                       "obligation": r.obligation}
                for r in REQUIREMENTS if p in r.products} for p in PRODUCTS}


def inventory(root) -> dict:
    """Capture all objects and TTree branch schemas without reading baskets."""
    if len(root.keys(recursive=True, cycle=True)) != len(root.keys(recursive=True, cycle=False)):
        raise ValueError("multiple ROOT key cycles need explicit handling")
    result = {}
    for name, cls in root.classnames(recursive=True, cycle=False).items():
        if cls in ("TDirectory", "TDirectoryFile"):
            continue
        if cls == "TTree":
            result[name] = {"class": cls, "rows": int(root[name].num_entries),
                            "branches": dict(root[name].typenames())}
        elif "RNTuple" in cls:
            raise ValueError(f"unsupported non-TTree input: {name}")
        else:
            result[name] = {"class": cls}
    return result


def coverage_report(source_inventory: Mapping, recovery_map: Mapping) -> dict:
    """Check total logical coverage, not whether a claimed mapping is truthful.

    Each source object or branch must name an internal target. Value validation
    remains mandatory; a path to an external native file cannot satisfy this.
    """
    expected = set()
    for name, obj in source_inventory.items():
        if obj["class"] == "TTree":
            expected.update((name, b) for b in obj["branches"])
        else:
            expected.add((name, None))
    missing, external, unverified = [], [], []
    for key in sorted(expected, key=str):
        record = recovery_map.get(key)
        if not record:
            missing.append(key)
            continue
        if record.get("storage") != "internal" or not record.get("target"):
            external.append(key)
        if record.get("readback_status") != "PASS_FULL_VALUES":
            unverified.append(key)
    extra = sorted(set(recovery_map) - expected, key=str)
    return {"status": "PASS_INFORMATION_COVERAGE" if not (missing or external or unverified or extra) else "INCOMPLETE",
            "expected_fields_objects": len(expected), "missing": missing,
            "external_dependencies": external, "unverified_values": unverified,
            "unrecognized_source_keys": extra,
            "scientific_acceptance": False, "retirement_authorized": False}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def publication_record(path, *, product: str, source_ids: Iterable[str],
                       validation_receipt_sha256: str) -> dict:
    """Measure once after close/readback; callers persist this with publication.

    No ROOT scan is needed for later totals. File mutation during hashing fails;
    source IDs name occurrence/range identities, never merely basenames.
    """
    if product not in PRODUCTS:
        raise ValueError("unknown product")
    if len(validation_receipt_sha256) != 64 or any(c not in "0123456789abcdef" for c in validation_receipt_sha256):
        raise ValueError("validated receipt SHA-256 required")
    source_ids = tuple(source_ids)
    if not source_ids or any(not isinstance(s, str) or not s for s in source_ids) or len(source_ids) != len(set(source_ids)):
        raise ValueError("nonempty unique source occurrence/range IDs required")
    path = Path(path).resolve(strict=True)
    before = path.stat()
    digest = sha256_file(path)
    after = path.stat()
    if (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns) != (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns):
        raise ValueError("published file changed during measurement")
    if not path.is_file() or after.st_size <= 0:
        raise ValueError("publication cannot be an empty or non-regular file")
    return {"schema": "CanonicalPartReceiptV1", "product": product,
            "path": str(path), "size_bytes": after.st_size,
            "mtime_ns": after.st_mtime_ns, "sha256": digest,
            "source_ids": list(source_ids),
            "validation_receipt_sha256": validation_receipt_sha256}


def storage_totals(records: Iterable[Mapping]) -> dict:
    """Aggregate frozen metadata only; do not recursively stat scientific data.

    Retry receipt repetitions are idempotent. Conflicting versions of the same
    path or repeated source ownership fail rather than inflate/shrink totals.
    Missing size contributes to an explicit unknown count, not an apparent zero.
    """
    parts, source_owners, path_owners = {}, {}, {}
    for record in records:
        key = (record["product"], record["path"])
        if key[0] not in PRODUCTS:
            raise ValueError("unknown product")
        if not isinstance(key[1], str) or not key[1]:
            raise ValueError("nonempty published path required")
        if key[1] in path_owners and path_owners[key[1]] != key[0]:
            raise ValueError("one physical path cannot belong to two products")
        path_owners[key[1]] = key[0]
        if key in parts:
            if dict(record) != dict(parts[key]):
                raise ValueError("conflicting receipt for one published path")
            continue
        size = record.get("size_bytes")
        if size is not None and (type(size) is not int or size <= 0):
            raise ValueError("unknown size must be null, not zero or a sentinel")
        sources = record.get("source_ids", [])
        if (not isinstance(sources, (tuple, list)) or
                any(not isinstance(s, str) or not s for s in sources) or
                len(sources) != len(set(sources))):
            raise ValueError("invalid or repeated source identity within one part")
        for source in sources:
            owner_key = (key[0], source)
            if owner_key in source_owners:
                raise ValueError("source occurrence/range appears in two parts")
            source_owners[owner_key] = key
        parts[key] = record
    summary = {}
    for product in PRODUCTS:
        selected = [r for (p, _), r in parts.items() if p == product]
        known = [r["size_bytes"] for r in selected if r.get("size_bytes") is not None]
        unknown = len(selected) - len(known)
        summary[product] = {"parts": len(selected), "measured_parts": len(known),
            "unknown_parts": unknown, "known_bytes": sum(known),
            "total_bytes": None if unknown or not selected else sum(known),
            "complete_size_measurement": bool(selected) and not unknown}
    return summary


def storage_breakdown(records: Iterable[Mapping], *, schema10_bytes=None) -> dict:
    """Four dataset totals from receipts, not a final-production forecast.

    A missing product is unknown, not zero. Logical native schema14 information
    inside the canonical package is not charged as a second physical copy.
    schema10 is protected and reported separately, never implicitly retired.
    """
    if schema10_bytes is not None and (type(schema10_bytes) is not int or schema10_bytes < 0):
        raise ValueError("schema10 bytes must be a measured nonnegative integer or null")
    products = storage_totals(records)
    groups = {"pp_data": ("pp_data",),
              "pp_sim": ("pp_photonjet_sim", "pp_inclusivejet_sim", "pp_di_sim"),
              "auau_data": ("auau_data",),
              "auau_sim": ("auau_photonjet_sim", "auau_inclusivejet_sim")}
    breakdown = {}
    for name, members in groups.items():
        complete = all(products[p]["complete_size_measurement"] for p in members)
        known = sum(products[p]["known_bytes"] for p in members)
        breakdown[name] = {"products": list(members), "known_bytes": known,
            "total_bytes": known if complete else None,
            "decimal_TB": known / 10**12 if complete else None,
            "missing_products": [p for p in members if not products[p]["parts"]],
            "unknown_parts": sum(products[p]["unknown_parts"] for p in members),
            "complete_size_measurement": complete}
    complete = all(row["complete_size_measurement"] for row in breakdown.values())
    known = sum(row["known_bytes"] for row in breakdown.values())
    return {"schema": "CanonicalStorageBreakdownV1", "products": products,
        "datasets": breakdown, "known_package_bytes": known,
        "total_package_bytes": known if complete else None,
        "protected_schema10_bytes": schema10_bytes,
        "package_plus_schema10_bytes": known + schema10_bytes if complete and schema10_bytes is not None else None,
        "scope": "SUPPLIED_RECEIPTS_ONLY_NOT_CAMPAIGN_COVERAGE_OR_FORECAST",
        "retirement_authorized": False}
