#!/usr/bin/env python3
"""Reduce parallel schema-10 ROOT histograms without reading any TTree.

This reducer intentionally touches only the producer-owned direct histogram
objects and ``analysis_config_yaml``.  It never opens, iterates, or projects a
schema-10 table.  Missing per-file histograms are legal because RecoilJets
books several objects lazily; the assembled campaign must establish that each
required observable exists at least once per system.
"""

from __future__ import annotations

import argparse
import base64
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import re
import sys
from typing import Any, Iterable


SCHEMA = "THE243Schema10DirectHistogramReductionV1"
PRODUCT_ROOT = Path(
    "/sphenix/tg/tg01/bulk/jbennett/thesisAna/recoiljets/outputs/"
    "the236_schema10_data_prod_20260818_813a2538_user04"
)
PP_ACCEPTED_PRODUCT_ROOTS = (
    PRODUCT_ROOT / "pp",
    Path(
        "/sphenix/tg/tg01/bulk/jbennett/thesisAna/recoiljets/outputs/"
        "the236_schema10_analyzable_recovery_20260820_v2/pp"
    ),
    Path(
        "/sphenix/tg/tg01/bulk/jbennett/thesisAna/recoiljets/outputs/"
        "the236_schema10_analyzable_recovery_exact611_s1_20260821/pp"
    ),
)

LEADING_ABCD_BASES = {
    "a": "h_xJpurityLead_isIsolated_isTight",
    "b": "h_xJpurityLead_notIsolated_isTight",
    "c": "h_xJpurityLead_isIsolated_notTight",
    "d": "h_xJpurityLead_notIsolated_notTight",
}
LEADING_ABCD_PT_BINS = ((15, 20), (20, 25), (25, 35))


def histogram_contract(system: str) -> tuple[str, dict[str, str]]:
    if system == "pp":
        top = "Photon_4_GeV_plus_MBD_NS_geq_1"
        suffix = ""
    elif system == "auau":
        # The producer books one observable family per configured scaled
        # trigger.  MBD-only is QA, not the nominal Photon10 population.
        top = "photon_10_plus_MBD_NS_geq_2_vtx_lt_150"
        suffix = "_cent_0_20"
    else:  # pragma: no cover - argparse prevents this path.
        raise ValueError(f"unsupported system: {system}")
    objects = {
        "recoil_xj_region_a": f"h2_unfoldReco_pTgamma_xJ_incl_r04{suffix}",
        "recoil_xj_region_c": f"h2_unfoldReco_pTgamma_xJ_incl_sidebandC_r04{suffix}",
        "recoil_xj_leading_check": f"h2_unfoldReco_pTgamma_xJ_lead_sam_compat_r04{suffix}",
        "recoil_dphi_region_a": f"h2_unfoldReco_pTgamma_dphi_incl_r04{suffix}",
        "photon_denominator_region_a": f"h_unfoldRecoPho_pTgamma{suffix}",
        "abcd_region_a": f"h_pTgamma_ABCD_A{suffix}",
        "abcd_region_b": f"h_pTgamma_ABCD_B{suffix}",
        "abcd_region_c": f"h_pTgamma_ABCD_C{suffix}",
        "abcd_region_d": f"h_pTgamma_ABCD_D{suffix}",
        "vertex_z": f"h_vertexZ{suffix}",
    }
    # These one-bin producer histograms are filled once per event using the
    # independently leading photon in each ABCD region.  They are the exact
    # direct-histogram counterpart of the compact reader's event-leading ABCD
    # projection; h_pTgamma_ABCD_* above is candidate-level and must not be
    # substituted for this parity check.
    for label, base in LEADING_ABCD_BASES.items():
        for low, high in LEADING_ABCD_PT_BINS:
            objects[f"abcd_leading_{label}_pt{low}_{high}"] = (
                f"{base}_pT_{low}_{high}{suffix}"
            )
    if system == "auau":
        objects["centrality"] = "h_centrality"
    return top, objects


def _axis_edges(axis: Any) -> list[float]:
    bins = int(axis.GetNbins())
    if bins <= 0:
        raise ValueError("positive histogram axis cardinality required")
    edges = [float(axis.GetBinLowEdge(index)) for index in range(1, bins + 1)]
    edges.append(float(axis.GetBinUpEdge(bins)))
    if not all(math.isfinite(value) for value in edges):
        raise ValueError("finite histogram axis edges required")
    if any(right <= left for left, right in zip(edges, edges[1:])):
        raise ValueError("strictly increasing histogram axis edges required")
    return edges


def _zeros(shape: list[int]) -> Any:
    if len(shape) == 1:
        return [0.0] * shape[0]
    return [_zeros(shape[1:]) for _ in range(shape[0])]


def _add_nested(target: Any, source: Any) -> None:
    if not isinstance(target, list) or not isinstance(source, list) or len(target) != len(source):
        raise ValueError("histogram payload shape differs")
    for index, value in enumerate(source):
        if isinstance(value, list):
            _add_nested(target[index], value)
        else:
            target[index] += float(value)


def extract_histogram(histogram: Any) -> dict[str, Any]:
    dimension = int(histogram.GetDimension())
    if dimension not in (1, 2):
        raise ValueError(f"only TH1/TH2 direct histograms are supported, got dimension {dimension}")
    axes = [_axis_edges(histogram.GetXaxis())]
    if dimension == 2:
        axes.append(_axis_edges(histogram.GetYaxis()))
    shape_with_flow = [len(edges) + 1 for edges in axes]
    sumw = _zeros(shape_with_flow)
    sumw2 = _zeros(shape_with_flow)
    if dimension == 1:
        for xbin in range(shape_with_flow[0]):
            value = float(histogram.GetBinContent(xbin))
            error = float(histogram.GetBinError(xbin))
            sumw[xbin] = value
            sumw2[xbin] = error * error
    else:
        for xbin in range(shape_with_flow[0]):
            for ybin in range(shape_with_flow[1]):
                value = float(histogram.GetBinContent(xbin, ybin))
                error = float(histogram.GetBinError(xbin, ybin))
                sumw[xbin][ybin] = value
                sumw2[xbin][ybin] = error * error
    return {
        "class_name": str(histogram.ClassName()),
        "dimension": dimension,
        "axis_edges": axes,
        "shape_with_flow": shape_with_flow,
        "sumw_with_flow": sumw,
        "sumw2_with_flow": sumw2,
        "entries": float(histogram.GetEntries()),
    }


def _same_edges(left: list[list[float]], right: list[list[float]]) -> bool:
    if len(left) != len(right):
        return False
    return all(
        len(a) == len(b)
        and all(math.isclose(x, y, rel_tol=0.0, abs_tol=1e-12) for x, y in zip(a, b))
        for a, b in zip(left, right)
    )


def merge_histogram(accumulator: dict[str, Any] | None, row: dict[str, Any]) -> dict[str, Any]:
    if accumulator is None:
        return {
            **row,
            "contributing_files": 1,
            "missing_files": 0,
        }
    for key in ("dimension", "shape_with_flow"):
        if accumulator[key] != row[key]:
            raise ValueError(f"direct histogram {key} differs")
    if not _same_edges(accumulator["axis_edges"], row["axis_edges"]):
        raise ValueError("direct histogram axis edges differ")
    _add_nested(accumulator["sumw_with_flow"], row["sumw_with_flow"])
    _add_nested(accumulator["sumw2_with_flow"], row["sumw2_with_flow"])
    accumulator["entries"] += row["entries"]
    accumulator["contributing_files"] += 1
    return accumulator


def _config_text(value: Any) -> str:
    if value is None:
        return ""
    if hasattr(value, "GetString"):
        string = value.GetString()
        if hasattr(string, "Data"):
            return str(string.Data())
        return str(string)
    if hasattr(value, "GetTitle"):
        return str(value.GetTitle())
    return str(value)


def _missing_root_object(value: Any) -> bool:
    """Recognize both Python ``None`` and PyROOT null-pointer proxies."""
    if value is None:
        return True
    try:
        return not bool(value)
    except (TypeError, ValueError):
        return False


def _approved_product_roots(system: str) -> tuple[Path, ...]:
    if system == "pp":
        return PP_ACCEPTED_PRODUCT_ROOTS
    if system == "auau":
        return (PRODUCT_ROOT / "auau",)
    raise ValueError(f"unsupported system: {system}")


def _require_approved_product_path(candidate: Path, system: str) -> None:
    for root in _approved_product_roots(system):
        try:
            candidate.relative_to(root)
            return
        except ValueError:
            continue
    raise ValueError(f"input is outside the frozen {system} product roots: {candidate}")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _normalize_manifest(payload: Any, system: str) -> list[dict[str, Any]]:
    rows = payload.get("inputs") if isinstance(payload, dict) else payload
    if not isinstance(rows, list) or not rows:
        raise ValueError("nonempty direct-histogram input manifest required")
    normalized: list[dict[str, Any]] = []
    for row in rows:
        if not isinstance(row, dict):
            raise ValueError("direct-histogram input rows must be objects")
        candidate = Path(str(row.get("path", row.get("root_path", ""))))
        size = int(row.get("size_bytes", row.get("root_size_bytes", -1)))
        digest = str(row.get("sha256", row.get("root_sha256", "")))
        if not re.fullmatch(r"[0-9a-f]{64}", digest) or size <= 0:
            raise ValueError("sealed path, SHA-256, and positive size required")
        _require_approved_product_path(candidate, system)
        normalized.append(
            {
                "path": str(candidate),
                "sha256": digest,
                "size_bytes": size,
                "cohort": row.get("cohort"),
            }
        )
    if len({row["path"] for row in normalized}) != len(normalized):
        raise ValueError("duplicate direct-histogram input path")
    return normalized


def _load_manifest(path: Path, system: str) -> list[dict[str, Any]]:
    return _normalize_manifest(json.loads(path.read_text(encoding="utf-8")), system)


def _load_manifest_b64_gzip(value: str, system: str) -> list[dict[str, Any]]:
    raw = gzip.decompress(base64.b64decode(value, validate=True))
    return _normalize_manifest(json.loads(raw), system)


def _direct_inputs(
    paths: Iterable[str],
    system: str,
    sha256s: Iterable[str] | None = None,
    sealed_sizes: Iterable[int] | None = None,
) -> list[dict[str, Any]]:
    path_rows = list(paths)
    digest_rows = list(sha256s or [])
    size_rows = list(sealed_sizes or [])
    if digest_rows and len(digest_rows) != len(path_rows):
        raise ValueError("one direct-input SHA-256 is required per input path")
    if size_rows and len(size_rows) != len(path_rows):
        raise ValueError("one direct-input sealed size is required per input path")
    rows = []
    for index, raw in enumerate(path_rows):
        path = Path(raw)
        _require_approved_product_path(path, system)
        size = path.stat().st_size
        sealed_size = size_rows[index] if size_rows else size
        if sealed_size != size:
            raise ValueError(f"direct-input sealed size changed: {path}")
        digest = digest_rows[index] if digest_rows else None
        if digest is not None and not re.fullmatch(r"[0-9a-f]{64}", digest):
            raise ValueError("direct-input SHA-256 must be lowercase hexadecimal")
        rows.append({"path": str(path), "sha256": digest, "size_bytes": sealed_size})
    return rows


def reduce_inputs(system: str, inputs: list[dict[str, Any]], root_module: Any) -> dict[str, Any]:
    topdir, contract = histogram_contract(system)
    accumulators: dict[str, dict[str, Any] | None] = {name: None for name in contract}
    missing: dict[str, int] = {name: 0 for name in contract}
    config_hash_counts: dict[str, int] = {}
    source_rows = []
    for source in inputs:
        path = Path(source["path"])
        stat_size = path.stat().st_size
        if stat_size != int(source["size_bytes"]):
            raise ValueError(f"source size changed: {path}")
        expected_sha256 = source.get("sha256")
        if expected_sha256 and _sha256(path) != expected_sha256:
            raise ValueError(f"source SHA-256 changed: {path}")
        root_file = root_module.TFile.Open(str(path), "READ")
        if root_file is None or not root_file.IsOpen() or root_file.IsZombie():
            raise ValueError(f"unreadable ROOT input: {path}")
        try:
            if bool(root_file.TestBit(root_module.TFile.kRecovered)):
                raise ValueError(f"recovered ROOT input is not accepted: {path}")
            config = _config_text(root_file.Get("analysis_config_yaml"))
            if not config:
                raise ValueError(f"analysis_config_yaml is absent: {path}")
            config_hash = hashlib.sha256(config.encode("utf-8")).hexdigest()
            config_hash_counts[config_hash] = config_hash_counts.get(config_hash, 0) + 1
            for logical_name, object_name in contract.items():
                obj = root_file.Get(f"{topdir}/{object_name}")
                if _missing_root_object(obj):
                    missing[logical_name] += 1
                    continue
                row = extract_histogram(obj)
                row["root_object"] = f"{topdir}/{object_name}"
                accumulators[logical_name] = merge_histogram(accumulators[logical_name], row)
            source_rows.append(
                {
                    "path": str(path),
                    "sha256": source.get("sha256"),
                    "size_bytes": stat_size,
                    "analysis_config_sha256": config_hash,
                    "cohort": source.get("cohort"),
                }
            )
        finally:
            root_file.Close()
    histograms: dict[str, Any] = {}
    for logical_name, accumulator in accumulators.items():
        if accumulator is None:
            histograms[logical_name] = {
                "root_object": f"{topdir}/{contract[logical_name]}",
                "contributing_files": 0,
                "missing_files": len(inputs),
                "status": "ABSENT_IN_THIS_REDUCTION_PARTITION",
            }
        else:
            accumulator["missing_files"] = missing[logical_name]
            accumulator["status"] = "PRESENT"
            histograms[logical_name] = accumulator
    return {
        "schema": SCHEMA,
        "status": "PASS",
        "system": system,
        "producer_histogram_contract": {
            "schema": "THE243Schema10ParallelDirectHistogramContractV1",
            "top_directory": topdir,
            "objects": contract,
            "ttree_reads": 0,
            "dst_reads": 0,
        },
        "input_count": len(inputs),
        "unique_input_count": len({row["path"] for row in source_rows}),
        "source_sha256_verification_count": sum(
            1 for row in source_rows if row.get("sha256")
        ),
        "input_size_bytes": sum(row["size_bytes"] for row in source_rows),
        "analysis_config_sha256_counts": dict(sorted(config_hash_counts.items())),
        "inputs": source_rows,
        "histograms": histograms,
    }


def write_json(path: str, payload: dict[str, Any]) -> None:
    text = json.dumps(payload, sort_keys=True, separators=(",", ":")) + "\n"
    if path == "-":
        sys.stdout.write(text)
        return
    output = Path(path)
    descriptor = os.open(output, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o640)
    with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
        stream.write(text)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--system", choices=("pp", "auau"), required=True)
    parser.add_argument("--input", action="append", default=[])
    parser.add_argument("--input-sha256", action="append", default=[])
    parser.add_argument("--input-size-bytes", action="append", type=int, default=[])
    parser.add_argument("--input-manifest", type=Path)
    parser.add_argument("--input-manifest-b64-gzip")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    input_modes = sum(
        (bool(args.input), bool(args.input_manifest), bool(args.input_manifest_b64_gzip))
    )
    if input_modes != 1:
        raise ValueError(
            "choose exactly one of --input, --input-manifest, or --input-manifest-b64-gzip"
        )
    if args.input_manifest:
        if args.input_sha256 or args.input_size_bytes:
            raise ValueError("sealed direct-input options cannot accompany --input-manifest")
        inputs = _load_manifest(args.input_manifest, args.system)
    elif args.input_manifest_b64_gzip:
        if args.input_sha256 or args.input_size_bytes:
            raise ValueError("sealed direct-input options cannot accompany encoded manifest")
        inputs = _load_manifest_b64_gzip(args.input_manifest_b64_gzip, args.system)
    else:
        inputs = _direct_inputs(
            args.input,
            args.system,
            sha256s=args.input_sha256,
            sealed_sizes=args.input_size_bytes,
        )
    import ROOT  # Imported only inside the validated sPHENIX runtime.

    ROOT.gROOT.SetBatch(True)
    payload = reduce_inputs(args.system, inputs, ROOT)
    write_json(args.output, payload)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(2)
