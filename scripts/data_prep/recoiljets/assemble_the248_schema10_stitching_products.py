#!/usr/bin/env python3
"""Assemble schema10 stitching histograms with explicit normalization channels."""

from __future__ import annotations

import argparse
import copy
import datetime as dt
import hashlib
import json
import math
import os
from collections import defaultdict
from pathlib import Path
from typing import Any

try:
    from scripts.data_prep.recoiljets.reduce_the248_schema10_stitching import (
        SAMPLE_CONTRACTS,
        SCHEMA as REDUCER_SCHEMA,
        merge_histograms as merge_reducer_histograms,
    )
    from scripts.data_prep.recoiljets.pp_schema10_weighting import (
        normalize_stitching_sample,
        normalization_contract,
    )
    from scripts.data_prep.recoiljets.auau_embedded_inclusive_schema10_weighting import (
        load_source_stitch_artifact,
    )
except ModuleNotFoundError:  # Support direct execution from this script's directory.
    from reduce_the248_schema10_stitching import (  # type: ignore[no-redef]
        SAMPLE_CONTRACTS,
        SCHEMA as REDUCER_SCHEMA,
        merge_histograms as merge_reducer_histograms,
    )
    from pp_schema10_weighting import (  # type: ignore[no-redef]
        normalize_stitching_sample,
        normalization_contract,
    )
    from auau_embedded_inclusive_schema10_weighting import (  # type: ignore[no-redef]
        load_source_stitch_artifact,
    )


REPO = Path(__file__).resolve().parents[3]
OUTPUT_ROOT = (
    REPO
    / "dataOutput/the243_golden_ppg_analysis_closure_deck_20260821"
    / "stitching_updates_current_20260820"
)
EMBEDDED_INCLUSIVE_ASSEMBLY = (
    OUTPUT_ROOT / "THE248_SCHEMA10_EMBEDDED_INCLUSIVE_STITCHING_ASSEMBLED_V1.json"
)
EMBEDDED_INCLUSIVE_RECEIPT = (
    OUTPUT_ROOT / "THE248_SCHEMA10_EMBEDDED_INCLUSIVE_STITCHING_ASSEMBLY_RECEIPT_V1.json"
)
ASSEMBLED = OUTPUT_ROOT / "THE248_SCHEMA10_STITCHING_ASSEMBLED_V4.json"
RECEIPT = OUTPUT_ROOT / "THE248_SCHEMA10_STITCHING_ASSEMBLY_RECEIPT_V4.json"
EXPECTED_INPUTS = 20006
EXPECTED_BULK_PRODUCTS = 408
EXPECTED_SAMPLES = 14
BIN_COUNT = 400
BIN_WIDTH_GEV = 0.25
EVENT_COUNT_FIELDS = (
    "generated_events",
    "events_with_observable",
    "events_in_histogram_range",
    "events_underflow",
    "events_overflow",
    "events_missing_observable",
    "events_with_finite_weight",
    "events_with_weight_fallback",
)
EVENT_WEIGHT_FIELDS = ("sum_event_weights", "sum_event_weights2")
TRUTH_COUNT_FIELDS = ("records_read", "eligible_records")
AUDIT_COUNT_FIELDS = (
    "entries_dropped",
    "events_in_window",
    "events_outside_window",
    "events_on_lower_edge",
    "events_on_upper_edge",
)


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"JSON object required: {path}")
    return value


def write_new(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o640)
    with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")


def _finite_number(value: Any, label: str, *, nonnegative: bool = False) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{label} is not numeric")
    try:
        numeric = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{label} is not numeric") from error
    if not math.isfinite(numeric) or (nonnegative and numeric < 0.0):
        raise ValueError(f"{label} is invalid")
    return numeric


def _nonnegative_int(value: Any, label: str) -> int:
    if type(value) is not int or value < 0:
        raise ValueError(f"{label} is not a non-negative integer")
    return value


def _expected_ownership(sample_id: str) -> dict[str, Any]:
    contract = SAMPLE_CONTRACTS[sample_id]
    high_is_unbounded = math.isinf(contract.ownership_high_gev)
    return {
        "low_gev": contract.ownership_low_gev,
        "high_gev": None if high_is_unbounded else contract.ownership_high_gev,
        "high_is_unbounded": high_is_unbounded,
        "lower_edge_inclusive": True,
        "upper_edge_inclusive": True,
    }


def _expected_reader_contract(sample_id: str) -> dict[str, Any]:
    photon = SAMPLE_CONTRACTS[sample_id].observable == "max_truth_photon_pt"
    return {
        "directory": "ReplayFoundationV1",
        "event_tree": "RJEventV1",
        "event_branches": ["event_id_hi", "event_id_lo", "event_weight"],
        "observable_tree": "RJTruthPhotonV1" if photon else "RJTruthJetV1",
        "observable_branches": (
            ["event_id_hi", "event_id_lo", "pt", "source_role"]
            if photon
            else ["event_id_hi", "event_id_lo", "algorithm", "radius", "pt"]
        ),
        "dst_reads": 0,
        "ttree_regeneration": False,
    }


def validate_payload(payload: dict[str, Any], sample_id: str, expected_inputs: list[str]) -> None:
    if sample_id not in SAMPLE_CONTRACTS:
        raise ValueError(f"unsupported stitching sample {sample_id}")
    contract = SAMPLE_CONTRACTS[sample_id]
    if (
        payload.get("schema") != REDUCER_SCHEMA
        or payload.get("status") != "PASS"
        or payload.get("sample_id") != sample_id
        or payload.get("system") != contract.system
        or payload.get("source_class") != contract.source_class
        or payload.get("observable") != contract.observable
        or payload.get("ownership_window") != _expected_ownership(sample_id)
        or payload.get("reader_contract") != _expected_reader_contract(sample_id)
        or payload.get("dst_reads") != 0
        or payload.get("ttree_regeneration") is not False
    ):
        raise ValueError(f"stitching payload contract differs for {sample_id}")
    expected_binning = {
        "minimum_gev": 0.0,
        "maximum_gev": 100.0,
        "width_gev": BIN_WIDTH_GEV,
        "edges_gev": [index * BIN_WIDTH_GEV for index in range(BIN_COUNT + 1)],
    }
    if payload.get("binning") != expected_binning:
        raise ValueError(f"fine-bin axis differs for {sample_id}")
    edges = payload["binning"]["edges_gev"]
    expected_edges = [index * 0.25 for index in range(401)]
    if edges != expected_edges:
        raise ValueError(f"fine-bin axis differs for {sample_id}")
    inputs = payload.get("inputs")
    if not isinstance(inputs, list):
        raise ValueError(f"input provenance differs for {sample_id}")
    reported_inputs: list[str] = []
    for index, row in enumerate(inputs):
        if not isinstance(row, dict) or not isinstance(row.get("path"), str):
            raise ValueError(f"input provenance differs for {sample_id}")
        _nonnegative_int(row.get("size_bytes"), f"inputs[{index}].size_bytes")
        if row["size_bytes"] == 0:
            raise ValueError(f"input provenance differs for {sample_id}")
        reported_inputs.append(row["path"])
    if len(reported_inputs) != len(set(reported_inputs)):
        raise ValueError(f"duplicate input provenance for {sample_id}")
    if reported_inputs != expected_inputs:
        raise ValueError(f"input provenance differs for {sample_id}")
    if payload.get("input_count") != len(expected_inputs):
        raise ValueError(f"input count differs for {sample_id}")
    histogram = payload.get("histogram", {})
    raw_counts = histogram.get("raw_counts")
    sumw = histogram.get("sumw")
    sumw2 = histogram.get("sumw2")
    if not isinstance(raw_counts, list) or len(raw_counts) != BIN_COUNT:
        raise ValueError(f"raw_counts shape differs for {sample_id}")
    if not isinstance(sumw, list) or len(sumw) != BIN_COUNT:
        raise ValueError(f"sumw shape differs for {sample_id}")
    if not isinstance(sumw2, list) or len(sumw2) != BIN_COUNT:
        raise ValueError(f"sumw2 shape differs for {sample_id}")
    for index, value in enumerate(raw_counts):
        _nonnegative_int(value, f"raw_counts[{index}]")
    for index, value in enumerate(sumw):
        _finite_number(value, f"sumw[{index}]")
    for index, value in enumerate(sumw2):
        _finite_number(value, f"sumw2[{index}]", nonnegative=True)

    flows: dict[str, dict[str, Any]] = {}
    for name in ("underflow", "overflow", "missing"):
        flow = histogram.get(name)
        if not isinstance(flow, dict):
            raise ValueError(f"{name} differs for {sample_id}")
        _nonnegative_int(flow.get("raw_count"), f"{name}.raw_count")
        _finite_number(flow.get("sumw"), f"{name}.sumw")
        _finite_number(flow.get("sumw2"), f"{name}.sumw2", nonnegative=True)
        flows[name] = flow

    totals = payload.get("event_totals")
    if not isinstance(totals, dict):
        raise ValueError(f"event totals differ for {sample_id}")
    for field in EVENT_COUNT_FIELDS:
        _nonnegative_int(totals.get(field), f"event_totals.{field}")
    _finite_number(totals.get("sum_event_weights"), "event_totals.sum_event_weights")
    _finite_number(
        totals.get("sum_event_weights2"),
        "event_totals.sum_event_weights2",
        nonnegative=True,
    )
    if totals["generated_events"] != totals["events_with_observable"] + totals["events_missing_observable"]:
        raise ValueError(f"event observable partition differs for {sample_id}")
    if totals["generated_events"] != totals["events_with_finite_weight"] + totals["events_with_weight_fallback"]:
        raise ValueError(f"event weight partition differs for {sample_id}")
    if totals["events_with_observable"] != (
        totals["events_in_histogram_range"] + totals["events_underflow"] + totals["events_overflow"]
    ):
        raise ValueError(f"histogram event partition differs for {sample_id}")
    if sum(raw_counts) != totals["events_in_histogram_range"]:
        raise ValueError(f"histogram raw count differs for {sample_id}")
    if flows["underflow"]["raw_count"] != totals["events_underflow"]:
        raise ValueError(f"underflow count differs for {sample_id}")
    if flows["overflow"]["raw_count"] != totals["events_overflow"]:
        raise ValueError(f"overflow count differs for {sample_id}")
    if flows["missing"]["raw_count"] != totals["events_missing_observable"]:
        raise ValueError(f"missing count differs for {sample_id}")
    summed_weight = sum(float(value) for value in sumw) + sum(float(flow["sumw"]) for flow in flows.values())
    summed_weight2 = sum(float(value) for value in sumw2) + sum(float(flow["sumw2"]) for flow in flows.values())
    if not math.isclose(summed_weight, float(totals["sum_event_weights"]), rel_tol=1e-12, abs_tol=1e-9):
        raise ValueError(f"event weight sum differs for {sample_id}")
    if not math.isclose(summed_weight2, float(totals["sum_event_weights2"]), rel_tol=1e-12, abs_tol=1e-9):
        raise ValueError(f"event weight-square sum differs for {sample_id}")

    audit = payload.get("ownership_audit")
    if not isinstance(audit, dict) or audit.get("application") != "audit_only_production_ownership_preserved":
        raise ValueError(f"ownership audit differs for {sample_id}")
    for field in AUDIT_COUNT_FIELDS:
        _nonnegative_int(audit.get(field), f"ownership_audit.{field}")
    if audit["entries_dropped"] != 0:
        raise ValueError(f"ownership audit dropped entries for {sample_id}")
    if audit["events_in_window"] + audit["events_outside_window"] != totals["events_with_observable"]:
        raise ValueError(f"ownership audit partition differs for {sample_id}")

    truth_totals = payload.get("truth_record_totals")
    if not isinstance(truth_totals, dict):
        raise ValueError(f"truth-record totals differ for {sample_id}")
    for field in TRUTH_COUNT_FIELDS:
        _nonnegative_int(truth_totals.get(field), f"truth_record_totals.{field}")
    if truth_totals["eligible_records"] > truth_totals["records_read"]:
        raise ValueError(f"truth-record partition differs for {sample_id}")


def merge_sample(sample_id: str, payloads: list[dict[str, Any]]) -> dict[str, Any]:
    if not payloads:
        raise ValueError(f"no products for {sample_id}")
    merged_payload = merge_reducer_histograms(payloads)
    merged_inputs = [str(row["path"]) for row in merged_payload["inputs"]]
    validate_payload(merged_payload, sample_id, merged_inputs)
    totals = copy.deepcopy(merged_payload["event_totals"])
    merged = {
        "sample_id": sample_id,
        "system": merged_payload["system"],
        "source_class": merged_payload["source_class"],
        "observable": merged_payload["observable"],
        "bin_edges_gev": list(merged_payload["binning"]["edges_gev"]),
        "raw_counts": list(merged_payload["histogram"]["raw_counts"]),
        "sumw": list(merged_payload["histogram"]["sumw"]),
        "sumw2": list(merged_payload["histogram"]["sumw2"]),
        "events_total": totals["generated_events"],
        "events_with_observable": totals["events_with_observable"],
        "events_missing_observable": totals["events_missing_observable"],
        "nonfinite_event_weight_fallbacks": totals["events_with_weight_fallback"],
        "event_totals": totals,
        "truth_record_totals": copy.deepcopy(merged_payload["truth_record_totals"]),
        "underflow": copy.deepcopy(merged_payload["histogram"]["underflow"]),
        "overflow": copy.deepcopy(merged_payload["histogram"]["overflow"]),
        "missing": copy.deepcopy(merged_payload["histogram"]["missing"]),
        "ownership_window": copy.deepcopy(merged_payload["ownership_window"]),
        "ownership_audit": copy.deepcopy(merged_payload["ownership_audit"]),
        "reader_contract": copy.deepcopy(merged_payload["reader_contract"]),
        "input_count": merged_payload["input_count"],
        "inputs": copy.deepcopy(merged_payload["inputs"]),
        "dst_reads": 0,
        "ttree_regeneration": False,
    }
    return merged


def apply_catalog_weight(merged: dict[str, Any], contract: dict[str, Any]) -> dict[str, Any]:
    """Compatibility entrypoint for canonical schema10 dual-channel scaling."""
    edges = merged.get("bin_edges_gev")
    if not isinstance(edges, list) or len(edges) != BIN_COUNT + 1:
        raise ValueError(f"{merged.get('sample_id')} bin edges differ")
    widths = [float(right) - float(left) for left, right in zip(edges, edges[1:])]
    if any(not math.isclose(width, BIN_WIDTH_GEV, rel_tol=0.0, abs_tol=1e-12) for width in widths):
        raise ValueError(f"{merged.get('sample_id')} bin width differs")
    return normalize_stitching_sample(merged, contract, bin_width_gev=BIN_WIDTH_GEV)


def assemble(
    *, catalog_path: Path, job_plan_path: Path, proof_path: Path,
    bulk_root: Path, acquisition_receipt_path: Path, output_path: Path, receipt_path: Path,
    embedded_inclusive_assembly_path: Path,
    embedded_inclusive_receipt_path: Path,
) -> dict[str, Any]:
    for path in (
        catalog_path,
        job_plan_path,
        proof_path,
        acquisition_receipt_path,
        embedded_inclusive_assembly_path,
        embedded_inclusive_receipt_path,
    ):
        if path.is_symlink() or not path.is_file():
            raise FileNotFoundError(path)
    if output_path.exists() or receipt_path.exists():
        raise FileExistsError("THE-248 assembled output identity already exists")
    catalog = load(catalog_path)
    plan = load(job_plan_path)
    proof = load(proof_path)
    acquisition = load(acquisition_receipt_path)
    embedded_source = load_source_stitch_artifact(
        embedded_inclusive_assembly_path,
        embedded_inclusive_receipt_path,
    )
    embedded_assembly = load(embedded_inclusive_assembly_path)
    if (
        catalog.get("schema") != "THE243Schema10AcceptedResponseCatalogV1"
        or catalog.get("status") != "PASS"
        or sha(catalog_path) != plan.get("accepted_catalog_sha256")
    ):
        raise ValueError("accepted catalog binding differs")
    if (
        plan.get("schema") != "THE248Schema10StitchingJobPlanV2"
        or plan.get("status") != "PASS"
        or plan.get("exact_jobs") != EXPECTED_BULK_PRODUCTS
    ):
        raise ValueError("THE-248 job plan differs")
    if (
        proof.get("schema") != "THE248Schema10StitchingForegroundPairProofV1"
        or proof.get("status") != "PASS"
        or proof.get("accepted_catalog_sha256") != sha(catalog_path)
    ):
        raise ValueError("foreground pair proof differs")
    if (
        acquisition.get("schema")
        != "THE248Schema10StitchingAcquisitionTerminalReceiptV1"
        or acquisition.get("status") != "PASS"
        or acquisition.get("job_plan_sha256") != sha(job_plan_path)
        or acquisition.get("local_root") != bulk_root.as_posix()
        or acquisition.get("member_count") != 3 * EXPECTED_BULK_PRODUCTS
        or acquisition.get("output_json_count") != EXPECTED_BULK_PRODUCTS
        or acquisition.get("measurement_json_count") != EXPECTED_BULK_PRODUCTS
        or acquisition.get("time_evidence_count") != EXPECTED_BULK_PRODUCTS
        or acquisition.get("bulk_input_root_count") != EXPECTED_INPUTS - 2
        or not isinstance(acquisition.get("products"), list)
        or len(acquisition["products"]) != EXPECTED_BULK_PRODUCTS
    ):
        raise ValueError("THE-248 acquisition terminal receipt differs")

    by_sample: dict[str, list[dict[str, Any]]] = defaultdict(list)
    input_union: set[str] = set()
    product_rows: list[dict[str, Any]] = []
    for job in plan["jobs"]:
        sample_id = str(job["sample_id"])
        input_list = job_plan_path.parent / str(job["input_list"])
        expected_inputs = [line.strip() for line in input_list.read_text().splitlines() if line.strip()]
        product = (
            bulk_root / "outputs" / sample_id / "shard_0000" / f"{job['job_id']}.json"
        )
        if product.is_symlink() or not product.is_file():
            raise FileNotFoundError(product)
        payload = load(product)
        validate_payload(payload, sample_id, expected_inputs)
        overlap = input_union.intersection(expected_inputs)
        if overlap:
            raise ValueError(f"bulk product input duplication: {sorted(overlap)[:1]}")
        input_union.update(expected_inputs)
        by_sample[sample_id].append(payload)
        product_rows.append({
            "kind": "bulk",
            "sample_id": sample_id,
            "job_id": job["job_id"],
            "path": product.as_posix(),
            "sha256": sha(product),
            "input_count": len(expected_inputs),
            "accepted_response_product_sha256": job["accepted_response_product_sha256"],
        })

    proof_products = proof.get("products", {})
    if set(proof_products) != {"pp_photon10", "auau_jet20"}:
        raise ValueError("foreground proof sample pair differs")
    for sample_id, payload in proof_products.items():
        expected_input = str(plan["foreground_proof_inputs"][sample_id])
        validate_payload(payload, sample_id, [expected_input])
        if expected_input in input_union:
            raise ValueError("foreground proof input leaked into bulk products")
        input_union.add(expected_input)
        by_sample[sample_id].append(payload)
        product_rows.append({
            "kind": "foreground_proof",
            "sample_id": sample_id,
            "path": proof_path.as_posix(),
            "sha256": sha(proof_path),
            "input_count": 1,
            "input_sha256": proof.get("input_sha256", {}).get(sample_id),
            "accepted_response_product_sha256": proof.get("accepted_response_product_sha256", {}).get(sample_id),
        })

    if len(product_rows) != EXPECTED_BULK_PRODUCTS + 2:
        raise ValueError("product count differs")
    if len(input_union) != EXPECTED_INPUTS:
        raise ValueError("complete accepted input union differs")
    if set(by_sample) != set(catalog["samples"]) or len(by_sample) != EXPECTED_SAMPLES:
        raise ValueError("sample inventory differs")

    samples: dict[str, Any] = {}
    for sample_id in sorted(by_sample):
        merged = merge_sample(sample_id, by_sample[sample_id])
        contract = catalog["samples"][sample_id]
        product_hashes = sorted({
            row["accepted_response_product_sha256"]
            for row in product_rows
            if row["sample_id"] == sample_id and row.get("accepted_response_product_sha256")
        })
        embedded_inclusive = (
            contract.get("system") == "auau"
            and contract.get("source_class") == "inclusive_background"
        )
        if embedded_inclusive:
            source_row = embedded_assembly["samples"].get(sample_id)
            if not isinstance(source_row, dict):
                raise ValueError(f"canonical embedded-inclusive sample missing: {sample_id}")
            embedded_source.sample(sample_id)
            for field in ("bin_edges_gev", "raw_counts", "sumw", "sumw2"):
                if source_row.get(field) != merged.get(field):
                    raise ValueError(
                        f"{sample_id}.{field} differs from the bound embedded-inclusive source artifact"
                    )
            if source_row.get("source_product_hashes") != product_hashes:
                raise ValueError(f"{sample_id} source-product hashes differ from canonical artifact")
            weighted = copy.deepcopy(source_row)
            weighted["analysis_weight_state"] = (
                "CANONICAL_SOURCE_STITCH_ONLY__CENTRALITY_NOT_APPLIED"
            )
            weighted["nominal_downstream_analysis_ready"] = False
        else:
            weighted = apply_catalog_weight(merged, contract)
            weighted["source_product_hashes"] = product_hashes
        samples[sample_id] = weighted

    assembled = {
        "schema": "THE248Schema10StitchingAssemblyV4",
        "status": "PASS",
        "created_at": utc_now(),
        "accepted_catalog_path": catalog_path.as_posix(),
        "accepted_catalog_sha256": sha(catalog_path),
        "job_plan_path": job_plan_path.as_posix(),
        "job_plan_sha256": sha(job_plan_path),
        "foreground_pair_proof_path": proof_path.as_posix(),
        "foreground_pair_proof_sha256": sha(proof_path),
        "acquisition_terminal_receipt_path": acquisition_receipt_path.as_posix(),
        "acquisition_terminal_receipt_sha256": sha(acquisition_receipt_path),
        "input_root_count": len(input_union),
        "unique_input_root_count": len(input_union),
        "bulk_product_count": EXPECTED_BULK_PRODUCTS,
        "foreground_product_count": 2,
        "sample_count": len(samples),
        "binning": {"lo_gev": 0.0, "hi_gev": 100.0, "width_gev": 0.25, "bins": 400},
        "normalization": {
            "schema": "THE248Schema10SplitStitchingNormalizationV4",
            "pp_and_noninclusive_contract": normalization_contract(),
            "auau_embedded_inclusive": {
                "formula": "raw_counts * ownership_effective_cross_section_pb / Npass_same_ownership_window / bin_width_gev",
                "source_artifact": embedded_source.binding_payload(),
                "centrality_reweighting_applied": False,
                "downstream_order": [
                    "existing_event_weight",
                    "canonical_source_stitch",
                    "auau_centrality_reweight",
                ],
            },
        },
        "source_raw_assembly": {
            "path": (OUTPUT_ROOT / "THE248_SCHEMA10_STITCHING_ASSEMBLED.json").as_posix(),
            "sha256": "631ec123d157125b880f0af00515b03e97ecc2cf5ce83febb5eed5e0b110e69d",
        },
        "supersedes": [
            {
                "path": (OUTPUT_ROOT / "THE248_SCHEMA10_STITCHING_ASSEMBLED.json").as_posix(),
                "sha256": "631ec123d157125b880f0af00515b03e97ecc2cf5ce83febb5eed5e0b110e69d",
                "reason": "V1 exposed an ambiguous weighted-density field with double-applied pp slice weighting",
            },
            {
                "path": (OUTPUT_ROOT / "THE248_SCHEMA10_STITCHING_ASSEMBLED_V2.json").as_posix(),
                "sha256": "5fcaa9baf48eb1418b9723f4105a78d081eba23c79480539af7f38a14c6c4ab3",
                "reason": "V2 used the downstream analysis-weight channel as the generator-stitching plot channel",
            },
            {
                "path": (OUTPUT_ROOT / "THE248_SCHEMA10_STITCHING_ASSEMBLED_V3.json").as_posix(),
                "sha256": "0006146ad2876bf7d0dfee7a7073138a250c43fbc56f31eb01b36e7d244c8705",
                "reason": "V3 normalized Au+Au embedded inclusive samples with cross section / generated events",
            },
        ],
        "reader_contract": {
            "source": "accepted ReplayFoundationV1 schema-10 TTrees",
            "dst_reads": 0,
            "ttree_regeneration": False,
            "photon_observable": "maximum finite positive source_role==1 RJTruthPhotonV1 pt per event",
            "jet_observable": "maximum finite positive antikt R=0.4 RJTruthJetV1 pt per event",
        },
        "products": sorted(product_rows, key=lambda row: (row["sample_id"], row["kind"], row.get("job_id", ""))),
        "samples": samples,
    }
    write_new(output_path, assembled)
    receipt = {
        "schema": "THE248Schema10StitchingAssemblyTerminalReceiptV4",
        "status": "PASS",
        "created_at": utc_now(),
        "assembled_path": output_path.as_posix(),
        "assembled_sha256": sha(output_path),
        "accepted_catalog_sha256": sha(catalog_path),
        "input_root_count": len(input_union),
        "sample_count": len(samples),
        "bulk_product_count": EXPECTED_BULK_PRODUCTS,
        "foreground_product_count": 2,
        "acquisition_terminal_receipt_sha256": sha(acquisition_receipt_path),
        "dst_reads": 0,
        "ttree_regeneration": False,
        "normalization_schema": "THE248Schema10SplitStitchingNormalizationV4",
        "auau_embedded_inclusive_source_stitch": embedded_source.binding_payload(),
        "source_raw_assembly_sha256": "631ec123d157125b880f0af00515b03e97ecc2cf5ce83febb5eed5e0b110e69d",
        "supersedes_assembled_sha256": [
            "631ec123d157125b880f0af00515b03e97ecc2cf5ce83febb5eed5e0b110e69d",
            "5fcaa9baf48eb1418b9723f4105a78d081eba23c79480539af7f38a14c6c4ab3",
            "0006146ad2876bf7d0dfee7a7073138a250c43fbc56f31eb01b36e7d244c8705",
        ],
    }
    write_new(receipt_path, receipt)
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", type=Path, required=True)
    parser.add_argument("--job-plan", type=Path, required=True)
    parser.add_argument("--proof", type=Path, required=True)
    parser.add_argument("--bulk-root", type=Path, required=True)
    parser.add_argument("--acquisition-receipt", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=ASSEMBLED)
    parser.add_argument("--receipt", type=Path, default=RECEIPT)
    parser.add_argument(
        "--embedded-inclusive-assembly",
        type=Path,
        default=EMBEDDED_INCLUSIVE_ASSEMBLY,
    )
    parser.add_argument(
        "--embedded-inclusive-receipt",
        type=Path,
        default=EMBEDDED_INCLUSIVE_RECEIPT,
    )
    args = parser.parse_args()
    result = assemble(
        catalog_path=args.catalog.resolve(), job_plan_path=args.job_plan.resolve(),
        proof_path=args.proof.resolve(), bulk_root=args.bulk_root.resolve(),
        acquisition_receipt_path=args.acquisition_receipt.resolve(),
        output_path=args.output.resolve(), receipt_path=args.receipt.resolve(),
        embedded_inclusive_assembly_path=args.embedded_inclusive_assembly.resolve(),
        embedded_inclusive_receipt_path=args.embedded_inclusive_receipt.resolve(),
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
