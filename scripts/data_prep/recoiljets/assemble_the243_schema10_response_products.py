#!/usr/bin/env python3
"""Validate and merge the accepted THE-243 response accumulators by sample."""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import math
import os
from collections import defaultdict
from pathlib import Path
import sys
from typing import Any, Mapping


REPO = Path(__file__).resolve().parents[3]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.data_prep.recoiljets.reduce_the243_schema10_response import merge_response_fanouts
from scripts.data_prep.recoiljets.reduce_the243_schema10_xjgamma import photon_id_contract
from scripts.data_prep.recoiljets.auau_embedded_inclusive_schema10_weighting import (
    SourceStitchArtifact,
    load_source_stitch_artifact,
)


EXPECTED_SCHEMA = "THE243Schema10InclusiveResponseThresholdFanoutV1"
EXPECTED_THRESHOLDS = [5.0, 7.0, 10.0, 12.0]
EXPECTED_CANARY_JOBS = 2
EXPECTED_FULL_JOBS = 290
EXPECTED_INPUTS = 7145
DIRECT_GROUP_PREFIX = "/direct/sphenix+tg+tg01/"
CANONICAL_GROUP_PREFIX = "/sphenix/tg/tg01/"


def canonical_input_path(value: str) -> str:
    """Map the SDCC direct-I/O alias back to the frozen group path identity."""
    if value.startswith(DIRECT_GROUP_PREFIX):
        return CANONICAL_GROUP_PREFIX + value[len(DIRECT_GROUP_PREFIX):]
    return value


def add_nested(*values: Any) -> Any:
    """Add equal-shaped numeric JSON arrays without flattening their axes."""
    if all(isinstance(value, list) for value in values):
        lengths = {len(value) for value in values}
        if len(lengths) != 1:
            raise ValueError("response fake-partition shapes differ")
        return [add_nested(*(value[index] for value in values))
                for index in range(next(iter(lengths)))]
    if any(isinstance(value, list) for value in values):
        raise ValueError("response fake-partition ranks differ")
    return sum(float(value) for value in values)


def equal_nested(left: Any, right: Any, tolerance: float = 1e-9) -> bool:
    if isinstance(left, list) and isinstance(right, list):
        return len(left) == len(right) and all(
            equal_nested(a, b, tolerance) for a, b in zip(left, right)
        )
    if isinstance(left, list) or isinstance(right, list):
        return False
    return abs(float(left) - float(right)) <= tolerance * max(1.0, abs(float(left)), abs(float(right)))


def validate_fake_partition(reduction: Mapping[str, Any], path: Path) -> None:
    """Require each reconstructed fake to have exactly one physical cause."""
    xj = reduction.get("xj")
    audit = reduction.get("audit")
    if not isinstance(xj, Mapping) or not isinstance(audit, Mapping):
        raise ValueError(f"response fake-partition metadata missing: {path}")
    fields = (
        "fakes", "fakes_sumw2",
        "fake_photon_reco", "fake_photon_reco_sumw2",
        "combinatoric_reco", "combinatoric_reco_sumw2",
        "detector_fakes_reco", "detector_fakes_reco_sumw2",
    )
    if any(field not in xj for field in fields):
        raise ValueError(f"response fake-partition fields missing: {path}")
    expected = add_nested(
        xj["fake_photon_reco"], xj["combinatoric_reco"], xj["detector_fakes_reco"]
    )
    expected_sumw2 = add_nested(
        xj["fake_photon_reco_sumw2"],
        xj["combinatoric_reco_sumw2"],
        xj["detector_fakes_reco_sumw2"],
    )
    if not equal_nested(xj["fakes"], expected) or not equal_nested(xj["fakes_sumw2"], expected_sumw2):
        raise ValueError(f"response fake partition is not additive: {path}")
    audit_total = int(audit.get("xj_reco_fake", 0))
    audit_parts = sum(int(audit.get(field, 0)) for field in (
        "xj_fake_photon_reco", "xj_combinatoric_reco", "xj_detector_fake_reco"
    ))
    if audit_parts != audit_total:
        raise ValueError(f"response fake audit partition differs: {path}")


def utc_now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")


def load(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"JSON object required: {path}")
    return value


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_exclusive(path: Path, value: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o640)
    with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2, sort_keys=True)
        stream.write("\n")


def validate_payload(payload: Mapping[str, Any], row: Mapping[str, Any], path: Path) -> None:
    if (payload.get("schema") != EXPECTED_SCHEMA
            or payload.get("result_status") != "RESPONSE_INPUT__NOT_UNFOLDED"
            or payload.get("system") != row.get("system")
            or payload.get("jet_pt_thresholds_gev") != EXPECTED_THRESHOLDS):
        raise ValueError(f"response contract differs: {path}")
    inputs = payload.get("inputs")
    expected_inputs = int(row.get("input_count", 1))
    if not isinstance(inputs, list) or len(inputs) != expected_inputs:
        raise ValueError(f"response input coverage differs: {path}")
    payload_paths = [
        canonical_input_path(str(item.get("path")))
        for item in inputs if isinstance(item, Mapping)
    ]
    if payload_paths != row.get("expected_input_paths"):
        raise ValueError(f"response input identity/order differs: {path}")
    provenance = payload.get("model_provenance", {})
    if (provenance.get("model_sha256") != row.get("model_sha256")
            or provenance.get("shower_definition") != "H70"
            or provenance.get("h0_fallbacks") != 0):
        raise ValueError(f"response H70 provenance differs: {path}")
    expected_photon_id = photon_id_contract(str(row.get("system")))
    if payload.get("photon_id_contract") != expected_photon_id:
        raise ValueError(f"response photon-ID contract differs: {path}")
    reductions = payload.get("reductions")
    if not isinstance(reductions, Mapping) or set(reductions) != {"pt5", "pt7", "pt10", "pt12"}:
        raise ValueError(f"response threshold reductions differ: {path}")
    for reduction in reductions.values():
        if (not isinstance(reduction, Mapping)
                or reduction.get("photon_id_contract") != expected_photon_id):
            raise ValueError(f"response threshold photon-ID contract differs: {path}")
        validate_fake_partition(reduction, path)


def rows_from(plan_path: Path, raw_root: Path, expected_jobs: int) -> list[tuple[dict[str, Any], Path]]:
    plan = load(plan_path)
    jobs = plan.get("jobs")
    if not isinstance(jobs, list) or len(jobs) != expected_jobs:
        raise ValueError(f"job-plan cardinality differs: {plan_path}")
    remote_base = Path(str(plan.get("remote_base")))
    result: list[tuple[dict[str, Any], Path]] = []
    for row in jobs:
        if not isinstance(row, dict):
            raise ValueError(f"job-plan row differs: {plan_path}")
        relative = Path(str(row["output_path"])).relative_to(remote_base)
        normalized = dict(row)
        normalized["input_count"] = int(row.get("input_count", 1))
        if row.get("input_list") is None:
            if normalized["input_count"] != 1 or not isinstance(row.get("input_path"), str):
                raise ValueError(f"response direct-input binding differs: {plan_path}")
            expected_input_paths = [str(row["input_path"])]
        else:
            input_list = plan_path.parent / str(row.get("input_list"))
            if (input_list.is_symlink() or not input_list.is_file()
                    or sha(input_list) != row.get("input_list_sha256")):
                raise ValueError(f"response input-list binding differs: {input_list}")
            expected_input_paths = [
                line.strip() for line in input_list.read_text(encoding="utf-8").splitlines()
                if line.strip()
            ]
            if (len(expected_input_paths) != normalized["input_count"]
                    or len(set(expected_input_paths)) != len(expected_input_paths)
                    or expected_input_paths[0] != row.get("first_input")
                    or expected_input_paths[-1] != row.get("last_input")):
                raise ValueError(f"response input-list content differs: {input_list}")
        normalized["expected_input_paths"] = expected_input_paths
        result.append((normalized, raw_root / relative))
    return result


def sample_contracts_from(plan_path: Path) -> dict[str, dict[str, Any]]:
    plan = load(plan_path)
    rows = plan.get("sample_contracts")
    if not isinstance(rows, list):
        raise ValueError(f"sample contracts missing: {plan_path}")
    result: dict[str, dict[str, Any]] = {}
    for row in rows:
        if not isinstance(row, dict) or not isinstance(row.get("sample_id"), str):
            raise ValueError(f"sample contract differs: {plan_path}")
        result[row["sample_id"]] = {
            "sample_id": row["sample_id"],
            "system": row["system"],
            "source_class": row["source_class"],
            "cross_section_pb": float(row["cross_section_pb"]),
            "production_campaign_tag": row["production_campaign_tag"],
        }
    return result


def bind_sample_contract(row: Mapping[str, Any], contracts: Mapping[str, Mapping[str, Any]]) -> dict[str, Any]:
    """Fill legacy canary labels from the frozen full-plan sample contract.

    The accepted two-row canary predates the full-plan bookkeeping fields but
    names the same sample IDs and production campaign.  The full plan is the
    immutable authority for source class and cross section; physics payloads
    are not modified.
    """
    sample_id = str(row.get("sample_id"))
    contract = contracts.get(sample_id)
    if not isinstance(contract, Mapping):
        raise ValueError(f"sample contract missing for {sample_id}")
    normalized = dict(row)
    for field in ("system", "source_class", "cross_section_pb", "production_campaign_tag"):
        existing = normalized.get(field)
        expected = contract[field]
        if existing is not None and existing != expected:
            raise ValueError(f"sample contract drift for {sample_id}: {field}")
        normalized[field] = expected
    return normalized


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--canary-job-plan", type=Path)
    parser.add_argument("--canary-root", type=Path)
    parser.add_argument("--full-job-plan", type=Path, required=True)
    parser.add_argument("--full-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument(
        "--auau-embedded-inclusive-stitching-assembly",
        type=Path,
        help="accepted sigma_eff/Npass assembly; required when Au+Au inclusive samples are present",
    )
    parser.add_argument(
        "--auau-embedded-inclusive-stitching-receipt",
        type=Path,
        help="validation receipt for the accepted sigma_eff/Npass assembly",
    )
    args = parser.parse_args()

    contracts = sample_contracts_from(args.full_job_plan.resolve())
    has_embedded_inclusive = any(
        row.get("system") == "auau" and row.get("source_class") == "inclusive_background"
        for row in contracts.values()
    )
    source_stitch: SourceStitchArtifact | None = None
    if has_embedded_inclusive:
        if (
            args.auau_embedded_inclusive_stitching_assembly is None
            or args.auau_embedded_inclusive_stitching_receipt is None
        ):
            raise ValueError(
                "Au+Au embedded inclusive response assembly requires the accepted "
                "hash-bound sigma_eff/Npass assembly and validation receipt"
            )
        source_stitch = load_source_stitch_artifact(
            args.auau_embedded_inclusive_stitching_assembly,
            args.auau_embedded_inclusive_stitching_receipt,
        )
    if (args.canary_job_plan is None) != (args.canary_root is None):
        raise ValueError("canary job plan and root must be supplied together")
    if args.canary_job_plan is None or args.canary_root is None:
        raise ValueError("accepted inclusive-z10 response canary products are required")
    planned_raw = [
        ("responsek_exclusiveleakage_canary_s10", *item)
        for item in rows_from(
            args.canary_job_plan.resolve(), args.canary_root.resolve(), EXPECTED_CANARY_JOBS
        )
    ] + [
        ("responsek_exclusiveleakage_complement_s12", *item)
        for item in rows_from(
            args.full_job_plan.resolve(), args.full_root.resolve(), EXPECTED_FULL_JOBS
        )
    ]
    planned = [(namespace, bind_sample_contract(row, contracts), path)
               for namespace, row, path in planned_raw]
    by_sample: dict[str, list[dict[str, Any]]] = defaultdict(list)
    sample_contracts: dict[str, dict[str, Any]] = {}
    catalog_rows: list[dict[str, Any]] = []
    input_paths: list[str] = []
    for namespace, row, path in planned:
        if path.is_symlink() or not path.is_file():
            raise ValueError(f"accepted response output missing: {path}")
        payload = load(path)
        validate_payload(payload, row, path)
        sample_id = str(row["sample_id"])
        contract = {
            "sample_id": sample_id,
            "system": row["system"],
            "source_class": row["source_class"],
            "cross_section_pb": float(row["cross_section_pb"]),
            "production_campaign_tag": row["production_campaign_tag"],
        }
        if sample_id in sample_contracts and sample_contracts[sample_id] != contract:
            raise ValueError(f"sample contract drift: {sample_id}")
        sample_contracts[sample_id] = contract
        by_sample[sample_id].append(payload)
        paths = [canonical_input_path(str(item.get("path"))) for item in payload["inputs"]]
        input_paths.extend(paths)
        catalog_rows.append({
            "namespace": namespace,
            "sample_id": sample_id,
            "job_id": row["job_id"],
            "system": row["system"],
            "source_class": row["source_class"],
            "input_count": len(paths),
            "path": path.as_posix(),
            "size_bytes": path.stat().st_size,
            "sha256": sha(path),
        })
    if len(catalog_rows) != EXPECTED_CANARY_JOBS + EXPECTED_FULL_JOBS:
        raise ValueError("accepted response output union differs")
    if len(input_paths) != EXPECTED_INPUTS or len(set(input_paths)) != EXPECTED_INPUTS:
        raise ValueError("accepted response input union is incomplete or duplicated")

    merged_receipts: dict[str, Any] = {}
    for sample_id in sorted(by_sample):
        merged = merge_response_fanouts(by_sample[sample_id])
        contract = sample_contracts[sample_id]
        generated_events = int(merged["reductions"]["pt5"]["record_counts"].get("events", 0))
        if generated_events <= 0:
            raise ValueError(f"sample has no generated events: {sample_id}")
        merged["sample_contract"] = contract
        merged["source_output_count"] = len(by_sample[sample_id])
        merged["generated_events"] = generated_events
        embedded_inclusive = (
            contract["system"] == "auau"
            and contract["source_class"] == "inclusive_background"
        )
        if embedded_inclusive:
            if source_stitch is None:
                raise ValueError("canonical Au+Au embedded inclusive source stitch is absent")
            stitch = source_stitch.sample(sample_id)
            if not math.isclose(
                float(contract["cross_section_pb"]),
                stitch.ownership_effective_cross_section_pb,
                rel_tol=0.0,
                abs_tol=1.0e-9,
            ):
                raise ValueError(f"{sample_id} effective cross section differs from source authority")
            merged.pop("cross_section_weight_pb_per_event", None)
            merged["canonical_source_stitch_weight_pb_per_owned_event"] = (
                stitch.stitching_weight_pb_per_owned_event
            )
            merged["canonical_source_stitch_denominator_events"] = (
                stitch.normalization_denominator_events
            )
            merged["canonical_source_stitch_ownership_window_gev"] = [
                stitch.ownership_low_gev,
                stitch.ownership_high_gev,
            ]
            merged["analysis_weight_state"] = (
                "SOURCE_STITCH_BOUND__CENTRALITY_NOT_APPLIED__NOT_NOMINAL_DOWNSTREAM"
            )
            merged["source_stitch_artifact"] = source_stitch.binding_payload()
        else:
            if contract["system"] == "auau" and contract["source_class"] != "photon_signal":
                raise ValueError(f"unsupported Au+Au response source class: {sample_id}")
            merged["cross_section_weight_pb_per_event"] = (
                contract["cross_section_pb"] / generated_events
            )
        destination = args.output_root.resolve() / "samples" / f"{sample_id}.json"
        write_exclusive(destination, merged)
        merged_receipts[sample_id] = {
            **contract,
            "generated_events": generated_events,
            "source_outputs": len(by_sample[sample_id]),
            "output_path": destination.as_posix(),
            "output_sha256": sha(destination),
            "output_size_bytes": destination.stat().st_size,
        }
        if embedded_inclusive:
            merged_receipts[sample_id].update({
                "canonical_source_stitch_weight_pb_per_owned_event": (
                    merged["canonical_source_stitch_weight_pb_per_owned_event"]
                ),
                "canonical_source_stitch_denominator_events": (
                    merged["canonical_source_stitch_denominator_events"]
                ),
                "analysis_weight_state": merged["analysis_weight_state"],
            })
        else:
            merged_receipts[sample_id]["cross_section_weight_pb_per_event"] = (
                merged["cross_section_weight_pb_per_event"]
            )

    catalog = {
        "schema": "THE243Schema10AcceptedResponseCatalogV2",
        "status": "PASS",
        "created_at": utc_now(),
        "output_json_count": len(catalog_rows),
        "input_root_count": len(input_paths),
        "unique_input_root_count": len(set(input_paths)),
        "sample_count": len(merged_receipts),
        "outputs": sorted(catalog_rows, key=lambda row: (row["sample_id"], row["job_id"])),
        "samples": merged_receipts,
        "auau_embedded_inclusive_source_stitch": (
            source_stitch.binding_payload() if source_stitch is not None else None
        ),
    }
    catalog_path = args.output_root.resolve() / "THE243_SCHEMA10_ACCEPTED_RESPONSE_CATALOG.json"
    write_exclusive(catalog_path, catalog)
    terminal = {
        "schema": "THE243Schema10AcceptedResponseAssemblyTerminalReceiptV2",
        "status": "PASS",
        "created_at": utc_now(),
        "catalog_path": catalog_path.as_posix(),
        "catalog_sha256": sha(catalog_path),
        "output_json_count": len(catalog_rows),
        "input_root_count": len(input_paths),
        "sample_count": len(merged_receipts),
        "jet_pt_thresholds_gev": EXPECTED_THRESHOLDS,
        "result_stage": "RESPONSE_INPUT__NOT_UNFOLDED",
        "dst_reads": 0,
        "accepted_canary_outputs": 2,
        "accepted_canary_inputs": 2,
        "samples": merged_receipts,
        "auau_embedded_inclusive_source_stitch": (
            source_stitch.binding_payload() if source_stitch is not None else None
        ),
    }
    terminal_path = args.output_root.resolve() / "THE243_SCHEMA10_ACCEPTED_RESPONSE_ASSEMBLY_TERMINAL_RECEIPT.json"
    write_exclusive(terminal_path, terminal)
    print(json.dumps({"status": "PASS", "terminal": terminal_path.as_posix(), "terminal_sha256": sha(terminal_path)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
