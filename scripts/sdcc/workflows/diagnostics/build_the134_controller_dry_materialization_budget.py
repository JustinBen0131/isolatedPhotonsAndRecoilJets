#!/usr/bin/env python3
"""Derive the THE-134 controller budget from the exact dry serializers.

This helper is deliberately non-submitting.  It reopens the pinned resolver
plan, resolver receipt, immutable bundle, materialization receipt, source
authority, and execution partition.  It then calls the production
``staged_row`` and ``staged_job`` serializers for every resolved row and every
one of the 18,577 partition records.

The strict budget receipt intentionally retains the exact schema consumed by
``project_the134_preextraction_storage_quota.py``.  A separately hash-linked
derivation receipt records the measurements and the manifest-envelope fixed
point without weakening that consumer contract.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import sys
from pathlib import Path, PurePosixPath
from typing import Any, Iterable, Mapping


HERE = Path(__file__).resolve().parent
MATERIALIZER_PATH = HERE / "materialize_the134_full_multiview_extraction.py"


def _load_local_module(name: str, path: Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load local module: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


materializer = _load_local_module(
    "the134_controller_budget_materializer", MATERIALIZER_PATH
)
resolver = materializer.resolver
projector = materializer.projector

BUDGET_SCHEMA = projector.CONTROLLER_BUDGET_SCHEMA
DERIVATION_SCHEMA = (
    "THE134_FULL_EXTRACTION_CONTROLLER_DRY_MATERIALIZATION_BUDGET_DERIVATION_V1"
)
PASS_STATUS = "PASS"
DERIVATION_PASS_STATUS = "PASS_EXACT_SERIALIZER_MEASUREMENT"
MAX_ARTIFACT_SIZE_BYTES = 9_999_999_999_999_999_999
MAX_FIXED_POINT_ITERATIONS = 32


class BudgetError(RuntimeError):
    """Malformed, inconsistent, or non-authoritative budget input."""


def canonical_json_bytes(payload: Any) -> bytes:
    return materializer.canonical_json_bytes(payload)


def file_sha256(path: Path) -> str:
    return materializer.file_sha256(path)


def artifact_record(path: Path) -> dict[str, Any]:
    absolute = path.absolute()
    if not absolute.is_file() or absolute.is_symlink():
        raise BudgetError(f"artifact must be a regular non-symlink file: {absolute}")
    return {
        "path": str(absolute),
        "sha256": file_sha256(absolute),
        "size_bytes": absolute.stat().st_size,
    }


def require_mapping(value: object, label: str) -> dict[str, Any]:
    try:
        return materializer.require_mapping(value, label)
    except materializer.ControllerError as exc:
        raise BudgetError(str(exc)) from exc


def require_sequence(value: object, label: str) -> list[Any]:
    try:
        return materializer.require_sequence(value, label)
    except materializer.ControllerError as exc:
        raise BudgetError(str(exc)) from exc


def require_sha256(value: object, label: str) -> str:
    try:
        return materializer.require_sha256(value, label)
    except materializer.ControllerError as exc:
        raise BudgetError(str(exc)) from exc


def load_pinned_json(
    path: Path, expected_sha256: str, label: str
) -> tuple[dict[str, Any], dict[str, Any]]:
    absolute = path.absolute()
    if absolute.is_symlink():
        raise BudgetError(f"{label} must not be a symlink: {absolute}")
    try:
        return materializer.load_pinned_json(absolute, expected_sha256, label)
    except materializer.ControllerError as exc:
        raise BudgetError(str(exc)) from exc


def validate_plan_surface(plan: Mapping[str, Any]) -> dict[str, str]:
    expected_keys = {
        "schema",
        "status",
        "execution_state",
        "submission_performed",
        "campaign",
        "training_period_si_contract",
        "training_period_si_contract_sha256",
        "authority",
        "artifact_profile",
        "closure_witness_boundary",
        "source_family_closure",
        "execution_partition",
        "input_manifests",
        "duplicate_contract",
        "duplicate_fingerprint_sha256",
        "execution_fingerprint_sha256",
        "rows",
    }
    try:
        materializer.require_exact_keys(plan, expected_keys, "resolved extraction plan")
        materializer.validate_authority_occurrences(plan, "resolved extraction plan")
        resolver.validate_preflight_authority_payload(plan, label="plan")
    except (materializer.ControllerError, resolver.ControllerError) as exc:
        raise BudgetError(str(exc)) from exc
    if (
        plan.get("schema") != resolver.PLAN_SCHEMA
        or plan.get("status") != "PREFLIGHT_PASS"
        or plan.get("execution_state") != "PREFLIGHT_ONLY_NO_CONDOR_MUTATION"
        or plan.get("submission_performed") is not False
        or plan.get("artifact_profile") != materializer.RESOLVED_ARTIFACT_PROFILE
    ):
        raise BudgetError("resolved plan schema/status/authority differs")
    authority = require_mapping(plan.get("authority"), "plan authority")
    if authority != {
        "requested_scope": "full",
        "full_training_authority": 0,
        "authority_state": resolver.PREFLIGHT_AUTHORITY_STATE,
    }:
        raise BudgetError("resolved plan authority differs")
    campaign_raw = require_mapping(plan.get("campaign"), "plan campaign")
    try:
        materializer.require_exact_keys(
            campaign_raw,
            {"tag", "output_root", "evidence_root", "submit_root"},
            "plan campaign",
        )
    except materializer.ControllerError as exc:
        raise BudgetError(str(exc)) from exc
    tag = campaign_raw.get("tag")
    if not isinstance(tag, str) or resolver.TAG_RE.fullmatch(tag) is None:
        raise BudgetError("campaign tag differs")
    campaign = {"tag": tag}
    for field in ("output_root", "evidence_root", "submit_root"):
        try:
            value = materializer.validate_remote_root(
                campaign_raw.get(field), f"campaign.{field}"
            )
        except materializer.ControllerError as exc:
            raise BudgetError(str(exc)) from exc
        if PurePosixPath(value).name != tag:
            raise BudgetError(f"campaign {field} basename differs from tag")
        campaign[field] = value
    if len(set(campaign.values())) != 4:
        raise BudgetError("campaign tag/namespaces must be distinct")
    return campaign


def validate_receipt_surface(
    receipt: Mapping[str, Any],
    *,
    plan_artifact: Mapping[str, Any],
) -> dict[str, Any]:
    if (
        receipt.get("schema") != resolver.RECEIPT_SCHEMA
        or receipt.get("status") != "PASS"
        or receipt.get("submission_performed") is not False
        or receipt.get("row_count") != materializer.EXPECTED_ROW_COUNT
        or receipt.get("requested_scope") != "full"
        or receipt.get("full_training_authority") != 0
        or receipt.get("authority_state") != resolver.PREFLIGHT_AUTHORITY_STATE
        or receipt.get("artifact_profile_sha256")
        != materializer.RESOLVED_ARTIFACT_PROFILE_SHA256
    ):
        raise BudgetError("resolver preflight receipt differs")
    try:
        materializer.validate_authority_occurrences(
            receipt, "resolver preflight receipt"
        )
    except materializer.ControllerError as exc:
        raise BudgetError(str(exc)) from exc
    artifacts = require_mapping(receipt.get("artifacts"), "preflight artifacts")
    plan_ref = require_mapping(artifacts.get("plan"), "preflight plan binding")
    if plan_ref.get("sha256") != plan_artifact["sha256"]:
        raise BudgetError("preflight receipt binds a different plan")
    return artifacts


def validate_input_chain(
    *,
    plan: Mapping[str, Any],
    plan_path: Path,
    plan_artifact: Mapping[str, Any],
    receipt: Mapping[str, Any],
    receipt_artifact: Mapping[str, Any],
) -> dict[str, Any]:
    campaign = validate_plan_surface(plan)
    receipt_artifacts = validate_receipt_surface(
        receipt, plan_artifact=plan_artifact
    )
    inputs = require_mapping(plan.get("input_manifests"), "plan input manifests")
    bundle_ref = require_mapping(inputs.get("bundle"), "plan bundle ref")
    materialization_ref = require_mapping(
        inputs.get("materialization"), "plan materialization ref"
    )
    source_ref = require_mapping(inputs.get("sources"), "plan source ref")
    bundle_path = Path(str(bundle_ref.get("path", ""))).absolute()
    materialization_path = Path(str(materialization_ref.get("path", ""))).absolute()
    source_path = Path(str(source_ref.get("path", ""))).absolute()
    bundle_payload, bundle_artifact = load_pinned_json(
        bundle_path, str(bundle_ref.get("sha256", "")), "immutable resolver bundle"
    )
    materialization_payload, materialization_artifact = load_pinned_json(
        materialization_path,
        str(materialization_ref.get("sha256", "")),
        "immutable materialization receipt",
    )
    source_payload, source_artifact = load_pinned_json(
        source_path,
        str(source_ref.get("sha256", "")),
        "source authority manifest",
    )
    try:
        bundle, by_role = materializer.validate_bundle(bundle_payload)
        materializer.validate_materialization_binding(
            materialization_payload,
            materialization_path=materialization_path,
            bundle_path=bundle_path,
            bundle_sha256=bundle_artifact["sha256"],
            bundle=bundle,
            plan_record=materialization_ref,
        )
    except materializer.ControllerError as exc:
        raise BudgetError(str(exc)) from exc
    if (
        source_payload.get("schema") != resolver.SOURCE_SCHEMA
        or source_payload.get("status") != "PASS"
    ):
        raise BudgetError("source authority schema/status differs")
    training_contract = require_mapping(
        plan.get("training_period_si_contract"), "training period/SI contract"
    )
    training_period = str(training_contract.get("period", ""))
    try:
        source_authority = resolver.validate_source_period_authority(
            source_payload, pp_period=training_period
        )
    except resolver.ControllerError as exc:
        raise BudgetError(str(exc)) from exc
    if source_authority != {
        "pp_period": training_period,
        "pp_si_di_role": training_contract.get("si_di_role"),
    }:
        raise BudgetError("source period/SI authority differs")
    if (
        receipt.get("bundle_manifest_sha256") != bundle_artifact["sha256"]
        or receipt.get("materialization_receipt_sha256")
        != materialization_artifact["sha256"]
        or receipt.get("source_manifest_sha256") != source_artifact["sha256"]
    ):
        raise BudgetError("preflight bundle/materialization/source binding differs")
    try:
        rows, rows_by_id = materializer.validate_rows(
            plan,
            campaign=campaign,
            bundle=bundle,
            by_role=by_role,
            source_manifest=source_payload,
        )
    except materializer.ControllerError as exc:
        raise BudgetError(str(exc)) from exc
    partition_ref = require_mapping(
        receipt_artifacts.get("partition"), "preflight partition binding"
    )
    partition_name = str(partition_ref.get("name", ""))
    plan_partition = require_mapping(
        require_mapping(plan.get("execution_partition"), "execution partition").get(
            "partition_artifact"
        ),
        "plan partition artifact",
    )
    if (
        not partition_name
        or Path(partition_name).name != partition_name
        or partition_name != plan_partition.get("name")
    ):
        raise BudgetError("partition artifact name differs")
    partition_path = plan_path.absolute().parent / partition_name
    partition_artifact = artifact_record(partition_path)
    if (
        partition_ref.get("sha256") != partition_artifact["sha256"]
        or partition_ref.get("record_count") != materializer.EXPECTED_JOB_COUNT
        or receipt.get("execution_partition_sha256")
        != materializer.canonical_sha256(plan["execution_partition"])
    ):
        raise BudgetError("preflight partition receipt differs")
    try:
        chunks, _ = materializer.validate_partition(
            partition_path, plan=plan, rows_by_id=rows_by_id
        )
    except materializer.ControllerError as exc:
        raise BudgetError(str(exc)) from exc
    duplicate = require_mapping(plan.get("duplicate_contract"), "duplicate contract")
    duplicate_sha = require_sha256(
        plan.get("duplicate_fingerprint_sha256"), "duplicate fingerprint"
    )
    if (
        duplicate.get("schema") != resolver.DUPLICATE_SCHEMA
        or duplicate.get("requested_scope") != "full"
        or duplicate.get("full_training_authority") != 0
        or duplicate.get("authority_state") != resolver.PREFLIGHT_AUTHORITY_STATE
        or duplicate_sha != materializer.canonical_sha256(duplicate)
        or receipt.get("duplicate_fingerprint_sha256") != duplicate_sha
    ):
        raise BudgetError("duplicate fingerprint differs")
    duplicate_sources = require_sequence(
        duplicate.get("sources"), "duplicate contract sources"
    )
    duplicate_ids = [
        record.get("row_id") if isinstance(record, dict) else None
        for record in duplicate_sources
    ]
    if (
        duplicate_ids != list(rows_by_id)
        or len(set(duplicate_ids)) != materializer.EXPECTED_ROW_COUNT
    ):
        raise BudgetError("duplicate contract rows are missing or duplicated")
    execution_sha = require_sha256(
        plan.get("execution_fingerprint_sha256"), "execution fingerprint"
    )
    expected_execution_sha = materializer.canonical_sha256(
        {
            "schema": resolver.EXECUTION_SCHEMA,
            "tag": campaign["tag"],
            "output_root": campaign["output_root"],
            "evidence_root": campaign["evidence_root"],
            "submit_root": campaign["submit_root"],
            "materialization_receipt_sha256": materialization_artifact["sha256"],
            "bundle_manifest_sha256": bundle_artifact["sha256"],
            "source_manifest_sha256": source_artifact["sha256"],
            "duplicate_fingerprint_sha256": duplicate_sha,
            "partition_artifact_sha256": partition_artifact["sha256"],
            "execution_partition_sha256": materializer.canonical_sha256(
                plan["execution_partition"]
            ),
            "row_fingerprints": [
                row["row_fingerprint_sha256"] for row in rows
            ],
        }
    )
    if (
        execution_sha != expected_execution_sha
        or receipt.get("execution_fingerprint_sha256") != execution_sha
    ):
        raise BudgetError("execution fingerprint receipt binding differs")
    return {
        "plan": plan,
        "plan_artifact": plan_artifact,
        "preflight_receipt_artifact": receipt_artifact,
        "bundle_artifact": bundle_artifact,
        "materialization_artifact": materialization_artifact,
        "source_artifact": source_artifact,
        "partition_artifact": partition_artifact,
        "campaign": campaign,
        "rows": rows,
        "rows_by_id": rows_by_id,
        "chunks": chunks,
        "duplicate_fingerprint_sha256": duplicate_sha,
        "execution_fingerprint_sha256": execution_sha,
    }


def envelope_artifact(path: Path) -> dict[str, Any]:
    """Return a maximum-width artifact record for future exact-path evidence."""

    return {
        "path": str(path.absolute()),
        "sha256": "f" * 64,
        "size_bytes": MAX_ARTIFACT_SIZE_BYTES,
    }


def measure_staged_job_records(
    chunks: Sequence[Mapping[str, Any]],
    rows_by_id: Mapping[str, Mapping[str, Any]],
) -> dict[str, int]:
    """Measure exact staged-job serialization without retaining every record."""

    record_count = 0
    exact_bytes = 0
    min_bytes: int | None = None
    max_bytes = 0
    for chunk in chunks:
        record_bytes = canonical_json_bytes(
            materializer.staged_job(chunk, rows_by_id[chunk["row_id"]])
        )
        size_bytes = len(record_bytes)
        record_count += 1
        exact_bytes += size_bytes
        min_bytes = (
            size_bytes if min_bytes is None else min(min_bytes, size_bytes)
        )
        max_bytes = max(max_bytes, size_bytes)
    if record_count == 0 or min_bytes is None or max_bytes <= 0:
        raise BudgetError("serialized per-job measurements are empty or invalid")
    return {
        "record_count": record_count,
        "exact_bytes": exact_bytes,
        "min_bytes": min_bytes,
        "max_bytes": max_bytes,
    }


def derive_budget(
    context: Mapping[str, Any],
    *,
    budget_output: Path,
    future_storage_certificate: Path,
) -> tuple[dict[str, Any], dict[str, Any]]:
    row_records = [
        materializer.staged_row(row) for row in context["rows"]
    ]
    job_measurements = measure_staged_job_records(
        context["chunks"], context["rows_by_id"]
    )
    if (
        len(row_records) != materializer.EXPECTED_ROW_COUNT
        or job_measurements["record_count"] != materializer.EXPECTED_JOB_COUNT
    ):
        raise BudgetError("exact serializer record counts differ")
    row_bytes = materializer.jsonl_bytes(row_records)
    description_bytes = materializer.submit_description_bytes(context)
    bytes_per_job = job_measurements["max_bytes"]

    fixed_inodes = len(materializer.OUTPUT_FILENAMES)
    inodes_per_job = 1
    fixed_bytes = len(row_bytes) + len(description_bytes)
    manifest_bytes = b""
    rendered: dict[str, bytes] = {}
    iterations = 0
    for iterations in range(1, MAX_FIXED_POINT_ITERATIONS + 1):
        projected_bytes = (
            fixed_bytes
            + materializer.EXPECTED_JOB_COUNT * bytes_per_job
        )
        projected_inodes = (
            fixed_inodes
            + materializer.EXPECTED_JOB_COUNT * inodes_per_job
        )
        envelope_context = {
            **context,
            "storage_certificate_artifact": envelope_artifact(
                future_storage_certificate
            ),
            "controller_budget_artifact": envelope_artifact(budget_output),
            "budget": {
                "projected_bytes": projected_bytes,
                "projected_inodes": projected_inodes,
            },
        }
        try:
            rendered = materializer.build_staged_artifacts(envelope_context)
        except materializer.ControllerError as exc:
            raise BudgetError(str(exc)) from exc
        manifest_bytes = rendered[materializer.OUTPUT_FILENAMES[3]]
        new_fixed = len(row_bytes) + len(description_bytes) + len(manifest_bytes)
        if new_fixed == fixed_bytes:
            break
        fixed_bytes = new_fixed
    else:
        raise BudgetError("budget manifest-envelope fixed point did not converge")

    projected_bytes = fixed_bytes + materializer.EXPECTED_JOB_COUNT * bytes_per_job
    projected_inodes = (
        fixed_inodes + materializer.EXPECTED_JOB_COUNT * inodes_per_job
    )
    actual_envelope_bytes = sum(len(payload) for payload in rendered.values())
    actual_envelope_inodes = len(rendered)
    if (
        job_measurements["exact_bytes"]
        > materializer.EXPECTED_JOB_COUNT * bytes_per_job
        or actual_envelope_bytes > projected_bytes
        or actual_envelope_inodes > projected_inodes
    ):
        raise BudgetError("derived controller budget does not cover exact serializer output")

    budget = {
        "schema": BUDGET_SCHEMA,
        "status": PASS_STATUS,
        "submission_performed": False,
        "expected_job_count": materializer.EXPECTED_JOB_COUNT,
        "storage_budget": {
            "fixed_bytes": fixed_bytes,
            "fixed_inodes": fixed_inodes,
            "bytes_per_job": bytes_per_job,
            "inodes_per_job": inodes_per_job,
        },
        "full_training_authority": 0,
        "full_extraction_authority": False,
    }
    budget_bytes = canonical_json_bytes(budget)
    budget_sha = hashlib.sha256(budget_bytes).hexdigest()
    derivation = {
        "schema": DERIVATION_SCHEMA,
        "status": DERIVATION_PASS_STATUS,
        "submission_performed": False,
        "authority": {
            "state": "NON_SUBMITTING_CONTROLLER_STORAGE_INPUT_ONLY",
            "full_training_authority": 0,
            "full_extraction_authority": False,
            "science_freeze_authority": False,
            "broad_production_authority": False,
            "canonical_promotion": False,
        },
        "serializer": {
            "path": str(MATERIALIZER_PATH),
            "sha256": file_sha256(MATERIALIZER_PATH),
            "row_function": "staged_row",
            "job_function": "staged_job",
            "manifest_function": "build_staged_artifacts",
        },
        "inputs": {
            "plan": context["plan_artifact"],
            "preflight_receipt": context["preflight_receipt_artifact"],
            "immutable_bundle": context["bundle_artifact"],
            "immutable_materialization": context["materialization_artifact"],
            "source_authority": context["source_artifact"],
            "source_partition": context["partition_artifact"],
            "duplicate_fingerprint_sha256": context[
                "duplicate_fingerprint_sha256"
            ],
            "execution_fingerprint_sha256": context[
                "execution_fingerprint_sha256"
            ],
        },
        "measurements": {
            "row_record_count": len(row_records),
            "row_manifest_bytes": len(row_bytes),
            "job_record_count": job_measurements["record_count"],
            "job_manifest_exact_bytes": job_measurements["exact_bytes"],
            "job_record_min_bytes": job_measurements["min_bytes"],
            "job_record_max_bytes": bytes_per_job,
            "job_record_ceiling_bytes": (
                materializer.EXPECTED_JOB_COUNT * bytes_per_job
            ),
            "submit_description_bytes": len(description_bytes),
            "manifest_envelope_bytes": len(manifest_bytes),
            "exact_envelope_total_bytes": actual_envelope_bytes,
            "exact_envelope_total_inodes": actual_envelope_inodes,
            "fixed_point_iterations": iterations,
        },
        "budget_derivation": {
            "fixed_bytes": {
                "value": fixed_bytes,
                "basis": (
                    "exact row manifest + exact submit description + "
                    "maximum-width future-binding materialization manifest"
                ),
            },
            "fixed_inodes": {
                "value": fixed_inodes,
                "basis": "exact dry-materialization fixed artifact inventory",
            },
            "bytes_per_job": {
                "value": bytes_per_job,
                "basis": "maximum exact staged_job canonical JSON record bytes",
            },
            "inodes_per_job": {
                "value": inodes_per_job,
                "basis": "conservative one logical controller record per job",
            },
            "expected_job_count": {
                "value": materializer.EXPECTED_JOB_COUNT,
                "basis": "validated resolver partition record count",
            },
            "projected_bytes": projected_bytes,
            "projected_inodes": projected_inodes,
        },
        "future_bindings": {
            "budget_output_path": str(budget_output.absolute()),
            "storage_certificate_path": str(
                future_storage_certificate.absolute()
            ),
            "artifact_size_width_ceiling": MAX_ARTIFACT_SIZE_BYTES,
        },
        "budget_receipt": {
            "path": str(budget_output.absolute()),
            "sha256": budget_sha,
            "size_bytes": len(budget_bytes),
        },
        "boundaries": [
            "No Condor submission or job control.",
            "No output, evidence, or submit namespace creation.",
            "No scientific, model, working-point, production, or CANONICAL authority.",
        ],
    }
    derivation["derivation_semantic_sha256"] = projector.semantic_sha256(
        derivation
    )
    return budget, derivation


def write_fresh(path: Path, payload: Mapping[str, Any], label: str) -> None:
    absolute = path.absolute()
    if absolute.exists() or absolute.is_symlink():
        raise BudgetError(f"{label} output must be fresh: {absolute}")
    if not absolute.parent.is_dir():
        raise BudgetError(f"{label} output parent is missing: {absolute.parent}")
    try:
        absolute.write_bytes(canonical_json_bytes(payload))
    except OSError as exc:
        raise BudgetError(f"cannot write {label}: {absolute}") from exc


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--plan-sha256", required=True)
    parser.add_argument("--preflight-receipt", type=Path, required=True)
    parser.add_argument("--preflight-receipt-sha256", required=True)
    parser.add_argument("--future-storage-certificate", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--derivation-output", type=Path, required=True)
    return parser.parse_args(argv)


def main(argv: Iterable[str] | None = None) -> int:
    try:
        args = parse_args(argv)
        output = args.output.absolute()
        derivation_output = args.derivation_output.absolute()
        future_storage = args.future_storage_certificate.absolute()
        if len({output, derivation_output, future_storage}) != 3:
            raise BudgetError("budget, derivation, and future storage paths must differ")
        for path, label in (
            (output, "budget"),
            (derivation_output, "derivation"),
        ):
            if path.exists() or path.is_symlink():
                raise BudgetError(f"{label} output must be fresh: {path}")
        plan, plan_artifact = load_pinned_json(
            args.plan, args.plan_sha256, "resolved extraction plan"
        )
        receipt, receipt_artifact = load_pinned_json(
            args.preflight_receipt,
            args.preflight_receipt_sha256,
            "resolver preflight receipt",
        )
        context = validate_input_chain(
            plan=plan,
            plan_path=args.plan,
            plan_artifact=plan_artifact,
            receipt=receipt,
            receipt_artifact=receipt_artifact,
        )
        budget, derivation = derive_budget(
            context,
            budget_output=output,
            future_storage_certificate=future_storage,
        )
        # Write the non-authoritative derivation first so a partial failure can
        # never leave a consumable strict budget without its audit receipt.
        write_fresh(derivation_output, derivation, "derivation")
        write_fresh(output, budget, "budget")
        budget_observed = artifact_record(output)
        if budget_observed != derivation["budget_receipt"]:
            raise BudgetError("written budget readback differs from derivation")
        result = {
            "schema": DERIVATION_SCHEMA,
            "status": DERIVATION_PASS_STATUS,
            "budget_receipt": budget_observed,
            "derivation_receipt": artifact_record(derivation_output),
            "expected_job_count": materializer.EXPECTED_JOB_COUNT,
            "projected_bytes": (
                budget["storage_budget"]["fixed_bytes"]
                + materializer.EXPECTED_JOB_COUNT
                * budget["storage_budget"]["bytes_per_job"]
            ),
            "projected_inodes": (
                budget["storage_budget"]["fixed_inodes"]
                + materializer.EXPECTED_JOB_COUNT
                * budget["storage_budget"]["inodes_per_job"]
            ),
            "submission_performed": False,
        }
        print(json.dumps(result, sort_keys=True))
        return 0
    except (BudgetError, FileNotFoundError) as exc:
        print(f"[THE134-CONTROLLER-BUDGET][ERROR] {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
