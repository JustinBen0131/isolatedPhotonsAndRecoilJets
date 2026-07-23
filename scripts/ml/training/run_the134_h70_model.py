#!/usr/bin/env python3
"""Plan, reuse, or run one frozen THE-134 factorial-view model job.

This is a focused local/worker entrypoint, not a Condor submitter.  It reuses
the established trainer with exactly one model spec and writes a separate
view-qualified metadata receipt.  The audited extraction matrix is never
modified; training operates on a private working copy because the established
trainer adds the global PPG12-exact weight column to its cache.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve()
_CONTRACTS = _HERE.parents[1] / "contracts"
sys.path.insert(0, str(_CONTRACTS))
from the134_h70_contract import (  # noqa: E402
    ALL_SHOWER_VIEWS,
    FEATURES_BY_SYSTEM,
    MODEL_DOMAIN_GEV,
    READER_ET_EDGES,
    VIEW_NAME,
    expected_model_origin,
    expected_model_split,
    expected_reuse_pinned_hashes,
    expected_sources,
    extraction_authority_binding,
    full_extraction_authority_checks,
    sha256_file,
    shower_semantic_sha256,
    valid_reuse_pinned_hashes,
)


SYSTEM_DEFAULTS = {
    "pp": {
        "campaign": "ppg12-sixpack",
        "model_id": "ppg12_base_v3E_bdt_noIso",
        "seed": 42,
        "test_size": 0.5,
        "row_cap_per_class": 2_000_000,
        "output_stem": "pp_tight_bdt_ppg12_base_v3E_bdt_noIso_tmva",
    },
    "auau": {
        "campaign": "corrected-baseline-shower-ladder",
        "model_id": "centAsFeatBase3x3_pt15to35",
        "seed": 13,
        "test_size": 0.1,
        "row_cap_per_class": 0,
        "output_stem": "auau_tight_bdt_centAsFeatBase3x3_pt15to35_tmva",
    },
}
AUTHORIZED_REUSE = {("pp", "H70"): "REUSED_THE116", ("auau", "H0"): "REUSED_THE111"}


def verify_reuse_audit_artifacts(reuse: dict, system: str, view_name: str) -> dict[str, bool]:
    """Recompute every artifact named by a reuse audit at consumption time."""

    artifacts = reuse.get("artifacts", {})
    artifact_pairs = {
        "xgboost": "xgboost_sha256",
        "tmva": "tmva_sha256",
        "metadata": "metadata_sha256",
        "authority_certificate": "authority_certificate_sha256",
    }
    checks = {}
    for path_key, sha_key in artifact_pairs.items():
        path = Path(str(artifacts.get(path_key, "")))
        checks[f"artifact_{path_key}"] = (
            path.is_file()
            and path.stat().st_size > 0
            and sha256_file(path) == artifacts.get(sha_key)
        )
    expected_pins = expected_reuse_pinned_hashes(system, view_name)
    checks["frozen_pins"] = valid_reuse_pinned_hashes(
        reuse.get("pinned_hashes"), system, view_name
    )
    pinned_artifacts = reuse.get("pinned_artifacts", {})
    if (system, view_name) == ("pp", "H70"):
        for key, expected_hash in expected_pins.items():
            path = Path(str(pinned_artifacts.get(key, "")))
            checks[f"pinned_artifact_{key}"] = (
                path.is_file() and sha256_file(path) == expected_hash
            )
    else:
        mapping = {
            "authority_certificate": "authority_certificate_sha256",
            "model_xgb": "model_xgb_sha256",
            "model_tmva": "model_tmva_sha256",
            "model_metadata": "model_metadata_sha256",
        }
        for path_key, pin_key in mapping.items():
            path = Path(str(pinned_artifacts.get(path_key, "")))
            checks[f"pinned_artifact_{path_key}"] = (
                path.is_file() and sha256_file(path) == expected_pins[pin_key]
            )
    checks["authority_gates"] = bool(reuse.get("authority_checks")) and all(
        reuse.get("authority_checks", {}).values()
    )
    checks["metadata_gates"] = bool(reuse.get("metadata_checks")) and all(
        reuse.get("metadata_checks", {}).values()
    )
    return checks


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--system", choices=sorted(SYSTEM_DEFAULTS), required=True)
    parser.add_argument("--view", choices=ALL_SHOWER_VIEWS, default=VIEW_NAME)
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--extraction-audit", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument(
        "--trainer",
        type=Path,
        default=_HERE.with_name("train_auau_photon_bdt.py"),
    )
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument(
        "--reuse-audit",
        type=Path,
        help="PASS certificate from audit_the134_model_reuse.py",
    )
    parser.add_argument("--execute", action="store_true")
    return parser.parse_args()


def deterministic_pp_cap(source: Path, destination: Path, cap: int, seed: int) -> dict:
    """Clone the established file-shuffle/within-file row-cap convention."""

    data = np.load(source, allow_pickle=True)
    columns = [str(item) for item in data["__columns__"].tolist()]
    arrays = {column: np.asarray(data[column]) for column in columns}
    labels = np.asarray(arrays["is_signal"], dtype=np.int8)
    file_index = np.asarray(arrays["input_file_index"], dtype=np.int64)
    rng = np.random.default_rng(seed)
    files = np.unique(file_index)
    rng.shuffle(files)
    # Preserve the established loader's exact row order: shuffled file order,
    # then class 0 and class 1 rows for each file, with selected in-file entry
    # indices sorted.  Sorting the final union would silently change XGBoost's
    # row order even though the selected set was identical.
    keep: list[np.ndarray] = []
    counts = {0: 0, 1: 0}
    for file_id in files:
        local_file = file_index == file_id
        local_keep = []
        for label in (0, 1):
            indices = np.flatnonzero(local_file & (labels == label))
            remaining = max(0, cap - counts[label])
            n_take = min(len(indices), remaining)
            if n_take <= 0:
                continue
            if n_take < len(indices):
                indices = rng.choice(indices, size=n_take, replace=False)
            local_keep.append(np.sort(indices))
            counts[label] += n_take
        if local_keep:
            keep.append(np.concatenate(local_keep))
        if all(counts[label] >= cap for label in (0, 1)):
            break
    selected = np.concatenate(keep) if keep else np.asarray([], dtype=np.int64)
    if not len(selected):
        raise SystemExit("deterministic pp row cap selected zero rows")
    payload = {column: values[selected] for column, values in arrays.items()}
    payload["__columns__"] = np.asarray(columns, dtype=object)
    destination.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(destination, **payload)
    return {
        "mode": "clone-established-file-shuffle-within-file-row-cap",
        "seed": seed,
        "cap_per_class": cap,
        "class_counts": {str(label): counts[label] for label in (0, 1)},
        "rows_before": int(len(labels)),
        "rows_after": int(len(selected)),
        "source_sha256": sha256_file(source),
        "output_sha256": sha256_file(destination),
    }


def enrich_weighted_cache(weighted: Path, snapshot: Path, output: Path) -> dict:
    """Reattach audited row provenance dropped by the generic trainer cache rewrite."""

    weighted_data = np.load(weighted, allow_pickle=True)
    snapshot_data = np.load(snapshot, allow_pickle=True)
    weighted_columns = [str(item) for item in weighted_data["__columns__"].tolist()]
    snapshot_columns = [str(item) for item in snapshot_data["__columns__"].tolist()]
    key_columns = ("source_sample", "input_file_index", "input_tree_entry")
    for name in key_columns:
        if name not in weighted_data.files or name not in snapshot_data.files:
            raise SystemExit(f"cannot enrich weighted cache without {name}")
    weighted_keys = list(
        zip(*(np.asarray(weighted_data[name]).astype(str).tolist() for name in key_columns))
    )
    snapshot_keys = list(
        zip(*(np.asarray(snapshot_data[name]).astype(str).tolist() for name in key_columns))
    )
    if len(snapshot_keys) != len(set(snapshot_keys)):
        raise SystemExit("preweight snapshot has duplicate source/file/entry row keys")
    index = {key: position for position, key in enumerate(snapshot_keys)}
    if any(key not in index for key in weighted_keys):
        raise SystemExit("weighted cache contains rows absent from preweight snapshot")
    aligned = np.asarray([index[key] for key in weighted_keys], dtype=np.int64)
    payload = {name: np.asarray(weighted_data[name]) for name in weighted_columns}
    for name in snapshot_columns:
        if name not in payload:
            payload[name] = np.asarray(snapshot_data[name])[aligned]
    payload["__columns__"] = np.asarray(sorted(payload), dtype=object)
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output, **payload)
    return {
        "schema": "THE134_FACTORIAL_VIEW_ENRICHED_WEIGHTED_CACHE_V1",
        "rows": len(weighted_keys),
        "training_cache": str(weighted),
        "training_cache_sha256": sha256_file(weighted),
        "preweight_snapshot": str(snapshot),
        "preweight_snapshot_sha256": sha256_file(snapshot),
        "output": str(output),
        "output_sha256": sha256_file(output),
        "reattached_columns": sorted(set(snapshot_columns) - set(weighted_columns)),
    }


def validate_reuse_weighted_matrix(path: Path) -> dict:
    """Require a reusable model to hand downstream an exact weighted matrix."""

    data = np.load(path, allow_pickle=True)
    if "__ppg12_exact_training_weight" not in data.files:
        raise SystemExit("reuse weight materialization omitted PPG12-exact weights")
    weights = np.asarray(data["__ppg12_exact_training_weight"], dtype=np.float64)
    if not len(weights) or not np.all(np.isfinite(weights) & (weights > 0.0)):
        raise SystemExit("reuse weight materialization produced invalid weights")
    return {
        "rows": int(len(weights)),
        "minimum_weight": float(np.min(weights)),
        "maximum_weight": float(np.max(weights)),
        "weighted_matrix_sha256": sha256_file(path),
    }


def build_training_completion_receipt(
    *,
    system: str,
    view_name: str,
    model_id: str,
    extraction_binding: dict,
    artifacts: dict,
    enrichment: dict,
) -> dict:
    """Build the immutable receipt consumed by every downstream model gate."""

    return {
        "schema": "THE134_FACTORIAL_VIEW_MODEL_TRAINING_COMPLETE_V1",
        "status": "PASS",
        "system": system,
        "model_id": model_id,
        "shower_definition": view_name,
        "shower_semantic_sha256": shower_semantic_sha256(view_name),
        "model_origin": "TRAINED_THE134",
        "extraction_authority": extraction_binding,
        "artifacts": artifacts,
        "weighted_cache_enrichment": enrichment,
    }


def main() -> int:
    args = parse_args()
    defaults = SYSTEM_DEFAULTS[args.system]
    for path in (args.matrix, args.extraction_audit, args.trainer):
        if not path.is_file():
            raise SystemExit(f"missing required input: {path}")
    extraction = json.loads(args.extraction_audit.read_text())
    extraction_checks = full_extraction_authority_checks(
        extraction,
        args.system,
        args.view,
        matrix_sha256=sha256_file(args.matrix),
    )
    if not all(extraction_checks.values()):
        raise SystemExit(
            "full extraction authority mismatch: "
            + json.dumps(extraction_checks, sort_keys=True)
        )
    extraction_binding = extraction_authority_binding(extraction)

    if args.reuse_audit is not None:
        if (args.system, args.view) not in AUTHORIZED_REUSE:
            raise SystemExit(
                f"reuse is not authorized for {args.system}/{args.view}; frozen contract "
                "requires a new matched model"
            )
        if not args.reuse_audit.is_file():
            raise SystemExit(f"missing reuse audit: {args.reuse_audit}")
        reuse = json.loads(args.reuse_audit.read_text())
        reuse_artifact_checks = verify_reuse_audit_artifacts(
            reuse, args.system, args.view
        )
        reuse_checks = {
            "schema": reuse.get("schema")
            == "THE134_FACTORIAL_VIEW_MODEL_REUSE_AUDIT_V1",
            "status": reuse.get("status") == "PASS",
            "system": reuse.get("system") == args.system,
            "view": reuse.get("shower_definition") == args.view,
            "semantic": reuse.get("shower_semantic_sha256")
            == shower_semantic_sha256(args.view),
            "feature_order": reuse.get("feature_order")
            == list(FEATURES_BY_SYSTEM[args.system]),
            "model_origin": reuse.get("model_origin")
            == AUTHORIZED_REUSE[(args.system, args.view)],
            "frozen_pins_exact": valid_reuse_pinned_hashes(
                reuse.get("pinned_hashes"), args.system, args.view
            ),
            "all_reuse_artifacts_reverified": all(
                reuse_artifact_checks.values()
            ),
        }
        if not all(reuse_checks.values()):
            raise SystemExit("reuse audit mismatch: " + json.dumps(reuse_checks, sort_keys=True))
        args.outdir.mkdir(parents=True, exist_ok=True)
        artifacts = reuse.get("artifacts", {})
        original_metadata_path = Path(str(artifacts.get("metadata", "")))
        if not original_metadata_path.is_file():
            raise SystemExit("reuse audit does not reference readable model metadata")
        reused_metadata = json.loads(original_metadata_path.read_text())
        reused_metadata["model_origin"] = reuse.get("model_origin")
        reused_metadata["reuse_pinned_hashes"] = reuse.get("pinned_hashes")
        reused_metadata["accepted_split_contract"] = reuse.get(
            "accepted_split_contract"
        )
        reused_metadata["the134_view_contract"] = {
            "schema": "THE134_FACTORIAL_VIEW_MODEL_VIEW_CONTRACT_V1",
            "system": args.system,
            "shower_definition": args.view,
            "shower_semantic_sha256": shower_semantic_sha256(args.view),
            "feature_order": list(FEATURES_BY_SYSTEM[args.system]),
            "model_domain_gev": list(MODEL_DOMAIN_GEV),
            "model_origin": reuse.get("model_origin"),
            "label_authority": "frozen source-role training_label materialized without relabeling",
            "extraction_audit": str(args.extraction_audit),
            "extraction_audit_sha256": sha256_file(args.extraction_audit),
            "reuse_audit": str(args.reuse_audit),
            "reuse_audit_sha256": sha256_file(args.reuse_audit),
            "pinned_hashes": reuse.get("pinned_hashes"),
            "accepted_split_contract": reuse.get("accepted_split_contract"),
            "extraction_authority": extraction_binding,
        }
        view_slug = args.view.lower()
        qualified_metadata = args.outdir / f"the134_{view_slug}_reused_model_metadata.json"
        qualified_metadata.write_text(
            json.dumps(reused_metadata, indent=2, sort_keys=True) + "\n"
        )
        working_matrix = args.outdir / f"the134_{view_slug}_reuse_weight_working.npz"
        preweight_snapshot = args.outdir / f"the134_{view_slug}_reuse_preweight_snapshot.npz"
        enriched_weighted_matrix = (
            args.outdir / f"the134_{view_slug}_reuse_weighted_enriched.npz"
        )
        weight_registry = args.outdir / f"the134_{view_slug}_reuse_weight_registry.json"
        weight_closure = args.outdir / f"the134_{view_slug}_reuse_weight_closure"
        weight_command = [
            args.python,
            str(args.trainer),
            "--task",
            "tight",
            "--campaign",
            str(defaults["campaign"]),
            "--campaign-spec-ids",
            str(defaults["model_id"]),
            "--cache-file",
            str(working_matrix),
            "--outdir",
            str(args.outdir),
            "--registry-output",
            str(weight_registry),
            "--label-contract",
            "extracted-is-signal",
            "--weight-mode",
            "ppg12-exact",
            "--no-event-weight",
            "--ppg12-exact-expected-samples",
            ",".join(expected_sources(args.system)),
            "--ppg12-exact-closure-dir",
            str(weight_closure),
            "--split-mode",
            str(expected_model_split(args.system, args.view)["mode"]),
            "--test-size",
            str(expected_model_split(args.system, args.view)["test_fraction_requested"]),
            "--random-seed",
            str(expected_model_split(args.system, args.view)["random_seed"]),
            "--cache-only",
        ]
        weight_materialization = {
            "status": "PLANNED_NOT_MATERIALIZED",
            "command": weight_command,
            "working_matrix": str(working_matrix),
            "preweight_snapshot": str(preweight_snapshot),
            "enriched_weighted_matrix": str(enriched_weighted_matrix),
        }
        if args.execute:
            if enriched_weighted_matrix.exists():
                raise SystemExit(
                    "reuse weighted completion artifact exists; refusing duplicate materialization: "
                    f"{enriched_weighted_matrix}"
                )
            if args.system == "pp" and defaults["row_cap_per_class"]:
                cap_report = deterministic_pp_cap(
                    args.matrix,
                    working_matrix,
                    int(defaults["row_cap_per_class"]),
                    int(defaults["seed"]),
                )
            else:
                shutil.copyfile(args.matrix, working_matrix)
                cap_report = {
                    "mode": "full-matrix-no-row-cap",
                    "source_sha256": sha256_file(args.matrix),
                    "output_sha256": sha256_file(working_matrix),
                }
            shutil.copyfile(working_matrix, preweight_snapshot)
            subprocess.run(weight_command, check=True)
            enrichment = enrich_weighted_cache(
                working_matrix, preweight_snapshot, enriched_weighted_matrix
            )
            weight_validation = validate_reuse_weighted_matrix(
                enriched_weighted_matrix
            )
            weight_materialization = {
                **weight_materialization,
                "status": "PASS",
                "row_cap": cap_report,
                "enrichment": enrichment,
                "enriched_weighted_matrix_sha256": sha256_file(
                    enriched_weighted_matrix
                ),
                "validation": weight_validation,
            }
        reuse_plan = {
            "schema": "THE134_FACTORIAL_VIEW_MODEL_REUSE_PLAN_V1",
            "status": "REUSE_READY",
            "system": args.system,
            "shower_definition": args.view,
            "shower_semantic_sha256": shower_semantic_sha256(args.view),
            "model_origin": reuse.get("model_origin", "REUSED_EXISTING"),
            "extraction_authority": extraction_binding,
            "pinned_hashes": reuse.get("pinned_hashes"),
            "reuse_audit": str(args.reuse_audit),
            "reuse_audit_sha256": sha256_file(args.reuse_audit),
            "artifacts": {
                **artifacts,
                "source_metadata": str(original_metadata_path),
                "source_metadata_sha256": sha256_file(original_metadata_path),
                "metadata": str(qualified_metadata),
                "metadata_sha256": sha256_file(qualified_metadata),
                "weighted_matrix": (
                    str(enriched_weighted_matrix) if args.execute else None
                ),
                "weighted_matrix_sha256": (
                    sha256_file(enriched_weighted_matrix) if args.execute else None
                ),
            },
            "weight_materialization": weight_materialization,
        }
        plan_path = args.outdir / f"the134_{view_slug}_model_reuse_plan.json"
        plan_path.write_text(json.dumps(reuse_plan, indent=2, sort_keys=True) + "\n")
        print(plan_path)
        if args.execute:
            receipt = {
                **reuse_plan,
                "schema": "THE134_FACTORIAL_VIEW_MODEL_REUSE_COMPLETE_V1",
                "status": "PASS",
            }
            receipt_path = args.outdir / f"the134_{view_slug}_model_reuse_complete.json"
            receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
            print(receipt_path)
        return 0
    if expected_model_origin(args.system, args.view).startswith("REUSED_"):
        raise SystemExit(
            f"duplicate-run guard: {args.system} {args.view} must reuse its accepted "
            "predecessor after exact semantic/hash audit; provide --reuse-audit. "
            "Retraining is not an implicit fallback."
        )
    if expected_model_origin(args.system, args.view) != "TRAINED_THE134":
        raise SystemExit("internal origin routing error: training branch is not authorized")

    args.outdir.mkdir(parents=True, exist_ok=True)
    view_slug = args.view.lower()
    completion = args.outdir / f"the134_{view_slug}_training_complete.json"
    if args.execute and completion.exists():
        raise SystemExit(f"completion marker exists; refusing duplicate training: {completion}")
    working_matrix = args.outdir / f"the134_{view_slug}_training_matrix_working.npz"
    preweight_snapshot = args.outdir / f"the134_{view_slug}_training_matrix_preweight_snapshot.npz"
    enriched_weighted_matrix = args.outdir / f"the134_{view_slug}_training_matrix_weighted_enriched.npz"
    if not args.execute:
        cap_report = {
            "mode": (
                "clone-established-file-shuffle-within-file-row-cap"
                if defaults["row_cap_per_class"]
                else "full-matrix-no-row-cap"
            ),
            "seed": defaults["seed"],
            "cap_per_class": defaults["row_cap_per_class"],
            "source_sha256": sha256_file(args.matrix),
            "status": "PLANNED_NOT_MATERIALIZED",
        }
    elif args.system == "pp" and defaults["row_cap_per_class"]:
        cap_report = deterministic_pp_cap(
            args.matrix,
            working_matrix,
            int(defaults["row_cap_per_class"]),
            int(defaults["seed"]),
        )
    else:
        shutil.copyfile(args.matrix, working_matrix)
        cap_report = {
            "mode": "full-matrix-no-row-cap",
            "source_sha256": sha256_file(args.matrix),
            "output_sha256": sha256_file(working_matrix),
        }
    if args.execute:
        shutil.copyfile(working_matrix, preweight_snapshot)

    registry = args.outdir / "model_registry.json"
    closure = args.outdir / "ppg12_exact_weight_closure"
    command = [
        args.python,
        str(args.trainer),
        "--task",
        "tight",
        "--campaign",
        str(defaults["campaign"]),
        "--campaign-spec-ids",
        str(defaults["model_id"]),
        "--cache-file",
        str(working_matrix),
        "--outdir",
        str(args.outdir),
        "--registry-output",
        str(registry),
        "--label-contract",
        "extracted-is-signal",
        "--weight-mode",
        "ppg12-exact",
        "--no-event-weight",
        "--ppg12-exact-expected-samples",
        ",".join(expected_sources(args.system)),
        "--ppg12-exact-closure-dir",
        str(closure),
        "--split-mode",
        "event50",
        "--require-global-event-key",
        "--test-size",
        str(defaults["test_size"]),
        "--random-seed",
        str(defaults["seed"]),
        "--pt-bins",
        ",".join(f"{edge:g}" for edge in READER_ET_EDGES),
        "--parallel-workers",
        "1",
        "--n-jobs",
        "2",
        "--n-estimators",
        "450",
        "--max-depth",
        "4",
        "--learning-rate",
        "0.035",
        "--subsample",
        "0.85",
        "--colsample-bytree",
        "0.85",
        "--reg-alpha",
        "5.0",
        "--reg-lambda",
        "0.3",
        "--grow-policy",
        "lossguide",
        "--max-bin",
        "256",
        "--tree-method",
        "hist",
    ]
    plan = {
        "schema": "THE134_FACTORIAL_VIEW_MODEL_TRAINING_PLAN_V1",
        "status": "EXECUTING" if args.execute else "DRY_RUN_READY",
        "system": args.system,
        "shower_definition": args.view,
        "shower_semantic_sha256": shower_semantic_sha256(args.view),
        "model_origin": "TRAINED_THE134",
        "extraction_authority": extraction_binding,
        "feature_order": list(FEATURES_BY_SYSTEM[args.system]),
        "model_domain_gev": list(MODEL_DOMAIN_GEV),
        "model_id": defaults["model_id"],
        "event_split": {
            "mode": "event50",
            "seed": defaults["seed"],
            "test_fraction": defaults["test_size"],
        },
        "row_cap": cap_report,
        "matrix": str(args.matrix),
        "matrix_sha256": sha256_file(args.matrix),
        "working_matrix": str(working_matrix),
        "preweight_snapshot": str(preweight_snapshot),
        "enriched_weighted_matrix": str(enriched_weighted_matrix),
        "extraction_audit": str(args.extraction_audit),
        "extraction_audit_sha256": sha256_file(args.extraction_audit),
        "trainer": str(args.trainer),
        "trainer_sha256": sha256_file(args.trainer),
        "command": command,
    }
    plan_path = args.outdir / f"the134_{view_slug}_training_plan.json"
    plan_path.write_text(json.dumps(plan, indent=2, sort_keys=True) + "\n")
    print(plan_path)
    if not args.execute:
        return 0

    subprocess.run(command, check=True)
    stem = str(defaults["output_stem"])
    source_tmva = args.outdir / f"{stem}.root"
    source_xgb = args.outdir / f"{stem}.xgb.json"
    source_metadata = args.outdir / f"{stem}.metadata.json"
    for path in (source_tmva, source_xgb, source_metadata, registry):
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"training did not produce required artifact: {path}")
    tmva = args.outdir / f"the134_{args.system}_{view_slug}_model.tmva.root"
    xgb = args.outdir / f"the134_{args.system}_{view_slug}_model.xgb.json"
    shutil.copyfile(source_tmva, tmva)
    shutil.copyfile(source_xgb, xgb)
    metadata = json.loads(source_metadata.read_text())
    metadata["model_origin"] = expected_model_origin(args.system, args.view)
    metadata["accepted_split_contract"] = expected_model_split(args.system, args.view)
    weighting = metadata.get("weighting", {})
    training_contract_gates = {
        "weight_mode": metadata.get("weight_mode") == "ppg12-exact",
        "event_weight_unused": weighting.get("event_weight_used") is False,
        "vertex_weight_unused": weighting.get("vertex_reweight") is False,
        "centrality_weight_unused": weighting.get("centrality_event_weight") is False,
        "cross_section_weight_unused": weighting.get("cross_section_weight_used_for_training")
        is False,
        "weights_computed_before_binning": weighting.get("weights_computed_before_binning")
        is True,
        "event_group_split": metadata.get("split", {}).get("mode") == "event50",
        "test_fraction": metadata.get("split", {}).get("test_fraction_requested")
        == defaults["test_size"],
        "old_low_calo_veto_disabled": not bool(
            metadata.get("event_quality_filter", {}).get("enabled", False)
        ),
    }
    if not all(training_contract_gates.values()):
        raise SystemExit(
            "trained metadata violates frozen THE-134 protocol: "
            + json.dumps(training_contract_gates, sort_keys=True)
        )
    metadata["the134_view_contract"] = {
        "schema": "THE134_FACTORIAL_VIEW_MODEL_VIEW_CONTRACT_V1",
        "system": args.system,
        "shower_definition": args.view,
        "shower_semantic_sha256": shower_semantic_sha256(args.view),
        "feature_order": list(FEATURES_BY_SYSTEM[args.system]),
        "model_domain_gev": list(MODEL_DOMAIN_GEV),
        "label_authority": "frozen source-role training_label materialized without relabeling",
        "extraction_audit": str(args.extraction_audit),
        "extraction_audit_sha256": sha256_file(args.extraction_audit),
        "source_metadata": str(source_metadata),
        "source_metadata_sha256": sha256_file(source_metadata),
        "training_protocol_gates": training_contract_gates,
        "extraction_authority": extraction_binding,
    }
    qualified_metadata = args.outdir / f"the134_{view_slug}_model_metadata.json"
    qualified_metadata.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
    enrichment = enrich_weighted_cache(
        working_matrix, preweight_snapshot, enriched_weighted_matrix
    )
    receipt = build_training_completion_receipt(
        system=args.system,
        view_name=args.view,
        model_id=defaults["model_id"],
        extraction_binding=extraction_binding,
        artifacts={
            "xgboost": str(xgb),
            "xgboost_sha256": sha256_file(xgb),
            "tmva": str(tmva),
            "tmva_sha256": sha256_file(tmva),
            "metadata": str(qualified_metadata),
            "metadata_sha256": sha256_file(qualified_metadata),
            "registry": str(registry),
            "registry_sha256": sha256_file(registry),
            "weighted_matrix": str(enriched_weighted_matrix),
            "weighted_matrix_sha256": sha256_file(enriched_weighted_matrix),
            "training_cache": str(working_matrix),
            "training_cache_sha256": sha256_file(working_matrix),
        },
        enrichment=enrichment,
    )
    completion.write_text(json.dumps(receipt, indent=2, sort_keys=True) + "\n")
    print(completion)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
