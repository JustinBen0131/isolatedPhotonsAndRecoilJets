#!/usr/bin/env python3
"""Validate THE-119 p+p archived-DI direct/writer canary pairs.

The direct ROOT remains the authoritative science output.  The writer ROOT
must be bit-for-bit neutral for every common histogram while also carrying a
complete, internally coherent ReplayFoundationV1 payload.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import sys
from pathlib import Path
from typing import Any

import ROOT


ROOT.gROOT.SetBatch(True)

EXPECTED_TREES = {
    "RJSourceOccurrenceV1",
    "RJEventV1",
    "RJPhotonCandidateV1",
    "RJModelEvaluationV1",
    "RJShowerCellV1",
    "RJIsolationConstituentV1",
    "RJIsolationWitnessV1",
    "RJJetV1",
    "RJJetConstituentV1",
    "RJPhotonJetPairV1",
    "RJTruthPhotonV1",
    "RJTruthJetV1",
    "RJRecoTruthLinkV1",
    "RJWeightComponentV1",
    "RJEventDisplaySnapshotV1",
}

HEX64 = re.compile(r"^[0-9a-f]{64}$")


def fail(errors: list[str], message: str) -> None:
    if len(errors) < 100:
        errors.append(message)


def open_root(path: str, errors: list[str]) -> Any:
    if not os.path.isfile(path):
        fail(errors, f"missing ROOT: {path}")
        return None
    if os.path.getsize(path) < 50_000:
        fail(errors, f"ROOT below 50 kB gate: {path}")
    root_file = ROOT.TFile.Open(path, "READ")
    if not root_file or root_file.IsZombie():
        fail(errors, f"unreadable or zombie ROOT: {path}")
        return None
    if root_file.TestBit(ROOT.TFile.kRecovered):
        fail(errors, f"recovered ROOT is forbidden: {path}")
    return root_file


def collect_histograms(directory: Any, prefix: str = "") -> dict[str, Any]:
    output: dict[str, Any] = {}
    keys = directory.GetListOfKeys()
    for key in keys:
        name = key.GetName()
        path = f"{prefix}/{name}" if prefix else name
        if path == "ReplayFoundationV1" or path.startswith("ReplayFoundationV1/"):
            continue
        obj = key.ReadObj()
        if obj.InheritsFrom("TDirectory"):
            output.update(collect_histograms(obj, path))
        elif obj.InheritsFrom("TH1"):
            output[path] = obj
    return output


def axis_signature(axis: Any) -> tuple[Any, ...]:
    bins = axis.GetXbins()
    variable = tuple(float(bins.At(i)) for i in range(bins.GetSize()))
    return (
        int(axis.GetNbins()),
        float(axis.GetXmin()),
        float(axis.GetXmax()),
        variable,
    )


def compare_histograms(direct: Any, writer: Any, errors: list[str]) -> dict[str, Any]:
    direct_hists = collect_histograms(direct)
    writer_hists = collect_histograms(writer)
    if set(direct_hists) != set(writer_hists):
        missing = sorted(set(direct_hists) - set(writer_hists))
        extra = sorted(set(writer_hists) - set(direct_hists))
        fail(errors, f"histogram key mismatch missing={missing[:8]} extra={extra[:8]}")

    compared = 0
    bins_compared = 0
    for path in sorted(set(direct_hists) & set(writer_hists)):
        left = direct_hists[path]
        right = writer_hists[path]
        if left.ClassName() != right.ClassName():
            fail(errors, f"histogram class mismatch {path}: {left.ClassName()} != {right.ClassName()}")
            continue
        shape_left = (
            int(left.GetDimension()),
            axis_signature(left.GetXaxis()),
            axis_signature(left.GetYaxis()),
            axis_signature(left.GetZaxis()),
            int(left.GetNcells()),
        )
        shape_right = (
            int(right.GetDimension()),
            axis_signature(right.GetXaxis()),
            axis_signature(right.GetYaxis()),
            axis_signature(right.GetZaxis()),
            int(right.GetNcells()),
        )
        if shape_left != shape_right:
            fail(errors, f"histogram binning mismatch: {path}")
            continue
        for index in range(left.GetNcells()):
            bins_compared += 1
            if float(left.GetBinContent(index)) != float(right.GetBinContent(index)):
                fail(errors, f"bin-content mismatch: {path} bin={index}")
                break
            if float(left.GetBinError(index)) != float(right.GetBinError(index)):
                fail(errors, f"bin-error mismatch: {path} bin={index}")
                break
        sumw2_left = left.GetSumw2()
        sumw2_right = right.GetSumw2()
        if sumw2_left.GetSize() != sumw2_right.GetSize():
            fail(errors, f"Sumw2 size mismatch: {path}")
        else:
            for index in range(sumw2_left.GetSize()):
                if float(sumw2_left.At(index)) != float(sumw2_right.At(index)):
                    fail(errors, f"Sumw2 mismatch: {path} bin={index}")
                    break
        if left.InheritsFrom("TProfile"):
            for index in range(left.GetNcells()):
                if float(left.GetBinEntries(index)) != float(right.GetBinEntries(index)):
                    fail(errors, f"profile-entry mismatch: {path} bin={index}")
                    break
        compared += 1
    return {"histograms_compared": compared, "histogram_bins_compared": bins_compared}


def named_title(directory: Any, name: str) -> str | None:
    obj = directory.Get(name)
    if not obj:
        return None
    return str(obj.GetTitle())


def identity(row: Any, stem: str) -> tuple[int, int]:
    return (int(getattr(row, f"{stem}_hi")), int(getattr(row, f"{stem}_lo")))


def finite_or_fail(value: float, errors: list[str], label: str) -> float:
    number = float(value)
    if not math.isfinite(number):
        fail(errors, f"nonfinite {label}")
    return number


def validate_writer(
    root_file: Any,
    errors: list[str],
    *,
    expected: dict[str, str],
    expected_period: str,
    candidate_model_hash: str,
) -> dict[str, Any]:
    replay = root_file.Get("ReplayFoundationV1")
    if not replay or not replay.InheritsFrom("TDirectory"):
        fail(errors, "missing ReplayFoundationV1 directory")
        return {}

    tree_names = {
        key.GetName()
        for key in replay.GetListOfKeys()
        if key.GetClassName() == "TTree"
    }
    if tree_names != EXPECTED_TREES:
        fail(errors, f"tree inventory mismatch: {sorted(tree_names)}")

    metadata = {
        key: named_title(replay, key)
        for key in (
            "rj_replay_schema",
            "rj_replay_schema_version",
            "rj_replay_complete",
            "schema_sha256",
            "semantic_sha256",
            "config_sha256",
            "model_sha256",
            "source_sha256",
            "code_sha256",
        )
    }
    if metadata["rj_replay_schema"] != "RJ_REPLAY_FOUNDATION_V1":
        fail(errors, "wrong replay schema identity")
    if metadata["rj_replay_schema_version"] != "1":
        fail(errors, "wrong replay schema version")
    if metadata["rj_replay_complete"] != "1":
        fail(errors, "writer completion marker is not 1")
    for key, wanted in expected.items():
        if metadata.get(key) != wanted:
            fail(errors, f"metadata mismatch {key}: {metadata.get(key)} != {wanted}")
    if not metadata["source_sha256"] or not HEX64.fullmatch(metadata["source_sha256"]):
        fail(errors, "source_sha256 is not lowercase SHA-256")

    trees = {name: replay.Get(name) for name in EXPECTED_TREES}
    entries = {name: int(tree.GetEntries()) for name, tree in trees.items()}

    source = trees["RJSourceOccurrenceV1"]
    if entries["RJSourceOccurrenceV1"] != 1:
        fail(errors, f"source occurrence count is {entries['RJSourceOccurrenceV1']}, expected 1")
    else:
        source.GetEntry(0)
        role = str(source.si_di_role).upper()
        period = str(source.period)
        if role != "DI":
            fail(errors, f"source SI/DI role is {role}, expected DI")
        if period != expected_period:
            fail(errors, f"source period is {period}, expected {expected_period}")
        for field in ("lane", "dataset", "sample", "ownership_state", "source_manifest_sha256"):
            if not str(getattr(source, field)):
                fail(errors, f"empty source field: {field}")

    events = trees["RJEventV1"]
    event_ids: set[tuple[int, int]] = set()
    candidate_count_sum = 0
    for event in events:
        event_id = identity(event, "event_id")
        if event_id == (0, 0) or event_id in event_ids:
            fail(errors, "null or duplicate event identity")
            break
        event_ids.add(event_id)
        candidate_count_sum += int(event.candidate_count)

    candidates = trees["RJPhotonCandidateV1"]
    candidate_et: dict[tuple[int, int], float] = {}
    for candidate in candidates:
        candidate_id = identity(candidate, "candidate_id")
        if candidate_id == (0, 0) or candidate_id in candidate_et:
            fail(errors, "null or duplicate candidate identity")
            break
        if identity(candidate, "event_id") not in event_ids:
            fail(errors, "candidate references unknown event")
            break
        et = finite_or_fail(candidate.cluster_et, errors, "candidate ET")
        if not (5.0 <= et < 40.0):
            fail(errors, f"candidate ET outside frozen capture domain: {et}")
        if (et < 15.0) != (int(candidate.below15_retention_state) == 1):
            fail(errors, "below-15 retention state disagrees with candidate ET")
        candidate_et[candidate_id] = et
    if candidate_count_sum != entries["RJPhotonCandidateV1"]:
        fail(errors, "event candidate-count sum does not close to candidate tree")

    models = trees["RJModelEvaluationV1"]
    below15_rows = 0
    validated_candidate_rows = 0
    model_keys: set[tuple[tuple[int, int], tuple[int, int]]] = set()
    for model in models:
        candidate_id = identity(model, "candidate_id")
        if candidate_id not in candidate_et:
            fail(errors, "model row references unknown candidate")
            break
        key = (candidate_id, identity(model, "model_id"))
        if key in model_keys:
            fail(errors, "duplicate candidate/model evaluation identity")
            break
        model_keys.add(key)
        et = candidate_et[candidate_id]
        state = int(model.applicability_state)
        witnesses = list(model.ordered_input_witnesses)
        if len(witnesses) != 11 or any(not math.isfinite(float(value)) for value in witnesses):
            fail(errors, "model row does not carry 11 finite ordered witnesses")
        if int(model.finite_score) != 1 or not math.isfinite(float(model.raw_score)):
            fail(errors, "model score is not finite")
        wp_values = [
            float(model.wp70),
            float(model.wp80),
            float(model.wp90),
            float(model.delta_wp70),
            float(model.delta_wp80),
            float(model.delta_wp90),
        ]
        if et < 15.0:
            below15_rows += 1
            if state != 1 or any(math.isfinite(value) for value in wp_values):
                fail(errors, "below-15 model row is not diagnostic-only with null WP state")
        elif et < 35.0:
            if state != 0:
                fail(errors, "15-35 model row is not VALIDATED_DOMAIN")
            if str(model.model_sha256) == candidate_model_hash:
                validated_candidate_rows += 1
                if any(not math.isfinite(value) for value in wp_values):
                    fail(errors, "candidate model lacks finite WP state in validated domain")
        else:
            if state != 1 or any(math.isfinite(value) for value in wp_values):
                fail(errors, "above-35 model row is not diagnostic-only with null WP state")

    weights = trees["RJWeightComponentV1"]
    if entries["RJWeightComponentV1"] != entries["RJEventV1"]:
        fail(errors, "event/weight population mismatch")
    weight_targets: set[tuple[int, int]] = set()
    for weight in weights:
        target = identity(weight, "target_id")
        if target not in event_ids or target in weight_targets:
            fail(errors, "weight target is unknown or duplicated")
            break
        weight_targets.add(target)
        if str(weight.component_type) != "event" or int(weight.application_count) != 1:
            fail(errors, "weight is not one event application")
        values = [
            float(weight.slice_weight),
            float(weight.cross_section_weight),
            float(weight.vertex_weight),
            float(weight.si_di_weight),
            float(weight.period_weight),
            float(weight.exposure_weight),
            float(weight.final_weight),
        ]
        if any(not math.isfinite(value) or value <= 0.0 for value in values):
            fail(errors, "weight contains a nonfinite or nonpositive component")
        if values[0] != values[1] or values[4] != values[5]:
            fail(errors, "weight witness duplication does not close")
        product = values[0] * values[2] * values[3] * values[4]
        if not math.isclose(product, values[6], rel_tol=2e-12, abs_tol=1e-15):
            fail(errors, f"final event weight does not close: {product} != {values[6]}")

    return {
        "metadata": metadata,
        "tree_entries": entries,
        "event_identity_count": len(event_ids),
        "candidate_identity_count": len(candidate_et),
        "model_identity_count": len(model_keys),
        "below15_model_rows": below15_rows,
        "validated_candidate_model_rows": validated_candidate_rows,
        "weight_once_rows": len(weight_targets),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--pair",
        action="append",
        nargs=4,
        metavar=("LABEL", "DIRECT_ROOT", "WRITER_ROOT", "EXPECTED_PERIOD"),
        required=True,
    )
    parser.add_argument("--output", type=Path)
    parser.add_argument("--schema-sha256", required=True)
    parser.add_argument("--semantic-sha256", required=True)
    parser.add_argument("--config-sha256", required=True)
    parser.add_argument("--code-sha256", required=True)
    parser.add_argument("--model-sha256", required=True)
    args = parser.parse_args()

    report: dict[str, Any] = {
        "schema": "THE119_PP_DI_CANARY_VALIDATION_V1",
        "status": "PASS",
        "pairs": {},
    }
    expected = {
        "schema_sha256": args.schema_sha256,
        "semantic_sha256": args.semantic_sha256,
        "config_sha256": args.config_sha256,
        "code_sha256": args.code_sha256,
        "model_sha256": args.model_sha256,
    }

    all_errors: list[str] = []
    for label, direct_path, writer_path, expected_period in args.pair:
        errors: list[str] = []
        direct = open_root(direct_path, errors)
        writer = open_root(writer_path, errors)
        pair_report: dict[str, Any] = {
            "direct_root": direct_path,
            "writer_root": writer_path,
            "direct_bytes": os.path.getsize(direct_path) if os.path.isfile(direct_path) else None,
            "writer_bytes": os.path.getsize(writer_path) if os.path.isfile(writer_path) else None,
            "expected_period": expected_period,
        }
        if direct and writer:
            pair_report.update(compare_histograms(direct, writer, errors))
            pair_report["writer_validation"] = validate_writer(
                writer,
                errors,
                expected=expected,
                expected_period=expected_period,
                candidate_model_hash=args.model_sha256,
            )
        if direct:
            direct.Close()
        if writer:
            writer.Close()
        pair_report["status"] = "PASS" if not errors else "FAIL"
        pair_report["errors"] = errors
        report["pairs"][label] = pair_report
        all_errors.extend(f"{label}: {error}" for error in errors)

    if all_errors:
        report["status"] = "FAIL"
    report["errors"] = all_errors
    payload = json.dumps(report, indent=2, sort_keys=True)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload + "\n")
    print(payload)
    return 0 if report["status"] == "PASS" else 2


if __name__ == "__main__":
    sys.exit(main())
