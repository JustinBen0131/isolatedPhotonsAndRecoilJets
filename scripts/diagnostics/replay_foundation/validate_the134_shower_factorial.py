#!/usr/bin/env python3
"""Certify the real THE-134 multi-view shower payload and active models.

The calculation below is intentionally independent of RJShowerFactorialV1.h.
It rebuilds every view from the absolute tower rows, checks the active
PhotonClusterBuilder feature vector against its declared view, and recomputes
the TMVA score.  It never chooses or retunes a working point or sideband.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np
import uproot

from validate_the119_jet_response import compare_histograms, _pyroot_health


TREE_NAMES = {
    "RJSourceOccurrenceV1", "RJEventV1", "RJPhotonCandidateV1",
    "RJModelEvaluationV1", "RJShowerCellV1", "RJShowerFeatureViewV1",
    "RJIsolationConstituentV1", "RJIsolationWitnessV1", "RJJetV1",
    "RJJetConstituentV1", "RJPhotonJetPairV1", "RJTruthPhotonV1",
    "RJTruthJetV1", "RJRecoTruthLinkV1", "RJWeightComponentV1",
    "RJEventDisplaySnapshotV1",
}

# energy source, rectangular membership, moment membership, floor [GeV]
DEFINITIONS: dict[str, tuple[int, int, int, float]] = {
    "H70": (0, 0, 1, 0.070), "H0": (0, 0, 1, 0.0),
    "G70": (0, 0, 0, 0.070), "G0": (0, 0, 0, 0.0),
    "O70": (0, 1, 1, 0.070), "O0": (0, 1, 1, 0.0),
    "R70": (1, 1, 1, 0.070),
}


def text(value: Any) -> str:
    return value.decode() if isinstance(value, bytes) else str(value)


def identity(hi: Any, lo: Any) -> tuple[int, int]:
    return int(hi), int(lo)


def semantic_text(name: str) -> str:
    energy, sums, moments, floor = DEFINITIONS[name]
    return (
        f"RJ_SHOWER_DEFINITION_FACTORIAL_V1|{name}|energy={energy}"
        f"|sums={sums}|moments={moments}"
        f"|floor_gev={'0.070000' if floor > 0 else '0.000000'}"
        "|grid=7x7|tower_quality=TowerInfo_get_isGood_only"
        "|center_excluded_from_cogx_numerator=1"
    )


def semantic_sha(name: str) -> str:
    return hashlib.sha256(semantic_text(name).encode()).hexdigest()


def close(left: float, right: float, tolerance: float) -> bool:
    if math.isnan(left) and math.isnan(right):
        return True
    return math.isclose(left, right, rel_tol=0.0, abs_tol=tolerance)


def view_from_cells(name: str, row: dict[str, Any], cells: list[dict[str, Any]]) -> dict[str, Any]:
    energy_source, sum_membership, moment_membership, floor = DEFINITIONS[name]
    center_eta = int(row["center_eta_index"])
    center_phi = int(row["center_phi_index"])
    cog_eta = 3.0 + (float(row["raw_center_eta"]) % 1.0 - 0.5)
    cog_phi = 3.0 + (float(row["raw_center_phi"]) % 1.0 - 0.5)
    sign_phi = 1 if cog_phi > 3.0 else -1
    out = {
        "e11": 0.0, "e33": 0.0, "e32": 0.0, "e35": 0.0,
        "moment_eta_numerator": 0.0, "moment_phi_numerator": 0.0,
        "moment_denominator": 0.0, "moment33_eta_numerator": 0.0,
        "moment33_phi_numerator": 0.0, "moment33_denominator": 0.0,
        "good_cell_count": 0, "owned_cell_count": 0,
        "active_sum_cell_count": 0, "active_moment_cell_count": 0,
        "exact_zero_count": 0, "negative_count": 0, "nonfinite_count": 0,
    }
    for cell in cells:
        i = int(cell["tower_eta_index"]) - center_eta + 3
        delta_phi = int(cell["tower_phi_index"]) - center_phi
        while delta_phi < -128:
            delta_phi += 256
        while delta_phi > 127:
            delta_phi -= 256
        j = delta_phi + 3
        if not (0 <= i <= 6 and 0 <= j <= 6):
            continue
        out["good_cell_count"] += int(cell["is_good"] != 0)
        out["owned_cell_count"] += int(cell["rawcluster_owned"] != 0)
        out["exact_zero_count"] += int(cell["is_good"] != 0 and cell["is_zero"] != 0)
        out["negative_count"] += int(cell["is_good"] != 0 and cell["is_negative"] != 0)
        out["nonfinite_count"] += int(cell["is_nonfinite"] != 0)
        if energy_source == 0:
            energy = float(cell["calibrated_energy"])
            selected = int(cell["is_good"]) != 0 and math.isfinite(energy) and energy > floor
        else:
            energy = float(cell["rawcluster_map_value"])
            selected = (
                int(cell["rawcluster_owned"]) != 0
                and int(cell["rawcluster_value_present"]) != 0
                and math.isfinite(energy) and energy > floor
            )
        if not selected:
            continue
        owned = int(cell["rawcluster_owned"]) != 0
        sum_member = sum_membership == 0 or owned
        moment_member = moment_membership == 0 or owned
        di, dj = abs(i - 3), abs(j - 3)
        if sum_member:
            out["active_sum_cell_count"] += 1
            if i == 3 and j == 3:
                out["e11"] += energy
            if di <= 1 and dj <= 1:
                out["e33"] += energy
            if di <= 1 and (j == 3 or j == 3 + sign_phi):
                out["e32"] += energy
            if di <= 1 and dj <= 2:
                out["e35"] += energy
        if moment_member:
            out["active_moment_cell_count"] += 1
            deta, dphi = i - cog_eta, j - cog_phi
            out["moment_denominator"] += energy
            if i != 3 or j != 3:
                out["moment_eta_numerator"] += energy * deta * deta
                out["moment_phi_numerator"] += energy * dphi * dphi
            if di <= 1 and dj <= 1:
                out["moment33_denominator"] += energy
                if i != 3 or j != 3:
                    out["moment33_eta_numerator"] += energy * deta * deta
                    out["moment33_phi_numerator"] += energy * dphi * dphi
    out["e11_over_e33"] = out["e11"] / out["e33"] if out["e33"] > 0 else math.nan
    out["e32_over_e35"] = out["e32"] / out["e35"] if out["e35"] > 0 else math.nan
    out["weta_cogx"] = out["moment_eta_numerator"] / out["moment_denominator"] if out["moment_denominator"] > 0 else math.nan
    out["wphi_cogx"] = out["moment_phi_numerator"] / out["moment_denominator"] if out["moment_denominator"] > 0 else math.nan
    out["weta33_cogx"] = out["moment33_eta_numerator"] / out["moment33_denominator"] if out["moment33_denominator"] > 0 else math.nan
    out["wphi33_cogx"] = out["moment33_phi_numerator"] / out["moment33_denominator"] if out["moment33_denominator"] > 0 else math.nan
    return out


def rows_by_candidate(tree: Any, branches: list[str]) -> dict[tuple[int, int], list[dict[str, Any]]]:
    arrays = tree.arrays(branches, library="np")
    result: dict[tuple[int, int], list[dict[str, Any]]] = defaultdict(list)
    for index in range(len(arrays[branches[0]])):
        key = identity(arrays["candidate_id_hi"][index], arrays["candidate_id_lo"][index])
        result[key].append({name: arrays[name][index] for name in branches})
    return result


def validate_writer(
    label: str, system: str, active_definition: str, path: Path,
    expected: dict[str, str], models: dict[str, tuple[Path, str]],
    minimum_bytes: int, feature_tolerance: float, score_tolerance: float,
) -> dict[str, Any]:
    failures: list[str] = []
    report: dict[str, Any] = {"label": label, "path": str(path), "failures": failures}
    report["size_bytes"] = path.stat().st_size if path.exists() else 0
    if not path.exists() or report["size_bytes"] < minimum_bytes:
        failures.append("missing_or_tiny_root")
        report["status"] = "FAIL"
        return report
    healthy, detail = _pyroot_health(path)
    report["root_health"] = detail
    if not healthy:
        failures.append(detail)
        report["status"] = "FAIL"
        return report

    with uproot.open(path) as root:
        if "ReplayFoundationV1" not in root:
            failures.append("missing_replay_directory")
            report["status"] = "FAIL"
            return report
        replay = root["ReplayFoundationV1"]
        version = text(replay["rj_replay_schema_version"].member("fTitle"))
        report["schema_version"] = version
        if version != "2":
            failures.append("schema_version_not_2")
        names = {key.split(";")[0] for key, cls in replay.classnames(recursive=False).items() if cls == "TTree"}
        report["tree_count"] = len(names)
        if names != TREE_NAMES:
            failures.append("tree_inventory_mismatch")
            report["missing_trees"] = sorted(TREE_NAMES - names)
            report["extra_trees"] = sorted(names - TREE_NAMES)
        metadata = {name: text(replay[name].member("fTitle")) for name in expected}
        report["metadata"] = metadata
        for name, value in expected.items():
            if metadata.get(name) != value:
                failures.append(f"{name}_mismatch")

        event_arrays = replay["RJEventV1"].arrays(
            ["event_id_hi", "event_id_lo", "vertex_z", "centrality"], library="np"
        )
        events = {
            identity(event_arrays["event_id_hi"][i], event_arrays["event_id_lo"][i]): (
                float(event_arrays["vertex_z"][i]), float(event_arrays["centrality"][i])
            ) for i in range(len(event_arrays["vertex_z"]))
        }
        candidate_arrays = replay["RJPhotonCandidateV1"].arrays(
            ["candidate_id_hi", "candidate_id_lo", "event_id_hi", "event_id_lo",
             "cluster_et", "eta", "ordered_features", "shower_definition_views"], library="np"
        )
        candidates: dict[tuple[int, int], dict[str, Any]] = {}
        for i in range(len(candidate_arrays["cluster_et"])):
            key = identity(candidate_arrays["candidate_id_hi"][i], candidate_arrays["candidate_id_lo"][i])
            candidates[key] = {
                "event": identity(candidate_arrays["event_id_hi"][i], candidate_arrays["event_id_lo"][i]),
                "et": float(candidate_arrays["cluster_et"][i]),
                "eta": float(candidate_arrays["eta"][i]),
                "features": np.asarray(candidate_arrays["ordered_features"][i], dtype=np.float32),
                "views": {text(value) for value in candidate_arrays["shower_definition_views"][i]},
            }

        cell_branches = [
            "candidate_id_hi", "candidate_id_lo", "tower_key", "tower_eta_index", "tower_phi_index",
            "calibrated_energy", "rawcluster_map_value", "is_good", "is_zero", "is_negative",
            "is_nonfinite", "rawcluster_owned", "rawcluster_value_present", "grid_membership_bitmask",
        ]
        cell_rows = rows_by_candidate(replay["RJShowerCellV1"], cell_branches)
        view_branches = [
            "candidate_id_hi", "candidate_id_lo", "definition_name", "semantic_sha256",
            "ordered_features", "floor_gev", "raw_center_eta", "raw_center_phi",
            "center_eta_index", "center_phi_index", "energy_source", "rectangular_membership",
            "moment_membership", "native_et1", "native_et2", "native_et3", "native_et4",
            "e11", "e33", "e32", "e35", "e11_over_e33", "e32_over_e35",
            "weta_cogx", "wphi_cogx", "weta33_cogx", "wphi33_cogx",
            "moment_eta_numerator", "moment_phi_numerator", "moment_denominator",
            "moment33_eta_numerator", "moment33_phi_numerator", "moment33_denominator",
            "good_cell_count", "owned_cell_count", "active_sum_cell_count",
            "active_moment_cell_count", "exact_zero_count", "negative_count", "nonfinite_count",
            "finite_feature_state",
        ]
        view_rows = rows_by_candidate(replay["RJShowerFeatureViewV1"], view_branches)
        scalar_fields = [
            "e11", "e33", "e32", "e35", "e11_over_e33", "e32_over_e35",
            "weta_cogx", "wphi_cogx", "weta33_cogx", "wphi33_cogx",
            "moment_eta_numerator", "moment_phi_numerator", "moment_denominator",
            "moment33_eta_numerator", "moment33_phi_numerator", "moment33_denominator",
        ]
        count_fields = [
            "good_cell_count", "owned_cell_count", "active_sum_cell_count",
            "active_moment_cell_count", "exact_zero_count", "negative_count", "nonfinite_count",
        ]
        feature_mismatches = 0
        replay_mismatches = 0
        bad_view_sets = 0
        view_lookup: dict[tuple[tuple[int, int], str], np.ndarray] = {}
        for key, candidate in candidates.items():
            candidate_views = view_rows.get(key, [])
            names_for_candidate = {text(row["definition_name"]) for row in candidate_views}
            if candidate["views"] != set(DEFINITIONS) or names_for_candidate != set(DEFINITIONS):
                bad_view_sets += 1
            coordinate_keys = [
                (int(row["tower_eta_index"]), int(row["tower_phi_index"])) for row in cell_rows.get(key, [])
            ]
            if len(coordinate_keys) != len(set(coordinate_keys)):
                failures.append("duplicate_absolute_tower_identity")
            if any(not (0 <= eta < 96 and 0 <= phi < 256) for eta, phi in coordinate_keys):
                failures.append("absolute_tower_coordinate_out_of_range")
            if any(int(row["grid_membership_bitmask"]) <= 0 for row in cell_rows.get(key, [])):
                failures.append("missing_grid_membership")
            vertex_z, centrality = events[candidate["event"]]
            for row in candidate_views:
                name = text(row["definition_name"])
                if name not in DEFINITIONS or text(row["semantic_sha256"]) != semantic_sha(name):
                    replay_mismatches += 1
                    continue
                expected_energy, expected_sum, expected_moment, expected_floor = DEFINITIONS[name]
                if (
                    int(row["energy_source"]) != expected_energy
                    or int(row["rectangular_membership"]) != expected_sum
                    or int(row["moment_membership"]) != expected_moment
                    or not close(float(row["floor_gev"]), expected_floor, 1.0e-15)
                ):
                    replay_mismatches += 1
                rebuilt = view_from_cells(name, row, cell_rows.get(key, []))
                replay_mismatches += sum(
                    not close(float(row[field]), float(rebuilt[field]), 2.0e-12) for field in scalar_fields
                )
                replay_mismatches += sum(int(row[field]) != int(rebuilt[field]) for field in count_fields)
                native = [float(row[f"native_et{i}"]) for i in range(1, 5)]
                expected_features = [
                    candidate["et"], rebuilt["weta_cogx"], rebuilt["wphi_cogx"], vertex_z,
                    candidate["eta"], rebuilt["e11_over_e33"], *native, rebuilt["e32_over_e35"],
                ]
                if system == "auau":
                    expected_features += [centrality, rebuilt["weta33_cogx"], rebuilt["wphi33_cogx"]]
                observed_features = np.asarray(row["ordered_features"], dtype=np.float32)
                expected_features_np = np.asarray(expected_features, dtype=np.float32)
                if observed_features.shape != expected_features_np.shape or not np.allclose(
                    observed_features, expected_features_np, rtol=0.0, atol=feature_tolerance, equal_nan=True
                ):
                    feature_mismatches += 1
                view_lookup[(key, name)] = observed_features
            active = view_lookup.get((key, active_definition))
            if active is None or candidate["features"].shape != active.shape or not np.allclose(
                candidate["features"], active, rtol=0.0, atol=feature_tolerance, equal_nan=True
            ):
                feature_mismatches += 1

        report.update({
            "candidate_count": len(candidates), "shower_cell_count": sum(map(len, cell_rows.values())),
            "shower_view_count": sum(map(len, view_rows.values())), "bad_view_sets": bad_view_sets,
            "replay_mismatches": replay_mismatches, "feature_mismatches": feature_mismatches,
        })
        if not candidates:
            failures.append("zero_candidate_witness")
        if bad_view_sets:
            failures.append("factorial_view_inventory_mismatch")
        if replay_mismatches:
            failures.append("independent_cell_replay_mismatch")
        if feature_mismatches:
            failures.append("active_or_view_feature_mismatch")

        model_branches = [
            "candidate_id_hi", "candidate_id_lo", "model_sha256", "shower_definition_id",
            "shower_semantic_sha256", "ordered_input_witnesses", "raw_score",
            "finite_score", "applicability_state", "wp70", "wp80", "wp90",
        ]
        model_arrays = replay["RJModelEvaluationV1"].arrays(model_branches, library="np")
        score_inputs: dict[str, list[np.ndarray]] = defaultdict(list)
        score_observed: dict[str, list[float]] = defaultdict(list)
        bad_model_semantics = 0
        bad_model_inputs = 0
        bad_below15 = 0
        for i in range(len(model_arrays["raw_score"])):
            key = identity(model_arrays["candidate_id_hi"][i], model_arrays["candidate_id_lo"][i])
            model_sha = text(model_arrays["model_sha256"][i])
            definition_name = text(model_arrays["shower_definition_id"][i])
            if model_sha not in models or models[model_sha][1] != definition_name:
                bad_model_semantics += 1
                continue
            if text(model_arrays["shower_semantic_sha256"][i]) != semantic_sha(definition_name):
                bad_model_semantics += 1
            expected_input = view_lookup.get((key, definition_name))
            observed_input = np.asarray(model_arrays["ordered_input_witnesses"][i], dtype=np.float32)
            if expected_input is None or expected_input.shape != observed_input.shape or not np.allclose(
                expected_input, observed_input, rtol=0.0, atol=feature_tolerance, equal_nan=True
            ):
                bad_model_inputs += 1
            applicability = int(model_arrays["applicability_state"][i])
            if candidates[key]["et"] < 15.0:
                null_wp = all(math.isnan(float(model_arrays[name][i])) for name in ("wp70", "wp80", "wp90"))
                if applicability != 1 or not null_wp:
                    bad_below15 += 1
            if applicability == 0 and int(model_arrays["finite_score"][i]) != 0 and expected_input is not None:
                score_inputs[model_sha].append(expected_input)
                score_observed[model_sha].append(float(model_arrays["raw_score"][i]))
        report.update({
            "model_rows": int(len(model_arrays["raw_score"])),
            "bad_model_semantics": bad_model_semantics,
            "bad_model_inputs": bad_model_inputs,
            "bad_below15_rows": bad_below15,
        })
        if bad_model_semantics:
            failures.append("model_shower_semantic_mismatch")
        if bad_model_inputs:
            failures.append("model_input_view_mismatch")
        if bad_below15:
            failures.append("below15_model_domain_violation")

    try:
        import ROOT  # type: ignore
        score_report: dict[str, Any] = {}
        for model_sha, inputs in score_inputs.items():
            matrix = np.ascontiguousarray(np.asarray(inputs, dtype=np.float32))
            runtime = ROOT.TMVA.Experimental.RBDT("myBDT", str(models[model_sha][0]))
            recomputed = np.asarray(runtime.Compute(matrix), dtype=np.float64).reshape(-1)
            observed = np.asarray(score_observed[model_sha], dtype=np.float64)
            difference = np.abs(recomputed - observed)
            max_abs = float(np.max(difference)) if len(difference) else math.nan
            score_report[model_sha] = {"rows": len(observed), "max_abs_difference": max_abs}
            if not len(observed) or not math.isfinite(max_abs) or max_abs > score_tolerance:
                failures.append("runtime_tmva_score_mismatch")
        report["score_parity"] = score_report
    except Exception as exc:  # pragma: no cover - runtime-dependent
        report["score_parity_error"] = str(exc)
        failures.append("runtime_tmva_score_validation_failed")

    report["failures"] = sorted(set(failures))
    report["status"] = "PASS" if not report["failures"] else "FAIL"
    return report


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pair", action="append", nargs=5, metavar=("LABEL", "SYSTEM", "ACTIVE_DEF", "DIRECT", "WRITER"), required=True)
    parser.add_argument("--model", action="append", nargs=3, metavar=("SHA256", "TMVA_ROOT", "DEFINITION"), required=True)
    parser.add_argument("--expected-code", required=True)
    parser.add_argument("--expected-schema", required=True)
    parser.add_argument("--expected-semantic", required=True)
    parser.add_argument("--minimum-bytes", type=int, default=50_000)
    parser.add_argument("--feature-tolerance", type=float, default=2.0e-6)
    parser.add_argument("--score-tolerance", type=float, default=2.0e-7)
    parser.add_argument("--output-json", type=Path)
    args = parser.parse_args()
    model_map: dict[str, tuple[Path, str]] = {}
    for sha, path, definition_name in args.model:
        model_path = Path(path)
        actual = hashlib.sha256(model_path.read_bytes()).hexdigest()
        if actual != sha:
            raise SystemExit(f"model hash mismatch: {model_path}: {actual} != {sha}")
        if definition_name not in DEFINITIONS:
            raise SystemExit(f"unknown model shower definition: {definition_name}")
        model_map[sha] = (model_path, definition_name)
    metadata = {
        "code_sha256": args.expected_code,
        "schema_sha256": args.expected_schema,
        "semantic_sha256": args.expected_semantic,
    }
    pairs = []
    writers = []
    for label, system, active_definition, direct, writer in args.pair:
        if system not in ("pp", "auau") or active_definition not in DEFINITIONS:
            raise SystemExit(f"invalid pair contract: {label} {system} {active_definition}")
        pairs.append(compare_histograms(label, Path(direct), Path(writer)))
        writers.append(validate_writer(
            label, system, active_definition, Path(writer), metadata, model_map,
            args.minimum_bytes, args.feature_tolerance, args.score_tolerance,
        ))
    total_writer_bytes = sum(item["size_bytes"] for item in writers)
    report = {
        "schema": "THE134_SHOWER_FACTORIAL_REAL_REPLAY_VALIDATION_V1",
        "pairs": pairs, "writers": writers,
        "storage": {"writer_count": len(writers), "total_writer_bytes": total_writer_bytes},
    }
    report["status"] = "PASS" if (
        pairs and writers and all(item["status"] == "PASS" for item in pairs + writers)
    ) else "FAIL"
    payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output_json:
        args.output_json.parent.mkdir(parents=True, exist_ok=True)
        args.output_json.write_text(payload)
    print(payload, end="")
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
