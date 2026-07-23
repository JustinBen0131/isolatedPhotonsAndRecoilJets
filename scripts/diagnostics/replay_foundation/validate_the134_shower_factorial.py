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

FROZEN_MINIMUM_WRITER_BYTES = 50_000
# Retained feature vectors and the independent reconstruction are both
# float32 PhotonClusterBuilder arithmetic.  Their contract is exact; unlike
# TMVA scoring and Au+Au isolation, no nonzero feature tolerance is authorized.
FROZEN_MAX_FEATURE_TOLERANCE = 0.0
FROZEN_MAX_SCORE_TOLERANCE = 2.0e-7


def valid_sha256(value: str) -> bool:
    return len(value) == 64 and all(character in "0123456789abcdef" for character in value)


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
        "|numeric=float32_PhotonClusterBuilder_row_major"
    )


def semantic_sha(name: str) -> str:
    return hashlib.sha256(semantic_text(name).encode()).hexdigest()


def close(left: float, right: float, tolerance: float) -> bool:
    if math.isnan(left) and math.isnan(right):
        return True
    return math.isclose(left, right, rel_tol=0.0, abs_tol=tolerance)


def _f32_add(left: np.float32, right: np.float32) -> np.float32:
    """Spell out one C++ ``float`` addition for the independent oracle."""

    return np.float32(np.float32(left) + np.float32(right))


def _f32_sub(left: np.float32, right: np.float32) -> np.float32:
    """Spell out one C++ ``float`` subtraction for the independent oracle."""

    return np.float32(np.float32(left) - np.float32(right))


def native_rawcluster_shape_from_cells(
    floor: float, cells: list[dict[str, Any]]
) -> dict[str, Any]:
    """Replay ``RawClusterv1::get_shower_shapes`` from retained map cells.

    The source implementation iterates the ordered ``RawCluster::TowerMap``
    and performs every accumulation in ``float``.  The replay rows carry
    absolute eta/phi indices and the original map value, so sorting by those
    indices reproduces the RawTowerDefs key order for CEMC.
    """

    owned = sorted(
        (
            cell
            for cell in cells
            if int(cell["rawcluster_owned"]) != 0
            and int(cell["rawcluster_value_present"]) != 0
        ),
        key=lambda cell: (int(cell["tower_eta_index"]), int(cell["tower_phi_index"])),
    )
    maximum_energy = np.float32(0.0)
    maximum_eta: int | None = None
    maximum_phi: int | None = None
    for cell in owned:
        energy = np.float32(cell["rawcluster_map_value"])
        if energy > maximum_energy:
            maximum_energy = energy
            maximum_eta = int(cell["tower_eta_index"])
            maximum_phi = int(cell["tower_phi_index"])
    if maximum_eta is None or maximum_phi is None:
        return {"valid": False}

    threshold = np.float32(floor)
    total_energy = np.float32(0.0)
    delta_eta_numerator = np.float32(0.0)
    delta_phi_numerator = np.float32(0.0)

    def wrapped_delta_phi(tower_phi: int, reference_phi: int) -> np.float32:
        delta = np.float32(tower_phi - reference_phi)
        wrapped = np.float32(256.0) - np.abs(delta)
        if np.abs(wrapped) < np.abs(delta):
            return np.float32(-wrapped if delta > 0 else wrapped)
        return delta

    for cell in owned:
        energy = np.float32(cell["rawcluster_map_value"])
        if not energy > threshold:
            continue
        delta_eta = np.float32(int(cell["tower_eta_index"]) - maximum_eta)
        delta_phi = wrapped_delta_phi(int(cell["tower_phi_index"]), maximum_phi)
        total_energy = _f32_add(total_energy, energy)
        delta_eta_numerator = _f32_add(
            delta_eta_numerator, np.float32(energy * delta_eta)
        )
        delta_phi_numerator = _f32_add(
            delta_phi_numerator, np.float32(energy * delta_phi)
        )
    if not np.isfinite(total_energy) or not total_energy > np.float32(0.0):
        return {"valid": False}

    delta_eta_mean = np.float32(delta_eta_numerator / total_energy)
    delta_phi_mean = np.float32(delta_phi_numerator / total_energy)
    local_center_eta = math.floor(float(_f32_add(delta_eta_mean, np.float32(0.5))))
    local_center_phi = math.floor(float(_f32_add(delta_phi_mean, np.float32(0.5))))
    eta_shift = -1 if _f32_sub(delta_eta_mean, np.float32(local_center_eta)) < 0 else 1
    phi_shift = -1 if _f32_sub(delta_phi_mean, np.float32(local_center_phi)) < 0 else 1
    eta1 = maximum_eta + local_center_eta
    eta2 = eta1 + eta_shift
    phi1 = (maximum_phi + local_center_phi) % 256
    phi2 = (maximum_phi + local_center_phi + phi_shift) % 256
    energy_by_coordinate = {
        (int(cell["tower_eta_index"]), int(cell["tower_phi_index"])): np.float32(
            cell["rawcluster_map_value"]
        )
        for cell in owned
    }

    def selected_energy(eta_index: int, phi_index: int) -> np.float32:
        energy = energy_by_coordinate.get((eta_index, phi_index), np.float32(0.0))
        return energy if energy > threshold else np.float32(0.0)

    e1 = selected_energy(eta1, phi1)
    e2 = selected_energy(eta1, phi2)
    e3 = selected_energy(eta2, phi2)
    e4 = selected_energy(eta2, phi1)
    numerator1 = _f32_add(_f32_add(_f32_add(e1, e2), e3), e4)
    numerator2 = _f32_sub(_f32_sub(_f32_add(e1, e2), e3), e4)
    numerator3 = _f32_add(_f32_sub(_f32_sub(e1, e2), e3), e4)
    raw_shape_eta = _f32_add(delta_eta_mean, np.float32(maximum_eta))
    raw_shape_phi = _f32_add(delta_phi_mean, np.float32(maximum_phi))
    raw_center_eta = float(raw_shape_eta) + 0.5
    raw_center_phi = float(raw_shape_phi) + 0.5
    return {
        "valid": True,
        "raw_center_eta": raw_center_eta,
        "raw_center_phi": raw_center_phi,
        "center_eta_index": int(math.floor(raw_center_eta)),
        "center_phi_index": int(math.floor(raw_center_phi)) % 256,
        "native_et1": np.float32(numerator1 / total_energy),
        "native_et2": np.float32(numerator2 / total_energy),
        "native_et3": np.float32(numerator3 / total_energy),
        "native_et4": np.float32(e3 / total_energy),
    }


def native_et_replay_mismatch(
    row: dict[str, Any], rebuilt: dict[str, Any], tolerance: float
) -> bool:
    """Return true when stored native shape values drift from cell replay."""

    return (not rebuilt.get("valid", False)) or any(
        not close(float(row[f"native_et{index}"]), float(rebuilt[f"native_et{index}"]), tolerance)
        for index in range(1, 5)
    )


def valid_grid_or_owned_provenance(cell: dict[str, Any]) -> bool:
    """Accept 7x7 members or complete owned-map provenance outside the grid."""

    bitmask = int(cell["grid_membership_bitmask"])
    if bitmask < 0 or bitmask > 3:
        return False
    if bitmask != 0:
        return True
    return (
        int(cell["rawcluster_owned"]) != 0
        and int(cell["rawcluster_value_present"]) != 0
    )


def view_from_cells(name: str, cells: list[dict[str, Any]]) -> dict[str, Any]:
    energy_source, sum_membership, moment_membership, floor = DEFINITIONS[name]
    native_shape = native_rawcluster_shape_from_cells(floor, cells)
    if not native_shape.get("valid", False):
        return native_shape
    center_eta = int(native_shape["center_eta_index"])
    center_phi = int(native_shape["center_phi_index"])
    raw_center_eta = np.float32(native_shape["raw_center_eta"])
    raw_center_phi = np.float32(native_shape["raw_center_phi"])
    cog_eta = np.float32(3.0) + np.float32(
        raw_center_eta - np.floor(raw_center_eta) - np.float32(0.5)
    )
    cog_phi = np.float32(3.0) + np.float32(
        raw_center_phi - np.floor(raw_center_phi) - np.float32(0.5)
    )
    sign_phi = 1 if cog_phi > 3.0 else -1
    out = {
        **native_shape,
        "e11": np.float32(0.0), "e33": np.float32(0.0),
        "e32": np.float32(0.0), "e35": np.float32(0.0),
        "moment_eta_numerator": np.float32(0.0),
        "moment_phi_numerator": np.float32(0.0),
        "moment_denominator": np.float32(0.0),
        "moment33_eta_numerator": np.float32(0.0),
        "moment33_phi_numerator": np.float32(0.0),
        "moment33_denominator": np.float32(0.0),
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
            energy = np.float32(cell["calibrated_energy"])
            selected = int(cell["is_good"]) != 0 and np.isfinite(energy) and energy > np.float32(floor)
        else:
            energy = np.float32(cell["rawcluster_map_value"])
            selected = (
                int(cell["rawcluster_owned"]) != 0
                and int(cell["rawcluster_value_present"]) != 0
                and np.isfinite(energy) and energy > np.float32(floor)
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
                out["e11"] = np.float32(out["e11"] + energy)
            if di <= 1 and dj <= 1:
                out["e33"] = np.float32(out["e33"] + energy)
            if di <= 1 and (j == 3 or j == 3 + sign_phi):
                out["e32"] = np.float32(out["e32"] + energy)
            if di <= 1 and dj <= 2:
                out["e35"] = np.float32(out["e35"] + energy)
        if moment_member:
            out["active_moment_cell_count"] += 1
            deta = np.float32(i) - cog_eta
            dphi = np.float32(j) - cog_phi
            out["moment_denominator"] = np.float32(out["moment_denominator"] + energy)
            if i != 3 or j != 3:
                out["moment_eta_numerator"] = np.float32(
                    out["moment_eta_numerator"] + np.float32(np.float32(energy * deta) * deta)
                )
                out["moment_phi_numerator"] = np.float32(
                    out["moment_phi_numerator"] + np.float32(np.float32(energy * dphi) * dphi)
                )
            if di <= 1 and dj <= 1:
                out["moment33_denominator"] = np.float32(out["moment33_denominator"] + energy)
                if i != 3 or j != 3:
                    out["moment33_eta_numerator"] = np.float32(
                        out["moment33_eta_numerator"] + np.float32(np.float32(energy * deta) * deta)
                    )
                    out["moment33_phi_numerator"] = np.float32(
                        out["moment33_phi_numerator"] + np.float32(np.float32(energy * dphi) * dphi)
                    )
    out["e11_over_e33"] = np.float32(out["e11"] / out["e33"]) if out["e33"] > 0 else math.nan
    out["e32_over_e35"] = np.float32(out["e32"] / out["e35"]) if out["e35"] > 0 else math.nan
    out["weta_cogx"] = np.float32(out["moment_eta_numerator"] / out["moment_denominator"]) if out["moment_denominator"] > 0 else math.nan
    out["wphi_cogx"] = np.float32(out["moment_phi_numerator"] / out["moment_denominator"]) if out["moment_denominator"] > 0 else math.nan
    out["weta33_cogx"] = np.float32(out["moment33_eta_numerator"] / out["moment33_denominator"]) if out["moment33_denominator"] > 0 else math.nan
    out["wphi33_cogx"] = np.float32(out["moment33_phi_numerator"] / out["moment33_denominator"]) if out["moment33_denominator"] > 0 else math.nan
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
    expected_model_shas: set[str],
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
        events: dict[tuple[int, int], tuple[float, float]] = {}
        for i in range(len(event_arrays["vertex_z"])):
            key = identity(event_arrays["event_id_hi"][i], event_arrays["event_id_lo"][i])
            if key in events:
                failures.append("duplicate_event_identity")
                continue
            events[key] = (
                float(event_arrays["vertex_z"][i]), float(event_arrays["centrality"][i])
            )
        candidate_arrays = replay["RJPhotonCandidateV1"].arrays(
            ["candidate_id_hi", "candidate_id_lo", "event_id_hi", "event_id_lo",
             "cluster_et", "eta", "ordered_features", "shower_definition_views"], library="np"
        )
        candidates: dict[tuple[int, int], dict[str, Any]] = {}
        for i in range(len(candidate_arrays["cluster_et"])):
            key = identity(candidate_arrays["candidate_id_hi"][i], candidate_arrays["candidate_id_lo"][i])
            if key in candidates:
                failures.append("duplicate_candidate_identity")
                continue
            candidate_view_names = [
                text(value) for value in candidate_arrays["shower_definition_views"][i]
            ]
            if (
                len(candidate_view_names) != len(DEFINITIONS)
                or len(set(candidate_view_names)) != len(candidate_view_names)
                or set(candidate_view_names) != set(DEFINITIONS)
            ):
                failures.append("candidate_declared_view_inventory_mismatch")
            candidates[key] = {
                "event": identity(candidate_arrays["event_id_hi"][i], candidate_arrays["event_id_lo"][i]),
                "et": float(candidate_arrays["cluster_et"][i]),
                "eta": float(candidate_arrays["eta"][i]),
                "features": np.asarray(candidate_arrays["ordered_features"][i], dtype=np.float32),
                "views": set(candidate_view_names),
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
        orphan_cell_candidates = set(cell_rows).difference(candidates)
        orphan_view_candidates = set(view_rows).difference(candidates)
        if orphan_cell_candidates:
            failures.append("shower_cell_candidate_foreign_key_missing")
        if orphan_view_candidates:
            failures.append("shower_view_candidate_foreign_key_missing")
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
            names_in_order = [text(row["definition_name"]) for row in candidate_views]
            names_for_candidate = set(names_in_order)
            if (
                candidate["views"] != set(DEFINITIONS)
                or names_for_candidate != set(DEFINITIONS)
                or len(candidate_views) != len(DEFINITIONS)
                or len(names_in_order) != len(names_for_candidate)
            ):
                bad_view_sets += 1
            coordinate_keys = [
                (int(row["tower_eta_index"]), int(row["tower_phi_index"])) for row in cell_rows.get(key, [])
            ]
            if len(coordinate_keys) != len(set(coordinate_keys)):
                failures.append("duplicate_absolute_tower_identity")
            if any(not (0 <= eta < 96 and 0 <= phi < 256) for eta, phi in coordinate_keys):
                failures.append("absolute_tower_coordinate_out_of_range")
            if any(
                not valid_grid_or_owned_provenance(row)
                for row in cell_rows.get(key, [])
            ):
                failures.append("invalid_grid_or_owned_provenance_membership")
            if candidate["event"] not in events:
                failures.append("candidate_event_foreign_key_missing")
                continue
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
                    or not close(float(row["floor_gev"]), expected_floor, 0.0)
                ):
                    replay_mismatches += 1
                rebuilt = view_from_cells(name, cell_rows.get(key, []))
                if not rebuilt.get("valid", False):
                    replay_mismatches += 1
                    continue
                replay_mismatches += int(
                    not close(
                        float(row["raw_center_eta"]),
                        float(rebuilt["raw_center_eta"]),
                        feature_tolerance,
                    )
                )
                replay_mismatches += int(
                    not close(
                        float(row["raw_center_phi"]),
                        float(rebuilt["raw_center_phi"]),
                        feature_tolerance,
                    )
                )
                replay_mismatches += int(
                    int(row["center_eta_index"]) != int(rebuilt["center_eta_index"])
                )
                replay_mismatches += int(
                    int(row["center_phi_index"]) != int(rebuilt["center_phi_index"])
                )
                replay_mismatches += int(
                    native_et_replay_mismatch(row, rebuilt, feature_tolerance)
                )
                replay_mismatches += sum(
                    not close(float(row[field]), float(rebuilt[field]), 0.0) for field in scalar_fields
                )
                replay_mismatches += sum(int(row[field]) != int(rebuilt[field]) for field in count_fields)
                native = [float(rebuilt[f"native_et{i}"]) for i in range(1, 5)]
                if system == "auau":
                    expected_features = [
                        candidate["et"], rebuilt["weta_cogx"], rebuilt["wphi_cogx"],
                        rebuilt["weta33_cogx"], rebuilt["wphi33_cogx"], vertex_z,
                        candidate["eta"], rebuilt["e11_over_e33"], *native,
                        rebuilt["e32_over_e35"], centrality,
                    ]
                else:
                    expected_features = [
                        candidate["et"], rebuilt["weta_cogx"], rebuilt["wphi_cogx"], vertex_z,
                        candidate["eta"], rebuilt["e11_over_e33"], *native,
                        rebuilt["e32_over_e35"],
                    ]
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
        observed_declared_models: set[str] = set()
        bad_model_semantics = 0
        bad_model_inputs = 0
        bad_below15 = 0
        for i in range(len(model_arrays["raw_score"])):
            key = identity(model_arrays["candidate_id_hi"][i], model_arrays["candidate_id_lo"][i])
            model_sha = text(model_arrays["model_sha256"][i])
            definition_name = text(model_arrays["shower_definition_id"][i])
            if key not in candidates:
                failures.append("model_candidate_foreign_key_missing")
                continue
            if model_sha not in models or models[model_sha][1] != definition_name:
                bad_model_semantics += 1
                continue
            observed_declared_models.add(model_sha)
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
        report["expected_model_sha256s"] = sorted(expected_model_shas)
        report["observed_model_sha256s"] = sorted(observed_declared_models)
        if observed_declared_models != expected_model_shas:
            failures.append("system_model_sha_inventory_mismatch")
        for model_sha in sorted(expected_model_shas):
            inputs = score_inputs.get(model_sha, [])
            matrix = np.ascontiguousarray(np.asarray(inputs, dtype=np.float32))
            observed = np.asarray(score_observed[model_sha], dtype=np.float64)
            if not len(inputs) or not len(observed):
                score_report[model_sha] = {"rows": 0, "max_abs_difference": math.nan}
                failures.append("runtime_tmva_zero_parity_rows")
                continue
            runtime = ROOT.TMVA.Experimental.RBDT("myBDT", str(models[model_sha][0]))
            recomputed = np.asarray(runtime.Compute(matrix), dtype=np.float64).reshape(-1)
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
    parser.add_argument(
        "--model", action="append", nargs=4,
        metavar=("SYSTEM", "SHA256", "TMVA_ROOT", "DEFINITION"), required=True,
        help="declare each exact system-scoped runtime model expected in its writer",
    )
    parser.add_argument("--expected-code", required=True)
    parser.add_argument("--expected-schema", required=True)
    parser.add_argument("--expected-semantic", required=True)
    parser.add_argument("--minimum-bytes", type=int, default=50_000)
    parser.add_argument("--feature-tolerance", type=float, default=0.0)
    parser.add_argument("--score-tolerance", type=float, default=2.0e-7)
    parser.add_argument("--output-json", type=Path)
    args = parser.parse_args()
    for label, value in (
        ("--expected-code", args.expected_code),
        ("--expected-schema", args.expected_schema),
        ("--expected-semantic", args.expected_semantic),
    ):
        if not valid_sha256(value):
            raise SystemExit(f"{label} must be an exact lowercase SHA-256")
    if args.minimum_bytes < FROZEN_MINIMUM_WRITER_BYTES:
        raise SystemExit(
            f"--minimum-bytes cannot be lower than the frozen {FROZEN_MINIMUM_WRITER_BYTES}-byte gate"
        )
    if not 0.0 <= args.feature_tolerance <= FROZEN_MAX_FEATURE_TOLERANCE:
        raise SystemExit(
            "--feature-tolerance must be nonnegative and cannot exceed the frozen "
            f"{FROZEN_MAX_FEATURE_TOLERANCE:.1e} gate"
        )
    if not 0.0 <= args.score_tolerance <= FROZEN_MAX_SCORE_TOLERANCE:
        raise SystemExit(
            "--score-tolerance must be nonnegative and cannot exceed the frozen "
            f"{FROZEN_MAX_SCORE_TOLERANCE:.1e} gate"
        )
    model_map: dict[str, tuple[Path, str]] = {}
    system_model_shas: dict[str, set[str]] = {"pp": set(), "auau": set()}
    for system, sha, path, definition_name in args.model:
        if system not in system_model_shas:
            raise SystemExit(f"unknown model system: {system}")
        if not valid_sha256(sha):
            raise SystemExit(f"model SHA-256 is malformed: {sha}")
        model_path = Path(path)
        actual = hashlib.sha256(model_path.read_bytes()).hexdigest()
        if actual != sha:
            raise SystemExit(f"model hash mismatch: {model_path}: {actual} != {sha}")
        if definition_name not in DEFINITIONS:
            raise SystemExit(f"unknown model shower definition: {definition_name}")
        if sha in model_map:
            raise SystemExit(f"duplicate declared model SHA-256: {sha}")
        model_map[sha] = (model_path, definition_name)
        system_model_shas[system].add(sha)
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
        if not system_model_shas[system]:
            raise SystemExit(f"no exact runtime model declared for pair system: {system}")
        pairs.append(compare_histograms(label, Path(direct), Path(writer)))
        writers.append(validate_writer(
            label, system, active_definition, Path(writer), metadata, model_map,
            system_model_shas[system],
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
