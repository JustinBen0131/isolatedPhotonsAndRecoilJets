#!/usr/bin/env python3
"""Certify THE-119 direct/writer neutrality and deterministic jet response.

The validator is intentionally selection-neutral.  It checks the normalized
replay payload and the common scientific histograms without deriving or
changing any photon, isolation, sideband, or recoil selection.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import struct
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import uproot


TREE_NAMES = {
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
HASH_NAMES = (
    "schema_sha256",
    "semantic_sha256",
    "source_sha256",
    "model_sha256",
    "config_sha256",
    "code_sha256",
)
LINK_MATCH = 0
LINK_FAKE = 1
LINK_MISS = 2
LINK_MATCH_CANDIDATE = 5
TYPE_NONE = 0
TYPE_PHOTON = 1
TYPE_JET = 2


def _id(hi: Any, lo: Any) -> tuple[int, int]:
    return int(hi), int(lo)


def _float_bits(value: Any) -> int:
    return struct.unpack("<Q", struct.pack("<d", float(value)))[0]


def _digest_rows(rows: Iterable[tuple[Any, ...]]) -> str:
    digest = hashlib.sha256()
    for row in sorted(rows):
        digest.update(repr(row).encode("utf-8"))
        digest.update(b"\n")
    return digest.hexdigest()


def _metadata(replay: uproot.ReadOnlyDirectory) -> dict[str, str]:
    result: dict[str, str] = {}
    for name in HASH_NAMES:
        result[name] = str(replay[name].member("fTitle"))
    return result


def _source_sha(path: Path) -> str:
    with uproot.open(path) as root:
        return _metadata(root["ReplayFoundationV1"])["source_sha256"]


def _pyroot_health(path: Path) -> tuple[bool, str]:
    try:
        import ROOT  # type: ignore
    except Exception as exc:  # pragma: no cover - depends on runtime
        return False, f"pyroot_unavailable:{exc}"
    root_file = ROOT.TFile.Open(str(path), "READ")
    if not root_file:
        return False, "pyroot_open_returned_null"
    try:
        if root_file.IsZombie():
            return False, "root_is_zombie"
        if root_file.TestBit(ROOT.TFile.kRecovered):
            return False, "root_is_recovered"
        return True, "healthy"
    finally:
        root_file.Close()


def _histograms(path: Path) -> dict[str, dict[str, Any]]:
    histograms: dict[str, dict[str, Any]] = {}
    with uproot.open(path) as root:
        for raw_key, class_name in root.classnames(recursive=True).items():
            key = raw_key.split(";")[0]
            if key.startswith("ReplayFoundationV1/"):
                continue
            if not class_name.startswith(("TH1", "TH2", "TH3", "TProfile")):
                continue
            obj = root[key]
            values = np.asarray(obj.values(flow=True))
            variances_raw = obj.variances(flow=True)
            variances = None if variances_raw is None else np.asarray(variances_raw)
            try:
                sumw2 = np.asarray(obj.member("fSumw2"))
            except Exception:
                sumw2 = np.asarray([], dtype=np.float64)
            axes = [np.asarray(axis.edges()) for axis in obj.axes]
            histograms[key] = {
                "class": class_name,
                "values": values,
                "variances": variances,
                "sumw2": sumw2,
                "axes": axes,
            }
    return histograms


def compare_histograms(label: str, direct: Path, writer: Path) -> dict[str, Any]:
    result: dict[str, Any] = {"label": label, "status": "FAIL", "failures": []}
    direct_hists = _histograms(direct)
    writer_hists = _histograms(writer)
    direct_keys, writer_keys = set(direct_hists), set(writer_hists)
    result["histogram_count"] = len(direct_keys)
    result["direct_only"] = sorted(direct_keys - writer_keys)
    result["writer_only"] = sorted(writer_keys - direct_keys)
    if not direct_keys:
        result["failures"].append("no_scientific_histograms")
    if direct_keys != writer_keys:
        result["failures"].append("histogram_inventory_mismatch")
    mismatches: list[str] = []
    for key in sorted(direct_keys & writer_keys):
        left, right = direct_hists[key], writer_hists[key]
        if left["class"] != right["class"]:
            mismatches.append(f"{key}:class")
            continue
        for field in ("values", "sumw2"):
            if not np.array_equal(left[field], right[field], equal_nan=True):
                mismatches.append(f"{key}:{field}")
        if (left["variances"] is None) != (right["variances"] is None):
            mismatches.append(f"{key}:variances_presence")
        elif left["variances"] is not None and not np.array_equal(
            left["variances"], right["variances"], equal_nan=True
        ):
            mismatches.append(f"{key}:variances")
        if len(left["axes"]) != len(right["axes"]) or any(
            not np.array_equal(a, b, equal_nan=True)
            for a, b in zip(left["axes"], right["axes"], strict=False)
        ):
            mismatches.append(f"{key}:axes")
    result["mismatches"] = mismatches
    if mismatches:
        result["failures"].append("histogram_content_or_sumw2_mismatch")
    result["failures"] = sorted(set(result["failures"]))
    result["status"] = "PASS" if not result["failures"] else "FAIL"
    return result


def _read_jets(replay: uproot.ReadOnlyDirectory, truth: bool) -> dict[tuple[int, int], dict[str, Any]]:
    tree_name = "RJTruthJetV1" if truth else "RJJetV1"
    prefix = "truth_jet_id" if truth else "jet_id"
    pt_branch = "pt" if truth else "corrected_pt"
    arrays = replay[tree_name].arrays(
        [
            f"{prefix}_hi",
            f"{prefix}_lo",
            "event_id_hi",
            "event_id_lo",
            "radius",
            pt_branch,
            "eta",
            "phi",
        ],
        library="np",
    )
    result: dict[tuple[int, int], dict[str, Any]] = {}
    for index in range(len(arrays[pt_branch])):
        key = _id(arrays[f"{prefix}_hi"][index], arrays[f"{prefix}_lo"][index])
        if key in result:
            raise ValueError(f"duplicate {tree_name} identity {key}")
        result[key] = {
            "event": _id(arrays["event_id_hi"][index], arrays["event_id_lo"][index]),
            "radius": float(arrays["radius"][index]),
            "pt": float(arrays[pt_branch][index]),
            "eta": float(arrays["eta"][index]),
            "phi": float(arrays["phi"][index]),
        }
    return result


def _read_links(replay: uproot.ReadOnlyDirectory) -> list[dict[str, Any]]:
    branches = [
        "link_id_hi",
        "link_id_lo",
        "reco_type",
        "reco_id_hi",
        "reco_id_lo",
        "truth_type",
        "truth_id_hi",
        "truth_id_lo",
        "match_metric",
        "link_class",
    ]
    arrays = replay["RJRecoTruthLinkV1"].arrays(branches, library="np")
    return [
        {
            "id": _id(arrays["link_id_hi"][i], arrays["link_id_lo"][i]),
            "reco_type": int(arrays["reco_type"][i]),
            "reco": _id(arrays["reco_id_hi"][i], arrays["reco_id_lo"][i]),
            "truth_type": int(arrays["truth_type"][i]),
            "truth": _id(arrays["truth_id_hi"][i], arrays["truth_id_lo"][i]),
            "metric": float(arrays["match_metric"][i]),
            "class": int(arrays["link_class"][i]),
        }
        for i in range(len(arrays["link_class"]))
    ]


def _delta_r(left: dict[str, Any], right: dict[str, Any]) -> float:
    dphi = math.atan2(math.sin(left["phi"] - right["phi"]), math.cos(left["phi"] - right["phi"]))
    return math.hypot(left["eta"] - right["eta"], dphi)


def validate_writer(
    label: str,
    path: Path,
    expected: dict[str, str],
    minimum_bytes: int,
    require_all_classes: bool,
) -> dict[str, Any]:
    result: dict[str, Any] = {
        "label": label,
        "path": str(path),
        "size_bytes": path.stat().st_size if path.exists() else 0,
        "status": "FAIL",
        "failures": [],
    }
    if not path.exists() or path.stat().st_size < minimum_bytes:
        result["failures"].append("missing_or_tiny_root")
        return result
    healthy, health_detail = _pyroot_health(path)
    result["root_health"] = health_detail
    if not healthy:
        result["failures"].append(health_detail)
        return result

    with uproot.open(path) as root:
        if "ReplayFoundationV1" not in root:
            result["failures"].append("missing_replay_directory")
            return result
        replay = root["ReplayFoundationV1"]
        classnames = replay.classnames(recursive=False)
        trees = {key.split(";")[0] for key, value in classnames.items() if value == "TTree"}
        result["tree_count"] = len(trees)
        result["missing_trees"] = sorted(TREE_NAMES - trees)
        result["extra_trees"] = sorted(trees - TREE_NAMES)
        if trees != TREE_NAMES:
            result["failures"].append("tree_inventory_mismatch")

        metadata = _metadata(replay)
        result["metadata"] = metadata
        for name, value in metadata.items():
            if len(value) != 64 or any(c not in "0123456789abcdef" for c in value):
                result["failures"].append(f"invalid_{name}")
        for name, value in expected.items():
            if metadata.get(name) != value:
                result["failures"].append(f"{name}_mismatch")

        source_tree = replay["RJSourceOccurrenceV1"]
        source_manifest = source_tree["source_manifest_sha256"].array(library="np")
        normalized_source = {
            value.decode() if isinstance(value, bytes) else str(value) for value in source_manifest
        }
        result["source_manifest_values"] = sorted(normalized_source)
        if normalized_source != {metadata["source_sha256"]}:
            result["failures"].append("source_manifest_metadata_mismatch")

        weights = replay["RJWeightComponentV1"].arrays(
            ["application_count", "final_weight"], library="np"
        )
        bad_weight_count = int(np.count_nonzero(weights["application_count"] != 1))
        nonfinite_weights = int(np.count_nonzero(~np.isfinite(weights["final_weight"])))
        result["weight_rows"] = int(len(weights["application_count"]))
        result["bad_application_count"] = bad_weight_count
        result["nonfinite_weights"] = nonfinite_weights
        if bad_weight_count:
            result["failures"].append("weight_not_applied_exactly_once")
        if nonfinite_weights:
            result["failures"].append("nonfinite_final_weight")

        candidates = replay["RJPhotonCandidateV1"].arrays(
            ["candidate_id_hi", "candidate_id_lo", "cluster_et"], library="np"
        )
        below_ids = {
            _id(hi, lo)
            for hi, lo, et in zip(
                candidates["candidate_id_hi"],
                candidates["candidate_id_lo"],
                candidates["cluster_et"],
                strict=True,
            )
            if float(et) < 15.0
        }
        evaluations = replay["RJModelEvaluationV1"].arrays(
            [
                "candidate_id_hi",
                "candidate_id_lo",
                "applicability_state",
                "wp70",
                "wp80",
                "wp90",
                "delta_wp70",
                "delta_wp80",
                "delta_wp90",
            ],
            library="np",
        )
        below_eval_rows = 0
        bad_below_rows = 0
        for i in range(len(evaluations["applicability_state"])):
            if _id(evaluations["candidate_id_hi"][i], evaluations["candidate_id_lo"][i]) not in below_ids:
                continue
            below_eval_rows += 1
            null_wp = all(
                math.isnan(float(evaluations[name][i]))
                for name in ("wp70", "wp80", "wp90", "delta_wp70", "delta_wp80", "delta_wp90")
            )
            if int(evaluations["applicability_state"][i]) != 1 or not null_wp:
                bad_below_rows += 1
        result["below15_candidate_count"] = len(below_ids)
        result["below15_evaluation_rows"] = below_eval_rows
        result["bad_below15_evaluation_rows"] = bad_below_rows
        if bad_below_rows:
            result["failures"].append("below15_model_domain_violation")

        reco = _read_jets(replay, truth=False)
        truth = _read_jets(replay, truth=True)
        links = _read_links(replay)
        jet_links = [row for row in links if row["reco_type"] == TYPE_JET or row["truth_type"] == TYPE_JET]
        photon_links = [row for row in links if row["reco_type"] == TYPE_PHOTON or row["truth_type"] == TYPE_PHOTON]
        result["reco_jet_count"] = len(reco)
        result["truth_jet_count"] = len(truth)
        result["jet_link_count"] = len(jet_links)
        result["photon_link_count"] = len(photon_links)
        result["photon_link_digest"] = _digest_rows(
            (
                row["reco_type"], row["reco"], row["truth_type"], row["truth"],
                row["class"], _float_bits(row["metric"]),
            )
            for row in photon_links
        )

        duplicate_link_ids = len(links) - len({row["id"] for row in links})
        result["duplicate_link_ids"] = duplicate_link_ids
        if duplicate_link_ids:
            result["failures"].append("duplicate_link_identities")

        candidate_ids = {
            _id(hi, lo)
            for hi, lo in zip(
                candidates["candidate_id_hi"], candidates["candidate_id_lo"], strict=True
            )
        }
        truth_photons = replay["RJTruthPhotonV1"].arrays(
            ["truth_photon_id_hi", "truth_photon_id_lo"], library="np"
        )
        truth_photon_ids = {
            _id(hi, lo)
            for hi, lo in zip(
                truth_photons["truth_photon_id_hi"],
                truth_photons["truth_photon_id_lo"],
                strict=True,
            )
        }
        photon_reco_final = [
            row for row in photon_links
            if row["reco_type"] == TYPE_PHOTON and row["class"] in (LINK_MATCH, LINK_FAKE)
        ]
        photon_truth_match = {
            row["truth"] for row in photon_links
            if row["truth_type"] == TYPE_PHOTON and row["class"] == LINK_MATCH
        }
        photon_truth_miss = [
            row["truth"] for row in photon_links
            if row["truth_type"] == TYPE_PHOTON and row["class"] == LINK_MISS
        ]
        if len(photon_reco_final) != len(candidate_ids) or {
            row["reco"] for row in photon_reco_final
        } != candidate_ids:
            result["failures"].append("photon_reco_link_partition_mismatch")
        if len(photon_truth_miss) != len(set(photon_truth_miss)):
            result["failures"].append("duplicate_photon_truth_miss")
        if photon_truth_match | set(photon_truth_miss) != truth_photon_ids:
            result["failures"].append("photon_truth_link_partition_mismatch")
        if photon_truth_match & set(photon_truth_miss):
            result["failures"].append("photon_truth_match_miss_overlap")

        by_class = {code: [row for row in jet_links if row["class"] == code] for code in (0, 1, 2, 5)}
        result["jet_link_class_counts"] = {str(code): len(rows) for code, rows in by_class.items()}
        if require_all_classes and any(not by_class[code] for code in (0, 1, 2, 5)):
            result["failures"].append("missing_required_jet_link_class_witness")

        reco_groups: dict[tuple[tuple[int, int], int], list[tuple[int, int]]] = defaultdict(list)
        truth_groups: dict[tuple[tuple[int, int], int], list[tuple[int, int]]] = defaultdict(list)
        for key, row in reco.items():
            reco_groups[(row["event"], int(round(row["radius"] * 1_000_000)))].append(key)
        for key, row in truth.items():
            truth_groups[(row["event"], int(round(row["radius"] * 1_000_000)))].append(key)

        expected_edges: dict[tuple[tuple[int, int], tuple[int, int]], float] = {}
        for group, reco_ids in reco_groups.items():
            for reco_id in reco_ids:
                for truth_id in truth_groups.get(group, []):
                    dr = _delta_r(reco[reco_id], truth[truth_id])
                    if dr < 0.3:
                        expected_edges[(reco_id, truth_id)] = dr

        candidate_edges: dict[tuple[tuple[int, int], tuple[int, int]], float] = {}
        malformed_jet_links = 0
        for row in by_class[LINK_MATCH_CANDIDATE]:
            if row["reco_type"] != TYPE_JET or row["truth_type"] != TYPE_JET:
                malformed_jet_links += 1
                continue
            key = (row["reco"], row["truth"])
            if key in candidate_edges:
                malformed_jet_links += 1
            candidate_edges[key] = row["metric"]
        result["malformed_or_duplicate_jet_links"] = malformed_jet_links
        if malformed_jet_links:
            result["failures"].append("malformed_or_duplicate_jet_links")
        if set(candidate_edges) != set(expected_edges):
            result["failures"].append("match_candidate_edge_set_mismatch")
        metric_mismatches = sum(
            not math.isclose(candidate_edges[key], expected_edges[key], rel_tol=0.0, abs_tol=1.0e-12)
            for key in set(candidate_edges) & set(expected_edges)
        )
        result["candidate_metric_mismatches"] = metric_mismatches
        if metric_mismatches:
            result["failures"].append("match_candidate_metric_mismatch")

        ordered_edges = sorted(
            expected_edges,
            key=lambda key: (
                expected_edges[key],
                -reco[key[0]]["pt"],
                key[0][0], key[0][1], key[1][0], key[1][1],
            ),
        )
        selected_reco: set[tuple[int, int]] = set()
        selected_truth: set[tuple[int, int]] = set()
        expected_matches: set[tuple[tuple[int, int], tuple[int, int]]] = set()
        for reco_id, truth_id in ordered_edges:
            if reco_id in selected_reco or truth_id in selected_truth:
                continue
            selected_reco.add(reco_id)
            selected_truth.add(truth_id)
            expected_matches.add((reco_id, truth_id))

        actual_matches = {(row["reco"], row["truth"]) for row in by_class[LINK_MATCH]}
        if actual_matches != expected_matches:
            result["failures"].append("deterministic_match_selection_mismatch")
        fake_ids = [row["reco"] for row in by_class[LINK_FAKE]]
        miss_ids = [row["truth"] for row in by_class[LINK_MISS]]
        if len(fake_ids) != len(set(fake_ids)) or len(miss_ids) != len(set(miss_ids)):
            result["failures"].append("duplicate_selected_jet_identity")
        if set(fake_ids) != set(reco) - selected_reco:
            result["failures"].append("reco_fake_partition_mismatch")
        if set(miss_ids) != set(truth) - selected_truth:
            result["failures"].append("truth_miss_partition_mismatch")

        jet_quality = replay["RJJetV1"]["quality_bitmask"].array(library="np")
        retained_by_bit = {
            key for key, bitmask in zip(reco, jet_quality, strict=True) if int(bitmask) & 1
        }
        expected_retained = {key for key, row in reco.items() if row["pt"] >= 5.0}
        if retained_by_bit != expected_retained:
            result["failures"].append("jet_constituent_ownership_bit_mismatch")
        constituents = replay["RJJetConstituentV1"].arrays(
            ["jet_id_hi", "jet_id_lo"], library="np"
        )
        constituent_jet_ids = {
            _id(hi, lo)
            for hi, lo in zip(constituents["jet_id_hi"], constituents["jet_id_lo"], strict=True)
        }
        below_boundary = sorted(
            key for key in constituent_jet_ids if key not in reco or reco[key]["pt"] < 5.0
        )
        result["constituent_jet_count"] = len(constituent_jet_ids)
        result["below_boundary_constituent_jet_count"] = len(below_boundary)
        if below_boundary:
            result["failures"].append("jet_constituent_storage_boundary_violation")

    result["failures"] = sorted(set(result["failures"]))
    result["status"] = "PASS" if not result["failures"] else "FAIL"
    return result


def photon_link_digest(path: Path) -> str:
    with uproot.open(path) as root:
        links = _read_links(root["ReplayFoundationV1"])
    return _digest_rows(
        (
            row["reco_type"], row["reco"], row["truth_type"], row["truth"],
            row["class"], _float_bits(row["metric"]),
        )
        for row in links
        if row["reco_type"] == TYPE_PHOTON or row["truth_type"] == TYPE_PHOTON
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pair", action="append", nargs=3, metavar=("LABEL", "DIRECT", "WRITER"), required=True)
    parser.add_argument(
        "--writer",
        action="append",
        nargs=5,
        metavar=("LABEL", "PATH", "CONFIG_SHA", "MODEL_SHA", "SOURCE_SHA"),
        required=True,
    )
    parser.add_argument("--photon-baseline", action="append", nargs=3, metavar=("LABEL", "WRITER", "BASELINE"))
    parser.add_argument("--expected-code", required=True)
    parser.add_argument("--expected-schema", required=True)
    parser.add_argument("--expected-semantic", required=True)
    parser.add_argument("--minimum-bytes", type=int, default=50_000)
    parser.add_argument("--require-all-link-classes", action="store_true")
    parser.add_argument("--output-json", type=Path)
    args = parser.parse_args()

    pairs = [compare_histograms(label, Path(direct), Path(writer)) for label, direct, writer in args.pair]
    writers = []
    for label, path, config_sha, model_sha, source_sha in args.writer:
        expected = {
            "code_sha256": args.expected_code,
            "schema_sha256": args.expected_schema,
            "semantic_sha256": args.expected_semantic,
            "config_sha256": config_sha,
            "model_sha256": model_sha,
            "source_sha256": source_sha,
        }
        writers.append(
            validate_writer(label, Path(path), expected, args.minimum_bytes, args.require_all_link_classes)
        )

    baselines = []
    for label, writer, baseline in args.photon_baseline or []:
        writer_digest = photon_link_digest(Path(writer))
        baseline_digest = photon_link_digest(Path(baseline))
        writer_source = _source_sha(Path(writer))
        baseline_source = _source_sha(Path(baseline))
        comparable = writer_source == baseline_source
        baselines.append(
            {
                "label": label,
                "writer": writer,
                "baseline": baseline,
                "writer_digest": writer_digest,
                "baseline_digest": baseline_digest,
                "writer_source_sha256": writer_source,
                "baseline_source_sha256": baseline_source,
                "status": (
                    "PASS" if comparable and writer_digest == baseline_digest
                    else "NOT_COMPARABLE_SOURCE_MISMATCH" if not comparable
                    else "FAIL"
                ),
            }
        )

    report = {
        "schema": "THE119_JET_RESPONSE_VALIDATION_V1",
        "pairs": pairs,
        "writers": writers,
        "photon_link_baselines": baselines,
    }
    report["status"] = "PASS" if (
        pairs
        and writers
        and all(item["status"] == "PASS" for item in pairs + writers)
        and all(item["status"] in ("PASS", "NOT_COMPARABLE_SOURCE_MISMATCH") for item in baselines)
    ) else "FAIL"
    payload = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output_json:
        args.output_json.parent.mkdir(parents=True, exist_ok=True)
        args.output_json.write_text(payload)
    print(payload, end="")
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
