#!/usr/bin/env python3
"""Build the fail-closed THE-255 AuAu H70 offline presentation package.

The builder never reads a data TTree or DST and never contacts SDCC.  When the
final AuAu membership gate is open it writes explicit status-only DATA rungs,
while materializing the already-accepted nominal H70 response, leakage and K
support, response-only closure diagnostics, and the historical comparison.
"""

from __future__ import annotations

import argparse
from array import array
import csv
from datetime import datetime, timezone
import hashlib
import io
import json
import math
import os
from pathlib import Path
import shutil
import sys
from typing import Any, Iterable, Mapping, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import ROOT


REPO = Path(__file__).resolve().parents[3]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.data_prep.recoiljets.the243_h70_corrections import (  # noqa: E402
    PT_EDGES,
    build_response_bundle,
    chi2_ndf,
    iterative_bayes,
)


SCHEMA = "THE255AuAuFullMoneyPlotPresentationPackageV1"
STATUS = "FOUNDATION_RUNNING"
MODEL_SHA256 = "8328af75235c8bcbf63659ba791ceff4e1165018de269c889472c6a80a26a41f"
RESPONSE_CATALOG_SHA256 = "99b67237e768871924128bb6fee0073cc0d20080f4ac406e17e9136dc6631922"
THRESHOLDS = (5, 7, 10, 12)
FINE_XJ_EDGES = np.asarray(
    [0.0, 0.30, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.10, 1.25, 1.50, 2.00, 3.00],
    dtype=float,
)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds").replace("+00:00", "Z")


def load(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"JSON object required: {path}")
    return value


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_text(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o640)
    with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
        stream.write(value)


def write_json(path: Path, value: Mapping[str, Any]) -> None:
    write_text(path, json.dumps(value, indent=2, sort_keys=True) + "\n")


def csv_payload(fieldnames: Sequence[str], rows: Iterable[Mapping[str, Any]]) -> str:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=fieldnames, lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue()


def write_csv(path: Path, fieldnames: Sequence[str], rows: Iterable[Mapping[str, Any]]) -> None:
    write_text(path, csv_payload(fieldnames, rows))


def set_read_only(path: Path) -> None:
    path.chmod(0o440)


def save_figure(fig: Any, stem: Path) -> None:
    for suffix in (".png", ".pdf"):
        destination = stem.with_suffix(suffix)
        fig.savefig(destination, dpi=200, bbox_inches="tight", facecolor="white")
        set_read_only(destination)
    plt.close(fig)


def status_figure(stem: Path, title: str, detail: str) -> None:
    fig, ax = plt.subplots(figsize=(10.8, 5.6))
    ax.axis("off")
    ax.text(0.5, 0.72, title, ha="center", va="center", fontsize=24, weight="bold")
    ax.text(
        0.5,
        0.45,
        "FOUNDATION RUNNING",
        ha="center",
        va="center",
        fontsize=22,
        color="#a33a2b",
        weight="bold",
    )
    ax.text(0.5, 0.25, detail, ha="center", va="center", fontsize=13, color="#333333")
    ax.text(
        0.5,
        0.07,
        "No incomplete-population value is present in this figure.",
        ha="center",
        va="center",
        fontsize=11,
        color="#6a6a6a",
    )
    save_figure(fig, stem)


def weighted_add(target: np.ndarray | None, value: Any, weight: float) -> np.ndarray:
    array_value = np.asarray(value, dtype=float) * weight
    return array_value.copy() if target is None else target + array_value


def merge_native_response(samples: Mapping[str, Mapping[str, Any]], threshold: int) -> dict[str, Any]:
    tag = f"pt{threshold}"
    required = ("auau_photon12", "auau_photon20")
    first = samples[required[0]]["reductions"][tag]
    merged: dict[str, Any] = {
        "reco_ptgamma_edges": list(map(float, first["reco_ptgamma_edges"])),
        "truth_ptgamma_edges": list(map(float, first["truth_ptgamma_edges"])),
        "xj_edges": list(map(float, first["xj_edges"])),
        "photon": {},
        "xj": {},
        "signal_leakage": {"truth_matched_abcd": {}, "truth_matched_abcd_sumw2": {}},
    }
    photon_sumw2 = {"reco_sumw2", "truth_sumw2", "fakes_sumw2", "misses_sumw2", "response_sumw2_truth_x_reco"}
    xj_sumw2 = {
        "reco_sumw2", "truth_sumw2", "fakes_sumw2", "misses_sumw2",
        "detector_fakes_reco_sumw2", "fake_photon_reco_sumw2", "combinatoric_reco_sumw2",
        "response_sumw2_truth_global_x_reco_global",
    }
    for sample_id in required:
        item = samples[sample_id]
        contract = item.get("sample_contract")
        if (
            not isinstance(contract, Mapping)
            or contract.get("system") != "auau"
            or contract.get("source_class") != "photon_signal"
        ):
            raise ValueError(
                f"{sample_id} is not an Au+Au photon-signal response source; "
                "embedded inclusive weights are forbidden in this merge"
            )
        reduction = item["reductions"][tag]
        if (
            reduction["reco_ptgamma_edges"] != first["reco_ptgamma_edges"]
            or reduction["truth_ptgamma_edges"] != first["truth_ptgamma_edges"]
            or reduction["xj_edges"] != first["xj_edges"]
        ):
            raise ValueError(f"native response axes differ: {sample_id}.{tag}")
        weight = float(item["cross_section_weight_pb_per_event"])
        for key, value in reduction["photon"].items():
            factor = weight * weight if key in photon_sumw2 else weight
            merged["photon"][key] = weighted_add(merged["photon"].get(key), value, factor)
        for key, value in reduction["xj"].items():
            factor = weight * weight if key in xj_sumw2 else weight
            merged["xj"][key] = weighted_add(merged["xj"].get(key), value, factor)
        leakage = reduction["signal_leakage"]
        for category in "ABCD":
            merged["signal_leakage"]["truth_matched_abcd"][category] = weighted_add(
                merged["signal_leakage"]["truth_matched_abcd"].get(category),
                leakage["truth_matched_abcd"][category],
                weight,
            )
            merged["signal_leakage"]["truth_matched_abcd_sumw2"][category] = weighted_add(
                merged["signal_leakage"]["truth_matched_abcd_sumw2"].get(category),
                leakage["truth_matched_abcd_sumw2"][category],
                weight * weight,
            )
    return merged


def th1(directory: Any, name: str, edges: Sequence[float], values: Any, variances: Any) -> None:
    directory.cd()
    histogram = ROOT.TH1D(name, name, len(edges) - 1, array("d", list(map(float, edges))))
    for index, (value, variance) in enumerate(zip(np.asarray(values), np.asarray(variances)), start=1):
        histogram.SetBinContent(index, float(value))
        histogram.SetBinError(index, math.sqrt(max(float(variance), 0.0)))
    histogram.Write()


def th2(
    directory: Any,
    name: str,
    xedges: Sequence[float],
    yedges: Sequence[float],
    values: Any,
    variances: Any,
) -> None:
    directory.cd()
    histogram = ROOT.TH2D(
        name,
        name,
        len(xedges) - 1,
        array("d", list(map(float, xedges))),
        len(yedges) - 1,
        array("d", list(map(float, yedges))),
    )
    matrix = np.asarray(values, dtype=float)
    variance_matrix = np.asarray(variances, dtype=float)
    if matrix.shape != (len(xedges) - 1, len(yedges) - 1):
        raise ValueError(f"TH2 shape differs for {name}: {matrix.shape}")
    for xindex in range(matrix.shape[0]):
        for yindex in range(matrix.shape[1]):
            histogram.SetBinContent(xindex + 1, yindex + 1, float(matrix[xindex, yindex]))
            histogram.SetBinError(
                xindex + 1,
                yindex + 1,
                math.sqrt(max(float(variance_matrix[xindex, yindex]), 0.0)),
            )
    histogram.Write()


def create_photon_response_root(path: Path, native: Mapping[int, Mapping[str, Any]], provenance: Mapping[str, Any]) -> None:
    ROOT.gROOT.SetBatch(True)
    root_file = ROOT.TFile.Open(str(path), "CREATE")
    if not root_file or root_file.IsZombie():
        raise ValueError(f"cannot create ROOT output: {path}")
    try:
        ROOT.TNamed("THE255_STATUS", "ACCEPTED_NOMINAL_H70_RESPONSE_INPUT__NOT_UNFOLDED").Write()
        ROOT.TObjString(json.dumps(provenance, sort_keys=True)).Write("THE255_PROVENANCE_JSON")
        for threshold, payload in native.items():
            directory = root_file.mkdir(f"pt{threshold}")
            photon = payload["photon"]
            truth_edges = payload["truth_ptgamma_edges"]
            reco_edges = payload["reco_ptgamma_edges"]
            th2(directory, "response_truth_x_reco", truth_edges, reco_edges,
                photon["response_truth_x_reco"], photon["response_sumw2_truth_x_reco"])
            th1(directory, "truth", truth_edges, photon["truth"], photon["truth_sumw2"])
            th1(directory, "misses", truth_edges, photon["misses"], photon["misses_sumw2"])
            th1(directory, "reco", reco_edges, photon["reco"], photon["reco_sumw2"])
            th1(directory, "fakes", reco_edges, photon["fakes"], photon["fakes_sumw2"])
            root_file.cd()
    finally:
        root_file.Close()
    set_read_only(path)


def create_xj_response_root(path: Path, native: Mapping[int, Mapping[str, Any]], provenance: Mapping[str, Any]) -> None:
    ROOT.gROOT.SetBatch(True)
    root_file = ROOT.TFile.Open(str(path), "CREATE")
    if not root_file or root_file.IsZombie():
        raise ValueError(f"cannot create ROOT output: {path}")
    try:
        ROOT.TNamed("THE255_STATUS", "ACCEPTED_NOMINAL_H70_RESPONSE_INPUT__NOT_UNFOLDED").Write()
        ROOT.TObjString(json.dumps(provenance, sort_keys=True)).Write("THE255_PROVENANCE_JSON")
        for threshold, payload in native.items():
            directory = root_file.mkdir(f"pt{threshold}")
            xj = payload["xj"]
            reco_pt = payload["reco_ptgamma_edges"]
            truth_pt = payload["truth_ptgamma_edges"]
            xj_edges = payload["xj_edges"]
            truth_global_bins = (len(truth_pt) - 1) * (len(xj_edges) - 1)
            reco_global_bins = (len(reco_pt) - 1) * (len(xj_edges) - 1)
            global_truth_edges = np.arange(truth_global_bins + 1, dtype=float)
            global_reco_edges = np.arange(reco_global_bins + 1, dtype=float)
            th2(
                directory,
                "response_truth_global_x_reco_global",
                global_truth_edges,
                global_reco_edges,
                xj["response_truth_global_x_reco_global"],
                xj["response_sumw2_truth_global_x_reco_global"],
            )
            for name, values, variances, pt_edges in (
                ("truth_ptgamma_x_xj", xj["truth"], xj["truth_sumw2"], truth_pt),
                ("misses_ptgamma_x_xj", xj["misses"], xj["misses_sumw2"], truth_pt),
                ("reco_ptgamma_x_xj", xj["reco"], xj["reco_sumw2"], reco_pt),
                ("fakes_ptgamma_x_xj", xj["fakes"], xj["fakes_sumw2"], reco_pt),
                ("detector_fakes_ptgamma_x_xj", xj["detector_fakes_reco"], xj["detector_fakes_reco_sumw2"], reco_pt),
                ("fake_photon_ptgamma_x_xj", xj["fake_photon_reco"], xj["fake_photon_reco_sumw2"], reco_pt),
                ("combinatoric_ptgamma_x_xj", xj["combinatoric_reco"], xj["combinatoric_reco_sumw2"], reco_pt),
            ):
                th2(directory, name, pt_edges, xj_edges, values, variances)
            ROOT.TObjString(json.dumps({
                "flattening": "global_bin = ptgamma_bin * n_xj + xj_bin",
                "truth_ptgamma_edges": truth_pt,
                "reco_ptgamma_edges": reco_pt,
                "xj_edges": xj_edges,
                "overflow_support": "2.0-3.0 retained",
            }, sort_keys=True)).Write("GLOBAL_BIN_MAP_JSON")
            root_file.cd()
    finally:
        root_file.Close()
    set_read_only(path)


def create_data_status_root(path: Path, observation: Mapping[str, Any]) -> None:
    root_file = ROOT.TFile.Open(str(path), "CREATE")
    if not root_file or root_file.IsZombie():
        raise ValueError(f"cannot create ROOT status output: {path}")
    try:
        ROOT.TNamed("THE255_STATUS", STATUS).Write()
        ROOT.TNamed("THE255_DATA_CONTENT", "NONE__FINAL_AUAU_MEMBERSHIP_UNSEALED").Write()
        ROOT.TObjString(json.dumps(observation, sort_keys=True)).Write("THE255_FOUNDATION_OBSERVATION_JSON")
    finally:
        root_file.Close()
    set_read_only(path)


def response_scan(bundle: Any) -> list[dict[str, Any]]:
    xj_truth = bundle.xj_truth.reshape(-1)
    xj_fakes = (bundle.xj_detector_fakes_reco + bundle.xj_boundary_fakes_reco).reshape(-1)
    xj_reco = bundle.xj_response.sum(axis=0) + xj_fakes
    photon_truth = bundle.photon_truth
    photon_fakes = bundle.photon_boundary_fakes_reco
    photon_reco = bundle.photon_response.sum(axis=0) + photon_fakes
    rows: list[dict[str, Any]] = []
    prior_xj: np.ndarray | None = None
    for iterations in range(2, 13):
        unfolded_xj, refolded_xj = iterative_bayes(
            xj_reco,
            bundle.xj_response,
            bundle.xj_misses,
            iterations,
            xj_truth,
            fakes=xj_fakes,
        )
        unfolded_photon, refolded_photon = iterative_bayes(
            photon_reco,
            bundle.photon_response,
            bundle.photon_misses,
            iterations,
            photon_truth,
            fakes=photon_fakes,
        )
        relative_change = (
            float(np.sum(np.abs(unfolded_xj - prior_xj)) / max(np.sum(np.abs(prior_xj)), 1.0e-30))
            if prior_xj is not None else None
        )
        rows.append({
            "iterations": iterations,
            "scope": "RESPONSE_ONLY_MC_TRUTH_CLOSURE__NOT_DATA_CANDIDATE_SELECTION",
            "xj_truth_closure_chi2_ndf": chi2_ndf(
                xj_truth, unfolded_xj, np.clip(xj_truth + unfolded_xj, 1.0e-30, None)
            ),
            "xj_refold_chi2_ndf": chi2_ndf(xj_reco, refolded_xj, np.clip(xj_reco, 1.0e-30, None)),
            "photon_truth_closure_chi2_ndf": chi2_ndf(
                photon_truth,
                unfolded_photon,
                np.clip(photon_truth + unfolded_photon, 1.0e-30, None),
            ),
            "photon_refold_chi2_ndf": chi2_ndf(
                photon_reco, refolded_photon, np.clip(photon_reco, 1.0e-30, None)
            ),
            "xj_relative_l1_change_from_previous": relative_change,
        })
        prior_xj = unfolded_xj
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=REPO)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--foundation-observation", type=Path, required=True)
    parser.add_argument("--recovery-snapshot", type=Path, required=True)
    parser.add_argument("--accounting-correction", type=Path, required=True)
    parser.add_argument("--event-gate-receipt", type=Path, required=True)
    parser.add_argument("--response-catalog", type=Path, required=True)
    parser.add_argument("--response-receipt", type=Path, required=True)
    args = parser.parse_args()
    repo = args.repo.resolve()
    output = args.output.resolve()
    if output.exists():
        raise ValueError(f"immutable output already exists: {output}")
    output.mkdir(parents=True)

    paths = {
        "foundation_observation": args.foundation_observation.resolve(),
        "recovery_snapshot": args.recovery_snapshot.resolve(),
        "accounting_correction": args.accounting_correction.resolve(),
        "event_gate_receipt": args.event_gate_receipt.resolve(),
        "response_catalog": args.response_catalog.resolve(),
        "response_receipt": args.response_receipt.resolve(),
        "reducer": repo / "scripts/data_prep/recoiljets/reduce_the243_schema10_xjgamma.py",
        "corrections": repo / "scripts/data_prep/recoiljets/the243_h70_corrections.py",
        "analysis_config": repo / "macros/analysis_config.yaml",
        "model": repo / "InputFiles/the134_auau_h70_model.tmva.root",
        "prior_k11": repo / "dataOutput/the243_golden_ppg_analysis_closure_deck_20260821/ttree_h70_z10_pt5_responsek_integratedpurity_s4/THE243_H70_COMMON1535_SELECTED_RESULT.json",
        "prior_k11_scan": repo / "dataOutput/the243_golden_ppg_analysis_closure_deck_20260821/ttree_h70_z10_pt5_responsek_integratedpurity_s4/THE243_H70_COMMON1535_OBJECTIVE_SCAN.json",
        "historical_points": repo / "dataOutput/the219_friday_ppg_20260814/the89_publication_xjgamma_overlay/the89_unfolded_xjgamma_pp_auau020_xjge0p5_points.json",
        "historical_covariance": repo / "dataOutput/the85_auau_xjgamma_unfolding_push/unfolded_xjgamma_firstpass/the85_unfolded_xjgamma_auau_0_20_responseKfix_iter5_covariance.npz",
    }
    for label, path in paths.items():
        if label == "model" and not path.exists():
            continue
        if not path.exists():
            raise ValueError(f"required source absent ({label}): {path}")

    observation = load(paths["foundation_observation"])
    if observation.get("status") != STATUS or observation["base_membership"].get("sealed") is not False:
        raise ValueError("foundation observation does not authorize fail-closed package")
    catalog = load(paths["response_catalog"])
    receipt = load(paths["response_receipt"])
    if sha256(paths["response_catalog"]) != RESPONSE_CATALOG_SHA256:
        raise ValueError("accepted response catalog SHA differs")
    if (
        receipt.get("status") != "PASS"
        or receipt.get("result_stage") != "RESPONSE_INPUT__NOT_UNFOLDED"
        or receipt.get("input_root_count") != 7145
        or receipt.get("dst_reads") != 0
    ):
        raise ValueError("accepted response terminal contract differs")

    samples: dict[str, dict[str, Any]] = {}
    for sample_id in ("auau_photon12", "auau_photon20"):
        binding = receipt["samples"][sample_id]
        path = Path(binding["output_path"])
        if not path.exists() or sha256(path) != binding["output_sha256"]:
            raise ValueError(f"accepted response sample binding differs: {sample_id}")
        sample = load(path)
        if sample.get("result_status") != "RESPONSE_INPUT__NOT_UNFOLDED":
            raise ValueError(f"response result stage differs: {sample_id}")
        sample_contract = sample.get("sample_contract")
        if (
            not isinstance(sample_contract, Mapping)
            or sample_contract.get("system") != "auau"
            or sample_contract.get("source_class") != "photon_signal"
        ):
            raise ValueError(
                f"{sample_id} is not an Au+Au photon-signal response source; "
                "embedded inclusive sigma_eff/Npass samples cannot enter this package"
            )
        samples[sample_id] = sample

    native = {threshold: merge_native_response(samples, threshold) for threshold in THRESHOLDS}
    bundles = {
        threshold: build_response_bundle(samples, "auau", threshold, native[threshold]["xj_edges"])
        for threshold in THRESHOLDS
    }
    fine_bundle = build_response_bundle(samples, "auau", 5, FINE_XJ_EDGES)
    iteration_scan = response_scan(fine_bundle)

    provenance = {
        "schema": "THE255AuAuProvenanceV1",
        "status": STATUS,
        "created_at": utc_now(),
        "data_membership_catalog": None,
        "data_membership_catalog_status": "FINAL_AUAU_BASE_CATALOG_NOT_YET_SEALED",
        "event_gate_uniform_coverage": False,
        "model_sha256": MODEL_SHA256,
        "analysis_config": {"path": str(paths["analysis_config"]), "sha256": sha256(paths["analysis_config"])},
        "reducer": {"path": str(paths["reducer"]), "sha256": sha256(paths["reducer"])},
        "corrections": {"path": str(paths["corrections"]), "sha256": sha256(paths["corrections"])},
        "response_catalog": {"path": str(paths["response_catalog"]), "sha256": RESPONSE_CATALOG_SHA256},
        "response_receipt": {"path": str(paths["response_receipt"]), "sha256": sha256(paths["response_receipt"])},
        "response_weights_pb_per_event": {
            sample_id: float(samples[sample_id]["cross_section_weight_pb_per_event"])
            for sample_id in samples
        },
        "embedded_inclusive_response_samples_used": 0,
        "embedded_inclusive_weight_contract": (
            "not applicable to this photon-signal-only response package; any future "
            "embedded-inclusive consumer must require the complete producer -> "
            "sigma_eff/Npass -> centrality artifact"
        ),
        "unfolding_implementation": "scripts/data_prep/recoiljets/the243_h70_corrections.py::iterative_bayes",
        "iteration_selection": "PENDING_FINAL_DATA; response-only k=2..12 closure scan included",
        "normalization": "PENDING_FINAL_DATA; sum unfolded jets / (sum unfolded photons * xJ bin width)",
        "data_ttree_reads": 0,
        "dst_reads": 0,
        "pp_inputs_used": 0,
        "historical_h0_inputs_used_as_nominal": 0,
    }
    write_json(output / "01_PROVENANCE.json", provenance)
    write_json(output / "02_DATA_COVERAGE.json", {
        "schema": "THE255AuAuDataCoverageV1",
        "status": STATUS,
        "final_auau_base_catalog_ready": False,
        "final_auau_base_root_count": None,
        "provisional_auau_base_root_count": observation["base_membership"]["current_provisional_roots"],
        "target_auau_base_root_count": observation["base_membership"]["target_if_all_recovery_rows_close"],
        "event_gate_uniform_coverage": False,
        "event_gate_covered_roots": observation["event_gate"]["existing_validated_companion_roots"],
        "event_gate_missing_roots": None,
        "event_gate_missing_roots_status": observation["event_gate"]["missing_root_count_status"],
        "base_catalog_path": None,
        "event_gate_catalog_or_receipt": str(paths["event_gate_receipt"]),
        "foundation_read_only_observation": observation,
        "data_reduction_performed": False,
    })

    nominal_contract = f"""# Nominal AuAu H70 contract

Status: **FROZEN; DATA EXECUTION WAITING ON FINAL MEMBERSHIP**

- Event: direct GL1 `ScaledVector` bit 22; `Photon 10 GeV + MBD N&S >= 2, trigger |z| < 150 cm`; `MinimumBiasInfo::isAuAuMinimumBias()==true`; offline `|z|<10 cm`.
- Centrality: `CentralityInfo::PROP::mbd_NS`; nominal 0-20%.
- Photon: 15-35 GeV, |eta|<0.7, event-leading category semantics.
- H70 model SHA-256: `{MODEL_SHA256}`; 14 features; 70 MeV affects shower-feature construction only.
- Tight: `score > 0.6927385742113304 + 0.0010476808935335573 * centrality`.
- NonTight: `WP90 < score < WP80`, using the frozen centrality-only surfaces.
- Isolation: UE-subtracted SUB1 CEMC+HCALIN+HCALOUT tower ET in R=0.4 minus photon ET; no core hole and no 70 MeV isolation-tower threshold. Thresholds: 7.241, 6.583, 5.925, 5.267, 4.609, 3.951, 2.964 GeV.
- Jets: anti-kT R=0.4 SUB1; corrected pT >=5 GeV response basis with 5/7/10/12 GeV fanout; |eta|<0.7; photon-jet DeltaR>=0.4; recoil DeltaPhi>=7pi/8; all selected recoil jets.

Authority: `scripts/data_prep/recoiljets/reduce_the243_schema10_xjgamma.py`, `macros/analysis_config.yaml`, and the user-frozen THE-255 contract. Historical H0 direct histograms are excluded from nominal input.
"""
    write_text(output / "03_NOMINAL_AUAU_CONTRACT.md", nominal_contract)

    pt_bins = [(15.0, 20.0), (20.0, 25.0), (25.0, 35.0)]
    status_note = "WITHHELD__FINAL_AUAU_BASE_CATALOG_AND_UNIFORM_EVENT_GATE_COVERAGE_OPEN"
    write_csv(output / "04_ABCD_COUNTS.csv",
              ("status", "ptgamma_low_gev", "ptgamma_high_gev", "region", "event_leading_count", "stat_error"),
              ({"status": status_note, "ptgamma_low_gev": lo, "ptgamma_high_gev": hi,
                "region": region, "event_leading_count": "", "stat_error": ""}
               for lo, hi in pt_bins for region in "ABCD"))

    leakage_rows = []
    for pindex, (lo, hi) in enumerate(pt_bins):
        a = fine_bundle.leakage[0, pindex]
        for cindex, category in enumerate("ABCD"):
            leakage_rows.append({
                "status": "RESPONSE_LEAKAGE_READY__DATA_PURITY_PENDING",
                "ptgamma_low_gev": lo,
                "ptgamma_high_gev": hi,
                "category": category,
                "weighted_prompt_leakage_yield": fine_bundle.leakage[cindex, pindex],
                "leakage_relative_to_A": fine_bundle.leakage[cindex, pindex] / a if a > 0 else "",
                "data_purity": "",
                "data_purity_stat_error": "",
                "alpha": "",
            })
    write_csv(output / "05_PURITY.csv",
              ("status", "ptgamma_low_gev", "ptgamma_high_gev", "category",
               "weighted_prompt_leakage_yield", "leakage_relative_to_A", "data_purity",
               "data_purity_stat_error", "alpha"), leakage_rows)
    recoil_fields = ("status", "ptgamma_low_gev", "ptgamma_high_gev", "xj_low", "xj_high", "yield", "stat_error")
    blank_recoil = [{"status": status_note, "ptgamma_low_gev": "", "ptgamma_high_gev": "",
                     "xj_low": "", "xj_high": "", "yield": "", "stat_error": ""}]
    write_csv(output / "06_A_RECOIL.csv", recoil_fields, blank_recoil)
    write_csv(output / "07_C_RECOIL.csv", recoil_fields, blank_recoil)
    write_csv(output / "08_AC_CORRECTED_RECOIL.csv", recoil_fields, blank_recoil)

    write_text(output / "09_K_COMBINATORIC_CONTRACT.md", """# Nominal AuAu K/combinatoric contract

Status: **SIM TEMPLATE READY; ABSOLUTE DATA SUBTRACTION PENDING FINAL MEMBERSHIP**

For each reconstructed photon-pT and xJ bin, the exact implementation defines

```text
K = Tcomb / Ngamma_SIM
absolute subtraction = K * Ngamma_DATA
```

`Tcomb` is the weighted inclusive unmatched-recoil template for truth-matched event-leading Region-A photons. A reconstructed recoil is combinatoric when no anti-kT R=0.4 truth jet above the active jet threshold lies within DeltaR<0.4. `Ngamma_SIM` is the weighted matched Region-A reconstructed-photon normalization in the same reconstructed pTgamma bin.

The absolute subtraction is applied after leakage-aware Region-C correction and before joint unfolding. `10_K_TEMPLATE.csv` is the accepted response-derived per-photon template. `11_K_SUBTRACTION.csv` and `12_K_CORRECTED_RECOIL.csv` deliberately contain no incomplete-population values.

Authority: `scripts/data_prep/recoiljets/the243_h70_corrections.py::correct_measured` and `build_response_bundle`.
""")
    k_rows = []
    xj_edges = fine_bundle.xj_edges
    for pindex, (lo, hi) in enumerate(pt_bins):
        denominator = fine_bundle.combinatoric_normalization_reco[pindex]
        denominator_var = fine_bundle.combinatoric_normalization_reco_sumw2[pindex]
        for xindex, (xlo, xhi) in enumerate(zip(xj_edges[:-1], xj_edges[1:])):
            numerator = fine_bundle.combinatoric_reco[pindex, xindex]
            numerator_var = fine_bundle.combinatoric_reco_sumw2[pindex, xindex]
            value = numerator / denominator if denominator > 0 else 0.0
            variance = (
                numerator_var / denominator**2 + numerator**2 * denominator_var / denominator**4
                if denominator > 0 else 0.0
            )
            k_rows.append({
                "status": "ACCEPTED_H70_RESPONSE_TEMPLATE",
                "jet_pt_threshold_gev": 5,
                "ptgamma_low_gev": lo,
                "ptgamma_high_gev": hi,
                "xj_low": xlo,
                "xj_high": xhi,
                "tcomb_weighted": numerator,
                "ngamma_sim_weighted": denominator,
                "k_per_photon": value,
                "k_stat_error": math.sqrt(max(variance, 0.0)),
            })
    write_csv(output / "10_K_TEMPLATE.csv",
              ("status", "jet_pt_threshold_gev", "ptgamma_low_gev", "ptgamma_high_gev",
               "xj_low", "xj_high", "tcomb_weighted", "ngamma_sim_weighted",
               "k_per_photon", "k_stat_error"), k_rows)
    k_pending_fields = ("status", "ptgamma_low_gev", "ptgamma_high_gev", "xj_low", "xj_high",
                        "ngamma_data", "absolute_k_subtraction", "stat_error")
    k_pending = [{"status": status_note, "ptgamma_low_gev": "", "ptgamma_high_gev": "",
                  "xj_low": "", "xj_high": "", "ngamma_data": "",
                  "absolute_k_subtraction": "", "stat_error": ""}]
    write_csv(output / "11_K_SUBTRACTION.csv", k_pending_fields, k_pending)
    write_csv(output / "12_K_CORRECTED_RECOIL.csv", recoil_fields, blank_recoil)

    write_text(output / "13_H70_RESPONSE_WEIGHTING_CONTRACT.md", f"""# Accepted nominal H70 response weighting contract

Status: **PASS — RESPONSE INPUT, NOT UNFOLDED**

- Catalog: `{paths['response_catalog']}`
- Catalog SHA-256: `{RESPONSE_CATALOG_SHA256}`
- Accepted inventory: 7,145 unique schema-10 ROOT inputs; `dst_reads=0`.
- Production campaign: `the121_the122_sim_prod_schema10_20260818eof_t1_photon_829eec8a`.
- Production-exclusive generator ownership: Photon12 owns generator photon [12,20) GeV and Photon20 owns >=20 GeV. All owned events are then summed so reco/truth migrations across 20 GeV remain represented once.
- Photon12: `2598.12425 / 9991946 = {samples['auau_photon12']['cross_section_weight_pb_per_event']:.17g}` pb/event.
- Photon20: `133.317866 / 9991946 = {samples['auau_photon20']['cross_section_weight_pb_per_event']:.17g}` pb/event.
- Histogram contents receive `w`; Sumw2 receives `w^2`; no later visual or physics rescaling is applied.
- Nominal model SHA-256: `{MODEL_SHA256}`.
- Response pT edges: reco [10,15,20,25,35,40] GeV; truth [5,10,15,20,25,35,40] GeV.
- Native xJ support: 0.05 bins through 2.0 plus a retained 2.0-3.0 overflow-support bin.
- Jet-threshold fanout: 5, 7, 10, 12 GeV; nominal response basis begins at corrected jet pT >=5 GeV.

There is no p+p SI/DI or 0/1.5 mrad machinery in this AuAu response. Historical H0 response products are not used.
""")

    create_data_status_root(output / "DATA_REDUCTION_STATUS.root", observation)
    create_photon_response_root(output / "14_PHOTON_RESPONSE.root", native, provenance)
    create_xj_response_root(output / "15_XJ_RESPONSE.root", native, provenance)

    bundle_audits = {}
    for threshold, bundle in bundles.items():
        xj_total = bundle.xj_response.sum(axis=1) + bundle.xj_misses
        photon_total = bundle.photon_response.sum(axis=1) + bundle.photon_misses
        bundle_audits[str(threshold)] = {
            "photon_response_sum": float(bundle.photon_response.sum()),
            "photon_miss_sum": float(bundle.photon_misses.sum()),
            "photon_boundary_fake_sum": float(bundle.photon_boundary_fakes_reco.sum()),
            "photon_min_supported_efficiency": float(np.min(np.divide(
                bundle.photon_response.sum(axis=1), photon_total,
                out=np.zeros_like(photon_total), where=photon_total > 0)[photon_total > 0])),
            "xj_response_sum": float(bundle.xj_response.sum()),
            "xj_miss_sum": float(bundle.xj_misses.sum()),
            "xj_detector_fake_sum": float(bundle.xj_detector_fakes_reco.sum()),
            "xj_boundary_fake_sum": float(bundle.xj_boundary_fakes_reco.sum()),
            "xj_min_supported_efficiency": float(np.min(np.divide(
                bundle.xj_response.sum(axis=1), xj_total,
                out=np.zeros_like(xj_total), where=xj_total > 0)[xj_total > 0])),
            "xj_reco_partition_max_abs_residual": float(np.max(np.abs(
                bundle.xj_reco - (
                    bundle.xj_response.sum(axis=0).reshape(bundle.xj_reco.shape)
                    + bundle.xj_fakes_reco + bundle.xj_boundary_fakes_reco
                )))),
            "xj_truth_partition_max_abs_residual": float(np.max(np.abs(
                bundle.xj_truth - (
                    bundle.xj_response.sum(axis=1).reshape(bundle.xj_truth.shape)
                    + bundle.xj_misses.reshape(bundle.xj_truth.shape)
                )))),
        }
    prior = load(paths["prior_k11"])
    write_json(output / "16_UNFOLDING_DIAGNOSTICS.json", {
        "schema": "THE255AuAuUnfoldingDiagnosticsV1",
        "status": "RESPONSE_DIAGNOSTICS_READY__DATA_UNFOLDING_PENDING_FINAL_MEMBERSHIP",
        "response_catalog_sha256": RESPONSE_CATALOG_SHA256,
        "response_input_status": "RESPONSE_INPUT__NOT_UNFOLDED",
        "response_bundle_audits": bundle_audits,
        "response_only_mc_closure_iteration_scan": iteration_scan,
        "data_iteration_scan_required": list(range(2, 13)),
        "data_iteration_selection_status": "NOT_RUN__FINAL_DATA_MEMBERSHIP_UNSEALED",
        "selection_rule": {
            "implementation": "scripts/slides/the243_h70_money_plot/make_the243_h70_final_overlay.py::candidate_score",
            "reference_shape_used_in_score": False,
            "requires": ["final data refolding", "ABCD physicality", "statistical precision",
                         "negative-input pressure", "MC truth closure", "response observability"],
            "may_not_be_evaluated_without_final_data": True,
        },
        "statistical_toys_required_for_selected_final_candidate": 6000,
        "statistical_toys_run_for_current_full_membership": 0,
        "prior_k11_checkpoint": {
            "label": "PRELIMINARY_STATISTICAL_CANDIDATE__DIRECT_PARITY_OPEN__INCOMPLETE_OLDER_DATA_POPULATION",
            "iterations": prior["selection"]["iterations"],
            "final_toys": prior["selection"]["final_toys"],
            "auau_refold_chi2_ndf": prior["auau"]["refold_chi2_ndf"],
            "auau_photon_refold_chi2_ndf": prior["auau"]["photon_refold_chi2_ndf"],
            "source": str(paths["prior_k11"]),
            "source_sha256": sha256(paths["prior_k11"]),
            "promoted": False,
        },
        "unfolded_data_present": False,
        "final_candidate_selected": False,
    })
    write_csv(output / "17_UNFOLDED_PHOTON_DENOMINATOR.csv",
              ("status", "ptgamma_low_gev", "ptgamma_high_gev", "unfolded_ngamma", "stat_error", "normalization_role"),
              [{"status": status_note, "ptgamma_low_gev": "", "ptgamma_high_gev": "",
                "unfolded_ngamma": "", "stat_error": "",
                "normalization_role": "sum over pTgamma after photon unfolding"}])
    write_csv(output / "18_FINAL_AUAU_POINTS.csv",
              ("status", "xj_low", "xj_high", "xj_center", "value", "stat_error", "covariance_reference", "normalization"),
              [{"status": status_note, "xj_low": "", "xj_high": "", "xj_center": "",
                "value": "", "stat_error": "", "covariance_reference": "",
                "normalization": "sum unfolded jets / (sum unfolded photons * delta_xj)"}])

    historical = load(paths["historical_points"])
    historical_rows = [{
        "contract_label": "THE-89 AuAu response/K statistical-correction candidate; diagnostic/noncanonical",
        "displayed": row["displayed"],
        "xj_low": row["xj_low"],
        "xj_high": row["xj_high"],
        "xj_center": row["xj_center"],
        "historical_value": row["auau_0_20"]["value"],
        "historical_stat_error": row["auau_0_20"]["stat_error"],
        "current_the255_value": "",
        "current_the255_stat_error": "",
        "current_status": STATUS,
    } for row in historical["bins"]]
    write_csv(output / "19_HISTORICAL_AUAU_COMPARISON.csv",
              ("contract_label", "displayed", "xj_low", "xj_high", "xj_center",
               "historical_value", "historical_stat_error", "current_the255_value",
               "current_the255_stat_error", "current_status"), historical_rows)

    plt.rcParams.update({
        "font.size": 11.5,
        "axes.titlesize": 15,
        "axes.labelsize": 12.5,
        "legend.fontsize": 9.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })
    status_figure(output / "ABCD_PURITY_STATUS", "H70 ABCD / purity",
                  "Final AuAu base catalog and uniform event-gate coverage are not sealed.")
    status_figure(output / "DETECTOR_LEVEL_RECOIL_CORRECTIONS_STATUS",
                  "A, C, ABCD and K detector-level recoil",
                  "K per-photon SIM template is ready; every DATA-dependent rung is withheld.")

    native5 = native[5]
    photon_matrix = np.asarray(native5["photon"]["response_truth_x_reco"], dtype=float)
    row_sum = photon_matrix.sum(axis=1, keepdims=True)
    photon_conditional = np.divide(photon_matrix, row_sum, out=np.zeros_like(photon_matrix), where=row_sum > 0)
    fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.4), gridspec_kw={"width_ratios": [1.35, 1.0]})
    image = axes[0].imshow(photon_conditional, origin="lower", aspect="auto", cmap="viridis", vmin=0)
    axes[0].set_xticks(range(photon_matrix.shape[1]))
    axes[0].set_xticklabels([f"{a:g}-{b:g}" for a, b in zip(native5["reco_ptgamma_edges"][:-1], native5["reco_ptgamma_edges"][1:])], rotation=35, ha="right")
    axes[0].set_yticks(range(photon_matrix.shape[0]))
    axes[0].set_yticklabels([f"{a:g}-{b:g}" for a, b in zip(native5["truth_ptgamma_edges"][:-1], native5["truth_ptgamma_edges"][1:])])
    axes[0].set_xlabel("reco pTgamma bin (GeV)")
    axes[0].set_ylabel("truth pTgamma bin (GeV)")
    axes[0].set_title("Photon response P(reco | matched truth)")
    fig.colorbar(image, ax=axes[0], label="row-normalized probability")
    bundle = bundles[5]
    totals = bundle.photon_response.sum(axis=1) + bundle.photon_misses
    efficiency = np.divide(bundle.photon_response.sum(axis=1), totals, out=np.zeros_like(totals), where=totals > 0)
    centers = 0.5 * (PT_EDGES[:-1] + PT_EDGES[1:])
    axes[1].plot(centers, efficiency, "o-", color="#184f78")
    axes[1].set_ylim(0, 1.05)
    axes[1].set_xlabel("truth pTgamma (GeV)")
    axes[1].set_ylabel("matched reconstruction efficiency")
    axes[1].set_title("15-35 GeV analysis support")
    axes[1].grid(alpha=0.25)
    fig.suptitle("Accepted nominal AuAu H70 photon response — jet pT >=5 GeV")
    save_figure(fig, output / "PHOTON_RESPONSE")

    matrix = bundle.xj_response
    conditional = np.divide(matrix, matrix.sum(axis=1, keepdims=True), out=np.zeros_like(matrix), where=matrix.sum(axis=1, keepdims=True) > 0)
    fig, ax = plt.subplots(figsize=(10.5, 8.0))
    image = ax.imshow(conditional, origin="lower", aspect="auto", cmap="magma", vmin=0,
                      vmax=float(np.quantile(conditional[conditional > 0], 0.995)))
    ax.set_xlabel("reco global bin: pTgamma x xJ")
    ax.set_ylabel("truth global bin: pTgamma x xJ")
    ax.set_title("Accepted nominal AuAu H70 joint response — native 0.05 xJ support")
    fig.colorbar(image, ax=ax, label="row-normalized probability")
    fig.text(0.01, 0.01, "2.0-3.0 overflow-support bin retained; ptjet >=5 GeV", fontsize=9, color="#555555")
    save_figure(fig, output / "XJ_RESPONSE")

    xj_total = bundle.xj_response.sum(axis=1) + bundle.xj_misses
    xj_eff = np.divide(bundle.xj_response.sum(axis=1), xj_total, out=np.zeros_like(xj_total), where=xj_total > 0).reshape(bundle.xj_truth.shape)
    xcenters = 0.5 * (bundle.xj_edges[:-1] + bundle.xj_edges[1:])
    fig, axes = plt.subplots(1, 2, figsize=(13.2, 5.3))
    for pindex, (lo, hi) in enumerate(pt_bins):
        axes[0].step(xcenters, xj_eff[pindex], where="mid", label=f"{lo:g}-{hi:g} GeV")
    axes[0].set_xlabel("truth xJgamma")
    axes[0].set_ylabel("matched response efficiency")
    axes[0].set_ylim(0, 1.05)
    axes[0].set_title("Efficiency / misses")
    axes[0].grid(alpha=0.25)
    axes[0].legend()
    matched = bundle.xj_response.sum(axis=0).reshape(bundle.xj_reco.shape).sum(axis=0)
    detector_fake = bundle.xj_detector_fakes_reco.sum(axis=0)
    boundary_fake = bundle.xj_boundary_fakes_reco.sum(axis=0)
    axes[1].step(xcenters, matched, where="mid", label="matched", color="#184f78")
    axes[1].step(xcenters, detector_fake, where="mid", label="detector/order fake", color="#a33a2b")
    axes[1].step(xcenters, boundary_fake, where="mid", label="boundary fake", color="#8a6d1f")
    axes[1].set_yscale("symlog", linthresh=1.0e-5)
    axes[1].set_xlabel("reco xJgamma")
    axes[1].set_ylabel("weighted yield")
    axes[1].set_title("Matched and fake support")
    axes[1].grid(alpha=0.25)
    axes[1].legend()
    fig.suptitle("Accepted H70 response support — misses, fakes and boundaries")
    save_figure(fig, output / "FAKE_MISS_SUPPORT")

    iterations = np.asarray([row["iterations"] for row in iteration_scan])
    fig, ax = plt.subplots(figsize=(9.7, 5.6))
    ax.semilogy(iterations, [max(row["xj_truth_closure_chi2_ndf"], 1.0e-32) for row in iteration_scan], "o-", label="xJ truth closure")
    ax.semilogy(iterations, [max(row["photon_truth_closure_chi2_ndf"], 1.0e-32) for row in iteration_scan], "s-", label="photon truth closure")
    ax.set_xlabel("Bayesian iterations k")
    ax.set_ylabel("response-only MC closure chi2/ndf")
    ax.set_title("Response-only truth closure, k=2...12")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.text(0.01, 0.01, "This is not the final DATA iteration selection.", fontsize=9, color="#a33a2b")
    save_figure(fig, output / "TRUTH_CLOSURE")

    fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.3))
    axes[0].plot(iterations, [row["xj_refold_chi2_ndf"] for row in iteration_scan], "o-", label="xJ")
    axes[0].plot(iterations, [row["photon_refold_chi2_ndf"] for row in iteration_scan], "s-", label="photon")
    axes[0].set_xlabel("Bayesian iterations k")
    axes[0].set_ylabel("response-only refold chi2/ndf")
    axes[0].set_title("MC refolding")
    axes[0].grid(alpha=0.25)
    axes[0].legend()
    changes = [np.nan if row["xj_relative_l1_change_from_previous"] is None else row["xj_relative_l1_change_from_previous"] for row in iteration_scan]
    axes[1].semilogy(iterations, np.maximum(changes, 1.0e-16), "o-", color="#6043aa")
    axes[1].set_xlabel("Bayesian iterations k")
    axes[1].set_ylabel("relative L1 change from k-1")
    axes[1].set_title("Response-only iteration stability")
    axes[1].grid(alpha=0.25)
    fig.suptitle("Accepted-response iteration diagnostics")
    fig.text(0.01, 0.01, "Final selection requires full DATA refolding and 6000 toys; k=11 is not forced.", fontsize=9, color="#a33a2b")
    save_figure(fig, output / "ITERATION_STABILITY")

    status_figure(output / "PHOTON_REFOLDING_STATUS", "Photon unfolding / refolding",
                  "Accepted response is offline; measured full-stat photon spectrum is not yet authorized.")
    status_figure(output / "XJ_REFOLDING_STATUS", "Joint pTgamma-xJ unfolding / refolding",
                  "K-corrected full-stat detector input does not yet exist.")
    status_figure(output / "FINAL_AUAU_XJGAMMA_STATUS", "Current nominal H70 AuAu xJgamma",
                  "Candidate points are withheld until exact full membership and all chain gates pass.")

    shown = [row for row in historical_rows if int(row["displayed"]) == 1]
    fig, ax = plt.subplots(figsize=(9.6, 5.8))
    hx = np.asarray([float(row["xj_center"]) for row in shown])
    hxe = np.asarray([(float(row["xj_high"]) - float(row["xj_low"])) / 2 for row in shown])
    hy = np.asarray([float(row["historical_value"]) for row in shown])
    hye = np.asarray([float(row["historical_stat_error"]) for row in shown])
    ax.errorbar(hx, hy, xerr=hxe, yerr=hye, fmt="o", color="#b3262e",
                label="THE-89 AuAu diagnostic/noncanonical")
    ax.axhline(0, color="#777777", linewidth=0.8)
    ax.set_xlabel("xJgamma")
    ax.set_ylabel("1/Ngamma dNjet/dxJgamma")
    ax.set_title("Historical AuAu regression comparison only")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.text(0.01, 0.01, "No THE-255 full-stat point is overlaid while the foundation is running.", fontsize=9, color="#a33a2b")
    save_figure(fig, output / "HISTORICAL_AUAU_COMPARISON")

    slide_map = f"""# Slide evidence map — AuAu side of the combined pp+AuAu story

## SLIDE 1 — measurement / selections

CLAIM=The nominal AuAu analysis is frozen to 15-35 GeV photons, 0-20% centrality, H70, SUB1 tower isolation, and anti-kT R=0.4 SUB1 recoil jets with corrected pT>=5 GeV.

NUMERICAL_EVIDENCE=Model SHA `{MODEL_SHA256}`; |z|<10 cm; |eta gamma|, |eta jet|<0.7; DeltaR>=0.4; DeltaPhi>=7pi/8; response fanout 5/7/10/12 GeV.

SOURCE_ARTIFACT=`03_NOMINAL_AUAU_CONTRACT.md`, `01_PROVENANCE.json`.

FIGURE=`PHOTON_RESPONSE.png` as the compact accepted-response selection witness.

CAVEAT=Final DATA membership is not sealed: 46,327 provisional of target 46,581 roots at the read-only observation.

DO_NOT_CLAIM=Do not claim a full-stat nominal DATA result, use H0 direct histograms, or imply the active recovery controller was modified.

## SLIDE 2 — H70 photon ID + ABCD purity

CLAIM=The response-derived exclusive prompt-photon leakage support is offline, but full-stat A/B/C/D and leakage-aware purity are withheld.

NUMERICAL_EVIDENCE=Accepted response catalog has 7,145 inputs and exact Photon12/Photon20 weights; event-gate companions currently cover 27,931 roots, not the unsealed final base.

SOURCE_ARTIFACT=`05_PURITY.csv`, `13_H70_RESPONSE_WEIGHTING_CONTRACT.md`, `02_DATA_COVERAGE.json`.

FIGURE=`ABCD_PURITY_STATUS.png`.

CAVEAT=No DATA count, alpha, purity, or corrected photon denominator in the package is full-stat nominal yet.

DO_NOT_CLAIM=Do not quote the older incomplete-population H70 purity or historical H0 purity as THE-255.

## SLIDE 3 — detector-level recoil corrections

CLAIM=The nominal K template per matched Region-A photon is available, while A, C, ABCD-corrected, absolute K subtraction, and the final detector input remain membership-gated.

NUMERICAL_EVIDENCE=`10_K_TEMPLATE.csv` preserves Tcomb, Ngamma_SIM, K and statistical support in three pTgamma bins and fine xJ bins; `11_K_SUBTRACTION.csv` contains no DATA value.

SOURCE_ARTIFACT=`09_K_COMBINATORIC_CONTRACT.md`, `10_K_TEMPLATE.csv`, `11_K_SUBTRACTION.csv`, `12_K_CORRECTED_RECOIL.csv`.

FIGURE=`DETECTOR_LEVEL_RECOIL_CORRECTIONS_STATUS.png`.

CAVEAT=Absolute K subtraction requires final leakage-corrected Ngamma_DATA.

DO_NOT_CLAIM=Do not show an older A/C/K chain as the nominal full-stat unfolding input.

## SLIDE 4 — response + unfolding

CLAIM=The accepted nominal H70 photon and joint pTgamma-xJ responses, misses, fakes, boundaries, and response-only k=2...12 closure diagnostics are offline.

NUMERICAL_EVIDENCE=Catalog SHA `{RESPONSE_CATALOG_SHA256}`; 7,145 unique ROOTs; `dst_reads=0`; reco photon edges [10,15,20,25,35,40], truth [5,10,15,20,25,35,40]; xJ 0.05 bins through 2.0 plus 2.0-3.0 support.

SOURCE_ARTIFACT=`14_PHOTON_RESPONSE.root`, `15_XJ_RESPONSE.root`, `16_UNFOLDING_DIAGNOSTICS.json`.

FIGURE=`PHOTON_RESPONSE.png`, `XJ_RESPONSE.png`, `FAKE_MISS_SUPPORT.png`, `TRUTH_CLOSURE.png`, `ITERATION_STABILITY.png`.

CAVEAT=DATA refolding and the 6,000-toy candidate scan are not run. The older k=11 result is recorded only as an incomplete-population checkpoint.

DO_NOT_CLAIM=Do not claim k=11 is selected for THE-255, or that response-only MC closure is DATA unfolding closure.

## SLIDE 5 — final AuAu xJgamma

CLAIM=No current full-stat nominal H70 AuAu candidate is released while the foundation is running; THE-89 is preserved only for diagnostic regression.

NUMERICAL_EVIDENCE=Current THE-255 point count = 0; historical displayed diagnostic bins = {len(shown)}.

SOURCE_ARTIFACT=`18_FINAL_AUAU_POINTS.csv`, `19_HISTORICAL_AUAU_COMPARISON.csv`.

FIGURE=`FINAL_AUAU_XJGAMMA_STATUS.png`; optionally `HISTORICAL_AUAU_COMPARISON.png` with its noncanonical label.

CAVEAT=The historical red result is `THE-89 AuAu response/K statistical-correction candidate; diagnostic/noncanonical`.

DO_NOT_CLAIM=Do not label THE-89 nominal/canonical, imply agreement with pp/ATLAS, or substitute it for THE-255.
"""
    write_text(output / "SLIDE_EVIDENCE_MAP.md", slide_map)

    tools = output / "offline_tools"
    tools.mkdir()
    for source in (paths["reducer"], paths["corrections"], Path(__file__).resolve()):
        shutil.copy2(source, tools / source.name)
    shutil.copy2(paths["prior_k11_scan"], tools / "PRIOR_K11_FIXED_CHECKPOINT_SCAN.json")
    shutil.copy2(paths["historical_covariance"], tools / paths["historical_covariance"].name)
    for path in tools.iterdir():
        if path.is_file():
            set_read_only(path)

    readme = f"""# THE-255 AuAu full money-plot offline package

`AUAU_MONEYPLOT_PRESENTATION_STATUS={STATUS}`

The final AuAu base catalog is not sealed and uniform `AuAuEventGateV1` coverage cannot be frozen until the final identity set exists. Therefore this package contains **no incomplete-population nominal DATA values or final points**.

Closed offline now:

- Accepted nominal H70 response catalog, hash and weights.
- Photon and joint pTgamma-xJ response ROOTs for pTjet thresholds 5/7/10/12 GeV, including Sumw2, misses, fakes, boundaries, fake-photon and combinatoric support.
- Nominal per-photon K template.
- Response-only truth-closure/refolding scan for k=2...12.
- Exact nominal contract, correction/unfolding code snapshot, historical THE-89 regression table, slide evidence map, PNG/PDF figures, and hashes.

Open by construction:

- Final catalog-bound TTree to compact-hist reduction.
- Full-stat A/B/C/D, purity/alpha, Region-C correction, absolute K subtraction, measured photon/xJ unfolding, DATA refolding, 6,000-toy iteration selection, photon normalization, and final candidate points.

The ROOT file `DATA_REDUCTION_STATUS.root` is metadata-only and explicitly contains no DATA histogram. The response ROOT files are physical accepted response products, not placeholders.

Offline reproduction sources are under `offline_tools/`. No SDCC access, DST read, pp input, replacement production, recovery mutation, H0 nominal input, or `ACTIVE_DEBUG_LEDGER.md` edit occurred.
"""
    write_text(output / "00_README.md", readme)

    manifest_rows = []
    for path in sorted(output.rglob("*")):
        if path.is_file() and path.name != "PACKAGE_MANIFEST.json":
            manifest_rows.append({
                "path": path.relative_to(output).as_posix(),
                "size_bytes": path.stat().st_size,
                "sha256": sha256(path),
            })
    write_json(output / "PACKAGE_MANIFEST.json", {
        "schema": "THE255OfflinePackageManifestV1",
        "status": STATUS,
        "created_at": utc_now(),
        "file_count_excluding_manifest": len(manifest_rows),
        "files": manifest_rows,
    })
    for path in output.rglob("*"):
        if path.is_file():
            set_read_only(path)
    for path in sorted((p for p in output.rglob("*") if p.is_dir()), reverse=True):
        path.chmod(0o550)
    output.chmod(0o550)
    print(json.dumps({
        "status": STATUS,
        "output": str(output),
        "provisional_auau_roots": observation["base_membership"]["current_provisional_roots"],
        "event_gate_covered_roots": observation["event_gate"]["existing_validated_companion_roots"],
        "response_input_roots": receipt["input_root_count"],
        "file_count": len(manifest_rows) + 1,
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"ERROR: {error}")
        raise SystemExit(2)
