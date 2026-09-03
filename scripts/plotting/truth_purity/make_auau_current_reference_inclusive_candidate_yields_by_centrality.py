#!/usr/bin/env python3
"""Plot inclusive Au+Au ABCD-candidate yields in three centrality panels.

This is an offline consumer of the lossless purity-factorial candidate blocks.
It never rereads Trees, performs event-leading arbitration, or derives/applies
centrality weights.  Each input block must already carry the complete canonical
producer -> source-stitch -> centrality analysis weight and its receipt binding.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile
from typing import Iterable, Mapping

import numpy as np


REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "scripts"))

from data_prep.recoiljets.auau_centrality_weight_contract import APPLICATION_ORDER  # noqa: E402
from data_prep.recoiljets.auau_embedded_inclusive_schema10_weighting import (  # noqa: E402
    load_analysis_weight_receipt,
    load_downstream_artifact_receipt,
    write_downstream_artifact_receipt,
)
from plotting.plot_label_contract import (  # noqa: E402
    FULL_ACCEPTED_SIMULATION_SCOPE,
    PlotLabelContract,
    canonical_cut_lines,
)
T70_INTERCEPT = 0.6927385742113304
T70_SLOPE = 0.0010476808935335573
T80_INTERCEPT = 0.6071724560222583
T80_SLOPE = 0.0013019695025559648
T90_INTERCEPT = 0.4634577166843198
T90_SLOPE = 0.0015466253667456042
ISO_INTERCEPT_GEV = 7.569999970933129
ISO_SLOPE_GEV = -0.06579999928964807

CENTRALITY_BINS = ((0.0, 20.0), (20.0, 50.0), (50.0, 80.0))
PT_EDGES = np.arange(15.0, 36.0, 1.0)
REGIONS = ("A", "B", "C", "D")
REGION_STYLES = {
    "A": ("#2774AE", "s", "Region A: tight, isolated"),
    "B": ("#E68613", "D", "Region B: tight, non-isolated"),
    "C": ("#2E9F5B", "o", "Region C: non-tight, isolated"),
    "D": ("#BF1E2E", "^", "Region D: non-tight, non-isolated"),
}
DISPLAY_OFFSETS_GEV = {"A": -0.18, "B": -0.06, "C": 0.06, "D": 0.18}
REQUIRED_POPULATION = "ALL_BASE_SELECTED_CANDIDATES_NO_EVENT_LEADER_ARBITRATION"
REQUIRED_WEIGHT_STATE = "COMPLETE_PRODUCER_STITCH_CENTRALITY"
CANDIDATE_BLOCK_SCHEMA = "THE259PurityFactorialCandidateBlockV1"
CANDIDATE_BLOCK_STATUS = "COMPLETE_EXISTING_RICH_TREE_CANDIDATE_BLOCK"


@dataclass(frozen=True)
class CandidateYieldRecord:
    """Validated candidate fields consumed by this downstream plot."""

    ptgamma: float
    centrality: float
    score: float
    isolation: float
    weight: float


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_candidate_block(
    path: Path,
) -> tuple[dict[str, object], tuple[CandidateYieldRecord, ...]]:
    """Read a candidate block without relying on private workspace modules.

    The full persisted layout is checked even though this plot consumes only
    five scalar fields. This keeps malformed or semantically incomplete blocks
    from being silently accepted by the public downstream consumer.
    """
    with np.load(path, allow_pickle=False) as payload:
        metadata = json.loads(payload["metadata_json_utf8"].tobytes().decode())
        if not isinstance(metadata, dict):
            raise ValueError("candidate-block metadata is not a JSON object")
        if (
            metadata.get("schema") != CANDIDATE_BLOCK_SCHEMA
            or metadata.get("status") != CANDIDATE_BLOCK_STATUS
        ):
            raise ValueError("candidate-block schema/status differs")
        count = int(metadata["candidate_count"])
        if count < 0:
            raise ValueError("candidate-block count is negative")
        one_d = (
            "prompt_valid", "ptgamma", "eta", "centrality", "vertex_z", "score",
            "isolation", "isolation_threshold", "isolation_sideband_threshold",
            "weight", "active_preselection_state", "prompt_match_metric",
        )
        if any(len(payload[name]) != count for name in one_d):
            raise ValueError("candidate-block vector length differs")
        if (
            payload["event_ids"].shape != (count, 2)
            or payload["candidate_ids"].shape != (count, 2)
            or payload["prompt_ids"].shape != (count, 2)
        ):
            raise ValueError("candidate identity shape differs")
        feature_count = int(metadata["feature_count"])
        if feature_count < 0 or payload["ordered_features"].shape != (count, feature_count):
            raise ValueError("feature matrix shape differs")
        offsets = payload["xj_offsets"]
        xj_values = payload["xj_values"]
        if (
            len(offsets) != count + 1
            or int(offsets[0]) != 0
            or int(offsets[-1]) != len(xj_values)
            or np.any(offsets[1:] < offsets[:-1])
        ):
            raise ValueError("xJ ragged offsets differ")
        if not np.all(np.isfinite(xj_values)) or not np.all(
            np.isfinite(payload["ordered_features"])
        ):
            raise ValueError("candidate block contains non-finite vector values")

        source_origin = metadata.get("source_origin")
        if source_origin not in {"data", "prompt", "inclusive"}:
            raise ValueError("candidate block has an unsupported source origin")
        if source_origin != "data" and (
            not metadata.get("sample_id") or not metadata.get("source_occurrence_id")
        ):
            raise ValueError("SIM candidate block lacks sample/source-occurrence identity")

        numeric_names = (
            "ptgamma", "eta", "centrality", "vertex_z", "score", "isolation",
            "isolation_threshold", "isolation_sideband_threshold", "weight",
        )
        if any(not np.all(np.isfinite(payload[name])) for name in numeric_names):
            raise ValueError("candidate block contains non-finite scalar values")
        if np.any(payload["weight"] < 0.0):
            raise ValueError("candidate weight must be nonnegative")
        if np.any(payload["isolation_sideband_threshold"] < payload["isolation_threshold"]):
            raise ValueError("isolation sideband threshold is below signal threshold")
        prompt_valid = payload["prompt_valid"].astype(bool)
        if np.any(prompt_valid & ~np.isfinite(payload["prompt_match_metric"])):
            raise ValueError("prompt-labelled candidate requires a finite match metric")

        rows = tuple(
            CandidateYieldRecord(
                ptgamma=float(payload["ptgamma"][index]),
                centrality=float(payload["centrality"][index]),
                score=float(payload["score"][index]),
                isolation=float(payload["isolation"][index]),
                weight=float(payload["weight"][index]),
            )
            for index in range(count)
        )
        return metadata, rows


def classify_current_reference(score: float, isolation: float, centrality: float) -> str | None:
    """Return the strict CURRENT_REFERENCE ABCD region for one candidate."""
    t70 = T70_INTERCEPT + T70_SLOPE * centrality
    t80 = T80_INTERCEPT + T80_SLOPE * centrality
    t90 = T90_INTERCEPT + T90_SLOPE * centrality
    iso = ISO_INTERCEPT_GEV + ISO_SLOPE_GEV * centrality
    tight = score > t70
    nontight = t90 < score < t80
    isolated = isolation < iso
    nonisolated = isolation > iso
    if tight and isolated:
        return "A"
    if tight and nonisolated:
        return "B"
    if nontight and isolated:
        return "C"
    if nontight and nonisolated:
        return "D"
    return None


def validate_candidate_block_metadata(
    metadata: Mapping[str, object],
    *,
    analysis_weight_receipt_sha256: str,
) -> None:
    """Fail closed on population, centrality support, or weight provenance drift."""
    if metadata.get("system") != "auau" or metadata.get("source_origin") != "inclusive":
        raise ValueError("plot accepts only Au+Au inclusive-source candidate blocks")
    if metadata.get("candidate_population") != REQUIRED_POPULATION:
        raise ValueError("candidate block is not the inclusive pre-ABCD candidate universe")
    support = metadata.get("candidate_centrality_range")
    if not isinstance(support, Mapping):
        raise ValueError("candidate block lacks explicit centrality support")
    if (
        float(support.get("minimum_inclusive", np.nan)) != 0.0
        or float(support.get("maximum_exclusive", np.nan)) != 80.0
    ):
        raise ValueError("candidate block must be materialized over exact [0,80) centrality")
    provenance = metadata.get("analysis_weight_provenance")
    if not isinstance(provenance, Mapping):
        raise ValueError("candidate block lacks upstream analysis-weight provenance")
    if provenance.get("state") != REQUIRED_WEIGHT_STATE:
        raise ValueError("candidate weights are not certified producer+stitch+centrality weights")
    if provenance.get("receipt_sha256") != analysis_weight_receipt_sha256:
        raise ValueError("candidate weight receipt differs from requested canonical receipt")
    if tuple(provenance.get("application_order", ())) != tuple(APPLICATION_ORDER):
        raise ValueError("candidate analysis-weight application order differs")


def centrality_bin_index(centrality: float) -> int | None:
    for index, (minimum, maximum) in enumerate(CENTRALITY_BINS):
        if minimum <= centrality < maximum:
            return index
    return None


def accumulate_candidate_blocks(
    paths: Iterable[Path],
    *,
    analysis_weight_receipt_sha256: str,
) -> tuple[np.ndarray, np.ndarray, dict[str, object]]:
    sumw = np.zeros((len(CENTRALITY_BINS), len(PT_EDGES) - 1, len(REGIONS)))
    sumw2 = np.zeros_like(sumw)
    occurrences: set[str] = set()
    candidate_count = 0
    classified_count = 0
    paths = tuple(sorted(paths))
    if not paths:
        raise ValueError("no candidate blocks supplied")
    for path in paths:
        metadata, records = read_candidate_block(path)
        validate_candidate_block_metadata(
            metadata,
            analysis_weight_receipt_sha256=analysis_weight_receipt_sha256,
        )
        occurrence = str(metadata["source_occurrence_id"])
        if occurrence in occurrences:
            raise ValueError(f"duplicate source occurrence: {occurrence}")
        occurrences.add(occurrence)
        for record in records:
            candidate_count += 1
            cent_index = centrality_bin_index(record.centrality)
            pt_index = int(np.searchsorted(PT_EDGES, record.ptgamma, side="right") - 1)
            if cent_index is None or not 0 <= pt_index < len(PT_EDGES) - 1:
                continue
            region = classify_current_reference(record.score, record.isolation, record.centrality)
            if region is None:
                continue
            region_index = REGIONS.index(region)
            sumw[cent_index, pt_index, region_index] += record.weight
            sumw2[cent_index, pt_index, region_index] += record.weight * record.weight
            classified_count += 1
    return sumw, sumw2, {
        "candidate_block_count": len(paths),
        "source_occurrence_count": len(occurrences),
        "candidate_count": candidate_count,
        "classified_candidate_count": classified_count,
        "duplicate_source_occurrences": 0,
    }


def write_csv(path: Path, sumw: np.ndarray, sumw2: np.ndarray) -> None:
    with path.open("x", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=(
            "centrality_min", "centrality_max", "pt_low", "pt_high", "region",
            "weighted_inclusive_candidate_yield", "statistical_uncertainty",
        ))
        writer.writeheader()
        for ci, (cent_min, cent_max) in enumerate(CENTRALITY_BINS):
            for pi, (pt_low, pt_high) in enumerate(zip(PT_EDGES[:-1], PT_EDGES[1:])):
                for ri, region in enumerate(REGIONS):
                    writer.writerow({
                        "centrality_min": cent_min,
                        "centrality_max": cent_max,
                        "pt_low": pt_low,
                        "pt_high": pt_high,
                        "region": region,
                        "weighted_inclusive_candidate_yield": sumw[ci, pi, ri],
                        "statistical_uncertainty": np.sqrt(sumw2[ci, pi, ri]),
                    })


def render(path: Path, sumw: np.ndarray, sumw2: np.ndarray) -> None:
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "font.family": "sans-serif", "font.size": 13, "axes.linewidth": 1.1,
        "xtick.direction": "in", "ytick.direction": "in",
        "xtick.top": True, "ytick.right": True,
    })
    fig, axes = plt.subplots(1, 3, figsize=(17.2, 6.4), sharex=True, sharey=True)
    fig.subplots_adjust(left=0.075, right=0.985, bottom=0.13, top=0.94, wspace=0.055)
    centers = 0.5 * (PT_EDGES[:-1] + PT_EDGES[1:])
    positive = sumw[sumw > 0.0]
    ymin = max(float(np.min(positive)) * 0.35, 1.0e-6) if positive.size else 1.0e-3
    ymax = float(np.max(sumw)) * 6.0 if positive.size else 1.0
    for ci, (axis, (cent_min, cent_max)) in enumerate(zip(axes, CENTRALITY_BINS)):
        for ri, region in enumerate(REGIONS):
            color, marker, label = REGION_STYLES[region]
            axis.errorbar(
                centers + DISPLAY_OFFSETS_GEV[region], sumw[ci, :, ri],
                yerr=np.sqrt(sumw2[ci, :, ri]), fmt=marker, linestyle="none",
                color=color, markerfacecolor="white" if region == "A" else color,
                markeredgecolor=color, markersize=5.8, capsize=2.0,
                label=label if ci == 0 else None,
            )
        axis.set_yscale("log")
        axis.set_xlim(15.0, 35.0)
        axis.set_ylim(ymin, ymax)
        axis.set_title(f"Au+Au {cent_min:.0f}–{cent_max:.0f}%", fontweight="bold")
        axis.set_xlabel(r"$E_T^{\gamma,\mathrm{reco}}$ [GeV]")
        axis.grid(axis="y", which="major", color="0.88", linewidth=0.6)
    axes[0].set_ylabel("Weighted inclusive photon-candidate yield per 1 GeV bin")
    axes[0].text(0.04, 0.95, r"$\it{sPHENIX}$  Internal", transform=axes[0].transAxes,
                 ha="left", va="top", fontsize=17)
    axes[0].text(0.04, 0.875, r"$\sqrt{s_{NN}}=200$ GeV", transform=axes[0].transAxes,
                 ha="left", va="top")
    axes[2].text(0.96, 0.95, "Embedded Inclusive+jet 12/20/30/40",
                 transform=axes[2].transAxes, ha="right", va="top", fontweight="bold")
    axes[0].legend(loc="lower left", frameon=False, fontsize=10.5)
    cut_text = (
        "Inclusive candidates; no event-leading arbitration\n"
        r"tight: $s_{\rm BDT}>0.692739+0.00104768c$" "\n"
        r"non-tight: $0.463458+0.00154663c<s_{\rm BDT}<0.607172+0.00130197c$" "\n"
        r"isolated: $E_T^{\rm iso}<7.57000-0.0658000c$ GeV; non-isolated: above same boundary"
    )
    axes[1].text(0.5, 0.03, cut_text, transform=axes[1].transAxes,
                 ha="center", va="bottom", fontsize=9.5,
                 bbox={"facecolor": "white", "edgecolor": "0.7", "alpha": 0.94})
    descriptor, temporary_name = tempfile.mkstemp(
        dir=path.parent, prefix=f".{path.name}.", suffix=".tmp.png"
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        fig.savefig(temporary, format="png", dpi=180, bbox_inches="tight")
        os.link(temporary, path)
    finally:
        plt.close(fig)
        temporary.unlink(missing_ok=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate-root", type=Path, required=True)
    parser.add_argument("--analysis-weight-receipt", type=Path, required=True)
    parser.add_argument("--expected-dependency-fingerprint", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    receipt = load_analysis_weight_receipt(
        args.analysis_weight_receipt,
        expected_dependency_fingerprint=args.expected_dependency_fingerprint,
        verify_files=True,
    )
    receipt_sha256 = sha256_file(args.analysis_weight_receipt)
    block_paths = tuple(sorted(args.candidate_root.rglob("*.npz")))
    cut_lines = canonical_cut_lines((
        "Inclusive candidates; no event-leading arbitration",
        "tight: s_BDT > 0.692739 + 0.00104768 c",
        "non-tight: 0.463458 + 0.00154663 c < s_BDT < 0.607172 + 0.00130197 c",
        "isolated: E_T^iso < 7.57000 - 0.0658000 c GeV; non-isolated: above same boundary",
    ))
    label_contract = PlotLabelContract(
        system="auau", energy_label=r"sqrt{s_{NN}} = 200 GeV",
        centrality_label="0–80% (panels: 0–20%, 20–50%, 50–80%)",
        sample_label="Embedded Inclusive+jet 12/20/30/40",
        sample_kind="simulation", plot_kind="physics",
        simulation_scope=FULL_ACCEPTED_SIMULATION_SCOPE,
        simulation_family="inclusive_background",
        cut_lines=cut_lines, input_paths=tuple(str(path) for path in block_paths),
        analysis_weight_receipt_path=str(args.analysis_weight_receipt),
        analysis_weight_dependency_fingerprint=args.expected_dependency_fingerprint,
    )
    sumw, sumw2, audit = accumulate_candidate_blocks(
        block_paths, analysis_weight_receipt_sha256=receipt_sha256,
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    png = args.output_dir / "auau_current_reference_inclusive_abcd_yields_by_centrality.png"
    csv_path = args.output_dir / "auau_current_reference_inclusive_abcd_yields_by_centrality.csv"
    audit_path = args.output_dir / "AUAU_CURRENT_REFERENCE_INCLUSIVE_ABCD_YIELDS_AUDIT.json"
    downstream_receipt_path = args.output_dir / (
        "AUAU_CURRENT_REFERENCE_INCLUSIVE_ABCD_YIELDS_DOWNSTREAM_RECEIPT.json"
    )
    targets = (png, csv_path, audit_path, downstream_receipt_path)
    existing = [str(path.resolve()) for path in targets if path.exists()]
    if existing:
        raise FileExistsError(
            "refusing to overwrite existing nominal output artifacts: {}".format(existing)
        )
    audit.update({
        "schema": "AuAuCurrentReferenceInclusiveABCDYieldsAuditV1",
        "status": "PASS",
        "candidate_population": REQUIRED_POPULATION,
        "analysis_weight_state": REQUIRED_WEIGHT_STATE,
        "analysis_weight_receipt_sha256": receipt_sha256,
        "analysis_weight_dependency_fingerprint": args.expected_dependency_fingerprint,
        "analysis_weight_contract_fingerprint": receipt["centrality_contract"][
            "contract_fingerprint"
        ],
        "candidate_blocks": [str(path.resolve()) for path in block_paths],
        "centrality_bins": [list(values) for values in CENTRALITY_BINS],
        "pt_edges": PT_EDGES.tolist(),
        "region_order": list(REGIONS),
        "display_offsets_gev": DISPLAY_OFFSETS_GEV,
        "output_png": str(png.resolve()),
        "output_csv": str(csv_path.resolve()),
        "downstream_artifact_receipt_path": str(downstream_receipt_path.resolve()),
        "generator": {
            "path": str(Path(__file__).resolve()),
            "sha256": sha256_file(Path(__file__).resolve()),
            "size_bytes": Path(__file__).resolve().stat().st_size,
        },
        "cut_constants": {
            "t70": [T70_INTERCEPT, T70_SLOPE],
            "t80": [T80_INTERCEPT, T80_SLOPE],
            "t90": [T90_INTERCEPT, T90_SLOPE],
            "isolation_gev": [ISO_INTERCEPT_GEV, ISO_SLOPE_GEV],
        },
    })
    label_contract.apply_to_audit(audit)
    created: list[Path] = []
    try:
        render(png, sumw, sumw2)
        created.append(png)
        write_csv(csv_path, sumw, sumw2)
        created.append(csv_path)
        with audit_path.open("x", encoding="utf-8") as stream:
            json.dump(audit, stream, indent=2, sort_keys=True, allow_nan=False)
            stream.write("\n")
        created.append(audit_path)
        write_downstream_artifact_receipt(
            downstream_receipt_path,
            analysis_weight_receipt_path=args.analysis_weight_receipt,
            expected_dependency_fingerprint=args.expected_dependency_fingerprint,
            audit_path=audit_path,
            artifacts={"plot_png": png, "yield_table_csv": csv_path},
            input_artifacts=block_paths,
        )
        created.append(downstream_receipt_path)
        load_downstream_artifact_receipt(
            downstream_receipt_path,
            expected_dependency_fingerprint=args.expected_dependency_fingerprint,
            expected_audit_path=audit_path,
            expected_artifacts={"plot_png": png, "yield_table_csv": csv_path},
            expected_inputs=block_paths,
            verify_files=True,
        )
    except Exception:
        for path in reversed(created):
            path.unlink(missing_ok=True)
        raise
    print(json.dumps({
        "status": "PASS",
        "png": str(png),
        "csv": str(csv_path),
        "audit": str(audit_path),
        "downstream_receipt": str(downstream_receipt_path),
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
