#!/usr/bin/env python3
"""Make a PPG12 Fig. 37-style SIM half-closure iteration-stability diagnostic.

This is intentionally a local plotting-only consumer of the canonical current
photon+jet SIM ROOT.  The response is trained with the deterministic first
half and unfolded against the statistically independent second half.  It is
therefore a closure diagnostic for the histogram family added for Fig. 37, not
a replacement for the final pp-data Fig. 37 result.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import ROOT


REPO = Path(__file__).resolve().parents[3]
CURRENT_POINTER = REPO / "dataOutput/current_recoiljets_artifacts/current/pp_sim_photonjet_merged/current.json"
DEFAULT_OUTDIR = (
    REPO
    / "dataOutput/ppg12Parity/the97_ppg12_final_parity_full_20260709_2230"
    / "fig37_sim_halfclosure_20260713"
)

H_RECO_HALF = "SIM/h_ppg12_fig37_input_reco_half_pT_0"
H_TRUTH_HALF = "SIM/h_ppg12_fig37_input_truth_half_pT_0"
H_RESPONSE_HALF = "SIM/h2_ppg12_fig37_input_response_half_reco_truth_0"
H_RECO_SECOND = "SIM/h_ppg12_fig37_input_reco_secondhalf_pT_0"
H_TRUTH_SECOND = "SIM/h_ppg12_fig37_input_truth_secondhalf_pT_0"

# The explicit Fig. 37 family is independently binned through 36 GeV.  Use
# every bin shared by its reco and truth spectra, excluding only truth's
# 8--10 support bin and 36--45 overflow tail.  Do not substitute the separate
# production-unfolding 26--35 binning here.
ANALYSIS_BINS = [(10, 12), (12, 14), (14, 16), (16, 18), (18, 20), (20, 22),
                 (22, 24), (24, 26), (26, 28), (28, 32), (32, 36)]


def detach(obj: ROOT.TObject, name: str) -> ROOT.TObject:
    clone = obj.Clone(name)
    if not clone:
        raise RuntimeError(f"failed to clone {obj.GetName()}")
    if hasattr(clone, "SetDirectory"):
        clone.SetDirectory(0)
    return clone


def find_matching_bins(hist: ROOT.TH1) -> list[int]:
    axis = hist.GetXaxis()
    selected: list[int] = []
    for lo, hi in ANALYSIS_BINS:
        match = None
        for ibin in range(1, axis.GetNbins() + 1):
            if abs(axis.GetBinLowEdge(ibin) - lo) < 1e-9 and abs(axis.GetBinUpEdge(ibin) - hi) < 1e-9:
                match = ibin
                break
        if match is None:
            raise RuntimeError(f"missing analysis bin {lo}-{hi} in {hist.GetName()}")
        selected.append(match)
    return selected


def relative_stat_and_change(
    unfolded: ROOT.TH1,
    previous: ROOT.TH1,
    covariance: ROOT.TMatrixD,
    selected_bins: list[int],
) -> tuple[float, float]:
    """Match BuildPhotonIterScan in AnalyzeRecoilJets_RooUnfoldPipeline.cpp."""
    sum_diag_cov = 0.0
    sum_v2_for_stat = 0.0
    nrows, ncols = covariance.GetNrows(), covariance.GetNcols()
    for hbin in selected_bins:
        index = hbin - 1
        if index < 0 or index >= nrows or index >= ncols:
            continue
        cii = float(covariance[index][index])
        if not (math.isfinite(cii) and cii > 0.0):
            continue
        value = float(unfolded.GetBinContent(hbin))
        if not math.isfinite(value):
            continue
        sum_diag_cov += cii
        sum_v2_for_stat += value * value
    rel_stat = math.sqrt(max(0.0, sum_diag_cov / sum_v2_for_stat)) if sum_v2_for_stat > 0.0 else 0.0

    num2 = 0.0
    den2 = 0.0
    for hbin in selected_bins:
        value = float(unfolded.GetBinContent(hbin))
        previous_value = float(previous.GetBinContent(hbin))
        num2 += (value - previous_value) ** 2
        den2 += value * value
    rel_change = math.sqrt(num2 / den2) if den2 > 0.0 else 0.0
    return (rel_stat if math.isfinite(rel_stat) else 0.0, rel_change if math.isfinite(rel_change) else 0.0)


def run_scan(reco_half: ROOT.TH1, truth_half: ROOT.TH1, response_matrix: ROOT.TH2, reco_second: ROOT.TH1,
             max_iterations: int, toys: int) -> tuple[list[dict[str, float]], int, str]:
    selected_bins = find_matching_bins(truth_half)
    response = ROOT.RooUnfoldResponse(reco_half, truth_half, response_matrix, "the97_fig37_half_response", "the97_fig37_half_response")
    previous = detach(truth_half, "the97_fig37_iteration_zero")
    previous.Reset("ICES")
    for ibin in range(1, previous.GetNbinsX() + 1):
        reco_bin = reco_second.GetXaxis().FindBin(previous.GetXaxis().GetBinCenter(ibin))
        if 1 <= reco_bin <= reco_second.GetNbinsX():
            previous.SetBinContent(ibin, reco_second.GetBinContent(reco_bin))
            previous.SetBinError(ibin, reco_second.GetBinError(reco_bin))

    rows: list[dict[str, float]] = []
    best_iteration = -1
    best_score = float("inf")
    best_reason = "minimum regularized stability score"
    stable_chosen = False
    for iteration in range(1, max_iterations + 1):
        unfolding = ROOT.RooUnfoldBayes(response, reco_second, iteration)
        unfolding.SetVerbose(0)
        unfolding.SetNToys(toys)
        unfolded = unfolding.Hreco(ROOT.RooUnfold.kCovToy)
        covariance = unfolding.Ereco(ROOT.RooUnfold.kCovToy)
        if not unfolded:
            raise RuntimeError(f"RooUnfold returned no spectrum at iteration {iteration}")
        unfolded = detach(unfolded, f"the97_fig37_unfolded_iter{iteration}")
        rel_stat, rel_change = relative_stat_and_change(unfolded, previous, covariance, selected_bins)
        stable_threshold = max(rel_stat, 0.025)
        stable_enough = iteration >= 3 and rel_change <= stable_threshold
        iteration_penalty = 0.010 * max(0, iteration - 3)
        score = math.sqrt(rel_stat * rel_stat + rel_change * rel_change) + iteration_penalty
        rows.append({
            "iteration": float(iteration),
            "total_relative_stat_uncertainty": rel_stat,
            "total_relative_deviation": rel_change,
            "stat_plus_dev_quadrature": score,
        })
        if stable_enough and not stable_chosen:
            best_iteration = iteration
            best_score = score
            best_reason = (
                f"earliest stable iteration: relChange {rel_change:.6g} <= "
                f"max(relStat {rel_stat:.6g}, 0.025)"
            )
            stable_chosen = True
        elif not stable_chosen and score < best_score:
            best_iteration = iteration
            best_score = score
        previous = unfolded
    return rows, best_iteration, best_reason


def draw(rows: list[dict[str, float]], out_png: Path, best_iteration: int) -> None:
    mpl.rcParams.update({
        "font.family": "DejaVu Sans",
        "mathtext.fontset": "dejavusans",
        "axes.linewidth": 1.1,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
    })
    iteration = np.array([r["iteration"] for r in rows])
    stat = np.array([r["total_relative_stat_uncertainty"] for r in rows])
    dev = np.array([r["total_relative_deviation"] for r in rows])
    quad = np.array([r["stat_plus_dev_quadrature"] for r in rows])
    ymax = max(0.12, 1.18 * float(np.max(np.concatenate([stat, dev, quad]))))

    fig, ax = plt.subplots(figsize=(6.25, 5.7), dpi=190)
    ax.plot(iteration, stat, "o", color="black", ms=4.5, label="total relative stat. uncertainty")
    ax.plot(iteration, dev, "o", color="#1f4ae0", ms=4.5, label="total relative deviation")
    ax.plot(iteration, quad, "o", color="#e72e25", ms=4.5, label="stat. + dev. + iteration penalty")
    ax.set(xlim=(0.0, 10.5), ylim=(0.0, ymax), xlabel="Iteration", ylabel=r"$\sqrt{\delta}$")
    ax.set_xticks(np.arange(0, 11, 1))
    ax.minorticks_on()
    ax.tick_params(which="major", length=6, labelsize=10)
    ax.tick_params(which="minor", length=3)
    ax.text(0.51, 0.965, r"$\it{\bf{sPHENIX}}$ Internal", transform=ax.transAxes, ha="left", va="top", fontsize=11)
    ax.text(0.51, 0.905, r"$p{+}p\ \sqrt{s}=200$ GeV", transform=ax.transAxes, ha="left", va="top", fontsize=10)
    ax.text(0.51, 0.850, r"SIM half-closure: first half $\rightarrow$ second half", transform=ax.transAxes, ha="left", va="top", fontsize=8.3)
    ax.text(0.51, 0.800, r"Fig. 37 inputs; diagnostic, not final data", transform=ax.transAxes, ha="left", va="top", fontsize=8.2, color="#8b1a1a")
    ax.legend(frameon=False, fontsize=8.7, loc="center right", bbox_to_anchor=(0.98, 0.49), handletextpad=0.4)
    if best_iteration > 0:
        ax.axvline(best_iteration, color="0.45", lw=0.9, ls="--", zorder=0)
        ax.text(best_iteration + 0.12, 0.08 * ymax, f"selected: {best_iteration}", color="0.3", fontsize=8.3)
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--iterations", type=int, default=10)
    parser.add_argument("--toys", type=int, default=120, help="Matches the existing RecoilJets Fig. 37 scan setting.")
    args = parser.parse_args()

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    if ROOT.gSystem.Load("libRooUnfold") < 0:
        raise RuntimeError("could not load libRooUnfold")
    current = json.loads(CURRENT_POINTER.read_text())
    roots = current.get("root_paths", [])
    if len(roots) != 1:
        raise RuntimeError(f"expected one canonical SIM root in {CURRENT_POINTER}, found {len(roots)}")
    root_path = Path(roots[0])
    root_sha256 = hashlib.sha256(root_path.read_bytes()).hexdigest()
    source = ROOT.TFile.Open(str(root_path), "READ")
    if not source or source.IsZombie():
        raise RuntimeError(f"could not read {root_path}")
    required = [H_RECO_HALF, H_TRUTH_HALF, H_RESPONSE_HALF, H_RECO_SECOND, H_TRUTH_SECOND]
    objects = {}
    for name in required:
        obj = source.Get(name)
        if not obj:
            raise KeyError(f"missing required histogram {name}")
        objects[name] = detach(obj, "the97_" + name.replace("/", "_"))
    source.Close()

    rows, best_iteration, best_reason = run_scan(
        objects[H_RECO_HALF], objects[H_TRUTH_HALF], objects[H_RESPONSE_HALF], objects[H_RECO_SECOND],
        args.iterations, args.toys,
    )
    outdir = args.out_dir
    png_path = outdir / "the97_fig37_sim_halfclosure_iteration_stability.png"
    csv_path = outdir / "the97_fig37_sim_halfclosure_iteration_stability.csv"
    manifest_path = outdir / "the97_fig37_sim_halfclosure_iteration_stability_manifest.json"
    draw(rows, png_path, best_iteration)
    outdir.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    manifest = {
        "plot": str(png_path),
        "csv": str(csv_path),
        "script": str(Path(__file__).resolve()),
        "current_pointer": str(CURRENT_POINTER),
        "campaign_tag": current.get("campaign_tag"),
        "canonical_source_root": str(root_path),
        "canonical_source_sha256": root_sha256,
        "required_input_histograms": required,
        "calculation": "Exact BuildPhotonIterScan formula from AnalyzeRecoilJets_RooUnfoldPipeline.cpp, with first-half response/reco/truth and second-half measured reco.",
        "analysis_bins_gev": [list(pair) for pair in ANALYSIS_BINS],
        "iterations": args.iterations,
        "toys_per_iteration": args.toys,
        "best_iteration": best_iteration,
        "best_iteration_reason": best_reason,
        "scope": "SIM half-closure diagnostic only; not final pp-data PPG12 Fig. 37.",
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(png_path)
    print(csv_path)
    print(manifest_path)
    print(json.dumps({"best_iteration": best_iteration, "reason": best_reason, "rows": rows}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
