#!/usr/bin/env python3
"""Build THE-219 Slide 1 as a Sam-parity leading-jet xJgamma check.

This is intentionally a reconstruction-level validation slide.  The current
ROOT artifacts contain the event-leading PPG12 photon and one highest-pT recoil
jet per event, but they do not contain a complete leading-jet unfolding
response.  The script therefore never labels the data as unfolded or
particle-level.

The display mirrors the consequential parts of Sam's PPG18 convention:
  * one maximum-pT recoil jet per photon;
  * the pT-dependent xJ turn-on is excluded from the displayed comparison;
  * each curve is shape-normalized over 0.5 < xJgamma < 2.0;
  * 0.1-wide xJ bins are used for the audience-facing view.

No smoothing is applied.  The script writes a full-slide PNG plus numeric,
layout, provenance, and speaker-note sidecars.
"""

from __future__ import annotations

import csv
import json
import math
import os
from pathlib import Path
from typing import Any

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import ROOT


ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

REPO = Path(__file__).resolve().parents[4]
SAM_REFERENCE_PDF = os.environ.get("RJ_SAM_PPG18_PDF")
OUT = REPO / "dataOutput/the219_friday_ppg_20260814/slide01_pp_leading_parity"
OUT.mkdir(parents=True, exist_ok=True)

PNG = OUT / "slide01_pp_leading_xjgamma_sam_parity_reco_candidate.png"
CSV = OUT / "slide01_pp_leading_xjgamma_sam_parity_points.csv"
POINTS = OUT / "slide01_pp_leading_xjgamma_sam_parity_points.json"
MANIFEST = OUT / "slide01_pp_leading_xjgamma_sam_parity_manifest.json"
LAYOUT = OUT / "slide01_pp_leading_xjgamma_sam_parity_layout_nodes.json"
NOTES = OUT / "slide01_pp_leading_xjgamma_sam_parity_speaker_notes.md"

DATA_POINTER = REPO / "dataOutput/current_recoiljets_artifacts/current/pp_data_merged/current.json"
SIM_POINTER = REPO / "dataOutput/current_recoiljets_artifacts/current/pp_sim_photonjet_merged/current.json"

DATA_TOP = "PPG12_scaledtrigger30"
SIM_TOP = "SIM"
RECO_OBJECT = "h_JES3_pT_xJ_alpha_r04"
TRUTH_OBJECT = "h_JES3TruthPure_pT_xJ_alpha_r04"

PT_GROUPS = [(16.0, 20.0), (20.0, 26.0), (26.0, 35.0)]
DISPLAY_XMIN = {(16.0, 20.0): 0.40, (20.0, 26.0): 0.30, (26.0, 35.0): 0.20}
XMAX = 2.0
SHAPE_NORM = (0.5, 2.0)
REBIN = 2  # native 0.05 bins -> audience-facing 0.10 bins


def open_root(path: Path) -> ROOT.TFile:
    source = ROOT.TFile.Open(str(path), "READ")
    if not source or source.IsZombie():
        raise RuntimeError(f"failed to open ROOT file: {path}")
    return source


def get_th3(source: ROOT.TFile, top: str, name: str) -> ROOT.TH3:
    obj = source.Get(f"{top}/{name}")
    if not obj or not obj.InheritsFrom("TH3"):
        raise KeyError(f"missing TH3 {top}/{name} in {source.GetName()}")
    clone = obj.Clone(f"{name}_{top}_the219_clone")
    clone.SetDirectory(0)
    clone.Sumw2()
    return clone


def exact_axis_bins(axis: ROOT.TAxis, lo: float, hi: float) -> list[int]:
    bins = []
    for index in range(1, axis.GetNbins() + 1):
        low = float(axis.GetBinLowEdge(index))
        high = float(axis.GetBinUpEdge(index))
        if low >= lo - 1e-9 and high <= hi + 1e-9:
            bins.append(index)
    width = sum(float(axis.GetBinWidth(index)) for index in bins)
    if not bins or not math.isclose(width, hi - lo, rel_tol=0.0, abs_tol=1e-8):
        raise RuntimeError(f"{lo:g}-{hi:g} GeV is not an exact union of native pT bins")
    if bins != list(range(bins[0], bins[-1] + 1)):
        raise RuntimeError(f"non-contiguous pT bins for {lo:g}-{hi:g} GeV")
    return bins


def project_xj(hist: ROOT.TH3, lo: float, hi: float, tag: str) -> ROOT.TH1:
    bins = exact_axis_bins(hist.GetXaxis(), lo, hi)
    projected = hist.ProjectionY(
        f"xj_{tag}_{int(lo)}_{int(hi)}",
        bins[0],
        bins[-1],
        1,
        hist.GetZaxis().GetNbins(),
        "e",
    )
    projected.SetDirectory(0)
    projected.Sumw2()
    if projected.GetNbinsX() % REBIN:
        raise RuntimeError("native xJ bin count is not divisible by requested rebin factor")
    rebinned = projected.Rebin(REBIN, f"{projected.GetName()}_rebin{REBIN}")
    rebinned.SetDirectory(0)
    return rebinned


def normalize_shape(hist: ROOT.TH1) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    axis = hist.GetXaxis()
    edges = np.array([float(axis.GetBinLowEdge(i)) for i in range(1, axis.GetNbins() + 2)])
    values = np.array([float(hist.GetBinContent(i)) for i in range(1, axis.GetNbins() + 1)])
    errors = np.array([float(hist.GetBinError(i)) for i in range(1, axis.GetNbins() + 1)])
    norm_mask = (edges[:-1] >= SHAPE_NORM[0] - 1e-9) & (edges[1:] <= SHAPE_NORM[1] + 1e-9)
    norm = float(values[norm_mask].sum())
    if not math.isfinite(norm) or norm <= 0:
        raise RuntimeError(f"non-positive shape normalization for {hist.GetName()}")
    # This matches Sam's per-bin shape scaling.  The common normalization
    # denominator is not propagated into the diagonal error bars.
    return edges, values / norm, errors / norm, norm


def panel_payload(
    data_hist: ROOT.TH3,
    reco_hist: ROOT.TH3,
    truth_hist: ROOT.TH3,
    lo: float,
    hi: float,
) -> dict[str, Any]:
    series = {}
    for label, source in (("data", data_hist), ("pythia_reco", reco_hist), ("pythia_truth", truth_hist)):
        projected = project_xj(source, lo, hi, label)
        edges, values, errors, norm = normalize_shape(projected)
        series[label] = {
            "values": values,
            "errors": errors,
            "normalization_0p5_to_2p0": norm,
            "raw_entries": float(projected.GetEntries()),
        }
    edges = normalize_shape(project_xj(data_hist, lo, hi, "data_edges"))[0]
    data = series["data"]["values"]
    data_err = series["data"]["errors"]
    reco = series["pythia_reco"]["values"]
    reco_err = series["pythia_reco"]["errors"]
    ratio = np.divide(data, reco, out=np.full_like(data, np.nan), where=reco > 0)
    ratio_var = np.zeros_like(data)
    valid = (data > 0) & (reco > 0)
    ratio_var[valid] = ratio[valid] ** 2 * (
        (data_err[valid] / data[valid]) ** 2 + (reco_err[valid] / reco[valid]) ** 2
    )
    turn_on = 5.0 / lo
    return {
        "pt": [lo, hi],
        "xj_edges": edges.tolist(),
        "display_xmin": DISPLAY_XMIN[(lo, hi)],
        "conservative_turn_on_5_over_pt_low": turn_on,
        "data": series["data"]["values"].tolist(),
        "data_stat": series["data"]["errors"].tolist(),
        "pythia_reco": series["pythia_reco"]["values"].tolist(),
        "pythia_reco_stat": series["pythia_reco"]["errors"].tolist(),
        "pythia_truth": series["pythia_truth"]["values"].tolist(),
        "pythia_truth_stat": series["pythia_truth"]["errors"].tolist(),
        "data_over_pythia_reco": ratio.tolist(),
        "data_over_pythia_reco_stat": np.sqrt(ratio_var).tolist(),
        "normalizations": {
            label: float(payload["normalization_0p5_to_2p0"])
            for label, payload in series.items()
        },
        "raw_entries": {label: float(payload["raw_entries"]) for label, payload in series.items()},
    }


def px_bbox(left: float, bottom: float, width: float, height: float) -> list[float]:
    return [2560 * left, 1440 * (1 - bottom - height), 2560 * (left + width), 1440 * (1 - bottom)]


def render(panels: list[dict[str, Any]]) -> None:
    mpl.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
        "mathtext.fontset": "stix",
        "axes.linewidth": 1.05,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
    })
    fig = plt.figure(figsize=(16, 9), dpi=160, facecolor="white")
    fig.text(
        0.047,
        0.956,
        r"Matching the leading-jet observable removes the apparent low-$x_{J\gamma}$ spike",
        ha="left",
        va="top",
        fontsize=27,
        fontweight="bold",
        color="#111111",
    )
    fig.text(
        0.049,
        0.896,
        r"Sam-parity detector-level check  |  one highest-$p_T$ recoil jet/event  |  shape normalized on $0.5<x_{J\gamma}<2.0$  |  no smoothing",
        ha="left",
        va="top",
        fontsize=15.2,
        color="#333333",
    )
    fig.text(
        0.953,
        0.913,
        "sPHENIX Internal",
        ha="right",
        va="top",
        fontsize=15.5,
        fontstyle="italic",
        fontweight="bold",
        color="#222222",
    )

    left0, gap, total_w = 0.073, 0.034, 0.871
    panel_w = (total_w - 2 * gap) / 3
    main_bottom, main_height = 0.327, 0.455
    ratio_bottom, ratio_height = 0.182, 0.118
    blue = "#1665a7"
    blue_fill = "#a9cee8"
    magenta = "#c2185b"

    for index, panel in enumerate(panels):
        left = left0 + index * (panel_w + gap)
        ax = fig.add_axes([left, main_bottom, panel_w, main_height])
        rax = fig.add_axes([left, ratio_bottom, panel_w, ratio_height], sharex=ax)
        edges = np.asarray(panel["xj_edges"], dtype=float)
        centers = 0.5 * (edges[:-1] + edges[1:])
        halfwidth = 0.5 * np.diff(edges)
        data = np.asarray(panel["data"], dtype=float)
        data_err = np.asarray(panel["data_stat"], dtype=float)
        reco = np.asarray(panel["pythia_reco"], dtype=float)
        reco_err = np.asarray(panel["pythia_reco_stat"], dtype=float)
        truth = np.asarray(panel["pythia_truth"], dtype=float)
        truth_err = np.asarray(panel["pythia_truth_stat"], dtype=float)
        ratio = np.asarray(panel["data_over_pythia_reco"], dtype=float)
        ratio_err = np.asarray(panel["data_over_pythia_reco_stat"], dtype=float)
        xmin = float(panel["display_xmin"])
        shown = (edges[:-1] >= xmin - 1e-9) & (edges[1:] <= XMAX + 1e-9)

        ax.fill_between(
            centers[shown],
            np.maximum(reco[shown] - reco_err[shown], 0.0),
            reco[shown] + reco_err[shown],
            step="mid",
            color=blue_fill,
            alpha=0.48,
            linewidth=0,
            zorder=1,
        )
        ax.plot(centers[shown], reco[shown], color=blue, linewidth=2.0, drawstyle="steps-mid", zorder=2)
        ax.errorbar(
            centers[shown],
            truth[shown],
            yerr=truth_err[shown],
            xerr=halfwidth[shown],
            fmt="o",
            markerfacecolor="white",
            markeredgecolor=magenta,
            ecolor=magenta,
            color=magenta,
            markersize=5.1,
            capsize=1.7,
            linewidth=1.0,
            zorder=3,
        )
        ax.plot(centers[shown], truth[shown], color=magenta, linewidth=1.35, alpha=0.82, zorder=2)
        ax.errorbar(
            centers[shown],
            data[shown],
            yerr=data_err[shown],
            xerr=halfwidth[shown],
            fmt="o",
            color="black",
            markersize=4.8,
            capsize=2.0,
            linewidth=1.2,
            zorder=4,
        )
        ax.set_xlim(xmin, XMAX)
        ymax = max(float(np.nanmax(data[shown] + data_err[shown])), float(np.nanmax(reco[shown] + reco_err[shown])), float(np.nanmax(truth[shown] + truth_err[shown])))
        ax.set_ylim(0.0, max(0.20, math.ceil(1.13 * ymax / 0.025) * 0.025))
        ax.grid(axis="y", color="#dedede", linewidth=0.65, alpha=0.78)
        ax.tick_params(labelsize=11.2, length=4)
        ax.tick_params(labelbottom=False)
        lo, hi = panel["pt"]
        ax.text(
            0.05,
            0.945,
            rf"${lo:.0f}<p_{{T}}^{{\gamma}}<{hi:.0f}$ GeV",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=15,
            fontweight="bold",
        )
        ax.text(
            0.05,
            0.850,
            rf"$N_{{data}}={panel['raw_entries']['data']:.0f}$  |  display starts at $x_{{J\gamma}}={xmin:.1f}$",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=12.2,
            color="#444444",
        )
        if index == 0:
            ax.set_ylabel(r"shape-normalized $\Delta N/\Delta x_{J\gamma}$", fontsize=16.2, labelpad=8)
            handles = [
                Line2D([0], [0], marker="o", color="black", linestyle="none", markersize=5.3, label="p+p data (reco)"),
                Line2D([0], [0], color=blue, linewidth=2.2, label="PYTHIA 8 (reco)"),
                Line2D([0], [0], marker="o", markerfacecolor="white", markeredgecolor=magenta, color=magenta, linewidth=1.2, markersize=5.4, label="PYTHIA 8 (particle)"),
            ]
            ax.legend(handles=handles, loc="upper right", frameon=False, fontsize=11.1, handletextpad=0.45, borderaxespad=0.3)
        else:
            ax.tick_params(labelleft=False)

        ratio_shown = shown & np.isfinite(ratio) & (reco > 0)
        rax.axhline(1.0, color="#666666", linewidth=1.0, linestyle="--")
        rax.errorbar(
            centers[ratio_shown],
            ratio[ratio_shown],
            yerr=ratio_err[ratio_shown],
            xerr=halfwidth[ratio_shown],
            fmt="o",
            color="black",
            markersize=4.2,
            capsize=1.6,
            linewidth=1.0,
        )
        rax.set_xlim(xmin, XMAX)
        rax.set_ylim(0.35, 1.65)
        rax.set_yticks([0.5, 1.0, 1.5])
        rax.grid(axis="y", color="#e2e2e2", linewidth=0.6)
        rax.tick_params(labelsize=10.4, length=3.5)
        rax.set_xlabel(r"$x_{J\gamma}=p_T^{jet1}/p_T^\gamma$", fontsize=14.2, labelpad=3)
        if index == 0:
            rax.set_ylabel("Data /\nMC reco", fontsize=11.8, labelpad=7)
        else:
            rax.tick_params(labelleft=False)

    callout = mpl.patches.FancyBboxPatch(
        (0.049, 0.014),
        0.902,
        0.075,
        transform=fig.transFigure,
        boxstyle="round,pad=0.007,rounding_size=0.007",
        facecolor="#eaf3fb",
        edgecolor="#6c9fc2",
        linewidth=1.0,
    )
    fig.patches.append(callout)
    fig.text(
        0.5,
        0.058,
        r"The retired candidate used every recoil jet; this view uses one leading jet, explicit turn-on edges, and Sam's shape normalization.",
        ha="center",
        va="center",
        fontsize=16.2,
        fontweight="bold",
        color="#153b57",
    )
    fig.text(
        0.5,
        0.033,
        r"Reco-level parity candidate; $26$-$35$ GeV is honestly statistics-limited. Leading response + covariance remains the particle-level gate.",
        ha="center",
        va="center",
        fontsize=13.2,
        color="#274c66",
    )
    fig.savefig(PNG, dpi=160, facecolor="white")
    plt.close(fig)


def write_layout() -> None:
    nodes: list[dict[str, Any]] = [
        {
            "name": "claim title",
            "kind": "text",
            "role": "title",
            "font_px": 60,
            "bbox": [120, 52, 2200, 118],
            "text": "Matching the leading-jet observable removes the apparent low-xJgamma spike",
        },
        {
            "name": "contract subtitle",
            "kind": "text",
            "role": "audience",
            "font_px": 34,
            "bbox": [125, 142, 2255, 192],
            "text": "Sam-parity detector-level check; one highest-pT recoil jet per event; shape normalized on 0.5<xJgamma<2.0; no smoothing",
        },
        {
            "name": "bottom takeaway",
            "kind": "text",
            "role": "audience",
            "font_px": 36,
            "bbox": [140, 1315, 2420, 1425],
            "text": "The retired candidate used every recoil jet; this view uses one leading jet, explicit turn-on edges, and Sam's shape normalization. The highest pT bin is statistics-limited.",
        },
    ]
    left0, gap, total_w = 0.073, 0.034, 0.871
    panel_w = (total_w - 2 * gap) / 3
    for index in range(3):
        left = left0 + index * (panel_w + gap)
        nodes.append({
            "name": f"result panel {index + 1}",
            "kind": "panel",
            "role": "audience",
            "symmetry_group": "three leading xJ panels",
            "bbox": px_bbox(left, 0.175, panel_w, 0.607),
        })
    LAYOUT.write_text(json.dumps({
        "schema": "slide_layout_nodes_v1",
        "slide": "the219_pp_leading_xjgamma_sam_parity_reco_candidate",
        "title_axis_x": 120,
        "minimum_audience_font_px": 29,
        "minimum_plot_annotation_font_px": 24,
        "minimum_title_font_px": 56,
        "nodes": nodes,
    }, indent=2) + "\n")


def finite_json(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: finite_json(item) for key, item in value.items()}
    if isinstance(value, list):
        return [finite_json(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def write_points(panels: list[dict[str, Any]]) -> None:
    POINTS.write_text(json.dumps(finite_json(panels), indent=2, allow_nan=False) + "\n")
    fields = [
        "pt_low", "pt_high", "xj_low", "xj_high", "data", "data_stat",
        "pythia_reco", "pythia_reco_stat", "pythia_truth", "pythia_truth_stat",
        "data_over_pythia_reco", "data_over_pythia_reco_stat", "displayed",
    ]
    with CSV.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for panel in panels:
            edges = panel["xj_edges"]
            for index, (low, high) in enumerate(zip(edges[:-1], edges[1:])):
                writer.writerow({
                    "pt_low": panel["pt"][0],
                    "pt_high": panel["pt"][1],
                    "xj_low": low,
                    "xj_high": high,
                    "data": panel["data"][index],
                    "data_stat": panel["data_stat"][index],
                    "pythia_reco": panel["pythia_reco"][index],
                    "pythia_reco_stat": panel["pythia_reco_stat"][index],
                    "pythia_truth": panel["pythia_truth"][index],
                    "pythia_truth_stat": panel["pythia_truth_stat"][index],
                    "data_over_pythia_reco": panel["data_over_pythia_reco"][index],
                    "data_over_pythia_reco_stat": panel["data_over_pythia_reco_stat"][index],
                    "displayed": low >= panel["display_xmin"] - 1e-9 and high <= XMAX + 1e-9,
                })


def main() -> None:
    data_pointer = json.loads(DATA_POINTER.read_text())
    sim_pointer = json.loads(SIM_POINTER.read_text())
    data_path = Path(data_pointer["root_paths"][0])
    sim_path = Path(sim_pointer["root_paths"][0])
    data_file = open_root(data_path)
    sim_file = open_root(sim_path)
    try:
        data_hist = get_th3(data_file, DATA_TOP, RECO_OBJECT)
        reco_hist = get_th3(sim_file, SIM_TOP, RECO_OBJECT)
        truth_hist = get_th3(sim_file, SIM_TOP, TRUTH_OBJECT)
        panels = [panel_payload(data_hist, reco_hist, truth_hist, lo, hi) for lo, hi in PT_GROUPS]
    finally:
        data_file.Close()
        sim_file.Close()

    render(panels)
    write_layout()
    write_points(panels)
    MANIFEST.write_text(json.dumps({
        "schema": "the219_pp_leading_xjgamma_sam_parity_v1",
        "status": "reco_level_parity_candidate",
        "slide_png": str(PNG),
        "data_pointer": str(DATA_POINTER),
        "sim_pointer": str(SIM_POINTER),
        "data_current_entry_id": data_pointer.get("current_entry_id"),
        "sim_current_entry_id": sim_pointer.get("current_entry_id"),
        "data_root": str(data_path),
        "sim_root": str(sim_path),
        "objects": {
            "data_reco": f"{DATA_TOP}/{RECO_OBJECT}",
            "pythia_reco": f"{SIM_TOP}/{RECO_OBJECT}",
            "pythia_particle": f"{SIM_TOP}/{TRUTH_OBJECT}",
        },
        "observable_contract": {
            "photon": "event-leading PPG12 isolated+tight reconstructed photon",
            "jet": "one global highest-pT non-overlapping jet, then eta and |dphi|>=7pi/8 veto; no fallback",
            "jet_radius": 0.4,
            "jet_pt_min_GeV": 5.0,
            "eta_abs_max_photon_and_jet": 0.7,
            "pt_groups_GeV": PT_GROUPS,
            "display_xmin": {f"{int(lo)}-{int(hi)}": DISPLAY_XMIN[(lo, hi)] for lo, hi in PT_GROUPS},
            "native_xj_width": 0.05,
            "display_xj_width": 0.10,
            "shape_normalization": list(SHAPE_NORM),
        },
        "claims_allowed": [
            "The selected p+p leading-jet reconstruction-level shapes are smooth after matching the observable contract.",
            "The low-xJ spike in the retired candidate was tied to an inclusive-recoil observable and turn-on display choice, not to a need for smoothing.",
        ],
        "claims_not_allowed": [
            "The p+p data points are unfolded or particle-level.",
            "Full systematic covariance is included.",
            "A data-to-MC JES correction has been applied in this bounded consumer.",
            "The leading-jet particle-level result is final.",
        ],
        "sam_reference": {
            "pdf": SAM_REFERENCE_PDF,
            "sdcc_repo_read_only_inspection": "/sphenix/user/samfred/projects/gammajet",
            "observed_contract": "getmaxjet -> one max-pT jet; shape scale on 0.5-2.0; Gaussian fit begins above minjet/low-photon-edge",
        },
        "local_code_evidence": {
            "leading_fill": "src/RecoilJets.cc: physics-output h_xJ uses recoil1Jet after global-leading selection",
            "inclusive_difference": "src/RecoilJets.cc: h2_unfoldReco_pTgamma_xJ_incl loops over all fiducial recoil jets",
        },
        "sidecars": [str(CSV), str(POINTS), str(LAYOUT), str(NOTES)],
    }, indent=2) + "\n")
    NOTES.write_text(
        "# Slide 1 speaker notes — p+p leading-jet parity\n\n"
        "**One-sentence message:** When we make the same leading-jet observable comparison as Sam, our fixed PPG12 photon definition gives the expected smooth p+p xJgamma shape; the earlier low-x spike belonged to a different inclusive-recoil observable and its turn-on region.\n\n"
        "- These black points are reconstruction-level data, not unfolded data.\n"
        "- We use exactly one highest-pT recoil jet per event, then require the R=0.4 eta and back-to-back cuts.\n"
        "- We rebin only from native 0.05 to 0.10 and shape-normalize on 0.5<xJgamma<2.0, matching Sam's audience-facing convention. No smoothing is applied.\n"
        "- The magenta particle-level PYTHIA curve shows the target direction, not an unfolded measurement.\n"
        "- The next physics gate is a genuine leading-jet response with data-to-MC JES handling and full systematic covariance. Until that exists, do not label this as the money particle-level result.\n"
        "- Keep the inclusive per-photon recoil yield as a separate Au+Au-target observable because it preserves recoil multiplicity for IAA; do not judge it by Gaussianity.\n"
    )
    print(PNG)


if __name__ == "__main__":
    main()
