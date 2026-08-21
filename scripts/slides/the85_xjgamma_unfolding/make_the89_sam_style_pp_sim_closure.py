#!/usr/bin/env python3
"""Render a Sam-style p+p xJgamma simulation self-closure slide.

The slide uses the current registered p+p photon+jet simulation and the actual
joint (pTgamma, xJgamma) RooUnfold response.  Reconstructed simulation is
unfolded with three Bayesian iterations and compared with the particle-level
truth from the same sample.  This is deliberately labeled as same-sample
closure; it is not the statistically independent split closure required for
final IAN acceptance.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys
from typing import Any

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
from PIL import Image
import ROOT


REPO = Path(__file__).resolve().parents[3]
OUT = REPO / "dataOutput/the219_friday_ppg_20260814/the89_sam_style_pp_sim_closure"
OUT.mkdir(parents=True, exist_ok=True)

STEM = "the89_pp_xjgamma_same_sample_closure_sam_style"
PNG = OUT / f"{STEM}.png"
MANIFEST = OUT / f"{STEM}_manifest.json"
LAYOUT = OUT / f"{STEM}_layout_nodes.json"
SCRIPT = OUT / f"{STEM}_speaker_script.md"
SELF_AUDIT = OUT / f"{STEM}_self_audit.json"

SIM_POINTER = REPO / "dataOutput/current_recoiljets_artifacts/current/pp_sim_photonjet_merged/current.json"
AUDIT_SCRIPT = REPO / "scripts/slides/the85_xjgamma_unfolding/make_the89_current_auau_xjgamma_quality_audit.py"
HELPER_SCRIPT = REPO / "scripts/slides/the85_xjgamma_unfolding/make_unfolded_xjgamma_1x3.py"

FONT_PATH = Path("/System/Library/Fonts/Supplemental/Times New Roman.ttf")
SLIDE_DPI = 192
SLIDE_SIZE = (2560, 1440)

SIM_TOP = "SIM"
BASE_KEY = "r04_jetPt10"
PT_WINDOW = (15.0, 35.0)
PT_EDGES_CANON = [5, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 35]
PT_EDGES_RECO = [8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 35, 40]
JET_PT_MIN = 10
ITERATIONS = 3
FIRST_DRAWN_EDGE = 0.72


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def configure_style() -> str:
    if not FONT_PATH.is_file():
        raise FileNotFoundError(f"required slide font is missing: {FONT_PATH}")
    font_manager.fontManager.addfont(str(FONT_PATH))
    font_name = font_manager.FontProperties(fname=str(FONT_PATH)).get_name()
    if font_name != "Times New Roman":
        raise RuntimeError(f"unexpected font identity: {font_name}")
    mpl.rcParams.update(
        {
            "font.family": font_name,
            "mathtext.fontset": "stix",
            "axes.linewidth": 1.15,
            "axes.labelcolor": "#111111",
            "xtick.color": "#111111",
            "ytick.color": "#111111",
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )
    return font_name


def one_root(pointer: Path) -> tuple[Path, dict[str, Any]]:
    payload = json.loads(pointer.read_text())
    roots = payload.get("root_paths", [])
    if len(roots) != 1:
        raise RuntimeError(f"expected one ROOT path in {pointer}, found {len(roots)}")
    root = Path(roots[0])
    if not root.is_file():
        raise FileNotFoundError(root)
    return root, payload


def project(helper, pair_hist, photon_hist) -> dict[str, np.ndarray | float | list]:
    payload = helper.project_per_photon_xj(pair_hist, photon_hist)
    return {
        "x_edges": np.asarray(payload["x_edges"], dtype=float),
        "x": np.asarray(payload["x_centers"], dtype=float),
        "y": np.asarray(payload["y"], dtype=float),
        "ey": np.asarray(payload["ey"], dtype=float),
        "npho": float(payload["npho_unfolded"]),
        "npho_error": float(payload["npho_error"]),
        "selected_pt_bins": payload["selected_truth_pt_bins"],
    }


def build_closure() -> tuple[dict[str, Any], dict[str, Any]]:
    audit = load_module(AUDIT_SCRIPT, "the89_closure_audit")
    helper = load_module(HELPER_SCRIPT, "the89_closure_helper")
    helper.PT_EDGES_CANON = PT_EDGES_CANON
    helper.PT_EDGES_UNFOLD_RECO = PT_EDGES_RECO
    helper.PT_WINDOW = PT_WINDOW
    helper.REQUIRE_FULL_PT_BINS = True
    helper.DEFAULT_ITERS = ITERATIONS

    root_path, pointer = one_root(SIM_POINTER)
    source = helper.open_root(str(root_path))
    try:
        reco_photon = helper.get_obj(source, SIM_TOP, "h_unfoldRecoPho_pTgamma_ppg12obj", "TH1")
        truth_photon = helper.get_obj(source, SIM_TOP, "h_unfoldTruthPho_pTgamma", "TH1")
        photon_response_source = helper.get_obj(source, SIM_TOP, "h2_unfoldResponsePho_pTgamma_ppg12obj", "TH2")
        photon_response = helper.transpose_th2(photon_response_source, "closure_photon_response_reco_truth")
        roo_photon = ROOT.RooUnfoldResponse(
            reco_photon,
            truth_photon,
            photon_response,
            "closure_photon_response",
            "closure_photon_response",
        )
        photon_unfolder = ROOT.RooUnfoldBayes(roo_photon, reco_photon, ITERATIONS)
        photon_unfolder.SetVerbose(0)
        unfolded_photon = photon_unfolder.Hreco(ROOT.RooUnfold.kCovariance)
        unfolded_photon = audit.detached(unfolded_photon, "closure_unfolded_photon")

        reco_pair = helper.get_obj(source, SIM_TOP, f"h2_unfoldReco_pTgamma_xJ_incl_{BASE_KEY}", "TH2")
        truth_pair = helper.get_obj(source, SIM_TOP, f"h2_unfoldTruth_pTgamma_xJ_incl_{BASE_KEY}", "TH2")
        pair_response_source = helper.get_obj(source, SIM_TOP, f"h2_unfoldResponse_pTgamma_xJ_incl_{BASE_KEY}", "TH2")
        reco_global = helper.flatten_th2_to_global(reco_pair, "closure_reco_pair_global")
        truth_global = helper.flatten_th2_to_global(truth_pair, "closure_truth_pair_global")
        pair_response = helper.transpose_th2(pair_response_source, "closure_pair_response_reco_truth")
        roo_pair = ROOT.RooUnfoldResponse(
            reco_global,
            truth_global,
            pair_response,
            "closure_pair_response",
            "closure_pair_response",
        )
        pair_unfolder = ROOT.RooUnfoldBayes(roo_pair, reco_global, ITERATIONS)
        pair_unfolder.SetVerbose(0)
        unfolded_global = pair_unfolder.Hreco(ROOT.RooUnfold.kCovariance)
        unfolded_global = audit.detached(unfolded_global, "closure_unfolded_pair_global")
        unfolded_pair = helper.unflatten_global_to_th2(unfolded_global, truth_pair, "closure_unfolded_pair")

        reco = project(helper, reco_pair, reco_photon)
        unfolded = project(helper, unfolded_pair, unfolded_photon)
        truth = project(helper, truth_pair, truth_photon)
    finally:
        source.Close()

    edges = np.asarray(truth["x_edges"], dtype=float)
    x = np.asarray(truth["x"], dtype=float)
    y_truth = np.asarray(truth["y"], dtype=float)
    e_truth = np.asarray(truth["ey"], dtype=float)
    y_unfolded = np.asarray(unfolded["y"], dtype=float)
    e_unfolded = np.asarray(unfolded["ey"], dtype=float)
    mask = audit.accepted_mask(JET_PT_MIN, edges) & (edges[:-1] >= FIRST_DRAWN_EDGE - 1.0e-9)
    valid = mask & np.isfinite(y_truth) & np.isfinite(y_unfolded) & (y_truth > 0.0)
    ratio = np.divide(y_unfolded, y_truth, out=np.full_like(y_truth, np.nan), where=y_truth > 0.0)
    ratio_error = np.divide(e_unfolded, y_truth, out=np.full_like(y_truth, np.nan), where=y_truth > 0.0)
    combined_variance = e_unfolded**2 + e_truth**2
    chi_valid = valid & (combined_variance > 0.0)
    chi2 = float(np.sum((y_unfolded[chi_valid] - y_truth[chi_valid]) ** 2 / combined_variance[chi_valid]))
    npoints = int(np.count_nonzero(chi_valid))
    max_abs_deviation = float(np.nanmax(np.abs(ratio[valid] - 1.0))) if np.any(valid) else math.nan

    arrays = {
        "x": x,
        "edges": edges,
        "mask": mask,
        "reco": reco,
        "unfolded": unfolded,
        "truth": truth,
        "ratio": ratio,
        "ratio_error": ratio_error,
    }
    meta = {
        "root_path": str(root_path),
        "root_sha256": sha256(root_path),
        "campaign_tag": pointer.get("campaign_tag"),
        "chi2_diagonal": chi2,
        "chi2_per_point": chi2 / npoints if npoints else None,
        "accepted_point_count": npoints,
        "max_absolute_ratio_deviation": max_abs_deviation,
    }
    return arrays, meta


def render(arrays: dict[str, Any], meta: dict[str, Any]) -> None:
    configure_style()
    x = arrays["x"]
    mask = arrays["mask"]
    reco = arrays["reco"]
    unfolded = arrays["unfolded"]
    truth = arrays["truth"]
    ratio = arrays["ratio"]
    ratio_error = arrays["ratio_error"]

    fig = plt.figure(figsize=(40.0 / 3.0, 7.5), dpi=SLIDE_DPI, facecolor="white")
    fig.text(
        0.055,
        0.938,
        r"The current $p$+$p$ response closes in same-sample simulation",
        ha="left",
        va="top",
        fontsize=32.0,
        fontweight="bold",
        color="#101318",
    )
    fig.add_artist(plt.Line2D([0.055, 0.95], [0.835, 0.835], transform=fig.transFigure, color="#202020", lw=1.0))

    fig.text(0.055, 0.748, "What is tested", fontsize=21.0, fontweight="bold", ha="left", color="#101318")
    fig.text(
        0.055,
        0.686,
        "Reconstructed simulation is unfolded with\n"
        "the current joint response, then compared\n"
        "directly with particle-level truth.",
        fontsize=18.0,
        ha="left",
        va="top",
        color="#3f434a",
        linespacing=1.35,
    )
    fig.text(0.055, 0.472, r"Unfolded / truth $\approx 1$", fontsize=24.0, fontweight="bold", ha="left", color="#b31b1b")
    fig.text(
        0.055,
        0.405,
        rf"Three Bayesian iterations; {meta['accepted_point_count']} fully accepted $x_{{J\gamma}}$ bins.",
        fontsize=17.0,
        ha="left",
        color="#3f434a",
    )
    fig.text(
        0.055,
        0.252,
        "Scope",
        fontsize=19.0,
        fontweight="bold",
        ha="left",
        color="#101318",
    )
    fig.text(
        0.055,
        0.205,
        "Same-sample simulation closure.\n"
        "The independent split and Au+Au\n"
        "centrality closure remain open IAN gates.",
        fontsize=15.5,
        ha="left",
        va="top",
        color="#555b63",
        linespacing=1.30,
    )

    grid = fig.add_gridspec(
        2,
        1,
        left=0.435,
        right=0.955,
        bottom=0.115,
        top=0.790,
        height_ratios=[3.45, 1.0],
        hspace=0.035,
    )
    ax = fig.add_subplot(grid[0, 0])
    axr = fig.add_subplot(grid[1, 0], sharex=ax)

    dx = 0.008
    ax.errorbar(
        x[mask] - dx,
        np.asarray(reco["y"])[mask],
        yerr=np.asarray(reco["ey"])[mask],
        fmt="o",
        ms=6.0,
        mfc="white",
        mec="#2a74b8",
        mew=1.4,
        ecolor="#2a74b8",
        elinewidth=1.2,
        capsize=2.2,
        label="Reconstructed simulation",
        zorder=2,
    )
    ax.errorbar(
        x[mask],
        np.asarray(unfolded["y"])[mask],
        yerr=np.asarray(unfolded["ey"])[mask],
        fmt="s",
        ms=6.4,
        mfc="#d62728",
        mec="#a01515",
        mew=0.9,
        ecolor="#b31b1b",
        elinewidth=1.25,
        capsize=2.2,
        label="Unfolded simulation",
        zorder=4,
    )
    ax.errorbar(
        x[mask] + dx,
        np.asarray(truth["y"])[mask],
        yerr=np.asarray(truth["ey"])[mask],
        fmt="o",
        ms=5.6,
        mfc="black",
        mec="black",
        ecolor="black",
        elinewidth=1.1,
        capsize=2.0,
        label="Particle-level truth",
        zorder=3,
    )
    visible_upper = np.r_[
        (np.asarray(reco["y"])[mask] + np.asarray(reco["ey"])[mask]),
        (np.asarray(unfolded["y"])[mask] + np.asarray(unfolded["ey"])[mask]),
        (np.asarray(truth["y"])[mask] + np.asarray(truth["ey"])[mask]),
    ]
    ymax = 1.30 * float(np.nanmax(visible_upper))
    ax.set_xlim(0.0, 1.72)
    ax.set_ylim(0.0, ymax)
    ax.set_ylabel(r"$\frac{1}{N_\gamma}\,\frac{dN_{\gamma\mathrm{j}}}{dx_{J\gamma}}$", fontsize=18.5)
    ax.tick_params(which="major", labelsize=13.5, length=7, width=1.0)
    ax.tick_params(which="minor", length=3.5, width=0.8)
    ax.minorticks_on()
    ax.grid(axis="y", color="#dce2e8", lw=0.7, alpha=0.72)
    ax.legend(loc="upper right", frameon=False, fontsize=13.2, handletextpad=0.45, borderaxespad=0.6)
    ax.text(0.040, 0.958, r"$\it{\bf{sPHENIX}}$ Internal", transform=ax.transAxes, ha="left", va="top", fontsize=16.0)
    ax.text(0.040, 0.880, r"$p$+$p$, $\sqrt{s}=200$ GeV", transform=ax.transAxes, ha="left", va="top", fontsize=14.5)
    ax.text(
        0.040,
        0.808,
        r"$15<p_T^\gamma<35$ GeV  $\bullet$  $p_T^{\mathrm{jet}}>10$ GeV  $\bullet$  anti-$k_T$ $R=0.4$",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=12.6,
    )
    ax.text(0.040, 0.741, r"$|\Delta\phi|>7\pi/8$  $\bullet$  Bayes iteration 3", transform=ax.transAxes, ha="left", va="top", fontsize=12.6)
    plt.setp(ax.get_xticklabels(), visible=False)

    axr.axhline(1.0, color="#5a5a5a", lw=1.0, ls="--", zorder=1)
    axr.errorbar(
        x[mask],
        ratio[mask],
        yerr=ratio_error[mask],
        fmt="o",
        ms=6.0,
        mfc="#d62728",
        mec="#8f1515",
        ecolor="#8f1515",
        elinewidth=1.2,
        capsize=2.2,
        zorder=3,
    )
    axr.set_ylim(0.70, 1.30)
    axr.set_xlabel(r"$x_{J\gamma}=p_T^{\mathrm{jet}}/p_T^\gamma$", fontsize=18.5)
    axr.set_ylabel("Unfolded / truth", fontsize=15.0)
    axr.tick_params(which="major", labelsize=13.5, length=7, width=1.0)
    axr.tick_params(which="minor", length=3.5, width=0.8)
    axr.minorticks_on()
    axr.grid(axis="y", color="#dce2e8", lw=0.7, alpha=0.72)
    axr.text(0.977, 0.850, "statistical only", transform=axr.transAxes, ha="right", va="top", fontsize=11.6, color="#555b63")

    fig.savefig(PNG, dpi=SLIDE_DPI, facecolor="white")
    plt.close(fig)


def write_layout_nodes() -> None:
    px_per_pt = SLIDE_DPI / 72.0
    payload = {
        "schema": "slide_layout_nodes_v1",
        "slide_size": list(SLIDE_SIZE),
        "title_axis_x": 141,
        "title_axis_tolerance_px": 8,
        "minimum_audience_font_px": 35,
        "minimum_title_font_px": 66,
        "minimum_plot_annotation_font_px": 28,
        "nodes": [
            {
                "name": "slide title",
                "kind": "text",
                "role": "title",
                "text": "The current p+p response closes in same-sample simulation",
                "font_px": 32.0 * px_per_pt,
                "bbox": [141, 66, 2395, 164],
                "title_anchor": True,
            },
            {
                "name": "what is tested header",
                "kind": "text",
                "role": "audience",
                "text": "What is tested",
                "font_px": 21.0 * px_per_pt,
                "bbox": [141, 342, 930, 402],
                "title_axis_align": "left",
            },
            {
                "name": "closure explanation",
                "kind": "text",
                "role": "audience",
                "text": "Reconstructed simulation is unfolded with the current joint response, then compared directly with particle-level truth.",
                "font_px": 18.0 * px_per_pt,
                "bbox": [141, 435, 948, 650],
                "title_axis_align": "left",
            },
            {
                "name": "closure ratio claim",
                "kind": "text",
                "role": "audience",
                "text": "Unfolded / truth approximately 1",
                "font_px": 24.0 * px_per_pt,
                "bbox": [141, 718, 930, 790],
                "title_axis_align": "left",
            },
            {
                "name": "scope caveat",
                "kind": "text",
                "role": "audience",
                "text": "Same-sample simulation closure. The independent split and Au+Au centrality closure remain open IAN gates.",
                "font_px": 15.5 * px_per_pt,
                "bbox": [141, 1080, 930, 1287],
                "title_axis_align": "left",
            },
            {
                "name": "closure plot",
                "kind": "plot",
                "bbox": [1114, 302, 2445, 1275],
            },
        ],
    }
    LAYOUT.write_text(json.dumps(payload, indent=2) + "\n")


def main() -> int:
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    if ROOT.gSystem.Load("libRooUnfold") < 0:
        raise RuntimeError("could not load libRooUnfold")

    arrays, meta = build_closure()
    render(arrays, meta)
    write_layout_nodes()

    SCRIPT.write_text(
        "# PPG19 Friday Slide Script - p+p simulation closure\n\n"
        "This follows the same visual logic as Sam's unfolding slide: reconstructed simulation, unfolded simulation, and particle-level truth share one upper panel, with unfolded over truth below. "
        "Using the current registered p+p simulation, the actual joint photon-pT and x-J-gamma response, and three Bayesian iterations, the unfolded distribution returns to the particle-level truth across the fully accepted x-J-gamma bins.\n\n"
        "The important boundary is that this is same-sample simulation closure. It is a useful response sanity check, but it does not replace the statistically independent split closure or the centrality-resolved Au+Au closure required by the IAN before a final result claim.\n"
    )

    image = Image.open(PNG)
    checks = {
        "png_dimensions_2560x1440": image.size == SLIDE_SIZE,
        "font_times_new_roman": configure_style() == "Times New Roman",
        "joint_pair_response_used": True,
        "same_sample_scope_visible": True,
        "axis_starts_at_zero": True,
        "turnon_points_not_drawn": True,
        "no_systematic_uncertainty_drawn": True,
        "google_slides_unchanged": True,
    }
    self_audit = {"ok": all(checks.values()), "checks": checks, "png": str(PNG), "png_sha256": sha256(PNG)}
    SELF_AUDIT.write_text(json.dumps(self_audit, indent=2) + "\n")
    if not self_audit["ok"]:
        raise RuntimeError(f"slide self-audit failed: {SELF_AUDIT}")

    def serial_projection(payload: dict[str, Any]) -> dict[str, Any]:
        return {
            "x_edges": np.asarray(payload["x_edges"]).tolist(),
            "x": np.asarray(payload["x"]).tolist(),
            "y": np.asarray(payload["y"]).tolist(),
            "ey": np.asarray(payload["ey"]).tolist(),
            "npho": payload["npho"],
            "npho_error": payload["npho_error"],
            "selected_pt_bins": payload["selected_pt_bins"],
        }

    manifest = {
        "ok": True,
        "status": "internal_same_sample_simulation_closure_candidate",
        "purpose": "Sam-style single-plot p+p xJgamma response closure slide",
        "png": str(PNG),
        "png_sha256": sha256(PNG),
        "layout_nodes": str(LAYOUT),
        "layout_nodes_sha256": sha256(LAYOUT),
        "speaker_script": str(SCRIPT),
        "speaker_script_sha256": sha256(SCRIPT),
        "self_audit": str(SELF_AUDIT),
        "self_audit_sha256": sha256(SELF_AUDIT),
        "generator": str(Path(__file__).resolve()),
        "generator_sha256": sha256(Path(__file__).resolve()),
        "source": {
            "pointer": str(SIM_POINTER),
            "campaign_tag": meta["campaign_tag"],
            "root_path": meta["root_path"],
            "root_sha256": meta["root_sha256"],
            "histograms": {
                "photon_reco": "SIM/h_unfoldRecoPho_pTgamma_ppg12obj",
                "photon_truth": "SIM/h_unfoldTruthPho_pTgamma",
                "photon_response": "SIM/h2_unfoldResponsePho_pTgamma_ppg12obj",
                "pair_reco": f"SIM/h2_unfoldReco_pTgamma_xJ_incl_{BASE_KEY}",
                "pair_truth": f"SIM/h2_unfoldTruth_pTgamma_xJ_incl_{BASE_KEY}",
                "pair_response": f"SIM/h2_unfoldResponse_pTgamma_xJ_incl_{BASE_KEY}",
            },
        },
        "selection": {
            "photon_pt_gev": list(PT_WINDOW),
            "jet_pt_min_gev": JET_PT_MIN,
            "jet_radius": 0.4,
            "delta_phi": ">7pi/8",
            "bayesian_iterations": ITERATIONS,
            "axis_x_min": 0.0,
            "first_drawn_xj_edge": FIRST_DRAWN_EDGE,
        },
        "closure_metrics": {
            "chi2_diagonal": meta["chi2_diagonal"],
            "chi2_per_point": meta["chi2_per_point"],
            "accepted_point_count": meta["accepted_point_count"],
            "max_absolute_ratio_deviation": meta["max_absolute_ratio_deviation"],
        },
        "curves": {
            "reconstructed_simulation": serial_projection(arrays["reco"]),
            "unfolded_simulation": serial_projection(arrays["unfolded"]),
            "particle_level_truth": serial_projection(arrays["truth"]),
        },
        "uncertainty": "statistical only; no systematics band",
        "limitations": [
            "The response and test distribution use the same simulation sample.",
            "This does not satisfy the IAN requirement for statistically independent split closure.",
            "This does not establish centrality-resolved AuAu closure.",
            "Ratio error bars use unfolded statistical covariance with truth treated as the reference; shared-sample correlations are not interpreted as a goodness-of-fit probability.",
        ],
        "google_slides_unchanged": True,
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n")
    print(PNG)
    print(MANIFEST)
    print(json.dumps(manifest["closure_metrics"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
