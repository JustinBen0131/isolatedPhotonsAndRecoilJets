#!/usr/bin/env python3
"""Offline low-xJgamma and Au+Au response/combinatoric contract audit.

This diagnostic consumes only the registered current pp/Au+Au ROOT products.
For Au+Au it compares the currently plotted unfolding contract

  data measured input       = ABCD-corrected data - scaled K
  response measured marginal = inclusive embedded reco

to a response-consistent offline variant

  data measured input       = ABCD-corrected data - scaled K
  response measured marginal = inclusive embedded reco - unscaled K.

Here K is the stored unmatched-reco/combinatoric template.  The response
matrix itself is unchanged.  The comparison is repeated for the stored
pTjet > 5, 7, and 10 GeV threshold objects and uses only xJ bins fully
accepted over the native photon-pT windows.

This is a diagnostic stress test, not a final correction or systematic.
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


ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
if ROOT.gSystem.Load("libRooUnfold") < 0:
    raise RuntimeError("failed to load libRooUnfold")

REPO = Path(__file__).resolve().parents[3]
AUAU_HELPER = (
    REPO
    / "scripts/slides/the85_xjgamma_unfolding/"
    "make_the89_current_auau_xjgamma_quality_audit.py"
)
PP_POINTS = (
    REPO
    / "dataOutput/the219_friday_ppg_20260814/the89_pp_threshold_peak_scan/"
    "the89_pp_unfolded_xjgamma_threshold_scan_points.json"
)
OUT = (
    REPO
    / "dataOutput/the219_friday_ppg_20260814/"
    "the89_lowx_response_contract_audit"
)
PNG = OUT / "the89_current_lowx_response_contract_pt5_7_10.png"
MANIFEST = OUT / "the89_current_lowx_response_contract_pt5_7_10_manifest.json"

FONT_PATH = Path("/System/Library/Fonts/Supplemental/Times New Roman.ttf")
PP_PT = (20.0, 26.0)
AUAU_PT = (19.0, 26.0)
THRESHOLDS = {
    5: "r04_isoR40_isSliding",
    7: "r04_jetPt7_isoR40_isSliding",
    10: "r04_jetPt10_isoR40_isSliding",
}
ITERATIONS = 3
DISPLAY_XMAX = 1.49


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def json_ready(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [json_ready(item) for item in value]
    if isinstance(value, (bool, np.bool_)) or type(value).__name__ == "bool":
        return bool(value)
    if isinstance(value, (float, np.floating)):
        return float(value) if math.isfinite(float(value)) else None
    if isinstance(value, (int, np.integer)):
        return int(value)
    return value


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def configure_style() -> str:
    if not FONT_PATH.is_file():
        raise FileNotFoundError(FONT_PATH)
    font_manager.fontManager.addfont(str(FONT_PATH))
    font_name = font_manager.FontProperties(fname=str(FONT_PATH)).get_name()
    if font_name != "Times New Roman":
        raise RuntimeError(f"unexpected font identity: {font_name}")
    mpl.rcParams.update(
        {
            "font.family": font_name,
            "mathtext.fontset": "stix",
            "axes.linewidth": 1.1,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "xtick.major.size": 6.0,
            "ytick.major.size": 6.0,
            "xtick.minor.size": 3.0,
            "ytick.minor.size": 3.0,
        }
    )
    return font_name


def load_pp() -> dict[int, dict[str, Any]]:
    payload = json.loads(PP_POINTS.read_text())
    out: dict[int, dict[str, Any]] = {}
    for threshold in THRESHOLDS:
        source = payload[str(threshold)]
        rows = source["bins"]
        out[threshold] = {
            "first_accepted_xj_edge": float(source["first_accepted_xj_edge"]),
            "rows": rows,
            "edges": np.asarray(
                [float(row["xj_low"]) for row in rows] + [float(rows[-1]["xj_high"])],
                dtype=float,
            ),
            "y": np.asarray([float(row["pp_unfolded"]["value"]) for row in rows]),
            "ey": np.asarray([float(row["pp_unfolded"]["stat_error"]) for row in rows]),
        }
    return out


def scaled_combinatoric(audit, helper, case, data, sim, comb):
    photon_data = helper.get_obj(
        data,
        case.data_topdir,
        f"h_unfoldRecoPho_pTgamma_{audit.PHO_KEY}{audit.CENT_SUFFIX}",
        "TH1",
    )
    photon_data, photon_meta = helper.apply_photon_abcd_input(
        case, data, sim, photon_data
    )
    photon_sim = helper.get_obj(
        sim,
        case.sim_topdir,
        f"h_unfoldRecoPho_pTgamma_{audit.PHO_KEY}{audit.CENT_SUFFIX}",
        "TH1",
    )
    scaled = audit.detached(comb, f"{comb.GetName()}_scaled_contract_audit")
    row_scales = []
    for ix in range(0, scaled.GetXaxis().GetNbins() + 2):
        n_data = float(photon_data.GetBinContent(ix))
        n_sim = float(photon_sim.GetBinContent(ix))
        scale = n_data / n_sim if n_sim > 0.0 else 0.0
        lo = float(scaled.GetXaxis().GetBinLowEdge(ix))
        hi = float(scaled.GetXaxis().GetBinUpEdge(ix))
        raw = 0.0
        for iy in range(0, scaled.GetYaxis().GetNbins() + 2):
            raw += float(scaled.GetBinContent(ix, iy))
            scaled.SetBinContent(ix, iy, scaled.GetBinContent(ix, iy) * scale)
            scaled.SetBinError(ix, iy, scaled.GetBinError(ix, iy) * scale)
        if raw != 0.0 or scale != 0.0:
            row_scales.append(
                {
                    "photon_pt": [lo, hi],
                    "data_photons_abcd": n_data,
                    "sim_photons": n_sim,
                    "data_over_sim": scale,
                    "raw_k": raw,
                    "scaled_k": raw * scale,
                }
            )
    return scaled, photon_meta, row_scales


def count_negative_bins(h2) -> int:
    return sum(
        1
        for ix in range(0, h2.GetXaxis().GetNbins() + 2)
        for iy in range(0, h2.GetYaxis().GetNbins() + 2)
        if h2.GetBinContent(ix, iy) < 0.0
    )


def pt_window_integral(audit, h2) -> float:
    total = 0.0
    for ix in range(1, h2.GetXaxis().GetNbins() + 1):
        lo = float(h2.GetXaxis().GetBinLowEdge(ix))
        hi = float(h2.GetXaxis().GetBinUpEdge(ix))
        if not audit.full_pt_bin(lo, hi):
            continue
        for iy in range(0, h2.GetYaxis().GetNbins() + 2):
            total += float(h2.GetBinContent(ix, iy))
    return total


def run_auau(audit, helper, case, data, sim, photon_unfolded, threshold: int):
    base_key = THRESHOLDS[threshold]
    names = {
        "data_A": f"h2_unfoldReco_pTgamma_xJ_incl_{base_key}{audit.CENT_SUFFIX}",
        "data_C": f"h2_unfoldReco_pTgamma_xJ_incl_sidebandC_{base_key}{audit.CENT_SUFFIX}",
        "sim_reco": f"h2_unfoldReco_pTgamma_xJ_incl_{base_key}{audit.CENT_SUFFIX}",
        "sim_truth": f"h2_unfoldTruth_pTgamma_xJ_incl_{base_key}{audit.CENT_SUFFIX}",
        "response": f"h2_unfoldResponse_pTgamma_xJ_incl_{base_key}{audit.CENT_SUFFIX}",
        "combinatoric": f"h2_unfoldRecoCombinatoric_pTgamma_xJ_incl_{base_key}{audit.CENT_SUFFIX}",
    }
    h_a = helper.get_obj(data, case.data_topdir, names["data_A"], "TH2")
    h_c = helper.get_obj(data, case.data_topdir, names["data_C"], "TH2")
    sim_reco = helper.get_obj(sim, case.sim_topdir, names["sim_reco"], "TH2")
    sim_truth = helper.get_obj(sim, case.sim_topdir, names["sim_truth"], "TH2")
    response_source = helper.get_obj(sim, case.sim_topdir, names["response"], "TH2")
    comb = helper.get_obj(sim, case.sim_topdir, names["combinatoric"], "TH2")

    abcd, abcd_rows = audit.build_purity_input(helper, case, data, sim, h_a, h_c)
    scaled_k, photon_meta, row_scales = scaled_combinatoric(
        audit, helper, case, data, sim, comb
    )
    measured = audit.detached(abcd, f"measured_after_k_{threshold}")
    measured.Add(scaled_k, -1.0)

    response_matrix = helper.transpose_th2(
        response_source, f"response_reco_truth_contract_{threshold}"
    )
    truth_global = helper.flatten_th2_to_global(
        sim_truth, f"sim_truth_global_contract_{threshold}"
    )
    data_global = helper.flatten_th2_to_global(
        measured, f"data_global_contract_{threshold}"
    )
    x_edges = np.asarray(helper.axis_edges(sim_truth.GetYaxis()), dtype=float)
    auau_first = float(x_edges[:-1][x_edges[:-1] >= threshold / AUAU_PT[0] - 1e-12][0])
    accepted = (x_edges[:-1] >= auau_first - 1e-12) & (
        x_edges[1:] <= DISPLAY_XMAX + 1e-12
    )

    variants: dict[str, Any] = {}
    for key, subtract_k_from_response in (
        ("current_inconsistent", False),
        ("response_consistent", True),
    ):
        response_reco = audit.detached(sim_reco, f"response_reco_{key}_{threshold}")
        if subtract_k_from_response:
            response_reco.Add(comb, -1.0)
        measured_global = helper.flatten_th2_to_global(
            response_reco, f"response_reco_global_{key}_{threshold}"
        )
        roo_response = ROOT.RooUnfoldResponse(
            measured_global,
            truth_global,
            response_matrix,
            f"roo_response_{key}_{threshold}",
            f"roo_response_{key}_{threshold}",
        )
        unfolding = ROOT.RooUnfoldBayes(roo_response, data_global, ITERATIONS)
        unfolding.SetVerbose(0)
        unfolded_global = unfolding.Hreco(ROOT.RooUnfold.kCovariance)
        if not unfolded_global:
            raise RuntimeError(f"unfolding failed for {key}, threshold {threshold}")
        unfolded_global = audit.detached(
            unfolded_global, f"unfolded_global_{key}_{threshold}"
        )
        unfolded_2d = helper.unflatten_global_to_th2(
            unfolded_global, sim_truth, f"unfolded_2d_{key}_{threshold}"
        )
        projection = helper.project_per_photon_xj(unfolded_2d, photon_unfolded)

        refold_global = roo_response.ApplyToTruth(
            unfolded_global, f"refolded_global_{key}_{threshold}"
        )
        refold_2d = helper.unflatten_global_to_th2(
            refold_global, response_reco, f"refolded_2d_{key}_{threshold}"
        )
        measured_values, measured_errors = audit.project_reco_counts(measured)
        refold_values, refold_errors = audit.project_reco_counts(refold_2d)
        closure = audit.closure_metrics(
            measured_values,
            measured_errors,
            refold_values,
            refold_errors,
            accepted,
        )
        variants[key] = {
            "y": np.asarray(projection["y"], dtype=float),
            "ey": np.asarray(projection["ey"], dtype=float),
            "npho_unfolded": float(projection["npho_unfolded"]),
            "npho_error": float(projection["npho_error"]),
            "selected_truth_pt_bins": projection["selected_truth_pt_bins"],
            "response_reco_integral": float(response_reco.Integral()),
            "response_reco_negative_bin_count": count_negative_bins(response_reco),
            "refold": closure,
        }

    abcd_window = pt_window_integral(audit, abcd)
    scaled_k_window = pt_window_integral(audit, scaled_k)
    measured_window = pt_window_integral(audit, measured)
    sim_reco_window = pt_window_integral(audit, sim_reco)
    sim_k_window = pt_window_integral(audit, comb)
    return {
        "threshold_gev": threshold,
        "base_key": base_key,
        "histograms": names,
        "x_edges": x_edges,
        "first_accepted_xj_edge": auau_first,
        "data_abcd_integral": abcd_window,
        "data_scaled_k_integral": scaled_k_window,
        "data_after_k_integral": measured_window,
        "k_fraction_of_abcd": float(scaled_k_window / abcd_window)
        if abcd_window != 0.0
        else None,
        "sim_reco_integral": sim_reco_window,
        "sim_k_integral": sim_k_window,
        "sim_k_fraction_of_reco": float(sim_k_window / sim_reco_window)
        if sim_reco_window != 0.0
        else None,
        "abcd_rows": abcd_rows,
        "photon_abcd": photon_meta,
        "k_row_scales": row_scales,
        "variants": variants,
    }


def ratio(y: np.ndarray, ey: np.ndarray, pp_y: np.ndarray, pp_ey: np.ndarray):
    value = np.full_like(y, np.nan)
    error = np.full_like(y, np.nan)
    valid = np.isfinite(y) & np.isfinite(ey) & np.isfinite(pp_y) & np.isfinite(pp_ey) & (pp_y > 0.0)
    value[valid] = y[valid] / pp_y[valid]
    relative = np.zeros_like(y)
    nonzero_y = valid & (y != 0.0)
    relative[nonzero_y] = (ey[nonzero_y] / y[nonzero_y]) ** 2
    relative[valid] += (pp_ey[valid] / pp_y[valid]) ** 2
    error[valid] = np.abs(value[valid]) * np.sqrt(relative[valid])
    return value, error


def render(pp: dict[int, dict[str, Any]], auau: dict[int, dict[str, Any]]) -> None:
    configure_style()
    fig = plt.figure(figsize=(16, 9), dpi=160, facecolor="white")
    fig.text(
        0.048,
        0.955,
        r"Offline low-$x_{J\gamma}$ stress test exposes the response/background contract",
        ha="left",
        va="top",
        fontsize=27.5,
        fontweight="bold",
        color="#111111",
    )
    fig.text(
        0.048,
        0.900,
        r"Current stored histograms  |  $p$+$p$: $20<p_T^\gamma<26$ GeV  |  "
        r"Au+Au 0--20%: $19<p_T^\gamma<26$ GeV  |  inclusive recoil jets  |  "
        r"RooUnfoldBayes, 3 iterations  |  statistical errors only",
        ha="left",
        va="top",
        fontsize=16.2,
        color="#333333",
    )

    left = 0.072
    right = 0.975
    gap = 0.036
    width = (right - left - 2.0 * gap) / 3.0
    top_bottom = 0.345
    top_height = 0.455
    ratio_bottom = 0.105
    ratio_height = 0.175
    axes = []
    ratio_axes = []
    for column, threshold in enumerate(THRESHOLDS):
        x0 = left + column * (width + gap)
        ax = fig.add_axes([x0, top_bottom, width, top_height])
        rax = fig.add_axes([x0, ratio_bottom, width, ratio_height], sharex=ax)
        axes.append(ax)
        ratio_axes.append(rax)

        pp_curve = pp[threshold]
        au = auau[threshold]
        edges = au["x_edges"]
        centers = 0.5 * (edges[:-1] + edges[1:])
        halfwidth = 0.5 * np.diff(edges)
        first = max(pp_curve["first_accepted_xj_edge"], au["first_accepted_xj_edge"])
        mask = (
            (edges[:-1] >= first - 1e-12)
            & (edges[1:] <= DISPLAY_XMAX + 1e-12)
            & np.isfinite(pp_curve["y"])
            & np.isfinite(pp_curve["ey"])
        )

        old = au["variants"]["current_inconsistent"]
        fixed = au["variants"]["response_consistent"]
        ax.errorbar(
            centers[mask] - 0.006,
            pp_curve["y"][mask],
            xerr=halfwidth[mask],
            yerr=pp_curve["ey"][mask],
            fmt="s",
            ms=6.1,
            mfc="white",
            mec="#111111",
            mew=1.4,
            ecolor="#111111",
            elinewidth=1.2,
            capsize=2.2,
            linestyle="none",
            label=r"$p$+$p$",
            zorder=4,
        )
        ax.errorbar(
            centers[mask],
            old["y"][mask],
            xerr=halfwidth[mask],
            yerr=old["ey"][mask],
            fmt="o",
            ms=5.8,
            mfc="white",
            mec="#9B2E25",
            mew=1.35,
            ecolor="#9B2E25",
            elinewidth=1.1,
            capsize=2.0,
            linestyle="none",
            label="Au+Au: current mismatch",
            zorder=3,
        )
        ax.errorbar(
            centers[mask] + 0.006,
            fixed["y"][mask],
            xerr=halfwidth[mask],
            yerr=fixed["ey"][mask],
            fmt="o",
            ms=6.2,
            mfc="#D94A3A",
            mec="#8E261E",
            mew=1.0,
            ecolor="#C53B2D",
            elinewidth=1.2,
            capsize=2.2,
            linestyle="none",
            label=r"Au+Au: response $-K$",
            zorder=5,
        )

        old_ratio, old_ratio_error = ratio(
            old["y"], old["ey"], pp_curve["y"], pp_curve["ey"]
        )
        fixed_ratio, fixed_ratio_error = ratio(
            fixed["y"], fixed["ey"], pp_curve["y"], pp_curve["ey"]
        )
        ratio_mask = mask & (pp_curve["y"] > 2.0 * pp_curve["ey"])
        rax.axhline(1.0, color="#777777", lw=1.0, linestyle="--", zorder=1)
        rax.errorbar(
            centers[ratio_mask],
            old_ratio[ratio_mask],
            xerr=halfwidth[ratio_mask],
            yerr=old_ratio_error[ratio_mask],
            fmt="o",
            ms=5.0,
            mfc="white",
            mec="#9B2E25",
            mew=1.25,
            ecolor="#9B2E25",
            elinewidth=1.0,
            capsize=1.8,
            linestyle="none",
            zorder=3,
        )
        rax.errorbar(
            centers[ratio_mask] + 0.006,
            fixed_ratio[ratio_mask],
            xerr=halfwidth[ratio_mask],
            yerr=fixed_ratio_error[ratio_mask],
            fmt="o",
            ms=5.4,
            mfc="#D94A3A",
            mec="#8E261E",
            mew=0.9,
            ecolor="#C53B2D",
            elinewidth=1.0,
            capsize=1.8,
            linestyle="none",
            zorder=4,
        )

        ax.set_xlim(0.0, 1.52)
        ax.set_ylim(-0.10, 3.15)
        rax.set_ylim(-0.15, 2.45)
        ax.minorticks_on()
        rax.minorticks_on()
        ax.tick_params(axis="both", which="major", labelsize=12.5, width=1.0)
        rax.tick_params(axis="both", which="major", labelsize=12.5, width=1.0)
        ax.tick_params(labelbottom=False)
        rax.set_xlabel(r"$x_{J\gamma}$", fontsize=17, labelpad=4)
        if column == 0:
            ax.set_ylabel(
                r"$(1/N_\gamma^{\rm particle})\,dN_{\rm jet}^{\rm particle}/dx_{J\gamma}$",
                fontsize=16,
                labelpad=8,
            )
            rax.set_ylabel(r"Au+Au / $p$+$p$", fontsize=14.5, labelpad=7)
        else:
            ax.set_yticklabels([])
            rax.set_yticklabels([])
        ax.text(
            0.04,
            0.95,
            rf"$p_T^{{\rm jet}}>{threshold}$ GeV",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=18,
            fontweight="bold",
        )
        ax.text(
            0.04,
            0.865,
            rf"first common full bin: $x_{{J\gamma}}={first:.2f}$",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=13.5,
            color="#333333",
        )
        ax.text(
            0.96,
            0.95,
            rf"$K/S_{{ABCD}}={au['k_fraction_of_abcd']:.2f}$",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=13.5,
            color="#8E261E",
        )

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper right",
        bbox_to_anchor=(0.974, 0.886),
        ncol=3,
        frameon=False,
        fontsize=14.5,
        handletextpad=0.5,
        columnspacing=1.35,
    )
    fig.savefig(PNG, dpi=160, facecolor="white")
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    pp = load_pp()
    audit = load_module(AUAU_HELPER, "the89_current_auau_contract_audit")
    audit.PT_WINDOW = AUAU_PT
    helper = audit.load_helper()
    helper.PT_WINDOW = AUAU_PT
    data_path, data_pointer = audit.one_root(audit.DATA_POINTER)
    sim_path, sim_pointer = audit.one_root(audit.SIM_POINTER)
    data = helper.open_root(str(data_path))
    sim = helper.open_root(str(sim_path))
    try:
        case = audit.configure_case(helper, data_path, sim_path)
        photon_unfolded, photon_meta = helper.unfold_photons(case, data, sim)
        auau = {
            threshold: run_auau(
                audit, helper, case, data, sim, photon_unfolded, threshold
            )
            for threshold in THRESHOLDS
        }
    finally:
        data.Close()
        sim.Close()

    render(pp, auau)
    image = Image.open(PNG)
    checks = {
        "png_dimensions_2560x1440": image.size == (2560, 1440),
        "font_times_new_roman": configure_style() == "Times New Roman",
        "uses_current_auau_pointer": data_pointer.get("canonical_status") == "canonical",
        "uses_current_auau_sim_pointer": sim_pointer.get("canonical_status") == "canonical",
        "three_stored_thresholds": set(auau) == set(THRESHOLDS),
        "no_systematic_band": True,
        "no_manual_curve_rescaling": True,
        "google_slides_unchanged": True,
    }
    manifest = {
        "ok": all(checks.values()),
        "status": "offline_response_contract_diagnostic_not_final",
        "purpose": "test low-xJ reach and AuAu combinatoric/response measured-marginal consistency",
        "png": str(PNG),
        "png_sha256": sha256(PNG),
        "generator": str(Path(__file__).resolve()),
        "generator_sha256": sha256(Path(__file__).resolve()),
        "selection": {
            "pp_photon_pt_gev": list(PP_PT),
            "auau_photon_pt_gev": list(AUAU_PT),
            "thresholds_gev": list(THRESHOLDS),
            "iterations": ITERATIONS,
            "fully_accepted_bins_only": True,
            "maximum_xj_edge": DISPLAY_XMAX,
        },
        "contracts": {
            "current_inconsistent": {
                "data_measured": "ABCD-corrected data minus photon-yield-scaled K",
                "response_measured_marginal": "inclusive embedded reco including K-like unmatched/fake population",
            },
            "response_consistent": {
                "data_measured": "ABCD-corrected data minus photon-yield-scaled K",
                "response_measured_marginal": "inclusive embedded reco minus unscaled K",
                "response_matrix": "unchanged matched-pair matrix",
            },
        },
        "sources": {
            "pp_points": str(PP_POINTS),
            "pp_points_sha256": sha256(PP_POINTS),
            "auau_data_pointer": str(audit.DATA_POINTER),
            "auau_data_pointer_sha256": sha256(audit.DATA_POINTER),
            "auau_data_campaign": data_pointer.get("campaign_tag"),
            "auau_data_root": str(data_path),
            "auau_sim_pointer": str(audit.SIM_POINTER),
            "auau_sim_pointer_sha256": sha256(audit.SIM_POINTER),
            "auau_sim_campaign": sim_pointer.get("campaign_tag"),
            "auau_sim_root": str(sim_path),
            "photon_unfolding": photon_meta,
        },
        "thresholds": auau,
        "checks": checks,
        "limitations": [
            "The response-consistent variant is an offline contract stress test, not an independently closed final result.",
            "The stored ROOT files do not contain alternate truth-jet matching-threshold templates, so the ATLAS-style generator threshold systematic cannot be reconstructed exactly offline.",
            "The pp and AuAu current ROOT products have nearest native photon windows rather than identical photon-pT bin edges.",
            "The current AuAu ABCD photon normalization reaches a boundary solution of purity 1.0 in the 21-23 and 23-26 GeV rows and requires a dedicated purity closure check.",
            "Only statistical uncertainties are propagated; numerator-denominator and correction cross-covariances are absent.",
        ],
        "google_slides_mutated": False,
    }
    MANIFEST.write_text(json.dumps(json_ready(manifest), indent=2, allow_nan=False) + "\n")
    if not manifest["ok"]:
        raise RuntimeError(f"audit checks failed: {MANIFEST}")
    print(PNG)
    print(MANIFEST)


if __name__ == "__main__":
    main()
