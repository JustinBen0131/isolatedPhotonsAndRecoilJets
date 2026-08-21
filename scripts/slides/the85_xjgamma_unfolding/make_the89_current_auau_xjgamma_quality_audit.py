#!/usr/bin/env python3
"""Audit and render the current Au+Au unfolded xJgamma candidate.

This is a local, no-production consumer of the canonical THE-88 Au+Au data and
embedded photon+jet ROOT files.  It deliberately follows the correction order
in AnalyzeRecoilJets_RooUnfoldPipeline.cpp:

  1. leakage-corrected event-leading ABCD photon normalization;
  2. region-C xJ background subtraction in each photon-pT row;
  3. photon-yield-scaled embedded combinatoric-jet subtraction;
  4. global-bin RooUnfoldBayes unfolding;
  5. refolding back to corrected reconstructed space.

The threshold/iteration scan is an audit.  It never rescales a physics curve
for appearance and it never draws a systematic band.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys
from typing import Any

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import ROOT


REPO = Path(__file__).resolve().parents[3]
HELPER_PATH = REPO / "scripts/slides/the85_xjgamma_unfolding/make_unfolded_xjgamma_1x3.py"
DATA_POINTER = REPO / "dataOutput/current_recoiljets_artifacts/current/auau_data_merged/current.json"
SIM_POINTER = REPO / "dataOutput/current_recoiljets_artifacts/current/auau_sim_photonjet_merged/current.json"
DEFAULT_OUTDIR = REPO / "dataOutput/the219_friday_ppg_20260814/the89_current_xjgamma_quality_audit"

DATA_TOP = "photon_10_plus_MBD_NS_geq_2_vtx_lt_150"
SIM_TOP = "SIM"
CENT_SUFFIX = "_cent_0_20"
PHO_KEY = "isoR40_isSliding"
PT_WINDOW = (15.0, 35.0)
CANONICAL_PT_EDGES = [15, 17, 19, 21, 23, 26, 35]
RECO_PT_EDGES = [10, 15, 17, 19, 21, 23, 26, 35, 40]
XJ_EDGES = np.array([0.0, 0.2, 0.24, 0.29, 0.35, 0.41, 0.5, 0.6, 0.72,
                     0.86, 1.03, 1.24, 1.49, 1.78, 2.14, 3.0], dtype=float)
THRESHOLDS = {
    5: "r04_isoR40_isSliding",
    7: "r04_jetPt7_isoR40_isSliding",
    10: "r04_jetPt10_isoR40_isSliding",
    12: "r04_jetPt12_isoR40_isSliding",
}
MAX_ITERATION = 6
SELECTED_THRESHOLD = 10


def load_helper():
    spec = importlib.util.spec_from_file_location("the85_unfold_helper", HELPER_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load helper {HELPER_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    module.PT_EDGES_CANON = CANONICAL_PT_EDGES
    module.PT_EDGES_UNFOLD_RECO = RECO_PT_EDGES
    module.PT_WINDOW = PT_WINDOW
    module.REQUIRE_FULL_PT_BINS = True
    module.PHO_KEY = PHO_KEY
    module.DEFAULT_ITERS = 3
    module.ERROR_MODE_NAME = "kCovariance"
    module.ERROR_MODE = ROOT.RooUnfold.kCovariance

    def current_suffix(lo: int, hi: int) -> str:
        if (lo, hi) == (10, 15):
            return "_pT_15_17"
        if (lo, hi) == (35, 40):
            return "_pT_26_35"
        return f"_pT_{lo}_{hi}"

    module.canonical_suffix_for_reco_bin = current_suffix
    return module


def one_root(pointer_path: Path) -> tuple[Path, dict[str, Any]]:
    payload = json.loads(pointer_path.read_text())
    paths = payload.get("root_paths", [])
    if len(paths) != 1:
        raise RuntimeError(f"expected one ROOT path in {pointer_path}, found {len(paths)}")
    root = Path(paths[0])
    if not root.is_file():
        raise FileNotFoundError(root)
    return root, payload


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def detached(obj, name: str):
    out = obj.Clone(name)
    if not out:
        raise RuntimeError(f"failed to clone {obj.GetName()}")
    if hasattr(out, "SetDirectory"):
        out.SetDirectory(0)
    if hasattr(out, "Sumw2"):
        out.Sumw2()
    return out


def full_pt_bin(lo: float, hi: float) -> bool:
    return lo >= PT_WINDOW[0] and hi <= PT_WINDOW[1]


def configure_case(helper, data_path: Path, sim_path: Path):
    return helper.Case(
        key="auau_current_0_20",
        title="Au+Au 0-20%",
        short_label="Au+Au 0-20%",
        data_file=str(data_path),
        data_topdir=DATA_TOP,
        sim_file=str(sim_path),
        sim_topdir=SIM_TOP,
        cent_suffix=CENT_SUFFIX,
        color="#d62728",
        marker="o",
        apply_abcd=True,
        apply_combinatoric_subtraction=True,
    )


def estimate_scale_variance(a: float, c_bkg: float, sa: float, esa: float, nbkg: float) -> float:
    if c_bkg <= 0.0 or nbkg <= 0.0:
        return 0.0
    e_nbkg2 = max(0.0, a) + esa * esa
    e_c2 = max(0.0, c_bkg)
    var = e_nbkg2 / (c_bkg * c_bkg) + (nbkg * nbkg * e_c2) / (c_bkg ** 4)
    return var if math.isfinite(var) and var > 0.0 else 0.0


def build_purity_input(helper, case, data_file, sim_file, h_a, h_c):
    """Current-bin equivalent of ApplyPurityCorrectionToRecoXJHist."""
    if h_c is None:
        raise KeyError("missing sideband-C xJ histogram")
    out = detached(h_a, f"{case.key}_purity_input")
    out.Reset("ICES")
    a_orig = detached(h_a, f"{case.key}_a_orig")
    reco_bins = helper.pt_bins_from_edges(RECO_PT_EDGES)
    leakage = helper.load_leakage_factors(sim_file, case.sim_topdir, case.cent_suffix)
    ny = out.GetYaxis().GetNbins()
    rows: list[dict[str, Any]] = []

    for lo, hi, _ in reco_bins:
        canon = helper.canonical_suffix_for_reco_bin(lo, hi)
        leak = leakage.get(canon, (0.0, 0.0, 0.0))
        suffix = helper.abcd_suffix(canon, case.cent_suffix)
        a, b, c, d, sa, esa, nbkg, meta = helper.compute_abcd_counts(
            data_file, case.data_topdir, suffix, leak
        )
        if a + b + c + d <= 0.0:
            continue
        f_c = max(0.0, leak[1])
        c_bkg = max(0.0, c - f_c * sa)
        scale_c = nbkg / c_bkg if c_bkg > 0.0 else 0.0
        denom = 1.0 - scale_c * f_c
        inv_denom = 1.0 / denom if math.isfinite(denom) and abs(denom) > 1e-6 else 1.0
        var_scale = estimate_scale_variance(a, c_bkg, sa, esa, nbkg)
        purity = sa / a if a > 0.0 else 0.0
        cen = 0.5 * (lo + hi)
        ix = out.GetXaxis().FindBin(cen)
        pooled = helper.pooled_sideband_c(h_c, ix, reco_bins, lo, hi, purity)
        for iy in range(0, ny + 2):
            value_a = a_orig.GetBinContent(ix, iy)
            error_a = a_orig.GetBinError(ix, iy)
            if pooled is None:
                value_c = h_c.GetBinContent(ix, iy)
                error_c = h_c.GetBinError(ix, iy)
            else:
                value_c = float(pooled[0][iy])
                error_c = float(pooled[1][iy])
            value = (value_a - scale_c * value_c) * inv_denom
            error2 = (error_a * error_a + scale_c * scale_c * error_c * error_c
                      + value_c * value_c * var_scale) * inv_denom * inv_denom
            out.SetBinContent(ix, iy, value if math.isfinite(value) else 0.0)
            out.SetBinError(ix, iy, math.sqrt(error2) if error2 > 0.0 and math.isfinite(error2) else 0.0)
        rows.append({
            "pt": [lo, hi], "A": a, "B": b, "C": c, "D": d, "SA": sa,
            "eSA": esa, "purity": purity, "scaleC": scale_c,
            "fCLeak": f_c, "pooled_C_shape": pooled is not None, **meta,
        })
    if not rows:
        raise RuntimeError("no current event-leading ABCD rows were found")
    return out, rows


def subtract_combinatoric(helper, case, data_file, sim_file, corrected, comb):
    photon_data = helper.get_obj(
        data_file, case.data_topdir,
        f"h_unfoldRecoPho_pTgamma_{PHO_KEY}{CENT_SUFFIX}", "TH1"
    )
    photon_data, photon_meta = helper.apply_photon_abcd_input(
        case, data_file, sim_file, photon_data
    )
    photon_sim = helper.get_obj(
        sim_file, case.sim_topdir,
        f"h_unfoldRecoPho_pTgamma_{PHO_KEY}{CENT_SUFFIX}", "TH1"
    )
    scaled = detached(comb, f"{comb.GetName()}_scaled")
    raw_integral = float(scaled.Integral())
    scaled_integral = 0.0
    for ix in range(0, scaled.GetXaxis().GetNbins() + 2):
        n_data = photon_data.GetBinContent(ix)
        n_sim = photon_sim.GetBinContent(ix)
        scale = n_data / n_sim if n_sim > 0.0 else 0.0
        row = 0.0
        for iy in range(0, scaled.GetYaxis().GetNbins() + 2):
            row += scaled.GetBinContent(ix, iy)
            scaled.SetBinContent(ix, iy, scaled.GetBinContent(ix, iy) * scale)
            scaled.SetBinError(ix, iy, scaled.GetBinError(ix, iy) * scale)
        scaled_integral += row * scale
    before = detached(corrected, f"{corrected.GetName()}_before_k")
    after = detached(corrected, f"{corrected.GetName()}_after_k")
    after.Add(scaled, -1.0)
    return after, {
        "raw_integral": raw_integral,
        "scaled_integral": scaled_integral,
        "purity_corrected_integral_before": float(before.Integral()),
        "integral_after": float(after.Integral()),
        "fraction_of_purity_corrected_input": (
            scaled_integral / float(before.Integral()) if before.Integral() != 0.0 else None
        ),
        "photon_abcd": photon_meta,
    }


def project_reco_counts(h2) -> tuple[np.ndarray, np.ndarray]:
    ny = h2.GetYaxis().GetNbins()
    values = np.zeros(ny)
    errors2 = np.zeros(ny)
    for ix in range(1, h2.GetXaxis().GetNbins() + 1):
        lo = h2.GetXaxis().GetBinLowEdge(ix)
        hi = h2.GetXaxis().GetBinUpEdge(ix)
        if not full_pt_bin(lo, hi):
            continue
        for iy in range(1, ny + 1):
            values[iy - 1] += h2.GetBinContent(ix, iy)
            errors2[iy - 1] += h2.GetBinError(ix, iy) ** 2
    return values, np.sqrt(errors2)


def accepted_mask(threshold: int, x_edges: np.ndarray) -> np.ndarray:
    uniform_turnon = threshold / PT_WINDOW[0]
    lows = x_edges[:-1]
    highs = x_edges[1:]
    return (lows + 1e-12 >= uniform_turnon) & (highs <= 1.49 + 1e-12)


def closure_metrics(measured: np.ndarray, measured_err: np.ndarray,
                    refold: np.ndarray, refold_err: np.ndarray, mask: np.ndarray) -> dict[str, Any]:
    variance = measured_err * measured_err + refold_err * refold_err
    valid = mask & np.isfinite(measured) & np.isfinite(refold) & (variance > 0.0)
    if not np.any(valid):
        return {"chi2": None, "ndf": 0, "chi2_ndf": None, "ratio_rms": None}
    chi2 = float(np.sum((refold[valid] - measured[valid]) ** 2 / variance[valid]))
    ratio_valid = valid & (measured != 0.0)
    ratios = refold[ratio_valid] / measured[ratio_valid]
    rms = float(np.sqrt(np.mean((ratios - 1.0) ** 2))) if ratios.size else None
    return {"chi2": chi2, "ndf": int(np.sum(valid)),
            "chi2_ndf": chi2 / int(np.sum(valid)), "ratio_rms": rms}


def run_threshold(helper, case, data_file, sim_file, photon_unfolded,
                  threshold: int, base_key: str) -> dict[str, Any]:
    helper.BASE_KEY = base_key
    names = {
        "data_A": f"h2_unfoldReco_pTgamma_xJ_incl_{base_key}{CENT_SUFFIX}",
        "data_C": f"h2_unfoldReco_pTgamma_xJ_incl_sidebandC_{base_key}{CENT_SUFFIX}",
        "sim_reco": f"h2_unfoldReco_pTgamma_xJ_incl_{base_key}{CENT_SUFFIX}",
        "sim_truth": f"h2_unfoldTruth_pTgamma_xJ_incl_{base_key}{CENT_SUFFIX}",
        "response": f"h2_unfoldResponse_pTgamma_xJ_incl_{base_key}{CENT_SUFFIX}",
        "combinatoric": f"h2_unfoldRecoCombinatoric_pTgamma_xJ_incl_{base_key}{CENT_SUFFIX}",
    }
    h_a = helper.get_obj(data_file, case.data_topdir, names["data_A"], "TH2")
    h_c = helper.get_obj(data_file, case.data_topdir, names["data_C"], "TH2")
    sim_reco = helper.get_obj(sim_file, case.sim_topdir, names["sim_reco"], "TH2")
    sim_truth = helper.get_obj(sim_file, case.sim_topdir, names["sim_truth"], "TH2")
    response_source = helper.get_obj(sim_file, case.sim_topdir, names["response"], "TH2")
    comb = helper.get_obj(sim_file, case.sim_topdir, names["combinatoric"], "TH2")
    corrected, abcd_rows = build_purity_input(helper, case, data_file, sim_file, h_a, h_c)
    corrected, comb_meta = subtract_combinatoric(
        helper, case, data_file, sim_file, corrected, comb
    )

    sim_reco_global = helper.flatten_th2_to_global(sim_reco, f"sim_reco_g_{threshold}")
    sim_truth_global = helper.flatten_th2_to_global(sim_truth, f"sim_truth_g_{threshold}")
    data_global = helper.flatten_th2_to_global(corrected, f"data_g_{threshold}")
    response = helper.transpose_th2(response_source, f"response_reco_truth_{threshold}")
    roo_response = ROOT.RooUnfoldResponse(
        sim_reco_global, sim_truth_global, response,
        f"response_threshold_{threshold}", f"response_threshold_{threshold}"
    )
    mask = accepted_mask(threshold, helper.axis_edges(sim_truth.GetYaxis()))
    iterations: dict[str, Any] = {}
    previous_values = None
    selected_iteration = 3
    selected_reason = "fallback to prior Au+Au and p+p stability choice"

    for iteration in range(1, MAX_ITERATION + 1):
        unfolding = ROOT.RooUnfoldBayes(roo_response, data_global, iteration)
        unfolding.SetVerbose(0)
        unfolded_global = unfolding.Hreco(ROOT.RooUnfold.kCovariance)
        if not unfolded_global:
            raise RuntimeError(f"unfolding failed for pT>{threshold}, iteration {iteration}")
        unfolded_global = detached(unfolded_global, f"unfolded_g_{threshold}_{iteration}")
        unfolded_2d = helper.unflatten_global_to_th2(
            unfolded_global, sim_truth, f"unfolded_2d_{threshold}_{iteration}"
        )
        projection = helper.project_per_photon_xj(unfolded_2d, photon_unfolded)
        values = np.asarray(projection["y"], dtype=float)
        errors = np.asarray(projection["ey"], dtype=float)
        accepted_values = values[mask]
        accepted_errors = errors[mask]
        rel_stat = (
            float(np.sqrt(np.sum(accepted_errors ** 2) / np.sum(accepted_values ** 2)))
            if np.sum(accepted_values ** 2) > 0.0 else None
        )
        rel_change = None
        if previous_values is not None and np.sum(accepted_values ** 2) > 0.0:
            rel_change = float(np.sqrt(
                np.sum((accepted_values - previous_values) ** 2) / np.sum(accepted_values ** 2)
            ))

        refold_global = roo_response.ApplyToTruth(
            unfolded_global, f"refolded_g_{threshold}_{iteration}"
        )
        if not refold_global:
            raise RuntimeError(f"refolding failed for pT>{threshold}, iteration {iteration}")
        refold_global = detached(refold_global, f"refolded_g_clone_{threshold}_{iteration}")
        refold_2d = helper.unflatten_global_to_th2(
            refold_global, sim_reco, f"refolded_2d_{threshold}_{iteration}"
        )
        measured, measured_err = project_reco_counts(corrected)
        refold, refold_err = project_reco_counts(refold_2d)
        closure = closure_metrics(measured, measured_err, refold, refold_err, mask)
        iterations[str(iteration)] = {
            "projection": {
                "x_edges": np.asarray(projection["x_edges"]).tolist(),
                "x": np.asarray(projection["x_centers"]).tolist(),
                "y": values.tolist(),
                "ey": errors.tolist(),
                "npho_unfolded": projection["npho_unfolded"],
                "npho_error": projection["npho_error"],
                "selected_truth_pt_bins": projection["selected_truth_pt_bins"],
            },
            "refold": {
                "measured": measured.tolist(), "measured_error": measured_err.tolist(),
                "prediction": refold.tolist(), "prediction_error": refold_err.tolist(),
            },
            "metrics": {
                **closure,
                "relative_stat": rel_stat,
                "relative_change_from_previous_iteration": rel_change,
                "accepted_bin_count": int(np.sum(mask)),
                "accepted_negative_bin_count": int(np.sum(accepted_values < 0.0)),
            },
        }
        if iteration >= 3 and rel_change is not None and rel_stat is not None:
            if rel_change <= max(rel_stat, 0.025) and selected_reason.startswith("fallback"):
                selected_iteration = iteration
                selected_reason = (
                    f"earliest stable iteration: relative change {rel_change:.4g} <= "
                    f"max(relative stat {rel_stat:.4g}, 0.025)"
                )
        previous_values = accepted_values.copy()

    return {
        "threshold_gev": threshold,
        "base_key": base_key,
        "histograms": names,
        "uniform_turnon": threshold / PT_WINDOW[0],
        "first_drawn_xj_edge": float(helper.axis_edges(sim_truth.GetYaxis())[:-1][mask][0]),
        "abcd_rows": abcd_rows,
        "combinatoric": comb_meta,
        "selected_iteration": selected_iteration,
        "selected_iteration_reason": selected_reason,
        "iterations": iterations,
    }


def write_scan_csv(results: dict[str, Any], path: Path) -> None:
    rows = []
    for threshold, result in results.items():
        for iteration, payload in result["iterations"].items():
            rows.append({
                "threshold_gev": threshold,
                "iteration": iteration,
                "selected": int(int(iteration) == result["selected_iteration"]),
                "uniform_turnon": result["uniform_turnon"],
                "first_drawn_xj_edge": result["first_drawn_xj_edge"],
                **payload["metrics"],
                "combinatoric_fraction": result["combinatoric"]["fraction_of_purity_corrected_input"],
            })
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def style() -> None:
    mpl.rcParams.update({
        "font.family": "DejaVu Sans",
        "mathtext.fontset": "dejavusans",
        "axes.linewidth": 1.25,
        "axes.labelsize": 16,
        "xtick.labelsize": 12.5,
        "ytick.labelsize": 12.5,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "legend.frameon": False,
    })


def draw_audit(results: dict[str, Any], path: Path) -> None:
    style()
    colors = {5: "#111111", 7: "#2b6cb0", 10: "#d62728", 12: "#6b46c1"}
    fig, axes = plt.subplots(2, 2, figsize=(13.6, 8.0), dpi=190, sharex=True)
    for ax, threshold in zip(axes.flat, [5, 7, 10, 12]):
        result = results[str(threshold)]
        selected = result["selected_iteration"]
        for iteration in [2, selected, 4]:
            payload = result["iterations"][str(iteration)]
            projection = payload["projection"]
            x = np.asarray(projection["x"])
            y = np.asarray(projection["y"])
            ey = np.asarray(projection["ey"])
            edges = np.asarray(projection["x_edges"])
            mask = accepted_mask(threshold, edges)
            alpha = 1.0 if iteration == selected else 0.34
            ax.errorbar(
                x[mask], y[mask], yerr=ey[mask], fmt="o-" if iteration == selected else "o--",
                ms=4.6 if iteration == selected else 3.8, lw=1.6 if iteration == selected else 1.0,
                color=colors[threshold], alpha=alpha,
                label=f"iteration {iteration}" + (" (selected)" if iteration == selected else ""),
            )
        metrics = result["iterations"][str(selected)]["metrics"]
        ax.set_xlim(0.0, 1.5)
        ax.set_ylim(bottom=-0.08)
        ax.minorticks_on()
        ax.set_title(rf"$p_T^{{\mathrm{{jet}}}}>{threshold}$ GeV", fontsize=15.5, fontweight="bold")
        ax.text(
            0.04, 0.94,
            rf"drawn $x_{{J\gamma}}\geq {result['first_drawn_xj_edge']:.2f}$" + "\n"
            + rf"refold $\chi^2/N={metrics['chi2_ndf']:.2f}$; stat={metrics['relative_stat']:.2f}",
            transform=ax.transAxes, va="top", fontsize=10.4,
        )
        ax.legend(loc="upper right", fontsize=9.2)
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$\frac{1}{N_\gamma}\,\frac{dN_{\mathrm{jet}}}{dx_{J\gamma}}$")
    for ax in axes[1, :]:
        ax.set_xlabel(r"$x_{J\gamma}=p_T^{\mathrm{jet}}/p_T^\gamma$")
    fig.suptitle("Current Au+Au unfolded $x_{J\gamma}$: threshold and iteration audit",
                 fontsize=20.5, fontweight="bold", y=0.982)
    fig.text(0.5, 0.912,
             r"0-20%, $15<p_T^\gamma<35$ GeV, anti-$k_T$ $R=0.4$, $\Delta\phi>7\pi/8$ • statistical uncertainties only",
             ha="center", fontsize=13.5)
    fig.subplots_adjust(left=0.09, right=0.985, bottom=0.10, top=0.845, wspace=0.17, hspace=0.20)
    fig.savefig(path, facecolor="white")
    plt.close(fig)


def draw_candidate(result: dict[str, Any], path: Path) -> None:
    style()
    iteration = result["selected_iteration"]
    payload = result["iterations"][str(iteration)]
    projection = payload["projection"]
    x = np.asarray(projection["x"])
    y = np.asarray(projection["y"])
    ey = np.asarray(projection["ey"])
    edges = np.asarray(projection["x_edges"])
    mask = accepted_mask(SELECTED_THRESHOLD, edges)
    measured = np.asarray(payload["refold"]["measured"])
    measured_err = np.asarray(payload["refold"]["measured_error"])
    predicted = np.asarray(payload["refold"]["prediction"])
    predicted_err = np.asarray(payload["refold"]["prediction_error"])
    ratio = np.divide(predicted, measured, out=np.full_like(predicted, np.nan), where=measured != 0.0)
    ratio_err = np.full_like(ratio, np.nan)
    valid = measured != 0.0
    ratio_err[valid] = np.sqrt(
        (predicted_err[valid] / measured[valid]) ** 2
        + (predicted[valid] * measured_err[valid] / (measured[valid] ** 2)) ** 2
    )
    metrics = payload["metrics"]

    fig = plt.figure(figsize=(40.0 / 3.0, 7.5), dpi=192)
    grid = fig.add_gridspec(1, 2, left=0.075, right=0.965, bottom=0.14, top=0.78, wspace=0.25)
    ax = fig.add_subplot(grid[0, 0])
    ax_refold = fig.add_subplot(grid[0, 1])

    ax.errorbar(x[mask], y[mask], yerr=ey[mask], fmt="o", ms=7.0, mfc="#d62728",
                mec="#8b1a1a", mew=0.9, ecolor="#8b1a1a", elinewidth=1.6,
                capsize=2.8, color="#d62728", label="Au+Au 0-20%")
    ax.set_xlim(0.0, 1.5)
    ymax = max(0.65, float(np.nanmax(y[mask] + ey[mask])) * 1.28)
    ax.set_ylim(0.0, ymax)
    ax.set_xlabel(r"$x_{J\gamma}=p_T^{\mathrm{jet}}/p_T^\gamma$")
    ax.set_ylabel(r"$\frac{1}{N_\gamma}\,\frac{dN_{\mathrm{jet}}}{dx_{J\gamma}}$")
    ax.minorticks_on()
    ax.legend(loc="upper right", fontsize=13.0, handletextpad=0.35)
    ax.text(0.05, 0.95, r"$\it{\bf{sPHENIX}}$ Internal", transform=ax.transAxes,
            ha="left", va="top", fontsize=15.0)
    ax.text(0.05, 0.875, r"Au+Au $\sqrt{s_{NN}}=200$ GeV", transform=ax.transAxes,
            ha="left", va="top", fontsize=13.0)
    ax.text(0.05, 0.81, rf"Bayes iteration {iteration} • stat. only", transform=ax.transAxes,
            ha="left", va="top", fontsize=12.0)

    ax_refold.axhline(1.0, color="0.25", lw=1.2, ls="--")
    ax_refold.errorbar(x[mask], ratio[mask], yerr=ratio_err[mask], fmt="o", ms=7.0,
                       mfc="#1f4e79", mec="#12324b", mew=0.9, ecolor="#12324b",
                       elinewidth=1.5, capsize=2.8)
    ax_refold.set_xlim(0.0, 1.5)
    finite = np.isfinite(ratio[mask]) & np.isfinite(ratio_err[mask])
    if np.any(finite):
        low = float(np.nanmin(ratio[mask][finite] - ratio_err[mask][finite]))
        high = float(np.nanmax(ratio[mask][finite] + ratio_err[mask][finite]))
        ax_refold.set_ylim(max(0.0, low - 0.18), min(2.1, max(1.35, high + 0.18)))
    else:
        ax_refold.set_ylim(0.4, 1.6)
    ax_refold.set_xlabel(r"reconstructed $x_{J\gamma}$")
    ax_refold.set_ylabel("refolded / corrected data")
    ax_refold.minorticks_on()
    ax_refold.text(0.05, 0.95, "Refolding validation", transform=ax_refold.transAxes,
                   va="top", fontsize=15.0, fontweight="bold")
    ax_refold.text(0.05, 0.865,
                   rf"$\chi^2/N={metrics['chi2_ndf']:.2f}$ (diagonal stat. diagnostic)",
                   transform=ax_refold.transAxes, va="top", fontsize=11.8)

    fig.suptitle(r"Unfolded $x_{J\gamma}$ with a reconstructed-space check",
                 x=0.075, y=0.955, ha="left", fontsize=23.0, fontweight="bold")
    fig.text(0.075, 0.865,
             r"$15<p_T^\gamma<35$ GeV • $p_T^{\mathrm{jet}}>10$ GeV • anti-$k_T$ $R=0.4$ • $\Delta\phi>7\pi/8$",
             ha="left", fontsize=15.5)
    fig.text(0.075, 0.815,
             rf"Axis begins at zero; points below the fully accepted region ($x_{{J\gamma}}<{result['first_drawn_xj_edge']:.2f}$) are not drawn.",
             ha="left", fontsize=12.6, color="#4b5563")
    fig.savefig(path, facecolor="white")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUTDIR)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    if ROOT.gSystem.Load("libRooUnfold") < 0:
        raise RuntimeError("could not load libRooUnfold")
    helper = load_helper()
    data_root, data_pointer = one_root(DATA_POINTER)
    sim_root, sim_pointer = one_root(SIM_POINTER)
    case = configure_case(helper, data_root, sim_root)
    data_file = helper.open_root(str(data_root))
    sim_file = helper.open_root(str(sim_root))

    helper.DEFAULT_ITERS = 3
    photon_unfolded, photon_meta = helper.unfold_photons(case, data_file, sim_file)
    results: dict[str, Any] = {}
    for threshold, base_key in THRESHOLDS.items():
        results[str(threshold)] = run_threshold(
            helper, case, data_file, sim_file, photon_unfolded, threshold, base_key
        )
    data_file.Close()
    sim_file.Close()

    scan_csv = args.out_dir / "the89_current_auau_xjgamma_threshold_iteration_scan.csv"
    audit_png = args.out_dir / "the89_current_auau_xjgamma_threshold_iteration_audit.png"
    candidate_png = args.out_dir / "the89_current_auau_xjgamma_pt10_unfold_refold_candidate.png"
    manifest_path = args.out_dir / "the89_current_auau_xjgamma_quality_audit_manifest.json"
    write_scan_csv(results, scan_csv)
    draw_audit(results, audit_png)
    draw_candidate(results[str(SELECTED_THRESHOLD)], candidate_png)

    manifest = {
        "status": "statistical_correction_candidate_not_final",
        "purpose": "current-input Au+Au xJgamma threshold/iteration audit and publication-style pTjet>10 refolding candidate",
        "script": str(Path(__file__).resolve()),
        "sources": {
            "data_pointer": str(DATA_POINTER),
            "data_campaign": data_pointer.get("campaign_tag"),
            "data_root": str(data_root),
            "data_root_sha256": sha256(data_root),
            "sim_pointer": str(SIM_POINTER),
            "sim_campaign": sim_pointer.get("campaign_tag"),
            "sim_root": str(sim_root),
            "sim_root_sha256": sha256(sim_root),
        },
        "selection": {
            "centrality": "0-20%",
            "photon_pt_gev": list(PT_WINDOW),
            "jet_radius": 0.4,
            "delta_phi": ">7pi/8",
            "photon_definition": PHO_KEY,
            "thresholds_audited_gev": list(THRESHOLDS),
            "display_threshold_gev": SELECTED_THRESHOLD,
            "display_points": "only stored xJ bins fully above pTjet_min / min(pTgamma)",
        },
        "method": {
            "correction_order": [
                "event-leading leakage-corrected ABCD photon normalization",
                "region-C xJ sideband subtraction with signed bins and scale variance",
                "photon-yield-scaled embedded combinatoric subtraction from measured reco input",
                "global-bin RooUnfoldBayes with kCovariance statistical errors",
                "ApplyToTruth refolding through the same RooUnfoldResponse",
            ],
            "iteration_rule": "earliest iteration >=3 whose change from the previous iteration is no larger than max(relative statistical uncertainty, 0.025)",
            "curve_rescaling": "none",
            "systematic_uncertainty": "not drawn; no systematic band is available",
            "closure_chi2_note": "diagonal statistical diagnostic; measured and refolded spectra are correlated",
        },
        "photon_unfolding": photon_meta,
        "results": results,
        "outputs": {
            "scan_csv": str(scan_csv),
            "audit_png": str(audit_png),
            "candidate_png": str(candidate_png),
        },
        "limitations": [
            "The Au+Au xJ-dependent background-shape and combinatoric-template closures are not independently final-approved.",
            "Only statistical covariance from the nominal unfolding is shown; systematic uncertainties are absent.",
            "The refolding chi2 is a correlated diagnostic and is not a formal goodness-of-fit probability.",
            "This local artifact does not constitute a production rerun, scientific acceptance, or publication approval.",
        ],
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(audit_png)
    print(candidate_png)
    print(scan_csv)
    print(manifest_path)
    for threshold in THRESHOLDS:
        result = results[str(threshold)]
        selected = result["selected_iteration"]
        print(json.dumps({
            "threshold": threshold,
            "selected_iteration": selected,
            "reason": result["selected_iteration_reason"],
            "first_drawn_xj_edge": result["first_drawn_xj_edge"],
            "combinatoric_fraction": result["combinatoric"]["fraction_of_purity_corrected_input"],
            "metrics": result["iterations"][str(selected)]["metrics"],
        }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
