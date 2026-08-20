#!/usr/bin/env python3
"""Build THE-219 slide 1: current-pointer p+p particle-level per-photon xJgamma.

This is a bounded, read-only consumer of the registered current p+p data and
photon+jet simulation ROOT files.  It mirrors the purity-corrected RooUnfold
path in AnalyzeRecoilJets_RooUnfoldPipeline.cpp and writes a full-slide PNG,
numeric sidecars, a provenance manifest, layout nodes, and speaker notes.

The output is deliberately a statistical candidate.  It does not claim the
full systematic covariance packet or numerator-denominator cross-covariance.
"""

from __future__ import annotations

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
from matplotlib.lines import Line2D
import numpy as np
import ROOT


ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
if ROOT.gSystem.Load("libRooUnfold") < 0:
    raise RuntimeError("failed to load libRooUnfold")

REPO = Path(__file__).resolve().parents[4]
OUT = REPO / "dataOutput/the219_friday_ppg_20260814/slide01_pp_money"
OUT.mkdir(parents=True, exist_ok=True)

PNG = OUT / "slide01_pp_particle_level_xjgamma_stat_candidate.png"
CSV = OUT / "slide01_pp_particle_level_xjgamma_points.csv"
POINTS = OUT / "slide01_pp_particle_level_xjgamma_points.json"
MANIFEST = OUT / "slide01_pp_particle_level_xjgamma_manifest.json"
LAYOUT = OUT / "slide01_pp_particle_level_xjgamma_layout_nodes.json"
NOTES = OUT / "slide01_pp_particle_level_xjgamma_speaker_notes.md"

DATA_POINTER = REPO / "dataOutput/current_recoiljets_artifacts/current/pp_data_merged/current.json"
SIM_POINTER = REPO / "dataOutput/current_recoiljets_artifacts/current/pp_sim_photonjet_merged/current.json"
HALFCLOSURE = REPO / "dataOutput/the219_friday_ppg_20260814/phase_a/iteration_halfclosure/the97_fig37_sim_halfclosure_iteration_stability_manifest.json"

DATA_TOP = "PPG12_scaledtrigger30"
SIM_TOP = "SIM"
ITERATIONS = 3
PT_GROUPS = [(16.0, 20.0), (20.0, 26.0), (26.0, 35.0)]
DISPLAY_XMIN = {(16.0, 20.0): 0.41, (20.0, 26.0): 0.29, (26.0, 35.0): 0.20}
CANON_PT_EDGES = [5, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 35]
RECO_PT_EDGES = [8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 35, 40]

OBJECTS = {
    "data_photon_reco": "h_unfoldRecoPho_pTgamma_ppg12obj",
    "sim_photon_reco": "h_unfoldRecoPho_pTgamma_ppg12obj",
    "sim_photon_truth": "h_unfoldTruthPho_pTgamma_ppg12obj",
    "sim_photon_response": "h2_unfoldResponsePho_pTgamma_ppg12obj",
    "data_xj_reco_a": "h2_unfoldReco_pTgamma_xJ_incl_r04",
    "data_xj_reco_c": "h2_unfoldReco_pTgamma_xJ_incl_sidebandC_r04",
    "sim_xj_reco": "h2_unfoldReco_pTgamma_xJ_incl_r04",
    "sim_xj_truth": "h2_unfoldTruth_pTgamma_xJ_incl_r04",
    "sim_xj_response": "h2_unfoldResponse_pTgamma_xJ_incl_r04",
}

LEGACY_PATH = REPO / "scripts/slides/the85_xjgamma_unfolding/make_unfolded_xjgamma_1x3.py"
spec = importlib.util.spec_from_file_location("the85_xj_helpers", LEGACY_PATH)
if spec is None or spec.loader is None:
    raise RuntimeError(f"cannot import helpers from {LEGACY_PATH}")
legacy = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = legacy
spec.loader.exec_module(legacy)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def json_ready(value: Any) -> Any:
    """Replace non-finite numeric sentinels with JSON null."""
    if isinstance(value, dict):
        return {key: json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_ready(item) for item in value]
    if isinstance(value, (float, np.floating)) and not math.isfinite(float(value)):
        return None
    return value


def open_root(path: Path) -> ROOT.TFile:
    f = ROOT.TFile.Open(str(path), "READ")
    if not f or f.IsZombie():
        raise RuntimeError(f"failed to open ROOT file: {path}")
    return f


def obj(f: ROOT.TFile, top: str, name: str, kind: str):
    source = f.Get(f"{top}/{name}")
    if not source:
        raise KeyError(f"missing {top}/{name} in {f.GetName()}")
    if kind == "TH1" and not source.InheritsFrom("TH1"):
        raise TypeError(f"{top}/{name} is not TH1")
    if kind == "TH2" and not source.InheritsFrom("TH2"):
        raise TypeError(f"{top}/{name} is not TH2")
    out = source.Clone(f"{name}_the219_clone")
    out.SetDirectory(0)
    out.Sumw2()
    return out


def scalar(f: ROOT.TFile, top: str, name: str) -> float:
    source = f.Get(f"{top}/{name}")
    return float(source.GetBinContent(1)) if source and source.InheritsFrom("TH1") else 0.0


def suffix_for_reco(lo: int, hi: int) -> str:
    if (lo, hi) == (8, 10):
        return "_pT_5_8"
    if (lo, hi) == (35, 40):
        return "_pT_26_35"
    return f"_pT_{lo}_{hi}"


def leakage(sim: ROOT.TFile, suffix: str) -> tuple[float, float, float, str]:
    preferred = f"h_xJpurityLead_sigABCD_MC{suffix}"
    fallback = f"h_sigABCD_MC{suffix}"
    source = sim.Get(f"{SIM_TOP}/{preferred}")
    used = preferred
    if not source:
        source = sim.Get(f"{SIM_TOP}/{fallback}")
        used = fallback
    if not source:
        return 0.0, 0.0, 0.0, "missing"
    a = float(source.GetBinContent(1))
    if a <= 0:
        return 0.0, 0.0, 0.0, used
    return (
        float(source.GetBinContent(2)) / a,
        float(source.GetBinContent(3)) / a,
        float(source.GetBinContent(4)) / a,
        used,
    )


def solve_leakage_sa(a: float, b: float, c: float, d: float, fb: float, fc: float, fd: float) -> tuple[bool, float]:
    if a <= 0:
        return True, 0.0
    s = min(max(a - b * c / d, 0.0), a) if d else a
    for _ in range(200):
        if fd > 0:
            s = min(s, max(0.0, 0.999 * d / fd))
        denom = d - fd * s
        if denom == 0:
            return False, 0.0
        fixed = a - (b - fb * s) * (c - fc * s) / denom
        if not math.isfinite(fixed):
            return False, 0.0
        next_s = min(max(0.75 * s + 0.25 * fixed, 0.0), a)
        if abs(next_s - s) < 1e-6:
            return True, next_s
        s = next_s
    return True, s


def abcd_row(data: ROOT.TFile, sim: ROOT.TFile, suffix: str) -> dict[str, Any]:
    a = scalar(data, DATA_TOP, f"h_xJpurityLead_isIsolated_isTight{suffix}")
    b = scalar(data, DATA_TOP, f"h_xJpurityLead_notIsolated_isTight{suffix}")
    c = scalar(data, DATA_TOP, f"h_xJpurityLead_isIsolated_notTight{suffix}")
    d = scalar(data, DATA_TOP, f"h_xJpurityLead_notIsolated_notTight{suffix}")
    fb, fc, fd, leak_source = leakage(sim, suffix)
    sa = min(max(a - b * c / d, 0.0), a) if a > 0 and d > 0 else 0.0
    ok, corrected = solve_leakage_sa(a, b, c, d, fb, fc, fd)
    if ok and (fb > 0 or fc > 0 or fd > 0):
        sa = min(max(corrected, 0.0), a)
    var_sa = max(a, 0.0)
    if d > 0:
        var_sa += (c / d) ** 2 * max(b, 0.0)
        var_sa += (b / d) ** 2 * max(c, 0.0)
        var_sa += (b * c / (d * d)) ** 2 * max(d, 0.0)
    esa = math.sqrt(max(var_sa, 0.0))
    nbkg_a = max(0.0, a - sa)
    c_bkg = max(0.0, c - fc * sa)
    scale_c = nbkg_a / c_bkg if c_bkg > 0 else 0.0
    denom = 1.0 - scale_c * fc
    inv_denom = 1.0 / denom if math.isfinite(denom) and abs(denom) > 1e-6 else 1.0
    var_scale = 0.0
    if c_bkg > 0 and nbkg_a > 0:
        var_nbkg = max(a, 0.0) + esa * esa
        var_scale = var_nbkg / (c_bkg * c_bkg) + (nbkg_a * nbkg_a * max(c_bkg, 0.0)) / (c_bkg ** 4)
    return {
        "A": a, "B": b, "C": c, "D": d, "SA": sa, "eSA": esa,
        "N_bkg_A": nbkg_a, "fB": fb, "fC": fc, "fD": fd,
        "leakage_source": leak_source, "C_bkg": c_bkg,
        "scaleC": scale_c, "varScaleC": max(var_scale, 0.0),
        "inverse_leakage_denominator": inv_denom,
    }


def pooled_c_shape(h_c, ix: int, lo: int, hi: int, purity: float) -> tuple[np.ndarray, np.ndarray] | None:
    sparse = lo >= 22 and hi <= 35 and (purity < 0.45 or hi >= 26)
    if not sparse:
        return None
    ny = h_c.GetNbinsY()
    row_raw = sum(h_c.GetBinContent(ix, iy) for iy in range(0, ny + 2))
    if row_raw <= 0:
        return None
    vals = np.zeros(ny + 2)
    err2 = np.zeros(ny + 2)
    total = 0.0
    for plo, phi in zip(RECO_PT_EDGES[:-1], RECO_PT_EDGES[1:]):
        if plo < 20 or phi > 35:
            continue
        ix_pool = h_c.GetXaxis().FindBin(0.5 * (plo + phi))
        for iy in range(0, ny + 2):
            v = h_c.GetBinContent(ix_pool, iy)
            e = h_c.GetBinError(ix_pool, iy)
            vals[iy] += v
            err2[iy] += e * e
            total += v
    if total <= 0:
        return None
    scale = row_raw / total
    return vals * scale, np.sqrt(err2) * scale


def purity_correct_photons(data: ROOT.TFile, sim: ROOT.TFile, reco) -> tuple[Any, list[dict[str, Any]]]:
    out = reco.Clone("the219_pp_photon_purity_input")
    out.SetDirectory(0)
    out.Reset("ICES")
    out.Sumw2()
    rows: list[dict[str, Any]] = []
    for lo, hi in zip(RECO_PT_EDGES[:-1], RECO_PT_EDGES[1:]):
        suffix = suffix_for_reco(lo, hi)
        row = abcd_row(data, sim, suffix)
        if row["A"] + row["B"] + row["C"] + row["D"] <= 0:
            continue
        ib = out.GetXaxis().FindBin(0.5 * (lo + hi))
        out.SetBinContent(ib, row["SA"])
        out.SetBinError(ib, row["eSA"])
        rows.append({"pt": [lo, hi], "suffix": suffix, **row})
    if not rows:
        raise RuntimeError("no event-leading ABCD counters available for photon unfolding input")
    return out, rows


def purity_correct_xj(data: ROOT.TFile, sim: ROOT.TFile, h_a, h_c) -> tuple[Any, list[dict[str, Any]]]:
    out = h_a.Clone("the219_pp_xj_purity_input")
    out.SetDirectory(0)
    out.Reset("ICES")
    out.Sumw2()
    ny = out.GetNbinsY()
    rows: list[dict[str, Any]] = []
    for lo, hi in zip(RECO_PT_EDGES[:-1], RECO_PT_EDGES[1:]):
        suffix = suffix_for_reco(lo, hi)
        row = abcd_row(data, sim, suffix)
        if row["A"] + row["B"] + row["C"] + row["D"] <= 0:
            continue
        ix = out.GetXaxis().FindBin(0.5 * (lo + hi))
        purity = row["SA"] / row["A"] if row["A"] > 0 else 0.0
        pooled = pooled_c_shape(h_c, ix, lo, hi, purity)
        for iy in range(0, ny + 2):
            av = h_a.GetBinContent(ix, iy)
            ae = h_a.GetBinError(ix, iy)
            cv = float(pooled[0][iy]) if pooled is not None else h_c.GetBinContent(ix, iy)
            ce = float(pooled[1][iy]) if pooled is not None else h_c.GetBinError(ix, iy)
            val = (av - row["scaleC"] * cv) * row["inverse_leakage_denominator"]
            var = ae * ae + row["scaleC"] ** 2 * ce * ce + cv * cv * row["varScaleC"]
            var *= row["inverse_leakage_denominator"] ** 2
            out.SetBinContent(ix, iy, val if math.isfinite(val) else 0.0)
            out.SetBinError(ix, iy, math.sqrt(max(var, 0.0)) if math.isfinite(var) else 0.0)
        rows.append({"pt": [lo, hi], "suffix": suffix, "purity": purity, "pooled_C_shape": pooled is not None, **row})
    if not rows:
        raise RuntimeError("no event-leading ABCD counters available for xJ unfolding input")
    return out, rows


def matrix_to_numpy(matrix, n: int) -> np.ndarray:
    try:
        return np.array([[float(matrix[i][j]) for j in range(n)] for i in range(n)], dtype=float)
    except Exception:
        return np.array([[float(matrix(i, j)) for j in range(n)] for i in range(n)], dtype=float)


def unfold_photons(data: ROOT.TFile, sim: ROOT.TFile) -> tuple[Any, np.ndarray, dict[str, Any]]:
    data_reco = obj(data, DATA_TOP, OBJECTS["data_photon_reco"], "TH1")
    sim_reco = obj(sim, SIM_TOP, OBJECTS["sim_photon_reco"], "TH1")
    sim_truth = obj(sim, SIM_TOP, OBJECTS["sim_photon_truth"], "TH1")
    response_raw = obj(sim, SIM_TOP, OBJECTS["sim_photon_response"], "TH2")
    corrected, rows = purity_correct_photons(data, sim, data_reco)
    response = legacy.transpose_th2(response_raw, "the219_photon_response_recoX_truthY")
    roo_response = ROOT.RooUnfoldResponse(sim_reco, sim_truth, response, "the219_photon_response", "the219_photon_response")
    unfold = ROOT.RooUnfoldBayes(roo_response, corrected, ITERATIONS)
    unfold.SetVerbose(0)
    result = unfold.Hreco(ROOT.RooUnfold.kCovariance)
    result.SetDirectory(0)
    cov = matrix_to_numpy(unfold.Ereco(ROOT.RooUnfold.kCovariance), result.GetNbinsX())
    return result, cov, {"abcd_rows": rows, "response_orientation": "truthX_recoY transposed to recoX_truthY"}


def unfold_xj(data: ROOT.TFile, sim: ROOT.TFile) -> tuple[Any, Any, np.ndarray, dict[str, Any]]:
    data_a = obj(data, DATA_TOP, OBJECTS["data_xj_reco_a"], "TH2")
    data_c = obj(data, DATA_TOP, OBJECTS["data_xj_reco_c"], "TH2")
    sim_reco = obj(sim, SIM_TOP, OBJECTS["sim_xj_reco"], "TH2")
    sim_truth = obj(sim, SIM_TOP, OBJECTS["sim_xj_truth"], "TH2")
    response_raw = obj(sim, SIM_TOP, OBJECTS["sim_xj_response"], "TH2")
    corrected, rows = purity_correct_xj(data, sim, data_a, data_c)
    data_global = legacy.flatten_th2_to_global(corrected, "the219_xj_data_global")
    reco_global = legacy.flatten_th2_to_global(sim_reco, "the219_xj_sim_reco_global")
    truth_global = legacy.flatten_th2_to_global(sim_truth, "the219_xj_sim_truth_global")
    response = legacy.transpose_th2(response_raw, "the219_xj_response_recoX_truthY")
    roo_response = ROOT.RooUnfoldResponse(reco_global, truth_global, response, "the219_xj_response", "the219_xj_response")
    unfold = ROOT.RooUnfoldBayes(roo_response, data_global, ITERATIONS)
    unfold.SetVerbose(0)
    result_global = unfold.Hreco(ROOT.RooUnfold.kCovariance)
    result_global.SetDirectory(0)
    cov = matrix_to_numpy(unfold.Ereco(ROOT.RooUnfold.kCovariance), result_global.GetNbinsX())
    result_2d = legacy.unflatten_global_to_th2(result_global, sim_truth, "the219_xj_truth_unfolded")
    return result_global, result_2d, cov, {"abcd_rows": rows, "response_orientation": "truthGlobalX_recoGlobalY transposed to recoX_truthY"}


def matching_bins(axis, lo: float, hi: float) -> list[int]:
    bins = []
    for ib in range(1, axis.GetNbins() + 1):
        blo = float(axis.GetBinLowEdge(ib))
        bhi = float(axis.GetBinUpEdge(ib))
        if blo >= lo - 1e-9 and bhi <= hi + 1e-9:
            bins.append(ib)
    covered = sum(float(axis.GetBinWidth(ib)) for ib in bins)
    if not bins or abs(covered - (hi - lo)) > 1e-6:
        raise RuntimeError(f"requested pT range {lo}-{hi} is not an exact union of native bins")
    return bins


def global_index(h2, ix: int, iy: int) -> int:
    return int(h2.GetBin(ix, iy))


def truth_panel(
    lo: float,
    hi: float,
    h_photon,
    cov_photon: np.ndarray,
    h_xj_global,
    h_xj_2d,
    cov_xj: np.ndarray,
    sim_truth_photon,
    sim_truth_xj,
) -> dict[str, Any]:
    pxbins = matching_bins(h_photon.GetXaxis(), lo, hi)
    xbins = matching_bins(h_xj_2d.GetXaxis(), lo, hi)
    sim_pxbins = matching_bins(sim_truth_photon.GetXaxis(), lo, hi)
    sim_xbins = matching_bins(sim_truth_xj.GetXaxis(), lo, hi)
    npho = sum(float(h_photon.GetBinContent(i)) for i in pxbins)
    pidx = [i - 1 for i in pxbins]
    var_npho = float(cov_photon[np.ix_(pidx, pidx)].sum())
    npho_truth = sum(float(sim_truth_photon.GetBinContent(i)) for i in sim_pxbins)
    yaxis = h_xj_2d.GetYaxis()
    edges = np.array([float(yaxis.GetBinLowEdge(i)) for i in range(1, yaxis.GetNbins() + 2)])
    widths = np.diff(edges)
    values, errors, truth, truth_errors = [], [], [], []
    raw_yield = 0.0
    raw_truth_yield = 0.0
    for iy in range(1, yaxis.GetNbins() + 1):
        gidx = [global_index(h_xj_2d, ix, iy) for ix in xbins]
        numerator = sum(float(h_xj_global.GetBinContent(g + 1)) for g in gidx)
        var_num = float(cov_xj[np.ix_(gidx, gidx)].sum())
        per_photon_bin = numerator / npho if npho > 0 else 0.0
        variance = var_num / (npho * npho) if npho > 0 else 0.0
        if npho > 0:
            variance += numerator * numerator * max(var_npho, 0.0) / (npho ** 4)
        raw_yield += per_photon_bin
        values.append(per_photon_bin / widths[iy - 1])
        errors.append(math.sqrt(max(variance, 0.0)) / widths[iy - 1])

        sim_num = sum(float(sim_truth_xj.GetBinContent(ix, iy)) for ix in sim_xbins)
        sim_var = sum(float(sim_truth_xj.GetBinError(ix, iy)) ** 2 for ix in sim_xbins)
        sim_bin = sim_num / npho_truth if npho_truth > 0 else 0.0
        raw_truth_yield += sim_bin
        truth.append(sim_bin / widths[iy - 1])
        truth_errors.append(math.sqrt(max(sim_var, 0.0)) / npho_truth / widths[iy - 1] if npho_truth > 0 else 0.0)
    values = np.asarray(values)
    errors = np.asarray(errors)
    truth = np.asarray(truth)
    truth_errors = np.asarray(truth_errors)
    ratio = np.divide(values, truth, out=np.full_like(values, np.nan), where=truth > 0)
    ratio_error = np.divide(errors, truth, out=np.full_like(errors, np.nan), where=truth > 0)
    return {
        "pt": [lo, hi], "xj_edges": edges.tolist(), "values": values.tolist(), "errors": errors.tolist(),
        "pythia_truth": truth.tolist(), "pythia_truth_errors": truth_errors.tolist(),
        "ratio_to_pythia": ratio.tolist(), "ratio_errors": ratio_error.tolist(),
        "unfolded_photons": npho, "unfolded_photon_variance": var_npho,
        "accepted_recoil_yield": raw_yield, "pythia_accepted_recoil_yield": raw_truth_yield,
        "truth_pt_bins": xbins, "photon_truth_bins": pxbins,
    }


def px_bbox(left: float, bottom: float, width: float, height: float) -> list[float]:
    return [2560 * left, 1440 * (1 - bottom - height), 2560 * (left + width), 1440 * (1 - bottom)]


def render(panels: list[dict[str, Any]]) -> None:
    mpl.rcParams.update({
        "font.family": "serif", "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
        "mathtext.fontset": "stix", "axes.linewidth": 1.0, "xtick.direction": "in", "ytick.direction": "in",
        "xtick.top": True, "ytick.right": True,
    })
    fig = plt.figure(figsize=(16, 9), dpi=160, facecolor="white")
    fig.text(0.045, 0.955, r"The PPG12 photon definition now reaches particle-level p+p $x_{J\gamma}$",
             ha="left", va="top", fontsize=28, fontweight="bold", color="#111111")
    fig.text(0.047, 0.895,
             r"Per particle-level photon  |  anti-$k_{T}$ $R=0.4$  |  $p_{T}^{jet}>5$ GeV  |  $|\eta^{\gamma,jet}|<0.7$  |  $|\Delta\phi|>7\pi/8$",
             ha="left", va="top", fontsize=16, color="#333333")
    fig.text(0.953, 0.910, "sPHENIX Internal", ha="right", va="top", fontsize=16,
             fontstyle="italic", fontweight="bold", color="#222222")

    left0, gap, total_w = 0.073, 0.035, 0.867
    panel_w = (total_w - 2 * gap) / 3
    main_bottom, main_height = 0.330, 0.460
    ratio_bottom, ratio_height = 0.190, 0.115
    magenta = "#c2185b"
    for ip, panel in enumerate(panels):
        left = left0 + ip * (panel_w + gap)
        ax = fig.add_axes([left, main_bottom, panel_w, main_height])
        rax = fig.add_axes([left, ratio_bottom, panel_w, ratio_height], sharex=ax)
        edges = np.asarray(panel["xj_edges"])
        centers = 0.5 * (edges[:-1] + edges[1:])
        halfwidth = 0.5 * np.diff(edges)
        values = np.asarray(panel["values"])
        errors = np.asarray(panel["errors"])
        truth = np.asarray(panel["pythia_truth"])
        truth_errors = np.asarray(panel["pythia_truth_errors"])
        ratio = np.asarray(panel["ratio_to_pythia"])
        ratio_error = np.asarray(panel["ratio_errors"])
        display_xmin = DISPLAY_XMIN[(float(panel["pt"][0]), float(panel["pt"][1]))]
        mask = (edges[:-1] >= display_xmin - 1e-9) & (centers <= 1.65)

        ax.errorbar(centers[mask], values[mask], yerr=errors[mask], xerr=halfwidth[mask], fmt="o",
                    color="black", markersize=4.8, capsize=2.0, linewidth=1.25, zorder=4)
        ax.errorbar(centers[mask], truth[mask], yerr=truth_errors[mask], xerr=halfwidth[mask], fmt="o",
                    mfc="white", mec=magenta, ecolor=magenta, color=magenta, markersize=5.1,
                    capsize=1.8, linewidth=1.05, zorder=3)
        ax.set_xlim(display_xmin, 1.65)
        panel_ymax = max(float(np.max(values + errors)), float(np.max(truth + truth_errors)))
        panel_ymax = max(0.8, math.ceil(1.12 * panel_ymax * 2.0) / 2.0)
        ax.set_ylim(0.0, panel_ymax)
        ax.grid(axis="y", color="#dedede", linewidth=0.65, alpha=0.75)
        ax.tick_params(labelsize=11, length=4)
        ax.tick_params(labelbottom=False)
        lo, hi = panel["pt"]
        ax.text(0.05, 0.94, rf"${lo:.0f}<p_{{T}}^{{\gamma}}<{hi:.0f}$ GeV", transform=ax.transAxes,
                ha="left", va="top", fontsize=15, fontweight="bold")
        ax.text(0.05, 0.83, rf"shown for $x_{{J\gamma}}\geq {display_xmin:.2f}$", transform=ax.transAxes,
                ha="left", va="top", fontsize=12.5, color="#333333")
        if ip == 0:
            ax.set_ylabel(r"$\frac{1}{N_{\gamma}^{particle}}\,\frac{dN_{jet}^{particle}}{dx_{J\gamma}}$",
                          fontsize=17, labelpad=9)
            handles = [
                Line2D([0], [0], marker="o", color="black", linestyle="none", markersize=5.5,
                       label="p+p data, unfolded"),
                Line2D([0], [0], marker="o", markerfacecolor="white", markeredgecolor=magenta,
                       color=magenta, linestyle="none", markersize=6, label="PYTHIA-8 truth"),
            ]
            ax.legend(handles=handles, loc="upper right", frameon=False, fontsize=11.8,
                      handletextpad=0.4, borderaxespad=0.4)
        else:
            ax.tick_params(labelleft=False)

        valid = mask & np.isfinite(ratio) & (truth > 0.02)
        rax.axhline(1.0, color="#666666", linewidth=1.0, linestyle="--")
        rax.errorbar(centers[valid], ratio[valid], yerr=ratio_error[valid], xerr=halfwidth[valid], fmt="o",
                     color="black", markersize=4.3, capsize=1.7, linewidth=1.0)
        rax.set_xlim(display_xmin, 1.65)
        rax.set_ylim(0.0, 5.2)
        rax.set_yticks([1.0, 3.0, 5.0])
        rax.grid(axis="y", color="#e4e4e4", linewidth=0.6)
        rax.tick_params(labelsize=10.5, length=3.5)
        rax.set_xlabel(r"$x_{J\gamma}=p_{T}^{jet}/p_{T}^{\gamma}$", fontsize=14, labelpad=3)
        if ip == 0:
            rax.set_ylabel("Data /\nPYTHIA", fontsize=12.5, labelpad=8)
        else:
            rax.tick_params(labelleft=False)

    fig.patches.extend([
        mpl.patches.FancyBboxPatch((0.047, 0.094), 0.906, 0.050, transform=fig.transFigure,
                                   boxstyle="round,pad=0.006,rounding_size=0.006", facecolor="#fff4cf",
                                   edgecolor="#d9a928", linewidth=1.0),
        mpl.patches.FancyBboxPatch((0.047, 0.022), 0.906, 0.052, transform=fig.transFigure,
                                   boxstyle="round,pad=0.006,rounding_size=0.006", facecolor="#eaf3fb",
                                   edgecolor="#7ca7c8", linewidth=1.0),
    ])
    fig.text(0.5, 0.119, "Statistical candidate — unfolding + ABCD input covariance shown; full systematic covariance not yet applied",
             ha="center", va="center", fontsize=15.5, fontweight="bold", color="#6b4f00")
    fig.text(0.5, 0.048,
             r"Why this is the Au+Au baseline: per-photon normalization preserves the accepted recoil yield needed for $I_{AA}(x_{J\gamma})$.",
             ha="center", va="center", fontsize=16, fontweight="bold", color="#153b57")
    fig.savefig(PNG, dpi=160, facecolor="white")
    plt.close(fig)


def write_layout() -> None:
    title = "The PPG12 photon definition now reaches particle-level p+p xJgamma"
    nodes: list[dict[str, Any]] = [
        {"name": "claim title", "kind": "text", "role": "title", "font_px": 28 * 160 / 72,
         "bbox": [115, 52, 2120, 112], "text": title},
        {"name": "selection line", "kind": "text", "role": "audience", "font_px": 16 * 160 / 72,
         "bbox": [120, 142, 2250, 190], "text": "Per particle-level photon; anti-kT R=0.4; jet pT > 5 GeV; eta and back-to-back cuts"},
        {"name": "status caveat", "kind": "text", "role": "audience", "font_px": 15.5 * 160 / 72,
         "bbox": [145, 1260, 2415, 1330], "text": "Statistical candidate - unfolding and ABCD input covariance shown; full systematic covariance not yet applied"},
        {"name": "AuAu bridge takeaway", "kind": "text", "role": "audience", "font_px": 16 * 160 / 72,
         "bbox": [145, 1340, 2415, 1425], "text": "Why this is the AuAu baseline - per-photon normalization preserves the accepted recoil yield needed for IAA(xJgamma)."},
    ]
    left0, gap, total_w = 0.073, 0.035, 0.867
    panel_w = (total_w - 2 * gap) / 3
    for idx in range(3):
        left = left0 + idx * (panel_w + gap)
        nodes.append({"name": f"result panel {idx + 1}", "kind": "panel", "role": "audience",
                      "symmetry_group": "three result panels", "bbox": px_bbox(left, 0.165, panel_w, 0.625)})
    LAYOUT.write_text(json.dumps({
        "schema": "slide_layout_nodes_v1", "slide": "the219_pp_particle_level_xjgamma_stat_candidate",
        "title_axis_x": 115, "minimum_audience_font_px": 33,
        "minimum_plot_annotation_font_px": 24, "minimum_title_font_px": 60,
        "nodes": nodes,
    }, indent=2) + "\n")


def main() -> None:
    data_pointer = json.loads(DATA_POINTER.read_text())
    sim_pointer = json.loads(SIM_POINTER.read_text())
    halfclosure = json.loads(HALFCLOSURE.read_text())
    if int(halfclosure["best_iteration"]) != ITERATIONS:
        raise RuntimeError("configured iteration does not match independent current-input half-closure")
    data_path = Path(data_pointer["root_paths"][0])
    sim_path = Path(sim_pointer["root_paths"][0])
    data = open_root(data_path)
    sim = open_root(sim_path)
    try:
        h_photon, cov_photon, photon_meta = unfold_photons(data, sim)
        h_xj_global, h_xj_2d, cov_xj, xj_meta = unfold_xj(data, sim)
        sim_truth_photon = obj(sim, SIM_TOP, OBJECTS["sim_photon_truth"], "TH1")
        sim_truth_xj = obj(sim, SIM_TOP, OBJECTS["sim_xj_truth"], "TH2")
        panels = [truth_panel(lo, hi, h_photon, cov_photon, h_xj_global, h_xj_2d, cov_xj,
                              sim_truth_photon, sim_truth_xj) for lo, hi in PT_GROUPS]
    finally:
        data.Close()
        sim.Close()

    render(panels)
    write_layout()
    POINTS.write_text(json.dumps(json_ready(panels), indent=2, allow_nan=False) + "\n")
    with CSV.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=["pt_low", "pt_high", "xj_low", "xj_high", "value", "stat_error", "pythia_truth", "pythia_truth_stat_error", "data_over_pythia", "ratio_stat_error"])
        writer.writeheader()
        for panel in panels:
            edges = panel["xj_edges"]
            for i in range(len(edges) - 1):
                writer.writerow({
                    "pt_low": panel["pt"][0], "pt_high": panel["pt"][1], "xj_low": edges[i], "xj_high": edges[i + 1],
                    "value": panel["values"][i], "stat_error": panel["errors"][i],
                    "pythia_truth": panel["pythia_truth"][i], "pythia_truth_stat_error": panel["pythia_truth_errors"][i],
                    "data_over_pythia": panel["ratio_to_pythia"][i], "ratio_stat_error": panel["ratio_errors"][i],
                })

    NOTES.write_text(
        "# Slide 1 speaker notes\n\n"
        "This is the first end-to-end p+p particle-level recoil result from our current qualified inputs. "
        "The photon definition is not being retuned: it is the PPG12-parity object already established in the purity work. "
        "The numerator is the R=0.4 joint pT-gamma versus xJ-gamma spectrum unfolded with three Bayesian iterations; "
        "the denominator is the independently unfolded particle-level photon count in the same exact native-bin union.\n\n"
        "The black points are unfolded p+p data and the open magenta points are the matched PYTHIA-8 truth reference. "
        "The integral printed in each panel is the accepted recoil yield per particle-level photon, so the spectrum is not unit-normalized. "
        "That is deliberate: the same observable carries the rate information required for the Au+Au IAA comparison.\n\n"
        "Do not call this final. The bars include RooUnfold statistical covariance for the joint numerator and photon denominator, "
        "plus the current binwise ABCD-input uncertainty. The numerator-denominator cross-covariance, xJ-dependent ABCD transfer covariance, "
        "and the complete detector, photon, jet, unfolding, and model systematic covariance are not yet included. "
        "The next result slide should therefore show the closure and covariance gates that turn this candidate into the final baseline.\n"
    )

    manifest = {
        "status": "statistical_candidate_not_final",
        "slide": str(PNG), "slide_sha256": sha256(PNG), "points_json": str(POINTS), "points_csv": str(CSV),
        "layout_nodes": str(LAYOUT), "speaker_notes": str(NOTES),
        "canonical_deck_reference": "https://docs.google.com/presentation/d/19pjXOc1CJ1Ed7xtw2ZrXffztpWkzCz4lN13L7CHPdQE/edit",
        "google_slides_mutated": False,
        "inputs": {
            "data_pointer": str(DATA_POINTER), "data_pointer_sha256": sha256(DATA_POINTER),
            "data_campaign": data_pointer["campaign_tag"], "data_status": data_pointer["canonical_status"],
            "data_root": str(data_path), "data_root_sha256": sha256(data_path),
            "sim_pointer": str(SIM_POINTER), "sim_pointer_sha256": sha256(SIM_POINTER),
            "sim_campaign": sim_pointer["campaign_tag"], "sim_status": sim_pointer["canonical_status"],
            "sim_root": str(sim_path), "sim_root_sha256": sha256(sim_path),
        },
        "objects": OBJECTS,
        "analysis": {
            "photon_definition": "PPG12 parity object (ppg12obj)", "jet_radius": 0.4, "jet_pt_min_gev": 5.0,
            "photon_eta_abs_max": 0.7, "jet_eta_abs_max": 0.7, "minimum_abs_delta_phi": "7pi/8",
            "pt_groups_gev": PT_GROUPS, "display_xj_min": {f"{lo:g}-{hi:g}": DISPLAY_XMIN[(lo, hi)] for lo, hi in PT_GROUPS},
            "display_xj_min_basis": "Sam-compatible pT-dependent fiducial display; closest exact current native xJ edges. Turn-on bins remain in sidecars.",
            "iterations": ITERATIONS, "iteration_source": str(HALFCLOSURE),
            "unfolding": "RooUnfoldBayes with kCovariance", "normalization": "per unfolded particle-level photon; not unit area",
            "purity": "event-leading ABCD counts with preferred event-leading SIM leakage, sideband-C subtraction, high-pT C-shape pooling",
            "numerator_covariance": "full RooUnfold joint-spectrum statistical covariance aggregated across exact truth-pT bins",
            "denominator_covariance": "full RooUnfold photon statistical covariance aggregated across exact truth-pT bins",
        },
        "caveats": [
            "Full systematic covariance is not applied.",
            "Numerator-denominator cross-covariance is not available and is not included.",
            "The ABCD scale variance is propagated binwise into the measured input, but a final xJ-dependent transfer covariance is not yet included.",
            "The p+p data current pointer is candidate_evidence/current-default, not final scientific canonicalization.",
            "Independent unfolded closure and complete systematic variations remain required before final-result language.",
        ],
        "photon_unfolding_meta": photon_meta, "xj_unfolding_meta": xj_meta, "panels": panels,
    }
    MANIFEST.write_text(json.dumps(json_ready(manifest), indent=2, allow_nan=False) + "\n")
    print(PNG)
    print(MANIFEST)
    print(LAYOUT)
    print(NOTES)


if __name__ == "__main__":
    main()
