#!/usr/bin/env python3
"""
THE-85 first-pass unfolded xJgamma slide.

This is a scoped plotting helper for the current AuAu default-BDT / pp baseline
outputs.  It deliberately mirrors the RooUnfold object semantics used by
macros/AnalyzeRecoilJets_RooUnfoldPipeline.cpp:

  - photon normalization: unfold N_gamma(pTgamma) with the SIM photon response
  - xJ: unfold the flattened (pTgamma, xJ) global-bin spectrum with the SIM
        response, then project truth xJ over 15 < pTgamma < 35 GeV
  - measured xJ input: region-A data with the same ABCD sideband-C subtraction
        shape used in the C++ pipeline when h_xJpurityLead_* counters are present

It is not a replacement for the full pipeline macro.  It produces a transparent
first-pass slide candidate and a manifest for review.
"""

from __future__ import annotations

import json
import math
import os
from dataclasses import dataclass, asdict, replace
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import ROOT


ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gSystem.Load("libRooUnfold")


REPO = Path(__file__).resolve().parents[3]
OUT_DIR = REPO / "dataOutput/the85_auau_xjgamma_unfolding_push/unfolded_xjgamma_firstpass"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# The nominal analysis back-to-back cut is 7pi/8.  In the RecoilJets key
# convention that keeps the legacy rKey name with no dphi suffix; _dphiPi2 is
# the looser pi/2 diagnostic variant.
BASE_KEY = "r04_isoR40_isSliding"
PHO_KEY = "isoR40_isSliding"
PT_EDGES_CANON = [5, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 35]
PT_EDGES_UNFOLD_RECO = [8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 35, 40]
PT_WINDOW = (15.0, 35.0)
REQUIRE_FULL_PT_BINS = os.environ.get("THE85_REQUIRE_FULL_PT_BINS", "0").strip() in {"1", "true", "TRUE", "yes", "YES"}
DEFAULT_ITERS = 5
ERROR_MODE_NAME = "kCovariance"
ERROR_MODE = ROOT.RooUnfold.kCovariance
NTOYS_FINAL = 600
OUTPUT_TAG = "iter5_covariance"
TAIL_SHAPE_XMIN = 0.5


def row_in_pt_window(lo: float, hi: float, cen: float) -> bool:
    """Select photon-pT rows for the nominal analysis window.

    The legacy center-based selector includes the [14,16] GeV row for a
    nominal 15-35 GeV analysis.  Set THE85_REQUIRE_FULL_PT_BINS=1 for
    ATLAS-style diagnostic plots that must avoid threshold-bin leakage.
    """
    if REQUIRE_FULL_PT_BINS:
        return lo >= PT_WINDOW[0] and hi <= PT_WINDOW[1]
    return PT_WINDOW[0] <= cen < PT_WINDOW[1]


@dataclass
class Case:
    key: str
    title: str
    short_label: str
    data_file: str
    data_topdir: str
    sim_file: str
    sim_topdir: str
    cent_suffix: str
    color: str
    marker: str
    apply_abcd: bool = True
    apply_combinatoric_subtraction: bool = False


CASES = [
    Case(
        key="auau_0_20",
        title="Au+Au 0-20%",
        short_label="0-20%",
        data_file=str(REPO / "InputFiles/the69_default_auau_physicsqa/RecoilJets_auau_ALL_preselectionNewPPG12_tightAuAuCentInputBase3x3BDT_nonTightAuAuBDTComplement_baseVariant.root"),
        data_topdir="photon_12_plus_MBD_NS_geq_2_vtx_lt_150",
        sim_file=str(REPO / "InputFiles/the69_leakageCentWP_fix/RecoilJets_embeddedPhoton12plus20_MERGED.root"),
        sim_topdir="SIM",
        cent_suffix="_cent_0_20",
        color="#1f77b4",
        marker="o",
        apply_combinatoric_subtraction=True,
    ),
    Case(
        key="auau_50_80",
        title="Au+Au 50-80%",
        short_label="50-80%",
        data_file=str(REPO / "InputFiles/the69_default_auau_physicsqa/RecoilJets_auau_ALL_preselectionNewPPG12_tightAuAuCentInputBase3x3BDT_nonTightAuAuBDTComplement_baseVariant.root"),
        data_topdir="photon_12_plus_MBD_NS_geq_2_vtx_lt_150",
        sim_file=str(REPO / "InputFiles/the69_leakageCentWP_fix/RecoilJets_embeddedPhoton12plus20_MERGED.root"),
        sim_topdir="SIM",
        cent_suffix="_cent_50_80",
        color="#2ca02c",
        marker="s",
        apply_combinatoric_subtraction=True,
    ),
    Case(
        key="pp_basev3e",
        title="p+p baseV3E",
        short_label="p+p raw-A",
        data_file=str(REPO / "dataOutput/ppg12TableQA/THE42_ppg12_tableqa_v1_basev3e_20260611/merged_roots/RecoilJets_pp_ALL_preselectionNewPPG12_tightNewPPG12_nonTightNewPPG12.root"),
        data_topdir="Photon_4_GeV_plus_MBD_NS_geq_1",
        sim_file=str(REPO / "dataOutput/ppg12TableQA/THE42_ppg12_tableqa_v1_basev3e_20260611/merged_roots/RecoilJets_photonjet5plus10plus20_MERGED.root"),
        sim_topdir="SIM",
        cent_suffix="",
        color="#d62728",
        marker="D",
        apply_abcd=False,
        apply_combinatoric_subtraction=False,
    ),
]


def open_root(path: str) -> ROOT.TFile:
    f = ROOT.TFile.Open(path, "READ")
    if not f or f.IsZombie():
        raise RuntimeError(f"failed to open ROOT file: {path}")
    return f


def get_obj(f: ROOT.TFile, topdir: str, name: str, cls_name: str):
    obj = f.Get(f"{topdir}/{name}")
    if not obj:
        raise KeyError(f"missing {topdir}/{name} in {f.GetName()}")
    if cls_name == "TH1" and not obj.InheritsFrom("TH1"):
        raise TypeError(f"{name} is not TH1")
    if cls_name == "TH2" and not obj.InheritsFrom("TH2"):
        raise TypeError(f"{name} is not TH2")
    out = obj.Clone(f"{name}_clone")
    out.SetDirectory(0)
    out.Sumw2()
    return out


def get_optional(f: ROOT.TFile, topdir: str, name: str):
    obj = f.Get(f"{topdir}/{name}")
    if not obj:
        return None
    out = obj.Clone(f"{name}_clone")
    out.SetDirectory(0)
    out.Sumw2()
    return out


def axis_edges(axis) -> np.ndarray:
    return np.array([axis.GetBinLowEdge(i) for i in range(1, axis.GetNbins() + 2)], dtype=float)


def axes_match(a, b, tol: float = 1e-8) -> bool:
    if a.GetNbins() != b.GetNbins():
        return False
    ea = axis_edges(a)
    eb = axis_edges(b)
    return bool(np.allclose(ea, eb, rtol=0, atol=tol))


def clone_th2_like(src, name: str):
    h = src.Clone(name)
    h.SetDirectory(0)
    h.Reset("ICES")
    h.Sumw2()
    return h


def transpose_th2(src, name: str):
    x = src.GetXaxis()
    y = src.GetYaxis()
    hx = axis_edges(x)
    hy = axis_edges(y)
    out = ROOT.TH2D(name, "", y.GetNbins(), hy, x.GetNbins(), hx)
    out.SetDirectory(0)
    out.Sumw2()
    for ix in range(0, x.GetNbins() + 2):
        for iy in range(0, y.GetNbins() + 2):
            out.SetBinContent(iy, ix, src.GetBinContent(ix, iy))
            out.SetBinError(iy, ix, src.GetBinError(ix, iy))
    return out


def orient_response(response, measured, truth, name: str):
    if axes_match(response.GetXaxis(), measured.GetXaxis()) and axes_match(response.GetYaxis(), truth.GetXaxis()):
        h = response.Clone(name)
        h.SetDirectory(0)
        return h, "as_is_measuredX_truthY"
    if axes_match(response.GetYaxis(), measured.GetXaxis()) and axes_match(response.GetXaxis(), truth.GetXaxis()):
        return transpose_th2(response, name), "transposed_to_measuredX_truthY"
    raise RuntimeError(
        "response axes do not match measured/truth axes: "
        f"rsp=({response.GetXaxis().GetNbins()},{response.GetYaxis().GetNbins()}) "
        f"meas={measured.GetXaxis().GetNbins()} truth={truth.GetXaxis().GetNbins()}"
    )


def flatten_th2_to_global(h2, name: str):
    nx = h2.GetXaxis().GetNbins()
    ny = h2.GetYaxis().GetNbins()
    n = (nx + 2) * (ny + 2)
    h = ROOT.TH1D(name, "", n, -0.5, n - 0.5)
    h.SetDirectory(0)
    h.Sumw2()
    for ix in range(0, nx + 2):
        for iy in range(0, ny + 2):
            g = h2.GetBin(ix, iy)
            h.SetBinContent(g + 1, h2.GetBinContent(ix, iy))
            h.SetBinError(g + 1, h2.GetBinError(ix, iy))
    return h


def unflatten_global_to_th2(hglob, tmpl, name: str):
    h2 = clone_th2_like(tmpl, name)
    nx = h2.GetXaxis().GetNbins()
    ny = h2.GetYaxis().GetNbins()
    for ix in range(0, nx + 2):
        for iy in range(0, ny + 2):
            g = h2.GetBin(ix, iy)
            b = g + 1
            if 1 <= b <= hglob.GetXaxis().GetNbins():
                h2.SetBinContent(ix, iy, hglob.GetBinContent(b))
                h2.SetBinError(ix, iy, hglob.GetBinError(b))
    return h2


def pt_bins_from_edges(edges: Sequence[float]) -> List[Tuple[int, int, str]]:
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        ilo, ihi = int(round(lo)), int(round(hi))
        out.append((ilo, ihi, f"_pT_{ilo}_{ihi}"))
    return out


def read_scalar_count(f: ROOT.TFile, topdir: str, name: str) -> float:
    obj = f.Get(f"{topdir}/{name}")
    if not obj or not obj.InheritsFrom("TH1"):
        return 0.0
    return float(obj.GetBinContent(1))


def solve_leakage_corrected_sa(A: float, B: float, C: float, D: float, fB: float, fC: float, fD: float) -> Tuple[bool, float]:
    if A <= 0.0:
        return True, 0.0
    s = A
    if D != 0.0:
        s = min(max(A - B * (C / D), 0.0), A)

    def fixed_point(x: float) -> float:
        denom = D - fD * x
        if denom == 0.0:
            return float("nan")
        return A - (B - fB * x) * (C - fC * x) / denom

    lam = 0.25
    for _ in range(200):
        if fD > 0.0:
            smax = (D / fD) * 0.999
            if math.isfinite(smax):
                s = min(s, max(0.0, smax))
        fp = fixed_point(s)
        if not math.isfinite(fp):
            return False, 0.0
        sn = (1.0 - lam) * s + lam * fp
        if not math.isfinite(sn):
            return False, 0.0
        sn = min(max(sn, 0.0), A)
        if abs(sn - s) < 1e-6:
            return True, sn
        s = sn
    return True, s


def load_leakage_factors(sim_f: ROOT.TFile, topdir: str, cent_suffix: str) -> Dict[str, Tuple[float, float, float]]:
    out: Dict[str, Tuple[float, float, float]] = {}
    for lo, hi, suffix in pt_bins_from_edges(PT_EDGES_CANON):
        name = f"h_sigABCD_MC_{PHO_KEY}{suffix}{cent_suffix}"
        h = sim_f.Get(f"{topdir}/{name}")
        if not h:
            out[suffix] = (0.0, 0.0, 0.0)
            continue
        A = float(h.GetBinContent(1))
        B = float(h.GetBinContent(2))
        C = float(h.GetBinContent(3))
        D = float(h.GetBinContent(4))
        out[suffix] = ((B / A) if A > 0 else 0.0, (C / A) if A > 0 else 0.0, (D / A) if A > 0 else 0.0)
    return out


def compute_abcd_counts(
    f: ROOT.TFile,
    topdir: str,
    suffix: str,
    leakage: Optional[Tuple[float, float, float]] = None,
) -> Tuple[float, float, float, float, float, float, float, Dict]:
    a = read_scalar_count(f, topdir, "h_xJpurityLead_isIsolated_isTight" + suffix)
    b = read_scalar_count(f, topdir, "h_xJpurityLead_notIsolated_isTight" + suffix)
    c = read_scalar_count(f, topdir, "h_xJpurityLead_isIsolated_notTight" + suffix)
    d = read_scalar_count(f, topdir, "h_xJpurityLead_notIsolated_notTight" + suffix)
    sa_raw = 0.0
    if a > 0.0 and d > 0.0:
        sa_raw = max(0.0, a - b * (c / d))
    sa = sa_raw
    leakage_applied = False
    leakage_ok = True
    fB, fC, fD = leakage if leakage is not None else (0.0, 0.0, 0.0)
    if leakage is not None and (fB > 0.0 or fC > 0.0 or fD > 0.0):
        leakage_ok, sa_corr = solve_leakage_corrected_sa(a, b, c, d, fB, fC, fD)
        if leakage_ok:
            sa = min(max(sa_corr, 0.0), a)
            leakage_applied = True
    var = max(0.0, a)
    if d > 0.0:
        var += (c / d) ** 2 * max(0.0, b)
        var += (b / d) ** 2 * max(0.0, c)
        var += ((b * c) / (d * d)) ** 2 * max(0.0, d)
    esa = math.sqrt(var) if var > 0.0 else 0.0
    nbkg = max(0.0, a - sa)
    meta = {"SA_raw": sa_raw, "SA": sa, "leakage_applied": leakage_applied, "leakage_ok": leakage_ok, "fB": fB, "fC": fC, "fD": fD}
    return a, b, c, d, sa, esa, nbkg, meta


def abcd_suffix(pt_suffix: str, cent_suffix: str) -> str:
    return f"_{PHO_KEY}{pt_suffix}{cent_suffix}"


def canonical_suffix_for_reco_bin(lo: int, hi: int) -> str:
    suffix = f"_pT_{lo}_{hi}"
    if (lo, hi) == (8, 10):
        return "_pT_5_8"
    if (lo, hi) == (35, 40):
        return "_pT_26_35"
    return suffix


def apply_photon_abcd_input(case: Case, data_file: ROOT.TFile, sim_file: ROOT.TFile, h_reco_pho) -> Tuple[object, Dict]:
    if not case.apply_abcd:
        return h_reco_pho, {"applied": False, "reason": "disabled for this first-pass pp baseline; local THE42 pp ABCD counters are sparse in the selected pT window"}
    out = h_reco_pho.Clone(f"{case.key}_hPhoRecoData_purityInput")
    out.SetDirectory(0)
    out.Reset("ICES")
    out.Sumw2()
    any_counts = False
    pt_rows = []
    leakage_factors = load_leakage_factors(sim_file, case.sim_topdir, case.cent_suffix)
    for lo, hi, suffix in pt_bins_from_edges(PT_EDGES_UNFOLD_RECO):
        # Match the C++ mapping for support bins.
        canon_suffix = canonical_suffix_for_reco_bin(lo, hi)
        a, b, c, d, sa, esa, nbkg, row_meta = compute_abcd_counts(
            data_file,
            case.data_topdir,
            abcd_suffix(canon_suffix, case.cent_suffix),
            leakage_factors.get(canon_suffix),
        )
        if a + b + c + d <= 0.0:
            continue
        any_counts = True
        cen = 0.5 * (lo + hi)
        ib = out.GetXaxis().FindBin(cen)
        out.SetBinContent(ib, sa)
        out.SetBinError(ib, esa)
        pt_rows.append({"pt": [lo, hi], "A": a, "B": b, "C": c, "D": d, "eSA": esa, "nbkgA": nbkg, **row_meta})
    if not any_counts:
        return h_reco_pho, {"applied": False, "reason": "missing h_xJpurityLead counters"}
    return out, {"applied": True, "method": "ABCD S_A per pT bin from h_xJpurityLead counters", "pt_rows": pt_rows}


def pooled_sideband_c(h_side_c, ix: int, reco_bins: Sequence[Tuple[int, int, str]], b_lo: int, b_hi: int, purity: float):
    sparse = b_lo >= 22 and b_hi <= 35 and (purity < 0.45 or b_hi >= 26)
    if not sparse:
        return None
    ny = h_side_c.GetYaxis().GetNbins()
    raw = sum(h_side_c.GetBinContent(ix, iy) for iy in range(0, ny + 2))
    if raw <= 0.0:
        return None
    vals = np.zeros(ny + 2)
    err2 = np.zeros(ny + 2)
    total = 0.0
    for lo, hi, _suffix in reco_bins:
        if lo < 20 or hi > 35:
            continue
        cen = 0.5 * (lo + hi)
        ixp = h_side_c.GetXaxis().FindBin(cen)
        for iy in range(0, ny + 2):
            v = h_side_c.GetBinContent(ixp, iy)
            e = h_side_c.GetBinError(ixp, iy)
            vals[iy] += v
            err2[iy] += e * e
            total += v
    if total <= 0.0:
        return None
    scale = raw / total
    return vals * scale, np.sqrt(err2) * scale


def apply_xj_abcd_input(case: Case, data_file: ROOT.TFile, sim_file: ROOT.TFile, h_reco_a, h_side_c) -> Tuple[object, Dict]:
    if not case.apply_abcd:
        return h_reco_a, {"applied": False, "reason": "disabled for this first-pass pp baseline; local THE42 pp ABCD counters are sparse in the selected pT window"}
    if h_side_c is None:
        return h_reco_a, {"applied": False, "reason": "missing sideband-C xJ histogram"}
    out = h_reco_a.Clone(f"{case.key}_h2RecoData_purityInput")
    out.SetDirectory(0)
    out.Reset("ICES")
    out.Sumw2()
    reco_bins = pt_bins_from_edges(PT_EDGES_UNFOLD_RECO)
    ny = out.GetYaxis().GetNbins()
    any_counts = False
    rows = []
    leakage_factors = load_leakage_factors(sim_file, case.sim_topdir, case.cent_suffix)
    for lo, hi, suffix in reco_bins:
        canon_suffix = canonical_suffix_for_reco_bin(lo, hi)
        leakage = leakage_factors.get(canon_suffix)
        a, b, c, d, sa, esa, nbkg, row_meta = compute_abcd_counts(
            data_file,
            case.data_topdir,
            abcd_suffix(canon_suffix, case.cent_suffix),
            leakage,
        )
        if a + b + c + d <= 0.0:
            continue
        any_counts = True
        purity = sa / a if a > 0 else 0.0
        fC_leak = max(0.0, leakage[1]) if leakage is not None else 0.0
        c_bkg = max(0.0, c - fC_leak * sa)
        scale_c = nbkg / c_bkg if c_bkg > 0 else 0.0
        denom_sig_leak = 1.0 - scale_c * fC_leak
        inv_denom = 1.0 / denom_sig_leak if math.isfinite(denom_sig_leak) and abs(denom_sig_leak) > 1e-6 else 1.0
        cen = 0.5 * (lo + hi)
        ix = out.GetXaxis().FindBin(cen)
        pooled = pooled_sideband_c(h_side_c, ix, reco_bins, lo, hi, purity)
        for iy in range(0, ny + 2):
            a_val = h_reco_a.GetBinContent(ix, iy)
            a_err = h_reco_a.GetBinError(ix, iy)
            if pooled is not None:
                c_val = float(pooled[0][iy])
                c_err = float(pooled[1][iy])
            else:
                c_val = h_side_c.GetBinContent(ix, iy)
                c_err = h_side_c.GetBinError(ix, iy)
            val = (a_val - scale_c * c_val) * inv_denom
            err = math.sqrt(max(0.0, a_err * a_err + scale_c * scale_c * c_err * c_err)) * abs(inv_denom)
            if not math.isfinite(val):
                val = 0.0
            if not math.isfinite(err):
                err = 0.0
            out.SetBinContent(ix, iy, val)
            out.SetBinError(ix, iy, err)
        rows.append({"pt": [lo, hi], "A": a, "B": b, "C": c, "D": d, "purity": purity, "scaleC": scale_c, "fCLeak": fC_leak, "pooled_C_shape": pooled is not None, **row_meta})
    if not any_counts:
        return h_reco_a, {"applied": False, "reason": "missing h_xJpurityLead counters"}
    return out, {"applied": True, "method": "H_A - (A-S_A)/C * H_C, with high-pT C-shape pooling as in the C++ pipeline", "pt_rows": rows}


def unfold_photons(case: Case, data_file: ROOT.TFile, sim_file: ROOT.TFile) -> Tuple[object, Dict]:
    reco_data = get_obj(data_file, case.data_topdir, f"h_unfoldRecoPho_pTgamma_{PHO_KEY}{case.cent_suffix}", "TH1")
    reco_sim = get_obj(sim_file, case.sim_topdir, f"h_unfoldRecoPho_pTgamma_{PHO_KEY}{case.cent_suffix}", "TH1")
    truth_sim = get_obj(sim_file, case.sim_topdir, f"h_unfoldTruthPho_pTgamma{case.cent_suffix}", "TH1")
    rsp_sim = get_obj(sim_file, case.sim_topdir, f"h2_unfoldResponsePho_pTgamma_{PHO_KEY}{case.cent_suffix}", "TH2")
    reco_input, purity_meta = apply_photon_abcd_input(case, data_file, sim_file, reco_data)
    rsp = transpose_th2(rsp_sim, f"{case.key}_pho_rsp_recoX_truthY")
    orientation = "transposed_truthX_recoY_to_recoX_truthY"
    resp = ROOT.RooUnfoldResponse(reco_sim, truth_sim, rsp, f"{case.key}_respPho", f"{case.key}_respPho")
    u = ROOT.RooUnfoldBayes(resp, reco_input, DEFAULT_ITERS)
    u.SetVerbose(0)
    if ERROR_MODE == ROOT.RooUnfold.kCovToy:
        u.SetNToys(NTOYS_FINAL)
    h_unfold = u.Hreco(ERROR_MODE)
    if not h_unfold:
        raise RuntimeError(f"photon RooUnfold returned null for {case.key}")
    h_unfold.SetDirectory(0)
    return h_unfold, {
        "reco_data": f"h_unfoldRecoPho_pTgamma_{PHO_KEY}{case.cent_suffix}",
        "reco_sim": f"h_unfoldRecoPho_pTgamma_{PHO_KEY}{case.cent_suffix}",
        "truth_sim": f"h_unfoldTruthPho_pTgamma{case.cent_suffix}",
        "response": f"h2_unfoldResponsePho_pTgamma_{PHO_KEY}{case.cent_suffix}",
        "response_orientation": orientation,
        "purity_input": purity_meta,
        "iterations": DEFAULT_ITERS,
        "error_mode": ERROR_MODE_NAME,
    }


def unfold_xj(case: Case, data_file: ROOT.TFile, sim_file: ROOT.TFile) -> Tuple[object, Dict]:
    reco_data_a = get_obj(data_file, case.data_topdir, f"h2_unfoldReco_pTgamma_xJ_incl_{BASE_KEY}{case.cent_suffix}", "TH2")
    reco_data_c = get_optional(data_file, case.data_topdir, f"h2_unfoldReco_pTgamma_xJ_incl_sidebandC_{BASE_KEY}{case.cent_suffix}")
    reco_sim = get_obj(sim_file, case.sim_topdir, f"h2_unfoldReco_pTgamma_xJ_incl_{BASE_KEY}{case.cent_suffix}", "TH2")
    truth_sim = get_obj(sim_file, case.sim_topdir, f"h2_unfoldTruth_pTgamma_xJ_incl_{BASE_KEY}{case.cent_suffix}", "TH2")
    rsp_sim = get_obj(sim_file, case.sim_topdir, f"h2_unfoldResponse_pTgamma_xJ_incl_{BASE_KEY}{case.cent_suffix}", "TH2")
    reco_data, purity_meta = apply_xj_abcd_input(case, data_file, sim_file, reco_data_a, reco_data_c)
    comb_meta = {"applied": False}
    if case.apply_combinatoric_subtraction:
        comb = get_optional(sim_file, case.sim_topdir, f"h2_unfoldRecoCombinatoric_pTgamma_xJ_incl_{BASE_KEY}{case.cent_suffix}")
        if comb is not None:
            pho_data_for_scale, _scale_meta = apply_photon_abcd_input(case, data_file, sim_file, get_obj(data_file, case.data_topdir, f"h_unfoldRecoPho_pTgamma_{PHO_KEY}{case.cent_suffix}", "TH1"))
            pho_sim_for_scale = get_obj(sim_file, case.sim_topdir, f"h_unfoldRecoPho_pTgamma_{PHO_KEY}{case.cent_suffix}", "TH1")
            raw_integral = float(comb.Integral())
            scaled_integral = 0.0
            for ix in range(0, comb.GetXaxis().GetNbins() + 2):
                n_data = pho_data_for_scale.GetBinContent(ix)
                n_sim = pho_sim_for_scale.GetBinContent(ix)
                scale = n_data / n_sim if n_sim > 0.0 else 0.0
                row = sum(comb.GetBinContent(ix, iy) for iy in range(0, comb.GetYaxis().GetNbins() + 2))
                scaled_integral += row * scale
                for iy in range(0, comb.GetYaxis().GetNbins() + 2):
                    comb.SetBinContent(ix, iy, comb.GetBinContent(ix, iy) * scale)
                    comb.SetBinError(ix, iy, comb.GetBinError(ix, iy) * scale)
            reco_data.Add(comb, -1.0)
            comb_meta = {
                "applied": True,
                "histogram": f"h2_unfoldRecoCombinatoric_pTgamma_xJ_incl_{BASE_KEY}{case.cent_suffix}",
                "raw_integral": raw_integral,
                "photon_yield_scaled_integral": scaled_integral,
            }
        else:
            comb_meta = {"applied": False, "reason": "missing combinatoric template"}
    reco_sim_glob = flatten_th2_to_global(reco_sim, f"{case.key}_reco_sim_global")
    truth_sim_glob = flatten_th2_to_global(truth_sim, f"{case.key}_truth_sim_global")
    data_glob = flatten_th2_to_global(reco_data, f"{case.key}_data_global")
    rsp = transpose_th2(rsp_sim, f"{case.key}_xj_rsp_recoX_truthY")
    orientation = "transposed_truthX_recoY_to_recoX_truthY"
    resp = ROOT.RooUnfoldResponse(reco_sim_glob, truth_sim_glob, rsp, f"{case.key}_respXJ", f"{case.key}_respXJ")
    u = ROOT.RooUnfoldBayes(resp, data_glob, DEFAULT_ITERS)
    u.SetVerbose(0)
    if ERROR_MODE == ROOT.RooUnfold.kCovToy:
        u.SetNToys(NTOYS_FINAL)
    h_unfold_glob = u.Hreco(ERROR_MODE)
    if not h_unfold_glob:
        raise RuntimeError(f"xJ RooUnfold returned null for {case.key}")
    h_unfold_glob.SetDirectory(0)
    h2_unfold = unflatten_global_to_th2(h_unfold_glob, truth_sim, f"{case.key}_h2_truth_unfolded")
    return h2_unfold, {
        "reco_data_A": f"h2_unfoldReco_pTgamma_xJ_incl_{BASE_KEY}{case.cent_suffix}",
        "reco_data_C": f"h2_unfoldReco_pTgamma_xJ_incl_sidebandC_{BASE_KEY}{case.cent_suffix}",
        "reco_sim": f"h2_unfoldReco_pTgamma_xJ_incl_{BASE_KEY}{case.cent_suffix}",
        "truth_sim": f"h2_unfoldTruth_pTgamma_xJ_incl_{BASE_KEY}{case.cent_suffix}",
        "response": f"h2_unfoldResponse_pTgamma_xJ_incl_{BASE_KEY}{case.cent_suffix}",
        "response_orientation": orientation,
        "purity_input": purity_meta,
        "combinatoric_subtraction": comb_meta,
        "iterations": DEFAULT_ITERS,
        "error_mode": ERROR_MODE_NAME,
    }


def project_per_photon_xj(h2_unfold, h_pho_unfold) -> Dict:
    xaxis = h2_unfold.GetXaxis()
    yaxis = h2_unfold.GetYaxis()
    x_edges = axis_edges(yaxis)
    vals = np.zeros(yaxis.GetNbins())
    err2 = np.zeros(yaxis.GetNbins())
    selected_pt_bins = []
    npho = 0.0
    npho_err2 = 0.0
    for ix in range(1, xaxis.GetNbins() + 1):
        lo = xaxis.GetBinLowEdge(ix)
        hi = xaxis.GetBinUpEdge(ix)
        cen = xaxis.GetBinCenter(ix)
        if not row_in_pt_window(lo, hi, cen):
            continue
        selected_pt_bins.append([float(xaxis.GetBinLowEdge(ix)), float(xaxis.GetBinUpEdge(ix)), float(cen)])
        for iy in range(1, yaxis.GetNbins() + 1):
            vals[iy - 1] += h2_unfold.GetBinContent(ix, iy)
            e = h2_unfold.GetBinError(ix, iy)
            err2[iy - 1] += e * e
    for ib in range(1, h_pho_unfold.GetXaxis().GetNbins() + 1):
        lo = h_pho_unfold.GetXaxis().GetBinLowEdge(ib)
        hi = h_pho_unfold.GetXaxis().GetBinUpEdge(ib)
        cen = h_pho_unfold.GetXaxis().GetBinCenter(ib)
        if not row_in_pt_window(lo, hi, cen):
            continue
        npho += h_pho_unfold.GetBinContent(ib)
        e = h_pho_unfold.GetBinError(ib)
        npho_err2 += e * e
    widths = np.diff(x_edges)
    y = np.zeros_like(vals)
    ey = np.zeros_like(vals)
    for i, (v, e2, w) in enumerate(zip(vals, err2, widths)):
        if npho <= 0.0 or w <= 0.0:
            continue
        y[i] = v / (npho * w)
        var = e2 / (npho * npho * w * w)
        var += (v * v * npho_err2) / (npho ** 4 * w * w)
        ey[i] = math.sqrt(var) if var > 0.0 and math.isfinite(var) else 0.0
    return {
        "x_edges": x_edges,
        "x_centers": 0.5 * (x_edges[:-1] + x_edges[1:]),
        "x_widths": widths,
        "y": y,
        "ey": ey,
        "npho_unfolded": npho,
        "npho_error": math.sqrt(npho_err2) if npho_err2 > 0.0 else 0.0,
        "selected_truth_pt_bins": selected_pt_bins,
    }


def draw_panel(result: Dict, case: Case, out_path: Path, standalone: bool = True):
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
        "mathtext.fontset": "dejavuserif",
        "axes.linewidth": 1.1,
    })
    fig, ax = plt.subplots(figsize=(7.8, 5.6), dpi=220)
    _draw_case_on_axis(ax, result, case)
    ax.text(
        0.03, 0.96,
        f"RooUnfoldBayes, {DEFAULT_ITERS} iterations\nABCD-subtracted input" if case.apply_abcd else f"RooUnfoldBayes, {DEFAULT_ITERS} iterations\nraw region-A input",
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=10.5,
        color="#2b2b2b",
    )
    fig.tight_layout(pad=0.35)
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def _draw_case_on_axis(ax, result: Dict, case: Case):
    x = result["x_centers"]
    y = result["y"]
    ey = result["ey"]
    finite = np.isfinite(y) & np.isfinite(ey)
    ax.errorbar(
        x[finite],
        y[finite],
        yerr=ey[finite],
        fmt=case.marker,
        ms=6.5,
        lw=1.3,
        capsize=2.5,
        color=case.color,
        mfc="white",
        mec=case.color,
        mew=1.5,
        label=case.short_label,
    )
    ax.axhline(0.0, color="#888888", lw=0.8)
    ax.set_title(case.title, fontsize=18, fontweight="bold", pad=8)
    ax.set_xlabel(r"$x_{J\gamma}$", fontsize=16)
    ax.set_ylabel(r"$(1/N_\gamma)\,dN/dx_{J\gamma}$", fontsize=15)
    ax.tick_params(axis="both", labelsize=12, direction="in", top=True, right=True)
    ax.grid(True, which="major", color="#e5e7eb", lw=0.8)
    ax.set_xlim(0.0, 1.2)
    ymax = np.nanmax(y[finite] + ey[finite]) if np.any(finite) else 1.0
    ymin = np.nanmin(y[finite] - ey[finite]) if np.any(finite) else 0.0
    if ymin < 0:
        ax.set_ylim(min(-0.05 * ymax, 1.15 * ymin), 1.28 * ymax)
    else:
        ax.set_ylim(0.0, 1.28 * ymax)
    ax.legend(loc="upper right", frameon=False, fontsize=12, handlelength=1.2)


def draw_slide(results: Dict[str, Dict], out_path: Path):
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
        "mathtext.fontset": "dejavuserif",
        "axes.linewidth": 1.1,
    })
    fig = plt.figure(figsize=(16, 9), dpi=160)
    fig.patch.set_facecolor("white")
    gs = fig.add_gridspec(
        2, 3,
        height_ratios=[0.18, 0.82],
        left=0.055,
        right=0.985,
        top=0.94,
        bottom=0.115,
        wspace=0.24,
        hspace=0.05,
    )
    title_ax = fig.add_subplot(gs[0, :])
    title_ax.axis("off")
    title_ax.text(
        0.0, 0.80,
        r"First-pass unfolded $x_{J\gamma}$ distributions",
        fontsize=31,
        fontweight="bold",
        ha="left",
        va="center",
        color="#111827",
    )
    title_ax.text(
        0.0, 0.31,
        r"$15 < E_T^\gamma < 35$ GeV, R = 0.4, $|\Delta\phi| > 7\pi/8$, sliding isolation; "
        r"default Au+Au BDT vs p+p baseV3E reference",
        fontsize=16.5,
        ha="left",
        va="center",
        color="#374151",
    )
    for idx, case in enumerate(CASES):
        ax = fig.add_subplot(gs[1, idx])
        _draw_case_on_axis(ax, results[case.key], case)
        if idx > 0:
            ax.set_ylabel("")
    fig.text(
        0.055, 0.040,
        f"RooUnfoldBayes, {DEFAULT_ITERS} iterations, {ERROR_MODE_NAME} errors. Inputs use the nominal 7pi/8 back-to-back object family.",
        fontsize=11.0,
        color="#4b5563",
        ha="left",
    )
    fig.text(
        0.055, 0.020,
        "p+p panel is response-unfolded raw region A: current THE42 baseV3E ABCD counters are too sparse/coarse in 15-35 GeV for a final purity-corrected p+p curve.",
        fontsize=10.6,
        color="#4b5563",
        ha="left",
    )
    fig.savefig(out_path, dpi=160)
    plt.close(fig)


def draw_auau020_pp_overlay(results: Dict[str, Dict], out_path: Path):
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
        "mathtext.fontset": "dejavuserif",
        "axes.linewidth": 1.15,
    })
    fig = plt.figure(figsize=(16, 9), dpi=160)
    fig.patch.set_facecolor("white")
    ax = fig.add_axes([0.075, 0.165, 0.64, 0.705])
    note_ax = fig.add_axes([0.745, 0.165, 0.205, 0.705])
    note_ax.axis("off")

    overlay = [
        ("auau_0_20", "Au+Au 0-20% corrected", "#1f77b4", "o", "full"),
        ("pp_basev3e", "p+p baseV3E raw-A diagnostic", "#d62728", "s", "none"),
    ]
    for key, label, color, marker, fill in overlay:
        r = results[key]
        x = r["x_centers"]
        y = r["y"]
        ey = r["ey"]
        finite = np.isfinite(y) & np.isfinite(ey)
        mfc = color if fill == "full" else "white"
        ax.errorbar(
            x[finite],
            y[finite],
            yerr=ey[finite],
            fmt=marker,
            ms=8.0,
            lw=1.55,
            capsize=3.0,
            color=color,
            mfc=mfc,
            mec=color,
            mew=1.6,
            label=label,
        )

    ax.axhspan(-0.01, 0.0, color="#f3f4f6", zorder=0)
    ax.axhline(0.0, color="#6b7280", lw=0.9)
    ax.axvspan(0.0, 0.25, color="#f3f4f6", alpha=0.65, zorder=0)
    ax.text(0.125, 0.695, "jet turn-on\nsensitive", ha="center", va="top", fontsize=11.5, color="#6b7280")
    ax.set_xlim(0.0, 1.45)
    ax.set_ylim(-0.035, 0.72)
    ax.set_xlabel(r"$x_{J\gamma}=p_T^{\mathrm{jet}}/p_T^\gamma$", fontsize=18)
    ax.set_ylabel(r"$(1/N_\gamma)\,dN/dx_{J\gamma}$", fontsize=18)
    ax.tick_params(axis="both", labelsize=14, direction="in", top=True, right=True, length=6)
    ax.grid(True, which="major", color="#e5e7eb", lw=0.85)
    ax.legend(loc="upper right", frameon=False, fontsize=13.0, handlelength=1.4)

    fig.text(
        0.075,
        0.935,
        r"Unfolded $x_{J\gamma}$: Au+Au 0-20% vs p+p",
        fontsize=31,
        fontweight="bold",
        ha="left",
        va="center",
        color="#111827",
    )
    fig.text(
        0.075,
        0.892,
        r"$15 < E_T^\gamma < 35$ GeV, anti-$k_T$ R = 0.4, $|\Delta\phi| > 7\pi/8$, sliding isolation",
        fontsize=16.5,
        ha="left",
        va="center",
        color="#374151",
    )

    y0 = 0.94
    note_ax.text(0.0, y0, "What is fixed", fontsize=17, fontweight="bold", color="#111827", va="top")
    note_ax.text(
        0.0,
        y0 - 0.085,
        "• Au+Au uses the default 14-feature BDT\n"
        "  and the rederived WP80 cut.\n"
        "• Au+Au input is ABCD-purity corrected\n"
        "  with the embedded combinatoric template.\n"
        "• Both curves use the nominal 7π/8\n"
        "  back-to-back object family.",
        fontsize=13.2,
        color="#374151",
        va="top",
        linespacing=1.32,
    )
    note_ax.text(0.0, 0.40, "Current limitation", fontsize=17, fontweight="bold", color="#111827", va="top")
    note_ax.text(
        0.0,
        0.315,
        "The current p+p baseV3E ABCD counters\n"
        "are too sparse/coarse in 15-35 GeV for\n"
        "a defensible final purity-corrected p+p\n"
        "overlay, so p+p is shown as a stable\n"
        "raw-region-A response-unfolded check.",
        fontsize=13.2,
        color="#374151",
        va="top",
        linespacing=1.32,
    )

    fig.text(
        0.075,
        0.060,
        f"RooUnfoldBayes, {DEFAULT_ITERS} iterations, {ERROR_MODE_NAME} errors. This is the cleanest current same-object-family overlay; final p+p purity correction is a separate THE-86 input fix.",
        fontsize=11.5,
        color="#4b5563",
        ha="left",
    )
    fig.savefig(out_path, dpi=160)
    plt.close(fig)


def _shape_scaled(result: Dict, xmin: float = 0.0) -> Tuple[np.ndarray, np.ndarray, float]:
    x = result["x_centers"]
    widths = result["x_widths"]
    mask = np.isfinite(result["y"]) & np.isfinite(result["ey"]) & (x >= xmin)
    area = float(np.sum(result["y"][mask] * widths[mask]))
    if not np.isfinite(area) or abs(area) < 1e-12:
        return np.zeros_like(result["y"]), np.zeros_like(result["ey"]), area
    return result["y"] / area, result["ey"] / abs(area), area


def draw_shape_normalized_overlay(results: Dict[str, Dict], out_path: Path):
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
        "mathtext.fontset": "dejavuserif",
        "axes.linewidth": 1.15,
    })
    fig = plt.figure(figsize=(16, 9), dpi=160)
    fig.patch.set_facecolor("white")
    ax = fig.add_axes([0.080, 0.155, 0.69, 0.715])
    note_ax = fig.add_axes([0.800, 0.165, 0.165, 0.705])
    note_ax.axis("off")

    overlay = [
        ("auau_0_20", "Au+Au 0-20% corrected", "#1f77b4", "o", "full"),
        ("pp_basev3e", "p+p baseV3E raw-A", "#d62728", "s", "none"),
    ]
    ymax = 0.0
    areas = {}
    for key, label, color, marker, fill in overlay:
        r = results[key]
        yshape, eyshape, area = _shape_scaled(r, TAIL_SHAPE_XMIN)
        areas[key] = area
        x = r["x_centers"]
        finite = np.isfinite(yshape) & np.isfinite(eyshape) & (x >= TAIL_SHAPE_XMIN)
        ymax = max(ymax, float(np.nanmax(yshape[finite] + eyshape[finite])) if np.any(finite) else 0.0)
        ax.errorbar(
            x[finite],
            yshape[finite],
            yerr=eyshape[finite],
            fmt=marker,
            ms=8.0,
            lw=1.55,
            capsize=3.0,
            color=color,
            mfc=(color if fill == "full" else "white"),
            mec=color,
            mew=1.6,
            label=label,
        )
    ax.axhline(0.0, color="#6b7280", lw=0.9)
    ax.axvspan(0.0, TAIL_SHAPE_XMIN, color="#f3f4f6", alpha=0.70, zorder=0)
    ax.axvline(TAIL_SHAPE_XMIN, color="#6b7280", lw=1.0, ls="--")
    ax.text(0.355, 2.92, "excluded from\nnormalization", ha="center", va="top", fontsize=11.5, color="#6b7280")
    ax.set_xlim(0.30, 1.45)
    ax.set_ylim(-0.08, max(3.1, 1.20 * ymax))
    ax.set_xlabel(r"$x_{J\gamma}=p_T^{\mathrm{jet}}/p_T^\gamma$", fontsize=18)
    ax.set_ylabel(r"tail-normalized shape: $(1/I_{x_J\geq0.5})\,(1/N_\gamma)\,dN/dx_{J\gamma}$", fontsize=16.2)
    ax.tick_params(axis="both", labelsize=14, direction="in", top=True, right=True, length=6)
    ax.grid(True, which="major", color="#e5e7eb", lw=0.85)
    ax.legend(loc="upper right", frameon=False, fontsize=13.0, handlelength=1.4)

    fig.text(0.080, 0.935, r"Tail-shape $x_{J\gamma}$ comparison: normalize for $x_{J\gamma}\geq0.5$", fontsize=30, fontweight="bold", ha="left", va="center", color="#111827")
    fig.text(0.080, 0.892, r"$15 < E_T^\gamma < 35$ GeV, anti-$k_T$ R = 0.4, $|\Delta\phi| > 7\pi/8$; shape-only diagnostic, not the physics normalization", fontsize=16.2, ha="left", va="center", color="#374151")

    note_ax.text(0.0, 0.94, "Why this view", fontsize=17, fontweight="bold", color="#111827", va="top")
    note_ax.text(
        0.0,
        0.85,
        "• Removes the different tail\n"
        "  per-photon yield integrals.\n"
        "• Tests whether the unfolded\n"
        "  xJ>0.5 shape is similar.\n"
        "• Does not replace the nominal\n"
        "  per-photon observable.",
        fontsize=13.0,
        color="#374151",
        va="top",
        linespacing=1.32,
    )
    note_ax.text(0.0, 0.38, r"$x_J\geq0.5$ areas", fontsize=17, fontweight="bold", color="#111827", va="top")
    note_ax.text(
        0.0,
        0.30,
        f"Au+Au 0-20: {areas['auau_0_20']:.3f}\n"
        f"p+p raw-A: {areas['pp_basev3e']:.3f}",
        fontsize=13.4,
        color="#374151",
        va="top",
        linespacing=1.35,
    )
    fig.text(0.080, 0.060, "Use this to compare tail shape only. The nominal physics plot stays normalized per unfolded photon and keeps the integral difference visible.", fontsize=11.5, color="#4b5563", ha="left")
    fig.savefig(out_path, dpi=160)
    plt.close(fig)


def draw_correction_step_diagnostics(
    correction_results: Dict[str, Dict],
    correction_cases: Sequence[Case],
    correction_meta: Dict[str, Dict],
    out_path: Path,
):
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
        "mathtext.fontset": "dejavuserif",
        "axes.linewidth": 1.15,
    })
    fig = plt.figure(figsize=(16, 9), dpi=160)
    fig.patch.set_facecolor("white")
    ax = fig.add_axes([0.075, 0.165, 0.61, 0.705])
    read_ax = fig.add_axes([0.720, 0.165, 0.245, 0.705])
    read_ax.axis("off")

    ymax = 0.0
    rows = []
    for case in correction_cases:
        r = correction_results[case.key]
        x = r["x_centers"]
        y = r["y"]
        ey = r["ey"]
        finite = np.isfinite(y) & np.isfinite(ey)
        ymax = max(ymax, float(np.nanmax(y[finite] + ey[finite])) if np.any(finite) else 0.0)
        ax.errorbar(
            x[finite],
            y[finite],
            yerr=ey[finite],
            fmt=case.marker,
            ms=7.4,
            lw=1.35,
            capsize=2.7,
            color=case.color,
            mfc=("white" if case.key != "auau_0_20_abcd_comb" else case.color),
            mec=case.color,
            mew=1.55,
            label=case.short_label,
        )
        widths = r["x_widths"]
        area = float(np.sum(y * widths))
        rows.append((case.short_label, correction_meta[case.key]["result"]["npho_unfolded"], area))

    ax.axhline(0.0, color="#6b7280", lw=0.9)
    ax.axvspan(0.0, 0.25, color="#f3f4f6", alpha=0.65, zorder=0)
    ax.set_xlim(0.0, 1.45)
    ax.set_ylim(-0.06, max(0.78, 1.20 * ymax))
    ax.set_xlabel(r"$x_{J\gamma}=p_T^{\mathrm{jet}}/p_T^\gamma$", fontsize=18)
    ax.set_ylabel(r"$(1/N_\gamma)\,dN/dx_{J\gamma}$", fontsize=18)
    ax.tick_params(axis="both", labelsize=14, direction="in", top=True, right=True, length=6)
    ax.grid(True, which="major", color="#e5e7eb", lw=0.85)
    ax.legend(loc="upper right", frameon=False, fontsize=13.0, handlelength=1.4)

    fig.text(0.075, 0.935, r"Au+Au 0-20% correction-step stress test", fontsize=31, fontweight="bold", ha="left", va="center", color="#111827")
    fig.text(0.075, 0.892, r"$15 < E_T^\gamma < 35$ GeV, $|\Delta\phi| > 7\pi/8$; same response, same photon window, different measured-input corrections", fontsize=16.2, ha="left", va="center", color="#374151")

    read_ax.text(0.0, 0.94, "Step readout", fontsize=18, fontweight="bold", color="#111827", va="top")
    y = 0.84
    for label, npho, area in rows:
        read_ax.text(0.0, y, label, fontsize=14.2, fontweight="bold", color="#111827", va="top")
        read_ax.text(0.0, y - 0.055, f"unfolded Nγ = {npho:,.0f}\n∫(1/Nγ)dN/dx = {area:.3f}", fontsize=12.7, color="#374151", va="top", linespacing=1.25)
        y -= 0.19
    comb = correction_meta["auau_0_20_abcd_comb"]["xj_unfolding"]["combinatoric_subtraction"]
    read_ax.text(0.0, 0.26, "Combinatoric template", fontsize=18, fontweight="bold", color="#111827", va="top")
    read_ax.text(
        0.0,
        0.18,
        f"raw SIM integral: {comb.get('raw_integral', 0.0):,.0f}\n"
        f"photon-yield scaled: {comb.get('photon_yield_scaled_integral', 0.0):.1f}\n"
        "subtracted before unfolding",
        fontsize=12.7,
        color="#374151",
        va="top",
        linespacing=1.28,
    )
    fig.text(0.075, 0.060, "This diagnoses the standard correction chain. An ML background remover would be a new estimator and needs closure/systematic validation before replacing this chain.", fontsize=11.5, color="#4b5563", ha="left")
    fig.savefig(out_path, dpi=160)
    plt.close(fig)


def main() -> None:
    manifest = {
        "script": str(Path(__file__).relative_to(REPO)),
        "base_key": BASE_KEY,
        "photon_key": PHO_KEY,
        "pt_window_gev": PT_WINDOW,
        "pt_selection_rule": "truth photon pT bins with center >=15 and <35 GeV",
        "roo_unfold": {
            "method": "RooUnfoldBayes",
            "iterations": DEFAULT_ITERS,
            "covariance": ERROR_MODE_NAME,
            "ntoys": NTOYS_FINAL if ERROR_MODE == ROOT.RooUnfold.kCovToy else 0,
            "toy_covariance_role": "stress-test diagnostic only; plotted first-pass statistical bars use kCovariance because kCovToy is unstable for sparse background-subtracted inputs",
        },
        "outputs": {},
        "cases": {},
    }
    results: Dict[str, Dict] = {}
    for case in CASES:
        data_f = open_root(case.data_file)
        sim_f = open_root(case.sim_file)
        h_pho_unf, pho_meta = unfold_photons(case, data_f, sim_f)
        h2_xj_unf, xj_meta = unfold_xj(case, data_f, sim_f)
        result = project_per_photon_xj(h2_xj_unf, h_pho_unf)
        results[case.key] = result
        panel_path = OUT_DIR / f"the85_unfolded_xjgamma_{case.key}_{OUTPUT_TAG}.png"
        draw_panel(result, case, panel_path)
        npz_path = OUT_DIR / f"the85_unfolded_xjgamma_{case.key}_{OUTPUT_TAG}.npz"
        np.savez(
            npz_path,
            x_edges=result["x_edges"],
            x_centers=result["x_centers"],
            x_widths=result["x_widths"],
            y=result["y"],
            ey=result["ey"],
            npho_unfolded=result["npho_unfolded"],
            npho_error=result["npho_error"],
            selected_truth_pt_bins=np.array(result["selected_truth_pt_bins"], dtype=float),
        )
        manifest["outputs"][case.key] = {
            "panel_png": str(panel_path),
            "npz": str(npz_path),
        }
        result_meta = {k: v for k, v in result.items() if k not in ("x_edges", "x_centers", "x_widths", "y", "ey")}
        manifest["cases"][case.key] = {
            **asdict(case),
            "photon_unfolding": pho_meta,
            "xj_unfolding": xj_meta,
            "result": result_meta,
            "integral_per_photon_dndx": float(np.sum(result["y"] * result["x_widths"])),
        }
        data_f.Close()
        sim_f.Close()
    correction_cases = [
        replace(
            CASES[0],
            key="auau_0_20_rawA",
            title="Au+Au 0-20% raw A",
            short_label="raw region A",
            color="#6b7280",
            marker="o",
            apply_abcd=False,
            apply_combinatoric_subtraction=False,
        ),
        replace(
            CASES[0],
            key="auau_0_20_abcd_only",
            title="Au+Au 0-20% ABCD only",
            short_label="ABCD purity",
            color="#9467bd",
            marker="s",
            apply_abcd=True,
            apply_combinatoric_subtraction=False,
        ),
        replace(
            CASES[0],
            key="auau_0_20_abcd_comb",
            title="Au+Au 0-20% ABCD + comb.",
            short_label="ABCD + comb.",
            color="#1f77b4",
            marker="o",
            apply_abcd=True,
            apply_combinatoric_subtraction=True,
        ),
    ]
    correction_results: Dict[str, Dict] = {}
    correction_meta: Dict[str, Dict] = {}
    for case in correction_cases:
        data_f = open_root(case.data_file)
        sim_f = open_root(case.sim_file)
        h_pho_unf, pho_meta = unfold_photons(case, data_f, sim_f)
        h2_xj_unf, xj_meta = unfold_xj(case, data_f, sim_f)
        result = project_per_photon_xj(h2_xj_unf, h_pho_unf)
        correction_results[case.key] = result
        result_meta = {k: v for k, v in result.items() if k not in ("x_edges", "x_centers", "x_widths", "y", "ey")}
        correction_meta[case.key] = {
            **asdict(case),
            "photon_unfolding": pho_meta,
            "xj_unfolding": xj_meta,
            "result": result_meta,
            "integral_per_photon_dndx": float(np.sum(result["y"] * result["x_widths"])),
        }
        npz_path = OUT_DIR / f"the85_unfolded_xjgamma_{case.key}_{OUTPUT_TAG}_correction_step.npz"
        np.savez(
            npz_path,
            x_edges=result["x_edges"],
            x_centers=result["x_centers"],
            x_widths=result["x_widths"],
            y=result["y"],
            ey=result["ey"],
            npho_unfolded=result["npho_unfolded"],
            npho_error=result["npho_error"],
            selected_truth_pt_bins=np.array(result["selected_truth_pt_bins"], dtype=float),
        )
        manifest["outputs"][case.key] = {"npz": str(npz_path)}
        data_f.Close()
        sim_f.Close()
    slide_path = OUT_DIR / f"slide02_unfolded_xjgamma_1x3_{OUTPUT_TAG}.png"
    draw_slide(results, slide_path)
    overlay_path = OUT_DIR / f"slide03_unfolded_xjgamma_auau020_vs_pp_{OUTPUT_TAG}.png"
    draw_auau020_pp_overlay(results, overlay_path)
    shape_path = OUT_DIR / f"slide04_unfolded_xjgamma_tail_shape_normalized_xjgt0p5_auau020_vs_pp_{OUTPUT_TAG}.png"
    draw_shape_normalized_overlay(results, shape_path)
    correction_path = OUT_DIR / f"slide05_unfolded_xjgamma_auau020_correction_steps_{OUTPUT_TAG}.png"
    draw_correction_step_diagnostics(correction_results, correction_cases, correction_meta, correction_path)
    manifest_path = OUT_DIR / f"slide02_unfolded_xjgamma_1x3_{OUTPUT_TAG}_manifest.json"
    manifest["outputs"]["slide_1x3_png"] = str(slide_path)
    manifest["outputs"]["slide_auau020_pp_overlay_png"] = str(overlay_path)
    manifest["outputs"]["slide_shape_normalized_auau020_pp_png"] = str(shape_path)
    manifest["outputs"]["slide_auau020_correction_steps_png"] = str(correction_path)
    manifest["correction_step_cases"] = correction_meta
    manifest["outputs"]["manifest_json"] = str(manifest_path)
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(json.dumps({"slide": str(slide_path), "manifest": str(manifest_path)}, indent=2))


if __name__ == "__main__":
    main()
