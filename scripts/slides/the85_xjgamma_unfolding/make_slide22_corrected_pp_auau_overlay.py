#!/usr/bin/env python3
"""Slide-22 candidate: corrected pp reference and AuAu response/K fix.

This helper is intentionally scoped to the urgent offline check:

  - pp uses the newer THE76 pp xJ region-A and sideband-C histograms,
    normalized with the PPG12 final-BDT leakage-corrected purity extraction.
  - pp is unfolded through the matching THE76 pp ppg12mix response.
  - AuAu is loaded from the response/K-contract candidate NPZ produced by the
    THE-89 audit.

The pp correction is a candidate, not a final closure claim.  The plotted error
bars are statistical RooUnfold kCovariance bars; the PPG12 purity uncertainty is
tracked in the manifest as a separate normalization/input systematic.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from typing import Dict, Tuple

import matplotlib.pyplot as plt
import numpy as np
import ROOT
from PIL import Image


REPO = Path(__file__).resolve().parents[3]
SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

import make_unfolded_xjgamma_1x3 as u  # noqa: E402
import make_pp_atlas_vs_sphenix_fig1_comparison as ppcomp  # noqa: E402
import make_atlas_vs_sphenix_unfolded_overlay_slide as atlas_overlay  # noqa: E402


OUT_DIR = REPO / "dataOutput/the85_auau_xjgamma_unfolding_push/unfolded_xjgamma_firstpass"
OUT_DIR.mkdir(parents=True, exist_ok=True)

PP_DATA = (
    REPO
    / "InputFiles/pp24/ppg12_photon_yield_v1_data_20260620/pp/"
    "RecoilJets_pp_ALL_jetMinPtScan_dphiScan_vz60_isoR40_isSliding_"
    "preselectionNewPPG12_tightNewPPG12_nonTightNewPPG12.root"
)
PP_SIM = (
    REPO
    / "InputFiles/pp24/ppg12_photon_yield_v1_signal_sim_ppg12mix_combined_20260625/sim/"
    "jetMinPtScan_dphiScan_vz60_isoR40_isSliding_preselectionNewPPG12_tightNewPPG12_nonTightNewPPG12/"
    "photonJet5and10and20ppg12mixmerged_SIM/RecoilJets_photonjet5plus10plus20_ppg12mix_MERGED.root"
)
PP_TOPDIR = "Photon_4_GeV_plus_MBD_NS_geq_1"
SIM_TOPDIR = "SIM"
BASE_KEY = "r04"

PP_NPZ = OUT_DIR / "the85_unfolded_xjgamma_pp_the76_ppg12purity_iter5_covariance_candidate.npz"
AUAU_NPZ = OUT_DIR / "the85_unfolded_xjgamma_auau_0_20_responseKfix_iter5_covariance.npz"
OLD_PP_NPZ = OUT_DIR / "the85_unfolded_xjgamma_pp_basev3e_iter5_covariance.npz"

OUT_PNG = OUT_DIR / "slide22_atlas_vs_sphenix_responseKfix_ppg12purity_candidate_v1.png"
OUT_MANIFEST = OUT_DIR / "slide22_atlas_vs_sphenix_responseKfix_ppg12purity_candidate_v1_manifest.json"
OUT_SCRIPT = OUT_DIR / "slide22_atlas_vs_sphenix_responseKfix_ppg12purity_candidate_v1_speaker_script.md"

XMIN_DISPLAY = 0.20
XMAX_DISPLAY = 1.80


def _load_npz(path: Path) -> Dict[str, np.ndarray]:
    if not path.exists():
        raise FileNotFoundError(path)
    z = np.load(path, allow_pickle=True)
    return {k: z[k] for k in z.files}


def _curve_summary(curve: Dict[str, np.ndarray]) -> Dict[str, float | list[float]]:
    widths = np.diff(curve["x_edges"])
    y = curve["y"]
    x = curve["x_centers"]
    ey = curve["ey"]
    finite = np.isfinite(y)
    peak_i = int(np.nanargmax(np.where(finite, y, -np.inf)))
    return {
        "integral": float(np.nansum(y * widths)),
        "tail_integral_xj_ge_0p4": float(np.nansum(y[x >= 0.4] * widths[x >= 0.4])),
        "tail_integral_xj_ge_0p5": float(np.nansum(y[x >= 0.5] * widths[x >= 0.5])),
        "peak_x": float(x[peak_i]),
        "peak_y": float(y[peak_i]),
        "peak_ey": float(ey[peak_i]),
        "npho_unfolded": float(curve.get("npho_unfolded", np.nan)),
    }


def _fitted_purity(lo: float, hi: float, rows: list[dict], fit_meta: dict) -> Tuple[float, float]:
    purity, purity_err = ppcomp.fitted_ppg12_purity(lo, hi, rows, fit_meta)
    return min(max(float(purity), 0.0), 0.995), max(float(purity_err), 0.0)


def _build_pp_xj_input(data_file: ROOT.TFile, purity_rows: list[dict], fit_meta: dict):
    h_a = u.get_obj(data_file, PP_TOPDIR, f"h2_unfoldReco_pTgamma_xJ_incl_{BASE_KEY}", "TH2")
    h_c = u.get_obj(data_file, PP_TOPDIR, f"h2_unfoldReco_pTgamma_xJ_incl_sidebandC_{BASE_KEY}", "TH2")
    out = h_a.Clone("pp_the76_ppg12purity_h2RecoData_input")
    out.SetDirectory(0)
    out.Reset("ICES")
    out.Sumw2()

    ny = out.GetYaxis().GetNbins()
    rows_meta = []
    for ix in range(1, h_a.GetXaxis().GetNbins() + 1):
        lo = float(h_a.GetXaxis().GetBinLowEdge(ix))
        hi = float(h_a.GetXaxis().GetBinUpEdge(ix))
        cen = float(h_a.GetXaxis().GetBinCenter(ix))
        purity, purity_err = _fitted_purity(lo, hi, purity_rows, fit_meta)
        row_a = sum(h_a.GetBinContent(ix, iy) for iy in range(1, ny + 1))
        row_c = sum(h_c.GetBinContent(ix, iy) for iy in range(1, ny + 1))
        bg_int = max(0.0, (1.0 - purity) * row_a)
        scale_c = bg_int / row_c if row_c > 0.0 else 0.0
        for iy in range(0, ny + 2):
            a = float(h_a.GetBinContent(ix, iy))
            ea = float(h_a.GetBinError(ix, iy))
            c = float(h_c.GetBinContent(ix, iy))
            ec = float(h_c.GetBinError(ix, iy))
            val = a - scale_c * c
            err = math.sqrt(max(0.0, ea * ea + scale_c * scale_c * ec * ec))
            out.SetBinContent(ix, iy, val if math.isfinite(val) else 0.0)
            out.SetBinError(ix, iy, err if math.isfinite(err) else 0.0)
        if row_a > 0.0 or row_c > 0.0:
            rows_meta.append(
                {
                    "pt": [lo, hi],
                    "center": cen,
                    "row_A_xJ_integral": float(row_a),
                    "row_C_xJ_integral": float(row_c),
                    "ppg12_leakage_corrected_purity": float(purity),
                    "ppg12_purity_err_recorded_not_in_stat_bars": float(purity_err),
                    "sideband_scale_to_C": float(scale_c),
                    "signal_xj_integral_after_subtraction": float(row_a - scale_c * row_c),
                }
            )
    return out, rows_meta


def _build_pp_photon_input(data_file: ROOT.TFile, purity_rows: list[dict], fit_meta: dict):
    h_raw = u.get_obj(data_file, PP_TOPDIR, "h_unfoldRecoPho_pTgamma_ppg12obj", "TH1")
    out = h_raw.Clone("pp_the76_ppg12purity_hRecoPho_input")
    out.SetDirectory(0)
    out.Reset("ICES")
    out.Sumw2()
    rows_meta = []
    for ib in range(1, h_raw.GetXaxis().GetNbins() + 1):
        lo = float(h_raw.GetXaxis().GetBinLowEdge(ib))
        hi = float(h_raw.GetXaxis().GetBinUpEdge(ib))
        cen = float(h_raw.GetXaxis().GetBinCenter(ib))
        raw = float(h_raw.GetBinContent(ib))
        eraw = float(h_raw.GetBinError(ib))
        purity, purity_err = _fitted_purity(lo, hi, purity_rows, fit_meta)
        val = raw * purity
        err = eraw * purity
        out.SetBinContent(ib, val if math.isfinite(val) else 0.0)
        out.SetBinError(ib, err if math.isfinite(err) else 0.0)
        if raw > 0.0:
            rows_meta.append(
                {
                    "pt": [lo, hi],
                    "center": cen,
                    "raw_reco_photon_count": raw,
                    "ppg12_leakage_corrected_purity": float(purity),
                    "ppg12_purity_err_recorded_not_in_stat_bars": float(purity_err),
                    "purity_corrected_reco_photon_count": float(val),
                }
            )
    return out, rows_meta


def build_pp_candidate() -> Dict:
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gSystem.Load("libRooUnfold")

    data_file = u.open_root(str(PP_DATA))
    sim_file = u.open_root(str(PP_SIM))
    purity_rows = ppcomp.load_ppg12_purity()
    purity_fit = ppcomp.fit_ppg12_purity(purity_rows)

    reco_data_xj, xj_rows = _build_pp_xj_input(data_file, purity_rows, purity_fit)
    reco_data_pho, pho_rows = _build_pp_photon_input(data_file, purity_rows, purity_fit)

    reco_sim_xj = u.get_obj(sim_file, SIM_TOPDIR, f"h2_unfoldReco_pTgamma_xJ_incl_{BASE_KEY}", "TH2")
    truth_sim_xj = u.get_obj(sim_file, SIM_TOPDIR, f"h2_unfoldTruth_pTgamma_xJ_incl_{BASE_KEY}", "TH2")
    rsp_sim_xj = u.get_obj(sim_file, SIM_TOPDIR, f"h2_unfoldResponse_pTgamma_xJ_incl_{BASE_KEY}", "TH2")

    reco_sim_glob = u.flatten_th2_to_global(reco_sim_xj, "pp_the76_reco_sim_global")
    truth_sim_glob = u.flatten_th2_to_global(truth_sim_xj, "pp_the76_truth_sim_global")
    data_glob = u.flatten_th2_to_global(reco_data_xj, "pp_the76_data_global")
    rsp_xj = u.transpose_th2(rsp_sim_xj, "pp_the76_xj_rsp_recoX_truthY")
    resp_xj = ROOT.RooUnfoldResponse(reco_sim_glob, truth_sim_glob, rsp_xj, "pp_the76_respXJ", "pp_the76_respXJ")
    unf_xj = ROOT.RooUnfoldBayes(resp_xj, data_glob, u.DEFAULT_ITERS)
    unf_xj.SetVerbose(0)
    h_unfold_glob = unf_xj.Hreco(u.ERROR_MODE)
    if not h_unfold_glob:
        raise RuntimeError("pp xJ RooUnfold returned null")
    h_unfold_glob.SetDirectory(0)
    h2_unfold = u.unflatten_global_to_th2(h_unfold_glob, truth_sim_xj, "pp_the76_h2_truth_unfolded")

    reco_sim_pho = u.get_obj(sim_file, SIM_TOPDIR, "h_unfoldRecoPho_pTgamma_ppg12obj", "TH1")
    truth_sim_pho = u.get_obj(sim_file, SIM_TOPDIR, "h_unfoldTruthPho_pTgamma_ppg12obj", "TH1")
    rsp_sim_pho = u.get_obj(sim_file, SIM_TOPDIR, "h2_unfoldResponsePho_pTgamma_ppg12obj", "TH2")
    rsp_pho = u.transpose_th2(rsp_sim_pho, "pp_the76_pho_rsp_recoX_truthY")
    resp_pho = ROOT.RooUnfoldResponse(reco_sim_pho, truth_sim_pho, rsp_pho, "pp_the76_respPho", "pp_the76_respPho")
    unf_pho = ROOT.RooUnfoldBayes(resp_pho, reco_data_pho, u.DEFAULT_ITERS)
    unf_pho.SetVerbose(0)
    h_pho_unfold = unf_pho.Hreco(u.ERROR_MODE)
    if not h_pho_unfold:
        raise RuntimeError("pp photon RooUnfold returned null")
    h_pho_unfold.SetDirectory(0)

    projected = u.project_per_photon_xj(h2_unfold, h_pho_unfold)
    np.savez(
        PP_NPZ,
        x_edges=projected["x_edges"],
        x_centers=projected["x_centers"],
        x_widths=projected["x_widths"],
        y=projected["y"],
        ey=projected["ey"],
        npho_unfolded=projected["npho_unfolded"],
        npho_error=projected["npho_error"],
    )

    data_file.Close()
    sim_file.Close()

    return {
        "npz": str(PP_NPZ),
        "purity_fit": purity_fit,
        "photon_rows": pho_rows,
        "xj_rows": xj_rows,
        "method": {
            "data_region_A": f"{PP_TOPDIR}/h2_unfoldReco_pTgamma_xJ_incl_{BASE_KEY}",
            "data_region_C": f"{PP_TOPDIR}/h2_unfoldReco_pTgamma_xJ_incl_sidebandC_{BASE_KEY}",
            "photon_denominator": f"{PP_TOPDIR}/h_unfoldRecoPho_pTgamma_ppg12obj times PPG12 leakage-corrected purity",
            "xj_numerator": "region A minus sideband C scaled row-by-row to (1-purity)*region-A xJ integral",
            "response": f"{SIM_TOPDIR}/h2_unfoldResponse_pTgamma_xJ_incl_{BASE_KEY}",
            "photon_response": f"{SIM_TOPDIR}/h2_unfoldResponsePho_pTgamma_ppg12obj",
            "iterations": u.DEFAULT_ITERS,
            "error_mode": u.ERROR_MODE_NAME,
            "stat_note": "purity uncertainties are recorded but not folded into plotted statistical bars",
        },
    }


def _draw_curve(ax, curve: Dict[str, np.ndarray], *, label: str, color: str, marker: str, open_marker: bool, zorder: int, alpha: float = 1.0) -> None:
    x = curve["x_centers"]
    edges = curve["x_edges"]
    y = curve["y"]
    ey = curve["ey"]
    xerr = 0.5 * np.diff(edges)
    mask = np.isfinite(y) & np.isfinite(ey) & (x >= XMIN_DISPLAY) & (x <= XMAX_DISPLAY)
    ax.errorbar(
        x[mask],
        y[mask],
        xerr=xerr[mask],
        yerr=ey[mask],
        fmt=marker,
        ms=7.4,
        lw=1.25,
        elinewidth=1.05,
        capsize=2.4,
        color=color,
        mfc="white" if open_marker else color,
        mec=color,
        mew=1.65,
        label=label,
        zorder=zorder,
        alpha=alpha,
    )


def _add_footer_column(fig, x: float, title: str, body: list[str], color: str) -> None:
    fig.text(x, 0.126, title, ha="left", va="top", fontsize=14.2, fontweight="bold", color=color)
    for i, line in enumerate(body):
        fig.text(x, 0.094 - 0.024 * i, line, ha="left", va="top", fontsize=10.8, color="#334155")


def draw_slide(pp_meta: Dict) -> Dict:
    atlas_crop = atlas_overlay.ensure_atlas_panel_crop()
    auau = _load_npz(AUAU_NPZ)
    pp = _load_npz(PP_NPZ)
    old_pp = _load_npz(OLD_PP_NPZ) if OLD_PP_NPZ.exists() else None

    pp_summary = _curve_summary(pp)
    auau_summary = _curve_summary(auau)
    old_pp_summary = _curve_summary(old_pp) if old_pp is not None else None

    y_max = max(
        1.2,
        float(np.nanmax(pp["y"] + pp["ey"])) if len(pp["y"]) else 0.0,
        float(np.nanmax(auau["y"] + auau["ey"])) if len(auau["y"]) else 0.0,
    )
    y_max = min(max(1.45, y_max * 1.18), 2.55)

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "dejavuserif",
            "axes.linewidth": 1.15,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.major.size": 6.5,
            "ytick.major.size": 6.5,
            "xtick.minor.size": 3.5,
            "ytick.minor.size": 3.5,
        }
    )

    title_color = "#121827"
    body_color = "#334155"
    blue = "#244cff"
    red = "#d62728"

    fig = plt.figure(figsize=(16, 9), dpi=170)
    fig.patch.set_facecolor("white")

    fig.text(
        0.045,
        0.94,
        r"Unfolded $x_{J\gamma}$ sanity check: pp quick fix and Au+Au response/K fix",
        ha="left",
        va="top",
        fontsize=29.0,
        fontweight="bold",
        color=title_color,
    )

    atlas_im = Image.open(atlas_crop).convert("RGB")
    ax_atlas = fig.add_axes([0.045, 0.285, 0.395, 0.465])
    ax_atlas.imshow(atlas_im)
    ax_atlas.axis("off")
    fig.text(
        0.055,
        0.855,
        "ATLAS reference: central Pb+Pb vs pp",
        ha="left",
        va="bottom",
        fontsize=18.5,
        fontweight="bold",
        color=body_color,
    )
    fig.text(
        0.055,
        0.820,
        r"Fig. 4, 0-10%; 5.02 TeV; $63.1 < p_T^\gamma < 79.6$ GeV",
        ha="left",
        va="bottom",
        fontsize=13.1,
        color="#526173",
    )
    fig.text(
        0.075,
        0.758,
        "ATLAS",
        ha="left",
        va="bottom",
        fontsize=16.5,
        fontweight="bold",
        color="#111827",
    )
    fig.text(
        0.170,
        0.758,
        "pp: blue open squares    Pb+Pb: red open squares",
        ha="left",
        va="bottom",
        fontsize=12.2,
        color="#334155",
    )

    fig.text(
        0.505,
        0.855,
        r"This analysis: offline sanity-check candidate",
        ha="left",
        va="bottom",
        fontsize=18.5,
        fontweight="bold",
        color=body_color,
    )
    fig.text(
        0.505,
        0.820,
        r"pp: THE76 + PPG12 purity candidate; Au+Au: response measured marginal corrected for K",
        ha="left",
        va="bottom",
        fontsize=12.9,
        color="#526173",
    )
    fig.text(
        0.505,
        0.790,
        r"$p$+$p$ and Au+Au 0-20%, $\sqrt{s_{NN}}=200$ GeV; $15<E_T^\gamma<35$ GeV, $|\Delta\phi|>7\pi/8$",
        ha="left",
        va="bottom",
        fontsize=12.1,
        color="#526173",
    )

    ax = fig.add_axes([0.505, 0.245, 0.435, 0.515])
    if old_pp is not None:
        _draw_curve(
            ax,
            old_pp,
            label=r"old p+p raw-A",
            color="#94a3b8",
            marker="s",
            open_marker=True,
            zorder=2,
            alpha=0.48,
        )
    _draw_curve(
        ax,
        pp,
        label=r"p+p PPG12-purity candidate",
        color=blue,
        marker="s",
        open_marker=True,
        zorder=5,
    )
    _draw_curve(
        ax,
        auau,
        label=r"Au+Au 0-20% response/K fix",
        color=red,
        marker="o",
        open_marker=False,
        zorder=6,
    )
    ax.set_xlim(XMIN_DISPLAY, XMAX_DISPLAY)
    ax.set_ylim(-0.06, y_max)
    ax.set_xlabel(r"$x_{J\gamma}=p_T^{jet}/p_T^\gamma$", fontsize=17.5)
    ax.set_ylabel(r"$(1/N_\gamma)\,dN/dx_{J\gamma}$", fontsize=15.1, labelpad=8)
    ax.grid(True, color="#dfe6ee", lw=0.8)
    ax.minorticks_on()
    ax.tick_params(labelsize=12.9, top=True, right=True)
    ax.axvline(0.4, color="#9ca3af", lw=1.0, ls="--", zorder=1)
    ax.legend(loc="upper right", frameon=False, fontsize=12.2, handlelength=1.8, borderpad=0.1, labelspacing=0.47)
    ax.text(
        0.04,
        0.95,
        r"$\it{\bf{sPHENIX}}$ Internal",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=15.1,
        color="#111827",
    )

    footer = plt.Rectangle((0.045, 0.040), 0.91, 0.125, transform=fig.transFigure, facecolor="#f8fafc", edgecolor="#d7e0ea", lw=1.0)
    fig.add_artist(footer)
    _add_footer_column(
        fig,
        0.065,
        "pp candidate",
        [
            "new THE76 xJ A/C + pp response",
            f"peak {pp_summary['peak_y']:.2f} at xJ={pp_summary['peak_x']:.2f}",
        ],
        blue,
    )
    _add_footer_column(
        fig,
        0.365,
        "Au+Au correction",
        [
            "subtract K from response measured marginal",
            f"tail xJ>=0.5: Au+Au {auau_summary['tail_integral_xj_ge_0p5']:.3f}, pp {pp_summary['tail_integral_xj_ge_0p5']:.3f}",
        ],
        red,
    )
    _add_footer_column(
        fig,
        0.675,
        "Status",
        [
            "no tail renormalization or hand scale",
            "pp scale not solved; closure pending",
        ],
        "#526173",
    )

    fig.savefig(OUT_PNG, dpi=170)
    plt.close(fig)

    manifest = {
        "slide": str(OUT_PNG),
        "pp_candidate_npz": str(PP_NPZ),
        "auau_responseK_npz": str(AUAU_NPZ),
        "old_pp_npz": str(OLD_PP_NPZ) if OLD_PP_NPZ.exists() else None,
        "atlas_reference": {
            "source_pdf": str(REPO / "usefulDocs/external_references/gamma_jet/ATLAS_xJgamma.pdf"),
            "published_reference": "ATLAS Phys. Lett. B 789 (2019), Fig. 4",
            "crop": str(atlas_crop),
            "comparison_role": "visual style/context reference, not same energy or photon-pT range",
        },
        "pp_method": pp_meta["method"],
        "pp_purity_fit": pp_meta["purity_fit"],
        "pp_photon_rows": pp_meta["photon_rows"],
        "pp_xj_rows": pp_meta["xj_rows"],
        "summaries": {
            "pp_candidate": pp_summary,
            "auau_responseK": auau_summary,
            "old_pp_rawA": old_pp_summary,
        },
        "caveats": [
            "No hand scaling or tail-shape renormalization was applied.",
            "The quick pp candidate does not produce the expected ATLAS-like pp peak height; it remains a sanity check, not a fixed final reference.",
            "The pp candidate uses PPG12 final-BDT leakage-corrected purity rather than the current raw THE76 ABCD counters.",
            "PPG12 purity uncertainty is not included in the plotted statistical bars; it is a separate input/systematic uncertainty.",
            "ATLAS panel is context only; energy and photon-pT selections differ.",
        ],
    }
    OUT_MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n")
    OUT_SCRIPT.write_text(
        "This is the fastest defensible offline sanity check for slide 22. The pp reference is no longer the old raw-A THE42 curve. "
        "It uses the newer THE76 pp region-A and sideband-C xJ histograms, scales the sideband with the PPG12 final-BDT leakage-corrected purity, "
        "unfolds through the matching pp ppg12mix response, and divides by a photon denominator corrected with the same PPG12 purity model.\n\n"
        "The AuAu points use the response/K fix from the combinatoric audit: the measured marginal in the response is made consistent with the data input where the embedded K template has already been subtracted. "
        "The important claim is not that this is the final closure-ready plot. It is that the correction is physics-motivated and offline-reproducible, with no tail renormalization or manual scale factor.\n\n"
        "The important sanity-check result is negative for pp: this consistent quick correction still does not raise the pp peak to the ATLAS-like 1.2 scale. "
        "The remaining caveat for the meeting is that pp closure still needs the current ABCD/purity counters and photon-normalization contract repaired or validated. "
        "The PPG12 purity uncertainty is recorded as an input uncertainty and is not included in these plotted statistical bars.\n"
    )

    return {
        "slide": str(OUT_PNG),
        "manifest": str(OUT_MANIFEST),
        "speaker_script": str(OUT_SCRIPT),
        "summaries": manifest["summaries"],
    }


def main() -> None:
    pp_meta = build_pp_candidate()
    result = draw_slide(pp_meta)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
