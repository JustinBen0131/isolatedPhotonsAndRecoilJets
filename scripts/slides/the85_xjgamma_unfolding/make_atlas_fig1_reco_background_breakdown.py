#!/usr/bin/env python3
"""ATLAS-Fig.-1-style reconstructed xJ background breakdown for THE-85.

This is a reconstructed-level diagnostic, not an unfolded final result.  It
shows the ingredients that feed the first-pass unfolding input:

  raw region A, the ABCD photon-ID sideband subtraction, combinatoric template,
  and the background-subtracted input.  The Au+Au subtraction matches the
  first-pass Python unfolding helper; pp is shown as the current baseV3E
  diagnostic.
"""

from __future__ import annotations

import json
import math
import os
import sys
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import make_unfolded_xjgamma_1x3 as u  # noqa: E402


OUT_DIR = REPO / "dataOutput/the85_auau_xjgamma_unfolding_push/unfolded_xjgamma_firstpass"
JET_PT_KEY = os.environ.get("THE85_JET_PT_KEY", "").strip()
OUT_SUFFIX = f"_{JET_PT_KEY}" if JET_PT_KEY else ""
OUT_PNG = OUT_DIR / f"slide07_atlas_fig1_style_reco_background_breakdown{OUT_SUFFIX}.png"
OUT_MANIFEST = OUT_DIR / f"slide07_atlas_fig1_style_reco_background_breakdown{OUT_SUFFIX}_manifest.json"
OUT_SCRIPT = OUT_DIR / f"slide07_atlas_fig1_style_reco_background_breakdown{OUT_SUFFIX}_speaker_script.md"
STANDARD_PURITY_CSV = REPO / (
    "dataOutput/the85_auau_xjgamma_unfolding_push/slides/"
    "slide02_purity_leakage_corrected_available_bins_1x3_v3_points.csv"
)
FIT_AUAU020_PURITY = os.environ.get("THE85_FIT_AUAU020_PURITY", "1").strip().lower() not in {
    "0",
    "false",
    "no",
}


@dataclass
class Panel:
    key: str
    label: str
    tag: str
    data_file: Path
    data_topdir: str
    sim_file: Path
    sim_topdir: str
    cent_suffix: str
    apply_comb: bool
    base_key: str = u.BASE_KEY
    pho_key: str = u.PHO_KEY
    use_leakage: bool = True
    y_label: str = "Entries"


PANELS = [
    Panel(
        key="pp",
        label="p+p",
        tag="p+p baseV3E",
        data_file=REPO
        / "InputFiles/pp24/ppg12_photon_yield_v1_data_20260620/pp/RecoilJets_pp_ALL_jetMinPtScan_dphiScan_vz60_isoR40_isSliding_preselectionNewPPG12_tightNewPPG12_nonTightNewPPG12.root",
        data_topdir="Photon_4_GeV_plus_MBD_NS_geq_1",
        sim_file=REPO
        / "InputFiles/pp24/ppg12_photon_yield_v1_signal_sim_ppg12mix_combined_20260625/sim/jetMinPtScan_dphiScan_vz60_isoR40_isSliding_preselectionNewPPG12_tightNewPPG12_nonTightNewPPG12/photonJet5and10and20ppg12mixmerged_SIM/RecoilJets_photonjet5plus10plus20_ppg12mix_MERGED.root",
        sim_topdir="SIM",
        cent_suffix="",
        apply_comb=False,
        base_key="r04",
        pho_key="",
        use_leakage=False,
    ),
    Panel(
        key="auau_50_80",
        label="Au+Au 50-80%",
        tag="Au+Au 50-80%",
        data_file=REPO / "InputFiles/the69_default_auau_physicsqa/RecoilJets_auau_ALL_preselectionNewPPG12_tightAuAuCentInputBase3x3BDT_nonTightAuAuBDTComplement_baseVariant.root",
        data_topdir="photon_12_plus_MBD_NS_geq_2_vtx_lt_150",
        sim_file=REPO / "InputFiles/the69_leakageCentWP_fix/RecoilJets_embeddedPhoton12plus20_MERGED.root",
        sim_topdir="SIM",
        cent_suffix="_cent_50_80",
        apply_comb=True,
    ),
    Panel(
        key="auau_0_20",
        label="Au+Au 0-20%",
        tag="Au+Au 0-20%",
        data_file=REPO / "InputFiles/the69_default_auau_physicsqa/RecoilJets_auau_ALL_preselectionNewPPG12_tightAuAuCentInputBase3x3BDT_nonTightAuAuBDTComplement_baseVariant.root",
        data_topdir="photon_12_plus_MBD_NS_geq_2_vtx_lt_150",
        sim_file=REPO / "InputFiles/the69_leakageCentWP_fix/RecoilJets_embeddedPhoton12plus20_MERGED.root",
        sim_topdir="SIM",
        cent_suffix="_cent_0_20",
        apply_comb=True,
    ),
]


def hist_key(prefix: str, panel: Panel) -> str:
    return f"{prefix}_{effective_base_key(panel)}{panel.cent_suffix}"


def effective_base_key(panel: Panel) -> str:
    if not JET_PT_KEY:
        return panel.base_key
    if f"_{JET_PT_KEY}" in panel.base_key:
        return panel.base_key
    if panel.base_key == "r04":
        return f"r04_{JET_PT_KEY}"
    if panel.base_key.startswith("r04_"):
        return f"r04_{JET_PT_KEY}_{panel.base_key[len('r04_'):]}"
    return f"{panel.base_key}_{JET_PT_KEY}"


def jet_pt_label() -> str:
    if not JET_PT_KEY:
        return r"$p_T^{jet}>5$ GeV"
    if JET_PT_KEY.startswith("jetPt"):
        return rf"$p_T^{{jet}}>{JET_PT_KEY[len('jetPt'):]}$ GeV"
    return JET_PT_KEY


def photon_pt_label() -> str:
    if u.REQUIRE_FULL_PT_BINS:
        return r"full bins in $15<E_T^\gamma<35$ GeV (16-35 GeV effective)"
    return r"$15<E_T^\gamma<35$ GeV"


def photon_hist_key(prefix: str, panel: Panel) -> str:
    if panel.pho_key:
        return f"{prefix}_{panel.pho_key}{panel.cent_suffix}"
    return f"{prefix}{panel.cent_suffix}"


def abcd_suffix(panel: Panel, pt_suffix: str) -> str:
    if panel.pho_key:
        return f"_{panel.pho_key}{pt_suffix}{panel.cent_suffix}"
    return f"{pt_suffix}{panel.cent_suffix}"


def load_leakage_factors(sim_f, panel: Panel) -> Dict[str, Tuple[float, float, float]]:
    if not panel.use_leakage:
        return {suffix: (0.0, 0.0, 0.0) for _lo, _hi, suffix in u.pt_bins_from_edges(u.PT_EDGES_CANON)}

    out: Dict[str, Tuple[float, float, float]] = {}
    for _lo, _hi, suffix in u.pt_bins_from_edges(u.PT_EDGES_CANON):
        if panel.pho_key:
            name = f"h_sigABCD_MC_{panel.pho_key}{suffix}{panel.cent_suffix}"
        else:
            name = f"h_sigABCD_MC{suffix}{panel.cent_suffix}"
        h = sim_f.Get(f"{panel.sim_topdir}/{name}")
        if not h:
            out[suffix] = (0.0, 0.0, 0.0)
            continue
        a = float(h.GetBinContent(1))
        b = float(h.GetBinContent(2))
        c = float(h.GetBinContent(3))
        d = float(h.GetBinContent(4))
        out[suffix] = ((b / a) if a > 0 else 0.0, (c / a) if a > 0 else 0.0, (d / a) if a > 0 else 0.0)
    return out


def use_fitted_purity(panel: Panel) -> bool:
    return FIT_AUAU020_PURITY and panel.key == "auau_0_20"


def pade11(x: np.ndarray, a: float, b: float, c: float) -> np.ndarray:
    return (a + b * x) / (1.0 + c * x)


def load_standard_purity_rows(label: str) -> Dict[Tuple[float, float], Dict]:
    """Load the photon-candidate ABCD purity table used by the purity slide."""
    if not STANDARD_PURITY_CSV.exists():
        raise FileNotFoundError(STANDARD_PURITY_CSV)
    out: Dict[Tuple[float, float], Dict] = {}
    with STANDARD_PURITY_CSV.open() as handle:
        for raw in csv.DictReader(handle):
            if raw.get("label") != label:
                continue
            row: Dict = {}
            for key, value in raw.items():
                if key in {"label", "system"}:
                    row[key] = value
                    continue
                if key in {"solver_ok", "all_abcd_found", "leak_found", "edge_bin_straddles_window"}:
                    row[key] = value == "True"
                    continue
                try:
                    row[key] = float(value)
                except ValueError:
                    row[key] = math.nan
            out[(row["pt_lo"], row["pt_hi"])] = row
    if not out:
        raise RuntimeError(f"no {label} purity rows in {STANDARD_PURITY_CSV}")
    return out


def fitted_purity_values(rows: List[Dict], panel: Panel) -> Tuple[Dict[int, float], Dict]:
    """Fit the standard leakage-corrected photon purity used to scale region C.

    ATLAS estimates photon purity from photon-candidate ABCD counts, then uses
    region C only as the xJ shape template.  Do not fit the event-leading xJ
    diagnostic counters here; those are a QA bookkeeping object, not the photon
    purity measurement.
    """
    label_by_panel = {"auau_0_20": "Au+Au 0-20%"}
    source_label = label_by_panel.get(panel.key)
    if source_label is None:
        return {}, {"enabled": False, "reason": f"no standard purity source configured for {panel.key}"}

    source_by_bin = load_standard_purity_rows(source_label)
    fit_rows: List[Tuple[Dict, Dict]] = []
    for r in rows:
        source = source_by_bin.get((float(r["lo"]), float(r["hi"])))
        if source is None:
            continue
        r["standard_raw_purity"] = source["raw_purity"]
        r["standard_raw_purity_err"] = source["raw_purity_err"]
        r["standard_corrected_purity"] = source["corrected_purity"]
        r["standard_corrected_purity_err"] = source["corrected_purity_err"]
        r["standard_purity_source"] = str(STANDARD_PURITY_CSV)
        if (
            source.get("all_abcd_found", False)
            and source.get("leak_found", False)
            and math.isfinite(source["corrected_purity"])
            and math.isfinite(source["corrected_purity_err"])
            and source["corrected_purity_err"] > 0.0
        ):
            fit_rows.append((r, source))
    if len(fit_rows) < 3:
        return {
            id(r): float(r["standard_corrected_purity"])
            for r, _source in fit_rows
            if math.isfinite(r.get("standard_corrected_purity", math.nan))
        }, {
            "enabled": False,
            "reason": "fewer than 3 valid standard photon-purity points",
            "source_csv": str(STANDARD_PURITY_CSV),
            "source_label": source_label,
        }

    x = np.array([source["pt_mid"] for _r, source in fit_rows], dtype=float)
    y = np.array([source["corrected_purity"] for _r, source in fit_rows], dtype=float)
    ey = np.array([max(0.04, min(0.55, source["corrected_purity_err"])) for _r, source in fit_rows], dtype=float)
    try:
        popt, _pcov = curve_fit(
            pade11,
            x,
            y,
            sigma=ey,
            absolute_sigma=True,
            p0=[0.35, 0.03, 0.02],
            bounds=([-2.0, -1.0, -0.09], [2.0, 1.0, 0.2]),
            maxfev=200000,
        )
        fitted = np.clip(pade11(x, *popt), 0.02, 0.98)
        residual = (y - fitted) / ey
        fit_meta = {
            "fit_model": "pade11",
            "fit_label": "Padé[1/1] fit to leakage-corrected photon purity",
            "pade11_parameters_a_b_c": [float(v) for v in popt],
            "chi2": float(np.sum(residual * residual)),
            "ndf": int(len(x) - len(popt)),
        }
    except Exception as exc:  # pragma: no cover - defensive for malformed input
        coeff = np.polyfit(x, y, 1, w=1.0 / ey)
        fitted = np.clip(np.polyval(coeff, x), 0.02, 0.98)
        fit_meta = {
            "fit_model": "weighted_linear_fallback",
            "fit_label": "weighted linear fallback fit to leakage-corrected photon purity",
            "coefficients_slope_intercept": [float(coeff[0]), float(coeff[1])],
            "fallback_reason": str(exc),
        }
    return (
        {id(r): float(v) for (r, _source), v in zip(fit_rows, fitted)},
        {
            "enabled": True,
            "method": fit_meta["fit_label"],
            "source_csv": str(STANDARD_PURITY_CSV),
            "source_label": source_label,
            "clip_range": [0.02, 0.98],
            **fit_meta,
            "input_points": [
                {
                    "pt": [source["pt_lo"], source["pt_hi"]],
                    "pt_center": source["pt_mid"],
                    "source_raw_purity": source["raw_purity"],
                    "source_raw_purity_err": source["raw_purity_err"],
                    "source_corrected_purity": source["corrected_purity"],
                    "source_corrected_purity_err": source["corrected_purity_err"],
                    "purity_fit": float(v),
                    "source_A": source["A"],
                    "source_C": source["C"],
                    "correction_row_A": r["A"],
                    "correction_row_C": r["C"],
                }
                for (r, source), v in zip(fit_rows, fitted)
            ],
        },
    )


def h2_to_row_arrays(h2) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    xedges = u.axis_edges(h2.GetYaxis())
    ny = h2.GetYaxis().GetNbins()
    vals = np.zeros(ny)
    err2 = np.zeros(ny)
    for ix in range(1, h2.GetXaxis().GetNbins() + 1):
        lo = h2.GetXaxis().GetBinLowEdge(ix)
        hi = h2.GetXaxis().GetBinUpEdge(ix)
        cen = h2.GetXaxis().GetBinCenter(ix)
        if not u.row_in_pt_window(lo, hi, cen):
            continue
        for iy in range(1, ny + 1):
            vals[iy - 1] += h2.GetBinContent(ix, iy)
            err2[iy - 1] += h2.GetBinError(ix, iy) ** 2
    return xedges, vals, np.sqrt(err2)


def build_components(panel: Panel) -> Dict:
    data_f = u.open_root(str(panel.data_file))
    sim_f = u.open_root(str(panel.sim_file))

    h_a = u.get_obj(data_f, panel.data_topdir, hist_key("h2_unfoldReco_pTgamma_xJ_incl", panel), "TH2")
    h_c = u.get_optional(data_f, panel.data_topdir, hist_key("h2_unfoldReco_pTgamma_xJ_incl_sidebandC", panel))
    if h_c is None:
        h_c = h_a.Clone(f"{panel.key}_empty_sideband_c")
        h_c.Reset("ICES")
        h_c.Sumw2()

    xedges, raw, raw_err = h2_to_row_arrays(h_a)
    side = np.zeros_like(raw)
    side_err2 = np.zeros_like(raw)
    side_direct = np.zeros_like(raw)
    side_direct_err2 = np.zeros_like(raw)
    corrected = np.zeros_like(raw)
    corrected_err2 = np.zeros_like(raw)
    pt_rows: List[Dict] = []
    leakage_factors = load_leakage_factors(sim_f, panel)
    reco_bins = u.pt_bins_from_edges(u.PT_EDGES_UNFOLD_RECO)
    ny = h_a.GetYaxis().GetNbins()

    row_work: List[Dict] = []
    for lo, hi, _suffix in reco_bins:
        cen = 0.5 * (lo + hi)
        if not u.row_in_pt_window(lo, hi, cen):
            continue
        canon_suffix = u.canonical_suffix_for_reco_bin(lo, hi)
        leakage = leakage_factors.get(canon_suffix)
        a, b, c, d, sa, esa, nbkg, row_meta = u.compute_abcd_counts(
            data_f,
            panel.data_topdir,
            abcd_suffix(panel, canon_suffix),
            leakage,
        )
        purity = sa / a if a > 0 else 0.0
        fC_leak = max(0.0, leakage[1]) if leakage is not None else 0.0
        row_work.append(
            {
                "lo": lo,
                "hi": hi,
                "pt": [lo, hi],
                "pt_center": cen,
                "ix": h_a.GetXaxis().FindBin(cen),
                "A": a,
                "B": b,
                "C": c,
                "D": d,
                "SA_raw_leakage_corrected": sa,
                "eSA": esa,
                "nbkg_raw": nbkg,
                "purity_raw": purity,
                "fC_leak": fC_leak,
                "leakage": leakage,
                "row_meta": row_meta,
            }
        )

    fit_map, purity_fit_meta = fitted_purity_values(row_work, panel) if use_fitted_purity(panel) else ({}, {"enabled": False})

    for row in row_work:
        lo = row["lo"]
        hi = row["hi"]
        ix = row["ix"]
        a = row["A"]
        b = row["B"]
        c = row["C"]
        d = row["D"]
        fC_leak = row["fC_leak"]
        purity_fit = fit_map.get(id(row))
        purity = purity_fit if purity_fit is not None else row["purity_raw"]
        sa_used = max(0.0, min(a, purity * a))
        nbkg = max(0.0, a - sa_used)
        c_bkg = max(0.0, c - fC_leak * sa_used)
        scale_c = nbkg / c_bkg if c_bkg > 0 else 0.0
        denom_sig_leak = 1.0 - scale_c * fC_leak
        inv_denom = 1.0 / denom_sig_leak if math.isfinite(denom_sig_leak) and abs(denom_sig_leak) > 1e-6 else 1.0
        pooled = u.pooled_sideband_c(h_c, ix, reco_bins, lo, hi, purity)
        for iy in range(1, ny + 1):
            a_val = h_a.GetBinContent(ix, iy)
            a_err = h_a.GetBinError(ix, iy)
            if pooled is not None:
                c_val = float(pooled[0][iy])
                c_err = float(pooled[1][iy])
            else:
                c_val = h_c.GetBinContent(ix, iy)
                c_err = h_c.GetBinError(ix, iy)
            direct_side_val = scale_c * c_val
            direct_side_err = scale_c * c_err
            abcd_val = (a_val - direct_side_val) * inv_denom
            abcd_err = math.sqrt(max(0.0, a_err * a_err + scale_c * scale_c * c_err * c_err)) * abs(inv_denom)
            side_bin = a_val - abcd_val
            side_direct[iy - 1] += direct_side_val
            side_direct_err2[iy - 1] += direct_side_err * direct_side_err
            side[iy - 1] += max(0.0, side_bin)
            side_err2[iy - 1] += (scale_c * c_err * abs(inv_denom)) ** 2
            corrected[iy - 1] += abcd_val
            corrected_err2[iy - 1] += abcd_err * abcd_err
        pt_rows.append({
            "pt": [lo, hi],
            "A": a,
            "B": b,
            "C": c,
            "D": d,
            "purity": purity,
            "purity_source": "standard_photon_abcd_pade11_fit" if purity_fit is not None else "event_leading_xj_abcd_diagnostic",
            "event_leading_xj_purity_diagnostic": row["purity_raw"],
            "standard_raw_purity": row.get("standard_raw_purity"),
            "standard_raw_purity_err": row.get("standard_raw_purity_err"),
            "standard_corrected_purity": row.get("standard_corrected_purity"),
            "standard_corrected_purity_err": row.get("standard_corrected_purity_err"),
            "standard_purity_source": row.get("standard_purity_source"),
            "purity_fit": purity_fit,
            "sideband_scale": scale_c,
            "c_background_after_prompt_leakage": c_bkg,
            "leakage_denominator": denom_sig_leak,
            "leakage_inverse_denominator": inv_denom,
            "pooled_C_shape": pooled is not None,
            "SA_used": sa_used,
            "SA_raw_leakage_corrected": row["SA_raw_leakage_corrected"],
            "eSA": row["eSA"],
            **row["row_meta"],
        })

    comb = np.zeros_like(raw)
    comb_err2 = np.zeros_like(raw)
    comb_meta = {"applied": False}
    if panel.apply_comb:
        h_comb = u.get_optional(sim_f, panel.sim_topdir, hist_key("h2_unfoldRecoCombinatoric_pTgamma_xJ_incl", panel))
        if h_comb is not None:
            h_pho_data = u.get_obj(data_f, panel.data_topdir, photon_hist_key("h_unfoldRecoPho_pTgamma", panel), "TH1")
            h_pho_data, _ = u.apply_photon_abcd_input(
                u.Case(panel.key, panel.label, panel.label, str(panel.data_file), panel.data_topdir, str(panel.sim_file), panel.sim_topdir, panel.cent_suffix, "#000000", "o", True, panel.apply_comb),
                data_f,
                sim_f,
                h_pho_data,
            )
            h_pho_sim = u.get_obj(sim_f, panel.sim_topdir, photon_hist_key("h_unfoldRecoPho_pTgamma", panel), "TH1")
            for ix in range(1, h_comb.GetXaxis().GetNbins() + 1):
                lo = h_comb.GetXaxis().GetBinLowEdge(ix)
                hi = h_comb.GetXaxis().GetBinUpEdge(ix)
                cen = h_comb.GetXaxis().GetBinCenter(ix)
                if not u.row_in_pt_window(lo, hi, cen):
                    continue
                n_data = h_pho_data.GetBinContent(ix)
                n_sim = h_pho_sim.GetBinContent(ix)
                scale = n_data / n_sim if n_sim > 0 else 0.0
                for iy in range(1, h_comb.GetYaxis().GetNbins() + 1):
                    comb[iy - 1] += scale * h_comb.GetBinContent(ix, iy)
                    comb_err2[iy - 1] += (scale * h_comb.GetBinError(ix, iy)) ** 2
            comb_meta = {"applied": True, "integral": float(np.sum(comb))}

    bkg_sub = corrected - comb
    bkg_sub_err = np.sqrt(corrected_err2 + comb_err2)
    data_f.Close()
    sim_f.Close()
    return {
        "x_edges": xedges,
        "x_centers": 0.5 * (xedges[:-1] + xedges[1:]),
        "raw": raw,
        "raw_err": raw_err,
        "sideband": side,
        "sideband_err": np.sqrt(side_err2),
        "sideband_direct": side_direct,
        "sideband_direct_err": np.sqrt(side_direct_err2),
        "comb": comb,
        "comb_err": np.sqrt(comb_err2),
        "bkg_sub": bkg_sub,
        "bkg_sub_err": bkg_sub_err,
        "pt_rows": pt_rows,
        "comb_meta": comb_meta,
        "purity_fit": purity_fit_meta,
        "integrals": {
            "raw": float(np.sum(raw)),
            "id_sideband_net_subtraction": float(np.sum(side)),
            "id_sideband_direct_scaled_C": float(np.sum(side_direct)),
            "comb": float(np.sum(comb)),
            "bkg_sub": float(np.sum(bkg_sub)),
        },
    }


def draw() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "mathtext.fontset": "dejavusans",
            "axes.linewidth": 1.1,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.major.size": 7,
            "ytick.major.size": 7,
            "xtick.minor.size": 4,
            "ytick.minor.size": 4,
        }
    )
    results = {panel.key: build_components(panel) for panel in PANELS}
    fig, axes = plt.subplots(1, 3, figsize=(17.2, 5.8), dpi=180, sharex=True)
    fig.patch.set_facecolor("white")

    for ax, panel in zip(axes, PANELS):
        r = results[panel.key]
        xedges = r["x_edges"]
        centers = r["x_centers"]
        widths = np.diff(xedges)
        mask = (centers >= 0.2) & (centers <= 1.85)

        ax.stairs(r["raw"], xedges, color="#8f8f8f", lw=2.0, label="Raw region A")
        if np.sum(r["comb"]) > 0:
            ax.stairs(r["comb"], xedges, color="#d62728", lw=2.0, linestyle=(0, (1, 1)), label="Comb. bkg.")
        if np.sum(r["sideband"]) > 0:
            ax.stairs(r["sideband"], xedges, color="#1f4cff", lw=2.0, linestyle=(0, (2, 2)), label="ABCD ID-sideband")
        ax.errorbar(
            centers[mask],
            r["bkg_sub"][mask],
            xerr=0.5 * widths[mask],
            yerr=r["bkg_sub_err"][mask],
            fmt="o",
            color="black",
            ms=4.6,
            elinewidth=1.05,
            capsize=0,
            label="Bkg.-sub. input",
            zorder=5,
        )
        ax.set_xlim(0.2, 1.85)
        ymax = max(1.0, float(np.nanmax(r["raw"][mask]) * 1.28))
        ax.set_ylim(0, ymax)
        ax.minorticks_on()
        ax.tick_params(labelsize=12.0, top=True, right=True)
        ax.set_xlabel(r"Reconstructed $x_{J\gamma}$", fontsize=14.5)
        ax.text(0.03, 0.95, r"$\bf{\it{sPHENIX}}$ Internal", transform=ax.transAxes, ha="left", va="top", fontsize=14.5)
        ax.text(0.03, 0.84, r"$15<E_T^\gamma<35$ GeV", transform=ax.transAxes, ha="left", va="top", fontsize=12.5)
        ax.text(0.03, 0.75, r"$|\Delta\phi|>7\pi/8$", transform=ax.transAxes, ha="left", va="top", fontsize=12.5)
        ax.text(0.03, 0.66, jet_pt_label(), transform=ax.transAxes, ha="left", va="top", fontsize=12.5)
        ax.text(0.60, 0.48, panel.label, transform=ax.transAxes, ha="left", va="center", fontsize=14.5, fontweight="bold")
        if panel.key == "pp":
            ax.text(0.03, 0.57, "pp: no comb. template", transform=ax.transAxes, ha="left", va="top", fontsize=10.5, color="#4b5563")

    axes[0].set_ylabel("Entries", fontsize=14.5)
    axes[2].legend(loc="upper right", frameon=False, fontsize=11.8, handlelength=2.6)
    fig.subplots_adjust(left=0.065, right=0.985, top=0.965, bottom=0.16, wspace=0.22)
    fig.savefig(OUT_PNG)
    plt.close(fig)

    manifest = {
        "slide": str(OUT_PNG),
        "source": "Reconstructed-level ATLAS Figure 1 analogue from current THE-85 inputs",
        "pt_window": list(u.PT_WINDOW),
        "require_full_pt_bins": u.REQUIRE_FULL_PT_BINS,
        "jet_pt_key": JET_PT_KEY or "nominal_key",
        "jet_pt_label": jet_pt_label(),
        "dphi": "nominal 7pi/8 object family; pp stores this configured row as r04, AuAu stores r04_isoR40_isSliding",
        "definition": {
            "raw": "region A reconstructed xJ histogram",
            "id_sideband_net_subtraction": "difference between raw region A and the leakage-aware ABCD-corrected signal input; this is the net photon-ID sideband subtraction shown in blue, not a standalone ATLAS dijet template",
            "id_sideband_direct_scaled_C": "literal scale_C * H_C shape recorded for QA; it is not drawn because leakage correction means raw is not decomposed as signal + direct C + comb",
            "combinatoric_background": "embedded combinatoric template scaled by reconstructed photon yield per photon-pT row",
            "bkg_sub_input": "ABCD-corrected input minus scaled combinatoric template, before RooUnfold",
        },
        "panels": {
            panel.key: {
                "label": panel.label,
                "data_file": str(panel.data_file),
                "sim_file": str(panel.sim_file),
                "base_key": panel.base_key,
                "effective_base_key": effective_base_key(panel),
                "photon_key": panel.pho_key or "(configured/default)",
                "leakage_correction": panel.use_leakage,
                "max_fC_leakage": max((row.get("fC", 0.0) for row in results[panel.key]["pt_rows"]), default=0.0),
                **results[panel.key]["integrals"],
                "comb": results[panel.key]["comb_meta"],
                "pt_rows": results[panel.key]["pt_rows"],
            }
            for panel in PANELS
        },
        "caveats": [
            "This is a reconstructed-input diagnostic, not an unfolded particle-level result.",
            "Blue is the ABCD photon-ID sideband subtraction used by this pipeline. It is not the same separately defined ATLAS 'Dijet Bkg.' template.",
            "pp combinatoric background is absent by construction; the pp blue curve is the scaled region-C photon-ID sideband, not a heavy-ion combinatoric template.",
            "AuAu centrality is 0-20 rather than ATLAS 0-10.",
        ],
    }
    OUT_MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n")
    OUT_SCRIPT.write_text(
        "This is the ATLAS Figure-1-style check for our current first-pass inputs. "
        "It is reconstructed xJ before unfolding: grey is raw region A, red is the embedded-AuAu combinatoric template where applicable, blue is the ABCD photon-ID sideband subtraction, and black points are the background-subtracted input that feeds unfolding.\n\n"
        "The key comparison is not the final particle-level shape; it is whether the subtraction ingredients look physically ordered and centrality-dependent in the way ATLAS shows. "
        "For pp, no combinatoric template is drawn because the pp RecoilJets/RooUnfold path has no heavy-ion combinatoric subtraction.\n"
    )
    print(json.dumps({"slide": str(OUT_PNG), "manifest": str(OUT_MANIFEST), "speaker_script": str(OUT_SCRIPT)}, indent=2))


if __name__ == "__main__":
    draw()
