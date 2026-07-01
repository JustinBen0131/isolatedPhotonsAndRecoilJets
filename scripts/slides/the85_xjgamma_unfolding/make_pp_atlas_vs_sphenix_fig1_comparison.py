#!/usr/bin/env python3
"""Side-by-side pp reconstructed-xJ background-subtraction comparison.

Left: ATLAS Figure 1 pp panel cropped from ATLAS_xJgamma.pdf.
Right: current sPHENIX/PPG12-baseV3E pp 15-35 GeV reconstructed-input
diagnostic from the THE-85 first-pass objects.
"""

from __future__ import annotations

import json
import csv
import sys
from pathlib import Path

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import make_atlas_fig1_reco_background_breakdown as fig1  # noqa: E402


OUT_DIR = REPO / "dataOutput/the85_auau_xjgamma_unfolding_push/unfolded_xjgamma_firstpass"
ATLAS_CROP = REPO / "dataOutput/the85_auau_xjgamma_unfolding_push/atlas_reference/atlas_fig1_pp_panel_crop_hi.png"
OUT_PNG = OUT_DIR / f"slide08_pp_atlas_vs_sphenix_fig1_comparison{fig1.OUT_SUFFIX}.png"
OUT_MANIFEST = OUT_DIR / f"slide08_pp_atlas_vs_sphenix_fig1_comparison{fig1.OUT_SUFFIX}_manifest.json"
OUT_SCRIPT = OUT_DIR / f"slide08_pp_atlas_vs_sphenix_fig1_comparison{fig1.OUT_SUFFIX}_speaker_script.md"
PPG12_PURITY_CSV = (
    REPO
    / "dataOutput/ppg12PhotonYield/ppg12_photon_yield_v1_data_20260620/purity_fig29_comparison/ppg12_ratio_diagnostic/ppg12_photon_final_bdt_nom_extract.csv"
)


def load_ppg12_purity() -> list[dict]:
    rows: list[dict] = []
    with PPG12_PURITY_CSV.open() as f:
        for row in csv.DictReader(f):
            if row["kind"] != "GRAPH" or row["name"] != "gpurity_leak" or row.get("source") != "data":
                continue
            x = float(row["x"])
            ex_low = float(row["ex_low"])
            ex_high = float(row["ex_high"])
            rows.append(
                {
                    "lo": x - ex_low,
                    "hi": x + ex_high,
                    "center": x,
                    "purity": float(row["y"]),
                    "err": 0.5 * (float(row["ey_low"]) + float(row["ey_high"])),
                }
            )
    if not rows:
        raise RuntimeError(f"no PPG12 data-source gpurity_leak rows found in {PPG12_PURITY_CSV}")
    return rows


def fit_ppg12_purity(purity_rows: list[dict]) -> dict:
    x = np.array([row["center"] for row in purity_rows], dtype=float)
    y = np.array([row["purity"] for row in purity_rows], dtype=float)
    ey = np.array([max(0.015, min(0.30, row["err"])) for row in purity_rows], dtype=float)
    try:
        popt, _pcov = curve_fit(
            fig1.pade11,
            x,
            y,
            sigma=ey,
            absolute_sigma=True,
            p0=[0.55, 0.012, 0.02],
            bounds=([-2.0, -1.0, -0.09], [2.0, 1.0, 0.2]),
            maxfev=200000,
        )
        fitted = np.clip(fig1.pade11(x, *popt), 0.02, 0.995)
        residual = (y - fitted) / ey
        return {
            "model": "pade11",
            "label": "Padé[1/1] fit to PPG12 leakage-corrected photon purity",
            "parameters_a_b_c": [float(v) for v in popt],
            "chi2": float(np.sum(residual * residual)),
            "ndf": int(len(x) - len(popt)),
            "source_points": [
                {
                    "pt": [row["lo"], row["hi"]],
                    "pt_center": row["center"],
                    "purity": row["purity"],
                    "purity_err": row["err"],
                    "fit_purity": float(v),
                }
                for row, v in zip(purity_rows, fitted)
            ],
        }
    except Exception as exc:
        coeff = np.polyfit(x, y, 1, w=1.0 / ey)
        fitted = np.clip(np.polyval(coeff, x), 0.02, 0.995)
        residual = (y - fitted) / ey
        return {
            "model": "weighted_linear_fallback",
            "label": "weighted linear fallback fit to PPG12 leakage-corrected photon purity",
            "coefficients_slope_intercept": [float(coeff[0]), float(coeff[1])],
            "fallback_reason": str(exc),
            "chi2": float(np.sum(residual * residual)),
            "ndf": int(len(x) - len(coeff)),
            "source_points": [
                {
                    "pt": [row["lo"], row["hi"]],
                    "pt_center": row["center"],
                    "purity": row["purity"],
                    "purity_err": row["err"],
                    "fit_purity": float(v),
                }
                for row, v in zip(purity_rows, fitted)
            ],
        }


def fitted_ppg12_purity(lo: float, hi: float, purity_rows: list[dict], fit_meta: dict) -> tuple[float, float]:
    center = 0.5 * (lo + hi)
    if fit_meta["model"] == "pade11":
        purity = float(fig1.pade11(np.array([center], dtype=float), *fit_meta["parameters_a_b_c"])[0])
    else:
        slope, intercept = fit_meta["coefficients_slope_intercept"]
        purity = float(slope * center + intercept)
    _raw_purity, raw_err = overlap_weighted_purity(lo, hi, purity_rows)
    return min(max(purity, 0.0), 0.995), raw_err


def overlap_weighted_purity(lo: float, hi: float, purity_rows: list[dict]) -> tuple[float, float]:
    num = 0.0
    err2 = 0.0
    den = 0.0
    for row in purity_rows:
        ov = max(0.0, min(hi, row["hi"]) - max(lo, row["lo"]))
        if ov <= 0.0:
            continue
        num += ov * row["purity"]
        err2 += (ov * row["err"]) ** 2
        den += ov
    if den <= 0.0:
        center = 0.5 * (lo + hi)
        row = min(purity_rows, key=lambda r: abs(center - 0.5 * (r["lo"] + r["hi"])))
        return min(max(row["purity"], 0.0), 0.995), row["err"]
    return min(max(num / den, 0.0), 0.995), (err2 ** 0.5) / den


def project_pp_with_ppg12_purity_norm() -> dict:
    """Build an ATLAS-like pp reco xJ input using current xJ shape + PPG12 purity.

    The current THE76 pp ROOT contains the right region-A and region-C xJ shapes,
    but its internal raw ABCD purity is still under review against PPG12 Fig.29.
    For this visual cross-check, normalize the C-region xJ shape so each photon
    pT row has the leakage-corrected purity from Shuhang's PPG12 final BDT
    graph. That gives a clean diagnostic without claiming the current counters
    have closed.
    """

    panel = fig1.PANELS[0]
    f = fig1.u.open_root(str(panel.data_file))
    h_a = fig1.u.get_obj(f, panel.data_topdir, fig1.hist_key("h2_unfoldReco_pTgamma_xJ_incl", panel), "TH2")
    h_c = fig1.u.get_obj(f, panel.data_topdir, fig1.hist_key("h2_unfoldReco_pTgamma_xJ_incl_sidebandC", panel), "TH2")
    xedges = fig1.u.axis_edges(h_a.GetYaxis())
    centers = 0.5 * (xedges[:-1] + xedges[1:])
    ny = h_a.GetYaxis().GetNbins()
    raw = np.zeros(ny)
    raw_err2 = np.zeros(ny)
    side = np.zeros(ny)
    side_err2 = np.zeros(ny)
    purity_rows = load_ppg12_purity()
    purity_fit = fit_ppg12_purity(purity_rows)
    pt_rows = []

    for ix in range(1, h_a.GetXaxis().GetNbins() + 1):
        lo = h_a.GetXaxis().GetBinLowEdge(ix)
        hi = h_a.GetXaxis().GetBinUpEdge(ix)
        cen = h_a.GetXaxis().GetBinCenter(ix)
        if not fig1.u.row_in_pt_window(lo, hi, cen):
            continue
        row_a = sum(h_a.GetBinContent(ix, iy) for iy in range(1, ny + 1))
        row_c = sum(h_c.GetBinContent(ix, iy) for iy in range(1, ny + 1))
        purity, purity_err = fitted_ppg12_purity(lo, hi, purity_rows, purity_fit)
        bg_int = max(0.0, (1.0 - purity) * row_a)
        scale_c = bg_int / row_c if row_c > 0.0 else 0.0
        scale_c_err = (row_a * purity_err / row_c) if row_c > 0.0 else 0.0
        for iy in range(1, ny + 1):
            a = h_a.GetBinContent(ix, iy)
            ea = h_a.GetBinError(ix, iy)
            c = h_c.GetBinContent(ix, iy)
            ec = h_c.GetBinError(ix, iy)
            raw[iy - 1] += a
            raw_err2[iy - 1] += ea * ea
            side_y = scale_c * c
            side[iy - 1] += side_y
            side_err2[iy - 1] += (scale_c * ec) ** 2 + (scale_c_err * c) ** 2
        pt_rows.append(
            {
                "pt": [lo, hi],
                "row_A_xJ_integral": row_a,
                "row_C_xJ_integral": row_c,
                "ppg12_leakage_corrected_purity": purity,
                "ppg12_purity_err": purity_err,
                "ppg12_purity_source": purity_fit["label"],
                "sideband_scale_to_C": scale_c,
            }
        )

    f.Close()
    bkg_sub = raw - side
    return {
        "x_edges": xedges,
        "x_centers": centers,
        "raw": raw,
        "raw_err": np.sqrt(raw_err2),
        "sideband": side,
        "sideband_err": np.sqrt(side_err2),
        "comb": np.zeros_like(raw),
        "bkg_sub": bkg_sub,
        "bkg_sub_err": np.sqrt(raw_err2 + side_err2),
        "pt_rows": pt_rows,
        "integrals": {
            "raw": float(np.sum(raw)),
            "sideband": float(np.sum(side)),
            "comb": 0.0,
            "bkg_sub": float(np.sum(bkg_sub)),
        },
        "purity_fit": purity_fit,
    }


def draw_sphenix_pp(ax) -> dict:
    r = project_pp_with_ppg12_purity_norm()
    xedges = r["x_edges"]
    centers = r["x_centers"]
    widths = np.diff(xedges)
    mask = (centers >= 0.2) & (centers <= 1.85)

    ax.stairs(r["raw"], xedges, color="#8f8f8f", lw=2.35, label="Raw data")
    ax.plot([xedges[0], xedges[-1]], [0.0, 0.0], color="#d62728", lw=2.1, linestyle=(0, (1, 1)), label="Comb. bkg. (pp = 0)")
    ax.stairs(
        r["sideband"],
        xedges,
        color="#1f4cff",
        lw=2.35,
        linestyle=(0, (2, 2)),
        label="ID sideband (PPG12 norm.)",
    )
    ax.errorbar(
        centers[mask],
        r["bkg_sub"][mask],
        xerr=0.5 * widths[mask],
        yerr=r["bkg_sub_err"][mask],
        fmt="o",
        color="black",
        ms=5.5,
        elinewidth=1.15,
        capsize=0,
        label="Bkg.-sub. input",
        zorder=5,
    )
    ax.set_xlim(0.2, 1.85)
    ymax = max(1.0, float(np.nanmax(r["raw"][mask]) * 1.28))
    ax.set_ylim(0.0, ymax)
    ax.set_xlabel(r"Reconstructed $x_{J\gamma}$", fontsize=18)
    ax.set_ylabel("Entries", fontsize=18)
    ax.minorticks_on()
    ax.tick_params(labelsize=14, top=True, right=True, direction="in", length=6)
    ax.tick_params(which="minor", top=True, right=True, direction="in", length=3)
    ax.legend(loc="upper right", frameon=False, fontsize=11.5, handlelength=2.8)
    ax.text(0.05, 0.94, r"$\bf{\it{sPHENIX}}$ Internal", transform=ax.transAxes, ha="left", va="top", fontsize=18)
    ax.text(0.05, 0.84, r"p+p, $15<E_T^\gamma<35$ GeV", transform=ax.transAxes, ha="left", va="top", fontsize=15)
    ax.text(0.05, 0.75, r"baseV3E, $|\Delta\phi|>7\pi/8$", transform=ax.transAxes, ha="left", va="top", fontsize=15)
    ax.text(0.05, 0.67, fig1.jet_pt_label(), transform=ax.transAxes, ha="left", va="top", fontsize=13.5, color="#334155")
    ax.text(0.05, 0.60, r"sideband norm: PPG12 Fig. 29", transform=ax.transAxes, ha="left", va="top", fontsize=12.5, color="#334155")
    ax.text(0.62, 0.46, "p+p", transform=ax.transAxes, ha="left", va="center", fontsize=18, fontstyle="italic")
    return r


def main() -> None:
    if not ATLAS_CROP.exists():
        raise FileNotFoundError(ATLAS_CROP)

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "mathtext.fontset": "dejavusans",
            "axes.linewidth": 1.25,
        }
    )

    title_color = "#121827"
    body_color = "#364153"
    fig = plt.figure(figsize=(16, 9), dpi=160)
    fig.patch.set_facecolor("white")
    fig.text(
        0.055,
        0.94,
        r"pp reconstructed $x_{J\gamma}$ input: ATLAS benchmark vs sPHENIX",
        ha="left",
        va="top",
        fontsize=26,
        fontweight="bold",
        color=title_color,
    )
    fig.text(
        0.055,
        0.885,
        "Same diagnostic role as ATLAS Fig. 1: raw region A, sideband background, and background-subtracted input before unfolding. Kinematics differ.",
        ha="left",
        va="top",
        fontsize=14.5,
        color=body_color,
    )

    ax_img = fig.add_axes([0.065, 0.145, 0.405, 0.675])
    img = mpimg.imread(ATLAS_CROP)
    ax_img.imshow(img)
    ax_img.axis("off")
    ax_img.set_title(r"ATLAS pp, $63.1<p_T^\gamma<79.6$ GeV", fontsize=17, pad=9, color=title_color)

    ax_spx = fig.add_axes([0.54, 0.17, 0.405, 0.63])
    result = draw_sphenix_pp(ax_spx)
    ax_spx.set_title(r"sPHENIX pp, $15<E_T^\gamma<35$ GeV", fontsize=17, pad=10, color=title_color)
    fig.savefig(OUT_PNG, dpi=160)
    plt.close(fig)

    manifest = {
        "slide": str(OUT_PNG),
        "atlas_source_pdf": str(REPO / "usefulDocs/external_references/gamma_jet/ATLAS_xJgamma.pdf"),
        "atlas_crop": str(ATLAS_CROP),
        "sphenix_source": {
            "pp_data": str(fig1.PANELS[0].data_file),
            "pp_sim": str(fig1.PANELS[0].sim_file),
            "selection": f"baseV3E p+p, 15<E_T^gamma<35 GeV, |Delta phi|>7pi/8, {fig1.jet_pt_label()}",
            "base_key": fig1.PANELS[0].base_key,
            "effective_base_key": fig1.effective_base_key(fig1.PANELS[0]),
            "jet_pt_key": fig1.JET_PT_KEY or "nominal_key",
        },
        "sphenix_integrals": result["integrals"],
        "sideband_normalization": {
            "method": "Current THE76 region-C xJ shape normalized pT-row-by-row to PPG12 final-BDT leakage-corrected purity.",
            "ppg12_purity_source": str(PPG12_PURITY_CSV),
            "ppg12_purity_fit": result["purity_fit"],
            "reason": "Current THE76 raw ABCD purity counters are known not to be PPG12-closed; using them directly makes the blue sideband visually and numerically too large.",
            "pt_rows": result["pt_rows"],
        },
        "comb_background": {
            "pp_value": 0.0,
            "reason": "The current RooUnfold macro applies the explicit combinatoric template only for embedded AuAu; p+p has no mixed-event heavy-ion combinatoric template.",
        },
        "caveat": "Side-by-side diagnostic comparison only; photon-pT ranges and experiment conditions differ.",
    }
    OUT_MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n")
    OUT_SCRIPT.write_text(
        "This slide compares the reconstructed-level pp input diagnostic from ATLAS Figure 1 with our current sPHENIX pp first-pass diagnostic. "
        "The purpose is to check whether the raw region-A spectrum, photon-ID sideband subtraction, and background-subtracted points have the same qualitative role before unfolding. "
        "For the sPHENIX panel, the xJ shapes come from the fuller THE76 pp output; the region-C sideband shape is normalized with the PPG12 final-BDT leakage-corrected purity because the current THE76 raw ABCD counters are still under closure review.\n\n"
        "It is not a same-kinematics comparison: ATLAS uses 63.1 to 79.6 GeV photons at 5.02 TeV, while the sPHENIX panel uses our current 15 to 35 GeV pp baseV3E diagnostic. "
        "The main useful feature is whether the sideband background is broad and subleading while the background-subtracted input remains peaked in a physically sensible xJ range. "
        "The red combinatoric line is zero by design for p+p in this pipeline; the explicit combinatoric template is an embedded-AuAu correction.\n"
    )
    print(json.dumps({"slide": str(OUT_PNG), "manifest": str(OUT_MANIFEST), "speaker_script": str(OUT_SCRIPT)}, indent=2))


if __name__ == "__main__":
    main()
