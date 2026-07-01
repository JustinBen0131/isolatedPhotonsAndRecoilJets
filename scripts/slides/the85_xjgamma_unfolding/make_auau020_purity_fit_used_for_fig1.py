#!/usr/bin/env python3
"""Show the Au+Au 0-20% purity fit used in the THE-85 Fig. 1 analogue."""

from __future__ import annotations

import json
import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


REPO = Path(__file__).resolve().parents[3]
OUT_DIR = REPO / "dataOutput/the85_auau_xjgamma_unfolding_push/unfolded_xjgamma_firstpass"
INPUT_MANIFEST = OUT_DIR / "slide09_atlas_vs_sphenix_fig1_2x2_reco_inputs_strictFullBins_manifest.json"
OUT_PNG = OUT_DIR / "slide10_auau020_purity_fit_used_for_fig1_correction.png"
OUT_MANIFEST = OUT_DIR / "slide10_auau020_purity_fit_used_for_fig1_correction_manifest.json"
OUT_SCRIPT = OUT_DIR / "slide10_auau020_purity_fit_used_for_fig1_correction_speaker_script.md"


def load_fit_payload() -> dict:
    if not INPUT_MANIFEST.exists():
        raise FileNotFoundError(INPUT_MANIFEST)
    manifest = json.loads(INPUT_MANIFEST.read_text())
    auau = manifest["this_analysis"]["auau_0_20"]
    fit = auau["purity_fit"]
    if not fit.get("enabled"):
        raise RuntimeError("input manifest does not contain an enabled AuAu 0-20 purity fit")
    points = fit["input_points"]
    rows = auau["pt_rows"]
    return {"source_manifest": manifest, "auau": auau, "fit": fit, "points": points, "rows": rows}


def load_source_purity_points(fit: dict) -> list[dict]:
    source_csv = fit.get("source_csv")
    source_label = fit.get("source_label")
    if not source_csv or not source_label:
        return []
    path = Path(source_csv)
    if not path.exists():
        return []
    rows: list[dict] = []
    with path.open() as handle:
        for raw in csv.DictReader(handle):
            if raw.get("label") != source_label:
                continue
            parsed = {}
            for key, value in raw.items():
                if key in {"label", "system"}:
                    parsed[key] = value
                    continue
                try:
                    parsed[key] = float(value)
                except ValueError:
                    parsed[key] = math.nan
            rows.append(parsed)
    return rows


def eval_fit(fit: dict, x: np.ndarray) -> np.ndarray:
    if fit.get("fit_model") == "pade11":
        a, b, c = fit["pade11_parameters_a_b_c"]
        return np.clip((a + b * x) / (1.0 + c * x), 0.02, 0.98)
    if "coefficients_slope_intercept" in fit:
        slope, intercept = fit["coefficients_slope_intercept"]
        return np.clip(slope * x + intercept, 0.02, 0.98)
    raise RuntimeError("purity fit manifest has no recognized fitted function")


def fit_readout(fit: dict) -> str:
    if fit.get("fit_model") == "pade11":
        a, b, c = fit["pade11_parameters_a_b_c"]
        return rf"Padé[1/1] fit: $P(E_T)=({a:.3f}{b:+.4f}E_T)/(1{c:+.4f}E_T)$"
    if "coefficients_slope_intercept" in fit:
        slope, intercept = fit["coefficients_slope_intercept"]
        return rf"Fit used in correction: $P(E_T)= {intercept:.3f} {slope:+.4f}\,E_T$"
    return "Fit used in correction: unrecognized function"


def draw() -> None:
    payload = load_fit_payload()
    points = payload["points"]
    rows = payload["rows"]
    fit = payload["fit"]
    source_points = load_source_purity_points(fit) or [
        {
            "pt_lo": p["pt"][0],
            "pt_hi": p["pt"][1],
            "pt_mid": p["pt_center"],
            "raw_purity": p["source_raw_purity"],
            "raw_purity_err": p["source_raw_purity_err"],
            "corrected_purity": p["source_corrected_purity"],
            "corrected_purity_err": p["source_corrected_purity_err"],
        }
        for p in points
    ]

    pt_lo = np.array([p["pt"][0] for p in points], dtype=float)
    pt_hi = np.array([p["pt"][1] for p in points], dtype=float)
    pt = np.array([p["pt_center"] for p in points], dtype=float)
    ex = 0.5 * (pt_hi - pt_lo)
    src_pt_lo = np.array([p["pt_lo"] for p in source_points], dtype=float)
    src_pt_hi = np.array([p["pt_hi"] for p in source_points], dtype=float)
    src_pt = np.array([p["pt_mid"] for p in source_points], dtype=float)
    src_ex = 0.5 * (src_pt_hi - src_pt_lo)
    purity_raw = np.array([p["raw_purity"] for p in source_points], dtype=float)
    purity_raw_err = np.array([p["raw_purity_err"] for p in source_points], dtype=float)
    purity_corr = np.array([p["corrected_purity"] for p in source_points], dtype=float)
    purity_corr_err = np.array([p["corrected_purity_err"] for p in source_points], dtype=float)
    purity_fit = np.array([p["purity_fit"] for p in points], dtype=float)

    beta_fit = np.array([r["sideband_scale"] for r in rows], dtype=float)
    beta_point = []
    for r in rows:
        a = float(r["A"])
        c = float(r["C"])
        f_c = float(r["fC"])
        p_point = float(r.get("standard_corrected_purity") or r["purity"])
        sa_raw = max(0.0, min(a, p_point * a))
        c_bkg = max(0.0, c - f_c * sa_raw)
        beta_point.append(max(0.0, a - sa_raw) / c_bkg if c_bkg > 0.0 else 0.0)
    beta_point = np.array(beta_point, dtype=float)

    plt.rcParams.update(
        {
            "font.family": "Times New Roman",
            "font.serif": ["Times New Roman", "Times", "Nimbus Roman", "DejaVu Serif"],
            "mathtext.fontset": "stix",
            "axes.linewidth": 1.35,
            "xtick.direction": "in",
            "ytick.direction": "in",
        }
    )

    fig = plt.figure(figsize=(16, 9), dpi=160)
    fig.patch.set_facecolor("white")
    title_color = "#111827"
    body_color = "#334155"
    blue = "#1f4cff"
    green = "#0f9f6e"
    orange = "#c96d1b"

    fig.text(
        0.055,
        0.948,
        r"Au+Au 0-20% photon purity used for the Fig. 1-style correction",
        ha="left",
        va="top",
        fontsize=25.5,
        fontweight="bold",
        color=title_color,
    )

    ax1 = fig.add_axes([0.075, 0.185, 0.405, 0.630])
    ax2 = fig.add_axes([0.555, 0.185, 0.375, 0.630])

    xgrid = np.linspace(15.5, 35.0, 300)
    ygrid = eval_fit(fit, xgrid)

    ax1.errorbar(
        src_pt,
        purity_raw,
        xerr=src_ex,
        yerr=purity_raw_err,
        fmt="o",
        ms=8.5,
        color="black",
        mfc="black",
        mec="black",
        mew=1.8,
        elinewidth=1.35,
        capsize=0,
        label="raw ABCD purity",
        zorder=5,
    )
    ax1.errorbar(
        src_pt,
        purity_corr,
        xerr=src_ex,
        yerr=purity_corr_err,
        fmt="D",
        ms=8.5,
        color="#0072b2",
        mfc="white",
        mec="#0072b2",
        mew=1.8,
        elinewidth=1.35,
        capsize=0,
        label="leakage-corrected purity",
        zorder=6,
    )
    if fit.get("fit_model") == "pade11":
        fit_label = "Padé[1/1] fit; scales region C"
    else:
        fit_label = f"{fit.get('fit_label', fit.get('method', 'fit'))}; scales region C"
    ax1.plot(xgrid, ygrid, color=blue, lw=2.0, label=fit_label)
    ax1.set_xlim(14.2, 35.0)
    ax1.set_ylim(0.0, 1.08)
    ax1.set_xlabel(r"photon $E_T$ bin center [GeV]", fontsize=16)
    ax1.set_ylabel(r"prompt-photon purity $P$", fontsize=16)
    ax1.tick_params(labelsize=13.5, top=True, right=True, length=7)
    ax1.tick_params(which="minor", top=True, right=True, length=4)
    ax1.minorticks_on()
    ax1.grid(True, color="#e5e7eb", lw=0.8)
    ax1.legend(loc="upper right", frameon=False, fontsize=12.6, handlelength=2.2)
    ax1.text(
        0.045,
        0.955,
        r"$\bf{\it{sPHENIX}}$ Internal",
        transform=ax1.transAxes,
        ha="left",
        va="top",
        fontsize=15.0,
    )
    ax1.text(
        0.045,
        0.865,
        r"Au+Au 0-20%, $|\Delta\phi|>7\pi/8$",
        transform=ax1.transAxes,
        ha="left",
        va="top",
        fontsize=13.0,
        color=body_color,
    )
    ax1.text(
        0.045,
        0.798,
        r"strict full-bin view: 16-35 GeV effective",
        transform=ax1.transAxes,
        ha="left",
        va="top",
        fontsize=12.2,
        color=body_color,
    )

    ax2.errorbar(
        pt,
        beta_point,
        xerr=ex,
        fmt="s",
        ms=8.0,
        color=orange,
        mfc="white",
        mec=orange,
        mew=1.8,
        elinewidth=1.25,
        capsize=0,
        label=r"using point-by-point leakage-corr. $P$",
    )
    ax2.errorbar(
        pt,
        beta_fit,
        xerr=ex,
        fmt="o",
        ms=8.0,
        color=green,
        mfc=green,
        mec=green,
        mew=1.5,
        elinewidth=1.25,
        capsize=0,
        label=r"using fitted $P$",
    )
    ax2.set_xlim(15.0, 35.0)
    ax2.set_ylim(0.0, max(1.18, 1.18 * np.nanmax([np.nanmax(beta_point), np.nanmax(beta_fit)])))
    ax2.set_xlabel(r"photon $E_T$ bin center [GeV]", fontsize=16)
    ax2.set_ylabel(r"region-C sideband scale $\beta$", fontsize=16)
    ax2.tick_params(labelsize=13.5, top=True, right=True, length=7)
    ax2.tick_params(which="minor", top=True, right=True, length=4)
    ax2.minorticks_on()
    ax2.grid(True, color="#e5e7eb", lw=0.8)
    ax2.legend(loc="upper right", frameon=False, fontsize=12.8, handlelength=2.0)
    ax2.text(
        0.045,
        0.955,
        r"$\beta_i=\frac{(1-P_i)A_i}{C_i-f_C P_i A_i}$",
        transform=ax2.transAxes,
        ha="left",
        va="top",
        fontsize=17.0,
        color=title_color,
    )
    ax2.text(
        0.045,
        0.855,
        r"fit stabilizes the amount of region-C subtraction",
        transform=ax2.transAxes,
        ha="left",
        va="top",
        fontsize=13.0,
        color=body_color,
    )

    readout = f"{fit_readout(fit)}; points are the standard photon-candidate ABCD purities from the 0-20% purity slide."
    fig.text(
        0.075,
        0.090,
        readout,
        ha="left",
        va="center",
        fontsize=15.5,
        color=title_color,
        bbox={
            "boxstyle": "round,pad=0.55,rounding_size=0.08",
            "facecolor": "white",
            "edgecolor": "#cbd5e1",
            "linewidth": 1.2,
        },
    )

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=160)
    plt.close(fig)

    manifest = {
        "slide": str(OUT_PNG),
        "source_manifest": str(INPUT_MANIFEST),
        "source_purity_csv": fit.get("source_csv"),
        "source_purity_label": fit.get("source_label"),
        "fit": fit,
        "beta_rows": [
            {
                "pt": r["pt"],
                "pt_center": float(p),
                "purity_raw": float(pr),
                "purity_corrected": float(pc),
                "purity_fit": float(pf),
                "beta_point_by_point_corrected_purity": float(br),
                "beta_fit": float(bf),
                "A": float(r["A"]),
                "C": float(r["C"]),
                "fC": float(r["fC"]),
            }
            for r, p, pr, pc, pf, br, bf in zip(rows, pt, purity_raw, purity_corr, purity_fit, beta_point, beta_fit)
        ],
        "status": "Companion diagnostic showing the standard photon-candidate ABCD purity fit used in slide09 strict-full-bin Fig.1-style subtraction input.",
        "caveat": "Strict-full-bin input is 16-35 GeV effective until the next pass has a true 15 GeV bin edge.",
    }
    OUT_MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n")
    OUT_SCRIPT.write_text(
        "This slide shows the standard photon-candidate ABCD purity used in the Au+Au 0-20 percent correction on the Fig. 1-style diagnostic. "
        "The left panel shows the raw ABCD purity, the signal-leakage-corrected purity, and the weighted linear fit used for the correction. "
        "The right panel shows how replacing noisy point-by-point leakage-corrected purity with the fitted purity changes the region-C sideband scale beta that controls the dijet or fake-photon subtraction. "
        "The important point is that the fit stabilizes the normalization of the region-C subtraction; it does not change the region-C xJ shape template.\n"
    )
    print(json.dumps({"slide": str(OUT_PNG), "manifest": str(OUT_MANIFEST), "speaker_script": str(OUT_SCRIPT)}, indent=2))


if __name__ == "__main__":
    draw()
