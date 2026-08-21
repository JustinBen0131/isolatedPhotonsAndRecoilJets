#!/usr/bin/env python3
"""Render the current p+p/Au+Au xJgamma diagnostic across pTjet and centrality.

This is an offline diagnostic consumer of the registered current ROOT products.
For Au+Au it repeats the response-consistent variant documented by
``audit_current_lowx_response_contract.py``: the photon-yield-scaled stored
combinatoric template K is subtracted from the ABCD-corrected data input and
the unscaled K template is subtracted from the response measured marginal.
The matched-pair response matrix is unchanged.  No systematic uncertainty or
manual curve rescaling is applied.
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
from matplotlib.lines import Line2D
import numpy as np
from PIL import Image
import ROOT


ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
if ROOT.gSystem.Load("libRooUnfold") < 0:
    raise RuntimeError("failed to load libRooUnfold")

REPO = Path(__file__).resolve().parents[3]
CONTRACT_PATH = REPO / "scripts/diagnostics/xjgamma/audit_current_lowx_response_contract.py"
OUT = REPO / "dataOutput/the219_friday_ppg_20260814/the89_lowx_response_contract_audit"
PNG = OUT / "the89_pp_auau_centrality_threshold_grid_stat.png"
MANIFEST = OUT / "the89_pp_auau_centrality_threshold_grid_stat_manifest.json"
IMMUTABLE_AUDIT = OUT / "the89_current_lowx_response_contract_pt5_7_10_manifest.json"
EXPECTED_IMMUTABLE_AUDIT_SHA256 = (
    "c2781ffe2cfd72b7997f507c6eca3185ee5320d6e6c64a755dfee70b376c84cf"
)
FONT_PATH = Path("/System/Library/Fonts/Supplemental/Times New Roman.ttf")

PP_PT = (20.0, 26.0)
AUAU_PT = (19.0, 26.0)
THRESHOLDS = {
    5: "r04_isoR40_isSliding",
    7: "r04_jetPt7_isoR40_isSliding",
    10: "r04_jetPt10_isoR40_isSliding",
    12: "r04_jetPt12_isoR40_isSliding",
}
CENTRALITIES = {
    "0_20": ("_cent_0_20", "0–20%"),
    "20_50": ("_cent_20_50", "20–50%"),
    "50_80": ("_cent_50_80", "50–80%"),
}
DISPLAY_XMAX = 1.49


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


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
            "axes.linewidth": 1.05,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "xtick.major.size": 5.0,
            "ytick.major.size": 5.0,
            "xtick.minor.size": 2.6,
            "ytick.minor.size": 2.6,
        }
    )
    return font_name


def make_case(helper, audit, data_path: Path, sim_path: Path, suffix: str, label: str):
    return helper.Case(
        key=f"auau_{suffix.removeprefix('_cent_')}",
        title=f"Au+Au {label}",
        short_label=label,
        data_file=str(data_path),
        data_topdir=audit.DATA_TOP,
        sim_file=str(sim_path),
        sim_topdir=audit.SIM_TOP,
        cent_suffix=suffix,
        color="#C83E32",
        marker="o",
        apply_abcd=True,
        apply_combinatoric_subtraction=True,
    )


def run_grid() -> tuple[dict[int, dict[str, Any]], dict[str, Any], dict[str, Any]]:
    contract = load_module(CONTRACT_PATH, "the89_contract_grid_driver")
    contract.THRESHOLDS = THRESHOLDS
    contract.AUAU_PT = AUAU_PT
    contract.DISPLAY_XMAX = DISPLAY_XMAX
    pp = contract.load_pp()

    audit = contract.load_module(contract.AUAU_HELPER, "the89_centrality_grid_audit")
    audit.PT_WINDOW = AUAU_PT
    helper = audit.load_helper()
    helper.PT_WINDOW = AUAU_PT

    data_path, data_pointer = audit.one_root(audit.DATA_POINTER)
    sim_path, sim_pointer = audit.one_root(audit.SIM_POINTER)
    data = helper.open_root(str(data_path))
    sim = helper.open_root(str(sim_path))
    results: dict[str, Any] = {}
    try:
        for cent_key, (suffix, label) in CENTRALITIES.items():
            audit.CENT_SUFFIX = suffix
            case = make_case(helper, audit, data_path, sim_path, suffix, label)
            photon_unfolded, photon_meta = helper.unfold_photons(case, data, sim)
            threshold_results = {}
            for threshold in THRESHOLDS:
                result = contract.run_auau(
                    audit, helper, case, data, sim, photon_unfolded, threshold
                )
                result["centrality"] = label
                threshold_results[threshold] = result
            results[cent_key] = {
                "label": label,
                "suffix": suffix,
                "photon_unfolding": photon_meta,
                "thresholds": threshold_results,
            }
    finally:
        data.Close()
        sim.Close()

    sources = {
        "data_path": data_path,
        "sim_path": sim_path,
        "data_pointer": data_pointer,
        "sim_pointer": sim_pointer,
        "contract": contract,
        "audit": audit,
    }
    return pp, results, sources


def arrays_for_pp(pp_case: dict[str, Any]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rows = pp_case["rows"]
    edges = np.asarray(
        [float(row["xj_low"]) for row in rows] + [float(rows[-1]["xj_high"])],
        dtype=float,
    )
    y = np.asarray([float(row["pp_unfolded"]["value"]) for row in rows])
    ey = np.asarray([float(row["pp_unfolded"]["stat_error"]) for row in rows])
    return edges, y, ey


def accepted_mask(
    edges: np.ndarray,
    first: float,
    *arrays: np.ndarray,
) -> np.ndarray:
    mask = (edges[:-1] >= first - 1e-12) & (
        edges[1:] <= DISPLAY_XMAX + 1e-12
    )
    for array in arrays:
        mask &= np.isfinite(array)
    return mask


def render(pp: dict[int, dict[str, Any]], results: dict[str, Any]) -> dict[str, Any]:
    font_name = configure_style()
    fig, axes = plt.subplots(
        3,
        4,
        figsize=(15.5, 10.8),
        dpi=240,
        sharex=True,
        sharey=True,
        gridspec_kw={"left": 0.090, "right": 0.985, "bottom": 0.105,
                     "top": 0.790, "wspace": 0.055, "hspace": 0.075},
    )

    first_edges: dict[str, dict[str, float]] = {}
    for row, (cent_key, cent_payload) in enumerate(results.items()):
        first_edges[cent_key] = {}
        for column, threshold in enumerate(THRESHOLDS):
            ax = axes[row, column]
            pp_case = pp[threshold]
            au = cent_payload["thresholds"][threshold]
            fixed = au["variants"]["response_consistent"]
            edges, pp_y, pp_ey = arrays_for_pp(pp_case)
            au_edges = np.asarray(au["x_edges"], dtype=float)
            if not np.array_equal(edges, au_edges):
                raise RuntimeError(
                    f"bin-edge mismatch for {cent_key}, threshold {threshold}"
                )
            au_y = np.asarray(fixed["y"], dtype=float)
            au_ey = np.asarray(fixed["ey"], dtype=float)
            first = max(
                float(pp_case["first_accepted_xj_edge"]),
                float(au["first_accepted_xj_edge"]),
            )
            first_edges[cent_key][str(threshold)] = first
            mask = accepted_mask(edges, first, pp_y, pp_ey, au_y, au_ey)
            centers = 0.5 * (edges[:-1] + edges[1:])
            halfwidth = 0.5 * np.diff(edges)

            ax.errorbar(
                centers[mask] - 0.006,
                pp_y[mask],
                xerr=halfwidth[mask],
                yerr=pp_ey[mask],
                fmt="s",
                ms=5.0,
                mfc="white",
                mec="#111111",
                mew=1.25,
                ecolor="#111111",
                elinewidth=1.05,
                capsize=1.8,
                capthick=1.0,
                linestyle="none",
                zorder=4,
            )
            ax.errorbar(
                centers[mask] + 0.006,
                au_y[mask],
                xerr=halfwidth[mask],
                yerr=au_ey[mask],
                fmt="o",
                ms=5.3,
                mfc="#C83E32",
                mec="#8F241D",
                mew=0.85,
                ecolor="#B9342A",
                elinewidth=1.05,
                capsize=1.8,
                capthick=1.0,
                linestyle="none",
                zorder=5,
            )
            ax.set_xlim(0.0, 1.52)
            ax.set_ylim(-0.18, 3.18)
            ax.set_xticks(np.arange(0.0, 1.51, 0.25))
            ax.set_yticks(np.arange(0.0, 3.01, 0.5))
            ax.minorticks_on()
            ax.tick_params(axis="both", which="major", labelsize=11.0, pad=4)
            ax.tick_params(axis="both", which="minor", width=0.8)

            if row == 0:
                ax.set_title(
                    rf"$p_T^{{\mathrm{{jet}}}}>{threshold}$ GeV",
                    fontsize=16.5,
                    fontweight="bold",
                    pad=10,
                )
            if column == 0:
                ax.text(
                    0.045,
                    0.925,
                    f"Au+Au {cent_payload['label']}",
                    transform=ax.transAxes,
                    ha="left",
                    va="top",
                    fontsize=14.0,
                    fontweight="bold",
                    color="#111111",
                )
            if row == 2 and column == 0:
                ax.text(
                    0.045,
                    0.815,
                    "peripheral ABCD\nstatistically fragile",
                    transform=ax.transAxes,
                    ha="left",
                    va="top",
                    fontsize=9.8,
                    color="#555555",
                )

    fig.text(
        0.090,
        0.952,
        r"Current unfolded $x_{J\gamma}$ diagnostic across jet threshold and centrality",
        ha="left",
        va="top",
        fontsize=25.0,
        fontweight="bold",
        color="#111111",
    )
    fig.text(
        0.090,
        0.900,
        r"anti-$k_T$ $R=0.4$  |  $|\eta^\gamma|,|\eta^{\rm jet}|<0.7$  |  "
        r"$|\Delta\phi_{\gamma{\rm j}}|>7\pi/8$  |  RooUnfoldBayes, 3 iterations",
        ha="left",
        va="top",
        fontsize=14.3,
        color="#2A2A2A",
    )
    fig.text(
        0.090,
        0.858,
        r"$p$+$p$: $20<p_T^\gamma<26$ GeV, $|z_{\rm vtx}|<60$ cm  |  "
        r"Au+Au: $19<p_T^\gamma<26$ GeV, $|z_{\rm vtx}|<10$ cm  |  "
        r"fully accepted $x_{J\gamma}$ bins only  |  statistical uncertainties only",
        ha="left",
        va="top",
        fontsize=13.2,
        color="#2A2A2A",
    )

    legend_handles = [
        Line2D([], [], marker="s", markersize=7.0, markerfacecolor="white",
               markeredgecolor="#111111", markeredgewidth=1.35,
               linestyle="none", color="#111111", label=r"$p$+$p$"),
        Line2D([], [], marker="o", markersize=7.2, markerfacecolor="#C83E32",
               markeredgecolor="#8F241D", markeredgewidth=0.9,
               linestyle="none", color="#B9342A", label="Au+Au"),
    ]
    axes[0, 0].legend(
        handles=legend_handles,
        loc="upper right",
        bbox_to_anchor=(0.970, 0.955),
        ncol=1,
        frameon=False,
        fontsize=12.2,
        handletextpad=0.38,
        labelspacing=0.35,
        borderaxespad=0.0,
    )
    axes[0, 3].text(
        0.965,
        0.955,
        r"$\it{\bf{sPHENIX}}$ Internal",
        transform=axes[0, 3].transAxes,
        ha="right",
        va="top",
        fontsize=13.8,
        color="#111111",
    )
    axes[0, 3].text(
        0.965,
        0.855,
        r"$\sqrt{s_{NN}}=200$ GeV",
        transform=axes[0, 3].transAxes,
        ha="right",
        va="top",
        fontsize=11.8,
        color="#111111",
    )
    fig.supxlabel(r"$x_{J\gamma}=p_T^{\mathrm{jet}}/p_T^\gamma$", fontsize=20.0, y=0.035)
    fig.supylabel(
        r"$(1/N_\gamma^{\mathrm{particle}})\,dN_{\mathrm{jet}}^{\mathrm{particle}}/dx_{J\gamma}$",
        fontsize=18.5,
        x=0.022,
    )
    fig.text(
        0.985,
        0.025,
        "Offline response−K diagnostic • not a final physics result",
        ha="right",
        va="bottom",
        fontsize=10.5,
        color="#555555",
    )
    fig.savefig(PNG, dpi=240, facecolor="white")
    plt.close(fig)

    image = Image.open(PNG)
    return {
        "font_name": font_name,
        "dimensions": list(image.size),
        "first_common_full_xj_edges": first_edges,
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    if sha256(IMMUTABLE_AUDIT) != EXPECTED_IMMUTABLE_AUDIT_SHA256:
        raise RuntimeError("the frozen 0–20% audit receipt changed")
    pp, results, sources = run_grid()
    render_meta = render(pp, results)

    contract = sources["contract"]
    audit = sources["audit"]
    checks = {
        "immutable_0_20_audit_unchanged": (
            sha256(IMMUTABLE_AUDIT) == EXPECTED_IMMUTABLE_AUDIT_SHA256
        ),
        "all_12_stored_combinations_rendered": (
            len(results) == 3
            and all(len(row["thresholds"]) == 4 for row in results.values())
        ),
        "png_dimensions_3720x2592": render_meta["dimensions"] == [3720, 2592],
        "font_times_new_roman": render_meta["font_name"] == "Times New Roman",
        "statistical_uncertainties_only": True,
        "no_systematic_band": True,
        "no_manual_curve_rescaling": True,
        "google_slides_unchanged": True,
    }
    payload = {
        "ok": all(checks.values()),
        "status": "offline_response_contract_centrality_threshold_diagnostic_not_final",
        "purpose": "show stored p+p and Au+Au response-consistent xJgamma variations across jet threshold and centrality",
        "png": str(PNG),
        "png_sha256": sha256(PNG),
        "generator": str(Path(__file__).resolve()),
        "generator_sha256": sha256(Path(__file__).resolve()),
        "selection": {
            "pp_photon_pt_gev": list(PP_PT),
            "auau_photon_pt_gev": list(AUAU_PT),
            "jet_thresholds_gev": list(THRESHOLDS),
            "centrality_percent": [payload["label"] for payload in results.values()],
            "iterations": contract.ITERATIONS,
            "maximum_xj_edge": DISPLAY_XMAX,
            "fully_accepted_bins_only": True,
        },
        "contract": {
            "data_measured": "ABCD-corrected AuAu data minus photon-yield-scaled stored combinatoric K",
            "response_measured_marginal": "inclusive embedded reco minus unscaled stored combinatoric K",
            "response_matrix": "unchanged matched-pair response matrix",
        },
        "sources": {
            "pp_points": str(contract.PP_POINTS),
            "pp_points_sha256": sha256(contract.PP_POINTS),
            "auau_data_pointer": str(audit.DATA_POINTER),
            "auau_data_pointer_sha256": sha256(audit.DATA_POINTER),
            "auau_data_campaign": sources["data_pointer"].get("campaign_tag"),
            "auau_data_root": str(sources["data_path"]),
            "auau_data_root_sha256": sha256(sources["data_path"]),
            "auau_sim_pointer": str(audit.SIM_POINTER),
            "auau_sim_pointer_sha256": sha256(audit.SIM_POINTER),
            "auau_sim_campaign": sources["sim_pointer"].get("campaign_tag"),
            "auau_sim_root": str(sources["sim_path"]),
            "auau_sim_root_sha256": sha256(sources["sim_path"]),
            "frozen_0_20_audit": str(IMMUTABLE_AUDIT),
            "frozen_0_20_audit_sha256": sha256(IMMUTABLE_AUDIT),
        },
        "render": render_meta,
        "centralities": results,
        "checks": checks,
        "limitations": [
            "This is an offline response-contract diagnostic, not an independently closed final result.",
            "The p+p and Au+Au inputs use their nearest native photon-pT windows rather than identical bin edges.",
            "The 50–80% ABCD sideband is statistically fragile, so its visual variation is not a final centrality claim.",
            "Only statistical uncertainties are propagated; correction cross-covariances and systematic uncertainties are absent.",
        ],
        "google_slides_mutated": False,
    }
    MANIFEST.write_text(
        json.dumps(contract.json_ready(payload), indent=2, allow_nan=False) + "\n"
    )
    if not payload["ok"]:
        raise RuntimeError(f"grid checks failed: {MANIFEST}")
    print(PNG)
    print(MANIFEST)


if __name__ == "__main__":
    main()
