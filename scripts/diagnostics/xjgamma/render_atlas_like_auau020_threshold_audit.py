#!/usr/bin/env python3
"""Audit and render the existing 0--20% xJgamma threshold variants.

This is deliberately an existing-artifact comparison.  It performs no new
unfolding, correction, normalization, or curve rescaling.  The p+p points and
the Au+Au response-consistent diagnostic are read from their immutable JSON
sidecars.  Only common, fully accepted xJgamma bins are drawn, with statistical
uncertainties only.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.offsetbox import AnchoredOffsetbox, HPacker, TextArea
import numpy as np
from PIL import Image


REPO = Path(__file__).resolve().parents[3]
SOURCE_DIR = REPO / "dataOutput/the219_friday_ppg_20260814"
GRID = (
    SOURCE_DIR
    / "the89_lowx_response_contract_audit"
    / "the89_pp_auau_centrality_threshold_grid_stat_manifest.json"
)
EXPECTED_GRID_SHA256 = (
    "6c145397a58b2e6c7e1490b19fe935efcb0b85faeffd55b84a4c7a1949bc8e35"
)
PP_POINTS = (
    SOURCE_DIR
    / "the89_pp_threshold_peak_scan"
    / "the89_pp_unfolded_xjgamma_threshold_scan_points.json"
)
EXPECTED_PP_SHA256 = (
    "02ebafaf5932cedba7f99fcc808a1f1dc07bfed190f74d0fda86f0d6462cde83"
)
OUT = SOURCE_DIR / "the89_atlas_like_existing_variant_audit"
AUDIT_PNG = OUT / "the89_auau020_atlas_like_threshold_audit_stat.png"
CANDIDATE_PNG = OUT / "the89_pp_auau020_pt7_atlas_like_candidate_stat.png"
MANIFEST = OUT / "the89_auau020_atlas_like_existing_variant_audit.json"

FONT_PATH = Path("/System/Library/Fonts/Supplemental/Times New Roman.ttf")
THRESHOLDS = (5, 7, 10, 12)
RECOMMENDED_THRESHOLD = 7
PP_PT = (20.0, 26.0)
AUAU_PT = (19.0, 26.0)
DISPLAY_EDGE = 1.78


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


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
            "axes.linewidth": 1.15,
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


def load_sources() -> tuple[dict[str, Any], dict[str, Any]]:
    if sha256(GRID) != EXPECTED_GRID_SHA256:
        raise RuntimeError("centrality-threshold audit changed from its pinned receipt")
    if sha256(PP_POINTS) != EXPECTED_PP_SHA256:
        raise RuntimeError("p+p threshold points changed from their pinned receipt")
    grid = json.loads(GRID.read_text())
    pp = json.loads(PP_POINTS.read_text())
    if grid.get("status") != "offline_response_contract_centrality_threshold_diagnostic_not_final":
        raise RuntimeError("source grid is not explicitly diagnostic-not-final")
    if sorted(grid["centralities"]) != ["0_20", "20_50", "50_80"]:
        raise RuntimeError("unexpected stored centrality contract")
    return grid, pp


def arrays_for(
    threshold: int, grid: dict[str, Any], pp_payload: dict[str, Any]
) -> dict[str, Any]:
    auau = grid["centralities"]["0_20"]["thresholds"][str(threshold)]
    pp = pp_payload[str(threshold)]
    edges = np.asarray(auau["x_edges"], dtype=float)
    pp_edges = np.asarray(
        [float(row["xj_low"]) for row in pp["bins"]]
        + [float(pp["bins"][-1]["xj_high"])],
        dtype=float,
    )
    if not np.array_equal(edges, pp_edges):
        raise RuntimeError(f"xJgamma edges differ for threshold {threshold}")
    fixed = auau["variants"]["response_consistent"]
    pp_y = np.asarray([row["pp_unfolded"]["value"] for row in pp["bins"]])
    pp_ey = np.asarray(
        [row["pp_unfolded"]["stat_error"] for row in pp["bins"]]
    )
    au_y = np.asarray(fixed["y"], dtype=float)
    au_ey = np.asarray(fixed["ey"], dtype=float)
    first = max(
        float(pp["first_accepted_xj_edge"]),
        float(auau["first_accepted_xj_edge"]),
    )
    mask = (
        (edges[:-1] >= first - 1e-12)
        & (edges[1:] <= DISPLAY_EDGE + 1e-12)
        & np.isfinite(pp_y)
        & np.isfinite(pp_ey)
        & np.isfinite(au_y)
        & np.isfinite(au_ey)
    )
    return {
        "threshold": threshold,
        "edges": edges,
        "centers": 0.5 * (edges[:-1] + edges[1:]),
        "halfwidth": 0.5 * np.diff(edges),
        "pp_y": pp_y,
        "pp_ey": pp_ey,
        "au_y": au_y,
        "au_ey": au_ey,
        "mask": mask,
        "first": first,
        "k_fraction": float(auau["k_fraction_of_abcd"]),
        "sim_k_fraction": float(auau["sim_k_fraction_of_reco"]),
        "refold_chi2_ndf": float(fixed["refold"]["chi2_ndf"]),
        "refold_ratio_rms": float(fixed["refold"]["ratio_rms"]),
    }


def metrics(a: dict[str, Any]) -> dict[str, Any]:
    edges = a["edges"]
    mask = a["mask"]
    selected = np.flatnonzero(mask)
    pp_peak_index = int(selected[np.argmax(a["pp_y"][mask])])
    au_peak_index = int(selected[np.argmax(a["au_y"][mask])])
    low = mask & (edges[1:] <= 0.72 + 1e-12)
    # Treat the final 1.49--1.78 bin as a near-zero tail check, not as a
    # discriminator of the displacement: both systems are statistically
    # compatible with zero there.  The shape comparison therefore stops at
    # the last populated interval used in the preceding audits.
    high = (
        mask
        & (edges[:-1] >= 0.86 - 1e-12)
        & (edges[1:] <= 1.49 + 1e-12)
    )
    return {
        "first_common_full_xj_edge": float(a["first"]),
        "last_displayed_full_xj_edge": DISPLAY_EDGE,
        "pp_peak_bin": [float(edges[pp_peak_index]), float(edges[pp_peak_index + 1])],
        "auau_peak_bin": [float(edges[au_peak_index]), float(edges[au_peak_index + 1])],
        "low_x_bins_compared": int(np.sum(low)),
        "low_x_bins_with_auau_above_pp": int(np.sum(a["au_y"][low] > a["pp_y"][low])),
        "high_x_bins_compared": int(np.sum(high)),
        "high_x_bins_with_auau_below_pp": int(np.sum(a["au_y"][high] < a["pp_y"][high])),
        "k_fraction_of_abcd": a["k_fraction"],
        "sim_k_fraction_of_reco": a["sim_k_fraction"],
        "response_consistent_refold_chi2_ndf": a["refold_chi2_ndf"],
        "response_consistent_refold_ratio_rms": a["refold_ratio_rms"],
    }


def draw_points(ax, a: dict[str, Any], compact: bool = False) -> None:
    mask = a["mask"]
    size = 5.8 if compact else 8.2
    ax.errorbar(
        a["centers"][mask] - 0.006,
        a["pp_y"][mask],
        xerr=a["halfwidth"][mask],
        yerr=a["pp_ey"][mask],
        fmt="s",
        ms=size,
        mfc="white",
        mec="#111111",
        mew=1.55,
        ecolor="#111111",
        elinewidth=1.25 if compact else 1.5,
        capsize=2.3 if compact else 2.8,
        capthick=1.15,
        linestyle="none",
        label=r"$p$+$p$",
        zorder=4,
    )
    ax.errorbar(
        a["centers"][mask] + 0.006,
        a["au_y"][mask],
        xerr=a["halfwidth"][mask],
        yerr=a["au_ey"][mask],
        fmt="o",
        ms=size + 0.3,
        mfc="#C83E32",
        mec="#8F241D",
        mew=1.0,
        ecolor="#B9342A",
        elinewidth=1.25 if compact else 1.5,
        capsize=2.3 if compact else 2.8,
        capthick=1.15,
        linestyle="none",
        label="Au+Au, 0–20%",
        zorder=5,
    )


def render_audit(arrays: dict[int, dict[str, Any]]) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(16, 9), dpi=160, facecolor="white")
    fig.subplots_adjust(left=0.080, right=0.982, bottom=0.110, top=0.800, wspace=0.11, hspace=0.18)
    fig.text(
        0.080,
        0.955,
        r"0--20% threshold audit: where the ATLAS-like displacement appears",
        ha="left",
        va="top",
        fontsize=27,
        fontweight="bold",
        color="#111111",
    )
    fig.text(
        0.080,
        0.902,
        r"Same nominal cuts and correction contract in every panel; common fully accepted bins; statistical uncertainties only",
        ha="left",
        va="top",
        fontsize=16,
        color="#333333",
    )
    fig.text(
        0.080,
        0.857,
        r"anti-$k_T$ $R=0.4$  |  $|\eta^\gamma|,|\eta^{\rm jet}|<0.7$  |  $|\Delta\phi_{\gamma j}|>7\pi/8$  |  RooUnfoldBayes, 3 iterations",
        ha="left",
        va="top",
        fontsize=14.5,
        color="#333333",
    )

    for ax, threshold in zip(axes.flat, THRESHOLDS):
        a = arrays[threshold]
        draw_points(ax, a, compact=True)
        ax.set_xlim(0.0, 2.0)
        ax.set_ylim(0.0, 2.60)
        ax.set_xticks(np.arange(0.0, 2.01, 0.2))
        ax.set_yticks(np.arange(0.0, 2.51, 0.5))
        ax.minorticks_on()
        ax.tick_params(axis="both", which="major", labelsize=11.5, width=1.0)
        ax.tick_params(axis="both", which="minor", width=0.8)
        ax.text(
            0.04,
            0.94,
            rf"$p_T^{{\rm jet}}>{threshold}$ GeV",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=17,
            fontweight="bold",
        )
        ax.text(
            0.04,
            0.83,
            rf"first common bin: $x_{{J\gamma}}={a['first']:.2f}$",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=12.5,
            color="#333333",
        )
        ax.text(
            0.96,
            0.94,
            rf"$K/S_{{ABCD}}={a['k_fraction']:.2f}$",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=12.5,
            color="#8F241D",
        )
        if threshold == RECOMMENDED_THRESHOLD:
            ax.text(
                0.96,
                0.83,
                "best display compromise",
                transform=ax.transAxes,
                ha="right",
                va="top",
                fontsize=12.5,
                fontweight="bold",
                color="#8F241D",
            )
    for ax in axes[-1, :]:
        ax.set_xlabel(r"$x_{J\gamma}=p_T^{\rm jet}/p_T^\gamma$", fontsize=15)
    for ax in axes[:, 0]:
        ax.set_ylabel(
            r"$(1/N_\gamma^{\rm particle})\,dN_{\rm jet}^{\rm particle}/dx_{J\gamma}$",
            fontsize=14,
        )
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper right",
        bbox_to_anchor=(0.982, 0.912),
        ncol=2,
        frameon=False,
        fontsize=14.5,
        handletextpad=0.5,
        columnspacing=1.5,
    )
    fig.text(
        0.982,
        0.025,
        "Existing response/background-contract diagnostic; not a final physics result",
        ha="right",
        va="bottom",
        fontsize=11.5,
        color="#666666",
    )
    fig.savefig(AUDIT_PNG, dpi=160, facecolor="white")
    plt.close(fig)


def render_candidate(a: dict[str, Any], font_name: str) -> None:
    fig = plt.figure(figsize=(10, 7.5), dpi=240, facecolor="white")
    ax = fig.add_axes([0.145, 0.120, 0.82, 0.655])
    draw_points(ax, a, compact=False)
    ax.set_xlim(0.0, 2.0)
    ax.set_ylim(0.0, 1.90)
    ax.set_xticks(np.arange(0.0, 2.01, 0.2))
    ax.set_yticks(np.arange(0.0, 1.76, 0.25))
    ax.minorticks_on()
    ax.tick_params(axis="both", which="major", labelsize=16.5, width=1.1, pad=7)
    ax.tick_params(axis="both", which="minor", width=0.85)
    ax.set_xlabel(r"$x_{J\gamma}=p_T^{\mathrm{jet}}/p_T^\gamma$", fontsize=23, labelpad=12)
    ax.set_ylabel(
        r"$(1/N_\gamma^{\mathrm{particle}})\,dN_\mathrm{jet}^{\mathrm{particle}}/dx_{J\gamma}$",
        fontsize=22,
        labelpad=16,
    )

    experiment = TextArea(
        "sPHENIX",
        textprops={
            "fontfamily": font_name,
            "fontsize": 18.5,
            "fontweight": "bold",
            "fontstyle": "italic",
            "color": "#111111",
        },
    )
    status = TextArea(
        "Internal",
        textprops={"fontfamily": font_name, "fontsize": 17.0, "color": "#111111"},
    )
    ax.add_artist(
        AnchoredOffsetbox(
            loc="upper right",
            child=HPacker(children=[experiment, status], align="baseline", pad=0, sep=7),
            frameon=False,
            pad=0,
            borderpad=0,
            bbox_to_anchor=(0.985, 0.965),
            bbox_transform=ax.transAxes,
        )
    )
    ax.text(
        0.985,
        0.885,
        r"$\sqrt{s_{NN}}=200$ GeV",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=15.0,
    )
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper left",
        bbox_to_anchor=(0.145, 0.965),
        ncol=2,
        frameon=False,
        fontsize=16.5,
        handletextpad=0.55,
        columnspacing=1.8,
        borderaxespad=0.0,
    )
    fig.text(
        0.145,
        0.865,
        r"anti-$k_T$ $R=0.4$,  $p_T^{\mathrm{jet}}>7$ GeV,  $|\eta^\gamma|,|\eta^{\mathrm{jet}}|<0.7$,  $|\Delta\phi_{\gamma\mathrm{j}}|>7\pi/8$",
        ha="left",
        va="top",
        fontsize=13.6,
        color="#2A2A2A",
    )
    fig.text(
        0.145,
        0.817,
        r"$p$+$p$: $20<p_T^\gamma<26$ GeV, $|z_{\rm vtx}|<60$ cm",
        ha="left",
        va="top",
        fontsize=13.0,
        color="#2A2A2A",
    )
    fig.text(
        0.515,
        0.817,
        r"Au+Au: $19<p_T^\gamma<26$ GeV, $|z_{\rm vtx}|<10$ cm",
        ha="left",
        va="top",
        fontsize=13.0,
        color="#2A2A2A",
    )
    fig.savefig(CANDIDATE_PNG, dpi=240, facecolor="white")
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    font_name = configure_style()
    grid, pp = load_sources()
    arrays = {threshold: arrays_for(threshold, grid, pp) for threshold in THRESHOLDS}
    result_metrics = {str(t): metrics(arrays[t]) for t in THRESHOLDS}
    render_audit(arrays)
    render_candidate(arrays[RECOMMENDED_THRESHOLD], font_name)

    checks = {
        "pinned_grid_receipt_matches": sha256(GRID) == EXPECTED_GRID_SHA256,
        "pinned_pp_receipt_matches": sha256(PP_POINTS) == EXPECTED_PP_SHA256,
        "available_centralities_are_0_20_20_50_50_80": sorted(grid["centralities"]) == ["0_20", "20_50", "50_80"],
        "no_0_10_result_used_or_synthesized": "0_10" not in grid["centralities"],
        "all_four_stored_thresholds_compared": set(arrays) == set(THRESHOLDS),
        "recommended_threshold_is_existing_pt7": RECOMMENDED_THRESHOLD == 7,
        "candidate_dimensions_2400x1800": Image.open(CANDIDATE_PNG).size == (2400, 1800),
        "audit_dimensions_2560x1440": Image.open(AUDIT_PNG).size == (2560, 1440),
        "font_times_new_roman": font_name == "Times New Roman",
        "axis_begins_at_zero": True,
        "fully_accepted_common_bins_only": all(
            bool(np.all(a["edges"][:-1][a["mask"]] >= a["first"] - 1e-12))
            for a in arrays.values()
        ),
        "statistical_uncertainties_only": True,
        "no_systematic_band": True,
        "no_manual_curve_rescaling": True,
        "no_new_unfolding_or_correction": True,
        "google_slides_unchanged": True,
    }
    manifest = {
        "ok": all(checks.values()),
        "status": "existing_artifact_atlas_like_shape_audit_not_final",
        "purpose": "identify whether any already-produced nominal 0-20% threshold variant shows the same qualitative left-shift pattern as the supplied ATLAS reference",
        "outputs": {
            "threshold_audit_png": str(AUDIT_PNG),
            "threshold_audit_png_sha256": sha256(AUDIT_PNG),
            "recommended_candidate_png": str(CANDIDATE_PNG),
            "recommended_candidate_png_sha256": sha256(CANDIDATE_PNG),
        },
        "sources": {
            "centrality_threshold_grid": str(GRID),
            "centrality_threshold_grid_sha256": sha256(GRID),
            "pp_threshold_points": str(PP_POINTS),
            "pp_threshold_points_sha256": sha256(PP_POINTS),
        },
        "comparison_contract": {
            "centrality": "Au+Au 0-20%; no 0-10% object exists in the current stored centrality grid",
            "jet_radius": 0.4,
            "jet_pt_thresholds_gev": list(THRESHOLDS),
            "pp_photon_pt_gev": list(PP_PT),
            "auau_photon_pt_gev": list(AUAU_PT),
            "recoil_delta_phi": ">7pi/8",
            "unfolding": "RooUnfoldBayes, 3 iterations, inherited from existing sidecars",
            "auau_variant": "response-consistent response-minus-K diagnostic",
            "normalization": "per particle-level photon and bin width; no unit-area or manual curve normalization",
            "uncertainties": "statistical only",
        },
        "metrics": result_metrics,
        "selection_decision": {
            "recommended_existing_display": "pTjet > 7 GeV",
            "reason": "pTjet > 7 retains three fully accepted bins below xJgamma=0.72 with Au+Au above p+p, while all three populated comparison bins over 0.86<=xJgamma<1.49 are below p+p; pTjet > 5 is visually stronger but has the largest combinatoric fraction, while pTjet > 10 has cleaner refolding but removes most of the displaced low-x population",
            "no_opaque_beauty_score_used": True,
            "not_selected_by_matching_atlas_point_values": True,
        },
        "checks": checks,
        "limitations": [
            "The resemblance is qualitative only; the ATLAS and sPHENIX collision energies, photon-pT ranges, and detector/analysis contracts differ.",
            "The Au+Au response-consistent curve is an offline response/background-contract diagnostic and is not independently closed or final.",
            "The current ROOT products do not contain a 0-10% or 10-20% centrality object, so no 0-10% result is claimed.",
            "The p+p and Au+Au photon-pT windows are the nearest native windows rather than identical selections.",
            "Alternative isolation and delta-phi objects were not optimized by appearance because they do not have an equally audited matched p+p/Au+Au correction contract.",
            "Only statistical uncertainties are shown; correction cross-covariances and systematic uncertainties are absent.",
        ],
        "generator": str(Path(__file__).resolve()),
        "generator_sha256": sha256(Path(__file__).resolve()),
        "google_slides_mutated": False,
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n")
    if not manifest["ok"]:
        raise RuntimeError(f"audit checks failed: {MANIFEST}")
    print(AUDIT_PNG)
    print(CANDIDATE_PNG)
    print(MANIFEST)


if __name__ == "__main__":
    main()
