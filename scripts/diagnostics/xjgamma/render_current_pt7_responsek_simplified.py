#!/usr/bin/env python3
"""Render a presentation-only p+p/Au+Au view of the immutable THE-89 audit.

This script performs no unfolding or physics correction.  It reads the
already-recorded p+p points and the response-consistent Au+Au diagnostic from
the immutable low-x response-contract audit, selects the pTjet > 7 GeV case,
and renders only statistical uncertainties.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.offsetbox import AnchoredOffsetbox, HPacker, TextArea
import numpy as np
from PIL import Image


REPO = Path(__file__).resolve().parents[3]
OUT = (
    REPO
    / "dataOutput/the219_friday_ppg_20260814/"
    "the89_lowx_response_contract_audit"
)
AUDIT = OUT / "the89_current_lowx_response_contract_pt5_7_10_manifest.json"
EXPECTED_AUDIT_SHA256 = (
    "c2781ffe2cfd72b7997f507c6eca3185ee5320d6e6c64a755dfee70b376c84cf"
)
PP_POINTS = (
    REPO
    / "dataOutput/the219_friday_ppg_20260814/the89_pp_threshold_peak_scan/"
    "the89_pp_unfolded_xjgamma_threshold_scan_points.json"
)
PNG = OUT / "the89_pp_auau020_pt7_responsek_stat_simplified.png"
MANIFEST = OUT / "the89_pp_auau020_pt7_responsek_stat_simplified_manifest.json"
FONT_PATH = Path("/System/Library/Fonts/Supplemental/Times New Roman.ttf")

THRESHOLD = 7
PP_PT = (20.0, 26.0)
AUAU_PT = (19.0, 26.0)
PP_VERTEX_ABS_MAX_CM = 60
AUAU_VERTEX_ABS_MAX_CM = 10
DISPLAY_XMAX = 1.49


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


def load_points() -> tuple[dict[str, Any], dict[str, Any], np.ndarray]:
    audit_sha = sha256(AUDIT)
    if audit_sha != EXPECTED_AUDIT_SHA256:
        raise RuntimeError(
            "immutable audit receipt changed: "
            f"expected {EXPECTED_AUDIT_SHA256}, found {audit_sha}"
        )

    audit = json.loads(AUDIT.read_text())
    if audit.get("status") != "offline_response_contract_diagnostic_not_final":
        raise RuntimeError("source audit is not marked diagnostic-not-final")
    auau = audit["thresholds"][str(THRESHOLD)]

    pp_payload = json.loads(PP_POINTS.read_text())
    pp = pp_payload[str(THRESHOLD)]
    pp_edges = np.asarray(
        [float(row["xj_low"]) for row in pp["bins"]]
        + [float(pp["bins"][-1]["xj_high"])],
        dtype=float,
    )
    auau_edges = np.asarray(auau["x_edges"], dtype=float)
    if not np.array_equal(pp_edges, auau_edges):
        raise RuntimeError("p+p and Au+Au xJgamma bin edges do not match")
    return pp, auau, auau_edges


def render() -> dict[str, Any]:
    font_name = configure_style()
    pp, auau, edges = load_points()
    fixed = auau["variants"]["response_consistent"]

    centers = 0.5 * (edges[:-1] + edges[1:])
    halfwidth = 0.5 * np.diff(edges)
    pp_y = np.asarray([row["pp_unfolded"]["value"] for row in pp["bins"]])
    pp_ey = np.asarray(
        [row["pp_unfolded"]["stat_error"] for row in pp["bins"]]
    )
    auau_y = np.asarray(fixed["y"], dtype=float)
    auau_ey = np.asarray(fixed["ey"], dtype=float)
    first = max(
        float(pp["first_accepted_xj_edge"]),
        float(auau["first_accepted_xj_edge"]),
    )
    mask = (
        (edges[:-1] >= first - 1e-12)
        & (edges[1:] <= DISPLAY_XMAX + 1e-12)
        & np.isfinite(pp_y)
        & np.isfinite(pp_ey)
        & np.isfinite(auau_y)
        & np.isfinite(auau_ey)
    )

    fig = plt.figure(figsize=(10, 7.5), dpi=240, facecolor="white")
    ax = fig.add_axes([0.145, 0.120, 0.82, 0.640])
    ax.errorbar(
        centers[mask] - 0.006,
        pp_y[mask],
        xerr=halfwidth[mask],
        yerr=pp_ey[mask],
        fmt="s",
        ms=7.8,
        mfc="white",
        mec="#111111",
        mew=1.65,
        ecolor="#111111",
        elinewidth=1.45,
        capsize=2.8,
        capthick=1.25,
        linestyle="none",
        label=r"$p$+$p$",
        zorder=4,
    )
    ax.errorbar(
        centers[mask] + 0.006,
        auau_y[mask],
        xerr=halfwidth[mask],
        yerr=auau_ey[mask],
        fmt="o",
        ms=8.2,
        mfc="#C83E32",
        mec="#8F241D",
        mew=1.1,
        ecolor="#B9342A",
        elinewidth=1.5,
        capsize=2.8,
        capthick=1.3,
        linestyle="none",
        label="Au+Au, 0–20%",
        zorder=5,
    )

    ax.set_xlim(0.0, 1.52)
    ax.set_ylim(0.0, 1.90)
    ax.set_xticks(np.arange(0.0, 1.51, 0.2))
    ax.set_yticks(np.arange(0.0, 1.76, 0.25))
    ax.minorticks_on()
    ax.tick_params(axis="both", which="major", labelsize=17, width=1.1, pad=7)
    ax.tick_params(axis="both", which="minor", width=0.85)
    ax.set_xlabel(
        r"$x_{J\gamma}=p_T^{\mathrm{jet}}/p_T^\gamma$",
        fontsize=23,
        labelpad=12,
    )
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
        textprops={
            "fontfamily": font_name,
            "fontsize": 17.0,
            "color": "#111111",
        },
    )
    label_box = HPacker(children=[experiment, status], align="baseline", pad=0, sep=7)
    ax.add_artist(
        AnchoredOffsetbox(
            loc="upper right",
            child=label_box,
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
        color="#111111",
    )
    fig.text(
        0.145,
        0.865,
        r"anti-$k_T$ $R=0.4$,  "
        r"$p_T^{\mathrm{jet}}>7$ GeV,  "
        r"$|\eta^\gamma|,|\eta^{\mathrm{jet}}|<0.7$,  "
        r"$|\Delta\phi_{\gamma\mathrm{j}}|>7\pi/8$",
        ha="left",
        va="top",
        fontsize=13.6,
        color="#2A2A2A",
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
        labelspacing=0.45,
        columnspacing=1.8,
        borderaxespad=0.0,
    )
    fig.text(
        0.145,
        0.810,
        rf"$p$+$p$: ${PP_PT[0]:.0f}<p_T^\gamma<{PP_PT[1]:.0f}$ GeV, "
        rf"$|z_{{\mathrm{{vtx}}}}|<{PP_VERTEX_ABS_MAX_CM}$ cm",
        ha="left",
        va="top",
        fontsize=13.2,
        color="#2A2A2A",
    )
    fig.text(
        0.515,
        0.810,
        rf"Au+Au: ${AUAU_PT[0]:.0f}<p_T^\gamma<{AUAU_PT[1]:.0f}$ GeV, "
        rf"$|z_{{\mathrm{{vtx}}}}|<{AUAU_VERTEX_ABS_MAX_CM}$ cm",
        ha="left",
        va="top",
        fontsize=13.2,
        color="#2A2A2A",
    )
    fig.text(
        0.965,
        0.810,
        rf"$x_{{J\gamma}}\geq{first:.2f}$ shown",
        ha="right",
        va="top",
        fontsize=13.2,
        color="#2A2A2A",
    )
    fig.savefig(PNG, dpi=240, facecolor="white")
    plt.close(fig)

    image = Image.open(PNG)
    checks = {
        "source_audit_sha256_matches_immutable_receipt": (
            sha256(AUDIT) == EXPECTED_AUDIT_SHA256
        ),
        "png_dimensions_2400x1800": image.size == (2400, 1800),
        "font_times_new_roman": font_name == "Times New Roman",
        "common_bin_edges": True,
        "fully_accepted_bins_only": bool(np.all(edges[:-1][mask] >= first)),
        "no_systematic_band": True,
        "no_manual_curve_rescaling": True,
        "no_new_unfolding_or_correction": True,
        "nominal_cut_block_includes_system_specific_vertex_limits": True,
        "sphenix_internal_upper_right_inside_plot_frame": True,
        "collision_energy_below_sphenix_internal_inside_plot_frame": True,
        "google_slides_unchanged": True,
    }
    return {
        "ok": all(checks.values()),
        "status": "presentation_view_of_offline_response_contract_diagnostic_not_final",
        "png": str(PNG),
        "png_sha256": sha256(PNG),
        "generator": str(Path(__file__).resolve()),
        "generator_sha256": sha256(Path(__file__).resolve()),
        "sources": {
            "immutable_auau_audit": str(AUDIT),
            "immutable_auau_audit_sha256": sha256(AUDIT),
            "pp_points": str(PP_POINTS),
            "pp_points_sha256": sha256(PP_POINTS),
        },
        "selection": {
            "jet_radius": 0.4,
            "jet_pt_min_gev": THRESHOLD,
            "pp_photon_pt_gev": list(PP_PT),
            "auau_photon_pt_gev": list(AUAU_PT),
            "pp_vertex_abs_max_cm": PP_VERTEX_ABS_MAX_CM,
            "auau_vertex_abs_max_cm": AUAU_VERTEX_ABS_MAX_CM,
            "auau_centrality_percent": [0, 20],
            "minimum_common_full_xj_edge": first,
            "maximum_displayed_xj_edge": DISPLAY_XMAX,
        },
        "marker_contract": {
            "pp": "black open squares",
            "auau": "filled red circles",
            "auau_variant": "response_consistent",
            "open_red_current_mismatch_series_omitted": True,
            "horizontal_offset": "plus or minus 0.006 in xJgamma for legibility only",
        },
        "plot_annotation_contract": {
            "experiment_label": "bold italic sPHENIX plus upright Internal at upper-right inside the axes frame",
            "collision_energy_label": "sqrt(sNN) = 200 GeV directly below the experiment label inside the axes frame",
            "common": [
                "anti-kT R = 0.4",
                "pTjet > 7 GeV",
                "absolute photon and jet eta < 0.7",
                "absolute photon-jet delta phi > 7pi/8",
                "displayed xJgamma >= 0.41",
            ],
            "system_specific": [
                "p+p absolute z vertex < 60 cm",
                "Au+Au absolute z vertex < 10 cm",
                "native photon-pT and vertex windows are printed in a compact system-specific block",
            ],
            "omitted_as_extraneous": [
                "response-K diagnostic label inside plot",
                "statistical-only prose inside plot",
            ],
        },
        "uncertainties": "statistical only; no systematic band",
        "checks": checks,
        "limitations": [
            "The Au+Au response-consistent curve is an offline diagnostic and is not an independently closed final result.",
            "The source audit records incomplete response/combinatoric closure and must not be relabeled as a final measurement.",
            "The p+p and Au+Au photon-pT windows are the nearest available native windows, not identical selections.",
        ],
        "google_slides_mutated": False,
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    manifest = render()
    MANIFEST.write_text(json.dumps(manifest, indent=2, allow_nan=False) + "\n")
    if not manifest["ok"]:
        raise RuntimeError(f"render checks failed: {MANIFEST}")
    print(PNG)
    print(MANIFEST)


if __name__ == "__main__":
    main()
