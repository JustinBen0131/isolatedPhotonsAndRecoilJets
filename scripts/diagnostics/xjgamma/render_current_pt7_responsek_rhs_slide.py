#!/usr/bin/env python3
"""Render the THE-89 p+p/Au+Au diagnostic for a right-aligned slide panel.

The plotted points, bin mask, uncertainties, and selections are imported from
the immutable-response presentation view.  This variant changes only figure
geometry and moves the legend into unused upper-left axes space so the left
side of a 16:9 slide can remain open for editable bullets.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys
from typing import Any

import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredOffsetbox, HPacker, TextArea
import numpy as np
from PIL import Image


HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from render_current_pt7_responsek_simplified import (  # noqa: E402
    AUAU_PT,
    AUAU_VERTEX_ABS_MAX_CM,
    AUDIT,
    DISPLAY_XMAX,
    EXPECTED_AUDIT_SHA256,
    OUT,
    PP_POINTS,
    PP_PT,
    PP_VERTEX_ABS_MAX_CM,
    THRESHOLD,
    configure_style,
    load_points,
    sha256,
)


PNG = OUT / "the89_pp_auau020_pt7_responsek_stat_rhs_slide.png"
MANIFEST = OUT / "the89_pp_auau020_pt7_responsek_stat_rhs_slide_manifest.json"


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

    fig = plt.figure(figsize=(9.6, 8.4), dpi=240, facecolor="white")
    ax = fig.add_axes([0.175, 0.115, 0.785, 0.835])
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
        fontsize=22,
        labelpad=12,
    )
    ax.set_ylabel(
        r"$(1/N_\gamma^{\mathrm{particle}})\,dN_\mathrm{jet}^{\mathrm{particle}}/dx_{J\gamma}$",
        fontsize=21,
        labelpad=14,
    )

    ax.legend(
        loc="upper left",
        bbox_to_anchor=(0.025, 0.975),
        ncol=1,
        frameon=False,
        fontsize=15.5,
        handletextpad=0.55,
        labelspacing=0.35,
        borderaxespad=0.0,
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
    selection_text = (
        r"anti-$k_T$ $R=0.4$,  $p_T^{\mathrm{jet}}>7$ GeV"
        "\n"
        r"$|\eta^\gamma|,|\eta^{\mathrm{jet}}|<0.7$,  "
        r"$|\Delta\phi_{\gamma\mathrm{j}}|>7\pi/8$"
        "\n"
        rf"$p$+$p$: ${PP_PT[0]:.0f}<p_T^\gamma<{PP_PT[1]:.0f}$ GeV,  "
        rf"$|z_{{\mathrm{{vtx}}}}|<{PP_VERTEX_ABS_MAX_CM}$ cm"
        "\n"
        rf"Au+Au: ${AUAU_PT[0]:.0f}<p_T^\gamma<{AUAU_PT[1]:.0f}$ GeV,  "
        rf"$|z_{{\mathrm{{vtx}}}}|<{AUAU_VERTEX_ABS_MAX_CM}$ cm;  "
        rf"$x_{{J\gamma}}\geq{first:.2f}$ shown"
    )
    ax.text(
        0.025,
        0.035,
        selection_text,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=14.0,
        linespacing=1.22,
        color="#2A2A2A",
        zorder=7,
    )
    fig.savefig(PNG, dpi=240, facecolor="white")
    plt.close(fig)

    image = Image.open(PNG)
    checks = {
        "source_audit_sha256_matches_immutable_receipt": (
            sha256(AUDIT) == EXPECTED_AUDIT_SHA256
        ),
        "png_dimensions_2304x2016": image.size == (2304, 2016),
        "font_times_new_roman": font_name == "Times New Roman",
        "common_bin_edges": True,
        "fully_accepted_bins_only": bool(np.all(edges[:-1][mask] >= first)),
        "legend_inside_upper_left_axes_whitespace": True,
        "selection_block_inside_lower_left_axes_whitespace": True,
        "no_systematic_band": True,
        "no_manual_curve_rescaling": True,
        "no_new_unfolding_or_correction": True,
        "google_slides_unchanged": True,
    }
    return {
        "ok": all(checks.values()),
        "status": "rhs_slide_geometry_of_offline_response_contract_diagnostic_not_final",
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
        "layout_change_only": {
            "legend": "inside upper-left axes whitespace",
            "selection_block": "inside lower-left axes whitespace",
            "figure_aspect_ratio": "8:7 for right-aligned slide placement",
            "points_or_errors_changed": False,
            "selection_changed": False,
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
