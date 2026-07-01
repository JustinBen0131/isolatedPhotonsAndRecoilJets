#!/usr/bin/env python3
"""Mechanics slide for the Au+Au combinatoric recoil-jet subtraction.

The slide intentionally reuses the same reconstructed-level component builder
as the ATLAS-Fig.-1-style diagnostic so that the explanatory plot is tied to
the current correction chain rather than a separate approximation.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch


REPO = Path(__file__).resolve().parents[3]
THIS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(THIS_DIR))

import make_atlas_fig1_reco_background_breakdown as fig1  # noqa: E402


OUT_DIR = REPO / "dataOutput/the85_auau_xjgamma_unfolding_push/combinatoric_validation"
OUT_PNG = OUT_DIR / "slide11_combinatoric_background_mechanics.png"
OUT_MANIFEST = OUT_DIR / "slide11_combinatoric_background_mechanics_manifest.json"
OUT_SCRIPT = OUT_DIR / "slide11_combinatoric_background_mechanics_speaker_script.md"


INK = "#111827"
MUTED = "#64748b"
GRID = "#e5e7eb"
GRAY = "#7a7a7a"
RED = "#d62728"
BLUE = "#1f77b4"
GREEN = "#0f9d58"
PANEL = "#f8fafc"
PANEL_EDGE = "#cbd5e1"


def style() -> None:
    plt.rcParams.update(
        {
            "font.family": "Times New Roman",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "stix",
            "axes.linewidth": 1.2,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "xtick.major.size": 7,
            "ytick.major.size": 7,
            "xtick.minor.size": 4,
            "ytick.minor.size": 4,
        }
    )


def add_box(fig, xywh, title, body, edge=PANEL_EDGE, accent=BLUE, body_size=14.2):
    x, y, w, h = xywh
    box = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.010,rounding_size=0.012",
        transform=fig.transFigure,
        linewidth=1.4,
        edgecolor=edge,
        facecolor=PANEL,
        zorder=1,
    )
    fig.patches.append(box)
    fig.patches.append(
        FancyBboxPatch(
            (x, y),
            0.007,
            h,
            boxstyle="round,pad=0.0,rounding_size=0.010",
            transform=fig.transFigure,
            linewidth=0,
            facecolor=accent,
            zorder=2,
        )
    )
    fig.text(x + 0.020, y + h * 0.74, title, ha="left", va="center", fontsize=19.6, fontweight="bold", color=INK)
    fig.text(x + 0.020, y + h * 0.34, body, ha="left", va="center", fontsize=body_size, color=INK, linespacing=1.08)


def finite_ylim(*arrays: np.ndarray) -> float:
    values = np.concatenate([np.asarray(a, dtype=float) for a in arrays])
    values = values[np.isfinite(values)]
    if values.size == 0:
        return 1.0
    return max(1.0, float(np.nanmax(values)) * 1.22)


def main() -> None:
    style()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    panel = next(p for p in fig1.PANELS if p.key == "auau_0_20")
    comp = fig1.build_components(panel)
    xedges = comp["x_edges"]
    xcenters = comp["x_centers"]
    widths = np.diff(xedges)
    pre_comb = comp["bkg_sub"] + comp["comb"]
    pre_comb_err = np.sqrt(np.square(comp["bkg_sub_err"]) + np.square(comp["comb_err"]))
    comb = comp["comb"]
    bkg_sub = comp["bkg_sub"]
    bkg_sub_err = comp["bkg_sub_err"]

    mask = (xcenters >= 0.25) & (xcenters <= 1.90)
    total_pre = float(np.sum(pre_comb[mask]))
    total_comb = float(np.sum(comb[mask]))
    total_after = float(np.sum(bkg_sub[mask]))
    comb_frac = total_comb / total_pre if total_pre > 0 else math.nan

    fig = plt.figure(figsize=(16, 9), dpi=160)
    fig.patch.set_facecolor("white")

    fig.text(
        0.045,
        0.935,
        "Combinatoric recoil-jet background is subtracted before unfolding",
        ha="left",
        va="top",
        fontsize=31.5,
        fontweight="bold",
        color=INK,
    )

    ax = fig.add_axes([0.066, 0.170, 0.565, 0.675])
    ax.set_facecolor("white")
    ax.grid(True, which="major", color=GRID, linewidth=0.8, alpha=0.75)
    ax.minorticks_on()
    ax.set_axisbelow(True)

    ax.stairs(pre_comb, xedges, color=GRAY, linewidth=2.3, label="Before combinatoric subtraction")
    ax.stairs(comb, xedges, color=RED, linewidth=2.2, linestyle=":", label="Scaled combinatoric template")
    ax.errorbar(
        xcenters[mask],
        bkg_sub[mask],
        yerr=bkg_sub_err[mask],
        xerr=0.5 * widths[mask],
        fmt="o",
        color=INK,
        markersize=5.2,
        capsize=2.2,
        elinewidth=1.2,
        label="After combinatoric subtraction",
        zorder=5,
    )

    ax.set_xlim(0.2, 1.9)
    ax.set_ylim(0.0, finite_ylim(pre_comb[mask], comb[mask], bkg_sub[mask]))
    ax.set_xlabel(r"reconstructed $x_{J\gamma}$", fontsize=20)
    ax.set_ylabel("Entries", fontsize=20)
    ax.tick_params(labelsize=16)
    ax.legend(
        loc="upper right",
        frameon=True,
        facecolor="white",
        edgecolor="#d1d5db",
        fontsize=13.8,
        borderpad=0.55,
        handlelength=2.3,
    )
    ax.text(
        0.035,
        0.955,
        r"$\it{\bf{sPHENIX}}$ Internal" "\n"
        r"Au+Au 0-20%, $\sqrt{s_{NN}}=200$ GeV" "\n"
        r"default Au+Au BDT WP80",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=14.5,
        color=INK,
    )
    ax.text(
        0.985,
        0.675,
        "Selection shown\n"
        r"$15<E_T^\gamma<35$ GeV target" "\n"
        r"current bins: 16-35 GeV effective" "\n"
        r"$p_T^{jet}>5$ GeV, $R=0.4$" "\n"
        r"$|\Delta\phi_{\gamma j}|>7\pi/8$" "\n"
        r"$|\eta_\gamma|<0.7$, $|\eta_{jet}|<0.7$" "\n"
        r"$|z_{vtx}|<150$ cm; MBD N/S $\geq2$" "\n"
        r"sliding isolation $R=0.4$",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=11.9,
        color=INK,
        linespacing=1.06,
        bbox=dict(boxstyle="round,pad=0.34", facecolor="white", edgecolor="#d1d5db", alpha=0.96),
    )
    ax.text(
        0.985,
        0.375,
        rf"Template/input area = {comb_frac:.2f}",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=16.5,
        fontweight="bold",
        color=RED,
        bbox=dict(boxstyle="round,pad=0.28", facecolor="white", edgecolor="#fecaca", alpha=0.96),
    )

    x0 = 0.665
    w = 0.287
    add_box(
        fig,
        (x0, 0.702, w, 0.143),
        "1. Fill in embedded simulation",
        "$H^{comb}_{i}(x)$ = reco away-side jets\n"
        "tagged as unrelated to the generator\n"
        "photon/recoil pair in photon-$E_T$ row $i$.",
        accent=BLUE,
        body_size=13.8,
    )
    add_box(
        fig,
        (x0, 0.536, w, 0.143),
        "2. Scale to the data photon yield",
        "$B^{comb}_{i}(x)=H^{comb}_{i}(x)\\,N^{data}_{\\gamma,i}/N^{sim}_{\\gamma,i}$.\n"
        "The data photon yield is after\n"
        "photon-ID/purity correction.",
        accent=GREEN,
        body_size=13.8,
    )
    add_box(
        fig,
        (x0, 0.370, w, 0.143),
        "3. Subtract at reconstructed level",
        "$S^{reco}_{i}(x)=S^{ABCD}_{i}(x)-B^{comb}_{i}(x)$.\n"
        "This is after region-C dijet\n"
        "subtraction and before unfolding.",
        accent=RED,
        body_size=13.8,
    )
    add_box(
        fig,
        (x0, 0.204, w, 0.143),
        "4. Validation target",
        "Template should be centrality ordered,\n"
        "low-$x_{J\\gamma}$ dominated, zero in p+p,\n"
        "and not dominate the balanced peak.",
        accent="#7c3aed",
        body_size=13.8,
    )

    fig.patches.append(
        FancyBboxPatch(
            (0.070, 0.055),
            0.882,
            0.074,
            boxstyle="round,pad=0.012,rounding_size=0.012",
            transform=fig.transFigure,
            linewidth=1.2,
            edgecolor="#dbeafe",
            facecolor="#eff6ff",
        )
    )
    fig.text(0.092, 0.093, "Takeaway", fontsize=22, fontweight="bold", color=INK, va="center", ha="left")
    fig.text(
        0.190,
        0.093,
        "Additive Au+Au recoil-jet correction to reconstructed $x_{J\\gamma}$; it is separate from photon-ID purity and from detector unfolding.",
        fontsize=17.2,
        color=INK,
        va="center",
        ha="left",
    )

    fig.savefig(OUT_PNG, dpi=160)
    plt.close(fig)

    manifest = {
        "slide": str(OUT_PNG),
        "purpose": "First THE-87 mechanics slide explaining combinatoric recoil-jet subtraction placement.",
        "panel": panel.key,
        "data_file": str(panel.data_file),
        "sim_file": str(panel.sim_file),
        "histogram": fig1.hist_key("h2_unfoldRecoCombinatoric_pTgamma_xJ_incl", panel),
        "selection_displayed": {
            "centrality": "Au+Au 0-20%",
            "photon_et_target": "15 < E_T^gamma < 35 GeV",
            "photon_et_effective_current_bins": "16-35 GeV because current histograms straddle the 15 GeV edge",
            "jet_pt": "p_T^jet > 5 GeV",
            "jet_radius": "R=0.4",
            "back_to_back": "|Delta phi_gammaj| > 7pi/8",
            "photon_eta": "|eta_gamma| < 0.7",
            "jet_eta": "|eta_jet| < 0.7",
            "vertex": "|z_vtx| < 150 cm",
            "event": "MBD N/S >= 2",
            "isolation": "sliding isolation R=0.4",
            "photon_id": "default AuAu BDT WP80",
        },
        "correction_chain": [
            "Region A raw xJ is corrected with leakage-aware region-C photon-ID/dijet sideband subtraction.",
            "The combinatoric template is taken from embedded signal MC and scaled per photon-E_T row by data/sim photon yield.",
            "The scaled combinatoric template is subtracted from the reconstructed xJ input before response-matrix unfolding.",
        ],
        "display_window": {"x_min": 0.25, "x_max": 1.90},
        "integrals_display_window": {
            "pre_combinatoric_S_ABCD": total_pre,
            "combinatoric_template": total_comb,
            "after_combinatoric_S_reco": total_after,
            "template_over_S_ABCD": comb_frac,
        },
        "source_component_integrals_full": comp["integrals"],
        "pt_rows": comp["pt_rows"],
        "purity_fit": comp["purity_fit"],
        "notes": [
            "Diagnostic is reconstructed-level and before response-matrix unfolding.",
            "pp combinatoric background is treated as zero for this ATLAS-style chain.",
            "Current strict-full-bin inputs are effectively 16-35 GeV because existing histograms straddle the 15 GeV edge.",
        ],
    }
    OUT_MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n")

    OUT_SCRIPT.write_text(
        "\n".join(
            [
                "# Speaker note",
                "",
                "This slide separates the Au+Au combinatoric recoil-jet correction from the photon-ID sideband correction.",
                "The red template is filled in embedded simulation from reco away-side jets not matched to the generator photon/recoil hard-scatter pair.",
                "Offline, it is scaled in each photon-E_T row by the corrected data photon yield over the corresponding signal-MC photon yield.",
                "That scaled template is subtracted after the ABCD photon-ID sideband correction and before the RooUnfold response-matrix step.",
                "For the current 0-20% input, the template is large enough that it deserves its own validation block before we interpret the final xJ shape.",
                "",
            ]
        )
    )

    print(OUT_PNG)
    print(OUT_MANIFEST)
    print(OUT_SCRIPT)


if __name__ == "__main__":
    main()
