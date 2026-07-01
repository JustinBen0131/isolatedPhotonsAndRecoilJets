#!/usr/bin/env python3
"""Full-slide comparison of reconstructed xJ inputs for jet pT thresholds.

This slide uses existing THE-85/THE-69 ROOT products.  It compares the current
analysis reconstructed inputs for the nominal jet threshold and the stored
`jetPt10` internal-scan threshold without rerunning Fun4All.
"""

from __future__ import annotations

import importlib
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.font_manager import FontProperties, findfont


REPO = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import make_atlas_fig1_reco_background_breakdown as fig1  # noqa: E402
import make_pp_atlas_vs_sphenix_fig1_comparison as ppcomp  # noqa: E402


OUT_DIR = REPO / "dataOutput/the85_auau_xjgamma_unfolding_push/unfolded_xjgamma_firstpass"
OUT_PNG = OUT_DIR / "slide10_sphenix_reco_inputs_ptjet5_vs10.png"
OUT_MANIFEST = OUT_DIR / "slide10_sphenix_reco_inputs_ptjet5_vs10_manifest.json"
OUT_SCRIPT = OUT_DIR / "slide10_sphenix_reco_inputs_ptjet5_vs10_speaker_script.md"

TITLE_COLOR = "#111827"
BODY_COLOR = "#334155"
GREY = "#8f8f8f"
RED = "#d62728"
BLUE = "#1f4cff"


def set_times() -> None:
    """Use Times New Roman when available; fall back to a Times-compatible serif."""
    try:
        findfont(FontProperties(family="Times New Roman"), fallback_to_default=False)
        family = "Times New Roman"
    except Exception:
        family = "Times"
    plt.rcParams.update(
        {
            "font.family": family,
            "font.serif": ["Times New Roman", "Times", "Nimbus Roman", "DejaVu Serif"],
            "mathtext.fontset": "stix",
            "axes.linewidth": 1.15,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "savefig.facecolor": "white",
        }
    )


def build_for_key(jet_key: str) -> tuple[dict, dict, fig1.Panel]:
    fig1.JET_PT_KEY = jet_key
    importlib.reload(ppcomp)
    pp_result = ppcomp.project_pp_with_ppg12_purity_norm()
    auau_panel = next(panel for panel in fig1.PANELS if panel.key == "auau_0_20")
    auau_result = fig1.build_components(auau_panel)
    return pp_result, auau_result, auau_panel


def jet_label(jet_key: str) -> str:
    return r"$p_T^{jet}>5$ GeV" if not jet_key else r"$p_T^{jet}>10$ GeV"


def draw_panel(ax, result: dict, *, sample_label: str, jet_key: str, show_ylabel: bool, show_xlabel: bool) -> None:
    xedges = result["x_edges"]
    centers = result["x_centers"]
    widths = np.diff(xedges)
    mask = (centers >= 0.2) & (centers <= 1.85)
    raw = result["raw"]
    side = result["sideband"]
    comb = result["comb"]
    bkg = result["bkg_sub"]
    bkg_err = result["bkg_sub_err"]

    ax.stairs(raw, xedges, color=GREY, lw=2.4, label="Raw region A")
    if float(np.sum(comb)) > 0:
        ax.stairs(comb, xedges, color=RED, lw=2.15, linestyle=(0, (1, 1)), label="Comb. bkg.")
    else:
        ax.plot([xedges[0], xedges[-1]], [0.0, 0.0], color=RED, lw=2.0, linestyle=(0, (1, 1)), label="Comb. bkg. (0)")
    if float(np.sum(side)) > 0:
        ax.stairs(side, xedges, color=BLUE, lw=2.3, linestyle=(0, (2, 2)), label="ABCD ID sideband")
    ax.errorbar(
        centers[mask],
        bkg[mask],
        xerr=0.5 * widths[mask],
        yerr=bkg_err[mask],
        fmt="o",
        color="black",
        ms=5.2,
        elinewidth=1.05,
        capsize=0,
        label="Bkg.-sub. input",
        zorder=5,
    )

    ymax = max(1.0, float(np.nanmax(raw[mask]) * 1.25))
    ax.set_xlim(0.2, 1.85)
    ax.set_ylim(0.0, ymax)
    ax.minorticks_on()
    ax.tick_params(axis="both", which="major", labelsize=16, top=True, right=True, length=7)
    ax.tick_params(axis="both", which="minor", top=True, right=True, length=3.5)
    if show_xlabel:
        ax.set_xlabel(r"Reconstructed $x_{J\gamma}$", fontsize=20, labelpad=6)
    else:
        ax.set_xticklabels([])
    if show_ylabel:
        ax.set_ylabel("Entries", fontsize=20, labelpad=7)
    else:
        ax.set_yticklabels([])

    ax.text(0.045, 0.925, r"$\bf{\it{sPHENIX}}$ Internal", transform=ax.transAxes, ha="left", va="top", fontsize=20)
    ax.text(0.045, 0.815, sample_label, transform=ax.transAxes, ha="left", va="top", fontsize=20)
    ax.text(0.045, 0.715, r"$15<E_T^\gamma<35$ GeV, $|\Delta\phi|>7\pi/8$", transform=ax.transAxes, ha="left", va="top", fontsize=16.5)
    ax.text(0.045, 0.625, jet_label(jet_key), transform=ax.transAxes, ha="left", va="top", fontsize=16.5, color=BODY_COLOR)


def draw() -> None:
    set_times()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    pp5, auau5, auau_panel = build_for_key("")
    pp10, auau10, _ = build_for_key("jetPt10")

    fig = plt.figure(figsize=(16, 9), dpi=160)
    fig.patch.set_facecolor("white")

    fig.text(
        0.055,
        0.955,
        r"Reconstructed $x_{J\gamma}$ inputs sharpen with a higher recoil-jet threshold",
        ha="left",
        va="top",
        fontsize=28,
        fontweight="bold",
        color=TITLE_COLOR,
    )

    left, right = 0.065, 0.955
    bottom, top = 0.105, 0.835
    wspace, hspace = 0.052, 0.092
    panel_w = (right - left - wspace) / 2.0
    panel_h = (top - bottom - hspace) / 2.0
    axes = {
        "pp5": fig.add_axes([left, bottom + panel_h + hspace, panel_w, panel_h]),
        "auau5": fig.add_axes([left + panel_w + wspace, bottom + panel_h + hspace, panel_w, panel_h]),
        "pp10": fig.add_axes([left, bottom, panel_w, panel_h]),
        "auau10": fig.add_axes([left + panel_w + wspace, bottom, panel_w, panel_h]),
    }

    draw_panel(axes["pp5"], pp5, sample_label="p+p", jet_key="", show_ylabel=True, show_xlabel=False)
    draw_panel(axes["auau5"], auau5, sample_label="Au+Au 0-20%", jet_key="", show_ylabel=False, show_xlabel=False)
    draw_panel(axes["pp10"], pp10, sample_label="p+p", jet_key="jetPt10", show_ylabel=True, show_xlabel=True)
    draw_panel(axes["auau10"], auau10, sample_label="Au+Au 0-20%", jet_key="jetPt10", show_ylabel=False, show_xlabel=True)

    handles, labels = axes["auau10"].get_legend_handles_labels()
    keep = []
    seen = set()
    for h, lab in zip(handles, labels):
        if lab not in seen:
            seen.add(lab)
            keep.append((h, lab))
    fig.legend(
        [h for h, _lab in keep],
        [lab for _h, lab in keep],
        loc="upper right",
        bbox_to_anchor=(0.955, 0.955),
        ncol=4,
        frameon=False,
        fontsize=15,
        handlelength=2.6,
        columnspacing=1.6,
    )

    manifest = {
        "slide": str(OUT_PNG),
        "source": "Current THE-85 reconstructed-input diagnostic using stored internal jet-pT scan histograms.",
        "selection_common": "15<E_T^gamma<35 GeV, |Delta phi|>7pi/8",
        "status": "Reconstructed input before unfolding; not a final particle-level result.",
        "panels": {
            "pp_jetPt5": {
                "effective_base_key": "r04",
                "integrals": pp5["integrals"],
                "pt_rows": pp5["pt_rows"],
            },
            "pp_jetPt10": {
                "effective_base_key": "r04_jetPt10",
                "integrals": pp10["integrals"],
                "pt_rows": pp10["pt_rows"],
            },
            "auau_0_20_jetPt5": {
                "effective_base_key": "r04_isoR40_isSliding",
                "integrals": auau5["integrals"],
                "pt_rows": auau5["pt_rows"],
                "data_file": str(auau_panel.data_file),
                "sim_file": str(auau_panel.sim_file),
            },
            "auau_0_20_jetPt10": {
                "effective_base_key": "r04_jetPt10_isoR40_isSliding",
                "integrals": auau10["integrals"],
                "pt_rows": auau10["pt_rows"],
                "data_file": str(auau_panel.data_file),
                "sim_file": str(auau_panel.sim_file),
            },
        },
        "caveats": [
            "Blue is the ABCD photon-ID sideband subtraction used by this pipeline, not a standalone ATLAS dijet template.",
            "The stricter jet threshold is available from existing internal-scan histograms; no new production was run.",
        ],
    }
    fig.savefig(OUT_PNG, dpi=160)
    plt.close(fig)

    OUT_MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n")
    OUT_SCRIPT.write_text(
        "This slide compares the same reconstructed xJgamma input at two recoil-jet thresholds. "
        "The top row uses the nominal five-GeV minimum jet pT. The bottom row uses the stored ten-GeV internal-scan histograms. "
        "The visual point is that the ten-GeV threshold removes much of the low-xJ turn-on region and makes the pp input cleaner, while the central AuAu background-subtracted input is still limited by the subtraction composition rather than the threshold alone.\n\n"
        "I would present this as a diagnostic choice, not a final physics result. "
        "If we want the ATLAS-style subtraction picture to be final-quality, the next fix is to separate the photon-ID sideband and heavy-ion combinatoric components more rigorously, not just tune the jet threshold.\n"
    )

    print(json.dumps({"slide": str(OUT_PNG), "manifest": str(OUT_MANIFEST), "speaker_script": str(OUT_SCRIPT)}, indent=2))


if __name__ == "__main__":
    draw()
