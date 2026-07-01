#!/usr/bin/env python3
"""Shape-normalized pp vs AuAu 0-20 reconstructed-xJ overlay for THE-85.

This uses the same current-analysis objects as the slide-14 ATLAS/sPHENIX
Fig.-1 analogue, but overlays pp and AuAu 0-20 directly.  It is a
reconstructed-input diagnostic before unfolding, not the final unfolded
per-photon xJgamma result.
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path
from typing import Dict, Tuple

import matplotlib.pyplot as plt
import numpy as np

# Match the strict-full-bin slide-14 object unless the caller deliberately
# overrides the environment before importing the existing helpers.
os.environ.setdefault("THE85_REQUIRE_FULL_PT_BINS", "1")

REPO = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import make_atlas_fig1_reco_background_breakdown as fig1  # noqa: E402
import make_pp_atlas_vs_sphenix_fig1_comparison as ppcomp  # noqa: E402


OUT_DIR = REPO / "dataOutput/the85_auau_xjgamma_unfolding_push/unfolded_xjgamma_firstpass"
OUT_PNG = OUT_DIR / "slide15_sphenix_pp_auau020_reco_shape_overlay.png"
OUT_MANIFEST = OUT_DIR / "slide15_sphenix_pp_auau020_reco_shape_overlay_manifest.json"
OUT_SCRIPT = OUT_DIR / "slide15_sphenix_pp_auau020_reco_shape_overlay_speaker_script.md"

XMIN = 0.2
XMAX = 1.85

PP_COLOR = "#d62728"
AUAU_COLOR = "#1f4cff"
TITLE_COLOR = "#111827"
BODY_COLOR = "#334155"


def normalize_density(
    values: np.ndarray,
    errors: np.ndarray,
    xedges: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, float]:
    """Return an area-normalized density over the displayed x range."""
    centers = 0.5 * (xedges[:-1] + xedges[1:])
    widths = np.diff(xedges)
    mask = (centers >= XMIN) & (centers <= XMAX)
    area = float(np.sum(values[mask] * widths[mask]))
    if not np.isfinite(area) or abs(area) < 1.0e-12:
        return np.zeros_like(values), np.zeros_like(errors), area
    return values / area, errors / abs(area), area


def build_inputs() -> Dict[str, Dict]:
    pp = ppcomp.project_pp_with_ppg12_purity_norm()
    auau_panel = next(panel for panel in fig1.PANELS if panel.key == "auau_0_20")
    auau = fig1.build_components(auau_panel)
    return {"pp": pp, "auau_0_20": auau}


def draw_curve(
    ax,
    result: Dict,
    key: str,
    mode: str,
    color: str,
    label: str,
    marker: str,
) -> Dict:
    xedges = result["x_edges"]
    centers = result["x_centers"]
    widths = np.diff(xedges)
    mask = (centers >= XMIN) & (centers <= XMAX)
    if mode == "raw":
        values = result["raw"]
        errors = result.get("raw_err", np.sqrt(np.clip(values, 0, None)))
    elif mode == "corrected":
        values = result["bkg_sub"]
        errors = result["bkg_sub_err"]
    else:
        raise ValueError(mode)

    density, density_err, area = normalize_density(values, errors, xedges)
    ax.stairs(density, xedges, color=color, lw=2.4, alpha=0.82)
    ax.errorbar(
        centers[mask],
        density[mask],
        xerr=0.5 * widths[mask],
        yerr=density_err[mask],
        fmt=marker,
        color=color,
        mfc="white" if key == "pp" else color,
        mec=color,
        mew=1.5,
        ms=7.0,
        elinewidth=1.2,
        capsize=0,
        label=label,
        zorder=5,
    )
    return {
        "area": area,
        "density_integral_check": float(np.sum(density[mask] * widths[mask])),
        "raw_integral": float(np.sum(values[mask])),
    }


def style_axis(ax, title: str, ylabel: bool = False) -> None:
    ax.set_xlim(XMIN, XMAX)
    ax.set_xlabel(r"Reconstructed $x_{J\gamma}$", fontsize=21, labelpad=6)
    if ylabel:
        ax.set_ylabel(r"Area-normalized $(1/N)\,dN/dx_{J\gamma}$", fontsize=20, labelpad=9)
    ax.set_title(title, fontsize=22, fontweight="bold", color=TITLE_COLOR, pad=12)
    ax.minorticks_on()
    ax.tick_params(labelsize=16.5, top=True, right=True, direction="in", length=7)
    ax.tick_params(which="minor", top=True, right=True, direction="in", length=3.5)
    ax.grid(True, which="major", color="#cbd5e1", alpha=0.38, lw=0.8)
    for spine in ax.spines.values():
        spine.set_linewidth(1.25)


def main() -> None:
    plt.rcParams.update(
        {
            "font.family": "Times New Roman",
            "font.serif": ["Times New Roman", "Times", "Nimbus Roman", "DejaVu Serif"],
            "mathtext.fontset": "stix",
            "axes.linewidth": 1.25,
        }
    )
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    inputs = build_inputs()

    fig = plt.figure(figsize=(16, 9), dpi=160)
    fig.patch.set_facecolor("white")

    fig.text(
        0.055,
        0.944,
        r"Shape-normalized pp vs Au+Au 0-20% reco $x_{J\gamma}$ inputs",
        ha="left",
        va="top",
        fontsize=30,
        fontweight="bold",
        color=TITLE_COLOR,
    )
    fig.text(
        0.055,
        0.888,
        (
            r"Same slide-14 objects; each curve has unit area over $0.2<x_{J\gamma}<1.85$. "
            "This is a shape check before unfolding, not a final yield comparison."
        ),
        ha="left",
        va="top",
        fontsize=18.0,
        color=BODY_COLOR,
    )
    fig.text(
        0.055,
        0.846,
        (
            r"Selection: effective 16-35 GeV photon binning, $|\Delta\phi|>7\pi/8$, $p_T^{jet}>5$ GeV. "
            "Corrected = photon-ID sideband plus Au+Au combinatoric subtraction."
        ),
        ha="left",
        va="top",
        fontsize=15.8,
        color=BODY_COLOR,
    )

    ax_raw = fig.add_axes([0.070, 0.105, 0.415, 0.660])
    ax_corr = fig.add_axes([0.545, 0.105, 0.415, 0.660])

    stats: Dict[str, Dict] = {"raw": {}, "corrected": {}}
    stats["raw"]["pp"] = draw_curve(ax_raw, inputs["pp"], "pp", "raw", PP_COLOR, "p+p", "s")
    stats["raw"]["auau_0_20"] = draw_curve(
        ax_raw, inputs["auau_0_20"], "auau_0_20", "raw", AUAU_COLOR, "Au+Au 0-20%", "o"
    )
    stats["corrected"]["pp"] = draw_curve(
        ax_corr, inputs["pp"], "pp", "corrected", PP_COLOR, "p+p", "s"
    )
    stats["corrected"]["auau_0_20"] = draw_curve(
        ax_corr, inputs["auau_0_20"], "auau_0_20", "corrected", AUAU_COLOR, "Au+Au 0-20%", "o"
    )

    style_axis(ax_raw, "Uncorrected: raw region A", ylabel=True)
    style_axis(ax_corr, "Corrected: input to unfolding", ylabel=False)

    for ax in (ax_raw, ax_corr):
        ax.legend(loc="upper right", frameon=False, fontsize=18, handlelength=2.5, labelspacing=0.55)
        ax.text(
            0.040,
            0.955,
            r"$\bf{\it{sPHENIX}}$ Internal",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=18.5,
        )
        ax.text(
            0.040,
            0.890,
            r"$\sqrt{s_{NN}}=200$ GeV (Au+Au), $\sqrt{s}=200$ GeV (p+p)",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=12.9,
            color=BODY_COLOR,
        )

    ymax = max(ax_raw.get_ylim()[1], ax_corr.get_ylim()[1])
    ax_raw.set_ylim(0, ymax * 1.03)
    ax_corr.set_ylim(0, ymax * 1.03)

    fig.savefig(OUT_PNG, dpi=160)
    plt.close(fig)

    manifest = {
        "slide": str(OUT_PNG),
        "source": "Same current-analysis objects used by slide09/slide14 ATLAS-vs-sPHENIX Fig.1 analogue.",
        "normalization": {
            "type": "area-normalized reconstructed shape",
            "formula": "density_i = counts_i / sum_j(counts_j * bin_width_j) over displayed x range",
            "x_range": [XMIN, XMAX],
            "reason": "shape-only overlay; not a final per-photon unfolded yield comparison",
        },
        "selection": {
            "photon_window": "strict full-bin mode for nominal 15<E_T<35 GeV; effective current bins are 16-35 GeV",
            "dphi": "|Delta phi| > 7pi/8",
            "jet_pt": fig1.jet_pt_label(),
        },
        "inputs": {
            "pp_data_file": str(fig1.PANELS[0].data_file),
            "pp_sideband_normalization": "PPG12 final-BDT leakage-corrected purity fit via project_pp_with_ppg12_purity_norm()",
            "auau_data_file": str(next(panel for panel in fig1.PANELS if panel.key == "auau_0_20").data_file),
            "auau_signal_mc_file": str(next(panel for panel in fig1.PANELS if panel.key == "auau_0_20").sim_file),
            "auau_correction": "fitted leakage-corrected purity, region-C sideband subtraction, then scaled combinatoric template subtraction",
        },
        "integrals_from_slide14_inputs": {
            "pp": inputs["pp"]["integrals"],
            "auau_0_20": inputs["auau_0_20"]["integrals"],
            "normalization_stats": stats,
        },
        "status": "Reconstructed-input diagnostic before unfolding; do not label as final unfolded per-photon result.",
    }
    OUT_MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n")
    OUT_SCRIPT.write_text(
        "This slide overlays p+p and Au+Au 0-20 percent using the same reconstructed xJ inputs shown on slide 14. "
        "The left panel is uncorrected region A, so it shows tight and isolated photon candidates before photon-ID sideband or combinatoric subtraction. "
        "The right panel is the corrected reconstructed input that would feed unfolding: the pp curve uses the PPG12-normalized region-C subtraction, while the Au+Au curve additionally subtracts the scaled embedded combinatoric template. "
        "Both panels are area-normalized over the visible x range, so this is a shape comparison, not a yield or final per-photon unfolded comparison.\n"
    )
    print(json.dumps({"slide": str(OUT_PNG), "manifest": str(OUT_MANIFEST), "speaker_script": str(OUT_SCRIPT)}, indent=2))


if __name__ == "__main__":
    main()
