#!/usr/bin/env python3
"""Clean THE-85 unfolded xJgamma overlay: Au+Au 0-20% vs p+p.

This intentionally keeps the slide to one plot.  It reads the first-pass
unfolded arrays already produced by make_unfolded_xjgamma_1x3.py and makes a
presentation-facing overlay without an I_AA panel.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict

import matplotlib.pyplot as plt
import numpy as np


REPO = Path(__file__).resolve().parents[3]
IN_DIR = REPO / "dataOutput/the85_auau_xjgamma_unfolding_push/unfolded_xjgamma_firstpass"
OUT_PNG = IN_DIR / "slide06_clean_overlay_auau020_vs_pp_noiaa.png"
OUT_MANIFEST = IN_DIR / "slide06_clean_overlay_auau020_vs_pp_noiaa_manifest.json"
OUT_SCRIPT = IN_DIR / "slide06_clean_overlay_auau020_vs_pp_noiaa_speaker_script.md"


def load_npz(name: str) -> Dict[str, np.ndarray]:
    path = IN_DIR / name
    if not path.exists():
        raise FileNotFoundError(path)
    z = np.load(path)
    return {k: z[k] for k in z.files}


def draw_curve(ax, curve: Dict[str, np.ndarray], *, label: str, color: str, marker: str, open_marker: bool) -> None:
    x = curve["x_centers"]
    edges = curve["x_edges"]
    y = curve["y"]
    ey = curve["ey"]
    xerr = 0.5 * np.diff(edges)
    mask = np.isfinite(y) & np.isfinite(ey) & (x >= 0.18) & (x <= 1.42)
    ax.errorbar(
        x[mask],
        y[mask],
        xerr=xerr[mask],
        yerr=ey[mask],
        fmt=marker,
        ms=9.0,
        lw=1.55,
        elinewidth=1.45,
        capsize=3.2,
        color=color,
        mfc="white" if open_marker else color,
        mec=color,
        mew=1.8,
        label=label,
        zorder=4 if not open_marker else 3,
    )


def main() -> None:
    auau = load_npz("the85_unfolded_xjgamma_auau_0_20_firstpass.npz")
    pp = load_npz("the85_unfolded_xjgamma_pp_basev3e_firstpass.npz")

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "dejavuserif",
            "axes.linewidth": 1.35,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.major.size": 7,
            "ytick.major.size": 7,
            "xtick.minor.size": 4,
            "ytick.minor.size": 4,
        }
    )

    title_color = "#121827"
    body_color = "#364153"
    blue = "#1f77b4"
    red = "#d62728"

    fig = plt.figure(figsize=(16, 9), dpi=160)
    fig.patch.set_facecolor("white")

    fig.text(
        0.075,
        0.93,
        r"First-pass unfolded $x_{J\gamma}$: central Au+Au vs p+p",
        ha="left",
        va="top",
        fontsize=33,
        fontweight="bold",
        color=title_color,
    )
    fig.text(
        0.075,
        0.868,
        r"$15 < E_T^\gamma < 35$ GeV, anti-$k_T$ $R=0.4$, $|\Delta\phi|>7\pi/8$",
        ha="left",
        va="top",
        fontsize=20,
        color=body_color,
    )

    ax = fig.add_axes([0.105, 0.17, 0.83, 0.62])
    ax.axhline(0.0, color="#6f7782", lw=1.0)
    draw_curve(
        ax,
        auau,
        label=r"Au+Au 0-20% corrected",
        color=blue,
        marker="o",
        open_marker=False,
    )
    draw_curve(
        ax,
        pp,
        label=r"p+p baseV3E reference",
        color=red,
        marker="s",
        open_marker=True,
    )

    ax.set_xlim(0.0, 1.45)
    ax.set_ylim(-0.08, 1.18)
    ax.set_xlabel(r"$x_{J\gamma}=p_T^{jet}/p_T^\gamma$", fontsize=23)
    ax.set_ylabel(r"$(1/N_\gamma)\,dN/dx_{J\gamma}$", fontsize=23)
    ax.grid(True, color="#e2e7ee", lw=0.85)
    ax.minorticks_on()
    ax.tick_params(labelsize=17, top=True, right=True)
    ax.legend(
        loc="upper right",
        frameon=False,
        fontsize=17.5,
        handlelength=2.1,
        borderpad=0.2,
        labelspacing=0.7,
    )

    fig.savefig(OUT_PNG, dpi=160)
    plt.close(fig)

    manifest = {
        "slide": str(OUT_PNG),
        "source_npz": {
            "auau_0_20": str(IN_DIR / "the85_unfolded_xjgamma_auau_0_20_firstpass.npz"),
            "pp_basev3e": str(IN_DIR / "the85_unfolded_xjgamma_pp_basev3e_firstpass.npz"),
        },
        "selection": {
            "photon_et": [15, 35],
            "jet": "anti-kT R=0.4",
            "dphi": "|Delta phi| > 7pi/8",
            "auau": "default 14-feature AuAu BDT, WP80, ABCD input, embedded combinatoric subtraction",
            "pp": "baseV3E p+p reference; current overlay uses first-pass pp diagnostic arrays",
        },
        "supersedes_for_review": "slide06_quenching_signature_auau020_5080_vs_pp_candidate.png",
        "caveat": "No I_AA panel. This is the clean overlay requested by Justin after rejecting the 50-80/ratio-heavy slide.",
    }
    OUT_MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n")

    OUT_SCRIPT.write_text(
        "This slide should be presented simply: it overlays the first-pass unfolded central Au+Au xJgamma distribution with the p+p reference using the same photon ET window and 7pi/8 recoil selection.\n\n"
        "The point of the slide is visual clarity, not a final modification-factor claim. The central Au+Au curve sits below the p+p reference through the balanced-recoil region, which is the quenching-sensitive pattern we are testing.\n\n"
        "Do not overstate it yet. The current p+p comparison is still a first-pass reference, and final pp ABCD cleanup plus Au+Au combinatoric-template closure remain the next validation steps.\n"
    )

    print(json.dumps({"slide": str(OUT_PNG), "manifest": str(OUT_MANIFEST), "speaker_script": str(OUT_SCRIPT)}, indent=2))


if __name__ == "__main__":
    main()
