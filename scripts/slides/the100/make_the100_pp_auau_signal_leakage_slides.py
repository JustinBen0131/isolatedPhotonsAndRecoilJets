#!/usr/bin/env python3
"""Build pp-versus-AuAu signal-leakage slide candidates, one per centrality.

The layout follows the IAN signal-leakage figure
(`figures/ppg12_reproduction/signal_leakage_inputs.tex`, rendered by
`scripts/plotting/pp_currentian/make_the97_corrected_si_candidate_signal_parity.py`):
one panel, three leakage coefficients drawn as black circles, red squares and
blue triangles, open markers for one sample and filled for the other, with the
sPHENIX annotation block upper left and the legend upper right.  The IAN's
current/PPG12 ratio subpanel is intentionally dropped.

Here the open/filled split carries pp versus Au+Au instead of PPG12 versus
current, and the inputs are exactly the points the THE-45 Region-A purity
comparison already derived, so the two slides share one provenance chain.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D


THIS_FILE = Path(__file__).resolve()
REPO = next((p for p in THIS_FILE.parents if (p / "AGENTS.md").exists()), THIS_FILE.parents[3])
POINTS_DIR = REPO / "dataOutput/slides/the45_jstg_20260720/the100_pp_regionA_purity_comparison"
OUTDIR = REPO / "dataOutput/slides/the45_jstg_20260720/the100_pp_auau_signal_leakage"

# Same coefficient styling as the IAN figure.
LEAKAGE = (
    ("fB", r"$c_B=B_{sig}/A_{sig}$", "black", "o"),
    ("fC", r"$c_C=C_{sig}/A_{sig}$", "#d62728", "s"),
    ("fD", r"$c_D=D_{sig}/A_{sig}$", "#1f77b4", "^"),
)
CENTRALITIES = (("0_20", "0–20%"), ("50_80", "50–80%"))
ET_LO, ET_HI = 15.0, 35.0
INK = "#142235"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_points(cent: str) -> tuple[list[dict[str, float]], list[dict[str, float]], Path]:
    """Return (pp rows, AuAu rows) for one centrality from the purity outputs."""
    source = POINTS_DIR / f"the100_pp_auau_{cent}_regionA_yield_and_purity_points.csv"
    if not source.is_file():
        raise FileNotFoundError(
            f"missing purity points CSV: {source}\n"
            "Run make_the100_pp_regionA_purity_comparison_slides.py first."
        )
    rows = list(csv.DictReader(source.open()))
    pp = [r for r in rows if r["system"] == "pp"]
    auau = [r for r in rows if r["system"] == "auau"]
    if not pp or not auau:
        raise RuntimeError(f"{source} is missing a pp or auau lane")
    for label, subset in (("auau", auau),):
        bad = {r["centrality"] for r in subset} - {cent}
        if bad:
            raise RuntimeError(f"{label} rows carry unexpected centralities {bad} in {source}")
    return pp, auau, source


def series(rows: list[dict[str, float]], key: str) -> tuple[np.ndarray, ...]:
    lo = np.array([float(r["pt_lo"]) for r in rows])
    hi = np.array([float(r["pt_hi"]) for r in rows])
    return (
        0.5 * (lo + hi),
        0.5 * (hi - lo),
        np.array([float(r[key]) for r in rows]),
        np.array([float(r[f"{key}_error"]) for r in rows]),
    )


def render(panels: list[tuple[str, str, list[dict], list[dict]]], output: Path) -> dict:
    """One slide, two stacked panels, sharing slide 16's grid geometry."""
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "dejavuserif",
            "axes.linewidth": 1.15,
        }
    )
    fig = plt.figure(figsize=(16, 9), dpi=160, facecolor="white")
    # Same split as the Region-A yield/purity slide.
    gs = fig.add_gridspec(
        # Same left/right/bottom/hspace/ratios as the Region-A yield/purity
        # slide; top pulled in from 0.780 to clear the second bullet, since
        # this slide carries two bullets above the panels rather than one.
        2, 1, left=0.090, right=0.972, top=0.726, bottom=0.088, hspace=0.28,
        height_ratios=(0.95, 1.05),
    )
    axes = [fig.add_subplot(gs[0]), fig.add_subplot(gs[1])]

    fig.text(
        0.052, 0.945,
        "Signal leakage into ABCD control regions: Au+Au and p+p",
        ha="left", va="top", fontsize=34.5, fontweight="bold", color=INK,
    )

    def bullet(y: float, label_text: str, body_text: str) -> None:
        """Blue arrowhead, bold 'Label:', then regular body on the same line."""
        fig.text(0.054, y, "\u25b6", ha="left", va="center", fontsize=15.0, color="#2468A8", fontfamily="DejaVu Sans")
        label = fig.text(0.075, y, label_text, fontsize=20.0, color=INK, ha="left", va="center", fontweight="bold")
        fig.canvas.draw()
        right = label.get_window_extent(renderer=fig.canvas.get_renderer()).x1
        body_x = fig.transFigure.inverted().transform((right, 0))[0] + 0.006
        fig.text(body_x, y, body_text, fontsize=20.0, color=INK, ha="left", va="center")

    # Evenly spaced between the slide title's ink bottom (y~146 px) and the top
    # panel's title (y~348 px): three equal ~41 px gaps around two 40 px rows.
    bullet(
        0.856,
        "Open markers:",
        "p+p; filled markers: Au+Au. Coefficients are signal-MC region ratios; each system keeps its native "
        r"$E_T$ bins.",
    )
    bullet(
        0.800,
        "Next step:",
        "tuning a non-tight sideband as is done in pp should reduce the region-C leakage discrepancy.",
    )

    # Draw both panels first so a common y-limit can be chosen.
    top_value = 0.0
    for ax, (cent, cent_label, pp, auau) in zip(axes, panels):
        for key, _, color, marker in LEAKAGE:
            for rows, filled, shift in ((pp, False, -0.13), (auau, True, +0.13)):
                x, xe, y, ye = series(rows, key)
                ax.errorbar(
                    x + shift, y, xerr=xe, yerr=ye, fmt=marker, ms=8.0,
                    mfc=color if filled else "white", mec=color, ecolor=color, color=color,
                    capsize=0, linestyle="none", elinewidth=1.2,
                )
                finite = y[np.isfinite(y)] + ye[np.isfinite(y)]
                if finite.size:
                    top_value = max(top_value, float(np.max(finite)))

    y_max = max(0.15, 1.75 * top_value)
    for index, (ax, (cent, cent_label, _, _)) in enumerate(zip(axes, panels)):
        ax.set_xlim(ET_LO, ET_HI)
        ax.set_ylim(0.0, y_max)
        ax.set_ylabel("Signal leakage", fontsize=18.0)
        ax.set_title(f"Au+Au {cent_label} vs p+p", fontsize=19.0, fontweight="bold", pad=8, color=INK)
        ax.tick_params(axis="both", labelsize=15.0, direction="in", top=True, right=True, length=6)
        ax.grid(True, color="#D7E0EA", linewidth=0.8, alpha=0.75)
        ax.set_axisbelow(True)
        if index == 0:
            ax.tick_params(labelbottom=False)
        else:
            ax.set_xlabel(r"$E_T^{\gamma,\mathrm{rec}}$ [GeV]", fontsize=19.0, labelpad=8)

    # Annotation block and sample composition live on the top panel.
    sphenix = axes[0].text(
        0.020, 0.955, "sPHENIX", transform=axes[0].transAxes, va="top", ha="left",
        fontsize=18, fontweight="bold", fontstyle="italic",
    )
    fig.canvas.draw()
    sphenix_right = sphenix.get_window_extent(renderer=fig.canvas.get_renderer()).x1
    internal_x = axes[0].transAxes.inverted().transform((sphenix_right, 0))[0] + 0.008
    axes[0].text(internal_x, 0.955, "Internal", transform=axes[0].transAxes, va="top", ha="left", fontsize=18)
    axes[0].text(0.020, 0.845, r"$\sqrt{s_{NN}}=200\ \mathrm{GeV}$,  $|\eta^\gamma|<0.7$", transform=axes[0].transAxes, va="top", ha="left", fontsize=16)
    axes[0].text(0.320, 0.955, "p+p: photon+jet 5/10/20 Pythia", transform=axes[0].transAxes, va="top", ha="left", fontsize=15.0, color=INK)
    axes[0].text(0.320, 0.845, "Au+Au: photon+jet 12/20 embedded Pythia", transform=axes[0].transAxes, va="top", ha="left", fontsize=15.0, color=INK)

    # Legend on the lower panel, which carries no annotation block.
    handles = []
    for key, text, color, marker in LEAKAGE:
        handles.append(Line2D([0], [0], marker=marker, color=color, mfc="white", mec=color, linestyle="none", ms=8.0, label=f"p+p {text}"))
    for key, text, color, marker in LEAKAGE:
        handles.append(Line2D([0], [0], marker=marker, color=color, mfc=color, mec=color, linestyle="none", ms=8.0, label=f"Au+Au {text}"))
    # Upper left: the leakage curves all sit below ~0.3, and the right side
    # carries the widest bins, so the left corner is the clear one.
    axes[1].legend(
        handles=handles, loc="upper left", bbox_to_anchor=(0.015, 0.995), ncol=2,
        frameon=False, fontsize=15.5, columnspacing=1.5, handletextpad=0.5, labelspacing=0.45,
    )

    fig.savefig(output, dpi=160, facecolor="white")
    plt.close(fig)
    return {"max_drawn_leakage": top_value, "y_max": y_max}


def main() -> int:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    panels = []
    sources = {}
    for cent, cent_label in CENTRALITIES:
        pp, auau, source = load_points(cent)
        panels.append((cent, cent_label, pp, auau))
        sources[cent] = {"source_points_csv": str(source), "source_points_sha256": sha256(source)}

    png = OUTDIR / "the100_pp_auau_signal_leakage_two_panel_slide.png"
    stats = render(panels, png)

    manifest = OUTDIR / "the100_pp_auau_signal_leakage_manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "schema": "THE45_THE100_PP_AUAU_SIGNAL_LEAKAGE_SLIDES_V2",
                "layout_reference": "IAN figures/ppg12_reproduction/signal_leakage_inputs.tex; ratio subpanel omitted; two-panel grid matches the THE-45 Region-A yield/purity slide",
                "marker_contract": "open = p+p, filled = Au+Au; black circle c_B, red square c_C, blue triangle c_D",
                "coefficient_definition": "c_R = N_R^sig / N_A^sig from the signal-MC ABCD histograms",
                "panels": {"top": "Au+Au 0-20% vs p+p", "bottom": "Au+Au 50-80% vs p+p"},
                "shared_y_limit": stats["y_max"],
                "et_range_GeV": [ET_LO, ET_HI],
                "binning_note": "pp and Au+Au keep their native, mutually offset E_T bins; only 26 and 35 GeV are shared edges",
                "upstream": "points derived by make_the100_pp_regionA_purity_comparison_slides.py from THE-100 complement data/signal and the current pp artifacts",
                "output": {"png": str(png), "png_sha256": sha256(png), **stats},
                "sources": sources,
            },
            indent=2,
        )
        + "\n"
    )
    print(json.dumps({"png": str(png), "manifest": str(manifest), **stats}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
