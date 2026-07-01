#!/usr/bin/env python3
"""Build a slide-ready AuAu isolation-stack diagnostic from current THE-69 outputs."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import patches

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "scripts/plotting/pp_currentian"))

from make_ppg12_style_isolation_stack import (  # noqa: E402
    DEFAULT_AUAU_DATA,
    DEFAULT_AUAU_SIGNAL,
    DEFAULT_PP_DATA,
    DEFAULT_PP_SIGNAL,
    SampleSpec,
    draw_step_band,
    old_ppg12_stack,
)


OUTDIR = REPO / "dataOutput/the85_auau_xjgamma_unfolding_push/isolation_stack_ppg12_style_20260628"
PNG = OUTDIR / "slide_isolation_stack_pp16to22_auau16to35_ppg12math_sourcecheck.png"
SCRIPT = OUTDIR / "slide_isolation_stack_pp16to22_auau16to35_ppg12math_sourcecheck_speaker_script.md"
MANIFEST = OUTDIR / "slide_isolation_stack_pp16to22_auau16to35_ppg12math_sourcecheck_manifest.json"


def auau_spec(cent: str) -> SampleSpec:
    cent_label = cent.replace("_", "-")
    return SampleSpec(
        tag=f"auau_default_bdt_cent_{cent}_clean_slide",
        collision_label=r"$\mathrm{Au+Au}\ \sqrt{s_{NN}}=200\ \mathrm{GeV}$",
        pt_label=r"$16 < E_T^\gamma < 35\ \mathrm{GeV}$",
        pt_bins=("16_18", "18_20", "20_22", "22_24", "24_26", "26_35"),
        data_file=DEFAULT_AUAU_DATA,
        signal_file=DEFAULT_AUAU_SIGNAL,
        data_dir="MBD_NS_geq_2_vtx_lt_150",
        signal_dir="SIM",
        tight_pattern=f"h_Eiso_tight_isoR40_pT_{{pt}}_cent_{cent}",
        nontight_pattern=f"h_Eiso_nonTight_isoR40_pT_{{pt}}_cent_{cent}",
        signal_pattern=f"h_Eiso_tight_isoR40_pT_{{pt}}_cent_{cent}",
        output_name=f"slide_auau_isolation_stack_16to35_cent{cent}_clean.png",
        sideband_note=f"Current AuAu non-tight is BDT complement; centrality {cent_label}%.",
    )


def pp_spec() -> SampleSpec:
    return SampleSpec(
        tag="pp_current_baseline_clean_slide",
        collision_label=r"$p{+}p\ \sqrt{s}=200\ \mathrm{GeV}$",
        pt_label=r"$16 < E_T^\gamma < 22\ \mathrm{GeV}$",
        pt_bins=("16_18", "18_20", "20_22"),
        data_file=DEFAULT_PP_DATA,
        signal_file=DEFAULT_PP_SIGNAL,
        data_dir="Photon_4_GeV_plus_MBD_NS_geq_1",
        signal_dir="SIM",
        tight_pattern="h_Eiso_tight_pT_{pt}",
        nontight_pattern="h_Eiso_nonTight_pT_{pt}",
        signal_pattern="h_Eiso_tight_pT_{pt}",
        output_name="slide_pp_isolation_stack_16to22_clean.png",
        sideband_note="Current pp ROOT has RecoilJets h_Eiso_nonTight, not PPG12 h_nontight_isoET; normalization follows PPG12 Fig. 3.",
    )


def add_arrow_bullet(fig: plt.Figure, x: float, y: float, text: str, *, fontsize: int = 24) -> None:
    tri = patches.RegularPolygon(
        (x, y + 0.003),
        numVertices=3,
        radius=0.008,
        orientation=-np.pi / 2,
        transform=fig.transFigure,
        facecolor="#1f4e79",
        edgecolor="none",
    )
    fig.patches.append(tri)
    fig.text(x + 0.018, y, text, ha="left", va="center", fontsize=fontsize, color="#111827")


def draw_panel(ax: plt.Axes, spec: SampleSpec, result: dict[str, object], cent_title: str) -> None:
    edges = result["edges"]
    tight = result["tight"]
    tight_err = result["tight_err"]
    bkg = result["background"]
    sig = result["signal_mc"]
    widths = np.diff(edges)
    centers = edges[:-1] + 0.5 * widths

    draw_step_band(ax, edges, np.zeros_like(bkg), bkg, "#f2aaa6", "#8a403c", "Scaled non-tight data", alpha=0.72, linewidth=1.0)
    draw_step_band(ax, edges, bkg, bkg + sig, "#999bf2", "#2f327d", "Signal MC stacked", alpha=0.64, linewidth=1.0)
    ax.errorbar(
        centers,
        tight,
        yerr=tight_err,
        fmt="o",
        color="black",
        ecolor="black",
        elinewidth=1.15,
        capsize=0,
        markersize=4.2,
        label="Tight data",
        zorder=5,
    )

    ymax = max(float(np.nanmax(tight + tight_err)), float(np.nanmax(bkg + sig)))
    ax.set_xlim(-2.0, 15.0)
    ax.set_xticks([-2, 0, 5, 10, 15])
    ax.set_ylim(0, ymax * 1.23 if ymax > 0 else 1.0)
    ax.set_title(cent_title, fontsize=18.5, fontweight="bold", pad=4)
    ax.set_xlabel(r"$E_T^{\mathrm{iso,reco}}\ \mathrm{[GeV]}$", fontsize=15)
    ax.set_ylabel("")
    ax.tick_params(axis="both", which="both", direction="in", top=True, right=True, labelsize=12.0, length=6)
    ax.tick_params(axis="both", which="minor", length=3.5)
    ax.minorticks_on()

    ax.text(0.94, 0.92, r"$\it{\bf{sPHENIX}}$ Internal", transform=ax.transAxes, ha="right", va="top", fontsize=12.2)
    ax.text(0.94, 0.835, spec.collision_label, transform=ax.transAxes, ha="right", va="top", fontsize=10.9)
    ax.text(0.94, 0.755, spec.pt_label, transform=ax.transAxes, ha="right", va="top", fontsize=10.9)


def main() -> int:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "dejavuserif",
            "axes.linewidth": 1.4,
            "xtick.major.width": 1.35,
            "ytick.major.width": 1.35,
            "xtick.minor.width": 1.0,
            "ytick.minor.width": 1.0,
        }
    )

    panel_specs: list[tuple[SampleSpec, str]] = [(pp_spec(), r"$p{+}p$")]
    panel_specs += [(auau_spec(cent), title) for cent, title in [("0_20", "Au+Au 0-20%"), ("20_50", "Au+Au 20-50%"), ("50_80", "Au+Au 50-80%")]]
    results = []
    for spec, title in panel_specs:
        if title == r"$p{+}p$":
            low_group, high_group = 1, 5
        else:
            low_group, high_group = 5, 10
        result = old_ppg12_stack(spec, tail_low=6.0, low_region_group=low_group, high_region_group=high_group)
        result["metadata"]["display_rebin_note"] = f"{low_group}/{high_group}"
        results.append((spec, result, title))

    fig = plt.figure(figsize=(12.8, 7.2), dpi=200)
    fig.patch.set_facecolor("white")
    fig.text(
        0.055,
        0.945,
        "Isolation-energy template comparison",
        ha="left",
        va="top",
        fontsize=30,
        fontweight="bold",
        color="#111827",
    )
    add_arrow_bullet(
        fig,
        0.062,
        0.865,
        r"Same PPG12 math: variable rebin, divide by bin width, then tail-match non-tight above 6 GeV.",
        fontsize=19.0,
    )
    add_arrow_bullet(
        fig,
        0.062,
        0.815,
        r"pp uses the 16-22 GeV Fig.-3-like window; Au+Au keeps the current 16-35 GeV baseline output.",
        fontsize=19.0,
    )

    fig.text(
        0.062,
        0.755,
        r"$|\eta^\gamma|<0.7$; Au+Au red is the broad BDT-complement sideband by construction.",
        ha="left",
        va="center",
        fontsize=18.0,
        color="#4b5563",
    )

    legend_y = 0.695
    legend_box = patches.FancyBboxPatch(
        (0.055, legend_y - 0.021),
        0.515,
        0.045,
        boxstyle="round,pad=0.006,rounding_size=0.006",
        transform=fig.transFigure,
        facecolor="#f8fafc",
        edgecolor="#d6dee8",
        linewidth=0.8,
    )
    fig.patches.append(legend_box)
    fig.text(0.070, legend_y, "components", ha="left", va="center", fontsize=12.8, color="#4b5563", fontweight="bold")
    fig.patches.append(patches.Rectangle((0.150, legend_y - 0.010), 0.019, 0.020, transform=fig.transFigure, facecolor="#f2aaa6", edgecolor="#8a403c", linewidth=0.9, alpha=0.74))
    fig.text(0.174, legend_y, "scaled non-tight data", ha="left", va="center", fontsize=14.3, color="#111827")
    fig.patches.append(patches.Rectangle((0.325, legend_y - 0.010), 0.019, 0.020, transform=fig.transFigure, facecolor="#8f91f0", edgecolor="#2f327d", linewidth=0.9, alpha=0.68))
    fig.text(0.349, legend_y, "signal MC", ha="left", va="center", fontsize=14.3, color="#111827")
    fig.text(0.433, legend_y, r"$\bullet$", ha="center", va="center", fontsize=18, color="black")
    fig.text(0.449, legend_y, "tight data", ha="left", va="center", fontsize=14.3, color="#111827")

    fig.text(0.019, 0.395, "Counts / Bin Width", ha="center", va="center", rotation=90, fontsize=17, color="#111827")
    grid = fig.add_gridspec(1, 4, left=0.062, right=0.974, bottom=0.135, top=0.625, wspace=0.23)
    axes = [fig.add_subplot(grid[0, i]) for i in range(4)]
    for ax, (spec, result, title) in zip(axes, results):
        draw_panel(ax, spec, result, title)

    fig.savefig(PNG)
    plt.close(fig)

    manifest = {
        "output_png": str(PNG),
        "data_file": str(DEFAULT_AUAU_DATA),
        "signal_file": str(DEFAULT_AUAU_SIGNAL),
        "pp_data_file": str(DEFAULT_PP_DATA),
        "pp_signal_file": str(DEFAULT_PP_SIGNAL),
        "pt_window": "pp: 16 < E_T^gamma < 22 GeV; AuAu: 16 < E_T^gamma < 35 GeV",
        "panels": [title for _, title in panel_specs],
        "normalization": "PPG12 Figure 3 style: non-tight data scaled to tight-data tail above Eiso=6 GeV after bin-width scaling; signal MC scaled to tight minus scaled non-tight integral",
        "render_rebin": "Hybrid display rebin: pp uses PPG12 1/5 source-bin grouping; AuAu uses 5/10 grouping for current statistics",
        "pp_source_caveat": "Current pp input lacks PPG12 h_tight_isoET/h_nontight_isoET. It contains RecoilJets h_Eiso_tight/h_Eiso_nonTight only, so the normalization math is PPG12-like but the non-tight source population is not proven identical to PPG12 Fig. 3.",
        "auau_source_caveat": "Current AuAu non-tight row is nonTightAuAuBDTComplement, a broad BDT-complement sideband; this can make the scaled red template sit higher or more signal-like than the bounded PPG12 non-tight sideband.",
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n")
    SCRIPT.write_text(
        "\n".join(
            [
                "# Speaker Script",
                "",
                "This slide keeps the PPG12 Figure 3 normalization logic but uses coarser display bins so the current lower-stat Au+Au panels are readable.",
                "The black points are tight selected data candidates. The red component is the scaled non-tight data sideband. The blue component is signal MC stacked on top of the red component, following the paper-style drawing convention.",
                "The key pp caveat is that the current pp ROOT does not contain PPG12's h_tight_isoET/h_nontight_isoET histograms. It contains RecoilJets h_Eiso_tight/h_Eiso_nonTight, so the normalization recipe matches PPG12 but the sideband source population is not the same object.",
                "The key Au+Au caveat is that the current non-tight row is the broad BDT-complement sideband. PPG12 used a bounded non-tight sideband, so a red template that sits high can reflect sideband definition, not just plotting.",
                "",
            ]
        )
    )
    print(PNG)
    print(SCRIPT)
    print(MANIFEST)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
