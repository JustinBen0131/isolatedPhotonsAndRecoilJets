#!/usr/bin/env python3
"""Build slide 48 with the pp isolation panel repaired to the final PPG12 pp contract."""

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
    SampleSpec,
    draw_step_band,
    old_ppg12_stack,
)


CAMPAIGN = "the76_ppg12_parity_full_20260701_003024"
CFG = "jetMinPtScan_dphiScan_vz60_isoR40_isSliding_preselectionNewPPG12_tightNewPPG12_nonTightNewPPG12"
PP_DATA = REPO / "dataOutput/current_recoiljets_artifacts/current/pp_data_merged/current.root"
PP_SIGNAL = (
    REPO
    / "dataOutput/ppg12Parity"
    / CAMPAIGN
    / "final_roots/photonjet"
    / CFG
    / "photonJet5and10and20merged_SIM/RecoilJets_photonjet5plus10plus20_MERGED.root"
)

OUTDIR = REPO / "dataOutput/ppg12Parity" / CAMPAIGN / "slide48_isolation_stack"
PNG = OUTDIR / "slide48_iso_template_comparison_ppfixed_finalpp_photonjet_sim.png"
MANIFEST = OUTDIR / "slide48_iso_template_comparison_ppfixed_finalpp_photonjet_sim_manifest.json"
SCRIPT = OUTDIR / "slide48_iso_template_comparison_ppfixed_finalpp_photonjet_sim_speaker_script.md"


def pp_spec() -> SampleSpec:
    return SampleSpec(
        tag="pp_final_hierarchical_scaledtrigger30_fig3_contract",
        collision_label=r"$p{+}p\ \sqrt{s}=200\ \mathrm{GeV}$",
        pt_label=r"$16 < E_T^\gamma < 22\ \mathrm{GeV}$",
        pt_bins=("3", "4", "5"),
        data_file=PP_DATA,
        signal_file=PP_SIGNAL,
        data_dir="PPG12_scaledtrigger30",
        signal_dir="SIM",
        tight_pattern="h_tight_isoET_0_{pt}",
        nontight_pattern="h_nontight_isoET_0_{pt}",
        signal_pattern="h_tight_isoET_0_{pt}",
        output_name="slide48_pp_fixed_finalpp.png",
        sideband_note=(
            "Final hierarchical July pp data. PPG12 Figure 3 contract: scaledtrigger30, "
            "h_tight_isoET_0_{3,4,5}, h_nontight_isoET_0_{3,4,5}; photon+jet SIM uses "
            "SIM/h_tight_isoET_0_{3,4,5}."
        ),
    )


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
        output_name=f"slide48_auau_isolation_stack_cent{cent}.png",
        sideband_note=f"Current AuAu non-tight is BDT complement; centrality {cent_label}%.",
    )


def add_arrow_bullet(fig: plt.Figure, x: float, y: float, text: str, *, fontsize: int = 20) -> None:
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


def draw_panel(ax: plt.Axes, spec: SampleSpec, result: dict[str, object], panel_title: str) -> None:
    edges = result["edges"]
    tight = result["tight"]
    tight_err = result["tight_err"]
    bkg = result["background"]
    sig = result["signal_mc"]
    widths = np.diff(edges)
    centers = edges[:-1] + 0.5 * widths

    draw_step_band(
        ax,
        edges,
        np.zeros_like(bkg),
        bkg,
        "#f2aaa6",
        "#8a403c",
        "scaled non-tight data",
        alpha=0.72,
        linewidth=1.0,
    )
    draw_step_band(
        ax,
        edges,
        bkg,
        bkg + sig,
        "#999bf2",
        "#2f327d",
        "signal MC",
        alpha=0.64,
        linewidth=1.0,
    )
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
        label="tight data",
        zorder=5,
    )

    ymax = max(float(np.nanmax(tight + tight_err)), float(np.nanmax(bkg + sig)))
    ax.set_xlim(-2.0, 15.0)
    ax.set_xticks([-2, 0, 5, 10, 15])
    ax.set_ylim(0, ymax * 1.23 if ymax > 0 else 1.0)
    ax.set_title(panel_title, fontsize=18.5, fontweight="bold", pad=4)
    ax.set_xlabel(r"$E_T^{\mathrm{iso,reco}}\ \mathrm{[GeV]}$", fontsize=15)
    ax.tick_params(axis="both", which="both", direction="in", top=True, right=True, labelsize=12.0, length=6)
    ax.tick_params(axis="both", which="minor", length=3.5)
    ax.minorticks_on()

    ax.text(0.94, 0.92, r"$\it{\bf{sPHENIX}}$ Internal", transform=ax.transAxes, ha="right", va="top", fontsize=12.2)
    ax.text(0.94, 0.835, spec.collision_label, transform=ax.transAxes, ha="right", va="top", fontsize=10.9)
    ax.text(0.94, 0.755, spec.pt_label, transform=ax.transAxes, ha="right", va="top", fontsize=10.9)


def write_speaker_script() -> None:
    SCRIPT.write_text(
        "\n".join(
            [
                "# Speaker Script",
                "",
                "This is the repaired slide-48 comparison. The pp panel now uses the final July pp data product and the current photon+jet merged SIM, not the older June pp inputs.",
                "For pp, the object contract matches PPG12 Figure 3: scaledtrigger30 data, pt bins 16-22 GeV via indices 3, 4, and 5, h_tight_isoET for the tight data, h_nontight_isoET for the red sideband, and photon+jet SIM h_tight_isoET for the blue signal component.",
                "The normalization follows Shuhang's CONF_plots.C: variable rebin, divide by bin width, scale non-tight to the tight tail above 6 GeV, then scale signal MC to the tight-minus-background integral.",
                "The AuAu panels are unchanged from the earlier slide-48 diagnostic; their red component is still the current broad BDT-complement sideband, so this slide is a template-shape comparison rather than an exact PPG12 sideband identity claim for AuAu.",
                "",
            ]
        )
    )


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

    panel_specs: list[tuple[SampleSpec, str, int, int]] = [
        (pp_spec(), r"$p{+}p$", 1, 5),
        (auau_spec("0_20"), "Au+Au 0-20%", 5, 10),
        (auau_spec("20_50"), "Au+Au 20-50%", 5, 10),
        (auau_spec("50_80"), "Au+Au 50-80%", 5, 10),
    ]
    results = []
    for spec, title, low_group, high_group in panel_specs:
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
        r"pp panel uses final July pp data plus current photon+jet SIM under the PPG12 Fig. 3 object contract.",
        fontsize=18.5,
    )
    add_arrow_bullet(
        fig,
        0.062,
        0.815,
        r"Same PPG12 math: variable rebin, bin-width scaling, tail-match above 6 GeV, then signal normalization.",
        fontsize=18.5,
    )
    fig.text(
        0.062,
        0.755,
        r"$|\eta^\gamma|<0.7$; pp uses 16-22 GeV, Au+Au keeps the current 16-35 GeV baseline output.",
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
    fig.patches.append(
        patches.Rectangle(
            (0.150, legend_y - 0.010),
            0.019,
            0.020,
            transform=fig.transFigure,
            facecolor="#f2aaa6",
            edgecolor="#8a403c",
            linewidth=0.9,
            alpha=0.74,
        )
    )
    fig.text(0.174, legend_y, "scaled non-tight data", ha="left", va="center", fontsize=14.3, color="#111827")
    fig.patches.append(
        patches.Rectangle(
            (0.325, legend_y - 0.010),
            0.019,
            0.020,
            transform=fig.transFigure,
            facecolor="#8f91f0",
            edgecolor="#2f327d",
            linewidth=0.9,
            alpha=0.68,
        )
    )
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
        "speaker_script": str(SCRIPT),
        "pp_data_file": str(PP_DATA),
        "pp_signal_file": str(PP_SIGNAL),
        "auau_data_file": str(DEFAULT_AUAU_DATA),
        "auau_signal_file": str(DEFAULT_AUAU_SIGNAL),
        "pp_reference_code": "ppg12codeGit/plotting/CONF_plots.C",
        "pp_pt_indices": [3, 4, 5],
        "normalization": "CONF_plots.C: variable rebin, Scale(width), normalize non-tight to tight tail above Eiso=6 GeV, normalize signal MC to tight-minus-background integral.",
        "panels": [title for _, _, title in results],
        "panel_metadata": [result["metadata"] for _, result, _ in results],
        "caveat": "AuAu sideband remains the current broad BDT-complement sideband from the earlier slide-48 diagnostic; only pp was repaired to the final PPG12 pp data/SIM contract.",
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n")
    write_speaker_script()

    print(PNG)
    print(MANIFEST)
    print(SCRIPT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
