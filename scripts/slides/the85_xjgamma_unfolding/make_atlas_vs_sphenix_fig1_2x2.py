#!/usr/bin/env python3
"""2x2 ATLAS-vs-sPHENIX reconstructed-xJ subtraction-input slide."""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.offsetbox import AnnotationBbox, HPacker, TextArea
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from PIL import Image, ImageChops


REPO = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import make_atlas_fig1_reco_background_breakdown as fig1  # noqa: E402
import make_pp_atlas_vs_sphenix_fig1_comparison as ppcomp  # noqa: E402


OUT_DIR = REPO / "dataOutput/the85_auau_xjgamma_unfolding_push/unfolded_xjgamma_firstpass"
ATLAS_DIR = REPO / "dataOutput/the85_auau_xjgamma_unfolding_push/atlas_reference"
ATLAS_PP = ATLAS_DIR / "atlas_fig1_pp_panel_crop_hi.png"
ATLAS_010 = ATLAS_DIR / "atlas_fig1_pbpb_0_10_panel_crop_hi.png"
ATLAS_PP_TIGHT = ATLAS_DIR / "atlas_fig1_pp_panel_crop_hi_tight.png"
ATLAS_010_TIGHT = ATLAS_DIR / "atlas_fig1_pbpb_0_10_panel_crop_hi_tight.png"
FORCED_OUT_SUFFIX = os.environ.get("THE85_OUT_SUFFIX", fig1.OUT_SUFFIX).strip()
OUT_PNG = OUT_DIR / f"slide09_atlas_vs_sphenix_fig1_2x2_reco_inputs{FORCED_OUT_SUFFIX}.png"
OUT_MANIFEST = OUT_DIR / f"slide09_atlas_vs_sphenix_fig1_2x2_reco_inputs{FORCED_OUT_SUFFIX}_manifest.json"
OUT_SCRIPT = OUT_DIR / f"slide09_atlas_vs_sphenix_fig1_2x2_reco_inputs{FORCED_OUT_SUFFIX}_speaker_script.md"


def ensure_atlas_010_crop() -> None:
    """Write the clean ATLAS Pb+Pb 0-10 crop from the high-res page render."""
    from PIL import Image

    full = ATLAS_DIR / "atlas_page8_figure1_full_6x.png"
    if not full.exists():
        raise FileNotFoundError(full)
    im = Image.open(full).convert("RGB")
    im.crop((2284, 470, 2284 + 1040, 470 + 1010)).save(ATLAS_010)


def write_tight_crop(src: Path, dst: Path) -> None:
    im = Image.open(src).convert("RGB")
    bg = Image.new("RGB", im.size, (255, 255, 255))
    diff = ImageChops.difference(im, bg).convert("L")
    mask = diff.point(lambda value: 255 if value > 8 else 0)
    box = mask.getbbox()
    if not box:
        im.save(dst)
        return
    pad = 4
    box = (
        max(0, box[0] - pad),
        max(0, box[1] - pad),
        min(im.size[0], box[2] + pad),
        min(im.size[1], box[3] + pad),
    )
    im.crop(box).save(dst)


def style_axis(ax, ylabel: bool = False, xlabel: bool = True) -> None:
    ax.set_xlim(0.2, 1.85)
    if xlabel:
        ax.set_xlabel(r"Reconstructed $x_{J\gamma}$", fontsize=15.8, labelpad=5)
    else:
        ax.set_xlabel("")
    if ylabel:
        ax.set_ylabel("Entries", fontsize=15.8, labelpad=5)
    ax.minorticks_on()
    ax.tick_params(labelsize=12.7, top=True, right=True, direction="in", length=6, labelbottom=xlabel)
    ax.tick_params(which="minor", top=True, right=True, direction="in", length=3)
    for spine in ax.spines.values():
        spine.set_linewidth(1.2)


def draw_component_panel(
    ax,
    result: dict,
    pp_panel: bool = False,
    xlabel: bool = True,
    canvas_label: str = "",
    selection_label: str = "",
) -> None:
    xedges = result["x_edges"]
    centers = result["x_centers"]
    widths = np.diff(xedges)
    mask = (centers >= 0.2) & (centers <= 1.85)
    raw = result["raw"]
    side = result["sideband"]
    comb = result["comb"]
    bkg = result["bkg_sub"]
    bkg_err = result["bkg_sub_err"]

    ax.stairs(raw, xedges, color="#8f8f8f", lw=2.2, label="Raw Data")
    if pp_panel:
        ax.plot(
            [xedges[0], xedges[-1]],
            [0.0, 0.0],
            color="#d62728",
            lw=2.0,
            linestyle=(0, (1, 1)),
            label="Comb. Bkg.",
        )
    elif float(np.sum(comb)) > 0:
        ax.stairs(comb, xedges, color="#d62728", lw=2.0, linestyle=(0, (1, 1)), label="Comb. Bkg.")
    if float(np.sum(side)) > 0:
        ax.stairs(side, xedges, color="#1f4cff", lw=2.1, linestyle=(0, (2, 2)), label="Dijet Bkg. / photon-ID sideband")
    ax.errorbar(
        centers[mask],
        bkg[mask],
        xerr=0.5 * widths[mask],
        yerr=bkg_err[mask],
        fmt="o",
        color="black",
        ms=4.5,
        elinewidth=1.0,
        capsize=0,
        label="Bkg-Sub Data",
        zorder=5,
    )
    # The sPHENIX panels carry in-canvas selection text.  Keep that TLatex-style
    # block in the upper-left open canvas, but reserve enough vertical headroom
    # that it does not sit on the high raw spectrum in the central Au+Au pad.
    y_headroom = 1.50 if selection_label else 1.22
    ymax = max(1.0, float(np.nanmax(raw[mask]) * y_headroom))
    ax.set_ylim(0.0, ymax)
    style_axis(ax, ylabel=True, xlabel=xlabel)
    ax.text(
        0.035,
        0.955,
        r"$\bf{\it{sPHENIX}}$ Internal",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=13.2,
    )
    if canvas_label:
        ax.text(
            0.035,
            0.885,
            canvas_label,
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=11.0,
            color="#111827",
        )
    if selection_label:
        ax.text(
            0.035,
            0.812,
            selection_label,
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=8.9,
            color="#111827",
            linespacing=1.08,
            bbox=dict(boxstyle="round,pad=0.20", facecolor="white", edgecolor="none", alpha=0.72),
        )


def main() -> None:
    ensure_atlas_010_crop()
    for path in (ATLAS_PP, ATLAS_010):
        if not path.exists():
            raise FileNotFoundError(path)
    write_tight_crop(ATLAS_PP, ATLAS_PP_TIGHT)
    write_tight_crop(ATLAS_010, ATLAS_010_TIGHT)

    plt.rcParams.update(
        {
            "font.family": "Times New Roman",
            "font.serif": ["Times New Roman", "Times", "Nimbus Roman", "DejaVu Serif"],
            "mathtext.fontset": "stix",
            "axes.linewidth": 1.2,
        }
    )
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    pp_result = ppcomp.project_pp_with_ppg12_purity_norm()
    auau_panel = next(panel for panel in fig1.PANELS if panel.key == "auau_0_20")
    auau_result = fig1.build_components(auau_panel)

    fig = plt.figure(figsize=(16, 9), dpi=240)
    fig.patch.set_facecolor("white")
    title_color = "#111827"
    body_color = "#334155"

    fig.text(
        0.055,
        0.955,
        r"Reconstructed $x_{J\gamma}$ subtraction inputs before unfolding",
        ha="left",
        va="top",
        fontsize=29.0,
        fontweight="bold",
        color=title_color,
    )
    # Four matched square plot slots.  The width is set from the height so the
    # rendered panels are square on a 16:9 canvas, not just equal in axes units.
    panel_h = 0.385
    panel_w = panel_h * 9.0 / 16.0
    left_x = 0.055
    right_x = 0.335
    bot_y = 0.070
    top_y = 0.490
    header_y = top_y + panel_h + 0.020
    atlas_x = left_x
    atlas_top_y = top_y
    atlas_bot_y = bot_y
    fig.text(
        left_x + 0.5 * panel_w,
        header_y + 0.006,
        "ATLAS published reference",
        ha="center",
        va="center",
        fontsize=18.7,
        fontweight="bold",
        color=body_color,
    )
    fig.text(
        left_x + 0.5 * panel_w,
        header_y - 0.018,
        "Phys. Lett. B 789 (2019), Fig. 1",
        ha="center",
        va="center",
        fontsize=10.3,
        color=body_color,
    )
    fig.text(right_x + 0.5 * panel_w, header_y, "This analysis", ha="center", va="center", fontsize=20.5, fontweight="bold", color=body_color)

    # Left column: exact ATLAS figure crops.
    ax_atlas_pp = fig.add_axes([atlas_x, atlas_top_y, panel_w, panel_h])
    ax_atlas_010 = fig.add_axes([atlas_x, atlas_bot_y, panel_w, panel_h])
    for ax, path in ((ax_atlas_pp, ATLAS_PP_TIGHT), (ax_atlas_010, ATLAS_010_TIGHT)):
        ax.imshow(mpimg.imread(path), interpolation="lanczos", aspect="auto")
        ax.axis("off")

    # Right column: regenerated current-analysis panels.
    ax_pp = fig.add_axes([right_x, top_y, panel_w, panel_h])
    ax_auau = fig.add_axes([right_x, bot_y, panel_w, panel_h])
    draw_component_panel(
        ax_pp,
        pp_result,
        pp_panel=True,
        xlabel=False,
        canvas_label=r"p+p, $\sqrt{s}=200$ GeV",
        selection_label=(
            r"$15<E_T^\gamma<35$ GeV, $|\eta^\gamma|<0.7$"
            "\n"
            r"$|\Delta\phi|>7\pi/8$, $p_T^{jet}>5$ GeV, $|z_{vtx}|<60$ cm"
        ),
    )
    draw_component_panel(
        ax_auau,
        auau_result,
        pp_panel=False,
        canvas_label=r"Au+Au 0-20%, $\sqrt{s_{NN}}=200$ GeV",
        selection_label=(
            r"$15<E_T^\gamma<35$ GeV, $|\eta^\gamma|<0.7$"
            "\n"
            r"$|\Delta\phi|>7\pi/8$, $p_T^{jet}>5$ GeV, $|z_{vtx}|<10$ cm"
        ),
    )
    handles, labels = ax_pp.get_legend_handles_labels()
    leg = fig.legend(
        handles,
        labels,
        loc="upper right",
        bbox_to_anchor=(0.982, 0.925),
        ncol=2,
        frameon=True,
        fancybox=True,
        fontsize=10.8,
        handlelength=1.85,
        columnspacing=0.72,
        borderpad=0.32,
        labelspacing=0.32,
        handletextpad=0.45,
    )
    leg.get_frame().set_facecolor("#f8fafc")
    leg.get_frame().set_edgecolor("#cbd5e1")
    leg.get_frame().set_linewidth(1.0)
    leg.get_frame().set_alpha(0.96)
    note_x = 0.585
    note_w = 0.402

    def add_flow_box_rich(
        x,
        y,
        w,
        h,
        title,
        bullets,
        edge="#cbd5e1",
        fill="#fbfdff",
        accent="#334155",
        title_color_local=None,
        title_size=13.0,
        bullet_size=10.4,
        bullet_gap=0.022,
        bullet_start=0.043,
        title_offset=0.014,
        left_pad=0.012,
        bullet_indent=0.030,
    ):
        box = FancyBboxPatch(
            (x, y - h),
            w,
            h,
            boxstyle="round,pad=0.010,rounding_size=0.010",
            transform=fig.transFigure,
            facecolor=fill,
            edgecolor=edge,
            linewidth=1.15,
            zorder=0,
        )
        fig.patches.append(box)
        fig.text(
            x + left_pad,
            y - title_offset,
            title,
            ha="left",
            va="top",
            fontsize=title_size,
            fontweight="bold",
            color=title_color_local or title_color,
        )
        yy = y - bullet_start
        for bullet in bullets:
            fig.text(
                x + left_pad + 0.005,
                yy + 0.001,
                "\u2022",
                ha="left",
                va="top",
                fontsize=bullet_size + 0.6,
                fontweight="bold",
                color=accent,
            )
            fig.text(
                x + bullet_indent,
                yy,
                bullet,
                ha="left",
                va="top",
                fontsize=bullet_size,
                color=body_color,
                linespacing=1.08,
            )
            yy -= bullet_gap * (bullet.count("\n") + 1)

    def add_down_arrow(x, y0, y1):
        fig.patches.append(
            FancyArrowPatch(
                (x, y0),
                (x, y1),
                transform=fig.transFigure,
                arrowstyle="-|>",
                mutation_scale=13,
                linewidth=1.4,
                color="#64748b",
                zorder=1,
            )
        )

    rhs_title_y = 0.846
    fig.text(
        note_x,
        rhs_title_y,
        "Correction chain before unfolding",
        ha="left",
        va="center",
        fontsize=19.2,
        fontweight="bold",
        color=title_color,
    )

    eq_y = 0.815
    eq_h = 0.158
    eq_box = FancyBboxPatch(
        (note_x, eq_y - eq_h),
        note_w,
        eq_h,
        boxstyle="round,pad=0.010,rounding_size=0.010",
        transform=fig.transFigure,
        facecolor="#ffffff",
        edgecolor="#cbd5e1",
        linewidth=1.2,
        zorder=0,
    )
    fig.patches.append(eq_box)
    fig.text(
        note_x + 0.014,
        eq_y - 0.014,
        r"Per photon-$p_T$ row:",
        ha="left",
        va="top",
        fontsize=12.6,
        fontweight="bold",
        color=body_color,
    )
    def eq_text(label, *, color=title_color, fill=None):
        props = {
            "fontsize": 13.2,
            "color": color,
            "fontfamily": "DejaVu Serif",
        }
        if fill:
            props["bbox"] = {
                "boxstyle": "round,pad=0.14,rounding_size=0.05",
                "facecolor": fill,
                "edgecolor": "none",
            }
        return TextArea(label, textprops=props)

    equation = HPacker(
        children=[
            eq_text(r"$S_i^{reco}(x)=$"),
            eq_text(r"$[$"),
            eq_text(r"$A_i(x)$", color="#334155", fill="#f1f5f9"),
            eq_text(r"$-$"),
            eq_text(r"$\beta_i C_i(x)$", color="#1d4ed8", fill="#dbeafe"),
            eq_text(r"$]$"),
            eq_text(r"$/(1-$"),
            eq_text(r"$\beta_i f_C$", color="#047857", fill="#dcfce7"),
            eq_text(r"$)-$"),
            eq_text(r"$\alpha_i K_i^{MC}(x)$", color="#991b1b", fill="#fee2e2"),
        ],
        align="baseline",
        pad=0,
        sep=1.0,
    )
    fig.add_artist(
        AnnotationBbox(
            equation,
            (note_x + 0.012, eq_y - 0.054),
            xycoords=fig.transFigure,
            box_alignment=(0.0, 0.5),
            frameon=False,
            pad=0,
        )
    )

    def key_text(label, color, fontsize=10.0):
        return TextArea(
            label,
            textprops={
                "fontsize": fontsize,
                "color": color,
                "fontfamily": "DejaVu Serif",
                "fontweight": "bold",
            },
        )

    alpha_row = HPacker(
        children=[
            key_text(r"$\alpha_i=N_{\gamma,i}^{data,corr}/N_{\gamma,i}^{sigMC}$", "#991b1b", fontsize=11.1),
        ],
        align="baseline",
        pad=0,
        sep=0,
    )
    fig.add_artist(
        AnnotationBbox(
            alpha_row,
            (note_x + 0.014, eq_y - 0.090),
            xycoords=fig.transFigure,
            box_alignment=(0.0, 0.5),
            frameon=False,
            pad=0,
        )
    )

    key_row = HPacker(
        children=[
            key_text(r"$A_i$: raw", "#334155", fontsize=10.5),
            key_text(r"$\beta_i C_i$: wrong photon", "#1d4ed8", fontsize=10.5),
            key_text(r"$f_C$: C leakage", "#047857", fontsize=10.5),
            key_text(r"$\alpha_iK_i^{MC}$: wrong recoil", "#991b1b", fontsize=10.5),
        ],
        align="baseline",
        pad=0,
        sep=9,
    )
    fig.add_artist(
        AnnotationBbox(
            key_row,
            (note_x + 0.014, eq_y - 0.117),
            xycoords=fig.transFigure,
            box_alignment=(0.0, 0.5),
            frameon=False,
            pad=0,
        )
    )

    bold_row = HPacker(
        children=[
            key_text("blue = wrong photon; red = right photon, wrong recoil jet", "#111827", fontsize=11.3),
        ],
        align="baseline",
        pad=0,
        sep=0,
    )
    fig.add_artist(
        AnnotationBbox(
            bold_row,
            (note_x + 0.014, eq_y - 0.143),
            xycoords=fig.transFigure,
            box_alignment=(0.0, 0.5),
            frameon=False,
            pad=0,
        )
    )

    box_x = note_x
    box_w = note_w
    add_flow_box_rich(
        box_x,
        0.646,
        box_w,
        0.086,
        r"1. Raw $A_i(x)$: gray",
        [
            r"tight + isolated photon candidates",
            r"paired with selected recoil jets",
        ],
        edge="#cbd5e1",
        fill="#f8fafc",
        accent="#475569",
        title_size=12.7,
        bullet_size=11.4,
        bullet_gap=0.023,
        bullet_start=0.044,
    )
    add_down_arrow(box_x + 0.5 * box_w, 0.556, 0.543)
    add_flow_box_rich(
        box_x,
        0.535,
        box_w,
        0.127,
        r"2. Subtract $\beta_i C_i(x)$: blue wrong-photon template",
        [
            r"$C_i(x)$: isolated non-tight $x_{J\gamma}$ shape",
            r"$\beta_i$ normalizes $C_i$ to fake yield from purity fit",
            r"physics: high-$z$ $\pi^0/\eta$ fragments from dijet events",
        ],
        edge="#bfdbfe",
        fill="#f8fbff",
        accent="#2563eb",
        title_color_local="#1d4ed8",
        title_size=12.0,
        bullet_size=11.0,
        bullet_gap=0.026,
        bullet_start=0.046,
    )
    add_down_arrow(box_x + 0.5 * box_w, 0.403, 0.390)
    add_flow_box_rich(
        box_x,
        0.382,
        box_w,
        0.170,
        r"3. Subtract $\alpha_i K_i^{MC}(x)$: red wrong-recoil template",
        [
            r"$K_i^{MC}(x)$: true prompt $\gamma$ + selected reco recoil jets with no truth-jet match",
            r"$\alpha_i=N_{\gamma,i}^{data,corr}/N_{\gamma,i}^{sigMC}$",
            r"meaning: random recoil jets per prompt $\gamma$ scaled to data prompt-$\gamma$ yield",
            r"physics: real $\gamma$ paired with unrelated HI recoil jet; pp term $\approx 0$",
        ],
        edge="#fecaca",
        fill="#fffafa",
        title_color_local="#991b1b",
        accent="#dc2626",
        title_size=12.1,
        bullet_size=10.6,
        bullet_gap=0.027,
        bullet_start=0.047,
    )
    add_down_arrow(box_x + 0.5 * box_w, 0.207, 0.194)
    add_flow_box_rich(
        box_x,
        0.186,
        box_w,
        0.090,
        r"4. Output: black points",
        [
            r"corrected reconstructed $x_{J\gamma}$ spectrum",
            r"input to response-matrix unfolding",
        ],
        edge="#ddd6fe",
        fill="#fdfcff",
        accent="#111827",
        title_size=12.8,
        bullet_size=11.4,
        bullet_gap=0.024,
        bullet_start=0.044,
    )

    # Keep the slide audience-facing; detailed provenance is in the manifest.
    fig.savefig(OUT_PNG, dpi=240)
    plt.close(fig)

    manifest = {
        "slide": str(OUT_PNG),
        "atlas": {
            "source_pdf": str(REPO / "usefulDocs/external_references/gamma_jet/ATLAS_xJgamma.pdf"),
            "pp_crop": str(ATLAS_PP),
            "pbpb_0_10_crop": str(ATLAS_010),
            "pp_tight_crop": str(ATLAS_PP_TIGHT),
            "pbpb_0_10_tight_crop": str(ATLAS_010_TIGHT),
            "note": "Exact visual crops from ATLAS Figure 1; not redrawn.",
        },
        "this_analysis": {
            "pp": {
                "data_file": str(fig1.PANELS[0].data_file),
                "sideband_normalization": "Current THE76 region-C xJ shape normalized by a Padé[1/1] fit to the PPG12 final-BDT leakage-corrected purity graph.",
                "sideband_label": "Dijet Bkg.; current pp visual uses region-C xJ shape normalized by the fitted PPG12 final-BDT leakage-corrected purity curve.",
                "integrals": pp_result["integrals"],
                "purity_fit": pp_result.get("purity_fit"),
                "pt_rows": pp_result["pt_rows"],
                "base_key": fig1.PANELS[0].base_key,
                "effective_base_key": fig1.effective_base_key(fig1.PANELS[0]),
            },
            "auau_0_20": {
                "data_file": str(auau_panel.data_file),
                "sim_file": str(auau_panel.sim_file),
                "base_key": auau_panel.base_key,
                "effective_base_key": fig1.effective_base_key(auau_panel),
                "photon_key": auau_panel.pho_key,
                "integrals": auau_result["integrals"],
                "comb_meta": auau_result["comb_meta"],
                "purity_fit": auau_result.get("purity_fit"),
                "pt_rows": auau_result["pt_rows"],
                "sideband_label": "Dijet Bkg.; net leakage-aware region-C subtraction using fitted leakage-corrected AuAu 0-20 purity. Direct scale_C * H_C is stored in the detailed slide07 manifest.",
            },
        },
        "selection": {
            "photon": fig1.photon_pt_label(),
            "photon_eta": "|eta_gamma| < 0.7",
            "back_to_back": "nominal |Delta phi| > 7pi/8 object family",
            "jet": fig1.jet_pt_label(),
            "pp_vertex": "|z_vtx| < 60 cm from source row name vz60",
            "auau_vertex": "|z_vtx| < 10 cm offline; source topdir name is trigger family vtx_lt_150",
        },
        "require_full_pt_bins": fig1.u.REQUIRE_FULL_PT_BINS,
        "jet_pt_key": fig1.JET_PT_KEY or "nominal_key",
        "status": "Reconstructed-input diagnostic before unfolding; not a final unfolded particle-level comparison.",
        "caveats": [
            "Generated from the existing first-pass 15-35 GeV response-input products; the visible slide intentionally omits old diagnostic caveat text.",
            "AuAu red combinatoric template is taken from embedded signal MC and normalized row-by-row to the data photon yield.",
            "Blue is the fitted-purity scaled region-C photon-ID sideband contribution in this analysis comparison.",
        ],
    }
    OUT_MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n")
    OUT_SCRIPT.write_text(
        "This slide lines up the ATLAS Figure 1 input-subtraction picture against our current first-pass inputs. "
        "The left column is the published ATLAS visual evidence: pp on top and central Pb+Pb on bottom. "
        "The right column is our corresponding reconstructed xJ input before unfolding: pp on top and Au+Au 0-20 percent on bottom.\n\n"
        "Read the colors by role. Grey is raw region A, red is the combinatoric background where applicable, blue is the region-C dijet or fake-photon background subtraction scaled by the fitted purity, and black points are the background-subtracted input. "
        "The pp scale uses a Padé fit to the PPG12 final-BDT purity graph; the Au+Au scale uses the Padé fit to the leakage-corrected Au+Au photon-candidate purity table. "
        "The main check is whether the subtraction ingredients are ordered sensibly and whether the central heavy-ion input shows the expected larger low-xJ background structure before unfolding.\n"
    )

    print(json.dumps({"slide": str(OUT_PNG), "manifest": str(OUT_MANIFEST), "speaker_script": str(OUT_SCRIPT)}, indent=2))


if __name__ == "__main__":
    main()
