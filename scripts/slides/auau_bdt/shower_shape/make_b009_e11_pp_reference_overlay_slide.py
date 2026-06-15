#!/usr/bin/env python3
"""Build an E11/E33 AuAu data/MC stage-flow slide with a pp reference row."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


THIS_FILE = Path(__file__).resolve()
REPO = next((p for p in THIS_FILE.parents if (p / "AGENTS.md").exists()), THIS_FILE.parents[5])
SCRIPTS_DIR = REPO / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.append(str(SCRIPTS_DIR))
if str(THIS_FILE.parent) not in sys.path:
    sys.path.append(str(THIS_FILE.parent))

from slides.common.slide_defaults import SLIDE_DPI, slide_figsize  # noqa: E402

import make_the42_b009_pre_vs_tight_shower_shape_slide as b009  # noqa: E402


PP_CAMPAIGN = REPO / "dataOutput/ppg12TableQA/THE42_ppg12_tableqa_v1_basev3e_20260611"
PP_ROOT_DIR = PP_CAMPAIGN / "merged_roots"
PP_INCLUSIVE_CACHE = PP_CAMPAIGN / "inclusive_sample_hist_cache/the42_ppg12_tableqa_v1_inclusive_sample_projectx_hists.json"
PP_PLOTTER_PATH = REPO / "scripts/plotting/pp_currentian/make_the42_ppg12_tableqa_v1_tables.py"
PP_DATA_ROOT = PP_ROOT_DIR / "RecoilJets_pp_ALL_preselectionNewPPG12_tightNewPPG12_nonTightNewPPG12.root"
PP_SIGNAL_ROOT = PP_ROOT_DIR / "RecoilJets_photonjet5plus10plus20_MERGED.root"
PP_STITCHED_INCLUSIVE_ROOT = PP_ROOT_DIR / "RecoilJets_jet5plus8plus12plus20plus30plus40_MERGED.root"
PP_TOPDIR_DATA = "Photon_4_GeV_plus_MBD_NS_geq_1"

OUTDIR = b009.OUTDIR / "e11_pp_reference_overlay"
VAR = "e11e33"
PP_VAR = "e11_to_e33"
PT_TOKEN = "1535"
PT_LABEL = "15 < E_T < 35 GeV"

STAGES = [
    ("Before preselection", "inclusive", "cut0"),
    ("After preselection", "pre", "cut1"),
    ("After tight ID", "tight", "cut2"),
]

AU_AU_SAMPLE_TAGS = {
    "inclusive": {
        "Data": "inclusive",
        "Signal MC": "inclusive_sig",
        "Inclusive MC": "inclusive_bkg",
    },
    "pre": {
        "Data": "pre",
        "Signal MC": "pre_sig",
        "Inclusive MC": "pre_bkg",
    },
    "tight": {
        "Data": "tight",
        "Signal MC": "tight_sig",
        "Inclusive MC": "tight_bkg",
    },
}

SAMPLE_COLORS = {
    "Data": "#111827",
    "Signal MC": "#C22F2F",
    "Inclusive MC": "#2E63D4",
}


def load_pp_plotter():
    spec = importlib.util.spec_from_file_location("pp_tableqa_plotter", PP_PLOTTER_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not import pp plotter from {PP_PLOTTER_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def add_sphenix_internal(ax, *, collision: str, details: str) -> None:
    ax.text(
        0.025,
        0.955,
        r"$\it{\bf{sPHENIX}}$ Internal",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8.8,
    )
    ax.text(
        0.025,
        0.805,
        f"{collision}\n{details}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=7.4,
        linespacing=1.05,
    )


def draw_b009_curve(ax, curve: b009.Curve, color: str) -> None:
    if curve.sample == "Data":
        ax.errorbar(
            curve.centers,
            curve.values,
            yerr=curve.errors,
            fmt="o",
            ms=3.2,
            color=color,
            mfc=color,
            mec="white",
            mew=0.55,
            elinewidth=0.75,
            capsize=1.5,
            capthick=0.75,
            alpha=0.95,
        )
    else:
        ax.step(curve.edges[:-1], curve.values, where="post", color=color, lw=1.75, alpha=0.96)


def draw_pp_shape(ax, arrays, *, label: str, color: str, marker: str | None = None, linewidth: float = 1.7) -> None:
    x, y, e = arrays
    if marker:
        ax.errorbar(x, y, yerr=e, fmt=marker, ms=3.0, color=color, mfc=color, mec="white", mew=0.55, elinewidth=0.75, capsize=1.4, label=label)
    else:
        ax.step(x, y, where="mid", color=color, lw=linewidth, label=label)
        ax.errorbar(x, y, yerr=e, fmt="none", ecolor=color, elinewidth=0.45, capsize=0, alpha=0.65)


def pp_arrays(plotter, files: dict, inclusive_cache: dict, stage: str) -> dict[str, tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray] | None]:
    xlim, rebin = plotter.ppg12_axis_settings(PP_VAR)
    out: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray] | None] = {}
    data = plotter.norm_arrays(plotter.get_hist(files["data"], PP_TOPDIR_DATA, PP_VAR, PT_TOKEN, stage), rebin, xlim)
    sig = plotter.norm_arrays(plotter.get_hist(files["signal_mc"], "SIM", PP_VAR, PT_TOKEN, stage), rebin, xlim)
    inc_payload = plotter.cache_hist(inclusive_cache, "current_ian_jet8to40", PP_VAR, PT_TOKEN, stage)
    inc = plotter.norm_payload(inc_payload, rebin, xlim)
    out["Data"] = data
    out["Signal MC"] = sig
    out["Inclusive MC"] = inc
    if stage == "cut0":
        npb_scale = plotter.npb_tail_scale(
            files,
            {
                "pt_token": PT_TOKEN,
                "pt_label": r"$15<E_T<35$ GeV",
                "cut": stage,
                "cut_label": "all candidates",
                "include_npb_template": True,
            },
        )
        npb = plotter.norm_arrays(plotter.get_hist(files["data"], PP_TOPDIR_DATA, PP_VAR, PT_TOKEN, "cut4"), rebin, xlim)
        out["NPB-tagged data"] = plotter.scale_arrays(npb, npb_scale)
    return out


def style_axis(ax, *, row: int, col: int, row_label: str, stage_label: str) -> None:
    ax.set_xlim(0.0, 1.02)
    ax.grid(True, axis="y", color="#E5E7EB", lw=0.55, alpha=0.78)
    ax.tick_params(labelsize=8.2, pad=1, direction="in", top=True, right=True)
    if row == 3:
        ax.set_xlabel(r"$E_{11}/E_{33}$", fontsize=10.8, labelpad=1)
    else:
        ax.set_xlabel("")
        ax.tick_params(labelbottom=False)
    if col == 0:
        ax.set_ylabel(row_label, fontsize=10.4, labelpad=8)
    else:
        ax.set_ylabel("")
    if row == 0:
        ax.set_title(stage_label, fontsize=15.5, fontweight="bold", pad=6, color="#173b63")


def main() -> int:
    b009.setup_style()
    OUTDIR.mkdir(parents=True, exist_ok=True)
    b009_indexes = {root: b009.build_index(root) for _, root, _, _ in b009.SAMPLES}
    pp_plotter = load_pp_plotter()
    pp_files = {
        "data": pp_plotter.open_root(PP_DATA_ROOT),
        "signal_mc": pp_plotter.open_root(PP_SIGNAL_ROOT),
    }
    pp_inclusive_cache = pp_plotter.load_inclusive_cache(PP_INCLUSIVE_CACHE, use_stitched_inclusive=False)

    fig, axes = plt.subplots(4, 3, figsize=slide_figsize(), constrained_layout=False)
    fig.subplots_adjust(left=0.087, right=0.985, top=0.765, bottom=0.087, wspace=0.135, hspace=0.260)
    fig.patch.set_facecolor("white")

    manifest: dict = {
        "schema": "B009_E11_PP_REFERENCE_OVERLAY_SLIDE_V1",
        "variable": "E11/E33",
        "pt_label": PT_LABEL,
        "auau_inputs": {
            "data": str(b009.DATA_ROOT),
            "signal_mc": str(b009.SIGNAL_ROOT),
            "inclusive_mc": str(b009.INCLUSIVE_ROOT),
        },
        "pp_inputs": {
            "data": str(PP_DATA_ROOT),
            "signal_mc": str(PP_SIGNAL_ROOT),
            "inclusive_sample_cache": str(PP_INCLUSIVE_CACHE),
            "stitched_inclusive_not_used_for_plot": str(PP_STITCHED_INCLUSIVE_ROOT),
        },
        "curves": [],
    }

    for row, (cent_label, cent_bins) in enumerate(b009.CENT_GROUPS):
        for col, (stage_label, auau_stage, _) in enumerate(STAGES):
            ax = axes[row, col]
            ymax = 0.0
            for sample, root, _, color in b009.SAMPLES:
                tag = AU_AU_SAMPLE_TAGS[auau_stage][sample]
                curve = b009.load_curve(
                    root,
                    b009_indexes[root],
                    sample,
                    stage_label,
                    tag,
                    VAR,
                    cent_bins,
                    (0.0, 1.02),
                    2,
                )
                draw_b009_curve(ax, curve, color)
                if len(curve.values):
                    ymax = max(ymax, float(np.max(curve.values + curve.errors)))
                manifest["curves"].append(
                    {
                        "system": "AuAu",
                        "centrality": cent_label,
                        "stage": stage_label,
                        "sample": sample,
                        "hist_tag": tag,
                        "integral": curve.integral,
                        "matches": len(curve.matches),
                        "missing": len(curve.missing),
                    }
                )
            ax.set_ylim(0, max(0.04, ymax * 1.18))
            style_axis(ax, row=row, col=col, row_label=f"AuAu\n{cent_label}", stage_label=stage_label)
            if row == 0 and col == 0:
                add_sphenix_internal(ax, collision=r"Au+Au $\sqrt{s_{NN}}=200$ GeV", details=PT_LABEL)

    pp_row = 3
    for col, (stage_label, _, pp_stage) in enumerate(STAGES):
        ax = axes[pp_row, col]
        arrays = pp_arrays(pp_plotter, pp_files, pp_inclusive_cache, pp_stage)
        ymax = 0.0
        for sample in ("Data", "Signal MC", "Inclusive MC"):
            arr = arrays[sample]
            if arr is None:
                continue
            x, y, e = arr
            ymax = max(ymax, float(np.max(y + e)) if len(y) else 0.0)
            draw_pp_shape(
                ax,
                arr,
                label=sample,
                color=SAMPLE_COLORS[sample],
                marker="o" if sample == "Data" else None,
            )
        if arrays.get("NPB-tagged data") is not None:
            arr = arrays["NPB-tagged data"]
            x, y, e = arr  # type: ignore[misc]
            ymax = max(ymax, float(np.max(y + e)) if len(y) else 0.0)
            draw_pp_shape(ax, arr, label="NPB-tagged data", color="#238b1e", linewidth=1.45)  # type: ignore[arg-type]
        ax.set_ylim(0, max(0.04, ymax * 1.18))
        style_axis(ax, row=pp_row, col=col, row_label="pp\nreference", stage_label=stage_label)
        if col == 0:
            add_sphenix_internal(ax, collision=r"$p$+$p$ $\sqrt{s}=200$ GeV", details=PT_LABEL)

    title = "E11/E33 overlay: AuAu data compared to signal and inclusive MC"
    fig.text(0.055, 0.955, title, fontsize=24.5, fontweight="bold", ha="left", va="top", color="#111827")
    fig.text(
        0.055,
        0.897,
        "   Each AuAu panel overlays data markers with embedded signal and inclusive MC lines across the photon-ID selection flow.",
        fontsize=14.6,
        ha="left",
        va="top",
        color="#172033",
    )
    fig.text(
        0.055,
        0.858,
        "   Bottom row is the repaired pp table-QA reference; active b002 data will be regenerated with matched table-QA MC when available.",
        fontsize=14.6,
        ha="left",
        va="top",
        color="#172033",
    )
    # Draw blue arrowheads with a font that contains the glyph; keep the body text in Times.
    for y in (0.897, 0.858):
        fig.text(0.055, y, "▶", fontsize=15.2, color="#2468a8", ha="left", va="top", fontfamily="DejaVu Sans")

    legend_items = [
        plt.Line2D([0], [0], color="#C22F2F", lw=2.2, label="Signal MC"),
        plt.Line2D([0], [0], color="#2E63D4", lw=2.2, label="Inclusive MC"),
        plt.Line2D([0], [0], color="#111827", marker="o", markersize=6.5, lw=0, markerfacecolor="#111827", markeredgecolor="white", label="Data"),
        plt.Line2D([0], [0], color="#238b1e", lw=1.8, label="NPB-tagged data (pp cut0 only)"),
    ]
    fig.legend(handles=legend_items, loc="upper right", bbox_to_anchor=(0.985, 0.825), frameon=False, ncol=4, fontsize=9.4, handlelength=1.6, columnspacing=1.2)

    out_png = OUTDIR / "b009_e11_to_e33_auau_data_mc_with_pp_reference_slide.png"
    out_manifest = OUTDIR / "b009_e11_to_e33_auau_data_mc_with_pp_reference_slide_manifest.json"
    out_script = OUTDIR / "b009_e11_to_e33_auau_data_mc_with_pp_reference_slide_speaker_script.md"
    fig.savefig(out_png, dpi=SLIDE_DPI)
    plt.close(fig)

    manifest["png"] = str(out_png)
    manifest["speaker_script"] = str(out_script)
    manifest["caveat"] = (
        "AuAu rows use completed b009 shower-shape products with existing embedded MC. "
        "This is the correct MC-overlaid shower-shape comparison, not the active b002 table-QA data-only subset."
    )
    out_manifest.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    out_script.write_text(
        "\n".join(
            [
                "# Speaker Script",
                "",
                "This slide focuses only on E11 over E33, because it is the cleanest shower-shape handle for seeing the photon-ID selection flow.",
                "Each AuAu row is a centrality bin. In every panel the black markers are AuAu data, the red curve is embedded signal MC, and the blue curve is embedded inclusive MC.",
                "The columns show the same progression as the pp reference: before preselection, after preselection, and after tight ID.",
                "The bottom row is the repaired pp table-QA reference, included to remind the audience what the validated photon-ID behavior looks like in pp.",
                "",
                "The visual readout is the rightmost column. After tight ID, the AuAu data move into the same high-E11/E33 region where signal MC is concentrated, while inclusive MC remains the broader comparison shape.",
                "This version uses the completed b009 shower-shape products for the AuAu data and embedded MC. The active b002 table-QA data will need matched table-QA MC before it can replace this overlay one-for-one.",
                "",
            ]
        )
    )
    print(out_png)
    print(out_manifest)
    print(out_script)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
