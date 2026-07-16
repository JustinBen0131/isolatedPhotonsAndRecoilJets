#!/usr/bin/env python3
"""Build the partial-data AuAu isolation-template slide for three centralities."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import subprocess
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


THIS_FILE = Path(__file__).resolve()
REPO = next((p for p in THIS_FILE.parents if (p / "AGENTS.md").exists()), THIS_FILE.parents[4])
SCRIPTS = REPO / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.append(str(SCRIPTS))
sys.path.append(str(REPO / "scripts/plotting/pp_currentian"))

from make_ppg12_style_isolation_stack import (  # noqa: E402
    add_hists,
    bin_width_scale,
    draw_step_band,
    ppg12_variable_rebin,
    tail_integral,
)
from slides.common.slide_defaults import (  # noqa: E402
    SLIDE_DPI,
    SLIDE_HEIGHT_PX,
    SLIDE_WIDTH_PX,
    slide_figsize,
)


DATA_ROOT = (
    REPO
    / "InputFiles/the88_auau_data_centrality_fixed_20260715T1113_partial_diagnostic"
    / "RecoilJets_auau_ALL_preselectionNewPPG12_tightAuAuCentInputBase3x3BDT_"
    "nonTightAuAuBDTSideband_baseVariant.root"
)
SIGNAL_ROOT = (
    REPO
    / "InputFiles/the88_bounded_nontight_20260708/simembedded"
    / "preselectionNewPPG12_tightAuAuCentInputBase3x3BDT_nonTightAuAuBDTSideband_baseVariant"
    / "photonJet12and20merged_SIM/RecoilJets_embeddedPhoton12plus20_MERGED.root"
)
OUT_DIR = (
    REPO
    / "dataOutput/auau/the88_auau_data_centrality_fixed_20260713"
    / "isolation_template_partial"
)
STEM = "auau_partial_isolation_template_15to35_by_centrality"
OUT_PNG = OUT_DIR / f"{STEM}.png"
OUT_CSV = OUT_DIR / f"{STEM}.csv"
OUT_MANIFEST = OUT_DIR / f"{STEM}.manifest.json"
OUT_LAYOUT = OUT_DIR / f"{STEM}.layout_nodes.json"
OUT_AUDIT = OUT_DIR / f"{STEM}.audit.json"
OUT_SCRIPT = OUT_DIR / f"{STEM}_speaker_script.md"

DATA_DIR = "MBD_NS_geq_2_vtx_lt_150"
SIGNAL_DIR = "SIM"
PT_BINS = ("15_17", "17_19", "19_21", "21_23", "23_26", "26_35")
CENTRALITIES = (("0_20", "0-20%"), ("20_50", "20-50%"), ("50_80", "50-80%"))
TAIL_LOW_GEV = 6.0

INK = "#111827"
MUTED = "#475569"
EDGE = "#cbd5e1"
RED_FILL = "#f4aaa7"
RED_EDGE = "#8a403c"
BLUE_FILL = "#9699ee"
BLUE_EDGE = "#2f327d"
AMBER = "#a16207"
AMBER_FILL = "#fffbeb"

TIMES_DIR = Path("/System/Library/Fonts/Supplemental")
TIMES_FONTS = (
    TIMES_DIR / "Times New Roman.ttf",
    TIMES_DIR / "Times New Roman Bold.ttf",
    TIMES_DIR / "Times New Roman Italic.ttf",
    TIMES_DIR / "Times New Roman Bold Italic.ttf",
)


def configure_fonts() -> None:
    for path in TIMES_FONTS:
        if path.exists():
            font_manager.fontManager.addfont(str(path))
    plt.rcParams.update(
        {
            "font.family": "Times New Roman",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "custom",
            "mathtext.rm": "Times New Roman",
            "mathtext.it": "Times New Roman:italic",
            "mathtext.bf": "Times New Roman:bold",
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "axes.unicode_minus": False,
            "axes.linewidth": 1.6,
        }
    )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def build_panel(cent_key: str) -> dict[str, object]:
    tight_pattern = f"h_Eiso_tight_isoR40_pT_{{pt}}_cent_{cent_key}"
    nontight_pattern = f"h_Eiso_nonTight_isoR40_pT_{{pt}}_cent_{cent_key}"
    signal_pattern = f"h_EisoReco_truthSigMatched_tight_isoR40_pT_{{pt}}_cent_{cent_key}"

    tight, tight_var, edges, tight_keys = add_hists(DATA_ROOT, DATA_DIR, tight_pattern, PT_BINS)
    nontight, nontight_var, nontight_edges, nontight_keys = add_hists(
        DATA_ROOT, DATA_DIR, nontight_pattern, PT_BINS
    )
    signal, signal_var, signal_edges, signal_keys = add_hists(
        SIGNAL_ROOT, SIGNAL_DIR, signal_pattern, PT_BINS
    )
    if not (np.array_equal(edges, nontight_edges) and np.array_equal(edges, signal_edges)):
        raise ValueError(f"Isolation binning mismatch in centrality {cent_key}")

    tight, tight_var, rebinned_edges = ppg12_variable_rebin(
        tight, tight_var, edges, low_region_group=1, high_region_group=5
    )
    nontight, nontight_var, _ = ppg12_variable_rebin(
        nontight, nontight_var, edges, low_region_group=1, high_region_group=5
    )
    signal, signal_var, _ = ppg12_variable_rebin(
        signal, signal_var, edges, low_region_group=1, high_region_group=5
    )

    tight, tight_var = bin_width_scale(tight, tight_var, rebinned_edges)
    nontight, nontight_var = bin_width_scale(nontight, nontight_var, rebinned_edges)
    signal, signal_var = bin_width_scale(signal, signal_var, rebinned_edges)

    tight_tail = tail_integral(tight, rebinned_edges, TAIL_LOW_GEV)
    nontight_tail = tail_integral(nontight, rebinned_edges, TAIL_LOW_GEV)
    if not (nontight_tail > 0.0 and tight_tail >= 0.0):
        raise ValueError(f"Invalid tail integrals for centrality {cent_key}")
    background_scale = tight_tail / nontight_tail
    background = nontight * background_scale
    background_var = nontight_var * background_scale * background_scale

    signed_remainder = float(np.sum(tight) - np.sum(background))
    signal_integral = float(np.sum(signal))
    physical_signal_normalization = signed_remainder > 0.0 and signal_integral > 0.0
    signal_scale = signed_remainder / signal_integral if physical_signal_normalization else 0.0
    signal_scaled = signal * signal_scale
    signal_scaled_var = signal_var * signal_scale * signal_scale

    return {
        "cent_key": cent_key,
        "edges": rebinned_edges,
        "tight": tight,
        "tight_var": tight_var,
        "background": background,
        "background_var": background_var,
        "signal": signal_scaled,
        "signal_var": signal_scaled_var,
        "tight_keys": tight_keys,
        "nontight_keys": nontight_keys,
        "signal_keys": signal_keys,
        "tight_tail": tight_tail,
        "nontight_tail": nontight_tail,
        "background_scale": background_scale,
        "signed_remainder": signed_remainder,
        "signal_integral_before_scale": signal_integral,
        "signal_scale": signal_scale,
        "physical_signal_normalization": physical_signal_normalization,
    }


def render(results: list[dict[str, object]]) -> None:
    configure_fonts()
    fig = plt.figure(figsize=slide_figsize(), dpi=SLIDE_DPI, facecolor="white")
    fig.text(
        0.050,
        0.952,
        "Isolation-energy templates by centrality",
        fontsize=34,
        fontweight="bold",
        color=INK,
        ha="left",
        va="top",
    )
    fig.text(
        0.950,
        0.943,
        "PRELIMINARY: PARTIAL DATA",
        fontsize=14.5,
        fontweight="bold",
        color=AMBER,
        ha="right",
        va="top",
    )

    fig.text(
        0.500,
        0.855,
        "Au+Au data and embedded photon+jet 12 + 20 simulation",
        fontsize=19.0,
        color=INK,
        ha="center",
        va="center",
    )
    fig.text(
        0.500,
        0.815,
        r"$15\leq E_T^\gamma<35$ GeV, $|\eta^\gamma|<0.7$, reconstructed isolation $R=0.4$",
        fontsize=17.5,
        color=MUTED,
        ha="center",
        va="center",
    )
    fig.text(
        0.500,
        0.775,
        r"Bounded non-tight data tail-matched to tight data above 6 GeV; truth-isolated, truth-matched tight signal MC",
        fontsize=15.7,
        color=MUTED,
        ha="center",
        va="center",
    )

    legend_handles = [
        Line2D([0], [0], marker="o", color="black", lw=0, markersize=8, label="tight data"),
        Patch(facecolor=RED_FILL, edgecolor=RED_EDGE, alpha=0.78, label="tail-matched bounded non-tight data"),
        Patch(facecolor=BLUE_FILL, edgecolor=BLUE_EDGE, alpha=0.66, label="truth-matched signal MC"),
    ]
    fig.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.50, 0.718),
        ncol=3,
        frameon=False,
        fontsize=15.0,
        handlelength=1.25,
        columnspacing=2.0,
    )

    grid = fig.add_gridspec(1, 3, left=0.073, right=0.965, bottom=0.195, top=0.595, wspace=0.205)
    axes = [fig.add_subplot(grid[0, i]) for i in range(3)]
    for ax, result, (_, cent_label) in zip(axes, results, CENTRALITIES):
        edges = np.asarray(result["edges"], dtype=float)
        tight = np.asarray(result["tight"], dtype=float)
        tight_err = np.sqrt(np.clip(np.asarray(result["tight_var"], dtype=float), 0.0, None))
        background = np.asarray(result["background"], dtype=float)
        signal = np.asarray(result["signal"], dtype=float)
        centers = edges[:-1] + 0.5 * np.diff(edges)

        draw_step_band(
            ax,
            edges,
            np.zeros_like(background),
            background,
            RED_FILL,
            RED_EDGE,
            "tail-matched non-tight data",
            alpha=0.78,
            linewidth=1.35,
        )
        if bool(result["physical_signal_normalization"]):
            draw_step_band(
                ax,
                edges,
                background,
                background + signal,
                BLUE_FILL,
                BLUE_EDGE,
                "truth-matched signal MC",
                alpha=0.66,
                linewidth=1.35,
            )
        ax.errorbar(
            centers,
            tight,
            yerr=tight_err,
            fmt="o",
            color="black",
            ecolor="black",
            elinewidth=1.35,
            markersize=5.2,
            capsize=0,
            zorder=8,
        )

        ymax = max(float(np.max(tight + tight_err)), float(np.max(background + signal)))
        ax.set_xlim(-2.0, 15.0)
        ax.set_xticks([-2, 0, 5, 10, 15])
        ax.set_ylim(0.0, ymax * 1.25 if ymax > 0.0 else 1.0)
        ax.set_title(f"Au+Au {cent_label}", fontsize=22, fontweight="bold", pad=10)
        ax.set_xlabel(r"$E_T^{\mathrm{iso,reco}}$ [GeV]", fontsize=17)
        ax.tick_params(axis="both", which="major", direction="in", top=True, right=True, labelsize=14, length=7)
        ax.tick_params(axis="both", which="minor", direction="in", top=True, right=True, length=4)
        ax.minorticks_on()
        ax.grid(axis="y", color="#e5e7eb", linewidth=0.7, alpha=0.70)

        ax.text(
            0.95,
            0.94,
            r"$\it{\bf{sPHENIX}}$ Internal",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=15.5,
        )
        ax.text(
            0.95,
            0.865,
            r"$\sqrt{s_{NN}}=200$ GeV",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=14.0,
        )
    axes[0].set_ylabel("Counts / bin width", fontsize=17)

    negative_labels = [
        label
        for result, (_, label) in zip(results, CENTRALITIES)
        if not bool(result["physical_signal_normalization"])
    ]
    if negative_labels:
        closure = "No positive signal remainder in " + " and ".join(negative_labels) + "."
    else:
        closure = "All centrality panels have a positive signal remainder."
    fig.text(
        0.500,
        0.045,
        closure + "  Regenerate after the corrected Au+Au data campaign reaches full coverage.",
        fontsize=14.8,
        color=AMBER,
        ha="center",
        va="center",
        fontweight="bold",
    )

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=SLIDE_DPI, facecolor="white")
    plt.close(fig)


def write_csv(results: list[dict[str, object]]) -> None:
    fieldnames = [
        "centrality",
        "eiso_low_gev",
        "eiso_high_gev",
        "tight_data_density",
        "tight_data_error",
        "scaled_nontight_density",
        "scaled_signal_mc_density",
        "background_tail_scale",
        "signed_tight_minus_background_integral",
        "physical_signal_normalization",
    ]
    with OUT_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for result, (_, label) in zip(results, CENTRALITIES):
            edges = np.asarray(result["edges"], dtype=float)
            tight = np.asarray(result["tight"], dtype=float)
            tight_err = np.sqrt(np.clip(np.asarray(result["tight_var"], dtype=float), 0.0, None))
            background = np.asarray(result["background"], dtype=float)
            signal = np.asarray(result["signal"], dtype=float)
            for i in range(len(tight)):
                writer.writerow(
                    {
                        "centrality": label,
                        "eiso_low_gev": edges[i],
                        "eiso_high_gev": edges[i + 1],
                        "tight_data_density": tight[i],
                        "tight_data_error": tight_err[i],
                        "scaled_nontight_density": background[i],
                        "scaled_signal_mc_density": signal[i],
                        "background_tail_scale": result["background_scale"],
                        "signed_tight_minus_background_integral": result["signed_remainder"],
                        "physical_signal_normalization": result["physical_signal_normalization"],
                    }
                )


def write_layout() -> None:
    payload = {
        "schema": "slide_layout_nodes_v1",
        "slide_size": [SLIDE_WIDTH_PX, SLIDE_HEIGHT_PX],
        "title_axis_x": 128,
        "minimum_title_font_px": 68,
        "minimum_audience_font_px": 27,
        "nodes": [
            {
                "kind": "text",
                "name": "title",
                "role": "title",
                "title_anchor": True,
                "text": "Isolation-energy templates by centrality",
                "bbox": [128, 65, 1850, 145],
                "font_px": 68,
            },
            {"kind": "panel", "name": "partial status", "bbox": [1960, 75, 2432, 130]},
            {
                "kind": "text",
                "name": "sample definition",
                "role": "audience",
                "text": "AuAu data and embedded photon+jet 12 + 20 simulation; 15 <= ETgamma < 35 GeV, R=0.4",
                "bbox": [250, 190, 2310, 240],
                "font_px": 38,
            },
            {
                "kind": "text",
                "name": "normalization definition",
                "role": "audience",
                "text": "Bounded non-tight tail match above 6 GeV; truth-isolated truth-matched tight signal MC",
                "bbox": [250, 245, 2310, 330],
                "font_px": 31,
            },
            {"kind": "panel", "name": "0-20 panel", "bbox": [187, 583, 846, 1159]},
            {"kind": "panel", "name": "20-50 panel", "bbox": [945, 583, 1604, 1159]},
            {"kind": "panel", "name": "50-80 panel", "bbox": [1703, 583, 2470, 1159]},
            {
                "kind": "text",
                "name": "partial data readout",
                "role": "audience",
                "text": "Partial-data tail-match closure status and required regeneration",
                "bbox": [180, 1340, 2380, 1395],
                "font_px": 30,
            },
        ],
    }
    OUT_LAYOUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def write_speaker_script(results: list[dict[str, object]]) -> None:
    scales = {label: float(result["background_scale"]) for result, (_, label) in zip(results, CENTRALITIES)}
    remainders = {label: float(result["signed_remainder"]) for result, (_, label) in zip(results, CENTRALITIES)}
    text = f"""# Au+Au isolation-energy template slide script

This slide shows the reconstructed isolation-energy template comparison separately for 0-20, 20-50, and 50-80 percent Au+Au centrality. It uses the current partial corrected Au+Au data merge and the historical Photon12 plus Photon20 embedded signal sample over 15 to 35 GeV.

The black points are BDT-tight data. The red distribution is the bounded BDT non-tight data sideband, independently scaled in each centrality interval so that its isolation tail above 6 GeV matches the tight-data tail. The blue component, where a physical normalization exists, is the truth-isolated, truth-matched, BDT-tight photon-plus-jet simulation normalized to the signed tight-minus-background remainder.

The non-tight tail scale factors are {scales['0-20%']:.3f}, {scales['20-50%']:.3f}, and {scales['50-80%']:.3f} from central to peripheral collisions. The corresponding signed tight-minus-background density integrals are {remainders['0-20%']:.1f}, {remainders['20-50%']:.1f}, and {remainders['50-80%']:.1f}.

With this partial data merge, the exact PPG12-style tail match leaves positive signal remainders in 0-20 and 20-50 percent. In 50-80 percent, the scaled sideband exceeds the total tight yield, so the slide deliberately does not draw an unphysical negative signal component. This is a diagnostic of template closure, not a final purity result.

The slide must be regenerated from the complete corrected Au+Au data merge before it is used as a final analysis statement.
"""
    OUT_SCRIPT.write_text(text, encoding="utf-8")


def run_audit() -> None:
    subprocess.run(
        [
            sys.executable,
            str(REPO / "scripts/slides/common/post_render_slide_audit.py"),
            "--png",
            str(OUT_PNG),
            "--layout-nodes",
            str(OUT_LAYOUT),
            "--output",
            str(OUT_AUDIT),
            "--require-pass",
        ],
        check=True,
    )


def main() -> int:
    for path in (DATA_ROOT, SIGNAL_ROOT):
        if not path.is_file():
            raise FileNotFoundError(path)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    results = [build_panel(key) for key, _ in CENTRALITIES]
    write_csv(results)
    render(results)
    write_layout()
    write_speaker_script(results)
    manifest = {
        "slide_png": str(OUT_PNG),
        "csv": str(OUT_CSV),
        "speaker_script": str(OUT_SCRIPT),
        "layout_nodes": str(OUT_LAYOUT),
        "audit": str(OUT_AUDIT),
        "data_root": str(DATA_ROOT),
        "data_root_sha256": sha256(DATA_ROOT),
        "data_root_size_bytes": DATA_ROOT.stat().st_size,
        "signal_root": str(SIGNAL_ROOT),
        "signal_root_sha256": sha256(SIGNAL_ROOT),
        "signal_root_size_bytes": SIGNAL_ROOT.stat().st_size,
        "data_status": "preliminary partial corrected AuAu merge; regenerate after exact full raw coverage and final merge",
        "pt_interval": "15 <= E_T_gamma < 35 GeV",
        "eta_selection": "|eta_gamma| < 0.7 from the production configuration",
        "isolation_cone_radius": 0.4,
        "data_tight_definition": "historical pre-MinBiasClassifier AuAuCentInputBase3x3BDT tight selection",
        "data_nontight_definition": "historical bounded AuAuBDTSideband selection",
        "signal_definition": "truth-isolated, truth-matched reconstructed photon satisfying the historical BDT-tight selection",
        "normalization": {
            "rebin": "PPG12 CONF_plots.C grouping: one source bin below 2.5 GeV and five source bins above 2.5 GeV",
            "bin_width": "all distributions divided by rebinned bin width",
            "background": "bounded non-tight data scaled to the tight-data integral for Eiso >= 6 GeV",
            "signal": "truth-matched tight signal MC scaled to the signed full-range tight-minus-scaled-background integral only when positive",
            "guard": "no signal component is drawn when the signed normalization is non-positive",
        },
        "panels": [
            {
                "centrality": label,
                "tight_keys": result["tight_keys"],
                "nontight_keys": result["nontight_keys"],
                "signal_keys": result["signal_keys"],
                "tight_tail": result["tight_tail"],
                "nontight_tail": result["nontight_tail"],
                "background_scale": result["background_scale"],
                "signed_remainder": result["signed_remainder"],
                "signal_scale": result["signal_scale"],
                "physical_signal_normalization": result["physical_signal_normalization"],
            }
            for result, (_, label) in zip(results, CENTRALITIES)
        ],
        "reference_code": "ppg12codeGit/plotting/CONF_plots.C",
        "slide_dimensions_px": [SLIDE_WIDTH_PX, SLIDE_HEIGHT_PX],
        "audience_canvas_label": "No campaign identifier is shown on the slide canvas.",
    }
    OUT_MANIFEST.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    run_audit()
    print(
        json.dumps(
            {
                "slide_png": str(OUT_PNG),
                "csv": str(OUT_CSV),
                "manifest": str(OUT_MANIFEST),
                "speaker_script": str(OUT_SCRIPT),
                "audit": str(OUT_AUDIT),
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
