#!/usr/bin/env python3
"""Build the THE-88 AuAu photon-efficiency stage slide."""

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
from matplotlib.patches import FancyBboxPatch


THIS_FILE = Path(__file__).resolve()
REPO = next((p for p in THIS_FILE.parents if (p / "AGENTS.md").exists()), THIS_FILE.parents[4])
SCRIPTS = REPO / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.append(str(SCRIPTS))

from slides.common.slide_defaults import (  # noqa: E402
    SLIDE_DPI,
    SLIDE_HEIGHT_PX,
    SLIDE_WIDTH_PX,
    slide_figsize,
)


ROOT_INPUT = (
    REPO
    / "InputFiles/the88_bounded_nontight_20260708/simembedded"
    / "preselectionNewPPG12_tightAuAuCentInputBase3x3BDT_nonTightAuAuBDTSideband_baseVariant"
    / "photonJet12and20merged_SIM/RecoilJets_embeddedPhoton12plus20_MERGED.root"
)
OUT_DIR = (
    REPO
    / "dataOutput/auau/the88_auau_data_centrality_fixed_20260713"
    / "photon_efficiency_stages"
)
STEM = "the88_auau_r40_cent_sliding_photon_efficiency_stages_15to35_by_centrality"
OUT_PNG = OUT_DIR / f"{STEM}.png"
OUT_CSV = OUT_DIR / f"{STEM}.csv"
OUT_MANIFEST = OUT_DIR / f"{STEM}.manifest.json"
OUT_LAYOUT = OUT_DIR / f"{STEM}.layout_nodes.json"
OUT_AUDIT = OUT_DIR / f"{STEM}.audit.json"
OUT_SCRIPT = OUT_DIR / f"{STEM}_speaker_script.md"

TOPDIR = "SIM"
ISO_TAG = "isoR40_isSliding"
CENTRALITIES = [
    ("0_20", "0-20%"),
    ("20_50", "20-50%"),
    ("50_80", "50-80%"),
]
STAGES = ["reco", "reco_id", "reco_id_iso"]

COLORS = {
    "reco": "#111111",
    "reco_id": "#d62d91",
    "reco_id_iso": "#238b45",
}
MARKERS = {"reco": "s", "reco_id": "o", "reco_id_iso": "P"}
LABELS = {
    "reco": r"$\varepsilon_{\mathrm{reco}}$",
    "reco_id": r"$\varepsilon_{\mathrm{reco}}\,\varepsilon_{\mathrm{ID}}$",
    "reco_id_iso": r"$\varepsilon_{\mathrm{reco}}\,\varepsilon_{\mathrm{ID}}\,\varepsilon_{\mathrm{iso}}$",
}

INK = "#111827"
MUTED = "#475569"
GRID = "#d9e2ec"
EDGE = "#cbd5e1"
PANEL_FILL = "#fbfdff"

TIMES_DIR = Path("/System/Library/Fonts/Supplemental")
TIMES_FONTS = [
    TIMES_DIR / "Times New Roman.ttf",
    TIMES_DIR / "Times New Roman Bold.ttf",
    TIMES_DIR / "Times New Roman Italic.ttf",
    TIMES_DIR / "Times New Roman Bold Italic.ttf",
]


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
        }
    )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def ratio_error(num: float, den: float, enum: float, eden: float) -> float:
    if den <= 0.0 or num < 0.0:
        return math.nan
    return math.sqrt((enum / den) ** 2 + ((num * eden) / (den * den)) ** 2)


def product_error(a: float, ea: float, b: float, eb: float) -> float:
    if not all(math.isfinite(value) for value in (a, ea, b, eb)):
        return math.nan
    return math.sqrt((b * ea) ** 2 + (a * eb) ** 2)


def require_hist(directory, name: str):
    hist = directory.Get(name) if directory else None
    if not hist:
        raise KeyError(f"Missing histogram: {TOPDIR}/{name}")
    hist.SetDirectory(0)
    return hist


def collect_rows() -> tuple[list[dict], dict]:
    import ROOT

    ROOT.gROOT.SetBatch(True)
    ROOT.TH1.AddDirectory(False)
    handle = ROOT.TFile.Open(str(ROOT_INPUT), "READ")
    if not handle or handle.IsZombie() or handle.TestBit(ROOT.TFile.kRecovered):
        raise OSError(f"Could not open clean ROOT file: {ROOT_INPUT}")
    sim = handle.Get(TOPDIR)
    if not sim:
        raise KeyError(f"Missing directory: {TOPDIR}")

    rows: list[dict] = []
    object_names: dict[str, dict[str, str]] = {}
    native_edges: list[float] | None = None
    try:
        for cent_key, cent_label in CENTRALITIES:
            names = {
                "truth": f"h_photonEffTruthDen_pTgamma_cent_{cent_key}",
                "reco": f"h_photonEffReco_pTgamma_cent_{cent_key}",
                "reco_iso": f"h_photonEffRecoIso_pTgamma_{ISO_TAG}_cent_{cent_key}",
                "tight_iso": f"h_photonEffRecoTightIso_pTgamma_{ISO_TAG}_cent_{cent_key}",
            }
            object_names[cent_key] = {key: f"{TOPDIR}/{value}" for key, value in names.items()}
            h_truth = require_hist(sim, names["truth"])
            h_reco = require_hist(sim, names["reco"])
            h_reco_iso = require_hist(sim, names["reco_iso"])
            h_tight_iso = require_hist(sim, names["tight_iso"])

            edges = [float(h_truth.GetXaxis().GetBinLowEdge(i)) for i in range(1, h_truth.GetNbinsX() + 2)]
            if native_edges is None:
                native_edges = edges
            elif edges != native_edges:
                raise ValueError(f"Truth-pT binning differs for centrality {cent_label}")

            for ib in range(1, h_truth.GetNbinsX() + 1):
                lo = float(h_truth.GetXaxis().GetBinLowEdge(ib))
                hi = float(h_truth.GetXaxis().GetBinUpEdge(ib))
                if lo < 15.0 or hi > 35.0:
                    continue
                mid = 0.5 * (lo + hi)

                def content(hist) -> tuple[float, float]:
                    jb = hist.GetXaxis().FindBin(mid)
                    return float(hist.GetBinContent(jb)), float(hist.GetBinError(jb))

                truth, etruth = content(h_truth)
                reco, ereco = content(h_reco)
                reco_iso, ereco_iso = content(h_reco_iso)
                tight_iso, etight_iso = content(h_tight_iso)
                if truth <= 0.0 or reco_iso <= 0.0:
                    continue

                reco_eff = reco / truth
                reco_eff_err = ratio_error(reco, truth, ereco, etruth)
                id_cond = tight_iso / reco_iso
                id_cond_err = ratio_error(tight_iso, reco_iso, etight_iso, ereco_iso)
                reco_id_eff = reco_eff * id_cond
                reco_id_eff_err = product_error(reco_eff, reco_eff_err, id_cond, id_cond_err)
                reco_id_iso_eff = tight_iso / truth
                reco_id_iso_eff_err = ratio_error(tight_iso, truth, etight_iso, etruth)

                values = {
                    "reco": (reco_eff, reco_eff_err, reco, truth),
                    "reco_id": (reco_id_eff, reco_id_eff_err, tight_iso, reco_iso),
                    "reco_id_iso": (reco_id_iso_eff, reco_id_iso_eff_err, tight_iso, truth),
                }
                for stage in STAGES:
                    value, error, numerator, denominator = values[stage]
                    if not (math.isfinite(value) and math.isfinite(error) and 0.0 <= value <= 1.2):
                        raise ValueError(
                            f"Invalid {stage} efficiency for {cent_label}, {lo:g}-{hi:g} GeV: {value} +/- {error}"
                        )
                    rows.append(
                        {
                            "centrality_key": cent_key,
                            "centrality_label": cent_label,
                            "stage": stage,
                            "pt_low_gev": lo,
                            "pt_high_gev": hi,
                            "pt_mid_gev": mid,
                            "numerator": numerator,
                            "denominator": denominator,
                            "efficiency": value,
                            "efficiency_error": error,
                        }
                    )
    finally:
        handle.Close()

    if len(rows) != len(CENTRALITIES) * len(STAGES) * 6:
        raise RuntimeError(f"Expected 54 efficiency rows, found {len(rows)}")

    metadata = {
        "root_input": str(ROOT_INPUT),
        "root_sha256": sha256(ROOT_INPUT),
        "root_size_bytes": ROOT_INPUT.stat().st_size,
        "top_directory": TOPDIR,
        "object_names": object_names,
        "native_truth_pt_edges_gev": native_edges,
        "display_window_gev": [15.0, 35.0],
        "display_interval": "15 <= pT_gamma_truth < 35 GeV",
        "reco_isolation": {
            "cone_radius": 0.4,
            "mode": "centrality-dependent sliding",
            "threshold": "Eiso < 7.57 - 0.0658*centrality_percentile GeV",
            "histogram_tag": ISO_TAG,
        },
        "truth_denominator": "prompt truth photons passing truth isolation Eiso_truth < 4 GeV",
        "photon_id": {
            "tight": "AuAuCentInputBase3x3BDT centlinear threshold 0.5387310379 + 0.0011102647*centrality",
            "nontight": "AuAuBDTSideband bounded relative-to-tight offsets [-0.20,-0.03]; not used by these tight-efficiency curves",
        },
        "embedded_minbias_classifier": {
            "required": False,
            "scope": "historical pre-MinBiasClassifier THE-88 embedded baseline",
            "note": "The standard setMinBiasClassifer YAML field does not apply the separate embedded-event filter.",
        },
        "stage_definitions": {
            "reco": "N_reco / N_truth",
            "reco_id": "(N_reco / N_truth) * (N_tight_and_iso / N_iso)",
            "reco_id_iso": "N_tight_and_iso / N_truth",
        },
        "production_artifact_exception": (
            "Explicit historical THE-88 embedded-signal diagnostic requested by Justin; "
            "the dated ROOT is used intentionally instead of a current pointer."
        ),
    }
    return rows, metadata


def stage_arrays(rows: list[dict], cent_key: str, stage: str) -> tuple[np.ndarray, ...]:
    selected = sorted(
        (row for row in rows if row["centrality_key"] == cent_key and row["stage"] == stage),
        key=lambda row: row["pt_mid_gev"],
    )
    return (
        np.array([row["pt_mid_gev"] for row in selected], dtype=float),
        np.array([0.5 * (row["pt_high_gev"] - row["pt_low_gev"]) for row in selected], dtype=float),
        np.array([row["efficiency"] for row in selected], dtype=float),
        np.array([row["efficiency_error"] for row in selected], dtype=float),
    )


def draw_card(fig: plt.Figure, bbox: tuple[float, float, float, float], stage: str, formula: str, body: str) -> None:
    x, y, w, h = bbox
    card = FancyBboxPatch(
        (x, y),
        w,
        h,
        transform=fig.transFigure,
        boxstyle="round,pad=0.006,rounding_size=0.012",
        facecolor="white",
        edgecolor=COLORS[stage],
        linewidth=2.5,
        zorder=2,
    )
    fig.add_artist(card)
    fig.text(x + 0.019, y + h * 0.76, LABELS[stage], fontsize=24, fontweight="bold", color=COLORS[stage], va="center")
    fig.text(x + 0.019, y + h * 0.44, formula, fontsize=17.6, color=INK, va="center")
    fig.text(x + 0.019, y + h * 0.12, body, fontsize=14.4, color=MUTED, va="center")


def render(rows: list[dict]) -> None:
    configure_fonts()
    fig = plt.figure(figsize=slide_figsize(), dpi=SLIDE_DPI, facecolor="white")
    fig.text(
        0.049,
        0.946,
        "Au+Au photon efficiency by centrality",
        fontsize=43,
        fontweight="bold",
        color=INK,
        ha="left",
        va="top",
    )
    fig.text(
        0.050,
        0.840,
        r"Embedded photon+jet 12 + 20 production slices; $15\leq p_T^{\gamma,\mathrm{truth}}<35$ GeV",
        fontsize=22.0,
        color=INK,
        ha="left",
        va="center",
    )

    card_y, card_h, card_w = 0.672, 0.140, 0.282
    draw_card(
        fig,
        (0.050, card_y, card_w, card_h),
        "reco",
        r"$N_{\mathrm{reco}}\,/\,N_{\mathrm{truth\ iso}}$",
        "Reco-matched / truth-isolated",
    )
    draw_card(
        fig,
        (0.359, card_y, card_w, card_h),
        "reco_id",
        r"$\varepsilon_{\mathrm{reco}}\times(N_{\mathrm{tight+iso}}/N_{\mathrm{iso}})$",
        "Historical tight-ID working point",
    )
    draw_card(
        fig,
        (0.668, card_y, card_w, card_h),
        "reco_id_iso",
        r"$N_{\mathrm{tight+iso}}\,/\,N_{\mathrm{truth\ iso}}$",
        "Tight + isolated / truth-isolated",
    )

    fig.text(
        0.500,
        0.628,
        r"Truth: $R=0.3$, $E_{T,\mathrm{truth}}^{\mathrm{iso}}<4$ GeV   |   Reco: $R=0.4$, $E_T^{\mathrm{iso}}<7.57-0.0658c$ GeV, $c=$ event centrality percentile",
        fontsize=17.2,
        color=MUTED,
        ha="center",
        va="center",
    )

    lefts = [0.063, 0.365, 0.667]
    axes = [fig.add_axes([left, 0.122, 0.278, 0.425]) for left in lefts]
    for index, ((cent_key, cent_label), ax) in enumerate(zip(CENTRALITIES, axes)):
        ax.set_facecolor(PANEL_FILL)
        for stage in STAGES:
            x, ex, y, ey = stage_arrays(rows, cent_key, stage)
            ax.errorbar(
                x,
                y,
                xerr=ex,
                yerr=ey,
                fmt=MARKERS[stage],
                color=COLORS[stage],
                markerfacecolor=COLORS[stage],
                markeredgecolor=COLORS[stage],
                markeredgewidth=0.9,
                markersize=7.2,
                elinewidth=1.15,
                capsize=0,
                linestyle="-",
                linewidth=1.45,
                zorder=4,
            )
        ax.set_title(f"Au+Au {cent_label}", fontsize=25, fontweight="bold", color=INK, pad=12)
        ax.set_xlim(15.0, 35.0)
        ax.set_ylim(0.0, 1.05)
        ax.set_xticks([15, 20, 25, 30, 35])
        ax.set_yticks(np.arange(0.0, 1.01, 0.2))
        ax.grid(axis="y", color=GRID, linewidth=1.0, zorder=0)
        ax.tick_params(direction="in", which="both", top=True, right=True, labelsize=17, length=6, width=1.05)
        ax.minorticks_on()
        ax.tick_params(which="minor", length=3.5, width=0.8)
        for spine in ax.spines.values():
            spine.set_color("#64748b")
            spine.set_linewidth(1.1)
        ax.set_xlabel(r"$p_T^{\gamma,\mathrm{truth}}$ [GeV]", fontsize=21, labelpad=8)
        if index == 0:
            ax.set_ylabel("Efficiency", fontsize=22, labelpad=12)
        else:
            ax.tick_params(labelleft=False)
        if index == 0:
            ax.text(
                0.055,
                0.305,
                r"$\it{\bf{sPHENIX}}$ Internal",
                transform=ax.transAxes,
                fontsize=16.5,
                color=INK,
                va="top",
            )
            ax.text(
                0.055,
                0.245,
                r"embedded PYTHIA8, $\sqrt{s_{NN}}=200$ GeV",
                transform=ax.transAxes,
                fontsize=13.4,
                color=INK,
                va="top",
            )
            ax.text(
                0.055,
                0.190,
                r"$|\eta^{\gamma,\mathrm{truth}}|<0.7$",
                transform=ax.transAxes,
                fontsize=13.4,
                color=INK,
                va="top",
            )

    legend_handles = [
        Line2D(
            [0],
            [0],
            color=COLORS[stage],
            marker=MARKERS[stage],
            markerfacecolor=COLORS[stage],
            lw=1.45,
            ms=7,
            label=LABELS[stage],
        )
        for stage in STAGES
    ]
    axes[1].legend(
        handles=legend_handles,
        loc="lower left",
        bbox_to_anchor=(0.03, 0.03),
        frameon=False,
        fontsize=15.5,
        handlelength=1.4,
        labelspacing=0.42,
    )

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_PNG, dpi=SLIDE_DPI, facecolor="white")
    plt.close(fig)


def write_csv(rows: list[dict]) -> None:
    with OUT_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_layout() -> None:
    payload = {
        "schema": "slide_layout_nodes_v1",
        "slide_size": [SLIDE_WIDTH_PX, SLIDE_HEIGHT_PX],
        "title_axis_x": 126,
        "minimum_title_font_px": 68,
        "minimum_audience_font_px": 27,
        "nodes": [
            {
                "kind": "text",
                "name": "title",
                "role": "title",
                "title_anchor": True,
                "text": "Au+Au photon efficiency by centrality",
                "bbox": [126, 66, 2200, 145],
                "font_px": 82,
            },
            {
                "kind": "text",
                "name": "sample line",
                "role": "audience",
                "title_axis_align": "left",
                "text": "Embedded photon+jet 12 + 20 production slices; 15 <= truth pT < 35 GeV",
                "bbox": [128, 210, 2380, 262],
                "font_px": 45,
            },
            {"kind": "panel", "name": "reco definition", "bbox": [128, 271, 850, 472]},
            {"kind": "panel", "name": "reco id definition", "bbox": [919, 271, 1640, 472]},
            {"kind": "panel", "name": "full definition", "bbox": [1710, 271, 2432, 472]},
            {
                "kind": "text",
                "name": "selection line",
                "role": "audience",
                "text": "Truth isolation, centrality-dependent reco isolation, and historical BDT-tight definitions",
                "bbox": [150, 504, 2410, 558],
                "font_px": 35,
            },
            {"kind": "panel", "name": "0-20 panel", "bbox": [161, 652, 873, 1269]},
            {"kind": "panel", "name": "20-50 panel", "bbox": [934, 652, 1646, 1269]},
            {"kind": "panel", "name": "50-80 panel", "bbox": [1708, 652, 2420, 1269]},
        ],
    }
    OUT_LAYOUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def write_speaker_script(rows: list[dict]) -> None:
    ranges: dict[str, dict[str, tuple[float, float]]] = {}
    for cent_key, cent_label in CENTRALITIES:
        ranges[cent_label] = {}
        for stage in STAGES:
            values = [
                row["efficiency"]
                for row in rows
                if row["centrality_key"] == cent_key and row["stage"] == stage
            ]
            ranges[cent_label][stage] = (min(values), max(values))

    text = f"""# Au+Au photon efficiency stages slide script

This slide shows the three cumulative photon-efficiency stages in the Au+Au embedded photon-plus-jet sample. From left to right, the panels are 0-20, 20-50, and 50-80 percent centrality, and every panel uses truth-photon bins from 15 GeV inclusive to 35 GeV exclusive.

The black curve is reconstruction efficiency: matched reconstructed photons divided by truth-isolated prompt photons. The magenta curve adds the nominal Au+Au BDT tight-identification efficiency using the same PPG12-style conditional construction. The green curve is the full tight-and-isolated numerator divided by the truth-isolated denominator.

The reconstructed isolation requirement is the nominal R equals 0.4 centrality-dependent window, E T iso less than 7.57 minus 0.0658 times the event centrality percentile. The truth denominator uses R equals 0.3 and a fixed truth-isolation requirement below 4 GeV. This historical THE-88 embedded sample predates the separate embedded MinimumBiasClassifier requirement; it uses the nominal AuAuCentInputBase3x3BDT centrality-dependent tight working point.

Across the displayed range, the full green efficiency spans {ranges['0-20%']['reco_id_iso'][0]:.3f} to {ranges['0-20%']['reco_id_iso'][1]:.3f} in 0-20 percent, {ranges['20-50%']['reco_id_iso'][0]:.3f} to {ranges['20-50%']['reco_id_iso'][1]:.3f} in 20-50 percent, and {ranges['50-80%']['reco_id_iso'][0]:.3f} to {ranges['50-80%']['reco_id_iso'][1]:.3f} in 50-80 percent.

The main point is that the same efficiency definition is evaluated independently in each centrality interval, so the separation among the panels directly shows the centrality dependence of reconstruction, tight identification, and isolation.
"""
    OUT_SCRIPT.write_text(text, encoding="utf-8")


def run_audit() -> None:
    command = [
        sys.executable,
        str(REPO / "scripts/slides/common/post_render_slide_audit.py"),
        "--png",
        str(OUT_PNG),
        "--layout-nodes",
        str(OUT_LAYOUT),
        "--output",
        str(OUT_AUDIT),
        "--require-pass",
    ]
    subprocess.run(command, check=True)


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows, metadata = collect_rows()
    write_csv(rows)
    render(rows)
    write_layout()
    write_speaker_script(rows)
    metadata.update(
        {
            "slide_png": str(OUT_PNG),
            "csv": str(OUT_CSV),
            "layout_nodes": str(OUT_LAYOUT),
            "audit": str(OUT_AUDIT),
            "speaker_script": str(OUT_SCRIPT),
            "slide_dimensions_px": [SLIDE_WIDTH_PX, SLIDE_HEIGHT_PX],
            "audience_canvas_label": "No THE-88 campaign identifier is shown on the slide canvas.",
        }
    )
    OUT_MANIFEST.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
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
