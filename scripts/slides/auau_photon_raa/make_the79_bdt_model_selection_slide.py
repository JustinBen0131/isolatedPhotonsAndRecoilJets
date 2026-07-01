#!/usr/bin/env python3
"""Build a THE-79 model-selection summary slide candidate."""

from __future__ import annotations

import csv
import json
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch, Rectangle


THIS_FILE = Path(__file__).resolve()
REPO = next((p for p in THIS_FILE.parents if (p / "AGENTS.md").exists()), THIS_FILE.parents[4])
SCRIPTS = REPO / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.append(str(SCRIPTS))

from slides.common.slide_defaults import SLIDE_DPI, SLIDE_HEIGHT_PX, SLIDE_WIDTH_PX, slide_figsize  # noqa: E402


OUT_DIR = REPO / "dataOutput/auauBDTModelSelection/THE79_5to40_candidate_summary_20260627"
OUT_PNG = OUT_DIR / "the79_bdt_model_selection_summary_slide.png"
OUT_SCRIPT = OUT_DIR / "the79_bdt_model_selection_summary_speaker_script.md"
OUT_CSV = OUT_DIR / "the79_bdt_model_selection_rankings.csv"
OUT_MANIFEST = OUT_DIR / "the79_bdt_model_selection_summary_manifest.json"
OUT_LAYOUT = OUT_DIR / "the79_bdt_model_selection_summary_layout_nodes.json"


@dataclass(frozen=True)
class Candidate:
    label: str
    product: str
    family: str
    auc: float
    auc_pt_mean: float
    auc_pt_cent_mean: float
    eligible_entries: int
    signal_mean: float
    background_mean: float
    selected: bool = False


CANDIDATES = [
    Candidate(
        label="pT × cent7",
        product="base14_perEtCent7",
        family="14-feature BDT, individual pT × 7 centrality bins",
        auc=0.8760020160331478,
        auc_pt_mean=0.8743841759379238,
        auc_pt_cent_mean=0.8723981333400556,
        eligible_entries=7_306_271,
        signal_mean=0.6652458906173706,
        background_mean=0.25469082593917847,
        selected=True,
    ),
    Candidate(
        label="pT × cent3",
        product="base14_perEtCent3",
        family="14-feature BDT, individual pT × 3 centrality bins",
        auc=0.8714949958139812,
        auc_pt_mean=0.8686390122739871,
        auc_pt_cent_mean=0.8666791675947387,
        eligible_entries=7_306_271,
        signal_mean=0.6615225076675415,
        background_mean=0.2570413649082184,
    ),
    Candidate(
        label="pT only",
        product="base14_perEt",
        family="14-feature BDT, individual pT bins",
        auc=0.8700753934466987,
        auc_pt_mean=0.8665956296706456,
        auc_pt_cent_mean=0.8624898887529775,
        eligible_entries=8_534_816,
        signal_mean=0.6608784794807434,
        background_mean=0.2565280497074127,
    ),
    Candidate(
        label="Global 15-40",
        product="centAsFeatBase3x3_pt15to40",
        family="single 14-feature BDT, centrality as feature",
        auc=0.8670403147378696,
        auc_pt_mean=0.8631294540607115,
        auc_pt_cent_mean=0.8590458289963268,
        eligible_entries=8_534_816,
        signal_mean=0.6557590365409851,
        background_mean=0.25613951683044434,
    ),
    Candidate(
        label="Global 5-40",
        product="centAsFeatBase3x3_pt5to40",
        family="single broad-range 14-feature BDT",
        auc=0.7891183693480203,
        auc_pt_mean=0.8662697895478007,
        auc_pt_cent_mean=0.8614777441331537,
        eligible_entries=27_424_648,
        signal_mean=0.6597945690155029,
        background_mean=0.3405023217201233,
    ),
]

SELECTED_HIGH_PT_AUC = [
    ("0-20%", "15-20", 0.866352022),
    ("0-20%", "20-25", 0.849490365),
    ("0-20%", "25-35", 0.859745556),
    ("20-50%", "15-20", 0.884749480),
    ("20-50%", "20-25", 0.863794231),
    ("20-50%", "25-35", 0.877582756),
    ("50-80%", "15-20", 0.895640230),
    ("50-80%", "20-25", 0.871040514),
    ("50-80%", "25-35", 0.883188046),
]


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "dejavuserif",
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "axes.facecolor": "white",
            "axes.edgecolor": "#1F2937",
            "axes.linewidth": 1.0,
            "axes.labelsize": 15,
            "axes.titlesize": 17,
            "xtick.labelsize": 12.5,
            "ytick.labelsize": 13.5,
        }
    )


def add_box(
    fig: plt.Figure,
    box: tuple[float, float, float, float],
    *,
    edge: str = "#CBD5E1",
    face: str = "white",
    lw: float = 1.4,
    radius: float = 0.018,
) -> plt.Axes:
    ax = fig.add_axes(box)
    ax.set_axis_off()
    patch = FancyBboxPatch(
        (0, 0),
        1,
        1,
        boxstyle=f"round,pad=0.012,rounding_size={radius}",
        transform=ax.transAxes,
        facecolor=face,
        edgecolor=edge,
        linewidth=lw,
    )
    ax.add_patch(patch)
    return ax


def add_metric(ax: plt.Axes, x: float, y: float, value: str, label: str, color: str = "#0F766E") -> None:
    ax.text(x, y, value, transform=ax.transAxes, fontsize=20.5, fontweight="bold", color=color, va="baseline")
    ax.text(x, y - 0.13, label, transform=ax.transAxes, fontsize=11.7, color="#334155", va="baseline")


def add_title_and_context(fig: plt.Figure) -> None:
    fig.text(
        0.047,
        0.942,
        "Binned pT × centrality BDT wins for AuAu data running",
        fontsize=28.5,
        fontweight="bold",
        color="#0F172A",
        va="top",
    )
    fig.text(
        0.047,
        0.883,
        "Full-stat embedded-MC validation on 27.4M scored rows; all candidates use the 14-feature AuAu baseline.",
        fontsize=15.6,
        color="#334155",
        va="top",
    )


def plot_auc_bars(fig: plt.Figure) -> None:
    ax = fig.add_axes([0.155, 0.565, 0.47, 0.255])
    candidates = list(reversed(CANDIDATES))
    y = np.arange(len(candidates))
    colors = ["#B91C1C" if c.product.endswith("pt5to40") else "#64748B" for c in candidates]
    colors[-1] = "#0F766E"
    ax.barh(y, [c.auc for c in candidates], color=colors, height=0.62)
    ax.set_yticks(y)
    ax.set_yticklabels([c.label for c in candidates])
    ax.set_xlim(0.775, 0.885)
    ax.set_xlabel("")
    ax.set_title("Selection metric - inclusive AUC ranking", loc="left", fontweight="bold", pad=8)
    ax.grid(axis="x", color="#E2E8F0", linewidth=1.0)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    for yi, cand in zip(y, candidates):
        ax.text(cand.auc + 0.0020, yi, f"{cand.auc:.3f}", va="center", fontsize=12.5, color="#0F172A")
    ax.text(
        0.876002016,
        y[-1] + 0.35,
        "selected",
        ha="center",
        va="bottom",
        fontsize=12.5,
        color="#0F766E",
        fontweight="bold",
    )


def plot_background_bars(fig: plt.Figure) -> None:
    ax = fig.add_axes([0.155, 0.295, 0.47, 0.205])
    ordered = [CANDIDATES[0], CANDIDATES[3], CANDIDATES[2], CANDIDATES[1], CANDIDATES[4]]
    y = np.arange(len(ordered))
    colors = ["#0F766E"] + ["#64748B", "#64748B", "#64748B", "#B91C1C"]
    ax.barh(y, [c.background_mean for c in reversed(ordered)], color=list(reversed(colors)), height=0.58)
    labels = [c.label for c in reversed(ordered)]
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlim(0.24, 0.35)
    ax.set_xlabel("Mean background BDT score; lower is cleaner")
    ax.set_title("Global 5-40 underperforms by lifting background scores", loc="left", fontweight="bold", pad=8)
    ax.grid(axis="x", color="#E2E8F0", linewidth=1.0)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    for yi, cand in zip(y, reversed(ordered)):
        ax.text(cand.background_mean + 0.0022, yi, f"{cand.background_mean:.3f}", va="center", fontsize=12.5, color="#0F172A")


def draw_result_card(fig: plt.Figure) -> None:
    ax = add_box(fig, [0.675, 0.565, 0.275, 0.255], edge="#99F6E4", face="#F8FFFD", lw=1.7)
    ax.add_patch(Rectangle((0.0, 0.0), 0.018, 1.0, transform=ax.transAxes, facecolor="#0F766E", edgecolor="none"))
    ax.text(0.065, 0.86, "What worked best", fontsize=20.5, fontweight="bold", color="#0F172A", va="top")
    ax.text(
        0.065,
        0.68,
        "Train the 14-feature baseline in\nindividual pT and centrality bins.",
        fontsize=15.5,
        color="#1E293B",
        va="top",
        linespacing=1.2,
    )
    ax.text(
        0.065,
        0.47,
        "Selected product - base14_perEtCent7",
        fontsize=14.0,
        fontweight="bold",
        color="#0F766E",
        va="top",
    )
    add_metric(ax, 0.065, 0.245, "0.876", "inclusive AUC")
    add_metric(ax, 0.405, 0.245, "0.872", "pT × cent AUC")
    add_metric(ax, 0.71, 0.245, "0.255", "bkg mean")


def draw_high_pt_matrix(fig: plt.Figure) -> None:
    ax = add_box(fig, [0.675, 0.295, 0.275, 0.205], edge="#CBD5E1", face="white", lw=1.3)
    ax.text(0.055, 0.88, "Selected model high-pT AUCs", transform=ax.transAxes, fontsize=16.5, fontweight="bold", color="#0F172A", va="top")
    ax.text(0.055, 0.735, "rows: centrality; columns: pT [GeV]", transform=ax.transAxes, fontsize=11.5, color="#475569", va="top")
    cents = ["0-20%", "20-50%", "50-80%"]
    pts = ["15-20", "20-25", "25-35"]
    values = {(c, p): v for c, p, v in SELECTED_HIGH_PT_AUC}
    x0, y0, cell_w, cell_h = 0.30, 0.13, 0.205, 0.14
    for j, pt in enumerate(pts):
        ax.text(x0 + j * cell_w + cell_w / 2, y0 + 3 * cell_h + 0.030, pt, ha="center", va="bottom", fontsize=10.8, color="#334155")
    for i, cent in enumerate(cents):
        y = y0 + (2 - i) * cell_h
        ax.text(0.06, y + cell_h / 2, cent, ha="left", va="center", fontsize=11.2, color="#334155")
        for j, pt in enumerate(pts):
            value = values[(cent, pt)]
            intensity = (value - 0.84) / (0.90 - 0.84)
            intensity = min(1.0, max(0.0, intensity))
            color = (
                231 - int(95 * intensity),
                245 - int(85 * intensity),
                241 - int(80 * intensity),
            )
            rect = Rectangle((x0 + j * cell_w, y), cell_w * 0.93, cell_h * 0.82, transform=ax.transAxes, facecolor=np.array(color) / 255.0, edgecolor="#CBD5E1", linewidth=0.8)
            ax.add_patch(rect)
            ax.text(x0 + j * cell_w + cell_w * 0.465, y + cell_h * 0.41, f"{value:.3f}", ha="center", va="center", fontsize=11.0, fontweight="bold", color="#0F172A")


def draw_bottom_decision(fig: plt.Figure) -> None:
    ax = add_box(fig, [0.047, 0.075, 0.903, 0.115], edge="#CBD5E1", face="#FFFFFF", lw=1.4, radius=0.014)
    ax.add_patch(Rectangle((0.023, 0.18), 0.018, 0.64, transform=ax.transAxes, facecolor="#0F766E", edgecolor="none"))
    ax.text(0.065, 0.65, "Decision", transform=ax.transAxes, fontsize=18.5, fontweight="bold", color="#0F172A", va="center")
    ax.text(
        0.19,
        0.66,
        "Use binned pT × cent7 as the primary AuAu data-running BDT.",
        transform=ax.transAxes,
        fontsize=17.5,
        color="#0F172A",
        va="center",
    )
    ax.text(
        0.19,
        0.32,
        "Keep global 15-40 as a simpler control; keep global 5-40 as a broad-range diagnostic only.",
        transform=ax.transAxes,
        fontsize=15.2,
        color="#334155",
        va="center",
    )


def write_rankings_csv() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with OUT_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "label",
                "product",
                "family",
                "inclusive_auc",
                "auc_pt_mean",
                "auc_pt_cent_mean",
                "eligible_entries",
                "signal_score_mean",
                "background_score_mean",
                "selected",
            ],
        )
        writer.writeheader()
        for candidate in CANDIDATES:
            writer.writerow(
                {
                    "label": candidate.label,
                    "product": candidate.product,
                    "family": candidate.family,
                    "inclusive_auc": candidate.auc,
                    "auc_pt_mean": candidate.auc_pt_mean,
                    "auc_pt_cent_mean": candidate.auc_pt_cent_mean,
                    "eligible_entries": candidate.eligible_entries,
                    "signal_score_mean": candidate.signal_mean,
                    "background_score_mean": candidate.background_mean,
                    "selected": candidate.selected,
                }
            )


def write_speaker_script() -> None:
    OUT_SCRIPT.write_text(
        """This slide is the punchline from the THE-79 model-selection validation.

The model type that worked best was not the single global 5-to-40 model. It was the binned fourteen-feature baseline trained in individual pT and centrality bins. The winning product is `base14_perEtCent7`.

The evidence is the full-stat embedded-MC validation. `base14_perEtCent7` has an inclusive AUC of 0.876, a pT-mean AUC of 0.874, and a pT-by-centrality mean AUC of 0.872. It also keeps the mean background score low, about 0.255.

The simpler global 15-to-40 model is close enough to keep as a control, with AUC 0.867. The global 5-to-40 model should not be the data-running choice. Its per-bin averages are not terrible, but inclusively it collapses to 0.789 AUC and pushes the background mean score up to 0.341. That is exactly the failure mode we were worried about from forcing one broad model to cover the whole 5-to-40 range.

So the next data-running proposal should use `base14_perEtCent7` as the primary BDT, keep global 15-to-40 as the simpler comparison, and leave global 5-to-40 as diagnostic evidence only.
""",
        encoding="utf-8",
    )


def write_manifest() -> None:
    selected = CANDIDATES[0]
    OUT_MANIFEST.write_text(
        json.dumps(
            {
                "status": "READY",
                "slide_png": str(OUT_PNG),
                "speaker_script": str(OUT_SCRIPT),
                "rankings_csv": str(OUT_CSV),
                "layout_nodes": str(OUT_LAYOUT),
                "selected_product": selected.product,
                "selected_model_type": selected.family,
                "selected_metrics": {
                    "inclusive_auc": selected.auc,
                    "auc_pt_mean": selected.auc_pt_mean,
                    "auc_pt_cent_mean": selected.auc_pt_cent_mean,
                    "finite_score_fraction": 1.0,
                    "eligible_entries": selected.eligible_entries,
                    "signal_score_mean": selected.signal_mean,
                    "background_score_mean": selected.background_mean,
                },
                "source_reports": {
                    "binned14": "/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_the79_true5to40_defaultsrc_20260626_0015/reports/model_validation_condor_the79_fullsrc_binned14_15to40_simval_20260627",
                    "global15": "/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_the79_true5to40_defaultsrc_20260626_0015/reports/model_validation_condor_the79_fullsrc_global14_pt15to40_simval_20260627",
                    "global5": "/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_the79_true5to40_defaultsrc_20260626_0015/reports/model_validation_condor_the79_fullsrc_global14_pt5to40_simval_20260627",
                },
                "validation_scope": {
                    "scored_entries": 27_424_648,
                    "source_tag": "the79_true5to40_defaultsrc_20260626_0015",
                    "feature_set": "baseV3E + centrality + weta33/wphi33",
                    "systematics_included": False,
                    "auau_data_running_launched": False,
                },
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )


def write_layout_nodes() -> None:
    def node(name: str, text: str, bbox: list[int], font_px: int, *, role: str = "audience", title_axis: bool = False) -> dict:
        item = {
            "kind": "text",
            "name": name,
            "role": role,
            "text": text,
            "bbox": bbox,
            "font_px": font_px,
        }
        if title_axis:
            item["title_axis_align"] = "left"
        if role == "title":
            item["title_anchor"] = True
        return item

    payload = {
        "slide_size": [SLIDE_WIDTH_PX, SLIDE_HEIGHT_PX],
        "title_axis_x": 120,
        "minimum_title_font_px": 66,
        "minimum_audience_font_px": 36,
        "nodes": [
            node(
                "title",
                "Binned pT × centrality BDT wins for AuAu data running",
                [120, 74, 1900, 130],
                79,
                role="title",
            ),
            node(
                "validation context",
                "Full-stat embedded-MC validation on 27.4M scored rows; all candidates use the 14-feature AuAu baseline.",
                [120, 168, 1820, 207],
                43,
                title_axis=True,
            ),
            node(
                "auc chart title",
                "Selection metric - inclusive AUC ranking",
                [396, 286, 1240, 330],
                47,
            ),
            node(
                "background chart title",
                "Global 5-40 underperforms by lifting background scores",
                [396, 713, 1300, 755],
                47,
            ),
            node(
                "what worked best",
                "What worked best",
                [1728, 322, 2320, 370],
                57,
            ),
            node(
                "method statement",
                "Train the 14-feature baseline in individual pT and centrality bins.",
                [1728, 432, 2385, 510],
                43,
            ),
            node(
                "selected product",
                "Selected product - base14_perEtCent7",
                [1728, 612, 2360, 650],
                39,
            ),
            node(
                "high pt matrix title",
                "Selected model high-pT AUCs",
                [1728, 742, 2380, 785],
                46,
            ),
            node(
                "decision label",
                "Decision",
                [166, 990, 440, 1038],
                51,
            ),
            node(
                "decision text",
                "Use binned pT × cent7 as the primary AuAu data-running BDT.",
                [486, 980, 1800, 1025],
                49,
            ),
            node(
                "fallback text",
                "Keep global 15-40 as a simpler control; keep global 5-40 as a broad-range diagnostic only.",
                [486, 1039, 2175, 1080],
                42,
            ),
        ],
    }
    OUT_LAYOUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def make_slide() -> None:
    setup_style()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=slide_figsize(), dpi=SLIDE_DPI)
    add_title_and_context(fig)
    plot_auc_bars(fig)
    plot_background_bars(fig)
    draw_result_card(fig)
    draw_high_pt_matrix(fig)
    draw_bottom_decision(fig)
    fig.savefig(OUT_PNG, dpi=SLIDE_DPI)
    plt.close(fig)


def main() -> int:
    make_slide()
    write_rankings_csv()
    write_speaker_script()
    write_manifest()
    write_layout_nodes()
    print(json.dumps({"slide_png": str(OUT_PNG), "speaker_script": str(OUT_SCRIPT), "manifest": str(OUT_MANIFEST)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
