#!/usr/bin/env python3
"""Opening scope slide for the current default-AuAu-BDT data sample.

The slide deliberately separates run/segment coverage, which comes from the
THE-69 production audit, from event-level counters, which come from the merged
DATA ROOT `cnt_*` histograms.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import uproot
from matplotlib.patches import FancyBboxPatch, Rectangle


REPO = Path(__file__).resolve().parents[3]
if str(REPO / "scripts") not in sys.path:
    sys.path.append(str(REPO / "scripts"))

from slides.common.slide_defaults import SLIDE_DPI, slide_figsize


DATA_ROOT = REPO / (
    "InputFiles/the69_default_auau_physicsqa/"
    "RecoilJets_auau_ALL_preselectionNewPPG12_tightAuAuCentInputBase3x3BDT_"
    "nonTightAuAuBDTComplement_baseVariant.root"
)
OUT_DIR = REPO / "dataOutput/the85_auau_xjgamma_unfolding_push/opening_scope_slide"
OUT_DIR.mkdir(parents=True, exist_ok=True)

PNG_OUT = OUT_DIR / "slide00_auau_data_scope_opening.png"
MANIFEST_OUT = OUT_DIR / "slide00_auau_data_scope_opening_manifest.json"
SCRIPT_OUT = OUT_DIR / "slide00_auau_data_scope_opening_speaker_script.md"


RUN_AUDIT = {
    "grl_runs_requested": 884,
    "analyzable_paired_runs": 882,
    "matched_segment_pairs": 176_577,
    "missing_both_stream_runs": [68491, 72588],
    "pairing_status": "OK_PARTIAL_RUN_COVERAGE / OK_FOR_AVAILABLE_OVERLAP",
    "raw_roots_total": 25_599,
    "raw_roots_non_tiny": 25_481,
    "raw_roots_tiny_or_low_size": 118,
    "raw_bytes": 243_933_291_369,
    "final_root_size_bytes": 136_496_289,
    "data_cluster": 5_552_123,
}

COUNTER_KEYS = {
    "MBD N/S >= 2, |z_vtx| < 150 cm": (
        "MBD_NS_geq_2_vtx_lt_150/cnt_MBD_NS_geq_2_vtx_lt_150"
    ),
    "photon-10 + MBD/vertex": (
        "photon_10_plus_MBD_NS_geq_2_vtx_lt_150/"
        "cnt_photon_10_plus_MBD_NS_geq_2_vtx_lt_150"
    ),
    "photon-12 + MBD/vertex": (
        "photon_12_plus_MBD_NS_geq_2_vtx_lt_150/"
        "cnt_photon_12_plus_MBD_NS_geq_2_vtx_lt_150"
    ),
}

VERTEX_KEYS = {
    "all": "MBD_NS_geq_2_vtx_lt_150/h_vertexZ",
    "0-20%": "MBD_NS_geq_2_vtx_lt_150/h_vertexZ_cent_0_20",
    "20-50%": "MBD_NS_geq_2_vtx_lt_150/h_vertexZ_cent_20_50",
    "50-80%": "MBD_NS_geq_2_vtx_lt_150/h_vertexZ_cent_50_80",
}
CENTRALITY_KEY = "MBD_NS_geq_2_vtx_lt_150/h_centrality"


COLORS = {
    "ink": "#111827",
    "muted": "#526071",
    "grid": "#d6dde6",
    "blue": "#2878a8",
    "green": "#50715f",
    "orange": "#9a6a22",
    "rose": "#cc4c6d",
    "panel": "#ffffff",
    "panel_edge": "#cbd5e1",
    "soft_blue": "#edf6ff",
    "soft_green": "#ecfdf3",
    "soft_rose": "#fff1f3",
}


def hist_values(f: uproot.ReadOnlyDirectory, key: str) -> tuple[np.ndarray, np.ndarray]:
    obj = f[key]
    values, edges = obj.to_numpy()
    return np.asarray(values, dtype=float), np.asarray(edges, dtype=float)


def hist_sum(f: uproot.ReadOnlyDirectory, key: str) -> float:
    values, _ = hist_values(f, key)
    return float(np.sum(values))


def hist_entries(f: uproot.ReadOnlyDirectory, key: str) -> float | None:
    obj = f[key]
    try:
        return float(obj.member("fEntries"))
    except Exception:
        return None


def hist_mean(f: uproot.ReadOnlyDirectory, key: str) -> float:
    values, edges = hist_values(f, key)
    centers = 0.5 * (edges[:-1] + edges[1:])
    denom = np.sum(values)
    return float(np.sum(values * centers) / denom) if denom > 0 else math.nan


def fmt_big(value: float, digits: int = 3) -> str:
    value = float(value)
    abs_value = abs(value)
    if abs_value >= 1e9:
        return f"{value / 1e9:.{digits}f}B"
    if abs_value >= 1e6:
        return f"{value / 1e6:.1f}M"
    if abs_value >= 1e3:
        return f"{value / 1e3:.1f}k"
    return f"{value:.0f}"


def add_box(
    ax: plt.Axes,
    x: float,
    y: float,
    w: float,
    h: float,
    accent: str,
    fill: str = "#ffffff",
    radius: float = 0.018,
) -> None:
    card = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle=f"round,pad=0.012,rounding_size={radius}",
        linewidth=1.35,
        edgecolor=COLORS["panel_edge"],
        facecolor=fill,
        transform=ax.transAxes,
        zorder=1,
    )
    ax.add_patch(card)
    ax.add_patch(
        Rectangle((x, y), 0.009, h, transform=ax.transAxes, linewidth=0, facecolor=accent, zorder=2)
    )


def add_stat_card(
    ax: plt.Axes,
    x: float,
    y: float,
    w: float,
    h: float,
    label: str,
    value: str,
    sublabel: str,
    accent: str,
    fill: str,
) -> None:
    add_box(ax, x, y, w, h, accent, fill)
    ax.text(
        x + 0.030,
        y + h - 0.036,
        label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=16.5,
        fontweight="bold",
        color=COLORS["ink"],
    )
    ax.text(
        x + 0.030,
        y + h * 0.36,
        value,
        transform=ax.transAxes,
        ha="left",
        va="center",
        fontsize=29,
        fontweight="bold",
        color=accent,
    )
    ax.text(
        x + 0.030,
        y + 0.012,
        sublabel,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=11.9,
        color=COLORS["muted"],
    )


def add_card(
    ax: plt.Axes,
    x: float,
    y: float,
    w: float,
    h: float,
    title: str,
    body: list[tuple[str, str]],
    accent: str,
    fill: str = "#ffffff",
    title_size: int = 19,
    body_size: int = 18,
    value_size: int = 30,
) -> None:
    add_box(ax, x, y, w, h, accent, fill)
    ax.text(
        x + 0.028,
        y + h - 0.040,
        title,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=title_size,
        fontweight="bold",
        color=COLORS["ink"],
    )
    if len(body) == 1:
        label, value = body[0]
        ax.text(
            x + 0.030,
            y + 0.070,
            value,
            transform=ax.transAxes,
            ha="left",
            va="bottom",
            fontsize=value_size,
            fontweight="bold",
            color=accent,
        )
        ax.text(
            x + 0.030,
            y + 0.045,
            label,
            transform=ax.transAxes,
            ha="left",
            va="bottom",
            fontsize=body_size,
            color=COLORS["muted"],
        )
        return

    top = y + h - 0.110
    row_step = min(0.064, (h - 0.150) / max(len(body), 1))
    for i, (label, value) in enumerate(body):
        yy = top - i * row_step
        ax.text(
            x + 0.030,
            yy,
            label,
            transform=ax.transAxes,
            ha="left",
            va="center",
            fontsize=body_size,
            color=COLORS["ink"],
        )
        ax.text(
            x + w - 0.026,
            yy,
            value,
            transform=ax.transAxes,
            ha="right",
            va="center",
            fontsize=body_size + 1,
            fontweight="bold",
            color=accent,
        )


def style_plot(ax: plt.Axes) -> None:
    ax.tick_params(axis="both", labelsize=17, width=1.0, length=6, color="#334155")
    for spine in ax.spines.values():
        spine.set_color("#334155")
        spine.set_linewidth(1.0)
    ax.grid(True, axis="y", color=COLORS["grid"], linewidth=0.8, alpha=0.75)
    ax.set_axisbelow(True)


def add_contract_panel(ax: plt.Axes, x: float, y: float, w: float, h: float) -> None:
    add_box(ax, x, y, w, h, COLORS["blue"], "#ffffff")
    ax.text(
        x + 0.030,
        y + h - 0.045,
        "Analysis contract used downstream",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=22,
        fontweight="bold",
        color=COLORS["ink"],
    )
    rows = [
        ("photon ET", "15-35 GeV"),
        ("tight ID", "14-feature Au+Au BDT WP80"),
        ("isolation", "sliding R=0.4"),
        ("recoil jets", "Delta phi > 7pi/8"),
    ]
    y0 = y + h - 0.092
    step = 0.030
    for i, (left, right) in enumerate(rows):
        yy = y0 - i * step
        ax.text(
            x + 0.038,
            yy,
            left,
            transform=ax.transAxes,
            ha="left",
            va="center",
            fontsize=14.7,
            color=COLORS["ink"],
        )
        ax.text(
            x + 0.165,
            yy,
            right,
            transform=ax.transAxes,
            ha="left",
            va="center",
            fontsize=14.7,
            fontweight="bold",
            color=COLORS["blue"],
        )


def make_slide() -> dict[str, Any]:
    with uproot.open(DATA_ROOT) as f:
        counters = {
            label: {
                "key": key,
                "sum": hist_sum(f, key),
                "entries": hist_entries(f, key),
            }
            for label, key in COUNTER_KEYS.items()
        }
        vertices = {
            label: {
                "key": key,
                "sum": hist_sum(f, key),
                "entries": hist_entries(f, key),
                "mean_cm": hist_mean(f, key),
            }
            for label, key in VERTEX_KEYS.items()
        }
        z_values, z_edges = hist_values(f, VERTEX_KEYS["all"])
        cent_values, cent_edges = hist_values(f, CENTRALITY_KEY)

    working_cent_labels = ["0-20%", "20-50%", "50-80%"]
    cent_counts = np.array([vertices[label]["sum"] for label in working_cent_labels], dtype=float)
    cent_total = float(np.sum(cent_counts))
    cent_centers = 0.5 * (cent_edges[:-1] + cent_edges[1:])
    centrality_0_80_sum = float(np.sum(cent_values[(cent_centers >= 0) & (cent_centers < 80)]))

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "dejavuserif",
        }
    )
    fig = plt.figure(figsize=slide_figsize(), dpi=SLIDE_DPI)
    fig.patch.set_facecolor("white")
    canvas = fig.add_axes([0, 0, 1, 1])
    canvas.axis("off")

    title = "Au+Au data sample for the default-BDT photon-jet pass"
    canvas.text(
        0.052,
        0.942,
        title,
        ha="left",
        va="top",
        fontsize=31,
        fontweight="bold",
        color=COLORS["ink"],
        transform=canvas.transAxes,
    )

    mbd_count = counters["MBD N/S >= 2, |z_vtx| < 150 cm"]["sum"]
    photon10 = counters["photon-10 + MBD/vertex"]["sum"]
    photon12 = counters["photon-12 + MBD/vertex"]["sum"]

    # Compact integrated scope text.
    stat_y = 0.807
    stat_items = [
        (fmt_big(mbd_count), "MBD+vertex events"),
        (f"{RUN_AUDIT['analyzable_paired_runs']} / {RUN_AUDIT['grl_runs_requested']}", "GRL runs with paired streams"),
        (f"{RUN_AUDIT['matched_segment_pairs']:,}", "matched segment pairs"),
        (fmt_big(photon12), "photon-12 events"),
    ]
    xs = [0.065, 0.300, 0.520, 0.760]
    for x, (value, label) in zip(xs, stat_items):
        canvas.text(
            x,
            stat_y,
            value,
            transform=canvas.transAxes,
            ha="left",
            va="baseline",
            fontsize=25,
            fontweight="bold",
            color=COLORS["ink"],
        )
        canvas.text(
            x,
            stat_y - 0.036,
            label,
            transform=canvas.transAxes,
            ha="left",
            va="baseline",
            fontsize=13.2,
            color=COLORS["muted"],
        )

    canvas.text(
        0.065,
        0.715,
        "Event gate: MBD N/S >= 2 and |z_vtx| < 150 cm. Counters are merged DATA ROOT cnt_* sums.",
        transform=canvas.transAxes,
        ha="left",
        va="baseline",
        fontsize=15.8,
        color=COLORS["ink"],
    )
    canvas.text(
        0.065,
        0.685,
        "Run coverage comes from the THE-69 GRL pairing audit; the centrality plot below is the stored h_centrality distribution.",
        transform=canvas.transAxes,
        ha="left",
        va="baseline",
        fontsize=15.8,
        color=COLORS["ink"],
    )

    # Vertex-z plot.
    ax_z = fig.add_axes([0.070, 0.225, 0.405, 0.405])
    z_centers = 0.5 * (z_edges[:-1] + z_edges[1:])
    dz = np.diff(z_edges)
    density = z_values / np.sum(z_values) / dz
    ax_z.step(z_centers, density, where="mid", color=COLORS["blue"], lw=2.4)
    ax_z.fill_between(z_centers, density, step="mid", color=COLORS["blue"], alpha=0.10)
    ax_z.set_xlim(-60, 60)
    ax_z.set_xlabel("z_vtx [cm]", fontsize=17)
    ax_z.set_ylabel("normalized event density", fontsize=17)
    ax_z.set_title("Au+Au z-vertex distribution", fontsize=21, fontweight="bold", pad=10)
    style_plot(ax_z)
    ax_z.text(
        0.055,
        0.925,
        f"mean = {vertices['all']['mean_cm']:.2f} cm",
        transform=ax_z.transAxes,
        ha="left",
        va="top",
        fontsize=14.5,
        color=COLORS["ink"],
    )

    # Centrality distribution.
    ax_c = fig.add_axes([0.555, 0.225, 0.375, 0.405])
    cent_density = cent_values / 1e6
    ax_c.step(cent_centers, cent_density, where="mid", color=COLORS["ink"], lw=2.2)
    ax_c.fill_between(cent_centers, cent_density, step="mid", color="#94a3b8", alpha=0.28)
    for lo, hi, label in [(0, 20, "0-20%"), (20, 50, "20-50%"), (50, 80, "50-80%")]:
        ax_c.axvspan(lo, hi, color=COLORS["blue"], alpha=0.055)
    ax_c.set_xlim(0, 100)
    ax_c.set_xlabel("centrality percentile", fontsize=17)
    ax_c.set_ylabel("events per bin [millions]", fontsize=17)
    ax_c.set_title("Au+Au centrality distribution", fontsize=21, fontweight="bold", pad=10)
    style_plot(ax_c)
    ymax = max(cent_density) * 1.18 if len(cent_density) else 1
    ax_c.set_ylim(0, ymax)
    for lo, hi, label in [(0, 20, "0-20%"), (20, 50, "20-50%"), (50, 80, "50-80%")]:
        ax_c.text(
            (lo + hi) / 2,
            ymax * 0.92,
            label,
            ha="center",
            va="center",
            fontsize=14,
            fontweight="bold",
            color=COLORS["ink"],
        )

    ax_c.text(
        0.98,
        0.06,
        f"0-80%: {fmt_big(centrality_0_80_sum)} entries",
        transform=ax_c.transAxes,
        ha="right",
        va="bottom",
        fontsize=13.5,
        color=COLORS["ink"],
    )

    manifest = {
        "png": str(PNG_OUT),
        "data_root": str(DATA_ROOT),
        "run_audit": RUN_AUDIT,
        "counter_histograms": counters,
        "vertex_histograms": vertices,
        "vertex_plot": {
            "shown_x_range_cm": [-60, 60],
            "stored_h_vertexZ_edges_cm": [float(z_edges[0]), float(z_edges[-1])],
            "note": "The stored h_vertexZ histogram in this merged ROOT is binned from -10 to 10 cm; the slide axis is widened to -60 to 60 cm per request.",
        },
        "working_centrality_sum": cent_total,
        "centrality_histogram": {
            "key": CENTRALITY_KEY,
            "sum": float(np.sum(cent_values)),
            "sum_0_80": centrality_0_80_sum,
        },
        "selection_note": {
            "event_gate": "MBD N/S >= 2 and |z_vtx| < 150 cm",
            "photon_analysis_window": "15 < E_T^gamma < 35 GeV",
            "tight_id": "centAsFeatBase3x3_pt15to35, WP80 T80(c)=0.53471108+0.0012284143*c",
            "isolation": "isoR40_isSliding",
            "back_to_back": "Delta phi > 7pi/8 for xJgamma slides",
        },
        "caveats": [
            "Run and segment coverage are from the THE-69 production audit, not inferred from ROOT run histograms.",
            "Photon-10 and photon-12 counts are cnt_* event-selection counters.",
            "Trigger max-cluster histograms are not used as event counts because their integrals count trigger-filling entries.",
        ],
    }

    fig.savefig(PNG_OUT, dpi=SLIDE_DPI)
    plt.close(fig)
    MANIFEST_OUT.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    SCRIPT_OUT.write_text(
        "\n".join(
            [
                "# Speaker script",
                "",
                "This opening slide is just the scope anchor for the Au+Au pass.",
                f"We start from {RUN_AUDIT['analyzable_paired_runs']} analyzable paired GRL runs, "
                f"with {RUN_AUDIT['matched_segment_pairs']:,} matched CALOFITTING/ZDC segment pairs.",
                f"After the MBD north/south and |z vertex| < 150 cm event gate, the merged ROOT has "
                f"{fmt_big(mbd_count)} event-level entries.",
                f"The stored centrality histogram has {fmt_big(centrality_0_80_sum)} entries in the 0-80 percent range.",
                f"The photon trigger streams are also large: {fmt_big(photon10)} photon-10 plus MBD/vertex entries "
                f"and {fmt_big(photon12)} photon-12 plus MBD/vertex entries.",
                "I am deliberately quoting the cnt histograms here, not the trigger max-cluster energy histograms, "
                "because those max-cluster histograms count trigger-filling entries and are not clean event counts.",
                "The z-vertex plot is drawn on a -60 to 60 cm axis, but the stored h_vertexZ histogram in this merged output is binned from -10 to 10 cm.",
                "",
            ]
        ),
        encoding="utf-8",
    )
    return manifest


def main() -> None:
    manifest = make_slide()
    print(json.dumps({"png": manifest["png"], "manifest": str(MANIFEST_OUT), "script": str(SCRIPT_OUT)}, indent=2))


if __name__ == "__main__":
    main()
