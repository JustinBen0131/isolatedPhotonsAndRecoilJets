#!/usr/bin/env python3
"""Make an audience-facing flow-chart slide for the total-calo threshold."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


OUTDIR = Path("dataOutput/auauTightBDTValidation/THE32_lowCaloDiagnosticClosure_20260603")
PNG = OUTDIR / "event_energy_veto_threshold_flow_slide_v8_20260606.png"
PHONE_PNG = OUTDIR / "event_energy_veto_threshold_flow_slide_v8_phone_refresh_20260606.png"
SCRIPT = OUTDIR / "event_energy_veto_threshold_flow_slide_script_20260606.md"
MANIFEST = OUTDIR / "event_energy_veto_threshold_flow_slide_manifest_20260606.json"

INK = "#172033"
MUTED = "#5D6B7A"
LIGHT_TEXT = "#F8FAFC"
BLUE = "#1F77B4"
BLUE_SOFT = "#EAF4FB"
RED = "#D62728"
RED_SOFT = "#FFF1F2"
GOLD = "#B7791F"
GOLD_SOFT = "#FFF7ED"
GREEN = "#1F7A4D"
GREEN_SOFT = "#ECFDF3"
GRAY_SOFT = "#F8FAFC"
GRAY_EDGE = "#CBD5E1"

plt.rcParams.update(
    {
        "font.family": "Times New Roman",
        "mathtext.fontset": "stix",
        "axes.unicode_minus": False,
    }
)


def rounded_box(
    ax,
    xy: tuple[float, float],
    wh: tuple[float, float],
    *,
    face: str = GRAY_SOFT,
    edge: str = GRAY_EDGE,
    lw: float = 1.4,
    radius: float = 0.018,
    zorder: int = 1,
) -> FancyBboxPatch:
    x, y = xy
    w, h = wh
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle=f"round,pad=0.012,rounding_size={radius}",
        facecolor=face,
        edgecolor=edge,
        linewidth=lw,
        transform=ax.transAxes,
        zorder=zorder,
    )
    ax.add_patch(patch)
    return patch


def text_box(
    ax,
    xy: tuple[float, float],
    wh: tuple[float, float],
    label: str,
    body_lines: list[str],
    *,
    face: str = GRAY_SOFT,
    edge: str = GRAY_EDGE,
    label_color: str = INK,
    body_color: str = MUTED,
    label_size: float = 16.0,
    body_size: float = 12.8,
    align: str = "left",
) -> None:
    x, y = xy
    w, h = wh
    rounded_box(ax, xy, wh, face=face, edge=edge)
    ha = "center" if align == "center" else "left"
    tx = x + w / 2 if align == "center" else x + 0.022
    ax.text(
        tx,
        y + h - 0.035,
        label,
        transform=ax.transAxes,
        fontsize=label_size,
        color=label_color,
        fontweight="bold",
        ha=ha,
        va="top",
        zorder=5,
    )
    ax.text(
        tx,
        y + h - 0.092,
        "\n".join(body_lines),
        transform=ax.transAxes,
        fontsize=body_size,
        color=body_color,
        ha=ha,
        va="top",
        linespacing=1.18,
        zorder=5,
    )


def arrow(
    ax,
    start: tuple[float, float],
    end: tuple[float, float],
    *,
    color: str = "#94A3B8",
    curve: float = 0.0,
    scale: float = 18.0,
) -> None:
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            transform=ax.transAxes,
            arrowstyle="-|>",
            connectionstyle=f"arc3,rad={curve}",
            mutation_scale=scale,
            linewidth=2.0,
            color=color,
            shrinkA=5,
            shrinkB=5,
            zorder=3,
        )
    )


def render() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(12.8, 7.2), dpi=200)
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    fig.patch.set_facecolor("white")
    ax.set_axis_off()

    ax.text(
        0.055,
        0.940,
        "How the event-energy veto threshold is derived",
        transform=ax.transAxes,
        fontsize=31,
        color=INK,
        fontweight="bold",
        va="top",
    )

    rounded_box(ax, (0.070, 0.805), (0.860, 0.065), face="#F1F5F9", edge="#D9E2EC")
    ax.text(
        0.500,
        0.837,
        r"One threshold is set in each 5% centrality bin using total event calorimeter energy only.",
        transform=ax.transAxes,
        fontsize=16.4,
        color=INK,
        fontweight="bold",
        ha="center",
        va="center",
    )

    rounded_box(ax, (0.095, 0.585), (0.810, 0.195), face=GOLD_SOFT, edge="#F1C27D", lw=2.0)
    ax.text(
        0.500,
        0.735,
        "Applied event-energy veto",
        transform=ax.transAxes,
        fontsize=18.5,
        color=GOLD,
        fontweight="bold",
        ha="center",
        va="center",
    )
    ax.text(
        0.500,
        0.672,
        r"$T_{\rm bin}=\max\left[\mathrm{median}(y)-5\cdot1.4826\,\mathrm{MAD}(y),\ q_{0.1\%}(y)\right]$",
        transform=ax.transAxes,
        fontsize=19.0,
        color=INK,
        fontweight="bold",
        ha="center",
        va="center",
    )
    ax.text(
        0.500,
        0.613,
        r"remove the event if $y<T_{\rm bin}$",
        transform=ax.transAxes,
        fontsize=16.4,
        color=RED,
        fontweight="bold",
        ha="center",
        va="center",
    )

    ax.text(
        0.070,
        0.525,
        "Derivation in each 5% centrality bin",
        transform=ax.transAxes,
        fontsize=19.5,
        color=INK,
        fontweight="bold",
        va="center",
    )
    ax.plot([0.070, 0.930], [0.500, 0.500], transform=ax.transAxes, color="#E2E8F0", lw=1.4)

    text_box(
        ax,
        (0.060, 0.292),
        (0.205, 0.186),
        "1. Input",
        [
            r"$y=\log_{10}(E_{\rm CEMC}+E_{\rm IHCal}$",
            r"$+E_{\rm OHCal}+1)$",
            "in one 5% bin",
        ],
        face=BLUE_SOFT,
        edge="#BBD7EA",
        label_color=BLUE,
        body_color=INK,
        label_size=16.0,
        body_size=12.0,
        align="center",
    )
    text_box(
        ax,
        (0.295, 0.292),
        (0.205, 0.186),
        "2. Center",
        [
            r"$\mathrm{median}(y)$",
            "marks the middle",
            "of the normal band",
        ],
        face=GRAY_SOFT,
        edge=GRAY_EDGE,
        label_color=BLUE,
        body_color=INK,
        label_size=16.0,
        body_size=12.9,
        align="center",
    )
    text_box(
        ax,
        (0.530, 0.292),
        (0.205, 0.186),
        "3. Width",
        [
            r"$\mathrm{MAD}$ measures spread",
            r"convert MAD into",
            r"a $\sigma$-like width",
            r"using the Gaussian rule below",
        ],
        face=GRAY_SOFT,
        edge=GRAY_EDGE,
        label_color=BLUE,
        body_color=INK,
        label_size=16.0,
        body_size=11.5,
        align="center",
    )
    text_box(
        ax,
        (0.765, 0.292),
        (0.205, 0.186),
        "4. Threshold",
        [
            "place the line",
            r"$5$ widths below center",
            "as a conservative edge",
            r"then apply $q_{0.1\%}$ floor",
        ],
        face=GRAY_SOFT,
        edge=GRAY_EDGE,
        label_color=BLUE,
        body_color=INK,
        label_size=16.0,
        body_size=11.6,
        align="center",
    )

    arrow(ax, (0.265, 0.385), (0.295, 0.385), color="#94A3B8")
    arrow(ax, (0.500, 0.385), (0.530, 0.385), color="#94A3B8")
    arrow(ax, (0.735, 0.385), (0.765, 0.385), color="#94A3B8")
    rounded_box(ax, (0.090, 0.122), (0.430, 0.130), face=GOLD_SOFT, edge="#F1C27D", lw=1.6)
    ax.text(
        0.125,
        0.214,
        "Where 1.4826 comes from",
        transform=ax.transAxes,
        fontsize=15.8,
        color=GOLD,
        fontweight="bold",
        va="center",
    )
    ax.text(
        0.125,
        0.172,
        r"For a normal band, $P(|Z|<0.6745)=0.5$.",
        transform=ax.transAxes,
        fontsize=12.8,
        color=INK,
        va="center",
    )
    ax.text(
        0.125,
        0.140,
        r"So $\mathrm{MAD}=0.6745\sigma$, and $\sigma=\mathrm{MAD}/0.6745=1.4826\,\mathrm{MAD}$.",
        transform=ax.transAxes,
        fontsize=12.3,
        color=INK,
        va="center",
    )

    rounded_box(ax, (0.550, 0.122), (0.360, 0.130), face=GREEN_SOFT, edge="#A7F3D0", lw=1.6)
    ax.text(
        0.585,
        0.214,
        "Why the max with the 0.1% quantile?",
        transform=ax.transAxes,
        fontsize=14.9,
        color=GREEN,
        fontweight="bold",
        va="center",
    )
    ax.text(
        0.585,
        0.170,
        r"$q_{0.1\%}(y)$ is a lower-percentile floor.",
        transform=ax.transAxes,
        fontsize=12.8,
        color=INK,
        va="center",
    )
    ax.text(
        0.585,
        0.140,
        r"It prevents rare extremes from pulling $T_{\rm bin}$ too low.",
        transform=ax.transAxes,
        fontsize=12.8,
        color=INK,
        va="center",
    )

    ax.text(
        0.500,
        0.068,
        "The derivation is blind to sample type, truth label, BDT score, and merge weights.",
        transform=ax.transAxes,
        fontsize=13.8,
        color=MUTED,
        ha="center",
        va="bottom",
    )

    fig.savefig(PNG, dpi=200)
    fig.savefig(PHONE_PNG, dpi=200)
    plt.close(fig)

    SCRIPT.write_text(
        "\n".join(
            [
                "# Event-Energy Veto Threshold Flow Slide Script",
                "",
                "This slide explains how the threshold was set before the BDT training.",
                "In each 5% centrality bin, I compute the total event calorimeter energy variable y.",
                "The median gives the center of the normal event-energy band.",
                "The median absolute deviation measures the width of that band without being pulled around by the low-energy tail.",
                "The factor 1.4826 converts the median absolute deviation into a sigma-like width for a Gaussian-shaped band.",
                "The threshold is then placed five of those widths below the median, which is intentionally conservative and keeps the normal band above the line.",
                "The lower 0.1% quantile is included as a floor so a few extreme events cannot push the threshold too low.",
                "The final veto is simple: remove events with total event-calo energy below T bin.",
                "",
            ]
        )
    )
    MANIFEST.write_text(
        json.dumps(
            {
                "png": str(PNG),
                "phone_refresh_png": str(PHONE_PNG),
                "script": str(SCRIPT),
                "purpose": "Audience-facing flow chart explaining the total-calo event-veto threshold derivation.",
                "formula": "T_bin = max[median(y) - 5*1.4826*MAD(y), q_0.1%(y)]",
                "cut_variable": "y = log10(E_CEMC + E_IHCal + E_OHCal + 1)",
                "visual_scope": "Explains the threshold derivation only; no Google Slides mutation.",
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )


if __name__ == "__main__":
    render()
