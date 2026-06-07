#!/usr/bin/env python3
"""Build a single THE-44 slide showing the connected fake-production path."""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.path import Path as MplPath
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, PathPatch, Rectangle


REPO = Path(__file__).resolve().parents[4]
BASE = REPO / "dataOutput/auauTightBDTValidation/THE8_branchA_ladder_scorecache_fullstat_20260527"
THE44_DIR = BASE / "diagnostics/the44_high_bdt_pythia_autopsy_20260605/slideReady/full_autopsy_20260605_2108"
DEFAULT_OUTDIR = THE44_DIR / "particle_intensive/connected_flow_master_slide_v6_20260607"

INK = "#111827"
MUTED = "#5B6572"
BG = "#F6F8FB"
CARD = "#FFFFFF"
GRID = "#D7DFEA"
YELLOW = "#FFF8DC"
YELLOW_EDGE = "#D9A600"
GREEN = "#2F8F52"
RED = "#C43C32"
ORANGE = "#D97706"
BLUE = "#2563A6"
PURPLE = "#7651A6"
SLATE = "#54616F"
TEAL = "#0F766E"
SOFT_BLUEGRAY = "#8EA7BF"
SOFT_GREEN = "#8FC7A6"
SOFT_TEAL = "#65B5AA"
SOFT_SLATE = "#9AA9B8"

SOURCE_COLORS = {"Jet20": ORANGE, "Jet30": BLUE, "Jet40": PURPLE}
TRUTH_COLORS = {"gamma truth": GREEN, "pi0/eta": ORANGE, "other": SLATE}
CONE_COLORS = {"photon-rich": GREEN, "mixed": TEAL, "non-photon-rich": RED}
COMP_COLORS = {
    "neutral_meson": ORANGE,
    "charged_hadron": BLUE,
    "photon": GREEN,
    "other": SLATE,
}

plt.rcParams.update(
    {
        "font.family": ["Times New Roman", "Times", "DejaVu Serif"],
        "axes.edgecolor": "#28313F",
        "axes.labelcolor": INK,
        "xtick.color": INK,
        "ytick.color": INK,
    }
)


@dataclass(frozen=True)
class Node:
    name: str
    x: float
    y0: float
    y1: float
    w: float
    count: int
    color: str

    @property
    def h(self) -> float:
        return self.y1 - self.y0

    @property
    def yc(self) -> float:
        return 0.5 * (self.y0 + self.y1)

    @property
    def left(self) -> float:
        return self.x

    @property
    def right(self) -> float:
        return self.x + self.w


def truth_bucket(name: str) -> str:
    if name in {"pi0", "eta"}:
        return "pi0/eta"
    if name == "gamma":
        return "gamma truth"
    return "other"


def load_data() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    cand = pd.read_csv(THE44_DIR / "the44_candidate_autopsy_rows.csv")
    part = pd.read_csv(THE44_DIR / "the44_particle_rows.csv")
    cand["source"] = cand["sample"].str.replace(" inclusive", "", regex=False).str.replace(" signal", "", regex=False)
    cand["truth_bucket"] = cand["truth_pid_name"].map(truth_bucket)
    cand["cone_bucket"] = np.select(
        [cand["final_state_photon_frac"] >= 0.80, cand["final_state_non_photon_frac"] >= 0.80],
        ["photon-rich", "non-photon-rich"],
        default="mixed",
    )
    inc = cand[(cand["label"].eq("high-BDT background")) & cand["sample"].str.contains("Jet", na=False)].copy()
    true = cand[(cand["label"].eq("middling true photon")) & cand["sample"].str.contains("Photon", na=False)].copy()
    return inc, true, part


def txt(
    ax: plt.Axes,
    x: float,
    y: float,
    text: str,
    *,
    size: float,
    color: str = INK,
    weight: str = "normal",
    ha: str = "left",
    va: str = "top",
    linespacing: float = 1.04,
    zorder: int = 30,
) -> None:
    ax.text(
        x,
        y,
        text,
        transform=ax.transAxes,
        fontsize=size,
        color=color,
        fontweight=weight,
        ha=ha,
        va=va,
        linespacing=linespacing,
        zorder=zorder,
    )


def rounded_box(
    ax: plt.Axes,
    x: float,
    y: float,
    w: float,
    h: float,
    *,
    face: str = CARD,
    edge: str = GRID,
    lw: float = 1.2,
    radius: float = 0.016,
    zorder: int = 5,
) -> None:
    ax.add_patch(
        FancyBboxPatch(
            (x, y),
            w,
            h,
            transform=ax.transAxes,
            boxstyle=f"round,pad=0.009,rounding_size={radius}",
            facecolor=face,
            edgecolor=edge,
            linewidth=lw,
            zorder=zorder,
        )
    )


def stack_nodes(
    order: list[str],
    counts: pd.Series,
    colors: dict[str, str],
    *,
    x: float,
    w: float,
    y0: float,
    y1: float,
    min_h: float = 0.042,
    gap: float = 0.015,
) -> dict[str, Node]:
    total = float(sum(counts.get(name, 0) for name in order))
    available = y1 - y0 - gap * (len(order) - 1)
    raw = {name: available * float(counts.get(name, 0)) / total for name in order}
    small = {name for name, h in raw.items() if 0 < h < min_h}
    fixed = min_h * len(small)
    remaining_raw = sum(raw[name] for name in order if name not in small)
    scale = max(available - fixed, 0.02) / max(remaining_raw, 1e-9)
    heights = {name: (min_h if name in small else raw[name] * scale) for name in order}

    nodes: dict[str, Node] = {}
    cursor = y1
    for name in order:
        h = heights[name]
        y_top = cursor
        y_bot = cursor - h
        nodes[name] = Node(name=name, x=x, y0=y_bot, y1=y_top, w=w, count=int(counts.get(name, 0)), color=colors[name])
        cursor = y_bot - gap
    return nodes


def ribbon(
    ax: plt.Axes,
    x0: float,
    y0a: float,
    y0b: float,
    x1: float,
    y1a: float,
    y1b: float,
    *,
    color: str,
    alpha: float,
    edge: str | None = None,
    lw: float = 0.0,
    zorder: int = 2,
) -> None:
    dx = x1 - x0
    c0 = x0 + 0.46 * dx
    c1 = x1 - 0.46 * dx
    verts = [
        (x0, y0a),
        (c0, y0a),
        (c1, y1a),
        (x1, y1a),
        (x1, y1b),
        (c1, y1b),
        (c0, y0b),
        (x0, y0b),
        (x0, y0a),
    ]
    codes = [
        MplPath.MOVETO,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.LINETO,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.CLOSEPOLY,
    ]
    ax.add_patch(
        PathPatch(
            MplPath(verts, codes),
            transform=ax.transAxes,
            facecolor=color,
            edgecolor=edge or color,
            linewidth=lw,
            alpha=alpha,
            zorder=zorder,
        )
    )


def segments(
    left_nodes: dict[str, Node],
    right_nodes: dict[str, Node],
    flows: list[tuple[str, str, int]],
    *,
    left_counts: pd.Series,
    right_counts: pd.Series,
) -> list[tuple[str, str, int, tuple[float, float], tuple[float, float]]]:
    left_cursor = {name: left_nodes[name].y0 for name in left_nodes}
    right_cursor = {name: right_nodes[name].y0 for name in right_nodes}
    out = []
    for a, b, count in flows:
        if count <= 0:
            continue
        left_node = left_nodes[a]
        right_node = right_nodes[b]
        left_h = left_node.h * count / max(float(left_counts.get(a, 0)), 1.0)
        right_h = right_node.h * count / max(float(right_counts.get(b, 0)), 1.0)
        l0 = left_cursor[a]
        l1 = l0 + left_h
        r0 = right_cursor[b]
        r1 = r0 + right_h
        left_cursor[a] = l1
        right_cursor[b] = r1
        out.append((a, b, count, (l0, l1), (r0, r1)))
    return out


def draw_node(ax: plt.Axes, node: Node, label: str, *, highlight: bool = False) -> None:
    edge = YELLOW_EDGE if highlight else "#FFFFFF"
    lw = 2.7 if highlight else 1.1
    rounded_box(ax, node.x, node.y0, node.w, node.h, face=node.color, edge=edge, lw=lw, radius=0.014, zorder=16)
    txt(
        ax,
        node.x + node.w / 2,
        node.yc,
        label,
        size=14.2 if node.h > 0.08 else 12.0,
        color="#FFFFFF",
        weight="bold",
        ha="center",
        va="center",
        linespacing=0.92,
        zorder=22,
    )


def read_key(ax: plt.Axes) -> None:
    rounded_box(ax, 0.055, 0.804, 0.890, 0.060, face="#FFFFFF", edge="#D8E0EA", lw=1.0, radius=0.012, zorder=24)
    txt(ax, 0.078, 0.845, "How to read", size=13.0, color=INK, weight="bold", va="center", zorder=31)
    x0 = 0.205
    ax.add_patch(Rectangle((x0, 0.830), 0.023, 0.022, transform=ax.transAxes, facecolor=RED, edgecolor=YELLOW_EDGE, linewidth=1.5, zorder=31))
    txt(ax, x0 + 0.033, 0.845, "block size = candidate count", size=11.8, color=INK, va="center", zorder=31)
    ribbon(ax, 0.405, 0.831, 0.852, 0.475, 0.831, 0.852, color=RED, alpha=0.42, edge=YELLOW_EDGE, lw=0.65, zorder=31)
    txt(ax, 0.492, 0.845, "ribbon width = flow count", size=11.8, color=INK, va="center", zorder=31)
    ax.add_patch(Rectangle((0.672, 0.830), 0.023, 0.022, transform=ax.transAxes, facecolor="#FFFFFF", edgecolor=YELLOW_EDGE, linewidth=1.7, zorder=31))
    txt(ax, 0.704, 0.845, "gold outline = dominant branch", size=11.8, color=INK, weight="bold", va="center", zorder=31)
    txt(ax, 0.704, 0.823, "node labels show n/160 and median BDT", size=10.4, color=MUTED, va="center", zorder=31)


def header_pill(ax: plt.Axes, x: float, label: str, w: float = 0.165) -> None:
    rounded_box(ax, x - w / 2, 0.750, w, 0.040, face=BG, edge=BG, lw=0.0, radius=0.010, zorder=20)
    txt(ax, x, 0.771, label, size=15.5, weight="bold", ha="center", va="center", color=INK, zorder=32)


def stage_arrow(ax: plt.Axes, x0: float, x1: float) -> None:
    ax.add_patch(
        FancyArrowPatch(
            (x0, 0.771),
            (x1, 0.771),
            transform=ax.transAxes,
            arrowstyle="-|>",
            mutation_scale=15,
            linewidth=1.8,
            color=SLATE,
            zorder=32,
        )
    )


def fmt_pct(n: int, total: int) -> str:
    return f"{100 * n / total:.0f}%"


def make_slide(inc: pd.DataFrame, true: pd.DataFrame, part: pd.DataFrame, outdir: Path) -> tuple[Path, dict[str, float]]:
    fig = plt.figure(figsize=(16, 9), dpi=160)
    fig.patch.set_facecolor(BG)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    total = len(inc)
    source_order = ["Jet40", "Jet30", "Jet20"]
    truth_order = ["gamma truth", "pi0/eta", "other"]
    cone_order = ["photon-rich", "mixed", "non-photon-rich"]

    source_counts = inc.groupby("source").size().reindex(source_order).fillna(0).astype(int)
    truth_counts = inc.groupby("truth_bucket").size().reindex(truth_order).fillna(0).astype(int)
    cone_counts = inc.groupby("cone_bucket").size().reindex(cone_order).fillna(0).astype(int)

    source_nodes = stack_nodes(source_order, source_counts, SOURCE_COLORS, x=0.060, w=0.130, y0=0.290, y1=0.735)
    truth_nodes = stack_nodes(truth_order, truth_counts, TRUTH_COLORS, x=0.305, w=0.135, y0=0.290, y1=0.735)
    cone_nodes = stack_nodes(cone_order, cone_counts, CONE_COLORS, x=0.550, w=0.155, y0=0.290, y1=0.735)
    endpoint_nodes = {"BDT tail": Node("BDT tail", x=0.815, y0=0.300, y1=0.725, w=0.130, count=total, color=RED)}

    txt(ax, 0.055, 0.945, "Jet30/40 reveal the high-BDT fake tail is a fragmentation path", size=28.0, weight="bold")
    txt(
        ax,
        0.055,
        0.902,
        "Ribbons propagate 160 high-score inclusive backgrounds candidate-by-candidate from source bin to truth seed to local cone to BDT endpoint.",
        size=15.8,
        color=MUTED,
    )
    read_key(ax)

    header_pill(ax, 0.125, "1. Pythia source bin", w=0.180)
    header_pill(ax, 0.372, "2. Truth seed", w=0.145)
    header_pill(ax, 0.630, "3. Local cone identity", w=0.190)
    header_pill(ax, 0.880, "4. BDT endpoint", w=0.165)
    stage_arrow(ax, 0.206, 0.258)
    stage_arrow(ax, 0.462, 0.514)
    stage_arrow(ax, 0.730, 0.782)

    rounded_box(ax, 0.043, 0.247, 0.912, 0.525, face="#FBFCFE", edge="#D6DEE9", lw=1.2, radius=0.016, zorder=0)

    source_truth = pd.crosstab(inc["source"], inc["truth_bucket"])
    truth_cone = pd.crosstab(inc["truth_bucket"], inc["cone_bucket"])
    cone_endpoint = pd.Series({"BDT tail": total})

    stage1_flows = [(s, t, int(source_truth.reindex(index=source_order, columns=truth_order).fillna(0).loc[s, t])) for s in source_order for t in truth_order]
    for a, b, count, left_seg, right_seg in segments(source_nodes, truth_nodes, stage1_flows, left_counts=source_counts, right_counts=truth_counts):
        highlight = a in {"Jet30", "Jet40"} and b == "pi0/eta"
        ribbon(
            ax,
            source_nodes[a].right,
            left_seg[0],
            left_seg[1],
            truth_nodes[b].left,
            right_seg[0],
            right_seg[1],
            color=SOURCE_COLORS[a] if highlight else SOFT_BLUEGRAY,
            alpha=0.68 if highlight else 0.27,
            edge=YELLOW_EDGE if highlight else None,
            lw=0.80 if highlight else 0.0,
            zorder=6 if highlight else 2,
        )

    stage2_flows = [(t, c, int(truth_cone.reindex(index=truth_order, columns=cone_order).fillna(0).loc[t, c])) for t in truth_order for c in cone_order]
    for a, b, count, left_seg, right_seg in segments(truth_nodes, cone_nodes, stage2_flows, left_counts=truth_counts, right_counts=cone_counts):
        highlight = a == "pi0/eta" and b == "non-photon-rich"
        ribbon(
            ax,
            truth_nodes[a].right,
            left_seg[0],
            left_seg[1],
            cone_nodes[b].left,
            right_seg[0],
            right_seg[1],
            color=RED if highlight else {"gamma truth": SOFT_GREEN, "pi0/eta": ORANGE, "other": SOFT_SLATE}[a],
            alpha=0.66 if highlight else 0.27,
            edge=YELLOW_EDGE if highlight else None,
            lw=0.80 if highlight else 0.0,
            zorder=7 if highlight else 3,
        )

    stage3_flows = [(c, "BDT tail", int(cone_counts.loc[c])) for c in cone_order]
    for a, b, count, left_seg, right_seg in segments(cone_nodes, endpoint_nodes, stage3_flows, left_counts=cone_counts, right_counts=cone_endpoint):
        highlight = a == "non-photon-rich"
        ribbon(
            ax,
            cone_nodes[a].right,
            left_seg[0],
            left_seg[1],
            endpoint_nodes[b].left,
            right_seg[0],
            right_seg[1],
            color=RED if highlight else {"photon-rich": SOFT_GREEN, "mixed": SOFT_TEAL, "non-photon-rich": RED}[a],
            alpha=0.62 if highlight else 0.24,
            edge=YELLOW_EDGE if highlight else None,
            lw=0.80 if highlight else 0.0,
            zorder=8 if highlight else 4,
        )

    bdt_by_source = inc.groupby("source")["auau_tight_bdt_score"].median()
    bdt_by_truth = inc.groupby("truth_bucket")["auau_tight_bdt_score"].median()
    bdt_by_cone = inc.groupby("cone_bucket")["auau_tight_bdt_score"].median()

    for name in source_order:
        n = int(source_counts.loc[name])
        draw_node(ax, source_nodes[name], f"{name}\n{n}/{total}\nBDT {bdt_by_source.loc[name]:.2f}", highlight=name in {"Jet30", "Jet40"})
    for name in truth_order:
        n = int(truth_counts.loc[name])
        draw_node(ax, truth_nodes[name], f"{name}\n{n}/{total}\nBDT {bdt_by_truth.loc[name]:.2f}", highlight=name == "pi0/eta")
    for name in cone_order:
        n = int(cone_counts.loc[name])
        label = f"{name.replace('-', ' ')}\n{n}/{total}\nBDT {bdt_by_cone.loc[name]:.2f}"
        draw_node(ax, cone_nodes[name], label, highlight=name == "non-photon-rich")

    endpoint = endpoint_nodes["BDT tail"]
    draw_node(ax, endpoint, f"high-BDT fake tail\nBDT >= 0.80\nn={total}\nmedian {inc['auau_tight_bdt_score'].median():.2f}", highlight=True)

    dom = inc[
        inc["source"].isin(["Jet30", "Jet40"])
        & inc["truth_bucket"].eq("pi0/eta")
        & inc["cone_bucket"].eq("non-photon-rich")
    ].copy()
    dom_ids = set(dom["candidate_uid"])
    pdom = part[part["candidate_uid"].isin(dom_ids) & part["status"].eq(1)]
    budget = pdom.groupby("group")["pt"].sum()
    budget = budget.reindex(["neutral_meson", "charged_hadron", "photon", "other"]).fillna(0.0)
    budget_frac = budget / max(float(budget.sum()), 1.0)

    rounded_box(ax, 0.055, 0.048, 0.890, 0.174, face=YELLOW, edge=YELLOW_EDGE, lw=1.35, radius=0.018, zorder=20)
    txt(
        ax,
        0.078,
        0.204,
        "Answer from the autopsy: high-BDT fakes are mostly hard neutral-meson fragmentation",
        size=15.2,
        weight="bold",
        color=INK,
        zorder=30,
    )
    rounded_box(ax, 0.073, 0.065, 0.285, 0.101, face="#FFFFFF", edge="#E7C75B", lw=0.9, radius=0.012, zorder=24)
    rounded_box(ax, 0.376, 0.065, 0.302, 0.101, face="#FFFFFF", edge="#E7C75B", lw=0.9, radius=0.012, zorder=24)
    rounded_box(ax, 0.696, 0.065, 0.225, 0.101, face="#FFFFFF", edge="#E7C75B", lw=0.9, radius=0.012, zorder=24)

    txt(ax, 0.091, 0.153, "Dominant answer", size=12.4, weight="bold", color=INK, zorder=30)
    txt(ax, 0.091, 0.128, f"{len(dom)}/{total}", size=22.5, weight="bold", color=RED, zorder=30)
    txt(ax, 0.196, 0.122, f"median BDT = {dom['auau_tight_bdt_score'].median():.2f}", size=13.2, weight="bold", color=INK, zorder=30)
    txt(
        ax,
        0.091,
        0.092,
        "Jet30/40 -> $\\pi^0/\\eta$ seed -> non-photon-rich cone -> high BDT",
        size=10.4,
        color=INK,
        zorder=30,
    )

    bar_x, bar_y, bar_w, bar_h = 0.398, 0.107, 0.250, 0.031
    txt(ax, 0.395, 0.154, "Inside the dominant branch's cone", size=12.4, weight="bold", color=INK, zorder=36)
    cursor = bar_x
    for group, frac in budget_frac.items():
        width = bar_w * float(frac)
        if width <= 0:
            continue
        ax.add_patch(
            Rectangle(
                (cursor, bar_y),
                width,
                bar_h,
                transform=ax.transAxes,
                facecolor=COMP_COLORS[group],
                edgecolor="white",
                linewidth=1.1,
                zorder=31,
            )
        )
        cursor += width
    labels = [
        (ORANGE, f"{100 * budget_frac['neutral_meson']:.0f}% neutral meson"),
        (BLUE, f"{100 * budget_frac['charged_hadron']:.0f}% charged"),
        (GREEN, f"~{100 * budget_frac['photon']:.0f}% photon"),
    ]
    for i, (color, label) in enumerate(labels):
        x = bar_x + 0.002 + 0.083 * i
        ax.scatter([x], [0.091], transform=ax.transAxes, s=34, marker="s", color=color, edgecolors="none", zorder=35)
        txt(ax, x + 0.011, 0.100, label, size=9.4, color=INK, va="center", zorder=35)
    txt(
        ax,
        0.395,
        0.078,
        "High BDT does not imply a photon-rich truth cone.",
        size=10.0,
        color=INK,
        zorder=30,
    )

    txt(ax, 0.714, 0.153, "Control check", size=12.4, weight="bold", color=GREEN, zorder=30)
    txt(
        ax,
        0.714,
        0.126,
        f"Photon-sample control:\nphoton-rich, median BDT {true['auau_tight_bdt_score'].median():.2f}",
        size=10.0,
        color=INK,
        linespacing=0.92,
        zorder=30,
    )
    txt(
        ax,
        0.714,
        0.087,
        "This tail is fragmentation that\nlooks photon-like to the BDT.",
        size=9.8,
        color=INK,
        weight="bold",
        linespacing=0.92,
        zorder=30,
    )

    out = outdir / "the44_connected_flow_master_slide_v6.png"
    fig.savefig(out, dpi=160)
    plt.close(fig)

    metrics = {
        "inclusive_high_bdt_backgrounds": float(total),
        "truth_photon_controls": float(len(true)),
        "jet30_40_candidates": float(inc["source"].isin(["Jet30", "Jet40"]).sum()),
        "pi0_eta_truth_candidates": float(inc["truth_bucket"].eq("pi0/eta").sum()),
        "non_photon_rich_cones": float(inc["cone_bucket"].eq("non-photon-rich").sum()),
        "dominant_path_candidates": float(len(dom)),
        "dominant_path_median_bdt": float(dom["auau_tight_bdt_score"].median()),
        "dominant_path_neutral_meson_pt_fraction": float(budget_frac["neutral_meson"]),
        "dominant_path_charged_hadron_pt_fraction": float(budget_frac["charged_hadron"]),
        "dominant_path_photon_pt_fraction": float(budget_frac["photon"]),
        "median_high_fake_bdt_score": float(inc["auau_tight_bdt_score"].median()),
        "median_high_fake_photon_fraction": float(inc["final_state_photon_frac"].median()),
        "median_truth_photon_bdt_score": float(true["auau_tight_bdt_score"].median()),
        "median_truth_photon_photon_fraction": float(true["final_state_photon_frac"].median()),
    }
    return out, metrics


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    inc, true, part = load_data()
    slide, metrics = make_slide(inc, true, part, args.outdir)

    script_path = args.outdir / "the44_connected_flow_master_slide_v6_script.md"
    script_path.write_text(
        "# THE-44 Connected Flow Master Slide Script\n\n"
        "This slide follows the high-BDT inclusive background candidates as a connected particle-production path, instead of asking the audience to decode isolated markers. "
        "The small guide at the top says how to read it: block height is the number of candidates at that stage, ribbon width is how many candidates flow between two categories, and the BDT number printed in each block is the median score for that block. "
        "Reading left to right, most of the high-score fake tail comes from Jet30 and Jet40 parent samples, then collapses into pi0/eta truth seeds, then lands in non-photon-rich status-1 local cones. "
        "The highlighted branch is the core result: Jet30/40 to pi0/eta to a non-photon cone accounts for 115 of the 160 high-BDT backgrounds, and that branch still has a median BDT score of about 0.84.\n\n"
        "The bottom verdict strip states why that path matters physically. The dominant branch is 115 of 160 candidates with median BDT about 0.84. Inside that branch's local cone, the status-1 pT budget is about 71 percent neutral meson and 29 percent charged hadron, with essentially no photon component. "
        "The photon-sample control is photon-rich but has a lower median BDT, so this tail is not clean photons; it is fragmentation that looks photon-like to the shower-shape BDT.\n"
    )
    manifest = {
        "status": "READY",
        "slide": str(slide),
        "script": str(script_path),
        "metrics": metrics,
        "inputs": {
            "candidate_rows": str(THE44_DIR / "the44_candidate_autopsy_rows.csv"),
            "particle_rows": str(THE44_DIR / "the44_particle_rows.csv"),
        },
        "caveat": "Internal qualitative diagnostic from saved THE-44 candidate and particle rows; not a production-weighted purity estimate.",
    }
    manifest_path = args.outdir / "the44_connected_flow_master_slide_v6_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
