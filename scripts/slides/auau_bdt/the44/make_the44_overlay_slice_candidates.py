#!/usr/bin/env python3
"""Build overlay/slice THE-44 slide candidates from saved BDT/Pythia tables."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Polygon, Rectangle


REPO = Path(__file__).resolve().parents[4]
BASE = REPO / "dataOutput/auauTightBDTValidation/THE8_branchA_ladder_scorecache_fullstat_20260527"
THE44_DIR = BASE / "diagnostics/the44_high_bdt_pythia_autopsy_20260605/slideReady/full_autopsy_20260605_2108"
DEFAULT_OUTDIR = THE44_DIR / "particle_intensive/deep_overlay_slice_candidates_20260606"

INK = "#111827"
MUTED = "#5D6675"
LIGHT = "#F6F7F9"
CARD = "#FFFFFF"
GRID = "#DCE3EC"
RED = "#C43C32"
GREEN = "#2F8F52"
ORANGE = "#D97706"
BLUE = "#2563A6"
PURPLE = "#7651A6"
TEAL = "#0F766E"
SLATE = "#475569"
YELLOW = "#FFF4C7"

SOURCE_COLORS = {"Jet20": "#F97316", "Jet30": GREEN, "Jet40": PURPLE}
SOURCE_MARKERS = {"Jet20": "o", "Jet30": "s", "Jet40": "^"}
TRUTH_COLORS = {"pi0/eta": ORANGE, "gamma truth": RED, "other": SLATE}
TRUTH_MARKERS = {"pi0/eta": "o", "gamma truth": "^", "other": "s"}
GROUP_COLORS = {
    "photon": RED,
    "neutral_meson": ORANGE,
    "charged_hadron": BLUE,
    "quark_gluon": PURPLE,
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


def canvas() -> tuple[plt.Figure, plt.Axes]:
    fig = plt.figure(figsize=(16, 9), dpi=160)
    fig.patch.set_facecolor(LIGHT)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    return fig, ax


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
    linespacing: float = 1.10,
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
        zorder=20,
    )


def card(
    ax: plt.Axes,
    x: float,
    y: float,
    w: float,
    h: float,
    *,
    face: str = CARD,
    edge: str = "#D4DCE7",
    lw: float = 1.2,
    radius: float = 0.018,
    zorder: int = 1,
) -> None:
    ax.add_patch(
        FancyBboxPatch(
            (x, y),
            w,
            h,
            transform=ax.transAxes,
            boxstyle=f"round,pad=0.010,rounding_size={radius}",
            facecolor=face,
            edgecolor=edge,
            linewidth=lw,
            zorder=zorder,
        )
    )


def title(ax: plt.Axes, claim: str, subtitle: str) -> None:
    txt(ax, 0.055, 0.935, claim, size=30, weight="bold")
    txt(ax, 0.055, 0.875, subtitle, size=18, color=MUTED)


def inset(fig: plt.Figure, x: float, y: float, w: float, h: float) -> plt.Axes:
    return fig.add_axes([x, y, w, h])


def finish(fig: plt.Figure, out: Path) -> Path:
    fig.savefig(out, dpi=160)
    plt.close(fig)
    return out


def truth_bucket(name: str) -> str:
    if name in {"pi0", "eta"}:
        return "pi0/eta"
    if name == "gamma":
        return "gamma truth"
    return "other"


def load_data() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    cand = pd.read_csv(THE44_DIR / "the44_candidate_autopsy_rows.csv")
    part = pd.read_csv(THE44_DIR / "the44_particle_rows.csv")
    examples = pd.read_csv(THE44_DIR / "the44_event_card_examples.csv")
    cand["source"] = cand["sample"].str.replace(" inclusive", "", regex=False).str.replace(" signal", "", regex=False)
    cand["truth_bucket"] = cand["truth_pid_name"].map(truth_bucket)
    cand["cone_bucket"] = np.select(
        [cand["final_state_photon_frac"] >= 0.80, cand["final_state_non_photon_frac"] >= 0.80],
        ["photon-rich", "non-photon-rich"],
        default="mixed",
    )
    inc = cand[(cand["label"].eq("high-BDT background")) & cand["sample"].str.contains("Jet", na=False)].copy()
    true = cand[(cand["label"].eq("middling true photon")) & cand["sample"].str.contains("Photon", na=False)].copy()
    return cand, inc, true, part.merge(cand[["candidate_uid", "label", "truth_bucket", "source"]], on="candidate_uid", how="left"), examples


def style_axis(ax: plt.Axes) -> None:
    ax.set_facecolor(CARD)
    ax.grid(True, color=GRID, linewidth=0.8, alpha=0.85)
    for spine in ax.spines.values():
        spine.set_color("#C8D2DF")
        spine.set_linewidth(1.0)


def legend_box(ax: plt.Axes, x: float, y: float, lines: list[tuple[str, str, str]]) -> None:
    card(ax, x, y - 0.035 * len(lines) - 0.010, 0.19, 0.035 * len(lines) + 0.030)
    for i, (label, color, marker) in enumerate(lines):
        yy = y - 0.035 * i
        ax.scatter([x + 0.022], [yy - 0.004], transform=ax.transAxes, s=70, marker=marker, color=color, edgecolors="white", linewidths=0.8, zorder=30)
        txt(ax, x + 0.045, yy + 0.005, label, size=12.5, va="center")


def scatter_by_source_truth(axp: plt.Axes, df: pd.DataFrame, x: str, y: str, *, alpha: float = 0.82, size_col: str | None = None) -> None:
    for source in ["Jet20", "Jet30", "Jet40"]:
        for truth in ["pi0/eta", "gamma truth", "other"]:
            sub = df[(df["source"].eq(source)) & (df["truth_bucket"].eq(truth))]
            if sub.empty:
                continue
            sizes = np.full(len(sub), 72.0)
            if size_col:
                sizes = np.clip(28 + 4.0 * sub[size_col].to_numpy(), 40, 150)
            axp.scatter(
                sub[x],
                sub[y],
                s=sizes,
                marker=SOURCE_MARKERS[source],
                facecolors=TRUTH_COLORS[truth],
                edgecolors=SOURCE_COLORS[source],
                linewidths=1.2,
                alpha=alpha,
                label=f"{source} {truth}",
            )


def slide_01_score_particle_overlay(inc: pd.DataFrame, true: pd.DataFrame, outdir: Path) -> Path:
    fig, ax = canvas()
    title(
        ax,
        "One overlay exposes the BDT's particle-level blind spot",
        "High-score inclusive backgrounds sit away from the photon-core control cloud.",
    )
    axp = inset(fig, 0.100, 0.180, 0.640, 0.630)
    style_axis(axp)
    axp.axvspan(0.80, 1.00, color="#FEE2E2", alpha=0.55, zorder=0)
    axp.axhspan(0.00, 0.20, color="#FFF7ED", alpha=0.75, zorder=0)
    axp.scatter(true["auau_tight_bdt_score"], true["final_state_photon_frac"], s=22, color=GREEN, alpha=0.20, edgecolors="none", label="truth photon controls")
    scatter_by_source_truth(axp, inc, "auau_tight_bdt_score", "final_state_photon_frac", size_col="cluster_Et")
    axp.axvline(0.80, color=RED, linestyle="--", linewidth=1.5)
    axp.axhline(0.80, color=GREEN, linestyle="--", linewidth=1.5)
    axp.set_xlim(0.45, 1.0)
    axp.set_ylim(-0.04, 1.04)
    axp.set_xlabel("BDT score", fontsize=15, fontweight="bold")
    axp.set_ylabel("final-state photon pT fraction in cone", fontsize=15, fontweight="bold")
    axp.tick_params(labelsize=12.5)
    non = int(inc["cone_bucket"].eq("non-photon-rich").sum())
    total = len(inc)
    med_true = float(true["auau_tight_bdt_score"].median())
    med_inc = float(inc["auau_tight_bdt_score"].median())
    card(ax, 0.775, 0.535, 0.175, 0.275)
    txt(ax, 0.795, 0.775, "Read the overlay", size=18, weight="bold")
    txt(ax, 0.795, 0.725, f"{non}/{total}", size=32, color=RED, weight="bold")
    txt(ax, 0.795, 0.675, "high-BDT inclusive\nbackgrounds are\nnon-photon-rich", size=15.5, linespacing=1.05)
    txt(ax, 0.795, 0.595, f"median BDT: fake {med_inc:.2f}\ncontrol {med_true:.2f}", size=15, color=MUTED)
    legend_box(
        ax,
        0.775,
        0.450,
        [("Jet20 marker", SOURCE_COLORS["Jet20"], SOURCE_MARKERS["Jet20"]), ("Jet30 marker", SOURCE_COLORS["Jet30"], SOURCE_MARKERS["Jet30"]), ("Jet40 marker", SOURCE_COLORS["Jet40"], SOURCE_MARKERS["Jet40"])],
    )
    txt(ax, 0.155, 0.125, "The confusing region is high BDT score but low photon-core fraction.", size=19, weight="bold", color=RED)
    return finish(fig, outdir / "overlay01_score_vs_particle_verdict.png")


def slide_02_energy_score_slices(inc: pd.DataFrame, true: pd.DataFrame, outdir: Path) -> Path:
    fig, ax = canvas()
    title(
        ax,
        "Energy slices show Jet30/40 supplying the hard high-score background",
        "Marker color is generator source; marker shape is truth source; size follows local non-photon fraction.",
    )
    axp = inset(fig, 0.090, 0.185, 0.750, 0.630)
    style_axis(axp)
    axp.scatter(true["cluster_Et"], true["auau_tight_bdt_score"], s=18, color="#9CA3AF", alpha=0.22, edgecolors="none", label="truth photon controls")
    for left, right, label in [(12, 18, "12-18"), (18, 24, "18-24"), (24, 35, "24-35")]:
        axp.axvspan(left, right, color="#EEF2F7" if label != "24-35" else "#FFF7ED", alpha=0.65, zorder=0)
        axp.text((left + right) / 2, 0.485, f"{label} GeV", ha="center", va="center", fontsize=12.5, color=MUTED)
    for source in ["Jet20", "Jet30", "Jet40"]:
        for truth in ["pi0/eta", "gamma truth", "other"]:
            sub = inc[(inc["source"].eq(source)) & (inc["truth_bucket"].eq(truth))]
            if sub.empty:
                continue
            sizes = np.clip(45 + 85 * sub["final_state_non_photon_frac"], 45, 140)
            axp.scatter(
                sub["cluster_Et"],
                sub["auau_tight_bdt_score"],
                s=sizes,
                marker=TRUTH_MARKERS[truth],
                facecolors=SOURCE_COLORS[source],
                edgecolors=TRUTH_COLORS[truth],
                linewidths=1.3,
                alpha=0.82,
            )
        src = inc[inc["source"].eq(source)]
        axp.scatter(
            [src["cluster_Et"].median()],
            [src["auau_tight_bdt_score"].median()],
            marker="D",
            s=230,
            color=SOURCE_COLORS[source],
            edgecolors="white",
            linewidths=1.8,
            zorder=10,
        )
    axp.set_xlim(12, 35)
    axp.set_ylim(0.48, 0.985)
    axp.set_xlabel("candidate cluster E$_T$ [GeV]", fontsize=15, fontweight="bold")
    axp.set_ylabel("BDT score", fontsize=15, fontweight="bold")
    axp.tick_params(labelsize=12.5)
    jet3040 = int(inc["source"].isin(["Jet30", "Jet40"]).sum())
    card(ax, 0.855, 0.525, 0.115, 0.290)
    txt(ax, 0.875, 0.775, f"{jet3040}/{len(inc)}", size=28, color=PURPLE, weight="bold")
    txt(ax, 0.875, 0.710, "high-BDT\ninclusive\nbackgrounds\ncome from\nJet30/40", size=15.5, linespacing=1.06)
    legend_box(
        ax,
        0.855,
        0.435,
        [("pi0/eta", TRUTH_COLORS["pi0/eta"], TRUTH_MARKERS["pi0/eta"]), ("gamma truth", TRUTH_COLORS["gamma truth"], TRUTH_MARKERS["gamma truth"]), ("other", TRUTH_COLORS["other"], TRUTH_MARKERS["other"])],
    )
    txt(ax, 0.145, 0.122, "The larger jet samples are not just more points; they populate the hard high-score background phase space.", size=18.5, weight="bold")
    return finish(fig, outdir / "overlay02_energy_score_slices.png")


def ternary_xy(df: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    p = df["final_state_photon_frac"].to_numpy(dtype=float)
    n = df["final_neutral_meson_frac"].to_numpy(dtype=float)
    c = df["final_charged_hadron_frac"].to_numpy(dtype=float)
    total = p + n + c
    total[total <= 0] = 1.0
    p, n, c = p / total, n / total, c / total
    neutral = np.array([0.10, 0.12])
    charged = np.array([0.90, 0.12])
    photon = np.array([0.50, 0.82])
    coords = n[:, None] * neutral + c[:, None] * charged + p[:, None] * photon
    return coords[:, 0], coords[:, 1]


def slide_03_composition_triangle(inc: pd.DataFrame, true: pd.DataFrame, outdir: Path) -> Path:
    fig, ax = canvas()
    title(
        ax,
        "The local cone composition separates photon cores from fragmentation cones",
        "Each marker is one candidate projected into photon / neutral-meson / charged-hadron pT fractions.",
    )
    tri_ax = inset(fig, 0.105, 0.120, 0.660, 0.700)
    tri_ax.set_axis_off()
    tri_ax.set_xlim(0, 1)
    tri_ax.set_ylim(0, 0.92)
    verts = np.array([[0.10, 0.12], [0.90, 0.12], [0.50, 0.82]])
    tri_ax.add_patch(Polygon(verts, closed=True, facecolor=CARD, edgecolor="#CBD5E1", linewidth=2.0, zorder=0))
    for frac in [0.25, 0.50, 0.75]:
        tri_ax.plot([0.10 + 0.40 * frac, 0.90 - 0.40 * frac], [0.12 + 0.70 * frac, 0.12 + 0.70 * frac], color=GRID, linewidth=0.9, zorder=1)
        tri_ax.plot([0.10 + 0.80 * frac, 0.50 + 0.40 * frac], [0.12, 0.82 - 0.70 * frac], color=GRID, linewidth=0.9, zorder=1)
        tri_ax.plot([0.90 - 0.80 * frac, 0.50 - 0.40 * frac], [0.12, 0.82 - 0.70 * frac], color=GRID, linewidth=0.9, zorder=1)
    x, y = ternary_xy(true)
    tri_ax.scatter(x, y, s=22, color=GREEN, alpha=0.22, edgecolors="none", zorder=2)
    for source in ["Jet20", "Jet30", "Jet40"]:
        sub = inc[inc["source"].eq(source)]
        sx, sy = ternary_xy(sub)
        sizes = np.clip(35 + 170 * (sub["auau_tight_bdt_score"] - 0.80), 45, 150)
        tri_ax.scatter(sx, sy, s=sizes, marker=SOURCE_MARKERS[source], color=SOURCE_COLORS[source], edgecolors="white", linewidths=1.0, alpha=0.86, zorder=5)
    tri_ax.text(0.50, 0.865, "photon core", ha="center", va="bottom", fontsize=18, fontweight="bold", color=GREEN)
    tri_ax.text(0.075, 0.075, "neutral meson", ha="left", va="top", fontsize=17, fontweight="bold", color=ORANGE)
    tri_ax.text(0.925, 0.075, "charged hadron", ha="right", va="top", fontsize=17, fontweight="bold", color=BLUE)
    med_x, med_y = ternary_xy(pd.DataFrame([inc.median(numeric_only=True)]))
    tri_ax.scatter(med_x, med_y, s=420, marker="*", color=RED, edgecolors="white", linewidths=1.5, zorder=9)
    tri_ax.text(med_x.item() + 0.025, med_y.item() + 0.018, "median high-BDT fake", fontsize=15, color=RED, fontweight="bold")
    card(ax, 0.790, 0.520, 0.165, 0.300)
    txt(ax, 0.812, 0.782, "Interpretation", size=19, weight="bold")
    txt(ax, 0.812, 0.725, "True photons pile up\nat the photon corner.", size=15.5, color=GREEN, linespacing=1.05)
    txt(ax, 0.812, 0.625, "High-BDT fakes sit\nalong the meson/\ncharged-fragment edge.", size=15.5, color=RED, linespacing=1.05)
    legend_box(ax, 0.790, 0.445, [("Jet20", SOURCE_COLORS["Jet20"], SOURCE_MARKERS["Jet20"]), ("Jet30", SOURCE_COLORS["Jet30"], SOURCE_MARKERS["Jet30"]), ("Jet40", SOURCE_COLORS["Jet40"], SOURCE_MARKERS["Jet40"])])
    txt(ax, 0.170, 0.060, "This plot turns the particle list into a visual fingerprint: photon-like detector score, non-photon local composition.", size=18.5, weight="bold")
    return finish(fig, outdir / "overlay03_cone_composition_triangle.png")


def slide_04_score_distribution_overlay(inc: pd.DataFrame, true: pd.DataFrame, outdir: Path) -> Path:
    fig, ax = canvas()
    title(
        ax,
        "The BDT score split is clear when the particle slices are overlaid",
        "Filled distributions compare truth-photon controls with the high-BDT inclusive background truth sources.",
    )
    axp = inset(fig, 0.105, 0.205, 0.735, 0.585)
    style_axis(axp)
    bins = np.linspace(0.45, 1.0, 36)
    axp.hist(true["auau_tight_bdt_score"], bins=bins, density=True, histtype="stepfilled", alpha=0.28, color=GREEN, edgecolor=GREEN, linewidth=2.0, label="truth photon control")
    for truth in ["pi0/eta", "gamma truth", "other"]:
        sub = inc[inc["truth_bucket"].eq(truth)]
        axp.hist(sub["auau_tight_bdt_score"], bins=bins, density=True, histtype="step", linewidth=3.0, color=TRUTH_COLORS[truth], label=f"high-BDT bg: {truth}")
        axp.axvline(sub["auau_tight_bdt_score"].median(), color=TRUTH_COLORS[truth], linestyle="--", linewidth=1.6)
    axp.axvline(true["auau_tight_bdt_score"].median(), color=GREEN, linestyle="--", linewidth=1.8)
    axp.set_xlim(0.45, 1.0)
    axp.set_xlabel("BDT score", fontsize=15, fontweight="bold")
    axp.set_ylabel("normalized candidate density", fontsize=15, fontweight="bold")
    axp.tick_params(labelsize=12.5)
    axp.legend(loc="upper left", frameon=True, facecolor=CARD, edgecolor="#CBD5E1", fontsize=12.5)
    card(ax, 0.855, 0.450, 0.120, 0.340)
    txt(ax, 0.875, 0.750, "Median markers", size=17, weight="bold")
    entries = [
        ("truth photons", true["auau_tight_bdt_score"].median(), GREEN),
        ("pi0/eta", inc[inc["truth_bucket"].eq("pi0/eta")]["auau_tight_bdt_score"].median(), ORANGE),
        ("gamma truth", inc[inc["truth_bucket"].eq("gamma truth")]["auau_tight_bdt_score"].median(), RED),
        ("other", inc[inc["truth_bucket"].eq("other")]["auau_tight_bdt_score"].median(), SLATE),
    ]
    for i, (label, val, color) in enumerate(entries):
        yy = 0.700 - 0.060 * i
        ax.add_patch(Rectangle((0.875, yy - 0.018), 0.020, 0.020, transform=ax.transAxes, facecolor=color, edgecolor="none", zorder=5))
        txt(ax, 0.905, yy, f"{label}\n{val:.2f}", size=13.5, va="center", linespacing=1.0)
    txt(ax, 0.150, 0.125, "This is why median score gap was more visually faithful than AUC for the slide-21 question.", size=18.5, weight="bold")
    return finish(fig, outdir / "overlay04_score_distribution_truth_slices.png")


def slide_05_particle_microscope_overlay(part: pd.DataFrame, examples: pd.DataFrame, outdir: Path) -> Path:
    fig, ax = canvas()
    title(
        ax,
        "A fine-grained particle slice shows what sits near the EM cluster",
        "Status-1 particles inside Delta R < 0.3: marker color is particle type, marker area follows pT.",
    )
    fake = examples[examples["example_role"].eq("highest-score non-photon background")].sort_values("auau_tight_bdt_score", ascending=False).iloc[0]
    true = examples[examples["example_role"].eq("busy true photon with lowered score")].sort_values("auau_tight_bdt_score", ascending=False).iloc[0]

    def draw_panel(x0: float, row: pd.Series, header: str, color: str) -> None:
        axp = inset(fig, x0, 0.240, 0.390, 0.520)
        style_axis(axp)
        sub = part[(part["candidate_uid"].eq(row["candidate_uid"])) & (part["status"].eq(1))].copy()
        for group, g in sub.groupby("group"):
            axp.scatter(g["dr"], g["pt"], s=np.clip(25 + 8 * g["pt"], 45, 360), color=GROUP_COLORS.get(group, SLATE), alpha=0.78, edgecolors="white", linewidths=1.0, label=group.replace("_", " "))
            for _, r in g.sort_values("pt", ascending=False).head(2).iterrows():
                axp.text(r["dr"] + 0.006, r["pt"] + 0.3, str(r["pdg_name"]), fontsize=10.5, color=GROUP_COLORS.get(group, SLATE), fontweight="bold")
        for d in [0.1, 0.2, 0.3]:
            axp.axvline(d, color=GRID, linestyle="--", linewidth=1.0)
        axp.set_xlim(0, 0.305)
        axp.set_ylim(0, max(16, sub["pt"].max() * 1.18))
        axp.set_xlabel("Delta R from cluster", fontsize=13.5, fontweight="bold")
        axp.set_ylabel("status-1 particle pT [GeV]", fontsize=13.5, fontweight="bold")
        axp.tick_params(labelsize=11.5)
        txt(ax, x0, 0.820, header, size=22, weight="bold", color=color)
        txt(ax, x0, 0.785, f"{row['sample']}  truth {row['truth_pid_name']}  BDT {row['auau_tight_bdt_score']:.3f}", size=14.5, color=INK, weight="bold")
        vals = [
            ("photon", row["final_state_photon_frac"], RED),
            ("neutral meson", row["final_neutral_meson_frac"], ORANGE),
            ("charged hadron", row["final_charged_hadron_frac"], BLUE),
        ]
        for i, (lab, val, col) in enumerate(vals):
            yy = 0.185 - 0.034 * i
            ax.add_patch(Rectangle((x0, yy - 0.015), 0.050 * float(val), 0.018, transform=ax.transAxes, facecolor=col, edgecolor="none", zorder=5))
            txt(ax, x0 + 0.060, yy, f"{lab}: {100 * float(val):.0f}%", size=12.8, va="center")

    draw_panel(0.085, fake, "High-BDT Jet40 fake", RED)
    draw_panel(0.535, true, "Truth-photon control", GREEN)
    txt(ax, 0.500, 0.070, "The fake panel is a local fragmentation cone; the control panel is a photon core with nearby activity.", size=18.5, weight="bold", ha="center")
    return finish(fig, outdir / "overlay05_dr_pt_particle_microscope.png")


def slide_06_source_truth_bubble_atlas(inc: pd.DataFrame, outdir: Path) -> Path:
    fig, ax = canvas()
    title(
        ax,
        "A source/truth bubble atlas compresses the high-BDT fake taxonomy",
        "Bubble area is candidate count; color intensity is median non-photon fraction; numbers are count and median BDT.",
    )
    axp = inset(fig, 0.120, 0.175, 0.660, 0.635)
    axp.set_facecolor(CARD)
    sources = ["Jet20", "Jet30", "Jet40"]
    truths = ["pi0/eta", "gamma truth", "other"]
    axp.set_xlim(-0.5, len(truths) - 0.5)
    axp.set_ylim(-0.5, len(sources) - 0.5)
    axp.set_xticks(range(len(truths)), truths, fontsize=15, fontweight="bold")
    axp.set_yticks(range(len(sources)), sources, fontsize=15, fontweight="bold")
    axp.grid(True, color=GRID, linewidth=1.0)
    for spine in axp.spines.values():
        spine.set_color("#CBD5E1")
    max_count = max(inc.groupby(["source", "truth_bucket"]).size().max(), 1)
    for yi, source in enumerate(sources):
        for xi, truth in enumerate(truths):
            sub = inc[(inc["source"].eq(source)) & (inc["truth_bucket"].eq(truth))]
            n = len(sub)
            if n == 0:
                continue
            med_non = float(sub["final_state_non_photon_frac"].median())
            med_score = float(sub["auau_tight_bdt_score"].median())
            size = 500 + 2300 * n / max_count
            color = TRUTH_COLORS[truth]
            axp.scatter([xi], [yi], s=size, color=color, alpha=0.32 + 0.55 * med_non, edgecolors=SOURCE_COLORS[source], linewidths=2.2, zorder=5)
            axp.text(xi, yi, f"{n}\n{med_score:.2f}", ha="center", va="center", fontsize=14.5, fontweight="bold", color=INK, zorder=8)
    axp.set_xlabel("truth source", fontsize=16, fontweight="bold")
    axp.set_ylabel("generator parent sample", fontsize=16, fontweight="bold")
    pi0eta = int(inc["truth_bucket"].eq("pi0/eta").sum())
    jet3040_pi0 = int(inc[inc["source"].isin(["Jet30", "Jet40"])]["truth_bucket"].eq("pi0/eta").sum())
    card(ax, 0.815, 0.545, 0.150, 0.270)
    txt(ax, 0.837, 0.775, "Dominant cell family", size=17, weight="bold")
    txt(ax, 0.837, 0.715, f"{jet3040_pi0}", size=34, color=ORANGE, weight="bold")
    txt(ax, 0.837, 0.655, "Jet30/40 plus\npi0/eta truth", size=15.5, linespacing=1.05)
    txt(ax, 0.837, 0.585, f"all pi0/eta:\n{pi0eta}/{len(inc)}", size=15, color=MUTED)
    card(ax, 0.815, 0.265, 0.150, 0.205, face=YELLOW, edge="#F2CB4C")
    txt(ax, 0.837, 0.425, "Slide answer", size=17, weight="bold")
    txt(ax, 0.837, 0.375, "The expanded\ntraining samples add\nspecific hard neutral-\nmeson backgrounds.", size=14.5, linespacing=1.04)
    txt(ax, 0.170, 0.095, "This is the compact taxonomy view: source energy on rows, truth identity on columns, BDT response printed inside.", size=18.5, weight="bold")
    return finish(fig, outdir / "overlay06_source_truth_bubble_atlas.png")


def contact_sheet(paths: list[Path], outdir: Path) -> Path:
    fig = plt.figure(figsize=(16, 9), dpi=160)
    fig.patch.set_facecolor(LIGHT)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    txt(ax, 0.045, 0.955, "THE-44 Overlay + Slice Candidate Set", size=30, weight="bold")
    txt(ax, 0.045, 0.905, "Six overlay-style slide candidates: markers, slices, distributions, and candidate-level particle microscopes.", size=17, color=MUTED)
    labels = [
        "1. Score vs particle verdict",
        "2. Energy score slices",
        "3. Cone composition triangle",
        "4. Score distribution overlay",
        "5. DeltaR-pT particle microscope",
        "6. Source/truth bubble atlas",
    ]
    positions = [(0.06, 0.525), (0.38, 0.525), (0.70, 0.525), (0.06, 0.115), (0.38, 0.115), (0.70, 0.115)]
    for p, label, (x, y) in zip(paths, labels, positions):
        txt(ax, x, y + 0.315, label, size=16.5, weight="bold")
        img_ax = fig.add_axes([x, y, 0.27, 0.255])
        img_ax.imshow(mpimg.imread(p))
        img_ax.set_xticks([])
        img_ax.set_yticks([])
        for spine in img_ax.spines.values():
            spine.set_color("#CBD5E1")
            spine.set_linewidth(1.0)
    out = outdir / "the44_overlay_slice_candidate_contact_sheet.png"
    fig.savefig(out, dpi=160)
    plt.close(fig)
    return out


def write_scripts(outdir: Path) -> dict[str, str]:
    scripts = {
        "overlay01_score_vs_particle_verdict_script.md": "This overlay puts BDT score on x and final-state photon fraction on y. The truth-photon controls cluster near photon fraction one, while high-BDT inclusive backgrounds often sit at high BDT but low photon fraction. That is the particle-level blind spot the autopsy exposes.\n",
        "overlay02_energy_score_slices_script.md": "This slide overlays candidate cluster energy and BDT score. Jet30 and Jet40 markers dominate the high-score inclusive background population, showing that the larger jet samples supply the hard background phase space the BDT must learn.\n",
        "overlay03_cone_composition_triangle_script.md": "The ternary view projects each candidate into photon, neutral-meson, and charged-hadron final-state pT fractions. Truth photons pile up at the photon corner; high-BDT fakes move along the meson/charged fragmentation edge.\n",
        "overlay04_score_distribution_truth_slices_script.md": "The score distribution overlay compares truth-photon controls with high-BDT background truth slices. The median markers make the qualitative score split visible in one axis, matching the median-gap argument better than AUC alone.\n",
        "overlay05_dr_pt_particle_microscope_script.md": "This candidate microscope uses status-1 particles inside Delta R less than 0.3. The Jet40 fake has neutral-meson and charged-hadron local content; the truth-photon control has a photon core. Marker size follows particle pT.\n",
        "overlay06_source_truth_bubble_atlas_script.md": "The bubble atlas compresses the taxonomy: rows are generator source, columns are truth source, bubble size is count, and printed values are count plus median BDT. The dominant lesson is Jet30/40 plus pi0/eta truth.\n",
    }
    out: dict[str, str] = {}
    for name, body in scripts.items():
        path = outdir / name
        path.write_text(body)
        out[name] = str(path)
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    _, inc, true, part, examples = load_data()
    paths = [
        slide_01_score_particle_overlay(inc, true, args.outdir),
        slide_02_energy_score_slices(inc, true, args.outdir),
        slide_03_composition_triangle(inc, true, args.outdir),
        slide_04_score_distribution_overlay(inc, true, args.outdir),
        slide_05_particle_microscope_overlay(part, examples, args.outdir),
        slide_06_source_truth_bubble_atlas(inc, args.outdir),
    ]
    contact = contact_sheet(paths, args.outdir)
    scripts = write_scripts(args.outdir)
    metrics = {
        "inclusive_high_bdt_backgrounds": float(len(inc)),
        "truth_photon_controls": float(len(true)),
        "jet30_40_candidates": float(inc["source"].isin(["Jet30", "Jet40"]).sum()),
        "pi0_eta_truth_candidates": float(inc["truth_bucket"].eq("pi0/eta").sum()),
        "non_photon_rich_cones": float(inc["cone_bucket"].eq("non-photon-rich").sum()),
        "median_high_fake_bdt_score": float(inc["auau_tight_bdt_score"].median()),
        "median_truth_photon_bdt_score": float(true["auau_tight_bdt_score"].median()),
    }
    manifest = {
        "status": "READY",
        "slides": [str(p) for p in paths],
        "contact_sheet": str(contact),
        "scripts": scripts,
        "metrics": metrics,
        "inputs": {
            "candidate_rows": str(THE44_DIR / "the44_candidate_autopsy_rows.csv"),
            "particle_rows": str(THE44_DIR / "the44_particle_rows.csv"),
            "event_cards": str(THE44_DIR / "the44_event_card_examples.csv"),
        },
        "caveat": "Internal qualitative overlay diagnostic from saved THE-44 rows; not a production-weighted purity estimate.",
    }
    manifest_path = args.outdir / "the44_overlay_slice_candidate_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
