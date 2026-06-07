#!/usr/bin/env python3
"""Build THE-44 triangle/evolution slides for local-cone composition."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Polygon, Rectangle


REPO = Path(__file__).resolve().parents[4]
BASE = REPO / "dataOutput/auauTightBDTValidation/THE8_branchA_ladder_scorecache_fullstat_20260527"
THE44_DIR = BASE / "diagnostics/the44_high_bdt_pythia_autopsy_20260605/slideReady/full_autopsy_20260605_2108"
DEFAULT_OUTDIR = THE44_DIR / "particle_intensive/triangle_evolution_sequence_20260606"

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
SLATE = "#475569"
YELLOW = "#FFF4C7"

SOURCE_COLORS = {"Jet20": "#F97316", "Jet30": GREEN, "Jet40": PURPLE}
SOURCE_MARKERS = {"Jet20": "o", "Jet30": "s", "Jet40": "^"}
GROUP_COLORS = {"photon": RED, "neutral_meson": ORANGE, "charged_hadron": BLUE, "quark_gluon": PURPLE, "other": SLATE}

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


def txt(ax: plt.Axes, x: float, y: float, text: str, *, size: float, color: str = INK, weight: str = "normal", ha: str = "left", va: str = "top", linespacing: float = 1.08) -> None:
    ax.text(x, y, text, transform=ax.transAxes, fontsize=size, color=color, fontweight=weight, ha=ha, va=va, linespacing=linespacing, zorder=30)


def card(ax: plt.Axes, x: float, y: float, w: float, h: float, *, face: str = CARD, edge: str = "#D4DCE7", lw: float = 1.2, radius: float = 0.018, zorder: int = 1) -> None:
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
    return cand, inc, true, part, examples


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


def draw_triangle(fig: plt.Figure, pos: tuple[float, float, float, float]) -> plt.Axes:
    ax = fig.add_axes(pos)
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 0.92)
    verts = np.array([[0.10, 0.12], [0.90, 0.12], [0.50, 0.82]])
    ax.add_patch(Polygon(verts, closed=True, facecolor=CARD, edgecolor="#CBD5E1", linewidth=2.0, zorder=0))
    for frac in [0.25, 0.50, 0.75]:
        ax.plot([0.10 + 0.40 * frac, 0.90 - 0.40 * frac], [0.12 + 0.70 * frac, 0.12 + 0.70 * frac], color=GRID, linewidth=0.9, zorder=1)
        ax.plot([0.10 + 0.80 * frac, 0.50 + 0.40 * frac], [0.12, 0.82 - 0.70 * frac], color=GRID, linewidth=0.9, zorder=1)
        ax.plot([0.90 - 0.80 * frac, 0.50 - 0.40 * frac], [0.12, 0.82 - 0.70 * frac], color=GRID, linewidth=0.9, zorder=1)
    ax.text(0.50, 0.860, "photon core", ha="center", va="bottom", fontsize=18, fontweight="bold", color=GREEN)
    ax.text(0.070, 0.075, "neutral meson", ha="left", va="top", fontsize=17, fontweight="bold", color=ORANGE)
    ax.text(0.930, 0.075, "charged hadron", ha="right", va="top", fontsize=17, fontweight="bold", color=BLUE)
    return ax


def med_point(df: pd.DataFrame) -> tuple[float, float]:
    vals = pd.DataFrame(
        [
            {
                "final_state_photon_frac": df["final_state_photon_frac"].median(),
                "final_neutral_meson_frac": df["final_neutral_meson_frac"].median(),
                "final_charged_hadron_frac": df["final_charged_hadron_frac"].median(),
            }
        ]
    )
    x, y = ternary_xy(vals)
    return x.item(), y.item()


def arrow(ax: plt.Axes, start: tuple[float, float], end: tuple[float, float], color: str, *, lw: float = 2.2, alpha: float = 0.90) -> None:
    ax.add_patch(FancyArrowPatch(start, end, transform=ax.transData, arrowstyle="-|>", mutation_scale=18, linewidth=lw, color=color, alpha=alpha, zorder=10))


def slide_01_phase_space(inc: pd.DataFrame, true: pd.DataFrame, outdir: Path) -> Path:
    fig, ax = canvas()
    title(
        ax,
        "The local cone has two physical attractors: photon core and fragmentation edge",
        "A ternary phase space turns the Pythia particle list into a direct visual fingerprint of what the BDT sees.",
    )
    tri = draw_triangle(fig, (0.055, 0.090, 0.670, 0.770))
    x, y = ternary_xy(true)
    tri.scatter(x, y, s=24, color=GREEN, alpha=0.20, edgecolors="none", zorder=2)
    for source in ["Jet20", "Jet30", "Jet40"]:
        sub = inc[inc["source"].eq(source)]
        sx, sy = ternary_xy(sub)
        sizes = np.clip(45 + 150 * (sub["auau_tight_bdt_score"] - 0.80), 45, 150)
        tri.scatter(sx, sy, s=sizes, marker=SOURCE_MARKERS[source], color=SOURCE_COLORS[source], edgecolors="white", linewidths=0.9, alpha=0.86, zorder=5)
    fake_x, fake_y = med_point(inc)
    true_x, true_y = med_point(true)
    tri.scatter([fake_x], [fake_y], s=440, marker="*", color=RED, edgecolors="white", linewidths=1.5, zorder=11)
    tri.scatter([true_x], [true_y], s=320, marker="*", color=GREEN, edgecolors="white", linewidths=1.5, zorder=11)
    tri.text(fake_x + 0.025, fake_y + 0.018, "median high-BDT fake", fontsize=15, color=RED, fontweight="bold")
    tri.text(true_x + 0.020, true_y + 0.022, "truth-photon median", fontsize=15, color=GREEN, fontweight="bold")
    card(ax, 0.760, 0.535, 0.185, 0.295)
    txt(ax, 0.785, 0.785, "What it says", size=19, weight="bold")
    txt(ax, 0.785, 0.730, "Truth photons collect\nat the photon-core\ncorner.", size=15.5, color=GREEN)
    txt(ax, 0.785, 0.620, "High-BDT fakes run\nalong the meson/\ncharged-fragment\nedge.", size=14.7, color=RED, linespacing=1.0)
    card(ax, 0.760, 0.330, 0.185, 0.155, face=YELLOW, edge="#F2CB4C")
    txt(ax, 0.785, 0.445, "Slide answer", size=18, weight="bold")
    txt(ax, 0.785, 0.395, "The BDT score is high,\nbut the particle cone\nis not photon-core-like.", size=14.8)
    return finish(fig, outdir / "triangle01_local_cone_phase_space.png")


def slide_02_source_vectors(inc: pd.DataFrame, true: pd.DataFrame, outdir: Path) -> Path:
    fig, ax = canvas()
    title(
        ax,
        "Jet30/40 do not move the fake class to photons; they fill the hard fragmentation edge",
        "Source medians are drawn as an evolution path across parent jet pThat bins inside the same cone-composition triangle.",
    )
    tri = draw_triangle(fig, (0.055, 0.090, 0.620, 0.760))
    x, y = ternary_xy(true)
    tri.scatter(x, y, s=16, color=GREEN, alpha=0.13, edgecolors="none", zorder=2)
    source_points: dict[str, tuple[float, float]] = {}
    for source in ["Jet20", "Jet30", "Jet40"]:
        sub = inc[inc["source"].eq(source)]
        sx, sy = ternary_xy(sub)
        tri.scatter(sx, sy, s=42, marker=SOURCE_MARKERS[source], color=SOURCE_COLORS[source], edgecolors="white", linewidths=0.7, alpha=0.50, zorder=4)
        mx, my = med_point(sub)
        source_points[source] = (mx, my)
        tri.scatter([mx], [my], s=380, marker=SOURCE_MARKERS[source], color=SOURCE_COLORS[source], edgecolors="white", linewidths=1.5, zorder=12)
        tri.text(mx + 0.018, my + 0.028, f"{source}\nmedian", fontsize=13.5, color=SOURCE_COLORS[source], fontweight="bold", zorder=13)
    arrow(tri, source_points["Jet20"], source_points["Jet30"], SOURCE_COLORS["Jet30"], lw=2.0)
    arrow(tri, source_points["Jet30"], source_points["Jet40"], SOURCE_COLORS["Jet40"], lw=2.0)
    rows = []
    for source in ["Jet20", "Jet30", "Jet40"]:
        sub = inc[inc["source"].eq(source)]
        rows.append((source, len(sub), sub["cluster_Et"].median(), sub["auau_tight_bdt_score"].median(), sub["final_state_non_photon_frac"].median(), SOURCE_COLORS[source]))
    card(ax, 0.700, 0.505, 0.240, 0.320)
    txt(ax, 0.725, 0.785, "Source medians", size=19, weight="bold")
    for i, (source, n, et, score, non, color) in enumerate(rows):
        yy = 0.720 - 0.075 * i
        ax.scatter([0.730], [yy - 0.012], transform=ax.transAxes, s=115, color=color, marker=SOURCE_MARKERS[source], edgecolors="white", linewidths=1.1, zorder=20)
        txt(ax, 0.760, yy, f"{source}: n={n}, E$_T$ {et:.1f} GeV\nmedian BDT {score:.2f}, non-photon {non:.0%}", size=13.8, va="top", linespacing=1.0)
    card(ax, 0.700, 0.285, 0.240, 0.145, face=YELLOW, edge="#F2CB4C")
    txt(ax, 0.725, 0.392, "Physical read", size=18, weight="bold")
    txt(ax, 0.725, 0.345, "The added samples broaden the\nhard-fragmentation examples the\nclassifier learns to separate.", size=14.6)
    return finish(fig, outdir / "triangle02_source_evolution_vectors.png")


def slide_03_fragmentation_tree(inc: pd.DataFrame, part: pd.DataFrame, outdir: Path) -> Path:
    fig, ax = canvas()
    title(
        ax,
        "The triangle can be read as a fragmentation history ending in a local cone",
        "HepMC ancestry supplies the upstream process; status-1 particles define the final cone point.",
    )
    tri = draw_triangle(fig, (0.520, 0.130, 0.405, 0.600))
    sub = inc[(inc["source"].isin(["Jet30", "Jet40"])) & inc["truth_bucket"].eq("pi0/eta")]
    sx, sy = ternary_xy(sub)
    tri.scatter(sx, sy, s=42, color=ORANGE, alpha=0.55, edgecolors="white", linewidths=0.6, zorder=4)
    mx, my = med_point(sub)
    tri.scatter([mx], [my], s=430, marker="*", color=RED, edgecolors="white", linewidths=1.5, zorder=10)
    tri.text(mx + 0.020, my + 0.020, "Jet30/40\npi0/eta median", fontsize=14, color=RED, fontweight="bold")
    psub = part[part["candidate_uid"].isin(sub["candidate_uid"])]
    stage_rows = [
        ("1", "hard source", "Jet30/40", f"{len(sub)} candidates", PURPLE),
        ("2", "truth seed", "pi0/eta", "neutral meson source", ORANGE),
        ("3", "final cone", "status-1", f"{int(sub['cone_bucket'].eq('non-photon-rich').sum())}/{len(sub)} non-photon", BLUE),
        ("4", "classifier", "BDT 0.84", "high-score tail", RED),
    ]
    x0, y0, w, h, gap = 0.070, 0.635, 0.095, 0.155, 0.018
    for i, (num, head, val, detail, color) in enumerate(stage_rows):
        xx = x0 + i * (w + gap)
        card(ax, xx, y0, w, h, face=CARD, edge=color, lw=1.6)
        txt(ax, xx + 0.014, y0 + 0.128, num, size=16, color=color, weight="bold")
        txt(ax, xx + 0.014, y0 + 0.100, head, size=12.3, color=color, weight="bold")
        txt(ax, xx + 0.014, y0 + 0.066, val, size=17, weight="bold")
        txt(ax, xx + 0.014, y0 + 0.030, detail, size=10.8, color=MUTED, linespacing=0.95)
        if i < len(stage_rows) - 1:
            ax.add_patch(FancyArrowPatch((xx + w + 0.003, y0 + h / 2), (xx + w + gap - 0.004, y0 + h / 2), transform=ax.transAxes, arrowstyle="-|>", mutation_scale=15, linewidth=1.8, color=SLATE, alpha=0.75, zorder=15))
    ax.add_patch(FancyArrowPatch((0.305, 0.620), (0.535, 0.365), transform=ax.transAxes, arrowstyle="-|>", mutation_scale=20, linewidth=2.2, color=SLATE, alpha=0.72, zorder=12))
    final = psub[psub["status"].eq(1)].groupby("group")["pt"].sum()
    total = float(final.sum())
    card(ax, 0.085, 0.165, 0.330, 0.225, face=YELLOW, edge="#F2CB4C")
    txt(ax, 0.110, 0.345, "Final status-1 pT budget", size=17, weight="bold")
    for i, group in enumerate(["neutral_meson", "charged_hadron", "photon"]):
        val = float(final.get(group, 0.0))
        yy = 0.295 - 0.055 * i
        ax.add_patch(Rectangle((0.110, yy - 0.016), 0.185 * val / total, 0.026, transform=ax.transAxes, facecolor=GROUP_COLORS[group], edgecolor="none", zorder=10))
        txt(ax, 0.315, yy, f"{group.replace('_', ' ')} {100 * val / total:.0f}%", size=13.5, va="center")
    return finish(fig, outdir / "triangle03_fragmentation_tree_to_cone.png")


def slide_04_candidate_time_path(examples: pd.DataFrame, outdir: Path) -> Path:
    fig, ax = canvas()
    title(
        ax,
        "One candidate path: neutral-meson fragmentation can reconstruct like a photon",
        "The ternary endpoint explains why the event is physically a fake even when detector variables look EM-like.",
    )
    fake = examples[examples["example_role"].eq("highest-score non-photon background")].sort_values("auau_tight_bdt_score", ascending=False).iloc[0]
    true = examples[examples["example_role"].eq("busy true photon with lowered score")].sort_values("auau_tight_bdt_score", ascending=False).iloc[0]
    tri = draw_triangle(fig, (0.560, 0.150, 0.350, 0.560))
    fx, fy = ternary_xy(pd.DataFrame([fake]))
    tx, ty = ternary_xy(pd.DataFrame([true]))
    tri.scatter([tx.item()], [ty.item()], s=300, marker="*", color=GREEN, edgecolors="white", linewidths=1.4, zorder=10)
    tri.scatter([fx.item()], [fy.item()], s=360, marker="*", color=RED, edgecolors="white", linewidths=1.4, zorder=11)
    tri.text(tx.item() + 0.020, ty.item() + 0.020, "truth photon\ncontrol", fontsize=13.5, color=GREEN, fontweight="bold")
    tri.text(fx.item() + 0.020, fy.item() + 0.020, "Jet40 pi0\nhigh-BDT fake", fontsize=13.5, color=RED, fontweight="bold")
    stages = [
        ("1. Parent sample", "Jet40 inclusive", "harder inclusive background bin", PURPLE),
        ("2. Truth seed", "pi0", "neutral meson near the cluster", ORANGE),
        ("3. Final cone", "0% photon\n92% neutral meson", "status-1 local particle endpoint", BLUE),
    ]
    for i, (head, val, detail, color) in enumerate(stages):
        yy = 0.650 - i * 0.130
        card(ax, 0.085, yy, 0.390, 0.100, face=CARD, edge=color, lw=1.5)
        txt(ax, 0.110, yy + 0.078, head, size=13.5, color=color, weight="bold")
        txt(ax, 0.110, yy + 0.046, val, size=16.8, weight="bold", linespacing=0.92)
        txt(ax, 0.300, yy + 0.046, detail, size=12.0, color=MUTED, linespacing=1.0)
        if i < len(stages) - 1:
            ax.add_patch(FancyArrowPatch((0.280, yy - 0.004), (0.280, yy - 0.028), transform=ax.transAxes, arrowstyle="-|>", mutation_scale=14, linewidth=1.6, color=SLATE, zorder=20))
    card(ax, 0.085, 0.160, 0.390, 0.220)
    txt(ax, 0.115, 0.342, f"BDT response {fake['auau_tight_bdt_score']:.3f}", size=18, color=RED, weight="bold")
    txt(ax, 0.115, 0.309, "Detector-level tension", size=16.5, weight="bold")
    rows = [
        ("BDT score", fake["auau_tight_bdt_score"], true["auau_tight_bdt_score"], 1.0),
        ("e11/e33", fake["e11_over_e33"], true["e11_over_e33"], 1.0),
        ("reco iso / 12", fake["reco_eiso"], true["reco_eiso"], 12.0),
        ("non-photon frac", fake["final_state_non_photon_frac"], true["final_state_non_photon_frac"], 1.0),
    ]
    for i, (label, fval, tval, scale) in enumerate(rows):
        yy = 0.270 - 0.038 * i
        txt(ax, 0.115, yy + 0.010, label, size=12.5, va="center")
        ax.add_patch(Rectangle((0.235, yy), 0.110, 0.018, transform=ax.transAxes, facecolor="#E8EDF3", edgecolor="none", zorder=2))
        ax.add_patch(Rectangle((0.235, yy), 0.110 * min(float(fval) / scale, 1.0), 0.018, transform=ax.transAxes, facecolor=RED, edgecolor="none", zorder=3))
        ax.add_patch(Rectangle((0.355, yy), 0.090, 0.018, transform=ax.transAxes, facecolor="#E8EDF3", edgecolor="none", zorder=2))
        ax.add_patch(Rectangle((0.355, yy), 0.090 * min(float(tval) / scale, 1.0), 0.018, transform=ax.transAxes, facecolor=GREEN, edgecolor="none", zorder=3))
    txt(ax, 0.110, 0.095, "The endpoint in the triangle is the missing physical context: high score does not mean photon-core origin.", size=18, weight="bold")
    return finish(fig, outdir / "triangle04_candidate_time_path.png")


def contact_sheet(paths: list[Path], outdir: Path) -> Path:
    fig = plt.figure(figsize=(16, 9), dpi=160)
    fig.patch.set_facecolor(LIGHT)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    txt(ax, 0.045, 0.955, "THE-44 Triangle Evolution Sequence", size=30, weight="bold")
    txt(ax, 0.045, 0.905, "Four triangle-centered candidates for local cone composition, source evolution, and HepMC fragmentation history.", size=17, color=MUTED)
    labels = [
        "1. Cone phase space",
        "2. Source evolution vectors",
        "3. Fragmentation tree to cone",
        "4. Candidate time path",
    ]
    positions = [(0.08, 0.515), (0.52, 0.515), (0.08, 0.105), (0.52, 0.105)]
    for p, label, (x, y) in zip(paths, labels, positions):
        txt(ax, x, y + 0.315, label, size=17, weight="bold")
        img_ax = fig.add_axes([x, y, 0.36, 0.260])
        img_ax.imshow(mpimg.imread(p))
        img_ax.set_xticks([])
        img_ax.set_yticks([])
        for spine in img_ax.spines.values():
            spine.set_color("#CBD5E1")
            spine.set_linewidth(1.0)
    out = outdir / "the44_triangle_evolution_contact_sheet.png"
    fig.savefig(out, dpi=160)
    plt.close(fig)
    return out


def write_scripts(outdir: Path) -> dict[str, str]:
    scripts = {
        "triangle01_local_cone_phase_space_script.md": "This slide introduces the ternary cone as a phase space. Truth photons cluster at the photon-core corner, while high-BDT inclusive backgrounds lie along the meson/charged fragmentation edge.\n",
        "triangle02_source_evolution_vectors_script.md": "This slide overlays Jet20, Jet30, and Jet40 source medians in the same triangle. The path shows that the added higher-jet samples fill the hard fragmentation edge rather than becoming photon-core objects.\n",
        "triangle03_fragmentation_tree_to_cone_script.md": "This slide reads the triangle as the endpoint of a HepMC production path: hard Jet30/40 source, pi0/eta truth seed, final status-1 local cone, and then BDT score.\n",
        "triangle04_candidate_time_path_script.md": "This slide makes the idea concrete with one Jet40 pi0 high-BDT fake. The detector variables can look photon-like, but the triangle endpoint says the final local particle cone is neutral-meson dominated.\n",
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
        slide_01_phase_space(inc, true, args.outdir),
        slide_02_source_vectors(inc, true, args.outdir),
        slide_03_fragmentation_tree(inc, part, args.outdir),
        slide_04_candidate_time_path(examples, args.outdir),
    ]
    contact = contact_sheet(paths, args.outdir)
    scripts = write_scripts(args.outdir)
    manifest = {
        "status": "READY",
        "slides": [str(p) for p in paths],
        "contact_sheet": str(contact),
        "scripts": scripts,
        "metrics": {
            "inclusive_high_bdt_backgrounds": float(len(inc)),
            "truth_photon_controls": float(len(true)),
            "jet30_40_candidates": float(inc["source"].isin(["Jet30", "Jet40"]).sum()),
            "pi0_eta_truth_candidates": float(inc["truth_bucket"].eq("pi0/eta").sum()),
            "non_photon_rich_cones": float(inc["cone_bucket"].eq("non-photon-rich").sum()),
        },
        "inputs": {
            "candidate_rows": str(THE44_DIR / "the44_candidate_autopsy_rows.csv"),
            "particle_rows": str(THE44_DIR / "the44_particle_rows.csv"),
            "event_cards": str(THE44_DIR / "the44_event_card_examples.csv"),
        },
        "caveat": "HepMC ancestry/status diagnostic from saved THE-44 particle rows; this is not a detector-time animation or a production-weighted purity estimate.",
    }
    manifest_path = args.outdir / "the44_triangle_evolution_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
