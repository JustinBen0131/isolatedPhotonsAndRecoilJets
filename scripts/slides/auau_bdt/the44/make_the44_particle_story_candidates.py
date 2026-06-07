#!/usr/bin/env python3
"""Build six THE-44 particle-history slide candidates for the BDT fake-tail story."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle


REPO = Path(__file__).resolve().parents[4]
BASE = REPO / "dataOutput/auauTightBDTValidation/THE8_branchA_ladder_scorecache_fullstat_20260527"
THE44_DIR = BASE / "diagnostics/the44_high_bdt_pythia_autopsy_20260605/slideReady/full_autopsy_20260605_2108"
DEFAULT_OUTDIR = THE44_DIR / "particle_intensive/deep_particle_story_candidates_20260606"

INK = "#111827"
MUTED = "#5D6675"
LIGHT = "#F6F7F9"
CARD = "#FFFFFF"
GRID = "#DEE5EE"
BLUE = "#2563A6"
RED = "#C43C32"
GREEN = "#2F8F52"
ORANGE = "#D97706"
PURPLE = "#7651A6"
TEAL = "#0F766E"
YELLOW = "#FFF4C7"
SLATE = "#475569"

SOURCE_COLORS = {"Jet20": "#F97316", "Jet30": GREEN, "Jet40": PURPLE}
TRUTH_COLORS = {"pi0/eta": ORANGE, "gamma truth": RED, "other": SLATE}
GROUP_COLORS = {"neutral_meson": ORANGE, "charged_hadron": BLUE, "photon": RED, "other": SLATE}
PARENT_COLORS = {
    "parton": "#7C3AED",
    "light meson resonance": "#CA8A04",
    "gamma parent": RED,
    "other parent": SLATE,
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
    txt(ax, 0.055, 0.935, claim, size=31, weight="bold")
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


def parent_category(pdg: int | float) -> str:
    try:
        ipdg = int(pdg)
    except (ValueError, TypeError):
        return "other parent"
    if ipdg in {1, 2, 3, 4, 5, 6, -1, -2, -3, -4, -5, -6, 21}:
        return "parton"
    if ipdg == 22:
        return "gamma parent"
    if abs(ipdg) in {113, 213, 223, 221, 331}:
        return "light meson resonance"
    return "other parent"


def load_data() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    cand = pd.read_csv(THE44_DIR / "the44_candidate_autopsy_rows.csv")
    part = pd.read_csv(THE44_DIR / "the44_particle_rows.csv")
    examples = pd.read_csv(THE44_DIR / "the44_event_card_examples.csv")
    inc = cand[(cand["label"] == "high-BDT background") & cand["sample"].str.contains("Jet", na=False)].copy()
    inc["source"] = inc["sample"].str.replace(" inclusive", "", regex=False)
    inc["truth_bucket"] = inc["truth_pid_name"].map(truth_bucket)
    inc["cone_bucket"] = np.select(
        [inc["final_state_photon_frac"] >= 0.80, inc["final_state_non_photon_frac"] >= 0.80],
        ["photon-rich", "non-photon-rich"],
        default="mixed",
    )
    return inc, part, examples


def pill(ax: plt.Axes, x: float, y: float, w: float, h: float, label: str, value: str, sub: str, color: str) -> None:
    card(ax, x, y, w, h, face=color, edge=color, radius=0.025)
    txt(ax, x + w / 2, y + h - 0.045, label, size=17, weight="bold", color="white", ha="center")
    txt(ax, x + w / 2, y + h * 0.52, value, size=32, weight="bold", color="white", ha="center", va="center")
    txt(ax, x + w / 2, y + 0.045, sub, size=14.5, color="white", ha="center", va="center", linespacing=1.0)


def stacked_hbar(
    ax: plt.Axes,
    x: float,
    y: float,
    w: float,
    h: float,
    values: dict[str, float],
    colors: dict[str, str],
    *,
    label_x: float | None = None,
    size: float = 13.5,
) -> None:
    total = sum(values.values()) or 1.0
    left = x
    for key, value in values.items():
        frac = value / total
        ax.add_patch(Rectangle((left, y), w * frac, h, transform=ax.transAxes, facecolor=colors[key], edgecolor="none", zorder=5))
        left += w * frac
    if label_x is not None:
        y0 = y - 0.030
        for i, (key, value) in enumerate(values.items()):
            yy = y0 - 0.032 * i
            ax.add_patch(Rectangle((label_x, yy - 0.010), 0.017, 0.017, transform=ax.transAxes, facecolor=colors[key], edgecolor="none", zorder=6))
            clean = key.replace("_", " ")
            txt(ax, label_x + 0.024, yy, f"{clean}: {100 * value / total:.0f}%", size=size, va="center")


def slide_01_coarse_flow(inc: pd.DataFrame, outdir: Path) -> Path:
    fig, ax = canvas()
    title(
        ax,
        "The high-BDT fake tail has a readable particle-production chain",
        "Candidate scores are joined to generator source, truth source, and local final-state particles.",
    )
    total = len(inc)
    jet3040 = int(inc["source"].isin(["Jet30", "Jet40"]).sum())
    meson = int(inc["truth_bucket"].eq("pi0/eta").sum())
    nonphoton = int(inc["cone_bucket"].eq("non-photon-rich").sum())
    med_score = float(inc["auau_tight_bdt_score"].median())
    steps = [
        ("Parent sample", f"{jet3040}/{total}", "Jet30/40 hard\ninclusive jets", GREEN),
        ("Truth source", f"{meson}/{total}", "pi0/eta starts\nthe fake cluster", ORANGE),
        ("Visible cone", f"{nonphoton}/{total}", "non-photon-rich\nstatus-1 particles", RED),
        ("BDT score", f"{med_score:.2f}", "still lands in the\nhigh-score tail", BLUE),
    ]
    xs = [0.070, 0.310, 0.550, 0.790]
    for i, (lab, val, sub, col) in enumerate(steps):
        pill(ax, xs[i], 0.52, 0.17, 0.22, lab, val, sub, col)
        if i < 3:
            ax.add_patch(
                FancyArrowPatch(
                    (xs[i] + 0.180, 0.63),
                    (xs[i + 1] - 0.012, 0.63),
                    transform=ax.transAxes,
                    arrowstyle="-|>",
                    mutation_scale=24,
                    linewidth=3,
                    color="#697386",
                    zorder=6,
                )
            )
    card(ax, 0.15, 0.255, 0.70, 0.17)
    txt(ax, 0.19, 0.382, "What this answers", size=20, weight="bold")
    txt(
        ax,
        0.19,
        0.335,
        "The confusing high-score region is not an undefined ML artifact. It is mostly hard-jet neutral-meson fragmentation that produces photon-like reconstructed clusters.",
        size=18,
        color=INK,
        linespacing=1.18,
    )
    card(ax, 0.17, 0.080, 0.66, 0.10, face=YELLOW, edge="#E5C65A")
    txt(ax, 0.50, 0.138, "This is the physical class Jet30/40 teach the BDT to push down.", size=21, weight="bold", ha="center", va="center")
    return finish(fig, outdir / "candidate01_coarse_fake_tail_chain.png")


def slide_02_source_truth_matrix(inc: pd.DataFrame, outdir: Path) -> Path:
    fig, ax = canvas()
    title(
        ax,
        "Jet30 and Jet40 do not add random fakes; they add pi0/eta fakes",
        "Rows are parent jet samples; colors are the truth source of high-BDT inclusive backgrounds.",
    )
    table = pd.crosstab(inc["source"], inc["truth_bucket"]).reindex(index=["Jet20", "Jet30", "Jet40"], columns=["pi0/eta", "gamma truth", "other"], fill_value=0)
    ax_bar = inset(fig, 0.095, 0.27, 0.56, 0.48)
    left = np.zeros(len(table))
    y = np.arange(len(table))
    for col in table.columns:
        vals = table[col].to_numpy()
        ax_bar.barh(y, vals, left=left, color=TRUTH_COLORS[col], height=0.66, label=col)
        for i, v in enumerate(vals):
            if v >= 6:
                ax_bar.text(left[i] + v / 2, i, f"{int(v)}", ha="center", va="center", fontsize=16, color="white", weight="bold")
        left += vals
    ax_bar.set_yticks(y)
    ax_bar.set_yticklabels(table.index, fontsize=17, fontweight="bold")
    ax_bar.invert_yaxis()
    ax_bar.set_xlabel("high-BDT inclusive candidates", fontsize=16)
    ax_bar.grid(True, axis="x", color=GRID)
    ax_bar.spines[["top", "right"]].set_visible(False)
    ax_bar.tick_params(axis="x", labelsize=13)

    card(ax, 0.710, 0.285, 0.225, 0.46)
    txt(ax, 0.740, 0.700, "Truth-source legend", size=19, weight="bold")
    for i, col in enumerate(table.columns):
        yy = 0.635 - 0.085 * i
        ax.add_patch(Rectangle((0.742, yy - 0.020), 0.026, 0.026, transform=ax.transAxes, facecolor=TRUTH_COLORS[col], edgecolor="none", zorder=5))
        txt(ax, 0.782, yy, col, size=17, weight="bold", va="center")
    txt(ax, 0.740, 0.405, "Key counts", size=19, weight="bold")
    txt(
        ax,
        0.740,
        0.358,
        f"Jet30 pi0/eta: {int(table.loc['Jet30','pi0/eta'])}\nJet40 pi0/eta: {int(table.loc['Jet40','pi0/eta'])}\nTotal pi0/eta: {int(table['pi0/eta'].sum())}/160",
        size=17,
        linespacing=1.15,
    )
    card(ax, 0.16, 0.080, 0.68, 0.10, face=YELLOW, edge="#E5C65A")
    txt(ax, 0.50, 0.138, "The dominant hard-background lesson is neutral-meson photon ID, not generic jet rejection.", size=20, weight="bold", ha="center", va="center")
    return finish(fig, outdir / "candidate02_source_truth_matrix.png")


def slide_03_status1_budget(inc: pd.DataFrame, part: pd.DataFrame, outdir: Path) -> Path:
    fig, ax = canvas()
    title(
        ax,
        "The local final-state pT budget is mostly mesons and charged fragments",
        "Status-1 Pythia particles within Delta R < 0.3 expose what actually surrounds the reconstructed cluster.",
    )
    fs = part.merge(inc[["candidate_uid"]], on="candidate_uid")
    fs = fs[fs["status"] == 1].copy()
    group_pt = fs.groupby("group")["pt"].sum()
    values = {
        "neutral_meson": group_pt.get("neutral_meson", 0.0),
        "charged_hadron": group_pt.get("charged_hadron", 0.0),
        "photon": group_pt.get("photon", 0.0),
        "other": max(group_pt.sum() - group_pt.reindex(["neutral_meson", "charged_hadron", "photon"]).fillna(0).sum(), 0.0),
    }
    card(ax, 0.085, 0.365, 0.45, 0.28)
    txt(ax, 0.125, 0.590, "Final-state pT composition", size=22, weight="bold")
    stacked_hbar(ax, 0.125, 0.505, 0.34, 0.060, values, GROUP_COLORS, label_x=0.125, size=15)

    top = fs.groupby(["pdg_name", "group"])["pt"].agg(["sum", "count"]).sort_values("sum", ascending=False).head(8).reset_index()
    ax_top = inset(fig, 0.610, 0.31, 0.30, 0.40)
    y = np.arange(len(top))
    colors = [GROUP_COLORS.get(g, SLATE) for g in top["group"]]
    ax_top.barh(y, top["sum"], color=colors)
    ax_top.set_yticks(y)
    ax_top.set_yticklabels(top["pdg_name"], fontsize=14)
    ax_top.invert_yaxis()
    ax_top.set_xlabel("summed status-1 pT", fontsize=14)
    ax_top.grid(True, axis="x", color=GRID)
    ax_top.spines[["top", "right"]].set_visible(False)
    ax_top.tick_params(axis="x", labelsize=12)
    txt(ax, 0.610, 0.755, "Dominant status-1 particles", size=21, weight="bold")

    card(ax, 0.15, 0.080, 0.70, 0.10, face=YELLOW, edge="#E5C65A")
    txt(ax, 0.50, 0.138, "At the particle-list level, the high-BDT fake tail is meson-rich, not photon-core dominated.", size=20, weight="bold", ha="center", va="center")
    return finish(fig, outdir / "candidate03_status1_particle_budget.png")


def slide_04_parentage(inc: pd.DataFrame, part: pd.DataFrame, outdir: Path) -> Path:
    fig, ax = canvas()
    title(
        ax,
        "The fine-grained HepMC view shows fragmentation ancestry",
        "Final-state particles are fed by partons and light-meson resonances before reaching the local cone.",
    )
    fs = part.merge(inc[["candidate_uid"]], on="candidate_uid")
    fs = fs[fs["status"] == 1].copy()
    fs["parent_cat"] = fs["parent0_pdg"].map(parent_category)
    parent_pt = fs.groupby("parent_cat")["pt"].sum().reindex(["parton", "light meson resonance", "gamma parent", "other parent"]).fillna(0)
    ax_parent = inset(fig, 0.09, 0.31, 0.36, 0.43)
    y = np.arange(len(parent_pt))
    ax_parent.barh(y, parent_pt.values, color=[PARENT_COLORS[k] for k in parent_pt.index])
    ax_parent.set_yticks(y)
    ax_parent.set_yticklabels(parent_pt.index, fontsize=14)
    ax_parent.invert_yaxis()
    ax_parent.set_xlabel("summed status-1 daughter pT", fontsize=14)
    ax_parent.grid(True, axis="x", color=GRID)
    ax_parent.spines[["top", "right"]].set_visible(False)
    txt(ax, 0.09, 0.775, "Parent category feeding the cone", size=21, weight="bold")

    group_parent = fs.groupby(["parent_cat", "group"])["pt"].sum().unstack(fill_value=0)
    card(ax, 0.545, 0.285, 0.38, 0.46)
    txt(ax, 0.575, 0.705, "How to read the ancestry", size=21, weight="bold")
    txt(
        ax,
        0.575,
        0.645,
        "Partons and light-meson resonances feed\nstatus-1 pi0/eta and charged fragments.\nThose fragments sit near the cluster and make\nthe reconstructed object photon-like enough\nto challenge the BDT.",
        size=17.2,
        linespacing=1.15,
    )
    rows = [
        ("parton -> neutral meson", group_parent.get("neutral_meson", pd.Series()).get("parton", 0.0), PURPLE),
        ("resonance -> neutral meson", group_parent.get("neutral_meson", pd.Series()).get("light meson resonance", 0.0), ORANGE),
        ("parton -> charged hadron", group_parent.get("charged_hadron", pd.Series()).get("parton", 0.0), BLUE),
    ]
    for i, (lab, val, col) in enumerate(rows):
        yy = 0.405 - 0.055 * i
        ax.add_patch(Rectangle((0.575, yy - 0.016), 0.020, 0.020, transform=ax.transAxes, facecolor=col, edgecolor="none", zorder=5))
        txt(ax, 0.607, yy, f"{lab}: {val:.0f} pT", size=15.5, va="center")

    card(ax, 0.15, 0.080, 0.70, 0.10, face=YELLOW, edge="#E5C65A")
    txt(ax, 0.50, 0.138, "The background is fragmentation history made visible in the photon-candidate cone.", size=20, weight="bold", ha="center", va="center")
    return finish(fig, outdir / "candidate04_hepmc_parentage_context.png")


def slide_05_detector_tension(inc: pd.DataFrame, examples: pd.DataFrame, outdir: Path) -> Path:
    fig, ax = canvas()
    title(
        ax,
        "Why the BDT can give a high score to a non-photon cone",
        "The same fake can have EM-like shower variables while the particle list says fragmentation.",
    )
    fake = examples[examples["example_role"].eq("highest-score non-photon background")].sort_values("auau_tight_bdt_score", ascending=False).iloc[0]
    true = examples[examples["example_role"].eq("busy true photon with lowered score")].sort_values("auau_tight_bdt_score", ascending=False).iloc[0]
    metrics = [
        ("BDT score", "auau_tight_bdt_score", 1.0, True),
        ("e11/e33", "e11_over_e33", 1.0, True),
        ("e32/e35", "e32_over_e35", 1.0, True),
        ("reco iso / 12", "reco_eiso", 12.0, True),
        ("non-photon frac", "final_state_non_photon_frac", 1.0, True),
    ]
    x0, y0 = 0.12, 0.660
    fake_x, true_x, bar_w = 0.295, 0.665, 0.285
    txt(ax, 0.12, 0.765, "Candidate-level variables", size=21, weight="bold")
    for i, (label, key, scale, _) in enumerate(metrics):
        yy = y0 - 0.090 * i
        txt(ax, x0, yy + 0.017, label, size=15.5, va="center")
        ax.add_patch(Rectangle((fake_x, yy), bar_w, 0.035, transform=ax.transAxes, facecolor="#E8EDF3", edgecolor="none", zorder=3))
        ax.add_patch(Rectangle((fake_x, yy), bar_w * min(float(fake[key]) / scale, 1.0), 0.035, transform=ax.transAxes, facecolor=RED, edgecolor="none", zorder=4))
        ax.add_patch(Rectangle((true_x, yy), bar_w, 0.035, transform=ax.transAxes, facecolor="#E8EDF3", edgecolor="none", zorder=3))
        ax.add_patch(Rectangle((true_x, yy), bar_w * min(float(true[key]) / scale, 1.0), 0.035, transform=ax.transAxes, facecolor=GREEN, edgecolor="none", zorder=4))
        txt(ax, fake_x + bar_w + 0.011, yy + 0.017, f"{float(fake[key]):.2f}", size=13.5, ha="left", va="center", color=RED, weight="bold")
        txt(ax, true_x + bar_w + 0.006, yy + 0.017, f"{float(true[key]):.2f}", size=13.5, ha="left", va="center", color=GREEN, weight="bold")
    txt(ax, fake_x + bar_w / 2, 0.713, "high-BDT fake", size=16, weight="bold", color=RED, ha="center")
    txt(ax, true_x + bar_w / 2, 0.713, "truth-photon control", size=16, weight="bold", color=GREEN, ha="center")

    card(ax, 0.15, 0.110, 0.70, 0.16)
    txt(ax, 0.185, 0.220, "Interpretation", size=20, weight="bold")
    txt(
        ax,
        0.185,
        0.175,
        "The fake is not absurd to the BDT: its compact EM-like variables can look signal-like.\nThe particle list adds the missing truth: the local cone is not photon-core dominated.",
        size=17.5,
        linespacing=1.15,
    )
    return finish(fig, outdir / "candidate05_detector_vs_particle_tension.png")


def slide_06_clean_microscope(examples: pd.DataFrame, outdir: Path) -> Path:
    fig, ax = canvas()
    title(
        ax,
        "The particle list turns one BDT score into an event-level diagnosis",
        "Two individual candidates show the distinction between photon core and fragmentation cone.",
    )
    fake = examples[examples["example_role"].eq("highest-score non-photon background")].sort_values("auau_tight_bdt_score", ascending=False).iloc[0]
    true = examples[examples["example_role"].eq("busy true photon with lowered score")].sort_values("auau_tight_bdt_score", ascending=False).iloc[0]

    def draw_case(x: float, row: pd.Series, header: str, color: str) -> None:
        card(ax, x, 0.225, 0.40, 0.55, face=CARD, edge=color, lw=1.8)
        txt(ax, x + 0.030, 0.725, header, size=22, color=color, weight="bold")
        txt(ax, x + 0.030, 0.665, f"{row['sample']}   truth {row['truth_pid_name']}", size=16, weight="bold")
        txt(ax, x + 0.030, 0.595, f"BDT {row['auau_tight_bdt_score']:.3f}", size=30, color=color, weight="bold")
        txt(ax, x + 0.250, 0.605, f"cluster E$_T$ {row['cluster_Et']:.1f} GeV\nreco iso {row['reco_eiso']:.1f} GeV", size=15.5, linespacing=1.15)
        vals = {
            "photon": float(row["final_state_photon_frac"]),
            "neutral_meson": float(row["final_neutral_meson_frac"]),
            "charged_hadron": float(row["final_charged_hadron_frac"]),
        }
        txt(ax, x + 0.030, 0.515, "Final-state cone", size=18, weight="bold")
        for i, (key, val) in enumerate(vals.items()):
            yy = 0.445 - 0.065 * i
            ax.add_patch(Rectangle((x + 0.030, yy), 0.245, 0.035, transform=ax.transAxes, facecolor="#E8EDF3", edgecolor="none", zorder=4))
            ax.add_patch(Rectangle((x + 0.030, yy), 0.245 * val, 0.035, transform=ax.transAxes, facecolor=GROUP_COLORS[key], edgecolor="none", zorder=5))
            txt(ax, x + 0.300, yy + 0.018, f"{key.replace('_', ' ')} {100 * val:.0f}%", size=14.5, va="center")
        tokens = str(row["top_particles"]).split()
        top = " ".join(tokens[:9]) + (" ..." if len(tokens) > 9 else "")
        txt(ax, x + 0.030, 0.265, f"Local list: {top}", size=15, color=MUTED)

    draw_case(0.075, fake, "High-BDT inclusive fake", RED)
    draw_case(0.525, true, "Truth-photon control", GREEN)
    card(ax, 0.18, 0.085, 0.64, 0.10, face=YELLOW, edge="#E5C65A")
    txt(ax, 0.50, 0.143, "The same score scale can now be interpreted as a physical local particle fingerprint.", size=20, weight="bold", ha="center", va="center")
    return finish(fig, outdir / "candidate06_clean_candidate_microscope.png")


def contact_sheet(paths: list[Path], outdir: Path) -> Path:
    fig = plt.figure(figsize=(16, 9), dpi=160)
    fig.patch.set_facecolor(LIGHT)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    txt(ax, 0.045, 0.945, "THE-44 Particle-Story Candidate Set", size=30, weight="bold")
    txt(ax, 0.045, 0.895, "Six alternate slide candidates, ordered from coarse chain to fine-grained candidate diagnosis.", size=16.5, color=MUTED)
    labels = [
        "1. Coarse fake-tail chain",
        "2. Source-to-truth matrix",
        "3. Status-1 particle budget",
        "4. HepMC parentage context",
        "5. Detector vs particle tension",
        "6. Clean candidate microscope",
    ]
    for i, (path, label) in enumerate(zip(paths, labels)):
        col = i % 3
        row = i // 3
        x = 0.045 + col * 0.315
        y = 0.515 - row * 0.390
        txt(ax, x, y + 0.295, label, size=14.5, weight="bold")
        img = mpimg.imread(path)
        iax = fig.add_axes([x, y, 0.285, 0.255])
        iax.imshow(img)
        iax.set_axis_off()
        iax.add_patch(Rectangle((0, 0), 1, 1, transform=iax.transAxes, facecolor="none", edgecolor="#C8D0DA", linewidth=2))
    return finish(fig, outdir / "the44_particle_story_candidate_contact_sheet.png")


def write_scripts(outdir: Path, metrics: dict[str, float]) -> dict[str, str]:
    scripts = {
        "candidate01_coarse_fake_tail_chain_script.md": f"""# Candidate 1 Script - Coarse Fake-Tail Chain

This slide is the one-sentence physical diagnosis of the high-BDT fake tail. Starting from 160 high-BDT inclusive backgrounds, 146 come from Jet30 or Jet40, 126 have pi0 or eta as the truth source, and 138 are non-photon-rich in the final-state cone.

The point is that this BDT region has a real particle-production identity. It is mostly hard neutral-meson fragmentation, not an undefined classifier artifact.
""",
        "candidate02_source_truth_matrix_script.md": """# Candidate 2 Script - Source-To-Truth Matrix

This slide breaks the same story down by parent jet sample. Jet30 and Jet40 are the main suppliers of the high-BDT inclusive tail, and their dominant truth source is pi0 or eta.

That means the extra samples are teaching the BDT the hard neutral-meson background class that was underrepresented in the lower-energy training.
""",
        "candidate03_status1_particle_budget_script.md": """# Candidate 3 Script - Status-1 Particle Budget

This slide goes one layer deeper and looks only at status-1 particles within Delta R less than 0.3 around the cluster. The summed pT is dominated by neutral mesons and charged hadrons, with a much smaller photon component.

So at the final-state particle level, these high-BDT fakes are not photon-core objects. They are fragmentation cones that reconstruct in a photon-like way.
""",
        "candidate04_hepmc_parentage_context_script.md": """# Candidate 4 Script - HepMC Parentage Context

This slide uses the parent information in the particle list. The local cone is fed mostly by parton fragmentation and light-meson resonance chains before producing the status-1 particles near the cluster.

That ancestry is the deep physical context: the background is a fragmentation process that ends in mesons and charged fragments close enough to the EM cluster to challenge photon ID.
""",
        "candidate05_detector_vs_particle_tension_script.md": """# Candidate 5 Script - Detector Versus Particle Tension

This slide explains why the BDT can assign a high score to something that is not a photon core. The fake can have compact EM-like shower variables, so it is not nonsensical that the classifier finds it signal-like.

The particle list adds the missing information: despite the EM-like detector variables, the local final-state cone is dominated by non-photon particles.
""",
        "candidate06_clean_candidate_microscope_script.md": """# Candidate 6 Script - Clean Candidate Microscope

This is the most concrete event-level version. The left candidate is a high-BDT inclusive fake with a pi0 truth source and a non-photon final-state cone. The right candidate is a truth-photon control with a photon-dominated core and some nearby activity.

The main message is that THE-44 lets us interpret the BDT score as a physical local particle fingerprint, candidate by candidate.
""",
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

    inc, part, examples = load_data()
    paths = [
        slide_01_coarse_flow(inc, args.outdir),
        slide_02_source_truth_matrix(inc, args.outdir),
        slide_03_status1_budget(inc, part, args.outdir),
        slide_04_parentage(inc, part, args.outdir),
        slide_05_detector_tension(inc, examples, args.outdir),
        slide_06_clean_microscope(examples, args.outdir),
    ]
    contact = contact_sheet(paths, args.outdir)
    metrics = {
        "high_bdt_inclusive_candidates": float(len(inc)),
        "jet30_40_candidates": float(inc["source"].isin(["Jet30", "Jet40"]).sum()),
        "pi0_eta_truth_candidates": float(inc["truth_bucket"].eq("pi0/eta").sum()),
        "non_photon_rich_cones": float(inc["cone_bucket"].eq("non-photon-rich").sum()),
        "median_high_fake_bdt_score": float(inc["auau_tight_bdt_score"].median()),
        "median_high_fake_photon_fraction": float(inc["final_state_photon_frac"].median()),
        "median_high_fake_neutral_meson_fraction": float(inc["final_neutral_meson_frac"].median()),
    }
    scripts = write_scripts(args.outdir, metrics)
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
        "caveat": "Internal qualitative diagnostic from THE-44 saved candidate/particle lists; not a production-weighted purity estimate.",
    }
    manifest_path = args.outdir / "the44_particle_story_candidate_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
