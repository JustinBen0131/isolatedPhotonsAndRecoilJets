#!/usr/bin/env python3
"""Build the THE-44 particle-list physics interpretation slide sequence."""

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
MECH_DIR = BASE / "diagnostics/hard_inclusive_jet3040_mechanism_20260605"
FIXED_DIR = BASE / "fixed_sample_controls"
DEFAULT_OUTDIR = THE44_DIR / "particle_intensive/deep_physics_sequence_20260606"

INK = "#111827"
MUTED = "#5B6472"
LIGHT = "#F6F7F9"
CARD = "#FFFFFF"
GRID = "#E0E5EC"
BLUE = "#2563A6"
RED = "#C43C32"
GREEN = "#2F8F52"
ORANGE = "#D97706"
PURPLE = "#7C5EA6"
TEAL = "#0F766E"
YELLOW = "#FFF4C7"
SLATE = "#475569"

SOURCE_COLORS = {
    "Jet12": "#8AA6C2",
    "Jet20": "#F97316",
    "Jet30": "#2F9E44",
    "Jet40": "#7651A6",
}
PARTICLE_COLORS = {
    "photon": RED,
    "neutral_meson": ORANGE,
    "charged_hadron": BLUE,
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
        zorder=10,
    )


def card(
    ax: plt.Axes,
    x: float,
    y: float,
    w: float,
    h: float,
    *,
    face: str = CARD,
    edge: str = "#D5DCE6",
    lw: float = 1.2,
    radius: float = 0.018,
    zorder: int = 1,
) -> FancyBboxPatch:
    patch = FancyBboxPatch(
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
    ax.add_patch(patch)
    return patch


def title(ax: plt.Axes, claim: str, subtitle: str) -> None:
    txt(ax, 0.055, 0.935, claim, size=32, weight="bold")
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


def load_inputs() -> dict[str, pd.DataFrame]:
    candidate = pd.read_csv(THE44_DIR / "the44_candidate_autopsy_rows.csv")
    particle = pd.read_csv(THE44_DIR / "the44_particle_rows.csv")
    event_cards = pd.read_csv(THE44_DIR / "the44_event_card_examples.csv")
    source = pd.read_csv(MECH_DIR / "source_fractions_holdout_cent020.csv")
    fixed_auc = pd.read_csv(MECH_DIR / "fixed_validation_training_effect_auc.csv")
    median_gap = pd.read_csv(
        FIXED_DIR
        / "the8_branchA_fixed_sample_holdout_3x3_direct_0to20_truthSignal_inclusiveJet_slide_v34_mediangap_ppbasev3e_clean.median_gap_summary.csv"
    )
    return {
        "candidate": candidate,
        "particle": particle,
        "event_cards": event_cards,
        "source": source,
        "fixed_auc": fixed_auc,
        "median_gap": median_gap,
    }


def high_inclusive(cand: pd.DataFrame) -> pd.DataFrame:
    inc = cand[(cand["label"] == "high-BDT background") & cand["sample"].str.contains("Jet", na=False)].copy()
    inc["source"] = inc["sample"].str.replace(" inclusive", "", regex=False)
    inc["truth_bucket"] = inc["truth_pid_name"].map(truth_bucket)
    inc["cone_bucket"] = np.select(
        [inc["final_state_photon_frac"] >= 0.80, inc["final_state_non_photon_frac"] >= 0.80],
        ["photon-rich", "non-photon-rich"],
        default="mixed",
    )
    return inc


def slide_score_metric(data: dict[str, pd.DataFrame], outdir: Path) -> Path:
    fig, ax = canvas()
    title(
        ax,
        "The visible BDT gain is score separation, not just AUC",
        "Holding the Jet12+20 validation row fixed, AUC barely moves while the median signal-background score gap opens.",
    )

    fixed = data["fixed_auc"]
    auc = (
        fixed[(fixed["scope"] == "direct_holdout_3x3") & (fixed["validation_sample"] == "Jet12+20")]
        .set_index("model_training_sample")
        .loc[["Jet12+20", "Jet12+20+30", "Jet12+20+30+40"], "auc_0_20"]
        .to_numpy()
    )
    gap = (
        data["median_gap"][data["median_gap"]["validation_sample"] == "Jet12+20"]
        .set_index("training_sample")
        .loc[["Jet12+20", "Jet12+20+30", "Jet12+20+30+40"], "median_bdt_score_gap_binned"]
        .to_numpy()
    )
    labels = ["Jet12+20", "+Jet30", "+Jet40"]
    x = np.arange(3)

    card(ax, 0.065, 0.20, 0.42, 0.62)
    txt(ax, 0.095, 0.775, "Global ranking: small movement", size=21, weight="bold")
    ax_auc = inset(fig, 0.115, 0.34, 0.30, 0.34)
    ax_auc.plot(x, auc, marker="o", ms=9, lw=4, color=BLUE)
    ax_auc.set_xticks(x)
    ax_auc.set_xticklabels(labels, fontsize=13)
    ax_auc.set_ylim(0.826, 0.838)
    ax_auc.set_ylabel("AUC, 0-20%", fontsize=15)
    ax_auc.grid(True, color=GRID)
    ax_auc.spines[["top", "right"]].set_visible(False)
    ax_auc.tick_params(labelsize=12)
    ax_auc.text(2, auc[-1] + 0.0007, f"+{auc[-1]-auc[0]:.3f}", ha="center", fontsize=17, weight="bold", color=BLUE)
    txt(ax, 0.095, 0.285, "AUC says the change is modest.", size=18, color=BLUE, weight="bold")

    card(ax, 0.535, 0.20, 0.40, 0.62)
    txt(ax, 0.565, 0.775, "Score polarization: visible movement", size=21, weight="bold")
    ax_gap = inset(fig, 0.595, 0.34, 0.27, 0.34)
    ax_gap.bar(x, gap, color=[BLUE, ORANGE, GREEN], width=0.68)
    ax_gap.set_xticks(x)
    ax_gap.set_xticklabels(labels, fontsize=13)
    ax_gap.set_ylim(0, 0.50)
    ax_gap.set_ylabel("median BDT gap", fontsize=15)
    ax_gap.grid(True, axis="y", color=GRID)
    ax_gap.spines[["top", "right"]].set_visible(False)
    ax_gap.tick_params(labelsize=12)
    for i, value in enumerate(gap):
        ax_gap.text(i, value + 0.018, f"{value:.2f}", ha="center", fontsize=19, weight="bold")
    txt(ax, 0.565, 0.285, "The typical photon moves farther from the typical jet.", size=18, color=GREEN, weight="bold")

    card(ax, 0.155, 0.075, 0.69, 0.10, face=YELLOW, edge="#E6C65D")
    txt(
        ax,
        0.50,
        0.132,
        "This is the clean metric for the qualitative top-row split:\nthe BDT score scale becomes more physically separated.",
        size=18.6,
        weight="bold",
        ha="center",
        va="center",
        linespacing=1.05,
    )
    return finish(fig, outdir / "01_score_gap_not_auc.png")


def slide_energy_coverage(data: dict[str, pd.DataFrame], inc: pd.DataFrame, outdir: Path) -> Path:
    fig, ax = canvas()
    title(
        ax,
        "Jet30/40 add the hard inclusive backgrounds that define the boundary",
        "The expanded training sample covers the parent-jet energy range that feeds the confusing photon-candidate E$_T$ region.",
    )

    source = data["source"]
    src = source[
        (source["branch"] == "Jet12+20+30+40")
        & (source["scope"] == "holdout10")
        & (source["centrality_bin"] == "0-20%")
    ]
    et_bins = ["15-20", "20-25", "25-30", "30-35"]
    x = np.arange(len(et_bins))
    ax_bar = inset(fig, 0.08, 0.24, 0.61, 0.51)
    bottom = np.zeros(len(et_bins))
    for sample in ["Jet12", "Jet20", "Jet30", "Jet40"]:
        vals = []
        for et in et_bins:
            row = src[(src["et_bin"] == et) & (src["source_bin"] == sample)]
            vals.append(float(row["source_fraction_weighted"].iloc[0]) if len(row) else 0.0)
        vals = np.asarray(vals) * 100.0
        ax_bar.bar(x, vals, bottom=bottom, color=SOURCE_COLORS[sample], width=0.72, label=sample)
        for i, value in enumerate(vals):
            if value >= 12:
                ax_bar.text(i, bottom[i] + value / 2, f"{value:.0f}%", ha="center", va="center", fontsize=14, color="white", weight="bold")
        bottom += vals
    ax_bar.set_ylim(0, 100)
    ax_bar.set_xticks(x)
    ax_bar.set_xticklabels(et_bins, fontsize=15)
    ax_bar.set_xlabel("candidate cluster E$_T$ [GeV]", fontsize=16)
    ax_bar.set_ylabel("weighted inclusive-background share", fontsize=16)
    ax_bar.grid(True, axis="y", color=GRID)
    ax_bar.spines[["top", "right"]].set_visible(False)
    ax_bar.tick_params(labelsize=13)

    card(ax, 0.735, 0.205, 0.20, 0.55)
    txt(ax, 0.765, 0.705, "Generator bins", size=20, weight="bold")
    rows = [("Jet20", "pThat 21-32 GeV"), ("Jet30", "pThat 32-42 GeV"), ("Jet40", "pThat 42-100 GeV")]
    for i, (name, desc) in enumerate(rows):
        y = 0.63 - 0.095 * i
        ax.add_patch(Rectangle((0.765, y - 0.018), 0.023, 0.023, transform=ax.transAxes, facecolor=SOURCE_COLORS[name], edgecolor="none", zorder=5))
        txt(ax, 0.803, y, name, size=17, weight="bold", va="center")
        txt(ax, 0.803, y - 0.030, desc, size=12.5, color=MUTED, va="center")
    med = inc.groupby("source")["cluster_Et"].median()
    txt(ax, 0.765, 0.365, "High-BDT fake\ncluster-E$_T$ medians:", size=17, weight="bold", linespacing=1.05)
    txt(
        ax,
        0.765,
        0.285,
        f"Jet20 {med.get('Jet20', np.nan):.1f} GeV\nJet30 {med.get('Jet30', np.nan):.1f} GeV\nJet40 {med.get('Jet40', np.nan):.1f} GeV",
        size=15.5,
        color=INK,
        linespacing=1.15,
    )

    card(ax, 0.12, 0.075, 0.76, 0.10, face=YELLOW, edge="#E6C65D")
    txt(
        ax,
        0.50,
        0.132,
        "Jet40 is not decorative extra data:\nit supplies the high-E$_T$ inclusive population the BDT must learn to reject.",
        size=18.5,
        weight="bold",
        ha="center",
        va="center",
        linespacing=1.05,
    )
    return finish(fig, outdir / "02_jet40_supplies_hard_background.png")


def slide_particle_pipeline(data: dict[str, pd.DataFrame], inc: pd.DataFrame, outdir: Path) -> Path:
    fig, ax = canvas()
    title(
        ax,
        "The hard high-BDT background is neutral-meson fragmentation",
        "THE-44 joins each candidate score to the local Pythia/HepMC particle list within Delta R < 0.3.",
    )

    total = len(inc)
    jet3040 = int(inc["source"].isin(["Jet30", "Jet40"]).sum())
    meson = int(inc["truth_bucket"].eq("pi0/eta").sum())
    nonphoton = int(inc["cone_bucket"].eq("non-photon-rich").sum())
    med_score = float(inc["auau_tight_bdt_score"].median())
    med_photon = float(inc["final_state_photon_frac"].median())

    steps = [
        ("Parent bin", f"{jet3040}/{total}", "Jet30/40\nhard generator bins", GREEN),
        ("Truth source", f"{meson}/{total}", "pi0/eta\ncandidate source", ORANGE),
        ("Visible cone", f"{nonphoton}/{total}", "non-photon-rich\nfinal state", RED),
        ("BDT problem", f"{med_score:.2f}", f"median score\nmedian photon frac {100 * med_photon:.0f}%", BLUE),
    ]
    xs = [0.07, 0.31, 0.55, 0.79]
    for i, (head, big, sub, color) in enumerate(steps):
        card(ax, xs[i], 0.51, 0.17, 0.22, face=color, edge=color, radius=0.025)
        txt(ax, xs[i] + 0.085, 0.685, head, size=18, weight="bold", color="white", ha="center")
        txt(ax, xs[i] + 0.085, 0.605, big, size=34, weight="bold", color="white", ha="center", va="center")
        txt(ax, xs[i] + 0.085, 0.545, sub, size=15.5, color="white", ha="center", va="center", linespacing=1.0)
        if i < len(steps) - 1:
            ax.add_patch(
                FancyArrowPatch(
                    (xs[i] + 0.18, 0.62),
                    (xs[i + 1] - 0.015, 0.62),
                    transform=ax.transAxes,
                    arrowstyle="-|>",
                    mutation_scale=24,
                    linewidth=3,
                    color="#697386",
                    zorder=6,
                )
            )

    particle = data["particle"]
    fs = particle.merge(inc[["candidate_uid"]], on="candidate_uid")
    fs = fs[fs["status"] == 1]
    group_pt = fs.groupby("group")["pt"].sum()
    shown = {
        "neutral_meson": group_pt.get("neutral_meson", 0.0),
        "charged_hadron": group_pt.get("charged_hadron", 0.0),
        "photon": group_pt.get("photon", 0.0),
    }
    shown["other"] = max(group_pt.sum() - sum(shown.values()), 0.0)
    total_pt = sum(shown.values())

    card(ax, 0.105, 0.22, 0.34, 0.22)
    txt(ax, 0.135, 0.395, "Final-state pT composition", size=20, weight="bold")
    x0, y0, width = 0.135, 0.315, 0.27
    left = x0
    for group in ["neutral_meson", "charged_hadron", "photon", "other"]:
        frac = shown[group] / total_pt if total_pt else 0.0
        ax.add_patch(Rectangle((left, y0), width * frac, 0.045, transform=ax.transAxes, facecolor=PARTICLE_COLORS[group], edgecolor="none", zorder=4))
        left += width * frac
    legend_rows = [("neutral meson", "neutral_meson"), ("charged hadron", "charged_hadron"), ("photon", "photon")]
    for i, (label, group) in enumerate(legend_rows):
        y = 0.280 - 0.035 * i
        ax.add_patch(Rectangle((0.135, y - 0.014), 0.018, 0.018, transform=ax.transAxes, facecolor=PARTICLE_COLORS[group], edgecolor="none", zorder=5))
        txt(ax, 0.163, y, f"{label}: {100 * shown[group] / total_pt:.0f}%", size=14, va="center")

    card(ax, 0.515, 0.22, 0.38, 0.22)
    txt(ax, 0.545, 0.395, "Interpretation", size=20, weight="bold")
    txt(
        ax,
        0.545,
        0.345,
        "These are photon-like reconstructed clusters,\nbut their local particle neighborhoods are\nmostly mesons and charged fragments.",
        size=17,
        color=INK,
        linespacing=1.12,
    )
    card(ax, 0.13, 0.075, 0.74, 0.10, face=YELLOW, edge="#E6C65D")
    txt(
        ax,
        0.50,
        0.132,
        "The expanded sample teaches the BDT the physical fake factory:\nhard pi0/eta fragmentation near the photon-candidate E$_T$ range.",
        size=18.0,
        weight="bold",
        ha="center",
        va="center",
        linespacing=1.05,
    )
    return finish(fig, outdir / "03_particle_pipeline_fake_factory.png")


def slide_candidate_microscope(data: dict[str, pd.DataFrame], inc: pd.DataFrame, outdir: Path) -> Path:
    fig, ax = canvas()
    title(
        ax,
        "At candidate level, the BDT is sorting local particle fingerprints",
        "A high-score fake can look EM-like in shower variables while its generator neighborhood is not photon-core dominated.",
    )

    cards = data["event_cards"]
    fake = cards[cards["example_role"].eq("highest-score non-photon background")].sort_values("auau_tight_bdt_score", ascending=False).iloc[0]
    true = cards[cards["example_role"].eq("busy true photon with lowered score")].sort_values("auau_tight_bdt_score", ascending=False).iloc[0]

    def draw_case(x: float, row: pd.Series, title_text: str, color: str) -> None:
        card(ax, x, 0.22, 0.40, 0.58, face=CARD, edge=color, lw=1.8)
        txt(ax, x + 0.025, 0.755, title_text, size=22, weight="bold", color=color)
        txt(ax, x + 0.025, 0.700, f"{row['sample']}   truth {row['truth_pid_name']}", size=16, weight="bold")
        txt(ax, x + 0.025, 0.642, f"BDT score {row['auau_tight_bdt_score']:.3f}", size=26, weight="bold", color=color)
        txt(ax, x + 0.245, 0.653, f"cluster E$_T$ {row['cluster_Et']:.1f} GeV\nreco iso {row['reco_eiso']:.1f} GeV", size=15, linespacing=1.15)
        vals = {
            "photon": float(row["final_state_photon_frac"]),
            "neutral_meson": float(row["final_neutral_meson_frac"]),
            "charged_hadron": float(row["final_charged_hadron_frac"]),
        }
        y0 = 0.485
        txt(ax, x + 0.025, y0 + 0.095, "Final-state cone composition", size=17, weight="bold")
        for i, (name, value) in enumerate(vals.items()):
            y = y0 - 0.060 * i
            ax.add_patch(Rectangle((x + 0.025, y), 0.27, 0.033, transform=ax.transAxes, facecolor="#E8EDF3", edgecolor="none", zorder=3))
            ax.add_patch(Rectangle((x + 0.025, y), 0.27 * value, 0.033, transform=ax.transAxes, facecolor=PARTICLE_COLORS[name], edgecolor="none", zorder=4))
            label = name.replace("_", " ")
            txt(ax, x + 0.315, y + 0.017, f"{label}: {100*value:.0f}%", size=14.5, va="center")
        txt(ax, x + 0.025, 0.260, f"Top local particles: {row['top_particles']}", size=15, color=MUTED)

    draw_case(0.075, fake, "High-BDT inclusive fake", RED)
    draw_case(0.525, true, "Truth-photon control", GREEN)

    card(ax, 0.13, 0.075, 0.74, 0.10, face=YELLOW, edge="#E6C65D")
    txt(
        ax,
        0.50,
        0.132,
        "This is why particle lists matter:\nthey show whether a score feature is a photon core or a fragmentation cone.",
        size=18.8,
        weight="bold",
        ha="center",
        va="center",
        linespacing=1.05,
    )
    return finish(fig, outdir / "04_candidate_level_particle_microscope.png")


def contact_sheet(paths: list[Path], outdir: Path) -> Path:
    fig = plt.figure(figsize=(16, 9), dpi=160)
    fig.patch.set_facecolor(LIGHT)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    txt(ax, 0.05, 0.945, "THE-44 Deep BDT Physics Sequence", size=30, weight="bold")
    txt(ax, 0.05, 0.895, "Use as four separate slides; this page is only the review contact sheet.", size=17, color=MUTED)
    labels = [
        "1. Score split beats AUC",
        "2. Jet40 supplies hard background",
        "3. Particle pipeline",
        "4. Candidate microscope",
    ]
    positions = [(0.045, 0.49), (0.515, 0.49), (0.045, 0.09), (0.515, 0.09)]
    for path, label, (x, y) in zip(paths, labels, positions):
        img = mpimg.imread(path)
        iax = fig.add_axes([x, y, 0.43, 0.305])
        iax.imshow(img)
        iax.set_axis_off()
        iax.add_patch(Rectangle((0, 0), 1, 1, transform=iax.transAxes, facecolor="none", edgecolor="#C9D0DA", linewidth=2))
        txt(ax, x, y + 0.335, label, size=16.5, weight="bold")
    return finish(fig, outdir / "the44_deep_physics_sequence_contact_sheet.png")


def write_scripts(outdir: Path, metrics: dict[str, float]) -> dict[str, str]:
    scripts: dict[str, str] = {}
    texts = {
        "01_score_gap_not_auc_script.md": f"""# THE-44 Slide 1 Script - Score Split Beats AUC

First I want to make the metric match what we see by eye. On the fixed Jet12 plus Jet20 validation row, the AUC only moves from {metrics['auc_start']:.3f} to {metrics['auc_end']:.3f}. That is a real change, but it is small because AUC averages over many easy signal-background pairs.

The more direct way to describe the visual separation is the median BDT score gap. That grows from {metrics['gap_start']:.2f} to {metrics['gap_end']:.2f}. In simple terms, the typical truth photon is being placed farther away from the typical inclusive-jet candidate in score space. That is the qualitative split we were trying to quantify.
""",
        "02_jet40_supplies_hard_background_script.md": """# THE-44 Slide 2 Script - Jet40 Supplies The Hard Background

The next question is why the expanded training sample changes that score split. This slide shows the source composition of the inclusive background in 0 to 20 percent centrality, broken out by candidate cluster energy.

The key feature is the Jet40 component. In the upper cluster-energy bins, Jet40 becomes the dominant source of inclusive background. So the extra samples are not just generic statistics. They add the hard parent-jet region that feeds the confusing photon-candidate energy range.
""",
        "03_particle_pipeline_fake_factory_script.md": f"""# THE-44 Slide 3 Script - Particle Pipeline Fake Factory

Now the particle-list autopsy tells us what that hard background physically is. Looking only at the high-BDT inclusive tail, {metrics['jet3040_count']:.0f} out of {metrics['high_inc_count']:.0f} candidates come from Jet30 or Jet40. At the truth-source level, {metrics['meson_count']:.0f} out of {metrics['high_inc_count']:.0f} are pi0 or eta. At the final-state cone level, {metrics['nonphoton_count']:.0f} out of {metrics['high_inc_count']:.0f} are non-photon rich.

So the BDT is not just learning an abstract classifier boundary. The expanded sample teaches it a specific physical background process: hard neutral-meson fragmentation that can create photon-like reconstructed clusters.
""",
        "04_candidate_level_particle_microscope_script.md": """# THE-44 Slide 4 Script - Candidate-Level Particle Microscope

This last slide is the event-level microscope. On the left is a high-BDT inclusive fake. It has a very signal-like score, but the local generator particle list shows a non-photon final-state cone dominated by neutral meson and charged-hadron fragments.

On the right is a truth-photon control with nearby activity. It still has a photon core, but the particle list explains why the BDT can be cautious: the local environment is not perfectly clean. This is the practical value of THE-44. We can now ask what physical object a BDT score corresponds to, candidate by candidate.
""",
    }
    for name, body in texts.items():
        path = outdir / name
        path.write_text(body)
        scripts[name] = str(path)
    return scripts


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    data = load_inputs()
    inc = high_inclusive(data["candidate"])
    paths = [
        slide_score_metric(data, args.outdir),
        slide_energy_coverage(data, inc, args.outdir),
        slide_particle_pipeline(data, inc, args.outdir),
        slide_candidate_microscope(data, inc, args.outdir),
    ]
    contact = contact_sheet(paths, args.outdir)

    auc_row = (
        data["fixed_auc"][(data["fixed_auc"]["scope"] == "direct_holdout_3x3") & (data["fixed_auc"]["validation_sample"] == "Jet12+20")]
        .set_index("model_training_sample")
        .loc[["Jet12+20", "Jet12+20+30+40"]]
    )
    gap_row = (
        data["median_gap"][data["median_gap"]["validation_sample"] == "Jet12+20"]
        .set_index("training_sample")
        .loc[["Jet12+20", "Jet12+20+30+40"]]
    )
    metrics = {
        "auc_start": float(auc_row["auc_0_20"].iloc[0]),
        "auc_end": float(auc_row["auc_0_20"].iloc[1]),
        "gap_start": float(gap_row["median_bdt_score_gap_binned"].iloc[0]),
        "gap_end": float(gap_row["median_bdt_score_gap_binned"].iloc[1]),
        "high_inc_count": float(len(inc)),
        "jet3040_count": float(inc["source"].isin(["Jet30", "Jet40"]).sum()),
        "meson_count": float(inc["truth_bucket"].eq("pi0/eta").sum()),
        "nonphoton_count": float(inc["cone_bucket"].eq("non-photon-rich").sum()),
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
            "source_fractions": str(MECH_DIR / "source_fractions_holdout_cent020.csv"),
            "fixed_validation_auc": str(MECH_DIR / "fixed_validation_training_effect_auc.csv"),
            "median_gap_summary": str(
                FIXED_DIR
                / "the8_branchA_fixed_sample_holdout_3x3_direct_0to20_truthSignal_inclusiveJet_slide_v34_mediangap_ppbasev3e_clean.median_gap_summary.csv"
            ),
        },
        "caveat": "Internal qualitative diagnostic from THE-44 saved candidate/particle lists; not a production-weighted purity estimate.",
    }
    manifest_path = args.outdir / "the44_deep_physics_sequence_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
