#!/usr/bin/env python3
"""Build a clean three-slide THE-44 Jet30/40 learning sequence."""

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
MECH_DIR = BASE / "diagnostics/hard_inclusive_jet3040_mechanism_20260605"
THE44_DIR = BASE / "diagnostics/the44_high_bdt_pythia_autopsy_20260605/slideReady/full_autopsy_20260605_2108"
DEFAULT_OUTDIR = THE44_DIR / "particle_intensive/jet3040_learning_sequence_20260606"

INK = "#1F2329"
MUTED = "#647083"
LIGHT = "#F6F7F9"
CARD = "#FFFFFF"
GRID = "#E3E7ED"
BLUE = "#2B6CB0"
ORANGE = "#F97316"
GREEN = "#2F9E44"
PURPLE = "#8E5A84"
RED = "#C2410C"
YELLOW = "#FFF3C4"

SAMPLE_COLORS = {
    "Jet12": "#6B8FB3",
    "Jet20": ORANGE,
    "Jet30": GREEN,
    "Jet40": PURPLE,
}

plt.rcParams.update(
    {
        "font.family": ["Times New Roman", "Times", "DejaVu Serif"],
        "axes.edgecolor": "#2B3038",
        "axes.labelcolor": INK,
        "xtick.color": INK,
        "ytick.color": INK,
    }
)


def canvas():
    fig = plt.figure(figsize=(16, 9), dpi=160)
    fig.patch.set_facecolor(LIGHT)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    return fig, ax


def txt(ax, x, y, s, size=18, color=INK, weight="normal", ha="left", va="top", **kwargs):
    ax.text(
        x,
        y,
        s,
        transform=ax.transAxes,
        fontsize=size,
        color=color,
        fontweight=weight,
        ha=ha,
        va=va,
        zorder=10,
        **kwargs,
    )


def title(ax, line1, line2):
    txt(ax, 0.055, 0.93, line1, size=34, weight="bold")
    txt(ax, 0.055, 0.875, line2, size=19, color=MUTED)


def card(ax, x, y, w, h, face=CARD, edge="#D8DDE5", radius=0.02, lw=1.2):
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        transform=ax.transAxes,
        boxstyle=f"round,pad=0.012,rounding_size={radius}",
        facecolor=face,
        edgecolor=edge,
        linewidth=lw,
        zorder=1,
    )
    ax.add_patch(patch)
    return patch


def inset(fig, x, y, w, h):
    return fig.add_axes([x, y, w, h])


def source_bucket(name: str) -> str:
    if name in {"pi0", "eta"}:
        return "pi0/eta"
    if name == "gamma":
        return "gamma"
    return "other"


def load_evidence():
    fixed = pd.read_csv(MECH_DIR / "fixed_validation_training_effect_auc.csv")
    wp80 = pd.read_csv(MECH_DIR / "wp80_full_scorecache_fake_rate_summary.csv")
    source = pd.read_csv(MECH_DIR / "source_fractions_holdout_cent020.csv")
    et_auc = pd.read_csv(MECH_DIR / "row_matched_holdout_cent020_et_metrics.csv")
    cand = pd.read_csv(THE44_DIR / "the44_candidate_autopsy_rows.csv")

    fixed = fixed[
        (fixed["scope"] == "direct_holdout_3x3")
        & (fixed["validation_sample"] == "Jet12+20")
    ].set_index("model_training_sample").loc[["Jet12+20", "Jet12+20+30", "Jet12+20+30+40"]]

    wp80 = wp80.set_index("sample").loc[["Jet12+20", "Jet12+20+30", "Jet12+20+30+40"]]

    source = source[
        (source["branch"] == "Jet12+20+30+40")
        & (source["scope"] == "holdout10")
        & (source["centrality_bin"] == "0-20%")
    ].copy()
    et_auc = et_auc[et_auc["sample"] == "Jet12+20+30+40"].copy()

    inc_high = cand[
        (cand["label"] == "high-BDT background")
        & cand["sample"].str.contains("Jet", na=False)
    ].copy()
    inc_high["sample_simple"] = inc_high["sample"].str.replace(" inclusive", "", regex=False)
    inc_high["truth_bucket"] = inc_high["truth_pid_name"].map(source_bucket)
    inc_high["particle_verdict"] = np.where(
        inc_high["final_state_photon_frac"] < 0.20,
        "non-photon cone",
        "photon-like/mixed",
    )
    return fixed, wp80, source, et_auc, inc_high


def finish(fig, out):
    fig.savefig(out, dpi=160)
    plt.close(fig)


def slide_auc_tail(fixed, wp80, outdir):
    fig, ax = canvas()
    title(
        ax,
        "AUC hides the part of the BDT gain we care about",
        "On the same Jet12+20 validation rows the AUC moves only slightly, while the high-score fake tail collapses.",
    )

    # Left: fixed-validation AUC.
    card(ax, 0.065, 0.205, 0.42, 0.60)
    txt(ax, 0.095, 0.755, "Global ranking metric", size=21, weight="bold")
    txt(ax, 0.095, 0.715, "same Jet12+20 validation rows", size=15, color=MUTED, weight="bold")
    auc_ax = inset(fig, 0.105, 0.34, 0.32, 0.31)
    labels = ["Jet12+20", "+Jet30", "+Jet40"]
    x = np.arange(3)
    auc = fixed["auc_all"].to_numpy()
    auc_ax.plot(x, auc, color=BLUE, marker="o", markersize=9, linewidth=4)
    auc_ax.set_xticks(x)
    auc_ax.set_xticklabels(labels, fontsize=13)
    auc_ax.set_ylim(0.848, 0.866)
    auc_ax.set_ylabel("AUC", fontsize=15)
    auc_ax.grid(True, color=GRID)
    auc_ax.spines[["top", "right"]].set_visible(False)
    auc_ax.tick_params(labelsize=12)
    auc_ax.text(2, auc[-1] + 0.001, f"+{auc[-1]-auc[0]:.3f}", ha="center", fontsize=17, color=BLUE, weight="bold")
    txt(ax, 0.095, 0.28, "AUC says: incremental", size=24, weight="bold", color=BLUE)

    # Right: WP80 fake rate.
    card(ax, 0.535, 0.205, 0.40, 0.60)
    txt(ax, 0.565, 0.755, "Tail-cleaning metric", size=21, weight="bold")
    txt(ax, 0.565, 0.715, "inclusive background passing WP80", size=15, color=MUTED, weight="bold")
    fake_ax = inset(fig, 0.595, 0.34, 0.27, 0.31)
    fake = wp80["wp80_background_fake_rate"].to_numpy() * 100.0
    fake_ax.bar(x, fake, color=[BLUE, ORANGE, GREEN], width=0.65)
    fake_ax.set_xticks(x)
    fake_ax.set_xticklabels(labels, fontsize=13)
    fake_ax.set_ylim(0, 48)
    fake_ax.set_ylabel("fake rate at WP80 (%)", fontsize=15)
    fake_ax.grid(True, axis="y", color=GRID)
    fake_ax.spines[["top", "right"]].set_visible(False)
    fake_ax.tick_params(labelsize=12)
    for i, v in enumerate(fake):
        fake_ax.text(i, v + 1.6, f"{v:.0f}%", ha="center", fontsize=18, weight="bold")
    txt(ax, 0.565, 0.28, "Score split says: large", size=24, weight="bold", color=GREEN)

    card(ax, 0.18, 0.08, 0.64, 0.105, face=YELLOW, edge="#E3C85C")
    txt(
        ax,
        0.50,
        0.142,
        "The meeting-visible improvement is a high-BDT tail cleanup, not a headline-AUC story.",
        size=22,
        weight="bold",
        ha="center",
        va="center",
    )
    out = outdir / "01_auc_misses_tail_cleanup.png"
    finish(fig, out)
    return out


def slide_source_energy(source, et_auc, outdir):
    fig, ax = canvas()
    title(
        ax,
        "Jet40 supplies the hard background the old model barely saw",
        "In 0-20% centrality, the expanded sample's high-E$_T$ inclusive background is dominated by Jet40.",
    )

    plot = inset(fig, 0.10, 0.24, 0.62, 0.50)
    et_bins = ["15-20", "20-25", "25-30", "30-35"]
    bottoms = np.zeros(len(et_bins))
    for sample in ["Jet12", "Jet20", "Jet30", "Jet40"]:
        vals = []
        for et in et_bins:
            row = source[(source["et_bin"] == et) & (source["source_bin"] == sample)]
            vals.append(float(row["source_fraction_weighted"].iloc[0]) if len(row) else 0.0)
        vals = np.asarray(vals) * 100.0
        plot.bar(np.arange(len(et_bins)), vals, bottom=bottoms, color=SAMPLE_COLORS[sample], width=0.72, label=sample)
        for i, (b, v) in enumerate(zip(bottoms, vals)):
            if v >= 15:
                plot.text(i, b + v / 2, f"{v:.0f}%", ha="center", va="center", fontsize=14, color="white", weight="bold")
        bottoms += vals
    plot.set_ylim(0, 100)
    plot.set_ylabel("weighted inclusive-background share (%)", fontsize=16)
    plot.set_xticks(np.arange(len(et_bins)))
    plot.set_xticklabels(et_bins, fontsize=15)
    plot.set_xlabel("cluster E$_T$ bin [GeV]", fontsize=16)
    plot.grid(True, axis="y", color=GRID)
    plot.spines[["top", "right"]].set_visible(False)
    plot.tick_params(axis="y", labelsize=14)

    gains = et_auc.set_index("et_label").loc[et_bins]["auc_delta_vs_jet12_20"].to_numpy()
    for i, gain in enumerate(gains):
        plot.text(i, 105, f"AUC gain +{gain:.2f}", ha="center", fontsize=14, color=GREEN, weight="bold", clip_on=False)

    # Right legend / explanation.
    card(ax, 0.765, 0.255, 0.18, 0.48)
    txt(ax, 0.79, 0.70, "Parent jet-energy bins", size=18, weight="bold")
    rows = [
        ("Jet20", "pThat 21-32 GeV", SAMPLE_COLORS["Jet20"]),
        ("Jet30", "pThat 32-42 GeV", SAMPLE_COLORS["Jet30"]),
        ("Jet40", "pThat 42-100 GeV", SAMPLE_COLORS["Jet40"]),
    ]
    for i, (lab, desc, col) in enumerate(rows):
        y = 0.63 - 0.105 * i
        ax.add_patch(Rectangle((0.795, y - 0.025), 0.026, 0.026, transform=ax.transAxes, facecolor=col, edgecolor="none", zorder=5))
        txt(ax, 0.835, y, lab, size=18, weight="bold", va="center")
        txt(ax, 0.835, y - 0.035, desc, size=13.5, color=MUTED, va="center")
    txt(ax, 0.79, 0.335, "By 25-35 GeV,\nJet40 is 76-81%\nof this background.", size=21, weight="bold", color=PURPLE, linespacing=1.05)

    card(ax, 0.14, 0.085, 0.72, 0.095, face=YELLOW, edge="#E3C85C")
    txt(
        ax,
        0.50,
        0.14,
        "The added samples are not generic statistics:\nthey populate the hard background region where the score boundary improves.",
        size=18.5,
        weight="bold",
        ha="center",
        va="center",
        linespacing=1.05,
    )
    out = outdir / "02_jet40_populates_hard_background.png"
    finish(fig, out)
    return out


def slide_particle_history(inc_high, outdir):
    fig, ax = canvas()
    title(
        ax,
        "The new hard background is neutral-meson fragmentation",
        "THE-44 particle lists identify what the BDT must learn to push down after Jet30/40 are included.",
    )

    total = len(inc_high)
    jet3040 = int(inc_high["sample_simple"].isin(["Jet30", "Jet40"]).sum())
    meson = int(inc_high["truth_bucket"].eq("pi0/eta").sum())
    non_photon = int(inc_high["particle_verdict"].eq("non-photon cone").sum())

    xs = [0.08, 0.32, 0.56, 0.80]
    y = 0.44
    w, h = 0.17, 0.20
    steps = [
        ("Parent bins", f"{jet3040}/{total}", "Jet30/40\npThat 32-100 GeV", GREEN),
        ("Truth source", f"{meson}/{total}", "mostly pi0/eta", ORANGE),
        ("Visible cone", f"{non_photon}/{total}", "non-photon rich", RED),
        ("BDT lesson", "tail lower", "not just AUC up", BLUE),
    ]
    for i, (label, big, sub, color) in enumerate(steps):
        card(ax, xs[i], y, w, h, face=color, edge=color, radius=0.025)
        txt(ax, xs[i] + w / 2, y + h - 0.045, label, size=18, color="white", weight="bold", ha="center")
        txt(ax, xs[i] + w / 2, y + 0.105, big, size=34, color="white", weight="bold", ha="center", va="center")
        txt(ax, xs[i] + w / 2, y + 0.045, sub, size=17, color="white", ha="center", va="center", linespacing=0.95)
        if i < 3:
            ax.add_patch(
                FancyArrowPatch(
                    (xs[i] + w + 0.015, y + h / 2),
                    (xs[i + 1] - 0.015, y + h / 2),
                    transform=ax.transAxes,
                    arrowstyle="-|>",
                    mutation_scale=24,
                    linewidth=3.0,
                    color="#6E7784",
                    zorder=5,
                )
            )

    # Simple readout cards.
    card(ax, 0.095, 0.24, 0.265, 0.11)
    txt(ax, 0.125, 0.315, "91%", size=34, weight="bold", color=GREEN, va="center")
    txt(ax, 0.215, 0.315, "of high-BDT inclusive fakes\ncome from Jet30/40", size=16.5, va="center", linespacing=1.05)
    card(ax, 0.385, 0.24, 0.265, 0.11)
    txt(ax, 0.415, 0.315, "79%", size=34, weight="bold", color=ORANGE, va="center")
    txt(ax, 0.505, 0.315, "have pi0/eta truth\nas the candidate source", size=16.5, va="center", linespacing=1.05)
    card(ax, 0.675, 0.24, 0.245, 0.11)
    txt(ax, 0.705, 0.315, "86%", size=34, weight="bold", color=RED, va="center")
    txt(ax, 0.795, 0.315, "have non-photon-rich\nfinal-state cones", size=16.5, va="center", linespacing=1.05)

    card(ax, 0.14, 0.085, 0.72, 0.095, face=YELLOW, edge="#E3C85C")
    txt(
        ax,
        0.50,
        0.14,
        "This is the physical class Jet30/40 teach the BDT:\nhard neutral-meson backgrounds that should move out of the signal-like tail.",
        size=18.5,
        weight="bold",
        ha="center",
        va="center",
        linespacing=1.05,
    )
    out = outdir / "03_particle_history_lesson.png"
    finish(fig, out)
    return out


def contact_sheet(paths, outdir):
    fig = plt.figure(figsize=(16, 9), dpi=160)
    fig.patch.set_facecolor(LIGHT)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    txt(ax, 0.05, 0.94, "Jet30/40 BDT Learning Sequence", size=30, weight="bold")
    labels = ["1. AUC misses tail cleanup", "2. Jet40 supplies hard background", "3. Particle-history lesson"]
    for i, path in enumerate(paths):
        img = mpimg.imread(path)
        x = 0.05 + i * 0.315
        iax = fig.add_axes([x, 0.23, 0.285, 0.50])
        iax.imshow(img)
        iax.set_axis_off()
        iax.add_patch(Rectangle((0, 0), 1, 1, transform=iax.transAxes, facecolor="none", edgecolor="#C9D0DA", linewidth=2))
        txt(ax, x, 0.78, labels[i], size=17, weight="bold")
    txt(ax, 0.05, 0.12, "Recommended delivery: use the three slides separately; the contact sheet is only for review.", size=17, color=MUTED)
    out = outdir / "jet3040_learning_sequence_contact_sheet.png"
    finish(fig, out)
    return out


def write_scripts(outdir, metrics):
    scripts = {}
    script1 = outdir / "01_auc_misses_tail_cleanup_script.md"
    script1.write_text(
        f"""# THE-44 Slide 1 Script - AUC Hides Tail Cleanup

First, I want to separate two different ideas. AUC is a global ranking metric, so it averages over a very large number of easy signal-background pairs. On the same Jet12 plus Jet20 validation rows, adding Jet30 and Jet40 only moves the AUC from {metrics['auc0']:.3f} to {metrics['auc2']:.3f}.

But the thing we actually see by eye in the score distributions is the high-score background tail. At the 80 percent signal-efficiency working point, the inclusive fake rate falls from {metrics['fake0']:.0f} percent to {metrics['fake2']:.0f} percent. So the score plot can look dramatically cleaner even when the headline AUC gain looks incremental.
"""
    )
    scripts["01_auc"] = str(script1)

    script2 = outdir / "02_jet40_populates_hard_background_script.md"
    script2.write_text(
        """# THE-44 Slide 2 Script - Jet40 Populates The Hard Background

Now the reason Jet30 and Jet40 matter is that they populate a different parent-jet energy region. This plot shows the source composition of inclusive background candidates in 0 to 20 percent centrality, broken out by cluster energy.

The key thing to notice is the purple Jet40 component. Above about 25 GeV cluster energy, Jet40 supplies most of the inclusive background. That is exactly the hard-background region where the old Jet12 plus Jet20 training had little support. So Jet30 and Jet40 are not just increasing the row count; they are adding the hard examples needed to define the boundary.
"""
    )
    scripts["02_source"] = str(script2)

    script3 = outdir / "03_particle_history_lesson_script.md"
    script3.write_text(
        f"""# THE-44 Slide 3 Script - Particle-History Lesson

The particle-list autopsy gives the physical identity of that hard background. In the final high-BDT inclusive tail, Jet30 and Jet40 supply {metrics['jet3040_frac']:.0f} percent of the candidates. The truth source is mostly pi0 or eta, and the visible final-state cone is usually non-photon rich.

That is the important interpretation. The BDT improvement is not only a numerical training improvement. The expanded sample teaches the model a specific physical class: hard neutral-meson and fragmentation backgrounds near the photon-candidate energy range. Once those examples are represented in training, the BDT can more consistently move that class out of the signal-like score tail.
"""
    )
    scripts["03_particle"] = str(script3)
    return scripts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    fixed, wp80, source, et_auc, inc_high = load_evidence()
    paths = [
        slide_auc_tail(fixed, wp80, args.outdir),
        slide_source_energy(source, et_auc, args.outdir),
        slide_particle_history(inc_high, args.outdir),
    ]
    contact = contact_sheet(paths, args.outdir)

    fake = wp80["wp80_background_fake_rate"].to_numpy() * 100.0
    total = len(inc_high)
    metrics = {
        "auc0": float(fixed["auc_all"].iloc[0]),
        "auc2": float(fixed["auc_all"].iloc[2]),
        "fake0": float(fake[0]),
        "fake2": float(fake[2]),
        "jet3040_frac": float(100.0 * inc_high["sample_simple"].isin(["Jet30", "Jet40"]).sum() / total),
        "meson_frac": float(100.0 * inc_high["truth_bucket"].eq("pi0/eta").sum() / total),
        "non_photon_frac": float(100.0 * inc_high["particle_verdict"].eq("non-photon cone").sum() / total),
    }
    scripts = write_scripts(args.outdir, metrics)

    manifest = {
        "status": "READY",
        "slides": [str(p) for p in paths],
        "contact_sheet": str(contact),
        "scripts": scripts,
        "metrics": metrics,
        "inputs": {
            "fixed_validation_training_effect_auc": str(MECH_DIR / "fixed_validation_training_effect_auc.csv"),
            "wp80_full_scorecache_fake_rate_summary": str(MECH_DIR / "wp80_full_scorecache_fake_rate_summary.csv"),
            "source_fractions_holdout_cent020": str(MECH_DIR / "source_fractions_holdout_cent020.csv"),
            "row_matched_holdout_cent020_et_metrics": str(MECH_DIR / "row_matched_holdout_cent020_et_metrics.csv"),
            "the44_candidate_autopsy_rows": str(THE44_DIR / "the44_candidate_autopsy_rows.csv"),
        },
        "caveat": "Particle identities are final-model THE-44 autopsy rows; paired old/new particle-history migration is a follow-up diagnostic.",
    }
    manifest_path = args.outdir / "jet3040_learning_sequence_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
