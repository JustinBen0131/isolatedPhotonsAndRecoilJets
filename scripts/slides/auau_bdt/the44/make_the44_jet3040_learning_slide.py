#!/usr/bin/env python3
"""Build the THE-44 Jet30/40 BDT-learning explanation slide."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch


REPO = Path(__file__).resolve().parents[4]
BASE = REPO / "dataOutput/auauTightBDTValidation/THE8_branchA_ladder_scorecache_fullstat_20260527"
MECH_DIR = BASE / "diagnostics/hard_inclusive_jet3040_mechanism_20260605"
THE44_DIR = (
    BASE
    / "diagnostics/the44_high_bdt_pythia_autopsy_20260605/slideReady/full_autopsy_20260605_2108"
)
DEFAULT_OUTDIR = THE44_DIR / "particle_intensive/jet3040_learning_story_20260606"

INK = "#20242A"
MUTED = "#5D6776"
GRID = "#E6E8EC"
BLUE = "#2B6CB0"
ORANGE = "#F97316"
GREEN = "#2F9E44"
PURPLE = "#7E3FB2"
RED = "#C2410C"
CARD = "#FFFFFF"
BG = "#F6F7F9"

SAMPLE_COLORS = {
    "Jet12": "#6B8FB3",
    "Jet20": "#F97316",
    "Jet30": "#2F9E44",
    "Jet40": "#8E5A84",
}


def add_card(ax, xy, wh, title=None, face=CARD, edge="#D7DCE2", radius=0.018, title_size=22):
    x, y = xy
    w, h = wh
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle=f"round,pad=0.012,rounding_size={radius}",
        transform=ax.transAxes,
        facecolor=face,
        edgecolor=edge,
        linewidth=1.1,
        zorder=1,
    )
    ax.add_patch(patch)
    if title:
        ax.text(
            x + 0.025,
            y + h - 0.045,
            title,
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=title_size,
            fontweight="bold",
            color=INK,
            zorder=2,
        )
    return patch


def text(ax, x, y, s, size=18, color=INK, weight="normal", ha="left", va="top", **kwargs):
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
        zorder=5,
        **kwargs,
    )


def inset(fig, bounds):
    return fig.add_axes(bounds)


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

    fixed_j1220 = fixed[
        (fixed["scope"] == "direct_holdout_3x3")
        & (fixed["validation_sample"] == "Jet12+20")
    ][["model_training_sample", "auc_all"]].copy()
    fixed_j1220 = fixed_j1220.set_index("model_training_sample").loc[
        ["Jet12+20", "Jet12+20+30", "Jet12+20+30+40"]
    ]

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

    return fixed_j1220, wp80, source, et_auc, inc_high


def draw_auc_vs_tail(fig, ax, fixed_j1220, wp80):
    add_card(ax, (0.045, 0.565), (0.385, 0.335), "1  AUC misses the tail cleanup", title_size=20)
    text(
        ax,
        0.07,
        0.822,
        "Same Jet12+20 validation rows:",
        size=15.5,
        color=MUTED,
        weight="bold",
    )

    auc_ax = inset(fig, [0.075, 0.655, 0.16, 0.125])
    labels = ["J12+20", "+J30", "+J40"]
    auc_vals = fixed_j1220["auc_all"].to_numpy()
    x = np.arange(len(labels))
    auc_ax.plot(x, auc_vals, color=BLUE, marker="o", linewidth=3.0)
    auc_ax.set_xticks(x)
    auc_ax.set_xticklabels(labels, fontsize=9)
    auc_ax.set_ylim(0.848, 0.866)
    auc_ax.set_ylabel("AUC", fontsize=10)
    auc_ax.grid(True, color=GRID, linewidth=0.8)
    auc_ax.spines[["top", "right"]].set_visible(False)
    auc_ax.tick_params(axis="y", labelsize=9)
    auc_ax.text(
        2,
        auc_vals[-1] + 0.0009,
        f"+{auc_vals[-1] - auc_vals[0]:.3f}",
        ha="center",
        fontsize=10.5,
        color=BLUE,
        fontweight="bold",
    )

    fake_ax = inset(fig, [0.265, 0.655, 0.14, 0.125])
    fake = wp80.set_index("sample").loc[["Jet12+20", "Jet12+20+30", "Jet12+20+30+40"]]
    fake_rates = fake["wp80_background_fake_rate"].to_numpy() * 100.0
    fake_ax.bar(x, fake_rates, color=[BLUE, ORANGE, GREEN], width=0.62)
    fake_ax.set_xticks(x)
    fake_ax.set_xticklabels(labels, fontsize=9)
    fake_ax.set_ylim(0, 48)
    fake_ax.set_ylabel("WP80 fake %", fontsize=10)
    fake_ax.grid(True, axis="y", color=GRID, linewidth=0.8)
    fake_ax.spines[["top", "right"]].set_visible(False)
    fake_ax.tick_params(axis="y", labelsize=9)
    for i, v in enumerate(fake_rates):
        fake_ax.text(i, v + 1.4, f"{v:.0f}%", ha="center", fontsize=10.5, fontweight="bold")

    text(
        ax,
        0.07,
        0.620,
        "WP80 high-score fake rate:\n42% -> 26% -> 17%",
        size=18.0,
        color=INK,
        weight="bold",
        linespacing=1.08,
    )
    text(
        ax,
        0.07,
        0.574,
        "Visual split improves more than AUC suggests.",
        size=12.8,
        color=MUTED,
    )


def draw_source_energy(fig, ax, source, et_auc):
    add_card(ax, (0.455, 0.565), (0.50, 0.335), "2  Jet30/40 add hard parent-jet bins", title_size=20)
    text(
        ax,
        0.48,
        0.822,
        "0-20% holdout; source mix by cluster E$_T$  (purple=Jet40, green=Jet30)",
        size=14.5,
        color=MUTED,
        weight="bold",
    )
    bar_ax = inset(fig, [0.505, 0.638, 0.39, 0.15])
    et_bins = ["15-20", "20-25", "25-30", "30-35"]
    bottoms = np.zeros(len(et_bins))
    for sample in ["Jet12", "Jet20", "Jet30", "Jet40"]:
        vals = []
        for et in et_bins:
            row = source[(source["et_bin"] == et) & (source["source_bin"] == sample)]
            vals.append(float(row["source_fraction_weighted"].iloc[0]) if len(row) else 0.0)
        bar_ax.bar(
            np.arange(len(et_bins)),
            np.asarray(vals) * 100.0,
            bottom=bottoms,
            color=SAMPLE_COLORS[sample],
            label=sample,
            width=0.72,
        )
        bottoms += np.asarray(vals) * 100.0

    bar_ax.set_ylim(0, 100)
    bar_ax.set_ylabel("background share", fontsize=10)
    bar_ax.set_xticks(np.arange(len(et_bins)))
    bar_ax.set_xticklabels(et_bins, fontsize=10)
    bar_ax.set_xlabel("")
    bar_ax.grid(True, axis="y", color=GRID, linewidth=0.8)
    bar_ax.spines[["top", "right"]].set_visible(False)
    bar_ax.tick_params(axis="y", labelsize=9)
    gains = et_auc.set_index("et_label").loc[et_bins]["auc_delta_vs_jet12_20"].to_numpy()
    for i, gain in enumerate(gains):
        bar_ax.text(
            i,
            95,
            f"AUC +{gain:.2f}",
            ha="center",
            va="center",
            fontsize=10.2,
            color="white",
            fontweight="bold",
            clip_on=True,
        )

    text(
        ax,
        0.48,
        0.607,
        "Jet40 is most of the 25-35 GeV background.\nThat hard phase space is missing from Jet12/20.",
        size=13.5,
        color=INK,
        linespacing=1.05,
    )


def draw_particle_chain(ax, inc_high):
    add_card(ax, (0.045, 0.155), (0.91, 0.39), None)
    text(ax, 0.07, 0.505, "3  Particle identity of the hard background tail", size=20.5, weight="bold", color=INK)
    total = len(inc_high)
    jet3040 = int(inc_high["sample_simple"].isin(["Jet30", "Jet40"]).sum())
    meson = int(inc_high["truth_bucket"].eq("pi0/eta").sum())
    gamma = int(inc_high["truth_bucket"].eq("gamma").sum())
    other = total - meson - gamma
    non_photon = int(inc_high["particle_verdict"].eq("non-photon cone").sum())

    x0, y = 0.075, 0.355
    box_w, box_h = 0.17, 0.13
    steps = [
        (x0, "Jet30/40 parent bins", f"{jet3040}/{total}", "pThat 32-100 GeV", GREEN),
        (x0 + 0.24, "Neutral-meson truth", f"{meson}/{total}", "mostly pi0/eta", ORANGE),
        (x0 + 0.48, "Visible cone", f"{non_photon}/{total}", "non-photon rich", RED),
        (x0 + 0.72, "BDT lesson", "tail lower", "not just AUC up", BLUE),
    ]
    for i, (x, title, big, sub, color) in enumerate(steps):
        patch = FancyBboxPatch(
            (x, y),
            box_w,
            box_h,
            boxstyle="round,pad=0.012,rounding_size=0.018",
            transform=ax.transAxes,
            facecolor=color,
            edgecolor="none",
            zorder=3,
        )
        ax.add_patch(patch)
        text(ax, x + box_w / 2, y + box_h - 0.028, title, size=12.8, color="white", weight="bold", ha="center")
        text(ax, x + box_w / 2, y + 0.070, big, size=28, color="white", weight="bold", ha="center", va="center")
        text(ax, x + box_w / 2, y + 0.027, sub, size=12.8, color="white", ha="center", va="center")
        if i < len(steps) - 1:
            ax.annotate(
                "",
                xy=(x + box_w + 0.045, y + box_h / 2),
                xytext=(x + box_w + 0.01, y + box_h / 2),
                xycoords=ax.transAxes,
                arrowprops=dict(arrowstyle="-|>", linewidth=2.8, color="#7A828E"),
                zorder=4,
            )

    text(
        ax,
        0.09,
        0.305,
        "Particle-list autopsy of the final high-BDT inclusive tail:",
        size=18,
        color=INK,
        weight="bold",
    )
    text(
        ax,
        0.09,
        0.267,
        (
            f"Jet30/40 supply {jet3040/total:.0%}; pi0/eta truth supplies {meson/total:.0%}; "
            f"non-photon final-state cones supply {non_photon/total:.0%}."
        ),
        size=17.5,
        color=INK,
        linespacing=1.13,
    )
    text(
        ax,
        0.09,
        0.220,
        (
            "Interpretation: Jet30/40 are not just more rows. They add hard neutral-meson / fragmentation backgrounds\n"
            "near photon-candidate E$_T$, so the BDT can define a cleaner low-score boundary."
        ),
        size=14.8,
        color=MUTED,
        linespacing=1.14,
    )


def write_script(outdir: Path, metrics: dict) -> Path:
    path = outdir / "the44_jet3040_learning_slide_script.md"
    body = f"""# THE-44 Slide Script - What Jet30/40 Teach The BDT

The point of this slide is that Jet30 and Jet40 are not just giving us more statistics. They are adding a new background phase space that the earlier training did not represent well.

On the upper left, the fixed Jet12 plus Jet20 validation AUC only moves from {metrics['fixed_auc0']:.3f} to {metrics['fixed_auc2']:.3f}. That is a small change, so AUC alone makes the gain look modest.

But the working-point behavior tells a different story. At the same 80 percent signal-efficiency target, the inclusive fake rate falls from {metrics['fake0']:.0f} percent to {metrics['fake2']:.0f} percent as Jet30 and Jet40 are added. That is the qualitative separation we see in the score plots.

The upper-right panel explains why. In the central 0 to 20 percent bin, the harder cluster-energy background is dominated by Jet30 and especially Jet40. Those samples are the parent jet-energy bins that populate the difficult hard-background tail.

The bottom panel connects that source story to the particle-list autopsy. In the high-BDT inclusive tail, Jet30 and Jet40 supply {metrics['jet3040_frac']:.0f} percent of the candidates, pi0 and eta supply {metrics['meson_frac']:.0f} percent of the truth sources, and {metrics['non_photon_frac']:.0f} percent of the visible cones are non-photon rich.

So the simple takeaway is: Jet30 and Jet40 improve the BDT because they teach it the hard neutral-meson and fragmentation backgrounds that sit near the photon-candidate energy range. AUC averages over the whole ranking problem, but the visual improvement is concentrated in cleaning this high-score background tail.
"""
    path.write_text(body)
    return path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    fixed_j1220, wp80, source, et_auc, inc_high = load_evidence()

    fig = plt.figure(figsize=(16, 9), dpi=160)
    fig.patch.set_facecolor(BG)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)

    text(ax, 0.045, 0.945, "What Jet30/40 Teach The BDT", size=29, weight="bold", color=INK)
    text(
        ax,
        0.045,
        0.905,
        "Hard Jet30/40 backgrounds enter training: neutral mesons and fragmentation near photon-candidate E$_T$.",
        size=16.3,
        color=MUTED,
    )
    text(
        ax,
        0.045,
        0.878,
        "The high-score fake tail cleans up; global AUC averages over that localized gain.",
        size=16.3,
        color=MUTED,
    )
    text(
        ax,
        0.865,
        0.943,
        "THE-44 + THE-8 evidence",
        size=13.5,
        color=MUTED,
        weight="bold",
        ha="right",
    )

    draw_auc_vs_tail(fig, ax, fixed_j1220, wp80)
    draw_source_energy(fig, ax, source, et_auc)
    draw_particle_chain(ax, inc_high)

    caveat = (
        "Evidence: THE-8 training-effect metrics plus THE-44 Pythia autopsy. "
        "Caveat: current particle rows do not yet store paired old/new model score migration for the same candidate."
    )
    text(ax, 0.047, 0.085, caveat, size=12.2, color="#6B7280")

    out_png = args.outdir / "the44_jet3040_bdt_learning_slide.png"
    fig.savefig(out_png, dpi=160)
    plt.close(fig)

    fake = wp80.set_index("sample").loc[["Jet12+20", "Jet12+20+30", "Jet12+20+30+40"]]
    total = len(inc_high)
    metrics = {
        "fixed_auc0": float(fixed_j1220.loc["Jet12+20", "auc_all"]),
        "fixed_auc2": float(fixed_j1220.loc["Jet12+20+30+40", "auc_all"]),
        "fake0": float(fake.loc["Jet12+20", "wp80_background_fake_rate"] * 100.0),
        "fake2": float(fake.loc["Jet12+20+30+40", "wp80_background_fake_rate"] * 100.0),
        "jet3040_frac": float(100.0 * inc_high["sample_simple"].isin(["Jet30", "Jet40"]).sum() / total),
        "meson_frac": float(100.0 * inc_high["truth_bucket"].eq("pi0/eta").sum() / total),
        "non_photon_frac": float(100.0 * inc_high["particle_verdict"].eq("non-photon cone").sum() / total),
    }
    script = write_script(args.outdir, metrics)

    manifest = {
        "status": "READY",
        "png": str(out_png),
        "script": str(script),
        "inputs": {
            "fixed_validation_training_effect_auc": str(MECH_DIR / "fixed_validation_training_effect_auc.csv"),
            "wp80_full_scorecache_fake_rate_summary": str(MECH_DIR / "wp80_full_scorecache_fake_rate_summary.csv"),
            "source_fractions_holdout_cent020": str(MECH_DIR / "source_fractions_holdout_cent020.csv"),
            "row_matched_holdout_cent020_et_metrics": str(MECH_DIR / "row_matched_holdout_cent020_et_metrics.csv"),
            "the44_candidate_autopsy_rows": str(THE44_DIR / "the44_candidate_autopsy_rows.csv"),
        },
        "metrics": metrics,
        "caveat": "Particle identities are final-model THE-44 autopsy rows; exact paired old/new particle-history score migration needs a follow-up paired diagnostic.",
    }
    manifest_path = args.outdir / "the44_jet3040_bdt_learning_slide_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
