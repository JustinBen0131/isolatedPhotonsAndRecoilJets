#!/usr/bin/env python3
"""Build particle-intensive THE-44 BDT autopsy story plots.

This script reads the full THE-44 candidate/particle CSV package and produces
slide-facing PNGs that explain what particles dominate high-BDT backgrounds,
how the BDT score relates to the particle neighborhood, and why the higher
energy inclusive samples are the cleanest diagnostic source.
"""

from __future__ import annotations

import argparse
import json
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
from matplotlib.patches import FancyBboxPatch, Patch, PathPatch
from matplotlib.path import Path as MplPath


DEFAULT_BASE = Path(
    "dataOutput/auauTightBDTValidation/"
    "THE8_branchA_ladder_scorecache_fullstat_20260527/diagnostics/"
    "the44_high_bdt_pythia_autopsy_20260605/slideReady/full_autopsy_20260605_2108"
)

GROUP_COLORS = {
    "photon": "#2F6FED",
    "neutral_meson": "#D95F02",
    "charged_hadron": "#1B9E77",
    "lepton": "#6A3D9A",
    "other": "#777777",
    "quark_gluon": "#7B6D8D",
}

GROUP_LABELS = {
    "photon": "Photon",
    "neutral_meson": "Neutral meson",
    "charged_hadron": "Charged hadron",
    "lepton": "Lepton",
    "other": "Other",
    "quark_gluon": "Quark/gluon",
}

SAMPLE_ORDER = ["Jet20 inclusive", "Jet30 inclusive", "Jet40 inclusive", "Photon12 signal", "Photon20 signal"]
INCLUSIVE_ORDER = ["Jet20 inclusive", "Jet30 inclusive", "Jet40 inclusive"]


def wrap(text: str, width: int = 42) -> str:
    return "\n".join(textwrap.wrap(text, width=width))


def draw_note(ax, x, y, w, h, title, body, face="#F7F7F7", edge="#D0D0D0", body_width=42):
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.012,rounding_size=0.014",
        transform=ax.transAxes,
        facecolor=face,
        edgecolor=edge,
        linewidth=1.1,
    )
    ax.add_patch(patch)
    ax.text(x + 0.025, y + h - 0.05, title, transform=ax.transAxes, va="top", fontsize=12.5, fontweight="bold")
    ax.text(
        x + 0.025,
        y + h - 0.12,
        wrap(body, body_width),
        transform=ax.transAxes,
        va="top",
        fontsize=8.8,
        linespacing=1.08,
    )


def source_bucket(name: str) -> str:
    if name in {"pi0", "eta"}:
        return "neutral meson"
    if name == "gamma":
        return "photon"
    if name in {"e+", "e-"}:
        return "electron"
    if name in {"pi+", "pi-", "K+", "K-", "p", "pbar"}:
        return "charged hadron"
    if name in {"D0", "D+", "D-", "521"}:
        return "heavy flavor"
    return "other"


def parent_bucket(row: pd.Series) -> str:
    parents = {int(row["parent0_pdg"]), int(row["parent1_pdg"])}
    abs_parents = {abs(v) for v in parents}
    if 22 in abs_parents:
        return "photon parent"
    if 11 in abs_parents:
        return "electron parent"
    if any(v in {1, 2, 3, 4, 5, 6, 21} for v in abs_parents):
        return "quark/gluon parent"
    if any(v in {111, 113, 211, 213, 221, 223, 313, 323, 331} for v in abs_parents):
        return "meson-resonance parent"
    return "other / unlisted parent"


def normalize_rows(frame: pd.DataFrame, value_cols: list[str]) -> pd.DataFrame:
    out = frame.copy()
    denom = out[value_cols].sum(axis=1).replace(0, np.nan)
    out[value_cols] = out[value_cols].div(denom, axis=0).fillna(0.0)
    return out


def load_tables(base: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    candidates = pd.read_csv(base / "the44_candidate_autopsy_rows.csv")
    particles = pd.read_csv(base / "the44_particle_rows.csv")
    particles = particles.merge(
        candidates[
            [
                "candidate_uid",
                "label",
                "truth_pid_name",
                "auau_tight_bdt_score",
                "cluster_Et",
                "reco_eiso",
                "final_state_non_photon_frac",
                "final_state_photon_frac",
            ]
        ],
        on="candidate_uid",
        how="left",
    )
    final_particles = particles[particles["status"] == 1].copy()
    final_particles["parent_bucket"] = final_particles.apply(parent_bucket, axis=1)
    return candidates, particles, final_particles


def make_particle_taxonomy_slide(candidates: pd.DataFrame, final_particles: pd.DataFrame, outdir: Path) -> Path:
    out = outdir / "the44_particle_taxonomy_slide.png"
    high = candidates[candidates["label"] == "high-BDT background"].copy()
    mid = candidates[candidates["label"] == "middling true photon"].copy()
    inc_high = high[high["sample"].str.contains("Jet", na=False)]

    fig = plt.figure(figsize=(16, 9), dpi=170)
    fig.patch.set_facecolor("white")
    gs = GridSpec(12, 16, figure=fig, left=0.09, right=0.97, top=0.825, bottom=0.08, wspace=1.25, hspace=1.65)

    fig.text(0.055, 0.95, "What the High-BDT Backgrounds Are Made Of", fontsize=25, fontweight="bold", family="serif")
    fig.text(
        0.055,
        0.912,
        "Candidate score is joined to final-state particles and HepMC parent context inside Delta R < 0.3",
        fontsize=14.5,
        color="#333333",
    )
    fig.legend(
        handles=[
            Patch(facecolor=GROUP_COLORS["photon"], label="Photon"),
            Patch(facecolor=GROUP_COLORS["neutral_meson"], label="Neutral meson"),
            Patch(facecolor=GROUP_COLORS["charged_hadron"], label="Charged hadron"),
            Patch(facecolor=GROUP_COLORS["lepton"], label="Lepton"),
            Patch(facecolor=GROUP_COLORS["other"], label="Other"),
        ],
        loc="upper left",
        bbox_to_anchor=(0.055, 0.885),
        ncol=5,
        frameon=False,
        fontsize=9.2,
        handlelength=1.2,
        columnspacing=1.2,
    )

    ax_comp = fig.add_subplot(gs[:5, :5])
    fp_high = final_particles[final_particles["label"] == "high-BDT background"].copy()
    comp = (
        fp_high.groupby(["sample", "group"])["pt"]
        .sum()
        .unstack(fill_value=0.0)
        .reindex(SAMPLE_ORDER)
        .fillna(0.0)
    )
    groups = ["photon", "neutral_meson", "charged_hadron", "lepton", "other"]
    comp = normalize_rows(comp.reset_index(), groups).set_index("sample")
    y = np.arange(len(comp))
    left = np.zeros(len(comp))
    for group in groups:
        vals = comp[group].to_numpy()
        ax_comp.barh(y, vals, left=left, color=GROUP_COLORS[group], height=0.58, label=GROUP_LABELS[group])
        left += vals
    short_index = [s.replace(" inclusive", "").replace(" signal", "") for s in comp.index]
    ax_comp.set_yticks(y)
    ax_comp.set_yticklabels(short_index, fontsize=9.7)
    ax_comp.invert_yaxis()
    ax_comp.set_xlim(0, 1)
    ax_comp.set_xlabel("fraction of final-state pT near candidate")
    ax_comp.set_title("High-BDT particle mix", fontsize=13, fontweight="bold")
    ax_comp.grid(True, axis="x", color="#E8E8E8")

    ax_truth = fig.add_subplot(gs[:5, 5:11])
    truth_counts = inc_high["truth_pid_name"].value_counts()
    truth_plot = truth_counts.head(7).sort_values()
    colors = ["#D95F02" if source_bucket(x) == "neutral meson" else "#2F6FED" if x == "gamma" else "#777777" for x in truth_plot.index]
    ax_truth.barh(truth_plot.index, truth_plot.values, color=colors)
    ax_truth.set_xlabel("candidate count")
    ax_truth.set_title("Inclusive high-BDT truth source", fontsize=13, fontweight="bold")
    ax_truth.grid(True, axis="x", color="#E8E8E8")
    ax_truth.text(
        0.98,
        0.05,
        f"pi0 + eta = {truth_counts.get('pi0', 0) + truth_counts.get('eta', 0)} / {len(inc_high)}",
        transform=ax_truth.transAxes,
        ha="right",
        fontsize=10.5,
        fontweight="bold",
        color="#7A3B00",
    )

    ax_counts = fig.add_subplot(gs[:5, 12:])
    counts = (
        candidates[candidates["sample"].isin(INCLUSIVE_ORDER)]
        .groupby(["sample", "label"])["candidate_uid"]
        .count()
        .unstack(fill_value=0)
        .reindex(INCLUSIVE_ORDER)
    )
    x = np.arange(len(counts))
    high_counts = counts.get("high-BDT background", pd.Series(0, index=counts.index))
    mid_counts = counts.get("middling true photon", pd.Series(0, index=counts.index))
    ax_counts.bar(x - 0.16, high_counts, width=0.32, color="#D95F02", label="high-BDT background")
    ax_counts.bar(x + 0.16, mid_counts, width=0.32, color="#2F6FED", label="middling true photon")
    ax_counts.set_xticks(x)
    ax_counts.set_xticklabels(["Jet20", "Jet30", "Jet40"])
    ax_counts.set_ylabel("saved candidates")
    ax_counts.set_title("Why Jet30/40 explain more", fontsize=13, fontweight="bold")
    ax_counts.grid(True, axis="y", color="#E8E8E8")
    ax_counts.legend(frameon=False, fontsize=8.5)
    ax_counts.annotate(
        "hard neutral-meson fakes\nin the 15-35 GeV window",
        xy=(1.5, max(high_counts.max(), 1)),
        xytext=(0.55, max(high_counts.max(), 1) * 0.72),
        arrowprops=dict(arrowstyle="->", color="#333333", lw=1.0),
        fontsize=9.2,
    )

    ax_map = fig.add_subplot(gs[6:, :8])
    for sample, marker in [("Jet20 inclusive", "o"), ("Jet30 inclusive", "s"), ("Jet40 inclusive", "^")]:
        sub = high[high["sample"] == sample]
        ax_map.scatter(
            sub["auau_tight_bdt_score"],
            sub["final_state_non_photon_frac"],
            s=34 + 1.4 * sub["cluster_Et"],
            alpha=0.78,
            marker=marker,
            edgecolor="white",
            linewidth=0.7,
            label=f"{sample.replace(' inclusive', '')} bkg n={len(sub)}",
        )
    ax_map.scatter(
        mid["auau_tight_bdt_score"],
        mid["final_state_non_photon_frac"],
        s=10,
        alpha=0.22,
        color="#2F6FED",
        label=f"middling true photons n={len(mid)}",
    )
    ax_map.axvline(0.8, color="#222222", linestyle=(0, (4, 3)), linewidth=1.0)
    ax_map.set_xlabel("BDT score")
    ax_map.set_ylabel("final-state non-photon pT fraction")
    ax_map.set_title("What the BDT is doing", fontsize=13, fontweight="bold")
    ax_map.set_xlim(0.43, 0.96)
    ax_map.set_ylim(-0.04, 1.04)
    ax_map.grid(True, color="#E8E8E8")
    ax_map.legend(frameon=False, fontsize=8.6, loc="upper left")

    ax_notes = fig.add_subplot(gs[6:, 8:])
    ax_notes.axis("off")
    inc_group = (
        final_particles[
            (final_particles["label"] == "high-BDT background")
            & (final_particles["sample"].str.contains("Jet", na=False))
        ]
        .groupby("group")["pt"]
        .sum()
    )
    inc_frac = inc_group / inc_group.sum()
    draw_note(
        ax_notes,
        0.02,
        0.68,
        0.96,
        0.26,
        "Particle answer",
        (
            f"In inclusive high-BDT backgrounds, final-state pT is "
            f"{100 * inc_frac.get('neutral_meson', 0):.1f}% neutral meson and "
            f"{100 * inc_frac.get('charged_hadron', 0):.1f}% charged hadron. "
            f"Only {100 * inc_frac.get('photon', 0):.1f}% is photon."
        ),
        face="#FFF3EA",
        edge="#E3B28D",
        body_width=47,
    )
    draw_note(
        ax_notes,
        0.02,
        0.38,
        0.96,
        0.24,
        "BDT answer",
        (
            "The BDT keeps true photons photon-core dominated, but high score is not truth. "
            "Pi0/eta fragments can look EM-like and become high-score false positives."
        ),
        face="#F2F6FF",
        edge="#B8C7E6",
        body_width=48,
    )
    draw_note(
        ax_notes,
        0.02,
        0.08,
        0.96,
        0.24,
        "Why Jet30/Jet40 matter",
        (
            "Jet20 is sparse. Jet30/Jet40 populate the hard neutral-meson phase space, "
            "so they reveal the false-positive taxonomy cleanly."
        ),
        face="#EEF8F4",
        edge="#A7CEC0",
        body_width=48,
    )

    fig.savefig(out)
    plt.close(fig)
    return out


def make_hepmc_trace_slide(candidates: pd.DataFrame, particles: pd.DataFrame, final_particles: pd.DataFrame, outdir: Path) -> Path:
    out = outdir / "the44_hepmc_trace_slide.png"
    inc_mask = lambda df: (df["label"] == "high-BDT background") & (df["sample"].str.contains("Jet", na=False))
    inc_final = final_particles[inc_mask(final_particles)].copy()
    inc_full = particles[inc_mask(particles)].copy()

    fig = plt.figure(figsize=(16, 9), dpi=170)
    fig.patch.set_facecolor("white")
    gs = GridSpec(12, 16, figure=fig, left=0.07, right=0.97, top=0.86, bottom=0.08, wspace=1.25, hspace=1.35)

    fig.text(0.055, 0.95, "HepMC Trace: Why the High-BDT Jet Background Is Physical", fontsize=24, fontweight="bold", family="serif")
    fig.text(
        0.055,
        0.912,
        "Final-state particles describe the visible neighborhood; full HepMC parents reveal the hard-process context",
        fontsize=14.2,
        color="#333333",
    )

    ax_pdg = fig.add_subplot(gs[:6, :6])
    pdg_pt = inc_final.groupby("pdg_name")["pt"].sum().sort_values(ascending=False).head(10)
    pdg_frac = (pdg_pt / inc_final["pt"].sum()).sort_values()
    y = np.arange(len(pdg_frac))
    bar_colors = ["#D95F02" if name in {"pi0", "eta"} else "#1B9E77" if name in {"pi+", "pi-", "K+", "K-"} else "#2F6FED" if name == "gamma" else "#777777" for name in pdg_frac.index]
    ax_pdg.barh(y, 100 * pdg_frac.values, color=bar_colors)
    ax_pdg.set_yticks(y)
    ax_pdg.set_yticklabels(pdg_frac.index, fontsize=10)
    ax_pdg.set_xlabel("share of final-state pT (%)")
    ax_pdg.set_title("Visible particles near inclusive high-BDT fakes", fontsize=13, fontweight="bold")
    ax_pdg.grid(True, axis="x", color="#E8E8E8")

    ax_parent = fig.add_subplot(gs[:6, 6:11])
    parent_pt = inc_final.groupby("parent_bucket")["pt"].sum().sort_values(ascending=False)
    parent_frac = (parent_pt / inc_final["pt"].sum()).sort_values()
    parent_colors = {
        "quark/gluon parent": "#7B6D8D",
        "meson-resonance parent": "#D95F02",
        "photon parent": "#2F6FED",
        "electron parent": "#6A3D9A",
        "other / unlisted parent": "#777777",
    }
    ax_parent.barh(parent_frac.index, 100 * parent_frac.values, color=[parent_colors.get(i, "#777777") for i in parent_frac.index])
    ax_parent.set_xlabel("share of final-state pT (%)")
    ax_parent.set_title("Immediate HepMC parent context", fontsize=13, fontweight="bold")
    ax_parent.grid(True, axis="x", color="#E8E8E8")

    ax_full = fig.add_subplot(gs[:6, 11:])
    full_group = inc_full.groupby("group")["pt"].sum()
    full_frac = full_group / full_group.sum()
    visible_group = inc_final.groupby("group")["pt"].sum()
    visible_frac = visible_group / visible_group.sum()
    compare = pd.DataFrame({"full record": full_frac, "final state": visible_frac}).fillna(0.0)
    compare = compare.reindex(["quark_gluon", "neutral_meson", "charged_hadron", "photon", "other"]).fillna(0.0)
    x = np.arange(len(compare))
    ax_full.bar(x - 0.17, 100 * compare["full record"], width=0.34, color="#B7A7C7", label="full HepMC record")
    ax_full.bar(x + 0.17, 100 * compare["final state"], width=0.34, color="#D95F02", label="final state")
    ax_full.set_xticks(x)
    ax_full.set_xticklabels([GROUP_LABELS.get(i, i) for i in compare.index], rotation=28, ha="right", fontsize=9)
    ax_full.set_ylabel("pT share (%)")
    ax_full.set_title("Do not mix ancestry with visible particles", fontsize=13, fontweight="bold")
    ax_full.grid(True, axis="y", color="#E8E8E8")
    ax_full.legend(frameon=False, fontsize=8.5)

    ax_flow = fig.add_subplot(gs[7:, :10])
    high = candidates[candidates["label"] == "high-BDT background"].copy()
    inc_high = high[high["sample"].str.contains("Jet", na=False)].copy()
    source = inc_high["truth_pid_name"].map(source_bucket).value_counts()
    source = source.reindex(["neutral meson", "photon", "charged hadron", "heavy flavor", "other"]).fillna(0)
    colors = ["#D95F02", "#2F6FED", "#1B9E77", "#9467BD", "#777777"]
    wedges, _ = ax_flow.pie(
        source.values,
        colors=colors,
        startangle=100,
        wedgeprops={"linewidth": 1, "edgecolor": "white"},
    )
    labels = [f"{idx}\n{int(val)}" for idx, val in source.items() if val > 0]
    ax_flow.legend(wedges, [f"{idx}: {int(val)}" for idx, val in source.items()], frameon=False, loc="center left", bbox_to_anchor=(0.88, 0.5), fontsize=10)
    ax_flow.set_title("Truth contributor type for inclusive high-BDT candidates", fontsize=13, fontweight="bold")
    ax_flow.text(0, -1.25, "Most false positives start as neutral-meson EM fragments", ha="center", fontsize=11, fontweight="bold")

    ax_notes = fig.add_subplot(gs[7:, 10:])
    ax_notes.axis("off")
    draw_note(
        ax_notes,
        0.02,
        0.66,
        0.96,
        0.27,
        "Trace in one sentence",
        (
            "Jet30/Jet40 create hard pi0/eta-rich fragments in the candidate pT window. "
            "They can look photon-like, even when the local final state is mostly non-photon."
        ),
        face="#FFF3EA",
        edge="#E3B28D",
        body_width=47,
    )
    draw_note(
        ax_notes,
        0.02,
        0.39,
        0.96,
        0.22,
        "How to read full HepMC",
        (
            "The full record is parton-heavy because it includes ancestry. "
            "Use status-1 particles for visible composition."
        ),
        face="#F2F6FF",
        edge="#B8C7E6",
        body_width=47,
    )
    draw_note(
        ax_notes,
        0.02,
        0.13,
        0.96,
        0.20,
        "BDT implication",
        (
            "The false positives are structured: neutral-meson and charged-hadron jet fragments "
            "mimic compact EM photon candidates."
        ),
        face="#EEF8F4",
        edge="#A7CEC0",
        body_width=49,
    )

    fig.savefig(out)
    plt.close(fig)
    return out


def write_particle_story_note(candidates: pd.DataFrame, particles: pd.DataFrame, final_particles: pd.DataFrame, outdir: Path) -> Path:
    out = outdir / "the44_particle_intensive_story.md"
    high = candidates[candidates["label"] == "high-BDT background"].copy()
    mid = candidates[candidates["label"] == "middling true photon"].copy()
    inc_high = high[high["sample"].str.contains("Jet", na=False)].copy()
    inc_final = final_particles[
        (final_particles["label"] == "high-BDT background")
        & (final_particles["sample"].str.contains("Jet", na=False))
    ].copy()
    inc_group_frac = inc_final.groupby("group")["pt"].sum()
    inc_group_frac = inc_group_frac / inc_group_frac.sum()
    truth_counts = inc_high["truth_pid_name"].value_counts()
    counts = (
        candidates[candidates["sample"].isin(INCLUSIVE_ORDER)]
        .groupby(["sample", "label"])["candidate_uid"]
        .count()
        .unstack(fill_value=0)
        .reindex(INCLUSIVE_ORDER)
    )
    top_final = inc_final.groupby("pdg_name")["pt"].sum().sort_values(ascending=False)
    top_final_frac = top_final / top_final.sum()
    high_particle_nonphoton = int((high["final_state_photon_frac"] < 0.20).sum())
    high_particle_photonish = int((high["final_state_photon_frac"] >= 0.20).sum())

    lines = [
        "# THE-44 Particle-Intensive BDT Autopsy",
        "",
        "## Simple answer",
        "",
        (
            "The high-BDT inclusive-background population is dominated by neutral-meson jet fragments, "
            "especially pi0, with substantial charged-hadron activity nearby. The BDT gives these rows "
            "high scores because the reconstructed shower can look photon-like even when the local "
            "generator final state is not photon dominated."
        ),
        "",
        "## Particle taxonomy",
        "",
        f"- Inclusive high-BDT background candidates: {len(inc_high)}.",
        f"- Truth contributors: pi0 {truth_counts.get('pi0', 0)}, eta {truth_counts.get('eta', 0)}, gamma {truth_counts.get('gamma', 0)}.",
        (
            f"- Final-state pT near inclusive high-BDT candidates: "
            f"neutral meson {100 * inc_group_frac.get('neutral_meson', 0):.1f}%, "
            f"charged hadron {100 * inc_group_frac.get('charged_hadron', 0):.1f}%, "
            f"photon {100 * inc_group_frac.get('photon', 0):.1f}%."
        ),
        (
            f"- Leading visible particle pT shares: pi0 {100 * top_final_frac.get('pi0', 0):.1f}%, "
            f"pi+ {100 * top_final_frac.get('pi+', 0):.1f}%, "
            f"gamma {100 * top_final_frac.get('gamma', 0):.1f}%, "
            f"pi- {100 * top_final_frac.get('pi-', 0):.1f}%, "
            f"eta {100 * top_final_frac.get('eta', 0):.1f}%."
        ),
        "",
        "## What the BDT is doing",
        "",
        (
            f"The middling true-photon population has {len(mid)} candidates and remains photon-core dominated "
            f"with mean final-state photon pT fraction {100 * mid['final_state_photon_frac'].mean():.1f}%. "
            f"Among all saved high-BDT background candidates, {high_particle_nonphoton} have less than 20% "
            f"final-state photon pT nearby, while {high_particle_photonish} are photon-like or mixed at particle level. "
            "The high-BDT inclusive false positives are therefore not an unstructured failure mode; the meeting concern "
            "is the lower-right population where the model score is photon-like but the visible particle neighborhood is not."
        ),
        "",
        "## Why the higher-energy inclusive samples explain more",
        "",
        (
            f"Saved high-BDT backgrounds by inclusive sample: Jet20 {int(counts.loc['Jet20 inclusive'].get('high-BDT background', 0))}, "
            f"Jet30 {int(counts.loc['Jet30 inclusive'].get('high-BDT background', 0))}, "
            f"Jet40 {int(counts.loc['Jet40 inclusive'].get('high-BDT background', 0))}. "
            "Jet20 is too sparse for a stable taxonomy, while Jet30 and Jet40 populate the hard neutral-meson-rich "
            "phase space that overlaps the 15-35 GeV photon-candidate window."
        ),
        "",
        "## Caveat",
        "",
        (
            "This is a qualitative diagnostic sample. The final-state particle fractions are the right view for "
            "visible local composition. The full HepMC record is parton-heavy by construction and should be used "
            "as ancestry/origin context, not as the visible-particle composition."
        ),
        "",
    ]
    out.write_text("\n".join(lines), encoding="utf-8")
    return out


def make_minimal_particle_summary(candidates: pd.DataFrame, outdir: Path) -> Path:
    out = outdir / "the44_minimal_particle_summary.png"
    high = candidates[candidates["label"] == "high-BDT background"].copy()
    mid = candidates[candidates["label"] == "middling true photon"].copy()
    high["source_bucket"] = high["truth_pid_name"].map(source_bucket)
    high["source_bucket"] = high["source_bucket"].replace(
        {
            "heavy flavor": "other",
            "charged hadron": "other",
        }
    )

    fig, (ax_map, ax_bar) = plt.subplots(
        1,
        2,
        figsize=(15.5, 6.8),
        dpi=180,
        gridspec_kw={"width_ratios": [1.42, 1.0], "wspace": 0.28},
    )
    fig.patch.set_facecolor("white")
    fig.subplots_adjust(left=0.075, right=0.98, bottom=0.12, top=0.80, wspace=0.28)
    fig.suptitle(
        "THE-44 BDT Particle Autopsy in Two Views",
        x=0.05,
        y=0.965,
        ha="left",
        fontsize=20,
        fontweight="bold",
        family="serif",
    )
    fig.text(
        0.05,
        0.895,
        "x-axis is the BDT decision; y-axis is whether the nearby final-state particle pT is photon-like or not",
        fontsize=11.5,
        color="#333333",
    )

    ax_map.scatter(
        mid["auau_tight_bdt_score"],
        mid["final_state_non_photon_frac"],
        s=13,
        color="#2F6FED",
        alpha=0.25,
        edgecolors="none",
        label=f"true photons, middling score (n={len(mid)})",
    )
    colors = {
        "neutral meson": "#D95F02",
        "photon": "#7B3FB2",
        "electron": "#1B9E77",
        "other": "#666666",
    }
    labels = {
        "neutral meson": "high-BDT bkg: pi0/eta",
        "photon": "high-BDT bkg: gamma",
        "electron": "high-BDT bkg: e+/e-",
        "other": "high-BDT bkg: other",
    }
    for bucket in ["neutral meson", "photon", "electron", "other"]:
        sub = high[high["source_bucket"] == bucket]
        if sub.empty:
            continue
        ax_map.scatter(
            sub["auau_tight_bdt_score"],
            sub["final_state_non_photon_frac"],
            s=40,
            color=colors[bucket],
            alpha=0.82,
            edgecolor="white",
            linewidth=0.6,
            label=f"{labels[bucket]} (n={len(sub)})",
        )
    ax_map.axvline(0.8, color="#222222", linestyle=(0, (4, 3)), linewidth=1.0)
    ax_map.axhspan(0.0, 0.15, color="#2F6FED", alpha=0.06)
    ax_map.axhspan(0.75, 1.02, color="#D95F02", alpha=0.06)
    ax_map.text(0.455, 0.105, "photon-core region", color="#1F4EA8", fontsize=10.5, fontweight="bold")
    ax_map.text(0.805, 0.94, "non-photon-rich high-score fakes", color="#8A3A00", fontsize=10.5, fontweight="bold")
    ax_map.set_xlim(0.43, 0.96)
    ax_map.set_ylim(-0.035, 1.035)
    ax_map.set_xlabel("BDT score")
    ax_map.set_ylabel("final-state non-photon pT fraction")
    ax_map.set_title("1. What the BDT score means physically", fontsize=13.5, fontweight="bold")
    ax_map.grid(True, color="#E8E8E8", linewidth=0.8)
    ax_map.legend(frameon=False, fontsize=8.4, loc="center left", bbox_to_anchor=(0.01, 0.55))

    inc = high[high["sample"].isin(INCLUSIVE_ORDER)].copy()
    inc["sample_short"] = inc["sample"].str.replace(" inclusive", "", regex=False)
    inc["source_bucket"] = inc["source_bucket"].replace({"electron": "other"})
    stack = (
        inc.groupby(["sample_short", "source_bucket"])["candidate_uid"]
        .count()
        .unstack(fill_value=0)
        .reindex(["Jet20", "Jet30", "Jet40"])
        .fillna(0)
    )
    bar_groups = ["neutral meson", "photon", "other"]
    x = np.arange(len(stack))
    bottom = np.zeros(len(stack))
    for bucket in bar_groups:
        vals = stack[bucket].to_numpy() if bucket in stack else np.zeros(len(stack))
        ax_bar.bar(
            x,
            vals,
            bottom=bottom,
            color=colors.get(bucket, "#666666"),
            width=0.58,
            label={"neutral meson": "pi0/eta", "photon": "gamma", "other": "other truth"}[bucket],
        )
        bottom += vals
    for xi, total in zip(x, bottom):
        ax_bar.text(xi, total + 2.0, f"{int(total)}", ha="center", va="bottom", fontsize=12, fontweight="bold")
    ax_bar.set_xticks(x)
    ax_bar.set_xticklabels(stack.index)
    ax_bar.set_ylabel("high-BDT background candidates")
    ax_bar.set_title("2. Why Jet30/Jet40 explain the background", fontsize=13.5, fontweight="bold")
    ax_bar.grid(True, axis="y", color="#E8E8E8", linewidth=0.8)
    ax_bar.legend(frameon=False, fontsize=9.5, loc="upper left")
    ax_bar.set_ylim(0, max(bottom) * 1.22 if len(bottom) else 1)
    fig.savefig(out)
    plt.close(fig)
    return out


def make_minimal_hepmc_trace_ladder(candidates: pd.DataFrame, particles: pd.DataFrame, final_particles: pd.DataFrame, outdir: Path) -> Path:
    out = outdir / "the44_minimal_hepmc_trace_ladder.png"
    inc_high = candidates[
        (candidates["label"] == "high-BDT background")
        & (candidates["sample"].str.contains("Jet", na=False))
    ].copy()
    inc_particles = particles[
        (particles["label"] == "high-BDT background")
        & (particles["sample"].str.contains("Jet", na=False))
    ].copy()
    inc_final = final_particles[
        (final_particles["label"] == "high-BDT background")
        & (final_particles["sample"].str.contains("Jet", na=False))
    ].copy()

    truth = inc_high["truth_pid_name"].map(source_bucket).replace({"heavy flavor": "other", "charged hadron": "other", "electron": "other"})
    truth_share = truth.value_counts(normalize=True).reindex(["neutral meson", "photon", "other"]).fillna(0.0)

    final_pt = inc_final.groupby("group")["pt"].sum()
    final_share = (final_pt / final_pt.sum()).reindex(["neutral_meson", "charged_hadron", "photon", "other"]).fillna(0.0)

    full_pt = inc_particles.groupby("group")["pt"].sum()
    full_share = (full_pt / full_pt.sum()).reindex(["quark_gluon", "neutral_meson", "charged_hadron", "photon", "other"]).fillna(0.0)

    rows = [
        ("candidate truth source", truth_share, {"neutral meson": "#D95F02", "photon": "#7B3FB2", "other": "#666666"}),
        ("visible status-1 particles", final_share, GROUP_COLORS),
        ("full HepMC record", full_share, GROUP_COLORS),
    ]
    label_map = {
        "neutral meson": "pi0/eta",
        "neutral_meson": "neutral meson",
        "charged_hadron": "charged hadron",
        "quark_gluon": "quark/gluon",
        "photon": "photon",
        "other": "other",
    }

    fig, ax = plt.subplots(figsize=(12.5, 5.7), dpi=180)
    fig.patch.set_facecolor("white")
    y_positions = np.arange(len(rows))[::-1]
    for y, (row_label, shares, color_map) in zip(y_positions, rows):
        left = 0.0
        for key, frac in shares.items():
            if frac <= 0:
                continue
            color = color_map.get(key, "#777777")
            ax.barh(y, 100 * frac, left=left, height=0.54, color=color, edgecolor="white", linewidth=1.0)
            if frac >= 0.055:
                ax.text(
                    left + 50 * frac,
                    y,
                    f"{label_map.get(key, key)}\n{100 * frac:.0f}%",
                    ha="center",
                    va="center",
                    fontsize=9.2,
                    color="white" if frac > 0.13 else "#111111",
                    fontweight="bold" if frac > 0.13 else "normal",
                )
            left += 100 * frac

    ax.set_yticks(y_positions)
    ax.set_yticklabels([r[0] for r in rows], fontsize=11.5)
    ax.set_xlim(0, 100)
    ax.set_xlabel("share (%)")
    ax.set_title("Inclusive High-BDT Background: One-Step HepMC Trace", fontsize=17, fontweight="bold", family="serif", pad=14)
    fig.text(
        0.19,
        0.085,
        "Truth source -> visible status-1 particles -> full HepMC ancestry.",
        fontsize=10.2,
        color="#333333",
    )
    fig.text(
        0.19,
        0.055,
        "Use the middle row for visible composition; the full record is parton-heavy by construction.",
        fontsize=10.2,
        color="#333333",
    )
    ax.grid(True, axis="x", color="#E8E8E8", linewidth=0.8)
    for spine in ["top", "right", "left"]:
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis="y", length=0)
    fig.subplots_adjust(left=0.19, right=0.98, top=0.84, bottom=0.22)
    fig.savefig(out)
    plt.close(fig)
    return out


def make_verdict_plane(candidates: pd.DataFrame, outdir: Path) -> Path:
    out = outdir / "the44_verdict_plane_bdt_vs_pythia.png"
    high = candidates[candidates["label"] == "high-BDT background"].copy()
    mid = candidates[candidates["label"] == "middling true photon"].copy()
    high["particle_verdict"] = np.where(
        high["final_state_photon_frac"] < 0.20,
        "particle non-photon",
        "particle photon-like/mixed",
    )

    colors = {
        "particle non-photon": "#D95F02",
        "particle photon-like/mixed": "#7B3FB2",
    }
    labels = {
        "particle non-photon": "high-BDT background: Pythia non-photon",
        "particle photon-like/mixed": "high-BDT background: Pythia photon-like/mixed",
    }

    fig, ax = plt.subplots(figsize=(11.6, 7.0), dpi=180)
    fig.patch.set_facecolor("white")
    ax.axvspan(0.8, 0.965, color="#F7E4D1", alpha=0.55, zorder=0)
    ax.axhspan(0.0, 0.20, color="#F7E4D1", alpha=0.48, zorder=0)
    ax.axhspan(0.86, 1.02, color="#E8F0FF", alpha=0.62, zorder=0)
    ax.axvline(0.8, color="#222222", linestyle=(0, (4, 3)), linewidth=1.1)

    ax.scatter(
        mid["auau_tight_bdt_score"],
        mid["final_state_photon_frac"],
        s=16,
        color="#2F6FED",
        alpha=0.27,
        edgecolors="none",
        label=f"middling true photons (n={len(mid)})",
    )
    for bucket in ["particle photon-like/mixed", "particle non-photon"]:
        sub = high[high["particle_verdict"] == bucket]
        if sub.empty:
            continue
        ax.scatter(
            sub["auau_tight_bdt_score"],
            sub["final_state_photon_frac"],
            s=44,
            color=colors[bucket],
            alpha=0.86,
            edgecolor="white",
            linewidth=0.7,
            label=f"{labels[bucket]} (n={len(sub)})",
        )

    ax.annotate(
        "false-positive\nmechanism",
        xy=(0.885, 0.025),
        xytext=(0.815, 0.15),
        color="#8A3A00",
        fontsize=13,
        fontweight="bold",
        arrowprops=dict(arrowstyle="-|>", color="#8A3A00", lw=1.2, shrinkA=4, shrinkB=4),
    )
    ax.text(0.465, 0.93, "true-photon control:\nPythia photon-core", color="#1F4EA8", fontsize=13, fontweight="bold")
    ax.text(0.823, 0.72, "BDT high,\nparticles still photon-rich", color="#4B2380", fontsize=11.8, fontweight="bold")
    ax.text(0.804, 1.015, "BDT > 0.8", color="#222222", fontsize=10.5, va="bottom")
    ax.set_xlim(0.43, 0.96)
    ax.set_ylim(-0.03, 1.04)
    ax.set_xlabel("BDT score: model says photon-like ->")
    ax.set_ylabel("Pythia status-1 photon pT fraction: particles say photon-like ->")
    ax.set_title("BDT Verdict vs Particle Verdict", fontsize=18, fontweight="bold", family="serif", pad=14)
    ax.grid(True, color="#E7E7E7", linewidth=0.8)
    ax.legend(frameon=False, fontsize=9.3, loc="center left", bbox_to_anchor=(0.02, 0.47))
    fig.text(
        0.08,
        0.045,
        "Read the lower-right cloud as the concern from the meeting: high BDT score, but nearby final-state particle pT is not photon dominated.",
        fontsize=10.6,
        color="#333333",
    )
    fig.subplots_adjust(left=0.11, right=0.97, top=0.86, bottom=0.13)
    fig.savefig(out)
    plt.close(fig)
    return out


def make_fake_factory_ladder(candidates: pd.DataFrame, particles: pd.DataFrame, final_particles: pd.DataFrame, outdir: Path) -> Path:
    out = outdir / "the44_fake_factory_ladder.png"
    inc_high = candidates[
        (candidates["label"] == "high-BDT background")
        & (candidates["sample"].str.contains("Jet", na=False))
    ].copy()
    inc_particles = particles[
        (particles["label"] == "high-BDT background")
        & (particles["sample"].str.contains("Jet", na=False))
    ].copy()
    inc_final = final_particles[
        (final_particles["label"] == "high-BDT background")
        & (final_particles["sample"].str.contains("Jet", na=False))
    ].copy()

    sample_share = (
        inc_high["sample"]
        .str.replace(" inclusive", "", regex=False)
        .value_counts(normalize=True)
        .reindex(["Jet20", "Jet30", "Jet40"])
        .fillna(0.0)
    )
    truth_share = (
        inc_high["truth_pid_name"]
        .map(source_bucket)
        .replace({"heavy flavor": "other", "charged hadron": "other", "electron": "other"})
        .value_counts(normalize=True)
        .reindex(["neutral meson", "photon", "other"])
        .fillna(0.0)
    )
    visible_pt = inc_final.groupby("group")["pt"].sum()
    visible_share = (
        visible_pt / visible_pt.sum()
    ).reindex(["neutral_meson", "charged_hadron", "photon", "other"]).fillna(0.0)
    full_pt = inc_particles.groupby("group")["pt"].sum()
    full_share = (
        full_pt / full_pt.sum()
    ).reindex(["quark_gluon", "neutral_meson", "charged_hadron", "photon", "other"]).fillna(0.0)

    rows = [
        ("which inclusive sample?", sample_share, {"Jet20": "#A6CEE3", "Jet30": "#1F78B4", "Jet40": "#08306B"}, "candidate count share"),
        ("candidate truth source", truth_share, {"neutral meson": "#D95F02", "photon": "#7B3FB2", "other": "#666666"}, "candidate count share"),
        ("visible status-1 particles", visible_share, GROUP_COLORS, "nearby pT share"),
        ("full HepMC record", full_share, GROUP_COLORS, "nearby pT share"),
    ]
    label_map = {
        "Jet20": "Jet20",
        "Jet30": "Jet30",
        "Jet40": "Jet40",
        "neutral meson": "pi0/eta",
        "neutral_meson": "neutral meson",
        "charged_hadron": "charged hadron",
        "quark_gluon": "quark/gluon",
        "photon": "photon",
        "other": "other",
    }

    fig, ax = plt.subplots(figsize=(12.8, 6.6), dpi=180)
    fig.patch.set_facecolor("white")
    y_positions = np.arange(len(rows))[::-1]
    for y, (row_label, shares, color_map, _) in zip(y_positions, rows):
        left = 0.0
        for key, frac in shares.items():
            if frac <= 0:
                continue
            color = color_map.get(key, "#777777")
            ax.barh(y, 100 * frac, left=left, height=0.55, color=color, edgecolor="white", linewidth=1.0)
            if frac >= 0.055:
                ax.text(
                    left + 50 * frac,
                    y,
                    f"{label_map.get(key, key)}\n{100 * frac:.0f}%",
                    ha="center",
                    va="center",
                    fontsize=9.2,
                    color="white" if frac > 0.13 and key not in {"Jet20"} else "#111111",
                    fontweight="bold" if frac > 0.13 else "normal",
                )
            left += 100 * frac
    ax.set_yticks(y_positions)
    ax.set_yticklabels([r[0] for r in rows], fontsize=11.2)
    ax.set_xlim(0, 100)
    ax.set_xlabel("share (%)")
    ax.set_title("Inclusive High-BDT Fake Factory", fontsize=18, fontweight="bold", family="serif", pad=15)
    ax.grid(True, axis="x", color="#E8E8E8", linewidth=0.8)
    for spine in ["top", "right", "left"]:
        ax.spines[spine].set_visible(False)
    ax.tick_params(axis="y", length=0)
    fig.text(
        0.21,
        0.075,
        "Read top to bottom. Jet30/40 supply the statistics; pi0/eta dominate truth source; status-1 particles are meson+hadron rich.",
        fontsize=10.4,
        color="#333333",
    )
    fig.text(
        0.21,
        0.045,
        "The full HepMC row is ancestry context, not visible composition.",
        fontsize=10.4,
        color="#333333",
    )
    fig.subplots_adjust(left=0.22, right=0.98, top=0.85, bottom=0.19)
    fig.savefig(out)
    plt.close(fig)
    return out


def make_fake_sorting_flow(candidates: pd.DataFrame, outdir: Path) -> Path:
    out = outdir / "the44_high_bdt_fake_sorting_flow.png"
    inc_high = candidates[
        (candidates["label"] == "high-BDT background")
        & (candidates["sample"].str.contains("Jet", na=False))
    ].copy()
    inc_high["sample_simple"] = inc_high["sample"].str.replace(" inclusive", "", regex=False)
    inc_high["truth_flow"] = inc_high["truth_pid_name"].map(source_bucket).replace(
        {"neutral meson": "pi0/eta truth", "photon": "gamma truth", "heavy flavor": "other truth", "charged hadron": "other truth", "electron": "other truth", "other": "other truth"}
    )
    inc_high["particle_verdict"] = np.where(
        inc_high["final_state_photon_frac"] < 0.20,
        "particle non-photon",
        "particle photon-like/mixed",
    )

    total = len(inc_high)
    sample_order = ["Jet20", "Jet30", "Jet40"]
    pthat_ranges = {
        "Jet20": "pThat 21-32",
        "Jet30": "pThat 32-42",
        "Jet40": "pThat 42-100",
    }
    sample_energy = (
        inc_high.groupby("sample_simple")["cluster_Et"]
        .agg(["median", lambda s: s.quantile(0.25), lambda s: s.quantile(0.75)])
        .rename(columns={"<lambda_0>": "q25", "<lambda_1>": "q75"})
        .reindex(sample_order)
    )
    truth_order = ["pi0/eta truth", "gamma truth", "other truth"]
    verdict_order = ["particle non-photon", "particle photon-like/mixed"]
    colors = {
        "Jet20": "#A6CEE3",
        "Jet30": "#1F78B4",
        "Jet40": "#08306B",
        "pi0/eta truth": "#D95F02",
        "gamma truth": "#7B3FB2",
        "other truth": "#666666",
        "particle non-photon": "#D95F02",
        "particle photon-like/mixed": "#7B3FB2",
    }
    column_counts = {
        "sample": inc_high["sample_simple"].value_counts().reindex(sample_order).fillna(0).astype(int),
        "truth": inc_high["truth_flow"].value_counts().reindex(truth_order).fillna(0).astype(int),
        "verdict": inc_high["particle_verdict"].value_counts().reindex(verdict_order).fillna(0).astype(int),
    }
    flow_sample_truth = (
        inc_high.groupby(["sample_simple", "truth_flow"]).size().reindex(
            pd.MultiIndex.from_product([sample_order, truth_order]), fill_value=0
        )
    )
    flow_truth_verdict = (
        inc_high.groupby(["truth_flow", "particle_verdict"]).size().reindex(
            pd.MultiIndex.from_product([truth_order, verdict_order]), fill_value=0
        )
    )

    def layout(counts: pd.Series, order: list[str], top: float = 0.80, bottom: float = 0.20, gap: float = 0.035) -> dict[str, tuple[float, float]]:
        usable = top - bottom - gap * (len(order) - 1)
        y = top
        intervals = {}
        for key in order:
            h = usable * (counts.get(key, 0) / total)
            intervals[key] = (y - h, y)
            y -= h + gap
        return intervals

    def take_segment(segments: dict[str, float], key: str, n: int, scale: float) -> tuple[float, float]:
        y0 = segments[key]
        y1 = y0 + n * scale
        segments[key] = y1
        return y0, y1

    def draw_ribbon(ax, x0: float, y0a: float, y0b: float, x1: float, y1a: float, y1b: float, color: str, alpha: float = 0.42) -> None:
        curve = 0.18
        verts = [
            (x0, y0a),
            (x0 + curve, y0a),
            (x1 - curve, y1a),
            (x1, y1a),
            (x1, y1b),
            (x1 - curve, y1b),
            (x0 + curve, y0b),
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
        ax.add_patch(PathPatch(MplPath(verts, codes), facecolor=color, edgecolor="none", alpha=alpha, zorder=1))

    def draw_node(
        ax,
        x: float,
        w: float,
        interval: tuple[float, float],
        key: str,
        count: int,
        label: str,
        fontsize: float = 10.2,
        linespacing: float = 0.95,
    ) -> None:
        y0, y1 = interval
        ax.add_patch(
            FancyBboxPatch(
                (x, y0),
                w,
                y1 - y0,
                boxstyle="round,pad=0.006,rounding_size=0.01",
                facecolor=colors.get(key, "#777777"),
                edgecolor="white",
                linewidth=1.2,
                zorder=3,
            )
        )
        text_color = "white" if key not in {"Jet20"} else "#111111"
        ax.text(
            x + w / 2,
            (y0 + y1) / 2,
            f"{label}\n{count} ({100 * count / total:.0f}%)",
            ha="center",
            va="center",
            fontsize=fontsize,
            linespacing=linespacing,
            color=text_color,
            fontweight="bold" if count / total > 0.12 else "normal",
            zorder=4,
        )

    sample_layout = layout(column_counts["sample"], sample_order)
    truth_layout = layout(column_counts["truth"], truth_order)
    verdict_layout = layout(column_counts["verdict"], verdict_order)
    sample_scale = {k: (sample_layout[k][1] - sample_layout[k][0]) / max(column_counts["sample"].get(k, 0), 1) for k in sample_order}
    truth_in_scale = {k: (truth_layout[k][1] - truth_layout[k][0]) / max(column_counts["truth"].get(k, 0), 1) for k in truth_order}
    truth_out_scale = truth_in_scale.copy()
    verdict_scale = {k: (verdict_layout[k][1] - verdict_layout[k][0]) / max(column_counts["verdict"].get(k, 0), 1) for k in verdict_order}

    sample_cursor = {k: sample_layout[k][0] for k in sample_order}
    truth_in_cursor = {k: truth_layout[k][0] for k in truth_order}
    truth_out_cursor = {k: truth_layout[k][0] for k in truth_order}
    verdict_cursor = {k: verdict_layout[k][0] for k in verdict_order}

    fig, ax = plt.subplots(figsize=(13.2, 7.0), dpi=180)
    fig.patch.set_facecolor("white")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    x_sample, x_truth, x_verdict = 0.055, 0.42, 0.76
    w_sample, w_truth, w_verdict = 0.19, 0.18, 0.21
    for sample in sample_order:
        for truth in truth_order:
            n = int(flow_sample_truth.loc[(sample, truth)])
            if n <= 0:
                continue
            y0a, y0b = take_segment(sample_cursor, sample, n, sample_scale[sample])
            y1a, y1b = take_segment(truth_in_cursor, truth, n, truth_in_scale[truth])
            draw_ribbon(ax, x_sample + w_sample, y0a, y0b, x_truth, y1a, y1b, colors[truth], alpha=0.34)

    for truth in truth_order:
        for verdict in verdict_order:
            n = int(flow_truth_verdict.loc[(truth, verdict)])
            if n <= 0:
                continue
            y0a, y0b = take_segment(truth_out_cursor, truth, n, truth_out_scale[truth])
            y1a, y1b = take_segment(verdict_cursor, verdict, n, verdict_scale[verdict])
            draw_ribbon(ax, x_truth + w_truth, y0a, y0b, x_verdict, y1a, y1b, colors[truth], alpha=0.44)

    sample_label_map = {}
    for key in sample_order:
        row = sample_energy.loc[key]
        if pd.isna(row["median"]):
            sample_label_map[key] = f"{key}\n{pthat_ranges[key]} GeV"
        else:
            sample_label_map[key] = (
                f"{key}\n{pthat_ranges[key]} GeV"
                f"\ncluster E_T {row['median']:.1f}"
            )

    label_map = {
        "pi0/eta truth": "pi0/eta",
        "gamma truth": "gamma",
        "other truth": "other",
        "particle non-photon": "non-photon\nparticles",
        "particle photon-like/mixed": "photon-like\nor mixed",
    }
    for key in sample_order:
        draw_node(
            ax,
            x_sample,
            w_sample,
            sample_layout[key],
            key,
            int(column_counts["sample"].get(key, 0)),
            sample_label_map[key],
            fontsize=8.6,
            linespacing=0.92,
        )
    for key in truth_order:
        draw_node(ax, x_truth, w_truth, truth_layout[key], key, int(column_counts["truth"].get(key, 0)), label_map[key])
    for key in verdict_order:
        draw_node(ax, x_verdict, w_verdict, verdict_layout[key], key, int(column_counts["verdict"].get(key, 0)), label_map[key])

    ax.text(x_sample + w_sample / 2, 0.87, "which Pythia jet-energy bin?", ha="center", va="center", fontsize=11.2, fontweight="bold")
    ax.text(x_truth + w_truth / 2, 0.87, "what particle made the candidate?", ha="center", va="center", fontsize=11.2, fontweight="bold")
    ax.text(x_verdict + w_verdict / 2, 0.87, "what did nearby particles say?", ha="center", va="center", fontsize=11.2, fontweight="bold")
    ax.text(0.04, 0.975, "High-BDT Fake Sorting Flow", ha="left", va="top", fontsize=19, fontweight="bold", family="serif")
    ax.text(0.04, 0.920, "Jet energy -> truth source -> particle-neighborhood verdict for inclusive backgrounds with BDT > 0.8", ha="left", va="top", fontsize=10.8, color="#333333")
    ax.text(0.04, 0.062, "Ribbon width is candidate count. Jet labels are generator pThat bins; cluster E_T is the reconstructed fake-candidate median.", ha="left", va="center", fontsize=9.9, color="#333333")
    ax.text(0.04, 0.032, "Key read: Jet30+Jet40 = 146/160 fakes; Jet40 lifts the median to 21.1 GeV; pi0/eta -> non-photon particles is the main high-BDT path.", ha="left", va="center", fontsize=9.9, color="#333333")
    fig.savefig(out)
    plt.close(fig)
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, default=DEFAULT_BASE, help="Full THE-44 slideReady package directory.")
    parser.add_argument("--outdir", type=Path, default=None, help="Output directory; defaults to <base>/particle_intensive.")
    args = parser.parse_args()

    outdir = args.outdir or args.base / "particle_intensive"
    outdir.mkdir(parents=True, exist_ok=True)
    candidates, particles, final_particles = load_tables(args.base)

    taxonomy = make_particle_taxonomy_slide(candidates, final_particles, outdir)
    trace = make_hepmc_trace_slide(candidates, particles, final_particles, outdir)
    minimal = make_minimal_particle_summary(candidates, outdir)
    trace_ladder = make_minimal_hepmc_trace_ladder(candidates, particles, final_particles, outdir)
    verdict = make_verdict_plane(candidates, outdir)
    factory = make_fake_factory_ladder(candidates, particles, final_particles, outdir)
    flow = make_fake_sorting_flow(candidates, outdir)
    note = write_particle_story_note(candidates, particles, final_particles, outdir)

    manifest = {
        "status": "READY",
        "inputs": {
            "base": str(args.base),
            "candidate_rows": str(args.base / "the44_candidate_autopsy_rows.csv"),
            "particle_rows": str(args.base / "the44_particle_rows.csv"),
        },
        "outputs": {
            "particle_taxonomy_slide": str(taxonomy),
            "hepmc_trace_slide": str(trace),
            "minimal_particle_summary": str(minimal),
            "minimal_hepmc_trace_ladder": str(trace_ladder),
            "verdict_plane_bdt_vs_pythia": str(verdict),
            "fake_factory_ladder": str(factory),
            "fake_sorting_flow": str(flow),
            "particle_story_note": str(note),
        },
        "n_candidates": int(len(candidates)),
        "n_particles": int(len(particles)),
    }
    manifest_path = outdir / "the44_particle_intensive_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
