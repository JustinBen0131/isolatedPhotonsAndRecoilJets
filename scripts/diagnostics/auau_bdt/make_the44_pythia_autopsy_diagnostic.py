#!/usr/bin/env python3
"""Build local THE-44 BDT/Pythia autopsy tables and slide diagnostics."""

from __future__ import annotations

import argparse
import csv
import json
import math
import textwrap
from collections import Counter, defaultdict
from pathlib import Path

import awkward as ak
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import uproot
from matplotlib.gridspec import GridSpec
from matplotlib.patches import FancyBboxPatch


DEFAULT_BASE = Path(
    "dataOutput/auauTightBDTValidation/"
    "THE8_branchA_ladder_scorecache_fullstat_20260527/"
    "diagnostics/the44_high_bdt_pythia_autopsy_20260605"
)

PARTICLE_GROUPS = [
    ("photon", "Photon", "#2F6FED"),
    ("neutral_meson", "Neutral meson", "#D95F02"),
    ("charged_hadron", "Charged hadron", "#1B9E77"),
    ("quark_gluon", "Quark/gluon", "#6A3D9A"),
    ("lepton", "Lepton", "#7570B3"),
    ("other", "Other", "#7F7F7F"),
]


def pdg_name(pid: int) -> str:
    names = {
        22: "gamma",
        111: "pi0",
        221: "eta",
        331: "eta'",
        211: "pi+",
        -211: "pi-",
        321: "K+",
        -321: "K-",
        311: "K0",
        -311: "K0bar",
        130: "K0L",
        310: "K0S",
        113: "rho0",
        213: "rho+",
        -213: "rho-",
        223: "omega",
        323: "K*+",
        -323: "K*-",
        411: "D+",
        -411: "D-",
        421: "D0",
        -421: "D0bar",
        431: "Ds+",
        -431: "Ds-",
        21: "g",
        1: "d",
        -1: "dbar",
        2: "u",
        -2: "ubar",
        3: "s",
        -3: "sbar",
        4: "c",
        -4: "cbar",
        5: "b",
        -5: "bbar",
        11: "e-",
        -11: "e+",
        13: "mu-",
        -13: "mu+",
        2212: "p",
        -2212: "pbar",
        2112: "n",
        -2112: "nbar",
    }
    return names.get(int(pid), str(int(pid)))


def particle_group(pid: int) -> str:
    apid = abs(int(pid))
    if int(pid) == 22:
        return "photon"
    if apid in {111, 221, 331, 113, 213, 223, 311, 130, 310, 313, 323}:
        return "neutral_meson"
    if apid in {211, 321, 411, 421, 431, 2212, 2112, 3122, 3222, 3312, 3334}:
        return "charged_hadron"
    if apid in {1, 2, 3, 4, 5, 6, 21}:
        return "quark_gluon"
    if apid in {11, 13, 15}:
        return "lepton"
    return "other"


def category_label(category: int, is_signal: int) -> str:
    if category == 1 and not is_signal:
        return "high-BDT background"
    if category == 2 and is_signal:
        return "middling true photon"
    if category == 1:
        return "high-BDT selected"
    if category == 2:
        return "middling selected"
    return "other"


def sample_label(path: Path) -> str:
    name = path.name
    if "embeddedJet40" in name:
        return "Jet40 inclusive"
    if "embeddedPhoton20" in name:
        return "Photon20 signal"
    if "embeddedPhoton12" in name:
        return "Photon12 signal"
    if "embeddedJet30" in name:
        return "Jet30 inclusive"
    if "embeddedJet20" in name:
        return "Jet20 inclusive"
    if "embeddedJet12" in name:
        return "Jet12 inclusive"
    return path.stem[:28]


def as_list(value) -> list:
    if value is None:
        return []
    try:
        return ak.to_list(value)
    except Exception:
        return list(value)


def load_roots(root_paths: list[Path]) -> tuple[pd.DataFrame, pd.DataFrame]:
    candidate_rows = []
    particle_rows = []
    scalar_branches = [
        "run",
        "evt",
        "category",
        "is_signal",
        "pt_bin",
        "cent_bin",
        "source_sample_code",
        "cluster_Et",
        "cluster_Eta",
        "cluster_Phi",
        "centrality",
        "reco_eiso",
        "auau_tight_bdt_score",
        "npb_score",
        "cluster_weta_cogx",
        "cluster_wphi_cogx",
        "e11_over_e33",
        "e32_over_e35",
        "cluster_truth_track_id",
        "cluster_truth_pid",
        "cluster_truth_barcode",
        "cluster_truth_econtrib",
        "matched_signal_barcode",
        "matched_signal_pt",
        "matched_signal_eta",
        "matched_signal_phi",
        "matched_signal_iso",
    ]
    vector_branches = [
        "part_pdg",
        "part_status",
        "part_barcode",
        "part_parent0_pdg",
        "part_parent1_pdg",
        "part_pt",
        "part_eta",
        "part_phi",
        "part_dr",
    ]

    for root_path in root_paths:
        with uproot.open(root_path) as root_file:
            tree = root_file["THE44PythiaAutopsyTree"]
            arrays = tree.arrays(scalar_branches + vector_branches, library="ak")
            sample = sample_label(root_path)
            for i in range(tree.num_entries):
                pdgs = [int(x) for x in as_list(arrays["part_pdg"][i])]
                pts = [float(x) for x in as_list(arrays["part_pt"][i])]
                drs = [float(x) for x in as_list(arrays["part_dr"][i])]
                parent0 = [int(x) for x in as_list(arrays["part_parent0_pdg"][i])]
                parent1 = [int(x) for x in as_list(arrays["part_parent1_pdg"][i])]
                statuses = [int(x) for x in as_list(arrays["part_status"][i])]
                barcodes = [int(x) for x in as_list(arrays["part_barcode"][i])]

                group_pt = defaultdict(float)
                group_n = Counter()
                final_group_pt = defaultdict(float)
                final_group_n = Counter()
                top_particles = []
                for j, pid in enumerate(pdgs):
                    pt = pts[j] if j < len(pts) else float("nan")
                    group = particle_group(pid)
                    status = statuses[j] if j < len(statuses) else None
                    if math.isfinite(pt):
                        group_pt[group] += pt
                        if status == 1:
                            final_group_pt[group] += pt
                    group_n[group] += 1
                    if status == 1:
                        final_group_n[group] += 1
                    if len(top_particles) < 8:
                        top_particles.append(pdg_name(pid))
                    particle_rows.append(
                        {
                            "candidate_uid": f"{root_path.stem}:{i}",
                            "sample": sample,
                            "row": i,
                            "pdg": pid,
                            "pdg_name": pdg_name(pid),
                            "group": group,
                            "status": statuses[j] if j < len(statuses) else None,
                            "barcode": barcodes[j] if j < len(barcodes) else None,
                            "parent0_pdg": parent0[j] if j < len(parent0) else None,
                            "parent1_pdg": parent1[j] if j < len(parent1) else None,
                            "pt": pt,
                            "dr": drs[j] if j < len(drs) else float("nan"),
                        }
                    )

                total_pt = sum(group_pt.values())
                non_photon_pt = total_pt - group_pt["photon"]
                final_total_pt = sum(final_group_pt.values())
                final_non_photon_pt = final_total_pt - final_group_pt["photon"]
                scalar = {name: arrays[name][i].item() for name in scalar_branches}
                cat = int(scalar["category"])
                sig = int(scalar["is_signal"])
                row = {
                    "candidate_uid": f"{root_path.stem}:{i}",
                    "sample": sample,
                    "root_file": root_path.name,
                    "row": i,
                    **scalar,
                    "label": category_label(cat, sig),
                    "truth_pid_name": pdg_name(int(scalar["cluster_truth_pid"])),
                    "n_particles": len(pdgs),
                    "top_particles": " ".join(top_particles),
                    "total_particle_pt": total_pt,
                    "non_photon_particle_pt": non_photon_pt,
                    "non_photon_pt_fraction": non_photon_pt / total_pt if total_pt > 0 else 0.0,
                    "photon_pt_fraction": group_pt["photon"] / total_pt if total_pt > 0 else 0.0,
                    "final_state_total_pt": final_total_pt,
                    "final_state_non_photon_pt": final_non_photon_pt,
                    "final_state_non_photon_frac": final_non_photon_pt / final_total_pt if final_total_pt > 0 else 0.0,
                    "final_state_photon_frac": final_group_pt["photon"] / final_total_pt if final_total_pt > 0 else 0.0,
                    "final_state_n_particles": sum(final_group_n.values()),
                }
                for group, _, _ in PARTICLE_GROUPS:
                    row[f"{group}_pt"] = group_pt[group]
                    row[f"{group}_n"] = group_n[group]
                    row[f"{group}_frac"] = group_pt[group] / total_pt if total_pt > 0 else 0.0
                    row[f"final_{group}_pt"] = final_group_pt[group]
                    row[f"final_{group}_n"] = final_group_n[group]
                    row[f"final_{group}_frac"] = final_group_pt[group] / final_total_pt if final_total_pt > 0 else 0.0
                candidate_rows.append(row)

    return pd.DataFrame(candidate_rows), pd.DataFrame(particle_rows)


def write_tables(candidates: pd.DataFrame, particles: pd.DataFrame, outdir: Path) -> dict:
    outdir.mkdir(parents=True, exist_ok=True)
    candidates_path = outdir / "the44_candidate_autopsy_rows.csv"
    particles_path = outdir / "the44_particle_rows.csv"
    summary_path = outdir / "the44_category_summary.csv"
    examples_path = outdir / "the44_event_card_examples.csv"
    manifest_path = outdir / "the44_pythia_autopsy_manifest.json"

    candidates.to_csv(candidates_path, index=False, quoting=csv.QUOTE_MINIMAL)
    particles.to_csv(particles_path, index=False, quoting=csv.QUOTE_MINIMAL)

    summary = (
        candidates.groupby(["sample", "label"], dropna=False)
        .agg(
            entries=("candidate_uid", "count"),
            mean_bdt=("auau_tight_bdt_score", "mean"),
            median_bdt=("auau_tight_bdt_score", "median"),
            mean_cluster_et=("cluster_Et", "mean"),
            mean_non_photon_frac=("non_photon_pt_fraction", "mean"),
            mean_photon_frac=("photon_pt_fraction", "mean"),
            mean_final_non_photon_frac=("final_state_non_photon_frac", "mean"),
            mean_final_photon_frac=("final_state_photon_frac", "mean"),
            mean_particles=("n_particles", "mean"),
            mean_final_particles=("final_state_n_particles", "mean"),
        )
        .reset_index()
    )
    summary.to_csv(summary_path, index=False)

    examples = select_examples(candidates)
    examples.to_csv(examples_path, index=False)

    manifest = {
        "status": "READY",
        "root_files": sorted(candidates["root_file"].unique().tolist()),
        "n_candidates": int(len(candidates)),
        "n_particles": int(len(particles)),
        "category_counts": candidates["label"].value_counts().to_dict(),
        "final_state_group_fractions": group_fraction_summary(candidates, prefix="final_"),
        "full_record_group_fractions": group_fraction_summary(candidates, prefix=""),
        "outputs": {
            "candidate_rows": str(candidates_path),
            "particle_rows": str(particles_path),
            "category_summary": str(summary_path),
            "event_examples": str(examples_path),
        },
        "interpretation_scope": (
            "Qualitative local autopsy of the THE-44 six-sample fanout "
            "(Photon12/20 plus Jet12/20/30/40); not a weighted production-yield estimate."
        ),
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    return {
        "candidates": candidates_path,
        "particles": particles_path,
        "summary": summary_path,
        "examples": examples_path,
        "manifest": manifest_path,
    }


def group_fraction_summary(candidates: pd.DataFrame, prefix: str) -> dict:
    out = {}
    total_col = "final_state_total_pt" if prefix == "final_" else "total_particle_pt"
    for (sample, label), sub in candidates.groupby(["sample", "label"]):
        total = float(sub[total_col].sum())
        key = f"{sample} | {label}"
        out[key] = {}
        for group, group_label, _ in PARTICLE_GROUPS:
            col = f"{prefix}{group}_pt" if prefix else f"{group}_pt"
            out[key][group_label] = float(sub[col].sum() / total) if total > 0 else 0.0
    return out


def select_examples(candidates: pd.DataFrame) -> pd.DataFrame:
    examples = []
    seen = set()

    def add(label: str, frame: pd.DataFrame, n: int, sort_cols: list[str], ascending: list[bool]):
        if frame.empty:
            return
        for _, row in frame.sort_values(sort_cols, ascending=ascending).head(n * 3).iterrows():
            uid = row["candidate_uid"]
            if uid in seen:
                continue
            copy = row.copy()
            copy["example_role"] = label
            examples.append(copy)
            seen.add(uid)
            if sum(1 for ex in examples if ex["example_role"] == label) >= n:
                break

    high_bkg = candidates[candidates["label"] == "high-BDT background"].copy()
    mid_sig = candidates[candidates["label"] == "middling true photon"].copy()
    if not high_bkg.empty:
        high_bkg["sample_preference"] = np.where(high_bkg["sample"].str.contains("Jet"), 0, 1)
    add(
        "highest-score non-photon background",
        high_bkg,
        3,
        ["sample_preference", "auau_tight_bdt_score", "final_state_non_photon_frac", "final_state_n_particles"],
        [True, False, False, False],
    )
    add(
        "busy true photon with lowered score",
        mid_sig,
        3,
        ["final_state_non_photon_frac", "final_state_n_particles", "auau_tight_bdt_score"],
        [False, False, True],
    )
    if not examples:
        return candidates.head(0)
    return pd.DataFrame(examples)


def card(ax, x, y, w, h, title, body, face="#F7F7F7", edge="#D0D0D0"):
    box = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.012,rounding_size=0.015",
        transform=ax.transAxes,
        facecolor=face,
        edgecolor=edge,
        linewidth=1.0,
    )
    ax.add_patch(box)
    ax.text(x + 0.025, y + h - 0.055, title, transform=ax.transAxes, fontsize=13.5, fontweight="bold", va="top")
    wrapped = "\n".join(textwrap.fill(line, width=42) for line in body.split("\n"))
    ax.text(x + 0.025, y + h - 0.12, wrapped, transform=ax.transAxes, fontsize=10.1, va="top", linespacing=1.18)


def make_raw_plot(candidates: pd.DataFrame, outdir: Path) -> Path:
    out = outdir / "the44_score_vs_particle_fingerprint.png"
    fig, ax = plt.subplots(figsize=(9.6, 6.2), dpi=180)
    colors = {
        "high-BDT background": "#D95F02",
        "middling true photon": "#2F6FED",
        "high-BDT selected": "#A6761D",
        "middling selected": "#66A61E",
    }
    for label, sub in candidates.groupby("label"):
        ax.scatter(
            sub["auau_tight_bdt_score"],
            sub["final_state_non_photon_frac"],
            s=35 + 3.2 * sub["final_state_n_particles"],
            alpha=0.76,
            color=colors.get(label, "#666666"),
            edgecolor="white",
            linewidth=0.7,
            label=f"{label} (n={len(sub)})",
        )
    ax.axvline(0.8, color="#333333", linestyle="--", linewidth=1.0)
    ax.set_xlabel("Au+Au tight BDT score")
    ax.set_ylabel("Final-state non-photon pT fraction")
    ax.set_title("THE-44 candidate score joined to same-event Pythia final-state particles")
    ax.grid(True, color="#E6E6E6", linewidth=0.8)
    ax.legend(frameon=False, loc="upper left")
    ax.text(
        0.99,
        0.02,
        "Marker size scales with final-state particles inside Delta R < 0.3",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=9,
        color="#555555",
    )
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)
    return out


def make_slide(candidates: pd.DataFrame, outdir: Path) -> Path:
    out = outdir / "the44_bdt_pythia_autopsy_slide.png"
    examples = select_examples(candidates)
    fig = plt.figure(figsize=(16, 9), dpi=160)
    fig.patch.set_facecolor("white")
    gs = GridSpec(12, 16, figure=fig, left=0.12, right=0.97, top=0.865, bottom=0.09, wspace=0.82, hspace=1.15)

    fig.text(0.055, 0.945, "BDT Autopsy: The Score Is Now Connected to Particles", fontsize=26, fontweight="bold", family="serif")
    fig.text(
        0.055,
        0.908,
        "Each suspicious photon candidate now carries its BDT score, truth tag, and same-event Pythia/HepMC neighborhood",
        fontsize=14.5,
        color="#333333",
    )

    ax_scatter = fig.add_subplot(gs[:6, :9])
    colors = {"high-BDT background": "#D95F02", "middling true photon": "#2F6FED"}
    markers = {"high-BDT background": "o", "middling true photon": "s"}
    for label, sub in candidates.groupby("label"):
        ax_scatter.scatter(
            sub["auau_tight_bdt_score"],
            sub["final_state_non_photon_frac"],
            s=42 + 3.0 * sub["final_state_n_particles"],
            alpha=0.78,
            marker=markers.get(label, "o"),
            color=colors.get(label, "#777777"),
            edgecolor="white",
            linewidth=0.8,
            label=f"{label}  n={len(sub)}",
        )
    ax_scatter.axvline(0.8, color="#222222", linestyle=(0, (4, 3)), linewidth=1.0)
    ax_scatter.text(0.805, 0.97, "high-BDT cut", transform=ax_scatter.get_xaxis_transform(), fontsize=9, va="top")
    ax_scatter.set_xlabel("BDT score", fontsize=12)
    ax_scatter.set_ylabel("final-state non-photon pT fraction", fontsize=12)
    ax_scatter.set_title("Candidate score vs final-state particle fingerprint", fontsize=14, fontweight="bold", pad=10)
    ax_scatter.set_xlim(0.42, max(0.93, candidates["auau_tight_bdt_score"].max() + 0.03))
    ax_scatter.set_ylim(-0.04, 1.04)
    ax_scatter.grid(True, color="#E8E8E8", linewidth=0.8)
    ax_scatter.legend(frameon=False, loc="upper left", fontsize=10)
    ax_scatter.text(
        0.98,
        0.05,
        "size = final-state multiplicity",
        transform=ax_scatter.transAxes,
        ha="right",
        fontsize=9.5,
        color="#555555",
    )

    ax_bar = fig.add_subplot(gs[7:, :10])
    if not examples.empty:
        examples = examples.reset_index(drop=True).head(6)
        y = np.arange(len(examples))
        left = np.zeros(len(examples))
        for group, label, color in PARTICLE_GROUPS:
            vals = examples[f"final_{group}_frac"].to_numpy()
            ax_bar.barh(y, vals, left=left, height=0.62, label=label, color=color)
            left += vals
        ylabels = []
        for idx, row in examples.iterrows():
            role = "Bkg" if row["label"] == "high-BDT background" else "Sig"
            ylabels.append(
                f"{chr(65 + idx)} {role} | {row['auau_tight_bdt_score']:.2f} | {row['truth_pid_name']}"
            )
        ax_bar.set_yticks(y)
        ax_bar.set_yticklabels(ylabels, fontsize=9.7)
        ax_bar.invert_yaxis()
        ax_bar.set_xlim(0, 1)
        ax_bar.set_xlabel("fraction of nearby final-state pT", fontsize=11)
        ax_bar.set_title("Selected event fingerprints inside Delta R < 0.3", fontsize=14, fontweight="bold", pad=10)
        ax_bar.grid(True, axis="x", color="#EAEAEA")
        ax_bar.legend(ncol=3, frameon=False, fontsize=8.7, loc="lower center", bbox_to_anchor=(0.5, -0.34))
    else:
        ax_bar.text(0.5, 0.5, "No event examples selected", ha="center", va="center")
        ax_bar.axis("off")

    ax_cards = fig.add_subplot(gs[:, 10:])
    ax_cards.axis("off")
    total = len(candidates)
    count_text = candidates["label"].value_counts().to_dict()
    card(
        ax_cards,
        0.02,
        0.74,
        0.96,
        0.22,
        "What changed",
        "Before: score plots had no clean event-level particle explanation.\n"
        "Now: score, truth contributor, and Pythia particles live in one tree.",
        face="#F2F6FF",
        edge="#B8C7E6",
    )
    card(
        ax_cards,
        0.02,
        0.49,
        0.96,
        0.21,
        "First physics answer",
        "High-BDT backgrounds are inspectable EM-like jet fragments.\n"
        "Their final-state neighborhood can be tested for photon-like vs non-photon pT.",
        face="#FFF3EA",
        edge="#E3B28D",
    )
    card(
        ax_cards,
        0.02,
        0.25,
        0.96,
        0.20,
        "Why some true photons score lower",
        "Middling-score signal rows are usually photon-core dominated.\n"
        "The exceptions expose nearby activity that makes the BDT less confident.",
        face="#EEF8F4",
        edge="#A7CEC0",
    )
    card(
        ax_cards,
        0.02,
        0.03,
        0.96,
        0.19,
        "Scope",
        f"Probe: {total} candidates.\n"
        f"High-BDT background: {count_text.get('high-BDT background', 0)}.\n"
        f"Middling true photon: {count_text.get('middling true photon', 0)}.\n"
        "Not a weighted purity result.",
        face="#F8F8F8",
        edge="#D5D5D5",
    )

    fig.savefig(out)
    plt.close(fig)
    return out


def write_script(outdir: Path) -> Path:
    path = outdir / "the44_bdt_pythia_autopsy_script.md"
    path.write_text(
        "# WP GammaJets Slide Script - BDT Autopsy: The Score Is Now Connected to Particles\n\n"
        "The point of this slide is that the high-BDT background question is no longer just a score-shape question. "
        "For these diagnostic rows, I can take the exact candidate that received a BDT score and look at the Pythia "
        "particles around that same reconstructed cluster.\n\n"
        "On the top-left plot, each point is one saved candidate. The horizontal axis is the BDT score, and the "
        "vertical axis is how much of the nearby final-state particle transverse momentum is not carried by photons. The high-score "
        "background rows are not random. They are cases where a non-signal candidate has a photon-like shower score, "
        "but the particle neighborhood tells us that it lives inside a hard jet-like environment.\n\n"
        "The lower panel turns individual rows into fingerprints. Each bar is one candidate, and the colors show what "
        "kind of final-state particles are inside Delta R less than 0.3. The high-BDT background examples contain non-photon truth "
        "contributors and meson or charged-hadron nearby activity. The middling-score true photons are usually photon-core "
        "dominated, but some have enough nearby activity that the BDT has a real reason to be less confident.\n\n"
        "So the new strategy is simple: when a BDT region looks suspicious, we do not only ask what the score distribution "
        "looks like. We pick candidates from that region and inspect the particle neighborhood that produced the score. "
        "That gives us an event-level physics explanation for the BDT behavior.\n",
        encoding="utf-8",
    )
    return path


def fmt_frac(value: float) -> str:
    return f"{100.0 * value:.1f}%"


def mean_or_zero(frame: pd.DataFrame, column: str) -> float:
    if frame.empty:
        return 0.0
    return float(frame[column].mean())


def write_interpretation_note(candidates: pd.DataFrame, outdir: Path) -> Path:
    path = outdir / "the44_full_autopsy_interpretation.md"
    high_bkg = candidates[candidates["label"] == "high-BDT background"]
    mid_sig = candidates[candidates["label"] == "middling true photon"]
    jet_high = high_bkg[high_bkg["sample"].str.contains("Jet", na=False)]
    photon_high = high_bkg[high_bkg["sample"].str.contains("Photon", na=False)]

    lines = [
        "# THE-44 Full Pythia Particle-List Autopsy",
        "",
        "## What is new",
        "",
        (
            "The candidate-level BDT score is now joined to same-event Pythia/HepMC "
            "particles within Delta R < 0.3 around the reconstructed cluster. "
            "This turns a suspicious score region into an inspectable event-level "
            "particle fingerprint."
        ),
        "",
        "## Main diagnostic answer",
        "",
        (
            f"- Full probe: {len(candidates)} saved candidates from 48 ROOT files "
            "(Photon12/20 plus Jet12/20/30/40)."
        ),
        (
            f"- High-BDT background rows: {len(high_bkg)} candidates, mean BDT "
            f"{mean_or_zero(high_bkg, 'auau_tight_bdt_score'):.3f}."
        ),
        (
            f"- Inclusive-jet high-BDT backgrounds: {len(jet_high)} candidates, "
            f"mean final-state non-photon pT fraction "
            f"{fmt_frac(mean_or_zero(jet_high, 'final_state_non_photon_frac'))}."
        ),
        (
            f"- Middling true photons: {len(mid_sig)} candidates, mean final-state "
            f"photon pT fraction {fmt_frac(mean_or_zero(mid_sig, 'final_state_photon_frac'))}."
        ),
        (
            f"- Photon-sample high-BDT backgrounds are mixed: {len(photon_high)} "
            f"candidates, mean final-state photon pT fraction "
            f"{fmt_frac(mean_or_zero(photon_high, 'final_state_photon_frac'))}; "
            "these need per-event inspection rather than a one-number label."
        ),
        "",
        "## Practical reading",
        "",
        (
            "The inclusive-jet high-BDT backgrounds are not mysterious score artifacts. "
            "They are photon-like reconstructed clusters whose nearby generator "
            "neighborhood is dominated by non-photon final-state pT, especially neutral "
            "mesons and charged hadrons. The middling true photons usually remain "
            "photon-core dominated, and the exceptions expose nearby activity that "
            "gives the BDT a physical reason to lower confidence."
        ),
        "",
        "## Caveat",
        "",
        (
            "This is a qualitative diagnostic sample, not a weighted production-yield "
            "or final purity estimate. Use the final-state fractions for the cleanest "
            "particle-composition statement; the full HepMC record is useful as "
            "ancestry context because it includes generator partons."
        ),
        "",
        "## Next-use recipe",
        "",
        (
            "Pick any suspicious BDT score band, select representative candidates, "
            "then inspect their final-state particle fingerprint and event-card rows. "
            "This directly answers whether a score feature is photon-core-like, "
            "neutral-meson-like, charged-hadron-rich, or busy because of nearby activity."
        ),
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root-dir",
        type=Path,
        default=DEFAULT_BASE / "roots",
        help="Directory containing pulled THE-44 ROOT files.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=DEFAULT_BASE / "slideReady",
        help="Output directory for tables and PNGs.",
    )
    parser.add_argument("--roots", nargs="*", type=Path, default=None, help="Explicit ROOT files.")
    args = parser.parse_args()

    roots = args.roots or sorted(args.root_dir.glob("*.root"))
    if not roots:
        raise SystemExit(f"No ROOT files found under {args.root_dir}")
    candidates, particles = load_roots(roots)
    tables = write_tables(candidates, particles, args.outdir)
    raw_plot = make_raw_plot(candidates, args.outdir)
    slide = make_slide(candidates, args.outdir)
    script = write_script(args.outdir)
    note = write_interpretation_note(candidates, args.outdir)

    print(json.dumps({
        "status": "READY",
        "n_roots": len(roots),
        "n_candidates": int(len(candidates)),
        "n_particles": int(len(particles)),
        "outputs": {k: str(v) for k, v in tables.items()},
        "raw_plot": str(raw_plot),
        "slide": str(slide),
        "script": str(script),
        "interpretation": str(note),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
