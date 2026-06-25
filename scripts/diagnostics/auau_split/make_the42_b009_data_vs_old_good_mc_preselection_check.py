#!/usr/bin/env python3
"""Diagnostic AuAu data-vs-older-good-MC E11/E33 side check.

This intentionally avoids tight-BDT conclusions.  It overlays b009 AuAu data
with older embedded MC that already passed the E11/E33 reconstruction-contract
sanity check, using only the common pre-BDT stages: all candidates and NPB pass.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import uproot


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")

DATA_ROOT = REPO / (
    "dataOutput/auauTightBDTValidation/THE42_wp80_centlinear_ss_20260604/"
    "b009_merged_data_roots_20260609/"
    "RecoilJets_auau_ALL_preselectionNewPPG12_tightAuAuCentInputBase3x3BDT_"
    "nonTightAuAuBDTComplement_baseVariant.root"
)
OLD_SIGNAL_ROOT = REPO / (
    "dataOutput/auau_widthstudy_pt1530_wp080/combinedSimOnlyEMBEDDED/"
    "preselectionNewPPG12_tightAuauCentInputBase3x3BDT_nonTightAuAuBDTComplement_baseVariant/"
    "photonJet12and20merged_SIM/RecoilJets_embeddedPhoton12plus20_MERGED.root"
)
OLD_INCLUSIVE_ROOT = REPO / (
    "dataOutput/auau_widthstudy_pt1530_wp080/combinedSimOnlyEMBEDDED/"
    "preselectionNewPPG12_tightAuauCentInputBase3x3BDT_nonTightAuAuBDTComplement_baseVariant/"
    "embeddedJet12and20merged_SIM/RecoilJets_embeddedJet12plus20_MERGED.root"
)

OUT_DIR = REPO / "dataOutput/auauTableQA/the42_b009_data_old_good_mc_preselection_check"
OUT_PNG = OUT_DIR / "b009_photon10_data_vs_old_good_mc_e11e33_preselection_check.png"
OUT_JSON = OUT_DIR / "b009_photon10_data_vs_old_good_mc_e11e33_preselection_check_manifest.json"

DATA_DIR = "photon_10_plus_MBD_NS_geq_2_vtx_lt_150"
PT_TOKENS = ("18_20", "20_22", "22_24", "24_26")
STAGES = (
    ("inclusive", "Before preselection"),
    ("npbPass", "After NPB preselection"),
)
CENTRALITIES = (
    ("0-20%", ("0_10", "10_20"), "0_20"),
    ("20-50%", ("20_30", "30_40", "40_50"), "20_50"),
    ("50-80%", ("50_60", "60_80"), "50_80"),
)


def sum_data_hist(path: Path, stage: str, fine_cent_tokens: tuple[str, ...]) -> tuple[np.ndarray, np.ndarray, list[str]]:
    values = None
    edges = None
    used: list[str] = []
    with uproot.open(path) as root_file:
        data_dir = root_file[DATA_DIR]
        for pt in PT_TOKENS:
            for cent in fine_cent_tokens:
                name = f"h_ss_e11e33_{stage}_pT_{pt}_cent_{cent}"
                if name not in data_dir:
                    continue
                hist_values, hist_edges = data_dir[name].to_numpy(flow=False)
                hist_values = np.asarray(hist_values, dtype=float)
                hist_edges = np.asarray(hist_edges, dtype=float)
                if values is None:
                    values = hist_values.copy()
                    edges = hist_edges.copy()
                else:
                    if values.shape != hist_values.shape or not np.allclose(edges, hist_edges):
                        raise ValueError(f"data binning mismatch for {name}")
                    values += hist_values
                used.append(f"{DATA_DIR}/{name}")
    if values is None or edges is None:
        raise KeyError(f"no data histograms found for stage={stage} cent={fine_cent_tokens}")
    return values, edges, used


def sum_old_mc_hist(path: Path, stage: str, coarse_cent_token: str) -> tuple[np.ndarray, np.ndarray, list[str]]:
    values = None
    edges = None
    used: list[str] = []
    with uproot.open(path) as root_file:
        sim_dir = root_file["SIM"]
        for pt in PT_TOKENS:
            name = f"h_ss_e11e33_{stage}_pT_{pt}_cent_{coarse_cent_token}"
            if name not in sim_dir:
                continue
            hist_values, hist_edges = sim_dir[name].to_numpy(flow=False)
            hist_values = np.asarray(hist_values, dtype=float)
            hist_edges = np.asarray(hist_edges, dtype=float)
            if values is None:
                values = hist_values.copy()
                edges = hist_edges.copy()
            else:
                if values.shape != hist_values.shape or not np.allclose(edges, hist_edges):
                    raise ValueError(f"MC binning mismatch for {name}")
                values += hist_values
            used.append(f"SIM/{name}")
    if values is None or edges is None:
        raise KeyError(f"no MC histograms found for {path} stage={stage} cent={coarse_cent_token}")
    return values, edges, used


def normalize(values: np.ndarray, edges: np.ndarray, xlim: tuple[float, float] = (0.0, 1.0)) -> tuple[np.ndarray, np.ndarray]:
    centers = 0.5 * (edges[:-1] + edges[1:])
    mask = (centers >= xlim[0]) & (centers <= xlim[1])
    out = values.astype(float).copy()
    total = float(np.sum(out[mask]))
    if total > 0:
        out /= total
    return out, mask


def metrics(values: np.ndarray, edges: np.ndarray) -> dict[str, float]:
    centers = 0.5 * (edges[:-1] + edges[1:])
    total = float(np.sum(values))
    sumw2 = float(np.sum(values * values))
    return {
        "sumw": total,
        "low_edge_fraction_x_lt_0p05": float(np.sum(values[centers < 0.05]) / total) if total else 0.0,
        "max_bin_fraction": float(np.max(values) / total) if total else 0.0,
        "effective_entries": float(total * total / sumw2) if sumw2 > 0 else 0.0,
        "mean": float(np.sum(values * centers) / total) if total else 0.0,
    }


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 11,
            "axes.linewidth": 1.1,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )

    fig, axes = plt.subplots(3, 2, figsize=(13.8, 10.4), constrained_layout=False, sharex=True)
    manifest: dict[str, object] = {
        "purpose": "diagnostic only: b009 AuAu data versus older-good embedded MC before tight BDT",
        "data_root": str(DATA_ROOT),
        "old_signal_root": str(OLD_SIGNAL_ROOT),
        "old_inclusive_root": str(OLD_INCLUSIVE_ROOT),
        "data_trigger_directory": DATA_DIR,
        "pt_tokens_summed": list(PT_TOKENS),
        "stages": [s for s, _ in STAGES],
        "centralities": {},
        "caveat": "Not a final THE-42 product and not a tight-BDT validation; data uses the photon_10_plus_MBD_NS_geq_2_vtx_lt_150 analysis trigger directory, old MC uses photon/jet 12+20 only, and the overlay is restricted to the common 18-26 GeV overlap.",
    }

    for row, (cent_label, data_cents, mc_cent) in enumerate(CENTRALITIES):
        manifest["centralities"][cent_label] = {}
        for col, (stage, stage_label) in enumerate(STAGES):
            ax = axes[row, col]
            data_values, edges, data_used = sum_data_hist(DATA_ROOT, stage, data_cents)
            sig_values, sig_edges, sig_used = sum_old_mc_hist(OLD_SIGNAL_ROOT, stage, mc_cent)
            inc_values, inc_edges, inc_used = sum_old_mc_hist(OLD_INCLUSIVE_ROOT, stage, mc_cent)
            if not (np.allclose(edges, sig_edges) and np.allclose(edges, inc_edges)):
                raise ValueError(f"binning mismatch for {cent_label} {stage}")

            centers = 0.5 * (edges[:-1] + edges[1:])
            data_norm, mask = normalize(data_values, edges)
            sig_norm, _ = normalize(sig_values, edges)
            inc_norm, _ = normalize(inc_values, edges)

            ax.step(edges[:-1], sig_norm, where="post", color="#d62728", linewidth=2.0, label="older signal MC")
            ax.step(edges[:-1], inc_norm, where="post", color="#1f77b4", linewidth=2.0, label="older inclusive MC")
            ax.errorbar(
                centers[mask],
                data_norm[mask],
                yerr=np.sqrt(np.maximum(data_values[mask], 0)) / max(np.sum(data_values[mask]), 1.0),
                fmt="o",
                ms=4.8,
                color="black",
                markerfacecolor="white",
                markeredgewidth=1.2,
                linewidth=0,
                capsize=0,
                label="b009 data",
            )
            ax.set_xlim(0, 1.0)
            ymax = max(float(np.max(data_norm[mask])), float(np.max(sig_norm[mask])), float(np.max(inc_norm[mask]))) * 1.22
            ax.set_ylim(0, max(0.055, ymax))
            ax.grid(True, axis="y", color="#d7dce2", linewidth=0.8, alpha=0.85)
            ax.set_title(f"{cent_label}  {stage_label}", loc="left", fontsize=12.5, fontweight="bold", pad=7)
            if row == 2:
                ax.set_xlabel(r"$E_{11}/E_{33}$")
            if col == 0:
                ax.set_ylabel("normalized counts")
            ax.text(
                0.035,
                0.955,
                f"data x<0.05 {metrics(data_values, edges)['low_edge_fraction_x_lt_0p05']:.3f}\n"
                f"sig x<0.05 {metrics(sig_values, edges)['low_edge_fraction_x_lt_0p05']:.3f}\n"
                f"inc x<0.05 {metrics(inc_values, edges)['low_edge_fraction_x_lt_0p05']:.3f}",
                transform=ax.transAxes,
                va="top",
                ha="left",
                fontsize=9.2,
                bbox={"boxstyle": "round,pad=0.22", "facecolor": "white", "edgecolor": "#c7c7c7", "alpha": 0.94},
            )
            if row == 0 and col == 1:
                ax.legend(loc="upper right", frameon=False, fontsize=10)

            manifest["centralities"][cent_label][stage] = {
                "data_histograms": data_used,
                "old_signal_histograms": sig_used,
                "old_inclusive_histograms": inc_used,
                "data_metrics": metrics(data_values, edges),
                "old_signal_metrics": metrics(sig_values, edges),
                "old_inclusive_metrics": metrics(inc_values, edges),
            }

    fig.suptitle("AuAu E11/E33 side check before tight BDT", fontsize=19, fontweight="bold", y=0.985)
    fig.text(
        0.5,
        0.948,
        "Photon10+MBD scaled-trigger data vs older embedded MC with validated shower-shape contract",
        ha="center",
        va="top",
        fontsize=12.5,
    )
    fig.text(
        0.5,
        0.924,
        "Common 18-26 GeV pT overlap only; columns stop before BDT-dependent interpretation.",
        ha="center",
        va="top",
        fontsize=11.4,
    )
    fig.text(
        0.5,
        0.018,
        "Diagnostic only: checks data shower-shape coherence before BDT-dependent interpretation. It is not the final THE-57/THE-42 baseline product.",
        ha="center",
        va="bottom",
        fontsize=10.8,
        color="#4d4d4d",
    )
    fig.subplots_adjust(left=0.07, right=0.985, top=0.885, bottom=0.07, hspace=0.25, wspace=0.16)
    fig.savefig(OUT_PNG, dpi=190)
    plt.close(fig)
    OUT_JSON.write_text(json.dumps(manifest, indent=2) + "\n")
    print(OUT_PNG)
    print(OUT_JSON)


if __name__ == "__main__":
    main()
