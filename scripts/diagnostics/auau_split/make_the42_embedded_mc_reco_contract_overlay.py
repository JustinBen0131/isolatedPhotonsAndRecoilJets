#!/usr/bin/env python3
"""Compare old working vs current embedded-AuAu MC E11/E33 reconstruction shapes.

This is an intentionally narrow diagnostic for the THE-42 table-QA failure mode:
the same legacy histogram names are read from an older working embedded-MC ROOT
and the current table-QA embedded-MC ROOTs, then normalized and overlaid.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import uproot


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")

OLD_SIGNAL = REPO / (
    "dataOutput/auau_widthstudy_pt1530_wp080/combinedSimOnlyEMBEDDED/"
    "preselectionNewPPG12_tightAuauCentInputBase3x3BDT_nonTightAuAuBDTComplement_baseVariant/"
    "photonJet12and20merged_SIM/RecoilJets_embeddedPhoton12plus20_MERGED.root"
)
OLD_INCLUSIVE = REPO / (
    "dataOutput/auau_widthstudy_pt1530_wp080/combinedSimOnlyEMBEDDED/"
    "preselectionNewPPG12_tightAuauCentInputBase3x3BDT_nonTightAuAuBDTComplement_baseVariant/"
    "embeddedJet12and20merged_SIM/RecoilJets_embeddedJet12plus20_MERGED.root"
)
CURRENT_SIGNAL = REPO / "InputFiles/the42_current_auau_tableqa/RecoilJets_embeddedPhoton12plus20_MERGED.root"
CURRENT_INCLUSIVE = REPO / "InputFiles/the42_current_auau_tableqa/RecoilJets_embeddedJet12plus20plus30plus40_MERGED.root"

OUT_DIR = REPO / "dataOutput/auauTableQA/the42_current_auau_reco_contract_debug"
OUT_PNG = OUT_DIR / "auau_embedded_mc_e11e33_old_vs_current_contract_overlay.png"
OUT_JSON = OUT_DIR / "auau_embedded_mc_e11e33_old_vs_current_contract_overlay_manifest.json"

PT_TOKENS = ("18_20", "20_22", "22_24", "24_26")
ROWS = [
    ("Signal MC, inclusive row", "h_ss_e11e33_inclusive_pT_", OLD_SIGNAL, CURRENT_SIGNAL),
    ("Inclusive MC, inclusive row", "h_ss_e11e33_inclusive_pT_", OLD_INCLUSIVE, CURRENT_INCLUSIVE),
    ("Signal MC, non-tight signal row", "h_ss_e11e33_nonTight_sig_pT_", OLD_SIGNAL, CURRENT_SIGNAL),
    ("Inclusive MC, non-tight background row", "h_ss_e11e33_nonTight_bkg_pT_", OLD_INCLUSIVE, CURRENT_INCLUSIVE),
]


def _sum_hist(path: Path, prefix: str) -> tuple[np.ndarray, np.ndarray, list[str]]:
    values = None
    edges = None
    names = [f"{prefix}{pt}" for pt in PT_TOKENS]
    used: list[str] = []

    with uproot.open(path) as root_file:
        if "SIM" not in root_file:
            raise KeyError(f"{path} has no SIM directory")
        sim_dir = root_file["SIM"]
        for name in names:
            if name not in sim_dir:
                raise KeyError(f"{path} missing SIM/{name}")
            hist_values, hist_edges = sim_dir[name].to_numpy(flow=False)
            hist_values = np.asarray(hist_values, dtype=float)
            hist_edges = np.asarray(hist_edges, dtype=float)
            if values is None:
                values = hist_values.copy()
                edges = hist_edges.copy()
            else:
                if values.shape != hist_values.shape or not np.allclose(edges, hist_edges):
                    raise ValueError(f"binning mismatch for {path}: {name}")
                values += hist_values
            used.append(f"SIM/{name}")

    if values is None or edges is None:
        raise RuntimeError(f"no histograms found for {path} prefix={prefix}")
    return values, edges, used


def _metrics(values: np.ndarray, edges: np.ndarray) -> dict[str, float]:
    total = float(np.sum(values))
    centers = 0.5 * (edges[:-1] + edges[1:])
    sumw2 = float(np.sum(values * values))
    return {
        "sumw": total,
        "low_edge_fraction_x_lt_0p05": float(np.sum(values[centers < 0.05]) / total) if total else 0.0,
        "max_bin_fraction": float(np.max(values) / total) if total else 0.0,
        "effective_entries": float(total * total / sumw2) if sumw2 > 0 else 0.0,
        "mean": float(np.sum(values * centers) / total) if total else 0.0,
    }


def _normalize(values: np.ndarray) -> np.ndarray:
    total = np.sum(values)
    if total <= 0:
        return values
    return values / total


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 11,
            "axes.linewidth": 1.2,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )

    fig, axes = plt.subplots(2, 2, figsize=(15.5, 8.7), constrained_layout=False)
    axes = axes.ravel()
    manifest: dict[str, object] = {
        "purpose": "THE-42 embedded AuAu MC reconstruction-contract old-vs-current E11/E33 diagnostic",
        "old_signal": str(OLD_SIGNAL),
        "old_inclusive": str(OLD_INCLUSIVE),
        "current_signal": str(CURRENT_SIGNAL),
        "current_inclusive": str(CURRENT_INCLUSIVE),
        "pt_tokens_summed": list(PT_TOKENS),
        "rows": [],
    }

    for ax, (title, prefix, old_path, current_path) in zip(axes, ROWS):
        old_values, edges, old_used = _sum_hist(old_path, prefix)
        current_values, current_edges, current_used = _sum_hist(current_path, prefix)
        if not np.allclose(edges, current_edges):
            raise ValueError(f"old/current binning mismatch for {title}")

        old_metrics = _metrics(old_values, edges)
        current_metrics = _metrics(current_values, edges)
        centers = 0.5 * (edges[:-1] + edges[1:])

        ax.step(edges[:-1], _normalize(old_values), where="post", color="#222222", linewidth=2.0, label="older working MC")
        ax.step(edges[:-1], _normalize(current_values), where="post", color="#c62828", linewidth=2.2, label="current table-QA MC")
        ax.set_xlim(0, 1)
        y_max = max(np.max(_normalize(old_values)), np.max(_normalize(current_values))) * 1.12
        ax.set_ylim(0, min(1.05, max(0.08, y_max)))
        ax.grid(True, axis="y", color="#d9d9d9", linewidth=0.8, alpha=0.8)
        ax.set_title(title, loc="left", fontsize=13, fontweight="bold", pad=8)
        ax.set_xlabel(r"$E_{11}/E_{33}$")
        ax.set_ylabel("normalized counts")
        ax.text(
            0.04,
            0.94,
            "same legacy ROOT keys\n"
            f"old x<0.05: {old_metrics['low_edge_fraction_x_lt_0p05']:.3f}\n"
            f"current x<0.05: {current_metrics['low_edge_fraction_x_lt_0p05']:.3f}",
            transform=ax.transAxes,
            va="top",
            ha="left",
            fontsize=10.5,
            bbox={"boxstyle": "round,pad=0.28", "facecolor": "white", "edgecolor": "#c9c9c9", "alpha": 0.95},
        )
        ax.axvspan(0.0, 0.05, color="#c62828", alpha=0.08, linewidth=0)
        ax.legend(loc="upper right", frameon=False, fontsize=10.5)

        manifest["rows"].append(
            {
                "title": title,
                "histogram_prefix": prefix,
                "old_path": str(old_path),
                "current_path": str(current_path),
                "old_histograms": old_used,
                "current_histograms": current_used,
                "old_metrics": old_metrics,
                "current_metrics": current_metrics,
                "bin_centers": centers.tolist(),
            }
        )

    fig.suptitle(
        "Embedded AuAu MC E11/E33 reconstruction-contract diagnostic",
        fontsize=19,
        fontweight="bold",
        y=0.985,
    )
    fig.text(
        0.5,
        0.935,
        "Current output is not a plotting/binning artifact: exact shared histograms move from near-zero low-edge content to dominant first-bin content.",
        ha="center",
        va="center",
        fontsize=13,
    )
    fig.subplots_adjust(left=0.065, right=0.985, top=0.875, bottom=0.08, wspace=0.18, hspace=0.33)
    fig.savefig(OUT_PNG, dpi=220)
    plt.close(fig)

    OUT_JSON.write_text(json.dumps(manifest, indent=2) + "\n")
    print(OUT_PNG)
    print(OUT_JSON)


if __name__ == "__main__":
    main()
