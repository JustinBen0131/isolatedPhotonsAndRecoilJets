#!/usr/bin/env python3
"""Plot BDT-score symptoms in the retired tainted AuAu table-QA MC output."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import uproot


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
SIGNAL_ROOT = REPO / "InputFiles/the42_current_auau_tableqa/RecoilJets_embeddedPhoton12plus20_MERGED.root"
INCLUSIVE_ROOT = REPO / "InputFiles/the42_current_auau_tableqa/RecoilJets_embeddedJet12plus20plus30plus40_MERGED.root"
OUT_DIR = REPO / "dataOutput/auauTableQA/the42_current_auau_reco_contract_debug"
OUT_PNG = OUT_DIR / "tainted_current_auau_bdt_score_symptom_diagnostic.png"
OUT_JSON = OUT_DIR / "tainted_current_auau_bdt_score_symptom_diagnostic_manifest.json"

CENTRALITIES = (("0-20%", "0_20"), ("20-50%", "20_50"), ("50-80%", "50_80"))
STAGES = (("cut0", "Before preselection"), ("cut1", "After preselection"), ("cut2", "After tight BDT"))
PT_TOKEN = "1535"


def load_hist(path: Path, cent: str, cut: str) -> tuple[np.ndarray, np.ndarray, str]:
    name = f"h1d_bdt_eta0_pt{PT_TOKEN}_cent{cent}_{cut}"
    with uproot.open(path) as f:
        d = f["SIM"]
        if name not in d:
            raise KeyError(f"{path} missing SIM/{name}")
        y, e = d[name].to_numpy(flow=False)
    return np.asarray(y, dtype=float), np.asarray(e, dtype=float), f"SIM/{name}"


def norm(y: np.ndarray) -> np.ndarray:
    total = float(np.sum(y))
    return y / total if total > 0 else y


def metrics(y: np.ndarray, e: np.ndarray) -> dict[str, float]:
    c = 0.5 * (e[:-1] + e[1:])
    total = float(np.sum(y))
    sumw2 = float(np.sum(y * y))
    return {
        "sumw": total,
        "mean": float(np.sum(y * c) / total) if total else 0.0,
        "frac_score_lt_0p05": float(np.sum(y[c < 0.05]) / total) if total else 0.0,
        "frac_score_lt_0p20": float(np.sum(y[c < 0.20]) / total) if total else 0.0,
        "frac_score_gt_0p80": float(np.sum(y[c > 0.80]) / total) if total else 0.0,
        "max_bin_fraction": float(np.max(y) / total) if total else 0.0,
        "effective_entries": float(total * total / sumw2) if sumw2 > 0 else 0.0,
    }


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 11,
        "axes.linewidth": 1.1,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
    })

    fig, axes = plt.subplots(3, 3, figsize=(14.4, 10.0), sharex=True, constrained_layout=False)
    manifest: dict[str, object] = {
        "purpose": "retired tainted AuAu table-QA BDT-score diagnostic",
        "signal_root": str(SIGNAL_ROOT),
        "inclusive_root": str(INCLUSIVE_ROOT),
        "pt_token": PT_TOKEN,
        "caveat": "This plots the retired current embedded-MC output known to have the E11/E33 reconstruction-contract failure; use only to illustrate symptoms.",
        "panels": {},
    }

    for row, (cent_label, cent_token) in enumerate(CENTRALITIES):
        manifest["panels"][cent_label] = {}
        for col, (cut, stage_label) in enumerate(STAGES):
            ax = axes[row, col]
            sig, edges, sig_name = load_hist(SIGNAL_ROOT, cent_token, cut)
            inc, inc_edges, inc_name = load_hist(INCLUSIVE_ROOT, cent_token, cut)
            if not np.allclose(edges, inc_edges):
                raise ValueError(f"bin mismatch {cent_label} {cut}")
            centers = 0.5 * (edges[:-1] + edges[1:])
            sig_n = norm(sig)
            inc_n = norm(inc)
            ax.step(edges[:-1], sig_n, where="post", color="#d62728", linewidth=2.0, label="tainted signal MC")
            ax.step(edges[:-1], inc_n, where="post", color="#1f77b4", linewidth=2.0, label="tainted inclusive MC")
            ax.set_xlim(0, 1)
            ymax = max(float(np.max(sig_n)), float(np.max(inc_n))) * 1.20
            ax.set_ylim(0, max(0.055, ymax))
            ax.grid(True, axis="y", color="#d8dde5", alpha=0.85)
            ax.set_title(f"{cent_label}  {stage_label}", loc="left", fontsize=12.2, fontweight="bold")
            if row == 2:
                ax.set_xlabel("BDT score")
            if col == 0:
                ax.set_ylabel("normalized counts")
            sm = metrics(sig, edges)
            im = metrics(inc, edges)
            ax.text(
                0.035,
                0.955,
                f"sig mean {sm['mean']:.3f}, <0.2 {sm['frac_score_lt_0p20']:.2f}\n"
                f"inc mean {im['mean']:.3f}, <0.2 {im['frac_score_lt_0p20']:.2f}\n"
                f"gap {sm['mean'] - im['mean']:.3f}",
                transform=ax.transAxes,
                va="top",
                ha="left",
                fontsize=8.9,
                bbox={"boxstyle": "round,pad=0.22", "facecolor": "white", "edgecolor": "#c7c7c7", "alpha": 0.94},
            )
            if row == 0 and col == 2:
                ax.legend(loc="upper right", frameon=False, fontsize=9.8)
            manifest["panels"][cent_label][cut] = {
                "signal_histogram": sig_name,
                "inclusive_histogram": inc_name,
                "signal_metrics": sm,
                "inclusive_metrics": im,
                "mean_score_gap_signal_minus_inclusive": sm["mean"] - im["mean"],
            }

    fig.suptitle("Retired AuAu MC BDT-score symptom: weak preselection-stage separation", fontsize=18, fontweight="bold", y=0.985)
    fig.text(
        0.5,
        0.948,
        "These are the current table-QA MC outputs with the known embedded-MC reconstruction-contract failure; cut2 is selection-biased by construction.",
        ha="center",
        va="top",
        fontsize=11.5,
    )
    fig.text(
        0.5,
        0.018,
        "Diagnostic only: the E11/E33 input failure is the primary proof; this shows why the score output should not be used as the AuAu baseline reference.",
        ha="center",
        va="bottom",
        fontsize=10.2,
        color="#4d4d4d",
    )
    fig.subplots_adjust(left=0.07, right=0.985, top=0.89, bottom=0.08, hspace=0.28, wspace=0.15)
    fig.savefig(OUT_PNG, dpi=190)
    plt.close(fig)
    OUT_JSON.write_text(json.dumps(manifest, indent=2) + "\n")
    print(OUT_PNG)
    print(OUT_JSON)


if __name__ == "__main__":
    main()
