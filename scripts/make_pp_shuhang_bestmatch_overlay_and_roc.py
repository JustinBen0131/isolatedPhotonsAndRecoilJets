#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import textwrap
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(
    "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/ppPhotonMLPipeline/"
    "ppg12_basev3E_currentIAN_finalShuhang_20260526_2115/validation/fullsim_shuhang_overlay"
)
SUMMARY = ROOT / "pp_currentian_basev3e_bdt_score_overlay_vs_shuhang_clean_summary.json"
VARIANTS = ROOT / "shuhang_projection_variant_scan.csv"


def load_variants() -> dict[tuple[str, str], np.ndarray]:
    rows: dict[tuple[str, str], list[tuple[float, float]]] = defaultdict(list)
    with VARIANTS.open() as handle:
        for row in csv.DictReader(handle):
            key = (row["sample"], f"pt{row['pt']} cut{row['cut']} {row['slice']}")
            rows[key].append((float(row["lo"]), float(row["value"])))
    return {key: np.asarray([v for _, v in sorted(vals)], dtype=float) for key, vals in rows.items()}


def auc_from_hists(signal: np.ndarray, background: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    sig = signal / signal.sum()
    bkg = background / background.sum()
    # Threshold moves from high score to low score.
    tpr = np.r_[0.0, np.cumsum(sig[::-1])]
    fpr = np.r_[0.0, np.cumsum(bkg[::-1])]
    order = np.argsort(fpr)
    auc = float(np.trapz(tpr[order], fpr[order]))
    return fpr, tpr, auc


def rms(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.sqrt(np.mean((a - b) ** 2)))


def main() -> None:
    data = json.loads(SUMMARY.read_text())
    bins = np.asarray(data["bins"], dtype=float)
    centers = 0.5 * (bins[:-1] + bins[1:])
    ours_sig = np.asarray(data["this_analysis_signal_hist"], dtype=float)
    ours_bkg = np.asarray(data["this_analysis_inclusive_hist"], dtype=float)
    variants = load_variants()

    common = sorted({v for _, v in variants if ("signal", v) in variants and ("inclusive", v) in variants})
    common_scores = []
    for variant in common:
        ds = rms(ours_sig, variants[("signal", variant)])
        db = rms(ours_bkg, variants[("inclusive", variant)])
        common_scores.append((ds + db, ds, db, variant))
    common_scores.sort()
    best_common = common_scores[0][3]

    best_signal = min(
        ((rms(ours_sig, hist), variant) for (sample, variant), hist in variants.items() if sample == "signal")
    )[1]
    best_bkg = min(
        ((rms(ours_bkg, hist), variant) for (sample, variant), hist in variants.items() if sample == "inclusive")
    )[1]

    fpr_ours, tpr_ours, auc_ours = auc_from_hists(ours_sig, ours_bkg)
    fpr_common, tpr_common, auc_common = auc_from_hists(
        variants[("signal", best_common)], variants[("inclusive", best_common)]
    )
    fpr_cut1, tpr_cut1, auc_cut1 = auc_from_hists(
        variants[("signal", "pt2 cut1 allY")], variants[("inclusive", "pt2 cut1 allY")]
    )

    fig, axes = plt.subplots(1, 2, figsize=(16, 9), dpi=180)
    ax = axes[0]
    ax.step(centers, ours_sig, where="mid", color="#d62728", linewidth=2.4, label="This analysis signal")
    ax.step(centers, ours_bkg, where="mid", color="#1f5eff", linewidth=2.4, label="This analysis inclusive")
    ax.step(
        centers,
        variants[("signal", best_common)],
        where="mid",
        color="#8b1a1a",
        linestyle="--",
        linewidth=2.0,
        label=f"Shuhang signal ({best_common})",
    )
    ax.step(
        centers,
        variants[("inclusive", best_common)],
        where="mid",
        color="#173a9a",
        linestyle="--",
        linewidth=2.0,
        label=f"Shuhang inclusive ({best_common})",
    )
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 0.56)
    ax.set_xlabel("BDT score")
    ax.set_ylabel("unit-normalized counts")
    ax.set_title("Closest common Shuhang histogram choice")
    ax.grid(True, alpha=0.22)
    ax.legend(loc="upper center", bbox_to_anchor=(0.52, 1.03), ncol=1, fontsize=10, frameon=False)

    ax = axes[1]
    ax.plot(fpr_ours, tpr_ours, color="black", linewidth=2.6, label=f"This analysis AUC={auc_ours:.3f}")
    ax.plot(fpr_common, tpr_common, color="#315fbd", linewidth=2.2, label=f"Shuhang best common AUC={auc_common:.3f}")
    ax.plot(fpr_cut1, tpr_cut1, color="#8b1a1a", linestyle="--", linewidth=2.0, label=f"Shuhang pt2 cut1 AUC={auc_cut1:.3f}")
    ax.plot([0, 1], [0, 1], color="#777777", linewidth=1.2, linestyle=":")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("False positive rate (inclusive accepted)")
    ax.set_ylabel("True positive rate (signal accepted)")
    ax.set_title("ROC from unit-normalized score histograms")
    ax.grid(True, alpha=0.22)
    ax.legend(loc="lower right", fontsize=10, frameon=False)

    fig.suptitle("pp current-IAN baseV3E BDT: Shuhang Overlay Diagnostic", fontsize=21, fontweight="bold")
    note = (
        f"Best common Shuhang choice by RMS(signal+inclusive): {best_common}. "
        f"Best individual matches are signal={best_signal} and inclusive={best_bkg}; those are diagnostic only "
        "because they do not use one common cut definition. The ROC curves here are derived from unit-normalized "
        "BDT-score histograms, not from event-level truth labels in Shuhang's ROOT files."
    )
    fig.text(
        0.055,
        0.055,
        "\n".join(textwrap.wrap(note, width=170)),
        fontsize=11.5,
        color="#333333",
        va="bottom",
    )
    fig.tight_layout(rect=[0.035, 0.13, 0.98, 0.91])
    out = ROOT / "pp_currentian_basev3e_bdt_shuhang_bestmatch_overlay_roc.png"
    fig.savefig(out)

    report = {
        "best_common_variant": best_common,
        "best_common_scores": [
            {"variant": v, "combined_rms": c, "signal_rms": s, "inclusive_rms": b}
            for c, s, b, v in common_scores
        ],
        "best_individual_signal_variant": best_signal,
        "best_individual_inclusive_variant": best_bkg,
        "auc_this_analysis": auc_ours,
        "auc_shuhang_best_common": auc_common,
        "auc_shuhang_pt2_cut1_allY": auc_cut1,
        "plot": str(out),
    }
    (ROOT / "pp_currentian_basev3e_bdt_shuhang_bestmatch_overlay_roc_summary.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n"
    )
    print(out)


if __name__ == "__main__":
    main()
