#!/usr/bin/env python3
"""Build a slide PNG comparing current pp BDT output with AuAu 50-80%."""

from __future__ import annotations

import json
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402


REPO = Path(__file__).resolve().parents[1]

PP_SUMMARY = (
    REPO
    / "dataOutput/ppPhotonMLPipeline/ppg12_basev3E_currentIAN_rawOverlayEnvFix_20260527_1630"
    / "validation/fullsim_shuhang_overlay_raw_inclusive_pt1535"
    / "pp_currentian_basev3e_bdt_score_overlay_vs_ppg12_pp_noCent_bdt_15_35_noNPB_eta0-pt3-cut0_summary.json"
)
AUAU_COMPACT_HIST = (
    REPO
    / "dataOutput/auauTightBDTValidation/THE8_branchA_ladder_scorecache_fullstat_20260527"
    / "the8_branch_a_ladder_compact_score_histograms.json"
)
AUAU_PANEL_SUMMARY = (
    REPO
    / "dataOutput/auauTightBDTValidation/THE8_branchA_ladder_scorecache_fullstat_20260527"
    / "slide23_candidates/the8_branchA_ladder_score_separation_slide23_style.json"
)
AUAU_BRANCH = "Jet12+20+30+40"
AUAU_CENTRALITY_KEY = "50_80"
AUAU_CENTRALITY_LABEL = "50-80%"

OUT_DIR = REPO / "dataOutput/jstg_slide_candidates/pp_vs_auau_50_80_bdt_check"
OUT_PNG = OUT_DIR / "pp_vs_auau_50_80_bdt_score_and_roc_slide.png"
OUT_JSON = OUT_DIR / "pp_vs_auau_50_80_bdt_score_and_roc_slide.json"
OUT_CSV = OUT_DIR / "pp_vs_auau_50_80_bdt_score_histograms.csv"


INK = "#111827"
MUTED = "#4b5563"
GRID = "#d9dee8"
BLUE = "#1f77b4"
ORANGE = "#f97316"
TEAL = "#0891b2"
PURPLE = "#7c3aed"
GREEN_BG = "#eef8f1"
GREEN_EDGE = "#8bd0a2"
BLUE_BG = "#eef5ff"
BLUE_EDGE = "#9dbcf5"
AMBER_BG = "#fff7e8"
AMBER_EDGE = "#e8b45b"


def density_to_prob(edges: np.ndarray, density: np.ndarray) -> np.ndarray:
    widths = np.diff(edges)
    probs = np.asarray(density, dtype=float) * widths
    total = float(np.sum(probs))
    if total <= 0:
        raise ValueError("non-positive histogram probability total")
    return probs / total


def auc_from_probs(sig_prob: np.ndarray, bkg_prob: np.ndarray) -> float:
    bkg_below = np.r_[0.0, np.cumsum(bkg_prob)[:-1]]
    favorable = float(np.sum(sig_prob * bkg_below))
    ties = 0.5 * float(np.sum(sig_prob * bkg_prob))
    return favorable + ties


def roc_from_probs(sig_prob: np.ndarray, bkg_prob: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    # Threshold scan from high score to low score.
    tpr = np.r_[0.0, np.cumsum(sig_prob[::-1])]
    fpr = np.r_[0.0, np.cumsum(bkg_prob[::-1])]
    return fpr, tpr


def step(ax: plt.Axes, edges: np.ndarray, density: np.ndarray, *, color: str, label: str, ls: str = "-") -> None:
    y = np.r_[density, density[-1]]
    ax.step(edges, y, where="post", color=color, lw=3.0, ls=ls, label=label)


def load_pp() -> dict[str, object]:
    with PP_SUMMARY.open() as handle:
        payload = json.load(handle)
    edges = np.asarray(payload["bins"], dtype=float)
    sig = np.asarray(payload["this_analysis_signal_hist"], dtype=float)
    bkg = np.asarray(payload["this_analysis_inclusive_hist"], dtype=float)
    sig_prob = density_to_prob(edges, sig)
    bkg_prob = density_to_prob(edges, bkg)
    fpr, tpr = roc_from_probs(sig_prob, bkg_prob)
    return {
        "label": "pp current-IAN",
        "score_note": "15 < E_T < 35 GeV, |eta| < 0.7, no NPB cut",
        "edges": edges,
        "signal_density": sig,
        "background_density": bkg,
        "signal_entries": int(payload["signal"]["rows_after_cuts"]),
        "background_entries": int(payload["inclusive"]["rows_after_cuts"]),
        "auc": auc_from_probs(sig_prob, bkg_prob),
        "fpr": fpr,
        "tpr": tpr,
        "source": str(PP_SUMMARY),
    }


def load_auau() -> dict[str, object]:
    compact = json.loads(AUAU_COMPACT_HIST.read_text())
    summary = json.loads(AUAU_PANEL_SUMMARY.read_text())
    branch = next((item for item in compact["branches"] if item["label"] == AUAU_BRANCH), None)
    if branch is None:
        raise SystemExit(f"No AuAu branch found for {AUAU_BRANCH}")
    cent = branch["by_centrality"].get(AUAU_CENTRALITY_KEY)
    if cent is None:
        raise SystemExit(f"No AuAu centrality block found for {AUAU_CENTRALITY_KEY}")
    row = next(
        (
            item
            for item in summary["rows"]
            if item["branch"] == AUAU_BRANCH and item["centrality"] == AUAU_CENTRALITY_LABEL
        ),
        None,
    )
    if row is None:
        raise SystemExit(f"No AuAu summary row found for {AUAU_BRANCH}, {AUAU_CENTRALITY_LABEL}")

    edges = np.asarray(branch["bin_edges"], dtype=float)
    sig = np.asarray(cent["signal"]["density"], dtype=float)
    bkg = np.asarray(cent["background"]["density"], dtype=float)
    sig_prob = density_to_prob(edges, sig)
    bkg_prob = density_to_prob(edges, bkg)
    fpr, tpr = roc_from_probs(sig_prob, bkg_prob)
    return {
        "label": f"AuAu {AUAU_CENTRALITY_LABEL}",
        "score_note": f"Branch A {AUAU_BRANCH}, baseV3E + centrality, 15 < E_T < 35 GeV",
        "edges": edges,
        "signal_density": sig,
        "background_density": bkg,
        "signal_entries": int(cent["signal"]["entries"]),
        "background_entries": int(cent["background"]["entries"]),
        "auc": float(row["auc_binned"]),
        "auc_binned": auc_from_probs(sig_prob, bkg_prob),
        "inclusive_auc": float(row["inclusive_auc"]),
        "fpr": fpr,
        "tpr": tpr,
        "source": str(AUAU_COMPACT_HIST),
        "summary_source": str(AUAU_PANEL_SUMMARY),
        "branch": AUAU_BRANCH,
        "centrality": AUAU_CENTRALITY_LABEL,
    }


def write_hist_csv(pp: dict[str, object], auau: dict[str, object]) -> None:
    rows: list[dict[str, object]] = []
    for block in (pp, auau):
        edges = np.asarray(block["edges"], dtype=float)
        for cls, vals in (
            ("signal", np.asarray(block["signal_density"], dtype=float)),
            ("background", np.asarray(block["background_density"], dtype=float)),
        ):
            for lo, hi, den in zip(edges[:-1], edges[1:], vals, strict=True):
                rows.append(
                    {
                        "sample": block["label"],
                        "class": cls,
                        "bin_low": lo,
                        "bin_high": hi,
                        "density": den,
                        "auc": block["auc"],
                        "signal_entries": block["signal_entries"],
                        "background_entries": block["background_entries"],
                    }
                )
    pd.DataFrame(rows).to_csv(OUT_CSV, index=False)


def make_slide(pp: dict[str, object], auau: dict[str, object]) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "font.family": ["Times New Roman", "Times", "DejaVu Serif"],
            "axes.labelsize": 20,
            "axes.titlesize": 23,
            "xtick.labelsize": 15,
            "ytick.labelsize": 15,
            "legend.fontsize": 14.5,
        }
    )

    fig = plt.figure(figsize=(16, 9), dpi=160)
    fig.patch.set_facecolor("#f8fafc")

    ax_score = fig.add_axes([0.055, 0.345, 0.520, 0.495])
    ax_roc = fig.add_axes([0.625, 0.505, 0.325, 0.335])

    fig.text(0.055, 0.942, "pp vs Peripheral AuAu BDT Comparison", fontsize=32, weight="bold", color=INK)
    fig.text(
        0.055,
        0.905,
        "Area-normalized BDT score shapes and ROC curves over the common 15-35 GeV photon-candidate window.",
        fontsize=17.5,
        color=MUTED,
    )

    step(ax_score, pp["edges"], pp["signal_density"], color=BLUE, label="pp signal")
    step(ax_score, pp["edges"], pp["background_density"], color=ORANGE, label="pp inclusive")
    step(ax_score, auau["edges"], auau["signal_density"], color=TEAL, label="AuAu 50-80 signal", ls="--")
    step(ax_score, auau["edges"], auau["background_density"], color=PURPLE, label="AuAu 50-80 background", ls="--")
    ax_score.set_xlim(0.0, 1.0)
    ymax = max(
        np.nanmax(pp["signal_density"]),
        np.nanmax(pp["background_density"]),
        np.nanmax(auau["signal_density"]),
        np.nanmax(auau["background_density"]),
    )
    positive = []
    for arr in (
        pp["signal_density"],
        pp["background_density"],
        auau["signal_density"],
        auau["background_density"],
    ):
        arr = np.asarray(arr, dtype=float)
        positive.extend(arr[arr > 0].tolist())
    ymin = max(min(positive) * 0.55, 2.0e-3)
    ax_score.set_yscale("log")
    ax_score.set_ylim(ymin, ymax * 2.65)
    ax_score.set_title("Signal/background BDT score shapes", loc="left", weight="bold", color=INK, pad=10)
    ax_score.set_xlabel("BDT score")
    ax_score.set_ylabel("Area-normalized density (log)")
    ax_score.grid(True, color=GRID, lw=0.8, alpha=0.72, which="both")
    ax_score.tick_params(direction="in", top=True, right=True, length=6)
    ax_score.legend(
        loc="upper right",
        bbox_to_anchor=(0.985, 0.985),
        frameon=True,
        facecolor="white",
        edgecolor="#d1d5db",
        framealpha=0.92,
        ncol=2,
        borderpad=0.55,
        columnspacing=1.15,
        handlelength=2.4,
    )

    ax_roc.plot(pp["fpr"], pp["tpr"], color=BLUE, lw=3.0, label=f"pp AUC {pp['auc']:.3f}")
    ax_roc.plot(auau["fpr"], auau["tpr"], color=PURPLE, lw=3.0, ls="--", label=f"AuAu 50-80 AUC {auau['auc']:.3f}")
    ax_roc.plot([0, 1], [0, 1], color="#9ca3af", lw=1.6, ls=":")
    ax_roc.set_xlim(0, 1)
    ax_roc.set_ylim(0, 1)
    ax_roc.set_xlabel("Background efficiency")
    ax_roc.set_ylabel("Signal efficiency")
    ax_roc.set_title("ROC comparison", loc="left", weight="bold", color=INK, pad=12)
    ax_roc.grid(True, color=GRID, lw=0.8, alpha=0.75)
    ax_roc.tick_params(direction="in", top=True, right=True, length=6)
    ax_roc.legend(loc="lower right", frameon=True, facecolor="white", edgecolor="#d1d5db")

    def box(
        x: float,
        y: float,
        w: float,
        h: float,
        title: str,
        body: str,
        fc: str,
        ec: str,
        *,
        wrap: int,
        body_size: float = 13.0,
    ) -> None:
        patch = matplotlib.patches.FancyBboxPatch(
            (x, y),
            w,
            h,
            boxstyle="round,pad=0.012,rounding_size=0.008",
            transform=fig.transFigure,
            facecolor=fc,
            edgecolor=ec,
            linewidth=1.2,
        )
        fig.add_artist(patch)
        fig.text(x + 0.014, y + h - 0.030, title, fontsize=18, weight="bold", color=INK, va="top")
        wrapped = "\n".join(textwrap.fill(part, width=wrap) for part in body.split("\n"))
        fig.text(
            x + 0.014,
            y + h - 0.066,
            wrapped,
            fontsize=body_size,
            color="#253247",
            va="top",
            linespacing=1.25,
        )

    box(
        0.635,
        0.290,
        0.315,
        0.135,
        "Inputs",
        "pp: 15-35 GeV, no NPB cut\nAuAu: Branch A Jet12+20+30+40, 50-80%, 15-35 GeV",
        BLUE_BG,
        BLUE_EDGE,
        wrap=40,
        body_size=13.2,
    )
    box(
        0.055,
        0.075,
        0.445,
        0.135,
        "Interpretation",
        "This checks the pp and peripheral AuAu photon-ID ranking shapes over the same 15-35 GeV window. The ROC values show both classifiers remain strongly separating.",
        GREEN_BG,
        GREEN_EDGE,
        wrap=68,
    )
    box(
        0.525,
        0.075,
        0.425,
        0.135,
        "Caveat",
        "This is now an apples-to-apples ET-window comparison. The remaining differences are sample environment and training contract: pp raw-inclusive full simulation versus AuAu embedded Branch A validation.",
        AMBER_BG,
        AMBER_EDGE,
        wrap=66,
    )

    for ax in (ax_score, ax_roc):
        for spine in ax.spines.values():
            spine.set_color(INK)
            spine.set_linewidth(1.2)

    fig.savefig(OUT_PNG)
    plt.close(fig)

    OUT_JSON.write_text(
        json.dumps(
            {
                "output_png": str(OUT_PNG),
                "canvas_px": [2560, 1440],
                "pp": {
                    "source": pp["source"],
                    "signal_entries": pp["signal_entries"],
                    "background_entries": pp["background_entries"],
                    "auc_binned_from_display_hist": pp["auc"],
                    "score_note": pp["score_note"],
                },
                "auau_50_80": {
                    "source": auau["source"],
                    "summary_source": auau["summary_source"],
                    "branch": auau["branch"],
                    "centrality": auau["centrality"],
                    "signal_entries": auau["signal_entries"],
                    "background_entries": auau["background_entries"],
                    "auc_from_branch_a_panel_row": auau["auc"],
                    "auc_binned_recomputed": auau["auc_binned"],
                    "inclusive_auc_all_centralities": auau["inclusive_auc"],
                    "score_note": auau["score_note"],
                },
            },
            indent=2,
        )
        + "\n"
    )


def main() -> None:
    pp = load_pp()
    auau = load_auau()
    make_slide(pp, auau)
    write_hist_csv(pp, auau)
    print(OUT_PNG)
    print(OUT_JSON)
    print(OUT_CSV)


if __name__ == "__main__":
    main()
