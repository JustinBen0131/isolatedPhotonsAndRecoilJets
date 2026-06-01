#!/usr/bin/env python3
from __future__ import annotations
# Keep purpose-folder helpers runnable when invoked directly.
import sys as _codex_sys
from pathlib import Path as _CodexPath
_CODEX_THIS_FILE = _CodexPath(__file__).resolve()
_CODEX_SCRIPTS_DIR = next((p for p in _CODEX_THIS_FILE.parents if p.name == "scripts"), _CODEX_THIS_FILE.parent)
_CODEX_SCRIPTS_DIR_STR = str(_CODEX_SCRIPTS_DIR)
if _CODEX_SCRIPTS_DIR_STR not in _codex_sys.path:
    _codex_sys.path.append(_CODEX_SCRIPTS_DIR_STR)
del _CODEX_THIS_FILE, _CODEX_SCRIPTS_DIR, _CODEX_SCRIPTS_DIR_STR

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch
from textwrap import fill


ROOT = Path(
    "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/auauMLDiagnosticRuns/"
    "global_etcent_inclusive3_sixpack_20260516_135439"
)
OUTDIR = ROOT / "slideReady" / "global32_feature_diagnostics"
SPLIT_CSV = OUTDIR / "global32_bdt_split_gain.csv"


BASE_V3E_PLUS_CENT = {
    "cluster_Et",
    "cluster_Eta",
    "vertexz",
    "centrality",
    "cluster_weta_cogx",
    "cluster_wphi_cogx",
    "cluster_et1",
    "cluster_et2",
    "cluster_et3",
    "cluster_et4",
    "e11_over_e33",
    "e32_over_e35",
}

GROUPS = [
    {
        "group": "base v3E + cent",
        "role": "shared reference",
        "features": sorted(BASE_V3E_PLUS_CENT),
        "color": "#64748B",
        "note": (
            "Grouped foundation: cluster kinematics, vertex, centrality, base widths, "
            "leading local tower energy sharing, E11/E33, E32/E35."
        ),
    },
    {
        "group": "E22 ratios",
        "role": "largest added handle",
        "features": ["e22_over_e33", "e22_over_e35", "e22_over_e37", "e22_over_e53"],
        "color": "#059669",
        "note": "E22/E33, E22/E35, E22/E37, E22/E53",
    },
    {
        "group": "width windows",
        "role": "wider shape information",
        "features": [
            "cluster_weta33_cogx",
            "cluster_wphi33_cogx",
            "cluster_weta35_cogx",
            "cluster_wphi53_cogx",
            "cluster_w32",
            "cluster_w52",
            "cluster_w72",
        ],
        "color": "#0284C7",
        "note": "w_eta33, w_phi33, w_eta35, w_phi53, w32, w52, w72",
    },
    {
        "group": "E11 ratios",
        "role": "additional compactness checks",
        "features": [
            "e11_over_e22",
            "e11_over_e13",
            "e11_over_e15",
            "e11_over_e17",
            "e11_over_e31",
            "e11_over_e51",
            "e11_over_e71",
        ],
        "color": "#F97316",
        "note": "E11/E22, E11/E13, E11/E15, E11/E17, E11/E31, E11/E51, E11/E71",
    },
    {
        "group": "width ratios",
        "role": "eta/phi asymmetry",
        "features": ["cluster_weta_over_wphi", "cluster_weta33_over_wphi33"],
        "color": "#7C3AED",
        "note": "w_eta/w_phi, w_eta33/w_phi33",
    },
]


def grouped_table(split: pd.DataFrame) -> pd.DataFrame:
    rows = []
    seen = set()
    for item in GROUPS:
        features = item["features"]
        seen.update(features)
        sub = split[split["feature"].isin(features)].copy()
        rows.append(
            {
                "group": item["group"],
                "role": item["role"],
                "features": "; ".join(features),
                "display_features": item["note"],
                "split_gain_fraction": float(sub["split_gain_fraction"].sum()),
                "split_gain_percent": 100.0 * float(sub["split_gain_fraction"].sum()),
                "split_count": int(sub["split_count"].sum()),
                "color": item["color"],
            }
        )
    missing = sorted(set(split["feature"]) - seen)
    if missing:
        raise SystemExit(f"Unhandled features in grouped split-gain plot: {missing}")
    return pd.DataFrame(rows)


DISPLAY = {
    "e22_over_e37": "E22/E37",
    "e22_over_e53": "E22/E53",
    "e22_over_e35": "E22/E35",
    "e22_over_e33": "E22/E33",
    "cluster_wphi53_cogx": "w_phi53",
    "cluster_weta35_cogx": "w_eta35",
    "cluster_wphi33_cogx": "w_phi33",
    "cluster_weta33_cogx": "w_eta33",
    "cluster_w72": "w72",
    "cluster_w52": "w52",
    "cluster_w32": "w32",
}


def draw(grouped: pd.DataFrame, split: pd.DataFrame, out: Path) -> None:
    plot_df = grouped.sort_values("split_gain_percent", ascending=False).reset_index(drop=True)
    fig = plt.figure(figsize=(16, 9), dpi=180)

    ax = fig.add_axes([0.125, 0.455, 0.48, 0.335])
    y = np.arange(len(plot_df))
    bars = ax.barh(
        y,
        plot_df["split_gain_percent"],
        height=0.62,
        color=plot_df["color"].tolist(),
        edgecolor="none",
        alpha=0.96,
    )
    ax.set_yticks(y)
    ax.set_yticklabels(plot_df["group"], fontsize=14.5)
    ax.invert_yaxis()
    ax.set_xlabel("Fraction of total split gain [%]", fontsize=14.2, labelpad=5)
    ax.set_xlim(0, 65)
    ax.grid(axis="x", color="#D1D5DB", lw=1.0, alpha=0.7)
    ax.set_axisbelow(True)
    ax.tick_params(axis="x", labelsize=13, direction="out", top=False, length=4)
    ax.tick_params(axis="y", labelsize=14.5, length=0, pad=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(1.1)
    ax.spines["bottom"].set_linewidth(1.1)

    for bar, (_, row) in zip(bars, plot_df.iterrows()):
        value = float(row["split_gain_percent"])
        if value > 38:
            x = value - 1.2
            ha = "right"
            color = "white"
        else:
            x = value + 1.1
            ha = "left"
            color = "#111827"
        ax.text(
            x,
            bar.get_y() + bar.get_height() / 2,
            f"{value:.1f}%  ({int(row['split_count'])})",
            ha=ha,
            va="center",
            fontsize=13.0,
            fontweight="bold",
            color=color,
        )

    text_x = 0.65
    fig.text(text_x, 0.795, "Interpretation", fontsize=17.5, fontweight="bold", color="#111827")
    fig.text(
        text_x,
        0.752,
        "Feature families replace the unreadable 32-row list.\n"
        "The bars show where the BDT actually spends split gain.",
        fontsize=11.0,
        color="#374151",
        linespacing=1.23,
    )

    grouped_idx = grouped.set_index("group")

    def card(y0: float, title: str, pct: float, body: str, color: str, height: float = 0.105) -> None:
        card_x = text_x - 0.012
        card_w = 0.31
        card_center = card_x + 0.5 * card_w
        patch = FancyBboxPatch(
            (card_x, y0 - height + 0.008),
            card_w,
            height,
            boxstyle="round,pad=0.005,rounding_size=0.012",
            linewidth=0,
            facecolor=color,
            alpha=0.10,
            transform=fig.transFigure,
            zorder=0,
        )
        fig.patches.append(patch)
        fig.text(card_center, y0 - 0.003, title, ha="center", fontsize=11.0, fontweight="bold", color=color)
        fig.text(
            card_center,
            y0 - 0.029,
            f"{pct:.1f}% of split gain",
            ha="center",
            fontsize=10.8,
            fontweight="bold",
            color="#111827",
        )
        fig.text(card_center, y0 - 0.051, body, ha="center", fontsize=8.8, color="#111827")

    card(
        0.692,
        "Main added handle: E22 geometry",
        grouped_idx.loc["E22 ratios", "split_gain_percent"],
        "E22/E37 = 44.3%; E22/E53 = 10.7%",
        "#059669",
        height=0.070,
    )
    card(
        0.595,
        "Secondary handle: wider windows",
        grouped_idx.loc["width windows", "split_gain_percent"],
        "3x3, 3x5, 5x3, w32, w52, w72",
        "#0284C7",
        height=0.070,
    )
    card(
        0.498,
        "Reference block remains important",
        grouped_idx.loc["base v3E + cent", "split_gain_percent"],
        "Compact PPG12-like base + centrality",
        "#64748B",
        height=0.070,
    )
    card(
        0.401,
        "Low-leverage additions",
        grouped_idx.loc["E11 ratios", "split_gain_percent"] + grouped_idx.loc["width ratios", "split_gain_percent"],
        "E11 ratios + width ratios are not drivers",
        "#F97316",
        height=0.070,
    )

    fig.text(
        0.06,
        0.330,
        "Exact 32 BDT inputs, grouped exactly as in the split-gain bars",
        ha="left",
        fontsize=13.8,
        fontweight="bold",
        color="#0F172A",
    )

    feature_cards = [
        (
            "base v3E + cent",
            "#64748B",
            "cluster_Et,  cluster_Eta,  vertexz,  centrality\n"
            "cluster_weta_cogx,  cluster_wphi_cogx\n"
            "cluster_et1,  cluster_et2\n"
            "cluster_et3,  cluster_et4\n"
            "e11_over_e33,  e32_over_e35",
            0.06,
            0.176,
            0.350,
            0.128,
        ),
        (
            "E22 ratios",
            "#059669",
            "e22_over_e33\n"
            "e22_over_e35\n"
            "e22_over_e37\n"
            "e22_over_e53",
            0.425,
            0.176,
            0.150,
            0.128,
        ),
        (
            "width windows",
            "#0284C7",
            "cluster_weta33_cogx,  cluster_wphi33_cogx\n"
            "cluster_weta35_cogx,  cluster_wphi53_cogx\n"
            "cluster_w32,  cluster_w52,  cluster_w72",
            0.590,
            0.176,
            0.350,
            0.128,
        ),
        (
            "E11 ratios",
            "#F97316",
            "e11_over_e22,  e11_over_e13,  e11_over_e15\n"
            "e11_over_e17,  e11_over_e31,  e11_over_e51\n"
            "e11_over_e71",
            0.060,
            0.045,
            0.435,
            0.115,
        ),
        (
            "width ratios",
            "#7C3AED",
            "cluster_weta_over_wphi\n"
            "cluster_weta33_over_wphi33",
            0.520,
            0.045,
            0.320,
            0.115,
        ),
    ]

    for title, color, body, x0, y0, w, h in feature_cards:
        patch = FancyBboxPatch(
            (x0, y0),
            w,
            h,
            boxstyle="round,pad=0.008,rounding_size=0.010",
            linewidth=0,
            facecolor=color,
            alpha=0.10,
            transform=fig.transFigure,
            zorder=0,
        )
        fig.patches.append(patch)
        cx = x0 + 0.5 * w
        fig.text(cx, y0 + h - 0.030, title, ha="center", fontsize=13.2, fontweight="bold", color=color)
        fig.text(
            cx,
            y0 + h - 0.052,
            body,
            ha="center",
            va="top",
            fontsize=10.0,
            color="#111827",
            linespacing=1.06,
        )

    fig.text(
        0.06,
        0.925,
        "32-feature BDT: where the added inputs matter",
        ha="left",
        fontsize=24,
        fontweight="bold",
        color="#0F172A",
    )
    fig.text(
        0.06,
        0.875,
        r"Embedded Photon12+20 signal vs Jet12+20+30 background; 15 < cluster $E_T$ < 35 GeV",
        ha="left",
        fontsize=13.2,
        color="#4B5563",
    )
    fig.text(0.845, 0.92, "sPHENIX", ha="right", fontsize=19, fontstyle="italic", fontweight="bold")
    fig.text(0.95, 0.92, "Internal", ha="right", fontsize=19)

    fig.savefig(out, facecolor="white")
    plt.close(fig)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    split = pd.read_csv(SPLIT_CSV)
    grouped = grouped_table(split)
    grouped.to_csv(OUTDIR / "global32_bdt_grouped_split_gain_slide32.csv", index=False)
    draw(grouped, split, OUTDIR / "global32_bdt_grouped_split_gain_slide32.png")
    print(OUTDIR / "global32_bdt_grouped_split_gain_slide32.png")
    print(OUTDIR / "global32_bdt_grouped_split_gain_slide32.csv")


if __name__ == "__main__":
    main()
