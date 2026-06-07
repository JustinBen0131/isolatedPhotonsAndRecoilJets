#!/usr/bin/env python3
"""Build the single cleanest THE-44 local-cone composition master slide."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Polygon, Rectangle


REPO = Path(__file__).resolve().parents[4]
BASE = REPO / "dataOutput/auauTightBDTValidation/THE8_branchA_ladder_scorecache_fullstat_20260527"
THE44_DIR = BASE / "diagnostics/the44_high_bdt_pythia_autopsy_20260605/slideReady/full_autopsy_20260605_2108"
DEFAULT_OUTDIR = THE44_DIR / "particle_intensive/triangle_master_slide_20260606"

INK = "#111827"
MUTED = "#5D6675"
LIGHT = "#F6F7F9"
CARD = "#FFFFFF"
GRID = "#DCE3EC"
RED = "#C43C32"
GREEN = "#2F8F52"
ORANGE = "#D97706"
BLUE = "#2563A6"
PURPLE = "#7651A6"
SLATE = "#475569"
YELLOW = "#FFF4C7"
TEAL = "#0F766E"

SOURCE_COLORS = {"Jet20": ORANGE, "Jet30": BLUE, "Jet40": PURPLE}
SOURCE_MARKERS = {"Jet20": "o", "Jet30": "s", "Jet40": "^"}
GROUP_COLORS = {"photon": RED, "neutral_meson": ORANGE, "charged_hadron": BLUE, "other": SLATE}

plt.rcParams.update(
    {
        "font.family": ["Times New Roman", "Times", "DejaVu Serif"],
        "axes.edgecolor": "#28313F",
        "axes.labelcolor": INK,
        "xtick.color": INK,
        "ytick.color": INK,
    }
)


def canvas() -> tuple[plt.Figure, plt.Axes]:
    fig = plt.figure(figsize=(16, 9), dpi=160)
    fig.patch.set_facecolor(LIGHT)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    return fig, ax


def txt(
    ax: plt.Axes,
    x: float,
    y: float,
    text: str,
    *,
    size: float,
    color: str = INK,
    weight: str = "normal",
    ha: str = "left",
    va: str = "top",
    linespacing: float = 1.08,
) -> None:
    ax.text(
        x,
        y,
        text,
        transform=ax.transAxes,
        fontsize=size,
        color=color,
        fontweight=weight,
        ha=ha,
        va=va,
        linespacing=linespacing,
        zorder=30,
    )


def card(
    ax: plt.Axes,
    x: float,
    y: float,
    w: float,
    h: float,
    *,
    face: str = CARD,
    edge: str = "#D4DCE7",
    lw: float = 1.2,
    radius: float = 0.018,
    zorder: int = 1,
) -> None:
    ax.add_patch(
        FancyBboxPatch(
            (x, y),
            w,
            h,
            transform=ax.transAxes,
            boxstyle=f"round,pad=0.010,rounding_size={radius}",
            facecolor=face,
            edgecolor=edge,
            linewidth=lw,
            zorder=zorder,
        )
    )


def truth_bucket(name: str) -> str:
    if name in {"pi0", "eta"}:
        return "pi0/eta"
    if name == "gamma":
        return "gamma truth"
    return "other"


def load_data() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    cand = pd.read_csv(THE44_DIR / "the44_candidate_autopsy_rows.csv")
    part = pd.read_csv(THE44_DIR / "the44_particle_rows.csv")
    cand["source"] = cand["sample"].str.replace(" inclusive", "", regex=False).str.replace(" signal", "", regex=False)
    cand["truth_bucket"] = cand["truth_pid_name"].map(truth_bucket)
    cand["cone_bucket"] = np.select(
        [cand["final_state_photon_frac"] >= 0.80, cand["final_state_non_photon_frac"] >= 0.80],
        ["photon-rich", "non-photon-rich"],
        default="mixed",
    )
    inc = cand[(cand["label"].eq("high-BDT background")) & cand["sample"].str.contains("Jet", na=False)].copy()
    true = cand[(cand["label"].eq("middling true photon")) & cand["sample"].str.contains("Photon", na=False)].copy()
    return inc, true, part


def ternary_xy(df: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    p = df["final_state_photon_frac"].to_numpy(dtype=float)
    n = df["final_neutral_meson_frac"].to_numpy(dtype=float)
    c = df["final_charged_hadron_frac"].to_numpy(dtype=float)
    total = p + n + c
    total[total <= 0] = 1.0
    p, n, c = p / total, n / total, c / total
    neutral = np.array([0.10, 0.12])
    charged = np.array([0.90, 0.12])
    photon = np.array([0.50, 0.82])
    coords = n[:, None] * neutral + c[:, None] * charged + p[:, None] * photon
    return coords[:, 0], coords[:, 1]


def med_point(df: pd.DataFrame) -> tuple[float, float]:
    vals = pd.DataFrame(
        [
            {
                "final_state_photon_frac": df["final_state_photon_frac"].median(),
                "final_neutral_meson_frac": df["final_neutral_meson_frac"].median(),
                "final_charged_hadron_frac": df["final_charged_hadron_frac"].median(),
            }
        ]
    )
    x, y = ternary_xy(vals)
    return x.item(), y.item()


def draw_triangle(fig: plt.Figure) -> plt.Axes:
    tri = fig.add_axes([0.055, 0.105, 0.660, 0.740])
    tri.set_axis_off()
    tri.set_xlim(0, 1)
    tri.set_ylim(0, 0.92)
    verts = np.array([[0.10, 0.12], [0.90, 0.12], [0.50, 0.82]])
    tri.add_patch(Polygon(verts, closed=True, facecolor=CARD, edgecolor="#CBD5E1", linewidth=2.1, zorder=0))
    for frac in [0.25, 0.50, 0.75]:
        tri.plot([0.10 + 0.40 * frac, 0.90 - 0.40 * frac], [0.12 + 0.70 * frac, 0.12 + 0.70 * frac], color=GRID, linewidth=0.9, zorder=1)
        tri.plot([0.10 + 0.80 * frac, 0.50 + 0.40 * frac], [0.12, 0.82 - 0.70 * frac], color=GRID, linewidth=0.9, zorder=1)
        tri.plot([0.90 - 0.80 * frac, 0.50 - 0.40 * frac], [0.12, 0.82 - 0.70 * frac], color=GRID, linewidth=0.9, zorder=1)
    tri.plot([0.10, 0.90], [0.12, 0.12], color="#EAB308", linewidth=6.0, alpha=0.22, zorder=2)
    tri.text(0.50, 0.860, "photon core", ha="center", va="bottom", fontsize=18, fontweight="bold", color=GREEN)
    tri.text(0.070, 0.075, "neutral meson", ha="left", va="top", fontsize=17, fontweight="bold", color=ORANGE)
    tri.text(0.930, 0.075, "charged hadron", ha="right", va="top", fontsize=17, fontweight="bold", color=BLUE)
    tri.text(0.50, 0.045, "fragmentation edge: status-1 cone is not photon-core dominated", ha="center", va="top", fontsize=14.2, color=SLATE, fontweight="bold")
    return tri


def evidence_pill(ax: plt.Axes, x: float, y: float, w: float, color: str, number: str, label: str, detail: str) -> None:
    card(ax, x, y, w, 0.115, face=CARD, edge=color, lw=1.5, radius=0.018)
    txt(ax, x + 0.018, y + 0.087, number, size=24, color=color, weight="bold")
    txt(ax, x + 0.018, y + 0.052, label, size=14.5, weight="bold")
    txt(ax, x + 0.018, y + 0.023, detail, size=11.8, color=MUTED, linespacing=0.95)


def make_slide(inc: pd.DataFrame, true: pd.DataFrame, part: pd.DataFrame, outdir: Path) -> Path:
    fig, ax = canvas()
    txt(ax, 0.055, 0.940, "The BDT's high-score fake tail is a fragmentation endpoint, not a photon core", size=28.5, weight="bold")
    txt(
        ax,
        0.055,
        0.882,
        "Each point is one candidate mapped to status-1 local-cone pT fractions within Delta R < 0.3.",
        size=17.2,
        color=MUTED,
    )

    tri = draw_triangle(fig)
    x_true, y_true = ternary_xy(true)
    tri.scatter(x_true, y_true, s=18, color=GREEN, alpha=0.16, edgecolors="none", zorder=3)
    for source in ["Jet20", "Jet30", "Jet40"]:
        sub = inc[inc["source"].eq(source)]
        x, y = ternary_xy(sub)
        size = np.clip(42 + 150 * (sub["auau_tight_bdt_score"] - 0.80), 44, 145)
        tri.scatter(
            x,
            y,
            s=size,
            marker=SOURCE_MARKERS[source],
            color=SOURCE_COLORS[source],
            edgecolors="white",
            linewidths=0.85,
            alpha=0.86,
            zorder=5,
        )
    fake_x, fake_y = med_point(inc)
    true_x, true_y = med_point(true)
    tri.scatter([true_x], [true_y], s=340, marker="*", color=GREEN, edgecolors="white", linewidths=1.6, zorder=12)
    tri.scatter([fake_x], [fake_y], s=470, marker="*", color=RED, edgecolors="white", linewidths=1.7, zorder=13)
    tri.text(true_x + 0.040, true_y - 0.012, "truth-photon\nmedian", fontsize=14.2, color=GREEN, fontweight="bold", va="top")
    tri.text(fake_x + 0.022, fake_y + 0.020, "high-BDT fake\nmedian", fontsize=14.5, color=RED, fontweight="bold")

    # Right-side conclusion panel.
    card(ax, 0.735, 0.610, 0.225, 0.275)
    txt(ax, 0.760, 0.832, "Definitive read", size=20, weight="bold")
    txt(ax, 0.760, 0.778, "Truth-photon controls", size=15, color=GREEN, weight="bold")
    txt(ax, 0.760, 0.747, "median photon fraction 100%\nmedian BDT 0.57", size=13.6, color=INK, linespacing=1.0)
    txt(ax, 0.760, 0.675, "High-BDT inclusive backgrounds", size=15, color=RED, weight="bold")
    txt(ax, 0.760, 0.642, "median photon fraction 0%\nmedian non-photon fraction 100%\nmedian BDT 0.84", size=13.2, color=INK, linespacing=0.98)

    card(ax, 0.735, 0.382, 0.225, 0.178, face=YELLOW, edge="#F2CB4C")
    txt(ax, 0.760, 0.522, "One-slide answer", size=18.5, weight="bold")
    txt(
        ax,
        0.760,
        0.478,
        "Dominant high-score path:\nJet30/40 -> pi0/eta ->\nnon-photon cone.\nEndpoint: 71% neutral meson,\n29% charged, 0% photon.",
        size=12.0,
        linespacing=0.94,
    )

    # Evidence chain below the conclusion panel.
    evidence_pill(ax, 0.735, 0.260, 0.103, PURPLE, "146/160", "Jet30/40", "hard parent sample")
    evidence_pill(ax, 0.850, 0.260, 0.110, ORANGE, "126/160", "pi0/eta", "truth seed")
    evidence_pill(ax, 0.735, 0.120, 0.103, BLUE, "138/160", "non-photon", "final local cone")
    evidence_pill(ax, 0.850, 0.120, 0.110, RED, "0.84", "BDT median", "high-score tail")
    ax.add_patch(FancyArrowPatch((0.842, 0.318), (0.850, 0.318), transform=ax.transAxes, arrowstyle="-|>", mutation_scale=14, linewidth=1.6, color=SLATE, zorder=20))
    ax.add_patch(FancyArrowPatch((0.902, 0.257), (0.902, 0.232), transform=ax.transAxes, arrowstyle="-|>", mutation_scale=14, linewidth=1.6, color=SLATE, zorder=20))
    ax.add_patch(FancyArrowPatch((0.842, 0.178), (0.850, 0.178), transform=ax.transAxes, arrowstyle="-|>", mutation_scale=14, linewidth=1.6, color=SLATE, zorder=20))

    for i, source in enumerate(["Jet20", "Jet30", "Jet40"]):
        ax.scatter([0.095 + 0.060 * i], [0.077], transform=ax.transAxes, marker=SOURCE_MARKERS[source], s=75, color=SOURCE_COLORS[source], edgecolors="white", linewidths=0.8, zorder=30)
        txt(ax, 0.112 + 0.060 * i, 0.087, source, size=12.5, color=INK, va="center")

    out = outdir / "the44_triangle_master_slide.png"
    fig.savefig(out, dpi=160)
    plt.close(fig)
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    inc, true, part = load_data()
    slide = make_slide(inc, true, part, args.outdir)
    script_path = args.outdir / "the44_triangle_master_slide_script.md"
    script_path.write_text(
        "The triangle is the whole argument. Truth-photon controls sit at the photon-core corner, while the high-BDT inclusive backgrounds sit on the neutral-meson/charged-hadron fragmentation edge. The dominant path is Jet30/40 to pi0/eta to a non-photon status-1 cone, even though the BDT score remains high.\n"
    )
    manifest = {
        "status": "READY",
        "slide": str(slide),
        "script": str(script_path),
        "metrics": {
            "inclusive_high_bdt_backgrounds": float(len(inc)),
            "truth_photon_controls": float(len(true)),
            "jet30_40_candidates": float(inc["source"].isin(["Jet30", "Jet40"]).sum()),
            "pi0_eta_truth_candidates": float(inc["truth_bucket"].eq("pi0/eta").sum()),
            "non_photon_rich_cones": float(inc["cone_bucket"].eq("non-photon-rich").sum()),
            "median_high_fake_bdt_score": float(inc["auau_tight_bdt_score"].median()),
            "median_high_fake_photon_fraction": float(inc["final_state_photon_frac"].median()),
            "median_truth_photon_bdt_score": float(true["auau_tight_bdt_score"].median()),
            "median_truth_photon_photon_fraction": float(true["final_state_photon_frac"].median()),
        },
        "inputs": {
            "candidate_rows": str(THE44_DIR / "the44_candidate_autopsy_rows.csv"),
            "particle_rows": str(THE44_DIR / "the44_particle_rows.csv"),
        },
        "caveat": "Internal qualitative diagnostic from saved THE-44 candidate and particle rows; not a production-weighted purity estimate.",
    }
    manifest_path = args.outdir / "the44_triangle_master_slide_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
