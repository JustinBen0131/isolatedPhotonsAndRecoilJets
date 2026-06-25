#!/usr/bin/env python3
"""Build a B-leakage focus slide for THE-69 AuAu photon-ID QA."""

from __future__ import annotations

import csv
import json
import math
import sys
import textwrap
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch


THIS_FILE = Path(__file__).resolve()
REPO = next((p for p in THIS_FILE.parents if (p / "AGENTS.md").exists()), THIS_FILE.parents[4])
SCRIPTS = REPO / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.append(str(SCRIPTS))

from slides.common.slide_defaults import SLIDE_DPI, slide_figsize  # noqa: E402


INPUT_ROOT = (
    REPO
    / "InputFiles/the69_default_auau_physicsqa/RecoilJets_embeddedPhoton12plus20_MERGED.root"
)
OUT_DIR = (
    REPO
    / "dataOutput/auauPhysicsQA/THE69_defaultAuAuBDT_physicsQA_20260620/b_leakage_focus_slide"
)
OUT_PNG = OUT_DIR / "the69_b_leakage_focus_slide.png"
OUT_CSV = OUT_DIR / "the69_b_leakage_focus_points.csv"
OUT_JSON = OUT_DIR / "the69_b_leakage_focus_manifest.json"
OUT_SCRIPT = OUT_DIR / "the69_b_leakage_focus_speaker_script.md"

TOPDIR = "SIM"
ISO_TAG = "isoR40_isSliding"
PT_BINS = [(14, 16), (16, 18), (18, 20), (20, 22), (22, 24), (24, 26), (26, 35)]
CENTRALITIES = [
    ("0_20", "0-20%", "#111827", "o", -0.08),
    ("20_50", "20-50%", "#2563EB", "s", 0.00),
    ("50_80", "50-80%", "#059669", "D", 0.08),
]


@dataclass(frozen=True)
class Count:
    value: float
    error: float


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "axes.facecolor": "white",
            "axes.edgecolor": "#1F2937",
            "axes.linewidth": 1.0,
            "axes.labelsize": 18,
            "axes.titlesize": 21,
            "xtick.labelsize": 15,
            "ytick.labelsize": 15,
            "legend.fontsize": 16,
            "mathtext.fontset": "dejavuserif",
        }
    )


def open_root(path: Path):
    import ROOT

    ROOT.gROOT.SetBatch(True)
    ROOT.TH1.AddDirectory(False)
    root_file = ROOT.TFile.Open(str(path), "READ")
    if not root_file or root_file.IsZombie():
        raise OSError(f"Could not open ROOT file: {path}")
    return root_file


def hist_path(lo: int, hi: int, cent_key: str) -> str:
    return f"{TOPDIR}/h_sigABCD_MC_{ISO_TAG}_pT_{lo}_{hi}_cent_{cent_key}"


def bin_count(hist, bin_index: int) -> Count:
    value = float(hist.GetBinContent(bin_index))
    error = float(hist.GetBinError(bin_index))
    if error <= 0.0 and value > 0.0:
        error = math.sqrt(value)
    return Count(value, error)


def ratio_with_error(num: Count, den: Count) -> tuple[float, float]:
    if den.value <= 0.0 or not np.isfinite(den.value):
        return float("nan"), float("nan")
    value = num.value / den.value
    if num.value > 0.0:
        rel2 = (num.error / num.value) ** 2 + (den.error / den.value) ** 2
        error = abs(value) * math.sqrt(max(0.0, rel2))
    else:
        error = num.error / den.value
    return value, error


def collect_points() -> list[dict]:
    root_file = open_root(INPUT_ROOT)
    rows: list[dict] = []
    try:
        for cent_key, cent_label, _, _, _ in CENTRALITIES:
            for lo, hi in PT_BINS:
                path = hist_path(lo, hi, cent_key)
                hist = root_file.Get(path)
                if not hist:
                    raise KeyError(f"Missing histogram: {path}")
                a_count = bin_count(hist, 1)
                b_count = bin_count(hist, 2)
                b_over_a, b_over_a_err = ratio_with_error(b_count, a_count)
                b_over_tight = b_count.value / (a_count.value + b_count.value) if (a_count.value + b_count.value) > 0 else float("nan")
                rows.append(
                    {
                        "centrality": cent_label,
                        "cent_key": cent_key,
                        "pt_lo": lo,
                        "pt_hi": hi,
                        "pt_mid": 0.5 * (lo + hi),
                        "a_signal": a_count.value,
                        "a_signal_err": a_count.error,
                        "b_signal": b_count.value,
                        "b_signal_err": b_count.error,
                        "b_over_a": b_over_a,
                        "b_over_a_err": b_over_a_err,
                        "b_over_a_plus_b": b_over_tight,
                        "histogram": path,
                    }
                )
    finally:
        root_file.Close()
    return rows


def write_points(rows: list[dict]) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "centrality",
        "cent_key",
        "pt_lo",
        "pt_hi",
        "pt_mid",
        "a_signal",
        "a_signal_err",
        "b_signal",
        "b_signal_err",
        "b_over_a",
        "b_over_a_err",
        "b_over_a_plus_b",
        "histogram",
    ]
    with OUT_CSV.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def add_card(fig, x: float, y: float, w: float, h: float, title: str, lines: list[str], accent: str) -> None:
    ax = fig.add_axes([x, y, w, h])
    ax.set_axis_off()
    box = FancyBboxPatch(
        (0.0, 0.0),
        1.0,
        1.0,
        boxstyle="round,pad=0.014,rounding_size=0.018",
        facecolor="white",
        edgecolor="#CBD5E1",
        linewidth=1.3,
        transform=ax.transAxes,
    )
    ax.add_patch(box)
    ax.add_patch(
        FancyBboxPatch(
            (0.0, 0.0),
            0.015,
            1.0,
            boxstyle="round,pad=0.014,rounding_size=0.018",
            facecolor=accent,
            edgecolor=accent,
            linewidth=0,
            transform=ax.transAxes,
        )
    )
    ax.text(0.05, 0.84, title, fontsize=18.5, fontweight="bold", color="#111827", va="center")
    y_text = 0.66
    for line in lines:
        ax.text(0.05, y_text, line, fontsize=13.8, color="#334155", va="top")
        y_text -= 0.135


def make_slide(rows: list[dict]) -> None:
    setup_style()
    fig = plt.figure(figsize=slide_figsize(), dpi=SLIDE_DPI)
    fig.text(
        0.045,
        0.94,
        "B leakage is the isolation-sideband signal component",
        fontsize=27,
        fontweight="bold",
        color="#111827",
        va="top",
    )
    fig.text(
        0.045,
        0.895,
        "Truth-isolated prompt signal, reconstructed as tight but non-isolated; pre-fix THE-69 output shown as a diagnostic.",
        fontsize=14.8,
        color="#334155",
        va="top",
    )

    ax = fig.add_axes([0.075, 0.205, 0.585, 0.575])
    for cent_key, cent_label, color, marker, offset in CENTRALITIES:
        subset = [row for row in rows if row["cent_key"] == cent_key]
        x = np.array([row["pt_mid"] + offset for row in subset], dtype=float)
        y = np.array([row["b_over_a"] for row in subset], dtype=float)
        yerr = np.array([row["b_over_a_err"] for row in subset], dtype=float)
        ax.errorbar(
            x,
            y,
            yerr=yerr,
            fmt=marker,
            linestyle="none",
            markersize=9.2,
            capsize=3.2,
            elinewidth=1.35,
            markerfacecolor=color,
            markeredgecolor="white",
            markeredgewidth=0.85,
            color=color,
            label=cent_label,
        )

    ax.set_title(r"B leakage: $f_B = B_{\rm sig}/A_{\rm sig}$", pad=10, fontweight="bold")
    ax.set_xlabel(r"reco cluster $E_T$ bin center [GeV]")
    ax.set_ylabel(r"truth-signal leakage into B / A")
    ax.set_xlim(13.5, 35.8)
    ax.set_ylim(0.43, 0.70)
    ax.set_xticks([15, 18, 20, 22, 24, 26, 30, 35])
    ax.grid(True, which="major", color="#DDE3EA", linewidth=0.8, alpha=0.8)
    ax.legend(
        loc="upper right",
        frameon=True,
        facecolor="white",
        edgecolor="#CBD5E1",
        borderpad=0.55,
        handletextpad=0.6,
    )
    ax.text(
        0.03,
        0.06,
        r"$A$: tight + isolated     $B$: tight + non-isolated",
        transform=ax.transAxes,
        fontsize=15,
        color="#334155",
        bbox=dict(boxstyle="round,pad=0.25", facecolor="white", edgecolor="#CBD5E1", alpha=0.94),
    )

    add_card(
        fig,
        0.70,
        0.59,
        0.255,
        0.245,
        "What this is",
        [
            r"$B$ is not fake purity by itself.",
            "It is true signal that passes",
            "tight ID but fails reco isolation.",
            r"Large $B/A$ means the isolation",
            "cut is rejecting real signal.",
        ],
        "#B45309",
    )
    add_card(
        fig,
        0.70,
        0.345,
        0.255,
        0.215,
        "Cut shown",
        [
            r"Reco isolation: $R=0.4$ sliding.",
            r"$E_{\rm iso} < 7.57 - 0.0658c$ GeV.",
            "Non-isolated starts at the",
            "same edge: sideGap = 0.",
        ],
        "#2563EB",
    )
    add_card(
        fig,
        0.70,
        0.12,
        0.255,
        0.205,
        "Rerun caveat",
        [
            "The default config is now fixed",
            "to use centrality WP80.",
            "This plotted output was made",
            "before that fix; rerun tagging.",
        ],
        "#DC2626",
    )

    fig.text(
        0.075,
        0.095,
        "Old slide reference: this is the same kind of B-leakage question, but that slide was an isolation-WP study, not final WP80 purity closure.",
        fontsize=13.8,
        color="#475569",
        va="center",
    )
    fig.savefig(OUT_PNG, dpi=SLIDE_DPI)
    plt.close(fig)


def write_manifest(rows: list[dict]) -> None:
    summary = {}
    for cent_key, cent_label, *_ in CENTRALITIES:
        vals = [row["b_over_a"] for row in rows if row["cent_key"] == cent_key and math.isfinite(row["b_over_a"])]
        summary[cent_label] = {
            "min_b_over_a": min(vals),
            "max_b_over_a": max(vals),
            "mean_b_over_a": sum(vals) / len(vals),
        }
    payload = {
        "input_root": str(INPUT_ROOT),
        "output_png": str(OUT_PNG),
        "output_csv": str(OUT_CSV),
        "histogram_pattern": f"{TOPDIR}/h_sigABCD_MC_{ISO_TAG}_pT_<lo>_<hi>_cent_<cent>",
        "definition": "B leakage = B_sig / A_sig; A=tight+isolated truth signal, B=tight+non-isolated truth signal.",
        "isolation": "isoR40_isSliding, Eiso < 7.57 - 0.0658*c GeV, sideGap=0.",
        "caveat": "The default local config has been fixed to centlinear, but the plotted THE-69 output was produced before that fix; tight/non-tight dependent outputs are diagnostic until rerun.",
        "summary": summary,
    }
    OUT_JSON.write_text(json.dumps(payload, indent=2) + "\n")


def write_speaker_script() -> None:
    OUT_SCRIPT.write_text(
        textwrap.dedent(
            """\
            # Speaker Script: B leakage focus

            This slide isolates only the B-side leakage question.  Here A is the truth-isolated prompt-signal population reconstructed as tight and isolated, while B is the same truth-signal population reconstructed as tight but non-isolated.  The plotted value is B over A.

            The important point is that B leakage is an isolation-sideband diagnostic, not fake purity by itself.  In this current output it is large, about one half to two thirds of A depending on centrality, because the non-isolated side uses the same sliding isolation boundary with no additional sideband gap.

            I would not turn this into a final purity conclusion yet.  The default config is now fixed to use the centrality-linear working-point keyword, but this plotted output was produced before that fix.  The tight/non-tight tagging must be rerun before final ABCD or purity statements.  This slide explains what the B number means and why it is visible.
            """
        )
    )


def main() -> None:
    rows = collect_points()
    write_points(rows)
    make_slide(rows)
    write_manifest(rows)
    write_speaker_script()
    print(OUT_PNG)
    print(OUT_CSV)
    print(OUT_JSON)
    print(OUT_SCRIPT)


if __name__ == "__main__":
    main()
