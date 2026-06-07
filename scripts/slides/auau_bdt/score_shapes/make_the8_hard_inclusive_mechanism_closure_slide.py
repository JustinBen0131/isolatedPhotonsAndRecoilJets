#!/usr/bin/env python3
"""Build the THE-8 Jet30/Jet40 mechanism + WP80 payoff closure slide."""

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

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch, Rectangle

from slides.common.slide_defaults import SLIDE_DPI, slide_figsize


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
AUDIT_DIR = (
    REPO
    / "dataOutput/auauTightBDTValidation/THE8_branchA_ladder_scorecache_fullstat_20260527"
    / "diagnostics/hard_inclusive_jet3040_mechanism_20260605"
)
DEFAULT_OUTDIR = (
    REPO
    / "dataOutput/auauTightBDTValidation/THE8_branchA_ladder_scorecache_fullstat_20260527"
    / "slideReady/hard_inclusive_jet3040_mechanism_closure_20260605"
)

TITLE = "Jet30/40 Gains: Further Diagnostic"
SUBTITLE = (
    r"The largest BDT gain appears exactly where Jet12+20 has the weakest hard-background coverage"
)

INK = "#121826"
MUTED = "#475467"
BORDER = "#C9D2DE"
PANEL = "#FFFFFF"
BG = "#FFFFFF"
SOFT_BLUE = "#F1F7FF"
SOFT_PURPLE = "#F7F2FC"
SOFT_ORANGE = "#FFF8EE"
GRID = "#D9DEE7"
BASELINE_FILL = "#DCEFE2"
EXPANDED_FILL = "#E7DDF7"

BRANCH_COLORS = {
    "Jet12+20": "#3E8E58",
    "Jet12+20+30": "#C99500",
    "Jet12+20+30+40": "#7B55B6",
}
SOURCE_COLORS = {
    "Jet12": "#4C78A8",
    "Jet20": "#F58518",
    "Jet30": "#54A24B",
    "Jet40": "#B279A2",
}
ET_LABELS = ["15-20", "20-25", "25-30", "30-35"]
SAMPLE_LABELS = ["Jet12+20", "Jet12+20+30", "Jet12+20+30+40"]
SOURCE_LABELS = ["Jet12", "Jet20", "Jet30", "Jet40"]


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def as_float(value: object, default: float = math.nan) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def as_int(value: object, default: int = 0) -> int:
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def fmt_count(n: float) -> str:
    if n >= 10000:
        return f"{n / 1000:.1f}k"
    if n >= 1000:
        return f"{n / 1000:.1f}k"
    return str(int(round(n)))


def add_round_box(
    fig: plt.Figure,
    xy: tuple[float, float],
    wh: tuple[float, float],
    face: str,
    *,
    edge: str = BORDER,
    lw: float = 1.2,
    radius: float = 0.012,
    zorder: int = -3,
) -> None:
    fig.patches.append(
        FancyBboxPatch(
            xy,
            wh[0],
            wh[1],
            boxstyle=f"round,pad=0.010,rounding_size={radius}",
            transform=fig.transFigure,
            facecolor=face,
            edgecolor=edge,
            linewidth=lw,
            zorder=zorder,
        )
    )


def fig_text(
    fig: plt.Figure,
    x: float,
    y: float,
    text: str,
    *,
    size: float,
    weight: str = "normal",
    color: str = INK,
    ha: str = "left",
    va: str = "top",
    linespacing: float = 1.18,
) -> None:
    fig.text(
        x,
        y,
        text,
        fontsize=size,
        fontweight=weight,
        color=color,
        ha=ha,
        va=va,
        linespacing=linespacing,
    )


def find_row(rows: list[dict[str, str]], **where: str) -> dict[str, str]:
    for row in rows:
        if all(row.get(key) == value for key, value in where.items()):
            return row
    raise KeyError(f"Missing row: {where}")


def load_payload(audit_dir: Path) -> dict[str, object]:
    row_metrics = read_rows(audit_dir / "row_matched_holdout_cent020_et_metrics.csv")
    source_rows = read_rows(audit_dir / "source_fractions_holdout_cent020.csv")
    fixed_rows = read_rows(audit_dir / "fixed_validation_training_effect_auc.csv")
    wp80_rows = read_rows(audit_dir / "wp80_full_scorecache_fake_rate_summary.csv")

    auc = {
        sample: [
            as_float(find_row(row_metrics, sample=sample, et_label=et)["auc_binned"])
            for et in ET_LABELS
        ]
        for sample in SAMPLE_LABELS
    }
    gain = [
        auc["Jet12+20+30+40"][idx] - auc["Jet12+20"][idx]
        for idx in range(len(ET_LABELS))
    ]
    baseline_counts = [
        as_int(find_row(row_metrics, sample="Jet12+20", et_label=et)["inclusive_entries"])
        for et in ET_LABELS
    ]
    expanded_counts = [
        as_int(find_row(row_metrics, sample="Jet12+20+30+40", et_label=et)["inclusive_entries"])
        for et in ET_LABELS
    ]
    source_counts = {
        sample: {
            et: {
                source: as_int(
                    find_row(source_rows, branch=sample, et_bin=et, source_bin=source)[
                        "inclusive_background_raw"
                    ]
                )
                for source in SOURCE_LABELS
            }
            for et in ET_LABELS
        }
        for sample in ["Jet12+20", "Jet12+20+30+40"]
    }
    source_fractions = {
        et: {
            source: as_float(
                find_row(source_rows, branch="Jet12+20+30+40", et_bin=et, source_bin=source)[
                    "source_fraction_raw"
                ]
            )
            for source in SOURCE_LABELS
        }
        for et in ET_LABELS
    }

    fixed = {
        "Jet12+20": [
            as_float(
                find_row(
                    fixed_rows,
                    scope="direct_holdout_3x3",
                    validation_sample="Jet12+20",
                    model_training_sample=sample,
                )["auc_0_20"]
            )
            for sample in SAMPLE_LABELS
        ],
        "Jet12+20+30+40": [
            as_float(
                find_row(
                    fixed_rows,
                    scope="direct_holdout_3x3",
                    validation_sample="Jet12+20+30+40",
                    model_training_sample=sample,
                )["auc_0_20"]
            )
            for sample in SAMPLE_LABELS
        ],
    }
    wp80 = {
        row["sample"]: {
            "auc": as_float(row["auc_all_full_scorecache"]),
            "fake_rate": as_float(row["wp80_background_fake_rate"]),
        }
        for row in wp80_rows
    }
    return {
        "auc": auc,
        "gain": gain,
        "baseline_counts": baseline_counts,
        "expanded_counts": expanded_counts,
        "source_counts": source_counts,
        "source_fractions": source_fractions,
        "fixed": fixed,
        "wp80": wp80,
        "inputs": {
            "row_metrics": str(audit_dir / "row_matched_holdout_cent020_et_metrics.csv"),
            "source_fractions": str(audit_dir / "source_fractions_holdout_cent020.csv"),
            "fixed_validation": str(audit_dir / "fixed_validation_training_effect_auc.csv"),
            "wp80": str(audit_dir / "wp80_full_scorecache_fake_rate_summary.csv"),
            "report": str(audit_dir / "hard_inclusive_jet3040_mechanism_report.md"),
        },
    }


def draw_auc_panel(fig: plt.Figure, payload: dict[str, object]) -> None:
    ax = fig.add_axes([0.075, 0.585, 0.565, 0.225])
    x = np.arange(len(ET_LABELS))
    ax.axvspan(1.70, 2.30, color="#F3E8FF", alpha=0.45, zorder=0)
    for sample, marker in zip(SAMPLE_LABELS, ["o", "s", "^"]):
        values = payload["auc"][sample]
        ax.plot(
            x,
            values,
            marker=marker,
            markersize=7.5,
            linewidth=2.4,
            color=BRANCH_COLORS[sample],
            label=sample,
            zorder=4,
        )
    for idx, delta in enumerate(payload["gain"]):
        y = payload["auc"]["Jet12+20+30+40"][idx]
        ax.text(
            idx,
            min(y + 0.006, 0.903),
            f"+{delta:.3f}",
            ha="center",
            va="bottom",
            fontsize=11.0,
            fontweight="bold" if idx == 2 else "normal",
            color=INK,
        )
    ax.set_title("Where the gain appears", fontsize=15.5, fontweight="bold", color=INK, pad=8)
    ax.set_xticks(x)
    ax.set_xticklabels(ET_LABELS, fontsize=11.2)
    ax.set_ylabel("AUC, 0-20% holdout", fontsize=11.5)
    ax.set_xlabel("")
    ax.set_ylim(0.68, 0.925)
    ax.set_yticks([0.70, 0.75, 0.80, 0.85, 0.90])
    ax.grid(True, color=GRID, alpha=0.62, linewidth=0.8)
    ax.set_axisbelow(True)
    ax.legend(
        frameon=False,
        fontsize=10.2,
        ncol=3,
        loc="lower left",
        bbox_to_anchor=(0.015, 0.01),
        handlelength=2.2,
    )
    for spine in ax.spines.values():
        spine.set_color("#344054")
        spine.set_linewidth(0.85)
    ax.tick_params(direction="in", top=True, right=True, labelsize=10.6)


def draw_support_panel(fig: plt.Figure, payload: dict[str, object]) -> None:
    fig_text(fig, 0.075, 0.545, "Why the gain appears", size=17.0, weight="bold", color=INK)
    fig_text(
        fig,
        0.075,
        0.516,
        "Counts show the coverage gap; fractions show which hard-jet samples supply the expanded pool.",
        size=10.2,
        color=MUTED,
    )
    fig_text(fig, 0.075, 0.492, r"1. Total support falls at high $E_T$", size=12.6, weight="bold", color=INK)
    fig_text(fig, 0.448, 0.492, "2. Expanded-pool source mix", size=12.6, weight="bold", color=INK)

    legend_y = 0.455
    count_legend = [
        (0.075, BASELINE_FILL, BRANCH_COLORS["Jet12+20"], "Jet12+20 baseline"),
        (0.245, EXPANDED_FILL, BRANCH_COLORS["Jet12+20+30+40"], "+Jet30/40 expanded"),
    ]
    for x0, fill, edge, label in count_legend:
        fig.patches.append(
            Rectangle(
                (x0, legend_y - 0.0075),
                0.014,
                0.015,
                transform=fig.transFigure,
                facecolor=fill,
                edgecolor=edge,
                linewidth=1.0,
                zorder=3,
            )
        )
        fig_text(fig, x0 + 0.019, legend_y, label, size=9.0, color=INK, va="center")

    for idx, source in enumerate(SOURCE_LABELS):
        x0 = 0.450 + idx * 0.044
        fig.patches.append(
            Rectangle(
                (x0, legend_y - 0.007),
                0.010,
                0.014,
                transform=fig.transFigure,
                facecolor=SOURCE_COLORS[source],
                edgecolor="none",
                zorder=3,
            )
        )
        fig_text(fig, x0 + 0.013, legend_y, source, size=8.8, color=INK, va="center")

    ax_counts = fig.add_axes([0.075, 0.215, 0.320, 0.205])
    x = np.arange(len(ET_LABELS))
    width = 0.32
    baseline = np.array(payload["baseline_counts"], dtype=float)
    expanded = np.array(payload["expanded_counts"], dtype=float)
    ax_counts.bar(
        x - width / 2,
        baseline,
        width=width,
        color=BASELINE_FILL,
        edgecolor=BRANCH_COLORS["Jet12+20"],
        linewidth=1.15,
        label="Jet12+20 baseline",
        zorder=3,
    )
    ax_counts.bar(
        x + width / 2,
        expanded,
        width=width,
        color=EXPANDED_FILL,
        edgecolor=BRANCH_COLORS["Jet12+20+30+40"],
        linewidth=1.15,
        label="+Jet30/40 expanded",
        zorder=3,
    )
    for idx, (b_val, e_val) in enumerate(zip(baseline, expanded)):
        ax_counts.text(
            idx - width / 2,
            max(b_val * 1.14, b_val + 3),
            fmt_count(b_val),
            ha="center",
            va="bottom",
            fontsize=8.7,
            color=BRANCH_COLORS["Jet12+20"],
            fontweight="bold" if idx >= 2 else "normal",
        )
        ax_counts.text(
            idx + width / 2,
            e_val * 1.10,
            fmt_count(e_val),
            ha="center",
            va="bottom",
            fontsize=8.7,
            color=BRANCH_COLORS["Jet12+20+30+40"],
            fontweight="bold" if idx >= 2 else "normal",
        )
    ax_counts.set_yscale("log")
    ax_counts.set_ylim(8, 40000)
    ax_counts.set_xticks(x)
    ax_counts.set_xticklabels(ET_LABELS, fontsize=10.2)
    ax_counts.set_ylabel("background candidates\n10% holdout, log", fontsize=10.3)
    ax_counts.set_xlabel(r"cluster $E_T$ bin [GeV]", fontsize=10.6, labelpad=1)
    ax_counts.grid(True, axis="y", which="both", color=GRID, alpha=0.58, linewidth=0.75)
    ax_counts.set_axisbelow(True)
    for spine in ax_counts.spines.values():
        spine.set_color("#344054")
        spine.set_linewidth(0.85)
    ax_counts.tick_params(direction="in", top=True, right=True, labelsize=9.8)

    ax_frac = fig.add_axes([0.448, 0.215, 0.190, 0.205])
    bottom = np.zeros(len(ET_LABELS))
    for source in SOURCE_LABELS:
        vals = np.array(
            [100.0 * payload["source_fractions"][et][source] for et in ET_LABELS],
            dtype=float,
        )
        ax_frac.bar(
            x,
            vals,
            bottom=bottom,
            width=0.62,
            color=SOURCE_COLORS[source],
            edgecolor="white",
            linewidth=0.62,
            label=source,
            zorder=3,
        )
        if source == "Jet40":
            for idx, value in enumerate(vals):
                ax_frac.text(
                    idx,
                    bottom[idx] + value / 2,
                    f"{value:.0f}%",
                    ha="center",
                    va="center",
                    fontsize=8.8,
                    fontweight="bold",
                    color="white",
                )
        if source == "Jet30":
            value = vals[3]
            ax_frac.text(
                3,
                bottom[3] + value / 2,
                f"{value:.0f}%",
                ha="center",
                va="center",
                fontsize=8.5,
                fontweight="bold",
                color="white",
            )
        bottom += vals
    ax_frac.set_ylim(0, 100)
    ax_frac.set_yticks([0, 50, 100])
    ax_frac.set_xticks(x)
    ax_frac.set_xticklabels(ET_LABELS, fontsize=10.2)
    ax_frac.set_ylabel("source fraction [%]", fontsize=10.3)
    ax_frac.grid(True, axis="y", color=GRID, alpha=0.58, linewidth=0.75)
    ax_frac.set_axisbelow(True)
    for spine in ax_frac.spines.values():
        spine.set_color("#344054")
        spine.set_linewidth(0.85)
    ax_frac.tick_params(direction="in", top=True, right=True, labelsize=9.8)


def arrow_line(
    fig: plt.Figure,
    x: float,
    y: float,
    label: str,
    values: list[float],
    *,
    fmt: str,
    size: float = 12.4,
    value_size: float = 13.2,
) -> None:
    fig_text(fig, x, y, label, size=size, weight="bold", color=INK)
    start = x + 0.068
    step = 0.057
    for idx, value in enumerate(values):
        fig_text(
            fig,
            start + idx * step,
            y,
            fmt.format(value),
            size=value_size,
            weight="bold",
            color=BRANCH_COLORS[SAMPLE_LABELS[idx]],
        )
        if idx < 2:
            fig_text(fig, start + idx * step + 0.038, y, "→", size=value_size, weight="bold", color=MUTED)


def draw_right_column(fig: plt.Figure, payload: dict[str, object]) -> None:
    add_round_box(fig, (0.730, 0.505), (0.240, 0.335), SOFT_BLUE, edge="#9CC7F5", radius=0.010)
    fig_text(fig, 0.750, 0.807, "Control: fixed validation", size=14.7, weight="bold", color="#1C57B7")
    fig_text(
        fig,
        0.750,
        0.762,
        "Same rows; only the trained BDT changes.",
        size=11.8,
        color=INK,
    )
    fixed = payload["fixed"]
    fig_text(fig, 0.750, 0.715, "0-20% AUC:", size=12.0, weight="bold")
    fig_text(
        fig,
        0.764,
        0.681,
        f"Jet12+20 val: {fixed['Jet12+20'][0]:.4f} → {fixed['Jet12+20'][1]:.4f} → {fixed['Jet12+20'][2]:.4f}",
        size=11.1,
        color=INK,
    )
    fig_text(
        fig,
        0.764,
        0.646,
        f"Full val: {fixed['Jet12+20+30+40'][0]:.4f} → {fixed['Jet12+20+30+40'][1]:.4f} → {fixed['Jet12+20+30+40'][2]:.4f}",
        size=11.1,
        color=INK,
    )
    fig_text(
        fig,
        0.750,
        0.595,
        "Coverage drives the big visual gain;\ntraining gain remains real.",
        size=11.7,
        color=MUTED,
        linespacing=1.15,
    )

    add_round_box(fig, (0.730, 0.205), (0.240, 0.260), SOFT_PURPLE, edge="#C4B5E4", radius=0.010)
    fig_text(fig, 0.750, 0.430, "Physics payoff at WP80", size=14.8, weight="bold", color="#6A44A0")
    wp80 = payload["wp80"]
    arrow_line(
        fig,
        0.750,
        0.386,
        "AUC",
        [wp80[sample]["auc"] for sample in SAMPLE_LABELS],
        fmt="{:.3f}",
        size=12.3,
        value_size=13.0,
    )
    arrow_line(
        fig,
        0.750,
        0.332,
        "Fake",
        [wp80[sample]["fake_rate"] for sample in SAMPLE_LABELS],
        fmt="{:.3f}",
        size=12.3,
        value_size=13.0,
    )
    fig_text(
        fig,
        0.750,
        0.278,
        "At fixed 80% signal efficiency,\ninclusive-background acceptance falls\nby more than a factor of two.",
        size=11.8,
        color=INK,
        linespacing=1.13,
    )


def draw_slide(payload: dict[str, object], out_png: Path) -> None:
    plt.rcParams.update(
        {
            "font.family": "Times New Roman",
            "mathtext.fontset": "stix",
            "axes.titlesize": 15,
            "axes.labelsize": 11,
        }
    )
    fig = plt.figure(figsize=slide_figsize(SLIDE_DPI), dpi=SLIDE_DPI, facecolor=BG)

    fig_text(fig, 0.035, 0.972, TITLE, size=26.0, weight="bold")
    fig_text(fig, 0.037, 0.902, SUBTITLE, size=15.4, color=MUTED)

    add_round_box(fig, (0.035, 0.165), (0.645, 0.695), PANEL, edge="#D5DCE6", radius=0.012)
    draw_auc_panel(fig, payload)
    draw_support_panel(fig, payload)
    draw_right_column(fig, payload)

    add_round_box(fig, (0.035, 0.055), (0.935, 0.088), SOFT_ORANGE, edge="#F0B36A", radius=0.010)
    fig_text(
        fig,
        0.055,
        0.111,
        "Jet30 and Jet40 do not just add more events.",
        size=15.5,
        weight="bold",
        color="#9A4B12",
    )
    fig_text(
        fig,
        0.055,
        0.082,
        "They close a missing hard-background phase-space hole, and the BDT gain appears exactly where that hole was largest.",
        size=14.0,
        color=INK,
    )

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=SLIDE_DPI, facecolor=fig.get_facecolor(), bbox_inches=None)
    plt.close(fig)


def write_script(out_path: Path) -> None:
    out_path.write_text(
        """# WP_GammaJets_6_3_26 Slide 11 Script - Jet30/40 Gains: Further Diagnostic

Now I want to close the loop on why the Jet30 and Jet40 samples matter. The key point is that Jet12 plus Jet20 is not wrong. It is just incomplete in the part of phase space where the BDT needs the most help: high cluster transverse energy.

The top panel shows where the gain appears. In the 0 to 20 percent centrality holdout, the AUC improves only moderately in the 15 to 20 GeV bin, but the gain grows as we move higher in cluster E_T. The largest improvement is in the 25 to 30 GeV bin. The 30 to 35 GeV point still improves, but that baseline row has very few inclusive-jet candidates, so I treat the exact number there more cautiously.

The bottom section separates the mechanism into two simpler pieces. On the left, I am only comparing total inclusive-background support: Jet12 plus Jet20 baseline versus the expanded Jet12 plus Jet20 plus Jet30 plus Jet40 pool. On the right, I am only looking inside the expanded pool to show which source samples fill it. That separation is important. The left plot shows the phase-space hole, and the right plot shows that Jet30 and especially Jet40 are what fill it. At low E_T, Jet12 plus Jet20 already has meaningful support. At higher E_T, that support collapses, and in the 30 to 35 GeV bin the expanded holdout is about eighty percent Jet40 and nineteen percent Jet30.

The control box on the right separates two effects. When I fix the validation sample and only change which BDT was trained, the gains are smaller but still real. Jet30 gives the main training improvement, and Jet40 is mostly a saturation and highest-tail coverage check. That tells me the big visual gain in the score-shape slides is mostly because the validation and training now include the hard inclusive-jet background that was missing before, not because the plot normalization changed.

The second box gives the practical payoff. At the WP80 working point, meaning fixed eighty percent signal efficiency, the inclusive-background fake rate drops from about 0.42 to 0.26 to 0.17 as the hard inclusive ladder is added. So the interpretation is: Jet30 and Jet40 are not just more events. They close a hard-background phase-space hole, and the BDT improvement appears exactly where that hole was largest.
""",
        encoding="utf-8",
    )


def write_manifest(out_path: Path, payload: dict[str, object], out_png: Path, script_path: Path) -> None:
    manifest = {
        "schema": "THE8_HARD_INCLUSIVE_MECHANISM_CLOSURE_SLIDE_V4",
        "created": "2026-06-05",
        "slide_title": "Jet30/40 Gains: Further Diagnostic",
        "purpose": "Polished full-slide PNG candidate to close the slides 8-10 Jet30/Jet40 BDT mechanism story, revised for shorter diagnostic title, white working-point-deck background, bottom-panel digestibility, gutter spacing, and non-overlapping legends/annotations.",
        "inputs": payload["inputs"],
        "outputs": {
            "png": str(out_png),
            "speaker_script": str(script_path),
            "manifest": str(out_path),
        },
        "claim_scope": {
            "supported": [
                "Jet12+20 is incomplete at high cluster ET in the 0-20% holdout.",
                "Jet30/Jet40 dominate the added high-ET inclusive-background support.",
                "Fixed-validation controls show smaller but real training gain.",
                "Full score-cache WP80 fake rate falls 0.419 -> 0.259 -> 0.172.",
            ],
            "caveats": [
                "30-35 GeV Jet12+20 baseline row is statistically limited.",
                "No event-by-event HEPMC ancestry was available in the ML matrix.",
                "Candidate PNG has not been inserted into Google Slides.",
            ],
        },
    }
    out_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit-dir", type=Path, default=AUDIT_DIR)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    payload = load_payload(args.audit_dir)
    args.outdir.mkdir(parents=True, exist_ok=True)
    out_png = args.outdir / "hard_inclusive_jet3040_mechanism_closure_slide.png"
    script_path = args.outdir / "hard_inclusive_jet3040_mechanism_closure_script.md"
    manifest_path = args.outdir / "hard_inclusive_jet3040_mechanism_closure_manifest.json"
    draw_slide(payload, out_png)
    write_script(script_path)
    write_manifest(manifest_path, payload, out_png, script_path)
    print(out_png)
    print(script_path)
    print(manifest_path)


if __name__ == "__main__":
    main()
