#!/usr/bin/env python3
"""Plot raw ABCD yields from the current THE-42 AuAu interim diagnostic CSV."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-codex")

import matplotlib.pyplot as plt


REPO = Path(__file__).resolve().parents[3]
DEFAULT_INPUT = (
    REPO
    / "dataOutput/ppg12TableQA/THE42_current_auau_interimEligibleDiagnostic_20260614"
    / "raw_abcd_purity/current_auau_interim_raw_abcd_purity_points.csv"
)
DEFAULT_OUTDIR = DEFAULT_INPUT.parents[1] / "raw_abcd_yields"

REGIONS = [
    ("A", "A: isolated + tight", "#111111", "o"),
    ("B", "B: non-isolated + tight", "#0072B2", "s"),
    ("C", "C: isolated + non-tight", "#D55E00", "^"),
    ("D", "D: non-isolated + non-tight", "#009E73", "D"),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT, help="Raw ABCD points CSV")
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR, help="Output directory")
    parser.add_argument("--tag", default="current_auau_interim_raw_abcd_yield_overlay")
    return parser.parse_args()


def read_rows(path: Path) -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = []
    with path.open() as handle:
        for row in csv.DictReader(handle):
            parsed: dict[str, float | str] = {}
            for key, value in row.items():
                if key in {"source", "centrality"}:
                    parsed[key] = value
                else:
                    parsed[key] = float(value) if value else math.nan
            rows.append(parsed)
    return rows


def select(rows: list[dict[str, float | str]], source_substr: str, centrality: str) -> list[dict[str, float | str]]:
    out = [
        row
        for row in rows
        if source_substr in str(row["source"]) and str(row["centrality"]) == centrality
    ]
    return sorted(out, key=lambda row: float(row["pt_mid"]))


def positive_points(rows: list[dict[str, float | str]], region: str) -> tuple[list[float], list[float], list[float], list[float]]:
    x: list[float] = []
    xerr: list[float] = []
    y: list[float] = []
    yerr: list[float] = []
    for row in rows:
        value = float(row[region])
        if value <= 0.0:
            continue
        lo = float(row["pt_low"])
        hi = float(row["pt_high"])
        x.append(float(row["pt_mid"]))
        xerr.append(0.5 * (hi - lo))
        y.append(value)
        err = float(row.get(f"{region}_error", math.nan))
        yerr.append(err if math.isfinite(err) and err > 0.0 else math.sqrt(value))
    return x, xerr, y, yerr


def max_positive(rows: list[dict[str, float | str]], regions: list[str]) -> float:
    vals = []
    for row in rows:
        for region in regions:
            val = float(row[region])
            if val > 0.0:
                vals.append(val)
    return max(vals) if vals else 1.0


def draw_region_yields(ax, rows: list[dict[str, float | str]], title: str) -> None:
    for region, label, color, marker in REGIONS:
        x, xerr, y, yerr = positive_points(rows, region)
        if not x:
            continue
        ax.errorbar(
            x,
            y,
            xerr=xerr,
            yerr=yerr,
            label=label,
            color=color,
            marker=marker,
            markersize=6.8,
            markerfacecolor="white" if region != "A" else color,
            markeredgewidth=1.7,
            linestyle="none",
            elinewidth=1.8,
            capsize=2.5,
        )
    ymax = max_positive(rows, ["A", "B", "C", "D"])
    ax.set_title(title, fontsize=15, fontweight="bold", pad=10)
    ax.set_yscale("log")
    ax.set_ylim(0.65, ymax * 2.8)
    ax.set_xlim(14.2, 35.8)
    ax.set_xlabel(r"Photon candidate $E_T$ [GeV]", fontsize=13)
    ax.set_ylabel("Raw candidate count", fontsize=13)
    ax.grid(True, which="major", axis="y", color="#d8d8d8", linewidth=0.9)
    ax.grid(True, which="minor", axis="y", color="#eeeeee", linewidth=0.5)
    ax.tick_params(labelsize=11, top=True, right=True)


def draw_a_overlay(
    ax,
    auau_rows: list[dict[str, float | str]],
    pp_rows: list[dict[str, float | str]],
) -> None:
    for rows, label, color, marker, open_marker in [
        (auau_rows, "AuAu 50-80%, A region", "#D55E00", "o", False),
        (pp_rows, "pp checkpoint, A region", "#0072B2", "s", True),
    ]:
        x, xerr, y, yerr = positive_points(rows, "A")
        ax.errorbar(
            x,
            y,
            xerr=xerr,
            yerr=yerr,
            label=label,
            color=color,
            marker=marker,
            markersize=7.2,
            markerfacecolor="white" if open_marker else color,
            markeredgewidth=1.8,
            linestyle="none",
            elinewidth=1.9,
            capsize=2.5,
        )
    ymax = max(max_positive(auau_rows, ["A"]), max_positive(pp_rows, ["A"]))
    ax.set_title("Tight A-region raw counts", fontsize=15, fontweight="bold", pad=10)
    ax.set_yscale("log")
    ax.set_ylim(0.65, ymax * 2.8)
    ax.set_xlim(14.2, 35.8)
    ax.set_xlabel(r"Photon candidate $E_T$ [GeV]", fontsize=13)
    ax.set_ylabel("Raw A-region count", fontsize=13)
    ax.grid(True, which="major", axis="y", color="#d8d8d8", linewidth=0.9)
    ax.grid(True, which="minor", axis="y", color="#eeeeee", linewidth=0.5)
    ax.tick_params(labelsize=11, top=True, right=True)
    ax.text(
        0.03,
        0.05,
        "Raw counts only\nnot luminosity-normalized",
        transform=ax.transAxes,
        fontsize=10.5,
        color="#444444",
        va="bottom",
    )


def make_plot(rows: list[dict[str, float | str]], out_png: Path, input_csv: Path) -> None:
    auau_020 = select(rows, "AuAu", "0_20")
    auau_5080 = select(rows, "AuAu", "50_80")
    pp = select(rows, "pp data", "pp")

    if not auau_020 or not auau_5080 or not pp:
        raise RuntimeError("Missing required 0-20, 50-80, or pp rows in input CSV")

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.linewidth": 1.05,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "savefig.facecolor": "white",
        }
    )
    fig, axes = plt.subplots(1, 3, figsize=(19.2, 6.4), constrained_layout=False)
    fig.patch.set_facecolor("white")

    draw_region_yields(axes[0], auau_020, "AuAu 0-20% raw ABCD yields")
    draw_region_yields(axes[1], auau_5080, "AuAu 50-80% raw ABCD yields")
    draw_a_overlay(axes[2], auau_5080, pp)

    for ax in axes:
        ax.legend(frameon=False, fontsize=10.8, loc="upper right")

    fig.text(0.028, 0.955, r"$\bf{\it{sPHENIX}}$ Internal", fontsize=17, va="top")
    fig.text(
        0.028,
        0.902,
        "Current THE-42 AuAu interim eligible merge: 331 per-run ROOTs; pp checkpoint shown for peripheral A-region context",
        fontsize=13.5,
        va="top",
    )
    fig.text(
        0.028,
        0.858,
        r"15 < $E_T$ < 35 GeV, $\Delta R=0.3$ sliding isolation; zero-count bins omitted on log axes",
        fontsize=12.2,
        color="#444444",
        va="top",
    )
    fig.text(0.99, 0.02, f"Input: {input_csv.name}", fontsize=8.5, color="#777777", ha="right")
    fig.subplots_adjust(left=0.06, right=0.985, bottom=0.14, top=0.78, wspace=0.22)

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=220)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    rows = read_rows(args.input)
    out_png = args.outdir / f"{args.tag}.png"
    make_plot(rows, out_png, args.input)

    manifest = {
        "schema": "THE42_CURRENT_AUAU_INTERIM_RAW_ABCD_YIELD_OVERLAY_V1",
        "input_csv": str(args.input),
        "output_png": str(out_png),
        "note": (
            "Raw count diagnostic. AuAu uses current THE-42 interim eligible per-run merges; "
            "pp overlay is raw pp checkpoint A-region count and is not luminosity-normalized. "
            "Yield points are marker-only by policy; no connecting lines are drawn."
        ),
        "zero_count_handling": "Omitted from log-y axes.",
        "panels": [
            "AuAu 0-20% ABCD raw counts",
            "AuAu 50-80% ABCD raw counts",
            "A-region raw counts: AuAu 50-80% vs pp checkpoint",
        ],
    }
    manifest_path = args.outdir / f"{args.tag}_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(out_png)
    print(manifest_path)


if __name__ == "__main__":
    main()
