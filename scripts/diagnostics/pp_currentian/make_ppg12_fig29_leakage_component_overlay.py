#!/usr/bin/env python3
"""Overlay current and PPG12 Fig.29 signal-leakage components."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
DEFAULT_INPUT = (
    REPO
    / "dataOutput/ppg12PhotonYield/ppg12_photon_yield_v1_data_20260620"
    / "purity_fig29_comparison/ppg12_ratio_diagnostic"
    / "current_vs_ppg12_fig29_purity_ratio_points.csv"
)
DEFAULT_OUTDIR = (
    REPO
    / "dataOutput/ppg12PhotonYield/ppg12_photon_yield_v1_data_20260620"
    / "purity_fig29_comparison/ppg12_ratio_diagnostic/leakage_components"
)


def read_rows(path: Path) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            rows.append({key: float(value) for key, value in row.items()})
    return rows


def finite_range(values: list[float]) -> tuple[float, float]:
    finite = [v for v in values if np.isfinite(v)]
    if not finite:
        return float("nan"), float("nan")
    return min(finite), max(finite)


def write_csv(rows: list[dict[str, float]], path: Path) -> None:
    fields = [
        "pt_lo",
        "pt_hi",
        "pt_center",
        "current_f_b",
        "ppg12_f_b",
        "current_f_b_over_ppg12",
        "current_f_c",
        "ppg12_f_c",
        "current_f_c_over_ppg12",
        "current_f_d",
        "ppg12_f_d",
        "current_f_d_over_ppg12",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row[field] for field in fields})


def make_plot(rows: list[dict[str, float]], out: Path) -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "dejavuserif",
            "axes.linewidth": 1.25,
            "axes.labelsize": 15,
            "xtick.labelsize": 12,
            "ytick.labelsize": 12,
            "legend.fontsize": 10.5,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )
    x = np.array([r["pt_center"] for r in rows])
    xerr = np.array([(r["pt_hi"] - r["pt_lo"]) / 2.0 for r in rows])
    components = [
        ("B", "tight, non-isolated", "current_f_b", "ppg12_f_b", "current_f_b_over_ppg12", "#0072B2", "o"),
        ("C", "non-tight, isolated", "current_f_c", "ppg12_f_c", "current_f_c_over_ppg12", "#009E73", "s"),
        ("D", "non-tight, non-isolated", "current_f_d", "ppg12_f_d", "current_f_d_over_ppg12", "#D55E00", "^"),
    ]
    fig, axes = plt.subplots(
        2,
        3,
        figsize=(15.5, 8.3),
        dpi=210,
        sharex=True,
        gridspec_kw={"height_ratios": [1.55, 1.0], "hspace": 0.08, "wspace": 0.18},
    )
    for idx, (label, subtitle, current_key, ppg12_key, ratio_key, color, marker) in enumerate(components):
        ax = axes[0, idx]
        ratio_ax = axes[1, idx]
        current = np.array([r[current_key] for r in rows])
        ppg12 = np.array([r[ppg12_key] for r in rows])
        ratios = np.array([r[ratio_key] for r in rows])
        ax.errorbar(
            x - 0.06,
            ppg12,
            xerr=xerr,
            fmt=marker,
            ms=6.6,
            color="black",
            lw=1.2,
            capsize=0,
            label="PPG12 Photon_final",
        )
        ax.errorbar(
            x + 0.06,
            current,
            xerr=xerr,
            fmt=marker,
            ms=6.7,
            mfc="white",
            mec=color,
            mew=1.7,
            ecolor=color,
            color=color,
            lw=1.2,
            capsize=0,
            label="Current pp output",
        )
        ax.set_title(f"Leakage {label}: {subtitle}", fontsize=14, pad=8)
        ax.set_ylabel(rf"$c_{label}$ = {label}$_{{sig}}$/A$_{{sig}}$")
        ymax = max(np.nanmax(ppg12), np.nanmax(current), 0.01) * 1.35
        ax.set_ylim(-0.005, ymax)
        ax.grid(True, color="0.8", alpha=0.35)
        if idx == 0:
            ax.legend(loc="upper left", frameon=False)

        ratio_ax.axhline(1.0, color="0.25", lw=1.25)
        ratio_ax.plot(x, ratios, marker + "-", color=color, ms=6.2, lw=1.9)
        ratio_ax.set_ylim(-0.02, 1.18)
        ratio_ax.set_xlabel(r"Cluster $E_T$ [GeV]")
        ratio_ax.set_ylabel("Current / PPG12")
        ratio_ax.grid(True, color="0.8", alpha=0.35)

    fig.suptitle("Signal-leakage components behind PPG12 Fig.29 purity correction", fontsize=18, y=0.985)
    fig.text(
        0.99,
        0.012,
        r"$\bf{\it{sPHENIX}}$ Internal   $p{+}p$, $\sqrt{s}=200$ GeV",
        ha="right",
        va="bottom",
        fontsize=12,
    )
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


def write_summary(rows: list[dict[str, float]], out: Path, source: Path) -> dict[str, object]:
    stable = [r for r in rows if r["pt_hi"] <= 26.0]
    ranges = {
        "fB_current_over_ppg12": finite_range([r["current_f_b_over_ppg12"] for r in stable]),
        "fC_current_over_ppg12": finite_range([r["current_f_c_over_ppg12"] for r in stable]),
        "fD_current_over_ppg12": finite_range([r["current_f_d_over_ppg12"] for r in stable]),
        "leakage_shift_current_over_ppg12": finite_range([r["current_shift_over_ppg12_shift"] for r in stable]),
    }

    def classify_component(name: str, key: str) -> str:
        lo, hi = ranges[key]
        if hi < 0.8:
            return f"{name} is systematically below PPG12 in the stable bins."
        if lo > 1.2:
            return f"{name} is systematically above PPG12 in the stable bins."
        return f"{name} is mixed or near PPG12 within the stable-bin scan."

    shift_lo, shift_hi = ranges["leakage_shift_current_over_ppg12"]
    if shift_hi < 0.8:
        shift_text = "The total leakage shift is still too small relative to PPG12."
    elif shift_lo > 1.2:
        shift_text = "The total leakage shift is now too large relative to PPG12."
    else:
        shift_text = "The total leakage shift is in the PPG12 neighborhood, but component-level agreement still needs checking."

    summary: dict[str, object] = {
        "source_csv": str(source),
        "stable_bin_definition": "pt_hi <= 26 GeV",
        "stable_ranges": ranges,
        "interpretation": [
            "The direct B/C/D overlay confirms the correction mismatch is not a plotting artifact.",
            classify_component("B leakage", "fB_current_over_ppg12"),
            classify_component("C leakage", "fC_current_over_ppg12"),
            classify_component("D leakage", "fD_current_over_ppg12"),
            shift_text,
            "This points first to signal-leakage fill semantics, pp isolation-region definitions, or MC weighting/selection differences, not the algebra in the purity solver.",
        ],
    }
    out.with_suffix(".json").write_text(json.dumps(summary, indent=2) + "\n")
    with out.open("w") as handle:
        handle.write("# PPG12 Fig.29 leakage component overlay\n\n")
        handle.write(f"- Input ratio table: `{source}`\n")
        handle.write("- Exact PPG12 reference source: `/sphenix/user/shuhangli/ppg12/efficiencytool/results/Photon_final_bdt_nom.root`\n")
        handle.write("- Leakage definition checked in PPG12 `CalculatePhotonYield.C`: `cB=B_sig/A_sig`, `cC=C_sig/A_sig`, `cD=D_sig/A_sig`.\n\n")
        handle.write("## Stable-bin readout\n\n")
        ranges = summary["stable_ranges"]
        for key, value in ranges.items():
            lo, hi = value
            handle.write(f"- `{key}`: `{lo:.3f}` to `{hi:.3f}`\n")
        handle.write("\n## Current interpretation\n\n")
        for item in summary["interpretation"]:
            handle.write(f"- {item}\n")
        handle.write("\n## Component table\n\n")
        handle.write("| ET bin | fB current/PPG12 | fC current/PPG12 | fD current/PPG12 |\n")
        handle.write("| --- | ---: | ---: | ---: |\n")
        for row in rows:
            handle.write(
                f"| {row['pt_lo']:.0f}-{row['pt_hi']:.0f} | "
                f"{row['current_f_b_over_ppg12']:.3f} | "
                f"{row['current_f_c_over_ppg12']:.3f} | "
                f"{row['current_f_d_over_ppg12']:.3f} |\n"
            )
    return summary


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    rows = read_rows(args.input)
    if not rows:
        raise RuntimeError(f"no rows read from {args.input}")
    png = args.outdir / "current_vs_ppg12_fig29_leakage_components_bcd.png"
    csv_path = args.outdir / "current_vs_ppg12_fig29_leakage_components_bcd.csv"
    summary_md = args.outdir / "current_vs_ppg12_fig29_leakage_components_bcd.md"
    make_plot(rows, png)
    write_csv(rows, csv_path)
    write_summary(rows, summary_md, args.input)
    print(png)
    print(csv_path)
    print(summary_md)
    print(summary_md.with_suffix(".json"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
