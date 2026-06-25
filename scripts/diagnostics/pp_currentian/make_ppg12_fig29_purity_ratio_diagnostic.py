#!/usr/bin/env python3
"""Compare current pp photon-yield purity points to PPG12 Photon_final graphs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
DEFAULT_OUTDIR = (
    REPO
    / "dataOutput/ppg12PhotonYield/ppg12_photon_yield_v1_data_20260620"
    / "purity_fig29_comparison/ppg12_ratio_diagnostic"
)
DEFAULT_CURRENT = (
    REPO
    / "dataOutput/ppg12PhotonYield/ppg12_photon_yield_v1_data_20260620"
    / "purity_fig29_comparison/current_pp_photon_yield_purity_points.csv"
)
DEFAULT_PPG12 = DEFAULT_OUTDIR / "ppg12_photon_final_bdt_nom_extract.csv"


def read_current(path: Path) -> dict[float, dict[str, float]]:
    out: dict[float, dict[str, float]] = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            x = 0.5 * (float(row["pt_lo"]) + float(row["pt_hi"]))
            out[x] = {
                "pt_lo": float(row["pt_lo"]),
                "pt_hi": float(row["pt_hi"]),
                "a": float(row["a"]),
                "b": float(row["b"]),
                "c": float(row["c"]),
                "d": float(row["d"]),
                "raw": float(row["raw"]),
                "raw_err": float(row["raw_err"]),
                "corrected": float(row["corrected"]),
                "corrected_err": float(row["corrected_err"]),
                "f_b": float(row["f_b"]),
                "f_c": float(row["f_c"]),
                "f_d": float(row["f_d"]),
            }
    return out


def read_ppg12(path: Path) -> tuple[dict[str, dict[float, dict[str, float]]], dict[str, dict[float, dict[str, float]]]]:
    graphs: dict[str, dict[float, dict[str, float]]] = {}
    hists: dict[str, dict[float, dict[str, float]]] = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            if row["source"] != "data":
                continue
            kind = row["kind"]
            name = row["name"]
            if kind == "GRAPH":
                x = float(row["x"])
                graphs.setdefault(name, {})[x] = {
                    "y": float(row["y"]),
                    "ex_low": float(row["ex_low"]),
                    "ex_high": float(row["ex_high"]),
                    "ey_low": float(row["ey_low"]),
                    "ey_high": float(row["ey_high"]),
                }
            elif kind == "HIST":
                x = float(row["x"])
                hists.setdefault(name, {})[x] = {
                    "bin_lo": float(row["bin_lo"]),
                    "bin_hi": float(row["bin_hi"]),
                    "content": float(row["content"]),
                    "error": float(row["error"]),
                }
    return graphs, hists


def ratio(num: float, den: float) -> float:
    return num / den if np.isfinite(num) and np.isfinite(den) and den != 0.0 else float("nan")


def build_rows(
    current: dict[float, dict[str, float]],
    graphs: dict[str, dict[float, dict[str, float]]],
    hists: dict[str, dict[float, dict[str, float]]],
) -> list[dict[str, float]]:
    raw = graphs.get("gpurity", {})
    corr = graphs.get("gpurity_leak", {})
    rows: list[dict[str, float]] = []
    for x in sorted(current):
        if x not in raw or x not in corr:
            continue
        c = current[x]
        p_raw = raw[x]
        p_corr = corr[x]
        this_shift = c["corrected"] - c["raw"]
        ppg12_shift = p_corr["y"] - p_raw["y"]
        ppg12_fb = hists.get("h_leak_B", {}).get(x, {}).get("content", float("nan"))
        ppg12_fc = hists.get("h_leak_C", {}).get(x, {}).get("content", float("nan"))
        ppg12_fd = hists.get("h_leak_D", {}).get(x, {}).get("content", float("nan"))
        rows.append(
            {
                "pt_center": x,
                "pt_lo": c["pt_lo"],
                "pt_hi": c["pt_hi"],
                "current_a": c["a"],
                "current_b": c["b"],
                "current_c": c["c"],
                "current_d": c["d"],
                "current_raw": c["raw"],
                "current_raw_err": c["raw_err"],
                "ppg12_raw": p_raw["y"],
                "ppg12_raw_err": 0.5 * (p_raw["ey_low"] + p_raw["ey_high"]),
                "current_raw_over_ppg12": ratio(c["raw"], p_raw["y"]),
                "current_corrected": c["corrected"],
                "current_corrected_err": c["corrected_err"],
                "ppg12_corrected": p_corr["y"],
                "ppg12_corrected_err": 0.5 * (p_corr["ey_low"] + p_corr["ey_high"]),
                "current_corrected_over_ppg12": ratio(c["corrected"], p_corr["y"]),
                "current_leakage_shift": this_shift,
                "ppg12_leakage_shift": ppg12_shift,
                "current_shift_over_ppg12_shift": ratio(this_shift, ppg12_shift),
                "current_f_b": c["f_b"] if "f_b" in c else float("nan"),
                "ppg12_f_b": ppg12_fb,
                "current_f_b_over_ppg12": ratio(c["f_b"] if "f_b" in c else float("nan"), ppg12_fb),
                "current_f_c": c["f_c"] if "f_c" in c else float("nan"),
                "ppg12_f_c": ppg12_fc,
                "current_f_c_over_ppg12": ratio(c["f_c"] if "f_c" in c else float("nan"), ppg12_fc),
                "current_f_d": c["f_d"] if "f_d" in c else float("nan"),
                "ppg12_f_d": ppg12_fd,
                "current_f_d_over_ppg12": ratio(c["f_d"] if "f_d" in c else float("nan"), ppg12_fd),
            }
        )
    return rows


def write_csv(rows: list[dict[str, float]], path: Path) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def make_plot(rows: list[dict[str, float]], out: Path) -> None:
    x = np.array([r["pt_center"] for r in rows])
    xerr = np.array([(r["pt_hi"] - r["pt_lo"]) / 2.0 for r in rows])
    raw_this = np.array([r["current_raw"] for r in rows])
    raw_ref = np.array([r["ppg12_raw"] for r in rows])
    corr_this = np.array([r["current_corrected"] for r in rows])
    corr_ref = np.array([r["ppg12_corrected"] for r in rows])
    raw_ratio = np.array([r["current_raw_over_ppg12"] for r in rows])
    corr_ratio = np.array([r["current_corrected_over_ppg12"] for r in rows])
    shift_ratio = np.array([r["current_shift_over_ppg12_shift"] for r in rows])

    fig, axes = plt.subplots(
        3,
        1,
        figsize=(10.5, 10.0),
        sharex=True,
        gridspec_kw={"height_ratios": [1.55, 1.0, 1.0], "hspace": 0.08},
    )
    ax = axes[0]
    ax.errorbar(
        x - 0.08,
        raw_ref,
        xerr=xerr,
        yerr=[r["ppg12_raw_err"] for r in rows],
        fmt="o",
        color="black",
        ms=6,
        label="PPG12 Photon_final raw",
    )
    ax.errorbar(
        x + 0.08,
        raw_this,
        xerr=xerr,
        yerr=[r["current_raw_err"] for r in rows],
        fmt="s",
        mfc="white",
        mec="#1f77b4",
        ecolor="#1f77b4",
        color="#1f77b4",
        ms=6,
        label="Current pp raw",
    )
    ax.errorbar(
        x - 0.08,
        corr_ref,
        xerr=xerr,
        yerr=[r["ppg12_corrected_err"] for r in rows],
        fmt="o",
        color="#d62728",
        ms=6,
        label="PPG12 Photon_final leakage corrected",
    )
    ax.errorbar(
        x + 0.08,
        corr_this,
        xerr=xerr,
        yerr=[r["current_corrected_err"] for r in rows],
        fmt="s",
        mfc="white",
        mec="#2ca02c",
        ecolor="#2ca02c",
        color="#2ca02c",
        ms=6,
        label="Current pp leakage corrected",
    )
    ax.set_ylabel("Purity")
    ax.set_ylim(0.0, 1.18)
    ax.grid(True, alpha=0.2)
    ax.legend(loc="upper left", frameon=False, ncol=2, fontsize=10)
    ax.text(0.985, 0.08, "sPHENIX Internal\npp, $\\sqrt{s}=200$ GeV", transform=ax.transAxes, ha="right", va="bottom", fontsize=12)

    axes[1].axhline(1.0, color="0.2", lw=1.2)
    axes[1].plot(x, raw_ratio, "o-", color="black", label="raw current / PPG12")
    axes[1].plot(x, corr_ratio, "s-", color="#d62728", label="corrected current / PPG12")
    axes[1].set_ylabel("Ratio")
    axes[1].set_ylim(0.2, 1.25)
    axes[1].grid(True, alpha=0.2)
    axes[1].legend(loc="lower left", frameon=False, fontsize=10)

    axes[2].axhline(1.0, color="0.2", lw=1.2)
    axes[2].plot(x, shift_ratio, "D-", color="#9467bd", label="leakage shift current / PPG12")
    axes[2].set_ylabel("Shift ratio")
    axes[2].set_xlabel(r"Cluster $E_T$ [GeV]")
    axes[2].set_ylim(0.0, 1.05)
    axes[2].grid(True, alpha=0.2)
    axes[2].legend(loc="upper left", frameon=False, fontsize=10)
    axes[2].set_xlim(9.4, 36.6)

    fig.suptitle("Current pp purity vs PPG12 Photon_final Fig. 29 source graphs", fontsize=16, y=0.985)
    fig.savefig(out, dpi=180, bbox_inches="tight")
    plt.close(fig)


def write_summary(rows: list[dict[str, float]], out: Path, ppg12_source: Path, current_source: Path) -> None:
    raw_ratios = [r["current_raw_over_ppg12"] for r in rows if np.isfinite(r["current_raw_over_ppg12"])]
    corr_ratios = [r["current_corrected_over_ppg12"] for r in rows if np.isfinite(r["current_corrected_over_ppg12"])]
    shift_ratios = [r["current_shift_over_ppg12_shift"] for r in rows if np.isfinite(r["current_shift_over_ppg12_shift"])]
    with out.open("w") as handle:
        handle.write("# PPG12 Fig.29 purity ratio diagnostic\n\n")
        handle.write(f"- PPG12 extracted graph source: `{ppg12_source}`\n")
        handle.write(f"- Current pp source: `{current_source}`\n")
        handle.write("- PPG12 SDCC ROOT source: `/sphenix/user/shuhangli/ppg12/efficiencytool/results/Photon_final_bdt_nom.root`\n\n")
        handle.write("## Readout\n\n")
        handle.write(f"- Raw current/PPG12 ratio range: `{min(raw_ratios):.3f}` to `{max(raw_ratios):.3f}`.\n")
        handle.write(f"- Leakage-corrected current/PPG12 ratio range: `{min(corr_ratios):.3f}` to `{max(corr_ratios):.3f}`.\n")
        handle.write(f"- Leakage-shift current/PPG12 ratio range: `{min(shift_ratios):.3f}` to `{max(shift_ratios):.3f}`.\n")
        handle.write("- Interpretation: the raw purity is already below PPG12 in most bins, so the discrepancy is not only the leakage correction.\n")
        handle.write("- The leakage correction in the current pp product is also much smaller than PPG12 in every stable bin, so leakage inputs/normalization are a second, independent difference.\n")
        handle.write("- High-ET current bins have small A-region counts, which explains the larger statistical errors and unstable top bins.\n\n")
        handle.write("## Bin table\n\n")
        handle.write("| ET bin | raw ratio | corrected ratio | leakage shift ratio | current A |\n")
        handle.write("| --- | ---: | ---: | ---: | ---: |\n")
        for r in rows:
            handle.write(
                f"| {r['pt_lo']:.0f}-{r['pt_hi']:.0f} | "
                f"{r['current_raw_over_ppg12']:.3f} | "
                f"{r['current_corrected_over_ppg12']:.3f} | "
                f"{r['current_shift_over_ppg12_shift']:.3f} | "
                f"{r['current_a']:.0f} |\n"
            )
        handle.write("\n## Leakage fraction comparison\n\n")
        handle.write("| ET bin | fB current/PPG12 | fC current/PPG12 | fD current/PPG12 |\n")
        handle.write("| --- | ---: | ---: | ---: |\n")
        for r in rows:
            handle.write(
                f"| {r['pt_lo']:.0f}-{r['pt_hi']:.0f} | "
                f"{r['current_f_b_over_ppg12']:.3g} | "
                f"{r['current_f_c_over_ppg12']:.3g} | "
                f"{r['current_f_d_over_ppg12']:.3g} |\n"
            )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--current", type=Path, default=DEFAULT_CURRENT)
    parser.add_argument("--ppg12", type=Path, default=DEFAULT_PPG12)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    current = read_current(args.current)
    graphs, hists = read_ppg12(args.ppg12)
    rows = build_rows(current, graphs, hists)
    if not rows:
        raise RuntimeError("no common current/PPG12 purity bins found")
    write_csv(rows, args.outdir / "current_vs_ppg12_fig29_purity_ratio_points.csv")
    make_plot(rows, args.outdir / "current_vs_ppg12_fig29_purity_ratio_diagnostic.png")
    write_summary(
        rows,
        args.outdir / "current_vs_ppg12_fig29_purity_ratio_diagnostic.md",
        args.ppg12,
        args.current,
    )


if __name__ == "__main__":
    main()
