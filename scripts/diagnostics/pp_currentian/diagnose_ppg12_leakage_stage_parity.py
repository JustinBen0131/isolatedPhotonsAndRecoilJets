#!/usr/bin/env python3
"""Compare current RecoilJets signal-leakage stage histograms to PPG12.

This diagnostic deliberately sits below the purity formula.  It compares the
raw signal ABCD histograms that feed PPG12 Fig.29 so the remaining C/D mismatch
can be localized to a tight/non-tight population issue, an isolation split, or
normalization/weighting.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import ROOT


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
DEFAULT_CURRENT_ROOT = (
    REPO
    / "InputFiles/pp24/ppg12_photon_yield_v1_signal_sim_ppg12mix_combined_20260625"
    / "sim/jetMinPtScan_dphiScan_vz60_isoR40_isSliding_preselectionNewPPG12_tightNewPPG12_nonTightNewPPG12"
    / "photonJet5and10and20ppg12mixmerged_SIM"
    / "RecoilJets_photonjet5plus10plus20_ppg12mix_MERGED.root"
)
DEFAULT_PPG12_STAGE = (
    REPO
    / "dataOutput/ppg12PhotonYield/ppg12_photon_yield_v1_data_20260620"
    / "purity_fig29_comparison/ppg12_ratio_diagnostic/leakage_components"
    / "ppg12_mc_efficiency_bdt_nom_stage_extract.csv"
)
DEFAULT_OUTDIR = (
    REPO
    / "dataOutput/ppg12PhotonYield/ppg12_photon_yield_v1_data_20260620"
    / "purity_fig29_comparison/ppg12_ratio_diagnostic/leakage_components"
)

HISTS = {
    "all": "h_all_cluster_signal_0",
    "tight": "h_tight_cluster_signal_0",
    "A": "h_tight_iso_cluster_signal_0",
    "B": "h_tight_noniso_cluster_signal_0",
    "C": "h_nontight_iso_cluster_signal_0",
    "D": "h_nontight_noniso_cluster_signal_0",
}


def safe_div(num: float, den: float) -> float:
    return num / den if den else float("nan")


def read_ppg12_stage(path: Path) -> dict[str, dict[tuple[float, float], tuple[float, float]]]:
    out: dict[str, dict[tuple[float, float], tuple[float, float]]] = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            hist = row["hist"]
            lo_hi = (float(row["pt_lo"]), float(row["pt_hi"]))
            out.setdefault(hist, {})[lo_hi] = (float(row["value"]), float(row["error"]))
    return out


def get_hist(file_handle: ROOT.TFile, hist_name: str) -> ROOT.TH1:
    for key in (f"SIM/{hist_name}", hist_name):
        hist = file_handle.Get(key)
        if hist:
            return hist
    raise RuntimeError(f"missing histogram {hist_name} in {file_handle.GetName()}")


def read_current_root(path: Path) -> dict[str, dict[tuple[float, float], tuple[float, float]]]:
    file_handle = ROOT.TFile.Open(str(path))
    if not file_handle or file_handle.IsZombie():
        raise RuntimeError(f"cannot open current ROOT file: {path}")
    out: dict[str, dict[tuple[float, float], tuple[float, float]]] = {}
    for hist_name in HISTS.values():
        hist = get_hist(file_handle, hist_name)
        bins: dict[tuple[float, float], tuple[float, float]] = {}
        axis = hist.GetXaxis()
        for ibin in range(1, hist.GetNbinsX() + 1):
            lo_hi = (float(axis.GetBinLowEdge(ibin)), float(axis.GetBinUpEdge(ibin)))
            bins[lo_hi] = (float(hist.GetBinContent(ibin)), float(hist.GetBinError(ibin)))
        out[hist_name] = bins
    file_handle.Close()
    return out


def build_rows(
    ppg12: dict[str, dict[tuple[float, float], tuple[float, float]]],
    current: dict[str, dict[tuple[float, float], tuple[float, float]]],
) -> list[dict[str, float]]:
    bins = sorted(ppg12[HISTS["A"]])
    rows: list[dict[str, float]] = []
    for lo, hi in bins:
        row: dict[str, float] = {
            "pt_lo": lo,
            "pt_hi": hi,
            "pt_center": 0.5 * (lo + hi),
        }
        for short, hist_name in HISTS.items():
            ppg12_value, ppg12_error = ppg12[hist_name][(lo, hi)]
            current_value, current_error = current[hist_name][(lo, hi)]
            row[f"ppg12_{short}"] = ppg12_value
            row[f"ppg12_{short}_err"] = ppg12_error
            row[f"current_{short}"] = current_value
            row[f"current_{short}_err"] = current_error
            row[f"current_over_ppg12_{short}"] = safe_div(current_value, ppg12_value)

        for prefix in ("ppg12", "current"):
            a = row[f"{prefix}_A"]
            row[f"{prefix}_B_over_A"] = safe_div(row[f"{prefix}_B"], a)
            row[f"{prefix}_C_over_A"] = safe_div(row[f"{prefix}_C"], a)
            row[f"{prefix}_D_over_A"] = safe_div(row[f"{prefix}_D"], a)
            row[f"{prefix}_nt_over_A"] = safe_div(row[f"{prefix}_C"] + row[f"{prefix}_D"], a)
            row[f"{prefix}_D_share_nt"] = safe_div(row[f"{prefix}_D"], row[f"{prefix}_C"] + row[f"{prefix}_D"])
            row[f"{prefix}_A_over_all"] = safe_div(row[f"{prefix}_A"], row[f"{prefix}_all"])
            row[f"{prefix}_nt_over_all"] = safe_div(row[f"{prefix}_C"] + row[f"{prefix}_D"], row[f"{prefix}_all"])

        for key in ("B_over_A", "C_over_A", "D_over_A", "nt_over_A", "D_share_nt", "A_over_all", "nt_over_all"):
            row[f"current_over_ppg12_{key}"] = safe_div(row[f"current_{key}"], row[f"ppg12_{key}"])
        rows.append(row)
    return rows


def finite_range(rows: list[dict[str, float]], key: str) -> tuple[float, float]:
    values = [row[key] for row in rows if math.isfinite(row[key])]
    return (min(values), max(values)) if values else (float("nan"), float("nan"))


def write_csv(rows: list[dict[str, float]], path: Path) -> None:
    fields = [
        "pt_lo",
        "pt_hi",
        "pt_center",
        "current_over_ppg12_all",
        "current_over_ppg12_tight",
        "current_over_ppg12_A",
        "current_over_ppg12_B",
        "current_over_ppg12_C",
        "current_over_ppg12_D",
        "current_over_ppg12_B_over_A",
        "current_over_ppg12_C_over_A",
        "current_over_ppg12_D_over_A",
        "current_over_ppg12_nt_over_A",
        "current_over_ppg12_D_share_nt",
        "current_over_ppg12_A_over_all",
        "current_over_ppg12_nt_over_all",
        "ppg12_A",
        "current_A",
        "ppg12_B",
        "current_B",
        "ppg12_C",
        "current_C",
        "ppg12_D",
        "current_D",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row[field] for field in fields})


def make_plot(rows: list[dict[str, float]], path: Path) -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "axes.linewidth": 1.2,
            "axes.labelsize": 13,
            "xtick.labelsize": 11,
            "ytick.labelsize": 11,
            "legend.fontsize": 10,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )
    x = [row["pt_center"] for row in rows]
    fig, axes = plt.subplots(2, 1, figsize=(9.2, 8.2), dpi=210, sharex=True)
    raw_specs = [
        ("A", "current_over_ppg12_A", "#0072B2", "o"),
        ("B", "current_over_ppg12_B", "#D55E00", "s"),
        ("C", "current_over_ppg12_C", "#009E73", "^"),
        ("D", "current_over_ppg12_D", "#CC79A7", "v"),
    ]
    for label, key, color, marker in raw_specs:
        axes[0].plot(x, [row[key] for row in rows], marker + "-", color=color, lw=1.8, ms=5.8, label=label)
    axes[0].axhline(1.0, color="0.25", lw=1.1)
    axes[0].set_ylabel("Raw stage: current / PPG12")
    axes[0].set_title("Signal-leakage stage parity against PPG12 MC_efficiency_bdt_nom", fontsize=14)
    axes[0].grid(True, color="0.82", alpha=0.4)
    axes[0].legend(frameon=False, ncol=4, loc="upper left")

    derived_specs = [
        ("B/A", "current_over_ppg12_B_over_A", "#D55E00", "s"),
        ("C/A", "current_over_ppg12_C_over_A", "#009E73", "^"),
        ("D/A", "current_over_ppg12_D_over_A", "#CC79A7", "v"),
        ("(C+D)/A", "current_over_ppg12_nt_over_A", "#000000", "o"),
    ]
    for label, key, color, marker in derived_specs:
        axes[1].plot(x, [row[key] for row in rows], marker + "-", color=color, lw=1.8, ms=5.8, label=label)
    axes[1].axhline(1.0, color="0.25", lw=1.1)
    axes[1].set_xlabel(r"Cluster $E_T$ [GeV]")
    axes[1].set_ylabel("Derived leakage ratio parity")
    axes[1].grid(True, color="0.82", alpha=0.4)
    axes[1].legend(frameon=False, ncol=4, loc="upper right")
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def write_summary(rows: list[dict[str, float]], path: Path, current_root: Path, ppg12_stage: Path) -> dict[str, object]:
    stable = [row for row in rows if row["pt_hi"] <= 26.0]
    ranges = {
        "raw_A_current_over_PPG12": finite_range(stable, "current_over_ppg12_A"),
        "raw_B_current_over_PPG12": finite_range(stable, "current_over_ppg12_B"),
        "raw_C_current_over_PPG12": finite_range(stable, "current_over_ppg12_C"),
        "raw_D_current_over_PPG12": finite_range(stable, "current_over_ppg12_D"),
        "B_over_A_current_over_PPG12": finite_range(stable, "current_over_ppg12_B_over_A"),
        "C_over_A_current_over_PPG12": finite_range(stable, "current_over_ppg12_C_over_A"),
        "D_over_A_current_over_PPG12": finite_range(stable, "current_over_ppg12_D_over_A"),
        "non_tight_over_A_current_over_PPG12": finite_range(stable, "current_over_ppg12_nt_over_A"),
        "D_share_nt_current_over_PPG12": finite_range(stable, "current_over_ppg12_D_share_nt"),
    }
    interpretation = [
        "B/A is aligned because raw B and raw A are inflated by nearly the same factor relative to PPG12.",
        "C/A and D/A are low mostly because the current tight-isolated signal denominator A is high relative to the non-tight signal population.",
        "The remaining target is tight/non-tight BDT-sideband parity, not another isolation-topology variation.",
        "The next code-level suspect is exact PPG12 BDT-score branch and sideband-selection equivalence in RecoilJets/PhotonClusterBuilder.",
    ]
    summary = {
        "current_root": str(current_root),
        "ppg12_stage_csv": str(ppg12_stage),
        "stable_bin_definition": "pt_hi <= 26 GeV",
        "stable_ranges": ranges,
        "interpretation": interpretation,
    }
    path.with_suffix(".json").write_text(json.dumps(summary, indent=2) + "\n")
    with path.open("w") as handle:
        handle.write("# PPG12 leakage stage-parity diagnostic\n\n")
        handle.write(f"- Current ROOT: `{current_root}`\n")
        handle.write(f"- PPG12 stage extract: `{ppg12_stage}`\n")
        handle.write("- PPG12 source ROOT: `/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_bdt_nom.root`\n\n")
        handle.write("## Stable-bin ranges\n\n")
        for key, (lo, hi) in ranges.items():
            handle.write(f"- `{key}`: `{lo:.3f}` to `{hi:.3f}`\n")
        handle.write("\n## Interpretation\n\n")
        for item in interpretation:
            handle.write(f"- {item}\n")
        handle.write("\n## Stage table\n\n")
        handle.write("| ET bin | raw A | raw B | raw C | raw D | B/A | C/A | D/A | (C+D)/A |\n")
        handle.write("| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |\n")
        for row in rows:
            handle.write(
                f"| {row['pt_lo']:.0f}-{row['pt_hi']:.0f} | "
                f"{row['current_over_ppg12_A']:.3f} | "
                f"{row['current_over_ppg12_B']:.3f} | "
                f"{row['current_over_ppg12_C']:.3f} | "
                f"{row['current_over_ppg12_D']:.3f} | "
                f"{row['current_over_ppg12_B_over_A']:.3f} | "
                f"{row['current_over_ppg12_C_over_A']:.3f} | "
                f"{row['current_over_ppg12_D_over_A']:.3f} | "
                f"{row['current_over_ppg12_nt_over_A']:.3f} |\n"
            )
    return summary


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--current-root", type=Path, default=DEFAULT_CURRENT_ROOT)
    parser.add_argument("--ppg12-stage-csv", type=Path, default=DEFAULT_PPG12_STAGE)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    ppg12 = read_ppg12_stage(args.ppg12_stage_csv)
    current = read_current_root(args.current_root)
    rows = build_rows(ppg12, current)
    stem = "current_vs_ppg12_fig29_leakage_stage_parity"
    csv_path = args.outdir / f"{stem}.csv"
    md_path = args.outdir / f"{stem}.md"
    png_path = args.outdir / f"{stem}.png"
    write_csv(rows, csv_path)
    make_plot(rows, png_path)
    write_summary(rows, md_path, args.current_root, args.ppg12_stage_csv)
    print(csv_path)
    print(md_path)
    print(md_path.with_suffix(".json"))
    print(png_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
