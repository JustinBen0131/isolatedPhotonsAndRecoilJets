#!/usr/bin/env python3
"""Render PPG12 IAN Fig. 11 directly from an SDCC read-only projection CSV."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
DEFAULT_OUTDIR = REPO / "dataOutput/ppg12Parity/ppg12_fig11_sdcc_reference"
DEFAULT_CSV = DEFAULT_OUTDIR / "fig11_sb_sdcc_verbatim_points.csv"
REMOTE_SIGNAL = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiencyshower_shape_signal_combined_showershape.root"
REMOTE_BACKGROUND = "/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiencyshower_shape_jet_inclusive_combined_showershape.root"
REMOTE_SIGNAL_SHA256 = "cd5b04a20b5588519ae9b65a3f1cc80c2fd74739b747ab806793320cb9fdbce8"
REMOTE_BACKGROUND_SHA256 = "b5115e793538c717e4efed31252e2db8def19520e84435c57241a5d6690bba84"
REMOTE_MACRO_SHA256 = "da21cd84af105385eebf9dc66fe77101a2f9e67039d9d8313deba39621bad8c5"


def read_points(path: Path) -> list[dict[str, float]]:
    with path.open() as handle:
        rows = [{key: float(value) for key, value in row.items()} for row in csv.DictReader(handle)]
    if len(rows) != 10:
        raise RuntimeError(f"Expected the ten nonzero Fig.11 points, found {len(rows)} in {path}")
    if any(row["s_over_b"] <= 0.0 or row["s_over_b_err"] < 0.0 for row in rows):
        raise RuntimeError("Fig.11 S/B inputs must be finite nonnegative values")
    return rows


def render(rows: list[dict[str, float]], output: Path) -> dict[str, float]:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "dejavuserif",
            "axes.linewidth": 1.25,
            "xtick.direction": "in",
            "ytick.direction": "in",
        }
    )
    fig, ax = plt.subplots(figsize=(6.0, 6.0), dpi=220)
    x = np.asarray([row["xcenter"] for row in rows])
    y = np.asarray([row["s_over_b"] for row in rows])
    yerr = np.asarray([row["s_over_b_err"] for row in rows])
    ax.errorbar(x, y, yerr=yerr, fmt="o", color="black", ms=5.4, elinewidth=1.0, capsize=0, zorder=3)
    ax.set_xlim(10.0, 32.0)
    ax.set_ylim(0.0, 1.0)
    ax.set_xticks(np.arange(10.0, 33.0, 2.0))
    ax.set_yticks(np.arange(0.0, 1.01, 0.1))
    ax.set_xlabel(r"$E_{\mathrm{T}}^{\gamma,\mathrm{rec}}$ [GeV]", fontsize=16, ha="right", x=1.0)
    ax.set_ylabel("S/B", fontsize=16)
    ax.minorticks_on()
    ax.tick_params(which="both", top=True, right=True, labelsize=12)
    ax.text(0.43, 0.93, "sPHENIX", transform=ax.transAxes, fontsize=13, fontstyle="italic", fontweight="bold")
    ax.text(0.66, 0.93, "Internal", transform=ax.transAxes, fontsize=13)
    ax.text(0.43, 0.855, r"$p$+$p$ $\sqrt{s}=200$ GeV", transform=ax.transAxes, fontsize=12)
    ax.text(0.43, 0.79, "PYTHIA8", transform=ax.transAxes, fontsize=12)
    fig.subplots_adjust(left=0.14, right=0.97, top=0.97, bottom=0.14)
    fig.savefig(output)
    plt.close(fig)
    return {"point_count": len(rows), "min_s_over_b": float(y.min()), "max_s_over_b": float(y.max())}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", type=Path, default=DEFAULT_CSV)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUTDIR)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    output = args.out_dir / "ppg12_ian_fig11_sb_sdcc_verbatim_reproduction.png"
    manifest = args.out_dir / "ppg12_ian_fig11_sb_sdcc_verbatim_reproduction_manifest.json"
    rows = read_points(args.csv)
    summary = render(rows, output)
    manifest.write_text(
        json.dumps(
            {
                "artifact": str(output),
                "source_csv": str(args.csv),
                "source_csv_sha256": hashlib.sha256(args.csv.read_bytes()).hexdigest(),
                "sdcc_source_macro": "/sphenix/user/shuhangli/ppg12/plotting/plot_SB.C",
                "sdcc_source_macro_sha256": REMOTE_MACRO_SHA256,
                "sdcc_signal_root": REMOTE_SIGNAL,
                "sdcc_signal_root_sha256": REMOTE_SIGNAL_SHA256,
                "sdcc_background_root": REMOTE_BACKGROUND,
                "sdcc_background_root_sha256": REMOTE_BACKGROUND_SHA256,
                "source_histogram": "h_ET_isoET_eta0 (TH2D in both roots)",
                "construction": "PPG12 plot_SB.C verbatim: RebinX(16); project -1 < isoET < 0.502095 + 0.0433036*ET; scale background by jet50cross/photon20cross = 7.3113/130.4461; TH1 Divide error propagation.",
                "no_digitization_or_fit": True,
                "summary": summary,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    print(output)
    print(manifest)


if __name__ == "__main__":
    main()
