#!/usr/bin/env python3
from __future__ import annotations

import csv
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-thesisanalysis")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm


REPORT = Path("dataOutput/auauTightBDTValidation/model_validation_condor_20260509_192942")
OUTDIR = Path("dataOutput/jstg_slide_candidates/slide17_ptcent_auc_gain")
TABLE = REPORT / "validation_auc_table.csv"

BASELINE = "centAsFeat_pt5to40"
ROWS = [
    ("ptCentDep3", r"$E_T \times$ 3 centrality-bin BDTs"),
    ("ptCentDepFine", r"$E_T \times$ 7 centrality-bin BDTs"),
]
CENTRALITIES = [("0_20", "0-20%"), ("20_50", "20-50%"), ("50_80", "50-80%")]
PT_BINS = [("6_10", "6-10"), ("10_15", "10-15"), ("15_20", "15-20"), ("20_25", "20-25"), ("25_35", "25-35")]


def read_auc_table() -> dict[tuple[str, str, str], dict[str, float]]:
    out: dict[tuple[str, str, str], dict[str, float]] = {}
    with TABLE.open() as handle:
        for row in csv.DictReader(handle):
            auc = np.nan if row["auc"].lower() == "nan" else float(row["auc"])
            out[(row["product"], row["centrality_bin"], row["pt_bin"])] = {
                "auc": auc,
                "entries": float(row["entries"]),
                "signal": float(row["signal_entries"]),
                "background": float(row["background_entries"]),
            }
    return out


def build_matrix(rows: dict[tuple[str, str, str], dict[str, float]]) -> np.ndarray:
    vals = np.full((len(ROWS), len(CENTRALITIES) * len(PT_BINS)), np.nan)
    for iy, (product, _) in enumerate(ROWS):
        ix = 0
        for cent, _ in CENTRALITIES:
            for pt, _ in PT_BINS:
                base = rows[(BASELINE, cent, pt)]["auc"]
                val = rows[(product, cent, pt)]["auc"]
                vals[iy, ix] = 100.0 * (val - base) / base
                ix += 1
    return vals


def sphenix_label(fig: plt.Figure) -> None:
    fig.text(0.205, 0.855, "sPHENIX", fontsize=18, fontstyle="italic", fontweight="bold", ha="left", va="baseline")
    fig.text(0.305, 0.855, "Internal", fontsize=18, ha="left", va="baseline")
    fig.text(0.205, 0.818, r"Pythia overlay, $\sqrt{s_{NN}} = 200$ GeV", fontsize=13, ha="left", va="baseline")


def write_csv(values: np.ndarray, auc_rows: dict[tuple[str, str, str], dict[str, float]]) -> None:
    path = OUTDIR / "ptcent_auc_gain_vs_centinput_table.csv"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "model",
                "centrality",
                "pt_bin",
                "baseline_auc",
                "dedicated_auc",
                "auc_gain_percent",
                "entries",
                "signal_entries",
                "background_entries",
            ]
        )
        for product, _ in ROWS:
            for cent, cent_label in CENTRALITIES:
                for pt, pt_label in PT_BINS:
                    base = auc_rows[(BASELINE, cent, pt)]
                    val = auc_rows[(product, cent, pt)]
                    gain = 100.0 * (val["auc"] - base["auc"]) / base["auc"]
                    writer.writerow(
                        [
                            product,
                            cent_label,
                            pt_label,
                            f"{base['auc']:.8f}",
                            f"{val['auc']:.8f}",
                            f"{gain:.5f}",
                            int(val["entries"]),
                            int(val["signal"]),
                            int(val["background"]),
                        ]
                    )


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    auc_rows = read_auc_table()
    values = build_matrix(auc_rows)
    write_csv(values, auc_rows)

    colors = [
        (0.18, 0.44, 0.70),
        (0.75, 0.84, 0.91),
        (0.96, 0.94, 0.86),
        (0.98, 0.70, 0.38),
        (0.88, 0.20, 0.20),
        (0.55, 0.00, 0.05),
    ]
    cmap = LinearSegmentedColormap.from_list("gain_map", colors, N=256)
    norm = TwoSlopeNorm(vmin=-1.0, vcenter=0.0, vmax=22.0)

    fig, ax = plt.subplots(figsize=(16.0, 7.2), dpi=180)
    fig.subplots_adjust(left=0.205, right=0.865, bottom=0.18, top=0.620)
    image = ax.imshow(values, cmap=cmap, norm=norm, aspect="auto")

    for iy in range(values.shape[0]):
        for ix in range(values.shape[1]):
            val = values[iy, ix]
            color = "white" if val >= 10.0 else "black"
            ax.text(ix, iy, f"{val:+.1f}%", ha="center", va="center", fontsize=12.5, color=color)

    ax.set_yticks(np.arange(len(ROWS)))
    ax.set_yticklabels([label for _, label in ROWS], fontsize=14)
    ax.set_xticks(np.arange(values.shape[1]))
    ax.set_xticklabels([pt_label for _cent, _cent_label in CENTRALITIES for _pt, pt_label in PT_BINS], fontsize=12)
    ax.set_xlabel(r"Photon candidate $E_T$ bin within each centrality group [GeV]", fontsize=17, labelpad=16)
    fig.suptitle(r"Percent gain in ROC-AUC from dedicated $E_T \times$ centrality BDTs", fontsize=22, y=0.965)

    ax.set_xticks(np.arange(-0.5, values.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-0.5, values.shape[0], 1), minor=True)
    ax.grid(which="minor", color="0.28", linewidth=0.75)
    ax.tick_params(which="both", length=0)

    cells_per_cent = len(PT_BINS)
    for boundary in (cells_per_cent - 0.5, 2 * cells_per_cent - 0.5):
        ax.axvline(boundary, color="0.35", linestyle="--", linewidth=2.0)

    for idx, (_cent, label) in enumerate(CENTRALITIES):
        center = idx * cells_per_cent + (cells_per_cent - 1) / 2
        ax.text(center, -0.83, label, ha="center", va="center", fontsize=17)
        x0 = idx * cells_per_cent - 0.42
        x1 = (idx + 1) * cells_per_cent - 0.58
        ax.plot([x0, x0, x1, x1], [-0.66, -0.74, -0.74, -0.66], color="0.25", linewidth=1.5, clip_on=False)

    for spine in ax.spines.values():
        spine.set_linewidth(1.2)
        spine.set_color("0.25")

    cax = fig.add_axes([0.888, 0.18, 0.020, 0.440])
    cb = fig.colorbar(image, cax=cax)
    cb.set_label("Percent gain wrt centrality-as-input BDT [%]", fontsize=16, labelpad=18)
    cb.ax.tick_params(labelsize=13)

    sphenix_label(fig)

    out = OUTDIR / "ptcent_auc_gain_vs_centinput_heatmap.png"
    fig.savefig(out)
    plt.close(fig)
    print(out)


if __name__ == "__main__":
    main()
