#!/usr/bin/env python3
"""Compose the current PPG12 Fig. 3 isolation template and data-shape ratio.

The upper panel uses the completed THE-97 p+p data and matched photon+jet
simulation under the PPG12 isolation-template normalization contract.  The
lower panel uses the separately audited, area-normalized tight-data comparison
against the PPG12 SDCC ROOT reference.  Inputs are the CSV/manifest products
of the two source generators; this script verifies their provenance before
composing the reader-facing PNG.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import numpy as np


REPO = Path(__file__).resolve().parents[3]
DEFAULT_OUTDIR = (
    REPO
    / "dataOutput/ppg12Parity"
    / "the97_ppg12_final_accepted_triple_full_20260714_1550"
    / "final_pp_data_canonical_20260717/isolation_consistency"
)
DEFAULT_TEMPLATE_CSV = DEFAULT_OUTDIR / "FULL_STAT_CURRENT.csv"
DEFAULT_TEMPLATE_MANIFEST = DEFAULT_OUTDIR / "FULL_STAT_CURRENT_manifest.json"
DEFAULT_RATIO_CSV = (
    DEFAULT_OUTDIR
    / "ppg12_fig3_iso_template_pt16_22_sdcc_vs_current_THE97_FULL_STAT_CURRENT_ratio.csv"
)
DEFAULT_RATIO_MANIFEST = (
    DEFAULT_OUTDIR
    / "ppg12_fig3_iso_template_pt16_22_sdcc_vs_current_THE97_FULL_STAT_CURRENT_ratio_manifest.json"
)
DEFAULT_DATA_POINTER = REPO / "dataOutput/current_recoiljets_artifacts/current/pp_data_merged/current.json"
DEFAULT_MC_POINTER = REPO / "dataOutput/current_recoiljets_artifacts/current/pp_sim_photonjet_merged/current.json"
DEFAULT_STEM = "ppg12_fig3_isolation_template_with_data_ratio_FULL_STAT_CURRENT"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def load_csv(path: Path) -> dict[str, np.ndarray]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise RuntimeError(f"empty source CSV: {path}")
    return {
        key: np.asarray([float(row[key]) for row in rows], dtype=float)
        for key in rows[0]
    }


def require_same_binning(template: dict[str, np.ndarray], ratio: dict[str, np.ndarray]) -> None:
    for key in ("bin_low", "bin_high", "bin_center"):
        if template[key].shape != ratio[key].shape or not np.allclose(
            template[key], ratio[key], rtol=0.0, atol=1.0e-9
        ):
            raise ValueError(f"template and ratio CSVs have incompatible {key}")


def validate_provenance(
    template_manifest: dict[str, Any],
    ratio_manifest: dict[str, Any],
    data_pointer: dict[str, Any],
    mc_pointer: dict[str, Any],
) -> None:
    coverage = template_manifest.get("coverage", {})
    if coverage.get("completed_files") != coverage.get("expected_files"):
        raise RuntimeError("template source is not full-stat")
    if template_manifest.get("status") != "FULL_STAT_CURRENT":
        raise RuntimeError("template source is not marked FULL_STAT_CURRENT")

    data_roots = [str(Path(path).resolve()) for path in template_manifest["inputs"]["data_roots"]]
    current_data_roots = [str(Path(path).resolve()) for path in data_pointer["root_paths"]]
    if data_roots != current_data_roots:
        raise RuntimeError("template data ROOT does not match the current pp-data pointer")

    mc_root = str(Path(template_manifest["inputs"]["mc_root"]).resolve())
    current_mc_roots = [str(Path(path).resolve()) for path in mc_pointer["root_paths"]]
    if mc_root not in current_mc_roots:
        raise RuntimeError("template MC ROOT does not match the current photon+jet pointer")

    if ratio_manifest["current"]["source_roots"] != data_roots:
        raise RuntimeError("ratio and template do not use the same current data ROOT")
    if ratio_manifest["processing_contract"]["ratio"] != "current density / PPG12 SDCC density":
        raise RuntimeError("unexpected ratio orientation")
    if not ratio_manifest["processing_contract"]["normalization"].startswith("current area scaled"):
        raise RuntimeError("tight-data ratio is not area-normalized to PPG12")

    expected_data_sha = template_manifest["inputs"]["data_sha256"][data_roots[0]]
    if sha256_file(Path(data_roots[0])) != expected_data_sha:
        raise RuntimeError("current pp-data ROOT SHA-256 changed")
    if sha256_file(Path(mc_root)) != template_manifest["inputs"]["mc_sha256"]:
        raise RuntimeError("current photon+jet ROOT SHA-256 changed")


def step_values(edges: np.ndarray, values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    return edges, np.r_[values, values[-1]]


def draw(
    output: Path,
    template: dict[str, np.ndarray],
    ratio: dict[str, np.ndarray],
    ratio_manifest: dict[str, Any],
) -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "mathtext.fontset": "dejavusans",
            "axes.linewidth": 1.35,
            "axes.labelsize": 18,
            "xtick.labelsize": 14,
            "ytick.labelsize": 14,
            "legend.fontsize": 13,
        }
    )

    edges = np.r_[template["bin_low"], template["bin_high"][-1]]
    centers = template["bin_center"]
    tight = template["data_signal_selection_counts_per_width"]
    tight_err = template["data_signal_selection_error"]
    background = template["data_background_counts_per_width"]
    signal_mc = template["signal_mc_counts_per_width"]
    ppg12_tight = ratio["ppg12_sdcc_counts_per_width"]
    ppg12_tight_err = ratio["ppg12_sdcc_error"]
    stack_total = background + signal_mc

    fig = plt.figure(figsize=(8.0, 7.2), dpi=200)
    grid = fig.add_gridspec(2, 1, height_ratios=(3.25, 1.0), hspace=0.035)
    ax = fig.add_subplot(grid[0])
    rax = fig.add_subplot(grid[1], sharex=ax)

    x_step, background_step = step_values(edges, background)
    _, stack_step = step_values(edges, stack_total)
    ax.fill_between(
        x_step,
        0.0,
        background_step,
        step="post",
        color="#f2a3a0",
        edgecolor="#c93b32",
        linewidth=1.25,
        alpha=0.72,
        label="Non-tight data background template",
        zorder=1,
    )
    ax.fill_between(
        x_step,
        background_step,
        stack_step,
        step="post",
        color="#aeb8ff",
        edgecolor="#4256c9",
        linewidth=1.25,
        alpha=0.72,
        label="Photon+jet signal simulation",
        zorder=2,
    )
    ax.errorbar(
        centers,
        tight,
        yerr=tight_err,
        fmt="o",
        markersize=3.5,
        markerfacecolor="black",
        markeredgecolor="black",
        color="black",
        ecolor="black",
        elinewidth=0.9,
        capsize=0,
        label="Tight photon candidates in data",
        zorder=5,
    )
    ax.errorbar(
        centers,
        ppg12_tight,
        yerr=ppg12_tight_err,
        fmt="o",
        markersize=6.0,
        markerfacecolor="none",
        markeredgecolor="#c54b7d",
        markeredgewidth=1.25,
        color="#c54b7d",
        ecolor="#d887a5",
        elinewidth=0.85,
        capsize=0,
        label="PPG12 SDCC tight photon candidates",
        zorder=6,
    )

    visible = (centers >= -1.0) & (centers <= 15.0)
    ymax = max(
        np.nanmax(tight[visible] + tight_err[visible]),
        np.nanmax(ppg12_tight[visible] + ppg12_tight_err[visible]),
        np.nanmax(stack_total[visible]),
    )
    ax.set_xlim(-1.0, 15.0)
    ax.set_ylim(0.0, 1.16 * ymax)
    ax.set_ylabel("Counts / bin width")
    ax.tick_params(which="both", direction="in", top=True, right=True, labelbottom=False)
    ax.tick_params(which="major", length=6)
    ax.tick_params(which="minor", length=3)
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))

    ax.text(
        0.965,
        0.955,
        r"$\it{\bf{sPHENIX}}$ Internal",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=17,
    )
    ax.text(
        0.965,
        0.895,
        r"$p{+}p\ \sqrt{s}=200\ \mathrm{GeV}$",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=15,
    )
    ax.text(
        0.965,
        0.838,
        r"$16<E_T^\gamma<22\ \mathrm{GeV},\quad |\eta^\gamma|<0.7$",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=14,
    )
    handles, labels = ax.get_legend_handles_labels()
    order = [3, 2, 0, 1]
    ax.legend(
        [handles[index] for index in order],
        [labels[index] for index in order],
        loc="upper right",
        bbox_to_anchor=(0.97, 0.78),
        frameon=False,
        handlelength=1.8,
        borderaxespad=0.0,
    )

    ratio_values = ratio["current_over_ppg12"]
    ratio_errors = ratio["current_over_ppg12_error"]
    ratio_mask = visible & np.isfinite(ratio_values) & np.isfinite(ratio_errors)
    rax.axhline(1.0, color="#666666", linewidth=1.0, linestyle=(0, (4, 4)))
    rax.errorbar(
        centers[ratio_mask],
        ratio_values[ratio_mask],
        yerr=ratio_errors[ratio_mask],
        fmt="o",
        markersize=3.6,
        markerfacecolor="#2878b5",
        markeredgecolor="#2878b5",
        color="#2878b5",
        ecolor="#5a9bd0",
        elinewidth=0.9,
        capsize=0,
    )
    positive_upper = ratio_values[ratio_mask] + ratio_errors[ratio_mask]
    ratio_ymax = max(1.2, float(np.nanmax(positive_upper)))
    ratio_ymax = float(np.ceil(ratio_ymax * 5.0) / 5.0)
    rax.set_ylim(0.0, ratio_ymax)
    rax.set_ylabel("This analysis / PPG12", fontsize=14)
    rax.set_xlabel(r"$E_T^{\mathrm{iso,reco}}\ \mathrm{[GeV]}$", loc="right")
    rax.tick_params(which="both", direction="in", top=True, right=True)
    rax.tick_params(which="major", length=6)
    rax.tick_params(which="minor", length=3)
    rax.xaxis.set_minor_locator(AutoMinorLocator(5))
    rax.yaxis.set_minor_locator(AutoMinorLocator(4))
    fig.subplots_adjust(left=0.13, right=0.975, bottom=0.105, top=0.975)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, facecolor="white")
    plt.close(fig)


def write_combined_csv(
    output: Path, template: dict[str, np.ndarray], ratio: dict[str, np.ndarray]
) -> None:
    columns = [
        "bin_low",
        "bin_high",
        "bin_center",
        "data_signal_selection_counts_per_width",
        "data_signal_selection_error",
        "data_background_counts_per_width",
        "data_background_error",
        "signal_mc_counts_per_width",
        "signal_mc_error",
        "stack_total_counts_per_width",
        "ppg12_sdcc_tight_data_counts_per_width",
        "ppg12_sdcc_tight_data_error",
        "current_tight_data_area_normalized_counts_per_width",
        "current_tight_data_area_normalized_error",
        "current_over_ppg12",
        "current_over_ppg12_error",
    ]
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(columns)
        for index in range(len(template["bin_center"])):
            writer.writerow(
                [
                    template["bin_low"][index],
                    template["bin_high"][index],
                    template["bin_center"][index],
                    template["data_signal_selection_counts_per_width"][index],
                    template["data_signal_selection_error"][index],
                    template["data_background_counts_per_width"][index],
                    template["data_background_error"][index],
                    template["signal_mc_counts_per_width"][index],
                    template["signal_mc_error"][index],
                    template["stack_total_counts_per_width"][index],
                    ratio["ppg12_sdcc_counts_per_width"][index],
                    ratio["ppg12_sdcc_error"][index],
                    ratio["current_counts_per_width"][index],
                    ratio["current_error"][index],
                    ratio["current_over_ppg12"][index],
                    ratio["current_over_ppg12_error"][index],
                ]
            )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--template-csv", type=Path, default=DEFAULT_TEMPLATE_CSV)
    parser.add_argument("--template-manifest", type=Path, default=DEFAULT_TEMPLATE_MANIFEST)
    parser.add_argument("--ratio-csv", type=Path, default=DEFAULT_RATIO_CSV)
    parser.add_argument("--ratio-manifest", type=Path, default=DEFAULT_RATIO_MANIFEST)
    parser.add_argument("--data-pointer", type=Path, default=DEFAULT_DATA_POINTER)
    parser.add_argument("--mc-pointer", type=Path, default=DEFAULT_MC_POINTER)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--stem", default=DEFAULT_STEM)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    template_manifest = load_json(args.template_manifest)
    ratio_manifest = load_json(args.ratio_manifest)
    data_pointer = load_json(args.data_pointer)
    mc_pointer = load_json(args.mc_pointer)
    validate_provenance(template_manifest, ratio_manifest, data_pointer, mc_pointer)

    template = load_csv(args.template_csv)
    ratio = load_csv(args.ratio_csv)
    require_same_binning(template, ratio)

    png_path = args.outdir / f"{args.stem}.png"
    csv_path = args.outdir / f"{args.stem}.csv"
    manifest_path = args.outdir / f"{args.stem}_manifest.json"
    draw(png_path, template, ratio, ratio_manifest)
    write_combined_csv(csv_path, template, ratio)

    manifest = {
        "artifact": "PPG12 Figure 3 isolation template with current-data shape ratio",
        "status": "FULL_STAT_CURRENT",
        "campaign": template_manifest["campaign"],
        "selection": template_manifest["selection"],
        "upper_panel": {
            "definition": "tight data over stacked non-tight data background template plus photon+jet signal simulation",
            "normalization": template_manifest["ppg12_contract"],
        },
        "lower_panel": {
            "definition": "current tight-data isolation density divided by PPG12 SDCC tight-data isolation density",
            "normalization_scale": ratio_manifest["current"]["normalization_scale"],
            "metrics": ratio_manifest["metrics"],
        },
        "inputs": {
            "template_csv": str(args.template_csv.resolve()),
            "template_csv_sha256": sha256_file(args.template_csv.resolve()),
            "template_manifest": str(args.template_manifest.resolve()),
            "template_manifest_sha256": sha256_file(args.template_manifest.resolve()),
            "ratio_csv": str(args.ratio_csv.resolve()),
            "ratio_csv_sha256": sha256_file(args.ratio_csv.resolve()),
            "ratio_manifest": str(args.ratio_manifest.resolve()),
            "ratio_manifest_sha256": sha256_file(args.ratio_manifest.resolve()),
            "data_pointer": str(args.data_pointer.resolve()),
            "data_pointer_sha256": sha256_file(args.data_pointer.resolve()),
            "mc_pointer": str(args.mc_pointer.resolve()),
            "mc_pointer_sha256": sha256_file(args.mc_pointer.resolve()),
        },
        "outputs": {
            "png": str(png_path.resolve()),
            "png_sha256": sha256_file(png_path),
            "csv": str(csv_path.resolve()),
            "csv_sha256": sha256_file(csv_path),
            "manifest": str(manifest_path.resolve()),
        },
        "limitations": [
            "The upper red component is a non-tight data background template, not background simulation.",
            "The lower panel is an area-normalized shape comparison and does not test absolute selected yield.",
        ],
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(manifest["outputs"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
