#!/usr/bin/env python3
"""Render the THE-240 photon-definition validation figure from offline artifacts.

The left panel reuses the frozen full-stat pp PPG12 purity-closure point table.
The Au+Au panels read the current registered THE-88 data, embedded-photon, and
embedded-inclusive ROOT files.  They integrate 15 < E_T^gamma < 35 GeV in
0--20% centrality and compare E11/E33 before and after tight photon ID.

No production, unfolding, calibration, systematic band, or Google Slides
mutation is performed by this script.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import uproot
from matplotlib.lines import Line2D


REPO = Path(__file__).resolve().parents[3]
CURRENT = REPO / "dataOutput/current_recoiljets_artifacts/current"
DEFAULT_OUT = REPO / "dataOutput/the219_friday_ppg_20260814/the240_photon_definition_validation"
PURITY_CSV = (
    REPO
    / "dataOutput/ppg12Parity/the97_ppg12_final_accepted_triple_full_20260714_1550"
    / "final_pp_data_canonical_20260717/purity_overlay_current"
    / "the97_pp_data_purity_ppg12_vs_current_fullstat_points.csv"
)
PURITY_MANIFEST = PURITY_CSV.with_name(
    "the97_pp_data_purity_ppg12_vs_current_fullstat_manifest.json"
)

PT_BINS = ("15_17", "17_19", "19_21", "21_23", "23_26", "26_35")
DATA_DIRECTORY = "photon_10_plus_MBD_NS_geq_2_vtx_lt_150"
CENTRALITY = "0_20"
DISPLAY_REBIN = 4
FILL_MULTIPLICITY = 4

INK = "#111827"
MUTED = "#526173"
GRID = "#D8E0E9"
PPG12_BLUE = "#2468A2"
PHOTON_RED = "#C9362B"
JET_BLUE = "#2864C7"


@dataclass(frozen=True)
class Curve:
    x: np.ndarray
    y: np.ndarray
    error: np.ndarray
    bin_width: float
    source_objects: tuple[str, ...]
    visible_sum: float
    exact_zero_sum: float
    effective_entries: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument(
        "--visual-qa-pass",
        action="store_true",
        help="Record PASS only after the rendered 2560x1440 PNG has been inspected.",
    )
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def resolve_pointer(sample_key: str) -> tuple[Path, Path, dict[str, Any]]:
    pointer_path = CURRENT / sample_key / "current.json"
    payload = json.loads(pointer_path.read_text(encoding="utf-8"))
    roots = [Path(item) for item in payload.get("root_paths", [])]
    if len(roots) != 1 or not roots[0].is_file():
        raise FileNotFoundError(f"Expected one readable ROOT path in {pointer_path}: {roots}")
    return roots[0], pointer_path, payload


def read_purity() -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    with PURITY_CSV.open(newline="", encoding="utf-8") as stream:
        for raw in csv.DictReader(stream):
            rows.append({key: float(value) for key, value in raw.items()})
    if not rows:
        raise RuntimeError(f"No purity rows found in {PURITY_CSV}")
    return rows


def histogram_name(stage: str, lane: str, pt_bin: str) -> str:
    if lane == "data":
        return f"{DATA_DIRECTORY}/h_ss_e11e33_{stage}_pT_{pt_bin}_cent_{CENTRALITY}"
    suffix = "sig" if lane == "photon" else "bkg"
    return f"SIM/h_ss_e11e33_{stage}_{suffix}_pT_{pt_bin}_cent_{CENTRALITY}"


def load_curve(root_path: Path, stage: str, lane: str) -> Curve:
    values: np.ndarray | None = None
    variances: np.ndarray | None = None
    edges: np.ndarray | None = None
    entries = 0.0
    names = tuple(histogram_name(stage, lane, pt_bin) for pt_bin in PT_BINS)
    with uproot.open(root_path) as root_file:
        for name in names:
            if name not in root_file:
                raise KeyError(f"Missing expected histogram {name} in {root_path}")
            hist = root_file[name]
            current_values = np.asarray(hist.values(flow=True), dtype=float)
            current_variances = hist.variances(flow=True)
            if current_variances is None:
                current_variances = np.abs(current_values)
            else:
                current_variances = np.asarray(current_variances, dtype=float)
            values = current_values.copy() if values is None else values + current_values
            variances = (
                current_variances.copy()
                if variances is None
                else variances + current_variances
            )
            current_edges = np.asarray(hist.axis().edges(), dtype=float)
            if edges is None:
                edges = current_edges
            elif not np.array_equal(edges, current_edges):
                raise RuntimeError(f"Inconsistent binning at {name}")
            entries += float(hist.member("fEntries"))

    assert values is not None and variances is not None and edges is not None
    # Flow layout is [underflow, ordinary bins..., overflow]. Bin 1 is the
    # exact/near-zero boundary component used by the established THE-88 display
    # contract; retain it in provenance but omit it from the continuous shape.
    exact_zero = float(values[1])
    visible_values = np.asarray(values[2:-1], dtype=float)
    visible_variances = np.asarray(variances[2:-1], dtype=float)
    centers = 0.5 * (edges[:-1] + edges[1:])
    visible_centers = np.asarray(centers[1:], dtype=float)
    if np.any(visible_values < 0):
        raise RuntimeError(f"Negative visible weights in {lane} {stage}")
    visible_sum = float(np.sum(visible_values))
    if visible_sum <= 0:
        raise RuntimeError(f"Nonpositive visible sum in {lane} {stage}")

    groups = [
        slice(start, min(start + DISPLAY_REBIN, len(visible_values)))
        for start in range(0, len(visible_values), DISPLAY_REBIN)
    ]
    rebinned_x = np.asarray([np.mean(visible_centers[group]) for group in groups])
    rebinned_y = np.asarray([np.sum(visible_values[group]) for group in groups])
    rebinned_variance = np.asarray(
        [np.sum(visible_variances[group]) for group in groups]
    )
    # These persisted pre/tight families contain four identical internal-
    # isolation-view fills. Multiplying the normalized errors by sqrt(4)
    # restores the unique-candidate statistical scale.
    rebinned_error = np.sqrt(rebinned_variance) * np.sqrt(FILL_MULTIPLICITY)
    return Curve(
        x=rebinned_x,
        y=rebinned_y / visible_sum,
        error=rebinned_error / visible_sum,
        bin_width=float((edges[1] - edges[0]) * DISPLAY_REBIN),
        source_objects=names,
        visible_sum=visible_sum,
        exact_zero_sum=exact_zero,
        effective_entries=entries / FILL_MULTIPLICITY,
    )


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "mathtext.fontset": "dejavusans",
            "axes.linewidth": 1.15,
            "axes.labelcolor": INK,
            "text.color": INK,
            "xtick.color": INK,
            "ytick.color": INK,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "savefig.facecolor": "white",
        }
    )


def style_axis(ax: plt.Axes) -> None:
    ax.grid(axis="y", color=GRID, lw=0.8, alpha=0.85, zorder=0)
    ax.minorticks_on()
    ax.tick_params(which="major", length=6, width=1.0, labelsize=12)
    ax.tick_params(which="minor", length=3, width=0.75)


def plot_purity(ax: plt.Axes, ratio_ax: plt.Axes, rows: list[dict[str, float]]) -> None:
    x = np.asarray([row["center"] for row in rows])
    xerr = np.asarray([row["half_width"] for row in rows])
    current = np.asarray([row["current_raw"] for row in rows])
    current_error = np.asarray([row["current_raw_error"] for row in rows])
    reference = np.asarray([row["ppg12_raw"] for row in rows])
    reference_error = np.asarray(
        [
            [row["ppg12_raw_error_low"] for row in rows],
            [row["ppg12_raw_error_high"] for row in rows],
        ]
    )
    ratio = np.asarray([row["current_over_ppg12_raw"] for row in rows])
    ratio_error = np.asarray(
        [
            [row["current_over_ppg12_raw_error_low"] for row in rows],
            [row["current_over_ppg12_raw_error_high"] for row in rows],
        ]
    )

    ax.errorbar(
        x,
        reference,
        xerr=xerr,
        yerr=reference_error,
        fmt="s",
        color=PPG12_BLUE,
        markerfacecolor="white",
        markeredgewidth=1.3,
        markersize=6.2,
        elinewidth=1.1,
        capsize=0,
        label="PPG12",
        zorder=3,
    )
    ax.errorbar(
        x,
        current,
        xerr=xerr,
        yerr=current_error,
        fmt="o",
        color=INK,
        markerfacecolor=INK,
        markersize=5.8,
        elinewidth=1.1,
        capsize=0,
        label="Current output",
        zorder=4,
    )
    ax.set_ylabel("Photon purity", fontsize=14)
    ax.set_ylim(0.25, 1.15)
    ax.set_xlim(9.2, 36.8)
    ax.tick_params(labelbottom=False)
    ax.legend(loc="lower right", frameon=False, fontsize=12, handletextpad=0.5)
    ax.set_title("(a) p+p: raw ABCD purity reproducibility", fontsize=16, loc="left", pad=12)
    style_axis(ax)

    ratio_ax.axhline(1.0, color=MUTED, lw=1.0, zorder=1)
    ratio_ax.errorbar(
        x,
        ratio,
        xerr=xerr,
        yerr=ratio_error,
        fmt="o",
        color=INK,
        markerfacecolor=INK,
        markersize=4.8,
        elinewidth=1.0,
        capsize=0,
        zorder=3,
    )
    ratio_ax.set_xlim(9.2, 36.8)
    ratio_ax.set_ylim(0.35, 1.65)
    ratio_ax.set_yticks([0.5, 1.0, 1.5])
    ratio_ax.set_xlabel(r"$E_{T}^{\gamma}$ [GeV]", fontsize=14)
    ratio_ax.set_ylabel("Current / PPG12", fontsize=11)
    style_axis(ratio_ax)


def plot_auau(
    ax: plt.Axes,
    curves: dict[str, Curve],
    *,
    title: str,
    show_y_label: bool,
    y_max: float,
) -> None:
    data = curves["data"]
    photon = curves["photon"]
    inclusive = curves["inclusive"]
    ax.step(
        photon.x,
        photon.y,
        where="mid",
        color=PHOTON_RED,
        lw=2.4,
        zorder=2,
    )
    ax.step(
        inclusive.x,
        inclusive.y,
        where="mid",
        color=JET_BLUE,
        lw=2.4,
        zorder=2,
    )
    ax.errorbar(
        data.x,
        data.y,
        yerr=data.error,
        fmt="o",
        color=INK,
        markerfacecolor=INK,
        markeredgecolor="white",
        markeredgewidth=0.45,
        markersize=5.5,
        elinewidth=1.05,
        capsize=0,
        zorder=4,
    )
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, y_max)
    ax.set_xlabel(r"$E_{1\times1}/E_{3\times3}$", fontsize=15)
    if show_y_label:
        ax.set_ylabel("Unit-normalized fraction / 0.04", fontsize=14)
    else:
        ax.tick_params(labelleft=False)
    ax.set_title(title, fontsize=16, loc="left", pad=12)
    style_axis(ax)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    setup_style()

    roots: dict[str, Path] = {}
    pointers: dict[str, Path] = {}
    pointer_payloads: dict[str, dict[str, Any]] = {}
    for sample_key in (
        "auau_data_merged",
        "auau_sim_photonjet_merged",
        "auau_sim_inclusivejet_merged",
    ):
        root, pointer, payload = resolve_pointer(sample_key)
        roots[sample_key] = root
        pointers[sample_key] = pointer
        pointer_payloads[sample_key] = payload

    purity_rows = read_purity()
    curves = {
        stage: {
            "data": load_curve(roots["auau_data_merged"], stage, "data"),
            "photon": load_curve(roots["auau_sim_photonjet_merged"], stage, "photon"),
            "inclusive": load_curve(roots["auau_sim_inclusivejet_merged"], stage, "inclusive"),
        }
        for stage in ("pre", "tight")
    }
    y_max = 1.18 * max(
        float(np.max(curve.y))
        for stage_curves in curves.values()
        for curve in stage_curves.values()
    )

    fig = plt.figure(figsize=(16, 9), dpi=160)
    grid = fig.add_gridspec(
        2,
        3,
        width_ratios=(1.12, 1.0, 1.0),
        height_ratios=(3.45, 1.08),
        left=0.064,
        right=0.985,
        bottom=0.155,
        top=0.745,
        wspace=0.30,
        hspace=0.05,
    )
    ax_purity = fig.add_subplot(grid[0, 0])
    ax_ratio = fig.add_subplot(grid[1, 0])
    ax_pre = fig.add_subplot(grid[:, 1])
    ax_tight = fig.add_subplot(grid[:, 2], sharey=ax_pre)

    plot_purity(ax_purity, ax_ratio, purity_rows)
    plot_auau(
        ax_pre,
        curves["pre"],
        title="(b) Au+Au 0–20%: before tight ID",
        show_y_label=True,
        y_max=y_max,
    )
    plot_auau(
        ax_tight,
        curves["tight"],
        title="(c) Au+Au 0–20%: after tight ID",
        show_y_label=False,
        y_max=y_max,
    )

    fig.text(
        0.064,
        0.957,
        "sPHENIX",
        fontsize=16,
        fontweight="bold",
        fontstyle="italic",
        ha="left",
        va="top",
    )
    fig.text(0.130, 0.957, "Internal", fontsize=16, ha="left", va="top")
    fig.text(
        0.064,
        0.912,
        "Photon definition: p+p reproducibility and central Au+Au morphology",
        fontsize=24,
        fontweight="bold",
        ha="left",
        va="top",
    )
    fig.text(
        0.064,
        0.866,
        r"Au+Au: $15 < E_T^{\gamma} < 35$ GeV, 0–20% centrality  •  statistical uncertainties only",
        fontsize=14,
        color=MUTED,
        ha="left",
        va="top",
    )
    legend_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color=INK,
            linestyle="none",
            markerfacecolor=INK,
            markersize=7,
            label="Au+Au data (stat.)",
        ),
        Line2D([0], [0], color=PHOTON_RED, lw=2.5, label="Embedded prompt photon"),
        Line2D([0], [0], color=JET_BLUE, lw=2.5, label="Embedded inclusive jet"),
    ]
    fig.legend(
        handles=legend_handles,
        loc="upper right",
        bbox_to_anchor=(0.985, 0.820),
        frameon=False,
        ncol=3,
        fontsize=12.5,
        columnspacing=1.3,
        handlelength=2.2,
    )
    fig.text(
        0.064,
        0.075,
        "Panel (a) shows raw ABCD purity before leakage correction; Au+Au shapes omit the exact-zero boundary bin. No systematic band is shown.",
        fontsize=11.5,
        color=MUTED,
        ha="left",
        va="center",
    )
    fig.text(
        0.064,
        0.045,
        r"$E_{1\times1}/E_{3\times3}$ enters photon ID: panels (b,c) are a morphology reasonableness check, not an independent purity measurement.",
        fontsize=11.5,
        color=MUTED,
        ha="left",
        va="center",
    )

    png = args.output_dir / "the240_photon_definition_validation.png"
    fig.savefig(png, dpi=160, facecolor="white")
    plt.close(fig)

    curve_csv = args.output_dir / "the240_photon_definition_validation_curves.csv"
    with curve_csv.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(("stage", "lane", "e11e33", "fraction", "stat_error"))
        for stage, stage_curves in curves.items():
            for lane, curve in stage_curves.items():
                for x, y, error in zip(curve.x, curve.y, curve.error):
                    writer.writerow((stage, lane, f"{x:.6g}", f"{y:.9g}", f"{error:.9g}"))

    manifest = {
        "schema_version": 1,
        "task": "THE-240",
        "figure_status": "local_candidate_pending_human_scientific_acceptance",
        "claim": "pp raw-ABCD-purity reproducibility plus qualitative central-AuAu E11/E33 morphology before and after tight photon ID",
        "not_claimed": [
            "independent photon-ID validation",
            "final AuAu purity",
            "systematic uncertainty",
            "unfolding",
            "jet calibration or jet-result validation",
        ],
        "purity": {
            "points_csv": str(PURITY_CSV),
            "points_csv_sha256": sha256(PURITY_CSV),
            "source_manifest": str(PURITY_MANIFEST),
            "source_manifest_sha256": sha256(PURITY_MANIFEST),
            "displayed_series": ["PPG12 raw", "Current output raw"],
            "uncertainty": "statistical/toy uncertainty stored in the frozen point table; no systematic band",
            "leakage_correction_displayed": False,
        },
        "auau": {
            "data_namespace": DATA_DIRECTORY,
            "centrality": CENTRALITY,
            "photon_et_gev": [15, 35],
            "pt_bins": list(PT_BINS),
            "variable": "E11/E33",
            "stages": ["pre", "tight"],
            "display_rebin": DISPLAY_REBIN,
            "display_bin_width": curves["pre"]["data"].bin_width,
            "normalization": "unit-normalized continuous visible shape after omitting exact-zero boundary bin",
            "fill_multiplicity_error_correction": FILL_MULTIPLICITY,
            "uncertainty": "statistical error bars only on data; simulation shown as weighted central-value steps",
            "sources": {
                sample_key: {
                    "pointer": str(pointers[sample_key]),
                    "pointer_sha256": sha256(pointers[sample_key]),
                    "root": str(roots[sample_key]),
                    "root_sha256": sha256(roots[sample_key]),
                    "pointer_payload": pointer_payloads[sample_key],
                }
                for sample_key in roots
            },
            "curves": {
                stage: {
                    lane: {
                        "objects": list(curve.source_objects),
                        "visible_sum_before_unit_normalization": curve.visible_sum,
                        "exact_zero_boundary_sum": curve.exact_zero_sum,
                        "effective_entries": curve.effective_entries,
                    }
                    for lane, curve in stage_curves.items()
                }
                for stage, stage_curves in curves.items()
            },
        },
        "outputs": {
            "png": str(png),
            "png_sha256": sha256(png),
            "curve_csv": str(curve_csv),
            "curve_csv_sha256": sha256(curve_csv),
        },
        "renderer": {
            "path": str(Path(__file__).resolve()),
            "sha256": sha256(Path(__file__).resolve()),
        },
        "slides_mutated": False,
        "production_or_remote_actions": False,
        "visual_qa": {
            "status": "pass" if args.visual_qa_pass else "pending_manual_inspection",
            "checks": (
                [
                    "full-canvas visual inspection",
                    "titles and legends do not overlap",
                    "all uncertainty bars remain inside axes",
                    "statistics-only labeling is explicit",
                    "no systematic band is present",
                ]
                if args.visual_qa_pass
                else []
            ),
        },
    }
    manifest_path = args.output_dir / "the240_photon_definition_validation_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(png)
    print(manifest_path)


if __name__ == "__main__":
    main()
