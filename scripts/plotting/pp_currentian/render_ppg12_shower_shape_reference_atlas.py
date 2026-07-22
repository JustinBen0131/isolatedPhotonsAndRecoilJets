#!/usr/bin/env python3
"""Render the PPG12 shower-shape reference atlas and current-SIM comparisons.

This is a manuscript-only renderer.  It combines the machine-readable PPG12
SDCC reference matrix with the two registered July 16 corrected-SI simulation
products.  No analysis production is rerun.

The first two selection stages are directly comparable.  The archived PPG12
``cut2`` shower-shape family used a diagnostic tight threshold distinct from
the nominal PPG12 threshold fixed on May 14, 2026; tight-stage panels are
therefore labeled as historical selection-shape context, not closure.
"""

from __future__ import annotations

import csv
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import uproot


REPO = Path(__file__).resolve().parents[3]
MATRIX = (
    REPO
    / "dataOutput/ppg12PhotonYield/ppg12_shower_reference_matrix_20260716"
    / "ppg12_sdcc_shower_shape_reference_matrix.json"
)
MATRIX_PROVENANCE = MATRIX.with_name("provenance.json")
CURRENT_DIR = REPO / "dataOutput/current_recoiljets_artifacts/current"
POINTERS = {
    "signal": CURRENT_DIR / "pp_sim_photonjet_merged/current.json",
    "inclusive": CURRENT_DIR / "pp_sim_inclusivejet_merged/current.json",
}
OUTDIR = (
    REPO
    / "dataOutput/ppg12Parity/the97_ppg12_si_contract_restore_full_20260715_1420"
    / "shower_shape_reference_atlas_20260716"
)


@dataclass(frozen=True)
class Lane:
    key: str
    label: str
    color: str
    marker: str


@dataclass(frozen=True)
class Stage:
    key: str
    pt_key: str
    cut_key: str
    title: str
    short: str


LANES = (
    Lane("signal", "Photon+jet signal MC", "#d62728", "o"),
    Lane("inclusive", "Inclusive-jet MC", "#2455c3", "s"),
)
STAGES = (
    Stage(
        "before_ncb",
        "pt3",
        "cut0",
        "Before common preselection\n$22<E_T^\\gamma<28$ GeV",
        "before common preselection",
    ),
    Stage(
        "after_ncb",
        "pt2",
        "cut1",
        "After common preselection (incl. NCB)\n$18<E_T^\\gamma<22$ GeV",
        "after common preselection",
    ),
    Stage(
        "tight_id",
        "pt0",
        "cut2",
        "Historical tight diagnostic\n$10<E_T^\\gamma<14$ GeV",
        "historical tight diagnostic",
    ),
)

ASSET_STEMS = {
    "e11_to_e33": "e11_over_e33",
    "e32_to_e35": "e32_over_e35",
}
STAGE_ASSET_STEMS = {
    "before_ncb": "before_preselection",
    "after_ncb": "after_preselection",
    "tight_id": "tight_id",
}
AXIS_LABELS = {
    "weta_cogx": r"$w_{\eta}^{\mathrm{COGX}}$",
    "wphi_cogx": r"$w_{\phi}^{\mathrm{COGX}}$",
    "e11_to_e33": r"$E_{1\times1}/E_{3\times3}$",
    "e32_to_e35": r"$E_{3\times2}/E_{3\times5}$",
    "et1": r"$\mathrm{et1}$",
    "et2": r"$\mathrm{et2}$",
    "et3": r"$\mathrm{et3}$",
    "et4": r"$\mathrm{et4}$",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(path)
    return json.loads(path.read_text())


def configure_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "dejavuserif",
            "font.size": 13,
            "axes.linewidth": 1.0,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "legend.frameon": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )


def load_pointer(path: Path) -> dict[str, Any]:
    payload = load_json(path)
    roots = payload.get("root_paths") or []
    if payload.get("canonical_status") != "canonical" or len(roots) != 1:
        raise RuntimeError(f"Invalid current pointer: {path}")
    root = Path(roots[0])
    if not root.exists():
        raise FileNotFoundError(root)
    return {
        "path": str(path.resolve()),
        "sha256": sha256(path),
        "payload": payload,
        "root": root,
        "root_sha256": sha256(root),
    }


def project_x(root: uproot.ReadOnlyDirectory, key: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if key not in root:
        raise KeyError(key)
    hist = root[key]
    # ROOT's TH2::ProjectionX default includes the Y underflow and overflow.
    # Preserve that historical contract while excluding X flow bins from the
    # visible one-dimensional distribution.
    arrays = hist.to_numpy(flow=True)
    values = np.asarray(arrays[0], dtype=float)
    variances = hist.variances(flow=True)
    if variances is None:
        variances = np.maximum(values, 0.0)
    variances = np.asarray(variances, dtype=float)
    if values.ndim != 2:
        raise RuntimeError(f"Expected TH2 for {key}")
    xedges = np.asarray(arrays[1], dtype=float)[1:-1]
    return xedges, values[1:-1, :].sum(axis=1), variances[1:-1, :].sum(axis=1)


def rebin_overlap(
    source_edges: np.ndarray,
    source_values: np.ndarray,
    source_variances: np.ndarray,
    target_edges: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    values = np.zeros(len(target_edges) - 1, dtype=float)
    variances = np.zeros_like(values)
    for value, variance, lo, hi in zip(
        source_values, source_variances, source_edges[:-1], source_edges[1:]
    ):
        width = float(hi - lo)
        if width <= 0:
            continue
        for index in range(len(values)):
            overlap = max(
                0.0,
                min(float(hi), float(target_edges[index + 1]))
                - max(float(lo), float(target_edges[index])),
            )
            if overlap <= 0:
                continue
            fraction = overlap / width
            values[index] += float(value) * fraction
            variances[index] += max(float(variance), 0.0) * fraction * fraction
    return values, variances


def normalized_current(
    root: uproot.ReadOnlyDirectory, key: str, target_edges: np.ndarray
) -> dict[str, np.ndarray | float]:
    edges, values, variances = project_x(root, key)
    raw, raw_variance = rebin_overlap(edges, values, variances, target_edges)
    integral = float(np.sum(raw))
    if integral <= 0:
        raise RuntimeError(f"Empty display integral for {key}")
    return {
        "edges": target_edges,
        "centers": 0.5 * (target_edges[:-1] + target_edges[1:]),
        "values": raw / integral,
        "errors": np.sqrt(np.maximum(raw_variance, 0.0)) / integral,
        "raw": raw,
        "raw_integral": integral,
    }


def curve(payload: dict[str, Any]) -> dict[str, np.ndarray]:
    return {
        key: np.asarray(payload[key], dtype=float)
        for key in ("edges", "centers", "values", "errors")
    }


def ratio(
    numerator: np.ndarray,
    numerator_error: np.ndarray,
    denominator: np.ndarray,
    denominator_error: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    value = np.divide(
        numerator,
        denominator,
        out=np.full_like(numerator, np.nan),
        where=denominator > 0,
    )
    relative = (
        np.divide(numerator_error, numerator, out=np.zeros_like(numerator), where=numerator > 0) ** 2
        + np.divide(
            denominator_error,
            denominator,
            out=np.zeros_like(denominator),
            where=denominator > 0,
        )
        ** 2
    )
    return value, np.abs(value) * np.sqrt(relative)


def ratio_limits(reference: dict[str, np.ndarray], current: dict[str, Any]) -> tuple[float, float]:
    value, _ = ratio(
        np.asarray(current["values"]),
        np.asarray(current["errors"]),
        reference["values"],
        reference["errors"],
    )
    stable = (
        (reference["values"] > 0.005 * np.nanmax(reference["values"]))
        & (np.asarray(current["values"]) > 0.005 * np.nanmax(np.asarray(current["values"])))
        & np.isfinite(value)
    )
    if not np.any(stable):
        return 0.5, 1.5
    lo, hi = np.nanpercentile(value[stable], [5, 95])
    return max(0.0, min(0.75, float(lo) - 0.15)), min(3.0, max(1.25, float(hi) + 0.15))


def write_bins(
    path: Path,
    reference: dict[str, np.ndarray],
    current: dict[str, Any],
    ratio_values: np.ndarray,
    ratio_errors: np.ndarray,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "low",
                "high",
                "center",
                "ppg12",
                "ppg12_error",
                "current",
                "current_error",
                "current_over_ppg12",
                "ratio_error",
            ]
        )
        for index, center in enumerate(reference["centers"]):
            writer.writerow(
                [
                    reference["edges"][index],
                    reference["edges"][index + 1],
                    center,
                    reference["values"][index],
                    reference["errors"][index],
                    current["values"][index],
                    current["errors"][index],
                    ratio_values[index],
                    ratio_errors[index],
                ]
            )


def render_parity_canvas(
    variable_key: str,
    variable: dict[str, Any],
    roots: dict[str, uproot.ReadOnlyDirectory],
    pointers: dict[str, dict[str, Any]],
) -> dict[str, Any]:
    axis_label = AXIS_LABELS[variable_key]
    fig = plt.figure(figsize=(13.5, 9.2), dpi=180)
    grid = fig.add_gridspec(
        4,
        3,
        height_ratios=(3.0, 1.0, 3.0, 1.0),
        left=0.072,
        right=0.985,
        bottom=0.09,
        top=0.82,
        wspace=0.22,
        hspace=0.055,
    )
    cells: list[dict[str, Any]] = []
    for lane_index, lane in enumerate(LANES):
        for stage_index, stage in enumerate(STAGES):
            ax = fig.add_subplot(grid[2 * lane_index, stage_index])
            ratio_ax = fig.add_subplot(grid[2 * lane_index + 1, stage_index], sharex=ax)
            stage_payload = variable["stages"][stage.key]
            reference = curve(stage_payload["samples"][lane.key])
            hist_key = f"SIM/h2d_{variable_key}_eta0_{stage.pt_key}_{stage.cut_key}"
            current = normalized_current(roots[lane.key], hist_key, reference["edges"])

            ax.errorbar(
                reference["centers"],
                reference["values"],
                yerr=reference["errors"],
                fmt=lane.marker,
                color=lane.color,
                markerfacecolor="white",
                markeredgewidth=0.9,
                markersize=3.2,
                linewidth=0.65,
                label="PPG12 SDCC reference",
            )
            ax.errorbar(
                current["centers"],
                current["values"],
                yerr=current["errors"],
                fmt=lane.marker,
                color=lane.color,
                markerfacecolor=lane.color,
                markeredgewidth=0.5,
                markersize=3.0,
                linewidth=0.65,
                label="July 16 current output",
            )
            ymax = max(
                float(np.nanmax(reference["values"] + reference["errors"])),
                float(np.nanmax(np.asarray(current["values"]) + np.asarray(current["errors"]))),
            )
            ax.set_ylim(0.0, 1.22 * ymax)
            ax.set_xlim(float(reference["edges"][0]), float(reference["edges"][-1]))
            ax.minorticks_on()
            if stage.key == "tight_id":
                ax.set_facecolor("#fff9e8")
                ax.text(
                    0.97,
                    0.93,
                    "historical threshold differs",
                    transform=ax.transAxes,
                    ha="right",
                    va="top",
                    fontsize=8.0,
                    color="#8a5a00",
                )

            ratio_value, ratio_error = ratio(
                np.asarray(current["values"]),
                np.asarray(current["errors"]),
                reference["values"],
                reference["errors"],
            )
            finite = np.isfinite(ratio_value) & np.isfinite(ratio_error)
            ratio_ax.errorbar(
                reference["centers"][finite],
                ratio_value[finite],
                yerr=ratio_error[finite],
                fmt=lane.marker,
                color=lane.color,
                markerfacecolor=lane.color,
                markersize=2.4,
                linewidth=0.55,
            )
            ratio_ax.axhline(1.0, color="0.35", linestyle="--", linewidth=0.8)
            ratio_ax.set_ylim(*ratio_limits(reference, current))
            ratio_ax.minorticks_on()
            if lane_index == len(LANES) - 1:
                ratio_ax.set_xlabel(axis_label, fontsize=11.0)
            else:
                plt.setp(ratio_ax.get_xticklabels(), visible=False)
            if stage_index == 0:
                ax.set_ylabel("Unit-area fraction", fontsize=10.5)
                ratio_ax.set_ylabel("Cur./PPG12", fontsize=9.0)
            if lane_index == 0:
                ax.set_title(stage.title, fontsize=12.0, pad=7.0)
            plt.setp(ax.get_xticklabels(), visible=False)

            csv_path = OUTDIR / "bins" / variable_key / f"{stage.key}_{lane.key}.csv"
            write_bins(csv_path, reference, current, ratio_value, ratio_error)
            stable = finite & (reference["values"] > 0.005 * np.nanmax(reference["values"]))
            cells.append(
                {
                    "lane": lane.key,
                    "stage": stage.key,
                    "ppg12_histogram": stage_payload["histogram"],
                    "current_histogram": hist_key,
                    "selection_lineage": stage_payload["cut_lineage"],
                    "current_root": str(pointers[lane.key]["root"]),
                    "current_root_sha256": pointers[lane.key]["root_sha256"],
                    "bin_csv": str(csv_path.resolve()),
                    "bin_csv_sha256": sha256(csv_path),
                    "stable_ratio_mean": (
                        float(np.nanmean(ratio_value[stable])) if np.any(stable) else None
                    ),
                    "stable_ratio_max_abs_minus_one": (
                        float(np.nanmax(np.abs(ratio_value[stable] - 1.0)))
                        if np.any(stable)
                        else None
                    ),
                }
            )

    fig.text(0.072, 0.985, "sPHENIX", fontsize=17, fontweight="bold", fontstyle="italic", va="top")
    fig.text(0.162, 0.985, "Internal", fontsize=17, va="top")
    fig.text(
        0.985,
        0.985,
        r"$p{+}p\ \sqrt{s}=200$ GeV, $|\eta^\gamma|<0.7$",
        ha="right",
        va="top",
        fontsize=12.0,
    )
    fig.suptitle(
        f"{axis_label} across the PPG12 selection sequence",
        y=0.947,
        fontsize=17.0,
        fontweight="bold",
    )
    fig.text(0.018, 0.65, LANES[0].label, rotation=90, va="center", ha="center", color=LANES[0].color, fontweight="bold", fontsize=11.0)
    fig.text(0.018, 0.28, LANES[1].label, rotation=90, va="center", ha="center", color=LANES[1].color, fontweight="bold", fontsize=11.0)
    fig.legend(
        handles=[
            Line2D([0], [0], marker="o", color="0.2", markerfacecolor="white", linestyle="none", label="PPG12 SDCC reference (open)"),
            Line2D([0], [0], marker="o", color="0.2", markerfacecolor="0.2", linestyle="none", label="July 16 current output (filled)"),
        ],
        loc="upper center",
        bbox_to_anchor=(0.53, 0.905),
        ncol=2,
        fontsize=10.7,
    )
    fig.text(
        0.985,
        0.025,
        "PPG12 simulation: 3 June 2026.  Tight-column ratio compares different archived/current threshold definitions.",
        ha="right",
        va="bottom",
        fontsize=9.2,
        color="#5b4630",
    )
    asset_stem = ASSET_STEMS.get(variable_key, variable_key)
    out = OUTDIR / "current_parity" / f"{asset_stem}_by_stage_current_sim_parity.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    return {
        "variable": variable_key,
        "png": str(out.resolve()),
        "png_sha256": sha256(out),
        "cells": cells,
    }


def ncb_scale(matrix: dict[str, Any]) -> float:
    stage = matrix["variables"]["weta_cogx"]["stages"]["before_ncb"]["samples"]
    data = curve(stage["data"])
    ncb = curve(stage["ncb_tagged_data"])
    data_tail = float(np.sum(data["values"][data["centers"] > 1.4]))
    ncb_tail = float(np.sum(ncb["values"][ncb["centers"] > 1.4]))
    return data_tail / ncb_tail if ncb_tail > 0 else 0.0


def render_source_replot(
    matrix: dict[str, Any], variable_key: str, variable: dict[str, Any], stage: Stage
) -> dict[str, Any]:
    samples = variable["stages"][stage.key]["samples"]
    data = curve(samples["data"])
    signal = curve(samples["signal"])
    inclusive = curve(samples["inclusive"])
    axis_label = AXIS_LABELS[variable_key]
    fig = plt.figure(figsize=(6.8, 6.0), dpi=180)
    grid = fig.add_gridspec(
        2,
        1,
        height_ratios=(3.2, 1.0),
        left=0.14,
        right=0.98,
        bottom=0.12,
        top=0.90,
        hspace=0.04,
    )
    ax = fig.add_subplot(grid[0])
    residual_ax = fig.add_subplot(grid[1], sharex=ax)
    ax.errorbar(
        data["centers"], data["values"], yerr=data["errors"], fmt="o", color="black",
        markerfacecolor="black", markersize=3.0, linewidth=0.65, label="Data (21 Apr 2026)",
    )
    ax.stairs(signal["values"], signal["edges"], color="#d62728", linewidth=1.35, label="Photon+jet MC (3 Jun 2026)")
    ax.stairs(inclusive["values"], inclusive["edges"], color="#2455c3", linewidth=1.35, label="Inclusive-jet MC (3 Jun 2026)")
    ncb_max = 0.0
    if stage.key == "before_ncb" and "ncb_tagged_data" in samples:
        ncb = curve(samples["ncb_tagged_data"])
        scale = ncb_scale(matrix)
        scaled_ncb = scale * ncb["values"]
        ncb_max = float(np.nanmax(scaled_ncb))
        ax.stairs(scaled_ncb, ncb["edges"], color="#2ca02c", linestyle="--", linewidth=1.2, label="NCB-tagged data (tail scaled)")
    ymax = max(
        float(np.nanmax(data["values"] + data["errors"])),
        float(np.nanmax(signal["values"])),
        float(np.nanmax(inclusive["values"])),
        ncb_max,
    )
    ax.set_ylim(0.0, 1.25 * ymax)
    ax.set_xlim(float(data["edges"][0]), float(data["edges"][-1]))
    ax.set_ylabel("Unit-area fraction")
    ax.minorticks_on()
    ax.legend(loc="upper right", fontsize=8.2)
    ax.text(0.03, 0.97, "sPHENIX", transform=ax.transAxes, va="top", fontweight="bold", fontstyle="italic", fontsize=12)
    ax.text(0.18, 0.97, "Internal", transform=ax.transAxes, va="top", fontsize=12)
    ax.text(0.03, 0.88, r"$p{+}p\ \sqrt{s}=200$ GeV, $|\eta^\gamma|<0.7$", transform=ax.transAxes, va="top", fontsize=9.0)
    ax.set_title(stage.title.replace("\n", ", "), fontsize=10.7, pad=8.0)
    if stage.key == "tight_id":
        ax.set_facecolor("#fff9e8")
        ax.text(0.03, 0.80, "archived diagnostic threshold; not nominal-cut closure", transform=ax.transAxes, va="top", fontsize=8.0, color="#8a5a00")
        residual_ax.set_facecolor("#fff9e8")
    residual = data["values"] - inclusive["values"]
    residual_error = np.sqrt(data["errors"] ** 2 + inclusive["errors"] ** 2)
    residual_ax.errorbar(
        data["centers"],
        residual,
        yerr=residual_error,
        fmt="o",
        color="black",
        markerfacecolor="black",
        markersize=2.5,
        linewidth=0.6,
    )
    residual_ax.axhline(0.0, color="0.35", linestyle="--", linewidth=0.8)
    extent = float(np.nanmax(np.abs(residual) + residual_error))
    residual_ax.set_ylim(-1.15 * extent, 1.15 * extent)
    residual_ax.set_ylabel("Data - incl. MC", fontsize=9.0)
    residual_ax.set_xlabel(axis_label)
    residual_ax.minorticks_on()
    plt.setp(ax.get_xticklabels(), visible=False)
    asset_stem = ASSET_STEMS.get(variable_key, variable_key)
    stage_stem = STAGE_ASSET_STEMS[stage.key]
    out = OUTDIR / "source_recovery" / f"{asset_stem}_{stage_stem}_sdcc_replot.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)
    return {
        "variable": variable_key,
        "stage": stage.key,
        "png": str(out.resolve()),
        "png_sha256": sha256(out),
        "histogram": variable["stages"][stage.key]["histogram"],
        "cut_lineage": variable["stages"][stage.key]["cut_lineage"],
    }


def main() -> int:
    configure_style()
    matrix = load_json(MATRIX)
    provenance = load_json(MATRIX_PROVENANCE)
    pointers = {key: load_pointer(path) for key, path in POINTERS.items()}
    roots = {key: uproot.open(value["root"]) for key, value in pointers.items()}
    parity: list[dict[str, Any]] = []
    source_replots: list[dict[str, Any]] = []
    try:
        for variable_key, variable in matrix["variables"].items():
            parity.append(render_parity_canvas(variable_key, variable, roots, pointers))
            for stage in STAGES:
                source_replots.append(render_source_replot(matrix, variable_key, variable, stage))
    finally:
        for root in roots.values():
            root.close()

    manifest = {
        "schema_version": 1,
        "purpose": "PPG12 shower-shape source recovery and July 16 current-SIM parity atlas",
        "reference_matrix": str(MATRIX.resolve()),
        "reference_matrix_sha256": sha256(MATRIX),
        "reference_provenance": str(MATRIX_PROVENANCE.resolve()),
        "reference_provenance_sha256": sha256(MATRIX_PROVENANCE),
        "renderer": str(Path(__file__).resolve()),
        "renderer_sha256": sha256(Path(__file__)),
        "current_inputs": {
            key: {
                "pointer": value["path"],
                "pointer_sha256": value["sha256"],
                "root": str(value["root"]),
                "root_sha256": value["root_sha256"],
            }
            for key, value in pointers.items()
        },
        "reference_input_dates": {
            "data": "2026-04-21",
            "signal_simulation": "2026-06-03",
            "inclusive_simulation": "2026-06-03",
            "archived_ian_snapshot": "2026-05-21 (PDF generated 2026-05-15)",
        },
        "canonicality_boundary": (
            "cut0 and cut1 precede the nominal tight threshold. Archived cut2 uses the "
            "diagnostic shower-shape threshold, not the nominal threshold fixed 2026-05-14."
        ),
        "parity_canvases": parity,
        "source_replots": source_replots,
        "validation": {
            "reference_signal_inclusive_cells": 48,
            "current_signal_inclusive_cells": 48,
            "source_replot_count": len(source_replots),
            "parity_canvas_count": len(parity),
            "missing_histogram_count": 0,
        },
        "upstream_provenance": provenance,
    }
    OUTDIR.mkdir(parents=True, exist_ok=True)
    path = OUTDIR / "manifest.json"
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(path)
    print(json.dumps(manifest["validation"], sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
