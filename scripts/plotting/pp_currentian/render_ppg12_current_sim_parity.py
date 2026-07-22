#!/usr/bin/env python3
"""Render current pp-SIM PPG12-parity diagnostics for the IAN.

The current photon+jet and inclusive-jet ROOT files are always resolved from
their registered ``current.json`` pointers.  This renderer covers the BDT
score and the eight pp baseV3E shower-shape inputs in the three PPG12 audit
stages used by the IAN:

* before the non-collision-background (NCB) preselection, 22--28 GeV;
* after the NCB preselection, 18--22 GeV;
* after the tight photon-ID selection, 10--14 GeV.

When a local machine-readable PPG12 signal/inclusive projection is available,
the output is a unit-area PPG12/current overlay with a Current/PPG12 ratio.
When it is not available, the current histogram is still rendered and the
panel says ``PPG12 numerical reference pending``.  Thus a reference-recovery
gap is never misreported as a missing current histogram.

The same command also regenerates the absolute signal-MC ABCD correction-input
comparison.  Its PPG12 side comes from the locally stored ROOT-extracted CSV;
its current side is read directly from the registered photon+jet ROOT.

Run with the repository analysis Python environment, which provides uproot::

    /Users/patsfan753/Desktop/analysis/env/bin/python3 \
      scripts/plotting/pp_currentian/render_ppg12_current_sim_parity.py
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import uproot


REPO = Path(__file__).resolve().parents[3]
CURRENT_DIR = REPO / "dataOutput/current_recoiljets_artifacts/current"
PHOTON_POINTER = CURRENT_DIR / "pp_sim_photonjet_merged/current.json"
INCLUSIVE_POINTER = CURRENT_DIR / "pp_sim_inclusivejet_merged/current.json"
REFERENCE_BASE = (
    REPO
    / "dataOutput/ppg12PhotonYield/ppg12_photon_yield_v1_data_20260620"
    / "shower_shape_reference_validation"
)
BDT_FIG13_19_JSON = (
    REFERENCE_BASE / "fig13_19_bdt/ppg12_sdcc_fig19_bdt_fourcurve_projections.json"
)
BDT_FIG20_JSON = (
    REFERENCE_BASE / "fig20_bdt/ppg12_sdcc_fig20_bdt_tight_fourcurve_projections.json"
)
E11_FIG13_JSON = (
    REFERENCE_BASE / "fig13_e11_e33/ppg12_sdcc_fig13_e11_to_e33_histograms.json"
)
WETA_FIG13_JSON = (
    REFERENCE_BASE / "fig13_weta_cogx/ppg12_sdcc_fig13_weta_cogx_histograms.json"
)
ABCD_REFERENCE_DIR = (
    REPO
    / "dataOutput/ppg12Parity/the97_ppg12_si_contract_restore_full_20260715_1420"
    / "corrected_si_fullstat_candidate"
)
ABCD_REFERENCE_CSV = ABCD_REFERENCE_DIR / "candidate_signal_abcd_points.csv"
ABCD_REFERENCE_MANIFEST = ABCD_REFERENCE_DIR / "candidate_signal_abcd_extract_manifest.json"
DEFAULT_OUTDIR = (
    REPO
    / "dataOutput/ppg12Parity/the97_ppg12_si_contract_restore_full_20260715_1420"
    / "ian_current_sim_refresh_20260716"
)


@dataclass(frozen=True)
class StageSpec:
    key: str
    pt_key: str
    cut_key: str
    selection: str
    short_label: str


@dataclass(frozen=True)
class LaneSpec:
    key: str
    label: str
    pointer: Path
    ppg12_curve: str
    color: str
    marker: str


@dataclass(frozen=True)
class VariableSpec:
    key: str
    label: str
    current_edges: tuple[float, ...]


STAGES = (
    StageSpec(
        "before_ncb",
        "pt3",
        "cut0",
        r"$22 < E_T^\gamma < 28$ GeV, before NCB preselection",
        "before NCB preselection",
    ),
    StageSpec(
        "after_ncb",
        "pt2",
        "cut1",
        r"$18 < E_T^\gamma < 22$ GeV, after NCB preselection",
        "after NCB preselection",
    ),
    StageSpec(
        "tight_id",
        "pt0",
        "cut2",
        r"$10 < E_T^\gamma < 14$ GeV, tight photon-ID selection",
        "tight photon-ID selection",
    ),
)

LANES = (
    LaneSpec(
        "photonjet_signal",
        "Photon+jet signal MC",
        PHOTON_POINTER,
        "signal_mc",
        "#d62728",
        "o",
    ),
    LaneSpec(
        "inclusivejet",
        "Inclusive-jet MC",
        INCLUSIVE_POINTER,
        "inclusive_mc",
        "#2455c3",
        "s",
    ),
)


def uniform_edges(start: float, stop: float, step: float) -> tuple[float, ...]:
    count = int(round((stop - start) / step))
    return tuple(float(x) for x in np.linspace(start, stop, count + 1))


SHOWER_VARIABLES = (
    VariableSpec("weta_cogx", r"$w_{\eta}^{\mathrm{cog}}$", uniform_edges(0.0, 2.0, 0.04)),
    VariableSpec("wphi_cogx", r"$w_{\phi}^{\mathrm{cog}}$", uniform_edges(0.0, 2.0, 0.04)),
    VariableSpec("e11_to_e33", r"$E_{1\times1}/E_{3\times3}$", uniform_edges(0.0, 1.0, 0.02)),
    VariableSpec("e32_to_e35", r"$E_{3\times2}/E_{3\times5}$", uniform_edges(0.0, 1.0, 0.02)),
    VariableSpec("et1", r"$E_{T,1}$ sharing variable", uniform_edges(0.0, 1.0, 0.02)),
    VariableSpec("et2", r"$E_{T,2}$ sharing variable", uniform_edges(-0.4, 1.0, 0.02)),
    VariableSpec("et3", r"$E_{T,3}$ sharing variable", uniform_edges(-0.4, 1.0, 0.02)),
    VariableSpec("et4", r"$E_{T,4}$ sharing variable", uniform_edges(0.0, 0.6, 0.02)),
)

ABCD_REGIONS = (
    ("A", "tight isolated", "black", "o", "SIM/h_tight_iso_cluster_signal_0"),
    ("B", "tight nonisolated", "#d62728", "s", "SIM/h_tight_noniso_cluster_signal_0"),
    ("C", "nontight isolated", "#1f77b4", "^", "SIM/h_nontight_iso_cluster_signal_0"),
    ("D", "nontight nonisolated", "#9467bd", "D", "SIM/h_nontight_noniso_cluster_signal_0"),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument(
        "--families",
        nargs="+",
        choices=("bdt", "shower", "abcd"),
        default=("bdt", "shower", "abcd"),
        help="Output families to render.",
    )
    return parser.parse_args()


def configure_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "dejavuserif",
            "font.size": 12,
            "axes.linewidth": 1.1,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "legend.frameon": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(path)
    return json.loads(path.read_text())


def load_pointer(path: Path) -> dict[str, Any]:
    payload = load_json(path)
    roots = payload.get("root_paths") or []
    if len(roots) != 1:
        raise RuntimeError(f"Expected one root_paths entry in {path}, found {roots}")
    root_path = Path(roots[0])
    if not root_path.exists():
        raise FileNotFoundError(root_path)
    if payload.get("canonical_status") != "canonical":
        raise RuntimeError(f"Pointer is not canonical: {path}")
    return {
        "path": str(path.resolve()),
        "sha256": sha256_file(path),
        "payload": payload,
        "root_path": str(root_path.resolve()),
        "root_sha256": sha256_file(root_path),
    }


def project_x(root_file: uproot.ReadOnlyDirectory, hist_key: str) -> dict[str, Any]:
    if hist_key not in root_file:
        raise KeyError(f"Current ROOT is missing required histogram {hist_key}")
    hist = root_file[hist_key]
    arrays = hist.to_numpy(flow=False)
    values = np.asarray(arrays[0], dtype=float)
    variances = hist.variances(flow=False)
    if variances is None:
        variances = np.maximum(values, 0.0)
    variances = np.asarray(variances, dtype=float)
    if values.ndim == 2:
        xedges = np.asarray(arrays[1], dtype=float)
        yedges = np.asarray(arrays[2], dtype=float)
        projected = values.sum(axis=1)
        projected_variance = variances.sum(axis=1)
        projection = "TH2 ProjectionX over all non-flow y bins"
        y_axis = [float(yedges[0]), float(yedges[-1])]
    elif values.ndim == 1:
        xedges = np.asarray(arrays[1], dtype=float)
        projected = values
        projected_variance = variances
        projection = "TH1 direct"
        y_axis = None
    else:
        raise RuntimeError(f"Unsupported dimension {values.ndim} for {hist_key}")
    return {
        "edges": xedges,
        "values": projected,
        "variances": projected_variance,
        "source_integral": float(np.sum(projected)),
        "entries": float(hist.member("fEntries")),
        "projection": projection,
        "source_x_axis": [float(xedges[0]), float(xedges[-1])],
        "source_y_axis": y_axis,
    }


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
        first = max(0, int(np.searchsorted(target_edges, lo, side="right") - 1))
        last = min(len(values) - 1, int(np.searchsorted(target_edges, hi, side="left")))
        for index in range(first, last + 1):
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


def normalized_current_curve(
    root_file: uproot.ReadOnlyDirectory,
    hist_key: str,
    target_edges: np.ndarray,
) -> dict[str, Any]:
    projected = project_x(root_file, hist_key)
    raw, raw_variance = rebin_overlap(
        projected["edges"],
        projected["values"],
        projected["variances"],
        target_edges,
    )
    integral = float(np.sum(raw))
    if integral <= 0:
        raise RuntimeError(f"Cannot normalize empty target range for {hist_key}")
    return {
        "edges": target_edges,
        "centers": 0.5 * (target_edges[:-1] + target_edges[1:]),
        "values": raw / integral,
        "errors": np.sqrt(np.maximum(raw_variance, 0.0)) / integral,
        "raw": raw,
        "raw_errors": np.sqrt(np.maximum(raw_variance, 0.0)),
        "target_integral": integral,
        **{key: value for key, value in projected.items() if key not in {"edges", "values", "variances"}},
    }


def reference_curve(payload: dict[str, Any]) -> dict[str, np.ndarray]:
    edges = np.asarray(payload["edges"], dtype=float)
    values = np.asarray(payload["values"], dtype=float)
    errors = np.asarray(payload["errors"], dtype=float)
    if len(edges) != len(values) + 1 or len(errors) != len(values):
        raise RuntimeError("Malformed PPG12 numerical reference")
    return {
        "edges": edges,
        "centers": 0.5 * (edges[:-1] + edges[1:]),
        "values": values,
        "errors": errors,
    }


def ratio_with_error(
    numerator: np.ndarray,
    numerator_error: np.ndarray,
    denominator: np.ndarray,
    denominator_error: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    ratio = np.divide(
        numerator,
        denominator,
        out=np.full_like(numerator, np.nan),
        where=denominator > 0,
    )
    relative_variance = (
        np.divide(
            numerator_error,
            numerator,
            out=np.zeros_like(numerator),
            where=numerator > 0,
        )
        ** 2
        + np.divide(
            denominator_error,
            denominator,
            out=np.zeros_like(denominator),
            where=denominator > 0,
        )
        ** 2
    )
    return ratio, np.abs(ratio) * np.sqrt(relative_variance)


def sphinx_label(ax: plt.Axes, x: float = 0.045, y: float = 0.955) -> None:
    ax.text(
        x,
        y,
        "sPHENIX",
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=14,
        fontweight="bold",
        fontstyle="italic",
    )
    ax.text(x + 0.205, y, "Internal", transform=ax.transAxes, va="top", ha="left", fontsize=14)


def figure_header(fig: plt.Figure, lane: LaneSpec, stage: StageSpec) -> None:
    """Place common physics labels in a reserved header above the data pad."""
    fig.text(
        0.135,
        0.978,
        "sPHENIX",
        va="top",
        ha="left",
        fontsize=14,
        fontweight="bold",
        fontstyle="italic",
    )
    fig.text(0.310, 0.978, "Internal", va="top", ha="left", fontsize=14)
    fig.text(
        0.135,
        0.943,
        r"$p{+}p\ \sqrt{s}=200$ GeV, $|\eta^\gamma|<0.7$",
        va="top",
        ha="left",
        fontsize=10.5,
    )
    fig.text(0.98, 0.943, stage.selection, va="top", ha="right", fontsize=10.0)
    fig.text(0.98, 0.910, lane.label, va="top", ha="right", fontsize=10.0)


def choose_ratio_axis(ratio: np.ndarray, stable: np.ndarray) -> tuple[str, tuple[float, float]]:
    populated = ratio[stable & np.isfinite(ratio) & (ratio > 0)]
    if not populated.size:
        return "linear", (0.0, 2.0)
    low = float(np.nanmin(populated))
    high = float(np.nanmax(populated))
    if high / max(low, 1e-12) > 15.0 or high > 5.0:
        return "log", (max(0.03, 0.75 * low), min(500.0, max(2.0, 1.35 * high)))
    return "linear", (max(0.0, min(0.75, 0.85 * low)), max(1.35, 1.15 * high))


def write_overlay_csv(
    path: Path,
    reference: dict[str, np.ndarray],
    current: dict[str, Any],
    ratio: np.ndarray,
    ratio_error: np.ndarray,
) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "bin",
                "low",
                "high",
                "center",
                "ppg12",
                "ppg12_error",
                "current",
                "current_error",
                "current_over_ppg12",
                "current_over_ppg12_error",
                "current_raw",
            ]
        )
        for index, center in enumerate(reference["centers"]):
            writer.writerow(
                [
                    index + 1,
                    f"{reference['edges'][index]:.12g}",
                    f"{reference['edges'][index + 1]:.12g}",
                    f"{center:.12g}",
                    f"{reference['values'][index]:.12g}",
                    f"{reference['errors'][index]:.12g}",
                    f"{current['values'][index]:.12g}",
                    f"{current['errors'][index]:.12g}",
                    f"{ratio[index]:.12g}",
                    f"{ratio_error[index]:.12g}",
                    f"{current['raw'][index]:.12g}",
                ]
            )


def draw_overlay_ratio(
    out_png: Path,
    out_csv: Path,
    variable_label: str,
    lane: LaneSpec,
    stage: StageSpec,
    reference: dict[str, np.ndarray],
    current: dict[str, Any],
) -> dict[str, Any]:
    ratio, ratio_error = ratio_with_error(
        current["values"],
        current["errors"],
        reference["values"],
        reference["errors"],
    )
    stable = (reference["values"] > 0.002) & (current["values"] > 0.002)
    ratio_scale, ratio_limits = choose_ratio_axis(ratio, stable)

    fig, (ax, ratio_ax) = plt.subplots(
        2,
        1,
        figsize=(7.4, 7.8),
        dpi=160,
        sharex=True,
        gridspec_kw={"height_ratios": [3.1, 1.0], "hspace": 0.04},
    )
    centers = reference["centers"]
    ax.errorbar(
        centers,
        reference["values"],
        yerr=reference["errors"],
        fmt=lane.marker,
        color=lane.color,
        markerfacecolor="white",
        markeredgewidth=1.0,
        markersize=4.2,
        linewidth=0.8,
        label="PPG12 SDCC reference",
    )
    ax.errorbar(
        centers,
        current["values"],
        yerr=current["errors"],
        fmt=lane.marker,
        color=lane.color,
        markerfacecolor=lane.color,
        markeredgewidth=0.7,
        markersize=3.3,
        linewidth=0.75,
        label="July 16 registered current output",
    )
    figure_header(fig, lane, stage)
    ymax = max(
        float(np.nanmax(reference["values"] + reference["errors"])),
        float(np.nanmax(current["values"] + current["errors"])),
    )
    ax.set_ylim(0.0, 1.20 * ymax)
    ax.set_ylabel("Unit-area fraction")
    ax.legend(loc="upper right", fontsize=9.6, handlelength=1.2)
    ax.tick_params(labelbottom=False)
    ax.minorticks_on()

    finite = np.isfinite(ratio) & np.isfinite(ratio_error)
    ratio_ax.errorbar(
        centers[finite],
        ratio[finite],
        yerr=ratio_error[finite],
        fmt=lane.marker,
        color="black",
        markerfacecolor="black",
        markersize=2.8,
        linewidth=0.65,
    )
    ratio_ax.axhline(1.0, color="0.4", linestyle="--", linewidth=0.9)
    ratio_ax.set_ylabel("Current / PPG12")
    ratio_ax.set_xlabel(variable_label)
    ratio_ax.set_yscale(ratio_scale)
    ratio_ax.set_ylim(*ratio_limits)
    ratio_ax.minorticks_on()
    if variable_label == "BDT score" and stage.key == "tight_id":
        # The PPG12 tight-stage JSON preserves the generic -1..2 storage axis,
        # but both samples are populated only over the selected score range.
        # Crop the display without changing the normalization or source bins.
        ax.set_xlim(0.50, 1.02)
    else:
        ax.set_xlim(float(reference["edges"][0]), float(reference["edges"][-1]))
    fig.subplots_adjust(left=0.135, right=0.98, top=0.875, bottom=0.10, hspace=0.04)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    write_overlay_csv(out_csv, reference, current, ratio, ratio_error)

    stable_ratio = ratio[stable & np.isfinite(ratio)]
    return {
        "display_mode": "ppg12_current_overlay_ratio",
        "normalization": "Each projection is normalized to unit area over the displayed PPG12 reference binning.",
        "ratio_definition": "Current / PPG12",
        "ratio_scale": ratio_scale,
        "stable_ratio_bin_count": int(stable_ratio.size),
        "stable_ratio_mean": float(np.mean(stable_ratio)) if stable_ratio.size else None,
        "stable_ratio_max_abs_minus_one": (
            float(np.max(np.abs(stable_ratio - 1.0))) if stable_ratio.size else None
        ),
    }


def write_current_only_csv(path: Path, current: dict[str, Any]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["bin", "low", "high", "center", "current", "current_error", "current_raw"])
        for index, center in enumerate(current["centers"]):
            writer.writerow(
                [
                    index + 1,
                    f"{current['edges'][index]:.12g}",
                    f"{current['edges'][index + 1]:.12g}",
                    f"{center:.12g}",
                    f"{current['values'][index]:.12g}",
                    f"{current['errors'][index]:.12g}",
                    f"{current['raw'][index]:.12g}",
                ]
            )


def draw_current_only(
    out_png: Path,
    out_csv: Path,
    variable_label: str,
    lane: LaneSpec,
    stage: StageSpec,
    current: dict[str, Any],
) -> dict[str, Any]:
    fig, ax = plt.subplots(figsize=(7.4, 5.8), dpi=160)
    ax.errorbar(
        current["centers"],
        current["values"],
        yerr=current["errors"],
        fmt=lane.marker,
        color=lane.color,
        markerfacecolor=lane.color,
        markeredgewidth=0.6,
        markersize=3.4,
        linewidth=0.75,
        label="July 16 registered current output",
    )
    figure_header(fig, lane, stage)
    ax.set_xlim(float(current["edges"][0]), float(current["edges"][-1]))
    ymax = float(np.nanmax(current["values"] + current["errors"]))
    ax.set_ylim(0.0, 1.20 * ymax)
    ax.set_ylabel("Unit-area fraction")
    ax.set_xlabel(variable_label)
    ax.legend(loc="upper right", fontsize=9.7)
    ax.minorticks_on()
    # Keep the pending-reference disclosure outside the data pad so it cannot
    # cover a feature peak for any of the eight variables.
    fig.text(
        0.98,
        0.035,
        "PPG12 numerical reference pending",
        ha="right",
        va="bottom",
        fontsize=9.7,
        color="#9a6700",
    )
    fig.subplots_adjust(left=0.125, right=0.98, top=0.875, bottom=0.19)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    write_current_only_csv(out_csv, current)
    return {
        "display_mode": "current_only_reference_pending",
        "normalization": "Current projection is normalized to unit area over the displayed range.",
        "reference_status": "PPG12 numerical reference pending",
        "reference_gap": (
            "No local machine-readable PPG12 signal/inclusive projection was found for this "
            "variable and selection stage. The current histogram exists and is rendered; "
            "reference recovery is the missing input."
        ),
    }


def bdt_reference_payloads() -> dict[str, dict[str, Any]]:
    fig13_19 = load_json(BDT_FIG13_19_JSON)
    fig20 = load_json(BDT_FIG20_JSON)
    return {
        "before_ncb": {
            "source_path": BDT_FIG13_19_JSON,
            "source_hash": sha256_file(BDT_FIG13_19_JSON),
            "panel": fig13_19["panels"]["bdt_no_npb_22_28"],
            "historical_source_key": "bdt_no_npb_22_28",
        },
        "after_ncb": {
            "source_path": BDT_FIG13_19_JSON,
            "source_hash": sha256_file(BDT_FIG13_19_JSON),
            "panel": fig13_19["panels"]["bdt_with_npb_18_22"],
            "historical_source_key": "bdt_with_npb_18_22",
        },
        "tight_id": {
            "source_path": BDT_FIG20_JSON,
            "source_hash": sha256_file(BDT_FIG20_JSON),
            "panel": fig20["panel"],
            "historical_source_key": "bdt_tight_10_14",
        },
    }


def shower_reference_payloads() -> dict[tuple[str, str], dict[str, Any]]:
    e11 = load_json(E11_FIG13_JSON)
    weta = load_json(WETA_FIG13_JSON)
    return {
        ("e11_to_e33", "before_ncb"): {
            "source_path": E11_FIG13_JSON,
            "source_hash": sha256_file(E11_FIG13_JSON),
            "payload": e11,
            "hist": e11.get("histname"),
            "source_roots": e11.get("source_paths"),
        },
        ("weta_cogx", "before_ncb"): {
            "source_path": WETA_FIG13_JSON,
            "source_hash": sha256_file(WETA_FIG13_JSON),
            "payload": weta,
            "hist": weta.get("histname"),
            "source_roots": weta.get("source_paths"),
        },
    }


def base_asset_manifest(
    lane: LaneSpec,
    stage: StageSpec,
    pointer: dict[str, Any],
    hist_key: str,
    current: dict[str, Any],
    png: Path,
    csv_path: Path,
) -> dict[str, Any]:
    return {
        "lane": lane.key,
        "lane_label": lane.label,
        "stage": stage.key,
        "selection": stage.selection,
        "current_pointer": pointer["path"],
        "current_pointer_sha256": pointer["sha256"],
        "current_entry_id": pointer["payload"].get("current_entry_id"),
        "current_campaign_tag": pointer["payload"].get("campaign_tag"),
        "current_root": pointer["root_path"],
        "current_root_sha256": pointer["root_sha256"],
        "current_histogram": hist_key,
        "projection": current["projection"],
        "source_x_axis": current["source_x_axis"],
        "source_y_axis": current["source_y_axis"],
        "source_integral": current["source_integral"],
        "target_integral": current["target_integral"],
        "entries": current["entries"],
        "png": str(png.resolve()),
        "png_sha256": sha256_file(png),
        "csv": str(csv_path.resolve()),
        "csv_sha256": sha256_file(csv_path),
    }


def render_bdt_family(
    outdir: Path,
    pointers: dict[str, dict[str, Any]],
    root_files: dict[str, uproot.ReadOnlyDirectory],
) -> list[dict[str, Any]]:
    references = bdt_reference_payloads()
    outputs: list[dict[str, Any]] = []
    for stage in STAGES:
        reference_meta = references[stage.key]
        panel = reference_meta["panel"]
        hist_key = f"SIM/h2d_bdt_eta0_{stage.pt_key}_{stage.cut_key}"
        for lane in LANES:
            reference = reference_curve(panel["curves"][lane.ppg12_curve])
            current = normalized_current_curve(root_files[lane.key], hist_key, reference["edges"])
            family_dir = outdir / "bdt_score" / stage.key
            png = family_dir / f"bdt_{stage.key}_{lane.key}_ppg12_vs_current_overlay_ratio.png"
            csv_path = family_dir / f"bdt_{stage.key}_{lane.key}_ppg12_vs_current_bins.csv"
            display = draw_overlay_ratio(
                png, csv_path, "BDT score", lane, stage, reference, current
            )
            source_roots = panel.get("summary", {}).get("source_files", {})
            output = base_asset_manifest(
                lane, stage, pointers[lane.key], hist_key, current, png, csv_path
            )
            output.update(
                {
                    "family": "bdt_score",
                    "variable": "bdt",
                    "ppg12_reference_json": str(reference_meta["source_path"].resolve()),
                    "ppg12_reference_json_sha256": reference_meta["source_hash"],
                    "ppg12_reference_curve": lane.ppg12_curve,
                    "ppg12_reference_histogram": (
                        panel.get("hist_name")
                        or panel.get("summary", {}).get("hist")
                        or f"h2d_bdt_eta0_{stage.pt_key}_{stage.cut_key}"
                    ),
                    "ppg12_source_root": source_roots.get(
                        "signal" if lane.key == "photonjet_signal" else "inclusive"
                    ),
                    "historical_source_key": reference_meta["historical_source_key"],
                    "terminology_note": (
                        "Historical file/key names use NPB; reader-facing labels use NCB "
                        "for the non-collision-background preselection."
                    ),
                    **display,
                }
            )
            outputs.append(output)
    return outputs


def render_shower_family(
    outdir: Path,
    pointers: dict[str, dict[str, Any]],
    root_files: dict[str, uproot.ReadOnlyDirectory],
) -> list[dict[str, Any]]:
    references = shower_reference_payloads()
    outputs: list[dict[str, Any]] = []
    for variable in SHOWER_VARIABLES:
        for stage in STAGES:
            reference_meta = references.get((variable.key, stage.key))
            hist_key = f"SIM/h2d_{variable.key}_eta0_{stage.pt_key}_{stage.cut_key}"
            for lane in LANES:
                family_dir = outdir / "shower_shapes" / variable.key / stage.key
                if reference_meta is not None:
                    curve_key = "signal" if lane.key == "photonjet_signal" else "inclusive"
                    reference = reference_curve(reference_meta["payload"][curve_key])
                    current = normalized_current_curve(
                        root_files[lane.key], hist_key, reference["edges"]
                    )
                    png = (
                        family_dir
                        / f"{variable.key}_{stage.key}_{lane.key}_ppg12_vs_current_overlay_ratio.png"
                    )
                    csv_path = (
                        family_dir
                        / f"{variable.key}_{stage.key}_{lane.key}_ppg12_vs_current_bins.csv"
                    )
                    display = draw_overlay_ratio(
                        png,
                        csv_path,
                        variable.label,
                        lane,
                        stage,
                        reference,
                        current,
                    )
                else:
                    target_edges = np.asarray(variable.current_edges, dtype=float)
                    current = normalized_current_curve(
                        root_files[lane.key], hist_key, target_edges
                    )
                    png = (
                        family_dir
                        / f"{variable.key}_{stage.key}_{lane.key}_current_reference_pending.png"
                    )
                    csv_path = (
                        family_dir
                        / f"{variable.key}_{stage.key}_{lane.key}_current_bins.csv"
                    )
                    display = draw_current_only(
                        png,
                        csv_path,
                        variable.label,
                        lane,
                        stage,
                        current,
                    )
                output = base_asset_manifest(
                    lane, stage, pointers[lane.key], hist_key, current, png, csv_path
                )
                output.update(
                    {
                        "family": "shower_shape",
                        "variable": variable.key,
                        "variable_label": variable.label,
                        "terminology_note": (
                            "Historical source naming may use NPB; reader-facing labels use NCB."
                        ),
                        **display,
                    }
                )
                if reference_meta is not None:
                    output.update(
                        {
                            "ppg12_reference_json": str(
                                reference_meta["source_path"].resolve()
                            ),
                            "ppg12_reference_json_sha256": reference_meta["source_hash"],
                            "ppg12_reference_curve": (
                                "signal" if lane.key == "photonjet_signal" else "inclusive"
                            ),
                            "ppg12_reference_histogram": reference_meta["hist"],
                            "ppg12_source_roots": reference_meta["source_roots"],
                        }
                    )
                outputs.append(output)
    return outputs


def _summary_reference(
    variable: VariableSpec | None,
    stage: StageSpec,
    lane: LaneSpec,
    bdt_references: dict[str, dict[str, Any]],
    shower_references: dict[tuple[str, str], dict[str, Any]],
) -> tuple[dict[str, np.ndarray] | None, dict[str, Any] | None]:
    """Return a local PPG12 curve and compact provenance for a summary cell."""
    if variable is None:
        reference_meta = bdt_references[stage.key]
        panel = reference_meta["panel"]
        source_roots = panel.get("summary", {}).get("source_files", {})
        return reference_curve(panel["curves"][lane.ppg12_curve]), {
            "ppg12_reference_json": str(reference_meta["source_path"].resolve()),
            "ppg12_reference_json_sha256": reference_meta["source_hash"],
            "ppg12_reference_curve": lane.ppg12_curve,
            "ppg12_reference_histogram": (
                panel.get("hist_name")
                or panel.get("summary", {}).get("hist")
                or f"h2d_bdt_eta0_{stage.pt_key}_{stage.cut_key}"
            ),
            "ppg12_source_root": source_roots.get(
                "signal" if lane.key == "photonjet_signal" else "inclusive"
            ),
        }

    reference_meta = shower_references.get((variable.key, stage.key))
    if reference_meta is None:
        return None, None
    curve_key = "signal" if lane.key == "photonjet_signal" else "inclusive"
    return reference_curve(reference_meta["payload"][curve_key]), {
        "ppg12_reference_json": str(reference_meta["source_path"].resolve()),
        "ppg12_reference_json_sha256": reference_meta["source_hash"],
        "ppg12_reference_curve": curve_key,
        "ppg12_reference_histogram": reference_meta["hist"],
        "ppg12_source_roots": reference_meta["source_roots"],
    }


def _draw_summary_main(
    ax: plt.Axes,
    lane: LaneSpec,
    reference: dict[str, np.ndarray] | None,
    current: dict[str, Any],
) -> None:
    if reference is not None:
        ax.errorbar(
            reference["centers"],
            reference["values"],
            yerr=reference["errors"],
            fmt=lane.marker,
            color=lane.color,
            markerfacecolor="white",
            markeredgewidth=0.85,
            markersize=3.1,
            linewidth=0.55,
            zorder=2,
        )
    ax.errorbar(
        current["centers"],
        current["values"],
        yerr=current["errors"],
        fmt=lane.marker,
        color=lane.color,
        markerfacecolor=lane.color,
        markeredgewidth=0.55,
        markersize=2.6,
        linewidth=0.55,
        zorder=3,
    )
    ymax = float(np.nanmax(current["values"] + current["errors"]))
    if reference is not None:
        ymax = max(
            ymax,
            float(np.nanmax(reference["values"] + reference["errors"])),
        )
    ax.set_ylim(0.0, 1.18 * ymax)
    ax.tick_params(labelsize=7.8, labelbottom=False)
    ax.minorticks_on()


def _draw_summary_ratio(
    ax: plt.Axes,
    lane: LaneSpec,
    reference: dict[str, np.ndarray],
    current: dict[str, Any],
) -> dict[str, Any]:
    ratio, ratio_error = ratio_with_error(
        current["values"],
        current["errors"],
        reference["values"],
        reference["errors"],
    )
    stable = (reference["values"] > 0.002) & (current["values"] > 0.002)
    ratio_scale, ratio_limits = choose_ratio_axis(ratio, stable)
    finite = np.isfinite(ratio) & np.isfinite(ratio_error)
    ax.errorbar(
        reference["centers"][finite],
        ratio[finite],
        yerr=ratio_error[finite],
        fmt=lane.marker,
        color="black",
        markerfacecolor="black",
        markersize=2.0,
        linewidth=0.45,
    )
    ax.axhline(1.0, color="0.45", linestyle="--", linewidth=0.75)
    ax.set_yscale(ratio_scale)
    ax.set_ylim(*ratio_limits)
    ax.tick_params(labelsize=7.2)
    ax.minorticks_on()
    stable_ratio = ratio[stable & np.isfinite(ratio)]
    return {
        "display_mode": "ppg12_current_overlay_ratio",
        "ratio_scale": ratio_scale,
        "stable_ratio_bin_count": int(stable_ratio.size),
        "stable_ratio_mean": float(np.mean(stable_ratio)) if stable_ratio.size else None,
        "stable_ratio_max_abs_minus_one": (
            float(np.max(np.abs(stable_ratio - 1.0))) if stable_ratio.size else None
        ),
    }


def _draw_summary_pending(ax: plt.Axes) -> dict[str, Any]:
    ax.set_facecolor("#fff7df")
    ax.text(
        0.5,
        0.52,
        "PPG12 numerical reference pending",
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=6.8,
        color="#8a5a00",
        fontweight="bold",
    )
    ax.set_ylim(0.0, 1.0)
    ax.set_yticks([])
    ax.tick_params(axis="x", labelsize=7.2)
    for side in ("left", "right", "top"):
        ax.spines[side].set_visible(False)
    return {
        "display_mode": "current_only_reference_pending",
        "reference_status": "PPG12 numerical reference pending",
        "reference_gap": (
            "The current histogram exists and is rendered. Only the local "
            "machine-readable PPG12 numerical projection is unavailable."
        ),
    }


def render_summary_canvas(
    outdir: Path,
    pointers: dict[str, dict[str, Any]],
    root_files: dict[str, uproot.ReadOnlyDirectory],
    variable: VariableSpec | None,
) -> dict[str, Any]:
    """Render one manuscript-ready two-lane by three-stage audit canvas."""
    bdt_references = bdt_reference_payloads()
    shower_references = shower_reference_payloads()
    variable_key = "bdt" if variable is None else variable.key
    variable_label = "BDT score" if variable is None else variable.label
    hist_variable = "bdt" if variable is None else variable.key

    if variable is None:
        family_dir = outdir / "bdt_score"
        png = family_dir / "bdt_lane_by_stage_2x3_ppg12_vs_current.png"
        title = "BDT-score parity across the PPG12 selection sequence"
    else:
        family_dir = outdir / "shower_shapes" / variable.key
        png = family_dir / f"{variable.key}_lane_by_stage_2x3_current_sim_parity.png"
        title = f"{variable.label} across the PPG12 selection sequence"
    family_dir.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=(15.0, 10.4), dpi=160)
    grid = fig.add_gridspec(
        4,
        3,
        height_ratios=(3.0, 1.0, 3.0, 1.0),
        left=0.075,
        right=0.985,
        bottom=0.085,
        top=0.815,
        wspace=0.24,
        hspace=0.055,
    )
    cells: list[dict[str, Any]] = []
    reference_cell_count = 0
    pending_cell_count = 0
    main_axes: list[list[plt.Axes]] = [[], []]

    for lane_index, lane in enumerate(LANES):
        main_row = 2 * lane_index
        ratio_row = main_row + 1
        for stage_index, stage in enumerate(STAGES):
            ax = fig.add_subplot(grid[main_row, stage_index])
            ratio_ax = fig.add_subplot(grid[ratio_row, stage_index], sharex=ax)
            main_axes[lane_index].append(ax)
            reference, reference_provenance = _summary_reference(
                variable,
                stage,
                lane,
                bdt_references,
                shower_references,
            )
            hist_key = f"SIM/h2d_{hist_variable}_eta0_{stage.pt_key}_{stage.cut_key}"
            if reference is not None:
                target_edges = reference["edges"]
            else:
                if variable is None:
                    raise RuntimeError(f"Missing mandatory BDT reference for {stage.key}")
                target_edges = np.asarray(variable.current_edges, dtype=float)
            current = normalized_current_curve(
                root_files[lane.key], hist_key, np.asarray(target_edges, dtype=float)
            )
            _draw_summary_main(ax, lane, reference, current)
            if reference is not None:
                display = _draw_summary_ratio(ratio_ax, lane, reference, current)
                reference_cell_count += 1
            else:
                display = _draw_summary_pending(ratio_ax)
                pending_cell_count += 1

            x_low = float(target_edges[0])
            x_high = float(target_edges[-1])
            if variable is None and stage.key == "tight_id":
                x_low, x_high = 0.50, 1.02
            ax.set_xlim(x_low, x_high)
            ratio_ax.set_xlim(x_low, x_high)
            if stage_index == 0:
                ax.set_ylabel("Unit-area fraction", fontsize=9.0)
                if reference is not None:
                    ratio_ax.set_ylabel("Cur./PPG12", fontsize=8.0)
            if lane_index == 0:
                stage_title = {
                    "before_ncb": "Before NCB preselection\n$22<E_T^\\gamma<28$ GeV",
                    "after_ncb": "After NCB preselection\n$18<E_T^\\gamma<22$ GeV",
                    "tight_id": "Tight photon-ID selection\n$10<E_T^\\gamma<14$ GeV",
                }[stage.key]
                ax.set_title(stage_title, fontsize=10.7, pad=7.0)

            cell = {
                "lane": lane.key,
                "lane_label": lane.label,
                "stage": stage.key,
                "selection": stage.selection,
                "current_pointer": pointers[lane.key]["path"],
                "current_pointer_sha256": pointers[lane.key]["sha256"],
                "current_root": pointers[lane.key]["root_path"],
                "current_root_sha256": pointers[lane.key]["root_sha256"],
                "current_histogram": hist_key,
                "projection": current["projection"],
                "source_integral": current["source_integral"],
                "target_integral": current["target_integral"],
                **display,
            }
            if reference_provenance is not None:
                cell.update(reference_provenance)
            cells.append(cell)

    fig.text(
        0.075,
        0.985,
        "sPHENIX",
        ha="left",
        va="top",
        fontsize=15,
        fontweight="bold",
        fontstyle="italic",
    )
    fig.text(0.155, 0.985, "Internal", ha="left", va="top", fontsize=15)
    fig.text(
        0.985,
        0.985,
        r"$p{+}p\ \sqrt{s}=200$ GeV, $|\eta^\gamma|<0.7$",
        ha="right",
        va="top",
        fontsize=10.5,
    )
    fig.suptitle(title, y=0.948, fontsize=16.0, fontweight="bold")
    fig.text(
        0.018,
        0.655,
        LANES[0].label,
        rotation=90,
        ha="center",
        va="center",
        fontsize=10.0,
        fontweight="bold",
        color=LANES[0].color,
    )
    fig.text(
        0.018,
        0.285,
        LANES[1].label,
        rotation=90,
        ha="center",
        va="center",
        fontsize=10.0,
        fontweight="bold",
        color=LANES[1].color,
    )
    legend_handles: list[Any] = []
    if reference_cell_count:
        legend_handles.append(
            Line2D(
                [0],
                [0],
                marker="o",
                color="0.2",
                markerfacecolor="white",
                linestyle="none",
                label="PPG12 SDCC numerical reference (open)",
            )
        )
    legend_handles.append(
        Line2D(
            [0],
            [0],
            marker="o",
            color="0.2",
            markerfacecolor="0.2",
            linestyle="none",
            label="July 16 registered current output (filled)",
        )
    )
    if pending_cell_count:
        legend_handles.append(
            Patch(
                facecolor="#fff7df",
                edgecolor="#c28a20",
                label="PPG12 numerical reference pending",
            )
        )
    fig.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.53, 0.912),
        ncol=len(legend_handles),
        fontsize=9.2,
        frameon=False,
        columnspacing=1.7,
    )
    fig.supxlabel(variable_label, y=0.025, fontsize=11.0)
    fig.savefig(png, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    return {
        "family": "bdt_summary_canvas" if variable is None else "shower_shape_summary_canvas",
        "variable": variable_key,
        "variable_label": variable_label,
        "layout": "two simulation lanes by three selection stages; each referenced cell has an overlay and ratio pad",
        "normalization": "Each current and PPG12 projection is independently normalized to unit area over the displayed cell range.",
        "reference_cell_count": reference_cell_count,
        "reference_pending_cell_count": pending_cell_count,
        "current_histogram_missing_count": 0,
        "png": str(png.resolve()),
        "png_sha256": sha256_file(png),
        "cells": cells,
    }


def render_summary_canvases(
    outdir: Path,
    pointers: dict[str, dict[str, Any]],
    root_files: dict[str, uproot.ReadOnlyDirectory],
    families: Iterable[str],
) -> list[dict[str, Any]]:
    summaries: list[dict[str, Any]] = []
    requested = set(families)
    if "bdt" in requested:
        summaries.append(render_summary_canvas(outdir, pointers, root_files, None))
    if "shower" in requested:
        summaries.extend(
            render_summary_canvas(outdir, pointers, root_files, variable)
            for variable in SHOWER_VARIABLES
        )
    return summaries


def load_abcd_reference() -> tuple[list[dict[str, str]], dict[str, Any]]:
    if not ABCD_REFERENCE_CSV.exists() or not ABCD_REFERENCE_MANIFEST.exists():
        raise FileNotFoundError(
            "Local machine-readable PPG12 signal-ABCD reference is unavailable: "
            f"{ABCD_REFERENCE_CSV} / {ABCD_REFERENCE_MANIFEST}"
        )
    with ABCD_REFERENCE_CSV.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    manifest = load_json(ABCD_REFERENCE_MANIFEST)
    if len(rows) != 11:
        raise RuntimeError(f"Expected 11 signal-ABCD reference bins, found {len(rows)}")
    return rows, manifest


def render_abcd_family(
    outdir: Path,
    pointer: dict[str, Any],
    root_file: uproot.ReadOnlyDirectory,
) -> list[dict[str, Any]]:
    reference_rows, reference_manifest = load_abcd_reference()
    edges = np.asarray(
        [float(reference_rows[0]["pt_lo"])]
        + [float(row["pt_hi"]) for row in reference_rows],
        dtype=float,
    )
    centers = 0.5 * (edges[:-1] + edges[1:])
    half_width = 0.5 * (edges[1:] - edges[:-1])
    curves: dict[str, dict[str, Any]] = {}
    for region, label, color, marker, hist_key in ABCD_REGIONS:
        projected = project_x(root_file, hist_key)
        values, variances = rebin_overlap(
            projected["edges"],
            projected["values"],
            projected["variances"],
            edges,
        )
        ppg12 = np.asarray([float(row[f"ppg12_{region}"]) for row in reference_rows])
        ppg12_error = np.asarray(
            [float(row[f"ppg12_{region}_err"]) for row in reference_rows]
        )
        current_error = np.sqrt(np.maximum(variances, 0.0))
        ratio, ratio_error = ratio_with_error(values, current_error, ppg12, ppg12_error)
        curves[region] = {
            "label": label,
            "color": color,
            "marker": marker,
            "hist_key": hist_key,
            "ppg12": ppg12,
            "ppg12_error": ppg12_error,
            "current": values,
            "current_error": current_error,
            "ratio": ratio,
            "ratio_error": ratio_error,
            "projected": projected,
        }

    family_dir = outdir / "raw_signal_abcd"
    family_dir.mkdir(parents=True, exist_ok=True)
    png = family_dir / "signal_abcd_ppg12_vs_current_photonjet_overlay_ratio.png"
    csv_path = family_dir / "signal_abcd_ppg12_vs_current_photonjet_bins.csv"
    fig, (ax, ratio_ax) = plt.subplots(
        2,
        1,
        figsize=(8.8, 7.7),
        dpi=160,
        sharex=True,
        gridspec_kw={"height_ratios": [3.1, 1.0], "hspace": 0.04},
    )
    for region, payload in curves.items():
        ax.errorbar(
            centers - 0.10,
            payload["ppg12"],
            xerr=half_width,
            yerr=payload["ppg12_error"],
            fmt=payload["marker"],
            color=payload["color"],
            markerfacecolor="white",
            markeredgewidth=1.0,
            markersize=4.8,
            linewidth=0.7,
            label=f"PPG12 {region}: {payload['label']}",
        )
        ax.errorbar(
            centers + 0.10,
            payload["current"],
            xerr=half_width,
            yerr=payload["current_error"],
            fmt=payload["marker"],
            color=payload["color"],
            markerfacecolor=payload["color"],
            markeredgewidth=0.7,
            markersize=4.2,
            linewidth=0.65,
            label=f"Current {region}: {payload['label']}",
        )
        finite = np.isfinite(payload["ratio"])
        ratio_ax.errorbar(
            centers[finite],
            payload["ratio"][finite],
            xerr=half_width[finite],
            yerr=payload["ratio_error"][finite],
            fmt=payload["marker"],
            color=payload["color"],
            markerfacecolor=payload["color"],
            markersize=3.5,
            linewidth=0.65,
            label=region,
        )
    positive = np.concatenate(
        [
            payload[key][payload[key] > 0]
            for payload in curves.values()
            for key in ("ppg12", "current")
        ]
    )
    ax.set_yscale("log")
    ax.set_ylim(float(np.min(positive)) / 3.0, float(np.max(positive)) * 4.0)
    fig.text(0.13, 0.978, "sPHENIX", va="top", ha="left", fontsize=14,
             fontweight="bold", fontstyle="italic")
    fig.text(0.285, 0.978, "Internal", va="top", ha="left", fontsize=14)
    fig.text(0.13, 0.943, r"$p{+}p\ \sqrt{s}=200$ GeV, $|\eta^\gamma|<0.7$",
             va="top", ha="left", fontsize=10.5)
    fig.text(0.985, 0.943,
             "Signal-MC ABCD correction inputs; absolute event-preweighted counts",
             va="top", ha="right", fontsize=9.7)
    ax.set_ylabel("Signal count")
    ax.legend(loc="upper right", ncol=2, fontsize=8.0, columnspacing=0.8, handletextpad=0.4)
    ax.tick_params(labelbottom=False)
    ratio_ax.axhline(1.0, color="0.4", linestyle="--", linewidth=0.9)
    all_ratio = np.concatenate(
        [payload["ratio"][np.isfinite(payload["ratio"])] for payload in curves.values()]
    )
    ratio_ax.set_ylim(max(0.0, 0.85 * float(np.min(all_ratio))), 1.15 * float(np.max(all_ratio)))
    ratio_ax.set_ylabel("Current / PPG12")
    ratio_ax.set_xlabel(r"$E_T^{\gamma,\mathrm{rec}}$ [GeV]")
    ratio_ax.set_xlim(float(edges[0]), float(edges[-1]))
    ratio_ax.minorticks_on()
    fig.subplots_adjust(left=0.13, right=0.985, top=0.885, bottom=0.10, hspace=0.04)
    fig.savefig(png, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)

    with csv_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "region",
                "pt_low",
                "pt_high",
                "ppg12",
                "ppg12_error",
                "current",
                "current_error",
                "current_over_ppg12",
                "current_over_ppg12_error",
            ]
        )
        for region, payload in curves.items():
            for index in range(len(centers)):
                writer.writerow(
                    [
                        region,
                        f"{edges[index]:.12g}",
                        f"{edges[index + 1]:.12g}",
                        f"{payload['ppg12'][index]:.12g}",
                        f"{payload['ppg12_error'][index]:.12g}",
                        f"{payload['current'][index]:.12g}",
                        f"{payload['current_error'][index]:.12g}",
                        f"{payload['ratio'][index]:.12g}",
                        f"{payload['ratio_error'][index]:.12g}",
                    ]
                )

    old_current_sha = reference_manifest.get("current_root_sha256")
    if old_current_sha and old_current_sha != pointer["root_sha256"]:
        raise RuntimeError(
            "The registered current photon ROOT differs from the ROOT used to extract the "
            "stored reference CSV. Regenerate the extraction before using it."
        )
    output = {
        "family": "raw_signal_abcd",
        "display_mode": "absolute_ppg12_current_overlay_ratio",
        "selection": r"10 <= E_T^gamma < 36 GeV, |eta^gamma| < 0.7",
        "normalization": "No normalization, fitted scale, or shape rescaling is applied.",
        "ratio_definition": "Current / PPG12",
        "current_pointer": pointer["path"],
        "current_pointer_sha256": pointer["sha256"],
        "current_entry_id": pointer["payload"].get("current_entry_id"),
        "current_root": pointer["root_path"],
        "current_root_sha256": pointer["root_sha256"],
        "ppg12_reference_csv": str(ABCD_REFERENCE_CSV.resolve()),
        "ppg12_reference_csv_sha256": sha256_file(ABCD_REFERENCE_CSV),
        "ppg12_reference_manifest": str(ABCD_REFERENCE_MANIFEST.resolve()),
        "ppg12_reference_manifest_sha256": sha256_file(ABCD_REFERENCE_MANIFEST),
        "ppg12_source_root": reference_manifest.get("ppg12_root"),
        "ppg12_source_root_sha256": reference_manifest.get("ppg12_root_sha256"),
        "histograms": {
            region: {
                "ppg12": reference_manifest.get("objects", {}).get("ppg12", {}).get(region),
                "current": payload["hist_key"],
                "current_source_integral": payload["projected"]["source_integral"],
                "current_entries": payload["projected"]["entries"],
            }
            for region, payload in curves.items()
        },
        "registered_root_matches_prior_extraction_sha256": bool(
            old_current_sha == pointer["root_sha256"]
        ),
        "png": str(png.resolve()),
        "png_sha256": sha256_file(png),
        "csv": str(csv_path.resolve()),
        "csv_sha256": sha256_file(csv_path),
    }
    return [output]


def verify_assets(outputs: Iterable[dict[str, Any]]) -> None:
    for output in outputs:
        for kind in ("png", "csv"):
            path = Path(output[kind])
            if not path.exists() or path.stat().st_size <= 0:
                raise RuntimeError(f"Missing or empty generated {kind}: {path}")
            expected = output[f"{kind}_sha256"]
            actual = sha256_file(path)
            if actual != expected:
                raise RuntimeError(f"Generated {kind} hash changed unexpectedly: {path}")


def verify_summary_canvases(outputs: Iterable[dict[str, Any]]) -> None:
    for output in outputs:
        path = Path(output["png"])
        if not path.exists() or path.stat().st_size <= 0:
            raise RuntimeError(f"Missing or empty generated summary canvas: {path}")
        if sha256_file(path) != output["png_sha256"]:
            raise RuntimeError(f"Generated summary canvas hash changed unexpectedly: {path}")


def main() -> int:
    args = parse_args()
    configure_style()
    args.outdir.mkdir(parents=True, exist_ok=True)

    pointers = {
        lane.key: load_pointer(lane.pointer)
        for lane in LANES
    }
    root_files = {
        lane.key: uproot.open(pointers[lane.key]["root_path"])
        for lane in LANES
    }
    outputs: list[dict[str, Any]] = []
    summary_canvases: list[dict[str, Any]] = []
    try:
        if "bdt" in args.families:
            outputs.extend(render_bdt_family(args.outdir, pointers, root_files))
        if "shower" in args.families:
            outputs.extend(render_shower_family(args.outdir, pointers, root_files))
        if "abcd" in args.families:
            outputs.extend(
                render_abcd_family(
                    args.outdir,
                    pointers["photonjet_signal"],
                    root_files["photonjet_signal"],
                )
            )
        summary_canvases.extend(
            render_summary_canvases(args.outdir, pointers, root_files, args.families)
        )
    finally:
        for root_file in root_files.values():
            root_file.close()

    verify_assets(outputs)
    verify_summary_canvases(summary_canvases)
    overlay_count = sum(
        output["display_mode"] in {
            "ppg12_current_overlay_ratio",
            "absolute_ppg12_current_overlay_ratio",
        }
        for output in outputs
    )
    pending_count = sum(
        output["display_mode"] == "current_only_reference_pending"
        for output in outputs
    )
    aggregate = {
        "schema": "ppg12_current_sim_ian_refresh_v1",
        "script": str(Path(__file__).resolve()),
        "script_sha256": sha256_file(Path(__file__).resolve()),
        "purpose": (
            "Current-pointer-resolved pp photon+jet and inclusive-jet SIM figures for "
            "the Stage-1 PPG12-parity IAN audit."
        ),
        "scientific_status": (
            "Current corrected-SI matched-pair diagnostic; current registration does not "
            "by itself establish historical PPG12 bin-by-bin identity."
        ),
        "families_requested": list(args.families),
        "terminology": (
            "NCB means non-collision background. Historical PPG12 paths and object-family "
            "labels may retain the older NPB spelling."
        ),
        "pointers": pointers,
        "asset_count": len(outputs),
        "summary_canvas_count": len(summary_canvases),
        "total_png_count": len(outputs) + len(summary_canvases),
        "overlay_ratio_asset_count": overlay_count,
        "current_only_reference_pending_asset_count": pending_count,
        "current_histogram_missing_count": 0,
        "reference_coverage_summary": (
            "All BDT stages have local numerical PPG12 signal and inclusive references. "
            "Among shower-shape panels, before-NCB weta_cogx and e11_to_e33 have local "
            "numerical PPG12 references; the remaining panels render the complete current "
            "histograms and explicitly mark reference recovery as pending."
        ),
        "summary_canvases": summary_canvases,
        "outputs": outputs,
    }
    manifest_path = args.outdir / "ppg12_current_sim_ian_refresh_manifest.json"
    manifest_path.write_text(json.dumps(aggregate, indent=2, sort_keys=True) + "\n")
    print(manifest_path)
    print(
        json.dumps(
            {
                "asset_count": len(outputs),
                "summary_canvas_count": len(summary_canvases),
                "total_png_count": len(outputs) + len(summary_canvases),
                "overlay_ratio_asset_count": overlay_count,
                "current_only_reference_pending_asset_count": pending_count,
                "current_histogram_missing_count": 0,
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
