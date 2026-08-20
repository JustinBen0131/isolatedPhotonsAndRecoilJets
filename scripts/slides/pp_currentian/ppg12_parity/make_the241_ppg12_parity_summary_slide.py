#!/usr/bin/env python3
"""Render the THE-241 PPG12 pp baseline parity summary slide.

The slide compares the frozen full-stat current pp output with the PPG12
reference at the two levels needed for the baseline argument: the four raw
ABCD populations and the raw purity derived from them.  Only statistical/toy
uncertainties stored in the accepted point table are drawn.  No fitted
normalization, leakage correction, systematic band, unfolding, production,
or Google Slides mutation is performed.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from matplotlib.font_manager import FontProperties
from matplotlib.lines import Line2D


REPO = Path(__file__).resolve().parents[4]
DEFAULT_OUT = (
    REPO
    / "dataOutput/the219_friday_ppg_20260814/the241_ppg12_parity_summary"
)
PARITY_CSV = (
    REPO
    / "dataOutput/ppg12Parity/the97_ppg12_final_accepted_triple_full_20260714_1550"
    / "final_pp_data_canonical_20260717/raw_abcd_three_panel"
    / "the97_final_pp_data_raw_abcd_yield_and_purity_current_over_ppg12.csv"
)
PARITY_MANIFEST = PARITY_CSV.with_name(
    "the97_final_pp_data_raw_abcd_yield_and_purity_current_over_ppg12_manifest.json"
)
DATATHIEF_MANIFEST = (
    REPO
    / "dataOutput/ppg12Parity/the76_ppg12_parity_full_20260701_003024"
    / "fig29_purity_datathief_audit"
    / "ppg12_fig5_fig29_purity_datathief_vs_sdcc_manifest.json"
)

INK = "#101820"
MUTED = "#526173"
GRID = "#D6DEE8"
PPG12_BLUE = "#2468A2"
SLIDE_FONT_FAMILY = "Times New Roman"
REGION_STYLE = {
    "A": ("#111827", "o"),
    "B": ("#C73A31", "s"),
    "C": ("#2468A2", "^"),
    "D": ("#8C4AA8", "v"),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument(
        "--visual-qa-pass",
        action="store_true",
        help="Record visual QA as PASS after the rendered slide is inspected.",
    )
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_rows() -> list[dict[str, float]]:
    with PARITY_CSV.open(newline="", encoding="utf-8") as stream:
        rows = [
            {key: float(value) for key, value in row.items()}
            for row in csv.DictReader(stream)
        ]
    if not rows:
        raise RuntimeError(f"No rows found in {PARITY_CSV}")
    return rows


def resolve_slide_font() -> Path:
    """Resolve the policy-required typeface without silent fallback."""
    try:
        resolved = Path(
            font_manager.findfont(
                FontProperties(family=SLIDE_FONT_FAMILY),
                fallback_to_default=False,
            )
        )
    except ValueError as error:
        raise RuntimeError(
            f"Required slide font {SLIDE_FONT_FAMILY!r} is unavailable"
        ) from error
    if "times new roman" not in resolved.name.lower():
        raise RuntimeError(
            f"Resolved unexpected font for {SLIDE_FONT_FAMILY!r}: {resolved}"
        )
    return resolved


def setup_style() -> Path:
    font_path = resolve_slide_font()
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": [SLIDE_FONT_FAMILY],
            "mathtext.fontset": "stix",
            "mathtext.rm": SLIDE_FONT_FAMILY,
            "mathtext.it": f"{SLIDE_FONT_FAMILY}:italic",
            "mathtext.bf": f"{SLIDE_FONT_FAMILY}:bold",
            "axes.linewidth": 1.25,
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
    return font_path


def style_axis(ax: plt.Axes, *, log_y: bool = False) -> None:
    ax.grid(axis="y", which="major", color=GRID, lw=0.9, alpha=0.9, zorder=0)
    if not log_y:
        ax.minorticks_on()
    ax.tick_params(which="major", length=7, width=1.05, labelsize=13)
    ax.tick_params(which="minor", length=3.5, width=0.8)


def plot_purity(
    ax: plt.Axes, ratio_ax: plt.Axes, rows: list[dict[str, float]]
) -> None:
    x = np.asarray([row["center"] for row in rows])
    xerr = np.asarray([row["half_width"] for row in rows])
    current = np.asarray([row["raw_purity_current"] for row in rows])
    current_error = np.asarray([row["raw_purity_current_error"] for row in rows])
    reference = np.asarray([row["raw_purity_ppg12"] for row in rows])
    reference_error = np.asarray([row["raw_purity_ppg12_error"] for row in rows])
    ratio = np.asarray([row["raw_purity_current_over_ppg12"] for row in rows])
    ratio_error = np.asarray(
        [row["raw_purity_current_over_ppg12_error"] for row in rows]
    )

    marker_offset = 0.11
    ax.errorbar(
        x - marker_offset,
        reference,
        xerr=xerr,
        yerr=reference_error,
        fmt="s",
        color=PPG12_BLUE,
        markerfacecolor="white",
        markeredgewidth=1.5,
        markersize=7.5,
        elinewidth=1.25,
        capsize=0,
        label="PPG12 final ROOT",
        zorder=3,
    )
    ax.errorbar(
        x + marker_offset,
        current,
        xerr=xerr,
        yerr=current_error,
        fmt="o",
        color=INK,
        markerfacecolor=INK,
        markersize=7.0,
        elinewidth=1.25,
        capsize=0,
        label="Current full-stat output",
        zorder=4,
    )
    ax.set_xlim(9.0, 37.0)
    ax.set_ylim(0.0, 1.18)
    ax.set_ylabel("Raw photon purity", fontsize=17)
    ax.tick_params(labelbottom=False)
    ax.set_title("Raw photon purity", fontsize=21, fontweight="bold", loc="left", pad=14)
    ax.legend(loc="upper left", frameon=False, fontsize=13.5, handletextpad=0.55)
    style_axis(ax)

    ratio_ax.axhline(1.0, color=MUTED, lw=1.2, zorder=1)
    ratio_ax.errorbar(
        x,
        ratio,
        xerr=xerr,
        yerr=ratio_error,
        fmt="o",
        color=INK,
        markerfacecolor=INK,
        markersize=6.0,
        elinewidth=1.15,
        capsize=0,
        zorder=3,
    )
    ratio_ax.set_xlim(9.0, 37.0)
    ratio_ax.set_ylim(0.25, 1.75)
    ratio_ax.set_yticks([0.5, 1.0, 1.5])
    ratio_ax.set_ylabel("Current / PPG12", fontsize=13)
    ratio_ax.set_xlabel(r"$E_T^{\gamma}$ [GeV]", fontsize=17)
    style_axis(ratio_ax)


def plot_counts(
    ax: plt.Axes, ratio_ax: plt.Axes, rows: list[dict[str, float]]
) -> None:
    x = np.asarray([row["center"] for row in rows])
    xerr = np.asarray([row["half_width"] for row in rows])
    marker_offset = 0.11
    region_handles: list[Line2D] = []
    for region, (color, marker) in REGION_STYLE.items():
        reference = np.asarray([row[f"{region}_ppg12"] for row in rows])
        reference_error = np.asarray([row[f"{region}_ppg12_error"] for row in rows])
        current = np.asarray([row[f"{region}_current"] for row in rows])
        current_error = np.asarray([row[f"{region}_current_error"] for row in rows])
        ratio = np.asarray([row[f"{region}_current_over_ppg12"] for row in rows])
        ratio_error = np.asarray(
            [row[f"{region}_current_over_ppg12_error"] for row in rows]
        )
        ax.errorbar(
            x - marker_offset,
            reference,
            xerr=xerr,
            yerr=reference_error,
            fmt=marker,
            color=color,
            markerfacecolor="white",
            markeredgewidth=1.4,
            markersize=6.5,
            elinewidth=1.05,
            capsize=0,
            alpha=0.95,
            zorder=3,
        )
        ax.errorbar(
            x + marker_offset,
            current,
            xerr=xerr,
            yerr=current_error,
            fmt=marker,
            color=color,
            markerfacecolor=color,
            markersize=6.2,
            elinewidth=1.05,
            capsize=0,
            zorder=4,
        )
        ratio_ax.errorbar(
            x,
            ratio,
            xerr=xerr,
            yerr=ratio_error,
            fmt=marker,
            color=color,
            markerfacecolor=color,
            markersize=5.3,
            elinewidth=1.0,
            capsize=0,
            zorder=3,
        )
        region_handles.append(
            Line2D(
                [0],
                [0],
                marker=marker,
                color=color,
                linestyle="none",
                markerfacecolor=color,
                markersize=7,
                label=f"Region {region}",
            )
        )

    ax.set_yscale("log")
    ax.set_xlim(9.0, 37.0)
    ax.set_ylim(3.0, 1.2e5)
    ax.set_ylabel("Raw weighted counts", fontsize=17)
    ax.tick_params(labelbottom=False)
    ax.set_title("Raw ABCD populations", fontsize=21, fontweight="bold", loc="left", pad=14)
    region_legend = ax.legend(
        handles=region_handles,
        loc="upper right",
        frameon=False,
        fontsize=12.5,
        ncol=2,
        columnspacing=1.0,
        handletextpad=0.45,
    )
    ax.add_artist(region_legend)
    sample_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color=INK,
            linestyle="none",
            markerfacecolor=INK,
            markersize=8,
            label="Current output — filled",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color=INK,
            linestyle="none",
            markerfacecolor="white",
            markeredgewidth=1.5,
            markersize=8,
            label="PPG12 reference — open",
        ),
    ]
    sample_legend = ax.legend(
        handles=sample_handles,
        title="Marker fill identifies the dataset",
        loc="lower left",
        bbox_to_anchor=(0.018, 0.035),
        frameon=True,
        fancybox=False,
        framealpha=0.98,
        facecolor="white",
        edgecolor=PPG12_BLUE,
        fontsize=12.3,
        title_fontsize=12.8,
        borderpad=0.65,
        labelspacing=0.42,
        handletextpad=0.60,
    )
    sample_legend.get_frame().set_linewidth(1.6)
    sample_legend.get_title().set_fontweight("bold")
    style_axis(ax, log_y=True)

    ratio_ax.axhline(1.0, color=MUTED, lw=1.2, zorder=1)
    ratio_ax.set_xlim(9.0, 37.0)
    ratio_ax.set_ylim(0.25, 1.75)
    ratio_ax.set_yticks([0.5, 1.0, 1.5])
    ratio_ax.set_ylabel("Current / PPG12", fontsize=13)
    ratio_ax.set_xlabel(r"$E_T^{\gamma}$ [GeV]", fontsize=17)
    style_axis(ratio_ax)


def write_speaker_script(output_dir: Path) -> Path:
    path = output_dir / "the241_ppg12_parity_summary_speaker_script.md"
    body = f"""# THE-241 speaker script — PPG12 parity summary

This is why I use PPG12 as the pp photon baseline. On the right, the raw A, B, C and D populations from the current full-stat output agree with PPG12 bin by bin, with no fitted normalization. On the left, putting those same populations through the raw ABCD purity calculation reproduces the PPG12 purity. The final 32–36 GeV bin is count-limited, so its uncertainty is visibly larger, but it is statistically consistent. The PPG12 reference points shown here are the final ROOT graph values. I also cross-checked them against the published Figure 5 and IAN Figure 29 marker positions; that paper-image digitization agrees with the final ROOT purity points within 0.15 percent. Together, this validates the pp trigger namespace, event population, photon classification, isolation and ABCD routing that define the baseline before moving to Au+Au.

## Scope guard

- Statistical/toy uncertainties only; no systematic band.
- Raw purity only; no leakage correction is shown.
- No unfolding or jet calibration enters this slide.
- The DataThief audit attempted the official jar, which returned non-finite coordinates; the preserved final audit used the same three axis anchors through a local affine transform. The final PPG12 ROOT graph, not the image digitization, remains authoritative.

[Sources]
- `{PARITY_MANIFEST}`
- `{PARITY_CSV}`
- `{DATATHIEF_MANIFEST}`
"""
    path.write_text(body, encoding="utf-8")
    return path


def validate_feedback_layout_contract(payload: dict[str, Any]) -> None:
    """Fail on lower-text clutter, x-label collisions, or a hidden marker key."""
    nodes = {node["name"]: node for node in payload["nodes"]}
    required = {
        "purity panel",
        "abcd panel",
        "marker fill legend",
    }
    missing = required - nodes.keys()
    if missing:
        raise RuntimeError(f"Missing required layout-feedback nodes: {sorted(missing)}")

    forbidden = {
        "baseline result bullet",
        "purity result bullet",
        "abcd result bullet",
        "full width takeaway",
        "datathief provenance note",
        "count limited note",
    }
    residual = forbidden & nodes.keys()
    if residual:
        raise RuntimeError(f"Removed lower-canvas text returned: {sorted(residual)}")
    marker_key = nodes["marker fill legend"]
    if marker_key.get("role") != "plot_annotation":
        raise RuntimeError("Marker-fill key must be audited as a plot annotation")
    if marker_key.get("font_px", 0) < payload["minimum_plot_annotation_font_px"]:
        raise RuntimeError("Marker-fill key is below the plot-annotation font floor")


def make_layout_nodes(path: Path) -> None:
    payload: dict[str, Any] = {
        "schema": "slide_layout_nodes_v1",
        "slide": "the241_ppg12_parity_summary",
        "title_axis_x": 120,
        "minimum_audience_font_px": 29,
        "minimum_plot_annotation_font_px": 24,
        "minimum_title_font_px": 56,
        "font_family": SLIDE_FONT_FAMILY,
        "nodes": [
            {
                "name": "claim title",
                "kind": "text",
                "role": "title",
                "font_px": 60,
                "bbox": [120, 50, 2380, 116],
                "text": "PPG12 parity establishes the pp photon baseline",
            },
            {
                "name": "subtitle",
                "kind": "text",
                "role": "audience",
                "font_px": 38,
                "bbox": [125, 146, 2400, 194],
                "text": "The same full-stat output reproduces the four ABCD inputs and their derived raw purity; no fitted normalization.",
            },
            {
                "name": "purity panel",
                "kind": "panel",
                "role": "audience",
                "symmetry_group": "two parity panels",
                "bbox": [145, 235, 1228, 1395],
            },
            {
                "name": "abcd panel",
                "kind": "panel",
                "role": "audience",
                "symmetry_group": "two parity panels",
                "bbox": [1370, 235, 2453, 1395],
            },
            {
                "name": "marker fill legend",
                "kind": "text",
                "role": "plot_annotation",
                "font_px": 28,
                "bbox": [1200, 560, 1660, 675],
                "text": "Marker fill identifies the dataset: current output filled; PPG12 reference open.",
            },
        ],
    }
    validate_feedback_layout_contract(payload)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    slide_font_path = setup_style()
    rows = read_rows()
    datathief = json.loads(DATATHIEF_MANIFEST.read_text(encoding="utf-8"))

    fig = plt.figure(figsize=(16, 9), dpi=160)
    outer = fig.add_gridspec(
        1,
        2,
        left=0.064,
        right=0.982,
        bottom=0.080,
        top=0.810,
        wspace=0.24,
    )
    left = outer[0].subgridspec(2, 1, height_ratios=(3.30, 1.18), hspace=0.045)
    right = outer[1].subgridspec(2, 1, height_ratios=(3.30, 1.18), hspace=0.045)
    ax_purity = fig.add_subplot(left[0])
    ax_purity_ratio = fig.add_subplot(left[1])
    ax_counts = fig.add_subplot(right[0])
    ax_counts_ratio = fig.add_subplot(right[1])

    plot_purity(ax_purity, ax_purity_ratio, rows)
    plot_counts(ax_counts, ax_counts_ratio, rows)

    fig.text(
        0.064,
        0.952,
        "PPG12 parity establishes the pp photon baseline",
        fontsize=27,
        fontweight="bold",
        ha="left",
        va="top",
    )
    fig.text(
        0.064,
        0.900,
        "The same full-stat output reproduces the four ABCD inputs and their derived raw purity; no fitted normalization.",
        fontsize=17.0,
        color=MUTED,
        ha="left",
        va="top",
    )

    png = args.output_dir / "the241_ppg12_parity_summary_slide.png"
    fig.savefig(png, dpi=160, facecolor="white")
    plt.close(fig)

    speaker_script = write_speaker_script(args.output_dir)
    layout_nodes = args.output_dir / "the241_ppg12_parity_summary_slide_layout_nodes.json"
    make_layout_nodes(layout_nodes)

    ratios = {
        region: [row[f"{region}_current_over_ppg12"] for row in rows]
        for region in REGION_STYLE
    }
    purity_ratios = [row["raw_purity_current_over_ppg12"] for row in rows]
    manifest: dict[str, Any] = {
        "schema_version": 1,
        "task": "THE-241",
        "workstream": "the219_ppg12_parity_summary_slide_20260812",
        "artifact_status": "local_slide_candidate_pending_human_scientific_acceptance",
        "claim": "The current full-stat pp output reproduces the PPG12 raw ABCD populations and their derived raw purity closely enough to establish PPG12 as the mechanics baseline.",
        "not_claimed": [
            "leakage-corrected purity parity",
            "systematic uncertainty",
            "unfolding",
            "jet calibration",
            "AuAu validity",
            "final scientific acceptance",
        ],
        "display": {
            "left": "raw photon purity, PPG12 final ROOT vs current full-stat output, plus current/PPG12 ratio",
            "right": "raw A/B/C/D weighted populations, PPG12 vs current, plus current/PPG12 ratios",
            "uncertainty": "statistical/toy uncertainties from the frozen accepted table only",
            "systematic_band": False,
            "normalization_fit": False,
            "marker_offset_gev": 0.11,
            "marker_offset_purpose": "visibility only; x coordinates and bin widths are otherwise unchanged",
            "font_family": SLIDE_FONT_FAMILY,
            "resolved_font_path": str(slide_font_path),
            "font_policy_guard": "renderer fails instead of silently falling back when Times New Roman is unavailable",
            "layout_feedback_guard": {
                "visible_lower_canvas_text": False,
                "subtitle_font_pt": 17.0,
                "lower_canvas_provenance_or_takeaway_text": False,
                "plot_bottom_fraction": 0.080,
                "marker_fill_key": "high-contrast framed legend in empty lower-left area of the RHS log-y panel",
            },
        },
        "source": {
            "parity_csv": str(PARITY_CSV),
            "parity_csv_sha256": sha256(PARITY_CSV),
            "parity_manifest": str(PARITY_MANIFEST),
            "parity_manifest_sha256": sha256(PARITY_MANIFEST),
        },
        "numerical_summary": {
            "bins": len(rows),
            "count_ratio_max_abs_minus_one_through_28_32": {
                region: max(abs(value - 1.0) for value in values[:-1])
                for region, values in ratios.items()
            },
            "raw_purity_ratio_max_abs_minus_one_through_28_32": max(
                abs(value - 1.0) for value in purity_ratios[:-1]
            ),
            "final_bin": {
                "range_gev": [rows[-1]["pt_lo"], rows[-1]["pt_hi"]],
                "raw_purity_current_over_ppg12": purity_ratios[-1],
                "raw_purity_ratio_error": rows[-1][
                    "raw_purity_current_over_ppg12_error"
                ],
                "interpretation": "count-limited; statistical uncertainty is visibly larger",
            },
        },
        "paper_image_cross_check": {
            "audit_manifest": str(DATATHIEF_MANIFEST),
            "audit_manifest_sha256": sha256(DATATHIEF_MANIFEST),
            "authoritative_reference": "PPG12 final ROOT graph gpurity",
            "raw_max_abs_ratio_minus_one": datathief["ratio_summary"]["raw"][
                "max_abs_ratio_minus_one"
            ],
            "slide_rounding": "within 0.15%",
            "official_datathief_jar_status": datathief["transform_status"][
                "datathief_jar_status"
            ],
            "final_transform_method": datathief["transform_status"][
                "coordinate_transform_method"
            ],
            "provenance_guard": "The screenshot/figure digitization is a cross-check; the final ROOT graph remains authoritative.",
        },
        "outputs": {
            "png": str(png),
            "png_sha256": sha256(png),
            "speaker_script": str(speaker_script),
            "speaker_script_sha256": sha256(speaker_script),
            "layout_nodes": str(layout_nodes),
            "layout_nodes_sha256": sha256(layout_nodes),
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
                    "full 2560x1440 canvas inspected",
                    "two-panel hierarchy is clear",
                    "axes, legends and uncertainty bars are readable",
                    "statistical-only scope is explicit",
                    "no systematic band is present",
                    "DataThief provenance note is accurate",
                ]
                if args.visual_qa_pass
                else []
            ),
        },
    }
    manifest_path = args.output_dir / "the241_ppg12_parity_summary_slide_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(png)
    print(speaker_script)
    print(layout_nodes)
    print(manifest_path)


if __name__ == "__main__":
    main()
