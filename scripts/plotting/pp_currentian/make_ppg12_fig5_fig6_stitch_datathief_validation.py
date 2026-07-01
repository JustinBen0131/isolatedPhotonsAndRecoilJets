#!/usr/bin/env python3
"""Validate PPG12 IAN Fig. 5/6 stitched spectra against SDCC ROOT pulls.

The figure PNGs are treated as the DataThief source. Colored marker centers are
detected from those PNGs, exported through the official DataThief coordinate
transform, then overlaid with the PPG12 SDCC ROOT extraction used by the local
current-IAN validation.
"""

from __future__ import annotations

import csv
import json
import math
import sys
import time
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from PIL import Image


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.append(str(SCRIPT_DIR))

from make_ppg12_datathief_validation_overlays import (  # noqa: E402
    DATATHIEF_JAR,
    jar_md5,
    load_datathief_csv,
    run_datathief_export,
)


CURRENT_IAN = (
    REPO
    / "dataOutput/ppPhotonMLPipeline/"
    / "ppg12_basev3E_currentIAN_fullsim_20260521_1811"
)
ROOTFIT_CSV = (
    CURRENT_IAN
    / "validation/currentIAN_stitching/"
    / "ppg12_currentian_rootfit_stitch_points.csv"
)
PHOTON_IMAGE = (
    CURRENT_IAN
    / "slide_assets/pp_currentIAN_photon_truth_stitch_ppg12_exact_style.png"
)
JET_IMAGE = (
    CURRENT_IAN
    / "slide_assets/pp_currentIAN_inclusive_jet_truth_stitch_ppg12_exact_style.png"
)

OUT_DIR = (
    REPO
    / "dataOutput/the85_auau_xjgamma_unfolding_push/"
    / "ppg12_stitching_fig5_fig6_validation"
)


SAMPLE_COLORS = {
    "photon5": "#d62aa0",
    "photon10": "#2ca02c",
    "photon20": "#1296f3",
    "jet8": "#d62aa0",
    "jet12": "#2ca02c",
    "jet20": "#1296f3",
    "jet30": "#ff6f00",
    "jet40": "#d62aa0",
}

FIGURE_SPECS = {
    "photon": {
        "image": PHOTON_IMAGE,
        "out": OUT_DIR / "ppg12_ian_fig5_photon_stitch_datathief_vs_sdcc_overlay_ratio.png",
        "manifest": OUT_DIR / "ppg12_ian_fig5_photon_stitch_datathief_vs_sdcc_manifest.json",
        "datathief_csv": OUT_DIR / "ppg12_ian_fig5_photon_stitch_datathief_points.csv",
        "samples": ["photon5", "photon10", "photon20"],
        "title": "PPG12 IAN Fig. 5: photon+jet stitching",
        "xlim": (10.0, 40.0),
        "visible_x_min": 10.0,
        "x_refs": ((160.0, 10.0), (1071.0, 40.0)),
        "ylim": (0.03, 1.0e5),
        "plot_ylim": (0.03, 3.5e4),
        "log_y_refs": ((804.0, -1.0), (164.0, 4.0)),
        "xlabel": r"Leading $E_T^\gamma$ [GeV]",
        "ylabel": r"$d\sigma/dE_T^\gamma$ [pb / GeV]",
        "source_note": "Figure 5 photon 5/10/20 GeV stitched leading-truth-photon spectrum",
        "sdcc_sources": [
            "/sphenix/user/shuhangli/ppg12/plotting/plot_combine_uncut.C",
            "/sphenix/user/shuhangli/ppg12/plotting/photon_max_pT_uncut.root",
        ],
    },
    "jet": {
        "image": JET_IMAGE,
        "out": OUT_DIR / "ppg12_ian_fig6_inclusivejet_stitch_datathief_vs_sdcc_overlay_ratio.png",
        "manifest": OUT_DIR / "ppg12_ian_fig6_inclusivejet_stitch_datathief_vs_sdcc_manifest.json",
        "datathief_csv": OUT_DIR / "ppg12_ian_fig6_inclusivejet_stitch_datathief_points.csv",
        "samples": ["jet8", "jet12", "jet20", "jet30", "jet40"],
        "title": "PPG12 IAN Fig. 6: inclusive-jet stitching",
        "xlim": (9.0, 50.0),
        "visible_x_min": 10.0,
        "x_refs": ((160.0, 9.0), (1071.0, 50.0)),
        "ylim": (1e4, 1e12),
        "plot_ylim": (1e4, 1e12),
        "log_y_refs": ((873.0, 4.0), (94.0, 12.0)),
        "xlabel": r"Leading $p_T^\mathrm{jet}$ [GeV]",
        "ylabel": "counts",
        "source_note": "Figure 6 jet 8/12/20/30/40 GeV stitched leading-truth-jet spectrum",
        "sdcc_sources": [
            "/sphenix/user/shuhangli/ppg12/efficiencytool/TruthSpectrumOverlay.C",
            "/sphenix/user/shuhangli/ppg12/efficiencytool/results/truth_spectrum_jet*.root",
        ],
    },
}

# Frame coordinates of the saved PPG12 ROOT-style PNGs. These are extracted from
# the black axis lines in the preserved figure images.
FRAME = {
    "note": "Axis calibration is per figure through x_refs and log_y_refs.",
}


def finite_export_count(raw: dict[str, list[tuple[float, float]]]) -> int:
    return sum(
        1
        for pts in raw.values()
        for x, y in pts
        if math.isfinite(float(x)) and math.isfinite(float(y))
    )


def hex_to_rgb(color: str) -> np.ndarray:
    color = color.lstrip("#")
    return np.asarray([int(color[i : i + 2], 16) for i in (0, 2, 4)], dtype=float)


def load_source_points(
    group: str,
    samples: list[str],
    xlim: tuple[float, float],
    visible_x_min: float | None = None,
) -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = []
    with ROOTFIT_CSV.open(newline="") as f:
        for row in csv.DictReader(f):
            if row["group"] != group or row["sample"] not in samples:
                continue
            if row["used_in_stitch"] != "1":
                continue
            x = float(row["bin_center"])
            y = float(row["value"])
            if not (xlim[0] <= x <= xlim[1]) or y <= 0:
                continue
            if visible_x_min is not None and x < visible_x_min:
                continue
            rows.append(
                {
                    "sample": row["sample"],
                    "x": x,
                    "y": y,
                    "ey": float(row["error"]),
                }
            )
    return rows


def data_to_pixel(
    x: float,
    y: float,
    x_refs: tuple[tuple[float, float], tuple[float, float]],
    log_y_refs: tuple[tuple[float, float], tuple[float, float]],
) -> tuple[float, float]:
    (px0, x0), (px1, x1) = x_refs
    px = px0 + (x - x0) / (x1 - x0) * (px1 - px0)
    logy = math.log10(y)
    (py0, log0), (py1, log1) = log_y_refs
    py = py0 + (logy - log0) / (log1 - log0) * (py1 - py0)
    return px, py


def find_marker_pixels(
    image: Path,
    source_rows: list[dict[str, float | str]],
    x_refs: tuple[tuple[float, float], tuple[float, float]],
    log_y_refs: tuple[tuple[float, float], tuple[float, float]],
) -> list[dict[str, float | str]]:
    arr = np.asarray(Image.open(image).convert("RGB"), dtype=float)
    out: list[dict[str, float | str]] = []
    for row in source_rows:
        sample = str(row["sample"])
        target = hex_to_rgb(SAMPLE_COLORS[sample])
        dist = np.sqrt(np.sum((arr - target) ** 2, axis=2))
        color_mask = dist < 105.0

        px0, py0 = data_to_pixel(float(row["x"]), float(row["y"]), x_refs, log_y_refs)
        xlo = max(0, int(round(px0)) - 8)
        xhi = min(arr.shape[1] - 1, int(round(px0)) + 8)
        ylo = max(0, int(round(py0)) - 8)
        yhi = min(arr.shape[0] - 1, int(round(py0)) + 8)
        sub = color_mask[ylo : yhi + 1, xlo : xhi + 1]
        yy, xx = np.where(sub)
        color_match_mode = "sample_color"
        if len(xx) < 2:
            fxlo = max(0, int(round(px0)) - 16)
            fxhi = min(arr.shape[1] - 1, int(round(px0)) + 16)
            fylo = max(0, int(round(py0)) - 16)
            fyhi = min(arr.shape[0] - 1, int(round(py0)) + 16)
            crop = arr[fylo : fyhi + 1, fxlo : fxhi + 1]
            local_max = crop.max(axis=2)
            local_min = crop.min(axis=2)
            saturated = (local_max - local_min > 45.0) & (local_max > 120.0)
            # The PPG12 inclusive-jet figure draws the red fit through the
            # highest jet8 markers. In those bins the visible marker pixels are
            # mostly fit-line red, so use saturated colored pixels in the same
            # local window rather than failing or grabbing black axis pixels.
            yy, xx = np.where(saturated)
            abs_x_offset = fxlo
            abs_y_offset = fylo
            color_match_mode = "saturated_local_fallback"
        else:
            abs_x_offset = xlo
            abs_y_offset = ylo
        if len(xx) < 2:
            raise RuntimeError(
                f"Could not find visible colored marker near {sample} x={row['x']} y={row['y']} "
                f"at predicted pixel ({px0:.1f}, {py0:.1f}) in {image}"
            )

        abs_x = abs_x_offset + xx.astype(float)
        abs_y = abs_y_offset + yy.astype(float)
        # Weight central pixels most strongly so anti-aliased error bars do not
        # pull the marker center away from the plotted point.
        r2 = (abs_x - px0) ** 2 + (abs_y - py0) ** 2
        weights = np.exp(-0.5 * r2 / (3.0**2))
        out.append(
            {
                **row,
                "pixel_x": float(np.average(abs_x, weights=weights)),
                "pixel_y": float(np.average(abs_y, weights=weights)),
                "pred_pixel_x": float(px0),
                "pred_pixel_y": float(py0),
                "pixels_used": int(len(xx)),
                "color_match_mode": color_match_mode,
            }
        )
    return out


def export_datathief(
    image: Path,
    pixel_rows: list[dict[str, float | str]],
    x_refs: tuple[tuple[float, float], tuple[float, float]],
    log_y_refs: tuple[tuple[float, float], tuple[float, float]],
    out_csv: Path,
) -> dict[str, list[tuple[float, float]]]:
    (py0, log0), (py1, log1) = log_y_refs
    # DataThief is asked to export a stable internal 0-100 y coordinate, then
    # Python converts that coordinate back to log10(y). The jar is fragile with
    # raw log-axis references on this old ROOT-style plot.
    internal_y0 = 0.0
    internal_y1 = 100.0
    (px0, x0), (px1, x1) = x_refs
    refs = [
        (0, px0, py0, x0, internal_y0),
        (1, px1, py0, x1, internal_y0),
        (2, px0, py1, x0, internal_y1),
    ]
    points = [(str(r["sample"]), float(r["pixel_x"]), float(r["pixel_y"])) for r in pixel_rows]

    expected_points = len(points)
    if out_csv.exists():
        existing = load_datathief_csv(out_csv)
        existing_points = sum(len(pts) for pts in existing.values())
        if existing_points == expected_points and finite_export_count(existing) == expected_points:
            raw = existing
        else:
            raw = {}
    else:
        raw = {}

    if not raw:
        last_error: Exception | None = None
        for attempt in range(1, 5):
            try:
                run_datathief_export(image, refs, points, out_csv)
                raw = load_datathief_csv(out_csv)
            except Exception as exc:  # DataThief can hang/fail in the legacy Java path.
                last_error = exc
                raw = {}
            if finite_export_count(raw) > 0:
                break
            if attempt < 4:
                time.sleep(0.5)
        if finite_export_count(raw) == 0 and last_error is not None:
            raise RuntimeError(f"DataThief export failed after retries for {image}") from last_error
    if finite_export_count(raw) == 0:
        raise RuntimeError(f"DataThief returned no finite points after retries for {image}")
    converted: dict[str, list[tuple[float, float]]] = {}
    for sample, pts in raw.items():
        converted[sample] = [
            (x, 10.0 ** (log0 + (internal_y - internal_y0) / (internal_y1 - internal_y0) * (log1 - log0)))
            for x, internal_y in pts
        ]
    return converted


def write_output(group: str) -> None:
    spec = FIGURE_SPECS[group]
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    source_rows = load_source_points(group, spec["samples"], spec["xlim"], spec.get("visible_x_min"))
    pixel_rows = find_marker_pixels(spec["image"], source_rows, spec["x_refs"], spec["log_y_refs"])
    exported = export_datathief(spec["image"], pixel_rows, spec["x_refs"], spec["log_y_refs"], spec["datathief_csv"])

    # Preserve row ordering while matching DataThief series order.
    sample_indices = {sample: 0 for sample in spec["samples"]}
    plotted_rows: list[dict[str, float | str]] = []
    for row in pixel_rows:
        sample = str(row["sample"])
        idx = sample_indices[sample]
        sample_indices[sample] += 1
        dt_x, dt_y = exported[sample][idx]
        ratio = dt_y / float(row["y"]) if float(row["y"]) > 0 else float("nan")
        plotted_rows.append(
            {
                **row,
                "dt_x": float(dt_x),
                "dt_y": float(dt_y),
                "ratio": float(ratio),
            }
        )

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "dejavuserif",
            "axes.linewidth": 1.1,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )
    fig, (ax, rax) = plt.subplots(
        2,
        1,
        figsize=(7.72, 9.98),
        dpi=200,
        sharex=True,
        gridspec_kw={"height_ratios": [3.25, 1.0], "hspace": 0.04},
    )
    ax.set_yscale("log")
    ax.set_xlim(*spec["xlim"])
    ax.set_ylim(*spec["plot_ylim"])
    for sample in spec["samples"]:
        rows = [r for r in plotted_rows if r["sample"] == sample]
        if not rows:
            continue
        color = SAMPLE_COLORS[sample]
        x = np.asarray([float(r["x"]) for r in rows])
        y = np.asarray([float(r["y"]) for r in rows])
        ey = np.asarray([float(r["ey"]) for r in rows])
        dtx = np.asarray([float(r["dt_x"]) for r in rows])
        dty = np.asarray([float(r["dt_y"]) for r in rows])
        rr = np.asarray([float(r["ratio"]) for r in rows])
        finite = np.isfinite(dtx) & np.isfinite(dty) & np.isfinite(rr)
        ax.errorbar(
            x,
            y,
            yerr=ey,
            fmt="o",
            color=color,
            markerfacecolor="white",
            markeredgewidth=1.2,
            ms=4.8,
            lw=0.9,
            zorder=3,
        )
        ax.plot(
            dtx[finite],
            dty[finite],
            "o",
            color=color,
            ms=3.4,
            zorder=4,
        )
        rax.plot(dtx[finite], rr[finite], "o", color=color, ms=3.4, zorder=3)

    ax.set_ylabel(spec["ylabel"], fontsize=19)
    rax.set_ylabel("DataThief / SDCC", fontsize=15)
    rax.set_xlabel(spec["xlabel"], fontsize=19)
    ax.tick_params(which="both", labelsize=14, length=6)
    rax.tick_params(which="both", labelsize=13, length=5)
    ax.minorticks_on()
    rax.minorticks_on()
    rax.axhline(1.0, color="0.25", lw=1.0, ls=(0, (4, 4)))
    rax.set_ylim(0.965, 1.035)

    ax.text(
        0.47,
        0.985,
        r"$\bf{\it{sPHENIX}}$ Internal" + "\n" + r"$p{+}p$ $\sqrt{s}=200$ GeV" + "\nPYTHIA8",
        transform=ax.transAxes,
        fontsize=15,
        va="top",
    )
    sample_handles = [
        Line2D([0], [0], marker="o", linestyle="none", color=SAMPLE_COLORS[s], label=s, markersize=6)
        for s in spec["samples"]
    ]
    marker_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            color="black",
            label="closed: DataThief",
            markersize=5,
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            color="black",
            markerfacecolor="white",
            markeredgewidth=1.2,
            label="open: SDCC ROOT",
            markersize=6,
        ),
    ]
    if group == "jet":
        sample_legend_kwargs = {
            "loc": "lower left",
            "bbox_to_anchor": (0.035, 0.045),
            "ncol": 2,
        }
        marker_legend_kwargs = {
            "loc": "upper right",
            "bbox_to_anchor": (0.985, 0.82),
            "ncol": 1,
        }
    else:
        sample_legend_kwargs = {
            "loc": "upper right",
            "bbox_to_anchor": (0.98, 0.80),
            "ncol": 1,
        }
        marker_legend_kwargs = {
            "loc": "lower left",
            "bbox_to_anchor": (0.03, 0.03),
            "ncol": 1,
        }
    leg1 = ax.legend(
        handles=sample_handles,
        **sample_legend_kwargs,
        frameon=False,
        fontsize=15,
        handlelength=1.1,
        handletextpad=0.4,
        columnspacing=0.9,
    )
    ax.add_artist(leg1)
    ax.legend(
        handles=marker_handles,
        **marker_legend_kwargs,
        frameon=False,
        fontsize=14.5,
        handlelength=1.1,
        handletextpad=0.5,
    )
    fig.subplots_adjust(left=0.16, right=0.98, top=0.975, bottom=0.08)
    fig.savefig(spec["out"])
    plt.close(fig)

    ratios = np.asarray([float(r["ratio"]) for r in plotted_rows], dtype=float)
    finite_ratios = ratios[np.isfinite(ratios)]
    manifest = {
        "artifact": str(spec["out"]),
        "source_note": spec["source_note"],
        "sdcc_sources": spec["sdcc_sources"],
        "source_image": str(spec["image"]),
        "sdcc_source_csv": str(ROOTFIT_CSV),
        "datathief_export_csv": str(spec["datathief_csv"]),
        "datathief_jar": str(DATATHIEF_JAR),
        "datathief_jar_md5": jar_md5(),
        "frame_pixels": FRAME,
        "xlim": list(spec["xlim"]),
        "visible_x_min": spec.get("visible_x_min"),
        "x_refs": [list(v) for v in spec["x_refs"]],
        "ylim": list(spec["plot_ylim"]),
        "log_y_refs": [list(v) for v in spec["log_y_refs"]],
        "n_points": len(plotted_rows),
        "mean_abs_ratio_minus_one": float(np.mean(np.abs(finite_ratios - 1.0))) if len(finite_ratios) else None,
        "max_abs_ratio_minus_one": float(np.max(np.abs(finite_ratios - 1.0))) if len(finite_ratios) else None,
        "samples": spec["samples"],
        "finite_datathief_points": int(sum(np.isfinite(float(r["dt_y"])) for r in plotted_rows)),
    }
    spec["manifest"].write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(spec["out"])
    print(spec["manifest"])
    mean = manifest["mean_abs_ratio_minus_one"]
    max_abs = manifest["max_abs_ratio_minus_one"]
    mean_text = "None" if mean is None else f"{mean:.4g}"
    max_text = "None" if max_abs is None else f"{max_abs:.4g}"
    print(
        f"{group}: n={len(plotted_rows)} finite={len(finite_ratios)} "
        f"mean|ratio-1|={mean_text} max|ratio-1|={max_text}"
    )


def main() -> None:
    write_output("photon")
    write_output("jet")


if __name__ == "__main__":
    main()
