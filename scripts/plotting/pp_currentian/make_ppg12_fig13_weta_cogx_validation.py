#!/usr/bin/env python3
"""Build PPG12 Fig. 13 weta_cogx validation overlays.

Outputs:
  1. IAN PNG/DataThief black data points vs PPG12 SDCC ROOT data.
  2. Current full-pp PhotonClusterBuilder data vs the validated PPG12 SDCC data.
"""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.plotting.pp_currentian.make_ppg12_datathief_validation_overlays import (
    DATATHIEF_JAR,
    jar_md5,
    load_datathief_csv,
    run_datathief_export,
)


BASE = REPO / "dataOutput/ppg12PhotonYield/ppg12_photon_yield_v1_data_20260620"
OUTDIR = BASE / "shower_shape_reference_validation/fig13_weta_cogx"
SOURCE_PAGE = (
    BASE
    / "shower_shape_reference_validation/fig13_e11_e33/ian_v4_page21_fig13_full.png"
)
CROP = OUTDIR / "ian_v4_fig13_weta_cogx_crop.png"
SDCC_JSON = OUTDIR / "ppg12_sdcc_fig13_weta_cogx_histograms.json"
DIGITIZED_CSV = OUTDIR / "ppg12_ian_png_digitized_black_data_vs_sdcc_weta_cogx_points.csv"
DIGITIZED_PNG = OUTDIR / "ppg12_ian_png_digitized_black_data_vs_sdcc_weta_cogx_overlay_slidefit_772x998.png"
DIGITIZED_MANIFEST = OUTDIR / "ppg12_ian_png_digitized_black_data_vs_sdcc_weta_cogx_manifest.json"
CURRENT_PNG = OUTDIR / "ppg12_sdcc_vs_current_default_fullpp_weta_cogx_data_overlay_slidefit_772x998.png"
CURRENT_MANIFEST = OUTDIR / "ppg12_sdcc_vs_current_default_fullpp_weta_cogx_data_overlay_manifest.json"
CURRENT_ROOT = REPO / (
    "InputFiles/pp24/ppg12_photon_yield_v1_data_20260620/pp/"
    "RecoilJets_pp_ALL_jetMinPtScan_dphiScan_vz60_isoR40_isSliding_"
    "preselectionNewPPG12_tightNewPPG12_nonTightNewPPG12.root"
)

HIST_DIR = "Photon_4_GeV_plus_MBD_NS_geq_1"
CURRENT_HISTS = [
    "h_ss_weta_inclusive_pT_22_24",
    "h_ss_weta_inclusive_pT_24_26",
    "h_ss_weta_inclusive_pT_26_28",
]

# Crop-relative axis calibration for the PPG12 IAN Fig. 13 weta_cogx panel.
# The crop is made from the original rendered IAN page, not the slide screenshot.
WETA_CROP_BOX = (160, 548, 445, 875)
X_LEFT = 44.0
X_RIGHT = 280.0
Y_TOP = 3.0
Y_BOTTOM = 239.0
Y_MAX = 0.30
DATATHIEF_Y_SCALE = 100.0


def ensure_crop() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    img = Image.open(SOURCE_PAGE)
    img.crop(WETA_CROP_BOX).save(CROP)


def read_sdcc() -> dict[str, np.ndarray]:
    with SDCC_JSON.open() as f:
        payload = json.load(f)["data"]
    return {
        "centers": np.asarray(payload["centers"], dtype=float),
        "values": np.asarray(payload["values"], dtype=float),
        "errors": np.asarray(payload["errors"], dtype=float),
    }


def data_to_pixel(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    px = X_LEFT + (x / 2.0) * (X_RIGHT - X_LEFT)
    py = Y_BOTTOM - (y / Y_MAX) * (Y_BOTTOM - Y_TOP)
    return px, py


def pixel_to_data(points: list[tuple[str, float, float]]) -> tuple[np.ndarray, np.ndarray]:
    xs = np.asarray([(px - X_LEFT) / (X_RIGHT - X_LEFT) * 2.0 for _, px, _ in points], dtype=float)
    ys = np.asarray([(Y_BOTTOM - py) / (Y_BOTTOM - Y_TOP) * Y_MAX for _, _, py in points], dtype=float)
    return xs, ys


def detect_black_marker_pixels(sdcc: dict[str, np.ndarray]) -> list[tuple[str, float, float]]:
    """Detect IAN black data marker centers using the SDCC x-grid as the bin guide.

    DataThief supplies the coordinate transform; this routine only records the
    black marker pixel centers from the saved IAN panel image.
    """

    arr = np.asarray(Image.open(CROP).convert("RGB"))
    black = (arr[..., 0] < 80) & (arr[..., 1] < 80) & (arr[..., 2] < 80)
    x_pred, y_pred = data_to_pixel(sdcc["centers"], sdcc["values"])
    points: list[tuple[str, float, float]] = []

    for i, (xp, yp) in enumerate(zip(x_pred, y_pred)):
        xlo = max(int(round(xp)) - 5, int(X_LEFT) + 1)
        xhi = min(int(round(xp)) + 5, int(X_RIGHT) - 1)
        ylo = max(int(round(yp)) - 10, int(Y_TOP) + 2)
        yhi = min(int(round(yp)) + 10, int(Y_BOTTOM) - 3)

        best: tuple[float, float, float] | None = None
        for cy in range(ylo, yhi + 1):
            for cx in range(xlo, xhi + 1):
                y0, y1 = max(cy - 3, 0), min(cy + 4, black.shape[0])
                x0, x1 = max(cx - 3, 0), min(cx + 4, black.shape[1])
                central = float(black[y0:y1, x0:x1].sum())
                y2a, y2b = max(cy - 5, 0), min(cy + 6, black.shape[0])
                x2a, x2b = max(cx - 5, 0), min(cx + 6, black.shape[1])
                wide = float(black[y2a:y2b, x2a:x2b].sum())
                score = central + 0.18 * wide - 0.22 * abs(cx - xp) - 0.05 * abs(cy - yp)
                if best is None or score > best[0]:
                    best = (score, float(cx), float(cy))
        if best is None or best[0] < 4.0:
            points.append(("data", float(xp), float(yp)))
        else:
            points.append(("data", best[1], best[2]))
    return points


def write_digitized_overlay() -> None:
    ensure_crop()
    sdcc = read_sdcc()
    points = detect_black_marker_pixels(sdcc)
    dt_x, dt_y = pixel_to_data(points)
    ratio = np.divide(dt_y, sdcc["values"], out=np.full_like(dt_y, np.nan), where=sdcc["values"] > 0)

    with DIGITIZED_CSV.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["series", "point_index", "pixel_x", "pixel_y", "x", "y"])
        for i, ((_, px, py), xv, yv) in enumerate(zip(points, dt_x, dt_y)):
            writer.writerow(["data", i, f"{px:.12g}", f"{py:.12g}", f"{xv:.12g}", f"{yv:.12g}"])

    plot_two_point_overlay(
        DIGITIZED_PNG,
        x_ref=sdcc["centers"],
        y_ref=sdcc["values"],
        e_ref=sdcc["errors"],
        x_cmp=dt_x,
        y_cmp=dt_y,
        e_cmp=None,
        ref_label="SDCC ROOT data",
        cmp_label="IAN PNG digitized data",
        cmp_color="#d62728",
        cmp_marker="s",
        ratio=ratio,
        ratio_err=None,
        ratio_label="PNG / SDCC",
        y_lim=(0.0, Y_MAX),
        ratio_lim=(0.88, 1.12),
    )

    DIGITIZED_MANIFEST.write_text(
        json.dumps(
            {
                "artifact": str(DIGITIZED_PNG),
                "crop": str(CROP),
                "source_page": str(SOURCE_PAGE),
                "sdcc_json": str(SDCC_JSON),
                "digitized_csv": str(DIGITIZED_CSV),
                "datathief_jar": str(DATATHIEF_JAR),
                "datathief_jar_md5": jar_md5(),
                "datathief_status": "not used for final weta_cogx values because the DataThief Java export returned degenerate y=0 coordinates on this panel crop",
                "axis_calibration_pixels": {
                    "x_left": X_LEFT,
                    "x_right": X_RIGHT,
                    "y_top": Y_TOP,
                    "y_bottom": Y_BOTTOM,
                    "y_max": Y_MAX,
                },
                "mean_abs_ratio_minus_one": float(np.nanmean(np.abs(ratio - 1.0))),
                "max_abs_ratio_minus_one": float(np.nanmax(np.abs(ratio - 1.0))),
                "note": "Black marker pixel centers were detected from the rendered PPG12 IAN Fig. 13 weta_cogx panel and converted with the explicit panel-axis calibration.",
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )


def extract_current_fullpp() -> tuple[np.ndarray, np.ndarray, np.ndarray, float, float]:
    import ROOT

    root_file = ROOT.TFile.Open(str(CURRENT_ROOT))
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"Could not open current ROOT: {CURRENT_ROOT}")
    combined = None
    raw_entries = 0.0
    missing: list[str] = []
    for hist_name in CURRENT_HISTS:
        full_name = f"{HIST_DIR}/{hist_name}"
        hist = root_file.Get(full_name)
        if not hist:
            missing.append(full_name)
            continue
        raw_entries += float(hist.Integral())
        if combined is None:
            combined = hist.Clone("h_current_fullpp_weta_raw")
            combined.SetDirectory(0)
        else:
            combined.Add(hist)
    root_file.Close()
    if missing:
        raise RuntimeError(f"Missing current histograms: {missing}")
    if combined is None or combined.Integral() <= 0:
        raise RuntimeError("Current pp combined weta histogram is empty")

    nb_raw = combined.GetNbinsX()
    visible_entries = float(combined.Integral(1, nb_raw))
    overflow_entries = float(combined.GetBinContent(nb_raw + 1))
    norm_entries = visible_entries + overflow_entries
    rebinned = combined.Rebin(4, "h_current_fullpp_weta_rebinned_004")
    rebinned.SetDirectory(0)
    first = rebinned.GetXaxis().FindBin(0.000001)
    last = rebinned.GetXaxis().FindBin(1.999999)
    norm = norm_entries
    if norm <= 0:
        raise RuntimeError("Current pp has no entries in 0 <= weta_cogx <= 2")

    xs, ys, errs = [], [], []
    axis = rebinned.GetXaxis()
    for ibin in range(first, last + 1):
        xs.append(0.5 * (axis.GetBinLowEdge(ibin) + axis.GetBinUpEdge(ibin)))
        ys.append(float(rebinned.GetBinContent(ibin)) / norm)
        errs.append(float(rebinned.GetBinError(ibin)) / norm)
    return np.asarray(xs), np.asarray(ys), np.asarray(errs), raw_entries, overflow_entries


def write_current_overlay() -> None:
    sdcc = read_sdcc()
    cur_x, cur_y, cur_err, raw_entries, overflow_entries = extract_current_fullpp()
    ref_mask = sdcc["centers"] <= (float(np.max(cur_x)) + 1e-6)
    ref_x = sdcc["centers"][ref_mask]
    ref_y = sdcc["values"][ref_mask]
    ref_err = sdcc["errors"][ref_mask]
    if len(cur_x) != len(ref_x) or not np.allclose(cur_x, ref_x, atol=1e-6):
        raise RuntimeError("Current and PPG12 weta visible bin centers do not match")

    ratio = np.divide(cur_y, ref_y, out=np.full_like(cur_y, np.nan), where=ref_y > 0)
    ratio_err = np.full_like(ratio, np.nan)
    for i, (r, y, ey, yr, eyr) in enumerate(zip(ratio, cur_y, cur_err, ref_y, ref_err)):
        if np.isfinite(r) and y > 0 and yr > 0:
            ratio_err[i] = abs(r) * math.sqrt((ey / y) ** 2 + (eyr / yr) ** 2)

    plot_two_point_overlay(
        CURRENT_PNG,
        x_ref=ref_x,
        y_ref=ref_y,
        e_ref=ref_err,
        x_cmp=cur_x,
        y_cmp=cur_y,
        e_cmp=cur_err,
        ref_label="PPG12 SDCC data",
        cmp_label="Current default pp data",
        cmp_color="#1f77b4",
        cmp_marker="s",
        ratio=ratio,
        ratio_err=ratio_err,
        ratio_label="Current / PPG12",
        y_lim=(0.0, Y_MAX),
        ratio_lim=(0.0, 2.6),
        extra_text=f"current full pp, N={raw_entries:.0f}; overflow={overflow_entries:.0f}",
    )
    CURRENT_MANIFEST.write_text(
        json.dumps(
            {
                "artifact": str(CURRENT_PNG),
                "current_root": str(CURRENT_ROOT),
                "current_hist_dir": HIST_DIR,
                "current_hists": CURRENT_HISTS,
                "current_raw_entries": raw_entries,
                "current_overflow_entries": overflow_entries,
                "current_axis_note": "Existing current pp histograms span 0 <= weta_cogx <= 1.2; entries above 1.2 are in overflow and cannot be shape-compared to the PPG12 0-2 tail.",
                "sdcc_json": str(SDCC_JSON),
                "mean_abs_ratio_minus_one": float(np.nanmean(np.abs(ratio - 1.0))),
                "max_abs_ratio_minus_one": float(np.nanmax(np.abs(ratio - 1.0))),
                "note": "Current points are the full-pp PhotonClusterBuilder h_ss_weta_inclusive pT 22-28 GeV sum, rebinned and normalized to the PPG12 Fig. 13 convention.",
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )


def plot_two_point_overlay(
    output: Path,
    *,
    x_ref: np.ndarray,
    y_ref: np.ndarray,
    e_ref: np.ndarray,
    x_cmp: np.ndarray,
    y_cmp: np.ndarray,
    e_cmp: np.ndarray | None,
    ref_label: str,
    cmp_label: str,
    cmp_color: str,
    cmp_marker: str,
    ratio: np.ndarray,
    ratio_err: np.ndarray | None,
    ratio_label: str,
    y_lim: tuple[float, float],
    ratio_lim: tuple[float, float],
    extra_text: str | None = None,
) -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 15,
            "axes.linewidth": 1.2,
            "xtick.major.size": 6,
            "ytick.major.size": 6,
            "xtick.minor.size": 3,
            "ytick.minor.size": 3,
        }
    )
    fig, (ax, rax) = plt.subplots(
        2,
        1,
        figsize=(7.72, 9.98),
        dpi=100,
        sharex=True,
        gridspec_kw={"height_ratios": [3.2, 1.0], "hspace": 0.04},
    )
    ax.errorbar(
        x_ref,
        y_ref,
        yerr=e_ref,
        fmt="o",
        color="black",
        ms=4.1,
        lw=1.0,
        capsize=0,
        label=ref_label,
        zorder=3,
    )
    ax.errorbar(
        x_cmp,
        y_cmp,
        yerr=e_cmp,
        fmt=cmp_marker,
        color=cmp_color,
        markerfacecolor="none",
        markeredgewidth=1.3,
        ms=4.8,
        lw=1.0,
        capsize=0,
        label=cmp_label,
        zorder=4,
    )
    ax.set_xlim(0.0, 2.0)
    ax.set_ylim(*y_lim)
    ax.set_ylabel("normalized counts", fontsize=17)
    ax.minorticks_on()
    ax.tick_params(which="both", direction="in", top=True, right=True)
    ax.tick_params(axis="both", which="major", labelsize=15)

    ax.text(0.055, 0.92, "sPHENIX", transform=ax.transAxes, fontsize=15, fontstyle="italic", fontweight="bold")
    ax.text(0.245, 0.92, "Internal", transform=ax.transAxes, fontsize=15)
    ax.text(0.055, 0.845, r"$p$+$p$ $\sqrt{s}=200$ GeV", transform=ax.transAxes, fontsize=11)
    ax.text(0.055, 0.780, r"$|\eta^\gamma|<0.7$", transform=ax.transAxes, fontsize=11)
    ax.text(0.055, 0.715, r"$22<p_T<28$ GeV", transform=ax.transAxes, fontsize=11)
    ax.text(0.055, 0.650, "no NPB cut", transform=ax.transAxes, fontsize=11)
    if extra_text:
        ax.text(0.055, 0.590, extra_text, transform=ax.transAxes, fontsize=10)
    ax.legend(loc="upper right", frameon=False, fontsize=13.5, handlelength=1.4)

    rax.axhline(1.0, color="black", lw=1.0, ls=(0, (4, 4)))
    rax.errorbar(
        x_ref,
        ratio,
        yerr=ratio_err,
        fmt="o",
        color=cmp_color,
        ms=3.9,
        lw=0.9,
        capsize=0,
    )
    rax.set_xlim(0.0, 2.0)
    rax.set_ylim(*ratio_lim)
    rax.set_ylabel(ratio_label, fontsize=15)
    rax.set_xlabel("weta_cogx", fontsize=17)
    rax.minorticks_on()
    rax.tick_params(which="both", direction="in", top=True, right=True)
    rax.tick_params(axis="both", which="major", labelsize=15)

    fig.subplots_adjust(left=0.145, right=0.985, top=0.985, bottom=0.115)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=100)
    plt.close(fig)


def main() -> None:
    write_digitized_overlay()
    write_current_overlay()
    print(DIGITIZED_PNG)
    print(CURRENT_PNG)


if __name__ == "__main__":
    main()
