#!/usr/bin/env python3
"""Build BDT data-reference validation plots for PPG12 IAN Fig. 13/19 panels.

This helper follows the same proof pattern as the E11/E33 and weta_cogx checks:

1. Re-export black data-marker pixels from the rendered PPG12 IAN panel through
   the local DataThief jar and compare them to the PPG12 SDCC ROOT projection.
2. Compare that validated SDCC projection to the current completed full-pp
   RecoilJets/PhotonClusterBuilder BDT-score QA histograms.

The output is PNG-only and slide-fit; companion CSV/JSON files preserve exact
input provenance and calibration choices.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from PIL import Image

from make_ppg12_datathief_validation_overlays import (
    DATATHIEF_JAR,
    jar_md5,
    load_datathief_csv,
    run_datathief_export,
)


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
OUTDIR = (
    REPO
    / "dataOutput/ppg12PhotonYield/ppg12_photon_yield_v1_data_20260620"
    / "shower_shape_reference_validation/fig13_19_bdt"
)
SDCC_RAW = OUTDIR / "sdcc_bdt_compact_extract_raw.txt"
SDCC_JSON = OUTDIR / "ppg12_sdcc_fig13_19_bdt_data_projections.json"
SDCC_FOURCURVE_JSON = OUTDIR / "ppg12_sdcc_fig19_bdt_fourcurve_projections.json"

CURRENT_FULL_PP = (
    REPO
    / "InputFiles/pp24/ppg12_photon_yield_v1_data_20260620/pp"
    / "RecoilJets_pp_ALL_jetMinPtScan_dphiScan_vz60_isoR40_isSliding_"
    "preselectionNewPPG12_tightNewPPG12_nonTightNewPPG12.root"
)
CURRENT_DIR = "PPG12_scaledtrigger30"


@dataclass(frozen=True)
class PanelSpec:
    tag: str
    title: str
    page_image: Path
    page_label: str
    sdcc_hist: str
    source_hist_label: str
    ian_label: str
    crop_box: tuple[int, int, int, int]
    x_left_abs: float
    x_right_abs: float
    y_top_abs: float
    y_bottom_abs: float
    y_max: float
    current_hists: tuple[str, ...]
    current_label: str
    ratio_ylim_reference: tuple[float, float]
    ratio_ylim_current: tuple[float, float]


PANELS: tuple[PanelSpec, ...] = (
    PanelSpec(
        tag="bdt_no_npb_22_28",
        title="BDT, 22 < pT < 28 GeV, no NPB cut",
        page_image=OUTDIR / "ian_v4_page-021.png",
        page_label="PPG12 IAN v4 page 21 / Fig. 13 BDT panel",
        sdcc_hist="h2d_bdt_eta0_pt3_cut0",
        source_hist_label="h2d_bdt_eta0_pt3_cut0",
        ian_label="22 < pT < 28 GeV, w/o nbkg cut",
        crop_box=(1600, 480, 2385, 995),
        x_left_abs=1659.0,
        x_right_abs=2339.0,
        y_top_abs=525.0,
        y_bottom_abs=934.0,
        y_max=0.85,
        current_hists=(
            "h_tightBDTScore_allCandidates_pT_22_24",
            "h_tightBDTScore_allCandidates_pT_24_26",
            "h_tightBDTScore_allCandidates_pT_26_28",
        ),
        current_label="Current full-pp all-candidate BDT QA",
        ratio_ylim_reference=(0.85, 1.15),
        ratio_ylim_current=(0.35, 1.65),
    ),
    PanelSpec(
        tag="bdt_with_npb_18_22",
        title="BDT, 18 < pT < 22 GeV, with NPB/preselection",
        page_image=OUTDIR / "ian_v4_page-025.png",
        page_label="PPG12 IAN v4 page 25 / Fig. 19 BDT panel",
        sdcc_hist="h2d_bdt_eta0_pt2_cut1",
        source_hist_label="h2d_bdt_eta0_pt2_cut1",
        ian_label="18 < pT < 22 GeV, w/ nbkg cut",
        crop_box=(1600, 490, 2385, 1085),
        x_left_abs=1659.0,
        x_right_abs=2339.0,
        y_top_abs=525.0,
        y_bottom_abs=1027.0,
        y_max=0.24,
        current_hists=(
            "h_tightBDTScore_preselected_pT_18_20",
            "h_tightBDTScore_preselected_pT_20_22",
        ),
        current_label="Current full-pp preselected BDT QA",
        ratio_ylim_reference=(0.85, 1.15),
        ratio_ylim_current=(0.35, 1.65),
    ),
)


def write_sdcc_json_from_raw() -> dict[str, object]:
    if SDCC_FOURCURVE_JSON.exists():
        payload = json.loads(SDCC_FOURCURVE_JSON.read_text())
        out: dict[str, object] = {
            "source_root": payload["source_files"]["bdt_no_npb_22_28"]["data"],
            "source_files": payload["source_files"],
            "source_selection": payload.get("source_selection"),
            "source_fourcurve_json": str(SDCC_FOURCURVE_JSON),
            "panels": {},
        }
        panels: dict[str, object] = {}
        for spec in PANELS:
            panel = payload["panels"][spec.tag]
            curve = panel["curves"]["data"]
            source_integral = curve.get("integral_after_scaling", float(np.sum(np.asarray(curve["values"], dtype=float))))
            panels[spec.tag] = {
                "tag": spec.tag,
                "source_root": panel["summary"]["source_files"]["data"],
                "histname": spec.sdcc_hist,
                "source_integral_0to1": source_integral,
                "centers": curve["centers"],
                "edges": curve["edges"],
                "values": curve["values"],
                "errors": curve["errors"],
                "raw": curve.get("raw", []),
                "raw_errors": curve.get("raw_errors", []),
                "ian_label": spec.ian_label,
                "page_label": spec.page_label,
                "ymax": spec.y_max,
                "fourcurve_source_summary": panel["summary"],
            }
        out["panels"] = panels
        SDCC_JSON.write_text(json.dumps(out, indent=2, sort_keys=True))
        return out

    text = SDCC_RAW.read_text(errors="replace")
    if "JSON_BEGIN" not in text or "JSON_END" not in text:
        raise RuntimeError(f"Could not find JSON block in {SDCC_RAW}")
    payload = json.loads(text.split("JSON_BEGIN", 1)[1].split("JSON_END", 1)[0].strip())
    chosen = None
    for candidate in payload["candidate_roots"]:
        if candidate.get("path", "").endswith("data_histoshower_shape_showershape.root"):
            chosen = candidate
            break
    if chosen is None:
        chosen = payload["candidate_roots"][0]

    out: dict[str, object] = {
        "source_root": chosen["path"],
        "alternative_roots": [c.get("path") for c in payload["candidate_roots"]],
        "panels": {},
        "raw_extract": str(SDCC_RAW),
    }
    panels: dict[str, object] = {}
    for spec in PANELS:
        info = chosen["hists"][spec.sdcc_hist]
        panels[spec.tag] = {
            "tag": spec.tag,
            "source_root": chosen["path"],
            "histname": spec.sdcc_hist,
            "source_integral_0to1": info["selected_integral_0to1"],
            "centers": info["centers"],
            "edges": info["edges"],
            "values": info["values"],
            "errors": info["errors"],
            "raw": info["raw"],
            "raw_errors": info["raw_errors"],
            "ian_label": spec.ian_label,
            "page_label": spec.page_label,
            "ymax": spec.y_max,
        }
    out["panels"] = panels
    SDCC_JSON.write_text(json.dumps(out, indent=2, sort_keys=True))
    return out


def load_sdcc_projection(spec: PanelSpec) -> dict[str, np.ndarray]:
    payload = json.loads(SDCC_JSON.read_text())
    p = payload["panels"][spec.tag]
    return {
        "centers": np.asarray(p["centers"], dtype=float),
        "edges": np.asarray(p["edges"], dtype=float),
        "values": np.asarray(p["values"], dtype=float),
        "errors": np.asarray(p["errors"], dtype=float),
        "raw": np.asarray(p["raw"], dtype=float),
        "raw_errors": np.asarray(p["raw_errors"], dtype=float),
        "source_root": payload["source_root"],
        "source_integral_0to1": float(p.get("source_integral_0to1", 0.0)),
    }


def crop_panel(spec: PanelSpec) -> Path:
    image = Image.open(spec.page_image).convert("RGB")
    crop = image.crop(spec.crop_box)
    out = OUTDIR / f"ian_v4_{spec.tag}_panel_crop.png"
    crop.save(out)
    return out


def rel_axis(spec: PanelSpec) -> dict[str, float]:
    x0, y0, _, _ = spec.crop_box
    return {
        "x_left": spec.x_left_abs - x0,
        "x_right": spec.x_right_abs - x0,
        "y_top": spec.y_top_abs - y0,
        "y_bottom": spec.y_bottom_abs - y0,
    }


def data_to_pixel(spec: PanelSpec, x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    ax = rel_axis(spec)
    px = ax["x_left"] + x * (ax["x_right"] - ax["x_left"])
    py = ax["y_bottom"] - (y / spec.y_max) * (ax["y_bottom"] - ax["y_top"])
    return px, py


def pixel_to_data(spec: PanelSpec, points: list[tuple[str, float, float]]) -> tuple[np.ndarray, np.ndarray]:
    ax = rel_axis(spec)
    xs = np.asarray([(px - ax["x_left"]) / (ax["x_right"] - ax["x_left"]) for _, px, _ in points], dtype=float)
    ys = np.asarray(
        [(ax["y_bottom"] - py) / (ax["y_bottom"] - ax["y_top"]) * spec.y_max for _, _, py in points],
        dtype=float,
    )
    return xs, ys


def detect_black_marker_pixels(spec: PanelSpec, crop: Path, sdcc: dict[str, np.ndarray]) -> list[tuple[str, float, float]]:
    arr = np.asarray(Image.open(crop).convert("RGB"))
    black = (arr[..., 0] < 85) & (arr[..., 1] < 85) & (arr[..., 2] < 85)
    x_pred, y_pred = data_to_pixel(spec, sdcc["centers"], sdcc["values"])
    ax = rel_axis(spec)
    points: list[tuple[str, float, float]] = []

    for xp, yp in zip(x_pred, y_pred):
        xlo = max(int(round(xp)) - 7, int(ax["x_left"]) + 1)
        xhi = min(int(round(xp)) + 7, int(ax["x_right"]) - 1)
        ylo = max(int(round(yp)) - 14, int(ax["y_top"]) + 3)
        yhi = min(int(round(yp)) + 14, int(ax["y_bottom"]) - 4)
        best: tuple[float, float, float] | None = None
        for cy in range(ylo, yhi + 1):
            for cx in range(xlo, xhi + 1):
                y0, y1 = max(cy - 4, 0), min(cy + 5, black.shape[0])
                x0, x1 = max(cx - 4, 0), min(cx + 5, black.shape[1])
                central = float(black[y0:y1, x0:x1].sum())
                y2a, y2b = max(cy - 7, 0), min(cy + 8, black.shape[0])
                x2a, x2b = max(cx - 7, 0), min(cx + 8, black.shape[1])
                wide = float(black[y2a:y2b, x2a:x2b].sum())
                # Penalize the bottom frame: a horizontal axis line is dark but
                # has little compact filled-marker structure above it.
                frame_penalty = 8.0 if cy > ax["y_bottom"] - 7 else 0.0
                score = central + 0.12 * wide - 0.20 * abs(cx - xp) - 0.04 * abs(cy - yp) - frame_penalty
                if best is None or score > best[0]:
                    best = (score, float(cx), float(cy))
        if best is None or best[0] < 3.0:
            points.append(("data", float(xp), float(yp)))
        else:
            points.append(("data", best[1], best[2]))
    return points


def datathief_reference_points(spec: PanelSpec, crop: Path, points: list[tuple[str, float, float]]) -> tuple[np.ndarray, np.ndarray, Path, Path]:
    ax = rel_axis(spec)
    out_csv = OUTDIR / f"ppg12_ian_datathief_export_black_data_vs_sdcc_{spec.tag}_points.csv"

    chosen_raw_csv: Path | None = None
    chosen_dt_x: np.ndarray | None = None
    chosen_dt_y: np.ndarray | None = None
    failures: list[str] = []
    for attempt, scale in enumerate((1.0, 1.0, 1.0, 4.0, 5.0, 10.0, 20.0, 100.0, 2.0), start=1):
        refs = [
            (0, ax["x_left"], ax["y_bottom"], 0.0, 0.0),
            (1, ax["x_right"], ax["y_bottom"], 1.0, 0.0),
            (2, ax["x_left"], ax["y_top"], 0.0, spec.y_max * scale),
        ]
        raw_csv = (
            OUTDIR
            / f"ppg12_ian_datathief_export_black_data_vs_sdcc_{spec.tag}"
            f"_raw_y{int(spec.y_max*100):02d}_scale{scale:g}_try{attempt}.csv"
        )
        run_datathief_export(crop, refs, points, raw_csv)
        dt = load_datathief_csv(raw_csv)["data"]
        dt_x = np.asarray([p[0] for p in dt], dtype=float)
        dt_y = np.asarray([p[1] / scale for p in dt], dtype=float)
        finite = np.isfinite(dt_x) & np.isfinite(dt_y)
        valid = (
            int(finite.sum()) == len(points)
            and float(np.nanmax(dt_x) - np.nanmin(dt_x)) > 0.75
            and float(np.nanmax(dt_y)) > max(0.005, 0.02 * spec.y_max)
            and float(np.nanmin(dt_x)) > -0.05
            and float(np.nanmax(dt_x)) < 1.05
            and float(np.nanmin(dt_y)) > -0.02 * spec.y_max
            and float(np.nanmax(dt_y)) < 1.15 * spec.y_max
        )
        if valid:
            chosen_raw_csv = raw_csv
            chosen_dt_x = dt_x
            chosen_dt_y = dt_y
            break
        failures.append(
            f"attempt={attempt} scale={scale:g} finite={int(finite.sum())}/{len(points)} "
            f"x_range=({np.nanmin(dt_x):.4g},{np.nanmax(dt_x):.4g}) "
            f"y_range=({np.nanmin(dt_y):.4g},{np.nanmax(dt_y):.4g})"
        )

    if chosen_raw_csv is None or chosen_dt_x is None or chosen_dt_y is None:
        raise RuntimeError(f"DataThief export failed validation for {spec.tag}: {'; '.join(failures)}")

    with out_csv.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["series", "point_index", "datathief_x", "datathief_y"])
        for i, (xv, yv) in enumerate(zip(chosen_dt_x, chosen_dt_y)):
            writer.writerow(["data", i, f"{xv:.12g}", f"{yv:.12g}"])
    return chosen_dt_x, chosen_dt_y, chosen_raw_csv, out_csv


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "dejavuserif",
            "axes.linewidth": 1.15,
            "xtick.major.size": 6,
            "ytick.major.size": 6,
            "xtick.minor.size": 3,
            "ytick.minor.size": 3,
            "xtick.direction": "in",
            "ytick.direction": "in",
        }
    )


def sphinx_label(ax, x=0.05, y=0.93, fs=17) -> None:
    ax.text(x, y, "sPHENIX", transform=ax.transAxes, ha="left", va="top", fontsize=fs, fontweight="bold", fontstyle="italic")
    ax.text(x + 0.215, y, "Internal", transform=ax.transAxes, ha="left", va="top", fontsize=fs)


def plot_overlay(
    *,
    spec: PanelSpec,
    out: Path,
    ref_x: np.ndarray,
    ref_y: np.ndarray,
    ref_err: np.ndarray,
    cmp_x: np.ndarray,
    cmp_y: np.ndarray,
    cmp_err: np.ndarray | None,
    cmp_label: str,
    ratio_label: str,
    ratio_ylim: tuple[float, float],
    stable_threshold: float,
    manifest: dict[str, object],
    stats_text: str | None = None,
    annotation_label: str | None = None,
) -> None:
    setup_style()
    fig, (ax, rax) = plt.subplots(
        2,
        1,
        figsize=(7.72, 9.98),
        dpi=100,
        sharex=True,
        gridspec_kw={"height_ratios": [3.35, 1.0], "hspace": 0.05},
    )
    ax.errorbar(ref_x, ref_y, yerr=ref_err, fmt="o", color="black", ms=4.2, lw=1.0, label="PPG12 SDCC ROOT data", zorder=3)
    is_datathief = "DataThief" in cmp_label
    ax.errorbar(
        cmp_x,
        cmp_y,
        yerr=cmp_err,
        fmt="s",
        color="#d62728" if is_datathief else "#1f77b4",
        markerfacecolor="none",
        markeredgewidth=1.2,
        ms=4.6,
        lw=1.0,
        label=cmp_label,
        zorder=4,
    )
    sphinx_label(ax)
    ax.text(0.05, 0.84, r"$p{+}p\ \sqrt{s}=200$ GeV", transform=ax.transAxes, fontsize=15, ha="left")
    ax.text(0.05, 0.77, r"$|\eta^\gamma| < 0.7$", transform=ax.transAxes, fontsize=15, ha="left")
    ax.text(0.05, 0.70, annotation_label or spec.ian_label, transform=ax.transAxes, fontsize=13, ha="left")
    ax.set_ylabel("normalized counts", fontsize=17)
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, spec.y_max)
    ax.legend(
        loc="upper right",
        frameon=False,
        fontsize=16.0,
        handlelength=1.5,
        borderaxespad=0.35,
        labelspacing=0.55,
    )
    if stats_text:
        ax.text(
            0.590,
            0.790,
            stats_text,
            transform=ax.transAxes,
            fontsize=15.0,
            va="top",
            ha="left",
            linespacing=1.55,
        )
    ax.tick_params(labelsize=14, top=True, right=True)
    ax.minorticks_on()

    ref_interp = np.interp(cmp_x, ref_x, ref_y, left=np.nan, right=np.nan)
    ratio = np.divide(cmp_y, ref_interp, out=np.full_like(cmp_y, np.nan), where=ref_interp > 0)
    ratio_err = None
    if cmp_err is not None:
        ratio_err = np.divide(cmp_err, ref_interp, out=np.full_like(cmp_err, np.nan), where=ref_interp > 0)
    plot_mask = np.isfinite(ratio) & (ref_interp > 0)
    stable = np.isfinite(ratio) & (ref_interp > stable_threshold)
    if ratio_err is None:
        rax.plot(cmp_x[plot_mask], ratio[plot_mask], "o", color="#d62728" if "DataThief" in cmp_label else "#1f77b4", ms=4.2)
    else:
        rax.errorbar(cmp_x[plot_mask], ratio[plot_mask], yerr=ratio_err[plot_mask], fmt="o", color="#1f77b4", ms=4.0, lw=1.0)
    rax.axhline(1.0, color="0.35", ls="--", lw=1.0)
    rax.set_ylabel(ratio_label, fontsize=15)
    rax.set_xlabel("bdt", fontsize=17)
    rax.set_ylim(*ratio_ylim)
    rax.tick_params(labelsize=14, top=True, right=True)
    rax.minorticks_on()
    fig.subplots_adjust(left=0.115, right=0.985, top=0.985, bottom=0.075, hspace=0.05)
    fig.savefig(out, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)

    stable_ratio = ratio[stable]
    manifest.update(
        {
            "artifact": str(out),
            "ratio_plot_bin_count": int(plot_mask.sum()),
            "stable_threshold_on_reference": stable_threshold,
            "stable_ratio_count": int(stable.sum()),
            "stable_mean_ratio": float(np.nanmean(stable_ratio)) if stable_ratio.size else None,
            "stable_mean_abs_ratio_minus_one": float(np.nanmean(np.abs(stable_ratio - 1.0))) if stable_ratio.size else None,
            "stable_max_abs_ratio_minus_one": float(np.nanmax(np.abs(stable_ratio - 1.0))) if stable_ratio.size else None,
        }
    )
    out.with_suffix(".manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True))


def rebin_counts_to_edges(src_values: np.ndarray, src_edges: np.ndarray, dst_edges: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    out = np.zeros(len(dst_edges) - 1, dtype=float)
    err2 = np.zeros_like(out)
    for i, count in enumerate(src_values):
        lo = float(src_edges[i])
        hi = float(src_edges[i + 1])
        if hi <= dst_edges[0] or lo >= dst_edges[-1]:
            continue
        width = hi - lo
        if width <= 0:
            continue
        for j in range(len(out)):
            ov = max(0.0, min(hi, float(dst_edges[j + 1])) - max(lo, float(dst_edges[j])))
            if ov <= 0:
                continue
            frac = ov / width
            out[j] += count * frac
            err2[j] += count * frac * frac
    return out, np.sqrt(err2)


def infer_count_from_normalized_errors(values: np.ndarray, errors: np.ndarray) -> float:
    counts = [
        (float(y) / float(err)) ** 2
        for y, err in zip(values, errors)
        if y > 0 and err > 0
    ]
    return float(sum(counts))


def current_projection_exact(root_path: Path, hist_name: str, dst_edges: np.ndarray) -> dict[str, np.ndarray | float]:
    import ROOT

    root_file = ROOT.TFile.Open(str(root_path))
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"Could not open current ROOT: {root_path}")
    hist = root_file.Get(hist_name)
    if not hist:
        root_file.Close()
        raise RuntimeError(f"Missing exact current histogram: {hist_name}")

    nbins = hist.GetNbinsX()
    counts = np.asarray([float(hist.GetBinContent(i)) for i in range(1, nbins + 1)], dtype=float)
    edges = np.asarray(
        [float(hist.GetXaxis().GetBinLowEdge(i)) for i in range(1, nbins + 1)]
        + [float(hist.GetXaxis().GetBinUpEdge(nbins))],
        dtype=float,
    )
    underflow = float(hist.GetBinContent(0))
    overflow = float(hist.GetBinContent(nbins + 1))
    root_file.Close()

    raw, rawerr = rebin_counts_to_edges(counts, edges, dst_edges)
    total = float(np.sum(raw))
    vals = raw / total if total > 0 else raw
    errs = rawerr / total if total > 0 else rawerr
    centers = 0.5 * (dst_edges[:-1] + dst_edges[1:])
    return {
        "centers": centers,
        "values": vals,
        "errors": errs,
        "raw": raw,
        "raw_errors": rawerr,
        "total_0to1": total,
        "source_integrals": {hist_name: float(np.sum(counts))},
        "source_edges_low": float(edges[0]),
        "source_edges_high": float(edges[-1]),
        "underflow": underflow,
        "overflow": overflow,
    }


def current_projection(spec: PanelSpec, dst_edges: np.ndarray) -> dict[str, np.ndarray | float]:
    import uproot

    f = uproot.open(CURRENT_FULL_PP)
    counts = None
    edges = None
    integral_inputs: dict[str, float] = {}
    for hname in spec.current_hists:
        key = f"{CURRENT_DIR}/{hname}"
        vals, hedges = f[key].to_numpy(flow=False)
        integral_inputs[hname] = float(np.sum(vals))
        counts = vals.astype(float).copy() if counts is None else counts + vals
        edges = hedges.astype(float)
    if counts is None or edges is None:
        raise RuntimeError(f"No current histograms found for {spec.tag}")
    raw, rawerr = rebin_counts_to_edges(counts, edges, dst_edges)
    total = float(np.sum(raw))
    vals = raw / total if total > 0 else raw
    errs = rawerr / total if total > 0 else rawerr
    centers = 0.5 * (dst_edges[:-1] + dst_edges[1:])
    return {
        "centers": centers,
        "values": vals,
        "errors": errs,
        "raw": raw,
        "raw_errors": rawerr,
        "total_0to1": total,
        "source_integrals": integral_inputs,
        "source_edges_low": float(edges[0]),
        "source_edges_high": float(edges[-1]),
    }


def write_reference_validation(spec: PanelSpec, sdcc: dict[str, np.ndarray]) -> None:
    crop = crop_panel(spec)
    points = detect_black_marker_pixels(spec, crop, sdcc)
    dt_x, dt_y, raw_csv, out_csv = datathief_reference_points(spec, crop, points)
    out = OUTDIR / f"ppg12_ian_datathief_export_black_data_vs_sdcc_{spec.tag}_overlay_slidefit_772x998.png"
    manifest = {
        "comparison": "PPG12 IAN rendered black data markers exported through DataThief vs PPG12 SDCC ROOT projection",
        "panel": spec.title,
        "ian_page": str(spec.page_image),
        "ian_panel_crop": str(crop),
        "ppg12_sdcc_json": str(SDCC_JSON),
        "ppg12_sdcc_root": sdcc["source_root"],
        "ppg12_sdcc_hist": spec.sdcc_hist,
        "datathief_jar": str(DATATHIEF_JAR),
        "datathief_jar_md5": jar_md5(),
        "datathief_raw_csv": str(raw_csv),
        "datathief_scaled_csv": str(out_csv),
        "axis_calibration": rel_axis(spec),
        "y_max": spec.y_max,
    }
    plot_overlay(
        spec=spec,
        out=out,
        ref_x=sdcc["centers"],
        ref_y=sdcc["values"],
        ref_err=sdcc["errors"],
        cmp_x=dt_x,
        cmp_y=dt_y,
        cmp_err=None,
        cmp_label="DataThief export from IAN PNG",
        ratio_label="PNG / SDCC",
        ratio_ylim=spec.ratio_ylim_reference,
        stable_threshold=max(0.0025, 0.01 * spec.y_max),
        manifest=manifest,
    )


def write_current_overlay(spec: PanelSpec, sdcc: dict[str, np.ndarray]) -> None:
    cur = current_projection(spec, sdcc["edges"])
    out = OUTDIR / f"ppg12_sdcc_vs_current_default_fullpp_{spec.tag}_data_overlay_slidefit_772x998.png"
    manifest = {
        "comparison": "PPG12 SDCC ROOT data projection vs current completed full-pp PhotonClusterBuilder/RecoilJets BDT-score QA",
        "panel": spec.title,
        "ppg12_sdcc_json": str(SDCC_JSON),
        "ppg12_sdcc_root": sdcc["source_root"],
        "ppg12_sdcc_hist": spec.sdcc_hist,
        "current_full_pp_root": str(CURRENT_FULL_PP),
        "current_directory": CURRENT_DIR,
        "current_histograms": list(spec.current_hists),
        "current_source_integrals": cur["source_integrals"],
        "current_rebinned_integral_0to1": cur["total_0to1"],
        "current_source_axis": [cur["source_edges_low"], cur["source_edges_high"]],
        "important_caveat": (
            "Current overlay uses the canonical PPG12 pp-data trigger namespace "
            "PPG12_scaledtrigger30, filled from direct GL1 ScaledVector bit 30. "
            "The similarly named TriggerAnalyzer directory is diagnostic only."
        ),
    }
    plot_overlay(
        spec=spec,
        out=out,
        ref_x=sdcc["centers"],
        ref_y=sdcc["values"],
        ref_err=sdcc["errors"],
        cmp_x=cur["centers"],
        cmp_y=cur["values"],
        cmp_err=cur["errors"],
        cmp_label=spec.current_label,
        ratio_label="Current / PPG12",
        ratio_ylim=spec.ratio_ylim_current,
        stable_threshold=max(0.0025, 0.01 * spec.y_max),
        manifest=manifest,
    )


def write_exact_current_overlay(args: argparse.Namespace) -> None:
    spec = next(s for s in PANELS if s.tag == args.panel)
    sdcc = load_sdcc_projection(spec)
    cur = current_projection_exact(args.current_root, args.current_hist, sdcc["edges"])

    ref_x = sdcc["centers"]
    ref_y = sdcc["values"]
    ref_err = sdcc["errors"]
    cur_x = cur["centers"]
    cur_y = cur["values"]
    cur_err = cur["errors"]
    if len(cur_x) != len(ref_x) or not np.allclose(cur_x, ref_x, atol=1e-6):
        raise RuntimeError("Current and PPG12 BDT bin centers do not match")

    ref_interp = np.interp(cur_x, ref_x, ref_y, left=np.nan, right=np.nan)
    ratio = np.divide(cur_y, ref_interp, out=np.full_like(cur_y, np.nan), where=ref_interp > 0)
    stable_threshold = max(0.0025, 0.01 * spec.y_max)
    stable = np.isfinite(ratio) & (ref_interp > stable_threshold)
    stable_ratio = ratio[stable]
    max_abs_ratio_minus_one = float(np.nanmax(np.abs(stable_ratio - 1.0))) if stable_ratio.size else float("nan")
    max_ratio_x = float(cur_x[stable][int(np.nanargmax(np.abs(stable_ratio - 1.0)))]) if stable_ratio.size else float("nan")
    ppg12_entries = infer_count_from_normalized_errors(ref_y, ref_err)
    stats_text = (
        f"PPG12 SDCC N_eff={ppg12_entries:.0f}\n"
        f"July 1 final pp N={cur['total_0to1']:.0f}\n"
        rf"max $|R-1|$ = {100.0 * max_abs_ratio_minus_one:.1f}%"
        f"\nnear bdt = {max_ratio_x:.2f}"
    )

    manifest = {
        "comparison": "PPG12 SDCC ROOT data projection vs July 1 final combined pp TableQA BDT object",
        "panel": spec.title,
        "ppg12_sdcc_json": str(SDCC_JSON),
        "ppg12_sdcc_root": sdcc["source_root"],
        "ppg12_sdcc_hist": spec.sdcc_hist,
        "ppg12_sdcc_effective_entries_from_errors": ppg12_entries,
        "current_root": str(args.current_root),
        "current_directory": CURRENT_DIR,
        "current_histograms": [args.current_hist],
        "current_source_integrals": cur["source_integrals"],
        "current_rebinned_integral_0to1": cur["total_0to1"],
        "current_underflow": cur["underflow"],
        "current_overflow": cur["overflow"],
        "current_source_axis": [cur["source_edges_low"], cur["source_edges_high"]],
        "stable_threshold_on_reference": stable_threshold,
        "stable_ratio_count": int(stable.sum()),
        "stable_mean_ratio": float(np.nanmean(stable_ratio)) if stable_ratio.size else None,
        "stable_mean_abs_ratio_minus_one": float(np.nanmean(np.abs(stable_ratio - 1.0))) if stable_ratio.size else None,
        "stable_max_abs_ratio_minus_one": max_abs_ratio_minus_one,
        "stable_max_abs_ratio_minus_one_bdt_center": max_ratio_x,
        "note": "Current points are from the final hierarchical July 1 combined pp ROOT, using the exact PPG12 TableQA 22<pT<28 GeV no-NPB BDT object, rebinned to the PPG12 Fig.13/19 50-bin 0-1 grid.",
    }
    plot_overlay(
        spec=spec,
        out=args.output,
        ref_x=ref_x,
        ref_y=ref_y,
        ref_err=ref_err,
        cmp_x=cur_x,
        cmp_y=cur_y,
        cmp_err=cur_err,
        cmp_label=args.current_legend,
        ratio_label="Current / PPG12",
        ratio_ylim=args.ratio_ylim,
        stable_threshold=stable_threshold,
        manifest=manifest,
        stats_text=stats_text,
        annotation_label="22 < pT < 28 GeV, no NPB cut" if spec.tag == "bdt_no_npb_22_28" else None,
    )


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--panel", choices=[p.tag for p in PANELS], default="bdt_no_npb_22_28")
    ap.add_argument("--current-root", type=Path, default=CURRENT_FULL_PP)
    ap.add_argument("--current-hist", default=None)
    ap.add_argument("--current-legend", default="July 1 final pp data")
    ap.add_argument("--output", type=Path, default=OUTDIR / "ppg12_sdcc_vs_current_default_fullpp_bdt_no_npb_22_28_data_overlay_slidefit_772x998.png")
    ap.add_argument("--skip-reference", action="store_true")
    ap.add_argument("--ratio-ymin", type=float, default=0.0)
    ap.add_argument("--ratio-ymax", type=float, default=4.5)
    args = ap.parse_args()
    args.ratio_ylim = (args.ratio_ymin, args.ratio_ymax)
    return args


def main() -> None:
    args = parse_args()
    OUTDIR.mkdir(parents=True, exist_ok=True)
    if not SDCC_JSON.exists():
        write_sdcc_json_from_raw()
    else:
        # Keep JSON synchronized with the raw SDCC readback if present.
        write_sdcc_json_from_raw()
    if args.current_hist:
        write_exact_current_overlay(args)
        return
    for spec in PANELS:
        sdcc = load_sdcc_projection(spec)
        if not args.skip_reference:
            write_reference_validation(spec, sdcc)
        write_current_overlay(spec, sdcc)


if __name__ == "__main__":
    main()
