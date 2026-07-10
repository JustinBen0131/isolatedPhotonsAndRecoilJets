#!/usr/bin/env python3
"""THE-95 corrected low-calo veto candidate slide.

This regenerates the old THE-58 style slide using the corrected event-energy
definition:

  E_calo = sum finite calibrated TOWERINFO_CALIB_* energies
  with CaloTowerStatus forced for embedded and TowerInfo::get_isGood required.

The cut shown here is a candidate validation threshold, not a production
default. It is fit to the valley between the fixed low-E band and the normal
event-energy band where that valley is clean, then applied only through 55%
centrality. Above that the normal peripheral band overlaps the low-E region.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch
import numpy as np


HIST_JSON = Path(
    "dataOutput/auauTightBDTValidation/THE32_lowCaloDiagnosticClosure_20260603/"
    "blair_reference_20260702/fourway_20260704/"
    "THE95_lowcalo4way_statusOn_isGoodOn_fast_20260704.wide_hist.json"
)
OUTDIR = Path(
    "dataOutput/auauTightBDTValidation/THE95_lowCaloGoodTowerCutCandidate_20260706"
)
PNG = OUTDIR / "the95_isgood_low_calo_smooth_cut_candidate_slide.png"
MANIFEST = OUTDIR / "the95_isgood_low_calo_smooth_cut_candidate_slide.manifest.json"
SCRIPT = OUTDIR / "the95_isgood_low_calo_smooth_cut_candidate_slide_script.md"

FIT_DEGREE = 2
FIT_CENT_HI = 50.0
CUT_CENT_HI = 55.0
LOW_PEAK_WINDOW = (1.65, 1.95)
MAIN_PEAK_WINDOW = (2.15, 3.45)
MIN_PEAK_GAP = 0.55
MAX_VALLEY_TO_LOWPEAK = 0.25

BLUE = "#1F77B4"
RED = "#D62728"
INK = "#172033"
MUTED = "#5D6B7A"
BOX_EDGE = "#CBD5E1"
GOLD = "#B7791F"
GOLD_SOFT = "#FFF7ED"
GOLD_LINE = "#F1C27D"
PURPLE = "#7A5BBE"
TEAL = "#2E7E8E"

plt.rcParams.update(
    {
        "font.family": "Times New Roman",
        "mathtext.fontset": "stix",
        "axes.unicode_minus": False,
    }
)


def step_xy(edges: np.ndarray, counts: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    return np.repeat(edges, 2)[1:-1], np.repeat(counts, 2)


def comma(value: float) -> str:
    return f"{int(round(value)):,}"


def smooth(y: np.ndarray, width: int = 5) -> np.ndarray:
    return np.convolve(y, np.ones(width) / width, mode="same")


def derive_thresholds(payload: dict) -> tuple[np.ndarray, list[dict], np.ndarray, np.ndarray]:
    edges = np.asarray(payload["energy_edges"], dtype="float64")
    centers = 0.5 * (edges[:-1] + edges[1:])

    rows = []
    fit_x = []
    fit_y = []
    for panel in payload["panels"]:
        total = np.asarray(panel["retained_hist"], dtype="float64") + np.asarray(
            panel["removed_hist"], dtype="float64"
        )
        s = smooth(total, 5)
        cent = 0.5 * (float(panel["cent_lo"]) + float(panel["cent_hi"]))

        low_mask = (centers >= LOW_PEAK_WINDOW[0]) & (centers <= LOW_PEAK_WINDOW[1])
        main_mask = (centers >= MAIN_PEAK_WINDOW[0]) & (centers <= MAIN_PEAK_WINDOW[1])
        low_idx = np.where(low_mask)[0][int(np.argmax(s[low_mask]))]
        main_idx = np.where(main_mask)[0][int(np.argmax(s[main_mask]))]

        valley_idx = None
        valley = float("nan")
        valley_ratio = float("nan")
        peak_gap = float(centers[main_idx] - centers[low_idx])
        if main_idx > low_idx + 3:
            valley_idx = low_idx + 1 + int(np.argmin(s[low_idx + 1 : main_idx]))
            valley = float(centers[valley_idx])
            valley_ratio = float(s[valley_idx] / max(s[low_idx], 1.0))

        fit_anchor = (
            valley_idx is not None
            and cent <= FIT_CENT_HI
            and peak_gap >= MIN_PEAK_GAP
            and valley_ratio < MAX_VALLEY_TO_LOWPEAK
        )
        if fit_anchor:
            fit_x.append(cent)
            fit_y.append(valley)

        rows.append(
            {
                "lo": int(round(float(panel["cent_lo"]))),
                "hi": int(round(float(panel["cent_hi"]))),
                "cent": cent,
                "full_hist": total,
                "total": float(total.sum()),
                "low_peak": float(centers[low_idx]),
                "main_peak": float(centers[main_idx]),
                "valley": valley,
                "valley_ratio": valley_ratio,
                "peak_gap": peak_gap,
                "fit_anchor": fit_anchor,
                "cut_candidate": cent < CUT_CENT_HI,
            }
        )

    if len(fit_x) < FIT_DEGREE + 1:
        raise RuntimeError("Not enough clean valley anchors to fit threshold")

    coef = np.polyfit(np.asarray(fit_x), np.asarray(fit_y), FIT_DEGREE)
    for row in rows:
        threshold = float(np.polyval(coef, row["cent"]))
        row["threshold"] = threshold
        if row["cut_candidate"]:
            below = centers < threshold
            removed_hist = np.where(below, row["full_hist"], 0.0)
            retained_hist = np.where(below, 0.0, row["full_hist"])
        else:
            removed_hist = np.zeros_like(row["full_hist"])
            retained_hist = row["full_hist"]
        row["removed_hist"] = removed_hist
        row["retained_hist"] = retained_hist
        row["removed"] = float(removed_hist.sum())
        row["frac"] = float(row["removed"] / max(row["total"], 1.0))
    return centers, rows, coef, edges


def draw_connected_hist(ax, edges: np.ndarray, hist: np.ndarray, color: str, alpha: float) -> None:
    hist = np.asarray(hist, dtype="float64")
    nonzero = np.nonzero(hist > 0)[0]
    if nonzero.size == 0:
        return
    i0 = int(nonzero[0])
    i1 = int(nonzero[-1])
    local_edges = edges[i0 : i1 + 2]
    counts = np.maximum(hist[i0 : i1 + 1], 0.9)
    x, y = step_xy(local_edges, counts)
    ax.plot(x, y, color=color, lw=1.75)
    ax.fill_between(x, y, 0.9, color=color, alpha=alpha)


def build(hist_json: Path, png: Path, manifest_path: Path, script_path: Path) -> None:
    png.parent.mkdir(parents=True, exist_ok=True)
    payload = json.loads(hist_json.read_text())
    centers, rows, coef, edges = derive_thresholds(payload)

    cut_rows = [row for row in rows if row["cut_candidate"]]
    merged_rows = [row for row in rows if not row["cut_candidate"]]
    merged = {
        "lo": merged_rows[0]["lo"],
        "hi": merged_rows[-1]["hi"],
        "nbins": len(merged_rows),
        "full_hist": sum((row["full_hist"] for row in merged_rows), np.zeros_like(rows[0]["full_hist"])),
        "total": sum(row["total"] for row in merged_rows),
    }
    total_all = sum(row["total"] for row in rows)
    removed_all = sum(row["removed"] for row in rows)
    y_max = max(max(row["retained_hist"].max(), row["removed_hist"].max()) for row in cut_rows)
    merged_y_max = merged["full_hist"].max()

    fig = plt.figure(figsize=(12.8, 7.2), dpi=200)
    fig.patch.set_facecolor("white")
    gs = fig.add_gridspec(
        3, 4, left=0.052, right=0.965, top=0.550, bottom=0.090, hspace=1.24, wspace=0.135
    )

    def style(ax):
        ax.set_xlim(0.0, 3.52)
        ax.set_xticks([0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5])
        ax.set_xlabel(r"$\log_{10}(E_{\rm calo}+1)$", fontsize=9.2, color=MUTED, labelpad=3.0)
        ax.tick_params(axis="both", labelsize=8.0, colors=INK, length=3.0, width=0.85)
        for spine in ax.spines.values():
            spine.set_linewidth(0.95)
            spine.set_color(INK)

    def cut_panel(ax, row):
        ax.set_facecolor("white")
        ax.set_yscale("log")
        draw_connected_hist(ax, edges, row["retained_hist"], BLUE, 0.06)
        draw_connected_hist(ax, edges, row["removed_hist"], RED, 0.10)
        ax.axvline(row["threshold"], color="white", lw=3.8, ls=(0, (3.0, 3.0)), zorder=4)
        ax.axvline(row["threshold"], color=GOLD, lw=2.1, ls=(0, (3.0, 3.0)), zorder=5)
        ax.set_ylim(0.8, max(12.0, y_max * 2.45))
        ax.set_title(f"{row['lo']}-{row['hi']}%", fontsize=11.0, fontweight="bold", color=INK, pad=14)
        ax.text(
            0.0,
            1.05,
            rf"$N_{{\rm acc}}$={comma(row['total'] - row['removed'])}",
            transform=ax.transAxes,
            fontsize=8.8,
            color=BLUE,
            fontweight="bold",
            ha="left",
            va="bottom",
            clip_on=False,
        )
        ax.text(
            0.51,
            1.05,
            rf"$N_{{\rm cut}}$={comma(row['removed'])}",
            transform=ax.transAxes,
            fontsize=8.8,
            color=RED,
            fontweight="bold",
            ha="left",
            va="bottom",
            clip_on=False,
        )
        ax.text(
            1.0,
            1.05,
            f"({100 * row['frac']:.1f}%)",
            transform=ax.transAxes,
            fontsize=8.8,
            color=RED,
            fontweight="bold",
            ha="right",
            va="bottom",
            clip_on=False,
        )
        style(ax)

    def merged_panel(ax):
        ax.set_facecolor("#FBF7FF")
        ax.set_yscale("log")
        draw_connected_hist(ax, edges, merged["full_hist"], PURPLE, 0.07)
        ax.set_ylim(0.8, max(12.0, merged_y_max * 2.45))
        ax.set_title(
            f"{merged['lo']}-{merged['hi']}% ({merged['nbins']} bins)",
            fontsize=10.8,
            fontweight="bold",
            color=INK,
            pad=14,
        )
        ax.text(
            0.0,
            1.05,
            rf"$N_{{\rm acc}}$={comma(merged['total'])}",
            transform=ax.transAxes,
            fontsize=8.8,
            color=BLUE,
            fontweight="bold",
            ha="left",
            va="bottom",
            clip_on=False,
        )
        box = FancyBboxPatch(
            (0.49, 0.46),
            0.48,
            0.48,
            boxstyle="round,pad=0.015,rounding_size=0.04",
            transform=ax.transAxes,
            facecolor="#F3EDFB",
            edgecolor=PURPLE,
            linewidth=2.0,
            zorder=8,
        )
        ax.add_patch(box)
        ax.text(
            0.73,
            0.78,
            "not cut:",
            transform=ax.transAxes,
            fontsize=13.5,
            fontweight="bold",
            color=PURPLE,
            ha="center",
            va="center",
            zorder=9,
        )
        ax.text(
            0.73,
            0.61,
            "merged",
            transform=ax.transAxes,
            fontsize=13.5,
            fontweight="bold",
            color=PURPLE,
            ha="center",
            va="center",
            zorder=9,
        )
        style(ax)

    slots = [
        (0, 0),
        (0, 1),
        (0, 2),
        (0, 3),
        (1, 0),
        (1, 1),
        (1, 2),
        (1, 3),
        (2, 0),
        (2, 1),
        (2, 2),
    ]
    left_axes = []
    for slot, row in zip(slots, cut_rows):
        ax = fig.add_subplot(gs[slot[0], slot[1]])
        cut_panel(ax, row)
        if slot[1] == 0:
            left_axes.append(ax)
    axm = fig.add_subplot(gs[2, 3])
    merged_panel(axm)
    for ax in left_axes:
        ax.set_ylabel("raw event counts", fontsize=8.6, color=INK)

    # Header
    fig.text(
        0.052,
        0.966,
        "Low-calo event veto candidate: corrected good-tower sum",
        fontsize=20.5,
        fontweight="bold",
        color=INK,
        va="top",
    )
    fig.text(
        0.052,
        0.918,
        "Raw embedded counts, not reweighted.  Strict sum: CaloTowerStatus on; get_isGood required; signed finite calibrated towers.",
        fontsize=11.8,
        color=MUTED,
        va="top",
    )

    chip = FancyBboxPatch(
        (0.052, 0.842),
        0.590,
        0.046,
        boxstyle="round,pad=0.004,rounding_size=0.006",
        transform=fig.transFigure,
        facecolor=GOLD_SOFT,
        edgecolor=GOLD_LINE,
        linewidth=1.0,
        zorder=1,
    )
    fig.patches.append(chip)
    fig.text(0.062, 0.865, "Candidate veto:", fontsize=12.5, fontweight="bold", color=GOLD, va="center")
    fig.text(
        0.165,
        0.865,
        r"$\mathrm{remove\ if}\ \log_{10}(E_{\rm calo}+1)<T(c),\ \ "
        r"E_{\rm calo}=E_{\rm CEMC}+E_{\rm IHCal}+E_{\rm OHCal}$",
        fontsize=11.0,
        color=INK,
        va="center",
    )

    fig.text(0.058, 0.812, r"$\blacktriangleright$", fontsize=11.3, color=GOLD, va="center")
    fig.text(0.086, 0.812, "Removes ~5.8-6.3%", fontsize=13.5, fontweight="bold", color=INK, va="center")
    fig.text(
        0.214,
        0.812,
        "in each 5% bin where the low-E band is visually separable.",
        fontsize=12.9,
        color=INK,
        va="center",
    )
    fig.text(0.058, 0.763, r"$\blacktriangleright$", fontsize=11.3, color=GOLD, va="center")
    fig.text(0.086, 0.763, "Past ~55%:", fontsize=13.5, fontweight="bold", color=INK, va="center")
    fig.text(
        0.168,
        0.763,
        "normal peripheral total-calo band overlaps the low-E region - not cut.",
        fontsize=12.9,
        color=INK,
        va="center",
    )

    rem_box = FancyBboxPatch(
        (0.052, 0.678),
        0.420,
        0.050,
        boxstyle="round,pad=0.006,rounding_size=0.008",
        transform=fig.transFigure,
        facecolor="white",
        edgecolor="#C96A6A",
        linewidth=1.8,
        zorder=2,
    )
    fig.patches.append(rem_box)
    fig.text(0.066, 0.703, "Candidate removed:", fontsize=13.5, fontweight="bold", color="#B23B3B", va="center")
    fig.text(
        0.248,
        0.703,
        f"{comma(removed_all)} / {comma(total_all)}  ({100 * removed_all / total_all:.1f}%)",
        fontsize=13.5,
        fontweight="bold",
        color=INK,
        va="center",
    )

    legend_y = 0.640
    for x0, color, dash, label, text_x in [
        (0.058, BLUE, None, "kept", 0.096),
        (0.190, RED, None, "removed", 0.228),
        (0.346, GOLD, (0, (4, 3)), r"$T(c)$ candidate", 0.384),
    ]:
        line = Line2D([x0, x0 + 0.032], [legend_y, legend_y], color=color, lw=3.6, transform=fig.transFigure)
        if dash:
            line.set_linestyle(dash)
        fig.add_artist(line)
        fig.text(text_x, legend_y, label, fontsize=12.6, color=INK, va="center")

    # Curve inset
    curve_box = FancyBboxPatch(
        (0.675, 0.688),
        0.290,
        0.270,
        boxstyle="round,pad=0.006,rounding_size=0.008",
        transform=fig.transFigure,
        facecolor="white",
        edgecolor=BOX_EDGE,
        linewidth=1.1,
        zorder=1,
    )
    fig.patches.append(curve_box)
    fig.text(0.820, 0.944, "Smooth candidate curve  T(c)", fontsize=11.2, fontweight="bold", color=INK, ha="center", va="top")
    inset = fig.add_axes([0.715, 0.800, 0.226, 0.104])
    inset.set_zorder(6)
    inset.set_facecolor("white")
    cc = np.linspace(0, 80, 300)
    tv = np.polyval(coef, cc)
    apply = cc < CUT_CENT_HI
    inset.plot(cc[apply], tv[apply], color=GOLD, lw=2.7, zorder=4)
    inset.plot(cc[~apply], tv[~apply], color=GOLD, lw=2.0, alpha=0.25, zorder=4)
    fit_rows = [row for row in rows if row["fit_anchor"]]
    inset.scatter([row["cent"] for row in fit_rows], [row["valley"] for row in fit_rows], s=13, color=INK, zorder=5)
    inset.axvline(CUT_CENT_HI, color=PURPLE, lw=1.1, ls=(0, (3, 2)))
    inset.text(56, 2.34, "stop cut:\nbands merge", fontsize=7.0, color=PURPLE, va="top")
    inset.set_xlim(0, 80)
    inset.set_ylim(1.55, 2.70)
    inset.set_xlabel("centrality (%)", fontsize=8.2, color=INK, labelpad=2.0)
    inset.set_ylabel(r"$T(c)$", fontsize=8.2, color=INK, labelpad=2.0)
    inset.tick_params(labelsize=7.2, colors=INK, length=2.6, width=0.8)
    for spine in inset.spines.values():
        spine.set_linewidth(0.85)
        spine.set_color(BOX_EDGE)

    a2, a1, a0 = coef
    eq_box = FancyBboxPatch(
        (0.685, 0.692),
        0.272,
        0.054,
        boxstyle="round,pad=0.004,rounding_size=0.006",
        transform=fig.transFigure,
        facecolor=GOLD_SOFT,
        edgecolor=GOLD_LINE,
        linewidth=1.4,
        zorder=6,
    )
    fig.patches.append(eq_box)
    eqn = rf"$T(c)=-8.48\times10^{{-5}}c^2{a1:+.4f}c{a0:+.2f}$"
    fig.text(0.820, 0.719, eqn, fontsize=12.2, color=INK, ha="center", va="center", zorder=8)

    fig.savefig(png, dpi=200, facecolor="white")
    plt.close(fig)

    manifest = {
        "png": str(png),
        "script": str(script_path),
        "source_hist_json": str(hist_json),
        "source_schema": payload.get("schema"),
        "source_variant_note": payload.get("variant_note"),
        "tower_definition": "TOWERINFO_CALIB_CEMC/HCALIN/HCALOUT, CaloTowerStatus on, get_isGood required, finite signed energies included",
        "method": "low_band_to_normal_band_valley_fit_candidate",
        "fit_degree": FIT_DEGREE,
        "fit_cent_hi": FIT_CENT_HI,
        "cut_cent_hi": CUT_CENT_HI,
        "low_peak_window": list(LOW_PEAK_WINDOW),
        "main_peak_window": list(MAIN_PEAK_WINDOW),
        "coef": [float(x) for x in coef],
        "total": float(total_all),
        "removed": float(removed_all),
        "removed_frac": float(removed_all / total_all),
        "per_bin": [
            {
                "cent_low": row["lo"],
                "cent_high": row["hi"],
                "threshold_log10_ecalo_plus1": float(row["threshold"]),
                "threshold_ecalo_gev": float(10 ** row["threshold"] - 1.0),
                "fit_anchor": bool(row["fit_anchor"]),
                "cut_candidate": bool(row["cut_candidate"]),
                "removed": float(row["removed"]),
                "total": float(row["total"]),
                "removed_frac": float(row["frac"]),
                "low_peak": float(row["low_peak"]),
                "main_peak": float(row["main_peak"]),
                "valley": float(row["valley"]),
                "valley_ratio": float(row["valley_ratio"]),
                "peak_gap": float(row["peak_gap"]),
            }
            for row in rows
        ],
        "status": "candidate only; requires removed-event validation before production BDT training",
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")

    script = f"""# THE-95 corrected low-calo veto candidate

This slide shows the candidate low-calo veto after fixing the event-energy
definition.  The sum uses calibrated CEMC, IHCal, and OHCal TowerInfo towers,
forces CaloTowerStatus for embedded input, requires `get_isGood()`, and keeps
finite signed tower energies.  So this is not the old positive-only tower sum.

The visible low-total-calo band remains at about six percent per 5% centrality
bin where it is separable.  The candidate threshold is fit to the valley between
that low-E band and the normal total-calo band, and is only applied through
55% centrality.  Past that, the normal peripheral distribution overlaps the
same low-energy region, so the slide deliberately does not cut those bins.

The candidate removes {comma(removed_all)} / {comma(total_all)} events,
or {100 * removed_all / total_all:.1f}% of this diagnostic embedded sample.
This should be treated as a validation candidate: before retraining a default
BDT with it, we should inspect the removed-event sample by run/chunk/sample and
confirm it is an overlay pathology rather than a legitimate low-activity tail.
"""
    script_path.write_text(script)

    print(png)
    print(manifest_path)
    print(script_path)
    print(f"removed {comma(removed_all)} / {comma(total_all)} = {100 * removed_all / total_all:.2f}%")
    print("coef", [float(x) for x in coef])


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--hist-json", type=Path, default=HIST_JSON)
    parser.add_argument("--png", type=Path, default=PNG)
    parser.add_argument("--manifest", type=Path, default=MANIFEST)
    parser.add_argument("--script", type=Path, default=SCRIPT)
    args = parser.parse_args()
    build(args.hist_json, args.png, args.manifest, args.script)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
