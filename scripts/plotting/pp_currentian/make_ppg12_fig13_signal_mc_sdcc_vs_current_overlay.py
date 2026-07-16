#!/usr/bin/env python3
"""Overlay PPG12 Fig. 13 Signal MC BDT template with current photon+jet output.

The PPG12 side is the SDCC-extracted Signal MC curve from the same
four-curve Fig. 13/19 readback used by the existing BDT data-validation plots.
The current side is the photon+jet combined RecoilJets TableQA BDT object for
the same Fig. 13 token: eta0, pt3 (22-28 GeV), cut0 (no NPB cut).
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import uproot


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
PPG12_JSON = (
    REPO
    / "dataOutput/ppg12PhotonYield/ppg12_photon_yield_v1_data_20260620"
    / "shower_shape_reference_validation/fig13_19_bdt"
    / "ppg12_sdcc_fig19_bdt_fourcurve_projections.json"
)
CURRENT_POINTER = REPO / "dataOutput/current_recoiljets_artifacts/current/pp_sim_photonjet_merged/current.json"
DEFAULT_CURRENT_ROOT = (
    REPO
    / "dataOutput/ppg12Parity/the76_ppg12_fig24_photonjet_fix_20260702_014217"
    / "final_roots/photonjet/RecoilJets_photonjet5plus10plus20_MERGED.root"
)
CURRENT_HIST = "SIM/h1d_bdt_eta0_pt3_cut0"
OUTDIR = REPO / "dataOutput/ppg12Parity/the76_ppg12_fig13_signal_mc_overlay_20260706"


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


def sphinx_label(ax, x: float = 0.05, y: float = 0.93, fs: float = 17) -> None:
    ax.text(
        x,
        y,
        "sPHENIX",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=fs,
        fontweight="bold",
        fontstyle="italic",
    )
    ax.text(x + 0.215, y, "Internal", transform=ax.transAxes, ha="left", va="top", fontsize=fs)


def load_current_root() -> tuple[Path, dict[str, object]]:
    if CURRENT_POINTER.exists():
        pointer = json.loads(CURRENT_POINTER.read_text())
        paths = pointer.get("root_paths") or []
        if paths:
            root = Path(paths[0])
            return root, pointer
    return DEFAULT_CURRENT_ROOT, {
        "note": "Fallback path used because current pointer was unavailable.",
        "root_paths": [str(DEFAULT_CURRENT_ROOT)],
    }


def load_ppg12_signal() -> dict[str, object]:
    payload = json.loads(PPG12_JSON.read_text())
    panel = payload["panels"]["bdt_no_npb_22_28"]
    curve = panel["curves"]["signal_mc"]
    return {
        "centers": np.asarray(curve["centers"], dtype=float),
        "edges": np.asarray(curve["edges"], dtype=float),
        "values": np.asarray(curve["values"], dtype=float),
        "errors": np.asarray(curve["errors"], dtype=float),
        "source_root": panel["summary"]["source_files"]["signal"],
        "source_hist": panel["hist_name"],
        "panel_summary": panel["summary"],
    }


def rebin_counts_to_edges(src_values: np.ndarray, src_edges: np.ndarray, dst_edges: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    out = np.zeros(len(dst_edges) - 1, dtype=float)
    err2 = np.zeros_like(out)
    for value, lo, hi in zip(src_values, src_edges[:-1], src_edges[1:]):
        width = float(hi - lo)
        if width <= 0:
            continue
        for j, (dst_lo, dst_hi) in enumerate(zip(dst_edges[:-1], dst_edges[1:])):
            overlap = max(0.0, min(float(hi), float(dst_hi)) - max(float(lo), float(dst_lo)))
            if overlap <= 0:
                continue
            frac = overlap / width
            out[j] += float(value) * frac
            err2[j] += max(float(value), 0.0) * frac * frac
    return out, np.sqrt(err2)


def load_current_signal(current_root: Path, dst_edges: np.ndarray) -> dict[str, object]:
    with uproot.open(current_root) as root_file:
        hist = root_file[CURRENT_HIST]
        values, edges = hist.to_numpy(flow=False)
        entries = float(hist.member("fEntries"))
    raw, raw_err = rebin_counts_to_edges(np.asarray(values, dtype=float), np.asarray(edges, dtype=float), dst_edges)
    total = float(np.sum(raw))
    norm = raw / total if total > 0 else raw
    err = raw_err / total if total > 0 else raw_err
    centers = 0.5 * (dst_edges[:-1] + dst_edges[1:])
    return {
        "centers": centers,
        "edges": dst_edges,
        "values": norm,
        "errors": err,
        "raw": raw,
        "raw_errors": raw_err,
        "raw_integral_0to1": total,
        "entries": entries,
    }


def write_outputs(ppg12: dict[str, object], current: dict[str, object], current_root: Path, pointer: dict[str, object]) -> tuple[Path, Path, Path]:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    out_png = OUTDIR / "ppg12_fig13_signal_mc_sdcc_vs_current_photonjet_overlay_ratio.png"
    out_csv = OUTDIR / "ppg12_fig13_signal_mc_sdcc_vs_current_photonjet_bins.csv"
    out_manifest = OUTDIR / "ppg12_fig13_signal_mc_sdcc_vs_current_photonjet_manifest.json"

    centers = np.asarray(ppg12["centers"], dtype=float)
    edges = np.asarray(ppg12["edges"], dtype=float)
    ref = np.asarray(ppg12["values"], dtype=float)
    ref_err = np.asarray(ppg12["errors"], dtype=float)
    cur = np.asarray(current["values"], dtype=float)
    cur_err = np.asarray(current["errors"], dtype=float)
    ratio = np.divide(ref, cur, out=np.full_like(ref, np.nan), where=cur > 0)
    ratio_err = ratio * np.sqrt(
        np.divide(ref_err, ref, out=np.zeros_like(ref), where=ref > 0) ** 2
        + np.divide(cur_err, cur, out=np.zeros_like(cur), where=cur > 0) ** 2
    )

    with out_csv.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "bin",
                "bdt_low",
                "bdt_high",
                "bdt_center",
                "ppg12_signal_mc",
                "ppg12_signal_mc_error",
                "current_signal_mc",
                "current_signal_mc_error",
                "ppg12_over_current",
                "ppg12_over_current_error",
                "current_raw_count",
            ]
        )
        for i, center in enumerate(centers):
            writer.writerow(
                [
                    i + 1,
                    f"{edges[i]:.12g}",
                    f"{edges[i + 1]:.12g}",
                    f"{center:.12g}",
                    f"{ref[i]:.12g}",
                    f"{ref_err[i]:.12g}",
                    f"{cur[i]:.12g}",
                    f"{cur_err[i]:.12g}",
                    f"{ratio[i]:.12g}",
                    f"{ratio_err[i]:.12g}",
                    f"{np.asarray(current['raw'], dtype=float)[i]:.12g}",
                ]
            )

    setup_style()
    fig, (ax, rax) = plt.subplots(
        2,
        1,
        figsize=(7.72, 9.98),
        dpi=100,
        sharex=True,
        gridspec_kw={"height_ratios": [3.35, 1.0], "hspace": 0.05},
    )
    ax.errorbar(
        centers,
        ref,
        yerr=ref_err,
        fmt="o",
        color="#d62728",
        markerfacecolor="none",
        markeredgewidth=1.15,
        ms=4.7,
        lw=1.0,
        label="PPG12 SDCC Signal MC",
        zorder=4,
    )
    ax.errorbar(
        centers,
        cur,
        yerr=cur_err,
        fmt="o",
        color="#d62728",
        markerfacecolor="#d62728",
        markeredgewidth=0.8,
        ms=3.8,
        lw=0.9,
        label="Current photon+jet Signal MC",
        zorder=3,
    )
    sphinx_label(ax)
    ax.text(0.05, 0.84, r"$p{+}p\ \sqrt{s}=200$ GeV", transform=ax.transAxes, fontsize=15, ha="left")
    ax.text(0.05, 0.77, r"$|\eta^\gamma| < 0.7$", transform=ax.transAxes, fontsize=15, ha="left")
    ax.text(0.05, 0.70, r"$22 < p_T < 28$ GeV, no NPB cut", transform=ax.transAxes, fontsize=13, ha="left")
    ax.text(
        0.58,
        0.62,
        "Signal MC only\nsame PPG12 50-bin BDT grid",
        transform=ax.transAxes,
        fontsize=14,
        ha="left",
        va="top",
        linespacing=1.35,
    )
    ax.set_ylabel("normalized counts", fontsize=17)
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, max(0.12, float(np.nanmax([np.nanmax(ref + ref_err), np.nanmax(cur + cur_err)])) * 1.18))
    ax.legend(
        loc="upper right",
        bbox_to_anchor=(0.985, 0.965),
        frameon=False,
        fontsize=15.0,
        handlelength=1.6,
        borderaxespad=0.25,
    )
    ax.tick_params(labelsize=14, top=True, right=True)
    ax.minorticks_on()

    finite = np.isfinite(ratio)
    rax.errorbar(centers[finite], ratio[finite], yerr=ratio_err[finite], fmt="o", color="black", ms=3.8, lw=0.9)
    rax.axhline(1.0, color="0.35", ls="--", lw=1.0)
    rax.set_ylabel("PPG12 / Current", fontsize=15)
    rax.set_xlabel("bdt", fontsize=17)
    finite_ratio = ratio[np.isfinite(ratio)]
    rax.set_yscale("log")
    ylo = max(0.08, min(0.55, float(np.nanmin(finite_ratio)) * 0.75)) if finite_ratio.size else 0.5
    yhi = min(300.0, max(1.45, float(np.nanmax(finite_ratio)) * 1.20)) if finite_ratio.size else 1.5
    rax.set_ylim(ylo, yhi)
    rax.tick_params(labelsize=14, top=True, right=True)
    rax.minorticks_on()
    fig.subplots_adjust(left=0.115, right=0.985, top=0.985, bottom=0.075, hspace=0.05)
    fig.savefig(out_png, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)

    stable = np.isfinite(ratio) & (cur > 0.002) & (ref > 0.002)
    manifest = {
        "schema": "the76_ppg12_fig13_signal_mc_sdcc_vs_current_overlay_v1",
        "comparison": "PPG12 SDCC Fig.13 Signal MC BDT projection vs current photon+jet RecoilJets TableQA Signal MC BDT object",
        "ppg12_source_json": str(PPG12_JSON),
        "ppg12_source_root": ppg12["source_root"],
        "ppg12_source_hist": ppg12["source_hist"],
        "ppg12_curve_key": "panels.bdt_no_npb_22_28.curves.signal_mc",
        "current_pointer": str(CURRENT_POINTER),
        "current_pointer_payload": pointer,
        "current_root": str(current_root),
        "current_hist": CURRENT_HIST,
        "current_hist_note": "TableQA BDT score for eta0/pt3(22-28 GeV)/cut0(no NPB), rebinned from 100 bins to the PPG12 50-bin BDT grid.",
        "normalization": "Each Signal MC projection is normalized to unit integral over 0 <= bdt <= 1 before ratio.",
        "output_png": str(out_png),
        "output_csv": str(out_csv),
        "stable_ratio_threshold": "ratio summary uses bins with both normalized contents > 0.002",
        "stable_ratio_bin_count": int(np.sum(stable)),
        "stable_ratio_mean_ppg12_over_current": float(np.nanmean(ratio[stable])) if np.any(stable) else None,
        "stable_ratio_max_abs_minus_one": float(np.nanmax(np.abs(ratio[stable] - 1.0))) if np.any(stable) else None,
        "all_bins_max_abs_ratio_minus_one": float(np.nanmax(np.abs(finite_ratio - 1.0))) if finite_ratio.size else None,
        "current_raw_integral_0to1": current["raw_integral_0to1"],
        "current_entries": current["entries"],
        "ppg12_panel_summary": ppg12["panel_summary"],
    }
    out_manifest.write_text(json.dumps(manifest, indent=2, sort_keys=True))
    return out_png, out_csv, out_manifest


def main() -> None:
    current_root, pointer = load_current_root()
    ppg12 = load_ppg12_signal()
    current = load_current_signal(current_root, np.asarray(ppg12["edges"], dtype=float))
    png, csv_path, manifest = write_outputs(ppg12, current, current_root, pointer)
    print(png)
    print(csv_path)
    print(manifest)


if __name__ == "__main__":
    main()
