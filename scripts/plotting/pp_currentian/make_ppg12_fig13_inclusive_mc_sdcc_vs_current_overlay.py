#!/usr/bin/env python3
"""Overlay PPG12 Fig. 13 Inclusive MC BDT template with inclusive-jet output.

The PPG12 side is the SDCC-extracted Inclusive MC curve from the Fig. 13/19
four-curve readback. The current side uses the available inclusive-jet slice
ROOTs, summing all-candidate tight-BDT score histograms over 22-28 GeV and
normalizing the shape on the same PPG12 50-bin BDT grid.
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
INCLUSIVE_ENTRY = (
    REPO
    / "dataOutput/current_recoiljets_artifacts/entries/pp_sim_inclusivejet_merged"
    / "2026-07-02T163931p0000_pp_sim_inclusivejet_merged_the76_ppg12_fig6_inclusive_strict_20260702_014757.json"
)
THE94_ROOT_DIR = (
    REPO
    / "dataOutput/ppg12Parity/the94_ppg12_inclusivejet_fig6_fixed_20260702_204032"
    / "final_roots/inclusivejet"
)
PT_HISTS = (
    "SIM/h_tightBDTScore_allCandidates_pT_22_24",
    "SIM/h_tightBDTScore_allCandidates_pT_24_26",
    "SIM/h_tightBDTScore_allCandidates_pT_26_28",
)
OUTDIR = REPO / "dataOutput/ppg12Parity/the76_ppg12_fig13_inclusive_mc_overlay_20260706"


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


def load_ppg12_inclusive() -> dict[str, object]:
    payload = json.loads(PPG12_JSON.read_text())
    panel = payload["panels"]["bdt_no_npb_22_28"]
    curve = panel["curves"]["inclusive_mc"]
    return {
        "centers": np.asarray(curve["centers"], dtype=float),
        "edges": np.asarray(curve["edges"], dtype=float),
        "values": np.asarray(curve["values"], dtype=float),
        "errors": np.asarray(curve["errors"], dtype=float),
        "source_root": panel["summary"]["source_files"]["inclusive"],
        "source_hist": panel["hist_name"],
        "panel_summary": panel["summary"],
    }


def registered_inclusive_roots() -> tuple[list[Path], dict[str, object]]:
    entry = json.loads(INCLUSIVE_ENTRY.read_text())
    roots = [Path(p) for p in entry.get("root_paths", [])]
    # The registered entry is paused and its older roots do not carry this BDT
    # family. Use the newer local THE-94 fixed roots when present, while keeping
    # the paused registry entry in the manifest.
    fixed = sorted(THE94_ROOT_DIR.glob("RecoilJets_jet*_ALL_*.root"))
    if fixed:
        return fixed, entry
    return roots, entry


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


def load_current_inclusive(dst_edges: np.ndarray) -> dict[str, object]:
    roots, entry = registered_inclusive_roots()
    raw = np.zeros(len(dst_edges) - 1, dtype=float)
    raw_err2 = np.zeros_like(raw)
    used: list[dict[str, object]] = []
    missing: list[dict[str, object]] = []

    for root_path in roots:
        with uproot.open(root_path) as root_file:
            for hist_name in PT_HISTS:
                key = hist_name + ";1"
                if key not in root_file:
                    missing.append({"root": str(root_path), "hist": hist_name})
                    continue
                values, edges = root_file[hist_name].to_numpy(flow=False)
                rebinned, rebinned_err = rebin_counts_to_edges(
                    np.asarray(values, dtype=float),
                    np.asarray(edges, dtype=float),
                    dst_edges,
                )
                raw += rebinned
                raw_err2 += rebinned_err * rebinned_err
                used.append(
                    {
                        "root": str(root_path),
                        "hist": hist_name,
                        "integral_native": float(np.sum(values)),
                        "integral_rebinned_0to1": float(np.sum(rebinned)),
                        "axis": [float(edges[0]), float(edges[-1])],
                    }
                )

    total = float(np.sum(raw))
    errors = np.sqrt(raw_err2)
    norm = raw / total if total > 0 else raw
    norm_err = errors / total if total > 0 else errors
    centers = 0.5 * (dst_edges[:-1] + dst_edges[1:])
    return {
        "centers": centers,
        "edges": dst_edges,
        "values": norm,
        "errors": norm_err,
        "raw": raw,
        "raw_errors": errors,
        "raw_integral_0to1": total,
        "used_inputs": used,
        "missing_inputs": missing,
        "registered_entry": entry,
        "root_paths": [str(p) for p in roots],
    }


def write_outputs(ppg12: dict[str, object], current: dict[str, object]) -> tuple[Path, Path, Path]:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    out_png = OUTDIR / "ppg12_fig13_inclusive_mc_sdcc_vs_current_inclusivejet_overlay_ratio.png"
    out_csv = OUTDIR / "ppg12_fig13_inclusive_mc_sdcc_vs_current_inclusivejet_bins.csv"
    out_manifest = OUTDIR / "ppg12_fig13_inclusive_mc_sdcc_vs_current_inclusivejet_manifest.json"

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
                "ppg12_inclusive_mc",
                "ppg12_inclusive_mc_error",
                "current_inclusive_mc",
                "current_inclusive_mc_error",
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
        color="#1f5fff",
        markerfacecolor="none",
        markeredgewidth=1.15,
        ms=4.7,
        lw=1.0,
        label="PPG12 SDCC Inclusive MC",
        zorder=4,
    )
    ax.errorbar(
        centers,
        cur,
        yerr=cur_err,
        fmt="o",
        color="#1f5fff",
        markerfacecolor="#1f5fff",
        markeredgewidth=0.8,
        ms=3.8,
        lw=0.9,
        label="Current inclusive-jet MC",
        zorder=3,
    )
    sphinx_label(ax)
    ax.text(0.05, 0.84, r"$p{+}p\ \sqrt{s}=200$ GeV", transform=ax.transAxes, fontsize=15, ha="left")
    ax.text(0.05, 0.77, r"$|\eta^\gamma| < 0.7$", transform=ax.transAxes, fontsize=15, ha="left")
    ax.text(0.05, 0.70, r"$22 < p_T < 28$ GeV, no NPB cut", transform=ax.transAxes, fontsize=13, ha="left")
    ax.text(
        0.55,
        0.62,
        "Inclusive MC only\nsame PPG12 50-bin BDT grid",
        transform=ax.transAxes,
        fontsize=14,
        ha="left",
        va="top",
        linespacing=1.35,
    )
    ax.set_ylabel("normalized counts", fontsize=17)
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, max(0.45, float(np.nanmax([np.nanmax(ref + ref_err), np.nanmax(cur + cur_err)])) * 1.16))
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
        "schema": "the76_ppg12_fig13_inclusive_mc_sdcc_vs_current_overlay_v1",
        "comparison": "PPG12 SDCC Fig.13 Inclusive MC BDT projection vs available inclusive-jet RecoilJets BDT score objects",
        "ppg12_source_json": str(PPG12_JSON),
        "ppg12_source_root": ppg12["source_root"],
        "ppg12_source_hist": ppg12["source_hist"],
        "ppg12_curve_key": "panels.bdt_no_npb_22_28.curves.inclusive_mc",
        "current_registry_entry": str(INCLUSIVE_ENTRY),
        "current_registry_status": current["registered_entry"].get("status"),
        "current_registry_canonical_status": current["registered_entry"].get("canonical_status"),
        "current_registry_plot_policy": current["registered_entry"].get("plot_policy"),
        "current_roots": current["root_paths"],
        "current_histograms_summed": list(PT_HISTS),
        "current_used_inputs": current["used_inputs"],
        "current_missing_inputs": current["missing_inputs"],
        "current_hist_note": "Summed all-candidate tight-BDT score histograms over 22-28 GeV across local THE-94 inclusive-jet slice roots, then rebinned to the PPG12 50-bin BDT grid.",
        "normalization": "Each Inclusive MC projection is normalized to unit integral over 0 <= bdt <= 1 before ratio.",
        "important_caveat": "The inclusive-jet registry entry is paused/not canonical; this is an available-artifact diagnostic, not a promoted current inclusive-jet baseline.",
        "output_png": str(out_png),
        "output_csv": str(out_csv),
        "stable_ratio_threshold": "ratio summary uses bins with both normalized contents > 0.002",
        "stable_ratio_bin_count": int(np.sum(stable)),
        "stable_ratio_mean_ppg12_over_current": float(np.nanmean(ratio[stable])) if np.any(stable) else None,
        "stable_ratio_max_abs_minus_one": float(np.nanmax(np.abs(ratio[stable] - 1.0))) if np.any(stable) else None,
        "all_bins_max_abs_ratio_minus_one": float(np.nanmax(np.abs(finite_ratio - 1.0))) if finite_ratio.size else None,
        "current_raw_integral_0to1": current["raw_integral_0to1"],
        "ppg12_panel_summary": ppg12["panel_summary"],
    }
    out_manifest.write_text(json.dumps(manifest, indent=2, sort_keys=True))
    return out_png, out_csv, out_manifest


def main() -> None:
    ppg12 = load_ppg12_inclusive()
    current = load_current_inclusive(np.asarray(ppg12["edges"], dtype=float))
    png, csv_path, manifest = write_outputs(ppg12, current)
    print(png)
    print(csv_path)
    print(manifest)


if __name__ == "__main__":
    main()
