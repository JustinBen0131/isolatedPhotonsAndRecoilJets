#!/usr/bin/env python3
"""Make THE-93 BDT before-NPB overlays against PPG12 Fig. 13/19 SDCC curves.

This uses the exact no-NPB 22-28 GeV PPG12 panel object family:
``h2d_bdt_eta0_pt3_cut0``.  Both PPG12 and current RecoilJets are projected
onto the BDT-score axis, rebinned to the same 50-bin 0..1 grid, then
unit-normalized.  The pp data current object is read from the
PPG12-scaled-trigger namespace, not the TriggerAnalyzer diagnostic namespace.
"""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import uproot


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
CAMPAIGN = "the93_ppg12_canonical_full_20260706_2145"
PPG12_JSON = (
    REPO
    / "dataOutput/ppg12PhotonYield/ppg12_photon_yield_v1_data_20260620"
    / "shower_shape_reference_validation/fig13_19_bdt"
    / "ppg12_sdcc_fig19_bdt_fourcurve_projections.json"
)
CURRENT_POINTER_DIR = REPO / "dataOutput/current_recoiljets_artifacts/current"
DEFAULT_OUTDIR = REPO / "dataOutput/ppg12Parity" / CAMPAIGN / "bdt_fig19_before_npb"


@dataclass(frozen=True)
class OverlaySpec:
    key: str
    title: str
    pointer_key: str
    current_hist: str
    ppg12_curve: str
    color: str
    ppg12_label: str
    current_label: str
    ratio_label: str
    ratio_direction: str
    y_minimum: float
    note: str


SPECS = (
    OverlaySpec(
        key="data",
        title="Data",
        pointer_key="pp_data_merged",
        current_hist="PPG12_scaledtrigger30/h2d_bdt_eta0_pt3_cut0",
        ppg12_curve="data",
        color="#1f77b4",
        ppg12_label="PPG12 SDCC ROOT data",
        current_label="THE-93 pp data",
        ratio_label="Current / PPG12",
        ratio_direction="current_over_ppg12",
        y_minimum=0.84,
        note="data uses direct GL1 ScaledVector bit 30 namespace PPG12_scaledtrigger30",
    ),
    OverlaySpec(
        key="photonjet_signal",
        title="Photon+jet Signal MC",
        pointer_key="pp_sim_photonjet_merged",
        current_hist="SIM/h2d_bdt_eta0_pt3_cut0",
        ppg12_curve="signal_mc",
        color="#e41a1c",
        ppg12_label="PPG12 SDCC Signal MC",
        current_label="THE-93 photon+jet Signal MC",
        ratio_label="PPG12 / Current",
        ratio_direction="ppg12_over_current",
        y_minimum=0.145,
        note="signal MC uses exact TableQA pt3 cut0 BDT object",
    ),
    OverlaySpec(
        key="inclusivejet",
        title="Inclusive-jet MC",
        pointer_key="pp_sim_inclusivejet_merged",
        current_hist="SIM/h2d_bdt_eta0_pt3_cut0",
        ppg12_curve="inclusive_mc",
        color="#2a6fdb",
        ppg12_label="PPG12 SDCC Inclusive MC",
        current_label="THE-93 inclusive-jet MC",
        ratio_label="PPG12 / Current",
        ratio_direction="ppg12_over_current",
        y_minimum=0.62,
        note="inclusive MC uses exact TableQA pt3 cut0 BDT object",
    ),
)


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


def load_pointer(pointer_key: str) -> dict[str, Any]:
    path = CURRENT_POINTER_DIR / pointer_key / "current.json"
    payload = json.loads(path.read_text())
    roots = payload.get("root_paths") or []
    if not roots:
        raise RuntimeError(f"{path} has no root_paths")
    payload["_pointer_path"] = str(path)
    payload["_root_path"] = roots[0]
    return payload


def load_ppg12_panel() -> dict[str, Any]:
    payload = json.loads(PPG12_JSON.read_text())
    return payload["panels"]["bdt_no_npb_22_28"]


def rebin_counts_to_edges(
    src_values: np.ndarray,
    src_edges: np.ndarray,
    dst_edges: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    return rebin_counts_and_variances_to_edges(
        src_values,
        np.maximum(src_values, 0.0),
        src_edges,
        dst_edges,
    )


def rebin_counts_and_variances_to_edges(
    src_values: np.ndarray,
    src_variances: np.ndarray,
    src_edges: np.ndarray,
    dst_edges: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    out = np.zeros(len(dst_edges) - 1, dtype=float)
    err2 = np.zeros_like(out)
    for count, variance, lo, hi in zip(src_values, src_variances, src_edges[:-1], src_edges[1:]):
        width = float(hi - lo)
        if width <= 0:
            continue
        for j, (dst_lo, dst_hi) in enumerate(zip(dst_edges[:-1], dst_edges[1:])):
            overlap = max(0.0, min(float(hi), float(dst_hi)) - max(float(lo), float(dst_lo)))
            if overlap <= 0:
                continue
            frac = overlap / width
            out[j] += float(count) * frac
            err2[j] += max(float(variance), 0.0) * frac * frac
    return out, np.sqrt(err2)


def normalize(raw: np.ndarray, raw_err: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    total = float(np.sum(raw))
    if total <= 0:
        return raw.copy(), raw_err.copy(), total
    return raw / total, raw_err / total, total


def load_current_curve(spec: OverlaySpec, dst_edges: np.ndarray) -> dict[str, Any]:
    pointer = load_pointer(spec.pointer_key)
    root_path = Path(pointer["_root_path"])
    with uproot.open(root_path) as root_file:
        if spec.current_hist not in root_file:
            raise KeyError(f"{root_path} is missing {spec.current_hist}")
        hist = root_file[spec.current_hist]
        hist_arrays = hist.to_numpy(flow=False)
        values = np.asarray(hist_arrays[0], dtype=float)
        variances = hist.variances(flow=False)
        if variances is None:
            variances = np.maximum(values, 0.0)
        variances = np.asarray(variances, dtype=float)
        if values.ndim == 2:
            edges = np.asarray(hist_arrays[1], dtype=float)
            yedges = np.asarray(hist_arrays[2], dtype=float)
            projected_values = values.sum(axis=1)
            projected_variances = variances.sum(axis=1)
            projection_note = (
                "TH2 projected onto x/BDT axis by summing all non-flow "
                "E_T^iso bins, matching PPG12 ProjectionX/RebinX macro mode"
            )
            source_axis = [float(edges[0]), float(edges[-1])]
            source_y_axis = [float(yedges[0]), float(yedges[-1])]
        elif values.ndim == 1:
            edges = np.asarray(hist_arrays[1], dtype=float)
            projected_values = values
            projected_variances = variances
            projection_note = "TH1 direct BDT histogram"
            source_axis = [float(edges[0]), float(edges[-1])]
            source_y_axis = None
        else:
            raise RuntimeError(f"Unsupported histogram dimension for {spec.current_hist}: {values.ndim}")
        entries = float(hist.member("fEntries"))
    raw, raw_err = rebin_counts_and_variances_to_edges(projected_values, projected_variances, edges, dst_edges)
    vals, errs, total = normalize(raw, raw_err)
    return {
        "values": vals,
        "errors": errs,
        "raw": raw,
        "raw_errors": raw_err,
        "raw_integral_0to1": total,
        "entries": entries,
        "root": str(root_path),
        "pointer": pointer,
        "hist": spec.current_hist,
        "projection_note": projection_note,
        "source_axis": source_axis,
        "source_y_axis": source_y_axis,
        "source_integral": float(np.sum(projected_values)),
    }


def ratio_values(spec: OverlaySpec, ref: np.ndarray, ref_err: np.ndarray, cur: np.ndarray, cur_err: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    if spec.ratio_direction == "current_over_ppg12":
        numerator, numerator_err = cur, cur_err
        denominator, denominator_err = ref, ref_err
    else:
        numerator, numerator_err = ref, ref_err
        denominator, denominator_err = cur, cur_err
    ratio = np.divide(numerator, denominator, out=np.full_like(numerator, np.nan), where=denominator > 0)
    rel2 = (
        np.divide(numerator_err, numerator, out=np.zeros_like(numerator), where=numerator > 0) ** 2
        + np.divide(denominator_err, denominator, out=np.zeros_like(denominator), where=denominator > 0) ** 2
    )
    return ratio, ratio * np.sqrt(rel2)


def infer_effective_entries(values: np.ndarray, errors: np.ndarray) -> float:
    terms = [(float(v) / float(e)) ** 2 for v, e in zip(values, errors) if v > 0 and e > 0]
    return float(sum(terms))


def sphinx_label(ax: plt.Axes, x: float = 0.055, y: float = 0.93, fs: float = 15.0) -> None:
    ax.text(x, y, "sPHENIX", transform=ax.transAxes, ha="left", va="top", fontsize=fs, fontweight="bold", fontstyle="italic")
    ax.text(x + 0.265, y, "Internal", transform=ax.transAxes, ha="left", va="top", fontsize=fs)


def draw_overlay(
    ax: plt.Axes,
    rax: plt.Axes,
    spec: OverlaySpec,
    centers: np.ndarray,
    ppg12: dict[str, Any],
    current: dict[str, Any],
    *,
    add_stats: bool,
) -> dict[str, Any]:
    ref = np.asarray(ppg12["values"], dtype=float)
    ref_err = np.asarray(ppg12["errors"], dtype=float)
    cur = np.asarray(current["values"], dtype=float)
    cur_err = np.asarray(current["errors"], dtype=float)
    ratio, ratio_err = ratio_values(spec, ref, ref_err, cur, cur_err)

    ax.errorbar(
        centers,
        ref,
        yerr=ref_err,
        fmt="o",
        color=spec.color,
        markerfacecolor="none",
        markeredgewidth=1.1,
        ms=4.2,
        lw=0.9,
        label=spec.ppg12_label,
        zorder=4,
    )
    ax.errorbar(
        centers,
        cur,
        yerr=cur_err,
        fmt="o",
        color=spec.color,
        markerfacecolor=spec.color,
        markeredgewidth=0.7,
        ms=3.2,
        lw=0.8,
        label=spec.current_label,
        zorder=3,
    )
    sphinx_label(ax)
    ax.text(0.055, 0.835, r"$p{+}p\ \sqrt{s}=200$ GeV", transform=ax.transAxes, fontsize=12.5, ha="left")
    ax.text(0.055, 0.755, r"$|\eta^\gamma| < 0.7$", transform=ax.transAxes, fontsize=12.5, ha="left")
    ax.text(0.055, 0.675, r"$22 < p_T < 28$ GeV, no NPB cut", transform=ax.transAxes, fontsize=10.8, ha="left")
    ax.set_xlim(0.0, 1.0)
    ymax = max(spec.y_minimum, float(np.nanmax([np.nanmax(ref + ref_err), np.nanmax(cur + cur_err)])) * 1.16)
    ax.set_ylim(0.0, ymax)
    ax.set_ylabel("normalized counts", fontsize=12.5)
    ax.legend(loc="upper right", frameon=False, fontsize=10.0, handlelength=1.4)
    ax.tick_params(labelsize=10.5, top=True, right=True)
    ax.minorticks_on()

    if spec.key != "data":
        ax.text(
            0.56,
            0.58,
            f"{spec.title} only\nsame PPG12 50-bin BDT grid",
            transform=ax.transAxes,
            fontsize=10.8,
            ha="left",
            va="top",
            linespacing=1.25,
        )
    elif add_stats:
        finite = np.isfinite(ratio) & (ref > 0.0085)
        stable_ratio = ratio[finite]
        max_abs = float(np.nanmax(np.abs(stable_ratio - 1.0))) if stable_ratio.size else float("nan")
        max_x = float(centers[finite][int(np.nanargmax(np.abs(stable_ratio - 1.0)))]) if stable_ratio.size else float("nan")
        stats = (
            f"PPG12 SDCC N_eff={infer_effective_entries(ref, ref_err):.0f}\n"
            f"THE-93 pp N={current['raw_integral_0to1']:.0f}\n"
            rf"max $|R-1|$ = {100.0 * max_abs:.1f}%"
            f"\nnear bdt = {max_x:.2f}"
        )
        ax.text(0.58, 0.76, stats, transform=ax.transAxes, fontsize=10.8, ha="left", va="top", linespacing=1.35)

    finite = np.isfinite(ratio)
    rax.errorbar(centers[finite], ratio[finite], yerr=ratio_err[finite], fmt="o", color="black", ms=2.6, lw=0.75)
    rax.axhline(1.0, color="0.35", ls="--", lw=0.8)
    rax.set_ylabel(spec.ratio_label, fontsize=10.8)
    rax.set_xlabel("bdt", fontsize=12.5)
    if spec.key == "data":
        rax.set_ylim(0.0, 1.8)
    else:
        positive = ratio[np.isfinite(ratio) & (ratio > 0)]
        rax.set_yscale("log")
        if positive.size:
            rax.set_ylim(max(0.05, float(np.nanmin(positive)) * 0.7), min(500.0, max(1.6, float(np.nanmax(positive)) * 1.25)))
    rax.tick_params(labelsize=10.5, top=True, right=True)
    rax.minorticks_on()

    stable = np.isfinite(ratio) & (ref > 0.002) & (cur > 0.002)
    stable_ratio = ratio[stable]
    return {
        "ratio_direction": spec.ratio_direction,
        "stable_bin_count": int(np.sum(stable)),
        "stable_ratio_mean": float(np.nanmean(stable_ratio)) if stable_ratio.size else None,
        "stable_ratio_max_abs_minus_one": float(np.nanmax(np.abs(stable_ratio - 1.0))) if stable_ratio.size else None,
        "raw_integral_0to1": current["raw_integral_0to1"],
        "ppg12_effective_entries": infer_effective_entries(ref, ref_err),
    }


def write_csv(path: Path, centers: np.ndarray, edges: np.ndarray, spec: OverlaySpec, ppg12: dict[str, Any], current: dict[str, Any]) -> None:
    ref = np.asarray(ppg12["values"], dtype=float)
    ref_err = np.asarray(ppg12["errors"], dtype=float)
    cur = np.asarray(current["values"], dtype=float)
    cur_err = np.asarray(current["errors"], dtype=float)
    ratio, ratio_err = ratio_values(spec, ref, ref_err, cur, cur_err)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "bin",
                "bdt_low",
                "bdt_high",
                "bdt_center",
                "ppg12",
                "ppg12_error",
                "current",
                "current_error",
                spec.ratio_direction,
                "ratio_error",
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


def make_plots(outdir: Path) -> dict[str, Any]:
    setup_style()
    outdir.mkdir(parents=True, exist_ok=True)
    panel = load_ppg12_panel()
    edges = np.asarray(panel["curves"]["data"]["edges"], dtype=float)
    centers = 0.5 * (edges[:-1] + edges[1:])
    manifest: dict[str, Any] = {
        "schema": "the93_fig19_bdt_before_npb_three_overlays_v1",
        "campaign": CAMPAIGN,
        "ppg12_json": str(PPG12_JSON),
        "ppg12_hist": panel["hist_name"],
        "ppg12_source_files": panel["summary"]["source_files"],
        "selection": "22 < pT < 28 GeV, |eta| < 0.7, no NPB cut, 50-bin BDT grid",
        "outputs": {},
    }

    combined = plt.figure(figsize=(15.9, 7.5), dpi=160)
    gs = combined.add_gridspec(2, 3, height_ratios=[3.35, 1.0], hspace=0.06, wspace=0.12)

    for col, spec in enumerate(SPECS):
        ppg12_curve = panel["curves"][spec.ppg12_curve]
        ppg12 = {
            "values": np.asarray(ppg12_curve["values"], dtype=float),
            "errors": np.asarray(ppg12_curve["errors"], dtype=float),
            "raw": np.asarray(ppg12_curve.get("raw", []), dtype=float),
            "raw_errors": np.asarray(ppg12_curve.get("raw_errors", []), dtype=float),
        }
        current = load_current_curve(spec, edges)

        fig, (ax, rax) = plt.subplots(
            2,
            1,
            figsize=(7.72, 9.98),
            dpi=100,
            sharex=True,
            gridspec_kw={"height_ratios": [3.35, 1.0], "hspace": 0.05},
        )
        summary = draw_overlay(ax, rax, spec, centers, ppg12, current, add_stats=True)
        png = outdir / f"fig19_bdt_before_npb_{spec.key}_sdcc_vs_current_overlay_ratio.png"
        csv_path = outdir / f"fig19_bdt_before_npb_{spec.key}_sdcc_vs_current_bins.csv"
        fig.subplots_adjust(left=0.12, right=0.985, top=0.985, bottom=0.075, hspace=0.05)
        fig.savefig(png, bbox_inches="tight", pad_inches=0.02)
        plt.close(fig)
        write_csv(csv_path, centers, edges, spec, ppg12, current)

        axc = combined.add_subplot(gs[0, col])
        raxc = combined.add_subplot(gs[1, col], sharex=axc)
        combined_summary = draw_overlay(axc, raxc, spec, centers, ppg12, current, add_stats=(spec.key == "data"))
        axc.set_title(spec.title, fontsize=14, pad=8)
        if col > 0:
            axc.set_ylabel("")
            raxc.set_ylabel("")
        plt.setp(axc.get_xticklabels(), visible=False)

        manifest["outputs"][spec.key] = {
            "png": str(png),
            "csv": str(csv_path),
            "ppg12_curve": spec.ppg12_curve,
            "ppg12_source_root": panel["summary"]["source_files"]["data" if spec.key == "data" else ("signal" if spec.key == "photonjet_signal" else "inclusive")],
            "current_pointer": current["pointer"]["_pointer_path"],
            "current_root": current["root"],
            "current_hist": spec.current_hist,
            "current_hist_note": spec.note,
            "current_source_axis": current["source_axis"],
            "current_source_integral": current["source_integral"],
            "current_entries": current["entries"],
            **summary,
        }
        # Preserve combined draw summary for debugging, without overwriting the
        # individual output numbers if matplotlib autoscaling ever changes.
        manifest["outputs"][spec.key]["combined_draw_summary"] = combined_summary

    combined.suptitle("BDT consistency check with PPG12 - before NPB", fontsize=21, fontweight="bold", y=0.965)
    combined_png = outdir / "fig19_bdt_before_npb_threepanel_sdcc_vs_current_overlay_ratio.png"
    combined.savefig(combined_png, bbox_inches="tight", pad_inches=0.03)
    plt.close(combined)
    manifest["combined_png"] = str(combined_png)
    manifest_path = outdir / "fig19_bdt_before_npb_threepanel_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True))
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    manifest = make_plots(args.outdir)
    print(manifest["combined_png"])
    for item in manifest["outputs"].values():
        print(item["png"])
    print(args.outdir / "fig19_bdt_before_npb_threepanel_manifest.json")


if __name__ == "__main__":
    main()
