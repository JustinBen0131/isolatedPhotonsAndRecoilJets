#!/usr/bin/env python3
"""Overlay PPG12 Fig. 19 after-NPB MC BDT curves with current outputs.

Produces separate Signal MC and Inclusive MC diagnostic overlays for the
18-22 GeV with-NPB/preselection BDT panel.  The PPG12 side is read from the
existing SDCC projection JSON.  The current signal side uses the exact
TableQA 1D object, while the current inclusive side sums the available
inclusive-jet preselected BDT-score slice histograms.
"""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass
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
PHOTONJET_CURRENT_POINTER = REPO / "dataOutput/current_recoiljets_artifacts/current/pp_sim_photonjet_merged/current.json"
PHOTONJET_FALLBACK_ROOT = (
    REPO
    / "dataOutput/ppg12Parity/the76_ppg12_fig24_photonjet_fix_20260702_014217"
    / "final_roots/photonjet/RecoilJets_photonjet5plus10plus20_MERGED.root"
)
INCLUSIVE_ENTRY = (
    REPO
    / "dataOutput/current_recoiljets_artifacts/entries/pp_sim_inclusivejet_merged"
    / "2026-07-02T163931p0000_pp_sim_inclusivejet_merged_the76_ppg12_fig6_inclusive_strict_20260702_014757.json"
)
THE94_INCLUSIVE_ROOT_DIR = (
    REPO
    / "dataOutput/ppg12Parity/the94_ppg12_inclusivejet_fig6_fixed_20260702_204032"
    / "final_roots/inclusivejet"
)
OUTDIR = REPO / "dataOutput/ppg12Parity/the76_ppg12_fig19_after_npb_mc_overlays_20260706"

PANEL_KEY = "bdt_with_npb_18_22"
PPG12_HIST = "h2d_bdt_eta0_pt2_cut1"
CURRENT_SIGNAL_HIST = "SIM/h1d_bdt_eta0_pt2_cut1"
CURRENT_INCLUSIVE_HISTS = (
    "SIM/h_tightBDTScore_preselected_pT_18_20",
    "SIM/h_tightBDTScore_preselected_pT_20_22",
)


@dataclass(frozen=True)
class PlotSpec:
    key: str
    label: str
    curve_key: str
    color: str
    current_kind: str
    output_stem: str


SPECS = (
    PlotSpec(
        key="signal",
        label="Signal MC",
        curve_key="signal_mc",
        color="#e41a1c",
        current_kind="photonjet_signal",
        output_stem="ppg12_fig19_after_npb_signal_mc_sdcc_vs_current_photonjet",
    ),
    PlotSpec(
        key="inclusive",
        label="Inclusive MC",
        curve_key="inclusive_mc",
        color="#1f5fff",
        current_kind="inclusivejet",
        output_stem="ppg12_fig19_after_npb_inclusive_mc_sdcc_vs_current_inclusivejet",
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


def sphinx_label(ax, x: float = 0.055, y: float = 0.925, fs: float = 17) -> None:
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


def load_ppg12_curve(spec: PlotSpec) -> dict[str, object]:
    payload = json.loads(PPG12_JSON.read_text())
    panel = payload["panels"][PANEL_KEY]
    curve = panel["curves"][spec.curve_key]
    source_key = "signal" if spec.key == "signal" else "inclusive"
    return {
        "centers": np.asarray(curve["centers"], dtype=float),
        "edges": np.asarray(curve["edges"], dtype=float),
        "values": np.asarray(curve["values"], dtype=float),
        "errors": np.asarray(curve["errors"], dtype=float),
        "source_root": panel["summary"]["source_files"][source_key],
        "source_hist": panel["hist_name"],
        "panel_summary": panel["summary"],
    }


def photonjet_current_root() -> tuple[Path, dict[str, object]]:
    if PHOTONJET_CURRENT_POINTER.exists():
        pointer = json.loads(PHOTONJET_CURRENT_POINTER.read_text())
        paths = pointer.get("root_paths") or []
        if paths:
            return Path(paths[0]), pointer
    return PHOTONJET_FALLBACK_ROOT, {
        "note": "Fallback path used because current pointer was unavailable.",
        "root_paths": [str(PHOTONJET_FALLBACK_ROOT)],
    }


def inclusive_roots() -> tuple[list[Path], dict[str, object]]:
    entry = json.loads(INCLUSIVE_ENTRY.read_text())
    fixed = sorted(THE94_INCLUSIVE_ROOT_DIR.glob("RecoilJets_jet*_ALL_*.root"))
    if fixed:
        return fixed, entry
    return [Path(p) for p in entry.get("root_paths", [])], entry


def rebin_counts_to_edges(
    src_values: np.ndarray,
    src_edges: np.ndarray,
    dst_edges: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
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


def load_current_signal(dst_edges: np.ndarray) -> dict[str, object]:
    root_path, pointer = photonjet_current_root()
    with uproot.open(root_path) as root_file:
        hist = root_file[CURRENT_SIGNAL_HIST]
        values, edges = hist.to_numpy(flow=False)
        entries = float(hist.member("fEntries"))
    raw, raw_err = rebin_counts_to_edges(np.asarray(values, dtype=float), np.asarray(edges, dtype=float), dst_edges)
    return {
        "raw": raw,
        "raw_errors": raw_err,
        "raw_integral_0to1": float(np.sum(raw)),
        "entries": entries,
        "root_paths": [str(root_path)],
        "histograms": [CURRENT_SIGNAL_HIST],
        "used_inputs": [
            {
                "root": str(root_path),
                "hist": CURRENT_SIGNAL_HIST,
                "integral_native": float(np.sum(values)),
                "integral_rebinned_0to1": float(np.sum(raw)),
                "entries": entries,
            }
        ],
        "missing_inputs": [],
        "registry_entry": pointer,
        "source_note": "Current photon+jet Signal MC uses exact TableQA h1d_bdt_eta0_pt2_cut1.",
    }


def load_current_inclusive(dst_edges: np.ndarray) -> dict[str, object]:
    roots, entry = inclusive_roots()
    raw = np.zeros(len(dst_edges) - 1, dtype=float)
    err2 = np.zeros_like(raw)
    used: list[dict[str, object]] = []
    missing: list[dict[str, object]] = []
    for root_path in roots:
        with uproot.open(root_path) as root_file:
            for hist_name in CURRENT_INCLUSIVE_HISTS:
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
                err2 += rebinned_err * rebinned_err
                used.append(
                    {
                        "root": str(root_path),
                        "hist": hist_name,
                        "integral_native": float(np.sum(values)),
                        "integral_rebinned_0to1": float(np.sum(rebinned)),
                        "axis": [float(edges[0]), float(edges[-1])],
                    }
                )
    return {
        "raw": raw,
        "raw_errors": np.sqrt(err2),
        "raw_integral_0to1": float(np.sum(raw)),
        "root_paths": [str(p) for p in roots],
        "histograms": list(CURRENT_INCLUSIVE_HISTS),
        "used_inputs": used,
        "missing_inputs": missing,
        "registry_entry": entry,
        "source_note": "Current Inclusive MC sums preselected BDT-score histograms over 18-22 GeV across available THE-94 inclusive-jet slice roots.",
    }


def normalize_current(current: dict[str, object], dst_edges: np.ndarray) -> dict[str, object]:
    raw = np.asarray(current["raw"], dtype=float)
    raw_err = np.asarray(current["raw_errors"], dtype=float)
    total = float(np.sum(raw))
    values = raw / total if total > 0 else raw
    errors = raw_err / total if total > 0 else raw_err
    return {
        **current,
        "centers": 0.5 * (dst_edges[:-1] + dst_edges[1:]),
        "edges": dst_edges,
        "values": values,
        "errors": errors,
    }


def draw_and_write(spec: PlotSpec) -> tuple[Path, Path, Path]:
    ppg12 = load_ppg12_curve(spec)
    dst_edges = np.asarray(ppg12["edges"], dtype=float)
    if spec.current_kind == "photonjet_signal":
        current = normalize_current(load_current_signal(dst_edges), dst_edges)
    elif spec.current_kind == "inclusivejet":
        current = normalize_current(load_current_inclusive(dst_edges), dst_edges)
    else:
        raise ValueError(spec.current_kind)

    OUTDIR.mkdir(parents=True, exist_ok=True)
    out_png = OUTDIR / f"{spec.output_stem}_overlay_ratio.png"
    out_csv = OUTDIR / f"{spec.output_stem}_bins.csv"
    out_manifest = OUTDIR / f"{spec.output_stem}_manifest.json"

    centers = np.asarray(ppg12["centers"], dtype=float)
    edges = np.asarray(ppg12["edges"], dtype=float)
    ref = np.asarray(ppg12["values"], dtype=float)
    ref_err = np.asarray(ppg12["errors"], dtype=float)
    cur = np.asarray(current["values"], dtype=float)
    cur_err = np.asarray(current["errors"], dtype=float)
    ratio = np.divide(cur, ref, out=np.full_like(cur, np.nan), where=ref > 0)
    ratio_err = ratio * np.sqrt(
        np.divide(cur_err, cur, out=np.zeros_like(cur), where=cur > 0) ** 2
        + np.divide(ref_err, ref, out=np.zeros_like(ref), where=ref > 0) ** 2
    )

    with out_csv.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "bin",
                "bdt_low",
                "bdt_high",
                "bdt_center",
                f"ppg12_{spec.key}_mc",
                f"ppg12_{spec.key}_mc_error",
                f"current_{spec.key}_mc",
                f"current_{spec.key}_mc_error",
                "current_over_ppg12",
                "current_over_ppg12_error",
                "current_raw_count",
            ]
        )
        raw = np.asarray(current["raw"], dtype=float)
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
                    f"{raw[i]:.12g}",
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
        color=spec.color,
        markerfacecolor="none",
        markeredgewidth=1.15,
        ms=4.7,
        lw=1.0,
        label=f"PPG12 SDCC {spec.label}",
        zorder=4,
    )
    ax.errorbar(
        centers,
        cur,
        yerr=cur_err,
        fmt="o",
        color=spec.color,
        markerfacecolor=spec.color,
        markeredgewidth=0.8,
        ms=3.8,
        lw=0.9,
        label=f"Current {spec.label}",
        zorder=3,
    )
    sphinx_label(ax)
    ax.text(0.055, 0.84, r"$p{+}p\ \sqrt{s}=200$ GeV", transform=ax.transAxes, fontsize=15, ha="left")
    ax.text(0.055, 0.77, r"$|\eta^\gamma| < 0.7$", transform=ax.transAxes, fontsize=15, ha="left")
    ax.text(0.055, 0.70, r"$18 < p_T < 22$ GeV, with NPB cut", transform=ax.transAxes, fontsize=13, ha="left")
    ax.text(
        0.59,
        0.62,
        f"{spec.label} only\nsame PPG12 50-bin BDT grid",
        transform=ax.transAxes,
        fontsize=14,
        ha="left",
        va="top",
        linespacing=1.35,
    )
    ax.set_ylabel("normalized counts", fontsize=17)
    ax.set_xlim(0.0, 1.0)
    ymax = max(0.24, float(np.nanmax([np.nanmax(ref + ref_err), np.nanmax(cur + cur_err)])) * 1.18)
    ax.set_ylim(0.0, ymax)
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
    rax.errorbar(
        centers[finite],
        ratio[finite],
        yerr=ratio_err[finite],
        fmt="o",
        color="black",
        ms=3.8,
        lw=0.9,
    )
    rax.axhline(1.0, color="0.35", ls="--", lw=1.0)
    rax.set_ylabel("Current / PPG12", fontsize=15)
    rax.set_xlabel("bdt", fontsize=17)
    finite_ratio = ratio[np.isfinite(ratio)]
    ylo = max(0.0, min(0.45, float(np.nanmin(finite_ratio)) * 0.78)) if finite_ratio.size else 0.5
    yhi = min(5.0, max(1.55, float(np.nanmax(finite_ratio)) * 1.18)) if finite_ratio.size else 1.5
    rax.set_ylim(ylo, yhi)
    rax.tick_params(labelsize=14, top=True, right=True)
    rax.minorticks_on()
    fig.subplots_adjust(left=0.115, right=0.985, top=0.985, bottom=0.075, hspace=0.05)
    fig.savefig(out_png, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)

    stable = np.isfinite(ratio) & (cur > 0.002) & (ref > 0.002)
    manifest = {
        "schema": "the76_ppg12_fig19_after_npb_mc_sdcc_vs_current_overlay_v1",
        "comparison": f"PPG12 SDCC Fig.19 after-NPB {spec.label} BDT projection vs current output",
        "panel_key": PANEL_KEY,
        "ppg12_source_json": str(PPG12_JSON),
        "ppg12_source_root": ppg12["source_root"],
        "ppg12_source_hist": PPG12_HIST,
        "ppg12_curve_key": f"panels.{PANEL_KEY}.curves.{spec.curve_key}",
        "current_kind": spec.current_kind,
        "current_roots": current["root_paths"],
        "current_histograms": current["histograms"],
        "current_used_inputs": current["used_inputs"],
        "current_missing_inputs": current["missing_inputs"],
        "current_registry_status": current["registry_entry"].get("status"),
        "current_registry_canonical_status": current["registry_entry"].get("canonical_status"),
        "current_registry_plot_policy": current["registry_entry"].get("plot_policy"),
        "current_source_note": current["source_note"],
        "normalization": "Each MC projection is normalized to unit integral over 0 <= bdt <= 1 before ratio.",
        "ratio_definition": "Current output / PPG12 SDCC",
        "important_caveat": (
            "Signal uses the registered photon+jet current ROOT. Inclusive uses available THE-94 inclusive-jet slice ROOTs; "
            "the inclusive-jet registry entry remains paused/not canonical, so the inclusive plot is diagnostic only."
        ),
        "output_png": str(out_png),
        "output_csv": str(out_csv),
        "stable_ratio_threshold": "ratio summary uses bins with both normalized contents > 0.002",
        "stable_ratio_bin_count": int(np.sum(stable)),
        "stable_ratio_mean_current_over_ppg12": float(np.nanmean(ratio[stable])) if np.any(stable) else None,
        "stable_ratio_max_abs_minus_one": float(np.nanmax(np.abs(ratio[stable] - 1.0))) if np.any(stable) else None,
        "all_bins_max_abs_ratio_minus_one": float(np.nanmax(np.abs(finite_ratio - 1.0))) if finite_ratio.size else None,
        "current_raw_integral_0to1": current["raw_integral_0to1"],
        "ppg12_panel_summary": ppg12["panel_summary"],
    }
    out_manifest.write_text(json.dumps(manifest, indent=2, sort_keys=True))
    return out_png, out_csv, out_manifest


def main() -> None:
    for spec in SPECS:
        for path in draw_and_write(spec):
            print(path)


if __name__ == "__main__":
    main()
