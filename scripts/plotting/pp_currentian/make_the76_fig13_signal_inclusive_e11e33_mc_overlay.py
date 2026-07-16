#!/usr/bin/env python3
"""Overlay PPG12 Fig.13 Signal/Inclusive MC E11/E33 templates with current output.

The PPG12 side comes from the SDCC ROOT-extracted Fig.13 JSON.  The current
signal side uses the promoted photon+jet merged ROOT pointer, and the current
inclusive side uses the latest local July-2 inclusive-jet final roots because
inclusive is not yet promoted as a `current` pointer.
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import uproot


REPO = Path(__file__).resolve().parents[3]
PPG12_JSON = REPO / (
    "dataOutput/ppg12PhotonYield/ppg12_photon_yield_v1_data_20260620/"
    "shower_shape_reference_validation/fig13_e11_e33/"
    "ppg12_sdcc_fig13_e11_to_e33_histograms.json"
)
PHOTONJET_POINTER = REPO / "dataOutput/current_recoiljets_artifacts/current/pp_sim_photonjet_merged/current.json"
INCLUSIVE_ROOT_GLOB = REPO / (
    "dataOutput/ppg12Parity/the76_ppg12_fig6_inclusive_strict_20260702_014757/"
    "final_roots/inclusivejet/*.root"
)
OUT_DIR = REPO / (
    "dataOutput/ppg12Parity/the76_ppg12_fig24_photonjet_fix_20260702_014217/"
    "fig13_e11e33_signal_inclusive_mc_overlay"
)
PT_HISTS = [
    "SIM/h_ss_e11e33_inclusive_pT_22_24",
    "SIM/h_ss_e11e33_inclusive_pT_24_26",
    "SIM/h_ss_e11e33_inclusive_pT_26_28",
]


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ppg12-json", type=Path, default=PPG12_JSON)
    ap.add_argument("--photonjet-pointer", type=Path, default=PHOTONJET_POINTER)
    ap.add_argument("--photonjet-root", type=Path, default=None)
    ap.add_argument("--inclusive-root-glob", default=str(INCLUSIVE_ROOT_GLOB))
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    return ap.parse_args()


def load_current_photonjet_root(pointer: Path) -> tuple[Path, dict]:
    with pointer.open() as f:
        meta = json.load(f)
    roots = [Path(p) for p in meta.get("root_paths", [])]
    if len(roots) != 1:
        raise RuntimeError(f"Expected exactly one photonjet current root in {pointer}, found {roots}")
    if not roots[0].exists():
        raise FileNotFoundError(roots[0])
    return roots[0], meta


def add_hist_from_root(root_path: Path, hist_names: list[str]) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, float]]:
    values = None
    variances = None
    edges = None
    integrals: dict[str, float] = {}
    with uproot.open(root_path) as f:
        for name in hist_names:
            if name not in f:
                continue
            h = f[name]
            h_values = h.values(flow=False).astype(float)
            h_vars = h.variances(flow=False)
            if h_vars is None:
                h_vars = np.abs(h_values)
            h_vars = h_vars.astype(float)
            h_edges = h.axis(0).edges().astype(float)
            if edges is None:
                edges = h_edges
                values = np.zeros_like(h_values, dtype=float)
                variances = np.zeros_like(h_vars, dtype=float)
            elif len(edges) != len(h_edges) or not np.allclose(edges, h_edges):
                raise RuntimeError(f"Histogram binning mismatch in {root_path}:{name}")
            values += h_values
            variances += h_vars
            integrals[name] = float(np.sum(h_values))
    if values is None or variances is None or edges is None:
        raise RuntimeError(f"No requested histograms found in {root_path}: {hist_names}")
    return edges, values, variances, integrals


def load_current_sum(
    root_paths: list[Path], hist_names: list[str]
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, dict[str, float]], list[str]]:
    total_values = None
    total_vars = None
    total_edges = None
    per_root: dict[str, dict[str, float]] = {}
    skipped: list[str] = []
    for root in root_paths:
        try:
            edges, values, variances, integrals = add_hist_from_root(root, hist_names)
        except RuntimeError as exc:
            if "No requested histograms found" in str(exc):
                skipped.append(str(root))
                continue
            raise
        per_root[str(root)] = integrals
        if total_edges is None:
            total_edges = edges
            total_values = np.zeros_like(values, dtype=float)
            total_vars = np.zeros_like(variances, dtype=float)
        elif len(total_edges) != len(edges) or not np.allclose(total_edges, edges):
            raise RuntimeError(f"Current root binning mismatch: {root}")
        total_values += values
        total_vars += variances
    if total_values is None or total_vars is None or total_edges is None:
        raise RuntimeError("No current roots loaded")
    return total_edges, total_values, total_vars, per_root, skipped


def rebin_to_edges(
    source_edges: np.ndarray,
    source_values: np.ndarray,
    source_variances: np.ndarray,
    target_edges: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    rebinned = np.zeros(len(target_edges) - 1, dtype=float)
    rebinned_vars = np.zeros(len(target_edges) - 1, dtype=float)
    for ibin, (lo, hi, value, var) in enumerate(
        zip(source_edges[:-1], source_edges[1:], source_values, source_variances)
    ):
        if hi <= target_edges[0] or lo >= target_edges[-1]:
            continue
        matches = np.where(np.isclose(target_edges[:-1], lo) & np.isclose(target_edges[1:], hi))[0]
        if len(matches) == 1:
            j = int(matches[0])
            rebinned[j] += value
            rebinned_vars[j] += var
            continue
        contained = np.where((target_edges[:-1] <= lo + 1e-9) & (target_edges[1:] >= hi - 1e-9))[0]
        if len(contained) != 1:
            raise RuntimeError(f"Cannot map source bin {ibin}: [{lo}, {hi}] into target edges")
        j = int(contained[0])
        rebinned[j] += value
        rebinned_vars[j] += var
    return rebinned, rebinned_vars


def normalize(values: np.ndarray, variances: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    norm = float(np.sum(values))
    if norm <= 0:
        raise RuntimeError("Cannot normalize empty histogram")
    return values / norm, np.sqrt(np.maximum(variances, 0.0)) / norm


def ratio(numer: np.ndarray, numer_err: np.ndarray, denom: np.ndarray, denom_err: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    out = np.divide(numer, denom, out=np.full_like(numer, np.nan), where=denom > 0)
    err = np.full_like(out, np.nan)
    for i, (r, n, ne, d, de) in enumerate(zip(out, numer, numer_err, denom, denom_err)):
        if not np.isfinite(r) or n <= 0 or d <= 0:
            continue
        err[i] = abs(r) * math.sqrt((ne / n) ** 2 + (de / d) ** 2)
    return out, err


def draw_plot(payload: dict, current: dict, out_png: Path) -> None:
    centers = np.asarray(payload["signal"]["centers"], dtype=float)
    ppg12_signal = np.asarray(payload["signal"]["values"], dtype=float)
    ppg12_signal_err = np.asarray(payload["signal"]["errors"], dtype=float)
    ppg12_incl = np.asarray(payload["inclusive"]["values"], dtype=float)
    ppg12_incl_err = np.asarray(payload["inclusive"]["errors"], dtype=float)
    cur_signal = current["signal_values"]
    cur_signal_err = current["signal_errors"]
    cur_incl = current["inclusive_values"]
    cur_incl_err = current["inclusive_errors"]

    r_sig, r_sig_err = ratio(cur_signal, cur_signal_err, ppg12_signal, ppg12_signal_err)
    r_inc, r_inc_err = ratio(cur_incl, cur_incl_err, ppg12_incl, ppg12_incl_err)

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 13,
            "axes.linewidth": 1.15,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )
    fig, (ax, rax) = plt.subplots(
        2,
        1,
        figsize=(7.6, 7.8),
        dpi=160,
        sharex=True,
        gridspec_kw={"height_ratios": [3.1, 1.0], "hspace": 0.04},
    )

    red = "#e74c3c"
    blue = "#2266dd"
    ax.errorbar(
        centers,
        ppg12_signal,
        yerr=ppg12_signal_err,
        fmt="o",
        mfc="white",
        mec=red,
        ecolor=red,
        color=red,
        ms=4.5,
        lw=0.9,
        label="PPG12 Signal MC",
    )
    ax.errorbar(
        centers,
        cur_signal,
        yerr=cur_signal_err,
        fmt="o",
        mfc=red,
        mec=red,
        ecolor=red,
        color=red,
        ms=3.8,
        lw=0.9,
        label="Current Signal MC",
    )
    ax.errorbar(
        centers,
        ppg12_incl,
        yerr=ppg12_incl_err,
        fmt="s",
        mfc="white",
        mec=blue,
        ecolor=blue,
        color=blue,
        ms=4.3,
        lw=0.9,
        label="PPG12 Inclusive MC",
    )
    ax.errorbar(
        centers,
        cur_incl,
        yerr=cur_incl_err,
        fmt="s",
        mfc=blue,
        mec=blue,
        ecolor=blue,
        color=blue,
        ms=3.7,
        lw=0.9,
        label="Current Inclusive MC",
    )

    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 0.18)
    ax.set_ylabel("normalized counts")
    ax.minorticks_on()
    ax.tick_params(which="both", top=True, right=True)
    ax.text(0.055, 0.925, r"$\it{\bf{sPHENIX}}$ Internal", transform=ax.transAxes, fontsize=14)
    ax.text(0.055, 0.865, r"$p$+$p$ $\sqrt{s}=200$ GeV", transform=ax.transAxes, fontsize=11)
    ax.text(0.055, 0.810, r"$|\eta^\gamma|<0.7$", transform=ax.transAxes, fontsize=11)
    ax.text(0.055, 0.755, r"$22 < p_T < 28$ GeV", transform=ax.transAxes, fontsize=11)
    ax.text(0.055, 0.700, "before NPB / preselection", transform=ax.transAxes, fontsize=10)
    ax.legend(loc="upper right", frameon=False, fontsize=10.5, ncol=1, handlelength=1.4)

    rax.axhline(1.0, color="0.35", ls=(0, (4, 4)), lw=1.0)
    rax.errorbar(centers, r_sig, yerr=r_sig_err, fmt="o", mfc=red, mec=red, ecolor=red, color=red, ms=3.8, lw=0.8, label="Signal")
    rax.errorbar(centers, r_inc, yerr=r_inc_err, fmt="s", mfc=blue, mec=blue, ecolor=blue, color=blue, ms=3.7, lw=0.8, label="Inclusive")
    finite = np.concatenate([r_sig[np.isfinite(r_sig)], r_inc[np.isfinite(r_inc)]])
    if finite.size:
        lo = max(0.0, min(0.55, float(np.nanmin(finite)) - 0.08))
        hi = min(1.8, max(1.35, float(np.nanmax(finite)) + 0.10))
    else:
        lo, hi = 0.6, 1.4
    rax.set_ylim(lo, hi)
    rax.set_ylabel("Current / PPG12")
    rax.set_xlabel("e11_to_e33")
    rax.minorticks_on()
    rax.tick_params(which="both", top=True, right=True)
    rax.legend(loc="upper left", frameon=False, fontsize=10, ncol=2)

    fig.subplots_adjust(left=0.12, right=0.97, top=0.97, bottom=0.10)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png)
    plt.close(fig)


def write_csv(out_csv: Path, payload: dict, current: dict) -> None:
    centers = np.asarray(payload["signal"]["centers"], dtype=float)
    rows = []
    r_sig, r_sig_err = ratio(
        current["signal_values"],
        current["signal_errors"],
        np.asarray(payload["signal"]["values"], dtype=float),
        np.asarray(payload["signal"]["errors"], dtype=float),
    )
    r_inc, r_inc_err = ratio(
        current["inclusive_values"],
        current["inclusive_errors"],
        np.asarray(payload["inclusive"]["values"], dtype=float),
        np.asarray(payload["inclusive"]["errors"], dtype=float),
    )
    for i, x in enumerate(centers):
        rows.append(
            {
                "center": x,
                "ppg12_signal": payload["signal"]["values"][i],
                "current_signal": current["signal_values"][i],
                "current_over_ppg12_signal": r_sig[i],
                "current_over_ppg12_signal_err": r_sig_err[i],
                "ppg12_inclusive": payload["inclusive"]["values"][i],
                "current_inclusive": current["inclusive_values"][i],
                "current_over_ppg12_inclusive": r_inc[i],
                "current_over_ppg12_inclusive_err": r_inc_err[i],
            }
        )
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    with args.ppg12_json.open() as f:
        payload = json.load(f)
    target_edges = np.asarray(payload["signal"]["edges"], dtype=float)
    centers = np.asarray(payload["signal"]["centers"], dtype=float)

    if args.photonjet_root is not None:
        photonjet_root = args.photonjet_root
        if not photonjet_root.exists():
            raise FileNotFoundError(photonjet_root)
        photonjet_meta = {
            "artifact_kind": "explicit_photonjet_root_override",
            "root": str(photonjet_root),
            "production_artifact_exception": "Explicit root override used for complete 22-28 diagnostic.",
        }
    else:
        photonjet_root, photonjet_meta = load_current_photonjet_root(args.photonjet_pointer)
    inclusive_roots = sorted(Path(p) for p in glob.glob(args.inclusive_root_glob) if Path(p).exists())
    if not inclusive_roots:
        raise RuntimeError(f"No inclusive roots matched {args.inclusive_root_glob}")

    sig_edges, sig_values_raw, sig_vars_raw, sig_integrals, sig_skipped = load_current_sum([photonjet_root], PT_HISTS)
    inc_edges, inc_values_raw, inc_vars_raw, inc_integrals, inc_skipped = load_current_sum(inclusive_roots, PT_HISTS)

    sig_rebin, sig_var_rebin = rebin_to_edges(sig_edges, sig_values_raw, sig_vars_raw, target_edges)
    inc_rebin, inc_var_rebin = rebin_to_edges(inc_edges, inc_values_raw, inc_vars_raw, target_edges)
    sig_norm, sig_err = normalize(sig_rebin, sig_var_rebin)
    inc_norm, inc_err = normalize(inc_rebin, inc_var_rebin)

    current = {
        "centers": centers,
        "signal_values": sig_norm,
        "signal_errors": sig_err,
        "inclusive_values": inc_norm,
        "inclusive_errors": inc_err,
    }

    out_png = args.out_dir / "fig13_e11e33_signal_inclusive_mc_ppg12_vs_current_overlay.png"
    out_csv = args.out_dir / "fig13_e11e33_signal_inclusive_mc_ppg12_vs_current_bins.csv"
    out_manifest = args.out_dir / "fig13_e11e33_signal_inclusive_mc_ppg12_vs_current_manifest.json"

    draw_plot(payload, current, out_png)
    write_csv(out_csv, payload, current)

    manifest = {
        "schema": "the76_fig13_e11e33_signal_inclusive_mc_overlay_v1",
        "script": str(Path(__file__).resolve()),
        "ppg12_json": str(args.ppg12_json),
        "ppg12_source_paths": payload.get("source_paths", {}),
        "photonjet_current_pointer": str(args.photonjet_pointer),
        "photonjet_current_meta": photonjet_meta,
        "photonjet_root_override": str(args.photonjet_root) if args.photonjet_root else None,
        "photonjet_root": str(photonjet_root),
        "inclusive_roots": [str(p) for p in inclusive_roots],
        "inclusive_current_exception": (
            "No promoted current inclusive pointer exists; using July-2 strict inclusive final roots. "
            "Replace with current pointer when inclusive rerun is promoted."
        ),
        "histograms_summed": PT_HISTS,
        "normalization": "unit area over PPG12 Fig.13 0.00-1.00 E11/E33 binning after current 0.01->0.04 rebin",
        "output_png": str(out_png),
        "output_csv": str(out_csv),
        "signal_raw_integrals": sig_integrals,
        "signal_skipped_roots": sig_skipped,
        "inclusive_raw_integrals": inc_integrals,
        "inclusive_skipped_roots": inc_skipped,
    }
    args.out_dir.mkdir(parents=True, exist_ok=True)
    with out_manifest.open("w") as f:
        json.dump(manifest, f, indent=2, sort_keys=True)

    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
