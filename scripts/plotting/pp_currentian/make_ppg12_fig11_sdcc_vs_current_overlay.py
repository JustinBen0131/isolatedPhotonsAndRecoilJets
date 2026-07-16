#!/usr/bin/env python3
"""Overlay PPG12 IAN Fig.11 SDCC values with the promoted current SIM output."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import ROOT


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
REFERENCE_DIR = REPO / "dataOutput/ppg12Parity/ppg12_fig11_sdcc_reference"
REFERENCE_CSV = REFERENCE_DIR / "fig11_sb_sdcc_verbatim_points.csv"
PHOTON_POINTER = REPO / "dataOutput/current_recoiljets_artifacts/current/pp_sim_photonjet_merged/current.json"
INCLUSIVE_POINTER = REPO / "dataOutput/current_recoiljets_artifacts/current/pp_sim_inclusivejet_merged/current.json"
ISO_A = 0.502095
ISO_B = 0.0433036
JET_TO_PHOTON20 = 7.3113 / 130.4461
SIGNAL_HIST = "h_ppg12_fig11_ET_isoET_eta0_signal_mc"
BACKGROUND_HIST = "h_ppg12_fig11_ET_isoET_eta0_jet_inclusive_mc"


def pointer_root(path: Path) -> tuple[Path, dict]:
    payload = json.loads(path.read_text())
    roots = payload.get("root_paths") or []
    if len(roots) != 1:
        raise RuntimeError(f"Expected exactly one ROOT path in {path}")
    return Path(roots[0]), payload


def read_reference(path: Path) -> list[dict[str, float]]:
    with path.open() as handle:
        rows = [{key: float(value) for key, value in row.items()} for row in csv.DictReader(handle)]
    if len(rows) != 10:
        raise RuntimeError(f"Expected ten Fig.11 SDCC rows, found {len(rows)}")
    return rows


def projected_sum(hist: ROOT.TH2, xbin: int, xcenter: float) -> tuple[float, float]:
    ylow = hist.GetYaxis().FindBin(-1.0)
    yhigh = hist.GetYaxis().FindBin(ISO_A + ISO_B * xcenter)
    value = sum(hist.GetBinContent(xbin, iy) for iy in range(ylow, yhigh + 1))
    variance = sum(hist.GetBinError(xbin, iy) ** 2 for iy in range(ylow, yhigh + 1))
    return value, math.sqrt(variance)


def current_rows(photon_root: Path, inclusive_root: Path, reference: list[dict[str, float]]) -> list[dict[str, float]]:
    photon_file = ROOT.TFile.Open(str(photon_root), "READ")
    inclusive_file = ROOT.TFile.Open(str(inclusive_root), "READ")
    if not photon_file or photon_file.IsZombie() or not inclusive_file or inclusive_file.IsZombie():
        raise RuntimeError("Could not open one of the promoted current SIM ROOTs")
    signal = photon_file.Get(f"SIM/{SIGNAL_HIST}")
    background = inclusive_file.Get(f"SIM/{BACKGROUND_HIST}")
    if not signal or not background:
        raise RuntimeError("Missing the dedicated PPG12 Fig.11 current histogram family")
    signal = signal.Clone("fig11_current_signal")
    background = background.Clone("fig11_current_background")
    signal.SetDirectory(0)
    background.SetDirectory(0)
    photon_file.Close()
    inclusive_file.Close()
    signal.RebinX(16)
    background.RebinX(16)

    rows = []
    for ref in reference:
        xcenter = ref["xcenter"]
        signal_bin = signal.GetXaxis().FindBin(xcenter)
        background_bin = background.GetXaxis().FindBin(xcenter)
        sig, sig_err = projected_sum(signal, signal_bin, xcenter)
        bkg_raw, bkg_raw_err = projected_sum(background, background_bin, xcenter)
        bkg = bkg_raw * JET_TO_PHOTON20
        bkg_err = bkg_raw_err * JET_TO_PHOTON20
        if sig <= 0.0 or bkg <= 0.0:
            raise RuntimeError(f"Nonpositive current Fig.11 projection at ET={xcenter:g} GeV")
        sb = sig / bkg
        sb_err = sb * math.sqrt((sig_err / sig) ** 2 + (bkg_err / bkg) ** 2)
        rows.append(
            {
                "xlow": ref["xlow"],
                "xhigh": ref["xhigh"],
                "xcenter": xcenter,
                "iso_max": ISO_A + ISO_B * xcenter,
                "signal": sig,
                "signal_err": sig_err,
                "background_raw": bkg_raw,
                "background_raw_err": bkg_raw_err,
                "background_scaled": bkg,
                "background_scaled_err": bkg_err,
                "s_over_b": sb,
                "s_over_b_err": sb_err,
            }
        )
    return rows


def write_points(path: Path, reference: list[dict[str, float]], current: list[dict[str, float]]) -> None:
    fields = [
        "xlow", "xhigh", "xcenter", "iso_max", "sdcc_s_over_b", "sdcc_s_over_b_err",
        "current_s_over_b", "current_s_over_b_err", "current_over_sdcc", "current_over_sdcc_err",
        "current_signal", "current_signal_err", "current_background_scaled", "current_background_scaled_err",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for ref, cur in zip(reference, current, strict=True):
            ratio = cur["s_over_b"] / ref["s_over_b"]
            ratio_err = ratio * math.sqrt(
                (cur["s_over_b_err"] / cur["s_over_b"]) ** 2
                + (ref["s_over_b_err"] / ref["s_over_b"]) ** 2
            )
            writer.writerow(
                {
                    "xlow": ref["xlow"], "xhigh": ref["xhigh"], "xcenter": ref["xcenter"],
                    "iso_max": ref["iso_max"], "sdcc_s_over_b": ref["s_over_b"],
                    "sdcc_s_over_b_err": ref["s_over_b_err"], "current_s_over_b": cur["s_over_b"],
                    "current_s_over_b_err": cur["s_over_b_err"], "current_over_sdcc": ratio,
                    "current_over_sdcc_err": ratio_err, "current_signal": cur["signal"],
                    "current_signal_err": cur["signal_err"], "current_background_scaled": cur["background_scaled"],
                    "current_background_scaled_err": cur["background_scaled_err"],
                }
            )


def draw(path: Path, reference: list[dict[str, float]], current: list[dict[str, float]]) -> dict[str, float]:
    plt.rcParams.update({"font.family": "serif", "font.serif": ["Times New Roman", "Times", "DejaVu Serif"], "mathtext.fontset": "dejavuserif", "axes.linewidth": 1.25})
    fig, (ax, rax) = plt.subplots(2, 1, figsize=(7.2, 8.7), dpi=220, sharex=True, gridspec_kw={"height_ratios": [3.3, 1.0], "hspace": 0.05})
    x = np.asarray([r["xcenter"] for r in reference])
    xerr = np.asarray([[r["xcenter"] - r["xlow"] for r in reference], [r["xhigh"] - r["xcenter"] for r in reference]])
    y_ref = np.asarray([r["s_over_b"] for r in reference])
    e_ref = np.asarray([r["s_over_b_err"] for r in reference])
    y_cur = np.asarray([r["s_over_b"] for r in current])
    e_cur = np.asarray([r["s_over_b_err"] for r in current])
    ratio = y_cur / y_ref
    ratio_err = ratio * np.hypot(e_cur / y_cur, e_ref / y_ref)
    ax.errorbar(x, y_ref, xerr=xerr, yerr=e_ref, fmt="o", ms=7.2, mfc="white", mec="black", mew=1.5, ecolor="black", elinewidth=1.0, capsize=0, label="PPG12 SDCC source", zorder=4)
    ax.errorbar(x, y_cur, xerr=xerr, yerr=e_cur, fmt="o", ms=5.2, mfc="black", mec="black", ecolor="black", elinewidth=1.0, capsize=0, label="Current output", zorder=5)
    rax.errorbar(x, ratio, xerr=xerr, yerr=ratio_err, fmt="o", ms=5.0, mfc="black", mec="black", ecolor="black", elinewidth=1.0, capsize=0, zorder=3)
    ax.set_xlim(10.0, 32.0)
    ax.set_ylim(0.0, 0.82)
    ax.set_ylabel("S/B", fontsize=17)
    ax.minorticks_on()
    ax.tick_params(which="both", direction="in", top=True, right=True, labelsize=12)
    ax.text(0.055, 0.93, "sPHENIX", transform=ax.transAxes, fontsize=14, fontstyle="italic", fontweight="bold")
    ax.text(0.275, 0.93, "Internal", transform=ax.transAxes, fontsize=14)
    ax.text(0.055, 0.855, r"$p$+$p$ $\sqrt{s}=200$ GeV", transform=ax.transAxes, fontsize=12)
    ax.text(0.055, 0.79, "PYTHIA8", transform=ax.transAxes, fontsize=12)
    ax.legend(handles=[Line2D([0], [0], marker="o", color="black", mfc="white", mec="black", mew=1.5, lw=0, ms=7, label="PPG12 SDCC source"), Line2D([0], [0], marker="o", color="black", mfc="black", mec="black", lw=0, ms=5.5, label="Current output")], loc="lower right", frameon=False, fontsize=12, handlelength=1.0)
    rax.axhline(1.0, color="0.4", lw=1.0, ls=(0, (4, 4)))
    rax.set_ylim(0.65, 1.95)
    rax.set_ylabel("Current /\nPPG12 SDCC", fontsize=12.5)
    rax.set_xlabel(r"$E_{\mathrm{T}}^{\gamma,\mathrm{rec}}$ [GeV]", fontsize=16, ha="right", x=1.0)
    rax.set_xticks(np.arange(10.0, 33.0, 2.0))
    rax.minorticks_on()
    rax.tick_params(which="both", direction="in", top=True, right=True, labelsize=11)
    fig.subplots_adjust(left=0.16, right=0.97, top=0.98, bottom=0.10)
    fig.savefig(path)
    plt.close(fig)
    return {"min_current_over_sdcc": float(ratio.min()), "max_current_over_sdcc": float(ratio.max()), "mean_current_over_sdcc": float(ratio.mean())}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-csv", type=Path, default=REFERENCE_CSV)
    parser.add_argument("--photon-root", type=Path, default=None)
    parser.add_argument("--inclusive-root", type=Path, default=None)
    parser.add_argument("--out-dir", type=Path, default=None)
    args = parser.parse_args()
    photon_root, photon_meta = pointer_root(PHOTON_POINTER) if args.photon_root is None else (args.photon_root, {"resolution": "explicit"})
    inclusive_root, inclusive_meta = pointer_root(INCLUSIVE_POINTER) if args.inclusive_root is None else (args.inclusive_root, {"resolution": "explicit"})
    campaign = photon_meta.get("campaign_tag", "current")
    out_dir = args.out_dir or REPO / "dataOutput/ppg12Parity" / campaign / "fig11_sb"
    out_dir.mkdir(parents=True, exist_ok=True)
    png = out_dir / "ppg12_ian_fig11_sb_sdcc_vs_current_overlay_ratio.png"
    csv_path = out_dir / "ppg12_ian_fig11_sb_sdcc_vs_current_overlay_points.csv"
    manifest = out_dir / "ppg12_ian_fig11_sb_sdcc_vs_current_overlay_manifest.json"
    reference = read_reference(args.reference_csv)
    current = current_rows(photon_root, inclusive_root, reference)
    write_points(csv_path, reference, current)
    summary = draw(png, reference, current)
    manifest.write_text(json.dumps({"artifact": str(png), "points_csv": str(csv_path), "ppg12_reference_csv": str(args.reference_csv), "ppg12_reference_csv_sha256": hashlib.sha256(args.reference_csv.read_bytes()).hexdigest(), "photon_pointer": str(PHOTON_POINTER) if args.photon_root is None else "explicit --photon-root", "inclusive_pointer": str(INCLUSIVE_POINTER) if args.inclusive_root is None else "explicit --inclusive-root", "photon_root": str(photon_root), "inclusive_root": str(inclusive_root), "photon_artifact_promotion_basis": photon_meta.get("promotion_basis", ""), "inclusive_artifact_promotion_basis": inclusive_meta.get("promotion_basis", ""), "current_histograms": {"signal": f"SIM/{SIGNAL_HIST}", "background": f"SIM/{BACKGROUND_HIST}"}, "construction": "Both source and current use PPG12 Fig.11 RebinX(16), -1 < isoET < 0.502095+0.0433036*ET, and jet50/photon20=7.3113/130.4461 background conversion. Current ratio errors propagate independent weighted projection errors.", "ratio_panel": "Current output / PPG12 SDCC source", "summary": summary}, indent=2, sort_keys=True) + "\n")
    print(png)
    print(csv_path)
    print(manifest)
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
