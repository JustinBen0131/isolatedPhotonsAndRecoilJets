#!/usr/bin/env python3
"""Make first-pass THE-76 PPG12 photon-yield SIM comparison plots."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import ROOT
from matplotlib.colors import LogNorm


ROOT.gROOT.SetBatch(True)

PPG12_RECO_BINS = [10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36]
PPG12_TRUTH_BINS = [8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36, 45]

DEFAULT_BASE = Path(
    "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/ppg12PhotonYield/"
    "THE76_ppg12_photon_yield_v1_sim_20260616"
)
DEFAULT_CFG = (
    "jetMinPtScan_dphiScan_vz60_isoR40_isSliding_"
    "preselectionNewPPG12_tightNewPPG12_nonTightNewPPG12"
)
DEFAULT_SIGNAL = (
    DEFAULT_BASE
    / "merged_roots"
    / DEFAULT_CFG
    / "photonJet5and10and20merged_SIM"
    / "RecoilJets_photonjet5plus10plus20_MERGED.root"
)
DEFAULT_INCLUSIVE = (
    DEFAULT_BASE
    / "merged_roots"
    / DEFAULT_CFG
    / "inclusiveJet5to40_SIM"
    / "RecoilJets_jet5plus8plus12plus20plus30plus40_MERGED.root"
)
DEFAULT_OUTDIR = DEFAULT_BASE / "ppg12_ian_sim_check_plots"


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "dejavuserif",
            "axes.linewidth": 1.2,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "legend.frameon": False,
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )


def open_root(path: Path) -> ROOT.TFile:
    f = ROOT.TFile.Open(str(path))
    if not f or f.IsZombie():
        raise RuntimeError(f"could not open ROOT file: {path}")
    return f


def get_hist(f: ROOT.TFile, name: str):
    hist = f.Get(name)
    if not hist:
        raise RuntimeError(f"missing histogram {name} in {f.GetName()}")
    return hist


def axis_edges(axis: ROOT.TAxis) -> np.ndarray:
    bins = axis.GetXbins()
    if bins.GetSize() == axis.GetNbins() + 1:
        return np.array([bins.At(i) for i in range(bins.GetSize())], dtype=float)
    return np.array(
        [axis.GetBinLowEdge(i) for i in range(1, axis.GetNbins() + 2)],
        dtype=float,
    )


def hist1_arrays(hist) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    edges = axis_edges(hist.GetXaxis())
    values = np.array([hist.GetBinContent(i) for i in range(1, hist.GetNbinsX() + 1)], dtype=float)
    errors = np.array([hist.GetBinError(i) for i in range(1, hist.GetNbinsX() + 1)], dtype=float)
    centers = 0.5 * (edges[:-1] + edges[1:])
    return centers, values, errors


def hist2_arrays(hist) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    xedges = axis_edges(hist.GetXaxis())
    yedges = axis_edges(hist.GetYaxis())
    z = np.zeros((hist.GetNbinsY(), hist.GetNbinsX()), dtype=float)
    for iy in range(1, hist.GetNbinsY() + 1):
        for ix in range(1, hist.GetNbinsX() + 1):
            z[iy - 1, ix - 1] = hist.GetBinContent(ix, iy)
    return xedges, yedges, z


def sphinx_label(ax, x=0.05, y=0.95, size=13) -> None:
    ax.text(
        x,
        y,
        r"$\it{\bf{sPHENIX}}$ Internal",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=size,
    )


def save(fig, out: Path) -> Path:
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return out


def plot_response(signal: ROOT.TFile, outdir: Path) -> dict:
    hist = get_hist(signal, "SIM/h_response_full_0")
    xedges, yedges, z = hist2_arrays(hist)
    positive = z[z > 0]
    fig, ax = plt.subplots(figsize=(7.0, 5.7))
    mesh = ax.pcolormesh(
        xedges,
        yedges,
        np.ma.masked_less_equal(z, 0),
        cmap="viridis",
        norm=LogNorm(vmin=max(float(positive.min()), 1e-12), vmax=float(positive.max())),
    )
    fig.colorbar(mesh, ax=ax, pad=0.02, label="matched photon counts")
    ax.plot([10, 36], [10, 36], color="white", lw=1.6, ls="--", alpha=0.9)
    ax.set_xlabel(r"reconstructed photon $E_T$ [GeV]")
    ax.set_ylabel(r"truth photon $E_T$ [GeV]")
    ax.set_title("THE-76 pp signal MC photon-yield response", fontsize=15, weight="bold")
    ax.set_xlim(10, 36)
    ax.set_ylim(8, 45)
    sphinx_label(ax, size=12.5)
    ax.text(
        0.98,
        0.05,
        "PPG12 binning\nsliding iso R=0.4",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=12,
    )
    path = save(fig, outdir / "the76_signal_response_matrix_h_response_full_0.png")
    return {
        "path": str(path),
        "hist": "SIM/h_response_full_0",
        "entries": float(hist.GetEntries()),
        "integral": float(hist.Integral()),
        "x_edges": list(map(float, xedges)),
        "y_edges": list(map(float, yedges)),
    }


def plot_response_spectra(signal: ROOT.TFile, outdir: Path) -> dict:
    specs = [
        ("SIM/h_truth_pT_0", "truth denominator", "#555555"),
        ("SIM/h_pT_truth_response_0", "truth matched", "#0072B2"),
        ("SIM/h_pT_reco_response_0", "reco matched", "#009E73"),
        ("SIM/h_pT_reco_fake_0", "reco fake", "#D55E00"),
    ]
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    meta = {}
    for name, label, color in specs:
        hist = get_hist(signal, name)
        x, y, err = hist1_arrays(hist)
        mask = y > 0
        ax.errorbar(x[mask], y[mask], yerr=err[mask], fmt="o", ms=4.8, lw=1.1, color=color, label=label)
        meta[name] = {"entries": float(hist.GetEntries()), "integral": float(hist.Integral())}
    ax.set_yscale("log")
    ax.set_xlabel(r"photon $E_T$ [GeV]")
    ax.set_ylabel("counts")
    ax.set_title("Signal MC response ingredients", fontsize=15, weight="bold")
    ax.set_xlim(8, 45)
    ax.grid(True, which="major", axis="y", alpha=0.18)
    sphinx_label(ax, size=12.5)
    ax.legend(loc="upper right")
    path = save(fig, outdir / "the76_signal_response_spectra.png")
    meta["path"] = str(path)
    return meta


def plot_abcd_inputs(signal: ROOT.TFile, inclusive: ROOT.TFile, outdir: Path) -> dict:
    regions = [
        ("A tight + isolated", "h_tight_iso_cluster_0", "#0072B2"),
        ("B tight + non-isolated", "h_tight_noniso_cluster_0", "#D55E00"),
        ("C non-tight + isolated", "h_nontight_iso_cluster_0", "#009E73"),
        ("D non-tight + non-isolated", "h_nontight_noniso_cluster_0", "#CC79A7"),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.9), sharey=False)
    meta = {}
    for ax, f, role in [(axes[0], signal, "signal MC"), (axes[1], inclusive, "inclusive MC")]:
        for label, base, color in regions:
            hist = get_hist(f, f"SIM/{base}")
            x, y, err = hist1_arrays(hist)
            mask = y > 0
            ax.errorbar(x[mask], y[mask], yerr=err[mask], fmt="o", ms=4.2, lw=1.0, color=color, label=label)
            meta[f"{role}:{base}"] = {"entries": float(hist.GetEntries()), "integral": float(hist.Integral())}
        ax.set_yscale("log")
        ax.set_xlabel(r"cluster $E_T$ [GeV]")
        ax.set_title(role, fontsize=13.5, weight="bold")
        ax.grid(True, which="major", axis="y", alpha=0.18)
    axes[0].set_ylabel("ABCD-region counts")
    axes[0].legend(loc="upper right", fontsize=10.5)
    fig.suptitle("THE-76 PPG12 photon-yield ABCD inputs", fontsize=16, weight="bold")
    path = save(fig, outdir / "the76_abcd_region_inputs_signal_vs_inclusive.png")
    meta["path"] = str(path)
    return meta


def plot_leakage(signal: ROOT.TFile, outdir: Path) -> dict:
    a = get_hist(signal, "SIM/h_tight_iso_cluster_signal_0")
    denominators = hist1_arrays(a)
    x = denominators[0]
    avals = denominators[1]
    curves = [
        ("B/A signal leakage", "SIM/h_tight_noniso_cluster_signal_0", "#D55E00"),
        ("C/A signal leakage", "SIM/h_nontight_iso_cluster_signal_0", "#009E73"),
        ("D/A signal leakage", "SIM/h_nontight_noniso_cluster_signal_0", "#CC79A7"),
    ]
    rows = []
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    for label, name, color in curves:
        hist = get_hist(signal, name)
        _, vals, _ = hist1_arrays(hist)
        frac = np.divide(vals, avals, out=np.full_like(vals, np.nan), where=avals > 0)
        mask = np.isfinite(frac)
        ax.plot(x[mask], frac[mask], marker="o", lw=2.0, ms=4.8, color=color, label=label)
        for xi, yi in zip(x[mask], frac[mask]):
            rows.append({"curve": label, "et_center": float(xi), "fraction": float(yi)})
    ax.set_xlabel(r"cluster $E_T$ [GeV]")
    ax.set_ylabel("signal leakage fraction relative to region A")
    ax.set_ylim(bottom=0)
    ax.set_title("Signal leakage inputs for ABCD purity correction", fontsize=15, weight="bold")
    ax.grid(True, axis="y", alpha=0.2)
    sphinx_label(ax, size=12.5)
    ax.legend(loc="upper right")
    png = save(fig, outdir / "the76_signal_leakage_fractions_relative_to_A.png")
    csv_path = outdir / "the76_signal_leakage_fractions_relative_to_A.csv"
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["curve", "et_center", "fraction"])
        writer.writeheader()
        writer.writerows(rows)
    return {"path": str(png), "csv": str(csv_path), "rows": len(rows)}


def normalize_hist(hist) -> tuple[np.ndarray, np.ndarray]:
    x, y, _ = hist1_arrays(hist)
    total = float(np.sum(y))
    if total > 0:
        y = y / total
    return x, y


def plot_shower_shapes(signal: ROOT.TFile, inclusive: ROOT.TFile, outdir: Path) -> dict:
    pt = "20_22"
    specs = [
        (signal, f"SIM/h_ss_e11e33_pre_sig_pT_{pt}", "signal MC preselection", "#D55E00", "-"),
        (inclusive, f"SIM/h_ss_e11e33_pre_bkg_pT_{pt}", "inclusive MC preselection", "#0072B2", "-"),
        (signal, f"SIM/h_ss_e11e33_tight_sig_pT_{pt}", "signal MC tight ID", "#D55E00", "--"),
        (inclusive, f"SIM/h_ss_e11e33_tight_bkg_pT_{pt}", "inclusive MC tight ID", "#0072B2", "--"),
    ]
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    meta = {}
    for f, name, label, color, ls in specs:
        hist = get_hist(f, name)
        x, y = normalize_hist(hist)
        ax.step(x, y, where="mid", lw=2.1, color=color, ls=ls, label=label)
        low_mask = x < 0.05
        meta[name] = {
            "entries": float(hist.GetEntries()),
            "normalized_low_e11e33_lt_0p05": float(np.sum(y[low_mask])),
        }
    ax.set_xlabel(r"$E_{1\times1}/E_{3\times3}$")
    ax.set_ylabel("unit-normalized candidates")
    ax.set_xlim(0, 1.08)
    ax.set_ylim(bottom=0)
    ax.set_title(r"Shower-shape sanity: $20<E_T<22$ GeV", fontsize=15, weight="bold")
    ax.grid(True, axis="y", alpha=0.18)
    ax.text(
        0.98,
        0.08,
        r"$\it{\bf{sPHENIX}}$ Internal",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=12.5,
    )
    ax.legend(loc="upper left", fontsize=10.5)
    png = save(fig, outdir / "the76_shower_shape_e11e33_signal_inclusive_pre_tight_20_22.png")
    meta["path"] = str(png)
    return meta


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--signal-root", type=Path, default=DEFAULT_SIGNAL)
    ap.add_argument("--inclusive-root", type=Path, default=DEFAULT_INCLUSIVE)
    ap.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = ap.parse_args()

    setup_style()
    args.outdir.mkdir(parents=True, exist_ok=True)
    signal = open_root(args.signal_root)
    inclusive = open_root(args.inclusive_root)

    manifest = {
        "signal_root": str(args.signal_root.resolve()),
        "inclusive_root": str(args.inclusive_root.resolve()),
        "ppg12_reco_bins": PPG12_RECO_BINS,
        "ppg12_truth_bins": PPG12_TRUTH_BINS,
        "plots": {
            "response": plot_response(signal, args.outdir),
            "response_spectra": plot_response_spectra(signal, args.outdir),
            "abcd_inputs": plot_abcd_inputs(signal, inclusive, args.outdir),
            "signal_leakage": plot_leakage(signal, args.outdir),
            "shower_shapes": plot_shower_shapes(signal, inclusive, args.outdir),
        },
        "truthfulness_note": (
            "These are THE-76 SIM-only PPG12 photon-yield contract plots from the final stitched "
            "signal and inclusive MC ROOTs; they are not pp data purity results."
        ),
    }
    manifest_path = args.outdir / "the76_ppg12_photon_yield_sim_check_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2))
    print(manifest_path)
    for plot in manifest["plots"].values():
        if "path" in plot:
            print(plot["path"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
