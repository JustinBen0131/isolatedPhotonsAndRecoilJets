#!/usr/bin/env python3
"""Render THE-42 top-10 isolation distributions on a linear y-axis."""

from __future__ import annotations

import argparse
import json
import re
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


REPO = next((p for p in Path(__file__).resolve().parents if (p / "AGENTS.md").exists()), Path.cwd())
DEFAULT_BASE = REPO / "dataOutput/auauTightBDTValidation/THE42_wp80_centlinear_ss_20260604"
DEFAULT_OUTDIR = DEFAULT_BASE / "top10_coherency_diagnostics_20260606"
DEFAULT_DATA = (
    DEFAULT_BASE
    / "top10_merged_data_roots_20260606"
    / "RecoilJets_auau_ALL_preselectionNewPPG12_tightAuAuCentInputBase3x3BDT_nonTightAuAuBDTComplement_baseVariant.root"
)
DEFAULT_SIGNAL = (
    DEFAULT_BASE
    / "sim_roots/preselectionNewPPG12_tightAuAuCentInputBase3x3BDT_nonTightAuAuBDTComplement_baseVariant"
    / "photonJet12and20merged_SIM/RecoilJets_embeddedPhoton12plus20_MERGED.root"
)
DEFAULT_INCLUSIVE = (
    DEFAULT_BASE
    / "sim_roots/preselectionNewPPG12_tightAuAuCentInputBase3x3BDT_nonTightAuAuBDTComplement_baseVariant"
    / "embeddedJet12and20and30and40merged_SIM/RecoilJets_embeddedJet12plus20plus30plus40_MERGED.root"
)

PT_BINS = ((22, 24), (24, 26), (26, 28))
CENT_BINS = ((0, 10), (10, 20), (20, 30), (30, 40), (40, 50), (50, 60), (60, 80))
ISO_VARS = (
    ("h_Eiso_isoR30", "Total isolation R=0.3"),
    ("h_Eiso_emcal_isoR30", "EMCal isolation R=0.3"),
    ("h_Eiso_hcalin_isoR30", "IHCal isolation R=0.3"),
    ("h_Eiso_hcalout_isoR30", "OHCal isolation R=0.3"),
)


@dataclass
class Curve:
    centers: np.ndarray
    edges: np.ndarray
    values: np.ndarray
    errors: np.ndarray
    raw_integral: float
    entries: float
    matches: int


def open_root(path: Path):
    import ROOT

    ROOT.gROOT.SetBatch(True)
    ROOT.TH1.AddDirectory(False)
    f = ROOT.TFile.Open(str(path), "READ")
    if not f or f.IsZombie():
        raise OSError(f"Could not open ROOT file: {path}")
    return f


def walk_root_dir(directory, prefix: str = ""):
    for key in directory.GetListOfKeys():
        obj = key.ReadObj()
        name = key.GetName()
        path = f"{prefix}/{name}" if prefix else name
        yield path, obj
        if obj.InheritsFrom("TDirectory"):
            yield from walk_root_dir(obj, path)


def add_hist(acc, hist):
    if acc is None:
        out = hist.Clone(f"{hist.GetName()}_sum")
        out.SetDirectory(0)
        return out
    acc.Add(hist)
    return acc


def matching_hist_names(base_name: str) -> set[str]:
    names = set()
    for pt_lo, pt_hi in PT_BINS:
        for c_lo, c_hi in CENT_BINS:
            names.add(f"{base_name}_pT_{pt_lo}_{pt_hi}_cent_{c_lo}_{c_hi}")
    return names


def sum_hists(root_path: Path, base_name: str, trigger_regex: re.Pattern[str] | None) -> tuple[object, int]:
    import ROOT  # noqa: F401

    target_names = matching_hist_names(base_name)
    acc = None
    matches = 0
    f = open_root(root_path)
    try:
        for object_path, obj in walk_root_dir(f):
            if not obj.InheritsFrom("TH1"):
                continue
            if Path(object_path).name not in target_names:
                continue
            if trigger_regex and not trigger_regex.search(object_path):
                continue
            h = obj.Clone(f"{base_name}_{matches}")
            h.SetDirectory(0)
            acc = add_hist(acc, h)
            matches += 1
    finally:
        f.Close()
    if acc is None:
        raise RuntimeError(f"No histograms found for {base_name} in {root_path}")
    return acc, matches


def curve_from_hist(hist, matches: int, x_min: float, x_max: float) -> Curve:
    nbins = hist.GetNbinsX()
    edges = np.array([hist.GetXaxis().GetBinLowEdge(i) for i in range(1, nbins + 2)], dtype=float)
    centers = np.array([hist.GetXaxis().GetBinCenter(i) for i in range(1, nbins + 1)], dtype=float)
    counts = np.array([hist.GetBinContent(i) for i in range(1, nbins + 1)], dtype=float)
    errs = np.array([hist.GetBinError(i) for i in range(1, nbins + 1)], dtype=float)
    mask = (centers >= x_min) & (centers < x_max)
    idx = np.flatnonzero(mask)
    if len(idx):
        edges = np.concatenate(([edges[idx[0]]], edges[idx + 1]))
        centers = centers[mask]
        counts = counts[mask]
        errs = errs[mask]
    raw_integral = float(np.sum(counts))
    values = np.zeros_like(counts)
    errors = np.zeros_like(errs)
    if raw_integral > 0:
        values = counts / raw_integral
        errors = errs / raw_integral
    return Curve(
        centers=centers,
        edges=edges,
        values=values,
        errors=errors,
        raw_integral=raw_integral,
        entries=float(hist.GetEntries()),
        matches=matches,
    )


def load_curves(args: argparse.Namespace) -> dict[str, dict[str, Curve]]:
    data_trigger_regex = re.compile(args.trigger_regex) if args.trigger_regex else None
    roots = {
        "Signal MC": Path(args.signal_root),
        "Inclusive MC": Path(args.inclusive_root),
        "Data top-10": Path(args.data_root),
    }
    curves: dict[str, dict[str, Curve]] = {}
    for base_name, _label in ISO_VARS:
        curves[base_name] = {}
        for sample, root_path in roots.items():
            trigger_regex = data_trigger_regex if sample == "Data top-10" else None
            hist, matches = sum_hists(root_path, base_name, trigger_regex)
            curves[base_name][sample] = curve_from_hist(hist, matches, args.x_min, args.x_max)
    return curves


def setup_style() -> None:
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "axes.edgecolor": "#111827",
        "axes.linewidth": 1.0,
    })


def render(curves: dict[str, dict[str, Curve]], outdir: Path) -> Path:
    colors = {"Signal MC": "#E32D26", "Inclusive MC": "#1F67D2", "Data top-10": "#111111"}
    fig, axes = plt.subplots(2, 2, figsize=(14.8, 8.8), dpi=170, sharex=True)
    for ax, (base_name, title) in zip(axes.flat, ISO_VARS):
        local_max = 0.0
        for sample in ("Signal MC", "Inclusive MC"):
            c = curves[base_name][sample]
            ax.stairs(c.values, c.edges, color=colors[sample], linewidth=1.9, label=sample)
            local_max = max(local_max, float(np.max(c.values)) if len(c.values) else 0.0)
        data = curves[base_name]["Data top-10"]
        ax.errorbar(
            data.centers,
            data.values,
            yerr=data.errors,
            fmt="o",
            markersize=2.7,
            color=colors["Data top-10"],
            elinewidth=0.7,
            capsize=0,
            label="Data top-10",
        )
        local_max = max(local_max, float(np.max(data.values)) if len(data.values) else 0.0)
        ax.set_title(title, fontsize=14.5, fontweight="bold")
        ax.set_ylim(0.0, local_max * 1.12 if local_max > 0 else 1.0)
        ax.grid(True, color="#D9DEE7", linewidth=0.7, alpha=0.8)
        ax.tick_params(direction="in", top=True, right=True, labelsize=10.5)
        ax.set_ylabel("unit-normalized counts", fontsize=12)
    axes.flat[0].legend(loc="upper right", frameon=False, fontsize=10.5)
    for ax in axes[1]:
        ax.set_xlabel("isolation energy", fontsize=12)
    fig.suptitle(
        "THE-42 top-10 output isolation distributions, linear scale: 22-28 GeV, 0-80% centrality",
        fontsize=18,
        fontweight="bold",
        y=0.985,
    )
    fig.text(
        0.015,
        0.012,
        "Merged top-10 AuAu data output compared to existing merged signal and inclusive MC; shapes are normalized in the displayed x-range.",
        fontsize=10.5,
        color="#374151",
    )
    fig.tight_layout(rect=(0, 0.035, 1, 0.955))
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / "top10_isolation_distributions_linear.png"
    fig.savefig(path)
    plt.close(fig)
    return path


def write_manifest(curves: dict[str, dict[str, Curve]], plot_path: Path, args: argparse.Namespace) -> Path:
    payload = {
        "schema": "THE42_TOP10_ISOLATION_LINEAR_DIAGNOSTIC_V1",
        "plot": str(plot_path.resolve()),
        "data_root": str(Path(args.data_root).resolve()),
        "signal_root": str(Path(args.signal_root).resolve()),
        "inclusive_root": str(Path(args.inclusive_root).resolve()),
        "pt_bins": PT_BINS,
        "centrality_bins": CENT_BINS,
        "x_range": [args.x_min, args.x_max],
        "y_scale": "linear",
        "curves": {
            base_name: {
                sample: {
                    "entries": curve.entries,
                    "displayed_integral": curve.raw_integral,
                    "matches": curve.matches,
                }
                for sample, curve in sample_curves.items()
            }
            for base_name, sample_curves in curves.items()
        },
    }
    manifest_path = plot_path.with_suffix(".manifest.json")
    manifest_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return manifest_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", default=str(DEFAULT_DATA))
    parser.add_argument("--signal-root", default=str(DEFAULT_SIGNAL))
    parser.add_argument("--inclusive-root", default=str(DEFAULT_INCLUSIVE))
    parser.add_argument("--outdir", default=str(DEFAULT_OUTDIR))
    parser.add_argument("--trigger-regex", default=r"MBD_NS_geq_2_vtx_lt_150")
    parser.add_argument("--x-min", type=float, default=-20.0)
    parser.add_argument("--x-max", type=float, default=52.0)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    setup_style()
    curves = load_curves(args)
    plot_path = render(curves, Path(args.outdir))
    manifest_path = write_manifest(curves, plot_path, args)
    print(f"Wrote {plot_path}")
    print(f"Wrote {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
