#!/usr/bin/env python3
"""Make reco-level xJgamma overlays for THE-42 b009 AuAu data vs pp."""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


THIS_FILE = Path(__file__).resolve()
REPO = next((p for p in THIS_FILE.parents if (p / "AGENTS.md").exists()), THIS_FILE.parents[3])
BASE = REPO / "dataOutput/auauTightBDTValidation/THE42_wp80_centlinear_ss_20260604"
CFG = "preselectionNewPPG12_tightAuAuCentInputBase3x3BDT_nonTightAuAuBDTComplement_baseVariant"
DEFAULT_AUAU = BASE / "b009_merged_data_roots_20260609" / f"RecoilJets_auau_ALL_{CFG}.root"
DEFAULT_PP = (
    REPO
    / "InputFiles/pp24/RecoilJets_pp_ALL_jetMinPt5_7pi_8_vz60_isoR40_fixedIso2GeV_preselectionReference_tightReference_nonTightReference.root"
)
DEFAULT_OUTDIR = BASE / "reco_xj_b009_20260610"

PT_RANGE = (15.0, 35.0)
CENT_GROUPS = {
    "0-20%": [(0, 10), (10, 20)],
    "20-40%": [(20, 30), (30, 40)],
    "40-80%": [(40, 50), (50, 60), (60, 80)],
}
CENT_COLORS = {
    "0-20%": "#D55E00",
    "20-40%": "#009E73",
    "40-80%": "#0072B2",
}


@dataclass
class Curve:
    label: str
    centers: np.ndarray
    values: np.ndarray
    errors: np.ndarray
    raw_integral: float
    raw_entries: float
    source_paths: list[str]
    missing_paths: list[str]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--auau", type=Path, default=DEFAULT_AUAU)
    parser.add_argument("--pp", type=Path, default=DEFAULT_PP)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--pt-min", type=float, default=PT_RANGE[0])
    parser.add_argument("--pt-max", type=float, default=PT_RANGE[1])
    parser.add_argument("--xj-min", type=float, default=0.0)
    parser.add_argument("--xj-max", type=float, default=2.0)
    return parser.parse_args()


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "axes.edgecolor": "#111827",
            "axes.linewidth": 1.1,
            "axes.labelcolor": "#111827",
            "xtick.color": "#111827",
            "ytick.color": "#111827",
            "legend.frameon": False,
        }
    )


def open_root(path: Path):
    import ROOT

    ROOT.gROOT.SetBatch(True)
    ROOT.TH1.AddDirectory(False)
    f = ROOT.TFile.Open(str(path), "READ")
    if not f or f.IsZombie():
        raise OSError(f"could not open ROOT file: {path}")
    return f


def walk(directory, prefix: str = ""):
    for key in directory.GetListOfKeys():
        obj = key.ReadObj()
        path = f"{prefix}/{key.GetName()}" if prefix else key.GetName()
        yield path, obj
        if obj.InheritsFrom("TDirectory"):
            yield from walk(obj, path)


def make_index(root_path: Path) -> dict[str, list[str]]:
    root_file = open_root(root_path)
    index: dict[str, list[str]] = {}
    try:
        for path, obj in walk(root_file):
            if obj.InheritsFrom("TH2"):
                index.setdefault(Path(path).name, []).append(path)
    finally:
        root_file.Close()
    return index


def clone_hist2(root_file, path: str, name: str):
    obj = root_file.Get(path)
    if not obj or not obj.InheritsFrom("TH2"):
        return None
    out = obj.Clone(name)
    out.SetDirectory(0)
    return out


def project_y(hist2, name: str, pt_min: float, pt_max: float):
    ix_lo = hist2.GetXaxis().FindBin(pt_min + 1e-6)
    ix_hi = hist2.GetXaxis().FindBin(pt_max - 1e-6)
    hist = hist2.ProjectionY(name, ix_lo, ix_hi, "e")
    hist.SetDirectory(0)
    return hist


def normalize_curve(hist, label: str, source_paths: list[str], missing_paths: list[str], xj_min: float, xj_max: float) -> Curve:
    nb = hist.GetNbinsX()
    centers = np.array([hist.GetXaxis().GetBinCenter(i) for i in range(1, nb + 1)], dtype=float)
    counts = np.array([hist.GetBinContent(i) for i in range(1, nb + 1)], dtype=float)
    errors = np.array([hist.GetBinError(i) for i in range(1, nb + 1)], dtype=float)
    mask = (centers >= xj_min) & (centers <= xj_max)
    centers = centers[mask]
    counts = counts[mask]
    errors = errors[mask]
    raw_integral = float(np.sum(counts))
    values = np.zeros_like(counts)
    scaled_errors = np.zeros_like(errors)
    if raw_integral > 0:
        values = counts / raw_integral
        scaled_errors = errors / raw_integral
    return Curve(
        label=label,
        centers=centers,
        values=values,
        errors=scaled_errors,
        raw_integral=raw_integral,
        raw_entries=float(hist.GetEntries()),
        source_paths=source_paths,
        missing_paths=missing_paths,
    )


def sum_auau_cent_group(
    root_path: Path,
    index: dict[str, list[str]],
    name_template: str,
    label: str,
    cent_bins: list[tuple[int, int]],
    pt_min: float,
    pt_max: float,
    xj_min: float,
    xj_max: float,
) -> Curve:
    root_file = open_root(root_path)
    acc = None
    source_paths: list[str] = []
    missing_paths: list[str] = []
    try:
        for cent_lo, cent_hi in cent_bins:
            hname = name_template.format(cent_lo=cent_lo, cent_hi=cent_hi)
            matches = index.get(hname, [])
            if not matches:
                missing_paths.append(hname)
                continue
            # Prefer the first exact basename match. The merged files use unique basenames.
            hist2 = clone_hist2(root_file, matches[0], f"{hname}_clone")
            if hist2 is None:
                missing_paths.append(matches[0])
                continue
            proj = project_y(hist2, f"{hname}_xj_{pt_min:g}_{pt_max:g}", pt_min, pt_max)
            source_paths.append(matches[0])
            if acc is None:
                acc = proj.Clone(f"{label}_sum")
                acc.SetDirectory(0)
            else:
                acc.Add(proj)
    finally:
        root_file.Close()
    if acc is None:
        raise RuntimeError(f"no AuAu histograms found for {label} using {name_template}")
    return normalize_curve(acc, label, source_paths, missing_paths, xj_min, xj_max)


def load_pp_curve(
    root_path: Path,
    index: dict[str, list[str]],
    pt_min: float,
    pt_max: float,
    xj_min: float,
    xj_max: float,
) -> Curve:
    target = "h2_unfoldReco_pTgamma_xJ_incl_r03"
    matches = index.get(target, [])
    if not matches:
        raise RuntimeError(f"missing pp reference histogram {target}")
    preferred = next((p for p in matches if p.startswith("Photon_4_GeV_plus_MBD_NS_geq_1/")), matches[0])
    root_file = open_root(root_path)
    try:
        hist2 = clone_hist2(root_file, preferred, "pp_h2_clone")
        if hist2 is None:
            raise RuntimeError(f"could not clone {preferred}")
        proj = project_y(hist2, f"pp_xj_{pt_min:g}_{pt_max:g}", pt_min, pt_max)
    finally:
        root_file.Close()
    return normalize_curve(proj, "pp Run24 reference", [preferred], [], xj_min, xj_max)


def plot_curves(curves: list[Curve], out_png: Path, title: str, subtitle: str) -> None:
    fig, ax = plt.subplots(figsize=(9.6, 6.2))
    fig.subplots_adjust(top=0.84, left=0.10, right=0.985, bottom=0.13)
    pp = curves[0]
    ax.step(pp.centers, pp.values, where="mid", color="#111827", lw=2.4, label=f"{pp.label} ({pp.raw_integral:.0f})")
    ax.errorbar(pp.centers, pp.values, yerr=pp.errors, fmt="o", ms=3.1, color="#111827", lw=1.2)

    for curve in curves[1:]:
        color = CENT_COLORS.get(curve.label, "#666666")
        ax.step(curve.centers, curve.values, where="mid", color=color, lw=2.2, label=f"AuAu {curve.label} ({curve.raw_integral:.0f})")
        ax.errorbar(curve.centers, curve.values, yerr=curve.errors, fmt="o", ms=4.0, color=color, lw=1.2, capsize=1.8)

    ax.set_xlim(0.0, 2.0)
    ymax = max((float(np.nanmax(c.values + c.errors)) if len(c.values) else 0.0) for c in curves)
    ax.set_ylim(0.0, max(0.055, ymax * 1.28))
    ax.set_xlabel(r"$x_{J\gamma} = p_T^{jet}/p_T^\gamma$")
    ax.set_ylabel("Unit-normalized candidates")
    fig.text(0.10, 0.955, title, ha="left", va="top", fontsize=17, fontweight="bold")
    fig.text(0.10, 0.912, subtitle, ha="left", va="top", fontsize=11.5)
    fig.text(
        0.10,
        0.875,
        "Reco-level projection from unfolding-input TH2s; not a completed unfolded result.",
        ha="left",
        va="top",
        fontsize=10,
        color="#4B5563",
    )
    ax.legend(loc="upper right", fontsize=10)
    ax.grid(True, color="#D1D5DB", linewidth=0.8, alpha=0.75)
    ax.axvline(1.0, color="#6B7280", lw=1.0, ls=":")
    fig.savefig(out_png, dpi=220)
    plt.close(fig)


def make_overlay(
    auau_path: Path,
    pp_path: Path,
    outdir: Path,
    name_template: str,
    tag: str,
    pt_min: float,
    pt_max: float,
    xj_min: float,
    xj_max: float,
) -> dict:
    auau_index = make_index(auau_path)
    pp_index = make_index(pp_path)
    curves = [load_pp_curve(pp_path, pp_index, pt_min, pt_max, xj_min, xj_max)]
    for label, cent_bins in CENT_GROUPS.items():
        curves.append(sum_auau_cent_group(auau_path, auau_index, name_template, label, cent_bins, pt_min, pt_max, xj_min, xj_max))

    out_png = outdir / f"the42_b009_reco_xj_centrality_vs_pp_{tag}.png"
    title = rf"Reco $x_{{J\gamma}}$, $\gamma$ $E_T$ {pt_min:g}-{pt_max:g} GeV"
    subtitle = f"AuAu b009 latest-production v008 BDT WP80 vs pp24 reference, anti-kT R=0.3, {tag}"
    plot_curves(curves, out_png, title, subtitle)

    return {
        "tag": tag,
        "png": str(out_png),
        "name_template": name_template,
        "pt_range": [pt_min, pt_max],
        "curves": [
            {
                "label": c.label,
                "raw_integral_in_projected_range": c.raw_integral,
                "raw_entries_from_root_projection": c.raw_entries,
                "source_paths": c.source_paths,
                "missing_paths": c.missing_paths,
            }
            for c in curves
        ],
    }


def main() -> int:
    args = parse_args()
    setup_style()
    args.outdir.mkdir(parents=True, exist_ok=True)

    overlays = [
        (
            "isoR40_fixedIso4GeV",
            "h2_unfoldReco_pTgamma_xJ_incl_r03_isoR40_fixedIso4GeV_cent_{cent_lo}_{cent_hi}",
        ),
        (
            "isoR40_isSliding",
            "h2_unfoldReco_pTgamma_xJ_incl_r03_isoR40_isSliding_cent_{cent_lo}_{cent_hi}",
        ),
    ]
    manifest = {
        "auau_root": str(args.auau),
        "pp_root": str(args.pp),
        "outdir": str(args.outdir),
        "note": "Reco-level xJ overlays from h2_unfoldReco input histograms. No completed validated b009 RooUnfold output product was found locally in the THE-42 output tree.",
        "unfolded_b009_product_found": False,
        "overlays": [],
    }
    for tag, template in overlays:
        manifest["overlays"].append(
            make_overlay(args.auau, args.pp, args.outdir, template, tag, args.pt_min, args.pt_max, args.xj_min, args.xj_max)
        )

    out_json = args.outdir / "the42_b009_reco_xj_overlay_manifest.json"
    out_json.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"manifest": str(out_json), "pngs": [o["png"] for o in manifest["overlays"]]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
