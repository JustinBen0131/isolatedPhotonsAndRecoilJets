#!/usr/bin/env python3
"""Clean one-variable b009 AuAu preselection vs tight-BDT shower-shape slide."""

from __future__ import annotations

import json
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


THIS_FILE = Path(__file__).resolve()
REPO = next((p for p in THIS_FILE.parents if (p / "AGENTS.md").exists()), THIS_FILE.parents[5])
SCRIPTS_DIR = REPO / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.append(str(SCRIPTS_DIR))

from slides.common.slide_defaults import SLIDE_DPI, slide_figsize  # noqa: E402


BASE = REPO / "dataOutput/auauTightBDTValidation/THE42_wp80_centlinear_ss_20260604"
CFG = "preselectionNewPPG12_tightAuAuCentInputBase3x3BDT_nonTightAuAuBDTComplement_baseVariant"
DATA_ROOT = BASE / "b009_merged_data_roots_20260609" / f"RecoilJets_auau_ALL_{CFG}.root"
SIGNAL_ROOT = BASE / "sim_roots" / CFG / "photonJet12and20merged_SIM" / "RecoilJets_embeddedPhoton12plus20_MERGED.root"
INCLUSIVE_ROOT = (
    BASE
    / "sim_roots"
    / CFG
    / "embeddedJet12and20and30and40merged_SIM"
    / "RecoilJets_embeddedJet12plus20plus30plus40_MERGED.root"
)
OUTDIR = BASE / "pre_vs_tight_shower_shape_b009_20260610"

PT_BINS = [(15, 18), (18, 20), (20, 22), (22, 24), (24, 26), (26, 28), (28, 30), (30, 35)]
CENT_GROUPS = [
    ("0-20%", [(0, 10), (10, 20)]),
    ("20-50%", [(20, 30), (30, 40), (40, 50)]),
    ("50-80%", [(50, 60), (60, 80)]),
]
VARIABLES = [
    ("weta", r"$w_{\eta}^{\mathrm{cogX}}$", "weta_cogx", (0.0, 0.35), 2),
    ("e11e33", r"$E_{11}/E_{33}$", "e11e33", (0.0, 1.08), 2),
]
STAGES = [
    ("After preselection", "pre"),
    ("After tight ID", "tight"),
]
SAMPLES = [
    ("Data", DATA_ROOT, {"pre": "pre", "tight": "tight"}, "#111827"),
    ("Signal MC", SIGNAL_ROOT, {"pre": "pre_sig", "tight": "tight_sig"}, "#C22F2F"),
    ("Inclusive MC", INCLUSIVE_ROOT, {"pre": "pre_bkg", "tight": "tight_bkg"}, "#2E63D4"),
]


@dataclass
class Curve:
    sample: str
    stage: str
    centers: np.ndarray
    edges: np.ndarray
    values: np.ndarray
    errors: np.ndarray
    integral: float
    matches: list[str]
    missing: list[str]


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "axes.edgecolor": "#111827",
            "axes.linewidth": 0.85,
            "xtick.color": "#111827",
            "ytick.color": "#111827",
        }
    )


def open_root(path: Path):
    import ROOT

    ROOT.gROOT.SetBatch(True)
    ROOT.TH1.AddDirectory(False)
    f = ROOT.TFile.Open(str(path), "READ")
    if not f or f.IsZombie():
        raise OSError(f"could not open {path}")
    return f


def walk(directory, prefix: str = ""):
    for key in directory.GetListOfKeys():
        obj = key.ReadObj()
        path = f"{prefix}/{key.GetName()}" if prefix else key.GetName()
        yield path, obj
        if obj.InheritsFrom("TDirectory"):
            yield from walk(obj, path)


def build_index(path: Path) -> dict[str, list[str]]:
    f = open_root(path)
    idx: dict[str, list[str]] = {}
    try:
        for obj_path, obj in walk(f):
            if obj.InheritsFrom("TH1"):
                idx.setdefault(Path(obj_path).name, []).append(obj_path)
    finally:
        f.Close()
    return idx


def target_name(var: str, tag: str, pt: tuple[int, int], cent: tuple[int, int]) -> str:
    return f"h_ss_{var}_{tag}_pT_{pt[0]}_{pt[1]}_cent_{cent[0]}_{cent[1]}"


def rebin(hist, factor: int):
    if factor <= 1:
        return hist
    if hist.GetNbinsX() % factor != 0:
        return hist
    out = hist.Rebin(factor, f"{hist.GetName()}_rebin{factor}")
    out.SetDirectory(0)
    return out


def load_curve(
    root_path: Path,
    index: dict[str, list[str]],
    sample: str,
    stage: str,
    tag: str,
    var: str,
    cent_bins: list[tuple[int, int]],
    x_range: tuple[float, float],
    rebin_factor: int,
) -> Curve:
    f = open_root(root_path)
    acc = None
    matches: list[str] = []
    missing: list[str] = []
    try:
        for pt in PT_BINS:
            for cent in cent_bins:
                name = target_name(var, tag, pt, cent)
                paths = index.get(name, [])
                if not paths:
                    missing.append(name)
                    continue
                obj = f.Get(paths[0])
                if not obj or not obj.InheritsFrom("TH1"):
                    missing.append(paths[0])
                    continue
                h = obj.Clone(f"{name}_{sample}_{stage}")
                h.SetDirectory(0)
                matches.append(paths[0])
                if acc is None:
                    acc = h.Clone(f"{var}_{sample}_{stage}_sum")
                    acc.SetDirectory(0)
                else:
                    acc.Add(h)
    finally:
        f.Close()
    if acc is None:
        raise RuntimeError(f"no histograms for {sample} {stage} {var}")
    acc = rebin(acc, rebin_factor)
    nb = acc.GetNbinsX()
    edges = np.array([acc.GetXaxis().GetBinLowEdge(i) for i in range(1, nb + 2)], dtype=float)
    centers = np.array([acc.GetXaxis().GetBinCenter(i) for i in range(1, nb + 1)], dtype=float)
    counts = np.array([acc.GetBinContent(i) for i in range(1, nb + 1)], dtype=float)
    errors = np.array([acc.GetBinError(i) for i in range(1, nb + 1)], dtype=float)
    mask = (centers >= x_range[0]) & (centers < x_range[1])
    keep = np.flatnonzero(mask)
    if len(keep):
        edges = np.concatenate(([edges[keep[0]]], edges[keep + 1]))
    centers = centers[mask]
    counts = counts[mask]
    errors = errors[mask]
    integral = float(np.sum(counts))
    values = np.zeros_like(counts)
    scaled_errors = np.zeros_like(errors)
    if integral > 0:
        values = counts / integral
        scaled_errors = errors / integral
    return Curve(sample, stage, centers, edges, values, scaled_errors, integral, matches, missing)


def draw_curve(ax, curve: Curve, color: str) -> None:
    if curve.sample != "Data":
        ax.step(curve.edges[:-1], curve.values, where="post", color=color, lw=2.8, alpha=0.96)
        return
    # Data is rendered as markers only so the MC line shapes remain easy to read.
    ax.errorbar(
        curve.centers,
        curve.values,
        yerr=curve.errors,
        fmt="o",
        ms=5.0,
        color=color,
        mfc=color,
        mec="white",
        mew=0.65,
        elinewidth=1.15,
        capsize=2.0,
        capthick=1.05,
        alpha=0.95,
    )


def add_sphenix_internal(ax) -> None:
    ax.text(
        0.025,
        0.955,
        "sPHENIX",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=10.8,
        fontstyle="italic",
        fontweight="bold",
    )
    ax.text(
        0.150,
        0.955,
        "Internal",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=10.8,
        fontstyle="normal",
        fontweight="normal",
    )


def render(variable: tuple[str, str, str, tuple[float, float], int]) -> tuple[Path, Path, Path]:
    setup_style()
    OUTDIR.mkdir(parents=True, exist_ok=True)
    indexes = {root: build_index(root) for _, root, _, _ in SAMPLES}

    fig, axes = plt.subplots(3, 2, figsize=slide_figsize(), constrained_layout=False)
    fig.subplots_adjust(left=0.080, right=0.985, top=0.762, bottom=0.090, wspace=0.120, hspace=0.285)
    manifest = {
        "data_root": str(DATA_ROOT),
        "signal_root": str(SIGNAL_ROOT),
        "inclusive_root": str(INCLUSIVE_ROOT),
        "pt_bins": PT_BINS,
        "centrality_groups": {label: bins for label, bins in CENT_GROUPS},
        "variable": variable[0],
        "curves": [],
    }
    var, var_label, slug, x_range, rebin_factor = variable

    for row, (cent_label, cent_bins) in enumerate(CENT_GROUPS):
        for col, (stage_label, stage_key) in enumerate(STAGES):
            ax = axes[row, col]
            ymax = 0.0
            for sample, root, tags, color in SAMPLES:
                curve = load_curve(
                    root,
                    indexes[root],
                    sample,
                    stage_label,
                    tags[stage_key],
                    var,
                    cent_bins,
                    x_range,
                    rebin_factor,
                )
                draw_curve(ax, curve, color)
                ymax = max(ymax, float(np.max(curve.values + curve.errors)) if len(curve.values) else 0.0)
                manifest["curves"].append(
                    {
                        "centrality": cent_label,
                        "variable": var,
                        "sample": sample,
                        "stage": stage_label,
                        "tag": tags[stage_key],
                        "integral": curve.integral,
                        "n_matches": len(curve.matches),
                        "n_missing": len(curve.missing),
                    }
                )
            ax.set_xlim(*x_range)
            ax.set_ylim(0.0, max(0.06, ymax * 1.18))
            ax.grid(True, color="#E5E7EB", lw=0.55, alpha=0.75)
            ax.tick_params(labelsize=9.8, pad=1)
            if row == 2:
                ax.set_xlabel(var_label, fontsize=13.2)
            if col == 0:
                ax.set_ylabel(f"{cent_label}\nunit-normalized", fontsize=11.2)
            else:
                ax.set_ylabel("")
            if row == 0:
                ax.set_title(stage_label, fontsize=17.0, fontweight="bold", pad=7)
                if col == 0:
                    add_sphenix_internal(ax)
            data_row = next(c for c in manifest["curves"] if c["centrality"] == cent_label and c["sample"] == "Data" and c["stage"] == stage_label)
            ax.text(
                0.98,
                0.91,
                f"data n={data_row['integral']:.0f}",
                transform=ax.transAxes,
                ha="right",
                va="top",
                fontsize=8.2,
                color="#374151",
            )

    fig.text(0.080, 0.958, rf"Tight ID pulls AuAu data toward signal-like {var_label}", fontsize=23, fontweight="bold", ha="left", va="top")
    fig.text(
        0.080,
        0.908,
        r"b009 latest-production v008, $15 \leq E_T^\gamma < 35$ GeV; MC shown as histogram lines, data as black markers",
        fontsize=14.0,
        ha="left",
        va="top",
    )

    legend_items = [
        plt.Line2D([0], [0], color="#C22F2F", lw=2.8, label="Signal MC"),
        plt.Line2D([0], [0], color="#2E63D4", lw=2.8, label="Inclusive MC"),
        plt.Line2D([0], [0], color="#111827", marker="o", markerfacecolor="#111827", markeredgecolor="white", markeredgewidth=0.7, markersize=8.0, lw=0, label="Data"),
    ]
    fig.legend(handles=legend_items, loc="upper right", bbox_to_anchor=(0.985, 0.850), ncol=3, frameon=False, fontsize=12.5)

    out_png = OUTDIR / f"the42_b009_pre_vs_tight_{slug}_summary_slide.png"
    out_json = OUTDIR / f"the42_b009_pre_vs_tight_{slug}_summary_manifest.json"
    out_script = OUTDIR / f"the42_b009_pre_vs_tight_{slug}_summary_speaker_script.md"
    fig.savefig(out_png, dpi=SLIDE_DPI)
    plt.close(fig)
    out_json.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    out_script.write_text(
        "\n".join(
            [
                "# Speaker Script",
                "",
                "This slide is a direct before-and-after check of the AuAu tight photon ID.",
                f"Each row is a centrality range, and the two columns show {var_label} after the analysis preselection and after the tight WP80 BDT cut.",
                "The red and blue histograms are the embedded signal and inclusive-jet MC shapes, while the black markers are the latest-production b009 AuAu data.",
                "",
                "The key point is the right column. After the tight ID, the data markers move much closer to the signal-like shape, especially in the central and mid-central bins.",
                "The inclusive MC remains the broader comparison shape, so the slide shows the BDT selection doing what it is supposed to do in data rather than only in simulation.",
                "",
                "This is still an internal validation plot, not a final corrected physics observable. It is a shape-closure check for the trained AuAu BDT selection.",
                "",
            ]
        )
    )
    return out_png, out_json, out_script


def main() -> int:
    outputs = []
    for variable in VARIABLES:
        out_png, out_json, out_script = render(variable)
        outputs.append({"variable": variable[0], "png": str(out_png), "manifest": str(out_json), "speaker_script": str(out_script)})
    print(json.dumps({"outputs": outputs}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
