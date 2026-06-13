#!/usr/bin/env python3
"""0-20% AuAu preselection/tight shower-shape overlay slide."""

from __future__ import annotations

import json
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle


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

AA_PT_BINS = [(15, 18), (18, 20), (20, 22), (22, 24), (24, 26), (26, 28), (28, 30), (30, 35)]
CENT_BINS_0_20 = [(0, 10), (10, 20)]

VARIABLES = [
    ("weta", r"$w_{\eta}^{\mathrm{cogX}}$", "weta_cogx", (0.0, 0.35), 2),
    ("wphi", r"$w_{\phi}^{\mathrm{cogX}}$", "wphi_cogx", (0.0, 0.35), 2),
    ("e11e33", r"$E_{11}/E_{33}$", "e11e33", (0.0, 1.08), 2),
]
STAGES = [
    ("After preselection", "pre"),
    ("After tight ID", "tight"),
]


@dataclass
class Curve:
    sample: str
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
            "figure.facecolor": "#F4F6F8",
            "savefig.facecolor": "#F4F6F8",
            "axes.edgecolor": "#1F2937",
            "axes.linewidth": 0.80,
            "xtick.color": "#111827",
            "ytick.color": "#111827",
            "xtick.labelsize": 9.0,
            "ytick.labelsize": 9.0,
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


def choose_path(paths: list[str], required_prefix: str | None) -> str | None:
    if not paths:
        return None
    if required_prefix is None:
        return paths[0]
    for path in paths:
        if path.startswith(required_prefix):
            return path
    return None


def aa_target_name(var: str, tag: str, pt: tuple[int, int], cent: tuple[int, int]) -> str:
    return f"h_ss_{var}_{tag}_pT_{pt[0]}_{pt[1]}_cent_{cent[0]}_{cent[1]}"


def rebin(hist, factor: int):
    if factor <= 1:
        return hist
    if hist.GetNbinsX() % factor != 0:
        return hist
    out = hist.Rebin(factor, f"{hist.GetName()}_rebin{factor}")
    out.SetDirectory(0)
    return out


def hist_to_curve(hist, sample: str, x_range: tuple[float, float], rebin_factor: int, matches: list[str], missing: list[str]) -> Curve:
    hist = rebin(hist, rebin_factor)
    nb = hist.GetNbinsX()
    edges = np.array([hist.GetXaxis().GetBinLowEdge(i) for i in range(1, nb + 2)], dtype=float)
    centers = np.array([hist.GetXaxis().GetBinCenter(i) for i in range(1, nb + 1)], dtype=float)
    counts = np.array([hist.GetBinContent(i) for i in range(1, nb + 1)], dtype=float)
    errors = np.array([hist.GetBinError(i) for i in range(1, nb + 1)], dtype=float)
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
    return Curve(sample, centers, edges, values, scaled_errors, integral, matches, missing)


def load_aa_curve(
    root_path: Path,
    index: dict[str, list[str]],
    sample: str,
    tag: str,
    var: str,
    x_range: tuple[float, float],
    rebin_factor: int,
) -> Curve:
    f = open_root(root_path)
    acc = None
    matches: list[str] = []
    missing: list[str] = []
    try:
        for pt in AA_PT_BINS:
            for cent in CENT_BINS_0_20:
                name = aa_target_name(var, tag, pt, cent)
                path = choose_path(index.get(name, []), None)
                if path is None:
                    missing.append(name)
                    continue
                obj = f.Get(path)
                if not obj or not obj.InheritsFrom("TH1"):
                    missing.append(path)
                    continue
                h = obj.Clone(f"{name}_{sample}")
                h.SetDirectory(0)
                matches.append(path)
                if acc is None:
                    acc = h.Clone(f"{var}_{sample}_{tag}_sum")
                    acc.SetDirectory(0)
                else:
                    acc.Add(h)
    finally:
        f.Close()
    if acc is None:
        raise RuntimeError(f"no AuAu histograms for {sample} {tag} {var}")
    return hist_to_curve(acc, sample, x_range, rebin_factor, matches, missing)


def add_sphenix_internal(ax) -> None:
    ax.text(
        0.026,
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
        0.151,
        0.955,
        "Internal",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=10.8,
        fontstyle="normal",
        fontweight="normal",
    )


def draw_line(ax, curve: Curve, color: str, linestyle: str = "-") -> None:
    ax.stairs(curve.values, curve.edges, color=color, lw=2.55, linestyle=linestyle, alpha=0.98)


def draw_open_markers(ax, curve: Curve, color: str, marker: str, size: float, zorder: int) -> None:
    ax.errorbar(
        curve.centers,
        curve.values,
        yerr=curve.errors,
        fmt=marker,
        ms=size,
        color=color,
        mfc="white",
        mec=color,
        mew=1.65,
        elinewidth=1.10,
        capsize=2.0,
        capthick=1.0,
        alpha=0.96,
        zorder=zorder,
    )


def apply_column_backdrops(fig, axes) -> None:
    for col, color in enumerate(("#EEF2F7", "#ECF7F1")):
        xs = [axes[row, col].get_position().x0 for row in range(3)]
        xe = [axes[row, col].get_position().x1 for row in range(3)]
        ys = [axes[row, col].get_position().y0 for row in range(3)]
        ye = [axes[row, col].get_position().y1 for row in range(3)]
        pad_x = 0.012
        pad_y = 0.018
        rect = Rectangle(
            (min(xs) - pad_x, min(ys) - pad_y),
            max(xe) - min(xs) + 2 * pad_x,
            max(ye) - min(ys) + 2 * pad_y,
            transform=fig.transFigure,
            facecolor=color,
            edgecolor="#D1D5DB",
            linewidth=0.8,
            zorder=-2,
        )
        fig.patches.append(rect)


def render() -> tuple[Path, Path, Path]:
    setup_style()
    OUTDIR.mkdir(parents=True, exist_ok=True)

    indexes = {
        "data": build_index(DATA_ROOT),
        "signal": build_index(SIGNAL_ROOT),
        "inclusive": build_index(INCLUSIVE_ROOT),
    }

    fig, axes = plt.subplots(3, 2, figsize=slide_figsize(), constrained_layout=False)
    fig.subplots_adjust(left=0.075, right=0.985, top=0.770, bottom=0.095, wspace=0.112, hspace=0.300)
    apply_column_backdrops(fig, axes)

    manifest = {
        "output": None,
        "data_root": str(DATA_ROOT),
        "signal_root": str(SIGNAL_ROOT),
        "inclusive_root": str(INCLUSIVE_ROOT),
        "auau_pt_bins": AA_PT_BINS,
        "centrality": "0-20%",
        "cfg_tag": CFG,
        "curves": [],
    }

    colors = {
        "signal": "#B91C1C",
        "inclusive": "#2E63D4",
        "auau": "#111827",
    }

    for row, (var, var_label, _slug, x_range, rebin_factor) in enumerate(VARIABLES):
        for col, (stage_label, stage_key) in enumerate(STAGES):
            ax = axes[row, col]
            ax.set_facecolor("#FFFFFF")
            ax.grid(True, axis="y", color="#E5E7EB", linewidth=0.75, alpha=0.95)
            ax.grid(True, axis="x", color="#EEF2F7", linewidth=0.55, alpha=0.70)

            curves = {
                "signal": load_aa_curve(SIGNAL_ROOT, indexes["signal"], "signal MC", f"{stage_key}_sig", var, x_range, rebin_factor),
                "inclusive": load_aa_curve(
                    INCLUSIVE_ROOT,
                    indexes["inclusive"],
                    "inclusive MC",
                    f"{stage_key}_bkg",
                    var,
                    x_range,
                    rebin_factor,
                ),
                "auau": load_aa_curve(DATA_ROOT, indexes["data"], "AuAu data", stage_key, var, x_range, rebin_factor),
            }

            draw_line(ax, curves["signal"], colors["signal"])
            draw_line(ax, curves["inclusive"], colors["inclusive"])
            draw_open_markers(ax, curves["auau"], colors["auau"], "o", 6.6, 6)

            max_y = max(float(np.nanmax(c.values)) if c.values.size else 0.0 for c in curves.values())
            ax.set_xlim(*x_range)
            ax.set_ylim(0.0, max_y * 1.28 if max_y > 0 else 1.0)
            ax.tick_params(direction="in", top=True, right=True, length=4.0, width=0.75)
            ax.set_xlabel(var_label, fontsize=10.8, labelpad=2.0)
            if col == 0:
                ax.set_ylabel("Unit-normalized counts", fontsize=10.3)
            else:
                ax.set_ylabel("")
            if row == 0:
                ax.set_title(stage_label, fontsize=15.2, fontweight="bold", pad=10)
            ax.text(
                0.985,
                0.900,
                f"{var_label}",
                transform=ax.transAxes,
                ha="right",
                va="top",
                fontsize=11.1,
                fontweight="bold",
                bbox={"boxstyle": "round,pad=0.22", "facecolor": "#F9FAFB", "edgecolor": "#E5E7EB", "linewidth": 0.6},
            )
            if row == 0 and col == 0:
                add_sphenix_internal(ax)

            for key, curve in curves.items():
                manifest["curves"].append(
                    {
                        "variable": var,
                        "stage": stage_key,
                        "sample": key,
                        "integral_before_normalization": curve.integral,
                        "matched_histograms": len(curve.matches),
                        "missing_histograms": curve.missing[:20],
                        "missing_histogram_count": len(curve.missing),
                    }
                )

    handles = [
        Line2D([0], [0], color=colors["signal"], lw=2.7, label="Signal MC"),
        Line2D([0], [0], color=colors["inclusive"], lw=2.7, label="Inclusive MC"),
        Line2D(
            [0],
            [0],
            marker="o",
            color=colors["auau"],
            mfc="white",
            mec=colors["auau"],
            mew=1.7,
            lw=0,
            ms=7.0,
            label="AuAu data",
        ),
    ]
    fig.legend(
        handles=handles,
        loc="upper right",
        bbox_to_anchor=(0.985, 0.905),
        ncol=3,
        frameon=True,
        fancybox=False,
        framealpha=0.92,
        edgecolor="#D1D5DB",
        facecolor="#F9FAFB",
        fontsize=11.7,
        handlelength=2.2,
        columnspacing=1.2,
        borderpad=0.55,
    )
    fig.text(
        0.075,
        0.955,
        "0-20% AuAu shower-shape overlays before and after tight photon ID",
        ha="left",
        va="top",
        fontsize=25.5,
        fontweight="bold",
        color="#111827",
    )
    fig.text(
        0.075,
        0.910,
        r"0-20% AuAu, photon $p_T$: 15-35 GeV; unit-normalized shapes.",
        ha="left",
        va="top",
        fontsize=13.5,
        color="#374151",
    )
    fig.text(
        0.075,
        0.875,
        "Latest-production data compared with embedded signal and inclusive MC.",
        ha="left",
        va="top",
        fontsize=12.2,
        color="#4B5563",
    )

    out = OUTDIR / "the42_b009_0_20_auau_only_top3_shapes_slide.png"
    manifest_path = OUTDIR / "the42_b009_0_20_auau_only_top3_shapes_manifest.json"
    script_path = OUTDIR / "the42_b009_0_20_auau_only_top3_shapes_speaker_script.md"
    fig.savefig(out, dpi=SLIDE_DPI)
    plt.close(fig)
    manifest["output"] = str(out)
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    script_path.write_text(
        "\n".join(
            [
                "# Speaker Notes",
                "",
                "This slide fixes the comparison to the 0-20% AuAu bin and shows three shower-shape variables in the same two-column structure: preselection on the left and tight photon ID on the right.",
                "",
                "The lines are embedded AuAu signal and inclusive MC. The open black circles are the latest-production AuAu data. Each panel is unit-normalized so the comparison is about shape rather than sample yield.",
                "",
                "The intended readout is whether the tight-ID column suppresses inclusive-like structure and pulls the AuAu data toward the signal-like shape without reproducing the old edge-spike pathology.",
                "",
            ]
        )
    )
    return out, manifest_path, script_path


def main() -> int:
    out, manifest, script = render()
    print(out)
    print(manifest)
    print(script)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
