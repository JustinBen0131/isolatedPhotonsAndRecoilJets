#!/usr/bin/env python3
"""Render THE-42 signed shape differences for energy-sum BDT inputs."""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np


THIS_FILE = Path(__file__).resolve()
REPO = next((p for p in THIS_FILE.parents if (p / "AGENTS.md").exists()), THIS_FILE.parents[3])
DEFAULT_BASE = REPO / "dataOutput/auauTightBDTValidation/THE42_wp80_centlinear_ss_20260604"
DEFAULT_ROOT_DIR = (
    DEFAULT_BASE
    / "sim_roots/preselectionNewPPG12_tightAuAuCentInputBase3x3BDT_nonTightAuAuBDTComplement_baseVariant"
)
DEFAULT_SIGNAL_ROOT = (
    DEFAULT_ROOT_DIR / "photonJet12and20merged_SIM/RecoilJets_embeddedPhoton12plus20_MERGED.root"
)
DEFAULT_INCLUSIVE_ROOT = (
    DEFAULT_ROOT_DIR / "embeddedJet12and20and30and40merged_SIM/RecoilJets_embeddedJet12plus20plus30plus40_MERGED.root"
)
DEFAULT_OUTDIR = DEFAULT_BASE / "efficiency_purity_qa"

ET_BINS = (
    (r"$22 \leq E_T < 24$ GeV", (22, 24)),
    (r"$24 \leq E_T < 26$ GeV", (24, 26)),
    (r"$26 \leq E_T < 28$ GeV", (26, 28)),
)
CENT_FOCUS_LABEL = "0-20%"
CENT_FOCUS_BINS = ((0, 10), (10, 20))
VARIABLES = (
    {
        "key": "et1",
        "feature": "cluster_et1",
        "label": r"$E_1/E_{\mathrm{cluster}}$",
        "plain": "cluster_et1",
        "xlim": (0.25, 1.02),
    },
    {
        "key": "e11e33",
        "feature": "e11_over_e33",
        "label": r"$E_{1\times1}/E_{3\times3}$",
        "plain": "E11/E33",
        "xlim": (0.00, 0.96),
    },
    {
        "key": "e32e35",
        "feature": "e32_over_e35",
        "label": r"$E_{3\times2}/E_{3\times5}$",
        "plain": "E32/E35",
        "xlim": (0.45, 1.02),
    },
)
SAMPLE_TAGS = {
    "Signal MC": {"before": "inclusive_sig", "tight": "tight_sig"},
    "Inclusive MC": {"before": "inclusive_bkg", "tight": "tight_bkg"},
}

INK = "#111827"
MUTED = "#566274"
GRID = "#E4EBF3"
PANEL_EDGE = "#C8D6E5"
SOFT_PANEL = "#F8FAFC"
SIGNAL = "#B91C1C"
INCLUSIVE = "#174EA6"
INCLUSIVE_REMOVED = "#6EA7DC"
INCLUSIVE_ADDED = "#174EA6"
TAKEAWAY = "#0F766E"


@dataclass(frozen=True)
class ShapeDiff:
    sample: str
    variable: str
    feature: str
    et_bin: str
    centrality: str
    edges: list[float]
    delta_density: list[float]
    before_mean: float
    tight_mean: float
    delta_mean: float
    before_integral: float
    tight_integral: float
    positive_area: float
    negative_area: float


def setup_style() -> None:
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
        "mathtext.fontset": "dejavuserif",
        "axes.edgecolor": INK,
        "axes.linewidth": 0.9,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
    })


def open_root(path: Path):
    import ROOT

    ROOT.gROOT.SetBatch(True)
    ROOT.TH1.AddDirectory(False)
    f = ROOT.TFile.Open(str(path), "READ")
    if not f or f.IsZombie():
        raise OSError(f"Could not open ROOT input: {path}")
    return f


def sum_hist(root_file, variable_key: str, tag: str, et_bin: tuple[int, int]):
    acc = None
    missing = []
    pt_lo, pt_hi = et_bin
    for c_lo, c_hi in CENT_FOCUS_BINS:
        name = f"SIM/h_ss_{variable_key}_{tag}_pT_{pt_lo}_{pt_hi}_cent_{c_lo}_{c_hi}"
        hist = root_file.Get(name)
        if not hist:
            missing.append(name)
            continue
        if acc is None:
            acc = hist.Clone(f"sum_{variable_key}_{tag}_{pt_lo}_{pt_hi}_{c_lo}_{c_hi}")
            acc.SetDirectory(0)
        else:
            acc.Add(hist)
    if acc is None:
        raise RuntimeError(f"No h_ss_{variable_key} histograms found for tag={tag}; first missing={missing[:3]}")
    return acc


def rebin_for_display(edges: np.ndarray, counts: np.ndarray, factor: int) -> tuple[np.ndarray, np.ndarray]:
    if factor <= 1:
        return edges, counts
    n = (len(counts) // factor) * factor
    if n <= 0:
        return edges, counts
    rebinned_counts = counts[:n].reshape(-1, factor).sum(axis=1)
    rebinned_edges = edges[: n + 1 : factor]
    if len(rebinned_edges) != len(rebinned_counts) + 1:
        rebinned_edges = np.append(rebinned_edges, edges[n])
    if n < len(counts):
        rebinned_counts = np.append(rebinned_counts, counts[n:].sum())
        rebinned_edges = np.append(rebinned_edges, edges[-1])
    return rebinned_edges, rebinned_counts


def hist_arrays(hist, variable: dict[str, object], rebin: int) -> tuple[np.ndarray, np.ndarray, float, float]:
    nbins = hist.GetNbinsX()
    raw_edges = np.array([hist.GetXaxis().GetBinLowEdge(i) for i in range(1, nbins + 2)], dtype=float)
    raw_centers = np.array([hist.GetXaxis().GetBinCenter(i) for i in range(1, nbins + 1)], dtype=float)
    raw_counts = np.array([hist.GetBinContent(i) for i in range(1, nbins + 1)], dtype=float)
    total = float(np.sum(raw_counts))
    mean = float(np.sum(raw_centers * raw_counts) / total) if total > 0 else math.nan

    edges, counts = rebin_for_display(raw_edges, raw_counts, rebin)
    left, right = variable["xlim"]
    centers = 0.5 * (edges[:-1] + edges[1:])
    keep = (centers >= left) & (centers <= right)
    if np.any(keep):
        first = int(np.where(keep)[0][0])
        last = int(np.where(keep)[0][-1])
        edges = edges[first:last + 2]
        counts = counts[first:last + 1]
    visible = float(np.sum(counts))
    widths = np.diff(edges)
    density = counts / (visible * widths) if visible > 0 else np.zeros_like(counts)
    return edges, density, mean, total


def make_diff(
    sample: str,
    variable: dict[str, object],
    et_label: str,
    before_hist,
    tight_hist,
    rebin: int,
) -> ShapeDiff:
    before_edges, before_density, before_mean, before_integral = hist_arrays(before_hist, variable, rebin)
    tight_edges, tight_density, tight_mean, tight_integral = hist_arrays(tight_hist, variable, rebin)
    if len(before_edges) != len(tight_edges) or np.max(np.abs(before_edges - tight_edges)) > 1e-9:
        raise RuntimeError(f"Histogram binning mismatch for {sample} {variable['plain']} {et_label}")
    delta = tight_density - before_density
    widths = np.diff(before_edges)
    positive_area = float(np.sum(np.clip(delta, 0, None) * widths))
    negative_area = float(np.sum(np.clip(delta, None, 0) * widths))
    return ShapeDiff(
        sample=sample,
        variable=str(variable["plain"]),
        feature=str(variable["feature"]),
        et_bin=et_label,
        centrality=CENT_FOCUS_LABEL,
        edges=[float(x) for x in before_edges],
        delta_density=[float(x) for x in delta],
        before_mean=before_mean,
        tight_mean=tight_mean,
        delta_mean=tight_mean - before_mean,
        before_integral=before_integral,
        tight_integral=tight_integral,
        positive_area=positive_area,
        negative_area=negative_area,
    )


def collect_diffs(signal_root: Path, inclusive_root: Path, rebin: int) -> list[ShapeDiff]:
    diffs: list[ShapeDiff] = []
    root_paths = {"Signal MC": signal_root, "Inclusive MC": inclusive_root}
    for sample, path in root_paths.items():
        root_file = open_root(path)
        try:
            for variable in VARIABLES:
                for et_label, et_bin in ET_BINS:
                    before = sum_hist(root_file, str(variable["key"]), SAMPLE_TAGS[sample]["before"], et_bin)
                    tight = sum_hist(root_file, str(variable["key"]), SAMPLE_TAGS[sample]["tight"], et_bin)
                    diffs.append(make_diff(sample, variable, et_label, before, tight, rebin))
        finally:
            root_file.Close()
    return diffs


def diff_lookup(diffs: list[ShapeDiff]) -> dict[tuple[str, str, str], ShapeDiff]:
    return {(d.sample, d.variable, d.et_bin): d for d in diffs}


def rounded_box(fig, xywh, face, edge=PANEL_EDGE, radius=0.018, lw=1.0):
    ax = fig.add_axes(xywh)
    ax.axis("off")
    patch = FancyBboxPatch(
        (0, 0), 1, 1,
        boxstyle=f"round,pad=0.006,rounding_size={radius}",
        facecolor=face,
        edgecolor=edge,
        linewidth=lw,
        transform=ax.transAxes,
        clip_on=False,
    )
    ax.add_patch(patch)
    return ax


def draw_legend(fig) -> None:
    ax = rounded_box(fig, [0.055, 0.812, 0.890, 0.065], SOFT_PANEL, edge=PANEL_EDGE, radius=0.014)
    ax.text(0.030, 0.63, "How to read each panel", ha="left", va="center",
            fontsize=15.2, fontweight="bold", color=INK, transform=ax.transAxes)
    ax.text(0.030, 0.27, "after tight WP80 minus before preselection",
            ha="left", va="center", fontsize=12.4, color=MUTED, transform=ax.transAxes)
    items = (
        ("less common after tight", INCLUSIVE_REMOVED, 0.375),
        ("more common after tight", INCLUSIVE_ADDED, 0.595),
        ("signal change", SIGNAL, 0.805),
    )
    for label, color, x in items:
        ax.plot([x, x + 0.055], [0.50, 0.50], color=color, linewidth=4.0, transform=ax.transAxes)
        ax.text(x + 0.066, 0.50, label, ha="left", va="center",
                fontsize=12.8, color=INK, transform=ax.transAxes)


def row_limits(diffs: list[ShapeDiff]) -> dict[str, float]:
    limits: dict[str, float] = {}
    for variable in [str(v["plain"]) for v in VARIABLES]:
        vals = []
        for d in diffs:
            if d.variable == variable:
                vals.extend(abs(x) for x in d.delta_density)
        limits[variable] = max(vals) * 1.18 if vals else 1.0
    return limits


def draw_panel(
    ax,
    lk: dict[tuple[str, str, str], ShapeDiff],
    variable: dict[str, object],
    et_label: str,
    ylim: float,
) -> None:
    variable_label = str(variable["plain"])
    inc = lk[("Inclusive MC", variable_label, et_label)]
    sig = lk[("Signal MC", variable_label, et_label)]
    edges = np.array(inc.edges, dtype=float)
    delta = np.array(inc.delta_density, dtype=float)

    ax.set_facecolor("white")
    for spine in ax.spines.values():
        spine.set_color(PANEL_EDGE)
        spine.set_linewidth(0.95)

    ax.axhline(0.0, color=INK, linewidth=1.1, zorder=2)
    ax.fill_between(
        edges[:-1], 0.0, np.clip(delta, None, 0.0),
        step="post", color=INCLUSIVE_REMOVED, alpha=0.58, linewidth=0, zorder=3,
    )
    ax.fill_between(
        edges[:-1], 0.0, np.clip(delta, 0.0, None),
        step="post", color=INCLUSIVE_ADDED, alpha=0.60, linewidth=0, zorder=3,
    )
    ax.stairs(delta, edges, color=INCLUSIVE, linewidth=2.4, zorder=4)

    sig_delta = np.array(sig.delta_density, dtype=float)
    ax.stairs(sig_delta, np.array(sig.edges, dtype=float), color=SIGNAL, linewidth=1.7, zorder=5)

    ax.set_xlim(*variable["xlim"])
    ax.set_ylim(-ylim, ylim)
    ax.grid(True, color=GRID, linewidth=0.75)
    ax.tick_params(labelsize=9.6, direction="in", top=True, right=True)
    ax.tick_params(axis="y", labelsize=8.5)
    ax.set_yticks([-ylim * 0.65, 0.0, ylim * 0.65])
    ax.set_yticklabels(["less", "0", "more"])


def render_slide(diffs: list[ShapeDiff], outdir: Path) -> dict[str, Path]:
    setup_style()
    outdir.mkdir(parents=True, exist_ok=True)
    lk = diff_lookup(diffs)
    ylimits = row_limits(diffs)

    fig = plt.figure(figsize=(16, 9), dpi=160)
    fig.text(0.045, 0.956, r"Tight WP80 removes background-like shower shapes",
             ha="left", va="top", fontsize=31.5, fontweight="bold", color=INK)
    fig.text(0.046, 0.904,
             r"0-20% centrality; blue shows how the inclusive MC shape changes after the tight BDT cut in each $E_T$ bin.",
             ha="left", va="top", fontsize=16.0, color=MUTED)

    draw_legend(fig)

    left = 0.150
    panel_w = 0.252
    panel_h = 0.178
    hgap = 0.030
    vgap = 0.049
    top = 0.768
    xs = [left + i * (panel_w + hgap) for i in range(3)]
    ys = [top - panel_h - i * (panel_h + vgap) for i in range(3)]

    for col_idx, (et_label, _) in enumerate(ET_BINS):
        fig.text(xs[col_idx] + panel_w / 2, top + 0.014, et_label,
                 ha="center", va="bottom", fontsize=17.0, fontweight="bold", color=INK)

    for row_idx, variable in enumerate(VARIABLES):
        label_ax = rounded_box(fig, [0.047, ys[row_idx] + 0.016, 0.078, panel_h - 0.032],
                               "#FFFFFF", edge=PANEL_EDGE, radius=0.016)
        label_ax.text(0.50, 0.62, str(variable["plain"]), ha="center", va="center",
                      fontsize=13.8, fontweight="bold", color=INK, transform=label_ax.transAxes)
        label_ax.text(0.50, 0.32, str(variable["feature"]), ha="center", va="center",
                      fontsize=8.5, color=MUTED, transform=label_ax.transAxes)
        for col_idx, (et_label, _) in enumerate(ET_BINS):
            ax = fig.add_axes([xs[col_idx], ys[row_idx], panel_w, panel_h])
            draw_panel(ax, lk, variable, et_label, ylimits[str(variable["plain"])])
            if row_idx == len(VARIABLES) - 1:
                ax.set_xlabel(str(variable["label"]), fontsize=11.8, labelpad=2)
            else:
                ax.tick_params(labelbottom=False)
            if col_idx == 0:
                ax.set_ylabel(r"$\Delta$ shape", fontsize=10.5, labelpad=5)
            else:
                ax.tick_params(labelleft=False)

    note = rounded_box(fig, [0.055, 0.036, 0.890, 0.052], "#ECFDF5", edge="#A8DCC6", radius=0.014)
    note.text(
        0.025,
        0.50,
        r"Key point: the blue background shape changes strongly, while the red signal change stays small. "
        r"The cut mainly removes background-like low-core or low-$E_1$ shoulders.",
        ha="left",
        va="center",
        fontsize=13.0,
        color=TAKEAWAY,
        transform=note.transAxes,
    )

    png = outdir / "the42_energy_sum_feature_shape_difference_grid_slide.png"
    manifest = outdir / "the42_energy_sum_feature_shape_difference_grid_manifest.json"
    script = outdir / "the42_energy_sum_feature_shape_difference_grid_script.md"
    fig.savefig(png, dpi=160)
    plt.close(fig)
    script.write_text(
        "# THE-42 Energy-Sum BDT Input Shape-Change Grid Script\n\n"
        "This slide shows the tight-WP80 selection effect directly. In each panel I subtract the before-preselection "
        "unit-normalized shape from the tight-WP80 unit-normalized shape. The columns are ET bins inside the 0 to 20 "
        "percent centrality interval, and the rows are the energy-sum BDT inputs written as THE-42 stage histograms: "
        "cluster_et1, E11/E33, and E32/E35.\n\n"
        "The audience-facing rule is simple: blue below zero means that shower-shape region is less common after the "
        "tight cut, and blue above zero means that shower-shape region is more common after the tight cut. The red "
        "line is the same shape change for the signal sample. The key point is that the inclusive MC changes strongly, "
        "while the signal line remains comparatively small.\n\n"
        "This is cleaner than the direct overlay because the audience no longer has to disentangle four overlapping "
        "histograms in each panel. The plotted object is the selection effect itself.\n",
        encoding="utf-8",
    )
    return {"png": png, "manifest": manifest, "speaker_script": script}


def write_manifest(path: Path, diffs: list[ShapeDiff], outputs: dict[str, Path], args: argparse.Namespace) -> None:
    payload = {
        "campaign": "THE-42 WP80 centrality-linear AuAu SS overlay",
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "script": str(THIS_FILE),
        "signal_root": str(args.signal_root.resolve()),
        "inclusive_root": str(args.inclusive_root.resolve()),
        "plot_mode": "signed shape difference grid: tight WP80 density minus before-preselection density",
        "display_rebin_factor": args.rebin,
        "normalization": "each before/tight histogram is unit-normalized in the shown x range before subtraction",
        "centrality_focus": {"label": CENT_FOCUS_LABEL, "fine_bins": CENT_FOCUS_BINS},
        "et_bins": [{"label": label, "bin": bin_edges} for label, bin_edges in ET_BINS],
        "variables": [
            {
                "hist_key": v["key"],
                "bdt_feature": v["feature"],
                "label": v["plain"],
                "xlim": v["xlim"],
            }
            for v in VARIABLES
        ],
        "stage_tags": SAMPLE_TAGS,
        "outputs": {k: str(v.resolve()) for k, v in outputs.items()},
        "shape_differences": [
            {
                key: value
                for key, value in asdict(d).items()
                if key not in {"edges", "delta_density"}
            }
            for d in diffs
        ],
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--signal-root", type=Path, default=DEFAULT_SIGNAL_ROOT)
    ap.add_argument("--inclusive-root", type=Path, default=DEFAULT_INCLUSIVE_ROOT)
    ap.add_argument("--output-dir", type=Path, default=DEFAULT_OUTDIR)
    ap.add_argument("--rebin", type=int, default=3, help="Display-only rebin factor applied before shape subtraction.")
    return ap.parse_args()


def main() -> int:
    args = parse_args()
    if args.rebin < 1:
        raise SystemExit("--rebin must be >= 1")
    for path in (args.signal_root, args.inclusive_root):
        if not path.exists():
            raise SystemExit(f"Missing ROOT input: {path}")
    diffs = collect_diffs(args.signal_root, args.inclusive_root, args.rebin)
    outputs = render_slide(diffs, args.output_dir)
    write_manifest(outputs["manifest"], diffs, outputs, args)
    for key, value in outputs.items():
        print(f"{key}={value}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
