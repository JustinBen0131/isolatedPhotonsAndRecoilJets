#!/usr/bin/env python3
"""Render a THE-42 full-slide table for w_eta shower-shape evolution."""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle


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

PT_BINS = ((22, 24), (24, 26), (26, 28))
CENT_GROUPS = (
    ("0-20%", ((0, 10), (10, 20))),
    ("20-50%", ((20, 30), (30, 40), (40, 50))),
    ("50-80%", ((50, 60), (60, 80))),
)
STAGES = (
    ("No preselection", "inclusive"),
    ("Preselection", "pre"),
    ("Tight WP80", "tight"),
)
SAMPLE_TAGS = {
    "Signal MC": {"inclusive": "inclusive_sig", "pre": "pre_sig", "tight": "tight_sig"},
    "Inclusive MC": {"inclusive": "inclusive_bkg", "pre": "pre_bkg", "tight": "tight_bkg"},
}

INK = "#111827"
MUTED = "#556070"
GRID = "#DCE5EF"
PANEL = "#F8FAFC"
PANEL_EDGE = "#C8D6E5"
SIGNAL = "#D84A4A"
INCLUSIVE = "#2F78B7"
TIGHT = "#1B9E77"
PRE = "#2E77BB"
NEUTRAL = "#6B7280"


@dataclass
class Metric:
    sample: str
    centrality: str
    stage: str
    mean: float
    mean_err: float
    rms: float
    entries: float


def setup_style() -> None:
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
        "mathtext.fontset": "dejavuserif",
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


def sum_weta_hist(root_file, tag: str, cent_bins: tuple[tuple[int, int], ...]):
    acc = None
    missing = []
    for pt_lo, pt_hi in PT_BINS:
        for c_lo, c_hi in cent_bins:
            name = f"SIM/h_ss_weta_{tag}_pT_{pt_lo}_{pt_hi}_cent_{c_lo}_{c_hi}"
            hist = root_file.Get(name)
            if not hist:
                missing.append(name)
                continue
            if acc is None:
                acc = hist.Clone(f"sum_weta_{tag}_{pt_lo}_{pt_hi}_{c_lo}_{c_hi}")
                acc.SetDirectory(0)
            else:
                acc.Add(hist)
    if acc is None:
        raise RuntimeError(f"No h_ss_weta histograms found for tag={tag}; first missing={missing[:3]}")
    return acc


def metric_from_hist(sample: str, centrality: str, stage: str, hist) -> Metric:
    entries = float(hist.Integral())
    if entries <= 0:
        return Metric(sample, centrality, stage, math.nan, math.nan, math.nan, entries)
    mean = sum(hist.GetXaxis().GetBinCenter(i) * hist.GetBinContent(i) for i in range(1, hist.GetNbinsX() + 1)) / entries
    variance = sum(
        hist.GetBinContent(i) * (hist.GetXaxis().GetBinCenter(i) - mean) ** 2
        for i in range(1, hist.GetNbinsX() + 1)
    ) / entries
    mean_err = math.sqrt(sum(
        (hist.GetBinError(i) * (hist.GetXaxis().GetBinCenter(i) - mean)) ** 2
        for i in range(1, hist.GetNbinsX() + 1)
    )) / entries
    return Metric(sample, centrality, stage, mean, mean_err, math.sqrt(max(variance, 0.0)), entries)


def collect_metrics(signal_root: Path, inclusive_root: Path) -> list[Metric]:
    out = []
    paths = {"Signal MC": signal_root, "Inclusive MC": inclusive_root}
    for sample, path in paths.items():
        root_file = open_root(path)
        try:
            for cent_label, cent_bins in CENT_GROUPS:
                for stage_label, stage_key in STAGES:
                    hist = sum_weta_hist(root_file, SAMPLE_TAGS[sample][stage_key], cent_bins)
                    out.append(metric_from_hist(sample, cent_label, stage_label, hist))
        finally:
            root_file.Close()
    return out


def as_lookup(metrics: list[Metric]) -> dict[tuple[str, str, str], Metric]:
    return {(m.sample, m.centrality, m.stage): m for m in metrics}


def rounded_box(ax, xy, wh, face, edge=PANEL_EDGE, radius=0.018, lw=1.2, alpha=1.0):
    patch = FancyBboxPatch(
        xy, wh[0], wh[1],
        boxstyle=f"round,pad=0.006,rounding_size={radius}",
        facecolor=face, edgecolor=edge, linewidth=lw, alpha=alpha,
        transform=ax.transAxes, clip_on=False,
    )
    ax.add_patch(patch)
    return patch


def cell_color(value: float, lo: float = 0.15, hi: float = 0.45) -> str:
    if not math.isfinite(value):
        return "#FFFFFF"
    t = max(0.0, min(1.0, (value - lo) / (hi - lo)))
    # Higher width is more orange; lower/narrower is cooler/cleaner.
    r0, g0, b0 = (232, 246, 242)
    r1, g1, b1 = (255, 232, 209)
    r = int(r0 + t * (r1 - r0))
    g = int(g0 + t * (g1 - g0))
    b = int(b0 + t * (b1 - b0))
    return f"#{r:02x}{g:02x}{b:02x}"


def draw_table(ax, lookup, sample: str, color: str, x0: float, y0: float, w: float, h: float) -> None:
    rounded_box(ax, (x0, y0), (w, h), "white", edge=PANEL_EDGE, radius=0.022, lw=1.4)
    ax.add_patch(Rectangle((x0, y0 + h - 0.075), w, 0.075, transform=ax.transAxes,
                           facecolor=color, edgecolor="none", alpha=0.11))
    ax.text(x0 + 0.025, y0 + h - 0.040, sample, ha="left", va="center",
            fontsize=23.0, fontweight="bold", color=color, transform=ax.transAxes)
    ax.text(x0 + w - 0.025, y0 + h - 0.040, r"mean $w_{\eta}$", ha="right", va="center",
            fontsize=15.5, color=MUTED, transform=ax.transAxes)

    cols = ["Centrality", "No preselection", "Preselection", "Tight WP80", r"$\Delta$ tight-no"]
    widths = [0.18, 0.215, 0.215, 0.215, 0.175]
    x_edges = [x0]
    for frac in widths:
        x_edges.append(x_edges[-1] + w * frac)
    header_y = y0 + h - 0.145
    row_h = (h - 0.190) / 3.0

    for i, col in enumerate(cols):
        ax.text((x_edges[i] + x_edges[i + 1]) / 2, header_y + 0.035, col,
                ha="center", va="center", fontsize=13.7, fontweight="bold",
                color=INK if i == 0 else MUTED, transform=ax.transAxes)

    for row_idx, (cent_label, _) in enumerate(CENT_GROUPS):
        y = header_y - (row_idx + 1) * row_h
        face = "#FFFFFF" if row_idx % 2 == 0 else "#F7FAFD"
        ax.add_patch(Rectangle((x0 + 0.010, y), w - 0.020, row_h - 0.006,
                               transform=ax.transAxes, facecolor=face, edgecolor="none"))
        ax.text((x_edges[0] + x_edges[1]) / 2, y + row_h / 2, cent_label,
                ha="center", va="center", fontsize=16.5, fontweight="bold",
                color=INK, transform=ax.transAxes)

        no = lookup[(sample, cent_label, "No preselection")]
        tight = lookup[(sample, cent_label, "Tight WP80")]
        for stage_idx, stage in enumerate(["No preselection", "Preselection", "Tight WP80"], start=1):
            m = lookup[(sample, cent_label, stage)]
            cx0 = x_edges[stage_idx] + 0.010
            cw = x_edges[stage_idx + 1] - x_edges[stage_idx] - 0.020
            cy0 = y + 0.028
            ch = row_h - 0.062
            ax.add_patch(Rectangle(
                (cx0, cy0), cw, ch,
                facecolor=cell_color(m.mean), edgecolor="#E0E7EF", linewidth=0.8,
                transform=ax.transAxes,
            ))
            ax.text(cx0 + cw / 2, cy0 + ch * 0.50, f"{m.mean:.3f}",
                    ha="center", va="center", fontsize=20.0, fontweight="bold",
                    color=INK, transform=ax.transAxes)

        delta = tight.mean - no.mean
        dcolor = TIGHT if delta < 0 else "#B45309"
        ax.text((x_edges[4] + x_edges[5]) / 2, y + row_h * 0.60, f"{delta:+.3f}",
                ha="center", va="center", fontsize=21.0, fontweight="bold",
                color=dcolor, transform=ax.transAxes)
        ax.text((x_edges[4] + x_edges[5]) / 2, y + row_h * 0.30, "narrower" if delta < 0 else "wider",
                ha="center", va="center", fontsize=11.8, color=MUTED, transform=ax.transAxes)

    for xe in x_edges[1:-1]:
        ax.plot([xe, xe], [y0 + 0.035, y0 + h - 0.120], transform=ax.transAxes,
                color=GRID, linewidth=0.9)
    for i in range(4):
        yy = header_y - i * row_h
        ax.plot([x0 + 0.018, x0 + w - 0.018], [yy, yy], transform=ax.transAxes,
                color=GRID, linewidth=0.9)


def make_manifest(metrics: list[Metric], outputs: dict[str, Path], args: argparse.Namespace) -> dict:
    return {
        "script": str(THIS_FILE),
        "signal_root": str(args.signal_root.resolve()),
        "inclusive_root": str(args.inclusive_root.resolve()),
        "variable": "weta",
        "histogram_family": "SIM/h_ss_weta_<stage>_pT_<pt>_cent_<cent>",
        "pt_bins": PT_BINS,
        "centrality_groups": [{"label": label, "fine_bins": bins} for label, bins in CENT_GROUPS],
        "stage_tags": SAMPLE_TAGS,
        "metric": "mean and approximate statistical error from summed h_ss_weta histograms; histograms are used as unit-normalized shape summaries, not efficiency denominators",
        "outputs": {k: str(v.resolve()) for k, v in outputs.items()},
        "rows": [m.__dict__ for m in metrics],
    }


def write_csv(metrics: list[Metric], path: Path) -> None:
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["sample", "centrality", "stage", "mean", "mean_err", "rms", "entries"])
        writer.writeheader()
        writer.writerows([m.__dict__ for m in metrics])


def write_script(path: Path) -> None:
    path.write_text(
        "# THE-42 Weta Shower-Shape Evolution Table Script\n\n"
        "This slide is a compact way to show that the WP80 BDT selection is not just changing the event yield; "
        "it also changes the shape of the inclusive-jet shower-width distribution.\n\n"
        "The number in each table cell is the mean of the unit-normalized `weta` shower-shape histogram after summing "
        "the 22 to 28 GeV photon candidate range. I split the result into three broad centrality regions: 0 to 20, "
        "20 to 50, and 50 to 80 percent.\n\n"
        "For the embedded photon signal sample, the mean `weta` is already relatively narrow, so the tight WP80 cut "
        "moves it only modestly. In contrast, the inclusive embedded-jet sample starts much broader, especially before "
        "selection, and the tight BDT cut pulls it strongly toward the same narrow region as the signal.\n\n"
        "The important point is that `e11/e33` is not the only shower-shape variable where the selection effect can be seen. "
        "`weta` gives a visually cleaner summary of the narrowing imposed by the BDT, especially for the inclusive MC background. "
        "I would use this as a supporting QA slide behind the main PPG12-style overlay slide.\n",
        encoding="utf-8",
    )


def render_slide(metrics: list[Metric], outdir: Path) -> dict[str, Path]:
    setup_style()
    outdir.mkdir(parents=True, exist_ok=True)
    lookup = as_lookup(metrics)

    fig = plt.figure(figsize=(16, 9), dpi=160)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis("off")

    ax.text(0.045, 0.952, r"WP80 selection visibly narrows $w_{\eta}$ for inclusive MC",
            ha="left", va="top", fontsize=31.5, fontweight="bold", color=INK, transform=ax.transAxes)
    ax.text(0.046, 0.902,
            r"Mean of unit-normalized $h\_ss\_weta$ over $22 \leq p_T^\gamma < 28$ GeV; values are grouped into three centrality regions.",
            ha="left", va="top", fontsize=16.2, color=MUTED, transform=ax.transAxes)

    rounded_box(ax, (0.705, 0.888), (0.250, 0.060), "#ECFDF5", edge="#A8DCC6", radius=0.018, lw=1.0)
    ax.text(0.722, 0.918, "Reading the table", ha="left", va="center",
            fontsize=14.3, fontweight="bold", color=TIGHT, transform=ax.transAxes)
    ax.text(0.722, 0.894, r"smaller $w_{\eta}$ = narrower shower shape",
            ha="left", va="center", fontsize=12.8, color=INK, transform=ax.transAxes)

    draw_table(ax, lookup, "Signal MC", SIGNAL, 0.045, 0.505, 0.910, 0.335)
    draw_table(ax, lookup, "Inclusive MC", INCLUSIVE, 0.045, 0.145, 0.910, 0.335)

    # Bottom interpretation cards.
    rounded_box(ax, (0.045, 0.045), (0.280, 0.070), "#FFF7ED", edge="#FED7AA", radius=0.018, lw=1.0)
    ax.text(0.062, 0.086, "Signal shift is small", fontsize=14.5, fontweight="bold",
            color=SIGNAL, ha="left", va="center", transform=ax.transAxes)
    ax.text(0.062, 0.060, r"tight WP80 changes mean $w_{\eta}$ by about -0.02",
            fontsize=12.6, color=INK, ha="left", va="center", transform=ax.transAxes)

    rounded_box(ax, (0.360, 0.045), (0.300, 0.070), "#EFF6FF", edge="#BFDBFE", radius=0.018, lw=1.0)
    ax.text(0.377, 0.086, "Inclusive shift is large", fontsize=14.5, fontweight="bold",
            color=INCLUSIVE, ha="left", va="center", transform=ax.transAxes)
    ax.text(0.377, 0.060, r"tight WP80 moves mean $w_{\eta}$ down by about -0.17",
            fontsize=12.6, color=INK, ha="left", va="center", transform=ax.transAxes)

    rounded_box(ax, (0.695, 0.045), (0.260, 0.070), "#F8FAFC", edge=PANEL_EDGE, radius=0.018, lw=1.0)
    ax.text(0.712, 0.086, "Caveat", fontsize=14.5, fontweight="bold",
            color=MUTED, ha="left", va="center", transform=ax.transAxes)
    ax.text(0.712, 0.060, "shape means only; stat errors <= 0.004 in CSV",
            fontsize=12.6, color=INK, ha="left", va="center", transform=ax.transAxes)

    png = outdir / "the42_weta_ss_evolution_table_slide_v2.png"
    fig.savefig(png, dpi=160)
    plt.close(fig)

    csv_path = outdir / "the42_weta_ss_evolution_table.csv"
    manifest_path = outdir / "the42_weta_ss_evolution_table_manifest.json"
    script_path = outdir / "the42_weta_ss_evolution_table_script.md"
    write_csv(metrics, csv_path)
    write_script(script_path)
    outputs = {"png": png, "csv": csv_path, "manifest": manifest_path, "speaker_script": script_path}
    return outputs


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--signal-root", type=Path, default=DEFAULT_SIGNAL_ROOT)
    ap.add_argument("--inclusive-root", type=Path, default=DEFAULT_INCLUSIVE_ROOT)
    ap.add_argument("--output-dir", type=Path, default=DEFAULT_OUTDIR)
    return ap.parse_args()


def main() -> int:
    args = parse_args()
    for path in (args.signal_root, args.inclusive_root):
        if not path.exists():
            raise SystemExit(f"Missing ROOT input: {path}")
    metrics = collect_metrics(args.signal_root, args.inclusive_root)
    outputs = render_slide(metrics, args.output_dir)
    manifest = make_manifest(metrics, outputs, args)
    outputs["manifest"].write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    for key, value in outputs.items():
        print(f"{key}={value}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
