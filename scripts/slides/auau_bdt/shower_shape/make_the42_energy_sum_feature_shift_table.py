#!/usr/bin/env python3
"""Render a THE-42 table for energy-sum BDT input shape shifts."""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import asdict, dataclass
from datetime import datetime
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
    ("Before preselection", "inclusive"),
    ("Tight WP80", "tight"),
)
SAMPLE_TAGS = {
    "Signal MC": {"inclusive": "inclusive_sig", "tight": "tight_sig"},
    "Inclusive MC": {"inclusive": "inclusive_bkg", "tight": "tight_bkg"},
}

VARIABLES = (
    {
        "hist_key": "et1",
        "feature": "cluster_et1",
        "label": "cluster_et1",
        "plain": "cluster_et1",
        "description": "leading tower energy sharing",
    },
    {
        "hist_key": "e11e33",
        "feature": "e11_over_e33",
        "label": "E11 / E33",
        "plain": "E11/E33",
        "description": "central core over 3x3 core",
    },
    {
        "hist_key": "e32e35",
        "feature": "e32_over_e35",
        "label": "E32 / E35",
        "plain": "E32/E35",
        "description": "3x2 strip over 3x5 region",
    },
)

UNWRITTEN_STAGE_FEATURES = (
    "cluster_et2",
    "cluster_et3",
    "cluster_et4",
    "e11_over_e22",
    "e11_over_e13",
    "e11_over_e15",
    "e11_over_e17",
    "e11_over_e31",
    "e11_over_e51",
    "e11_over_e71",
    "e22_over_e33",
    "e22_over_e35",
    "e22_over_e37",
    "e22_over_e53",
)

INK = "#111827"
MUTED = "#566274"
LIGHT_TEXT = "#6B7280"
GRID = "#DCE5EF"
PANEL = "#F8FAFC"
PANEL_EDGE = "#C8D6E5"
SIGNAL = "#D84A4A"
INCLUSIVE = "#2F78B7"
POS = "#126F55"
NEG = "#B45309"


@dataclass(frozen=True)
class Metric:
    sample: str
    variable: str
    feature: str
    centrality: str
    stage: str
    mean: float
    mean_stat_error: float
    rms: float
    entries: float


@dataclass(frozen=True)
class SummaryRow:
    sample: str
    variable: str
    feature: str
    centrality: str
    before_mean: float
    before_mean_stat_error: float
    before_rms: float
    before_entries: float
    tight_mean: float
    tight_mean_stat_error: float
    tight_rms: float
    tight_entries: float
    delta_tight_minus_before: float
    relative_delta_percent: float


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


def sum_hist(root_file, hist_key: str, tag: str, cent_bins: tuple[tuple[int, int], ...]):
    acc = None
    missing = []
    for pt_lo, pt_hi in PT_BINS:
        for c_lo, c_hi in cent_bins:
            name = f"SIM/h_ss_{hist_key}_{tag}_pT_{pt_lo}_{pt_hi}_cent_{c_lo}_{c_hi}"
            hist = root_file.Get(name)
            if not hist:
                missing.append(name)
                continue
            if acc is None:
                acc = hist.Clone(f"sum_{hist_key}_{tag}_{pt_lo}_{pt_hi}_{c_lo}_{c_hi}")
                acc.SetDirectory(0)
            else:
                acc.Add(hist)
    if acc is None:
        raise RuntimeError(f"No h_ss_{hist_key} histograms found for tag={tag}; first missing={missing[:3]}")
    return acc


def metric_from_hist(sample: str, variable: dict[str, str], centrality: str, stage: str, hist) -> Metric:
    entries = float(hist.Integral())
    if entries <= 0:
        return Metric(
            sample=sample,
            variable=variable["plain"],
            feature=variable["feature"],
            centrality=centrality,
            stage=stage,
            mean=math.nan,
            mean_stat_error=math.nan,
            rms=math.nan,
            entries=entries,
        )
    mean = sum(hist.GetXaxis().GetBinCenter(i) * hist.GetBinContent(i) for i in range(1, hist.GetNbinsX() + 1)) / entries
    variance = sum(
        hist.GetBinContent(i) * (hist.GetXaxis().GetBinCenter(i) - mean) ** 2
        for i in range(1, hist.GetNbinsX() + 1)
    ) / entries
    # Propagate the stored per-bin statistical uncertainties into the weighted mean.
    mean_err = math.sqrt(sum(
        (hist.GetBinError(i) * (hist.GetXaxis().GetBinCenter(i) - mean)) ** 2
        for i in range(1, hist.GetNbinsX() + 1)
    )) / entries
    return Metric(
        sample=sample,
        variable=variable["plain"],
        feature=variable["feature"],
        centrality=centrality,
        stage=stage,
        mean=mean,
        mean_stat_error=mean_err,
        rms=math.sqrt(max(variance, 0.0)),
        entries=entries,
    )


def collect_metrics(signal_root: Path, inclusive_root: Path) -> list[Metric]:
    out: list[Metric] = []
    for sample, path in {"Signal MC": signal_root, "Inclusive MC": inclusive_root}.items():
        root_file = open_root(path)
        try:
            for variable in VARIABLES:
                for cent_label, cent_bins in CENT_GROUPS:
                    for stage_label, stage_key in STAGES:
                        tag = SAMPLE_TAGS[sample][stage_key]
                        hist = sum_hist(root_file, variable["hist_key"], tag, cent_bins)
                        out.append(metric_from_hist(sample, variable, cent_label, stage_label, hist))
        finally:
            root_file.Close()
    return out


def summarize(metrics: list[Metric]) -> list[SummaryRow]:
    lookup = {(m.sample, m.variable, m.centrality, m.stage): m for m in metrics}
    rows: list[SummaryRow] = []
    for sample in SAMPLE_TAGS:
        for variable in VARIABLES:
            for cent_label, _ in CENT_GROUPS:
                before = lookup[(sample, variable["plain"], cent_label, "Before preselection")]
                tight = lookup[(sample, variable["plain"], cent_label, "Tight WP80")]
                delta = tight.mean - before.mean
                relative = 100.0 * delta / before.mean if before.mean and math.isfinite(before.mean) else math.nan
                rows.append(SummaryRow(
                    sample=sample,
                    variable=variable["plain"],
                    feature=variable["feature"],
                    centrality=cent_label,
                    before_mean=before.mean,
                    before_mean_stat_error=before.mean_stat_error,
                    before_rms=before.rms,
                    before_entries=before.entries,
                    tight_mean=tight.mean,
                    tight_mean_stat_error=tight.mean_stat_error,
                    tight_rms=tight.rms,
                    tight_entries=tight.entries,
                    delta_tight_minus_before=delta,
                    relative_delta_percent=relative,
                ))
    return rows


def rounded_box(ax, xy, wh, face, edge=PANEL_EDGE, radius=0.018, lw=1.0, alpha=1.0):
    patch = FancyBboxPatch(
        xy, wh[0], wh[1],
        boxstyle=f"round,pad=0.006,rounding_size={radius}",
        facecolor=face,
        edgecolor=edge,
        linewidth=lw,
        alpha=alpha,
        transform=ax.transAxes,
        clip_on=False,
    )
    ax.add_patch(patch)
    return patch


def fmt(value: float, digits: int = 3) -> str:
    if not math.isfinite(value):
        return "n/a"
    return f"{value:.{digits}f}"


def fmt_delta(value: float) -> str:
    if not math.isfinite(value):
        return "n/a"
    return f"{value:+.3f}"


def row_band_color(sample: str, delta: float) -> str:
    if sample == "Signal MC":
        return "#FFF8F7" if abs(delta) > 0.05 else "#FFFFFF"
    return "#EFF6FF" if abs(delta) > 0.05 else "#FFFFFF"


def draw_sample_table(ax, rows: list[SummaryRow], sample: str, color: str, x0: float, y0: float, w: float, h: float) -> None:
    sample_rows = [r for r in rows if r.sample == sample]
    rounded_box(ax, (x0, y0), (w, h), "white", edge=PANEL_EDGE, radius=0.020, lw=1.3)
    ax.add_patch(Rectangle((x0, y0 + h - 0.073), w, 0.073, transform=ax.transAxes,
                           facecolor=color, edgecolor="none", alpha=0.10))
    ax.text(x0 + 0.020, y0 + h - 0.038, sample, ha="left", va="center",
            fontsize=22.5, fontweight="bold", color=color, transform=ax.transAxes)
    ax.text(x0 + w - 0.020, y0 + h - 0.038, "mean of stage histogram",
            ha="right", va="center", fontsize=13.8, color=MUTED, transform=ax.transAxes)

    cols = [
        ("BDT input", 0.225),
        ("Cent.", 0.135),
        ("Before", 0.185),
        ("Tight", 0.165),
        ("Delta", 0.155),
        ("Rel.", 0.135),
    ]
    x_edges = [x0]
    for _, frac in cols:
        x_edges.append(x_edges[-1] + w * frac)
    header_y = y0 + h - 0.132
    usable_h = h - 0.160
    row_h = usable_h / len(sample_rows)

    for i, (label, _) in enumerate(cols):
        ax.text((x_edges[i] + x_edges[i + 1]) / 2, header_y + 0.028, label,
                ha="center", va="center", fontsize=13.6, fontweight="bold",
                color=INK if i < 2 else MUTED, transform=ax.transAxes)

    prev_variable = None
    for idx, row in enumerate(sample_rows):
        y = header_y - (idx + 1) * row_h
        face = row_band_color(sample, row.delta_tight_minus_before)
        if idx % 2 == 1 and face == "#FFFFFF":
            face = "#F8FAFC"
        ax.add_patch(Rectangle((x0 + 0.010, y), w - 0.020, row_h - 0.004,
                               transform=ax.transAxes, facecolor=face, edgecolor="none"))

        show_variable = row.variable != prev_variable
        if show_variable:
            ax.text(x_edges[0] + 0.015, y + row_h * 0.62, row.variable,
                    ha="left", va="center", fontsize=15.2, fontweight="bold",
                    color=INK, transform=ax.transAxes)
            desc = next(v["description"] for v in VARIABLES if v["plain"] == row.variable)
            ax.text(x_edges[0] + 0.015, y + row_h * 0.30, desc,
                    ha="left", va="center", fontsize=9.6, color=LIGHT_TEXT,
                    transform=ax.transAxes)
            prev_variable = row.variable

        ax.text((x_edges[1] + x_edges[2]) / 2, y + row_h * 0.50, row.centrality,
                ha="center", va="center", fontsize=13.8, fontweight="bold",
                color=INK, transform=ax.transAxes)
        ax.text((x_edges[2] + x_edges[3]) / 2, y + row_h * 0.50, fmt(row.before_mean),
                ha="center", va="center", fontsize=15.6, fontweight="bold",
                color=INK, transform=ax.transAxes)
        ax.text((x_edges[3] + x_edges[4]) / 2, y + row_h * 0.50, fmt(row.tight_mean),
                ha="center", va="center", fontsize=15.6, fontweight="bold",
                color=INK, transform=ax.transAxes)
        delta_color = POS if row.delta_tight_minus_before >= 0 else NEG
        ax.text((x_edges[4] + x_edges[5]) / 2, y + row_h * 0.50, fmt_delta(row.delta_tight_minus_before),
                ha="center", va="center", fontsize=15.8, fontweight="bold",
                color=delta_color, transform=ax.transAxes)
        ax.text((x_edges[5] + x_edges[6]) / 2, y + row_h * 0.50, f"{row.relative_delta_percent:+.1f}%",
                ha="center", va="center", fontsize=13.2, fontweight="bold",
                color=delta_color, transform=ax.transAxes)

        if idx in (2, 5):
            ax.plot([x0 + 0.012, x0 + w - 0.012], [y, y], transform=ax.transAxes,
                    color=color, linewidth=1.1, alpha=0.30)

    for xe in x_edges[1:-1]:
        ax.plot([xe, xe], [y0 + 0.025, y0 + h - 0.110], transform=ax.transAxes,
                color=GRID, linewidth=0.8)
    for i in range(len(sample_rows) + 1):
        yy = header_y - i * row_h
        ax.plot([x0 + 0.012, x0 + w - 0.012], [yy, yy], transform=ax.transAxes,
                color=GRID, linewidth=0.7)


def draw_feature_note(ax) -> None:
    rounded_box(ax, (0.045, 0.815), (0.910, 0.062), "#F3F8FF", edge="#BED3EE", radius=0.018, lw=1.0)
    ax.text(0.065, 0.848,
            "Energy-sum stage histograms available for the BDT inputs: cluster_et1, E11/E33, and E32/E35.",
            ha="left", va="center", fontsize=15.8, color=INK, transform=ax.transAxes)
    ax.text(0.065, 0.824,
            "Means are computed before preselection and after the centrality-linear WP80 tight cut.",
            ha="left", va="center", fontsize=13.5, color=MUTED, transform=ax.transAxes)


def draw_bottom_notes(ax, rows: list[SummaryRow]) -> None:
    inclusive = [r for r in rows if r.sample == "Inclusive MC"]
    largest = max(inclusive, key=lambda r: abs(r.delta_tight_minus_before))
    e11 = [r for r in inclusive if r.variable == "E11/E33"]
    e32 = [r for r in inclusive if r.variable == "E32/E35"]

    rounded_box(ax, (0.045, 0.035), (0.275, 0.096), "#FFF7ED", edge="#FED7AA", radius=0.018, lw=1.0)
    ax.text(0.062, 0.097, "Largest inclusive response", fontsize=14.8,
            fontweight="bold", color=NEG, ha="left", va="center", transform=ax.transAxes)
    ax.text(0.062, 0.069, f"{largest.variable} shifts by {fmt_delta(largest.delta_tight_minus_before)}",
            fontsize=13.3, color=INK, ha="left", va="center", transform=ax.transAxes)

    rounded_box(ax, (0.355, 0.035), (0.285, 0.096), "#EFF6FF", edge="#BFDBFE", radius=0.018, lw=1.0)
    ax.text(0.372, 0.097, "Core-ratio movement", fontsize=14.8,
            fontweight="bold", color=INCLUSIVE, ha="left", va="center", transform=ax.transAxes)
    ax.text(0.372, 0.071,
            f"E11/E33: {fmt_delta(min(r.delta_tight_minus_before for r in e11))} to {fmt_delta(max(r.delta_tight_minus_before for r in e11))}",
            fontsize=13.3, color=INK, ha="left", va="center", transform=ax.transAxes)
    ax.text(0.372, 0.052,
            f"E32/E35: {fmt_delta(min(r.delta_tight_minus_before for r in e32))} to {fmt_delta(max(r.delta_tight_minus_before for r in e32))}",
            fontsize=12.4, color=MUTED, ha="left", va="center", transform=ax.transAxes)

    rounded_box(ax, (0.675, 0.035), (0.280, 0.096), "#F8FAFC", edge=PANEL_EDGE, radius=0.018, lw=1.0)
    ax.text(0.692, 0.097, "Output limitation", fontsize=14.8,
            fontweight="bold", color=MUTED, ha="left", va="center", transform=ax.transAxes)
    ax.text(0.692, 0.071, "cluster_et2-4 and extended ratios are",
            fontsize=12.7, color=INK, ha="left", va="center", transform=ax.transAxes)
    ax.text(0.692, 0.052, "training inputs but not stage histograms here.",
            fontsize=12.7, color=INK, ha="left", va="center", transform=ax.transAxes)


def write_csv(rows: list[SummaryRow], path: Path) -> None:
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(asdict(rows[0]).keys()))
        writer.writeheader()
        writer.writerows([asdict(row) for row in rows])


def write_script(path: Path, rows: list[SummaryRow]) -> None:
    inclusive = [r for r in rows if r.sample == "Inclusive MC"]
    signal = [r for r in rows if r.sample == "Signal MC"]
    strongest = max(inclusive, key=lambda r: abs(r.delta_tight_minus_before))
    path.write_text(
        "# THE-42 Energy-Sum BDT Input Shift Table Script\n\n"
        "This slide summarizes the stage-evolution histograms that are available for energy-sum BDT inputs in "
        "the THE-42 merged simulation output. The columns compare the mean before preselection to the mean after "
        "the centrality-linear WP80 tight selection, summed over 22 to 28 GeV and grouped into 0-20, 20-50, and "
        "50-80 percent centrality regions.\n\n"
        "The signal sample is stable: the largest signal shifts are small compared with the inclusive-jet shifts. "
        "The inclusive sample moves much more strongly. The largest inclusive change in this table is "
        f"{strongest.variable} in {strongest.centrality}, with delta mean {strongest.delta_tight_minus_before:+.3f}. "
        "This is the clearest compact evidence that the tight WP80 selection is changing the background shower-shape "
        "composition rather than just reducing yield.\n\n"
        "A limitation is explicit on the slide: `cluster_et2`, `cluster_et3`, `cluster_et4`, and the extended "
        "ratio features from the larger training list are not written as no-preselection/tight stage histograms "
        "in this THE-42 overlay output. They would require either a training-tree-level table or an additional "
        "histogram-writing pass.\n\n"
        f"Signal rows summarized: {len(signal)}. Inclusive rows summarized: {len(inclusive)}.\n",
        encoding="utf-8",
    )


def make_manifest(metrics: list[Metric], rows: list[SummaryRow], outputs: dict[str, Path], args: argparse.Namespace) -> dict:
    return {
        "campaign": "THE-42 WP80 centrality-linear AuAu SS overlay",
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "script": str(THIS_FILE),
        "signal_root": str(args.signal_root.resolve()),
        "inclusive_root": str(args.inclusive_root.resolve()),
        "pt_bins": PT_BINS,
        "centrality_groups": [{"label": label, "fine_bins": bins} for label, bins in CENT_GROUPS],
        "stage_comparison": ["Before preselection", "Tight WP80"],
        "stage_tags": SAMPLE_TAGS,
        "variables_in_table": [
            {
                "hist_key": v["hist_key"],
                "bdt_feature": v["feature"],
                "label": v["plain"],
                "description": v["description"],
            }
            for v in VARIABLES
        ],
        "bdt_energy_sum_features_not_written_as_stage_histograms": list(UNWRITTEN_STAGE_FEATURES),
        "metric": (
            "Means, RMS values, integrals, and propagated mean statistical uncertainties are computed "
            "from summed SIM/h_ss_<var>_<stage>_pT_<lo>_<hi>_cent_<lo>_<hi> histograms."
        ),
        "outputs": {k: str(v.resolve()) for k, v in outputs.items()},
        "metrics_by_stage": [asdict(m) for m in metrics],
        "summary_rows": [asdict(row) for row in rows],
    }


def render_slide(rows: list[SummaryRow], outdir: Path) -> dict[str, Path]:
    setup_style()
    outdir.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=(16, 9), dpi=160)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis("off")

    ax.text(0.045, 0.955, "Energy-sum BDT inputs shift under tight WP80",
            ha="left", va="top", fontsize=31.5, fontweight="bold", color=INK, transform=ax.transAxes)
    ax.text(0.046, 0.898,
            "THE-42 simulation, Photon+Jet embedded 12+20 and Inclusive Jet embedded 12+20+30+40; "
            "22 <= cluster pT < 28 GeV.",
            ha="left", va="top", fontsize=15.8, color=MUTED, transform=ax.transAxes)

    draw_feature_note(ax)
    draw_sample_table(ax, rows, "Signal MC", SIGNAL, 0.045, 0.175, 0.438, 0.610)
    draw_sample_table(ax, rows, "Inclusive MC", INCLUSIVE, 0.517, 0.175, 0.438, 0.610)
    draw_bottom_notes(ax, rows)

    png = outdir / "the42_energy_sum_feature_shift_table_slide.png"
    fig.savefig(png, dpi=160)
    plt.close(fig)

    csv_path = outdir / "the42_energy_sum_feature_shift_table.csv"
    manifest_path = outdir / "the42_energy_sum_feature_shift_table_manifest.json"
    script_path = outdir / "the42_energy_sum_feature_shift_table_script.md"
    write_csv(rows, csv_path)
    write_script(script_path, rows)
    return {"png": png, "csv": csv_path, "manifest": manifest_path, "speaker_script": script_path}


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
    rows = summarize(metrics)
    outputs = render_slide(rows, args.output_dir)
    manifest = make_manifest(metrics, rows, outputs, args)
    outputs["manifest"].write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    for key, value in outputs.items():
        print(f"{key}={value}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
