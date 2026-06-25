#!/usr/bin/env python3
"""Render THE-42 energy-sum BDT input distributions before and after tight WP80."""

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

PT_GROUPS = (
    (r"Low $p_T$: 15-22 GeV", ((15, 18), (18, 20), (20, 22))),
    (r"Mid $p_T$: 22-28 GeV", ((22, 24), (24, 26), (26, 28))),
    (r"High $p_T$: 28-35 GeV", ((28, 30), (30, 35))),
)
PT_COLUMN_STYLES = (
    {"face": "#F1F8FF", "edge": "#8DBBEA", "accent": "#2563EB"},
    {"face": "#FFF8E6", "edge": "#E7C772", "accent": "#B7791F"},
    {"face": "#F7F1FF", "edge": "#C4A8F3", "accent": "#7C3AED"},
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
STAGES = (
    ("Before preselection", "inclusive", "--"),
    ("Tight WP80", "tight", "-"),
)
SAMPLE_TAGS = {
    "Signal MC": {"inclusive": "inclusive_sig", "tight": "tight_sig"},
    "Inclusive MC": {"inclusive": "inclusive_bkg", "tight": "tight_bkg"},
}
SAMPLE_COLORS = {
    "Signal MC": {"before": "#E4574F", "tight": "#B91C1C", "soft": "#FFF4F2"},
    "Inclusive MC": {"before": "#2F7FC8", "tight": "#174EA6", "soft": "#F0F6FF"},
}

INK = "#111827"
MUTED = "#566274"
GRID = "#E4EBF3"
PANEL_EDGE = "#C8D6E5"
SOFT_PANEL = "#F8FAFC"
TAKEAWAY = "#0F766E"
OVERLAP = "#7C3AED"


@dataclass(frozen=True)
class Curve:
    sample: str
    variable: str
    feature: str
    pt_group: str
    centrality: str
    stage: str
    edges: list[float]
    density: list[float]
    visible_integral: float
    total_integral: float
    mean: float
    rms: float


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


def sum_hist(
    root_file,
    variable_key: str,
    tag: str,
    pt_bins: tuple[tuple[int, int], ...],
    cent_bins: tuple[tuple[int, int], ...],
):
    acc = None
    missing = []
    for pt_lo, pt_hi in pt_bins:
        for c_lo, c_hi in cent_bins:
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


def curve_from_hist(
    sample: str,
    variable: dict[str, object],
    pt_label: str,
    stage: str,
    hist,
    rebin: int,
) -> Curve:
    nbins = hist.GetNbinsX()
    raw_edges = np.array([hist.GetXaxis().GetBinLowEdge(i) for i in range(1, nbins + 2)], dtype=float)
    raw_centers = np.array([hist.GetXaxis().GetBinCenter(i) for i in range(1, nbins + 1)], dtype=float)
    raw_counts = np.array([hist.GetBinContent(i) for i in range(1, nbins + 1)], dtype=float)
    total_integral = float(np.sum(raw_counts))
    mean = float(np.sum(raw_centers * raw_counts) / total_integral) if total_integral > 0 else math.nan
    rms = (
        float(np.sqrt(np.sum(raw_counts * (raw_centers - mean) ** 2) / total_integral))
        if total_integral > 0
        else math.nan
    )

    edges, counts = rebin_for_display(raw_edges, raw_counts, rebin)
    left, right = variable["xlim"]
    centers = 0.5 * (edges[:-1] + edges[1:])
    keep = (centers >= left) & (centers <= right)
    edges_kept = edges[np.r_[keep, False]]
    if len(edges_kept) == 0:
        edges_kept = edges
        counts_kept = counts
    else:
        first = int(np.where(keep)[0][0])
        last = int(np.where(keep)[0][-1])
        edges_kept = edges[first:last + 2]
        counts_kept = counts[first:last + 1]

    visible = float(np.sum(counts_kept))
    widths = np.diff(edges_kept)
    density = counts_kept / (visible * widths) if visible > 0 else np.zeros_like(counts_kept)
    return Curve(
        sample=sample,
        variable=str(variable["plain"]),
        feature=str(variable["feature"]),
        pt_group=pt_label,
        centrality=CENT_FOCUS_LABEL,
        stage=stage,
        edges=[float(x) for x in edges_kept],
        density=[float(y) for y in density],
        visible_integral=visible,
        total_integral=total_integral,
        mean=mean,
        rms=rms,
    )


def collect_curves(signal_root: Path, inclusive_root: Path, rebin: int) -> list[Curve]:
    curves: list[Curve] = []
    root_paths = {"Signal MC": signal_root, "Inclusive MC": inclusive_root}
    for sample, path in root_paths.items():
        f = open_root(path)
        try:
            for variable in VARIABLES:
                for pt_label, pt_bins in PT_GROUPS:
                    for stage_label, stage_key, _ in STAGES:
                        tag = SAMPLE_TAGS[sample][stage_key]
                        hist = sum_hist(f, str(variable["key"]), tag, pt_bins, CENT_FOCUS_BINS)
                        curves.append(curve_from_hist(sample, variable, pt_label, stage_label, hist, rebin))
        finally:
            f.Close()
    return curves


def load_curves_json(path: Path) -> tuple[list[Curve], dict[str, object]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    curves = [Curve(**item) for item in payload.get("curves", [])]
    if not curves:
        raise RuntimeError(f"No curves found in {path}")
    return curves, payload


def curve_lookup(curves: list[Curve]) -> dict[tuple[str, str, str, str], Curve]:
    return {(c.sample, c.variable, c.pt_group, c.stage): c for c in curves}


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


def draw_legend(fig, *, source_payload: dict[str, object] | None = None) -> None:
    before_label = "top = before tight-ID" if source_payload else "top = no preselection"
    after_label = "bottom = after tight-ID" if source_payload else "bottom = tight WP80"
    guide = rounded_box(fig, [0.055, 0.795, 0.250, 0.065], "#F5F8FC", edge=PANEL_EDGE, radius=0.014)
    guide.text(0.50, 0.64, "Read vertically", ha="center", va="center",
               fontsize=15.4, fontweight="bold", color=INK, transform=guide.transAxes)
    guide.text(0.28, 0.28, before_label,
               fontsize=12.8, color=MUTED, ha="center", va="center", transform=guide.transAxes)
    guide.text(0.74, 0.28, after_label,
               fontsize=12.8, color=TAKEAWAY, ha="center", va="center", transform=guide.transAxes)

    ax = rounded_box(fig, [0.323, 0.795, 0.622, 0.065], SOFT_PANEL, edge=PANEL_EDGE, radius=0.014)
    tight_band_label = "green band = after tight-ID" if source_payload else "green band = tight WP80"
    items = (
        ("Signal MC", SAMPLE_COLORS["Signal MC"]["tight"], "-", 0.125),
        ("Inclusive MC", SAMPLE_COLORS["Inclusive MC"]["tight"], "-", 0.375),
        ("shared shape region", OVERLAP, "-", 0.625),
        (tight_band_label, TAKEAWAY, "-", 0.875),
    )
    for label, color, style, x_center in items:
        ax.plot([x_center - 0.050, x_center + 0.050], [0.64, 0.64], color=color, linewidth=3.2,
                linestyle=style, transform=ax.transAxes)
        ax.text(x_center, 0.30, label, ha="center", va="center",
                fontsize=12.8, color=INK, transform=ax.transAxes)


def draw_column_guides(
    fig,
    xs: list[float],
    ys: list[float],
    panel_w: float,
    panel_h: float,
    top: float,
) -> None:
    bottom = ys[-1] - 0.014
    guide_h = top - bottom + 0.021
    for idx, x in enumerate(xs):
        style = PT_COLUMN_STYLES[idx]
        rail = rounded_box(
            fig,
            [x - 0.012, bottom, panel_w + 0.024, guide_h],
            style["face"],
            edge=style["edge"],
            radius=0.014,
            lw=0.9,
        )
        rail.patches[0].set_alpha(0.36)
        rail.plot([0.040, 0.960], [0.982, 0.982], color=style["accent"],
                  linewidth=2.0, alpha=0.72, transform=rail.transAxes)


def draw_panel(
    ax,
    lk: dict[tuple[str, str, str, str], Curve],
    variable: dict[str, object],
    pt_label: str,
    column_style: dict[str, str],
    *,
    source_payload: dict[str, object] | None = None,
) -> None:
    variable_label = str(variable["plain"])
    plotted = []
    for sample in ("Signal MC", "Inclusive MC"):
        for stage_label, _, style in STAGES:
            plotted.append(lk[(sample, variable_label, pt_label, stage_label)])
    ymax = max(max(c.density) if c.density else 0.0 for c in plotted)

    ax.set_facecolor("white")
    for spine in ax.spines.values():
        spine.set_color(column_style["edge"])
        spine.set_linewidth(1.05)

    scale = 0.39 / ymax if ymax > 0 else 1.0
    lanes = (
        ("Before preselection", 1.02, 0.22, 2.45),
        ("Tight WP80", 0.14, 0.15, 2.35),
    )

    ax.axhspan(0.96, 1.46, color="#F5F8FC", zorder=0)
    ax.axhspan(0.08, 0.58, color="#F0FDF7", zorder=0)
    ax.hlines([1.02, 0.14], variable["xlim"][0], variable["xlim"][1],
              color="#D6E0EA", linewidth=0.95, zorder=1)

    for stage_label, base, fill_alpha, lw in lanes:
        inc_curve = lk[("Inclusive MC", variable_label, pt_label, stage_label)]
        sig_curve = lk[("Signal MC", variable_label, pt_label, stage_label)]
        inc_y = base + np.array(inc_curve.density, dtype=float) * scale
        sig_y = base + np.array(sig_curve.density, dtype=float) * scale
        edges = np.array(inc_curve.edges, dtype=float)
        overlap_y = base + np.minimum(
            np.array(inc_curve.density, dtype=float),
            np.array(sig_curve.density, dtype=float),
        ) * scale

        for sample in ("Inclusive MC", "Signal MC"):
            colors = SAMPLE_COLORS[sample]
            curve = lk[(sample, variable_label, pt_label, stage_label)]
            edges = np.array(curve.edges, dtype=float)
            y = base + np.array(curve.density, dtype=float) * scale
            color = colors["before"] if stage_label == "Before preselection" else colors["tight"]
            ax.fill_between(
                edges[:-1], base, y, step="post", color=color,
                alpha=fill_alpha, linewidth=0, zorder=2 if sample == "Inclusive MC" else 3,
            )
            ax.stairs(
                y, edges, color=color, linewidth=lw,
                zorder=4 if sample == "Inclusive MC" else 5,
            )
        ax.fill_between(
            edges[:-1],
            base,
            overlap_y,
            step="post",
            color=OVERLAP,
            alpha=0.16 if stage_label == "Before preselection" else 0.25,
            linewidth=0,
            zorder=3,
        )
        ax.stairs(overlap_y, edges, color=OVERLAP, linewidth=1.0, alpha=0.70, zorder=3.5)

    x_left = variable["xlim"][0] + 0.02 * (variable["xlim"][1] - variable["xlim"][0])
    before_label = "Before tight-ID" if source_payload else "Before (no preselection)"
    after_label = "After tight-ID" if source_payload else "After tight selection WP80"
    ax.text(x_left, 1.45, before_label, ha="left", va="top", fontsize=16.4,
            fontweight="bold", color="#465366",
            bbox=dict(boxstyle="round,pad=0.075", facecolor="white", edgecolor="none", alpha=0.78))
    ax.text(x_left, 0.57, after_label, ha="left", va="top", fontsize=16.4,
            fontweight="bold", color=TAKEAWAY,
            bbox=dict(boxstyle="round,pad=0.075", facecolor="white", edgecolor="none", alpha=0.78))

    ax.set_xlim(*variable["xlim"])
    ax.set_ylim(0.04, 1.52)
    ax.grid(True, color=GRID, linewidth=0.75)
    ax.tick_params(labelsize=9.6, direction="in", top=True, right=True)
    ax.set_yticklabels([])
    ax.tick_params(axis="y", length=0)


def render_slide(
    curves: list[Curve],
    outdir: Path,
    *,
    source_payload: dict[str, object] | None = None,
) -> dict[str, Path]:
    setup_style()
    outdir.mkdir(parents=True, exist_ok=True)
    lk = curve_lookup(curves)

    fig = plt.figure(figsize=(16, 9), dpi=160)
    title = r"0-20% only: tight WP80 pulls background shapes toward signal"
    if source_payload:
        title = r"0-20% only: default tight-ID pulls background shapes toward signal"
    fig.text(0.045, 0.956, title,
             ha="left", va="top", fontsize=31.0, fontweight="bold", color=INK)
    subtitle = (
        r"Shower-shape distributions across $15 \leq p_T < 35$ GeV; top lanes are before selection, "
        r"bottom lanes are after tight WP80 selection."
    )
    if source_payload:
        wp80 = source_payload.get("wp80_formula", {})
        if isinstance(wp80, dict) and "intercept" in wp80 and "slope" in wp80:
            formula = f"T80(c) = {float(wp80['intercept']):.4f} + {float(wp80['slope']):.6f} c"
        elif isinstance(wp80, dict):
            formula = str(wp80.get("expression", "default WP80"))
        else:
            formula = "default WP80"
        subtitle = (
            r"Full weighted signal/inclusive MC sample; top lanes are before tight selection, "
            rf"bottom lanes use default {formula}."
        )
    fig.text(0.046, 0.904, subtitle, ha="left", va="top", fontsize=19.0, color=MUTED)

    draw_legend(fig, source_payload=source_payload)

    left = 0.165
    panel_w = 0.245
    panel_h = 0.166
    hgap = 0.030
    vgap = 0.043
    top = 0.736
    xs = [left + i * (panel_w + hgap) for i in range(3)]
    ys = [top - panel_h - i * (panel_h + vgap) for i in range(3)]

    draw_column_guides(fig, xs, ys, panel_w, panel_h, top)

    for col_idx, (pt_label, _) in enumerate(PT_GROUPS):
        style = PT_COLUMN_STYLES[col_idx]
        header = rounded_box(
            fig,
            [xs[col_idx] + 0.022, top + 0.004, panel_w - 0.044, 0.033],
            style["face"],
            edge=style["edge"],
            radius=0.012,
            lw=1.1,
        )
        header.text(0.50, 0.53, pt_label, ha="center", va="center",
                    fontsize=14.1, fontweight="bold", color=style["accent"],
                    transform=header.transAxes)

    for row_idx, variable in enumerate(VARIABLES):
        label_ax = rounded_box(fig, [0.042, ys[row_idx] + 0.016, 0.098, panel_h - 0.032],
                               "#FFFFFF", edge=PANEL_EDGE, radius=0.016)
        label_ax.text(0.50, 0.62, str(variable["label"]), ha="center", va="center",
                      fontsize=17.2, fontweight="bold", color=INK, transform=label_ax.transAxes)
        label_ax.text(0.50, 0.32, str(variable["feature"]), ha="center", va="center",
                      fontsize=10.4, color=MUTED, transform=label_ax.transAxes)
        for col_idx, (pt_label, _) in enumerate(PT_GROUPS):
            ax = fig.add_axes([xs[col_idx], ys[row_idx], panel_w, panel_h])
            draw_panel(ax, lk, variable, pt_label, PT_COLUMN_STYLES[col_idx], source_payload=source_payload)
            if row_idx == len(VARIABLES) - 1:
                ax.set_xlabel(str(variable["label"]), fontsize=11.8, labelpad=2)
            else:
                ax.tick_params(labelbottom=False)
            if col_idx == 0:
                ax.set_ylabel("normalized shape", fontsize=10.5, labelpad=5)

    note = rounded_box(fig, [0.055, 0.036, 0.890, 0.052], "#FFF3BF", edge="#EAB308", radius=0.014)
    note_text = (
        r"Key point: with no preselection, blue background and red signal are visibly separated; "
        r"after tight WP80, blue contracts into the red signal-like region."
    )
    if source_payload:
        note_text = (
            r"Key point: before the default tight-ID split, inclusive MC is broader; after tight-ID, "
            r"the selected inclusive shape contracts into the signal-like region."
        )
    note.text(0.025, 0.50, note_text,
              ha="left", va="center", fontsize=13.0, color="#7A4B00", transform=note.transAxes)

    png = outdir / "the42_energy_sum_feature_distribution_grid_slide.png"
    manifest = outdir / "the42_energy_sum_feature_distribution_grid_manifest.json"
    script = outdir / "the42_energy_sum_feature_distribution_grid_script.md"
    fig.savefig(png, dpi=160)
    plt.close(fig)
    source_sentence = (
        "Rows are the energy-sum BDT inputs written as THE-42 stage histograms: cluster_et1, E11/E33, and E32/E35."
    )
    if source_payload:
        source_sentence = (
            "Rows are the energy-sum BDT inputs rebuilt from the THE-57 full weighted scored matrix: "
            "cluster_et1, E11/E33, and E32/E35."
        )
    script.write_text(
        "# Energy-Sum BDT Input Distribution Grid Script\n\n"
        "This slide shows the actual shower-shape distributions, not a subtraction or a table. The centrality "
        f"scope is only 0-20 percent. {source_sentence} Columns are broad cluster-pT groups spanning the full "
        "15-35 GeV working range. Each panel has two horizontal stage lanes: the top lane is before tight-ID "
        "and the bottom lane is after the tight centrality-linear WP80 BDT cut. Signal MC is always red and "
        "inclusive MC is always blue.\n\n"
        "The point to emphasize is visual. Before tight-ID, the blue inclusive-MC shape has a broader "
        "background-like shoulder and is visibly separated from the red signal shape. After tight-ID, the blue "
        "shape contracts into the red signal-like region. The subtle purple fill marks the shared shape region, "
        "which is easier to see in the tight lane. The vertical offsets are visual only; each lane is independently "
        "normalized in the shown x range.\n",
        encoding="utf-8",
    )
    return {"png": png, "manifest": manifest, "speaker_script": script}


def write_manifest(
    path: Path,
    curves: list[Curve],
    outputs: dict[str, Path],
    args: argparse.Namespace,
    source_payload: dict[str, object] | None = None,
) -> None:
    payload = {
        "campaign": (
            "THE-57 default 14-feature AuAu baseline WP80 energy-sum overlay"
            if source_payload else
            "THE-42 WP80 centrality-linear AuAu SS overlay"
        ),
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "script": str(THIS_FILE),
        "signal_root": str(args.signal_root.resolve()) if args.signal_root else None,
        "inclusive_root": str(args.inclusive_root.resolve()) if args.inclusive_root else None,
        "curves_json": str(args.curves_json.resolve()) if args.curves_json else None,
        "source_payload_summary": {
            key: source_payload.get(key)
            for key in (
                "schema",
                "matrix",
                "model",
                "full_matrix_rows",
                "rows_loaded",
                "weight_mode",
                "wp80_formula",
                "source_label",
                "model_label",
            )
        } if source_payload else None,
        "plot_mode": (
            "distribution grid, before default tight-ID versus after default tight-ID; "
            "curves unit-normalized in shown x range"
            if source_payload else
            "distribution grid, before preselection versus tight WP80; curves unit-normalized in shown x range"
        ),
        "display_rebin_factor": args.rebin,
        "pt_groups": [{"label": label, "fine_bins": bins} for label, bins in PT_GROUPS],
        "centrality_focus": {"label": CENT_FOCUS_LABEL, "fine_bins": CENT_FOCUS_BINS},
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
        "curve_summaries": [
            {
                key: value
                for key, value in asdict(c).items()
                if key not in {"edges", "density"}
            }
            for c in curves
        ],
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--signal-root", type=Path, default=DEFAULT_SIGNAL_ROOT)
    ap.add_argument("--inclusive-root", type=Path, default=DEFAULT_INCLUSIVE_ROOT)
    ap.add_argument("--curves-json", type=Path, default=None, help="Optional precomputed Curve payload; skips ROOT inputs.")
    ap.add_argument("--output-dir", type=Path, default=DEFAULT_OUTDIR)
    ap.add_argument("--rebin", type=int, default=3, help="Display-only rebin factor applied after summing ROOT histograms.")
    return ap.parse_args()


def main() -> int:
    args = parse_args()
    if args.rebin < 1:
        raise SystemExit("--rebin must be >= 1")
    source_payload = None
    if args.curves_json:
        if not args.curves_json.exists():
            raise SystemExit(f"Missing curves JSON input: {args.curves_json}")
        curves, source_payload = load_curves_json(args.curves_json)
        args.signal_root = None
        args.inclusive_root = None
    else:
        for path in (args.signal_root, args.inclusive_root):
            if not path.exists():
                raise SystemExit(f"Missing ROOT input: {path}")
        curves = collect_curves(args.signal_root, args.inclusive_root, args.rebin)
    outputs = render_slide(curves, args.output_dir, source_payload=source_payload)
    write_manifest(outputs["manifest"], curves, outputs, args, source_payload=source_payload)
    for key, value in outputs.items():
        print(f"{key}={value}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
