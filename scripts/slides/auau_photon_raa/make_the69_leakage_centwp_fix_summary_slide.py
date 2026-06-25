#!/usr/bin/env python3
"""Build the THE-69 corrected AuAu ABCD leakage summary slide."""

from __future__ import annotations

import csv
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch, Polygon


THIS_FILE = Path(__file__).resolve()
REPO = next((p for p in THIS_FILE.parents if (p / "AGENTS.md").exists()), THIS_FILE.parents[4])
SCRIPTS = REPO / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.append(str(SCRIPTS))

from slides.common.slide_defaults import SLIDE_DPI, slide_figsize  # noqa: E402


CORRECTED_AUAU_ROOT = (
    REPO
    / "InputFiles/the69_leakageCentWP_fix/RecoilJets_embeddedPhoton12plus20_MERGED.root"
)
PP_REFERENCE_ROOT = (
    REPO
    / "dataOutput/ppg12PhotonYield/THE76_ppg12_photon_yield_v1_sim_20260616"
    / "merged_roots/jetMinPtScan_dphiScan_vz60_isoR40_isSliding_preselectionNewPPG12_tightNewPPG12_nonTightNewPPG12"
    / "photonJet5and10and20merged_SIM/RecoilJets_photonjet5plus10plus20_MERGED.root"
)
OUT_DIR = REPO / "dataOutput/auauPhysicsQA/THE69_leakageCentWP_fix_20260620"
OUT_PNG = OUT_DIR / "the69_leakage_centwp_fix_summary_slide.png"
OUT_CSV = OUT_DIR / "the69_leakage_centwp_fix_points.csv"
OUT_JSON = OUT_DIR / "the69_leakage_centwp_fix_manifest.json"
OUT_SCRIPT = OUT_DIR / "the69_leakage_centwp_fix_speaker_script.md"

TOPDIR = "SIM"
ISO_TAG = "isoR40_isSliding"
PT_BINS = [(16, 18), (18, 20), (20, 22), (22, 24), (24, 26), (26, 35)]
CENTRALITIES = [
    ("0_20", "0-20%", "#111827", "o", -0.18),
    ("20_50", "20-50%", "#2563EB", "s", 0.00),
    ("50_80", "50-80%", "#059669", "D", 0.18),
]
REGIONS = [
    ("B", 2, r"$f_B = B_{\rm sig}/A_{\rm sig}$", "tight + non-isolated", (-0.006, 0.145)),
    ("C", 3, r"$f_C = C_{\rm sig}/A_{\rm sig}$", "isolated + non-tight", (0.0, 0.430)),
    ("D", 4, r"$f_D = D_{\rm sig}/A_{\rm sig}$", "non-isolated + non-tight", (-0.0025, 0.050)),
]
PP_SOURCE = "pp_ppg12_sliding_reference"
PPG12_SIGNAL_HISTS = {
    "A": "h_tight_iso_cluster_signal_0",
    "B": "h_tight_noniso_cluster_signal_0",
    "C": "h_nontight_iso_cluster_signal_0",
    "D": "h_nontight_noniso_cluster_signal_0",
}


@dataclass(frozen=True)
class Count:
    value: float
    error: float


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "axes.facecolor": "white",
            "axes.edgecolor": "#1F2937",
            "axes.linewidth": 1.0,
            "axes.labelsize": 15.0,
            "axes.titlesize": 19.0,
            "xtick.labelsize": 12.2,
            "ytick.labelsize": 12.2,
            "legend.fontsize": 12.0,
            "mathtext.fontset": "dejavuserif",
        }
    )


def add_arrow_bullet(
    fig,
    x: float,
    y: float,
    text: str,
    *,
    fontsize: float = 13.0,
    text_color: str = "#334155",
    arrow_color: str = "#0F172A",
) -> None:
    """Draw a slide-style arrowhead bullet without relying on font glyphs."""
    triangle = Polygon(
        [[x, y - 0.0046], [x, y + 0.0046], [x + 0.0075, y]],
        closed=True,
        transform=fig.transFigure,
        facecolor=arrow_color,
        edgecolor=arrow_color,
        linewidth=0,
        clip_on=False,
    )
    fig.patches.append(triangle)
    fig.text(
        x + 0.0125,
        y,
        text,
        fontsize=fontsize,
        color=text_color,
        va="center",
    )


def open_root(path: Path):
    import ROOT

    ROOT.gROOT.SetBatch(True)
    ROOT.TH1.AddDirectory(False)
    root_file = ROOT.TFile.Open(str(path), "READ")
    if not root_file or root_file.IsZombie():
        raise OSError(f"Could not open ROOT file: {path}")
    return root_file


def auau_hist_path(lo: int, hi: int, cent_key: str) -> str:
    return f"{TOPDIR}/h_sigABCD_MC_{ISO_TAG}_pT_{lo}_{hi}_cent_{cent_key}"


def pp_hist_path(region: str) -> str:
    return f"{TOPDIR}/{PPG12_SIGNAL_HISTS[region]}"


def bin_count(hist, bin_index: int) -> Count:
    value = float(hist.GetBinContent(bin_index))
    error = float(hist.GetBinError(bin_index))
    if error <= 0.0 and value > 0.0:
        error = math.sqrt(value)
    return Count(value, error)


def integrated_count(hist, lo: int, hi: int) -> Count:
    """Integrate a PPG12 photon-yield histogram into the AuAu comparison bin."""
    value = 0.0
    err2 = 0.0
    for bin_index in range(1, hist.GetNbinsX() + 1):
        center = float(hist.GetXaxis().GetBinCenter(bin_index))
        if lo <= center < hi:
            value += float(hist.GetBinContent(bin_index))
            err = float(hist.GetBinError(bin_index))
            if err <= 0.0 and hist.GetBinContent(bin_index) > 0.0:
                err = math.sqrt(float(hist.GetBinContent(bin_index)))
            err2 += err * err
    return Count(value, math.sqrt(err2))


def ratio_with_error(num: Count, den: Count) -> tuple[float, float]:
    if den.value <= 0.0:
        return float("nan"), float("nan")
    value = num.value / den.value
    if num.value > 0.0:
        rel2 = (num.error / num.value) ** 2 + (den.error / den.value) ** 2
        error = abs(value) * math.sqrt(max(0.0, rel2))
    else:
        error = num.error / den.value
    return value, error


def collect_auau_points(root_path: Path) -> list[dict]:
    root_file = open_root(root_path)
    rows: list[dict] = []
    try:
        for cent_key, cent_label, _, _, _ in CENTRALITIES:
            for lo, hi in PT_BINS:
                path = auau_hist_path(lo, hi, cent_key)
                hist = root_file.Get(path)
                if not hist:
                    raise KeyError(f"Missing histogram: {path}")
                a_count = bin_count(hist, 1)
                for region, bin_index, _, _, _ in REGIONS:
                    side_count = bin_count(hist, bin_index)
                    ratio, ratio_err = ratio_with_error(side_count, a_count)
                    rows.append(
                        {
                            "source": "corrected_auau",
                            "region": region,
                            "centrality": cent_label,
                            "cent_key": cent_key,
                            "pt_lo": lo,
                            "pt_hi": hi,
                            "pt_mid": 0.5 * (lo + hi),
                            "a_signal": a_count.value,
                            "a_signal_err": a_count.error,
                            "sideband_signal": side_count.value,
                            "sideband_signal_err": side_count.error,
                            "leakage_fraction": ratio,
                            "leakage_fraction_err": ratio_err,
                            "histogram": path,
                        }
                    )
    finally:
        root_file.Close()
    return rows


def collect_pp_points(root_path: Path) -> list[dict]:
    root_file = open_root(root_path)
    rows: list[dict] = []
    try:
        a_hist = root_file.Get(pp_hist_path("A"))
        if not a_hist:
            raise KeyError(f"Missing histogram: {pp_hist_path('A')}")
        region_hists = {}
        for region, *_ in REGIONS:
            path = pp_hist_path(region)
            hist = root_file.Get(path)
            if not hist:
                raise KeyError(f"Missing histogram: {path}")
            region_hists[region] = hist
        for lo, hi in PT_BINS:
            a_count = integrated_count(a_hist, lo, hi)
            for region, _, _, _, _ in REGIONS:
                side_count = integrated_count(region_hists[region], lo, hi)
                ratio, ratio_err = ratio_with_error(side_count, a_count)
                rows.append(
                    {
                        "source": PP_SOURCE,
                        "region": region,
                        "centrality": "pp",
                        "cent_key": "pp",
                        "pt_lo": lo,
                        "pt_hi": hi,
                        "pt_mid": 0.5 * (lo + hi),
                        "a_signal": a_count.value,
                        "a_signal_err": a_count.error,
                        "sideband_signal": side_count.value,
                        "sideband_signal_err": side_count.error,
                        "leakage_fraction": ratio,
                        "leakage_fraction_err": ratio_err,
                        "histogram": f"{pp_hist_path('A')} :: {pp_hist_path(region)} integrated {lo}-{hi}",
                    }
                )
    finally:
        root_file.Close()
    return rows


def weighted_mean(rows: list[dict], source: str, region: str, cent_key: str) -> float:
    subset = [
        r
        for r in rows
        if r["source"] == source and r["region"] == region and r["cent_key"] == cent_key
    ]
    den = sum(float(r["a_signal"]) for r in subset)
    num = sum(float(r["sideband_signal"]) for r in subset)
    return num / den if den > 0.0 else float("nan")


def format_leakage(value: float) -> str:
    if not math.isfinite(value):
        return "n/a"
    if 0.0 < abs(value) < 0.001:
        return f"{value:.1e}"
    return f"{value:.3f}"


def write_points(rows: list[dict]) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "source",
        "region",
        "centrality",
        "cent_key",
        "pt_lo",
        "pt_hi",
        "pt_mid",
        "a_signal",
        "a_signal_err",
        "sideband_signal",
        "sideband_signal_err",
        "leakage_fraction",
        "leakage_fraction_err",
        "histogram",
    ]
    with OUT_CSV.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def plot_panel(
    ax,
    rows: list[dict],
    region: str,
    title: str,
    subtitle: str,
    ylim: tuple[float, float],
    show_ylabel: bool,
) -> None:
    for cent_key, cent_label, color, marker, offset in CENTRALITIES:
        subset = [
            r
            for r in rows
            if r["source"] == "corrected_auau" and r["region"] == region and r["cent_key"] == cent_key
        ]
        x = np.array([r["pt_mid"] + offset for r in subset], dtype=float)
        y = np.array([r["leakage_fraction"] for r in subset], dtype=float)
        yerr = np.array([r["leakage_fraction_err"] for r in subset], dtype=float)
        ax.errorbar(
            x,
            y,
            yerr=yerr,
            fmt=marker,
            linestyle="none",
            markersize=8.0,
            capsize=2.8,
            elinewidth=1.15,
            markerfacecolor=color,
            markeredgecolor=color,
            markeredgewidth=1.25,
            color=color,
            alpha=0.96,
            label=cent_label,
        )

    pp_rows = [
        r for r in rows if r["source"] == PP_SOURCE and r["region"] == region
    ]
    x = np.array([r["pt_mid"] + 0.36 for r in pp_rows], dtype=float)
    y = np.array([r["leakage_fraction"] for r in pp_rows], dtype=float)
    yerr = np.array([r["leakage_fraction_err"] for r in pp_rows], dtype=float)
    ax.errorbar(
        x,
        y,
        yerr=yerr,
        fmt="o",
        linestyle="none",
        markersize=8.4,
        capsize=2.8,
        elinewidth=1.15,
        markerfacecolor="white",
        markeredgecolor="#DC2626",
        markeredgewidth=2.0,
        color="#DC2626",
        alpha=0.98,
        label="pp PPG12 baseline",
    )

    ax.set_title(title, fontweight="bold", pad=8)
    ax.text(
        0.025,
        0.930,
        subtitle,
        transform=ax.transAxes,
        fontsize=12.4,
        fontweight="bold",
        color="#475569",
        va="top",
    )
    ax.set_xlabel(r"reco cluster $E_T$ [GeV]", labelpad=7)
    ax.set_ylabel("truth-signal leakage fraction" if show_ylabel else "", labelpad=8)
    ax.set_xlim(15.2, 35.8)
    ax.set_ylim(*ylim)
    ax.set_xticks([16, 18, 20, 22, 24, 26, 30, 35])
    ax.set_xticklabels(["16", "18", "20", "22", "24", "26", "30", "35"])
    ax.tick_params(axis="x", pad=3)
    ax.tick_params(axis="y", pad=4)
    ax.grid(True, which="major", color="#DDE3EA", linewidth=0.8, alpha=0.85)


def add_summary_table(fig, rows: list[dict]) -> None:
    ax = fig.add_axes([0.075, 0.062, 0.87, 0.230])
    ax.set_axis_off()
    box = FancyBboxPatch(
        (0.0, 0.0),
        1.0,
        1.0,
        boxstyle="round,pad=0.010,rounding_size=0.015",
        facecolor="#F8FAFC",
        edgecolor="#CBD5E1",
        linewidth=1.0,
        transform=ax.transAxes,
    )
    ax.add_patch(box)
    ax.text(
        0.018,
        0.885,
        r"A-weighted mean leakage over 16-35 GeV",
        fontsize=15.2,
        fontweight="bold",
        color="#111827",
        va="center",
        transform=ax.transAxes,
    )
    ax.text(
        0.018,
        0.760,
        "Rows also define the plot marker key. B/D are non-isolated leakage checks; C is isolated, non-tight signal.",
        fontsize=12.2,
        color="#475569",
        va="center",
        transform=ax.transAxes,
    )
    headers = [
        "marker / sample",
        "B: tight, non-isolated",
        "C: isolated, non-tight",
        "D: non-isolated, non-tight",
    ]
    col_x = [0.030, 0.285, 0.555, 0.815]
    ax.add_patch(
        FancyBboxPatch(
            (0.014, 0.565),
            0.972,
            0.115,
            boxstyle="round,pad=0.004,rounding_size=0.010",
            facecolor="#EFF6FF",
            edgecolor="none",
            transform=ax.transAxes,
            zorder=0,
        )
    )
    for x, header in zip(col_x, headers):
        ax.text(
            x,
            0.622,
            header,
            fontsize=12.4,
            fontweight="bold",
            color="#111827",
            va="center",
            transform=ax.transAxes,
        )
    y0 = 0.440
    for idx, (cent_key, cent_label, color, _, _) in enumerate(CENTRALITIES):
        y = y0 - idx * 0.122
        marker = CENTRALITIES[idx][3]
        ax.plot(
            [col_x[0] + 0.012],
            [y],
            marker=marker,
            markersize=7.6,
            markerfacecolor=color,
            markeredgecolor=color,
            linestyle="none",
            transform=ax.transAxes,
            clip_on=False,
        )
        ax.text(
            col_x[0] + 0.032,
            y,
            f"Au+Au {cent_label}",
            fontsize=12.8,
            fontweight="bold",
            color=color,
            va="center",
            transform=ax.transAxes,
        )
        for col_idx, (region, _, _, _, _) in enumerate(REGIONS, start=1):
            value = weighted_mean(rows, "corrected_auau", region, cent_key)
            text = format_leakage(value)
            ax.text(
                col_x[col_idx],
                y,
                text,
                fontsize=12.6,
                color="#334155",
                va="center",
                transform=ax.transAxes,
            )
    pp_y = y0 - len(CENTRALITIES) * 0.122 - 0.020
    ax.plot(
        [0.012, 0.988],
        [pp_y + 0.067, pp_y + 0.067],
        color="#E2E8F0",
        linewidth=1.0,
        transform=ax.transAxes,
        clip_on=False,
    )
    ax.plot(
        [col_x[0] + 0.012],
        [pp_y],
        marker="o",
        markersize=8.0,
        markerfacecolor="white",
        markeredgecolor="#DC2626",
        markeredgewidth=1.8,
        linestyle="none",
        transform=ax.transAxes,
        clip_on=False,
    )
    ax.text(
        col_x[0] + 0.032,
        pp_y,
        "pp PPG12",
        fontsize=12.8,
        fontweight="bold",
        color="#DC2626",
        va="center",
        transform=ax.transAxes,
    )
    for col_idx, (region, _, _, _, _) in enumerate(REGIONS, start=1):
        value = weighted_mean(rows, PP_SOURCE, region, "pp")
        ax.text(
            col_x[col_idx],
            pp_y,
            format_leakage(value),
            fontsize=12.6,
            color="#334155",
            va="center",
            transform=ax.transAxes,
        )


def make_slide(rows: list[dict]) -> None:
    setup_style()
    fig = plt.figure(figsize=slide_figsize(), dpi=SLIDE_DPI)
    fig.text(
        0.045,
        0.945,
        "Corrected Au+Au leakage compared with exact PPG12 pp",
        fontsize=25.0,
        fontweight="bold",
        color="#111827",
        va="top",
    )
    fig.text(
        0.045,
        0.900,
        (
        "Leakage is truth-signal sideband content normalized to A = tight + isolated signal."
        ),
        fontsize=14.0,
        color="#334155",
        va="top",
    )
    add_arrow_bullet(
        fig,
        0.045,
        0.861,
        (
            r"Red open circles: exact PPG12 pp sliding baseline "
            r"($R=0.4$, $E_T^{iso}<0.490+0.037E_T$, gap=0.8)."
        ),
        fontsize=14.4,
    )
    add_arrow_bullet(
        fig,
        0.045,
        0.832,
        r"Au+Au: corrected centrality-dependent WP80 with current sideband setting (sideGap=0); B/D are the non-isolated leakage checks.",
        fontsize=14.4,
    )

    axes = [
        fig.add_axes([0.080, 0.330, 0.265, 0.355]),
        fig.add_axes([0.382, 0.330, 0.265, 0.355]),
        fig.add_axes([0.684, 0.330, 0.265, 0.355]),
    ]
    for idx, (ax, (region, _, title, subtitle, ylim)) in enumerate(zip(axes, REGIONS)):
        plot_panel(ax, rows, region, title, subtitle, ylim, show_ylabel=(idx == 0))

    add_summary_table(fig, rows)
    fig.savefig(OUT_PNG, dpi=SLIDE_DPI)
    plt.close(fig)


def write_manifest(rows: list[dict]) -> None:
    summary: dict[str, dict[str, dict[str, float]]] = {}
    for cent_key, cent_label, *_ in CENTRALITIES:
        summary[cent_label] = {}
        for region, *_ in REGIONS:
            after = weighted_mean(rows, "corrected_auau", region, cent_key)
            summary[cent_label][region] = {
                "corrected": after,
            }
    summary["pp PPG12 baseline"] = {
        region: {"reference": weighted_mean(rows, PP_SOURCE, region, "pp")}
        for region, *_ in REGIONS
    }
    manifest = {
        "campaign": "THE-69 leakage-centWP signal-MC proof",
        "corrected_campaign_tag": "the69_leakageCentWP_fix_signalMC_20260620",
        "corrected_input_root": str(CORRECTED_AUAU_ROOT),
        "pp_reference_root": str(PP_REFERENCE_ROOT),
        "output_png": str(OUT_PNG),
        "points_csv": str(OUT_CSV),
        "model_id": "centAsFeatBase3x3_pt15to35",
        "wp80": {
            "row": "auauCentInputBase3x3BDT|centlinear|0.53471108|0.0012284143|15|35|1",
            "mode": "centlinear",
            "formula": "T80(c) = 0.53471108 + 0.0012284143*c",
        },
        "isolation": {
            "view": "isoR40_isSliding",
            "formula": "Eiso(R=0.4) < 7.57 - 0.0658*c",
            "sideGap": 0,
        },
        "pp_reference_isolation": {
            "view": "isoR40_isSliding",
            "formula": "Eiso(R=0.4) < 0.490 + 0.037*E_T",
            "sideband_gap_GeV": 0.8,
            "source_note": "Exact current PPG12 photon-yield pp baseline output, not the old fixed-isolation reference and not the table-QA proxy h_sigABCD_MC_isoR40_isSliding family.",
        },
        "leakage_definition": "f_X = N_sig(X) / N_sig(A), with A=tight+isolated truth signal",
        "shown_et_bins_GeV": PT_BINS,
        "excluded_boundary_note": "14-16 GeV straddles the 15 GeV BDT/WP80 lower edge and is excluded from this summary slide.",
        "summary_weighted_means_16_35_GeV": summary,
        "interpretation": (
            "This audience-facing slide shows only the corrected AuAu output and the exact PPG12 pp sliding-isolation baseline overlay. "
            "B and D are small but nonzero in AuAu with sideGap=0; pp B and D are also nonzero but much smaller and are printed in scientific notation in the table."
        ),
    }
    OUT_JSON.write_text(json.dumps(manifest, indent=2) + "\n")


def write_speaker_script() -> None:
    OUT_SCRIPT.write_text(
        """# THE-69 corrected leakage comparison

This slide is a focused signal-leakage check. A is tight and isolated truth signal, and each panel shows sideband truth signal normalized to A.

The red open circles are the exact current PPG12 pp sliding-isolation baseline, not a fixed-isolation or table-QA proxy. The Au+Au points use the corrected centrality-dependent WP80 and the current sideband setting, sideGap=0. The bottom table doubles as the marker key and gives the A-weighted mean leakage over 16-35 GeV.

B and D are the non-isolated leakage checks. pp is near zero on this scale, while Au+Au remains visibly finite. C is isolated signal that fails the tight ID, so it is a working-point and composition check rather than an isolation-sideband failure.

The audience takeaway is that the corrected centrality-dependent BDT cut is now applied, the pp comparison is the exact PPG12 baseline, and any remaining non-isolated leakage question is now isolated to the Au+Au sideband definition rather than a stale WP80 bug.
""",
        encoding="utf-8",
    )


def main() -> None:
    rows = collect_auau_points(CORRECTED_AUAU_ROOT) + collect_pp_points(PP_REFERENCE_ROOT)
    write_points(rows)
    make_slide(rows)
    write_manifest(rows)
    write_speaker_script()
    print(f"wrote {OUT_PNG}")
    print(f"wrote {OUT_JSON}")
    print(f"wrote {OUT_SCRIPT}")


if __name__ == "__main__":
    main()
