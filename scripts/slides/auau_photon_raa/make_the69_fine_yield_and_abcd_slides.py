#!/usr/bin/env python3
"""Build fine-ET tight-yield and ABCD-purity QA slides for THE-69."""

from __future__ import annotations

import csv
import json
import math
import sys
import textwrap
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


THIS_FILE = Path(__file__).resolve()
REPO = next((p for p in THIS_FILE.parents if (p / "AGENTS.md").exists()), THIS_FILE.parents[4])
SCRIPTS = REPO / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.append(str(SCRIPTS))

from slides.common.slide_defaults import SLIDE_DPI, slide_figsize  # noqa: E402


INPUT_DIR = REPO / "InputFiles/the69_default_auau_physicsqa"
DATA_ROOT = INPUT_DIR / "RecoilJets_auau_ALL_preselectionNewPPG12_tightAuAuCentInputBase3x3BDT_nonTightAuAuBDTComplement_baseVariant.root"
SIGNAL_ROOT = INPUT_DIR / "RecoilJets_embeddedPhoton12plus20_MERGED.root"
INCLUSIVE_ROOT = INPUT_DIR / "RecoilJets_embeddedJet12plus20plus30plus40_MERGED.root"
OUT_DIR = REPO / "dataOutput/auauPhysicsQA/THE69_defaultAuAuBDT_physicsQA_20260620/fine_yield_abcd_companion_slides"

DATA_TOP = "MBD_NS_geq_2_vtx_lt_150"
SIM_TOP = "SIM"
CENT = [
    ("0_20", "0-20%", 10.0, "#1F77B4"),
    ("20_50", "20-50%", 35.0, "#2CA02C"),
    ("50_80", "50-80%", 65.0, "#7A5AA6"),
]
PT_BINS = [(14, 16), (16, 18), (18, 20), (20, 22), (22, 24), (24, 26), (26, 35)]
SAMPLES = [
    ("Data", DATA_ROOT, DATA_TOP, "#111827", "o", 0.0),
    ("Signal MC", SIGNAL_ROOT, SIM_TOP, "#C43C39", "s", -0.18),
    ("Inclusive MC", INCLUSIVE_ROOT, SIM_TOP, "#2F63C6", "^", 0.18),
]
TRUTH_SIGNAL_SAMPLE = ("Truth-tagged signal MC", SIGNAL_ROOT, SIM_TOP, "#059669", "D", 0.12)
REGIONS = [
    ("A", "#111827", "o"),
    ("B", "#B45309", "s"),
    ("C", "#7C3AED", "^"),
    ("D", "#0F766E", "D"),
]


@dataclass
class Point:
    value: float
    error: float
    source: str


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
            "axes.labelsize": 11,
            "axes.titlesize": 14,
            "xtick.labelsize": 9.5,
            "ytick.labelsize": 9.5,
            "legend.fontsize": 9.2,
        }
    )


def open_root(path: Path):
    import ROOT

    ROOT.gROOT.SetBatch(True)
    ROOT.TH1.AddDirectory(False)
    f = ROOT.TFile.Open(str(path), "READ")
    if not f or f.IsZombie():
        raise OSError(f"Could not open ROOT file: {path}")
    return f


def integral_and_error_from_file(root_file, obj_path: str) -> Point:
    h = root_file.Get(obj_path)
    if not h or not hasattr(h, "GetNbinsX"):
        return Point(float("nan"), float("nan"), obj_path)
    value = 0.0
    err2 = 0.0
    for i in range(1, h.GetNbinsX() + 1):
        value += float(h.GetBinContent(i))
        err2 += float(h.GetBinError(i)) ** 2
    return Point(value, math.sqrt(err2), obj_path)


def tight_iso_path(top: str, lo: int, hi: int, cent: str) -> str:
    return f"{top}/h_Eiso_tight_isoR40_pT_{lo}_{hi}_cent_{cent}"


def truth_signal_tight_iso_path(top: str, lo: int, hi: int, cent: str) -> str:
    return f"{top}/h_EisoReco_truthSigMatched_tight_isoR40_pT_{lo}_{hi}_cent_{cent}"


def abcd_path(top: str, region: str, lo: int, hi: int, cent: str) -> str:
    return f"{top}/h_Eiso_ABCD_{region}_isoR40_isSliding_pT_{lo}_{hi}_cent_{cent}"


def purity_from_abcd(a: Point, b: Point, c: Point, d: Point) -> tuple[float, float]:
    if not all(np.isfinite([a.value, b.value, c.value, d.value])) or a.value <= 0 or d.value <= 0:
        return float("nan"), float("nan")
    q = b.value * c.value / d.value
    raw = (a.value - q) / a.value
    purity = max(0.0, min(1.0, raw))
    var_q = 0.0
    if b.value > 0:
        var_q += (q * b.error / b.value) ** 2
    if c.value > 0:
        var_q += (q * c.error / c.value) ** 2
    if d.value > 0:
        var_q += (q * d.error / d.value) ** 2
    dp_da = q / (a.value * a.value)
    dp_dq = -1.0 / a.value
    var = (dp_da * a.error) ** 2 + (dp_dq**2) * var_q
    return purity, math.sqrt(max(0.0, var))


def add_title(fig, title: str, subtitle: str) -> None:
    fig.text(0.045, 0.945, title, ha="left", va="top", fontsize=24, fontweight="bold", color="#111827")
    fig.text(0.045, 0.895, subtitle, ha="left", va="top", fontsize=13.5, color="#374151")


def add_takeaway(fig, text: str) -> None:
    fig.text(
        0.055,
        0.047,
        textwrap.fill(text, 132),
        ha="left",
        va="bottom",
        fontsize=13.2,
        color="#111827",
        bbox=dict(facecolor="#F8FAFC", edgecolor="#CBD5E1", linewidth=1.1, boxstyle="round,pad=0.45"),
    )


def collect_points() -> tuple[list[dict], list[dict]]:
    yield_rows: list[dict] = []
    abcd_rows: list[dict] = []
    root_files = {
        str(DATA_ROOT): open_root(DATA_ROOT),
        str(SIGNAL_ROOT): open_root(SIGNAL_ROOT),
        str(INCLUSIVE_ROOT): open_root(INCLUSIVE_ROOT),
    }
    try:
        for cent, cent_label, _, _ in CENT:
            for lo, hi in PT_BINS:
                width = hi - lo
                mid = 0.5 * (lo + hi)
                for sample, root, top, _, _, _ in SAMPLES:
                    p = integral_and_error_from_file(root_files[str(root)], tight_iso_path(top, lo, hi, cent))
                    yield_rows.append(
                        {
                            "sample": sample,
                            "centrality": cent_label,
                            "cent_key": cent,
                            "pt_lo": lo,
                            "pt_hi": hi,
                            "pt_mid": mid,
                            "width": width,
                            "yield_per_gev": p.value / width if np.isfinite(p.value) else float("nan"),
                            "yield_per_gev_err": p.error / width if np.isfinite(p.error) else float("nan"),
                            "raw_integral": p.value,
                            "raw_error": p.error,
                            "source": p.source,
                        }
                    )
                sample, root, top, _, _, _ = TRUTH_SIGNAL_SAMPLE
                p = integral_and_error_from_file(root_files[str(root)], truth_signal_tight_iso_path(top, lo, hi, cent))
                yield_rows.append(
                    {
                        "sample": sample,
                        "centrality": cent_label,
                        "cent_key": cent,
                        "pt_lo": lo,
                        "pt_hi": hi,
                        "pt_mid": mid,
                        "width": width,
                        "yield_per_gev": p.value / width if np.isfinite(p.value) else float("nan"),
                        "yield_per_gev_err": p.error / width if np.isfinite(p.error) else float("nan"),
                        "raw_integral": p.value,
                        "raw_error": p.error,
                        "source": p.source,
                    }
                )
                regions = {
                    region: integral_and_error_from_file(root_files[str(DATA_ROOT)], abcd_path(DATA_TOP, region, lo, hi, cent))
                    for region, _, _ in REGIONS
                }
                purity, purity_err = purity_from_abcd(regions["A"], regions["B"], regions["C"], regions["D"])
                row = {
                    "centrality": cent_label,
                    "cent_key": cent,
                    "pt_lo": lo,
                    "pt_hi": hi,
                    "pt_mid": mid,
                    "width": width,
                    "raw_purity": purity,
                    "raw_purity_err": purity_err,
                }
                for region, _, _ in REGIONS:
                    p = regions[region]
                    row[f"{region}"] = p.value
                    row[f"{region}_err"] = p.error
                    row[f"{region}_per_gev"] = p.value / width if np.isfinite(p.value) else float("nan")
                    row[f"{region}_source"] = p.source
                abcd_rows.append(row)
    finally:
        for f in root_files.values():
            f.Close()
    return yield_rows, abcd_rows


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        return
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_text_artifact(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(textwrap.dedent(text).strip() + "\n")
    return path


def make_yield_slide(rows: list[dict]) -> Path:
    fig, axes = plt.subplots(1, 3, figsize=slide_figsize(), dpi=SLIDE_DPI)
    fig.subplots_adjust(left=0.070, right=0.965, top=0.755, bottom=0.190, wspace=0.235)
    add_title(
        fig,
        "Tight-ID yield-shape QA: data compared with signal and inclusive MC",
        "Fine cluster-ET bins, isolated tight WP80 candidates. Each sample is area-normalized within a centrality panel; this is shape QA, not MC/data closure.",
    )
    for ax, (cent, cent_label, _, _) in zip(axes, CENT):
        cent_rows = [r for r in rows if r["cent_key"] == cent]
        for sample, _, _, color, marker, offset in SAMPLES:
            pts = [r for r in cent_rows if r["sample"] == sample]
            total = sum(float(r["yield_per_gev"]) * float(r["width"]) for r in pts if np.isfinite(float(r["yield_per_gev"])))
            x = np.array([float(r["pt_mid"]) + offset for r in pts])
            y = np.array([float(r["yield_per_gev"]) / total if total > 0 else float("nan") for r in pts])
            e = np.array([float(r["yield_per_gev_err"]) / total if total > 0 else float("nan") for r in pts])
            ax.errorbar(x, y, yerr=e, fmt=marker, color=color, markersize=6.0, capsize=2.2, linewidth=1.05, label=sample)
        ax.set_title(f"AuAu {cent_label}")
        ax.set_xlabel("cluster ET [GeV]")
        ax.set_yscale("log")
        ax.grid(True, axis="y", which="both", alpha=0.25)
        ax.set_xlim(13.2, 35.8)
        ax.set_ylim(1.0e-4, 3.0e-1)
        if ax is axes[0]:
            ax.set_ylabel("area-normalized tight yield / GeV")
        else:
            ax.set_yticklabels([])
    handles = [
        plt.Line2D([0], [0], color=color, marker=marker, lw=0, markersize=7.0, label=sample)
        for sample, _, _, color, marker, _ in SAMPLES
    ]
    leg = fig.legend(
        handles=handles,
        loc="upper right",
        bbox_to_anchor=(0.935, 0.855),
        ncol=3,
        frameon=True,
        fancybox=False,
        edgecolor="#CBD5E1",
        facecolor="white",
        framealpha=0.96,
        borderpad=0.55,
        handletextpad=0.55,
        columnspacing=1.25,
        title="Shape-normalized samples",
    )
    leg.get_title().set_fontweight("bold")
    add_takeaway(
        fig,
        "Use this as a sanity check on the tight-candidate ET shape only. The absolute data/MC rates are intentionally not compared here; the next normalization step is exposure/luminosity and efficiency.",
    )
    out = OUT_DIR / "the69_tight_yield_shape_data_signal_inclusive_mc.png"
    fig.savefig(out, dpi=SLIDE_DPI, facecolor="white")
    plt.close(fig)
    return out


def make_abcd_slide(rows: list[dict]) -> Path:
    fig, axes = plt.subplots(2, 3, figsize=slide_figsize(), dpi=SLIDE_DPI)
    fig.subplots_adjust(left=0.070, right=0.965, top=0.705, bottom=0.185, wspace=0.235, hspace=0.340)
    add_title(
        fig,
        "Raw ABCD purity companion: same fine ET bins as the yield-shape check",
        "Data-only ABCD sideband estimate using isoR40 sliding isolation. Purity = max(0, A - BC/D) / A; this is a diagnostic before final purity systematics.",
    )
    for col, (cent, cent_label, _, _) in enumerate(CENT):
        pts = [r for r in rows if r["cent_key"] == cent]
        ax = axes[0][col]
        x = np.array([float(r["pt_mid"]) for r in pts])
        p = np.array([float(r["raw_purity"]) for r in pts])
        pe = np.array([float(r["raw_purity_err"]) for r in pts])
        ax.errorbar(x, p, yerr=pe, fmt="o", color="#111827", markersize=6.0, capsize=2.2, linewidth=1.05)
        ax.set_title(f"AuAu {cent_label}")
        ax.set_ylim(-0.05, 1.10)
        ax.set_xlim(13.2, 35.8)
        ax.grid(True, axis="y", alpha=0.28)
        if col == 0:
            ax.set_ylabel("raw ABCD purity")
        else:
            ax.set_yticklabels([])
        ax.tick_params(labelbottom=False)

        axb = axes[1][col]
        for region, color, marker in REGIONS:
            y = np.array([float(r[f"{region}_per_gev"]) for r in pts])
            ye = np.array([float(r[f"{region}_err"]) / float(r["width"]) for r in pts])
            axb.errorbar(x, y, yerr=ye, fmt=marker, color=color, markersize=4.9, capsize=1.8, linewidth=0.9, label=region)
        axb.set_yscale("log")
        axb.set_xlim(13.2, 35.8)
        axb.set_xlabel("cluster ET [GeV]")
        axb.grid(True, axis="y", which="both", alpha=0.24)
        if col == 0:
            axb.set_ylabel("ABCD region counts / GeV")
        else:
            axb.set_yticklabels([])
    handles = [
        plt.Line2D([0], [0], color=color, marker=marker, lw=0, markersize=6.0, label=region)
        for region, color, marker in REGIONS
    ]
    leg = fig.legend(
        handles=handles,
        loc="upper right",
        bbox_to_anchor=(0.935, 0.832),
        ncol=4,
        frameon=True,
        fancybox=False,
        edgecolor="#CBD5E1",
        facecolor="white",
        framealpha=0.96,
        borderpad=0.50,
        handletextpad=0.45,
        columnspacing=0.95,
        title="ABCD regions",
    )
    leg.get_title().set_fontweight("bold")
    add_takeaway(
        fig,
        "This catches composition drift that a yield-only plot can hide. The peripheral and high-ET bins are visibly sparse, so these bins need careful treatment before a final corrected-yield result.",
    )
    out = OUT_DIR / "the69_raw_abcd_purity_companion_fine_et.png"
    fig.savefig(out, dpi=SLIDE_DPI, facecolor="white")
    plt.close(fig)
    return out


def make_ratio_slide(rows: list[dict]) -> Path:
    fig, axes = plt.subplots(1, 3, figsize=slide_figsize(), dpi=SLIDE_DPI)
    fig.subplots_adjust(left=0.070, right=0.965, top=0.755, bottom=0.190, wspace=0.235)
    add_title(
        fig,
        "Tight-ID yield-shape ratios: data compared with signal and inclusive MC",
        "Same fine ET bins as the yield slide. Ratios use area-normalized tight+isolated spectra, so they test shape agreement only, not absolute MC/data closure.",
    )

    ratio_styles = [
        ("Data / Signal MC", "Signal MC", "#C43C39", "o", -0.10),
        ("Data / Inclusive MC", "Inclusive MC", "#2F63C6", "s", 0.10),
    ]
    all_ratio_values: list[float] = []
    for ax, (cent, cent_label, _, _) in zip(axes, CENT):
        cent_rows = [r for r in rows if r["cent_key"] == cent]
        by_sample = {sample: [r for r in cent_rows if r["sample"] == sample] for sample, *_ in SAMPLES}
        totals = {
            sample: sum(float(r["yield_per_gev"]) * float(r["width"]) for r in pts if np.isfinite(float(r["yield_per_gev"])))
            for sample, pts in by_sample.items()
        }
        norm: dict[str, dict[tuple[float, float], tuple[float, float, float]]] = {}
        for sample, pts in by_sample.items():
            sample_points: dict[tuple[float, float], tuple[float, float, float]] = {}
            total = totals.get(sample, 0.0)
            for r in pts:
                key = (float(r["pt_lo"]), float(r["pt_hi"]))
                y = float(r["yield_per_gev"])
                e = float(r["yield_per_gev_err"])
                if total > 0:
                    sample_points[key] = (float(r["pt_mid"]), y / total, e / total)
            norm[sample] = sample_points

        for label, denom_sample, color, marker, offset in ratio_styles:
            xs, ratios, errs = [], [], []
            for lo, hi in PT_BINS:
                key = (float(lo), float(hi))
                if key not in norm.get("Data", {}) or key not in norm.get(denom_sample, {}):
                    continue
                x_mid, data_y, data_e = norm["Data"][key]
                _, den_y, den_e = norm[denom_sample][key]
                if not (np.isfinite(data_y) and np.isfinite(den_y)) or data_y <= 0 or den_y <= 0:
                    continue
                ratio = data_y / den_y
                err = ratio * math.sqrt((data_e / data_y) ** 2 + (den_e / den_y) ** 2)
                xs.append(x_mid + offset)
                ratios.append(ratio)
                errs.append(err)
                all_ratio_values.append(ratio)
            ax.errorbar(xs, ratios, yerr=errs, fmt=marker, color=color, markersize=6.0, capsize=2.2, linewidth=1.05, label=label)
        ax.axhline(1.0, color="#6B7280", lw=1.0, ls=(0, (4, 3)), alpha=0.75)
        ax.set_title(f"AuAu {cent_label}")
        ax.set_xlabel("cluster ET [GeV]")
        ax.set_xlim(13.2, 35.8)
        ax.grid(True, axis="y", alpha=0.26)
        if ax is axes[0]:
            ax.set_ylabel("shape ratio")
        else:
            ax.set_yticklabels([])

    ymax = max([2.0] + [v for v in all_ratio_values if np.isfinite(v)])
    for ax in axes:
        ax.set_ylim(0.0, min(4.0, max(1.8, ymax * 1.28)))

    handles = [
        plt.Line2D([0], [0], color=color, marker=marker, lw=0, markersize=7.0, label=label)
        for label, _, color, marker, _ in ratio_styles
    ]
    handles.append(plt.Line2D([0], [0], color="#6B7280", lw=1.0, ls=(0, (4, 3)), label="ratio = 1"))
    leg = fig.legend(
        handles=handles,
        loc="upper right",
        bbox_to_anchor=(0.935, 0.855),
        ncol=3,
        frameon=True,
        fancybox=False,
        edgecolor="#CBD5E1",
        facecolor="white",
        framealpha=0.96,
        borderpad=0.55,
        handletextpad=0.55,
        columnspacing=1.25,
        title="Shape-normalized ratios",
    )
    leg.get_title().set_fontweight("bold")
    add_takeaway(
        fig,
        "Ratios near unity mean the data ET shape follows that MC template after normalizing away total rate. Coherent shape ratios support the QA; absolute normalization still needs exposure, purity, and efficiency.",
    )
    out = OUT_DIR / "the69_tight_yield_shape_ratios_data_over_mc.png"
    fig.savefig(out, dpi=SLIDE_DPI, facecolor="white")
    plt.close(fig)
    return out


def normalized_points(rows: list[dict], cent: str, sample: str) -> dict[tuple[float, float], tuple[float, float, float]]:
    pts = [r for r in rows if r["cent_key"] == cent and r["sample"] == sample]
    total = sum(float(r["yield_per_gev"]) * float(r["width"]) for r in pts if np.isfinite(float(r["yield_per_gev"])))
    out: dict[tuple[float, float], tuple[float, float, float]] = {}
    if total <= 0:
        return out
    for r in pts:
        y = float(r["yield_per_gev"])
        e = float(r["yield_per_gev_err"])
        if not np.isfinite(y) or y <= 0:
            continue
        out[(float(r["pt_lo"]), float(r["pt_hi"]))] = (float(r["pt_mid"]), y / total, e / total)
    return out


def make_data_over_reco_signal_ratio_slide(rows: list[dict]) -> Path:
    fig, axes = plt.subplots(1, 3, figsize=slide_figsize(), dpi=SLIDE_DPI)
    fig.subplots_adjust(left=0.070, right=0.965, top=0.755, bottom=0.190, wspace=0.235)
    add_title(
        fig,
        "Tight+isolated yield-shape ratio: data / reco signal MC",
        "Fine cluster-ET bins. Data and reco signal MC spectra are each area-normalized within centrality, so this tests shape agreement, not absolute MC/data closure.",
    )

    all_ratios: list[float] = []
    for ax, (cent, cent_label, _, _) in zip(axes, CENT):
        data = normalized_points(rows, cent, "Data")
        sig = normalized_points(rows, cent, "Signal MC")
        xs, ratios, errs = [], [], []
        for lo, hi in PT_BINS:
            key = (float(lo), float(hi))
            if key not in data or key not in sig:
                continue
            x_mid, data_y, data_e = data[key]
            _, sig_y, sig_e = sig[key]
            if not (np.isfinite(data_y) and np.isfinite(sig_y)) or data_y <= 0 or sig_y <= 0:
                continue
            ratio = data_y / sig_y
            err = ratio * math.sqrt((data_e / data_y) ** 2 + (sig_e / sig_y) ** 2)
            xs.append(x_mid)
            ratios.append(ratio)
            errs.append(err)
            all_ratios.append(ratio)
        ax.errorbar(xs, ratios, yerr=errs, fmt="o", color="#111827", markersize=6.4, capsize=2.4, linewidth=1.05)
        ax.axhline(1.0, color="#6B7280", lw=1.0, ls=(0, (4, 3)), alpha=0.78)
        ax.set_title(f"AuAu {cent_label}")
        ax.set_xlabel("cluster ET [GeV]")
        ax.set_xlim(13.2, 35.8)
        ax.grid(True, axis="y", alpha=0.28)
        if ax is axes[0]:
            ax.set_ylabel("data / reco signal MC shape")
        else:
            ax.set_yticklabels([])

    finite = [v for v in all_ratios if np.isfinite(v)]
    ymin = min([0.72] + finite)
    ymax = max([1.28] + finite)
    pad = max(0.12, 0.18 * (ymax - ymin))
    for ax in axes:
        ax.set_ylim(max(0.0, ymin - pad), min(3.0, ymax + pad))

    handles = [
        plt.Line2D([0], [0], color="#111827", marker="o", lw=0, markersize=7.0, label="Data / reco signal MC"),
        plt.Line2D([0], [0], color="#6B7280", lw=1.0, ls=(0, (4, 3)), label="ratio = 1"),
    ]
    leg = fig.legend(
        handles=handles,
        loc="upper right",
        bbox_to_anchor=(0.935, 0.855),
        ncol=2,
        frameon=True,
        fancybox=False,
        edgecolor="#CBD5E1",
        facecolor="white",
        framealpha=0.96,
        borderpad=0.55,
        handletextpad=0.55,
        columnspacing=1.25,
        title="Shape-normalized comparison",
    )
    leg.get_title().set_fontweight("bold")
    add_takeaway(
        fig,
        "Use this as the first tight+isolated data/MC shape closure check. A final physics comparison still needs exposure, purity, and efficiency corrections.",
    )
    out = OUT_DIR / "the69_tight_iso_data_over_reco_signal_shape_ratio.png"
    fig.savefig(out, dpi=SLIDE_DPI, facecolor="white")
    plt.close(fig)
    return out


def make_data_vs_truth_signal_yield_slide(rows: list[dict]) -> Path:
    fig, axes = plt.subplots(1, 3, figsize=slide_figsize(), dpi=SLIDE_DPI)
    fig.subplots_adjust(left=0.070, right=0.965, top=0.755, bottom=0.190, wspace=0.235)
    add_title(
        fig,
        "Tight+isolated yield-shape QA: data vs truth-tagged signal MC",
        "Fine cluster-ET bins. Truth-tagged MC uses h_EisoReco_truthSigMatched_tight; both spectra are area-normalized within centrality.",
    )

    plot_specs = [
        ("Data", "#111827", "o", 0.0),
        ("Truth-tagged signal MC", "#059669", "D", 0.14),
    ]
    for ax, (cent, cent_label, _, _) in zip(axes, CENT):
        for sample, color, marker, offset in plot_specs:
            norm = normalized_points(rows, cent, sample)
            xs, ys, es = [], [], []
            for lo, hi in PT_BINS:
                key = (float(lo), float(hi))
                if key not in norm:
                    continue
                x_mid, y, e = norm[key]
                xs.append(x_mid + offset)
                ys.append(y)
                es.append(e)
            ax.errorbar(xs, ys, yerr=es, fmt=marker, color=color, markersize=6.0, capsize=2.2, linewidth=1.05, label=sample)
        ax.set_title(f"AuAu {cent_label}")
        ax.set_xlabel("cluster ET [GeV]")
        ax.set_yscale("log")
        ax.grid(True, axis="y", which="both", alpha=0.25)
        ax.set_xlim(13.2, 35.8)
        ax.set_ylim(1.0e-4, 3.0e-1)
        if ax is axes[0]:
            ax.set_ylabel("area-normalized tight yield / GeV")
        else:
            ax.set_yticklabels([])

    handles = [
        plt.Line2D([0], [0], color=color, marker=marker, lw=0, markersize=7.0, label=sample)
        for sample, color, marker, _ in plot_specs
    ]
    leg = fig.legend(
        handles=handles,
        loc="upper right",
        bbox_to_anchor=(0.935, 0.855),
        ncol=2,
        frameon=True,
        fancybox=False,
        edgecolor="#CBD5E1",
        facecolor="white",
        framealpha=0.96,
        borderpad=0.55,
        handletextpad=0.55,
        columnspacing=1.25,
        title="Shape-normalized yields",
    )
    leg.get_title().set_fontweight("bold")
    add_takeaway(
        fig,
        "This isolates the signal-like MC reference after truth matching. The plotted object is a shape QA view; absolute yield interpretation waits for the exposure and correction chain.",
    )
    out = OUT_DIR / "the69_tight_iso_data_vs_truth_tagged_signal_yield_shape.png"
    fig.savefig(out, dpi=SLIDE_DPI, facecolor="white")
    plt.close(fig)
    return out


def write_requested_speaker_scripts() -> tuple[Path, Path]:
    ratio_script = write_text_artifact(
        OUT_DIR / "the69_tight_iso_data_over_reco_signal_shape_ratio_speaker_script.md",
        """
        This slide is the first shape-closure check after the tight BDT and isolation selections.
        In each centrality bin I take the tight-plus-isolated data spectrum and the reco-level signal MC spectrum, normalize each spectrum to unit area within that centrality, and plot data divided by signal MC in fine cluster-ET bins.
        Because both spectra are area-normalized, this is not an absolute yield or luminosity comparison.
        What I want to see here is whether the data ET shape is broadly consistent with the reco signal MC template after the photon-ID selection.
        The low-ET bins sit close to unity, while the higher-ET tail tends to fall below the reco signal MC shape, with the peripheral bin becoming statistically sparse.
        The next step is not to over-interpret this as physics closure yet; it tells us where the exposure, purity, and efficiency corrections need to be checked carefully.
        """,
    )
    truth_script = write_text_artifact(
        OUT_DIR / "the69_tight_iso_data_vs_truth_tagged_signal_yield_shape_speaker_script.md",
        """
        This slide compares the same tight-plus-isolated data spectrum with the truth-tagged signal component in the embedded photon MC.
        The MC object is h_EisoReco_truthSigMatched_tight, so it is still the reconstructed tight and isolated candidate yield, but restricted to clusters matched to the signal photon truth.
        Each spectrum is area-normalized within the centrality panel and plotted on log-y to make the tail visible.
        The point of this slide is to isolate the signal-like MC reference from the inclusive MC mixture and check whether the selected data shape is sane against that reference.
        The data and truth-tagged signal MC agree best at low ET, while the data tail is lower than the truth-tagged signal shape in these current QA plots.
        That is a useful diagnostic, but the absolute interpretation waits for luminosity or sampled-event exposure, purity, and efficiency corrections.
        """,
    )
    return ratio_script, truth_script


def main() -> int:
    setup_style()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    yield_rows, abcd_rows = collect_points()
    write_csv(OUT_DIR / "the69_tight_yield_shape_points.csv", yield_rows)
    write_csv(OUT_DIR / "the69_raw_abcd_purity_points.csv", abcd_rows)
    yield_png = make_yield_slide(yield_rows)
    ratio_png = make_ratio_slide(yield_rows)
    data_reco_ratio_png = make_data_over_reco_signal_ratio_slide(yield_rows)
    truth_yield_png = make_data_vs_truth_signal_yield_slide(yield_rows)
    data_reco_ratio_script, truth_yield_script = write_requested_speaker_scripts()
    abcd_png = make_abcd_slide(abcd_rows)
    manifest = {
        "schema": "THE69_FINE_YIELD_ABCD_COMPANION_SLIDES_V1",
        "data_root": str(DATA_ROOT),
        "signal_root": str(SIGNAL_ROOT),
        "inclusive_root": str(INCLUSIVE_ROOT),
        "data_top": DATA_TOP,
        "sim_top": SIM_TOP,
        "pt_bins": PT_BINS,
        "yield_definition": "h_Eiso_tight_isoR40 fine-ET integrals, divided by bin width, then area-normalized per centrality/sample for the slide.",
        "truth_signal_definition": "h_EisoReco_truthSigMatched_tight_isoR40 fine-ET integrals from signal MC, divided by bin width, then area-normalized per centrality for the slide.",
        "abcd_definition": "h_Eiso_ABCD_{A,B,C,D}_isoR40_isSliding fine-ET integrals; raw purity=max(0,A-BC/D)/A.",
        "yield_png": str(yield_png),
        "ratio_definition": "Data/Signal MC and Data/Inclusive MC ratios after area-normalizing each tight+isolated ET spectrum within each centrality/sample.",
        "ratio_png": str(ratio_png),
        "data_over_reco_signal_ratio_definition": "Data/Signal MC ratio after area-normalizing each tight+isolated ET spectrum within each centrality/sample.",
        "data_over_reco_signal_ratio_png": str(data_reco_ratio_png),
        "data_over_reco_signal_ratio_speaker_script": str(data_reco_ratio_script),
        "data_vs_truth_signal_yield_png": str(truth_yield_png),
        "data_vs_truth_signal_yield_speaker_script": str(truth_yield_script),
        "abcd_png": str(abcd_png),
        "yield_csv": str(OUT_DIR / "the69_tight_yield_shape_points.csv"),
        "abcd_csv": str(OUT_DIR / "the69_raw_abcd_purity_points.csv"),
    }
    manifest_path = OUT_DIR / "the69_fine_yield_abcd_companion_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(yield_png)
    print(ratio_png)
    print(data_reco_ratio_png)
    print(truth_yield_png)
    print(abcd_png)
    print(manifest_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
