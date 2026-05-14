#!/usr/bin/env python3
"""Slide-candidate plots for the AuAu target80 BDT MC campaigns.

This script reads the final merged RecoilJets ROOT outputs from the completed
target80 campaign and produces compact, presentation-oriented comparisons:

  * truth-photon recovery vs E_T
  * tight-ID fraction among reconstructed photon candidates
  * inclusive-jet tight fraction
  * signal-recovery vs jet-leakage tradeoff

It intentionally keeps the plot text scientific-only. Slide titles and
takeaway text belong in Google Slides, not inside the PNG.
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import ROOT


ROOT.gROOT.SetBatch(True)


REPO = Path(__file__).resolve().parents[1]
TARGET_ROOT = REPO / "dataOutput/target80_all_available/bdt_target80_gated_20260512_001012"
FLAT050_ROOT = REPO / "dataOutput/combinedSimOnlyEMBEDDED"
OUT_DIR = TARGET_ROOT / "deep_dive_plots"

SIGNAL_REL = Path("photonJet12and20merged_SIM/RecoilJets_embeddedPhoton12plus20_MERGED.root")
BKG_REL = Path("embeddedJet12and20merged_SIM/RecoilJets_embeddedJet12plus20_MERGED.root")

CENTS = [
    ("0_20", "0-20% central"),
    ("20_50", "20-50% mid-central"),
    ("50_80", "50-80% peripheral"),
]

PT_BINS = [
    (5, 8),
    (8, 10),
    (10, 12),
    (12, 14),
    (14, 16),
    (16, 18),
    (18, 20),
    (20, 22),
    (22, 24),
    (24, 26),
    (26, 35),
]

COLORS = {
    "reference": "#111111",
    "blue": "#0072B2",
    "orange": "#E69F00",
    "green": "#009E73",
    "vermillion": "#D55E00",
    "purple": "#6F3FB5",
    "magenta": "#CC79A7",
    "sky": "#56B4E9",
    "gray": "#666666",
}


@dataclass(frozen=True)
class Variant:
    key: str
    label: str
    config_dir: str
    cfg: str
    pt_min: float
    pt_max: float
    color: str
    marker: str
    target80: bool = True
    cfg_signal: str | None = None
    cfg_bkg: str | None = None

    def signal_path(self) -> Path:
        cfg = self.cfg_signal or self.cfg
        if self.target80:
            return TARGET_ROOT / self.config_dir / "simembedded" / cfg / SIGNAL_REL
        return FLAT050_ROOT / cfg / SIGNAL_REL

    def bkg_path(self) -> Path:
        cfg = self.cfg_bkg or self.cfg
        if self.target80:
            return TARGET_ROOT / self.config_dir / "simembeddedinclusive" / cfg / BKG_REL
        return FLAT050_ROOT / cfg / BKG_REL


REFERENCE_CFG = "preselectionNewPPG12_tightReference_nonTightReference_baseVariant"
NEWPPG12_CFG = "preselectionNewPPG12_tightNewPPG12_nonTightReference_baseVariant"
CENTINPUT_CFG = "preselectionNewPPG12_tightAuAuCentInputBDT_nonTightAuAuBDTComplement_baseVariant"
CENTINPUT_3X3_CFG = "preselectionNewPPG12_tightAuAuCentInput3x3BDT_nonTightAuAuBDTComplement_baseVariant"
CENTINPUT_BASE3X3_CFG = "preselectionNewPPG12_tightAuAuCentInputBase3x3BDT_nonTightAuAuBDTComplement_baseVariant"
CENTINPUT_BASE3X3_LOWERCASE_SIGNAL = (
    "preselectionNewPPG12_tightAuauCentInputBase3x3BDT_nonTightAuAuBDTComplement_baseVariant"
)


REFERENCE = Variant(
    "reference_box",
    "Box-cuts",
    "analysis_config_expanded_5to40_target80",
    REFERENCE_CFG,
    5.0,
    35.0,
    COLORS["reference"],
    "o",
)

FLAT050_CENTINPUT = Variant(
    "flat050_centinput_5to40",
    "Cent. as input, score > 0.50",
    "",
    CENTINPUT_CFG,
    5.0,
    35.0,
    COLORS["gray"],
    "^",
    target80=False,
)


RECOVERY_CACHE: dict[tuple[str, str, float, float], list[dict]] = {}
TIGHT_CACHE: dict[tuple[str, str, float, float], list[dict]] = {}


def expanded_variants() -> list[Variant]:
    d = "analysis_config_expanded_5to40_target80"
    return [
        REFERENCE,
        Variant(
            "expanded_no_cent_5to40_t80",
            "No centrality input",
            d,
            "preselectionNewPPG12_tightAuAuNoCentBDT_nonTightAuAuBDTComplement_baseVariant",
            5.0,
            35.0,
            COLORS["gray"],
            "s",
        ),
        Variant(
            "expanded_centinput_5to40_t80",
            "Centrality as input",
            d,
            CENTINPUT_CFG,
            5.0,
            35.0,
            COLORS["blue"],
            "o",
        ),
        Variant(
            "expanded_centinput3x3_5to40_t80",
            "Centrality as input + local widths",
            d,
            CENTINPUT_3X3_CFG,
            5.0,
            35.0,
            COLORS["sky"],
            "D",
        ),
        Variant(
            "expanded_minopt_5to40_t80",
            "Minority-balanced",
            d,
            "preselectionNewPPG12_tightAuAuCentInputMinOptBDT_nonTightAuAuBDTComplement_baseVariant",
            5.0,
            35.0,
            COLORS["green"],
            "P",
        ),
        Variant(
            "expanded_cent3_5to40_t80",
            "3 centrality-bin BDTs",
            d,
            "preselectionNewPPG12_tightAuAuCent3BDT_nonTightAuAuBDTComplement_baseVariant",
            5.0,
            35.0,
            COLORS["orange"],
            "v",
        ),
        Variant(
            "expanded_cent7_5to40_t80",
            "7 centrality-bin BDTs",
            d,
            "preselectionNewPPG12_tightAuAuCent7BDT_nonTightAuAuBDTComplement_baseVariant",
            5.0,
            35.0,
            COLORS["purple"],
            "<",
        ),
        Variant(
            "expanded_ptbin_centinput_5to40_t80",
            r"$E_T$ bins + centrality input",
            d,
            "preselectionNewPPG12_tightAuAuPtBinCentInputBDT_nonTightAuAuBDTComplement_baseVariant",
            5.0,
            35.0,
            COLORS["magenta"],
            ">",
        ),
        Variant(
            "expanded_ptcent3_5to40_t80",
            r"$E_T$ x 3 centrality-bin BDTs",
            d,
            "preselectionNewPPG12_tightAuAuPtCent3BDT_nonTightAuAuBDTComplement_baseVariant",
            5.0,
            35.0,
            COLORS["vermillion"],
            "*",
        ),
        Variant(
            "expanded_ptcent7_5to40_t80",
            r"$E_T$ x 7 centrality-bin BDTs",
            d,
            "preselectionNewPPG12_tightAuAuPtCent7BDT_nonTightAuAuBDTComplement_baseVariant",
            5.0,
            35.0,
            COLORS["purple"],
            "X",
        ),
    ]


def width_variants(window: str, pt_min: float, pt_max: float) -> list[Variant]:
    d = f"analysis_config_widthstudy_{window}_target80"
    cfg_sig = None
    if window == "pt10to35":
        cfg_sig = CENTINPUT_BASE3X3_LOWERCASE_SIGNAL
    return [
        REFERENCE,
        Variant(
            f"width_{window}_base_t80",
            f"Base widths, {window_label(window)}",
            d,
            CENTINPUT_CFG,
            pt_min,
            pt_max,
            COLORS["blue"],
            "o",
        ),
        Variant(
            f"width_{window}_3x3_t80",
            f"Local 3x3 widths, {window_label(window)}",
            d,
            CENTINPUT_3X3_CFG,
            pt_min,
            pt_max,
            COLORS["orange"],
            "s",
        ),
        Variant(
            f"width_{window}_base3x3_t80",
            f"Base + local widths, {window_label(window)}",
            d,
            CENTINPUT_BASE3X3_CFG,
            pt_min,
            pt_max,
            COLORS["vermillion"],
            "D",
            cfg_signal=cfg_sig,
        ),
    ]


def etfine_variants() -> list[Variant]:
    d = "analysis_config_etfine_15to35_target80"
    return [
        REFERENCE,
        Variant(
            "etfine_global_base3x3_t80",
            "Global 15-35 BDT",
            d,
            CENTINPUT_BASE3X3_CFG,
            15.0,
            35.0,
            COLORS["blue"],
            "o",
        ),
        Variant(
            "etfine_centinput_t80",
            r"Fine $E_T$ bins + centrality input",
            d,
            "preselectionNewPPG12_tightAuAuEtFineCentInputBDT_nonTightAuAuBDTComplement_baseVariant",
            15.0,
            35.0,
            COLORS["green"],
            "s",
        ),
        Variant(
            "etfine_cent3_t80",
            r"Fine $E_T$ x 3 centrality-bin BDTs",
            d,
            "preselectionNewPPG12_tightAuAuEtFineCent3BDT_nonTightAuAuBDTComplement_baseVariant",
            15.0,
            35.0,
            COLORS["orange"],
            "D",
        ),
        Variant(
            "etfine_cent7_t80",
            r"Fine $E_T$ x 7 centrality-bin BDTs",
            d,
            "preselectionNewPPG12_tightAuAuEtFineCent7BDT_nonTightAuAuBDTComplement_baseVariant",
            15.0,
            35.0,
            COLORS["purple"],
            "^",
        ),
    ]


def window_label(window: str) -> str:
    return {
        "pt5to35": "5-35 GeV",
        "pt10to35": "10-35 GeV",
        "pt1530": "15-30 GeV",
        "pt15to35": "15-35 GeV",
    }.get(window, window)


def open_sim(path: Path) -> tuple[ROOT.TFile | None, ROOT.TDirectory | None]:
    if not path.exists():
        print(f"[WARN] missing file: {path}")
        return None, None
    f = ROOT.TFile.Open(str(path), "READ")
    if not f or f.IsZombie():
        print(f"[WARN] could not open: {path}")
        return None, None
    d = f.Get("SIM")
    if not d:
        print(f"[WARN] missing SIM directory: {path}")
        f.Close()
        return None, None
    return f, d


def get_hist(d: ROOT.TDirectory, name: str):
    h = d.Get(name)
    if not h:
        return None
    h.SetStats(False)
    return h


def sum_hist(h) -> tuple[float, float]:
    if not h:
        return 0.0, 0.0
    total = 0.0
    err2 = 0.0
    for ib in range(0, h.GetNbinsX() + 2):
        total += h.GetBinContent(ib)
        err2 += h.GetBinError(ib) ** 2
    return total, math.sqrt(err2)


def ratio_point(num: float, den_part: float, enum: float, eden_part: float) -> tuple[float, float, float]:
    den = num + den_part
    if den <= 0:
        return math.nan, math.nan, den
    value = num / den
    err = math.sqrt((den_part * enum) ** 2 + (num * eden_part) ** 2) / (den * den)
    if not math.isfinite(err):
        err = 0.0
    return value, err, den


def recovery_points(path: Path, cent: str, pt_min: float, pt_max: float) -> list[dict]:
    key = (str(path), cent, float(pt_min), float(pt_max))
    if key in RECOVERY_CACHE:
        return RECOVERY_CACHE[key]
    f, d = open_sim(path)
    if not d:
        RECOVERY_CACHE[key] = []
        return RECOVERY_CACHE[key]
    h_truth = get_hist(d, f"h_unfoldTruthPho_pTgamma_cent_{cent}")
    h_miss = get_hist(d, f"h_unfoldTruthPhoMisses_pTgamma_isoR30_isSliding_cent_{cent}")
    points: list[dict] = []
    if h_truth and h_miss:
        for ib in range(1, h_truth.GetNbinsX() + 1):
            lo = h_truth.GetXaxis().GetBinLowEdge(ib)
            hi = h_truth.GetXaxis().GetBinUpEdge(ib)
            mid = 0.5 * (lo + hi)
            if mid < pt_min or mid > pt_max:
                continue
            im = h_miss.GetXaxis().FindBin(mid)
            if im < 1 or im > h_miss.GetNbinsX():
                continue
            truth = h_truth.GetBinContent(ib)
            miss = h_miss.GetBinContent(im)
            if truth <= 0:
                continue
            value = 1.0 - miss / truth
            ey = 0.0
            if miss > 0:
                e_truth = h_truth.GetBinError(ib)
                e_miss = h_miss.GetBinError(im)
                miss_frac = miss / truth
                ey = miss_frac * math.hypot(e_miss / miss, e_truth / truth)
            points.append({"x": mid, "y": value, "ey": ey, "num": truth - miss, "den": truth})
    f.Close()
    RECOVERY_CACHE[key] = points
    return RECOVERY_CACHE[key]


def integrated_recovery(path: Path, cent: str, pt_min: float, pt_max: float) -> tuple[float, float, float]:
    points = recovery_points(path, cent, pt_min, pt_max)
    num = sum(p["num"] for p in points)
    den = sum(p["den"] for p in points)
    if den <= 0:
        return math.nan, math.nan, 0.0
    value = num / den
    err = math.sqrt(max(value * (1.0 - value), 0.0) / den) if den > 0 else math.nan
    return value, err, den


def tight_fraction_points(path: Path, cent: str, pt_min: float, pt_max: float) -> list[dict]:
    key = (str(path), cent, float(pt_min), float(pt_max))
    if key in TIGHT_CACHE:
        return TIGHT_CACHE[key]
    f, d = open_sim(path)
    if not d:
        TIGHT_CACHE[key] = []
        return TIGHT_CACHE[key]
    points: list[dict] = []
    for lo, hi in PT_BINS:
        mid = 0.5 * (lo + hi)
        if mid < pt_min or mid > pt_max:
            continue
        tag = f"isoR30_pT_{lo}_{hi}_cent_{cent}"
        h_tight = get_hist(d, f"h_Eiso_tight_{tag}")
        h_nontight = get_hist(d, f"h_Eiso_nonTight_{tag}")
        tight, e_tight = sum_hist(h_tight)
        nontight, e_nontight = sum_hist(h_nontight)
        value, err, den = ratio_point(tight, nontight, e_tight, e_nontight)
        if math.isfinite(value):
            points.append({"x": mid, "y": value, "ey": err, "num": tight, "den": den})
    f.Close()
    TIGHT_CACHE[key] = points
    return TIGHT_CACHE[key]


def integrated_tight_fraction(path: Path, cent: str, pt_min: float, pt_max: float) -> tuple[float, float, float]:
    points = tight_fraction_points(path, cent, pt_min, pt_max)
    num = sum(p["num"] for p in points)
    den = sum(p["den"] for p in points)
    if den <= 0:
        return math.nan, math.nan, 0.0
    value = num / den
    err = math.sqrt(max(value * (1.0 - value), 0.0) / den)
    return value, err, den


def mean_y(points: list[dict]) -> float:
    vals = [p["y"] for p in points if math.isfinite(p["y"])]
    return float(np.mean(vals)) if vals else math.nan


def smoothness(points: list[dict]) -> float:
    vals = [p["y"] for p in sorted(points, key=lambda p: p["x"]) if math.isfinite(p["y"])]
    if len(vals) < 3:
        return math.nan
    diffs = np.diff(vals)
    return float(np.sqrt(np.mean(diffs * diffs)))


def style_axes(ax, ylabel: str | None = None, xlabel: str | None = None):
    ax.tick_params(direction="in", which="both", top=True, right=True, length=6)
    ax.minorticks_on()
    ax.grid(True, color="#D0D0D0", alpha=0.45, linewidth=0.7)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=15)
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=15)
    ax.tick_params(labelsize=12)


def draw_sphenix_header(ax, y: float = 0.88, line2: str = "Pythia overlay, $\\sqrt{s_{NN}} = 200$ GeV"):
    ax.text(0.06, y, "sPHENIX", transform=ax.transAxes, ha="left", va="top",
            fontsize=13, fontweight="bold", fontstyle="italic")
    ax.text(0.265, y, "Internal", transform=ax.transAxes, ha="left", va="top", fontsize=13)
    ax.text(0.06, y - 0.10, line2, transform=ax.transAxes, ha="left", va="top", fontsize=11)
    ax.text(0.06, y - 0.17, "$R = 0.3$, sliding isolation", transform=ax.transAxes,
            ha="left", va="top", fontsize=11)


def errorbar_points(ax, pts: list[dict], variant: Variant, shift: float = 0.0, label: str | None = None):
    if not pts:
        return
    x = np.array([p["x"] + shift for p in pts], dtype=float)
    y = np.array([p["y"] for p in pts], dtype=float)
    ey = np.array([p["ey"] for p in pts], dtype=float)
    ax.errorbar(
        x,
        y,
        yerr=ey,
        fmt=variant.marker,
        linestyle="none",
        color=variant.color,
        markerfacecolor=variant.color,
        markeredgecolor=variant.color,
        markersize=6.4,
        elinewidth=1.4,
        capsize=2.4,
        label=label or variant.label,
        zorder=3,
    )


def plot_1x3(
    variants: list[Variant],
    out_name: str,
    ylabel: str,
    y_lim: tuple[float, float],
    point_reader: Callable[[Path, str, float, float], list[dict]],
    path_getter: Callable[[Variant], Path],
    pt_min: float,
    pt_max: float,
    shifts: Iterable[float] | None = None,
    header_y: float = 0.89,
):
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(17.0, 5.3), sharey=True)
    shifts_list = list(shifts if shifts is not None else np.linspace(-0.18, 0.18, len(variants)))
    for ic, (cent, title) in enumerate(CENTS):
        ax = axes[ic]
        for iv, variant in enumerate(variants):
            pts = point_reader(path_getter(variant), cent, max(pt_min, variant.pt_min), min(pt_max, variant.pt_max))
            errorbar_points(ax, pts, variant, shifts_list[iv] if len(variants) > 1 else 0.0)
        ax.set_title(title, fontsize=17, pad=12)
        ax.set_xlim(pt_min - 0.5, pt_max + 0.8)
        ax.set_ylim(*y_lim)
        x_title = r"Truth photon $E_T$ [GeV]" if "Truth-photon" in ylabel else r"Photon candidate $E_T$ [GeV]"
        style_axes(ax, ylabel if ic == 0 else None, x_title if ic == 1 else None)
        if ic == 0:
            draw_sphenix_header(ax, y=header_y)
        if ic == 2:
            leg = ax.legend(
                loc="upper right",
                bbox_to_anchor=(0.98, 0.97),
                frameon=True,
                framealpha=0.92,
                fontsize=10.5,
                ncol=1 if len(variants) <= 4 else 2,
                borderpad=0.55,
                handletextpad=0.5,
                columnspacing=0.9,
            )
            leg.get_frame().set_linewidth(0.0)
    fig.tight_layout(w_pad=1.5)
    out = OUT_DIR / out_name
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return out


def write_metric_summary(variants: list[Variant], out_csv: Path, pt_min: float = 15.0, pt_max: float = 35.0):
    seen: set[str] = set()
    rows = []
    for v in variants:
        if v.key in seen:
            continue
        seen.add(v.key)
        if v.key.startswith("reference"):
            continue
        sig_path = v.signal_path()
        bkg_path = v.bkg_path()
        if not sig_path.exists() or not bkg_path.exists():
            continue
        for cent, _ in CENTS:
            lo = max(pt_min, v.pt_min)
            hi = min(pt_max, v.pt_max)
            rec, rec_err, rec_den = integrated_recovery(sig_path, cent, lo, hi)
            sig_tight, sig_tight_err, sig_tight_den = integrated_tight_fraction(sig_path, cent, lo, hi)
            leak, leak_err, leak_den = integrated_tight_fraction(bkg_path, cent, lo, hi)
            rec_pts = recovery_points(sig_path, cent, lo, hi)
            rows.append(
                {
                    "model_key": v.key,
                    "model_label": v.label,
                    "centrality": cent,
                    "pt_min": lo,
                    "pt_max": hi,
                    "signal_recovery": rec,
                    "signal_recovery_err": rec_err,
                    "signal_recovery_den": rec_den,
                    "signal_candidate_tight_fraction": sig_tight,
                    "signal_candidate_tight_fraction_err": sig_tight_err,
                    "signal_candidate_tight_fraction_den": sig_tight_den,
                    "inclusive_jet_tight_fraction": leak,
                    "inclusive_jet_tight_fraction_err": leak_err,
                    "inclusive_jet_tight_fraction_den": leak_den,
                    "recovery_smoothness_rms_adjacent": smoothness(rec_pts),
                    "mean_recovery_points": mean_y(rec_pts),
                    "source_signal": str(sig_path.relative_to(REPO)),
                    "source_background": str(bkg_path.relative_to(REPO)),
                }
            )
    with out_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else [])
        if rows:
            writer.writeheader()
            writer.writerows(rows)
    return rows


def aggregate_rows(rows: list[dict]) -> list[dict]:
    by_key: dict[str, list[dict]] = {}
    for row in rows:
        by_key.setdefault(row["model_key"], []).append(row)
    out = []
    for key, group in by_key.items():
        label = group[0]["model_label"]
        # Central collisions are the priority. Penalize high inclusive-jet tight
        # fraction and visibly jagged recovery curves, but keep this as a
        # heuristic for plot selection rather than a physics selection rule.
        central = [r for r in group if r["centrality"] == "0_20"]
        rec_c = float(np.nanmean([r["signal_recovery"] for r in central]))
        leak_c = float(np.nanmean([r["inclusive_jet_tight_fraction"] for r in central]))
        smooth_c = float(np.nanmean([r["recovery_smoothness_rms_adjacent"] for r in central]))
        rec_all = float(np.nanmean([r["signal_recovery"] for r in group]))
        leak_all = float(np.nanmean([r["inclusive_jet_tight_fraction"] for r in group]))
        score = rec_c - 0.35 * leak_c - 0.40 * smooth_c
        out.append(
            {
                "model_key": key,
                "model_label": label,
                "central_signal_recovery": rec_c,
                "central_inclusive_jet_tight_fraction": leak_c,
                "central_recovery_smoothness": smooth_c,
                "mean_signal_recovery": rec_all,
                "mean_inclusive_jet_tight_fraction": leak_all,
                "plot_selection_score": score,
            }
        )
    out.sort(key=lambda r: r["plot_selection_score"], reverse=True)
    return out


def write_aggregate_csv(rows: list[dict], out_csv: Path):
    with out_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else [])
        if rows:
            writer.writeheader()
            writer.writerows(rows)


def plot_tradeoff(variants: list[Variant], out_name: str, pt_min: float = 15.0, pt_max: float = 35.0):
    fig, axes = plt.subplots(1, 3, figsize=(16.7, 5.2), sharey=True)
    for ic, (cent, title) in enumerate(CENTS):
        ax = axes[ic]
        for v in variants:
            if v.key.startswith("reference"):
                continue
            lo = max(pt_min, v.pt_min)
            hi = min(pt_max, v.pt_max)
            rec, rec_err, _ = integrated_recovery(v.signal_path(), cent, lo, hi)
            leak, leak_err, _ = integrated_tight_fraction(v.bkg_path(), cent, lo, hi)
            if not (math.isfinite(rec) and math.isfinite(leak)):
                continue
            ax.errorbar(
                leak,
                rec,
                xerr=leak_err,
                yerr=rec_err,
                fmt=v.marker,
                linestyle="none",
                color=v.color,
                markersize=7.0,
                capsize=2.2,
                elinewidth=1.2,
                label=v.label,
            )
        ax.set_title(title, fontsize=17, pad=12)
        ax.set_xlim(0.0, 0.62)
        ax.set_ylim(0.0, 0.62)
        style_axes(ax, "Truth-photon recovery" if ic == 0 else None, "Inclusive-jet tight fraction" if ic == 1 else None)
        if ic == 0:
            draw_sphenix_header(ax, y=0.92, line2="Pythia overlay, $15 < E_T^\\gamma < 35$ GeV")
        if ic == 2:
            leg = ax.legend(
                loc="lower right",
                frameon=True,
                framealpha=0.92,
                fontsize=10.0,
                ncol=1,
                borderpad=0.55,
            )
            leg.get_frame().set_linewidth(0.0)
    fig.tight_layout(w_pad=1.2)
    out = OUT_DIR / out_name
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return out


def plot_rank_bar(rows: list[dict], out_name: str):
    top = rows[:10]
    labels = [r["model_label"] for r in top][::-1]
    recovery = [r["central_signal_recovery"] for r in top][::-1]
    leakage = [r["central_inclusive_jet_tight_fraction"] for r in top][::-1]
    y = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(10.5, 6.2))
    ax.barh(y + 0.18, recovery, height=0.32, color=COLORS["blue"], label="Truth-photon recovery")
    ax.barh(y - 0.18, leakage, height=0.32, color=COLORS["vermillion"], label="Inclusive-jet tight fraction")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=10)
    ax.set_xlim(0.0, 0.62)
    style_axes(ax, None, "Integrated fraction, 0-20% central, 15-35 GeV")
    draw_sphenix_header(ax, y=0.98, line2="Pythia overlay, target-80 BDT working points")
    leg = ax.legend(loc="lower right", frameon=True, framealpha=0.92, fontsize=11)
    leg.get_frame().set_linewidth(0.0)
    fig.tight_layout()
    out = OUT_DIR / out_name
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return out


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    expanded = expanded_variants()
    width_1535 = width_variants("pt15to35", 15.0, 35.0)
    width_1530 = width_variants("pt1530", 15.0, 30.0)
    width_1035 = width_variants("pt10to35", 10.0, 35.0)
    width_0535 = width_variants("pt5to35", 5.0, 35.0)
    etfine = etfine_variants()

    all_for_ranking = []
    for group in [expanded, width_0535, width_1035, width_1530, width_1535, etfine]:
        for v in group:
            if v.key != REFERENCE.key:
                all_for_ranking.append(v)
    summary_rows = write_metric_summary(all_for_ranking, OUT_DIR / "target80_deep_dive_metric_summary.csv")
    aggregate = aggregate_rows(summary_rows)
    write_aggregate_csv(aggregate, OUT_DIR / "target80_deep_dive_model_ranking.csv")

    best_key = aggregate[0]["model_key"] if aggregate else "expanded_centinput3x3_5to40_t80"
    variant_by_key = {v.key: v for v in all_for_ranking}
    best = variant_by_key.get(best_key, expanded[3])

    # Simple story plots.
    plot_1x3(
        [REFERENCE, best],
        "best_target80_vs_boxcuts_truth_photon_recovery_1x3.png",
        "Truth-photon recovery",
        (0.0, 0.64),
        recovery_points,
        lambda v: v.signal_path(),
        5.0,
        35.0,
        shifts=[-0.08, 0.08],
        header_y=0.89,
    )
    plot_1x3(
        [REFERENCE, best],
        "best_target80_vs_boxcuts_signal_tight_fraction_1x3.png",
        "Tight-ID fraction among reco candidates",
        (0.0, 1.05),
        tight_fraction_points,
        lambda v: v.signal_path(),
        5.0,
        35.0,
        shifts=[-0.08, 0.08],
        header_y=0.89,
    )
    plot_1x3(
        [REFERENCE, best],
        "best_target80_vs_boxcuts_inclusive_jet_tight_fraction_1x3.png",
        "Inclusive-jet tight fraction",
        (0.0, 0.75),
        tight_fraction_points,
        lambda v: v.bkg_path(),
        5.0,
        35.0,
        shifts=[-0.08, 0.08],
        header_y=0.89,
    )

    # 5-40 expanded family: readable subset.
    expanded_core = [
        REFERENCE,
        expanded[2],
        expanded[3],
        expanded[5],
        expanded[6],
        expanded[8],
        expanded[9],
    ]
    plot_1x3(
        expanded_core,
        "expanded5to40_target80_core_family_truth_photon_recovery_1x3.png",
        "Truth-photon recovery",
        (0.0, 0.64),
        recovery_points,
        lambda v: v.signal_path(),
        5.0,
        35.0,
        header_y=0.89,
    )
    plot_1x3(
        expanded_core,
        "expanded5to40_target80_core_family_inclusive_jet_tight_fraction_1x3.png",
        "Inclusive-jet tight fraction",
        (0.0, 0.75),
        tight_fraction_points,
        lambda v: v.bkg_path(),
        5.0,
        35.0,
        header_y=0.89,
    )

    # Flat score >0.50 vs validation-derived target80 for the same 5-40 centrality-input model.
    plot_1x3(
        [REFERENCE, FLAT050_CENTINPUT, expanded[2], expanded[3]],
        "centinput_5to40_flat050_vs_target80_truth_photon_recovery_1x3.png",
        "Truth-photon recovery",
        (0.0, 0.64),
        recovery_points,
        lambda v: v.signal_path(),
        5.0,
        35.0,
        shifts=[-0.15, -0.05, 0.05, 0.15],
        header_y=0.89,
    )

    # Width study and Et-fine comparisons in the trusted high-pT window.
    plot_1x3(
        width_1535,
        "widthstudy_15to35_target80_truth_photon_recovery_1x3.png",
        "Truth-photon recovery",
        (0.0, 0.64),
        recovery_points,
        lambda v: v.signal_path(),
        15.0,
        35.0,
        shifts=[-0.15, -0.05, 0.05, 0.15],
        header_y=0.89,
    )
    plot_1x3(
        width_1535,
        "widthstudy_15to35_target80_signal_tight_fraction_1x3.png",
        "Tight-ID fraction among reco candidates",
        (0.0, 1.05),
        tight_fraction_points,
        lambda v: v.signal_path(),
        15.0,
        35.0,
        shifts=[-0.15, -0.05, 0.05, 0.15],
        header_y=0.89,
    )
    plot_1x3(
        [REFERENCE] + [v for v in etfine if not v.key.startswith("reference")],
        "etfine_15to35_target80_truth_photon_recovery_1x3.png",
        "Truth-photon recovery",
        (0.0, 0.64),
        recovery_points,
        lambda v: v.signal_path(),
        15.0,
        35.0,
        header_y=0.89,
    )

    plot_tradeoff(expanded_core, "expanded5to40_target80_signal_recovery_vs_jet_leakage_1x3.png")
    plot_rank_bar(aggregate, "target80_top_model_recovery_leakage_rank_0to20.png")

    print(f"[DONE] wrote deep-dive target80 plots to {OUT_DIR}")
    if aggregate:
        print("[TOP] plot-selection heuristic:")
        for row in aggregate[:6]:
            print(
                f"  {row['plot_selection_score']:.4f}  "
                f"recovery={row['central_signal_recovery']:.4f}  "
                f"jetTight={row['central_inclusive_jet_tight_fraction']:.4f}  "
                f"smooth={row['central_recovery_smoothness']:.4f}  "
                f"{row['model_label']} ({row['model_key']})"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
