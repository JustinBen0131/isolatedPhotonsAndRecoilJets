#!/usr/bin/env python3
"""Mine new target80 AuAu photon-ID money plots.

This script deliberately avoids regenerating the already-used first-look plots.
It uses the completed target80 RecoilJets MC outputs plus MLP validation CSVs to
make PPG12/ATLAS-inspired diagnostics:

  * all-variant signal-recovery vs background-control Pareto frontier
  * ABCD R-factor and truth-signal leakage diagnostics
  * truth-photon food-chain survival
  * shower-shape central/peripheral templates
  * MLP-vs-BDT validation challenger plots
"""

from __future__ import annotations

import csv
import math
import sys
import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import ROOT
from matplotlib.lines import Line2D


ROOT.gROOT.SetBatch(True)

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

import plot_target80_deep_dive as deep  # noqa: E402


TARGET_ROOT = deep.TARGET_ROOT
OUT = TARGET_ROOT / "money_plot_mining"
MLP_REPORT = (
    REPO
    / "dataOutput/auauTightMLPValidation/"
    / "mlp_model_validation_condor_deep_primary_ratios_nostat_fullval_20260512_145449"
)

CENTS = deep.CENTS
PT_BINS = [
    (15, 16),
    (16, 18),
    (18, 20),
    (20, 22),
    (22, 24),
    (24, 26),
    (26, 30),
    (30, 35),
]

REGIONS = [
    ("A", "isIsolated_isTight", 1),
    ("B", "notIsolated_isTight", 2),
    ("C", "isIsolated_notTight", 3),
    ("D", "notIsolated_notTight", 4),
]

FILE_CACHE: dict[str, tuple[ROOT.TFile, ROOT.TDirectory]] = {}

COLORS = {
    "box": "#111111",
    "bdt": "#0072B2",
    "bdt2": "#D55E00",
    "mlp": "#009E73",
    "gray": "#8A8A8A",
    "light": "#D6D6D6",
    "purple": "#6F3FB5",
    "orange": "#E69F00",
    "sky": "#56B4E9",
    "magenta": "#CC79A7",
}


@dataclass(frozen=True)
class NamedVariant:
    variant: deep.Variant
    short: str
    family: str


def all_target80_variants() -> list[NamedVariant]:
    items: list[NamedVariant] = []
    seen: set[str] = set()

    groups = [
        ("expanded 5-40", deep.expanded_variants()),
        ("width 5-35", deep.width_variants("pt5to35", 5.0, 35.0)),
        ("width 10-35", deep.width_variants("pt10to35", 10.0, 35.0)),
        ("width 15-30", deep.width_variants("pt1530", 15.0, 30.0)),
        ("width 15-35", deep.width_variants("pt15to35", 15.0, 35.0)),
        ("fine ET 15-35", deep.etfine_variants()),
    ]

    short_names = {
        "reference_box": "Box-cuts",
        "expanded_no_cent_5to40_t80": "No cent.",
        "expanded_centinput_5to40_t80": "Cent. input",
        "expanded_centinput3x3_5to40_t80": "Cent. input + 3x3",
        "expanded_minopt_5to40_t80": "Minority-balanced",
        "expanded_cent3_5to40_t80": "3 cent. BDTs",
        "expanded_cent7_5to40_t80": "7 cent. BDTs",
        "expanded_ptbin_centinput_5to40_t80": "ET bins + cent.",
        "expanded_ptcent3_5to40_t80": "ET x 3 cent.",
        "expanded_ptcent7_5to40_t80": "ET x 7 cent.",
        "width_pt5to35_base_t80": "5-35 base widths",
        "width_pt5to35_3x3_t80": "5-35 3x3 widths",
        "width_pt5to35_base3x3_t80": "5-35 base+3x3",
        "width_pt10to35_base_t80": "10-35 base widths",
        "width_pt10to35_3x3_t80": "10-35 3x3 widths",
        "width_pt10to35_base3x3_t80": "10-35 base+3x3",
        "width_pt1530_base_t80": "15-30 base widths",
        "width_pt1530_3x3_t80": "15-30 3x3 widths",
        "width_pt1530_base3x3_t80": "15-30 base+3x3",
        "width_pt15to35_base_t80": "15-35 base widths",
        "width_pt15to35_3x3_t80": "15-35 3x3 widths",
        "width_pt15to35_base3x3_t80": "15-35 base+3x3",
        "etfine_global_base3x3_t80": "15-35 global",
        "etfine_centinput_t80": "Fine ET + cent.",
        "etfine_cent3_t80": "Fine ET x 3 cent.",
        "etfine_cent7_t80": "Fine ET x 7 cent.",
    }

    for family, variants in groups:
        for v in variants:
            if v.key in seen:
                continue
            seen.add(v.key)
            items.append(NamedVariant(v, short_names.get(v.key, v.label), family))
    return items


def set_style(ax, xlabel: str | None = None, ylabel: str | None = None):
    ax.tick_params(direction="in", which="both", top=True, right=True, length=6)
    ax.minorticks_on()
    ax.grid(True, color="#D0D0D0", alpha=0.45, linewidth=0.8)
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=14)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=14)
    ax.tick_params(labelsize=11)


def header(ax, y: float = 0.93, line2: str = r"Pythia overlay, $\sqrt{s_{NN}} = 200$ GeV"):
    ax.text(
        0.055,
        y,
        "sPHENIX",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=13,
        fontweight="bold",
        fontstyle="italic",
    )
    ax.annotate(
        "Internal",
        xy=(0.055, y),
        xycoords=ax.transAxes,
        xytext=(73, 0),
        textcoords="offset points",
        ha="left",
        va="top",
        fontsize=13,
    )
    ax.text(0.055, y - 0.095, line2, transform=ax.transAxes, ha="left", va="top", fontsize=10.8)
    ax.text(0.055, y - 0.165, r"$R = 0.3$, sliding isolation", transform=ax.transAxes, ha="left", va="top", fontsize=10.8)


def save(fig, name: str) -> Path:
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / name
    fig.savefig(path, dpi=190)
    plt.close(fig)
    return path


def finite(x: float) -> bool:
    return math.isfinite(float(x))


def cached_dir(path: Path) -> ROOT.TDirectory | None:
    key = str(path)
    if key in FILE_CACHE:
        return FILE_CACHE[key][1]
    if not path.exists():
        print(f"[WARN] missing file: {path}")
        return None
    f = ROOT.TFile.Open(key, "READ")
    if not f or f.IsZombie():
        print(f"[WARN] could not open ROOT file: {path}")
        return None
    d = f.Get("SIM")
    if not d:
        print(f"[WARN] missing SIM directory: {path}")
        f.Close()
        return None
    FILE_CACHE[key] = (f, d)
    return d


def close_cached_files():
    for f, _ in FILE_CACHE.values():
        f.Close()
    FILE_CACHE.clear()


def hist_sum(h) -> tuple[float, float]:
    if not h:
        return 0.0, 0.0
    total = 0.0
    err2 = 0.0
    for ib in range(0, h.GetNbinsX() + 2):
        total += h.GetBinContent(ib)
        err2 += h.GetBinError(ib) ** 2
    return total, math.sqrt(err2)


def cached_integrated_recovery(path: Path, cent: str, pt_min: float, pt_max: float) -> tuple[float, float, float]:
    d = cached_dir(path)
    if not d:
        return math.nan, math.nan, 0.0
    h_truth = deep.get_hist(d, f"h_unfoldTruthPho_pTgamma_cent_{cent}")
    h_miss = deep.get_hist(d, f"h_unfoldTruthPhoMisses_pTgamma_isoR30_isSliding_cent_{cent}")
    if not h_truth or not h_miss:
        return math.nan, math.nan, 0.0
    truth_sum = 0.0
    miss_sum = 0.0
    truth_err2 = 0.0
    miss_err2 = 0.0
    for ib in range(1, h_truth.GetNbinsX() + 1):
        mid = h_truth.GetXaxis().GetBinCenter(ib)
        if not (pt_min <= mid <= pt_max):
            continue
        im = h_miss.GetXaxis().FindBin(mid)
        if im < 1 or im > h_miss.GetNbinsX():
            continue
        truth_sum += h_truth.GetBinContent(ib)
        truth_err2 += h_truth.GetBinError(ib) ** 2
        miss_sum += h_miss.GetBinContent(im)
        miss_err2 += h_miss.GetBinError(im) ** 2
    if truth_sum <= 0:
        return math.nan, math.nan, 0.0
    value = 1.0 - miss_sum / truth_sum
    if miss_sum > 0:
        err = (miss_sum / truth_sum) * math.hypot(math.sqrt(miss_err2) / miss_sum, math.sqrt(truth_err2) / truth_sum)
    else:
        err = 0.0
    return value, err, truth_sum


def cached_integrated_tight_fraction(path: Path, cent: str, pt_min: float, pt_max: float) -> tuple[float, float, float]:
    d = cached_dir(path)
    if not d:
        return math.nan, math.nan, 0.0
    tight_sum = 0.0
    nontight_sum = 0.0
    tight_err2 = 0.0
    nontight_err2 = 0.0
    for lo, hi in deep.PT_BINS:
        mid = 0.5 * (lo + hi)
        if not (pt_min <= mid <= pt_max):
            continue
        tag = f"isoR30_pT_{lo}_{hi}_cent_{cent}"
        h_tight = deep.get_hist(d, f"h_Eiso_tight_{tag}")
        h_nontight = deep.get_hist(d, f"h_Eiso_nonTight_{tag}")
        t, et = hist_sum(h_tight)
        nt, ent = hist_sum(h_nontight)
        tight_sum += t
        nontight_sum += nt
        tight_err2 += et * et
        nontight_err2 += ent * ent
    den = tight_sum + nontight_sum
    if den <= 0:
        return math.nan, math.nan, 0.0
    value = tight_sum / den
    err = math.sqrt((nontight_sum * math.sqrt(tight_err2)) ** 2 + (tight_sum * math.sqrt(nontight_err2)) ** 2) / (den * den)
    return value, err, den


def pareto_indices(points: list[dict]) -> list[int]:
    winners = []
    for i, p in enumerate(points):
        dominated = False
        for j, q in enumerate(points):
            if i == j:
                continue
            better_or_equal = q["bkg"] <= p["bkg"] and q["rec"] >= p["rec"]
            strictly_better = q["bkg"] < p["bkg"] or q["rec"] > p["rec"]
            if better_or_equal and strictly_better:
                dominated = True
                break
        if not dominated:
            winners.append(i)
    return winners


def build_frontier_rows() -> list[dict]:
    rows: list[dict] = []
    items = all_target80_variants()
    for idx, item in enumerate(items, start=1):
        v = item.variant
        print(f"[frontier] {idx}/{len(items)} {item.short} ({item.family})", flush=True)
        if not v.signal_path().exists() or not v.bkg_path().exists():
            continue
        for cent, cent_label in CENTS:
            lo = max(15.0, v.pt_min)
            hi = min(35.0, v.pt_max)
            if hi <= lo:
                continue
            rec, rec_err, rec_den = cached_integrated_recovery(v.signal_path(), cent, lo, hi)
            bkg, bkg_err, bkg_den = cached_integrated_tight_fraction(v.bkg_path(), cent, lo, hi)
            sig_tight, sig_tight_err, sig_den = cached_integrated_tight_fraction(v.signal_path(), cent, lo, hi)
            if finite(rec) and finite(bkg):
                rows.append(
                    {
                        "key": v.key,
                        "label": item.short,
                        "family": item.family,
                        "centrality": cent,
                        "centrality_label": cent_label,
                        "pt_min": lo,
                        "pt_max": hi,
                        "rec": rec,
                        "rec_err": rec_err,
                        "rec_den": rec_den,
                        "bkg": bkg,
                        "bkg_err": bkg_err,
                        "bkg_den": bkg_den,
                        "sig_tight": sig_tight,
                        "sig_tight_err": sig_tight_err,
                        "sig_tight_den": sig_den,
                    }
                )
    return rows


def build_frontier_rows_from_existing_expanded_csv() -> list[dict]:
    csv_path = (
        TARGET_ROOT
        / "analysis_config_expanded_5to40_target80"
        / "first_look_plots"
        / "target80_integrated_metric_summary_15to35.csv"
    )
    if not csv_path.exists():
        return []
    rows = []
    label_map = {
        "box": "Box-cuts",
        "noCent": "No cent.",
        "centInput": "Cent. input",
        "minOpt": "Minority-balanced",
        "cent3": "3 cent. BDTs",
        "cent7": "7 cent. BDTs",
        "ptCent3": "ET x 3 cent.",
    }
    key_map = {
        "box": "reference_box",
        "noCent": "expanded_no_cent_5to40_t80",
        "centInput": "expanded_centinput_5to40_t80",
        "minOpt": "expanded_minopt_5to40_t80",
        "cent3": "expanded_cent3_5to40_t80",
        "cent7": "expanded_cent7_5to40_t80",
        "ptCent3": "expanded_ptcent3_5to40_t80",
    }
    for r in csv.DictReader(csv_path.open()):
        key = r["model_key"]
        if key not in key_map:
            continue
        rows.append(
            {
                "key": key_map[key],
                "label": label_map.get(key, r["model_label"]),
                "family": "expanded 5-40 target80",
                "centrality": r["centrality"],
                "centrality_label": dict(CENTS).get(r["centrality"], r["centrality"]),
                "pt_min": 15.0,
                "pt_max": 35.0,
                "rec": float(r["signal_recovery"]),
                "rec_err": float(r["signal_recovery_err"]),
                "rec_den": math.nan,
                "bkg": float(r["background_tight_fraction"]),
                "bkg_err": float(r["background_tight_fraction_err"]),
                "bkg_den": math.nan,
                "sig_tight": float(r["signal_candidate_tight_fraction"]),
                "sig_tight_err": float(r["signal_candidate_tight_fraction_err"]),
                "sig_tight_den": math.nan,
            }
        )
    return rows


def write_csv(rows: list[dict], path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        return
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def plot_pareto_frontier(rows: list[dict]) -> Path:
    fig, axes = plt.subplots(1, 3, figsize=(17.5, 5.5), sharey=True)
    for ic, (cent, title) in enumerate(CENTS):
        ax = axes[ic]
        pts = [r for r in rows if r["centrality"] == cent]
        winners = pareto_indices(pts)
        winner_set = {pts[i]["key"] for i in winners}

        ax.errorbar(
            [p["bkg"] for p in pts],
            [p["rec"] for p in pts],
            xerr=[p["bkg_err"] for p in pts],
            yerr=[p["rec_err"] for p in pts],
            fmt="o",
            linestyle="none",
            markersize=4.2,
            color=COLORS["light"],
            ecolor="#C0C0C0",
            elinewidth=0.8,
            capsize=1.5,
            zorder=1,
        )

        # Draw the Pareto front in increasing background-fraction order.
        front = sorted([p for p in pts if p["key"] in winner_set], key=lambda r: r["bkg"])
        if front:
            ax.plot([p["bkg"] for p in front], [p["rec"] for p in front], color="#444444", lw=1.4, zorder=2)
        best_bdt = None
        for p in front:
            color = COLORS["box"] if p["key"] == "reference_box" else COLORS["bdt"]
            marker = "o" if p["key"] == "reference_box" else "D"
            ax.errorbar(
                p["bkg"],
                p["rec"],
                xerr=p["bkg_err"],
                yerr=p["rec_err"],
                fmt=marker,
                color=color,
                markersize=7.5,
                capsize=2.2,
                elinewidth=1.2,
                zorder=4,
            )
            if p["key"] != "reference_box" and (best_bdt is None or p["rec"] > best_bdt["rec"]):
                best_bdt = p

        box = next((p for p in pts if p["key"] == "reference_box"), None)
        if box:
            ax.errorbar(
                box["bkg"],
                box["rec"],
                xerr=box["bkg_err"],
                yerr=box["rec_err"],
                fmt="o",
                color=COLORS["box"],
                markersize=7.5,
                capsize=2.2,
                elinewidth=1.2,
                zorder=5,
            )
            ax.annotate(
                "Box-cuts",
                (box["bkg"], box["rec"]),
                xytext=(7, -12),
                textcoords="offset points",
                fontsize=9.4,
                ha="left",
                va="top",
            )
        if best_bdt:
            ax.annotate(
                "BDT frontier",
                (best_bdt["bkg"], best_bdt["rec"]),
                xytext=(-10, 9),
                textcoords="offset points",
                fontsize=9.4,
                ha="right",
                va="bottom",
                color=COLORS["bdt"],
            )

        ax.set_title(title, fontsize=17, pad=12)
        ax.set_xlim(0.05, 0.56)
        ax.set_ylim(0.10, 0.52)
        set_style(
            ax,
            "Inclusive-jet tight fraction" if ic == 1 else None,
            "Truth-photon recovery" if ic == 0 else None,
        )
        if ic == 0:
            header(ax, y=0.94, line2=r"Pythia overlay, $15 < E_T^\gamma < 35$ GeV")
        if ic == 2:
            handles = [
                Line2D([0], [0], marker="o", color="none", markerfacecolor=COLORS["light"], markeredgecolor=COLORS["light"], markersize=6, label="tested target-80 variants"),
                Line2D([0], [0], marker="D", color=COLORS["bdt"], markerfacecolor=COLORS["bdt"], markersize=7, label="Pareto-winning BDT"),
                Line2D([0], [0], marker="o", color=COLORS["box"], markerfacecolor=COLORS["box"], markersize=7, label="box-cuts"),
                Line2D([0], [0], color="#444444", lw=1.4, label="Pareto frontier"),
            ]
            leg = ax.legend(handles=handles, loc="lower right", frameon=True, framealpha=0.92, fontsize=10.0)
            leg.get_frame().set_linewidth(0)

    fig.tight_layout(w_pad=1.15)
    return save(fig, "money_pareto_signal_recovery_vs_inclusive_jet_tight_fraction.png")


def open_sim(path: Path):
    f, d = deep.open_sim(path)
    return f, d


def hsum(h) -> tuple[float, float]:
    if not h:
        return 0.0, 0.0
    total = 0.0
    err2 = 0.0
    for ib in range(0, h.GetNbinsX() + 2):
        total += h.GetBinContent(ib)
        err2 += h.GetBinError(ib) ** 2
    return total, math.sqrt(err2)


def count_region(d, base: str, region_token: str, pt: tuple[int, int], cent: str) -> tuple[float, float]:
    lo, hi = pt
    name = f"{base}_{region_token}_isoR30_isSliding_pT_{lo}_{hi}_cent_{cent}"
    return hsum(deep.get_hist(d, name))


def signal_abcd_counts(d, pt: tuple[int, int], cent: str, base: str = "h_sigABCD_MC") -> dict[str, tuple[float, float]]:
    lo, hi = pt
    h = deep.get_hist(d, f"{base}_isoR30_isSliding_pT_{lo}_{hi}_cent_{cent}")
    out = {}
    for reg, _, ibin in REGIONS:
        if h:
            out[reg] = (h.GetBinContent(ibin), h.GetBinError(ibin))
        else:
            out[reg] = (0.0, 0.0)
    return out


def subtract_counts(total: tuple[float, float], sig: tuple[float, float]) -> tuple[float, float]:
    val = max(total[0] - sig[0], 0.0)
    err = math.hypot(total[1], sig[1])
    return val, err


def ratio_with_error(num: tuple[float, float], den: tuple[float, float]) -> tuple[float, float]:
    if den[0] <= 0 or num[0] < 0:
        return math.nan, math.nan
    value = num[0] / den[0]
    err = value * math.hypot(num[1] / num[0] if num[0] > 0 else 0.0, den[1] / den[0])
    return value, err


def rbkg_with_error(a, b, c, d) -> tuple[float, float]:
    if min(a[0], b[0], c[0], d[0]) <= 0:
        return math.nan, math.nan
    value = (a[0] * d[0]) / (b[0] * c[0])
    rel2 = 0.0
    for val, err in (a, b, c, d):
        rel2 += (err / val) ** 2 if val > 0 else 0.0
    return value, abs(value) * math.sqrt(rel2)


def build_abcd_rows(variant: deep.Variant, label: str) -> list[dict]:
    rows: list[dict] = []
    f_bkg, d_bkg = open_sim(variant.bkg_path())
    f_sig, d_sig = open_sim(variant.signal_path())
    if not d_bkg or not d_sig:
        if f_bkg:
            f_bkg.Close()
        if f_sig:
            f_sig.Close()
        return rows

    for cent, _ in CENTS:
        for pt in PT_BINS:
            mid = 0.5 * (pt[0] + pt[1])
            totals = {
                reg: count_region(d_bkg, "h_xJpurityLead", token, pt, cent)
                for reg, token, _ in REGIONS
            }
            sig_lead = signal_abcd_counts(d_bkg, pt, cent, "h_xJpurityLead_sigABCD_MC")
            bkg_only = {reg: subtract_counts(totals[reg], sig_lead[reg]) for reg, _, _ in REGIONS}
            r_raw, r_raw_err = rbkg_with_error(totals["A"], totals["B"], totals["C"], totals["D"])
            r_bkg, r_bkg_err = rbkg_with_error(bkg_only["A"], bkg_only["B"], bkg_only["C"], bkg_only["D"])

            sig_abcd = signal_abcd_counts(d_sig, pt, cent, "h_sigABCD_MC")
            a_sig = sig_abcd["A"]
            for reg in ("B", "C", "D"):
                frac, frac_err = ratio_with_error(sig_abcd[reg], a_sig)
                rows.append(
                    {
                        "model": label,
                        "centrality": cent,
                        "pt_low": pt[0],
                        "pt_high": pt[1],
                        "pt_mid": mid,
                        "quantity": f"c{reg}",
                        "value": frac,
                        "error": frac_err,
                    }
                )
            rows.append(
                {
                    "model": label,
                    "centrality": cent,
                    "pt_low": pt[0],
                    "pt_high": pt[1],
                    "pt_mid": mid,
                    "quantity": "R_raw_inclusive_jet",
                    "value": r_raw,
                    "error": r_raw_err,
                }
            )
            rows.append(
                {
                    "model": label,
                    "centrality": cent,
                    "pt_low": pt[0],
                    "pt_high": pt[1],
                    "pt_mid": mid,
                    "quantity": "R_bkg_truth_subtracted",
                    "value": r_bkg,
                    "error": r_bkg_err,
                }
            )
    f_bkg.Close()
    f_sig.Close()
    return rows


def plot_abcd_stack(rows: list[dict]) -> list[Path]:
    paths: list[Path] = []
    selected_models = ["Box-cuts", "15-35 base+3x3 BDT"]
    fig, axes = plt.subplots(1, 3, figsize=(17.4, 5.2), sharey=True)
    for ic, (cent, title) in enumerate(CENTS):
        ax = axes[ic]
        for model, color, marker in [
            ("Box-cuts", COLORS["box"], "o"),
            ("15-35 base+3x3 BDT", COLORS["bdt"], "s"),
        ]:
            pts = [
                r
                for r in rows
                if r["model"] == model
                and r["centrality"] == cent
                and r["quantity"] == "R_bkg_truth_subtracted"
                and finite(r["value"])
            ]
            pts = sorted(pts, key=lambda r: r["pt_mid"])
            ax.errorbar(
                [r["pt_mid"] for r in pts],
                [r["value"] for r in pts],
                yerr=[r["error"] for r in pts],
                fmt=marker,
                linestyle="none",
                color=color,
                markersize=6.5,
                elinewidth=1.3,
                capsize=2.2,
                label=model,
            )
        ax.axhline(1.0, color="#555555", lw=1.2, ls="--")
        ax.set_title(title, fontsize=17, pad=12)
        ax.set_xlim(14.3, 35.8)
        ax.set_ylim(0.0, 3.4)
        set_style(ax, r"Photon candidate $E_T$ [GeV]" if ic == 1 else None, r"$R_{\mathrm{bkg}} = A_{\mathrm{bkg}}D_{\mathrm{bkg}}/(B_{\mathrm{bkg}}C_{\mathrm{bkg}})$" if ic == 0 else None)
        if ic == 0:
            header(ax, y=0.94, line2=r"Pythia overlay, truth-subtracted inclusive-jet background")
        if ic == 2:
            leg = ax.legend(loc="upper right", frameon=True, framealpha=0.93, fontsize=10.5)
            leg.get_frame().set_linewidth(0)
    fig.tight_layout(w_pad=1.2)
    paths.append(save(fig, "money_abcd_truth_subtracted_rfactor_box_vs_bdt.png"))

    fig, axes = plt.subplots(1, 3, figsize=(17.4, 5.2), sharey=True)
    qstyles = {
        "cB": (COLORS["orange"], "o", r"$c_B$: tight, non-isolated"),
        "cC": (COLORS["bdt"], "s", r"$c_C$: non-tight, isolated"),
        "cD": (COLORS["purple"], "D", r"$c_D$: non-tight, non-isolated"),
    }
    for ic, (cent, title) in enumerate(CENTS):
        ax = axes[ic]
        for q, (color, marker, lab) in qstyles.items():
            pts = [
                r
                for r in rows
                if r["model"] == "15-35 base+3x3 BDT"
                and r["centrality"] == cent
                and r["quantity"] == q
                and finite(r["value"])
            ]
            pts = sorted(pts, key=lambda r: r["pt_mid"])
            ax.errorbar(
                [r["pt_mid"] for r in pts],
                [r["value"] for r in pts],
                yerr=[r["error"] for r in pts],
                fmt=marker,
                linestyle="none",
                color=color,
                markersize=6.5,
                elinewidth=1.3,
                capsize=2.2,
                label=lab,
            )
        ax.set_title(title, fontsize=17, pad=12)
        ax.set_xlim(14.3, 35.8)
        ax.set_ylim(0.0, 1.10)
        set_style(ax, r"Photon candidate $E_T$ [GeV]" if ic == 1 else None, "Truth-signal leakage / Region A" if ic == 0 else None)
        if ic == 0:
            header(ax, y=0.94, line2=r"Pythia overlay, truth-matched signal photons")
        if ic == 2:
            leg = ax.legend(loc="upper right", frameon=True, framealpha=0.93, fontsize=10.5)
            leg.get_frame().set_linewidth(0)
    fig.tight_layout(w_pad=1.2)
    paths.append(save(fig, "money_abcd_signal_leakage_factors_bdt.png"))
    return paths


def get_1d_sum(d, names: Iterable[str]) -> tuple[float, float]:
    val = 0.0
    err2 = 0.0
    for name in names:
        h = deep.get_hist(d, name)
        if not h:
            continue
        s, e = hsum(h)
        val += s
        err2 += e * e
    return val, math.sqrt(err2)


def food_chain_counts(variant: deep.Variant, cent: str, pt_min: float = 15.0, pt_max: float = 35.0) -> dict[str, tuple[float, float]]:
    f, d = open_sim(variant.signal_path())
    if not d:
        if f:
            f.Close()
        return {}

    # Denominator: truth-isolated prompt photons from h_unfoldTruthPho.
    truth = 0.0
    truth_err2 = 0.0
    h_truth = deep.get_hist(d, f"h_unfoldTruthPho_pTgamma_cent_{cent}")
    if h_truth:
        for ib in range(1, h_truth.GetNbinsX() + 1):
            mid = h_truth.GetXaxis().GetBinCenter(ib)
            if pt_min <= mid <= pt_max:
                truth += h_truth.GetBinContent(ib)
                truth_err2 += h_truth.GetBinError(ib) ** 2

    sig_abcd_total = {"A": (0.0, 0.0), "B": (0.0, 0.0), "C": (0.0, 0.0), "D": (0.0, 0.0)}
    lead_abcd_total = {"A": (0.0, 0.0), "B": (0.0, 0.0), "C": (0.0, 0.0), "D": (0.0, 0.0)}
    for pt in PT_BINS:
        mid = 0.5 * (pt[0] + pt[1])
        if not (pt_min <= mid <= pt_max):
            continue
        sig = signal_abcd_counts(d, pt, cent, "h_sigABCD_MC")
        lead = signal_abcd_counts(d, pt, cent, "h_xJpurityLead_sigABCD_MC")
        for reg in sig_abcd_total:
            sig_abcd_total[reg] = (
                sig_abcd_total[reg][0] + sig[reg][0],
                math.hypot(sig_abcd_total[reg][1], sig[reg][1]),
            )
            lead_abcd_total[reg] = (
                lead_abcd_total[reg][0] + lead[reg][0],
                math.hypot(lead_abcd_total[reg][1], lead[reg][1]),
            )

    a = sig_abcd_total["A"]
    b = sig_abcd_total["B"]
    c = sig_abcd_total["C"]
    dreg = sig_abcd_total["D"]
    matched_pre = (a[0] + b[0] + c[0] + dreg[0], math.sqrt(a[1] ** 2 + b[1] ** 2 + c[1] ** 2 + dreg[1] ** 2))
    tight = (a[0] + b[0], math.hypot(a[1], b[1]))
    isolated = (a[0] + c[0], math.hypot(a[1], c[1]))
    region_a = a
    xj_lead = lead_abcd_total["A"]
    f.Close()
    return {
        "truth signal": (truth, math.sqrt(truth_err2)),
        "reco + preselection": matched_pre,
        "tight ID": tight,
        "isolated": isolated,
        "Region A": region_a,
        "event-leading Region A": xj_lead,
    }


def plot_food_chain(box: deep.Variant, bdt: deep.Variant) -> Path:
    stages = ["truth signal", "reco + preselection", "tight ID", "isolated", "Region A", "event-leading Region A"]
    fig, axes = plt.subplots(1, 3, figsize=(17.6, 5.45), sharey=True)
    for ic, (cent, title) in enumerate(CENTS):
        ax = axes[ic]
        for variant, color, marker, label, offset in [
            (box, COLORS["box"], "o", "Box-cuts", -0.055),
            (bdt, COLORS["bdt"], "s", "15-35 base+3x3 BDT", 0.055),
        ]:
            counts = food_chain_counts(variant, cent)
            den = counts.get("truth signal", (0.0, 0.0))[0]
            xs = np.arange(len(stages), dtype=float) + offset
            ys = []
            yerr = []
            for st in stages:
                val, err = counts.get(st, (0.0, 0.0))
                frac = val / den if den > 0 else math.nan
                frac_err = frac * math.hypot(err / val if val > 0 else 0.0, counts["truth signal"][1] / den if den > 0 else 0.0) if den > 0 else math.nan
                ys.append(frac)
                yerr.append(frac_err)
            ax.errorbar(xs, ys, yerr=yerr, fmt=marker, linestyle="-", color=color, markersize=6.2, lw=1.5, capsize=2.3, label=label)
        ax.set_title(title, fontsize=17, pad=12)
        ax.set_xticks(np.arange(len(stages)))
        ax.set_xticklabels(["truth", "reco+\npre", "tight", "iso", "A", "lead\nA"], fontsize=10.5)
        ax.set_ylim(0.0, 1.08)
        set_style(ax, None, "Fraction of truth-isolated photons" if ic == 0 else None)
        if ic == 0:
            header(ax, y=0.94, line2=r"Pythia overlay, $15 < E_T^\gamma < 35$ GeV")
        if ic == 2:
            leg = ax.legend(loc="upper right", frameon=True, framealpha=0.93, fontsize=10.5)
            leg.get_frame().set_linewidth(0)
    fig.tight_layout(w_pad=1.2)
    return save(fig, "money_truth_photon_food_chain_box_vs_bdt.png")


def clone_norm(h):
    if not h:
        return None
    out = h.Clone(f"{h.GetName()}_clone")
    out.SetDirectory(0)
    s = out.Integral(0, out.GetNbinsX() + 1)
    if s > 0:
        out.Scale(1.0 / s)
    return out


def hist_to_step(h):
    xs = []
    ys = []
    if not h:
        return np.array([]), np.array([])
    for ib in range(1, h.GetNbinsX() + 1):
        lo = h.GetXaxis().GetBinLowEdge(ib)
        hi = h.GetXaxis().GetBinUpEdge(ib)
        y = h.GetBinContent(ib)
        xs.extend([lo, hi])
        ys.extend([y, y])
    return np.array(xs), np.array(ys)


def add_hists(hists):
    acc = None
    for h in hists:
        if not h:
            continue
        if acc is None:
            acc = h.Clone(f"{h.GetName()}_sum")
            acc.SetDirectory(0)
        else:
            acc.Add(h)
    return acc


def fetch_ss_sum(d, var: str, tag: str, cent: str, pt_bins: list[tuple[int, int]]):
    hists = []
    for lo, hi in pt_bins:
        for name in [
            f"h_ss_{var}_{tag}_pT_{lo}_{hi}_cent_{cent}",
            f"h_ss_{var}_{tag}_pT_{lo}_{hi}",
        ]:
            h = deep.get_hist(d, name)
            if h:
                hists.append(h)
                break
    return clone_norm(add_hists(hists))


def plot_shower_shape_templates(bdt: deep.Variant) -> Path:
    f_sig, d_sig = open_sim(bdt.signal_path())
    f_bkg, d_bkg = open_sim(bdt.bkg_path())
    if not d_sig or not d_bkg:
        if f_sig:
            f_sig.Close()
        if f_bkg:
            f_bkg.Close()
        raise RuntimeError("Missing signal/background ROOT input for shower-shape plot")

    variables = [
        ("weta", r"base $\eta$ width"),
        ("weta33", r"local 3x3 $\eta$ width"),
        ("wphi", r"base $\phi$ width"),
        ("wphi33", r"local 3x3 $\phi$ width"),
    ]
    cent_pairs = [("0_20", "0-20% central"), ("50_80", "50-80% peripheral")]
    fig, axes = plt.subplots(2, 4, figsize=(17.7, 7.8), sharey="row")
    for ir, (cent, cent_label) in enumerate(cent_pairs):
        for ic, (var, label) in enumerate(variables):
            ax = axes[ir, ic]
            h_sig = fetch_ss_sum(d_sig, var, "isIsolated_isTight_isoR30_isSliding", cent, PT_BINS)
            h_bkg = fetch_ss_sum(d_bkg, var, "isIsolated_isTight_isoR30_isSliding", cent, PT_BINS)
            xs, ys = hist_to_step(h_sig)
            xb, yb = hist_to_step(h_bkg)
            if len(xs):
                ax.plot(xs, ys, color=COLORS["bdt"], lw=2.0, label="truth-matched photons")
            if len(xb):
                ax.plot(xb, yb, color=COLORS["bdt2"], lw=2.0, label="inclusive-jet candidates")
            ax.set_title(label if ir == 0 else "", fontsize=14, pad=10)
            ax.text(0.05, 0.90, cent_label, transform=ax.transAxes, fontsize=11.5, ha="left", va="top")
            set_style(ax, label if ir == 1 else None, "Normalized candidates" if ic == 0 else None)
            if ir == 0 and ic == 3:
                leg = ax.legend(loc="upper right", frameon=True, framealpha=0.93, fontsize=10.0)
                leg.get_frame().set_linewidth(0)
    fig.text(0.035, 0.975, "sPHENIX", ha="left", va="top", fontsize=14, fontweight="bold", fontstyle="italic")
    fig.text(0.105, 0.975, "Internal", ha="left", va="top", fontsize=14)
    fig.text(0.035, 0.938, r"Pythia overlay, BDT Region A; $15 < E_T^\gamma < 35$ GeV", ha="left", va="top", fontsize=11.5)
    fig.tight_layout(rect=(0, 0, 1, 0.90), w_pad=1.1, h_pad=1.0)
    f_sig.Close()
    f_bkg.Close()
    return save(fig, "money_shower_width_templates_regionA_central_vs_peripheral.png")


def plot_mlp_challenger() -> list[Path]:
    paths: list[Path] = []
    summary = MLP_REPORT / "slideReady/mlp_vs_bdt_base3x3_pt15to30_summary.csv"
    scan = MLP_REPORT / "slideReady/explicit_validation/mlp_validation_wp_scan_by_centrality_pt.csv"
    if not summary.exists() or not scan.exists():
        return paths

    rows = list(csv.DictReader(summary.open()))
    scopes = ["cent_0_20", "cent_20_50", "cent_50_80", "pt_15_20", "pt_20_25", "pt_25_30"]
    labels = ["0-20%", "20-50%", "50-80%", "15-20", "20-25", "25-30"]
    bdt = {r["scope"]: r for r in rows if r["model"].startswith("BDT")}
    mlp = {r["scope"]: r for r in rows if r["model"].startswith("MLP")}
    x = np.arange(len(scopes))

    fig, ax = plt.subplots(figsize=(11.0, 5.6))
    bdt_fake = [float(bdt[s]["wp80_fake_rate"]) for s in scopes]
    mlp_fake = [float(mlp[s]["wp80_fake_rate"]) for s in scopes]
    ax.bar(x - 0.18, bdt_fake, width=0.34, color=COLORS["bdt"], label="BDT validation")
    ax.bar(x + 0.18, mlp_fake, width=0.34, color=COLORS["mlp"], label="MLP validation")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylim(0.0, 0.62)
    set_style(ax, "Validation slice", "Background acceptance at 80% signal efficiency")
    header(ax, y=0.94, line2=r"Pythia overlay validation sample, $15 < E_T^\gamma < 30$ GeV")
    leg = ax.legend(loc="upper right", frameon=True, framealpha=0.93, fontsize=11.0)
    leg.get_frame().set_linewidth(0)
    fig.tight_layout()
    paths.append(save(fig, "money_mlp_vs_bdt_fake_rate_at_80sig_validation.png"))

    fig, ax = plt.subplots(figsize=(11.0, 5.6))
    bdt_auc = [float(bdt[s]["auc"]) for s in scopes]
    mlp_auc = [float(mlp[s]["auc"]) for s in scopes]
    ax.plot(x, bdt_auc, "o", color=COLORS["bdt"], markersize=7.0, label="BDT")
    ax.plot(x, mlp_auc, "s", color=COLORS["mlp"], markersize=7.0, label="MLP")
    for xi, ya, yb in zip(x, bdt_auc, mlp_auc):
        ax.plot([xi, xi], [ya, yb], color="#B0B0B0", lw=1.2, zorder=0)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylim(0.66, 0.85)
    set_style(ax, "Validation slice", "ROC AUC")
    header(ax, y=0.94, line2=r"Pythia overlay validation sample, $15 < E_T^\gamma < 30$ GeV")
    leg = ax.legend(loc="upper right", frameon=True, framealpha=0.93, fontsize=11.0)
    leg.get_frame().set_linewidth(0)
    fig.tight_layout()
    paths.append(save(fig, "money_mlp_vs_bdt_auc_validation.png"))
    return paths


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--frontier-source",
        choices=("existing-expanded-csv", "full-root"),
        default="existing-expanded-csv",
        help="Use existing expanded-family CSV for quick frontier, or scan all target80 ROOT outputs.",
    )
    args = parser.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    produced: list[Path] = []

    rows = build_frontier_rows_from_existing_expanded_csv()
    if args.frontier_source == "full-root" or not rows:
        rows = build_frontier_rows()
    write_csv(rows, OUT / "money_pareto_signal_recovery_vs_background.csv")
    produced.append(plot_pareto_frontier(rows))

    reference = deep.REFERENCE
    bdt_1535 = deep.width_variants("pt15to35", 15.0, 35.0)[3]

    abcd_rows = build_abcd_rows(reference, "Box-cuts")
    abcd_rows.extend(build_abcd_rows(bdt_1535, "15-35 base+3x3 BDT"))
    write_csv(abcd_rows, OUT / "money_abcd_trust_stack_values.csv")
    produced.extend(plot_abcd_stack(abcd_rows))

    produced.append(plot_food_chain(reference, bdt_1535))
    produced.append(plot_shower_shape_templates(bdt_1535))
    produced.extend(plot_mlp_challenger())

    print("[DONE] money plot mining outputs:")
    for path in produced:
        print(path)
    close_cached_files()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
