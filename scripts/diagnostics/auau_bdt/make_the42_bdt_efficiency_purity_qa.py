#!/usr/bin/env python3
"""Make THE-42 simulation-only BDT efficiency and shower-shape QA plots."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


THIS_FILE = Path(__file__).resolve()
REPO = next((p for p in THIS_FILE.parents if (p / "AGENTS.md").exists()), THIS_FILE.parents[3])
DEFAULT_BASE = REPO / "dataOutput/auauTightBDTValidation/THE42_wp80_centlinear_ss_20260604"
DEFAULT_SIGNAL = (
    DEFAULT_BASE
    / "sim_roots/preselectionNewPPG12_tightAuAuCentInputBase3x3BDT_nonTightAuAuBDTComplement_baseVariant"
    / "photonJet12and20merged_SIM/RecoilJets_embeddedPhoton12plus20_MERGED.root"
)
DEFAULT_INCLUSIVE = (
    DEFAULT_BASE
    / "sim_roots/preselectionNewPPG12_tightAuAuCentInputBase3x3BDT_nonTightAuAuBDTComplement_baseVariant"
    / "embeddedJet12and20and30and40merged_SIM/RecoilJets_embeddedJet12plus20plus30plus40_MERGED.root"
)
DEFAULT_WP_JSON = (
    REPO
    / "dataOutput/auauTightBDTValidation/THE8_branchA_ladder_scorecache_fullstat_20260527"
    / "the41_centdep_wp_slides_20260604/the41_centdep_bdt_wp_exact_thresholds.json"
)

VARS = ("weta", "wphi", "weta33", "wphi33", "weta35", "wphi53", "et1", "e11e33", "e32e35", "npbScore")
VAR_LABELS = {
    "weta": r"$w_{\eta}$",
    "wphi": r"$w_{\phi}$",
    "weta33": r"$w_{\eta,3\times3}$",
    "wphi33": r"$w_{\phi,3\times3}$",
    "weta35": r"$w_{\eta,3\times5}$",
    "wphi53": r"$w_{\phi,5\times3}$",
    "et1": r"$E_{1}/E_{\mathrm{cluster}}$",
    "e11e33": r"$e_{11}/e_{33}$",
    "e32e35": r"$e_{32}/e_{35}$",
    "npbScore": "NPB score",
}
PT_BINS = ((22, 24), (24, 26), (26, 28))
CENT_BINS = ((0, 10), (10, 20), (20, 30), (30, 40), (40, 50), (50, 60), (60, 80))
STAGE_TAGS = {
    "no preselection": ("inclusive_sig", "inclusive_bkg"),
    "NPB preselection": ("pre_sig", "pre_bkg"),
    "tight WP80 BDT": ("tight_sig", "tight_bkg"),
}


@dataclass
class HistCurve:
    edges: np.ndarray
    centers: np.ndarray
    values: np.ndarray
    errors: np.ndarray
    raw_integral: float
    raw_error: float


def open_root(path: Path):
    import ROOT

    ROOT.gROOT.SetBatch(True)
    ROOT.TH1.AddDirectory(False)
    f = ROOT.TFile.Open(str(path), "READ")
    if not f or f.IsZombie():
        raise OSError(f"Could not open ROOT file: {path}")
    return f


def get_sum_hist(root_file, var: str, tag: str, pt_bins=PT_BINS, cent_bins=CENT_BINS):
    acc = None
    missing: list[str] = []
    for pt_lo, pt_hi in pt_bins:
        for c_lo, c_hi in cent_bins:
            name = f"SIM/h_ss_{var}_{tag}_pT_{pt_lo}_{pt_hi}_cent_{c_lo}_{c_hi}"
            h = root_file.Get(name)
            if not h:
                missing.append(name)
                continue
            if acc is None:
                acc = h.Clone(f"{var}_{tag}_sum")
                acc.SetDirectory(0)
            else:
                acc.Add(h)
    if acc is None:
        raise RuntimeError(f"No histograms found for var={var}, tag={tag}; first missing={missing[:3]}")
    return acc, missing


def curve_from_hist(hist, x_min: float | None = None, x_max: float | None = None) -> HistCurve:
    nb = hist.GetNbinsX()
    edges = np.array([hist.GetXaxis().GetBinLowEdge(i) for i in range(1, nb + 2)], dtype=float)
    centers = np.array([hist.GetXaxis().GetBinCenter(i) for i in range(1, nb + 1)], dtype=float)
    counts = np.array([hist.GetBinContent(i) for i in range(1, nb + 1)], dtype=float)
    errs = np.array([hist.GetBinError(i) for i in range(1, nb + 1)], dtype=float)
    mask = np.ones_like(centers, dtype=bool)
    if x_min is not None:
        mask &= centers >= x_min
    if x_max is not None:
        mask &= centers < x_max
    centers = centers[mask]
    counts = counts[mask]
    errs = errs[mask]
    idx = np.flatnonzero(mask)
    edges = np.concatenate(([edges[idx[0]]], edges[idx + 1])) if len(idx) else edges
    integral = float(np.sum(counts))
    raw_error = float(np.sqrt(np.sum(errs * errs)))
    values = np.zeros_like(counts)
    errors = np.zeros_like(errs)
    if integral > 0:
        values = counts / integral
        errors = errs / integral
    return HistCurve(edges, centers, values, errors, integral, raw_error)


def tvd(a: np.ndarray, b: np.ndarray) -> float:
    if len(a) != len(b):
        raise ValueError("TVD inputs must have the same binning")
    return float(0.5 * np.sum(np.abs(a - b)))


def ratio(num: float, num_err: float, den: float, den_err: float) -> tuple[float, float]:
    if den <= 0:
        return math.nan, math.nan
    val = num / den
    rel2 = 0.0
    if num > 0:
        rel2 += (num_err / num) ** 2
    if den > 0:
        rel2 += (den_err / den) ** 2
    return val, abs(val) * math.sqrt(rel2)


def setup_style() -> None:
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
        "axes.edgecolor": "#111827",
        "axes.linewidth": 1.0,
        "xtick.color": "#111827",
        "ytick.color": "#111827",
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
    })


def draw_internal(ax) -> None:
    ax.text(0.035, 0.955, r"$\it{\bf{sPHENIX}}$ Internal", transform=ax.transAxes,
            ha="left", va="top", fontsize=11.5)
    ax.text(0.035, 0.885, r"Au+Au embedded, $22 \leq p_T^\gamma < 28$ GeV",
            transform=ax.transAxes, ha="left", va="top", fontsize=10.5)
    ax.text(0.035, 0.820, "0-80% centrality", transform=ax.transAxes,
            ha="left", va="top", fontsize=10.5)


def load_curves(signal_root: Path, inclusive_root: Path):
    signal_file = open_root(signal_root)
    inclusive_file = open_root(inclusive_root)
    curves: dict[str, dict[str, dict[str, HistCurve]]] = {}
    skipped: dict[str, str] = {}
    try:
        for var in VARS:
            var_curves = {}
            try:
                for stage, (sig_tag, bkg_tag) in STAGE_TAGS.items():
                    sig_hist, _ = get_sum_hist(signal_file, var, sig_tag)
                    bkg_hist, _ = get_sum_hist(inclusive_file, var, bkg_tag)
                    var_curves[stage] = {
                        "signal": curve_from_hist(sig_hist),
                        "inclusive": curve_from_hist(bkg_hist),
                    }
            except RuntimeError as exc:
                skipped[var] = str(exc)
                continue
            curves[var] = var_curves
    finally:
        signal_file.Close()
        inclusive_file.Close()
    if not curves:
        raise RuntimeError("No compatible shower-shape variables were found")
    return curves, skipped


def make_shift_summary(curves) -> list[dict[str, float | str]]:
    rows = []
    for var in curves:
        b0 = curves[var]["no preselection"]["inclusive"].values
        bp = curves[var]["NPB preselection"]["inclusive"].values
        bt = curves[var]["tight WP80 BDT"]["inclusive"].values
        st = curves[var]["tight WP80 BDT"]["signal"].values
        rows.append({
            "variable": var,
            "label": VAR_LABELS[var],
            "bkg_no_to_pre_tvd": tvd(b0, bp),
            "bkg_pre_to_tight_tvd": tvd(bp, bt),
            "bkg_no_to_tight_tvd": tvd(b0, bt),
            "tight_signal_vs_bkg_tvd": tvd(st, bt),
            "bkg_no_integral": curves[var]["no preselection"]["inclusive"].raw_integral,
            "bkg_tight_integral": curves[var]["tight WP80 BDT"]["inclusive"].raw_integral,
        })
    rows.sort(key=lambda r: float(r["bkg_no_to_tight_tvd"]), reverse=True)
    return rows


def render_shift_ranking(rows: list[dict[str, float | str]], outdir: Path) -> Path:
    labels = [str(r["label"]) for r in rows]
    y = np.arange(len(rows))
    no_tight = np.array([float(r["bkg_no_to_tight_tvd"]) for r in rows])
    pre_tight = np.array([float(r["bkg_pre_to_tight_tvd"]) for r in rows])
    sig_bkg = np.array([float(r["tight_signal_vs_bkg_tvd"]) for r in rows])

    fig, ax = plt.subplots(figsize=(12.8, 7.2), dpi=160)
    ax.barh(y + 0.24, no_tight, height=0.22, color="#2E77BB", label="Inclusive MC: no preselection -> tight")
    ax.barh(y, pre_tight, height=0.22, color="#1B9E77", label="Inclusive MC: NPB -> tight")
    ax.barh(y - 0.24, sig_bkg, height=0.22, color="#D33F49", label="Tight panel: Signal MC vs Inclusive MC")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=12)
    ax.invert_yaxis()
    ax.set_xlabel("total variation distance between unit-normalized shapes", fontsize=13)
    ax.set_title("Which shower-shape variables move most under the WP80 selection?", fontsize=19, fontweight="bold", pad=14)
    ax.grid(axis="x", color="#DDE3EA", linewidth=0.8)
    ax.legend(loc="lower right", fontsize=11.5, frameon=False)
    ax.text(0.01, 0.015, r"$22 \leq p_T^\gamma < 28$ GeV, 0-80% centrality; THE-42 simulation only",
            transform=ax.transAxes, fontsize=11.5, color="#4B5563")
    fig.tight_layout()
    path = outdir / "the42_shower_shape_variable_shift_ranking.png"
    fig.savefig(path)
    plt.close(fig)
    return path


def render_top_variable_grid(curves, rows: list[dict[str, float | str]], outdir: Path, n: int = 6) -> Path:
    top = [str(r["variable"]) for r in rows[:n]]
    fig, axes = plt.subplots(2, 3, figsize=(14.8, 8.4), dpi=160)
    colors = {
        "no preselection": "#6B7280",
        "NPB preselection": "#2364C8",
        "tight WP80 BDT": "#1B9E77",
    }
    for ax, var in zip(axes.flat, top):
        for stage in STAGE_TAGS:
            c = curves[var][stage]["inclusive"]
            ax.stairs(c.values, c.edges, color=colors[stage], linewidth=1.9, label=stage)
        ax.set_title(VAR_LABELS[var], fontsize=15.5, fontweight="bold")
        ax.set_ylabel("unit-normalized Inclusive MC", fontsize=10.5)
        ax.grid(True, color="#E5EAF0", linewidth=0.7)
        ax.tick_params(labelsize=10, direction="in", top=True, right=True)
    axes.flat[0].legend(loc="upper right", fontsize=10.5, frameon=False)
    fig.suptitle("Inclusive-MC shower-shape evolution for the most responsive variables",
                 fontsize=20, fontweight="bold", y=0.985)
    fig.text(0.055, 0.025,
             r"THE-42 simulation only; curves are unit-normalized after summing $22 \leq p_T^\gamma < 28$ GeV and 0-80% centrality.",
             fontsize=11.5, color="#4B5563")
    fig.tight_layout(rect=[0, 0.04, 1, 0.95])
    path = outdir / "the42_top_shower_shape_background_evolution_grid.png"
    fig.savefig(path)
    plt.close(fig)
    return path


def centrality_integrals(root_file, var: str, tag: str):
    out = []
    for cent in CENT_BINS:
        acc = None
        for pt in PT_BINS:
            name = f"SIM/h_ss_{var}_{tag}_pT_{pt[0]}_{pt[1]}_cent_{cent[0]}_{cent[1]}"
            h = root_file.Get(name)
            if not h:
                raise RuntimeError(f"Missing {name}")
            if acc is None:
                acc = h.Clone(f"{var}_{tag}_{cent[0]}_{cent[1]}")
                acc.SetDirectory(0)
            else:
                acc.Add(h)
        err2 = 0.0
        for i in range(1, acc.GetNbinsX() + 1):
            err2 += acc.GetBinError(i) ** 2
        out.append((float(acc.Integral()), math.sqrt(err2)))
    return out


def binom_err(eff: float, n: float) -> float:
    if not math.isfinite(eff) or n <= 0:
        return math.nan
    eff = min(max(eff, 0.0), 1.0)
    return math.sqrt(eff * (1.0 - eff) / n)


def purity_err(signal_pass: float, signal_pass_err: float, bkg_pass: float, bkg_pass_err: float) -> float:
    total = signal_pass + bkg_pass
    if total <= 0:
        return math.nan
    d_ds = bkg_pass / (total * total)
    d_db = -signal_pass / (total * total)
    return math.sqrt((d_ds * signal_pass_err) ** 2 + (d_db * bkg_pass_err) ** 2)


def interp_scan(scan_rows: list[dict], cent_label: str, threshold: float) -> tuple[float, float]:
    rows = [r for r in scan_rows if r["centrality_label"] == cent_label]
    rows.sort(key=lambda r: r["threshold"])
    xs = np.array([float(r["threshold"]) for r in rows], dtype=float)
    sig = np.array([float(r["signal_efficiency"]) for r in rows], dtype=float)
    bkg = np.array([float(r["background_fake_rate"]) for r in rows], dtype=float)
    threshold = float(np.clip(threshold, xs.min(), xs.max()))
    return float(np.interp(threshold, xs, sig)), float(np.interp(threshold, xs, bkg))


def compute_scorecache_wp_qa(wp_json: Path) -> tuple[list[dict], dict]:
    payload = json.loads(wp_json.read_text())
    row_by_key = {
        (r["wp_label"], r["centrality_label"]): r
        for r in payload["rows"]
    }
    scan = payload["efficiency_scan"]
    fit_rows = {}
    for r in payload["rows"]:
        label = r["wp_label"]
        fit_rows.setdefault(label, {
            "intercept": float(r.get("fit_intercept", math.nan)),
            "slope": float(r.get("fit_slope", math.nan)),
        })
    if any(not math.isfinite(v["intercept"]) for v in fit_rows.values()):
        # Fall back to the exact fitted values used in THE-41 if the compact JSON
        # predates those columns.
        summary_csv = wp_json.with_name("the41_centdep_bdt_wp_summary.csv")
        if summary_csv.exists():
            with summary_csv.open() as f:
                for row in csv.DictReader(f):
                    if row.get("row_type") == "flat_fit_constant":
                        fit_rows[row["wp_label"]] = {
                            "intercept": float(row["fit_intercept"]),
                            "slope": float(row["fit_slope"]),
                        }

    rows = []
    for wp_label in ("WP90", "WP80", "WP70"):
        target = float(wp_label.replace("WP", "")) / 100.0
        fit = fit_rows[wp_label]
        for c_lo, c_hi in CENT_BINS:
            c = 0.5 * (c_lo + c_hi)
            cent_label = f"{c_lo}-{c_hi}%"
            base = row_by_key[(wp_label, cent_label)]
            signal_n = float(base["signal_entries"])
            bkg_n = float(base["background_entries"])
            threshold = fit["intercept"] + fit["slope"] * c
            sig_eff, bkg_fake = interp_scan(scan, cent_label, threshold)
            sig_eff_err = binom_err(sig_eff, signal_n)
            bkg_fake_err = binom_err(bkg_fake, bkg_n)
            signal_pass = signal_n * sig_eff
            bkg_pass = bkg_n * bkg_fake
            signal_pass_err = signal_n * sig_eff_err
            bkg_pass_err = bkg_n * bkg_fake_err
            purity = signal_pass / (signal_pass + bkg_pass) if (signal_pass + bkg_pass) > 0 else math.nan
            purity_stat = purity_err(signal_pass, signal_pass_err, bkg_pass, bkg_pass_err)
            sb_before = signal_n / bkg_n if bkg_n > 0 else math.nan
            sb_after = signal_pass / bkg_pass if bkg_pass > 0 else math.nan
            rows.append({
                "wp_label": wp_label,
                "target_signal_efficiency": target,
                "cent_low": c_lo,
                "cent_high": c_hi,
                "cent_center": c,
                "cent_label": cent_label,
                "linear_threshold": threshold,
                "signal_entries": signal_n,
                "background_entries": bkg_n,
                "signal_efficiency": sig_eff,
                "signal_efficiency_stat_err": sig_eff_err,
                "background_fake_rate": bkg_fake,
                "background_fake_rate_stat_err": bkg_fake_err,
                "mc_purity_proxy": purity,
                "mc_purity_proxy_stat_err": purity_stat,
                "sb_before": sb_before,
                "sb_after": sb_after,
                "sb_enrichment": sb_after / sb_before if sb_before > 0 and math.isfinite(sb_after) else math.nan,
            })
    meta = {
        "wp_json": str(wp_json.resolve()),
        "score_key": payload.get("score_key"),
        "source_label": payload.get("source_label"),
        "training_sample": payload.get("training_sample"),
        "rows_loaded": payload.get("rows_loaded"),
    }
    return rows, meta


def render_centrality_qa(rows, outdir: Path) -> Path:
    colors = {"WP90": "#D95F9F", "WP80": "#1B9E77", "WP70": "#2E77BB"}
    labels = {"WP90": "WP90 context", "WP80": "WP80 nominal", "WP70": "WP70 context"}

    fig, axes = plt.subplots(2, 2, figsize=(13.8, 8.2), dpi=160)
    for wp_label in ("WP90", "WP80", "WP70"):
        sub = [r for r in rows if r["wp_label"] == wp_label]
        x = np.array([r["cent_center"] for r in sub])
        xerr = np.array([(r["cent_high"] - r["cent_low"]) / 2 for r in sub])
        lw = 2.8 if wp_label == "WP80" else 1.8
        alpha = 1.0 if wp_label == "WP80" else 0.72
        color = colors[wp_label]
        axes[0, 0].errorbar(x, [r["signal_efficiency"] for r in sub],
                            yerr=[r["signal_efficiency_stat_err"] for r in sub],
                            xerr=xerr, fmt="o-", linewidth=lw, color=color, alpha=alpha,
                            label=labels[wp_label])
        axes[0, 0].axhline(sub[0]["target_signal_efficiency"], color=color, alpha=0.18, linewidth=1.6)
        axes[0, 1].errorbar(x, [r["background_fake_rate"] for r in sub],
                            yerr=[r["background_fake_rate_stat_err"] for r in sub],
                            xerr=xerr, fmt="o-", linewidth=lw, color=color, alpha=alpha,
                            label=labels[wp_label])
        axes[1, 0].errorbar(x, [r["mc_purity_proxy"] for r in sub],
                            yerr=[r["mc_purity_proxy_stat_err"] for r in sub],
                            xerr=xerr, fmt="o-", linewidth=lw, color=color, alpha=alpha,
                            label=labels[wp_label])
        axes[1, 1].plot(x, [r["sb_enrichment"] for r in sub], "o-",
                        linewidth=lw, color=color, alpha=alpha, label=labels[wp_label])

    axes[0, 0].set_ylabel("signal efficiency", fontsize=12)
    axes[0, 0].set_ylim(0.64, 0.93)
    axes[0, 1].set_ylabel("inclusive-MC background acceptance", fontsize=12)
    axes[0, 1].set_yscale("log")
    axes[0, 1].legend(fontsize=9.8, frameon=False, loc="upper right")
    axes[1, 0].set_ylabel("MC purity proxy S/(S+B)", fontsize=12)
    axes[1, 0].set_xlabel("centrality percentile", fontsize=12)
    axes[1, 0].set_ylim(0, 1)
    axes[1, 1].axhline(1, color="#6B7280", linestyle=":", linewidth=1.0)
    axes[1, 1].set_ylabel("S/B enrichment vs before BDT cut", fontsize=12)
    axes[1, 1].set_xlabel("centrality percentile", fontsize=12)

    for ax in axes.flat:
        ax.set_xlim(0, 80)
        ax.set_xticks([0, 10, 20, 30, 40, 50, 60, 70, 80])
        ax.grid(True, color="#E5EAF0", linewidth=0.75)
        ax.tick_params(labelsize=10.5, direction="in", top=True, right=True)

    fig.suptitle("Centrality QA for the WP80 linear BDT working point", fontsize=19.5, fontweight="bold", y=0.985)
    fig.text(0.055, 0.020,
             "THE41 full-stat score-threshold scan evaluated at fitted linear cuts; binomial statistical errors shown. "
             "Purity is S/(S+B) from MC score-cache entries, not data-driven purity.",
             fontsize=10.2, color="#4B5563")
    fig.tight_layout(rect=[0, 0.070, 1, 0.94])
    path = outdir / "the42_wp80_efficiency_purity_centrality_qa.png"
    fig.savefig(path)
    plt.close(fig)
    return path


def write_csv(rows, path: Path) -> None:
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def build(args: argparse.Namespace) -> dict[str, str]:
    setup_style()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    curves, skipped_variables = load_curves(args.signal_root, args.inclusive_root)
    shift_rows = make_shift_summary(curves)
    qa_rows, qa_meta = compute_scorecache_wp_qa(args.wp_json)

    ranking_png = render_shift_ranking(shift_rows, args.output_dir)
    top_grid_png = render_top_variable_grid(curves, shift_rows, args.output_dir)
    qa_png = render_centrality_qa(qa_rows, args.output_dir)
    shift_csv = args.output_dir / "the42_shower_shape_shift_metrics.csv"
    qa_csv = args.output_dir / "the42_wp80_efficiency_purity_centrality_qa.csv"
    write_csv(shift_rows, shift_csv)
    write_csv(qa_rows, qa_csv)
    manifest = {
        "script": str(THIS_FILE),
        "signal_root": str(args.signal_root.resolve()),
        "inclusive_root": str(args.inclusive_root.resolve()),
        "pt_bins": PT_BINS,
        "centrality_bins": CENT_BINS,
        "variables_requested": VARS,
        "variables_plotted": tuple(curves.keys()),
        "variables_skipped": skipped_variables,
        "outputs": {
            "shift_ranking_png": str(ranking_png.resolve()),
            "top_variable_grid_png": str(top_grid_png.resolve()),
            "centrality_qa_png": str(qa_png.resolve()),
            "shift_metrics_csv": str(shift_csv.resolve()),
            "centrality_qa_csv": str(qa_csv.resolve()),
        },
        "centrality_qa_source": qa_meta,
        "definitions": {
            "shape_shift_metric": "total variation distance between unit-normalized histograms",
            "signal_efficiency": "THE41 score-threshold scan signal efficiency evaluated at the fitted linear WP threshold in each centrality bin",
            "background_fake_rate": "THE41 score-threshold scan Inclusive MC acceptance evaluated at the fitted linear WP threshold in each centrality bin",
            "purity_proxy": "score-cache Signal MC pass count / (Signal MC pass count + Inclusive MC pass count); simulation proxy, not data-driven purity",
        },
        "top_background_shift_variables": [
            {"variable": r["variable"], "bkg_no_to_tight_tvd": r["bkg_no_to_tight_tvd"]}
            for r in shift_rows[:6]
        ],
    }
    manifest_path = args.output_dir / "the42_efficiency_purity_qa_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return {k: str(v) for k, v in manifest["outputs"].items()} | {"manifest": str(manifest_path)}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--signal-root", type=Path, default=DEFAULT_SIGNAL)
    parser.add_argument("--inclusive-root", type=Path, default=DEFAULT_INCLUSIVE)
    parser.add_argument("--wp-json", type=Path, default=DEFAULT_WP_JSON)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_BASE / "efficiency_purity_qa")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    for path in (args.signal_root, args.inclusive_root, args.wp_json):
        if not path.exists():
            raise SystemExit(f"Missing input: {path}")
    outputs = build(args)
    for key, path in outputs.items():
        print(f"{key}={path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
