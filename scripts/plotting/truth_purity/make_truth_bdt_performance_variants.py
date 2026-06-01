#!/usr/bin/env python3
from __future__ import annotations
# Keep purpose-folder helpers runnable when invoked directly.
import sys as _codex_sys
from pathlib import Path as _CodexPath
_CODEX_THIS_FILE = _CodexPath(__file__).resolve()
_CODEX_SCRIPTS_DIR = next((p for p in _CODEX_THIS_FILE.parents if p.name == "scripts"), _CODEX_THIS_FILE.parent)
_CODEX_SCRIPTS_DIR_STR = str(_CODEX_SCRIPTS_DIR)
if _CODEX_SCRIPTS_DIR_STR not in _codex_sys.path:
    _codex_sys.path.append(_CODEX_SCRIPTS_DIR_STR)
del _CODEX_THIS_FILE, _CODEX_SCRIPTS_DIR, _CODEX_SCRIPTS_DIR_STR

from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch


ROOT = Path(
    "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/auauMLDiagnosticRuns/"
    "global_etcent_inclusive3_sixpack_20260516_135439"
)
VALID_BDT = ROOT / "validation" / "bdt_finished_only_20260516_180817"
VALID_BINNED = ROOT / "validation" / "bdt_binned_sidecars_fullstat_20260517_2152"
OUTDIR = ROOT / "slideReady" / "truth_bdt_performance_variants"

GLOBAL_HIST = VALID_BDT / "fine7x8_score_histograms_globalEtCent1535_bdt_iso_noIso.csv"
BINNED_HIST = VALID_BINNED / "coarse_centrality_score_histograms_binned_noIso.csv"
BASEV3E_FINE_HIST = ROOT / "validation" / "basev3e_controls_20260518_1110" / "basev3e_with_centrality_score_histograms.csv"
BASEV3E_DIAG = ROOT / "validation" / "basev3e_controls_20260518_1110" / "validation_deep_diagnostics.json"

CENT_ORDER = ["0_20", "20_50", "50_80"]
CENT_LABEL = {"0_20": "0-20%", "20_50": "20-50%", "50_80": "50-80%"}
FINE_CENT_ORDER = ["0-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-80"]
FINE_CENT_COLORS = ["#111827", "#374151", "#0072B2", "#009E73", "#D55E00", "#CC79A7", "#E69F00"]
ET3_ORDER = ["15-20", "20-25", "25-35"]
ET3_COLORS = {"15-20": "#0072B2", "20-25": "#009E73", "25-35": "#D55E00"}
ET8_ORDER = ["15-17", "17-19", "19-21", "21-23", "23-25", "25-27", "27-30", "30-35"]

MODELS = {
    "global32": {
        "label": "Global 32-feature BDT",
        "short": "global",
        "color": "#6B7280",
        "source_product": "globalEtCent1535_bdt_noIso",
    },
    "ptcent3": {
        "label": r"8 $p_T$ x 3 centrality BDT",
        "short": "8x3",
        "color": "#0072B2",
        "source_product": "globalEtCent1535_bdt_noIso_ptCent3",
    },
    "ptcent7": {
        "label": r"8 $p_T$ x 7 centrality BDT",
        "short": "8x7",
        "color": "#009E73",
        "source_product": "globalEtCent1535_bdt_noIso_ptCent7",
    },
}


def coarse_cent_from_label(label: str) -> str:
    low = int(str(label).split("-")[0])
    if low < 20:
        return "0_20"
    if low < 50:
        return "20_50"
    return "50_80"


def style_ax(ax, grid_axis: str = "both") -> None:
    ax.tick_params(direction="in", top=True, right=True, labelsize=10.5)
    ax.grid(True, axis=grid_axis, color="#D1D5DB", lw=0.7, alpha=0.75)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_linewidth(1.0)


def sphinx_label(fig, x=0.845, y=0.94, size=16) -> None:
    fig.text(x, y, "sPHENIX", ha="right", fontsize=size, fontstyle="italic", fontweight="bold")
    fig.text(x + 0.095, y, "Internal", ha="right", fontsize=size)


def load_global_hist() -> pd.DataFrame:
    df = pd.read_csv(GLOBAL_HIST)
    df = df[df["product"] == "globalEtCent1535_bdt_noIso"].copy()
    df["model"] = "global32"
    df["cent_bin"] = df["centrality_bin"].map(coarse_cent_from_label)
    g = (
        df.groupby(["model", "cent_bin", "score_low", "score_high"], as_index=False)[
            ["signal_count", "background_count"]
        ]
        .sum()
        .sort_values(["model", "cent_bin", "score_low"])
    )
    return g


def load_global_fine_hist() -> pd.DataFrame:
    df = pd.read_csv(GLOBAL_HIST)
    df = df[df["product"] == "globalEtCent1535_bdt_noIso"].copy()
    df["model"] = "global32"
    return df[
        [
            "model",
            "centrality_bin",
            "et_bin",
            "score_low",
            "score_high",
            "signal_count",
            "background_count",
        ]
    ].copy()


def load_basev3e_fine_hist() -> pd.DataFrame:
    df = pd.read_csv(BASEV3E_FINE_HIST)
    df = df[df["product"] == "baseBDT_v3E_withCentrality"].copy()
    df["model"] = "basev3e_cent"
    return df[
        [
            "model",
            "centrality_bin",
            "et_bin",
            "score_low",
            "score_high",
            "signal_count",
            "background_count",
        ]
    ].copy()


def load_binned_hist() -> pd.DataFrame:
    df = pd.read_csv(BINNED_HIST)
    rows = []
    for _, r in df.iterrows():
        model = None
        for key, meta in MODELS.items():
            if r["product"] == meta["source_product"]:
                model = key
                break
        if model is None:
            continue
        width = float(r["bin_hi"] - r["bin_lo"])
        rows.append(
            {
                "model": model,
                "cent_bin": r["cent_bin"],
                "score_low": float(r["bin_lo"]),
                "score_high": float(r["bin_hi"]),
                "signal_count": float(r["signal_density"]) * float(r["signal_entries"]) * width,
                "background_count": float(r["background_density"]) * float(r["background_entries"]) * width,
            }
        )
    return pd.DataFrame(rows)


def add_curves(hist: pd.DataFrame) -> pd.DataFrame:
    pieces = []
    for (model, cent), g in hist.groupby(["model", "cent_bin"], sort=False):
        g = g.sort_values("score_low").copy()
        sig = g["signal_count"].to_numpy(dtype=float)
        bkg = g["background_count"].to_numpy(dtype=float)
        sig_total = sig.sum()
        bkg_total = bkg.sum()
        sig_above = sig[::-1].cumsum()[::-1]
        bkg_above = bkg[::-1].cumsum()[::-1]
        total_above = sig_above + bkg_above
        g["score_center"] = 0.5 * (g["score_low"] + g["score_high"])
        g["truth_purity_bin"] = sig / np.where(sig + bkg > 0, sig + bkg, np.nan)
        g["signal_above"] = sig_above
        g["background_above"] = bkg_above
        g["total_above"] = total_above
        g["truth_purity_above"] = sig_above / np.where(total_above > 0, total_above, np.nan)
        g["signal_efficiency"] = sig_above / sig_total
        g["background_fake_rate"] = bkg_above / bkg_total
        g["background_rejection"] = 1.0 - g["background_fake_rate"]
        g["raw_truth_purity"] = sig_total / (sig_total + bkg_total)
        g["signal_entries_total"] = sig_total
        g["background_entries_total"] = bkg_total
        pieces.append(g)
    curves = pd.concat(pieces, ignore_index=True)
    curves = curves[curves["cent_bin"].isin(CENT_ORDER)].copy()
    return curves


def cumulative_curves(hist: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    pieces = []
    for keys, g in hist.groupby(group_cols, sort=False):
        if not isinstance(keys, tuple):
            keys = (keys,)
        g = g.sort_values("score_low").copy()
        sig = g["signal_count"].to_numpy(dtype=float)
        bkg = g["background_count"].to_numpy(dtype=float)
        sig_total = sig.sum()
        bkg_total = bkg.sum()
        sig_above = sig[::-1].cumsum()[::-1]
        bkg_above = bkg[::-1].cumsum()[::-1]
        total_above = sig_above + bkg_above
        out = g[group_cols + ["score_low", "score_high"]].copy()
        out["score_center"] = 0.5 * (g["score_low"].to_numpy() + g["score_high"].to_numpy())
        out["truth_purity_above"] = sig_above / np.where(total_above > 0, total_above, np.nan)
        out["signal_efficiency"] = sig_above / sig_total if sig_total > 0 else np.nan
        out["background_fake_rate"] = bkg_above / bkg_total if bkg_total > 0 else np.nan
        out["raw_truth_purity"] = sig_total / (sig_total + bkg_total) if (sig_total + bkg_total) > 0 else np.nan
        out["signal_above"] = sig_above
        out["background_above"] = bkg_above
        out["total_above"] = total_above
        out["signal_entries_total"] = sig_total
        out["background_entries_total"] = bkg_total
        pieces.append(out)
    return pd.concat(pieces, ignore_index=True)


def clean_turnon_tail(g: pd.DataFrame, min_signal_eff: float = 0.025, min_selected: float = 500.0) -> pd.DataFrame:
    """Drop low-stat extreme-score endpoints that produce non-physics visual spikes."""
    out = g.copy()
    if "signal_efficiency" in out:
        out = out[out["signal_efficiency"] >= min_signal_eff]
    if "total_above" in out:
        out = out[out["total_above"] >= min_selected]
    return out


def op_points(curves: pd.DataFrame, targets=(0.70, 0.80, 0.90)) -> pd.DataFrame:
    rows = []
    for (model, cent), g in curves.groupby(["model", "cent_bin"]):
        g = g.sort_values("signal_efficiency")
        for target in targets:
            idx = (g["signal_efficiency"] - target).abs().idxmin()
            r = g.loc[idx]
            rows.append(
                {
                    "model": model,
                    "model_label": MODELS[model]["label"],
                    "centrality": cent,
                    "centrality_label": CENT_LABEL[cent],
                    "target_signal_efficiency": target,
                    "score_threshold": r["score_low"],
                    "actual_signal_efficiency": r["signal_efficiency"],
                    "truth_purity": r["truth_purity_above"],
                    "background_fake_rate": r["background_fake_rate"],
                    "background_rejection": r["background_rejection"],
                    "raw_truth_purity": r["raw_truth_purity"],
                }
            )
    return pd.DataFrame(rows)


def plot_variant_a(curves: pd.DataFrame, ops: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15.6, 4.85), dpi=180, sharey=True)
    for ax, cent in zip(axes, CENT_ORDER):
        cent_curves = curves[curves["cent_bin"] == cent]
        for model in ["global32", "ptcent3", "ptcent7"]:
            meta = MODELS[model]
            g = cent_curves[cent_curves["model"] == model].sort_values("signal_efficiency")
            ax.plot(
                g["signal_efficiency"],
                g["truth_purity_above"],
                color=meta["color"],
                lw=2.5,
                label=meta["label"] if cent == "0_20" else None,
            )
            o = ops[(ops["model"] == model) & (ops["centrality"] == cent) & (ops["target_signal_efficiency"] == 0.80)]
            if len(o):
                ax.scatter(
                    o["actual_signal_efficiency"],
                    o["truth_purity"],
                    s=55,
                    color=meta["color"],
                    edgecolor="black",
                    linewidth=0.5,
                    zorder=5,
                )
        raw = cent_curves[cent_curves["model"] == "ptcent7"]["raw_truth_purity"].iloc[0]
        ax.axhline(raw, color="#7A7A7A", ls=":", lw=1.4)
        ax.text(0.965, raw + 0.008, "pre-BDT truth mix", ha="right", va="bottom", fontsize=8.8, color="#4B5563")
        ax.set_title(f"{CENT_LABEL[cent]} centrality", fontsize=13.5, fontweight="bold")
        ax.set_xlim(0.50, 1.0)
        ax.set_ylim(0.80, 1.005)
        ax.set_xlabel("Signal efficiency kept")
        if ax is axes[0]:
            ax.set_ylabel("Truth purity after BDT threshold")
            ax.legend(loc="lower left", frameon=False, fontsize=9.0)
        style_ax(ax)
    fig.text(0.055, 0.955, "Truth-labeled BDT performance: purity vs signal efficiency", fontsize=17.5, fontweight="bold")
    fig.text(0.055, 0.907, r"Photon12+20 truth signal vs Jet12+20+30 inclusive background; dots mark 80% signal-efficiency operating point", fontsize=11.7, color="#4B5563")
    sphinx_label(fig, x=0.855, y=0.948)
    fig.tight_layout(rect=(0.02, 0.02, 0.99, 0.88))
    fig.savefig(OUTDIR / "variant_a_truth_purity_vs_signal_efficiency_model_comparison_1x3.png")
    plt.close(fig)


def plot_variant_b(ops: pd.DataFrame) -> None:
    data = ops[ops["model"] == "ptcent7"].copy()
    fig, axes = plt.subplots(2, 3, figsize=(15.5, 7.0), dpi=180, sharex=True)
    targets = [0.70, 0.80, 0.90]
    colors = {0.70: "#0072B2", 0.80: "#009E73", 0.90: "#D55E00"}
    for col, cent in enumerate(CENT_ORDER):
        sub = data[data["centrality"] == cent]
        xs = np.arange(len(targets))
        top = axes[0, col]
        bot = axes[1, col]
        purity = [float(sub[sub["target_signal_efficiency"] == t]["truth_purity"].iloc[0]) for t in targets]
        fake = [float(sub[sub["target_signal_efficiency"] == t]["background_fake_rate"].iloc[0]) for t in targets]
        top.bar(xs, purity, color=[colors[t] for t in targets], width=0.62)
        bot.bar(xs, fake, color=[colors[t] for t in targets], width=0.62)
        raw = float(sub["raw_truth_purity"].iloc[0])
        top.axhline(raw, color="#6B7280", ls=":", lw=1.5)
        top.text(1.95, raw + 0.006, "pre-BDT truth mix", ha="right", va="bottom", fontsize=8.5, color="#4B5563")
        for ax, vals in [(top, purity), (bot, fake)]:
            for x, v in zip(xs, vals):
                ax.text(x, v + 0.012, f"{v:.2f}", ha="center", va="bottom", fontsize=10, fontweight="bold")
            ax.set_xticks(xs, [f"{int(t*100)}%" for t in targets])
            style_ax(ax, "y")
        top.set_ylim(0.80, 1.01)
        bot.set_ylim(0, 0.42)
        top.set_title(f"{CENT_LABEL[cent]} centrality", fontsize=13.5, fontweight="bold")
        if col == 0:
            top.set_ylabel("Truth purity")
            bot.set_ylabel("Background fake rate")
        bot.set_xlabel("Target signal efficiency")
    fig.text(0.055, 0.965, r"Truth-label operating points for the 8 $p_T$ x 7 centrality BDT", fontsize=17.5, fontweight="bold")
    fig.text(0.055, 0.918, r"At fixed signal efficiency, the selected candidate sample is substantially purer than the pre-BDT truth mixture", fontsize=11.7, color="#4B5563")
    sphinx_label(fig, x=0.855, y=0.957)
    fig.tight_layout(rect=(0.02, 0.02, 0.99, 0.89))
    fig.savefig(OUTDIR / "variant_b_ptcent7_truth_purity_fake_rate_operating_points_2x3.png")
    plt.close(fig)


def plot_variant_c(curves: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15.8, 4.9), dpi=180, sharey=True)
    for ax, cent in zip(axes, CENT_ORDER):
        g = curves[(curves["model"] == "ptcent7") & (curves["cent_bin"] == cent)].sort_values("score_low")
        ax.plot(g["score_low"], g["truth_purity_above"], color="#009E73", lw=2.8, label="truth purity")
        ax.plot(g["score_low"], g["signal_efficiency"], color="#0072B2", lw=2.2, ls="--", label="signal efficiency")
        ax.plot(g["score_low"], g["background_fake_rate"], color="#D55E00", lw=2.2, ls=":", label="background fake rate")
        for target in [0.70, 0.80, 0.90]:
            idx = (g["signal_efficiency"] - target).abs().idxmin()
            r = g.loc[idx]
            ax.scatter([r["score_low"]], [r["truth_purity_above"]], color="#009E73", edgecolor="black", s=42, zorder=5)
        ax.set_title(f"{CENT_LABEL[cent]} centrality", fontsize=13.5, fontweight="bold")
        ax.set_xlabel("Minimum BDT score")
        ax.set_xlim(0.0, 1.0)
        ax.set_ylim(0.0, 1.03)
        if ax is axes[0]:
            ax.set_ylabel("Fraction")
            ax.legend(loc="center left", frameon=False, fontsize=9.5)
        style_ax(ax)
    fig.text(0.055, 0.955, r"Truth-side threshold scan for the 8 $p_T$ x 7 centrality BDT", fontsize=17.5, fontweight="bold")
    fig.text(0.055, 0.907, "The same threshold axis simultaneously raises purity, lowers background fake rate, and trades signal efficiency.", fontsize=11.7, color="#4B5563")
    sphinx_label(fig, x=0.855, y=0.948)
    fig.tight_layout(rect=(0.02, 0.02, 0.99, 0.88))
    fig.savefig(OUTDIR / "variant_c_ptcent7_threshold_scan_truth_purity_eff_fake_1x3.png")
    plt.close(fig)


def plot_variant_d(curves: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15.8, 4.9), dpi=180, sharey=True)
    for ax, cent in zip(axes, CENT_ORDER):
        g = curves[(curves["model"] == "ptcent7") & (curves["cent_bin"] == cent)].sort_values("score_center")
        # Use every other bin for readability but preserve shape.
        g = g.iloc[::2].copy()
        x = g["score_center"].to_numpy()
        width = 0.034
        purity = g["truth_purity_bin"].to_numpy()
        ax.bar(x, 1.0 - purity, width=width, color="#D55E00", alpha=0.72, label="truth background" if cent == "0_20" else None)
        ax.bar(x, purity, bottom=1.0 - purity, width=width, color="#0072B2", alpha=0.82, label="truth signal" if cent == "0_20" else None)
        ax.plot(x, purity, color="#0F172A", lw=2.0)
        ax.set_title(f"{CENT_LABEL[cent]} centrality", fontsize=13.5, fontweight="bold")
        ax.set_xlabel("BDT score bin")
        ax.set_xlim(0.0, 1.0)
        ax.set_ylim(0.0, 1.02)
        if ax is axes[0]:
            ax.set_ylabel("Truth composition in score bin")
            ax.legend(loc="lower right", frameon=False, fontsize=9.5)
        style_ax(ax)
    fig.text(0.055, 0.955, "Truth composition by BDT score bin", fontsize=17.5, fontweight="bold")
    fig.text(0.055, 0.907, r"High-score bins are signal-dominated in every coarse centrality region for the 8 $p_T$ x 7 centrality BDT.", fontsize=11.7, color="#4B5563")
    sphinx_label(fig, x=0.855, y=0.948)
    fig.tight_layout(rect=(0.02, 0.02, 0.99, 0.88))
    fig.savefig(OUTDIR / "variant_d_ptcent7_truth_composition_by_score_bin_1x3.png")
    plt.close(fig)


def plot_variant_e(ops: pd.DataFrame) -> None:
    sub = ops[ops["target_signal_efficiency"] == 0.80].copy()
    fig, ax = plt.subplots(figsize=(13.6, 6.2), dpi=180)
    ax.axis("off")
    fig.text(0.055, 0.94, "Truth-label WP80 summary before reco/ABCD complications", fontsize=17.2, fontweight="bold")
    fig.text(0.055, 0.895, "Each cell uses the truth-labeled validation mixture only: selected signal / selected signal+background.", fontsize=11.8, color="#4B5563")
    sphinx_label(fig, x=0.855, y=0.935, size=14.5)

    x0, y0 = 0.07, 0.72
    col_w, row_h = 0.205, 0.145
    headers = ["Model", "0-20%", "20-50%", "50-80%", "Meaning"]
    widths = [0.25, 0.16, 0.16, 0.16, 0.25]
    xs = [x0]
    for w in widths[:-1]:
        xs.append(xs[-1] + w)

    def cell(x, y, w, h, text, face="#FFFFFF", color="#111827", weight="normal", size=11.0):
        patch = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.006,rounding_size=0.006", linewidth=0.8, edgecolor="#CBD5E1", facecolor=face, transform=fig.transFigure)
        fig.patches.append(patch)
        fig.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=size, color=color, fontweight=weight, linespacing=1.18)

    for i, h in enumerate(headers):
        cell(xs[i], y0, widths[i], 0.075, h, face="#E5E7EB", weight="bold", size=11.5)

    rows = [
        ("global32", "Global 32-feature BDT", "one model over all pT and centrality"),
        ("ptcent3", r"8 $p_T$ x 3 centrality BDT", "local pT bins and coarse centrality routing"),
        ("ptcent7", r"8 $p_T$ x 7 centrality BDT", "local pT bins and fine centrality routing"),
    ]
    for ir, (model, label, meaning) in enumerate(rows):
        y = y0 - (ir + 1) * row_h
        cell(xs[0], y, widths[0], row_h - 0.010, label, face="#F8FAFC", weight="bold", size=11.0)
        for ic, cent in enumerate(CENT_ORDER, start=1):
            r = sub[(sub["model"] == model) & (sub["centrality"] == cent)].iloc[0]
            text = f"purity {r.truth_purity:.3f}\nfake {r.background_fake_rate:.3f}\ncut {r.score_threshold:.2f}"
            face = "#ECFDF5" if model == "ptcent7" else "#F8FAFC"
            cell(xs[ic], y, widths[ic], row_h - 0.010, text, face=face, size=10.7)
        cell(xs[4], y, widths[4], row_h - 0.010, meaning, face="#F8FAFC", size=10.5)

    fig.text(0.075, 0.085, "Readout: the fine-routed BDT gives high truth purity at the nominal 80% signal-efficiency point, without relying on reco ABCD raw purity.", fontsize=11.4, color="#111827")
    fig.savefig(OUTDIR / "variant_e_wp80_truth_purity_fake_summary_table.png", bbox_inches="tight")
    plt.close(fig)


def plot_variant_f(curves: pd.DataFrame, ops: pd.DataFrame) -> None:
    model = "ptcent7"
    op80 = ops[(ops["model"] == model) & (ops["target_signal_efficiency"] == 0.80)]
    fig, axes = plt.subplots(1, 3, figsize=(15.2, 4.9), dpi=180)
    for ax, metric, ylabel, color, ylim in [
        (axes[0], "raw_truth_purity", "Pre-BDT truth mix", "#6B7280", (0.82, 0.88)),
        (axes[1], "truth_purity", "WP80 truth purity", "#009E73", (0.90, 0.98)),
        (axes[2], "background_rejection", "WP80 background rejection", "#D55E00", (0.68, 0.82)),
    ]:
        vals = [float(op80[op80["centrality"] == c][metric].iloc[0]) for c in CENT_ORDER]
        xs = np.arange(len(CENT_ORDER))
        ax.plot(xs, vals, marker="o", ms=7.5, lw=2.6, color=color)
        for x, v in zip(xs, vals):
            ax.text(x, v + 0.006, f"{v:.3f}", ha="center", va="bottom", fontsize=10, fontweight="bold")
        ax.set_xticks(xs, [CENT_LABEL[c] for c in CENT_ORDER])
        ax.set_ylabel(ylabel)
        ax.set_ylim(*ylim)
        style_ax(ax, "y")
    fig.text(0.055, 0.955, r"Clean truth-side readout for the 8 $p_T$ x 7 centrality BDT", fontsize=17.5, fontweight="bold")
    fig.text(0.055, 0.907, "The validation truth labels show the expected behavior: WP80 keeps signal while selecting a substantially purer sample.", fontsize=11.7, color="#4B5563")
    sphinx_label(fig, x=0.855, y=0.948)
    fig.tight_layout(rect=(0.02, 0.02, 0.99, 0.88))
    fig.savefig(OUTDIR / "variant_f_ptcent7_truth_side_readout_1x3.png")
    plt.close(fig)


def plot_turnon_coarse_overlay(curves: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(15.2, 5.2), dpi=180, sharey=True)
    products = [("global32", "Global 32-feature BDT"), ("ptcent7", r"8 $p_T$ x 7 centrality BDT")]
    colors = {"0_20": "#111827", "20_50": "#0072B2", "50_80": "#D55E00"}
    for ax, (model, title) in zip(axes, products):
        for cent in CENT_ORDER:
            g = curves[(curves["model"] == model) & (curves["cent_bin"] == cent)].sort_values("score_low")
            g = clean_turnon_tail(g)
            ax.plot(g["score_low"], g["truth_purity_above"], lw=2.7, color=colors[cent], label=CENT_LABEL[cent])
            raw = float(g["raw_truth_purity"].iloc[0])
            ax.axhline(raw, color=colors[cent], lw=1.1, ls=":", alpha=0.55)
        ax.set_title(title, fontsize=13.5, fontweight="bold")
        ax.set_xlabel("Minimum BDT score")
        ax.set_xlim(0, 1)
        ax.set_ylim(0.82, 1.005)
        style_ax(ax)
    axes[0].set_ylabel("Truth purity above threshold")
    axes[0].legend(frameon=False, fontsize=10, loc="lower right", title="Centrality")
    fig.text(0.055, 0.955, "Truth-purity turn-on curves: coarse centrality overlay", fontsize=17.5, fontweight="bold")
    fig.text(0.055, 0.907, "Dotted lines show the pre-BDT truth mixture; solid curves show selected purity after requiring score > threshold.", fontsize=11.7, color="#4B5563")
    sphinx_label(fig, x=0.855, y=0.948)
    fig.tight_layout(rect=(0.02, 0.02, 0.99, 0.88))
    fig.savefig(OUTDIR / "variant_g_truth_purity_turnon_coarse_cent_overlay_global_vs_ptcent7.png")
    plt.close(fig)


def plot_turnon_fine_cent_integrated(global_fine: pd.DataFrame, base_fine: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(15.5, 5.3), dpi=180, sharey=True)
    specs = [(global_fine, "global32", "Global 32-feature BDT"), (base_fine, "basev3e_cent", "Base v3E + centrality")]
    for ax, (hist, model, title) in zip(axes, specs):
        curves = cumulative_curves(hist, ["model", "centrality_bin"])
        for cent, color in zip(FINE_CENT_ORDER, FINE_CENT_COLORS):
            g = curves[(curves["model"] == model) & (curves["centrality_bin"] == cent)].sort_values("score_low")
            if len(g) == 0:
                continue
            g = clean_turnon_tail(g)
            ax.plot(g["score_low"], g["truth_purity_above"], lw=2.0, color=color, label=f"{cent}%")
        ax.set_title(title, fontsize=13.5, fontweight="bold")
        ax.set_xlabel("Minimum BDT score")
        ax.set_xlim(0, 1)
        ax.set_ylim(0.78, 1.005)
        style_ax(ax)
    axes[0].set_ylabel("Truth purity above threshold")
    axes[1].legend(frameon=False, fontsize=8.4, loc="lower right", ncol=1, title="Fine centrality")
    fig.text(0.055, 0.955, "Fine-centrality truth-purity turn-on, integrated over photon $E_T$", fontsize=17.5, fontweight="bold")
    fig.text(0.055, 0.907, "This is the clean truth-side behavior: high BDT thresholds select signal-dominated candidates across all fine centrality bins.", fontsize=11.7, color="#4B5563")
    sphinx_label(fig, x=0.855, y=0.948)
    fig.tight_layout(rect=(0.02, 0.02, 0.99, 0.88))
    fig.savefig(OUTDIR / "variant_h_truth_purity_turnon_fine_cent_overlay_integrated_et.png")
    plt.close(fig)


def plot_turnon_et_overlay(global_fine: pd.DataFrame, base_fine: pd.DataFrame) -> None:
    curves_global = cumulative_curves(global_fine, ["model", "et_bin"])
    curves_base = cumulative_curves(base_fine, ["model", "et_bin"])
    fig, axes = plt.subplots(1, 2, figsize=(15.5, 5.3), dpi=180, sharey=True)
    cmap = plt.get_cmap("viridis")
    for ax, curves, model, et_order, title in [
        (axes[0], curves_global, "global32", ET8_ORDER, "Global 32-feature BDT: 8 $E_T$ bins"),
        (axes[1], curves_base, "basev3e_cent", ET3_ORDER, "Base v3E + centrality: 3 $E_T$ bins"),
    ]:
        for i, et in enumerate(et_order):
            g = curves[(curves["model"] == model) & (curves["et_bin"] == et)].sort_values("score_low")
            if len(g) == 0:
                continue
            g = clean_turnon_tail(g)
            color = cmap(0.08 + 0.86 * i / max(len(et_order) - 1, 1))
            ax.plot(g["score_low"], g["truth_purity_above"], lw=2.1, color=color, label=f"{et} GeV")
        ax.set_title(title, fontsize=13.5, fontweight="bold")
        ax.set_xlabel("Minimum BDT score")
        ax.set_xlim(0, 1)
        ax.set_ylim(0.78, 1.005)
        ax.legend(frameon=False, fontsize=8.0, loc="lower right", ncol=1)
        style_ax(ax)
    axes[0].set_ylabel("Truth purity above threshold")
    fig.text(0.055, 0.955, "$E_T$-resolved truth-purity turn-on, integrated over centrality", fontsize=17.5, fontweight="bold")
    fig.text(0.055, 0.907, "The turn-on is visible within each photon-energy slice, not only after integrating all candidates.", fontsize=11.7, color="#4B5563")
    sphinx_label(fig, x=0.855, y=0.948)
    fig.tight_layout(rect=(0.02, 0.02, 0.99, 0.88))
    fig.savefig(OUTDIR / "variant_i_truth_purity_turnon_et_overlay_integrated_centrality.png")
    plt.close(fig)


def plot_turnon_coarse_by_et_grid(base_fine: pd.DataFrame) -> None:
    tmp = base_fine.copy()
    tmp["coarse_cent"] = tmp["centrality_bin"].map(lambda x: coarse_cent_from_label(str(x).replace("_", "-")))
    curves = cumulative_curves(tmp, ["model", "coarse_cent", "et_bin"])
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 5.0), dpi=180, sharey=True)
    for ax, cent in zip(axes, CENT_ORDER):
        for et in ET3_ORDER:
            g = curves[(curves["model"] == "basev3e_cent") & (curves["coarse_cent"] == cent) & (curves["et_bin"] == et)].sort_values("score_low")
            if len(g) == 0:
                continue
            g = clean_turnon_tail(g)
            ax.plot(g["score_low"], g["truth_purity_above"], color=ET3_COLORS[et], lw=2.5, label=f"{et} GeV")
        ax.set_title(f"{CENT_LABEL[cent]} centrality", fontsize=13.5, fontweight="bold")
        ax.set_xlabel("Minimum BDT score")
        ax.set_xlim(0, 1)
        ax.set_ylim(0.75, 1.005)
        style_ax(ax)
    axes[0].set_ylabel("Truth purity above threshold")
    axes[0].legend(frameon=False, fontsize=9.6, loc="lower right", title="$E_T$ bin")
    fig.text(0.055, 0.955, "Base v3E + centrality purity turn-on by coarse centrality and $E_T$", fontsize=17.5, fontweight="bold")
    fig.text(0.055, 0.907, "Three photon-energy slices per centrality bin show the same monotonic truth-purity turn-on.", fontsize=11.7, color="#4B5563")
    sphinx_label(fig, x=0.855, y=0.948)
    fig.tight_layout(rect=(0.02, 0.02, 0.99, 0.88))
    fig.savefig(OUTDIR / "variant_j_basev3e_cent_truth_purity_turnon_coarse_cent_by_et_1x3.png")
    plt.close(fig)


def plot_integrated_centrality_feature_points() -> None:
    import json

    diag = json.load(open(BASEV3E_DIAG))
    products = [
        ("baseBDT_v3E_withOutCentraltiy", "base v3E", "#6B7280"),
        ("baseBDT_v3E_withCentrality", "base v3E + centrality", "#009E73"),
        ("baseBDT_v3E_withOutCentraltiy_w33", "base v3E + 3x3 widths", "#0072B2"),
        ("baseBDT_v3E_withCentrality_w33", "base v3E + centrality + 3x3 widths", "#D55E00"),
    ]
    total_sig = diag["products"]["baseBDT_v3E_withCentrality"]["score_summary"]["signal"]["entries"]
    total_bkg = diag["products"]["baseBDT_v3E_withCentrality"]["score_summary"]["background"]["entries"]
    rows = []
    for product, label, color in products:
        for item in diag["products"][product]["thresholds_inclusive"]["by_signal_efficiency"]:
            eff = float(item["signal_efficiency"])
            fake = float(item["background_fake_rate"])
            purity = eff * total_sig / (eff * total_sig + fake * total_bkg)
            rows.append(
                {
                    "product": product,
                    "label": label,
                    "color": color,
                    "target_signal_efficiency": item["target_signal_efficiency"],
                    "actual_signal_efficiency": eff,
                    "background_fake_rate": fake,
                    "truth_purity": purity,
                    "threshold": item["threshold"],
                }
            )
    df = pd.DataFrame(rows)
    df.to_csv(OUTDIR / "basev3e_cent_vs_no_cent_inclusive_truth_purity_points.csv", index=False)

    fig, ax = plt.subplots(figsize=(12.6, 6.4), dpi=180)
    for product, label, color in products:
        g = df[df["product"] == product].sort_values("actual_signal_efficiency")
        ax.plot(g["actual_signal_efficiency"], g["truth_purity"], marker="o", lw=2.4, ms=6.5, color=color, label=label)
    raw = total_sig / (total_sig + total_bkg)
    ax.axhline(raw, color="#111827", ls=":", lw=1.3)
    ax.text(0.515, raw + 0.004, "pre-BDT truth mix", ha="left", va="bottom", fontsize=9.0)
    ax.set_xlabel("Signal efficiency kept")
    ax.set_ylabel("Truth purity after threshold")
    ax.set_xlim(0.46, 0.97)
    ax.set_ylim(0.86, 0.965)
    ax.legend(frameon=False, fontsize=9.2, loc="lower left")
    style_ax(ax)
    fig.text(0.06, 0.955, "Inclusive truth purity: centrality feature vs no centrality", fontsize=16.2, fontweight="bold")
    fig.text(
        0.06,
        0.908,
        "Exact inclusive validation operating points. Per-bin no-centrality purity curves need score histograms that are not present locally.",
        fontsize=10.4,
        color="#4B5563",
    )
    sphinx_label(fig, x=0.88, y=0.948, size=14.0)
    fig.tight_layout(rect=(0.03, 0.02, 0.99, 0.87))
    fig.savefig(OUTDIR / "variant_k_basev3e_cent_vs_no_cent_inclusive_truth_purity_points.png")
    plt.close(fig)


def plot_integrated_et_zoom(curves: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(11.4, 6.2), dpi=180)
    for model in ["global32", "ptcent3", "ptcent7"]:
        meta = MODELS[model]
        g = cumulative_curves(
            curves[curves["model"] == model][["model", "score_low", "score_high", "signal_count", "background_count"]],
            ["model"],
        )
        g = clean_turnon_tail(g, min_signal_eff=0.015, min_selected=1000.0)
        ax.plot(g["score_low"], g["truth_purity_above"], color=meta["color"], lw=2.8, label=meta["label"])
    ax.set_xlabel("Minimum BDT score")
    ax.set_ylabel("Truth purity above threshold")
    ax.set_xlim(0.25, 1.0)
    ax.set_ylim(0.86, 1.005)
    ax.legend(frameon=False, fontsize=10, loc="lower right")
    style_ax(ax)
    fig.text(0.06, 0.955, "Integrated truth-purity turn-on across all centralities and $E_T$", fontsize=15.8, fontweight="bold")
    fig.text(0.06, 0.908, "Clean single-panel readout: increasing BDT score monotonically enriches true photons.", fontsize=10.5, color="#4B5563")
    sphinx_label(fig, x=0.88, y=0.948, size=14.0)
    fig.tight_layout(rect=(0.03, 0.02, 0.99, 0.88))
    fig.savefig(OUTDIR / "variant_l_truth_purity_turnon_integrated_all_cent_et_model_overlay.png")
    plt.close(fig)


def main() -> int:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    hist = pd.concat([load_global_hist(), load_binned_hist()], ignore_index=True)
    curves = add_curves(hist)
    ops = op_points(curves)
    global_fine = load_global_fine_hist()
    base_fine = load_basev3e_fine_hist()
    curves.to_csv(OUTDIR / "truth_bdt_score_threshold_curves.csv", index=False)
    ops.to_csv(OUTDIR / "truth_bdt_target_efficiency_operating_points.csv", index=False)

    plot_variant_a(curves, ops)
    plot_variant_b(ops)
    plot_variant_c(curves)
    plot_variant_d(curves)
    plot_variant_e(ops)
    plot_variant_f(curves, ops)
    plot_turnon_coarse_overlay(curves)
    plot_turnon_fine_cent_integrated(global_fine, base_fine)
    plot_turnon_et_overlay(global_fine, base_fine)
    plot_turnon_coarse_by_et_grid(base_fine)
    plot_integrated_centrality_feature_points()
    plot_integrated_et_zoom(hist)

    print(OUTDIR)
    for p in sorted(OUTDIR.glob("variant_*.png")):
        print(p)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
