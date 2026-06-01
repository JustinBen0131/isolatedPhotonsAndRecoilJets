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

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = Path(
    "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/auauMLDiagnosticRuns/"
    "global_etcent_inclusive3_sixpack_20260516_135439"
)
OUTDIR = ROOT / "slideReady" / "truth_bdt_performance_variants"

BASEV3E_HIST = ROOT / "validation" / "basev3e_controls_20260518_1110" / "basev3e_with_centrality_score_histograms.csv"
REFERENCE_EFF = Path(
    "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/target80_first_offline/"
    "bdt_target80_gated_20260512_001012/analysis_config_etfine_15to35_target80/"
    "efficiency_qa/truth_photon_recovery_box_cuts_points.csv"
)
REFERENCE_PUR = (
    ROOT
    / "slideReady"
    / "current_best_target80_recoiljets"
    / "reference_vs_ptcent7bdt_target80_reco_truth_raw_abcd_purity_pt15to35_points.csv"
)
BASEV3E_TARGET80_PUR_8 = (
    ROOT
    / "slideReady"
    / "basev3e_reco_efficiency_overlay"
    / "basev3e_target80_raw_sim_purity_recoPt15to35_points.csv"
)
BASEV3E_TARGET80_EFF_8 = (
    ROOT
    / "slideReady"
    / "basev3e_reco_efficiency_overlay"
    / "basev3e_target80_vs_reference_boxcuts_reco_efficiency_pt15to35_points.csv"
)

CENT_ORDER = ["0_20", "20_50", "50_80"]
CENT_LABEL = {"0_20": "0-20%", "20_50": "20-50%", "50_80": "50-80%"}
ET_GROUPS = {
    "15-20": {"label": "15-20", "center": 17.5, "low": 15.0, "high": 20.0},
    "20-25": {"label": "20-25", "center": 22.5, "low": 20.0, "high": 25.0},
    "25-35": {"label": "25-35", "center": 30.0, "low": 25.0, "high": 35.0},
}

REFERENCE_EFF_BINS = {
    16.0: (15.0, 17.0),
    18.0: (17.0, 19.0),
    20.0: (19.0, 21.0),
    22.0: (21.0, 23.0),
    24.0: (23.0, 25.0),
    26.0: (25.0, 27.0),
    28.5: (27.0, 30.0),
    32.5: (30.0, 35.0),
}
ET8_GROUPS = {
    "15-17": {"label": "15-17", "center": 16.0, "low": 15.0, "high": 17.0},
    "17-19": {"label": "17-19", "center": 18.0, "low": 17.0, "high": 19.0},
    "19-21": {"label": "19-21", "center": 20.0, "low": 19.0, "high": 21.0},
    "21-23": {"label": "21-23", "center": 22.0, "low": 21.0, "high": 23.0},
    "23-25": {"label": "23-25", "center": 24.0, "low": 23.0, "high": 25.0},
    "25-27": {"label": "25-27", "center": 26.0, "low": 25.0, "high": 27.0},
    "27-30": {"label": "27-30", "center": 28.5, "low": 27.0, "high": 30.0},
    "30-35": {"label": "30-35", "center": 32.5, "low": 30.0, "high": 35.0},
}


def style_ax(ax, grid_axis: str = "both") -> None:
    ax.tick_params(direction="in", top=True, right=True, labelsize=10.0)
    ax.grid(True, axis=grid_axis, color="#D1D5DB", lw=0.7, alpha=0.75)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_linewidth(1.0)


def sphinx_label(fig, x=0.865, y=0.958, size=14.5) -> None:
    fig.text(x, y, "sPHENIX", ha="right", fontsize=size, fontstyle="italic", fontweight="bold")
    fig.text(x + 0.095, y, "Internal", ha="right", fontsize=size)


def coarse_cent(label: str) -> str:
    low = int(str(label).split("-")[0])
    if low < 20:
        return "0_20"
    if low < 50:
        return "20_50"
    return "50_80"


def overlap_fraction(low: float, high: float, target_low: float, target_high: float) -> float:
    overlap = max(0.0, min(high, target_high) - max(low, target_low))
    width = max(high - low, 1e-12)
    return overlap / width


def cumulative_curve(hist: pd.DataFrame) -> pd.DataFrame:
    g = hist.sort_values("score_low").copy()
    sig = g["signal_count"].to_numpy(float)
    bkg = g["background_count"].to_numpy(float)
    sig_total = sig.sum()
    bkg_total = bkg.sum()
    sig_above = sig[::-1].cumsum()[::-1]
    bkg_above = bkg[::-1].cumsum()[::-1]
    total_above = sig_above + bkg_above
    out = g[["score_low", "score_high"]].copy()
    out["signal_efficiency"] = sig_above / sig_total if sig_total > 0 else np.nan
    out["truth_purity"] = sig_above / np.where(total_above > 0, total_above, np.nan)
    out["background_fake_rate"] = bkg_above / bkg_total if bkg_total > 0 else np.nan
    return out


def interp_at_eff(curve: pd.DataFrame, target_eff: float) -> dict[str, float]:
    curve = curve.replace([np.inf, -np.inf], np.nan).dropna(subset=["signal_efficiency", "truth_purity"])
    curve = curve.sort_values("signal_efficiency").copy()
    x = curve["signal_efficiency"].to_numpy(float)
    keep = np.r_[True, np.diff(x) > 1e-12]
    curve = curve.iloc[keep].copy()
    x = curve["signal_efficiency"].to_numpy(float)
    target = float(np.clip(target_eff, np.nanmin(x), np.nanmax(x)))
    return {
        "efficiency": target,
        "truth_purity": float(np.interp(target, x, curve["truth_purity"].to_numpy(float))),
        "score_threshold": float(np.interp(target, x, curve["score_low"].to_numpy(float))),
        "background_fake_rate": float(np.interp(target, x, curve["background_fake_rate"].to_numpy(float))),
    }


def load_reference_efficiency() -> pd.DataFrame:
    df = pd.read_csv(REFERENCE_EFF)
    df = df[df["label"].str.contains("Box", case=False, na=False)].copy()
    df["centrality"] = df["cent"].astype(str)
    df["pt_low"] = df["pt_mid"].map(lambda x: REFERENCE_EFF_BINS.get(float(x), (np.nan, np.nan))[0])
    df["pt_high"] = df["pt_mid"].map(lambda x: REFERENCE_EFF_BINS.get(float(x), (np.nan, np.nan))[1])
    rows = []
    for cent in CENT_ORDER:
        for et, meta in ET_GROUPS.items():
            g = df[df["centrality"] == cent].copy()
            g["overlap_fraction"] = g.apply(
                lambda row: overlap_fraction(row["pt_low"], row["pt_high"], meta["low"], meta["high"]), axis=1
            )
            g = g[g["overlap_fraction"] > 0].copy()
            if len(g) == 0:
                continue
            # Prefer inverse-variance weighting for the published points, with
            # fractional overlap for bins that straddle the coarse ET boundary.
            w = g["overlap_fraction"].to_numpy(float) / np.maximum(g["error"].to_numpy(float), 1e-9) ** 2
            rows.append(
                {
                    "centrality": cent,
                    "et_bin": et,
                    "pt_center": meta["center"],
                    "reference_efficiency": float(np.average(g["value"], weights=w)),
                    "reference_efficiency_error": float(np.sqrt(1.0 / w.sum())),
                }
            )
    return pd.DataFrame(rows)


def load_reference_purity() -> pd.DataFrame:
    df = pd.read_csv(REFERENCE_PUR)
    df = df[(df["model"] == "Reference box cuts") & (df["quantity"] == "truth_labeled_A_purity")].copy()
    df["centrality"] = df["centrality"].astype(str)
    rows = []
    for cent in CENT_ORDER:
        for et, meta in ET_GROUPS.items():
            g = df[df["centrality"] == cent].copy()
            g["overlap_fraction"] = g.apply(
                lambda row: overlap_fraction(row["pt_low"], row["pt_high"], meta["low"], meta["high"]), axis=1
            )
            g = g[g["overlap_fraction"] > 0].copy()
            if len(g) == 0:
                continue
            sig = (g["signal_tight"] * g["overlap_fraction"]).sum()
            bkg = (g["background_tight"] * g["overlap_fraction"]).sum()
            purity = sig / (sig + bkg)
            err = np.sqrt(purity * (1.0 - purity) / max(sig + bkg, 1.0))
            rows.append(
                {
                    "centrality": cent,
                    "et_bin": et,
                    "pt_center": meta["center"],
                    "reference_truth_purity": float(purity),
                    "reference_truth_purity_error": float(err),
                    "reference_signal_tight": float(sig),
                    "reference_background_tight": float(bkg),
                }
            )
    return pd.DataFrame(rows)


def load_reference_purity_for_bins(bin_map: dict[str, dict]) -> pd.DataFrame:
    df = pd.read_csv(REFERENCE_PUR)
    df = df[(df["model"] == "Reference box cuts") & (df["quantity"] == "truth_labeled_A_purity")].copy()
    df["centrality"] = df["centrality"].astype(str)
    rows = []
    for cent in CENT_ORDER:
        for et, meta in bin_map.items():
            g = df[df["centrality"] == cent].copy()
            g["overlap_fraction"] = g.apply(
                lambda row: overlap_fraction(row["pt_low"], row["pt_high"], meta["low"], meta["high"]), axis=1
            )
            g = g[g["overlap_fraction"] > 0].copy()
            if len(g) == 0:
                continue
            sig = (g["signal_tight"] * g["overlap_fraction"]).sum()
            bkg = (g["background_tight"] * g["overlap_fraction"]).sum()
            purity = sig / (sig + bkg)
            err = np.sqrt(purity * (1.0 - purity) / max(sig + bkg, 1.0))
            rows.append(
                {
                    "centrality": cent,
                    "et_bin": et,
                    "pt_center": meta["center"],
                    "reference_truth_purity": float(purity),
                    "reference_truth_purity_error": float(err),
                    "reference_signal_tight": float(sig),
                    "reference_background_tight": float(bkg),
                }
            )
    return pd.DataFrame(rows)


def load_basev3e_matched_points(ref_eff: pd.DataFrame) -> pd.DataFrame:
    hist = pd.read_csv(BASEV3E_HIST)
    hist = hist[hist["product"] == "baseBDT_v3E_withCentrality"].copy()
    hist["centrality"] = hist["centrality_bin"].map(coarse_cent)
    hist = (
        hist.groupby(["centrality", "et_bin", "score_low", "score_high"], as_index=False)[
            ["signal_count", "background_count"]
        ]
        .sum()
        .sort_values(["centrality", "et_bin", "score_low"])
    )
    rows = []
    for _, target in ref_eff.iterrows():
        h = hist[(hist["centrality"] == target["centrality"]) & (hist["et_bin"] == target["et_bin"])]
        if len(h) == 0:
            continue
        curve = cumulative_curve(h)
        r = interp_at_eff(curve, float(target["reference_efficiency"]))
        rows.append(
            {
                "centrality": target["centrality"],
                "et_bin": target["et_bin"],
                "pt_center": target["pt_center"],
                "bdt_truth_purity_at_reference_eff": r["truth_purity"],
                "bdt_efficiency": r["efficiency"],
                "bdt_score_threshold": r["score_threshold"],
                "bdt_background_fake_rate": r["background_fake_rate"],
            }
        )
    return pd.DataFrame(rows)


def build_comparison() -> pd.DataFrame:
    ref_eff = load_reference_efficiency()
    ref_pur = load_reference_purity()
    bdt = load_basev3e_matched_points(ref_eff)
    df = ref_eff.merge(ref_pur, on=["centrality", "et_bin", "pt_center"]).merge(
        bdt, on=["centrality", "et_bin", "pt_center"]
    )
    df["centrality_label"] = df["centrality"].map(CENT_LABEL)
    df["purity_gain_absolute"] = df["bdt_truth_purity_at_reference_eff"] - df["reference_truth_purity"]
    df["purity_gain_relative_percent"] = 100.0 * df["purity_gain_absolute"] / df["reference_truth_purity"]
    df["efficiency_delta"] = df["bdt_efficiency"] - df["reference_efficiency"]
    return df


def plot_digestible(df: pd.DataFrame) -> None:
    fig, axes = plt.subplots(
        2,
        3,
        figsize=(15.8, 7.25),
        dpi=180,
        sharex=True,
        gridspec_kw={"height_ratios": [3.1, 1.0], "hspace": 0.12, "wspace": 0.16},
    )
    xlabels = list(ET_GROUPS.keys())
    x = np.arange(len(xlabels))
    for col, cent in enumerate(CENT_ORDER):
        g = df[df["centrality"] == cent].set_index("et_bin").loc[xlabels].reset_index()
        ax = axes[0, col]
        ax.errorbar(
            x,
            g["reference_truth_purity"],
            yerr=g["reference_truth_purity_error"],
            fmt="o",
            ms=7.0,
            lw=1.8,
            capsize=3.0,
            color="#6B7280",
            label="Reference box cuts",
            zorder=4,
        )
        ax.plot(
            x,
            g["bdt_truth_purity_at_reference_eff"],
            "s-",
            ms=7.2,
            lw=2.6,
            color="#009E73",
            label="Base v3E + centrality BDT\nthreshold matched to reference efficiency",
            zorder=5,
        )
        for i, row in g.iterrows():
            y0 = row["reference_truth_purity"]
            y1 = row["bdt_truth_purity_at_reference_eff"]
            ax.annotate(
                "",
                xy=(i, y1 - 0.0025),
                xytext=(i, y0 + 0.004),
                arrowprops=dict(arrowstyle="-|>", lw=1.5, color="#10B981", alpha=0.55),
            )
            ax.text(
                i + (0.06 if i == 0 else -0.02 if i == len(g) - 1 else 0.0),
                min(y1 + 0.006, 0.985),
                f"+{row['purity_gain_relative_percent']:.1f}%",
                ha="center",
                va="bottom",
                fontsize=9.4,
                color="#047857",
                fontweight="bold",
            )
        ax.set_xlim(-0.16, len(xlabels) - 0.84)
        ax.set_ylim(0.84, 0.995)
        ax.set_title(f"{CENT_LABEL[cent]} centrality", fontsize=13.5, fontweight="bold")
        if col == 0:
            ax.set_ylabel("Truth purity of selected sample")
            ax.legend(frameon=False, fontsize=8.8, loc="lower left")
        style_ax(ax, "y")

        ax2 = axes[1, col]
        ax2.errorbar(
            x,
            g["reference_efficiency"],
            yerr=g["reference_efficiency_error"],
            fmt="o",
            ms=5.6,
            lw=1.5,
            capsize=2.5,
            color="#6B7280",
        )
        ax2.plot(x, g["bdt_efficiency"], "s--", ms=5.8, lw=2.1, color="#009E73")
        ax2.set_ylim(0.11, 0.46)
        ax2.set_xticks(x, xlabels)
        ax2.set_xlabel(r"Cluster $E_T$ [GeV]")
        if col == 0:
            ax2.set_ylabel("Truth eff.")
        style_ax(ax2, "y")
    fig.text(0.055, 0.965, "BDT improves truth purity at the same reference-cut efficiency", fontsize=17.4, fontweight="bold")
    fig.text(
        0.055,
        0.918,
        "Base v3E + centrality BDT; threshold chosen in each centrality and $E_T$ bin to match the reference box-cut truth-photon efficiency.",
        fontsize=11.2,
        color="#4B5563",
    )
    sphinx_label(fig, x=0.855, y=0.957)
    fig.tight_layout(rect=(0.02, 0.02, 0.99, 0.89))
    fig.savefig(OUTDIR / "basev3e_cent_reference_matched_truth_purity_efficiency_2x3.png")
    plt.close(fig)


def plot_single_panel_summary(df: pd.DataFrame) -> None:
    # Compact slide-inset version: each centrality is a paired purity bar with
    # the efficiency closure as text below.
    rows = []
    for cent, g in df.groupby("centrality", sort=False):
        w = 1.0 / np.maximum(g["reference_truth_purity_error"].to_numpy(float), 1e-9) ** 2
        rows.append(
            {
                "centrality": cent,
                "centrality_label": CENT_LABEL[cent],
                "reference_purity": np.average(g["reference_truth_purity"], weights=w),
                "bdt_purity": np.average(g["bdt_truth_purity_at_reference_eff"], weights=w),
                "reference_eff": np.average(g["reference_efficiency"], weights=w),
                "bdt_eff": np.average(g["bdt_efficiency"], weights=w),
            }
        )
    s = pd.DataFrame(rows)
    s["gain_percent"] = 100.0 * (s["bdt_purity"] - s["reference_purity"]) / s["reference_purity"]
    s.to_csv(OUTDIR / "basev3e_cent_reference_matched_truth_purity_summary.csv", index=False)

    fig, ax = plt.subplots(figsize=(11.2, 5.9), dpi=180)
    x = np.arange(len(s))
    width = 0.33
    ax.bar(x - width / 2, s["reference_purity"], width, color="#9CA3AF", label="Reference box cuts")
    ax.bar(x + width / 2, s["bdt_purity"], width, color="#009E73", label="Base v3E + centrality BDT")
    for i, row in s.iterrows():
        ax.text(i + width / 2, row["bdt_purity"] + 0.006, f"+{row['gain_percent']:.1f}%", ha="center", va="bottom", color="#047857", fontsize=12.5, fontweight="bold")
        ax.text(i, 0.846, f"truth eff. fixed at {row['reference_eff']:.2f}", ha="center", va="bottom", fontsize=9.8, color="#374151")
    ax.set_xticks(x, s["centrality_label"])
    ax.set_ylim(0.84, 1.0)
    ax.set_ylabel("ET-averaged truth purity")
    ax.legend(frameon=False, fontsize=10.5, loc="upper left")
    style_ax(ax, "y")
    fig.text(0.065, 0.955, "Cleaner photons at fixed truth efficiency", fontsize=16.4, fontweight="bold")
    fig.text(
        0.065,
        0.905,
        "Base v3E + centrality BDT compared to reference box cuts; thresholds are retuned to match reference truth-photon efficiency.",
        fontsize=10.4,
        color="#4B5563",
    )
    sphinx_label(fig, x=0.94, y=0.953, size=13.8)
    fig.tight_layout(rect=(0.03, 0.02, 0.99, 0.875))
    fig.savefig(OUTDIR / "basev3e_cent_reference_matched_truth_purity_summary.png")
    plt.close(fig)


def build_target80_8bin_comparison() -> pd.DataFrame:
    ref_pur = load_reference_purity_for_bins(ET8_GROUPS)

    pur = pd.read_csv(BASEV3E_TARGET80_PUR_8)
    pur = pur[pur["source"] == "recoiljets_final_root"].copy()
    pur["centrality"] = pur["centrality"].astype(str)
    pur["et_bin"] = pur.apply(lambda r: f"{int(r['pt_low'])}-{int(r['pt_high'])}", axis=1)
    pur = pur.rename(
        columns={
            "pt_mid": "pt_center",
            "value": "bdt_truth_purity",
            "error": "bdt_truth_purity_error",
            "signal_tight": "bdt_signal_tight",
            "background_tight": "bdt_background_tight",
        }
    )[
        [
            "centrality",
            "et_bin",
            "pt_center",
            "bdt_truth_purity",
            "bdt_truth_purity_error",
            "bdt_signal_tight",
            "bdt_background_tight",
        ]
    ]

    eff = pd.read_csv(BASEV3E_TARGET80_EFF_8)
    ref_eff = eff[eff["source"] == "reference_box_cuts"].rename(
        columns={"value": "reference_efficiency", "error": "reference_efficiency_error"}
    )[["centrality", "pt_mid", "reference_efficiency", "reference_efficiency_error"]]
    bdt_eff = eff[eff["source"] == "basev3e_centrality_target80"].rename(
        columns={"value": "bdt_efficiency", "error": "bdt_efficiency_error"}
    )[["centrality", "pt_mid", "bdt_efficiency", "bdt_efficiency_error"]]
    for frame in (ref_eff, bdt_eff):
        frame["centrality"] = frame["centrality"].astype(str)
        frame["pt_center"] = frame["pt_mid"].astype(float)
        frame.drop(columns=["pt_mid"], inplace=True)

    df = ref_pur.merge(pur, on=["centrality", "et_bin", "pt_center"])
    df = df.merge(ref_eff, on=["centrality", "pt_center"]).merge(bdt_eff, on=["centrality", "pt_center"])
    df["centrality_label"] = df["centrality"].map(CENT_LABEL)
    df["purity_gain_absolute"] = df["bdt_truth_purity"] - df["reference_truth_purity"]
    df["purity_gain_relative_percent"] = 100.0 * df["purity_gain_absolute"] / df["reference_truth_purity"]
    df["efficiency_gain_absolute"] = df["bdt_efficiency"] - df["reference_efficiency"]
    return df


def plot_target80_8bin(df: pd.DataFrame) -> None:
    fig, axes = plt.subplots(
        2,
        3,
        figsize=(17.5, 7.8),
        dpi=180,
        sharex=True,
        gridspec_kw={"height_ratios": [2.25, 1.35], "hspace": 0.13, "wspace": 0.17},
    )
    labels = list(ET8_GROUPS.keys())
    x = np.arange(len(labels))
    for col, cent in enumerate(CENT_ORDER):
        g = df[df["centrality"] == cent].set_index("et_bin").loc[labels].reset_index()
        ax = axes[0, col]
        ax.errorbar(
            x,
            g["reference_truth_purity"],
            yerr=g["reference_truth_purity_error"],
            fmt="o-",
            ms=4.6,
            lw=1.6,
            capsize=2.4,
            color="#6B7280",
            label="Reference box cuts",
            zorder=3,
        )
        ax.errorbar(
            x,
            g["bdt_truth_purity"],
            yerr=g["bdt_truth_purity_error"],
            fmt="s-",
            ms=4.8,
            lw=2.0,
            capsize=2.4,
            color="#009E73",
            label="Base v3E + centrality BDT target-80",
            zorder=4,
        )
        for idx in [0, 2, 4, 7]:
            row = g.iloc[idx]
            ax.text(
                idx,
                max(row["reference_truth_purity"], row["bdt_truth_purity"]) + 0.010,
                f"{row['purity_gain_relative_percent']:+.1f}%",
                ha="center",
                va="bottom",
                fontsize=8.2,
                color="#047857" if row["purity_gain_relative_percent"] >= 0 else "#B91C1C",
                fontweight="bold",
            )
        ax.set_ylim(0.82, 0.94)
        ax.set_title(f"{CENT_LABEL[cent]} centrality", fontsize=12.8, fontweight="bold")
        if col == 0:
            ax.set_ylabel("Truth purity of selected sample")
            ax.legend(frameon=False, fontsize=8.6, loc="lower left")
        style_ax(ax, "y")

        ax2 = axes[1, col]
        ax2.errorbar(
            x,
            g["reference_efficiency"],
            yerr=g["reference_efficiency_error"],
            fmt="o-",
            ms=4.3,
            lw=1.5,
            capsize=2.2,
            color="#6B7280",
        )
        ax2.errorbar(
            x,
            g["bdt_efficiency"],
            yerr=g["bdt_efficiency_error"],
            fmt="s-",
            ms=4.5,
            lw=2.0,
            capsize=2.2,
            color="#009E73",
        )
        ax2.set_ylim(0.06, 1.02)
        ax2.set_xticks(x, labels, rotation=35, ha="right")
        ax2.set_xlabel(r"Cluster $E_T$ bin [GeV]")
        if col == 0:
            ax2.set_ylabel("Truth-photon recovery")
        style_ax(ax2, "y")

    fig.text(0.055, 0.965, "8-bin production view: BDT purity and efficiency versus reference box cuts", fontsize=17.0, fontweight="bold")
    fig.text(
        0.055,
        0.918,
        "Base v3E + centrality BDT target-80 applied in RecoilJets; this is the 8-bin production comparison, not the reference-efficiency-matched validation scan.",
        fontsize=10.6,
        color="#4B5563",
    )
    sphinx_label(fig, x=0.86, y=0.958, size=14.5)
    fig.tight_layout(rect=(0.02, 0.02, 0.99, 0.89))
    fig.savefig(OUTDIR / "basev3e_cent_reference_vs_target80_truth_purity_efficiency_8et_2x3.png")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Make clean truth-purity comparisons for base v3E + centrality BDT."
    )
    parser.add_argument(
        "--include-production-8bin",
        action="store_true",
        help=(
            "Also write the 8-bin RecoilJets target-80 production comparison. "
            "This is not the fixed-reference-efficiency validation scan."
        ),
    )
    args = parser.parse_args()

    OUTDIR.mkdir(parents=True, exist_ok=True)
    df = build_comparison()
    df.to_csv(OUTDIR / "basev3e_cent_reference_matched_truth_purity_efficiency.csv", index=False)
    plot_digestible(df)
    plot_single_panel_summary(df)
    print(OUTDIR / "basev3e_cent_reference_matched_truth_purity_efficiency_2x3.png")
    print(OUTDIR / "basev3e_cent_reference_matched_truth_purity_summary.png")
    if args.include_production_8bin:
        df8 = build_target80_8bin_comparison()
        df8.to_csv(OUTDIR / "basev3e_cent_reference_vs_target80_truth_purity_efficiency_8et.csv", index=False)
        plot_target80_8bin(df8)
        print(OUTDIR / "basev3e_cent_reference_vs_target80_truth_purity_efficiency_8et_2x3.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
