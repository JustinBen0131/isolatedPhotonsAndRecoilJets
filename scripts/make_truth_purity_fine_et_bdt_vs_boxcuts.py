#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


ROOT = Path(
    "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/auauMLDiagnosticRuns/"
    "global_etcent_inclusive3_sixpack_20260516_135439"
)
OUTDIR = ROOT / "slideReady" / "truth_bdt_performance_variants"

GLOBAL_HIST = ROOT / "validation" / "bdt_finished_only_20260516_180817" / "fine7x8_score_histograms_globalEtCent1535_bdt_iso_noIso.csv"
REFERENCE_EFF = Path(
    "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/target80_first_offline/"
    "bdt_target80_gated_20260512_001012/analysis_config_etfine_15to35_target80/"
    "efficiency_qa/truth_photon_recovery_box_cuts_points.csv"
)
REFERENCE_COUNTS = (
    ROOT
    / "slideReady"
    / "abcd_raw_purity_reference_vs_basev3e"
    / "reference_boxcuts_vs_basev3e_centrality_target80_reco_abcd_raw_purity_pt15to35_points.csv"
)

CSV_OUT = OUTDIR / "basev3e_cent_bdt_vs_boxcuts_truth_purity_fine_et_effmatched.csv"
LINE_PNG = OUTDIR / "basev3e_cent_bdt_vs_boxcuts_truth_purity_fine_et_effmatched_1x3.png"
GAIN_PNG = OUTDIR / "basev3e_cent_bdt_vs_boxcuts_truth_purity_gain_fine_et_heatmap.png"

CENT_ORDER = ["0_20", "20_50", "50_80"]
CENT_LABEL = {"0_20": "0-20%", "20_50": "20-50%", "50_80": "50-80%"}
ET_BINS = [
    ("15-17", 15.0, 17.0, 16.0),
    ("17-19", 17.0, 19.0, 18.0),
    ("19-21", 19.0, 21.0, 20.0),
    ("21-23", 21.0, 23.0, 22.0),
    ("23-25", 23.0, 25.0, 24.0),
    ("25-27", 25.0, 27.0, 26.0),
    ("27-30", 27.0, 30.0, 28.5),
    ("30-35", 30.0, 35.0, 32.5),
]


def sphinx_label(fig, x=0.86, y=0.955, size=14.0) -> None:
    fig.text(x, y, "sPHENIX", ha="right", fontsize=size, fontstyle="italic", fontweight="bold")
    fig.text(x + 0.088, y, "Internal", ha="right", fontsize=size)


def style_axis(ax, grid_axis="y") -> None:
    ax.tick_params(direction="in", top=True, right=True, labelsize=10.5)
    if grid_axis:
        ax.grid(True, axis=grid_axis, color="#D1D5DB", lw=0.75, alpha=0.82)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_linewidth(1.05)


def coarse_cent(fine_label: str) -> str:
    low = int(str(fine_label).split("-")[0])
    if low < 20:
        return "0_20"
    if low < 50:
        return "20_50"
    return "50_80"


def overlap_fraction(source_low: float, source_high: float, target_low: float, target_high: float) -> float:
    overlap = max(0.0, min(source_high, target_high) - max(source_low, target_low))
    return overlap / max(source_high - source_low, 1e-12)


def ratio_error(signal: float, background: float) -> tuple[float, float]:
    total = signal + background
    if total <= 0.0:
        return np.nan, np.nan
    value = signal / total
    return float(value), float(np.sqrt(max(value * (1.0 - value) / total, 0.0)))


def load_reference_purity() -> pd.DataFrame:
    ref = pd.read_csv(REFERENCE_COUNTS)
    ref = ref[(ref["model"] == "Reference box cuts") & (ref["quantity"] == "truth_labeled_A_purity")].copy()
    rows = []
    for cent in CENT_ORDER:
        gcent = ref[ref["centrality"].astype(str) == cent]
        for et_bin, et_low, et_high, et_mid in ET_BINS:
            signal = 0.0
            background = 0.0
            source_bins = []
            for _, row in gcent.iterrows():
                frac = overlap_fraction(float(row["pt_low"]), float(row["pt_high"]), et_low, et_high)
                if frac <= 0.0:
                    continue
                signal += float(row["signal_tight"]) * frac
                background += float(row["background_tight"]) * frac
                source_bins.append(f"{row['pt_low']:g}-{row['pt_high']:g}:{frac:.3g}")
            purity, err = ratio_error(signal, background)
            rows.append(
                {
                    "centrality": cent,
                    "centrality_label": CENT_LABEL[cent],
                    "et_bin": et_bin,
                    "et_low": et_low,
                    "et_high": et_high,
                    "pt_mid": et_mid,
                    "reference_truth_purity": purity,
                    "reference_truth_purity_error": err,
                    "reference_purity_source_bins": ";".join(source_bins),
                }
            )
    return pd.DataFrame(rows)


def load_reference_efficiency() -> pd.DataFrame:
    eff = pd.read_csv(REFERENCE_EFF)
    eff = eff[eff["label"].str.contains("Box", case=False, na=False)].copy()
    eff["centrality"] = eff["cent"].astype(str)
    eff = eff[eff["centrality"].isin(CENT_ORDER)].copy()
    rows = []
    for _, row in eff.iterrows():
        pt_mid = float(row["pt_mid"])
        match = [b for b in ET_BINS if abs(b[3] - pt_mid) < 1e-9]
        if not match:
            continue
        et_bin, et_low, et_high, _ = match[0]
        rows.append(
            {
                "centrality": str(row["centrality"]),
                "et_bin": et_bin,
                "reference_efficiency": float(row["value"]),
                "reference_efficiency_error": float(row["error"]),
            }
        )
    return pd.DataFrame(rows)


def bdt_curve_for_bin(hist: pd.DataFrame) -> pd.DataFrame:
    g = hist.sort_values("score_low").copy()
    sig = g["signal_count"].to_numpy(float)
    bkg = g["background_count"].to_numpy(float)
    sig_total = float(sig.sum())
    bkg_total = float(bkg.sum())
    sig_above = sig[::-1].cumsum()[::-1]
    bkg_above = bkg[::-1].cumsum()[::-1]
    total_above = sig_above + bkg_above
    out = g[["score_low", "score_high"]].copy()
    out["signal_efficiency"] = sig_above / sig_total if sig_total > 0 else np.nan
    out["truth_purity"] = sig_above / np.where(total_above > 0, total_above, np.nan)
    out["background_fake_rate"] = bkg_above / bkg_total if bkg_total > 0 else np.nan
    out["signal_above"] = sig_above
    out["background_above"] = bkg_above
    out["signal_total"] = sig_total
    out["background_total"] = bkg_total
    return out


def interp_curve(curve: pd.DataFrame, target_eff: float) -> dict[str, float]:
    c = curve.replace([np.inf, -np.inf], np.nan).dropna(subset=["signal_efficiency", "truth_purity"])
    c = c.sort_values("signal_efficiency").copy()
    x = c["signal_efficiency"].to_numpy(float)
    keep = np.r_[True, np.diff(x) > 1e-12]
    c = c.iloc[keep].copy()
    x = c["signal_efficiency"].to_numpy(float)
    target = float(np.clip(target_eff, np.nanmin(x), np.nanmax(x)))
    sig_above = float(np.interp(target, x, c["signal_above"].to_numpy(float)))
    bkg_above = float(np.interp(target, x, c["background_above"].to_numpy(float)))
    purity, purity_err = ratio_error(sig_above, bkg_above)
    return {
        "actual_signal_efficiency": target,
        "truth_purity": purity,
        "truth_purity_error": purity_err,
        "score_threshold": float(np.interp(target, x, c["score_low"].to_numpy(float))),
        "background_fake_rate": float(np.interp(target, x, c["background_fake_rate"].to_numpy(float))),
        "signal_above": sig_above,
        "background_above": bkg_above,
    }


def load_bdt_matched() -> pd.DataFrame:
    hist = pd.read_csv(GLOBAL_HIST)
    hist = hist[hist["product"] == "globalEtCent1535_bdt_noIso"].copy()
    hist["centrality"] = hist["centrality_bin"].map(coarse_cent)
    hist = (
        hist.groupby(["centrality", "et_bin", "score_low", "score_high"], as_index=False)[
            ["signal_count", "background_count"]
        ]
        .sum()
        .sort_values(["centrality", "et_bin", "score_low"])
    )
    eff = load_reference_efficiency()
    rows = []
    for cent in CENT_ORDER:
        for et_bin, _, _, et_mid in ET_BINS:
            e = eff[(eff["centrality"] == cent) & (eff["et_bin"] == et_bin)]
            h = hist[(hist["centrality"] == cent) & (hist["et_bin"] == et_bin)]
            if len(e) == 0 or len(h) == 0:
                raise SystemExit(f"Missing BDT hist/reference efficiency for {cent} {et_bin}")
            target_eff = float(e.iloc[0]["reference_efficiency"])
            r = interp_curve(bdt_curve_for_bin(h), target_eff)
            rows.append(
                {
                    "centrality": cent,
                    "et_bin": et_bin,
                    "pt_mid": et_mid,
                    "bdt_truth_purity_at_reference_eff": r["truth_purity"],
                    "bdt_truth_purity_error": r["truth_purity_error"],
                    "bdt_actual_efficiency": r["actual_signal_efficiency"],
                    "bdt_score_threshold": r["score_threshold"],
                    "bdt_background_fake_rate": r["background_fake_rate"],
                    "bdt_signal_above_interp": r["signal_above"],
                    "bdt_background_above_interp": r["background_above"],
                }
            )
    return pd.DataFrame(rows)


def build_table() -> pd.DataFrame:
    table = load_reference_purity().merge(load_reference_efficiency(), on=["centrality", "et_bin"], how="inner")
    table = table.merge(load_bdt_matched(), on=["centrality", "et_bin", "pt_mid"], how="inner")
    table["purity_gain_absolute"] = table["bdt_truth_purity_at_reference_eff"] - table["reference_truth_purity"]
    table["purity_gain_pp"] = 100.0 * table["purity_gain_absolute"]
    table["purity_gain_relative_percent"] = 100.0 * table["purity_gain_absolute"] / table["reference_truth_purity"]
    expected = len(CENT_ORDER) * len(ET_BINS)
    if len(table) != expected:
        raise SystemExit(f"Expected {expected} rows, found {len(table)}")
    return table.sort_values(["centrality", "et_low"]).reset_index(drop=True)


def plot_line(table: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(16.4, 5.35), dpi=190, sharey=True, gridspec_kw={"wspace": 0.10})
    for col, cent in enumerate(CENT_ORDER):
        ax = axes[col]
        g = table[table["centrality"] == cent].sort_values("et_low").reset_index(drop=True)
        x = g["pt_mid"].to_numpy(float)
        xerr = 0.5 * (g["et_high"].to_numpy(float) - g["et_low"].to_numpy(float))
        ax.errorbar(
            x,
            g["reference_truth_purity"],
            xerr=xerr,
            yerr=g["reference_truth_purity_error"],
            fmt="o",
            ms=5.2,
            lw=0.0,
            elinewidth=1.4,
            capsize=2.5,
            color="#6B7280",
            label="Reference box cuts" if col == 0 else None,
            zorder=3,
        )
        ax.errorbar(
            x,
            g["bdt_truth_purity_at_reference_eff"],
            xerr=xerr,
            yerr=g["bdt_truth_purity_error"],
            fmt="s",
            ms=5.4,
            lw=0.0,
            elinewidth=1.7,
            capsize=2.5,
            color="#009E73",
            label="BDT at fixed reference efficiency" if col == 0 else None,
            zorder=4,
        )
        ax.set_title(f"{CENT_LABEL[cent]} centrality", fontsize=13.5, fontweight="bold", pad=10)
        ax.set_xticks([15, 20, 25, 30, 35])
        ax.set_xlabel(r"Cluster $E_T$ [GeV]", fontsize=11.6)
        ax.set_ylim(0.835, 0.986)
        ax.set_xlim(14.2, 35.8)
        if col == 0:
            ax.set_ylabel("Truth purity of selected photons", fontsize=12.6)
            ax.legend(frameon=False, fontsize=10.4, loc="lower left")
        style_axis(ax, grid_axis=None)

    fig.text(
        0.055,
        0.956,
        "Truth purity at fixed reference-cut efficiency",
        fontsize=17.0,
        fontweight="bold",
    )
    sphinx_label(fig, x=0.855, y=0.955)
    fig.subplots_adjust(left=0.065, right=0.985, bottom=0.165, top=0.865, wspace=0.10)
    fig.savefig(LINE_PNG)
    plt.close(fig)


def plot_gain_heatmap(table: pd.DataFrame) -> None:
    labels = [b[0] for b in ET_BINS]
    mat = np.zeros((len(CENT_ORDER), len(labels)))
    for iy, cent in enumerate(CENT_ORDER):
        g = table[table["centrality"] == cent].set_index("et_bin").loc[labels]
        mat[iy, :] = g["purity_gain_pp"].to_numpy(float)

    cmap = LinearSegmentedColormap.from_list("gain", ["#ECFDF5", "#A7F3D0", "#10B981", "#047857"])
    fig, ax = plt.subplots(figsize=(12.8, 4.7), dpi=190)
    im = ax.imshow(mat, aspect="auto", cmap=cmap, vmin=0.0, vmax=max(12.5, float(np.nanmax(mat))))
    ax.set_xticks(np.arange(len(labels)), labels, rotation=35, ha="right")
    ax.set_yticks(np.arange(len(CENT_ORDER)), [CENT_LABEL[c] for c in CENT_ORDER])
    ax.set_xlabel(r"Cluster $E_T$ bin [GeV]", fontsize=12.2)
    ax.set_ylabel("Centrality", fontsize=12.2)
    ax.tick_params(labelsize=11.0)
    for iy in range(mat.shape[0]):
        for ix in range(mat.shape[1]):
            val = mat[iy, ix]
            ax.text(ix, iy, f"+{val:.1f}", ha="center", va="center", fontsize=11.2, fontweight="bold", color="#052E16")
    for spine in ax.spines.values():
        spine.set_linewidth(1.0)
    cbar = fig.colorbar(im, ax=ax, fraction=0.032, pad=0.018)
    cbar.set_label("Truth-purity gain [percentage points]", fontsize=11.3)
    cbar.ax.tick_params(labelsize=10.0)
    fig.text(
        0.055,
        0.952,
        "BDT truth-purity gain at matched truth recovery",
        fontsize=16.6,
        fontweight="bold",
    )
    fig.text(
        0.055,
        0.902,
        "Cell value = BDT purity minus box-cut purity, using the same truth-photon recovery in that bin.",
        fontsize=11.2,
        color="#4B5563",
    )
    sphinx_label(fig, x=0.835, y=0.952, size=13.5)
    fig.subplots_adjust(left=0.105, right=0.93, bottom=0.19, top=0.82)
    fig.savefig(GAIN_PNG)
    plt.close(fig)


def main() -> int:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    table = build_table()
    table.to_csv(CSV_OUT, index=False)
    plot_line(table)
    plot_gain_heatmap(table)
    print(CSV_OUT)
    print(LINE_PNG)
    print(GAIN_PNG)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
