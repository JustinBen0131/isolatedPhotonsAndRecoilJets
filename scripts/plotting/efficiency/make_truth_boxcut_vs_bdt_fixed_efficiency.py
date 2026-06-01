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
REFERENCE_PUR = (
    ROOT
    / "slideReady"
    / "current_best_target80_recoiljets"
    / "reference_vs_ptcent7bdt_target80_reco_truth_raw_abcd_purity_pt15to35_points.csv"
)
PROD_EFF = (
    ROOT
    / "slideReady"
    / "current_best_target80_recoiljets"
    / "reference_vs_ptcent7bdt_target80_reco_efficiency_pt15to35_points.csv"
)
PROD_PUR = REFERENCE_PUR

CENT_ORDER = ["0_20", "20_50", "50_80"]
CENT_LABEL = {"0_20": "0-20%", "20_50": "20-50%", "50_80": "50-80%"}
ET_BINS = [
    ("15-17", 16.0),
    ("17-19", 18.0),
    ("19-21", 20.0),
    ("21-23", 22.0),
    ("23-25", 24.0),
    ("25-27", 26.0),
    ("27-30", 28.5),
    ("30-35", 32.5),
]


def style_ax(ax, grid_axis: str = "both") -> None:
    ax.tick_params(direction="in", top=True, right=True, labelsize=10.5)
    ax.grid(True, axis=grid_axis, color="#D1D5DB", lw=0.7, alpha=0.75)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_linewidth(1.0)


def sphinx_label(fig, x=0.855, y=0.955, size=15) -> None:
    fig.text(x, y, "sPHENIX", ha="right", fontsize=size, fontstyle="italic", fontweight="bold")
    fig.text(x + 0.095, y, "Internal", ha="right", fontsize=size)


def coarse_cent(label: str) -> str:
    low = int(str(label).split("-")[0])
    if low < 20:
        return "0_20"
    if low < 50:
        return "20_50"
    return "50_80"


def bdt_curve_for_bin(hist: pd.DataFrame) -> pd.DataFrame:
    g = hist.sort_values("score_low").copy()
    sig = g["signal_count"].to_numpy(float)
    bkg = g["background_count"].to_numpy(float)
    sig_total = sig.sum()
    bkg_total = bkg.sum()
    sig_above = sig[::-1].cumsum()[::-1]
    bkg_above = bkg[::-1].cumsum()[::-1]
    total_above = sig_above + bkg_above
    out = g[["score_low", "score_high"]].copy()
    out["score_center"] = 0.5 * (g["score_low"].to_numpy() + g["score_high"].to_numpy())
    out["signal_efficiency"] = sig_above / sig_total if sig_total > 0 else np.nan
    out["truth_purity"] = sig_above / np.where(total_above > 0, total_above, np.nan)
    out["background_fake_rate"] = bkg_above / bkg_total if bkg_total > 0 else np.nan
    out["signal_above"] = sig_above
    out["background_above"] = bkg_above
    out["signal_total"] = sig_total
    out["background_total"] = bkg_total
    return out


def interpolate_at_efficiency(curve: pd.DataFrame, target_eff: float) -> dict[str, float]:
    curve = curve.replace([np.inf, -np.inf], np.nan).dropna(subset=["signal_efficiency", "truth_purity"])
    curve = curve.sort_values("signal_efficiency").copy()
    x = curve["signal_efficiency"].to_numpy(float)
    # Drop duplicate x values so np.interp sees a well-defined monotonic coordinate.
    keep = np.r_[True, np.diff(x) > 1e-12]
    curve = curve.iloc[keep].copy()
    x = curve["signal_efficiency"].to_numpy(float)
    target = float(np.clip(target_eff, np.nanmin(x), np.nanmax(x)))
    return {
        "signal_efficiency": target,
        "truth_purity": float(np.interp(target, x, curve["truth_purity"].to_numpy(float))),
        "background_fake_rate": float(np.interp(target, x, curve["background_fake_rate"].to_numpy(float))),
        "score_threshold": float(np.interp(target, x, curve["score_low"].to_numpy(float))),
        "signal_total": float(curve["signal_total"].iloc[0]),
        "background_total": float(curve["background_total"].iloc[0]),
    }


def load_bdt_fixed_efficiency_points() -> pd.DataFrame:
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

    ref_eff = pd.read_csv(REFERENCE_EFF)
    ref_eff = ref_eff[ref_eff["label"].str.contains("Box", case=False, na=False)].copy()
    ref_eff["centrality"] = ref_eff["cent"].astype(str)
    ref_eff["pt_mid"] = ref_eff["pt_mid"].astype(float)
    ref_eff = ref_eff[ref_eff["centrality"].isin(CENT_ORDER)].copy()

    rows = []
    for cent in CENT_ORDER:
        for et_label, pt_mid in ET_BINS:
            eff_rows = ref_eff[(ref_eff["centrality"] == cent) & (np.isclose(ref_eff["pt_mid"], pt_mid))]
            h = hist[(hist["centrality"] == cent) & (hist["et_bin"] == et_label)]
            if len(eff_rows) == 0 or len(h) == 0:
                continue
            target_eff = float(eff_rows.iloc[0]["value"])
            curve = bdt_curve_for_bin(h)
            r = interpolate_at_efficiency(curve, target_eff)
            rows.append(
                {
                    "model": "global32_bdt_matched_to_reference_eff",
                    "label": "BDT threshold at reference efficiency",
                    "centrality": cent,
                    "centrality_label": CENT_LABEL[cent],
                    "et_bin": et_label,
                    "pt_mid": pt_mid,
                    "reference_efficiency": target_eff,
                    "bdt_actual_efficiency": float(r["signal_efficiency"]),
                    "bdt_truth_purity": float(r["truth_purity"]),
                    "bdt_background_fake_rate": float(r["background_fake_rate"]),
                    "bdt_score_threshold": float(r["score_threshold"]),
                    "signal_total": float(r["signal_total"]),
                    "background_total": float(r["background_total"]),
                    "definition": "BDT score threshold chosen in this centrality x ET bin to match reference box-cut truth-photon recovery",
                }
            )
    return pd.DataFrame(rows)


def load_reference_purity_points() -> pd.DataFrame:
    ref = pd.read_csv(REFERENCE_PUR)
    ref = ref[(ref["model"] == "Reference box cuts") & (ref["quantity"] == "truth_labeled_A_purity")].copy()
    ref["centrality"] = ref["centrality"].astype(str)
    ref["pt_mid"] = ref["pt_mid"].astype(float)
    # The existing reference-purity production uses 14-16, 16-18, ..., 26-35 bins.
    # Map those centers onto the reference-efficiency/validation ET centers.
    # The final wide 26-35 bin is compared to the 27-30 validation center and
    # the 30-35 point is left out of the fixed-efficiency overlay.
    center_map = {15.0: 16.0, 17.0: 18.0, 19.0: 20.0, 21.0: 22.0, 23.0: 24.0, 25.0: 26.0, 30.5: 28.5}
    ref["matched_pt_mid"] = ref["pt_mid"].map(center_map)
    ref = ref.dropna(subset=["matched_pt_mid"]).copy()
    return ref


def build_fixed_efficiency_comparison() -> pd.DataFrame:
    bdt = load_bdt_fixed_efficiency_points()
    ref = load_reference_purity_points()
    eff = pd.read_csv(REFERENCE_EFF)
    eff = eff[eff["label"].str.contains("Box", case=False, na=False)].copy()
    eff["centrality"] = eff["cent"].astype(str)
    eff["pt_mid"] = eff["pt_mid"].astype(float)
    eff = eff[["centrality", "pt_mid", "value", "error"]].rename(
        columns={"value": "reference_efficiency", "error": "reference_efficiency_error"}
    )

    rows = []
    for _, b in bdt.iterrows():
        rref = ref[(ref["centrality"] == b["centrality"]) & (np.isclose(ref["matched_pt_mid"], b["pt_mid"]))]
        reff = eff[(eff["centrality"] == b["centrality"]) & (np.isclose(eff["pt_mid"], b["pt_mid"]))]
        if len(rref) == 0 or len(reff) == 0:
            continue
        rr = rref.iloc[0]
        ee = reff.iloc[0]
        rows.append(
            {
                "centrality": b["centrality"],
                "centrality_label": b["centrality_label"],
                "et_bin": b["et_bin"],
                "pt_mid": b["pt_mid"],
                "reference_truth_purity": float(rr["value"]),
                "reference_truth_purity_error": float(rr["error"]),
                "reference_efficiency": float(ee["reference_efficiency"]),
                "reference_efficiency_error": float(ee["reference_efficiency_error"]),
                "bdt_truth_purity_at_same_efficiency": float(b["bdt_truth_purity"]),
                "bdt_actual_efficiency": float(b["bdt_actual_efficiency"]),
                "bdt_score_threshold": float(b["bdt_score_threshold"]),
                "purity_gain_absolute": float(b["bdt_truth_purity"] - rr["value"]),
                "purity_gain_relative_percent": 100.0 * float((b["bdt_truth_purity"] - rr["value"]) / rr["value"]),
                "note": "BDT is rethresholded to match the reference box-cut truth efficiency in the same coarse-centrality and ET bin.",
            }
        )
    return pd.DataFrame(rows)


def weighted_integrated(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for cent, g in df.groupby("centrality", sort=False):
        w = 1.0 / np.maximum(g["reference_truth_purity_error"].to_numpy(float), 1e-6) ** 2
        ref_p = np.average(g["reference_truth_purity"], weights=w)
        bdt_p = np.average(g["bdt_truth_purity_at_same_efficiency"], weights=w)
        ref_e = np.average(g["reference_efficiency"], weights=w)
        bdt_e = np.average(g["bdt_actual_efficiency"], weights=w)
        rows.append(
            {
                "centrality": cent,
                "centrality_label": CENT_LABEL[cent],
                "reference_truth_purity": ref_p,
                "bdt_truth_purity_at_same_efficiency": bdt_p,
                "reference_efficiency": ref_e,
                "bdt_actual_efficiency": bdt_e,
                "purity_gain_absolute": bdt_p - ref_p,
                "purity_gain_relative_percent": 100.0 * (bdt_p - ref_p) / ref_p,
            }
        )
    return pd.DataFrame(rows)


def plot_fixed_efficiency(df: pd.DataFrame) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(16.0, 8.0), dpi=180, sharex=True)
    for col, cent in enumerate(CENT_ORDER):
        g = df[df["centrality"] == cent].sort_values("pt_mid")
        ax = axes[0, col]
        ax.errorbar(
            g["pt_mid"],
            g["reference_truth_purity"],
            yerr=g["reference_truth_purity_error"],
            fmt="o",
            ms=5.5,
            lw=1.7,
            capsize=2.8,
            color="#6B7280",
            label="Reference box cuts",
        )
        ax.plot(g["pt_mid"], g["bdt_truth_purity_at_same_efficiency"], "s-", ms=5.0, lw=2.2, color="#009E73", label="BDT, same truth efficiency")
        for x, y0, y1 in zip(g["pt_mid"], g["reference_truth_purity"], g["bdt_truth_purity_at_same_efficiency"]):
            ax.plot([x, x], [y0, y1], color="#009E73", lw=1.1, alpha=0.35)
        ax.set_ylim(0.84, 0.985)
        ax.set_title(f"{CENT_LABEL[cent]} centrality", fontsize=13.5, fontweight="bold")
        if col == 0:
            ax.set_ylabel("Truth purity of selected sample")
            ax.legend(frameon=False, fontsize=9.2, loc="lower left")
        style_ax(ax, "y")

        ax2 = axes[1, col]
        ax2.errorbar(
            g["pt_mid"],
            g["reference_efficiency"],
            yerr=g["reference_efficiency_error"],
            fmt="o",
            ms=5.5,
            lw=1.7,
            capsize=2.8,
            color="#6B7280",
            label="Reference box cuts",
        )
        ax2.plot(g["pt_mid"], g["bdt_actual_efficiency"], "s-", ms=5.0, lw=2.2, color="#009E73", label="BDT threshold match")
        ax2.set_ylim(0.08, 0.48)
        ax2.set_xlabel(r"Cluster $E_T$ bin center [GeV]")
        if col == 0:
            ax2.set_ylabel("Truth signal efficiency kept")
        style_ax(ax2, "y")
    fig.text(0.055, 0.965, "Truth purity gain at fixed signal efficiency", fontsize=18.0, fontweight="bold")
    fig.text(
        0.055,
        0.918,
        "For each centrality and photon-$E_T$ bin, the BDT threshold is chosen to match the reference box-cut truth efficiency.",
        fontsize=11.7,
        color="#4B5563",
    )
    sphinx_label(fig, x=0.855, y=0.957)
    fig.tight_layout(rect=(0.02, 0.02, 0.99, 0.89))
    fig.savefig(OUTDIR / "variant_m_reference_boxcuts_vs_bdt_fixed_efficiency_truth_purity_2x3.png")
    plt.close(fig)


def plot_integrated_summary(df: pd.DataFrame) -> None:
    integ = weighted_integrated(df)
    integ.to_csv(OUTDIR / "variant_n_reference_boxcuts_vs_bdt_fixed_efficiency_integrated_summary.csv", index=False)
    fig, axes = plt.subplots(1, 2, figsize=(13.8, 5.8), dpi=180)
    x = np.arange(len(integ))
    width = 0.34
    axes[0].bar(x - width / 2, integ["reference_truth_purity"], width, color="#9CA3AF", label="Reference box cuts")
    axes[0].bar(x + width / 2, integ["bdt_truth_purity_at_same_efficiency"], width, color="#009E73", label="BDT, same efficiency")
    for i, row in integ.iterrows():
        axes[0].text(i + width / 2, row["bdt_truth_purity_at_same_efficiency"] + 0.004, f"+{row['purity_gain_relative_percent']:.1f}%", ha="center", va="bottom", fontsize=10.5, color="#047857", fontweight="bold")
    axes[0].set_xticks(x, integ["centrality_label"])
    axes[0].set_ylim(0.84, 1.005)
    axes[0].set_ylabel("ET-averaged truth purity")
    axes[0].legend(frameon=False, fontsize=9.5, loc="lower right")
    style_ax(axes[0], "y")

    axes[1].plot(x, integ["reference_efficiency"], "o-", color="#6B7280", lw=2.2, label="Reference")
    axes[1].plot(x, integ["bdt_actual_efficiency"], "s--", color="#009E73", lw=2.2, label="BDT threshold match")
    axes[1].set_xticks(x, integ["centrality_label"])
    axes[1].set_ylim(0.22, 0.40)
    axes[1].set_ylabel("ET-averaged signal efficiency")
    axes[1].legend(frameon=False, fontsize=9.5, loc="lower right")
    style_ax(axes[1], "y")

    fig.text(0.055, 0.955, "Integrated fixed-efficiency truth-purity check", fontsize=15.8, fontweight="bold")
    fig.text(0.055, 0.908, "BDT threshold is tuned to the reference efficiency; any vertical gain is purity, not signal loss.", fontsize=10.8, color="#4B5563")
    sphinx_label(fig, x=0.86, y=0.948, size=14.5)
    fig.tight_layout(rect=(0.03, 0.02, 0.99, 0.875))
    fig.savefig(OUTDIR / "variant_n_reference_boxcuts_vs_bdt_fixed_efficiency_integrated_summary.png")
    plt.close(fig)


def plot_production_selected_sample() -> None:
    pur = pd.read_csv(PROD_PUR)
    pur = pur[pur["quantity"] == "truth_labeled_A_purity"].copy()
    eff = pd.read_csv(PROD_EFF)
    # This is the production-style selected-region comparison. It is not
    # efficiency-matched, so keep it separate from the fixed-efficiency proof.
    fig, axes = plt.subplots(1, 3, figsize=(15.8, 4.9), dpi=180, sharey=True)
    colors = {"Reference box cuts": "#6B7280", "8x7 routed BDT target-80": "#009E73"}
    markers = {"Reference box cuts": "o", "8x7 routed BDT target-80": "s"}
    for ax, cent in zip(axes, CENT_ORDER):
        for model in ["Reference box cuts", "8x7 routed BDT target-80"]:
            g = pur[(pur["model"] == model) & (pur["centrality"] == cent)].sort_values("pt_mid")
            ax.errorbar(
                g["pt_mid"],
                g["value"],
                yerr=g["error"],
                fmt=markers[model] + "-",
                color=colors[model],
                lw=2.0,
                ms=5.2,
                capsize=2.4,
                label=model if cent == "0_20" else None,
            )
        ax.set_title(f"{CENT_LABEL[cent]} centrality", fontsize=13.5, fontweight="bold")
        ax.set_xlabel(r"Cluster $E_T$ bin center [GeV]")
        ax.set_ylim(0.84, 0.94)
        if ax is axes[0]:
            ax.set_ylabel("Truth-labeled purity in selected region A")
            ax.legend(frameon=False, fontsize=9.2, loc="lower left")
        style_ax(ax, "y")
    fig.text(0.055, 0.955, "Production selected-region truth purity: reference cuts vs BDT", fontsize=17.2, fontweight="bold")
    fig.text(0.055, 0.907, "This uses the actual RecoilJets selected A region; use the fixed-efficiency plot for the apples-to-apples proof.", fontsize=11.3, color="#4B5563")
    sphinx_label(fig, x=0.855, y=0.948)
    fig.tight_layout(rect=(0.02, 0.02, 0.99, 0.88))
    fig.savefig(OUTDIR / "variant_o_reference_boxcuts_vs_bdt_production_truth_purity_1x3.png")
    plt.close(fig)


def main() -> int:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    df = build_fixed_efficiency_comparison()
    df.to_csv(OUTDIR / "variant_m_reference_boxcuts_vs_bdt_fixed_efficiency_truth_purity.csv", index=False)
    plot_fixed_efficiency(df)
    plot_integrated_summary(df)
    plot_production_selected_sample()
    print(OUTDIR)
    for name in [
        "variant_m_reference_boxcuts_vs_bdt_fixed_efficiency_truth_purity_2x3.png",
        "variant_n_reference_boxcuts_vs_bdt_fixed_efficiency_integrated_summary.png",
        "variant_o_reference_boxcuts_vs_bdt_production_truth_purity_1x3.png",
    ]:
        print(OUTDIR / name)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
