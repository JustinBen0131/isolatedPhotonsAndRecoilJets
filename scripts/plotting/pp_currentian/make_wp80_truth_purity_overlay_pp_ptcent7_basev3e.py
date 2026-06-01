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
import io

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
ROOT = REPO / "dataOutput/auauMLDiagnosticRuns/global_etcent_inclusive3_sixpack_20260516_135439"
OUTDIR = ROOT / "slideReady/truth_bdt_performance_variants"

PP_HIST = (
    REPO
    / "dataOutput/ppPhotonMLPipeline/ppg12_fix4_score_overlays/score_separation/"
    / "pp_baseline_bdt_score_histograms_compact.csv"
)
PTCENT7_WP80_CELLS = ROOT / "slideReady/ml_bdt_wp80_cent3_fits/binned_bdt_noiso_ptcent7_wp80_cent3_threshold_cells.csv"
BASEV3E_HIST = ROOT / "validation/basev3e_controls_20260518_1110/basev3e_with_centrality_score_histograms.csv"

OUT_CSV = OUTDIR / "wp80_truth_purity_overlay_pp_ptcent7_basev3e_points.csv"
OUT_PNG = OUTDIR / "wp80_truth_purity_overlay_pp_ptcent7_basev3e_fine_et_1x3.png"

CENT_ORDER = ["0_20", "20_50", "50_80"]
CENT_LABEL = {"0_20": "0-20%", "20_50": "20-50%", "50_80": "50-80%"}
DISPLAY_ET_BINS = [
    (15.0, 17.0, "15-17"),
    (17.0, 19.0, "17-19"),
    (19.0, 21.0, "19-21"),
    (21.0, 23.0, "21-23"),
    (23.0, 25.0, "23-25"),
    (25.0, 30.0, "25-30"),
    (30.0, 35.0, "30-35"),
]

MODEL_STYLE = {
    "pp_baseline": {
        "label": "pp baseline BDT",
        "color": "#6B7280",
        "marker": "D",
        "linestyle": "None",
        "zorder": 5,
    },
    "ptcent7_best": {
        "label": "8 ET x 7 cent 32-feature BDT",
        "color": "#009E73",
        "marker": "o",
        "linestyle": "None",
        "zorder": 4,
    },
    "basev3e_cent": {
        "label": "base v3E + centrality BDT",
        "color": "#0072B2",
        "marker": "s",
        "linestyle": "None",
        "zorder": 3,
    },
}


def coarse_cent(label: str) -> str:
    low = int(str(label).replace("%", "").replace("_", "-").split("-")[0])
    if low < 20:
        return "0_20"
    if low < 50:
        return "20_50"
    return "50_80"


def et_mid(label: str) -> float:
    lo, hi = str(label).split("-")
    return 0.5 * (float(lo) + float(hi))


def interval_overlap(a0: float, a1: float, b0: float, b1: float) -> float:
    return max(0.0, min(a1, b1) - max(a0, b0))


def read_csv_after_header(path: Path, header_prefix: str) -> pd.DataFrame:
    text = path.read_text()
    lines = text.splitlines()
    idx = next(i for i, line in enumerate(lines) if line.startswith(header_prefix))
    return pd.read_csv(io.StringIO("\n".join(lines[idx:])))


def aggregate_native_selected_to_display(selected: pd.DataFrame, model: str, native_binning: str, source: Path) -> pd.DataFrame:
    rows = []
    for cent in CENT_ORDER:
        cent_rows = selected[selected["centrality"] == cent]
        for et_low, et_high, et_bin in DISPLAY_ET_BINS:
            sig = bkg = sig_entries = bkg_entries = 0.0
            thresholds = []
            for _, r in cent_rows.iterrows():
                width = float(r["et_high"] - r["et_low"])
                if width <= 0:
                    continue
                frac = interval_overlap(float(r["et_low"]), float(r["et_high"]), et_low, et_high) / width
                if frac <= 0:
                    continue
                sig += frac * float(r["signal_pass_wp80"])
                bkg += frac * float(r["background_pass_wp80"])
                sig_entries += frac * float(r["signal_entries"])
                bkg_entries += frac * float(r["background_entries"])
                if np.isfinite(float(r["score_threshold"])):
                    thresholds.append(float(r["score_threshold"]))
            total = sig + bkg
            purity = sig / total if total > 0 else np.nan
            err_lo, err_hi = wilson_err(sig, total)
            rows.append(
                {
                    "model": model,
                    "model_label": MODEL_STYLE[model]["label"],
                    "centrality": cent,
                    "centrality_label": CENT_LABEL[cent],
                    "et_bin": et_bin,
                    "pt_mid": 0.5 * (et_low + et_high),
                    "truth_purity_wp80": purity,
                    "truth_purity_err_low": err_lo,
                    "truth_purity_err_high": err_hi,
                    "actual_signal_efficiency": sig / sig_entries if sig_entries > 0 else np.nan,
                    "signal_pass_wp80": sig,
                    "background_pass_wp80": bkg,
                    "native_binning": native_binning,
                    "score_threshold": float(np.nanmean(thresholds)) if thresholds else np.nan,
                    "source": str(source.relative_to(REPO)),
                }
            )
    return pd.DataFrame(rows)


def wilson_err(k: float, n: float) -> tuple[float, float]:
    if n <= 0:
        return np.nan, np.nan
    z = 1.0
    p = k / n
    denom = 1.0 + z * z / n
    center = (p + z * z / (2.0 * n)) / denom
    half = z * np.sqrt((p * (1.0 - p) + z * z / (4.0 * n)) / n) / denom
    lo = max(0.0, center - half)
    hi = min(1.0, center + half)
    return p - lo, hi - p


def wp80_for_group(g: pd.DataFrame) -> dict:
    g = g.sort_values("score_low").copy()
    sig = g["signal_count"].to_numpy(dtype=float)
    bkg = g["background_count"].to_numpy(dtype=float)
    sig_total = sig.sum()
    bkg_total = bkg.sum()
    sig_above = sig[::-1].cumsum()[::-1]
    bkg_above = bkg[::-1].cumsum()[::-1]
    total_above = sig_above + bkg_above
    if sig_total <= 0:
        return {
            "score_threshold": np.nan,
            "actual_signal_efficiency": np.nan,
            "truth_purity": np.nan,
            "truth_purity_err_low": np.nan,
            "truth_purity_err_high": np.nan,
            "signal_above": np.nan,
            "background_above": np.nan,
            "signal_entries": sig_total,
            "background_entries": bkg_total,
        }
    eff = sig_above / sig_total
    idx = int(np.nanargmin(np.abs(eff - 0.80)))
    selected = total_above[idx]
    purity = sig_above[idx] / selected if selected > 0 else np.nan
    err_lo, err_hi = wilson_err(sig_above[idx], selected)
    return {
        "score_threshold": float(g["score_low"].iloc[idx]),
        "actual_signal_efficiency": float(eff[idx]),
        "truth_purity": float(purity),
        "truth_purity_err_low": float(err_lo),
        "truth_purity_err_high": float(err_hi),
        "signal_above": float(sig_above[idx]),
        "background_above": float(bkg_above[idx]),
        "signal_entries": float(sig_total),
        "background_entries": float(bkg_total),
    }


def aggregate_selected(rows: list[dict], model: str, model_label: str, cent: str, et_bin: str, native_binning: str) -> dict:
    sig = float(np.nansum([r["signal_above"] for r in rows]))
    bkg = float(np.nansum([r["background_above"] for r in rows]))
    total = sig + bkg
    purity = sig / total if total > 0 else np.nan
    err_lo, err_hi = wilson_err(sig, total)
    sig_entries = float(np.nansum([r["signal_entries"] for r in rows]))
    eff = sig / sig_entries if sig_entries > 0 else np.nan
    return {
        "model": model,
        "model_label": model_label,
        "centrality": cent,
        "centrality_label": CENT_LABEL[cent],
        "et_bin": et_bin,
        "pt_mid": et_mid(et_bin),
        "truth_purity_wp80": purity,
        "truth_purity_err_low": err_lo,
        "truth_purity_err_high": err_hi,
        "actual_signal_efficiency": eff,
        "signal_pass_wp80": sig,
        "background_pass_wp80": bkg,
        "native_binning": native_binning,
        "source": "",
    }


def load_pp_points() -> pd.DataFrame:
    df = read_csv_after_header(PP_HIST, "scope,")
    df = df[(df["scope"] != "integrated") & (df["et_low"] >= 15.0) & (df["et_high"] <= 35.0)].copy()
    selected_rows = []
    for scope, g in df.groupby("scope", sort=False):
        piv = g.pivot_table(
            index=["bin_low", "bin_high"],
            columns="class",
            values="count",
            aggfunc="sum",
            fill_value=0,
        ).reset_index()
        piv = piv.rename(
            columns={
                "bin_low": "score_low",
                "bin_high": "score_high",
                "signal": "signal_count",
                "background": "background_count",
            }
        )
        r = wp80_for_group(piv)
        selected_rows.append(
            {
                "centrality": "0_20",
                "et_low": float(g["et_low"].iloc[0]),
                "et_high": float(g["et_high"].iloc[0]),
                "signal_pass_wp80": r["signal_above"],
                "background_pass_wp80": r["background_above"],
                "signal_entries": r["signal_entries"],
                "background_entries": r["background_entries"],
                "score_threshold": r["score_threshold"],
            }
        )
    pp_native = pd.DataFrame(selected_rows)
    repeated = []
    for cent in CENT_ORDER:
        tmp = pp_native.copy()
        tmp["centrality"] = cent
        repeated.append(tmp)
    selected = pd.concat(repeated, ignore_index=True)
    return aggregate_native_selected_to_display(
        selected,
        "pp_baseline",
        "overlap-rebinned from pp native compact ET bins: 15-18, 18-22, 22-28, 28-35",
        PP_HIST,
    )


def load_ptcent7_best_points() -> pd.DataFrame:
    df = pd.read_csv(PTCENT7_WP80_CELLS)
    selected_rows = []
    for _, r in df.iterrows():
        cent = f"{int(r['centrality_min'])}_{int(r['centrality_max'])}"
        signal_pass = float(r["signal_efficiency"]) * float(r["signal_entries"])
        background_pass = float(r["background_fake_rate"]) * float(r["background_entries"])
        selected_rows.append(
            {
                "centrality": cent,
                "et_low": float(r["pt_min"]),
                "et_high": float(r["pt_max"]),
                "actual_signal_efficiency": float(r["signal_efficiency"]),
                "signal_pass_wp80": signal_pass,
                "background_pass_wp80": background_pass,
                "signal_entries": float(r["signal_entries"]),
                "background_entries": float(r["background_entries"]),
                "score_threshold": float(r["threshold"]),
            }
        )
    selected = pd.DataFrame(selected_rows)
    return aggregate_native_selected_to_display(
        selected,
        "ptcent7_best",
        "native 8-bin ptCent7 WP80 cells; 25-27 and 27-30 combined into display 25-30",
        PTCENT7_WP80_CELLS,
    )


def load_basev3e_points() -> pd.DataFrame:
    df = pd.read_csv(BASEV3E_HIST)
    df = df[df["product"] == "baseBDT_v3E_withCentrality"].copy()
    selected_rows = []
    for (cent_fine, et_bin), g in df.groupby(["centrality_bin", "et_bin"], sort=False):
        r = wp80_for_group(g)
        selected_rows.append(
            {
                "centrality": coarse_cent(cent_fine),
                "et_low": float(g["et_low"].iloc[0]),
                "et_high": float(g["et_high"].iloc[0]),
                "signal_pass_wp80": r["signal_above"],
                "background_pass_wp80": r["background_above"],
                "signal_entries": r["signal_entries"],
                "background_entries": r["background_entries"],
                "score_threshold": r["score_threshold"],
            }
        )
    selected = pd.DataFrame(selected_rows)
    return aggregate_native_selected_to_display(
        selected,
        "basev3e_cent",
        "overlap-rebinned from base v3E native ET bins: 15-20, 20-25, 25-35",
        BASEV3E_HIST,
    )


def style_axis(ax) -> None:
    ax.grid(True, color="#D1D5DB", lw=0.8, alpha=0.75)
    ax.set_axisbelow(True)
    ax.tick_params(direction="in", top=True, right=True, labelsize=11.0)
    for spine in ax.spines.values():
        spine.set_linewidth(1.0)


def sphinx_label(fig, x=0.90, y=0.952, size=14.5) -> None:
    fig.text(x, y, "sPHENIX", ha="right", fontsize=size, fontstyle="italic", fontweight="bold")
    fig.text(x + 0.085, y, "Internal", ha="right", fontsize=size)


def plot(points: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15.7, 5.2), dpi=190, sharey=True)
    for ax, cent in zip(axes, CENT_ORDER):
        sub = points[points["centrality"] == cent].copy()
        for model in ["pp_baseline", "ptcent7_best", "basev3e_cent"]:
            g = sub[sub["model"] == model].sort_values("pt_mid")
            if g.empty:
                continue
            st = MODEL_STYLE[model]
            ax.errorbar(
                g["pt_mid"],
                g["truth_purity_wp80"],
                yerr=[g["truth_purity_err_low"], g["truth_purity_err_high"]],
                xerr=[
                    g["pt_mid"] - g["et_bin"].str.split("-").str[0].astype(float),
                    g["et_bin"].str.split("-").str[1].astype(float) - g["pt_mid"],
                ],
                fmt=st["marker"],
                ms=5.8 if model != "ptcent7_best" else 5.2,
                lw=0.0,
                elinewidth=1.8,
                capsize=2.0,
                color=st["color"],
                linestyle="None",
                label=st["label"] if cent == "0_20" else None,
                zorder=st["zorder"],
                alpha=0.96,
            )
        ax.set_title(f"{CENT_LABEL[cent]} centrality", fontsize=13.8, fontweight="bold", pad=9)
        ax.set_xlim(14.4, 35.6)
        ax.set_ylim(0.82, 1.005)
        ax.set_xticks([0.5 * (lo + hi) for lo, hi, _ in DISPLAY_ET_BINS])
        ax.set_xticklabels([label for _, _, label in DISPLAY_ET_BINS], rotation=35, ha="right")
        ax.set_xlabel(r"Cluster $E_T$ bin [GeV]", fontsize=12.0)
        if ax is axes[0]:
            ax.set_ylabel("Truth purity at WP80", fontsize=12.4)
            ax.legend(frameon=False, fontsize=9.0, loc="lower left")
        style_axis(ax)

    fig.text(
        0.055,
        0.965,
        "Truth purity at the 80% signal-efficiency BDT working point",
        fontsize=17.2,
        fontweight="bold",
    )
    fig.text(
        0.055,
        0.918,
        "All models are displayed in seven common ET bins from 15-35 GeV; pp is repeated as the no-centrality reference.",
        fontsize=10.8,
        color="#4B5563",
    )
    sphinx_label(fig)
    fig.tight_layout(rect=(0.02, 0.02, 0.99, 0.89))
    fig.savefig(OUT_PNG)
    plt.close(fig)


def main() -> int:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    points = pd.concat([load_pp_points(), load_ptcent7_best_points(), load_basev3e_points()], ignore_index=True)
    points = points.replace([np.inf, -np.inf], np.nan).sort_values(["centrality", "model", "pt_mid"])
    points.to_csv(OUT_CSV, index=False)
    plot(points)
    print(OUT_PNG)
    print(OUT_CSV)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
