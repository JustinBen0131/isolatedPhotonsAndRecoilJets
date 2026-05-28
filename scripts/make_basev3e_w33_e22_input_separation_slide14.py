#!/usr/bin/env python3
"""Build a full-slide replacement PNG for the slide-14 input diagnostic."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager
from matplotlib.patches import Patch


ROOT = Path(
    "dataOutput/auauMLDiagnosticRuns/"
    "global_etcent_inclusive3_sixpack_20260516_135439"
)
FEATURE_SUMMARY = (
    ROOT
    / "validation/basev3e_w33_e22ratio_20260518_192630/validation/validation_feature_summary.csv"
)
OUT_DIR = ROOT / "slideReady/basev3e_w33_feature_diagnostics"

PT_SCOPES = ("pt_15_20", "pt_20_25", "pt_25_35")

FEATURE_ORDER = [
    "e22_over_e37",
    "e22_over_e53",
    "cluster_weta_cogx",
    "cluster_wphi_cogx",
    "cluster_wphi33_cogx",
    "cluster_weta33_cogx",
    "cluster_et1",
    "e32_over_e35",
    "e11_over_e33",
    "cluster_et3",
    "cluster_et2",
    "cluster_et4",
    "cluster_Et",
    "centrality",
    "cluster_Eta",
    "vertexz",
]

FEATURE_LABELS = {
    "centrality": "centrality",
    "cluster_Et": r"cluster $E_T$",
    "cluster_Eta": r"cluster $\eta$",
    "cluster_et1": "cluster et1",
    "cluster_et2": "cluster et2",
    "cluster_et3": "cluster et3",
    "cluster_et4": "cluster et4",
    "cluster_weta_cogx": r"base $w_{\eta}$",
    "cluster_wphi_cogx": r"base $w_{\phi}$",
    "cluster_weta33_cogx": r"local $w_{\eta}^{3\times3}$",
    "cluster_wphi33_cogx": r"local $w_{\phi}^{3\times3}$",
    "e11_over_e33": r"$E_{11}/E_{33}$",
    "e32_over_e35": r"$E_{32}/E_{35}$",
    "e22_over_e37": r"$E_{22}/E_{37}$",
    "e22_over_e53": r"$E_{22}/E_{53}$",
    "vertexz": r"$z_{\rm vtx}$",
}

HIGHLIGHT_FEATURES = {
    "e22_over_e37",
    "e22_over_e53",
    "cluster_weta33_cogx",
    "cluster_wphi33_cogx",
}

BLUE = "#1F77B4"
ORANGE = "#D95F02"
HIGHLIGHT_YELLOW = "#FFF2A8"
INK = "#111827"
MUTED = "#374151"
GRID = "#D1D5DB"
TIMES_FONT_FILES = [
    Path("/System/Library/Fonts/Supplemental/Times New Roman.ttf"),
    Path("/System/Library/Fonts/Supplemental/Times New Roman Bold.ttf"),
    Path("/System/Library/Fonts/Supplemental/Times New Roman Italic.ttf"),
    Path("/System/Library/Fonts/Supplemental/Times New Roman Bold Italic.ttf"),
]
TIMES_FAMILY = "Times New Roman"


def register_times_new_roman() -> None:
    missing = [path for path in TIMES_FONT_FILES if not path.exists()]
    if missing:
        raise RuntimeError(f"Missing Times New Roman font files: {missing}")
    for path in TIMES_FONT_FILES:
        font_manager.fontManager.addfont(str(path))


def combine_stats(part: pd.DataFrame) -> tuple[int, float, float]:
    entries = int(part["entries"].sum())
    if entries <= 0:
        return 0, np.nan, np.nan
    mean = float((part["entries"] * part["mean"]).sum() / entries)
    second_moment = float((part["entries"] * (part["std"] ** 2 + part["mean"] ** 2)).sum() / entries)
    variance = max(0.0, second_moment - mean**2)
    return entries, mean, float(np.sqrt(variance))


def build_table() -> pd.DataFrame:
    raw = pd.read_csv(FEATURE_SUMMARY)
    rows = []
    for feature in FEATURE_ORDER:
        sig = raw[
            raw["feature"].eq(feature)
            & raw["scope"].isin(PT_SCOPES)
            & raw["class"].eq("signal")
        ]
        bkg = raw[
            raw["feature"].eq(feature)
            & raw["scope"].isin(PT_SCOPES)
            & raw["class"].eq("background")
        ]
        if sig.empty or bkg.empty:
            raise RuntimeError(f"Missing signal/background feature summary rows for {feature}")
        sig_entries, sig_mean, sig_std = combine_stats(sig)
        bkg_entries, bkg_mean, bkg_std = combine_stats(bkg)
        pooled_width = float(np.sqrt(max(1.0e-12, 0.5 * (sig_std**2 + bkg_std**2))))
        separation = (sig_mean - bkg_mean) / pooled_width
        rows.append(
            {
                "feature": feature,
                "label": FEATURE_LABELS[feature],
                "signal_entries": sig_entries,
                "background_entries": bkg_entries,
                "signal_mean": sig_mean,
                "background_mean": bkg_mean,
                "signal_std": sig_std,
                "background_std": bkg_std,
                "standardized_mean_difference": separation,
                "abs_standardized_mean_difference": abs(separation),
                "signal_higher": separation >= 0,
                "highlighted_from_prior_slide": feature in HIGHLIGHT_FEATURES,
            }
        )
    return pd.DataFrame(rows)


def add_slide_text(fig: plt.Figure) -> None:
    title_font = {"fontfamily": TIMES_FAMILY, "color": INK}
    body_font = {"fontfamily": TIMES_FAMILY, "color": "black"}

    fig.text(
        0.012,
        0.966,
        "Shower shape drives separation",
        ha="left",
        va="top",
        fontsize=36.0,
        fontweight="bold",
        **title_font,
    )

    fig.text(
        0.018,
        0.842,
        "Pre-BDT input diagnostic: for each model input,\n"
        "ask how different signal and background look\n"
        "before the BDT combines variables.",
        ha="left",
        va="top",
        fontsize=19.2,
        linespacing=1.16,
        **body_font,
    )

    fig.text(
        0.018,
        0.692,
        "It compares the signal mean to the background mean,\n"
        "normalized by the combined distribution width.",
        ha="left",
        va="top",
        fontsize=19.2,
        linespacing=1.18,
        **body_font,
    )

    fig.text(
        0.018,
        0.570,
        r"$d = (\mu_{\rm sig}-\mu_{\rm bkg})/\sqrt{(\sigma_{\rm sig}^{2}+\sigma_{\rm bkg}^{2})/2}$",
        ha="left",
        va="top",
        fontsize=18.2,
        **body_font,
    )

    fig.text(
        0.018,
        0.472,
        "Blue: signal average is higher.\n"
        "Orange: background average is higher.\n"
        r"Longer bar: stronger one-variable separation $|d|$.",
        ha="left",
        va="top",
        fontsize=18.6,
        linespacing=1.20,
        **body_font,
    )

    fig.text(
        0.018,
        0.326,
        "Full input contract shown",
        ha="left",
        va="top",
        fontsize=18.0,
        fontweight="bold",
        **title_font,
    )
    fig.text(
        0.018,
        0.286,
        r"base v3E + centrality + $w_{\eta,\phi}^{3\times3}$ + E22 ratios"
        "\n"
        r"Yellow highlight: local 3x3 widths and E22 ratios",
        ha="left",
        va="top",
        fontsize=16.6,
        linespacing=1.25,
        **body_font,
    )
    fig.text(
        0.018,
        0.200,
        r"$E_T$, centrality, $\eta$, and $z_{\rm vtx}$ are context inputs;"
        "\n"
        "their bars are unweighted validation mean differences.",
        ha="left",
        va="top",
        fontsize=15.6,
        linespacing=1.22,
        **body_font,
    )

    fig.text(
        0.908,
        0.940,
        "sPHENIX",
        ha="right",
        va="top",
        fontsize=18.5,
        fontstyle="italic",
        fontweight="bold",
        color=INK,
    )
    fig.text(
        0.912,
        0.940,
        " Internal",
        ha="left",
        va="top",
        fontsize=18.5,
        color=INK,
    )


def add_plot(fig: plt.Figure, df: pd.DataFrame) -> None:
    ax = fig.add_axes([0.505, 0.155, 0.470, 0.705])

    y = np.arange(len(df))
    signed = df["standardized_mean_difference"].to_numpy()
    strength = df["abs_standardized_mean_difference"].to_numpy()
    colors = np.where(signed >= 0, BLUE, ORANGE)

    for idx in y:
        if idx % 2 == 0:
            ax.axhspan(idx - 0.43, idx + 0.43, color="#F8FAFC", zorder=0)

    bars = ax.barh(
        y,
        strength,
        height=0.58,
        color=colors,
        edgecolor="white",
        linewidth=0.8,
        zorder=2,
    )

    ax.set_xlim(0.0, 1.12)
    ax.set_ylim(-0.55, len(df) - 0.45)
    ax.invert_yaxis()
    ax.set_yticks(y)
    tick_labels = ax.set_yticklabels(df["label"], fontsize=12.4, fontfamily=TIMES_FAMILY)
    for tick_label, highlighted in zip(tick_labels, df["highlighted_from_prior_slide"]):
        if highlighted:
            tick_label.set_bbox(
                {
                    "facecolor": HIGHLIGHT_YELLOW,
                    "edgecolor": "none",
                    "boxstyle": "round,pad=0.16",
                    "alpha": 0.95,
                }
            )
    ax.tick_params(axis="y", length=0, pad=5)
    ax.tick_params(axis="x", labelsize=12.0)
    ax.set_xlabel(r"One-variable separation strength $|d|$", fontsize=14.0, labelpad=5)
    ax.grid(axis="x", color=GRID, linewidth=1.0, zorder=1)
    ax.set_axisbelow(True)

    for spine in ["top", "right", "left"]:
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_color(INK)
    ax.spines["bottom"].set_linewidth(1.1)

    for bar, (_, row) in zip(bars, df.iterrows()):
        value = float(row["abs_standardized_mean_difference"])
        ax.text(
            min(value + 0.025, 1.085),
            bar.get_y() + bar.get_height() / 2,
            rf"{value:.2f}",
            ha="left",
            va="center",
            fontsize=9.9,
            fontweight="bold",
            color=INK,
            zorder=4,
        )

    ax.text(
        0.0,
        1.065,
        r"Full input set: base v3E + centrality + local 3x3 widths + E22 ratios",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=11.2,
        color=MUTED,
        fontfamily=TIMES_FAMILY,
    )
    ax.text(
        0.0,
        1.020,
        r"Photon12+20 signal vs Jet12+20+30 inclusive background, $15<E_T<35$ GeV",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=12.0,
        color=MUTED,
        fontfamily=TIMES_FAMILY,
    )

    legend_handles = [
        Patch(facecolor=BLUE, edgecolor="white", label="signal mean higher"),
        Patch(facecolor=ORANGE, edgecolor="white", label="background mean higher"),
    ]
    ax.legend(
        handles=legend_handles,
        loc="upper right",
        bbox_to_anchor=(1.0, -0.105),
        ncol=1,
        frameon=False,
        fontsize=10.7,
        handlelength=1.5,
        borderaxespad=0.0,
    )


def plot(df: pd.DataFrame) -> Path:
    register_times_new_roman()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_png = OUT_DIR / "basev3e_w33_e22_input_separation_slide14_clean.png"
    out_csv = OUT_DIR / "basev3e_w33_e22_input_separation_slide14_clean.csv"
    df.to_csv(out_csv, index=False)

    plt.rcParams.update(
        {
            "font.family": TIMES_FAMILY,
            "font.serif": [TIMES_FAMILY],
            "mathtext.fontset": "custom",
            "mathtext.rm": TIMES_FAMILY,
            "mathtext.it": f"{TIMES_FAMILY}:italic",
            "mathtext.bf": f"{TIMES_FAMILY}:bold",
            "axes.unicode_minus": False,
        }
    )

    fig = plt.figure(figsize=(16.0, 9.0), dpi=200)
    fig.patch.set_facecolor("white")

    add_slide_text(fig)
    add_plot(fig, df)

    fig.savefig(out_png, facecolor="white")
    plt.close(fig)
    print(out_png)
    print(out_csv)
    return out_png


def main() -> None:
    df = build_table()
    plot(df)


if __name__ == "__main__":
    main()
