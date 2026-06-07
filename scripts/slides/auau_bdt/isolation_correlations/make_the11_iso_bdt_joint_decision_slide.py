#!/usr/bin/env python3
"""Build a simple slide joining R=0.4 isolation, BDT WP80, and correlation."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib import font_manager  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from matplotlib.patches import FancyBboxPatch, Rectangle  # noqa: E402


REPO = Path(__file__).resolve().parents[4]
ISO_FINE_CSV = REPO / (
    "dataOutput/auauMLDiagnosticRuns/global_etcent_inclusive3_sixpack_20260516_135439/"
    "slideReady/isolation_wp_checks/"
    "auau_sliding_isolation_wp_r03_r04_eff70_80_90_fine7_no_numbers_coefficients.csv"
)
ISO_LINEAR_CSV = REPO / (
    "dataOutput/auauMLDiagnosticRuns/global_etcent_inclusive3_sixpack_20260516_135439/"
    "slideReady/isolation_wp_checks/"
    "auau_sliding_isolation_wp_r03_r04_eff70_80_90_no_numbers_coefficients.csv"
)
BDT_SUMMARY_CSV = REPO / (
    "dataOutput/auauTightBDTValidation/THE8_branchA_ladder_scorecache_fullstat_20260527/"
    "the41_centdep_wp_slides_20260604/the41_centdep_bdt_wp_summary.csv"
)
BDT_MANIFEST_JSON = REPO / (
    "dataOutput/auauTightBDTValidation/THE8_branchA_ladder_scorecache_fullstat_20260527/"
    "the41_centdep_wp_slides_20260604/the41_centdep_bdt_wp_slide_manifest.json"
)
CORR_CSV = REPO / (
    "dataOutput/auauMLDiagnosticRuns/global_etcent_inclusive3_sixpack_20260516_135439/"
    "slideReady/isolation_feature_correlations/noiso_score_diagnostic/"
    "isolation_feature_score_correlations.csv"
)
OUT_DIR = REPO / (
    "dataOutput/auauMLDiagnosticRuns/global_etcent_inclusive3_sixpack_20260516_135439/"
    "slideReady/isolation_feature_correlations/noiso_score_diagnostic/"
    "the11_iso_bdt_joint_decision_20260604"
)

TARGET_SCORE = "score_globalEtCent1535_bdt_noIso_ptCent7"
CENT_FINE = ["0-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-80"]
ET_ROWS = ["15-20", "20-25", "25-30", "30-35"]
CENT_SCOPES = [("cent_0_20", "0-20%"), ("cent_20_50", "20-50%"), ("cent_50_80", "50-80%")]
CLASS_ORDER = [
    ("all", "all candidates", "#7057a8"),
    ("signal", "truth photons", "#c94b4b"),
    ("background", "inclusive jets", "#3168a8"),
]

W, H, DPI = 2560, 1440, 200
INK = "#111827"
MUTED = "#59616f"
GRID = "#dfe5ec"
PANEL_EDGE = "#d8dee8"
ORANGE = "#d87306"
GREEN = "#1fa77a"
PURPLE = "#7057a8"
BLUE = "#3168a8"
RED = "#c94b4b"
TEAL = "#168a83"


def setup_style() -> None:
    available = {f.name for f in font_manager.fontManager.ttflist}
    for name in ["Times New Roman", "Times", "DejaVu Serif"]:
        if name in available:
            plt.rcParams["font.family"] = name
            break
    plt.rcParams.update(
        {
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "axes.facecolor": "white",
            "axes.edgecolor": INK,
            "axes.linewidth": 1.1,
            "axes.labelcolor": INK,
            "xtick.color": INK,
            "ytick.color": INK,
            "mathtext.fontset": "dejavuserif",
        }
    )


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle))


def load_inputs() -> dict:
    iso_rows = [
        r
        for r in read_csv(ISO_FINE_CSV)
        if r["cone"] == "R0.4" and abs(float(r["efficiency"]) - 0.9) < 1e-9
    ]
    iso_linear = next(
        r
        for r in read_csv(ISO_LINEAR_CSV)
        if r["cone"] == "R0.4" and abs(float(r["efficiency"]) - 0.9) < 1e-9
    )
    bdt_rows = [
        r
        for r in read_csv(BDT_SUMMARY_CSV)
        if r["row_type"] == "flat_fit_constant" and abs(float(r["target_signal_efficiency"]) - 0.8) < 1e-9
    ]
    bdt_fit = json.loads(BDT_MANIFEST_JSON.read_text())["fits"]["WP80"]
    return {
        "iso_rows": iso_rows,
        "iso_linear": iso_linear,
        "bdt_rows": bdt_rows,
        "bdt_fit": bdt_fit,
        "corr_rows": read_csv(CORR_CSV),
    }


def corr_value(rows: list[dict[str, str]], scope: str, klass: str) -> float:
    matches = [
        r
        for r in rows
        if r["scope"] == scope
        and r["class"] == klass
        and r["iso_variable"] == "reco_eiso_r40"
        and r["target_variable"] == TARGET_SCORE
    ]
    if not matches:
        raise KeyError(f"missing R=0.4 correlation for {scope} {klass}")
    return float(matches[0]["pearson"])


def rounded(fig, xywh, facecolor="white", edgecolor=PANEL_EDGE, lw=1.2, radius=0.010) -> None:
    fig.add_artist(
        FancyBboxPatch(
            (xywh[0], xywh[1]),
            xywh[2],
            xywh[3],
            boxstyle=f"round,pad=0.008,rounding_size={radius}",
            transform=fig.transFigure,
            facecolor=facecolor,
            edgecolor=edgecolor,
            linewidth=lw,
            zorder=0.1,
        )
    )


def add_header(fig) -> None:
    fig.text(
        0.045,
        0.955,
        r"Cuts vary with centrality; no $E_T$ slope is applied in 15-35 GeV",
        fontsize=27.5,
        fontweight="bold",
        color=INK,
        ha="left",
        va="top",
    )
    fig.text(
        0.047,
        0.906,
        r"Each row is a cluster-$E_T$ bin. Repeated colors show the same cut is used across $E_T$ after the window check.",
        fontsize=15.3,
        color=MUTED,
        ha="left",
        va="top",
    )
    fig.add_artist(Rectangle((0.045, 0.874), 0.91, 0.003, transform=fig.transFigure, color="#d8dee8", lw=0))


def threshold_map(ax, values, title, formula, cmap, text_color, value_fmt, show_ylabel=True):
    data = np.tile(np.asarray(values, dtype=float), (len(ET_ROWS), 1))
    im = ax.imshow(data, aspect="auto", cmap=cmap)
    ax.set_title(title, fontsize=15.8, fontweight="bold", pad=9)
    ax.text(
        0.025,
        0.925,
        formula,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=11.8,
        color=text_color,
        fontweight="bold",
        bbox={
            "boxstyle": "round,pad=0.24",
            "facecolor": "white",
            "edgecolor": "#e5e7eb",
            "alpha": 0.92,
        },
    )
    ax.set_xticks(np.arange(len(CENT_FINE)))
    ax.set_xticklabels([f"{x}%" for x in CENT_FINE], fontsize=10.2)
    ax.set_yticks(np.arange(len(ET_ROWS)))
    ax.set_yticklabels([f"{x}" for x in ET_ROWS], fontsize=10.7)
    ax.set_xlabel("")
    if show_ylabel:
        ax.set_ylabel(r"cluster $E_T$ bin [GeV]", fontsize=11.8)
    else:
        ax.set_ylabel("")
    ax.set_xticks(np.arange(-0.5, len(CENT_FINE), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(ET_ROWS), 1), minor=True)
    ax.grid(which="minor", color="white", lw=1.4)
    ax.tick_params(which="minor", bottom=False, left=False)
    for col, value in enumerate(values):
        ax.text(
            col,
            1.5,
            value_fmt(value),
            ha="center",
            va="center",
            fontsize=13.3,
            fontweight="bold",
            color="white",
            bbox={"boxstyle": "round,pad=0.18", "facecolor": text_color, "edgecolor": "none", "alpha": 0.88},
        )
    return im


def add_map_note(fig) -> None:
    rounded(fig, [0.365, 0.494, 0.270, 0.058], facecolor="#f7fafc")
    fig.text(
        0.500,
        0.523,
        r"No extra $E_T$ slope is applied in this working window: the derived thresholds are centrality-only.",
        fontsize=11.3,
        color=MUTED,
        ha="center",
        va="center",
    )


def draw_correlation(ax, corr_rows: list[dict[str, str]]) -> dict:
    x = np.arange(len(CENT_SCOPES), dtype=float)
    width = 0.22
    stats = {}
    for idx, (klass, label, color) in enumerate(CLASS_ORDER):
        vals = [corr_value(corr_rows, scope, klass) for scope, _ in CENT_SCOPES]
        stats[klass] = vals
        bars = ax.bar(x + (idx - 1) * width, vals, width=width, color=color, label=label, zorder=3)
        for bar, val in zip(bars, vals):
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                val - 0.028 if val < -0.08 else val + 0.010,
                f"{val:.2f}",
                ha="center",
                va="top" if val < -0.08 else "bottom",
                color="white" if val < -0.08 else INK,
                fontsize=10.8,
                fontweight="bold",
            )
    ax.axhline(0, color="#6f7682", lw=1.1)
    ax.axhspan(-0.07, 0.03, color="#f5e6e6", zorder=0)
    ax.set_ylim(-0.56, 0.06)
    ax.set_xticks(x)
    ax.set_xticklabels([label for _, label in CENT_SCOPES], fontsize=11.6)
    ax.set_ylabel(r"Pearson $\rho$ with BDT score", fontsize=12.2)
    ax.set_title(r"R=0.4 isolation correlation with BDT score", fontsize=15.0, fontweight="bold", pad=8)
    ax.grid(True, axis="y", color=GRID, lw=0.8, zorder=0)
    ax.legend(loc="lower left", ncol=3, frameon=False, fontsize=10.7)
    return stats


def add_readout(fig, iso_values, bdt_values, corr_stats) -> None:
    rounded(fig, [0.700, 0.105, 0.245, 0.310], facecolor="#f7fafc")
    fig.text(0.724, 0.376, "Fast read", fontsize=18.2, fontweight="bold", color=INK, ha="left", va="top")
    y = 0.326
    rows = [
        (ORANGE, rf"Isolation cut falls: {iso_values[0]:.2f} $\rightarrow$ {iso_values[-1]:.2f} GeV"),
        (GREEN, rf"BDT cut rises: {bdt_values[0]:.3f} $\rightarrow$ {bdt_values[-1]:.3f}"),
        (INK, r"$E_T$: repeated rows = no applied slope"),
        (RED, r"Truth photons: $|\rho|<0.06$"),
        (BLUE, rf"Inclusive jets: $\rho$ {max(corr_stats['background']):.2f} to {min(corr_stats['background']):.2f}"),
        (PURPLE, rf"Mixture: $\rho$ reaches {min(corr_stats['all']):.2f}"),
    ]
    for color, text in rows:
        fig.text(0.724, y, text, fontsize=12.4, color=color, ha="left", va="top", fontweight="bold" if color in [ORANGE, GREEN] else "normal")
        y -= 0.041


def add_footer(fig) -> None:
    cards = [
        ("Isolation", r"R=0.4 configured WP90: pass if reco $E_T^{iso}<I_{90}(c)$.", ORANGE),
        ("BDT", r"THE-41 WP80: pass if BDT score $>T_{80}(c)$.", GREEN),
        ("Interpretation", "The BDT-isolation trend is visible for inclusive jets, not truth photons.", PURPLE),
    ]
    x0, y, h, gap = 0.055, 0.045, 0.110, 0.018
    w = (0.89 - 2 * gap) / 3
    for idx, (title, body, color) in enumerate(cards):
        x = x0 + idx * (w + gap)
        rounded(fig, [x, y, w, h])
        fig.text(x + 0.018, y + h - 0.028, title, fontsize=14.8, fontweight="bold", color=color, ha="left", va="top")
        fig.text(x + 0.018, y + h - 0.064, body, fontsize=11.2, color=INK, ha="left", va="top")


def write_script(out: Path, iso_a: float, iso_b: float, bdt_a: float, bdt_b: float, corr_stats: dict) -> Path:
    path = out / "the11_iso_bdt_joint_decision_script.md"
    text = f"""# THE-11 Script - Simple Isolation and BDT Decision View

This slide is a compact way to read the isolation and BDT choices together.

The two heatmaps use centrality on the horizontal axis and the cluster E T bins from 15 to 35 GeV on the vertical axis. The repeated color down each column is intentional: after the E T-bin checks, the proposed thresholds are centrality-only in this working window.

On the left, the R equals 0.4 isolation cut decreases with centrality percentile. The fitted threshold is I90 of c equals {iso_a:.2f} {iso_b:+.4f} c GeV. On the right, the THE-41 BDT threshold increases with centrality percentile. The fitted threshold is T80 of c equals {bdt_a:.4f} {bdt_b:+.5f} c.

The bottom strip is the physical insight. For R equals 0.4, truth photons stay close to zero correlation with the no-isolation BDT score, while inclusive jets carry the visible negative correlation. So the BDT is not simply replaying isolation; isolation is mostly exposing a background-tail axis, while the BDT remains a photon-ID score.
"""
    path.write_text(text)
    return path


def render() -> dict:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    setup_style()
    inputs = load_inputs()

    iso_values = [float(r["threshold_gev"]) for r in inputs["iso_rows"]]
    bdt_values = [float(r["threshold"]) for r in inputs["bdt_rows"]]
    iso_a = float(inputs["iso_linear"]["aGeV"])
    iso_b = float(inputs["iso_linear"]["bPerGeV"])
    bdt_a = float(inputs["bdt_fit"]["intercept"])
    bdt_b = float(inputs["bdt_fit"]["slope"])

    orange_cmap = LinearSegmentedColormap.from_list("iso_orange", ["#fff7ed", "#fdba74", "#c2410c"])
    green_cmap = LinearSegmentedColormap.from_list("bdt_green", ["#ecfdf5", "#6ee7b7", "#047857"])

    fig = plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI)
    add_header(fig)

    ax_iso = fig.add_axes([0.060, 0.500, 0.430, 0.290])
    threshold_map(
        ax_iso,
        iso_values,
        r"R=0.4 isolation cut by centrality",
        rf"$I_{{90}}(c)={iso_a:.2f}{iso_b:+.4f}c$ GeV",
        orange_cmap,
        ORANGE,
        lambda v: f"{v:.2f}",
    )

    ax_bdt = fig.add_axes([0.530, 0.500, 0.430, 0.290])
    threshold_map(
        ax_bdt,
        bdt_values,
        "BDT score cut by centrality",
        rf"$T_{{80}}(c)={bdt_a:.3f}{bdt_b:+.5f}c$",
        green_cmap,
        GREEN,
        lambda v: f"{v:.3f}",
        show_ylabel=False,
    )

    ax_corr = fig.add_axes([0.070, 0.105, 0.610, 0.310])
    corr_stats = draw_correlation(ax_corr, inputs["corr_rows"])
    add_readout(fig, iso_values, bdt_values, corr_stats)

    png = OUT_DIR / "the11_iso_bdt_joint_decision_slide.png"
    fig.savefig(png, dpi=DPI)
    plt.close(fig)

    script = write_script(OUT_DIR, iso_a, iso_b, bdt_a, bdt_b, corr_stats)
    manifest = {
        "schema": "THE11_ISO_BDT_JOINT_DECISION_V2",
        "png": str(png),
        "script": str(script),
        "google_slides_mutated": False,
        "sources": {
            "r04_sliding_isolation_fine_csv": str(ISO_FINE_CSV),
            "r04_sliding_isolation_linear_csv": str(ISO_LINEAR_CSV),
            "bdt_wp_summary_csv": str(BDT_SUMMARY_CSV),
            "bdt_wp_manifest_json": str(BDT_MANIFEST_JSON),
            "r04_correlation_csv": str(CORR_CSV),
        },
        "definitions": {
            "isolation": {
                "cone": "R0.4",
                "target_signal_efficiency": 0.90,
                "pass_rule": "reco E_T^iso(Delta R < 0.4) < I90(c)",
                "intercept_GeV": iso_a,
                "slope_GeV_per_percent": iso_b,
                "fine_thresholds_GeV": iso_values,
            },
            "bdt": {
                "target_signal_efficiency": 0.80,
                "pass_rule": "BDT score > T80(c)",
                "intercept": bdt_a,
                "slope_per_percent": bdt_b,
                "rms_residual": float(inputs["bdt_fit"]["rms_residual"]),
                "max_abs_residual": float(inputs["bdt_fit"]["max_abs_residual"]),
                "fine_thresholds": bdt_values,
            },
        },
        "r04_pearson": corr_stats,
        "caveat": "The heatmap rows show the checked 15-35 GeV cluster-E_T bins; the final displayed thresholds are centrality-only in this window. R=0.4 relationship uses exact aggregate correlations because local row-level R=0.4 density shards were not present.",
    }
    manifest_path = OUT_DIR / "the11_iso_bdt_joint_decision_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return manifest


if __name__ == "__main__":
    render()
