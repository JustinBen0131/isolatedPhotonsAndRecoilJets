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
import csv
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize
from matplotlib.patches import FancyBboxPatch, Rectangle


FEATURES = [
    ("e22_over_e37", r"$E_{22}/E_{37}$"),
    ("e22_over_e53", r"$E_{22}/E_{53}$"),
    ("e11_over_e33", r"$E_{11}/E_{33}$"),
    ("e32_over_e35", r"$E_{32}/E_{35}$"),
]
ISO_COLUMNS = [
    ("reco_eiso_r30", r"$R=0.3\ E_T^{iso}$"),
    ("reco_eiso_r40", r"$R=0.4\ E_T^{iso}$"),
]


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input-json", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument(
        "--output-name",
        default="slide20_eiso_shower_shape_th2_correlation.png",
    )
    ap.add_argument(
        "--qualitative",
        action="store_true",
        help="Render as a qualitative TH2 density slide without Pearson-r panel labels.",
    )
    return ap.parse_args()


def add_text(
    fig: plt.Figure,
    x: float,
    y: float,
    text: str,
    *,
    size: float,
    weight: str = "normal",
    color: str = "#111827",
    ha: str = "left",
    va: str = "top",
    linespacing: float = 1.08,
) -> None:
    fig.text(
        x,
        y,
        text,
        fontsize=size,
        fontweight=weight,
        color=color,
        ha=ha,
        va=va,
        linespacing=linespacing,
    )


def add_box(
    fig: plt.Figure,
    xy: tuple[float, float],
    wh: tuple[float, float],
    face: str,
    *,
    edge: str = "none",
    lw: float = 0.0,
    radius: float = 0.018,
    zorder: int = -1,
) -> None:
    fig.patches.append(
        FancyBboxPatch(
            xy,
            wh[0],
            wh[1],
            boxstyle=f"round,pad=0.010,rounding_size={radius}",
            transform=fig.transFigure,
            facecolor=face,
            edgecolor=edge,
            linewidth=lw,
            zorder=zorder,
        )
    )


def write_summary_csv(payload: dict, path: Path) -> None:
    rows = []
    for feature, feature_label in FEATURES:
        for iso_col, iso_label in ISO_COLUMNS:
            panel = payload["panels"][iso_col][feature]
            corr = panel["correlations"]
            rows.append(
                {
                    "iso_variable": iso_col,
                    "iso_label": iso_label,
                    "feature": feature,
                    "feature_label": feature_label,
                    "pearson_all": corr["all"]["pearson"],
                    "pearson_signal": corr["signal"]["pearson"],
                    "pearson_background": corr["background"]["pearson"],
                    "n_all": corr["all"]["n"],
                    "n_signal": corr["signal"]["n"],
                    "n_background": corr["background"]["n"],
                    "in_range_fraction": panel["in_range_fraction"],
                }
            )
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def fmt_corr(value: float | str | None) -> str:
    try:
        val = float(value)
    except Exception:
        return "nan"
    return f"{val:+.2f}" if math.isfinite(val) else "nan"


def fmt_count(value: int | float | str) -> str:
    try:
        val = int(float(value))
    except Exception:
        return "0"
    if val >= 1_000_000:
        return f"{val / 1_000_000:.2f}M"
    if val >= 1_000:
        return f"{val / 1_000:.0f}k"
    return str(val)


def render(payload: dict, out_png: Path, *, qualitative: bool = False) -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "dejavuserif",
            "axes.unicode_minus": False,
        }
    )

    navy = "#101828"
    gray = "#475467"
    border = "#D0D5DD"
    pale_blue = "#EAF3FF"
    pale_orange = "#FFF1E7"
    yellow = "#FFF2B8"

    fig = plt.figure(figsize=(16, 9), dpi=180, facecolor="white")
    fig.patches.append(Rectangle((0, 0), 1, 1, transform=fig.transFigure, facecolor="white", zorder=-5))

    add_text(
        fig,
        0.035,
        0.955,
        r"Raw cone isolation versus shower-shape variables",
        size=27.5,
        weight="bold",
        color=navy,
    )
    add_text(
        fig,
        0.035,
        0.912,
        r"Centrality-integrated SIM candidates, $15 < cluster\ E_T < 35$ GeV; each row shows the same-cluster TH2 density.",
        size=14.2,
        color=gray,
    )

    add_box(fig, (0.045, 0.158), (0.910, 0.710), "white", edge=border, lw=1.0, radius=0.012)
    add_box(fig, (0.118, 0.812), (0.366, 0.040), pale_blue, edge="#BBD7F5", lw=0.8, radius=0.010)
    add_box(fig, (0.532, 0.812), (0.366, 0.040), pale_orange, edge="#F3C7A5", lw=0.8, radius=0.010)
    add_text(fig, 0.301, 0.841, r"$R=0.3\ E_T^{iso}$", size=17, weight="bold", color=navy, ha="center")
    add_text(fig, 0.715, 0.841, r"$R=0.4\ E_T^{iso}$", size=17, weight="bold", color=navy, ha="center")
    add_text(
        fig,
        0.500,
        0.792,
        (
            r"Each panel is a TH2 occupancy map; color is $\log_{10}(N_{\mathrm{bin}}+1)$."
            if qualitative
            else r"Panel labels show Pearson $r$ for all / signal / background; heatmap color is $\log_{10}(N_{\mathrm{bin}}+1)$."
        ),
        size=11.5,
        weight="bold",
        color=gray,
        ha="center",
    )

    left_label_x = 0.062
    left, right = 0.127, 0.925
    bottom, top = 0.245, 0.760
    col_gap = 0.055
    row_gap = 0.026
    ncols, nrows = 2, 4
    ax_w = (right - left - col_gap) / 2
    ax_h = (top - bottom - row_gap * (nrows - 1)) / nrows
    vmax = max(
        float(np.nanmax(np.asarray(payload["panels"][iso][feature]["log_counts"], dtype=float)))
        for iso, _ in ISO_COLUMNS
        for feature, _ in FEATURES
    )
    norm = Normalize(vmin=0.0, vmax=max(vmax, 1.0))
    cmap = "viridis"
    first_im = None

    for r, (feature, feature_label) in enumerate(FEATURES):
        y0 = top - (r + 1) * ax_h - r * row_gap
        add_text(
            fig,
            left_label_x,
            y0 + ax_h * 0.60,
            feature_label,
            size=17,
            weight="bold",
            color=navy,
            ha="center",
            va="center",
        )
        for c, (iso_col, _iso_label) in enumerate(ISO_COLUMNS):
            x0 = left + c * (ax_w + col_gap)
            ax = fig.add_axes([x0, y0, ax_w, ax_h])
            panel = payload["panels"][iso_col][feature]
            xedges = np.asarray(panel["x_edges"], dtype=float)
            yedges = np.asarray(panel["y_edges"], dtype=float)
            z = np.asarray(panel["log_counts"], dtype=float).T
            im = ax.pcolormesh(xedges, yedges, z, shading="auto", cmap=cmap, norm=norm)
            first_im = im
            ax.grid(False)
            ax.tick_params(labelsize=8.2, direction="in", top=True, right=True, length=3)
            if r == nrows - 1:
                ax.set_xlabel(r"$E_T^{iso}$ [GeV]", fontsize=10.5)
            else:
                ax.set_xticklabels([])
            if c == 0:
                ax.set_ylabel("feature value", fontsize=10.5)
            else:
                ax.set_yticklabels([])
            corr = panel["correlations"]
            label = (
                f"N={fmt_count(corr['all']['n'])}"
                if qualitative
                else (
                    f"r all {fmt_corr(corr['all']['pearson'])}   "
                    f"sig {fmt_corr(corr['signal']['pearson'])}   "
                    f"bkg {fmt_corr(corr['background']['pearson'])}\n"
                    f"N={fmt_count(corr['all']['n'])}"
                )
            )
            ax.text(
                0.025,
                0.925,
                label,
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=9.2 if qualitative else 8.8,
                fontweight="bold",
                color=navy,
                bbox=dict(facecolor="white", edgecolor="#EAECF0", alpha=0.92, boxstyle="round,pad=0.22"),
            )

    cax = fig.add_axes([0.935, bottom, 0.012, top - bottom])
    cbar = fig.colorbar(first_im, cax=cax)
    cbar.set_label(r"$\log_{10}(N_{\mathrm{bin}}+1)$", fontsize=10.5)
    cbar.ax.tick_params(labelsize=8.5)

    all_abs = [
        abs(float(payload["panels"][iso][feature]["correlations"]["all"]["pearson"]))
        for iso, _ in ISO_COLUMNS
        for feature, _ in FEATURES
    ]
    max_abs = max(all_abs) if all_abs else math.nan
    takeaway = (
        r"Takeaway    The isolation cones broaden the same shower-shape bands; this is a visual same-cluster TH2 check, not a coefficient test."
        if qualitative
        else (
            r"Takeaway    Isolation-shape coupling is modest and mostly background-driven; "
            rf"signal $r \simeq 0$ while max all-candidate $|r|={max_abs:.2f}$. "
            r"Raw $E_T^{iso}$ is not simply duplicating the shower-shape inputs."
        )
    )
    add_box(fig, (0.065, 0.055), (0.870, 0.065), yellow, edge="none", radius=0.014)
    add_text(fig, 0.092, 0.098, takeaway, size=12.7, weight="bold", color=navy)

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=180)
    plt.close(fig)


def main() -> int:
    args = parse_args()
    payload = json.loads(args.input_json.read_text())
    args.outdir.mkdir(parents=True, exist_ok=True)
    out_png = args.outdir / args.output_name
    render(payload, out_png, qualitative=args.qualitative)
    write_summary_csv(payload, args.outdir / "slide20_eiso_shower_shape_th2_correlations.csv")
    summary_out = args.outdir / "slide20_eiso_shower_shape_th2_summary.json"
    summary_out.write_text(json.dumps(payload["summary"], indent=2, sort_keys=True) + "\n")
    print(out_png)
    print(summary_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
