#!/usr/bin/env python3
"""Make slide-6-style AuAu isolation WP/leakage panels from compact CSVs."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from matplotlib.container import ErrorbarContainer
from matplotlib.legend_handler import HandlerErrorbar


COLORS = {
    0.70: "#1f77b4",
    0.80: "#2ca02c",
    0.90: "#d627b5",
}
MARKERS = {
    0.70: "s",
    0.80: "o",
    0.90: "^",
}


def add_sphenix_label(ax, x=0.58, y=0.94):
    ax.text(
        x,
        y,
        "sPHENIX",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=15,
        fontweight="bold",
        fontstyle="italic",
        family="serif",
    )
    ax.annotate(
        " Internal",
        xy=(x, y),
        xycoords=ax.transAxes,
        xytext=(75, 0),
        textcoords="offset points",
        ha="left",
        va="top",
        fontsize=15,
        family="serif",
    )


def style_axes(ax):
    ax.tick_params(which="both", direction="in", top=True, right=True, length=6, width=1.2)
    ax.tick_params(which="minor", length=3)
    for spine in ax.spines.values():
        spine.set_linewidth(1.2)
    ax.minorticks_on()


def draw_wp_panel(ax, wp: pd.DataFrame, cone_label: str):
    wp = wp.copy()
    wp["plot_value"] = pd.to_numeric(wp["value"], errors="coerce")
    wp["plot_value"] = wp["plot_value"].fillna(pd.to_numeric(wp["cut_gev"], errors="coerce"))
    for eff in (0.70, 0.80, 0.90):
        part = wp[np.isclose(wp["efficiency"].astype(float), eff)].sort_values("cent_center")
        x = part["cent_center"].astype(float).to_numpy()
        y = part["plot_value"].astype(float).to_numpy()
        if len(x) == 0:
            continue

        coeff = np.polyfit(x, y, deg=1)
        xfit = np.linspace(0, 80, 200)
        yfit = np.polyval(coeff, xfit)
        ax.plot(xfit, yfit, color=COLORS[eff], linestyle="--", lw=2.0, alpha=0.95)
        ax.plot(
            x,
            y,
            linestyle="None",
            marker=MARKERS[eff],
            ms=5.6,
            color=COLORS[eff],
            label=f"{int(round(eff * 100))}% Efficiency",
        )

    ax.set_xlim(0, 80)
    ymax = max(6.8, wp["plot_value"].astype(float).max() * 1.16)
    ax.set_ylim(0, ymax)
    ax.set_xlabel("Centrality [%]", fontsize=14, family="serif")
    ax.set_ylabel(r"$E_T^{iso}$ Cutoff [GeV]", fontsize=15, family="serif")
    ax.set_xticks(np.arange(0, 81, 10))
    ax.text(0.08, 0.89, "Photon+Jet 12+20 Embedded SIM", transform=ax.transAxes, fontsize=13, family="serif")
    ax.text(0.08, 0.80, r"$|v_z| < 10$ cm, $|\eta^\gamma| < 0.7$", transform=ax.transAxes, fontsize=13, family="serif")
    ax.text(0.08, 0.71, rf"$\Delta R_{{cone}} < {cone_label}$", transform=ax.transAxes, fontsize=13, family="serif")
    add_sphenix_label(ax, x=0.55, y=0.94)
    ax.legend(loc="upper right", bbox_to_anchor=(0.98, 0.68), frameon=False, fontsize=10)
    style_axes(ax)


def draw_leakage_panel(ax, leak: pd.DataFrame, cone_label: str, wp: pd.DataFrame | None = None):
    cent_order = ["0-20%", "20-50%", "50-80%"]
    cent_colors = {"0-20%": "#111111", "20-50%": "#0000cc", "50-80%": "#ff7f0e"}
    # The compact B-leakage CSV stores the two isolation branches with the
    # historical mode names, but the slide-6 visual convention is what matters
    # here: the upper branch is the fixed Eiso comparison and the lower branch
    # is the centrality-dependent isolation comparison.
    mode_styles = {
        "isSliding": {
            "markerface": "color",
            "label": r"Fixed $E_T^{iso}<4$ GeV",
        },
        "fixedIso4GeV": {
            "markerface": "none",
            "label": r"Centrality-dependent $E_T^{iso}$ cut",
        },
    }
    fit_label = mode_styles["isSliding"]["label"]
    if wp is not None and not wp.empty:
        fit_part = wp[np.isclose(wp["efficiency"].astype(float), 0.90)].copy()
        fit_part["plot_value"] = pd.to_numeric(fit_part["value"], errors="coerce")
        fit_part["plot_value"] = fit_part["plot_value"].fillna(pd.to_numeric(fit_part["cut_gev"], errors="coerce"))
        fit_part = fit_part.dropna(subset=["plot_value", "cent_center"])
        if len(fit_part) >= 2:
            slope, intercept = np.polyfit(
                fit_part["cent_center"].astype(float).to_numpy(),
                fit_part["plot_value"].astype(float).to_numpy(),
                deg=1,
            )
            sign = "-" if slope < 0 else "+"
            fit_label = (
                rf"$E_T^{{iso}}(\mathrm{{cent}}) < {intercept:.2f} "
                rf"{sign} {abs(slope):.4f}[\mathrm{{cent}}\ \%]$"
            )
    centrality_handles = []
    for cent in cent_order:
        cent_handle = ax.errorbar(
            [],
            [],
            yerr=[],
            color=cent_colors[cent],
            marker="o",
            markersize=4.8,
            linestyle="None",
            elinewidth=1.25,
            capsize=3.0,
            label=cent,
        )
        centrality_handles.append(cent_handle)
        for mode in ("isSliding", "fixedIso4GeV"):
            style = mode_styles[mode]
            part = leak[(leak["centrality"] == cent) & (leak["mode"] == mode)].sort_values("pt_center")
            if part.empty:
                continue
            ax.errorbar(
                part["pt_center"].astype(float),
                part["value"].astype(float),
                yerr=part["error"].astype(float),
                color=cent_colors[cent],
                marker="o",
                markersize=(4.8 if style["markerface"] == "color" else 5.2),
                markerfacecolor=(cent_colors[cent] if style["markerface"] == "color" else "white"),
                markeredgecolor=cent_colors[cent],
                linestyle="None",
                linewidth=0.0,
                elinewidth=1.35,
                capsize=3.2,
                capthick=1.35,
                markeredgewidth=(1.25 if style["markerface"] == "color" else 1.65),
                label=f"{cent} {style['label']}",
            )

    ax.set_xlim(15, 35)
    ax.set_ylim(-0.025, 1.0)
    ax.set_xlabel(r"$p_T^\gamma$ [GeV]", fontsize=14, family="serif", loc="right")
    ax.set_ylabel(r"Region B signal leakage, $f_B = B_\mathrm{sig}/A_\mathrm{sig}$", fontsize=14, family="serif")
    ax.set_xticks([15, 20, 25, 30, 35])
    ax.set_yticks(np.arange(0.0, 1.01, 0.1))
    ax.yaxis.set_minor_locator(MultipleLocator(0.05))
    ax.text(
        0.01,
        1.01,
        "B-leakage, BDT selection, Photon+Jet Embedded Pythia (12 + 20) GeV",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=9.0,
        family="serif",
    )
    cent_leg = ax.legend(
        handles=centrality_handles,
        loc="upper left",
        frameon=False,
        fontsize=9.8,
        bbox_to_anchor=(0.055, 1.0),
        handlelength=1.8,
        labelspacing=0.42,
        handler_map={ErrorbarContainer: HandlerErrorbar(xerr_size=0.0, yerr_size=0.55)},
    )
    ax.add_artist(cent_leg)

    iso_handles = []
    for mode in ("isSliding", "fixedIso4GeV"):
        style = mode_styles[mode]
        handle = ax.errorbar(
            [],
            [],
            yerr=[],
            color="#111111",
            marker="o",
            markersize=(4.8 if style["markerface"] == "color" else 5.2),
            markerfacecolor=("#111111" if style["markerface"] == "color" else "white"),
            markeredgecolor="#111111",
            linestyle="None",
            linewidth=0.0,
            elinewidth=1.25,
            capsize=3.0,
            markeredgewidth=(1.25 if style["markerface"] == "color" else 1.65),
            label=(style["label"] if mode == "isSliding" else fit_label),
        )
        iso_handles.append(handle)
    ax.legend(
        handles=iso_handles,
        loc="upper center",
        frameon=False,
        fontsize=9.4,
        bbox_to_anchor=(0.58, 1.0),
        handlelength=2.5,
        labelspacing=0.50,
        handler_map={ErrorbarContainer: HandlerErrorbar(xerr_size=0.0, yerr_size=0.55)},
    )
    ax.text(0.70, 0.77, "sPHENIX", transform=ax.transAxes, fontsize=11.0, fontweight="bold", family="serif")
    ax.text(0.83, 0.77, "Internal", transform=ax.transAxes, fontsize=11.0, fontstyle="italic", family="serif")
    ax.text(0.58, 0.70, r"Pythia Overlay   $\sqrt{s_{NN}} = 200$ GeV", transform=ax.transAxes, fontsize=10.0, family="serif")
    style_axes(ax)


def make_plot(csv_path: Path, output_path: Path, cone_label: str):
    df = pd.read_csv(csv_path)
    wp = df[df["kind"] == "wp"].copy()
    leak = df[df["kind"] == "b_leakage"].copy()

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "dejavuserif",
            "axes.linewidth": 1.2,
        }
    )
    fig, axes = plt.subplots(1, 2, figsize=(14.4, 5.2), dpi=180)
    draw_wp_panel(axes[0], wp, cone_label)
    draw_leakage_panel(axes[1], leak, cone_label, wp=wp)
    fig.subplots_adjust(left=0.07, right=0.985, bottom=0.15, top=0.95, wspace=0.19)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--csv", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--cone-label", required=True)
    args = parser.parse_args()
    make_plot(args.csv, args.out, args.cone_label)
    print(args.out)


if __name__ == "__main__":
    main()
