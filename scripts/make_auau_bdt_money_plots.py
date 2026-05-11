#!/usr/bin/env python3
"""Make slide-candidate plots from expanded AuAu tight-BDT validation outputs."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


KEY_PRODUCTS = [
    "centINDcontrol_allRange",
    "centAsFeat_allRange",
    "centAsFeatMinOpt_pt5to40",
    "centDepBDTs_allRange",
    "centDepFineBDTs_allRange",
]

GRANULAR_PRODUCTS = [
    "ptBinCentAsFeat",
    "ptCentDep3",
    "ptCentDepFine",
]

ALL_PRODUCTS = KEY_PRODUCTS + GRANULAR_PRODUCTS

MODEL_COMPLEXITY = {
    "centINDcontrol_allRange": 1,
    "centAsFeat_allRange": 1,
    "centAsFeatMinOpt_pt5to40": 1,
    "centDepBDTs_allRange": 3,
    "centDepFineBDTs_allRange": 7,
    "ptBinCentAsFeat": 11,
    "ptCentDep3": 33,
    "ptCentDepFine": 77,
}

SHORT_LABELS = {
    "centINDcontrol_allRange": "no cent.",
    "centAsFeat_allRange": "cent input",
    "centAsFeatMinOpt_pt5to40": "minority opt.",
    "centDepBDTs_allRange": "3 cent.",
    "centDepFineBDTs_allRange": "7 cent.",
    "ptBinCentAsFeat": r"$E_T$ bin",
    "ptCentDep3": r"$E_T$x3c",
    "ptCentDepFine": r"$E_T$xfine",
}

CELL_LABELS = {
    "centINDcontrol_allRange": "No cent.\ninput",
    "centAsFeat_allRange": "Cent. as\ninput",
    "centAsFeatMinOpt_pt5to40": "Minority\nbalanced",
    "centDepBDTs_allRange": "3 cent.\nBDTs",
    "centDepFineBDTs_allRange": "7 cent.\nBDTs",
    "ptBinCentAsFeat": "$E_T$ bin\n+ cent.",
    "ptCentDep3": "$E_T$\n+ 3 cent.",
    "ptCentDepFine": "$E_T$\n+ 7 cent.",
}

LABELS = {
    "centINDcontrol_allRange": "No centrality input",
    "centAsFeat_allRange": "Centrality input",
    "centAsFeatMinOpt_pt5to40": "Minority optimized",
    "centDepBDTs_allRange": "3 centrality BDTs",
    "centDepFineBDTs_allRange": "7 centrality BDTs",
    "ptBinCentAsFeat": r"$E_T$-binned, cent input",
    "ptCentDep3": r"$E_T$ x 3 centrality",
    "ptCentDepFine": r"$E_T$ x fine centrality",
}

COLORS = {
    "centINDcontrol_allRange": "#0072B2",
    "centAsFeat_allRange": "#D55E00",
    "centAsFeatMinOpt_pt5to40": "#CC79A7",
    "centDepBDTs_allRange": "#009E73",
    "centDepFineBDTs_allRange": "#E45756",
    "ptBinCentAsFeat": "#56B4E9",
    "ptCentDep3": "#117733",
    "ptCentDepFine": "#332288",
}

WINNER_MAP_COLORS = {
    "centINDcontrol_allRange": "#4E79A7",
    "centAsFeat_allRange": "#F28E2B",
    "centAsFeatMinOpt_pt5to40": "#8CD17D",
    "centDepBDTs_allRange": "#59A14F",
    "centDepFineBDTs_allRange": "#E15759",
    "ptBinCentAsFeat": "#EDC948",
    "ptCentDep3": "#76B7B2",
    "ptCentDepFine": "#1F4E79",
}

PRESENTATION_CMAP = "cividis"
FEATURE_CMAP = "YlGnBu"


def label(product: str) -> str:
    return LABELS.get(product, product.replace("_", " "))


def color(product: str) -> str:
    return COLORS.get(product, "#777777")


def short_label(product: str) -> str:
    return SHORT_LABELS.get(product, label(product))


def cell_label(product: str) -> str:
    return CELL_LABELS.get(product, short_label(product))


def sphenix_label(ax, y=0.98):
    ax.text(
        0.02,
        y,
        "sPHENIX",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=15,
        fontstyle="italic",
        fontweight="bold",
    )
    ax.text(
        0.315,
        y,
        "Internal",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=15,
    )
    ax.text(
        0.02,
        y - 0.075,
        r"Pythia overlay, $\sqrt{s_{NN}}=200$ GeV",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=12,
    )


def sphenix_fig_label(fig, x=0.12, y=0.96):
    fig.text(
        x,
        y,
        "sPHENIX",
        ha="left",
        va="top",
        fontsize=15,
        fontstyle="italic",
        fontweight="bold",
    )
    fig.text(
        x + 0.145,
        y,
        "Internal",
        ha="left",
        va="top",
        fontsize=15,
    )
    fig.text(
        x,
        y - 0.045,
        r"Pythia overlay, $\sqrt{s_{NN}}=200$ GeV",
        ha="left",
        va="top",
        fontsize=12,
    )


def polish(ax):
    ax.grid(True, axis="y", alpha=0.25, linewidth=0.8)
    ax.tick_params(labelsize=11)
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)


def heatmap_text_color(value: float, vmin: float, vmax: float, cmap_name: str) -> str:
    if not math.isfinite(value) or vmax <= vmin:
        return "black"
    rgba = matplotlib.colormaps[cmap_name]((value - vmin) / (vmax - vmin))
    r, g, b = rgba[:3]
    luminance = 0.2126 * r + 0.7152 * g + 0.0722 * b
    return "black" if luminance > 0.58 else "white"


def text_color_for_hex(hex_color: str) -> str:
    r, g, b = matplotlib.colors.to_rgb(hex_color)
    luminance = 0.2126 * r + 0.7152 * g + 0.0722 * b
    return "black" if luminance > 0.52 else "white"


def roc_auc(y_true: np.ndarray, score: np.ndarray) -> float:
    good = np.isfinite(score)
    y = np.asarray(y_true[good], dtype=np.int8)
    s = np.asarray(score[good], dtype=float)
    if y.size < 2 or len(np.unique(y)) < 2:
        return np.nan
    order = np.argsort(s)[::-1]
    y = y[order]
    pos = np.sum(y == 1)
    neg = np.sum(y == 0)
    if pos == 0 or neg == 0:
        return np.nan
    tpr = np.r_[0.0, np.cumsum(y == 1) / pos, 1.0]
    fpr = np.r_[0.0, np.cumsum(y == 0) / neg, 1.0]
    return float(np.trapezoid(tpr, fpr))


def parse_range_key(key: str) -> tuple[float, float]:
    lo, hi = key.split("_", 1)
    return float(lo), float(hi)


def range_label(key: str, suffix: str = "") -> str:
    lo, hi = parse_range_key(key)
    return f"{lo:g}-{hi:g}{suffix}"


def save(fig, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(path)


def plot_ranking(report: Path, out: Path):
    df = pd.read_csv(report / "validation_model_rankings.csv")
    products = [p for p in KEY_PRODUCTS + GRANULAR_PRODUCTS if p in set(df["product"])]
    df = df.set_index("product").loc[products].reset_index()

    fig, ax = plt.subplots(figsize=(9.4, 5.9))
    fig.subplots_adjust(top=0.82)
    fig.text(0.125, 0.965, "Expanded Au+Au tight-BDT model comparison", ha="left", va="top", fontsize=16)
    fig.text(0.78, 0.965, "sPHENIX", ha="left", va="top", fontsize=15, fontstyle="italic", fontweight="bold")
    fig.text(0.91, 0.965, "Internal", ha="left", va="top", fontsize=15)
    fig.text(0.985, 0.92, r"Pythia overlay, $\sqrt{s_{NN}}=200$ GeV", ha="right", va="top", fontsize=12)
    y = np.arange(len(df))
    h = 0.34
    ax.barh(y - h / 2, df["auc_inclusive"], height=h, color="#4C78A8", label="Inclusive AUC")
    ax.barh(y + h / 2, df["auc_pt_cent_mean"], height=h, color="#F58518", label=r"Mean $E_T\times$centrality AUC")
    ax.set_yticks(y)
    ax.set_yticklabels([label(p) for p in df["product"]])
    ax.invert_yaxis()
    ax.set_xlim(0.74, 0.90)
    ax.set_xlabel("ROC AUC", fontsize=13)
    ax.legend(loc="lower right", frameon=False, fontsize=11)
    for yy, val in zip(y, df["auc_inclusive"]):
        ax.text(val + 0.002, yy - h / 2, f"{val:.3f}", va="center", fontsize=9)
    polish(ax)
    save(fig, out / "money_model_ranking_auc.png")


def plot_centrality(report: Path, out: Path):
    diag = json.loads((report / "validation_deep_diagnostics.json").read_text())
    fig, ax = plt.subplots(figsize=(8.7, 5.8))
    for product in KEY_PRODUCTS:
        pdata = diag["products"].get(product)
        if not pdata:
            continue
        items = sorted(pdata["auc_by_centrality"].items(), key=lambda kv: parse_range_key(kv[0])[0])
        x = np.arange(len(items))
        y = [cell["auc"] for _, cell in items]
        ax.plot(x, y, marker="o", linewidth=2.4, color=color(product), label=label(product))
    ax.axvspan(-0.35, 0.35, color="#EEEEEE", zorder=-10)
    ax.text(0, 0.868, "central Au+Au focus", ha="center", va="bottom", fontsize=10, color="#555555")
    ax.set_xticks(np.arange(3))
    ax.set_xticklabels(["0-20%", "20-50%", "50-80%"])
    ax.set_ylim(0.86, 0.895)
    ax.set_ylabel("ROC AUC", fontsize=13)
    ax.set_title("Tight-BDT separation versus centrality", fontsize=16, pad=14)
    ax.legend(loc="lower right", frameon=False, fontsize=10)
    polish(ax)
    sphenix_label(ax, y=0.98)
    save(fig, out / "money_auc_by_centrality.png")


def plot_pt(report: Path, out: Path):
    diag = json.loads((report / "validation_deep_diagnostics.json").read_text())
    fig, ax = plt.subplots(figsize=(9.2, 5.8))
    fig.subplots_adjust(top=0.78)
    for product in KEY_PRODUCTS:
        pdata = diag["products"].get(product)
        if not pdata:
            continue
        items = sorted(pdata["auc_by_pt"].items(), key=lambda kv: parse_range_key(kv[0])[0])
        keys = [k for k, _ in items if k != "35_40"]
        y = [pdata["auc_by_pt"][k]["auc"] for k in keys]
        ax.plot(np.arange(len(keys)), y, marker="o", linewidth=2.4, color=color(product), label=label(product))
    ax.set_xticks(np.arange(5))
    ax.set_xticklabels(["6-10", "10-15", "15-20", "20-25", "25-35"])
    ax.set_ylim(0.60, 0.91)
    ax.set_xlabel(r"Cluster $E_T$ validation bin [GeV]", fontsize=13)
    ax.set_ylabel("ROC AUC", fontsize=13)
    ax.set_title(r"Tight-BDT separation versus cluster $E_T$", fontsize=16, pad=14)
    ax.text(
        0.98,
        0.08,
        r"Dedicated $35$-$40$ GeV pT-binned BDTs skipped;"
        "\n"
        r"tail covered by inclusive $5$-$40$ GeV models.",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=10,
        color="#444444",
    )
    ax.legend(loc="upper right", frameon=False, fontsize=9)
    polish(ax)
    sphenix_fig_label(fig, x=0.13, y=0.94)
    save(fig, out / "money_auc_by_cluster_et.png")


def plot_heatmaps(report: Path, out: Path):
    diag = json.loads((report / "validation_deep_diagnostics.json").read_text())
    products = ["centINDcontrol_allRange", "centAsFeat_allRange", "centDepBDTs_allRange", "centDepFineBDTs_allRange"]
    cent_keys = ["0_20", "20_50", "50_80"]
    pt_keys = ["6_10", "10_15", "15_20", "20_25", "25_35"]
    fig, axes = plt.subplots(2, 2, figsize=(10.2, 7.4), sharex=True, sharey=True)
    fig.subplots_adjust(top=0.80, bottom=0.14, left=0.10, right=0.88, hspace=0.42, wspace=0.20)
    vmin, vmax = 0.60, 0.90
    last = None
    for ax, product in zip(axes.ravel(), products):
        pdata = diag["products"][product]["auc_by_centrality_and_pt"]
        z = np.full((len(cent_keys), len(pt_keys)), np.nan)
        for ic, ckey in enumerate(cent_keys):
            for ip, pkey in enumerate(pt_keys):
                cell = pdata.get(ckey, {}).get(pkey)
                if cell:
                    z[ic, ip] = cell["auc"]
        last = ax.imshow(z, vmin=vmin, vmax=vmax, cmap=PRESENTATION_CMAP, aspect="auto")
        ax.set_title(label(product), fontsize=13)
        ax.set_xticks(np.arange(len(pt_keys)))
        ax.set_xticklabels([range_label(k) for k in pt_keys], rotation=35, ha="right")
        ax.set_yticks(np.arange(len(cent_keys)))
        ax.set_yticklabels([range_label(k, "%") for k in cent_keys])
        for ic in range(len(cent_keys)):
            for ip in range(len(pt_keys)):
                if math.isfinite(z[ic, ip]):
                    txt = heatmap_text_color(z[ic, ip], vmin, vmax, PRESENTATION_CMAP)
                    ax.text(ip, ic, f"{z[ic, ip]:.2f}", ha="center", va="center", color=txt, fontsize=9)
    fig.suptitle(r"AUC map in centrality and cluster $E_T$", fontsize=17, y=0.99)
    fig.text(0.10, 0.925, "sPHENIX", ha="left", va="top", fontsize=13, fontstyle="italic", fontweight="bold")
    fig.text(0.19, 0.925, "Internal", ha="left", va="top", fontsize=13)
    fig.text(0.10, 0.895, r"Pythia overlay, $\sqrt{s_{NN}}=200$ GeV", ha="left", va="top", fontsize=10)
    fig.text(0.5, 0.04, r"Cluster $E_T$ validation bin [GeV]", ha="center", fontsize=12)
    fig.text(0.04, 0.5, "Centrality", va="center", rotation=90, fontsize=12)
    cbar = fig.colorbar(last, ax=axes.ravel().tolist(), shrink=0.88, pad=0.02)
    cbar.set_label("ROC AUC")
    save(fig, out / "money_auc_heatmaps_cent_pt.png")


def plot_score_overlay(report: Path, out: Path):
    data = json.loads((report / "score_histograms.json").read_text())
    edges = np.asarray(data["bin_edges"], dtype=float)
    centers = 0.5 * (edges[:-1] + edges[1:])
    product = "centDepFineBDTs_allRange"
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.8), sharey=True)
    fig.subplots_adjust(top=0.74)
    for ax, ckey, title in zip(axes, ["0_20", "50_80"], ["0-20% central", "50-80% peripheral"]):
        cdata = data["products"][product]["by_centrality"][ckey]
        ax.step(centers, cdata["background"]["density"], where="mid", color="#4C78A8", linewidth=2.2, label="background")
        ax.step(centers, cdata["signal"]["density"], where="mid", color="#E45756", linewidth=2.5, label="signal")
        ax.set_yscale("log")
        ax.set_xlim(0, 1)
        ax.set_ylim(0.008, 6)
        ax.set_xlabel("BDT score", fontsize=12)
        ax.set_title(title, fontsize=13)
        polish(ax)
    axes[0].set_ylabel("Unit-normalized candidates", fontsize=12)
    axes[1].legend(loc="upper left", frameon=False, fontsize=11)
    sphenix_fig_label(fig, x=0.12, y=0.93)
    fig.suptitle("Signal/background score separation for best inclusive candidate", fontsize=16, y=0.99)
    save(fig, out / "money_score_overlay_best_central_peripheral.png")


def plot_thresholds(report: Path, out: Path):
    diag = json.loads((report / "validation_deep_diagnostics.json").read_text())
    products = KEY_PRODUCTS
    rows = []
    for product in products:
        thresh = diag["products"][product]["thresholds_inclusive"]["by_signal_efficiency"]
        match = min(thresh, key=lambda item: abs(item["target_signal_efficiency"] - 0.8))
        rows.append((product, match["background_rejection"], match["threshold"]))
    fig, ax = plt.subplots(figsize=(8.6, 5.2))
    x = np.arange(len(rows))
    vals = [r[1] for r in rows]
    ax.bar(x, vals, color=[color(r[0]) for r in rows])
    ax.set_xticks(x)
    ax.set_xticklabels([label(r[0]) for r in rows], rotation=25, ha="right")
    ax.set_ylim(0.80, 0.86)
    ax.set_ylabel("Background rejection at 80% signal efficiency", fontsize=12)
    ax.set_title("Working-point tradeoff for top model families", fontsize=15, pad=12)
    for xx, val, thr in zip(x, vals, [r[2] for r in rows]):
        ax.text(xx, val + 0.002, f"{val:.3f}\nscore>{thr:.2f}", ha="center", va="bottom", fontsize=9)
    polish(ax)
    sphenix_label(ax, y=0.98)
    save(fig, out / "money_background_rejection_at_80sig.png")


def plot_feature_separation(report: Path, out: Path):
    df = pd.read_csv(report / "money_tables" / "money_feature_separation.csv")
    df = df[df["feature"] != "reco_eiso"]
    df = df.head(10).iloc[::-1]
    fig, ax = plt.subplots(figsize=(8.6, 5.8))
    y = np.arange(len(df))
    ax.barh(y, df["separation"], color="#4C78A8")
    ax.set_yticks(y)
    ax.set_yticklabels(df["feature"])
    ax.set_xlabel("Signal/background separation of input distribution", fontsize=12)
    ax.set_title("Most separating input variables before the BDT", fontsize=15, pad=12)
    polish(ax)
    ax.text(
        0.84,
        0.20,
        "sPHENIX",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=15,
        fontstyle="italic",
        fontweight="bold",
    )
    ax.text(
        0.98,
        0.20,
        "Internal",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=15,
    )
    ax.text(
        0.98,
        0.125,
        r"Pythia overlay, $\sqrt{s_{NN}}=200$ GeV",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=12,
    )
    save(fig, out / "money_feature_separation_inputs.png")


def plot_auc_gain_maps(report: Path, out: Path):
    diag = json.loads((report / "validation_deep_diagnostics.json").read_text())
    baseline = diag["products"]["centINDcontrol_allRange"]["auc_by_centrality_and_pt"]
    products = ["centAsFeat_allRange", "centDepBDTs_allRange", "centDepFineBDTs_allRange"]
    cent_keys = ["0_20", "20_50", "50_80"]
    pt_keys = ["6_10", "10_15", "15_20", "20_25", "25_35"]
    fig, axes = plt.subplots(1, 3, figsize=(12.0, 4.2), sharex=True, sharey=True)
    fig.subplots_adjust(top=0.76, bottom=0.20, left=0.08, right=0.90, wspace=0.22)
    vmax = 0.025
    last = None
    for ax, product in zip(axes, products):
        pdata = diag["products"][product]["auc_by_centrality_and_pt"]
        z = np.full((len(cent_keys), len(pt_keys)), np.nan)
        for ic, ckey in enumerate(cent_keys):
            for ip, pkey in enumerate(pt_keys):
                cell = pdata.get(ckey, {}).get(pkey)
                base = baseline.get(ckey, {}).get(pkey)
                if cell and base:
                    z[ic, ip] = cell["auc"] - base["auc"]
        last = ax.imshow(z, vmin=-vmax, vmax=vmax, cmap="RdBu_r", aspect="auto")
        ax.set_title(label(product), fontsize=12)
        ax.set_xticks(np.arange(len(pt_keys)))
        ax.set_xticklabels([range_label(k) for k in pt_keys], rotation=35, ha="right")
        ax.set_yticks(np.arange(len(cent_keys)))
        ax.set_yticklabels([range_label(k, "%") for k in cent_keys])
        for ic in range(len(cent_keys)):
            for ip in range(len(pt_keys)):
                if math.isfinite(z[ic, ip]):
                    txt = heatmap_text_color(z[ic, ip], -vmax, vmax, "RdBu_r")
                    ax.text(ip, ic, f"{z[ic, ip]:+.2f}", ha="center", va="center", color=txt, fontsize=8)
    fig.suptitle(r"AUC gain over no-centrality baseline", fontsize=15, y=0.98)
    fig.text(0.08, 0.90, "sPHENIX", ha="left", va="top", fontsize=13, fontstyle="italic", fontweight="bold")
    fig.text(0.17, 0.90, "Internal", ha="left", va="top", fontsize=13)
    fig.text(0.08, 0.855, r"Pythia overlay, $\sqrt{s_{NN}}=200$ GeV", ha="left", va="top", fontsize=10)
    fig.text(0.5, 0.055, r"Cluster $E_T$ validation bin [GeV]", ha="center", fontsize=11)
    fig.text(0.02, 0.5, "Centrality", va="center", rotation=90, fontsize=11)
    cbar = fig.colorbar(last, ax=axes.ravel().tolist(), shrink=0.82, pad=0.02)
    cbar.set_label(r"$\Delta$ROC AUC")
    save(fig, out / "money_auc_gain_over_no_centrality_baseline.png")


def plot_working_point_curves(report: Path, out: Path):
    th = pd.read_csv(report / "validation_threshold_table.csv")
    products = KEY_PRODUCTS
    fig, ax = plt.subplots(figsize=(8.8, 5.6))
    fig.subplots_adjust(top=0.80)
    for product in products:
        rows = th[
            (th["product"] == product)
            & (th["scope"] == "inclusive")
            & (th["target_type"] == "signal_efficiency")
        ].sort_values("signal_efficiency")
        ax.plot(
            rows["signal_efficiency"],
            rows["background_rejection"],
            marker="o",
            linewidth=2.3,
            color=color(product),
            label=label(product),
        )
    ax.set_xlim(0.48, 0.965)
    ax.set_ylim(0.62, 0.93)
    ax.set_xlabel("Signal efficiency retained by BDT threshold", fontsize=12)
    ax.set_ylabel("Background rejection", fontsize=12)
    ax.set_title("Efficiency-rejection working-point scan", fontsize=15, pad=12)
    ax.legend(loc="lower left", frameon=False, fontsize=9)
    polish(ax)
    sphenix_fig_label(fig, x=0.13, y=0.95)
    save(fig, out / "money_working_point_scan.png")


def plot_signal_eff_at_fixed_rejection(report: Path, out: Path):
    th = pd.read_csv(report / "validation_threshold_table.csv")
    products = KEY_PRODUCTS
    scopes_pt = ["pt_6_10", "pt_10_15", "pt_15_20", "pt_20_25", "pt_25_35"]
    scopes_cent = ["cent_0_20", "cent_20_50", "cent_50_80"]
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.7))
    for product in products:
        rows = th[
            (th["product"] == product)
            & (th["target_type"] == "background_fake_rate")
            & (th["target"] == 0.1)
        ]
        by_scope = rows.set_index("scope")
        ypt = [by_scope.loc[s, "signal_efficiency"] if s in by_scope.index else np.nan for s in scopes_pt]
        ycent = [by_scope.loc[s, "signal_efficiency"] if s in by_scope.index else np.nan for s in scopes_cent]
        axes[0].plot(np.arange(len(scopes_pt)), ypt, marker="o", linewidth=2.0, color=color(product), label=label(product))
        axes[1].plot(np.arange(len(scopes_cent)), ycent, marker="o", linewidth=2.0, color=color(product), label=label(product))
    axes[0].set_xticks(np.arange(len(scopes_pt)))
    axes[0].set_xticklabels(["6-10", "10-15", "15-20", "20-25", "25-35"], rotation=25, ha="right")
    axes[1].set_xticks(np.arange(len(scopes_cent)))
    axes[1].set_xticklabels(["0-20%", "20-50%", "50-80%"])
    for ax in axes:
        ax.set_ylim(0.0, 0.78)
        ax.set_ylabel("Signal efficiency at 90% background rejection", fontsize=11)
        polish(ax)
    axes[0].set_xlabel(r"Cluster $E_T$ validation bin [GeV]", fontsize=11)
    axes[1].set_xlabel("Centrality", fontsize=11)
    axes[1].legend(loc="lower right", frameon=False, fontsize=8)
    sphenix_label(axes[0], y=0.98)
    fig.suptitle("Candidate-level BDT efficiency at a fixed rejection point", fontsize=15, y=1.02)
    save(fig, out / "money_signal_eff_at_90rej.png")


def plot_feature_separation_by_centrality(report: Path, out: Path):
    df = pd.read_csv(report / "validation_feature_summary.csv")
    df = df[df["feature"] != "reco_eiso"]
    piv = {}
    for scope in ["inclusive", "cent_0_20", "cent_20_50", "cent_50_80"]:
        sub = df[df["scope"] == scope]
        sig = sub[sub["class"] == "signal"].set_index("feature")
        bkg = sub[sub["class"] == "background"].set_index("feature")
        common = sig.index.intersection(bkg.index)
        sep = (sig.loc[common, "mean"] - bkg.loc[common, "mean"]).abs()
        denom = np.sqrt(0.5 * (sig.loc[common, "std"] ** 2 + bkg.loc[common, "std"] ** 2)).replace(0, np.nan)
        piv[scope] = (sep / denom).replace([np.inf, -np.inf], np.nan)
    top = piv["inclusive"].sort_values(ascending=False).head(10).index.tolist()
    mat = np.array([[piv[scope].get(feat, np.nan) for scope in ["inclusive", "cent_0_20", "cent_20_50", "cent_50_80"]] for feat in top])
    fig, ax = plt.subplots(figsize=(8.8, 6.1))
    fig.subplots_adjust(top=0.78)
    vmin = 0.0
    vmax = np.nanmax(mat) * 1.05
    im = ax.imshow(mat, cmap=FEATURE_CMAP, aspect="auto", vmin=vmin, vmax=vmax)
    ax.set_yticks(np.arange(len(top)))
    ax.set_yticklabels(top)
    ax.set_xticks(np.arange(4))
    ax.set_xticklabels(["inclusive", "0-20%", "20-50%", "50-80%"])
    for iy in range(mat.shape[0]):
        for ix in range(mat.shape[1]):
            if math.isfinite(mat[iy, ix]):
                txt = heatmap_text_color(mat[iy, ix], vmin, vmax, FEATURE_CMAP)
                ax.text(ix, iy, f"{mat[iy, ix]:.2f}", ha="center", va="center", color=txt, fontsize=8)
    ax.set_title("Input-variable separation by centrality", fontsize=15, pad=20)
    cbar = fig.colorbar(im, ax=ax, pad=0.02)
    cbar.set_label("Signal/background separation")
    fig.text(0.125, 0.955, "sPHENIX", ha="left", va="top", fontsize=14, fontstyle="italic", fontweight="bold")
    fig.text(0.245, 0.955, "Internal", ha="left", va="top", fontsize=14)
    fig.text(0.125, 0.91, r"Pythia overlay, $\sqrt{s_{NN}}=200$ GeV", ha="left", va="top", fontsize=11)
    save(fig, out / "money_feature_separation_by_centrality.png")


def load_score_cache_columns(report: Path, products: list[str]):
    files = sorted((report / "score_caches").glob("score_cache_*.npz"))
    cols = {"is_signal": [], "cluster_Et": [], "centrality": [], "reco_eiso": []}
    for product in products:
        cols[f"score_{product}"] = []
    for path in files:
        with np.load(path) as z:
            for key in ["is_signal", "cluster_Et", "centrality", "reco_eiso"]:
                cols[key].append(np.asarray(z[key]))
            for product in products:
                cols[f"score_{product}"].append(np.asarray(z[f"score_{product}"]))
    return {key: np.concatenate(vals) for key, vals in cols.items()}


def plot_score_response_maps(report: Path, out: Path):
    product = "centDepFineBDTs_allRange"
    tab = load_score_cache_columns(report, [product])
    score = tab[f"score_{product}"]
    sig = tab["is_signal"] == 1
    bkg = ~sig
    et_edges = np.array([5, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 35], dtype=float)
    eiso_edges = np.array([-10, -5, -2, 0, 2, 5, 10, 20, 40], dtype=float)
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.9), sharey=True)
    for ax, mask, title, col in [
        (axes[0], sig, "truth-matched signal", "#E45756"),
        (axes[1], bkg, "embedded-jet background", "#4C78A8"),
    ]:
        means_et = []
        centers_et = []
        for lo, hi in zip(et_edges[:-1], et_edges[1:]):
            m = mask & (tab["cluster_Et"] >= lo) & (tab["cluster_Et"] < hi) & np.isfinite(score)
            if m.sum() > 0:
                means_et.append(np.mean(score[m]))
                centers_et.append(0.5 * (lo + hi))
        ax.plot(centers_et, means_et, marker="o", color=col, linewidth=2.2, label=r"vs $E_T$")
        ax.set_xlabel(r"Cluster $E_T$ [GeV]", fontsize=11)
        ax.set_title(title, fontsize=12)
        ax.set_ylim(0.0, 0.9)
        polish(ax)
    axes[0].set_ylabel("Mean BDT score", fontsize=12)
    sphenix_label(axes[0], y=0.98)
    fig.suptitle("BDT response versus candidate energy", fontsize=15, y=1.02)
    save(fig, out / "money_score_response_vs_et.png")

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.9), sharey=True)
    for ax, mask, title, col in [
        (axes[0], sig, "truth-matched signal", "#E45756"),
        (axes[1], bkg, "embedded-jet background", "#4C78A8"),
    ]:
        means = []
        centers = []
        for lo, hi in zip(eiso_edges[:-1], eiso_edges[1:]):
            m = mask & (tab["reco_eiso"] >= lo) & (tab["reco_eiso"] < hi) & np.isfinite(score)
            if m.sum() > 20:
                means.append(np.mean(score[m]))
                centers.append(0.5 * (lo + hi))
        ax.plot(centers, means, marker="o", color=col, linewidth=2.2)
        ax.set_xlabel(r"Reco isolation energy [GeV]", fontsize=11)
        ax.set_title(title, fontsize=12)
        ax.set_ylim(0.0, 0.9)
        polish(ax)
    axes[0].set_ylabel("Mean BDT score", fontsize=12)
    sphenix_label(axes[0], y=0.98)
    fig.suptitle("BDT response versus isolation environment", fontsize=15, y=1.02)
    save(fig, out / "money_score_response_vs_eiso.png")


def plot_isolation_threshold_scan(report: Path, out: Path):
    tab = load_score_cache_columns(report, ["centDepFineBDTs_allRange"])
    eiso = tab["reco_eiso"]
    is_signal = tab["is_signal"] == 1
    cent = tab["centrality"]
    thresholds = np.linspace(-8.0, 16.0, 80)
    cent_bins = [(0, 20, "0-20%"), (20, 50, "20-50%"), (50, 80, "50-80%")]
    cent_cols = ["#0072B2", "#009E73", "#D55E00"]

    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.9), sharex=True)
    for (lo, hi, lab), col in zip(cent_bins, cent_cols):
        cmask = (cent >= lo) & (cent < hi)
        sig = cmask & is_signal & np.isfinite(eiso)
        bkg = cmask & (~is_signal) & np.isfinite(eiso)
        if sig.sum() == 0 or bkg.sum() == 0:
            continue
        sig_eff = [np.mean(eiso[sig] < thr) for thr in thresholds]
        bkg_acc = [np.mean(eiso[bkg] < thr) for thr in thresholds]
        axes[0].plot(thresholds, sig_eff, color=col, linewidth=2.4, label=lab)
        axes[1].plot(thresholds, bkg_acc, color=col, linewidth=2.4, label=lab)

    for ax in axes:
        ax.set_xlabel(r"Reco isolation requirement, $E_{\mathrm{iso}}<$ cut [GeV]", fontsize=11)
        ax.set_ylim(0, 1.02)
        polish(ax)
    axes[0].set_ylabel("Truth-matched photon retention", fontsize=11)
    axes[1].set_ylabel("Embedded-jet background acceptance", fontsize=11)
    axes[1].legend(loc="lower right", frameon=False, fontsize=10)
    sphenix_label(axes[0], y=0.94)
    fig.suptitle("Isolation threshold scan on validation candidates", fontsize=15, y=1.02)
    save(fig, out / "money_isolation_threshold_scan_by_centrality.png")


def plot_auc_by_isolation_region(report: Path, out: Path):
    products = KEY_PRODUCTS
    tab = load_score_cache_columns(report, products)
    eiso = tab["reco_eiso"]
    y = tab["is_signal"]
    regions = [
        (-np.inf, 0.0, r"$E_{\mathrm{iso}}<0$"),
        (0.0, 2.0, "0-2"),
        (2.0, 5.0, "2-5"),
        (5.0, 10.0, "5-10"),
        (10.0, np.inf, r"$>10$"),
    ]
    x = np.arange(len(regions))

    fig, ax = plt.subplots(figsize=(9.2, 5.6))
    for product in products:
        scores = tab[f"score_{product}"]
        vals = []
        for lo, hi, _ in regions:
            mask = (eiso >= lo) & (eiso < hi) & np.isfinite(eiso)
            vals.append(roc_auc(y[mask], scores[mask]) if mask.sum() > 200 else np.nan)
        ax.plot(x, vals, marker="o", linewidth=2.3, color=color(product), label=label(product))
    ax.set_xticks(x)
    ax.set_xticklabels([lab for _, _, lab in regions])
    ax.set_ylim(0.80, 0.93)
    ax.set_xlabel(r"Reco isolation region [GeV]", fontsize=12)
    ax.set_ylabel("ROC AUC inside isolation region", fontsize=12)
    ax.set_title("BDT separation stability versus isolation region", fontsize=15, pad=12)
    ax.text(
        0.02,
        0.08,
        r"$E_{\mathrm{iso}}$ is a validation split, not a BDT input",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=10,
        color="#444444",
    )
    ax.legend(loc="lower right", frameon=False, fontsize=8.5)
    polish(ax)
    sphenix_label(ax, y=0.98)
    save(fig, out / "money_bdt_auc_by_isolation_region.png")


def plot_score_overlay_by_isolation_region(report: Path, out: Path):
    product = "centDepFineBDTs_allRange"
    tab = load_score_cache_columns(report, [product])
    score = tab[f"score_{product}"]
    sig = tab["is_signal"] == 1
    eiso = tab["reco_eiso"]
    cent = tab["centrality"]
    central = (cent >= 0) & (cent < 20)
    regions = [
        (sig & central & (eiso < 2.0), "signal, isolated", "#D55E00", "-"),
        (sig & central & (eiso >= 5.0), "signal, non-isolated", "#D55E00", "--"),
        ((~sig) & central & (eiso < 2.0), "background, isolated", "#0072B2", "-"),
        ((~sig) & central & (eiso >= 5.0), "background, non-isolated", "#0072B2", "--"),
    ]
    bins = np.linspace(0.0, 1.0, 51)
    centers = 0.5 * (bins[:-1] + bins[1:])

    fig, ax = plt.subplots(figsize=(8.8, 5.6))
    fig.subplots_adjust(top=0.78)
    for mask, lab, col, ls in regions:
        mask = mask & np.isfinite(score)
        if mask.sum() < 100:
            continue
        hist, _ = np.histogram(score[mask], bins=bins, density=True)
        ax.step(centers, hist, where="mid", color=col, linestyle=ls, linewidth=2.3, label=lab)
    ax.set_yscale("log")
    ax.set_xlim(0, 1)
    ax.set_ylim(0.01, 8)
    ax.set_xlabel("BDT score", fontsize=12)
    ax.set_ylabel("Unit-normalized candidates", fontsize=12)
    ax.set_title("Score separation split by isolation in 0-20% Au+Au", fontsize=15, pad=12)
    ax.legend(loc="lower left", frameon=False, fontsize=9)
    polish(ax)
    sphenix_fig_label(fig, x=0.13, y=0.95)
    save(fig, out / "money_score_overlay_by_isolation_region.png")


def plot_performance_vs_complexity(report: Path, out: Path):
    df = pd.read_csv(report / "validation_model_rankings.csv")
    products = [p for p in ALL_PRODUCTS if p in set(df["product"])]
    df = df.set_index("product").loc[products].reset_index()
    fig, ax = plt.subplots(figsize=(8.8, 5.8))
    offset = {
        "centINDcontrol_allRange": (1.08, -0.0020),
        "centAsFeat_allRange": (1.08, 0.0003),
        "centAsFeatMinOpt_pt5to40": (1.08, 0.0026),
    }
    for _, row in df.iterrows():
        product = row["product"]
        x = MODEL_COMPLEXITY.get(product, np.nan)
        ax.scatter(
            x,
            row["auc_pt_cent_mean"],
            s=120,
            color=color(product),
            edgecolor="black",
            linewidth=0.8,
            zorder=3,
        )
        xmul, dy = offset.get(product, (1.06, 0.0))
        ax.text(x * xmul, row["auc_pt_cent_mean"] + dy, label(product), va="center", fontsize=9)
    ax.set_xscale("log")
    ax.set_xlim(0.75, 120)
    ax.set_ylim(0.76, 0.835)
    ax.set_xlabel("Number of separately trained BDTs in model family", fontsize=12)
    ax.set_ylabel(r"Mean AUC across $E_T\times$centrality cells", fontsize=12)
    ax.set_title("Performance gained versus model complexity", fontsize=15, pad=12)
    ax.text(
        0.03,
        0.04,
        "Useful complexity should raise the difficult-bin average,\nnot only the inclusive score.",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=10,
        color="#444444",
    )
    polish(ax)
    sphenix_label(ax, y=0.98)
    save(fig, out / "money_performance_vs_complexity.png")


def plot_auc_stability_bands(report: Path, out: Path):
    auc = pd.read_csv(report / "validation_auc_table.csv")
    products = [p for p in ALL_PRODUCTS if p in set(auc["product"])]
    rows = []
    for product in products:
        vals = auc.loc[auc["product"] == product, "auc"].dropna().to_numpy()
        if vals.size:
            rows.append((product, np.nanmin(vals), np.nanpercentile(vals, 25), np.nanmedian(vals), np.nanpercentile(vals, 75), np.nanmax(vals)))
    fig, ax = plt.subplots(figsize=(9.4, 6.0))
    fig.subplots_adjust(top=0.78)
    y = np.arange(len(rows))
    for yy, (product, vmin, q25, med, q75, vmax) in zip(y, rows):
        ax.hlines(yy, vmin, vmax, color=color(product), linewidth=2.0, alpha=0.35)
        ax.hlines(yy, q25, q75, color=color(product), linewidth=8.0, alpha=0.85)
        ax.plot(med, yy, marker="o", color="white", markeredgecolor="black", markersize=7, zorder=4)
        ax.text(vmax + 0.003, yy, f"worst {vmin:.2f}", va="center", fontsize=8, color="#444444")
    ax.set_yticks(y)
    ax.set_yticklabels([label(r[0]) for r in rows])
    ax.invert_yaxis()
    ax.set_xlim(0.58, 0.94)
    ax.set_xlabel(r"AUC across centrality and $E_T$ validation cells", fontsize=12)
    ax.set_title("Model robustness across the full validation phase space", fontsize=15, pad=12)
    ax.text(0.02, 0.05, "Thick band = middle 50%; dot = median; thin line = full range", transform=ax.transAxes, fontsize=9)
    polish(ax)
    sphenix_fig_label(fig, x=0.19, y=0.94)
    save(fig, out / "money_auc_stability_bands.png")


def plot_phase_space_winner_map(report: Path, out: Path):
    diag = json.loads((report / "validation_deep_diagnostics.json").read_text())
    products = [p for p in ALL_PRODUCTS if p in diag["products"]]
    cent_keys = ["0_20", "20_50", "50_80"]
    pt_keys = ["6_10", "10_15", "15_20", "20_25", "25_35"]
    family_colors = {p: i for i, p in enumerate(products)}
    z = np.full((len(cent_keys), len(pt_keys)), np.nan)
    winners: list[list[str | None]] = [[None for _ in pt_keys] for _ in cent_keys]
    best_auc = np.full((len(cent_keys), len(pt_keys)), np.nan)
    for ic, ckey in enumerate(cent_keys):
        for ip, pkey in enumerate(pt_keys):
            best_product = None
            best_val = -np.inf
            for product in products:
                cell = diag["products"][product]["auc_by_centrality_and_pt"].get(ckey, {}).get(pkey)
                if cell and math.isfinite(cell["auc"]) and cell["auc"] > best_val:
                    best_product = product
                    best_val = cell["auc"]
            if best_product:
                z[ic, ip] = family_colors[best_product]
                winners[ic][ip] = best_product
                best_auc[ic, ip] = best_val
    winner_colors = {p: WINNER_MAP_COLORS.get(p, color(p)) for p in products}
    cmap = matplotlib.colors.ListedColormap([winner_colors[p] for p in products])
    fig, ax = plt.subplots(figsize=(9.4, 5.25))
    fig.subplots_adjust(top=0.73, bottom=0.23, left=0.13, right=0.98)
    ax.imshow(z, cmap=cmap, vmin=-0.5, vmax=len(products) - 0.5, aspect="auto")
    ax.set_xticks(np.arange(len(pt_keys)))
    ax.set_xticklabels([range_label(k) for k in pt_keys], rotation=25, ha="right")
    ax.set_yticks(np.arange(len(cent_keys)))
    ax.set_yticklabels([range_label(k, "%") for k in cent_keys])
    for ic in range(len(cent_keys)):
        for ip in range(len(pt_keys)):
            product = winners[ic][ip]
            if product:
                txt = text_color_for_hex(winner_colors[product])
                ax.text(
                    ip,
                    ic,
                    f"{cell_label(product)}\nAUC={best_auc[ic, ip]:.2f}",
                    ha="center",
                    va="center",
                    fontsize=8.5,
                    color=txt,
                    linespacing=0.95,
                    fontweight="bold",
                )
    ax.set_xlabel(r"Cluster $E_T$ validation bin [GeV]", fontsize=12)
    ax.set_ylabel("Centrality", fontsize=12)
    fig.suptitle("Winning model family by BDT separation quality", fontsize=15, y=0.985)
    fig.text(
        0.13,
        0.905,
        "Each cell shows the model family with the highest AUC; AUC=1 is perfect separation, AUC=0.5 is random.",
        ha="left",
        va="top",
        fontsize=10,
        color="#333333",
    )
    fig.text(0.13, 0.855, "sPHENIX", ha="left", va="top", fontsize=13, fontstyle="italic", fontweight="bold")
    fig.text(0.26, 0.855, "Internal", ha="left", va="top", fontsize=13)
    fig.text(0.13, 0.805, r"Pythia overlay, $\sqrt{s_{NN}}=200$ GeV", ha="left", va="top", fontsize=10)
    save(fig, out / "money_phase_space_winner_map.png")


def plot_phase_space_winner_margin_map(report: Path, out: Path):
    diag = json.loads((report / "validation_deep_diagnostics.json").read_text())
    products = [p for p in ALL_PRODUCTS if p in diag["products"]]
    cent_keys = ["0_20", "20_50", "50_80"]
    pt_keys = ["6_10", "10_15", "15_20", "20_25", "25_35"]
    margin = np.full((len(cent_keys), len(pt_keys)), np.nan)
    winners: list[list[str | None]] = [[None for _ in pt_keys] for _ in cent_keys]
    for ic, ckey in enumerate(cent_keys):
        for ip, pkey in enumerate(pt_keys):
            vals: list[tuple[float, str]] = []
            for product in products:
                cell = diag["products"][product]["auc_by_centrality_and_pt"].get(ckey, {}).get(pkey)
                if cell and math.isfinite(cell["auc"]):
                    vals.append((float(cell["auc"]), product))
            vals.sort(reverse=True)
            if len(vals) >= 2 and vals[1][0] > 0:
                best_auc, best_product = vals[0]
                second_auc, _ = vals[1]
                winners[ic][ip] = best_product
                margin[ic, ip] = 100.0 * (best_auc - second_auc) / second_auc

    fig, ax = plt.subplots(figsize=(9.4, 5.25))
    fig.subplots_adjust(top=0.73, bottom=0.23, left=0.13, right=0.98)
    vmax = max(1.0, float(np.nanmax(margin)) if np.isfinite(margin).any() else 1.0)
    im = ax.imshow(margin, cmap="YlGnBu", vmin=0.0, vmax=vmax, aspect="auto")
    ax.set_xticks(np.arange(len(pt_keys)))
    ax.set_xticklabels([range_label(k) for k in pt_keys], rotation=25, ha="right")
    ax.set_yticks(np.arange(len(cent_keys)))
    ax.set_yticklabels([range_label(k, "%") for k in cent_keys])
    for ic in range(len(cent_keys)):
        for ip in range(len(pt_keys)):
            val = margin[ic, ip]
            product = winners[ic][ip]
            if math.isfinite(val) and product:
                txt = "white" if val > 0.55 * vmax else "black"
                ax.text(
                    ip,
                    ic,
                    f"+{val:.1f}%\n{cell_label(product)}",
                    ha="center",
                    va="center",
                    fontsize=8.5,
                    color=txt,
                    linespacing=0.95,
                    fontweight="bold",
                )
    cb = fig.colorbar(im, ax=ax, fraction=0.040, pad=0.018)
    cb.set_label("AUC gain over next-best model [%]", fontsize=10)
    ax.set_xlabel(r"Cluster $E_T$ validation bin [GeV]", fontsize=12)
    ax.set_ylabel("Centrality", fontsize=12)
    fig.suptitle("How decisive is each winning model?", fontsize=15, y=0.985)
    fig.text(
        0.13,
        0.905,
        "Percent difference = 100 x (best AUC - next-best AUC) / next-best AUC in the same bin.",
        ha="left",
        va="top",
        fontsize=10,
        color="#333333",
    )
    fig.text(0.13, 0.855, "sPHENIX", ha="left", va="top", fontsize=13, fontstyle="italic", fontweight="bold")
    fig.text(0.26, 0.855, "Internal", ha="left", va="top", fontsize=13)
    fig.text(0.13, 0.805, r"Pythia overlay, $\sqrt{s_{NN}}=200$ GeV", ha="left", va="top", fontsize=10)
    save(fig, out / "money_phase_space_winner_margin_pct.png")


def plot_working_points_by_centrality(report: Path, out: Path):
    th = pd.read_csv(report / "validation_threshold_table.csv")
    products = KEY_PRODUCTS
    scopes = ["cent_0_20", "cent_20_50", "cent_50_80"]
    targets = [0.7, 0.8, 0.9]
    fig, axes = plt.subplots(1, 3, figsize=(13.0, 4.8), sharey=True)
    fig.subplots_adjust(top=0.74, wspace=0.20)
    for ax, target in zip(axes, targets):
        for product in products:
            rows = th[
                (th["product"] == product)
                & (th["target_type"] == "signal_efficiency")
                & (np.isclose(th["target"], target))
            ].set_index("scope")
            vals = [rows.loc[s, "background_rejection"] if s in rows.index else np.nan for s in scopes]
            ax.plot(np.arange(len(scopes)), vals, marker="o", linewidth=2.1, color=color(product), label=label(product))
        ax.set_title(f"{int(target * 100)}% signal retained", fontsize=12)
        ax.set_xticks(np.arange(len(scopes)))
        ax.set_xticklabels(["0-20%", "20-50%", "50-80%"])
        ax.set_ylim(0.55, 0.93)
        ax.set_xlabel("Centrality", fontsize=11)
        polish(ax)
    axes[0].set_ylabel("Background rejection at fixed signal efficiency", fontsize=11)
    axes[-1].legend(loc="lower right", frameon=False, fontsize=8)
    fig.text(0.12, 0.91, "sPHENIX", ha="left", va="top", fontsize=13, fontstyle="italic", fontweight="bold")
    fig.text(0.24, 0.91, "Internal", ha="left", va="top", fontsize=13)
    fig.text(0.12, 0.86, r"Pythia overlay, $\sqrt{s_{NN}}=200$ GeV", ha="left", va="top", fontsize=10)
    fig.suptitle("Working-point stability versus centrality", fontsize=15, y=0.99)
    save(fig, out / "money_working_points_by_centrality.png")


def plot_score_correlation_heatmap(report: Path, out: Path):
    diag = json.loads((report / "validation_deep_diagnostics.json").read_text())
    products = KEY_PRODUCTS
    features = [
        "e32_over_e35",
        "cluster_et1",
        "cluster_weta_cogx",
        "cluster_wphi_cogx",
        "cluster_Et",
        "e11_over_e33",
        "cluster_et2",
        "cluster_et3",
        "cluster_et4",
        "centrality",
        "reco_eiso",
    ]
    mat = np.full((len(features), len(products)), np.nan)
    for ip, product in enumerate(products):
        corr = diag["products"][product]["score_correlations"]
        for ifeat, feat in enumerate(features):
            mat[ifeat, ip] = corr.get(feat, np.nan)
    fig, ax = plt.subplots(figsize=(9.8, 6.6))
    fig.subplots_adjust(top=0.80, left=0.20)
    im = ax.imshow(mat, vmin=-0.8, vmax=0.8, cmap="RdBu_r", aspect="auto")
    ax.set_yticks(np.arange(len(features)))
    ylabels = [f"{f}*" if f == "reco_eiso" else f for f in features]
    ax.set_yticklabels(ylabels)
    ax.set_xticks(np.arange(len(products)))
    ax.set_xticklabels([label(p) for p in products], rotation=25, ha="right")
    for iy in range(mat.shape[0]):
        for ix in range(mat.shape[1]):
            if math.isfinite(mat[iy, ix]):
                txt = heatmap_text_color(mat[iy, ix], -0.8, 0.8, "RdBu_r")
                ax.text(ix, iy, f"{mat[iy, ix]:+.2f}", ha="center", va="center", fontsize=8, color=txt)
    ax.set_title("BDT score correlations with candidate observables", fontsize=15, pad=16)
    ax.text(0.0, -0.22, "* isolation is diagnostic only, not a BDT input", transform=ax.transAxes, fontsize=9)
    cbar = fig.colorbar(im, ax=ax, pad=0.02)
    cbar.set_label("Pearson correlation with BDT score")
    fig.text(0.20, 0.96, "sPHENIX", ha="left", va="top", fontsize=13, fontstyle="italic", fontweight="bold")
    fig.text(0.32, 0.96, "Internal", ha="left", va="top", fontsize=13)
    fig.text(0.20, 0.915, r"Pythia overlay, $\sqrt{s_{NN}}=200$ GeV", ha="left", va="top", fontsize=10)
    save(fig, out / "money_score_correlations_observables.png")


def plot_isolation_operating_map(report: Path, out: Path):
    tab = load_score_cache_columns(report, ["centDepFineBDTs_allRange"])
    eiso = tab["reco_eiso"]
    sig = tab["is_signal"] == 1
    cent = tab["centrality"]
    et = tab["cluster_Et"]
    cent_bins = [(0, 20), (20, 50), (50, 80)]
    pt_bins = [(6, 10), (10, 15), (15, 20), (20, 25), (25, 35)]
    cut = 2.0
    mats = [np.full((len(cent_bins), len(pt_bins)), np.nan), np.full((len(cent_bins), len(pt_bins)), np.nan)]
    for ic, (clo, chi) in enumerate(cent_bins):
        for ip, (plo, phi) in enumerate(pt_bins):
            base = (cent >= clo) & (cent < chi) & (et >= plo) & (et < phi) & np.isfinite(eiso)
            sm = base & sig
            bm = base & (~sig)
            if sm.sum() > 0:
                mats[0][ic, ip] = np.mean(eiso[sm] < cut)
            if bm.sum() > 0:
                mats[1][ic, ip] = np.mean(eiso[bm] < cut)
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.8), sharex=True, sharey=True)
    fig.subplots_adjust(top=0.74, bottom=0.22, left=0.10, right=0.88, wspace=0.18)
    titles = [r"Truth-matched photon retention", r"Embedded-jet background acceptance"]
    for ax, mat, title in zip(axes, mats, titles):
        im = ax.imshow(mat, vmin=0, vmax=1, cmap=PRESENTATION_CMAP, aspect="auto")
        ax.set_title(title, fontsize=12)
        ax.set_xticks(np.arange(len(pt_bins)))
        ax.set_xticklabels([f"{lo}-{hi}" for lo, hi in pt_bins], rotation=25, ha="right")
        ax.set_yticks(np.arange(len(cent_bins)))
        ax.set_yticklabels([f"{lo}-{hi}%" for lo, hi in cent_bins])
        for ic in range(len(cent_bins)):
            for ip in range(len(pt_bins)):
                if math.isfinite(mat[ic, ip]):
                    txt = heatmap_text_color(mat[ic, ip], 0, 1, PRESENTATION_CMAP)
                    ax.text(ip, ic, f"{mat[ic, ip]:.2f}", ha="center", va="center", color=txt, fontsize=9)
    fig.suptitle(rf"Isolation operating map for $E_{{\mathrm{{iso}}}}<{cut:g}$ GeV", fontsize=15, y=0.99)
    fig.text(0.5, 0.055, r"Cluster $E_T$ validation bin [GeV]", ha="center", fontsize=11)
    fig.text(0.04, 0.5, "Centrality", va="center", rotation=90, fontsize=11)
    cbar = fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.85, pad=0.02)
    cbar.set_label("Fraction passing isolation requirement")
    fig.text(0.10, 0.91, "sPHENIX", ha="left", va="top", fontsize=13, fontstyle="italic", fontweight="bold")
    fig.text(0.23, 0.91, "Internal", ha="left", va="top", fontsize=13)
    fig.text(0.10, 0.86, r"Pythia overlay, $\sqrt{s_{NN}}=200$ GeV", ha="left", va="top", fontsize=10)
    save(fig, out / "money_isolation_operating_map_eiso_lt2.png")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("report_dir", type=Path)
    parser.add_argument("--out-dir", type=Path, default=None)
    args = parser.parse_args()
    report = args.report_dir
    out = args.out_dir or report / "money_plots"
    plot_ranking(report, out)
    plot_centrality(report, out)
    plot_pt(report, out)
    plot_heatmaps(report, out)
    plot_score_overlay(report, out)
    plot_thresholds(report, out)
    plot_feature_separation(report, out)
    plot_auc_gain_maps(report, out)
    plot_working_point_curves(report, out)
    plot_signal_eff_at_fixed_rejection(report, out)
    plot_feature_separation_by_centrality(report, out)
    plot_score_response_maps(report, out)
    plot_isolation_threshold_scan(report, out)
    plot_auc_by_isolation_region(report, out)
    plot_score_overlay_by_isolation_region(report, out)
    plot_performance_vs_complexity(report, out)
    plot_auc_stability_bands(report, out)
    plot_phase_space_winner_map(report, out)
    plot_phase_space_winner_margin_map(report, out)
    plot_working_points_by_centrality(report, out)
    plot_score_correlation_heatmap(report, out)
    plot_isolation_operating_map(report, out)
    print(f"[OK] wrote plots to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
