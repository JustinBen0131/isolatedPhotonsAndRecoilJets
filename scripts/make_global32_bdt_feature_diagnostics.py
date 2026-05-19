#!/usr/bin/env python3
from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(
    "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/auauMLDiagnosticRuns/"
    "global_etcent_inclusive3_sixpack_20260516_135439"
)
OUTDIR = ROOT / "slideReady" / "global32_feature_diagnostics"
META_PATH = OUTDIR / "auau_tight_bdt_globalEtCent1535_bdt_noIso_tmva.metadata.json"
XGB_PATH = OUTDIR / "auau_tight_bdt_globalEtCent1535_bdt_noIso_tmva.xgb.json"
FEATURE_SUMMARY = ROOT / "validation" / "bdt_finished_only_20260516_180817" / "validation_feature_summary.csv"


LABELS = {
    "cluster_Et": r"cluster $E_T$",
    "cluster_Eta": r"cluster $\eta$",
    "vertexz": r"$z_{vtx}$",
    "centrality": "centrality",
    "cluster_weta_cogx": r"$w_{\eta}$",
    "cluster_wphi_cogx": r"$w_{\phi}$",
    "cluster_weta33_cogx": r"$w_{\eta}^{3x3}$",
    "cluster_wphi33_cogx": r"$w_{\phi}^{3x3}$",
    "cluster_weta35_cogx": r"$w_{\eta}^{3x5}$",
    "cluster_wphi53_cogx": r"$w_{\phi}^{5x3}$",
    "cluster_weta_over_wphi": r"$w_{\eta}/w_{\phi}$",
    "cluster_weta33_over_wphi33": r"$w_{\eta}^{3x3}/w_{\phi}^{3x3}$",
    "cluster_et1": "cluster et1",
    "cluster_et2": "cluster et2",
    "cluster_et3": "cluster et3",
    "cluster_et4": "cluster et4",
    "cluster_w32": r"$w_{32}$",
    "cluster_w52": r"$w_{52}$",
    "cluster_w72": r"$w_{72}$",
    "e11_over_e33": r"$E_{11}/E_{33}$",
    "e32_over_e35": r"$E_{32}/E_{35}$",
    "e11_over_e22": r"$E_{11}/E_{22}$",
    "e11_over_e13": r"$E_{11}/E_{13}$",
    "e11_over_e15": r"$E_{11}/E_{15}$",
    "e11_over_e17": r"$E_{11}/E_{17}$",
    "e11_over_e31": r"$E_{11}/E_{31}$",
    "e11_over_e51": r"$E_{11}/E_{51}$",
    "e11_over_e71": r"$E_{11}/E_{71}$",
    "e22_over_e33": r"$E_{22}/E_{33}$",
    "e22_over_e35": r"$E_{22}/E_{35}$",
    "e22_over_e37": r"$E_{22}/E_{37}$",
    "e22_over_e53": r"$E_{22}/E_{53}$",
}


def read_json_with_shell_preamble(path: Path) -> dict:
    text = path.read_text()
    start = text.find("{")
    if start < 0:
        raise ValueError(f"No JSON object found in {path}")
    return json.loads(text[start:])


def extract_split_gain(xgb_path: Path, features: list[str]) -> pd.DataFrame:
    data = read_json_with_shell_preamble(xgb_path)
    trees = data["learner"]["gradient_booster"]["model"]["trees"]
    gain = defaultdict(float)
    count = defaultdict(int)
    for tree in trees:
        for left, right, idx, loss in zip(
            tree["left_children"],
            tree["right_children"],
            tree["split_indices"],
            tree["loss_changes"],
        ):
            if left < 0 and right < 0:
                continue
            if idx < len(features):
                feature = features[int(idx)]
                gain[feature] += max(float(loss), 0.0)
                count[feature] += 1
    total = sum(gain.values())
    rows = []
    for feature in features:
        rows.append(
            {
                "feature": feature,
                "split_count": count.get(feature, 0),
                "split_gain": gain.get(feature, 0.0),
                "split_gain_fraction": gain.get(feature, 0.0) / total if total else 0.0,
            }
        )
    return pd.DataFrame(rows).sort_values("split_gain_fraction", ascending=False)


def build_separation(features: list[str]) -> pd.DataFrame:
    raw = pd.read_csv(FEATURE_SUMMARY)
    raw = raw[raw["scope"].eq("inclusive") & raw["feature"].isin(features)]
    sig = raw[raw["class"].eq("signal")].set_index("feature")
    bkg = raw[raw["class"].eq("background")].set_index("feature")
    rows = []
    for feature in features:
        s = sig.loc[feature]
        b = bkg.loc[feature]
        denom = np.sqrt(max(1.0e-12, 0.5 * (float(s["std"]) ** 2 + float(b["std"]) ** 2)))
        sep = (float(s["mean"]) - float(b["mean"])) / denom
        rows.append(
            {
                "feature": feature,
                "signal_entries": int(s["entries"]),
                "background_entries": int(b["entries"]),
                "signal_mean": float(s["mean"]),
                "background_mean": float(b["mean"]),
                "signal_std": float(s["std"]),
                "background_std": float(b["std"]),
                "standardized_mean_difference": sep,
                "abs_standardized_mean_difference": abs(sep),
            }
        )
    return pd.DataFrame(rows).sort_values("abs_standardized_mean_difference", ascending=False)


def style_axis(ax):
    ax.tick_params(axis="x", direction="in", top=True, length=6, width=1.2, labelsize=14)
    ax.tick_params(axis="y", direction="in", right=True, length=6, width=1.2, labelsize=13)
    for spine in ax.spines.values():
        spine.set_linewidth(1.25)
        spine.set_color("0.12")
    ax.grid(axis="x", color="0.87", linewidth=1.0)
    ax.set_axisbelow(True)


def add_header(fig, title: str, subtitle: str):
    fig.text(0.145, 0.920, title, ha="left", va="bottom", fontsize=22, fontweight="bold")
    fig.text(0.145, 0.885, subtitle, ha="left", va="bottom", fontsize=15.5, color="0.30")
    fig.text(
        0.145,
        0.853,
        r"Embedded Photon12+20 signal vs embedded inclusive Jet12+20+30 background, 15 < cluster $E_T$ < 35 GeV",
        ha="left",
        va="bottom",
        fontsize=12.8,
        color="0.22",
    )
    fig.text(
        0.865,
        0.942,
        "sPHENIX",
        ha="right",
        va="bottom",
        fontsize=22,
        fontstyle="italic",
        fontweight="bold",
    )
    fig.text(0.965, 0.942, "Internal", ha="right", va="bottom", fontsize=22)


def plot_split_gain(df: pd.DataFrame, out: Path) -> None:
    plot_df = df.sort_values("split_gain_fraction", ascending=True).copy()
    plot_df["fraction_percent"] = 100.0 * plot_df["split_gain_fraction"]
    fig, ax = plt.subplots(figsize=(14.4, 17.0), dpi=180)
    fig.subplots_adjust(left=0.22, right=0.975, top=0.805, bottom=0.075)
    y = np.arange(len(plot_df))
    colors = plt.cm.Blues(np.linspace(0.35, 0.9, len(plot_df)))
    bars = ax.barh(y, plot_df["fraction_percent"], color=colors, edgecolor="none", height=0.56)
    ax.set_yticks(y)
    ax.set_yticklabels([LABELS.get(v, v) for v in plot_df["feature"]])
    ax.tick_params(axis="y", labelsize=17)
    ax.set_xlabel("Fraction of total split gain [%]", fontsize=22, labelpad=12)
    xmax = max(10.0, float(plot_df["fraction_percent"].max()) * 1.18)
    ax.set_xlim(0, xmax)
    style_axis(ax)
    ax.tick_params(axis="x", labelsize=17)
    ax.tick_params(axis="y", labelsize=17)
    xmax = ax.get_xlim()[1]
    for bar, (_, row) in zip(bars, plot_df.iterrows()):
        value = float(row["fraction_percent"])
        count = int(row["split_count"])
        x = min(value + 0.35, xmax - 4.8)
        ha = "left"
        if value > xmax - 5.8:
            x = xmax - 0.4
            ha = "right"
        ax.text(
            x,
            bar.get_y() + bar.get_height() / 2,
            f"{value:.1f}% ({count})",
            ha=ha,
            va="center",
            fontsize=16.0,
            fontweight="bold",
            color="0.05",
            clip_on=True,
        )
    add_header(
        fig,
        "XGBoost split-gain usage: 32-feature global BDT",
        "gain summed over all tree splits; labels give gain fraction and split count",
    )
    fig.savefig(out, facecolor="white")
    plt.close(fig)


def plot_separation(df: pd.DataFrame, out: Path) -> None:
    plot_df = df.sort_values("standardized_mean_difference", ascending=True).copy()
    fig, ax = plt.subplots(figsize=(14.4, 10.2), dpi=180)
    fig.subplots_adjust(left=0.18, right=0.97, top=0.78, bottom=0.11)
    y = np.arange(len(plot_df))
    vals = plot_df["standardized_mean_difference"].to_numpy()
    colors = np.where(vals >= 0, "#1f77b4", "#d95f02")
    bars = ax.barh(y, vals, color=colors, edgecolor="none", height=0.72)
    ax.axvline(0, color="0.12", linewidth=1.5)
    ax.set_yticks(y)
    ax.set_yticklabels([LABELS.get(v, v) for v in plot_df["feature"]])
    ax.set_xlabel("Standardized signal-background mean difference", fontsize=18, labelpad=8)
    lim = max(1.9, float(np.nanmax(np.abs(vals))) * 1.18)
    ax.set_xlim(-lim, lim)
    style_axis(ax)
    for bar, val in zip(bars, vals):
        if val >= 0:
            x = val + 0.035 * lim
            ha = "left"
        else:
            x = val - 0.035 * lim
            ha = "right"
        ax.text(
            x,
            bar.get_y() + bar.get_height() / 2,
            f"{val:+.2f}",
            ha=ha,
            va="center",
            fontsize=12.5,
            fontweight="bold",
            color="0.05",
            clip_on=True,
        )
    add_header(
        fig,
        "Input-variable separation: 32-feature global BDT",
        "blue = signal mean higher; orange = background mean higher; longer bars separate more strongly",
    )
    fig.savefig(out, facecolor="white")
    plt.close(fig)


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    meta = read_json_with_shell_preamble(META_PATH)
    features = list(meta["features"])
    if len(features) != 32:
        raise SystemExit(f"Expected 32 features, found {len(features)}")
    split = extract_split_gain(XGB_PATH, features)
    sep = build_separation(features)
    split.to_csv(OUTDIR / "global32_bdt_split_gain.csv", index=False)
    sep.to_csv(OUTDIR / "global32_bdt_feature_separation.csv", index=False)
    plot_split_gain(split, OUTDIR / "global32_bdt_split_gain_usage_all_features.png")
    plot_separation(sep, OUTDIR / "global32_bdt_input_separation_all_features.png")
    print(OUTDIR / "global32_bdt_split_gain_usage_all_features.png")
    print(OUTDIR / "global32_bdt_input_separation_all_features.png")


if __name__ == "__main__":
    main()
