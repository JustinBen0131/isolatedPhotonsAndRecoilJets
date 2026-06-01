#!/usr/bin/env python3
"""Summarize raw-isolation routed-BDT split gain into slide-sized groups."""

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
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


LABELS = {
    "reco_eiso_r30": r"raw R=0.3 $E_T^{iso}$",
    "reco_eiso_r40": r"raw R=0.4 $E_T^{iso}$",
    "e22_over_e37": r"$E_{22}/E_{37}$",
    "e22_over_e53": r"$E_{22}/E_{53}$",
    "e22_over_e35": r"$E_{22}/E_{35}$",
    "e22_over_e33": r"$E_{22}/E_{33}$",
}

COLORS = {
    "isolation": "#CC79A7",
    "core energy ratios": "#009E73",
    "other energy sharing": "#6EE7B7",
    "shower widths": "#0072B2",
    "kinematics/context": "#6B7280",
    "other small inputs": "#CBD5E1",
}

SHOWER_WIDTHS = {
    "cluster_weta_cogx",
    "cluster_wphi_cogx",
    "cluster_weta33_cogx",
    "cluster_wphi33_cogx",
    "cluster_weta35_cogx",
    "cluster_wphi53_cogx",
    "cluster_w32",
    "cluster_w52",
    "cluster_w72",
    "cluster_weta_over_wphi",
    "cluster_weta33_over_wphi33",
}

KINEMATICS = {"cluster_Et", "cluster_Eta", "vertexz", "centrality"}
E22_RATIOS = {"e22_over_e53", "e22_over_e35", "e22_over_e33"}


def read_rows(path: Path) -> list[dict]:
    with path.open() as handle:
        rows = list(csv.DictReader(handle))
    for row in rows:
        row["gain_value"] = float(row.get("split_gain") or row.get("route_gain") or 0.0)
        row["split_count"] = int(row.get("split_count") or row.get("route_split_count") or 0)
    return rows


def add_group(groups: dict[str, dict], key: str, label: str, family: str, frac: float, count: int) -> None:
    item = groups.setdefault(key, {"label": label, "family": family, "fraction": 0.0, "split_count": 0})
    item["fraction"] += frac
    item["split_count"] += count


def summarize(rows: list[dict], *, keep_r30_only: bool = False) -> list[dict]:
    groups: dict[str, dict] = {}
    total_gain = sum(float(row["gain_value"]) for row in rows)
    for row in rows:
        feature = row["feature"]
        frac = float(row["gain_value"]) / total_gain if total_gain > 0 else 0.0
        count = row["split_count"]
        if feature in {"reco_eiso_r30", "reco_eiso_r40", "e22_over_e37"}:
            add_group(groups, feature, LABELS[feature], "isolation" if feature.startswith("reco_eiso") else "core energy ratios", frac, count)
        elif feature in E22_RATIOS:
            add_group(groups, "other_e22_ratios", r"other $E_{22}$ ratios", "core energy ratios", frac, count)
        elif feature in SHOWER_WIDTHS:
            add_group(groups, "shower_widths", "all shower-width inputs", "shower widths", frac, count)
        elif feature in KINEMATICS:
            add_group(groups, "kinematics", r"$E_T$, $\eta$, $z_{vtx}$, centrality", "kinematics/context", frac, count)
        elif row["feature_family"] == "energy sharing":
            add_group(groups, "other_energy_sharing", "other tower-energy sharing", "other energy sharing", frac, count)
        else:
            add_group(groups, "other_small_inputs", "other small inputs", "other small inputs", frac, count)

    preferred_order = [
        "reco_eiso_r30",
        "reco_eiso_r40",
        "e22_over_e37",
        "other_e22_ratios",
        "shower_widths",
        "other_energy_sharing",
        "kinematics",
        "other_small_inputs",
    ]
    out = [groups[key] for key in preferred_order if key in groups and groups[key]["fraction"] > 0]
    if keep_r30_only:
        return out
    return out


def write_summary_csv(rows: list[dict], path: Path) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["label", "family", "fraction", "percent", "split_count"])
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "label": row["label"],
                    "family": row["family"],
                    "fraction": f"{row['fraction']:.9f}",
                    "percent": f"{100.0 * row['fraction']:.4f}",
                    "split_count": row["split_count"],
                }
            )


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.linewidth": 1.1,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )


def draw_single(rows: list[dict], *, title: str, subtitle: str, out: Path) -> None:
    setup_style()
    fig, ax = plt.subplots(figsize=(14.8, 8.4), dpi=180)
    fig.patch.set_facecolor("white")
    y = list(range(len(rows)))
    vals = [100.0 * row["fraction"] for row in rows]
    colors = [COLORS[row["family"]] for row in rows]
    labels = [row["label"] for row in rows]

    ax.barh(y, vals, color=colors, height=0.58)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=14)
    ax.invert_yaxis()
    ax.set_xlabel("Fraction of total split gain [%]", fontsize=15)
    ax.tick_params(axis="x", labelsize=12.5)
    ax.grid(axis="x", color="#CBD5E1", linewidth=0.8, alpha=0.8)
    ax.set_axisbelow(True)
    xmax = max(vals) * 1.18
    ax.set_xlim(0, xmax)
    for yy, val, row in zip(y, vals, rows):
        ax.text(val + xmax * 0.012, yy, f"{val:.1f}% ({row['split_count']})", ha="left", va="center", fontsize=13, fontweight="bold" if val >= 5 else "normal")

    fig.text(0.065, 0.965, title, ha="left", va="top", fontsize=23, fontweight="bold")
    fig.text(0.065, 0.920, subtitle, ha="left", va="top", fontsize=13.2, color="#374151")
    fig.text(0.065, 0.877, "sPHENIX", ha="left", va="top", fontsize=17, fontstyle="italic", fontweight="bold")
    fig.text(0.139, 0.877, " Internal", ha="left", va="top", fontsize=17)
    fig.text(0.065, 0.842, r"Photon12+20 signal vs Jet12+20+30 background, $15 < E_T < 35$ GeV", ha="left", va="top", fontsize=12.5, color="#111827")
    fig.text(0.065, 0.060, "Long tail of sub-percent variables is grouped by physics family; raw cone-isolation inputs are shown explicitly.", ha="left", va="bottom", fontsize=11.0, color="#6B7280")
    fig.subplots_adjust(left=0.285, right=0.960, top=0.790, bottom=0.135)
    fig.savefig(out)
    plt.close(fig)


def draw_comparison(r30_rows: list[dict], r40_rows: list[dict], both_rows: list[dict], out: Path) -> None:
    setup_style()
    fig, axes = plt.subplots(1, 3, figsize=(23.5, 8.4), dpi=180, sharex=True)
    fig.patch.set_facecolor("white")
    for ax, rows, title in zip(
        axes,
        [r30_rows, r40_rows, both_rows],
        [r"raw R=0.3 $E_T^{iso}$", r"raw R=0.4 $E_T^{iso}$", r"raw R=0.3 + R=0.4 $E_T^{iso}$"],
    ):
        y = list(range(len(rows)))
        vals = [100.0 * row["fraction"] for row in rows]
        colors = [COLORS[row["family"]] for row in rows]
        labels = [row["label"] for row in rows]
        ax.barh(y, vals, color=colors, height=0.58)
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=11.5)
        ax.tick_params(axis="y", pad=2)
        ax.invert_yaxis()
        ax.set_title(title, fontsize=15.5, fontweight="bold", pad=12)
        ax.grid(axis="x", color="#CBD5E1", linewidth=0.8, alpha=0.8)
        ax.set_axisbelow(True)
        ax.set_xlim(0, 72)
        for yy, val, row in zip(y, vals, rows):
            ax.text(val + 0.8, yy, f"{val:.1f}%", ha="left", va="center", fontsize=11.5, fontweight="bold" if val >= 5 else "normal")
    for ax in axes:
        ax.set_xlabel("Fraction of total split gain [%]", fontsize=12.5)
    fig.text(0.055, 0.965, "Split-gain summary: raw isolation routed BDTs", ha="left", va="top", fontsize=22, fontweight="bold")
    fig.text(0.055, 0.922, r"8 $E_T$ x 7 centrality routing; 56 route-specific XGBoost BDTs per column", ha="left", va="top", fontsize=13, color="#374151")
    fig.text(0.790, 0.950, "sPHENIX", ha="left", va="top", fontsize=17, fontstyle="italic", fontweight="bold")
    fig.text(0.865, 0.950, " Internal", ha="left", va="top", fontsize=17)
    fig.text(0.055, 0.060, "Interpretation: compare R=0.3-only and R=0.4-only directly; the both-cone model shows whether the two radii carry complementary split gain.", ha="left", va="bottom", fontsize=11.0, color="#6B7280")
    fig.subplots_adjust(left=0.115, right=0.985, top=0.820, bottom=0.135, wspace=0.62)
    fig.savefig(out)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--r30-csv", type=Path, required=True)
    parser.add_argument("--r40-csv", type=Path, required=True)
    parser.add_argument("--r30r40-csv", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    r30 = summarize(read_rows(args.r30_csv))
    r40 = summarize(read_rows(args.r40_csv))
    r30r40 = summarize(read_rows(args.r30r40_csv))
    write_summary_csv(r30, args.outdir / "raw_eisoR30_ptCent7_split_gain_summary.csv")
    write_summary_csv(r40, args.outdir / "raw_eisoR40_ptCent7_split_gain_summary.csv")
    write_summary_csv(r30r40, args.outdir / "raw_eisoR30R40_ptCent7_split_gain_summary.csv")
    draw_single(
        r30,
        title=r"Split-gain summary: raw R=0.3 $E_T^{iso}$ BDT",
        subtitle=r"Grouped view of 8 $E_T$ x 7 centrality routed split gain; labels give gain fraction and split count",
        out=args.outdir / "raw_eisoR30_ptCent7_split_gain_summary.png",
    )
    draw_single(
        r40,
        title=r"Split-gain summary: raw R=0.4 $E_T^{iso}$ BDT",
        subtitle=r"Grouped view of 8 $E_T$ x 7 centrality routed split gain; labels give gain fraction and split count",
        out=args.outdir / "raw_eisoR40_ptCent7_split_gain_summary.png",
    )
    draw_single(
        r30r40,
        title=r"Split-gain summary: raw R=0.3 + R=0.4 $E_T^{iso}$ BDT",
        subtitle=r"Grouped view of 8 $E_T$ x 7 centrality routed split gain; labels give gain fraction and split count",
        out=args.outdir / "raw_eisoR30R40_ptCent7_split_gain_summary.png",
    )
    draw_comparison(r30, r40, r30r40, args.outdir / "raw_eiso_split_gain_summary_comparison_1x3.png")
    print(args.outdir / "raw_eisoR30_ptCent7_split_gain_summary.png")
    print(args.outdir / "raw_eisoR40_ptCent7_split_gain_summary.png")
    print(args.outdir / "raw_eisoR30R40_ptCent7_split_gain_summary.png")
    print(args.outdir / "raw_eiso_split_gain_summary_comparison_1x3.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
