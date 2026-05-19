#!/usr/bin/env python3
"""Make a slide-sized split-gain summary for the isolation-input routed BDT."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


FAMILY_COLORS = {
    "isolation inputs": "#C65F9B",
    "energy sharing": "#009E73",
    "shower widths": "#0072B2",
    "kinematic/context": "#6B7280",
    "other inputs": "#BDBDBD",
}

SLIDE_LABELS = {
    "reco_eiso_clip30": r"R=0.4 $E_T^{iso}$, clipped",
    "reco_eiso_signed_log1p": r"R=0.4 signed log$(1+|E_T^{iso}|)$",
    "reco_eiso_over_cluster_Et": r"R=0.4 $E_T^{iso}/E_T^{cluster}$",
    "reco_eiso_r30": r"R=0.3 $E_T^{iso}$",
    "reco_eiso_r40": r"R=0.4 $E_T^{iso}$",
    "e22_over_e37": r"$E_{22}/E_{37}$",
    "e22_over_e53": r"$E_{22}/E_{53}$",
    "cluster_wphi53_cogx": r"$w_{\phi}^{5x3}$",
    "e11_over_e33": r"$E_{11}/E_{33}$",
    "cluster_et3": "cluster et3",
    "e22_over_e35": r"$E_{22}/E_{35}$",
    "cluster_Eta": r"cluster $\eta$",
    "cluster_weta35_cogx": r"$w_{\eta}^{3x5}$",
    "cluster_et2": "cluster et2",
    "cluster_weta_cogx": r"$w_{\eta}$",
    "cluster_wphi33_cogx": r"$w_{\phi}^{3x3}$",
    "cluster_w72": r"$w_{72}$",
}


def read_rows(path: Path) -> list[dict]:
    with path.open() as f:
        rows = list(csv.DictReader(f))
    for row in rows:
        row["split_gain_fraction"] = float(row["split_gain_fraction"])
        row["split_count"] = int(row["split_count"])
    rows.sort(key=lambda item: item["split_gain_fraction"], reverse=True)
    return rows


def summarize(rows: list[dict], top_n: int) -> list[dict]:
    top = [dict(row) for row in rows[:top_n]]
    other = rows[top_n:]
    if other:
        top.append(
            {
                "feature": f"all_other_{len(other)}",
                "feature_label": f"All other {len(other)} inputs",
                "feature_family": "other inputs",
                "split_gain_fraction": sum(float(row["split_gain_fraction"]) for row in other),
                "split_count": sum(int(row["split_count"]) for row in other),
            }
        )
    return top


def draw(rows: list[dict], out: Path, *, top_n: int) -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.linewidth": 1.15,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )
    fig, ax = plt.subplots(figsize=(15.6, 8.75), dpi=220)
    fig.patch.set_facecolor("white")
    y = list(range(len(rows)))
    values = [100.0 * float(row["split_gain_fraction"]) for row in rows]
    labels = [SLIDE_LABELS.get(str(row["feature"]), str(row["feature_label"])) for row in rows]
    colors = [FAMILY_COLORS.get(str(row["feature_family"]), FAMILY_COLORS["other inputs"]) for row in rows]

    ax.barh(y, values, color=colors, height=0.70)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=15.2)
    ax.invert_yaxis()
    ax.set_xlim(0.0, 31.0)
    ax.set_xlabel("Fraction of total split gain [%]", fontsize=17)
    ax.tick_params(axis="x", labelsize=14.2)
    ax.grid(axis="x", color="#D1D5DB", linewidth=1.0, alpha=0.75)
    ax.set_axisbelow(True)
    for yy, val, row in zip(y, values, rows):
        ax.text(
            val + 0.22,
            yy,
            f"{val:.1f}% ({int(row['split_count']):,})",
            ha="left",
            va="center",
            fontsize=14.4,
            fontweight="bold" if val >= 2.0 else "normal",
        )

    legend_items = [
        Line2D([0], [0], color=FAMILY_COLORS[name], lw=9, label=name)
        for name in ("isolation inputs", "energy sharing", "shower widths", "kinematic/context", "other inputs")
    ]
    ax.legend(
        handles=legend_items,
        loc="lower right",
        frameon=True,
        facecolor="white",
        edgecolor="#E5E7EB",
        fontsize=13.1,
        title="Feature family",
        title_fontsize=13.8,
    )

    fig.text(0.052, 0.962, "Split-gain usage for the isolation-input binned BDT", ha="left", va="top", fontsize=25, fontweight="bold")
    fig.text(
        0.052,
        0.919,
        rf"8 $E_T$ x 7 centrality XGBoost BDTs; gain summed over 56 models; top {top_n} inputs plus remainder",
        ha="left",
        va="top",
        fontsize=15.2,
        color="#374151",
    )
    fig.text(
        0.052,
        0.885,
        r"Isolation inputs are transforms of reconstructed $E_T^{iso}$ from the R=0.4 cone training tree, not binary isolation pass/fail flags",
        ha="left",
        va="top",
        fontsize=13.4,
        color="#111827",
    )
    fig.text(0.835, 0.958, "sPHENIX", ha="left", va="top", fontsize=19, fontstyle="italic", fontweight="bold")
    fig.text(0.918, 0.958, " Internal", ha="left", va="top", fontsize=19)

    fig.tight_layout(rect=[0.045, 0.055, 0.985, 0.825])
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out)
    plt.close(fig)


def write_summary_csv(rows: list[dict], out: Path) -> None:
    fieldnames = ["feature", "feature_label", "feature_family", "split_gain_fraction", "split_count"]
    with out.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-csv", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--top-n", type=int, default=14)
    parser.add_argument("--tag", default="binned_bdt_iso_ptcent7_split_gain_slide31_summary")
    args = parser.parse_args()

    rows = summarize(read_rows(args.input_csv), args.top_n)
    png = args.outdir / f"{args.tag}.png"
    csv_path = args.outdir / f"{args.tag}.csv"
    draw(rows, png, top_n=args.top_n)
    write_summary_csv(rows, csv_path)
    print(png)
    print(csv_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
