#!/usr/bin/env python3
"""Reco-efficiency overlay for target-80 RecoilJets reruns."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


CENTS = [
    ("0_20", "0-20%"),
    ("20_50", "20-50%"),
    ("50_80", "50-80%"),
]


def load_rows(path: Path, *, source: str, label: str) -> list[dict]:
    rows: list[dict] = []
    with path.open() as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            cent = row.get("cent") or row.get("centrality")
            if cent not in {c for c, _ in CENTS}:
                continue
            x = row.get("pt_mid")
            y = row.get("value")
            ey = row.get("error")
            if x is None or y is None or ey is None:
                continue
            rows.append(
                {
                    "source": source,
                    "label": label,
                    "centrality": cent,
                    "pt_mid": float(x),
                    "value": float(y),
                    "error": float(ey),
                    "numerator": row.get("numerator", row.get("numerator_or_truth", "")),
                    "denominator": row.get("denominator", ""),
                    "input_csv": str(path),
                }
            )
    return rows


def setup_axes(ax, ylabel: str | None, xlabel: str | None) -> None:
    ax.set_xlim(14.5, 35.5)
    ax.set_ylim(0.0, 0.62)
    ax.grid(True, color="#dfdfdf", linewidth=0.8, alpha=0.75)
    ax.tick_params(direction="in", top=True, right=True, labelsize=12, length=6)
    ax.minorticks_on()
    if ylabel:
        ax.set_ylabel("Reco efficiency", fontsize=18)
    if xlabel:
        ax.set_xlabel(r"truth photon $p_T$ [GeV]", fontsize=16)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reference-csv", required=True, type=Path)
    ap.add_argument("--variant-csv", "--basev3e-csv", action="append", required=True, type=Path)
    ap.add_argument("--variant-source", action="append", default=[])
    ap.add_argument("--variant-label", action="append", default=[])
    ap.add_argument("--stem", default="target80_vs_reference_boxcuts_reco_efficiency_pt15to35")
    ap.add_argument("--out-dir", required=True, type=Path)
    args = ap.parse_args()
    if args.variant_label and len(args.variant_label) != len(args.variant_csv):
        raise ValueError("--variant-label count must match --variant-csv count")
    if args.variant_source and len(args.variant_source) != len(args.variant_csv):
        raise ValueError("--variant-source count must match --variant-csv count")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    rows.extend(load_rows(args.reference_csv, source="reference_box_cuts", label="Reference box cuts"))
    for idx, path in enumerate(args.variant_csv):
        label = args.variant_label[idx] if args.variant_label else f"Target-80 ML model {idx + 1}"
        source = args.variant_source[idx] if args.variant_source else f"target80_variant_{idx + 1}"
        rows.extend(load_rows(path, source=source, label=label))

    fig, axes = plt.subplots(1, 3, figsize=(17.2, 5.2), sharey=True)
    def style_for(label: str) -> dict:
        if label == "Reference box cuts":
            return dict(color="#252525", marker="s", linestyle="none", lw=0.0, ms=7.4)
        if "Global" in label or "global" in label:
            return dict(color="#CC79A7", marker="D", linestyle="none", lw=0.0, ms=7.8)
        return dict(color="#0072B2", marker="o", linestyle="none", lw=0.0, ms=8.2)

    labels = []
    for row in rows:
        if row["label"] not in labels:
            labels.append(row["label"])

    for ic, (cent, title) in enumerate(CENTS):
        ax = axes[ic]
        for label in labels:
            style = style_for(label)
            pts = sorted([r for r in rows if r["centrality"] == cent and r["label"] == label], key=lambda r: r["pt_mid"])
            if not pts:
                continue
            x = np.array([p["pt_mid"] for p in pts], dtype=float)
            y = np.array([p["value"] for p in pts], dtype=float)
            ey = np.array([p["error"] for p in pts], dtype=float)
            ax.errorbar(
                x,
                y,
                yerr=ey,
                label=label,
                capsize=2.2,
                elinewidth=1.2,
                markerfacecolor=style["color"],
                markeredgecolor=style["color"],
                **style,
            )
        ax.set_title(title, fontsize=20, fontweight="bold", pad=10)
        setup_axes(ax, "Reco efficiency" if ic == 0 else None, r"truth photon $p_T$ [GeV]" if ic == 1 else None)

    axes[2].legend(loc="lower right", fontsize=13, frameon=True, framealpha=0.92, borderpad=0.55).get_frame().set_linewidth(0)
    fig.text(0.023, 0.972, "sPHENIX", fontsize=18, fontweight="bold", fontstyle="italic", va="top")
    fig.text(0.111, 0.972, "Internal", fontsize=18, va="top")
    fig.text(0.023, 0.902, r"Embedded Photon12+20 signal, $15 < p_{T} < 35$ GeV", fontsize=14, color="#333333")
    fig.text(0.023, 0.855, "Reco efficiency = 1 - truth photons missed after reconstruction and ID selection", fontsize=14, color="#333333")
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.81), w_pad=1.2)

    png = args.out_dir / f"{args.stem}_1x3.png"
    out_csv = args.out_dir / f"{args.stem}_points.csv"
    fig.savefig(png, dpi=190)
    plt.close(fig)

    with out_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["source", "label", "centrality", "pt_mid", "value", "error", "numerator", "denominator", "input_csv"],
        )
        writer.writeheader()
        writer.writerows(rows)

    print(png)
    print(out_csv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
