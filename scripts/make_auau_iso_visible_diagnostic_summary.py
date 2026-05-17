#!/usr/bin/env python3
"""Summarize isolation-visible AuAu photon-ID diagnostic validation outputs."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


CENT_BINS = [("0_20", "0-20%"), ("20_50", "20-50%"), ("50_80", "50-80%")]
WARNING = "uses isolation-derived inputs; diagnostic ceiling test only; not ABCD-safe photon ID"


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bdt-report", type=Path, default=None)
    ap.add_argument("--mlp-report", type=Path, default=None)
    ap.add_argument("--iso-stack-dir", type=Path, default=None)
    ap.add_argument("--clean-context-stack-dir", type=Path, default=None)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--best-clean-bdt-auc", type=float, default=0.886)
    ap.add_argument("--clean-mlp-auc", type=float, default=0.871)
    ap.add_argument("--clean-nnstack-auc", type=float, default=0.800)
    return ap.parse_args()


def load_json(path: Path | None):
    if path is None or not path.is_file():
        return None
    return json.loads(path.read_text())


def as_float(value) -> float:
    try:
        out = float(value)
    except Exception:
        return math.nan
    return out if math.isfinite(out) else math.nan


def add_bdt_rows(rows: list[dict], report: Path | None) -> None:
    payload = load_json(None if report is None else report / "validation_metrics.json")
    if not payload:
        return
    for product, metrics in payload.get("products", {}).items():
        row = {
            "family": "BDT + isolation inputs",
            "model": product,
            "split": "validation",
            "diagnostic_warning": WARNING,
            "auc": metrics.get("auc_inclusive"),
            "wp80_fake": "",
            "finite_fraction": metrics.get("finite_score_fraction"),
            "score_vs_eiso_corr": metrics.get("score_eiso_pearson"),
        }
        by_cent = metrics.get("auc_by_centrality", {})
        for key, label in CENT_BINS:
            row[f"auc_cent_{label}"] = by_cent.get(key)
        rows.append(row)


def add_mlp_rows(rows: list[dict], report: Path | None) -> None:
    payload = load_json(None if report is None else report / "validation_metrics.json")
    if not payload:
        return
    for product, metrics in payload.get("products", {}).items():
        row = {
            "family": "MLP + isolation inputs",
            "model": product,
            "split": "validation",
            "diagnostic_warning": WARNING,
            "auc": metrics.get("auc"),
            "wp80_fake": (metrics.get("thresholds", {}).get("wp080") or {}).get("background_fake_rate"),
            "finite_fraction": metrics.get("finite_fraction"),
            "score_vs_eiso_corr": (metrics.get("correlations", {}) or {}).get("score_vs_reco_eiso"),
        }
        for item in metrics.get("centrality_bins", []):
            label = f"{item.get('lo'):g}-{item.get('hi'):g}%"
            row[f"auc_cent_{label}"] = item.get("auc")
        rows.append(row)


def read_stack_rank(path: Path) -> list[dict]:
    rank = path / "stacked_sweep_rank_table.csv"
    if not rank.is_file():
        return []
    with rank.open(newline="") as handle:
        return list(csv.DictReader(handle))


def add_stack_rows(rows: list[dict], stack_dir: Path | None, family: str) -> None:
    if stack_dir is None:
        return
    candidates = [row for row in read_stack_rank(stack_dir) if row.get("split") == "test"]
    if not candidates:
        return

    def key(row):
        fake = as_float(row.get("highpt_wp80_fake_20_35"))
        auc = as_float(row.get("auc"))
        return (math.inf if not math.isfinite(fake) else fake, -auc if math.isfinite(auc) else math.inf)

    best = sorted(candidates, key=key)[0]
    row = {
        "family": family,
        "model": best.get("model"),
        "split": "held-out test",
        "diagnostic_warning": WARNING,
        "auc": best.get("auc"),
        "wp80_fake": best.get("wp80_fake"),
        "finite_fraction": best.get("finite_fraction"),
        "score_vs_eiso_corr": best.get("score_vs_eiso_corr"),
        "highpt_auc_20_35": best.get("highpt_auc_20_35"),
        "highpt_wp80_fake_20_35": best.get("highpt_wp80_fake_20_35"),
    }
    for _key, label in CENT_BINS:
        row[f"auc_cent_{label}"] = best.get(f"cent_{_key}_auc")
    rows.append(row)


def write_csv(path: Path, rows: list[dict]) -> None:
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys or ["family", "model"])
        writer.writeheader()
        writer.writerows(rows)


def make_auc_gain_plot(path: Path, rows: list[dict], anchors: dict[str, float]) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
    except Exception:
        return

    selected = [row for row in rows if any(math.isfinite(as_float(row.get(f"auc_cent_{label}"))) for _, label in CENT_BINS)]
    if not selected:
        return
    labels = [label for _, label in CENT_BINS]
    x = np.arange(len(labels), dtype=float)
    width = min(0.22, 0.8 / max(1, len(selected)))
    fig, ax = plt.subplots(figsize=(8.8, 5.2))
    for idx, row in enumerate(selected):
        offset = (idx - 0.5 * (len(selected) - 1)) * width
        values = [as_float(row.get(f"auc_cent_{label}")) for label in labels]
        ax.bar(x + offset, values, width=width, label=row.get("family", row.get("model", "model")))
    for label, anchor in anchors.items():
        if math.isfinite(anchor):
            ax.axhline(anchor, linestyle="--", linewidth=1.2, alpha=0.65, label=label)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylim(0.55, 1.0)
    ax.set_ylabel("AUC")
    ax.set_title("Isolation-visible diagnostic AUC by centrality")
    ax.text(
        0.01,
        0.02,
        "Isolation-derived inputs included; diagnostic ceiling test only, not ABCD-safe.",
        transform=ax.transAxes,
        fontsize=9,
        color="#7a1f1f",
    )
    ax.grid(axis="y", alpha=0.22)
    ax.legend(frameon=False, fontsize=8, ncol=2)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main() -> int:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []
    add_bdt_rows(rows, args.bdt_report)
    add_mlp_rows(rows, args.mlp_report)
    add_stack_rows(rows, args.iso_stack_dir, "NN stack: iso BDT + iso MLP + isolation context")
    add_stack_rows(rows, args.clean_context_stack_dir, "NN stack: clean scores + isolation context")
    write_csv(args.outdir / "iso_visible_diagnostic_rank_table.csv", rows)
    anchors = {
        "best clean BDT anchor": args.best_clean_bdt_auc,
        "clean MLP anchor": args.clean_mlp_auc,
        "clean NN stack anchor": args.clean_nnstack_auc,
    }
    make_auc_gain_plot(args.outdir / "iso_visible_auc_by_centrality.png", rows, anchors)
    summary = [
        "RECOILJETS_AUAU_ISOLATION_VISIBLE_DIAGNOSTIC_SUMMARY_V1",
        f"status={'READY' if rows else 'CHECK'}",
        f"warning={WARNING}",
        f"rank_table={args.outdir / 'iso_visible_diagnostic_rank_table.csv'}",
        f"auc_plot={args.outdir / 'iso_visible_auc_by_centrality.png'}",
    ]
    if rows:
        best = max(rows, key=lambda row: as_float(row.get("auc")))
        summary.append(f"best_by_auc={best.get('family')}::{best.get('model')}::{best.get('auc')}")
    (args.outdir / "iso_visible_diagnostic_summary.txt").write_text("\n".join(summary) + "\n")
    print("\n".join(summary))
    return 0 if rows else 3


if __name__ == "__main__":
    raise SystemExit(main())
