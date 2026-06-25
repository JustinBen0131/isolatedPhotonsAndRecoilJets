#!/usr/bin/env python3
"""Regenerate slide-8-style validation stability check for THE-57 floor veto."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


REPO = Path(__file__).resolve().parents[4]
BASE = REPO / "dataOutput/auauTightBDTValidation/THE57_baseline_validation_20260615"
A1 = BASE / "a1_withcut"
A2 = BASE / "a2_nocut"
A1_SCORE_CACHE_DIR = BASE / "a1_withcut_score_caches/score_caches"
A2_SCORE_CACHE_DIR = BASE / "a2_nocut_score_caches/score_caches"
OUT = BASE / "slide8_floor_veto_validation"
PRODUCT = "centAsFeatBase3x3_pt15to35"
SCORE_KEY = f"score_{PRODUCT}"
CENT_KEYS = ["0_20", "20_50", "50_80"]
CENT_LABELS = ["0-20%", "20-50%", "50-80%"]
CENT_CENTERS = np.array([10.0, 35.0, 65.0])
BROAD_CENT_BINS = [(0.0, 20.0), (20.0, 50.0), (50.0, 80.0)]
FINE_CENT_BINS = [(float(lo), float(lo + 5)) for lo in range(0, 80, 5)]


plt.rcParams.update(
    {
        "font.family": ["Times New Roman", "DejaVu Serif"],
        "axes.titlesize": 16.5,
        "axes.labelsize": 13.5,
        "xtick.labelsize": 11.8,
        "ytick.labelsize": 11.8,
        "legend.fontsize": 11.5,
    }
)


def read_json(path: Path) -> dict:
    with path.open() as f:
        return json.load(f)


def parse_summary(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    for line in path.read_text().splitlines():
        if "=" in line:
            k, v = line.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def metric_float(summary: dict[str, str], key: str) -> float:
    return float(summary[key])


def centrality_metrics(deep: dict) -> dict[str, dict]:
    return deep["products"][PRODUCT]["auc_by_centrality"]


def wp80_fake(row: dict) -> float:
    for item in row["thresholds"]["by_signal_efficiency"]:
        if abs(float(item["target_signal_efficiency"]) - 0.8) < 1e-9:
            return float(item["background_fake_rate"])
    raise KeyError("No target_signal_efficiency=0.8 row")


def load_score_rows(cache_dir: Path) -> dict[str, np.ndarray]:
    paths = sorted(cache_dir.glob("score_cache_*.npz"))
    if not paths:
        raise SystemExit(f"No score caches found under {cache_dir}")
    chunks: dict[str, list[np.ndarray]] = {
        "score": [],
        "is_signal": [],
        "centrality": [],
        "cluster_Et": [],
        "weight": [],
    }
    for path in paths:
        with np.load(path, allow_pickle=True) as data:
            score = np.asarray(data[SCORE_KEY], dtype=float)
            cent = np.asarray(data["centrality"], dtype=float)
            pt = np.asarray(data["cluster_Et"], dtype=float)
            is_sig = np.asarray(data["is_signal"], dtype=bool)
            weight = np.asarray(data["event_weight"], dtype=float) if "event_weight" in data.files else np.ones_like(score, dtype=float)
            keep = np.isfinite(score) & np.isfinite(cent) & np.isfinite(pt) & (pt >= 15.0) & (pt < 35.0) & (cent >= 0.0) & (cent < 80.0)
            chunks["score"].append(score[keep])
            chunks["is_signal"].append(is_sig[keep])
            chunks["centrality"].append(cent[keep])
            chunks["cluster_Et"].append(pt[keep])
            chunks["weight"].append(weight[keep])
    return {key: np.concatenate(vals) for key, vals in chunks.items()}


def weighted_quantile(values: np.ndarray, weights: np.ndarray, q: float) -> float:
    if values.size == 0:
        return float("nan")
    order = np.argsort(values)
    vals = values[order]
    w = weights[order]
    total = float(np.sum(w))
    if total <= 0:
        return float("nan")
    cdf = np.cumsum(w) / total
    return float(vals[min(np.searchsorted(cdf, q, side="left"), vals.size - 1)])


def wp80_fake_for_mask(rows: dict[str, np.ndarray], mask: np.ndarray) -> dict[str, float | int]:
    score = rows["score"]
    is_sig = rows["is_signal"].astype(bool)
    weight = rows["weight"]
    sig = mask & is_sig
    bkg = mask & ~is_sig
    threshold = weighted_quantile(score[sig], weight[sig], 0.20)
    if not np.isfinite(threshold):
        return {
            "threshold": float("nan"),
            "signal_efficiency": float("nan"),
            "background_fake_rate": float("nan"),
            "signal_entries": int(np.count_nonzero(sig)),
            "background_entries": int(np.count_nonzero(bkg)),
        }
    sig_den = float(np.sum(weight[sig]))
    bkg_den = float(np.sum(weight[bkg]))
    sig_pass = float(np.sum(weight[sig & (score > threshold)]))
    bkg_pass = float(np.sum(weight[bkg & (score > threshold)]))
    return {
        "threshold": threshold,
        "signal_efficiency": sig_pass / sig_den if sig_den > 0 else float("nan"),
        "background_fake_rate": bkg_pass / bkg_den if bkg_den > 0 else float("nan"),
        "signal_entries": int(np.count_nonzero(sig)),
        "background_entries": int(np.count_nonzero(bkg)),
        "signal_weight_sum": sig_den,
        "background_weight_sum": bkg_den,
    }


def finite_row_counts(rows: dict[str, np.ndarray], cent_bins: list[tuple[float, float]]) -> list[int]:
    cent = rows["centrality"]
    return [int(np.count_nonzero((cent >= lo) & (cent < hi))) for lo, hi in cent_bins]


def fine_fake_ratio(a1_rows: dict[str, np.ndarray], a2_rows: dict[str, np.ndarray]) -> list[dict]:
    out = []
    for lo, hi in FINE_CENT_BINS:
        a1_mask = (a1_rows["centrality"] >= lo) & (a1_rows["centrality"] < hi)
        a2_mask = (a2_rows["centrality"] >= lo) & (a2_rows["centrality"] < hi)
        a1_wp = wp80_fake_for_mask(a1_rows, a1_mask)
        a2_wp = wp80_fake_for_mask(a2_rows, a2_mask)
        ratio = float(a1_wp["background_fake_rate"]) / float(a2_wp["background_fake_rate"]) if float(a2_wp["background_fake_rate"]) > 0 else float("nan")
        out.append(
            {
                "centrality_min": lo,
                "centrality_max": hi,
                "centrality_center": 0.5 * (lo + hi),
                "with_cut": a1_wp,
                "no_cut": a2_wp,
                "wp80_fake_rate_ratio": ratio,
            }
        )
    return out


def collect_payload() -> dict:
    a1_summary = parse_summary(A1 / "validation_summary.txt")
    a2_summary = parse_summary(A2 / "validation_summary.txt")
    a1_metrics = read_json(A1 / "validation_metrics.json")
    a2_metrics = read_json(A2 / "validation_metrics.json")
    a1_deep = read_json(A1 / "validation_deep_diagnostics.json")
    a2_deep = read_json(A2 / "validation_deep_diagnostics.json")
    a1_rows = load_score_rows(A1_SCORE_CACHE_DIR)
    a2_rows = load_score_rows(A2_SCORE_CACHE_DIR)

    a1_cent = centrality_metrics(a1_deep)
    a2_cent = centrality_metrics(a2_deep)
    a1_broad_counts = finite_row_counts(a1_rows, BROAD_CENT_BINS)
    a2_broad_counts = finite_row_counts(a2_rows, BROAD_CENT_BINS)
    cells = []
    for idx, (key, label) in enumerate(zip(CENT_KEYS, CENT_LABELS)):
        no_cut_entries = a2_broad_counts[idx]
        cut_entries = a1_broad_counts[idx]
        removed = no_cut_entries - cut_entries
        no_cut_fake = wp80_fake(a2_cent[key])
        cut_fake = wp80_fake(a1_cent[key])
        cells.append(
            {
                "centrality_key": key,
                "centrality_label": label,
                "no_cut_auc": float(a2_cent[key]["auc"]),
                "with_cut_auc": float(a1_cent[key]["auc"]),
                "auc_delta": float(a1_cent[key]["auc"]) - float(a2_cent[key]["auc"]),
                "no_cut_wp80_fake_rate": no_cut_fake,
                "with_cut_wp80_fake_rate": cut_fake,
                "wp80_fake_rate_ratio": cut_fake / no_cut_fake if no_cut_fake > 0 else float("nan"),
                "no_cut_validation_entries": no_cut_entries,
                "with_cut_validation_entries": cut_entries,
                "removed_validation_entries": removed,
                "removed_fraction": removed / no_cut_entries if no_cut_entries else float("nan"),
            }
        )

    no_cut_total = sum(c["no_cut_validation_entries"] for c in cells)
    cut_total = sum(c["with_cut_validation_entries"] for c in cells)
    removed_total = no_cut_total - cut_total
    fine = fine_fake_ratio(a1_rows, a2_rows)
    return {
        "schema": "THE57_SLIDE8_FLOOR_VETO_VALIDATION_V1",
        "product": PRODUCT,
        "baseline_model": "diagnostic no-cut model",
        "default_model": "THE-58 floor-veto model",
        "source_dirs": {"with_cut": str(A1), "no_cut": str(A2)},
        "dataset": "local pulled validation outputs; 15 <= cluster_Et < 35 GeV, 0 <= centrality < 80",
        "score_cache_dirs": {"with_cut": str(A1_SCORE_CACHE_DIR), "no_cut": str(A2_SCORE_CACHE_DIR)},
        "floor_veto": {
            "enabled_for_default": True,
            "cut_json_path": a1_metrics["counts"]["event_quality_filter"]["cut_json_path"],
            "cut_json_sha256": a1_metrics["counts"]["event_quality_filter"]["cut_json_sha256"],
        },
        "inclusive": {
            "no_cut_auc": metric_float(a2_summary, f"{PRODUCT}_auc"),
            "with_cut_auc": metric_float(a1_summary, f"{PRODUCT}_auc"),
            "auc_delta": metric_float(a1_summary, f"{PRODUCT}_auc") - metric_float(a2_summary, f"{PRODUCT}_auc"),
            "no_cut_score_mean_signal": metric_float(a2_summary, f"{PRODUCT}_signal_score_mean"),
            "with_cut_score_mean_signal": metric_float(a1_summary, f"{PRODUCT}_signal_score_mean"),
            "no_cut_score_mean_background": metric_float(a2_summary, f"{PRODUCT}_background_score_mean"),
            "with_cut_score_mean_background": metric_float(a1_summary, f"{PRODUCT}_background_score_mean"),
        },
        "validation_rows": {
            "no_cut": no_cut_total,
            "with_cut": cut_total,
            "removed": removed_total,
            "removed_fraction": removed_total / no_cut_total if no_cut_total else float("nan"),
        },
        "cells": cells,
        "fine_fake_ratio": fine,
    }


def card(fig, x, y, w, h, *, title, value, subtitle, face, edge, value_color):
    box = matplotlib.patches.FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.012,rounding_size=0.018",
        transform=fig.transFigure,
        facecolor=face,
        edgecolor=edge,
        linewidth=1.2,
    )
    fig.add_artist(box)
    fig.text(x + 0.018, y + h - 0.023, title, fontsize=12.5, fontweight="bold", color="#14213d", ha="left", va="top")
    fig.text(x + 0.018, y + 0.045, value, fontsize=20, fontweight="bold", color=value_color, ha="left", va="bottom")
    fig.text(x + 0.018, y + 0.018, subtitle, fontsize=10.5, color="#64748b", ha="left", va="bottom")


def render(payload: dict) -> tuple[Path, Path]:
    OUT.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(16, 9), dpi=160)
    fig.patch.set_facecolor("white")
    fig.text(0.047, 0.952, "sPHENIX", fontsize=14.5, fontstyle="italic", fontweight="bold", ha="left", va="top")
    fig.text(0.132, 0.952, "Internal", fontsize=14.5, ha="left", va="top")
    fig.text(0.047, 0.908, "BDT validation is stable after the total-calo floor veto", fontsize=25, fontweight="bold", ha="left", va="top", color="#111827")

    top = matplotlib.patches.FancyBboxPatch(
        (0.055, 0.805),
        0.89,
        0.054,
        boxstyle="round,pad=0.010,rounding_size=0.015",
        transform=fig.transFigure,
        facecolor="#eff6ff",
        edgecolor="#c7d7eb",
        linewidth=1.1,
    )
    fig.add_artist(top)
    fig.text(
        0.5,
        0.833,
        "Validation split, 15-35 GeV photons: no-cut diagnostic vs retrained default after removing low total-calo-energy events.",
        fontsize=12.7,
        fontweight="bold",
        ha="center",
        va="center",
        color="#1f2937",
    )

    rows = payload["validation_rows"]
    inc = payload["inclusive"]
    fine = payload["fine_fake_ratio"]
    ratios = np.array([c["wp80_fake_rate_ratio"] for c in fine], dtype=float)
    card(
        fig,
        0.055,
        0.674,
        0.29,
        0.112,
        title="Validation rows removed by floor veto",
        value=f"{rows['removed']:,}",
        subtitle=f"{100*rows['removed_fraction']:.2f}% of no-cut validation rows",
        face="#fff1f2",
        edge="#fecdd3",
        value_color="#dc2626",
    )
    card(
        fig,
        0.360,
        0.674,
        0.29,
        0.112,
        title="Inclusive AUC",
        value=f"{inc['with_cut_auc']:.6f}",
        subtitle=f"no-cut {inc['no_cut_auc']:.6f}; change {inc['auc_delta']:+.6f}",
        face="#ecfdf5",
        edge="#bbf7d0",
        value_color="#059669",
    )
    card(
        fig,
        0.665,
        0.674,
        0.28,
        0.112,
        title="Inclusive WP80 fake-rate ratio",
        value=f"{np.nanmean(ratios):.3f}",
        subtitle=f"with cut / no cut in 5% bins; range {np.nanmin(ratios):.3f}-{np.nanmax(ratios):.3f}",
        face="#eff6ff",
        edge="#bfdbfe",
        value_color="#2563eb",
    )

    ax1 = fig.add_axes([0.075, 0.365, 0.405, 0.245])
    ax2 = fig.add_axes([0.555, 0.365, 0.37, 0.245])
    x = np.arange(len(CENT_LABELS))
    a2_auc = np.array([c["no_cut_auc"] for c in payload["cells"]])
    a1_auc = np.array([c["with_cut_auc"] for c in payload["cells"]])
    width = 0.34
    ax1.bar(x - width / 2, a2_auc, width, color="#b9c9da", edgecolor="#64748b", label="No-cut diagnostic")
    ax1.bar(x + width / 2, a1_auc, width, color="#1f77b4", edgecolor="#174b73", label="With floor veto + retrain")
    for i, delta in enumerate(a1_auc - a2_auc):
        ax1.text(i, max(a1_auc[i], a2_auc[i]) + 0.0032, f"{delta:+.4f}", ha="center", va="bottom", fontsize=11.5, fontweight="bold", color="#334155")
    ax1.set_xticks(x, CENT_LABELS)
    ax1.set_ylabel("AUC")
    ax1.set_title("AUC by centrality", fontweight="bold")
    ax1.set_ylim(min(np.min(a1_auc), np.min(a2_auc)) - 0.012, max(np.max(a1_auc), np.max(a2_auc)) + 0.012)
    ax1.grid(axis="y", alpha=0.28)
    ax1.legend(loc="upper left", frameon=False)

    fine_centers = np.array([c["centrality_center"] for c in fine], dtype=float)
    ax2.axhline(1.0, color="#94a3b8", lw=1.4, linestyle="--", label="No change")
    ax2.plot(fine_centers, ratios, color="#d17a16", marker="o", markersize=5.5, lw=2.0, label="With floor veto / no cut")
    idx_hi = int(np.nanargmax(ratios))
    idx_lo = int(np.nanargmin(ratios))
    for idx in sorted(set([idx_hi, idx_lo])):
        xval = fine_centers[idx]
        ratio = ratios[idx]
        ax2.text(xval, ratio + (0.004 if ratio <= 1 else 0.003), f"{(ratio-1)*100:+.2f}%", ha="center", va="bottom", fontsize=10.4, fontweight="bold", color="#334155")
    ax2.set_title("Fine-bin fake-rate ratio with WP80 rederived per 5%", fontweight="bold")
    ax2.set_xlabel("Centrality percentile")
    ax2.set_ylabel("Fake-rate ratio")
    ax2.set_xlim(0, 80)
    pad = max(0.015, float(np.nanmax(np.abs(ratios - 1.0))) + 0.008)
    ax2.set_ylim(1.0 - pad, 1.0 + pad)
    ax2.grid(alpha=0.28)
    ax2.legend(loc="upper right", frameon=False)

    bottom = matplotlib.patches.FancyBboxPatch(
        (0.055, 0.075),
        0.89,
        0.205,
        boxstyle="round,pad=0.010,rounding_size=0.018",
        transform=fig.transFigure,
        facecolor="#fffafa",
        edgecolor="#fecaca",
        linewidth=1.1,
    )
    fig.add_artist(bottom)
    fig.text(0.075, 0.245, "15-35 GeV validation rows removed in each centrality bin", fontsize=15.2, fontweight="bold", ha="left", va="top", color="#111827")
    card_w = 0.225
    for idx, cell in enumerate(payload["cells"]):
        x0 = 0.085 + idx * 0.31
        inner = matplotlib.patches.FancyBboxPatch(
            (x0, 0.105),
            card_w,
            0.095,
            boxstyle="round,pad=0.010,rounding_size=0.010",
            transform=fig.transFigure,
            facecolor="white",
            edgecolor="#dbeafe",
            linewidth=1.0,
        )
        fig.add_artist(inner)
        fig.text(x0 + card_w / 2, 0.178, cell["centrality_label"], fontsize=16, fontweight="bold", ha="center", va="top", color="#1f2937")
        fig.text(
            x0 + card_w / 2,
            0.142,
            f"{cell['removed_validation_entries']:,} / {cell['no_cut_validation_entries']:,} rows removed",
            fontsize=11.2,
            fontweight="bold",
            ha="center",
            va="top",
            color="#dc2626" if cell["removed_fraction"] > 0.02 else "#64748b",
        )
        fig.text(x0 + card_w / 2, 0.116, f"{100*cell['removed_fraction']:.2f}% of that validation bin", fontsize=10.2, ha="center", va="top", color="#64748b")

    png = OUT / "the57_slide8_floor_veto_validation.png"
    script = OUT / "the57_slide8_floor_veto_validation_script.md"
    fig.savefig(png, facecolor="white")
    plt.close(fig)
    script.write_text(
        "# Slide 8 spoken script\n\n"
        "This slide answers whether the THE-58 event-level floor veto destabilized the BDT validation.\n"
        "The comparison is between the no-cut diagnostic training and the retrained default model that uses the floor veto before training.\n"
        "The inclusive AUC changes only at the few-times-10^-4 level, and the centrality-binned AUC bars remain essentially overlapping.\n"
        "The fake-rate ratio is near unity, so the veto is not producing a large artificial gain or loss in the validation working point.\n"
        "The bottom row records how many validation rows are removed in each broad centrality bin for the 15-35 GeV scored validation set.\n"
    )
    return png, script


def main() -> None:
    payload = collect_payload()
    png, script = render(payload)
    payload_path = OUT / "the57_slide8_floor_veto_validation_payload.json"
    manifest_path = OUT / "the57_slide8_floor_veto_validation_manifest.json"
    payload_path.write_text(json.dumps(payload, indent=2) + "\n")
    manifest = {
        "schema": "THE57_SLIDE8_FLOOR_VETO_VALIDATION_MANIFEST_V1",
        "status": "READY",
        "png": str(png),
        "speaker_script": str(script),
        "payload": str(payload_path),
        "source_dirs": payload["source_dirs"],
        "score_cache_dirs": payload["score_cache_dirs"],
        "dataset": payload["dataset"],
        "baseline_model": payload["baseline_model"],
        "default_model": payload["default_model"],
        "floor_veto": payload["floor_veto"],
        "inclusive": payload["inclusive"],
        "validation_rows": payload["validation_rows"],
        "fine_fake_ratio": payload["fine_fake_ratio"],
        "google_slides_mutated": False,
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print("THE57_SLIDE8_FLOOR_VETO_VALIDATION_READY")
    print(f"png={png}")
    print(f"script={script}")
    print(f"payload={payload_path}")
    print(f"manifest={manifest_path}")


if __name__ == "__main__":
    main()
