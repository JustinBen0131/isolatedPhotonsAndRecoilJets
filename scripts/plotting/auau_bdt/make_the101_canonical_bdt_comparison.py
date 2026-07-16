#!/usr/bin/env python3
"""Compare the THE-101 canonical AuAu BDT candidate with THE-95."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


PRODUCT = "centAsFeatBase3x3_pt15to35"
CENTRALITY_KEYS = ("0_20", "20_50", "50_80")
BLUE = "#1565C0"
ORANGE = "#E07A1F"
GRAY = "#5F6368"
GREEN = "#1B8A5A"


def load_json(path: Path) -> dict:
    with path.open() as stream:
        return json.load(stream)


def product_metrics(payload: dict) -> dict:
    if payload.get("status") != "READY":
        raise ValueError(f"validation is not READY: {payload.get('status')}")
    return payload["products"][PRODUCT]


def wp80_rows(payload: dict, key: str) -> list[dict]:
    rows = [row for row in payload[key] if row.get("wp_label") == "WP80"]
    if len(rows) != 7:
        raise ValueError(f"expected seven WP80 {key}, found {len(rows)}")
    return sorted(rows, key=lambda row: float(row["centrality_center"]))


def weighted_line(rows: list[dict]) -> dict:
    x = np.asarray([row["centrality_center"] for row in rows], dtype=float)
    y = np.asarray([row["threshold"] for row in rows], dtype=float)
    err = np.asarray([row["threshold_stat_err"] for row in rows], dtype=float)
    good = np.isfinite(err) & (err > 0)
    if not good.all():
        raise ValueError("WP80 flat-fit rows contain invalid uncertainties")
    slope, intercept = np.polyfit(x, y, 1, w=1.0 / err)
    residual = y - (intercept + slope * x)
    return {
        "intercept": float(intercept),
        "slope_per_percent": float(slope),
        "max_abs_residual": float(np.max(np.abs(residual))),
        "rms_residual": float(np.sqrt(np.mean(residual**2))),
    }


def direct_fake_rates(payload: dict) -> tuple[np.ndarray, np.ndarray]:
    rows = wp80_rows(payload, "rows")
    centers = np.asarray([row["centrality_center"] for row in rows], dtype=float)
    rates = np.asarray([row["background_fake_rate"] for row in rows], dtype=float)
    return centers, rates


def normalized_training_contract(metadata: dict) -> dict:
    split = metadata["split"]
    closure = metadata["ppg12_exact_closure"]
    weighting = closure["weighting"]
    return {
        "campaign": metadata["campaign"],
        "product": metadata["product"],
        "features": metadata["features"],
        "pt_range": metadata["pt_range"],
        "label_branch": metadata["label_branch"],
        "task": metadata["task"],
        "weight_mode": metadata["weight_mode"],
        "xgboost": metadata["xgboost"],
        "background_subsampling": metadata["background_subsampling"],
        "majority_class_optimization": metadata["majority_class_optimization"],
        "split_mode": split["mode"],
        "test_fraction_requested": split["test_fraction_requested"],
        "source_samples": closure["sample_validation"]["observed_samples"],
        "weighting_policy": {
            key: weighting[key]
            for key in (
                "centrality_event_weight",
                "cross_section_weight_used_for_training",
                "event_weight_used",
                "source",
                "vertex_reweight",
                "weight_mode",
                "weights_computed_before_binning",
            )
        },
    }


def build_summary(
    candidate_validation: dict,
    baseline_validation: dict,
    candidate_wp: dict,
    baseline_same_source_wp: dict,
    baseline_published_wp: dict,
    candidate_metadata: dict,
    baseline_metadata: dict,
) -> dict:
    cand = product_metrics(candidate_validation)
    base = product_metrics(baseline_validation)
    if candidate_validation["counts"] != baseline_validation["counts"]:
        raise ValueError("candidate and baseline validation counts differ")

    cand_flat = wp80_rows(candidate_wp, "flat_rows")
    base_flat = wp80_rows(baseline_same_source_wp, "flat_rows")
    published_flat = wp80_rows(baseline_published_wp, "flat_rows")
    cand_fit = weighted_line(cand_flat)
    base_fit = weighted_line(base_flat)
    published_fit = weighted_line(published_flat)
    candidate_contract = normalized_training_contract(candidate_metadata)
    baseline_contract = normalized_training_contract(baseline_metadata)

    auc_delta = float(cand["auc_inclusive"] - base["auc_inclusive"])
    cent_deltas = {
        key: float(cand["auc_by_centrality"][key] - base["auc_by_centrality"][key])
        for key in CENTRALITY_KEYS
    }
    cand_cent, cand_fake = direct_fake_rates(candidate_wp)
    base_cent, base_fake = direct_fake_rates(baseline_same_source_wp)
    if not np.array_equal(cand_cent, base_cent):
        raise ValueError("candidate and baseline WP80 centrality grids differ")
    fake_delta = cand_fake - base_fake
    score_behavior = {
        "candidate_signal_mean": float(cand["signal_score_mean"]),
        "the95_same_source_signal_mean": float(base["signal_score_mean"]),
        "signal_mean_delta": float(cand["signal_score_mean"] - base["signal_score_mean"]),
        "candidate_background_mean": float(cand["background_score_mean"]),
        "the95_same_source_background_mean": float(base["background_score_mean"]),
        "background_mean_delta": float(
            cand["background_score_mean"] - base["background_score_mean"]
        ),
        "candidate_score_eiso_pearson": float(cand["score_eiso_pearson"]),
        "the95_same_source_score_eiso_pearson": float(base["score_eiso_pearson"]),
        "score_eiso_pearson_delta": float(
            cand["score_eiso_pearson"] - base["score_eiso_pearson"]
        ),
    }

    gates = {
        "identical_validation_counts": candidate_validation["counts"] == baseline_validation["counts"],
        "training_contract_unchanged": candidate_contract == baseline_contract,
        "finite_score_fraction_is_one": cand["finite_score_fraction"] == 1.0
        and base["finite_score_fraction"] == 1.0,
        "inclusive_auc_abs_delta_le_0p005": abs(auc_delta) <= 0.005,
        "centrality_auc_max_abs_delta_le_0p005": max(abs(v) for v in cent_deltas.values()) <= 0.005,
        "score_behavior_max_abs_delta_le_0p03": max(
            abs(score_behavior[key])
            for key in (
                "signal_mean_delta",
                "background_mean_delta",
                "score_eiso_pearson_delta",
            )
        )
        <= 0.03,
        "wp80_fit_max_abs_delta_le_0p03": max(
            abs(
                (cand_fit["intercept"] + cand_fit["slope_per_percent"] * c)
                - (base_fit["intercept"] + base_fit["slope_per_percent"] * c)
            )
            for c in np.linspace(0.0, 80.0, 81)
        )
        <= 0.03,
        "wp80_published_fit_max_abs_delta_le_0p03": max(
            abs(
                (cand_fit["intercept"] + cand_fit["slope_per_percent"] * c)
                - (published_fit["intercept"] + published_fit["slope_per_percent"] * c)
            )
            for c in np.linspace(0.0, 80.0, 81)
        )
        <= 0.03,
        "wp80_background_fake_rate_max_increase_le_0p02": float(np.max(fake_delta)) <= 0.02,
        "candidate_wp80_fit_max_residual_le_0p03": cand_fit["max_abs_residual"] <= 0.03,
    }
    gates = {name: bool(passed) for name, passed in gates.items()}

    return {
        "schema": "THE101_CANONICAL_AUAU_BDT_COMPARISON_V1",
        "status": "PASS" if all(gates.values()) else "FAIL",
        "model_product": PRODUCT,
        "comparison_scope": "same classifier-pass six-sample extraction",
        "validation_counts": candidate_validation["counts"],
        "training_contract": {
            "unchanged": candidate_contract == baseline_contract,
            "candidate": candidate_contract,
            "the95": baseline_contract,
            "candidate_rows": candidate_metadata["n_rows"],
            "the95_rows": baseline_metadata["n_rows"],
            "row_count_change_reason": "embedded MinimumBiasClassifier pass requirement",
        },
        "auc": {
            "candidate_inclusive": cand["auc_inclusive"],
            "the95_same_source_inclusive": base["auc_inclusive"],
            "inclusive_delta": auc_delta,
            "candidate_by_centrality": cand["auc_by_centrality"],
            "the95_same_source_by_centrality": base["auc_by_centrality"],
            "delta_by_centrality": cent_deltas,
        },
        "score_behavior": score_behavior,
        "wp80": {
            "candidate_fit": cand_fit,
            "the95_same_source_fit": base_fit,
            "the95_published_fit": published_fit,
            "candidate_background_fake_rate": cand_fake.tolist(),
            "the95_same_source_background_fake_rate": base_fake.tolist(),
            "background_fake_rate_delta": fake_delta.tolist(),
            "centrality_centers": cand_cent.tolist(),
        },
        "gates": gates,
    }


def write_csv(summary: dict, path: Path) -> None:
    auc = summary["auc"]
    wp = summary["wp80"]
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["metric", "bin", "candidate", "the95_same_source", "delta"])
        writer.writerow(
            [
                "auc",
                "inclusive",
                auc["candidate_inclusive"],
                auc["the95_same_source_inclusive"],
                auc["inclusive_delta"],
            ]
        )
        for key in CENTRALITY_KEYS:
            writer.writerow(
                [
                    "auc",
                    key,
                    auc["candidate_by_centrality"][key],
                    auc["the95_same_source_by_centrality"][key],
                    auc["delta_by_centrality"][key],
                ]
            )
        for metric, candidate_key, baseline_key, delta_key in (
            (
                "signal_score_mean",
                "candidate_signal_mean",
                "the95_same_source_signal_mean",
                "signal_mean_delta",
            ),
            (
                "background_score_mean",
                "candidate_background_mean",
                "the95_same_source_background_mean",
                "background_mean_delta",
            ),
            (
                "score_eiso_pearson",
                "candidate_score_eiso_pearson",
                "the95_same_source_score_eiso_pearson",
                "score_eiso_pearson_delta",
            ),
        ):
            score = summary["score_behavior"]
            writer.writerow(
                [metric, "inclusive", score[candidate_key], score[baseline_key], score[delta_key]]
            )
        for center, cand, base, delta in zip(
            wp["centrality_centers"],
            wp["candidate_background_fake_rate"],
            wp["the95_same_source_background_fake_rate"],
            wp["background_fake_rate_delta"],
        ):
            writer.writerow(["wp80_background_fake_rate", center, cand, base, delta])


def add_sphenix_label(ax: plt.Axes) -> None:
    ax.text(
        0.03,
        0.96,
        "sPHENIX",
        transform=ax.transAxes,
        va="top",
        fontsize=13,
        fontweight="bold",
        fontstyle="italic",
    )
    ax.text(0.285, 0.96, "Internal", transform=ax.transAxes, va="top", fontsize=13)


def plot(summary: dict, candidate_wp: dict, baseline_wp: dict, published_wp: dict, path: Path) -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.size": 11,
            "axes.linewidth": 1.0,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 5.0), constrained_layout=True)
    auc = summary["auc"]
    labels = ["Inclusive", "0-20%", "20-50%", "50-80%"]
    cand_auc = [auc["candidate_inclusive"]] + [auc["candidate_by_centrality"][k] for k in CENTRALITY_KEYS]
    base_auc = [auc["the95_same_source_inclusive"]] + [auc["the95_same_source_by_centrality"][k] for k in CENTRALITY_KEYS]
    x = np.arange(len(labels))
    axes[0].plot(x, base_auc, "o", color=GRAY, ms=7, label="THE-95 on same events")
    axes[0].plot(x, cand_auc, "s", color=BLUE, ms=7, label="Classifier-pass retrain")
    axes[0].set_xticks(x, labels, rotation=20, ha="right")
    axes[0].set_ylabel("ROC AUC")
    axes[0].set_ylim(0.835, 0.885)
    axes[0].set_title("Same-source discrimination")
    axes[0].legend(frameon=False, loc="lower right", fontsize=10)
    add_sphenix_label(axes[0])

    for payload, color, marker, label in (
        (published_wp, ORANGE, "^", "Published THE-95"),
        (baseline_wp, GRAY, "o", "THE-95 on classifier-pass events"),
        (candidate_wp, BLUE, "s", "Classifier-pass retrain"),
    ):
        rows = wp80_rows(payload, "flat_rows")
        centers = np.asarray([row["centrality_center"] for row in rows], dtype=float)
        thresholds = np.asarray([row["threshold"] for row in rows], dtype=float)
        errors = np.asarray([row["threshold_stat_err"] for row in rows], dtype=float)
        fit = weighted_line(rows)
        axes[1].errorbar(centers, thresholds, yerr=errors, fmt=marker, color=color, ms=6, capsize=2, label=label)
        grid = np.linspace(0.0, 80.0, 161)
        axes[1].plot(grid, fit["intercept"] + fit["slope_per_percent"] * grid, color=color, lw=1.8)
    axes[1].set_xlabel("Centrality percentile [%]")
    axes[1].set_ylabel("BDT threshold at 80% signal efficiency")
    axes[1].set_xlim(0, 80)
    axes[1].set_title("Centrality-dependent WP80")
    axes[1].legend(frameon=False, fontsize=9, loc="lower right")

    centers = np.asarray(summary["wp80"]["centrality_centers"])
    axes[2].plot(
        centers,
        summary["wp80"]["the95_same_source_background_fake_rate"],
        "o-",
        color=GRAY,
        lw=1.5,
        ms=6,
        label="THE-95 on same events",
    )
    axes[2].plot(
        centers,
        summary["wp80"]["candidate_background_fake_rate"],
        "s-",
        color=GREEN,
        lw=1.5,
        ms=6,
        label="Classifier-pass retrain",
    )
    axes[2].set_xlabel("Centrality percentile [%]")
    axes[2].set_ylabel("Inclusive-jet acceptance at local WP80")
    axes[2].set_xlim(0, 80)
    axes[2].set_title("Background acceptance")
    axes[2].legend(frameon=False, fontsize=10, loc="upper right")

    fig.suptitle("Canonical Au+Au photon-ID BDT validation", fontsize=18, fontweight="bold")
    fig.text(
        0.5,
        0.925,
        "Embedded photon+jet 12+20 and inclusive-jet 12+20+30+40; identical classifier-pass source",
        ha="center",
        fontsize=11.5,
        color="#3C4043",
    )
    fig.savefig(path, dpi=180, facecolor="white")
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate-validation", type=Path, required=True)
    parser.add_argument("--baseline-validation", type=Path, required=True)
    parser.add_argument("--candidate-wp", type=Path, required=True)
    parser.add_argument("--baseline-same-source-wp", type=Path, required=True)
    parser.add_argument("--baseline-published-wp", type=Path, required=True)
    parser.add_argument("--candidate-metadata", type=Path, required=True)
    parser.add_argument("--baseline-metadata", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    candidate_validation = load_json(args.candidate_validation)
    baseline_validation = load_json(args.baseline_validation)
    candidate_wp = load_json(args.candidate_wp)
    baseline_same_source_wp = load_json(args.baseline_same_source_wp)
    baseline_published_wp = load_json(args.baseline_published_wp)
    candidate_metadata = load_json(args.candidate_metadata)
    baseline_metadata = load_json(args.baseline_metadata)
    summary = build_summary(
        candidate_validation,
        baseline_validation,
        candidate_wp,
        baseline_same_source_wp,
        baseline_published_wp,
        candidate_metadata,
        baseline_metadata,
    )
    json_path = args.out_dir / "the101_canonical_bdt_acceptance.json"
    csv_path = args.out_dir / "the101_canonical_bdt_comparison.csv"
    png_path = args.out_dir / "the101_canonical_bdt_comparison.png"
    json_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    write_csv(summary, csv_path)
    plot(summary, candidate_wp, baseline_same_source_wp, baseline_published_wp, png_path)
    print(f"status={summary['status']}")
    print(f"json={json_path}")
    print(f"csv={csv_path}")
    print(f"png={png_path}")


if __name__ == "__main__":
    main()
