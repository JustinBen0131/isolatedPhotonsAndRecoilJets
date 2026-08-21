#!/usr/bin/env python3
"""Compare corrected and historical Au+Au BDTs on identical corrected rows."""

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
NAVY = "#163A5F"
GREEN = "#1B8A5A"
ORANGE = "#D97706"
GRAY = "#667085"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text())


def metrics(payload: dict) -> dict:
    if payload.get("status") != "READY":
        raise ValueError(f"validation is not READY: {payload.get('status')}")
    return payload["products"][PRODUCT]


def wp80_rows(payload: dict, key: str) -> list[dict]:
    rows = [row for row in payload[key] if row.get("wp_label") == "WP80"]
    if len(rows) != 7:
        raise ValueError(f"expected seven WP80 rows in {key}; found {len(rows)}")
    return sorted(rows, key=lambda row: float(row["centrality_center"]))


def weighted_line(rows: list[dict]) -> dict:
    x = np.asarray([row["centrality_center"] for row in rows], dtype=float)
    y = np.asarray([row["threshold"] for row in rows], dtype=float)
    err = np.asarray([row["threshold_stat_err"] for row in rows], dtype=float)
    if not (np.isfinite(err) & (err > 0)).all():
        raise ValueError("working-point rows contain invalid uncertainties")
    slope, intercept = np.polyfit(x, y, 1, w=1.0 / err)
    residual = y - (intercept + slope * x)
    return {
        "intercept": float(intercept),
        "slope_per_percent": float(slope),
        "rms_residual": float(np.sqrt(np.mean(residual**2))),
        "max_abs_residual": float(np.max(np.abs(residual))),
    }


def training_contract(metadata: dict) -> dict:
    split = metadata["split"]
    closure = metadata["ppg12_exact_closure"]
    weighting = closure["weighting"]
    # THE-101's trainer used these executable defaults but its older metadata
    # did not serialize all of them.  The corrected trainer records them
    # explicitly.  Normalize the two metadata generations before testing the
    # frozen training contract so the gate measures parameter values rather
    # than provenance verbosity.
    raw_xgboost = metadata["xgboost"]
    xgboost_defaults = {
        "objective": "binary:logistic",
        "eval_metric": ["auc", "logloss"],
        "random_state": 13,
        "n_jobs": 1,
    }
    xgboost_keys = (
        "n_estimators",
        "max_depth",
        "learning_rate",
        "subsample",
        "colsample_bytree",
        "reg_alpha",
        "reg_lambda",
        "grow_policy",
        "max_bin",
        "tree_method",
        "objective",
        "eval_metric",
        "random_state",
        "n_jobs",
    )
    normalized_xgboost = {
        key: raw_xgboost.get(key, xgboost_defaults.get(key)) for key in xgboost_keys
    }
    return {
        "campaign": metadata["campaign"],
        "product": metadata["product"],
        "features": metadata["features"],
        "pt_range": metadata["pt_range"],
        "label_branch": metadata["label_branch"],
        "task": metadata["task"],
        "weight_mode": metadata["weight_mode"],
        "xgboost": normalized_xgboost,
        "background_subsampling": metadata["background_subsampling"],
        "majority_class_optimization": metadata["majority_class_optimization"],
        "split_mode": split["mode"],
        "test_fraction_requested": split["test_fraction_requested"],
        "split_random_seed": split.get("random_seed", 13),
        "split_stratified_by_class": split.get("stratified_by_class", True),
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


def fake_rates(payload: dict) -> tuple[np.ndarray, np.ndarray]:
    rows = wp80_rows(payload, "rows")
    return (
        np.asarray([row["centrality_center"] for row in rows], dtype=float),
        np.asarray([row["background_fake_rate"] for row in rows], dtype=float),
    )


def build_summary(args: argparse.Namespace) -> dict:
    new_validation = load_json(args.corrected_validation)
    old_validation = load_json(args.historical_validation)
    new_wp = load_json(args.corrected_wp)
    old_same_wp = load_json(args.historical_same_source_wp)
    old_published_wp = load_json(args.historical_published_wp)
    new_metadata = load_json(args.corrected_metadata)
    old_metadata = load_json(args.historical_metadata)
    extraction_audit = load_json(args.extraction_audit)

    new = metrics(new_validation)
    old = metrics(old_validation)
    if new_validation["counts"] != old_validation["counts"]:
        raise ValueError("models were not validated on identical corrected rows")

    new_contract = training_contract(new_metadata)
    old_contract = training_contract(old_metadata)
    auc_delta = float(new["auc_inclusive"] - old["auc_inclusive"])
    cent_delta = {
        key: float(new["auc_by_centrality"][key] - old["auc_by_centrality"][key])
        for key in CENTRALITY_KEYS
    }
    new_cent, new_fake = fake_rates(new_wp)
    old_cent, old_fake = fake_rates(old_same_wp)
    if not np.array_equal(new_cent, old_cent):
        raise ValueError("working-point centrality grids differ")

    new_fit = weighted_line(wp80_rows(new_wp, "flat_rows"))
    old_same_fit = weighted_line(wp80_rows(old_same_wp, "flat_rows"))
    old_published_fit = weighted_line(wp80_rows(old_published_wp, "flat_rows"))
    score_deltas = {
        "signal_mean": float(new["signal_score_mean"] - old["signal_score_mean"]),
        "background_mean": float(new["background_score_mean"] - old["background_score_mean"]),
        "score_eiso_pearson": float(
            new["score_eiso_pearson"] - old["score_eiso_pearson"]
        ),
    }
    grid = np.linspace(0.0, 80.0, 81)
    fit_delta = (
        new_fit["intercept"]
        + new_fit["slope_per_percent"] * grid
        - old_same_fit["intercept"]
        - old_same_fit["slope_per_percent"] * grid
    )
    gates = {
        "corrected_extraction_audit_passed": extraction_audit.get("status")
        in {"PASS", "PASSED"},
        "identical_validation_counts": new_validation["counts"]
        == old_validation["counts"],
        "training_contract_unchanged": new_contract == old_contract,
        "finite_score_fraction_is_one": new["finite_score_fraction"] == 1.0
        and old["finite_score_fraction"] == 1.0,
        "new_wp80_fit_residual_le_0p03": new_fit["max_abs_residual"] <= 0.03,
        "metrics_are_finite": all(
            math.isfinite(value)
            for value in [
                new["auc_inclusive"],
                old["auc_inclusive"],
                *cent_delta.values(),
                *score_deltas.values(),
                *new_fake,
                *old_fake,
            ]
        ),
    }
    gates = {key: bool(value) for key, value in gates.items()}
    return {
        "schema": "CORRECTED_AUAU_SHOWER_CONTRACT_BDT_COMPARISON_V1",
        "status": "PASS" if all(gates.values()) else "FAIL",
        "promotion_status": "NOT_PROMOTED",
        "comparison_scope": "both models evaluated on identical corrected-feature rows",
        "validation_counts": new_validation["counts"],
        "training_contract": {
            "unchanged": new_contract == old_contract,
            "corrected": new_contract,
            "historical": old_contract,
            "corrected_rows": new_metadata["n_rows"],
            "historical_rows": old_metadata["n_rows"],
        },
        "auc": {
            "corrected_inclusive": new["auc_inclusive"],
            "historical_inclusive": old["auc_inclusive"],
            "delta": auc_delta,
            "corrected_by_centrality": new["auc_by_centrality"],
            "historical_by_centrality": old["auc_by_centrality"],
            "delta_by_centrality": cent_delta,
        },
        "score_behavior": {
            "corrected_signal_mean": new["signal_score_mean"],
            "historical_signal_mean": old["signal_score_mean"],
            "corrected_background_mean": new["background_score_mean"],
            "historical_background_mean": old["background_score_mean"],
            "corrected_score_eiso_pearson": new["score_eiso_pearson"],
            "historical_score_eiso_pearson": old["score_eiso_pearson"],
            "deltas": score_deltas,
        },
        "wp80": {
            "corrected_fit": new_fit,
            "historical_same_source_fit": old_same_fit,
            "historical_published_fit": old_published_fit,
            "centrality_centers": new_cent.tolist(),
            "corrected_background_fake_rate": new_fake.tolist(),
            "historical_background_fake_rate": old_fake.tolist(),
            "background_fake_rate_delta": (new_fake - old_fake).tolist(),
            "max_abs_same_source_fit_delta": float(np.max(np.abs(fit_delta))),
        },
        "gates": gates,
    }


def write_csv(summary: dict, path: Path) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["metric", "bin", "corrected", "historical", "delta"])
        writer.writerow(
            [
                "auc",
                "inclusive",
                summary["auc"]["corrected_inclusive"],
                summary["auc"]["historical_inclusive"],
                summary["auc"]["delta"],
            ]
        )
        for key in CENTRALITY_KEYS:
            writer.writerow(
                [
                    "auc",
                    key,
                    summary["auc"]["corrected_by_centrality"][key],
                    summary["auc"]["historical_by_centrality"][key],
                    summary["auc"]["delta_by_centrality"][key],
                ]
            )
        wp = summary["wp80"]
        for center, new, old, delta in zip(
            wp["centrality_centers"],
            wp["corrected_background_fake_rate"],
            wp["historical_background_fake_rate"],
            wp["background_fake_rate_delta"],
            strict=True,
        ):
            writer.writerow(["wp80_background_acceptance", center, new, old, delta])


def sphenix_label(ax: plt.Axes) -> None:
    ax.text(
        0.04,
        0.96,
        r"$\bf{\it{sPHENIX}}$ Internal",
        transform=ax.transAxes,
        va="top",
        fontsize=12,
    )


def render(summary: dict, corrected_wp: dict, historical_wp: dict, path: Path) -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.size": 10.5,
            "axes.linewidth": 1.1,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.7), constrained_layout=True)
    auc = summary["auc"]
    labels = ["Inclusive", "0–20%", "20–50%", "50–80%"]
    new_auc = [auc["corrected_inclusive"]] + [
        auc["corrected_by_centrality"][key] for key in CENTRALITY_KEYS
    ]
    old_auc = [auc["historical_inclusive"]] + [
        auc["historical_by_centrality"][key] for key in CENTRALITY_KEYS
    ]
    x = np.arange(len(labels))
    axes[0].plot(x, old_auc, "o", color=GRAY, ms=6, label="Historical model")
    axes[0].plot(x, new_auc, "s", color=NAVY, ms=6, label="Corrected-input retrain")
    axes[0].set_xticks(x, labels, rotation=18, ha="right")
    axes[0].set_ylabel("ROC AUC")
    axes[0].set_title("Held-out discrimination")
    axes[0].legend(frameon=False, loc="lower right")
    sphenix_label(axes[0])

    for payload, color, marker, label in (
        (historical_wp, GRAY, "o", "Historical model on corrected rows"),
        (corrected_wp, GREEN, "s", "Corrected-input retrain"),
    ):
        rows = wp80_rows(payload, "flat_rows")
        centers = np.asarray([row["centrality_center"] for row in rows])
        thresholds = np.asarray([row["threshold"] for row in rows])
        errors = np.asarray([row["threshold_stat_err"] for row in rows])
        fit = weighted_line(rows)
        axes[1].errorbar(
            centers,
            thresholds,
            yerr=errors,
            fmt=marker,
            color=color,
            ms=5.5,
            capsize=2,
            label=label,
        )
        grid = np.linspace(0, 80, 161)
        axes[1].plot(
            grid,
            fit["intercept"] + fit["slope_per_percent"] * grid,
            color=color,
            lw=1.7,
        )
    axes[1].set_xlabel("Centrality percentile [%]")
    axes[1].set_ylabel("BDT threshold at 80% signal efficiency")
    axes[1].set_title("Re-derived WP80")
    axes[1].legend(frameon=False, fontsize=8.5)

    centers = np.asarray(summary["wp80"]["centrality_centers"])
    axes[2].plot(
        centers,
        summary["wp80"]["historical_background_fake_rate"],
        "o-",
        color=GRAY,
        ms=5.5,
        label="Historical model",
    )
    axes[2].plot(
        centers,
        summary["wp80"]["corrected_background_fake_rate"],
        "s-",
        color=ORANGE,
        ms=5.5,
        label="Corrected-input retrain",
    )
    axes[2].set_xlabel("Centrality percentile [%]")
    axes[2].set_ylabel("Inclusive-jet acceptance at local WP80")
    axes[2].set_title("Background acceptance")
    axes[2].legend(frameon=False, fontsize=8.5)

    fig.suptitle(
        "Au+Au photon-ID response to the corrected shower-feature contract",
        fontsize=15,
        fontweight="bold",
    )
    fig.savefig(path, dpi=200, facecolor="white")
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--corrected-validation", type=Path, required=True)
    parser.add_argument("--historical-validation", type=Path, required=True)
    parser.add_argument("--corrected-wp", type=Path, required=True)
    parser.add_argument("--historical-same-source-wp", type=Path, required=True)
    parser.add_argument("--historical-published-wp", type=Path, required=True)
    parser.add_argument("--corrected-metadata", type=Path, required=True)
    parser.add_argument("--historical-metadata", type=Path, required=True)
    parser.add_argument("--extraction-audit", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    summary = build_summary(args)
    json_path = args.out_dir / "corrected_shower_contract_bdt_comparison.json"
    csv_path = args.out_dir / "corrected_shower_contract_bdt_comparison.csv"
    png_path = args.out_dir / "corrected_shower_contract_bdt_comparison.png"
    json_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    write_csv(summary, csv_path)
    render(
        summary,
        load_json(args.corrected_wp),
        load_json(args.historical_same_source_wp),
        png_path,
    )
    print(f"status={summary['status']}")
    print(f"json={json_path}")
    print(f"csv={csv_path}")
    print(f"png={png_path}")


if __name__ == "__main__":
    main()
