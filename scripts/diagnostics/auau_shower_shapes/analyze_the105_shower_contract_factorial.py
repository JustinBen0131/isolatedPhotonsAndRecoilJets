#!/usr/bin/env python3
"""Analyze the bounded THE-105 AuAu shower-contract factorial canary.

The three reconstruction arms are deliberately factorized:

  historical  : data TowerInfo+70 MeV; embedded RawCluster+70 MeV
  towerinfo70 : all populations TowerInfo+70 MeV
  canonical   : all populations TowerInfo+0 MeV

The script reads one diagnostic candidate row per reconstructed photon and
compares the same three selection stages: before preselection, after the full
NCB preselection, and after the frozen tight classifier.  It never promotes an
artifact or modifies a canonical current pointer.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import uproot


VARIANTS = ("historical", "towerinfo70", "canonical")
POPULATIONS = ("data", "signal", "inclusive")
STAGES = (
    ("before", "Before preselection"),
    ("preselection", "After NCB preselection"),
    ("tight", "After tight ID"),
)
CENTRALITIES = (
    ("cent0_20", "0--20%", 0.0, 20.0),
    ("cent20_50", "20--50%", 20.0, 50.0),
    ("cent50_80", "50--80%", 50.0, 80.0),
)


@dataclass(frozen=True)
class Variable:
    branch: str
    slug: str
    label: str
    xmin: float
    xmax: float
    bins: int
    primary: bool = False


VARIABLES = (
    Variable("cluster_weta_cogx", "weta_cogx", r"$w_{\eta}^{\mathrm{COGX}}$", 0.0, 1.25, 80, True),
    Variable("cluster_wphi_cogx", "wphi_cogx", r"$w_{\phi}^{\mathrm{COGX}}$", 0.0, 1.25, 80, True),
    Variable("e11_over_e33", "e11_over_e33", r"$E_{1\times1}/E_{3\times3}$", 0.0, 1.2, 60, True),
    Variable("e32_over_e35", "e32_over_e35", r"$E_{3\times2}/E_{3\times5}$", 0.0, 1.2, 60, True),
    Variable("cluster_et1", "et1", r"$e_{T,1}$", 0.0, 1.2, 60),
    Variable("cluster_et2", "et2", r"$e_{T,2}$", 0.0, 1.2, 60),
    Variable("cluster_et3", "et3", r"$e_{T,3}$", 0.0, 1.2, 60),
    Variable("cluster_et4", "et4", r"$e_{T,4}$", 0.0, 0.3, 60),
)

BASE_BRANCHES = (
    "run",
    "evt",
    "event_count",
    "is_sim",
    "is_sim_embedded",
    "sample_code",
    "cluster_Et",
    "cluster_Eta",
    "cluster_Phi",
    "centrality",
    "event_weight",
    "event_calo_total_energy",
    "reco_eiso",
    "preselection_pass",
    "baseline_bdt_tight",
    "active_tight_tag",
    "npb_pass",
    "auau_tight_bdt_score",
    "baseline_wp80_threshold",
    "cluster_truth_barcode",
)
BRANCHES = BASE_BRANCHES + tuple(v.branch for v in VARIABLES)

COLORS = {
    "historical": "#4D4D4D",
    "towerinfo70": "#D55E00",
    "canonical": "#0072B2",
}
LINESTYLES = {"historical": "-", "towerinfo70": "--", "canonical": "-"}
LABELS = {
    "historical": "Historical routing, 70 MeV",
    "towerinfo70": "TowerInfo, 70 MeV",
    "canonical": "TowerInfo, 0 MeV",
}
POP_LABELS = {
    "data": "Au+Au data",
    "signal": "Embedded photon simulation",
    "inclusive": "Embedded inclusive-jet simulation",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def config_text(root_file: uproot.ReadOnlyDirectory) -> str:
    if "analysis_config_yaml" not in root_file:
        return ""
    obj = root_file["analysis_config_yaml"]
    try:
        return str(obj)
    except Exception:
        try:
            return str(obj.member("fString"))
        except Exception:
            return ""


def stage_mask(frame: pd.DataFrame, stage: str) -> np.ndarray:
    if stage == "before":
        return np.ones(len(frame), dtype=bool)
    if stage == "preselection":
        return frame["preselection_pass"].to_numpy(dtype=np.int64) != 0
    if stage == "tight":
        return frame["baseline_bdt_tight"].to_numpy(dtype=np.int64) != 0
    raise KeyError(stage)


def cent_mask(frame: pd.DataFrame, lo: float, hi: float) -> np.ndarray:
    values = frame["centrality"].to_numpy(dtype=float)
    return np.isfinite(values) & (values >= lo) & (values < hi)


def event_weights(frame: pd.DataFrame, population: str) -> np.ndarray:
    if population == "data":
        return np.ones(len(frame), dtype=float)
    values = frame["event_weight"].to_numpy(dtype=float)
    return np.where(np.isfinite(values) & (values > 0.0), values, 0.0)


def candidate_keys(frame: pd.DataFrame) -> pd.DataFrame:
    stable = pd.DataFrame(
        {
            "sample_code": frame["sample_code"].astype("int64"),
            "run": frame["run"].astype("int64"),
            "evt": frame["evt"].astype("int64"),
            "cluster_Et": frame["cluster_Et"].round(5),
            "cluster_Eta": frame["cluster_Eta"].round(5),
            "cluster_Phi": frame["cluster_Phi"].round(5),
            "centrality": frame["centrality"].round(4),
            "event_calo_total_energy": frame["event_calo_total_energy"].round(3),
            "reco_eiso": frame["reco_eiso"].round(4),
            "cluster_truth_barcode": frame["cluster_truth_barcode"].astype("int64"),
        }
    )
    frame = frame.copy()
    frame["candidate_key"] = pd.util.hash_pandas_object(stable, index=False).astype("uint64")
    frame["candidate_dup"] = frame.groupby("candidate_key", sort=False).cumcount().astype("int32")
    return frame


def discover_and_load(
    input_root: Path,
    variant: str,
    population: str,
    input_rows: list[dict[str, object]],
    failures: list[str],
) -> pd.DataFrame:
    lane = input_root / variant / population
    files = sorted(lane.rglob("*.root")) if lane.exists() else []
    if not files:
        failures.append(f"missing ROOT files for {variant}/{population}: {lane}")
        return pd.DataFrame(columns=BRANCHES)

    frames: list[pd.DataFrame] = []
    for path in files:
        row: dict[str, object] = {
            "variant": variant,
            "population": population,
            "path": str(path.resolve()),
            "size": path.stat().st_size,
            "sha256": sha256(path),
            "tree_entries": 0,
            "variant_stamp_ok": False,
        }
        try:
            with uproot.open(path) as root_file:
                cfg = config_text(root_file)
                row["variant_stamp_ok"] = (
                    f"cemc_shower_shape_diagnostic_variant: {variant}" in cfg
                    or f"cemc_shower_shape_diagnostic_variant: '{variant}'" in cfg
                    or f'cemc_shower_shape_diagnostic_variant: "{variant}"' in cfg
                )
                if "AuAuPhotonCandidateSkim" not in root_file:
                    failures.append(f"missing AuAuPhotonCandidateSkim: {path}")
                    input_rows.append(row)
                    continue
                tree = root_file["AuAuPhotonCandidateSkim"]
                missing = [branch for branch in BRANCHES if branch not in tree.keys()]
                if missing:
                    failures.append(f"missing skim branches in {path}: {missing}")
                    input_rows.append(row)
                    continue
                arrays = tree.arrays(BRANCHES, library="np")
                frame = pd.DataFrame({name: arrays[name] for name in BRANCHES})
                row["tree_entries"] = len(frame)
                frame["source_root"] = str(path.resolve())
                frames.append(frame)
        except Exception as exc:
            failures.append(f"failed to read {path}: {exc}")
        input_rows.append(row)

    if not frames:
        return pd.DataFrame(columns=BRANCHES)
    frame = pd.concat(frames, ignore_index=True)
    base = (
        np.isfinite(frame["cluster_Et"].to_numpy(dtype=float))
        & (frame["cluster_Et"].to_numpy(dtype=float) >= 15.0)
        & (frame["cluster_Et"].to_numpy(dtype=float) < 35.0)
        & np.isfinite(frame["cluster_Eta"].to_numpy(dtype=float))
        & (np.abs(frame["cluster_Eta"].to_numpy(dtype=float)) < 0.7)
    )
    return candidate_keys(frame.loc[base].reset_index(drop=True))


def weighted_sum(mask: np.ndarray, weights: np.ndarray) -> float:
    return float(np.sum(weights[mask], dtype=np.float64))


def compute_metrics(
    frames: dict[tuple[str, str], pd.DataFrame],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    boundary_rows: list[dict[str, object]] = []
    acceptance_rows: list[dict[str, object]] = []
    for population in POPULATIONS:
        for cent_key, cent_label, cent_lo, cent_hi in CENTRALITIES:
            for variant in VARIANTS:
                frame = frames[(variant, population)]
                if frame.empty:
                    continue
                weights = event_weights(frame, population)
                cmask = cent_mask(frame, cent_lo, cent_hi)
                before_weight = weighted_sum(cmask, weights)
                for stage, stage_label in STAGES:
                    smask = cmask & stage_mask(frame, stage)
                    total_weight = weighted_sum(smask, weights)
                    acceptance_rows.append(
                        {
                            "population": population,
                            "centrality": cent_key,
                            "centrality_label": cent_label,
                            "variant": variant,
                            "stage": stage,
                            "candidate_rows": int(np.count_nonzero(smask)),
                            "weighted_candidates": total_weight,
                            "acceptance_from_before": total_weight / before_weight if before_weight > 0 else math.nan,
                        }
                    )
                    for variable in VARIABLES:
                        values = frame[variable.branch].to_numpy(dtype=float)
                        finite = np.isfinite(values)
                        zero = finite & np.isclose(values, 0.0, rtol=0.0, atol=1.0e-12)
                        under = finite & (values < variable.xmin)
                        over = finite & (values >= variable.xmax)
                        visible_nonzero = finite & (values > variable.xmin) & (values < variable.xmax)
                        finite_weight = weighted_sum(smask & finite, weights)
                        row = {
                            "population": population,
                            "centrality": cent_key,
                            "centrality_label": cent_label,
                            "variant": variant,
                            "stage": stage,
                            "variable": variable.branch,
                            "candidate_rows": int(np.count_nonzero(smask)),
                            "weighted_candidates": total_weight,
                            "finite_fraction": finite_weight / total_weight if total_weight > 0 else math.nan,
                            "nonfinite_fraction": weighted_sum(smask & ~finite, weights) / total_weight if total_weight > 0 else math.nan,
                            "exact_zero_fraction": weighted_sum(smask & zero, weights) / total_weight if total_weight > 0 else math.nan,
                            "underflow_fraction": weighted_sum(smask & under, weights) / total_weight if total_weight > 0 else math.nan,
                            "overflow_fraction": weighted_sum(smask & over, weights) / total_weight if total_weight > 0 else math.nan,
                            "visible_nonzero_fraction": weighted_sum(smask & visible_nonzero, weights) / total_weight if total_weight > 0 else math.nan,
                            "finite_mean": math.nan,
                            "finite_rms": math.nan,
                        }
                        fmask = smask & finite
                        fweight = weights[fmask]
                        if fweight.size and np.sum(fweight) > 0:
                            fvalue = values[fmask]
                            mean = float(np.average(fvalue, weights=fweight))
                            row["finite_mean"] = mean
                            row["finite_rms"] = float(np.sqrt(np.average((fvalue - mean) ** 2, weights=fweight)))
                        boundary_rows.append(row)
    return pd.DataFrame(boundary_rows), pd.DataFrame(acceptance_rows)


def normalized_histogram(
    frame: pd.DataFrame,
    population: str,
    variable: Variable,
    mask: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    values = frame[variable.branch].to_numpy(dtype=float)
    weights = event_weights(frame, population)
    denom = weighted_sum(mask, weights)
    bins = np.linspace(variable.xmin, variable.xmax, variable.bins + 1)
    continuous = mask & np.isfinite(values) & (values > variable.xmin) & (values < variable.xmax)
    hist, edges = np.histogram(values[continuous], bins=bins, weights=weights[continuous])
    if denom > 0:
        hist = hist.astype(float) / denom
    else:
        hist = hist.astype(float) * math.nan
    return hist, edges


def step_xy(hist: np.ndarray, edges: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    return edges, np.r_[hist, hist[-1] if hist.size else math.nan]


def render_overlay_matrix(
    output_dir: Path,
    frames: dict[tuple[str, str], pd.DataFrame],
    variable: Variable,
    cent_key: str,
    cent_label: str,
    cent_lo: float,
    cent_hi: float,
) -> Path:
    fig, axes = plt.subplots(3, 3, figsize=(15.6, 11.0), sharex=True, constrained_layout=False)
    fig.subplots_adjust(top=0.86, bottom=0.09, left=0.09, right=0.985, hspace=0.23, wspace=0.20)
    for row, population in enumerate(POPULATIONS):
        for col, (stage, stage_label) in enumerate(STAGES):
            ax = axes[row, col]
            for variant in VARIANTS:
                frame = frames[(variant, population)]
                mask = cent_mask(frame, cent_lo, cent_hi) & stage_mask(frame, stage)
                hist, edges = normalized_histogram(frame, population, variable, mask)
                x, y = step_xy(hist, edges)
                ax.step(
                    x,
                    y,
                    where="post",
                    color=COLORS[variant],
                    linestyle=LINESTYLES[variant],
                    linewidth=2.0 if variant != "historical" else 1.8,
                    label=LABELS[variant],
                )
            ax.set_xlim(variable.xmin, variable.xmax)
            ax.set_ylim(bottom=0.0)
            ax.grid(axis="y", color="#D7DEE8", linewidth=0.6, alpha=0.8)
            ax.tick_params(direction="in", top=True, right=True, labelsize=10)
            for spine in ax.spines.values():
                spine.set_linewidth(1.0)
            if row == 0:
                ax.set_title(stage_label, fontsize=13, fontweight="bold", color="#17365D", pad=8)
            if col == 0:
                ax.set_ylabel(f"{POP_LABELS[population]}\nweighted candidate fraction / bin", fontsize=10.5)
            if row == 2:
                ax.set_xlabel(variable.label, fontsize=12)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(0.59, 0.925), fontsize=11)
    fig.suptitle(
        f"Au+Au shower-contract factorial: {variable.label}, {cent_label}",
        x=0.09,
        y=0.975,
        ha="left",
        fontsize=19,
        fontweight="bold",
        color="#14213D",
    )
    fig.text(0.09, 0.935, r"$15\leq E_T^\gamma<35$ GeV, $|\eta^\gamma|<0.7$", fontsize=12, color="#314E6E")
    fig.text(0.09, 0.018, "sPHENIX Internal", fontsize=11.5, fontweight="bold")
    fig.text(
        0.985,
        0.018,
        "Continuous nonzero bins are divided by the full stage weight; zero, non-finite, underflow, and overflow fractions remain in the denominator and are tabulated separately.",
        fontsize=8.8,
        ha="right",
        color="#4B5563",
    )
    path = output_dir / f"the105_{variable.slug}_{cent_key}_stage_overlays.png"
    fig.savefig(path, dpi=220, facecolor="white")
    plt.close(fig)
    return path


def render_ratio_matrix(
    output_dir: Path,
    frames: dict[tuple[str, str], pd.DataFrame],
    variable: Variable,
    cent_key: str,
    cent_label: str,
    cent_lo: float,
    cent_hi: float,
) -> Path:
    fig, axes = plt.subplots(3, 3, figsize=(15.6, 10.5), sharex=True, constrained_layout=False)
    fig.subplots_adjust(top=0.84, bottom=0.09, left=0.09, right=0.985, hspace=0.22, wspace=0.20)
    contrasts = (
        ("towerinfo70", "historical", "TowerInfo70 / historical", "#D55E00"),
        ("canonical", "towerinfo70", "TowerInfo0 / TowerInfo70", "#0072B2"),
    )
    for row, population in enumerate(POPULATIONS):
        for col, (stage, stage_label) in enumerate(STAGES):
            ax = axes[row, col]
            cached: dict[str, tuple[np.ndarray, np.ndarray]] = {}
            for variant in VARIANTS:
                frame = frames[(variant, population)]
                mask = cent_mask(frame, cent_lo, cent_hi) & stage_mask(frame, stage)
                cached[variant] = normalized_histogram(frame, population, variable, mask)
            for numerator, denominator, label, color in contrasts:
                num, edges = cached[numerator]
                den, _ = cached[denominator]
                ratio = np.full_like(num, np.nan, dtype=float)
                np.divide(num, den, out=ratio, where=den > 0)
                centers = 0.5 * (edges[:-1] + edges[1:])
                ax.plot(centers, ratio, marker="o", markersize=2.4, linewidth=1.25, color=color, label=label)
            ax.axhline(1.0, color="#6B7280", linewidth=1.0, linestyle="--")
            ax.set_xlim(variable.xmin, variable.xmax)
            ax.set_ylim(0.45, 1.55)
            ax.grid(axis="y", color="#D7DEE8", linewidth=0.6, alpha=0.8)
            ax.tick_params(direction="in", top=True, right=True, labelsize=10)
            if row == 0:
                ax.set_title(stage_label, fontsize=13, fontweight="bold", color="#17365D", pad=8)
            if col == 0:
                ax.set_ylabel(f"{POP_LABELS[population]}\nfactorial ratio", fontsize=10.5)
            if row == 2:
                ax.set_xlabel(variable.label, fontsize=12)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.65, 0.91), fontsize=11)
    fig.suptitle(
        f"Factorized reconstruction effects: {variable.label}, {cent_label}",
        x=0.09,
        y=0.97,
        ha="left",
        fontsize=19,
        fontweight="bold",
        color="#14213D",
    )
    fig.text(0.09, 0.925, "Source effect at fixed 70 MeV; threshold effect at fixed TowerInfo source", fontsize=12, color="#314E6E")
    fig.text(0.09, 0.018, "sPHENIX Internal", fontsize=11.5, fontweight="bold")
    path = output_dir / f"the105_{variable.slug}_{cent_key}_factorial_ratios.png"
    fig.savefig(path, dpi=220, facecolor="white")
    plt.close(fig)
    return path


def compute_migrations(
    frames: dict[tuple[str, str], pd.DataFrame],
    failures: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[dict[str, object]] = []
    delta_rows: list[dict[str, object]] = []
    comparisons = (
        ("historical", "towerinfo70", "source_at_70mev"),
        ("towerinfo70", "canonical", "floor_at_towerinfo"),
        ("historical", "canonical", "combined_historical_to_canonical"),
    )
    for population in POPULATIONS:
        key_cols = ["candidate_key", "candidate_dup"]
        keep = key_cols + [
            "centrality",
            "event_weight",
            "preselection_pass",
            "baseline_bdt_tight",
            "auau_tight_bdt_score",
        ] + [v.branch for v in VARIABLES]
        for source_variant, target_variant, comparison in comparisons:
            source = frames[(source_variant, population)][keep].copy()
            target = frames[(target_variant, population)][keep].copy()
            merged = source.merge(
                target,
                on=key_cols,
                how="inner",
                suffixes=(f"_{source_variant}", f"_{target_variant}"),
            )
            overlap = len(merged) / max(len(source), len(target), 1)
            if overlap < 0.995:
                failures.append(
                    f"candidate overlap below 99.5% for {population}/{comparison}: {overlap:.6f}"
                )
            centrality_column = f"centrality_{target_variant}"
            weight_column = f"event_weight_{target_variant}"
            for cent_key, cent_label, cent_lo, cent_hi in CENTRALITIES:
                centrality = merged[centrality_column].to_numpy(dtype=float)
                cmask = np.isfinite(centrality) & (centrality >= cent_lo) & (centrality < cent_hi)
                subset = merged.loc[cmask]
                weights = (
                    np.ones(len(subset), dtype=float)
                    if population == "data"
                    else np.where(
                        np.isfinite(subset[weight_column].to_numpy(dtype=float))
                        & (subset[weight_column].to_numpy(dtype=float) > 0),
                        subset[weight_column].to_numpy(dtype=float),
                        0.0,
                    )
                )
                for stage, branch in (("preselection", "preselection_pass"), ("tight", "baseline_bdt_tight")):
                    source_pass = subset[f"{branch}_{source_variant}"].to_numpy(dtype=int) != 0
                    target_pass = subset[f"{branch}_{target_variant}"].to_numpy(dtype=int) != 0
                    for transition, mask in (
                        ("fail_to_fail", ~source_pass & ~target_pass),
                        ("fail_to_pass", ~source_pass & target_pass),
                        ("pass_to_fail", source_pass & ~target_pass),
                        ("pass_to_pass", source_pass & target_pass),
                    ):
                        rows.append(
                            {
                                "population": population,
                                "centrality": cent_key,
                                "centrality_label": cent_label,
                                "comparison": comparison,
                                "source_variant": source_variant,
                                "target_variant": target_variant,
                                "stage": stage,
                                "transition_source_to_target": transition,
                                "candidate_rows": int(np.count_nonzero(mask)),
                                "weighted_candidates": float(np.sum(weights[mask])),
                                "matched_rows": len(subset),
                                "global_overlap_fraction": overlap,
                            }
                        )
                for variable in VARIABLES:
                    source_values = subset[f"{variable.branch}_{source_variant}"].to_numpy(dtype=float)
                    target_values = subset[f"{variable.branch}_{target_variant}"].to_numpy(dtype=float)
                    finite = np.isfinite(source_values) & np.isfinite(target_values)
                    delta = target_values[finite] - source_values[finite]
                    delta_rows.append(
                        {
                            "population": population,
                            "centrality": cent_key,
                            "comparison": comparison,
                            "source_variant": source_variant,
                            "target_variant": target_variant,
                            "variable": variable.branch,
                            "matched_finite_rows": int(delta.size),
                            "mean_target_minus_source": float(np.mean(delta)) if delta.size else math.nan,
                            "rms_target_minus_source": float(np.sqrt(np.mean(delta**2))) if delta.size else math.nan,
                            "median_target_minus_source": float(np.median(delta)) if delta.size else math.nan,
                            "p95_abs_delta": float(np.quantile(np.abs(delta), 0.95)) if delta.size else math.nan,
                            "max_abs_delta": float(np.max(np.abs(delta))) if delta.size else math.nan,
                        }
                    )

        # Historical and TowerInfo70 must be identical in data by construction.
        if population == "data":
            left = frames[("historical", "data")]
            right = frames[("towerinfo70", "data")]
            check = left.merge(right, on=key_cols, how="inner", suffixes=("_historical", "_towerinfo70"))
            checksum_overlap = len(check) / max(len(left), len(right), 1)
            if checksum_overlap < 0.999999:
                failures.append(f"data historical/TowerInfo70 checksum overlap is {checksum_overlap:.8f}")
            for variable in VARIABLES:
                a = check[f"{variable.branch}_historical"].to_numpy(dtype=float)
                b = check[f"{variable.branch}_towerinfo70"].to_numpy(dtype=float)
                finite = np.isfinite(a) & np.isfinite(b)
                if np.any(np.abs(a[finite] - b[finite]) > 1.0e-7):
                    failures.append(f"data source-control checksum changed {variable.branch}")
            for branch in ("preselection_pass", "baseline_bdt_tight"):
                if not np.array_equal(
                    check[f"{branch}_historical"].to_numpy(),
                    check[f"{branch}_towerinfo70"].to_numpy(),
                ):
                    failures.append(f"data source-control checksum changed {branch}")
    return pd.DataFrame(rows), pd.DataFrame(delta_rows)


def render_acceptance(output_dir: Path, acceptance: pd.DataFrame) -> Path:
    fig, axes = plt.subplots(3, 3, figsize=(14.8, 10.2), sharey=True)
    fig.subplots_adjust(top=0.84, bottom=0.09, left=0.09, right=0.985, hspace=0.26, wspace=0.15)
    x = np.arange(2, dtype=float)
    width = 0.23
    for row, population in enumerate(POPULATIONS):
        for col, (cent_key, cent_label, _, _) in enumerate(CENTRALITIES):
            ax = axes[row, col]
            for offset, variant in enumerate(VARIANTS):
                values = []
                for stage in ("preselection", "tight"):
                    query = acceptance[
                        (acceptance.population == population)
                        & (acceptance.centrality == cent_key)
                        & (acceptance.variant == variant)
                        & (acceptance.stage == stage)
                    ]
                    values.append(float(query.acceptance_from_before.iloc[0]) if len(query) else math.nan)
                ax.bar(x + (offset - 1) * width, values, width=width, color=COLORS[variant], label=LABELS[variant], alpha=0.9)
            ax.set_xticks(x, ["NCB preselection", "tight ID"], fontsize=9)
            ax.set_ylim(0.0, 1.0)
            ax.grid(axis="y", color="#D7DEE8", linewidth=0.6, alpha=0.8)
            ax.tick_params(direction="in", top=True, right=True)
            if row == 0:
                ax.set_title(cent_label, fontsize=13, fontweight="bold", color="#17365D")
            if col == 0:
                ax.set_ylabel(f"{POP_LABELS[population]}\nweighted acceptance", fontsize=10.5)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(0.62, 0.91), fontsize=10.5)
    fig.suptitle("Selection migration under the shower-shape reconstruction contract", x=0.09, y=0.97, ha="left", fontsize=18, fontweight="bold", color="#14213D")
    fig.text(0.09, 0.925, r"$15\leq E_T^\gamma<35$ GeV, $|\eta^\gamma|<0.7$; acceptances are relative to the same before-preselection candidate population", fontsize=11.5, color="#314E6E")
    fig.text(0.09, 0.018, "sPHENIX Internal", fontsize=11.5, fontweight="bold")
    path = output_dir / "the105_selection_acceptance_by_variant.png"
    fig.savefig(path, dpi=220, facecolor="white")
    plt.close(fig)
    return path


def write_csv(path: Path, frame: pd.DataFrame) -> None:
    frame.to_csv(path, index=False, quoting=csv.QUOTE_MINIMAL)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--all-variables", action="store_true", help="Render all eight baseV3E shower variables; default renders the four routing-sensitive variables.")
    parser.add_argument("--strict", action="store_true", help="Return nonzero if any provenance, candidate-overlap, or data-control checksum fails.")
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    failures: list[str] = []
    input_rows: list[dict[str, object]] = []
    frames: dict[tuple[str, str], pd.DataFrame] = {}
    for variant in VARIANTS:
        for population in POPULATIONS:
            frames[(variant, population)] = discover_and_load(
                args.input_root, variant, population, input_rows, failures
            )

    for row in input_rows:
        if not row.get("variant_stamp_ok", False):
            failures.append(f"missing or incorrect variant stamp: {row['path']}")

    boundary, acceptance = compute_metrics(frames)
    migrations, deltas = compute_migrations(frames, failures)
    write_csv(args.output_dir / "input_roots.csv", pd.DataFrame(input_rows))
    write_csv(args.output_dir / "boundary_and_range_fractions.csv", boundary)
    write_csv(args.output_dir / "selection_acceptance.csv", acceptance)
    write_csv(args.output_dir / "selection_migration.csv", migrations)
    write_csv(args.output_dir / "matched_feature_deltas.csv", deltas)

    plots: list[str] = []
    selected = VARIABLES if args.all_variables else tuple(v for v in VARIABLES if v.primary)
    if all(not frame.empty for frame in frames.values()):
        for variable in selected:
            for cent_key, cent_label, cent_lo, cent_hi in CENTRALITIES:
                plots.append(str(render_overlay_matrix(args.output_dir, frames, variable, cent_key, cent_label, cent_lo, cent_hi)))
                plots.append(str(render_ratio_matrix(args.output_dir, frames, variable, cent_key, cent_label, cent_lo, cent_hi)))
        plots.append(str(render_acceptance(args.output_dir, acceptance)))

    summary = {
        "status": "pass" if not failures else "fail",
        "input_root": str(args.input_root.resolve()),
        "output_dir": str(args.output_dir.resolve()),
        "variants": list(VARIANTS),
        "populations": list(POPULATIONS),
        "stages": [stage for stage, _ in STAGES],
        "centralities": [key for key, _, _, _ in CENTRALITIES],
        "rendered_variables": [variable.branch for variable in selected],
        "failures": failures,
        "plots": plots,
        "interpretation_contract": {
            "towerinfo70_over_historical": "energy-source effect at fixed 70 MeV; data is an exact control checksum",
            "canonical_over_towerinfo70": "tower-floor effect at fixed complete good-TowerInfo source",
            "after_tight": "frozen old-model diagnostic; not retrained-model performance",
            "normalization": "continuous nonzero visible bins divided by the full stage candidate weight",
        },
    }
    (args.output_dir / "validation_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps({"status": summary["status"], "failures": failures, "plots": len(plots)}, indent=2))
    return 2 if failures and args.strict else 0


if __name__ == "__main__":
    raise SystemExit(main())
