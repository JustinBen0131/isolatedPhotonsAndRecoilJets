#!/usr/bin/env python3
"""Build the final THE-100 bounded-versus-complement comparison artifacts."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-the100")

import numpy as np
import uproot


CENTRALITIES = ("0_20", "20_50", "50_80")
CENTRALITY_LABELS = {"0_20": "0-20%", "20_50": "20-50%", "50_80": "50-80%"}
ISOLATIONS = (
    "isoR30_fixedIso4GeV",
    "isoR40_fixedIso4GeV",
    "isoR30_isSliding",
    "isoR40_isSliding",
)
CANONICAL_ISOLATION = "isoR40_isSliding"
PT_BINS = ((15, 17), (17, 19), (19, 21), (21, 23), (23, 26), (26, 35))
PT_MIN = 15.0
PT_MAX = 35.0
DATA_TOPDIR = "photon_12_plus_MBD_NS_geq_2_vtx_lt_150"
SIM_TOPDIR = "SIM"
XJ_KEY = "r04"
VIEWS = ("bounded", "complement")
ROLES = ("data", "signal", "inclusive")


def finite_or_none(value: float | None) -> float | None:
    if value is None or not math.isfinite(value):
        return None
    return float(value)


def ratio_with_error(
    numerator: float,
    numerator_variance: float,
    denominator: float,
    denominator_variance: float,
) -> tuple[float | None, float | None]:
    if denominator == 0.0:
        return None, None
    value = numerator / denominator
    variance = numerator_variance / denominator**2
    variance += numerator**2 * denominator_variance / denominator**4
    return finite_or_none(value), finite_or_none(math.sqrt(max(variance, 0.0)))


def hist_sum_1d(hist: Any, low: float, high: float) -> tuple[float, float]:
    values = np.asarray(hist.values(flow=False), dtype=np.float64)
    variances = hist.variances(flow=False)
    if variances is None:
        raise RuntimeError(f"histogram has no variances: {hist.name}")
    variances = np.asarray(variances, dtype=np.float64)
    edges = np.asarray(hist.axis().edges(), dtype=np.float64)
    mask = (edges[:-1] >= low - 1e-9) & (edges[1:] <= high + 1e-9)
    if not np.any(mask):
        raise RuntimeError(f"no bins in [{low}, {high}] for {hist.name}")
    if not np.all(np.isfinite(values[mask])) or not np.all(np.isfinite(variances[mask])):
        raise RuntimeError(f"non-finite bins in {hist.name}")
    return float(np.sum(values[mask])), float(np.sum(variances[mask]))


def hist_sum_2d_pt_window(hist: Any, low: float, high: float) -> tuple[float, float]:
    values = np.asarray(hist.values(flow=False), dtype=np.float64)
    variances = hist.variances(flow=False)
    if variances is None:
        raise RuntimeError(f"histogram has no variances: {hist.name}")
    variances = np.asarray(variances, dtype=np.float64)
    edges = np.asarray(hist.axis(0).edges(), dtype=np.float64)
    mask = (edges[:-1] >= low - 1e-9) & (edges[1:] <= high + 1e-9)
    if not np.any(mask):
        raise RuntimeError(f"no pT bins in [{low}, {high}] for {hist.name}")
    selected_values = values[mask, ...]
    selected_variances = variances[mask, ...]
    if not np.all(np.isfinite(selected_values)) or not np.all(np.isfinite(selected_variances)):
        raise RuntimeError(f"non-finite bins in {hist.name}")
    return float(np.sum(selected_values)), float(np.sum(selected_variances))


def abcd_metrics(counts: dict[str, float], variances: dict[str, float]) -> dict[str, Any]:
    a, b, c, d = (counts[region] for region in "ABCD")
    va, vb, vc, vd = (variances[region] for region in "ABCD")

    ad_over_bc = None
    ad_over_bc_error = None
    if b != 0.0 and c != 0.0:
        ad_over_bc = a * d / (b * c)
        relative_variance = 0.0
        for value, variance in ((a, va), (b, vb), (c, vc), (d, vd)):
            if value != 0.0:
                relative_variance += variance / value**2
        ad_over_bc_error = abs(ad_over_bc) * math.sqrt(max(relative_variance, 0.0))

    background = None
    background_variance = None
    if d != 0.0:
        background = b * c / d
        background_variance = (
            (c / d) ** 2 * vb
            + (b / d) ** 2 * vc
            + (b * c / d**2) ** 2 * vd
        )

    subtracted = None if background is None else a - background
    subtracted_fraction = None if subtracted is None or a == 0.0 else subtracted / a
    subtracted_variance = None if background_variance is None else va + background_variance
    relative_precision = None
    if subtracted is not None and subtracted != 0.0 and subtracted_variance is not None:
        relative_precision = math.sqrt(max(subtracted_variance, 0.0)) / abs(subtracted)

    return {
        "ad_over_bc": finite_or_none(ad_over_bc),
        "ad_over_bc_error": finite_or_none(ad_over_bc_error),
        "abcd_background_A": finite_or_none(background),
        "abcd_background_A_error": finite_or_none(
            None if background_variance is None else math.sqrt(max(background_variance, 0.0))
        ),
        "subtracted_A": finite_or_none(subtracted),
        "subtracted_fraction": finite_or_none(subtracted_fraction),
        "subtracted_relative_stat_uncertainty": finite_or_none(relative_precision),
    }


def prompt_leakage(
    root_file: uproot.ReadOnlyDirectory,
    topdir: str,
    isolation: str,
    centrality: str,
) -> dict[str, Any]:
    counts = np.zeros(4, dtype=np.float64)
    variances = np.zeros(4, dtype=np.float64)
    keys: list[str] = []
    for low, high in PT_BINS:
        key = (
            f"{topdir}/h_xJpurityLead_sigABCD_MC_{isolation}_"
            f"pT_{low}_{high}_cent_{centrality}"
        )
        hist = root_file[key]
        values = np.asarray(hist.values(flow=False), dtype=np.float64)
        hist_variances = hist.variances(flow=False)
        if values.shape != (4,) or hist_variances is None:
            raise RuntimeError(f"unexpected prompt-leakage histogram: {key}")
        counts += values
        variances += np.asarray(hist_variances, dtype=np.float64)
        keys.append(key)

    total = float(np.sum(counts))
    total_variance = float(np.sum(variances))
    leaked = float(np.sum(counts[1:]))
    leaked_variance = float(np.sum(variances[1:]))
    leakage, leakage_error = ratio_with_error(leaked, leaked_variance, total, total_variance)
    return {
        "prompt_A": float(counts[0]),
        "prompt_B": float(counts[1]),
        "prompt_C": float(counts[2]),
        "prompt_D": float(counts[3]),
        "prompt_leakage_fraction": leakage,
        "prompt_leakage_fraction_error": leakage_error,
        "prompt_leakage_keys": keys,
    }


def collect_sample_row(
    root_file: uproot.ReadOnlyDirectory,
    role: str,
    view: str,
    topdir: str,
    isolation: str,
    centrality: str,
) -> dict[str, Any]:
    counts: dict[str, float] = {}
    variances: dict[str, float] = {}
    abcd_keys: dict[str, str] = {}
    for region in "ABCD":
        key = f"{topdir}/h_pTgamma_ABCD_{region}_{isolation}_cent_{centrality}"
        counts[region], variances[region] = hist_sum_1d(root_file[key], PT_MIN, PT_MAX)
        abcd_keys[region] = key

    xj_a_key = f"{topdir}/h2_unfoldReco_pTgamma_xJ_incl_{XJ_KEY}_{isolation}_cent_{centrality}"
    xj_c_key = (
        f"{topdir}/h2_unfoldReco_pTgamma_xJ_incl_sidebandC_"
        f"{XJ_KEY}_{isolation}_cent_{centrality}"
    )
    xj_a, xj_a_variance = hist_sum_2d_pt_window(root_file[xj_a_key], PT_MIN, PT_MAX)
    xj_c, xj_c_variance = hist_sum_2d_pt_window(root_file[xj_c_key], PT_MIN, PT_MAX)
    xj_c_over_a, xj_c_over_a_error = ratio_with_error(
        xj_c, xj_c_variance, xj_a, xj_a_variance
    )

    row: dict[str, Any] = {
        "row_type": "sample",
        "role": role,
        "view": view,
        "isolation": isolation,
        "centrality": centrality,
        "pt_min_GeV": PT_MIN,
        "pt_max_GeV": PT_MAX,
        **counts,
        **{f"{region}_error": math.sqrt(max(variances[region], 0.0)) for region in "ABCD"},
        **abcd_metrics(counts, variances),
        "inclusive_background_closure_ratio": None,
        "inclusive_background_closure_residual": None,
        "prompt_A": None,
        "prompt_B": None,
        "prompt_C": None,
        "prompt_D": None,
        "prompt_leakage_fraction": None,
        "prompt_leakage_fraction_error": None,
        "xj_region_A_yield": xj_a,
        "xj_region_A_error": math.sqrt(max(xj_a_variance, 0.0)),
        "xj_region_C_yield": xj_c,
        "xj_region_C_error": math.sqrt(max(xj_c_variance, 0.0)),
        "xj_region_C_over_A": xj_c_over_a,
        "xj_region_C_over_A_error": xj_c_over_a_error,
        "xjgamma_structural_ready": bool(xj_a > 0.0 and xj_c > 0.0),
        "subtraction_closure_ratio": None,
        "subtraction_closure_residual_fraction": None,
        "subtraction_closure_residual_uncertainty": None,
        "source_keys": {"abcd": abcd_keys, "xj_region_A": xj_a_key, "xj_region_C": xj_c_key},
    }
    if role == "inclusive" and row["abcd_background_A"] is not None and counts["A"] != 0.0:
        row["inclusive_background_closure_ratio"] = row["abcd_background_A"] / counts["A"]
        row["inclusive_background_closure_residual"] = 1.0 - row["inclusive_background_closure_ratio"]
    if role == "signal":
        leakage = prompt_leakage(root_file, topdir, isolation, centrality)
        row.update({key: value for key, value in leakage.items() if key != "prompt_leakage_keys"})
        row["source_keys"]["prompt_leakage"] = leakage["prompt_leakage_keys"]
    row["_counts"] = counts
    row["_variances"] = variances
    return row


def closure_row(signal: dict[str, Any], inclusive: dict[str, Any]) -> dict[str, Any]:
    signal_counts = signal["_counts"]
    inclusive_counts = inclusive["_counts"]
    signal_variances = signal["_variances"]
    inclusive_variances = inclusive["_variances"]
    combined_counts = {
        region: signal_counts[region] + inclusive_counts[region] for region in "ABCD"
    }
    combined_variances = {
        region: signal_variances[region] + inclusive_variances[region] for region in "ABCD"
    }
    metrics = abcd_metrics(combined_counts, combined_variances)
    predicted_signal = metrics["subtracted_A"]
    known_signal = signal_counts["A"]
    closure_ratio = None if predicted_signal is None or known_signal == 0.0 else predicted_signal / known_signal
    residual = None if closure_ratio is None else closure_ratio - 1.0

    predicted_background = metrics["abcd_background_A"]
    predicted_background_error = metrics["abcd_background_A_error"]
    residual_uncertainty = None
    if predicted_background_error is not None and known_signal != 0.0:
        residual_uncertainty = math.sqrt(
            inclusive_variances["A"] + predicted_background_error**2
        ) / abs(known_signal)

    return {
        "row_type": "closure",
        "role": "signal_plus_inclusive",
        "view": signal["view"],
        "isolation": signal["isolation"],
        "centrality": signal["centrality"],
        "pt_min_GeV": PT_MIN,
        "pt_max_GeV": PT_MAX,
        **combined_counts,
        **{f"{region}_error": math.sqrt(max(combined_variances[region], 0.0)) for region in "ABCD"},
        **metrics,
        "inclusive_background_closure_ratio": None,
        "inclusive_background_closure_residual": None,
        "prompt_A": None,
        "prompt_B": None,
        "prompt_C": None,
        "prompt_D": None,
        "prompt_leakage_fraction": None,
        "prompt_leakage_fraction_error": None,
        "xj_region_A_yield": None,
        "xj_region_A_error": None,
        "xj_region_C_yield": None,
        "xj_region_C_error": None,
        "xj_region_C_over_A": None,
        "xj_region_C_over_A_error": None,
        "xjgamma_structural_ready": None,
        "subtraction_closure_ratio": finite_or_none(closure_ratio),
        "subtraction_closure_residual_fraction": finite_or_none(residual),
        "subtraction_closure_residual_uncertainty": finite_or_none(residual_uncertainty),
        "source_keys": {
            "signal": signal["source_keys"]["abcd"],
            "inclusive": inclusive["source_keys"]["abcd"],
            "definition": "(signal+inclusive A - (signal+inclusive B)(signal+inclusive C)/(signal+inclusive D))/signal A",
        },
    }


CSV_FIELDS = (
    "row_type", "role", "view", "isolation", "centrality", "pt_min_GeV", "pt_max_GeV",
    "A", "A_error", "B", "B_error", "C", "C_error", "D", "D_error",
    "ad_over_bc", "ad_over_bc_error", "abcd_background_A", "abcd_background_A_error",
    "subtracted_A", "subtracted_fraction", "subtracted_relative_stat_uncertainty",
    "inclusive_background_closure_ratio", "inclusive_background_closure_residual",
    "prompt_A", "prompt_B", "prompt_C", "prompt_D", "prompt_leakage_fraction",
    "prompt_leakage_fraction_error", "xj_region_A_yield", "xj_region_A_error",
    "xj_region_C_yield", "xj_region_C_error", "xj_region_C_over_A",
    "xj_region_C_over_A_error", "xjgamma_structural_ready", "subtraction_closure_ratio",
    "subtraction_closure_residual_fraction", "subtraction_closure_residual_uncertainty",
)


def clean_row(row: dict[str, Any]) -> dict[str, Any]:
    return {key: value for key, value in row.items() if not key.startswith("_")}


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_FIELDS)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key) for key in CSV_FIELDS})


def canonical_rows(rows: list[dict[str, Any]], row_type: str, role: str) -> list[dict[str, Any]]:
    return [
        row for row in rows
        if row["row_type"] == row_type
        and row["role"] == role
        and row["isolation"] == CANONICAL_ISOLATION
    ]


def plot_comparison(path: Path, rows: list[dict[str, Any]]) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    colors = {"bounded": "#007C91", "complement": "#D55E00"}
    region_colors = {"A": "#0072B2", "B": "#E69F00", "C": "#009E73", "D": "#777777"}
    centers = np.arange(len(CENTRALITIES), dtype=float)
    width = 0.34

    fig, axes = plt.subplots(3, 2, figsize=(13.2, 12.0))
    fig.patch.set_facecolor("white")
    for ax in axes.flat:
        ax.set_facecolor("white")
        ax.grid(axis="y", color="#d9d9d9", linewidth=0.8, alpha=0.7)
        ax.set_axisbelow(True)
        ax.tick_params(direction="in", top=True, right=True)

    data_rows = canonical_rows(rows, "sample", "data")
    signal_rows = canonical_rows(rows, "sample", "signal")
    closure_rows = canonical_rows(rows, "closure", "signal_plus_inclusive")

    ax = axes[0, 0]
    for view_index, view in enumerate(VIEWS):
        bottom = np.zeros(len(CENTRALITIES), dtype=float)
        subset = {row["centrality"]: row for row in data_rows if row["view"] == view}
        totals = np.asarray([sum(subset[cent][region] for region in "ABCD") for cent in CENTRALITIES])
        x = centers + (view_index - 0.5) * width
        for region in "ABCD":
            values = np.asarray([subset[cent][region] for cent in CENTRALITIES]) / totals
            ax.bar(
                x,
                values,
                width,
                bottom=bottom,
                color=region_colors[region],
                edgecolor="#333333" if view == "complement" else "white",
                linewidth=0.55,
                hatch="//" if view == "complement" else None,
            )
            bottom += values
    ax.set_title("Data ABCD composition")
    ax.set_ylabel("Fraction of A+B+C+D")
    ax.set_ylim(0.0, 1.0)
    region_legend = ax.legend(
        handles=[Patch(facecolor=region_colors[region], label=region) for region in "ABCD"],
        ncol=4,
        fontsize=9,
        loc="upper center",
    )
    ax.add_artist(region_legend)
    ax.legend(
        handles=[
            Patch(facecolor="white", edgecolor="#333333", label="Bounded"),
            Patch(facecolor="white", edgecolor="#333333", hatch="//", label="Complement"),
        ],
        ncol=2,
        fontsize=8,
        loc="lower center",
    )

    def grouped_bars(ax: Any, source_rows: list[dict[str, Any]], metric: str, ylabel: str, title: str, error: str | None = None) -> None:
        for view_index, view in enumerate(VIEWS):
            subset = {row["centrality"]: row for row in source_rows if row["view"] == view}
            values = [subset[cent][metric] for cent in CENTRALITIES]
            errors = None if error is None else [subset[cent][error] for cent in CENTRALITIES]
            ax.bar(
                centers + (view_index - 0.5) * width,
                values,
                width,
                yerr=errors,
                capsize=3,
                color=colors[view],
                label=view.capitalize(),
            )
        ax.set_title(title)
        ax.set_ylabel(ylabel)
        ax.legend(fontsize=9)

    grouped_bars(axes[0, 1], data_rows, "ad_over_bc", "AD / BC", "Data ABCD factorization", "ad_over_bc_error")
    axes[0, 1].axhline(1.0, color="#222222", linestyle="--", linewidth=1.1)
    grouped_bars(
        axes[1, 0], signal_rows, "prompt_leakage_fraction", "Leakage fraction",
        "Truth-signal leakage outside A", "prompt_leakage_fraction_error",
    )
    grouped_bars(
        axes[1, 1], data_rows, "subtracted_relative_stat_uncertainty", "Relative statistical uncertainty",
        "Data ABCD subtraction precision",
    )
    grouped_bars(
        axes[2, 0], closure_rows, "subtraction_closure_residual_fraction", "(predicted - true signal) / true signal",
        "Signal+inclusive subtraction closure", "subtraction_closure_residual_uncertainty",
    )
    axes[2, 0].axhline(0.0, color="#222222", linestyle="--", linewidth=1.1)
    grouped_bars(
        axes[2, 1], data_rows, "xj_region_C_over_A", "Region C / Region A recoil yield",
        "xJgamma structural readiness", "xj_region_C_over_A_error",
    )
    ready_count = sum(bool(row["xjgamma_structural_ready"]) for row in data_rows)
    axes[2, 1].text(
        0.98, 0.95, f"finite, populated A/C: {ready_count}/{len(data_rows)}",
        transform=axes[2, 1].transAxes, ha="right", va="top", fontsize=9,
        bbox={"facecolor": "white", "edgecolor": "#aaaaaa", "pad": 4},
    )

    for ax in axes.flat:
        ax.set_xticks(centers, [CENTRALITY_LABELS[cent] for cent in CENTRALITIES])
    fig.text(0.055, 0.985, r"$\bf{\it{sPHENIX}}$ Internal", ha="left", va="top", fontsize=15)
    fig.suptitle(
        "THE-100 Au+Au bounded sideband vs unrestricted complement\n"
        r"$15<E_T^\gamma<35$ GeV, $R=0.4$ sliding isolation, photon-12 data stream",
        fontsize=15,
        y=0.985,
    )
    fig.tight_layout(rect=(0.03, 0.03, 0.99, 0.94), h_pad=2.0, w_pad=1.5)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    for role in ROLES:
        for view in VIEWS:
            parser.add_argument(f"--{role}-{view}", type=Path)
    parser.add_argument("--qa-report", action="append", type=Path, default=[])
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--prefix", default="the100_dualview_final_comparison")
    parser.add_argument("--skip-png", action="store_true")
    parser.add_argument(
        "--plot-json",
        type=Path,
        help="Render the PNG from an existing comparison JSON without reopening ROOT files.",
    )
    args = parser.parse_args()

    if args.plot_json:
        report = json.loads(args.plot_json.read_text())
        png_path = args.output_dir / f"{args.prefix}.png"
        plot_comparison(png_path, report["rows"])
        report["outputs"]["png"] = str(png_path)
        report["plot_rendered_from_json"] = True
        args.plot_json.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n")
        print(json.dumps({"status": report["status"], "png": str(png_path)}, indent=2))
        return 0 if report["status"] == "PASS" else 1

    qa_reports: list[dict[str, Any]] = []
    failures: list[str] = []
    for path in args.qa_report:
        payload = json.loads(path.read_text())
        qa_reports.append({"path": str(path), "status": payload.get("status"), "failures": payload.get("failures", [])})
        if payload.get("status") != "PASS" or payload.get("failures"):
            failures.append(f"strict QA is not PASS: {path}")

    missing_args = [
        f"--{role}-{view}"
        for role in ROLES for view in VIEWS
        if getattr(args, f"{role}_{view}") is None
    ]
    if missing_args:
        parser.error("ROOT extraction requires " + ", ".join(missing_args))
    sources = {
        (role, view): Path(getattr(args, f"{role}_{view}"))
        for role in ROLES for view in VIEWS
    }
    root_files: dict[tuple[str, str], uproot.ReadOnlyDirectory] = {}
    rows: list[dict[str, Any]] = []
    try:
        for key, path in sources.items():
            if not path.is_file() or path.stat().st_size < 50_000:
                raise RuntimeError(f"missing or tiny final ROOT: {path}")
            root_files[key] = uproot.open(path)

        sample_index: dict[tuple[str, str, str, str], dict[str, Any]] = {}
        for role in ROLES:
            topdir = DATA_TOPDIR if role == "data" else SIM_TOPDIR
            for view in VIEWS:
                root_file = root_files[(role, view)]
                for isolation in ISOLATIONS:
                    for centrality in CENTRALITIES:
                        row = collect_sample_row(root_file, role, view, topdir, isolation, centrality)
                        rows.append(row)
                        sample_index[(role, view, isolation, centrality)] = row

        for view in VIEWS:
            for isolation in ISOLATIONS:
                for centrality in CENTRALITIES:
                    rows.append(
                        closure_row(
                            sample_index[("signal", view, isolation, centrality)],
                            sample_index[("inclusive", view, isolation, centrality)],
                        )
                    )

        for role in ROLES:
            for isolation in ISOLATIONS:
                for centrality in CENTRALITIES:
                    bounded = sample_index[(role, "bounded", isolation, centrality)]
                    complement = sample_index[(role, "complement", isolation, centrality)]
                    for region in "AB":
                        if not math.isclose(bounded[region], complement[region], rel_tol=0.0, abs_tol=1e-6):
                            failures.append(f"tight region {region} differs: {role} {isolation} {centrality}")
                    for region in "CD":
                        if bounded[region] > complement[region] + 1e-6:
                            failures.append(f"bounded {region} exceeds complement: {role} {isolation} {centrality}")
                    if bounded["xj_region_A_yield"] != complement["xj_region_A_yield"]:
                        failures.append(f"xJ region A differs: {role} {isolation} {centrality}")
                    if bounded["xj_region_C_yield"] > complement["xj_region_C_yield"] + 1e-6:
                        failures.append(f"bounded xJ region C exceeds complement: {role} {isolation} {centrality}")

        expected_rows = len(ROLES) * len(VIEWS) * len(ISOLATIONS) * len(CENTRALITIES)
        expected_rows += len(VIEWS) * len(ISOLATIONS) * len(CENTRALITIES)
        if len(rows) != expected_rows:
            failures.append(f"wrong row count: {len(rows)} != {expected_rows}")
        for row in rows:
            for key in ("A", "B", "C", "D", "ad_over_bc", "abcd_background_A"):
                if row.get(key) is None or not math.isfinite(float(row[key])):
                    failures.append(f"non-finite {key}: {row['row_type']} {row['role']} {row['view']} {row['isolation']} {row['centrality']}")
        data_ready = [
            row for row in rows
            if row["row_type"] == "sample" and row["role"] == "data"
        ]
        if not all(row["xjgamma_structural_ready"] for row in data_ready):
            failures.append("one or more data xJgamma A/C surfaces are empty")

        args.output_dir.mkdir(parents=True, exist_ok=True)
        csv_path = args.output_dir / f"{args.prefix}.csv"
        json_path = args.output_dir / f"{args.prefix}.json"
        png_path = args.output_dir / f"{args.prefix}.png"
        write_csv(csv_path, rows)
        if not args.skip_png:
            plot_comparison(png_path, rows)

        report = {
            "schema": "THE100_AUAU_DUALVIEW_FINAL_COMPARISON_V1",
            "status": "PASS" if not failures else "FAIL",
            "campaign": "the100_auau_dualview_20260714",
            "definitions": {
                "ABCD": "A isolated-tight, B nonisolated-tight, C isolated-nontight, D nonisolated-nontight; 15-35 GeV integrals",
                "AD_over_BC": "A*D/(B*C); unity is the factorization reference",
                "prompt_leakage": "truth-signal leading photons in (B+C+D)/(A+B+C+D)",
                "precision": "sqrt(var(A)+var(B*C/D))/abs(A-B*C/D)",
                "subtraction_closure": "((signal+inclusive A)-BC/D)/signal A, reported as ratio and ratio-minus-one residual",
                "xJgamma_readiness": "region-A and region-C r04 recoil histograms are finite and populated in 15-35 GeV",
            },
            "canonical_png_selection": {
                "data_topdir": DATA_TOPDIR,
                "isolation": CANONICAL_ISOLATION,
                "xj_key": XJ_KEY,
                "photon_pt_GeV": [PT_MIN, PT_MAX],
                "centralities": list(CENTRALITIES),
            },
            "source_roots": {f"{role}_{view}": str(path) for (role, view), path in sources.items()},
            "strict_qa_reports": qa_reports,
            "row_count": len(rows),
            "failures": failures,
            "outputs": {
                "csv": str(csv_path),
                "json": str(json_path),
                "png": None if args.skip_png else str(png_path),
            },
            "rows": [clean_row(row) for row in rows],
        }
        json_path.write_text(json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n")
        print(json.dumps({key: report[key] for key in ("status", "row_count", "failures", "outputs")}, indent=2))
        return 0 if not failures else 1
    finally:
        for root_file in root_files.values():
            root_file.close()


if __name__ == "__main__":
    raise SystemExit(main())
