#!/usr/bin/env python3
"""Render corrected-SI photon-SIM signal leakage and ABCD parity diagnostics."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


REGIONS = (
    ("A", "tight isolated", "black", "o"),
    ("B", "tight nonisolated", "#d62728", "s"),
    ("C", "nontight isolated", "#1f77b4", "^"),
    ("D", "nontight nonisolated", "#9467bd", "D"),
)
LEAKAGE = (
    ("B", r"$c_B=B_{sig}/A_{sig}$", "black", "o"),
    ("C", r"$c_C=C_{sig}/A_{sig}$", "#d62728", "s"),
    ("D", r"$c_D=D_{sig}/A_{sig}$", "#1f77b4", "^"),
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def value_error_ratio(n: float, ne: float, d: float, de: float) -> tuple[float, float]:
    if d == 0:
        return math.nan, math.nan
    ratio = n / d
    variance = (ne / d) ** 2 + ((n * de) / (d * d)) ** 2
    return ratio, math.sqrt(max(0.0, variance))


def load_rows(path: Path) -> list[dict[str, float]]:
    with path.open(newline="") as handle:
        raw = list(csv.DictReader(handle))
    if len(raw) != 11:
        raise RuntimeError(f"expected 11 bins, found {len(raw)}")
    rows: list[dict[str, float]] = []
    required = {"bin_index", "pt_lo", "pt_hi"}
    for source in ("ppg12", "current"):
        for region, *_ in REGIONS:
            required.update({f"{source}_{region}", f"{source}_{region}_err"})
    missing = required - set(raw[0])
    if missing:
        raise RuntimeError(f"missing fields: {sorted(missing)}")
    for item in raw:
        row = {key: float(value) for key, value in item.items()}
        row["bin_index"] = int(float(item["bin_index"]))
        row["pt_center"] = 0.5 * (row["pt_lo"] + row["pt_hi"])
        row["pt_half_width"] = 0.5 * (row["pt_hi"] - row["pt_lo"])
        rows.append(row)
    expected = [(10, 12), (12, 14), (14, 16), (16, 18), (18, 20), (20, 22), (22, 24), (24, 26), (26, 28), (28, 32), (32, 36)]
    got = [(row["pt_lo"], row["pt_hi"]) for row in rows]
    if got != expected:
        raise RuntimeError(f"unexpected binning: {got}")
    return rows


def configure() -> None:
    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
        "font.size": 12,
        "axes.linewidth": 1.0,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "legend.frameon": False,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
    })


def annotate(ax: plt.Axes, subtitle: str, *, compact_top: bool = False) -> None:
    if compact_top:
        y_internal, y_system, y_eta, y_subtitle = 0.220, 0.170, 0.120, 0.050
    else:
        y_internal, y_system, y_eta, y_subtitle = 0.965, 0.895, 0.835, 0.775
    ax.text(0.025, y_internal, r"$\bf{\it{sPHENIX}}$ Internal", transform=ax.transAxes, va="top", ha="left", fontsize=13)
    ax.text(0.025, y_system, r"$p+p\ \sqrt{s}=200\ \mathrm{GeV}$", transform=ax.transAxes, va="top", ha="left")
    ax.text(0.025, y_eta, r"$|\eta^\gamma|<0.7$", transform=ax.transAxes, va="top", ha="left")
    ax.text(0.025, y_subtitle, subtitle, transform=ax.transAxes, va="top", ha="left", fontsize=10.5)


def write_csv(path: Path, rows: list[dict[str, float | int | str]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def leakage_products(rows: list[dict[str, float]], outdir: Path) -> tuple[Path, Path, dict]:
    derived: list[dict[str, float | int | str]] = []
    for row in rows:
        for region, label, *_ in LEAKAGE:
            ppg, ppge = value_error_ratio(row[f"ppg12_{region}"], row[f"ppg12_{region}_err"], row["ppg12_A"], row["ppg12_A_err"])
            cur, cure = value_error_ratio(row[f"current_{region}"], row[f"current_{region}_err"], row["current_A"], row["current_A_err"])
            ratio, ratioe = value_error_ratio(cur, cure, ppg, ppge)
            derived.append({
                "bin_index": row["bin_index"], "pt_lo": row["pt_lo"], "pt_hi": row["pt_hi"],
                "region": region, "definition": label,
                "ppg12_leakage": ppg, "ppg12_error": ppge,
                "current_leakage": cur, "current_error": cure,
                "current_over_ppg12": ratio, "current_over_ppg12_error": ratioe,
            })
    csv_path = outdir / "corrected_si_candidate_signal_leakage_points.csv"
    write_csv(csv_path, derived)

    fig, (ax, ratio_ax) = plt.subplots(2, 1, figsize=(9.2, 8.0), gridspec_kw={"height_ratios": [3.1, 1.0], "hspace": 0.03}, sharex=True)
    for region, label, color, marker in LEAKAGE:
        vals = [item for item in derived if item["region"] == region]
        x = np.array([(item["pt_lo"] + item["pt_hi"]) / 2 for item in vals])
        xe = np.array([(item["pt_hi"] - item["pt_lo"]) / 2 for item in vals])
        p = np.array([item["ppg12_leakage"] for item in vals])
        pe = np.array([item["ppg12_error"] for item in vals])
        c = np.array([item["current_leakage"] for item in vals])
        ce = np.array([item["current_error"] for item in vals])
        r = np.array([item["current_over_ppg12"] for item in vals])
        re = np.array([item["current_over_ppg12_error"] for item in vals])
        ax.errorbar(x - 0.12, p, xerr=xe, yerr=pe, fmt=marker, ms=5.5, mfc="white", mec=color, ecolor=color, color=color, capsize=0, linestyle="none", label=f"PPG12 {label}")
        ax.errorbar(x + 0.12, c, xerr=xe, yerr=ce, fmt=marker, ms=5.5, mfc=color, mec=color, ecolor=color, color=color, capsize=0, linestyle="none", label=f"Current {label}")
        ratio_ax.errorbar(x, r, xerr=xe, yerr=re, fmt=marker, ms=5.2, mfc=color, mec=color, ecolor=color, color=color, capsize=0, linestyle="none", label=label)
    finite_top = [item[key] for item in derived for key in ("ppg12_leakage", "current_leakage") if math.isfinite(item[key])]
    ax.set_ylim(0, max(0.15, 1.28 * max(finite_top)))
    ax.set_ylabel("Signal leakage")
    annotate(ax, "Corrected-SI full-stat candidate (SI+DI)")
    ax.legend(loc="upper center", bbox_to_anchor=(0.67, 0.99), ncol=2, fontsize=9.2, columnspacing=1.0, handletextpad=0.45)
    ax.tick_params(labelbottom=False)
    finite_ratio = [item["current_over_ppg12"] for item in derived if math.isfinite(item["current_over_ppg12"])]
    lo, hi = min(finite_ratio), max(finite_ratio)
    span = max(0.15, hi - lo)
    ratio_ax.set_ylim(max(0, lo - 0.25 * span), hi + 0.25 * span)
    ratio_ax.axhline(1.0, color="0.45", linewidth=1.0, linestyle="--")
    ratio_ax.set_ylabel("Current / PPG12")
    ratio_ax.set_xlabel(r"$E_T^{\gamma,\mathrm{rec}}$ [GeV]")
    ratio_ax.set_xlim(10, 36)
    fig.subplots_adjust(left=0.12, right=0.98, top=0.98, bottom=0.11)
    png = outdir / "corrected_si_candidate_signal_leakage_ppg12_overlay_ratio.png"
    fig.savefig(png, dpi=180)
    plt.close(fig)
    return png, csv_path, {"ratio_min": lo, "ratio_max": hi}


def abcd_products(rows: list[dict[str, float]], outdir: Path) -> tuple[Path, Path, dict]:
    derived: list[dict[str, float | int | str]] = []
    for row in rows:
        for region, label, *_ in REGIONS:
            ratio, ratioe = value_error_ratio(row[f"current_{region}"], row[f"current_{region}_err"], row[f"ppg12_{region}"], row[f"ppg12_{region}_err"])
            derived.append({
                "bin_index": row["bin_index"], "pt_lo": row["pt_lo"], "pt_hi": row["pt_hi"],
                "region": region, "definition": label,
                "ppg12_count": row[f"ppg12_{region}"], "ppg12_error": row[f"ppg12_{region}_err"],
                "current_count": row[f"current_{region}"], "current_error": row[f"current_{region}_err"],
                "current_over_ppg12": ratio, "current_over_ppg12_error": ratioe,
            })
    csv_path = outdir / "corrected_si_candidate_signal_abcd_counts_points.csv"
    write_csv(csv_path, derived)

    fig, (ax, ratio_ax) = plt.subplots(2, 1, figsize=(9.2, 8.0), gridspec_kw={"height_ratios": [3.1, 1.0], "hspace": 0.03}, sharex=True)
    for region, label, color, marker in REGIONS:
        vals = [item for item in derived if item["region"] == region]
        x = np.array([(item["pt_lo"] + item["pt_hi"]) / 2 for item in vals])
        xe = np.array([(item["pt_hi"] - item["pt_lo"]) / 2 for item in vals])
        p = np.array([item["ppg12_count"] for item in vals])
        pe = np.array([item["ppg12_error"] for item in vals])
        c = np.array([item["current_count"] for item in vals])
        ce = np.array([item["current_error"] for item in vals])
        r = np.array([item["current_over_ppg12"] for item in vals])
        re = np.array([item["current_over_ppg12_error"] for item in vals])
        ax.errorbar(x - 0.12, p, xerr=xe, yerr=pe, fmt=marker, ms=5.2, mfc="white", mec=color, ecolor=color, color=color, capsize=0, linestyle="none", label=f"PPG12 {region}: {label}")
        ax.errorbar(x + 0.12, c, xerr=xe, yerr=ce, fmt=marker, ms=5.2, mfc=color, mec=color, ecolor=color, color=color, capsize=0, linestyle="none", label=f"Current {region}: {label}")
        ratio_ax.errorbar(x, r, xerr=xe, yerr=re, fmt=marker, ms=5.2, mfc=color, mec=color, ecolor=color, color=color, capsize=0, linestyle="none", label=region)
    positive = [item[key] for item in derived for key in ("ppg12_count", "current_count") if item[key] > 0]
    ax.set_yscale("log")
    ax.set_ylim(max(1e-3, min(positive) / 3), max(positive) * 6)
    ax.set_ylabel("Event-preweighted signal count")
    annotate(ax, "Corrected-SI full-stat candidate (SI+DI), no rescaling", compact_top=True)
    ax.legend(loc="upper center", bbox_to_anchor=(0.66, 0.99), ncol=2, fontsize=8.8, columnspacing=0.9, handletextpad=0.4)
    ax.tick_params(labelbottom=False)
    finite_ratio = [item["current_over_ppg12"] for item in derived if math.isfinite(item["current_over_ppg12"])]
    lo, hi = min(finite_ratio), max(finite_ratio)
    span = max(0.15, hi - lo)
    ratio_ax.set_ylim(max(0, lo - 0.25 * span), hi + 0.25 * span)
    ratio_ax.axhline(1.0, color="0.45", linewidth=1.0, linestyle="--")
    ratio_ax.set_ylabel("Current / PPG12")
    ratio_ax.set_xlabel(r"$E_T^{\gamma,\mathrm{rec}}$ [GeV]")
    ratio_ax.set_xlim(10, 36)
    fig.subplots_adjust(left=0.13, right=0.98, top=0.98, bottom=0.11)
    png = outdir / "corrected_si_candidate_signal_abcd_counts_ppg12_overlay_ratio.png"
    fig.savefig(png, dpi=180)
    plt.close(fig)
    return png, csv_path, {"ratio_min": lo, "ratio_max": hi}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-csv", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--candidate-root", required=True)
    parser.add_argument("--candidate-root-sha256", required=True)
    parser.add_argument("--ppg12-root", required=True)
    parser.add_argument("--ppg12-root-sha256", required=True)
    parser.add_argument("--merge-audit-json", required=True)
    args = parser.parse_args()
    configure()
    rows = load_rows(args.input_csv)
    args.outdir.mkdir(parents=True, exist_ok=True)
    leak_png, leak_csv, leak_metrics = leakage_products(rows, args.outdir)
    abcd_png, abcd_csv, abcd_metrics = abcd_products(rows, args.outdir)
    manifest = {
        "schema": "THE97_CORRECTED_SI_FULLSTAT_PHOTON_SIGNAL_PARITY_V1",
        "status": "corrected-SI full-stat candidate diagnostic; neutral and noncanonical",
        "candidate_root": args.candidate_root,
        "candidate_root_sha256": args.candidate_root_sha256,
        "ppg12_root": args.ppg12_root,
        "ppg12_root_sha256": args.ppg12_root_sha256,
        "source_points_csv": str(args.input_csv),
        "source_points_csv_sha256": sha256(args.input_csv),
        "object_paths": {
            "ppg12": {
                "A": "h_tight_iso_cluster_signal_0",
                "B": "h_tight_noniso_cluster_signal_0",
                "C": "h_nontight_iso_cluster_signal_0",
                "D": "h_nontight_noniso_cluster_signal_0",
            },
            "current": {
                "A": "SIM/h_tight_iso_cluster_signal_0",
                "B": "SIM/h_tight_noniso_cluster_signal_0",
                "C": "SIM/h_nontight_iso_cluster_signal_0",
                "D": "SIM/h_nontight_noniso_cluster_signal_0",
            },
        },
        "bin_edges_gev": [10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 32, 36],
        "comparison": "absolute event-preweighted signal counts and cB/cC/cD; no normalization, rescaling, fitting, reweighting, or retuning",
        "ratio_definition": "Current / PPG12 with propagated independent statistical errors",
        "merge_audit_json": args.merge_audit_json,
        "leakage": {"png": str(leak_png), "csv": str(leak_csv), **leak_metrics},
        "signal_abcd": {"png": str(abcd_png), "csv": str(abcd_csv), **abcd_metrics},
    }
    manifest_path = args.outdir / "corrected_si_candidate_signal_parity_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(leak_png)
    print(abcd_png)
    print(manifest_path)


if __name__ == "__main__":
    main()
