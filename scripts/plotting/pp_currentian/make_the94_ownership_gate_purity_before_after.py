#!/usr/bin/env python3
"""Compare the pre-gate and ownership-gated inclusive-SIM purity results.

Both inputs use the same PPG12 SDCC reference and the same unsuffixed
ABCD estimator contract.  Their difference therefore isolates the change in
the registered inclusive-simulation production artifact.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


REPO = Path(__file__).resolve().parents[3]
DEFAULT_PRIOR = (
    REPO
    / "dataOutput/ppg12Parity/the97_ppg12_si_contract_restore_full_20260715_1420"
    / "ian_current_sim_refresh_20260716/inclusive_purity_unsuffixed_fix"
    / "ppg12_ian_fig3_purity_sim_sdcc_vs_current_unsuffixed_abcd_overlay_ratio_points.csv"
)
DEFAULT_CURRENT = (
    REPO
    / "dataOutput/ppg12Parity/the94_ppg12_inclusive_ownership_gate_final_20260717_1349EDT"
    / "purity_comparison"
    / "ppg12_ian_fig3_purity_sim_sdcc_vs_current_unsuffixed_abcd_overlay_ratio_points.csv"
)
DEFAULT_OUTDIR = DEFAULT_CURRENT.parent

SERIES = (
    ("truth", "Truth purity", "#d62728"),
    ("raw", "Raw ABCD purity", "#222222"),
    ("corrected", "Leakage-corrected purity", "#2457e6"),
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_points(path: Path) -> dict[str, list[dict[str, float]]]:
    rows: dict[str, list[dict[str, float]]] = {key: [] for key, _, _ in SERIES}
    with path.open(newline="") as handle:
        for raw in csv.DictReader(handle):
            key = raw["series"]
            if key not in rows:
                continue
            rows[key].append({name: float(value) for name, value in raw.items() if name != "series"})
    for key in rows:
        rows[key].sort(key=lambda row: row["x_gev"])
    return rows


def arrays(rows: list[dict[str, float]], prefix: str) -> tuple[np.ndarray, ...]:
    x = np.array([row["x_gev"] for row in rows])
    y = np.array([row[prefix] for row in rows])
    if prefix == "ppg12_sdcc":
        low = np.array([row["ppg12_error_low"] for row in rows])
        high = np.array([row["ppg12_error_high"] for row in rows])
    else:
        low = high = np.array([row["current_error"] for row in rows])
    return x, y, low, high


def require_common_contract(
    prior: dict[str, list[dict[str, float]]],
    current: dict[str, list[dict[str, float]]],
) -> None:
    for key, _, _ in SERIES:
        if len(prior[key]) != len(current[key]):
            raise RuntimeError(f"{key}: point-count mismatch")
        for old, new in zip(prior[key], current[key]):
            for field in ("x_gev", "ppg12_sdcc", "ppg12_error_low", "ppg12_error_high"):
                if abs(old[field] - new[field]) > 1.0e-12:
                    raise RuntimeError(f"{key}: PPG12 reference mismatch in {field}")


def render(
    prior: dict[str, list[dict[str, float]]],
    current: dict[str, list[dict[str, float]]],
    output: Path,
    prior_label: str,
    current_label: str,
    title: str,
) -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.linewidth": 1.0,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )
    fig = plt.figure(figsize=(13.2, 5.1), constrained_layout=False)
    grid = fig.add_gridspec(2, 3, height_ratios=(3.2, 1.25), hspace=0.05, wspace=0.25)

    for column, (key, title, color) in enumerate(SERIES):
        top = fig.add_subplot(grid[0, column])
        bottom = fig.add_subplot(grid[1, column], sharex=top)

        x, ref, ref_lo, ref_hi = arrays(current[key], "ppg12_sdcc")
        _, old, old_lo, old_hi = arrays(prior[key], "current")
        _, new, new_lo, new_hi = arrays(current[key], "current")

        top.errorbar(
            x,
            ref,
            yerr=np.vstack((ref_lo, ref_hi)),
            xerr=0.78,
            fmt="o",
            ms=4.4,
            mfc="white",
            mec=color,
            ecolor=color,
            elinewidth=0.9,
            capsize=2,
            label="PPG12 SDCC",
        )
        top.errorbar(
            x - 0.12,
            old,
            yerr=np.vstack((old_lo, old_hi)),
            fmt="s",
            ms=3.8,
            color="#777777",
            ecolor="#777777",
            elinewidth=0.8,
            capsize=1.8,
            label=prior_label,
        )
        top.errorbar(
            x + 0.12,
            new,
            yerr=np.vstack((new_lo, new_hi)),
            fmt="o",
            ms=4.3,
            color=color,
            ecolor=color,
            elinewidth=0.9,
            capsize=2,
            label=current_label,
        )
        top.set_title(title, fontsize=11, fontweight="bold", pad=7)
        top.set_xlim(9.8, 36.2)
        all_low = np.concatenate((ref - ref_lo, old - old_lo, new - new_lo))
        all_high = np.concatenate((ref + ref_hi, old + old_hi, new + new_hi))
        pad = max(0.04, 0.10 * (all_high.max() - all_low.min()))
        top.set_ylim(max(0.0, all_low.min() - pad), min(1.2, all_high.max() + pad))
        top.grid(axis="y", color="#dddddd", linewidth=0.55, alpha=0.8)
        top.tick_params(labelbottom=False)
        if column == 0:
            top.set_ylabel("Purity")
            top.text(0.04, 0.95, r"$\bf{sPHENIX}$ Internal", transform=top.transAxes, va="top")
            top.text(0.04, 0.86, r"$p+p$, $\sqrt{s}=200$ GeV", transform=top.transAxes, va="top")
        if column == 2:
            top.legend(loc="lower right", fontsize=8.2, frameon=False, handletextpad=0.5)

        safe = ref != 0.0
        old_ratio = np.divide(old, ref, out=np.full_like(old, np.nan), where=safe)
        new_ratio = np.divide(new, ref, out=np.full_like(new, np.nan), where=safe)
        old_ratio_err = old_ratio * np.sqrt((old_lo / old) ** 2 + (ref_lo / ref) ** 2)
        new_ratio_err = new_ratio * np.sqrt((new_lo / new) ** 2 + (ref_lo / ref) ** 2)
        bottom.axhline(1.0, color="#555555", linestyle="--", linewidth=0.85)
        bottom.errorbar(x - 0.12, old_ratio, yerr=old_ratio_err, fmt="s", ms=3.5, color="#777777", capsize=1.5)
        bottom.errorbar(x + 0.12, new_ratio, yerr=new_ratio_err, fmt="o", ms=3.9, color=color, capsize=1.5)
        envelope = np.concatenate((old_ratio - old_ratio_err, old_ratio + old_ratio_err, new_ratio - new_ratio_err, new_ratio + new_ratio_err))
        lower = min(0.95, float(np.nanmin(envelope)))
        upper = max(1.05, float(np.nanmax(envelope)))
        span = max(upper - lower, 0.15)
        bottom.set_ylim(max(0.0, lower - 0.10 * span), upper + 0.10 * span)
        bottom.set_xlabel(r"$E_{T}^{\gamma,\mathrm{rec}}$ [GeV]")
        bottom.grid(axis="y", color="#e4e4e4", linewidth=0.5, alpha=0.8)
        if column == 0:
            bottom.set_ylabel("Current / PPG12", fontsize=9)
        else:
            bottom.tick_params(labelleft=False)

    fig.suptitle(
        title,
        fontsize=13,
        fontweight="bold",
        y=0.985,
    )
    fig.text(
        0.5,
        0.012,
        "Same PPG12 reference and unsuffixed A/B/C/D estimator in both RecoilJets constructions",
        ha="center",
        fontsize=9,
    )
    fig.subplots_adjust(left=0.065, right=0.99, top=0.88, bottom=0.12)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180, facecolor="white")
    plt.close(fig)


def write_comparison_csv(
    output: Path,
    prior: dict[str, list[dict[str, float]]],
    current: dict[str, list[dict[str, float]]],
) -> None:
    fields = (
        "series",
        "x_gev",
        "ppg12_sdcc",
        "prior_current",
        "ownership_gated_current",
        "prior_over_ppg12",
        "ownership_gated_over_ppg12",
        "ownership_gated_minus_prior",
    )
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for key, _, _ in SERIES:
            for old, new in zip(prior[key], current[key]):
                reference = new["ppg12_sdcc"]
                writer.writerow(
                    {
                        "series": key,
                        "x_gev": new["x_gev"],
                        "ppg12_sdcc": reference,
                        "prior_current": old["current"],
                        "ownership_gated_current": new["current"],
                        "prior_over_ppg12": old["current"] / reference,
                        "ownership_gated_over_ppg12": new["current"] / reference,
                        "ownership_gated_minus_prior": new["current"] - old["current"],
                    }
                )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prior-points", type=Path, default=DEFAULT_PRIOR)
    parser.add_argument("--current-points", type=Path, default=DEFAULT_CURRENT)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--output-stem", default="ownership_gate_purity_before_after")
    parser.add_argument("--prior-label", default="Prior registered output")
    parser.add_argument("--current-label", default="Ownership-gated output")
    parser.add_argument(
        "--title",
        default="Inclusive-simulation purity: effect of the truth-jet ownership contract",
    )
    parser.add_argument(
        "--interpretation",
        default=(
            "The comparison isolates the change from the registered pre-gate inclusive output "
            "to the ownership-gated output. It does not attribute the earlier above-unity "
            "corrected-purity artifact to the ownership gate; that artifact was removed by "
            "using the PPG12-equivalent unsuffixed estimator populations."
        ),
    )
    args = parser.parse_args()

    prior = read_points(args.prior_points)
    current = read_points(args.current_points)
    require_common_contract(prior, current)

    png = args.outdir / f"{args.output_stem}.png"
    table = args.outdir / f"{args.output_stem}.csv"
    manifest = args.outdir / f"{args.output_stem}_manifest.json"
    render(prior, current, png, args.prior_label, args.current_label, args.title)
    write_comparison_csv(table, prior, current)
    payload = {
        "status": "production_contract_comparison",
        "observable": "inclusive-simulation truth, raw-ABCD, and leakage-corrected purity",
        "ppg12_reference": "identical stored SDCC arrays in both inputs",
        "estimator_contract": "identical unsuffixed A/B/C/D PPG12 estimator in both inputs",
        "prior_points": str(args.prior_points),
        "prior_points_sha256": sha256(args.prior_points),
        "ownership_gated_points": str(args.current_points),
        "ownership_gated_points_sha256": sha256(args.current_points),
        "prior_label": args.prior_label,
        "current_label": args.current_label,
        "comparison_csv": str(table),
        "comparison_csv_sha256": sha256(table),
        "png": str(png),
        "png_sha256": sha256(png),
        "interpretation": args.interpretation,
    }
    manifest.write_text(json.dumps(payload, indent=2) + "\n")
    print(png)
    print(table)
    print(manifest)


if __name__ == "__main__":
    main()
