#!/usr/bin/env python3
"""Diagnose whether PPG12 leakage C/D residuals are population or isolation split.

Consumes the Fig.29 current-vs-PPG12 leakage ratio table and writes a compact
sideband-population diagnostic next to the existing B/C/D overlay.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
DEFAULT_INPUT = (
    REPO
    / "dataOutput/ppg12PhotonYield/ppg12_photon_yield_v1_data_20260620"
    / "purity_fig29_comparison/ppg12_ratio_diagnostic"
    / "leakage_components/current_vs_ppg12_fig29_leakage_components_bcd.csv"
)
DEFAULT_OUTDIR = DEFAULT_INPUT.parent


def read_rows(path: Path) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            rows.append({key: float(value) for key, value in row.items()})
    return rows


def safe_div(num: float, den: float) -> float:
    return num / den if den else float("nan")


def augment_rows(rows: list[dict[str, float]]) -> list[dict[str, float]]:
    out: list[dict[str, float]] = []
    for row in rows:
        current_nt = row["current_f_c"] + row["current_f_d"]
        ppg12_nt = row["ppg12_f_c"] + row["ppg12_f_d"]
        current_d_share = safe_div(row["current_f_d"], current_nt)
        ppg12_d_share = safe_div(row["ppg12_f_d"], ppg12_nt)
        current_c_share = safe_div(row["current_f_c"], current_nt)
        ppg12_c_share = safe_div(row["ppg12_f_c"], ppg12_nt)
        out.append(
            {
                **row,
                "current_f_nt": current_nt,
                "ppg12_f_nt": ppg12_nt,
                "current_f_nt_over_ppg12": safe_div(current_nt, ppg12_nt),
                "current_d_share_of_nt": current_d_share,
                "ppg12_d_share_of_nt": ppg12_d_share,
                "current_d_share_over_ppg12": safe_div(current_d_share, ppg12_d_share),
                "current_c_share_of_nt": current_c_share,
                "ppg12_c_share_of_nt": ppg12_c_share,
                "current_c_share_over_ppg12": safe_div(current_c_share, ppg12_c_share),
            }
        )
    return out


def finite_range(rows: list[dict[str, float]], key: str) -> tuple[float, float]:
    vals = [row[key] for row in rows if row[key] == row[key]]
    if not vals:
        return float("nan"), float("nan")
    return min(vals), max(vals)


def write_csv(rows: list[dict[str, float]], path: Path) -> None:
    fields = [
        "pt_lo",
        "pt_hi",
        "pt_center",
        "current_f_b_over_ppg12",
        "current_f_c_over_ppg12",
        "current_f_d_over_ppg12",
        "current_f_nt",
        "ppg12_f_nt",
        "current_f_nt_over_ppg12",
        "current_c_share_of_nt",
        "ppg12_c_share_of_nt",
        "current_c_share_over_ppg12",
        "current_d_share_of_nt",
        "ppg12_d_share_of_nt",
        "current_d_share_over_ppg12",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row[field] for field in fields})


def write_summary(rows: list[dict[str, float]], path: Path, source: Path) -> dict[str, object]:
    stable = [row for row in rows if row["pt_hi"] <= 26.0]
    summary = {
        "source_csv": str(source),
        "stable_bin_definition": "pt_hi <= 26 GeV",
        "stable_ranges": {
            "B_current_over_PPG12": finite_range(stable, "current_f_b_over_ppg12"),
            "C_current_over_PPG12": finite_range(stable, "current_f_c_over_ppg12"),
            "D_current_over_PPG12": finite_range(stable, "current_f_d_over_ppg12"),
            "non_tight_total_current_over_PPG12": finite_range(stable, "current_f_nt_over_ppg12"),
            "D_share_within_non_tight_current_over_PPG12": finite_range(stable, "current_d_share_over_ppg12"),
        },
        "interpretation": [
            "B/A is now aligned with PPG12 in the stable bins.",
            "C and D are both low because the total non-tight signal population (C+D)/A is low.",
            "D/(C+D) is much closer to PPG12 than D/A, so D is not primarily a separate isolation split failure.",
            "The next fix should target non-tight candidate/BDT-sideband parity, not another broad topo-isolation rerun.",
        ],
    }
    path.with_suffix(".json").write_text(json.dumps(summary, indent=2) + "\n")
    with path.open("w") as handle:
        handle.write("# PPG12 Fig.29 non-tight leakage population diagnostic\n\n")
        handle.write(f"- Source: `{source}`\n")
        handle.write("- Purpose: separate missing non-tight signal population from wrong C/D isolation split.\n\n")
        handle.write("## Stable-bin ranges\n\n")
        for key, value in summary["stable_ranges"].items():
            lo, hi = value
            handle.write(f"- `{key}`: `{lo:.3f}` to `{hi:.3f}`\n")
        handle.write("\n## Interpretation\n\n")
        for item in summary["interpretation"]:
            handle.write(f"- {item}\n")
        handle.write("\n## Table\n\n")
        handle.write("| ET bin | B/PPG12 | C/PPG12 | D/PPG12 | (C+D)/PPG12 | current D/(C+D) | PPG12 D/(C+D) | D-share ratio |\n")
        handle.write("| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |\n")
        for row in rows:
            handle.write(
                f"| {row['pt_lo']:.0f}-{row['pt_hi']:.0f} | "
                f"{row['current_f_b_over_ppg12']:.3f} | "
                f"{row['current_f_c_over_ppg12']:.3f} | "
                f"{row['current_f_d_over_ppg12']:.3f} | "
                f"{row['current_f_nt_over_ppg12']:.3f} | "
                f"{row['current_d_share_of_nt']:.3f} | "
                f"{row['ppg12_d_share_of_nt']:.3f} | "
                f"{row['current_d_share_over_ppg12']:.3f} |\n"
            )
    return summary


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    rows = augment_rows(read_rows(args.input))
    if not rows:
        raise RuntimeError(f"no rows read from {args.input}")

    stem = "current_vs_ppg12_fig29_nontight_population_diagnostic"
    csv_path = args.outdir / f"{stem}.csv"
    md_path = args.outdir / f"{stem}.md"
    write_csv(rows, csv_path)
    write_summary(rows, md_path, args.input)
    print(csv_path)
    print(md_path)
    print(md_path.with_suffix(".json"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
