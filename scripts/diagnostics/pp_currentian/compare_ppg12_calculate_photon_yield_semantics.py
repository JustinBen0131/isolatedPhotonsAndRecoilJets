#!/usr/bin/env python3
"""Decompose PPG12 Fig. 29 purity parity into ABCD vs leakage vs solver terms.

This is a local forensic helper.  It does not read or mutate SDCC state.  The
inputs are the already-extracted PPG12/RecoilJets CSV diagnostics for the
current THE-76 pp purity campaign.
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path


BASE = Path(
    "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/ppg12PhotonYield/"
    "ppg12_photon_yield_v1_data_20260620/purity_fig29_comparison"
)
DEFAULT_OUTDIR = BASE / "globalmbd_mbddigi_componentmix_20260629_current_pp"
DEFAULT_ABCD_CSV = DEFAULT_OUTDIR / "the76_current_vs_ppg12_data_abcd_region_ratios_20260629.csv"
DEFAULT_ROOT_LEVEL_CSV = DEFAULT_OUTDIR / "the76_globalmbd_current_vs_ppg12_fig29_purity_root_level.csv"
DEFAULT_PPG12_EXTRACT_CSV = BASE / "ppg12_ratio_diagnostic/ppg12_photon_final_bdt_nom_extract.csv"


@dataclass(frozen=True)
class BinKey:
    lo: float
    hi: float

    @property
    def label(self) -> str:
        return f"{self.lo:g}-{self.hi:g}"


def ffloat(row: dict[str, str], key: str, default: float = 0.0) -> float:
    value = row.get(key, "")
    if value == "":
        return default
    return float(value)


def load_csv_by_bin(path: Path) -> dict[BinKey, dict[str, str]]:
    out: dict[BinKey, dict[str, str]] = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            out[BinKey(float(row["pt_lo"]), float(row["pt_hi"]))] = row
    return out


def load_ppg12_extract(path: Path) -> dict[tuple[str, str], dict[float, float]]:
    out: dict[tuple[str, str], dict[float, float]] = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            if row["source"] != "data":
                continue
            key = (row["kind"], row["name"])
            x = ffloat(row, "x")
            y = ffloat(row, "y")
            if row["kind"] == "HIST":
                y = ffloat(row, "content")
                x = 0.5 * (ffloat(row, "bin_lo") + ffloat(row, "bin_hi"))
            out.setdefault(key, {})[round(x, 6)] = y
    return out


def by_center(values: dict[float, float], lo: float, hi: float) -> float:
    center = round(0.5 * (lo + hi), 6)
    if center in values:
        return values[center]
    # PPG12 stores 28-32 and 32-36 at centers 30 and 34 exactly, but keep a
    # nearest-neighbor fallback for CSV floating point noise.
    nearest = min(values, key=lambda x: abs(x - center))
    if abs(nearest - center) > 1.0e-3:
        raise KeyError(f"missing value for bin center {center}")
    return values[nearest]


def myfunc_residual(x: float, a: float, b: float, c: float, d: float, cb: float, cc: float, cd: float) -> float:
    denom = d - cd * x
    if abs(denom) < 1.0e-12:
        denom = 1.0e-12 if denom >= 0 else -1.0e-12
    return x - (a - (b - cb * x) * (c - cc * x) / denom)


def solve_ppg12_root(a: float, b: float, c: float, d: float, cb: float = 0.0, cc: float = 0.0, cd: float = 0.0) -> float:
    """Solve PPG12 myfunc(x)=0 for R=1 over the same nominal [0, 2A] range."""
    if a <= 0.0 or d <= 0.0:
        return math.nan

    if abs(cb) < 1.0e-14 and abs(cc) < 1.0e-14 and abs(cd) < 1.0e-14:
        return a - b * c / d

    # (A-x)(D-cD*x) - (B-cB*x)(C-cC*x) = 0
    q2 = cd - cb * cc
    q1 = -a * cd - d + b * cc + c * cb
    q0 = a * d - b * c

    candidates: list[float] = []
    if abs(q2) < 1.0e-14:
        if abs(q1) > 1.0e-14:
            candidates.append(-q0 / q1)
    else:
        disc = q1 * q1 - 4.0 * q2 * q0
        if disc >= -1.0e-9:
            disc = max(disc, 0.0)
            sqrt_disc = math.sqrt(disc)
            candidates.append((-q1 - sqrt_disc) / (2.0 * q2))
            candidates.append((-q1 + sqrt_disc) / (2.0 * q2))

    valid = []
    for x in candidates:
        if not math.isfinite(x):
            continue
        if -1.0e-9 <= x <= 2.0 * a + 1.0e-9 and abs(d - cd * x) > 1.0e-9:
            valid.append(max(0.0, min(2.0 * a, x)))
    if not valid:
        return math.nan
    return min(valid, key=lambda x: abs(myfunc_residual(x, a, b, c, d, cb, cc, cd)))


def purity(a: float, root: float) -> float:
    if a <= 0.0 or not math.isfinite(root):
        return math.nan
    return root / a


def rms(values: list[float]) -> float:
    finite = [v for v in values if math.isfinite(v)]
    if not finite:
        return math.nan
    return math.sqrt(sum(v * v for v in finite) / len(finite))


def median(values: list[float]) -> float:
    finite = sorted(v for v in values if math.isfinite(v))
    if not finite:
        return math.nan
    n = len(finite)
    if n % 2:
        return finite[n // 2]
    return 0.5 * (finite[n // 2 - 1] + finite[n // 2])


def fmt(value: float) -> str:
    if value is None or not math.isfinite(value):
        return "nan"
    return f"{value:.6g}"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--abcd-csv", type=Path, default=DEFAULT_ABCD_CSV)
    parser.add_argument("--root-level-csv", type=Path, default=DEFAULT_ROOT_LEVEL_CSV)
    parser.add_argument("--ppg12-extract-csv", type=Path, default=DEFAULT_PPG12_EXTRACT_CSV)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--tag", default="20260629")
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)

    abcd = load_csv_by_bin(args.abcd_csv)
    root_level = load_csv_by_bin(args.root_level_csv)
    extract = load_ppg12_extract(args.ppg12_extract_csv)

    ppg12_g_raw = extract[("GRAPH", "gpurity")]
    ppg12_g_leak = extract[("GRAPH", "gpurity_leak")]
    ppg12_leak_b = extract[("HIST", "h_leak_B")]
    ppg12_leak_c = extract[("HIST", "h_leak_C")]
    ppg12_leak_d = extract[("HIST", "h_leak_D")]

    rows: list[dict[str, float | str]] = []
    for key in sorted(abcd, key=lambda k: k.lo):
        arow = abcd[key]
        rrow = root_level[key]

        cur_a = ffloat(arow, "current_A")
        cur_b = ffloat(arow, "current_B")
        cur_c = ffloat(arow, "current_C")
        cur_d = ffloat(arow, "current_D")
        ppg_a = ffloat(arow, "ppg12_A")
        ppg_b = ffloat(arow, "ppg12_B")
        ppg_c = ffloat(arow, "ppg12_C")
        ppg_d = ffloat(arow, "ppg12_D")

        cur_cb = ffloat(rrow, "current_fB")
        cur_cc = ffloat(rrow, "current_fC")
        cur_cd = ffloat(rrow, "current_fD")
        ppg_cb = by_center(ppg12_leak_b, key.lo, key.hi)
        ppg_cc = by_center(ppg12_leak_c, key.lo, key.hi)
        ppg_cd = by_center(ppg12_leak_d, key.lo, key.hi)

        ppg_graph_raw = by_center(ppg12_g_raw, key.lo, key.hi)
        ppg_graph_leak = by_center(ppg12_g_leak, key.lo, key.hi)

        cur_raw_formula = purity(cur_a, solve_ppg12_root(cur_a, cur_b, cur_c, cur_d))
        cur_leak_formula = purity(cur_a, solve_ppg12_root(cur_a, cur_b, cur_c, cur_d, cur_cb, cur_cc, cur_cd))
        cur_abcd_ppg_leak = purity(cur_a, solve_ppg12_root(cur_a, cur_b, cur_c, cur_d, ppg_cb, ppg_cc, ppg_cd))

        ppg_raw_formula = purity(ppg_a, solve_ppg12_root(ppg_a, ppg_b, ppg_c, ppg_d))
        ppg_leak_formula = purity(ppg_a, solve_ppg12_root(ppg_a, ppg_b, ppg_c, ppg_d, ppg_cb, ppg_cc, ppg_cd))
        ppg_abcd_cur_leak = purity(ppg_a, solve_ppg12_root(ppg_a, ppg_b, ppg_c, ppg_d, cur_cb, cur_cc, cur_cd))

        d_needed_raw = math.nan
        if cur_a > 0.0 and (1.0 - ppg_graph_raw) > 0.0:
            d_needed_raw = cur_b * cur_c / (cur_a * (1.0 - ppg_graph_raw))

        rows.append(
            {
                "pt_lo": key.lo,
                "pt_hi": key.hi,
                "bin": key.label,
                "current_B_over_A": cur_b / cur_a if cur_a else math.nan,
                "current_C_over_A": cur_c / cur_a if cur_a else math.nan,
                "current_D_over_A": cur_d / cur_a if cur_a else math.nan,
                "ppg12_B_over_A": ppg_b / ppg_a if ppg_a else math.nan,
                "ppg12_C_over_A": ppg_c / ppg_a if ppg_a else math.nan,
                "ppg12_D_over_A": ppg_d / ppg_a if ppg_a else math.nan,
                "current_raw_formula": cur_raw_formula,
                "current_leak_formula": cur_leak_formula,
                "current_abcd_ppg12_leak_formula": cur_abcd_ppg_leak,
                "ppg12_raw_formula": ppg_raw_formula,
                "ppg12_leak_formula": ppg_leak_formula,
                "ppg12_abcd_current_leak_formula": ppg_abcd_cur_leak,
                "ppg12_graph_raw": ppg_graph_raw,
                "ppg12_graph_leak": ppg_graph_leak,
                "ppg12_raw_formula_minus_graph": ppg_raw_formula - ppg_graph_raw,
                "ppg12_leak_formula_minus_graph": ppg_leak_formula - ppg_graph_leak,
                "current_raw_minus_ppg12_graph": cur_raw_formula - ppg_graph_raw,
                "current_leak_minus_ppg12_graph": cur_leak_formula - ppg_graph_leak,
                "current_abcd_ppg12_leak_minus_ppg12_graph": cur_abcd_ppg_leak - ppg_graph_leak,
                "ppg12_abcd_current_leak_minus_ppg12_graph": ppg_abcd_cur_leak - ppg_graph_leak,
                "current_D_over_A_div_ppg12": (cur_d / cur_a) / (ppg_d / ppg_a) if cur_a and ppg_a and ppg_d else math.nan,
                "current_D_over_A_needed_for_ppg12_raw": d_needed_raw / cur_a if cur_a else math.nan,
                "current_D_scale_needed_for_ppg12_raw": d_needed_raw / cur_d if cur_d else math.nan,
                "current_fB": cur_cb,
                "current_fC": cur_cc,
                "current_fD": cur_cd,
                "ppg12_fB": ppg_cb,
                "ppg12_fC": ppg_cc,
                "ppg12_fD": ppg_cd,
            }
        )

    csv_path = args.outdir / f"the76_calculate_photon_yield_semantics_decomp_{args.tag}.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    raw_err = [float(row["current_raw_minus_ppg12_graph"]) for row in rows]
    leak_err = [float(row["current_leak_minus_ppg12_graph"]) for row in rows]
    cur_abcd_ppg_leak_err = [float(row["current_abcd_ppg12_leak_minus_ppg12_graph"]) for row in rows]
    ppg_abcd_cur_leak_err = [float(row["ppg12_abcd_current_leak_minus_ppg12_graph"]) for row in rows]
    ppg_formula_err = [float(row["ppg12_leak_formula_minus_graph"]) for row in rows]
    d_ratio = [float(row["current_D_over_A_div_ppg12"]) for row in rows]

    md_path = args.outdir / f"the76_calculate_photon_yield_semantics_decomp_{args.tag}.md"
    lines = [
        "# THE-76 PPG12 CalculatePhotonYield semantics decomposition",
        "",
        f"ABCD CSV: `{args.abcd_csv}`",
        f"Root-level purity CSV: `{args.root_level_csv}`",
        f"PPG12 extract CSV: `{args.ppg12_extract_csv}`",
        "",
        "## Executive result",
        "",
        (
            "The PPG12 `CalculatePhotonYield.C` equation and graph convention are not the dominant "
            "remaining discrepancy.  Re-evaluating PPG12's own ABCD and leakage inputs with the "
            "same R=1 root equation reproduces the PPG12 Fig.29 leakage-corrected graph to within "
            f"RMS {fmt(rms(ppg_formula_err))} purity."
        ),
        "",
        (
            "Replacing only the leakage fractions does not repair the current points: current ABCD "
            f"+ PPG12 leakage still differs from PPG12 by RMS {fmt(rms(cur_abcd_ppg_leak_err))}.  "
            "Replacing the data ABCD side with PPG12 while keeping current leakage is much closer "
            f"(RMS {fmt(rms(ppg_abcd_cur_leak_err))})."
        ),
        "",
        (
            "The strongest data-side symptom is D/A.  Median current D/A divided by PPG12 D/A is "
            f"{fmt(median(d_ratio))}; in the high-ET bins it is far below one."
        ),
        "",
        "## Error scales",
        "",
        "| comparison | RMS purity error vs PPG12 graph | median signed error |",
        "| --- | ---: | ---: |",
        f"| current raw formula vs PPG12 raw graph | {fmt(rms(raw_err))} | {fmt(median(raw_err))} |",
        f"| current leakage formula vs PPG12 leakage graph | {fmt(rms(leak_err))} | {fmt(median(leak_err))} |",
        f"| current ABCD + PPG12 leakage vs PPG12 leakage graph | {fmt(rms(cur_abcd_ppg_leak_err))} | {fmt(median(cur_abcd_ppg_leak_err))} |",
        f"| PPG12 ABCD + current leakage vs PPG12 leakage graph | {fmt(rms(ppg_abcd_cur_leak_err))} | {fmt(median(ppg_abcd_cur_leak_err))} |",
        f"| PPG12 ABCD + PPG12 leakage formula vs PPG12 leakage graph | {fmt(rms(ppg_formula_err))} | {fmt(median(ppg_formula_err))} |",
        "",
        "## Bin details",
        "",
        (
            "| bin | current D/A | PPG12 D/A | D/A ratio | current raw | PPG12 raw graph | "
            "current corr | PPG12 corr graph | current ABCD + PPG12 leak | PPG12 ABCD + current leak |"
        ),
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in rows:
        lines.append(
            f"| {row['bin']} | {fmt(float(row['current_D_over_A']))} | {fmt(float(row['ppg12_D_over_A']))} | "
            f"{fmt(float(row['current_D_over_A_div_ppg12']))} | {fmt(float(row['current_raw_formula']))} | "
            f"{fmt(float(row['ppg12_graph_raw']))} | {fmt(float(row['current_leak_formula']))} | "
            f"{fmt(float(row['ppg12_graph_leak']))} | {fmt(float(row['current_abcd_ppg12_leak_formula']))} | "
            f"{fmt(float(row['ppg12_abcd_current_leak_formula']))} |"
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "- The signal-only leakage campaign is no longer the first blocker for Fig.29 parity.",
            "- The remaining parity break is upstream in the pp data ABCD population, especially the non-tight non-isolated D region.",
            "- The next local test should compare PPG12 `MakeDataHisto`/CaloAna24 data-side classification against RecoilJets pp data region fills for the same event or, at minimum, the same per-bin BDT/isolation gates.",
            "- A pp data rerun is only justified after that data-side region-fill difference is identified and patched in a tiny foreground/local test.",
            "",
        ]
    )
    md_path.write_text("\n".join(lines))

    print(f"wrote {csv_path}")
    print(f"wrote {md_path}")
    print(f"ppg12_formula_rms={fmt(rms(ppg_formula_err))}")
    print(f"current_abcd_ppg12_leak_rms={fmt(rms(cur_abcd_ppg_leak_err))}")
    print(f"ppg12_abcd_current_leak_rms={fmt(rms(ppg_abcd_cur_leak_err))}")
    print(f"median_current_DA_over_ppg12_DA={fmt(median(d_ratio))}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
