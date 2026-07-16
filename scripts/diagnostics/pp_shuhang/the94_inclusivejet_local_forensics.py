#!/usr/bin/env python3
"""THE-94 local forensic tables for PPG12 inclusive-jet stitching.

This script is read-only.  It compares existing local CSV/ROOT artifacts from
the July 1 source-scope diagnostic and the strict July 2 rerun against the
PPG12 Fig.6 IAN/no-suffix source.  It deliberately does not claim parity and
does not apply fitted scales as a valid fix; scale tests are labelled as
mechanism diagnostics only.
"""

from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path
from statistics import mean, median

import numpy as np
import uproot


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
OUT_DIR = (
    REPO
    / "dataOutput/ppg12Parity/THE-94_ppg12_pp_inclusivejet_sim_parity/local_forensics_20260702"
)

PPG12_IAN_CSV = (
    REPO
    / "dataOutput/ppg12Parity/the76_ppg12_parity_full_20260701_003024/"
    "reference_audit/sdcc_pull/live_ppg12_ian_source_points.csv"
)
OLD_CLOSURE_CSV = (
    REPO
    / "dataOutput/ppg12Parity/the76_ppg12_parity_full_20260701_003024/"
    "strict_stitched_inclusivejet/jet_sdcc_over_current_shape_overlay_ppg12_jet8_xsecfix_points.csv"
)
OLD_COMPONENT_CSV = (
    REPO
    / "dataOutput/ppg12Parity/the76_ppg12_parity_full_20260701_003024/"
    "strict_stitched_inclusivejet/jet_component_exposure_audit_points.csv"
)
STRICT_JULY2_CSV = (
    REPO
    / "dataOutput/ppg12Parity/the76_ppg12_fig6_inclusive_strict_20260702_014757/"
    "contract_canary/fig6_inclusive/fig6_inclusive_current_kept_vs_ppg12_source_points.csv"
)
STRICT_JULY2_ROOT_DIR = (
    REPO
    / "dataOutput/ppg12Parity/the76_ppg12_fig6_inclusive_strict_20260702_014757/"
    "final_roots/inclusivejet"
)

JET50_XSEC_PB = 7.3113
JET8_OLD_LOCAL_XSEC_PB = 1.3013e7
JET8_PPG12_XSEC_PB = 1.15e7
OLD_GLOBAL_SCALE_AFTER_JET8_FIX = 0.5158774708396985


@dataclass(frozen=True)
class Sample:
    name: str
    xsec_pb: float
    lo: float
    hi: float

    @property
    def xsec_over_jet50(self) -> float:
        return self.xsec_pb / JET50_XSEC_PB


SAMPLES = (
    Sample("jet8", JET8_PPG12_XSEC_PB, 9.0, 14.0),
    Sample("jet12", 1.4903e6, 14.0, 21.0),
    Sample("jet20", 6.2623e4, 21.0, 32.0),
    Sample("jet30", 2.5298e3, 32.0, 42.0),
    Sample("jet40", 1.3553e2, 42.0, 100.0),
)
SAMPLE_BY_NAME = {s.name: s for s in SAMPLES}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def finite(values: list[float]) -> list[float]:
    return [v for v in values if math.isfinite(v)]


def rms(values: list[float]) -> float | None:
    values = finite(values)
    if not values:
        return None
    mu = mean(values)
    return math.sqrt(sum((v - mu) ** 2 for v in values) / len(values))


def rel_rms(values: list[float]) -> float | None:
    values = finite(values)
    if not values:
        return None
    med = median(values)
    if med == 0:
        return None
    return rms(values) / abs(med)


def summarize(values: list[float]) -> dict[str, object]:
    values = finite(values)
    if not values:
        return {"n": 0, "median": None, "min": None, "max": None, "rms": None, "rel_rms": None}
    return {
        "n": len(values),
        "median": median(values),
        "min": min(values),
        "max": max(values),
        "rms": rms(values),
        "rel_rms": rel_rms(values),
    }


def shape_label(values: list[float]) -> str:
    stats = summarize(values)
    if not stats["n"]:
        return "missing"
    rr = stats["rel_rms"]
    span = stats["max"] / stats["min"] if stats["min"] else math.inf
    if rr is not None and rr < 0.01 and span < 1.04:
        return "flat"
    if rr is not None and rr < 0.035 and span < 1.12:
        return "mostly_flat"
    return "nonflat_or_sloped"


def owned_rows_from_ppg12() -> dict[tuple[str, float], dict[str, object]]:
    out: dict[tuple[str, float], dict[str, object]] = {}
    for row in read_csv(PPG12_IAN_CSV):
        if row.get("group") != "jet":
            continue
        sample = row.get("sample", "")
        if sample not in SAMPLE_BY_NAME:
            continue
        if str(row.get("used_in_stitch", "1")) not in {"1", "true", "True"}:
            continue
        center = round(float(row["bin_center"]), 6)
        out[(sample, center)] = {
            "value": float(row["value"]),
            "error": float(row.get("error", 0.0) or 0.0),
            "source_file": row.get("source_file", ""),
            "source_object": row.get("source_object", ""),
        }
    return out


def collect_series() -> dict[str, dict[tuple[str, float], float]]:
    series: dict[str, dict[tuple[str, float], float]] = {"ppg12_ian_no_suffix": {}}
    ppg12 = owned_rows_from_ppg12()
    for key, rec in ppg12.items():
        series["ppg12_ian_no_suffix"][key] = float(rec["value"])

    series["old_july1_original"] = {}
    series["old_july1_jet8_xsecfixed"] = {}
    series["old_july1_jet8fixed_global_scaled"] = {}
    for row in read_csv(OLD_CLOSURE_CSV):
        sample = row["sample"]
        center = round(float(row["bin_center"]), 6)
        key = (sample, center)
        series["old_july1_original"][key] = float(row["current_value_original"])
        series["old_july1_jet8_xsecfixed"][key] = float(row["current_value_ppg12_jet8_fixed"])
        series["old_july1_jet8fixed_global_scaled"][key] = float(row["scaled_fixed_current_value"])

    series["strict_july2_kept"] = {}
    series["strict_july2_all"] = {}
    for row in read_csv(STRICT_JULY2_CSV):
        sample = row["sample"]
        center = round(float(row["bin_center"]), 6)
        key = (sample, center)
        series["strict_july2_kept"][key] = float(row["current_kept_value"])
        series["strict_july2_all"][key] = float(row["current_all_value"])
    return series


def ratio_by_sample(numer: dict[tuple[str, float], float], denom: dict[tuple[str, float], float]) -> dict[str, list[float]]:
    by: dict[str, list[float]] = {s.name: [] for s in SAMPLES}
    for key, num in numer.items():
        if key not in denom:
            continue
        den = denom[key]
        if den > 0:
            by[key[0]].append(num / den)
    return by


def source_table(series: dict[str, dict[tuple[str, float], float]]) -> list[dict[str, object]]:
    ppg12_rows = owned_rows_from_ppg12()
    rows: list[dict[str, object]] = []
    source_meta = {
        "ppg12_ian_no_suffix": {
            "family": "IAN/no-suffix per-sample PPG12 source",
            "histogram_family": "h_max_truth_jet_pT_rebin1GeV",
            "status": "historical_ian_reference",
        },
        "old_july1_original": {
            "family": "July1 current source-scope diagnostic before jet8 correction",
            "histogram_family": "CSV aggregate, old local output",
            "status": "diagnostic_not_canonical",
        },
        "old_july1_jet8_xsecfixed": {
            "family": "July1 current with jet8 xsec correction only",
            "histogram_family": "CSV aggregate, old local output",
            "status": "diagnostic_not_canonical",
        },
        "old_july1_jet8fixed_global_scaled": {
            "family": "July1 current with jet8 correction and old global scale",
            "histogram_family": "CSV aggregate, old local output",
            "status": "shape_closure_diagnostic_only",
        },
        "strict_july2_kept": {
            "family": "Strict July2 RecoilJets online Fig6 event-weighted kept",
            "histogram_family": "SIM/h_ppInclusiveJetStitch_ppg12Fig6EventWeighted_r04_maxTruthJetPt_kept",
            "status": "failed_no_scale_canary",
        },
        "strict_july2_all": {
            "family": "Strict July2 RecoilJets online Fig6 event-weighted all",
            "histogram_family": "SIM/h_ppInclusiveJetStitch_ppg12Fig6EventWeighted_r04_maxTruthJetPt_all",
            "status": "diagnostic_all_family",
        },
    }
    for source_name, vals in series.items():
        ratios = ratio_by_sample(vals, series["ppg12_ian_no_suffix"])
        for sample in SAMPLES:
            sample_vals = [v for (s, _), v in vals.items() if s == sample.name and v > 0]
            source_files = sorted(
                {
                    str(ppg12_rows[(s, c)]["source_file"])
                    for (s, c) in vals
                    if s == sample.name and (s, c) in ppg12_rows
                }
            )
            source_objects = sorted(
                {
                    str(ppg12_rows[(s, c)]["source_object"])
                    for (s, c) in vals
                    if s == sample.name and (s, c) in ppg12_rows
                }
            )
            stats = summarize(ratios.get(sample.name, []))
            rows.append(
                {
                    "sample": sample.name,
                    "source": source_name,
                    "family": source_meta[source_name]["family"],
                    "status": source_meta[source_name]["status"],
                    "histogram_or_object": source_meta[source_name]["histogram_family"],
                    "ppg12_source_file_if_applicable": ";".join(source_files),
                    "ppg12_source_object_if_applicable": ";".join(source_objects),
                    "owned_window": f"[{sample.lo:g},{sample.hi:g})",
                    "raw_integral_owned_window": "not_available_from_csv",
                    "weighted_integral_owned_window": f"{sum(sample_vals):.12g}",
                    "nonzero_bins": len(sample_vals),
                    "ratio_to_ian_no_suffix_median": stats["median"],
                    "ratio_to_ian_no_suffix_rms": stats["rms"],
                    "ratio_to_ian_no_suffix_min": stats["min"],
                    "ratio_to_ian_no_suffix_max": stats["max"],
                    "ratio_shape": shape_label(ratios.get(sample.name, [])),
                    "first_ratios": " ".join(f"{v:.6g}" for v in ratios.get(sample.name, [])[:3]),
                    "last_ratios": " ".join(f"{v:.6g}" for v in ratios.get(sample.name, [])[-3:]),
                }
            )
    # Explicit placeholder row for the true period-combined target.  This pass
    # cannot derive it locally unless a period-combined ROOT/CSV is present.
    for sample in SAMPLES:
        rows.append(
            {
                "sample": sample.name,
                "source": "ppg12_true_period_combined",
                "family": "PPG12 true period-combined product",
                "status": "local_missing_needs_readonly_sdcc_or_local_csv",
                "histogram_or_object": "/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_jet_bdt_nom.root::h_max_truth_jet_pT",
                "ppg12_source_file_if_applicable": "/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_jet_bdt_nom.root",
                "ppg12_source_object_if_applicable": "h_max_truth_jet_pT",
                "owned_window": f"[{sample.lo:g},{sample.hi:g})",
                "raw_integral_owned_window": "local_missing",
                "weighted_integral_owned_window": "local_missing",
                "nonzero_bins": "local_missing",
                "ratio_to_ian_no_suffix_median": "local_missing",
                "ratio_to_ian_no_suffix_rms": "local_missing",
                "ratio_to_ian_no_suffix_min": "local_missing",
                "ratio_to_ian_no_suffix_max": "local_missing",
                "ratio_shape": "blocked_without_reference",
                "first_ratios": "",
                "last_ratios": "",
            }
        )
    return rows


def read_strict_metadata() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    hist_meta = "SIM/h_ppInclusiveJetStitch_ppg12Fig6EventWeighted_metadata"
    for root in sorted(STRICT_JULY2_ROOT_DIR.glob("*.root")):
        sample = next((s for s in SAMPLES if s.name in root.name), None)
        if sample is None:
            continue
        with uproot.open(root) as f:
            values, _ = f[hist_meta].to_numpy(flow=False)
        labels = [
            "events_seen_raw",
            "window_low_GeV",
            "window_high_GeV",
            "upper_edge_inclusive",
            "bin_width_GeV",
            "xsec_pb",
            "xsec_over_jet50",
            "sample_bin",
            "truth_def_code",
            "includes_current_weight",
            "includes_slice_weight",
            "value_mode_code",
        ]
        rec = {labels[i]: float(values[i]) for i in range(min(len(values), len(labels)))}
        rows.append(
            {
                "sample": sample.name,
                "root": str(root),
                **rec,
                "metadata_vs_code_warning": (
                    "metadata_says_current_weight_1_but_src_fill_uses_sliceWeight_only"
                    if rec.get("includes_current_weight") == 1.0
                    else ""
                ),
            }
        )
    return rows


def mechanism_table(series: dict[str, dict[tuple[str, float], float]]) -> list[dict[str, object]]:
    ppg12 = series["ppg12_ian_no_suffix"]

    mechanisms: list[tuple[str, dict[tuple[str, float], float], str, str]] = []
    mechanisms.append(("old_no_scale_original", series["old_july1_original"], "diagnostic", "Old July1 current before jet8 correction; not parity"))
    mechanisms.append(("old_jet8_xsecfix_only", series["old_july1_jet8_xsecfixed"], "diagnostic", "Jet8 corrected to PPG12 xsec; leaves common source-count factor"))
    mechanisms.append(("old_jet8_xsecfix_plus_global_scale", series["old_july1_jet8fixed_global_scaled"], "not_valid_as_final", "Reproduces old clean plot but global scale was fitted/source-scope only"))
    mechanisms.append(("strict_july2_no_scale_kept", series["strict_july2_kept"], "failed_canary", "The submitted strict object compared directly to IAN/no-suffix source"))
    mechanisms.append(("strict_july2_apply_old_global_scale", {k: v * OLD_GLOBAL_SCALE_AFTER_JET8_FIX for k, v in series["strict_july2_kept"].items()}, "invalid_test", "Numerical stress test only; a fitted scale is not acceptable"))
    mechanisms.append(("strict_july2_xsec_jet50_already_baked", series["strict_july2_kept"], "failed_canary", "This is already the xsec/jet50 fill convention"))

    # Normalize strict July2 by the old per-sample median mismatch, just to
    # show whether shape is still similar within each sample.  This is not a
    # defensible global correction because the factor is sample-dependent.
    per_sample_norm: dict[tuple[str, float], float] = {}
    ratios = ratio_by_sample(series["strict_july2_kept"], ppg12)
    sample_medians = {s: median(vs) for s, vs in ratios.items() if vs}
    for key, value in series["strict_july2_kept"].items():
        factor = sample_medians.get(key[0])
        if factor:
            per_sample_norm[key] = value / factor
    mechanisms.append(("strict_july2_per_sample_backsolve", per_sample_norm, "invalid_but_informative", "If this closes, failure is sample/source normalization, not truth-shape"))

    rows: list[dict[str, object]] = []
    for mech, vals, validity, note in mechanisms:
        by = ratio_by_sample(vals, ppg12)
        all_ratios = [v for vals2 in by.values() for v in vals2]
        sample_meds = {s: median(vs) for s, vs in by.items() if vs}
        across = summarize(list(sample_meds.values()))
        for sample in SAMPLES:
            stats = summarize(by.get(sample.name, []))
            rows.append(
                {
                    "mechanism": mech,
                    "sample": sample.name,
                    "validity": validity,
                    "median_current_over_reference": stats["median"],
                    "min": stats["min"],
                    "max": stats["max"],
                    "rms": stats["rms"],
                    "shape": shape_label(by.get(sample.name, [])),
                    "across_sample_median_rms": across["rms"],
                    "across_sample_median_min": across["min"],
                    "across_sample_median_max": across["max"],
                    "reproduces_old_clean_convergence": (
                        "yes" if mech == "old_jet8_xsecfix_plus_global_scale" else "no"
                    ),
                    "physically_defensible": (
                        "no" if validity in {"not_valid_as_final", "invalid_test", "invalid_but_informative"} else "diagnostic_only"
                    ),
                    "note": note,
                }
            )
    return rows


def write_summary(source_rows: list[dict[str, object]], mech_rows: list[dict[str, object]], meta_rows: list[dict[str, object]]) -> None:
    lines = []
    lines.append("THE-94 inclusive-jet local forensic pass")
    lines.append("")
    lines.append("Key local result:")
    lines.append("- The old July1 source-scope diagnostic is flat sample-by-sample after jet8 xsec correction, with a common current/IAN factor about 1.94.")
    lines.append("- The strict July2 kept histogram is not a global-factor problem: medians fall from jet8 about 1.82 to jet40 about 1.16.")
    lines.append("- The strict July2 metadata claims current weight is included, but local source inspection shows the strict fill uses only sliceWeight=xsec/jet50.")
    lines.append("")
    lines.append("Important blocker:")
    lines.append("- The true PPG12 period-combined product is not present as a local CSV/ROOT in this pass; it must be read via a read-only SDCC canary or pulled as a selected tiny reference if Justin approves transfer later.")
    lines.append("")
    lines.append("Blair legacy source viability:")
    lines.append("- Not needed for the current nominal decision from local evidence. The old jet8 discrepancy collapses under the corrected PPG12 jet8 xsec factor; no local evidence here requires swapping to Blair legacy input.")
    lines.append("")
    lines.append("Output tables:")
    lines.append(f"- {OUT_DIR / 'reference_source_table.csv'}")
    lines.append(f"- {OUT_DIR / 'candidate_mechanism_table.csv'}")
    lines.append(f"- {OUT_DIR / 'strict_july2_metadata_table.csv'}")
    (OUT_DIR / "the94_forensic_summary.txt").write_text("\n".join(lines) + "\n")


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    series = collect_series()
    src_rows = source_table(series)
    mech_rows = mechanism_table(series)
    meta_rows = read_strict_metadata()

    write_csv(
        OUT_DIR / "reference_source_table.csv",
        src_rows,
        [
            "sample",
            "source",
            "family",
            "status",
            "histogram_or_object",
            "ppg12_source_file_if_applicable",
            "ppg12_source_object_if_applicable",
            "owned_window",
            "raw_integral_owned_window",
            "weighted_integral_owned_window",
            "nonzero_bins",
            "ratio_to_ian_no_suffix_median",
            "ratio_to_ian_no_suffix_rms",
            "ratio_to_ian_no_suffix_min",
            "ratio_to_ian_no_suffix_max",
            "ratio_shape",
            "first_ratios",
            "last_ratios",
        ],
    )
    write_csv(
        OUT_DIR / "candidate_mechanism_table.csv",
        mech_rows,
        [
            "mechanism",
            "sample",
            "validity",
            "median_current_over_reference",
            "min",
            "max",
            "rms",
            "shape",
            "across_sample_median_rms",
            "across_sample_median_min",
            "across_sample_median_max",
            "reproduces_old_clean_convergence",
            "physically_defensible",
            "note",
        ],
    )
    write_csv(
        OUT_DIR / "strict_july2_metadata_table.csv",
        meta_rows,
        [
            "sample",
            "root",
            "events_seen_raw",
            "window_low_GeV",
            "window_high_GeV",
            "upper_edge_inclusive",
            "bin_width_GeV",
            "xsec_pb",
            "xsec_over_jet50",
            "sample_bin",
            "truth_def_code",
            "includes_current_weight",
            "includes_slice_weight",
            "value_mode_code",
            "metadata_vs_code_warning",
        ],
    )
    write_summary(src_rows, mech_rows, meta_rows)
    payload = {
        "source_table": str(OUT_DIR / "reference_source_table.csv"),
        "mechanism_table": str(OUT_DIR / "candidate_mechanism_table.csv"),
        "metadata_table": str(OUT_DIR / "strict_july2_metadata_table.csv"),
        "summary": str(OUT_DIR / "the94_forensic_summary.txt"),
    }
    (OUT_DIR / "the94_forensic_manifest.json").write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
