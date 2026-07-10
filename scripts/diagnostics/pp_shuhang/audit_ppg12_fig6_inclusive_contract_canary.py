#!/usr/bin/env python3
"""PPG12 Fig.6 inclusive-jet contract canary.

This is a pre-resubmission diagnostic. It compares the current RecoilJets
inclusive-jet Fig.6 event-weighted histograms against the exact PPG12
efficiency-tool source convention without applying any fitted/global scale.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path
from statistics import median

import numpy as np

try:
    import uproot
except ModuleNotFoundError as exc:  # pragma: no cover - environment guard
    raise SystemExit(
        "This canary requires uproot. Run with:\n"
        "/Users/patsfan753/Desktop/analysis/env/bin/python "
        "scripts/diagnostics/pp_shuhang/audit_ppg12_fig6_inclusive_contract_canary.py"
    ) from exc


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
DEFAULT_TAG = "the76_ppg12_fig6_inclusive_strict_20260702_014757"
DEFAULT_REFERENCE_TAG = "the76_ppg12_parity_full_20260701_003024"
DEFAULT_CURRENT_ROOT_DIR = (
    REPO
    / "dataOutput/ppg12Parity"
    / DEFAULT_TAG
    / "final_roots/inclusivejet"
)
DEFAULT_PPG12_SOURCE_CSV = (
    REPO
    / "dataOutput/ppg12Parity"
    / DEFAULT_REFERENCE_TAG
    / "reference_audit/sdcc_pull/live_ppg12_ian_source_points.csv"
)
DEFAULT_OUT_DIR = (
    REPO
    / "dataOutput/ppg12Parity"
    / DEFAULT_TAG
    / "contract_canary/fig6_inclusive"
)


JET50_XSEC_PB = 7.3113


@dataclass(frozen=True)
class JetSample:
    sample: str
    root_token: str
    xsec_pb: float
    window_low: float
    window_high: float
    sample_bin: int

    @property
    def xsec_over_jet50(self) -> float:
        return self.xsec_pb / JET50_XSEC_PB


SAMPLES: tuple[JetSample, ...] = (
    JetSample("jet8", "jet8", 1.15e7, 9.0, 14.0, 2),
    JetSample("jet12", "jet12", 1.4903e6, 14.0, 21.0, 3),
    JetSample("jet20", "jet20", 6.2623e4, 21.0, 32.0, 4),
    JetSample("jet30", "jet30", 2.5298e3, 32.0, 42.0, 5),
    JetSample("jet40", "jet40", 1.3553e2, 42.0, 100.0, 6),
)

SAMPLE_BY_NAME = {s.sample: s for s in SAMPLES}

HIST_ALL = "SIM/h_ppInclusiveJetStitch_ppg12Fig6EventWeighted_r04_maxTruthJetPt_all"
HIST_KEPT = "SIM/h_ppInclusiveJetStitch_ppg12Fig6EventWeighted_r04_maxTruthJetPt_kept"
HIST_REJECTED = "SIM/h_ppInclusiveJetStitch_ppg12Fig6EventWeighted_r04_maxTruthJetPt_rejected"
HIST_META = "SIM/h_ppInclusiveJetStitch_ppg12Fig6EventWeighted_metadata"


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(path)
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def detect_sample(path: Path) -> JetSample | None:
    name = path.name.lower()
    for sample in SAMPLES:
        if sample.root_token in name:
            return sample
    return None


def load_ppg12_source_points(path: Path) -> dict[tuple[str, float], dict[str, float]]:
    out: dict[tuple[str, float], dict[str, float]] = {}
    for row in read_csv(path):
        if row.get("group") != "jet":
            continue
        sample = row.get("sample", "")
        if sample not in SAMPLE_BY_NAME:
            continue
        used = row.get("ppg12_bin_center_window", row.get("used_in_stitch", "1"))
        if str(used) not in {"1", "True", "true"}:
            continue
        x = round(float(row["bin_center"]), 6)
        value_key = "value" if "value" in row else "ppg12_source_value"
        error_key = "error" if "error" in row else "ppg12_source_error"
        value = float(row[value_key])
        if value <= 0:
            continue
        out[(sample, x)] = {
            "value": value,
            "error": float(row.get(error_key, "0") or 0.0),
        }
    if not out:
        raise RuntimeError(f"No PPG12 Fig.6 jet source points loaded from {path}")
    return out


def hist_to_numpy(root_file: uproot.ReadOnlyDirectory, hist_name: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if hist_name not in root_file:
        raise KeyError(f"missing histogram {hist_name}")
    hist = root_file[hist_name]
    values, edges = hist.to_numpy(flow=False)
    variances = hist.variances(flow=False)
    if variances is None:
        variances = np.clip(values, 0.0, None)
    return values.astype(float), variances.astype(float), edges.astype(float)


def metadata_values(root_file: uproot.ReadOnlyDirectory) -> dict[str, float] | None:
    if HIST_META not in root_file:
        return None
    values, _ = root_file[HIST_META].to_numpy(flow=False)
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
    return {label: float(values[i]) for i, label in enumerate(labels[: len(values)])}


def in_ppg12_window(centers: np.ndarray, sample: JetSample) -> np.ndarray:
    # PPG12 code rejects maxjetpT > upper and maxjetpT < lower, so an exact
    # upper-edge event is accepted. The 1 GeV comparison uses bin centers, so
    # this equals the displayed Fig.6 row selection.
    return (centers >= sample.window_low) & (centers <= sample.window_high)


def summarize_ratio(values: list[float]) -> dict[str, object]:
    if not values:
        return {
            "n_bins": 0,
            "median": None,
            "min": None,
            "max": None,
            "mean": None,
            "closure_status": "missing",
        }
    mean = sum(values) / len(values)
    med = median(values)
    status = "pass_2pct" if 0.98 <= med <= 1.02 and min(values) >= 0.95 and max(values) <= 1.05 else "fail"
    return {
        "n_bins": len(values),
        "median": med,
        "min": min(values),
        "max": max(values),
        "mean": mean,
        "closure_status": status,
    }


def audit_current_roots(current_root_dir: Path, source_points: dict[tuple[str, float], dict[str, float]]) -> tuple[list[dict[str, object]], list[dict[str, object]], list[str]]:
    rows: list[dict[str, object]] = []
    summaries: list[dict[str, object]] = []
    warnings: list[str] = []
    root_paths = sorted(current_root_dir.glob("*.root"))
    if not root_paths:
        raise FileNotFoundError(f"No ROOT files found in {current_root_dir}")

    found_samples: set[str] = set()
    for root_path in root_paths:
        sample = detect_sample(root_path)
        if sample is None:
            warnings.append(f"skipped unrecognized ROOT file: {root_path}")
            continue
        found_samples.add(sample.sample)
        with uproot.open(root_path) as f:
            values_all, variances_all, edges = hist_to_numpy(f, HIST_ALL)
            values_kept, variances_kept, kept_edges = hist_to_numpy(f, HIST_KEPT)
            values_rejected, _, rejected_edges = hist_to_numpy(f, HIST_REJECTED)
            if not np.allclose(edges, kept_edges) or not np.allclose(edges, rejected_edges):
                raise RuntimeError(f"{root_path}: all/kept/rejected bin edges differ")
            meta = metadata_values(f)

        centers = 0.5 * (edges[:-1] + edges[1:])
        mask = in_ppg12_window(centers, sample)
        sample_ratios_kept: list[float] = []
        sample_ratios_all: list[float] = []
        for i, center in enumerate(centers):
            if not mask[i]:
                continue
            key = (sample.sample, round(float(center), 6))
            if key not in source_points:
                continue
            ppg12_value = source_points[key]["value"]
            ppg12_error = source_points[key]["error"]
            current_kept = float(values_kept[i])
            current_all = float(values_all[i])
            current_rejected = float(values_rejected[i])
            current_kept_error = math.sqrt(max(float(variances_kept[i]), 0.0))
            current_all_error = math.sqrt(max(float(variances_all[i]), 0.0))
            ratio_kept = current_kept / ppg12_value if ppg12_value else math.nan
            ratio_all = current_all / ppg12_value if ppg12_value else math.nan
            sample_ratios_kept.append(ratio_kept)
            sample_ratios_all.append(ratio_all)
            rows.append(
                {
                    "sample": sample.sample,
                    "bin_low": f"{edges[i]:.8g}",
                    "bin_high": f"{edges[i + 1]:.8g}",
                    "bin_center": f"{center:.8g}",
                    "ppg12_source_value": f"{ppg12_value:.12g}",
                    "ppg12_source_error": f"{ppg12_error:.12g}",
                    "current_kept_value": f"{current_kept:.12g}",
                    "current_kept_error": f"{current_kept_error:.12g}",
                    "current_all_value": f"{current_all:.12g}",
                    "current_all_error": f"{current_all_error:.12g}",
                    "current_rejected_value": f"{current_rejected:.12g}",
                    "current_kept_over_ppg12": f"{ratio_kept:.12g}",
                    "current_all_over_ppg12": f"{ratio_all:.12g}",
                    "all_minus_kept": f"{current_all - current_kept:.12g}",
                    "xsec_pb_contract": f"{sample.xsec_pb:.12g}",
                    "xsec_over_jet50_contract": f"{sample.xsec_over_jet50:.12g}",
                    "window_low_contract": f"{sample.window_low:.8g}",
                    "window_high_contract": f"{sample.window_high:.8g}",
                    "current_root": str(root_path),
                    "parity_object": HIST_KEPT,
                    "diagnostic_all_object": HIST_ALL,
                }
            )

        meta_warning = None
        if meta:
            if not math.isclose(meta.get("window_low_GeV", sample.window_low), sample.window_low, rel_tol=0.0, abs_tol=1e-6):
                meta_warning = "merged_scalar_metadata_not_single_value"
            elif not math.isclose(meta.get("xsec_over_jet50", sample.xsec_over_jet50), sample.xsec_over_jet50, rel_tol=1e-6, abs_tol=1e-9):
                meta_warning = "merged_scalar_metadata_not_single_value"
        summaries.append(
            {
                "sample": sample.sample,
                "root": str(root_path),
                "events_seen_raw_metadata": meta.get("events_seen_raw") if meta else None,
                "kept_window_integral": float(np.sum(values_kept[mask])),
                "all_window_integral": float(np.sum(values_all[mask])),
                "rejected_window_integral": float(np.sum(values_rejected[mask])),
                "all_minus_kept_window_integral": float(np.sum(values_all[mask] - values_kept[mask])),
                "ratio_kept_over_ppg12": summarize_ratio(sample_ratios_kept),
                "ratio_all_over_ppg12": summarize_ratio(sample_ratios_all),
                "metadata_scalar_warning": meta_warning,
                "metadata_raw": meta,
                "contract": {
                    "hist_fill_weight": "jetXcross/jet50cross",
                    "event_weight": "mix_weight * period_lumi/lumi_target * truth_vertex_weight(z_hard) [* truth_vertex_weight(z_mb) for DI]",
                    "window": [sample.window_low, sample.window_high],
                    "xsec_pb": sample.xsec_pb,
                    "xsec_over_jet50": sample.xsec_over_jet50,
                    "parity_histogram": HIST_KEPT,
                    "no_global_scale": True,
                },
            }
        )

    missing = sorted(set(SAMPLE_BY_NAME) - found_samples)
    if missing:
        warnings.append(f"missing current ROOT samples: {', '.join(missing)}")
    return rows, summaries, warnings


def build_contract_map() -> dict[str, object]:
    return {
        "ppg12_top_down_map": [
            {
                "step": "input",
                "ppg12_file": "ppg12codeGit/efficiencytool/config_bdt_nom_{0rad,1p5mrad}.yaml",
                "ppg12_contract": "slimtree inputs from /sphenix/user/shuhangli/ppg12/FunWithxgboost/...; truth_vertex_reweight_on=1; per-period lumi values are pre-scaled to lumi_target=64.3718",
                "our_contract": "RecoilJets Fun4All reads the pp SIM DST/G4 lanes, rebuilds relevant objects, then writes the Fig.6 histogram directly.",
            },
            {
                "step": "sample constants",
                "ppg12_file": "ppg12codeGit/efficiencytool/CrossSectionWeights.h",
                "ppg12_contract": "jet8/12/20/30/40 windows 9-14, 14-21, 21-32, 32-42, 42-100; weights are jetXcross/jet50cross.",
                "our_contract": "src/RecoilJets.cc ppg12InclusiveJetSliceWindow and ppg12InclusiveJetSliceXsecPb must match these values.",
            },
            {
                "step": "single/double blend",
                "ppg12_file": "ppg12codeGit/efficiencytool/oneforall_tree_double.sh",
                "ppg12_contract": "single and double are filled separately with mix_weight=(1-fDI) and fDI, then hadded.",
                "our_contract": "strict Fig.6 campaigns must preserve period/component weights before hadd; generic finalStitch scaling must not reweight this object.",
            },
            {
                "step": "event weight",
                "ppg12_file": "ppg12codeGit/efficiencytool/RecoEffCalculator_TTreeReader.C",
                "ppg12_contract": "weight=(jetXcross/jet50cross)*mix_weight*(lumi/lumi_target)*truthVertexWeight(z_hard[, z_mb]). reco vertex reweight is bypassed when truth_vertex_reweight_on=1.",
                "our_contract": "RJMCWeighting::CurrentWeight must be set to mix*lumi*truth-vtx before fill, while the explicit TH1 fill weight is jetXcross/jet50cross.",
            },
            {
                "step": "truth observable and gate",
                "ppg12_file": "ppg12codeGit/efficiencytool/RecoEffCalculator_TTreeReader.C",
                "ppg12_contract": "max truth jet pT is computed, events outside the owned window are skipped, then h_max_truth_jet_pT is filled.",
                "our_contract": "Use max R=0.4 truth jet pT; parity object is *_kept, not *_all.",
            },
            {
                "step": "plot/reference",
                "ppg12_file": "ppg12codeGit/plotting/plot_sample_combining_background.C",
                "ppg12_contract": "Fig.6 source uses MC_efficiency_jet*_bdt_nom.root::h_max_truth_jet_pT, rebinned to 1 GeV weighted counts.",
                "our_contract": "Compare current/PPG12 directly, no fitted/global scale.",
            },
        ],
        "historical_reference_warning": "PPG12 reports document an old jet8 SI source exception in the historical run28 tree; current ppg12codeGit jet8 Fun4All macro is corrected to g4hits.list. Keep historical IAN-source reproduction separate from corrected PPG12 policy.",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--current-root-dir", type=Path, default=DEFAULT_CURRENT_ROOT_DIR)
    parser.add_argument("--ppg12-source-csv", type=Path, default=DEFAULT_PPG12_SOURCE_CSV)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()

    source_points = load_ppg12_source_points(args.ppg12_source_csv)
    rows, summaries, warnings = audit_current_roots(args.current_root_dir, source_points)

    out_dir = args.out_dir
    points_csv = out_dir / "fig6_inclusive_current_kept_vs_ppg12_source_points.csv"
    summary_json = out_dir / "fig6_inclusive_contract_canary_summary.json"
    note_txt = out_dir / "fig6_inclusive_contract_canary_note.txt"

    fields = [
        "sample",
        "bin_low",
        "bin_high",
        "bin_center",
        "ppg12_source_value",
        "ppg12_source_error",
        "current_kept_value",
        "current_kept_error",
        "current_all_value",
        "current_all_error",
        "current_rejected_value",
        "current_kept_over_ppg12",
        "current_all_over_ppg12",
        "all_minus_kept",
        "xsec_pb_contract",
        "xsec_over_jet50_contract",
        "window_low_contract",
        "window_high_contract",
        "current_root",
        "parity_object",
        "diagnostic_all_object",
    ]
    write_csv(points_csv, rows, fields)

    all_ratios = [
        float(row["current_kept_over_ppg12"])
        for row in rows
        if math.isfinite(float(row["current_kept_over_ppg12"]))
    ]
    payload = {
        "schema": "THE76_PPG12_FIG6_INCLUSIVE_CONTRACT_CANARY_V1",
        "scale_policy": "none",
        "current_root_dir": str(args.current_root_dir),
        "ppg12_source_csv": str(args.ppg12_source_csv),
        "points_csv": str(points_csv),
        "contract_map": build_contract_map(),
        "global_kept_over_ppg12": summarize_ratio(all_ratios),
        "sample_summaries": summaries,
        "warnings": warnings,
        "component_decomposition_status": "not_available_from_local_final_roots; requires period/component ROOTs or per-event canary output before another submission",
        "decision": "pass" if summarize_ratio(all_ratios)["closure_status"] == "pass_2pct" else "fail_stop_before_resubmission",
    }
    write_json(summary_json, payload)

    lines = [
        "PPG12 Fig.6 inclusive-jet contract canary",
        f"decision: {payload['decision']}",
        f"scale_policy: {payload['scale_policy']}",
        f"current_root_dir: {args.current_root_dir}",
        f"ppg12_source_csv: {args.ppg12_source_csv}",
        "",
        "Sample medians, current kept / PPG12 source:",
    ]
    for summary in summaries:
        ratio = summary["ratio_kept_over_ppg12"]
        lines.append(
            f"  {summary['sample']}: median={ratio['median']} "
            f"min={ratio['min']} max={ratio['max']} status={ratio['closure_status']}"
        )
    lines.extend(
        [
            "",
            "Interpretation:",
            "  The exact parity object is *_kept. The current local final roots do not close to PPG12 Fig.6 without a scale.",
            "  The failure is sample-dependent, so it is not a single harmless plotting-scale choice.",
            "  Component decomposition is still required before another Condor resubmission.",
        ]
    )
    if warnings:
        lines.append("")
        lines.append("Warnings:")
        lines.extend(f"  {w}" for w in warnings)
    note_txt.parent.mkdir(parents=True, exist_ok=True)
    note_txt.write_text("\n".join(lines) + "\n")

    print(note_txt)
    print(summary_json)
    print(points_csv)
    for summary in summaries:
        ratio = summary["ratio_kept_over_ppg12"]
        print(
            f"{summary['sample']} median={ratio['median']} "
            f"min={ratio['min']} max={ratio['max']} status={ratio['closure_status']}"
        )
    print(f"decision={payload['decision']}")
    return 0 if payload["decision"] == "pass" else 2


if __name__ == "__main__":
    raise SystemExit(main())
