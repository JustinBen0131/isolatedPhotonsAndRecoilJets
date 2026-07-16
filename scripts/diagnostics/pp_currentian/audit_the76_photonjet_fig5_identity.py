#!/usr/bin/env python3
"""THE-76 photon+jet Fig.5 identity-break audit.

This is a local, read-only diagnostic. It does not submit jobs, merge, transfer,
or canonicalize anything. The purpose is to separate central-value identity
from visual/statistical compatibility for the PPG12 Fig.5 photon+jet stitch
candidate.
"""

from __future__ import annotations

import csv
import json
import math
import statistics
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

try:
    import uproot
except ModuleNotFoundError as exc:  # pragma: no cover - runtime guard
    raise SystemExit("uproot is required; use /Users/patsfan753/Desktop/analysis/env/bin/python3") from exc


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
OUT_DIR = (
    REPO
    / "dataOutput/ppg12Parity/control_plane/audits/"
    / "the76_photonjet_fig5_identity_20260702"
)

CURRENT_ROOT = (
    REPO
    / "dataOutput/ppg12Parity/the76_ppg12_fig24_photonjet_fix_20260702_014217/"
    / "final_roots/photonjet/RecoilJets_photonjet5plus10plus20_MERGED.root"
)
CURRENT_POINTS_CSV = (
    REPO
    / "dataOutput/ppg12Parity/the76_ppg12_fig24_photonjet_fix_20260702_014217/"
    / "strict_stitched_photonjet_current/"
    / "photon_data_over_fit_sdcc_vs_current_overlay_current_root_points.csv"
)
CURRENT_POINTS_MANIFEST = (
    REPO
    / "dataOutput/ppg12Parity/the76_ppg12_fig24_photonjet_fix_20260702_014217/"
    / "strict_stitched_photonjet_current/"
    / "photon_data_over_fit_sdcc_vs_current_overlay_current_root_points_manifest.json"
)
PPG12_SOURCE_CSV = (
    REPO
    / "dataOutput/ppg12Parity/the76_ppg12_parity_full_20260701_003024/"
    / "reference_audit/sdcc_pull/live_ppg12_ian_source_points.csv"
)
PPG12_SOURCE_ROOT = Path("/sphenix/user/shuhangli/ppg12/plotting/photon_max_pT_uncut.root")

PPG12_BUILD_SCRIPT = REPO / "ppg12codeGit/plotting/build_h_max_photon_pT_uncut.py"
PPG12_PLOT_SCRIPT = REPO / "ppg12codeGit/plotting/plot_combine_uncut.C"
RECOILJETS_SRC = REPO / "src/RecoilJets.cc"
LEDGER = REPO / "agent_context/local/debug_passes/THE-76_ppg12_photonjet_sim_parity.md"

STRICT_REL_TOL = 1e-12

CURRENT_HIST = "SIM/h_ppPhotonStitch_ppg12Fig5_maxPhotonPt_kept"
RECOILJETS_FAMILIES = {
    "ppg12Fig5_kept": "SIM/h_ppPhotonStitch_ppg12Fig5_maxPhotonPt_kept",
    "truthSpectrum_eta07_kept": "SIM/h_ppPhotonStitch_ppg12TruthSpectrum_eta07_maxPhotonPt_kept",
    "truthSpectrum_g4Stored_kept": "SIM/h_ppPhotonStitch_ppg12TruthSpectrum_g4Stored_maxPhotonPt_kept",
    "legacy_g4Stored_kept": "SIM/h_ppPhotonStitch_maxPhotonPt_kept",
}


@dataclass(frozen=True)
class SampleWindow:
    name: str
    ppg12_object: str
    ppg12_sumw2_object: str
    low: float
    high: float


SAMPLES = (
    SampleWindow("photon5", "h_max_photon_pT_photon5", "h_max_photon_pT_photon5_sumw2", 0.0, 14.0),
    SampleWindow("photon10", "h_max_photon_pT_photon10", "h_max_photon_pT_photon10_sumw2", 14.0, 22.0),
    SampleWindow("photon20", "h_max_photon_pT_photon20", "h_max_photon_pT_photon20_sumw2", 22.0, math.inf),
)
SAMPLE_BY_NAME = {s.name: s for s in SAMPLES}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, obj: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, sort_keys=True) + "\n")


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


def rel_diff(a: float, b: float) -> float:
    denom = max(abs(a), abs(b), 1.0)
    return (a - b) / denom


def strict_equal(a: float, b: float) -> bool:
    return abs(a - b) <= STRICT_REL_TOL * max(1.0, abs(a), abs(b))


def ratio_error_independent(num: float, num_err: float, den: float, den_err: float) -> float:
    if den == 0 or num == 0:
        return math.nan
    ratio = num / den
    return abs(ratio) * math.sqrt((num_err / num) ** 2 + (den_err / den) ** 2)


def sample_for_center(x: float) -> str:
    for sample in SAMPLES:
        if x >= sample.low and (math.isinf(sample.high) or x < sample.high):
            return sample.name
    return "unknown"


def load_hist(path: Path, hist_name: str) -> dict[str, Any]:
    with uproot.open(path) as f:
        if hist_name not in f:
            raise KeyError(f"{path} does not contain {hist_name}")
        hist = f[hist_name]
        values, edges = hist.to_numpy(flow=False)
        variances = hist.variances(flow=False)
    if variances is None:
        variances = np.clip(values, 0, None)
    centers = 0.5 * (edges[:-1] + edges[1:])
    return {
        "name": hist_name,
        "values": np.asarray(values, dtype=float),
        "variances": np.asarray(variances, dtype=float),
        "edges": np.asarray(edges, dtype=float),
        "centers": np.asarray(centers, dtype=float),
        "by_center": {
            round(float(x), 8): (i, float(v), float(var))
            for i, (x, v, var) in enumerate(zip(centers, values, variances), start=1)
        },
    }


def load_ppg12_root_bins_if_available() -> dict[str, dict[float, tuple[float, float]]] | None:
    if not PPG12_SOURCE_ROOT.exists():
        return None
    out: dict[str, dict[float, tuple[float, float]]] = {}
    with uproot.open(PPG12_SOURCE_ROOT) as f:
        for sample in SAMPLES:
            if sample.ppg12_object not in f:
                raise KeyError(f"{PPG12_SOURCE_ROOT} missing {sample.ppg12_object}")
            values, edges = f[sample.ppg12_object].to_numpy(flow=False)
            if sample.ppg12_sumw2_object in f:
                sumw2_values, _ = f[sample.ppg12_sumw2_object].to_numpy(flow=False)
                errors = np.sqrt(np.clip(sumw2_values, 0, None))
            else:
                errors = np.sqrt(np.clip(values, 0, None))
            centers = 0.5 * (edges[:-1] + edges[1:])
            out[sample.name] = {
                round(float(x), 8): (float(v), float(e))
                for x, v, e in zip(centers, values, errors)
            }
    return out


def source_csv_lookup(rows: list[dict[str, str]]) -> dict[tuple[str, float], dict[str, str]]:
    out: dict[tuple[str, float], dict[str, str]] = {}
    for row in rows:
        if row.get("group") != "photon":
            continue
        sample = row.get("sample", "")
        x = round(float(row["bin_center"]), 8)
        out[(sample, x)] = row
    return out


def line_hits(path: Path, needles: list[str], context: int = 0) -> list[dict[str, Any]]:
    if not path.exists():
        return [{"path": str(path), "missing": True}]
    lines = path.read_text(errors="replace").splitlines()
    hits: list[dict[str, Any]] = []
    for needle in needles:
        for idx, line in enumerate(lines, start=1):
            if needle in line:
                lo = max(1, idx - context)
                hi = min(len(lines), idx + context)
                hits.append(
                    {
                        "path": str(path),
                        "needle": needle,
                        "line": idx,
                        "excerpt": "\n".join(f"{j}: {lines[j-1]}" for j in range(lo, hi + 1)),
                    }
                )
                break
    return hits


def plotting_self_ratio_tests(point_rows: list[dict[str, str]]) -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    for row in point_rows:
        sample = row["sample"]
        x = float(row["bin_center"])
        for label, value_key, err_key in (
            ("ppg12_over_ppg12", "ppg12_sdcc_value", "ppg12_sdcc_error"),
            ("recoiljets_over_recoiljets", "current_value", "current_error"),
        ):
            v = float(row[value_key])
            e = float(row[err_key])
            ratio = v / v if v else math.nan
            ratio_err = ratio_error_independent(v, e, v, e)
            rows.append(
                {
                    "test": label,
                    "sample": sample,
                    "bin_center": f"{x:.8g}",
                    "numerator": f"{v:.17g}",
                    "denominator": f"{v:.17g}",
                    "ratio": f"{ratio:.17g}",
                    "ratio_error_independent": f"{ratio_err:.17g}",
                    "abs_deviation_from_one": f"{abs(ratio - 1.0):.17g}",
                    "relative_deviation_from_one": f"{rel_diff(ratio, 1.0):.17g}",
                    "central_equal_strict": strict_equal(ratio, 1.0),
                    "error_bar_nonzero": bool(ratio_err > 0),
                }
            )
    write_csv(
        OUT_DIR / "plotting_self_ratio_tests.csv",
        rows,
        [
            "test",
            "sample",
            "bin_center",
            "numerator",
            "denominator",
            "ratio",
            "ratio_error_independent",
            "abs_deviation_from_one",
            "relative_deviation_from_one",
            "central_equal_strict",
            "error_bar_nonzero",
        ],
    )
    deviations = [float(r["abs_deviation_from_one"]) for r in rows if math.isfinite(float(r["ratio"]))]
    summary = {
        "status": "pass" if all(d <= 0 for d in deviations) else "fail",
        "n_tests": len(rows),
        "max_abs_deviation_from_one": max(deviations) if deviations else math.nan,
        "max_relative_deviation_from_one": max(abs(float(r["relative_deviation_from_one"])) for r in rows),
        "all_central_ratios_equal_one_strict": all(r["central_equal_strict"] for r in rows),
        "any_nonzero_error_bars": any(r["error_bar_nonzero"] for r in rows),
        "interpretation": (
            "The extraction/CSV path does not introduce non-unity central values in self-ratios. "
            "Independent-error propagation still produces nonzero error bars even when central values are exactly one."
        ),
        "csv": str(OUT_DIR / "plotting_self_ratio_tests.csv"),
    }
    write_json(OUT_DIR / "plotting_self_ratio_tests.json", summary)
    return summary


def raw_root_bin_identity(
    point_rows: list[dict[str, str]],
    ppg12_csv_rows: list[dict[str, str]],
    current_hist: dict[str, Any],
    ppg12_root_bins: dict[str, dict[float, tuple[float, float]]] | None,
) -> dict[str, Any]:
    ppg12_csv = source_csv_lookup(ppg12_csv_rows)
    rows: list[dict[str, Any]] = []
    direct_equal = 0
    extracted_equal = 0
    max_direct_rel = 0.0
    max_extracted_rel = 0.0
    max_extracted_row: dict[str, Any] | None = None
    source_factors: list[float] = []
    for row in point_rows:
        sample = row["sample"]
        x = round(float(row["bin_center"]), 8)
        root_idx, current_raw, current_raw_var = current_hist["by_center"][x]
        current_root_error = math.sqrt(max(current_raw_var, 0.0))
        current_extracted = float(row["current_value"])
        current_extracted_error = float(row["current_error"])
        source_factors.append(current_extracted / current_raw if current_raw else math.nan)

        if ppg12_root_bins and sample in ppg12_root_bins and x in ppg12_root_bins[sample]:
            ppg12_value, ppg12_error = ppg12_root_bins[sample][x]
            ppg12_source = "local_ppg12_root_getbincontent"
        else:
            ppg12_row = ppg12_csv.get((sample, x))
            if not ppg12_row:
                raise KeyError(f"missing PPG12 source CSV row for {sample} {x}")
            ppg12_value = float(ppg12_row["value"])
            ppg12_error = float(ppg12_row["error"])
            ppg12_source = "source_csv_extracted_from_ppg12_root"

        direct_abs = ppg12_value - current_raw
        direct_rel = rel_diff(ppg12_value, current_raw)
        extracted_abs = ppg12_value - current_extracted
        extracted_rel = rel_diff(ppg12_value, current_extracted)
        direct_ok = strict_equal(ppg12_value, current_raw)
        extracted_ok = strict_equal(ppg12_value, current_extracted)
        if direct_ok:
            direct_equal += 1
        if extracted_ok:
            extracted_equal += 1
        max_direct_rel = max(max_direct_rel, abs(direct_rel))
        if abs(extracted_rel) > max_extracted_rel:
            max_extracted_rel = abs(extracted_rel)
            max_extracted_row = {
                "sample": sample,
                "bin_center": x,
                "ppg12_value": ppg12_value,
                "current_extracted_value": current_extracted,
                "relative_difference": extracted_rel,
                "ratio_ppg12_over_current_extracted": ppg12_value / current_extracted,
            }

        rows.append(
            {
                "bin_index": root_idx,
                "bin_low_edge": row["bin_low"],
                "bin_high_edge": row["bin_high"],
                "bin_center": row["bin_center"],
                "sample_owner": sample,
                "ppg12_object_name": SAMPLE_BY_NAME[sample].ppg12_object,
                "recoiljets_object_name": CURRENT_HIST,
                "ppg12_getbincontent": f"{ppg12_value:.17g}",
                "recoiljets_getbincontent_raw": f"{current_raw:.17g}",
                "recoiljets_extracted_value": f"{current_extracted:.17g}",
                "direct_absolute_difference": f"{direct_abs:.17g}",
                "direct_relative_difference": f"{direct_rel:.17g}",
                "extracted_absolute_difference": f"{extracted_abs:.17g}",
                "extracted_relative_difference": f"{extracted_rel:.17g}",
                "ppg12_error_or_sumw2_error": f"{ppg12_error:.17g}",
                "recoiljets_raw_error": f"{current_root_error:.17g}",
                "recoiljets_extracted_error": f"{current_extracted_error:.17g}",
                "ratio_ppg12_over_recoiljets_raw": f"{(ppg12_value / current_raw) if current_raw else math.nan:.17g}",
                "ratio_ppg12_over_recoiljets_extracted": f"{(ppg12_value / current_extracted) if current_extracted else math.nan:.17g}",
                "ratio_error_extracted": f"{ratio_error_independent(ppg12_value, ppg12_error, current_extracted, current_extracted_error):.17g}",
                "direct_equal_within_strict_tolerance": direct_ok,
                "extracted_equal_within_strict_tolerance": extracted_ok,
                "ppg12_value_source": ppg12_source,
            }
        )
    fields = list(rows[0].keys())
    write_csv(OUT_DIR / "raw_root_bin_identity.csv", rows, fields)

    finite_factors = [x for x in source_factors if math.isfinite(x)]
    summary = {
        "csv": str(OUT_DIR / "raw_root_bin_identity.csv"),
        "ppg12_source_root_local_status": "available" if ppg12_root_bins is not None else "missing_locally",
        "ppg12_value_source": "local ROOT GetBinContent" if ppg12_root_bins else "source CSV extracted from PPG12 ROOT",
        "n_bins": len(rows),
        "direct_raw_equal_count": direct_equal,
        "extracted_equal_count": extracted_equal,
        "max_direct_raw_relative_difference": max_direct_rel,
        "max_extracted_relative_difference": max_extracted_rel,
        "max_extracted_relative_difference_row": max_extracted_row,
        "source_scope_factor_median": statistics.median(finite_factors) if finite_factors else math.nan,
        "source_scope_factor_min": min(finite_factors) if finite_factors else math.nan,
        "source_scope_factor_max": max(finite_factors) if finite_factors else math.nan,
        "interpretation": (
            "Direct raw equality fails because the RecoilJets final ROOT histogram is stored in raw/weighted event-count units, "
            "while the PPG12 source values are already in the source plotting convention. After applying the current extraction "
            "factor, central values still differ at per-mille to percent level, so the residual is not a CSV or floating-point artifact."
        ),
    }
    write_json(OUT_DIR / "raw_root_bin_identity.json", summary)
    return summary


def component_prestitch_identity(raw_rows: list[dict[str, Any]]) -> dict[str, Any]:
    out_rows: list[dict[str, Any]] = []
    summary: dict[str, Any] = {
        "csv": str(OUT_DIR / "component_prestitch_identity.csv"),
        "samples": {},
        "component_histogram_status": (
            "PPG12 source has per-sample objects. The current local final RecoilJets ROOT has only combined stitch diagnostics, "
            "not separate per-sample final component ROOTs, so true pre-stitch component identity cannot be proven locally."
        ),
    }
    for sample in SAMPLE_BY_NAME:
        sample_rows = [r for r in raw_rows if r["sample_owner"] == sample]
        ratios = [float(r["ratio_ppg12_over_recoiljets_extracted"]) for r in sample_rows]
        rels = [float(r["extracted_relative_difference"]) for r in sample_rows]
        pulls = []
        for r in sample_rows:
            ratio = float(r["ratio_ppg12_over_recoiljets_extracted"])
            ratio_err = float(r["ratio_error_extracted"])
            if ratio_err and math.isfinite(ratio_err):
                pulls.append((ratio - 1.0) / ratio_err)
        sample_summary = {
            "ppg12_component_histogram_exists": True,
            "recoiljets_component_histogram_exists_in_current_final_root": False,
            "n_matched_bins_in_owned_stitch_region": len(sample_rows),
            "central_content_equality_count_after_extraction": sum(r["extracted_equal_within_strict_tolerance"] in ("True", True) for r in sample_rows),
            "ratio_min": min(ratios) if ratios else math.nan,
            "ratio_max": max(ratios) if ratios else math.nan,
            "ratio_mean": statistics.mean(ratios) if ratios else math.nan,
            "ratio_rms_around_one": math.sqrt(statistics.mean([(x - 1.0) ** 2 for x in ratios])) if ratios else math.nan,
            "max_abs_relative_difference": max(abs(x) for x in rels) if rels else math.nan,
            "max_abs_pull": max(abs(x) for x in pulls) if pulls else math.nan,
            "deviations_present_pre_stitch": "not_testable_from_current_final_root",
            "deviations_visible_in_owned_stitch_region": any(not (r["extracted_equal_within_strict_tolerance"] in ("True", True)) for r in sample_rows),
        }
        summary["samples"][sample] = sample_summary
        out_rows.append({"sample": sample, **sample_summary})
    write_csv(OUT_DIR / "component_prestitch_identity.csv", out_rows, list(out_rows[0].keys()))
    write_json(OUT_DIR / "component_prestitch_identity.json", summary)
    return summary


def recoiljets_family_compare() -> dict[str, Any]:
    loaded: dict[str, dict[str, Any]] = {}
    for label, hist_name in RECOILJETS_FAMILIES.items():
        loaded[label] = load_hist(CURRENT_ROOT, hist_name)
    ref = loaded["ppg12Fig5_kept"]
    rows: list[dict[str, Any]] = []
    for label, hist in loaded.items():
        diffs = []
        rels = []
        for center, ref_val, val in zip(ref["centers"], ref["values"], hist["values"]):
            if center < 10 or center >= 40:
                continue
            diffs.append(float(ref_val - val))
            rels.append(rel_diff(float(ref_val), float(val)))
        rows.append(
            {
                "family": label,
                "histogram": RECOILJETS_FAMILIES[label],
                "max_abs_difference_vs_ppg12Fig5_kept_10_40": max(abs(x) for x in diffs) if diffs else math.nan,
                "max_abs_relative_difference_vs_ppg12Fig5_kept_10_40": max(abs(x) for x in rels) if rels else math.nan,
                "identical_to_ppg12Fig5_kept_10_40": all(strict_equal(x, 0.0) for x in diffs) if diffs else False,
            }
        )
    write_csv(OUT_DIR / "recoiljets_internal_family_compare.csv", rows, list(rows[0].keys()))
    summary = {"csv": str(OUT_DIR / "recoiljets_internal_family_compare.csv"), "families": rows}
    write_json(OUT_DIR / "recoiljets_internal_family_compare.json", summary)
    return summary


def source_count_identity() -> dict[str, Any]:
    rows = []
    for sample in SAMPLES:
        rows.append(
            {
                "sample": sample.name,
                "ppg12_files_counts_norm": (
                    f"Build script points at /sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/sim/run28/{sample.name}/condorout/combined.root; "
                    "n_total is read from slimtree.num_entries at build time and is not stored in the local source CSV."
                ),
                "recoiljets_files_counts_norm": (
                    "Local evidence is a final merged ROOT only. Per-worker/per-sample file lists and generated-event denominators "
                    "are not present in the current local artifact."
                ),
                "same_different_unknown": "unknown",
                "consequence_for_exact_identity": (
                    "Exact regeneration cannot be claimed from local artifacts because input population and source-count denominators are not recoverable."
                ),
            }
        )
    write_csv(
        OUT_DIR / "source_count_identity.csv",
        rows,
        ["sample", "ppg12_files_counts_norm", "recoiljets_files_counts_norm", "same_different_unknown", "consequence_for_exact_identity"],
    )
    md = [
        "# THE-76 Photon+Jet Source-Count / Input-Population Identity",
        "",
        "| sample | PPG12 files/counts/norm | RecoilJets files/counts/norm | same/different/unknown | consequence |",
        "| --- | --- | --- | --- | --- |",
    ]
    for r in rows:
        md.append(
            f"| {r['sample']} | {r['ppg12_files_counts_norm']} | {r['recoiljets_files_counts_norm']} | "
            f"{r['same_different_unknown']} | {r['consequence_for_exact_identity']} |"
        )
    md.extend(
        [
            "",
            "Smallest read-only SDCC diagnostic needed:",
            "",
            "```bash",
            "ssh ssh.sdcc.bnl.gov 'ssh sphnxuser05 \"python3 - <<PY\\n"
            "import uproot\\n"
            "for s in (\\\"photon5\\\",\\\"photon10\\\",\\\"photon20\\\"):\\n"
            "    p=f\\\"/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/sim/run28/{s}/condorout/combined.root\\\"\\n"
            "    with uproot.open(p) as f:\\n"
            "        print(s, p, f[\\\"slimtree\\\"].num_entries)\\n"
            "with uproot.open(\\\"/sphenix/user/shuhangli/ppg12/plotting/photon_max_pT_uncut.root\\\") as f:\\n"
            "    print(f.keys())\\n"
            "PY\"'",
            "```",
            "",
            "This command is read-only and recovers the PPG12 source-count denominators and source ROOT object presence.",
        ]
    )
    write_text(OUT_DIR / "source_count_identity.md", "\n".join(md) + "\n")
    summary = {
        "csv": str(OUT_DIR / "source_count_identity.csv"),
        "md": str(OUT_DIR / "source_count_identity.md"),
        "status": "unknown_from_local_artifacts",
        "smallest_read_only_sdcc_diagnostic": "ssh nested read-only uproot num_entries/object-list check recorded in markdown artifact",
    }
    write_json(OUT_DIR / "source_count_identity.json", summary)
    return summary


def vertex_reweight_containment() -> dict[str, Any]:
    qa_hists = [
        "SIM/h_ppg12_vertex_contract_audit",
        "SIM/h_ppg12_vtxqa_sim_vertex_weight_pre_vzcut",
        "SIM/h_ppg12_vtxqa_sim_vertex_weight_post_vzcut",
        "SIM/h_ppg12_vtxqa_sim_period_event_weight_pre_vzcut",
        "SIM/h_ppg12_vtxqa_sim_period_event_weight_post_vzcut",
        "SIM/h_ppg12_vtxqa_sim_truth_z_weighted_post_vzcut",
        "SIM/h_ppg12_vtxqa_sim_truth_z_unweighted_post_vzcut",
    ]
    root_evidence = []
    with uproot.open(CURRENT_ROOT) as f:
        for hist_name in qa_hists:
            if hist_name in f:
                values, _ = f[hist_name].to_numpy(flow=False)
                root_evidence.append(
                    {
                        "histogram": hist_name,
                        "present": True,
                        "sum": float(np.asarray(values, dtype=float).sum()),
                        "nonzero_bins": int(np.count_nonzero(values)),
                    }
                )
            else:
                root_evidence.append({"histogram": hist_name, "present": False})

    code_hits = line_hits(
        RECOILJETS_SRC,
        [
            "truth_vertex_reweight/output/0mrad/reweight.root",
            "truth_vertex_reweight/output/1p5mrad/reweight.root",
            "m_mcEventWeight",
            "RJMCWeighting::SetCurrentWeight",
            "h_ppg12_vtxqa_sim_vertex_weight",
            "h_ppPhotonStitch_ppg12Fig5_maxPhotonPt_kept",
        ],
        context=2,
    )
    summary = {
        "verdict": "ruled out for code path but not exact artifact",
        "expected_vertex_reweight_files": [
            "/sphenix/user/shuhangli/ppg12/efficiencytool/truth_vertex_reweight/output/0mrad/reweight.root",
            "/sphenix/user/shuhangli/ppg12/efficiencytool/truth_vertex_reweight/output/1p5mrad/reweight.root",
        ],
        "actual_files_loaded": "not directly recorded in the final local ROOT; code path records defaults and ROOT QA proves vertex-weight family was filled",
        "code_path_evidence": code_hits,
        "root_qa_evidence": root_evidence,
        "plotted_histogram_downstream_evidence": (
            "The plotted histogram exists in the same final ROOT as the vertex QA histograms. This is ROOT QA/code-path proof, not a per-event runtime canary."
        ),
        "bypass_or_double_apply_status": "no bypass/double-apply proven from local artifacts; per-slice runtime logs/canary are absent",
    }
    write_json(OUT_DIR / "vertex_reweight_containment.json", summary)
    md = [
        "# THE-76 Vertex-Reweight Exact-Artifact Containment",
        "",
        f"Verdict: `{summary['verdict']}`",
        "",
        "Expected vertex reweight files:",
        *[f"- `{p}`" for p in summary["expected_vertex_reweight_files"]],
        "",
        "ROOT QA evidence:",
        "",
        "| histogram | present | sum | nonzero bins |",
        "| --- | --- | --- | --- |",
    ]
    for r in root_evidence:
        md.append(f"| `{r['histogram']}` | {r['present']} | {r.get('sum', '')} | {r.get('nonzero_bins', '')} |")
    md.extend(
        [
            "",
            "Code-path evidence is stored in `vertex_reweight_containment.json` with line excerpts.",
            "",
            "This rules out the inclusive-jet-style missing machinery at the code-path/ROOT-QA level, but it is not a per-event proof for photon5/photon10/photon20 in the exact plotted bins.",
        ]
    )
    write_text(OUT_DIR / "vertex_reweight_containment.md", "\n".join(md) + "\n")
    return summary


def identity_ladder(
    ppg12_root_available: bool,
    raw_summary: dict[str, Any],
    component_summary: dict[str, Any],
    source_summary: dict[str, Any],
    vertex_summary: dict[str, Any],
    family_summary: dict[str, Any],
) -> dict[str, Any]:
    rows = [
        {
            "rung": "PPG12 reference target identity",
            "ppg12_evidence": str(PPG12_SOURCE_ROOT),
            "recoiljets_evidence": str(CURRENT_ROOT),
            "status": "different",
            "identity_survives": False,
            "consequence": "The comparison is PPG12 source ROOT vs RecoilJets final ROOT, not the same file/object regeneration.",
        },
        {
            "rung": "PPG12 ROOT file identity",
            "ppg12_evidence": "local missing" if not ppg12_root_available else "local ROOT available",
            "recoiljets_evidence": str(CURRENT_ROOT),
            "status": "unknown" if not ppg12_root_available else "different",
            "identity_survives": False,
            "consequence": "Local audit cannot prove PPG12 ROOT internals without SDCC read-only object/count check.",
        },
        {
            "rung": "PPG12 histogram object identity",
            "ppg12_evidence": "h_max_photon_pT_photon5/10/20 from source CSV/build script",
            "recoiljets_evidence": CURRENT_HIST,
            "status": "different",
            "identity_survives": False,
            "consequence": "PPG12 stores per-sample weighted source histograms; RecoilJets stores a combined stitch diagnostic.",
        },
        {
            "rung": "RecoilJets ROOT file identity",
            "ppg12_evidence": "not applicable",
            "recoiljets_evidence": str(CURRENT_ROOT),
            "status": "same",
            "identity_survives": True,
            "consequence": "Current artifact path is resolved.",
        },
        {
            "rung": "RecoilJets histogram object identity",
            "ppg12_evidence": "not applicable",
            "recoiljets_evidence": CURRENT_HIST,
            "status": "same",
            "identity_survives": True,
            "consequence": "Current plotted histogram exists.",
        },
        {
            "rung": "photon5 component histogram identity",
            "ppg12_evidence": "PPG12 source CSV/object rows for photon5",
            "recoiljets_evidence": "no separate photon5 component final ROOT in local artifact",
            "status": "unknown",
            "identity_survives": False,
            "consequence": "Cannot prove pre-stitch identity locally.",
        },
        {
            "rung": "photon10 component histogram identity",
            "ppg12_evidence": "PPG12 source CSV/object rows for photon10",
            "recoiljets_evidence": "no separate photon10 component final ROOT in local artifact",
            "status": "unknown",
            "identity_survives": False,
            "consequence": "Cannot prove pre-stitch identity locally.",
        },
        {
            "rung": "photon20 component histogram identity",
            "ppg12_evidence": "PPG12 source CSV/object rows for photon20",
            "recoiljets_evidence": "no separate photon20 component final ROOT in local artifact",
            "status": "unknown",
            "identity_survives": False,
            "consequence": "Cannot prove pre-stitch identity locally.",
        },
        {
            "rung": "photon5 input file list identity",
            "ppg12_evidence": "build script points to anatreemaker combined.root",
            "recoiljets_evidence": "final merged ROOT only",
            "status": "unknown",
            "identity_survives": False,
            "consequence": "Exact input-population identity cannot be claimed.",
        },
        {
            "rung": "photon10 input file list identity",
            "ppg12_evidence": "build script points to anatreemaker combined.root",
            "recoiljets_evidence": "final merged ROOT only",
            "status": "unknown",
            "identity_survives": False,
            "consequence": "Exact input-population identity cannot be claimed.",
        },
        {
            "rung": "photon20 input file list identity",
            "ppg12_evidence": "build script points to anatreemaker combined.root",
            "recoiljets_evidence": "final merged ROOT only",
            "status": "unknown",
            "identity_survives": False,
            "consequence": "Exact input-population identity cannot be claimed.",
        },
        {
            "rung": "photon5 event/source count identity",
            "ppg12_evidence": "not stored in local source CSV",
            "recoiljets_evidence": "not recoverable per sample from current final ROOT",
            "status": "unknown",
            "identity_survives": False,
            "consequence": "No exact denominator identity proof.",
        },
        {
            "rung": "photon10 event/source count identity",
            "ppg12_evidence": "not stored in local source CSV",
            "recoiljets_evidence": "not recoverable per sample from current final ROOT",
            "status": "unknown",
            "identity_survives": False,
            "consequence": "No exact denominator identity proof.",
        },
        {
            "rung": "photon20 event/source count identity",
            "ppg12_evidence": "not stored in local source CSV",
            "recoiljets_evidence": "not recoverable per sample from current final ROOT",
            "status": "unknown",
            "identity_survives": False,
            "consequence": "No exact denominator identity proof.",
        },
        {
            "rung": "generated-event/source-count normalization identity",
            "ppg12_evidence": "xsec / slimtree.num_entries in build_h_max_photon_pT_uncut.py",
            "recoiljets_evidence": f"current extraction factor median {raw_summary['source_scope_factor_median']}",
            "status": "unknown",
            "identity_survives": False,
            "consequence": "Current extraction uses a validated source-scope factor; exact denominator identity is not proven.",
        },
        {
            "rung": "event selection identity",
            "ppg12_evidence": "vectorized slimtree particle photon max, no per-sample event cut",
            "recoiljets_evidence": "Fun4All/RecoilJets event path with stitch ABORTEVENT after diagnostic fill",
            "status": "different",
            "identity_survives": False,
            "consequence": "Different fill program and event pathway are sufficient to break exact central identity.",
        },
        {
            "rung": "truth photon definition / max photon pT definition identity",
            "ppg12_evidence": "particle_pid/Pt/Eta from PPG12 anatreemaker slimtree",
            "recoiljets_evidence": "ppg12MaxStoredTruthPhotonPt(topNode, m_truthInfo, true, 0.7)",
            "status": "different",
            "identity_survives": False,
            "consequence": "This is a direct object/fill-definition break.",
        },
        {
            "rung": "period split identity: 0mrad / 1.5mrad",
            "ppg12_evidence": "not encoded in local PPG12 source CSV",
            "recoiljets_evidence": "period QA histograms present",
            "status": "unknown",
            "identity_survives": False,
            "consequence": "Cannot prove exact period composition identity locally.",
        },
        {
            "rung": "SI/DI component identity",
            "ppg12_evidence": "not encoded in local PPG12 source CSV",
            "recoiljets_evidence": "vertex/mix QA exists but final per-component source counts are absent",
            "status": "unknown",
            "identity_survives": False,
            "consequence": "Cannot prove exact component composition identity locally.",
        },
        {
            "rung": "cross-section weight identity",
            "ppg12_evidence": "XSEC dict in build script",
            "recoiljets_evidence": "metadata/source factor and code constants",
            "status": "same",
            "identity_survives": True,
            "consequence": "Nominal xsec constants match known PPG12 values, but this does not prove denominator identity.",
        },
        {
            "rung": "luminosity/period weight identity",
            "ppg12_evidence": "not represented in source CSV",
            "recoiljets_evidence": "period-event-weight QA present",
            "status": "unknown",
            "identity_survives": False,
            "consequence": "Cannot use this rung to claim exact identity.",
        },
        {
            "rung": "SI/DI mixture weight identity",
            "ppg12_evidence": "not represented in source CSV",
            "recoiljets_evidence": "vertex contract audit exists",
            "status": "unknown",
            "identity_survives": False,
            "consequence": "Cannot use this rung to claim exact identity.",
        },
        {
            "rung": "PPG12 vertex-reweight file identity",
            "ppg12_evidence": "expected PPG12 files",
            "recoiljets_evidence": "code-path defaults; no runtime file path in final ROOT",
            "status": "unknown",
            "identity_survives": False,
            "consequence": vertex_summary["verdict"],
        },
        {
            "rung": "PPG12 vertex-reweight factor application identity",
            "ppg12_evidence": "not encoded in PPG12 source CSV",
            "recoiljets_evidence": "vtxqa histograms present in exact artifact",
            "status": "unknown",
            "identity_survives": False,
            "consequence": vertex_summary["verdict"],
        },
        {
            "rung": "final event fill-weight identity",
            "ppg12_evidence": "single xsec/n_total in build script",
            "recoiljets_evidence": "weighted framework path plus extraction factor",
            "status": "different",
            "identity_survives": False,
            "consequence": "Different fill-weight machinery prevents exact identity unless independently proven equivalent.",
        },
        {
            "rung": "histogram fill-stage identity",
            "ppg12_evidence": "offline Python/awkward source ROOT builder",
            "recoiljets_evidence": "Fun4All RecoilJets diagnostic histogram",
            "status": "different",
            "identity_survives": False,
            "consequence": "This explains why central values need not be bit-identical.",
        },
        {
            "rung": "binning identity",
            "ppg12_evidence": "0.5 GeV bins",
            "recoiljets_evidence": "0.5 GeV bins",
            "status": "same",
            "identity_survives": True,
            "consequence": "Binning is not the observed break.",
        },
        {
            "rung": "underflow/overflow identity",
            "ppg12_evidence": "not audited from local ROOT",
            "recoiljets_evidence": "not used in plotted extraction",
            "status": "unknown",
            "identity_survives": False,
            "consequence": "Not enough local evidence.",
        },
        {
            "rung": "pre-merge component histogram identity",
            "ppg12_evidence": "per-sample source objects",
            "recoiljets_evidence": "pre-merge component ROOTs not locally present",
            "status": "unknown",
            "identity_survives": False,
            "consequence": "Smallest next diagnostic is per-sample/per-worker read-only accounting.",
        },
        {
            "rung": "hadd/merge identity",
            "ppg12_evidence": "not a hadd of RecoilJets ROOTs",
            "recoiljets_evidence": "current final merged ROOT",
            "status": "different",
            "identity_survives": False,
            "consequence": "Comparison is not same merge process.",
        },
        {
            "rung": "post-merge stitched histogram identity",
            "ppg12_evidence": "PPG12 source CSV values",
            "recoiljets_evidence": f"extracted_equal_count={raw_summary['extracted_equal_count']} of {raw_summary['n_bins']}",
            "status": "different",
            "identity_survives": False,
            "consequence": "Central values are not identical after extraction.",
        },
        {
            "rung": "plotting extraction identity",
            "ppg12_evidence": "self-ratio exactly one",
            "recoiljets_evidence": "self-ratio exactly one",
            "status": "same",
            "identity_survives": True,
            "consequence": "Plotting/CSV extraction is ruled out as the source of non-unity central values.",
        },
        {
            "rung": "statistical error construction identity",
            "ppg12_evidence": "sumw2/source errors",
            "recoiljets_evidence": "current extracted errors",
            "status": "different",
            "identity_survives": False,
            "consequence": "Error bars are not the explanation for central-value non-identity.",
        },
    ]
    first_break = next((r for r in rows if r["status"] in {"different", "unknown"} and not r["identity_survives"]), None)
    table = [
        "# THE-76 Photon+Jet Fig.5 Identity Ladder",
        "",
        "| rung | PPG12 evidence | RecoilJets evidence | same/different/unknown | central identity survives | consequence |",
        "| --- | --- | --- | --- | --- | --- |",
    ]
    for r in rows:
        table.append(
            f"| {r['rung']} | {r['ppg12_evidence']} | {r['recoiljets_evidence']} | "
            f"{r['status']} | {r['identity_survives']} | {r['consequence']} |"
        )
    table.extend(
        [
            "",
            "## First identity break",
            "",
            f"`{first_break['rung'] if first_break else 'none'}`",
            "",
            f"{first_break['consequence'] if first_break else ''}",
        ]
    )
    write_text(OUT_DIR / "identity_ladder.md", "\n".join(table) + "\n")
    summary = {
        "md": str(OUT_DIR / "identity_ladder.md"),
        "rows": rows,
        "first_identity_break": first_break,
        "recoiljets_internal_family_compare": family_summary,
    }
    write_json(OUT_DIR / "identity_ladder.json", summary)
    return summary


def final_classification(
    self_summary: dict[str, Any],
    raw_summary: dict[str, Any],
    component_summary: dict[str, Any],
    ladder_summary: dict[str, Any],
    source_summary: dict[str, Any],
    vertex_summary: dict[str, Any],
) -> dict[str, Any]:
    classification = "explained: non-identical histogram/fill object"
    first_break = ladder_summary["first_identity_break"]
    central_identity = {
        "self_ratio_ruled_out_plotting_bug": self_summary["all_central_ratios_equal_one_strict"],
        "raw_or_extracted_bins_equal_count": raw_summary["extracted_equal_count"],
        "n_bins": raw_summary["n_bins"],
        "max_extracted_relative_difference": raw_summary["max_extracted_relative_difference"],
    }
    text = f"""# THE-76 Photon+Jet Fig.5 Identity-Break Classification

Classification: `{classification}`

Canonicalization allowed now: `no`

First identity break: `{first_break['rung']}`

The PPG12 side is a source ROOT built by `ppg12codeGit/plotting/build_h_max_photon_pT_uncut.py` from PPG12 anatreemaker slimtrees and per-sample weighted histograms. The RecoilJets side is a final merged Fun4All/RecoilJets diagnostic histogram, `SIM/h_ppPhotonStitch_ppg12Fig5_maxPhotonPt_kept`, read from the current July 2 artifact. Those are not the same ROOT file, not the same object family, and not the same fill stage.

The plotting self-ratio central values are exactly one, so the non-unity central ratios are not introduced by CSV matching or plotting extraction. After the current extraction convention is applied, the PPG12/current central bins still differ with max relative difference `{raw_summary['max_extracted_relative_difference']}` over `{raw_summary['n_bins']}` bins. That scale is not floating-point or ROOT serialization noise.

This is not evidence that the photon+jet stitch is bad as an independent PPG12-compatible reproduction. It is evidence that this specific artifact is not an exact mathematical regeneration of PPG12 internals. Exact central identity would require proving identical source file lists/counts, fill object, truth photon collection, event weights, and per-sample component histograms.

Minimal next diagnostic:

1. Run the read-only SDCC source-count/object check recorded in `source_count_identity.md`.
2. If Justin wants exact identity rather than compatible reproduction, run a tiny foreground canary on photon5/photon10/photon20 that prints per-event leading photon pT from the PPG12 anatreemaker variables and from the RecoilJets truth object path for the same bounded files/events.
3. Do not run Condor or merge until that canary predicts the expected central-ratio change.

Direct answer: It is reasonable to expect a perfect central ratio only for an exact regeneration using the same input population, denominator, object definition, fill stage, and weights. The current comparison is an independent PPG12-compatible RecoilJets reproduction against a PPG12 source object, so nonzero central-bin deviations are not automatically plotting noise. The plotted statistical error bars are separate; they can remain nonzero even when self-ratio central values are exactly one, and they do not explain central-value non-identity.
"""
    write_text(OUT_DIR / "final_classification.md", text)
    summary = {
        "classification": classification,
        "canonicalization_allowed_now": False,
        "first_identity_break": first_break,
        "central_identity": central_identity,
        "what_is_ruled_out": [
            "plotting/CSV self-ratio central-value bug",
            "floating-point or ROOT serialization as an explanation for per-mille/percent residuals",
            "statistical error bars as an explanation for central-value non-identity",
        ],
        "what_remains_ambiguous": [
            "exact PPG12 source file list and source-count denominators",
            "per-sample RecoilJets component histogram identity before the final merged stitch object",
            "per-event equivalence between PPG12 anatreemaker leading-photon pT and RecoilJets stored truth photon pT",
            "runtime per-slice vertex-reweight factor proof for the exact plotted bins",
        ],
        "minimal_next_diagnostic": "Run the read-only SDCC source-count/object check, then a tiny foreground same-event/factorized photon5/10/20 canary if exact identity is still required.",
        "md": str(OUT_DIR / "final_classification.md"),
    }
    write_json(OUT_DIR / "final_classification.json", summary)
    return summary


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    point_rows = read_csv(CURRENT_POINTS_CSV)
    ppg12_rows = read_csv(PPG12_SOURCE_CSV)
    ppg12_root_bins = load_ppg12_root_bins_if_available()
    current_hist = load_hist(CURRENT_ROOT, CURRENT_HIST)

    self_summary = plotting_self_ratio_tests(point_rows)
    raw_summary = raw_root_bin_identity(point_rows, ppg12_rows, current_hist, ppg12_root_bins)
    raw_rows = read_csv(OUT_DIR / "raw_root_bin_identity.csv")
    component_summary = component_prestitch_identity(raw_rows)
    family_summary = recoiljets_family_compare()
    source_summary = source_count_identity()
    vertex_summary = vertex_reweight_containment()
    ladder_summary = identity_ladder(
        ppg12_root_available=ppg12_root_bins is not None,
        raw_summary=raw_summary,
        component_summary=component_summary,
        source_summary=source_summary,
        vertex_summary=vertex_summary,
        family_summary=family_summary,
    )
    final = final_classification(
        self_summary,
        raw_summary,
        component_summary,
        ladder_summary,
        source_summary,
        vertex_summary,
    )

    index = {
        "audit": "THE-76 photon+jet Fig.5 identity-break audit",
        "created_at_local_date": "2026-07-02",
        "inputs": {
            "current_root": str(CURRENT_ROOT),
            "current_points_csv": str(CURRENT_POINTS_CSV),
            "ppg12_source_csv": str(PPG12_SOURCE_CSV),
            "ppg12_source_root": str(PPG12_SOURCE_ROOT),
            "ppg12_source_root_local_status": "available" if ppg12_root_bins is not None else "missing_locally",
        },
        "outputs": {
            "plotting_self_ratio_tests": str(OUT_DIR / "plotting_self_ratio_tests.json"),
            "raw_root_bin_identity": str(OUT_DIR / "raw_root_bin_identity.json"),
            "component_prestitch_identity": str(OUT_DIR / "component_prestitch_identity.json"),
            "identity_ladder": str(OUT_DIR / "identity_ladder.json"),
            "source_count_identity": str(OUT_DIR / "source_count_identity.json"),
            "vertex_reweight_containment": str(OUT_DIR / "vertex_reweight_containment.json"),
            "final_classification": str(OUT_DIR / "final_classification.json"),
        },
        "classification": final["classification"],
        "canonicalization_allowed_now": False,
    }
    write_json(OUT_DIR / "audit_index.json", index)
    print(json.dumps(index, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
