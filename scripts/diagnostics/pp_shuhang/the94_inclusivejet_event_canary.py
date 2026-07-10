#!/usr/bin/env python3
"""THE-94 inclusive-jet event-level canary.

This is intentionally read-only. It inspects existing local ROOT products for
event-level provenance needed to compare the RecoilJets strict Fig.6 fill
against the PPG12 RecoEffCalculator formula. If the local products only contain
merged histograms, it emits an explicit blocker and the minimal histogram-level
summary that can still be computed.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import uproot


REPO = Path(__file__).resolve().parents[3]
STRICT_ROOT_DIR = (
    REPO
    / "dataOutput/ppg12Parity/the76_ppg12_fig6_inclusive_strict_20260702_014757"
    / "final_roots/inclusivejet"
)
FORENSIC_DIR = (
    REPO
    / "dataOutput/ppg12Parity/THE-94_ppg12_pp_inclusivejet_sim_parity"
    / "local_forensics_20260702"
)
OUT_DIR = (
    REPO
    / "dataOutput/ppg12Parity/THE-94_ppg12_pp_inclusivejet_sim_parity"
    / "event_canary_20260702"
)

JET50_XSEC = 7.3113
SAMPLES: dict[str, dict[str, float]] = {
    "jet8": {"lo": 9.0, "hi": 14.0, "xsec": 1.15e7},
    "jet12": {"lo": 14.0, "hi": 21.0, "xsec": 1.4903e6},
    "jet20": {"lo": 21.0, "hi": 32.0, "xsec": 6.2623e4},
    "jet30": {"lo": 32.0, "hi": 42.0, "xsec": 2.5298e3},
    "jet40": {"lo": 42.0, "hi": 100.0, "xsec": 1.3553e2},
}

STRICT_KEPT = "SIM/h_ppInclusiveJetStitch_ppg12Fig6EventWeighted_r04_maxTruthJetPt_kept"
STRICT_ALL = "SIM/h_ppInclusiveJetStitch_ppg12Fig6EventWeighted_r04_maxTruthJetPt_all"
RAW_KEPT = "SIM/h_ppInclusiveJetStitch_ppg12TruthSpectrum_r04_maxTruthJetPt_kept"
META = "SIM/h_ppInclusiveJetStitch_ppg12Fig6EventWeighted_metadata"
EVENT_WEIGHT_HIST = "SIM/h_ppg12_vtxqa_sim_period_event_weight_pre_vzcut"
VERTEX_WEIGHT_HIST = "SIM/h_ppg12_vtxqa_sim_vertex_weight_pre_vzcut"
PERIOD_VALUES = "SIM/h_ppg12_period_contract_values"


def read_reference_integrals() -> tuple[dict[str, float], dict[str, float]]:
    path = FORENSIC_DIR / "sdcc_ppg12_reference_integrals.csv"
    ian: dict[str, float] = {}
    period: dict[str, float] = {}
    if not path.exists():
        return ian, period
    with path.open() as f:
        for row in csv.DictReader(f):
            sample = row["sample"]
            if row["reference"] == "ian_no_suffix":
                ian[sample] = float(row["integral"])
            elif row["reference"] == "period_combined_all":
                period[sample] = float(row["integral"])
    return ian, period


def find_root(sample: str) -> Path | None:
    matches = sorted(STRICT_ROOT_DIR.glob(f"RecoilJets_{sample}_ALL_*.root"))
    return matches[0] if matches else None


def classnames(root_path: Path) -> dict[str, str]:
    with uproot.open(root_path) as f:
        return dict(f.classnames(recursive=True))


def has_event_tree(root_path: Path) -> bool:
    return any(("TTree" in cls or "TNtuple" in cls) for cls in classnames(root_path).values())


def hist_values(root_path: Path, key: str) -> tuple[np.ndarray, np.ndarray] | None:
    with uproot.open(root_path) as f:
        if key not in f:
            return None
        obj = f[key]
        values, edges = obj.to_numpy(flow=False)
        return np.asarray(values, dtype=float), np.asarray(edges, dtype=float)


def integral_window(root_path: Path, key: str, lo: float, hi: float) -> tuple[float | None, int | None]:
    hv = hist_values(root_path, key)
    if hv is None:
        return None, None
    values, edges = hv
    centers = 0.5 * (edges[:-1] + edges[1:])
    mask = (centers >= lo) & (centers < hi)
    return float(values[mask].sum()), int(np.count_nonzero(values[mask]))


def hist_mean(root_path: Path, key: str) -> tuple[float | None, float | None]:
    hv = hist_values(root_path, key)
    if hv is None:
        return None, None
    values, edges = hv
    total = float(values.sum())
    if total <= 0.0:
        return 0.0, total
    centers = 0.5 * (edges[:-1] + edges[1:])
    return float((values * centers).sum() / total), total


def metadata_values(root_path: Path) -> dict[str, float | None]:
    hv = hist_values(root_path, META)
    if hv is None:
        return {}
    values, _ = hv
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
    return {label: float(values[i]) if i < len(values) else None for i, label in enumerate(labels)}


def period_values(root_path: Path) -> dict[str, float | None]:
    hv = hist_values(root_path, PERIOD_VALUES)
    if hv is None:
        return {}
    values, _ = hv
    labels = [
        "run_min",
        "run_max_excl",
        "lumi_pb_inv",
        "lumi_target_pb_inv",
        "lumi_weight",
        "f_single",
        "f_double",
        "z_closure_cm",
        "component_double",
        "mix_weight",
        "mix_auto",
        "vertex_file_auto",
        "data_period_filter",
    ]
    return {label: float(values[i]) if i < len(values) else None for i, label in enumerate(labels)}


def ratio(num: float | None, den: float | None) -> float | None:
    if num is None or den is None or den == 0.0 or not math.isfinite(den):
        return None
    return num / den


def fmt(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        if not math.isfinite(value):
            return ""
        return f"{value:.10g}"
    return str(value)


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: fmt(row.get(field)) for field in fields})


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    ian_ref, period_ref = read_reference_integrals()
    event_rows: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []
    inventory: dict[str, Any] = {}
    hypothesis: list[dict[str, str]] = []

    for sample, cfg in SAMPLES.items():
        root_path = find_root(sample)
        if root_path is None:
            event_rows.append({"sample": sample, "status": "blocked_missing_root"})
            summary_rows.append({"sample": sample, "status": "blocked_missing_root"})
            continue

        classes = classnames(root_path)
        tree_keys = [key for key, cls in classes.items() if "TTree" in cls or "TNtuple" in cls]
        inventory[sample] = {
            "root": str(root_path),
            "n_keys": len(classes),
            "tree_keys": tree_keys,
            "has_event_tree": bool(tree_keys),
        }

        xsec = cfg["xsec"]
        xsec_over_jet50 = xsec / JET50_XSEC
        lo = cfg["lo"]
        hi = cfg["hi"]
        strict_kept, kept_nonzero = integral_window(root_path, STRICT_KEPT, lo, hi)
        strict_all, _ = integral_window(root_path, STRICT_ALL, lo, hi)
        raw_kept, raw_nonzero = integral_window(root_path, RAW_KEPT, lo, hi)
        event_weight_mean, event_weight_entries = hist_mean(root_path, EVENT_WEIGHT_HIST)
        vertex_weight_mean, vertex_weight_entries = hist_mean(root_path, VERTEX_WEIGHT_HIST)
        meta = metadata_values(root_path)
        period = period_values(root_path)
        r_ian = ratio(strict_kept, ian_ref.get(sample))
        r_period = ratio(strict_kept, period_ref.get(sample))

        event_rows.append(
            {
                "sample": sample,
                "input_family_path": str(root_path),
                "component": "ALL_hadded_unknown",
                "period": "ALL_hadded_unknown",
                "event_index": "",
                "max_r04_truth_jet_pt": "",
                "owned_window_pass": "",
                "histogram_family_filled": "strict_july2_recoiljets_fig6",
                "fill_target": "_kept",
                "raw_fill_count": "",
                "xsec": xsec,
                "jet50_xsec": JET50_XSEC,
                "xsec_over_jet50": xsec_over_jet50,
                "period_lumi_factor": "",
                "si_di_mix_factor": "",
                "hard_truth_vertex": "",
                "mb_truth_vertex_if_di": "",
                "hard_truth_vertex_weight": "",
                "mb_truth_vertex_weight_if_di": "",
                "final_event_weight": "",
                "final_histogram_fill_weight": "",
                "equivalent_ppg12_formula_weight": "",
                "ratio_our_fill_over_ppg12_formula": "",
                "status": "blocked_missing_event_tree",
                "blocker": "local final ROOT contains histograms only; no TTree/TNtuple/event audit rows",
            }
        )

        summary_rows.append(
            {
                "sample": sample,
                "component": "ALL_hadded_unknown",
                "period": "ALL_hadded_unknown",
                "root": str(root_path),
                "has_event_tree": bool(tree_keys),
                "events_seen_metadata": meta.get("events_seen_raw"),
                "events_in_owned_pt_window": "",
                "raw_integral": raw_kept,
                "weighted_integral": strict_kept,
                "strict_all_integral_owned_window": strict_all,
                "expected_ppg12_weighted_integral_if_computable": "",
                "ian_no_suffix_integral": ian_ref.get(sample),
                "period_combined_integral": period_ref.get(sample),
                "ratio_to_ian_no_suffix_reference": r_ian,
                "ratio_to_true_period_combined_reference": r_period,
                "nonzero_bins_kept": kept_nonzero,
                "nonzero_bins_raw": raw_nonzero,
                "event_weight_mean_from_hist": event_weight_mean,
                "event_weight_hist_entries": event_weight_entries,
                "vertex_weight_mean_from_hist": vertex_weight_mean,
                "vertex_weight_hist_entries": vertex_weight_entries,
                "metadata_xsec_over_jet50_sum": meta.get("xsec_over_jet50"),
                "metadata_includes_current_weight_sum": meta.get("includes_current_weight"),
                "period_lumi_weight_sum": period.get("lumi_weight"),
                "period_mix_weight_sum": period.get("mix_weight"),
                "movement_vs_old_flat_1p94": "away" if r_ian is not None and abs(r_ian - 1.94) > 0.15 else "near_old_flat",
                "status": "hist_summary_only",
            }
        )

    all_have_trees = all(item.get("has_event_tree") for item in inventory.values())
    if not all_have_trees:
        hypothesis.extend(
            [
                {
                    "hypothesis": "H1_wrong_configured_input_component_population",
                    "status": "open_moved_up",
                    "evidence": "final roots are ALL-hadded histogram products with no event/component rows",
                },
                {
                    "hypothesis": "H2_wrong_SI_DI_period_composition",
                    "status": "open_moved_up",
                    "evidence": "period/component information is hadded into scalar histograms and cannot be factorized event-by-event",
                },
                {
                    "hypothesis": "H3_wrong_strict_run_event_weight_factor",
                    "status": "blocked",
                    "evidence": "event-weight distribution histograms exist, but no per-event truth-z/mix/lumi/final-weight rows exist",
                },
                {
                    "hypothesis": "H4_DST_population_mismatch_vs_PPG12_period_combined",
                    "status": "open_moved_up",
                    "evidence": "strict/true-period-combined ratios remain sample-dependent in histogram summaries",
                },
                {
                    "hypothesis": "H5_merge_hadd_composition_issue",
                    "status": "open_moved_up",
                    "evidence": "only ALL-hadded local final roots are available, so component-level hadd composition cannot be audited from event rows",
                },
                {
                    "hypothesis": "H6_jet8_legacy_exception",
                    "status": "down_non_nominal",
                    "evidence": "Blair legacy source is documented as an upstream historical culprit, not nominal policy",
                },
            ]
        )

    event_fields = [
        "sample",
        "input_family_path",
        "component",
        "period",
        "event_index",
        "max_r04_truth_jet_pt",
        "owned_window_pass",
        "histogram_family_filled",
        "fill_target",
        "raw_fill_count",
        "xsec",
        "jet50_xsec",
        "xsec_over_jet50",
        "period_lumi_factor",
        "si_di_mix_factor",
        "hard_truth_vertex",
        "mb_truth_vertex_if_di",
        "hard_truth_vertex_weight",
        "mb_truth_vertex_weight_if_di",
        "final_event_weight",
        "final_histogram_fill_weight",
        "equivalent_ppg12_formula_weight",
        "ratio_our_fill_over_ppg12_formula",
        "status",
        "blocker",
    ]
    summary_fields = [
        "sample",
        "component",
        "period",
        "root",
        "has_event_tree",
        "events_seen_metadata",
        "events_in_owned_pt_window",
        "raw_integral",
        "weighted_integral",
        "strict_all_integral_owned_window",
        "expected_ppg12_weighted_integral_if_computable",
        "ian_no_suffix_integral",
        "period_combined_integral",
        "ratio_to_ian_no_suffix_reference",
        "ratio_to_true_period_combined_reference",
        "nonzero_bins_kept",
        "nonzero_bins_raw",
        "event_weight_mean_from_hist",
        "event_weight_hist_entries",
        "vertex_weight_mean_from_hist",
        "vertex_weight_hist_entries",
        "metadata_xsec_over_jet50_sum",
        "metadata_includes_current_weight_sum",
        "period_lumi_weight_sum",
        "period_mix_weight_sum",
        "movement_vs_old_flat_1p94",
        "status",
    ]
    hyp_fields = ["hypothesis", "status", "evidence"]

    write_csv(OUT_DIR / "event_factor_rows.csv", event_rows, event_fields)
    write_csv(OUT_DIR / "component_summary.csv", summary_rows, summary_fields)
    write_csv(OUT_DIR / "hypothesis_scoreboard.csv", hypothesis, hyp_fields)

    manifest = {
        "task": "THE-94 ppg12_pp_inclusivejet_sim_parity",
        "mode": "local_read_only_event_canary",
        "strict_root_dir": str(STRICT_ROOT_DIR),
        "event_level_rows_available": all_have_trees,
        "root_inventory": inventory,
        "outputs": {
            "event_factor_rows": str(OUT_DIR / "event_factor_rows.csv"),
            "component_summary": str(OUT_DIR / "component_summary.csv"),
            "hypothesis_scoreboard": str(OUT_DIR / "hypothesis_scoreboard.csv"),
            "manifest": str(OUT_DIR / "event_canary_manifest.json"),
            "blocker_note": str(OUT_DIR / "event_canary_blocker.md"),
        },
        "blocked": not all_have_trees,
        "blocker": (
            "Existing local strict July 2 final ROOTs contain histograms but no "
            "TTree/TNtuple/event audit rows, so per-event PPG12 factorization "
            "cannot be reconstructed offline."
            if not all_have_trees
            else ""
        ),
        "minimal_next_instrumentation": (
            "Run a foreground RecoilJets canary over 1k-2k events per "
            "sample/component with a gated diagnostic TTree/CSV containing "
            "sample, input path, SI/DI, period, event index, max R=0.4 truth "
            "jet pT, window pass, xsec/jet50, lumi weight, mix weight, truth "
            "vertices, truth-vertex weights, final event weight, fill weight, "
            "and PPG12 formula ratio."
            if not all_have_trees
            else ""
        ),
    }
    (OUT_DIR / "event_canary_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")

    blocker_lines = [
        "# THE-94 Event Canary Blocker",
        "",
        "The requested event-level factorization cannot be reconstructed from the local strict July 2 final ROOTs.",
        "",
        "Reason: the ROOTs contain merged histograms only; no TTree/TNtuple/event-audit object is present.",
        "",
        "Minimal temporary instrumentation needed:",
        "",
        "- Add a gated foreground-only diagnostic output, not a promoted production object.",
        "- Emit one row per event before the Fig.6 fill for 1k-2k events per sample/component.",
        "- Include sample, input path, SI/DI, period, event index, max R=0.4 truth jet pT, owned-window pass, xsec/jet50, lumi weight, mix weight, truth z values, truth-vertex weights, final event weight, final fill weight, PPG12-equivalent formula weight, and ratio.",
        "- Run locally/foreground on the smallest accessible inputs before any Condor rerun.",
        "",
        "Histogram-level summaries were still written to `component_summary.csv`.",
    ]
    (OUT_DIR / "event_canary_blocker.md").write_text("\n".join(blocker_lines) + "\n")

    print(f"Wrote {OUT_DIR / 'event_factor_rows.csv'}")
    print(f"Wrote {OUT_DIR / 'component_summary.csv'}")
    print(f"Wrote {OUT_DIR / 'hypothesis_scoreboard.csv'}")
    print(f"Wrote {OUT_DIR / 'event_canary_manifest.json'}")
    if not all_have_trees:
        print("BLOCKED: no event-level TTree/TNtuple/event audit rows in local final ROOTs")
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
