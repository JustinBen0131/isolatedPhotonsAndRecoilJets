#!/usr/bin/env python3
"""THE-76 Fig24/Fig3 Eiso path canary from existing artifacts.

This is a deliberately narrow diagnostic wrapper.  The current merged
RecoilJets ROOTs do not persist the per-cluster isolation path choice
(`stored ppg12_topo_raw_eiso_04` vs recomputed topo cone).  The script emits
the canary-shaped CSV/JSON artifacts requested for the Fig24/Fig3 high-Eiso
tail investigation and records the exact blocker instead of silently treating
the final mixed histogram as path-resolved evidence.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[3]
DEFAULT_PRIOR = REPO / "dataOutput/ppg12Parity/control_plane/audits/the76_fig24_fig3_eiso_tail_20260706"
DEFAULT_OUT = REPO / "dataOutput/ppg12Parity/control_plane/audits/the76_eiso_path_canary_20260706"
DEFAULT_ROOT = (
    REPO
    / "dataOutput/ppg12Parity/the76_ppg12_fig24_photonjet_fix_20260702_014217"
    / "final_roots/photonjet/RecoilJets_photonjet5plus10plus20_MERGED.root"
)


PER_CLUSTER_FIELDS = [
    "sample_component_source_file",
    "event_id",
    "cluster_id",
    "cluster_pt",
    "cluster_eta",
    "truth_photon_pt",
    "truth_match_status",
    "ppg12_cluster_iso_topo_04",
    "recoiljets_stored_ppg12_topo_raw_eiso_04",
    "stored_value_valid",
    "vertex_compatibility_passes",
    "recoiljets_recomputed_topo_cone_raw_isolation",
    "final_raw_isolation_chosen_by_recoiljets",
    "path_choice",
    "corrected_eiso_1p2_raw_plus_0p1",
    "fig24_pt_bin",
    "fig24_eiso_bin",
    "event_cluster_weight",
    "mismatch_label",
]


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def write_md(path: Path, lines: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines).rstrip() + "\n")


def f(row: dict[str, str], key: str) -> float | None:
    try:
        return float(row[key])
    except Exception:
        return None


def high_pt_final_rows(prior: Path) -> list[dict[str, Any]]:
    rows = read_csv(prior / "fig24_eiso_shape_comparison_by_pt.csv")
    out: list[dict[str, Any]] = []
    for row in rows:
        pt_low = f(row, "pt_low")
        pt_high = f(row, "pt_high")
        if pt_low is None or pt_high is None or pt_low < 30.0:
            continue
        out.append(
            {
                "path": "final_mixed_current",
                "pt_low": pt_low,
                "pt_high": pt_high,
                "entries": f(row, "current_entries"),
                "mean": f(row, "current_mean"),
                "median": None,
                "rms": f(row, "current_rms"),
                "q70": f(row, "current_cut70"),
                "q80": f(row, "current_cut80"),
                "q90": f(row, "current_cut90"),
                "tail_fraction_eiso_gt_2": (
                    None
                    if f(row, "current_cdf_le_2") is None
                    else 1.0 - float(row["current_cdf_le_2"])
                ),
                "tail_fraction_eiso_gt_3": (
                    None
                    if f(row, "current_cdf_le_3") is None
                    else 1.0 - float(row["current_cdf_le_3"])
                ),
                "tail_fraction_eiso_gt_4": (
                    None
                    if f(row, "current_cdf_le_4") is None
                    else 1.0 - float(row["current_cdf_le_4"])
                ),
                "ppg12_q70": f(row, "ppg12_cut70"),
                "ppg12_q80": f(row, "ppg12_cut80"),
                "ppg12_q90": f(row, "ppg12_cut90"),
                "sdcc_over_current_q70": f(row, "sdcc_over_current_cut70"),
                "sdcc_over_current_q80": f(row, "sdcc_over_current_cut80"),
                "sdcc_over_current_q90": f(row, "sdcc_over_current_cut90"),
                "status": "available_from_final_mixed_histogram",
            }
        )
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--prior", type=Path, default=DEFAULT_PRIOR)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--current-root", type=Path, default=DEFAULT_ROOT)
    args = ap.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)

    blocker = {
        "sample_component_source_file": "photon20/high-pT/current THE-76 merged and inspected worker artifacts",
        "event_id": "",
        "cluster_id": "",
        "cluster_pt": "",
        "cluster_eta": "",
        "truth_photon_pt": "",
        "truth_match_status": "",
        "ppg12_cluster_iso_topo_04": "",
        "recoiljets_stored_ppg12_topo_raw_eiso_04": "",
        "stored_value_valid": "",
        "vertex_compatibility_passes": "",
        "recoiljets_recomputed_topo_cone_raw_isolation": "",
        "final_raw_isolation_chosen_by_recoiljets": "",
        "path_choice": "blocked_missing_path_instrumentation",
        "corrected_eiso_1p2_raw_plus_0p1": "",
        "fig24_pt_bin": "30-36 GeV requested",
        "fig24_eiso_bin": "",
        "event_cluster_weight": "",
        "mismatch_label": (
            "current artifacts do not preserve per-cluster stored/recomputed raw-isolation "
            "path choice; exact same-cluster PPG12 comparison is blocked"
        ),
    }
    per_cluster_path = args.out / "per_cluster_eiso_path_canary.csv"
    write_csv(per_cluster_path, [blocker], PER_CLUSTER_FIELDS)

    usage_rows = [
        {
            "path": "stored",
            "clusters": "",
            "fraction": "",
            "pt_scope": "30-36 GeV photon20 signal clusters",
            "status": "blocked_missing_path_instrumentation",
            "evidence": "No AuAuPhotonIDTrainingTree and no stored/recomputed split histogram in current inspected artifacts.",
        },
        {
            "path": "recomputed",
            "clusters": "",
            "fraction": "",
            "pt_scope": "30-36 GeV photon20 signal clusters",
            "status": "blocked_missing_path_instrumentation",
            "evidence": "ppg12PhotonYieldRawEiso returns one raw value; current outputs do not persist the fallback decision.",
        },
        {
            "path": "final_mixed",
            "clusters": "",
            "fraction": "1.0 of persisted Fig24 shape",
            "pt_scope": "30-36 GeV high-pT summary",
            "status": "available",
            "evidence": str(args.prior / "fig24_eiso_shape_comparison_by_pt.csv"),
        },
    ]
    usage_path = args.out / "stored_vs_recomputed_usage_fraction.csv"
    write_csv(usage_path, usage_rows, ["path", "clusters", "fraction", "pt_scope", "status", "evidence"])

    quantile_rows = [
        {
            "path": "stored",
            "pt_low": 30.0,
            "pt_high": 36.0,
            "entries": "",
            "mean": "",
            "median": "",
            "rms": "",
            "q70": "",
            "q80": "",
            "q90": "",
            "tail_fraction_eiso_gt_2": "",
            "tail_fraction_eiso_gt_3": "",
            "tail_fraction_eiso_gt_4": "",
            "status": "blocked_missing_path_instrumentation",
        },
        {
            "path": "recomputed",
            "pt_low": 30.0,
            "pt_high": 36.0,
            "entries": "",
            "mean": "",
            "median": "",
            "rms": "",
            "q70": "",
            "q80": "",
            "q90": "",
            "tail_fraction_eiso_gt_2": "",
            "tail_fraction_eiso_gt_3": "",
            "tail_fraction_eiso_gt_4": "",
            "status": "blocked_missing_path_instrumentation",
        },
        *high_pt_final_rows(args.prior),
    ]
    quantile_fields = [
        "path",
        "pt_low",
        "pt_high",
        "entries",
        "mean",
        "median",
        "rms",
        "q70",
        "q80",
        "q90",
        "tail_fraction_eiso_gt_2",
        "tail_fraction_eiso_gt_3",
        "tail_fraction_eiso_gt_4",
        "ppg12_q70",
        "ppg12_q80",
        "ppg12_q90",
        "sdcc_over_current_q70",
        "sdcc_over_current_q80",
        "sdcc_over_current_q90",
        "status",
    ]
    quantile_path = args.out / "stored_vs_recomputed_eiso_quantiles.csv"
    write_csv(quantile_path, quantile_rows, quantile_fields)

    summary = {
        "current_state": (
            "Fig24/Fig3 high-Eiso tail is real in final mixed histograms, but the current ROOTs "
            "do not preserve the per-cluster raw-isolation path needed to decide stored vs recomputed."
        ),
        "diagnostic_scope": "photon20 high-pT 30-36 GeV signal-cluster Eiso path canary from existing artifacts",
        "current_root": str(args.current_root),
        "per_cluster_csv": str(per_cluster_path),
        "usage_fraction_csv": str(usage_path),
        "quantiles_csv": str(quantile_path),
        "stored_vs_recomputed_usage_fraction": {
            "stored": "blocked_missing_path_instrumentation",
            "recomputed": "blocked_missing_path_instrumentation",
            "final_mixed": "available",
        },
        "high_tail_driven_by_fallback_recompute": "not proven from current artifacts",
        "same_cluster_ppg12_comparison": {
            "status": "blocked",
            "blocker": (
                "PPG12 slimTree has cluster_iso_topo_04, but the current RecoilJets outputs do not "
                "preserve per-cluster stored raw Eiso, recomputed raw Eiso, or path choice."
            ),
        },
        "first_proven_mechanism": "not proven; first proven issue is missing path instrumentation",
        "exact_fix_rerun_needed": [
            "Add a diagnostic-only foreground canary around RecoilJets::ppg12PhotonYieldRawEiso.",
            "For Fig24-matched photon20 clusters, write sample/component/source file, event id, cluster id, cluster pT/eta, truth match, stored valid flag, vertex compatibility, stored raw Eiso, recomputed raw Eiso, final raw Eiso, path choice, corrected Eiso, pT/Eiso bins, and weight.",
            "Run it on a bounded photon20 high-pT foreground subset before any broad production.",
            "If stored-only matches PPG12 and recomputed is tail-heavy, fix by forcing/repairing the stored PPG12 topo-isolation branch for this parity path.",
            "If stored-only also disagrees, repair PhotonClusterBuilder stored branch construction.",
            "If raw Eiso agrees but cluster set differs, align the signal-cluster truth-match selection.",
        ],
        "canonicalization_allowed_now": False,
    }
    summary_path = args.out / "eiso_path_canary_summary.json"
    write_json(summary_path, summary)
    write_md(
        args.out / "eiso_path_canary_summary.md",
        [
            "# THE-76 Eiso Path Canary From Existing Artifacts",
            "",
            f"- Current state: {summary['current_state']}",
            f"- Scope: {summary['diagnostic_scope']}",
            f"- Per-cluster CSV: `{per_cluster_path}`",
            f"- Stored/recomputed usage: `{usage_path}`",
            f"- Quantiles: `{quantile_path}`",
            "- Stored-vs-recomputed usage fraction: blocked for stored and recomputed paths; final mixed shape is available.",
            "- High-Eiso tail driven by fallback/recompute: not proven from current artifacts.",
            "- Same-cluster PPG12 comparison: blocked by missing RecoilJets per-cluster path fields.",
            "- First proven mechanism: not a physics mechanism; the proven blocker is missing path instrumentation.",
            "- Canonicalization allowed now: no.",
            "",
            "## Exact Fix / Rerun Needed",
            *[f"- {x}" for x in summary["exact_fix_rerun_needed"]],
        ],
    )
    write_json(
        args.out / "audit_index.json",
        {
            "audit": "the76_eiso_path_canary_20260706",
            "outputs": {
                "per_cluster_csv": str(per_cluster_path),
                "usage_fraction_csv": str(usage_path),
                "quantiles_csv": str(quantile_path),
                "summary_json": str(summary_path),
                "summary_md": str(args.out / "eiso_path_canary_summary.md"),
            },
        },
    )
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
