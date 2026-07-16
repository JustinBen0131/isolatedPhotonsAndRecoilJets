#!/usr/bin/env python3
"""Close the THE-76 Fig24/Fig3 Eiso contract audit from existing artifacts.

This is intentionally read-only with respect to analysis outputs.  It
consolidates the prior Fig24/Fig3 Eiso-tail audit into a narrower closure
record and records where exact same-cluster/stored-vs-recomputed proof is
blocked by missing diagnostic fields in the current final ROOT.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[3]
PRIOR = ROOT / "dataOutput/ppg12Parity/control_plane/audits/the76_fig24_fig3_eiso_tail_20260706"
OUT = ROOT / "dataOutput/ppg12Parity/control_plane/audits/the76_eiso_contract_closure_20260706"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fieldnames = sorted({k for row in rows for k in row})
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_md(path: Path, lines: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines).rstrip() + "\n")


def as_float(row: dict[str, str], key: str) -> float:
    try:
        return float(row[key])
    except Exception:
        return float("nan")


def make_contract_outputs() -> tuple[dict[str, Any], dict[str, Any]]:
    ppg12_rows = read_csv(PRIOR / "ppg12_source_chain_table.csv")
    current_rows = read_csv(PRIOR / "current_source_chain_table.csv")

    ppg12 = {
        "status": "frozen_from_prior_audit",
        "rows": ppg12_rows,
        "summary": [
            "Fig24 PPG12 source is /sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_bdt_nom.root:h_singal_reco_isoET_0.",
            "Fig24 observable is corrected_Eiso = 1.2 * cluster_iso_topo_04 + 0.1 for MC nominal config.",
            "Fig24 cutoff macro FindETCut.C scans the TH2; prior audit ruled out cutoff extraction as the cause.",
            "Fig3 data templates use the same efficiency-tool isolation convention, without the MC scale/shift branch.",
        ],
        "evidence": [
            "ppg12codeGit/efficiencytool/config_bdt_nom.yaml:61,81-83",
            "ppg12codeGit/efficiencytool/RecoEffCalculator_TTreeReader.C:690,1213-1215,2151-2164,2261-2266,2916",
            "ppg12codeGit/efficiencytool/FindETCut.C:7-105",
            str(PRIOR / "ppg12_source_chain_table.csv"),
        ],
    }
    current = {
        "status": "frozen_from_prior_audit_plus_current_root_tree_check",
        "rows": current_rows,
        "summary": [
            "Fig24 current source is the current THE-76 photonjet merged ROOT: SIM/h_singal_reco_isoET_0.",
            "RecoilJets calls ppg12PhotonYieldRawEiso(recoMatch, topNode), then ppg12PhotonYieldEiso(raw).",
            "Raw path first attempts PhotonClusterv1 ppg12_topo_valid_04 / ppg12_topo_raw_eiso_04 with vertex compatibility, then recomputes a topo-cone sum if that fails.",
            "The current final merged ROOT has no AuAuPhotonIDTrainingTree, so the exact cluster-level raw path cannot be recovered from this artifact.",
        ],
        "evidence": [
            "src/RecoilJets.cc:9157,9180-9244",
            "src/RecoilJets.cc:13487-13627",
            "src/RecoilJets.h:1672-1674",
            "ROOT readback: final merged ROOT has only SIM directory and no AuAuPhotonIDTrainingTree",
            str(PRIOR / "current_source_chain_table.csv"),
        ],
    }

    write_json(OUT / "ppg12_eiso_contract.json", ppg12)
    write_md(
        OUT / "ppg12_eiso_contract.md",
        [
            "# PPG12 Eiso Contract",
            "",
            *[f"- {x}" for x in ppg12["summary"]],
            "",
            "## Evidence",
            *[f"- `{x}`" for x in ppg12["evidence"]],
        ],
    )
    write_json(OUT / "recoiljets_eiso_contract.json", current)
    write_md(
        OUT / "recoiljets_eiso_contract.md",
        [
            "# RecoilJets Eiso Contract",
            "",
            *[f"- {x}" for x in current["summary"]],
            "",
            "## Evidence",
            *[f"- `{x}`" for x in current["evidence"]],
        ],
    )
    return ppg12, current


def make_path_usage_outputs() -> dict[str, Any]:
    rows = [
        {
            "question": "clusters using stored ppg12_topo_raw_eiso_04",
            "status": "blocked_missing_path_label",
            "available_from_current_root": "no",
            "evidence": "current final ROOT has no AuAuPhotonIDTrainingTree and no stored/recomputed split histogram",
            "consequence": "cannot prove stored-only agreement or fallback-driven tail from the merged artifact",
        },
        {
            "question": "clusters using recomputed topo-cone fallback",
            "status": "blocked_missing_path_label",
            "available_from_current_root": "no",
            "evidence": "ppg12PhotonYieldRawEiso returns only a double; no per-cluster enum/counter is persisted in this ROOT",
            "consequence": "fallback/recompute remains the leading controllable suspect, not a proven cause",
        },
        {
            "question": "Eiso shape from final mixed path",
            "status": "available",
            "available_from_current_root": "yes",
            "evidence": str(PRIOR / "fig24_eiso_shape_comparison_by_pt.csv"),
            "consequence": "final mixed current distribution is broader/harder at high pT",
        },
    ]
    write_csv(OUT / "stored_vs_recomputed_eiso_canary.csv", rows)
    payload = {
        "classification": "blocked_missing_path_label",
        "rows": rows,
        "minimal_unblock": [
            "Add diagnostic-only stored/recomputed path counters or a compact per-cluster canary output.",
            "For each Fig24-matched cluster record stored valid flag, vertex compatibility, stored raw Eiso, recomputed raw Eiso, final raw Eiso used, corrected Eiso, pT bin, sample, and event key.",
            "Then rerun a tiny foreground photon20 high-pT canary before any broad production.",
        ],
    }
    write_json(OUT / "stored_vs_recomputed_eiso_canary.json", payload)
    return payload


def make_same_cluster_outputs() -> dict[str, Any]:
    rows = [
        {
            "test": "same-cluster PPG12 cluster_iso_topo_04 vs RecoilJets stored ppg12_topo_raw_eiso_04",
            "status": "blocked",
            "blocker": "current merged ROOT has no per-cluster tree and no stored raw Eiso branch",
            "smallest_next_check": "tiny foreground RecoilJets canary on shared photon20 high-pT events with diagnostic tree/csv; compare to PPG12 slimTree cluster_iso_topo_04 by event, truth track, eta/phi",
        },
        {
            "test": "same-cluster PPG12 cluster_iso_topo_04 vs RecoilJets recomputed topo-cone isolation",
            "status": "blocked",
            "blocker": "recomputed raw value is not persisted; current code only returns final raw double to the histogram fill",
            "smallest_next_check": "instrument standalone canary or diagnostic-only helper to print recomputed topo sum before correction",
        },
        {
            "test": "same-cluster corrected Eiso comparison",
            "status": "blocked",
            "blocker": "raw branch/source is blocked; corrected comparison alone would not identify the mechanism",
            "smallest_next_check": "same canary must print raw and corrected values plus correction application count",
        },
    ]
    write_csv(OUT / "same_cluster_eiso_canary.csv", rows)
    payload = {
        "classification": "blocked_missing_same_cluster_artifact",
        "rows": rows,
        "available_current_root_tree_status": "NO_TREE for AuAuPhotonIDTrainingTree in current final merged ROOT",
    }
    write_json(OUT / "same_cluster_eiso_canary.json", payload)
    return payload


def make_distribution_surrogate() -> dict[str, Any]:
    comparison = read_csv(PRIOR / "fig24_eiso_shape_comparison_by_pt.csv")
    raw_rows = read_csv(PRIOR / "fig24_current_component_and_raw_localization.csv")
    surrogate_rows: list[dict[str, Any]] = []
    for row in comparison:
        pt_low = as_float(row, "pt_low")
        pt_high = as_float(row, "pt_high")
        surrogate_rows.append(
            {
                "source": "fig24_final_mixed_ppg12_vs_current",
                "pt_low": pt_low,
                "pt_high": pt_high,
                "ppg12_mean": as_float(row, "ppg12_mean"),
                "current_mean": as_float(row, "current_mean"),
                "current_minus_ppg12_mean": as_float(row, "current_minus_ppg12_mean"),
                "ppg12_q70": as_float(row, "ppg12_cut70"),
                "current_q70": as_float(row, "current_cut70"),
                "sdcc_over_current_q70": as_float(row, "sdcc_over_current_cut70"),
                "ppg12_q80": as_float(row, "ppg12_cut80"),
                "current_q80": as_float(row, "current_cut80"),
                "sdcc_over_current_q80": as_float(row, "sdcc_over_current_cut80"),
                "ppg12_q90": as_float(row, "ppg12_cut90"),
                "current_q90": as_float(row, "current_cut90"),
                "sdcc_over_current_q90": as_float(row, "sdcc_over_current_cut90"),
                "current_minus_ppg12_cdf_le_1": as_float(row, "current_minus_ppg12_cdf_le_1"),
                "current_minus_ppg12_cdf_le_2": as_float(row, "current_minus_ppg12_cdf_le_2"),
                "interpretation": "current is harder/broader; fewer events pass low Eiso thresholds at high pT",
            }
        )

    for row in raw_rows:
        if row.get("source") in {"raw_all", "raw_tight", "raw_nonTight"}:
            surrogate_rows.append(
                {
                    "source": row.get("source"),
                    "pt_low": 26.0,
                    "pt_high": 35.0,
                    "entries": as_float(row, "entries"),
                    "mean_raw": as_float(row, "mean"),
                    "q70_raw": as_float(row, "q70"),
                    "q80_raw": as_float(row, "q80"),
                    "q90_raw": as_float(row, "q90"),
                    "q70_after_1p2_plus_0p1": 1.2 * as_float(row, "q70") + 0.1,
                    "q80_after_1p2_plus_0p1": 1.2 * as_float(row, "q80") + 0.1,
                    "q90_after_1p2_plus_0p1": 1.2 * as_float(row, "q90") + 0.1,
                    "interpretation": "raw current distribution exists but is not split by stored vs recomputed source",
                }
            )

    write_csv(OUT / "eiso_distribution_surrogate.csv", surrogate_rows)
    payload = {
        "classification": "distribution_surrogate_confirms_wrong_final_shape_but_not_mechanism",
        "rows": surrogate_rows,
        "high_pt_summary": [
            r for r in surrogate_rows
            if r.get("source") == "fig24_final_mixed_ppg12_vs_current" and r.get("pt_low") == 34.0
        ],
    }
    write_json(OUT / "eiso_distribution_surrogate.json", payload)
    return payload


def make_fig3_outputs() -> dict[str, Any]:
    rows = read_csv(PRIOR / "fig3_data_tail_summary.csv")
    out_rows = []
    for row in rows:
        out_rows.append(
            {
                **row,
                "path_split_available": "no",
                "interpretation": (
                    "Fig3 data tail is also high in current/PPG12, supporting a shared isolation-path issue; "
                    "it still does not identify stored-vs-recomputed without path labels."
                ),
            }
        )
    write_csv(OUT / "fig3_data_eiso_path_canary.csv", out_rows)
    payload = {
        "classification": "data_tail_supports_shared_isolation_path_issue_but_path_split_blocked",
        "rows": out_rows,
        "evidence": str(PRIOR / "fig3_data_tail_summary.csv"),
    }
    write_json(OUT / "fig3_data_eiso_path_canary.json", payload)
    return payload


def make_final(ppg12: dict[str, Any], current: dict[str, Any], path_usage: dict[str, Any], same_cluster: dict[str, Any], surrogate: dict[str, Any], fig3: dict[str, Any]) -> dict[str, Any]:
    final = {
        "current_state": "Fig24/Fig3 high-isolation discrepancy is real in filled Eiso shapes; exact mechanism is not proven from current final artifacts.",
        "classification": "unexplained: exact same-cluster evidence blocked",
        "first_causal_mechanism": "not proven; leading controllable suspect is recomputed topo-isolation fallback versus PPG12 cluster_iso_topo_04",
        "ruled_out": [
            "FindETCut.C / cutoff extraction as the cause",
            "documented linear MC correction mismatch as the leading cause",
            "no-vtx/coarse variant as sufficient explanation",
            "pure photon+jet stitching/weighting as the leading proven cause for the shared Fig3/Fig24 tail",
        ],
        "ambiguous": [
            "Whether stored ppg12_topo_raw_eiso_04 alone matches PPG12 cluster_iso_topo_04",
            "Whether recomputed topo-cone fallback is tail-heavy and grows at high pT",
            "Whether PhotonClusterBuilder stored branch itself differs from PPG12 cluster_iso_topo_04",
            "Whether cluster selection differences remain after raw Eiso identity is checked",
        ],
        "exact_controllable_fix_or_rerun_needed": [
            "Do not broad-rerun yet.",
            "Add a diagnostic-only foreground canary or temporary split histograms that record stored-vs-recomputed path choice and both raw values for Fig24-matched clusters.",
            "If stored-only matches PPG12 and recomputed is tail-heavy, fix by enforcing/repairing the PhotonClusterBuilder PPG12 stored topo isolation branch and avoiding silent recompute for PPG12 parity outputs.",
            "If stored differs from PPG12, fix PhotonClusterBuilder branch construction or upstream topo-cluster contract.",
            "If raw values match but selected clusters differ, fix truth-match/cluster-selection contract.",
        ],
        "validation_target_after_fix": [
            "stored_vs_recomputed_eiso_canary shows which subpath caused the high-pT tail",
            "same_cluster_eiso_canary has matched PPG12 cluster_iso_topo_04 and RecoilJets raw isolation values for photon20 high-pT clusters",
            "Fig24 h_singal_reco_isoET_0 cutoff ratios close bin-by-bin without the 34-36 GeV 0.87-0.89 SDCC/current drop",
            "Fig3 data template tail ratio no longer has the same high-isolation excess if the same path is causal",
        ],
        "canonicalization_allowed_now": False,
        "artifacts": {
            "ppg12_eiso_contract": str(OUT / "ppg12_eiso_contract.md"),
            "recoiljets_eiso_contract": str(OUT / "recoiljets_eiso_contract.md"),
            "stored_vs_recomputed": str(OUT / "stored_vs_recomputed_eiso_canary.csv"),
            "same_cluster": str(OUT / "same_cluster_eiso_canary.csv"),
            "distribution_surrogate": str(OUT / "eiso_distribution_surrogate.csv"),
            "fig3_data_linkage": str(OUT / "fig3_data_eiso_path_canary.csv"),
        },
    }
    write_json(OUT / "final_classification.json", final)
    write_md(
        OUT / "final_classification.md",
        [
            "# THE-76 Eiso Contract Closure",
            "",
            f"- Classification: `{final['classification']}`",
            f"- First causal mechanism: {final['first_causal_mechanism']}",
            f"- Canonicalization allowed now: `{final['canonicalization_allowed_now']}`",
            "",
            "## Ruled Out",
            *[f"- {x}" for x in final["ruled_out"]],
            "",
            "## Ambiguous",
            *[f"- {x}" for x in final["ambiguous"]],
            "",
            "## Exact Controllable Fix / Rerun Needed",
            *[f"- {x}" for x in final["exact_controllable_fix_or_rerun_needed"]],
            "",
            "## Validation Target After Fix",
            *[f"- {x}" for x in final["validation_target_after_fix"]],
        ],
    )
    return final


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    ppg12, current = make_contract_outputs()
    path_usage = make_path_usage_outputs()
    same_cluster = make_same_cluster_outputs()
    surrogate = make_distribution_surrogate()
    fig3 = make_fig3_outputs()
    final = make_final(ppg12, current, path_usage, same_cluster, surrogate, fig3)
    write_json(
        OUT / "audit_index.json",
        {
            "audit": "the76_eiso_contract_closure_20260706",
            "prior_audit": str(PRIOR),
            "outputs": final["artifacts"] | {
                "final_classification": str(OUT / "final_classification.md"),
                "audit_index": str(OUT / "audit_index.json"),
            },
        },
    )
    print(json.dumps({"out": str(OUT), "classification": final["classification"]}, indent=2))


if __name__ == "__main__":
    main()
