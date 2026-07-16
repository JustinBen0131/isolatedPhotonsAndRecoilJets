#!/usr/bin/env python3
"""Probe PPG12 BDT model-selection effects on RecoilJets PPG12-mode rows.

This is a local forensic helper for THE-76.  It reads an existing
AuAuPhotonIDTrainingTree produced by RecoilJets with PPG12 diagnostic branches,
re-evaluates the PPG12 base_E/base_v3E TMVA models from the stored BDT inputs,
and asks whether PPG12's ET-binned score branch choice can explain a deficit in
the non-tight sideband population.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any


PT_BINS = [10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 32.0, 36.0]


@dataclass(frozen=True)
class ThresholdSpec:
    name: str
    tight_intercept: float
    tight_slope: float
    nontight_min_intercept: float
    nontight_min_slope: float
    nontight_max_intercept: float
    nontight_max_slope: float

    def tight_min(self, et: float) -> float:
        return self.tight_intercept + self.tight_slope * et

    def nontight_min(self, et: float) -> float:
        return self.nontight_min_intercept + self.nontight_min_slope * et

    def nontight_max(self, et: float) -> float:
        return self.nontight_max_intercept + self.nontight_max_slope * et


THRESHOLD_SPECS = {
    # Constants currently mirrored from ppg12codeGit/efficiencytool/config_bdt_nom.yaml.
    "code_config": ThresholdSpec(
        name="code_config",
        tight_intercept=0.815625,
        tight_slope=-0.0015625,
        nontight_min_intercept=0.7333333333333333,
        nontight_min_slope=-0.01333333333333333,
        nontight_max_intercept=0.684375,
        nontight_max_slope=0.0015625,
    ),
    # Canonical PPG12 paper-era cuts; keep the exact fractions used by RecoilJets.
    "paper_canonical": ThresholdSpec(
        name="paper_canonical",
        tight_intercept=0.815625,
        tight_slope=-0.0015625,
        nontight_min_intercept=0.7333333333333333,
        nontight_min_slope=-0.01333333333333333,
        nontight_max_intercept=0.684375,
        nontight_max_slope=0.0015625,
    ),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rj-root", required=True, help="RecoilJets diagnostic ROOT file")
    parser.add_argument("--tree", default="AuAuPhotonIDTrainingTree", help="Tree name")
    parser.add_argument("--base-e-model", required=True, help="PPG12 base_E TMVA ROOT model")
    parser.add_argument("--base-v3e-model", required=True, help="PPG12 base_v3E TMVA ROOT model")
    parser.add_argument("--out-md", required=True, help="Markdown report path")
    parser.add_argument("--out-csv", required=True, help="CSV report path")
    parser.add_argument(
        "--max-entries",
        type=int,
        default=-1,
        help="Optional maximum tree entries to process",
    )
    return parser.parse_args()


def pt_bin_index(pt: float) -> int | None:
    for idx in range(len(PT_BINS) - 1):
        if PT_BINS[idx] <= pt < PT_BINS[idx + 1]:
            return idx
    return None


def in_range(value: float, low: float, high: float) -> bool:
    return math.isfinite(value) and low < value < high


def classify_with_score(row: Any, score: float, thresholds: ThresholdSpec) -> str:
    """Replicate the PPG12 tight/non-tight tag using stored RJ shape inputs.

    The foreground diagnostic tree does not store cluster_prob.  This helper
    only evaluates rows that already have ppg12_common_pass=true, so the common
    probability gate has already passed and the prob nfail leg is inert.
    """

    if not bool(getattr(row, "ppg12_common_pass")):
        return "not_common"

    et = float(getattr(row, "cluster_Et"))
    weta = float(getattr(row, "cluster_weta_cogx"))
    wphi = float(getattr(row, "cluster_wphi_cogx"))
    et1 = float(getattr(row, "cluster_et1"))
    et2 = float(getattr(row, "cluster_et2"))
    et3 = float(getattr(row, "cluster_et3"))
    et4 = float(getattr(row, "cluster_et4"))
    e11_e33 = float(getattr(row, "e11_over_e33"))
    e32_e35 = float(getattr(row, "e32_over_e35"))

    tight_shape = (
        in_range(weta, 0.0, 1.0)
        and in_range(wphi, 0.0, 1.0)
        and in_range(et1, 0.5, 1.0)
        and in_range(et2, 0.0, 1.0)
        and in_range(et3, 0.0, 1.0)
        and in_range(et4, 0.0, 1.0)
        and in_range(e11_e33, 0.0, 1.0)
        and in_range(e32_e35, 0.8, 1.0)
    )
    tight_bdt = score > thresholds.tight_min(et)
    if tight_shape and tight_bdt:
        return "tight"

    nontight_shape = (
        in_range(weta, 0.0, 1.0)
        and in_range(wphi, 0.0, 1.0)
        and in_range(et1, 0.6, 1.0)
        and in_range(et4, 0.0, 1.0)
        and in_range(e11_e33, 0.0, 1.0)
        and in_range(e32_e35, 0.8, 1.0)
    )
    low, high = thresholds.nontight_min(et), thresholds.nontight_max(et)
    nontight_bdt = low < score < high
    nfail = 0
    if not in_range(weta, 0.0, 1.0):
        nfail += 1
    if not tight_bdt:
        nfail += 1
    if nontight_shape and nontight_bdt and nfail > 0:
        return "nontight"
    return "neither"


def load_rbdt(model_path: str):
    import ROOT  # type: ignore

    return ROOT.TMVA.Experimental.RBDT("myBDT", model_path)


def eval_rbdt(model: Any, feature_names: list[str], values: dict[str, float]) -> float:
    import ROOT  # type: ignore

    vec = ROOT.std.vector("float")()
    for name in feature_names:
        vec.push_back(float(values[name]))
    result = model.Compute(vec)
    return float(result[0]) if len(result) else float("nan")


def ratio(num: float, den: float) -> float:
    return num / den if den else float("nan")


def main() -> int:
    args = parse_args()

    import ROOT  # type: ignore

    ROOT.gROOT.SetBatch(True)
    ROOT.TMVA.Tools.Instance()

    base_e_tail_features = [
        "vertexz",
        "cluster_Eta",
        "e11_over_e33",
        "cluster_et1",
        "cluster_et2",
        "cluster_et3",
        "cluster_et4",
    ]
    base_v3e_tail_features = [
        "cluster_weta_cogx",
        "cluster_wphi_cogx",
        "vertexz",
        "cluster_Eta",
        "e11_over_e33",
        "cluster_et1",
        "cluster_et2",
        "cluster_et3",
        "cluster_et4",
        "e32_over_e35",
    ]

    base_e_model = load_rbdt(args.base_e_model)
    base_v3e_model = load_rbdt(args.base_v3e_model)

    root_file = ROOT.TFile.Open(args.rj_root)
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"Could not open {args.rj_root}")
    tree = root_file.Get(args.tree)
    if not tree:
        raise RuntimeError(f"Could not find tree {args.tree} in {args.rj_root}")

    branches = {branch.GetName() for branch in tree.GetListOfBranches()}
    score_et_branch = "cluster_Et_score_input" if "cluster_Et_score_input" in branches else "cluster_Et"
    base_e_features = [score_et_branch] + base_e_tail_features
    base_v3e_features = [score_et_branch] + base_v3e_tail_features
    needed = set(base_e_features + base_v3e_features)
    needed.update(["cluster_Et", "is_signal", "ppg12_common_pass", "ppg12_tight_tag", "tight_bdt_score"])
    missing = sorted(needed - branches)
    if missing:
        raise RuntimeError(f"Missing required branches: {missing}")

    rows: dict[int, dict[str, float]] = defaultdict(lambda: defaultdict(float))
    flips: dict[int, dict[str, float]] = defaultdict(lambda: defaultdict(float))
    max_entries = tree.GetEntries() if args.max_entries < 0 else min(args.max_entries, tree.GetEntries())

    for i, row in enumerate(tree):
        if i >= max_entries:
            break
        if hasattr(row, "is_signal") and not bool(getattr(row, "is_signal")):
            continue
        et = float(getattr(row, "cluster_Et"))
        bin_idx = pt_bin_index(et)
        if bin_idx is None:
            continue

        values = {name: float(getattr(row, name)) for name in needed if hasattr(row, name)}
        base_e_score = eval_rbdt(base_e_model, base_e_features, values)
        base_v3e_score = eval_rbdt(base_v3e_model, base_v3e_features, values)
        stored_score = float(getattr(row, "tight_bdt_score"))
        selected_model = "base_v3E" if 8.0 <= et < 35.0 else "base_E"
        selected_score = base_v3e_score if selected_model == "base_v3E" else base_e_score

        stored_tag_code = int(getattr(row, "ppg12_tight_tag"))
        stored_tag = {1: "tight", 2: "nontight", 3: "neither"}.get(stored_tag_code, f"tag{stored_tag_code}")
        classified = {
            f"stored_score_{name}": classify_with_score(row, stored_score, spec)
            for name, spec in THRESHOLD_SPECS.items()
        }
        classified.update(
            {
                f"selected_score_{name}": classify_with_score(row, selected_score, spec)
                for name, spec in THRESHOLD_SPECS.items()
            }
        )

        bucket = rows[bin_idx]
        bucket["signal"] += 1.0
        bucket["score_stored_sum"] += stored_score
        bucket["score_base_e_sum"] += base_e_score
        bucket["score_base_v3e_sum"] += base_v3e_score
        bucket["score_selected_sum"] += selected_score
        bucket[f"selected_model_{selected_model}"] += 1.0
        if bool(getattr(row, "ppg12_common_pass")):
            bucket["common"] += 1.0
            bucket[f"stored_{stored_tag}"] += 1.0
            for label, tag in classified.items():
                bucket[f"{label}_{tag}"] += 1.0
                if tag != stored_tag:
                    flips[bin_idx][f"stored_{stored_tag}_vs_{label}_{tag}"] += 1.0

    records: list[dict[str, Any]] = []
    for bin_idx in range(len(PT_BINS) - 1):
        bucket = rows[bin_idx]
        common = bucket["common"]
        record: dict[str, Any] = {
            "pt_bin": f"{PT_BINS[bin_idx]:.0f}-{PT_BINS[bin_idx + 1]:.0f}",
            "signal": int(bucket["signal"]),
            "common": int(common),
            "stored_tight_common": ratio(bucket["stored_tight"], common),
            "stored_nontight_common": ratio(bucket["stored_nontight"], common),
            "stored_neither_common": ratio(bucket["stored_neither"], common),
            "stored_score_code_config_tight_common": ratio(bucket["stored_score_code_config_tight"], common),
            "stored_score_code_config_nontight_common": ratio(bucket["stored_score_code_config_nontight"], common),
            "stored_score_code_config_neither_common": ratio(bucket["stored_score_code_config_neither"], common),
            "stored_score_paper_canonical_tight_common": ratio(bucket["stored_score_paper_canonical_tight"], common),
            "stored_score_paper_canonical_nontight_common": ratio(bucket["stored_score_paper_canonical_nontight"], common),
            "stored_score_paper_canonical_neither_common": ratio(bucket["stored_score_paper_canonical_neither"], common),
            "selected_score_code_config_tight_common": ratio(bucket["selected_score_code_config_tight"], common),
            "selected_score_code_config_nontight_common": ratio(bucket["selected_score_code_config_nontight"], common),
            "selected_score_code_config_neither_common": ratio(bucket["selected_score_code_config_neither"], common),
            "selected_score_paper_canonical_tight_common": ratio(bucket["selected_score_paper_canonical_tight"], common),
            "selected_score_paper_canonical_nontight_common": ratio(bucket["selected_score_paper_canonical_nontight"], common),
            "selected_score_paper_canonical_neither_common": ratio(bucket["selected_score_paper_canonical_neither"], common),
            "mean_stored_score": ratio(bucket["score_stored_sum"], bucket["signal"]),
            "mean_base_e_score": ratio(bucket["score_base_e_sum"], bucket["signal"]),
            "mean_base_v3e_score": ratio(bucket["score_base_v3e_sum"], bucket["signal"]),
            "mean_selected_score": ratio(bucket["score_selected_sum"], bucket["signal"]),
            "selected_base_e_rows": int(bucket["selected_model_base_E"]),
            "selected_base_v3e_rows": int(bucket["selected_model_base_v3E"]),
        }
        for key, value in sorted(flips[bin_idx].items()):
            record[key] = int(value)
        records.append(record)

    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = sorted({key for record in records for key in record})
    preferred = [
        "pt_bin",
        "signal",
        "common",
        "stored_nontight_common",
        "stored_score_code_config_nontight_common",
        "stored_score_paper_canonical_nontight_common",
        "selected_score_code_config_nontight_common",
        "selected_score_paper_canonical_nontight_common",
        "stored_tight_common",
        "stored_score_code_config_tight_common",
        "stored_score_paper_canonical_tight_common",
        "selected_score_code_config_tight_common",
        "selected_score_paper_canonical_tight_common",
        "stored_neither_common",
        "stored_score_code_config_neither_common",
        "stored_score_paper_canonical_neither_common",
        "selected_score_code_config_neither_common",
        "selected_score_paper_canonical_neither_common",
        "mean_stored_score",
        "mean_base_e_score",
        "mean_base_v3e_score",
        "mean_selected_score",
        "selected_base_e_rows",
        "selected_base_v3e_rows",
    ]
    fieldnames = preferred + [name for name in fieldnames if name not in preferred]
    with out_csv.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)

    out_md = Path(args.out_md)
    out_md.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# THE-76 RJ PPG12 BDT Model-Selection Probe",
        "",
        f"- RecoilJets ROOT: `{args.rj_root}`",
        f"- Processed entries: `{max_entries}`",
        f"- base_E model: `{args.base_e_model}`",
        f"- base_v3E model: `{args.base_v3e_model}`",
        f"- Score ET branch used for TMVA recomputation: `{score_et_branch}`",
        "- Cut/model-selection ET branch: `cluster_Et`",
        "",
        "This re-evaluates PPG12 `base_E` and `base_v3E` scores from the stored RecoilJets feature branches.",
        "The selected-score tag uses PPG12's RecoEff ET-bin rule on the cut ET: `base_v3E` for `8 <= ET < 35`, otherwise `base_E`.",
        "It then verifies the stored and recomputed tags against the canonical PPG12 paper-era BDT thresholds.",
        "",
        "| pT bin | common | actual RJ NT/common | stored score + code NT/common | stored score + canonical NT/common | recomputed score + code NT/common | recomputed score + canonical NT/common | actual RJ tight/common | stored score + canonical tight/common | recomputed score + canonical tight/common | selected base_E rows |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for record in records:
        lines.append(
            "| {pt_bin} | {common} | {stored_nt:.4f} | {stored_code_nt:.4f} | {stored_ian_nt:.4f} | {sel_code_nt:.4f} | {sel_ian_nt:.4f} | {stored_t:.4f} | {stored_ian_t:.4f} | {sel_ian_t:.4f} | {base_e} |".format(
                pt_bin=record["pt_bin"],
                common=record["common"],
                stored_nt=record["stored_nontight_common"],
                stored_code_nt=record["stored_score_code_config_nontight_common"],
                stored_ian_nt=record["stored_score_paper_canonical_nontight_common"],
                sel_code_nt=record["selected_score_code_config_nontight_common"],
                sel_ian_nt=record["selected_score_paper_canonical_nontight_common"],
                stored_t=record["stored_tight_common"],
                stored_ian_t=record["stored_score_paper_canonical_tight_common"],
                sel_ian_t=record["selected_score_paper_canonical_tight_common"],
                base_e=record["selected_base_e_rows"],
            )
        )
    lines.extend(
        [
            "",
            "## Interpretation Guard",
            "",
            "- This is a local foreground probe, not final production evidence.",
            "- RBDT recomputation uses raw score-input ET when the diagnostic tree stores it; threshold and model-selection cuts use `cluster_Et`.",
            "- Because the probe starts from `ppg12_common_pass`, the probability leg is treated as already passed for common rows.",
            "- Movement from actual RJ/stored-score code thresholds to stored-score IAN thresholds isolates the threshold formula.",
            "- Additional movement from stored-score to recomputed-score columns isolates score-input/model-evaluation differences.",
        ]
    )
    out_md.write_text("\n".join(lines) + "\n")

    print(f"Wrote {out_md}")
    print(f"Wrote {out_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
