#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
from pathlib import Path
import argparse

import numpy as np


REPORT = Path(
    "/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/mlp_models/"
    "global_etcent_inclusive3_sixpack_20260516_135439/bdt_ablation_sidecars/"
    "basev3e_w33_e22ratio_20260518_192630/validation"
)
MODEL_DIR = REPORT.parent / "bdt_models"


def read_json_with_preamble(path: Path) -> dict:
    text = path.read_text()
    return json.loads(text[text.find("{") :])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--product", default="baseBDT_v3E_withCentrality_w33_E22E37")
    args = parser.parse_args()

    meta = json.loads((MODEL_DIR / f"auau_tight_bdt_{args.product}_tmva.metadata.json").read_text())
    model = read_json_with_preamble(MODEL_DIR / f"auau_tight_bdt_{args.product}_tmva.xgb.json")
    features = list(meta["features"])
    ratio_features = [f for f in features if f.startswith("e") and "_over_" in f]

    gain = {f: 0.0 for f in features}
    count = {f: 0 for f in features}
    for tree in model["learner"]["gradient_booster"]["model"]["trees"]:
        for left, right, idx, loss in zip(
            tree["left_children"],
            tree["right_children"],
            tree["split_indices"],
            tree["loss_changes"],
        ):
            if left < 0 and right < 0:
                continue
            idx = int(idx)
            if 0 <= idx < len(features):
                feature = features[idx]
                gain[feature] += max(float(loss), 0.0)
                count[feature] += 1

    ratio_features.sort(key=lambda f: gain.get(f, 0.0), reverse=True)
    cent_bins = [(0.0, 20.0, "0-20%"), (20.0, 50.0, "20-50%"), (50.0, 80.0, "50-80%")]
    edges = np.linspace(0.0, 1.2, 61)
    rows = []
    cache_paths = [
        Path(line.strip())
        for line in (REPORT / "score_caches.list").read_text().splitlines()
        if line.strip()
    ]

    for feature in ratio_features:
        sums = {
            (cent_label, cls): np.zeros(len(edges) - 1, dtype=np.float64)
            for _, _, cent_label in cent_bins
            for cls in (0, 1)
        }
        totals = {(cent_label, cls): 0 for _, _, cent_label in cent_bins for cls in (0, 1)}
        for path in cache_paths:
            with np.load(path, allow_pickle=True) as data:
                if feature not in data:
                    continue
                values = data[feature].astype(np.float64)
                centrality = data["centrality"].astype(np.float64)
                labels = data["is_signal"].astype(np.int32)
                finite = np.isfinite(values) & np.isfinite(centrality) & np.isin(labels, [0, 1])
                for c_lo, c_hi, cent_label in cent_bins:
                    cent_mask = finite & (centrality >= c_lo) & (centrality < c_hi)
                    for cls in (0, 1):
                        selected = values[cent_mask & (labels == cls)]
                        hist, _ = np.histogram(selected, bins=edges)
                        sums[(cent_label, cls)] += hist
                        totals[(cent_label, cls)] += int(selected.size)

        for cent_label in [label for _, _, label in cent_bins]:
            for cls in (0, 1):
                hist = sums[(cent_label, cls)]
                total = totals[(cent_label, cls)]
                density = hist / max(1, total) / np.diff(edges)
                for low, high, den, cnt in zip(edges[:-1], edges[1:], density, hist):
                    rows.append(
                        {
                            "feature": feature,
                            "feature_order_by_gain": ratio_features.index(feature) + 1,
                            "split_gain": gain.get(feature, 0.0),
                            "split_count": count.get(feature, 0),
                            "cent_bin": cent_label,
                            "class": "signal" if cls == 1 else "background",
                            "bin_low": float(low),
                            "bin_high": float(high),
                            "density": float(den),
                            "count": int(cnt),
                            "total_entries": int(total),
                        }
                    )

    suffix = args.product.removeprefix("baseBDT_v3E_withCentrality_w33_")
    out = REPORT / f"energy_ratio_histograms_{suffix}_fullstat.csv"
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(out)


if __name__ == "__main__":
    main()
