#!/usr/bin/env python3
"""Train JetML residual-pT regressors from RecoilJets_AuAu row trees.

The expected input tree is ``JetResidualMLTrainingTree``.  One row is one
inclusive recoil-jet candidate from a selected gamma+jet event.  The regression
target is

    delta_pt = truth_pt - reco_areaSub_pt

using only truth-matched recoil jets for training.  Unmatched rows remain useful
in the ROOT files for fake/combinatoric QA, but they are intentionally excluded
from the residual regressor.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Iterable


CORE_FEATURES = [
    "reco_areaSub_pt",
    "raw_pt",
    "jet_area",
    "jet_eta",
    "centrality",
    "R",
]


def expand_input_paths(items: list[Path]) -> list[Path]:
    paths: list[Path] = []
    for item in items:
        text = str(item)
        if text.startswith("@"):
            manifest = Path(text[1:])
            if not manifest.is_file():
                raise SystemExit(f"Input manifest does not exist: {manifest}")
            for raw in manifest.read_text().splitlines():
                line = raw.strip()
                if line and not line.startswith("#"):
                    paths.append(Path(line))
        else:
            paths.append(item)
    if not paths:
        raise SystemExit("No input ROOT files supplied")
    missing = [str(path) for path in paths if not path.is_file()]
    if missing:
        preview = "\n  ".join(missing[:20])
        extra = "" if len(missing) <= 20 else f"\n  ... {len(missing) - 20} more"
        raise SystemExit(f"Input ROOT files are missing:\n  {preview}{extra}")
    return paths


def finite_mask(frame, columns: Iterable[str]):
    import numpy as np

    mask = np.ones(len(frame), dtype=bool)
    for col in columns:
        mask &= np.isfinite(frame[col].to_numpy())
    return mask


def load_frame(paths: list[Path], tree_name: str, columns: list[str]):
    try:
        import uproot
        import pandas as pd
    except Exception as exc:
        raise SystemExit(
            "This trainer needs uproot and pandas in the active environment. "
            "Run it from the same Python env used for analysis training."
        ) from exc

    frames = []
    for path in paths:
        with uproot.open(path) as f:
            if tree_name not in f:
                raise SystemExit(f"{path} does not contain tree {tree_name}")
            tree = f[tree_name]
            missing = [col for col in columns if col not in tree.keys()]
            if missing:
                raise SystemExit(f"{path}:{tree_name} is missing columns: {', '.join(missing)}")
            frame = tree.arrays(columns, library="pd")
            frame["source_file"] = str(path)
            frames.append(frame)

    if not frames:
        raise SystemExit("No rows loaded")
    return pd.concat(frames, ignore_index=True)


def split_by_event(frame, seed: int):
    import numpy as np

    keys = (frame["run"].astype(str) + ":" + frame["evt"].astype(str)).to_numpy()
    unique = np.unique(keys)
    rng = np.random.default_rng(seed)
    rng.shuffle(unique)

    n = len(unique)
    n_train = int(0.6 * n)
    n_val = int(0.2 * n)
    train_keys = set(unique[:n_train])
    val_keys = set(unique[n_train:n_train + n_val])

    is_train = np.array([key in train_keys for key in keys], dtype=bool)
    is_val = np.array([key in val_keys for key in keys], dtype=bool)
    is_test = ~(is_train | is_val)
    return is_train, is_val, is_test


def train_model(frame, features: list[str], seed: int):
    import numpy as np
    from xgboost import XGBRegressor

    train_mask, val_mask, test_mask = split_by_event(frame, seed)
    x_train = frame.loc[train_mask, features].to_numpy(dtype=np.float32)
    y_train = frame.loc[train_mask, "delta_pt"].to_numpy(dtype=np.float32)
    x_val = frame.loc[val_mask, features].to_numpy(dtype=np.float32)
    y_val = frame.loc[val_mask, "delta_pt"].to_numpy(dtype=np.float32)

    model = XGBRegressor(
        n_estimators=500,
        max_depth=3,
        learning_rate=0.035,
        subsample=0.8,
        colsample_bytree=0.9,
        objective="reg:squarederror",
        reg_lambda=2.0,
        random_state=seed,
        n_jobs=4,
    )
    model.fit(
        x_train,
        y_train,
        eval_set=[(x_val, y_val)],
        verbose=False,
    )

    def rms(values):
        arr = np.asarray(values, dtype=np.float64)
        return float(math.sqrt(np.mean(arr * arr))) if len(arr) else float("nan")

    x_all = frame[features].to_numpy(dtype=np.float32)
    pred = model.predict(x_all)
    reco_resid = frame["reco_areaSub_pt"].to_numpy(dtype=np.float64) - frame["truth_pt"].to_numpy(dtype=np.float64)
    ml_resid = (frame["reco_areaSub_pt"].to_numpy(dtype=np.float64) + pred) - frame["truth_pt"].to_numpy(dtype=np.float64)

    metrics = {
        "rows_train": int(train_mask.sum()),
        "rows_val": int(val_mask.sum()),
        "rows_test": int(test_mask.sum()),
        "area_residual_mean": float(np.mean(reco_resid)),
        "area_residual_rms": rms(reco_resid),
        "ml_residual_mean": float(np.mean(ml_resid)),
        "ml_residual_rms": rms(ml_resid),
    }
    return model, metrics


def save_tmva(model, output: Path, features: list[str]):
    try:
        import ROOT
    except Exception as exc:
        raise SystemExit(
            "Could not import ROOT Python bindings. Use the analysis ROOT environment "
            "when exporting the TMVA/RBDT model."
        ) from exc

    output.parent.mkdir(parents=True, exist_ok=True)
    booster = model.get_booster()
    booster.feature_names = [f"f{i}" for i in range(len(features))]
    try:
        ROOT.TMVA.Experimental.SaveXGBoost(model, "myBDT", str(output), num_inputs=len(features))
    except (AttributeError, TypeError):
        ROOT.TMVA.Experimental.SaveXGBoost(booster, "myBDT", str(output), num_inputs=len(features))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", nargs="+", type=Path, required=True,
                        help="Input ROOT files, or @manifest.txt")
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--tree", default="JetResidualMLTrainingTree")
    parser.add_argument("--features", default=",".join(CORE_FEATURES),
                        help="Comma-separated feature order. Must match YAML inference order.")
    parser.add_argument("--model-name", default="jetResidualBDT_r03_allCent_tmva.root")
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--min-truth-pt", type=float, default=0.0)
    parser.add_argument("--cent-min", type=float, default=None)
    parser.add_argument("--cent-max", type=float, default=None)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    features = [item.strip() for item in args.features.split(",") if item.strip()]
    if not features:
        raise SystemExit("Feature list is empty")

    required = sorted(set(features + [
        "run", "evt", "is_recoil", "is_matched", "centrality",
        "reco_areaSub_pt", "truth_pt", "delta_pt",
    ]))
    paths = expand_input_paths(args.input)
    frame = load_frame(paths, args.tree, required)

    mask = (frame["is_recoil"].to_numpy() == 1) & (frame["is_matched"].to_numpy() == 1)
    mask &= frame["truth_pt"].to_numpy() > args.min_truth_pt
    mask &= finite_mask(frame, features + ["reco_areaSub_pt", "truth_pt", "delta_pt", "centrality"])
    if args.cent_min is not None:
        mask &= frame["centrality"].to_numpy() >= args.cent_min
    if args.cent_max is not None:
        mask &= frame["centrality"].to_numpy() < args.cent_max

    train_frame = frame.loc[mask].copy()
    if len(train_frame) < 100:
        raise SystemExit(f"Too few matched recoil rows after cuts: {len(train_frame)}")

    model, metrics = train_model(train_frame, features, args.seed)

    args.outdir.mkdir(parents=True, exist_ok=True)
    tmva_path = args.outdir / args.model_name
    json_model = tmva_path.with_suffix(".xgb.json")
    model.get_booster().save_model(json_model)
    save_tmva(model, tmva_path, features)

    (args.outdir / "featureOrder.txt").write_text("\n".join(features) + "\n")
    summary = {
        "input_files": [str(path) for path in paths],
        "tree": args.tree,
        "features": features,
        "target": "delta_pt = truth_pt - reco_areaSub_pt",
        "rows_loaded": int(len(frame)),
        "rows_used": int(len(train_frame)),
        "selection": {
            "is_recoil": 1,
            "is_matched": 1,
            "min_truth_pt": args.min_truth_pt,
            "cent_min": args.cent_min,
            "cent_max": args.cent_max,
        },
        "metrics": metrics,
        "tmva_model": str(tmva_path),
        "xgb_json_model": str(json_model),
    }
    (args.outdir / "trainingSummary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    sys.exit(main())
