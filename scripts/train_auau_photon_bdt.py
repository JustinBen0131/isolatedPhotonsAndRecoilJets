#!/usr/bin/env python3
"""Train AuAu photon BDT TMVA files from RecoilJets_AuAu training trees.

The expected input tree is ``AuAuPhotonIDTrainingTree``.  For tight photon-ID,
``is_signal`` labels truth-matched isolated prompt photons from embedded
photon+jet signal samples and non-signal clusters from embedded inclusive jet
background samples.  For NPB, this follows the PPG12 convention:
``npb_label=1`` means a physics-like cluster and ``npb_label=0`` means a
timing-tagged non-physics background cluster from data.  The resulting score is
therefore a physics-like score, so the analysis cut remains
``auau_npb_score > 0.5``.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Iterable


PPG12_TIGHT_FEATURES = [
    "cluster_Et",
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

PPG12_NPB_FEATURES = [
    "cluster_Et",
    "cluster_Eta",
    "vertexz",
    "e11_over_e33",
    "e32_over_e35",
    "e11_over_e22",
    "e11_over_e13",
    "e11_over_e15",
    "e11_over_e17",
    "e11_over_e31",
    "e11_over_e51",
    "e11_over_e71",
    "e22_over_e33",
    "e22_over_e35",
    "e22_over_e37",
    "e22_over_e53",
    "cluster_weta_cogx",
    "cluster_wphi_cogx",
    "cluster_et1",
    "cluster_et2",
    "cluster_et3",
    "cluster_et4",
    "cluster_w32",
    "cluster_w52",
    "cluster_w72",
]


def parse_cent_bins(text: str) -> list[tuple[float, float]]:
    bins: list[tuple[float, float]] = []
    if not text:
        return bins
    for item in text.split(","):
        lo_s, hi_s = item.split(":", 1)
        bins.append((float(lo_s), float(hi_s)))
    return bins


def cent_tag(lo: float, hi: float) -> str:
    return f"cent_{int(round(lo)):03d}_{int(round(hi)):03d}"


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


def load_frame(
    paths: list[Path],
    tree_name: str,
    required_columns: list[str],
    optional_columns: list[str],
    missing_label_branch: str | None,
    missing_label_value: int | None,
):
    try:
        import pandas as pd
        import uproot
    except ImportError as exc:
        raise SystemExit(
            "Missing dependency. Use an environment with uproot, pandas, numpy, "
            "scikit-learn, xgboost, and PyROOT for TMVA export."
        ) from exc

    frames = []
    seen_optional: set[str] = set()
    for path in paths:
        with uproot.open(path) as root_file:
            if tree_name not in root_file:
                raise SystemExit(f"{path} does not contain tree {tree_name}")
            tree = root_file[tree_name]
            keys = set(tree.keys())
            missing = [col for col in required_columns if col not in keys]
            allow_missing_label = (
                missing_label_branch is not None
                and missing == [missing_label_branch]
                and missing_label_value is not None
            )
            if allow_missing_label:
                missing = []
            if missing:
                raise SystemExit(f"{path}:{tree_name} missing required branches: {', '.join(missing)}")
            present_optional = [col for col in optional_columns if col in keys]
            seen_optional.update(present_optional)
            read_columns = [col for col in required_columns if col in keys] + present_optional
            frame = tree.arrays(read_columns, library="pd")
            if allow_missing_label:
                frame[missing_label_branch] = int(missing_label_value)
            frames.append(frame)

    frame = pd.concat(frames, ignore_index=True)
    for col in optional_columns:
        if col not in frame.columns:
            frame[col] = 1.0
    return frame, sorted(seen_optional)


def normalize_mean_one(values):
    import numpy as np

    values = np.asarray(values, dtype="float64")
    finite = np.isfinite(values) & (values > 0.0)
    if not finite.any():
        return np.ones_like(values, dtype="float64")
    mean = float(values[finite].mean())
    if mean <= 0.0 or not math.isfinite(mean):
        return np.ones_like(values, dtype="float64")
    values = np.where(finite, values / mean, 1.0)
    return values


def inverse_pdf_factors(values, labels, nbins: int, max_factor: float):
    import numpy as np

    values = np.asarray(values, dtype="float64")
    labels = np.asarray(labels, dtype="int32")
    factors = np.ones(len(values), dtype="float64")
    finite = np.isfinite(values)
    if finite.sum() < max(10, nbins):
        return factors

    lo, hi = np.nanpercentile(values[finite], [1.0, 99.0])
    if not math.isfinite(float(lo)) or not math.isfinite(float(hi)) or hi <= lo:
        return factors
    edges = np.linspace(lo, hi, nbins + 1)

    for cls in (0, 1):
        cls_mask = finite & (labels == cls) & (values >= lo) & (values <= hi)
        if cls_mask.sum() < nbins:
            continue
        counts, _ = np.histogram(values[cls_mask], bins=edges)
        counts = counts.astype("float64")
        nonzero = counts > 0.0
        if not nonzero.any():
            continue
        target = float(np.mean(counts[nonzero]))
        idx = np.clip(np.searchsorted(edges, values[cls_mask], side="right") - 1, 0, nbins - 1)
        local = np.ones(cls_mask.sum(), dtype="float64")
        valid = counts[idx] > 0.0
        local[valid] = target / counts[idx[valid]]
        local = np.clip(local, 1.0 / max_factor, max_factor)
        factors[cls_mask] *= normalize_mean_one(local)

    return factors


def compute_weights(frame, label_branch: str, args) -> tuple[object, dict]:
    import numpy as np

    labels = frame[label_branch].to_numpy(dtype="int32")
    weights = np.ones(len(frame), dtype="float64")
    diagnostics: dict[str, object] = {
        "use_event_weight": bool(args.use_event_weight),
        "et_reweight": bool(args.et_reweight),
        "eta_reweight": bool(args.eta_reweight),
    }

    if args.use_event_weight and args.weight_branch in frame.columns:
        ew = frame[args.weight_branch].to_numpy(dtype="float64")
        weights *= np.where(np.isfinite(ew) & (ew > 0.0), ew, 1.0)

    class_sums: dict[str, float] = {}
    for cls in (0, 1):
        mask = labels == cls
        class_sums[str(cls)] = float(weights[mask].sum())
    target = 0.5 * sum(class_sums.values())
    if target > 0.0:
        for cls in (0, 1):
            mask = labels == cls
            denom = float(weights[mask].sum())
            if denom > 0.0:
                weights[mask] *= target / denom

    if args.et_reweight and "cluster_Et" in frame.columns:
        weights *= inverse_pdf_factors(frame["cluster_Et"].to_numpy(), labels, args.flatten_bins, args.max_flatten_factor)
    if args.eta_reweight and "cluster_Eta" in frame.columns:
        weights *= inverse_pdf_factors(frame["cluster_Eta"].to_numpy(), labels, args.flatten_bins, args.max_flatten_factor)

    finite_positive = np.isfinite(weights) & (weights > 0.0)
    if finite_positive.any():
        median = float(np.median(weights[finite_positive]))
        cap = max(args.max_total_weight_factor * median, 1.0e-12)
        weights = np.where(finite_positive, np.clip(weights, 1.0e-12, cap), 1.0)
    else:
        weights = np.ones(len(frame), dtype="float64")

    diagnostics["sum_weight_class0"] = float(weights[labels == 0].sum())
    diagnostics["sum_weight_class1"] = float(weights[labels == 1].sum())
    diagnostics["min_weight"] = float(np.min(weights)) if len(weights) else math.nan
    diagnostics["max_weight"] = float(np.max(weights)) if len(weights) else math.nan
    diagnostics["mean_weight"] = float(np.mean(weights)) if len(weights) else math.nan
    return weights, diagnostics


def train_one(frame, features: list[str], label_branch: str, output: Path, metadata: dict, args) -> dict | None:
    import numpy as np
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import train_test_split
    from xgboost import XGBClassifier

    cols = features + [label_branch]
    frame = frame.loc[finite_mask(frame, cols)].copy()
    frame[label_branch] = frame[label_branch].astype(int)
    frame = frame[frame[label_branch].isin([0, 1])].copy()

    n_sig = int((frame[label_branch] == 1).sum())
    n_bkg = int((frame[label_branch] == 0).sum())
    enough = n_sig >= args.min_rows_per_class and n_bkg >= args.min_rows_per_class
    if not enough:
        msg = (
            f"Skipping {output.name}: class counts below minimum "
            f"(signal={n_sig}, background={n_bkg}, minimum={args.min_rows_per_class})"
        )
        if metadata.get("cent_range") == "all":
            raise SystemExit(msg)
        print(f"[WARN] {msg}")
        return None

    x = frame[features].to_numpy(dtype="float32")
    y = frame[label_branch].to_numpy(dtype="int32")
    weights, weight_report = compute_weights(frame, label_branch, args)

    stratify = y if min(n_sig, n_bkg) >= 2 else None
    x_train, x_test, y_train, y_test, w_train, w_test = train_test_split(
        x, y, weights, test_size=args.test_size, random_state=args.random_seed, stratify=stratify
    )

    model = XGBClassifier(
        n_estimators=args.n_estimators,
        max_depth=args.max_depth,
        learning_rate=args.learning_rate,
        subsample=args.subsample,
        colsample_bytree=args.colsample_bytree,
        objective="binary:logistic",
        eval_metric="auc",
        tree_method=args.tree_method,
        random_state=args.random_seed,
    )
    model.fit(x_train, y_train, sample_weight=w_train)

    pred = model.predict_proba(x_test)[:, 1]
    auc = float(roc_auc_score(y_test, pred, sample_weight=w_test)) if len(np.unique(y_test)) == 2 else math.nan

    output.parent.mkdir(parents=True, exist_ok=True)
    json_model = output.with_suffix(".xgb.json")
    model.get_booster().save_model(json_model)

    export_status = "not_attempted"
    export_error = ""
    try:
        import ROOT

        booster = model.get_booster()
        booster.feature_names = [f"f{i}" for i in range(len(features))]
        # ROOT 6.32's SaveXGBoost cannot parse XGBoost >=3 vector-valued
        # base_score strings like "[5E-1]". Patch only the Python config view
        # handed to TMVA; the saved XGBoost JSON remains untouched for audit.
        original_save_config = booster.save_config

        def tmva_save_config():
            import re

            text = original_save_config()
            return re.sub(r'"base_score":"\[([0-9eE+\-.]+)\]"', r'"base_score":"\1"', text)

        booster.save_config = tmva_save_config  # type: ignore[method-assign]
        ROOT.TMVA.Experimental.SaveXGBoost(model, "myBDT", str(output), num_inputs=len(features))
        export_status = "ok"
    except Exception as exc:  # noqa: BLE001
        export_status = f"failed: {exc}"
        export_error = str(exc)

    report = {
        **metadata,
        "output_tmva": str(output),
        "output_xgb_json": str(json_model),
        "features": features,
        "label_branch": label_branch,
        "n_rows": int(len(frame)),
        "n_signal": n_sig,
        "n_background": n_bkg,
        "auc": auc,
        "tmva_export": export_status,
        "weighting": weight_report,
        "xgboost": {
            "n_estimators": args.n_estimators,
            "max_depth": args.max_depth,
            "learning_rate": args.learning_rate,
            "subsample": args.subsample,
            "colsample_bytree": args.colsample_bytree,
            "tree_method": args.tree_method,
        },
    }
    output.with_suffix(".metadata.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    if export_error and not args.allow_tmva_export_failure:
        raise SystemExit(
            f"TMVA export failed for {output}: {export_error}. "
            "The XGBoost JSON was written for diagnostics, but this pipeline needs "
            "the TMVA/RBDT ROOT file for analysis consumption."
        )
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task", choices=["tight", "npb"], required=True)
    parser.add_argument("--input", nargs="+", type=Path, required=True, help="ROOT files or @manifest.list files")
    parser.add_argument("--tree", default="AuAuPhotonIDTrainingTree")
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--prefix", default="")
    parser.add_argument("--cent-bins", default="0:10,10:20,20:40,40:60,60:80,80:100")
    parser.add_argument("--label-branch", default=None)
    parser.add_argument(
        "--missing-label-value",
        type=int,
        choices=[0, 1],
        default=None,
        help="For files missing only the label branch, fill that label with this value. For PPG12-style NPB, use 1 for embedded physics-side sim.",
    )
    parser.add_argument("--features", default=None, help="Comma-separated override feature order")
    parser.add_argument("--weight-branch", default="event_weight")
    parser.add_argument("--use-event-weight", dest="use_event_weight", action="store_true", default=True)
    parser.add_argument("--no-event-weight", dest="use_event_weight", action="store_false")
    parser.add_argument("--et-reweight", dest="et_reweight", action="store_true", default=True)
    parser.add_argument("--no-et-reweight", dest="et_reweight", action="store_false")
    parser.add_argument("--eta-reweight", dest="eta_reweight", action="store_true", default=True)
    parser.add_argument("--no-eta-reweight", dest="eta_reweight", action="store_false")
    parser.add_argument("--flatten-bins", type=int, default=20)
    parser.add_argument("--max-flatten-factor", type=float, default=8.0)
    parser.add_argument("--max-total-weight-factor", type=float, default=50.0)
    parser.add_argument("--min-rows-per-class", type=int, default=10)
    parser.add_argument("--test-size", type=float, default=0.25)
    parser.add_argument("--random-seed", type=int, default=13)
    parser.add_argument("--n-estimators", type=int, default=450)
    parser.add_argument("--max-depth", type=int, default=4)
    parser.add_argument("--learning-rate", type=float, default=0.035)
    parser.add_argument("--subsample", type=float, default=0.85)
    parser.add_argument("--colsample-bytree", type=float, default=0.85)
    parser.add_argument("--tree-method", default="hist")
    parser.add_argument("--allow-tmva-export-failure", action="store_true",
                        help="Write diagnostics but do not fail if TMVA export fails. Not recommended for production pipeline tests.")
    args = parser.parse_args()

    features = [item.strip() for item in args.features.split(",")] if args.features else (
        PPG12_TIGHT_FEATURES if args.task == "tight" else PPG12_NPB_FEATURES
    )
    label_branch = args.label_branch or ("is_signal" if args.task == "tight" else "npb_label")
    if args.task == "npb" and label_branch == "is_signal":
        raise SystemExit("NPB training needs a real NPB label branch, not is_signal.")
    if args.task == "npb" and label_branch == "is_npb":
        raise SystemExit(
            "Use npb_label for PPG12-style NPB training: 1=physics-like, 0=NPB. "
            "The is_npb branch is kept only as an audit inverse."
        )

    paths = expand_input_paths(args.input)
    required_columns = sorted(set(features + [label_branch, "centrality"]))
    optional_columns = [args.weight_branch]
    frame, optional_seen = load_frame(
        paths,
        args.tree,
        required_columns,
        optional_columns,
        label_branch,
        args.missing_label_value,
    )

    args.outdir.mkdir(parents=True, exist_ok=True)
    prefix = args.prefix or f"auau_{args.task}_bdt"
    reports = []
    common_metadata = {
        "task": args.task,
        "input_files": [str(path) for path in paths],
        "tree": args.tree,
        "optional_branches_seen": optional_seen,
        "python": sys.version,
    }

    all_report = train_one(
        frame,
        features,
        label_branch,
        args.outdir / f"{prefix}_allCent_tmva.root",
        {**common_metadata, "cent_range": "all"},
        args,
    )
    if all_report is not None:
        reports.append(all_report)

    for lo, hi in parse_cent_bins(args.cent_bins):
        sub = frame[(frame["centrality"] >= lo) & (frame["centrality"] < hi)]
        if len(sub) == 0:
            print(f"[WARN] Skipping {cent_tag(lo, hi)}: no rows")
            continue
        report = train_one(
            sub,
            features,
            label_branch,
            args.outdir / f"{prefix}_{cent_tag(lo, hi)}_tmva.root",
            {**common_metadata, "cent_range": [lo, hi]},
            args,
        )
        if report is not None:
            reports.append(report)

    summary = args.outdir / f"{prefix}_summary.json"
    summary.write_text(json.dumps(reports, indent=2, sort_keys=True) + "\n")
    print(f"[OK] wrote {summary}")
    for report in reports:
        print(
            f"{Path(report['output_tmva']).name}: rows={report['n_rows']} "
            f"sig={report['n_signal']} bkg={report['n_background']} "
            f"auc={report['auc']:.4f} tmva={report['tmva_export']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
