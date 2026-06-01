#!/usr/bin/env python3
"""Train compact AuAu photon-ID logistic-regression artifacts.

This is a direct photon-feature baseline, not a BDT/MLP stacker.  It consumes
the same AuAuPhotonIDTrainingTree extraction products as the BDT/MLP trainers
and exports deterministic JSON artifacts for C++ runtime evaluation.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from train_auau_photon_mlp import (  # noqa: E402
    BASE_AND_3X3_FEATURES,
    EXTENDED_SHOWER_FEATURES,
    WIDTH_RATIO_FEATURES,
    add_derived_features,
    auc_score,
    calibration_report,
    compute_weights,
    expand_input_paths,
    expand_required_columns,
    json_ready,
    limit_paths_by_group,
    load_frame,
    parse_bin_spec,
    parse_range,
    select_training_rows,
    sigmoid,
    split_masks,
    threshold_for_signal_efficiency,
)


SCHEMA = "RJ_AUAU_TIGHT_LOGREG_V1"
REGISTRY_SCHEMA = "RJ_AUAU_TIGHT_LOGREG_REGISTRY_V1"
DEFAULT_PT_BINS = "15,17,19,21,23,25,27,30,35"
DEFAULT_CENT3 = "0,20,50,80"
DEFAULT_CENT7 = "0,10,20,30,40,50,60,80"

FULL_FEATURES = []
for _feature in BASE_AND_3X3_FEATURES + EXTENDED_SHOWER_FEATURES + WIDTH_RATIO_FEATURES + ["centrality"]:
    if _feature not in FULL_FEATURES:
        FULL_FEATURES.append(_feature)

PRODUCT_SPECS = {
    "centEtFullLogReg_pt1535": {
        "variant": "auauTightLogReg",
        "label": "Global full-feature logistic regression with E_T and centrality",
        "filename": "auau_tight_logreg_centEtFull_pt1535.json",
        "features": FULL_FEATURES,
        "pt_range": (15.0, 35.0),
        "routes": [],
    },
    "etFineCent7FullLogReg_pt1535": {
        "variant": "auauTightLogReg",
        "label": "Fine E_T x 7-centrality full-feature routed logistic regression",
        "filename": "auau_tight_logreg_etFineCent7Full_pt1535.json",
        "features": FULL_FEATURES,
        "pt_range": (15.0, 35.0),
        "pt_edges": [15.0, 17.0, 19.0, 21.0, 23.0, 25.0, 27.0, 30.0, 35.0],
        "cent_edges": [0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0],
    },
    "etFineCent3FullLogReg_pt1535": {
        "variant": "auauTightLogReg",
        "label": "Fine E_T x 3-centrality full-feature routed logistic regression",
        "filename": "auau_tight_logreg_etFineCent3Full_pt1535.json",
        "features": FULL_FEATURES,
        "pt_range": (15.0, 35.0),
        "pt_edges": [15.0, 17.0, 19.0, 21.0, 23.0, 25.0, 27.0, 30.0, 35.0],
        "cent_edges": [0.0, 20.0, 50.0, 80.0],
    },
}


@dataclass
class LinearModel:
    feature_names: list[str]
    impute: np.ndarray
    mean: np.ndarray
    scale: np.ndarray
    coef: np.ndarray
    intercept: float


def write_json(path: Path, payload) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(json_ready(payload), indent=2, sort_keys=True) + "\n")


def load_artifact(path: Path | str) -> dict:
    return json.loads(Path(path).read_text())


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def parse_products(text: str) -> list[str]:
    key = text.strip().lower()
    if key in ("all", "default"):
        return list(PRODUCT_SPECS)
    if key in ("primary", "global"):
        return ["centEtFullLogReg_pt1535"]
    if key in ("routed", "etfine-cent7", "cent7"):
        return ["etFineCent7FullLogReg_pt1535"]
    products = [item.strip() for item in text.split(",") if item.strip()]
    unknown = [item for item in products if item not in PRODUCT_SPECS]
    if unknown:
        raise SystemExit("Unknown logreg product(s): " + ", ".join(unknown))
    return products


def stable_route_seed(base: int, product: str, route: str) -> int:
    import hashlib

    raw = f"{base}|{product}|{route}".encode("utf-8")
    return int(hashlib.sha256(raw).hexdigest()[:8], 16)


def finite_design(frame, features: list[str]) -> tuple[np.ndarray, np.ndarray]:
    x = frame[features].to_numpy(dtype="float64")
    finite = np.isfinite(x)
    return x, finite.all(axis=1)


def standardize_with_impute(x: np.ndarray, train_mask: np.ndarray, clip: float) -> tuple[np.ndarray, LinearModel]:
    x_train = x[train_mask]
    impute = np.nanmedian(np.where(np.isfinite(x_train), x_train, np.nan), axis=0)
    impute = np.where(np.isfinite(impute), impute, 0.0)
    clean = np.where(np.isfinite(x), x, impute)
    mean = np.mean(clean[train_mask], axis=0)
    scale = np.std(clean[train_mask], axis=0)
    scale = np.where(np.isfinite(scale) & (scale > 1.0e-9), scale, 1.0)
    mean = np.where(np.isfinite(mean), mean, 0.0)
    z = (clean - mean) / scale
    if clip > 0.0:
        z = np.clip(z, -clip, clip)
    dummy = LinearModel([], impute.astype("float64"), mean.astype("float64"), scale.astype("float64"), np.zeros(x.shape[1]), 0.0)
    return z.astype("float64"), dummy


def weighted_bce(y: np.ndarray, p: np.ndarray, w: np.ndarray, coef: np.ndarray, l2: float) -> float:
    p = np.clip(p, 1.0e-7, 1.0 - 1.0e-7)
    denom = max(float(np.sum(w)), 1.0)
    loss = -float(np.sum(w * (y * np.log(p) + (1.0 - y) * np.log(1.0 - p))) / denom)
    return loss + 0.5 * float(l2) * float(np.sum(coef * coef))


def fit_logreg(
    name: str,
    features: list[str],
    selected,
    weights: np.ndarray,
    train_mask: np.ndarray,
    val_mask: np.ndarray,
    args,
    seed: int,
) -> tuple[LinearModel, list[dict]]:
    x_raw, finite_rows = finite_design(selected, features)
    if not finite_rows.any():
        raise SystemExit(f"{name}: no finite rows for requested features")
    labels = selected[args.label_branch].to_numpy(dtype="int32")
    z, base = standardize_with_impute(x_raw, train_mask, args.input_clip)
    y_train = labels[train_mask].astype("float64")
    w_train = weights[train_mask].astype("float64")
    if not (np.any(y_train == 0) and np.any(y_train == 1)):
        raise SystemExit(f"{name}: train split is missing a class")

    rng = np.random.default_rng(seed)
    coef = rng.normal(0.0, 0.01, size=z.shape[1])
    intercept = 0.0
    m_coef = np.zeros_like(coef)
    v_coef = np.zeros_like(coef)
    m_int = 0.0
    v_int = 0.0
    beta1 = 0.9
    beta2 = 0.999
    step = 0
    best_loss = math.inf
    best_coef = coef.copy()
    best_intercept = intercept
    stale = 0
    history: list[dict] = []
    train_idx = np.flatnonzero(train_mask)
    val_idx = np.flatnonzero(val_mask)
    batch_size = min(max(1024, int(args.batch_size)), max(1, len(train_idx)))
    steps_per_epoch = max(1, int(math.ceil(len(train_idx) / batch_size)))
    for epoch in range(1, args.epochs + 1):
        rng.shuffle(train_idx)
        for start in range(0, len(train_idx), batch_size):
            step += 1
            idx = train_idx[start : start + batch_size]
            zz = z[idx]
            yy = labels[idx].astype("float64")
            ww = weights[idx].astype("float64")
            prob = sigmoid(zz @ coef + intercept)
            err = (prob - yy) * ww / max(float(np.mean(ww)), 1.0e-12)
            grad_coef = (zz.T @ err) / max(1, len(idx)) + args.l2 * coef
            grad_int = float(np.mean(err))
            m_coef = beta1 * m_coef + (1.0 - beta1) * grad_coef
            v_coef = beta2 * v_coef + (1.0 - beta2) * grad_coef * grad_coef
            m_int = beta1 * m_int + (1.0 - beta1) * grad_int
            v_int = beta2 * v_int + (1.0 - beta2) * grad_int * grad_int
            coef -= args.learning_rate * (m_coef / (1.0 - beta1**step)) / (np.sqrt(v_coef / (1.0 - beta2**step)) + 1.0e-8)
            intercept -= args.learning_rate * (m_int / (1.0 - beta1**step)) / (math.sqrt(v_int / (1.0 - beta2**step)) + 1.0e-8)

        if epoch == 1 or epoch == args.epochs or epoch % args.progress_every == 0:
            train_prob = sigmoid(z[train_idx] @ coef + intercept)
            row = {
                "model": name,
                "epoch": int(epoch),
                "step": int(step),
                "train_loss": weighted_bce(labels[train_idx], train_prob, weights[train_idx], coef, args.l2),
                "train_auc": auc_score(labels[train_idx], train_prob),
            }
            if len(val_idx):
                val_prob = sigmoid(z[val_idx] @ coef + intercept)
                row["val_loss"] = weighted_bce(labels[val_idx], val_prob, weights[val_idx], coef, args.l2)
                row["val_auc"] = auc_score(labels[val_idx], val_prob)
                monitor = float(row["val_loss"])
            else:
                row["val_loss"] = math.nan
                row["val_auc"] = math.nan
                monitor = float(row["train_loss"])
            history.append(row)
            print(
                f"[trainAuAuPhotonLogReg] model={name} epoch={epoch} step={step}/{args.epochs * steps_per_epoch} "
                f"train_loss={row['train_loss']:.6g} val_loss={row['val_loss']:.6g} "
                f"train_auc={row['train_auc']:.5f} val_auc={row['val_auc']:.5f}",
                flush=True,
            )
            if monitor + args.min_delta < best_loss:
                best_loss = monitor
                best_coef = coef.copy()
                best_intercept = float(intercept)
                stale = 0
            else:
                stale += args.progress_every
                if stale >= args.patience:
                    print(f"[trainAuAuPhotonLogReg] model={name} early_stop epoch={epoch}", flush=True)
                    break

    return LinearModel(features, base.impute, base.mean, base.scale, best_coef.astype("float64"), best_intercept), history


def score_model(model: LinearModel, frame) -> np.ndarray:
    x = frame[model.feature_names].to_numpy(dtype="float64")
    clean = np.where(np.isfinite(x), x, model.impute)
    z = np.nan_to_num((clean - model.mean) / model.scale, nan=0.0, posinf=0.0, neginf=0.0)
    return sigmoid(z @ model.coef + model.intercept)


def model_payload(model: LinearModel, name: str, metadata: dict) -> dict:
    return {
        "kind": "logistic",
        "name": name,
        "feature_names": model.feature_names,
        "impute": model.impute,
        "mean": model.mean,
        "scale": model.scale,
        "coef": model.coef,
        "intercept": float(model.intercept),
        "metadata": metadata,
    }


def artifact_score(artifact: dict, frame) -> np.ndarray:
    if artifact.get("routes"):
        out = np.full(len(frame), np.nan, dtype="float64")
        et = frame["cluster_Et"].to_numpy(dtype="float64")
        cent = frame["centrality"].to_numpy(dtype="float64")
        for route in artifact["routes"]:
            mask = np.isfinite(et) & (et >= float(route["pt_lo"])) & (et < float(route["pt_hi"]))
            if route.get("cent_lo") is not None and route.get("cent_hi") is not None:
                mask &= np.isfinite(cent) & (cent >= float(route["cent_lo"])) & (cent < float(route["cent_hi"]))
            if not mask.any():
                continue
            m = payload_to_model(route["model"])
            out[mask] = score_model(m, frame.loc[mask])
        return out
    return score_model(payload_to_model(artifact["model"]), frame)


def payload_to_model(payload: dict) -> LinearModel:
    return LinearModel(
        list(payload["feature_names"]),
        np.asarray(payload["impute"], dtype="float64"),
        np.asarray(payload["mean"], dtype="float64"),
        np.asarray(payload["scale"], dtype="float64"),
        np.asarray(payload["coef"], dtype="float64"),
        float(payload["intercept"]),
    )


def summarize_scores(product: str, frame, score: np.ndarray, target: float) -> dict:
    y = frame["is_signal"].to_numpy(dtype="int32")
    finite = np.isfinite(score) & np.isin(y, [0, 1])
    wp = threshold_for_signal_efficiency(y[finite], score[finite], target)
    return {
        "product": product,
        "entries": int(len(frame)),
        "scored_entries": int(finite.sum()),
        "finite_fraction": float(finite.mean()) if len(finite) else math.nan,
        "auc": auc_score(y[finite], score[finite]) if finite.any() else math.nan,
        "wp80_fake": None if wp is None else float(wp["background_fake_rate"]),
        "wp80_threshold": None if wp is None else float(wp["threshold"]),
        "calibration": calibration_report(y[finite], score[finite]) if finite.any() else {},
    }


def train_product(product: str, frame, args, outdir: Path) -> dict:
    spec = PRODUCT_SPECS[product]
    features = list(spec["features"])
    selected = select_training_rows(frame, features, args)
    labels = selected[args.label_branch].to_numpy(dtype="int32")
    counts = {str(cls): int((labels == cls).sum()) for cls in (0, 1)}
    if min(counts.values()) < args.min_rows_per_class:
        raise SystemExit(f"{product}: not enough rows per class after cuts: {counts}")
    weights, weight_report = compute_weights(selected, args.label_branch, args)
    train_mask, val_mask, test_mask = split_masks(selected, args.label_branch, args)
    print(
        f"[trainAuAuPhotonLogReg] product={product} selected_rows={len(selected)} "
        f"class_counts={counts} train={int(train_mask.sum())} val={int(val_mask.sum())} test={int(test_mask.sum())}",
        flush=True,
    )

    routes: list[dict] = []
    histories: list[dict] = []
    if "pt_edges" in spec and "cent_edges" in spec:
        pt_edges = list(spec["pt_edges"])
        cent_edges = list(spec["cent_edges"])
        for ip, (pt_lo, pt_hi) in enumerate(zip(pt_edges[:-1], pt_edges[1:])):
            for ic, (cent_lo, cent_hi) in enumerate(zip(cent_edges[:-1], cent_edges[1:])):
                route_label = f"pt{pt_lo:g}to{pt_hi:g}_cent{cent_lo:g}to{cent_hi:g}".replace(".", "p")
                mask = (
                    selected["cluster_Et"].to_numpy(dtype="float64") >= pt_lo
                ) & (
                    selected["cluster_Et"].to_numpy(dtype="float64") < pt_hi
                ) & (
                    selected["centrality"].to_numpy(dtype="float64") >= cent_lo
                ) & (
                    selected["centrality"].to_numpy(dtype="float64") < cent_hi
                )
                route_labels = labels[mask]
                if int((route_labels == 0).sum()) < args.min_rows_per_route_class or int((route_labels == 1).sum()) < args.min_rows_per_route_class:
                    print(f"[trainAuAuPhotonLogReg][WARN] skip route {route_label}: class counts too small", flush=True)
                    continue
                local = selected.loc[mask].copy()
                local_weights = weights[mask]
                local_train, local_val, local_test = split_masks(local, args.label_branch, args)
                model, history = fit_logreg(
                    f"{product}:{route_label}",
                    features,
                    local,
                    local_weights,
                    local_train,
                    local_val,
                    args,
                    stable_route_seed(args.random_seed, product, route_label),
                )
                histories.extend(history)
                routes.append(
                    {
                        "label": route_label,
                        "pt_lo": float(pt_lo),
                        "pt_hi": float(pt_hi),
                        "cent_lo": float(cent_lo),
                        "cent_hi": float(cent_hi),
                        "model": model_payload(model, f"{product}:{route_label}", {"route": route_label}),
                    }
                )
    else:
        model, history = fit_logreg(product, features, selected, weights, train_mask, val_mask, args, args.random_seed)
        histories.extend(history)
        routes = []

    artifact = {
        "schema": SCHEMA,
        "name": product,
        "variant": spec["variant"],
        "label": spec["label"],
        "features": features,
        "pt_range": [float(spec["pt_range"][0]), float(spec["pt_range"][1])],
        "score_definition": "sigmoid(standardized_features dot coef + intercept)",
        "metadata": {
            "abcd_safe": True,
            "no_isolation_inputs": True,
            "no_bdt_or_mlp_score_inputs": True,
            "weight_report": weight_report,
            "training_args": {
                "epochs": int(args.epochs),
                "patience": int(args.patience),
                "learning_rate": float(args.learning_rate),
                "l2": float(args.l2),
                "batch_size": int(args.batch_size),
                "pt_range": args.pt_range,
            },
        },
        "routes": routes,
    }
    if not routes:
        artifact["model"] = model_payload(model, product, {"global": True})
    if "pt_edges" in spec:
        artifact["pt_edges"] = spec["pt_edges"]
        artifact["cent_edges"] = spec["cent_edges"]

    out_json = outdir / spec["filename"]
    write_json(out_json, artifact)
    score = artifact_score(artifact, selected)
    report = summarize_scores(product, selected, score, args.selection_target_signal_efficiency)
    report["class_counts"] = counts
    report["output_json"] = str(out_json)
    report["routes"] = len(routes)
    write_csv(outdir / f"{product}.training_history.csv", histories)
    write_json(outdir / f"{product}.training_report.json", report)
    print(
        f"[trainAuAuPhotonLogReg] product={product} auc={report['auc']:.5f} "
        f"wp80_fake={report['wp80_fake']} routes={len(routes)} output={out_json}",
        flush=True,
    )
    return report


def load_training_frame(args, products: list[str]):
    required_columns = {args.label_branch, "cluster_Et", "cluster_Eta", "centrality", "vertexz", "run", "evt"}
    for product in products:
        required_columns.update(expand_required_columns(PRODUCT_SPECS[product]["features"]))
    optional_columns = [args.weight_branch]
    all_paths = expand_input_paths(args.input, args.source, args.check_input_files)
    paths, manifest_report = limit_paths_by_group(all_paths, args.max_files_per_sample, args.random_seed)
    frame, seen_optional, read_report = load_frame(paths, args.tree, sorted(required_columns), optional_columns)
    frame = add_derived_features(frame)
    write_json(
        args.outdir / "training_manifest_summary.json",
        {
            "schema": "RJ_AUAU_TIGHT_LOGREG_TRAINING_MANIFEST_SUMMARY_V1",
            "source": str(args.source) if args.source else None,
            "products": products,
            "manifest_report": manifest_report,
            "seen_optional": seen_optional,
            "read_report": read_report,
        },
    )
    return frame


def self_test() -> int:
    class A:
        pass

    args = A()
    args.label_branch = "is_signal"
    args.input_clip = 8.0
    args.batch_size = 128
    args.epochs = 6
    args.patience = 8
    args.progress_every = 2
    args.learning_rate = 0.04
    args.l2 = 1.0e-4
    args.min_delta = 0.0
    rng = np.random.default_rng(123)
    import pandas as pd

    n = 2000
    x1 = rng.normal(size=n)
    x2 = rng.normal(size=n)
    p = sigmoid(1.8 * x1 - 0.6 * x2 + 0.2)
    y = rng.binomial(1, p)
    frame = pd.DataFrame({"is_signal": y, "x1": x1, "x2": x2})
    weights = np.ones(n)
    train = np.zeros(n, dtype=bool)
    train[:1400] = True
    val = np.zeros(n, dtype=bool)
    val[1400:1700] = True
    model, _ = fit_logreg("selftest", ["x1", "x2"], frame, weights, train, val, args, 123)
    score = score_model(model, frame)
    auc = auc_score(y, score)
    if auc < 0.75:
        raise SystemExit(f"self-test AUC too low: {auc}")
    artifact = {"schema": SCHEMA, "name": "selftest", "model": model_payload(model, "selftest", {})}
    roundtrip = score_model(payload_to_model(artifact["model"]), frame)
    if float(np.max(np.abs(score - roundtrip))) > 1.0e-12:
        raise SystemExit("artifact roundtrip mismatch")
    print(f"[trainAuAuPhotonLogReg] self-test passed auc={auc:.5f}")
    return 0


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--input", type=Path, nargs="*", default=[])
    ap.add_argument("--source", type=Path, default=None)
    ap.add_argument("--check-input-files", action="store_true")
    ap.add_argument("--tree", default="AuAuPhotonIDTrainingTree")
    ap.add_argument("--outdir", type=Path, required=False)
    ap.add_argument("--products", default="all")
    ap.add_argument("--pt-range", default="15:35")
    ap.add_argument("--centrality-range", default="")
    ap.add_argument("--train-pt-bins", default=DEFAULT_PT_BINS)
    ap.add_argument("--label-branch", default="is_signal")
    ap.add_argument("--weight-branch", default="event_weight")
    ap.add_argument("--use-event-weight", dest="use_event_weight", action="store_true", default=True)
    ap.add_argument("--no-event-weight", dest="use_event_weight", action="store_false")
    ap.add_argument("--et-reweight", dest="et_reweight", action="store_true", default=True)
    ap.add_argument("--no-et-reweight", dest="et_reweight", action="store_false")
    ap.add_argument("--eta-reweight", dest="eta_reweight", action="store_true", default=True)
    ap.add_argument("--no-eta-reweight", dest="eta_reweight", action="store_false")
    ap.add_argument("--flatten-bins", type=int, default=20)
    ap.add_argument("--max-flatten-factor", type=float, default=8.0)
    ap.add_argument("--max-total-weight-factor", type=float, default=50.0)
    ap.add_argument("--pt-bin-weight-mode", choices=["none", "equal", "highpt"], default="equal")
    ap.add_argument("--pt-bin-weight-spec", default="")
    ap.add_argument("--hard-example-branch", default="")
    ap.add_argument("--hard-background-factor", type=float, default=0.0)
    ap.add_argument("--hard-signal-factor", type=float, default=0.0)
    ap.add_argument("--hard-example-power", type=float, default=2.0)
    ap.add_argument("--max-files-per-sample", type=int, default=int(os.environ.get("RJ_AUAU_LOGREG_TRAIN_MAX_FILES_PER_SAMPLE", "0")))
    ap.add_argument("--max-rows", type=int, default=int(os.environ.get("RJ_AUAU_LOGREG_TRAIN_MAX_ROWS", "0")))
    ap.add_argument("--max-rows-per-class", type=int, default=int(os.environ.get("RJ_AUAU_LOGREG_TRAIN_MAX_ROWS_PER_CLASS", "0")))
    ap.add_argument("--max-rows-per-pt-bin-class", type=int, default=int(os.environ.get("RJ_AUAU_LOGREG_TRAIN_MAX_ROWS_PER_PT_BIN_CLASS", "160000")))
    ap.add_argument("--min-rows-per-class", type=int, default=100)
    ap.add_argument("--min-rows-per-route-class", type=int, default=100)
    ap.add_argument("--validation-fraction", type=float, default=0.20)
    ap.add_argument("--test-fraction", type=float, default=0.10)
    ap.add_argument("--random-seed", type=int, default=98713)
    ap.add_argument("--epochs", type=int, default=int(os.environ.get("RJ_AUAU_LOGREG_TRAIN_EPOCHS", "60")))
    ap.add_argument("--patience", type=int, default=int(os.environ.get("RJ_AUAU_LOGREG_TRAIN_PATIENCE", "14")))
    ap.add_argument("--progress-every", type=int, default=int(os.environ.get("RJ_AUAU_LOGREG_TRAIN_PROGRESS_EVERY", "2")))
    ap.add_argument("--batch-size", type=int, default=int(os.environ.get("RJ_AUAU_LOGREG_TRAIN_BATCH_SIZE", "262144")))
    ap.add_argument("--learning-rate", type=float, default=float(os.environ.get("RJ_AUAU_LOGREG_TRAIN_LEARNING_RATE", "0.018")))
    ap.add_argument("--l2", type=float, default=float(os.environ.get("RJ_AUAU_LOGREG_TRAIN_L2", "2.0e-3")))
    ap.add_argument("--min-delta", type=float, default=1.0e-5)
    ap.add_argument("--input-clip", type=float, default=8.0)
    ap.add_argument("--selection-target-signal-efficiency", type=float, default=0.80)
    return ap.parse_args()


def main() -> int:
    args = parse_args()
    if args.self_test:
        return self_test()
    if args.outdir is None:
        raise SystemExit("--outdir is required unless --self-test is used")
    args.outdir.mkdir(parents=True, exist_ok=True)
    parse_range(args.pt_range)
    if args.centrality_range:
        parse_range(args.centrality_range)
    if args.train_pt_bins:
        parse_bin_spec(args.train_pt_bins)
    products = parse_products(args.products)
    frame = load_training_frame(args, products)
    reports = []
    for product in products:
        reports.append(train_product(product, frame, args, args.outdir))
    registry = {
        "schema": REGISTRY_SCHEMA,
        "models": reports,
        "product_specs": {product: PRODUCT_SPECS[product] for product in products},
    }
    write_json(args.outdir / "model_registry.json", registry)
    write_csv(args.outdir / "logreg_training_summary.csv", reports)
    print(f"[trainAuAuPhotonLogReg] wrote registry: {args.outdir / 'model_registry.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
