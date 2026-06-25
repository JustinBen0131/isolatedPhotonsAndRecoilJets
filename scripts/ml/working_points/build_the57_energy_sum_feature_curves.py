#!/usr/bin/env python3
"""Build slide-12 energy-sum feature curves from the full THE-57 matrix."""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

FEATURES = [
    "cluster_Et",
    "cluster_weta_cogx",
    "cluster_wphi_cogx",
    "cluster_weta33_cogx",
    "cluster_wphi33_cogx",
    "vertexz",
    "cluster_Eta",
    "e11_over_e33",
    "cluster_et1",
    "cluster_et2",
    "cluster_et3",
    "cluster_et4",
    "e32_over_e35",
    "centrality",
]

VARIABLES = (
    {
        "key": "et1",
        "feature": "cluster_et1",
        "variable": "cluster_et1",
        "xlim": (0.25, 1.02),
        "edges": np.linspace(0.0, 1.02, 103),
    },
    {
        "key": "e11e33",
        "feature": "e11_over_e33",
        "variable": "E11/E33",
        "xlim": (0.00, 0.96),
        "edges": np.linspace(0.0, 1.00, 101),
    },
    {
        "key": "e32e35",
        "feature": "e32_over_e35",
        "variable": "E32/E35",
        "xlim": (0.45, 1.02),
        "edges": np.linspace(0.0, 1.02, 103),
    },
)

PT_GROUPS = (
    (r"Low $p_T$: 15-22 GeV", ((15.0, 18.0), (18.0, 20.0), (20.0, 22.0))),
    (r"Mid $p_T$: 22-28 GeV", ((22.0, 24.0), (24.0, 26.0), (26.0, 28.0))),
    (r"High $p_T$: 28-35 GeV", ((28.0, 30.0), (30.0, 35.0))),
)

CENT_FOCUS = (0.0, 20.0)
PPG12_EXACT_ETA_RANGE = (-0.7, 0.7)
PPG12_EXACT_N_BINS = 20
PPG12_EXACT_ET_WEIGHT_CAP = 800.0


@dataclass(frozen=True)
class Curve:
    sample: str
    variable: str
    feature: str
    pt_group: str
    centrality: str
    stage: str
    edges: list[float]
    density: list[float]
    visible_integral: float
    total_integral: float
    mean: float
    rms: float


def normalize_mean_one(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype="float64")
    finite = np.isfinite(values)
    if not finite.any():
        return values
    mean = float(np.mean(values[finite]))
    if mean <= 0.0 or not math.isfinite(mean):
        return values
    out = values.copy()
    out[finite] = out[finite] / mean
    return out


def inverse_pdf_weights(values: np.ndarray, *, n_bins: int, fixed_range=None, weight_cap: float | None = None):
    from scipy.interpolate import UnivariateSpline

    values = np.asarray(values, dtype="float64")
    weights = np.ones(len(values), dtype="float64")
    finite = np.isfinite(values)
    report = {
        "n_entries": int(len(values)),
        "n_finite": int(finite.sum()),
        "n_bins": int(n_bins),
        "fixed_range": list(fixed_range) if fixed_range is not None else None,
        "weight_cap": float(weight_cap) if weight_cap is not None else None,
        "status": "ok",
    }
    if finite.sum() < max(10, n_bins):
        report["status"] = "insufficient_finite_values"
        return weights, report
    vals = values[finite]
    lo, hi = (float(np.min(vals)), float(np.max(vals))) if fixed_range is None else (float(fixed_range[0]), float(fixed_range[1]))
    report["range"] = [lo, hi]
    if not math.isfinite(lo) or not math.isfinite(hi) or hi <= lo:
        report["status"] = "invalid_range"
        return weights, report
    hist, bin_edges = np.histogram(vals, bins=np.linspace(lo, hi, n_bins + 1), density=True)
    hist = hist.astype("float64") * float(n_bins)
    centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    good_hist = np.isfinite(hist)
    if good_hist.sum() < 4:
        report["status"] = "insufficient_histogram_support"
        return weights, report
    spline = UnivariateSpline(centers[good_hist], hist[good_hist], s=0.0)
    pdf = spline(vals)
    report["pdf_min_before_clip"] = float(np.nanmin(pdf)) if len(pdf) else math.nan
    report["pdf_max_before_clip"] = float(np.nanmax(pdf)) if len(pdf) else math.nan
    report["pdf_negative_fraction_before_clip"] = float(np.mean(pdf < 0.0)) if len(pdf) else math.nan
    pdf = np.clip(pdf, a_min=1.0e-3, a_max=None)
    local = 1.0 / pdf
    if weight_cap is not None:
        local = np.clip(local, a_min=None, a_max=float(weight_cap))
    local = normalize_mean_one(local)
    weights[finite] = local
    report["min_weight"] = float(np.min(local)) if len(local) else math.nan
    report["max_weight"] = float(np.max(local)) if len(local) else math.nan
    report["mean_weight"] = float(np.mean(local)) if len(local) else math.nan
    return weights, report


def compute_ppg12_exact_weights(labels: np.ndarray, et: np.ndarray, eta: np.ndarray):
    labels = np.asarray(labels, dtype="int32")
    weights = np.ones(len(labels), dtype="float64")
    report: dict[str, object] = {
        "weight_mode": "ppg12-exact",
        "event_weight_used": False,
        "vertex_reweight": False,
        "centrality_event_weight": False,
        "cross_section_weight_used_for_training": False,
        "weights_computed_before_binning": True,
        "eta_range": list(PPG12_EXACT_ETA_RANGE),
        "eta_bins": PPG12_EXACT_N_BINS,
        "et_bins": PPG12_EXACT_N_BINS,
        "et_weight_cap": PPG12_EXACT_ET_WEIGHT_CAP,
    }
    class_counts = {str(cls): int((labels == cls).sum()) for cls in (0, 1)}
    if class_counts["0"] <= 0 or class_counts["1"] <= 0:
        raise SystemExit(f"PPG12-exact weights need both classes; observed counts={class_counts}")
    n_total = class_counts["0"] + class_counts["1"]
    class_factors = {}
    eta_reports = {}
    et_reports = {}
    for cls in (0, 1):
        mask = labels == cls
        factor = float(n_total) / (2.0 * float(class_counts[str(cls)]))
        class_factors[str(cls)] = factor
        weights[mask] *= factor
        eta_w, eta_report = inverse_pdf_weights(
            eta[mask],
            n_bins=PPG12_EXACT_N_BINS,
            fixed_range=PPG12_EXACT_ETA_RANGE,
            weight_cap=None,
        )
        weights[mask] *= eta_w
        eta_reports[str(cls)] = eta_report
        et_w, et_report = inverse_pdf_weights(
            et[mask],
            n_bins=PPG12_EXACT_N_BINS,
            fixed_range=None,
            weight_cap=PPG12_EXACT_ET_WEIGHT_CAP,
        )
        weights[mask] *= et_w
        et_reports[str(cls)] = et_report
    finite_positive = np.isfinite(weights) & (weights > 0.0)
    if not finite_positive.all():
        raise SystemExit(f"PPG12-exact weights produced {int((~finite_positive).sum())} bad rows")
    report["class_counts"] = class_counts
    report["class_weight_factors"] = class_factors
    report["eta_reweight"] = eta_reports
    report["et_reweight"] = et_reports
    report["sum_weight_class0"] = float(weights[labels == 0].sum())
    report["sum_weight_class1"] = float(weights[labels == 1].sum())
    report["min_weight"] = float(np.min(weights)) if len(weights) else math.nan
    report["max_weight"] = float(np.max(weights)) if len(weights) else math.nan
    report["mean_weight"] = float(np.mean(weights)) if len(weights) else math.nan
    return weights, report


def score_matrix(z, model: Path, chunk_size: int) -> np.ndarray:
    from xgboost import XGBClassifier

    clf = XGBClassifier()
    clf.load_model(str(model))
    n = len(z["is_signal"])
    score = np.empty(n, dtype=np.float32)
    for start in range(0, n, chunk_size):
        stop = min(n, start + chunk_size)
        x = np.column_stack([np.asarray(z[name][start:stop], dtype=np.float32) for name in FEATURES])
        score[start:stop] = clf.predict_proba(x)[:, 1].astype(np.float32, copy=False)
        print(f"[score] {stop}/{n}", flush=True)
    return score


def crop_and_density(edges: np.ndarray, counts: np.ndarray, xlim: tuple[float, float]):
    centers = 0.5 * (edges[:-1] + edges[1:])
    keep = (centers >= xlim[0]) & (centers <= xlim[1])
    if keep.any():
        first = int(np.where(keep)[0][0])
        last = int(np.where(keep)[0][-1])
        edges = edges[first : last + 2]
        counts = counts[first : last + 1]
    visible = float(np.sum(counts))
    widths = np.diff(edges)
    density = counts / (visible * widths) if visible > 0.0 else np.zeros_like(counts)
    return edges, density, visible


def make_curve(sample: str, variable: dict, pt_label: str, stage: str, values, weights) -> Curve:
    values = np.asarray(values, dtype="float64")
    weights = np.asarray(weights, dtype="float64")
    good = np.isfinite(values) & np.isfinite(weights) & (weights > 0.0)
    values = values[good]
    weights = weights[good]
    edges = np.asarray(variable["edges"], dtype="float64")
    counts, edges = np.histogram(values, bins=edges, weights=weights)
    total = float(np.sum(counts))
    mean = float(np.sum(values * weights) / np.sum(weights)) if weights.size and np.sum(weights) > 0.0 else math.nan
    rms = (
        float(np.sqrt(np.sum(weights * (values - mean) ** 2) / np.sum(weights)))
        if weights.size and np.sum(weights) > 0.0 and math.isfinite(mean)
        else math.nan
    )
    shown_edges, density, visible = crop_and_density(edges, counts, tuple(variable["xlim"]))
    return Curve(
        sample=sample,
        variable=str(variable["variable"]),
        feature=str(variable["feature"]),
        pt_group=pt_label,
        centrality="0-20%",
        stage=stage,
        edges=[float(x) for x in shown_edges],
        density=[float(x) for x in density],
        visible_integral=visible,
        total_integral=total,
        mean=mean,
        rms=rms,
    )


def build_payload(args: argparse.Namespace) -> dict:
    z = np.load(args.matrix, allow_pickle=True)
    missing = [key for key in [*FEATURES, "is_signal", "cluster_Et", "cluster_Eta", "centrality"] if key not in z.files]
    if missing:
        raise SystemExit(f"{args.matrix} missing keys: {missing}")
    labels = np.asarray(z["is_signal"], dtype=np.int8)
    et = np.asarray(z["cluster_Et"], dtype=np.float32)
    eta = np.asarray(z["cluster_Eta"], dtype=np.float32)
    cent = np.asarray(z["centrality"], dtype=np.float32)
    weights, weight_report = compute_ppg12_exact_weights(labels, et, eta)
    score = score_matrix(z, args.model, args.chunk_size)
    base = (
        np.isfinite(score)
        & np.isfinite(et)
        & np.isfinite(cent)
        & np.isin(labels, [0, 1])
        & (et >= 15.0)
        & (et < 35.0)
        & (cent >= CENT_FOCUS[0])
        & (cent < CENT_FOCUS[1])
    )
    tight = score > (args.wp_intercept + args.wp_slope * cent)
    curves = []
    for sample, label_value in (("Signal MC", 1), ("Inclusive MC", 0)):
        sample_mask = base & (labels == label_value)
        for variable in VARIABLES:
            values = np.asarray(z[variable["feature"]], dtype=np.float32)
            for pt_label, pt_bins in PT_GROUPS:
                pt_mask = np.zeros(len(labels), dtype=bool)
                for lo, hi in pt_bins:
                    pt_mask |= (et >= lo) & (et < hi)
                before = sample_mask & pt_mask
                after = before & tight
                curves.append(make_curve(sample, variable, pt_label, "Before preselection", values[before], weights[before]))
                curves.append(make_curve(sample, variable, pt_label, "Tight WP80", values[after], weights[after]))
    return {
        "schema": "THE57_FULL_WEIGHTED_ENERGY_SUM_FEATURE_CURVES_V1",
        "matrix": str(args.matrix),
        "model": str(args.model),
        "full_matrix_rows": int(len(labels)),
        "rows_loaded": int(base.sum()),
        "weight_mode": "ppg12-exact",
        "weighting": weight_report,
        "wp80_formula": {
            "type": "linear",
            "intercept": float(args.wp_intercept),
            "slope": float(args.wp_slope),
            "expression": f"score > {args.wp_intercept:.8f} + {args.wp_slope:.10f} * centrality_percentile",
        },
        "score_key": "score_centAsFeatBase3x3_pt15to35",
        "source_label": "THE-57 full combined weighted training matrix",
        "model_label": "default 14-feature Au+Au baseline with THE-58 cut",
        "centrality_focus": "0-20%",
        "pt_groups": [{"label": label, "fine_bins": bins} for label, bins in PT_GROUPS],
        "curves": [asdict(c) for c in curves],
    }


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--matrix", type=Path, required=True)
    ap.add_argument("--model", type=Path, required=True)
    ap.add_argument("--json-out", type=Path, required=True)
    ap.add_argument("--chunk-size", type=int, default=500_000)
    ap.add_argument("--wp-intercept", type=float, default=0.53471108)
    ap.add_argument("--wp-slope", type=float, default=0.0012284143)
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    payload = build_payload(args)
    if str(args.json_out) == "-":
        print("__THE57_ENERGY_SUM_JSON_BEGIN__")
        print(json.dumps(payload, separators=(",", ":")))
        print("__THE57_ENERGY_SUM_JSON_END__")
    else:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(args.json_out)


if __name__ == "__main__":
    main()
