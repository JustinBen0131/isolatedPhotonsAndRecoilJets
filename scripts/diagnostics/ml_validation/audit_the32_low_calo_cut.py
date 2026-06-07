#!/usr/bin/env python3
"""Forensic audit for the THE-32 low-calo event-quality cut."""

from __future__ import annotations

import argparse
import csv
import json
import math
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np


PRODUCT = "globalEtCent1535_bdt_noIso"
SCORE_COL = f"score_{PRODUCT}"
CENT_BINS = [
    ("0_20", "0-20%", 0.0, 20.0),
    ("20_50", "20-50%", 20.0, 50.0),
    ("50_80", "50-80%", 50.0, 80.0),
]
CENT_SLICES = [(float(lo), float(lo + 5)) for lo in range(0, 80, 5)]
SAMPLES = [
    ("run28_embeddedPhoton12", "Photon12"),
    ("run28_embeddedPhoton20", "Photon20"),
    ("run28_embeddedJet12", "Jet12"),
    ("run28_embeddedJet20", "Jet20"),
    ("run28_embeddedJet30", "Jet30"),
    ("run28_embeddedJet40", "Jet40"),
]
SAMPLE_TO_CODE = {name: idx + 1 for idx, (name, _) in enumerate(SAMPLES)}
CODE_TO_SAMPLE = {idx + 1: name for idx, (name, _) in enumerate(SAMPLES)}
CODE_TO_LABEL = {idx + 1: label for idx, (_, label) in enumerate(SAMPLES)}


@dataclass
class ScoreData:
    score: np.ndarray
    y: np.ndarray
    source_code: np.ndarray
    cent: np.ndarray
    et: np.ndarray
    eta: np.ndarray
    log_calo: np.ndarray
    cemc: np.ndarray
    ihcal: np.ndarray
    ohcal: np.ndarray
    run: np.ndarray
    evt: np.ndarray
    event_source_code: np.ndarray
    event_cent: np.ndarray
    event_log_calo: np.ndarray
    event_cemc: np.ndarray
    event_ihcal: np.ndarray
    event_ohcal: np.ndarray
    event_run: np.ndarray
    event_evt: np.ndarray
    cache_paths: list[Path]


def json_ready(obj):
    if isinstance(obj, dict):
        return {k: json_ready(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [json_ready(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return json_ready(obj.tolist())
    if isinstance(obj, np.generic):
        return obj.item()
    if isinstance(obj, float) and (math.isnan(obj) or math.isinf(obj)):
        return None
    return obj


def sample_codes(source: np.ndarray) -> np.ndarray:
    text = source.astype(str)
    out = np.zeros(len(text), dtype="uint8")
    for name, code in SAMPLE_TO_CODE.items():
        out[np.char.find(text, name) >= 0] = code
    for name, code in SAMPLE_TO_CODE.items():
        short = name.replace("run28_", "")
        out[(out == 0) & (np.char.find(text, short) >= 0)] = code
    return out


def read_cache_list(path: Path) -> list[Path]:
    if path.is_dir():
        manifest = path / "score_caches.list"
        if manifest.exists():
            path = manifest
        else:
            return sorted((path / "score_caches").glob("score_cache_*.npz"))
    paths = [Path(line.strip()) for line in path.read_text().splitlines() if line.strip()]
    if not paths:
        raise SystemExit(f"No score caches listed in {path}")
    return paths


def require_columns(cache: Path, data, cols: list[str]) -> None:
    missing = [col for col in cols if col not in data.files]
    if missing:
        raise SystemExit(f"{cache} is missing required columns: {missing}")


def load_score_data(cache_list: Path, score_col: str) -> tuple[ScoreData, dict]:
    cache_paths = read_cache_list(cache_list)
    required = [
        score_col,
        "is_signal",
        "source_sample",
        "run",
        "evt",
        "centrality",
        "cluster_Et",
        "cluster_Eta",
        "event_calo_cemc_energy",
        "event_calo_ihcal_energy",
        "event_calo_ohcal_energy",
        "event_calo_log10_total_energy_plus1",
    ]
    parts = {name: [] for name in ["score", "y", "src", "cent", "et", "eta", "log", "cemc", "ihcal", "ohcal", "run", "evt"]}
    ev_parts = {name: [] for name in ["src", "cent", "log", "cemc", "ihcal", "ohcal", "run", "evt"]}
    source_rows = {name: 0 for name, _ in SAMPLES}
    source_events_per_cache = {name: 0 for name, _ in SAMPLES}
    raw_rows = 0
    good_rows = 0
    finite_score_rows = 0
    unknown_source_rows = 0

    for idx, cache in enumerate(cache_paths, 1):
        if idx == 1 or idx % 10 == 0 or idx == len(cache_paths):
            print(f"[audit] loading {idx}/{len(cache_paths)} {cache}", flush=True, file=sys.stderr)
        with np.load(cache, allow_pickle=True) as z:
            require_columns(cache, z, required)
            score = z[score_col].astype("float32", copy=False)
            y = z["is_signal"].astype("int8", copy=False)
            src = sample_codes(z["source_sample"])
            cent = z["centrality"].astype("float32", copy=False)
            et = z["cluster_Et"].astype("float32", copy=False)
            eta = z["cluster_Eta"].astype("float32", copy=False)
            log_calo = z["event_calo_log10_total_energy_plus1"].astype("float32", copy=False)
            cemc = z["event_calo_cemc_energy"].astype("float32", copy=False)
            ihcal = z["event_calo_ihcal_energy"].astype("float32", copy=False)
            ohcal = z["event_calo_ohcal_energy"].astype("float32", copy=False)
            run = z["run"].astype("int32", copy=False)
            evt = z["evt"].astype("int64", copy=False)
            raw_rows += len(score)
            finite_score_rows += int(np.isfinite(score).sum())
            unknown_source_rows += int((src == 0).sum())
            good = np.isfinite(cent) & np.isfinite(et) & np.isfinite(eta) & np.isfinite(log_calo) & (src > 0)
            good_rows += int(good.sum())
            if not np.any(good):
                continue
            for code, sample in CODE_TO_SAMPLE.items():
                source_rows[sample] += int(np.sum(good & (src == code)))
            parts["score"].append(score[good])
            parts["y"].append(y[good])
            parts["src"].append(src[good])
            parts["cent"].append(cent[good])
            parts["et"].append(et[good])
            parts["eta"].append(eta[good])
            parts["log"].append(log_calo[good])
            parts["cemc"].append(cemc[good])
            parts["ihcal"].append(ihcal[good])
            parts["ohcal"].append(ohcal[good])
            parts["run"].append(run[good])
            parts["evt"].append(evt[good])

            key = np.empty(int(good.sum()), dtype=[("src", "u1"), ("run", "i4"), ("evt", "i8")])
            key["src"] = src[good]
            key["run"] = run[good]
            key["evt"] = evt[good]
            unique_idx = np.unique(key, return_index=True)[1]
            ev_src = src[good][unique_idx]
            for code, sample in CODE_TO_SAMPLE.items():
                source_events_per_cache[sample] += int(np.sum(ev_src == code))
            ev_parts["src"].append(ev_src)
            ev_parts["cent"].append(cent[good][unique_idx])
            ev_parts["log"].append(log_calo[good][unique_idx])
            ev_parts["cemc"].append(cemc[good][unique_idx])
            ev_parts["ihcal"].append(ihcal[good][unique_idx])
            ev_parts["ohcal"].append(ohcal[good][unique_idx])
            ev_parts["run"].append(run[good][unique_idx])
            ev_parts["evt"].append(evt[good][unique_idx])

    data = ScoreData(
        score=np.concatenate(parts["score"]),
        y=np.concatenate(parts["y"]),
        source_code=np.concatenate(parts["src"]),
        cent=np.concatenate(parts["cent"]),
        et=np.concatenate(parts["et"]),
        eta=np.concatenate(parts["eta"]),
        log_calo=np.concatenate(parts["log"]),
        cemc=np.concatenate(parts["cemc"]),
        ihcal=np.concatenate(parts["ihcal"]),
        ohcal=np.concatenate(parts["ohcal"]),
        run=np.concatenate(parts["run"]),
        evt=np.concatenate(parts["evt"]),
        event_source_code=np.concatenate(ev_parts["src"]),
        event_cent=np.concatenate(ev_parts["cent"]),
        event_log_calo=np.concatenate(ev_parts["log"]),
        event_cemc=np.concatenate(ev_parts["cemc"]),
        event_ihcal=np.concatenate(ev_parts["ihcal"]),
        event_ohcal=np.concatenate(ev_parts["ohcal"]),
        event_run=np.concatenate(ev_parts["run"]),
        event_evt=np.concatenate(ev_parts["evt"]),
        cache_paths=cache_paths,
    )
    ev_key = np.empty(len(data.event_cent), dtype=[("src", "u1"), ("run", "i4"), ("evt", "i8")])
    ev_key["src"] = data.event_source_code
    ev_key["run"] = data.event_run
    ev_key["evt"] = data.event_evt
    global_unique_events = int(len(np.unique(ev_key)))
    audit = {
        "raw_candidate_rows": raw_rows,
        "good_candidate_rows_after_finite_filter": good_rows,
        "finite_score_rows": finite_score_rows,
        "unknown_source_rows": unknown_source_rows,
        "event_rows_after_per_cache_dedup": int(len(data.event_cent)),
        "global_unique_event_keys_after_per_cache_dedup": global_unique_events,
        "duplicate_event_keys_after_per_cache_dedup": int(len(data.event_cent) - global_unique_events),
        "source_candidate_rows_after_finite_filter": source_rows,
        "source_event_rows_after_per_cache_dedup": source_events_per_cache,
        "centrality_min": float(np.nanmin(data.event_cent)),
        "centrality_max": float(np.nanmax(data.event_cent)),
    }
    return data, audit


def robust_envelope(events_cent: np.ndarray, events_log: np.ndarray, *, mad_scale: float, quantile_floor: float) -> list[dict]:
    out = []
    finite = np.isfinite(events_cent) & np.isfinite(events_log) & (events_cent >= 0.0) & (events_cent < 80.0)
    for lo, hi in CENT_SLICES:
        mask = finite & (events_cent >= lo) & (events_cent < hi)
        vals = events_log[mask].astype("float64", copy=False)
        med = float(np.nanmedian(vals))
        mad = float(np.nanmedian(np.abs(vals - med)))
        sigma = 1.4826 * mad
        median_minus = med - mad_scale * sigma
        q = float(np.nanquantile(vals, quantile_floor))
        threshold = max(median_minus, q)
        out.append(
            {
                "cent_lo": lo,
                "cent_hi": hi,
                "n_events": int(len(vals)),
                "median": med,
                "mad_sigma": sigma,
                "median_minus_scaled_MAD": median_minus,
                "quantile_floor": q,
                "threshold": float(threshold),
                "winning_term": "quantile_floor" if q >= median_minus else "median_minus_scaled_MAD",
            }
        )
    return out


def threshold_for_cent(cent: np.ndarray, envelope: list[dict]) -> np.ndarray:
    out = np.full(len(cent), np.nan, dtype="float32")
    for row in envelope:
        mask = (cent >= float(row["cent_lo"])) & (cent < float(row["cent_hi"]))
        out[mask] = float(row["threshold"])
    return out


def cut_masks(data: ScoreData, envelope: list[dict]) -> tuple[np.ndarray, np.ndarray]:
    cand_t = threshold_for_cent(data.cent, envelope)
    ev_t = threshold_for_cent(data.event_cent, envelope)
    cand = np.isfinite(cand_t) & np.isfinite(data.log_calo) & (data.log_calo < cand_t)
    ev = np.isfinite(ev_t) & np.isfinite(data.event_log_calo) & (data.event_log_calo < ev_t)
    return cand, ev


def centrality_mask(values: np.ndarray, lo: float, hi: float) -> np.ndarray:
    return np.isfinite(values) & (values >= lo) & (values < hi)


def normalize_mean_one(values: np.ndarray) -> np.ndarray:
    finite = np.isfinite(values) & (values > 0.0)
    if not np.any(finite):
        return np.ones(len(values), dtype="float64")
    mean = float(values[finite].mean())
    return np.where(finite, values / mean, 1.0) if mean > 0.0 and math.isfinite(mean) else np.ones(len(values), dtype="float64")


def inverse_pdf_weights(values: np.ndarray, *, n_bins: int, fixed_range=None, weight_cap: float | None = None) -> np.ndarray:
    values = np.asarray(values, dtype="float64")
    weights = np.ones(len(values), dtype="float64")
    finite = np.isfinite(values)
    if finite.sum() < max(10, n_bins):
        return weights
    vals = values[finite]
    lo, hi = (float(np.min(vals)), float(np.max(vals))) if fixed_range is None else (float(fixed_range[0]), float(fixed_range[1]))
    if not math.isfinite(lo) or not math.isfinite(hi) or hi <= lo:
        return weights
    edges = np.linspace(lo, hi, n_bins + 1)
    hist, _ = np.histogram(vals, bins=edges, density=True)
    hist = hist.astype("float64") * float(n_bins)
    idx = np.clip(np.searchsorted(edges, vals, side="right") - 1, 0, n_bins - 1)
    local = 1.0 / np.clip(hist[idx], 1.0e-3, None)
    if weight_cap is not None:
        local = np.clip(local, a_min=None, a_max=float(weight_cap))
    weights[finite] = normalize_mean_one(local)
    return weights


def ppg12_exact_weights(data: ScoreData) -> np.ndarray:
    y = data.y.astype("int32", copy=False)
    weights = np.ones(len(y), dtype="float64")
    valid = np.isin(y, [0, 1]) & np.isfinite(data.score)
    counts = {cls: int(np.sum(valid & (y == cls))) for cls in (0, 1)}
    for cls in (0, 1):
        mask = valid & (y == cls)
        weights[mask] *= float(counts[0] + counts[1]) / (2.0 * float(counts[cls]))
        weights[mask] *= inverse_pdf_weights(data.eta[mask], n_bins=20, fixed_range=(-0.7, 0.7))
        weights[mask] *= inverse_pdf_weights(data.et[mask], n_bins=20, weight_cap=800.0)
    return weights


def weighted_hist(values: np.ndarray, weights: np.ndarray, bins: np.ndarray) -> np.ndarray:
    good = np.isfinite(values) & np.isfinite(weights) & (weights > 0.0)
    return np.histogram(values[good], bins=bins, weights=weights[good])[0].astype("float64")


def binned_auc(sig_counts: np.ndarray, bkg_counts: np.ndarray) -> float:
    sig_total = float(np.sum(sig_counts))
    bkg_total = float(np.sum(bkg_counts))
    if sig_total <= 0.0 or bkg_total <= 0.0:
        return math.nan
    bkg_below = np.cumsum(bkg_counts) - bkg_counts
    return float((np.sum(sig_counts * bkg_below) + 0.5 * np.sum(sig_counts * bkg_counts)) / (sig_total * bkg_total))


def weighted_wp_threshold(score: np.ndarray, weight: np.ndarray, target_eff: float = 0.80) -> float:
    good = np.isfinite(score) & np.isfinite(weight) & (weight > 0.0)
    if not np.any(good):
        return math.nan
    s = score[good].astype("float64", copy=False)
    w = weight[good].astype("float64", copy=False)
    order = np.argsort(s)
    cdf = np.cumsum(w[order]) / float(np.sum(w))
    q = 1.0 - target_eff
    return float(s[order][min(len(s) - 1, int(np.searchsorted(cdf, q, side="left")))])


def weighted_rate(score: np.ndarray, weight: np.ndarray, threshold: float) -> float:
    good = np.isfinite(score) & np.isfinite(weight) & (weight > 0.0)
    return float(np.sum(weight[good & (score >= threshold)]) / np.sum(weight[good])) if np.any(good) and math.isfinite(threshold) else math.nan


def closure_rows(data: ScoreData, cand_reject: np.ndarray, weights: np.ndarray) -> list[dict]:
    rows = []
    score_bins = np.linspace(0.0, 1.0, 81)
    is_scored = np.isfinite(data.score)
    is_signal = np.isin(data.source_code, [SAMPLE_TO_CODE["run28_embeddedPhoton12"], SAMPLE_TO_CODE["run28_embeddedPhoton20"]]) & (data.y == 1)
    is_inclusive = np.isin(data.source_code, [SAMPLE_TO_CODE[s] for s in ["run28_embeddedJet12", "run28_embeddedJet20", "run28_embeddedJet30", "run28_embeddedJet40"]])
    for key, label, lo, hi in CENT_BINS:
        before = centrality_mask(data.cent, lo, hi) & is_scored
        after = before & ~cand_reject
        sig_before = before & is_signal
        inc_before = before & is_inclusive
        sig_after = after & is_signal
        inc_after = after & is_inclusive
        sig_h_before = weighted_hist(data.score[sig_before], weights[sig_before], score_bins)
        inc_h_before = weighted_hist(data.score[inc_before], weights[inc_before], score_bins)
        sig_h_after = weighted_hist(data.score[sig_after], weights[sig_after], score_bins)
        inc_h_after = weighted_hist(data.score[inc_after], weights[inc_after], score_bins)
        wp_before = weighted_wp_threshold(data.score[sig_before], weights[sig_before], 0.80)
        wp_after = weighted_wp_threshold(data.score[sig_after], weights[sig_after], 0.80)
        rows.append(
            {
                "centrality_key": key,
                "centrality_bin": label,
                "candidate_total": int(before.sum()),
                "candidate_rejected": int((before & cand_reject).sum()),
                "candidate_rejected_fraction": float((before & cand_reject).sum() / max(1, before.sum())),
                "signal_before": int(sig_before.sum()),
                "signal_after": int(sig_after.sum()),
                "inclusive_before": int(inc_before.sum()),
                "inclusive_after": int(inc_after.sum()),
                "auc_before": binned_auc(sig_h_before, inc_h_before),
                "auc_after": binned_auc(sig_h_after, inc_h_after),
                "auc_shift": binned_auc(sig_h_after, inc_h_after) - binned_auc(sig_h_before, inc_h_before),
                "wp80_threshold_before": wp_before,
                "wp80_threshold_after": wp_after,
                "wp80_inclusive_rate_before": weighted_rate(data.score[inc_before], weights[inc_before], wp_before),
                "wp80_inclusive_rate_after": weighted_rate(data.score[inc_after], weights[inc_after], wp_after),
                "wp80_inclusive_rate_shift": weighted_rate(data.score[inc_after], weights[inc_after], wp_after)
                - weighted_rate(data.score[inc_before], weights[inc_before], wp_before),
            }
        )
    return rows


def rejection_rows(data: ScoreData, cand_reject: np.ndarray, event_reject: np.ndarray) -> list[dict]:
    rows = []
    for code, sample in CODE_TO_SAMPLE.items():
        for _, label, lo, hi in CENT_BINS:
            ev_mask = (data.event_source_code == code) & centrality_mask(data.event_cent, lo, hi)
            ca_mask = (data.source_code == code) & centrality_mask(data.cent, lo, hi)
            rows.append(
                {
                    "source_sample": sample,
                    "source_label": CODE_TO_LABEL[code],
                    "centrality_bin": label,
                    "event_total": int(ev_mask.sum()),
                    "event_rejected": int((ev_mask & event_reject).sum()),
                    "event_rejected_fraction": float((ev_mask & event_reject).sum() / max(1, ev_mask.sum())),
                    "candidate_total": int(ca_mask.sum()),
                    "candidate_rejected": int((ca_mask & cand_reject).sum()),
                    "candidate_rejected_fraction": float((ca_mask & cand_reject).sum() / max(1, ca_mask.sum())),
                }
            )
    return rows


def component_medians(data: ScoreData, event_reject: np.ndarray) -> list[dict]:
    rows = []
    total = data.event_cemc + data.event_ihcal + data.event_ohcal
    for _, label, lo, hi in CENT_BINS:
        base = centrality_mask(data.event_cent, lo, hi)
        for name, arr in [
            ("CEMC", data.event_cemc),
            ("IHCal", data.event_ihcal),
            ("OHCal", data.event_ohcal),
            ("total_calo", total),
            ("log10_total_plus1", data.event_log_calo),
        ]:
            ret = base & ~event_reject & np.isfinite(arr)
            rej = base & event_reject & np.isfinite(arr)
            retained = float(np.nanmedian(arr[ret])) if np.any(ret) else math.nan
            rejected = float(np.nanmedian(arr[rej])) if np.any(rej) else math.nan
            rows.append(
                {
                    "centrality_bin": label,
                    "component": name,
                    "retained_median": retained,
                    "rejected_median": rejected,
                    "rejected_over_retained": rejected / retained if retained > 0 and math.isfinite(rejected) else math.nan,
                    "retained_events": int(ret.sum()),
                    "rejected_events": int(rej.sum()),
                }
            )
    return rows


def sensitivity(data: ScoreData, weights: np.ndarray, mad_scales: list[float], floors: list[float]) -> list[dict]:
    rows = []
    for mad_scale in mad_scales:
        for floor in floors:
            env = robust_envelope(data.event_cent, data.event_log_calo, mad_scale=mad_scale, quantile_floor=floor)
            cand_reject, event_reject = cut_masks(data, env)
            closures = closure_rows(data, cand_reject, weights)
            event_total_080 = int(np.sum(centrality_mask(data.event_cent, 0, 80)))
            event_rej_080 = int(np.sum(centrality_mask(data.event_cent, 0, 80) & event_reject))
            row = {
                "mad_scale": mad_scale,
                "quantile_floor": floor,
                "total_rejected_events": event_rej_080,
                "total_event_rejected_fraction_0_80": event_rej_080 / max(1, event_total_080),
            }
            for _, label, lo, hi in CENT_BINS:
                m = centrality_mask(data.event_cent, lo, hi)
                row[f"event_rejected_fraction_{label}"] = float(np.sum(m & event_reject) / max(1, np.sum(m)))
            for code, sample in CODE_TO_SAMPLE.items():
                m = (data.event_source_code == code) & centrality_mask(data.event_cent, 0, 80)
                row[f"event_rejected_fraction_{CODE_TO_LABEL[code]}"] = float(np.sum(m & event_reject) / max(1, np.sum(m)))
            for c in closures:
                key = c["centrality_key"]
                row[f"auc_shift_{key}"] = c["auc_shift"]
                row[f"wp80_inclusive_rate_shift_{key}"] = c["wp80_inclusive_rate_shift"]
            rows.append(row)
    return rows


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def git_status() -> str:
    try:
        return subprocess.check_output(["git", "status", "--short"], text=True, stderr=subprocess.STDOUT)
    except Exception as exc:
        return f"git status unavailable: {exc}"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--score-cache-list", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--score-col", default=SCORE_COL)
    ap.add_argument("--mad-scale", type=float, default=5.0)
    ap.add_argument("--quantile-floor", type=float, default=0.001)
    ap.add_argument("--emit-only", action="store_true", help="Print JSON only; do not write output files.")
    args = ap.parse_args()
    if not args.emit_only:
        args.outdir.mkdir(parents=True, exist_ok=True)

    data, row_audit = load_score_data(args.score_cache_list, args.score_col)
    envelope = robust_envelope(data.event_cent, data.event_log_calo, mad_scale=args.mad_scale, quantile_floor=args.quantile_floor)
    cand_reject, event_reject = cut_masks(data, envelope)
    weights = ppg12_exact_weights(data)
    rejection = rejection_rows(data, cand_reject, event_reject)
    components = component_medians(data, event_reject)
    closure = closure_rows(data, cand_reject, weights)
    sens = sensitivity(data, weights, [4.0, 5.0, 6.0], [0.0005, 0.001, 0.005])

    paths = {
        "audit_json": args.outdir / "the32_low_calo_forensic_audit_v1.json",
        "envelope_csv": args.outdir / "the32_low_calo_envelope_audit_v1.csv",
        "component_csv": args.outdir / "the32_low_calo_component_medians_v1.csv",
        "sensitivity_csv": args.outdir / "the32_low_calo_sensitivity_v1.csv",
        "rejection_csv": args.outdir / "the32_low_calo_rejection_audit_v1.csv",
        "closure_csv": args.outdir / "the32_low_calo_closure_audit_v1.csv",
    }
    audit = {
        "schema": "THE32_LOW_CALO_FORENSIC_AUDIT_V1",
        "score_cache_list": str(args.score_cache_list),
        "score_column": args.score_col,
        "cut_contract": {
            "direct_threshold_columns": [
                "centrality",
                "event_calo_log10_total_energy_plus1",
            ],
            "event_calo_components_for_interpretation": [
                "event_calo_cemc_energy",
                "event_calo_ihcal_energy",
                "event_calo_ohcal_energy",
            ],
            "excluded_from_threshold": [
                "source_sample",
                "is_signal",
                "truth photon label",
                "truth isolation label",
                "BDT score",
                "candidate pT / cluster_Et",
                "candidate shower-shape variables",
                "train/test split",
                "sample weight",
            ],
            "formula": "threshold per 5% centrality slice = max(median(log_calo) - mad_scale * 1.4826*MAD(log_calo), quantile_floor(log_calo))",
            "mad_scale": args.mad_scale,
            "quantile_floor": args.quantile_floor,
        },
        "row_audit": row_audit,
        "envelope": envelope,
        "rejection_rows": rejection,
        "component_medians": components,
        "closure_rows": closure,
        "sensitivity_rows": sens,
        "outputs": {key: str(path) for key, path in paths.items()},
        "commands": {
            "local_script": "scripts/diagnostics/ml_validation/audit_the32_low_calo_cut.py",
            "example": f"python3 scripts/diagnostics/ml_validation/audit_the32_low_calo_cut.py --score-cache-list {args.score_cache_list} --outdir {args.outdir}",
        },
        "git_status_short_at_runtime": git_status(),
    }
    if not args.emit_only:
        write_csv(paths["envelope_csv"], envelope)
        write_csv(paths["component_csv"], components)
        write_csv(paths["sensitivity_csv"], sens)
        write_csv(paths["rejection_csv"], rejection)
        write_csv(paths["closure_csv"], closure)
        paths["audit_json"].write_text(json.dumps(json_ready(audit), indent=2, sort_keys=True))
    print(json.dumps(json_ready(audit), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
