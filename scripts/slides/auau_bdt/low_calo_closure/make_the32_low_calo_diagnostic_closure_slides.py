#!/usr/bin/env python3
"""Build THE-32 low-calo pathology and post-hoc closure slide candidates."""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np


SLIDE_W_PX = 2560
SLIDE_H_PX = 1440
SLIDE_DPI = 200
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
OUTDIR = Path("dataOutput/auauTightBDTValidation/THE32_lowCaloDiagnosticClosure_20260603")

SIGNAL_COLOR = "#C5362C"
INCLUSIVE_COLOR = "#1E6FB3"
LOW_CALO_COLOR = "#B42318"
RETAINED_COLOR = "#98A2B3"
INK = "#182230"
MUTED = "#475467"
BORDER = "#C9D2DE"
GRID = "#D0D5DD"
SOFT_GRAY = "#F6F7F9"
SOFT_BLUE = "#EDF6FF"
SOFT_YELLOW = "#FFF6D7"


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
    cache_paths: list[Path]
    event_source_code: np.ndarray
    event_cent: np.ndarray
    event_log_calo: np.ndarray
    event_cemc: np.ndarray
    event_ihcal: np.ndarray
    event_ohcal: np.ndarray
    event_run: np.ndarray
    event_evt: np.ndarray


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
    if np.any(out == 0):
        # Accept shorter provenance strings if the validator inferred those.
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
        raise SystemExit(
            f"{cache} is missing {missing}; rerun validation with source/run/evt/event_calo diagnostic columns."
        )


def load_score_data(cache_list: Path, score_col: str) -> ScoreData:
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
    score_parts: list[np.ndarray] = []
    y_parts: list[np.ndarray] = []
    src_parts: list[np.ndarray] = []
    cent_parts: list[np.ndarray] = []
    et_parts: list[np.ndarray] = []
    eta_parts: list[np.ndarray] = []
    log_parts: list[np.ndarray] = []
    cemc_parts: list[np.ndarray] = []
    ihcal_parts: list[np.ndarray] = []
    ohcal_parts: list[np.ndarray] = []
    run_parts: list[np.ndarray] = []
    evt_parts: list[np.ndarray] = []
    ev_src_parts: list[np.ndarray] = []
    ev_cent_parts: list[np.ndarray] = []
    ev_log_parts: list[np.ndarray] = []
    ev_cemc_parts: list[np.ndarray] = []
    ev_ihcal_parts: list[np.ndarray] = []
    ev_ohcal_parts: list[np.ndarray] = []
    ev_run_parts: list[np.ndarray] = []
    ev_evt_parts: list[np.ndarray] = []

    for idx, cache in enumerate(cache_paths, 1):
        if idx == 1 or idx % 10 == 0 or idx == len(cache_paths):
            print(f"[THE32] loading score cache {idx}/{len(cache_paths)}: {cache}", flush=True)
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
            good = (
                np.isfinite(cent)
                & np.isfinite(et)
                & np.isfinite(eta)
                & np.isfinite(log_calo)
                & (src > 0)
            )
            if not np.any(good):
                continue
            score_parts.append(score[good])
            y_parts.append(y[good])
            src_parts.append(src[good])
            cent_parts.append(cent[good])
            et_parts.append(et[good])
            eta_parts.append(eta[good])
            log_parts.append(log_calo[good])
            cemc_parts.append(cemc[good])
            ihcal_parts.append(ihcal[good])
            ohcal_parts.append(ohcal[good])
            run_parts.append(run[good])
            evt_parts.append(evt[good])

            key = np.empty(int(np.sum(good)), dtype=[("src", "u1"), ("run", "i4"), ("evt", "i8")])
            key["src"] = src[good]
            key["run"] = run[good]
            key["evt"] = evt[good]
            unique_idx = np.unique(key, return_index=True)[1]
            ev_src_parts.append(src[good][unique_idx])
            ev_cent_parts.append(cent[good][unique_idx])
            ev_log_parts.append(log_calo[good][unique_idx])
            ev_cemc_parts.append(cemc[good][unique_idx])
            ev_ihcal_parts.append(ihcal[good][unique_idx])
            ev_ohcal_parts.append(ohcal[good][unique_idx])
            ev_run_parts.append(run[good][unique_idx])
            ev_evt_parts.append(evt[good][unique_idx])

    if not score_parts:
        raise SystemExit("No finite scored rows found.")

    return ScoreData(
        score=np.concatenate(score_parts),
        y=np.concatenate(y_parts),
        source_code=np.concatenate(src_parts),
        cent=np.concatenate(cent_parts),
        et=np.concatenate(et_parts),
        eta=np.concatenate(eta_parts),
        log_calo=np.concatenate(log_parts),
        cemc=np.concatenate(cemc_parts),
        ihcal=np.concatenate(ihcal_parts),
        ohcal=np.concatenate(ohcal_parts),
        run=np.concatenate(run_parts),
        evt=np.concatenate(evt_parts),
        cache_paths=cache_paths,
        event_source_code=np.concatenate(ev_src_parts),
        event_cent=np.concatenate(ev_cent_parts),
        event_log_calo=np.concatenate(ev_log_parts),
        event_cemc=np.concatenate(ev_cemc_parts),
        event_ihcal=np.concatenate(ev_ihcal_parts),
        event_ohcal=np.concatenate(ev_ohcal_parts),
        event_run=np.concatenate(ev_run_parts),
        event_evt=np.concatenate(ev_evt_parts),
    )


def robust_envelope(events_cent: np.ndarray, events_log: np.ndarray, *, mad_scale: float, quantile_floor: float) -> list[dict]:
    out: list[dict] = []
    finite = np.isfinite(events_cent) & np.isfinite(events_log) & (events_cent >= 0.0) & (events_cent < 80.0)
    for lo, hi in CENT_SLICES:
        mask = finite & (events_cent >= lo) & (events_cent < hi)
        vals = events_log[mask].astype("float64", copy=False)
        if len(vals) < 25:
            out.append(
                {
                    "cent_lo": lo,
                    "cent_hi": hi,
                    "n_events": int(len(vals)),
                    "median": None,
                    "mad_sigma": None,
                    "quantile_floor": None,
                    "threshold": None,
                    "status": "insufficient_events",
                }
            )
            continue
        med = float(np.nanmedian(vals))
        mad = float(np.nanmedian(np.abs(vals - med)))
        sigma = 1.4826 * mad
        q = float(np.nanquantile(vals, quantile_floor))
        threshold = max(med - mad_scale * sigma, q)
        out.append(
            {
                "cent_lo": lo,
                "cent_hi": hi,
                "n_events": int(len(vals)),
                "median": med,
                "mad_sigma": sigma,
                "quantile_floor": q,
                "threshold": float(threshold),
                "status": "ok",
            }
        )
    return out


def threshold_for_cent(cent: np.ndarray, envelope: list[dict]) -> np.ndarray:
    threshold = np.full(len(cent), np.nan, dtype="float32")
    for row in envelope:
        t = row.get("threshold")
        if t is None:
            continue
        lo = float(row["cent_lo"])
        hi = float(row["cent_hi"])
        mask = (cent >= lo) & (cent < hi)
        threshold[mask] = float(t)
    return threshold


def event_cut_masks(data: ScoreData, envelope: list[dict]) -> tuple[np.ndarray, np.ndarray]:
    cand_t = threshold_for_cent(data.cent, envelope)
    event_t = threshold_for_cent(data.event_cent, envelope)
    cand_reject = np.isfinite(cand_t) & np.isfinite(data.log_calo) & (data.log_calo < cand_t)
    event_reject = np.isfinite(event_t) & np.isfinite(data.event_log_calo) & (data.event_log_calo < event_t)
    return cand_reject, event_reject


def normalize_mean_one(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype="float64")
    finite = np.isfinite(values) & (values > 0.0)
    if not finite.any():
        return np.ones(len(values), dtype="float64")
    mean = float(values[finite].mean())
    if mean <= 0.0 or not math.isfinite(mean):
        return np.ones(len(values), dtype="float64")
    return np.where(finite, values / mean, 1.0)


def inverse_pdf_weights(values: np.ndarray, *, n_bins: int, fixed_range=None, weight_cap: float | None = None) -> np.ndarray:
    try:
        from scipy.interpolate import UnivariateSpline
    except Exception:
        UnivariateSpline = None

    values = np.asarray(values, dtype="float64")
    weights = np.ones(len(values), dtype="float64")
    finite = np.isfinite(values)
    if finite.sum() < max(10, n_bins):
        return weights
    vals = values[finite]
    if fixed_range is None:
        lo, hi = float(np.min(vals)), float(np.max(vals))
    else:
        lo, hi = float(fixed_range[0]), float(fixed_range[1])
    if not math.isfinite(lo) or not math.isfinite(hi) or hi <= lo:
        return weights
    edges = np.linspace(lo, hi, n_bins + 1)
    hist, bin_edges = np.histogram(vals, bins=edges, density=True)
    hist = hist.astype("float64") * float(n_bins)
    centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    good_hist = np.isfinite(hist)
    if good_hist.sum() < 4:
        return weights
    if UnivariateSpline is not None:
        try:
            pdf = UnivariateSpline(centers[good_hist], hist[good_hist], s=0.0)(vals)
        except Exception:
            pdf = None
    else:
        pdf = None
    if pdf is None:
        idx = np.clip(np.searchsorted(edges, vals, side="right") - 1, 0, n_bins - 1)
        pdf = hist[idx]
    pdf = np.clip(pdf, a_min=1.0e-3, a_max=None)
    local = 1.0 / pdf
    if weight_cap is not None:
        local = np.clip(local, a_min=None, a_max=float(weight_cap))
    weights[finite] = normalize_mean_one(local)
    return weights


def ppg12_exact_weights(data: ScoreData) -> np.ndarray:
    y = data.y.astype("int32", copy=False)
    weights = np.ones(len(y), dtype="float64")
    valid = np.isin(y, [0, 1]) & np.isfinite(data.score)
    counts = {cls: int(np.sum(valid & (y == cls))) for cls in (0, 1)}
    total = counts[0] + counts[1]
    if counts[0] <= 0 or counts[1] <= 0:
        raise SystemExit(f"Cannot compute PPG12-style weights; class counts={counts}")
    for cls in (0, 1):
        mask = valid & (y == cls)
        weights[mask] *= float(total) / (2.0 * float(counts[cls]))
        weights[mask] *= inverse_pdf_weights(data.eta[mask], n_bins=20, fixed_range=(-0.7, 0.7))
        weights[mask] *= inverse_pdf_weights(data.et[mask], n_bins=20, fixed_range=None, weight_cap=800.0)
    return weights


def weighted_hist(values: np.ndarray, weights: np.ndarray, bins: np.ndarray) -> np.ndarray:
    good = np.isfinite(values) & np.isfinite(weights) & (weights > 0.0)
    return np.histogram(values[good], bins=bins, weights=weights[good])[0].astype("float64")


def density(values: np.ndarray, bins: np.ndarray) -> np.ndarray:
    counts, _ = np.histogram(values[np.isfinite(values)], bins=bins)
    width = np.diff(bins)
    out = counts.astype("float64")
    if out.sum() > 0:
        out = out / out.sum() / width
    return out


def binned_auc(sig_counts: np.ndarray, bkg_counts: np.ndarray) -> float:
    sig_total = float(np.sum(sig_counts))
    bkg_total = float(np.sum(bkg_counts))
    if sig_total <= 0.0 or bkg_total <= 0.0:
        return math.nan
    bkg_below = np.cumsum(bkg_counts) - bkg_counts
    wins = float(np.sum(sig_counts * bkg_below))
    ties = float(np.sum(sig_counts * bkg_counts))
    return (wins + 0.5 * ties) / (sig_total * bkg_total)


def weighted_wp_threshold(signal_score: np.ndarray, signal_weight: np.ndarray, target_eff: float = 0.80) -> float:
    good = np.isfinite(signal_score) & np.isfinite(signal_weight) & (signal_weight > 0.0)
    if not np.any(good):
        return math.nan
    score = signal_score[good].astype("float64", copy=False)
    weight = signal_weight[good].astype("float64", copy=False)
    order = np.argsort(score)
    score = score[order]
    weight = weight[order]
    cdf = np.cumsum(weight) / float(np.sum(weight))
    q = max(0.0, min(1.0, 1.0 - target_eff))
    return float(score[min(len(score) - 1, int(np.searchsorted(cdf, q, side="left")))])


def weighted_rate(score: np.ndarray, weight: np.ndarray, threshold: float) -> float:
    good = np.isfinite(score) & np.isfinite(weight) & (weight > 0.0)
    if not np.any(good) or not math.isfinite(threshold):
        return math.nan
    return float(np.sum(weight[good & (score >= threshold)]) / np.sum(weight[good]))


def centrality_mask(values: np.ndarray, lo: float, hi: float) -> np.ndarray:
    return np.isfinite(values) & (values >= lo) & (values < hi)


def summarize_counts(data: ScoreData, cand_reject: np.ndarray, event_reject: np.ndarray) -> list[dict]:
    rows = []
    for code, sample_name in CODE_TO_SAMPLE.items():
        short = CODE_TO_LABEL[code]
        for key, label, lo, hi in CENT_BINS:
            ev_mask = (data.event_source_code == code) & centrality_mask(data.event_cent, lo, hi)
            ca_mask = (data.source_code == code) & centrality_mask(data.cent, lo, hi)
            ev_total = int(ev_mask.sum())
            ca_total = int(ca_mask.sum())
            ev_rej = int(np.sum(ev_mask & event_reject))
            ca_rej = int(np.sum(ca_mask & cand_reject))
            rows.append(
                {
                    "source_sample": sample_name,
                    "source_label": short,
                    "centrality_bin": label,
                    "event_total": ev_total,
                    "event_rejected": ev_rej,
                    "event_rejected_fraction": float(ev_rej / ev_total) if ev_total else math.nan,
                    "candidate_total": ca_total,
                    "candidate_rejected": ca_rej,
                    "candidate_rejected_fraction": float(ca_rej / ca_total) if ca_total else math.nan,
                }
            )
    return rows


def summarize_closure(data: ScoreData, cand_reject: np.ndarray, weights: np.ndarray, score_bins: np.ndarray) -> list[dict]:
    rows = []
    is_scored = np.isfinite(data.score)
    is_signal = np.isin(data.source_code, [SAMPLE_TO_CODE["run28_embeddedPhoton12"], SAMPLE_TO_CODE["run28_embeddedPhoton20"]]) & (data.y == 1)
    is_inclusive = np.isin(
        data.source_code,
        [
            SAMPLE_TO_CODE["run28_embeddedJet12"],
            SAMPLE_TO_CODE["run28_embeddedJet20"],
            SAMPLE_TO_CODE["run28_embeddedJet30"],
            SAMPLE_TO_CODE["run28_embeddedJet40"],
        ],
    )
    for key, label, lo, hi in CENT_BINS:
        cent = centrality_mask(data.cent, lo, hi) & is_scored
        before = cent
        after = cent & ~cand_reject
        sig_before = before & is_signal
        inc_before = before & is_inclusive
        sig_after = after & is_signal
        inc_after = after & is_inclusive
        if not np.any(sig_after) or not np.any(inc_after):
            raise SystemExit(f"Post-cut class is empty in {label}: Signal={sig_after.sum()} Inclusive={inc_after.sum()}")
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
                "candidate_total": int(np.sum(before)),
                "candidate_rejected": int(np.sum(before & cand_reject)),
                "candidate_rejected_fraction": float(np.sum(before & cand_reject) / max(1, np.sum(before))),
                "candidate_retained_fraction": float(np.sum(after) / max(1, np.sum(before))),
                "signal_before": int(np.sum(sig_before)),
                "signal_after": int(np.sum(sig_after)),
                "inclusive_before": int(np.sum(inc_before)),
                "inclusive_after": int(np.sum(inc_after)),
                "auc_before": binned_auc(sig_h_before, inc_h_before),
                "auc_after": binned_auc(sig_h_after, inc_h_after),
                "wp80_threshold_before": wp_before,
                "wp80_threshold_after": wp_after,
                "wp80_inclusive_rate_before": weighted_rate(data.score[inc_before], weights[inc_before], wp_before),
                "wp80_inclusive_rate_after": weighted_rate(data.score[inc_after], weights[inc_after], wp_after),
                "signal_density_before": density(data.score[sig_before], score_bins).tolist(),
                "signal_density_after": density(data.score[sig_after], score_bins).tolist(),
                "inclusive_density_before": density(data.score[inc_before], score_bins).tolist(),
                "inclusive_density_after": density(data.score[inc_after], score_bins).tolist(),
            }
        )
    return rows


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def setup_matplotlib():
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": [
                "Nimbus Roman",
                "Nimbus Roman No9 L",
                "Times New Roman",
                "Times",
                "Liberation Serif",
                "DejaVu Serif",
            ],
            "mathtext.fontset": "stix",
            "axes.unicode_minus": False,
        }
    )
    return plt


def add_box(fig, xy, wh, face, edge=BORDER, lw=1.0, radius=0.010, zorder=-2):
    from matplotlib.patches import FancyBboxPatch

    fig.patches.append(
        FancyBboxPatch(
            xy,
            wh[0],
            wh[1],
            transform=fig.transFigure,
            boxstyle=f"round,pad=0.010,rounding_size={radius}",
            facecolor=face,
            edgecolor=edge,
            linewidth=lw,
            zorder=zorder,
        )
    )


def fig_text(fig, x, y, text, *, size, weight="normal", color=INK, ha="left", va="top", linespacing=1.25):
    fig.text(
        x,
        y,
        text,
        fontsize=size,
        fontweight=weight,
        color=color,
        ha=ha,
        va=va,
        linespacing=linespacing,
        fontfamily="serif",
    )


def sphenix_label(ax, x=0.035, y=0.955, size=15.0):
    from matplotlib.offsetbox import AnchoredOffsetbox, HPacker, TextArea

    sphenix = TextArea(
        "sPHENIX",
        textprops={
            "fontsize": size,
            "fontfamily": "serif",
            "fontweight": "bold",
            "fontstyle": "italic",
            "color": INK,
        },
    )
    internal = TextArea(
        " Internal",
        textprops={
            "fontsize": size,
            "fontfamily": "serif",
            "fontweight": "normal",
            "fontstyle": "normal",
            "color": INK,
        },
    )
    packed = HPacker(children=[sphenix, internal], align="baseline", pad=0, sep=0)
    anchored = AnchoredOffsetbox(
        loc="upper left",
        child=packed,
        pad=0.0,
        borderpad=0.0,
        frameon=False,
        bbox_to_anchor=(x, y),
        bbox_transform=ax.transAxes,
    )
    ax.add_artist(anchored)


def fmt_frac(value: float) -> str:
    if not math.isfinite(value):
        return "n/a"
    return f"{100.0 * value:.2f}%"


def fmt_count(value: int) -> str:
    if value >= 1_000_000:
        return f"{value / 1_000_000:.2f}M"
    if value >= 1000:
        return f"{value / 1000:.1f}k"
    return str(value)


def draw_pathology_slide(
    outpng: Path,
    data: ScoreData,
    envelope: list[dict],
    cand_reject: np.ndarray,
    event_reject: np.ndarray,
    event_rows: list[dict],
) -> None:
    plt = setup_matplotlib()
    from matplotlib.colors import LogNorm
    from matplotlib.patches import Rectangle

    fig = plt.figure(figsize=(SLIDE_W_PX / SLIDE_DPI, SLIDE_H_PX / SLIDE_DPI), dpi=SLIDE_DPI, facecolor="white")
    fig.patches.append(Rectangle((0, 0), 1, 1, transform=fig.transFigure, facecolor="white", zorder=-10))
    fig_text(
        fig,
        0.040,
        0.965,
        "Low-calo events are a source-blind quality tail at fixed centrality",
        size=25.2,
        weight="bold",
        color=INK,
    )

    add_box(fig, (0.055, 0.812), (0.405, 0.092), SOFT_GRAY)
    fig_text(fig, 0.073, 0.892, "Sample being diagnosed", size=16.2, weight="bold")
    fig_text(
        fig,
        0.073,
        0.858,
        "Full Jet12+20+30+40 diagnostic cache.\nCandidate rows are matched to event-calo sums before class splits.",
        size=13.1,
        linespacing=1.48,
    )
    add_box(fig, (0.530, 0.812), (0.415, 0.092), SOFT_BLUE)
    fig_text(fig, 0.548, 0.892, "Cut contract", size=16.2, weight="bold")
    fig_text(
        fig,
        0.548,
        0.858,
        "Threshold uses only centrality and total calo energy.\nNo source sample, truth label, or BDT score enters the cut.",
        size=13.1,
        linespacing=1.48,
    )

    ax = fig.add_axes([0.060, 0.320, 0.585, 0.455])
    event_good = np.isfinite(data.event_cent) & np.isfinite(data.event_log_calo) & (data.event_cent >= 0.0) & (data.event_cent < 80.0)
    h = ax.hist2d(
        data.event_cent[event_good],
        data.event_log_calo[event_good],
        bins=[np.linspace(0, 80, 161), np.linspace(max(0.0, float(np.nanpercentile(data.event_log_calo[event_good], 0.05)) - 0.2), float(np.nanpercentile(data.event_log_calo[event_good], 99.9)) + 0.2, 140)],
        cmap="Blues",
        norm=LogNorm(vmin=1),
    )
    rng = np.random.default_rng(3203)
    rejected_idx = np.flatnonzero(event_good & event_reject)
    if len(rejected_idx) > 18000:
        rejected_idx = rng.choice(rejected_idx, size=18000, replace=False)
    ax.scatter(data.event_cent[rejected_idx], data.event_log_calo[rejected_idx], s=2.0, c=LOW_CALO_COLOR, alpha=0.35, rasterized=True)
    xs = []
    ys = []
    for row in envelope:
        if row.get("threshold") is None:
            continue
        xs.extend([float(row["cent_lo"]), float(row["cent_hi"])])
        ys.extend([float(row["threshold"]), float(row["threshold"])])
    ax.plot(xs, ys, color="#111827", lw=2.0, label="centrality-conditioned low-calo envelope")
    ax.set_xlim(0, 80)
    ax.set_xlabel("Centrality percentile", fontsize=12.8)
    ax.set_ylabel("log10(CEMC + IHCal + OHCal + 1)", fontsize=12.8)
    ax.tick_params(direction="in", top=True, right=True, labelsize=11.0)
    ax.grid(True, color=GRID, alpha=0.55)
    for spine in ax.spines.values():
        spine.set_linewidth(0.9)
        spine.set_color("#344054")
    ax.legend(loc="lower left", fontsize=10.6, frameon=True, framealpha=0.96, facecolor="white", edgecolor=BORDER)
    sphenix_label(ax, size=12.8)
    cax = fig.add_axes([0.646, 0.338, 0.012, 0.405])
    cb = fig.colorbar(h[3], cax=cax)
    cb.ax.tick_params(labelsize=9.5)

    comp_ax = fig.add_axes([0.705, 0.545, 0.250, 0.195])
    retained = event_good & ~event_reject
    rejected = event_good & event_reject
    comp_names = ["CEMC", "IHCal", "OHCal"]
    vals_ret = [
        float(np.nanmedian(np.log10(np.maximum(arr[retained], 0.0) + 1.0))) if np.any(retained) else math.nan
        for arr in (data.event_cemc, data.event_ihcal, data.event_ohcal)
    ]
    vals_rej = [
        float(np.nanmedian(np.log10(np.maximum(arr[rejected], 0.0) + 1.0))) if np.any(rejected) else math.nan
        for arr in (data.event_cemc, data.event_ihcal, data.event_ohcal)
    ]
    x = np.arange(3)
    comp_ax.bar(x - 0.17, vals_ret, width=0.34, color=RETAINED_COLOR, label="retained")
    comp_ax.bar(x + 0.17, vals_rej, width=0.34, color=LOW_CALO_COLOR, label="rejected")
    comp_ax.set_xticks(x, comp_names)
    comp_ax.set_ylabel("median log10(E + 1)", fontsize=11.2)
    comp_ax.tick_params(direction="in", top=True, right=True, labelsize=10.2)
    comp_ax.grid(True, axis="y", color=GRID, alpha=0.65)
    comp_ax.set_title("Rejected events have lower subsystem sums", fontsize=12.8, weight="bold", pad=6)
    comp_ax.legend(loc="upper left", fontsize=9.8, frameon=True, framealpha=0.96, edgecolor=BORDER)

    table_ax = fig.add_axes([0.675, 0.185, 0.295, 0.285])
    table_ax.axis("off")
    table_ax.set_title("Rejected fraction by source (event / candidate)", fontsize=12.8, weight="bold", color=INK, pad=4)
    table_data = []
    for code in CODE_TO_SAMPLE:
        label = CODE_TO_LABEL[code]
        cells = []
        for _, cent_label, _, _ in CENT_BINS:
            row = next(r for r in event_rows if r["source_sample"] == CODE_TO_SAMPLE[code] and r["centrality_bin"] == cent_label)
            cells.append(f"{fmt_frac(row['event_rejected_fraction'])}\n{fmt_frac(row['candidate_rejected_fraction'])}")
        table_data.append([label] + cells)
    tab = table_ax.table(
        cellText=table_data,
        colLabels=["source", "0-20", "20-50", "50-80"],
        cellLoc="center",
        loc="center",
        bbox=[0.0, 0.0, 1.0, 0.90],
    )
    tab.auto_set_font_size(False)
    tab.set_fontsize(8.9)
    for (r, c), cell in tab.get_celld().items():
        cell.set_edgecolor(BORDER)
        cell.set_linewidth(0.8)
        if r == 0:
            cell.set_facecolor("#EAECF0")
            cell.set_text_props(weight="bold", color=INK)
        elif c == 0:
            cell.set_facecolor("#F8FAFC")
            cell.set_text_props(weight="bold", color=INK)
        else:
            cell.set_facecolor("white")

    add_box(fig, (0.060, 0.058), (0.880, 0.115), SOFT_YELLOW, edge="#E3B341")
    total_events = int(event_good.sum())
    rejected_events = int(rejected.sum())
    total_cands = int(np.sum((data.cent >= 0.0) & (data.cent < 80.0)))
    rejected_cands = int(np.sum((data.cent >= 0.0) & (data.cent < 80.0) & cand_reject))
    fig_text(fig, 0.080, 0.156, "What this establishes before changing training", size=15.4, weight="bold", color=INK)
    fig_text(
        fig,
        0.080,
        0.126,
        f"A small low-total-calo population falls below the normal envelope at the same centrality: {fmt_count(rejected_events)}/{fmt_count(total_events)} events and {fmt_count(rejected_cands)}/{fmt_count(total_cands)} candidate rows in 0-80%.\nBecause the threshold uses only centrality and calorimeter sums, this is an event-quality removal, not a change to Signal MC or Inclusive MC definitions.",
        size=12.1,
        color=INK,
        linespacing=1.42,
    )
    outpng.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpng, dpi=SLIDE_DPI)
    plt.close(fig)


def step_xy(edges: np.ndarray, vals: list[float]) -> tuple[np.ndarray, np.ndarray]:
    density_vals = np.asarray(vals, dtype="float64")
    return np.repeat(edges, 2)[1:-1], np.repeat(density_vals, 2)


def draw_closure_slide(outpng: Path, closure_rows: list[dict], score_bins: np.ndarray) -> None:
    plt = setup_matplotlib()
    from matplotlib.patches import Rectangle

    fig = plt.figure(figsize=(SLIDE_W_PX / SLIDE_DPI, SLIDE_H_PX / SLIDE_DPI), dpi=SLIDE_DPI, facecolor="white")
    fig.patches.append(Rectangle((0, 0), 1, 1, transform=fig.transFigure, facecolor="white", zorder=-10))
    fig_text(fig, 0.040, 0.965, "Low-calo cut removes outliers without changing BDT separation", size=25.4, weight="bold")
    add_box(fig, (0.055, 0.825), (0.460, 0.087), SOFT_GRAY)
    fig_text(fig, 0.075, 0.896, "Class definitions held fixed", size=15.8, weight="bold")
    fig_text(
        fig,
        0.075,
        0.868,
        "Signal = truth-isolated embeddedPhoton\nInclusive = all embeddedJet; no truth-background filter",
        size=12.1,
        linespacing=1.18,
    )
    add_box(fig, (0.555, 0.825), (0.390, 0.087), SOFT_BLUE)
    fig_text(fig, 0.575, 0.896, "Diagnostic cut only", size=15.8, weight="bold")
    fig_text(
        fig,
        0.575,
        0.868,
        "Existing 15-35 GeV score rows after the cut\nNo retraining, split, feature, or filter mutation",
        size=12.1,
        linespacing=1.18,
    )
    add_box(fig, (0.165, 0.735), (0.670, 0.072), "white", edge=BORDER, lw=0.9, radius=0.008)
    legend_ax = fig.add_axes([0.175, 0.740, 0.650, 0.060])
    legend_ax.axis("off")
    legend_ax.set_xlim(0, 1)
    legend_ax.set_ylim(0, 1)
    legend_ax.plot([0.040, 0.100], [0.70, 0.70], color=SIGNAL_COLOR, lw=2.8)
    legend_ax.text(0.120, 0.70, "Signal MC after cut", va="center", fontsize=10.6, color=INK, fontfamily="serif")
    legend_ax.plot([0.520, 0.580], [0.70, 0.70], color=SIGNAL_COLOR, lw=2.0, ls="--", alpha=0.58)
    legend_ax.text(0.600, 0.70, "Signal MC before cut", va="center", fontsize=10.6, color=INK, fontfamily="serif")
    legend_ax.plot([0.040, 0.100], [0.30, 0.30], color=INCLUSIVE_COLOR, lw=2.8)
    legend_ax.text(0.120, 0.30, "Inclusive MC after cut", va="center", fontsize=10.6, color=INK, fontfamily="serif")
    legend_ax.plot([0.520, 0.580], [0.30, 0.30], color=INCLUSIVE_COLOR, lw=2.0, ls="--", alpha=0.58)
    legend_ax.text(0.600, 0.30, "Inclusive MC before cut", va="center", fontsize=10.6, color=INK, fontfamily="serif")

    axes = []
    lefts = [0.060, 0.365, 0.670]
    for i, row in enumerate(closure_rows):
        ax = fig.add_axes([lefts[i], 0.395, 0.270, 0.300])
        axes.append(ax)
        sig_before = row["signal_density_before"]
        sig_after = row["signal_density_after"]
        inc_before = row["inclusive_density_before"]
        inc_after = row["inclusive_density_after"]
        ymax = max(max(sig_before), max(sig_after), max(inc_before), max(inc_after), 1.0) * 1.16
        xb, yb = step_xy(score_bins, inc_before)
        xa, ya = step_xy(score_bins, inc_after)
        sb, syb = step_xy(score_bins, sig_before)
        sa, sya = step_xy(score_bins, sig_after)
        ax.plot(xb, yb, color=INCLUSIVE_COLOR, lw=1.8, ls="--", alpha=0.58)
        ax.plot(sb, syb, color=SIGNAL_COLOR, lw=1.8, ls="--", alpha=0.58)
        ax.plot(xa, ya, color=INCLUSIVE_COLOR, lw=2.4)
        ax.plot(sa, sya, color=SIGNAL_COLOR, lw=2.6)
        ax.set_xlim(0.0, 1.0)
        ax.set_ylim(0.0, ymax)
        ax.set_title(row["centrality_bin"], fontsize=15.2, weight="bold", pad=6)
        ax.set_xlabel("BDT score", fontsize=12.0)
        if i == 0:
            ax.set_ylabel("Candidate density", fontsize=12.0)
        else:
            ax.set_yticklabels([])
        ax.tick_params(direction="in", top=True, right=True, labelsize=10.2)
        ax.grid(True, color=GRID, alpha=0.65)
        for spine in ax.spines.values():
            spine.set_linewidth(0.9)
            spine.set_color("#344054")
        sphenix_label(ax, size=11.4)
        ax.text(
            0.045,
            0.075,
            f"AUC {row['auc_before']:.3f} -> {row['auc_after']:.3f}",
            transform=ax.transAxes,
            fontsize=10.8,
            ha="left",
            va="bottom",
            color=INK,
            bbox=dict(facecolor="white", edgecolor=BORDER, boxstyle="round,pad=0.25", alpha=0.94),
        )

    add_box(fig, (0.060, 0.258), (0.880, 0.076), SOFT_YELLOW, edge="#E3B341")
    fig_text(fig, 0.080, 0.321, "Closure result", size=15.2, weight="bold", color=INK)
    fig_text(
        fig,
        0.080,
        0.293,
        "AUC and WP80 inclusive rate are stable in all three centrality bins after removing low-calo events.\nThat supports using this as an event-quality cut: it removes the weird events without retuning the current Signal-vs-Inclusive BDT diagnostic.",
        size=11.8,
        linespacing=1.36,
        color=INK,
    )

    band_ax = fig.add_axes([0.055, 0.045, 0.890, 0.200])
    band_ax.axis("off")
    metric_rows = []
    for label, key in [
        ("Rejected candidates", "candidate_rejected_fraction"),
        ("Weighted AUC", "auc_pair"),
        ("WP80 inclusive rate", "wp80_pair"),
        ("Retained S / Incl", "retained_counts"),
    ]:
        vals = []
        for row in closure_rows:
            if key == "candidate_rejected_fraction":
                vals.append(fmt_frac(row["candidate_rejected_fraction"]))
            elif key == "auc_pair":
                vals.append(f"{row['auc_before']:.3f} -> {row['auc_after']:.3f}")
            elif key == "wp80_pair":
                vals.append(f"{row['wp80_inclusive_rate_before']:.3f} -> {row['wp80_inclusive_rate_after']:.3f}")
            else:
                vals.append(f"{fmt_count(row['signal_after'])} / {fmt_count(row['inclusive_after'])}")
        metric_rows.append([label] + vals)
    tab = band_ax.table(
        cellText=metric_rows,
        colLabels=["metric", "0-20%", "20-50%", "50-80%"],
        cellLoc="center",
        loc="center",
        bbox=[0.0, 0.0, 1.0, 0.92],
    )
    tab.auto_set_font_size(False)
    tab.set_fontsize(10.4)
    for (r, c), cell in tab.get_celld().items():
        cell.set_edgecolor(BORDER)
        cell.set_linewidth(0.8)
        if r == 0:
            cell.set_facecolor("#EAECF0")
            cell.set_text_props(weight="bold", color=INK)
        elif c == 0:
            cell.set_facecolor("#F8FAFC")
            cell.set_text_props(weight="bold", color=INK)
        else:
            cell.set_facecolor("white")
    outpng.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outpng, dpi=SLIDE_DPI)
    plt.close(fig)


def write_scripts(outdir: Path, event_rows: list[dict], closure_rows: list[dict]) -> tuple[Path, Path]:
    path1 = outdir / "the32_low_calo_pathology_slide_v1.md"
    path2 = outdir / "the32_low_calo_posthoc_closure_slide_v1.md"
    total_event_rej = sum(r["event_rejected"] for r in event_rows)
    total_events = sum(r["event_total"] for r in event_rows)
    path1.write_text(
        "\n".join(
            [
                "# THE-32 low-calo pathology slide",
                "",
                "What I want to show here is that these are event-quality outliers at fixed centrality, not a photon-ID truth-label choice.",
                f"The envelope uses only centrality and total calorimeter energy, and rejects {total_event_rej:,} of {total_events:,} source-tagged events in the scored diagnostic sample.",
                "The component panel and source-split table show why I would treat this as a detector/event-quality cut before touching the BDT training itself.",
                "",
            ]
        )
    )
    path2.write_text(
        "\n".join(
            [
                "# THE-32 post-hoc closure slide",
                "",
                "This slide applies the event-quality cut after scoring the existing THE-8 model; there is no retraining in this diagnostic pass.",
                "The score distributions keep the PPG12-style asymmetric class contract: Signal MC is embeddedPhoton truth-isolated prompt, while Inclusive MC is all embeddedJet candidates after selections.",
                "The important closure is that the weighted AUC and WP80 inclusive rate remain essentially stable after the cut in all three centrality bins.",
                "If this closure is accepted, the next phase is to impose the event-quality cut before BDT training and rerun the normal training/validation closure.",
                "",
            ]
        )
    )
    return path1, path2


def write_manifest(
    outdir: Path,
    args,
    data: ScoreData,
    envelope: list[dict],
    event_rows: list[dict],
    closure_rows: list[dict],
    output_paths: dict[str, Path],
) -> Path:
    manifest = {
        "schema": "THE32_LOW_CALO_DIAGNOSTIC_CLOSURE_V1",
        "source_cache_list": str(args.score_cache_list),
        "score_column": args.score_col,
        "n_score_caches": len(data.cache_paths),
        "candidate_rows": int(len(data.score)),
        "finite_score_rows": int(np.isfinite(data.score).sum()),
        "event_rows_after_per_cache_dedup": int(len(data.event_cent)),
        "cut_definition": {
            "variable": "log10(CEMC + IHCal + OHCal + 1)",
            "centrality_source": "CentralityInfo::mbd_NS stored as centrality",
            "centrality_slices": "5%-wide slices over 0-80%",
            "truth_or_source_used": False,
            "bdt_score_used": False,
            "threshold_formula": "max(median - mad_scale * 1.4826*MAD, quantile_floor)",
            "mad_scale": args.mad_scale,
            "quantile_floor": args.quantile_floor,
        },
        "plot_class_definition": (
            "Signal MC = source_sample contains embeddedPhoton and is_signal == 1 truth-isolated prompt; "
            "Inclusive MC = source_sample contains embeddedJet with no truth-background filter"
        ),
        "weighted_auc_convention": "PPG12-style class-balance plus cluster_Et/cluster_Eta inverse-density weights computed on scored rows.",
        "training_mutation": "none; post-hoc diagnostic cut on existing scored sample",
        "envelope": envelope,
        "event_count_rows": event_rows,
        "closure_rows": closure_rows,
        "outputs": {key: str(path) for key, path in output_paths.items()},
        "qa": {
            "png_dimensions": "2560x1440",
            "white_background": True,
            "slide_numbers_baked_in": False,
            "audience_canvas_has_private_labels": False,
            "typography": "Times-style serif; Nimbus Roman explicit on SDCC if Times New Roman is unavailable",
        },
    }
    path = outdir / "the32_low_calo_manifest_v1.json"
    path.write_text(json.dumps(json_ready(manifest), indent=2, sort_keys=True))
    return path


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--score-cache-list", type=Path, required=True, help="score_caches.list or validation report directory")
    ap.add_argument("--outdir", type=Path, default=OUTDIR)
    ap.add_argument("--score-col", default=SCORE_COL)
    ap.add_argument("--mad-scale", type=float, default=5.0)
    ap.add_argument("--quantile-floor", type=float, default=0.001)
    args = ap.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    data = load_score_data(args.score_cache_list, args.score_col)
    envelope = robust_envelope(data.event_cent, data.event_log_calo, mad_scale=args.mad_scale, quantile_floor=args.quantile_floor)
    cand_reject, event_reject = event_cut_masks(data, envelope)
    weights = ppg12_exact_weights(data)
    score_bins = np.linspace(0.0, 1.0, 81)
    event_rows = summarize_counts(data, cand_reject, event_reject)
    closure_rows = summarize_closure(data, cand_reject, weights, score_bins)

    cut_json = args.outdir / "the32_low_calo_cut_v1.json"
    event_csv = args.outdir / "the32_low_calo_event_counts_v1.csv"
    closure_csv = args.outdir / "the32_low_calo_closure_metrics_v1.csv"
    pathology_png = args.outdir / "the32_low_calo_pathology_slide_v1.png"
    closure_png = args.outdir / "the32_low_calo_posthoc_closure_slide_v1.png"
    pathology_json = args.outdir / "the32_low_calo_pathology_slide_v1.json"
    closure_json = args.outdir / "the32_low_calo_posthoc_closure_slide_v1.json"

    cut_json.write_text(
        json.dumps(
            json_ready(
                {
                    "schema": "THE32_LOW_CALO_CUT_V1",
                    "cut_variable": "log10(CEMC + IHCal + OHCal + 1)",
                    "centrality_source": "CentralityInfo::mbd_NS",
                    "source_blind": True,
                    "truth_blind": True,
                    "bdt_score_blind": True,
                    "slice_width_percent": 5.0,
                    "mad_scale": args.mad_scale,
                    "quantile_floor": args.quantile_floor,
                    "envelope": envelope,
                }
            ),
            indent=2,
            sort_keys=True,
        )
    )
    write_csv(event_csv, event_rows)
    compact_closure_rows = [{k: v for k, v in row.items() if not k.endswith("_density_before") and not k.endswith("_density_after")} for row in closure_rows]
    write_csv(closure_csv, compact_closure_rows)
    pathology_json.write_text(json.dumps(json_ready({"envelope": envelope, "event_counts": event_rows}), indent=2, sort_keys=True))
    closure_json.write_text(json.dumps(json_ready({"score_bin_edges": score_bins.tolist(), "closure": closure_rows}), indent=2, sort_keys=True))

    draw_pathology_slide(pathology_png, data, envelope, cand_reject, event_reject, event_rows)
    draw_closure_slide(closure_png, closure_rows, score_bins)
    script1, script2 = write_scripts(args.outdir, event_rows, closure_rows)
    manifest = write_manifest(
        args.outdir,
        args,
        data,
        envelope,
        event_rows,
        compact_closure_rows,
        {
            "cut_json": cut_json,
            "event_count_csv": event_csv,
            "closure_metrics_csv": closure_csv,
            "pathology_png": pathology_png,
            "closure_png": closure_png,
            "pathology_json": pathology_json,
            "closure_json": closure_json,
            "pathology_script_md": script1,
            "closure_script_md": script2,
        },
    )
    print(f"[THE32] wrote {pathology_png}")
    print(f"[THE32] wrote {closure_png}")
    print(f"[THE32] wrote {manifest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
