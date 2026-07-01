#!/usr/bin/env python3
"""Derive and render the Au+Au BDT bounded non-tight sideband slide.

The slide is intentionally score-space first: it uses the exact THE-57 default
Au+Au BDT validation cache and defines the non-tight sideband relative to the
centrality-dependent WP80 cut, T80(c).  This makes the sideband centrality aware
without inventing separate hard-coded windows for each centrality bin.
"""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch


REPO = Path(__file__).resolve().parents[3]
CACHE_DIR = REPO / "dataOutput/auauTightBDTValidation/THE57_baseline_validation_20260615/a1_withcut_score_caches/score_caches"
OUT_DIR = REPO / "dataOutput/the88_dijet_combinatoric_root_cause/sideband_derivation"
OUT_PNG = OUT_DIR / "slide02_auau_bdt_nontight_sideband_derivation.png"
OUT_MANIFEST = OUT_DIR / "slide02_auau_bdt_nontight_sideband_derivation_manifest.json"
OUT_CSV = OUT_DIR / "slide02_auau_bdt_nontight_sideband_scan_metrics.csv"
OUT_NOTES = OUT_DIR / "slide02_auau_bdt_nontight_sideband_derivation_speaker_script.md"

MODEL_ID = "centAsFeatBase3x3_pt15to35"
SCORE_KEY = "score_centAsFeatBase3x3_pt15to35"
T80_INTERCEPT = 0.53471108
T80_SLOPE = 0.0012284143
SIDEBAND_GAP = 0.03
SIDEBAND_WIDTH = 0.20
CENT_BINS = [(0, 20), (20, 50), (50, 80)]


@dataclass
class ScoreData:
    rel_score: np.ndarray
    centrality: np.ndarray
    label: np.ndarray
    weight: np.ndarray
    et: np.ndarray
    n_raw_rows: int


def t80(c: np.ndarray | float) -> np.ndarray | float:
    return T80_INTERCEPT + T80_SLOPE * c


def load_score_data() -> ScoreData:
    files = sorted(CACHE_DIR.glob("*.npz"))
    if not files:
        raise FileNotFoundError(f"No score cache files found in {CACHE_DIR}")

    rels: list[np.ndarray] = []
    cents: list[np.ndarray] = []
    labels: list[np.ndarray] = []
    weights: list[np.ndarray] = []
    ets: list[np.ndarray] = []
    n_raw = 0

    for path in files:
        with np.load(path, allow_pickle=True) as z:
            score = z[SCORE_KEY].astype(float)
            cent = z["centrality"].astype(float)
            et = z["cluster_Et"].astype(float)
            label = z["is_signal"].astype(np.int8)
            weight = z["event_weight"].astype(float) if "event_weight" in z.files else np.ones_like(score)
            n_raw += int(score.size)
            mask = (
                (et >= 15.0)
                & (et < 35.0)
                & (cent >= 0.0)
                & (cent < 80.0)
                & np.isfinite(score)
                & np.isfinite(cent)
                & np.isfinite(weight)
            )
            if not np.any(mask):
                continue
            rels.append(score[mask] - t80(cent[mask]))
            cents.append(cent[mask])
            labels.append(label[mask])
            weights.append(weight[mask])
            ets.append(et[mask])

    return ScoreData(
        rel_score=np.concatenate(rels),
        centrality=np.concatenate(cents),
        label=np.concatenate(labels),
        weight=np.concatenate(weights),
        et=np.concatenate(ets),
        n_raw_rows=n_raw,
    )


def ratio(num: float, den: float) -> float:
    return float(num / den) if den else float("nan")


def compute_metrics(data: ScoreData, gap: float, width: float) -> dict[str, dict[str, float]]:
    rel = data.rel_score
    lab = data.label
    cent = data.centrality
    w = data.weight

    tight = rel > 0.0
    complement = rel <= 0.0
    sideband = (rel > -width) & (rel < -gap)
    gap_region = (rel >= -gap) & (rel <= 0.0)
    far_tail = rel <= -width

    metrics: dict[str, dict[str, float]] = {}
    for lo, hi in CENT_BINS:
        key = f"{lo}-{hi}%"
        cm = (cent >= lo) & (cent < hi)
        sig = lab == 1
        bkg = lab == 0
        sig_t = float(w[cm & sig & tight].sum())
        sig_c = float(w[cm & sig & complement].sum())
        sig_s = float(w[cm & sig & sideband].sum())
        sig_gap = float(w[cm & sig & gap_region].sum())
        sig_far = float(w[cm & sig & far_tail].sum())
        bkg_t = float(w[cm & bkg & tight].sum())
        bkg_c = float(w[cm & bkg & complement].sum())
        bkg_s = float(w[cm & bkg & sideband].sum())
        bkg_gap = float(w[cm & bkg & gap_region].sum())
        bkg_far = float(w[cm & bkg & far_tail].sum())
        metrics[key] = {
            "signal_tight": sig_t,
            "signal_complement": sig_c,
            "signal_sideband": sig_s,
            "signal_gap": sig_gap,
            "signal_far_tail": sig_far,
            "background_tight": bkg_t,
            "background_complement": bkg_c,
            "background_sideband": bkg_s,
            "background_gap": bkg_gap,
            "background_far_tail": bkg_far,
            "signal_sideband_over_tight": ratio(sig_s, sig_t),
            "signal_complement_over_tight": ratio(sig_c, sig_t),
            "background_sideband_over_tight": ratio(bkg_s, bkg_t),
            "background_complement_over_tight": ratio(bkg_c, bkg_t),
            "sideband_rows": float(np.count_nonzero(cm & sideband)),
            "gap_rows": float(np.count_nonzero(cm & gap_region)),
            "far_tail_rows": float(np.count_nonzero(cm & far_tail)),
        }
    return metrics


def scan_windows(data: ScoreData) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    for gap in [0.02, 0.03, 0.04, 0.05, 0.06]:
        for width in [0.12, 0.16, 0.20, 0.24, 0.28, 0.32]:
            if width <= gap:
                continue
            metrics = compute_metrics(data, gap, width)
            leaks = [metrics[f"{lo}-{hi}%"]["signal_sideband_over_tight"] for lo, hi in CENT_BINS]
            bkgs = [metrics[f"{lo}-{hi}%"]["background_sideband_over_tight"] for lo, hi in CENT_BINS]
            side_rows = [metrics[f"{lo}-{hi}%"]["sideband_rows"] for lo, hi in CENT_BINS]
            row = {
                "gap": gap,
                "width": width,
                "signal_sideband_over_tight_mean": float(np.nanmean(leaks)),
                "signal_sideband_over_tight_spread": float(np.nanmax(leaks) - np.nanmin(leaks)),
                "background_sideband_over_tight_mean": float(np.nanmean(bkgs)),
                "background_sideband_over_tight_spread": float(np.nanmax(bkgs) - np.nanmin(bkgs)),
                "min_sideband_rows": float(np.min(side_rows)),
                "selected": float(abs(gap - SIDEBAND_GAP) < 1e-9 and abs(width - SIDEBAND_WIDTH) < 1e-9),
            }
            rows.append(row)
    return rows


def write_scan_csv(rows: list[dict[str, float]]) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0].keys())
    with OUT_CSV.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "dejavuserif",
            "axes.linewidth": 1.1,
            "axes.labelsize": 13,
            "axes.titlesize": 15,
            "xtick.labelsize": 11,
            "ytick.labelsize": 11,
            "legend.fontsize": 11,
        }
    )


def add_box(fig: plt.Figure, xywh: tuple[float, float, float, float], *, fc: str, ec: str = "#cbd5e1", lw: float = 1.2) -> None:
    box = FancyBboxPatch(
        (xywh[0], xywh[1]),
        xywh[2],
        xywh[3],
        boxstyle="round,pad=0.012,rounding_size=0.012",
        linewidth=lw,
        edgecolor=ec,
        facecolor=fc,
        transform=fig.transFigure,
        zorder=-1,
    )
    fig.add_artist(box)


def render_slide(data: ScoreData, metrics: dict[str, dict[str, float]], scan_rows: list[dict[str, float]]) -> None:
    setup_style()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    fig = plt.figure(figsize=(12.8, 7.2), dpi=200)
    fig.patch.set_facecolor("white")

    fig.text(
        0.045,
        0.947,
        "Deriving a bounded Au+Au BDT non-tight sideband",
        fontsize=24,
        fontweight="bold",
        color="#0f172a",
        va="top",
    )
    fig.text(0.055, 0.865, "signal MC", fontsize=11.2, color="#0b8f65", fontweight="bold")
    fig.text(0.135, 0.865, "inclusive MC", fontsize=11.2, color="#dc7f23", fontweight="bold")
    fig.text(0.250, 0.865, "blue band = proposed sideband", fontsize=11.2, color="#1d4ed8", fontweight="bold")
    fig.text(0.455, 0.865, "black line = tight T80(c)", fontsize=11.2, color="#111827", fontweight="bold")

    # Score distributions by centrality, using relative score to make the
    # centrality-dependent sideband visible in one common coordinate.
    panel_y = 0.535
    panel_h = 0.265
    x0 = 0.055
    panel_w = 0.176
    dx = 0.192
    bins = np.linspace(-0.55, 0.30, 70)
    for idx, (lo, hi) in enumerate(CENT_BINS):
        ax = fig.add_axes([x0 + idx * dx, panel_y, panel_w, panel_h])
        cm = (data.centrality >= lo) & (data.centrality < hi)
        sig = cm & (data.label == 1)
        bkg = cm & (data.label == 0)
        ax.axvspan(-SIDEBAND_WIDTH, -SIDEBAND_GAP, color="#dbeafe", alpha=0.95)
        ax.axvspan(-SIDEBAND_GAP, 0.0, color="#e5e7eb", alpha=0.65)
        ax.axvline(0.0, color="#111827", lw=1.6, ls="-")
        ax.axvline(-SIDEBAND_GAP, color="#2563eb", lw=1.2, ls="--")
        ax.axvline(-SIDEBAND_WIDTH, color="#2563eb", lw=1.2, ls="--")
        ax.hist(
            data.rel_score[bkg],
            bins=bins,
            density=True,
            histtype="step",
            lw=1.8,
            color="#dc7f23",
        )
        ax.hist(
            data.rel_score[sig],
            bins=bins,
            density=True,
            histtype="step",
            lw=1.8,
            color="#0b8f65",
        )
        ax.set_title(f"Au+Au {lo}-{hi}%", pad=6, fontweight="bold")
        ax.set_xlim(-0.55, 0.30)
        ax.set_ylim(bottom=0)
        ax.grid(True, color="#cbd5e1", alpha=0.35, lw=0.7)
        if idx == 0:
            ax.set_ylabel("normalized entries")
        else:
            ax.set_yticklabels([])
        ax.set_xlabel(r"score $-$ $T_{80}(c)$")

    # Centrality-dependent cut geometry.
    ax_cut = fig.add_axes([0.675, 0.520, 0.290, 0.295])
    c = np.linspace(0, 80, 240)
    upper = t80(c) - SIDEBAND_GAP
    lower = t80(c) - SIDEBAND_WIDTH
    ax_cut.fill_between(c, lower, upper, color="#dbeafe", alpha=0.9)
    ax_cut.plot(c, t80(c), color="#111827", lw=2.2)
    ax_cut.plot(c, upper, color="#2563eb", lw=1.8, ls="--")
    ax_cut.plot(c, lower, color="#2563eb", lw=1.8, ls=":")
    ax_cut.set_title("score thresholds vs centrality", fontweight="bold", pad=6)
    ax_cut.set_xlabel("centrality percentile c", labelpad=1)
    ax_cut.set_ylabel("BDT score")
    ax_cut.set_xlim(0, 80)
    ax_cut.set_ylim(0.30, 0.66)
    ax_cut.set_xticks([0, 20, 50, 80])
    ax_cut.grid(True, color="#cbd5e1", alpha=0.45, lw=0.7)
    label_box = dict(facecolor="white", edgecolor="none", alpha=0.82, pad=1.5)
    ax_cut.text(77, t80(77) + 0.006, "tight T80(c)", color="#111827", fontsize=10.5, ha="right", va="bottom", bbox=label_box)
    ax_cut.text(77, upper[-1] - 0.004, "sideband upper", color="#1d4ed8", fontsize=10.5, ha="right", va="top", bbox=label_box)
    ax_cut.text(77, lower[-1] + 0.004, "sideband lower", color="#1d4ed8", fontsize=10.5, ha="right", va="bottom", bbox=label_box)

    # Metrics table.
    ax_tbl = fig.add_axes([0.055, 0.205, 0.590, 0.275])
    ax_tbl.axis("off")
    rows = []
    for lo, hi in CENT_BINS:
        key = f"{lo}-{hi}%"
        m = metrics[key]
        rows.append(
            [
                key,
                f"{100*m['signal_sideband_over_tight']:.1f}%",
                f"{m['background_sideband_over_tight']:.2f}",
                f"{m['background_complement_over_tight']:.2f}",
                f"{int(m['sideband_rows']):,}",
            ]
        )
    col_labels = [
        "centrality\nbin",
        "prompt γ leakage\n(sideband / tight)",
        "bounded sideband\n(background / tight)",
        "old broad complement\n(background / tight)",
        "sideband\nentries",
    ]
    table = ax_tbl.table(
        cellText=rows,
        colLabels=col_labels,
        cellLoc="center",
        colLoc="center",
        loc="center",
        colWidths=[0.13, 0.24, 0.24, 0.24, 0.15],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(11.3)
    table.scale(1.00, 2.16)
    for (r, cidx), cell in table.get_celld().items():
        cell.set_edgecolor("#cbd5e1")
        cell.set_linewidth(0.9)
        if r == 0:
            cell.set_facecolor("#eef2ff")
            cell.set_text_props(fontweight="bold", color="#0f172a")
        elif cidx == 3:
            cell.set_facecolor("#fff1f2")
        elif cidx in (1, 2):
            cell.set_facecolor("#f0fdf4")
        else:
            cell.set_facecolor("white")

    # Right-side interpretation card.
    add_box(fig, (0.705, 0.215, 0.250, 0.218), fc="#f8fafc", ec="#cbd5e1")
    cx = 0.830
    fig.text(cx, 0.404, "Why this window", fontsize=15.4, fontweight="bold", color="#0f172a", ha="center")
    fig.text(cx, 0.363, "same centrality slope", fontsize=13.3, fontweight="bold", color="#334155", ha="center")
    fig.text(cx, 0.327, "0.03 gap avoids tight-boundary migration", fontsize=12.6, color="#334155", ha="center")
    fig.text(cx, 0.292, "0.20 width keeps sideband statistics", fontsize=12.6, color="#334155", ha="center")
    fig.text(cx, 0.257, "rejects the broad low-score complement", fontsize=12.6, color="#334155", ha="center")

    # Bottom decision band.  Scan details stay in the manifest/speaker notes so
    # the audience-facing slide stays focused on the decision.
    add_box(fig, (0.055, 0.055, 0.900, 0.120), fc="#eff6ff", ec="#bfdbfe", lw=1.4)
    fig.text(
        0.075,
        0.140,
        "Nominal non-tight sideband to validate",
        fontsize=14.3,
        fontweight="bold",
        color="#0f172a",
        va="center",
    )
    fig.text(
        0.075,
        0.096,
        r"$T_{80}(c)-0.20 < \mathrm{score} < T_{80}(c)-0.03$",
        fontsize=19,
        fontweight="bold",
        color="#1d4ed8",
        va="center",
    )
    fig.text(
        0.570,
        0.123,
        "Use this as the bounded Region C/D sideband.",
        fontsize=12.0,
        fontweight="bold",
        color="#334155",
        va="center",
    )
    fig.text(
        0.570,
        0.087,
        "Scores below the lower edge are neither tight nor sideband.",
        fontsize=11.4,
        color="#334155",
        va="center",
    )

    fig.savefig(OUT_PNG, dpi=200)
    plt.close(fig)


def write_manifest(data: ScoreData, metrics: dict[str, dict[str, float]], scan_rows: list[dict[str, float]]) -> None:
    payload = {
        "artifact": str(OUT_PNG),
        "model_id": MODEL_ID,
        "score_key": SCORE_KEY,
        "cache_dir": str(CACHE_DIR),
        "cache_files": len(list(CACHE_DIR.glob("*.npz"))),
        "raw_rows_seen": data.n_raw_rows,
        "selected_rows_15_35_cent_0_80": int(data.rel_score.size),
        "tight_wp80": {
            "mode": "centlinear",
            "formula": "T80(c) = 0.53471108 + 0.0012284143*c",
            "intercept": T80_INTERCEPT,
            "slope": T80_SLOPE,
        },
        "sideband_candidate": {
            "relative_score_variable": "score - T80(c)",
            "lower": "T80(c) - 0.20",
            "upper": "T80(c) - 0.03",
            "gap": SIDEBAND_GAP,
            "width": SIDEBAND_WIDTH,
            "rule": "T80(c)-0.20 < score < T80(c)-0.03",
            "below_lower_edge": "neither tight nor region C",
        },
        "metrics_by_centrality": metrics,
        "scan_csv": str(OUT_CSV),
        "validation_status": "score-level derivation only; needs bounded-sideband RecoilJets canary before default production freeze",
    }
    OUT_MANIFEST.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def write_speaker_notes(metrics: dict[str, dict[str, float]]) -> None:
    lines = [
        "# Au+Au BDT non-tight sideband derivation",
        "",
        "The sideband is defined relative to the default Au+Au WP80 cut, so centrality is handled by construction:",
        "`T80(c)-0.20 < score < T80(c)-0.03`.",
        "",
        "The 0.03 gap avoids using photons right at the tight boundary, where migration and prompt-photon leakage are most dangerous. "
        "The 0.20 width preserves statistics while removing the far below-threshold complement tail.",
        "",
        "The evidence on the slide is the sideband-to-tight behavior in the exact THE57 default-model validation score cache. "
        "The bounded window gives stable signal leakage across centrality while the old broad complement gives a very large inclusive-background control region.",
        "",
        "Key numbers:",
    ]
    for cent, vals in metrics.items():
        lines.append(
            f"- {cent}: signal side/tight={100*vals['signal_sideband_over_tight']:.1f}%, "
            f"bounded inclusive side/tight={vals['background_sideband_over_tight']:.2f}, "
            f"broad-complement inclusive side/tight={vals['background_complement_over_tight']:.2f}."
        )
    lines.extend(
        [
            "",
            "This does not yet make the sideband final. The next validation is a small RecoilJets canary using `nonTightAuAuBDTSideband`, followed by leakage/purity/xJ closure.",
        ]
    )
    OUT_NOTES.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    data = load_score_data()
    metrics = compute_metrics(data, SIDEBAND_GAP, SIDEBAND_WIDTH)
    scan_rows = scan_windows(data)
    write_scan_csv(scan_rows)
    render_slide(data, metrics, scan_rows)
    write_manifest(data, metrics, scan_rows)
    write_speaker_notes(metrics)
    print(OUT_PNG)
    print(OUT_MANIFEST)
    print(OUT_CSV)
    print(OUT_NOTES)


if __name__ == "__main__":
    main()
