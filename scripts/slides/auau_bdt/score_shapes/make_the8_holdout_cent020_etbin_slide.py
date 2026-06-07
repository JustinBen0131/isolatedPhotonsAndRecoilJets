#!/usr/bin/env python3
"""Build THE-8 own-holdout 0-20% centrality ET-bin score-shape slide."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np


BRANCHES = [
    ("jet12_20", "Jet12+20"),
    ("jet12_20_30", "Jet12+20+30"),
    ("jet12_20_30_40", "Jet12+20+30+40"),
]
ET_BINS = [
    ("15-35", 15.0, 35.0),
    ("15-20", 15.0, 20.0),
    ("20-25", 20.0, 25.0),
    ("25-30", 25.0, 30.0),
    ("30-35", 30.0, 35.0),
]
SAMPLE_ACCENTS = (
    {"fill": "#ECF7EE", "edge": "#6EA77B", "strip": "#3F8F5A"},
    {"fill": "#FFF4CC", "edge": "#D6A84A", "strip": "#B88400"},
    {"fill": "#F3ECFA", "edge": "#9A7CC3", "strip": "#7651A6"},
)
PRODUCT = "globalEtCent1535_bdt_noIso"
SCORE_COL = f"score_{PRODUCT}"
WEIGHT_COL = "__ppg12_exact_training_weight"
PLOT_CLASS_DEFINITION = (
    "Signal MC = source_sample contains embeddedPhoton and is_signal == 1; "
    "Signal MC label is truth-isolated prompt; "
    "Inclusive MC = source_sample contains embeddedJet with no truth-background filter"
)
DEFAULT_SOURCE_BASE = Path("/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining")
DEFAULT_MODEL_BASE = Path("/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/bdt_models")
DEFAULT_OUTDIR = Path("dataOutput/auauTightBDTValidation/THE8_branchA_ladder_scorecache_fullstat_20260527/slide23_candidates")


def safe_ratio(frame: dict[str, np.ndarray], num: str, den: str) -> np.ndarray:
    n = frame[num].astype("float64", copy=False)
    d = frame[den].astype("float64", copy=False)
    out = np.full(len(n), np.nan, dtype="float64")
    good = np.isfinite(n) & np.isfinite(d) & (np.abs(d) > 1.0e-9)
    out[good] = n[good] / d[good]
    return out


def add_needed_derived(frame: dict[str, np.ndarray]) -> None:
    if "cluster_weta_over_wphi" not in frame:
        frame["cluster_weta_over_wphi"] = safe_ratio(frame, "cluster_weta_cogx", "cluster_wphi_cogx")
    if "cluster_weta33_over_wphi33" not in frame:
        frame["cluster_weta33_over_wphi33"] = safe_ratio(frame, "cluster_weta33_cogx", "cluster_wphi33_cogx")


def load_registry(path: Path) -> dict[str, object]:
    data = json.loads(path.read_text())
    model = data["models"][0]
    report = model["report"]
    return {
        "features": list(model["features"]),
        "pt_range": model.get("pt_range") or report.get("pt_range"),
        "cent_range": model.get("cent_range") or report.get("cent_range"),
        "reported_holdout_rows": int(report["overfit_diagnostics"]["holdout_rows"]),
        "reported_train_rows": int(report["overfit_diagnostics"]["train_rows"]),
        "reported_holdout_auc": float(report["overfit_diagnostics"]["holdout_auc"]),
    }


def load_matrix(path: Path, features: list[str]) -> dict[str, np.ndarray]:
    z = np.load(path, allow_pickle=True)
    needed = set(features) | {"is_signal", "centrality", "cluster_Et", "source_sample", WEIGHT_COL}
    base_needed = set()
    for name in needed:
        if name == "cluster_weta_over_wphi":
            base_needed.update(["cluster_weta_cogx", "cluster_wphi_cogx"])
        elif name == "cluster_weta33_over_wphi33":
            base_needed.update(["cluster_weta33_cogx", "cluster_wphi33_cogx"])
        else:
            base_needed.add(name)
    missing = sorted(col for col in base_needed if col not in z.files)
    if missing:
        raise SystemExit(f"{path} is missing columns: {missing}")
    frame = {name: z[name] for name in base_needed}
    add_needed_derived(frame)
    return frame


def filtered_indices(frame: dict[str, np.ndarray], features: list[str], pt_range, cent_range) -> np.ndarray:
    mask = np.ones(len(frame["is_signal"]), dtype=bool)
    if pt_range is not None:
        lo, hi = float(pt_range[0]), float(pt_range[1])
        et = frame["cluster_Et"]
        mask &= np.isfinite(et) & (et >= lo) & (et < hi)
    if cent_range is not None:
        lo, hi = float(cent_range[0]), float(cent_range[1])
        cent = frame["centrality"]
        mask &= np.isfinite(cent) & (cent >= lo) & (cent < hi)
    y = frame["is_signal"].astype("int32", copy=False)
    mask &= np.isin(y, [0, 1])
    for name in features:
        mask &= np.isfinite(frame[name])
    return np.flatnonzero(mask)


def own_holdout_indices(frame: dict[str, np.ndarray], registry: dict[str, object], seed: int) -> np.ndarray:
    from sklearn.model_selection import train_test_split

    idx = filtered_indices(frame, registry["features"], registry["pt_range"], registry["cent_range"])
    y = frame["is_signal"][idx].astype("int32", copy=False)
    _, test_idx = train_test_split(idx, test_size=0.10, random_state=seed, stratify=y)
    if len(test_idx) != registry["reported_holdout_rows"]:
        raise SystemExit(f"Holdout mismatch: got {len(test_idx)}, expected {registry['reported_holdout_rows']}")
    return np.asarray(test_idx, dtype="int64")


def density_hist(values: np.ndarray, bins: np.ndarray) -> dict[str, object]:
    counts, _ = np.histogram(values[np.isfinite(values)], bins=bins)
    width = np.diff(bins)
    density = counts.astype("float64")
    if density.sum() > 0:
        density = density / density.sum() / width
    return {"counts": counts.astype(int).tolist(), "density": density.astype(float).tolist(), "entries": int(counts.sum())}


def binned_auc(signal_counts: np.ndarray, background_counts: np.ndarray) -> float:
    sig_total = float(np.sum(signal_counts))
    bkg_total = float(np.sum(background_counts))
    if sig_total <= 0.0 or bkg_total <= 0.0:
        return math.nan
    bkg_below = np.cumsum(background_counts) - background_counts
    wins = float(np.sum(signal_counts * bkg_below))
    ties = float(np.sum(signal_counts * background_counts))
    return (wins + 0.5 * ties) / (sig_total * bkg_total)


def cache_paths(report_dir: Path) -> list[Path]:
    manifest = report_dir / "score_caches.list"
    if manifest.exists():
        return [Path(line.strip()) for line in manifest.read_text().splitlines() if line.strip()]
    return sorted((report_dir / "score_caches").glob("score_cache_*.npz"))


def extract_branch(tag: str, label: str, source_base: Path, model_base: Path, seed: int, score_bins: np.ndarray) -> dict[str, object]:
    model_dir = model_base / f"THE8_branchA_{tag}_global_noiso_bdt_mem24_20260527"
    source_dir = source_base / f"THE8_branchA_ladder_{tag}_20260527"
    report_dir = source_dir / "reports" / f"model_validation_condor_THE8_branchA_{tag}_scorecache_fullstat_20260527"
    registry = load_registry(model_dir / "model_registry.json")
    frame = load_matrix(model_dir / "training_matrix.npz", registry["features"])
    holdout = own_holdout_indices(frame, registry, seed)
    holdout_mask = np.zeros(len(frame["is_signal"]), dtype=bool)
    holdout_mask[holdout] = True

    score_parts = []
    y_parts = []
    source_parts = []
    cent_parts = []
    et_parts = []
    offset = 0
    scored_holdout_rows = 0
    for cache in cache_paths(report_dir):
        z = np.load(cache, allow_pickle=True)
        score = z[SCORE_COL].astype("float32", copy=False)
        n = len(score)
        stop = offset + n
        if stop > len(holdout_mask):
            raise SystemExit(f"Score cache exceeds matrix rows for {label}: {cache}")
        y_cache = z["is_signal"].astype("int32", copy=False)
        if not np.array_equal(y_cache, frame["is_signal"][offset:stop].astype("int32", copy=False)):
            raise SystemExit(f"Score-cache row order mismatch at {cache}")
        keep = holdout_mask[offset:stop] & np.isfinite(score)
        if np.any(keep):
            local = slice(offset, stop)
            score_parts.append(score[keep])
            y_parts.append(frame["is_signal"][local][keep].astype("int32", copy=False))
            source_parts.append(frame["source_sample"][local][keep].astype(str))
            cent_parts.append(frame["centrality"][local][keep].astype("float32", copy=False))
            et_parts.append(frame["cluster_Et"][local][keep].astype("float32", copy=False))
            scored_holdout_rows += int(np.sum(keep))
        offset = stop
    if offset != len(holdout_mask):
        raise SystemExit(f"Score caches do not cover matrix for {label}: cache={offset} matrix={len(holdout_mask)}")

    score = np.concatenate(score_parts)
    y = np.concatenate(y_parts)
    source = np.concatenate(source_parts)
    cent = np.concatenate(cent_parts)
    et = np.concatenate(et_parts)
    base = np.isfinite(score) & np.isfinite(cent) & np.isfinite(et) & (cent >= 0.0) & (cent < 20.0) & (et >= 15.0) & (et < 35.0)
    signal_base = base & (np.char.find(source, "embeddedPhoton") >= 0) & (y == 1)
    inclusive_base = base & (np.char.find(source, "embeddedJet") >= 0)

    et_payload = []
    for et_label, lo, hi in ET_BINS:
        emask = (et >= lo) & (et < hi)
        sig = signal_base & emask
        inc = inclusive_base & emask
        sig_hist = density_hist(score[sig], score_bins)
        inc_hist = density_hist(score[inc], score_bins)
        auc = binned_auc(np.asarray(sig_hist["counts"], dtype=float), np.asarray(inc_hist["counts"], dtype=float))
        et_payload.append(
            {
                "et_label": et_label,
                "et_lo": lo,
                "et_hi": hi,
                "auc_binned": auc,
                "signal": sig_hist,
                "inclusive": inc_hist,
            }
        )
        if sig_hist["entries"] <= 0 or inc_hist["entries"] <= 0:
            raise SystemExit(f"Empty class for {label} {et_label}: S={sig_hist['entries']} I={inc_hist['entries']}")

    return {
        "label": label,
        "model_dir": str(model_dir),
        "source": str(source_dir),
        "report_dir": str(report_dir),
        "model_registry": str(model_dir / "model_registry.json"),
        "reported_holdout_rows": registry["reported_holdout_rows"],
        "reported_train_rows": registry["reported_train_rows"],
        "scored_holdout_rows": scored_holdout_rows,
        "centrality": "0-20%",
        "et_bins": et_payload,
    }


def extract_payload(args: argparse.Namespace) -> dict[str, object]:
    score_bins = np.linspace(0.0, 1.0, args.score_bins + 1)
    return {
        "schema": "THE8_BRANCH_A_OWN_HOLDOUT_CENT020_ETBIN_SCORE_SHAPES_V1",
        "description": "Own 10% validation holdout score shapes in 0-20% centrality, using slide-7 truth-isolated signal vs inclusive-jet class contract.",
        "plot_class_definition": PLOT_CLASS_DEFINITION,
        "validation_scope": "Each row uses its own 10% training holdout; 0 <= centrality < 20; 15 <= cluster_Et < 35.",
        "columns": [{"label": label, "et_lo": lo, "et_hi": hi} for label, lo, hi in ET_BINS],
        "score_edges": score_bins.astype(float).tolist(),
        "branches": [extract_branch(tag, label, args.source_base, args.model_base, args.random_seed, score_bins) for tag, label in BRANCHES],
    }


def fmt_count(value: int) -> str:
    if value >= 1_000_000:
        return f"{value / 1_000_000:.2f}M"
    if value >= 100_000:
        return f"{value / 1_000:.0f}k"
    if value >= 1_000:
        return f"{value / 1_000:.1f}k"
    return str(value)


def add_card(
    fig,
    *,
    xy,
    wh,
    title,
    body,
    face,
    edge,
    title_color="#111827",
    title_size=16.6,
    body_size=13.6,
    body_linespacing=1.45,
    vertical="top",
    align="left",
):
    import matplotlib.patches as patches

    x, y = xy
    w, h = wh
    rect = patches.FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.010,rounding_size=0.010", transform=fig.transFigure, linewidth=1.05, edgecolor=edge, facecolor=face)
    fig.add_artist(rect)
    text_x = x + 0.5 * w if align == "center" else x + 0.016
    ha = "center" if align == "center" else "left"
    if vertical == "center":
        title_y = y + h - 0.024
        body_y = title_y - 0.033
    else:
        title_y = y + h - 0.021
        body_y = y + h - 0.052
    fig.text(text_x, title_y, title, ha=ha, va="top", fontsize=title_size, fontweight="bold", color=title_color)
    fig.text(text_x, body_y, body, ha=ha, va="top", fontsize=body_size, color="#334155", linespacing=body_linespacing)


def add_text_with_et_subscript(fig, *, x, y, prefix, suffix="", fontsize=16.0, fontweight="normal", color="#111827"):
    from matplotlib.font_manager import FontProperties
    from matplotlib.textpath import TextPath

    prop = FontProperties(family="serif", weight=fontweight, size=fontsize)
    fig_width_pts = fig.get_figwidth() * 72.0

    def width_frac(text: str) -> float:
        if not text:
            return 0.0
        return TextPath((0, 0), text, prop=prop).get_extents().width / fig_width_pts

    fig.text(x, y, prefix, ha="left", va="top", fontsize=fontsize, fontweight=fontweight, color=color)
    x_et = x + width_frac(prefix) + 0.004
    fig.text(x_et, y, "E", ha="left", va="top", fontsize=fontsize, fontweight=fontweight, color=color)
    x_sub = x_et + width_frac("E") + 0.001
    fig.text(x_sub, y - 0.014, "T", ha="left", va="top", fontsize=fontsize * 0.58, fontweight=fontweight, color=color)
    x_suffix = x_sub + width_frac("T") * 0.58 + 0.006
    fig.text(x_suffix, y, suffix, ha="left", va="top", fontsize=fontsize, fontweight=fontweight, color=color)


def draw_density(ax, edges: np.ndarray, density: np.ndarray, color: str, linestyle: str) -> None:
    y = np.r_[density, density[-1]]
    ax.fill_between(edges, y, step="post", color=color, alpha=0.10, linewidth=0)
    ax.step(edges, y, where="post", color=color, lw=1.95, linestyle=linestyle)


def formatted_sample_label(label: str) -> str:
    if label == "Jet12+20+30":
        return "Jet12+20\n+30"
    if label == "Jet12+20+30+40":
        return "Jet12+20\n+30+40"
    return label


def add_column_header(fig, *, xy, wh, label, ink, highlight=False):
    import matplotlib.patches as patches

    x, y = xy
    w, h = wh
    face = "#FDF2F8" if highlight else "#F8FAFC"
    edge = "#F9A8D4" if highlight else "#C9D2DE"
    rect = patches.FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.006,rounding_size=0.008",
        transform=fig.transFigure,
        linewidth=1.25 if highlight else 1.0,
        edgecolor=edge,
        facecolor=face,
    )
    fig.add_artist(rect)
    fig.text(
        x + 0.5 * w,
        y + 0.5 * h,
        f"{label} GeV",
        transform=fig.transFigure,
        ha="center",
        va="center",
        fontsize=12.6,
        fontweight="bold",
        color=ink,
    )


def add_inclusive_column_band(fig, axes, *, color="#FFF7FB", edge="#FBCFE8"):
    import matplotlib.patches as patches

    col_axes = [axes[i, 0] for i in range(axes.shape[0])]
    x0 = min(ax.get_position().x0 for ax in col_axes) - 0.010
    x1 = max(ax.get_position().x1 for ax in col_axes) + 0.010
    y0 = min(ax.get_position().y0 for ax in col_axes) - 0.030
    y1 = max(ax.get_position().y1 for ax in col_axes) + 0.062
    band = patches.FancyBboxPatch(
        (x0, y0),
        x1 - x0,
        y1 - y0,
        boxstyle="round,pad=0.004,rounding_size=0.010",
        transform=fig.transFigure,
        linewidth=1.0,
        edgecolor=edge,
        facecolor=color,
        alpha=0.62,
        zorder=-0.5,
    )
    fig.add_artist(band)


def add_row_header_card(fig, *, xy, wh, label, accent, ink):
    import matplotlib.patches as patches

    x, y = xy
    w, h = wh
    rect = patches.FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.010,rounding_size=0.008",
        transform=fig.transFigure,
        linewidth=1.1,
        edgecolor=accent["edge"],
        facecolor=accent["fill"],
    )
    fig.add_artist(rect)
    fig.add_artist(
        patches.Rectangle(
            (x + 0.004, y + 0.010),
            0.006,
            max(0.001, h - 0.020),
            transform=fig.transFigure,
            facecolor=accent["strip"],
            edgecolor="none",
        )
    )
    sample_label = formatted_sample_label(label)
    fig.text(
        x + 0.5 * w + 0.006,
        y + 0.5 * h + 0.028,
        "VALIDATION",
        transform=fig.transFigure,
        ha="center",
        va="center",
        fontsize=12.0,
        fontweight="bold",
        color=ink,
    )
    fig.text(
        x + 0.5 * w + 0.006,
        y + 0.5 * h - 0.008,
        sample_label,
        transform=fig.transFigure,
        ha="center",
        va="center",
        fontsize=13.1 if "\n" in sample_label else 13.8,
        fontweight="bold",
        color=ink,
        linespacing=1.18,
    )


def render_slide(payload: dict[str, object], outdir: Path, tag: str) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as patches

    outdir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "axes.linewidth": 0.95,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )
    branches = payload["branches"]
    edges = np.asarray(payload["score_edges"], dtype=float)
    sig_color = "#DC2626"
    inc_color = "#1D4ED8"
    ink = "#111827"
    muted = "#475569"
    grid = "#E5E7EB"
    ymax = 0.0
    for branch in branches:
        for item in branch["et_bins"]:
            ymax = max(ymax, float(np.max(item["signal"]["density"])), float(np.max(item["inclusive"]["density"])))
    ymax = 1.18 * ymax if ymax > 0 else 1.0

    fig, axes = plt.subplots(
        3,
        5,
        figsize=(16, 9),
        dpi=160,
        sharex=True,
        sharey=True,
        gridspec_kw={"left": 0.180, "right": 0.982, "bottom": 0.190, "top": 0.630, "hspace": 0.34, "wspace": 0.135},
    )
    fig.patch.set_facecolor("white")
    add_text_with_et_subscript(
        fig,
        x=0.035,
        y=0.967,
        prefix="0-20% validation holdout: score shapes by",
        suffix="bin",
        fontsize=27.0,
        fontweight="bold",
        color=ink,
    )
    add_card(
        fig,
        xy=(0.045, 0.802),
        wh=(0.590, 0.115),
        title="Class definition fixed; only energy bin changes",
        body="Signal MC: truth-isolated prompt embeddedPhoton rows.\nInclusive MC: all embeddedJet candidates, no truth-background filter.\nOwn 10% validation holdout; 0-20% centrality only.",
        face="#F8FAFC",
        edge="#CBD5E1",
        title_size=17.2,
        body_size=12.9,
        body_linespacing=1.16,
        vertical="center",
        align="left",
    )
    add_card(
        fig,
        xy=(0.670, 0.802),
        wh=(0.285, 0.115),
        title="Columns",
        body="15-35 is the inclusive validation window;\nremaining columns are 5 GeV sub-bins.",
        face="#EEF6FF",
        edge="#93C5FD",
        title_color="#1D4ED8",
        body_size=13.8,
    )
    legend = patches.FancyBboxPatch((0.210, 0.704), 0.590, 0.055, boxstyle="round,pad=0.006,rounding_size=0.010", transform=fig.transFigure, linewidth=1.05, edgecolor="#D1D5DB", facecolor="#FFFFFF")
    fig.add_artist(legend)
    yc = 0.7315
    fig.add_artist(matplotlib.lines.Line2D([0.270, 0.330], [yc, yc], transform=fig.transFigure, color=sig_color, lw=3.4))
    fig.text(0.340, yc, "Signal MC (truth-isolated prompt)", transform=fig.transFigure, ha="left", va="center", fontsize=12.9, color=ink)
    fig.add_artist(matplotlib.lines.Line2D([0.515, 0.575], [yc, yc], transform=fig.transFigure, color=inc_color, lw=3.4, linestyle="--"))
    fig.text(0.585, yc, "Inclusive MC (embedded jet; no truth filter)", transform=fig.transFigure, ha="left", va="center", fontsize=12.9, color=ink)

    add_inclusive_column_band(fig, axes)
    summary_rows = []
    for irow, branch in enumerate(branches):
        for icol, item in enumerate(branch["et_bins"]):
            ax = axes[irow, icol]
            if icol == 0:
                ax.set_facecolor("#FFFCFD")
            draw_density(ax, edges, np.asarray(item["signal"]["density"], dtype=float), sig_color, "-")
            draw_density(ax, edges, np.asarray(item["inclusive"]["density"], dtype=float), inc_color, "--")
            ax.set_xlim(0.0, 1.0)
            ax.set_ylim(0.0, ymax)
            ax.grid(True, color=grid, lw=0.50, alpha=0.85)
            ax.tick_params(labelsize=10.0, pad=2.5)
            for spine in ax.spines.values():
                spine.set_color(ink)
                spine.set_linewidth(0.95)
            if irow == 0:
                pos = ax.get_position()
                add_column_header(
                    fig,
                    xy=(pos.x0 + 0.004, pos.y1 + 0.014),
                    wh=(pos.width - 0.008, 0.034),
                    label=item["et_label"],
                    ink=ink,
                    highlight=(icol == 0),
                )
            if icol == 0:
                ax.set_ylabel("Unit-area\ndensity", fontsize=10.5, color=ink, labelpad=6)
                pos = ax.get_position()
                add_row_header_card(
                    fig,
                    xy=(0.030, pos.y0 + 0.010),
                    wh=(0.098, pos.height - 0.020),
                    label=branch["label"],
                    accent=SAMPLE_ACCENTS[irow % len(SAMPLE_ACCENTS)],
                    ink=ink,
                )
            if irow == 2:
                ax.set_xlabel("BDT score", fontsize=11.9, color=ink)
            if irow == 0 and icol == 0:
                ax.text(0.035, 0.900, r"$\it{\bf{sPHENIX}}$ Internal", transform=ax.transAxes, ha="left", va="top", fontsize=10.3, color=ink)
                ax.text(0.035, 0.750, "0-20% holdout", transform=ax.transAxes, ha="left", va="top", fontsize=9.4, color=ink)
            sig_entries = int(item["signal"]["entries"])
            inc_entries = int(item["inclusive"]["entries"])
            auc = float(item["auc_binned"])
            ax.text(0.965, 0.910, f"AUC {auc:.3f}", transform=ax.transAxes, ha="right", va="top", fontsize=12.2, fontweight="bold", color=ink, bbox=dict(boxstyle="round,pad=0.20", fc="white", ec="#CBD5E1", lw=0.75, alpha=0.96))
            ax.text(0.965, 0.675, f"S {fmt_count(sig_entries)}  I {fmt_count(inc_entries)}", transform=ax.transAxes, ha="right", va="top", fontsize=9.6, fontweight="bold", color=muted)
            summary_rows.append(
                {
                    "sample": branch["label"],
                    "et_label": item["et_label"],
                    "et_lo": item["et_lo"],
                    "et_hi": item["et_hi"],
                    "auc_binned": auc,
                    "signal_entries": sig_entries,
                    "inclusive_entries": inc_entries,
                }
            )

    add_card(
        fig,
        xy=(0.038, 0.036),
        wh=(0.924, 0.085),
        title="Read as a validation-sample shape table",
        body="Rows compare the three independently trained Branch A samples; columns compare the inclusive 15-35 GeV validation window and its four 5 GeV slices in central 0-20% Au+Au.",
        face="#FFF7ED",
        edge="#FDBA74",
        title_color="#9A3412",
        body_size=13.3,
    )

    png_path = outdir / f"{tag}.png"
    csv_path = outdir / f"{tag}.csv"
    json_path = outdir / f"{tag}.json"
    md_path = outdir / f"{tag}.md"
    fig.savefig(png_path)
    plt.close(fig)
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0]))
        writer.writeheader()
        writer.writerows(summary_rows)
    plot_class_definition = str(payload.get("plot_class_definition", ""))
    plot_class_definition = plot_class_definition.replace(
        "Signal MC = source_sample contains embeddedPhoton and is_signal == 1;",
        "Signal MC = source_sample contains embeddedPhoton and is_signal == 1 (truth-isolated prompt label);",
    )
    json_path.write_text(
        json.dumps(
            {
                "schema": "THE8_BRANCH_A_OWN_HOLDOUT_CENT020_ETBIN_SLIDE_V1",
                "input_schema": payload.get("schema", ""),
                "plot_class_definition": plot_class_definition,
                "validation_scope": payload.get("validation_scope", ""),
                "png": str(png_path),
                "summary_csv": str(csv_path),
                "rows": summary_rows,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    md_path.write_text(make_script_text() + "\n")
    print(png_path)
    print(csv_path)
    print(json_path)
    print(md_path)


def make_script_text() -> str:
    return """# WP GammaJets Slide Script - 0-20% Holdout ET-Bin Score Shapes

This slide is the same score-shape diagnostic as the slide-7 validation comparison, but now restricted to the most central 0-20 percent Au+Au bin and sliced by photon transverse energy.

The red curve is Signal MC: embedded-photon candidates that are truth-isolated prompt photons. The blue curve is Inclusive MC: all embedded-jet candidates, without applying a truth-background filter. So the legend keeps the same intentionally asymmetric contract as the slide-7 diagnostic.

Each row corresponds to one independently trained Branch A sample: Jet12 plus Jet20, then adding Jet30, then adding Jet40. Each row is evaluated only on its own ten percent validation holdout, separate from training.

The first column shows the full 15 to 35 GeV validation window. The next four columns split that same window into 15-20, 20-25, 25-30, and 30-35 GeV. This lets us see whether the apparent signal versus inclusive-jet separation is stable across energy inside the central Au+Au region.

The point of this slide is not to introduce a new validation definition. It is to make the 0-20 percent validation output more granular by energy, while keeping the signal-MC versus inclusive-MC meaning exactly consistent with the previous slide.
"""


def main() -> int:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)
    extract = sub.add_parser("extract")
    extract.add_argument("--source-base", type=Path, default=DEFAULT_SOURCE_BASE)
    extract.add_argument("--model-base", type=Path, default=DEFAULT_MODEL_BASE)
    extract.add_argument("--random-seed", type=int, default=13)
    extract.add_argument("--score-bins", type=int, default=50)
    render = sub.add_parser("render")
    render.add_argument("--input", type=Path, required=True)
    render.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    render.add_argument("--tag", default="the8_branchA_ladder_holdout_cent020_etbins_truthSignal_inclusiveJet_v1")
    args = parser.parse_args()
    if args.command == "extract":
        print(json.dumps(extract_payload(args), indent=2, sort_keys=True))
        return 0
    payload = json.loads(args.input.read_text())
    render_slide(payload, args.outdir, args.tag)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
