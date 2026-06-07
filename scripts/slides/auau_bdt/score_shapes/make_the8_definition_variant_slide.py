#!/usr/bin/env python3
"""Build a THE-8 Jet12+20+30+40 score-definition diagnostic slide."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np


CENTRALITY = [("0_20", "0-20%", 0.0, 20.0), ("20_50", "20-50%", 20.0, 50.0), ("50_80", "50-80%", 50.0, 80.0)]
DEFINITION_ACCENTS = (
    {"fill": "#F5EFE7", "edge": "#C8A886", "strip": "#8B5E34"},
    {"fill": "#F2F3F0", "edge": "#A8AC9A", "strip": "#6B705C"},
    {"fill": "#F1F1F1", "edge": "#A7A7A7", "strip": "#4A4A4A"},
)
PRODUCT = "globalEtCent1535_bdt_noIso"
SCORE_COL = f"score_{PRODUCT}"
WEIGHT_COL = "__ppg12_exact_training_weight"
BRANCH_LABEL = "Jet12+20+30+40"
DEFAULT_MODEL_DIR = Path("/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/bdt_models/THE8_branchA_jet12_20_30_40_global_noiso_bdt_mem24_20260527")
DEFAULT_REPORT_DIR = Path("/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/THE8_branchA_ladder_jet12_20_30_40_20260527/reports/model_validation_condor_THE8_branchA_jet12_20_30_40_scorecache_fullstat_20260527")
DEFAULT_OUTDIR = Path("dataOutput/auauTightBDTValidation/THE8_branchA_ladder_scorecache_fullstat_20260527/slide23_candidates")


VARIANTS = [
    {
        "key": "current",
        "label": "PPG12 equivalent",
        "row_label": "PPG12 equivalent\nS = truth-isolated\nprompt\nB = all jet",
        "row_title": "PPG12 equivalent",
        "signal_label": "truth-isolated prompt",
        "background_label": "all jet",
        "signal_definition": "source_sample contains embeddedPhoton and is_signal == 1 (truth-isolated prompt)",
        "background_definition": "source_sample contains embeddedJet, no is_signal filter",
    },
    {
        "key": "truth_tagged",
        "label": "Truth-tagged classes",
        "row_label": "Truth-tagged\nS = truth-isolated\nprompt\nB = jet truth bkg",
        "row_title": "Truth-tagged",
        "signal_label": "truth-isolated prompt",
        "background_label": "jet truth bkg",
        "signal_definition": "source_sample contains embeddedPhoton and is_signal == 1 (truth-isolated prompt)",
        "background_definition": "source_sample contains embeddedJet and is_signal == 0",
    },
    {
        "key": "source_only",
        "label": "Source-only samples",
        "row_label": "Source-only\nS = all photon\nB = all jet",
        "row_title": "Source-only",
        "signal_label": "all photon",
        "background_label": "all jet",
        "signal_definition": "source_sample contains embeddedPhoton, no is_signal filter",
        "background_definition": "source_sample contains embeddedJet, no is_signal filter",
    },
]


def density_hist(values: np.ndarray, bins: np.ndarray) -> dict[str, object]:
    counts, _ = np.histogram(values[np.isfinite(values)], bins=bins)
    width = np.diff(bins)
    density = counts.astype("float64")
    if density.sum() > 0:
        density = density / density.sum() / width
    return {"counts": counts.astype(int).tolist(), "density": density.astype(float).tolist(), "entries": int(counts.sum())}


def weighted_hist(values: np.ndarray, weights: np.ndarray, bins: np.ndarray) -> np.ndarray:
    good = np.isfinite(values) & np.isfinite(weights) & (weights > 0.0)
    return np.histogram(values[good], bins=bins, weights=weights[good])[0].astype("float64")


def binned_auc_from_weighted_counts(signal_counts: np.ndarray, background_counts: np.ndarray) -> float:
    sig_total = float(np.sum(signal_counts))
    bkg_total = float(np.sum(background_counts))
    if sig_total <= 0.0 or bkg_total <= 0.0:
        return math.nan
    bkg_below = np.cumsum(background_counts) - background_counts
    wins = float(np.sum(signal_counts * bkg_below))
    ties = float(np.sum(signal_counts * background_counts))
    return (wins + 0.5 * ties) / (sig_total * bkg_total)
    return wins / (pos_total * neg_total)


def variant_masks(variant_key: str, is_photon: np.ndarray, is_jet: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    if variant_key == "current":
        return is_photon & (y == 1), is_jet
    if variant_key == "truth_tagged":
        return is_photon & (y == 1), is_jet & (y == 0)
    if variant_key == "source_only":
        return is_photon, is_jet
    raise ValueError(f"Unknown variant {variant_key}")


def extract_payload(model_dir: Path, report_dir: Path) -> dict[str, object]:
    matrix_path = model_dir / "training_matrix.npz"
    cache_list = report_dir / "score_caches.list"
    cache_paths = [Path(line.strip()) for line in cache_list.read_text().splitlines() if line.strip()] if cache_list.exists() else sorted((report_dir / "score_caches").glob("score_cache_*.npz"))
    if not cache_paths:
        raise SystemExit(f"No score caches found under {report_dir}")

    matrix = np.load(matrix_path, allow_pickle=True)
    missing = [name for name in ["source_sample", "is_signal", "centrality", WEIGHT_COL] if name not in matrix.files]
    if missing:
        raise SystemExit(f"{matrix_path} is missing required columns: {missing}")

    source_all = matrix["source_sample"].astype(str)
    y_all = matrix["is_signal"].astype("int8", copy=False)
    cent_all = matrix["centrality"].astype("float32", copy=False)
    weight_all = matrix[WEIGHT_COL].astype("float64", copy=False)

    source_code_parts = []
    y_parts = []
    cent_parts = []
    weight_parts = []
    score_parts = []
    total_entries = 0
    scored_entries = 0
    offset = 0
    for cache_path in cache_paths:
        cache = np.load(cache_path, allow_pickle=True)
        if SCORE_COL not in cache.files:
            raise SystemExit(f"{cache_path} is missing {SCORE_COL}")
        score = cache[SCORE_COL].astype("float32", copy=False)
        n = len(score)
        stop = offset + n
        if stop > len(y_all):
            raise SystemExit(f"Score caches exceed matrix rows: stop={stop} matrix={len(y_all)}")
        y_cache = cache["is_signal"].astype("int8", copy=False)
        if not np.array_equal(y_cache, y_all[offset:stop]):
            raise SystemExit(f"Score-cache row order mismatch at {cache_path}")
        finite = np.isfinite(score)
        source = source_all[offset:stop]
        source_code = np.zeros(n, dtype="uint8")
        source_code[np.char.find(source, "embeddedPhoton") >= 0] = 1
        source_code[np.char.find(source, "embeddedJet") >= 0] = 2
        keep = finite & (source_code > 0)
        source_code_parts.append(source_code[keep])
        y_parts.append(y_cache[keep])
        cent_parts.append(cent_all[offset:stop][keep])
        weight_parts.append(weight_all[offset:stop][keep])
        score_parts.append(score[keep])
        total_entries += n
        scored_entries += int(np.sum(finite))
        offset = stop
    if offset != len(y_all):
        raise SystemExit(f"Score caches do not cover the matrix: cache={offset} matrix={len(y_all)}")

    source_code = np.concatenate(source_code_parts)
    y = np.concatenate(y_parts)
    cent = np.concatenate(cent_parts)
    weight = np.concatenate(weight_parts)
    score = np.concatenate(score_parts)
    is_photon = source_code == 1
    is_jet = source_code == 2
    bins = np.linspace(0.0, 1.0, 51)

    rows = []
    for variant in VARIANTS:
        signal_mask, background_mask = variant_masks(variant["key"], is_photon, is_jet, y)
        if not np.any(signal_mask) or not np.any(background_mask):
            raise SystemExit(f"{variant['key']} has empty classes: S={int(signal_mask.sum())} B={int(background_mask.sum())}")
        sig_weighted = weighted_hist(score[signal_mask], weight[signal_mask], bins)
        bkg_weighted = weighted_hist(score[background_mask], weight[background_mask], bins)
        variant_row = {
            **variant,
            "global_auc": binned_auc_from_weighted_counts(sig_weighted, bkg_weighted),
            "signal_entries": int(signal_mask.sum()),
            "background_entries": int(background_mask.sum()),
            "by_centrality": {},
        }
        for key, label, lo, hi in CENTRALITY:
            cmask = np.isfinite(cent) & (cent >= lo) & (cent < hi)
            sig = signal_mask & cmask
            bkg = background_mask & cmask
            sig_cent_weighted = weighted_hist(score[sig], weight[sig], bins)
            bkg_cent_weighted = weighted_hist(score[bkg], weight[bkg], bins)
            variant_row["by_centrality"][key] = {
                "label": label,
                "auc": binned_auc_from_weighted_counts(sig_cent_weighted, bkg_cent_weighted),
                "signal": density_hist(score[sig], bins),
                "background": density_hist(score[bkg], bins),
            }
        rows.append(variant_row)

    return {
        "schema": "THE8_BRANCH_A_JET1234_DEFINITION_VARIANT_SCORE_SHAPES_V1",
        "branch": BRANCH_LABEL,
        "model_dir": str(model_dir),
        "matrix": str(matrix_path),
        "report_dir": str(report_dir),
        "score_cache_list": str(cache_list),
        "product": PRODUCT,
        "row_scope": "full allotted sample: training rows plus 10% holdout",
        "total_entries": total_entries,
        "scored_entries": scored_entries,
        "finite_score_fraction": float(scored_entries / total_entries) if total_entries else math.nan,
        "bin_edges": bins.astype(float).tolist(),
        "variants": rows,
        "class_contract_note": "All rows use the same scored Jet12+20+30+40 full sample; only the plotted red/blue class definitions change.",
        "auc_convention": "Weighted AUC computed from 50-bin weighted score histograms.",
    }


def fmt_count(value: int) -> str:
    if value >= 1_000_000:
        return f"{value / 1_000_000:.2f}M"
    if value >= 100_000:
        return f"{value / 1_000:.0f}k"
    if value >= 1_000:
        return f"{value / 1_000:.1f}k"
    return str(value)


def add_card(fig, *, xy, wh, title, body, face, edge, title_color, title_size=16.8, body_size=13.5, line_spacing=1.35):
    import matplotlib.patches as patches

    x, y = xy
    w, h = wh
    rect = patches.FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.010,rounding_size=0.010",
        transform=fig.transFigure,
        linewidth=1.05,
        edgecolor=edge,
        facecolor=face,
    )
    fig.add_artist(rect)
    fig.text(x + 0.018, y + h - 0.021, title, ha="left", va="top", fontsize=title_size, fontweight="bold", color=title_color)
    fig.text(x + 0.018, y + h - 0.052, body, ha="left", va="top", fontsize=body_size, color="#334155", linespacing=line_spacing)


def add_row_header_card(fig, *, xy, wh, variant, face, edge, sig_color, bkg_color, ink):
    import matplotlib.patches as patches

    x, y = xy
    w, h = wh
    rect = patches.FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.010,rounding_size=0.010",
        transform=fig.transFigure,
        linewidth=1.15,
        edgecolor=edge,
        facecolor=face,
    )
    fig.add_artist(rect)
    fig.add_artist(
        patches.Rectangle(
            (x + 0.004, y + 0.010),
            0.006,
            max(0.001, h - 0.020),
            transform=fig.transFigure,
            facecolor=variant["accent_strip"],
            edgecolor="none",
        )
    )
    title_text = str(variant["row_title"])
    title_has_break = "\n" in title_text
    title_size = 15.8 if title_text == "PPG12 equivalent" else (16.2 if title_has_break else 18.4)
    title_spacing = 0.90 if title_has_break else 1.00
    fig.text(
        x + 0.5 * w + 0.006,
        y + h - 0.020,
        variant["row_title"],
        ha="center",
        va="top",
        fontsize=title_size,
        fontweight="bold",
        color=ink,
        linespacing=title_spacing,
    )
    fig.text(
        x + 0.020,
        y + 0.060,
        f"S = {variant['signal_label']}",
        ha="left",
        va="center",
        fontsize=13.9,
        fontweight="bold",
        color="#B91C1C",
    )
    fig.text(
        x + 0.020,
        y + 0.033,
        f"B = {variant['background_label']}",
        ha="left",
        va="center",
        fontsize=13.9,
        fontweight="bold",
        color="#1E40AF",
    )


def add_column_header(fig, *, xy, wh, label, ink):
    import matplotlib.patches as patches

    x, y = xy
    w, h = wh
    rect = patches.FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.006,rounding_size=0.008",
        transform=fig.transFigure,
        linewidth=1.0,
        edgecolor="#C9D2DE",
        facecolor="#F8FAFC",
    )
    fig.add_artist(rect)
    fig.text(
        x + 0.5 * w,
        y + 0.5 * h,
        f"{label} centrality",
        transform=fig.transFigure,
        ha="center",
        va="center",
        fontsize=14.3,
        fontweight="bold",
        color=ink,
    )


def draw_density(ax, edges: np.ndarray, density: np.ndarray, color: str, linestyle: str) -> None:
    y = np.r_[density, density[-1]]
    ax.fill_between(edges, y, step="post", color=color, alpha=0.10, linewidth=0)
    ax.step(edges, y, where="post", color=color, lw=2.25, linestyle=linestyle)


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
            "axes.linewidth": 1.05,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )

    variant_display = {variant["key"]: variant for variant in VARIANTS}
    variants = []
    for variant in payload["variants"]:
        merged = dict(variant)
        display = variant_display.get(str(variant.get("key", "")), {})
        for field in ("label", "row_label", "row_title", "signal_label", "background_label", "signal_definition", "background_definition"):
            if field in display:
                merged[field] = display[field]
        variants.append(merged)
    edges = np.asarray(payload["bin_edges"], dtype=float)
    sig_color = "#DC2626"
    bkg_color = "#1D4ED8"
    ink = "#111827"
    muted = "#475569"
    light_grid = "#E5E7EB"
    ymax = 0.0
    for variant in variants:
        for cent_key, _, _, _ in CENTRALITY:
            item = variant["by_centrality"][cent_key]
            ymax = max(ymax, float(np.max(item["signal"]["density"])), float(np.max(item["background"]["density"])))
    ymax = 1.18 * ymax if ymax > 0 else 1.0

    fig, axes = plt.subplots(
        3,
        3,
        figsize=(16, 9),
        dpi=160,
        sharex=True,
        sharey=True,
        gridspec_kw={"left": 0.245, "right": 0.968, "bottom": 0.205, "top": 0.648, "hspace": 0.35, "wspace": 0.115},
    )
    fig.patch.set_facecolor("white")

    fig.text(0.035, 0.967, "Definition check: what counts as signal and background?", ha="left", va="top", fontsize=27.5, fontweight="bold", color=ink)
    add_card(
        fig,
        xy=(0.045, 0.805),
        wh=(0.910, 0.112),
        title="Same scored embedded Photon+Jet 12+20 and embedded Inclusive Jet 12+20+30+40",
        body=(
            "BDT training and AUC use PPG12-exact cluster ET/eta weights; only the red/blue class definitions change by row.\n"
            "Curves are raw-entry unit-area score densities; S/B labels are raw entry counts."
        ),
        face="#F8FAFC",
        edge="#CBD5E1",
        title_color=ink,
        title_size=17.2,
        body_size=14.2,
        line_spacing=1.45,
    )

    legend = patches.FancyBboxPatch(
        (0.190, 0.717),
        0.620,
        0.050,
        boxstyle="round,pad=0.006,rounding_size=0.010",
        transform=fig.transFigure,
        linewidth=1.05,
        edgecolor="#D1D5DB",
        facecolor="#FFFFFF",
    )
    fig.add_artist(legend)
    yc = 0.742
    fig.add_artist(matplotlib.lines.Line2D([0.230, 0.290], [yc, yc], transform=fig.transFigure, color=sig_color, lw=3.5))
    fig.text(0.300, yc, "Red curve uses row S definition", transform=fig.transFigure, ha="left", va="center", fontsize=13.4, color=ink)
    fig.add_artist(matplotlib.lines.Line2D([0.575, 0.635], [yc, yc], transform=fig.transFigure, color=bkg_color, lw=3.5, linestyle="--"))
    fig.text(0.645, yc, "Blue curve uses row B definition", transform=fig.transFigure, ha="left", va="center", fontsize=13.4, color=ink)

    summary_rows = []
    for irow, variant in enumerate(variants):
        for icol, (cent_key, cent_label, _, _) in enumerate(CENTRALITY):
            ax = axes[irow, icol]
            item = variant["by_centrality"][cent_key]
            sig_density = np.asarray(item["signal"]["density"], dtype=float)
            bkg_density = np.asarray(item["background"]["density"], dtype=float)
            draw_density(ax, edges, sig_density, sig_color, "-")
            draw_density(ax, edges, bkg_density, bkg_color, "--")
            ax.set_xlim(0.0, 1.0)
            ax.set_ylim(0.0, ymax)
            ax.grid(True, color=light_grid, lw=0.55, alpha=0.85)
            ax.tick_params(labelsize=11.4, pad=2)
            ax.set_facecolor("#FFFFFF")
            for spine in ax.spines.values():
                spine.set_color(ink)
                spine.set_linewidth(1.05)
            if irow == 0:
                pos = ax.get_position()
                add_column_header(
                    fig,
                    xy=(pos.x0 + 0.004, pos.y1 + 0.016),
                    wh=(pos.width - 0.008, 0.034),
                    label=cent_label,
                    ink=ink,
                )
            if icol == 0:
                ax.set_ylabel("Unit-area\ndensity", fontsize=12.0, color=ink, labelpad=8)
                pos = ax.get_position()
                accent = DEFINITION_ACCENTS[irow % len(DEFINITION_ACCENTS)]
                variant["accent_strip"] = accent["strip"]
                add_row_header_card(
                    fig,
                    xy=(0.038, pos.y0 - 0.004),
                    wh=(0.148, pos.height + 0.008),
                    variant=variant,
                    face=accent["fill"],
                    edge=accent["edge"],
                    sig_color=sig_color,
                    bkg_color=bkg_color,
                    ink=ink,
                )
            if irow == 2:
                ax.set_xlabel("BDT score", fontsize=13.8, color=ink)
            if irow == 0 and icol == 0:
                ax.text(0.035, 0.900, "sPHENIX", transform=ax.transAxes, ha="left", va="top", fontsize=12.3, fontstyle="italic", fontweight="bold", color=ink)
                ax.text(0.238, 0.900, "Internal", transform=ax.transAxes, ha="left", va="top", fontsize=12.3, color=ink)
                ax.text(0.035, 0.755, "full sample", transform=ax.transAxes, ha="left", va="top", fontsize=11.3, color=ink)
            sig_entries = int(item["signal"]["entries"])
            bkg_entries = int(item["background"]["entries"])
            auc = float(item["auc"])
            ax.text(
                0.965,
                0.920,
                f"AUC {auc:.3f}",
                transform=ax.transAxes,
                ha="right",
                va="top",
                fontsize=15.3,
                fontweight="bold",
                color=ink,
                bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#CBD5E1", lw=0.9, alpha=0.96),
            )
            ax.text(0.965, 0.675, f"S {fmt_count(sig_entries)}   B {fmt_count(bkg_entries)}", transform=ax.transAxes, ha="right", va="top", fontsize=12.7, fontweight="bold", color=muted)
            summary_rows.append(
                {
                    "variant": variant["label"],
                    "variant_key": variant["key"],
                    "centrality": cent_label,
                    "auc_weighted_binned": auc,
                    "signal_entries": sig_entries,
                    "background_entries": bkg_entries,
                    "signal_definition": variant["signal_definition"],
                    "background_definition": variant["background_definition"],
                }
            )

    global_aucs = {str(v["key"]): float(v["global_auc"]) for v in variants}
    auc_min = min(global_aucs.values())
    auc_max = max(global_aucs.values())
    auc_spread = auc_max - auc_min
    current_auc = global_aucs["current"]
    truth_delta = global_aucs["truth_tagged"] - current_auc
    source_delta = global_aucs["source_only"] - current_auc
    add_card(
        fig,
        xy=(0.038, 0.026),
        wh=(0.924, 0.110),
        title="Definition choice has a small effect on the score-separation metric",
        body=(
            f"Weighted full-sample AUC changes only {auc_min:.3f}-{auc_max:.3f} (spread {auc_spread:.3f}).\n"
            f"Relative to PPG12 equivalent: truth-tagged {truth_delta:+.3f}; source-only {source_delta:+.3f}."
        ),
        face="#FFF7ED",
        edge="#FDBA74",
        title_color="#9A3412",
        title_size=18.4,
        body_size=15.6,
        line_spacing=1.30,
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
    json_path.write_text(
        json.dumps(
            {
                "schema": "THE8_BRANCH_A_JET1234_DEFINITION_VARIANT_SLIDE_V1",
                "input_schema": payload.get("schema", ""),
                "branch": payload.get("branch", ""),
                "row_scope": payload.get("row_scope", ""),
                "class_contract_note": payload.get("class_contract_note", ""),
                "auc_convention": payload.get("auc_convention", ""),
                "product": payload.get("product", ""),
                "model_dir": payload.get("model_dir", ""),
                "matrix": payload.get("matrix", ""),
                "report_dir": payload.get("report_dir", ""),
                "score_cache_list": payload.get("score_cache_list", ""),
                "png": str(png_path),
                "summary_csv": str(csv_path),
                "variants": [
                    {
                        "key": v["key"],
                        "label": v["label"],
                        "signal_definition": v["signal_definition"],
                        "background_definition": v["background_definition"],
                        "global_auc": v["global_auc"],
                        "signal_entries": v["signal_entries"],
                        "background_entries": v["background_entries"],
                    }
                    for v in variants
                ],
                "rows": summary_rows,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    md_path.write_text(make_script_text(variants) + "\n")
    print(png_path)
    print(csv_path)
    print(json_path)
    print(md_path)


def make_script_text(variants: list[dict[str, object]]) -> str:
    values = "; ".join(f"{v['label']} AUC {float(v['global_auc']):.3f}" for v in variants)
    return f"""# WP GammaJets Slide Script - Score Definition Variants

Now I want to pause on a subtle but important point in the score-shape plots: what exactly we mean by the red signal curve and the blue background curve.

Every panel on this slide uses the same Jet12+20+30+40 full allotted simulation sample and the same BDT score. The only thing changing by row is how I define the candidates that get drawn as red and blue.

The top row is the PPG12-equivalent diagnostic definition. The red curve is truth-isolated prompt candidates from embedded-photon Monte Carlo, while the blue curve is all candidates from embedded inclusive-jet Monte Carlo. This matches the signal-template versus inclusive-candidate-population question.

The middle row shows what happens if I also truth-tag the inclusive-jet side as background. In that case the blue curve excludes truth-positive jet-source candidates, so it is a more symmetric truth-class diagnostic but no longer represents the full inclusive-jet candidate population.

The bottom row removes truth tagging from both source samples. That shows the raw source-sample comparison: all embedded-photon candidates against all embedded-jet candidates. This is useful as a source-composition check, but the red curve is no longer a clean prompt-photon signal template.

The binned weighted full-sample AUCs are {values}. The point is not to replace the held-out validation claim, but to make clear which score-shape definition is being used and why the PPG12-equivalent asymmetric definition is the one that matches the signal-template versus inclusive-MC diagnostic.
"""


def main() -> int:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)
    extract = sub.add_parser("extract")
    extract.add_argument("--model-dir", type=Path, default=DEFAULT_MODEL_DIR)
    extract.add_argument("--report-dir", type=Path, default=DEFAULT_REPORT_DIR)
    render = sub.add_parser("render")
    render.add_argument("--input", type=Path, required=True)
    render.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    render.add_argument("--tag", default="the8_branchA_jet12_20_30_40_definition_variants_fullAllotted_v1")
    args = parser.parse_args()

    if args.command == "extract":
        print(json.dumps(extract_payload(args.model_dir, args.report_dir), indent=2, sort_keys=True))
        return 0
    payload = json.loads(args.input.read_text())
    render_slide(payload, args.outdir, args.tag)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
