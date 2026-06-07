#!/usr/bin/env python3
"""Make THE-8 Branch A ladder score-separation slide from compact validation histograms."""

from __future__ import annotations
# Keep purpose-folder helpers runnable when invoked directly.
import sys as _codex_sys
from pathlib import Path as _CodexPath
_CODEX_THIS_FILE = _CodexPath(__file__).resolve()
_CODEX_SCRIPTS_DIR = next((p for p in _CODEX_THIS_FILE.parents if p.name == "scripts"), _CODEX_THIS_FILE.parent)
_CODEX_SCRIPTS_DIR_STR = str(_CODEX_SCRIPTS_DIR)
if _CODEX_SCRIPTS_DIR_STR not in _codex_sys.path:
    _codex_sys.path.append(_CODEX_SCRIPTS_DIR_STR)
del _CODEX_THIS_FILE, _CODEX_SCRIPTS_DIR, _CODEX_SCRIPTS_DIR_STR

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np


CENTRALITY = [("0_20", "0-20%"), ("20_50", "20-50%"), ("50_80", "50-80%")]
SAMPLE_ACCENTS = (
    {"fill": "#ECF7EE", "edge": "#6EA77B", "strip": "#3F8F5A"},
    {"fill": "#FFF4CC", "edge": "#D6A84A", "strip": "#B88400"},
    {"fill": "#F3ECFA", "edge": "#9A7CC3", "strip": "#7651A6"},
)


def binned_auc(signal_counts: np.ndarray, background_counts: np.ndarray) -> float:
    sig_total = float(np.sum(signal_counts))
    bkg_total = float(np.sum(background_counts))
    if sig_total <= 0 or bkg_total <= 0:
        return float("nan")
    bkg_below = np.cumsum(background_counts) - background_counts
    wins = float(np.sum(signal_counts * bkg_below))
    ties = float(np.sum(signal_counts * background_counts))
    return (wins + 0.5 * ties) / (sig_total * bkg_total)


def wp_fake_rate_at_signal_eff(
    signal_counts: np.ndarray,
    background_counts: np.ndarray,
    *,
    edges: np.ndarray,
    target_efficiency: float = 0.80,
) -> dict[str, float]:
    """Approximate WP fake rate from the same low-to-high score-bin counts."""
    sig_total = float(np.sum(signal_counts))
    bkg_total = float(np.sum(background_counts))
    if sig_total <= 0 or bkg_total <= 0:
        return {"threshold": float("nan"), "signal_efficiency": float("nan"), "background_fake_rate": float("nan")}
    signal_above = np.cumsum(signal_counts[::-1])[::-1] / sig_total
    background_above = np.cumsum(background_counts[::-1])[::-1] / bkg_total
    idx = int(np.nanargmin(np.abs(signal_above - target_efficiency)))
    return {
        "threshold": float(edges[idx]),
        "signal_efficiency": float(signal_above[idx]),
        "background_fake_rate": float(background_above[idx]),
    }


def step(ax: plt.Axes, edges: np.ndarray, density: np.ndarray, color: str, label: str, linestyle: str) -> None:
    y = np.r_[density, density[-1]]
    ax.step(edges, y, where="post", color=color, lw=1.9, linestyle=linestyle, label=label)


def draw_density(ax: plt.Axes, edges: np.ndarray, density: np.ndarray, color: str, label: str, *, linestyle: str) -> None:
    y = np.r_[density, density[-1]]
    ax.fill_between(edges, y, step="post", color=color, alpha=0.12, linewidth=0)
    ax.step(edges, y, where="post", color=color, lw=2.25, linestyle=linestyle, label=label)


def fmt_millions(value: int) -> str:
    if value >= 1_000_000:
        return f"{value / 1_000_000:.2f}M"
    if value >= 1_000:
        return f"{value / 1_000:.0f}k"
    return str(value)


def auc_color(auc: float) -> str:
    if not np.isfinite(auc):
        return "#64748B"
    if auc >= 0.88:
        return "#047857"
    if auc >= 0.82:
        return "#2563EB"
    return "#B45309"


def add_card(
    fig: plt.Figure,
    *,
    xy: tuple[float, float],
    wh: tuple[float, float],
    title: str,
    body: str,
    face: str,
    edge: str,
    title_color: str,
    title_size: float = 18.6,
    body_size: float = 15.4,
    title_offset: float = 0.024,
    body_offset: float = 0.066,
    line_spacing: float = 1.10,
    align: str = "left",
) -> None:
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
    text_x = x + (0.5 * w if align == "center" else 0.018)
    fig.text(text_x, y + h - title_offset, title, ha=align, va="top", fontsize=title_size, fontweight="bold", color=title_color)
    fig.text(text_x, y + h - body_offset, body, ha=align, va="top", fontsize=body_size, color="#334155", linespacing=line_spacing)


def add_legend_band(fig: plt.Figure, *, x: float, y: float, w: float, h: float, sig_color: str, bkg_color: str, ink: str) -> None:
    rect = patches.FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.006,rounding_size=0.010",
        transform=fig.transFigure,
        linewidth=1.05,
        edgecolor="#D1D5DB",
        facecolor="#FFFFFF",
    )
    fig.add_artist(rect)
    yc = y + 0.5 * h
    fig.add_artist(matplotlib.lines.Line2D([x + 0.026, x + 0.078], [yc, yc], transform=fig.transFigure, color=sig_color, lw=3.5))
    fig.text(x + 0.086, yc, "Signal MC (truth-isolated prompt)", transform=fig.transFigure, ha="left", va="center", fontsize=13.2, color=ink)
    fig.add_artist(
        matplotlib.lines.Line2D(
            [x + 0.382, x + 0.434],
            [yc, yc],
            transform=fig.transFigure,
            color=bkg_color,
            lw=3.5,
            linestyle="--",
        )
    )
    fig.text(
        x + 0.442,
        yc,
        "Inclusive MC (embedded jet)",
        transform=fig.transFigure,
        ha="left",
        va="center",
        fontsize=13.2,
        color=ink,
    )


def formatted_sample_label(label: str) -> str:
    if label == "Jet12+20+30":
        return "Jet12+20\n+30"
    if label == "Jet12+20+30+40":
        return "Jet12+20\n+30+40"
    return label


def add_column_header(fig: plt.Figure, *, xy: tuple[float, float], wh: tuple[float, float], label: str, ink: str) -> None:
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


def add_row_header_card(
    fig: plt.Figure,
    *,
    xy: tuple[float, float],
    wh: tuple[float, float],
    role: str,
    label: str,
    accent: dict[str, str],
    ink: str,
) -> None:
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
    fig.text(
        x + 0.5 * w + 0.006,
        y + 0.5 * h + 0.028,
        role,
        transform=fig.transFigure,
        ha="center",
        va="center",
        fontsize=12.4,
        fontweight="bold",
        color=ink,
    )
    fig.text(
        x + 0.5 * w + 0.006,
        y + 0.5 * h - 0.008,
        formatted_sample_label(label),
        transform=fig.transFigure,
        ha="center",
        va="center",
        fontsize=13.6 if "\n" in formatted_sample_label(label) else 14.3,
        fontweight="bold",
        color=ink,
        linespacing=1.18,
    )


def slide_copy(branches: list[dict]) -> dict[str, str]:
    modes = {str(b.get("summary", {}).get("validation_mode", "")) for b in branches}
    if modes == {"own_full_allotted_sample_truth_signal_inclusive_jet"}:
        return {
            "title": "Full-sample diagnostic: training + holdout score separation",
            "subtitle": "",
            "left_title": "Rows shown in each pad",
            "left_body": "Red is truth-isolated prompt candidates from embedded-photon MC;\nblue is all candidates from embedded inclusive-jet MC.",
            "right_title": "Training + holdout diagnostic",
            "right_body": "Each BDT scores its full matching sample:\nJet12+20, then +Jet30, then +Jet40.",
            "bottom_title": "Weighted full-sample AUC by row",
            "bottom_suffix": "Full sample includes training rows plus the 10% holdout; use prior slide for held-out validation.",
            "plot_sample": "training+holdout",
        }
    if modes == {"own_10pct_training_holdout_truth_signal_inclusive_jet"}:
        return {
            "title": "Held-out validation: truth-isolated signal vs inclusive-jet MC",
            "subtitle": r"Global no-isolation photon-ID BDT; each row uses only its own 10% training holdout, $15 < E_T < 35$ GeV.",
            "left_title": "Validation rows in each pad",
            "left_body": "Red is truth-isolated prompt candidates from embedded-photon MC;\nblue is all candidates from embedded inclusive-jet MC.",
            "right_title": "Rows use independent holdout samples",
            "right_body": "Each BDT is validated on its matching training sample:\nJet12+20, then +Jet30, then +Jet40.",
            "bottom_title": "Weighted holdout AUC by row",
            "bottom_suffix": "Pad AUC/S/B labels use truth-isolated signal and unfiltered inclusive-jet MC rows.",
            "plot_sample": "10% row holdout",
        }
    if modes == {"common_jet12_20_10pct_training_holdout_truth_signal_inclusive_jet"}:
        return {
            "title": "Common Jet12+20 holdout: truth-isolated signal vs inclusive-jet MC",
            "subtitle": r"All three independently trained BDTs are scored on the same Jet12+20 10% holdout, $15 < E_T < 35$ GeV.",
            "left_title": "Same validation rows in every pad",
            "left_body": "Red is truth-isolated prompt candidates from embedded-photon MC;\nblue is all candidates from the same inclusive-jet MC split.",
            "right_title": "This is the same-sample control",
            "right_body": "Any row-to-row change here comes from the trained model,\nnot from adding Jet30 or Jet40 rows to validation.",
            "bottom_title": "Weighted AUC on common Jet12+20 holdout",
            "bottom_suffix": "Pad AUC/S/B labels use the same truth-isolated signal and inclusive-jet rows in every row.",
            "plot_sample": "Jet12+20 holdout",
        }
    if modes == {"own_10pct_training_holdout"}:
        return {
            "title": "Held-out validation confirms the Branch A ladder separation gain",
            "subtitle": r"Global no-isolation photon-ID BDT; each row uses only its own 10% training holdout, $15 < E_T < 35$ GeV.",
            "left_title": "Validation rows in each pad",
            "left_body": "Red uses embedded-photon rows; blue uses all embedded-jet rows\nfrom the held-out 10% split named on that row.",
            "right_title": "Rows use independent holdout samples",
            "right_body": "Each BDT is validated on its matching training sample:\nJet12+20, then +Jet30, then +Jet40.",
            "bottom_title": "Weighted holdout AUC by row",
            "bottom_suffix": "Pad AUC/S/B labels use source-defined sample rows within each centrality bin.",
            "plot_sample": "10% row holdout",
        }
    if modes == {"common_jet12_20_10pct_training_holdout"}:
        return {
            "title": "Common Jet12+20 holdout isolates the model-comparison effect",
            "subtitle": r"All three independently trained BDTs are scored on the same Jet12+20 10% holdout, $15 < E_T < 35$ GeV.",
            "left_title": "Same validation rows in every pad",
            "left_body": "Red uses embedded-photon rows; blue uses all embedded-jet rows\nfrom the same Jet12+20 held-out split.",
            "right_title": "This is the same-sample control",
            "right_body": "Any row-to-row change here comes from the trained model,\nnot from adding Jet30 or Jet40 rows to validation.",
            "bottom_title": "Weighted AUC on common Jet12+20 holdout",
            "bottom_suffix": "Pad AUC/S/B labels use source-defined sample rows in every row.",
            "plot_sample": "Jet12+20 holdout",
        }
    return {
        "title": "Jet40-inclusive validation background improves BDT score separation",
        "subtitle": r"Global no-isolation photon-ID BDT; embedded prompt-photon signal vs embedded inclusive-jet background, $15 < E_T < 35$ GeV.",
        "left_title": "Validation sample in each pad",
        "left_body": "Embedded prompt-photon signal is compared to embedded inclusive-jet\nbackground from the combined sample named on that row.",
        "right_title": "Rows are not one common inclusive pool",
        "right_body": "Each row validates the matching BDT against its matching combined\ninclusive pool: Jet12+20, then +Jet30, then +Jet40.",
        "bottom_title": "Inclusive validation AUC by row",
        "bottom_suffix": "Pad AUC/S/B labels are computed within each centrality bin.",
        "plot_sample": "Au+Au embedded validation",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--tag", default="the8_branchA_ladder_score_separation_slide23_style")
    parser.add_argument(
        "--annotate-wp80-fake",
        action="store_true",
        help="Print the binned WP80 background fake rate below each panel AUC.",
    )
    args = parser.parse_args()

    payload = json.loads(args.input.read_text())
    branches = payload["branches"]
    copy = slide_copy(branches)
    args.outdir.mkdir(parents=True, exist_ok=True)

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

    fig, axes = plt.subplots(
        3,
        3,
        figsize=(16, 9),
        sharex=True,
        sharey=True,
        dpi=160,
        gridspec_kw={"left": 0.190, "right": 0.968, "bottom": 0.205, "top": 0.590, "hspace": 0.34, "wspace": 0.115},
    )
    fig.patch.set_facecolor("white")
    sig_color = "#DC2626"
    bkg_color = "#1D4ED8"
    ink = "#111827"
    muted = "#475569"
    light_grid = "#E5E7EB"
    summary_rows: list[dict[str, object]] = []

    ymax = 0.0
    cache: list[dict[str, object]] = []
    for branch in branches:
        edges = np.asarray(branch["bin_edges"], dtype=float)
        for key, label in CENTRALITY:
            cent = branch["by_centrality"][key]
            sig_density = np.asarray(cent["signal"]["density"], dtype=float)
            bkg_density = np.asarray(cent["background"]["density"], dtype=float)
            sig_counts = np.asarray(cent["signal"]["counts"], dtype=float)
            bkg_counts = np.asarray(cent["background"]["counts"], dtype=float)
            item = {
                "branch": branch,
                "cent_key": key,
                "cent_label": label,
                "edges": edges,
                "sig_density": sig_density,
                "bkg_density": bkg_density,
                "sig_counts": sig_counts,
                "bkg_counts": bkg_counts,
                "auc": binned_auc(sig_counts, bkg_counts),
                "wp80": wp_fake_rate_at_signal_eff(sig_counts, bkg_counts, edges=edges),
            }
            cache.append(item)
            ymax = max(ymax, float(np.nanmax(sig_density)), float(np.nanmax(bkg_density)))
    ymax = 1.18 * ymax if ymax > 0 else 1.0

    for irow, branch in enumerate(branches):
        for icol, (_, cent_label) in enumerate(CENTRALITY):
            ax = axes[irow, icol]
            item = next(x for x in cache if x["branch"] is branch and x["cent_label"] == cent_label)
            draw_density(ax, item["edges"], item["sig_density"], sig_color, "Signal", linestyle="-")
            draw_density(ax, item["edges"], item["bkg_density"], bkg_color, "Background", linestyle="--")
            ax.set_xlim(0.0, 1.0)
            ax.set_ylim(0.0, ymax)
            ax.grid(True, color=light_grid, lw=0.55, alpha=0.85)
            ax.tick_params(labelsize=12.0, pad=3)
            ax.set_facecolor("#FFFFFF")
            for spine in ax.spines.values():
                spine.set_color("#111827")
                spine.set_linewidth(1.05)
            if irow == 0:
                pos = ax.get_position()
                add_column_header(
                    fig,
                    xy=(pos.x0 + 0.004, pos.y1 + 0.017),
                    wh=(pos.width - 0.008, 0.034),
                    label=cent_label,
                    ink=ink,
                )
            if icol == 0:
                ax.set_ylabel("Unit-area\ndensity", fontsize=11.8, color=ink, labelpad=8)
                pos = ax.get_position()
                add_row_header_card(
                    fig,
                    xy=(0.028, pos.y0 + 0.010),
                    wh=(0.098, pos.height - 0.020),
                    role="SAMPLE" if copy["plot_sample"] == "training+holdout" else "VALIDATION",
                    label=branch["label"],
                    accent=SAMPLE_ACCENTS[irow % len(SAMPLE_ACCENTS)],
                    ink=ink,
                )
            if irow == 2:
                ax.set_xlabel("BDT score", fontsize=14.2, color=ink)
            if irow == 0 and icol == 0:
                ax.text(
                    0.035,
                    0.900,
                    "sPHENIX",
                    transform=ax.transAxes,
                    ha="left",
                    va="top",
                    fontsize=12.8,
                    fontstyle="italic",
                    fontweight="bold",
                    color=ink,
                )
                ax.text(
                    0.240,
                    0.900,
                    "Internal",
                    transform=ax.transAxes,
                    ha="left",
                    va="top",
                    fontsize=12.8,
                    color=ink,
                )
                ax.text(
                    0.035,
                    0.760,
                    copy["plot_sample"],
                    transform=ax.transAxes,
                    ha="left",
                    va="top",
                    fontsize=11.8,
                    color=ink,
                )
            signal_entries = int(np.sum(item["sig_counts"]))
            background_entries = int(np.sum(item["bkg_counts"]))
            auc = float(item["auc"])
            wp80 = item["wp80"]
            ax.text(
                0.965,
                0.935,
                f"AUC {auc:.3f}",
                transform=ax.transAxes,
                ha="right",
                va="top",
                fontsize=16.0,
                fontweight="bold",
                color=ink,
            )
            s_b_y = 0.665
            if args.annotate_wp80_fake:
                wp80_fake = float(wp80["background_fake_rate"])
                ax.text(
                    0.965,
                    0.735,
                    f"WP80 fake {100.0 * wp80_fake:.1f}%",
                    transform=ax.transAxes,
                    ha="right",
                    va="top",
                    fontsize=15.0,
                    fontweight="bold",
                    color="#C2410C",
                    bbox=dict(boxstyle="round,pad=0.20", fc="#FFF7ED", ec="#FDBA74", lw=0.85, alpha=0.96),
                )
                s_b_y = 0.545
            ax.text(
                0.965,
                s_b_y,
                f"S {fmt_millions(signal_entries)}   B {fmt_millions(background_entries)}",
                transform=ax.transAxes,
                ha="right",
                va="top",
                fontsize=13.4,
                fontweight="bold",
                color=muted,
            )
            summary_rows.append(
                {
                    "branch": branch["label"],
                    "centrality": cent_label,
                    "auc_binned": item["auc"],
                    "wp80_threshold_binned": wp80["threshold"],
                    "wp80_signal_efficiency_binned": wp80["signal_efficiency"],
                    "wp80_background_fake_rate_binned": wp80["background_fake_rate"],
                    "signal_entries": signal_entries,
                    "background_entries": background_entries,
                    "inclusive_auc": branch["summary"].get("globalEtCent1535_bdt_noIso_auc", ""),
                    "finite_score_fraction": branch["summary"].get("finite_score_fraction", ""),
                    "total_entries": branch["summary"].get("total_entries", ""),
                    "scored_entries": branch["summary"].get("scored_entries", ""),
                }
            )

    jet12_auc = float(branches[0]["summary"].get("globalEtCent1535_bdt_noIso_auc", "nan"))
    jet123_auc = float(branches[1]["summary"].get("globalEtCent1535_bdt_noIso_auc", "nan"))
    jet1234_auc = float(branches[2]["summary"].get("globalEtCent1535_bdt_noIso_auc", "nan"))
    bottom_suffix = copy["bottom_suffix"]
    if args.annotate_wp80_fake:
        bottom_suffix = "Orange pad labels show background accepted at binned WP80."
    fig.text(
        0.035,
        0.967,
        copy["title"],
        ha="left",
        va="top",
        fontsize=28.0,
        fontweight="bold",
        color=ink,
    )
    add_card(
        fig,
        xy=(0.041, 0.732),
        wh=(0.425, 0.150),
        title=copy["left_title"],
        body=copy["left_body"],
        face="#F8FAFC",
        edge="#CBD5E1",
        title_color=ink,
        title_size=19.0,
        body_size=15.7,
        body_offset=0.064,
        line_spacing=1.50,
    )
    add_card(
        fig,
        xy=(0.534, 0.732),
        wh=(0.425, 0.150),
        title=copy["right_title"],
        body=copy["right_body"],
        face="#EEF6FF",
        edge="#93C5FD",
        title_color="#1D4ED8",
        title_size=19.0,
        body_size=15.7,
        body_offset=0.064,
        line_spacing=1.50,
    )
    add_legend_band(fig, x=0.168, y=0.660, w=0.704, h=0.052, sig_color=sig_color, bkg_color=bkg_color, ink=ink)
    add_card(
        fig,
        xy=(0.038, 0.036),
        wh=(0.924, 0.072),
        title=copy["bottom_title"],
        body=(
            f"Jet12+20: {jet12_auc:.3f}   |   Jet12+20+30: {jet123_auc:.3f}   |   "
            f"Jet12+20+30+40: {jet1234_auc:.3f}      "
            f"{bottom_suffix}"
        ),
        face="#FFF7ED",
        edge="#FDBA74",
        title_color="#9A3412",
        title_offset=0.018,
        body_offset=0.047,
        line_spacing=1.20,
    )

    png = args.outdir / f"{args.tag}.png"
    csv_path = args.outdir / f"{args.tag}.csv"
    json_path = args.outdir / f"{args.tag}.json"
    fig.savefig(png)
    plt.close(fig)

    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary_rows[0]))
        writer.writeheader()
        writer.writerows(summary_rows)
    plot_class_definition = str(payload.get("plot_class_definition", ""))
    plot_class_definition = plot_class_definition.replace(
        "source_sample contains embeddedPhoton and is_signal == 1;",
        "source_sample contains embeddedPhoton and is_signal == 1 (truth-isolated prompt label);",
    )
    json_path.write_text(
        json.dumps(
            {
                "schema": "THE8_BRANCH_A_LADDER_SCORE_SEPARATION_PANEL_V1",
                "input": str(args.input),
                "input_schema": payload.get("schema", ""),
                "description": payload.get("description", ""),
                "plot_class_definition": plot_class_definition,
                "png": str(png),
                "summary_csv": str(csv_path),
                "rows": summary_rows,
                "note": (
                    "Panel AUC and WP80 fake-rate values are computed from compact binned histograms; "
                    "row-summary AUC values are read from the input branch summaries."
                ),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    print(png)
    print(csv_path)
    print(json_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
