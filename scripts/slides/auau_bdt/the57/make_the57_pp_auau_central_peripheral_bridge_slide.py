#!/usr/bin/env python3
"""Regenerate the working-point deck pp/AuAu BDT bridge slide with THE-57 output."""

from __future__ import annotations

import csv
import json
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredOffsetbox, HPacker, TextArea
import numpy as np


def find_repo() -> Path:
    here = Path(__file__).resolve()
    for parent in here.parents:
        if (parent / "AGENTS.md").exists() and (parent / "agent_context").exists():
            return parent
    raise RuntimeError("Could not resolve ThesisAnalysis repo root")


REPO = find_repo()
PP_SUMMARY = (
    REPO
    / "dataOutput/ppPhotonMLPipeline/ppg12_basev3E_currentIAN_rawOverlayEnvFix_20260527_1630"
    / "validation/fullsim_shuhang_overlay_raw_inclusive_pt1535"
    / "pp_currentian_basev3e_bdt_score_overlay_vs_ppg12_pp_noCent_bdt_15_35_noNPB_eta0-pt3-cut0_summary.json"
)
AUAU_SCORE_CACHE_DIR = (
    REPO
    / "dataOutput/auauTightBDTValidation/THE57_baseline_validation_20260615"
    / "a1_withcut_score_caches/score_caches"
)
PRODUCT = "centAsFeatBase3x3_pt15to35"
SCORE_KEY = f"score_{PRODUCT}"
CENTRALITIES = [
    ("0_20", "0-20%", "Central Au+Au", 0.0, 20.0),
    ("50_80", "50-80%", "Peripheral Au+Au", 50.0, 80.0),
]

OUT_DIR = REPO / "dataOutput/auauTightBDTValidation/THE57_baseline_validation_20260615/slide9_pp_auau_bridge"
OUT_PNG = OUT_DIR / "the57_pp_auau_central_peripheral_bdt_bridge.png"
OUT_JSON = OUT_DIR / "the57_pp_auau_central_peripheral_bdt_bridge_manifest.json"
OUT_CSV = OUT_DIR / "the57_pp_auau_central_peripheral_bdt_bridge_histograms.csv"
OUT_SCRIPT = OUT_DIR / "the57_pp_auau_central_peripheral_bdt_bridge_script.md"

INK = "#111827"
MUTED = "#4b5563"
GRID = "#d9dee8"
SIGNAL_RED = "#d62728"
BACKGROUND_BLUE = "#1f77b4"
CENTRAL_EDGE = "#d95f5f"
CENTRAL_HEADER = "#fff0f0"
PERIPHERAL_EDGE = "#2f8f67"
PERIPHERAL_HEADER = "#edf8f0"


def density_to_prob(edges: np.ndarray, density: np.ndarray) -> np.ndarray:
    widths = np.diff(edges)
    probs = np.asarray(density, dtype=float) * widths
    total = float(np.sum(probs))
    if total <= 0:
        raise ValueError("non-positive histogram probability total")
    return probs / total


def auc_from_probs(sig_prob: np.ndarray, bkg_prob: np.ndarray) -> float:
    bkg_below = np.r_[0.0, np.cumsum(bkg_prob)[:-1]]
    favorable = float(np.sum(sig_prob * bkg_below))
    ties = 0.5 * float(np.sum(sig_prob * bkg_prob))
    return favorable + ties


def weighted_auc_exact(signal_scores: np.ndarray, background_scores: np.ndarray, signal_w: np.ndarray, background_w: np.ndarray) -> float:
    scores = np.concatenate([signal_scores, background_scores])
    labels = np.concatenate([np.ones_like(signal_scores, dtype=bool), np.zeros_like(background_scores, dtype=bool)])
    weights = np.concatenate([signal_w, background_w]).astype(float)
    order = np.argsort(scores, kind="mergesort")
    scores = scores[order]
    labels = labels[order]
    weights = weights[order]
    total_sig = float(np.sum(weights[labels]))
    total_bkg = float(np.sum(weights[~labels]))
    if total_sig <= 0 or total_bkg <= 0:
        return float("nan")
    auc_num = 0.0
    bkg_below = 0.0
    start = 0
    n = scores.size
    while start < n:
        end = start + 1
        while end < n and scores[end] == scores[start]:
            end += 1
        group_labels = labels[start:end]
        group_weights = weights[start:end]
        sig_w = float(np.sum(group_weights[group_labels]))
        bkg_w = float(np.sum(group_weights[~group_labels]))
        auc_num += sig_w * (bkg_below + 0.5 * bkg_w)
        bkg_below += bkg_w
        start = end
    return auc_num / (total_sig * total_bkg)


def roc_from_probs(sig_prob: np.ndarray, bkg_prob: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    tpr = np.r_[0.0, np.cumsum(sig_prob[::-1])]
    fpr = np.r_[0.0, np.cumsum(bkg_prob[::-1])]
    return fpr, tpr


def load_pp() -> dict[str, object]:
    payload = json.loads(PP_SUMMARY.read_text())
    edges = np.asarray(payload["bins"], dtype=float)
    sig = np.asarray(payload["this_analysis_signal_hist"], dtype=float)
    bkg = np.asarray(payload["this_analysis_inclusive_hist"], dtype=float)
    sig_prob = density_to_prob(edges, sig)
    bkg_prob = density_to_prob(edges, bkg)
    fpr, tpr = roc_from_probs(sig_prob, bkg_prob)
    return {
        "label": "pp",
        "edges": edges,
        "signal_density": sig,
        "background_density": bkg,
        "signal_entries": int(payload["signal"]["rows_after_cuts"]),
        "background_entries": int(payload["inclusive"]["rows_after_cuts"]),
        "auc": auc_from_probs(sig_prob, bkg_prob),
        "fpr": fpr,
        "tpr": tpr,
        "source": str(PP_SUMMARY),
        "note": "current-IAN baseV3E validation source, raw inclusive MC, no NPB cut, 15 < E_T < 35 GeV",
    }


def load_auau_raw() -> dict[str, np.ndarray]:
    paths = sorted(AUAU_SCORE_CACHE_DIR.glob("score_cache_*.npz"))
    if not paths:
        raise SystemExit(f"No THE-57 score caches found under {AUAU_SCORE_CACHE_DIR}")
    chunks: dict[str, list[np.ndarray]] = {k: [] for k in ("score", "is_signal", "centrality", "cluster_Et", "weight")}
    for path in paths:
        with np.load(path, allow_pickle=True) as data:
            score = np.asarray(data[SCORE_KEY], dtype=float)
            cent = np.asarray(data["centrality"], dtype=float)
            pt = np.asarray(data["cluster_Et"], dtype=float)
            is_sig = np.asarray(data["is_signal"], dtype=bool)
            weight = np.asarray(data["event_weight"], dtype=float) if "event_weight" in data.files else np.ones_like(score)
            keep = np.isfinite(score) & np.isfinite(cent) & np.isfinite(pt) & (pt >= 15.0) & (pt < 35.0) & (cent >= 0.0) & (cent < 80.0)
            chunks["score"].append(score[keep])
            chunks["is_signal"].append(is_sig[keep])
            chunks["centrality"].append(cent[keep])
            chunks["cluster_Et"].append(pt[keep])
            chunks["weight"].append(weight[keep])
    return {k: np.concatenate(v) for k, v in chunks.items()}


def make_density(scores: np.ndarray, weights: np.ndarray, edges: np.ndarray) -> np.ndarray:
    counts, _ = np.histogram(scores, bins=edges, weights=weights)
    widths = np.diff(edges)
    integral = float(np.sum(counts))
    if integral <= 0:
        raise ValueError("empty score histogram")
    return counts / integral / widths


def load_auau_blocks(edges: np.ndarray) -> list[dict[str, object]]:
    raw = load_auau_raw()
    out: list[dict[str, object]] = []
    for key, label, title, clo, chi in CENTRALITIES:
        in_cent = (raw["centrality"] >= clo) & (raw["centrality"] < chi)
        sig_mask = in_cent & raw["is_signal"]
        bkg_mask = in_cent & ~raw["is_signal"]
        sig_density = make_density(raw["score"][sig_mask], raw["weight"][sig_mask], edges)
        bkg_density = make_density(raw["score"][bkg_mask], raw["weight"][bkg_mask], edges)
        sig_prob = density_to_prob(edges, sig_density)
        bkg_prob = density_to_prob(edges, bkg_density)
        fpr, tpr = roc_from_probs(sig_prob, bkg_prob)
        auc_exact = weighted_auc_exact(raw["score"][sig_mask], raw["score"][bkg_mask], raw["weight"][sig_mask], raw["weight"][bkg_mask])
        out.append(
            {
                "key": key,
                "label": label,
                "title": title,
                "edges": edges,
                "signal_density": sig_density,
                "background_density": bkg_density,
                "signal_entries": int(np.count_nonzero(sig_mask)),
                "background_entries": int(np.count_nonzero(bkg_mask)),
                "signal_weight_sum": float(np.sum(raw["weight"][sig_mask])),
                "background_weight_sum": float(np.sum(raw["weight"][bkg_mask])),
                "auc": auc_exact,
                "auc_binned": auc_from_probs(sig_prob, bkg_prob),
                "fpr": fpr,
                "tpr": tpr,
                "source": str(AUAU_SCORE_CACHE_DIR),
                "note": "THE-57 default 14-feature baseline validation score caches, THE-58 cleaned, 15 < E_T < 35 GeV",
            }
        )
    return out


def step(ax: plt.Axes, edges: np.ndarray, density: np.ndarray, *, color: str, label: str, ls: str = "-", lw: float = 2.55, alpha: float = 1.0) -> None:
    y = np.r_[density, density[-1]]
    ax.step(edges, y, where="post", color=color, lw=lw, ls=ls, label=label, alpha=alpha)


def fmt_count(n: int) -> str:
    if n >= 1_000_000:
        return f"{n / 1_000_000:.2f}M"
    if n >= 1_000:
        return f"{n / 1_000:.0f}k"
    return str(n)


def decorate_axis(ax: plt.Axes) -> None:
    ax.grid(True, color=GRID, lw=0.75, alpha=0.75, which="both")
    ax.tick_params(direction="in", top=True, right=True, length=5)
    for spine in ax.spines.values():
        spine.set_color(INK)
        spine.set_linewidth(1.0)


def add_plot_label(ax: plt.Axes) -> None:
    sphenix = TextArea("sPHENIX", textprops={"fontsize": 11.4, "fontstyle": "italic", "fontweight": "bold", "fontfamily": "Times New Roman", "color": INK})
    internal = TextArea(" Internal", textprops={"fontsize": 11.4, "fontfamily": "Times New Roman", "color": INK})
    packed_label = HPacker(children=[sphenix, internal], align="baseline", pad=0, sep=0)
    anchored_label = AnchoredOffsetbox(loc="lower right", child=packed_label, frameon=False, pad=0, borderpad=0, bbox_to_anchor=(0.995, 1.006), bbox_transform=ax.transAxes)
    ax.add_artist(anchored_label)


def add_roc_readout(ax: plt.Axes, *, pp_auc: float, auau_auc: float, central_gap: float, closure_fraction: float | None) -> None:
    auc_gap = pp_auc - auau_auc
    lines = [f"AUC gap to pp: {auc_gap:.3f}"]
    if closure_fraction is not None:
        lines.append(f"{closure_fraction:.0%} of central gap recovered")
    else:
        lines.append("largest ranking penalty")
    ax.text(
        0.965,
        0.080,
        "\n".join(lines),
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=12.0,
        color=INK,
        linespacing=1.12,
        bbox={"boxstyle": "round,pad=0.26,rounding_size=0.02", "facecolor": "white", "edgecolor": "#cbd5e1", "linewidth": 0.8, "alpha": 0.94},
    )


def make_slide(pp: dict[str, object], auau_blocks: list[dict[str, object]]) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "font.family": ["Times New Roman", "Times", "DejaVu Serif"],
            "axes.labelsize": 15.5,
            "axes.titlesize": 20,
            "xtick.labelsize": 12.8,
            "ytick.labelsize": 12.8,
            "legend.fontsize": 14.5,
        }
    )

    fig = plt.figure(figsize=(16, 9), dpi=160)
    fig.patch.set_facecolor("white")
    central = auau_blocks[0]
    peripheral = auau_blocks[1]
    central_gap = float(pp["auc"]) - float(central["auc"])
    peripheral_gap = float(pp["auc"]) - float(peripheral["auc"])
    closure_fraction = (float(peripheral["auc"]) - float(central["auc"])) / central_gap

    fig.text(0.055, 0.935, "Peripheral Au+Au BDT separation moves back toward the pp baseline", fontsize=28.0, weight="bold", color=INK)
    insight_box = matplotlib.patches.FancyBboxPatch((0.055, 0.833), 0.890, 0.068, boxstyle="round,pad=0.008,rounding_size=0.004", transform=fig.transFigure, facecolor="#f8fbff", edgecolor="#b9cce8", linewidth=1.1)
    fig.add_artist(insight_box)
    insight_lead = TextArea("Expected centrality trend:", textprops={"fontsize": 17.0, "fontweight": "bold", "fontfamily": "Times New Roman", "color": INK})
    insight_body = TextArea(f" less underlying event \u2192 more pp-like; 50-80% is {peripheral_gap:.3f} AUC from pp, vs 0-20% at {central_gap:.3f}.", textprops={"fontsize": 17.0, "fontfamily": "Times New Roman", "color": INK})
    insight_text = HPacker(children=[insight_lead, insight_body], align="baseline", pad=0, sep=0)
    fig.add_artist(AnchoredOffsetbox(loc="center left", child=insight_text, frameon=False, pad=0, borderpad=0, bbox_to_anchor=(0.073, 0.867), bbox_transform=fig.transFigure))
    legend_handles = [
        matplotlib.lines.Line2D([0], [0], color=SIGNAL_RED, lw=3.0, ls="-", label="pp signal"),
        matplotlib.lines.Line2D([0], [0], color=BACKGROUND_BLUE, lw=3.0, ls="-", label="pp inclusive"),
        matplotlib.lines.Line2D([0], [0], color=SIGNAL_RED, lw=3.2, ls="--", label="Au+Au signal"),
        matplotlib.lines.Line2D([0], [0], color=BACKGROUND_BLUE, lw=3.2, ls="--", label="Au+Au background"),
    ]
    fig.legend(handles=legend_handles, loc="upper center", bbox_to_anchor=(0.530, 0.810), ncol=4, frameon=False, columnspacing=1.9, handlelength=2.8, fontsize=14.8)

    card_specs = [
        {
            "block": central,
            "x": 0.055,
            "edge": CENTRAL_EDGE,
            "header": CENTRAL_HEADER,
            "title": "0-20% central Au+Au vs pp",
            "subtitle": f"AUC {central['auc']:.3f} | gap to pp {central_gap:.3f} | S/B {fmt_count(central['signal_entries'])}/{fmt_count(central['background_entries'])}",
        },
        {
            "block": peripheral,
            "x": 0.525,
            "edge": PERIPHERAL_EDGE,
            "header": PERIPHERAL_HEADER,
            "title": "50-80% peripheral Au+Au vs pp",
            "subtitle": f"AUC {peripheral['auc']:.3f} | gap to pp {peripheral_gap:.3f} | {closure_fraction:.0%} of central-to-pp gap recovered",
        },
    ]
    for spec in card_specs:
        fig.add_artist(matplotlib.patches.FancyBboxPatch((spec["x"], 0.080), 0.420, 0.665, boxstyle="round,pad=0.010,rounding_size=0.006", transform=fig.transFigure, facecolor="white", edgecolor=spec["edge"], linewidth=1.5, zorder=-10))
        fig.add_artist(matplotlib.patches.FancyBboxPatch((spec["x"], 0.668), 0.420, 0.078, boxstyle="round,pad=0.010,rounding_size=0.006", transform=fig.transFigure, facecolor=spec["header"], edgecolor=spec["edge"], linewidth=1.4, zorder=-9))
        fig.text(spec["x"] + 0.020, 0.719, spec["title"], fontsize=18.0, weight="bold", color=INK, va="center")
        fig.text(spec["x"] + 0.020, 0.686, spec["subtitle"], fontsize=14.0, color=MUTED, va="center")

    score_axes = [fig.add_axes([0.102, 0.404, 0.340, 0.198]), fig.add_axes([0.572, 0.404, 0.340, 0.198])]
    roc_axes = [fig.add_axes([0.102, 0.151, 0.340, 0.205]), fig.add_axes([0.572, 0.151, 0.340, 0.205])]
    for x in (0.102, 0.572):
        fig.text(x, 0.614, "Score-shape overlay (BDT score)", fontsize=15.0, weight="bold", color=INK)
        fig.text(x, 0.365, "ROC: ranking performance", fontsize=15.0, weight="bold", color=INK)

    all_positive: list[float] = []
    for arr in [pp["signal_density"], pp["background_density"]]:
        vals = np.asarray(arr, dtype=float)
        all_positive.extend(vals[vals > 0])
    for block in auau_blocks:
        for arr in [block["signal_density"], block["background_density"]]:
            vals = np.asarray(arr, dtype=float)
            all_positive.extend(vals[vals > 0])
    ymin = max(min(all_positive) * 0.65, 1.8e-3)
    ymax = max(all_positive) * 1.95

    for ax_s, ax_r, block in zip(score_axes, roc_axes, auau_blocks, strict=True):
        step(ax_s, pp["edges"], pp["signal_density"], color=SIGNAL_RED, label="pp signal", alpha=0.72, lw=2.7)
        step(ax_s, pp["edges"], pp["background_density"], color=BACKGROUND_BLUE, label="pp inclusive", alpha=0.72, lw=2.7)
        step(ax_s, block["edges"], block["signal_density"], color=SIGNAL_RED, label=f"Au+Au {block['label']} signal", ls="--", lw=3.2)
        step(ax_s, block["edges"], block["background_density"], color=BACKGROUND_BLUE, label=f"Au+Au {block['label']} background", ls="--", lw=3.2)
        ax_s.set_yscale("log")
        ax_s.set_xlim(0.0, 1.0)
        ax_s.set_ylim(ymin, ymax)
        ax_s.set_ylabel("Area density (log)", labelpad=9)
        ax_s.tick_params(labelbottom=False)
        add_plot_label(ax_s)
        decorate_axis(ax_s)

        ax_r.plot(pp["fpr"], pp["tpr"], color=INK, lw=2.8, alpha=0.72, label=f"pp AUC {pp['auc']:.3f}")
        ax_r.plot(block["fpr"], block["tpr"], color=BACKGROUND_BLUE, lw=3.2, ls="--", label=f"Au+Au {block['label']} AUC {block['auc']:.3f}")
        ax_r.plot([0, 1], [0, 1], color="#9ca3af", lw=1.4, ls=":")
        ax_r.set_xlim(0.0, 1.0)
        ax_r.set_ylim(0.0, 1.0)
        ax_r.set_ylabel("Signal efficiency", labelpad=9)
        add_plot_label(ax_r)
        block_gap = float(pp["auc"]) - float(block["auc"])
        block_recovery = None if str(block["label"]) == "0-20%" else (central_gap - block_gap) / central_gap
        ax_r.legend(loc="lower right", bbox_to_anchor=(0.965, 0.405), frameon=False, fontsize=10.6, handlelength=2.3, borderaxespad=0.0)
        add_roc_readout(ax_r, pp_auc=float(pp["auc"]), auau_auc=float(block["auc"]), central_gap=central_gap, closure_fraction=block_recovery)
        decorate_axis(ax_r)

    for ax in roc_axes:
        ax.set_xlabel("Background efficiency")

    fig.savefig(OUT_PNG)
    plt.close(fig)


def write_csv(pp: dict[str, object], auau_blocks: list[dict[str, object]]) -> None:
    rows: list[dict[str, object]] = []

    def append_block(sample: str, centrality: str, block: dict[str, object]) -> None:
        edges = np.asarray(block["edges"], dtype=float)
        for cls, values in [("signal", np.asarray(block["signal_density"], dtype=float)), ("background", np.asarray(block["background_density"], dtype=float))]:
            for lo, hi, density in zip(edges[:-1], edges[1:], values, strict=True):
                rows.append(
                    {
                        "sample": sample,
                        "centrality": centrality,
                        "class": cls,
                        "bin_low": f"{lo:.6g}",
                        "bin_high": f"{hi:.6g}",
                        "density": f"{float(density):.9g}",
                        "auc": f"{float(block['auc']):.9g}",
                        "auc_binned": f"{float(block.get('auc_binned', block['auc'])):.9g}",
                        "signal_entries": block["signal_entries"],
                        "background_entries": block["background_entries"],
                    }
                )

    append_block("pp", "none", pp)
    for block in auau_blocks:
        append_block("AuAu", str(block["label"]), block)
    with OUT_CSV.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def write_manifest(pp: dict[str, object], auau_blocks: list[dict[str, object]]) -> None:
    OUT_JSON.write_text(
        json.dumps(
            {
                "schema": "THE57_PP_AUAU_CENTRAL_PERIPHERAL_BDT_BRIDGE_V1",
                "output_png": str(OUT_PNG),
                "speaker_script": str(OUT_SCRIPT),
                "histogram_csv": str(OUT_CSV),
                "canvas_px": [2560, 1440],
                "deck_context": "Working-point deck slide 9 regeneration, not the THE-57 WP80 heatmap slide.",
                "pp": {
                    "source": pp["source"],
                    "signal_entries": pp["signal_entries"],
                    "background_entries": pp["background_entries"],
                    "auc_binned_from_display_hist": pp["auc"],
                    "note": pp["note"],
                },
                "auau": [
                    {
                        "source": block["source"],
                        "product": PRODUCT,
                        "centrality": block["label"],
                        "signal_entries": block["signal_entries"],
                        "background_entries": block["background_entries"],
                        "signal_weight_sum": block["signal_weight_sum"],
                        "background_weight_sum": block["background_weight_sum"],
                        "auc_exact_weighted": block["auc"],
                        "auc_binned_recomputed": block["auc_binned"],
                        "note": block["note"],
                    }
                    for block in auau_blocks
                ],
                "interpretation": [
                    "This slide is the working-point deck slide-9 bridge from pp/PPG12 consistency to AuAu centrality dependence.",
                    "The pp reference is the existing validated current-IAN pp BDT source used by the original slide.",
                    "The AuAu comparison is regenerated from the THE-57 default 14-feature baseline validation score caches.",
                    "Central AuAu remains farther from pp; peripheral AuAu moves back toward the pp score/ROC baseline.",
                ],
                "google_slides_mutated": False,
            },
            indent=2,
        )
        + "\n"
    )


def write_speaker_script(pp: dict[str, object], auau_blocks: list[dict[str, object]]) -> None:
    central = auau_blocks[0]
    peripheral = auau_blocks[1]
    central_gap = float(pp["auc"]) - float(central["auc"])
    peripheral_gap = float(pp["auc"]) - float(peripheral["auc"])
    closure_fraction = (float(peripheral["auc"]) - float(central["auc"])) / central_gap
    OUT_SCRIPT.write_text(
        "\n".join(
            [
                "# Working-point slide 9 script - pp to AuAu BDT score bridge",
                "",
                "This is the same slide-9 comparison layout from the working-point deck, regenerated with the current THE-57 default AuAu BDT validation output.",
                "",
                "The pp curves are the existing current-IAN pp reference used by the original slide. The AuAu curves are the default 14-feature baseline with THE-58 cleaning, evaluated on validation score-cache rows in the same 15 to 35 GeV candidate window.",
                "",
                f"The central bin is still the harder environment: AUC {central['auc']:.3f}, which is {central_gap:.3f} below the pp reference AUC of {pp['auc']:.3f}. The peripheral bin has AUC {peripheral['auc']:.3f}, leaving a {peripheral_gap:.3f} gap to pp and recovering about {closure_fraction:.0%} of the central-to-pp gap.",
                "",
                "The point is the same as before but now tied to the corrected default model: the BDT ranking becomes more pp-like as the AuAu underlying event gets smaller.",
                "",
            ]
        )
    )


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    pp = load_pp()
    auau_blocks = load_auau_blocks(np.asarray(pp["edges"], dtype=float))
    make_slide(pp, auau_blocks)
    write_csv(pp, auau_blocks)
    write_manifest(pp, auau_blocks)
    write_speaker_script(pp, auau_blocks)
    print(OUT_PNG)
    print(OUT_JSON)
    print(OUT_CSV)
    print(OUT_SCRIPT)


if __name__ == "__main__":
    main()
