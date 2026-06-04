#!/usr/bin/env python3
"""Build a pp-to-AuAu BDT score/ROC bridge slide candidate."""

from __future__ import annotations

import csv
import json
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.offsetbox import AnchoredOffsetbox, HPacker, TextArea  # noqa: E402
import numpy as np  # noqa: E402


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
AUAU_COMPACT_HIST = (
    REPO
    / "dataOutput/auauTightBDTValidation/THE8_branchA_ladder_scorecache_fullstat_20260527"
    / "the8_branch_a_ladder_compact_score_histograms.json"
)
AUAU_PANEL_SUMMARY = (
    REPO
    / "dataOutput/auauTightBDTValidation/THE8_branchA_ladder_scorecache_fullstat_20260527"
    / "slide23_candidates/the8_branchA_ladder_score_separation_slide23_style.json"
)

AUAU_BRANCH = "Jet12+20+30+40"
CENTRALITIES = [
    ("0_20", "0-20%", "Central Au+Au"),
    ("50_80", "50-80%", "Peripheral Au+Au"),
]

OUT_DIR = REPO / "dataOutput/slides/wp_gammajets_6_3_26/pp_auau_bridge_20260603"
OUT_PNG = OUT_DIR / "pp_auau_central_peripheral_bdt_bridge.png"
OUT_JSON = OUT_DIR / "pp_auau_central_peripheral_bdt_bridge_manifest.json"
OUT_CSV = OUT_DIR / "pp_auau_central_peripheral_bdt_bridge_histograms.csv"
OUT_SCRIPT = OUT_DIR / "pp_auau_central_peripheral_bdt_bridge_script.md"

INK = "#111827"
MUTED = "#4b5563"
GRID = "#d9dee8"
SIGNAL_RED = "#d62728"
BACKGROUND_BLUE = "#1f77b4"
GREEN_BG = "#eef8f1"
GREEN_EDGE = "#8bd0a2"
BLUE_BG = "#eef5ff"
BLUE_EDGE = "#9dbcf5"
AMBER_BG = "#fff7e8"
AMBER_EDGE = "#e8b45b"
CENTRAL_EDGE = "#d95f5f"
CENTRAL_HEADER = "#fff0f0"
PERIPHERAL_EDGE = "#2f8f67"
PERIPHERAL_HEADER = "#edf8f0"
PP_HEADER = "#f3f4f6"


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
        "note": "current-IAN baseV3E, raw inclusive MC, no NPB cut, 15 < E_T < 35 GeV",
    }


def load_auau_blocks() -> list[dict[str, object]]:
    compact = json.loads(AUAU_COMPACT_HIST.read_text())
    panel = json.loads(AUAU_PANEL_SUMMARY.read_text())
    branch = next((item for item in compact["branches"] if item["label"] == AUAU_BRANCH), None)
    if branch is None:
        raise RuntimeError(f"No compact AuAu branch found for {AUAU_BRANCH}")

    out: list[dict[str, object]] = []
    for key, label, title in CENTRALITIES:
        cent = branch["by_centrality"].get(key)
        if cent is None:
            raise RuntimeError(f"No AuAu centrality block found for {key}")
        row = next(
            (
                item
                for item in panel["rows"]
                if item["branch"] == AUAU_BRANCH and item["centrality"] == label
            ),
            None,
        )
        if row is None:
            raise RuntimeError(f"No AuAu panel summary row found for {AUAU_BRANCH} {label}")
        edges = np.asarray(branch["bin_edges"], dtype=float)
        sig = np.asarray(cent["signal"]["density"], dtype=float)
        bkg = np.asarray(cent["background"]["density"], dtype=float)
        sig_prob = density_to_prob(edges, sig)
        bkg_prob = density_to_prob(edges, bkg)
        fpr, tpr = roc_from_probs(sig_prob, bkg_prob)
        out.append(
            {
                "label": label,
                "title": title,
                "edges": edges,
                "signal_density": sig,
                "background_density": bkg,
                "signal_entries": int(cent["signal"]["entries"]),
                "background_entries": int(cent["background"]["entries"]),
                "auc": float(row["auc_binned"]),
                "auc_recomputed": auc_from_probs(sig_prob, bkg_prob),
                "fpr": fpr,
                "tpr": tpr,
                "source": str(AUAU_COMPACT_HIST),
                "summary_source": str(AUAU_PANEL_SUMMARY),
                "note": f"Branch A {AUAU_BRANCH} global no-isolation BDT, 15 < E_T < 35 GeV",
            }
        )
    return out


def step(
    ax: plt.Axes,
    edges: np.ndarray,
    density: np.ndarray,
    *,
    color: str,
    label: str,
    ls: str = "-",
    lw: float = 2.55,
    alpha: float = 1.0,
) -> None:
    y = np.r_[density, density[-1]]
    ax.step(edges, y, where="post", color=color, lw=lw, ls=ls, label=label, alpha=alpha)


def fmt_count(n: int) -> str:
    if n >= 1_000_000:
        return f"{n / 1_000_000:.2f}M"
    if n >= 1_000:
        return f"{n / 1_000:.0f}k"
    return str(n)


def add_box(
    fig: plt.Figure,
    x: float,
    y: float,
    w: float,
    h: float,
    title: str,
    body: str,
    fc: str,
    ec: str,
    *,
    wrap: int,
    body_size: float = 14.0,
    title_size: float = 17.0,
    xpad: float = 0.014,
    title_offset: float = 0.032,
    body_offset: float = 0.068,
    linespacing: float = 1.18,
) -> None:
    patch = matplotlib.patches.FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.010,rounding_size=0.006",
        transform=fig.transFigure,
        facecolor=fc,
        edgecolor=ec,
        linewidth=1.1,
    )
    fig.add_artist(patch)
    fig.text(x + xpad, y + h - title_offset, title, fontsize=title_size, weight="bold", color=INK, va="top")
    wrapped = "\n".join(textwrap.fill(part, width=wrap) for part in body.split("\n"))
    fig.text(
        x + xpad,
        y + h - body_offset,
        wrapped,
        fontsize=body_size,
        color="#253247",
        va="top",
        linespacing=linespacing,
    )


def decorate_axis(ax: plt.Axes) -> None:
    ax.grid(True, color=GRID, lw=0.75, alpha=0.75, which="both")
    ax.tick_params(direction="in", top=True, right=True, length=5)
    for spine in ax.spines.values():
        spine.set_color(INK)
        spine.set_linewidth(1.0)


def add_plot_label(ax: plt.Axes, *, loc: str) -> None:
    sphenix = TextArea(
        "sPHENIX",
        textprops={
            "fontsize": 11.4,
            "fontstyle": "italic",
            "fontweight": "bold",
            "fontfamily": "Times New Roman",
            "color": INK,
        },
    )
    internal = TextArea(
        " Internal",
        textprops={
            "fontsize": 11.4,
            "fontfamily": "Times New Roman",
            "color": INK,
        },
    )
    packed_label = HPacker(children=[sphenix, internal], align="baseline", pad=0, sep=0)
    anchored_label = AnchoredOffsetbox(
        loc="lower right",
        child=packed_label,
        frameon=False,
        pad=0,
        borderpad=0,
        bbox_to_anchor=(0.995, 1.006),
        bbox_transform=ax.transAxes,
    )
    ax.add_artist(anchored_label)


def add_roc_readout(
    ax: plt.Axes,
    *,
    pp_auc: float,
    auau_auc: float,
    central_gap: float,
    closure_fraction: float | None,
) -> None:
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
        bbox={
            "boxstyle": "round,pad=0.26,rounding_size=0.02",
            "facecolor": "white",
            "edgecolor": "#cbd5e1",
            "linewidth": 0.8,
            "alpha": 0.94,
        },
    )
    if central_gap <= 0:
        return


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

    fig.text(
        0.055,
        0.935,
        "Peripheral Au+Au BDT separation moves back toward the pp baseline",
        fontsize=28.0,
        weight="bold",
        color=INK,
    )
    insight_box = matplotlib.patches.FancyBboxPatch(
        (0.055, 0.833),
        0.890,
        0.068,
        boxstyle="round,pad=0.008,rounding_size=0.004",
        transform=fig.transFigure,
        facecolor="#f8fbff",
        edgecolor="#b9cce8",
        linewidth=1.1,
    )
    fig.add_artist(insight_box)
    insight_lead = TextArea(
        "Expected centrality trend:",
        textprops={
            "fontsize": 17.0,
            "fontweight": "bold",
            "fontfamily": "Times New Roman",
            "color": INK,
        },
    )
    insight_body = TextArea(
        f" less underlying event \u2192 more pp-like; 50-80% is {peripheral_gap:.3f} AUC from pp, vs 0-20% at {central_gap:.3f}.",
        textprops={
            "fontsize": 17.0,
            "fontfamily": "Times New Roman",
            "color": INK,
        },
    )
    insight_text = HPacker(children=[insight_lead, insight_body], align="baseline", pad=0, sep=0)
    insight_anchor = AnchoredOffsetbox(
        loc="center left",
        child=insight_text,
        frameon=False,
        pad=0,
        borderpad=0,
        bbox_to_anchor=(0.073, 0.867),
        bbox_transform=fig.transFigure,
    )
    fig.add_artist(insight_anchor)
    legend_handles = [
        matplotlib.lines.Line2D([0], [0], color=SIGNAL_RED, lw=3.0, ls="-", label="pp signal"),
        matplotlib.lines.Line2D([0], [0], color=BACKGROUND_BLUE, lw=3.0, ls="-", label="pp inclusive"),
        matplotlib.lines.Line2D([0], [0], color=SIGNAL_RED, lw=3.2, ls="--", label="Au+Au signal"),
        matplotlib.lines.Line2D([0], [0], color=BACKGROUND_BLUE, lw=3.2, ls="--", label="Au+Au background"),
    ]
    fig.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.530, 0.810),
        ncol=4,
        frameon=False,
        columnspacing=1.9,
        handlelength=2.8,
        fontsize=14.8,
    )

    card_specs = [
        {
            "block": central,
            "x": 0.055,
            "edge": CENTRAL_EDGE,
            "header": CENTRAL_HEADER,
            "title": "0-20% central Au+Au vs pp",
            "subtitle": (
                f"AUC {central['auc']:.3f} | gap to pp {central_gap:.3f} | "
                f"S/B {fmt_count(central['signal_entries'])}/{fmt_count(central['background_entries'])}"
            ),
        },
        {
            "block": peripheral,
            "x": 0.525,
            "edge": PERIPHERAL_EDGE,
            "header": PERIPHERAL_HEADER,
            "title": "50-80% peripheral Au+Au vs pp",
            "subtitle": (
                f"AUC {peripheral['auc']:.3f} | gap to pp {peripheral_gap:.3f} | "
                f"{closure_fraction:.0%} of central-to-pp gap recovered"
            ),
        },
    ]
    for spec in card_specs:
        card = matplotlib.patches.FancyBboxPatch(
            (spec["x"], 0.080),
            0.420,
            0.665,
            boxstyle="round,pad=0.010,rounding_size=0.006",
            transform=fig.transFigure,
            facecolor="white",
            edgecolor=spec["edge"],
            linewidth=1.5,
            zorder=-10,
        )
        header = matplotlib.patches.FancyBboxPatch(
            (spec["x"], 0.668),
            0.420,
            0.078,
            boxstyle="round,pad=0.010,rounding_size=0.006",
            transform=fig.transFigure,
            facecolor=spec["header"],
            edgecolor=spec["edge"],
            linewidth=1.4,
            zorder=-9,
        )
        fig.add_artist(card)
        fig.add_artist(header)
        fig.text(spec["x"] + 0.020, 0.719, spec["title"], fontsize=18.0, weight="bold", color=INK, va="center")
        fig.text(spec["x"] + 0.020, 0.686, spec["subtitle"], fontsize=14.0, color=MUTED, va="center")

    score_axes = [
        fig.add_axes([0.102, 0.404, 0.340, 0.198]),
        fig.add_axes([0.572, 0.404, 0.340, 0.198]),
    ]
    roc_axes = [
        fig.add_axes([0.102, 0.151, 0.340, 0.205]),
        fig.add_axes([0.572, 0.151, 0.340, 0.205]),
    ]

    fig.text(0.102, 0.614, "Score-shape overlay (BDT score)", fontsize=15.0, weight="bold", color=INK)
    fig.text(0.102, 0.365, "ROC: ranking performance", fontsize=15.0, weight="bold", color=INK)
    fig.text(0.572, 0.614, "Score-shape overlay (BDT score)", fontsize=15.0, weight="bold", color=INK)
    fig.text(0.572, 0.365, "ROC: ranking performance", fontsize=15.0, weight="bold", color=INK)

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
        step(
            ax_s,
            pp["edges"],
            pp["background_density"],
            color=BACKGROUND_BLUE,
            label="pp inclusive",
            alpha=0.72,
            lw=2.7,
        )
        step(
            ax_s,
            block["edges"],
            block["signal_density"],
            color=SIGNAL_RED,
            label=f"Au+Au {block['label']} signal",
            ls="--",
            lw=3.2,
        )
        step(
            ax_s,
            block["edges"],
            block["background_density"],
            color=BACKGROUND_BLUE,
            label=f"Au+Au {block['label']} background",
            ls="--",
            lw=3.2,
        )
        ax_s.set_yscale("log")
        ax_s.set_xlim(0.0, 1.0)
        ax_s.set_ylim(ymin, ymax)
        ax_s.set_ylabel("Area density (log)", labelpad=9)
        ax_s.tick_params(labelbottom=False)
        add_plot_label(ax_s, loc="upper_left")
        decorate_axis(ax_s)

        ax_r.plot(pp["fpr"], pp["tpr"], color=INK, lw=2.8, alpha=0.72, label=f"pp AUC {pp['auc']:.3f}")
        ax_r.plot(
            block["fpr"],
            block["tpr"],
            color=BACKGROUND_BLUE,
            lw=3.2,
            ls="--",
            label=f"Au+Au {block['label']} AUC {block['auc']:.3f}",
        )
        ax_r.plot([0, 1], [0, 1], color="#9ca3af", lw=1.4, ls=":")
        ax_r.set_xlim(0.0, 1.0)
        ax_r.set_ylim(0.0, 1.0)
        ax_r.set_ylabel("Signal efficiency", labelpad=9)
        add_plot_label(ax_r, loc="lower_right")
        block_gap = float(pp["auc"]) - float(block["auc"])
        block_recovery = None if str(block["label"]) == "0-20%" else (central_gap - block_gap) / central_gap
        ax_r.legend(
            loc="lower right",
            bbox_to_anchor=(0.965, 0.405),
            frameon=False,
            fontsize=10.6,
            handlelength=2.3,
            borderaxespad=0.0,
        )
        add_roc_readout(
            ax_r,
            pp_auc=float(pp["auc"]),
            auau_auc=float(block["auc"]),
            central_gap=central_gap,
            closure_fraction=block_recovery,
        )
        decorate_axis(ax_r)

    for ax in roc_axes:
        ax.set_xlabel("Background efficiency")

    fig.savefig(OUT_PNG)
    plt.close(fig)


def write_csv(pp: dict[str, object], auau_blocks: list[dict[str, object]]) -> None:
    rows: list[dict[str, object]] = []

    def append_block(sample: str, centrality: str, block: dict[str, object]) -> None:
        edges = np.asarray(block["edges"], dtype=float)
        for cls, values in [
            ("signal", np.asarray(block["signal_density"], dtype=float)),
            ("background", np.asarray(block["background_density"], dtype=float)),
        ]:
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
                "output_png": str(OUT_PNG),
                "speaker_script": str(OUT_SCRIPT),
                "histogram_csv": str(OUT_CSV),
                "canvas_px": [2560, 1440],
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
                        "summary_source": block["summary_source"],
                        "branch": AUAU_BRANCH,
                        "centrality": block["label"],
                        "signal_entries": block["signal_entries"],
                        "background_entries": block["background_entries"],
                        "auc_from_panel_summary": block["auc"],
                        "auc_binned_recomputed": block["auc_recomputed"],
                        "note": block["note"],
                    }
                    for block in auau_blocks
                ],
                "interpretation": [
                    "This slide is a two-column bridge from pp/PPG12 consistency to Au+Au centrality dependence.",
                    "The pp reference uses raw inclusive MC with no NPB cut in the same 15-35 GeV window.",
                    "The Au+Au comparison uses the highest-AUC Branch A Jet12+20+30+40 global no-isolation BDT validation product.",
                    "Central Au+Au is farther from the pp ROC baseline; peripheral Au+Au recovers most of the central-to-pp AUC gap.",
                ],
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
                "# WP GammaJets Slide Script - pp to Au+Au BDT score bridge",
                "",
                "Now that the pp baseline is internally consistent with the PPG12 reference, the next thing I want to check is whether the Au+Au classifier behavior moves in the direction we expect with centrality.",
                "",
                f"The comparison is organized left to right. The left card is the 0 to 20 percent central Au+Au bin overlaid against the pp reference, and the right card is the 50 to 80 percent peripheral bin overlaid against the same pp reference. In both cards the top plot is the signal and background score shape, and the bottom plot is the ROC curve.",
                "",
                f"The central bin is visibly the harder environment. Its AUC is about {central['auc']:.3f}, which is {central_gap:.3f} below the pp reference AUC of {pp['auc']:.3f}. In the peripheral bin the Au+Au ROC moves much closer to pp, with an AUC of about {peripheral['auc']:.3f}; that leaves only a {peripheral_gap:.3f} AUC gap to pp and recovers about {closure_fraction:.0%} of the central-to-pp gap.",
                "",
                "So the point is not just that peripheral looks better. The point is that the centrality trend is coherent with the underlying-event picture: as the environment becomes less central, the photon-ID ranking problem becomes more pp-like.",
                "",
                "That is the bridge into the working-point decision: after pp closure, the remaining question is how much centrality dependence we need to carry into the Au+Au fake-rate, composition, isolation, and ABCD closure checks.",
                "",
            ]
        )
    )


def main() -> None:
    pp = load_pp()
    auau_blocks = load_auau_blocks()
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
