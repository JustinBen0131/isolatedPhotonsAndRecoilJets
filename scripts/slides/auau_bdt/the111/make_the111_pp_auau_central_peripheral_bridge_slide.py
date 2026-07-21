#!/usr/bin/env python3
"""Build the JSTG pp-to-AuAu BDT bridge slide from validated score products."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.offsetbox import AnchoredOffsetbox, HPacker, TextArea
from matplotlib.patches import FancyBboxPatch
import numpy as np


def find_repo() -> Path:
    here = Path(__file__).resolve()
    for parent in here.parents:
        if (parent / "AGENTS.md").exists() and (parent / "agent_context").exists():
            return parent
    raise RuntimeError("Could not resolve ThesisAnalysis repository root")


REPO = find_repo()
PP_SOURCE = (
    REPO
    / "dataOutput/ppPhotonMLPipeline/ppg12_basev3E_currentIAN_rawOverlayEnvFix_20260527_1630"
    / "validation/fullsim_shuhang_overlay_raw_inclusive_pt1535"
    / "pp_currentian_basev3e_bdt_score_overlay_vs_ppg12_pp_noCent_bdt_15_35_noNPB_eta0-pt3-cut0_summary.json"
)
THE111_SOURCE = (
    REPO
    / "dataOutput/slides/the45_jstg_20260720/the111_pp_bridge/source"
    / "the111_holdout_score_curves.json"
)
THE111_VALIDATION = (
    REPO
    / "dataOutput/auauTightBDTValidation/the111_combined_corrected_shower_ppg12_labels_20260719_1618"
    / "paired/the111_combined_corrected_shower_ppg12_labels_validation.json"
)
OUT_DIR = REPO / "dataOutput/slides/the45_jstg_20260720/the111_pp_bridge"
OUT_PNG = OUT_DIR / "the111_pp_auau_central_peripheral_bdt_bridge.png"
OUT_CSV = OUT_DIR / "the111_pp_auau_central_peripheral_bdt_bridge.csv"
OUT_MANIFEST = OUT_DIR / "the111_pp_auau_central_peripheral_bdt_bridge_manifest.json"
OUT_SCRIPT = OUT_DIR / "the111_pp_auau_central_peripheral_bdt_bridge_speaker_notes.md"

INK = "#152238"
MUTED = "#56657a"
GRID = "#d7dee8"
SIGNAL = "#c93b30"
BACKGROUND = "#2f67c7"
CENTRAL_EDGE = "#d56060"
CENTRAL_FILL = "#fff6f4"
PERIPHERAL_EDGE = "#2f8f68"
PERIPHERAL_FILL = "#f2faf5"
ACCENT = "#2474b8"
# Shared JSTG subtitle-arrowhead convention (THE-111 slides 11/12).
BLUE_BULLET = "#2468A8"
PP_AUC_EXPECTED = 0.934245502817205

# Text drawn after the italic bold "sPHENIX" tag on every plotted canvas.
# Default preserves the original slide exactly; --status-tag overrides it.
STATUS_TAG = " Simulation"

# Shaded header panel behind the subtitle block.  The subtitle text is drawn
# independently, so clearing this leaves the wording in place.
DRAW_HEADER_RIBBON = True


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def density_to_prob(edges: np.ndarray, density: np.ndarray) -> np.ndarray:
    prob = np.asarray(density, dtype=float) * np.diff(edges)
    total = float(np.sum(prob))
    if total <= 0:
        raise ValueError("Histogram has non-positive total probability")
    return prob / total


def roc_from_probabilities(signal: np.ndarray, background: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    return np.r_[0.0, np.cumsum(background[::-1])], np.r_[0.0, np.cumsum(signal[::-1])]


def auc_from_probabilities(signal: np.ndarray, background: np.ndarray) -> float:
    background_below = np.r_[0.0, np.cumsum(background)[:-1]]
    return float(np.sum(signal * background_below) + 0.5 * np.sum(signal * background))


def load_pp() -> dict[str, object]:
    payload = json.loads(PP_SOURCE.read_text())
    edges = np.asarray(payload["bins"], dtype=float)
    # The accepted pp summary stores unit-sum bin probabilities, whereas the
    # THE-111 extraction stores probability density.  Convert the pp arrays to
    # the same density convention before drawing the cross-system overlay.
    signal_prob = np.asarray(payload["this_analysis_signal_hist"], dtype=float)
    background_prob = np.asarray(payload["this_analysis_inclusive_hist"], dtype=float)
    signal_prob = signal_prob / np.sum(signal_prob)
    background_prob = background_prob / np.sum(background_prob)
    signal_density = signal_prob / np.diff(edges)
    background_density = background_prob / np.diff(edges)
    fpr, tpr = roc_from_probabilities(signal_prob, background_prob)
    auc = auc_from_probabilities(signal_prob, background_prob)
    if abs(auc - PP_AUC_EXPECTED) > 5e-4:
        raise RuntimeError(f"Unexpected pp AUC {auc:.8f}; expected {PP_AUC_EXPECTED:.8f}")
    return {
        "edges": edges,
        "signal_density": signal_density,
        "background_density": background_density,
        "fpr": fpr,
        "tpr": tpr,
        "auc": auc,
        "signal_rows": int(payload["signal"]["rows_after_cuts"]),
        "background_rows": int(payload["inclusive"]["rows_after_cuts"]),
    }


def load_auau() -> list[dict[str, object]]:
    payload = json.loads(THE111_SOURCE.read_text())
    if payload.get("schema") != "THE111_HOLDOUT_SCORE_CURVES_V1":
        raise RuntimeError("Unexpected THE-111 score-curve schema")
    edges = np.asarray(payload["bin_edges"], dtype=float)
    validation = json.loads(THE111_VALIDATION.read_text())
    validated = {
        tuple(block["centrality"]): block
        for block in validation["native"]["ppg12_variant"]["centrality_bins"]
    }
    out: list[dict[str, object]] = []
    for block in payload["blocks"]:
        key = tuple(block["centrality"])
        if key not in validated:
            raise RuntimeError(f"Missing validation block for centrality {key}")
        authoritative = validated[key]
        if abs(float(block["weighted_auc"]) - float(authoritative["weighted_auc"])) > 1e-12:
            raise RuntimeError(f"THE-111 AUC mismatch for centrality {key}")
        out.append(
            {
                "centrality": list(key),
                "label": f"{int(key[0])}\N{EN DASH}{int(key[1])}%",
                "edges": edges,
                "signal_density": np.asarray(block["signal_density"], dtype=float),
                "background_density": np.asarray(block["background_density"], dtype=float),
                "fpr": np.asarray(block["roc_fpr"], dtype=float),
                "tpr": np.asarray(block["roc_tpr"], dtype=float),
                "auc": float(authoritative["weighted_auc"]),
                "wp80_background_acceptance": float(authoritative["wp80"]["background_acceptance"]),
                "wp80_threshold": float(authoritative["wp80"]["threshold"]),
                "signal_rows": int(block["signal_rows"]),
                "background_rows": int(block["background_rows"]),
            }
        )
    return out


def step(ax: plt.Axes, edges: np.ndarray, density: np.ndarray, *, color: str, linestyle: str, linewidth: float, alpha: float = 1.0) -> None:
    ax.step(edges, np.r_[density, density[-1]], where="post", color=color, ls=linestyle, lw=linewidth, alpha=alpha)


def add_sphenix_label(ax: plt.Axes) -> None:
    experiment = TextArea(
        "sPHENIX",
        textprops={"fontsize": 10.2, "fontstyle": "italic", "fontweight": "bold", "fontfamily": "Times New Roman", "color": INK},
    )
    status = TextArea(
        STATUS_TAG,
        textprops={"fontsize": 10.2, "fontfamily": "Times New Roman", "color": INK},
    )
    label = HPacker(children=[experiment, status], align="baseline", pad=0, sep=0)
    ax.add_artist(
        AnchoredOffsetbox(
            loc="upper right",
            child=label,
            frameon=False,
            pad=0,
            borderpad=0,
            bbox_to_anchor=(0.985, 0.965),
            bbox_transform=ax.transAxes,
        )
    )


def decorate(ax: plt.Axes) -> None:
    ax.grid(True, color=GRID, lw=0.65, alpha=0.72, which="both")
    ax.tick_params(direction="in", top=True, right=True, length=4.5, width=0.9)
    for spine in ax.spines.values():
        spine.set_color(INK)
        spine.set_linewidth(0.95)


def make_slide(pp: dict[str, object], auau: list[dict[str, object]]) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "font.family": ["Times New Roman", "Times", "DejaVu Serif"],
            "axes.labelsize": 15.0,
            "xtick.labelsize": 12.3,
            "ytick.labelsize": 12.3,
        }
    )
    central, peripheral = auau
    pp_auc = float(pp["auc"])
    central_gap = pp_auc - float(central["auc"])
    peripheral_gap = pp_auc - float(peripheral["auc"])
    recovered = (float(peripheral["auc"]) - float(central["auc"])) / central_gap

    fig = plt.figure(figsize=(16, 9), dpi=160, facecolor="white")
    fig.text(
        0.050,
        0.936,
        "Au+Au vs pp photon BDT: 0–20% and 50–80% simulation validation",
        fontsize=29.5,
        weight="bold",
        color=INK,
        va="top",
    )

    if DRAW_HEADER_RIBBON:
        ribbon = FancyBboxPatch(
            (0.050, 0.828),
            0.900,
            0.070,
            transform=fig.transFigure,
            boxstyle="round,pad=0.006,rounding_size=0.004",
            facecolor="#f8fbff",
            edgecolor="#b7cce5",
            linewidth=1.05,
        )
        fig.add_artist(ribbon)
    # Two arrowhead bullets, centred in the band between the title ink and the
    # legend, using the shared JSTG convention (THE-111 slides 11/12).
    bullets = [
        (
            0.858,
            r"Au+Au is pp's 11 baseV3E features + $w_{\eta|\phi}^{3\times3}$ + centrality",
        ),
        (
            0.8055,
            r"Both pp/Au+Au trained using personal analysis code, 15–35 GeV, $|\eta|<0.7$ (not PPG12 BDT)",
        ),
    ]
    for bullet_y, bullet_text in bullets:
        fig.text(
            0.050,
            bullet_y,
            "▶",
            fontsize=15.0,
            color=BLUE_BULLET,
            ha="left",
            va="center",
            fontfamily="DejaVu Sans",
        )
        fig.text(
            0.071,
            bullet_y,
            bullet_text,
            fontsize=17.5,
            color=INK,
            va="center",
        )

    legend = [
        Line2D([0], [0], color=SIGNAL, lw=3.0, ls="-", label="pp signal"),
        Line2D([0], [0], color=BACKGROUND, lw=3.0, ls="-", label="pp inclusive"),
        Line2D([0], [0], color=SIGNAL, lw=3.1, ls="--", label="Au+Au signal"),
        Line2D([0], [0], color=BACKGROUND, lw=3.1, ls="--", label="Au+Au background"),
    ]
    fig.legend(
        handles=legend,
        loc="upper center",
        bbox_to_anchor=(0.535, 0.792),
        ncol=4,
        frameon=False,
        fontsize=14.4,
        handlelength=2.8,
        columnspacing=1.85,
    )

    cards = [
        (central, 0.050, CENTRAL_EDGE, CENTRAL_FILL, "Central Au+Au"),
        (peripheral, 0.525, PERIPHERAL_EDGE, PERIPHERAL_FILL, "Peripheral Au+Au"),
    ]
    for block, x, edge, fill, label in cards:
        fig.add_artist(
            FancyBboxPatch(
                (x, 0.067),
                0.425,
                0.669,
                transform=fig.transFigure,
                boxstyle="round,pad=0.009,rounding_size=0.006",
                facecolor="white",
                edgecolor=edge,
                linewidth=1.45,
                zorder=-10,
            )
        )
        # Title band holds the label alone; the AUC / WP80 numbers already
        # appear in the ROC legend and gap box below.  Band shortened and the
        # label enlarged, centred with a small even buffer top and bottom.
        band_y0, band_h = 0.674, 0.062
        fig.add_artist(
            FancyBboxPatch(
                (x, band_y0),
                0.425,
                band_h,
                transform=fig.transFigure,
                boxstyle="round,pad=0.009,rounding_size=0.006",
                facecolor=fill,
                edgecolor=edge,
                linewidth=1.25,
                zorder=-9,
            )
        )
        fig.text(
            x + 0.020,
            band_y0 + band_h / 2.0,
            f"{block['label']} {label} vs pp",
            fontsize=21.5,
            weight="bold",
            color=INK,
            va="center",
        )

    score_axes = [fig.add_axes([0.099, 0.437, 0.345, 0.160]), fig.add_axes([0.574, 0.437, 0.345, 0.160])]
    roc_axes = [fig.add_axes([0.099, 0.132, 0.345, 0.190]), fig.add_axes([0.574, 0.132, 0.345, 0.190])]
    for x in (0.099, 0.574):
        fig.text(x, 0.612, "Native BDT-score distributions", fontsize=15.0, weight="bold", color=INK)
        fig.text(x, 0.347, "ROC: held-out ranking performance", fontsize=15.0, weight="bold", color=INK)

    positive = []
    for density in (pp["signal_density"], pp["background_density"]):
        arr = np.asarray(density, dtype=float)
        positive.extend(arr[arr > 0])
    for block in auau:
        for density in (block["signal_density"], block["background_density"]):
            arr = np.asarray(density, dtype=float)
            positive.extend(arr[arr > 0])
    ymin = max(float(np.min(positive)) * 0.62, 1.2e-3)
    ymax = float(np.max(positive)) * 1.75

    for score_ax, roc_ax, block in zip(score_axes, roc_axes, auau, strict=True):
        step(score_ax, np.asarray(pp["edges"]), np.asarray(pp["signal_density"]), color=SIGNAL, linestyle="-", linewidth=2.65, alpha=0.72)
        step(score_ax, np.asarray(pp["edges"]), np.asarray(pp["background_density"]), color=BACKGROUND, linestyle="-", linewidth=2.65, alpha=0.72)
        step(score_ax, np.asarray(block["edges"]), np.asarray(block["signal_density"]), color=SIGNAL, linestyle="--", linewidth=3.05)
        step(score_ax, np.asarray(block["edges"]), np.asarray(block["background_density"]), color=BACKGROUND, linestyle="--", linewidth=3.05)
        score_ax.set_yscale("log")
        score_ax.set_xlim(0.0, 1.0)
        score_ax.set_ylim(ymin, ymax)
        score_ax.set_ylabel("Area density")
        score_ax.set_xlabel("Native BDT score")
        add_sphenix_label(score_ax)
        decorate(score_ax)

        roc_ax.plot(np.asarray(pp["fpr"]), np.asarray(pp["tpr"]), color=INK, lw=2.8, alpha=0.72, label=f"pp baseV3E  AUC {pp_auc:.3f}")
        roc_ax.plot(np.asarray(block["fpr"]), np.asarray(block["tpr"]), color=ACCENT, lw=3.05, ls="--", label=f"Au+Au {block['label']}  AUC {block['auc']:.3f}")
        roc_ax.plot([0, 1], [0, 1], color="#9ca3af", lw=1.15, ls=":")
        roc_ax.set_xlim(0.0, 1.0)
        roc_ax.set_ylim(0.0, 1.0)
        roc_ax.set_xlabel("Background efficiency")
        roc_ax.set_ylabel("Signal efficiency")
        roc_ax.legend(loc="lower right", bbox_to_anchor=(0.965, 0.38), frameon=False, fontsize=10.5, handlelength=2.4)
        add_sphenix_label(roc_ax)
        decorate(roc_ax)

    roc_axes[0].text(
        0.965,
        0.075,
        f"AUC gap to pp: {central_gap:.3f}",
        transform=roc_axes[0].transAxes,
        ha="right",
        va="bottom",
        fontsize=11.6,
        color=INK,
        bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "edgecolor": "#cbd5e1", "linewidth": 0.8, "alpha": 0.95},
    )
    roc_axes[1].text(
        0.965,
        0.075,
        f"AUC gap to pp: {peripheral_gap:.3f}\n{100.0 * recovered:.0f}% of central gap recovered",
        transform=roc_axes[1].transAxes,
        ha="right",
        va="bottom",
        fontsize=11.6,
        color=INK,
        bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "edgecolor": "#cbd5e1", "linewidth": 0.8, "alpha": 0.95},
    )

    fig.savefig(OUT_PNG, dpi=160, facecolor="white")
    plt.close(fig)


def write_csv(pp: dict[str, object], auau: list[dict[str, object]]) -> None:
    rows: list[dict[str, object]] = []
    for sample, centrality, block in [("pp", "none", pp)] + [("AuAu", str(b["label"]), b) for b in auau]:
        edges = np.asarray(block["edges"], dtype=float)
        for class_name, density in (("signal", block["signal_density"]), ("background", block["background_density"])):
            for lo, hi, value in zip(edges[:-1], edges[1:], np.asarray(density, dtype=float), strict=True):
                rows.append(
                    {
                        "record": "score_histogram",
                        "sample": sample,
                        "centrality": centrality,
                        "class": class_name,
                        "x": "",
                        "y": "",
                        "bin_low": f"{lo:.8g}",
                        "bin_high": f"{hi:.8g}",
                        "density": f"{value:.10g}",
                        "auc": f"{float(block['auc']):.10g}",
                    }
                )
        for x, y in zip(np.asarray(block["fpr"]), np.asarray(block["tpr"]), strict=True):
            rows.append(
                {
                    "record": "roc",
                    "sample": sample,
                    "centrality": centrality,
                    "class": "",
                    "x": f"{x:.10g}",
                    "y": f"{y:.10g}",
                    "bin_low": "",
                    "bin_high": "",
                    "density": "",
                    "auc": f"{float(block['auc']):.10g}",
                }
            )
    with OUT_CSV.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def write_manifest(pp: dict[str, object], auau: list[dict[str, object]]) -> None:
    central, peripheral = auau
    central_gap = float(pp["auc"]) - float(central["auc"])
    peripheral_gap = float(pp["auc"]) - float(peripheral["auc"])
    recovered = (float(peripheral["auc"]) - float(central["auc"])) / central_gap
    manifest = {
        "schema": "THE45_THE111_PP_AUAU_BDT_BRIDGE_SLIDE_V1",
        "campaign": "THE-45 JSTG July 22 practice deck",
        "reference_slide": {
            "presentation_id": "167x-He2rOOBO2i4nNS6Pdcqu7Wv03GeFMuWH9tRRx-8",
            "slide_object_id": "g3eebd15c25a_0_54",
            "reference_image_title": "the57_pp_auau_central_peripheral_bdt_bridge.png",
        },
        "outputs": {"png": str(OUT_PNG), "csv": str(OUT_CSV), "speaker_notes": str(OUT_SCRIPT)},
        "canvas_status_tag": f"sPHENIX{STATUS_TAG}",
        "pp": {
            "identity": "accepted PPG12 baseV3E 11-feature pp model-validation baseline",
            "source": str(PP_SOURCE),
            "source_sha256": sha256(PP_SOURCE),
            "auc": float(pp["auc"]),
            "selection": "15 <= reconstructed photon ET < 35 GeV, |eta| < 0.7, no NCB preselection",
            "source_choice": "The July production ROOTs are newer parity products but do not define a newer pp model or a complete exact wide-bin classifier-validation ROC. The accepted baseV3E model-validation source remains the current pp classifier baseline for this slide.",
        },
        "auau": {
            "identity": "THE-111 combined corrected-shower plus PPG12 source-role 14-feature candidate",
            "source": str(THE111_SOURCE),
            "source_sha256": sha256(THE111_SOURCE),
            "validation": str(THE111_VALIDATION),
            "validation_sha256": sha256(THE111_VALIDATION),
            "holdout_sha256": "4663dd99e9dc84812c086b175eb2bc8054194075743512fd2cd63c5e44ea45e4",
            "selection": "15 <= reconstructed photon ET < 35 GeV, |eta| < 0.7, |zvtx| < 10 cm, candidate-row holdout",
            "central_auc": float(central["auc"]),
            "peripheral_auc": float(peripheral["auc"]),
            "central_wp80_background_acceptance": float(central["wp80_background_acceptance"]),
            "peripheral_wp80_background_acceptance": float(peripheral["wp80_background_acceptance"]),
        },
        "derived": {
            "central_auc_gap_to_pp": central_gap,
            "peripheral_auc_gap_to_pp": peripheral_gap,
            "fraction_of_central_gap_recovered": recovered,
        },
        "boundaries": [
            "simulation-only comparison",
            "candidate-row holdout is not event-grouped",
            "native pp and AuAu score values are model-specific; ROC/AUC is the common ranking metric",
            "no data scoring or production rerun",
            "no Google Slides mutation",
        ],
    }
    OUT_MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n")


def write_speaker_notes(pp: dict[str, object], auau: list[dict[str, object]]) -> None:
    central, peripheral = auau
    central_gap = float(pp["auc"]) - float(central["auc"])
    peripheral_gap = float(pp["auc"]) - float(peripheral["auc"])
    recovered = (float(peripheral["auc"]) - float(central["auc"])) / central_gap
    OUT_SCRIPT.write_text(
        "\n".join(
            [
                "# THE-111 pp-to-AuAu BDT bridge - speaker notes",
                "",
                "This updates the earlier central-versus-peripheral bridge with the combined Au+Au candidate that includes both corrections: calibrated TowerInfo shower shapes on the same good 7x7 contract and PPG12 source-role signal/background labels.",
                "",
                f"The accepted pp baseV3E reference has AUC {float(pp['auc']):.3f}. The corrected Au+Au candidate reaches {float(central['auc']):.3f} in 0-20 percent and {float(peripheral['auc']):.3f} in 50-80 percent. The peripheral gap to pp is {peripheral_gap:.3f}, compared with {central_gap:.3f} in central events, so about {100.0 * recovered:.0f} percent of the central-to-pp ranking gap is recovered.",
                "",
                f"At WP80, weighted background acceptance falls from {100.0 * float(central['wp80_background_acceptance']):.1f} percent in 0-20 percent to {100.0 * float(peripheral['wp80_background_acceptance']):.1f} percent in 50-80 percent.",
                "",
                "The score curves use each model's native response and are shown only to make the signal/background separation visible. The ROC and AUC are the common cross-system comparison. This remains simulation-only candidate-row validation; the matched data and embedding production is separate.",
                "",
            ]
        )
    )


def apply_variant(status_tag: str, suffix: str, header_ribbon: bool = True) -> None:
    """Rebind the canvas status tag, header styling, and output paths."""
    global STATUS_TAG, DRAW_HEADER_RIBBON, OUT_PNG, OUT_CSV, OUT_MANIFEST, OUT_SCRIPT
    STATUS_TAG = status_tag
    DRAW_HEADER_RIBBON = header_ribbon
    if not suffix:
        return
    OUT_PNG = OUT_PNG.with_name(f"{OUT_PNG.stem}{suffix}{OUT_PNG.suffix}")
    OUT_CSV = OUT_CSV.with_name(f"{OUT_CSV.stem}{suffix}{OUT_CSV.suffix}")
    OUT_MANIFEST = OUT_MANIFEST.with_name(f"{OUT_MANIFEST.stem}{suffix}{OUT_MANIFEST.suffix}")
    OUT_SCRIPT = OUT_SCRIPT.with_name(f"{OUT_SCRIPT.stem}{suffix}{OUT_SCRIPT.suffix}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--status-tag",
        default=" Simulation",
        help="text drawn after the bold italic sPHENIX tag on each canvas (default: ' Simulation')",
    )
    parser.add_argument(
        "--variant-suffix",
        default="",
        help="suffix appended to every output filename, e.g. '_internal'",
    )
    parser.add_argument(
        "--no-header-ribbon",
        action="store_true",
        help="drop the shaded panel behind the subtitle block; the subtitle text is kept",
    )
    args = parser.parse_args()
    apply_variant(args.status_tag, args.variant_suffix, header_ribbon=not args.no_header_ribbon)

    for source in (PP_SOURCE, THE111_SOURCE, THE111_VALIDATION):
        if not source.exists():
            raise SystemExit(f"Missing required source: {source}")
    pp = load_pp()
    auau = load_auau()
    make_slide(pp, auau)
    write_csv(pp, auau)
    write_manifest(pp, auau)
    write_speaker_notes(pp, auau)
    print(OUT_PNG)
    print(OUT_CSV)
    print(OUT_MANIFEST)
    print(OUT_SCRIPT)


if __name__ == "__main__":
    main()
