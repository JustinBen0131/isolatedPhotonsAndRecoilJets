#!/usr/bin/env python3
"""Render the THE95 signed-tower low-calo conclusion slide.

This intentionally writes the same slide PNG path that used to hold the THE58
smooth-veto slide. The point of THE95 is to replace that old cut decision with
the signed-tower diagnostic outcome.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch


DEFAULT_HIST_JSON = Path(
    "dataOutput/auauTightBDTValidation/THE32_lowCaloDiagnosticClosure_20260603/"
    "blair_reference_20260702/the95_signed_tower_full_count_histograms_v1.json"
)
DEFAULT_KBIRD = Path(
    "dataOutput/auauTightBDTValidation/THE32_lowCaloDiagnosticClosure_20260603/"
    "blair_reference_20260702/blair_cent_vs_total_caloE_GeV_ROOT_kBird_clean_noline_20260702.png"
)
DEFAULT_PNG = Path(
    "dataOutput/auauTightBDTValidation/THE58_lowCaloSmoothCut_20260614/"
    "the58_low_calo_smooth_cut_slide5_v3.png"
)
DEFAULT_MANIFEST = Path(
    "dataOutput/auauTightBDTValidation/THE58_lowCaloSmoothCut_20260614/"
    "the58_low_calo_smooth_cut_slide5_v3.manifest.json"
)
DEFAULT_OLD_VETO_MANIFEST = Path(
    "dataOutput/auauTightBDTValidation/THE58_lowCaloSmoothCut_20260614/"
    "the95_signedtower_old_the58_veto_manifest.json"
)

INK = "#1E293B"
MUTED = "#536173"
BLUE = "#1F77B4"
RED = "#C2410C"
GREEN = "#16794C"
ORANGE = "#B45309"
PALE_BLUE = "#EFF6FF"
PALE_GREEN = "#ECFDF5"
PALE_ORANGE = "#FFF7ED"
EDGE = "#CBD5E1"


def pct(x: float) -> str:
    return f"{100.0 * x:.2f}%"


def card(ax, xy, wh, face, edge=EDGE, lw=1.2, radius=0.018):
    patch = FancyBboxPatch(
        xy,
        wh[0],
        wh[1],
        boxstyle=f"round,pad=0.012,rounding_size={radius}",
        transform=ax.transAxes,
        facecolor=face,
        edgecolor=edge,
        linewidth=lw,
        zorder=1,
    )
    ax.add_patch(patch)
    return patch


def render(
    hist_json: Path,
    kbird_png: Path,
    out_png: Path,
    manifest_path: Path,
    old_veto_manifest_path: Path,
) -> None:
    out_png.parent.mkdir(parents=True, exist_ok=True)
    hist = json.loads(hist_json.read_text())
    schema = str(hist.get("schema", ""))
    is_good_tower = "ISGOOD" in schema.upper()
    source_tag = hist.get("source_tag") or ("THE95_isgoodtower_lowcalo_diag_extract_20260703" if is_good_tower else "THE95_signedtower_lowcalo_diag_extract_20260702")
    title = "Good-tower total-calo diagnostic" if is_good_tower else "Signed tower-sum low-calo diagnostic"
    subtitle = (
        "THE95 rerun sums finite calibrated energies from TowerInfo::get_isGood towers only."
        if is_good_tower
        else "THE95 rerun includes finite negative calibrated tower energies in "
        "E_calo = E_CEMC + E_IHCal + E_OHCal."
    )
    primary = (
        "Good-tower selection does not remove the low-E population."
        if is_good_tower
        else "The old broad THE58 pretraining veto is stale for signed tower sums."
    )
    implication = (
        "Do not reintroduce the old low-calo training cut from this plot. "
        "The remaining low-E population is an overlay/QA question, not a "
        "hot-tower artifact in the event-calo diagnostic."
        if is_good_tower
        else "Do not train a new default AuAu BDT with the old low-calo cut. "
        "First freeze the signed-tower diagnostic, then decide whether a "
        "separate overlay-quality gate is needed."
    )
    conclusion = (
        "The get_isGood event-calo rerun leaves the low-E rates effectively "
        "unchanged from the signed-tower diagnostic; the old THE58 broad "
        "pretraining veto remains stale."
        if is_good_tower
        else "Signed tower-sum diagnostic makes the old broad THE58 pretraining "
        "low-calo veto stale; low-E population remains a QA issue, not an "
        "automatic BDT-training cut."
    )
    low = hist["low_energy_rate_0_5"]
    stats = {
        "cache_count": hist.get("cache_count"),
        "raw_candidate_rows": hist.get("raw_candidate_rows"),
        "event_rows_after_per_cache_dedup": hist.get("event_rows_after_per_cache_dedup"),
        "low_700": low["700"],
        "low_725": low["725"],
        "low_750": low["750"],
    }

    # This comes from rerunning the old smooth-veto decision code on the signed
    # histogram. It should be tiny if the old THE58 broad veto is no longer
    # supported.
    previous_manifest = {}
    if old_veto_manifest_path.exists():
        try:
            previous_manifest = json.loads(old_veto_manifest_path.read_text())
        except json.JSONDecodeError:
            previous_manifest = {}
    recomputed_removed = previous_manifest.get("removed")
    recomputed_total = previous_manifest.get("total")
    recomputed_frac = previous_manifest.get("removed_frac")

    img = mpimg.imread(kbird_png)

    plt.rcParams.update(
        {
            "font.family": "Times New Roman",
            "mathtext.fontset": "stix",
            "axes.unicode_minus": False,
        }
    )
    fig = plt.figure(figsize=(12.8, 7.2), dpi=200)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis("off")
    fig.patch.set_facecolor("white")

    ax.text(
        0.045,
        0.94,
        title,
        fontsize=25,
        fontweight="bold",
        color=INK,
        ha="left",
        va="top",
    )
    ax.text(
        0.045,
        0.895,
        subtitle,
        fontsize=13.5,
        color=MUTED,
        ha="left",
        va="top",
    )

    # Main plot
    plot_ax = fig.add_axes([0.060, 0.300, 0.560, 0.500])
    plot_ax.imshow(img)
    plot_ax.axis("off")
    ax.text(
        0.060,
        0.813,
        "Blair reference view: centrality vs total calo energy (GeV)",
        fontsize=13.5,
        fontweight="bold",
        color=INK,
        ha="left",
        va="bottom",
    )

    # Result cards
    card(ax, (0.655, 0.665), (0.300, 0.165), PALE_GREEN, edge="#86EFAC")
    ax.text(0.675, 0.792, "Primary outcome", fontsize=14.5, fontweight="bold", color=GREEN, va="top")
    ax.text(
        0.675,
        0.758,
        primary,
        fontsize=12.5,
        color=INK,
        va="top",
        wrap=True,
    )
    if recomputed_frac is not None:
        ax.text(
            0.675,
            0.705,
            f"Old smooth-veto logic removes {int(round(recomputed_removed)):,} / "
            f"{int(round(recomputed_total)):,} = {pct(float(recomputed_frac))}.",
            fontsize=12.0,
            color=INK,
            va="top",
            wrap=True,
        )

    card(ax, (0.655, 0.450), (0.300, 0.185), PALE_BLUE, edge="#93C5FD")
    ax.text(0.675, 0.593, "0-5% low-E check", fontsize=14.5, fontweight="bold", color=BLUE, va="top")
    lines = [
        f"E_calo < 700 GeV: {pct(stats['low_700']['fraction'])}",
        f"E_calo < 725 GeV: {pct(stats['low_725']['fraction'])}",
        f"E_calo < 750 GeV: {pct(stats['low_750']['fraction'])}",
    ]
    ax.text(0.675, 0.552, "\n".join(lines), fontsize=12.0, color=INK, va="top", linespacing=1.20)
    card(ax, (0.655, 0.235), (0.300, 0.185), PALE_ORANGE, edge="#FDBA74")
    ax.text(0.675, 0.380, "Production implication", fontsize=14.5, fontweight="bold", color=ORANGE, va="top")
    ax.text(
        0.675,
        0.342,
        implication,
        fontsize=12.2,
        color=INK,
        va="top",
        wrap=True,
    )

    # Bottom provenance band
    card(ax, (0.045, 0.078), (0.910, 0.145), "white", edge=EDGE)
    ax.text(0.065, 0.190, "Provenance", fontsize=13.5, fontweight="bold", color=INK, va="top")
    prov = (
        f"8574 ROOTs from {source_tag}; "
        f"{stats['raw_candidate_rows']:,} candidate rows; "
        f"{stats['event_rows_after_per_cache_dedup']:,} event rows after per-file dedup.\n"
        + (
            "Tower-sum code requires tower->get_isGood() and skips only non-finite energies."
            if is_good_tower
            else "Tower-sum code skips only non-finite tower energies, not E <= 0."
        )
    )
    ax.text(0.065, 0.157, prov, fontsize=10.6, color=MUTED, va="top", linespacing=1.15)
    ax.text(
        0.065,
        0.107,
        "Interpretation: this fixes the diagnostic definition, but it does not by itself prove the "
        "overlay low-energy population is gone. Decide from tower-sum QA before any new AuAu BDT training.",
        fontsize=10.8,
        color=RED,
        va="top",
        wrap=True,
    )

    fig.savefig(out_png, dpi=200, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)

    manifest = {
        "schema": "THE95_ISGOOD_TOWER_LOW_CALO_CONCLUSION_SLIDE_V1"
        if is_good_tower
        else "THE95_SIGNED_TOWER_LOW_CALO_CONCLUSION_SLIDE_V1",
        "hist_json": str(hist_json),
        "kbird_png": str(kbird_png),
        "out_png": str(out_png),
        "old_smooth_veto_recomputed": {
            "manifest": str(old_veto_manifest_path),
            "removed": recomputed_removed,
            "total": recomputed_total,
            "removed_frac": recomputed_frac,
        },
        "low_energy_rate_0_5": low,
        "cache_count": stats["cache_count"],
        "raw_candidate_rows": stats["raw_candidate_rows"],
        "event_rows_after_per_cache_dedup": stats["event_rows_after_per_cache_dedup"],
        "conclusion": conclusion,
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--hist-json", type=Path, default=DEFAULT_HIST_JSON)
    parser.add_argument("--kbird-png", type=Path, default=DEFAULT_KBIRD)
    parser.add_argument("--out-png", type=Path, default=DEFAULT_PNG)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--old-veto-manifest", type=Path, default=DEFAULT_OLD_VETO_MANIFEST)
    args = parser.parse_args()
    render(args.hist_json, args.kbird_png, args.out_png, args.manifest, args.old_veto_manifest)


if __name__ == "__main__":
    main()
