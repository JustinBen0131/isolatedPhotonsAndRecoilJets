#!/usr/bin/env python3
"""Build the THE-42 Phase-1 PPG12 pp BDT-score reproduction packet.

This script is intentionally local-only.  It consumes existing pp BDT-score
summary JSONs/projections and records, in one manifest, which PPG12-style
panels can be made from reusable products and which need a targeted score-QA
rerun.
"""

from __future__ import annotations

import argparse
import json
import math
import shutil
import textwrap
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np


REPO = Path(__file__).resolve().parents[3]
DEFAULT_OUTDIR = (
    REPO
    / "dataOutput/ppPhotonMLPipeline/the42_ppg12_bdt_score_reproduction_20260611"
)

FIG13_SUMMARY = (
    REPO
    / "dataOutput/ppPhotonMLPipeline/ppg12_basev3E_currentIAN_truthWindowOverlayFix_20260528_1305"
    / "validation/fullsim_shuhang_overlay_truthwindow_ppg12_weighted"
    / "pp_currentian_basev3e_bdt_score_overlay_vs_ppg12_pp_noCent_bdt_22_28_noNPB_eta0-pt3-cut0_summary.json"
)
FIG13_DATA_PROJECTION = (
    REPO
    / "dataOutput/ppPhotonMLPipeline/ppg12_basev3E_currentIAN_finalShuhang_20260526_2115"
    / "validation/fullsim_shuhang_overlay"
    / "ppg12_data_h2d_bdt_eta0_pt3_cut0_projection.json"
)
FIG19_SUMMARY = (
    REPO
    / "dataOutput/ppPhotonMLPipeline/ppg12_basev3E_currentIAN_fullsim_20260521_1811"
    / "validation/fullsim_shuhang_overlay"
    / "ppg12_fig19_equiv_pp_noCent_bdt_18_22_summary.json"
)
PT1535_SUMMARY = (
    REPO
    / "dataOutput/ppPhotonMLPipeline/ppg12_basev3E_currentIAN_rawOverlayEnvFix_20260527_1630"
    / "validation/fullsim_shuhang_overlay_raw_inclusive_pt1535"
    / "pp_currentian_basev3e_bdt_score_overlay_vs_ppg12_pp_noCent_bdt_15_35_noNPB_eta0-pt3-cut0_summary.json"
)


@dataclass
class Curve:
    label: str
    values: np.ndarray
    edges: np.ndarray
    color: str
    source: str
    draw: str = "line"
    errors: np.ndarray | None = None

    @property
    def centers(self) -> np.ndarray:
        return 0.5 * (self.edges[:-1] + self.edges[1:])


def read_json_with_optional_prefix(path: Path) -> dict[str, Any]:
    text = path.read_text()
    start = text.find("{")
    if start < 0:
        raise ValueError(f"{path} does not contain a JSON object")
    return json.loads(text[start:])


def unit(values: np.ndarray) -> np.ndarray:
    total = float(np.sum(values))
    if total <= 0.0 or not math.isfinite(total):
        return values
    return values / total


def curve_from_summary(
    data: dict[str, Any],
    key_candidates: list[str],
    label: str,
    color: str,
    source: str,
    draw: str = "line",
) -> Curve:
    for key in key_candidates:
        if key in data:
            values = np.asarray(data[key], dtype=float)
            break
    else:
        raise KeyError(f"none of {key_candidates} are present in {source}")
    edges = np.asarray(data["bins"], dtype=float)
    if len(edges) != len(values) + 1:
        raise ValueError(f"{source}: bins length {len(edges)} does not match values {len(values)}")
    return Curve(label=label, values=unit(values), edges=edges, color=color, source=source, draw=draw)


def curve_from_projection(path: Path, label: str, color: str) -> Curve:
    data = read_json_with_optional_prefix(path)
    values = np.asarray(data.get("hist") or data.get("raw_hist"), dtype=float)
    edges = np.asarray(data["bins"], dtype=float)
    raw = np.asarray(data.get("raw_hist", values), dtype=float)
    errors = np.sqrt(np.clip(raw, 0.0, None))
    total = float(np.sum(raw))
    if total > 0:
        errors = errors / total
    return Curve(
        label=label,
        values=unit(values),
        edges=edges,
        color=color,
        source=str(path),
        draw="markers",
        errors=errors,
    )


def draw_sphenix_label(ax, *, y: float = 0.955) -> None:
    ax.text(
        0.045,
        y,
        r"$\bf{\it{sPHENIX}}$ Internal",
        transform=ax.transAxes,
        fontsize=13.0,
        ha="left",
        va="top",
    )


def draw_panel(
    out: Path,
    title_lines: list[str],
    curves: list[Curve],
    *,
    note: str | None = None,
    unavailable: bool = False,
) -> dict[str, Any]:
    out.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(5.4, 6.9), dpi=190)
    ax.set_facecolor("white")
    draw_sphenix_label(ax)
    ax.text(
        0.045,
        0.875,
        "\n".join(title_lines),
        transform=ax.transAxes,
        fontsize=11.0,
        ha="left",
        va="top",
        linespacing=1.18,
    )

    if unavailable:
        ax.axis("off")
        wrapped_note = "\n".join(textwrap.wrap(note or "Unavailable from current local products", width=54))
        ax.text(
            0.045,
            0.54,
            wrapped_note,
            transform=ax.transAxes,
            fontsize=10.8,
            ha="left",
            va="top",
            linespacing=1.25,
            bbox={"boxstyle": "round,pad=0.45", "facecolor": "#F8FAFC", "edgecolor": "#CBD5E1"},
        )
    else:
        ymax = 0.0
        for curve in curves:
            centers = curve.centers
            if curve.draw == "markers":
                ax.errorbar(
                    centers,
                    curve.values,
                    yerr=curve.errors,
                    fmt="o",
                    ms=4.8,
                    mfc=curve.color,
                    mec=curve.color,
                    mew=0.8,
                    lw=0.9,
                    capsize=1.8,
                    color=curve.color,
                    label=curve.label,
                    zorder=5,
                )
            else:
                ax.stairs(curve.values, curve.edges, color=curve.color, linewidth=1.8, label=curve.label)
            if curve.values.size:
                ymax = max(ymax, float(np.nanmax(curve.values)))
        ax.set_xlim(0.0, 1.0)
        ax.set_ylim(0.0, ymax * 1.27 if ymax > 0 else 1.0)
        ax.set_xlabel("BDT score", fontsize=13.0, loc="right")
        ax.set_ylabel("normalized counts", fontsize=13.0)
        ax.tick_params(direction="in", which="both", top=True, right=True, labelsize=11.5)
        ax.minorticks_on()
        ax.legend(loc="upper right", frameon=False, fontsize=10.8, handlelength=2.0)
        if note:
            wrapped_note = "\n".join(textwrap.wrap(note, width=58))
            ax.text(
                0.045,
                0.045,
                wrapped_note,
                transform=ax.transAxes,
                fontsize=8.2,
                ha="left",
                va="bottom",
                color="#334155",
                bbox={
                    "boxstyle": "round,pad=0.22",
                    "facecolor": "white",
                    "edgecolor": "none",
                    "alpha": 0.88,
                },
            )
    fig.subplots_adjust(left=0.16, right=0.96, top=0.96, bottom=0.12)
    fig.savefig(out)
    plt.close(fig)
    return {
        "output": str(out),
        "status": "blocked_missing_product" if unavailable else "generated",
        "curves": [
            {
                "label": curve.label,
                "source": curve.source,
                "integral_after_unit_norm": float(np.sum(curve.values)),
                "max": float(np.max(curve.values)) if curve.values.size else 0.0,
            }
            for curve in curves
        ],
        "note": note,
    }


def copy_source_png(summary: dict[str, Any], out: Path) -> str | None:
    plot = summary.get("plot")
    if not plot:
        return None
    source = Path(str(plot))
    if not source.is_file():
        local_guess = REPO / str(plot).lstrip("/")
        if local_guess.is_file():
            source = local_guess
    if not source.is_file():
        return None
    out.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, out)
    return str(out)


def write_report(outdir: Path, manifest: dict[str, Any]) -> Path:
    report = outdir / "decision_report.md"
    fig20_status = next((x for x in manifest["outputs"] if x["panel"] == "fig20_style_10_14_tight"), {})
    lines = [
        "# THE-42 pp BDT-score Phase-1 decision report",
        "",
        "## Decision",
        "",
        "Do not start a new broad AuAu campaign from this QA alone. Existing local pp products can be reused for Fig.13-like and Fig.19-like BDT-score sanity checks, but they do not fully satisfy the requested PPG12 table contract because the full data/NPB-tagged projections and final tight-stage score histograms are not all present locally.",
        "",
        "The next production action, if Justin approves it, should be a targeted pp score-QA pass, not a broad pp rerun: pp data, pp signal MC, and pp inclusive MC should persist true continuous BDT-score histograms for all candidates, preselection, NPB-pass, tight-ID pass, non-tight, and isolation pass/fail stages.",
        "",
        "## Evidence",
        "",
        "- Fig.13-style 22-28 GeV no-preselection plot uses the strict truth-window/PPG12-weighted existing summary plus the saved PPG12 data projection where available.",
        "- Fig.19-style 18-22 GeV preselection+NPB plot uses the existing current-IAN fullsim summary with the May-21 convention that `cluster_weta_cogx` is not cut in preselection.",
        "- Fig.20-style 10-14 GeV tight plot is blocked as a matched data/signal/inclusive reproduction from current local products.",
        "- Local RecoilJets ROOTs contain `h_tightBDTScore_allCandidates` and `h_tightBDTScore_preselected`, but not a final persisted `h_tightBDTScore_tight` distribution for the pp/AuAu score-QA contract.",
        "- The available PPG12 data projection found locally is only `h2d_bdt_eta0_pt3_cut0`, so data-marker coverage for Fig.19/Fig.20 is incomplete without restored PPG12 result ROOTs or a targeted pp score-QA pass.",
        "",
        "## Fig.20 blocker",
        "",
        str(fig20_status.get("reason", "missing final tight-stage score product")),
        "",
        "## Deliverables",
        "",
    ]
    for item in manifest["outputs"]:
        lines.append(f"- `{item['panel']}`: `{item['output']}` ({item['status']})")
    lines.extend(
        [
            f"- Manifest: `{manifest['manifest']}`",
            "",
            "## Targeted pp score-QA product required if we want exact PPG12 table reproduction",
            "",
            "Required persisted objects:",
            "",
            "- true photon-ID BDT score, not NPB score",
            "- pp data, signal MC, inclusive MC, and NPB-tagged data projections",
            "- stages: all candidates, no-preselection/cut0, preselection/cut1, NPB-pass, tight/cut2, non-tight/cut3, isolation pass/fail",
            "- pT bins: PPG12 bins including 10-14, 18-22, 22-28, plus 15-35 trained-BDT diagnostic",
            "- manifest recording weights, cuts, input branches/histograms, model registry, and source files",
            "",
        ]
    )
    report.write_text("\n".join(lines))
    return report


def build(outdir: Path) -> dict[str, Any]:
    outdir.mkdir(parents=True, exist_ok=True)
    fig13 = read_json_with_optional_prefix(FIG13_SUMMARY)
    fig19 = read_json_with_optional_prefix(FIG19_SUMMARY)
    pt1535 = read_json_with_optional_prefix(PT1535_SUMMARY)

    outputs: list[dict[str, Any]] = []
    outputs.append(
        {
            "panel": "fig13_style_22_28_no_preselection",
            **draw_panel(
                outdir / "fig13_style_22_28_no_preselection_bdt_score.png",
                [
                    r"$p$+$p$ $\sqrt{s}=200$ GeV",
                    r"$|\eta^\gamma| < 0.7$",
                    r"$22 < E_T^\gamma < 28$ GeV, w/o NPB cut",
                ],
                [
                    curve_from_projection(FIG13_DATA_PROJECTION, "Data", "black"),
                    curve_from_summary(fig13, ["this_analysis_signal_hist", "signal_raw_hist"], "Signal MC", "#E43D30", str(FIG13_SUMMARY)),
                    curve_from_summary(fig13, ["this_analysis_inclusive_hist", "inclusive_raw_hist"], "Inclusive MC", "#3B5BFF", str(FIG13_SUMMARY)),
                ],
                note="Strict existing pp summary: truth-window inclusive MC and PPG12 shower-shape weights.",
            ),
            "cuts": fig13.get("cuts"),
            "source_summary": str(FIG13_SUMMARY),
            "data_projection": str(FIG13_DATA_PROJECTION),
        }
    )
    outputs.append(
        {
            "panel": "fig19_style_18_22_preselection_npb",
            **draw_panel(
                outdir / "fig19_style_18_22_preselection_npb_bdt_score.png",
                [
                    r"$p$+$p$ $\sqrt{s}=200$ GeV",
                    r"$|\eta^\gamma| < 0.7$",
                    r"$18 < E_T^\gamma < 22$ GeV, preselection + NPB $>0.5$",
                ],
                [
                    curve_from_summary(fig19, ["signal_raw_hist"], "Signal MC", "#E43D30", str(FIG19_SUMMARY)),
                    curve_from_summary(fig19, ["inclusive_raw_hist"], "Inclusive MC", "#3B5BFF", str(FIG19_SUMMARY)),
                ],
                note="Existing summary lacks local data/NPB projection for this bin; MC-only diagnostic.",
            ),
            "cuts": fig19.get("cuts"),
            "source_summary": str(FIG19_SUMMARY),
            "data_projection": None,
        }
    )
    outputs.append(
        {
            "panel": "fig20_style_10_14_tight",
            **draw_panel(
                outdir / "fig20_style_10_14_tight_bdt_score_blocked.png",
                [
                    r"$p$+$p$ $\sqrt{s}=200$ GeV",
                    r"$|\eta^\gamma| < 0.7$",
                    r"$10 < E_T^\gamma < 14$ GeV, tight BDT selection",
                ],
                [],
                note=(
                    "Blocked from existing local products.\n\n"
                    "Need final tight-stage true BDT-score histograms or source score trees for pp data, "
                    "signal MC, and inclusive MC. Current RecoilJets histograms store all-candidate and "
                    "preselected score views, not a final tight/cut2 BDT-score view."
                ),
                unavailable=True,
            ),
            "cuts": {
                "pt_range": "10:14",
                "stage": "tight/cut2",
                "required_variable": "true photon-ID BDT score",
            },
            "reason": "missing final tight-stage pp data/signal/inclusive score product in local existing outputs",
            "source_summary": None,
            "data_projection": None,
        }
    )
    outputs.append(
        {
            "panel": "trained_bdt_15_35_no_preselection",
            **draw_panel(
                outdir / "trained_bdt_15_35_no_preselection_bdt_score.png",
                [
                    r"$p$+$p$ $\sqrt{s}=200$ GeV",
                    r"$|\eta^\gamma| < 0.7$",
                    r"$15 < E_T^\gamma < 35$ GeV, w/o NPB cut",
                ],
                [
                    curve_from_summary(pt1535, ["this_analysis_signal_hist", "signal_raw_hist"], "Signal MC", "#E43D30", str(PT1535_SUMMARY)),
                    curve_from_summary(pt1535, ["this_analysis_inclusive_hist", "inclusive_raw_hist"], "Inclusive MC", "#3B5BFF", str(PT1535_SUMMARY)),
                ],
                note="Existing raw-inclusive diagnostic; not a PPG12 figure bin.",
            ),
            "cuts": pt1535.get("cuts"),
            "source_summary": str(PT1535_SUMMARY),
            "data_projection": None,
        }
    )

    copied_sources: list[dict[str, str]] = []
    for label, summary in [("fig13_existing", fig13), ("fig19_existing", fig19), ("pt1535_existing", pt1535)]:
        copied = copy_source_png(summary, outdir / f"source_{label}.png")
        if copied:
            copied_sources.append({"label": label, "copy": copied, "source": str(summary.get("plot"))})

    manifest = {
        "schema": "THE42_PPG12_BDT_SCORE_REPRODUCTION_V1",
        "created_by": str(Path(__file__).relative_to(REPO)),
        "intent": "Phase 1 pp-only reproduction gate before any new broad AuAu submission",
        "ppg12_reference": "PPG12 current IAN figures 13, 19, 20",
        "status": "partial_reproduction_targeted_pp_score_qa_needed",
        "no_remote_or_condor_actions": True,
        "outputs": outputs,
        "copied_existing_source_pngs": copied_sources,
        "product_decision": {
            "reuse_existing_pp_for_fig13_fig19_sanity": True,
            "full_pp_table_matched": False,
            "targeted_pp_score_qa_needed": True,
            "broad_pp_rerun_needed_now": False,
            "new_auau_submission_allowed_from_this_gate": False,
        },
    }
    manifest_path = outdir / "manifest.json"
    manifest["manifest"] = str(manifest_path)
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    report = write_report(outdir, manifest)
    manifest["decision_report"] = str(report)
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    return ap.parse_args()


def main() -> int:
    args = parse_args()
    manifest = build(args.outdir)
    print(manifest["manifest"])
    print(manifest["decision_report"])
    for item in manifest["outputs"]:
        print(f"{item['panel']}: {item['output']} [{item['status']}]")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
