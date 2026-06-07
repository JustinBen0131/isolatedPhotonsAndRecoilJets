#!/usr/bin/env python3
"""Build local THE-32 upstream low-calo closure artifacts from confirmed evidence."""

from __future__ import annotations

import csv
import json
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


OUTDIR = Path(
    "dataOutput/auauTightBDTValidation/"
    "THE32_lowCaloDiagnosticClosure_20260603/upstream_closure_20260605"
)

BASELINE_REPORT_ROOT = (
    "/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/"
    "auauTightBDT_THE32_lowcalo_diag_extract_fixcalo_20260604/reports/"
    "model_validation_condor_THE32_lowcalo_diag_scorecache_fixcalo_fullstat_20260604"
)
UPSTREAM_REPORT_ROOT = (
    "/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/"
    "auauTightBDT_THE32_lowcalo_diag_extract_fixcalo_20260604/reports/"
    "model_validation_condor_THE32_lowcalo_upstreamcut_weightfix_noiso_fullstat_20260605"
)
TRAINING_TAG = (
    "THE32_lowcalo_upstreamcut_jet12_20_30_40_global_noiso_bdt_"
    "frommatrix_fastaudit_j40fix_weightfix_noiso_20260605"
)
PRODUCT = "globalEtCent1535_bdt_noIso"


EVIDENCE = {
    "baseline": {
        "tag": "THE32_lowcalo_diag_scorecache_fixcalo_fullstat_20260604",
        "status": "READY",
        "report_root": BASELINE_REPORT_ROOT,
        "model_dir": (
            "/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/bdt_models/"
            "THE8_branchA_jet12_20_30_40_global_noiso_bdt_mem24_20260527"
        ),
        "total_entries": 20369392,
        "signal_entries": 8807166,
        "background_entries": 11562226,
        "scored_entries": 20369392,
        "auc": 0.894693,
        "finite_score_fraction": 1.0,
        "notes": "none",
    },
    "upstream": {
        "tag": "THE32_lowcalo_upstreamcut_weightfix_noiso_fullstat_20260605",
        "status": "READY",
        "report_root": UPSTREAM_REPORT_ROOT,
        "model_dir": (
            "/gpfs/mnt/gpfs02/sphenix/user/patsfan753/thesisAnalysis/bdt_models/"
            f"{TRAINING_TAG}"
        ),
        "total_entries": 19788523,
        "signal_entries": 8550678,
        "background_entries": 11237845,
        "scored_entries": 19788523,
        "auc": 0.895321,
        "finite_score_fraction": 1.0,
        "notes": "none",
        "retained_below_envelope_events": 0,
        "retained_below_envelope_candidates": 0,
        "training_auc": 0.8892111484693804,
        "non_placeholder_ppg12_exact_weights": True,
    },
}


def pct(value: float) -> str:
    return f"{100.0 * value:.2f}%"


def write_audit_csv(path: Path) -> list[dict[str, object]]:
    baseline = EVIDENCE["baseline"]
    upstream = EVIDENCE["upstream"]
    removed_total = baseline["total_entries"] - upstream["total_entries"]
    removed_signal = baseline["signal_entries"] - upstream["signal_entries"]
    removed_background = baseline["background_entries"] - upstream["background_entries"]
    rows = [
        {
            "check": "validation_status",
            "baseline": baseline["status"],
            "upstream_cut": upstream["status"],
            "closure_interpretation": "both validations READY",
        },
        {
            "check": "full_stat_entries",
            "baseline": baseline["total_entries"],
            "upstream_cut": upstream["total_entries"],
            "closure_interpretation": f"{removed_total} candidate rows removed ({pct(removed_total / baseline['total_entries'])})",
        },
        {
            "check": "signal_entries",
            "baseline": baseline["signal_entries"],
            "upstream_cut": upstream["signal_entries"],
            "closure_interpretation": f"{removed_signal} signal rows removed ({pct(removed_signal / baseline['signal_entries'])})",
        },
        {
            "check": "background_entries",
            "baseline": baseline["background_entries"],
            "upstream_cut": upstream["background_entries"],
            "closure_interpretation": f"{removed_background} background rows removed ({pct(removed_background / baseline['background_entries'])})",
        },
        {
            "check": "validation_auc",
            "baseline": baseline["auc"],
            "upstream_cut": upstream["auc"],
            "closure_interpretation": f"AUC change {upstream['auc'] - baseline['auc']:+.6f}",
        },
        {
            "check": "finite_score_fraction",
            "baseline": baseline["finite_score_fraction"],
            "upstream_cut": upstream["finite_score_fraction"],
            "closure_interpretation": "all scored entries finite",
        },
        {
            "check": "retained_below_threshold_events",
            "baseline": "not applied upstream",
            "upstream_cut": upstream["retained_below_envelope_events"],
            "closure_interpretation": "zero retained events below the event-energy threshold",
        },
        {
            "check": "retained_below_threshold_candidates",
            "baseline": "not applied upstream",
            "upstream_cut": upstream["retained_below_envelope_candidates"],
            "closure_interpretation": "zero retained candidates below the event-energy threshold",
        },
        {
            "check": "training_weight_guard",
            "baseline": "legacy model",
            "upstream_cut": "non-placeholder PPG12 exact weights",
            "closure_interpretation": "training handoff rejected all-unit placeholder weights",
        },
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return rows


def draw_summary(path: Path) -> None:
    baseline = EVIDENCE["baseline"]
    upstream = EVIDENCE["upstream"]
    removed_total = baseline["total_entries"] - upstream["total_entries"]
    removal_fraction = removed_total / baseline["total_entries"]
    auc_delta = upstream["auc"] - baseline["auc"]

    plt.rcParams.update(
        {
            "font.family": "Times New Roman",
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "axes.edgecolor": "#344054",
            "axes.labelcolor": "#101828",
            "xtick.color": "#344054",
            "ytick.color": "#344054",
            "text.color": "#101828",
        }
    )
    fig = plt.figure(figsize=(12.8, 7.2), dpi=200)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis("off")

    ax.text(0.055, 0.925, "THE-32 upstream event-energy filter closure", fontsize=26, fontweight="bold", va="top")
    ax.text(
        0.055,
        0.875,
        "noIso BDT retrained and full-stat validated after applying the total-calo-energy event filter",
        fontsize=14.5,
        color="#475467",
        va="top",
    )

    cards = [
        {
            "x": 0.055,
            "title": "1. Filter applied before training",
            "accent": "#C4320A",
            "face": "#FFF4ED",
            "big": "0 retained",
            "big_sub": "events/candidates below threshold",
            "line1": f"{removed_total:,} rows removed vs baseline",
            "line2": f"{pct(removal_fraction)} of corrected diagnostic validation",
        },
        {
            "x": 0.365,
            "title": "2. Validation stayed stable",
            "accent": "#026AA2",
            "face": "#F0F9FF",
            "big": f"{upstream['auc']:.6f}",
            "big_sub": f"AUC after cut ({auc_delta:+.6f})",
            "line1": f"baseline AUC {baseline['auc']:.6f}",
            "line2": "finite score fraction = 1",
        },
        {
            "x": 0.675,
            "title": "3. Training handoff guarded",
            "accent": "#027A48",
            "face": "#F0FDF4",
            "big": f"{upstream['training_auc']:.6f}",
            "big_sub": "training AUC",
            "line1": "PPG12 exact weights are non-placeholder",
            "line2": "model_registry.json status = READY",
        },
    ]
    for card in cards:
        x = card["x"]
        ax.add_patch(Rectangle((x, 0.395), 0.27, 0.315, facecolor=card["face"], edgecolor="#D0D5DD", linewidth=1.1))
        ax.add_patch(Rectangle((x, 0.690), 0.27, 0.020, facecolor=card["accent"], edgecolor=card["accent"], linewidth=0))
        ax.text(x + 0.016, 0.655, card["title"], fontsize=13.7, fontweight="bold", va="top")
        ax.text(x + 0.016, 0.596, card["big"], fontsize=20.0, fontweight="bold", color=card["accent"], va="top")
        ax.text(x + 0.016, 0.548, card["big_sub"], fontsize=11.5, color="#344054", va="top")
        ax.text(x + 0.016, 0.492, card["line1"], fontsize=11.2, color="#101828", va="top")
        ax.text(x + 0.016, 0.455, card["line2"], fontsize=11.2, color="#101828", va="top")

    ax.add_patch(Rectangle((0.055, 0.185), 0.890, 0.135, facecolor="#F8FAFC", edgecolor="#D0D5DD", linewidth=1.1))
    ax.text(0.075, 0.292, "Evidence used", fontsize=13.5, fontweight="bold", va="top")
    ax.text(0.075, 0.255, "Corrected baseline", fontsize=11.0, fontweight="bold", va="top")
    ax.text(0.245, 0.255, "READY, 20,369,392 scored rows, AUC 0.894693", fontsize=11.0, color="#344054", va="top")
    ax.text(0.075, 0.216, "Upstream-cut validation", fontsize=11.0, fontweight="bold", va="top")
    ax.text(0.245, 0.216, "READY, 19,788,523 scored rows, AUC 0.895321, notes=none", fontsize=11.0, color="#344054", va="top")
    ax.text(0.705, 0.255, "Product", fontsize=11.0, fontweight="bold", va="top")
    ax.text(0.785, 0.255, PRODUCT, fontsize=11.0, color="#344054", va="top")

    ax.text(0.055, 0.105, "Conclusion", fontsize=14.5, fontweight="bold", va="top")
    ax.text(
        0.155,
        0.105,
        "The low-event-energy pathology is removed upstream and the full-stat noIso BDT validation remains stable.",
        fontsize=13.5,
        color="#101828",
        va="top",
    )
    fig.savefig(path, bbox_inches="tight", pad_inches=0.0)
    plt.close(fig)


def write_report(path: Path, audit_rows: list[dict[str, object]], manifest_path: Path, slide_path: Path) -> None:
    baseline = EVIDENCE["baseline"]
    upstream = EVIDENCE["upstream"]
    removed_total = baseline["total_entries"] - upstream["total_entries"]
    text = f"""# THE-32 Upstream Low-Calo Closure Report

## Verdict

**CLOSED for the noIso validation path.**

The event-energy filter has now been applied upstream before BDT training. The
resulting full-stat validation is READY, retains zero events/candidates below
the stored threshold, and has stable validation performance relative to the
corrected diagnostic baseline.

## Evidence

- Training tag: `{TRAINING_TAG}`
- Product: `{PRODUCT}`
- Upstream training READY: one model, training AUC `{upstream['training_auc']}`,
  non-placeholder PPG12 exact weights from `__ppg12_exact_training_weight`.
- Upstream validation READY: `{upstream['total_entries']}` total/scored rows,
  AUC `{upstream['auc']}`, finite score fraction `{upstream['finite_score_fraction']}`,
  notes `{upstream['notes']}`.
- Upstream event-quality audit: retained below threshold events
  `{upstream['retained_below_envelope_events']}`, retained below threshold
  candidates `{upstream['retained_below_envelope_candidates']}`.
- Baseline validation READY: `{baseline['total_entries']}` total/scored rows,
  AUC `{baseline['auc']}`, finite score fraction `{baseline['finite_score_fraction']}`,
  notes `{baseline['notes']}`.
- Rows removed by upstream filter relative to corrected baseline:
  `{removed_total}` candidate rows (`{pct(removed_total / baseline['total_entries'])}`).

## Remote Evidence Roots

- Baseline report root: `{BASELINE_REPORT_ROOT}`
- Upstream-cut report root: `{UPSTREAM_REPORT_ROOT}`

## Local Package

- Audit CSV: `{OUTDIR / 'the32_upstream_closure_audit_20260605.csv'}`
- Manifest: `{manifest_path}`
- Summary slide PNG: `{slide_path}`

## Audit Table

| check | baseline | upstream cut | interpretation |
| --- | --- | --- | --- |
"""
    for row in audit_rows:
        text += (
            f"| {row['check']} | {row['baseline']} | {row['upstream_cut']} | "
            f"{row['closure_interpretation']} |\n"
        )
    path.write_text(text)


def main() -> int:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    audit_path = OUTDIR / "the32_upstream_closure_audit_20260605.csv"
    manifest_path = OUTDIR / "the32_upstream_closure_manifest_20260605.json"
    report_path = OUTDIR / "the32_upstream_closure_report_20260605.md"
    slide_path = OUTDIR / "the32_upstream_closure_summary_slide_20260605.png"
    audit_rows = write_audit_csv(audit_path)
    draw_summary(slide_path)
    manifest = {
        "created_at_local": datetime.now().astimezone().isoformat(timespec="seconds"),
        "campaign": "THE-32 low-calo upstream closure",
        "verdict": "CLOSED_FOR_NOISO_VALIDATION_PATH",
        "product": PRODUCT,
        "training_tag": TRAINING_TAG,
        "evidence": EVIDENCE,
        "artifacts": {
            "audit_csv": str(audit_path),
            "manifest_json": str(manifest_path),
            "report_md": str(report_path),
            "summary_slide_png": str(slide_path),
            "blair_marker_png": (
                "dataOutput/auauTightBDTValidation/THE32_lowCaloDiagnosticClosure_20260603/"
                "the32_event_total_calo_energy_threshold_readable_20260605_223429.png"
            ),
        },
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    write_report(report_path, audit_rows, manifest_path, slide_path)
    print(f"wrote {audit_path}")
    print(f"wrote {manifest_path}")
    print(f"wrote {report_path}")
    print(f"wrote {slide_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
