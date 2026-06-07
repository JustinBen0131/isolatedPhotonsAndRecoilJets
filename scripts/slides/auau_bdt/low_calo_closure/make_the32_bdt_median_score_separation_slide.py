#!/usr/bin/env python3
"""Make a slide for BDT median signal/background score separation after the floor veto."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch


ROOT = Path("dataOutput/auauTightBDTValidation/THE32_lowCaloDiagnosticClosure_20260603")
OUTDIR = ROOT / "bdt_effect_20260606"
SUMMARY_JSON = OUTDIR / "the32_bdt_effect_summary_20260606.json"
PNG = OUTDIR / "event_energy_veto_bdt_median_score_separation_20260606.png"
PHONE_PNG = OUTDIR / "event_energy_veto_bdt_median_score_separation_phone_refresh_20260606.png"
CSV = OUTDIR / "event_energy_veto_bdt_median_score_separation_20260606.csv"
SCRIPT = OUTDIR / "event_energy_veto_bdt_median_score_separation_script_20260606.md"
MANIFEST = OUTDIR / "event_energy_veto_bdt_median_score_separation_manifest_20260606.json"

INK = "#172033"
MUTED = "#64748B"
BLUE = "#1F77B4"
LIGHT_BLUE = "#B7C9DD"
ORANGE = "#C47A1C"
GREEN = "#1F7A4D"
RED = "#D62728"
PANEL = "#F8FAFC"

CENT_KEYS = ["0_20", "20_50", "50_80"]
CENT_LABELS = ["0-20%", "20-50%", "50-80%"]

plt.rcParams.update(
    {
        "font.family": "Times New Roman",
        "mathtext.fontset": "stix",
        "axes.unicode_minus": False,
    }
)


def rounded(ax, xy, wh, face, edge="#CBD5E1", lw=1.2, radius=0.018):
    patch = FancyBboxPatch(
        xy,
        wh[0],
        wh[1],
        boxstyle=f"round,pad=0.012,rounding_size={radius}",
        facecolor=face,
        edgecolor=edge,
        linewidth=lw,
        transform=ax.transAxes,
        zorder=1,
    )
    ax.add_patch(patch)
    return patch


def load_rows() -> tuple[dict, list[dict[str, float | int | str]]]:
    data = json.loads(SUMMARY_JSON.read_text())
    rows: list[dict[str, float | int | str]] = []
    for key, label in zip(CENT_KEYS, CENT_LABELS):
        b = data["baseline"]["centrality"][key]
        c = data["upstream_cut"]["centrality"][key]
        base_gap = float(b["signal_score_median"] - b["background_score_median"])
        cut_gap = float(c["signal_score_median"] - c["background_score_median"])
        rows.append(
            {
                "centrality": label,
                "centrality_key": key,
                "baseline_signal_median": float(b["signal_score_median"]),
                "baseline_background_median": float(b["background_score_median"]),
                "baseline_median_gap": base_gap,
                "after_signal_median": float(c["signal_score_median"]),
                "after_background_median": float(c["background_score_median"]),
                "after_median_gap": cut_gap,
                "delta_median_gap": cut_gap - base_gap,
                "effect": "helps" if cut_gap > base_gap else "hurts" if cut_gap < base_gap else "unchanged",
                "baseline_entries": int(b["entries"]),
                "after_entries": int(c["entries"]),
            }
        )
    return data, rows


def write_csv(rows):
    CSV.parent.mkdir(parents=True, exist_ok=True)
    with CSV.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def metric_card(ax, x, title, value, sub, face, value_color=INK):
    rounded(ax, (x, 0.705), (0.285, 0.115), face, edge="#D0D5DD")
    ax.text(x + 0.018, 0.797, title, transform=ax.transAxes, fontsize=12.8, color=INK, fontweight="bold", va="top")
    ax.text(x + 0.018, 0.758, value, transform=ax.transAxes, fontsize=20.5, color=value_color, fontweight="bold", va="top")
    ax.text(x + 0.018, 0.718, sub, transform=ax.transAxes, fontsize=10.8, color=MUTED, va="top")


def draw_slide(data, rows):
    fig = plt.figure(figsize=(12.8, 7.2), dpi=200)
    fig.patch.set_facecolor("white")

    canvas = fig.add_axes([0, 0, 1, 1])
    canvas.axis("off")

    canvas.text(
        0.055,
        0.962,
        "Median BDT score separation changes by only about 0.002",
        fontsize=25.5,
        color=INK,
        fontweight="bold",
        va="top",
    )
    rounded(canvas, (0.055, 0.840), (0.890, 0.052), "#F1F5F9", edge="#D9E2EC")
    canvas.text(
        0.500,
        0.866,
        "Definition: median score separation = median(truth-photon BDT score) - median(inclusive-jet BDT score).",
        fontsize=13.1,
        color=INK,
        fontweight="bold",
        ha="center",
        va="center",
    )

    global_base_gap = float(data["baseline"]["signal_score_mean"] - data["baseline"]["background_score_mean"])
    global_cut_gap = float(data["upstream_cut"]["signal_score_mean"] - data["upstream_cut"]["background_score_mean"])
    deltas = [float(r["delta_median_gap"]) for r in rows]
    max_abs_delta = max(abs(d) for d in deltas)
    help_count = sum(1 for d in deltas if d > 0)

    metric_card(
        canvas,
        0.055,
        "Largest median-gap shift",
        f"{max_abs_delta:.4f}",
        "absolute score units across centrality bins",
        "#FFF8F8",
        value_color=RED,
    )
    metric_card(
        canvas,
        0.365,
        "Bins where gap increases",
        f"{help_count} / {len(rows)}",
        "larger gap means more median separation",
        "#ECFDF3",
        value_color=GREEN,
    )
    metric_card(
        canvas,
        0.675,
        "Global mean-score gap",
        f"{global_cut_gap:.3f}",
        f"baseline {global_base_gap:.3f}; change {global_cut_gap - global_base_gap:+.4f}",
        "#EAF4FB",
        value_color=BLUE,
    )

    x = np.arange(len(rows))
    width = 0.34
    base_gap = np.array([float(r["baseline_median_gap"]) for r in rows])
    cut_gap = np.array([float(r["after_median_gap"]) for r in rows])
    delta = np.array([float(r["delta_median_gap"]) for r in rows])

    ax_gap = fig.add_axes([0.075, 0.315, 0.420, 0.320])
    ax_gap.bar(x - width / 2, base_gap, width, color=LIGHT_BLUE, edgecolor="#5D6B7A", linewidth=0.7, label="Baseline")
    ax_gap.bar(x + width / 2, cut_gap, width, color=BLUE, edgecolor="#174A76", linewidth=0.7, label="After floor veto + retrain")
    ax_gap.set_title("Median signal-background score gap", fontsize=15.0, fontweight="bold", color=INK, pad=8)
    ax_gap.set_ylabel("Median score gap", fontsize=12.8)
    ax_gap.set_xticks(x)
    ax_gap.set_xticklabels([r["centrality"] for r in rows], fontsize=12.2)
    ax_gap.set_ylim(0.52, 0.81)
    ax_gap.grid(axis="y", color="#E2E8F0", linewidth=0.8)
    ax_gap.tick_params(axis="y", labelsize=11)
    ax_gap.legend(loc="upper left", frameon=False, fontsize=10.6)
    for i, d in enumerate(delta):
        effect = "helps" if d > 0 else "hurts" if d < 0 else "flat"
        color = GREEN if d > 0 else RED if d < 0 else MUTED
        ax_gap.text(
            i,
            max(base_gap[i], cut_gap[i]) + 0.010,
            f"{effect} {d:+.4f}",
            ha="center",
            va="bottom",
            fontsize=11.3,
            color=color,
            fontweight="bold",
        )

    ax_med = fig.add_axes([0.580, 0.315, 0.355, 0.320])
    sig_base = np.array([float(r["baseline_signal_median"]) for r in rows])
    sig_cut = np.array([float(r["after_signal_median"]) for r in rows])
    bkg_base = np.array([float(r["baseline_background_median"]) for r in rows])
    bkg_cut = np.array([float(r["after_background_median"]) for r in rows])
    ax_med.plot(x, sig_base, color="#8492A6", marker="o", lw=2.0, label="Signal median, baseline")
    ax_med.plot(x, sig_cut, color=BLUE, marker="o", lw=2.2, label="Signal median, after")
    ax_med.plot(x, bkg_base, color="#D6A15D", marker="s", lw=2.0, label="Background median, baseline")
    ax_med.plot(x, bkg_cut, color=ORANGE, marker="s", lw=2.2, label="Background median, after")
    ax_med.set_title("What moved: class medians", fontsize=15.0, fontweight="bold", color=INK, pad=8)
    ax_med.set_ylabel("Median BDT score", fontsize=12.8)
    ax_med.set_xticks(x)
    ax_med.set_xticklabels([r["centrality"] for r in rows], fontsize=12.2)
    ax_med.set_ylim(0.00, 0.86)
    ax_med.grid(axis="y", color="#E2E8F0", linewidth=0.8)
    ax_med.tick_params(axis="y", labelsize=11)
    ax_med.legend(loc="center right", frameon=False, fontsize=9.7)

    rounded(canvas, (0.055, 0.055), (0.890, 0.190), "#F8FAFC", edge="#D0D5DD")
    canvas.text(0.080, 0.224, "Per-centrality readout", transform=canvas.transAxes, fontsize=15.3, color=INK, fontweight="bold", va="top")
    col_x = [0.185, 0.500, 0.815]
    for xx, r in zip(col_x, rows):
        d = float(r["delta_median_gap"])
        color = GREEN if d > 0 else RED if d < 0 else MUTED
        rounded(canvas, (xx - 0.118, 0.078), (0.236, 0.105), "#FFFFFF", edge="#E2E8F0", lw=0.9, radius=0.010)
        canvas.text(xx, 0.166, str(r["centrality"]), transform=canvas.transAxes, fontsize=14.9, color=INK, fontweight="bold", ha="center", va="top")
        canvas.text(
            xx,
            0.130,
            f"{r['baseline_median_gap']:.3f} -> {r['after_median_gap']:.3f}",
            transform=canvas.transAxes,
            fontsize=12.1,
            color=INK,
            fontweight="bold",
            ha="center",
            va="top",
        )
        canvas.text(
            xx,
            0.100,
            f"{r['effect']} by {d:+.4f}",
            transform=canvas.transAxes,
            fontsize=11.3,
            color=color,
            fontweight="bold",
            ha="center",
            va="top",
        )

    for ax in [ax_gap, ax_med]:
        for spine in ax.spines.values():
            spine.set_color("#CBD5E1")
            spine.set_linewidth(1.0)

    fig.savefig(PNG, dpi=200)
    fig.savefig(PHONE_PNG, dpi=200)
    plt.close(fig)


def write_script(data, rows):
    SCRIPT.write_text(
        "\n".join(
            [
                "# BDT Median Score Separation Slide Script",
                "",
                "This slide checks whether the upstream total-calo floor veto materially changes the median score separation between truth photons and inclusive-jet background.",
                "The metric is median(truth-photon BDT score) minus median(inclusive-jet BDT score), computed in each broad centrality bin from the full-stat embedded validation summaries.",
                "The changes are tiny: about -0.002 in 0-20%, -0.002 in 20-50%, and +0.002 in 50-80%.",
                "So by median score separation, the floor veto does not strongly help or hurt the BDT. It removes the event-quality pathology while leaving the class separation essentially stable.",
                "",
            ]
        )
    )


def write_manifest(data, rows):
    MANIFEST.write_text(
        json.dumps(
            {
                "png": str(PNG),
                "phone_refresh_png": str(PHONE_PNG),
                "csv": str(CSV),
                "script": str(SCRIPT),
                "summary_json": str(SUMMARY_JSON),
                "product": data["product"],
                "metric": "median signal/background BDT score gap = signal median score - background median score",
                "comparison": "corrected diagnostic baseline vs retrained upstream-filtered total-calo floor veto model",
                "rows": rows,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )


def main() -> int:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    data, rows = load_rows()
    write_csv(rows)
    draw_slide(data, rows)
    write_script(data, rows)
    write_manifest(data, rows)
    print(PNG)
    print(PHONE_PNG)
    print(CSV)
    print(SCRIPT)
    print(MANIFEST)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
