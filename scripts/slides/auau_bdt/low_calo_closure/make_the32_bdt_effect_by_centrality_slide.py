#!/usr/bin/env python3
"""Make a slide showing BDT validation behavior after the total-calo floor veto."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch


ROOT = Path("dataOutput/auauTightBDTValidation/THE32_lowCaloDiagnosticClosure_20260603")
OUTDIR = ROOT / "bdt_effect_20260606"
SUMMARY_JSON = OUTDIR / "the32_bdt_effect_summary_20260606.json"
FINE5_JSON = OUTDIR / "the32_bdt_effect_fine5_wp80_summary_20260606.json"
PNG = OUTDIR / "event_energy_veto_bdt_effect_by_centrality_20260606.png"
PHONE_PNG = OUTDIR / "event_energy_veto_bdt_effect_by_centrality_phone_refresh_20260606.png"
CSV = OUTDIR / "event_energy_veto_bdt_effect_by_centrality_20260606.csv"
FINE5_CSV = OUTDIR / "event_energy_veto_bdt_effect_fine5_wp80_ratio_20260606.csv"
SCRIPT = OUTDIR / "event_energy_veto_bdt_effect_by_centrality_script_20260606.md"
MANIFEST = OUTDIR / "event_energy_veto_bdt_effect_by_centrality_manifest_20260606.json"

INK = "#172033"
MUTED = "#5D6B7A"
BLUE = "#1F77B4"
ORANGE = "#C47A1C"
RED = "#D62728"
GREEN = "#1F7A4D"
GRAY = "#E2E8F0"
PANEL = "#F8FAFC"
YELLOW = "#FFF7ED"
CARD_BLUE = "#EAF4FB"
CARD_GREEN = "#ECFDF3"
CARD_RED = "#FFF1F2"

CENT_LABELS = ["0-20%", "20-50%", "50-80%"]
CENT_KEYS = ["0_20", "20_50", "50_80"]

plt.rcParams.update(
    {
        "font.family": "Times New Roman",
        "mathtext.fontset": "stix",
        "axes.unicode_minus": False,
    }
)


def pct(x: float, digits: int = 1) -> str:
    return f"{100.0 * x:.{digits}f}%"


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


def card(ax, x, title, value, sub, *, color, face):
    rounded(ax, (x, 0.665), (0.285, 0.140), face, edge="#D0D5DD")
    ax.text(x + 0.018, 0.783, title, transform=ax.transAxes, fontsize=13.5, color=INK, fontweight="bold", va="top")
    ax.text(x + 0.018, 0.742, value, transform=ax.transAxes, fontsize=22.0, color=color, fontweight="bold", va="top")
    ax.text(x + 0.018, 0.696, sub, transform=ax.transAxes, fontsize=11.8, color=MUTED, va="top")


def load_rows() -> tuple[dict, list[dict[str, float | int | str]], list[dict[str, float | int | str]]]:
    data = json.loads(SUMMARY_JSON.read_text())
    fine = json.loads(FINE5_JSON.read_text())
    base = data["baseline"]
    cut = data["upstream_cut"]
    rows = []
    for key, label in zip(CENT_KEYS, CENT_LABELS):
        b = base["centrality"][key]
        c = cut["centrality"][key]
        removed = int(b["entries"] - c["entries"])
        rows.append(
            {
                "centrality": label,
                "centrality_key": key,
                "baseline_auc": float(b["auc"]),
                "cut_auc": float(c["auc"]),
                "delta_auc": float(c["auc"] - b["auc"]),
                "baseline_wp80_fake": float(b["wp80_background_fake_rate"]),
                "cut_wp80_fake": float(c["wp80_background_fake_rate"]),
                "delta_wp80_fake": float(c["wp80_background_fake_rate"] - b["wp80_background_fake_rate"]),
                "baseline_threshold": float(b["wp80_threshold"]),
                "cut_threshold": float(c["wp80_threshold"]),
                "baseline_entries": int(b["entries"]),
                "cut_entries": int(c["entries"]),
                "removed_entries": removed,
                "removed_fraction": float(removed / b["entries"]),
            }
        )
    fine_rows = []
    for lo in range(0, 80, 5):
        key = f"{lo}_{lo + 5}"
        b = fine["fine5"]["baseline"][key]
        c = fine["fine5"]["upstream_cut"][key]
        base_fake = float(b["wp80_background_fake_rate"])
        cut_fake = float(c["wp80_background_fake_rate"])
        fine_rows.append(
            {
                "centrality": f"{lo}-{lo + 5}%",
                "centrality_mid": lo + 2.5,
                "centrality_key": key,
                "baseline_wp80_fake": base_fake,
                "cut_wp80_fake": cut_fake,
                "fake_ratio": cut_fake / base_fake,
                "fake_ratio_percent_change": 100.0 * (cut_fake / base_fake - 1.0),
                "baseline_entries": int(b["entries"]),
                "cut_entries": int(c["entries"]),
                "removed_entries": int(b["entries"]) - int(c["entries"]),
            }
        )
    return data, rows, fine_rows


def write_csv(rows):
    CSV.parent.mkdir(parents=True, exist_ok=True)
    with CSV.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_fine5_csv(fine_rows):
    with FINE5_CSV.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fine_rows[0]))
        writer.writeheader()
        writer.writerows(fine_rows)


def draw_slide(data, rows, fine_rows):
    fig = plt.figure(figsize=(12.8, 7.2), dpi=200)
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    fig.patch.set_facecolor("white")

    canvas = fig.add_axes([0, 0, 1, 1])
    canvas.axis("off")
    canvas.text(
        0.055,
        0.962,
        "BDT validation is stable after the total-calo floor veto",
        fontsize=27,
        color=INK,
        fontweight="bold",
        va="top",
    )
    rounded(canvas, (0.055, 0.835), (0.890, 0.050), "#F1F5F9", edge="#D9E2EC")
    canvas.text(
        0.500,
        0.860,
        "Full-stat embedded validation, 15-35 GeV photons; baseline vs retrained model after removing low total-calo-energy events.",
        fontsize=13.3,
        color=INK,
        ha="center",
        va="center",
        fontweight="bold",
    )

    base = data["baseline"]
    cut = data["upstream_cut"]
    global_removed = base["counts"]["total_entries"] - cut["counts"]["total_entries"]
    global_removed_frac = global_removed / base["counts"]["total_entries"]
    auc_delta = cut["inclusive_auc"] - base["inclusive_auc"]
    fake_delta = cut["threshold_wp80"]["background_fake_rate"] - base["threshold_wp80"]["background_fake_rate"]

    card(
        canvas,
        0.055,
        "Validation rows removed by floor veto",
        f"{global_removed:,}",
        f"{pct(global_removed_frac, 2)} of corrected baseline validation rows",
        color=RED,
        face=CARD_RED,
    )
    card(
        canvas,
        0.365,
        "Inclusive AUC",
        f"{cut['inclusive_auc']:.6f}",
        f"baseline {base['inclusive_auc']:.6f}; change {auc_delta:+.6f}",
        color=GREEN,
        face=CARD_GREEN,
    )
    fake_ratio = cut["threshold_wp80"]["background_fake_rate"] / base["threshold_wp80"]["background_fake_rate"]
    card(
        canvas,
        0.675,
        "Inclusive WP80 fake-rate ratio",
        f"{fake_ratio:.3f}",
        f"after floor veto / baseline; change {100*fake_delta:+.2f} pp",
        color=BLUE,
        face=CARD_BLUE,
    )

    # Main AUC comparison
    ax_auc = fig.add_axes([0.075, 0.322, 0.415, 0.265])
    ax_auc.set_facecolor("white")
    x = np.arange(len(rows))
    width = 0.34
    base_auc = np.array([r["baseline_auc"] for r in rows], dtype=float)
    cut_auc = np.array([r["cut_auc"] for r in rows], dtype=float)
    ax_auc.bar(x - width / 2, base_auc, width=width, color="#B7C9DD", edgecolor="#5D6B7A", linewidth=0.7, label="Baseline")
    ax_auc.bar(x + width / 2, cut_auc, width=width, color=BLUE, edgecolor="#174A76", linewidth=0.7, label="After floor veto + retrain")
    ax_auc.set_title("AUC by centrality", fontsize=15, color=INK, fontweight="bold", pad=8)
    ax_auc.set_ylabel("AUC", fontsize=13)
    ax_auc.set_xticks(x)
    ax_auc.set_xticklabels([r["centrality"] for r in rows], fontsize=12)
    ax_auc.set_ylim(0.84, 0.915)
    ax_auc.grid(axis="y", color="#E2E8F0", linewidth=0.8)
    ax_auc.tick_params(axis="y", labelsize=11)
    ax_auc.legend(loc="upper left", fontsize=10.5, frameon=False)
    for i, r in enumerate(rows):
        ax_auc.text(i, max(r["baseline_auc"], r["cut_auc"]) + 0.0018, f"{r['delta_auc']:+.4f}", ha="center", va="bottom", fontsize=12.2, color=INK, fontweight="bold")

    # WP80 fake-rate ratio comparison.
    ax_fake = fig.add_axes([0.565, 0.322, 0.375, 0.265])
    fine_x = np.array([float(r["centrality_mid"]) for r in fine_rows], dtype=float)
    fake_ratio_by_cent = np.array([float(r["fake_ratio"]) for r in fine_rows], dtype=float)
    ax_fake.axhline(1.0, color="#8492A6", linewidth=1.35, linestyle="--", label="No change")
    ax_fake.plot(fine_x, fake_ratio_by_cent, marker="o", markersize=4.8, linewidth=2.0, color=ORANGE, label="After floor veto / baseline")
    ax_fake.set_title("Fake-rate ratio with centrality-dependent BDT cut", fontsize=14.2, color=INK, fontweight="bold", pad=8)
    ax_fake.set_ylabel("Fake-rate ratio\n(after floor veto / baseline)", fontsize=11.5)
    ax_fake.set_xlabel("Centrality percentile", fontsize=11.5)
    ax_fake.set_xlim(0, 80)
    ax_fake.set_xticks(np.arange(0, 81, 10))
    ax_fake.set_xticklabels([f"{int(v)}" for v in np.arange(0, 81, 10)], fontsize=10.5)
    ax_fake.set_ylim(0.990, 1.011)
    ax_fake.grid(axis="y", color="#E2E8F0", linewidth=0.8)
    ax_fake.tick_params(axis="y", labelsize=11)
    ax_fake.legend(loc="upper right", fontsize=9.8, frameon=False)
    min_idx = int(np.nanargmin(fake_ratio_by_cent))
    max_idx = int(np.nanargmax(fake_ratio_by_cent))
    for idx in sorted({min_idx, max_idx}):
        ratio = fake_ratio_by_cent[idx]
        dy = 100.0 * (ratio - 1.0)
        yoff = 0.00065 if idx == max_idx else -0.00075
        va = "bottom" if idx == max_idx else "top"
        ax_fake.text(fine_x[idx], ratio + yoff, f"{dy:+.2f}%", ha="center", va=va, fontsize=10.2, color=INK, fontweight="bold")

    rounded(canvas, (0.055, 0.045), (0.890, 0.180), "#FFF8F8", edge="#F3D4D4")
    canvas.text(0.080, 0.206, "15-35 GeV validation rows removed in each centrality bin", transform=canvas.transAxes, fontsize=15.2, color=INK, fontweight="bold", va="top")
    col_x = [0.185, 0.500, 0.815]
    for xx, row in zip(col_x, rows):
        rounded(canvas, (xx - 0.108, 0.066), (0.216, 0.094), "#FFFFFF", edge="#E2E8F0", lw=0.9, radius=0.010)
        canvas.text(xx, 0.149, row["centrality"], transform=canvas.transAxes, fontsize=14.8, color=INK, fontweight="bold", ha="center", va="top")
        canvas.text(
            xx,
            0.116,
            f"{row['removed_entries']:,} / {row['baseline_entries']:,} rows removed",
            transform=canvas.transAxes,
            fontsize=11.2,
            color=RED if row["removed_fraction"] > 0.01 else MUTED,
            fontweight="bold",
            ha="center",
            va="top",
        )
        canvas.text(
            xx,
            0.087,
            f"{pct(row['removed_fraction'], 2)} of that validation bin",
            transform=canvas.transAxes,
            fontsize=10.7,
            color=MUTED,
            ha="center",
            va="top",
        )
    for ax in [ax_auc, ax_fake]:
        for spine in ax.spines.values():
            spine.set_color("#CBD5E1")
            spine.set_linewidth(1.0)

    fig.savefig(PNG, dpi=200)
    fig.savefig(PHONE_PNG, dpi=200)
    plt.close(fig)


def write_script(data, rows):
    base = data["baseline"]
    cut = data["upstream_cut"]
    removed = base["counts"]["total_entries"] - cut["counts"]["total_entries"]
    SCRIPT.write_text(
        "\n".join(
            [
                "# Total-Calo Floor Veto BDT Effect Slide Script",
                "",
                "This slide checks whether applying the total-calo floor veto upstream changes the BDT behavior in a meaningful way.",
                f"The comparison is the corrected diagnostic baseline against the retrained upstream-filtered model, using the same BDT product without isolation inputs and the full-stat embedded validation.",
                f"Globally, the filter removes {removed:,} validation rows, which is {pct(removed / base['counts']['total_entries'], 2)} of the corrected baseline sample.",
                f"The inclusive AUC changes from {base['inclusive_auc']:.6f} to {cut['inclusive_auc']:.6f}, so the net change is only {cut['inclusive_auc'] - base['inclusive_auc']:+.6f}.",
                "The left plot shows the same stability by centrality: the AUC shifts are at the few-times-ten-to-the-minus-four level.",
                "The right plot checks a more operational quantity: the ratio of the WP80 background fake rate after the floor veto to the baseline fake rate in 5% centrality bins.",
                "That ratio stays close to one across centrality; the largest visible excursion is below one percent.",
                "So the conclusion is that the event-energy pathology is removed upstream, but the BDT separation itself is not being artificially driven by that cut.",
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
                "fine5_csv": str(FINE5_CSV),
                "script": str(SCRIPT),
                "summary_json": str(SUMMARY_JSON),
                "product": data["product"],
                "comparison": "corrected diagnostic baseline vs retrained upstream-filtered total-calo floor veto model",
                "validation_scope": "full-stat embedded validation, 15-35 GeV BDT product, centrality bins 0-20, 20-50, 50-80",
                "metrics_plotted": [
                    "AUC by centrality",
                    "ratio of background fake rate at per-centrality 80% signal-efficiency threshold after floor veto versus baseline, shown in 5% centrality bins",
                    "rows removed by centrality",
                ],
                "baseline_root": data["baseline"]["root"],
                "upstream_cut_root": data["upstream_cut"]["root"],
                "rows": rows,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )


def main() -> int:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    data, rows, fine_rows = load_rows()
    write_csv(rows)
    write_fine5_csv(fine_rows)
    draw_slide(data, rows, fine_rows)
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
