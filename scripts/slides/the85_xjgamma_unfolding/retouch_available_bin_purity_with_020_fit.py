#!/usr/bin/env python3
"""Retouch the available-bin ABCD purity slide with the 0-20% fit overlaid."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch


REPO = Path(__file__).resolve().parents[3]
SLIDE_DIR = REPO / "dataOutput/the85_auau_xjgamma_unfolding_push/slides"
POINTS_CSV = SLIDE_DIR / "slide02_purity_leakage_corrected_available_bins_1x3_v3_points.csv"
BASE_MANIFEST = SLIDE_DIR / "slide02_purity_leakage_corrected_available_bins_1x3_v6_manifest.json"
FIT_MANIFEST = (
    REPO
    / "dataOutput/the85_auau_xjgamma_unfolding_push/unfolded_xjgamma_firstpass"
    / "slide10_auau020_purity_fit_used_for_fig1_correction_manifest.json"
)

OUT_PNG = SLIDE_DIR / "slide02_purity_leakage_corrected_available_bins_1x3_v7_fit020.png"
OUT_MANIFEST = SLIDE_DIR / "slide02_purity_leakage_corrected_available_bins_1x3_v7_fit020_manifest.json"
OUT_SCRIPT = SLIDE_DIR / "slide02_purity_leakage_corrected_available_bins_1x3_v7_fit020_speaker_script.md"

SYSTEM_ORDER = ["Au+Au 0-20%", "Au+Au 50-80%", "pp baseV3E"]


def load_rows() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with POINTS_CSV.open() as handle:
        for row in csv.DictReader(handle):
            parsed: dict[str, object] = {}
            for key, value in row.items():
                if key in {"label", "system"}:
                    parsed[key] = value
                elif value in {"", "nan", "NaN"}:
                    parsed[key] = math.nan
                elif value in {"True", "False"}:
                    parsed[key] = value == "True"
                else:
                    try:
                        parsed[key] = float(value)
                    except ValueError:
                        parsed[key] = value
            rows.append(parsed)
    return rows


def valid_rows(rows: list[dict[str, object]], label: str) -> list[dict[str, object]]:
    out = []
    for row in rows:
        if row.get("label") != label:
            continue
        needed = ["A", "B", "C", "D", "raw_purity", "raw_purity_err", "corrected_purity", "corrected_purity_err"]
        if row.get("all_abcd_found") is not True:
            continue
        if all(math.isfinite(float(row[k])) for k in needed):
            out.append(row)
    return out


def load_fit() -> dict[str, object]:
    manifest = json.loads(FIT_MANIFEST.read_text())
    fit = manifest["fit"]
    if fit.get("fit_model") != "pade11":
        raise RuntimeError(f"expected Padé[1/1] fit in {FIT_MANIFEST}")
    return fit


def eval_pade11(fit: dict[str, object], x: np.ndarray) -> np.ndarray:
    a, b, c = [float(v) for v in fit["pade11_parameters_a_b_c"]]
    return np.clip((a + b * x) / (1.0 + c * x), 0.02, 0.98)


def summary(rows: list[dict[str, object]]) -> dict[str, float]:
    return {
        "valid_bins": float(len(rows)),
        "A": float(np.nansum([float(r["A"]) for r in rows])),
        "B": float(np.nansum([float(r["B"]) for r in rows])),
        "C": float(np.nansum([float(r["C"]) for r in rows])),
        "D": float(np.nansum([float(r["D"]) for r in rows])),
        "mean_raw": float(np.nanmean([float(r["raw_purity"]) for r in rows])),
        "mean_corr": float(np.nanmean([float(r["corrected_purity"]) for r in rows])),
    }


def draw_card(fig, x: float, y: float, w: float, h: float, title: str, metric: str, line: str) -> None:
    card = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.012,rounding_size=0.010",
        transform=fig.transFigure,
        facecolor="white",
        edgecolor="#cbd5e1",
        linewidth=1.05,
        zorder=1,
    )
    fig.patches.append(card)
    fig.text(x + 0.020, y + h - 0.032, title, ha="left", va="top", fontsize=23.5, fontweight="bold", color="#111827")
    fig.text(x + 0.020, y + h - 0.076, metric, ha="left", va="top", fontsize=22.5, fontweight="bold", color="#0f4c81")
    fig.text(x + 0.020, y + h - 0.118, line, ha="left", va="top", fontsize=15.2, color="#475569")


def draw_panel(ax, rows: list[dict[str, object]], label: str, fit: dict[str, object] | None) -> None:
    x = np.array([float(r["pt_mid"]) for r in rows])
    ex = np.array([float(r["pt_width"]) / 2.0 for r in rows])
    raw = np.array([float(r["raw_purity"]) for r in rows])
    raw_e = np.array([float(r["raw_purity_err"]) for r in rows])
    corr = np.array([float(r["corrected_purity"]) for r in rows])
    corr_e = np.array([float(r["corrected_purity_err"]) for r in rows])

    ax.errorbar(
        x,
        raw,
        xerr=ex,
        yerr=raw_e,
        fmt="o",
        ms=6.0,
        mfc="black",
        mec="black",
        ecolor="black",
        elinewidth=1.1,
        capsize=0,
        zorder=4,
    )
    ax.errorbar(
        x,
        corr,
        xerr=ex,
        yerr=corr_e,
        fmt="D",
        ms=6.2,
        mfc="white",
        mec="#0072B2",
        mew=1.55,
        ecolor="#0072B2",
        elinewidth=1.15,
        capsize=0,
        zorder=5,
    )

    ax.set_title(label, fontsize=18, fontweight="bold", pad=8)
    ax.set_xlim(14.0, 35.0)
    ax.set_ylim(0.0, 1.18)
    ax.set_xticks([15, 18, 21, 24, 28, 32, 35])
    ax.set_yticks([0.0, 0.25, 0.5, 0.75, 1.0])
    ax.grid(True, color="#e2e8f0", linewidth=0.75, alpha=0.86)
    ax.tick_params(which="both", direction="in", top=True, right=True, labelsize=13.2, length=6)
    ax.tick_params(which="minor", length=3.5)
    ax.minorticks_on()
    for spine in ax.spines.values():
        spine.set_linewidth(1.1)

    if fit is not None:
        xx = np.linspace(15.0, 35.0, 260)
        ax.plot(xx, eval_pade11(fit, xx), color="#1f4cff", lw=2.15, zorder=3)
        handles = [
            Line2D([0], [0], color="#1f4cff", lw=2.15, label="Padé[1/1] fit"),
            Line2D([0], [0], marker="o", color="black", lw=0, markersize=6.5, label="raw ABCD"),
            Line2D([0], [0], marker="D", color="#0072B2", lw=0, markerfacecolor="white", markeredgewidth=1.5, markersize=6.5, label="leakage-corrected"),
        ]
        leg = ax.legend(
            handles=handles,
            loc="upper right",
            bbox_to_anchor=(0.985, 0.985),
            frameon=True,
            facecolor="white",
            edgecolor="#cbd5e1",
            framealpha=0.94,
            fontsize=10.8,
            handlelength=1.7,
            borderpad=0.45,
            labelspacing=0.35,
        )
        leg.get_frame().set_linewidth(0.8)


def main() -> None:
    rows = load_rows()
    by_label = {label: valid_rows(rows, label) for label in SYSTEM_ORDER}
    fit = load_fit()

    plt.rcParams.update(
        {
            "font.family": "Times New Roman",
            "font.serif": ["Times New Roman", "Times", "Nimbus Roman", "DejaVu Serif"],
            "mathtext.fontset": "stix",
            "axes.linewidth": 1.1,
        }
    )

    fig = plt.figure(figsize=(16, 9), dpi=160, facecolor="white")
    fig.text(0.052, 0.948, "Available-bin ABCD purity QA for 15-35 GeV xJgamma", fontsize=33.5, fontweight="bold", ha="left", va="top", color="#111827")
    card_y, card_h, gap = 0.685, 0.162, 0.018
    card_w = (0.916 - 2 * gap) / 3.0
    card_x0 = 0.052
    counts_summary: dict[str, dict[str, float]] = {}
    for i, label in enumerate(SYSTEM_ORDER):
        s = summary(by_label[label])
        counts_summary[label] = s
        available_rows = len([r for r in rows if r.get("label") == label])
        dropped = available_rows - int(s["valid_bins"])
        middle = f"{int(s['valid_bins'])} valid bins"
        if dropped:
            middle += f"; {dropped} incomplete excluded"
        draw_card(
            fig,
            card_x0 + i * (card_w + gap),
            card_y,
            card_w,
            card_h,
            label,
            f"raw {s['mean_raw']:.2f}  ->  corrected {s['mean_corr']:.2f}",
            middle,
        )

    axes = [
        fig.add_axes([0.052, 0.045, 0.293, 0.500]),
        fig.add_axes([0.365, 0.045, 0.293, 0.500]),
        fig.add_axes([0.678, 0.045, 0.293, 0.500]),
    ]
    for idx, (ax, label) in enumerate(zip(axes, SYSTEM_ORDER)):
        draw_panel(ax, by_label[label], label, fit if idx == 0 else None)
        ax.axvline(15.0, color="#94a3b8", lw=1.0, ls=":", zorder=1)
        ax.axvline(35.0, color="#94a3b8", lw=1.0, ls=":", zorder=1)
        ax.set_xlabel(r"cluster $E_T$ bin [GeV]", fontsize=14.2)
        if idx == 0:
            ax.set_ylabel("Purity estimate", fontsize=14.8)
        else:
            ax.set_yticklabels([])

    base_manifest = json.loads(BASE_MANIFEST.read_text())
    out_manifest = {
        "schema": "THE85_SLIDE02_AVAILABLE_BIN_PURITY_QA_1X3_V7_FIT020",
        "output_png": str(OUT_PNG),
        "points_csv": str(POINTS_CSV),
        "base_manifest": str(BASE_MANIFEST),
        "fit_manifest": str(FIT_MANIFEST),
        "speaker_script": str(OUT_SCRIPT),
        "definitions": base_manifest.get("definitions", {}),
        "fit_overlay": {
            "applied_to": "Au+Au 0-20%",
            "model": "Padé[1/1]",
            "parameters_a_b_c": fit.get("pade11_parameters_a_b_c"),
            "source": str(FIT_MANIFEST),
            "purpose": "In-place purity fit used to scale region C in the Au+Au 0-20 Fig.1-style correction.",
        },
        "counts_summary": counts_summary,
    }
    OUT_MANIFEST.write_text(json.dumps(out_manifest, indent=2) + "\n")
    OUT_SCRIPT.write_text(
        "This retouched version keeps the available-bin purity QA layout and overlays the Au+Au 0-20 percent Padé purity fit directly on the 0-20 panel. "
        "That makes the separate fit-only slide redundant while preserving the raw ABCD and signal-leakage-corrected points for all three systems.\n",
        encoding="utf-8",
    )
    fig.savefig(OUT_PNG, dpi=160)
    plt.close(fig)
    print(json.dumps({"slide": str(OUT_PNG), "manifest": str(OUT_MANIFEST)}, indent=2))


if __name__ == "__main__":
    main()
