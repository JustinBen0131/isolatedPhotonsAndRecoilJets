#!/usr/bin/env python3
"""Make simple signed-tower total-calo 5% centrality panels."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


DEFAULT_HIST_JSON = Path(
    "dataOutput/auauTightBDTValidation/THE32_lowCaloDiagnosticClosure_20260603/"
    "blair_reference_20260702/the95_signed_tower_full_count_histograms_v1.json"
)
DEFAULT_OUT = Path(
    "dataOutput/auauTightBDTValidation/THE32_lowCaloDiagnosticClosure_20260603/"
    "blair_reference_20260702/blair_5pct_log10_total_caloE_panels_signedtower_blueonly_20260702.png"
)
DEFAULT_MANIFEST = Path(
    "dataOutput/auauTightBDTValidation/THE32_lowCaloDiagnosticClosure_20260603/"
    "blair_reference_20260702/blair_5pct_log10_total_caloE_panels_signedtower_blueonly_20260702.manifest.json"
)

BLUE = "#1F77B4"
INK = "#1E293B"
MUTED = "#536173"


def step_xy(edges: np.ndarray, counts: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    return np.repeat(edges, 2)[1:-1], np.repeat(counts, 2)


def draw_hist(ax, edges: np.ndarray, counts: np.ndarray, color: str, alpha: float) -> None:
    nz = np.nonzero(counts > 0)[0]
    if nz.size == 0:
        return
    i0, i1 = int(nz[0]), int(nz[-1])
    x, y = step_xy(edges[i0 : i1 + 2], counts[i0 : i1 + 1])
    y = np.maximum(y, 0.85)
    ax.plot(x, y, color=color, lw=1.6)
    ax.fill_between(x, y, 0.85, color=color, alpha=alpha)


def source_labels(hist: dict) -> dict[str, str]:
    schema = str(hist.get("schema", ""))
    is_good = "ISGOOD" in schema.upper()
    status_on = "STATUSON" in schema.upper() or "STATUS_ON" in schema.upper()
    if is_good:
        return {
            "schema": "THE95_ISGOOD_TOWER_5PCT_LOW_CALO_PANELS_V1",
            "title_prefix": "Good-tower",
            "selection": "TowerInfo::get_isGood() required; finite negative good-tower energies included",
            "subtitle_suffix": (
                "CaloTowerStatus on; get_isGood required"
                if status_on
                else "CaloTowerStatus off; get_isGood required"
            ),
        }
    return {
        "schema": "THE95_SIGNED_TOWER_5PCT_LOW_CALO_PANELS_V1",
        "title_prefix": "Signed tower-sum",
        "selection": "finite tower energies included; no positive-energy-only cut",
        "subtitle_suffix": (
            "CaloTowerStatus on; all finite signed towers"
            if status_on
            else "CaloTowerStatus off; all finite signed towers"
        ),
    }


def render(
    hist_json: Path,
    out_png: Path,
    manifest_path: Path,
    *,
    x_axis: str = "log",
    y_scale: str = "log",
    title: str | None = None,
    subtitle: str | None = None,
    tower_selection: str | None = None,
) -> None:
    out_png.parent.mkdir(parents=True, exist_ok=True)
    hist = json.loads(hist_json.read_text())
    labels = source_labels(hist)
    log_edges = np.asarray(hist["energy_edges"], dtype=float)
    if x_axis == "gev":
        edges = np.power(10.0, log_edges) - 1.0
        xlabel = r"$E_{\rm calo}$ [GeV]"
        title_quantity = "total calo energy"
    else:
        edges = log_edges
        xlabel = r"$\log_{10}(E_{\rm calo}+1)$"
        title_quantity = "log total calo energy"

    panels = []
    ymax = 1.0
    for p in hist["panels"]:
        total = np.asarray(p["retained_hist"], dtype=float) + np.asarray(p["removed_hist"], dtype=float)
        total_n = float(total.sum())
        ymax = max(ymax, float(total.max()))
        panels.append(
            {
                "cent_lo": p["cent_lo"],
                "cent_hi": p["cent_hi"],
                "total": total,
                "total_n": total_n,
            }
        )

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "mathtext.fontset": "dejavusans",
            "axes.unicode_minus": False,
        }
    )
    fig, axes = plt.subplots(4, 4, figsize=(14.5, 9.2), dpi=180, sharex=True, sharey=True)
    fig.patch.set_facecolor("white")

    for ax, p in zip(axes.flat, panels):
        ax.set_yscale(y_scale)
        draw_hist(ax, edges, p["total"], BLUE, 0.12)
        ax.set_title(f"{int(p['cent_lo'])}-{int(p['cent_hi'])}%", fontsize=12, fontweight="bold", pad=5)
        ax.set_xlim(float(edges[0]), float(edges[-1]))
        if y_scale == "log":
            ax.set_ylim(0.85, ymax * 2.2)
        else:
            ax.set_ylim(0.0, ymax * 1.12)
        ax.grid(True, which="major", axis="y", color="#E2E8F0", lw=0.45)
        ax.tick_params(axis="both", labelsize=9, colors=INK)
        for spine in ax.spines.values():
            spine.set_color("#334155")
            spine.set_linewidth(0.8)

    for ax in axes[:, 0]:
        ax.set_ylabel("raw event count", fontsize=10.5, color=INK)
    for ax in axes[-1, :]:
        ax.set_xlabel(xlabel, fontsize=10.5, color=INK)

    fig.suptitle(
        title or f"{labels['title_prefix']} {title_quantity} by 5% centrality bin",
        fontsize=24,
        fontweight="bold",
        color=INK,
        y=0.988,
    )
    fig.text(
        0.5,
        0.935,
        subtitle
        or (
            "Embedded Photon12/20 + Jet12/20/30/40; raw event counts, no reweighting; "
            f"{labels['subtitle_suffix']}"
        ),
        ha="center",
        va="center",
        fontsize=14.5,
        color=MUTED,
    )
    fig.tight_layout(rect=[0.035, 0.045, 0.995, 0.925], h_pad=1.1, w_pad=0.7)
    fig.savefig(out_png, dpi=180, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)

    manifest = {
        "schema": labels["schema"],
        "hist_json": str(hist_json),
        "out_png": str(out_png),
        "x_axis": "E_calo [GeV]" if x_axis == "gev" else "log10(E_calo + 1)",
        "y_axis_scale": y_scale,
        "dataset": "embedded Photon12/20 + Jet12/20/30/40",
        "tower_selection": tower_selection or labels["selection"],
        "source_hist_schema": str(hist.get("schema", "")),
        "source_variant_note": str(hist.get("variant_note", "")),
        "weights": "none; raw unweighted event counts",
        "plot_filters": [
            "histogram JSON builder deduplicates per file by (run, evt)",
            "requires finite centrality and finite log10(E_calo + 1)",
            "requires 0 <= centrality < 80",
            "no BDT, photon-ID, low-E, sample-weight, reweighting, or random subsampling cut is applied by this panel renderer",
        ],
        "panels": [
            {
                "cent_lo": p["cent_lo"],
                "cent_hi": p["cent_hi"],
                "total_count": p["total_n"],
            }
            for p in panels
        ],
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--hist-json", type=Path, default=DEFAULT_HIST_JSON)
    parser.add_argument("--out-png", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--x-axis", choices=["log", "gev"], default="log")
    parser.add_argument("--y-scale", choices=["log", "linear"], default="log")
    parser.add_argument("--title", default=None)
    parser.add_argument("--subtitle", default=None)
    parser.add_argument("--tower-selection", default=None)
    args = parser.parse_args()
    render(
        args.hist_json,
        args.out_png,
        args.manifest,
        x_axis=args.x_axis,
        y_scale=args.y_scale,
        title=args.title,
        subtitle=args.subtitle,
        tower_selection=args.tower_selection,
    )


if __name__ == "__main__":
    main()
