#!/usr/bin/env python3
"""Overlay PPG12 IAN Fig. 2 DataThief points with current RecoilJets output."""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import uproot


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
DEFAULT_DATATHIEF = (
    REPO
    / "dataOutput/ppg12PhotonYield/ppg12_fig2_datathief_overlay_20260630"
    / "ppg12_fig2_datathief_points.csv"
)
CURRENT_ARTIFACT_POINTER = (
    REPO
    / "dataOutput/current_recoiljets_artifacts/current/pp_sim_photonjet_merged/current.json"
)
LEGACY_DEFAULT_ROOT = (
    REPO
    / "dataOutput/ppg12Parity/the76_ppg12_parity_full_20260701_003024"
    / "final_roots/photonjet"
    / "jetMinPtScan_dphiScan_vz60_isoR40_isSliding_preselectionNewPPG12_tightNewPPG12_nonTightNewPPG12"
    / "photonJet5and10and20merged_SIM/RecoilJets_photonjet5plus10plus20_MERGED.root"
)


@dataclass(frozen=True)
class Series:
    key: str
    label: str
    color: str
    marker: str
    hist_label: str


SERIES = [
    Series("pt10_15", r"$10 < p_T^\gamma < 15$ GeV", "#d948a2", "o", "10_15"),
    Series("pt15_20", r"$15 < p_T^\gamma < 20$ GeV", "#2ca02c", "s", "15_20"),
    Series("pt25_30", r"$25 < p_T^\gamma < 30$ GeV", "#1f77ff", "^", "25_30"),
]
PADS = {
    "direct": {
        "hist_class": "direct",
        "title": "Direct photons",
        "ylim": (0.90, 1.05),
        "ratio_ylim": (0.985, 1.012),
    },
    "fragmentation": {
        "hist_class": "frag",
        "title": "Fragmentation photons",
        "ylim": (0.80, 1.10),
        "ratio_ylim": (0.93, 1.05),
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--datathief", type=Path, default=DEFAULT_DATATHIEF)
    parser.add_argument(
        "--root",
        type=Path,
        default=None,
        help="Current RecoilJets photon+jet ROOT. Defaults to the pp_sim_photonjet_merged current-artifact pointer.",
    )
    parser.add_argument(
        "--outdir",
        type=Path,
        default=None,
        help="Output directory. Defaults to a fig2_truth_iso_fraction_overlay directory next to the resolved campaign output.",
    )
    args = parser.parse_args()
    if args.root is None:
        args.root = resolve_current_root()
    if args.outdir is None:
        args.outdir = default_outdir_for_root(args.root)
    return args


def resolve_current_root() -> Path:
    if CURRENT_ARTIFACT_POINTER.exists():
        payload = json.loads(CURRENT_ARTIFACT_POINTER.read_text())
        roots = payload.get("root_paths") or []
        if roots:
            return Path(roots[0])
    return LEGACY_DEFAULT_ROOT


def default_outdir_for_root(root_path: Path) -> Path:
    parts = root_path.resolve().parts
    try:
        idx = parts.index("ppg12Parity")
    except ValueError:
        return root_path.parent / "fig2_truth_iso_fraction_overlay"
    if idx + 1 >= len(parts):
        return root_path.parent / "fig2_truth_iso_fraction_overlay"
    return Path(*parts[: idx + 2]) / "fig2_truth_iso_fraction_overlay"


def load_datathief(path: Path) -> dict[tuple[str, str], list[dict[str, float | str]]]:
    rows: dict[tuple[str, str], list[dict[str, float | str]]] = {}
    with path.open() as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            key = (row["source_pad"], row["series"])
            rows.setdefault(key, []).append(
                {
                    "cutoff": float(row["cutoff_gev"]),
                    "ppg12_cutoff_upper": float(int(row["point_index"]) + 1),
                    "point_index": int(row["point_index"]),
                    "fraction": float(row["fraction"]),
                    "note": row.get("note", ""),
                }
            )
    for key in rows:
        rows[key].sort(key=lambda item: float(item["cutoff"]))
    return rows


def cumulative_fraction(hist: uproot.behaviors.TH1.Histogram, cutoff: float) -> float:
    values = hist.values(flow=True)
    edges = hist.axis().edges()
    underflow = float(values[0])
    normal = np.asarray(values[1:-1], dtype=float)
    overflow = float(values[-1])
    denom = underflow + float(np.sum(normal)) + overflow
    if denom <= 0.0:
        return float("nan")
    selected = edges[1:] <= cutoff + 1.0e-9
    numer = underflow + float(np.sum(normal[selected]))
    return numer / denom


def load_current(
    root_path: Path,
    datathief: dict[tuple[str, str], list[dict[str, float | str]]],
) -> dict[tuple[str, str], list[dict[str, float | str]]]:
    current: dict[tuple[str, str], list[dict[str, float | str]]] = {}
    with uproot.open(root_path) as root_file:
        for pad in PADS:
            for series in SERIES:
                hist_name = (
                    f"SIM/h_ppg12_fig2_truthIso_{PADS[pad]['hist_class']}_pT_"
                    f"{series.hist_label}_eta07_vz30_r03"
                )
                if hist_name not in root_file:
                    raise KeyError(f"Missing current-analysis histogram: {hist_name}")
                hist = root_file[hist_name]
                rows = []
                for dt_row in datathief[(pad, series.key)]:
                    cutoff = float(dt_row["ppg12_cutoff_upper"])
                    rows.append(
                        {
                            "cutoff": float(dt_row["cutoff"]),
                            "ppg12_cutoff_upper": cutoff,
                            "fraction": cumulative_fraction(hist, cutoff),
                            "hist_name": hist_name,
                        }
                    )
                current[(pad, series.key)] = rows
    return current


def write_csv(
    path: Path,
    datathief: dict[tuple[str, str], list[dict[str, float | str]]],
    current: dict[tuple[str, str], list[dict[str, float | str]]],
) -> None:
    with path.open("w", newline="") as handle:
        fields = [
            "pad",
            "series",
            "label",
            "point_index",
            "cutoff_gev",
            "ppg12_cutoff_upper_gev",
            "ppg12_datathief_fraction",
            "current_fraction",
            "datathief_over_current",
            "current_hist",
            "datathief_note",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for pad in PADS:
            for series in SERIES:
                dt_rows = datathief[(pad, series.key)]
                cur_rows = current[(pad, series.key)]
                for idx, (dt_row, cur_row) in enumerate(zip(dt_rows, cur_rows)):
                    dt_frac = float(dt_row["fraction"])
                    cur_frac = float(cur_row["fraction"])
                    ratio = dt_frac / cur_frac if cur_frac > 0.0 else float("nan")
                    writer.writerow(
                        {
                            "pad": pad,
                            "series": series.key,
                            "label": series.label.replace("$", ""),
                            "point_index": idx,
                            "cutoff_gev": float(dt_row["cutoff"]),
                            "ppg12_cutoff_upper_gev": float(dt_row["ppg12_cutoff_upper"]),
                            "ppg12_datathief_fraction": dt_frac,
                            "current_fraction": cur_frac,
                            "datathief_over_current": ratio,
                            "current_hist": cur_row["hist_name"],
                            "datathief_note": dt_row.get("note", ""),
                        }
                    )


def draw_plot(
    path: Path,
    datathief: dict[tuple[str, str], list[dict[str, float | str]]],
    current: dict[tuple[str, str], list[dict[str, float | str]]],
) -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "dejavuserif",
            "axes.linewidth": 1.2,
        }
    )
    fig = plt.figure(figsize=(14.8, 8.6), dpi=160)
    gs = fig.add_gridspec(
        2,
        2,
        height_ratios=[3.2, 1.15],
        left=0.075,
        right=0.985,
        bottom=0.105,
        top=0.875,
        wspace=0.18,
        hspace=0.045,
    )
    top_axes = [fig.add_subplot(gs[0, idx]) for idx in range(2)]
    ratio_axes = [fig.add_subplot(gs[1, idx], sharex=top_axes[idx]) for idx in range(2)]

    for col, pad in enumerate(("direct", "fragmentation")):
        ax = top_axes[col]
        rax = ratio_axes[col]
        for series in SERIES:
            key = (pad, series.key)
            dt = datathief[key]
            cur = current[key]
            x = np.asarray([float(row["cutoff"]) for row in dt])
            y_dt = np.asarray([float(row["fraction"]) for row in dt])
            y_cur = np.asarray([float(row["fraction"]) for row in cur])
            ratio = np.divide(y_dt, y_cur, out=np.full_like(y_dt, np.nan), where=y_cur > 0)

            ax.scatter(
                x,
                y_dt,
                marker=series.marker,
                s=36,
                facecolors="none",
                edgecolors=series.color,
                linewidths=1.35,
                label=f"PPG12 {series.label}",
                zorder=4,
            )
            ax.scatter(
                x,
                y_cur,
                marker=series.marker,
                s=24,
                color=series.color,
                linewidths=0.9,
                label=f"Current {series.label}",
                zorder=3,
            )
            rax.scatter(
                x,
                ratio,
                marker=series.marker,
                s=22,
                color=series.color,
                linewidths=0.8,
                zorder=3,
            )

        ax.set_xlim(0.0, 20.0)
        ax.set_ylim(*PADS[pad]["ylim"])
        ax.set_title(PADS[pad]["title"], loc="left", fontsize=20, pad=10)
        ax.text(
            0.03,
            0.92,
            r"$R=0.3,\ |\eta^\gamma|<0.7,\ |z_{vtx}|<30$ cm",
            transform=ax.transAxes,
            fontsize=13,
        )
        ax.text(
            0.03,
            0.855,
            "open = PPG12 IAN DataThief, filled = current RecoilJets output",
            transform=ax.transAxes,
            fontsize=12.5,
        )
        ax.set_ylabel("Fraction of Events", fontsize=17)
        ax.tick_params(which="both", direction="in", top=True, right=True, labelsize=13)
        ax.minorticks_on()
        ax.grid(False)
        ax.tick_params(labelbottom=False)

        rax.axhline(1.0, color="0.45", linestyle=(0, (4, 4)), linewidth=1.1, zorder=1)
        rax.set_ylim(*PADS[pad]["ratio_ylim"])
        rax.set_xlabel(r"$E_T^{iso}$ Cutoff [GeV]", fontsize=17)
        rax.set_ylabel("IAN / current", fontsize=15)
        rax.tick_params(which="both", direction="in", top=True, right=True, labelsize=12)
        rax.minorticks_on()
        rax.grid(False)

    sample_handles = []
    source_handles = []
    for series in SERIES:
        sample_handles.append(
            plt.Line2D(
                [0],
                [0],
                marker=series.marker,
                color="none",
                markerfacecolor=series.color,
                markeredgecolor=series.color,
                markersize=7.5,
                linestyle="None",
                label=series.label,
            )
        )
    source_handles = [
        plt.Line2D(
            [0],
            [0],
            marker="o",
            color="black",
            markerfacecolor="none",
            markersize=7.5,
            linestyle="None",
            label="PPG12 IAN DataThief",
        ),
        plt.Line2D(
            [0],
            [0],
            marker="o",
            color="black",
            markerfacecolor="black",
            markersize=6.5,
            linestyle="None",
            label="Current output",
        ),
    ]
    leg1 = top_axes[1].legend(
        handles=sample_handles,
        title="pT window",
        loc="lower right",
        bbox_to_anchor=(0.58, 0.045),
        frameon=False,
        fontsize=12,
        title_fontsize=12,
    )
    top_axes[1].add_artist(leg1)
    top_axes[1].legend(
        handles=source_handles,
        title="source",
        loc="lower right",
        bbox_to_anchor=(0.98, 0.045),
        frameon=False,
        fontsize=12,
        title_fontsize=12,
    )

    fig.text(0.075, 0.972, r"PPG12 Fig. 2 truth-isolation fraction parity", fontsize=22, weight="bold", va="top")
    fig.text(0.075, 0.934, r"Photon+jet SIM: direct and fragmentation photons, three truth-$p_T$ windows", fontsize=14, va="top")
    fig.savefig(path)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    if not args.datathief.exists():
        raise FileNotFoundError(args.datathief)
    if not args.root.exists():
        raise FileNotFoundError(args.root)
    args.outdir.mkdir(parents=True, exist_ok=True)

    datathief = load_datathief(args.datathief)
    missing = [(pad, s.key) for pad in PADS for s in SERIES if (pad, s.key) not in datathief]
    if missing:
        raise KeyError(f"Missing DataThief series: {missing}")
    current = load_current(args.root, datathief)

    csv_path = args.outdir / "ppg12_fig2_datathief_vs_current_points.csv"
    png_path = args.outdir / "ppg12_fig2_datathief_vs_current_overlay_ratio.png"
    manifest_path = args.outdir / "ppg12_fig2_datathief_vs_current_manifest.json"
    write_csv(csv_path, datathief, current)
    draw_plot(png_path, datathief, current)

    ratios: dict[str, dict[str, float]] = {}
    with csv_path.open() as handle:
        for row in csv.DictReader(handle):
            key = f"{row['pad']}_{row['series']}"
            ratios.setdefault(key, {"min": float("inf"), "max": float("-inf")})
            value = float(row["datathief_over_current"])
            ratios[key]["min"] = min(ratios[key]["min"], value)
            ratios[key]["max"] = max(ratios[key]["max"], value)

    manifest = {
        "plot": str(png_path),
        "points_csv": str(csv_path),
        "ppg12_datathief_csv": str(args.datathief),
        "current_root": str(args.root),
        "current_histogram_pattern": "SIM/h_ppg12_fig2_truthIso_{direct,frag}_pT_{10_15,15_20,25_30}_eta07_vz30_r03",
        "fraction_definition": "PPG12 Fig. 2 convention: marker x is the drawn TH1 bin center, but the fraction is evaluated at the TH1 bin upper edge; current side uses point_index+1 GeV as the cutoff upper edge.",
        "ratio_definition": "PPG12 IAN DataThief fraction / current RecoilJets fraction",
        "ratio_ranges": ratios,
        "caveat": "Direct PPG12 IAN curves are visually degenerate in the DataThief extraction and share the common visible curve, as recorded in the source DataThief manifest.",
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(manifest, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
