#!/usr/bin/env python3
"""Make a PPG12 Fig. 3-style truth-pT overlay for photon+jet samples.

The PPG12 reference macro is:
  ppg12codeGit/efficiencytool/FindTruthETCut.C

It reads h_direct_pT_truth_isoET_0 and h_frag_pT_truth_isoET_0 from
MC_efficiency_noiso.root, projects E_T^iso <= 4 GeV, rebins by 10, and draws
total/direct/fragmentation spectra with draw_1D_multiple_plot_ratio().

This helper keeps the PPG12 canvas/pad/margin contract, but replaces the lower
pad with PPG12 SDCC / current-output ratios for all three spectra.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import uproot


REPO = Path(__file__).resolve().parents[3]
CURRENT_POINTER = (
    REPO
    / "dataOutput/current_recoiljets_artifacts/current/pp_sim_photonjet_merged/current.json"
)
DEFAULT_PPG12_ROOT = (
    REPO
    / "dataOutput/ppg12Parity/reference_roots/ppg12/efficiencytool/MC_efficiency_noiso.root"
)
DEFAULT_PPG12_POINTS_CSV = (
    REPO
    / "dataOutput/ppg12Parity/the76_ppg12_fig3_sdcc_reference/ppg12_fig3_sdcc_root_points.csv"
)
DEFAULT_OUT_DIR = (
    REPO
    / "dataOutput/ppg12Parity/the76_ppg12_fig24_photonjet_fix_20260702_014217"
    / "fig3_truth_pt_overlay"
)

PPG12_DIRECT_H2 = "h_direct_pT_truth_isoET_0"
PPG12_FRAG_H2 = "h_frag_pT_truth_isoET_0"
CURRENT_HISTS = {
    "total": "SIM/h_ppg12_fig3_truthPt_total_eta07_vz30_r03_iso4",
    "direct": "SIM/h_ppg12_fig3_truthPt_direct_eta07_vz30_r03_iso4",
    "frag": "SIM/h_ppg12_fig3_truthPt_frag_eta07_vz30_r03_iso4",
}


@dataclass(frozen=True)
class Series:
    key: str
    label: str
    color: str
    marker_open: str
    marker_closed: str
    zorder: int


SERIES = [
    Series("total", "Total", "#e15b9a", "o", "o", 8),
    Series("direct", "Direct", "#2c7a19", "s", "s", 7),
    Series("frag", "Fragmentation", "#2b83ff", "^", "^", 6),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ppg12-root", type=Path, default=DEFAULT_PPG12_ROOT)
    parser.add_argument(
        "--ppg12-points-csv",
        type=Path,
        help="CSV extracted from the PPG12 SDCC ROOT source with series/xcenter_gev/count columns.",
    )
    parser.add_argument("--ppg12-datathief-csv", type=Path)
    parser.add_argument("--current-root", type=Path)
    parser.add_argument("--current-pointer", type=Path, default=CURRENT_POINTER)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--current-only-audit", action="store_true")
    parser.add_argument("--ratio-ymin", type=float, default=0.80)
    parser.add_argument("--ratio-ymax", type=float, default=1.20)
    return parser.parse_args()


def current_root_from_pointer(path: Path) -> Path:
    meta = json.loads(path.read_text())
    roots = meta.get("root_paths") or []
    if len(roots) != 1:
        raise RuntimeError(f"{path} must contain exactly one root_paths entry")
    return Path(roots[0])


def root_like_find_bin(edges: np.ndarray, value: float) -> int:
    """ROOT-like normal-bin index, 1-based and clamped to normal bins."""
    nbins = len(edges) - 1
    if value < edges[0]:
        return 1
    if value >= edges[-1]:
        return nbins
    return int(np.searchsorted(edges, value, side="right"))


def project_iso_leq4(hist: uproot.behaviors.TH2.Histogram) -> tuple[np.ndarray, np.ndarray]:
    values, xedges, yedges = hist.to_numpy(flow=False)
    values = values.astype(float)
    ymask = yedges[1:] <= 4.0 + 1.0e-12
    if not np.any(ymask):
        raise RuntimeError(f"{hist.name} has no normal y bins with upper edge <= 4 GeV")
    projected = np.sum(values[:, ymask], axis=1)
    return projected, xedges.astype(float)


def rebin_like_ppg12(values: np.ndarray, edges: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Apply PPG12's Rebin(10) when compatible; otherwise preserve 1 GeV bins."""
    widths = np.diff(edges)
    if len(values) % 10 == 0 and np.allclose(widths, widths[0]) and math.isclose(widths[0] * 10.0, 1.0, rel_tol=0, abs_tol=1e-6):
        rebinned = values.reshape(-1, 10).sum(axis=1)
        new_edges = edges[::10]
        if len(new_edges) != len(rebinned) + 1:
            new_edges = np.append(new_edges, edges[-1])
        return rebinned, new_edges
    if np.allclose(widths, 1.0, rtol=0, atol=1e-6):
        return values, edges
    raise RuntimeError(
        "Cannot reproduce PPG12 Fig.3 binning: expected 0.1 GeV bins for Rebin(10) "
        f"or already-1 GeV bins, found nbins={len(values)} first_width={widths[0]:.6g}"
    )


def load_ppg12_reference(path: Path) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    if not path.exists():
        raise FileNotFoundError(
            f"Missing PPG12 Fig.3 source ROOT: {path}\n"
            "Expected local copy of SDCC /sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_noiso.root"
        )
    with uproot.open(path) as handle:
        missing = [name for name in (PPG12_DIRECT_H2, PPG12_FRAG_H2) if name not in handle]
        if missing:
            raise KeyError(f"{path} is missing PPG12 Fig.3 objects: {missing}")
        direct, edges = rebin_like_ppg12(*project_iso_leq4(handle[PPG12_DIRECT_H2]))
        frag, frag_edges = rebin_like_ppg12(*project_iso_leq4(handle[PPG12_FRAG_H2]))
    if not np.allclose(edges, frag_edges):
        raise RuntimeError("PPG12 direct and fragmentation pT axes differ")
    return {
        "total": (direct + frag, edges),
        "direct": (direct, edges),
        "frag": (frag, edges),
    }


def load_ppg12_datathief(path: Path) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    if not path.exists():
        raise FileNotFoundError(f"Missing PPG12 Fig.3 DataThief CSV: {path}")
    by_series: dict[str, list[tuple[float, float]]] = {series.key: [] for series in SERIES}
    with path.open() as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            key = row["series"]
            if key not in by_series:
                continue
            by_series[key].append((float(row["expected_xcenter_gev"]), float(row["counts"])))
    out: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for series in SERIES:
        rows = sorted(by_series[series.key])
        if not rows:
            raise RuntimeError(f"{path} has no DataThief rows for {series.key}")
        centers_arr = np.asarray([row[0] for row in rows], dtype=float)
        values = np.asarray([row[1] for row in rows], dtype=float)
        if not np.allclose(np.diff(centers_arr), 1.0, atol=1e-6):
            raise RuntimeError(f"{series.key} DataThief centers are not contiguous 1 GeV bins: {centers_arr}")
        edges = np.concatenate([[centers_arr[0] - 0.5], centers_arr + 0.5])
        out[series.key] = (values, edges)
    return out


def load_ppg12_points_csv(path: Path) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    if not path.exists():
        raise FileNotFoundError(f"Missing PPG12 Fig.3 SDCC points CSV: {path}")
    by_series: dict[str, list[tuple[float, float]]] = {series.key: [] for series in SERIES}
    with path.open() as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            key = row.get("series", "")
            if key not in by_series:
                continue
            x_raw = row.get("xcenter_gev") or row.get("xcenter") or row.get("expected_xcenter_gev")
            y_raw = row.get("count") or row.get("counts")
            if x_raw is None or y_raw is None:
                raise RuntimeError(f"{path} row is missing x/count columns: {row}")
            by_series[key].append((float(x_raw), float(y_raw)))
    out: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for series in SERIES:
        rows = sorted(by_series[series.key])
        if not rows:
            raise RuntimeError(f"{path} has no rows for {series.key}")
        centers_arr = np.asarray([row[0] for row in rows], dtype=float)
        values = np.asarray([row[1] for row in rows], dtype=float)
        if len(centers_arr) < 2:
            raise RuntimeError(f"{path} has too few rows for {series.key}")
        widths = np.diff(centers_arr)
        if not np.allclose(widths, widths[0], rtol=0, atol=1e-6):
            raise RuntimeError(f"{series.key} SDCC ROOT centers are not uniform: {centers_arr}")
        half_width = 0.5 * float(widths[0])
        edges = np.concatenate([[centers_arr[0] - half_width], centers_arr + half_width])
        out[series.key] = (values, edges)
    return out


def load_current(path: Path) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    if not path.exists():
        raise FileNotFoundError(f"Missing current photon+jet ROOT: {path}")
    out: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    with uproot.open(path) as handle:
        for key, hist_name in CURRENT_HISTS.items():
            if hist_name not in handle:
                raise KeyError(f"{path} is missing current Fig.3 histogram {hist_name}")
            hist = handle[hist_name]
            values, edges = hist.to_numpy(flow=False)
            out[key] = (values.astype(float), edges.astype(float))
    return out


def centers(edges: np.ndarray) -> np.ndarray:
    return 0.5 * (edges[:-1] + edges[1:])


def select_range(values: np.ndarray, edges: np.ndarray, xmin: float = 10.0, xmax: float = 35.0) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x = centers(edges)
    mask = (x >= xmin) & (x <= xmax)
    return x[mask], values[mask], np.diff(edges)[mask]


def series_rows(source: str, spectra: dict[str, tuple[np.ndarray, np.ndarray]]) -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = []
    for series in SERIES:
        values, edges = spectra[series.key]
        x, y, width = select_range(values, edges)
        for bin_index, (xc, yc, bw) in enumerate(zip(x, y, width), start=1):
            rows.append(
                {
                    "source": source,
                    "series": series.key,
                    "label": series.label,
                    "bin_index_in_range": bin_index,
                    "xcenter": float(xc),
                    "bin_width": float(bw),
                    "count": float(yc),
                    "stat_err": math.sqrt(yc) if yc >= 0 else math.nan,
                }
            )
    return rows


def ratio_rows(ppg12: dict[str, tuple[np.ndarray, np.ndarray]], current: dict[str, tuple[np.ndarray, np.ndarray]]) -> list[dict[str, float | str]]:
    rows: list[dict[str, float | str]] = []
    for series in SERIES:
        ref_values, ref_edges = ppg12[series.key]
        cur_values, cur_edges = current[series.key]
        ref_x, ref_y, _ = select_range(ref_values, ref_edges)
        cur_x, cur_y, _ = select_range(cur_values, cur_edges)
        ref_by_x = {round(float(x), 6): float(y) for x, y in zip(ref_x, ref_y)}
        cur_by_x = {round(float(x), 6): float(y) for x, y in zip(cur_x, cur_y)}
        common_x = sorted(set(ref_by_x) & set(cur_by_x))
        if not common_x:
            raise RuntimeError(f"No common PPG12/current bin centers for {series.key}")
        for idx, xkey in enumerate(common_x, start=1):
            xc = xkey
            rv = ref_by_x[xkey]
            cv = cur_by_x[xkey]
            ratio = rv / cv if cv > 0 else math.nan
            ratio_err = ratio * math.sqrt((1.0 / rv if rv > 0 else 0.0) + (1.0 / cv if cv > 0 else 0.0)) if math.isfinite(ratio) else math.nan
            rows.append(
                {
                    "series": series.key,
                    "label": series.label,
                    "bin_index_in_range": idx,
                    "xcenter": float(xc),
                    "ppg12_count": float(rv),
                    "current_count": float(cv),
                    "sdcc_over_current": float(ratio),
                    "ratio_stat_err": float(ratio_err),
                }
            )
    return rows


def write_csv(path: Path, rows: list[dict[str, float | str]]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def draw(
    path: Path,
    ppg12: dict[str, tuple[np.ndarray, np.ndarray]] | None,
    current: dict[str, tuple[np.ndarray, np.ndarray]],
    ratios: list[dict[str, float | str]],
    ratio_y: tuple[float, float],
    reference_label: str = "PPG12 SDCC",
    ratio_label: str = "SDCC / Current",
) -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 18,
            "axes.linewidth": 1.35,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.major.size": 8,
            "ytick.major.size": 8,
            "xtick.minor.size": 4,
            "ytick.minor.size": 4,
            "xtick.major.width": 1.1,
            "ytick.major.width": 1.1,
            "xtick.minor.width": 0.9,
            "ytick.minor.width": 0.9,
        }
    )
    fig = plt.figure(figsize=(8.0, 8.89), dpi=100)
    ax = fig.add_axes([0.13, 0.421, 0.79, 0.507])
    rax = fig.add_axes([0.13, 0.100, 0.79, 0.292], sharex=ax)

    for axis in (ax, rax):
        axis.tick_params(which="both", top=True, right=True, labelsize=15)
        axis.minorticks_on()

    ax.set_yscale("log")
    ax.set_xlim(10.0, 35.0)
    ax.set_ylim(7.0e3, 1.3e8)
    rax.set_ylim(*ratio_y)
    ax.tick_params(labelbottom=False)

    for series in SERIES:
        if ppg12 is not None:
            ref_values, ref_edges = ppg12[series.key]
            x, y, _ = select_range(ref_values, ref_edges)
            ax.errorbar(
                x,
                y,
                yerr=np.sqrt(np.clip(y, 0.0, None)),
                fmt=series.marker_open,
                markersize=5.5,
                markerfacecolor="white",
                markeredgecolor=series.color,
                markeredgewidth=1.15,
                ecolor=series.color,
                elinewidth=0.9,
                capsize=0,
                linestyle="none",
                zorder=series.zorder + 4,
            )
        cur_values, cur_edges = current[series.key]
        x, y, _ = select_range(cur_values, cur_edges)
        ax.errorbar(
            x,
            y,
            yerr=np.sqrt(np.clip(y, 0.0, None)),
            fmt=series.marker_closed,
            markersize=4.8,
            markerfacecolor=series.color,
            markeredgecolor=series.color,
            ecolor=series.color,
            elinewidth=0.75,
            capsize=0,
            linestyle="none",
            zorder=series.zorder,
        )

    for series in SERIES:
        rows = [row for row in ratios if row["series"] == series.key]
        if not rows:
            continue
        rax.errorbar(
            [float(row["xcenter"]) for row in rows],
            [float(row["sdcc_over_current"]) for row in rows],
            yerr=[float(row["ratio_stat_err"]) for row in rows],
            fmt=series.marker_closed,
            markersize=4.8,
            markerfacecolor=series.color,
            markeredgecolor=series.color,
            ecolor=series.color,
            elinewidth=0.8,
            capsize=0,
            linestyle="none",
            zorder=series.zorder,
        )

    rax.axhline(1.0, color="0.35", linestyle=(0, (4, 4)), linewidth=1.0)

    ax.set_ylabel("Counts", fontsize=18, labelpad=13)
    rax.set_ylabel(ratio_label, fontsize=18, labelpad=13)
    rax.set_xlabel(r"$p_T$ [GeV]", fontsize=18, labelpad=8)

    ax.text(0.065, 0.97, r"$\bf{\it{sPHENIX}}$ Simulation", transform=ax.transAxes, fontsize=13, va="top")
    ax.text(0.065, 0.92, "Photon Jet Samples", transform=ax.transAxes, fontsize=13, va="top")
    ax.text(0.25, 0.85, r"Pythia, $\sqrt{s}$=200 GeV", transform=ax.transAxes, fontsize=12, va="top")
    ax.text(0.25, 0.80, r"$|\eta^\gamma| < 0.7$", transform=ax.transAxes, fontsize=12, va="top")
    ax.text(0.25, 0.75, r"$R = 0.3,\ E_T^{iso} < 4$ GeV", transform=ax.transAxes, fontsize=12, va="top")

    sample_handles = [
        plt.Line2D([0], [0], marker=s.marker_closed, color=s.color, markerfacecolor=s.color, linestyle="none", markersize=6)
        for s in SERIES
    ]
    leg_sample = ax.legend(
        sample_handles,
        [s.label for s in SERIES],
        loc="upper left",
        bbox_to_anchor=(0.58, 0.82),
        frameon=False,
        fontsize=12,
        borderaxespad=0.0,
        labelspacing=0.35,
        handletextpad=0.5,
    )
    ax.add_artist(leg_sample)

    source_handles = [
        plt.Line2D([0], [0], marker="o", color="black", markerfacecolor="white", linestyle="none", markersize=6),
        plt.Line2D([0], [0], marker="o", color="black", markerfacecolor="black", linestyle="none", markersize=6),
    ]
    ax.legend(
        source_handles,
        [reference_label, "Current output"],
        loc="upper left",
        bbox_to_anchor=(0.58, 0.60),
        frameon=False,
        fontsize=11.5,
        borderaxespad=0.0,
        labelspacing=0.35,
        handletextpad=0.5,
    )

    fig.savefig(path)
    plt.close(fig)


def main() -> int:
    args = parse_args()
    current_root = args.current_root or current_root_from_pointer(args.current_pointer)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    current = load_current(current_root)
    if args.current_only_audit:
        ppg12 = None
        ppg12_source_kind = None
        ppg12_source_path = None
    elif args.ppg12_points_csv:
        ppg12 = load_ppg12_points_csv(args.ppg12_points_csv)
        ppg12_source_kind = "sdcc_root_projected_csv"
        ppg12_source_path = args.ppg12_points_csv
    elif args.ppg12_datathief_csv:
        ppg12 = load_ppg12_datathief(args.ppg12_datathief_csv)
        ppg12_source_kind = "datathief_csv"
        ppg12_source_path = args.ppg12_datathief_csv
    else:
        ppg12 = load_ppg12_reference(args.ppg12_root)
        ppg12_source_kind = "sdcc_root"
        ppg12_source_path = args.ppg12_root
    ratios = [] if ppg12 is None else ratio_rows(ppg12, current)

    points_csv = args.out_dir / "ppg12_fig3_truth_pt_sdcc_vs_current_points.csv"
    ratio_csv = args.out_dir / "ppg12_fig3_truth_pt_sdcc_over_current_ratios.csv" if ppg12 is not None else None
    manifest_path = args.out_dir / "ppg12_fig3_truth_pt_sdcc_vs_current_manifest.json"
    png = args.out_dir / "ppg12_fig3_truth_pt_sdcc_vs_current_overlay.png"

    rows = series_rows("current_output", current)
    if ppg12 is not None:
        rows = series_rows("ppg12_sdcc", ppg12) + rows
    write_csv(points_csv, rows)
    if ratio_csv is not None:
        write_csv(ratio_csv, ratios)
    reference_label = "PPG12 DataThief" if ppg12_source_kind == "datathief_csv" else "PPG12 SDCC"
    ratio_label = "DataThief / Current" if ppg12_source_kind == "datathief_csv" else "PPG12 SDCC / Current"
    draw(png, ppg12, current, ratios, (args.ratio_ymin, args.ratio_ymax), reference_label, ratio_label)

    manifest = {
        "artifact": "PPG12 Fig.3 truth-pT direct/fragmentation parity overlay",
        "ppg12_macro": "ppg12codeGit/efficiencytool/FindTruthETCut.C",
        "ppg12_draw_helper": "ppg12codeGit/efficiencytool/draw.C::draw_1D_multiple_plot_ratio",
        "ppg12_source_kind": ppg12_source_kind,
        "ppg12_source_root": str(args.ppg12_root) if ppg12_source_kind == "sdcc_root" else None,
        "ppg12_source_datathief_csv": str(ppg12_source_path) if ppg12_source_kind == "datathief_csv" else None,
        "ppg12_source_projected_csv": str(ppg12_source_path) if ppg12_source_kind == "sdcc_root_projected_csv" else None,
        "ppg12_source_remote": "/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_noiso.root",
        "ppg12_source_objects": [PPG12_DIRECT_H2, PPG12_FRAG_H2],
        "current_pointer": str(args.current_pointer),
        "current_root": str(current_root),
        "current_objects": CURRENT_HISTS,
        "current_merge_contract": "Resolved via pp_sim_photonjet_merged/current.json: registered photon+jet pp SIM current artifact produced by the PPG12 photon-only policy merge. This is the current downstream merged photon+jet output; the strict Fig.5 pre-vz source-stage stitch is a separate canonical source-diagnostic contract.",
        "algorithm": "PPG12 reference projects h_direct/h_frag over E_T^iso bins whose upper edge <= 4 GeV, applies Rebin(10) when the x-axis is 0.1 GeV, builds total=direct+frag, and draws 10<=pT<=35. Current output uses dedicated 1 GeV Fig.3 spectra from the registered current RecoilJets merged photon+jet artifact.",
        "marker_convention": "same colors/shapes by spectrum; open markers are PPG12 SDCC, filled markers are current output; lower pad is PPG12 SDCC / current output.",
        "plot_reference_label": reference_label,
        "plot_ratio_label": ratio_label,
        "png": str(png),
        "points_csv": str(points_csv),
        "ratio_csv": str(ratio_csv) if ratio_csv is not None else None,
        "canvas_contract": {
            "source": "PPG12 draw.C",
            "size_px_at_root": [800, 889],
            "top_pad": [0, 0.4, 1, 1],
            "bottom_pad": [0, 0, 1, 0.4],
            "margins": {
                "top_pad": {"top": 0.12, "left": 0.13, "bottom": 0.035, "right": 0.08},
                "bottom_pad": {"top": 0.02, "left": 0.13, "bottom": 0.25, "right": 0.08},
            },
        },
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "png": str(png),
                "points_csv": str(points_csv),
                "ratio_csv": str(ratio_csv) if ratio_csv is not None else None,
                "manifest": str(manifest_path),
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
