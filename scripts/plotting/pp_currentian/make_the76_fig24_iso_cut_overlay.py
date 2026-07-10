#!/usr/bin/env python3
"""Build the THE-76 PPG12 Fig.24 isolation-cut overlay.

This reads the dedicated current-analysis Fig.24 TH2 filled from the
truth-matched signal loop and reproduces the PPG12 FindETCut.C cutoff scan.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import uproot


REPO = Path(__file__).resolve().parents[3]
DEFAULT_REFERENCE = (
    REPO
    / "dataOutput/ppg12Parity/the76_ppg12_fig24_photonjet_fix_20260702_014217"
    / "fig24_iso_efficiency/fig24_ppg12_cutoff_points_from_sdcc_nom_root.csv"
)

PPG12_FITS = {
    0.70: (0.312, 0.0311),
    0.80: (0.490, 0.0370),
    0.90: (0.919, 0.0448),
}

COLORS = {
    0.70: "#169df2",
    0.80: "#2ca02c",
    0.90: "#d6279f",
}

MARKERS = {
    0.70: "s",
    0.80: "^",
    0.90: "o",
}


def root_find_bin(edges: np.ndarray, value: float) -> int:
    """Return a 1-based ROOT-like bin number, clamped to normal bins."""
    nbins = len(edges) - 1
    if value < edges[0]:
        return 1
    if value >= edges[-1]:
        return nbins
    return int(np.searchsorted(edges, value, side="right"))


def fit_value(eff: float, x: float) -> float:
    a, b = PPG12_FITS[eff]
    return a + b * x


def find_eff_cutoff(
    values: np.ndarray,
    xedges: np.ndarray,
    yedges: np.ndarray,
    eff: float,
    pt_low: float,
    pt_high: float,
    iso_low: float = -1.0,
) -> tuple[float, float, float]:
    """Replicate PPG12 FindETCut.C for one pT interval."""
    xlo = root_find_bin(xedges, pt_low)
    xhi = root_find_bin(xedges, pt_high)
    xlo = max(1, min(values.shape[0], xlo))
    xhi = max(1, min(values.shape[0], xhi))
    if xhi < xlo:
        xlo, xhi = xhi, xlo

    block = values[xlo - 1 : xhi, :]
    total = float(np.nansum(block))
    yerr = float((yedges[1] - yedges[0]) / 2.0)
    if total <= 0:
        return math.nan, yerr, total

    target = total * eff
    ylo = root_find_bin(yedges, iso_low)
    ylo = max(1, min(values.shape[1], ylo))
    current = 0.0
    found_y = values.shape[1]
    for ybin in range(ylo, values.shape[1] + 1):
        current += float(np.nansum(block[:, ybin - 1]))
        if current >= target:
            found_y = ybin
            break
    return float(yedges[found_y]), yerr, total


def load_ppg12_points(path: Path) -> list[dict[str, float]]:
    points: list[dict[str, float]] = []
    lines = path.read_text().splitlines()
    try:
        start = lines.index("POINTS_BEGIN") + 1
    except ValueError as exc:
        raise RuntimeError(f"{path} does not contain POINTS_BEGIN") from exc
    reader = csv.DictReader(lines[start:])
    for row in reader:
        if row.get("kind") != "POINT":
            continue
        cut = float(row["cut"])
        if cut <= 0:
            continue
        eff = float(row["eff"])
        xcenter = float(row["xcenter"])
        points.append(
            {
                "source": "ppg12_sdcc_points",
                "eff": eff,
                "bin": int(row["bin"]),
                "xlow": float(row["xlow"]),
                "xhigh": float(row["xhigh"]),
                "xcenter": xcenter,
                "cut": cut,
                "yerr": float(row["yerr"]),
                "ratio_to_ppg12_fit": cut / fit_value(eff, xcenter),
                "integral": math.nan,
            }
        )
    return points


def load_current_hist(root_path: Path, hist_names: list[str]) -> tuple[str, np.ndarray, np.ndarray, np.ndarray]:
    with uproot.open(root_path) as handle:
        available = set(handle.keys(cycle=False))
        for hist_name in hist_names:
            if hist_name in available:
                hist = handle[hist_name]
                values, xedges, yedges = hist.to_numpy(flow=False)
                return hist_name, values.astype(float), xedges.astype(float), yedges.astype(float)
    raise RuntimeError(
        "None of the required Fig24 histograms were found in "
        f"{root_path}: {', '.join(hist_names)}"
    )


def derive_current_points(values: np.ndarray, xedges: np.ndarray, yedges: np.ndarray) -> list[dict[str, float]]:
    points: list[dict[str, float]] = []
    pt_edges = np.linspace(10.0, 36.0, 14)
    for eff in (0.70, 0.80, 0.90):
        for idx in range(13):
            xlow = float(pt_edges[idx])
            xhigh = float(pt_edges[idx + 1])
            xcenter = 0.5 * (xlow + xhigh)
            cut, yerr, integral = find_eff_cutoff(values, xedges, yedges, eff, xlow, xhigh)
            points.append(
                {
                    "source": "current_dedicated_fig24_th2",
                    "eff": eff,
                    "bin": idx + 1,
                    "xlow": xlow,
                    "xhigh": xhigh,
                    "xcenter": xcenter,
                    "cut": cut,
                    "yerr": yerr,
                    "ratio_to_ppg12_fit": cut / fit_value(eff, xcenter) if math.isfinite(cut) else math.nan,
                    "integral": integral,
                }
            )
    return points


def write_points(path: Path, rows: list[dict[str, float]]) -> None:
    fieldnames = [
        "source",
        "eff",
        "bin",
        "xlow",
        "xhigh",
        "xcenter",
        "cut",
        "yerr",
        "ratio_to_ppg12_fit",
        "integral",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def draw_plot(path: Path, ppg12_points: list[dict[str, float]], current_points: list[dict[str, float]]) -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Serif",
            "font.size": 18,
            "axes.linewidth": 1.4,
            "xtick.major.size": 8,
            "xtick.minor.size": 4,
            "ytick.major.size": 8,
            "ytick.minor.size": 4,
            "xtick.direction": "in",
            "ytick.direction": "in",
        }
    )
    fig, (ax, rax) = plt.subplots(
        2,
        1,
        figsize=(7.4, 7.2),
        gridspec_kw={"height_ratios": [3.0, 1.35], "hspace": 0.03},
        sharex=True,
    )
    xfit = np.linspace(10.0, 36.0, 300)
    ratio_rows: list[dict[str, float]] = []

    for eff in (0.70, 0.80, 0.90):
        color = COLORS[eff]
        marker = MARKERS[eff]
        ref = [p for p in ppg12_points if abs(p["eff"] - eff) < 1e-9]
        cur = [p for p in current_points if abs(p["eff"] - eff) < 1e-9 and math.isfinite(p["cut"])]
        ax.plot(xfit, [fit_value(eff, x) for x in xfit], color=color, linestyle="--", linewidth=1.6, zorder=1)
        if ref:
            ax.errorbar(
                [p["xcenter"] for p in ref],
                [p["cut"] for p in ref],
                yerr=[p["yerr"] for p in ref],
                fmt=marker,
                markersize=6.2,
                markerfacecolor="none",
                markeredgecolor=color,
                markeredgewidth=1.35,
                ecolor=color,
                elinewidth=0.95,
                capsize=0,
                linestyle="none",
                zorder=7,
            )
        if cur:
            ax.errorbar(
                [p["xcenter"] for p in cur],
                [p["cut"] for p in cur],
                yerr=[p["yerr"] for p in cur],
                fmt=marker,
                markersize=5.2,
                markerfacecolor=color,
                markeredgecolor=color,
                ecolor=color,
                elinewidth=1.0,
                capsize=0,
                linestyle="none",
                zorder=4,
            )
        ref_by_x = {p["xcenter"]: p for p in ref}
        cur_by_x = {p["xcenter"]: p for p in cur}
        for x in sorted(set(ref_by_x) & set(cur_by_x)):
            rp = ref_by_x[x]
            cp = cur_by_x[x]
            if cp["cut"] <= 0 or not math.isfinite(cp["cut"]):
                continue
            ratio = rp["cut"] / cp["cut"]
            ratio_err = ratio * math.sqrt(
                (rp["yerr"] / rp["cut"]) ** 2 + (cp["yerr"] / cp["cut"]) ** 2
            )
            ratio_rows.append(
                {
                    "eff": eff,
                    "xcenter": x,
                    "ratio": ratio,
                    "ratio_err": ratio_err,
                    "marker": marker,
                    "color": color,
                }
            )

    for eff in (0.70, 0.80, 0.90):
        rows = [r for r in ratio_rows if abs(r["eff"] - eff) < 1e-9]
        if not rows:
            continue
        rax.errorbar(
            [r["xcenter"] for r in rows],
            [r["ratio"] for r in rows],
            yerr=[r["ratio_err"] for r in rows],
            fmt=MARKERS[eff],
            markersize=5.2,
            markerfacecolor=COLORS[eff],
            markeredgecolor=COLORS[eff],
            ecolor=COLORS[eff],
            elinewidth=0.9,
            capsize=0,
            linestyle="none",
            zorder=5,
        )

    ax.set_xlim(10.0, 36.0)
    ax.set_ylim(0.0, 4.5)
    rax.set_ylim(0.84, 1.12)
    ax.set_ylabel(r"$E_T^{iso}$ cutoff [GeV]", labelpad=8)
    rax.set_ylabel("SDCC / Current", labelpad=8)
    rax.set_xlabel(r"Cluster $p_T$ [GeV]")
    rax.axhline(1.0, color="0.45", linestyle=(0, (4, 4)), linewidth=1.0, zorder=0)

    for axis in (ax, rax):
        axis.tick_params(which="both", top=True, right=True)
        axis.minorticks_on()

    ax.text(0.10, 0.93, r"$\it{\bf{sPHENIX}}$ Internal", transform=ax.transAxes, fontsize=16)
    ax.text(0.10, 0.86, r"$p{+}p\ \sqrt{s}=200$ GeV", transform=ax.transAxes, fontsize=15)
    ax.text(0.10, 0.79, r"$|\eta^\gamma| < 0.7$", transform=ax.transAxes, fontsize=15)
    ax.text(0.88, 0.93, "PYTHIA Signal", transform=ax.transAxes, ha="right", fontsize=15)

    legend_handles = []
    legend_labels = []
    for eff in (0.70, 0.80, 0.90):
        a, b = PPG12_FITS[eff]
        handle = plt.Line2D(
            [0],
            [0],
            color=COLORS[eff],
            marker=MARKERS[eff],
            markerfacecolor=COLORS[eff],
            markeredgecolor=COLORS[eff],
            linestyle="--",
            linewidth=1.4,
            markersize=7,
        )
        legend_handles.append(handle)
        legend_labels.append(f"{int(100*eff)}%: {a:.3f} + {b:.4f}  $p_T$")
    leg1 = ax.legend(
        legend_handles,
        legend_labels,
        loc="upper left",
        bbox_to_anchor=(0.10, 0.72),
        frameon=False,
        fontsize=11.5,
        handlelength=2.2,
        borderaxespad=0.0,
        labelspacing=0.35,
    )
    ax.add_artist(leg1)

    source_handles = [
        plt.Line2D(
            [0],
            [0],
            color="black",
            marker="o",
            markerfacecolor="white",
            markeredgecolor="black",
            linestyle="none",
            markersize=6,
        ),
        plt.Line2D(
            [0],
            [0],
            color="black",
            marker="o",
            markerfacecolor="black",
            markeredgecolor="black",
            linestyle="none",
            markersize=6,
        ),
    ]
    ax.legend(
        source_handles,
        ["PPG12 SDCC source", "Current output"],
        loc="upper right",
        bbox_to_anchor=(0.985, 0.72),
        frameon=False,
        fontsize=10.8,
        borderaxespad=0.0,
        labelspacing=0.35,
        handletextpad=0.5,
    )

    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--current-root", required=True, type=Path)
    parser.add_argument("--campaign-tag", default="the76_ppg12_fig24_photonjet_fix_20260702_014217")
    parser.add_argument("--reference-csv", default=DEFAULT_REFERENCE, type=Path)
    parser.add_argument("--out-dir", type=Path)
    parser.add_argument(
        "--hist",
        action="append",
        default=["SIM/h_singal_reco_isoET_0", "SIM/h_ppg12_fig24_signal_reco_isoET_eta0"],
        help="Current ROOT histogram to use. Can be supplied multiple times.",
    )
    args = parser.parse_args()

    out_dir = args.out_dir or REPO / "dataOutput/ppg12Parity" / args.campaign_tag / "fig24_iso_efficiency"
    out_dir.mkdir(parents=True, exist_ok=True)
    png = out_dir / "fig24_cutoff_ppg12_sdcc_vs_current_dedicated_fig24_overlay.png"
    points_csv = out_dir / "fig24_cutoff_ppg12_sdcc_vs_current_dedicated_fig24_overlay_points.csv"
    manifest_path = out_dir / "fig24_cutoff_ppg12_sdcc_vs_current_dedicated_fig24_overlay_manifest.json"

    hist_name, values, xedges, yedges = load_current_hist(args.current_root, args.hist)
    ppg12_points = load_ppg12_points(args.reference_csv)
    current_points = derive_current_points(values, xedges, yedges)
    rows = ppg12_points + current_points
    write_points(points_csv, rows)
    draw_plot(png, ppg12_points, current_points)

    missing_bins = [
        {
            "eff": p["eff"],
            "bin": p["bin"],
            "xlow": p["xlow"],
            "xhigh": p["xhigh"],
        }
        for p in current_points
        if not math.isfinite(p["cut"]) or p["integral"] <= 0
    ]
    manifest = {
        "artifact": "THE-76 dedicated PPG12 Fig.24 isolation cutoff overlay",
        "png": str(png),
        "points_csv": str(points_csv),
        "current_root": str(args.current_root),
        "current_histogram": hist_name,
        "reference_csv": str(args.reference_csv),
        "algorithm": "Replicates PPG12 efficiencytool/FindETCut.C with 13 pT bins from 10 to 36 GeV and isolow=-1.0.",
        "ppg12_fit_equations": {f"{eff:.2f}": {"a": a, "b": b} for eff, (a, b) in PPG12_FITS.items()},
        "marker_convention": "open markers = PPG12 SDCC source points; filled markers = current output from dedicated Fig24 TH2.",
        "ratio_panel": "PPG12 SDCC source cutoff divided by current-output cutoff for bins where both sources exist.",
        "missing_current_bins": missing_bins,
        "current_integral_summary": {
            f"{eff:.2f}": {
                "min_integral": min(p["integral"] for p in current_points if abs(p["eff"] - eff) < 1e-9),
                "max_integral": max(p["integral"] for p in current_points if abs(p["eff"] - eff) < 1e-9),
            }
            for eff in PPG12_FITS
        },
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"png": str(png), "points_csv": str(points_csv), "manifest": str(manifest_path)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
