#!/usr/bin/env python3
"""Plot an offline pp photon+jet topo-vs-tower isolation diagnostic.

The locally available RecoilJets products are histogram-only.  The current
canonical pp photon+jet file contains the PPG12 raw topo-isolation spectra,
while the preserved May tower-cone product contains the tower-isolation
spectra.  They are not candidate-matched, so this script intentionally makes
distribution-level overlays and pT trends rather than a fake 2D correlation.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import uproot


REPO = Path(__file__).resolve().parents[3]
DEFAULT_TOPO_ROOT = (
    REPO
    / "dataOutput/current_recoiljets_artifacts/current/pp_sim_photonjet_merged/current.root"
)
DEFAULT_TOPO_POINTER = (
    REPO
    / "dataOutput/current_recoiljets_artifacts/current/pp_sim_photonjet_merged/current.json"
)
DEFAULT_TOWER_ROOT = (
    REPO
    / "dataOutput/combinedSimOnly"
    / "jetMinPt5_7pi_8_vz60_isoR40_fixedIso2GeV_"
    "preselectionReference_tightReference_nonTightReference"
    / "photonJet5and10and20merged_SIM"
    / "RecoilJets_photonjet5plus10plus20_MERGED.root"
)
DEFAULT_OUTDIR = (
    REPO / "dataOutput/ppg12Parity/the150_topocluster_vs_tower_offline_20260817"
)

PT_BINS = [
    (5.0, 8.0),
    (8.0, 10.0),
    (10.0, 12.0),
    (12.0, 14.0),
    (14.0, 16.0),
    (16.0, 18.0),
    (18.0, 20.0),
    (20.0, 22.0),
    (22.0, 24.0),
    (24.0, 26.0),
    (26.0, 35.0),
]
OVERLAY_BINS = {(5.0, 8.0), (10.0, 12.0), (16.0, 18.0), (26.0, 35.0)}


@dataclass(frozen=True)
class HistSummary:
    edges: np.ndarray
    values: np.ndarray
    variances: np.ndarray
    underflow: float
    overflow: float
    mean: float
    median: float
    std: float
    integral: float


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def pt_label(lo: float, hi: float) -> str:
    return f"{int(lo)}_{int(hi)}"


def binned_median(edges: np.ndarray, values: np.ndarray) -> float:
    total = float(np.sum(values))
    if not np.isfinite(total) or total <= 0:
        return float("nan")
    cumulative = np.cumsum(values)
    idx = int(np.searchsorted(cumulative, 0.5 * total, side="left"))
    idx = min(max(idx, 0), len(values) - 1)
    before = float(cumulative[idx - 1]) if idx > 0 else 0.0
    content = float(values[idx])
    fraction = 0.5 if content <= 0 else (0.5 * total - before) / content
    fraction = min(max(fraction, 0.0), 1.0)
    return float(edges[idx] + fraction * (edges[idx + 1] - edges[idx]))


def summarize_hist(hist: uproot.behaviors.TH1.Histogram) -> HistSummary:
    edges = np.asarray(hist.axis().edges(), dtype=float)
    values_flow = np.asarray(hist.values(flow=True), dtype=float)
    variances_flow = hist.variances(flow=True)
    if variances_flow is None:
        variances_flow = np.maximum(values_flow, 0.0)
    variances_flow = np.asarray(variances_flow, dtype=float)
    values = values_flow[1:-1]
    variances = variances_flow[1:-1]
    centers = 0.5 * (edges[:-1] + edges[1:])
    integral = float(np.sum(values))
    if integral > 0:
        mean = float(np.sum(values * centers) / integral)
        variance = float(np.sum(values * (centers - mean) ** 2) / integral)
        std = float(np.sqrt(max(variance, 0.0)))
    else:
        mean = std = float("nan")
    return HistSummary(
        edges=edges,
        values=values,
        variances=variances,
        underflow=float(values_flow[0]),
        overflow=float(values_flow[-1]),
        mean=mean,
        median=binned_median(edges, values),
        std=std,
        integral=integral,
    )


def load_spectra(
    root_path: Path,
    base: str,
) -> dict[tuple[float, float], HistSummary]:
    spectra: dict[tuple[float, float], HistSummary] = {}
    with uproot.open(root_path) as root_file:
        for lo, hi in PT_BINS:
            key = f"SIM/{base}_pT_{pt_label(lo, hi)}"
            if key not in root_file:
                raise KeyError(f"missing required histogram {key} in {root_path}")
            spectra[(lo, hi)] = summarize_hist(root_file[key])
    return spectra


def configure_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "mathtext.fontset": "dejavusans",
            "axes.linewidth": 1.15,
            "axes.labelsize": 13,
            "axes.titlesize": 13,
            "xtick.labelsize": 10.5,
            "ytick.labelsize": 10.5,
            "legend.fontsize": 10.5,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )


def density(summary: HistSummary) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    widths = np.diff(summary.edges)
    norm = summary.integral
    if norm <= 0:
        return summary.values.copy(), summary.values.copy(), widths
    values = summary.values / (norm * widths)
    errors = np.sqrt(np.maximum(summary.variances, 0.0)) / (norm * widths)
    return values, errors, widths


def make_overlay_plot(
    tower: dict[tuple[float, float], HistSummary],
    topo: dict[tuple[float, float], HistSummary],
    output: Path,
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12.6, 9.3), dpi=190, sharex=True, sharey=True)
    for ax, pt_bin in zip(axes.flat, sorted(OVERLAY_BINS)):
        lo, hi = pt_bin
        tower_hist = tower[pt_bin]
        topo_hist = topo[pt_bin]
        if not np.allclose(tower_hist.edges, topo_hist.edges):
            raise ValueError(f"incompatible isolation binning for pT {lo:g}-{hi:g} GeV")
        centers = 0.5 * (tower_hist.edges[:-1] + tower_hist.edges[1:])
        tower_density, tower_error, _ = density(tower_hist)
        topo_density, topo_error, _ = density(topo_hist)
        ax.step(
            centers,
            tower_density,
            where="mid",
            color="#0072B2",
            linewidth=1.8,
            label="Tower cone, R=0.4",
        )
        ax.fill_between(
            centers,
            np.maximum(tower_density - tower_error, 1.0e-8),
            tower_density + tower_error,
            step="mid",
            color="#0072B2",
            alpha=0.16,
            linewidth=0,
        )
        ax.step(
            centers,
            topo_density,
            where="mid",
            color="#D55E00",
            linewidth=1.8,
            label="Raw topo clusters, R=0.4",
        )
        ax.fill_between(
            centers,
            np.maximum(topo_density - topo_error, 1.0e-8),
            topo_density + topo_error,
            step="mid",
            color="#D55E00",
            alpha=0.16,
            linewidth=0,
        )
        ax.axvline(0.0, color="0.35", linestyle=":", linewidth=1.0)
        ax.set_yscale("log")
        ax.set_xlim(-5.0, 12.0)
        ax.set_ylim(2.0e-5, 3.0)
        ax.set_title(rf"${lo:g} < p_T^{{\gamma}} < {hi:g}$ GeV", loc="left")
        ax.grid(axis="y", which="both", alpha=0.16)
    axes[0, 0].legend(frameon=False, loc="upper right")
    for ax in axes[:, 0]:
        ax.set_ylabel(r"Unit-normalized candidates / GeV")
    for ax in axes[1, :]:
        ax.set_xlabel(r"Raw $E_T^{\mathrm{iso}}(R=0.4)$ [GeV]")
    fig.suptitle(
        r"$p{+}p$ photon+jet simulation: tower-cone and topo-cluster isolation",
        fontsize=17,
        y=0.985,
    )
    fig.text(
        0.017,
        0.014,
        r"$\bf{\it{sPHENIX}}$ Internal   $\sqrt{s}=200$ GeV",
        ha="left",
        va="bottom",
        fontsize=11.5,
    )
    fig.text(
        0.983,
        0.014,
        "Offline diagnostic: different production snapshots; not candidate-matched",
        ha="right",
        va="bottom",
        fontsize=10.5,
        color="0.30",
    )
    fig.subplots_adjust(left=0.09, right=0.985, bottom=0.09, top=0.92, hspace=0.10, wspace=0.08)
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def make_summary_plot(
    tower: dict[tuple[float, float], HistSummary],
    topo: dict[tuple[float, float], HistSummary],
    output: Path,
) -> None:
    bins = PT_BINS
    centers = np.asarray([(lo + hi) / 2.0 for lo, hi in bins])
    xerr = np.asarray([(hi - lo) / 2.0 for lo, hi in bins])
    tower_mean = np.asarray([tower[b].mean for b in bins])
    topo_mean = np.asarray([topo[b].mean for b in bins])
    tower_median = np.asarray([tower[b].median for b in bins])
    topo_median = np.asarray([topo[b].median for b in bins])

    fig, (ax, dax) = plt.subplots(
        2,
        1,
        figsize=(11.4, 8.1),
        dpi=190,
        sharex=True,
        gridspec_kw={"height_ratios": [2.0, 1.0], "hspace": 0.06},
    )
    ax.errorbar(
        centers - 0.08,
        tower_mean,
        xerr=xerr,
        fmt="o-",
        markersize=5.5,
        linewidth=1.65,
        color="#0072B2",
        label="Tower mean",
    )
    ax.errorbar(
        centers + 0.08,
        topo_mean,
        xerr=xerr,
        fmt="o-",
        markersize=5.5,
        linewidth=1.65,
        color="#D55E00",
        label="Topo mean",
    )
    ax.errorbar(
        centers - 0.08,
        tower_median,
        xerr=xerr,
        fmt="s--",
        markersize=5.0,
        linewidth=1.45,
        color="#0072B2",
        markerfacecolor="white",
        label="Tower median",
    )
    ax.errorbar(
        centers + 0.08,
        topo_median,
        xerr=xerr,
        fmt="s--",
        markersize=5.0,
        linewidth=1.45,
        color="#D55E00",
        markerfacecolor="white",
        label="Topo median",
    )
    ax.set_ylabel(r"Isolation summary [GeV]")
    ax.set_ylim(-1.0, 4.2)
    ax.grid(axis="y", alpha=0.22)
    ax.legend(frameon=False, ncol=2, loc="upper left")
    ax.text(
        0.985,
        0.96,
        r"$\bf{\it{sPHENIX}}$ Internal" + "\n" + r"$p{+}p$, $\sqrt{s}=200$ GeV",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=11.5,
    )

    dax.axhline(0.0, color="0.25", linewidth=1.1)
    dax.errorbar(
        centers,
        topo_mean - tower_mean,
        xerr=xerr,
        fmt="o-",
        color="#6A3D9A",
        markersize=5.3,
        linewidth=1.6,
        label="Mean difference",
    )
    dax.errorbar(
        centers,
        topo_median - tower_median,
        xerr=xerr,
        fmt="s--",
        color="#009E73",
        markerfacecolor="white",
        markersize=5.0,
        linewidth=1.45,
        label="Median difference",
    )
    dax.set_ylabel("Topo - tower [GeV]")
    dax.set_xlabel(r"Photon $p_T^{\gamma}$ [GeV]")
    dax.set_ylim(-4.0, 2.0)
    dax.grid(axis="y", alpha=0.22)
    dax.legend(frameon=False, ncol=2, loc="lower right")
    dax.text(
        0.015,
        0.08,
        "Difference of distribution summaries, not event-level $\\Delta E_T^{iso}$",
        transform=dax.transAxes,
        ha="left",
        va="bottom",
        fontsize=10.2,
        color="0.30",
    )
    fig.suptitle(
        r"Raw $R=0.4$ isolation versus photon $p_T$: offline pp photon+jet diagnostic",
        fontsize=16.5,
        y=0.985,
    )
    fig.subplots_adjust(left=0.11, right=0.985, bottom=0.10, top=0.91)
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def write_csv(
    tower: dict[tuple[float, float], HistSummary],
    topo: dict[tuple[float, float], HistSummary],
    output: Path,
) -> None:
    fields = [
        "pt_lo",
        "pt_hi",
        "tower_mean",
        "tower_median",
        "tower_std",
        "tower_integral_inrange",
        "tower_underflow",
        "tower_overflow",
        "topo_mean",
        "topo_median",
        "topo_std",
        "topo_integral_inrange",
        "topo_underflow",
        "topo_overflow",
        "topo_minus_tower_mean",
        "topo_minus_tower_median",
    ]
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for lo, hi in PT_BINS:
            t = tower[(lo, hi)]
            c = topo[(lo, hi)]
            writer.writerow(
                {
                    "pt_lo": lo,
                    "pt_hi": hi,
                    "tower_mean": t.mean,
                    "tower_median": t.median,
                    "tower_std": t.std,
                    "tower_integral_inrange": t.integral,
                    "tower_underflow": t.underflow,
                    "tower_overflow": t.overflow,
                    "topo_mean": c.mean,
                    "topo_median": c.median,
                    "topo_std": c.std,
                    "topo_integral_inrange": c.integral,
                    "topo_underflow": c.underflow,
                    "topo_overflow": c.overflow,
                    "topo_minus_tower_mean": c.mean - t.mean,
                    "topo_minus_tower_median": c.median - t.median,
                }
            )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--topo-root", type=Path, default=DEFAULT_TOPO_ROOT)
    parser.add_argument("--topo-pointer", type=Path, default=DEFAULT_TOPO_POINTER)
    parser.add_argument("--tower-root", type=Path, default=DEFAULT_TOWER_ROOT)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = parser.parse_args()

    for path in (args.topo_root, args.topo_pointer, args.tower_root):
        if not path.exists():
            raise FileNotFoundError(path)
    args.outdir.mkdir(parents=True, exist_ok=True)
    configure_style()

    tower = load_spectra(args.tower_root, "h_Eiso")
    topo = load_spectra(args.topo_root, "h_Eiso_ppg12_topo_raw")

    overlay_png = args.outdir / "pp_photonjet_topo_vs_tower_iso_ptbin_overlays.png"
    summary_png = args.outdir / "pp_photonjet_topo_vs_tower_iso_summary_vs_pt.png"
    csv_path = args.outdir / "pp_photonjet_topo_vs_tower_iso_summary.csv"
    manifest_path = args.outdir / "pp_photonjet_topo_vs_tower_iso_manifest.json"

    make_overlay_plot(tower, topo, overlay_png)
    make_summary_plot(tower, topo, summary_png)
    write_csv(tower, topo, csv_path)

    pointer = json.loads(args.topo_pointer.read_text())
    manifest = {
        "ok": True,
        "artifact_status": "diagnostic_only_not_candidate_matched",
        "comparison": "raw R=0.4 topo-cluster versus tower-cone isolation",
        "sample": "pp photon+jet simulation, photon5+10+20",
        "topo_input": {
            "pointer": str(args.topo_pointer),
            "root": str(args.topo_root.resolve()),
            "sha256": sha256(args.topo_root.resolve()),
            "campaign_tag": pointer.get("campaign_tag"),
            "histogram": "SIM/h_Eiso_ppg12_topo_raw_pT_<lo>_<hi>",
        },
        "tower_input": {
            "root": str(args.tower_root.resolve()),
            "sha256": sha256(args.tower_root.resolve()),
            "histogram": "SIM/h_Eiso_pT_<lo>_<hi>",
            "production_artifact_exception": (
                "The current topo-mode ROOT does not retain tower-cone spectra. "
                "This preserved historical pp photon+jet tower-mode product is used "
                "for a shape-only offline diagnostic."
            ),
        },
        "selection_scope": (
            "Reco photon candidates with pT>=5 GeV and |eta| inside the configured "
            "acceptance, before isolation and tight-ID requirements."
        ),
        "normalization": "unit area independently in each pT bin for overlays",
        "known_limitations": [
            "The two histograms come from different production snapshots.",
            "The final ROOT files are histogram-only and do not permit same-candidate 2D correlation or event-level delta isolation.",
            "Mean and median differences are differences between marginal distributions, not paired-candidate differences.",
            "The histogram-visible isolation range is -5 to 12 GeV; underflow and overflow are recorded in the CSV.",
        ],
        "outputs": [str(overlay_png), str(summary_png), str(csv_path)],
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")

    for path in (overlay_png, summary_png, csv_path, manifest_path):
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
