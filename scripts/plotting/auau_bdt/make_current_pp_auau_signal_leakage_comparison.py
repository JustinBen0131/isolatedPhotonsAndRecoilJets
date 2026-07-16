#!/usr/bin/env python3
"""Compare current pp and AuAu truth-signal ABCD leakage fractions."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


REPO = Path(__file__).resolve().parents[3]
DEFAULT_PP_ROOT = (
    REPO
    / "InputFiles/the97_ppg12_final_accepted_triple_full_20260714_1550"
    / "final_merged_roots/final_combined_canonical_20260715"
    / "RecoilJets_photonjet5plus10plus20_si_di_period_combined_MERGED.root"
)
DEFAULT_AUAU_ROOT = (
    REPO
    / "InputFiles/the88_bounded_nontight_20260708/simembedded"
    / "preselectionNewPPG12_tightAuAuCentInputBase3x3BDT_nonTightAuAuBDTSideband_baseVariant"
    / "photonJet12and20merged_SIM/RecoilJets_embeddedPhoton12plus20_MERGED.root"
)
DEFAULT_OUT_DIR = (
    REPO
    / "dataOutput/auau/the88_auau_data_centrality_fixed_20260713"
    / "signal_leakage_current_pp_comparison"
)

PP_BINS = [(16, 18), (18, 20), (20, 22), (22, 24), (24, 26), (26, 35)]
AUAU_BINS = [(15, 17), (17, 19), (19, 21), (21, 23), (23, 26), (26, 35)]
CENTRALITIES = [("0_20", "0–20%", "#0284c7", "s"), ("20_50", "20–50%", "#f59e0b", "^"), ("50_80", "50–80%", "#15803d", "D")]
REGIONS = [
    ("B", 2, r"$N_B^{\mathrm{sig}}/N_A^{\mathrm{sig}}$", "tight, non-isolated"),
    ("C", 3, r"$N_C^{\mathrm{sig}}/N_A^{\mathrm{sig}}$", "non-tight, isolated"),
    ("D", 4, r"$N_D^{\mathrm{sig}}/N_A^{\mathrm{sig}}$", "non-tight, non-isolated"),
]

PP_COLOR = "#7c3aed"
INK = "#111827"
MUTED = "#475569"
GRID = "#dbe4ef"
PANEL = "#fbfdff"


@dataclass(frozen=True)
class Point:
    sample: str
    sample_label: str
    centrality: str
    region: str
    et_low_gev: float
    et_high_gev: float
    n_a_signal: float
    n_a_error: float
    n_side_signal: float
    n_side_error: float
    leakage_ratio: float
    leakage_error: float
    histogram: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pp-root", type=Path, default=DEFAULT_PP_ROOT)
    parser.add_argument("--auau-root", type=Path, default=DEFAULT_AUAU_ROOT)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for chunk in iter(lambda: source.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def open_root(path: Path):
    import ROOT

    ROOT.gROOT.SetBatch(True)
    ROOT.TH1.AddDirectory(False)
    root_file = ROOT.TFile.Open(str(path), "READ")
    if not root_file or root_file.IsZombie() or root_file.TestBit(ROOT.TFile.kRecovered):
        raise OSError(f"ROOT integrity failure: {path}")
    return root_file


def ratio(num: float, num_error: float, den: float, den_error: float) -> tuple[float, float]:
    if den <= 0.0:
        raise ValueError("Region-A signal denominator is not positive")
    value = num / den
    variance = (num_error / den) ** 2 + ((num * den_error) / (den * den)) ** 2
    error = math.sqrt(max(0.0, variance))
    if not (math.isfinite(value) and math.isfinite(error) and value >= 0.0):
        raise ValueError(f"Invalid leakage ratio: {value} +/- {error}")
    return value, error


def extract_hist_point(
    root_file,
    *,
    path: str,
    sample: str,
    sample_label: str,
    centrality: str,
    region: str,
    region_bin: int,
    et_low: int,
    et_high: int,
) -> Point:
    hist = root_file.Get(path)
    if not hist:
        raise KeyError(f"Missing histogram: {path}")
    if hist.GetNbinsX() < 4:
        raise ValueError(f"ABCD histogram has fewer than four bins: {path}")
    n_a = float(hist.GetBinContent(1))
    n_a_error = float(hist.GetBinError(1))
    n_side = float(hist.GetBinContent(region_bin))
    n_side_error = float(hist.GetBinError(region_bin))
    if n_a_error <= 0.0 and n_a > 0.0:
        n_a_error = math.sqrt(n_a)
    if n_side_error <= 0.0 and n_side > 0.0:
        n_side_error = math.sqrt(n_side)
    leakage, leakage_error = ratio(n_side, n_side_error, n_a, n_a_error)
    return Point(
        sample=sample,
        sample_label=sample_label,
        centrality=centrality,
        region=region,
        et_low_gev=float(et_low),
        et_high_gev=float(et_high),
        n_a_signal=n_a,
        n_a_error=n_a_error,
        n_side_signal=n_side,
        n_side_error=n_side_error,
        leakage_ratio=leakage,
        leakage_error=leakage_error,
        histogram=path,
    )


def collect(pp_root: Path, auau_root: Path) -> list[Point]:
    pp_file = open_root(pp_root)
    auau_file = open_root(auau_root)
    points: list[Point] = []
    try:
        for region, region_bin, _ylabel, _description in REGIONS:
            for low, high in PP_BINS:
                points.append(
                    extract_hist_point(
                        pp_file,
                        path=f"SIM/h_sigABCD_MC_pT_{low}_{high}",
                        sample="pp",
                        sample_label="p+p",
                        centrality="inclusive",
                        region=region,
                        region_bin=region_bin,
                        et_low=low,
                        et_high=high,
                    )
                )
            for cent_key, cent_label, _color, _marker in CENTRALITIES:
                for low, high in AUAU_BINS:
                    points.append(
                        extract_hist_point(
                            auau_file,
                            path=(
                                "SIM/h_sigABCD_MC_isoR40_isSliding_"
                                f"pT_{low}_{high}_cent_{cent_key}"
                            ),
                            sample="auau",
                            sample_label="Au+Au embedded photon+jet 12 + 20",
                            centrality=cent_label,
                            region=region,
                            region_bin=region_bin,
                            et_low=low,
                            et_high=high,
                        )
                    )
    finally:
        pp_file.Close()
        auau_file.Close()
    return points


def series(points: list[Point], region: str, sample: str, centrality: str) -> list[Point]:
    return sorted(
        (
            point
            for point in points
            if point.region == region and point.sample == sample and point.centrality == centrality
        ),
        key=lambda point: point.et_low_gev,
    )


def draw(points: list[Point], output: Path) -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "mathtext.fontset": "dejavusans",
            "axes.unicode_minus": False,
        }
    )
    fig, axes = plt.subplots(3, 1, figsize=(14.5, 10.3), sharex=True)
    fig.subplots_adjust(left=0.105, right=0.975, bottom=0.105, top=0.79, hspace=0.30)
    fig.suptitle(
        "Truth-signal ABCD leakage fractions in photon+jet simulation",
        fontsize=25,
        fontweight="bold",
        y=0.965,
        color=INK,
    )
    fig.text(
        0.105,
        0.915,
        "sPHENIX",
        fontsize=18,
        fontstyle="italic",
        fontweight="bold",
        color=INK,
        ha="left",
        va="center",
    )
    fig.text(
        0.188,
        0.915,
        "Internal",
        fontsize=18,
        color=INK,
        ha="left",
        va="center",
    )
    fig.text(
        0.975,
        0.915,
        "sideband truth signal relative to Region A",
        fontsize=15.5,
        color=MUTED,
        ha="right",
        va="center",
    )
    fig.text(
        0.5,
        0.872,
        r"p+p photon+jet 5 + 10 + 20 and Au+Au embedded photon+jet 12 + 20, $\sqrt{s_{NN}}=200$ GeV",
        fontsize=15.5,
        color=MUTED,
        ha="center",
        va="center",
    )
    fig.text(
        0.5,
        0.837,
        r"$R=0.4$ sliding isolation; horizontal error bars show each production's native stored $E_T$ interval",
        fontsize=14.5,
        color=MUTED,
        ha="center",
        va="center",
    )

    for axis, (region, _region_bin, ylabel, description) in zip(axes, REGIONS):
        axis.set_facecolor(PANEL)
        axis.set_title(description, loc="left", fontsize=16.5, fontweight="bold", pad=7, color=INK)
        axis.set_ylabel(ylabel, fontsize=17, color=INK)
        axis.grid(axis="y", color=GRID, linewidth=1.0, zorder=0)
        axis.tick_params(axis="both", labelsize=13.5, colors=INK, length=5, width=1.0)
        for spine in axis.spines.values():
            spine.set_color("#c8d3e2")
            spine.set_linewidth(1.0)

        pp_points = series(points, region, "pp", "inclusive")
        axis.errorbar(
            [(point.et_low_gev + point.et_high_gev) / 2.0 for point in pp_points],
            [point.leakage_ratio for point in pp_points],
            xerr=[(point.et_high_gev - point.et_low_gev) / 2.0 for point in pp_points],
            yerr=[point.leakage_error for point in pp_points],
            fmt="o",
            color=PP_COLOR,
            markeredgecolor="white",
            markeredgewidth=0.8,
            markersize=7.5,
            elinewidth=1.2,
            capsize=2.5,
            label="p+p",
            zorder=5,
        )
        for cent_key, cent_label, color, marker in CENTRALITIES:
            del cent_key
            cent_points = series(points, region, "auau", cent_label)
            axis.errorbar(
                [(point.et_low_gev + point.et_high_gev) / 2.0 for point in cent_points],
                [point.leakage_ratio for point in cent_points],
                xerr=[(point.et_high_gev - point.et_low_gev) / 2.0 for point in cent_points],
                yerr=[point.leakage_error for point in cent_points],
                fmt=marker,
                color=color,
                markeredgecolor="white",
                markeredgewidth=0.8,
                markersize=7.2,
                elinewidth=1.15,
                capsize=2.3,
                label=f"Au+Au {cent_label}",
                zorder=4,
            )

        maximum = max(point.leakage_ratio + point.leakage_error for point in points if point.region == region)
        axis.set_ylim(0.0, maximum * 1.28)
        axis.set_xlim(15.0, 35.0)

    axes[0].legend(
        loc="upper left",
        bbox_to_anchor=(0.18, 1.02),
        ncol=4,
        frameon=False,
        fontsize=14.2,
        handletextpad=0.4,
        columnspacing=1.3,
    )
    axes[-1].set_xlabel(r"reconstructed photon-cluster $E_T$ [GeV]", fontsize=17, color=INK)
    axes[-1].set_xticks([15, 17, 19, 21, 23, 26, 30, 35])
    fig.savefig(output, dpi=180, facecolor="white")
    plt.close(fig)


def write_outputs(points: list[Point], pp_root: Path, auau_root: Path, out_dir: Path) -> tuple[Path, Path, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    png = out_dir / "current_pp_vs_the88_auau_signal_leakage_bcd_15to35.png"
    csv_path = out_dir / "current_pp_vs_the88_auau_signal_leakage_bcd_15to35.csv"
    manifest_path = out_dir / "current_pp_vs_the88_auau_signal_leakage_bcd_15to35_manifest.json"

    rows = [asdict(point) for point in points]
    with csv_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    draw(points, png)
    manifest = {
        "schema_version": 1,
        "plot": str(png),
        "csv": str(csv_path),
        "quantity": "truth-signal sideband leakage N_X^sig / N_A^sig for X=B,C,D",
        "selection": {
            "pp": "canonical current PPG12 pp photon+jet output; R=0.4 sliding isolation",
            "auau": "THE-88 embedded Photon12+20 weighted signal stitch; R=0.4 sliding isolation",
            "centrality_percent": [[0, 20], [20, 50], [50, 80]],
            "cluster_et_range_gev": [15, 35],
        },
        "source_native_binning": {
            "pp": PP_BINS,
            "auau": AUAU_BINS,
            "policy": "No fractional splitting of stored counters. Horizontal error bars show each source interval.",
        },
        "inputs": {
            "pp_root": str(pp_root),
            "pp_sha256": sha256(pp_root),
            "auau_root": str(auau_root),
            "auau_sha256": sha256(auau_root),
        },
        "production_artifact_exceptions": [
            "The validated 2026-07-15 THE-97 accepted-triple photon+jet ROOT is newer than the registered pp current pointer and is used because the user requested the most up-to-date pp output.",
            "The AuAu input is the validated historical THE-88 signal product matched to the corrected THE-88 data configuration; no generic AuAu current pointer exists.",
        ],
        "uncertainty": "Independent numerator/denominator propagation using ROOT bin errors (Sumw2 where stored).",
        "point_count": len(points),
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    return png, csv_path, manifest_path


def main() -> int:
    args = parse_args()
    for path in (args.pp_root, args.auau_root):
        if not path.is_file():
            raise FileNotFoundError(path)
    points = collect(args.pp_root, args.auau_root)
    png, csv_path, manifest_path = write_outputs(points, args.pp_root, args.auau_root, args.out_dir)
    print(f"wrote {png}")
    print(f"wrote {csv_path}")
    print(f"wrote {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
