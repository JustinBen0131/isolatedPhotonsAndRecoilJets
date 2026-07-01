#!/usr/bin/env python3
"""Audit the two pp stitching contracts against PPG12 Fig. 5/6 sources.

This intentionally separates two concepts that should not be mixed:

1. clean_density
   The in-situ RecoilJets QA contract used for the smooth slide-34/35 style
   closure: raw counts are normalized by processed events and bin width, so
   the result is a physical-looking density.

2. ppg12_verbatim
   The source-histogram convention used by the current PPG12 IAN Fig. 5/6
   extraction. Photon+jet Fig. 5 uses per-bin weighted contents
   raw * xsec / events_processed (no bin-width division). Inclusive+jet Fig. 6
   uses PPG12 efficiency-tool weighted counts raw * (xsec / jet50), then needs
   an exposure normalization if our local processed event count is not the
   same as Shuhang's source ROOTs.

The output is a local diagnostic artifact, not a production merge.
"""

from __future__ import annotations

# Keep purpose-folder helpers runnable when invoked directly.
import sys as _codex_sys
from pathlib import Path as _CodexPath

_CODEX_THIS_FILE = _CodexPath(__file__).resolve()
_CODEX_SCRIPTS_DIR = next(
    (p for p in _CODEX_THIS_FILE.parents if p.name == "scripts"),
    _CODEX_THIS_FILE.parent,
)
_CODEX_SCRIPTS_DIR_STR = str(_CODEX_SCRIPTS_DIR)
if _CODEX_SCRIPTS_DIR_STR not in _codex_sys.path:
    _codex_sys.path.append(_CODEX_SCRIPTS_DIR_STR)
del _CODEX_THIS_FILE, _CODEX_SCRIPTS_DIR, _CODEX_SCRIPTS_DIR_STR

import argparse
import csv
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from statistics import median
from typing import Iterable

try:
    import matplotlib.pyplot as plt
    import numpy as np
except ModuleNotFoundError as exc:
    if "-h" in sys.argv or "--help" in sys.argv:
        plt = None
        np = None
    else:
        raise SystemExit(
            "This diagnostic needs the ThesisAnalysis Python environment. Run:\n"
            "/Users/patsfan753/Desktop/analysis/env/bin/python "
            "scripts/diagnostics/pp_shuhang/audit_pp_stitching_dual_contracts.py"
        ) from exc


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
DEFAULT_REF_CSV = (
    REPO
    / "dataOutput/ppPhotonMLPipeline/ppg12_basev3E_currentIAN_fullsim_20260521_1811/"
    / "validation/currentIAN_stitching/ppg12_currentian_rootfit_stitch_points.csv"
)
DEFAULT_PHOTON_CSV = (
    REPO
    / "dataOutput/ppPhotonMLPipeline/ppg12_basev3E_currentIAN_exactStitchPhoton0p5_20260526_1307/"
    / "validation/insitu_stitching/pp_currentian_photon0p5_exactstitch_contract_points.csv"
)
DEFAULT_JET_CSV = (
    REPO
    / "dataOutput/ppPhotonMLPipeline/ppg12_basev3E_currentIAN_exactStitch_20260523_0915/"
    / "validation/insitu_stitching/pp_currentian_exactstitch_contract_points.csv"
)
DEFAULT_OUT_DIR = (
    REPO
    / "dataOutput/the85_auau_xjgamma_unfolding_push/"
    / "ppg12_stitching_contract_debug_20260630"
)

PHOTON_XSECS_PB = {
    "run28_photonjet5": 146359.3,
    "run28_photonjet10": 6944.675,
    "run28_photonjet20": 130.4461,
}
JET_XSECS_CURRENT_PB = {
    "run28_jet8": 1.3013e7,
    "run28_jet12": 1.4903e6,
    "run28_jet20": 6.2623e4,
    "run28_jet30": 2.5298e3,
    "run28_jet40": 1.3553e2,
}
JET_XSECS_PPG12_LEGACY_PB = {
    "run28_jet8": 1.15e7,
    "run28_jet12": 1.4903e6,
    "run28_jet20": 6.2623e4,
    "run28_jet30": 2.5298e3,
    "run28_jet40": 1.3553e2,
}
JET50_XSEC_PB = 7.3113

SAMPLE_TO_REF = {
    "run28_photonjet5": "photon5",
    "run28_photonjet10": "photon10",
    "run28_photonjet20": "photon20",
    "run28_jet8": "jet8",
    "run28_jet12": "jet12",
    "run28_jet20": "jet20",
    "run28_jet30": "jet30",
    "run28_jet40": "jet40",
}
REF_TO_CURRENT = {v: k for k, v in SAMPLE_TO_REF.items()}

COLORS = {
    "photon5": "#d62aa0",
    "photon10": "#2ca02c",
    "photon20": "#1296f3",
    "jet8": "#d62aa0",
    "jet12": "#2ca02c",
    "jet20": "#1296f3",
    "jet30": "#ff7f00",
    "jet40": "#d62aa0",
}


@dataclass(frozen=True)
class Point:
    group: str
    sample: str
    x: float
    y: float
    ey: float


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        raise FileNotFoundError(path)
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def is_window_row(row: dict[str, str]) -> bool:
    if "ppg12_bin_center_window" in row:
        return row["ppg12_bin_center_window"] == "1"
    if "used_in_stitch" in row:
        return row["used_in_stitch"] == "1"
    return False


def load_reference_points(ref_csv: Path, group: str) -> dict[tuple[str, float], Point]:
    points: dict[tuple[str, float], Point] = {}
    for row in read_csv(ref_csv):
        if row["group"] != group or row["sample"] == "stitched":
            continue
        if not is_window_row(row):
            continue
        y = float(row["value"])
        if y <= 0:
            continue
        sample = row["sample"]
        x = round(float(row["bin_center"]), 6)
        points[(sample, x)] = Point(group, sample, x, y, float(row["error"]))
    if not points:
        raise RuntimeError(f"No PPG12 source points loaded for group={group}")
    return points


def load_reference_fit(ref_csv: Path, group: str) -> tuple[str, str]:
    for row in read_csv(ref_csv):
        if row["group"] == group and row["sample"] == "stitched" and row["root_fit_params"]:
            return row["root_fit_params"], row["ratio_label"]
    raise RuntimeError(f"No PPG12 fit contract loaded for group={group}")


def eval_root_fit(x: np.ndarray, params: str) -> np.ndarray:
    p = np.array([float(v) for v in params.split(";")], dtype=float)
    return p[0] * np.power(p[1] / x, p[2] + p[3] * np.log(x / p[1]) + p[4] * x)


def current_photon_points(rows: Iterable[dict[str, str]], mode: str) -> dict[tuple[str, float], Point]:
    points: dict[tuple[str, float], Point] = {}
    for row in rows:
        if row["group"] != "photon" or not is_window_row(row):
            continue
        sample = row["sample"]
        ref_sample = SAMPLE_TO_REF.get(sample)
        if ref_sample is None:
            continue
        x = round(float(row["bin_center"]), 6)
        raw = float(row["raw_all_events"])
        xsec = float(row["xsec_pb"])
        events = float(row["events_processed_metadata"])
        width = float(row["bin_high"]) - float(row["bin_low"])
        if raw <= 0 or events <= 0 or width <= 0:
            continue
        if mode == "clean_density":
            y = float(row["density_pb_per_gev"])
            ey = float(row["density_err_pb_per_gev"])
        elif mode == "ppg12_verbatim":
            y = raw * xsec / events
            ey = math.sqrt(raw) * xsec / events
        else:
            raise ValueError(mode)
        points[(ref_sample, x)] = Point("photon", ref_sample, x, y, ey)
    return points


def current_jet_points(
    rows: Iterable[dict[str, str]],
    mode: str,
    exposure_scales: dict[str, float] | None = None,
) -> dict[tuple[str, float], Point]:
    points: dict[tuple[str, float], Point] = {}
    for row in rows:
        if row["group"] != "jet" or not is_window_row(row):
            continue
        sample = row["sample"]
        ref_sample = SAMPLE_TO_REF.get(sample)
        if ref_sample is None:
            continue
        x = round(float(row["bin_center"]), 6)
        raw = float(row["raw_all_events"])
        width = float(row["bin_high"]) - float(row["bin_low"])
        if raw <= 0 or width <= 0:
            continue
        if mode == "clean_density":
            y = float(row["density_pb_per_gev"])
            ey = float(row["density_err_pb_per_gev"])
        elif mode in {"ppg12_weighted_counts", "ppg12_weighted_counts_scaled"}:
            xsec = JET_XSECS_PPG12_LEGACY_PB[sample]
            scale = 1.0 if exposure_scales is None else exposure_scales.get(ref_sample, 1.0)
            y = raw * xsec / JET50_XSEC_PB * scale
            ey = math.sqrt(raw) * xsec / JET50_XSEC_PB * scale
        else:
            raise ValueError(mode)
        points[(ref_sample, x)] = Point("jet", ref_sample, x, y, ey)
    return points


def ratios(current: dict[tuple[str, float], Point], reference: dict[tuple[str, float], Point]) -> dict[str, list[float]]:
    out: dict[str, list[float]] = {}
    for key, cur in current.items():
        ref = reference.get(key)
        if ref is None or ref.y <= 0 or cur.y <= 0:
            continue
        out.setdefault(cur.sample, []).append(cur.y / ref.y)
    return out


def summarize_ratio_values(values: list[float]) -> dict[str, float | int]:
    finite = [v for v in values if math.isfinite(v)]
    if not finite:
        return {"n": 0}
    return {
        "n": len(finite),
        "mean": float(np.mean(finite)),
        "median": float(np.median(finite)),
        "rms_to_one": float(math.sqrt(np.mean([(v - 1.0) ** 2 for v in finite]))),
        "min": float(min(finite)),
        "max": float(max(finite)),
    }


def summarize_by_sample(
    name: str,
    current: dict[tuple[str, float], Point],
    reference: dict[tuple[str, float], Point],
) -> dict[str, object]:
    sample_ratios = ratios(current, reference)
    all_ratios = [v for values in sample_ratios.values() for v in values]
    return {
        "comparison": name,
        "all": summarize_ratio_values(all_ratios),
        "by_sample": {sample: summarize_ratio_values(values) for sample, values in sorted(sample_ratios.items())},
    }


def exposure_scales_from_reference(
    current: dict[tuple[str, float], Point],
    reference: dict[tuple[str, float], Point],
    samples: Iterable[str] | None = None,
) -> dict[str, float]:
    wanted = set(samples) if samples is not None else None
    grouped: dict[str, list[float]] = {}
    for key, cur in current.items():
        ref = reference.get(key)
        if ref is None or cur.y <= 0 or ref.y <= 0:
            continue
        if wanted is not None and cur.sample not in wanted:
            continue
        grouped.setdefault(cur.sample, []).append(ref.y / cur.y)
    return {sample: float(median(values)) for sample, values in grouped.items() if values}


def one_global_scale_from_samples(
    current: dict[tuple[str, float], Point],
    reference: dict[tuple[str, float], Point],
    samples: Iterable[str],
) -> float:
    wanted = set(samples)
    vals: list[float] = []
    for key, cur in current.items():
        ref = reference.get(key)
        if ref is None or cur.y <= 0 or ref.y <= 0 or cur.sample not in wanted:
            continue
        vals.append(ref.y / cur.y)
    if not vals:
        raise RuntimeError("Cannot compute global scale")
    return float(median(vals))


def apply_sample_scales(
    current: dict[tuple[str, float], Point],
    scales: dict[str, float],
) -> dict[tuple[str, float], Point]:
    out: dict[tuple[str, float], Point] = {}
    for key, p in current.items():
        scale = scales.get(p.sample, 1.0)
        out[key] = Point(p.group, p.sample, p.x, p.y * scale, p.ey * scale)
    return out


def plot_overlay(
    *,
    group: str,
    reference: dict[tuple[str, float], Point],
    current: dict[tuple[str, float], Point],
    fit_params: str,
    ratio_label: str,
    out_png: Path,
    xlim: tuple[float, float],
    ylim: tuple[float, float],
    xlabel: str,
    ylabel: str,
    note: str,
) -> None:
    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "mathtext.fontset": "dejavusans",
            "axes.linewidth": 1.15,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )
    fig = plt.figure(figsize=(7.0, 8.8), dpi=220)
    gs = fig.add_gridspec(2, 1, height_ratios=(3.2, 1.0), hspace=0.035, left=0.16, right=0.96, top=0.965, bottom=0.105)
    ax = fig.add_subplot(gs[0])
    rax = fig.add_subplot(gs[1], sharex=ax)
    ax.set_yscale("log")
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_ylabel(ylabel, fontsize=16)
    rax.set_xlabel(xlabel, fontsize=16)
    rax.set_ylabel("current / PPG12", fontsize=13)
    ax.tick_params(which="both", labelsize=13.5, labelbottom=False, length=6)
    rax.tick_params(which="both", labelsize=13.5, length=6)

    grid = np.linspace(max(xlim[0], 1e-3), xlim[1], 800)
    fit = eval_root_fit(grid, fit_params)
    ax.plot(grid, fit, color="0.25", lw=1.2, zorder=1, label="PPG12 fit")

    samples = sorted({k[0] for k in reference.keys()}, key=lambda s: (s.rstrip("0123456789"), float("".join(c for c in s if c.isdigit()) or 0)))
    for sample in samples:
        color = COLORS.get(sample, "black")
        ref_pts = [p for (s, _), p in reference.items() if s == sample and xlim[0] <= p.x <= xlim[1]]
        cur_pts = [current[(s, x)] for (s, x), p in reference.items() if s == sample and (s, x) in current and xlim[0] <= p.x <= xlim[1]]
        if ref_pts:
            rx = np.array([p.x for p in ref_pts])
            ry = np.array([p.y for p in ref_pts])
            rey = np.array([p.ey for p in ref_pts])
            ax.errorbar(rx, ry, yerr=rey, fmt="s", ms=4.2, mfc="white", mec=color, mew=1.0, ecolor=color, elinewidth=0.65, linestyle="none", zorder=3)
        if cur_pts:
            cx = np.array([p.x for p in cur_pts])
            cy = np.array([p.y for p in cur_pts])
            cey = np.array([p.ey for p in cur_pts])
            ratios_to_ref = []
            ratio_err = []
            for p in cur_pts:
                ref = reference[(p.sample, p.x)]
                ratios_to_ref.append(p.y / ref.y)
                ratio_err.append(p.ey / ref.y)
            ax.errorbar(cx, cy, yerr=cey, fmt="o", ms=4.0, mfc=color, mec=color, mew=0.75, ecolor=color, elinewidth=0.65, linestyle="none", zorder=4)
            rax.errorbar(cx, ratios_to_ref, yerr=ratio_err, fmt="o", ms=3.5, mfc=color, mec=color, mew=0.65, ecolor=color, elinewidth=0.6, linestyle="none", zorder=4)

    from matplotlib.lines import Line2D

    ax.legend(
        handles=[
            Line2D([0], [0], marker="s", color="black", mfc="white", mec="black", lw=0, label="PPG12 SDCC"),
            Line2D([0], [0], marker="o", color="black", mfc="black", mec="black", lw=0, label="RecoilJets contract"),
            Line2D([0], [0], color="0.25", lw=1.2, label="PPG12 fit"),
        ],
        loc="lower left",
        frameon=False,
        fontsize=11.5,
        handletextpad=0.45,
        borderpad=0.2,
    )
    ax.text(0.97, 0.94, r"$\it{\bf{sPHENIX}}$ Internal", transform=ax.transAxes, ha="right", va="top", fontsize=15)
    ax.text(0.97, 0.875, note, transform=ax.transAxes, ha="right", va="top", fontsize=11.2)
    rax.axhline(1.0, color="0.45", lw=1.0, ls="--")
    rax.set_ylim(0.82, 1.18)
    if group == "jet":
        rax.set_ylim(0.84, 1.28)
    fig.savefig(out_png)
    plt.close(fig)


def write_metrics_csv(path: Path, summaries: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    for summary in summaries:
        comparison = str(summary["comparison"])
        for sample, metrics in summary["by_sample"].items():
            row = {"comparison": comparison, "sample": sample}
            row.update(metrics)
            rows.append(row)
        row = {"comparison": comparison, "sample": "all"}
        row.update(summary["all"])
        rows.append(row)
    fields = ["comparison", "sample", "n", "mean", "median", "rms_to_one", "min", "max"]
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-csv", type=Path, default=DEFAULT_REF_CSV)
    parser.add_argument("--photon-csv", type=Path, default=DEFAULT_PHOTON_CSV)
    parser.add_argument("--jet-csv", type=Path, default=DEFAULT_JET_CSV)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    args = parser.parse_args()

    photon_rows = read_csv(args.photon_csv)
    jet_rows = read_csv(args.jet_csv)
    photon_ref = load_reference_points(args.reference_csv, "photon")
    jet_ref = load_reference_points(args.reference_csv, "jet")
    photon_fit, photon_ratio_label = load_reference_fit(args.reference_csv, "photon")
    jet_fit, jet_ratio_label = load_reference_fit(args.reference_csv, "jet")

    photon_density = current_photon_points(photon_rows, "clean_density")
    photon_verbatim = current_photon_points(photon_rows, "ppg12_verbatim")
    jet_density = current_jet_points(jet_rows, "clean_density")
    jet_weighted = current_jet_points(jet_rows, "ppg12_weighted_counts")

    jet_global_nonjet_scale = one_global_scale_from_samples(
        jet_weighted,
        jet_ref,
        samples=["jet12", "jet20", "jet30", "jet40"],
    )
    jet_global_scaled = apply_sample_scales(
        jet_weighted,
        {sample: jet_global_nonjet_scale for sample in ["jet8", "jet12", "jet20", "jet30", "jet40"]},
    )
    jet_sample_exposure_scales = exposure_scales_from_reference(jet_weighted, jet_ref)
    jet_sample_scaled = apply_sample_scales(jet_weighted, jet_sample_exposure_scales)

    summaries = [
        summarize_by_sample("photon_clean_density_vs_ppg12_source", photon_density, photon_ref),
        summarize_by_sample("photon_ppg12_verbatim_vs_ppg12_source", photon_verbatim, photon_ref),
        summarize_by_sample("jet_clean_density_vs_ppg12_source_not_a_valid_contract", jet_density, jet_ref),
        summarize_by_sample("jet_ppg12_legacy_weighted_counts_unscaled_vs_ppg12_source", jet_weighted, jet_ref),
        summarize_by_sample("jet_ppg12_legacy_weighted_counts_global_nonjet_exposure_scaled", jet_global_scaled, jet_ref),
        summarize_by_sample("jet_ppg12_legacy_weighted_counts_per_sample_exposure_scaled", jet_sample_scaled, jet_ref),
    ]

    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    plot_overlay(
        group="photon",
        reference=photon_ref,
        current=photon_verbatim,
        fit_params=photon_fit,
        ratio_label=photon_ratio_label,
        out_png=out_dir / "photon_ppg12_verbatim_current_vs_sdcc_ratio.png",
        xlim=(10.0, 40.0),
        ylim=(2.0e-2, 2.0e3),
        xlabel=r"Leading $p_T^\gamma$ [GeV]",
        ylabel=r"weighted entries / bin",
        note="photon+jet Fig. 5 contract",
    )
    plot_overlay(
        group="jet",
        reference=jet_ref,
        current=jet_global_scaled,
        fit_params=jet_fit,
        ratio_label=jet_ratio_label,
        out_png=out_dir / "inclusivejet_ppg12_legacy_global_nonjet_scaled_vs_sdcc_ratio.png",
        xlim=(9.0, 50.0),
        ylim=(1.0e4, 1.0e12),
        xlabel=r"Leading truth-jet $p_T$ [GeV]",
        ylabel="weighted entries",
        note="inclusive+jet Fig. 6, one non-jet global exposure scale",
    )
    plot_overlay(
        group="jet",
        reference=jet_ref,
        current=jet_sample_scaled,
        fit_params=jet_fit,
        ratio_label=jet_ratio_label,
        out_png=out_dir / "inclusivejet_ppg12_legacy_per_sample_scaled_vs_sdcc_ratio.png",
        xlim=(9.0, 50.0),
        ylim=(1.0e4, 1.0e12),
        xlabel=r"Leading truth-jet $p_T$ [GeV]",
        ylabel="weighted entries",
        note="inclusive+jet Fig. 6, per-sample exposure matched",
    )

    metrics_csv = out_dir / "ppg12_stitching_dual_contract_metrics.csv"
    write_metrics_csv(metrics_csv, summaries)
    manifest = {
        "inputs": {
            "ppg12_sdcc_reference_csv": str(args.reference_csv),
            "current_photon_insitu_csv": str(args.photon_csv),
            "current_jet_insitu_csv": str(args.jet_csv),
        },
        "contracts": {
            "clean_density": "RecoilJets smooth QA: raw * xsec / events_processed / bin_width.",
            "photon_ppg12_verbatim": "raw * xsec / events_processed; no bin-width division. This fixes the factor-of-two photon discrepancy from 0.5 GeV bins.",
            "jet_ppg12_verbatim": "raw * (legacy PPG12 xsec / jet50); compare shapes after exposure normalization because PPG12 source ROOTs are weighted-count histograms, not physical cross-section densities.",
        },
        "ppg12_legacy_constants": {
            "jet50_xsec_pb": JET50_XSEC_PB,
            "jet_xsecs_pb": JET_XSECS_PPG12_LEGACY_PB,
            "current_clean_jet_xsecs_pb": JET_XSECS_CURRENT_PB,
            "photon_xsecs_pb": PHOTON_XSECS_PB,
        },
        "derived_scales": {
            "jet_global_nonjet_exposure_scale": jet_global_nonjet_scale,
            "jet_per_sample_exposure_scales": jet_sample_exposure_scales,
            "interpretation": "The non-jet samples share one exposure scale; jet8 has a separate residual. This is a source-exposure/legacy-contract issue, not a reason to use clean density as PPG12 Fig. 6.",
        },
        "outputs": {
            "metrics_csv": str(metrics_csv),
            "photon_png": str(out_dir / "photon_ppg12_verbatim_current_vs_sdcc_ratio.png"),
            "inclusivejet_global_nonjet_scaled_png": str(out_dir / "inclusivejet_ppg12_legacy_global_nonjet_scaled_vs_sdcc_ratio.png"),
            "inclusivejet_per_sample_scaled_png": str(out_dir / "inclusivejet_ppg12_legacy_per_sample_scaled_vs_sdcc_ratio.png"),
        },
        "summaries": summaries,
    }
    manifest_path = out_dir / "ppg12_stitching_dual_contract_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"manifest": str(manifest_path), "metrics_csv": str(metrics_csv)}, indent=2))


if __name__ == "__main__":
    main()
