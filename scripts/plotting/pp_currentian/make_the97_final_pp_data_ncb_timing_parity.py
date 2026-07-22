#!/usr/bin/env python3
"""Render full-stat PPG12/RecoilJets NCB timing comparisons for the IAN.

The PPG12 references are the locally recovered SDCC ROOT objects used for
IAN Figures 14 and 15.  The current side is resolved from the registered
``pp_data_merged/current.json`` pointer.  Historical ROOT object names retain
``npb``; reader-facing labels use NCB (non-collision background).
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import uproot


REPO = Path(__file__).resolve().parents[3]
CURRENT_POINTER = REPO / "dataOutput/current_recoiljets_artifacts/current/pp_data_merged/current.json"
REFERENCE_TIMING = (
    REPO
    / "dataOutput/ppg12PhotonYield/ppg12_photon_yield_v1_data_20260620"
    / "reference_roots/fig14_15_truth_photon/cluster_time_analysis_data.root"
)
REFERENCE_SCORE = (
    REPO
    / "dataOutput/ppg12PhotonYield/ppg12_photon_yield_v1_data_20260620"
    / "reference_roots/fig14_15_truth_photon/data_histoshower_shape_.root"
)
DEFAULT_OUTDIR = (
    REPO
    / "dataOutput/ppg12Parity/the97_ppg12_final_accepted_triple_full_20260714_1550"
    / "final_pp_data_canonical_20260717/ncb_timing_parity"
)

CURRENT_DIR = "PPG12_scaledtrigger30"
TIMING_KEY = "h_all_delta_t_mbd_vs_eta"
REFERENCE_NCB_TIMING_KEY = "h_npb_delta_t_mbd_vs_eta"
SCORE_TIME_KEY = "h_npb_score_vs_time_eta0_pt0"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def configure_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "dejavuserif",
            "font.size": 12,
            "axes.linewidth": 1.15,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "legend.frameon": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )


def load_current() -> tuple[Path, dict[str, Any]]:
    pointer = json.loads(CURRENT_POINTER.read_text())
    roots = pointer.get("root_paths") or []
    if len(roots) != 1:
        raise RuntimeError(f"expected one registered pp-data ROOT, found {roots}")
    root = Path(roots[0])
    if not root.exists():
        raise FileNotFoundError(root)
    expected = "5854abfa21629482ca3aef87c421b790945c093ca97dcee6093a7d6af53f9a6f"
    actual = sha256(root)
    if actual != expected:
        raise RuntimeError(f"registered pp-data SHA mismatch: {actual} != {expected}")
    return root, pointer


def read_th2(path: Path, key: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    with uproot.open(path) as handle:
        hist = handle[key]
        values, xedges, yedges = hist.to_numpy(flow=False)
        variances = hist.variances(flow=False)
    if variances is None:
        variances = np.clip(values, 0.0, None)
    return (
        np.asarray(values, dtype=float),
        np.asarray(variances, dtype=float),
        np.asarray(xedges, dtype=float),
        np.asarray(yedges, dtype=float),
    )


def profile_y(values: np.ndarray, variances: np.ndarray, yedges: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    centers = 0.5 * (yedges[:-1] + yedges[1:])
    counts = np.sum(values, axis=1)
    weighted = np.sum(values * centers[None, :], axis=1)
    mean = np.divide(weighted, counts, out=np.full_like(counts, np.nan), where=counts > 0)
    second = np.divide(
        np.sum(values * centers[None, :] ** 2, axis=1),
        counts,
        out=np.full_like(counts, np.nan),
        where=counts > 0,
    )
    spread = np.sqrt(np.clip(second - mean**2, 0.0, None))
    # ROOT-profile-like uncertainty on the mean.  Counts are unweighted here;
    # Sumw2 is retained in the source manifest for provenance.
    error = np.divide(spread, np.sqrt(counts), out=np.full_like(spread, np.nan), where=counts > 0)
    _ = variances
    return mean, error


def render_timing_profile(
    out: Path,
    current_root: Path,
    current: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray],
    reference_all: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray],
    reference_ncb: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray],
) -> list[dict[str, float]]:
    cvals, cvar, xedges, yedges = current
    rvals, rvar, rxedges, ryedges = reference_all
    nvals, nvar, nxedges, nyedges = reference_ncb
    for candidate in (rxedges, nxedges):
        if not np.allclose(candidate, xedges):
            raise RuntimeError("timing-versus-eta x binning mismatch")
    for candidate in (ryedges, nyedges):
        if not np.allclose(candidate, yedges):
            raise RuntimeError("timing-versus-eta y binning mismatch")
    cmean, cerr = profile_y(cvals, cvar, yedges)
    rmean, rerr = profile_y(rvals, rvar, yedges)
    nmean, nerr = profile_y(nvals, nvar, yedges)
    x = 0.5 * (xedges[:-1] + xedges[1:])
    residual = cmean - rmean
    residual_err = np.sqrt(cerr**2 + rerr**2)

    fig, (ax, rax) = plt.subplots(
        2,
        1,
        figsize=(8.4, 7.2),
        sharex=True,
        gridspec_kw={"height_ratios": [3.0, 1.05], "hspace": 0.03},
    )
    ax.errorbar(x, rmean, yerr=rerr, fmt="o", ms=4.0, mfc="white", mec="black", color="black", label="PPG12 all clusters")
    ax.errorbar(x, cmean, yerr=cerr, fmt="o", ms=3.5, color="#d62728", label="This analysis")
    ax.errorbar(x, nmean, yerr=nerr, fmt="s", ms=3.0, mfc="white", mec="#1f77b4", color="#1f77b4", label="PPG12 NCB-tagged reference")
    ax.axhline(0.0, color="0.6", lw=0.9, ls="--")
    ax.set_ylabel(r"Mean cluster$-$MBD time [ns]")
    ax.set_ylim(-5.2, 2.2)
    ax.grid(axis="y", color="0.9", lw=0.7)
    ax.legend(loc="lower left", fontsize=10.5, ncol=1)
    ax.text(0.97, 0.95, r"$\bf{\it{sPHENIX}}$ Internal", ha="right", va="top", transform=ax.transAxes, fontsize=14)
    ax.text(0.97, 0.865, r"$p+p\ \sqrt{s}=200$ GeV", ha="right", va="top", transform=ax.transAxes)
    ax.text(0.97, 0.80, "Photon-4-GeV trigger", ha="right", va="top", transform=ax.transAxes, fontsize=10.5)
    rax.axhline(0.0, color="0.45", lw=1.0, ls="--")
    rax.errorbar(x, residual, yerr=residual_err, fmt="o", ms=3.5, color="#d62728")
    rax.set_ylabel("This analysis $-$\nPPG12 [ns]", fontsize=10.5)
    rax.set_xlabel(r"Cluster $\eta$")
    rax.set_xlim(-1.0, 1.0)
    finite_residual = residual[np.isfinite(residual)]
    if finite_residual.size:
        span = max(0.35, 0.12 * float(np.ptp(finite_residual)))
        rax.set_ylim(float(np.min(finite_residual)) - span, float(np.max(finite_residual)) + span)
    else:
        rax.set_ylim(-1.0, 1.0)
    rax.grid(axis="y", color="0.9", lw=0.7)
    fig.subplots_adjust(left=0.14, right=0.97, bottom=0.11, top=0.97)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=220)
    plt.close(fig)

    return [
        {
            "eta": float(xx),
            "ppg12_all_mean_ns": float(rm),
            "ppg12_all_error_ns": float(re),
            "current_all_mean_ns": float(cm),
            "current_all_error_ns": float(ce),
            "current_minus_ppg12_ns": float(dr),
            "current_minus_ppg12_error_ns": float(de),
            "ppg12_ncb_mean_ns": float(nm),
            "ppg12_ncb_error_ns": float(ne),
        }
        for xx, rm, re, cm, ce, dr, de, nm, ne in zip(
            x, rmean, rerr, cmean, cerr, residual, residual_err, nmean, nerr
        )
    ]


def normalize_score_rows(values: np.ndarray) -> np.ndarray:
    # ROOT axes are time (x) and NCB score (y); normalize each score bin over time.
    sums = np.sum(values, axis=0, keepdims=True)
    return np.divide(values, sums, out=np.zeros_like(values), where=sums > 0)


def render_score_timing(
    out: Path,
    current: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray],
    reference: tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray],
) -> list[dict[str, float]]:
    cvals, _, xedges, yedges = current
    rvals, _, rxedges, ryedges = reference
    if not np.allclose(xedges, rxedges) or not np.allclose(yedges, ryedges):
        raise RuntimeError("NCB-score timing binning mismatch")
    cnorm = normalize_score_rows(cvals)
    rnorm = normalize_score_rows(rvals)
    diff = cnorm - rnorm
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    common_max = float(max(np.max(cnorm), np.max(rnorm)))
    diff_max = float(np.max(np.abs(diff)))

    fig, axes = plt.subplots(1, 3, figsize=(12.0, 4.15), constrained_layout=True)
    for ax, array, title in zip(
        axes[:2],
        (rnorm, cnorm),
        ("PPG12 SDCC reference", "This analysis"),
    ):
        image = ax.imshow(array.T, origin="lower", aspect="auto", extent=extent, cmap="viridis", vmin=0.0, vmax=common_max)
        ax.set_title(title, fontsize=12.5, fontweight="bold")
        ax.set_xlabel(r"Cluster$-$MBD time [ns]")
        ax.set_ylabel("NCB score")
        fig.colorbar(image, ax=ax, pad=0.015, fraction=0.047, label="Conditional fraction")
    image = axes[2].imshow(diff.T, origin="lower", aspect="auto", extent=extent, cmap="RdBu_r", vmin=-diff_max, vmax=diff_max)
    axes[2].set_title("This analysis $-$ PPG12", fontsize=12.5, fontweight="bold")
    axes[2].set_xlabel(r"Cluster$-$MBD time [ns]")
    axes[2].set_ylabel("NCB score")
    fig.colorbar(image, ax=axes[2], pad=0.015, fraction=0.047, label="Conditional-fraction difference")
    fig.suptitle(r"$\bf{\it{sPHENIX}}$ Internal   $p+p\ \sqrt{s}=200$ GeV   $10<E_T^{\mathrm{cluster}}<14$ GeV, $|\eta|<0.7$", fontsize=13.5)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=220)
    plt.close(fig)

    time_centers = 0.5 * (xedges[:-1] + xedges[1:])
    score_centers = 0.5 * (yedges[:-1] + yedges[1:])
    rows: list[dict[str, float]] = []
    for ix, time in enumerate(time_centers):
        for iy, score in enumerate(score_centers):
            rows.append(
                {
                    "time_ns": float(time),
                    "ncb_score": float(score),
                    "ppg12_conditional_fraction": float(rnorm[ix, iy]),
                    "current_conditional_fraction": float(cnorm[ix, iy]),
                    "current_minus_ppg12": float(diff[ix, iy]),
                }
            )
    return rows


def write_csv(path: Path, rows: list[dict[str, float]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = parser.parse_args()
    configure_style()
    current_root, pointer = load_current()
    args.outdir.mkdir(parents=True, exist_ok=True)

    current_timing = read_th2(current_root, f"{CURRENT_DIR}/{TIMING_KEY}")
    reference_all = read_th2(REFERENCE_TIMING, TIMING_KEY)
    reference_ncb = read_th2(REFERENCE_TIMING, REFERENCE_NCB_TIMING_KEY)
    timing_png = args.outdir / "the97_fullstat_timing_mean_vs_eta_ppg12_vs_current.png"
    timing_csv = args.outdir / "the97_fullstat_timing_mean_vs_eta_ppg12_vs_current.csv"
    write_csv(
        timing_csv,
        render_timing_profile(timing_png, current_root, current_timing, reference_all, reference_ncb),
    )

    current_score = read_th2(current_root, f"{CURRENT_DIR}/{SCORE_TIME_KEY}")
    reference_score = read_th2(REFERENCE_SCORE, SCORE_TIME_KEY)
    score_png = args.outdir / "the97_fullstat_ncb_score_vs_timing_ppg12_vs_current.png"
    score_csv = args.outdir / "the97_fullstat_ncb_score_vs_timing_ppg12_vs_current.csv"
    write_csv(score_csv, render_score_timing(score_png, current_score, reference_score))

    manifest = {
        "schema": "THE97_FINAL_PP_DATA_NCB_TIMING_PARITY_V1",
        "status": "FULL_STAT_CURRENT",
        "reader_facing_term": "NCB (non-collision background)",
        "historical_object_name_note": "ROOT keys retain npb for compatibility with PPG12 source naming.",
        "current": {
            "pointer": str(CURRENT_POINTER),
            "pointer_sha256": sha256(CURRENT_POINTER),
            "current_entry_id": pointer.get("current_entry_id"),
            "root": str(current_root),
            "root_sha256": sha256(current_root),
            "directory": CURRENT_DIR,
            "objects": [TIMING_KEY, SCORE_TIME_KEY],
            "missing_reference_contract_object": REFERENCE_NCB_TIMING_KEY,
        },
        "ppg12_reference": {
            "timing_root": str(REFERENCE_TIMING),
            "timing_root_sha256": sha256(REFERENCE_TIMING),
            "score_root": str(REFERENCE_SCORE),
            "score_root_sha256": sha256(REFERENCE_SCORE),
            "objects": [TIMING_KEY, REFERENCE_NCB_TIMING_KEY, SCORE_TIME_KEY],
        },
        "processing": {
            "timing_profile": "mean time and standard error in each eta bin; residual is current minus PPG12",
            "score_timing": "each NCB-score bin normalized independently over cluster-MBD time; third panel is current minus PPG12 conditional fraction",
            "normalization": "no absolute normalization comparison",
        },
        "outputs": {
            "timing_png": str(timing_png),
            "timing_png_sha256": sha256(timing_png),
            "timing_csv": str(timing_csv),
            "timing_csv_sha256": sha256(timing_csv),
            "score_png": str(score_png),
            "score_png_sha256": sha256(score_png),
            "score_csv": str(score_csv),
            "score_csv_sha256": sha256(score_csv),
        },
        "limitations": [
            "The current ROOT does not contain h_npb_delta_t_mbd_vs_eta, so no current NCB-tagged eta-slice timing overlay is claimed.",
            "The score-timing comparison is conditional-shape parity, not an absolute-yield comparison.",
            "The PPG12 Figure 17 threshold scan remains screenshot-only in the local source set and is not numerically compared here.",
        ],
    }
    manifest_path = args.outdir / "the97_fullstat_ncb_timing_parity_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"timing_png": str(timing_png), "score_png": str(score_png), "manifest": str(manifest_path)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
