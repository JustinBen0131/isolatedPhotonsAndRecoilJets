#!/usr/bin/env python3
"""Overlay PPG12 SDCC stitch spectra with current RecoilJets pp outputs.

The PPG12 points come from the validated SDCC ROOT extraction CSV. The current
pp side reads diagnostic histograms from locally available RecoilJets outputs.

Important: final ``*plus*MERGED.root`` products are not a PPG12 stitch-contract
source. They are useful for old one-off shape diagnostics, but they silently
mix downstream merge conventions with the Fig. 5/6 source-histogram contract.
By default this helper now refuses those paths. Set
``RJ_ALLOW_FINAL_MERGED_STITCH_PROXY=1`` only to reproduce the historical proxy
plots, and use ``audit_pp_stitching_dual_contracts.py`` for the current
contract/parity audit.
"""

from __future__ import annotations

import csv
import json
import math
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import uproot


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
OUT_DIR = (
    REPO
    / "dataOutput/the85_auau_xjgamma_unfolding_push/"
    / "ppg12_stitching_fig5_fig6_validation"
)
PPG12_SDCC_CSV = (
    REPO
    / "dataOutput/ppPhotonMLPipeline/"
    / "ppg12_basev3E_currentIAN_fullsim_20260521_1811/"
    / "validation/currentIAN_stitching/ppg12_currentian_rootfit_stitch_points.csv"
)


PANELS = {
    "photonjet": {
        "group": "photon",
        "out_png": OUT_DIR / "ppg12_sdcc_vs_current_pp_photonjet_stitch_overlay_ratio.png",
        "manifest": OUT_DIR / "ppg12_sdcc_vs_current_pp_photonjet_stitch_overlay_manifest.json",
        "current_only_png": OUT_DIR / "current_pp_photonjet_ppg12_color_stitch_shape_ratio.png",
        "current_only_manifest": OUT_DIR / "current_pp_photonjet_ppg12_color_stitch_shape_manifest.json",
        "current_root": (
            REPO
            / "InputFiles/pp24/ppg12_globalmbd_mbddigi_componentmix_20260629/"
            / "jetMinPtScan_dphiScan_vz60_isoR40_isSliding_preselectionNewPPG12_tightNewPPG12_nonTightNewPPG12/"
            / "photonJet5and10and20componentmix_SIM/"
            / "RecoilJets_photonjet5plus10plus20_componentmix_MERGED.root"
        ),
        "current_hist": "SIM/h_ppPhotonStitch_ppg12Fig5_maxPhotonPt_kept",
        "xlim": (10.0, 40.0),
        "ylim": (2.0e-2, 2.0e3),
        "xlabel": r"Leading $E_T^\gamma$ [GeV]",
        "ylabel": r"$d\sigma/dE_T^\gamma$ [pb / GeV]",
        "ratio_ylim": (0.85, 1.15),
        "fit_range": (10.0, 36.0),
        "label": "photon+jet",
        "reference_note": "PPG12 SDCC Fig. 5 photon 5/10/20 GeV stitched leading-truth-photon spectrum",
        "current_note": "latest local current pp photon+jet component-mix output found on 2026-06-30",
        "samples": [
            {"name": "photon 5", "low": 10.0, "high": 14.0, "color": "#d62aa0"},
            {"name": "photon 10", "low": 14.0, "high": 22.0, "color": "#2ca02c"},
            {"name": "photon 20", "low": 22.0, "high": 40.0, "color": "#1296f3"},
        ],
    },
    "inclusivejet": {
        "group": "jet",
        "out_png": OUT_DIR / "ppg12_sdcc_vs_current_pp_inclusivejet_stitch_overlay_ratio.png",
        "manifest": OUT_DIR / "ppg12_sdcc_vs_current_pp_inclusivejet_stitch_overlay_manifest.json",
        "current_only_png": OUT_DIR / "current_pp_inclusivejet_ppg12_color_stitch_shape_ratio.png",
        "current_only_manifest": OUT_DIR / "current_pp_inclusivejet_ppg12_color_stitch_shape_manifest.json",
        "current_root": (
            REPO
            / "dataOutput/ppg12PhotonYield/THE76_ppg12_photon_yield_v1_sim_20260616/"
            / "merged_roots/jetMinPtScan_dphiScan_vz60_isoR40_isSliding_preselectionNewPPG12_tightNewPPG12_nonTightNewPPG12/"
            / "inclusiveJet5to40_SIM/"
            / "RecoilJets_jet5plus8plus12plus20plus30plus40_MERGED.root"
        ),
        "current_hist": "SIM/h_ppInclusiveJetStitch_ppg12TruthSpectrum_r04_maxTruthJetPt_kept",
        "xlim": (9.0, 50.0),
        "ylim": (1.0e4, 1.0e12),
        "xlabel": r"Leading $p_T^\mathrm{jet}$ [GeV]",
        "ylabel": "counts",
        "ratio_ylim": (0.85, 1.15),
        "fit_range": (10.0, 50.0),
        "label": "inclusive+jet",
        "reference_note": "PPG12 SDCC Fig. 6 jet 8/12/20/30/40 GeV stitched leading-truth-jet spectrum",
        "current_note": "latest local current pp inclusive-jet output found on 2026-06-30",
        "samples": [
            {"name": "jet 8", "low": 9.0, "high": 14.0, "color": "#d62aa0"},
            {"name": "jet 12", "low": 14.0, "high": 21.0, "color": "#2ca02c"},
            {"name": "jet 20", "low": 21.0, "high": 32.0, "color": "#1296f3"},
            {"name": "jet 30", "low": 32.0, "high": 42.0, "color": "#ff6f00"},
            {"name": "jet 40", "low": 42.0, "high": 50.0, "color": "#d62aa0"},
        ],
    },
}

CANVAS_PX = (1544, 1996)
DPI = 220
FIGSIZE = (CANVAS_PX[0] / DPI, CANVAS_PX[1] / DPI)
ALLOW_FINAL_MERGED_PROXY = os.environ.get("RJ_ALLOW_FINAL_MERGED_STITCH_PROXY", "0") in {
    "1",
    "true",
    "TRUE",
    "yes",
    "YES",
}


def reject_final_merged_proxy(path: Path) -> None:
    text = str(path)
    if ALLOW_FINAL_MERGED_PROXY:
        return
    if "plus" in path.name.lower() or path.name.endswith("_MERGED.root"):
        raise RuntimeError(
            "Refusing final merged pp stitch proxy ROOT for a PPG12 parity plot: "
            f"{path}. Use the in-situ contract CSV audit instead, or set "
            "RJ_ALLOW_FINAL_MERGED_STITCH_PROXY=1 only to reproduce the old diagnostic."
        )
    if "plus" in text and "MERGED.root" in text:
        raise RuntimeError(
            "Refusing final merged pp stitch proxy ROOT for a PPG12 parity plot: "
            f"{path}. Use scripts/diagnostics/pp_shuhang/audit_pp_stitching_dual_contracts.py."
        )


def load_reference(group: str, xlim: tuple[float, float]) -> dict[float, tuple[float, float]]:
    points: dict[float, tuple[float, float]] = {}
    with PPG12_SDCC_CSV.open(newline="") as f:
        for row in csv.DictReader(f):
            if row["group"] != group or row["sample"] != "stitched" or row["used_in_stitch"] != "1":
                continue
            x = float(row["bin_center"])
            y = float(row["value"])
            if xlim[0] <= x <= xlim[1] and y > 0:
                points[round(x, 6)] = (y, float(row["error"]))
    if not points:
        raise RuntimeError(f"No PPG12 SDCC reference points found for group={group}")
    return points


def load_reference_rows(group: str, xlim: tuple[float, float]) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    with PPG12_SDCC_CSV.open(newline="") as f:
        for row in csv.DictReader(f):
            if row["group"] != group or row["used_in_stitch"] != "1":
                continue
            x = float(row["bin_center"])
            y = float(row["value"])
            if xlim[0] <= x <= xlim[1] and y > 0:
                rows.append(
                    {
                        "sample": row["sample"],
                        "x": x,
                        "y": y,
                        "ey": float(row["error"]),
                        "fit": float(row["root_fit_value"] or "nan"),
                        "fit_params": row["root_fit_params"],
                        "ratio_label": row["ratio_label"],
                    }
                )
    if not rows:
        raise RuntimeError(f"No PPG12 SDCC reference rows found for group={group}")
    return rows


def load_fit_contract(group: str) -> tuple[str, str]:
    with PPG12_SDCC_CSV.open(newline="") as f:
        for row in csv.DictReader(f):
            if row["group"] == group and row["sample"] == "stitched" and row["root_fit_params"]:
                return row["root_fit_params"], row["ratio_label"]
    raise RuntimeError(f"No stitched fit parameters found for group={group}")


def eval_root_hagedorn(grid: np.ndarray, params: str) -> np.ndarray:
    p = np.array([float(x) for x in params.split(";")], dtype=float)
    return p[0] * np.power(p[1] / grid, p[2] + p[3] * np.log(grid / p[1]) + p[4] * grid)


def load_current(path: Path, hist_name: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    reject_final_merged_proxy(path)
    if not path.exists():
        raise FileNotFoundError(path)
    with uproot.open(path) as f:
        if hist_name not in f:
            raise KeyError(f"{path} does not contain {hist_name}")
        hist = f[hist_name]
        values, edges = hist.to_numpy(flow=False)
        variances = hist.variances(flow=False)
    centers = 0.5 * (edges[:-1] + edges[1:])
    if variances is None:
        variances = np.clip(values, 0, None)
    return centers.astype(float), values.astype(float), np.asarray(variances, dtype=float)


def matched_arrays(panel: dict[str, object]) -> dict[str, np.ndarray]:
    group = str(panel["group"])
    xlim = panel["xlim"]
    assert isinstance(xlim, tuple)
    ref_points = load_reference(group, xlim)
    centers, current, current_var = load_current(Path(panel["current_root"]), str(panel["current_hist"]))

    xs: list[float] = []
    ref: list[float] = []
    ref_err: list[float] = []
    cur: list[float] = []
    cur_err: list[float] = []
    for x, y, v in zip(centers, current, current_var):
        key = round(float(x), 6)
        if key not in ref_points or y <= 0:
            continue
        xs.append(float(x))
        ref_y, ref_ey = ref_points[key]
        ref.append(ref_y)
        ref_err.append(ref_ey)
        cur.append(float(y))
        cur_err.append(math.sqrt(max(float(v), 0.0)))

    if not xs:
        raise RuntimeError(f"No matching bins between {panel['current_root']} and PPG12 reference")

    order = np.argsort(xs)
    return {
        "x": np.asarray(xs, dtype=float)[order],
        "ref": np.asarray(ref, dtype=float)[order],
        "ref_err": np.asarray(ref_err, dtype=float)[order],
        "current": np.asarray(cur, dtype=float)[order],
        "current_err": np.asarray(cur_err, dtype=float)[order],
    }


def robust_scale(ref: np.ndarray, current: np.ndarray) -> float:
    mask = (ref > 0) & (current > 0) & np.isfinite(ref) & np.isfinite(current)
    if not np.any(mask):
        raise RuntimeError("Cannot compute current-to-reference scale")
    return float(np.median(ref[mask] / current[mask]))


def sample_mask(x: np.ndarray, sample: dict[str, object]) -> np.ndarray:
    low = float(sample["low"])
    high = float(sample["high"])
    return (x >= low) & (x < high)


def ratio_limits(ratio: np.ndarray) -> tuple[float, float]:
    finite = ratio[np.isfinite(ratio)]
    if len(finite) == 0:
        return (0.75, 1.25)
    lo = max(0.65, float(np.nanpercentile(finite, 1)) - 0.06)
    hi = min(1.35, float(np.nanpercentile(finite, 99)) + 0.06)
    if hi - lo < 0.18:
        mid = 0.5 * (hi + lo)
        lo, hi = mid - 0.09, mid + 0.09
    return (lo, hi)


def panel_ratio_ylim(panel: dict[str, object]) -> tuple[float, float]:
    ratio_ylim = panel.get("ratio_ylim", (0.85, 1.15))
    assert isinstance(ratio_ylim, tuple)
    return ratio_ylim


def configure_axes(fig: plt.Figure, panel: dict[str, object]) -> tuple[plt.Axes, plt.Axes]:
    gs = fig.add_gridspec(
        2,
        1,
        height_ratios=(3.35, 1.05),
        hspace=0.035,
        left=0.18,
        right=0.955,
        top=0.965,
        bottom=0.105,
    )
    ax = fig.add_subplot(gs[0])
    rax = fig.add_subplot(gs[1], sharex=ax)
    ax.set_yscale("log")
    ax.set_xlim(*panel["xlim"])
    ax.set_ylim(*panel["ylim"])
    ax.tick_params(labelbottom=False, which="both", length=7, labelsize=15)
    ax.tick_params(which="minor", length=3.5)
    rax.tick_params(which="both", length=7, labelsize=15)
    rax.tick_params(which="minor", length=3.5)
    ax.set_ylabel(str(panel["ylabel"]), fontsize=18)
    rax.set_xlabel(str(panel["xlabel"]), fontsize=18)
    return ax, rax


def set_plot_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "mathtext.fontset": "dejavusans",
            "axes.linewidth": 1.25,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )


def draw_panel(name: str, panel: dict[str, object]) -> dict[str, object]:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    arrays = matched_arrays(panel)
    ref_rows = load_reference_rows(str(panel["group"]), panel["xlim"])
    fit_params, ratio_label = load_fit_contract(str(panel["group"]))
    x = arrays["x"]
    ref = arrays["ref"]
    cur = arrays["current"]
    cur_err = arrays["current_err"]
    scale = robust_scale(ref, cur)
    cur_scaled = cur * scale
    cur_scaled_err = cur_err * scale
    fit_at_current = eval_root_hagedorn(x, fit_params)
    ratio = cur_scaled / fit_at_current
    ratio_err = cur_scaled_err / fit_at_current

    set_plot_style()

    fig = plt.figure(figsize=FIGSIZE, dpi=DPI)
    ax, rax = configure_axes(fig, panel)

    fit_grid = np.linspace(float(panel["fit_range"][0]), float(panel["xlim"][1]), 900)
    fit_curve = eval_root_hagedorn(fit_grid, fit_params)
    ax.plot(fit_grid, fit_curve, color="red", lw=1.7, zorder=2)

    samples = list(panel["samples"])
    for sample in samples:
        color = str(sample["color"])
        sample_rows = [
            r
            for r in ref_rows
            if r["sample"] != "stitched" and sample_mask(np.asarray([float(r["x"])]), sample)[0]
        ]
        if sample_rows:
            sx = np.asarray([float(r["x"]) for r in sample_rows], dtype=float)
            sy = np.asarray([float(r["y"]) for r in sample_rows], dtype=float)
            sey = np.asarray([float(r["ey"]) for r in sample_rows], dtype=float)
            sfit = eval_root_hagedorn(sx, fit_params)
            ax.errorbar(
                sx,
                sy,
                yerr=sey,
                fmt="o",
                ms=4.5,
                mfc="white",
                mec=color,
                mew=1.05,
                ecolor=color,
                elinewidth=0.75,
                capsize=0,
                linestyle="none",
                zorder=4,
            )
            rax.errorbar(
                sx,
                sy / sfit,
                yerr=sey / sfit,
                fmt="o",
                ms=4.1,
                mfc="white",
                mec=color,
                mew=1.0,
                ecolor=color,
                elinewidth=0.7,
                capsize=0,
                linestyle="none",
                zorder=4,
            )

        cmask = sample_mask(x, sample)
        if np.any(cmask):
            ax.errorbar(
                x[cmask],
                cur_scaled[cmask],
                yerr=cur_scaled_err[cmask],
                fmt="s",
                ms=4.1,
                mfc=color,
                mec=color,
                mew=0.8,
                ecolor=color,
                elinewidth=0.7,
                capsize=0,
                linestyle="none",
                zorder=5,
            )
            rax.errorbar(
                x[cmask],
                ratio[cmask],
                yerr=ratio_err[cmask],
                fmt="s",
                ms=3.8,
                mfc=color,
                mec=color,
                mew=0.7,
                ecolor=color,
                elinewidth=0.65,
                capsize=0,
                linestyle="none",
                zorder=5,
            )

    ax.text(
        0.97,
        0.92,
        r"$\it{\bf{sPHENIX}}$ Internal",
        transform=ax.transAxes,
        fontsize=18,
        ha="right",
        va="top",
    )
    ax.text(
        0.97,
        0.84,
        f"{panel['label']} stitch comparison",
        transform=ax.transAxes,
        fontsize=15,
        ha="right",
        va="top",
    )
    ax.text(
        0.97,
        0.77,
        "current output scaled by one constant",
        transform=ax.transAxes,
        fontsize=13.5,
        ha="right",
        va="top",
    )
    from matplotlib.lines import Line2D

    ax.legend(
        handles=[
            Line2D([0], [0], marker="o", color="black", mfc="white", mec="black", lw=0, label="PPG12 SDCC source"),
            Line2D([0], [0], marker="s", color="black", mfc="black", mec="black", lw=0, label="Current pp output"),
            Line2D([0], [0], color="red", lw=1.7, label="PPG12 fit"),
        ],
        loc="lower left",
        frameon=False,
        fontsize=12.5,
        handletextpad=0.55,
        borderpad=0.2,
    )
    rax.axhline(1.0, color="0.45", lw=1.0, ls="--", zorder=1)
    rax.set_ylim(*panel_ratio_ylim(panel))
    rax.set_ylabel(ratio_label, fontsize=14.5)

    out_png = Path(panel["out_png"])
    fig.savefig(out_png)
    plt.close(fig)

    summary = {
        "name": name,
        "output_png": str(out_png),
        "ppg12_sdcc_csv": str(PPG12_SDCC_CSV),
        "current_root": str(panel["current_root"]),
        "current_hist": str(panel["current_hist"]),
        "reference_note": panel["reference_note"],
        "current_note": panel["current_note"],
        "normalization": {
            "mode": "single robust constant scale",
            "scale_applied_to_current": scale,
            "definition": "median(PPG12_SDCC / current_ROOT) over matched displayed bins",
        },
        "matched_bins": int(len(x)),
        "x_range": list(panel["xlim"]),
        "fit": {
            "params": fit_params,
            "ratio_label": ratio_label,
        },
        "ratio_current_scaled_over_fit": {
            "median": float(np.nanmedian(ratio)),
            "mean": float(np.nanmean(ratio)),
            "min": float(np.nanmin(ratio)),
            "max": float(np.nanmax(ratio)),
            "max_abs_deviation_from_unity": float(np.nanmax(np.abs(ratio - 1.0))),
        },
    }
    manifest = Path(panel["manifest"])
    manifest.write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def draw_current_only_panel(name: str, panel: dict[str, object]) -> dict[str, object]:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    arrays = matched_arrays(panel)
    fit_params, ratio_label = load_fit_contract(str(panel["group"]))
    x = arrays["x"]
    ref = arrays["ref"]
    cur = arrays["current"]
    cur_err = arrays["current_err"]
    scale = robust_scale(ref, cur)
    cur_scaled = cur * scale
    cur_scaled_err = cur_err * scale
    fit_at_current = eval_root_hagedorn(x, fit_params)
    ratio = cur_scaled / fit_at_current
    ratio_err = cur_scaled_err / fit_at_current

    set_plot_style()
    fig = plt.figure(figsize=FIGSIZE, dpi=DPI)
    ax, rax = configure_axes(fig, panel)

    fit_grid = np.linspace(float(panel["fit_range"][0]), float(panel["xlim"][1]), 900)
    fit_curve = eval_root_hagedorn(fit_grid, fit_params)
    ax.plot(fit_grid, fit_curve, color="red", lw=1.7, zorder=2)

    samples = list(panel["samples"])
    for sample in samples:
        mask = sample_mask(x, sample)
        if not np.any(mask):
            continue
        color = str(sample["color"])
        ax.errorbar(
            x[mask],
            cur_scaled[mask],
            yerr=cur_scaled_err[mask],
            fmt="o",
            ms=4.9,
            mfc=color,
            mec=color,
            mew=0.8,
            ecolor=color,
            elinewidth=0.8,
            capsize=0,
            linestyle="none",
            label=str(sample["name"]),
            zorder=3,
        )
        rax.errorbar(
            x[mask],
            ratio[mask],
            yerr=ratio_err[mask],
            fmt="o",
            ms=4.4,
            mfc=color,
            mec=color,
            mew=0.7,
            ecolor=color,
            elinewidth=0.75,
            capsize=0,
            linestyle="none",
            zorder=3,
        )

    ax.text(
        0.97,
        0.92,
        r"$\it{\bf{sPHENIX}}$ Internal",
        transform=ax.transAxes,
        fontsize=18,
        ha="right",
        va="top",
    )
    ax.text(
        0.97,
        0.84,
        f"current pp {panel['label']} output",
        transform=ax.transAxes,
        fontsize=15,
        ha="right",
        va="top",
    )
    ax.text(
        0.97,
        0.77,
        "PPG12 stitch-window colors",
        transform=ax.transAxes,
        fontsize=13.5,
        ha="right",
        va="top",
    )
    legend_kwargs = {
        "frameon": False,
        "fontsize": 12.5,
        "handletextpad": 0.45,
        "borderpad": 0.2,
        "columnspacing": 0.9,
    }
    if name == "inclusivejet":
        ax.legend(loc="lower left", ncol=2, **legend_kwargs)
    else:
        ax.legend(loc="lower left", **legend_kwargs)

    rax.axhline(1.0, color="0.45", lw=1.0, ls="--", zorder=1)
    rax.set_ylim(*panel_ratio_ylim(panel))
    rax.set_ylabel(ratio_label, fontsize=14.5)

    out_png = Path(panel["current_only_png"])
    fig.savefig(out_png)
    plt.close(fig)

    sample_summaries = []
    for sample in samples:
        mask = sample_mask(x, sample)
        if not np.any(mask):
            continue
        sample_ratio = ratio[mask]
        sample_summaries.append(
            {
                "sample": sample["name"],
                "window": [sample["low"], sample["high"]],
                "color": sample["color"],
                "bins": int(np.sum(mask)),
                "ratio_median": float(np.nanmedian(sample_ratio)),
                "ratio_min": float(np.nanmin(sample_ratio)),
                "ratio_max": float(np.nanmax(sample_ratio)),
            }
        )

    summary = {
        "name": name,
        "output_png": str(out_png),
        "ppg12_sdcc_csv": str(PPG12_SDCC_CSV),
        "current_root": str(panel["current_root"]),
        "current_hist": str(panel["current_hist"]),
        "normalization": {
            "mode": "single robust constant scale",
            "scale_applied_to_current": scale,
            "definition": "median(PPG12_SDCC / current_ROOT) over matched displayed bins",
        },
        "matched_bins": int(len(x)),
        "x_range": list(panel["xlim"]),
        "fit": {
            "params": fit_params,
            "ratio_label": ratio_label,
        },
        "sample_windows": sample_summaries,
        "ratio_current_scaled_over_fit": {
            "median": float(np.nanmedian(ratio)),
            "mean": float(np.nanmean(ratio)),
            "min": float(np.nanmin(ratio)),
            "max": float(np.nanmax(ratio)),
            "max_abs_deviation_from_unity": float(np.nanmax(np.abs(ratio - 1.0))),
        },
    }
    Path(panel["current_only_manifest"]).write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def main() -> None:
    summaries = []
    for name, panel in PANELS.items():
        summaries.append(draw_panel(name, panel))
        summaries.append(draw_current_only_panel(name, panel))
    print(json.dumps(summaries, indent=2))


if __name__ == "__main__":
    main()
