#!/usr/bin/env python3
"""Overlay PPG12 Fig.81 data-vertex reference with current RecoilJets output.

The PPG12 side is a compact JSON extract of
truth_vertex_reweight/output/{period}/reweight.root:h_D from SDCC.  The current
side is the explicit Fig.81 h_D reference-contract histogram family written by
RecoilJets before the analysis vz cut.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import ROOT  # type: ignore  # noqa: E402


REPO = Path(__file__).resolve().parents[3]
CAMPAIGN = "the93_ppg12_canonical_full_20260706_2145"
CURRENT_JSON = REPO / "dataOutput/current_recoiljets_artifacts/current/pp_data_merged/current.json"
REF_JSON = (
    REPO
    / "dataOutput/ppg12Parity"
    / CAMPAIGN
    / "data_vertex_fig81/sdcc_reference_extract_20260709"
    / "ppg12_fig81_current_reweight_hD_extract.clean.json"
)
OUT_DIR = (
    REPO
    / "dataOutput/ppg12Parity"
    / CAMPAIGN
    / "data_vertex_fig81/current_ppg12_20260422"
)

PERIODS = {
    "1p5mrad": {
        "title": "1.5 mrad reco vertex",
        "color": "#1f77b4",
        "hist": "PPG12_scaledtrigger30/h_ppg12_fig81_data_reco_z_slimtree_contract_1p5mrad",
        "z_cut": 83.0,
    },
    "0mrad": {
        "title": "0 mrad reco vertex",
        "color": "#d62728",
        "hist": "PPG12_scaledtrigger30/h_ppg12_fig81_data_reco_z_slimtree_contract_0mrad",
        "z_cut": 145.0,
    },
}

HIST_FAMILIES = {
    "hD_reference_contract": {
        "label": "This analysis",
        "note": "RecoilJets PPG12_scaledtrigger30/h_ppg12_fig81_data_reco_z_hD_reference_contract_{period}",
        "hist_template": "PPG12_scaledtrigger30/h_ppg12_fig81_data_reco_z_hD_reference_contract_{period}",
    },
    "slimtree_contract": {
        "label": "Current full-stat legacy post-vz slimtree proxy",
        "note": "RecoilJets PPG12_scaledtrigger30/h_ppg12_fig81_data_reco_z_slimtree_contract_{period}",
        "hist_template": "PPG12_scaledtrigger30/h_ppg12_fig81_data_reco_z_slimtree_contract_{period}",
    },
    "vtxqa_pre_vzcut": {
        "label": "Current full-stat bit30 pre-vzcut QA",
        "note": "RecoilJets PPG12_scaledtrigger30/h_ppg12_vtxqa_data_reco_z_triggered_{period}_pre_vzcut",
        "hist_template": "PPG12_scaledtrigger30/h_ppg12_vtxqa_data_reco_z_triggered_{period}_pre_vzcut",
    },
    "vtxqa_post_vzcut": {
        "label": "Current full-stat bit30 post-vzcut QA",
        "note": "RecoilJets PPG12_scaledtrigger30/h_ppg12_vtxqa_data_reco_z_triggered_{period}_post_vzcut",
        "hist_template": "PPG12_scaledtrigger30/h_ppg12_vtxqa_data_reco_z_triggered_{period}_post_vzcut",
    },
}


def load_current_root(path: Path | None) -> Path:
    if path is not None:
        return path
    with CURRENT_JSON.open() as f:
        current = json.load(f)
    roots = current.get("root_paths") or []
    if not roots:
        raise RuntimeError(f"no root_paths in {CURRENT_JSON}")
    return Path(roots[0])


def read_root_hist(root_path: Path, hist_path: str) -> dict[str, np.ndarray | float | int]:
    f = ROOT.TFile.Open(str(root_path), "READ")
    if not f or f.IsZombie() or f.TestBit(ROOT.TFile.kRecovered):
        raise RuntimeError(f"invalid current ROOT: {root_path}")
    h = f.Get(hist_path)
    if not h:
        raise RuntimeError(f"missing current histogram {hist_path} in {root_path}")
    nb = h.GetNbinsX()
    edges = np.array([h.GetXaxis().GetBinLowEdge(i) for i in range(1, nb + 2)], dtype=float)
    content = np.array([h.GetBinContent(i) for i in range(1, nb + 1)], dtype=float)
    error = np.array([h.GetBinError(i) for i in range(1, nb + 1)], dtype=float)
    entries = float(h.GetEntries())
    integral = float(h.Integral(1, nb))
    f.Close()
    return {
        "edges": edges,
        "content": content,
        "error": error,
        "entries": entries,
        "integral": integral,
        "nbins": nb,
    }


def read_ppg12_ref(ref: dict[str, Any], period: str) -> dict[str, Any]:
    rec = ref["periods"][period]
    bins = rec["bins"]
    edges = np.array([bins[0]["low"]] + [b["high"] for b in bins], dtype=float)
    content = np.array([b["content"] for b in bins], dtype=float)
    stored_err = np.array([b["error"] for b in bins], dtype=float)
    # The PPG12 h_D reference does not store useful bin errors in this ROOT.
    # For shape-ratio diagnostics, use Poisson raw-count errors.
    error = np.where(stored_err > 0.0, stored_err, np.sqrt(np.clip(content, 0.0, None)))
    return {
        "edges": edges,
        "content": content,
        "error": error,
        "integral": float(rec["integral_bins"]),
        "entries": float(rec["entries"]),
        "metadata": rec.get("metadata", {}),
        "path": rec["path"],
        "stat_mtime_local": rec.get("stat_mtime_local", ""),
    }


def shape(content: np.ndarray, error: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    total = float(np.sum(content))
    if total <= 0:
        raise RuntimeError("cannot normalize empty histogram")
    return content / total, error / total


def align_current_to_reference(current: dict[str, Any], reference_edges: np.ndarray) -> dict[str, Any]:
    edges = current["edges"]
    content = current["content"]
    error = current["error"]
    if len(edges) == len(reference_edges) and np.allclose(edges, reference_edges):
        return current
    if not (np.isclose(edges[0], reference_edges[0]) and np.isclose(edges[-1], reference_edges[-1])):
        raise RuntimeError("current/reference histogram ranges do not match")
    ref_n = len(reference_edges) - 1
    cur_n = len(edges) - 1
    if cur_n % ref_n != 0:
        raise RuntimeError(f"cannot rebin current {cur_n} bins onto reference {ref_n} bins")
    factor = cur_n // ref_n
    if not np.allclose(edges[::factor], reference_edges):
        raise RuntimeError("current bin edges are not an integer subdivision of reference edges")
    rebinned_content = content.reshape(ref_n, factor).sum(axis=1)
    rebinned_error = np.sqrt((error.reshape(ref_n, factor) ** 2).sum(axis=1))
    out = dict(current)
    out["edges"] = reference_edges.copy()
    out["content"] = rebinned_content
    out["error"] = rebinned_error
    out["rebin_factor_to_reference"] = factor
    out["nbins"] = ref_n
    return out


def ratio(num: np.ndarray, num_err: np.ndarray, den: np.ndarray, den_err: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    r = np.full_like(num, np.nan, dtype=float)
    e = np.full_like(num, np.nan, dtype=float)
    ok = (num > 0) & (den > 0)
    r[ok] = num[ok] / den[ok]
    e[ok] = r[ok] * np.sqrt((num_err[ok] / num[ok]) ** 2 + (den_err[ok] / den[ok]) ** 2)
    return r, e


def draw_panel(
    ax_top: plt.Axes,
    ax_bot: plt.Axes,
    period: str,
    ppg12: dict[str, Any],
    current: dict[str, Any],
    family_key: str,
) -> dict[str, Any]:
    cfg = PERIODS[period]
    family = HIST_FAMILIES[family_key]
    centers = 0.5 * (ppg12["edges"][:-1] + ppg12["edges"][1:])
    widths = ppg12["edges"][1:] - ppg12["edges"][:-1]

    p_shape, p_err = shape(ppg12["content"], ppg12["error"])
    c_shape, c_err = shape(current["content"], current["error"])
    r, r_err = ratio(c_shape, c_err, p_shape, p_err)

    color = cfg["color"]
    ax_top.set_yscale("log")
    ax_top.errorbar(
        centers,
        p_shape,
        yerr=p_err,
        xerr=0.5 * widths,
        fmt="o",
        ms=3.0,
        lw=0.8,
        color="black",
        mfc="white",
        mec="black",
        label="PPG12 reference",
    )
    ax_top.errorbar(
        centers,
        c_shape,
        yerr=c_err,
        xerr=0.0,
        fmt="o",
        ms=2.8,
        lw=0.8,
        color=color,
        label=family["label"],
    )
    z_cut = float(cfg["z_cut"])
    ax_top.axvline(-z_cut, color=color, ls=":", lw=0.8, alpha=0.7)
    ax_top.axvline(z_cut, color=color, ls=":", lw=0.8, alpha=0.7)
    ax_top.set_xlim(-200, 200)
    ax_top.set_ylim(1e-6, 3e-1)
    ax_top.grid(True, alpha=0.25)
    ax_top.text(0.06, 0.91, cfg["title"], transform=ax_top.transAxes, fontsize=10)
    ax_top.text(0.06, 0.82, r"$\it{\bf{sPHENIX}}$ Internal", transform=ax_top.transAxes, fontsize=8.5)
    ax_top.text(0.06, 0.74, r"$p+p$ $\sqrt{s}=200$ GeV", transform=ax_top.transAxes, fontsize=8.5)
    ax_top.text(0.06, 0.66, "Photon-4-GeV trigger; before vertex selection", transform=ax_top.transAxes, fontsize=7.8)
    ax_top.legend(loc="upper right", frameon=False, fontsize=8)

    ax_bot.axhline(1.0, color="0.35", lw=0.9, ls="--")
    ax_bot.errorbar(centers, r, yerr=r_err, fmt="o", ms=2.8, lw=0.8, color=color)
    ax_bot.set_xlim(-200, 200)
    ax_bot.set_ylim(0.0, 6.2)
    ax_bot.grid(True, alpha=0.25)
    ax_bot.set_xlabel(r"$z_{\rm reco}$ (cm)")
    ax_bot.set_ylabel("This analysis / PPG12")

    return {
        "period": period,
        "ppg12_raw_integral": ppg12["integral"],
        "current_raw_integral": current["integral"],
        "raw_population_ratio_current_over_ppg12": current["integral"] / ppg12["integral"],
        "ppg12_path": ppg12["path"],
        "ppg12_mtime": ppg12["stat_mtime_local"],
        "ppg12_metadata": ppg12["metadata"],
        "current_hist": str(family["hist_template"]).format(period=period),
        "ratio_min": float(np.nanmin(r)),
        "ratio_max": float(np.nanmax(r)),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--current-root", type=Path, default=None)
    parser.add_argument("--reference-json", type=Path, default=REF_JSON)
    parser.add_argument("--out-dir", type=Path, default=OUT_DIR)
    parser.add_argument("--current-family", choices=sorted(HIST_FAMILIES), default="hD_reference_contract")
    args = parser.parse_args()

    family = HIST_FAMILIES[args.current_family]
    out_dir = args.out_dir
    if args.out_dir == OUT_DIR:
        out_dir = args.out_dir / args.current_family
    out_dir.mkdir(parents=True, exist_ok=True)

    current_root = load_current_root(args.current_root)
    with args.reference_json.open() as f:
        ref = json.load(f)

    fig, axes = plt.subplots(
        2,
        2,
        figsize=(12.0, 7.2),
        sharex="col",
        gridspec_kw={"height_ratios": [3.2, 1.15], "hspace": 0.05, "wspace": 0.24},
    )
    summaries = []
    csv_rows = ["period,bin,low,high,center,ppg12_raw,current_raw,ppg12_shape,current_shape,current_over_ppg12_shape"]

    for col, period in enumerate(["1p5mrad", "0mrad"]):
        ppg12 = read_ppg12_ref(ref, period)
        current_hist = str(family["hist_template"]).format(period=period)
        current = align_current_to_reference(read_root_hist(current_root, current_hist), ppg12["edges"])
        summary = draw_panel(axes[0, col], axes[1, col], period, ppg12, current, args.current_family)
        summaries.append(summary)

        centers = 0.5 * (ppg12["edges"][:-1] + ppg12["edges"][1:])
        p_shape, p_err = shape(ppg12["content"], ppg12["error"])
        c_shape, c_err = shape(current["content"], current["error"])
        r, _ = ratio(c_shape, c_err, p_shape, p_err)
        for ib, center in enumerate(centers):
            csv_rows.append(
                f"{period},{ib+1},{ppg12['edges'][ib]:.6g},{ppg12['edges'][ib+1]:.6g},{center:.6g},"
                f"{ppg12['content'][ib]:.12g},{current['content'][ib]:.12g},"
                f"{p_shape[ib]:.12g},{c_shape[ib]:.12g},{r[ib]:.12g}"
            )

    axes[0, 0].set_ylabel("shape-normalized events")
    axes[1, 0].set_ylabel("This analysis / PPG12")
    for ax in axes.ravel():
        ax.tick_params(direction="in", top=True, right=True)

    fig.suptitle("Reconstructed-vertex consistency", y=0.985, fontsize=15, fontweight="bold")
    stem = f"fig81_data_reco_vertex_sdcc_vs_current_fullstat_{args.current_family}_ratio"
    png = out_dir / f"{stem}.png"
    csv = out_dir / f"{stem}.csv"
    manifest = out_dir / f"{stem}_manifest.json"
    fig.savefig(png, dpi=180, bbox_inches="tight")
    csv.write_text("\n".join(csv_rows) + "\n")
    manifest.write_text(
        json.dumps(
            {
                "artifact": f"PPG12 Fig.81 vertex overlay using current SDCC reweight.root h_D and current full-stat {args.current_family} hists",
                "current_root": str(current_root),
                "current_artifact_pointer": str(CURRENT_JSON),
                "reference_json": str(args.reference_json),
                "ratio_definition": "shape-normalized current output / shape-normalized PPG12 current h_D",
                "ppg12_reference_note": "SDCC truth_vertex_reweight/output/{period}/reweight.root:h_D, not local reweight.root.*_backup",
                "current_family": args.current_family,
                "current_note": family["note"],
                "periods": summaries,
                "outputs": {"png": str(png), "csv": str(csv), "manifest": str(manifest)},
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    print(png)
    print(manifest)


if __name__ == "__main__":
    main()
