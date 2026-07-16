#!/usr/bin/env python3
"""Draw a PPG12 Fig. 36-style photon response matrix from the current pp photon+jet ROOT.

The current RecoilJets output stores the response matrix after normal pp MC
weights, vertex weights, and SI/DI mix weights. PPG12's Fig. 36 caption refers
to an additional truth-pT prior reweight applied while filling the response.
This script can draw the raw stored matrix or an offline bin-center application
of the current PPG12 truth-prior formula for a visual parity check.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import uproot
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import LogNorm


REPO = Path(__file__).resolve().parents[3]
CURRENT_POINTER = REPO / "dataOutput/current_recoiljets_artifacts/current/pp_sim_photonjet_merged/current.json"
OUT_DIR = REPO / "dataOutput/ppg12Parity/the76_ppg12_fig24_photonjet_fix_20260702_014217/fig36_response_matrix"

# PPG12 config_bdt_nom.yaml analysis.unfold.truth_reweight, re-derived 2026-04-27.
TRUTH_PRIOR_PARAMS = [4.05592, -0.984728, -0.478818, 0.0723232, 0.0522681]

# ROOT kBird-like blue -> cyan/green -> yellow palette. This is intentionally
# close to the PPG12/ROOT heatmap look without requiring a ROOT plotting macro.
KBIRD_LIKE = LinearSegmentedColormap.from_list(
    "root_kbird_like",
    [
        (0.00, "#00004f"),
        (0.18, "#0033a0"),
        (0.36, "#1178bd"),
        (0.55, "#2fb7b5"),
        (0.74, "#8bd3a8"),
        (0.90, "#e7ef8a"),
        (1.00, "#ffffcc"),
    ],
)


def ppg12_truth_prior_weight(x: np.ndarray) -> np.ndarray:
    p0, p1, p2, p3, p4 = TRUTH_PRIOR_PARAMS
    return (p0 + p1 * x + p3 * x * x) / (1.0 + p2 * x + p4 * x * x)


def resolve_current_root() -> tuple[Path, dict]:
    with CURRENT_POINTER.open() as f:
        meta = json.load(f)
    roots = [Path(p) for p in meta.get("root_paths", [])]
    if not roots:
        raise RuntimeError(f"no root_paths in {CURRENT_POINTER}")
    root = roots[0]
    if not root.exists():
        raise FileNotFoundError(root)
    return root, meta


def load_response(root: Path, object_name: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    with uproot.open(root) as f:
        if object_name not in f:
            raise KeyError(f"{object_name} not found in {root}")
        h = f[object_name]
        values = h.values(flow=False).astype(float)
        xedges = h.axis(0).edges()
        yedges = h.axis(1).edges()
    return values, xedges, yedges


def apply_prior(values: np.ndarray, yedges: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    ycenters = 0.5 * (yedges[:-1] + yedges[1:])
    weights = ppg12_truth_prior_weight(ycenters)
    # uproot values shape is [xbin, ybin]. Apply truth-y prior per truth bin.
    return values * weights[np.newaxis, :], weights


def draw_matrix(values: np.ndarray, xedges: np.ndarray, yedges: np.ndarray, out_png: Path, *, reweighted: bool) -> None:
    mpl.rcParams.update({
        "font.family": "DejaVu Sans",
        "mathtext.fontset": "dejavusans",
        "axes.linewidth": 1.1,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
    })

    positive = values[values > 0]
    vmax = max(float(np.nanmax(positive)) if positive.size else 1.0, 10.0)
    fig = plt.figure(figsize=(4.6, 4.6), dpi=180)
    ax = fig.add_axes([0.18, 0.15, 0.66, 0.76])
    cax = fig.add_axes([0.86, 0.15, 0.035, 0.76])

    mesh = ax.pcolormesh(
        xedges,
        yedges,
        values.T,
        cmap=KBIRD_LIKE,
        norm=LogNorm(vmin=1.0, vmax=max(vmax, 1.0e7)),
        shading="flat",
    )
    cb = fig.colorbar(mesh, cax=cax)
    cb.set_ticks([1, 10, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7])
    cb.ax.tick_params(labelsize=8, direction="in")

    ax.set_xlim(10, 36)
    ax.set_ylim(8, 45)
    ax.set_xticks([10, 15, 20, 25, 30, 35])
    ax.set_yticks([10, 15, 20, 25, 30, 35, 40, 45])
    ax.minorticks_on()
    ax.tick_params(which="major", length=6, width=1.0, labelsize=9)
    ax.tick_params(which="minor", length=3, width=0.8)
    ax.set_xlabel(r"$E_T^{\gamma,\mathrm{rec}}\ \mathrm{[GeV]}$", fontsize=11, loc="right")
    ax.set_ylabel(r"$E_T^{\gamma,\mathrm{truth}}\ \mathrm{[GeV]}$", fontsize=11)

    ax.text(0.08, 0.93, r"$\it{\bf{sPHENIX}}$ Internal", transform=ax.transAxes,
            ha="left", va="top", fontsize=8.5)
    ax.text(0.08, 0.86, r"$p{+}p\ \sqrt{s}=200\ \mathrm{GeV}$", transform=ax.transAxes,
            ha="left", va="top", fontsize=7.5)
    label = "response matrix"
    if reweighted:
        label += "\noffline truth-prior reweighted"
    ax.text(0.08, 0.76, label, transform=ax.transAxes, ha="left", va="top", fontsize=7)

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=180)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--object", default="SIM/h_response_full_0")
    parser.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = parser.parse_args()

    root, meta = resolve_current_root()
    values, xedges, yedges = load_response(root, args.object)
    prior_values, prior_weights = apply_prior(values, yedges)

    raw_png = args.out_dir / "fig36_current_h_response_full_raw_ppg12_style.png"
    prior_png = args.out_dir / "fig36_current_h_response_full_offline_truth_prior_reweighted_ppg12_style.png"
    draw_matrix(values, xedges, yedges, raw_png, reweighted=False)
    draw_matrix(prior_values, xedges, yedges, prior_png, reweighted=True)

    manifest = {
        "script": str(Path(__file__).resolve()),
        "current_pointer": str(CURRENT_POINTER),
        "campaign_tag": meta.get("campaign_tag"),
        "root": str(root),
        "object": args.object,
        "raw_png": str(raw_png),
        "offline_prior_reweighted_png": str(prior_png),
        "x_edges": xedges.tolist(),
        "y_edges": yedges.tolist(),
        "raw_sum": float(np.sum(values)),
        "prior_reweighted_sum": float(np.sum(prior_values)),
        "prior_formula": "([0] + [1]*x + [3]*x*x) / (1 + [2]*x + [4]*x*x)",
        "prior_params": TRUTH_PRIOR_PARAMS,
        "prior_weight_by_truth_bin_center": prior_weights.tolist(),
        "caveat": (
            "The current ROOT stores h_response_full_0 from RecoilJets. This script applies the "
            "PPG12 truth-pT prior formula offline at truth-bin centers for a Fig.36-style visual. "
            "Exact per-candidate continuous prior reweighting would need to be applied during filling."
        ),
    }
    args.out_dir.mkdir(parents=True, exist_ok=True)
    with (args.out_dir / "fig36_response_matrix_manifest.json").open("w") as f:
        json.dump(manifest, f, indent=2, sort_keys=True)
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
