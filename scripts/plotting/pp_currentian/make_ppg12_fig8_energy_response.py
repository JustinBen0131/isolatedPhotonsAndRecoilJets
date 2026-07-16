#!/usr/bin/env python3
"""Render a PPG12 IAN Fig. 8-style photon energy-response plot.

Input is a RecoilJets ROOT file written with RJ_PPG12_PHOTON_YIELD=1.  The
target histogram is filled in pp SIM from truth-matched signal photons:

  x = truth photon E_T
  y = reconstructed cluster E_T / truth photon E_T

The input histogram is finer and can contain far more weighted entries than
the historical IAN panel.  Display-only rebinning and a clearly recorded
global shape normalization are available to reproduce the IAN presentation
without changing the underlying response distribution.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np

try:
    import matplotlib.pyplot as plt
    import uproot
    from matplotlib import colors
except ModuleNotFoundError as exc:
    raise SystemExit(
        "Missing plotting dependency. Run with "
        "/Users/patsfan753/Desktop/analysis/env/bin/python3 or another environment "
        "that provides matplotlib and uproot."
    ) from exc


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
DEFAULT_OUTDIR = REPO / "dataOutput/ppg12PhotonYield/ppg12_fig8_energy_response"
DEFAULT_HIST = "h2_ppg12_fig8_respET_vs_truthET_reco_vz30_eta07"
DEFAULT_AUDIT = "h_ppg12_fig8_audit_flow"


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("input_root", type=Path, help="RecoilJets ROOT output containing the Fig. 8 histogram")
    ap.add_argument("--trigger", default=None, help="Trigger/directory key, for example SIM. Auto-detects if omitted.")
    ap.add_argument("--hist", default=DEFAULT_HIST, help="Histogram name inside the trigger directory")
    ap.add_argument("--audit-hist", default=DEFAULT_AUDIT, help="Audit-flow histogram name")
    ap.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    ap.add_argument("--output-name", default=None)
    ap.add_argument("--cluster-label", default="Cluster source not specified")
    ap.add_argument("--title", default="Pythia signal MC")
    ap.add_argument("--logz", action="store_true", help="Use logarithmic color scale")
    ap.add_argument("--rebin-x", type=int, default=1, help="Display-only truth-ET rebin factor")
    ap.add_argument("--rebin-y", type=int, default=1, help="Display-only response rebin factor")
    ap.add_argument(
        "--shape-normalize-max",
        type=float,
        default=None,
        help="Globally scale the visible histogram maximum to this value (shape-only display)",
    )
    ap.add_argument("--vmax", type=float, default=None, help="Fixed color-scale maximum")
    ap.add_argument("--cmap", default="viridis", help="Matplotlib color map")
    ap.add_argument("--xlim", nargs=2, type=float, default=(10.0, 30.0))
    ap.add_argument("--ylim", nargs=2, type=float, default=(0.0, 2.0))
    ap.add_argument("--dpi", type=int, default=180)
    return ap.parse_args()


def strip_cycle(key: str) -> str:
    return key.split(";", 1)[0]


def find_hist_path(root_file: uproot.ReadOnlyDirectory, trigger: str | None, hist_name: str) -> str:
    if trigger:
        path = f"{trigger}/{hist_name}"
        if path in root_file:
            return path
        raise KeyError(f"Missing histogram {path}")

    matches = []
    for key in root_file.keys(recursive=True):
        clean = strip_cycle(key)
        if clean.endswith("/" + hist_name) or clean == hist_name:
            matches.append(clean)
    if not matches:
        raise KeyError(f"Missing histogram named {hist_name}")
    preferred = [m for m in matches if m.startswith("SIM/")]
    return (preferred or matches)[0]


def read_th2(root_file: uproot.ReadOnlyDirectory, path: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    obj = root_file[path]
    values, x_edges, y_edges = obj.to_numpy(flow=False)
    return np.asarray(values, dtype=float), np.asarray(x_edges, dtype=float), np.asarray(y_edges, dtype=float)


def read_audit(root_file: uproot.ReadOnlyDirectory, hist_path: str, audit_name: str) -> dict[str, Any]:
    trigger = hist_path.rsplit("/", 1)[0] if "/" in hist_path else ""
    audit_path = f"{trigger}/{audit_name}" if trigger else audit_name
    if audit_path not in root_file:
        return {"path": audit_path, "present": False}

    obj = root_file[audit_path]
    values, edges = obj.to_numpy(flow=False)
    axis = obj.axis()
    try:
        labels = list(axis.labels())
    except Exception:
        labels = [f"bin_{idx + 1}" for idx in range(len(values))]
    labels = [label or f"bin_{idx + 1}" for idx, label in enumerate(labels)]
    return {
        "path": audit_path,
        "present": True,
        "edges": [float(x) for x in edges],
        "counts": {labels[i]: float(values[i]) for i in range(len(values))},
    }


def rebin_2d(
    values: np.ndarray,
    x_edges: np.ndarray,
    y_edges: np.ndarray,
    factor_x: int,
    factor_y: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if factor_x < 1 or factor_y < 1:
        raise ValueError("Rebin factors must be positive integers")
    nx, ny = values.shape
    if nx % factor_x or ny % factor_y:
        raise ValueError(
            f"Histogram shape {values.shape} is not divisible by rebin factors "
            f"({factor_x}, {factor_y})"
        )
    rebinned = values.reshape(nx // factor_x, factor_x, ny // factor_y, factor_y).sum(axis=(1, 3))
    return rebinned, x_edges[::factor_x], y_edges[::factor_y]


def visible_maximum(
    values: np.ndarray,
    x_edges: np.ndarray,
    y_edges: np.ndarray,
    xlim: tuple[float, float],
    ylim: tuple[float, float],
) -> float:
    xmask = (x_edges[:-1] < xlim[1]) & (x_edges[1:] > xlim[0])
    ymask = (y_edges[:-1] < ylim[1]) & (y_edges[1:] > ylim[0])
    visible = values[np.ix_(xmask, ymask)]
    return float(np.max(visible)) if visible.size else 0.0


def render(
    values: np.ndarray,
    x_edges: np.ndarray,
    y_edges: np.ndarray,
    output_png: Path,
    *,
    title: str,
    cluster_label: str,
    xlim: tuple[float, float],
    ylim: tuple[float, float],
    logz: bool,
    vmax: float | None,
    cmap_name: str,
    colorbar_label: str,
    dpi: int,
) -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "axes.linewidth": 1.2,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
    })

    fig, ax = plt.subplots(figsize=(6.4, 5.0), dpi=dpi)
    masked = np.ma.masked_where(values.T <= 0.0, values.T)
    data_max = float(np.nanmax(masked)) if masked.count() else 1.0
    if logz:
        norm = colors.LogNorm(
            vmin=max(1.0, float(np.nanmin(masked)) if masked.count() else 1.0),
            vmax=vmax or data_max,
        )
    else:
        norm = colors.Normalize(vmin=0.0, vmax=vmax or data_max)
    cmap = plt.get_cmap(cmap_name).copy()
    cmap.set_bad("white")
    mesh = ax.pcolormesh(x_edges, y_edges, masked, shading="auto", cmap=cmap, norm=norm)

    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xlabel(r"$E_T^{\gamma,\mathrm{truth}}$ [GeV]", fontsize=12)
    ax.set_ylabel(r"$E_T^{\mathrm{cluster}} / E_T^{\gamma,\mathrm{truth}}$", fontsize=12)
    ax.text(0.05, 0.95, r"$\bf{\it{sPHENIX}}$ Internal", transform=ax.transAxes,
            ha="left", va="top", fontsize=12)
    ax.text(0.05, 0.88, title, transform=ax.transAxes, ha="left", va="top", fontsize=10.5)
    ax.text(0.05, 0.81, cluster_label, transform=ax.transAxes, ha="left", va="top", fontsize=10.5)
    ax.text(0.96, 0.95, r"$|z_{\mathrm{vertex}}| < 30$ cm",
            transform=ax.transAxes, ha="right", va="top", fontsize=10.5)
    ax.text(0.96, 0.88, r"$|\eta^\gamma_{\mathrm{truth}}| < 0.7$",
            transform=ax.transAxes, ha="right", va="top", fontsize=10.5)

    from matplotlib.ticker import AutoMinorLocator, MultipleLocator
    ax.xaxis.set_major_locator(MultipleLocator(2.0))
    ax.yaxis.set_major_locator(MultipleLocator(0.2))
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))

    cbar = fig.colorbar(mesh, ax=ax, pad=0.02)
    cbar.set_label(colorbar_label, fontsize=10.5)
    cbar.ax.tick_params(labelsize=9.5)
    ax.tick_params(labelsize=10)

    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(output_png)
    plt.close(fig)


def main() -> int:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    output_name = args.output_name or f"{args.hist}.png"
    output_png = args.outdir / output_name
    manifest = args.outdir / (Path(output_name).stem + ".manifest.json")

    with uproot.open(args.input_root) as root_file:
        hist_path = find_hist_path(root_file, args.trigger, args.hist)
        values, x_edges, y_edges = read_th2(root_file, hist_path)
        audit = read_audit(root_file, hist_path, args.audit_hist)

    raw_shape = list(values.shape)
    raw_visible_max = visible_maximum(values, x_edges, y_edges, tuple(args.xlim), tuple(args.ylim))
    values, x_edges, y_edges = rebin_2d(values, x_edges, y_edges, args.rebin_x, args.rebin_y)
    rebinned_visible_max = visible_maximum(values, x_edges, y_edges, tuple(args.xlim), tuple(args.ylim))
    display_scale = 1.0
    if args.shape_normalize_max is not None:
        if args.shape_normalize_max <= 0 or rebinned_visible_max <= 0:
            raise ValueError("Shape-normalization target and visible histogram maximum must be positive")
        display_scale = args.shape_normalize_max / rebinned_visible_max
        values = values * display_scale
    display_vmax = args.vmax or args.shape_normalize_max
    colorbar_label = (
        "Display-normalized entries / bin"
        if args.shape_normalize_max is not None
        else "Weighted entries / bin"
    )

    render(
        values,
        x_edges,
        y_edges,
        output_png,
        title=args.title,
        cluster_label=args.cluster_label,
        xlim=tuple(args.xlim),
        ylim=tuple(args.ylim),
        logz=args.logz,
        vmax=display_vmax,
        cmap_name=args.cmap,
        colorbar_label=colorbar_label,
        dpi=args.dpi,
    )

    manifest.write_text(
        json.dumps(
            {
                "input_root": str(args.input_root),
                "hist_path": hist_path,
                "output_png": str(output_png),
                "plot_definition": {
                    "x": "truth photon E_T [GeV]",
                    "y": "cluster E_T / truth photon E_T",
                    "selection": "pp SIM, PPG12 truth signal photon, |zvertex| < 30 cm, |eta_truth| < 0.7, reco cluster E_T >= 5 GeV",
                    "cluster_target": args.cluster_label,
                    "ian_reference": "PPG12 current IAN Fig. 8; source macro compare_energy_response.C books h_respET_* as clusterET / truth pT vs truth pT",
                },
                "display_transform": {
                    "raw_histogram_shape": raw_shape,
                    "rebin_x": args.rebin_x,
                    "rebin_y": args.rebin_y,
                    "raw_visible_maximum": raw_visible_max,
                    "rebinned_visible_maximum": rebinned_visible_max,
                    "shape_normalize_max": args.shape_normalize_max,
                    "global_display_scale": display_scale,
                    "color_scale_maximum": display_vmax,
                    "color_map": args.cmap,
                    "log_color_scale": args.logz,
                    "normalization_scope": "one global multiplicative display scale; no per-bin or per-slice tuning",
                },
                "audit": audit,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    print(output_png)
    print(manifest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
