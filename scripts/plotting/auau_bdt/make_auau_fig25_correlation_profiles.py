#!/usr/bin/env python3
"""Render AuAu embedded-background Fig. 25 style correlation profiles."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import uproot


CENTRALITIES = (("0-20%", "cent_0_20"), ("20-50%", "cent_20_50"), ("50-80%", "cent_50_80"))
FAMILIES = {
    "e11e33": ("E11/E33", r"$E_{1\times1}/E_{3\times3}$"),
    "bdtScore": ("Photon-ID BDT score", "Au+Au photon-ID BDT score"),
}
BASE = "h2_auauFig25_{axis}_vs_Eiso_background_pT15to35"


def clean_keys(root_file: uproot.ReadOnlyDirectory) -> list[str]:
    return sorted(str(key).split(";", 1)[0] for key in root_file.keys(recursive=True))


def find_surface(keys: list[str], axis: str, cent_token: str) -> str:
    token = BASE.format(axis=axis)
    matches = [key for key in keys if token in key and cent_token in key]
    if len(matches) != 1:
        raise RuntimeError(f"expected one {axis}/{cent_token} surface, found {matches}")
    return matches[0]


def profile(hist: uproot.behaviors.TH2.Histogram) -> dict[str, np.ndarray | float]:
    weights = np.asarray(hist.values(flow=False), dtype=np.float64)
    sumw2 = np.asarray(hist.variances(flow=False), dtype=np.float64)
    x_edges = np.asarray(hist.axis(0).edges(), dtype=np.float64)
    y_edges = np.asarray(hist.axis(1).edges(), dtype=np.float64)
    x = 0.5 * (x_edges[:-1] + x_edges[1:])
    y = 0.5 * (y_edges[:-1] + y_edges[1:])
    sumw = np.sum(weights, axis=1)
    total_sumw2 = np.sum(sumw2, axis=1)
    weighted_y = np.sum(weights * y[np.newaxis, :], axis=1)
    mean = np.divide(weighted_y, sumw, out=np.full_like(sumw, np.nan), where=sumw > 0)
    variance = np.divide(
        np.sum(weights * (y[np.newaxis, :] - mean[:, np.newaxis]) ** 2, axis=1),
        sumw,
        out=np.full_like(sumw, np.nan),
        where=sumw > 0,
    )
    neff = np.divide(sumw**2, total_sumw2, out=np.zeros_like(sumw), where=total_sumw2 > 0)
    error = np.sqrt(np.divide(variance, neff, out=np.full_like(sumw, np.nan), where=neff > 0))

    xx, yy = np.meshgrid(x, y, indexing="ij")
    total = float(np.sum(weights))
    if total > 0:
        mean_x = float(np.sum(weights * xx) / total)
        mean_y = float(np.sum(weights * yy) / total)
        cov = float(np.sum(weights * (xx - mean_x) * (yy - mean_y)) / total)
        var_x = float(np.sum(weights * (xx - mean_x) ** 2) / total)
        var_y = float(np.sum(weights * (yy - mean_y) ** 2) / total)
        correlation = cov / np.sqrt(var_x * var_y) if var_x > 0 and var_y > 0 else float("nan")
    else:
        correlation = float("nan")
    return {"x": x, "mean": mean, "error": error, "sumw": sumw, "neff": neff,
            "correlation": correlation}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    root_file = uproot.open(args.root)
    keys = clean_keys(root_file)
    csv_rows: list[dict[str, object]] = []
    manifest: dict[str, object] = {
        "schema": "AUAU_EMBEDDED_BACKGROUND_FIG25_PROFILES_V1",
        "input_root": str(args.root.resolve()),
        "sample": "canonical weighted embedded Jet12+20+30+40 background complement",
        "photon_pt_gev": [15.0, 35.0],
        "centrality_percent": [[0, 20], [20, 50], [50, 80]],
        "embedded_minbias_classifier_filter": False,
        "profiles": {},
    }

    for axis, (title_axis, xlabel) in FAMILIES.items():
        fig, panels = plt.subplots(1, 3, figsize=(15.8, 5.2), sharey=True)
        family_manifest: dict[str, object] = {}
        for panel, (cent_label, cent_token) in zip(panels, CENTRALITIES):
            key = find_surface(keys, axis, cent_token)
            stats = profile(root_file[key])
            x = np.asarray(stats["x"])
            mean = np.asarray(stats["mean"])
            error = np.asarray(stats["error"])
            sumw = np.asarray(stats["sumw"])
            neff = np.asarray(stats["neff"])
            mask = np.isfinite(mean) & np.isfinite(error) & (sumw > 0)
            panel.errorbar(x[mask], mean[mask], yerr=error[mask], fmt="o-", color="#2367d1",
                           markersize=3.2, linewidth=1.5, capsize=1.8)
            panel.set_xlim(0.0, 1.0)
            panel.grid(alpha=0.18, linewidth=0.7)
            panel.set_title(cent_label, fontsize=15, fontweight="bold")
            panel.set_xlabel(xlabel, fontsize=13)
            panel.tick_params(labelsize=11)
            panel.text(0.04, 0.05, f"Correlation: {float(stats['correlation']):+.3f}",
                       transform=panel.transAxes, ha="left", va="bottom", fontsize=10.5,
                       bbox={"facecolor": "white", "edgecolor": "#999999", "alpha": 0.92,
                             "boxstyle": "square,pad=0.35", "linewidth": 0.8})
            family_manifest[cent_label] = {"object": key, "correlation": float(stats["correlation"])}
            for i in np.flatnonzero(mask):
                csv_rows.append({
                    "family": axis,
                    "centrality": cent_label,
                    "x_center": float(x[i]),
                    "mean_eiso_gev": float(mean[i]),
                    "mean_eiso_error_gev": float(error[i]),
                    "sum_event_weight": float(sumw[i]),
                    "effective_entries": float(neff[i]),
                })
        panels[0].set_ylabel(r"$\langle E_T^{iso,reco}\rangle$ [GeV]", fontsize=13)
        fig.suptitle("Au+Au Embedded Inclusive-Jet Simulation, 15-35 GeV", fontsize=20,
                     fontweight="bold", y=0.995)
        fig.tight_layout(rect=(0.02, 0.02, 0.99, 0.94))
        output = args.output_dir / f"auau_embedded_background_{axis}_isolation_profiles_cent3.png"
        fig.savefig(output, dpi=180, facecolor="white")
        plt.close(fig)
        manifest["profiles"][axis] = {"png": str(output.resolve()), "centralities": family_manifest}

    csv_path = args.output_dir / "auau_embedded_background_fig25_profiles.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(csv_rows[0]))
        writer.writeheader()
        writer.writerows(csv_rows)
    manifest["csv"] = str(csv_path.resolve())
    json_path = args.output_dir / "auau_embedded_background_fig25_profiles_manifest.json"
    json_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": "PASS", "outputs": manifest["profiles"], "csv": str(csv_path),
                      "manifest": str(json_path)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
