#!/usr/bin/env python3
"""Plot signal-only BDT score split by current AuAu isolation decision."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


CENT3_BINS = [(0.0, 20.0, "0-20%"), (20.0, 50.0, "20-50%"), (50.0, 80.0, "50-80%")]
FINE_CENT_EDGES = np.array([0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0])
R30_WP90_GEV = np.array([5.7165, 5.2095, 4.7025, 4.1955, 3.6885, 3.1815, 2.4210])


def score_cache_paths(path: Path) -> list[Path]:
    if path.is_dir():
        manifest = path / "score_caches.local.list"
        if not manifest.is_file():
            manifest = path / "score_caches.list"
    else:
        manifest = path
    if manifest.suffix == ".npz" and manifest.is_file():
        return [manifest]
    if not manifest.is_file():
        raise SystemExit(f"Missing score-cache manifest: {manifest}")
    paths = [Path(line.strip()) for line in manifest.read_text().splitlines() if line.strip()]
    if not paths:
        raise SystemExit(f"No score-cache paths in {manifest}")
    return paths


def load_arrays(manifest: Path, score_key: str) -> dict[str, np.ndarray]:
    columns = ["is_signal", "cluster_Et", "centrality", "reco_eiso", score_key]
    pieces = {key: [] for key in columns}
    for path in score_cache_paths(manifest):
        with np.load(path, allow_pickle=True) as z:
            missing = [key for key in columns if key not in z.files]
            if missing:
                raise SystemExit(f"{path} is missing {missing}; available={z.files}")
            for key in columns:
                pieces[key].append(np.asarray(z[key]))
    return {key: np.concatenate(value) for key, value in pieces.items()}


def r30_cent_dependent_wp(centrality: np.ndarray) -> np.ndarray:
    idx = np.searchsorted(FINE_CENT_EDGES, centrality, side="right") - 1
    idx = np.clip(idx, 0, len(R30_WP90_GEV) - 1)
    return R30_WP90_GEV[idx]


def unit_density(values: np.ndarray, bins: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    hist, edges = np.histogram(values, bins=bins)
    width = np.diff(edges)
    total = hist.sum()
    if total <= 0:
        return np.full(len(width), np.nan), 0.5 * (edges[:-1] + edges[1:])
    return hist / (total * width), 0.5 * (edges[:-1] + edges[1:])


def add_sphenix_label(fig: plt.Figure) -> None:
    fig.text(0.046, 0.835, "sPHENIX", ha="left", va="center", fontsize=18, fontstyle="italic", fontweight="bold")
    fig.text(0.122, 0.835, " Internal", ha="left", va="center", fontsize=18)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", type=Path, required=True)
    ap.add_argument("--score-key", default="score_globalEtCent1535_bdt_iso_ptCent7")
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--tag", default="binned_bdt_iso_ptcent7_signal_iso_split")
    ap.add_argument("--model-label", default="Binned BDT, 8 E_T x 7 centrality + isolation inputs")
    args = ap.parse_args()

    data = load_arrays(args.manifest, args.score_key)
    y = data["is_signal"].astype(np.int8)
    et = data["cluster_Et"].astype(float)
    cent = data["centrality"].astype(float)
    eiso = data["reco_eiso"].astype(float)
    score = data[args.score_key].astype(float)

    valid = (
        (y == 1)
        & np.isfinite(et)
        & np.isfinite(cent)
        & np.isfinite(eiso)
        & np.isfinite(score)
        & (et >= 15.0)
        & (et < 35.0)
        & (cent >= 0.0)
        & (cent < 80.0)
    )
    wp = r30_cent_dependent_wp(cent)
    iso_pass = eiso < wp

    bins = np.linspace(0.0, 1.0, 51)
    colors = {"pass": "#0072B2", "fail": "#D55E00"}
    fig, axes = plt.subplots(1, 3, figsize=(17.2, 6.0), sharex=True, sharey=False, dpi=185)
    fig.patch.set_facecolor("white")
    rows: list[dict] = []

    for ax, (clo, chi, label) in zip(axes, CENT3_BINS):
        cmask = valid & (cent >= clo) & (cent < chi)
        pass_mask = cmask & iso_pass
        fail_mask = cmask & (~iso_pass)
        max_density = 0.0
        for mask, split, style, color in [
            (pass_mask, "passes R=0.3 cent-dependent isolation", "-", colors["pass"]),
            (fail_mask, "fails R=0.3 cent-dependent isolation", "--", colors["fail"]),
        ]:
            values = score[mask]
            density, centers = unit_density(values, bins)
            if np.isfinite(density).any():
                max_density = max(max_density, float(np.nanmax(density)))
            ax.step(centers, density, where="mid", color=color, linestyle=style, linewidth=2.8, label=split)
            rows.append(
                {
                    "centrality_bin": label,
                    "split": split,
                    "entries": int(values.size),
                    "fraction_of_signal_bin": float(values.size / max(1, cmask.sum())),
                    "mean_score": float(np.mean(values)) if values.size else float("nan"),
                    "median_score": float(np.median(values)) if values.size else float("nan"),
                    "mean_reco_eiso": float(np.mean(eiso[mask])) if values.size else float("nan"),
                    "mean_wp": float(np.mean(wp[mask])) if values.size else float("nan"),
                }
            )

        ax.set_title(f"{label} centrality", fontsize=16, fontweight="bold")
        ax.set_xlim(0.0, 1.0)
        ax.set_ylim(0.0, max(2.2, 1.22 * max_density))
        ax.grid(True, color="#D1D5DB", linewidth=1.0, alpha=0.75)
        ax.tick_params(direction="in", top=True, right=True, labelsize=12)
        ax.minorticks_on()
        ax.set_xlabel("BDT score", fontsize=14)
        total = int(cmask.sum())
        pass_frac = pass_mask.sum() / max(1, total)
        ax.text(
            0.04,
            0.93,
            f"signal entries: {total:,}\niso-pass fraction: {pass_frac:.3f}",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=10.8,
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.84, "pad": 3},
        )
    axes[0].set_ylabel("Unit-normalized signal candidates", fontsize=14)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        bbox_to_anchor=(0.52, 0.025),
        ncol=2,
        frameon=False,
        fontsize=11.6,
        handlelength=2.6,
        columnspacing=2.2,
    )

    fig.text(0.045, 0.94, "Signal BDT-score split by isolation decision", ha="left", va="top", fontsize=25, fontweight="bold")
    fig.text(
        0.045,
        0.898,
        f"{args.model_label}; embedded Photon12+20 signal, 15 < cluster $E_T$ < 35 GeV",
        ha="left",
        va="top",
        fontsize=13.8,
        color="#374151",
    )
    add_sphenix_label(fig)
    fig.tight_layout(rect=[0.04, 0.135, 0.995, 0.80])

    args.outdir.mkdir(parents=True, exist_ok=True)
    png = args.outdir / f"{args.tag}.png"
    csv_path = args.outdir / f"{args.tag}.csv"
    json_path = args.outdir / f"{args.tag}.json"
    fig.savefig(png)
    plt.close(fig)

    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    json_path.write_text(
        json.dumps(
            {
                "schema": "AUAU_SIGNAL_ISO_SPLIT_BDT_SCORE_V1",
                "manifest": str(args.manifest),
                "score_key": args.score_key,
                "model_label": args.model_label,
                "sample": "embedded Photon12+20 signal only",
                "pt_range_GeV": [15.0, 35.0],
                "centrality_bins": CENT3_BINS,
                "isolation_definition": {
                    "cone": "R=0.3",
                    "type": "current AuAu fine-centrality-dependent WP90 reco isolation",
                    "fine_centrality_edges": FINE_CENT_EDGES.tolist(),
                    "wp_GeV": R30_WP90_GEV.tolist(),
                    "pass_rule": "reco_eiso < wp_GeV(centrality)",
                },
                "rows": rows,
                "outputs": {"png": str(png), "csv": str(csv_path), "json": str(json_path)},
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    print(png)
    print(csv_path)
    print(json_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
