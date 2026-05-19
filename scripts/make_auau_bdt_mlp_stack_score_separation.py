#!/usr/bin/env python3
"""Make a 3x3 score-separation panel aligned to the exact stack ROC inputs."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import make_auau_bdt_mlp_stack_roc_overlay as rocutil  # noqa: E402
import train_auau_stacked_bdt_mlp_sweep as stack  # noqa: E402


def write_json(path: Path, payload) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(rocutil.json_ready(payload), indent=2, sort_keys=True) + "\n")


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mlp-cache", type=Path, required=True)
    ap.add_argument("--bdt-cache", type=Path, required=True)
    ap.add_argument("--artifact", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--mlp-score", default=stack.DEFAULT_MLP_SCORE)
    ap.add_argument("--bdt-score", default=stack.DEFAULT_BDT_SCORE)
    ap.add_argument("--pt-min", type=float, default=15.0)
    ap.add_argument("--pt-max", type=float, default=35.0)
    ap.add_argument("--cent-bins", default="0,20,50,80")
    ap.add_argument("--split", choices=["val", "test", "all"], default="test")
    ap.add_argument("--train-fraction", type=float, default=0.60)
    ap.add_argument("--val-fraction", type=float, default=0.20)
    ap.add_argument("--random-seed", type=int, default=24681357)
    ap.add_argument("--max-shards", type=int, default=0)
    ap.add_argument("--max-rows", type=int, default=0)
    ap.add_argument("--bins", type=int, default=34)
    return ap.parse_args()


def density_hist(ax, signal, background, bins, color_sig, color_bkg):
    edges = np.linspace(0.0, 1.0, bins + 1)
    ax.hist(signal, bins=edges, density=True, histtype="step", lw=1.8, color=color_sig, label="Signal")
    ax.hist(background, bins=edges, density=True, histtype="step", lw=1.8, color=color_bkg, label="Background")


def main() -> None:
    args = parse_args()
    if args.pt_min >= args.pt_max:
        raise SystemExit("--pt-min must be smaller than --pt-max")
    args.outdir.mkdir(parents=True, exist_ok=True)
    artifact = json.loads(args.artifact.read_text())
    load_args = argparse.Namespace(
        mlp_cache=args.mlp_cache,
        bdt_cache=args.bdt_cache,
        logreg_cache=None,
        mlp_score=args.mlp_score,
        bdt_score=args.bdt_score,
        logreg_score="",
        max_shards=args.max_shards,
        max_rows=args.max_rows,
        require_full_stat=False,
        expected_shards=0,
        random_seed=args.random_seed,
    )
    frame, cache_info = stack.load_aligned_caches(load_args)
    y = np.asarray(frame["is_signal"], dtype="int32")
    masks = stack.stratified_split(y, args.random_seed, args.train_fraction, args.val_fraction)
    split_mask = np.ones(len(y), dtype=bool) if args.split == "all" else masks[args.split]
    stack_score = rocutil.route_score(artifact, frame)

    et = np.asarray(frame["cluster_Et"], dtype="float64")
    cent = np.asarray(frame["centrality"], dtype="float64")
    cent_edges = rocutil.parse_edges(args.cent_bins)
    base = split_mask & np.isfinite(et) & np.isfinite(cent) & (et >= args.pt_min) & (et < args.pt_max)
    rows = [
        ("Input BDT\n8 $E_{T}$ x 7 cent", "Input BDT, 8 E_T x 7 centrality", np.asarray(frame["bdt_score"], dtype="float64")),
        ("Input MLP\nsingle small NN", "Input MLP, single small NN", np.asarray(frame["mlp_score"], dtype="float64")),
        ("BDT+MLP stack\nNN combiner", "BDT+MLP stack, NN combiner", stack_score),
    ]

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "axes.linewidth": 1.0,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )
    fig, axes = plt.subplots(3, len(cent_edges) - 1, figsize=(15.2, 8.9), sharex=True, sharey=True)
    summary: list[dict] = []
    ymax_global = 0.0
    hist_payload: list[dict] = []
    for irow, (row_label, model_label, score) in enumerate(rows):
        for icol, (clo, chi) in enumerate(zip(cent_edges[:-1], cent_edges[1:], strict=True)):
            mask = base & (cent >= clo) & (cent < chi) & np.isfinite(score)
            sig = score[mask & (y == 1)]
            bkg = score[mask & (y == 0)]
            _, _, auc = rocutil.roc_curve_auc(y[mask], score[mask])
            edges = np.linspace(0.0, 1.0, args.bins + 1)
            sig_density = np.histogram(sig, bins=edges, density=True)[0] if sig.size else np.zeros(args.bins)
            bkg_density = np.histogram(bkg, bins=edges, density=True)[0] if bkg.size else np.zeros(args.bins)
            ymax_global = max(ymax_global, float(np.nanmax(sig_density)), float(np.nanmax(bkg_density)))
            hist_payload.append(
                {
                    "irow": irow,
                    "icol": icol,
                    "row_label": row_label,
                    "model_label": model_label,
                    "cent_label": f"{clo:g}-{chi:g}% centrality",
                    "cent_lo": clo,
                    "cent_hi": chi,
                    "sig": sig,
                    "bkg": bkg,
                    "auc": auc,
                    "entries": int(mask.sum()),
                    "signal_entries": int(sig.size),
                    "background_entries": int(bkg.size),
                }
            )
    ymax_global = ymax_global * 1.18 if ymax_global > 0 else 1.0

    for payload in hist_payload:
            irow = payload["irow"]
            icol = payload["icol"]
            ax = axes[irow, icol]
            density_hist(ax, payload["sig"], payload["bkg"], args.bins, "#1769aa", "#c43c4e")
            ax.set_ylim(0.0, ymax_global)
            ax.text(
                0.045,
                0.86,
                f"AUC = {payload['auc']:.3f}",
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=11.0,
                fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.22", fc="white", ec="0.78", alpha=0.92),
            )
            if irow == 0:
                ax.set_title(payload["cent_label"], fontsize=13.2, fontweight="bold")
            if icol == 0:
                ax.set_ylabel("Area-normalized density", fontsize=11.0)
                ax.text(
                    -0.35,
                    0.5,
                    payload["row_label"],
                    transform=ax.transAxes,
                    rotation=90,
                    ha="center",
                    va="center",
                    fontsize=12.0,
                    fontweight="bold",
                )
            if irow == 2:
                ax.set_xlabel("Classifier score", fontsize=10.5)
            if irow == 0 and icol == len(cent_edges) - 2:
                ax.legend(loc="upper right", fontsize=9.5, frameon=False)
            ax.grid(True, color="0.90", lw=0.5)
            ax.set_xlim(0, 1)
            summary.append(
                {
                    "model": payload["model_label"],
                    "cent_lo": payload["cent_lo"],
                    "cent_hi": payload["cent_hi"],
                    "split": args.split,
                    "entries": payload["entries"],
                    "signal_entries": payload["signal_entries"],
                    "background_entries": payload["background_entries"],
                    "auc": payload["auc"],
                }
            )

    fig.text(0.035, 0.975, r"$\bf{\it{sPHENIX}}$ Internal", ha="left", va="top", fontsize=14)
    fig.text(0.035, 0.944, "Au+Au embedded validation", ha="left", va="top", fontsize=10.5)
    fig.text(0.50, 0.975, r"Signal/background score separation, $15 < E_{T} < 35$ GeV", ha="center", va="top", fontsize=17, fontweight="bold")
    split_label = {"val": "held-out validation split", "test": "held-out test split", "all": "all scored rows"}[args.split]
    fig.text(0.50, 0.944, f"{split_label}; same exact stack artifact/cache inputs as ROC overlay", ha="center", va="top", fontsize=10.3, color="0.35")
    fig.text(0.50, 0.017, "BDT+MLP stack uses the exported JSON artifact evaluated on aligned enriched validation caches.", ha="center", va="bottom", fontsize=8.6, color="0.35")
    fig.tight_layout(rect=[0.045, 0.055, 0.995, 0.91], h_pad=1.0, w_pad=1.8)

    pt_tag = f"{args.pt_min:g}to{args.pt_max:g}".replace(".", "p")
    png = args.outdir / f"score_separation_3x3_exact_stack_pt{pt_tag}_{args.split}.png"
    fig.savefig(png, dpi=220)
    plt.close(fig)

    csv_path = args.outdir / "score_separation_3x3_summary.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary[0]))
        writer.writeheader()
        writer.writerows(summary)
    metadata = {
        "schema": "RECOILJETS_AUAU_BDT_MLP_STACK_SCORE_SEPARATION_V1",
        "plot": str(png),
        "summary_csv": str(csv_path),
        "artifact": str(args.artifact),
        "artifact_name": artifact.get("name", ""),
        "mlp_cache": str(args.mlp_cache),
        "bdt_cache": str(args.bdt_cache),
        "mlp_score": cache_info["mlp_score_key"],
        "bdt_score": cache_info["bdt_score_key"],
        "rows_loaded": cache_info["rows"],
        "selection": f"{args.pt_min:g} < cluster_Et < {args.pt_max:g} GeV",
        "centrality_bins": cent_edges,
        "split": args.split,
        "note": "This panel is aligned to the corrected ROC overlay; do not compare with older all-row/sweep-artifact panels.",
    }
    write_json(args.outdir / "score_separation_3x3_metadata.json", metadata)
    print(f"[stackSep] plot={png}", flush=True)
    print(f"[stackSep] summary={csv_path}", flush=True)
    print(f"[stackSep] metadata={args.outdir / 'score_separation_3x3_metadata.json'}", flush=True)


if __name__ == "__main__":
    main()
