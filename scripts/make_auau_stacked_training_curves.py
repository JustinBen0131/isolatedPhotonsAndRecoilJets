#!/usr/bin/env python3
"""Plot training histories from the AuAu BDT+MLP stacker sweep."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path


def read_csv(path: Path) -> list[dict]:
    if not path.is_file():
        raise SystemExit(f"Missing CSV: {path}")
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def to_float(value, default=math.nan) -> float:
    try:
        out = float(value)
    except Exception:
        return default
    return out if math.isfinite(out) else default


def progress_value(row: dict) -> tuple[str, float]:
    for key, label in (("epoch", "epoch"), ("iteration", "boosting iteration"), ("step", "optimizer step")):
        if key in row and str(row[key]).strip():
            return label, to_float(row[key])
    return "training step", math.nan


def ranked_models(rank_table: Path | None, top_n: int) -> list[str]:
    if rank_table is None or not rank_table.is_file():
        return []
    rows = [row for row in read_csv(rank_table) if row.get("split") == "test"]
    rows.sort(
        key=lambda row: (
            to_float(row.get("highpt_wp80_fake_20_35"), 999.0),
            -to_float(row.get("highpt_auc_20_35"), -999.0),
            -to_float(row.get("auc"), -999.0),
        )
    )
    out = []
    for row in rows:
        model = row.get("model", "")
        if model and model not in out:
            out.append(model)
        if len(out) >= top_n:
            break
    return out


def aggregate_history(rows: list[dict]) -> dict[str, list[dict]]:
    grouped = defaultdict(list)
    for row in rows:
        model = row.get("model", "")
        label, progress = progress_value(row)
        if not model or not math.isfinite(progress):
            continue
        key = (model, label, progress)
        grouped[key].append(row)

    by_model = defaultdict(list)
    for (model, label, progress), items in grouped.items():
        out = {"model": model, "progress_label": label, "progress": progress, "routes": len(items)}
        for field in ("train_loss", "val_loss", "train_auc", "val_auc"):
            vals = [to_float(item.get(field)) for item in items]
            vals = [val for val in vals if math.isfinite(val)]
            out[field] = sum(vals) / len(vals) if vals else math.nan
        by_model[model].append(out)
    for model in by_model:
        by_model[model].sort(key=lambda row: row["progress"])
    return by_model


def write_summary(path: Path, models: list[str], by_model: dict[str, list[dict]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = ["model", "progress_label", "progress", "routes", "train_loss", "val_loss", "train_auc", "val_auc"]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for model in models:
            for row in by_model.get(model, []):
                writer.writerow(row)


def plot_model_curves(outdir: Path, models: list[str], by_model: dict[str, list[dict]]) -> list[Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    outdir.mkdir(parents=True, exist_ok=True)
    outputs = []
    for model in models:
        rows = by_model.get(model, [])
        if not rows:
            continue
        x = [row["progress"] for row in rows]
        xlabel = rows[0]["progress_label"]
        fig, axes = plt.subplots(1, 2, figsize=(12.2, 4.6), dpi=170)
        ax = axes[0]
        ax.plot(x, [row["train_loss"] for row in rows], color="#1f77b4", lw=2.0, label="train loss")
        ax.plot(x, [row["val_loss"] for row in rows], color="#d62728", lw=2.0, label="validation loss")
        ax.set_xlabel(xlabel)
        ax.set_ylabel("binary cross-entropy")
        ax.grid(alpha=0.28)
        ax.legend(frameon=False, fontsize=9)

        ax = axes[1]
        ax.plot(x, [row["train_auc"] for row in rows], color="#1f77b4", lw=2.0, label="train AUC")
        ax.plot(x, [row["val_auc"] for row in rows], color="#d62728", lw=2.0, label="validation AUC")
        ax.set_xlabel(xlabel)
        ax.set_ylabel("AUC")
        ax.set_ylim(0.45, 1.01)
        ax.grid(alpha=0.28)
        ax.legend(frameon=False, fontsize=9)

        fig.text(0.025, 0.94, "sPHENIX", fontsize=14, fontstyle="italic", fontweight="bold", ha="left")
        fig.text(0.122, 0.94, "Internal", fontsize=14, ha="left")
        fig.suptitle(model.replace("_", " "), fontsize=14, fontweight="bold", y=0.965)
        fig.text(0.5, 0.02, "For routed models, curves are route-averaged at each training step.", ha="center", color="0.35", fontsize=8.5)
        fig.tight_layout(rect=(0, 0.05, 1, 0.90))
        safe = "".join(ch if ch.isalnum() or ch in ("_", "-") else "_" for ch in model)
        png = outdir / f"{safe}_training_curves.png"
        fig.savefig(png)
        plt.close(fig)
        outputs.append(png)
    return outputs


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--history-csv", type=Path, required=True)
    ap.add_argument("--rank-table", type=Path, default=None)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--models", default="")
    ap.add_argument("--top-n", type=int, default=4)
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    rows = read_csv(args.history_csv)
    requested = [item.strip() for item in args.models.split(",") if item.strip()]
    models = requested or ranked_models(args.rank_table, args.top_n)
    if not models:
        models = []
        for row in rows:
            model = row.get("model", "")
            if model and model not in models:
                models.append(model)
            if len(models) >= args.top_n:
                break
    by_model = aggregate_history(rows)
    summary = args.outdir / "stacked_training_curve_summary.csv"
    write_summary(summary, models, by_model)
    plots = plot_model_curves(args.outdir, models, by_model)
    print(f"[stackTrainingCurves] summary={summary}", flush=True)
    for plot in plots:
        print(f"[stackTrainingCurves] plot={plot}", flush=True)


if __name__ == "__main__":
    main()
