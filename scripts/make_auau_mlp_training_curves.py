#!/usr/bin/env python3
"""Plot AuAu MLP training and validation loss versus epoch from real histories."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path


LOG_RE = re.compile(
    r"candidate=(?P<label>\S+)\s+epoch=(?P<epoch>\d+)\s+"
    r"train_loss=(?P<train>[-+0-9.eE]+)\s+val_loss=(?P<val>[-+0-9.eE]+)"
    r"(?:\s+truth_val_loss=(?P<truth>[-+0-9.eE]+))?.*?"
    r"(?:val_auc=(?P<auc>[-+0-9.eE]+))?"
)


def parse_label_path(text: str) -> tuple[str, Path]:
    if "=" in text:
        label, path = text.split("=", 1)
        return label.strip(), Path(path.strip())
    path = Path(text)
    return path.stem, path


def as_float(value) -> float:
    try:
        out = float(value)
    except Exception:
        return math.nan
    return out if math.isfinite(out) else math.nan


def read_history_csv(path: Path, label: str) -> list[dict]:
    rows = []
    with path.open() as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            epoch = int(float(row.get("epoch", len(rows) + 1)))
            rows.append(
                {
                    "label": label,
                    "epoch": epoch,
                    "train_loss": as_float(row.get("train_loss")),
                    "validation_loss": as_float(row.get("validation_loss", row.get("val_loss"))),
                    "truth_validation_loss": as_float(row.get("truth_validation_loss", row.get("truth_val_loss"))),
                    "validation_auc": as_float(row.get("validation_auc", row.get("val_auc"))),
                    "source": str(path),
                }
            )
    return rows


def read_log(path: Path, wanted_label: str | None = None) -> list[dict]:
    rows = []
    with path.open(errors="replace") as handle:
        for line in handle:
            match = LOG_RE.search(line)
            if not match:
                continue
            label = match.group("label")
            if wanted_label and wanted_label not in label:
                continue
            rows.append(
                {
                    "label": wanted_label or label,
                    "epoch": int(match.group("epoch")),
                    "train_loss": as_float(match.group("train")),
                    "validation_loss": as_float(match.group("val")),
                    "truth_validation_loss": as_float(match.group("truth")),
                    "validation_auc": as_float(match.group("auc")),
                    "source": str(path),
                }
            )
    return rows


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = ["label", "epoch", "train_loss", "validation_loss", "truth_validation_loss", "validation_auc", "source"]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, payload) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--history", action="append", default=[], help="label=/path/to/history.csv or /path")
    ap.add_argument("--log", action="append", default=[], help="label=/path/to/log or /path")
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--title", default="MLP training curves")
    args = ap.parse_args()

    all_rows: list[dict] = []
    for item in args.history:
        label, path = parse_label_path(item)
        all_rows.extend(read_history_csv(path, label))
    for item in args.log:
        label, path = parse_label_path(item)
        all_rows.extend(read_log(path, label))
    if not all_rows:
        raise SystemExit("No training-history rows found")

    args.outdir.mkdir(parents=True, exist_ok=True)
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

    labels = []
    for row in all_rows:
        if row["label"] not in labels:
            labels.append(row["label"])
    fig, axes = plt.subplots(len(labels), 1, figsize=(8.2, max(3.6, 3.1 * len(labels))), sharex=False)
    if len(labels) == 1:
        axes = [axes]
    colors = {"train": "#1f77b4", "val": "#d95f02", "truth": "#7570b3"}
    summary = []
    for ax, label in zip(axes, labels, strict=True):
        rows = sorted([r for r in all_rows if r["label"] == label], key=lambda r: r["epoch"])
        epochs = [r["epoch"] for r in rows]
        train = [r["train_loss"] for r in rows]
        val = [r["validation_loss"] for r in rows]
        truth = [r["truth_validation_loss"] for r in rows]
        ax.plot(epochs, train, color=colors["train"], lw=1.9, marker="o", ms=2.7, label="Training loss")
        ax.plot(epochs, val, color=colors["val"], lw=1.9, marker="s", ms=2.7, label="Validation loss")
        if any(math.isfinite(x) for x in truth):
            ax.plot(epochs, truth, color=colors["truth"], lw=1.5, ls="--", label="Truth-only validation loss")
        best = min((r for r in rows if math.isfinite(r["validation_loss"])), key=lambda r: r["validation_loss"], default=None)
        if best:
            ax.axvline(best["epoch"], color="0.35", lw=1.0, ls=":")
            ax.text(0.985, 0.90, f"best val epoch {best['epoch']}", transform=ax.transAxes, ha="right", va="top", fontsize=9.0, color="0.25")
        ax.set_title(label, fontsize=12.5, fontweight="bold")
        ax.set_ylabel("Loss", fontsize=11)
        ax.grid(True, color="0.90", lw=0.6)
        ax.legend(loc="best", fontsize=9, frameon=False)
        ax.set_xlabel("Epoch", fontsize=11)
        summary.append(
            {
                "label": label,
                "epochs": len(rows),
                "last_epoch": rows[-1]["epoch"],
                "best_validation_epoch": None if best is None else best["epoch"],
                "best_validation_loss": None if best is None else best["validation_loss"],
                "last_train_loss": rows[-1]["train_loss"],
                "last_validation_loss": rows[-1]["validation_loss"],
                "source": rows[0]["source"],
            }
        )
    fig.text(0.02, 0.985, r"$\bf{\it{sPHENIX}}$ Internal", ha="left", va="top", fontsize=11.5)
    fig.suptitle(args.title, fontsize=13.5, fontweight="bold", x=0.60, y=0.985)
    fig.tight_layout(rect=[0.02, 0.02, 1.0, 0.925], h_pad=1.25)
    png = args.outdir / "mlp_training_validation_loss_vs_epoch.png"
    fig.savefig(png, dpi=220)
    plt.close(fig)

    csv_path = args.outdir / "mlp_training_history_combined.csv"
    write_csv(csv_path, all_rows)
    write_json(args.outdir / "mlp_training_curve_summary.json", {"plot": str(png), "history_csv": str(csv_path), "models": summary})
    print(f"[mlpCurves] plot={png}", flush=True)
    print(f"[mlpCurves] history={csv_path}", flush=True)
    print(f"[mlpCurves] summary={args.outdir / 'mlp_training_curve_summary.json'}", flush=True)


if __name__ == "__main__":
    main()
