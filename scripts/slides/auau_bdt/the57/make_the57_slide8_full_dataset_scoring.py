#!/usr/bin/env python3
"""Make the THE-57 slide-8 floor-veto check from full-dataset scoring caches.

This intentionally avoids the word "validation" for the source product:
these score caches are an application/scoring pass over the full extracted
sample, not the 90/10 validation split used during training.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch


try:
    REPO = Path(__file__).resolve().parents[4]
except IndexError:
    REPO = Path.cwd()
BASE = REPO / "dataOutput/auauTightBDTFullDatasetScoring/THE57_baseline_full_dataset_scoring_20260615"
OUT = BASE / "slide8_full_dataset_scoring"
PRODUCT = "centAsFeatBase3x3_pt15to35"
SCORE_KEY = f"score_{PRODUCT}"
CENT_KEYS = ["0_20", "20_50", "50_80"]
CENT_LABELS = ["0-20%", "20-50%", "50-80%"]
BROAD_CENT_BINS = [(0.0, 20.0), (20.0, 50.0), (50.0, 80.0)]
FINE_CENT_BINS = [(float(lo), float(lo + 5)) for lo in range(0, 80, 5)]

INK = "#172033"
MUTED = "#5D6B7A"
BLUE = "#1F77B4"
ORANGE = "#C47A1C"
RED = "#D62728"
GREEN = "#1F7A4D"
CARD_BLUE = "#EAF4FB"
CARD_GREEN = "#ECFDF3"
CARD_RED = "#FFF1F2"

plt.rcParams.update(
    {
        "font.family": "Times New Roman",
        "mathtext.fontset": "stix",
        "axes.unicode_minus": False,
    }
)


def pct(x: float, digits: int = 2) -> str:
    return f"{100.0 * x:.{digits}f}%"


def safe_float(value) -> float:
    try:
        return float(value)
    except Exception:
        return math.nan


def read_cache_paths(cache_dir: Path) -> list[Path]:
    paths = sorted(cache_dir.glob("score_cache_*.npz"))
    if not paths:
        raise SystemExit(f"No score_cache_*.npz files under {cache_dir}")
    return paths


def load_scored_rows(cache_dir: Path, *, product: str = PRODUCT) -> dict:
    """Load only the columns needed for the slide from score-cache shards."""

    score_key = f"score_{product}"
    paths = read_cache_paths(cache_dir)
    pieces = {name: [] for name in ["score", "is_signal", "centrality"]}
    counts = {
        "score_cache_files": len(paths),
        "files": 0,
        "total_entries": 0,
        "signal_entries": 0,
        "background_entries": 0,
        "scored_entries": 0,
        "selected_entries": 0,
        "selected_signal_entries": 0,
        "selected_background_entries": 0,
    }
    for idx, path in enumerate(paths, 1):
        if idx % 25 == 0 or idx == 1 or idx == len(paths):
            print(f"[slide8-full] reading {idx}/{len(paths)} {path}", flush=True)
        with np.load(path, allow_pickle=True) as data:
            if score_key not in data.files:
                raise SystemExit(f"Missing {score_key} in {path}")
            score = np.asarray(data[score_key], dtype="float32")
            y = np.asarray(data["is_signal"], dtype="int8")
            cent = np.asarray(data["centrality"], dtype="float32")
            pt = np.asarray(data["cluster_Et"], dtype="float32")
            if "counts_json" in data.files:
                cached_counts = json.loads(str(data["counts_json"].item()))
                for key in ["files", "total_entries", "signal_entries", "background_entries", "scored_entries"]:
                    counts[key] += int(cached_counts.get(key, 0) or 0)
            else:
                counts["total_entries"] += int(len(score))
                counts["scored_entries"] += int(len(score))
                counts["signal_entries"] += int(np.count_nonzero(y == 1))
                counts["background_entries"] += int(np.count_nonzero(y == 0))

            keep = (
                np.isfinite(score)
                & np.isfinite(cent)
                & np.isfinite(pt)
                & (pt >= 15.0)
                & (pt < 35.0)
                & (cent >= 0.0)
                & (cent < 80.0)
                & np.isin(y, [0, 1])
            )
            pieces["score"].append(score[keep])
            pieces["is_signal"].append(y[keep].astype("int8", copy=False))
            pieces["centrality"].append(cent[keep])

    rows = {
        key: np.concatenate(vals) if vals else np.asarray([], dtype="float32")
        for key, vals in pieces.items()
    }
    counts["selected_entries"] = int(len(rows["score"]))
    counts["selected_signal_entries"] = int(np.count_nonzero(rows["is_signal"] == 1))
    counts["selected_background_entries"] = int(np.count_nonzero(rows["is_signal"] == 0))
    return {"rows": rows, "counts": counts, "cache_dir": str(cache_dir)}


def auc_for(y: np.ndarray, score: np.ndarray) -> float:
    from sklearn.metrics import roc_auc_score

    mask = np.isfinite(score) & np.isin(y, [0, 1])
    if mask.sum() < 2 or len(np.unique(y[mask])) < 2:
        return math.nan
    return float(roc_auc_score(y[mask], score[mask]))


def wp_for_signal_efficiency(y: np.ndarray, score: np.ndarray, target: float = 0.80) -> dict:
    y = np.asarray(y)
    score = np.asarray(score)
    mask = np.isfinite(score) & np.isin(y, [0, 1])
    sig = score[mask & (y == 1)]
    bkg = score[mask & (y == 0)]
    if sig.size <= 0 or bkg.size <= 0:
        return {
            "threshold": math.nan,
            "signal_efficiency": math.nan,
            "background_fake_rate": math.nan,
            "signal_entries": int(sig.size),
            "background_entries": int(bkg.size),
        }
    threshold = float(np.quantile(sig, max(0.0, min(1.0, 1.0 - target))))
    return {
        "threshold": threshold,
        "signal_efficiency": float(np.mean(sig > threshold)),
        "background_fake_rate": float(np.mean(bkg > threshold)),
        "signal_entries": int(sig.size),
        "background_entries": int(bkg.size),
    }


def row_block(dataset: dict, lo: float, hi: float) -> dict:
    rows = dataset["rows"]
    cent = rows["centrality"]
    mask = (cent >= lo) & (cent < hi)
    y = rows["is_signal"][mask]
    score = rows["score"][mask]
    wp = wp_for_signal_efficiency(y, score, 0.80)
    return {
        "centrality_low": lo,
        "centrality_high": hi,
        "centrality_key": f"{lo:g}_{hi:g}",
        "entries": int(mask.sum()),
        "signal_entries": int(np.count_nonzero(y == 1)),
        "background_entries": int(np.count_nonzero(y == 0)),
        "auc": auc_for(y, score),
        "wp80_threshold": wp["threshold"],
        "wp80_signal_efficiency": wp["signal_efficiency"],
        "wp80_background_fake_rate": wp["background_fake_rate"],
    }


def compute_payload(a1_cache_dir: Path, a2_cache_dir: Path) -> dict:
    print("[slide8-full] loading A1 with-cut full-dataset scoring caches", flush=True)
    with_cut = load_scored_rows(a1_cache_dir)
    print("[slide8-full] loading A2 no-cut full-dataset scoring caches", flush=True)
    no_cut = load_scored_rows(a2_cache_dir)

    y_cut = with_cut["rows"]["is_signal"]
    score_cut = with_cut["rows"]["score"]
    y_base = no_cut["rows"]["is_signal"]
    score_base = no_cut["rows"]["score"]

    baseline = {
        "counts": no_cut["counts"],
        "inclusive_auc": auc_for(y_base, score_base),
        "threshold_wp80": wp_for_signal_efficiency(y_base, score_base, 0.80),
        "centrality": {},
    }
    cut = {
        "counts": with_cut["counts"],
        "inclusive_auc": auc_for(y_cut, score_cut),
        "threshold_wp80": wp_for_signal_efficiency(y_cut, score_cut, 0.80),
        "centrality": {},
    }
    for (lo, hi), key in zip(BROAD_CENT_BINS, CENT_KEYS):
        baseline["centrality"][key] = row_block(no_cut, lo, hi)
        cut["centrality"][key] = row_block(with_cut, lo, hi)

    fine5 = {"baseline": {}, "upstream_cut": {}}
    for lo, hi in FINE_CENT_BINS:
        key = f"{lo:g}_{hi:g}"
        fine5["baseline"][key] = row_block(no_cut, lo, hi)
        fine5["upstream_cut"][key] = row_block(with_cut, lo, hi)

    return {
        "schema": "THE57_SLIDE8_FULL_DATASET_SCORING_V1",
        "product": PRODUCT,
        "source_product_name": "full-dataset scoring/application, not 90/10 validation split",
        "a1_with_cut_cache_dir": str(a1_cache_dir),
        "a2_no_cut_cache_dir": str(a2_cache_dir),
        "selection": "15 <= cluster_Et < 35 GeV, 0 <= centrality < 80 for AUC/WP80 cells; global top card uses score-cache counts",
        "baseline": baseline,
        "upstream_cut": cut,
        "fine5": fine5,
    }


def rounded(ax, xy, wh, face, edge="#CBD5E1", lw=1.2, radius=0.018):
    patch = FancyBboxPatch(
        xy,
        wh[0],
        wh[1],
        boxstyle=f"round,pad=0.012,rounding_size={radius}",
        facecolor=face,
        edgecolor=edge,
        linewidth=lw,
        transform=ax.transAxes,
        zorder=1,
    )
    ax.add_patch(patch)
    return patch


def card(ax, x, title, value, sub, *, color, face):
    card_y = 0.645
    rounded(ax, (x, card_y), (0.285, 0.140), face, edge="#D0D5DD")
    ax.text(x + 0.018, card_y + 0.118, title, transform=ax.transAxes, fontsize=13.0, color=INK, fontweight="bold", va="top")
    ax.text(x + 0.018, card_y + 0.077, value, transform=ax.transAxes, fontsize=21.0, color=color, fontweight="bold", va="top")
    ax.text(x + 0.018, card_y + 0.031, sub, transform=ax.transAxes, fontsize=11.0, color=MUTED, va="top")


def slide_rows(payload: dict) -> tuple[list[dict], list[dict]]:
    base = payload["baseline"]
    cut = payload["upstream_cut"]
    rows = []
    for key, label in zip(CENT_KEYS, CENT_LABELS):
        b = base["centrality"][key]
        c = cut["centrality"][key]
        base_fake = safe_float(b["wp80_background_fake_rate"])
        cut_fake = safe_float(c["wp80_background_fake_rate"])
        rows.append(
            {
                "centrality": label,
                "centrality_key": key,
                "baseline_auc": safe_float(b["auc"]),
                "cut_auc": safe_float(c["auc"]),
                "delta_auc": safe_float(c["auc"]) - safe_float(b["auc"]),
                "baseline_wp80_fake": base_fake,
                "cut_wp80_fake": cut_fake,
                "delta_wp80_fake": cut_fake - base_fake,
                "baseline_threshold": safe_float(b["wp80_threshold"]),
                "cut_threshold": safe_float(c["wp80_threshold"]),
                "baseline_entries": int(b["entries"]),
                "cut_entries": int(c["entries"]),
                "removed_entries": int(b["entries"]) - int(c["entries"]),
                "removed_fraction": (int(b["entries"]) - int(c["entries"])) / int(b["entries"]) if int(b["entries"]) else math.nan,
            }
        )
    fine_rows = []
    for lo in range(0, 80, 5):
        key = f"{lo}_{lo + 5}"
        b = payload["fine5"]["baseline"][key]
        c = payload["fine5"]["upstream_cut"][key]
        base_fake = safe_float(b["wp80_background_fake_rate"])
        cut_fake = safe_float(c["wp80_background_fake_rate"])
        ratio = cut_fake / base_fake if base_fake > 0 else math.nan
        fine_rows.append(
            {
                "centrality": f"{lo}-{lo + 5}%",
                "centrality_mid": lo + 2.5,
                "centrality_key": key,
                "baseline_wp80_fake": base_fake,
                "cut_wp80_fake": cut_fake,
                "fake_ratio": ratio,
                "fake_ratio_percent_change": 100.0 * (ratio - 1.0) if np.isfinite(ratio) else math.nan,
                "baseline_entries": int(b["entries"]),
                "cut_entries": int(c["entries"]),
                "removed_entries": int(b["entries"]) - int(c["entries"]),
            }
        )
    return rows, fine_rows


def render_slide(payload: dict, png: Path, csv_path: Path, fine_csv_path: Path, script_path: Path, manifest_path: Path) -> None:
    rows, fine_rows = slide_rows(payload)
    png.parent.mkdir(parents=True, exist_ok=True)

    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    with fine_csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(fine_rows[0]))
        writer.writeheader()
        writer.writerows(fine_rows)

    fig = plt.figure(figsize=(12.8, 7.2), dpi=200)
    fig.subplots_adjust(left=0, right=1, bottom=0, top=1)
    fig.patch.set_facecolor("white")

    canvas = fig.add_axes([0, 0, 1, 1])
    canvas.axis("off")
    canvas.text(
        0.055,
        0.970,
        "BDT response is stable after the total-calo floor veto",
        fontsize=26.0,
        color=INK,
        fontweight="bold",
        va="top",
    )
    rounded(canvas, (0.055, 0.820), (0.890, 0.076), "#F1F5F9", edge="#D9E2EC")
    canvas.text(
        0.500,
        0.872,
        "Full-dataset scoring/application, 15-35 GeV photons; no-cut diagnostic vs default model with the event floor veto.",
        fontsize=12.8,
        color=INK,
        ha="center",
        va="center",
        fontweight="bold",
    )
    canvas.text(
        0.500,
        0.838,
        "Takeaway: the veto removes a real low-calo population while leaving BDT discrimination essentially unchanged.",
        fontsize=12.4,
        color=GREEN,
        ha="center",
        va="center",
        fontweight="bold",
    )

    base = payload["baseline"]
    cut = payload["upstream_cut"]
    base_total = int(base["counts"].get("scored_entries") or base["counts"].get("total_entries") or 0)
    cut_total = int(cut["counts"].get("scored_entries") or cut["counts"].get("total_entries") or 0)
    global_removed = base_total - cut_total
    global_removed_frac = global_removed / base_total if base_total else math.nan
    auc_delta = safe_float(cut["inclusive_auc"]) - safe_float(base["inclusive_auc"])
    fake_ratio = safe_float(cut["threshold_wp80"]["background_fake_rate"]) / safe_float(base["threshold_wp80"]["background_fake_rate"])
    fake_delta = safe_float(cut["threshold_wp80"]["background_fake_rate"]) - safe_float(base["threshold_wp80"]["background_fake_rate"])

    card(
        canvas,
        0.055,
        "Full-dataset rows removed by floor veto",
        f"{global_removed:,}",
        f"{pct(global_removed_frac, 2)} of no-cut scored rows",
        color=RED,
        face=CARD_RED,
    )
    card(
        canvas,
        0.365,
        "Inclusive AUC",
        f"{safe_float(cut['inclusive_auc']):.6f}",
        f"no-cut {safe_float(base['inclusive_auc']):.6f}; change {auc_delta:+.6f}",
        color=GREEN,
        face=CARD_GREEN,
    )
    card(
        canvas,
        0.675,
        "Inclusive WP80 fake-rate ratio",
        f"{fake_ratio:.3f}",
        f"with cut / no cut; change {100*fake_delta:+.2f} pp",
        color=BLUE,
        face=CARD_BLUE,
    )

    ax_auc = fig.add_axes([0.075, 0.322, 0.415, 0.265])
    ax_auc.set_facecolor("white")
    x = np.arange(len(rows))
    width = 0.34
    base_auc = np.array([r["baseline_auc"] for r in rows], dtype=float)
    cut_auc = np.array([r["cut_auc"] for r in rows], dtype=float)
    y_min = max(0.0, min(np.nanmin(base_auc), np.nanmin(cut_auc)) - 0.012)
    y_max = min(1.0, max(np.nanmax(base_auc), np.nanmax(cut_auc)) + 0.012)
    ax_auc.bar(x - width / 2, base_auc, width=width, color="#B7C9DD", edgecolor="#5D6B7A", linewidth=0.7, label="No-cut diagnostic")
    ax_auc.bar(x + width / 2, cut_auc, width=width, color=BLUE, edgecolor="#174A76", linewidth=0.7, label="With floor veto + retrain")
    ax_auc.set_title("AUC by centrality", fontsize=15, color=INK, fontweight="bold", pad=8)
    ax_auc.set_ylabel("AUC", fontsize=13)
    ax_auc.set_xticks(x)
    ax_auc.set_xticklabels([r["centrality"] for r in rows], fontsize=12)
    ax_auc.set_ylim(y_min, y_max)
    ax_auc.grid(axis="y", color="#E2E8F0", linewidth=0.8)
    ax_auc.tick_params(axis="y", labelsize=11)
    ax_auc.legend(loc="upper left", fontsize=10.5, frameon=False)
    for i, row in enumerate(rows):
        ax_auc.text(i, max(row["baseline_auc"], row["cut_auc"]) + 0.002, f"{row['delta_auc']:+.4f}", ha="center", va="bottom", fontsize=12.2, color=INK, fontweight="bold")

    ax_fake = fig.add_axes([0.565, 0.322, 0.375, 0.265])
    fine_x = np.array([float(r["centrality_mid"]) for r in fine_rows], dtype=float)
    fake_ratio_by_cent = np.array([float(r["fake_ratio"]) for r in fine_rows], dtype=float)
    ax_fake.axhline(1.0, color="#8492A6", linewidth=1.35, linestyle="--", label="No change")
    ax_fake.plot(fine_x, fake_ratio_by_cent, marker="o", markersize=4.8, linewidth=2.0, color=ORANGE, label="With cut / no cut")
    ax_fake.set_title("Fine-bin fake-rate ratio with WP80 rederived per 5%", fontsize=13.5, color=INK, fontweight="bold", pad=8)
    ax_fake.set_ylabel("Fake-rate ratio", fontsize=11.8)
    ax_fake.set_xlabel("Centrality percentile", fontsize=11.5)
    ax_fake.set_xlim(0, 80)
    ax_fake.set_xticks(np.arange(0, 81, 10))
    ax_fake.set_ylim(float(np.nanmin(fake_ratio_by_cent)) - 0.006, float(np.nanmax(fake_ratio_by_cent)) + 0.006)
    ax_fake.grid(axis="y", color="#E2E8F0", linewidth=0.8)
    ax_fake.tick_params(axis="both", labelsize=10.8)
    ax_fake.legend(loc="upper right", fontsize=9.8, frameon=False)
    finite = np.flatnonzero(np.isfinite(fake_ratio_by_cent))
    if finite.size:
        min_idx = int(finite[np.nanargmin(fake_ratio_by_cent[finite])])
        max_idx = int(finite[np.nanargmax(fake_ratio_by_cent[finite])])
        for idx in sorted({min_idx, max_idx}):
            ratio = fake_ratio_by_cent[idx]
            dy = 100.0 * (ratio - 1.0)
            yoff = 0.002 if idx == max_idx else -0.002
            va = "bottom" if idx == max_idx else "top"
            ax_fake.text(fine_x[idx], ratio + yoff, f"{dy:+.2f}%", ha="center", va=va, fontsize=10.2, color=INK, fontweight="bold")

    rounded(canvas, (0.055, 0.045), (0.890, 0.180), "#FFF8F8", edge="#F3D4D4")
    canvas.text(0.080, 0.206, "15-35 GeV full-dataset scored rows removed in each centrality bin", transform=canvas.transAxes, fontsize=15.0, color=INK, fontweight="bold", va="top")
    for xx, row in zip([0.185, 0.500, 0.815], rows):
        rounded(canvas, (xx - 0.108, 0.066), (0.216, 0.094), "#FFFFFF", edge="#E2E8F0", lw=0.9, radius=0.010)
        canvas.text(xx, 0.149, row["centrality"], transform=canvas.transAxes, fontsize=14.8, color=INK, fontweight="bold", ha="center", va="top")
        canvas.text(
            xx,
            0.111,
            f"{row['removed_entries']:,} / {row['baseline_entries']:,} rows removed",
            transform=canvas.transAxes,
            fontsize=10.5,
            color=RED if row["removed_fraction"] > 0.01 else MUTED,
            fontweight="bold",
            ha="center",
            va="top",
        )
        canvas.text(xx, 0.083, f"{pct(row['removed_fraction'], 2)} of that bin", transform=canvas.transAxes, fontsize=10.2, color=MUTED, ha="center", va="top")

    script_path.write_text(
        "\n".join(
            [
                "# Speaker Notes",
                "",
                "- This remake uses full-dataset scoring/application caches, not the capped validation-split cache.",
                "- The comparison is no-cut diagnostic model versus the default model trained with the total-calo floor veto.",
                "- The bottom row and fine-bin fake-rate curve use 15 <= cluster ET < 35 GeV and centrality 0-80%.",
            ]
        )
        + "\n"
    )
    manifest = {
        "schema": "THE57_SLIDE8_FULL_DATASET_SCORING_MANIFEST_V1",
        "png": str(png),
        "payload": str(manifest_path.with_name("the57_slide8_full_dataset_scoring_payload.json")),
        "csv": str(csv_path),
        "fine5_csv": str(fine_csv_path),
        "script": str(script_path),
        "product": PRODUCT,
        "source_product_name": payload.get("source_product_name"),
        "global_removed_rows": global_removed,
        "global_removed_fraction": global_removed_frac,
        "inclusive_auc_with_cut": safe_float(cut["inclusive_auc"]),
        "inclusive_auc_no_cut": safe_float(base["inclusive_auc"]),
        "inclusive_wp80_fake_ratio": fake_ratio,
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    fig.savefig(png)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--a1-cache-dir", type=Path, help="A1 with-cut full-dataset score_caches directory")
    parser.add_argument("--a2-cache-dir", type=Path, help="A2 no-cut full-dataset score_caches directory")
    parser.add_argument("--payload", type=Path, default=OUT / "the57_slide8_full_dataset_scoring_payload.json")
    parser.add_argument("--png", type=Path, default=OUT / "the57_slide8_full_dataset_scoring.png")
    parser.add_argument("--compute", action="store_true", help="Stream caches and write payload JSON")
    parser.add_argument("--render", action="store_true", help="Render slide PNG from payload JSON")
    parser.add_argument("--print-payload", action="store_true", help="Print payload JSON to stdout after compute")
    args = parser.parse_args()

    if not args.compute and not args.render:
        args.compute = bool(args.a1_cache_dir and args.a2_cache_dir)
        args.render = True

    if args.compute:
        if args.a1_cache_dir is None or args.a2_cache_dir is None:
            raise SystemExit("--compute requires --a1-cache-dir and --a2-cache-dir")
        payload = compute_payload(args.a1_cache_dir, args.a2_cache_dir)
        args.payload.parent.mkdir(parents=True, exist_ok=True)
        args.payload.write_text(json.dumps(payload, indent=2, sort_keys=True, allow_nan=True) + "\n")
        if args.print_payload:
            print("THE57_SLIDE8_PAYLOAD_JSON_BEGIN")
            print(json.dumps(payload, sort_keys=True, allow_nan=True))
            print("THE57_SLIDE8_PAYLOAD_JSON_END")

    if args.render:
        payload = json.loads(args.payload.read_text())
        render_slide(
            payload,
            args.png,
            args.png.with_name("the57_slide8_full_dataset_scoring.csv"),
            args.png.with_name("the57_slide8_full_dataset_scoring_fine5.csv"),
            args.png.with_name("the57_slide8_full_dataset_scoring_script.md"),
            args.png.with_name("the57_slide8_full_dataset_scoring_manifest.json"),
        )
        print(f"Wrote {args.png}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
