#!/usr/bin/env python3
"""Reduce and render the THE-11 weighted full-stat BDT-isolation slide."""

from __future__ import annotations

import argparse
import json
import math
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import colors, font_manager  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from matplotlib.patches import FancyBboxPatch  # noqa: E402
import numpy as np  # noqa: E402


REPO = Path(__file__).resolve().parents[4]
DEFAULT_LOCAL_OUT = (
    REPO
    / "dataOutput/auauTightBDTValidation/the11_weighted_bdt_iso_fullstat_basev3e_20260606/slideReady"
)
DEFAULT_REMOTE_REPORT = Path(
    "/sphenix/tg/tg01/bulk/jbennett/thesisAnaTraining/auauTightBDT_eiso_cone_raw_20260518_2220/"
    "reports/model_validation_condor_the11_weighted_bdt_iso_fullstat_basev3e_20260606_1645"
)

SCORE_COLUMN = "score_centInput_pt1535"
ISO_COLUMN = "reco_eiso_r30"
PT_RANGE = (15.0, 35.0)
CENT_BINS = [(0.0, 20.0, "0-20%"), (20.0, 50.0, "20-50%"), (50.0, 80.0, "50-80%")]
CLASS_ROWS = [(1, "truth photons", "#a32035"), (0, "inclusive background", "#1764a8")]
X_BINS = np.linspace(-10.0, 18.0, 57)
Y_BINS = np.linspace(0.0, 1.0, 51)
W, H, DPI = 2560, 1440, 200
WHITE = "#ffffff"
INK = "#151515"
MUTED = "#5c6470"
GRID = "#dfe5ec"
SUBTITLE_BG = "#f4f7fb"
SUBTITLE_EDGE = "#cad6e4"
CARD_BG = "#fff8df"
CARD_EDGE = "#d8c67f"


def kbird_cmap() -> LinearSegmentedColormap:
    stops = np.linspace(0.0, 1.0, 9)
    red = [0.2082, 0.0592, 0.0780, 0.0232, 0.1802, 0.5301, 0.8186, 0.9956, 0.9764]
    green = [0.1664, 0.3599, 0.5041, 0.6419, 0.7178, 0.7492, 0.7328, 0.7862, 0.9832]
    blue = [0.5293, 0.8684, 0.8385, 0.8385, 0.6425, 0.4662, 0.3499, 0.1968, 0.0539]
    return LinearSegmentedColormap.from_list("root_kbird_like", list(zip(stops, zip(red, green, blue))), N=256)


def setup_style() -> None:
    available = {f.name for f in font_manager.fontManager.ttflist}
    for family in ("Times New Roman", "Times", "DejaVu Serif"):
        if family in available:
            plt.rcParams["font.family"] = family
            break
    plt.rcParams.update(
        {
            "figure.facecolor": WHITE,
            "savefig.facecolor": WHITE,
            "axes.facecolor": WHITE,
            "axes.edgecolor": INK,
            "axes.linewidth": 0.9,
            "axes.labelcolor": INK,
            "xtick.color": INK,
            "ytick.color": INK,
            "text.color": INK,
            "mathtext.default": "regular",
        }
    )


def read_cache_paths(report_dir: Path, cache_manifest: Path | None = None) -> list[Path]:
    manifest = cache_manifest or (report_dir / "score_caches.list")
    if manifest.exists():
        paths = [Path(line.strip()) for line in manifest.read_text().splitlines() if line.strip()]
    else:
        paths = sorted((report_dir / "score_caches").glob("score_cache_*.npz"))
    if not paths:
        raise FileNotFoundError(f"No score caches found from {manifest} or {report_dir / 'score_caches'}")
    return paths


def weighted_quantile(values: np.ndarray, weights: np.ndarray, quantile: float) -> float:
    mask = np.isfinite(values) & np.isfinite(weights) & (weights > 0.0)
    if not np.any(mask):
        return float("nan")
    v = values[mask].astype("float64", copy=False)
    w = weights[mask].astype("float64", copy=False)
    order = np.argsort(v)
    v = v[order]
    w = w[order]
    cdf = np.cumsum(w)
    target = quantile * cdf[-1]
    return float(np.interp(target, cdf, v))


def weighted_corr(x: np.ndarray, y: np.ndarray, weights: np.ndarray) -> float:
    mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(weights) & (weights > 0.0)
    if np.count_nonzero(mask) < 3:
        return float("nan")
    x = x[mask].astype("float64", copy=False)
    y = y[mask].astype("float64", copy=False)
    w = weights[mask].astype("float64", copy=False)
    wsum = np.sum(w)
    mx = np.sum(w * x) / wsum
    my = np.sum(w * y) / wsum
    cov = np.sum(w * (x - mx) * (y - my)) / wsum
    vx = np.sum(w * (x - mx) ** 2) / wsum
    vy = np.sum(w * (y - my) ** 2) / wsum
    if vx <= 0.0 or vy <= 0.0:
        return float("nan")
    return float(cov / math.sqrt(vx * vy))


def reduce_caches(report_dir: Path, out_dir: Path, cache_manifest: Path | None = None) -> tuple[Path, Path]:
    cache_paths = read_cache_paths(report_dir, cache_manifest)
    if len(cache_paths) != 25:
        raise RuntimeError(f"Expected 25 full-stat score caches; found {len(cache_paths)}")

    hist = np.zeros((len(CLASS_ROWS), len(CENT_BINS), len(X_BINS) - 1, len(Y_BINS) - 1), dtype="float64")
    cells: dict[tuple[int, int], dict[str, list[np.ndarray] | float | int]] = {}
    for class_idx, _ in enumerate(CLASS_ROWS):
        for cent_idx, _ in enumerate(CENT_BINS):
            cells[(class_idx, cent_idx)] = {"eiso": [], "score": [], "weight": [], "n": 0, "wsum": 0.0}

    required = [
        "is_signal",
        "cluster_Et",
        "centrality",
        ISO_COLUMN,
        SCORE_COLUMN,
        "event_weight",
        "source_sample",
        "input_file_index",
        "input_tree_entry",
        "reco_eiso_r40",
    ]
    for cache_path in cache_paths:
        with np.load(cache_path, allow_pickle=True) as data:
            missing = [name for name in required if name not in data.files]
            if missing:
                raise KeyError(f"{cache_path} missing required columns: {missing}")
            is_signal = data["is_signal"].astype("int8", copy=False)
            et = data["cluster_Et"].astype("float32", copy=False)
            cent = data["centrality"].astype("float32", copy=False)
            eiso = data[ISO_COLUMN].astype("float32", copy=False)
            score = data[SCORE_COLUMN].astype("float32", copy=False)
            weight = data["event_weight"].astype("float64", copy=False)

            finite = (
                np.isfinite(et)
                & np.isfinite(cent)
                & np.isfinite(eiso)
                & np.isfinite(score)
                & np.isfinite(weight)
                & (weight > 0.0)
                & (et >= PT_RANGE[0])
                & (et < PT_RANGE[1])
                & (cent >= 0.0)
                & (cent < 80.0)
            )
            for class_idx, (label_value, _, _) in enumerate(CLASS_ROWS):
                class_mask = finite & (is_signal == label_value)
                for cent_idx, (cent_lo, cent_hi, _) in enumerate(CENT_BINS):
                    mask = class_mask & (cent >= cent_lo) & (cent < cent_hi)
                    if not np.any(mask):
                        continue
                    h, _, _ = np.histogram2d(eiso[mask], score[mask], bins=[X_BINS, Y_BINS], weights=weight[mask])
                    hist[class_idx, cent_idx] += h
                    cell = cells[(class_idx, cent_idx)]
                    cell["eiso"].append(eiso[mask])
                    cell["score"].append(score[mask])
                    cell["weight"].append(weight[mask])
                    cell["n"] = int(cell["n"]) + int(np.count_nonzero(mask))
                    cell["wsum"] = float(cell["wsum"]) + float(np.sum(weight[mask]))

    summaries: dict[str, object] = {}
    for class_idx, (_, class_label, _) in enumerate(CLASS_ROWS):
        summaries[class_label] = {}
        for cent_idx, (_, _, cent_label) in enumerate(CENT_BINS):
            cell = cells[(class_idx, cent_idx)]
            eiso = np.concatenate(cell["eiso"]) if cell["eiso"] else np.array([], dtype="float32")
            score = np.concatenate(cell["score"]) if cell["score"] else np.array([], dtype="float32")
            weight = np.concatenate(cell["weight"]) if cell["weight"] else np.array([], dtype="float64")
            positive_busy = (eiso > 5.0) & (score >= 0.60)
            clean_high = (eiso < 2.0) & (score >= 0.60)
            summaries[class_label][cent_label] = {
                "n": int(cell["n"]),
                "weight_sum": float(cell["wsum"]),
                "median_eiso": weighted_quantile(eiso, weight, 0.5),
                "q25_eiso": weighted_quantile(eiso, weight, 0.25),
                "q75_eiso": weighted_quantile(eiso, weight, 0.75),
                "median_score": weighted_quantile(score, weight, 0.5),
                "weighted_pearson_score_eiso": weighted_corr(score, eiso, weight),
                "weighted_fraction_score_ge_0p60_eiso_lt_2": float(np.sum(weight[clean_high]) / np.sum(weight))
                if np.sum(weight) > 0.0
                else float("nan"),
                "weighted_fraction_score_ge_0p60_eiso_gt_5": float(np.sum(weight[positive_busy]) / np.sum(weight))
                if np.sum(weight) > 0.0
                else float("nan"),
            }

    out_dir.mkdir(parents=True, exist_ok=True)
    npz_path = out_dir / "weighted_bdt_iso_fullstat_summary_r30.npz"
    json_path = out_dir / "weighted_bdt_iso_fullstat_summary_r30.json"
    np.savez_compressed(
        npz_path,
        hist=hist.astype("float64"),
        x_bins=X_BINS.astype("float32"),
        y_bins=Y_BINS.astype("float32"),
        class_labels=np.asarray([row[1] for row in CLASS_ROWS], dtype=object),
        centrality_labels=np.asarray([row[2] for row in CENT_BINS], dtype=object),
    )
    payload = {
        "source_report": str(report_dir),
        "cache_count": len(cache_paths),
        "score_column": SCORE_COLUMN,
        "isolation_column": ISO_COLUMN,
        "pt_range_gev": list(PT_RANGE),
        "normalization": "Histograms are event-weighted; plotting normalizes each class/centrality panel to unit area.",
        "summaries": summaries,
    }
    json_path.write_text(json.dumps(payload, indent=2, sort_keys=True))
    return npz_path, json_path


def add_box(fig: plt.Figure, xywh: tuple[float, float, float, float], face: str, edge: str, zorder: float = 0.2) -> None:
    fig.add_artist(
        FancyBboxPatch(
            (xywh[0], xywh[1]),
            xywh[2],
            xywh[3],
            boxstyle="round,pad=0.010,rounding_size=0.010",
            transform=fig.transFigure,
            facecolor=face,
            edgecolor=edge,
            linewidth=1.0,
            zorder=zorder,
        )
    )


def write_text(fig: plt.Figure, x: float, y: float, text: str, width: int, **kwargs) -> None:
    fig.text(x, y, textwrap.fill(text, width=width), ha="left", va="top", linespacing=1.12, **kwargs)


def plot_slide(summary_npz: Path, summary_json: Path, out_dir: Path) -> tuple[Path, Path, Path]:
    setup_style()
    with np.load(summary_npz, allow_pickle=True) as data:
        hist = data["hist"].astype("float64")
        x_bins = data["x_bins"].astype("float64")
        y_bins = data["y_bins"].astype("float64")
        class_labels = [str(x) for x in data["class_labels"]]
        cent_labels = [str(x) for x in data["centrality_labels"]]
    meta = json.loads(summary_json.read_text())

    panel_prob = np.zeros_like(hist)
    for i in range(hist.shape[0]):
        for j in range(hist.shape[1]):
            total = float(np.sum(hist[i, j]))
            if total > 0.0:
                panel_prob[i, j] = hist[i, j] / total
    positive = panel_prob[panel_prob > 0.0]
    vmax = float(np.quantile(positive, 0.995)) if positive.size else 1.0
    vmin = max(vmax / 1500.0, 1e-8)
    norm = colors.LogNorm(vmin=vmin, vmax=vmax)
    cmap = kbird_cmap()

    fig = plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI)
    fig.patch.set_facecolor(WHITE)
    fig.text(
        0.055,
        0.955,
        "Full-stat weighted sample confirms distinct BDT-isolation regions",
        ha="left",
        va="top",
        fontsize=23,
        fontweight="bold",
    )
    add_box(fig, (0.055, 0.866, 0.885, 0.064), SUBTITLE_BG, SUBTITLE_EDGE)
    write_text(
        fig,
        0.071,
        0.914,
        "Au+Au embedded Photon+Jet + Inclusive+Jet, 15-35 GeV clusters. Color is weighted candidates per bin on a log scale; each panel is normalized within its own class and centrality bin, so the shape comparison is not driven by raw sample size.",
        142,
        fontsize=12.4,
        color="#253243",
    )

    left0, bottom0 = 0.125, 0.300
    ax_w, ax_h = 0.238, 0.232
    x_gap, y_gap = 0.028, 0.058
    mesh = None
    for row_idx, (_, class_name, class_color) in enumerate(CLASS_ROWS):
        for col_idx, (_, _, cent_label) in enumerate(CENT_BINS):
            ax = fig.add_axes([left0 + col_idx * (ax_w + x_gap), bottom0 + (1 - row_idx) * (ax_h + y_gap), ax_w, ax_h])
            values = np.ma.masked_less_equal(panel_prob[row_idx, col_idx].T, 0.0)
            mesh = ax.pcolormesh(x_bins, y_bins, values, cmap=cmap, norm=norm, shading="auto")
            ax.set_xlim(-10, 18)
            ax.set_ylim(0, 1)
            ax.grid(True, color=GRID, linewidth=0.45, alpha=0.75)
            ax.tick_params(axis="both", labelsize=9.5, length=3)
            if row_idx == 1:
                ax.set_xlabel(r"reco $E_T^{iso}$, $\Delta R < 0.3$ [GeV]", fontsize=10.4)
            else:
                ax.set_xticklabels([])
            if col_idx == 0:
                ax.set_ylabel("BDT score", fontsize=11)
            else:
                ax.set_yticklabels([])
            if row_idx == 0:
                ax.set_title(cent_label + " centrality", fontsize=13, fontweight="bold", pad=6)
            summary = meta["summaries"][class_name][cent_label]
            ax.text(
                0.02,
                0.965,
                f"N={summary['n']:,}\n"
                f"med iso={summary['median_eiso']:.2f} GeV\n"
                f"rho={summary['weighted_pearson_score_eiso']:+.2f}",
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=8.8,
                bbox=dict(boxstyle="round,pad=0.22", facecolor="white", edgecolor="#b8c2ce", alpha=0.90),
            )
        fig.text(
            0.063,
            bottom0 + (1 - row_idx) * (ax_h + y_gap) + ax_h / 2,
            class_name,
            ha="center",
            va="center",
            rotation=90,
            fontsize=15,
            fontweight="bold",
            color=class_color,
        )

    if mesh is not None:
        cax = fig.add_axes([0.925, 0.333, 0.018, 0.470])
        cb = fig.colorbar(mesh, cax=cax)
        cb.ax.tick_params(labelsize=9)
        cb.set_label("weighted fraction / bin\nlog scale", fontsize=9.5)

    add_box(fig, (0.074, 0.038, 0.833, 0.130), CARD_BG, CARD_EDGE)
    write_text(
        fig,
        0.097,
        0.145,
        "Why this replaces the sampled plot: this uses the full embedded source, event weights, and explicit R=0.3 cone isolation rather than a small validation sample.",
        46,
        fontsize=10.8,
        color=INK,
    )
    write_text(
        fig,
        0.376,
        0.145,
        "How to read it: truth photons concentrate at high BDT score and low cone energy. Inclusive background spreads into lower-score and positive-isolation regions.",
        50,
        fontsize=10.8,
        color=INK,
    )
    write_text(
        fig,
        0.671,
        0.145,
        "Physics implication: BDT and isolation are not independent background axes; the key follow-up is the correlated background surviving both working points.",
        45,
        fontsize=10.8,
        color=INK,
    )

    out_dir.mkdir(parents=True, exist_ok=True)
    png = out_dir / "weighted_bdt_isolation_fullstat_slide.png"
    script = out_dir / "weighted_bdt_isolation_fullstat_speaker_script.md"
    manifest = out_dir / "weighted_bdt_isolation_fullstat_manifest.json"
    fig.savefig(png, dpi=DPI)
    plt.close(fig)

    script.write_text(
        "# THE-11 Weighted BDT-Isolation Full-Stat Slide Script\n\n"
        "This slide is the full-stat version of the BDT versus isolation picture. "
        "The important change is that I am no longer showing a small sampled validation cache. "
        "This is the embedded Photon+Jet and Inclusive+Jet source, with the event weights carried through, "
        "and with the R equals 0.3 reconstructed isolation branch written explicitly in the score caches.\n\n"
        "The top row is truth photons. They concentrate at high BDT score and low reconstructed cone energy, "
        "which is the behavior we want from an isolated photon population. The bottom row is the inclusive "
        "background. It fills a broader positive-isolation tail and lower BDT score region, which is the "
        "hadronic activity around the EM cluster showing up in both variables.\n\n"
        "The main takeaway is that isolation and BDT are not independent axes for the background. The BDT is "
        "doing useful photon-likeness classification, but the remaining background structure is correlated "
        "with cone activity, so this is exactly the place where the ABCD assumption needs a quantitative "
        "closure or systematic treatment rather than a visual assumption of factorization.\n"
    )
    manifest.write_text(
        json.dumps(
            {
                "png": str(png),
                "speaker_script": str(script),
                "summary_npz": str(summary_npz),
                "summary_json": str(summary_json),
                "render_size_px": [W, H],
                "score_column": SCORE_COLUMN,
                "isolation_column": ISO_COLUMN,
                "normalization": meta["normalization"],
            },
            indent=2,
            sort_keys=True,
        )
    )
    return png, script, manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["reduce", "plot", "both"], default="both")
    parser.add_argument("--report-dir", type=Path, default=DEFAULT_REMOTE_REPORT)
    parser.add_argument("--cache-manifest", type=Path, default=None)
    parser.add_argument("--summary-npz", type=Path, default=None)
    parser.add_argument("--summary-json", type=Path, default=None)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_LOCAL_OUT)
    args = parser.parse_args()

    summary_npz = args.summary_npz
    summary_json = args.summary_json
    if args.mode in {"reduce", "both"}:
        summary_npz, summary_json = reduce_caches(args.report_dir, args.out_dir, args.cache_manifest)
        print(f"WROTE_SUMMARY_NPZ={summary_npz}")
        print(f"WROTE_SUMMARY_JSON={summary_json}")
    if args.mode in {"plot", "both"}:
        if summary_npz is None or summary_json is None:
            raise SystemExit("--summary-npz and --summary-json are required for plot mode")
        png, script, manifest = plot_slide(summary_npz, summary_json, args.out_dir)
        print(f"WROTE_PNG={png}")
        print(f"WROTE_SCRIPT={script}")
        print(f"WROTE_MANIFEST={manifest}")


if __name__ == "__main__":
    main()
