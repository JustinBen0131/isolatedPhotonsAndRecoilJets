#!/usr/bin/env python3
"""Make quick plot-only variants for default isolation vs BDT score."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402
from PIL import Image, ImageDraw, ImageFont  # noqa: E402
import numpy as np  # noqa: E402


REPO = Path(__file__).resolve().parents[3]
REPORT = REPO / "dataOutput/auauTightBDTValidation/model_validation_condor_20260511_194832"
MANIFEST = REPORT / "score_caches.local.list"
OUT_DIR = REPORT / "slideReady/default_eiso_vs_bdt_score_20260604/plot_variants"

SCORE = "score_ptFine_cent7"
EISO = "reco_eiso"
PT_RANGE = (15.0, 35.0)
CENT_BINS = [(0.0, 20.0, "0-20%"), (20.0, 50.0, "20-50%"), (50.0, 80.0, "50-80%")]

RED = "#c84c4c"
BLUE = "#315f9c"
INK = "#171717"
MUTED = "#5b6472"


def load_arrays() -> dict[str, np.ndarray]:
    arrays: dict[str, list[np.ndarray]] = {k: [] for k in ["eiso", "score", "is_signal", "et", "cent"]}
    for line in MANIFEST.read_text().splitlines():
        if not line.strip():
            continue
        path = REPO / line.strip()
        data = np.load(path, allow_pickle=True)
        arrays["eiso"].append(data[EISO].astype("float32", copy=False))
        arrays["score"].append(data[SCORE].astype("float32", copy=False))
        arrays["is_signal"].append(data["is_signal"].astype("int8", copy=False))
        arrays["et"].append(data["cluster_Et"].astype("float32", copy=False))
        arrays["cent"].append(data["centrality"].astype("float32", copy=False))
    out = {k: np.concatenate(v) for k, v in arrays.items()}
    selected = (
        np.isfinite(out["eiso"])
        & np.isfinite(out["score"])
        & (out["et"] >= PT_RANGE[0])
        & (out["et"] < PT_RANGE[1])
        & (out["cent"] >= 0.0)
        & (out["cent"] < 80.0)
    )
    return {k: v[selected] for k, v in out.items()}


def setup_ax(ax, title: str, xlim: tuple[float, float], ylim: tuple[float, float] = (0.0, 1.0)) -> None:
    ax.set_title(title, fontsize=15, weight="bold")
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xlabel(r"default $E_T^{iso}$ ($reco\_eiso$) [GeV]", fontsize=12)
    ax.set_ylabel("BDT score", fontsize=12)
    ax.grid(True, color="#dfe4ea", lw=0.7)
    ax.tick_params(labelsize=10)


def hist2d(ax, x, y, xlim, bins=(90, 70), title="") -> None:
    h = ax.hist2d(x, y, bins=bins, range=[xlim, [0.0, 1.0]], norm=LogNorm(vmin=1), cmap="viridis")
    plt.colorbar(h[3], ax=ax, pad=0.012, label="candidates / bin")
    setup_ax(ax, title, xlim)


def quantiles_by_bin(x: np.ndarray, y: np.ndarray, bins: np.ndarray, min_n: int = 30):
    centers, q16, q25, q50, q75, q84 = [], [], [], [], [], []
    for lo, hi in zip(bins[:-1], bins[1:]):
        mask = (x >= lo) & (x < hi)
        if mask.sum() < min_n:
            continue
        vals = y[mask]
        centers.append(0.5 * (lo + hi))
        q16.append(np.nanpercentile(vals, 16))
        q25.append(np.nanpercentile(vals, 25))
        q50.append(np.nanpercentile(vals, 50))
        q75.append(np.nanpercentile(vals, 75))
        q84.append(np.nanpercentile(vals, 84))
    return map(np.asarray, (centers, q16, q25, q50, q75, q84))


def smooth_hist(h: np.ndarray, passes: int = 2) -> np.ndarray:
    kernel = np.array([1.0, 2.0, 1.0], dtype="float64")
    kernel /= kernel.sum()
    out = h.astype("float64", copy=True)
    for _ in range(passes):
        out = np.apply_along_axis(lambda v: np.convolve(v, kernel, mode="same"), 0, out)
        out = np.apply_along_axis(lambda v: np.convolve(v, kernel, mode="same"), 1, out)
    return out


def save_log_density(data: dict[str, np.ndarray], xlim: tuple[float, float]) -> Path:
    fig, ax = plt.subplots(figsize=(8.4, 6.2), dpi=180)
    hist2d(ax, data["eiso"], data["score"], xlim, title="Log-density view: default isolation vs BDT score")
    out = OUT_DIR / "01_log_density_default_eiso_vs_bdt.png"
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)
    return out


def save_rank_density(data: dict[str, np.ndarray]) -> Path:
    x = data["eiso"]
    y = data["score"]
    order = np.argsort(x)
    rank = np.empty_like(x, dtype="float32")
    rank[order] = np.linspace(0.0, 100.0, len(x), dtype="float32")
    fig, ax = plt.subplots(figsize=(8.4, 6.2), dpi=180)
    h = ax.hist2d(rank, y, bins=(90, 70), range=[[0, 100], [0, 1]], norm=LogNorm(vmin=1), cmap="magma")
    plt.colorbar(h[3], ax=ax, pad=0.012, label="candidates / bin")
    ax.set_title("Percentile view: isolation rank vs BDT score", fontsize=15, weight="bold")
    ax.set_xlabel("default isolation percentile", fontsize=12)
    ax.set_ylabel("BDT score", fontsize=12)
    ax.grid(True, color="#dfe4ea", lw=0.7)
    out = OUT_DIR / "02_isolation_percentile_density.png"
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)
    return out


def save_contours(data: dict[str, np.ndarray], xlim: tuple[float, float]) -> Path:
    fig, ax = plt.subplots(figsize=(8.4, 6.2), dpi=180)
    xedges = np.linspace(xlim[0], xlim[1], 90)
    yedges = np.linspace(0.0, 1.0, 70)
    for mask, color, label in [
        (data["is_signal"].astype(bool), RED, "truth photons"),
        (~data["is_signal"].astype(bool), BLUE, "inclusive jets"),
    ]:
        h, _, _ = np.histogram2d(data["eiso"][mask], data["score"][mask], bins=[xedges, yedges])
        h = smooth_hist(h)
        h = h / max(float(h.sum()), 1.0)
        levels = np.nanpercentile(h[h > 0], [55, 75, 88, 96])
        xx = 0.5 * (xedges[:-1] + xedges[1:])
        yy = 0.5 * (yedges[:-1] + yedges[1:])
        ax.contour(xx, yy, h.T, levels=np.unique(levels), colors=color, linewidths=1.6)
        ax.plot([], [], color=color, lw=2.5, label=label)
    setup_ax(ax, "Class-normalized density contours", xlim)
    ax.legend(frameon=False, loc="upper right")
    out = OUT_DIR / "03_class_density_contours.png"
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)
    return out


def save_quantile_ribbons(data: dict[str, np.ndarray], xlim: tuple[float, float]) -> Path:
    fig, ax = plt.subplots(figsize=(8.4, 6.2), dpi=180)
    bins = np.linspace(xlim[0], xlim[1], 24)
    for mask, color, label in [
        (np.ones_like(data["score"], dtype=bool), INK, "all candidates"),
        (data["is_signal"].astype(bool), RED, "truth photons"),
        (~data["is_signal"].astype(bool), BLUE, "inclusive jets"),
    ]:
        xc, q16, q25, q50, q75, q84 = quantiles_by_bin(data["eiso"][mask], data["score"][mask], bins)
        ax.fill_between(xc, q16, q84, color=color, alpha=0.10, lw=0)
        ax.fill_between(xc, q25, q75, color=color, alpha=0.18, lw=0)
        ax.plot(xc, q50, color=color, lw=2.3, label=label)
    setup_ax(ax, "Binned quantile ribbons: score trend without point clutter", xlim)
    ax.legend(frameon=False, loc="upper right")
    out = OUT_DIR / "04_binned_quantile_ribbons.png"
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)
    return out


def save_centrality_facets(data: dict[str, np.ndarray], xlim: tuple[float, float]) -> Path:
    fig, axes = plt.subplots(1, 3, figsize=(12.2, 4.3), dpi=180, sharex=True, sharey=True)
    for ax, (lo, hi, label) in zip(axes, CENT_BINS):
        mask = (data["cent"] >= lo) & (data["cent"] < hi)
        h = ax.hist2d(
            data["eiso"][mask],
            data["score"][mask],
            bins=(70, 55),
            range=[xlim, [0, 1]],
            norm=LogNorm(vmin=1),
            cmap="viridis",
        )
        setup_ax(ax, label, xlim)
        ax.set_xlabel(r"$reco\_eiso$ [GeV]", fontsize=11)
    axes[0].set_ylabel("BDT score", fontsize=12)
    fig.colorbar(h[3], ax=axes.ravel().tolist(), pad=0.012, label="candidates / bin")
    fig.suptitle("Centrality facets: same relationship, split by event activity", fontsize=15, weight="bold")
    out = OUT_DIR / "05_centrality_faceted_log_density.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


def save_pass_fraction(data: dict[str, np.ndarray], xlim: tuple[float, float]) -> Path:
    sig = data["is_signal"].astype(bool)
    threshold = float(np.nanpercentile(data["score"][sig], 20.0))
    bins = np.linspace(xlim[0], xlim[1], 26)
    fig, ax = plt.subplots(figsize=(8.4, 5.5), dpi=180)
    for mask, color, label in [(sig, RED, "truth photons"), (~sig, BLUE, "inclusive jets")]:
        xs, vals, errs = [], [], []
        for lo, hi in zip(bins[:-1], bins[1:]):
            b = mask & (data["eiso"] >= lo) & (data["eiso"] < hi)
            n = int(b.sum())
            if n < 25:
                continue
            p = float((data["score"][b] >= threshold).mean())
            xs.append(0.5 * (lo + hi))
            vals.append(p)
            errs.append((p * (1.0 - p) / n) ** 0.5)
        ax.errorbar(xs, vals, yerr=errs, color=color, marker="o", ms=4, lw=2, capsize=2, label=label)
    ax.axhline(0.80, color=RED, ls="--", lw=1.2, alpha=0.65)
    ax.set_title("BDT pass fraction vs default isolation", fontsize=15, weight="bold")
    ax.set_xlabel(r"default $E_T^{iso}$ ($reco\_eiso$) [GeV]", fontsize=12)
    ax.set_ylabel(f"fraction with score >= {threshold:.3f}", fontsize=12)
    ax.set_xlim(*xlim)
    ax.set_ylim(0, 1.02)
    ax.grid(True, color="#dfe4ea", lw=0.7)
    ax.legend(frameon=False, loc="upper right")
    out = OUT_DIR / "06_bdt_pass_fraction_vs_default_eiso.png"
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)
    return out


def save_signal_probability(data: dict[str, np.ndarray], xlim: tuple[float, float]) -> Path:
    sig = data["is_signal"].astype(bool)
    xedges = np.linspace(xlim[0], xlim[1], 70)
    yedges = np.linspace(0.0, 1.0, 58)
    all_h, _, _ = np.histogram2d(data["eiso"], data["score"], bins=[xedges, yedges])
    sig_h, _, _ = np.histogram2d(data["eiso"][sig], data["score"][sig], bins=[xedges, yedges])
    prob = np.divide(sig_h, all_h, out=np.full_like(sig_h, np.nan, dtype="float64"), where=all_h >= 8)
    fig, ax = plt.subplots(figsize=(8.4, 6.2), dpi=180)
    mesh = ax.pcolormesh(xedges, yedges, prob.T, cmap="RdYlBu_r", vmin=0.0, vmax=1.0, shading="auto")
    plt.colorbar(mesh, ax=ax, pad=0.012, label="truth-photon fraction in bin")
    setup_ax(ax, "Signal-probability map: where the classes separate", xlim)
    out = OUT_DIR / "07_signal_probability_map.png"
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)
    return out


def save_side_by_side_class_density(data: dict[str, np.ndarray], xlim: tuple[float, float]) -> Path:
    sig = data["is_signal"].astype(bool)
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.8), dpi=180, sharex=True, sharey=True)
    last = None
    for ax, mask, title, cmap in [
        (axes[0], sig, "Truth photons", "Reds"),
        (axes[1], ~sig, "Inclusive jets", "Blues"),
    ]:
        last = ax.hist2d(
            data["eiso"][mask],
            data["score"][mask],
            bins=(70, 55),
            range=[xlim, [0, 1]],
            norm=LogNorm(vmin=1),
            cmap=cmap,
        )
        setup_ax(ax, title, xlim)
        ax.set_xlabel(r"$reco\_eiso$ [GeV]", fontsize=11)
    axes[0].set_ylabel("BDT score", fontsize=12)
    fig.colorbar(last[3], ax=axes.ravel().tolist(), pad=0.012, label="candidates / bin")
    fig.suptitle("Side-by-side class densities: avoid overplotting", fontsize=15, weight="bold")
    out = OUT_DIR / "08_side_by_side_class_densities.png"
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


def save_score_distributions_by_eiso_quantile(data: dict[str, np.ndarray]) -> Path:
    q = np.nanpercentile(data["eiso"], [0, 20, 40, 60, 80, 100])
    score_bins = np.linspace(0.0, 1.0, 55)
    centers = 0.5 * (score_bins[:-1] + score_bins[1:])
    colors = ["#2c7bb6", "#00a6ca", "#00ccbc", "#fdae61", "#d7191c"]
    fig, ax = plt.subplots(figsize=(8.4, 5.7), dpi=180)
    for i, (lo, hi) in enumerate(zip(q[:-1], q[1:])):
        mask = (data["eiso"] >= lo) & (data["eiso"] <= hi if i == 4 else data["eiso"] < hi)
        hist, _ = np.histogram(data["score"][mask], bins=score_bins, density=True)
        label = f"q{i + 1}: {lo:.1f} to {hi:.1f} GeV"
        ax.plot(centers, hist, lw=2.2, color=colors[i], label=label)
    ax.set_title("Score shape by isolation quantile", fontsize=15, weight="bold")
    ax.set_xlabel("BDT score", fontsize=12)
    ax.set_ylabel("normalized candidates", fontsize=12)
    ax.grid(True, color="#dfe4ea", lw=0.7)
    ax.legend(frameon=False, fontsize=9.5)
    out = OUT_DIR / "09_score_distribution_by_isolation_quantile.png"
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)
    return out


def save_median_score_cent_eiso(data: dict[str, np.ndarray], xlim: tuple[float, float]) -> Path:
    xedges = np.linspace(xlim[0], xlim[1], 42)
    yedges = np.linspace(0.0, 80.0, 17)
    med = np.full((len(xedges) - 1, len(yedges) - 1), np.nan)
    for ix, (x0, x1) in enumerate(zip(xedges[:-1], xedges[1:])):
        for iy, (y0, y1) in enumerate(zip(yedges[:-1], yedges[1:])):
            mask = (data["eiso"] >= x0) & (data["eiso"] < x1) & (data["cent"] >= y0) & (data["cent"] < y1)
            if mask.sum() >= 20:
                med[ix, iy] = np.nanmedian(data["score"][mask])
    fig, ax = plt.subplots(figsize=(8.4, 5.7), dpi=180)
    mesh = ax.pcolormesh(xedges, yedges, med.T, cmap="viridis", vmin=0.25, vmax=0.75, shading="auto")
    plt.colorbar(mesh, ax=ax, pad=0.012, label="median BDT score")
    ax.set_title("Median score across isolation and centrality", fontsize=15, weight="bold")
    ax.set_xlabel(r"default $E_T^{iso}$ ($reco\_eiso$) [GeV]", fontsize=12)
    ax.set_ylabel("centrality [%]", fontsize=12)
    ax.grid(False)
    out = OUT_DIR / "10_median_score_eiso_centrality_heatmap.png"
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)
    return out


def make_contact_sheet(paths: list[Path]) -> Path:
    thumbs = []
    for p in paths:
        img = Image.open(p).convert("RGB")
        img.thumbnail((760, 520), Image.Resampling.LANCZOS)
        canvas = Image.new("RGB", (800, 590), "white")
        canvas.paste(img, ((800 - img.width) // 2, 28))
        draw = ImageDraw.Draw(canvas)
        try:
            font = ImageFont.truetype("Arial.ttf", 22)
        except OSError:
            font = ImageFont.load_default()
        draw.text((24, 555), p.name, fill=(20, 20, 20), font=font)
        thumbs.append(canvas)
    rows = (len(thumbs) + 1) // 2
    sheet = Image.new("RGB", (1600, 590 * rows), (246, 246, 246))
    for i, img in enumerate(thumbs):
        sheet.paste(img, ((i % 2) * 800, (i // 2) * 590))
    out = OUT_DIR / "00_contact_sheet_default_eiso_bdt_variants.png"
    sheet.save(out)
    return out


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    data = load_arrays()
    xlim = tuple(float(v) for v in np.nanpercentile(data["eiso"], [0.5, 99.5]))
    paths = [
        save_log_density(data, xlim),
        save_rank_density(data),
        save_contours(data, xlim),
        save_quantile_ribbons(data, xlim),
        save_centrality_facets(data, xlim),
        save_pass_fraction(data, xlim),
        save_signal_probability(data, xlim),
        save_side_by_side_class_density(data, xlim),
        save_score_distributions_by_eiso_quantile(data),
        save_median_score_cent_eiso(data, xlim),
    ]
    sheet = make_contact_sheet(paths)
    manifest = OUT_DIR / "plot_variants_manifest.json"
    manifest.write_text(
        json.dumps(
            {
                "source_report": str(REPORT),
                "manifest": str(MANIFEST),
                "score": SCORE,
                "eiso": EISO,
                "selected_entries": int(len(data["score"])),
                "signal_entries": int(data["is_signal"].sum()),
                "background_entries": int((~data["is_signal"].astype(bool)).sum()),
                "xlim_percentile_0p5_99p5": list(xlim),
                "outputs": [str(p) for p in [sheet, *paths]],
            },
            indent=2,
        )
        + "\n"
    )
    print(f"wrote {sheet}")
    for path in paths:
        print(f"wrote {path}")
    print(f"wrote {manifest}")


if __name__ == "__main__":
    main()
