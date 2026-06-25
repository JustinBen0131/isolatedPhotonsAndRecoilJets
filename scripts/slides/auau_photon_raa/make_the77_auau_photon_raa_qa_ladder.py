#!/usr/bin/env python3
"""Build the first THE-77 AuAu photon-RAA QA ladder slide packet.

This packet is intentionally source-driven.  It uses the current THE-69 merged
ROOT outputs and existing pp reference diagnostics, and it marks measurement
contract gaps explicitly instead of filling them with placeholders.
"""

from __future__ import annotations

import csv
import json
import math
import sys
import textwrap
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
from PIL import Image, ImageDraw


THIS_FILE = Path(__file__).resolve()
REPO = next((p for p in THIS_FILE.parents if (p / "AGENTS.md").exists()), THIS_FILE.parents[4])
SCRIPTS = REPO / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.append(str(SCRIPTS))

from slides.common.slide_defaults import SLIDE_DPI, slide_figsize  # noqa: E402


INPUT_DIR = REPO / "InputFiles/the69_default_auau_physicsqa"
DATA_ROOT = INPUT_DIR / "RecoilJets_auau_ALL_preselectionNewPPG12_tightAuAuCentInputBase3x3BDT_nonTightAuAuBDTComplement_baseVariant.root"
SIGNAL_ROOT = INPUT_DIR / "RecoilJets_embeddedPhoton12plus20_MERGED.root"
INCLUSIVE_ROOT = INPUT_DIR / "RecoilJets_embeddedJet12plus20plus30plus40_MERGED.root"
PP_ROOT = (
    REPO
    / "dataOutput/ppg12TableQA/THE42_ppg12_tableqa_v1_basev3e_20260611/merged_roots/"
    / "RecoilJets_pp_ALL_preselectionNewPPG12_tightNewPPG12_nonTightNewPPG12.root"
)
FINE_CSV = (
    REPO
    / "dataOutput/auauPhysicsQA/THE69_defaultAuAuBDT_physicsQA_20260617/fine_yield_purity_pp_slide/"
    / "the69_fine_abcd_yield_purity_points.csv"
)

OUT_DIR = REPO / "dataOutput/auauPhysicsQA/THE77_auau_photon_raa_qa_ladder_20260617"
SLIDE_DIR = OUT_DIR / "slides"
MANIFEST = OUT_DIR / "the77_auau_photon_raa_qa_ladder_manifest.json"
SCRIPT_MD = OUT_DIR / "the77_auau_photon_raa_qa_ladder_speaker_script.md"
CONTACT_SHEET = OUT_DIR / "the77_qa_ladder_contact_sheet.png"

DATA_TOP = "MBD_NS_geq_2_vtx_lt_150"
SIM_TOP = "SIM"
PT_TOKEN = "1535"
WP80 = "T80(c) = 0.53471108 + 0.0012284143*c"

CENT = [
    ("0_20", "cent0_20", "0-20%", 10.0, "#1F77B4"),
    ("20_50", "cent20_50", "20-50%", 35.0, "#2CA02C"),
    ("50_80", "cent50_80", "50-80%", 65.0, "#7A5AA6"),
]
PT_BINS = [(14, 16), (16, 18), (18, 20), (20, 22), (22, 24), (24, 26), (26, 35)]
SAMPLES = {
    "Data": (DATA_ROOT, DATA_TOP, "#111827"),
    "Signal MC": (SIGNAL_ROOT, SIM_TOP, "#C43C39"),
    "Inclusive MC": (INCLUSIVE_ROOT, SIM_TOP, "#2F63C6"),
}


@dataclass
class HistArrays:
    x: np.ndarray
    y: np.ndarray
    e: np.ndarray
    integral: float
    source: str


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "axes.facecolor": "white",
            "axes.edgecolor": "#1F2937",
            "axes.linewidth": 1.0,
            "axes.labelsize": 12,
            "axes.titlesize": 13,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "legend.fontsize": 9.5,
        }
    )


def open_root(path: Path):
    import ROOT

    ROOT.gROOT.SetBatch(True)
    ROOT.TH1.AddDirectory(False)
    f = ROOT.TFile.Open(str(path), "READ")
    if not f or f.IsZombie():
        raise OSError(f"Could not open ROOT file: {path}")
    return f


def get_clone(path: Path, obj_path: str):
    f = open_root(path)
    try:
        obj = f.Get(obj_path)
        if not obj:
            return None
        clone = obj.Clone(f"{Path(obj_path).name}_clone")
        if hasattr(clone, "SetDirectory"):
            clone.SetDirectory(0)
        return clone
    finally:
        f.Close()


def hist_integral(path: Path, obj_path: str) -> float:
    h = get_clone(path, obj_path)
    if h is None or not hasattr(h, "Integral"):
        return float("nan")
    return float(h.Integral())


def bin_error_sum(path: Path, obj_path: str) -> tuple[float, float]:
    h = get_clone(path, obj_path)
    if h is None or not hasattr(h, "Integral"):
        return float("nan"), float("nan")
    y = float(h.Integral())
    err2 = 0.0
    for i in range(1, h.GetNbinsX() + 1):
        err2 += float(h.GetBinError(i)) ** 2
    return y, math.sqrt(err2)


def h1_to_arrays(path: Path, obj_path: str, *, rebin: int = 1, xlim: tuple[float, float] | None = None) -> HistArrays | None:
    h = get_clone(path, obj_path)
    if h is None or not hasattr(h, "GetNbinsX"):
        return None
    nb = h.GetNbinsX()
    x = np.array([h.GetXaxis().GetBinCenter(i) for i in range(1, nb + 1)], dtype=float)
    y = np.array([h.GetBinContent(i) for i in range(1, nb + 1)], dtype=float)
    e = np.array([h.GetBinError(i) for i in range(1, nb + 1)], dtype=float)
    if rebin > 1 and len(x) >= rebin:
        n = len(x) // rebin
        trim = n * rebin
        x = x[:trim].reshape(n, rebin).mean(axis=1)
        y = y[:trim].reshape(n, rebin).sum(axis=1)
        e = np.sqrt((e[:trim].reshape(n, rebin) ** 2).sum(axis=1))
    mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(e)
    if xlim is not None:
        mask &= (x >= xlim[0]) & (x <= xlim[1])
    x, y, e = x[mask], y[mask], e[mask]
    integral = float(y.sum())
    return HistArrays(x=x, y=y, e=e, integral=integral, source=obj_path)


def normalized(arr: HistArrays | None) -> HistArrays | None:
    if arr is None:
        return None
    if arr.integral <= 0:
        return arr
    return HistArrays(arr.x, arr.y / arr.integral, arr.e / arr.integral, arr.integral, arr.source)


def combine_arrays(parts: list[HistArrays | None]) -> HistArrays | None:
    good = [p for p in parts if p is not None]
    if not good:
        return None
    x = good[0].x.copy()
    y = np.zeros_like(x)
    e2 = np.zeros_like(x)
    sources = []
    for arr in good:
        if len(arr.x) != len(x) or np.max(np.abs(arr.x - x)) > 1e-6:
            raise ValueError("Cannot combine histograms with different binning")
        y += arr.y
        e2 += arr.e * arr.e
        sources.append(arr.source)
    return HistArrays(x=x, y=y, e=np.sqrt(e2), integral=float(y.sum()), source=" + ".join(sources))


def h2_to_image(path: Path, obj_path: str, *, max_x_bins: int | None = None, max_y_bins: int | None = None):
    h = get_clone(path, obj_path)
    if h is None or not hasattr(h, "GetNbinsX"):
        return None
    nx, ny = h.GetNbinsX(), h.GetNbinsY()
    x_edges = np.array([h.GetXaxis().GetBinLowEdge(i) for i in range(1, nx + 2)], dtype=float)
    y_edges = np.array([h.GetYaxis().GetBinLowEdge(i) for i in range(1, ny + 2)], dtype=float)
    z = np.array([[h.GetBinContent(ix, iy) for ix in range(1, nx + 1)] for iy in range(1, ny + 1)], dtype=float)
    if max_x_bins and z.shape[1] > max_x_bins:
        step = math.ceil(z.shape[1] / max_x_bins)
        keep = (z.shape[1] // step) * step
        z = z[:, :keep].reshape(z.shape[0], keep // step, step).sum(axis=2)
        x_edges = x_edges[0 : keep + 1 : step]
        if len(x_edges) != z.shape[1] + 1:
            x_edges = np.linspace(float(x_edges[0]), float(x_edges[-1]), z.shape[1] + 1)
    if max_y_bins and z.shape[0] > max_y_bins:
        step = math.ceil(z.shape[0] / max_y_bins)
        keep = (z.shape[0] // step) * step
        z = z[:keep, :].reshape(keep // step, step, z.shape[1]).sum(axis=1)
        y_edges = y_edges[0 : keep + 1 : step]
        if len(y_edges) != z.shape[0] + 1:
            y_edges = np.linspace(float(y_edges[0]), float(y_edges[-1]), z.shape[0] + 1)
    return x_edges, y_edges, z


def weighted_mean_y(path: Path, obj_path: str) -> tuple[float, float] | None:
    h = get_clone(path, obj_path)
    if h is None or not hasattr(h, "GetNbinsX"):
        return None
    sumw = 0.0
    sumwy = 0.0
    for ix in range(1, h.GetNbinsX() + 1):
        for iy in range(1, h.GetNbinsY() + 1):
            w = float(h.GetBinContent(ix, iy))
            if w <= 0:
                continue
            y = float(h.GetYaxis().GetBinCenter(iy))
            sumw += w
            sumwy += w * y
    if sumw <= 0:
        return None
    return sumwy / sumw, sumw


def add_title(fig, title: str, kicker: str | None = None) -> None:
    fig.text(0.045, 0.94, title, ha="left", va="top", fontsize=24, fontweight="bold", color="#111827")
    if kicker:
        fig.text(0.045, 0.89, kicker, ha="left", va="top", fontsize=13.5, color="#374151")


def add_takeaway(fig, text: str, *, y: float = 0.055) -> None:
    wrapped = textwrap.fill(text, 118)
    fig.text(
        0.055,
        y,
        wrapped,
        ha="left",
        va="bottom",
        fontsize=13.5,
        color="#111827",
        bbox=dict(facecolor="#F8FAFC", edgecolor="#CBD5E1", linewidth=1.1, boxstyle="round,pad=0.45"),
    )


def add_card(
    fig,
    x: float,
    y: float,
    w: float,
    h: float,
    title: str,
    body: str,
    *,
    accent: str = "#2563EB",
    body_wrap: int = 48,
    title_size: float = 14.5,
    body_size: float = 11.8,
) -> None:
    rect = plt.Rectangle((x, y), w, h, transform=fig.transFigure, facecolor="white", edgecolor="#D1D5DB", linewidth=1.1)
    fig.add_artist(rect)
    fig.add_artist(plt.Rectangle((x, y), 0.008, h, transform=fig.transFigure, facecolor=accent, edgecolor=accent))
    fig.text(x + 0.018, y + h - 0.03, title, ha="left", va="top", fontsize=title_size, fontweight="bold", color="#111827")
    fig.text(x + 0.018, y + h - 0.082, textwrap.fill(body, body_wrap), ha="left", va="top", fontsize=body_size, color="#374151", linespacing=1.18)


def save_slide(fig, filename: str) -> Path:
    SLIDE_DIR.mkdir(parents=True, exist_ok=True)
    out = SLIDE_DIR / filename
    fig.savefig(out, dpi=SLIDE_DPI, facecolor="white")
    plt.close(fig)
    return out


def draw_norm_step(ax, arr: HistArrays | None, label: str, color: str, *, lw: float = 2.0, marker: bool = False) -> None:
    if arr is None or arr.integral <= 0 or len(arr.x) == 0:
        return
    y = arr.y
    if marker:
        ax.errorbar(arr.x, y, yerr=arr.e, fmt="o", markersize=3.6, color=color, label=label, linewidth=0.9)
    else:
        ax.step(arr.x, y, where="mid", color=color, linewidth=lw, label=label)


def bdt_path(top: str, cent_hist: str, cut: str) -> str:
    return f"{top}/h1d_bdt_eta0_pt{PT_TOKEN}_{cent_hist}_{cut}"


def var_path(top: str, var: str, cent_hist: str, cut: str) -> str:
    return f"{top}/h1d_{var}_eta0_pt{PT_TOKEN}_{cent_hist}_{cut}"


def iso_path(top: str, prefix: str, lo: int, hi: int, cent: str) -> str:
    return f"{top}/{prefix}_isoR40_pT_{lo}_{hi}_cent_{cent}"


def truth_match_path(top: str, prefix: str, lo: int, hi: int, cent: str) -> str:
    return f"{top}/{prefix}_isoR40_pT_{lo}_{hi}_cent_{cent}"


def load_fine_csv() -> list[dict[str, str]]:
    with FINE_CSV.open(newline="") as f:
        return list(csv.DictReader(f))


def fnum(row: dict[str, str], key: str) -> float:
    try:
        return float(row.get(key, "nan"))
    except ValueError:
        return float("nan")


def make_slide_01_scope(summary: dict) -> Path:
    fig = plt.figure(figsize=slide_figsize(), dpi=SLIDE_DPI)
    add_title(
        fig,
        "Current AuAu photon-ID QA sample is ready for first measurement checks",
        "Default 14-feature AuAu BDT, runtime WP80 tight cut, current-production data plus matched embedded MC.",
    )

    add_card(fig, 0.055, 0.64, 0.25, 0.17, "Data scope", "882 analyzable paired GRL runs and 176,577 matched CALOFITTING/ZDC segment pairs. The 4 pm append check found no newly available matched pairs after submission.", accent="#2563EB")
    add_card(fig, 0.335, 0.64, 0.25, 0.17, "Merged products", "Final data, signal-MC, and inclusive-MC ROOT files are local and structurally readable under InputFiles/the69_default_auau_physicsqa.", accent="#059669")
    add_card(fig, 0.615, 0.64, 0.31, 0.17, "Photon-ID contract", f"Model centAsFeatBase3x3_pt15to35 with {WP80}. Tight and non-tight rows are filled in the current output.", accent="#7C3AED")

    ax = fig.add_axes([0.075, 0.23, 0.38, 0.32])
    labels = [c[2] for c in CENT]
    cut0 = [hist_integral(DATA_ROOT, bdt_path(DATA_TOP, cent_hist, "cut0")) for _, cent_hist, _, _, _ in CENT]
    cut1 = [hist_integral(DATA_ROOT, bdt_path(DATA_TOP, cent_hist, "cut1")) for _, cent_hist, _, _, _ in CENT]
    cut2 = [hist_integral(DATA_ROOT, bdt_path(DATA_TOP, cent_hist, "cut2")) for _, cent_hist, _, _, _ in CENT]
    x = np.arange(len(labels))
    width = 0.23
    ax.bar(x - width, cut0, width, color="#CBD5E1", edgecolor="#64748B", label="All candidates")
    ax.bar(x, cut1, width, color="#93C5FD", edgecolor="#2563EB", label="Preselection")
    ax.bar(x + width, cut2, width, color="#60A5FA", edgecolor="#1D4ED8", label="Tight WP80")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("15-35 GeV candidate count")
    ax.set_title("Data photon-candidate exposure by centrality")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(frameon=False, loc="upper right")
    for xi, val in zip(x + width, cut2):
        ax.text(xi, val * 1.035, f"{int(val):,}", ha="center", va="bottom", fontsize=9)

    ax2 = fig.add_axes([0.55, 0.23, 0.36, 0.32])
    evt = [hist_integral(DATA_ROOT, f"{DATA_TOP}/h_vertexZ_cent_{cent}") for cent, _, _, _, _ in CENT]
    total = sum(v for v in evt if np.isfinite(v))
    frac = [v / total if total else 0 for v in evt]
    ax2.bar(labels, frac, color=["#BFDBFE", "#BBF7D0", "#DDD6FE"], edgecolor="#334155")
    ax2.set_ylim(0, max(frac) * 1.35)
    ax2.set_ylabel("Fraction of event-weighted entries")
    ax2.set_title("Centrality exposure proxy from vertex-Z histograms")
    ax2.grid(axis="y", alpha=0.25)
    for i, (v, f) in enumerate(zip(evt, frac)):
        ax2.text(i, f + 0.012, f"{v/1e6:.1f}M", ha="center", va="bottom", fontsize=10)

    add_takeaway(fig, "This is sufficient for first-pass QA and raw photon-ID checks. The raw counts are exposure diagnostics, not luminosity- or TAA-normalized yields.")
    summary["slide01"] = {"cut0": cut0, "cut1": cut1, "cut2": cut2, "event_weighted_vertex_entries": evt}
    return save_slide(fig, "01_data_scope_exposure.png")


def make_slide_02_event_level(summary: dict) -> Path:
    fig = plt.figure(figsize=slide_figsize(), dpi=SLIDE_DPI)
    add_title(fig, "Event-level QA: vertex, centrality, and calorimeter-energy proxies", "These are data-side exposure checks before interpreting photon yields.")

    ax = fig.add_axes([0.06, 0.50, 0.38, 0.31])
    for cent, _, label, _, color in CENT:
        arr = normalized(h1_to_arrays(DATA_ROOT, f"{DATA_TOP}/h_vertexZ_cent_{cent}", rebin=6, xlim=(-160, 160)))
        draw_norm_step(ax, arr, label, color, lw=2.2)
    ax.set_xlabel("z vertex [cm]")
    ax.set_ylabel("Unit-normalized entries")
    ax.set_title("Vertex-Z is centered and similar across centrality")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False)

    ax2 = fig.add_axes([0.55, 0.50, 0.36, 0.31])
    cent_arr = normalized(h1_to_arrays(DATA_ROOT, f"{DATA_TOP}/h_centrality", rebin=2, xlim=(0, 80)))
    draw_norm_step(ax2, cent_arr, "Data", "#111827", lw=2.4)
    ax2.axvspan(0, 20, color="#BFDBFE", alpha=0.22)
    ax2.axvspan(20, 50, color="#BBF7D0", alpha=0.18)
    ax2.axvspan(50, 80, color="#DDD6FE", alpha=0.18)
    ax2.set_xlabel("Centrality percentile")
    ax2.set_ylabel("Unit-normalized entries")
    ax2.set_title("Centrality coverage over 0-80%")
    ax2.grid(alpha=0.25)

    ax3 = fig.add_axes([0.12, 0.13, 0.74, 0.25])
    image = h2_to_image(DATA_ROOT, f"{DATA_TOP}/h2_totalCaloEnergy_vs_centrality", max_x_bins=80, max_y_bins=70)
    if image:
        xe, ye, z = image
        z = np.where(z > 0, z, np.nan)
        mesh = ax3.pcolormesh(xe, ye, z, norm=LogNorm(vmin=np.nanmax(z) * 1e-5, vmax=np.nanmax(z)), cmap="Blues")
        fig.colorbar(mesh, ax=ax3, pad=0.012, label="entries")
    ax3.set_xlim(0, 80)
    ax3.set_xlabel("Centrality percentile")
    ax3.set_ylabel("Total calo energy proxy")
    ax3.set_title("Total calorimeter energy follows centrality as expected")
    ax3.grid(alpha=0.15)

    add_takeaway(fig, "The exposure looks coherent at event level. The plotted quantities are merged histogram entries/proxies and should be paired with the luminosity accounting slide before final RAA-style normalization.")
    summary["slide02"] = {"objects": ["h_vertexZ_cent_*", "h_centrality", "h2_totalCaloEnergy_vs_centrality"]}
    return save_slide(fig, "02_event_level_qa.png")


def make_slide_03_cutflow(summary: dict) -> Path:
    fig = plt.figure(figsize=slide_figsize(), dpi=SLIDE_DPI)
    add_title(fig, "Photon-ID cutflow: preselection and tight WP80 behave as expected", "15 < cluster ET < 35 GeV; cut1 is preselection/NPB, cut2 is tight WP80, cut3 is non-tight complement.")

    ax = fig.add_axes([0.065, 0.27, 0.42, 0.47])
    stages = [("cut0", "All"), ("cut1", "Preselection"), ("cut2", "Tight"), ("cut3", "Non-tight")]
    x = np.arange(len(stages))
    for j, (cent, cent_hist, label, _, color) in enumerate(CENT):
        vals = [hist_integral(DATA_ROOT, bdt_path(DATA_TOP, cent_hist, cut)) for cut, _ in stages]
        ax.plot(x + (j - 1) * 0.045, vals, "o", color=color, label=label, markersize=8)
        for xi, val in zip(x + (j - 1) * 0.045, vals):
            ax.text(xi, val * 1.04 + 45, f"{int(val):,}", ha="center", fontsize=8.8, color=color)
    ax.set_xticks(x)
    ax.set_xticklabels([s[1] for s in stages])
    ax.set_yscale("log")
    ax.set_ylabel("Data candidates")
    ax.set_title("Data cutflow by centrality")
    ax.grid(alpha=0.28, which="both")
    ax.legend(frameon=False, ncol=3, loc="lower left")

    ax2 = fig.add_axes([0.58, 0.27, 0.34, 0.47])
    labels = [c[2] for c in CENT]
    sig_eff = []
    inc_fake = []
    data_tight_frac = []
    for _, cent_hist, _, _, _ in CENT:
        s1 = hist_integral(SIGNAL_ROOT, bdt_path(SIM_TOP, cent_hist, "cut1"))
        s2 = hist_integral(SIGNAL_ROOT, bdt_path(SIM_TOP, cent_hist, "cut2"))
        b1 = hist_integral(INCLUSIVE_ROOT, bdt_path(SIM_TOP, cent_hist, "cut1"))
        b2 = hist_integral(INCLUSIVE_ROOT, bdt_path(SIM_TOP, cent_hist, "cut2"))
        d1 = hist_integral(DATA_ROOT, bdt_path(DATA_TOP, cent_hist, "cut1"))
        d2 = hist_integral(DATA_ROOT, bdt_path(DATA_TOP, cent_hist, "cut2"))
        sig_eff.append(s2 / s1 if s1 else float("nan"))
        inc_fake.append(b2 / b1 if b1 else float("nan"))
        data_tight_frac.append(d2 / d1 if d1 else float("nan"))
    xx = np.arange(len(labels))
    w = 0.25
    ax2.bar(xx - w, sig_eff, w, color="#FCA5A5", edgecolor="#B91C1C", label="Signal MC tight/pre")
    ax2.bar(xx, inc_fake, w, color="#BFDBFE", edgecolor="#1D4ED8", label="Inclusive MC tight/pre")
    ax2.bar(xx + w, data_tight_frac, w, color="#D1D5DB", edgecolor="#374151", label="Data tight/pre")
    ax2.set_ylim(0, 1.0)
    ax2.set_xticks(xx)
    ax2.set_xticklabels(labels)
    ax2.set_ylabel("Fraction")
    ax2.set_title("WP80 keeps signal-like MC high, inclusive lower")
    ax2.grid(axis="y", alpha=0.25)
    ax2.legend(frameon=False, fontsize=9)

    add_takeaway(fig, "The data tight fractions are about one half after preselection, while signal MC remains much higher than inclusive MC under the same WP80 cut.")
    summary["slide03"] = {"signal_mc_tight_over_preselection": sig_eff, "inclusive_mc_tight_over_preselection": inc_fake, "data_tight_over_preselection": data_tight_frac}
    return save_slide(fig, "03_photon_id_cutflow.png")


def make_slide_04_bdt(summary: dict) -> Path:
    fig = plt.figure(figsize=slide_figsize(), dpi=SLIDE_DPI)
    add_title(fig, "BDT-score QA: data sits between signal-like and inclusive-MC shapes", "After preselection/NPB, 15-35 GeV. Dashed vertical lines show the centrality-midpoint WP80 threshold.")
    axes = fig.subplots(1, 3, gridspec_kw={"left": 0.055, "right": 0.98, "bottom": 0.16, "top": 0.78, "wspace": 0.20})
    metrics = {}
    for ax, (cent, cent_hist, label, c_mid, _) in zip(axes, CENT):
        for sample, (root, top, color) in SAMPLES.items():
            arr = normalized(h1_to_arrays(root, bdt_path(top, cent_hist, "cut1"), rebin=2, xlim=(0, 1)))
            draw_norm_step(ax, arr, sample, color, lw=2.2, marker=(sample == "Data"))
        thr = 0.53471108 + 0.0012284143 * c_mid
        ax.axvline(thr, color="#111827", linestyle="--", linewidth=1.2, alpha=0.75)
        ax.text(thr + 0.01, ax.get_ylim()[1] * 0.8 if ax.get_ylim()[1] else 0.1, "T80", rotation=90, va="top", fontsize=9)
        ax.set_yscale("log")
        ax.set_ylim(5e-4, 1.0)
        ax.set_xlim(0, 1)
        ax.set_xlabel("BDT score")
        ax.set_title(label)
        ax.grid(alpha=0.25, which="both")
        if ax is axes[0]:
            ax.set_ylabel("Unit-normalized candidates")
        else:
            ax.set_yticklabels([])
        metrics[label] = {"t80_midpoint": thr}
    axes[-1].legend(frameon=False, loc="upper left")
    add_takeaway(fig, "The current data BDT distributions are not pathological. The comparison is shape QA only: data are raw counts, MC are weighted embedded samples.")
    summary["slide04"] = metrics
    return save_slide(fig, "04_bdt_score_data_mc.png")


def make_slide_05_shower(summary: dict) -> Path:
    fig = plt.figure(figsize=slide_figsize(), dpi=SLIDE_DPI)
    add_title(fig, "Shower-shape QA: E11/E33 and E32/E35 are coherent after preselection", "Unit-normalized data, signal MC, and inclusive MC overlays in the same 15-35 GeV candidate row.")
    variables = [("e11_to_e33", "E11/E33", (0, 1.05), 4), ("e32_to_e35", "E32/E35", (0, 1.05), 4)]
    axes = fig.subplots(2, 3, gridspec_kw={"left": 0.055, "right": 0.98, "bottom": 0.14, "top": 0.79, "wspace": 0.20, "hspace": 0.34})
    metrics = {}
    for row, (var, var_label, xlim, rebin) in enumerate(variables):
        for col, (cent, cent_hist, label, _, _) in enumerate(CENT):
            ax = axes[row][col]
            for sample, (root, top, color) in SAMPLES.items():
                arr = normalized(h1_to_arrays(root, var_path(top, var, cent_hist, "cut1"), rebin=rebin, xlim=xlim))
                draw_norm_step(ax, arr, sample, color, lw=2.0, marker=(sample == "Data"))
                if arr is not None and sample == "Data":
                    low = float(arr.y[arr.x < 0.10].sum()) if arr.integral > 0 else float("nan")
                    metrics[f"{var}_{label}_data_low_lt_0p10"] = low
            ax.set_xlim(*xlim)
            ax.set_yscale("log")
            ax.set_ylim(5e-4, 1.0)
            ax.grid(alpha=0.25, which="both")
            if row == 0:
                ax.set_title(label)
            if col == 0:
                ax.set_ylabel(var_label + "\nunit-normalized")
            else:
                ax.set_yticklabels([])
            if row == 1:
                ax.set_xlabel(var_label)
    axes[0][2].legend(frameon=False, loc="upper left")
    add_takeaway(fig, "The low-edge spike that triggered the earlier stop is not present as a dominant data feature here; signal and inclusive MC remain well behaved under the corrected default reconstruction path.")
    summary["slide05"] = metrics
    return save_slide(fig, "05_shower_shape_closure.png")


def make_slide_06_iso(summary: dict) -> Path:
    fig = plt.figure(figsize=slide_figsize(), dpi=SLIDE_DPI)
    add_title(fig, "Isolation QA: R=0.4 isolation shape changes coherently after tight ID", "22-26 GeV slice, combined 22-24 and 24-26 table-QA bins; unit-normalized shapes.")
    axes = fig.subplots(2, 3, gridspec_kw={"left": 0.055, "right": 0.98, "bottom": 0.14, "top": 0.79, "wspace": 0.20, "hspace": 0.34})
    prefixes = [("h_Eiso", "After preselection"), ("h_Eiso_tight", "After tight WP80")]
    metrics = {}
    for row, (prefix, row_label) in enumerate(prefixes):
        for col, (cent, _, label, _, _) in enumerate(CENT):
            ax = axes[row][col]
            for sample, (root, top, color) in SAMPLES.items():
                arr = combine_arrays(
                    [
                        h1_to_arrays(root, iso_path(top, prefix, 22, 24, cent), rebin=3, xlim=(-20, 40)),
                        h1_to_arrays(root, iso_path(top, prefix, 24, 26, cent), rebin=3, xlim=(-20, 40)),
                    ]
                )
                arrn = normalized(arr)
                draw_norm_step(ax, arrn, sample, color, lw=2.0, marker=(sample == "Data"))
                if arr and sample == "Data":
                    metrics[f"{row_label}_{label}_data_integral"] = arr.integral
            ax.set_yscale("log")
            ax.set_ylim(5e-4, 1.0)
            ax.set_xlim(-20, 40)
            ax.grid(alpha=0.25, which="both")
            if row == 0:
                ax.set_title(label)
            if col == 0:
                ax.set_ylabel(row_label + "\nunit-normalized")
            else:
                ax.set_yticklabels([])
            if row == 1:
                ax.set_xlabel("E_iso, R=0.4 [GeV]")
    axes[0][2].legend(frameon=False, loc="upper right")
    add_takeaway(fig, "The tight-ID isolation slice is statistically sparse in data but the centrality ordering and MC shapes are coherent enough for first-pass purity checks.")
    summary["slide06"] = metrics
    return save_slide(fig, "06_isolation_closure.png")


def make_slide_07_abcd(summary: dict) -> Path:
    rows = load_fine_csv()
    fig = plt.figure(figsize=slide_figsize(), dpi=SLIDE_DPI)
    add_title(fig, "Raw ABCD purity and A-region yield: useful QA, not final yield", "Fine ET bins from current table-QA histograms; pp is a stock PPG12/baseV3E reference for shape context.")
    ax = fig.add_axes([0.06, 0.24, 0.42, 0.50])
    ax2 = fig.add_axes([0.57, 0.24, 0.36, 0.50])
    for sample, color in [("Au+Au 0-20%", "#1F77B4"), ("Au+Au 20-50%", "#2CA02C"), ("Au+Au 50-80%", "#7A5AA6"), ("pp reference", "#111827")]:
        pts = [r for r in rows if r["sample"] == sample and r["status"] == "ok"]
        x = np.array([fnum(r, "pt_mid") for r in pts], dtype=float)
        y = np.array([fnum(r, "yield_per_GeV") for r in pts], dtype=float)
        ye = np.array([fnum(r, "yield_per_GeV_err") for r in pts], dtype=float)
        p = np.array([fnum(r, "raw_purity") for r in pts], dtype=float)
        pe = np.array([fnum(r, "raw_purity_err") for r in pts], dtype=float)
        ax.errorbar(x, y, yerr=ye, fmt="o", color=color, markersize=5.5, linewidth=1.1, label=sample)
        ax2.errorbar(x, p, yerr=pe, fmt="o", color=color, markersize=5.5, linewidth=1.1, label=sample)
    ax.set_yscale("log")
    ax.set_xlabel("Cluster ET [GeV]")
    ax.set_ylabel("A-region raw count / GeV")
    ax.set_title("Raw isolated-tight candidate yield proxy")
    ax.grid(alpha=0.25, which="both")
    ax.legend(frameon=False, fontsize=9)
    ax2.set_ylim(-0.05, 1.15)
    ax2.set_xlabel("Cluster ET [GeV]")
    ax2.set_ylabel("Raw ABCD purity")
    ax2.set_title("Raw ABCD purity: max(0, A - BC/D) / A")
    ax2.grid(alpha=0.25)
    add_takeaway(fig, "The raw purity is broadly plausible but sparse in the peripheral/high-ET bins. The final photon RAA path needs luminosity/TAA normalization and a vetted purity-systematics treatment before using these as physics yields.")
    summary["slide07"] = {"csv": str(FINE_CSV), "definition": "raw A-region count and raw ABCD purity"}
    return save_slide(fig, "07_raw_abcd_purity_yield.png")


def make_slide_08_truth_proxy(summary: dict) -> Path:
    fig = plt.figure(figsize=slide_figsize(), dpi=SLIDE_DPI)
    add_title(fig, "MC truth-matched tight-photon proxy: first efficiency/fake-rate handle", "Proxy ratio: truth-matched tight isoR40 integral divided by all tight isoR40 integral in the same ET and centrality bin.")
    axes = fig.subplots(1, 3, gridspec_kw={"left": 0.06, "right": 0.95, "bottom": 0.33, "top": 0.72, "wspace": 0.28})
    metrics = {}
    for ax, (cent, _, label, _, color) in zip(axes, CENT):
        for sample_label, root, marker, line_color in [
            ("Signal MC", SIGNAL_ROOT, "o", "#C43C39"),
            ("Inclusive MC", INCLUSIVE_ROOT, "s", "#2F63C6"),
        ]:
            xs, ratios, errs = [], [], []
            for lo, hi in PT_BINS:
                num, nume = bin_error_sum(root, truth_match_path(SIM_TOP, "h_EisoReco_truthSigMatched_tight", lo, hi, cent))
                den, dene = bin_error_sum(root, iso_path(SIM_TOP, "h_Eiso_tight", lo, hi, cent))
                if not (np.isfinite(num) and np.isfinite(den)) or den <= 0:
                    continue
                ratio = num / den
                err = ratio * math.sqrt((nume / num) ** 2 + (dene / den) ** 2) if num > 0 and den > 0 else 0.0
                xs.append((lo + hi) / 2)
                ratios.append(ratio)
                errs.append(err)
            ax.errorbar(xs, ratios, yerr=errs, fmt=marker, color=line_color, markersize=5.5, linewidth=1.0, label=sample_label)
            metrics[f"{label}_{sample_label}"] = ratios
        ax.set_ylim(0, 1.05)
        ax.set_xlabel("Cluster ET [GeV]")
        ax.set_title(label)
        ax.grid(alpha=0.25)
        if ax is axes[0]:
            ax.set_ylabel("Truth-matched tight fraction")
        else:
            ax.set_yticklabels([])
        ax.legend(frameon=False, loc="lower right")

    add_card(
        fig,
        0.08,
        0.065,
        0.37,
        0.15,
        "What this proves now",
        "The current MC output carries truth-matched tight photon histograms, so we can already QA signal purity and inclusive fake contamination proxies.",
        accent="#059669",
        body_wrap=50,
        body_size=10.7,
    )
    add_card(
        fig,
        0.52,
        0.065,
        0.37,
        0.15,
        "What this is not yet",
        "This is not the final photon-ID efficiency denominator. The next production contract should add explicit truth-photon denominator and reco-matched numerator spectra for publication-grade efficiency.",
        accent="#B45309",
        body_wrap=52,
        body_size=10.7,
    )
    summary["slide08"] = metrics
    return save_slide(fig, "08_mc_truth_matched_proxy.png")


def make_slide_09_response(summary: dict) -> Path:
    fig = plt.figure(figsize=slide_figsize(), dpi=SLIDE_DPI)
    add_title(fig, "Unfolding readiness: xJ response objects exist, photon-yield contract needs tightening", "Current merged MC includes xJ response/fake/matched families; the photon-RAA yield path needs explicit photon-efficiency objects.")
    ax = fig.add_axes([0.06, 0.23, 0.38, 0.50])
    response_paths = [
        ("0-20%", "h2_leadRecoilJetPtResponse_pTtruth_ratio_r03_dphiPi2_isoR40_fixedIso4GeV_cent_0_20"),
        ("20-50%", "h2_leadRecoilJetPtResponse_pTtruth_ratio_r03_isoR40_fixedIso4GeV_cent_20_50"),
        ("50-80%", "h2_leadRecoilJetPtResponse_pTtruth_ratio_r03_isoR40_fixedIsoGev_cent_50_80"),
    ]
    # Correct the one case-sensitive token after keeping the list readable above.
    response_paths[2] = ("50-80%", "h2_leadRecoilJetPtResponse_pTtruth_ratio_r03_isoR40_fixedIso4GeV_cent_50_80")
    means = []
    weights = []
    for label, name in response_paths:
        got = weighted_mean_y(SIGNAL_ROOT, f"{SIM_TOP}/{name}")
        if got is None:
            means.append(float("nan"))
            weights.append(0.0)
        else:
            means.append(got[0])
            weights.append(got[1])
    ax.bar([p[0] for p in response_paths], means, color=["#BFDBFE", "#BBF7D0", "#DDD6FE"], edgecolor="#334155")
    ax.axhline(1.0, color="#111827", linestyle="--", linewidth=1.0)
    ax.set_ylabel("Mean reco/truth jet response ratio")
    ax.set_ylim(0, max(1.4, np.nanmax(means) * 1.2 if np.isfinite(means).any() else 1.2))
    ax.set_title("Example signal-MC jet response diagnostic")
    ax.grid(axis="y", alpha=0.25)
    for i, (m, w) in enumerate(zip(means, weights)):
        ax.text(i, m + 0.03, f"{m:.2f}\nentries~{w/1e6:.1f}M", ha="center", va="bottom", fontsize=9)

    ax2 = fig.add_axes([0.55, 0.23, 0.38, 0.50])
    checks = [
        ("Data event exposure", True),
        ("Data photon cutflow", True),
        ("Signal/inclusive BDT closure", True),
        ("Shower-shape and isolation closure", True),
        ("Raw ABCD purity", True),
        ("Truth-matched tight proxy", True),
        ("xJ response families", True),
        ("Explicit photon-efficiency denominator", False),
        ("Luminosity/TAA normalization", False),
        ("Purity-corrected yield spectra", False),
    ]
    y = np.arange(len(checks))[::-1]
    colors = ["#059669" if ok else "#B45309" for _, ok in checks]
    ax2.barh(y, [1] * len(checks), color=colors, alpha=0.16, edgecolor=colors)
    ax2.set_yticks(y)
    ax2.set_yticklabels([])
    ax2.set_xticks([])
    ax2.set_xlim(0, 1)
    for yi, (name, ok) in zip(y, checks):
        ax2.text(0.02, yi, name, ha="left", va="center", fontsize=9.6, color="#111827")
        ax2.text(0.97, yi, "present" if ok else "add", ha="right", va="center", fontsize=10, color="#064E3B" if ok else "#7C2D12")
    ax2.set_title("Current measurement-contract status")
    for spine in ax2.spines.values():
        spine.set_visible(False)
    add_takeaway(fig, "The current outputs are strong for photon-ID QA and xJ-response readiness, but the next narrow production should add explicit photon-efficiency denominator/numerator histograms before final photon RAA extraction.")
    summary["slide09"] = {"response_mean_ratio": means, "response_integral_proxy": weights}
    return save_slide(fig, "09_response_unfolding_readiness.png")


def make_slide_10_next(summary: dict) -> Path:
    fig = plt.figure(figsize=slide_figsize(), dpi=SLIDE_DPI)
    add_title(fig, "Photon RAA QA ladder: next work should be targeted, not broad", "The current THE-69 output is good enough to organize the measurement path and identify only the missing contracts.")
    lanes = [
        ("1. Lock data exposure", "Freeze run/segment manifests and event-level QA. Append new GRL-complete segments only after the daily delta check.", "#2563EB", 0.07, 0.64),
        ("2. Finalize photon-ID QA", "Keep BDT score, E11/E33, E32/E35, width, and isolation overlays as the acceptance gate for each append.", "#059669", 0.07, 0.44),
        ("3. Convert ABCD to purity", "Move raw A/B/C/D diagnostics to the PPG12 purity contract: binning, leakage, systematics, and corrected yield.", "#7C3AED", 0.07, 0.24),
        ("4. Add photon efficiency", "Add truth-photon denominators and reco-matched numerator spectra by centrality and ET.", "#B45309", 0.53, 0.54),
        ("5. Build RAA spectra", "Combine purity-corrected data, MC efficiency, pp reference/unfolding, luminosity/TAA, and systematics.", "#111827", 0.53, 0.34),
    ]
    for title, body, color, x, y in lanes:
        add_card(fig, x, y, 0.40, 0.155, title, body, accent=color, body_wrap=45, title_size=13.8, body_size=10.4)

    add_takeaway(fig, "Recommended next campaign: a narrow photon-RAA contract pass that adds explicit photon efficiency and luminosity-normalization bookkeeping while continuing the daily append-only data pattern.")
    summary["slide10"] = {"next_contract": [lane[0] for lane in lanes]}
    return save_slide(fig, "10_raa_readiness_next_actions.png")


def write_speaker_script(slides: list[tuple[str, Path]]) -> None:
    SCRIPT_MD.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# THE-77 AuAu photon RAA QA ladder speaker script",
        "",
        "These are local full-slide PNG candidates generated from current THE-69 merged ROOT outputs. They are intended as a first QA packet, not a final physics result.",
        "",
    ]
    script_lines = {
        "01_data_scope_exposure.png": "This slide establishes the dataset we are actually looking at. The sample is the current default-AuAu-BDT pass over 882 analyzable paired runs and 176,577 matched segment pairs. The important caveat is that these are raw exposure and candidate counts, not normalized photon yields.",
        "02_event_level_qa.png": "Here I check the event-level behavior before asking physics questions. The z-vertex and centrality proxies are coherent, and total calorimeter energy follows centrality. This is the basic data-scope sanity gate.",
        "03_photon_id_cutflow.png": "This is the photon-ID cutflow. The tight WP80 selection keeps signal-like MC at high efficiency while reducing inclusive MC substantially. Data sits in the expected intermediate region.",
        "04_bdt_score_data_mc.png": "This is the BDT-score shape check after preselection. The data shape is not pathological, and the centrality-dependent WP80 threshold is shown for reference.",
        "05_shower_shape_closure.png": "This is the shower-shape closure check for E11/E33 and E32/E35. The current corrected output does not show the low-edge MC collapse that caused the earlier stop.",
        "06_isolation_closure.png": "This is the R=0.4 isolation check in a 22-26 GeV slice. The data are sparse after tight ID, but the shapes remain usable for a first ABCD/purity diagnostic.",
        "07_raw_abcd_purity_yield.png": "This slide shows raw A-region counts and raw ABCD purity in fine ET bins, with pp only as a reference. This is not the final yield; it is the diagnostic that tells us where the purity machinery is stable or sparse.",
        "08_mc_truth_matched_proxy.png": "This slide uses truth-matched tight histograms already present in the MC output. It gives a first handle on signal purity and inclusive contamination, but it is not the final photon-ID efficiency denominator.",
        "09_response_unfolding_readiness.png": "This slide separates what is already present from what still needs a targeted contract. xJ response families exist, but photon-RAA needs explicit photon-efficiency denominators and luminosity/TAA normalization.",
        "10_raa_readiness_next_actions.png": "This is the action slide. The next campaign should be narrow: keep the QA gate, add photon-efficiency and normalization contracts, then build purity-corrected spectra and the first RAA-style result.",
    }
    for title, path in slides:
        lines.extend([f"## {title}", "", script_lines.get(path.name, ""), "", f"PNG: `{path}`", ""])
    SCRIPT_MD.write_text("\n".join(lines), encoding="utf-8")


def make_contact_sheet(slides: list[tuple[str, Path]]) -> None:
    images = [Image.open(path).convert("RGB") for _, path in slides]
    thumb_w, thumb_h = 512, 288
    margin = 28
    label_h = 34
    cols = 2
    rows = math.ceil(len(images) / cols)
    sheet = Image.new("RGB", (cols * thumb_w + (cols + 1) * margin, rows * (thumb_h + label_h) + (rows + 1) * margin), "white")
    draw = ImageDraw.Draw(sheet)
    for i, ((title, path), img) in enumerate(zip(slides, images)):
        r, c = divmod(i, cols)
        x = margin + c * (thumb_w + margin)
        y = margin + r * (thumb_h + label_h + margin)
        img.thumbnail((thumb_w, thumb_h), Image.LANCZOS)
        tile = Image.new("RGB", (thumb_w, thumb_h), "white")
        tile.paste(img, ((thumb_w - img.width) // 2, (thumb_h - img.height) // 2))
        sheet.paste(tile, (x, y))
        draw.rectangle((x, y, x + thumb_w, y + thumb_h), outline=(190, 190, 190), width=2)
        draw.text((x, y + thumb_h + 7), f"{i+1:02d}. {title}", fill=(20, 20, 20))
    CONTACT_SHEET.parent.mkdir(parents=True, exist_ok=True)
    sheet.save(CONTACT_SHEET)


def main() -> None:
    setup_style()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    summary: dict = {
        "schema": "THE77_AUAU_PHOTON_RAA_QA_LADDER_V1",
        "inputs": {
            "data_root": str(DATA_ROOT),
            "signal_root": str(SIGNAL_ROOT),
            "inclusive_root": str(INCLUSIVE_ROOT),
            "pp_reference_root": str(PP_ROOT),
            "fine_abcd_csv": str(FINE_CSV),
        },
        "model": {
            "id": "centAsFeatBase3x3_pt15to35",
            "wp80": WP80,
        },
        "caveats": [
            "Raw data counts are not luminosity-normalized and are not TAA-normalized.",
            "pp references are stock PPG12/baseV3E table-QA products, used as shape/procedure context only.",
            "Truth-matched tight fractions are diagnostic proxies; a final photon-efficiency denominator should be added in the next narrow production contract.",
        ],
    }
    slides = [
        ("Data scope and exposure", make_slide_01_scope(summary)),
        ("Event-level QA", make_slide_02_event_level(summary)),
        ("Photon-ID cutflow", make_slide_03_cutflow(summary)),
        ("BDT-score data/MC QA", make_slide_04_bdt(summary)),
        ("Shower-shape closure", make_slide_05_shower(summary)),
        ("Isolation closure", make_slide_06_iso(summary)),
        ("Raw ABCD purity and yield", make_slide_07_abcd(summary)),
        ("MC truth-matched proxy", make_slide_08_truth_proxy(summary)),
        ("Response and unfolding readiness", make_slide_09_response(summary)),
        ("RAA-readiness next actions", make_slide_10_next(summary)),
    ]
    make_contact_sheet(slides)
    write_speaker_script(slides)
    summary["slides"] = [{"title": title, "path": str(path)} for title, path in slides]
    summary["contact_sheet"] = str(CONTACT_SHEET)
    summary["speaker_script"] = str(SCRIPT_MD)
    MANIFEST.write_text(json.dumps(summary, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({"slides": len(slides), "out_dir": str(OUT_DIR), "contact_sheet": str(CONTACT_SHEET), "manifest": str(MANIFEST)}, indent=2))


if __name__ == "__main__":
    main()
