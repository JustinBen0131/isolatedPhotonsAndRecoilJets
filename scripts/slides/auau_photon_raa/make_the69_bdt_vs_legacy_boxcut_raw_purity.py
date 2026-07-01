#!/usr/bin/env python3
"""Compare THE-69 AuAu BDT raw ABCD purity to legacy AuAu box-cut purity."""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


THIS_FILE = Path(__file__).resolve()
REPO = next((p for p in THIS_FILE.parents if (p / "AGENTS.md").exists()), THIS_FILE.parents[4])
SCRIPTS = REPO / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.append(str(SCRIPTS))

from slides.common.slide_defaults import SLIDE_DPI, slide_figsize  # noqa: E402


CURRENT_BDT_CSV = (
    REPO
    / "dataOutput/auauPhysicsQA/THE69_defaultAuAuBDT_physicsQA_20260620/"
    / "fine_yield_abcd_companion_slides/the69_raw_abcd_purity_points.csv"
)
LEGACY_BOX_ROOT = (
    REPO
    / "InputFiles/auau25/"
    / "RecoilJets_auau_ALL_jetMinPt5_7pi_8_vz60_isoR40_isSliding_baseVariant_"
    / "preselectionReference_tightReference_nonTightReference.root"
)
if not LEGACY_BOX_ROOT.exists():
    LEGACY_BOX_ROOT = (
        REPO
        / "InputFiles/auau25/"
        / "RecoilJets_auau_ALL_jetMinPt5_7pi_8_vz60_isoR40_isSliding_baseVariant_preselectionReference_tightReference_nonTightReference.root"
    )

OUT_DIR = (
    REPO
    / "dataOutput/auauPhysicsQA/THE69_defaultAuAuBDT_physicsQA_20260620/"
    / "bdt_vs_legacy_boxcut_raw_purity"
)
OUT_PNG = OUT_DIR / "the69_bdt_vs_legacy_boxcut_raw_abcd_purity.png"
OUT_CSV = OUT_DIR / "the69_bdt_vs_legacy_boxcut_raw_abcd_purity_points.csv"
OUT_MANIFEST = OUT_DIR / "the69_bdt_vs_legacy_boxcut_raw_abcd_purity_manifest.json"

TOP_DIR = "MBD_NS_geq_2_vtx_lt_150"
CENTRALITIES = [("0_20", "0-20%"), ("20_50", "20-50%"), ("50_80", "50-80%")]
PT_BINS = [(14, 16), (16, 18), (18, 20), (20, 22), (22, 24), (24, 26), (26, 35)]


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "axes.facecolor": "white",
            "axes.edgecolor": "#1F2937",
            "axes.linewidth": 1.05,
            "axes.labelsize": 14.5,
            "axes.titlesize": 16.0,
            "xtick.labelsize": 12.0,
            "ytick.labelsize": 12.0,
            "legend.fontsize": 11.4,
        }
    )


def fnum(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def raw_purity(a: float, b: float, c: float, d: float) -> float:
    if a <= 0.0 or d <= 0.0:
        return float("nan")
    return max(0.0, min(1.0, (a - b * c / d) / a))


def raw_purity_error(a: float, b: float, c: float, d: float) -> float:
    if a <= 0.0 or d <= 0.0:
        return float("nan")
    d_p_da = (b * c) / (a * a * d)
    d_p_db = -c / (a * d)
    d_p_dc = -b / (a * d)
    d_p_dd = (b * c) / (a * d * d)
    var = d_p_da * d_p_da * a + d_p_db * d_p_db * b + d_p_dc * d_p_dc * c + d_p_dd * d_p_dd * d
    return math.sqrt(var) if var > 0.0 else 0.0


def hist_integral_and_error(hist) -> tuple[float, float]:
    if not hist:
        return float("nan"), float("nan")
    total = 0.0
    err2 = 0.0
    for idx in range(1, hist.GetNbinsX() + 1):
        val = float(hist.GetBinContent(idx))
        err = float(hist.GetBinError(idx))
        if err <= 0.0 and val > 0.0:
            err = math.sqrt(val)
        total += val
        err2 += err * err
    return total, math.sqrt(err2)


def load_current_bdt() -> dict[tuple[str, int, int], dict]:
    out: dict[tuple[str, int, int], dict] = {}
    with CURRENT_BDT_CSV.open() as f:
        for row in csv.DictReader(f):
            key = (row["cent_key"], int(row["pt_lo"]), int(row["pt_hi"]))
            out[key] = {
                "purity": fnum(row["raw_purity"]),
                "purity_err": fnum(row["raw_purity_err"]),
                "A": fnum(row["A"]),
                "B": fnum(row["B"]),
                "C": fnum(row["C"]),
                "D": fnum(row["D"]),
                "source": row["A_source"],
            }
    return out


def load_legacy_box() -> dict[tuple[str, int, int], dict]:
    import ROOT

    ROOT.gROOT.SetBatch(True)
    ROOT.TH1.AddDirectory(False)
    root_file = ROOT.TFile.Open(str(LEGACY_BOX_ROOT), "READ")
    if not root_file or root_file.IsZombie():
        raise OSError(f"Could not open legacy box-cut ROOT: {LEGACY_BOX_ROOT}")
    directory = root_file.Get(TOP_DIR)
    if not directory:
        raise KeyError(f"Missing {TOP_DIR} in {LEGACY_BOX_ROOT}")

    out: dict[tuple[str, int, int], dict] = {}
    try:
        for cent_key, _ in CENTRALITIES:
            for lo, hi in PT_BINS:
                vals: dict[str, float] = {}
                errs: dict[str, float] = {}
                sources: dict[str, str] = {}
                for region in ("A", "B", "C", "D"):
                    hist_name = f"h_Eiso_ABCD_{region}_isoR40_isSliding_pT_{lo}_{hi}_cent_{cent_key}"
                    hist = directory.Get(hist_name)
                    val, err = hist_integral_and_error(hist)
                    vals[region] = val
                    errs[region] = err
                    sources[region] = f"{TOP_DIR}/{hist_name}"
                a, b, c, d = vals["A"], vals["B"], vals["C"], vals["D"]
                out[(cent_key, lo, hi)] = {
                    "purity": raw_purity(a, b, c, d),
                    "purity_err": raw_purity_error(a, b, c, d),
                    "A": a,
                    "B": b,
                    "C": c,
                    "D": d,
                    "source": sources["A"],
                }
    finally:
        root_file.Close()
    return out


def build_points() -> list[dict]:
    current = load_current_bdt()
    legacy = load_legacy_box()
    points: list[dict] = []
    for cent_key, cent_label in CENTRALITIES:
        for lo, hi in PT_BINS:
            width = hi - lo
            mid = 0.5 * (lo + hi)
            for method, label, data in (
                ("current_bdt", "THE-69 BDT WP80", current.get((cent_key, lo, hi), {})),
                ("legacy_box", "Legacy box cuts", legacy.get((cent_key, lo, hi), {})),
            ):
                points.append(
                    {
                        "method": method,
                        "label": label,
                        "centrality": cent_label,
                        "cent_key": cent_key,
                        "pt_lo": lo,
                        "pt_hi": hi,
                        "pt_mid": mid,
                        "pt_width": width,
                        "raw_purity": data.get("purity", float("nan")),
                        "raw_purity_err": data.get("purity_err", float("nan")),
                        "A": data.get("A", float("nan")),
                        "B": data.get("B", float("nan")),
                        "C": data.get("C", float("nan")),
                        "D": data.get("D", float("nan")),
                        "source": data.get("source", ""),
                    }
                )
    return points


def write_points(points: list[dict]) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with OUT_CSV.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(points[0].keys()))
        writer.writeheader()
        writer.writerows(points)


def draw(points: list[dict]) -> None:
    setup_style()
    fig = plt.figure(figsize=slide_figsize(), dpi=SLIDE_DPI)
    gs = fig.add_gridspec(
        2,
        3,
        height_ratios=[0.22, 1.0],
        left=0.065,
        right=0.975,
        top=0.91,
        bottom=0.105,
        hspace=0.18,
        wspace=0.18,
    )
    title_ax = fig.add_subplot(gs[0, :])
    title_ax.axis("off")
    title_ax.text(
        0.0,
        0.92,
        "Raw ABCD purity: current Au+Au BDT vs legacy box cuts",
        ha="left",
        va="top",
        fontsize=26.5,
        fontweight="bold",
        color="#111827",
    )
    title_ax.text(
        0.0,
        0.26,
        r"Historical overlay only: legacy box cuts use an older Au+Au ROOT, not the current THE-69 run list.",
        ha="left",
        va="center",
        fontsize=15.5,
        color="#374151",
    )

    axes = [fig.add_subplot(gs[1, idx]) for idx in range(3)]
    colors = {"current_bdt": "#111827", "legacy_box": "#D62728"}
    markers = {"current_bdt": "o", "legacy_box": "s"}
    labels = {"current_bdt": "Current BDT WP80", "legacy_box": "Legacy box cuts"}

    for ax, (cent_key, cent_label) in zip(axes, CENTRALITIES):
        for method in ("current_bdt", "legacy_box"):
            rows = [p for p in points if p["cent_key"] == cent_key and p["method"] == method]
            x = np.array([p["pt_mid"] for p in rows], dtype=float)
            xerr = np.array([0.5 * p["pt_width"] for p in rows], dtype=float)
            y = np.array([p["raw_purity"] for p in rows], dtype=float)
            yerr = np.array([p["raw_purity_err"] for p in rows], dtype=float)
            good = np.isfinite(y)
            ax.errorbar(
                x[good],
                y[good],
                xerr=xerr[good],
                yerr=yerr[good],
                fmt=markers[method],
                ms=8.0 if method == "current_bdt" else 8.4,
                mfc=colors[method] if method == "current_bdt" else "white",
                mec="white" if method == "current_bdt" else colors[method],
                mew=0.8 if method == "current_bdt" else 2.0,
                color=colors[method],
                ecolor=colors[method],
                elinewidth=1.1,
                capsize=3.0,
                label=labels[method],
            )
        ax.set_title(cent_label, fontweight="bold", pad=8)
        ax.set_xlim(13.4, 35.8)
        ax.set_ylim(-0.05, 1.08)
        ax.set_xticks([14, 16, 18, 20, 22, 24, 26, 30, 34])
        ax.set_yticks(np.arange(0.0, 1.01, 0.2))
        ax.grid(True, which="major", color="#D9E2EC", alpha=0.9, linewidth=0.85)
        ax.set_xlabel(r"cluster $E_T$ [GeV]")
        if ax is axes[0]:
            ax.set_ylabel("Raw ABCD purity")
            ax.legend(
                loc="lower left",
                frameon=True,
                framealpha=0.96,
                facecolor="white",
                edgecolor="#CBD5E1",
                handlelength=1.2,
                borderpad=0.45,
            )
        else:
            ax.tick_params(labelleft=False)

    fig.savefig(OUT_PNG, dpi=SLIDE_DPI, facecolor="white")
    plt.close(fig)


def write_manifest(points: list[dict]) -> None:
    manifest = {
        "schema": "THE69_BDT_VS_LEGACY_BOX_RAW_ABCD_PURITY_V1",
        "output_png": str(OUT_PNG),
        "points_csv": str(OUT_CSV),
        "current_bdt_csv": str(CURRENT_BDT_CSV),
        "legacy_box_root": str(LEGACY_BOX_ROOT),
        "top_dir": TOP_DIR,
        "centralities": CENTRALITIES,
        "pt_bins": PT_BINS,
        "raw_purity_definition": "P_raw=(A-B*C/D)/A using h_Eiso_ABCD_{A,B,C,D}_isoR40_isSliding integrals.",
        "current_bdt_definition": "THE-69 default AuAu BDT, centAsFeatBase3x3_pt15to35, WP80 T80(c)=0.53471108+0.0012284143*c.",
        "legacy_box_definition": {
            "variant": "preselectionReference_tightReference_nonTightReference",
            "preselection": {
                "e11e33_max": 0.98,
                "et1_min": 0.60,
                "et1_max": 1.00,
                "e32e35_min": 0.80,
                "e32e35_max": 1.00,
                "weta_max": 0.60,
            },
            "tight": {
                "w_lo": 0.0,
                "w_hi": "0.15 + 0.006 * cluster_Et",
                "e11e33_min": 0.40,
                "e11e33_max": 0.98,
                "et1_min": 0.90,
                "et1_max": 1.00,
                "e32e35_min": 0.92,
                "e32e35_max": 1.00,
            },
            "source": "src_AuAu/RecoilJets_AuAu.h constants used by src_AuAu/RecoilJets_AuAu.cc ABCD tight-axis logic.",
        },
        "comparison_caveat": "Current BDT points use the THE-69 current-production data merge; legacy box points use the existing old AuAu box-cut ROOT, not a same-segment current-data rerun. A fair cut-to-cut comparison requires rerunning the current THE-69 run list with a legacy box-cut photon-ID row.",
        "bins": len(points),
    }
    OUT_MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n")


def main() -> int:
    points = build_points()
    write_points(points)
    write_manifest(points)
    draw(points)
    print(OUT_PNG)
    print(OUT_CSV)
    print(OUT_MANIFEST)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
