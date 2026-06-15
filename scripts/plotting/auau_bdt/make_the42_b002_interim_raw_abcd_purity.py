#!/usr/bin/env python3
"""Plot interim raw ABCD purity from current THE-42 AuAu data cache."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-codex")

import matplotlib.pyplot as plt
import ROOT


ROOT.gROOT.SetBatch(True)

REPO = Path(__file__).resolve().parents[3]
DEFAULT_AUAU_CACHE = (
    REPO
    / "dataOutput/ppg12TableQA/THE42_ppg12_tableqa_v1_auauDataNewV008_b002_20260612"
    / "interim_complete_runs"
    / "the42_b002_interim_complete_run_raw_abcd_counts.json"
)
DEFAULT_PP_ROOT = (
    REPO
    / "dataOutput/ppg12TableQA/THE42_ppg12_tableqa_v1_basev3e_20260611"
    / "merged_roots"
    / "RecoilJets_pp_ALL_preselectionNewPPG12_tightNewPPG12_nonTightNewPPG12.root"
)
DEFAULT_OUTDIR = DEFAULT_AUAU_CACHE.parent / "raw_abcd_purity"

REGION_PREFIX = {
    "A": "h_isIsolated_isTight",
    "B": "h_notIsolated_isTight",
    "C": "h_isIsolated_notTight",
    "D": "h_notIsolated_notTight",
}
CENT_LABELS = {
    "0_20": "AuAu 0-20%",
    "20_50": "AuAu 20-50%",
    "50_80": "AuAu 50-80%",
}
CENT_COLORS = {
    "0_20": "#111111",
    "20_50": "#1f77b4",
    "50_80": "#d95f02",
}


def hist_sum_and_error(hist) -> tuple[float, float]:
    total = 0.0
    err2 = 0.0
    for ibin in range(0, hist.GetNbinsX() + 2):
        total += float(hist.GetBinContent(ibin))
        err2 += float(hist.GetBinError(ibin)) ** 2
    if err2 <= 0.0 and total > 0.0:
        err2 = total
    return total, math.sqrt(err2)


def raw_abcd_purity(counts: dict[str, tuple[float, float]]) -> tuple[float, float, float]:
    a, ea = counts["A"]
    b, eb = counts["B"]
    c, ec = counts["C"]
    d, ed = counts["D"]
    if a <= 0.0 or d <= 0.0:
        return math.nan, math.nan, math.nan
    bkg_est = b * c / d
    raw = 1.0 - bkg_est / a
    value = max(0.0, raw)
    dp_da = (b * c) / (a * a * d)
    dp_db = -c / (a * d)
    dp_dc = -b / (a * d)
    dp_dd = (b * c) / (a * d * d)
    var = (
        dp_da * dp_da * ea * ea
        + dp_db * dp_db * eb * eb
        + dp_dc * dp_dc * ec * ec
        + dp_dd * dp_dd * ed * ed
    )
    return value, math.sqrt(max(var, 0.0)), bkg_est


def load_auau_rows(path: Path) -> tuple[list[dict], dict]:
    payload = json.loads(path.read_text())
    allowed_schemas = {
        "THE42_B002_INTERIM_RAW_ABCD_COUNTS_V1",
        "THE42_CURRENT_AUAU_INTERIM_ELIGIBLE_RAW_ABCD_COUNTS_V1",
    }
    if payload.get("schema") not in allowed_schemas:
        raise RuntimeError(f"Unexpected cache schema in {path}")
    rows = []
    for row in payload["rows"]:
        counts = {r: (float(row[r]), float(row[f"{r}_error"])) for r in "ABCD"}
        value, error, bkg_est = raw_abcd_purity(counts)
        rows.append(
            {
                "source": "current THE-42 AuAu interim data",
                "centrality": row["centrality"],
                "pt_low": float(row["pt_low"]),
                "pt_high": float(row["pt_high"]),
                "pt_mid": 0.5 * (float(row["pt_low"]) + float(row["pt_high"])),
                "value": value,
                "error": error,
                "abcd_background_estimate": bkg_est,
                **{k: row[k] for k in row if k in {"A", "B", "C", "D", "A_error", "B_error", "C_error", "D_error"}},
            }
        )
    return rows, payload.get("metadata", {})


def read_pp_rows(path: Path, topdir: str, cone: str, pt_min: float, pt_max: float) -> list[dict]:
    handle = ROOT.TFile.Open(str(path), "READ")
    if not handle or handle.IsZombie():
        raise RuntimeError(f"Could not open {path}")
    try:
        directory = handle.Get(topdir)
        if not directory:
            raise RuntimeError(f"Missing topdir {topdir} in {path}")
        pattern = re.compile(
            rf"^({'|'.join(re.escape(v) for v in REGION_PREFIX.values())})_{re.escape(cone)}_pT_([0-9]+)_([0-9]+)$"
        )
        counts: dict[tuple[float, float], dict[str, tuple[float, float]]] = {}
        for key in directory.GetListOfKeys():
            name = key.GetName()
            match = pattern.match(name)
            if not match:
                continue
            prefix, lo_s, hi_s = match.groups()
            lo, hi = float(lo_s), float(hi_s)
            mid = 0.5 * (lo + hi)
            if mid < pt_min or mid > pt_max:
                continue
            region = next(r for r, p in REGION_PREFIX.items() if p == prefix)
            hist = key.ReadObj()
            counts.setdefault((lo, hi), {})[region] = hist_sum_and_error(hist)

        rows = []
        for (lo, hi), region_counts in sorted(counts.items()):
            if any(r not in region_counts for r in "ABCD"):
                continue
            value, error, bkg_est = raw_abcd_purity(region_counts)
            out = {
                "source": "pp data checkpoint",
                "centrality": "pp",
                "pt_low": lo,
                "pt_high": hi,
                "pt_mid": 0.5 * (lo + hi),
                "value": value,
                "error": error,
                "abcd_background_estimate": bkg_est,
            }
            for region in "ABCD":
                out[region] = region_counts[region][0]
                out[f"{region}_error"] = region_counts[region][1]
            rows.append(out)
        return rows
    finally:
        handle.Close()


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "source",
        "centrality",
        "pt_low",
        "pt_high",
        "pt_mid",
        "value",
        "error",
        "abcd_background_estimate",
        "A",
        "A_error",
        "B",
        "B_error",
        "C",
        "C_error",
        "D",
        "D_error",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def plot(path: Path, auau_rows: list[dict], pp_rows: list[dict], meta: dict, min_a: float) -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "axes.linewidth": 1.05,
            "xtick.direction": "in",
            "ytick.direction": "in",
        }
    )
    fig, axes = plt.subplots(1, 2, figsize=(15.8, 6.2), sharey=True)
    fig.patch.set_facecolor("white")

    def draw_points(ax, rows, label, color, marker="o", open_marker=False, require_min_a=False):
        good = [r for r in rows if math.isfinite(float(r["value"]))]
        if require_min_a:
            good = [r for r in good if float(r.get("A", 0.0)) >= min_a]
        if not good:
            return
        x = [float(r["pt_mid"]) for r in good]
        y = [float(r["value"]) for r in good]
        yerr = [float(r["error"]) for r in good]
        xerr = [0.5 * (float(r["pt_high"]) - float(r["pt_low"])) for r in good]
        ax.errorbar(
            x,
            y,
            xerr=xerr,
            yerr=yerr,
            linestyle="none",
            marker=marker,
            markersize=7.0,
            markerfacecolor="white" if open_marker else color,
            markeredgecolor=color,
            markeredgewidth=1.7 if open_marker else 0.8,
            color=color,
            ecolor=color,
            elinewidth=1.4,
            capsize=3,
            label=label,
        )

    ax = axes[0]
    for cent in ("0_20", "20_50", "50_80"):
        draw_points(
            ax,
            [r for r in auau_rows if r["centrality"] == cent],
            CENT_LABELS[cent],
            CENT_COLORS[cent],
            require_min_a=True,
        )
    ax.set_title("AuAu centrality split", fontsize=16, fontweight="bold")
    ax.set_ylabel(r"raw ABCD purity  $\max(0, A-BC/D)/A$", fontsize=13)
    ax.legend(frameon=False, fontsize=11, loc="lower right")

    ax = axes[1]
    draw_points(
        ax,
        [r for r in auau_rows if r["centrality"] == "50_80"],
        "AuAu 50-80%",
        CENT_COLORS["50_80"],
        require_min_a=True,
    )
    draw_points(ax, pp_rows, "pp checkpoint data", "#0072B2", marker="s", open_marker=True)
    ax.set_title("Peripheral AuAu compared with pp", fontsize=16, fontweight="bold")
    ax.legend(frameon=False, fontsize=11, loc="lower right")

    for ax in axes:
        ax.set_xlim(14.0, 36.0)
        ax.set_ylim(-0.05, 1.05)
        ax.grid(True, axis="y", color="#dfe5ee", linewidth=0.9)
        ax.tick_params(axis="both", which="major", labelsize=11, top=True, right=True, length=6)
        ax.minorticks_on()
        ax.tick_params(axis="both", which="minor", top=True, right=True, length=3)
        ax.set_xlabel(r"photon candidate $E_T$ [GeV]", fontsize=13)

    fig.text(0.045, 0.965, r"$\bf{\it{sPHENIX}}$ Internal", ha="left", va="top", fontsize=16)
    fig.text(
        0.045,
        0.912,
        "Interim current-campaign data-only ABCD check",
        ha="left",
        va="top",
        fontsize=13,
        color="#333333",
    )
    detail = (
        f"Current THE-42 AuAu: {meta.get('files_opened', meta.get('complete_root_files_used', '?'))} "
        f"merged per-run ROOTs; plotted AuAu bins require A >= {min_a:g}; "
        r"$15<E_T<35$ GeV, $\Delta R=0.3$ sliding isolation"
    )
    fig.text(0.045, 0.868, detail, ha="left", va="top", fontsize=12.2, color="#555555")
    fig.tight_layout(rect=(0.035, 0.055, 0.995, 0.835), w_pad=2.0)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=220)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--auau-cache", type=Path, default=DEFAULT_AUAU_CACHE)
    parser.add_argument("--pp-root", type=Path, default=DEFAULT_PP_ROOT)
    parser.add_argument("--pp-topdir", default="Photon_4_GeV_plus_MBD_NS_geq_1")
    parser.add_argument("--cone", default="isoR30_isSliding")
    parser.add_argument("--pt-min", type=float, default=15.0)
    parser.add_argument("--pt-max", type=float, default=35.0)
    parser.add_argument("--min-a", type=float, default=5.0)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = parser.parse_args()

    auau_rows, meta = load_auau_rows(args.auau_cache)
    pp_rows = read_pp_rows(args.pp_root, args.pp_topdir, args.cone, args.pt_min, args.pt_max)
    rows = auau_rows + pp_rows
    out_csv = args.outdir / "current_auau_interim_raw_abcd_purity_points.csv"
    out_png = args.outdir / "current_auau_interim_raw_abcd_purity_auau_cent_and_pp50_80_overlay.png"
    manifest = args.outdir / "current_auau_interim_raw_abcd_purity_manifest.json"
    write_csv(out_csv, rows)
    plot(out_png, auau_rows, pp_rows, meta, args.min_a)
    manifest.write_text(
        json.dumps(
            {
                "schema": "THE42_CURRENT_AUAU_INTERIM_RAW_ABCD_PURITY_PLOT_V1",
                "png": str(out_png),
                "csv": str(out_csv),
                "auau_cache": str(args.auau_cache),
                "pp_root": str(args.pp_root),
                "pp_topdir": args.pp_topdir,
                "cone": args.cone,
                "pt_range": [args.pt_min, args.pt_max],
                "auau_plot_min_A": args.min_a,
                "caveat": "Interim current THE-42 AuAu data-only ABCD purity from eligible merged outputs; regenerate after the full current campaign drains/merges. Not a matched-sim efficiency or truth-purity result.",
                "metadata": meta,
            },
            indent=2,
        )
        + "\n"
    )
    print(out_png)
    print(out_csv)
    print(manifest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
