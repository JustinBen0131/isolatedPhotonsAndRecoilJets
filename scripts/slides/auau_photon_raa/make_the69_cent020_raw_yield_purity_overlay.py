#!/usr/bin/env python3
"""Build a 0-20% raw-yield plus raw/corrected ABCD purity overlay for THE-69."""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path
from typing import Iterable

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


RAW_ABCD_CSV = (
    REPO
    / "dataOutput/auauPhysicsQA/THE69_defaultAuAuBDT_physicsQA_20260620/"
    / "fine_yield_abcd_companion_slides/the69_raw_abcd_purity_points.csv"
)
LEAKAGE_CSV = (
    REPO
    / "dataOutput/auauPhysicsQA/THE69_leakageCentWP_fix_20260620/"
    / "the69_leakage_centwp_fix_points.csv"
)
PP_EXACT_BIN_DATA_ROOT = (
    REPO
    / "dataOutput/ppg12TableQA/THE42_ppg12_tableqa_v1_basev3e_20260611/"
    / "merged_roots/RecoilJets_pp_ALL_preselectionNewPPG12_tightNewPPG12_nonTightNewPPG12.root"
)
OUT_DIR = (
    REPO
    / "dataOutput/auauPhysicsQA/THE69_defaultAuAuBDT_physicsQA_20260620/"
    / "cent020_raw_yield_purity_overlay"
)
OUT_PNG = OUT_DIR / "the69_cent020_raw_yield_raw_vs_leakage_corrected_purity.png"
OUT_CSV = OUT_DIR / "the69_cent020_raw_yield_raw_vs_leakage_corrected_purity_points.csv"
OUT_MANIFEST = OUT_DIR / "the69_cent020_raw_yield_raw_vs_leakage_corrected_purity_manifest.json"

CENT_KEY = "0_20"
CENT_LABEL = "0-20%"
PP_TRIGGER_DIR = "Photon_4_GeV_plus_MBD_NS_geq_1"
VALID_PT_LO = 16
VALID_PT_HI = 35


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
            "axes.labelsize": 15,
            "axes.titlesize": 17,
            "xtick.labelsize": 12.5,
            "ytick.labelsize": 12.5,
            "legend.fontsize": 12.0,
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


def raw_signal_yield(a: float, b: float, c: float, d: float) -> float:
    if a <= 0.0 or d <= 0.0:
        return float("nan")
    return max(a - b * c / d, 0.0)


def raw_signal_yield_error(a: float, b: float, c: float, d: float) -> float:
    if not all(math.isfinite(v) for v in (a, b, c, d)) or d <= 0.0:
        return float("nan")
    q = b * c / d
    var_q = 0.0
    if b > 0.0:
        var_q += (q / b) ** 2 * b
    if c > 0.0:
        var_q += (q / c) ** 2 * c
    if d > 0.0:
        var_q += (q / d) ** 2 * d
    return math.sqrt(max(a + var_q, 0.0))


def solve_leakage_corrected_sa(
    a: float, b: float, c: float, d: float, f_b: float, f_c: float, f_d: float
) -> tuple[float, bool]:
    """Python copy of the RecoilJets fixed-point leakage correction."""
    if a <= 0.0:
        return 0.0, True

    s = a
    if d != 0.0:
        s = min(max(a - b * (c / d), 0.0), a)

    def fixed_point(value: float) -> float:
        denom = d - f_d * value
        if denom == 0.0:
            return float("nan")
        return a - (b - f_b * value) * (c - f_c * value) / denom

    lam = 0.25
    for iteration in range(200):
        if f_d > 0.0:
            s_max = (d / f_d) * 0.999
            if math.isfinite(s_max):
                s = min(s, max(0.0, s_max))

        f_val = fixed_point(s)
        if not math.isfinite(f_val):
            return s, False

        s_new = (1.0 - lam) * s + lam * f_val
        if not math.isfinite(s_new):
            return s, False

        s_new = min(max(s_new, 0.0), a)
        delta = abs(s_new - s)
        s = s_new
        if delta < 1.0e-6:
            return s, True
        if iteration > 10 and delta > 0.5 * a and lam > 0.05:
            lam *= 0.5

    return s, False


def corrected_purity_value(
    a: float, b: float, c: float, d: float, f_b: float, f_c: float, f_d: float
) -> tuple[float, bool]:
    sa, ok = solve_leakage_corrected_sa(a, b, c, d, f_b, f_c, f_d)
    if ok and a > 0.0:
        return max(0.0, min(1.0, sa / a)), True
    return raw_purity(a, b, c, d), False


def corrected_purity_error(
    a: float,
    b: float,
    c: float,
    d: float,
    f_b: float,
    f_c: float,
    f_d: float,
    f_b_err: float,
    f_c_err: float,
    f_d_err: float,
    has_correction: bool,
) -> float:
    if not has_correction:
        return raw_purity_error(a, b, c, d)
    if a <= 0.0:
        return float("nan")

    base = [a, b, c, d, f_b, f_c, f_d]
    widths = [math.sqrt(max(v, 1.0)) for v in (a, b, c, d)] + [
        max(f_b_err, 0.0),
        max(f_c_err, 0.0),
        max(f_d_err, 0.0),
    ]

    def value(vals: Iterable[float]) -> float:
        aa, bb, cc, dd, ff_b, ff_c, ff_d = vals
        pur, _ = corrected_purity_value(aa, bb, cc, dd, ff_b, ff_c, ff_d)
        return pur

    var = 0.0
    for idx, sigma in enumerate(widths):
        if sigma <= 0.0:
            continue
        up = list(base)
        down = list(base)
        up[idx] = max(0.0, up[idx] + sigma)
        down[idx] = max(0.0, down[idx] - sigma)
        denom = up[idx] - down[idx]
        if denom <= 0.0:
            continue
        derivative = (value(up) - value(down)) / denom
        var += derivative * derivative * sigma * sigma
    return math.sqrt(var) if var > 0.0 else 0.0


def load_raw_rows() -> list[dict]:
    rows: list[dict] = []
    with RAW_ABCD_CSV.open() as f:
        for row in csv.DictReader(f):
            if row["cent_key"] != CENT_KEY:
                continue
            lo = int(row["pt_lo"])
            hi = int(row["pt_hi"])
            if lo < VALID_PT_LO or hi > VALID_PT_HI:
                continue
            rows.append(row)
    rows.sort(key=lambda r: int(r["pt_lo"]))
    return rows


def load_leakage_map() -> dict[tuple[int, int], dict[str, tuple[float, float]]]:
    leakage: dict[tuple[int, int], dict[str, tuple[float, float]]] = {}
    with LEAKAGE_CSV.open() as f:
        for row in csv.DictReader(f):
            if row["source"] != "corrected_auau" or row["cent_key"] != CENT_KEY:
                continue
            lo = int(row["pt_lo"])
            hi = int(row["pt_hi"])
            if lo < VALID_PT_LO or hi > VALID_PT_HI:
                continue
            leakage.setdefault((lo, hi), {})[row["region"]] = (
                fnum(row["leakage_fraction"]),
                fnum(row["leakage_fraction_err"]),
            )
    return leakage


def load_pp_exact_bin_abcd() -> dict[tuple[int, int], dict[str, tuple[float, float, str]]]:
    import ROOT

    ROOT.gROOT.SetBatch(True)
    ROOT.TH1.AddDirectory(False)
    root_file = ROOT.TFile.Open(str(PP_EXACT_BIN_DATA_ROOT), "READ")
    if not root_file or root_file.IsZombie():
        raise OSError(f"Could not open pp exact-bin ROOT: {PP_EXACT_BIN_DATA_ROOT}")
    directory = root_file.Get(PP_TRIGGER_DIR)
    if not directory:
        raise KeyError(f"Missing pp trigger directory {PP_TRIGGER_DIR} in {PP_EXACT_BIN_DATA_ROOT}")
    out: dict[tuple[int, int], dict[str, tuple[float, float, str]]] = {}
    try:
        for lo, hi in [(16, 18), (18, 20), (20, 22), (22, 24), (24, 26), (26, 35)]:
            out[(lo, hi)] = {}
            for region in ("A", "B", "C", "D"):
                hist_name = f"h_Eiso_ABCD_{region}_isoR40_isSliding_pT_{lo}_{hi}"
                hist_path = f"{PP_TRIGGER_DIR}/{hist_name}"
                hist = directory.Get(hist_name)
                if not hist:
                    raise KeyError(f"Missing pp ABCD histogram {hist_path}")
                value = 0.0
                err2 = 0.0
                for bin_idx in range(1, hist.GetNbinsX() + 1):
                    value += float(hist.GetBinContent(bin_idx))
                    err = float(hist.GetBinError(bin_idx))
                    if err <= 0.0 and hist.GetBinContent(bin_idx) > 0.0:
                        err = math.sqrt(float(hist.GetBinContent(bin_idx)))
                    err2 += err * err
                out[(lo, hi)][region] = (value, math.sqrt(err2), hist_path)
    finally:
        root_file.Close()
    return out


def build_points() -> list[dict]:
    raw_rows = load_raw_rows()
    leakage = load_leakage_map()
    pp_abcd = load_pp_exact_bin_abcd()
    points: list[dict] = []
    for row in raw_rows:
        lo = int(row["pt_lo"])
        hi = int(row["pt_hi"])
        factors = leakage.get((lo, hi), {})
        if not all(region in factors for region in ("B", "C", "D")):
            continue

        a = fnum(row["A"])
        b = fnum(row["B"])
        c = fnum(row["C"])
        d = fnum(row["D"])
        f_b, f_b_err = factors["B"]
        f_c, f_c_err = factors["C"]
        f_d, f_d_err = factors["D"]
        width = fnum(row["width"])
        auau_signal = raw_signal_yield(a, b, c, d)
        auau_signal_err = raw_signal_yield_error(a, b, c, d)
        pp_regions = pp_abcd[(lo, hi)]
        pp_a, pp_a_err, pp_a_source = pp_regions["A"]
        pp_b, pp_b_err, pp_b_source = pp_regions["B"]
        pp_c, pp_c_err, pp_c_source = pp_regions["C"]
        pp_d, pp_d_err, pp_d_source = pp_regions["D"]
        pp_signal = raw_signal_yield(pp_a, pp_b, pp_c, pp_d)
        pp_signal_err = raw_signal_yield_error(pp_a, pp_b, pp_c, pp_d)
        corrected, ok = corrected_purity_value(a, b, c, d, f_b, f_c, f_d)
        points.append(
            {
                "centrality": row["centrality"],
                "cent_key": row["cent_key"],
                "pt_lo": lo,
                "pt_hi": hi,
                "pt_mid": fnum(row["pt_mid"]),
                "width": width,
                "A_data": a,
                "A_data_err": fnum(row["A_err"]),
                "A_data_per_gev": fnum(row["A_per_gev"]),
                "A_data_per_gev_err": fnum(row["A_err"]) / width,
                "B_data": b,
                "C_data": c,
                "D_data": d,
                "raw_signal_yield": auau_signal,
                "raw_signal_yield_err": auau_signal_err,
                "raw_signal_yield_per_gev": auau_signal / width,
                "raw_signal_yield_per_gev_err": auau_signal_err / width,
                "pp_A_data": pp_a,
                "pp_A_data_err": pp_a_err,
                "pp_A_data_per_gev": pp_a / width,
                "pp_A_data_per_gev_err": pp_a_err / width,
                "pp_B_data": pp_b,
                "pp_B_data_err": pp_b_err,
                "pp_C_data": pp_c,
                "pp_C_data_err": pp_c_err,
                "pp_D_data": pp_d,
                "pp_D_data_err": pp_d_err,
                "pp_raw_signal_yield": pp_signal,
                "pp_raw_signal_yield_err": pp_signal_err,
                "pp_raw_signal_yield_per_gev": pp_signal / width,
                "pp_raw_signal_yield_per_gev_err": pp_signal_err / width,
                "pp_A_data_source": pp_a_source,
                "pp_B_data_source": pp_b_source,
                "pp_C_data_source": pp_c_source,
                "pp_D_data_source": pp_d_source,
                "raw_purity": raw_purity(a, b, c, d),
                "raw_purity_err": raw_purity_error(a, b, c, d),
                "corrected_purity": corrected,
                "corrected_purity_err": corrected_purity_error(
                    a, b, c, d, f_b, f_c, f_d, f_b_err, f_c_err, f_d_err, ok
                ),
                "correction_converged": ok,
                "fB": f_b,
                "fB_err": f_b_err,
                "fC": f_c,
                "fC_err": f_c_err,
                "fD": f_d,
                "fD_err": f_d_err,
            }
        )
    return points


def write_points(points: list[dict]) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fields = list(points[0].keys()) if points else []
    with OUT_CSV.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(points)


def draw(points: list[dict]) -> None:
    setup_style()
    x = np.array([p["pt_mid"] for p in points], dtype=float)
    xerr = np.array([0.5 * p["width"] for p in points], dtype=float)
    raw_y = np.array([p["A_data_per_gev"] for p in points], dtype=float)
    raw_yerr = np.array([p["A_data_per_gev_err"] for p in points], dtype=float)
    pp_y = np.array([p["pp_A_data_per_gev"] for p in points], dtype=float)
    pp_yerr = np.array([p["pp_A_data_per_gev_err"] for p in points], dtype=float)
    raw_p = np.array([p["raw_purity"] for p in points], dtype=float)
    raw_perr = np.array([p["raw_purity_err"] for p in points], dtype=float)
    corr_p = np.array([p["corrected_purity"] for p in points], dtype=float)
    corr_perr = np.array([p["corrected_purity_err"] for p in points], dtype=float)

    fig = plt.figure(figsize=slide_figsize(), dpi=SLIDE_DPI)
    gs = fig.add_gridspec(
        2,
        1,
        left=0.095,
        right=0.965,
        top=0.83,
        bottom=0.105,
        height_ratios=[1.0, 1.05],
        hspace=0.23,
    )
    ax_yield = fig.add_subplot(gs[0])
    ax_purity = fig.add_subplot(gs[1], sharex=ax_yield)

    fig.text(
        0.055,
        0.94,
        r"0-20% Au+Au vs pp: Region-A yield and ABCD purity",
        ha="left",
        va="top",
        fontsize=27,
        fontweight="bold",
        color="#111827",
    )
    ax_yield.errorbar(
        x,
        raw_y,
        xerr=xerr,
        yerr=raw_yerr,
        fmt="o",
        ms=8.0,
        mfc="#111827",
        mec="white",
        mew=0.8,
        ecolor="#374151",
        elinewidth=1.15,
        capsize=3.0,
        color="#111827",
        label=r"Au+Au Region A ($A$)",
    )
    ax_yield.errorbar(
        x,
        pp_y,
        xerr=xerr,
        yerr=pp_yerr,
        fmt="s",
        ms=8.3,
        mfc="white",
        mec="#DC2626",
        mew=2.0,
        ecolor="#F87171",
        elinewidth=1.15,
        capsize=3.0,
        color="#DC2626",
        label=r"pp PPG12/baseV3E Region A ($A$)",
    )
    ax_yield.set_yscale("log")
    ax_yield.set_ylabel(r"Raw Region-A yield $A$ / GeV")
    ax_yield.set_title("Raw tight+isolated Region-A yield", pad=8, fontweight="bold")
    ax_yield.grid(True, which="major", color="#D9E2EC", alpha=0.85, linewidth=0.85)
    ax_yield.grid(True, which="minor", axis="y", color="#E5E7EB", alpha=0.45, linewidth=0.55)
    ax_yield.legend(
        loc="upper right",
        frameon=True,
        framealpha=0.97,
        facecolor="white",
        edgecolor="#CBD5E1",
        handlelength=1.2,
        borderpad=0.45,
    )
    plt.setp(ax_yield.get_xticklabels(), visible=False)

    ax_purity.errorbar(
        x - 0.08,
        raw_p,
        xerr=xerr,
        yerr=raw_perr,
        fmt="o",
        ms=8.0,
        mfc="#111827",
        mec="white",
        mew=0.8,
        ecolor="#4B5563",
        elinewidth=1.1,
        capsize=3.0,
        color="#111827",
        label="Raw ABCD: (A - BC/D) / A",
    )
    ax_purity.errorbar(
        x + 0.08,
        corr_p,
        xerr=xerr,
        yerr=corr_perr,
        fmt="o",
        ms=8.5,
        mfc="white",
        mec="#2563EB",
        mew=2.0,
        ecolor="#60A5FA",
        elinewidth=1.1,
        capsize=3.0,
        color="#2563EB",
        label="Leakage-corrected ABCD",
    )
    ax_purity.set_title(r"ABCD purity: raw $(A-BC/D)/A$ vs leakage-corrected", pad=8, fontweight="bold")
    ax_purity.set_ylabel("Purity")
    ax_purity.set_xlabel(r"cluster $E_T$ [GeV]", labelpad=8)
    ax_purity.set_ylim(-0.04, 1.08)
    ax_purity.set_xlim(15.5, 35.5)
    ax_purity.set_xticks([16, 18, 20, 22, 24, 26, 28, 30, 32, 34])
    ax_purity.set_yticks(np.arange(0.0, 1.01, 0.2))
    ax_purity.grid(True, which="major", color="#D9E2EC", alpha=0.9, linewidth=0.85)
    ax_purity.legend(
        loc="upper right",
        ncol=1,
        frameon=True,
        framealpha=0.97,
        facecolor="white",
        edgecolor="#CBD5E1",
        handlelength=1.2,
        borderpad=0.45,
    )

    fig.savefig(OUT_PNG, dpi=SLIDE_DPI, facecolor="white")
    plt.close(fig)


def write_manifest(points: list[dict]) -> None:
    manifest = {
        "schema": "THE69_CENT020_RAW_YIELD_RAW_VS_LEAKAGE_CORRECTED_PURITY_V1",
        "output_png": str(OUT_PNG),
        "points_csv": str(OUT_CSV),
        "raw_abcd_csv": str(RAW_ABCD_CSV),
        "leakage_csv": str(LEAKAGE_CSV),
        "pp_exact_bin_data_root": str(PP_EXACT_BIN_DATA_ROOT),
        "pp_trigger_dir": PP_TRIGGER_DIR,
        "centrality": CENT_LABEL,
        "pt_range": [VALID_PT_LO, VALID_PT_HI],
        "yield_definition": "Top panel is raw Region-A yield A per ET bin divided by bin width, where A is tight+isolated. AuAu uses default THE-69 BDT/WP80 and sliding isolation; pp overlay uses exact-bin PPG12/baseV3E table-QA data histograms.",
        "raw_purity_definition": "P_raw=(A - B*C/D)/A using data sideband counts.",
        "corrected_purity_definition": "P_corr=S_A/A where S_A solves S_A=A-(B-fB*S_A)*(C-fC*S_A)/(D-fD*S_A).",
        "default_model": "centAsFeatBase3x3_pt15to35",
        "wp80": "T80(c)=0.53471108+0.0012284143*c, centlinear",
        "isolation": "isoR40_isSliding, Eiso < 7.57 - 0.0658*c, sideGap=0",
        "bins": len(points),
    }
    OUT_MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n")


def main() -> int:
    points = build_points()
    if not points:
        raise RuntimeError("No 0-20% points available for corrected leakage overlay")
    write_points(points)
    write_manifest(points)
    draw(points)
    print(OUT_PNG)
    print(OUT_CSV)
    print(OUT_MANIFEST)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
