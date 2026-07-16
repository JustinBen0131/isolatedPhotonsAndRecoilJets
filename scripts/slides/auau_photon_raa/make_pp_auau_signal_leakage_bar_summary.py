#!/usr/bin/env python3
"""Build a pp-vs-AuAu signal-leakage grouped-bar slide candidate."""

from __future__ import annotations

import csv
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager
from matplotlib.patches import Patch


THIS_FILE = Path(__file__).resolve()
REPO = next((p for p in THIS_FILE.parents if (p / "AGENTS.md").exists()), THIS_FILE.parents[4])
SCRIPTS = REPO / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.append(str(SCRIPTS))

from slides.common.slide_defaults import SLIDE_DPI, slide_figsize  # noqa: E402


PP_ROOT = (
    REPO
    / "dataOutput/ppg12Parity/the76_ppg12_fig24_photonjet_fix_20260702_014217"
    / "final_roots/photonjet/RecoilJets_photonjet5plus10plus20_MERGED.root"
)
AUAU_ROOT = REPO / "InputFiles/the69_leakageCentWP_fix/RecoilJets_embeddedPhoton12plus20_MERGED.root"
OUT_DIR = REPO / "dataOutput/auauPhysicsQA/THE69_leakageCentWP_fix_20260706"
OUT_PNG = OUT_DIR / "pp_auau_signal_leakage_bar_summary.png"
OUT_CSV = OUT_DIR / "pp_auau_signal_leakage_bar_summary.csv"
OUT_JSON = OUT_DIR / "pp_auau_signal_leakage_bar_summary_manifest.json"

TOPDIR = "SIM"
AUAU_ISO_TAG = "isoR40_isSliding"

# Four broad bins aligned to the stored pT-counter edges. The 10-14 and 14-16
# bins are excluded because the current AuAu BDT/leakage counters are not a
# stable populated comparison surface there.
ET_GROUPS = [
    ("16-20", [(16, 18), (18, 20)]),
    ("20-24", [(20, 22), (22, 24)]),
    ("24-26", [(24, 26)]),
    ("26-35", [(26, 35)]),
]

SAMPLES = [
    ("pp", "p+p", "#7c3aed", "pp"),
    ("auau_0_20", "0-20%", "#0284c7", "0_20"),
    ("auau_20_50", "20-50%", "#f59e0b", "20_50"),
    ("auau_50_80", "50-80%", "#15803d", "50_80"),
]

REGIONS = [
    ("B", 2, r"$N_{B}^{sig}/N_{A}^{sig}$", "tight non-isolated leakage"),
    ("C", 3, r"$N_{C}^{sig}/N_{A}^{sig}$", "non-tight isolated leakage"),
    ("D", 4, r"$N_{D}^{sig}/N_{A}^{sig}$", "non-tight non-isolated leakage"),
]

TIMES_DIR = Path("/System/Library/Fonts/Supplemental")
TIMES_FONTS = [
    TIMES_DIR / "Times New Roman.ttf",
    TIMES_DIR / "Times New Roman Bold.ttf",
    TIMES_DIR / "Times New Roman Italic.ttf",
    TIMES_DIR / "Times New Roman Bold Italic.ttf",
]

INK = "#111827"
MUTED = "#475569"
GRID = "#dbe4ef"
PANEL_FILL = "#fbfdff"
EDGE = "#c8d3e2"


@dataclass(frozen=True)
class Count:
    value: float
    err2: float


def configure_fonts() -> None:
    for font in TIMES_FONTS:
        if font.exists():
            font_manager.fontManager.addfont(str(font))
    plt.rcParams.update(
        {
            "font.family": "Times New Roman",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "custom",
            "mathtext.rm": "Times New Roman",
            "mathtext.it": "Times New Roman:italic",
            "mathtext.bf": "Times New Roman:bold",
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "axes.unicode_minus": False,
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


def hist_path(sample_key: str, cent_key: str, lo: int, hi: int) -> str:
    if sample_key == "pp":
        return f"{TOPDIR}/h_sigABCD_MC_pT_{lo}_{hi}"
    return f"{TOPDIR}/h_sigABCD_MC_{AUAU_ISO_TAG}_pT_{lo}_{hi}_cent_{cent_key}"


def get_count(hist, bin_index: int) -> Count:
    value = float(hist.GetBinContent(bin_index))
    error = float(hist.GetBinError(bin_index))
    if error <= 0.0 and value > 0.0:
        error = math.sqrt(value)
    return Count(value=value, err2=error * error)


def sum_counts(counts: list[Count]) -> Count:
    return Count(sum(c.value for c in counts), sum(c.err2 for c in counts))


def ratio(num: Count, den: Count) -> tuple[float, float]:
    if den.value <= 0.0:
        return float("nan"), float("nan")
    r = num.value / den.value
    if num.value <= 0.0:
        return r, math.sqrt(num.err2) / den.value
    rel2 = num.err2 / (num.value * num.value) + den.err2 / (den.value * den.value)
    return r, abs(r) * math.sqrt(max(0.0, rel2))


def collect_rows() -> list[dict]:
    pp_file = open_root(PP_ROOT)
    auau_file = open_root(AUAU_ROOT)
    rows: list[dict] = []
    try:
        for sample_key, sample_label, _color, cent_key in SAMPLES:
            source_file = pp_file if sample_key == "pp" else auau_file
            for group_label, bins in ET_GROUPS:
                paths: list[str] = []
                a_counts: list[Count] = []
                region_counts: dict[str, list[Count]] = {r[0]: [] for r in REGIONS}
                for lo, hi in bins:
                    path = hist_path(sample_key, cent_key, lo, hi)
                    hist = source_file.Get(path)
                    if not hist:
                        raise KeyError(f"Missing histogram: {path}")
                    paths.append(path)
                    a_counts.append(get_count(hist, 1))
                    for region, bin_index, _ylabel, _desc in REGIONS:
                        region_counts[region].append(get_count(hist, bin_index))

                a_total = sum_counts(a_counts)
                for region, _bin_index, ylabel, description in REGIONS:
                    n_total = sum_counts(region_counts[region])
                    value, error = ratio(n_total, a_total)
                    rows.append(
                        {
                            "sample_key": sample_key,
                            "sample_label": sample_label,
                            "et_group": group_label,
                            "region": region,
                            "region_ratio": ylabel,
                            "description": description,
                            "n_a_signal": a_total.value,
                            "n_side_signal": n_total.value,
                            "ratio": value,
                            "ratio_error": error,
                            "histograms": ";".join(paths),
                        }
                    )
    finally:
        pp_file.Close()
        auau_file.Close()
    return rows


def values_for(rows: list[dict], region: str) -> tuple[np.ndarray, np.ndarray]:
    values = np.full((len(SAMPLES), len(ET_GROUPS)), np.nan, dtype=float)
    errors = np.full_like(values, np.nan)
    lookup = {(r["sample_key"], r["et_group"], r["region"]): r for r in rows}
    for i, (sample_key, _sample_label, _color, _cent_key) in enumerate(SAMPLES):
        for j, (group_label, _bins) in enumerate(ET_GROUPS):
            row = lookup[(sample_key, group_label, region)]
            values[i, j] = float(row["ratio"])
            errors[i, j] = float(row["ratio_error"])
    return values, errors


def style_axis(ax: plt.Axes, title: str, ylabel: str, show_x: bool) -> None:
    ax.set_facecolor(PANEL_FILL)
    ax.set_title(title, loc="left", fontsize=17.5, fontweight="bold", pad=7, color=INK)
    ax.set_ylabel(ylabel, fontsize=15.5, color=INK)
    ax.grid(axis="y", color=GRID, linewidth=1.05, zorder=1)
    ax.tick_params(axis="both", labelsize=13.5, colors=INK, length=5, width=1.0)
    ax.tick_params(axis="x", labelbottom=show_x)
    for spine in ax.spines.values():
        spine.set_edgecolor(EDGE)
        spine.set_linewidth(1.0)


def draw_grouped_bars(ax: plt.Axes, rows: list[dict], region: str, title: str, ylabel: str, show_x: bool) -> None:
    vals, errs = values_for(rows, region)
    style_axis(ax, title, ylabel, show_x)
    x = np.arange(len(ET_GROUPS), dtype=float)
    n = len(SAMPLES)
    width = 0.17
    offsets = (np.arange(n) - (n - 1) / 2.0) * width

    for j in range(len(ET_GROUPS)):
        if j % 2 == 0:
            ax.axvspan(j - 0.5, j + 0.5, color="#eef4fb", alpha=0.58, zorder=0)

    for i, (_sample_key, sample_label, color, _cent_key) in enumerate(SAMPLES):
        ax.bar(
            x + offsets[i],
            vals[i],
            width=width * 0.9,
            color=color,
            edgecolor="white",
            linewidth=0.7,
            label=sample_label,
            zorder=3,
        )
        ax.errorbar(
            x + offsets[i],
            vals[i],
            yerr=errs[i],
            fmt="none",
            ecolor="#0f172a",
            elinewidth=0.65,
            capsize=1.8,
            zorder=4,
        )

    finite_max = float(np.nanmax(vals + errs))
    top = max(0.05, finite_max * 1.22)
    if region == "B":
        top = max(top, 0.15)
    elif region == "C":
        top = max(top, 0.36)
    elif region == "D":
        top = max(top, 0.055)
    ax.set_ylim(0.0, top)
    ax.set_xlim(-0.55, len(ET_GROUPS) - 0.45)
    ax.set_xticks(x)
    ax.set_xticklabels([g[0] for g in ET_GROUPS], fontsize=15.0, fontweight="bold")
    if show_x:
        ax.set_xlabel(r"cluster $E_T$ bin (GeV)", fontsize=16.5, fontweight="bold", labelpad=2)


def write_outputs(rows: list[dict]) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "sample_key",
        "sample_label",
        "et_group",
        "region",
        "region_ratio",
        "description",
        "n_a_signal",
        "n_side_signal",
        "ratio",
        "ratio_error",
        "histograms",
    ]
    with OUT_CSV.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    manifest = {
        "schema_version": 1,
        "plot": str(OUT_PNG),
        "csv": str(OUT_CSV),
        "pp_root": str(PP_ROOT),
        "auau_root": str(AUAU_ROOT),
        "histogram_family": {
            "pp": "SIM/h_sigABCD_MC_pT_<lo>_<hi>",
            "auau": f"SIM/h_sigABCD_MC_{AUAU_ISO_TAG}_pT_<lo>_<hi>_cent_<cent>",
            "bins": "1=A, 2=B, 3=C, 4=D",
        },
        "et_groups": [{"label": label, "source_bins": bins} for label, bins in ET_GROUPS],
        "samples": [
            {"key": key, "label": label, "cent_key": cent if key != "pp" else None}
            for key, label, _color, cent in SAMPLES
        ],
        "notes": [
            "Bar heights are sideband truth-signal leakage fractions N_X^sig / N_A^sig.",
            "The plot uses existing local ROOT outputs only; no production, merge, transfer, or Slides mutation was performed.",
            "The 10-14 and 14-16 GeV source bins are excluded because the current AuAu tight-BDT signal ABCD counters are threshold-edge dominated there.",
        ],
    }
    OUT_JSON.write_text(json.dumps(manifest, indent=2) + "\n")


def draw_slide(rows: list[dict]) -> None:
    configure_fonts()
    fig = plt.figure(figsize=slide_figsize(), dpi=SLIDE_DPI, facecolor="white")
    fig.text(
        0.055,
        0.945,
        r"Signal leakage fractions: pp versus AuAu",
        ha="left",
        va="top",
        fontsize=35.0,
        fontweight="bold",
        color=INK,
    )
    fig.text(
        0.055,
        0.892,
        r"truth-signal sideband leakage relative to Region A; current pp photon+jet and THE-69 AuAu centrality output",
        ha="left",
        va="top",
        fontsize=15.3,
        color=MUTED,
    )

    handles = [Patch(facecolor=color, edgecolor="none", label=label) for _key, label, color, _cent in SAMPLES]
    fig.legend(
        handles=handles,
        loc="upper right",
        bbox_to_anchor=(0.945, 0.862),
        frameon=False,
        ncol=4,
        fontsize=17.0,
        columnspacing=1.2,
        handlelength=1.4,
        handletextpad=0.45,
    )

    positions = [
        [0.075, 0.600, 0.875, 0.180],
        [0.075, 0.360, 0.875, 0.180],
        [0.075, 0.120, 0.875, 0.180],
    ]
    for idx, ((region, _bin, ylabel, desc), pos) in enumerate(zip(REGIONS, positions)):
        ax = fig.add_axes(pos)
        draw_grouped_bars(ax, rows, region, desc, ylabel, show_x=(idx == 2))

    fig.savefig(OUT_PNG, dpi=SLIDE_DPI)
    plt.close(fig)


def main() -> int:
    rows = collect_rows()
    write_outputs(rows)
    draw_slide(rows)
    print(f"wrote {OUT_PNG}")
    print(f"wrote {OUT_CSV}")
    print(f"wrote {OUT_JSON}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
