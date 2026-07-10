#!/usr/bin/env python3
"""Replot PPG12 IAN Fig. 20 BDT panel directly from SDCC ROOT projections."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FixedLocator, FuncFormatter, MultipleLocator


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
OUTDIR = (
    REPO
    / "dataOutput/ppg12PhotonYield/ppg12_photon_yield_v1_data_20260620"
    / "shower_shape_reference_validation/fig20_bdt"
)
RAW_EXTRACT = OUTDIR / "ppg12_sdcc_fig20_bdt_tight_macro_mode_scan_raw.txt"
JSON_EXTRACT = OUTDIR / "ppg12_sdcc_fig20_bdt_tight_fourcurve_projections.json"
PANEL_SOURCE_KEY = "bdt_tight_10_14__showershape__mixed"


def load_payload() -> dict:
    text = RAW_EXTRACT.read_text(errors="replace")
    if "JSON_BEGIN" not in text or "JSON_END" not in text:
        raise RuntimeError(f"Could not find JSON block in {RAW_EXTRACT}")
    scan_payload = json.loads(text.split("JSON_BEGIN", 1)[1].split("JSON_END", 1)[0].strip())
    panel = dict(scan_payload["panels"][PANEL_SOURCE_KEY])
    panel["center_note"] = "With tight cut"
    panel["label"] = r"$10<p_T<14$ GeV, w/ tight cut"
    panel["ymin"] = 0.0
    panel["ymax"] = 0.285
    panel["diff_ymin"] = -0.04
    panel["diff_ymax"] = 0.04
    payload = {
        "normalization": (
            "Exact plot_showershapes_variations.C macro mode: showershape suffix, "
            "mixed-mode combined signal/inclusive ROOT files, h2d_bdt_eta0_pt0_cut2, "
            "ProjectionX after RebinX(2), x-axis 0..1, and unit-normalized curves."
        ),
        "panel": panel,
        "scan_raw": str(RAW_EXTRACT),
        "source_selection": PANEL_SOURCE_KEY,
        "target": scan_payload["target"],
    }
    JSON_EXTRACT.write_text(json.dumps(payload, indent=2, sort_keys=True))
    return payload


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "dejavuserif",
            "axes.linewidth": 1.2,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.major.size": 7,
            "ytick.major.size": 7,
            "xtick.minor.size": 3.5,
            "ytick.minor.size": 3.5,
        }
    )


def as_arrays(curve: dict) -> dict[str, np.ndarray]:
    return {
        "edges": np.asarray(curve["edges"], dtype=float),
        "centers": np.asarray(curve["centers"], dtype=float),
        "values": np.asarray(curve["values"], dtype=float),
        "errors": np.asarray(curve["errors"], dtype=float),
    }


def draw_step(ax, curve: dict[str, np.ndarray], *, color: str, label: str) -> None:
    vals = curve["values"]
    ax.step(curve["edges"], np.r_[vals, vals[-1]], where="post", color=color, lw=1.7, label=label)


def rootish_tick(value: float, _pos: int) -> str:
    if abs(value) < 1e-12:
        return "0"
    return f"{value:.2f}".rstrip("0").rstrip(".")


def draw_panel(payload: dict) -> Path:
    panel = payload["panel"]
    curves = {name: as_arrays(curve) for name, curve in panel["curves"].items()}

    setup_style()
    out = OUTDIR / "ppg12_sdcc_root_replot_bdt_tight_10_14_fourcurve_slidefit_772x998.png"
    fig, (ax, rax) = plt.subplots(
        2,
        1,
        figsize=(7.72, 9.98),
        dpi=200,
        sharex=True,
        gridspec_kw={"height_ratios": [3.35, 1.0], "hspace": 0.035},
    )

    data = curves["data"]
    signal = curves["signal_mc"]
    inclusive = curves["inclusive_mc"]

    draw_step(ax, signal, color="red", label="Signal MC")
    draw_step(ax, inclusive, color="blue", label="Inclusive MC")
    ax.axhline(
        float(np.nanmax(signal["values"])),
        color="black",
        lw=0.9,
        zorder=2,
    )
    ax.errorbar(
        data["centers"],
        data["values"],
        yerr=data["errors"],
        fmt="o",
        color="black",
        ms=4.4,
        lw=1.0,
        label="Data",
        zorder=4,
    )

    ax.text(
        0.05,
        0.96,
        "sPHENIX",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=18,
        fontweight="bold",
        fontstyle="italic",
    )
    ax.text(0.265, 0.96, "Internal", transform=ax.transAxes, ha="left", va="top", fontsize=18)
    ax.text(0.05, 0.885, r"$p{+}p\ \sqrt{s}=200$ GeV", transform=ax.transAxes, fontsize=17, ha="left", va="top")
    ax.text(0.05, 0.815, r"$|\eta^\gamma| < 0.7$", transform=ax.transAxes, fontsize=17, ha="left", va="top")
    ax.text(0.05, 0.748, panel["label"], transform=ax.transAxes, fontsize=14.5, ha="left", va="top")
    ax.text(0.50, 0.48, panel["center_note"], transform=ax.transAxes, fontsize=20, fontweight="bold", ha="center", va="center")
    ax.text(
        0.50,
        0.40,
        "Replotted directly from PPG12 SDCC ROOT",
        transform=ax.transAxes,
        fontsize=12.5,
        ha="center",
        va="center",
    )

    ax.set_ylabel("normalized counts", fontsize=20)
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(panel["ymin"], panel["ymax"])
    ax.yaxis.set_major_locator(FixedLocator(np.arange(0.0, 0.2501, 0.05)))
    ax.yaxis.set_minor_locator(MultipleLocator(0.01))
    ax.yaxis.set_major_formatter(FuncFormatter(rootish_tick))
    ax.tick_params(labelsize=16, top=True, right=True)
    ax.minorticks_on()

    handles, labels = ax.get_legend_handles_labels()
    order = [labels.index(name) for name in ["Data", "Signal MC", "Inclusive MC"] if name in labels]
    ax.legend(
        [handles[i] for i in order],
        [labels[i] for i in order],
        loc="upper right",
        frameon=False,
        fontsize=16,
        handlelength=1.6,
        borderpad=0.2,
        labelspacing=0.38,
    )

    diff = data["values"] - inclusive["values"]
    diff_err = np.sqrt(data["errors"] ** 2 + inclusive["errors"] ** 2)
    rax.errorbar(data["centers"], diff, yerr=diff_err, fmt="o", color="black", ms=3.6, lw=0.9)
    rax.axhline(0.0, color="black", ls="--", lw=1.0)
    rax.set_ylabel("Data - Incl. MC", fontsize=17)
    rax.set_xlabel("bdt", fontsize=20)
    rax.set_xlim(0.0, 1.0)
    rax.set_ylim(panel["diff_ymin"], panel["diff_ymax"])
    rax.yaxis.set_major_locator(FixedLocator([-0.04, -0.02, 0.0, 0.02, 0.04]))
    rax.yaxis.set_minor_locator(MultipleLocator(0.01))
    rax.yaxis.set_major_formatter(FuncFormatter(rootish_tick))
    rax.xaxis.set_major_locator(FixedLocator(np.arange(0.0, 1.0001, 0.2)))
    rax.xaxis.set_minor_locator(MultipleLocator(0.05))
    rax.xaxis.set_major_formatter(FuncFormatter(rootish_tick))
    rax.tick_params(labelsize=16, top=True, right=True)
    rax.minorticks_on()

    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)

    manifest = {
        "artifact": str(out),
        "compact_json": str(JSON_EXTRACT),
        "raw_extract": str(RAW_EXTRACT),
        "source_selection": PANEL_SOURCE_KEY,
        "panel": panel,
        "note": (
            "Direct PPG12 SDCC ROOT replot of IAN Fig.20 BDT tight panel. "
            "The selected source key matches the IAN chi2/ndf fingerprint: "
            f"{panel['chi2']:.2f}/{panel['ndf']} = {panel['chi2'] / panel['ndf']:.2f}."
        ),
    }
    out.with_suffix(".manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True))
    return out


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    payload = load_payload()
    print(draw_panel(payload))


if __name__ == "__main__":
    main()
