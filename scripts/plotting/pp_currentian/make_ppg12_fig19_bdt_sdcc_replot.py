#!/usr/bin/env python3
"""Replot PPG12 BDT shower-shape panels directly from SDCC ROOT projections.

This avoids using screenshot digitization as the primary proof for BDT panels,
where the data points are compressed near zero and DataThief is fragile.
The input is a compact JSON block produced by a read-only SDCC ROOT extraction.
"""

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
    / "shower_shape_reference_validation/fig13_19_bdt"
)
RAW_EXTRACT = OUTDIR / "ppg12_sdcc_bdt_macro_mode_scan_raw.txt"
JSON_EXTRACT = OUTDIR / "ppg12_sdcc_fig19_bdt_fourcurve_projections.json"
SHOW_CHI2_ON_PNG = False

PANEL_SOURCE_KEYS = {
    "bdt_no_npb_22_28": "bdt_no_npb_22_28__showershape__mixed",
    "bdt_with_npb_18_22": "bdt_with_npb_18_22__showershape__mixed",
}


def load_payload() -> dict:
    text = RAW_EXTRACT.read_text(errors="replace")
    if "JSON_BEGIN" not in text or "JSON_END" not in text:
        raise RuntimeError(f"Could not find JSON block in {RAW_EXTRACT}")
    scan_payload = json.loads(text.split("JSON_BEGIN", 1)[1].split("JSON_END", 1)[0].strip())
    panels = {}
    for tag, scan_key in PANEL_SOURCE_KEYS.items():
        panel = scan_payload["panels"][scan_key]
        panel = dict(panel)
        panel["center_note"] = "no NPB cut" if tag == "bdt_no_npb_22_28" else "With NPB cut"
        panel["label"] = pretty_label(tag, "")
        panels[tag] = panel
    payload = {
        "normalization": (
            "Exact plot_showershapes_variations.C macro mode: showershape suffix, mixed-mode "
            "combined signal/inclusive ROOT files, ProjectionX after RebinX(2), x-axis 0..1, "
            "unit-normalized curves; no-NPB NPB overlay tail-scaled with weta_cogx > 1.4."
        ),
        "panels": panels,
        "source_files": {
            tag: panels[tag]["summary"]["source_files"] for tag in panels
        },
        "scan_raw": str(RAW_EXTRACT),
        "source_selection": PANEL_SOURCE_KEYS,
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


def draw_step(ax, curve: dict[str, np.ndarray], *, color: str, label: str, lw: float = 1.7) -> None:
    vals = curve["values"]
    y = np.r_[vals, vals[-1]]
    ax.step(curve["edges"], y, where="post", color=color, lw=lw, label=label)


def pretty_label(tag: str, raw: str) -> str:
    if tag == "bdt_no_npb_22_28":
        return r"$22<p_T<28$ GeV, w/o nbkg cut"
    if tag == "bdt_with_npb_18_22":
        return r"$18<p_T<22$ GeV, w/ nbkg cut"
    return raw.replace("p_{T}", r"$p_T$")


def rootish_tick(value: float, _pos: int) -> str:
    if abs(value) < 1e-12:
        return "0"
    return f"{value:.2f}".rstrip("0").rstrip(".")


def apply_ppg12_axis_ranges(tag: str, ax, rax, panel: dict) -> None:
    if tag == "bdt_no_npb_22_28":
        ax.set_ylim(0.0, 0.85)
        ax.yaxis.set_major_locator(FixedLocator(np.arange(0.0, 0.8001, 0.1)))
        ax.yaxis.set_minor_locator(MultipleLocator(0.05))
        ax.yaxis.set_major_formatter(FuncFormatter(rootish_tick))

        rax.set_ylim(-0.5, 0.5)
        rax.yaxis.set_major_locator(FixedLocator([-0.5, 0.0, 0.5]))
        rax.yaxis.set_minor_locator(MultipleLocator(0.1))
        rax.yaxis.set_major_formatter(FuncFormatter(rootish_tick))
        return

    if tag == "bdt_with_npb_18_22":
        ax.set_ylim(0.0, 0.235)
        ax.yaxis.set_major_locator(FixedLocator(np.arange(0.0, 0.2201, 0.02)))
        ax.yaxis.set_minor_locator(MultipleLocator(0.01))
        ax.yaxis.set_major_formatter(FuncFormatter(rootish_tick))

        rax.set_ylim(-0.03, 0.03)
        rax.yaxis.set_major_locator(FixedLocator([-0.02, 0.0, 0.02]))
        rax.yaxis.set_minor_locator(MultipleLocator(0.01))
        rax.yaxis.set_major_formatter(FuncFormatter(rootish_tick))
        return

    ax.set_ylim(panel["ymin"], panel["ymax"])
    rax.set_ylim(panel["diff_ymin"], panel["diff_ymax"])


def draw_panel(payload: dict, tag: str) -> Path:
    panel = payload["panels"][tag]
    curves = {name: as_arrays(curve) for name, curve in panel["curves"].items()}

    setup_style()
    out = OUTDIR / f"ppg12_sdcc_root_replot_{tag}_fourcurve_slidefit_772x998.png"
    fig, (ax, rax) = plt.subplots(
        2,
        1,
        figsize=(7.72, 9.98),
        dpi=200,
        sharex=True,
        gridspec_kw={"height_ratios": [3.35, 1.0], "hspace": 0.035},
    )

    data = curves["data"]
    sig = curves["signal_mc"]
    incl = curves["inclusive_mc"]
    npb = curves.get("npb_tagged_data")

    draw_step(ax, sig, color="red", label="Signal MC")
    draw_step(ax, incl, color="blue", label="Inclusive MC")
    if npb is not None:
        draw_step(ax, npb, color="green", label="NPB-tagged data")
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

    ax.text(0.05, 0.96, "sPHENIX", transform=ax.transAxes, ha="left", va="top", fontsize=18, fontweight="bold", fontstyle="italic")
    ax.text(0.265, 0.96, "Internal", transform=ax.transAxes, ha="left", va="top", fontsize=18)
    ax.text(0.05, 0.885, r"$p{+}p\ \sqrt{s}=200$ GeV", transform=ax.transAxes, fontsize=17, ha="left", va="top")
    ax.text(0.05, 0.815, r"$|\eta^\gamma| < 0.7$", transform=ax.transAxes, fontsize=17, ha="left", va="top")
    ax.text(0.05, 0.748, pretty_label(tag, panel["label"]), transform=ax.transAxes, fontsize=14.5, ha="left", va="top")
    ax.text(0.48, 0.48, panel["center_note"], transform=ax.transAxes, fontsize=20, fontweight="bold", ha="center", va="center")
    ax.text(
        0.48,
        0.40,
        "Replotted directly from PPG12 SDCC ROOT",
        transform=ax.transAxes,
        fontsize=12.5,
        ha="center",
        va="center",
    )

    chi2 = panel.get("chi2")
    ndf = panel.get("ndf")
    pvalue = panel.get("pvalue")
    if SHOW_CHI2_ON_PNG and chi2 is not None and ndf:
        ax.text(0.05, 0.675, rf"$\chi^2$/ndf = {chi2:.1f}/{ndf:d} = {chi2 / ndf:.2f}", transform=ax.transAxes, fontsize=12.5, ha="left", va="top")
    if SHOW_CHI2_ON_PNG and pvalue is not None:
        ax.text(0.05, 0.62, f"p-value = {pvalue:.4f}", transform=ax.transAxes, fontsize=12.5, ha="left", va="top")

    ax.set_ylabel("normalized counts", fontsize=20)
    ax.set_xlim(0.0, 1.0)
    ax.tick_params(labelsize=16, top=True, right=True)
    ax.minorticks_on()
    handles, labels = ax.get_legend_handles_labels()
    order = [labels.index(name) for name in ["Data", "Signal MC", "Inclusive MC"] if name in labels]
    if "NPB-tagged data" in labels:
        order.append(labels.index("NPB-tagged data"))
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

    diff = data["values"] - incl["values"]
    diff_err = np.sqrt(data["errors"] ** 2 + incl["errors"] ** 2)
    rax.errorbar(data["centers"], diff, yerr=diff_err, fmt="o", color="black", ms=3.6, lw=0.9)
    rax.axhline(0.0, color="black", ls="--", lw=1.0)
    rax.set_ylabel("Data - Incl. MC", fontsize=17)
    rax.set_xlabel("bdt", fontsize=20)
    rax.set_xlim(0.0, 1.0)
    apply_ppg12_axis_ranges(tag, ax, rax, panel)
    rax.tick_params(labelsize=16, top=True, right=True)
    rax.minorticks_on()

    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)

    manifest = {
        "artifact": str(out),
        "compact_json": str(JSON_EXTRACT),
        "raw_extract": str(RAW_EXTRACT),
        "panel": panel,
        "source_files": payload["source_files"],
        "source_selection": payload["source_selection"],
        "note": (
            "Direct PPG12 SDCC ROOT replot using the exact plot_showershapes_variations.C "
            "mixed-mode source family selected by IAN chi2/ndf and visible peak fingerprints; "
            "no screenshot digitization used for plotted curves."
        ),
    }
    out.with_suffix(".manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True))
    return out


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    payload = load_payload()
    outputs = [draw_panel(payload, tag) for tag in payload["panels"]]
    print("\n".join(str(p) for p in outputs))


if __name__ == "__main__":
    main()
