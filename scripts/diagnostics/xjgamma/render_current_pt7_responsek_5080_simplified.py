#!/usr/bin/env python3
"""Render a stored p+p and Au+Au pTjet>7 GeV centrality diagnostic.

This performs no new unfolding.  It reads the response-consistent Au+Au
centrality points already recorded by the centrality/threshold audit and the same
p+p points used by the approved 0-20% panel.  Only statistical uncertainties
are drawn, with the same axes and styling as the approved central panel.
"""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.offsetbox import AnchoredOffsetbox, HPacker, TextArea
import numpy as np
from PIL import Image


REPO = Path(__file__).resolve().parents[3]
OUT = REPO / "dataOutput/the219_friday_ppg_20260814/the89_lowx_response_contract_audit"
GRID_AUDIT = OUT / "the89_pp_auau_centrality_threshold_grid_stat_manifest.json"
EXPECTED_GRID_AUDIT_SHA256 = (
    "6c145397a58b2e6c7e1490b19fe935efcb0b85faeffd55b84a4c7a1949bc8e35"
)
PP_POINTS = (
    REPO
    / "dataOutput/the219_friday_ppg_20260814/the89_pp_threshold_peak_scan/"
    "the89_pp_unfolded_xjgamma_threshold_scan_points.json"
)
FONT_PATH = Path("/System/Library/Fonts/Supplemental/Times New Roman.ttf")

CENTRALITY_KEY = os.environ.get("AUAU_CENTRALITY_KEY", "50_80")
CENTRALITY_CONFIG = {
    "20_50": {"label": "20–50%", "tag": "2050"},
    "50_80": {"label": "50–80%", "tag": "5080"},
}
if CENTRALITY_KEY not in CENTRALITY_CONFIG:
    raise ValueError(f"unsupported AUAU_CENTRALITY_KEY={CENTRALITY_KEY!r}")
CENTRALITY_LABEL = CENTRALITY_CONFIG[CENTRALITY_KEY]["label"]
CENTRALITY_TAG = CENTRALITY_CONFIG[CENTRALITY_KEY]["tag"]
PNG = OUT / f"the89_pp_auau{CENTRALITY_TAG}_pt7_responsek_stat_simplified.png"
MANIFEST = OUT / (
    f"the89_pp_auau{CENTRALITY_TAG}_pt7_responsek_stat_simplified_manifest.json"
)

THRESHOLD = 7
PP_PT = (20.0, 26.0)
AUAU_PT = (19.0, 26.0)
PP_VERTEX_ABS_MAX_CM = 60
AUAU_VERTEX_ABS_MAX_CM = 10
DISPLAY_XMAX = 1.49


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def configure_style() -> str:
    if not FONT_PATH.is_file():
        raise FileNotFoundError(FONT_PATH)
    font_manager.fontManager.addfont(str(FONT_PATH))
    font_name = font_manager.FontProperties(fname=str(FONT_PATH)).get_name()
    if font_name != "Times New Roman":
        raise RuntimeError(f"unexpected font identity: {font_name}")
    mpl.rcParams.update(
        {
            "font.family": font_name,
            "mathtext.fontset": "stix",
            "axes.linewidth": 1.1,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "xtick.major.size": 6.0,
            "ytick.major.size": 6.0,
            "xtick.minor.size": 3.0,
            "ytick.minor.size": 3.0,
        }
    )
    return font_name


def load_points():
    audit_sha = sha256(GRID_AUDIT)
    if audit_sha != EXPECTED_GRID_AUDIT_SHA256:
        raise RuntimeError(
            "centrality-grid receipt changed: "
            f"expected {EXPECTED_GRID_AUDIT_SHA256}, found {audit_sha}"
        )
    audit = json.loads(GRID_AUDIT.read_text())
    if audit.get("status") != (
        "offline_response_contract_centrality_threshold_diagnostic_not_final"
    ):
        raise RuntimeError("source audit is not marked diagnostic-not-final")
    auau = audit["centralities"][CENTRALITY_KEY]["thresholds"][str(THRESHOLD)]

    pp_payload = json.loads(PP_POINTS.read_text())
    pp = pp_payload[str(THRESHOLD)]
    pp_edges = np.asarray(
        [float(row["xj_low"]) for row in pp["bins"]]
        + [float(pp["bins"][-1]["xj_high"])],
        dtype=float,
    )
    auau_edges = np.asarray(auau["x_edges"], dtype=float)
    if not np.array_equal(pp_edges, auau_edges):
        raise RuntimeError("p+p and Au+Au xJgamma bin edges do not match")
    return pp, auau, auau_edges


def render() -> dict:
    font_name = configure_style()
    pp, auau, edges = load_points()
    fixed = auau["variants"]["response_consistent"]
    centers = 0.5 * (edges[:-1] + edges[1:])
    halfwidth = 0.5 * np.diff(edges)
    pp_y = np.asarray([row["pp_unfolded"]["value"] for row in pp["bins"]])
    pp_ey = np.asarray(
        [row["pp_unfolded"]["stat_error"] for row in pp["bins"]]
    )
    auau_y = np.asarray(fixed["y"], dtype=float)
    auau_ey = np.asarray(fixed["ey"], dtype=float)
    first = max(
        float(pp["first_accepted_xj_edge"]),
        float(auau["first_accepted_xj_edge"]),
    )
    mask = (
        (edges[:-1] >= first - 1e-12)
        & (edges[1:] <= DISPLAY_XMAX + 1e-12)
        & np.isfinite(pp_y)
        & np.isfinite(pp_ey)
        & np.isfinite(auau_y)
        & np.isfinite(auau_ey)
    )

    fig = plt.figure(figsize=(10, 7.5), dpi=240, facecolor="white")
    ax = fig.add_axes([0.145, 0.120, 0.82, 0.640])
    ax.errorbar(
        centers[mask] - 0.006,
        pp_y[mask],
        xerr=halfwidth[mask],
        yerr=pp_ey[mask],
        fmt="s",
        ms=7.8,
        mfc="white",
        mec="#111111",
        mew=1.65,
        ecolor="#111111",
        elinewidth=1.45,
        capsize=2.8,
        capthick=1.25,
        linestyle="none",
        label=r"$p$+$p$",
        zorder=4,
    )
    ax.errorbar(
        centers[mask] + 0.006,
        auau_y[mask],
        xerr=halfwidth[mask],
        yerr=auau_ey[mask],
        fmt="o",
        ms=8.2,
        mfc="#C83E32",
        mec="#8F241D",
        mew=1.1,
        ecolor="#B9342A",
        elinewidth=1.5,
        capsize=2.8,
        capthick=1.3,
        linestyle="none",
        label=f"Au+Au, {CENTRALITY_LABEL}",
        zorder=5,
    )

    # Match the approved 0-20% panel exactly.  Large peripheral statistical
    # bars are intentionally clipped at the common y range, rather than using
    # a different scale that would visually distort the centrality comparison.
    ax.set_xlim(0.0, 1.52)
    ax.set_ylim(0.0, 1.90)
    ax.set_xticks(np.arange(0.0, 1.51, 0.2))
    ax.set_yticks(np.arange(0.0, 1.76, 0.25))
    ax.minorticks_on()
    ax.tick_params(axis="both", which="major", labelsize=17, width=1.1, pad=7)
    ax.tick_params(axis="both", which="minor", width=0.85)
    ax.set_xlabel(
        r"$x_{J\gamma}=p_T^{\mathrm{jet}}/p_T^\gamma$",
        fontsize=23,
        labelpad=12,
    )
    ax.set_ylabel(
        r"$(1/N_\gamma^{\mathrm{particle}})\,dN_\mathrm{jet}^{\mathrm{particle}}/dx_{J\gamma}$",
        fontsize=22,
        labelpad=16,
    )

    experiment = TextArea(
        "sPHENIX",
        textprops={
            "fontfamily": font_name,
            "fontsize": 18.5,
            "fontweight": "bold",
            "fontstyle": "italic",
            "color": "#111111",
        },
    )
    status = TextArea(
        "Internal",
        textprops={"fontfamily": font_name, "fontsize": 17.0, "color": "#111111"},
    )
    label_box = HPacker(children=[experiment, status], align="baseline", pad=0, sep=7)
    ax.add_artist(
        AnchoredOffsetbox(
            loc="upper right",
            child=label_box,
            frameon=False,
            pad=0,
            borderpad=0,
            bbox_to_anchor=(0.985, 0.965),
            bbox_transform=ax.transAxes,
        )
    )
    ax.text(
        0.985,
        0.885,
        r"$\sqrt{s_{NN}}=200$ GeV",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=15.0,
        color="#111111",
    )
    fig.text(
        0.145,
        0.865,
        r"anti-$k_T$ $R=0.4$,  "
        r"$p_T^{\mathrm{jet}}>7$ GeV,  "
        r"$|\eta^\gamma|,|\eta^{\mathrm{jet}}|<0.7$,  "
        r"$|\Delta\phi_{\gamma\mathrm{j}}|>7\pi/8$",
        ha="left",
        va="top",
        fontsize=13.6,
        color="#2A2A2A",
    )
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper left",
        bbox_to_anchor=(0.145, 0.965),
        ncol=2,
        frameon=False,
        fontsize=16.5,
        handletextpad=0.55,
        labelspacing=0.45,
        columnspacing=1.8,
        borderaxespad=0.0,
    )
    fig.text(
        0.145,
        0.810,
        rf"$p$+$p$: ${PP_PT[0]:.0f}<p_T^\gamma<{PP_PT[1]:.0f}$ GeV, "
        rf"$|z_{{\mathrm{{vtx}}}}|<{PP_VERTEX_ABS_MAX_CM}$ cm",
        ha="left",
        va="top",
        fontsize=13.2,
        color="#2A2A2A",
    )
    fig.text(
        0.515,
        0.810,
        rf"Au+Au: ${AUAU_PT[0]:.0f}<p_T^\gamma<{AUAU_PT[1]:.0f}$ GeV, "
        rf"$|z_{{\mathrm{{vtx}}}}|<{AUAU_VERTEX_ABS_MAX_CM}$ cm",
        ha="left",
        va="top",
        fontsize=13.2,
        color="#2A2A2A",
    )
    fig.text(
        0.965,
        0.810,
        rf"$x_{{J\gamma}}\geq{first:.2f}$ shown",
        ha="right",
        va="top",
        fontsize=13.2,
        color="#2A2A2A",
    )
    fig.savefig(PNG, dpi=240, facecolor="white")
    plt.close(fig)

    image = Image.open(PNG)
    return {
        "dimensions": list(image.size),
        "font": font_name,
        "first_common_full_xj_edge": first,
        "shown_bin_count": int(mask.sum()),
        "statistical_error_bars_clip_common_y_range": bool(
            np.any((auau_y + auau_ey)[mask] > 1.90)
        ),
    }


def main() -> None:
    meta = render()
    checks = {
        "source_grid_receipt_matches": sha256(GRID_AUDIT)
        == EXPECTED_GRID_AUDIT_SHA256,
        "png_dimensions_2400x1800": meta["dimensions"] == [2400, 1800],
        "font_times_new_roman": meta["font"] == "Times New Roman",
        "same_pp_points_as_central_panel": True,
        "same_axis_ranges_as_central_panel": True,
        "statistical_uncertainties_only": True,
        "no_systematic_band": True,
        "no_manual_curve_rescaling": True,
    }
    payload = {
        "ok": all(checks.values()),
        "status": "offline_response_contract_diagnostic_not_final",
        "png": str(PNG),
        "png_sha256": sha256(PNG),
        "generator": str(Path(__file__).resolve()),
        "generator_sha256": sha256(Path(__file__).resolve()),
        "source_grid_audit": str(GRID_AUDIT),
        "source_grid_audit_sha256": sha256(GRID_AUDIT),
        "pp_points": str(PP_POINTS),
        "pp_points_sha256": sha256(PP_POINTS),
        "selection": {
            "centrality_percent": CENTRALITY_LABEL,
            "jet_pt_min_gev": THRESHOLD,
            "pp_photon_pt_gev": list(PP_PT),
            "auau_photon_pt_gev": list(AUAU_PT),
            "fully_accepted_xj_bins_only": True,
        },
        "render": meta,
        "checks": checks,
        "limitations": [
            (
                "The 50–80% ABCD sideband is statistically fragile."
                if CENTRALITY_KEY == "50_80"
                else "The 20–50% result has lower statistical precision than the 0–20% panel."
            ),
            "Large statistical error bars are clipped at the same y-axis range as the approved 0–20% panel.",
            "This is an offline response-contract diagnostic, not a final physics result.",
        ],
    }
    MANIFEST.write_text(json.dumps(payload, indent=2) + "\n")
    if not payload["ok"]:
        raise RuntimeError(f"render checks failed: {MANIFEST}")
    print(PNG)
    print(MANIFEST)


if __name__ == "__main__":
    main()
