#!/usr/bin/env python3
"""Render an ATLAS Figure-1-style reconstruction/subtraction comparison.

The ATLAS panels are exact raster crops from the published paper.  The
sPHENIX panels are rebuilt from the current registered p+p and Au+Au ROOT
inputs with the same correction helpers and selections associated with the
linked THE-89 diagnostic.  This slide deliberately stops *before* unfolding:
it compares the high-level construction of the nominal corrected data input,
not collision energies, absolute yields, or final physics results.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys
from typing import Any

import matplotlib as mpl
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import ROOT


ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

REPO = Path(__file__).resolve().parents[3]
OUT = REPO / "dataOutput/the219_friday_ppg_20260814/the89_atlas_sphenix_reco_subtraction_pt7"
OUT.mkdir(parents=True, exist_ok=True)

PNG = OUT / "the89_atlas_sphenix_reco_subtraction_pt7_slide_atlas_selection.png"
MANIFEST = OUT / "the89_atlas_sphenix_reco_subtraction_pt7_atlas_selection_manifest.json"
POINTS = OUT / "the89_atlas_sphenix_reco_subtraction_pt7_atlas_selection_points.json"
LAYOUT = OUT / "the89_atlas_sphenix_reco_subtraction_pt7_atlas_selection_layout_nodes.json"
NOTES = OUT / "the89_atlas_sphenix_reco_subtraction_pt7_atlas_selection_speaker_script.md"

ATLAS_DIR = REPO / "dataOutput/the85_auau_xjgamma_unfolding_push/atlas_reference"
ATLAS_PP = ATLAS_DIR / "atlas_fig1_pp_panel_crop_hi_tight.png"
ATLAS_AUAU = ATLAS_DIR / "atlas_fig1_pbpb_0_10_panel_crop_hi_tight.png"
ATLAS_PDF = REPO / "usefulDocs/external_references/gamma_jet/ATLAS_xJgamma.pdf"
ATLAS_DOI = "https://doi.org/10.1016/j.physletb.2018.12.023"

LINKED_MANIFEST = (
    REPO
    / "dataOutput/the219_friday_ppg_20260814/the89_lowx_response_contract_audit"
    / "the89_pp_auau020_pt7_responsek_stat_simplified_manifest.json"
)
IMMUTABLE_AUDIT = (
    REPO
    / "dataOutput/the219_friday_ppg_20260814/the89_lowx_response_contract_audit"
    / "the89_current_lowx_response_contract_pt5_7_10_manifest.json"
)

PP_DATA_POINTER = REPO / "dataOutput/current_recoiljets_artifacts/current/pp_data_merged/current.json"
PP_SIM_POINTER = REPO / "dataOutput/current_recoiljets_artifacts/current/pp_sim_photonjet_merged/current.json"
AUAU_DATA_POINTER = REPO / "dataOutput/current_recoiljets_artifacts/current/auau_data_merged/current.json"
AUAU_SIM_POINTER = REPO / "dataOutput/current_recoiljets_artifacts/current/auau_sim_photonjet_merged/current.json"

PP_HELPER_PATH = REPO / "scripts/slides/pp_currentian/xjgamma/make_the219_pp_money_slide.py"
AUAU_AUDIT_PATH = REPO / "scripts/slides/the85_xjgamma_unfolding/make_the89_current_auau_xjgamma_quality_audit.py"

PP_PT = (20.0, 26.0)
AUAU_PT = (19.0, 26.0)
JET_PT_MIN = 7
DISPLAY_X = (0.2, 1.8)
EXPECTED_IMMUTABLE_AUDIT_SHA = "c2781ffe2cfd72b7997f507c6eca3185ee5320d6e6c64a755dfee70b376c84cf"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def current_root(pointer_path: Path) -> tuple[Path, dict[str, Any]]:
    payload = json.loads(pointer_path.read_text())
    roots = payload.get("root_paths", [])
    if len(roots) != 1:
        raise RuntimeError(f"expected one root path in {pointer_path}, found {len(roots)}")
    root = Path(roots[0])
    if not root.is_file():
        raise FileNotFoundError(root)
    return root, payload


def open_root(path: Path):
    handle = ROOT.TFile.Open(str(path), "READ")
    if not handle or handle.IsZombie():
        raise RuntimeError(f"failed to open ROOT file: {path}")
    return handle


def axis_edges(axis) -> np.ndarray:
    return np.asarray(
        [float(axis.GetBinLowEdge(i)) for i in range(1, axis.GetNbins() + 2)],
        dtype=float,
    )


def exact_pt_bins(axis, lo: float, hi: float) -> list[int]:
    bins: list[int] = []
    for ib in range(1, axis.GetNbins() + 1):
        blo = float(axis.GetBinLowEdge(ib))
        bhi = float(axis.GetBinUpEdge(ib))
        if blo >= lo - 1.0e-9 and bhi <= hi + 1.0e-9:
            bins.append(ib)
    covered = sum(float(axis.GetBinWidth(ib)) for ib in bins)
    if not bins or abs(covered - (hi - lo)) > 1.0e-6:
        raise RuntimeError(f"pT window {lo:g}-{hi:g} is not an exact union of native bins")
    return bins


def project_counts(h2, pt_window: tuple[float, float]) -> dict[str, np.ndarray]:
    bins = exact_pt_bins(h2.GetXaxis(), *pt_window)
    edges = axis_edges(h2.GetYaxis())
    values = np.zeros(len(edges) - 1, dtype=float)
    errors2 = np.zeros(len(edges) - 1, dtype=float)
    for ix in bins:
        for iy in range(1, h2.GetYaxis().GetNbins() + 1):
            values[iy - 1] += float(h2.GetBinContent(ix, iy))
            errors2[iy - 1] += float(h2.GetBinError(ix, iy)) ** 2
    return {"edges": edges, "values": values, "errors": np.sqrt(errors2), "pt_bins": np.asarray(bins)}


def ensure_matching_edges(*items: dict[str, np.ndarray]) -> np.ndarray:
    reference = items[0]["edges"]
    for item in items[1:]:
        if not np.array_equal(reference, item["edges"]):
            raise RuntimeError("component xJ bin edges do not match")
    return reference


def build_pp() -> tuple[dict[str, Any], dict[str, Any]]:
    pp = load_module("the89_pp_current_helper", PP_HELPER_PATH)
    data_path, data_pointer = current_root(PP_DATA_POINTER)
    sim_path, sim_pointer = current_root(PP_SIM_POINTER)
    data = open_root(data_path)
    sim = open_root(sim_path)
    try:
        base = "r04_jetPt7"
        h_a = pp.obj(data, pp.DATA_TOP, f"h2_unfoldReco_pTgamma_xJ_incl_{base}", "TH2")
        h_c = pp.obj(data, pp.DATA_TOP, f"h2_unfoldReco_pTgamma_xJ_incl_sidebandC_{base}", "TH2")
        purity, rows = pp.purity_correct_xj(data, sim, h_a, h_c)
        raw = project_counts(h_a, PP_PT)
        corrected = project_counts(purity, PP_PT)
    finally:
        data.Close()
        sim.Close()

    edges = ensure_matching_edges(raw, corrected)
    photon_background = raw["values"] - corrected["values"]
    combinatoric = np.zeros_like(photon_background)
    closure = raw["values"] - photon_background - combinatoric - corrected["values"]
    selected_rows = [row for row in rows if row["pt"][0] >= PP_PT[0] and row["pt"][1] <= PP_PT[1]]
    result = {
        "key": "pp",
        "collision": "p+p",
        "pt_window": list(PP_PT),
        "edges": edges,
        "raw": raw["values"],
        "raw_error": raw["errors"],
        "photon_background": photon_background,
        "combinatoric": combinatoric,
        "corrected": corrected["values"],
        "corrected_error": corrected["errors"],
        "additivity_residual": closure,
        "abcd_rows": selected_rows,
    }
    provenance = {
        "data_pointer": str(PP_DATA_POINTER),
        "data_pointer_sha256": sha256(PP_DATA_POINTER),
        "data_campaign": data_pointer.get("campaign_tag"),
        "data_root": str(data_path),
        "sim_pointer": str(PP_SIM_POINTER),
        "sim_pointer_sha256": sha256(PP_SIM_POINTER),
        "sim_campaign": sim_pointer.get("campaign_tag"),
        "sim_root": str(sim_path),
        "histograms": {
            "region_A": f"{pp.DATA_TOP}/h2_unfoldReco_pTgamma_xJ_incl_{base}",
            "region_C": f"{pp.DATA_TOP}/h2_unfoldReco_pTgamma_xJ_incl_sidebandC_{base}",
        },
    }
    return result, provenance


def build_auau() -> tuple[dict[str, Any], dict[str, Any]]:
    audit = load_module("the89_auau_current_helper", AUAU_AUDIT_PATH)
    helper = audit.load_helper()
    helper.PT_WINDOW = AUAU_PT
    audit.PT_WINDOW = AUAU_PT
    data_path, data_pointer = current_root(AUAU_DATA_POINTER)
    sim_path, sim_pointer = current_root(AUAU_SIM_POINTER)
    case = audit.configure_case(helper, data_path, sim_path)
    data = open_root(data_path)
    sim = open_root(sim_path)
    base = audit.THRESHOLDS[JET_PT_MIN]
    names = {
        "region_A": f"h2_unfoldReco_pTgamma_xJ_incl_{base}{audit.CENT_SUFFIX}",
        "region_C": f"h2_unfoldReco_pTgamma_xJ_incl_sidebandC_{base}{audit.CENT_SUFFIX}",
        "combinatoric": f"h2_unfoldRecoCombinatoric_pTgamma_xJ_incl_{base}{audit.CENT_SUFFIX}",
    }
    try:
        h_a = helper.get_obj(data, case.data_topdir, names["region_A"], "TH2")
        h_c = helper.get_obj(data, case.data_topdir, names["region_C"], "TH2")
        h_k = helper.get_obj(sim, case.sim_topdir, names["combinatoric"], "TH2")
        purity, rows = audit.build_purity_input(helper, case, data, sim, h_a, h_c)
        after_k, k_meta = audit.subtract_combinatoric(helper, case, data, sim, purity, h_k)
        raw = project_counts(h_a, AUAU_PT)
        purity_projected = project_counts(purity, AUAU_PT)
        corrected = project_counts(after_k, AUAU_PT)
    finally:
        data.Close()
        sim.Close()

    edges = ensure_matching_edges(raw, purity_projected, corrected)
    photon_background = raw["values"] - purity_projected["values"]
    combinatoric = purity_projected["values"] - corrected["values"]
    closure = raw["values"] - photon_background - combinatoric - corrected["values"]
    selected_rows = [row for row in rows if row["pt"][0] >= AUAU_PT[0] and row["pt"][1] <= AUAU_PT[1]]
    result = {
        "key": "auau_0_20",
        "collision": "Au+Au 0-20%",
        "pt_window": list(AUAU_PT),
        "edges": edges,
        "raw": raw["values"],
        "raw_error": raw["errors"],
        "photon_background": photon_background,
        "combinatoric": combinatoric,
        "corrected": corrected["values"],
        "corrected_error": corrected["errors"],
        "additivity_residual": closure,
        "abcd_rows": selected_rows,
        "combinatoric_meta": k_meta,
    }
    provenance = {
        "data_pointer": str(AUAU_DATA_POINTER),
        "data_pointer_sha256": sha256(AUAU_DATA_POINTER),
        "data_campaign": data_pointer.get("campaign_tag"),
        "data_root": str(data_path),
        "sim_pointer": str(AUAU_SIM_POINTER),
        "sim_pointer_sha256": sha256(AUAU_SIM_POINTER),
        "sim_campaign": sim_pointer.get("campaign_tag"),
        "sim_root": str(sim_path),
        "histograms": {key: f"{case.data_topdir if key != 'combinatoric' else case.sim_topdir}/{value}" for key, value in names.items()},
    }
    return result, provenance


def finite_json(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: finite_json(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [finite_json(item) for item in value]
    if isinstance(value, (np.floating, float)):
        return float(value) if math.isfinite(float(value)) else None
    if isinstance(value, (np.integer, int)):
        return int(value)
    return value


def configure_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "Nimbus Roman", "DejaVu Serif"],
            "mathtext.fontset": "stix",
            "axes.linewidth": 1.15,
            "axes.labelsize": 13.8,
            "xtick.labelsize": 10.7,
            "ytick.labelsize": 10.7,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "legend.frameon": False,
        }
    )


def draw_sphenix_label(ax) -> None:
    ax.text(
        0.035,
        0.965,
        r"$\bf{\it{sPHENIX}}$ Internal",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=13.8,
        color="#111111",
    )


def draw_sphenix_panel(ax, result: dict[str, Any], *, auau: bool) -> None:
    edges = result["edges"]
    centers = 0.5 * (edges[:-1] + edges[1:])
    widths = np.diff(edges)
    mask = (edges[:-1] >= DISPLAY_X[0] - 1.0e-12) & (edges[1:] <= DISPLAY_X[1] + 1.0e-12)
    indices = np.flatnonzero(mask)
    if not len(indices):
        raise RuntimeError("no bins in display range")
    lo, hi = int(indices[0]), int(indices[-1]) + 1
    display_edges = edges[lo : hi + 1]
    raw = result["raw"][lo:hi]
    photon_bkg = result["photon_background"][lo:hi]
    comb = result["combinatoric"][lo:hi]
    corrected = result["corrected"][lo:hi]
    corrected_error = result["corrected_error"][lo:hi]
    display_centers = centers[lo:hi]
    display_widths = widths[lo:hi]

    ax.stairs(raw, display_edges, color="#949494", lw=2.0, label="Raw pairs")
    ax.stairs(
        photon_bkg,
        display_edges,
        color="#3255d9",
        lw=1.9,
        linestyle=(0, (2, 2)),
        label="ABCD photon bkg. (net)",
    )
    if auau:
        ax.stairs(
            comb,
            display_edges,
            color="#e33d32",
            lw=2.0,
            linestyle=(0, (1, 1)),
            label="Comb. recoil bkg.",
        )
    else:
        ax.plot(
            [DISPLAY_X[0], DISPLAY_X[1]],
            [0.0, 0.0],
            color="#e33d32",
            lw=2.0,
            linestyle=(0, (1, 1)),
            label="Comb. recoil bkg. (Au+Au only)",
        )
    ax.errorbar(
        display_centers,
        corrected,
        xerr=0.5 * display_widths,
        yerr=corrected_error,
        fmt="o",
        color="#111111",
        mfc="#111111",
        mec="#111111",
        ms=4.2,
        elinewidth=1.0,
        capsize=0,
        label="Corrected reco input",
        zorder=5,
    )

    candidates = np.concatenate(
        [raw, photon_bkg, comb, corrected - corrected_error, corrected + corrected_error]
    )
    finite = candidates[np.isfinite(candidates)]
    ymin = min(0.0, float(np.min(finite)))
    ymax = max(1.0, float(np.max(finite)))
    span = ymax - ymin
    ax.set_ylim(ymin - 0.07 * span, ymax + 0.52 * span)
    ax.set_xlim(*DISPLAY_X)
    ax.set_xticks(np.arange(0.2, 1.81, 0.2))
    ax.set_xlabel(r"Reconstructed $x_{J\gamma}$", labelpad=5)
    ax.set_ylabel("Photon-jet pairs / bin", labelpad=5)
    ax.minorticks_on()
    ax.tick_params(which="major", length=5.5, width=1.05)
    ax.tick_params(which="minor", length=3.0, width=0.9)
    ax.set_box_aspect(1.0)
    ax.set_anchor("C")
    draw_sphenix_label(ax)
    energy = r"Au+Au $\sqrt{s_{NN}}=200$ GeV" if auau else r"p+p $\sqrt{s}=200$ GeV"
    ax.text(0.035, 0.885, energy, transform=ax.transAxes, ha="left", va="top", fontsize=11.2)


def panel_bbox(left: float, bottom: float, width: float, height: float) -> list[float]:
    return [
        round(2560 * left, 1),
        round(1440 * (1.0 - bottom - height), 1),
        round(2560 * (left + width), 1),
        round(1440 * (1.0 - bottom), 1),
    ]


def render(pp: dict[str, Any], auau: dict[str, Any]) -> None:
    configure_style()
    fig = plt.figure(figsize=(16, 9), dpi=160, facecolor="white")
    title_color = "#111111"
    body_color = "#30343b"
    fig.text(
        0.050,
        0.958,
        "Our reconstructed input uses the same high-level subtraction logic as ATLAS",
        ha="left",
        va="top",
        fontsize=27.5,
        fontweight="bold",
        color=title_color,
    )
    fig.text(
        0.052,
        0.902,
        "Published ATLAS reference above; current sPHENIX reconstructed inputs below",
        ha="left",
        va="top",
        fontsize=16.2,
        color=body_color,
    )

    # All four visual pads use the same physical square.  Wide logical slots
    # preserve column alignment while Matplotlib centers each square pad.
    panel_w = 0.395
    panel_h = 0.340
    # Place the four equal squares on quarter points.  This keeps the two
    # columns exactly symmetric while opening a calmer center gutter for the
    # shared selection contract.
    left_x, right_x = 0.0525, 0.5525
    top_y, bottom_y = 0.485, 0.060
    positions = {
        "atlas_pp": [left_x, top_y, panel_w, panel_h],
        "atlas_auau": [right_x, top_y, panel_w, panel_h],
        "sphenix_pp": [left_x, bottom_y, panel_w, panel_h],
        "sphenix_auau": [right_x, bottom_y, panel_w, panel_h],
    }
    fig.text(
        0.5,
        0.850,
        "ATLAS published Figure 1",
        ha="center",
        va="bottom",
        fontsize=18.2,
        fontweight="bold",
        color=title_color,
    )
    fig.text(
        0.5,
        0.455,
        "sPHENIX current reconstructed inputs",
        ha="center",
        va="bottom",
        fontsize=18.2,
        fontweight="bold",
        color=title_color,
    )

    ax_atlas_pp = fig.add_axes(positions["atlas_pp"])
    ax_atlas_auau = fig.add_axes(positions["atlas_auau"])
    for ax, image in ((ax_atlas_pp, ATLAS_PP), (ax_atlas_auau, ATLAS_AUAU)):
        ax.set_box_aspect(1.0)
        ax.set_anchor("C")
        ax.imshow(mpimg.imread(image), interpolation="lanczos", aspect="equal")
        ax.axis("off")

    # Published Figure 1 selection, verified against the ATLAS paper rather
    # than inferred only from the raster labels.  The hierarchy mirrors the
    # sPHENIX selection block below so the threshold and acceptance differences
    # are immediately comparable.
    fig.text(
        0.5,
        0.785,
        "ATLAS Figure 1 selection",
        ha="center",
        va="center",
        fontsize=14.2,
        fontweight="bold",
        color=title_color,
    )
    for y, label in (
        (0.750, r"$p_T^\gamma=63.1$--$79.6$ GeV"),
        (0.715, r"anti-$k_T$ $R=0.4$  |  $p_T^{jet}>31.6$ GeV"),
        (0.680, r"$|\eta^\gamma|<2.37$ (excl. 1.37--1.52)  |  $|\eta^{jet}|<2.8$"),
        (0.645, r"$|\Delta\phi_{\gamma j}|>7\pi/8$"),
    ):
        fig.text(
            0.5,
            y,
            label,
            ha="center",
            va="center",
            fontsize=10.8,
            color=body_color,
        )
    atlas_system_x = {"pp": 0.438, "pbpb": 0.538}
    fig.text(
        atlas_system_x["pp"],
        0.595,
        r"$p+p$",
        ha="center",
        va="center",
        fontsize=11.6,
        fontweight="bold",
        color=title_color,
    )
    fig.text(
        atlas_system_x["pbpb"],
        0.595,
        r"Pb+Pb 0--10%",
        ha="center",
        va="center",
        fontsize=11.6,
        fontweight="bold",
        color=title_color,
    )
    for x, energy, luminosity, isolation in (
        (atlas_system_x["pp"], r"$\sqrt{s}=5.02$ TeV", r"25 pb$^{-1}$", r"$E_T^{iso}<3$ GeV"),
        (atlas_system_x["pbpb"], r"$\sqrt{s_{NN}}=5.02$ TeV", r"0.49 nb$^{-1}$", r"$E_T^{iso}<8$ GeV"),
    ):
        fig.text(
            x,
            0.565,
            energy,
            ha="center",
            va="center",
            fontsize=10.3,
            color=body_color,
        )
        fig.text(
            x,
            0.535,
            luminosity,
            ha="center",
            va="center",
            fontsize=10.3,
            color=body_color,
        )
        fig.text(
            x,
            0.505,
            isolation,
            ha="center",
            va="center",
            fontsize=10.3,
            color=body_color,
        )
    fig.add_artist(
        Line2D(
            [0.488, 0.488],
            [0.493, 0.607],
            transform=fig.transFigure,
            color="#d5d7da",
            linewidth=0.9,
        )
    )

    draw_sphenix_panel(fig.add_axes(positions["sphenix_pp"]), pp, auau=False)
    draw_sphenix_panel(fig.add_axes(positions["sphenix_auau"]), auau, auau=True)

    shared_handles = [
        Line2D([0], [0], color="#949494", lw=2.4, label="Raw pairs"),
        Line2D([0], [0], color="#3255d9", lw=2.2, linestyle=(0, (2, 2)), label="ABCD photon background (net)"),
        Line2D([0], [0], color="#e33d32", lw=2.2, linestyle=(0, (1, 1)), label="Combinatoric recoil (Au+Au only)"),
        Line2D([0], [0], color="#111111", marker="o", linestyle="none", markersize=5.2, label="Corrected reco input"),
    ]
    fig.legend(
        handles=shared_handles,
        loc="center",
        bbox_to_anchor=(0.5, 0.427),
        ncol=4,
        frameon=False,
        fontsize=11.2,
        handlelength=2.1,
        handletextpad=0.55,
        columnspacing=1.35,
        borderaxespad=0.0,
    )
    fig.text(
        0.5,
        0.325,
        "Current selection",
        ha="center",
        va="center",
        fontsize=14.2,
        fontweight="bold",
        color=title_color,
    )
    fig.text(
        0.5,
        0.290,
        r"anti-$k_T$ $R=0.4$  |  $p_T^{jet}>7$ GeV",
        ha="center",
        va="center",
        fontsize=11.5,
        color=body_color,
    )
    fig.text(
        0.5,
        0.257,
        r"$|\eta^{\gamma,jet}|<0.7$  |  $|\Delta\phi|>7\pi/8$",
        ha="center",
        va="center",
        fontsize=11.5,
        color=body_color,
    )
    # The system-specific cuts are intentionally split into two compact
    # columns.  Long single-line cut strings made the lower row feel compressed
    # and competed with the Au+Au y-axis title.
    system_x = {"pp": 0.438, "auau": 0.538}
    fig.text(
        system_x["pp"],
        0.210,
        r"$p+p$",
        ha="center",
        va="center",
        fontsize=11.8,
        fontweight="bold",
        color=title_color,
    )
    fig.text(
        system_x["auau"],
        0.210,
        r"Au+Au 0--20%",
        ha="center",
        va="center",
        fontsize=11.8,
        fontweight="bold",
        color=title_color,
    )
    for x, photon_cut, vertex_cut in (
        (system_x["pp"], r"$20<p_T^\gamma<26$ GeV", r"$|z_{vtx}|<60$ cm"),
        (system_x["auau"], r"$19<p_T^\gamma<26$ GeV", r"$|z_{vtx}|<10$ cm"),
    ):
        fig.text(
            x,
            0.180,
            photon_cut,
            ha="center",
            va="center",
            fontsize=10.8,
            color=body_color,
        )
        fig.text(
            x,
            0.151,
            vertex_cut,
            ha="center",
            va="center",
            fontsize=10.8,
            color=body_color,
        )
    fig.add_artist(
        Line2D(
            [0.488, 0.488],
            [0.142, 0.220],
            transform=fig.transFigure,
            color="#d5d7da",
            linewidth=0.9,
        )
    )
    fig.savefig(PNG, dpi=160, facecolor="white")
    plt.close(fig)

    nodes = [
        {"type": "text", "name": "title", "bbox": [128, 52, 2420, 126], "font_px": 61},
        {"type": "text", "name": "logic_line", "bbox": [133, 132, 2390, 185], "font_px": 36},
    ]
    square_w = panel_h * 9.0 / 16.0
    square_positions = {
        key: [left + 0.5 * (width - square_w), bottom, square_w, height]
        for key, (left, bottom, width, height) in positions.items()
    }
    nodes.extend(
        [
            {"type": "text", "name": "atlas_row_header", "bbox": [760, 205, 1800, 245], "font_px": 40},
            {"type": "text", "name": "atlas_selection_header", "bbox": [930, 285, 1630, 330], "font_px": 32},
            {"type": "text", "name": "atlas_selection_photon_pt", "bbox": [995, 340, 1565, 380], "font_px": 24},
            {"type": "text", "name": "atlas_selection_jet", "bbox": [930, 390, 1630, 430], "font_px": 24},
            {"type": "text", "name": "atlas_selection_eta", "bbox": [895, 440, 1665, 485], "font_px": 24},
            {"type": "text", "name": "atlas_selection_dphi", "bbox": [1030, 490, 1530, 530], "font_px": 24},
            {"type": "text", "name": "atlas_pp_header", "bbox": [1040, 560, 1205, 595], "font_px": 26},
            {"type": "text", "name": "atlas_pp_energy", "bbox": [1000, 605, 1245, 640], "font_px": 23},
            {"type": "text", "name": "atlas_pp_luminosity", "bbox": [1015, 648, 1230, 683], "font_px": 23},
            {"type": "text", "name": "atlas_pp_iso", "bbox": [985, 690, 1260, 725], "font_px": 23},
            {"type": "text", "name": "atlas_pbpb_header", "bbox": [1270, 560, 1490, 595], "font_px": 26},
            {"type": "text", "name": "atlas_pbpb_energy", "bbox": [1250, 605, 1515, 640], "font_px": 23},
            {"type": "text", "name": "atlas_pbpb_luminosity", "bbox": [1270, 648, 1495, 683], "font_px": 23},
            {"type": "text", "name": "atlas_pbpb_iso", "bbox": [1240, 690, 1515, 725], "font_px": 23},
            {"type": "text", "name": "sphenix_row_header", "bbox": [690, 765, 1870, 810], "font_px": 40},
            {"type": "text", "name": "shared_sphenix_legend", "bbox": [300, 810, 2260, 855], "font_px": 25},
            {"type": "text", "name": "selection_header", "bbox": [980, 945, 1580, 990], "font_px": 32},
            {"type": "text", "name": "selection_common_1", "bbox": [950, 995, 1610, 1035], "font_px": 26},
            {"type": "text", "name": "selection_common_2", "bbox": [950, 1040, 1610, 1080], "font_px": 26},
            {"type": "text", "name": "selection_pp_header", "bbox": [1040, 1110, 1205, 1145], "font_px": 26},
            {"type": "text", "name": "selection_pp_pt", "bbox": [955, 1150, 1290, 1185], "font_px": 24},
            {"type": "text", "name": "selection_pp_vertex", "bbox": [980, 1195, 1265, 1230], "font_px": 24},
            {"type": "text", "name": "selection_auau_header", "bbox": [1270, 1110, 1490, 1145], "font_px": 26},
            {"type": "text", "name": "selection_auau_pt", "bbox": [1210, 1150, 1545, 1185], "font_px": 24},
            {"type": "text", "name": "selection_auau_vertex", "bbox": [1235, 1195, 1520, 1230], "font_px": 24},
        ]
    )
    for key, position in square_positions.items():
        nodes.append({"type": "panel", "name": key, "bbox": panel_bbox(*position)})
        panel_box = panel_bbox(*position)
    LAYOUT.write_text(json.dumps(nodes, indent=2) + "\n")


def validate_contract() -> dict[str, Any]:
    for path in (ATLAS_PP, ATLAS_AUAU, ATLAS_PDF, LINKED_MANIFEST, IMMUTABLE_AUDIT):
        if not path.is_file():
            raise FileNotFoundError(path)
    linked = json.loads(LINKED_MANIFEST.read_text())
    selection = linked.get("selection", {})
    required = {
        "jet_radius": 0.4,
        "jet_pt_min_gev": 7,
        "pp_photon_pt_gev": [20.0, 26.0],
        "auau_photon_pt_gev": [19.0, 26.0],
        "pp_vertex_abs_max_cm": 60,
        "auau_vertex_abs_max_cm": 10,
        "auau_centrality_percent": [0, 20],
    }
    mismatches = {key: {"expected": value, "observed": selection.get(key)} for key, value in required.items() if selection.get(key) != value}
    if mismatches:
        raise RuntimeError(f"linked selection contract mismatch: {mismatches}")
    observed_audit_sha = sha256(IMMUTABLE_AUDIT)
    if observed_audit_sha != EXPECTED_IMMUTABLE_AUDIT_SHA:
        raise RuntimeError(
            f"immutable audit changed: expected {EXPECTED_IMMUTABLE_AUDIT_SHA}, observed {observed_audit_sha}"
        )
    if linked.get("sources", {}).get("immutable_auau_audit_sha256") != observed_audit_sha:
        raise RuntimeError("linked result and immutable audit SHA do not agree")
    return {
        "linked_manifest": str(LINKED_MANIFEST),
        "linked_manifest_sha256": sha256(LINKED_MANIFEST),
        "linked_status": linked.get("status"),
        "immutable_audit": str(IMMUTABLE_AUDIT),
        "immutable_audit_sha256": observed_audit_sha,
        "selection": required,
        "limitations": linked.get("limitations", []),
    }


def write_notes() -> None:
    NOTES.write_text(
        """# Speaker script — reconstructed subtraction chain

The top row is the published ATLAS Figure 1 reference: raw photon–jet pairs, estimated backgrounds, and the corrected reconstructed distribution before unfolding. The center block makes its kinematic contract explicit: 63.1–79.6 GeV photons, anti-kT R=0.4 jets above 31.6 GeV, |eta_gamma|<2.37 excluding the barrel–endcap transition, |eta_jet|<2.8, and |Delta phi|>7pi/8. The p+p and 0–10% Pb+Pb columns also state the collision energy, integrated luminosity, and system-specific photon-isolation threshold.

The bottom row is the equivalent *construction* from our current registered inputs. In p+p, the ABCD photon-ID correction removes the net sideband-estimated photon background. In Au+Au, we then subtract the photon-yield-scaled embedded combinatoric recoil estimate. The black points are the nominal corrected reconstructed data inputs produced by the unfolding pipeline.

The key claim is deliberately narrow: the backend follows the same high-level sequence as ATLAS. This is not a numerical ATLAS comparison—the collision energy, photon-pT windows, centrality, detector, and background estimators differ. The linked simplified Au+Au overlay uses the separately labeled response-consistent stress-test variant; it must not be equated with these nominal data-side black points. This is also not a final physics plot: the Au+Au response/combinatoric closure remains incomplete, only statistical uncertainties are shown, and the exact JES payload provenance was not re-audited for this slide.

[Sources]
- ATLAS Collaboration, Phys. Lett. B 789 (2019) 167–190, https://doi.org/10.1016/j.physletb.2018.12.023, especially Sections 1, 4, 5 and Figure 1.
- Local preserved paper: usefulDocs/external_references/gamma_jet/ATLAS_xJgamma.pdf.
- Exact ATLAS Figure 1 panel crops: dataOutput/the85_auau_xjgamma_unfolding_push/atlas_reference/.
"""
    )


def main() -> int:
    contract = validate_contract()
    pp, pp_sources = build_pp()
    auau, auau_sources = build_auau()
    if not np.array_equal(pp["edges"], auau["edges"]):
        raise RuntimeError("p+p and Au+Au xJ axes differ")
    render(pp, auau)
    write_notes()

    points = finite_json({"pp": pp, "auau_0_20": auau})
    POINTS.write_text(json.dumps(points, indent=2) + "\n")
    from PIL import Image

    with Image.open(PNG) as image:
        dimensions = list(image.size)
    additivity = {
        "pp_max_abs_bin_residual": float(np.max(np.abs(pp["additivity_residual"]))),
        "auau_max_abs_bin_residual": float(np.max(np.abs(auau["additivity_residual"]))),
    }
    manifest = {
        "status": "slide_ready_reconstructed_input_method_comparison_not_final_physics_result",
        "claim": "high_level_subtraction_chain_matches_atlas_roles_before_unfolding",
        "not_claimed": [
            "numerical agreement with ATLAS",
            "final Au+Au measurement or quenching significance",
            "closed systematic uncertainty evaluation",
            "new JES calibration validation",
        ],
        "slide_png": str(PNG),
        "slide_png_sha256": sha256(PNG),
        "dimensions_px": dimensions,
        "font_policy": "Times New Roman family; mathematical text uses STIX",
        "uncertainties": "statistical propagation only on corrected reconstructed input",
        "unfolding_on_this_slide": False,
        "unfolding_relationship": (
            "the black points are the nominal data-side corrected reconstructed inputs built by the unfolding pipeline; "
            "the linked simplified Au+Au overlay instead displays the separately labeled response-consistent stress-test variant"
        ),
        "construction": {
            "pp": "raw Region A minus net leakage-corrected ABCD photon background",
            "auau": "raw Region A minus net leakage-corrected ABCD photon background minus photon-yield-scaled embedded combinatoric recoil background",
            "component_identity": "background components are defined by exact binwise differences so raw = backgrounds + corrected input",
            "additivity_checks": additivity,
        },
        "selection": contract["selection"],
        "atlas_figure1_selection": {
            "photon_pt_gev": [63.1, 79.6],
            "jet_algorithm": "anti-kT R=0.4",
            "jet_pt_min_gev": 31.6,
            "photon_abs_eta_max": 2.37,
            "photon_eta_excluded": [1.37, 1.52],
            "jet_abs_eta_max": 2.8,
            "delta_phi_min": "7pi/8",
            "pp": {
                "sqrt_s_tev": 5.02,
                "integrated_luminosity": "25 pb^-1",
                "photon_isolation_max_gev": 3.0,
            },
            "pbpb_0_10": {
                "sqrt_snn_tev": 5.02,
                "integrated_luminosity": "0.49 nb^-1",
                "photon_isolation_max_gev": 8.0,
            },
            "source": ATLAS_DOI,
        },
        "sources": {
            "contract": contract,
            "atlas_pdf": str(ATLAS_PDF),
            "atlas_pdf_sha256": sha256(ATLAS_PDF),
            "atlas_doi": ATLAS_DOI,
            "atlas_pp_crop": str(ATLAS_PP),
            "atlas_pp_crop_sha256": sha256(ATLAS_PP),
            "atlas_pbpb_0_10_crop": str(ATLAS_AUAU),
            "atlas_pbpb_0_10_crop_sha256": sha256(ATLAS_AUAU),
            "pp": pp_sources,
            "auau": auau_sources,
            "generator": str(Path(__file__).resolve()),
            "generator_sha256": sha256(Path(__file__).resolve()),
        },
        "sidecars": {
            "points": str(POINTS),
            "layout_nodes": str(LAYOUT),
            "speaker_script": str(NOTES),
        },
        "caveats": [
            "ATLAS and sPHENIX panels have different collision energies, photon-pT windows, centralities, detectors, and background estimators.",
            "The linked Au+Au result is an offline response-contract diagnostic with incomplete response/combinatoric closure and must not be relabeled final.",
            "The linked simplified Au+Au red curve is the response-consistent stress-test variant and is not identical to the nominal data-side black points shown here.",
            "The p+p and Au+Au photon-pT windows are the nearest exact native windows and are not identical.",
            "Only statistical propagation is drawn; no systematic uncertainty band is present.",
            "Exact JES CDB payload provenance was not newly re-established for this slide; current registered ROOT pointers were used unchanged.",
        ],
    }
    MANIFEST.write_text(json.dumps(finite_json(manifest), indent=2) + "\n")
    print(PNG)
    print(MANIFEST)
    print(POINTS)
    print(NOTES)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
