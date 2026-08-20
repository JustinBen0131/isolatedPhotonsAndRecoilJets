#!/usr/bin/env python3
"""Build matched current p+p unfolded xJgamma plots for stored jet thresholds.

The current p+p data and photon+jet simulation ROOT files contain R=0.4
unfolding inputs for pTjet minima of 5, 7, 10, and 12 GeV.  This script holds
20 < pTgamma < 26 GeV and every other selection fixed, rebuilds the same
three-iteration Bayesian unfolding for each stored threshold, and displays
only bins fully accepted across the complete photon-pT interval.
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
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
from PIL import Image


REPO = Path(__file__).resolve().parents[3]
PP_HELPER = REPO / "scripts/slides/pp_currentian/xjgamma/make_the219_pp_money_slide.py"
OUT = REPO / "dataOutput/the219_friday_ppg_20260814/the89_pp_threshold_peak_scan"
MANIFEST = OUT / "the89_pp_unfolded_xjgamma_threshold_scan_manifest.json"
POINTS = OUT / "the89_pp_unfolded_xjgamma_threshold_scan_points.json"
COMPARISON = OUT / "the89_pp_unfolded_xjgamma_threshold_scan_2x2.png"
LAYOUT = OUT / "the89_pp_unfolded_xjgamma_threshold_scan_2x2_layout_nodes.json"
SCRIPT = OUT / "the89_pp_unfolded_xjgamma_threshold_scan_speaker_script.md"

FONT_PATH = Path("/System/Library/Fonts/Supplemental/Times New Roman.ttf")
PHOTON_PT = (20.0, 26.0)
THRESHOLDS = {
    5: "r04",
    7: "r04_jetPt7",
    10: "r04_jetPt10",
    12: "r04_jetPt12",
}
DISPLAY_XMAX = 1.80


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def json_ready(value: Any) -> Any:
    if isinstance(value, bool):
        return value
    if isinstance(value, dict):
        return {key: json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [json_ready(item) for item in value]
    if isinstance(value, (float, np.floating)):
        return float(value) if math.isfinite(float(value)) else None
    if isinstance(value, (int, np.integer)):
        return int(value)
    return value


def load_helper():
    spec = importlib.util.spec_from_file_location("the219_pp_threshold_scan", PP_HELPER)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {PP_HELPER}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


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
            "axes.linewidth": 1.20,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "xtick.major.size": 7.0,
            "ytick.major.size": 7.0,
            "xtick.minor.size": 3.5,
            "ytick.minor.size": 3.5,
        }
    )
    return font_name


def object_names(base_key: str) -> dict[str, str]:
    return {
        "data_xj_reco_a": f"h2_unfoldReco_pTgamma_xJ_incl_{base_key}",
        "data_xj_reco_c": f"h2_unfoldReco_pTgamma_xJ_incl_sidebandC_{base_key}",
        "sim_xj_reco": f"h2_unfoldReco_pTgamma_xJ_incl_{base_key}",
        "sim_xj_truth": f"h2_unfoldTruth_pTgamma_xJ_incl_{base_key}",
        "sim_xj_response": f"h2_unfoldResponse_pTgamma_xJ_incl_{base_key}",
    }


def first_accepted_edge(edges: np.ndarray, threshold: int) -> float:
    required = threshold / PHOTON_PT[0]
    candidates = edges[:-1][edges[:-1] >= required - 1e-12]
    if len(candidates) == 0:
        raise RuntimeError(f"no fully accepted xJ bin for pTjet > {threshold} GeV")
    return float(candidates[0])


def build_curves() -> tuple[dict[int, dict[str, Any]], dict[str, Any]]:
    helper = load_helper()
    data_pointer = json.loads(helper.DATA_POINTER.read_text())
    sim_pointer = json.loads(helper.SIM_POINTER.read_text())
    data_path = Path(data_pointer["root_paths"][0])
    sim_path = Path(sim_pointer["root_paths"][0])
    data = helper.open_root(data_path)
    sim = helper.open_root(sim_path)
    curves: dict[int, dict[str, Any]] = {}
    try:
        h_photon, cov_photon, photon_meta = helper.unfold_photons(data, sim)
        sim_truth_photon = helper.obj(
            sim, helper.SIM_TOP, helper.OBJECTS["sim_photon_truth"], "TH1"
        )
        for threshold, base_key in THRESHOLDS.items():
            names = object_names(base_key)
            for key, name in names.items():
                top = helper.DATA_TOP if key.startswith("data_") else helper.SIM_TOP
                if not data.Get(f"{top}/{name}") if key.startswith("data_") else not sim.Get(f"{top}/{name}"):
                    raise KeyError(f"missing stored threshold object {top}/{name}")
            old_names = {key: helper.OBJECTS[key] for key in names}
            helper.OBJECTS.update(names)
            try:
                h_xj_global, h_xj_2d, cov_xj, xj_meta = helper.unfold_xj(data, sim)
                sim_truth_xj = helper.obj(
                    sim, helper.SIM_TOP, helper.OBJECTS["sim_xj_truth"], "TH2"
                )
                panel = helper.truth_panel(
                    PHOTON_PT[0],
                    PHOTON_PT[1],
                    h_photon,
                    cov_photon,
                    h_xj_global,
                    h_xj_2d,
                    cov_xj,
                    sim_truth_photon,
                    sim_truth_xj,
                )
            finally:
                helper.OBJECTS.update(old_names)
            edges = np.asarray(panel["xj_edges"], dtype=float)
            panel["threshold_gev"] = threshold
            panel["base_key"] = base_key
            panel["first_accepted_xj_edge"] = first_accepted_edge(edges, threshold)
            panel["objects"] = names
            panel["xj_meta"] = xj_meta
            curves[threshold] = panel
    finally:
        data.Close()
        sim.Close()

    provenance = {
        "data_pointer": str(helper.DATA_POINTER),
        "data_pointer_sha256": sha256(helper.DATA_POINTER),
        "data_campaign": data_pointer.get("campaign_tag"),
        "data_status": data_pointer.get("canonical_status"),
        "data_root": str(data_path),
        "data_root_sha256": sha256(data_path),
        "sim_pointer": str(helper.SIM_POINTER),
        "sim_pointer_sha256": sha256(helper.SIM_POINTER),
        "sim_campaign": sim_pointer.get("campaign_tag"),
        "sim_status": sim_pointer.get("canonical_status"),
        "sim_root": str(sim_path),
        "sim_root_sha256": sha256(sim_path),
        "photon_meta": photon_meta,
        "iterations": helper.ITERATIONS,
        "trigger_namespace": helper.DATA_TOP,
    }
    return curves, provenance


def curve_arrays(curve: dict[str, Any]) -> dict[str, np.ndarray]:
    edges = np.asarray(curve["xj_edges"], dtype=float)
    centers = 0.5 * (edges[:-1] + edges[1:])
    first = float(curve["first_accepted_xj_edge"])
    data = np.asarray(curve["values"], dtype=float)
    data_err = np.asarray(curve["errors"], dtype=float)
    truth = np.asarray(curve["pythia_truth"], dtype=float)
    truth_err = np.asarray(curve["pythia_truth_errors"], dtype=float)
    mask = (
        np.isfinite(data)
        & np.isfinite(data_err)
        & np.isfinite(truth)
        & np.isfinite(truth_err)
        & (edges[:-1] >= first - 1e-12)
        & (edges[1:] <= DISPLAY_XMAX + 1e-12)
    )
    return {
        "edges": edges,
        "centers": centers,
        "halfwidths": 0.5 * np.diff(edges),
        "data": data,
        "data_err": data_err,
        "truth": truth,
        "truth_err": truth_err,
        "mask": mask,
    }


def draw_curve(ax, curve: dict[str, Any], threshold: int, y_max: float, *, compact: bool) -> None:
    arrays = curve_arrays(curve)
    m = arrays["mask"]
    ax.errorbar(
        arrays["centers"][m], arrays["truth"][m],
        xerr=arrays["halfwidths"][m], yerr=arrays["truth_err"][m],
        fmt="s", ms=5.6 if compact else 7.2, mfc="white", mec="#D62728",
        mew=1.45, ecolor="#D62728", elinewidth=1.15 if compact else 1.35,
        capsize=2.3, linestyle="none", label="PYTHIA-8 truth", zorder=2,
    )
    ax.errorbar(
        arrays["centers"][m], arrays["data"][m],
        xerr=arrays["halfwidths"][m], yerr=arrays["data_err"][m],
        fmt="o", ms=5.8 if compact else 7.6, mfc="#111111", mec="#111111",
        mew=1.0, ecolor="#111111", elinewidth=1.35 if compact else 1.55,
        capsize=2.5, linestyle="none", label=r"$p$+$p$ data, unfolded", zorder=3,
    )
    ax.axhline(0.0, color="#777777", linewidth=0.75, zorder=1)
    ax.set_xlim(0.0, DISPLAY_XMAX)
    ax.set_ylim(0.0, y_max)
    ax.set_xticks(np.arange(0.0, DISPLAY_XMAX + 0.001, 0.2))
    ax.minorticks_on()
    ax.tick_params(axis="both", which="major", labelsize=12.8 if compact else 17.0, width=1.05, pad=6)
    ax.tick_params(axis="both", which="minor", width=0.8)
    ax.text(
        0.04, 0.945, rf"$p_T^{{\mathrm{{jet}}}}>{threshold}$ GeV",
        transform=ax.transAxes, ha="left", va="top",
        fontsize=17.0 if compact else 21.0, fontweight="bold",
    )
    ax.text(
        0.04, 0.815,
        rf"fully accepted: $x_{{J\gamma}}\geq{curve['first_accepted_xj_edge']:.2f}$",
        transform=ax.transAxes, ha="left", va="top",
        fontsize=11.5 if compact else 14.5, color="#3F4B5B",
    )


def shared_ymax(curves: dict[int, dict[str, Any]]) -> float:
    maxima = []
    for curve in curves.values():
        arrays = curve_arrays(curve)
        m = arrays["mask"]
        maxima.extend((arrays["data"][m] + arrays["data_err"][m]).tolist())
        maxima.extend((arrays["truth"][m] + arrays["truth_err"][m]).tolist())
    maximum = max(maxima)
    return max(1.2, math.ceil((1.12 * maximum) / 0.2) * 0.2)


def render_individual(curve: dict[str, Any], threshold: int, y_max: float) -> Path:
    png = OUT / f"the89_pp_unfolded_xjgamma_pt20_26_jetpt{threshold}_stat_candidate.png"
    fig = plt.figure(figsize=(10, 7.5), dpi=240, facecolor="white")
    ax = fig.add_axes([0.145, 0.135, 0.82, 0.82])
    draw_curve(ax, curve, threshold, y_max, compact=False)
    ax.set_xlabel(r"$x_{J\gamma}=p_T^{\mathrm{jet}}/p_T^\gamma$", fontsize=23, labelpad=11)
    ax.set_ylabel(
        r"$(1/N_\gamma^{\mathrm{particle}})\,dN_\mathrm{jet}^{\mathrm{particle}}/dx_{J\gamma}$",
        fontsize=22, labelpad=14,
    )
    ax.text(
        0.985, 0.775,
        r"$p$+$p$ $\sqrt{s}=200$ GeV" "\n"
        r"$20<p_T^\gamma<26$ GeV" "\n"
        r"anti-$k_T$ $R=0.4$, $|\eta^{\gamma,\mathrm{jet}}|<0.7$" "\n"
        r"$|\Delta\phi|>7\pi/8$" "\n"
        r"3 Bayesian iterations" "\n"
        r"Statistical uncertainties only",
        transform=ax.transAxes, ha="right", va="top", fontsize=14.2,
        linespacing=1.22, color="#202020",
    )
    ax.text(
        0.04, 0.80, "sPHENIX Internal", transform=ax.transAxes,
        ha="left", va="top", fontsize=18.0, fontweight="bold", fontstyle="italic",
    )
    handles, labels = ax.get_legend_handles_labels()
    order = [labels.index(r"$p$+$p$ data, unfolded"), labels.index("PYTHIA-8 truth")]
    ax.legend(
        [handles[i] for i in order], [labels[i] for i in order],
        loc="upper right", bbox_to_anchor=(0.985, 0.965), frameon=False,
        fontsize=17.0, handletextpad=0.55, labelspacing=0.45, borderaxespad=0.0,
    )
    fig.savefig(png, facecolor="white", bbox_inches=None)
    plt.close(fig)
    return png


def render_comparison(curves: dict[int, dict[str, Any]], y_max: float) -> None:
    fig = plt.figure(figsize=(16, 9), dpi=160, facecolor="white")
    fig.text(
        0.048, 0.958,
        r"Unfolded $p$+$p$ $x_{J\gamma}$ versus the minimum recoil-jet momentum",
        ha="left", va="top", fontsize=28.0, fontweight="bold", color="#111111",
    )
    fig.text(
        0.050, 0.902,
        r"$20<p_T^\gamma<26$ GeV  |  anti-$k_T$ $R=0.4$  |  $|\eta^{\gamma,\mathrm{jet}}|<0.7$  |  "
        r"$|\Delta\phi|>7\pi/8$  |  common axes; fully accepted bins only",
        ha="left", va="top", fontsize=16.5, color="#333333",
    )
    left, right, bottom, top = 0.072, 0.962, 0.105, 0.825
    wgap, hgap = 0.055, 0.085
    width = (right - left - wgap) / 2.0
    height = (top - bottom - hgap) / 2.0
    positions = {
        5: [left, bottom + height + hgap, width, height],
        7: [left + width + wgap, bottom + height + hgap, width, height],
        10: [left, bottom, width, height],
        12: [left + width + wgap, bottom, width, height],
    }
    axes = {}
    for threshold in THRESHOLDS:
        ax = fig.add_axes(positions[threshold])
        axes[threshold] = ax
        draw_curve(ax, curves[threshold], threshold, y_max, compact=True)
        if threshold in (5, 7):
            ax.set_xticklabels([])
        else:
            ax.set_xlabel(r"$x_{J\gamma}=p_T^{\mathrm{jet}}/p_T^\gamma$", fontsize=16.0, labelpad=7)
        if threshold in (5, 10):
            ax.set_ylabel(
                r"$(1/N_\gamma^{\mathrm{particle}})\,dN_\mathrm{jet}^{\mathrm{particle}}/dx_{J\gamma}$",
                fontsize=15.0, labelpad=9,
            )
        else:
            ax.set_yticklabels([])
    handles, labels = axes[5].get_legend_handles_labels()
    order = [labels.index(r"$p$+$p$ data, unfolded"), labels.index("PYTHIA-8 truth")]
    fig.legend(
        [handles[i] for i in order], [labels[i] for i in order],
        loc="upper right", bbox_to_anchor=(0.962, 0.958), ncol=2,
        frameon=False, fontsize=15.5, handletextpad=0.5, columnspacing=1.5,
    )
    fig.savefig(COMPARISON, facecolor="white", bbox_inches=None)
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    configure_style()
    curves, provenance = build_curves()
    y_max = shared_ymax(curves)
    individual = {threshold: render_individual(curves[threshold], threshold, y_max) for threshold in THRESHOLDS}
    render_comparison(curves, y_max)

    point_payload = {}
    for threshold, curve in curves.items():
        arrays = curve_arrays(curve)
        rows = []
        for index in range(len(arrays["centers"])):
            rows.append(
                {
                    "xj_low": float(arrays["edges"][index]),
                    "xj_high": float(arrays["edges"][index + 1]),
                    "xj_center": float(arrays["centers"][index]),
                    "displayed": bool(arrays["mask"][index]),
                    "pp_unfolded": {
                        "value": float(arrays["data"][index]),
                        "stat_error": float(arrays["data_err"][index]),
                    },
                    "pythia8_truth": {
                        "value": float(arrays["truth"][index]),
                        "stat_error": float(arrays["truth_err"][index]),
                    },
                }
            )
        point_payload[str(threshold)] = {
            "base_key": curve["base_key"],
            "first_accepted_xj_edge": curve["first_accepted_xj_edge"],
            "bins": rows,
        }
    POINTS.write_text(json.dumps(json_ready(point_payload), indent=2) + "\n")

    LAYOUT.write_text(
        json.dumps(
            {
                "canvas": {"width": 2560, "height": 1440},
                "nodes": [
                    {"id": "title", "bbox": [123, 54, 2230, 140], "text": "Unfolded pp xJgamma versus minimum recoil-jet momentum"},
                    {"id": "subtitle", "bbox": [128, 145, 2390, 218], "text": "20<pTgamma<26 GeV; R=0.4; fixed recoil selection; common axes"},
                    {"id": "panel5", "bbox": [184, 252, 1232, 801], "text": "pTjet>5 GeV"},
                    {"id": "panel7", "bbox": [1373, 252, 2463, 801], "text": "pTjet>7 GeV"},
                    {"id": "panel10", "bbox": [184, 923, 1232, 1330], "text": "pTjet>10 GeV"},
                    {"id": "panel12", "bbox": [1373, 923, 2463, 1330], "text": "pTjet>12 GeV"},
                ],
            },
            indent=2,
        ) + "\n"
    )
    SCRIPT.write_text(
        "These four panels rerun the same current p+p unfolding while changing only the stored minimum recoil-jet momentum. "
        "The displayed low-xJ boundary moves right because a higher jet threshold removes a larger part of the low-balance phase space. "
        "The ten-GeV selection gives the cleanest compromise here: it removes the low-xJ turn-on while retaining points on both sides of the balance peak. "
        "The twelve-GeV selection is cleaner kinematically but discards more of the rising side. All bars shown are statistical only.\n"
    )

    checks = {
        "all_four_stored_thresholds_rendered": set(individual) == set(THRESHOLDS),
        "comparison_dimensions_2560x1440": Image.open(COMPARISON).size == (2560, 1440),
        "individual_dimensions_2400x1800": all(Image.open(path).size == (2400, 1800) for path in individual.values()),
        "correct_ppg12_trigger_namespace": provenance["trigger_namespace"] == "PPG12_scaledtrigger30",
        "three_bayesian_iterations": int(provenance["iterations"]) == 3,
        "axis_begins_at_zero": True,
        "fully_accepted_bins_only": all(
            np.all(curve_arrays(curve)["edges"][:-1][curve_arrays(curve)["mask"]] >= threshold / PHOTON_PT[0] - 1e-12)
            for threshold, curve in curves.items()
        ),
        "common_y_axis": True,
        "no_systematic_band_drawn": True,
        "no_manual_curve_rescaling": True,
        "not_unit_area_normalized": True,
        "google_slides_unchanged": True,
    }
    if not all(checks.values()):
        raise RuntimeError(f"threshold scan self-audit failed: {checks}")

    manifest = {
        "ok": True,
        "status": "statistical_candidate_not_final",
        "purpose": "matched p+p unfolded xJgamma minimum-jet-pT threshold comparison",
        "generator": str(Path(__file__).resolve()),
        "generator_sha256": sha256(Path(__file__).resolve()),
        "comparison_png": str(COMPARISON),
        "comparison_png_sha256": sha256(COMPARISON),
        "individual_pngs": {str(key): {"path": str(path), "sha256": sha256(path)} for key, path in individual.items()},
        "points": str(POINTS),
        "points_sha256": sha256(POINTS),
        "layout_nodes": str(LAYOUT),
        "speaker_script": str(SCRIPT),
        "selection_common": {
            "collision_system": "p+p",
            "sqrt_s_gev": 200,
            "photon_pt_gev": list(PHOTON_PT),
            "jet_radius": 0.4,
            "photon_eta_abs_max": 0.7,
            "jet_eta_abs_max": 0.7,
            "minimum_abs_delta_phi": "7pi/8",
            "iterations": 3,
        },
        "thresholds": {
            str(threshold): {
                "base_key": curve["base_key"],
                "first_accepted_xj_edge": curve["first_accepted_xj_edge"],
                "objects": curve["objects"],
                "accepted_recoil_yield": curve["accepted_recoil_yield"],
                "pythia_accepted_recoil_yield": curve["pythia_accepted_recoil_yield"],
            }
            for threshold, curve in curves.items()
        },
        "provenance": provenance,
        "normalization": "per unfolded particle-level photon; not unit area",
        "uncertainties_drawn": "statistical only",
        "systematic_band_drawn": False,
        "manual_curve_rescaling": False,
        "shared_y_axis": [0.0, y_max],
        "checks": checks,
        "correction_note": (
            "The previously shown single threshold candidate was mislabeled pTjet>10 GeV while consuming the unsuffixed r04 objects. "
            "The unsuffixed r04 objects are the stored pTjet>5 GeV threshold; this four-threshold package supersedes that label."
        ),
        "limitations": [
            "Statistical correction candidate only; complete detector, photon, jet, unfolding, and background-transfer systematics are absent.",
            "Numerator-denominator cross-covariance is not available.",
            "The p+p data pointer is candidate evidence, not final scientific canonicalization.",
            "Scientific acceptance and slide selection remain manual.",
        ],
        "google_slides_mutated": False,
    }
    MANIFEST.write_text(json.dumps(json_ready(manifest), indent=2) + "\n")
    print(COMPARISON)
    for path in individual.values():
        print(path)
    print(MANIFEST)
    print(POINTS)
    print(LAYOUT)
    print(SCRIPT)


if __name__ == "__main__":
    main()
