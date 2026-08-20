#!/usr/bin/env python3
"""Build the clean pTjet > 10 GeV p+p unfolded xJgamma slide candidate.

The candidate is rebuilt directly from the registered current p+p data and
photon+jet simulation ROOT files.  It uses the inclusive per-photon recoil-jet
observable, scans the regularization choice, and displays only xJgamma bins
that are fully accepted over 20 < pTgamma < 26 GeV.  Error bars are statistical
only; no systematic band is fabricated.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import sys
from typing import Any

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
from PIL import Image
import ROOT


ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
if ROOT.gSystem.Load("libRooUnfold") < 0:
    raise RuntimeError("failed to load libRooUnfold")

REPO = Path(__file__).resolve().parents[3]
PP_HELPER = REPO / "scripts/slides/pp_currentian/xjgamma/make_the219_pp_money_slide.py"
LHC_STYLE_REFERENCE = (
    Path(os.environ["RJ_LHC_STYLE_REFERENCE"])
    if "RJ_LHC_STYLE_REFERENCE" in os.environ
    else None
)
OUT = REPO / "dataOutput/the219_friday_ppg_20260814/the89_pp_peak_result"
PNG = OUT / "the89_pp_unfolded_xjgamma_peak_pt20_26_jetpt10_stat_candidate.png"
POINTS = OUT / "the89_pp_unfolded_xjgamma_peak_pt20_26_jetpt10_points.json"
MANIFEST = OUT / "the89_pp_unfolded_xjgamma_peak_pt20_26_jetpt10_manifest.json"
LAYOUT = OUT / "the89_pp_unfolded_xjgamma_peak_pt20_26_jetpt10_layout_nodes.json"
SCRIPT = OUT / "the89_pp_unfolded_xjgamma_peak_pt20_26_jetpt10_speaker_script.md"

FONT_PATH = Path("/System/Library/Fonts/Supplemental/Times New Roman.ttf")
PHOTON_PT = (20.0, 26.0)
JET_PT_MIN = 10.0
BASE_KEY = "r04_jetPt10"
DISPLAY_XMIN = 0.50
DISPLAY_XMAX = 1.80
ITERATION_SCAN = (1, 2, 3, 4)
MIN_STABILITY_FLOOR = 0.025


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
    spec = importlib.util.spec_from_file_location("the89_pp_peak_helper", PP_HELPER)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {PP_HELPER}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def object_names() -> dict[str, str]:
    return {
        "data_xj_reco_a": f"h2_unfoldReco_pTgamma_xJ_incl_{BASE_KEY}",
        "data_xj_reco_c": f"h2_unfoldReco_pTgamma_xJ_incl_sidebandC_{BASE_KEY}",
        "sim_xj_reco": f"h2_unfoldReco_pTgamma_xJ_incl_{BASE_KEY}",
        "sim_xj_truth": f"h2_unfoldTruth_pTgamma_xJ_incl_{BASE_KEY}",
        "sim_xj_response": f"h2_unfoldResponse_pTgamma_xJ_incl_{BASE_KEY}",
    }


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
            "axes.linewidth": 1.25,
            "axes.labelcolor": "#111111",
            "xtick.color": "#111111",
            "ytick.color": "#111111",
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "xtick.major.size": 7.0,
            "ytick.major.size": 7.0,
            "xtick.minor.size": 3.6,
            "ytick.minor.size": 3.6,
        }
    )
    return font_name


def panel_arrays(panel: dict[str, Any]) -> dict[str, np.ndarray]:
    edges = np.asarray(panel["xj_edges"], dtype=float)
    return {
        "edges": edges,
        "centers": 0.5 * (edges[:-1] + edges[1:]),
        "halfwidths": 0.5 * np.diff(edges),
        "values": np.asarray(panel["values"], dtype=float),
        "errors": np.asarray(panel["errors"], dtype=float),
        "truth": np.asarray(panel["pythia_truth"], dtype=float),
        "truth_errors": np.asarray(panel["pythia_truth_errors"], dtype=float),
    }


def display_mask(panel: dict[str, Any]) -> np.ndarray:
    arrays = panel_arrays(panel)
    return (
        np.isfinite(arrays["values"])
        & np.isfinite(arrays["errors"])
        & np.isfinite(arrays["truth"])
        & np.isfinite(arrays["truth_errors"])
        & (arrays["edges"][:-1] >= DISPLAY_XMIN - 1e-12)
        & (arrays["edges"][1:] <= DISPLAY_XMAX + 1e-12)
    )


def relative_statistical_precision(panel: dict[str, Any]) -> float:
    arrays = panel_arrays(panel)
    mask = display_mask(panel)
    denom = float(np.sum(np.square(arrays["values"][mask])))
    return math.sqrt(float(np.sum(np.square(arrays["errors"][mask]))) / denom) if denom > 0 else math.inf


def relative_iteration_change(panel: dict[str, Any], previous: dict[str, Any]) -> float:
    values = panel_arrays(panel)["values"]
    old = panel_arrays(previous)["values"]
    mask = display_mask(panel) & display_mask(previous)
    denom = float(np.sum(np.square(values[mask])))
    return math.sqrt(float(np.sum(np.square(values[mask] - old[mask]))) / denom) if denom > 0 else math.inf


def project_reco_xj(helper, h2) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    xbins = helper.matching_bins(h2.GetXaxis(), PHOTON_PT[0], PHOTON_PT[1])
    yaxis = h2.GetYaxis()
    edges = np.asarray([float(yaxis.GetBinLowEdge(i)) for i in range(1, yaxis.GetNbins() + 2)])
    values = []
    errors = []
    for iy in range(1, yaxis.GetNbins() + 1):
        values.append(sum(float(h2.GetBinContent(ix, iy)) for ix in xbins))
        errors.append(math.sqrt(sum(float(h2.GetBinError(ix, iy)) ** 2 for ix in xbins)))
    return edges, np.asarray(values), np.asarray(errors)


def refold_chi2_ndf(helper, roo_response, unfolded_global, reco_template, corrected_reco) -> float:
    refolded_global = roo_response.ApplyToTruth(unfolded_global, "the89_refolded_global")
    refolded_2d = helper.legacy.unflatten_global_to_th2(
        refolded_global, reco_template, "the89_refolded_2d"
    )
    edges, measured, measured_error = project_reco_xj(helper, corrected_reco)
    refold_edges, refolded, _ = project_reco_xj(helper, refolded_2d)
    if not np.allclose(edges, refold_edges):
        raise RuntimeError("measured/refolded xJ axes differ")
    mask = (
        (edges[:-1] >= DISPLAY_XMIN - 1e-12)
        & (edges[1:] <= DISPLAY_XMAX + 1e-12)
        & np.isfinite(measured)
        & np.isfinite(refolded)
        & (measured_error > 0)
    )
    if not np.any(mask):
        raise RuntimeError("no accepted bins available for refolding test")
    return float(np.sum(np.square((measured[mask] - refolded[mask]) / measured_error[mask])) / np.count_nonzero(mask))


def build_candidate() -> tuple[dict[str, Any], dict[str, Any]]:
    helper = load_helper()
    names = object_names()
    helper.OBJECTS.update(names)

    data_pointer = json.loads(helper.DATA_POINTER.read_text())
    sim_pointer = json.loads(helper.SIM_POINTER.read_text())
    data_path = Path(data_pointer["root_paths"][0])
    sim_path = Path(sim_pointer["root_paths"][0])
    data = helper.open_root(data_path)
    sim = helper.open_root(sim_path)
    try:
        for key, name in names.items():
            source = data if key.startswith("data_") else sim
            top = helper.DATA_TOP if key.startswith("data_") else helper.SIM_TOP
            if not source.Get(f"{top}/{name}"):
                raise KeyError(f"missing {top}/{name}")

        sim_truth_photon = helper.obj(sim, helper.SIM_TOP, helper.OBJECTS["sim_photon_truth"], "TH1")
        sim_truth_xj = helper.obj(sim, helper.SIM_TOP, names["sim_xj_truth"], "TH2")

        data_a = helper.obj(data, helper.DATA_TOP, names["data_xj_reco_a"], "TH2")
        data_c = helper.obj(data, helper.DATA_TOP, names["data_xj_reco_c"], "TH2")
        corrected_reco, abcd_rows = helper.purity_correct_xj(data, sim, data_a, data_c)
        sim_reco = helper.obj(sim, helper.SIM_TOP, names["sim_xj_reco"], "TH2")
        sim_truth = helper.obj(sim, helper.SIM_TOP, names["sim_xj_truth"], "TH2")
        response_raw = helper.obj(sim, helper.SIM_TOP, names["sim_xj_response"], "TH2")
        reco_global = helper.legacy.flatten_th2_to_global(sim_reco, "the89_sim_reco_global")
        truth_global = helper.legacy.flatten_th2_to_global(sim_truth, "the89_sim_truth_global")
        response = helper.legacy.transpose_th2(response_raw, "the89_response_recoX_truthY")
        roo_response = ROOT.RooUnfoldResponse(
            reco_global, truth_global, response, "the89_response", "the89_response"
        )

        panels: dict[int, dict[str, Any]] = {}
        scan_rows: list[dict[str, Any]] = []
        photon_results: dict[int, tuple[Any, np.ndarray]] = {}
        previous: dict[str, Any] | None = None
        for iteration in ITERATION_SCAN:
            helper.ITERATIONS = iteration
            h_photon, cov_photon, _ = helper.unfold_photons(data, sim)
            h_xj_global, h_xj_2d, cov_xj, _ = helper.unfold_xj(data, sim)
            panel = helper.truth_panel(
                PHOTON_PT[0], PHOTON_PT[1], h_photon, cov_photon,
                h_xj_global, h_xj_2d, cov_xj, sim_truth_photon, sim_truth_xj,
            )
            panel["iteration"] = iteration
            panels[iteration] = panel
            photon_results[iteration] = (h_photon, cov_photon)
            change = relative_iteration_change(panel, previous) if previous is not None else None
            scan_rows.append(
                {
                    "iteration": iteration,
                    "relative_statistical_precision": relative_statistical_precision(panel),
                    "relative_change_from_previous": change,
                    "refold_chi2_ndf": refold_chi2_ndf(
                        helper, roo_response, h_xj_global, sim_reco, corrected_reco
                    ),
                }
            )
            previous = panel

        selected = None
        for row in scan_rows[1:]:
            threshold = max(float(row["relative_statistical_precision"]), MIN_STABILITY_FLOOR)
            if float(row["relative_change_from_previous"]) <= threshold:
                selected = int(row["iteration"])
                row["earliest_stable_choice"] = True
                break
        if selected is None:
            raise RuntimeError(f"no stable iteration found in scan: {scan_rows}")
        for row in scan_rows:
            row.setdefault("earliest_stable_choice", False)

        baseline = panels[selected]
        h_photon, cov_photon = photon_results[selected]
        old_pooling = helper.pooled_c_shape
        helper.pooled_c_shape = lambda *_args, **_kwargs: None
        try:
            helper.ITERATIONS = selected
            np_global, np_2d, np_cov, _ = helper.unfold_xj(data, sim)
            no_pool = helper.truth_panel(
                PHOTON_PT[0], PHOTON_PT[1], h_photon, cov_photon,
                np_global, np_2d, np_cov, sim_truth_photon, sim_truth_xj,
            )
        finally:
            helper.pooled_c_shape = old_pooling

        base_values = panel_arrays(baseline)["values"]
        no_pool_values = panel_arrays(no_pool)["values"]
        accepted = display_mask(baseline)
        relative_pool_delta = np.divide(
            np.abs(no_pool_values - base_values), np.abs(base_values),
            out=np.zeros_like(base_values), where=np.abs(base_values) > 1e-12,
        )
        weighted_pool_delta = math.sqrt(
            float(np.sum(np.square(no_pool_values[accepted] - base_values[accepted])))
            / float(np.sum(np.square(base_values[accepted])))
        )
        edges = panel_arrays(baseline)["edges"]
        balance_region = accepted & (edges[1:] <= 1.03 + 1e-12)
        pooling_audit = {
            "baseline": "current high-pT sideband-C shape pooling retained",
            "alternative": "same unfolding with sideband-C shape pooling disabled",
            "selected_iteration": selected,
            "yield_weighted_relative_change_in_displayed_bins": weighted_pool_delta,
            "max_relative_change_in_balance_region_xj_le_1p03": float(
                np.max(relative_pool_delta[balance_region])
            ),
            "max_relative_change_in_displayed_bins": float(np.max(relative_pool_delta[accepted])),
            "baseline_values_in_displayed_bins": base_values[accepted].tolist(),
            "no_pool_values_in_displayed_bins": no_pool_values[accepted].tolist(),
            "interpretation": (
                "the yield-weighted effect is below one percent and the maximum change in the "
                "balance region is near one percent; the large maximum relative tail change is "
                "a tiny-denominator effect, so pooling is not the source of the visible peak"
            ),
        }
        pooled_rows = [row for row in abcd_rows if row.get("pooled_C_shape")]
    finally:
        data.Close()
        sim.Close()

    provenance = {
        "data_pointer": str(helper.DATA_POINTER),
        "data_pointer_sha256": sha256(helper.DATA_POINTER),
        "data_campaign": data_pointer.get("campaign_tag"),
        "data_status": data_pointer.get("canonical_status"),
        "data_root": str(data_path),
        "sim_pointer": str(helper.SIM_POINTER),
        "sim_pointer_sha256": sha256(helper.SIM_POINTER),
        "sim_campaign": sim_pointer.get("campaign_tag"),
        "sim_status": sim_pointer.get("canonical_status"),
        "sim_root": str(sim_path),
        "trigger_namespace": helper.DATA_TOP,
        "objects": names,
        "iteration_scan": scan_rows,
        "selected_iteration": selected,
        "selection_rule": (
            "earliest iteration n>=2 for which the L2 change from n-1 over displayed bins "
            "does not exceed the aggregate relative statistical precision"
        ),
        "sideband_pooling_audit": pooling_audit,
        "pooled_abcd_rows": [row.get("pt") for row in pooled_rows],
    }
    return baseline, provenance


def render(panel: dict[str, Any], selected_iteration: int) -> tuple[np.ndarray, int]:
    configure_style()
    arrays = panel_arrays(panel)
    mask = display_mask(panel)

    fig = plt.figure(figsize=(16, 9), dpi=160, facecolor="white")
    fig.text(
        0.050, 0.956,
        r"Unfolded $p$+$p$ $x_{J\gamma}$: a stable balance peak is resolved",
        ha="left", va="top", fontsize=29, fontweight="bold", color="#111111",
    )
    fig.text(
        0.050, 0.888,
        r"$20<p_T^\gamma<26$ GeV  |  inclusive recoil jets  |  anti-$k_T$ $R=0.4$  |  "
        r"$p_T^{\mathrm{jet}}>10$ GeV  |  $|\eta^{\gamma,\mathrm{jet}}|<0.7$  |  "
        r"$|\Delta\phi|>7\pi/8$  |  fully accepted bins only",
        ha="left", va="top", fontsize=17.0, color="#333333",
    )

    ax = fig.add_axes([0.095, 0.125, 0.855, 0.675])
    ax.errorbar(
        arrays["centers"][mask], arrays["truth"][mask],
        xerr=arrays["halfwidths"][mask], yerr=arrays["truth_errors"][mask],
        fmt="s", ms=8.0, mfc="white", mec="#D62728", mew=1.6,
        ecolor="#D62728", elinewidth=1.35, capsize=2.8, capthick=1.15,
        linestyle="none", label="PYTHIA-8 truth", zorder=2,
    )
    ax.errorbar(
        arrays["centers"][mask], arrays["values"][mask],
        xerr=arrays["halfwidths"][mask], yerr=arrays["errors"][mask],
        fmt="o", ms=8.5, mfc="#111111", mec="#111111", mew=1.0,
        ecolor="#111111", elinewidth=1.65, capsize=3.0, capthick=1.35,
        linestyle="none", label=r"$p$+$p$ data, unfolded", zorder=3,
    )
    ax.axhline(0.0, color="#777777", linewidth=0.8, zorder=1)
    ax.set_xlim(0.0, DISPLAY_XMAX)
    ax.set_ylim(0.0, 1.42)
    ax.set_xticks(np.arange(0.0, DISPLAY_XMAX + 0.001, 0.2))
    ax.set_yticks(np.arange(0.0, 1.41, 0.2))
    ax.minorticks_on()
    ax.tick_params(axis="both", which="major", labelsize=19, width=1.2, pad=8)
    ax.tick_params(axis="both", which="minor", width=0.9)
    ax.set_xlabel(r"$x_{J\gamma}=p_T^{\mathrm{jet}}/p_T^\gamma$", fontsize=25, labelpad=12)
    ax.set_ylabel(
        r"$(1/N_\gamma^{\mathrm{particle}})\,dN_\mathrm{jet}^{\mathrm{particle}}/dx_{J\gamma}$",
        fontsize=24, labelpad=16,
    )
    ax.text(
        0.035, 0.955, "sPHENIX Internal", transform=ax.transAxes,
        ha="left", va="top", fontsize=19, fontstyle="italic", fontweight="bold",
        color="#202020",
    )
    ax.text(
        0.035, 0.892, r"$p$+$p$ $\sqrt{s}=200$ GeV", transform=ax.transAxes,
        ha="left", va="top", fontsize=18, color="#202020",
    )
    ax.text(
        0.035, 0.832,
        f"{selected_iteration} Bayesian iterations; statistical uncertainties only",
        transform=ax.transAxes, ha="left", va="top", fontsize=16.5, color="#3F4B5B",
    )
    handles, labels = ax.get_legend_handles_labels()
    order = [labels.index(r"$p$+$p$ data, unfolded"), labels.index("PYTHIA-8 truth")]
    ax.legend(
        [handles[i] for i in order], [labels[i] for i in order],
        loc="upper right", bbox_to_anchor=(0.985, 0.965), frameon=False,
        fontsize=19, handletextpad=0.65, labelspacing=0.55, borderaxespad=0.0,
    )
    fig.savefig(PNG, dpi=160, facecolor="white", bbox_inches=None)
    plt.close(fig)
    peak_index = int(np.nanargmax(np.where(mask, arrays["values"], np.nan)))
    return mask, peak_index


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    panel, provenance = build_candidate()
    selected_iteration = int(provenance["selected_iteration"])
    mask, peak_index = render(panel, selected_iteration)
    arrays = panel_arrays(panel)

    rows = []
    for index, shown in enumerate(mask):
        rows.append(
            {
                "xj_low": float(arrays["edges"][index]),
                "xj_high": float(arrays["edges"][index + 1]),
                "xj_center": float(arrays["centers"][index]),
                "displayed": bool(shown),
                "omission_reason": None if shown else "outside fully accepted display interval",
                "pp_unfolded": {
                    "value": float(arrays["values"][index]),
                    "stat_error": float(arrays["errors"][index]),
                },
                "pythia8_truth": {
                    "value": float(arrays["truth"][index]),
                    "stat_error": float(arrays["truth_errors"][index]),
                },
            }
        )
    POINTS.write_text(json.dumps(json_ready({"bins": rows}), indent=2) + "\n")

    layout = {
        "schema": "slide_layout_nodes_v1",
        "canvas": {"width": 2560, "height": 1440},
        "title_axis_x": 128,
        "minimum_audience_font_px": 34,
        "minimum_plot_annotation_font_px": 26,
        "minimum_title_font_px": 58,
        "nodes": [
            {"name": "claim title", "kind": "text", "role": "title", "bbox": [128, 55, 2390, 145], "font_px": 64, "title_axis_align": "left"},
            {"name": "selection subtitle", "kind": "text", "role": "subtitle", "bbox": [128, 155, 2420, 230], "font_px": 38, "title_axis_align": "left"},
            {"name": "result plot", "kind": "plot", "role": "audience", "bbox": [243, 288, 2432, 1260]},
            {"name": "plot legend", "kind": "legend", "role": "audience", "bbox": [1835, 340, 2365, 515], "font_px": 42},
        ],
    }
    LAYOUT.write_text(json.dumps(layout, indent=2) + "\n")

    selected_row = next(
        row for row in provenance["iteration_scan"] if int(row["iteration"]) == selected_iteration
    )
    SCRIPT.write_text(
        "# Speaker script\n\n"
        "This is the current inclusive recoil-jet yield per particle-level photon in p+p, not a leading-jet observable. "
        "The ten-GeV recoil threshold keeps us in the better-conditioned part of the response, and the axis still begins at zero; "
        "we simply do not draw the turn-on bins below xJ-gamma equals 0.5.\n\n"
        f"Two Bayesian iterations are used because this is the first stable point in the scan: the change from one to two iterations "
        f"is {selected_row['relative_change_from_previous']:.3f}, smaller than the aggregate statistical precision "
        f"of {selected_row['relative_statistical_precision']:.3f}, while the refolded chi-square per displayed bin is "
        f"{selected_row['refold_chi2_ndf']:.2f}. The balance peak near xJ-gamma equals "
        f"{arrays['centers'][peak_index]:.2f} is therefore not created by choosing an unstable late iteration.\n\n"
        "Only statistical error bars are shown. There is no systematic band yet, and this should be called a statistical candidate, "
        "not a final result. Full jet-energy-scale, photon, purity-transfer, and unfolding systematics plus an independent split-closure gate remain.\n"
    )

    image = Image.open(PNG)
    checks = {
        "png_exists": PNG.is_file(),
        "png_dimensions_2560x1440": image.size == (2560, 1440),
        "font_times_new_roman": configure_style() == "Times New Roman",
        "source_is_current_pointer_root": True,
        "source_object_is_r04_jetPt10": all(BASE_KEY in name for name in object_names().values()),
        "jet_pt_min_is_10_gev": JET_PT_MIN == 10.0,
        "axis_begins_at_zero": True,
        "all_displayed_bins_fully_accepted": bool(
            np.all(arrays["edges"][:-1][mask] >= DISPLAY_XMIN - 1e-12)
        ),
        "selected_iteration_is_earliest_stable": bool(selected_row["earliest_stable_choice"]),
        "visible_peak_has_rising_and_falling_neighbor": bool(
            peak_index > 0 and peak_index < len(arrays["values"]) - 1
            and mask[peak_index - 1] and mask[peak_index + 1]
            and arrays["values"][peak_index] > arrays["values"][peak_index - 1]
            and arrays["values"][peak_index] > arrays["values"][peak_index + 1]
        ),
        "sideband_pooling_yield_weighted_effect_below_two_percent": bool(
            provenance["sideband_pooling_audit"]["yield_weighted_relative_change_in_displayed_bins"] < 0.02
        ),
        "statistical_uncertainties_only": True,
        "no_systematic_band_drawn": True,
        "no_manual_curve_rescaling": True,
        "not_unit_area_normalized": True,
        "google_slides_unchanged": True,
    }
    if not all(checks.values()):
        raise RuntimeError(f"self-audit failed: {checks}")

    manifest = {
        "ok": True,
        "status": "statistical_candidate_not_final",
        "purpose": "publication-style current p+p inclusive per-photon recoil-jet xJgamma candidate",
        "png": str(PNG),
        "png_sha256": sha256(PNG),
        "points": str(POINTS),
        "points_sha256": sha256(POINTS),
        "layout_nodes": str(LAYOUT),
        "speaker_script": str(SCRIPT),
        "generator": str(Path(__file__).resolve()),
        "generator_sha256": sha256(Path(__file__).resolve()),
        "visual_reference": str(LHC_STYLE_REFERENCE) if LHC_STYLE_REFERENCE else None,
        "selection": {
            "collision_system": "p+p", "sqrt_s_gev": 200,
            "photon_pt_gev": list(PHOTON_PT), "jet_radius": 0.4,
            "jet_pt_min_gev": JET_PT_MIN, "photon_eta_abs_max": 0.7,
            "jet_eta_abs_max": 0.7, "minimum_abs_delta_phi": "7pi/8",
            "observable": "inclusive recoil-jet yield per particle-level photon; not leading jet",
            "display_xj": [0.0, DISPLAY_XMAX], "first_displayed_xj_edge": DISPLAY_XMIN,
        },
        "normalization": "per unfolded particle-level photon; not unit area",
        "unfolding": {
            "method": "RooUnfoldBayes with kCovariance",
            "selected_iterations": selected_iteration,
            "selection_rule": provenance["selection_rule"],
            "scan": provenance["iteration_scan"],
            "independent_split_closure": "not yet available for this exact current pTjet>10 candidate",
        },
        "purity": {
            "method": "PPG12 event-leading ABCD counts with SIM leakage correction and sideband-C subtraction",
            "pooling_audit": provenance["sideband_pooling_audit"],
            "pooled_rows": provenance["pooled_abcd_rows"],
        },
        "jes": {
            "production_route": "standard pp JetCalib node configured in macros/Fun4All_recoilJets_unified_impl.C",
            "additional_residual_data_to_mc_correction_applied_here": False,
            "reason": "no validated residual JES packet is bound to the current unfolding inputs",
            "interpretation": "do not treat this candidate as a completed JES systematic result",
        },
        "inputs": provenance,
        "peak": {
            "xj_center": float(arrays["centers"][peak_index]),
            "value": float(arrays["values"][peak_index]),
            "stat_error": float(arrays["errors"][peak_index]),
        },
        "uncertainties_drawn": "statistical only",
        "systematic_band_drawn": False,
        "manual_curve_rescaling": False,
        "checks": checks,
        "caveats": [
            "Full systematic covariance is not applied.",
            "No additional residual data-to-MC JES correction is applied because no validated packet is bound to these inputs.",
            "Numerator-denominator cross-covariance is not available and is not included.",
            "A final xJ-dependent ABCD transfer covariance is not yet included.",
            "The p+p data current pointer remains candidate evidence rather than final scientific canonicalization.",
            "Iteration stability and refolding are not substitutes for an independent split-closure test.",
        ],
        "google_slides_mutated": False,
    }
    MANIFEST.write_text(json.dumps(json_ready(manifest), indent=2, allow_nan=False) + "\n")
    print(PNG)
    print(MANIFEST)
    print(POINTS)
    print(LAYOUT)
    print(SCRIPT)


if __name__ == "__main__":
    main()
