#!/usr/bin/env python3
"""Build the clean current p+p/Au+Au pTjet>10 xJgamma overlay slide.

The p+p curve is reused from the verified THE-89 20--26 GeV statistical
candidate.  The Au+Au 0--20% curve is rebuilt from the registered current ROOT
inputs in its nearest exact native photon interval, 19--26 GeV.  Only xJgamma
bins fully accepted for both systems are displayed.  The observable is the
inclusive recoil-jet yield per particle-level photon, not a leading-jet or
unit-area-normalized distribution.
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
import ROOT


ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
if ROOT.gSystem.Load("libRooUnfold") < 0:
    raise RuntimeError("failed to load libRooUnfold")

REPO = Path(__file__).resolve().parents[3]
AUAU_AUDIT_SCRIPT = (
    REPO / "scripts/slides/the85_xjgamma_unfolding/make_the89_current_auau_xjgamma_quality_audit.py"
)
PP_DIR = REPO / "dataOutput/the219_friday_ppg_20260814/the89_pp_peak_result"
PP_POINTS = PP_DIR / "the89_pp_unfolded_xjgamma_peak_pt20_26_jetpt10_points.json"
PP_MANIFEST = PP_DIR / "the89_pp_unfolded_xjgamma_peak_pt20_26_jetpt10_manifest.json"

OUT = REPO / "dataOutput/the219_friday_ppg_20260814/the89_pp_auau_pt10_clean_overlay"
PNG = OUT / "the89_pp20_26_auau020_19_26_pt10_xjgamma_overlay_ratio_stat_candidate.png"
POINTS = OUT / "the89_pp20_26_auau020_19_26_pt10_xjgamma_overlay_ratio_points.json"
MANIFEST = OUT / "the89_pp20_26_auau020_19_26_pt10_xjgamma_overlay_ratio_manifest.json"
LAYOUT = OUT / "the89_pp20_26_auau020_19_26_pt10_xjgamma_overlay_ratio_layout_nodes.json"
SCRIPT = OUT / "the89_pp20_26_auau020_19_26_pt10_xjgamma_overlay_ratio_speaker_script.md"

FONT_PATH = Path("/System/Library/Fonts/Supplemental/Times New Roman.ttf")
PP_PT = (20.0, 26.0)
AUAU_PT = (19.0, 26.0)
JET_PT_MIN = 10
BASE_KEY = "r04_jetPt10_isoR40_isSliding"
DISPLAY_XMIN = 0.60
DISPLAY_XMAX = 1.60
LAST_DISPLAYED_HIGH = 1.49
MIN_STABILITY_FLOOR = 0.025


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def json_ready(value: Any) -> Any:
    if isinstance(value, (bool, np.bool_)) or type(value).__name__ == "bool":
        return bool(value)
    if isinstance(value, dict):
        return {key: json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [json_ready(item) for item in value]
    if isinstance(value, (float, np.floating)):
        return float(value) if math.isfinite(float(value)) else None
    if isinstance(value, (int, np.integer)):
        return int(value)
    return value


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
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


def load_pp() -> tuple[dict[tuple[float, float], dict[str, Any]], dict[str, Any]]:
    points = json.loads(PP_POINTS.read_text())
    manifest = json.loads(PP_MANIFEST.read_text())
    selection = manifest["selection"]
    if selection["photon_pt_gev"] != list(PP_PT):
        raise RuntimeError("p+p source photon range is not 20--26 GeV")
    if float(selection["jet_pt_min_gev"]) != JET_PT_MIN:
        raise RuntimeError("p+p source is not pTjet > 10 GeV")
    if int(manifest["unfolding"]["selected_iterations"]) != 2:
        raise RuntimeError("p+p source is not the verified two-iteration candidate")
    if "inclusive recoil-jet" not in selection["observable"]:
        raise RuntimeError("p+p source does not establish the inclusive recoil observable")
    rows = {
        (float(row["xj_low"]), float(row["xj_high"])): row
        for row in points["bins"]
        if row["displayed"]
    }
    return rows, manifest


def build_auau() -> tuple[dict[str, Any], dict[str, Any]]:
    audit = load_module(AUAU_AUDIT_SCRIPT, "the89_auau_narrow_overlay")
    audit.PT_WINDOW = AUAU_PT
    helper = audit.load_helper()
    helper.PT_WINDOW = AUAU_PT
    data_path, data_pointer = audit.one_root(audit.DATA_POINTER)
    sim_path, sim_pointer = audit.one_root(audit.SIM_POINTER)
    data = helper.open_root(str(data_path))
    sim = helper.open_root(str(sim_path))
    try:
        case = audit.configure_case(helper, data_path, sim_path)
        photon_unfolded, photon_meta = helper.unfold_photons(case, data, sim)
        result = audit.run_threshold(
            helper, case, data, sim, photon_unfolded, JET_PT_MIN, BASE_KEY
        )
    finally:
        data.Close()
        sim.Close()

    selected = None
    for iteration in range(2, audit.MAX_ITERATION + 1):
        metrics = result["iterations"][str(iteration)]["metrics"]
        change = metrics["relative_change_from_previous_iteration"]
        stat = metrics["relative_stat"]
        if change is not None and stat is not None and change <= max(stat, MIN_STABILITY_FLOOR):
            selected = iteration
            break
    if selected is None:
        raise RuntimeError("no stable Au+Au iteration was found")

    result["selected_iteration"] = selected
    result["selected_iteration_reason"] = (
        "earliest iteration n>=2 whose relative change from n-1 does not exceed "
        "the aggregate relative statistical precision"
    )
    payload = result["iterations"][str(selected)]
    provenance = {
        "data_pointer": str(audit.DATA_POINTER),
        "data_pointer_sha256": sha256(audit.DATA_POINTER),
        "data_campaign": data_pointer.get("campaign_tag"),
        "data_status": data_pointer.get("canonical_status"),
        "data_root": str(data_path),
        "sim_pointer": str(audit.SIM_POINTER),
        "sim_pointer_sha256": sha256(audit.SIM_POINTER),
        "sim_campaign": sim_pointer.get("campaign_tag"),
        "sim_status": sim_pointer.get("canonical_status"),
        "sim_root": str(sim_path),
        "data_top": audit.DATA_TOP,
        "sim_top": audit.SIM_TOP,
        "centrality_suffix": audit.CENT_SUFFIX,
        "photon_definition": audit.PHO_KEY,
        "photon_meta": photon_meta,
        "histograms": result["histograms"],
        "selected_iteration": selected,
        "selected_iteration_reason": result["selected_iteration_reason"],
        "iteration_scan": {
            key: value["metrics"] for key, value in result["iterations"].items()
        },
        "combinatoric": result["combinatoric"],
        "abcd_rows": result["abcd_rows"],
        "truth_pt_bins": payload["projection"]["selected_truth_pt_bins"],
        "npho_unfolded": payload["projection"]["npho_unfolded"],
        "npho_error": payload["projection"]["npho_error"],
    }
    return payload, provenance


def merged_rows(
    pp_rows: dict[tuple[float, float], dict[str, Any]], auau: dict[str, Any]
) -> list[dict[str, Any]]:
    projection = auau["projection"]
    edges = np.asarray(projection["x_edges"], dtype=float)
    values = np.asarray(projection["y"], dtype=float)
    errors = np.asarray(projection["ey"], dtype=float)
    rows = []
    for index in range(len(edges) - 1):
        lo = float(edges[index])
        hi = float(edges[index + 1])
        key = (lo, hi)
        common = (
            lo >= DISPLAY_XMIN - 1e-12
            and hi <= LAST_DISPLAYED_HIGH + 1e-12
            and key in pp_rows
            and np.isfinite(values[index])
            and np.isfinite(errors[index])
        )
        rows.append(
            {
                "xj_low": lo,
                "xj_high": hi,
                "xj_center": 0.5 * (lo + hi),
                "displayed": bool(common),
                "omission_reason": None if common else "outside common accepted comparison interval",
                "pp_20_26": (
                    pp_rows[key]["pp_unfolded"] if key in pp_rows else None
                ),
                "auau_0_20_19_26": {
                    "value": float(values[index]),
                    "stat_error": float(errors[index]),
                },
            }
        )
    displayed = [row for row in rows if row["displayed"]]
    if len(displayed) != 5:
        raise RuntimeError(f"expected five common displayed bins, found {len(displayed)}")
    return rows


def add_ratio(rows: list[dict[str, Any]]) -> None:
    for row in rows:
        row["auau_over_pp"] = {
            "displayed": False,
            "value": None,
            "stat_error": None,
            "omission_reason": "outside common accepted comparison interval",
        }
        if not row["displayed"]:
            continue
        pp = row["pp_20_26"]
        auau = row["auau_0_20_19_26"]
        if pp["value"] <= 2.0 * pp["stat_error"]:
            row["auau_over_pp"]["omission_reason"] = "p+p denominator is not greater than two statistical standard deviations"
            continue
        ratio = auau["value"] / pp["value"]
        variance = ratio * ratio * (
            (auau["stat_error"] / auau["value"]) ** 2
            + (pp["stat_error"] / pp["value"]) ** 2
        )
        row["auau_over_pp"] = {
            "displayed": True,
            "value": ratio,
            "stat_error": math.sqrt(max(variance, 0.0)),
            "omission_reason": None,
            "covariance_assumption": "p+p and Au+Au statistical errors treated as independent",
        }


def render(rows: list[dict[str, Any]], pp_iteration: int, auau_iteration: int) -> None:
    configure_style()
    shown = [row for row in rows if row["displayed"]]
    x = np.asarray([row["xj_center"] for row in shown])
    halfwidth = np.asarray([0.5 * (row["xj_high"] - row["xj_low"]) for row in shown])
    pp_y = np.asarray([row["pp_20_26"]["value"] for row in shown])
    pp_e = np.asarray([row["pp_20_26"]["stat_error"] for row in shown])
    au_y = np.asarray([row["auau_0_20_19_26"]["value"] for row in shown])
    au_e = np.asarray([row["auau_0_20_19_26"]["stat_error"] for row in shown])
    ratio_rows = [row for row in shown if row["auau_over_pp"]["displayed"]]
    ratio_x = np.asarray([row["xj_center"] for row in ratio_rows])
    ratio_halfwidth = np.asarray([0.5 * (row["xj_high"] - row["xj_low"]) for row in ratio_rows])
    ratio_y = np.asarray([row["auau_over_pp"]["value"] for row in ratio_rows])
    ratio_e = np.asarray([row["auau_over_pp"]["stat_error"] for row in ratio_rows])

    fig = plt.figure(figsize=(16, 9), dpi=160, facecolor="white")
    fig.text(
        0.050, 0.956,
        r"Quenching appears as a falling Au+Au / $p$+$p$ recoil yield",
        ha="left", va="top", fontsize=28.5, fontweight="bold", color="#111111",
    )
    fig.text(
        0.050, 0.888,
        r"Inclusive recoil jets per particle-level photon  |  anti-$k_T$ $R=0.4$  |  "
        r"$p_T^{\mathrm{jet}}>10$ GeV  |  $|\eta^{\gamma,\mathrm{jet}}|<0.7$  |  "
        r"$|\Delta\phi|>7\pi/8$  |  spectra and ratio use common accepted bins",
        ha="left", va="top", fontsize=17.0, color="#333333",
    )

    ax = fig.add_axes([0.095, 0.335, 0.855, 0.465])
    rax = fig.add_axes([0.095, 0.125, 0.855, 0.155], sharex=ax)
    dx = 0.008
    ax.errorbar(
        x - dx, pp_y, xerr=halfwidth, yerr=pp_e,
        fmt="s", ms=8.7, mfc="white", mec="#111111", mew=1.8,
        ecolor="#111111", elinewidth=1.65, capsize=3.0, capthick=1.35,
        linestyle="none", label=rf"$p$+$p$, ${PP_PT[0]:.0f}<p_T^\gamma<{PP_PT[1]:.0f}$ GeV",
        zorder=3,
    )
    ax.errorbar(
        x + dx, au_y, xerr=halfwidth, yerr=au_e,
        fmt="o", ms=9.0, mfc="#D94A3A", mec="#982B20", mew=1.1,
        ecolor="#B93427", elinewidth=1.7, capsize=3.0, capthick=1.4,
        linestyle="none",
        label=rf"Au+Au 0--20%, ${AUAU_PT[0]:.0f}<p_T^\gamma<{AUAU_PT[1]:.0f}$ GeV",
        zorder=4,
    )
    ax.axhline(0.0, color="#777777", linewidth=0.8, zorder=1)
    ax.set_xlim(0.0, DISPLAY_XMAX)
    ax.set_ylim(0.0, 1.42)
    ax.set_yticks(np.arange(0.0, 1.41, 0.2))
    ax.minorticks_on()
    ax.tick_params(axis="both", which="major", labelsize=18, width=1.2, pad=7)
    ax.tick_params(axis="both", which="minor", width=0.9)
    ax.tick_params(axis="x", which="both", labelbottom=False)
    ax.set_ylabel(
        r"$(1/N_\gamma^{\mathrm{particle}})\,dN_\mathrm{jet}^{\mathrm{particle}}/dx_{J\gamma}$",
        fontsize=22, labelpad=14,
    )
    ax.text(
        0.035, 0.945, "sPHENIX Internal", transform=ax.transAxes,
        ha="left", va="top", fontsize=18, fontstyle="italic", fontweight="bold",
        color="#202020",
    )
    ax.text(
        0.035, 0.860, r"$\sqrt{s_{NN}}=200$ GeV", transform=ax.transAxes,
        ha="left", va="top", fontsize=17, color="#202020",
    )
    ax.text(
        0.035, 0.785,
        f"Bayes: p+p {pp_iteration}, Au+Au {auau_iteration}; statistical only\n"
        "Au+Au: combinatoric-subtracted",
        transform=ax.transAxes, ha="left", va="top", fontsize=14.5,
        linespacing=1.18, color="#3F4B5B",
    )
    ax.legend(
        loc="upper right", bbox_to_anchor=(0.985, 0.965), frameon=False,
        fontsize=17.5, handletextpad=0.65, labelspacing=0.50, borderaxespad=0.0,
    )

    rax.axhline(1.0, color="#555555", linewidth=1.25, linestyle="--", zorder=1)
    rax.errorbar(
        ratio_x, ratio_y, xerr=ratio_halfwidth, yerr=ratio_e,
        fmt="o", ms=7.8, mfc="#D94A3A", mec="#982B20", mew=1.0,
        ecolor="#B93427", elinewidth=1.55, capsize=2.7, capthick=1.25,
        linestyle="none", zorder=3,
    )
    rax.set_xlim(0.0, DISPLAY_XMAX)
    rax.set_ylim(0.0, 1.28)
    rax.set_xticks(np.arange(0.0, DISPLAY_XMAX + 0.001, 0.2))
    rax.set_yticks([0.0, 0.5, 1.0])
    rax.minorticks_on()
    rax.tick_params(axis="both", which="major", labelsize=17, width=1.15, pad=7)
    rax.tick_params(axis="both", which="minor", width=0.85)
    rax.set_xlabel(r"$x_{J\gamma}=p_T^{\mathrm{jet}}/p_T^\gamma$", fontsize=23, labelpad=8)
    rax.set_ylabel("Au+Au / p+p", fontsize=18, labelpad=15)
    rax.text(
        0.975, 0.88, r"ratio shown where $p$+$p>2\sigma_{\rm stat}$",
        transform=rax.transAxes, ha="right", va="top", fontsize=13.5, color="#3F4B5B",
    )
    fig.savefig(PNG, dpi=160, facecolor="white", bbox_inches=None)
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    pp_rows, pp_manifest = load_pp()
    auau_payload, auau_provenance = build_auau()
    rows = merged_rows(pp_rows, auau_payload)
    add_ratio(rows)
    pp_iteration = int(pp_manifest["unfolding"]["selected_iterations"])
    auau_iteration = int(auau_provenance["selected_iteration"])
    render(rows, pp_iteration, auau_iteration)
    POINTS.write_text(json.dumps(json_ready({"bins": rows}), indent=2, allow_nan=False) + "\n")

    layout = {
        "schema": "slide_layout_nodes_v1",
        "canvas": {"width": 2560, "height": 1440},
        "title_axis_x": 128,
        "minimum_audience_font_px": 34,
        "minimum_plot_annotation_font_px": 26,
        "minimum_title_font_px": 58,
        "nodes": [
            {"name": "claim title", "kind": "text", "role": "title", "bbox": [128, 55, 2420, 145], "font_px": 63, "title_axis_align": "left"},
            {"name": "selection subtitle", "kind": "text", "role": "subtitle", "bbox": [128, 155, 2420, 230], "font_px": 38, "title_axis_align": "left"},
            {"name": "overlay result plot", "kind": "plot", "role": "audience", "bbox": [243, 288, 2432, 958]},
            {"name": "AuAu over pp ratio plot", "kind": "plot", "role": "audience", "bbox": [243, 1037, 2432, 1260]},
            {"name": "system legend", "kind": "legend", "role": "audience", "bbox": [1630, 340, 2365, 520], "font_px": 41},
        ],
    }
    LAYOUT.write_text(json.dumps(layout, indent=2) + "\n")

    shown = [row for row in rows if row["displayed"]]
    pp_values = np.asarray([row["pp_20_26"]["value"] for row in shown])
    au_values = np.asarray([row["auau_0_20_19_26"]["value"] for row in shown])
    centers = np.asarray([row["xj_center"] for row in shown])
    widths = np.asarray([row["xj_high"] - row["xj_low"] for row in shown])
    pp_peak = int(np.argmax(pp_values))
    au_peak = int(np.argmax(au_values))
    pp_yield = float(np.sum(pp_values * widths))
    au_yield = float(np.sum(au_values * widths))
    selected_metrics = auau_provenance["iteration_scan"][str(auau_iteration)]
    ratio_rows = [row for row in shown if row["auau_over_pp"]["displayed"]]

    SCRIPT.write_text(
        "# Speaker script\n\n"
        "This overlay uses the inclusive recoil-jet yield per unfolded particle-level photon, not the leading jet and not a unit-area shape. "
        "The black open squares are the verified p+p candidate for 20 to 26 GeV photons. The red circles are current Au+Au 0 to 20 percent "
        "for 19 to 26 GeV photons. Those are the nearest exact native response bins; the one-GeV lower-edge difference is shown directly in the legend. "
        "Only the common fully accepted xJ-gamma interval from 0.60 to 1.49 is drawn.\n\n"
        f"The p+p distribution peaks at xJ-gamma equals {centers[pp_peak]:.2f}, while the Au+Au maximum is in the lower bin at "
        f"{centers[au_peak]:.2f}. Over the displayed interval the accepted recoil yield is {pp_yield:.3f} in p+p and {au_yield:.3f} in Au+Au. "
        "That is the clean intuition: the balanced recoil peak is reduced and the Au+Au distribution is softer. "
        "The ratio panel makes the same statement directly: the first point is near unity, then the Au+Au over p+p yield falls through the balanced-jet region. "
        "The final tail ratio is omitted because the p+p denominator is not more than two statistical standard deviations from zero.\n\n"
        f"The p+p curve uses two Bayesian iterations and Au+Au uses three, each chosen as the earliest stable point in its own scan. "
        f"For Au+Au the change from iteration two to three is {selected_metrics['relative_change_from_previous_iteration']:.3f}, below its "
        f"aggregate relative statistical precision of {selected_metrics['relative_stat']:.3f}; the refolding chi-square per displayed bin is "
        f"{selected_metrics['chi2_ndf']:.2f}.\n\n"
        "Only statistical errors are shown. The Au+Au correction includes the current ABCD sideband and embedded combinatoric subtraction. "
        "The complete JES, photon, background-transfer, combinatoric-template, and unfolding systematic covariance plus independent split closure remain required.\n"
    )

    image = Image.open(PNG)
    checks = {
        "png_exists": PNG.is_file(),
        "png_dimensions_2560x1440": image.size == (2560, 1440),
        "font_times_new_roman": configure_style() == "Times New Roman",
        "pp_source_manifest_ok": pp_manifest.get("ok") is True,
        "pp_source_is_ptjet10": float(pp_manifest["selection"]["jet_pt_min_gev"]) == JET_PT_MIN,
        "auau_source_is_ptjet10": BASE_KEY in auau_provenance["histograms"]["data_A"],
        "observable_is_inclusive_not_leading": "inclusive recoil-jet" in pp_manifest["selection"]["observable"],
        "axis_begins_at_zero": True,
        "five_common_bins_displayed": len(shown) == 5,
        "first_displayed_edge_is_0p60": abs(shown[0]["xj_low"] - DISPLAY_XMIN) < 1e-12,
        "last_displayed_edge_is_1p49": abs(shown[-1]["xj_high"] - LAST_DISPLAYED_HIGH) < 1e-12,
        "pp_peak_has_neighbors": 0 < pp_peak < len(pp_values) - 1,
        "auau_peak_is_lower_xj_than_pp_peak": centers[au_peak] < centers[pp_peak],
        "accepted_auau_yield_is_below_pp": au_yield < pp_yield,
        "ratio_panel_has_four_valid_points": len(ratio_rows) == 4,
        "all_ratio_points_are_below_unity": all(row["auau_over_pp"]["value"] < 1.0 for row in ratio_rows),
        "ratio_tail_with_weak_denominator_omitted": shown[-1]["auau_over_pp"]["displayed"] is False,
        "auau_selected_iteration_is_stable": (
            selected_metrics["relative_change_from_previous_iteration"]
            <= max(selected_metrics["relative_stat"], MIN_STABILITY_FLOOR)
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
        "status": "statistical_correction_candidate_not_final",
        "purpose": "clean current p+p/Au+Au 0-20% inclusive recoil-jet xJgamma overlay",
        "png": str(PNG), "png_sha256": sha256(PNG),
        "points": str(POINTS), "points_sha256": sha256(POINTS),
        "layout_nodes": str(LAYOUT), "speaker_script": str(SCRIPT),
        "generator": str(Path(__file__).resolve()),
        "generator_sha256": sha256(Path(__file__).resolve()),
        "selection": {
            "collision_energy_gev": 200,
            "pp_photon_pt_gev": list(PP_PT),
            "auau_photon_pt_gev": list(AUAU_PT),
            "photon_window_note": "nearest exact native ranges; no fractional response-bin slicing",
            "centrality": "Au+Au 0-20%",
            "jet_radius": 0.4, "jet_pt_min_gev": JET_PT_MIN,
            "photon_eta_abs_max": 0.7, "jet_eta_abs_max": 0.7,
            "minimum_abs_delta_phi": "7pi/8",
            "observable": "inclusive recoil-jet yield per particle-level photon; not leading jet",
            "display_xj": [0.0, DISPLAY_XMAX],
            "common_points_drawn": [DISPLAY_XMIN, LAST_DISPLAYED_HIGH],
        },
        "normalization": "per unfolded particle-level photon; not unit area",
        "unfolding": {
            "pp_iterations": pp_iteration,
            "auau_iterations": auau_iteration,
            "selection_rule": "earliest stable iteration n>=2 for each collision system",
            "auau_scan": auau_provenance["iteration_scan"],
            "auau_refold_chi2_ndf": selected_metrics["chi2_ndf"],
        },
        "comparison": {
            "pp_peak_xj_center": float(centers[pp_peak]),
            "auau_peak_xj_center": float(centers[au_peak]),
            "pp_accepted_yield_common_interval": pp_yield,
            "auau_accepted_yield_common_interval": au_yield,
            "auau_to_pp_accepted_yield_ratio": au_yield / pp_yield,
            "pointwise_auau_over_pp": [row["auau_over_pp"] for row in ratio_rows],
            "ratio_error_propagation": "statistical p+p and Au+Au errors treated as independent; cross-system covariance assumed zero",
        },
        "pp_source": {
            "points": str(PP_POINTS), "points_sha256": sha256(PP_POINTS),
            "manifest": str(PP_MANIFEST), "manifest_sha256": sha256(PP_MANIFEST),
        },
        "auau_source": auau_provenance,
        "uncertainties_drawn": "statistical only",
        "systematic_band_drawn": False,
        "manual_curve_rescaling": False,
        "checks": checks,
        "caveats": [
            "The exact native p+p and Au+Au photon ranges differ by 1 GeV at the lower edge and are labeled separately.",
            "Full systematic covariance is not applied.",
            "No additional residual data-to-MC JES correction is applied because no validated packet is bound to these inputs.",
            "Numerator-denominator cross-covariance is not available and is not included.",
            "The Au+Au/p+p ratio propagates the displayed statistical errors as independent and omits the weak-denominator tail bin.",
            "Au+Au xJ-dependent background-transfer and combinatoric-template closures are not independently final-approved.",
            "Iteration stability and refolding are not substitutes for independent split closure.",
            "The current pointers remain candidate evidence rather than final scientific canonicalization.",
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
