#!/usr/bin/env python3
"""Current detector-level leading-xJ overlay for THE-219.

Reproduce the visual logic of the supplied pp/Au+Au leading-jet reference with
the current registered RecoilJets products.  The current Au+Au products expose
0-20%, 20-50%, and 50-80% centrality.  The data projection uses the
photon-10-plus-MBD trigger directory selected by AnalyzeRecoilJets.h, not the
MBD-only directory.

The plotted object is the event-leading JES3 reconstruction-level family.  It
is not background-subtracted or unfolded.  Each curve is shape-normalized on
0.5 < xJgamma < 2.0, matching the presentation convention used for the
existing THE-219 leading-jet parity candidate.  The blue MC reference uses
the dedicated reconstructed photon-plus-jet pair truth-tagged JES3 family;
no pure-truth curve is drawn.
"""

from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import math
from pathlib import Path
from typing import Any

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


REPO = Path(__file__).resolve().parents[3]
HELPER_PATH = REPO / "scripts/slides/pp_currentian/xjgamma/make_the221_pp_auau_leading_match_slide.py"
OUT_DIR = REPO / "dataOutput/the219_friday_ppg_20260814/current_leading_xj_centrality_overlay"
OUT_DIR.mkdir(parents=True, exist_ok=True)
SCAN_DIR = OUT_DIR / "jet_pt_min_scan"
SCAN_DIR.mkdir(parents=True, exist_ok=True)
JET_PT_MIN_SCAN = (5.0, 7.0, 10.0, 12.0)

PP_PT_RANGE = (16.0, 18.0)
AUAU_PT_RANGE = (15.0, 17.0)
DISPLAY_RANGE = (0.20, 2.0)
CENTRALITY_KEYS = ("auau_50_80", "auau_0_20")
AUAU_DATA_TRIGGER_DIR = "photon_10_plus_MBD_NS_geq_2_vtx_lt_150"
PP_TAGGED_RECO_OBJECT = "SIM/h_JES3RecoTruthTagged_pT_xJ_alpha_r04"
AUAU_TAGGED_RECO_OBJECT = "SIM/h_JES3RecoTruthTagged_pT_xJ_alpha_r04_isoR40_isSliding_cent_{centrality}"


def artifact_paths(jet_pt_min: float) -> dict[str, Path | str]:
    threshold_tag = f"{jet_pt_min:g}".replace(".", "p")
    stem = f"the219_clean_leading_xj_pp16_18_auau15_17_jetptmin{threshold_tag}"
    return {
        "stem": stem,
        "png": SCAN_DIR / f"{stem}.png",
        "csv": SCAN_DIR / f"{stem}_points.csv",
        "json_points": SCAN_DIR / f"{stem}_points.json",
        "manifest": SCAN_DIR / f"{stem}_manifest.json",
        "speaker_script": SCAN_DIR / f"{stem}_speaker_script.md",
        "layout": SCAN_DIR / f"{stem}_layout_nodes.json",
    }


def rkey_for_threshold(jet_pt_min: float) -> str:
    if math.isclose(jet_pt_min, 5.0, rel_tol=0.0, abs_tol=1e-9):
        return "r04"
    return f"r04_jetPt{jet_pt_min:g}".replace(".", "p")


def load_helper():
    spec = importlib.util.spec_from_file_location("the221_leading_helper", HELPER_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load helper: {HELPER_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def finite(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: finite(item) for key, item in value.items()}
    if isinstance(value, list):
        return [finite(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def px_bbox(left: float, bottom: float, width: float, height: float) -> list[float]:
    return [2560 * left, 1440 * (1 - bottom - height), 2560 * (left + width), 1440 * (1 - bottom)]


def build_payload(jet_pt_min: float) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    helper = load_helper()
    pointers = {key: json.loads(path.read_text()) for key, path in helper.POINTERS.items()}
    files = {
        key: helper.open_root(Path(pointer["root_paths"][0]))
        for key, pointer in pointers.items()
    }
    try:
        specs = {spec["key"]: dict(spec) for spec in helper.PANEL_SPECS}
        for key in ("auau_0_20", "auau_20_50", "auau_50_80"):
            object_name = specs[key]["data_object"].split("/", 1)[1]
            specs[key]["data_object"] = f"{AUAU_DATA_TRIGGER_DIR}/{object_name}"
            specs[key]["reco_object"] = AUAU_TAGGED_RECO_OBJECT.format(centrality=key.removeprefix("auau_"))
        specs["pp"]["reco_object"] = PP_TAGGED_RECO_OBJECT
        rkey = rkey_for_threshold(jet_pt_min)
        for spec in specs.values():
            for object_key in ("data_object", "reco_object", "truth_object"):
                spec[object_key] = spec[object_key].replace("r04", rkey, 1)
        helper.PT_RANGE = PP_PT_RANGE
        panels = {"pp": helper.build_panel(specs["pp"], files)}
        helper.PT_RANGE = AUAU_PT_RANGE
        panels.update({key: helper.build_panel(specs[key], files) for key in CENTRALITY_KEYS})
    finally:
        for source in files.values():
            source.Close()
    return panels["pp"], {key: panels[key] for key in CENTRALITY_KEYS}, pointers


def visible_arrays(panel: dict[str, Any]) -> tuple[np.ndarray, ...]:
    edges = np.asarray(panel["xj_edges"], dtype=float)
    centers = 0.5 * (edges[:-1] + edges[1:])
    halfwidth = 0.5 * np.diff(edges)
    shown = (edges[:-1] >= DISPLAY_RANGE[0] - 1e-9) & (edges[1:] <= DISPLAY_RANGE[1] + 1e-9)
    return (
        edges,
        centers,
        halfwidth,
        shown,
        np.asarray(panel["data"], dtype=float),
        np.asarray(panel["data_stat"], dtype=float),
        np.asarray(panel["reco"], dtype=float),
        np.asarray(panel["reco_stat"], dtype=float),
    )


def render(
    pp: dict[str, Any],
    auau: dict[str, dict[str, Any]],
    jet_pt_min: float,
    png: Path,
) -> None:
    turn_on_xj = max(jet_pt_min / PP_PT_RANGE[0], jet_pt_min / AUAU_PT_RANGE[0])
    mpl.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "stix",
            "axes.linewidth": 1.15,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )
    fig = plt.figure(figsize=(16, 9), dpi=160, facecolor="white")
    fig.text(
        0.050,
        0.965,
        rf"Current leading-jet xJγ at $p_T^{{jet}}>{jet_pt_min:g}$ GeV: $p+p$ and Au+Au",
        ha="left",
        va="top",
        fontsize=27.0,
        fontweight="bold",
        color="#111111",
    )
    fig.text(
        0.052,
        0.888,
        r"Nearest native bins: $p+p$ $16<p_T^\gamma<18$ GeV; Au+Au $15<p_T^\gamma<17$ GeV  |  anti-$k_T$ $R=0.4$",
        ha="left",
        va="top",
        fontsize=15.5,
        color="#333333",
    )

    specs = (
        ("auau_50_80", "Peripheral: Au+Au 50–80%", 0.070),
        ("auau_0_20", "Central: Au+Au 0–20%", 0.535),
    )
    pp_arrays = visible_arrays(pp)
    y_candidates: list[float] = []
    for key, _title, _left in specs:
        au_arrays = visible_arrays(auau[key])
        for values, errors, shown in (
            (pp_arrays[4], pp_arrays[5], pp_arrays[3]),
            (au_arrays[4], au_arrays[5], au_arrays[3]),
            (au_arrays[6], au_arrays[7], au_arrays[3]),
        ):
            y_candidates.append(float(np.nanmax(values[shown] + errors[shown])))
    ymax = max(2.0, math.ceil(1.13 * max(y_candidates) / 0.25) * 0.25)

    red = "#d62728"
    blue = "#2358c6"
    black = "#111111"
    axes = []
    for index, (key, title, left) in enumerate(specs):
        ax = fig.add_axes([left, 0.215, 0.405, 0.545])
        axes.append(ax)
        edges, centers, halfwidth, shown, au_data, au_err, au_reco, au_reco_err = visible_arrays(auau[key])
        pp_edges, pp_centers, pp_halfwidth, pp_shown, pp_data, pp_err, _pp_reco, _pp_reco_err = pp_arrays
        if not np.allclose(edges, pp_edges, atol=1e-12, rtol=0.0):
            raise RuntimeError(f"xJ axes differ between pp and {key}")

        ax.axvspan(DISPLAY_RANGE[0], turn_on_xj, color="#e5e7eb", alpha=0.70, zorder=0)
        ax.axvline(turn_on_xj, color="#6b7280", linestyle="--", linewidth=1.2, zorder=1)
        ax.errorbar(
            pp_centers[pp_shown],
            pp_data[pp_shown],
            xerr=pp_halfwidth[pp_shown],
            yerr=pp_err[pp_shown],
            fmt="D",
            color=red,
            markerfacecolor="white",
            markeredgecolor=red,
            markersize=5.8,
            capsize=1.8,
            linewidth=1.1,
            zorder=4,
        )
        ax.errorbar(
            centers[shown],
            au_data[shown],
            xerr=halfwidth[shown],
            yerr=au_err[shown],
            fmt="o",
            color=black,
            markerfacecolor=black,
            markeredgecolor=black,
            markersize=5.6,
            capsize=1.9,
            linewidth=1.15,
            zorder=5,
        )
        ax.errorbar(
            centers[shown],
            au_reco[shown],
            xerr=halfwidth[shown],
            yerr=au_reco_err[shown],
            fmt="d",
            color=blue,
            markerfacecolor="white",
            markeredgecolor=blue,
            markersize=5.3,
            capsize=1.6,
            linewidth=1.0,
            zorder=3,
        )
        ax.set_xlim(*DISPLAY_RANGE)
        ax.set_ylim(0.0, ymax)
        ax.set_title(title, fontsize=21.0, fontweight="bold", pad=12)
        ax.set_xlabel(r"Leading-jet $x_{J\gamma}=p_T^{jet1}/p_T^\gamma$", fontsize=18.0, labelpad=8)
        if index == 0:
            ax.set_ylabel(r"Shape-normalized density", fontsize=18.0, labelpad=8)
        else:
            ax.tick_params(labelleft=False)
        ax.tick_params(labelsize=14.2, length=6)
        ax.tick_params(which="minor", length=3.2)
        ax.minorticks_on()
        ax.grid(axis="y", color="#d9dde3", linewidth=0.7, alpha=0.70)
        ax.text(
            0.035,
            0.955,
            r"$\bf{\it{sPHENIX}}$ Internal",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=16.5,
        )
        ax.text(
            0.035,
            0.885,
            rf"$p_T^{{jet}}>{jet_pt_min:g}$ GeV, $|\Delta\phi|\geq7\pi/8$",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=13.1,
            color="#333333",
        )
        ax.text(
            0.035,
            0.820,
            rf"shape norm entries: $p+p={pp['normalization_yields']['data']:.0f}$, Au+Au $={auau[key]['normalization_yields']['data']:.0f}$",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=12.2,
            color="#4b5563",
        )
        ax.text(
            turn_on_xj + 0.015,
            0.06 * ymax,
            rf"conservative ${jet_pt_min:g}/15$ turn-on",
            ha="left",
            va="bottom",
            fontsize=11.6,
            rotation=90,
            color="#555b65",
        )

    handles = [
        Line2D([0], [0], marker="D", linestyle="none", markerfacecolor="white", markeredgecolor=red, color=red, markersize=6.5, label=r"$p+p$ data"),
        Line2D([0], [0], marker="o", linestyle="none", markerfacecolor=black, markeredgecolor=black, color=black, markersize=6.5, label="Au+Au data"),
        Line2D([0], [0], marker="d", linestyle="none", markerfacecolor="white", markeredgecolor=blue, color=blue, markersize=6.5, label="MC reco pair, truth tagged (stat. only)"),
    ]
    axes[0].legend(handles=handles, loc="upper right", frameon=False, fontsize=13.0, handletextpad=0.45, labelspacing=0.55)

    fig.text(
        0.065,
        0.074,
        rf"• Central Au+Au retains a broad high-$x_{{J\gamma}}$ population after requiring $p_T^{{jet}}>{jet_pt_min:g}$ GeV.",
        ha="left",
        va="center",
        fontsize=15.6,
        fontweight="bold",
        color="#181818",
    )
    fig.text(
        0.065,
        0.038,
        "• Detector-level diagnostic: statistical errors only; no background subtraction, unfolding, residual correction, or systematics.",
        ha="left",
        va="center",
        fontsize=13.5,
        color="#333333",
    )
    fig.savefig(png, dpi=160, facecolor="white")
    plt.close(fig)


def write_sidecars(
    pp: dict[str, Any],
    auau: dict[str, dict[str, Any]],
    pointers: dict[str, Any],
    jet_pt_min: float,
    paths: dict[str, Path | str],
) -> None:
    stem = str(paths["stem"])
    png = Path(paths["png"])
    csv_path = Path(paths["csv"])
    json_points = Path(paths["json_points"])
    manifest_path = Path(paths["manifest"])
    speaker_script = Path(paths["speaker_script"])
    layout_path = Path(paths["layout"])
    turn_on_xj = max(jet_pt_min / PP_PT_RANGE[0], jet_pt_min / AUAU_PT_RANGE[0])
    payload = {"pp": pp, **auau}
    json_points.write_text(json.dumps(finite(payload), indent=2, allow_nan=False) + "\n")
    with csv_path.open("w", newline="") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=("sample", "xj_low", "xj_high", "data", "data_stat", "signal_mc_reco", "signal_mc_reco_stat", "displayed"),
        )
        writer.writeheader()
        for key, panel in payload.items():
            edges = panel["xj_edges"]
            for index, (low, high) in enumerate(zip(edges[:-1], edges[1:])):
                writer.writerow(
                    {
                        "sample": key,
                        "xj_low": low,
                        "xj_high": high,
                        "data": panel["data"][index],
                        "data_stat": panel["data_stat"][index],
                        "signal_mc_reco": panel["reco"][index],
                        "signal_mc_reco_stat": panel["reco_stat"][index],
                        "displayed": low >= DISPLAY_RANGE[0] - 1e-9 and high <= DISPLAY_RANGE[1] + 1e-9,
                    }
                )

    layout_path.write_text(
        json.dumps(
            {
                "schema": "slide_layout_nodes_v1",
                "slide": stem,
                "title_axis_x": 128,
                "minimum_audience_font_px": 29,
                "minimum_plot_annotation_font_px": 24,
                "minimum_title_font_px": 56,
                "nodes": [
                    {"name": "claim title", "kind": "text", "role": "title", "font_px": 60, "bbox": [128, 52, 2430, 122]},
                    {"name": "selection line", "kind": "text", "role": "audience", "font_px": 34, "bbox": [132, 150, 2360, 200]},
                    {"name": "peripheral panel", "kind": "panel", "role": "audience", "symmetry_group": "centrality contrast", "bbox": px_bbox(0.070, 0.215, 0.405, 0.545)},
                    {"name": "central panel", "kind": "panel", "role": "audience", "symmetry_group": "centrality contrast", "bbox": px_bbox(0.535, 0.215, 0.405, 0.545)},
                    {"name": "bottom bullet one", "kind": "text", "role": "audience", "font_px": 35, "bbox": [166, 1305, 2390, 1358]},
                    {"name": "bottom bullet two", "kind": "text", "role": "audience", "font_px": 30, "bbox": [166, 1362, 2390, 1414]},
                ],
            },
            indent=2,
        )
        + "\n"
    )

    manifest = {
        "schema": "the219_current_leading_xj_pp_auau_centrality_overlay_v1",
        "status": "detector_level_existing_artifact_diagnostic",
        "png": str(png),
        "script": str(Path(__file__).relative_to(REPO)),
        "script_sha256": sha256(Path(__file__)),
        "helper": str(HELPER_PATH.relative_to(REPO)),
        "helper_sha256": sha256(HELPER_PATH),
        "pointers": {
            key: {
                "path": str(path),
                "current_entry_id": pointers[key].get("current_entry_id"),
                "canonical_status": pointers[key].get("canonical_status"),
                "root": pointers[key]["root_paths"][0],
                "root_sha256": sha256(Path(pointers[key]["root_paths"][0])),
            }
            for key, path in load_helper().POINTERS.items()
        },
        "selection": {
            "pp_photon_pt_GeV": list(PP_PT_RANGE),
            "auau_photon_pt_GeV": list(AUAU_PT_RANGE),
            "photon_pt_binning_note": "nearest high-stat native bins; not an exact matched-pT comparison",
            "jet_radius": 0.4,
            "jet_pt_min_GeV": jet_pt_min,
            "recoil_delta_phi": ">=7pi/8",
            "jet_order": "one global highest-pT non-overlapping jet, then eta and delta-phi vetoes; no fallback",
            "jet_energy_treatment": "standard JetCalib applied to pp AntiKt_Tower_r04 and AuAu UE-subtracted AntiKt_Tower_r04_Sub1; no residual response correction or unfolding",
            "mc_reference": "reconstructed photon matched to prompt truth photon plus reconstructed leading recoil jet matched to truth leading recoil jet; h_JES3RecoTruthTagged family; not pure truth",
            "pp_photon": "PPG12 isolated+tight reconstructed photon",
            "auau_photon": "AuAu baseline BDT plus centrality-aware isolation",
            "shape_normalization": "each curve independently normalized on 0.5<xJgamma<2.0, then divided by bin width",
            "uncertainty_model": "vertical data/MC errors come from projected ROOT bin errors (Sumw2), divided by the fixed plotted normalization and bin width; normalization-denominator uncertainty and bin-to-bin covariance are not propagated; horizontal bars are bin half-widths; no uncertainty band and no systematic uncertainty are drawn",
            "display_range": list(DISPLAY_RANGE),
            "conservative_turn_on_xj": turn_on_xj,
        },
        "objects": {key: panel["objects"] for key, panel in payload.items()},
        "normalization_yields": {key: panel["normalization_yields"] for key, panel in payload.items()},
        "centrality_mapping": {
            "reference_screenshot": ["40-60% peripheral", "0-10% central mentioned by user"],
            "current_available": ["50-80% peripheral", "0-20% central", "20-50% mid-central"],
            "plotted": ["50-80% peripheral", "0-20% central"],
            "mapping_note": "Current bins do not reproduce the reference's 40-60% and 0-10% centralities exactly.",
        },
        "auau_data_trigger_directory": AUAU_DATA_TRIGGER_DIR,
        "claims_allowed": [
            "The plot compares current photon-triggered AuAu leading-jet detector-level shapes across centrality with current pp.",
            "The peripheral current-bin data retain a recoil-like peak while the central current-bin data contain a broad high-xJgamma tail.",
            "The blue MC curve is reconstructed and pair truth tagged; no pure-truth curve is drawn.",
        ],
        "claims_not_allowed": [
            "The central excess is a measured combinatorial fraction.",
            "The curves are background-subtracted, unfolded, residual-response corrected, or equipped with systematic covariance.",
            "The plot is an exact reproduction of the screenshot centralities or photon/trigger contract.",
            "The pp and AuAu curves use exactly matched photon-pT bins.",
            "The centrality contrast alone demonstrates jet quenching.",
        ],
        "sidecars": [str(csv_path), str(json_points), str(layout_path), str(speaker_script)],
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    speaker_script.write_text(
        f"# Speaker script — leading-jet centrality contrast, pTjet > {jet_pt_min:g} GeV\n\n"
        "This is the closest statistically useful version of the reference plot that we can make from the current registered outputs. "
        "Our plotted Au+Au bins are 50 to 80 percent and 0 to 20 percent, rather than 40 to 60 and 0 to 10. "
        "The data come from the photon-10-plus-MBD trigger directory used by the current AnalyzeRecoilJets configuration.\n\n"
        f"Every curve uses one highest-pT recoil jet per event and the exact internal-scan histogram for pTjet above {jet_pt_min:g} GeV. The pp JES3 object uses 16 to 18 GeV while the AuAu JES3 object uses 15 to 17 GeV. "
        "These are the nearest high-stat native bins, not an exact matched-pT comparison; reproducing the reference's 12 to 14 GeV AuAu bin would require event-level reprocessing. "
        "The shapes are normalized independently, so the comparison is about shape rather than absolute yield. "
        "The blue curve is the reconstructed pair truth-tagged subset: the reconstructed photon is matched to the prompt truth photon and the reconstructed leading recoil jet is matched to the truth leading recoil jet. It is not a pure-truth distribution. "
        "The 50-to-80 and 0-to-20 percent panels show the current centrality dependence directly.\n\n"
        "Any central mismatch must be treated as a detector-level diagnostic: before unfolding or interpreting quenching, central Au+Au needs an explicit combinatorial-jet treatment. "
        "The vertical error bars are statistical only; the horizontal bars show bin widths. No uncertainty band is drawn. "
        "Do not call this a measured background fraction or a final corrected result. It is a detector-level diagnostic from existing outputs.\n"
    )


def main() -> None:
    outputs = []
    for jet_pt_min in JET_PT_MIN_SCAN:
        paths = artifact_paths(jet_pt_min)
        pp, auau, pointers = build_payload(jet_pt_min)
        render(pp, auau, jet_pt_min, Path(paths["png"]))
        write_sidecars(pp, auau, pointers, jet_pt_min, paths)
        outputs.append(
            {
                "jet_pt_min_GeV": jet_pt_min,
                "png": str(paths["png"]),
                "manifest": str(paths["manifest"]),
                "speaker_script": str(paths["speaker_script"]),
            }
        )
    index_path = SCAN_DIR / "the219_leading_xj_jet_pt_min_scan_index.json"
    index_path.write_text(json.dumps({"thresholds": outputs}, indent=2) + "\n")
    print(json.dumps({"index": str(index_path), "thresholds": outputs}, indent=2))


if __name__ == "__main__":
    main()
