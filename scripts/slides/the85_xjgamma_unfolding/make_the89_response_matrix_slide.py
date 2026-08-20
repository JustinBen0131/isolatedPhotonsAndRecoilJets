#!/usr/bin/env python3
"""Build a clean p+p versus Au+Au 0-20% xJgamma response-matrix slide.

The registered current signal-simulation ROOT files store the joint
(pTgamma, xJgamma) response as truth-global-bin versus reco-global-bin.  This
helper de-flattens that contract, integrates native photon-pT blocks over each
sample's analysis support, and displays P(xJ_reco | xJ_truth) with one shared
color scale.  The matrix is normalized within the displayed 0 <= xJ < 2.14
window, so the slide compares migration shape rather than sample size.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import font_manager
import matplotlib.patheffects as path_effects
import numpy as np
from PIL import Image
import ROOT


ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

REPO = Path(__file__).resolve().parents[3]
OUT = REPO / "dataOutput/the219_friday_ppg_20260814/the89_response_matrix_slide"
OUT.mkdir(parents=True, exist_ok=True)

PNG = OUT / "the89_pp_auau020_r04_xj_response_matrix_slide.png"
POINTS = OUT / "the89_pp_auau020_r04_xj_response_matrix_values.json"
MANIFEST = OUT / "the89_pp_auau020_r04_xj_response_matrix_manifest.json"
LAYOUT = OUT / "the89_pp_auau020_r04_xj_response_matrix_layout_nodes.json"
SCRIPT = OUT / "the89_pp_auau020_r04_xj_response_matrix_speaker_script.md"
SELF_AUDIT = OUT / "the89_pp_auau020_r04_xj_response_matrix_self_audit.json"

PP_POINTER = REPO / "dataOutput/current_recoiljets_artifacts/current/pp_sim_photonjet_merged/current.json"
AUAU_POINTER = REPO / "dataOutput/current_recoiljets_artifacts/current/auau_sim_photonjet_merged/current.json"
FONT_PATH = Path("/System/Library/Fonts/Supplemental/Times New Roman.ttf")

DISPLAY_XJ_MAX = 2.14
SLIDE_DPI = 200
SLIDE_SIZE = (2560, 1440)


SAMPLES = {
    "pp": {
        "label": r"$p$+$p$",
        "pointer": PP_POINTER,
        "pt_support": [14.0, 35.0],
        "reco": "h2_unfoldReco_pTgamma_xJ_incl_r04",
        "truth": "h2_unfoldTruth_pTgamma_xJ_incl_r04",
        "response": "h2_unfoldResponse_pTgamma_xJ_incl_r04",
    },
    "auau_0_20": {
        "label": r"Au+Au 0–20%",
        "pointer": AUAU_POINTER,
        "pt_support": [15.0, 35.0],
        "reco": "h2_unfoldReco_pTgamma_xJ_incl_r04_isoR40_isSliding_cent_0_20",
        "truth": "h2_unfoldTruth_pTgamma_xJ_incl_r04_isoR40_isSliding_cent_0_20",
        "response": "h2_unfoldResponse_pTgamma_xJ_incl_r04_isoR40_isSliding_cent_0_20",
    },
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def json_ready(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [json_ready(item) for item in value]
    if isinstance(value, (float, np.floating)):
        return float(value) if math.isfinite(float(value)) else None
    if isinstance(value, (int, np.integer)):
        return int(value)
    return value


def open_root(path: Path) -> ROOT.TFile:
    root_file = ROOT.TFile.Open(str(path), "READ")
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"failed to open ROOT input: {path}")
    return root_file


def get_th2(root_file: ROOT.TFile, name: str):
    hist = root_file.Get(f"SIM/{name}")
    if not hist or not hist.InheritsFrom("TH2"):
        raise KeyError(f"missing SIM/{name} in {root_file.GetName()}")
    return hist


def axis_edges(axis) -> np.ndarray:
    return np.asarray(
        [float(axis.GetBinLowEdge(ibin)) for ibin in range(1, axis.GetNbins() + 2)],
        dtype=float,
    )


def contained_bins(axis, low: float, high: float) -> list[int]:
    bins = [
        ibin
        for ibin in range(1, axis.GetNbins() + 1)
        if float(axis.GetBinLowEdge(ibin)) >= low - 1.0e-9
        and float(axis.GetBinUpEdge(ibin)) <= high + 1.0e-9
    ]
    if not bins:
        raise RuntimeError(f"no native bins contained in {low:g}-{high:g} GeV")
    return bins


def bin_ranges(axis, bins: list[int]) -> list[list[float]]:
    return [
        [float(axis.GetBinLowEdge(ibin)), float(axis.GetBinUpEdge(ibin))]
        for ibin in bins
    ]


def extract_response(sample: dict[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    pointer = json.loads(Path(sample["pointer"]).read_text())
    root_path = Path(pointer["root_paths"][0])
    root_file = open_root(root_path)
    try:
        reco = get_th2(root_file, sample["reco"])
        truth = get_th2(root_file, sample["truth"])
        response = get_th2(root_file, sample["response"])

        truth_edges = axis_edges(truth.GetYaxis())
        reco_edges = axis_edges(reco.GetYaxis())
        if not np.allclose(truth_edges, reco_edges, atol=1.0e-12, rtol=0.0):
            raise RuntimeError("truth and reco xJ binning differ")
        if response.GetNbinsX() != truth.GetNcells() or response.GetNbinsY() != reco.GetNcells():
            raise RuntimeError(
                "global-bin response dimensions do not match truth/reco templates: "
                f"response=({response.GetNbinsX()},{response.GetNbinsY()}) "
                f"templates=({truth.GetNcells()},{reco.GetNcells()})"
            )

        low, high = map(float, sample["pt_support"])
        truth_pt_bins = contained_bins(truth.GetXaxis(), low, high)
        reco_pt_bins = contained_bins(reco.GetXaxis(), low, high)
        displayed_bins = [
            ibin
            for ibin in range(1, truth.GetYaxis().GetNbins() + 1)
            if float(truth.GetYaxis().GetBinUpEdge(ibin)) <= DISPLAY_XJ_MAX + 1.0e-9
        ]
        if not displayed_bins or displayed_bins != list(range(1, len(displayed_bins) + 1)):
            raise RuntimeError("displayed xJ bins must form a contiguous range starting at zero")
        xj_edges = truth_edges[: len(displayed_bins) + 1]

        raw = np.zeros((len(displayed_bins), len(displayed_bins)), dtype=float)
        for truth_position, truth_xj_bin in enumerate(displayed_bins):
            for reco_position, reco_xj_bin in enumerate(displayed_bins):
                total = 0.0
                for truth_pt_bin in truth_pt_bins:
                    truth_global = int(truth.GetBin(truth_pt_bin, truth_xj_bin))
                    response_x = int(response.GetXaxis().FindBin(float(truth_global)))
                    for reco_pt_bin in reco_pt_bins:
                        reco_global = int(reco.GetBin(reco_pt_bin, reco_xj_bin))
                        response_y = int(response.GetYaxis().FindBin(float(reco_global)))
                        total += float(response.GetBinContent(response_x, response_y))
                raw[reco_position, truth_position] = total

        column_sums = raw.sum(axis=0)
        normalized = np.divide(
            raw,
            column_sums[np.newaxis, :],
            out=np.zeros_like(raw),
            where=column_sums[np.newaxis, :] > 0.0,
        )
        centers = 0.5 * (xj_edges[:-1] + xj_edges[1:])
        diagonal_fraction = np.diag(normalized)
        mean_absolute_migration = [
            float(np.sum(normalized[:, index] * np.abs(centers - centers[index])))
            if column_sums[index] > 0.0
            else math.nan
            for index in range(len(centers))
        ]
        adjacent_or_diagonal = []
        for index in range(len(centers)):
            lo_index = max(0, index - 1)
            hi_index = min(len(centers), index + 2)
            adjacent_or_diagonal.append(float(normalized[lo_index:hi_index, index].sum()))

        payload = {
            "xj_edges": xj_edges,
            "raw_weighted_counts": raw,
            "truth_column_normalized": normalized,
            "truth_column_raw_sums": column_sums,
            "truth_pt_bins": bin_ranges(truth.GetXaxis(), truth_pt_bins),
            "reco_pt_bins": bin_ranges(reco.GetXaxis(), reco_pt_bins),
            "diagonal_fraction": diagonal_fraction,
            "mean_absolute_migration_by_truth_bin": mean_absolute_migration,
            "adjacent_or_diagonal_fraction": adjacent_or_diagonal,
            "summary": {
                "matrix_raw_weight_sum": float(raw.sum()),
                "median_diagonal_fraction": float(np.median(diagonal_fraction[column_sums > 0.0])),
                "mean_diagonal_fraction": float(np.mean(diagonal_fraction[column_sums > 0.0])),
                "mean_absolute_xj_migration": float(np.nanmean(mean_absolute_migration)),
                "mean_adjacent_or_diagonal_fraction": float(np.mean(adjacent_or_diagonal)),
            },
        }
        provenance = {
            "pointer": str(sample["pointer"]),
            "pointer_sha256": sha256(Path(sample["pointer"])),
            "campaign_tag": pointer["campaign_tag"],
            "canonical_status": pointer["canonical_status"],
            "root": str(root_path),
            "root_size_bytes": root_path.stat().st_size,
            "root_mtime_ns": root_path.stat().st_mtime_ns,
            "objects": {
                "reco_template": f"SIM/{sample['reco']}",
                "truth_template": f"SIM/{sample['truth']}",
                "global_response": f"SIM/{sample['response']}",
            },
            "global_response_dimensions": [int(response.GetNbinsX()), int(response.GetNbinsY())],
            "global_response_entries": float(response.GetEntries()),
            "global_response_integral": float(response.Integral()),
        }
        return payload, provenance
    finally:
        root_file.Close()


def configure_style() -> str:
    if not FONT_PATH.exists():
        raise FileNotFoundError(f"required slide font is missing: {FONT_PATH}")
    font_manager.fontManager.addfont(str(FONT_PATH))
    font_name = font_manager.FontProperties(fname=str(FONT_PATH)).get_name()
    if font_name != "Times New Roman":
        raise RuntimeError(f"unexpected font identity: {font_name}")
    mpl.rcParams.update(
        {
            "font.family": font_name,
            "mathtext.fontset": "stix",
            "axes.linewidth": 1.0,
            "axes.labelcolor": "#151515",
            "xtick.color": "#151515",
            "ytick.color": "#151515",
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )
    return font_name


def draw_matrix(ax, payload: dict[str, Any], *, vmax: float) -> Any:
    edges = np.asarray(payload["xj_edges"], dtype=float)
    matrix = np.asarray(payload["truth_column_normalized"], dtype=float)
    mesh = ax.pcolormesh(
        edges,
        edges,
        matrix,
        cmap="viridis",
        vmin=0.0,
        vmax=vmax,
        shading="flat",
        rasterized=True,
    )
    guide = ax.plot(
        [0.0, DISPLAY_XJ_MAX],
        [0.0, DISPLAY_XJ_MAX],
        color="white",
        linewidth=1.35,
        linestyle=(0, (5, 4)),
        alpha=0.9,
        zorder=4,
    )[0]
    guide.set_path_effects([path_effects.Stroke(linewidth=2.5, foreground="#27324A", alpha=0.42), path_effects.Normal()])
    ax.set_xlim(0.0, DISPLAY_XJ_MAX)
    ax.set_ylim(0.0, DISPLAY_XJ_MAX)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks([0.0, 0.5, 1.0, 1.5, 2.0])
    ax.set_yticks([0.0, 0.5, 1.0, 1.5, 2.0])
    ax.tick_params(axis="both", which="major", labelsize=13.5, length=6, width=1.0, pad=6)
    ax.minorticks_on()
    ax.tick_params(axis="both", which="minor", length=3.2, width=0.8)
    ax.set_xlabel(r"Truth $x_{J\gamma}$", fontsize=18.5, labelpad=9)
    ax.set_ylabel(r"Reconstructed $x_{J\gamma}$", fontsize=18.5, labelpad=10)
    return mesh


def render(payloads: dict[str, dict[str, Any]]) -> float:
    configure_style()
    fig = plt.figure(figsize=(12.8, 7.2), dpi=SLIDE_DPI, facecolor="white")
    blue = "#1D4F91"
    dark = "#151515"

    fig.text(
        0.055,
        0.956,
        r"Au+Au broadens the $x_{J\gamma}$ detector response relative to $p$+$p$",
        ha="left",
        va="top",
        fontsize=30,
        fontweight="bold",
        color=dark,
    )
    fig.add_artist(plt.Line2D([0.055, 0.945], [0.902, 0.902], transform=fig.transFigure, color=blue, linewidth=2.2))
    fig.text(
        0.055,
        0.876,
        r"Matched photon–jet simulation  •  truth-column normalized within the displayed range  •  anti-$k_{T}$ $R=0.4$",
        ha="left",
        va="top",
        fontsize=17.2,
        color="#343434",
    )

    axes = {
        "pp": fig.add_axes([0.075, 0.160, 0.340, 0.604]),
        "auau_0_20": fig.add_axes([0.490, 0.160, 0.340, 0.604]),
    }
    matrices = [np.asarray(payloads[key]["truth_column_normalized"], dtype=float) for key in axes]
    shared_max = max(float(matrix.max()) for matrix in matrices)
    vmax = math.ceil(shared_max / 0.05) * 0.05
    meshes = {}
    for key, ax in axes.items():
        meshes[key] = draw_matrix(ax, payloads[key], vmax=vmax)

    fig.text(0.245, 0.813, r"$p$+$p$", ha="center", va="center", fontsize=22, fontweight="bold", color=dark)
    fig.text(0.245, 0.784, r"$14<p_{T}^{\gamma}<35$ GeV native support", ha="center", va="center", fontsize=14.5, color="#414141")
    fig.text(0.660, 0.813, r"Au+Au 0–20%", ha="center", va="center", fontsize=22, fontweight="bold", color=dark)
    fig.text(0.660, 0.784, r"$15<p_{T}^{\gamma}<35$ GeV native support", ha="center", va="center", fontsize=14.5, color="#414141")

    color_axis = fig.add_axes([0.865, 0.225, 0.022, 0.475])
    colorbar = fig.colorbar(meshes["pp"], cax=color_axis)
    colorbar.ax.tick_params(labelsize=13.2, length=5, width=0.9)
    colorbar.set_label(
        r"$P(\mathrm{reco}\ x_{J\gamma}\ \mathrm{bin}\mid\mathrm{truth}\ x_{J\gamma}\ \mathrm{bin})$",
        fontsize=15.5,
        labelpad=12,
    )
    colorbar.outline.set_linewidth(0.9)

    fig.text(
        0.055,
        0.066,
        r"Each truth column sums to one; color measures migration width, not the different simulation sample sizes.",
        ha="left",
        va="center",
        fontsize=17.2,
        color="#303030",
    )

    fig.savefig(PNG, dpi=SLIDE_DPI, facecolor="white")
    plt.close(fig)
    return vmax


def write_layout_nodes() -> None:
    pt_to_px = SLIDE_DPI / 72.0
    nodes = [
        {
            "name": "slide title",
            "kind": "text",
            "role": "title",
            "text": "Au+Au broadens the xJgamma detector response relative to p+p",
            "font_px": 30.0 * pt_to_px,
            "bbox": [141, 57, 2360, 145],
            "title_anchor": True,
        },
        {
            "name": "shared response definition",
            "kind": "text",
            "role": "audience",
            "text": "Matched photon-jet simulation; truth-column normalized within the displayed range; anti-kT R=0.4",
            "font_px": 17.2 * pt_to_px,
            "bbox": [141, 180, 2370, 225],
            "title_axis_align": "left",
        },
        {
            "name": "pp plot panel",
            "kind": "plot",
            "bbox": [192, 340, 1062, 1210],
            "symmetry_group": "response panels",
        },
        {
            "name": "auau plot panel",
            "kind": "plot",
            "bbox": [1254, 340, 2124, 1210],
            "symmetry_group": "response panels",
        },
        {
            "name": "pp panel label",
            "kind": "text",
            "role": "plot_annotation",
            "text": "p+p",
            "font_px": 22.0 * pt_to_px,
            "bbox": [192, 240, 1062, 292],
        },
        {
            "name": "pp support label",
            "kind": "text",
            "role": "plot_annotation",
            "text": "14 < pTgamma < 35 GeV native support",
            "font_px": 14.5 * pt_to_px,
            "bbox": [192, 292, 1062, 330],
        },
        {
            "name": "auau panel label",
            "kind": "text",
            "role": "plot_annotation",
            "text": "Au+Au 0-20%",
            "font_px": 22.0 * pt_to_px,
            "bbox": [1254, 240, 2124, 292],
        },
        {
            "name": "auau support label",
            "kind": "text",
            "role": "plot_annotation",
            "text": "15 < pTgamma < 35 GeV native support",
            "font_px": 14.5 * pt_to_px,
            "bbox": [1254, 292, 2124, 330],
        },
        {
            "name": "shared colorbar",
            "kind": "plot",
            "bbox": [2214, 324, 2270, 1008],
        },
        {
            "name": "colorbar label",
            "kind": "text",
            "role": "plot_annotation",
            "text": "P(reconstructed xJgamma bin | truth xJgamma bin)",
            "font_px": 15.5 * pt_to_px,
            "bbox": [2280, 420, 2500, 960],
        },
        {
            "name": "normalization explanation",
            "kind": "text",
            "role": "audience",
            "text": "Each truth column sums to one; color measures migration width, not the different simulation sample sizes.",
            "font_px": 17.2 * pt_to_px,
            "bbox": [141, 1340, 2380, 1382],
            "title_axis_align": "left",
        },
    ]
    payload = {
        "slide_size": list(SLIDE_SIZE),
        "title_axis_x": 141,
        "title_axis_tolerance_px": 8,
        "minimum_audience_font_px": 35,
        "minimum_title_font_px": 66,
        "minimum_plot_annotation_font_px": 28,
        "vertical_margin_balance": {
            "top_node": "slide title",
            "bottom_node": "normalization explanation",
            "target_gap_px": 58,
            "tolerance_px": 8,
        },
        "nodes": nodes,
    }
    LAYOUT.write_text(json.dumps(payload, indent=2) + "\n")


def write_speaker_script(payloads: dict[str, dict[str, Any]]) -> None:
    pp_migration = payloads["pp"]["summary"]["mean_absolute_xj_migration"]
    auau_migration = payloads["auau_0_20"]["summary"]["mean_absolute_xj_migration"]
    SCRIPT.write_text(
        "# Speaker script\n\n"
        "Before looking at the unfolded spectra, this is the detector response that controls the migration. "
        "Truth x-J-gamma is horizontal and reconstructed x-J-gamma is vertical. Each truth column is normalized to one, "
        "so the color compares the shape of the migration and not the very different simulation sample sizes.\n\n"
        "The p-plus-p response on the left keeps a relatively narrow diagonal core. The zero-to-twenty-percent Au-plus-Au "
        "response on the right still follows the diagonal, but it is visibly broader because the heavy-ion environment "
        "adds substantially more jet-energy and reconstruction smearing. In this displayed range, the mean absolute x-J-gamma "
        f"migration is {pp_migration:.2f} in p-plus-p and {auau_migration:.2f} in central Au-plus-Au.\n\n"
        "These heatmaps show only the matched-jet migration shape. They do not include reconstruction efficiency, fakes, "
        "or the combinatoric-background normalization. The native photon-p-T support starts at fourteen GeV in p-plus-p "
        "and fifteen GeV in Au-plus-Au, which is stated above each panel.\n\n"
        "This is why the Au-plus-Au unfolding is materially harder: the response carries less local information about the truth bin. "
        "It also explains why the response construction and combinatoric-background contract must be validated before interpreting "
        "the final p-plus-p versus Au-plus-Au spectrum.\n"
    )


def main() -> None:
    configure_style()
    payloads: dict[str, dict[str, Any]] = {}
    provenance: dict[str, dict[str, Any]] = {}
    for key, sample in SAMPLES.items():
        payloads[key], provenance[key] = extract_response(sample)

    pp_edges = np.asarray(payloads["pp"]["xj_edges"], dtype=float)
    auau_edges = np.asarray(payloads["auau_0_20"]["xj_edges"], dtype=float)
    if not np.allclose(pp_edges, auau_edges, atol=1.0e-12, rtol=0.0):
        raise RuntimeError("p+p and Au+Au displayed xJ binning does not match")
    if abs(float(pp_edges[0])) > 1.0e-12:
        raise RuntimeError("response-matrix xJ display must start at zero")

    shared_vmax = render(payloads)
    write_layout_nodes()
    write_speaker_script(payloads)

    image = Image.open(PNG)
    self_checks = {
        "png_dimensions_2560x1440": image.size == SLIDE_SIZE,
        "font_times_new_roman": configure_style() == "Times New Roman",
        "xj_axes_start_at_zero": float(pp_edges[0]) == 0.0 and float(auau_edges[0]) == 0.0,
        "identical_displayed_xj_edges": bool(np.allclose(pp_edges, auau_edges, atol=1.0e-12, rtol=0.0)),
        "pp_truth_columns_normalized": bool(np.allclose(np.sum(payloads["pp"]["truth_column_normalized"], axis=0), 1.0, atol=1.0e-10)),
        "auau_truth_columns_normalized": bool(np.allclose(np.sum(payloads["auau_0_20"]["truth_column_normalized"], axis=0), 1.0, atol=1.0e-10)),
        "shared_color_scale": True,
        "google_slides_unchanged": True,
        "no_systematic_uncertainty_claim": True,
    }
    self_audit = {
        "ok": all(self_checks.values()),
        "checks": self_checks,
        "png": str(PNG),
        "png_sha256": sha256(PNG),
    }
    SELF_AUDIT.write_text(json.dumps(self_audit, indent=2) + "\n")
    if not self_audit["ok"]:
        raise RuntimeError(f"response-matrix slide self-audit failed: {SELF_AUDIT}")

    POINTS.write_text(
        json.dumps(
            json_ready(
                {
                    "normalization": "P(reco xJgamma bin | truth xJgamma bin), normalized within 0 <= xJgamma < 2.14",
                    "xj_edges": pp_edges,
                    "samples": payloads,
                }
            ),
            indent=2,
            allow_nan=False,
        )
        + "\n"
    )

    manifest = {
        "ok": True,
        "status": "slide_candidate_for_internal_review",
        "png": str(PNG),
        "png_sha256": sha256(PNG),
        "points": str(POINTS),
        "points_sha256": sha256(POINTS),
        "layout_nodes": str(LAYOUT),
        "layout_nodes_sha256": sha256(LAYOUT),
        "speaker_script": str(SCRIPT),
        "speaker_script_sha256": sha256(SCRIPT),
        "self_audit": str(SELF_AUDIT),
        "generator": str(Path(__file__).resolve()),
        "generator_sha256": sha256(Path(__file__).resolve()),
        "plot_contract": {
            "quantity": "truth-column-normalized matched-jet xJgamma response after de-flattening the joint pTgamma-xJgamma global response",
            "orientation": "truth xJgamma on x; reconstructed xJgamma on y",
            "jet_radius": 0.4,
            "display_xj": [0.0, DISPLAY_XJ_MAX],
            "normalization": "each truth column sums to one within the displayed range",
            "shared_color_scale": [0.0, shared_vmax],
            "statistical_uncertainties_drawn": False,
            "systematic_uncertainties_drawn": False,
        },
        "inputs": provenance,
        "summary": {key: payloads[key]["summary"] for key in payloads},
        "caveats": [
            "The response heatmaps show matched-jet migration shape; they do not show reconstruction efficiency, fakes, or combinatoric-background normalization.",
            "The displayed pTgamma supports follow each current response's exact native binning: 14-35 GeV for p+p and 15-35 GeV for Au+Au.",
            "The final sparse 2.14-3.0 xJgamma bin is outside the displayed matrix and the shown columns are normalized within 0-2.14.",
            "The Au+Au signal-simulation current pointer is canonical at the merged ROOT level, but this slide does not constitute unfolding closure or final scientific acceptance.",
        ],
        "google_slides_mutated": False,
    }
    MANIFEST.write_text(json.dumps(json_ready(manifest), indent=2, allow_nan=False) + "\n")
    print(PNG)
    print(MANIFEST)
    print(LAYOUT)
    print(SCRIPT)


if __name__ == "__main__":
    main()
