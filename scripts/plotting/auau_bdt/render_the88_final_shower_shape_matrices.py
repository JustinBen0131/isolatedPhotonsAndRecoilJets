#!/usr/bin/env python3
"""Render the current pp/Au+Au shower-shape matrices for the IAN.

All Au+Au lanes resolve through the registered final THE-88 pointers.  The
exact-zero boundary bin and histogram under/overflow are not drawn, but they
remain in the normalization denominator.  The visible curves therefore show
the fraction of all finite candidates in each continuous in-range bin rather
than a unit-area conditional interior shape.  This avoids visually amplifying
the sparse high-width component in the raw data while preserving it without
an arbitrary width cut.  The renderer also corrects the statistical errors of
the pre/tight histograms for their four identical internal-isolation-view
fills.  It deliberately leaves et2--et4 as explicit unavailable cells because
those histogram families were not persisted in the completed production.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = None
OPEN_FILES: dict[Path, object] = {}


CURRENT_ROOT = Path(
    "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/"
    "current_recoiljets_artifacts/current"
)
AUAU_DATA_POINTER = CURRENT_ROOT / "auau_data_merged/current.json"
AUAU_SIGNAL_POINTER = CURRENT_ROOT / "auau_sim_photonjet_merged/current.json"
AUAU_INCLUSIVE_POINTER = CURRENT_ROOT / "auau_sim_inclusivejet_merged/current.json"
PP_PHOTONJET_POINTER = Path(
    "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/current_recoiljets_artifacts/"
    "current/pp_sim_photonjet_merged/current.json"
)
PP_INCLUSIVE_POINTER = Path(
    "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/current_recoiljets_artifacts/"
    "current/pp_sim_inclusivejet_merged/current.json"
)
PP_INCLUSIVE_COMPAT_ROOT = Path(
    "/Users/patsfan753/Desktop/ThesisAnalysis/InputFiles/"
    "the97_ppg12_si_contract_restore_full_20260715_1420/final_merged_roots/"
    "final_combined_canonical_20260716/"
    "RecoilJets_jet8plus12plus20plus30plus40_correctedSI_validDI_period_combined_MERGED.root"
)
DATA_DIRECTORY = "photon_12_plus_MBD_NS_geq_2_vtx_lt_150"

PT_SLICES = ("15_17", "17_19", "19_21", "21_23", "23_26", "26_35")
CENTRALITIES = (
    ("0_20", "0--20%", "0_20"),
    ("20_50", "20--50%", "20_50"),
    ("50_80", "50--80%", "50_80"),
)
STAGES = (("inclusive", "Before preselection"), ("pre", "After preselection"), ("tight", "After tight ID"))
AUAU_STAGE_FILL_MULTIPLICITY = {"inclusive": 1, "pre": 4, "tight": 4}
VARIABLES = {
    "weta": r"$w_{\eta}^{\mathrm{COGX}}$",
    "wphi": r"$w_{\phi}^{\mathrm{COGX}}$",
    "e11e33": r"$E_{11}/E_{33}$",
    "e32e35": r"$E_{32}/E_{35}$",
    "et1": r"$e_{T,1}$",
}
UNAVAILABLE = {
    "et2": r"$e_{T,2}$",
    "et3": r"$e_{T,3}$",
    "et4": r"$e_{T,4}$",
}
PP_H2D_VARIABLES = {
    "weta": "weta_cogx",
    "wphi": "wphi_cogx",
    "e11e33": "e11_to_e33",
    "e32e35": "e32_to_e35",
    "et1": "et1",
    "et2": "et2",
    "et3": "et3",
    "et4": "et4",
}
PP_STAGES = (("0", "Before preselection"), ("1", "After preselection"), ("2", "After tight ID"))
COLORS = {"Data": "#161B27", "Signal MC": "#C2352B", "Inclusive MC": "#2D64D7"}


@dataclass
class Curve:
    x: np.ndarray
    y: np.ndarray
    e: np.ndarray
    raw_entries: float
    effective_entries: float
    fill_multiplicity: int
    total_weight: float
    retained_weight: float
    retained_fraction: float
    underflow_weight: float
    first_bin_weight: float
    overflow_weight: float
    first_bin_low: float
    first_bin_high: float
    object_names: list[str]


def require_root():
    global ROOT
    if ROOT is None:
        import ROOT as root

        root.gROOT.SetBatch(True)
        root.TH1.AddDirectory(False)
        ROOT = root
    return ROOT


def open_root(path: Path):
    if path in OPEN_FILES:
        return OPEN_FILES[path]
    root = require_root()
    f = root.TFile.Open(str(path), "READ")
    if not f or f.IsZombie():
        raise OSError(f"Cannot open ROOT input: {path}")
    # The completed data file has tens of thousands of top-level keys.  This
    # local ROOT build can terminate while destructing that file after a
    # close(), so retain one read-only handle per finite render process.  The
    # process exits immediately after the render and releases the descriptors.
    OPEN_FILES[path] = f
    return f


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def summed_histogram(file_handle, object_names: list[str], label: str):
    total = None
    for name in object_names:
        obj = file_handle.Get(name)
        if not obj or not obj.InheritsFrom("TH1"):
            raise KeyError(f"Missing expected histogram {name}")
        if total is None:
            total = obj.Clone(f"{label}_sum")
            total.SetDirectory(0)
        else:
            total.Add(obj)
    return total


def curve_from_histogram(
    hist,
    object_names: list[str],
    *,
    fill_multiplicity: int = 1,
) -> Curve:
    """Return the visible in-range contribution to the full candidate sample.

    The first ordinary bin is the exact/near-zero boundary bin for every
    persisted shower-shape family used here.  It is a distinct unresolved or
    degenerate component in current Au+Au data, not part of the continuous
    shape compared against embedding.  Underflow, that first bin, and overflow
    are not drawn, but they remain in the normalization denominator.  The area
    of the visible curve is therefore the retained fraction, and no discarded
    component is promoted by renormalization.

    The completed Au+Au production evaluated four internal isolation views.
    The canonical inclusive histograms were guarded to one view, while the
    pre/tight shower-shape histograms were filled identically in all four.
    Repetition leaves the normalized central values unchanged but makes naive
    ROOT errors too small by sqrt(4); ``fill_multiplicity`` restores the
    unique-candidate statistical scale.
    """
    if fill_multiplicity < 1:
        raise ValueError(f"Invalid fill multiplicity: {fill_multiplicity}")
    axis = hist.GetXaxis()
    count = hist.GetNbinsX()
    # Bin 1 is the zero-boundary component; bins 2..N are the resolved shape.
    x = np.asarray([axis.GetBinCenter(i) for i in range(2, count + 1)], dtype=float)
    y = np.asarray([hist.GetBinContent(i) for i in range(2, count + 1)], dtype=float)
    e = np.asarray([hist.GetBinError(i) for i in range(2, count + 1)], dtype=float)
    underflow = float(hist.GetBinContent(0))
    first_bin = float(hist.GetBinContent(1))
    overflow = float(hist.GetBinContent(count + 1))
    retained = float(np.sum(y))
    total = retained + underflow + first_bin + overflow
    if retained <= 0.0:
        raise ValueError(f"Nonpositive integral for {object_names[0]}")
    if total <= 0.0:
        raise ValueError(f"Nonpositive full finite-candidate weight for {object_names[0]}")
    raw_entries = float(hist.GetEntries())
    return Curve(
        x=x,
        y=y / total,
        e=(e * np.sqrt(float(fill_multiplicity))) / total,
        raw_entries=raw_entries,
        effective_entries=raw_entries / float(fill_multiplicity),
        fill_multiplicity=fill_multiplicity,
        total_weight=total,
        retained_weight=retained,
        retained_fraction=(retained / total) if total > 0.0 else 0.0,
        underflow_weight=underflow,
        first_bin_weight=first_bin,
        overflow_weight=overflow,
        first_bin_low=float(axis.GetBinLowEdge(1)),
        first_bin_high=float(axis.GetBinUpEdge(1)),
        object_names=object_names,
    )


def data_names(variable: str, stage: str, centrality: str) -> list[str]:
    return [
        f"{DATA_DIRECTORY}/h_ss_{variable}_{stage}_pT_{pt}_cent_{centrality}"
        for pt in PT_SLICES
    ]


def sim_names(variable: str, stage: str, category: str, centrality: str) -> list[str]:
    return [
        f"SIM/h_ss_{variable}_{stage}_{category}_pT_{pt}_cent_{cent}"
        for pt in PT_SLICES
        for cent in (centrality,)
    ]


def load_curve(
    root_path: Path,
    names: list[str],
    label: str,
    *,
    fill_multiplicity: int,
) -> Curve:
    handle = open_root(root_path)
    histogram = summed_histogram(handle, names, label)
    return curve_from_histogram(
        histogram,
        names,
        fill_multiplicity=fill_multiplicity,
    )


def root_from_pointer(pointer: Path) -> tuple[Path, dict]:
    payload = json.loads(pointer.read_text(encoding="utf-8"))
    roots = [Path(item) for item in payload.get("root_paths", [])]
    if len(roots) != 1 or not roots[0].is_file():
        raise FileNotFoundError(f"Expected one readable ROOT path in {pointer}: {roots}")
    return roots[0], payload


def root_has_object(root_path: Path, object_name: str) -> bool:
    handle = open_root(root_path)
    obj = handle.Get(object_name)
    return bool(obj and obj.InheritsFrom("TH1"))


def pp_names(variable: str, stage: str) -> list[str]:
    return [f"SIM/h2d_{PP_H2D_VARIABLES[variable]}_eta0_pt1535_cut{stage}"]


def load_pp_curve(root_path: Path, variable: str, stage: str, label: str) -> Curve:
    names = pp_names(variable, stage)
    handle = open_root(root_path)
    source = handle.Get(names[0])
    if not source or not source.InheritsFrom("TH2"):
        raise KeyError(f"Missing expected p+p 2D histogram {names[0]}")
    projection = source.ProjectionX(f"{label}_{variable}_{stage}_projection")
    projection.SetDirectory(0)
    return curve_from_histogram(projection, names)


def curve_record(curve: Curve, **identity) -> dict:
    return {
        **identity,
        "raw_entries": curve.raw_entries,
        "effective_entries": curve.effective_entries,
        "fill_multiplicity": curve.fill_multiplicity,
        "total_weight_including_boundaries": curve.total_weight,
        "retained_weight": curve.retained_weight,
        "retained_fraction": curve.retained_fraction,
        "excluded_fraction": 1.0 - curve.retained_fraction,
        "underflow_weight": curve.underflow_weight,
        "zero_boundary_bin_weight": curve.first_bin_weight,
        "overflow_weight": curve.overflow_weight,
        "zero_boundary_bin_edges": [curve.first_bin_low, curve.first_bin_high],
        "objects": curve.object_names,
    }


def style_axes(ax):
    ax.grid(True, axis="y", color="#E6E9EF", linewidth=0.55)
    ax.tick_params(direction="in", top=True, right=True, labelsize=7, pad=1)
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)
        spine.set_color("#1F2937")


def draw_curve(ax, curve: Curve, label: str):
    color = COLORS[label]
    if label == "Data":
        ax.errorbar(curve.x, curve.y, yerr=curve.e, fmt="o", color=color, ms=3.0, mfc=color, mec="white", mew=0.45, lw=0.65, capsize=1.0, label=label, zorder=4)
    else:
        ax.step(curve.x, curve.y, where="mid", color=color, lw=1.65, label=label, zorder=2)
        ax.errorbar(curve.x, curve.y, yerr=curve.e, fmt="none", ecolor=color, elinewidth=0.35, alpha=0.42, zorder=1)


def draw_pp_row(ax, variable: str, stage: str, photonjet_root: Path, inclusive_root: Path, manifest: dict):
    ax.set_facecolor("#F3F7FC")
    curves = {
        "Signal MC": load_pp_curve(photonjet_root, variable, stage, "pp_photonjet"),
        "Inclusive MC": load_pp_curve(inclusive_root, variable, stage, "pp_inclusive"),
    }
    ymax = 0.0
    for label, curve in curves.items():
        draw_curve(ax, curve, label)
        ymax = max(ymax, float(np.max(curve.y + curve.e)))
        manifest["curves"].append(
            curve_record(
                curve,
                variable=variable,
                centrality="pp_reference",
                stage=f"cut{stage}",
                lane="photonjet" if label == "Signal MC" else "inclusivejet",
            )
        )
    ax.set_ylim(0.0, max(0.025, 1.16 * ymax))
    style_axes(ax)
    ax.set_ylabel("p+p\nreference", fontsize=8.4, labelpad=7, color="#183B5B")
    ax.text(0.035, 0.94, r"$\it{\bf{sPHENIX}}$ Internal", transform=ax.transAxes, ha="left", va="top", fontsize=7.5)
    ax.text(0.035, 0.80, r"$p$+$p$ $\sqrt{s}=200$ GeV" + "\n" + r"$15<E_T^\gamma<35$ GeV, $|\eta|<0.7$", transform=ax.transAxes, ha="left", va="top", fontsize=6.5, color="#183B5B")


def plot_matrix(variable: str, outdir: Path, manifest: dict, data_root: Path, signal_root: Path, inclusive_root: Path, pp_photonjet_root: Path, pp_inclusive_root: Path):
    figure, axes = plt.subplots(4, 3, figsize=(7.35, 8.25), constrained_layout=False)
    figure.subplots_adjust(left=0.115, right=0.982, top=0.865, bottom=0.085, wspace=0.16, hspace=0.18)
    figure.patch.set_facecolor("white")
    for col, (stage, stage_label) in enumerate(PP_STAGES):
        ax = axes[0, col]
        draw_pp_row(ax, variable, stage, pp_photonjet_root, pp_inclusive_root, manifest)
        ax.set_title(stage_label, fontsize=10, fontweight="bold", color="#173B63", pad=5)
        if col != 0:
            ax.set_ylabel("")
        ax.tick_params(labelbottom=False)
    for auau_row, (cent_token, cent_label, sim_cent) in enumerate(CENTRALITIES, start=1):
        for col, (stage, stage_label) in enumerate(STAGES):
            ax = axes[auau_row, col]
            fill_multiplicity = AUAU_STAGE_FILL_MULTIPLICITY[stage]
            curves = {
                "Data": load_curve(
                    data_root,
                    data_names(variable, stage, cent_token),
                    f"data_{variable}_{cent_token}_{stage}",
                    fill_multiplicity=fill_multiplicity,
                ),
                "Signal MC": load_curve(
                    signal_root,
                    sim_names(variable, stage, "sig", sim_cent),
                    f"signal_{variable}_{cent_token}_{stage}",
                    fill_multiplicity=fill_multiplicity,
                ),
                "Inclusive MC": load_curve(
                    inclusive_root,
                    sim_names(variable, stage, "bkg", sim_cent),
                    f"inclusive_{variable}_{cent_token}_{stage}",
                    fill_multiplicity=fill_multiplicity,
                ),
            }
            ymax = 0.0
            for label, curve in curves.items():
                draw_curve(ax, curve, label)
                ymax = max(ymax, float(np.max(curve.y + curve.e)))
                manifest["curves"].append(
                    curve_record(
                        curve,
                        variable=variable,
                        centrality=cent_token,
                        stage=stage,
                        lane=label,
                    )
                )
            ax.set_ylim(0.0, max(0.025, 1.16 * ymax))
            style_axes(ax)
            if col == 0:
                ax.set_ylabel(f"Au+Au\n{cent_label}", fontsize=8.4, labelpad=7)
                ax.text(0.035, 0.95, r"$\it{\bf{sPHENIX}}$ Internal", transform=ax.transAxes, ha="left", va="top", fontsize=7.6)
                ax.text(0.035, 0.80, r"Au+Au $\sqrt{s_{NN}}=200$ GeV" + "\n" + r"$15<E_T^\gamma<35$ GeV, $|\eta|<0.7$", transform=ax.transAxes, ha="left", va="top", fontsize=6.5)
            data_coverage = curves["Data"].retained_fraction
            if data_coverage < 0.995:
                ax.text(
                    0.965,
                    0.955,
                    f"displayed data: {100.0 * data_coverage:.0f}%",
                    transform=ax.transAxes,
                    ha="right",
                    va="top",
                    fontsize=5.4,
                    color="#4B5563",
                    bbox={"boxstyle": "round,pad=0.18", "facecolor": "white", "edgecolor": "#CBD5E1", "linewidth": 0.45, "alpha": 0.92},
                    zorder=6,
                )
            if auau_row == 3:
                ax.set_xlabel(VARIABLES[variable], fontsize=9)
            else:
                ax.tick_params(labelbottom=False)
    figure.text(0.115, 0.975, f"p+p and Au+Au shower-shape comparison: {VARIABLES[variable]}", ha="left", va="top", fontsize=13.0, fontweight="bold", color="#111827")
    figure.text(0.115, 0.946, "Shaded p+p row: matched photon+jet and inclusive-jet simulation.", ha="left", va="top", fontsize=7.25, color="#52657A")
    figure.text(0.115, 0.930, "Au+Au rows: current data and matched embedding; fractions of all finite candidates.", ha="left", va="top", fontsize=7.25, color="#52657A")
    handles = [
        plt.Line2D([0], [0], color=COLORS["Signal MC"], lw=2.0, label="Photon simulation"),
        plt.Line2D([0], [0], color=COLORS["Inclusive MC"], lw=2.0, label="Inclusive-jet simulation"),
        plt.Line2D([0], [0], color=COLORS["Data"], marker="o", lw=0, ms=5, mfc=COLORS["Data"], mec="white", label="Au+Au data"),
    ]
    figure.legend(handles=handles, loc="upper right", bbox_to_anchor=(0.982, 0.914), frameon=False, ncol=3, fontsize=6.65, handlelength=1.45, columnspacing=0.7)
    figure.text(
        0.115,
        0.030,
        'Exact-zero and under/overflow components remain in the denominator; "displayed data" is the visible fraction.',
        ha="left",
        va="bottom",
        fontsize=6.15,
        color="#52657A",
    )
    output = outdir / f"pp_auau_{variable}_stage_matrix_the88_final.png"
    figure.savefig(output, dpi=230)
    plt.close(figure)
    manifest["assets"][output.name] = sha256(output)


def plot_unavailable(variable: str, outdir: Path, manifest: dict, pp_photonjet_root: Path, pp_inclusive_root: Path):
    figure, axes = plt.subplots(4, 3, figsize=(7.35, 8.25), constrained_layout=False)
    figure.subplots_adjust(left=0.115, right=0.982, top=0.865, bottom=0.085, wspace=0.16, hspace=0.18)
    figure.patch.set_facecolor("white")
    for col, (stage, stage_label) in enumerate(PP_STAGES):
        ax = axes[0, col]
        draw_pp_row(ax, variable, stage, pp_photonjet_root, pp_inclusive_root, manifest)
        ax.set_title(stage_label, fontsize=10, fontweight="bold", color="#173B63", pad=5)
        if col != 0:
            ax.set_ylabel("")
        ax.tick_params(labelbottom=False)
    for auau_row, (_, cent_label, _) in enumerate(CENTRALITIES, start=1):
        for col, (_, stage_label) in enumerate(STAGES):
            ax = axes[auau_row, col]
            ax.set_facecolor("#FFF9ED")
            ax.set_xticks([]); ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_color("#C47A14"); spine.set_linewidth(0.95)
            if col == 0:
                ax.set_ylabel(f"Au+Au\n{cent_label}", fontsize=8.4, labelpad=7)
            ax.text(0.5, 0.60, "SOURCE HISTOGRAM\nNOT PERSISTED", transform=ax.transAxes, ha="center", va="center", fontsize=9, fontweight="bold", color="#A65C00", linespacing=1.18)
            ax.text(0.5, 0.30, "The current Au+Au products contain no\ncommon-contract histogram for this cell.", transform=ax.transAxes, ha="center", va="center", fontsize=6.7, color="#6B4A1D")
    figure.text(0.115, 0.975, f"p+p and Au+Au shower-shape availability: {UNAVAILABLE[variable]}", ha="left", va="top", fontsize=13.0, fontweight="bold", color="#111827")
    figure.text(0.115, 0.946, "Shaded p+p row: matched photon+jet and inclusive-jet simulation.", ha="left", va="top", fontsize=7.25, color="#52657A")
    figure.text(0.115, 0.930, "No Au+Au curve is inferred from another observable or from simulation.", ha="left", va="top", fontsize=7.25, color="#52657A")
    handles = [
        plt.Line2D([0], [0], color=COLORS["Signal MC"], lw=2.0, label="Photon simulation"),
        plt.Line2D([0], [0], color=COLORS["Inclusive MC"], lw=2.0, label="Inclusive-jet simulation"),
    ]
    figure.legend(handles=handles, loc="upper right", bbox_to_anchor=(0.982, 0.914), frameon=False, ncol=2, fontsize=6.65, handlelength=1.45, columnspacing=0.7)
    output = outdir / f"pp_auau_{variable}_stage_matrix_the88_final.png"
    figure.savefig(output, dpi=230)
    plt.close(figure)
    manifest["assets"][output.name] = sha256(output)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    data_root, data_meta = root_from_pointer(AUAU_DATA_POINTER)
    signal_root, signal_meta = root_from_pointer(AUAU_SIGNAL_POINTER)
    inclusive_root, inclusive_meta = root_from_pointer(AUAU_INCLUSIVE_POINTER)
    pp_photonjet_root, pp_photonjet_meta = root_from_pointer(PP_PHOTONJET_POINTER)
    pp_inclusive_pointer_root, pp_inclusive_meta = root_from_pointer(PP_INCLUSIVE_POINTER)
    pp_inclusive_root = pp_inclusive_pointer_root
    pp_inclusive_resolution = "current pointer"
    required_pp_shape = "SIM/h2d_weta_cogx_eta0_pt1535_cut0"
    if not root_has_object(pp_inclusive_root, required_pp_shape):
        if not PP_INCLUSIVE_COMPAT_ROOT.is_file() or not root_has_object(PP_INCLUSIVE_COMPAT_ROOT, required_pp_shape):
            raise KeyError(
                f"Current p+p inclusive ROOT lacks {required_pp_shape}, and no compatible fallback is available"
            )
        pp_inclusive_root = PP_INCLUSIVE_COMPAT_ROOT
        pp_inclusive_resolution = (
            "last compatible registered shower-shape source; the newer ownership-gated current ROOT "
            "does not persist SIM/h2d_<variable>_eta0_pt1535_cut{0,1,2}"
        )
    for path in (data_root, signal_root, inclusive_root, pp_photonjet_root, pp_inclusive_root):
        if not path.is_file():
            raise FileNotFoundError(path)
    args.outdir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "serif", "font.serif": ["Times New Roman", "Times", "DejaVu Serif"], "axes.facecolor": "white", "figure.facecolor": "white"})
    manifest = {
        "schema": "IAN_AUAU_SHOWER_SHAPE_MATRICES_V4",
        "status": "current_final_the88_triplet",
        "renderer": {
            "path": str(Path(__file__).resolve()),
            "sha256": sha256(Path(__file__).resolve()),
        },
        "data": {"pointer": str(AUAU_DATA_POINTER), "current_entry_id": data_meta.get("current_entry_id"), "path": str(data_root), "sha256": sha256(data_root), "coverage": "50460/50460 jobs", "selection_directory": DATA_DIRECTORY, "role": "completed corrected-centrality THE-88 data merge"},
        "signal": {"pointer": str(AUAU_SIGNAL_POINTER), "current_entry_id": signal_meta.get("current_entry_id"), "path": str(signal_root), "sha256": sha256(signal_root), "role": "completed THE-88 embedded Photon12+20 signal"},
        "inclusive": {"pointer": str(AUAU_INCLUSIVE_POINTER), "current_entry_id": inclusive_meta.get("current_entry_id"), "path": str(inclusive_root), "sha256": sha256(inclusive_root), "role": "completed THE-88 embedded Jet12+20+30+40 background"},
        "pp_reference": {
            "photonjet_pointer": str(PP_PHOTONJET_POINTER),
            "photonjet_root": str(pp_photonjet_root),
            "photonjet_root_sha256": sha256(pp_photonjet_root),
            "photonjet_current_entry_id": pp_photonjet_meta.get("current_entry_id"),
            "inclusive_pointer": str(PP_INCLUSIVE_POINTER),
            "inclusive_pointer_current_root": str(pp_inclusive_pointer_root),
            "inclusive_pointer_current_entry_id": pp_inclusive_meta.get("current_entry_id"),
            "inclusive_root": str(pp_inclusive_root),
            "inclusive_root_sha256": sha256(pp_inclusive_root),
            "inclusive_resolution": pp_inclusive_resolution,
            "selection": "SIM/h2d_<variable>_eta0_pt1535_cut0, cut1, cut2; ProjectionX over the stored isolation axis",
            "role": "last compatible p+p shower-shape simulation reference; no p+p data curve is included",
        },
        "contract": {
            "photon_et_gev": "15 < E_T^gamma < 35",
            "photon_abs_eta": "< 0.7",
            "centrality": ["0-20%", "20-50%", "50-80%"],
            "stages": ["inclusive", "pre", "tight"],
            "normalization": "visible bins divided by the full finite-candidate histogram weight, including underflow, the exact-zero first ordinary bin, and overflow",
            "interpretation": "fraction of all finite candidates in each visible continuous in-range shower-shape bin; the visible curve area equals the retained fraction",
        },
        "quality_contract": {
            "boundary_bin_treatment": "The first ordinary bin is not drawn in any lane, but remains in the normalization denominator. Its edges and weight are recorded per curve; this is a symmetric display-domain definition, not point-by-point editing.",
            "overflow_treatment": "Histogram underflow and overflow are not drawn, but remain in the normalization denominator and are recorded per curve.",
            "auau_stage_fill_multiplicity": AUAU_STAGE_FILL_MULTIPLICITY,
            "stage_error_correction": "Pre and tight histograms were filled identically for four internal isolation views. Central values are unchanged; normalized ROOT errors are multiplied by sqrt(4), and raw entries are divided by 4 for effective-entry reporting.",
            "source_code_evidence": [
                "src_AuAu/RecoilJets_AuAu.cc: processCandidates loops over m_internalIsoViews",
                "src_AuAu/RecoilJets_AuAu.cc: inclusive fillSSSpectra is guarded by fillCanonicalThisView",
                "src_AuAu/RecoilJets_AuAu.cc: fillSSPPG12 pre/tight block is not guarded by fillCanonicalThisView",
                "scripts/sdcc/runtime/condor/RecoilJets_Condor_submit.sh: RJ_INTERNAL_ISO_VIEWS contains four views",
            ],
        },
        "variables_populated": list(VARIABLES),
        "variables_unavailable": list(UNAVAILABLE),
        "weta_raw_tail_audit": {
            "scope": "data, photon_12_plus_MBD_NS_geq_2_vtx_lt_150, 0-20%, 21<E_T^gamma<23 GeV, 0.9<=w_eta^COGX<1.0",
            "raw_entries": 5270,
            "effective_after_preselection": 94,
            "effective_after_tight_id": 80,
            "interpretation": "sparse/degenerate raw-candidate component; retained in the figure, strongly reduced by common shower-quality cuts, and not a selected-photon closure result",
            "source_resolution_limit": "The merged ROOT contains no TTree/RNTuple, so run/event localization requires per-run inputs or a future candidate-skim diagnostic.",
        },
        "assets": {},
        "curves": [],
        "replacement_condition": "Regenerate if the registered Au+Au triplet changes. Replace the p+p inclusive compatibility source when a newer accepted ROOT persists the required h2d shower-shape families. et2--et4 require new persisted Au+Au histogram families before their placeholder cells can be populated.",
    }
    for variable in VARIABLES:
        plot_matrix(variable, args.outdir, manifest, data_root, signal_root, inclusive_root, pp_photonjet_root, pp_inclusive_root)
    for variable in UNAVAILABLE:
        plot_unavailable(variable, args.outdir, manifest, pp_photonjet_root, pp_inclusive_root)
    (args.outdir / "provenance.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"outdir": str(args.outdir), "assets": manifest["assets"], "status": manifest["status"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
