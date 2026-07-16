#!/usr/bin/env python3
"""Build common-contract AuAu/pp photon-ID stage-flow matrices.

The renderer is shared by slide and manuscript outputs so the pp reference is
not redrawn through a second, drifting implementation. Missing AuAu data
histograms are rendered as explicit non-data slots; available simulation is
never substituted for data.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import re
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


THIS_FILE = Path(__file__).resolve()
REPO = next((p for p in THIS_FILE.parents if (p / "AGENTS.md").exists()), THIS_FILE.parents[4])
SCRIPTS_DIR = REPO / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.append(str(SCRIPTS_DIR))

from slides.common.slide_defaults import SLIDE_DPI, slide_figsize  # noqa: E402


CFG = "preselectionNewPPG12_tightAuAuCentInputBase3x3BDT_nonTightAuAuBDTComplement_baseVariant"
LOCAL_INPUT = REPO / "InputFiles/the42_current_auau_tableqa"
DEFAULT_DATA_ROOT = LOCAL_INPUT / f"RecoilJets_auau_ALL_{CFG}.root"
DEFAULT_SIGNAL_ROOT = LOCAL_INPUT / "RecoilJets_embeddedPhoton12plus20_MERGED.root"
DEFAULT_INCLUSIVE_ROOT = LOCAL_INPUT / "RecoilJets_embeddedJet12plus20plus30plus40_MERGED.root"
DEFAULT_INTERIM_DATA_CACHE = (
    REPO
    / "dataOutput/ppg12TableQA/THE42_ppg12_tableqa_v1_auauDataNewV008_b002_20260612"
    / "interim_complete_runs/the42_b002_interim_complete_run_shower_shape_hists.json"
)

PP_PLOTTER_PATH = REPO / "scripts/plotting/pp_currentian/make_the42_ppg12_tableqa_v1_tables.py"


def registered_root(sample_key: str) -> tuple[Path, Path]:
    pointer = REPO / "dataOutput/current_recoiljets_artifacts/current" / sample_key / "current.json"
    payload = json.loads(pointer.read_text())
    roots = payload.get("root_paths", [])
    if len(roots) != 1:
        raise RuntimeError(f"Expected one registered ROOT in {pointer}, found {len(roots)}")
    return Path(roots[0]), pointer


PP_DATA_ROOT, PP_DATA_POINTER = registered_root("pp_data_merged")
PP_SIGNAL_ROOT, PP_SIGNAL_POINTER = registered_root("pp_sim_photonjet_merged")
PP_INCLUSIVE_ROOT, PP_INCLUSIVE_POINTER = registered_root("pp_sim_inclusivejet_merged")
PP_TOPDIR_DATA = "PPG12_scaledtrigger30"

DEFAULT_OUTDIR = REPO / "dataOutput/ppg12TableQA/THE42_current_auau_tableqa_20260614/e11_overlay_slide"
VAR = "e11_to_e33"
PT_TOKEN = "1535"
PT_LABEL = r"$15<E_T<35$ GeV"
VARIABLE_CONFIG = {
    "weta_cogx": {
        "axis": r"$w_{\eta}^{\mathrm{COGX}}$",
        "slug": "weta_cogx",
        "title": r"$w_{\eta}^{\mathrm{COGX}}$ shower-shape overlay: AuAu centrality and pp reference",
        "lead": "Rows show AuAu centrality bins and the pp reference; columns follow the photon-ID selection flow.",
        "takeaway": "The width evolution tests whether the selection retains the same compact-shower ordering across occupancy and collision system.",
    },
    "wphi_cogx": {
        "axis": r"$w_{\phi}^{\mathrm{COGX}}$",
        "slug": "wphi_cogx",
        "title": r"$w_{\phi}^{\mathrm{COGX}}$ shower-shape overlay: AuAu centrality and pp reference",
        "lead": "Rows show AuAu centrality bins and the pp reference; columns follow the photon-ID selection flow.",
        "takeaway": "The azimuthal width provides the complementary lateral-shape check to the pseudorapidity width.",
    },
    "e11_to_e33": {
        "axis": r"$E_{11}/E_{33}$",
        "slug": "e11_to_e33",
        "title": r"$E_{11}/E_{33}$ shower-shape overlay: AuAu centrality and pp reference",
        "lead": "Rows show AuAu centrality bins and the validated pp reference; columns follow the photon-ID selection flow.",
        "takeaway": "The tight-ID column should move the data toward the signal-like shower-shape region while retaining a coherent inclusive-MC comparison.",
    },
    "et1": {
        "axis": r"$\mathrm{et1}$",
        "slug": "et1",
        "title": r"$\mathrm{et1}$ energy-sharing overlay: AuAu centrality and pp reference",
        "lead": "Rows show AuAu centrality bins and the pp reference; columns follow the photon-ID selection flow.",
        "takeaway": "The first ordered energy-sharing input checks whether the same shower-core partition is populated in data and simulation.",
    },
    "et2": {
        "axis": r"$\mathrm{et2}$",
        "slug": "et2",
        "title": r"$\mathrm{et2}$ energy-sharing overlay: AuAu centrality and pp reference",
        "lead": "Rows show AuAu centrality bins and the pp reference; columns follow the photon-ID selection flow.",
        "takeaway": "The second ordered energy-sharing input tests a complementary component of the cluster topology.",
    },
    "et3": {
        "axis": r"$\mathrm{et3}$",
        "slug": "et3",
        "title": r"$\mathrm{et3}$ energy-sharing overlay: AuAu centrality and pp reference",
        "lead": "Rows show AuAu centrality bins and the pp reference; columns follow the photon-ID selection flow.",
        "takeaway": "The third ordered energy-sharing input tests whether broader shower energy is modeled consistently.",
    },
    "et4": {
        "axis": r"$\mathrm{et4}$",
        "slug": "et4",
        "title": r"$\mathrm{et4}$ energy-sharing overlay: AuAu centrality and pp reference",
        "lead": "Rows show AuAu centrality bins and the pp reference; columns follow the photon-ID selection flow.",
        "takeaway": "The fourth ordered energy-sharing input is especially sensitive to small outlying components of the cluster energy pattern.",
    },
    "e32_to_e35": {
        "axis": r"$E_{3\times2}/E_{3\times5}$",
        "slug": "e32_to_e35",
        "title": r"$E_{3\times2}/E_{3\times5}$ containment overlay: AuAu centrality and pp reference",
        "lead": "Rows show AuAu centrality bins and the pp reference; columns follow the photon-ID selection flow.",
        "takeaway": "The elongated-core containment ratio tests compactness in a geometry complementary to $E_{11}/E_{33}$.",
    },
    "bdt": {
        "axis": "BDT score",
        "slug": "bdt_score",
        "title": "BDT-score separation overlay: current AuAu default BDT data with matched MC",
        "lead": "Rows show AuAu centrality bins and the pp reference; columns follow the photon-ID selection flow.",
        "takeaway": "Signal MC should concentrate at higher score than inclusive MC; data should sit between the two without pathological pileups.",
    },
}

STAGES = [
    ("cut0", "Before preselection"),
    ("cut1", "After preselection"),
    ("cut2", "After tight ID"),
]
CENTRALITIES = [
    ("cent0_20", "AuAu 0-20%"),
    ("cent20_50", "AuAu 20-50%"),
    ("cent50_80", "AuAu 50-80%"),
]
CENTRALITY_MIDPOINTS = {
    "cent0_20": 10.0,
    "cent20_50": 35.0,
    "cent50_80": 65.0,
}
WP80_INTERCEPT = 0.53471108
WP80_SLOPE = 0.0012284143
SAMPLE_COLORS = {
    "Data": "#111827",
    "Signal MC": "#C22F2F",
    "Inclusive MC": "#2E63D4",
    "NPB-tagged data": "#238B1E",
}


@dataclass
class Arrays:
    x: np.ndarray
    y: np.ndarray
    e: np.ndarray
    integral: float
    source: str


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
            "axes.facecolor": "white",
            "axes.edgecolor": "#111827",
            "axes.linewidth": 0.85,
            "xtick.color": "#111827",
            "ytick.color": "#111827",
        }
    )


def open_root(path: Path):
    import ROOT

    ROOT.gROOT.SetBatch(True)
    ROOT.TH1.AddDirectory(False)
    f = ROOT.TFile.Open(str(path), "READ")
    if not f or f.IsZombie():
        raise OSError(f"Could not open ROOT file: {path}")
    return f


def walk(directory, prefix: str = ""):
    for key in directory.GetListOfKeys():
        obj = key.ReadObj()
        path = f"{prefix}/{key.GetName()}" if prefix else key.GetName()
        yield path, obj
        if obj.InheritsFrom("TDirectory"):
            yield from walk(obj, path)


def build_index(path: Path) -> dict[str, list[str]]:
    f = open_root(path)
    index: dict[str, list[str]] = {}
    try:
        for obj_path, obj in walk(f):
            if obj.InheritsFrom("TH1"):
                index.setdefault(Path(obj_path).name, []).append(obj_path)
    finally:
        f.Close()
    return index


def tableqa_name(cent: str, stage: str) -> str:
    return f"h1d_{VAR}_eta0_pt{PT_TOKEN}_{cent}_{stage}"


def rebin_arrays(x: np.ndarray, y: np.ndarray, e: np.ndarray, factor: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if factor <= 1 or x.size < factor:
        return x, y, e
    n = x.size // factor
    trim = n * factor
    xr = x[:trim].reshape(n, factor)
    yr = y[:trim].reshape(n, factor)
    er = e[:trim].reshape(n, factor)
    return xr.mean(axis=1), yr.sum(axis=1), np.sqrt(np.sum(er * er, axis=1))


def coarsen_payload(payload: tuple[np.ndarray, np.ndarray, np.ndarray], factor: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x, y, e = payload
    return rebin_arrays(np.asarray(x, dtype=float), np.asarray(y, dtype=float), np.asarray(e, dtype=float), factor)


def coarsen_arrays(arr: Arrays, factor: int) -> Arrays:
    x, y, e = rebin_arrays(arr.x, arr.y, arr.e, factor)
    suffix = f"; display_coarsen={factor}" if factor > 1 else ""
    return Arrays(x=x, y=y, e=e, integral=arr.integral, source=f"{arr.source}{suffix}")


def compress_sideband_for_display(arr: Arrays, target_peak: float) -> Arrays:
    peak = float(np.nanmax(arr.y)) if arr.y.size else 0.0
    if peak <= 0.0 or target_peak <= 0.0:
        return Arrays(x=arr.x, y=np.zeros_like(arr.y), e=np.zeros_like(arr.e), integral=arr.integral, source=f"{arr.source}; compressed_sideband_display")
    scale = target_peak / peak
    return Arrays(
        x=arr.x,
        y=arr.y * scale,
        e=arr.e * scale,
        integral=arr.integral,
        source=f"{arr.source}; compressed_sideband_display_scale={scale:.6g}",
    )


def hist_to_arrays(hist, *, source: str, xlim: tuple[float, float] = (0.0, 1.0), rebin: int = 4) -> Arrays:
    nb = hist.GetNbinsX()
    x = np.array([hist.GetXaxis().GetBinCenter(i) for i in range(1, nb + 1)], dtype=float)
    y = np.array([hist.GetBinContent(i) for i in range(1, nb + 1)], dtype=float)
    e = np.array([hist.GetBinError(i) for i in range(1, nb + 1)], dtype=float)
    x, y, e = rebin_arrays(x, y, e, rebin)
    mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(e) & (x >= xlim[0]) & (x <= xlim[1])
    x, y, e = x[mask], y[mask], e[mask]
    integral = float(np.sum(y))
    if integral > 0:
        y = y / integral
        e = e / integral
    return Arrays(x=x, y=y, e=e, integral=integral, source=source)


def load_tableqa_arrays(
    root: Path,
    index: dict[str, list[str]],
    cent: str,
    stage: str,
    *,
    path_regex: re.Pattern[str] | None = None,
) -> Arrays:
    name = tableqa_name(cent, stage)
    paths = index.get(name, [])
    if not paths:
        raise KeyError(f"Missing {name} in {root}")
    if path_regex is not None:
        matched = [p for p in paths if path_regex.search(p)]
        if matched:
            paths = matched
    f = open_root(root)
    try:
        obj = f.Get(paths[0])
        if not obj or not obj.InheritsFrom("TH1"):
            raise TypeError(f"{paths[0]} is not a TH1 in {root}")
        h = obj.Clone(f"{name}_clone")
        h.SetDirectory(0)
    finally:
        f.Close()
    xlim, rebin = load_pp_plotter().ppg12_axis_settings(VAR)
    return hist_to_arrays(h, source=paths[0], xlim=xlim, rebin=rebin)


def load_interim_cache(path: Path) -> dict:
    payload = json.loads(path.read_text())
    if payload.get("schema") != "THE42_B002_INTERIM_COMPLETE_RUN_SHOWER_SHAPES_V1":
        raise RuntimeError(f"Unexpected interim data cache schema in {path}")
    return payload


def load_cache_arrays(cache: dict, cache_path: Path, cent: str, stage: str) -> Arrays | None:
    variable_payload = cache.get("hists", {}).get(VAR)
    if variable_payload is None:
        return None
    item = variable_payload.get(cent, {}).get(stage)
    if item is None:
        return None
    x = np.asarray(item["x"], dtype=float)
    y = np.asarray(item["y"], dtype=float)
    e = np.asarray(item["e"], dtype=float)
    xlim, rebin = load_pp_plotter().ppg12_axis_settings(VAR)
    x, y, e = rebin_arrays(x, y, e, rebin)
    mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(e) & (x >= xlim[0]) & (x <= xlim[1])
    x, y, e = x[mask], y[mask], e[mask]
    integral = float(np.sum(y))
    if integral > 0:
        y = y / integral
        e = e / integral
    return Arrays(x=x, y=y, e=e, integral=integral, source=f"{cache_path}:{VAR}:{cent}:{stage}")


def load_pp_plotter():
    spec = importlib.util.spec_from_file_location("pp_tableqa_plotter", PP_PLOTTER_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not import pp plotter from {PP_PLOTTER_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_pp_arrays(stage: str, *, npb_display_coarsen: int = 1, include_npb: bool = True) -> dict[str, Arrays]:
    plotter = load_pp_plotter()
    files = {
        "data": plotter.open_root(PP_DATA_ROOT),
        "signal_mc": plotter.open_root(PP_SIGNAL_ROOT),
        "inclusive_mc": plotter.open_root(PP_INCLUSIVE_ROOT),
    }
    xlim, rebin = plotter.ppg12_axis_settings(VAR)

    out: dict[str, Arrays] = {}
    data = plotter.norm_arrays(plotter.get_hist(files["data"], PP_TOPDIR_DATA, VAR, PT_TOKEN, stage), rebin, xlim)
    sig = plotter.norm_arrays(plotter.get_hist(files["signal_mc"], "SIM", VAR, PT_TOKEN, stage), rebin, xlim)
    inc = plotter.norm_arrays(plotter.get_hist(files["inclusive_mc"], "SIM", VAR, PT_TOKEN, stage), rebin, xlim)
    if data is None or sig is None or inc is None:
        raise RuntimeError(f"Missing registered pp histogram for {VAR}, {stage}")
    out["Data"] = Arrays(*data, integral=float(np.sum(data[1])), source=f"{PP_DATA_ROOT}:{stage}")
    out["Signal MC"] = Arrays(*sig, integral=float(np.sum(sig[1])), source=f"{PP_SIGNAL_ROOT}:{stage}")
    out["Inclusive MC"] = Arrays(*inc, integral=float(np.sum(inc[1])), source=f"{PP_INCLUSIVE_ROOT}:{stage}")
    if include_npb and stage == "cut0":
        npb_hist = plotter.get_hist(files["data"], PP_TOPDIR_DATA, VAR, PT_TOKEN, "cut4")
        npb_raw_entries = float(npb_hist.GetEntries())
        npb = plotter.norm_arrays(npb_hist, rebin, xlim)
        npb = coarsen_payload(npb, npb_display_coarsen)
        out["NPB-tagged data"] = Arrays(
            *npb,
            integral=npb_raw_entries,
            source=f"{PP_DATA_ROOT}:cut4 raw-sideband-shape; display_coarsen={npb_display_coarsen}",
        )
    for root_file in files.values():
        root_file.Close()
    return out


def draw_curve(ax, arrays: Arrays, label: str) -> None:
    color = SAMPLE_COLORS[label]
    if label == "Data":
        ax.errorbar(
            arrays.x,
            arrays.y,
            yerr=arrays.e,
            fmt="o",
            ms=3.5,
            color=color,
            mfc=color,
            mec="white",
            mew=0.55,
            elinewidth=0.75,
            capsize=1.3,
            capthick=0.7,
            alpha=0.96,
            label=label,
        )
    elif label == "NPB-tagged data":
        ax.fill_between(arrays.x, 0.0, arrays.y, step="mid", color=color, alpha=0.115, linewidth=0)
        ax.step(arrays.x, arrays.y, where="mid", color=color, lw=3.0, alpha=0.98, label=label)
        ax.plot(
            arrays.x,
            arrays.y,
            "s",
            ms=2.8,
            color=color,
            mfc="white",
            mec=color,
            mew=0.8,
            alpha=0.98,
        )
    else:
        ax.step(arrays.x, arrays.y, where="mid", color=color, lw=1.85, alpha=0.97, label=label)
        ax.errorbar(arrays.x, arrays.y, yerr=arrays.e, fmt="none", ecolor=color, elinewidth=0.38, alpha=0.42)


def add_sphenix_label(ax, *, system: str) -> None:
    if system == "pp":
        collision = r"$p$+$p$ $\sqrt{s}=200$ GeV"
    else:
        collision = r"Au+Au $\sqrt{s_{NN}}=200$ GeV"
    ax.text(0.035, 0.945, r"$\it{\bf{sPHENIX}}$ Internal", transform=ax.transAxes, ha="left", va="top", fontsize=7.7)
    ax.text(0.035, 0.790, f"{collision}\n{PT_LABEL}, $|\\eta|<0.7$", transform=ax.transAxes, ha="left", va="top", fontsize=6.7, linespacing=1.03)


def draw_arrow_text(fig, x: float, y: float, text: str) -> None:
    fig.text(x, y, "▶", fontsize=14.5, color="#2468A8", ha="left", va="top", fontfamily="DejaVu Sans")
    fig.text(x + 0.021, y, text, fontsize=14.2, color="#172033", ha="left", va="top")


def draw_missing_data_slot(ax, *, row_label: str, stage_label: str) -> None:
    """Render an unmistakable non-data slot when the bounded cache lacks a variable."""
    ax.set_facecolor("#FFF9ED")
    for spine in ax.spines.values():
        spine.set_color("#C47A14")
        spine.set_linewidth(1.15)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.text(
        0.5,
        0.59,
        "AUAU DATA SLOT",
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=8.4,
        fontweight="bold",
        color="#A65C00",
    )
    ax.text(
        0.5,
        0.38,
        "No common-contract data histogram\nin the bounded cache",
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=6.7,
        color="#6B4A1D",
        linespacing=1.15,
    )
    ax.text(
        0.5,
        0.14,
        f"{row_label}; {stage_label}",
        transform=ax.transAxes,
        ha="center",
        va="center",
        fontsize=6.2,
        color="#6B4A1D",
    )


def render(args: argparse.Namespace) -> dict:
    global VAR
    VAR = args.var
    var_cfg = VARIABLE_CONFIG[VAR]
    include_npb = bool(args.include_npb_sideband and VAR != "bdt")
    setup_style()
    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    data_root = args.data_root
    signal_root = args.signal_root
    inclusive_root = args.inclusive_root
    data_cache_path = args.data_cache
    for path in (signal_root, inclusive_root):
        if not path.exists():
            raise FileNotFoundError(path)
    if data_cache_path is None and not data_root.exists():
        raise FileNotFoundError(data_root)
    if data_cache_path is not None and not data_cache_path.exists():
        raise FileNotFoundError(data_cache_path)

    data_cache = load_interim_cache(data_cache_path) if data_cache_path is not None else None
    data_index = build_index(data_root) if data_cache is None else {}
    signal_index = build_index(signal_root)
    inclusive_index = build_index(inclusive_root)
    data_trigger_regex = re.compile(args.data_trigger_regex) if args.data_trigger_regex else None

    manuscript_layout = args.layout == "manuscript"
    fig_size = (7.35, 8.7) if manuscript_layout else slide_figsize()
    fig, axes = plt.subplots(4, 3, figsize=fig_size, constrained_layout=False)
    fig.patch.set_facecolor("white")
    if manuscript_layout:
        fig.subplots_adjust(left=0.105, right=0.985, top=0.875, bottom=0.065, wspace=0.135, hspace=0.215)
    else:
        fig.subplots_adjust(left=0.070, right=0.990, top=0.720, bottom=0.070, wspace=0.120, hspace=0.220)

    manifest: dict = {
        "schema": "CURRENT_AUAU_TABLEQA_OVERLAY_SLIDE_V2",
        "variable": VAR,
        "pt_token": PT_TOKEN,
        "pt_label": "15 < E_T < 35 GeV",
        "data_root": str(data_root) if data_cache is None else None,
        "data_interim_cache": str(data_cache_path) if data_cache is not None else None,
        "data_interim_cache_metadata": data_cache.get("metadata", {}) if data_cache is not None else None,
        "signal_root": str(signal_root),
        "inclusive_root": str(inclusive_root),
        "pp_reference": {
            "data": str(PP_DATA_ROOT),
            "data_pointer": str(PP_DATA_POINTER),
            "signal": str(PP_SIGNAL_ROOT),
            "signal_pointer": str(PP_SIGNAL_POINTER),
            "inclusive": str(PP_INCLUSIVE_ROOT),
            "inclusive_pointer": str(PP_INCLUSIVE_POINTER),
        },
        "npb_sideband_note": (
            "Green is raw cut4 NPB-tagged data in both pp and AuAu. "
            "It is drawn as a compressed diagnostic strip and does not set the unit-area data/MC y-axis."
            if include_npb
            else None
        ),
        "npb_display_coarsen": args.npb_display_coarsen,
        "include_npb_sideband": include_npb,
        "layout": args.layout,
        "curves": [],
        "missing_data_slots": [],
    }

    rows = CENTRALITIES + [("pp", "pp reference")]
    for row, (cent_key, row_label) in enumerate(rows):
        for col, (stage, stage_label) in enumerate(STAGES):
            ax = axes[row, col]
            plotted: list[Arrays] = []
            if cent_key == "pp":
                curves = load_pp_arrays(stage, npb_display_coarsen=args.npb_display_coarsen, include_npb=include_npb)
            else:
                data = (
                    load_cache_arrays(data_cache, data_cache_path, cent_key, stage)
                    if data_cache is not None
                    else load_tableqa_arrays(data_root, data_index, cent_key, stage, path_regex=data_trigger_regex)
                )
                if data is None:
                    draw_missing_data_slot(ax, row_label=row_label, stage_label=stage_label)
                    if row == 0:
                        ax.set_title(
                            stage_label,
                            fontsize=10.2 if manuscript_layout else 14.7,
                            fontweight="bold",
                            pad=5,
                            color="#173B63",
                        )
                    if col == 0:
                        ax.set_ylabel(
                            row_label.replace(" ", "\n", 1),
                            fontsize=8.4 if manuscript_layout else 10.4,
                            labelpad=7,
                        )
                    manifest["missing_data_slots"].append(
                        {"row": row_label, "stage": stage, "reason": f"{VAR} absent from bounded AuAu data cache"}
                    )
                    continue
                curves = {
                    "Data": data,
                    "Signal MC": load_tableqa_arrays(signal_root, signal_index, cent_key, stage),
                    "Inclusive MC": load_tableqa_arrays(inclusive_root, inclusive_index, cent_key, stage),
                }
                if include_npb and stage == "cut0":
                    npb = (
                        load_cache_arrays(data_cache, data_cache_path, cent_key, "cut4")
                        if data_cache is not None
                        else load_tableqa_arrays(data_root, data_index, cent_key, "cut4", path_regex=data_trigger_regex)
                    )
                    if npb is not None:
                        npb = coarsen_arrays(npb, args.npb_display_coarsen)
                        curves["NPB-tagged data"] = npb
            for sample in ("Data", "Signal MC", "Inclusive MC"):
                arr = curves.get(sample)
                if arr is None:
                    continue
                draw_curve(ax, arr, sample)
                plotted.append(arr)
                manifest["curves"].append(
                    {
                        "row": row_label,
                        "stage": stage,
                        "sample": sample,
                        "integral_before_normalization": arr.integral,
                        "source": arr.source,
                    }
                )
            ymax = max((float(np.nanmax(a.y + a.e)) for a in plotted if a.y.size), default=0.02)
            npb = curves.get("NPB-tagged data") if include_npb else None
            if npb is not None:
                npb_display = compress_sideband_for_display(npb, target_peak=max(0.012, ymax * 0.185))
                draw_curve(ax, npb_display, "NPB-tagged data")
                manifest["curves"].append(
                    {
                        "row": row_label,
                        "stage": stage,
                        "sample": "NPB-tagged data",
                        "integral_before_normalization": npb.integral,
                        "source": npb.source,
                        "display_source": npb_display.source,
                        "display_mode": "compressed_sideband_strip",
                    }
                )
            if include_npb and stage == "cut0":
                npb = curves.get("NPB-tagged data")
                if npb is not None:
                    note = "cut4 NPB: N=0" if npb.integral <= 0 else f"cut4 NPB: N={npb.integral:.0f}"
                    ax.text(
                        0.975,
                        0.905,
                        note,
                        transform=ax.transAxes,
                        ha="right",
                        va="top",
                        fontsize=10.2,
                        color=SAMPLE_COLORS["NPB-tagged data"],
                        fontweight="bold",
                        bbox={"boxstyle": "round,pad=0.16", "facecolor": "white", "edgecolor": "none", "alpha": 0.78},
                    )
            xlim, _ = load_pp_plotter().ppg12_axis_settings(VAR)
            ax.set_xlim(*xlim)
            ax.set_ylim(0.0, max(0.025, ymax * 1.16))
            if VAR == "bdt" and cent_key in CENTRALITY_MIDPOINTS:
                wp80 = WP80_INTERCEPT + WP80_SLOPE * CENTRALITY_MIDPOINTS[cent_key]
                ax.axvline(wp80, color="#4B5563", lw=1.05, ls=(0, (3.2, 2.6)), alpha=0.72, zorder=0)
            ax.grid(True, axis="y", color="#E5E7EB", lw=0.52, alpha=0.78)
            ax.tick_params(labelsize=6.4 if manuscript_layout else 7.8, pad=1, direction="in", top=True, right=True)
            if col == 0:
                ax.set_ylabel(row_label.replace(" ", "\n", 1), fontsize=8.4 if manuscript_layout else 10.4, labelpad=7)
                add_sphenix_label(ax, system="pp" if cent_key == "pp" else "AuAu")
            if row == 0:
                ax.set_title(stage_label, fontsize=10.2 if manuscript_layout else 14.7, fontweight="bold", pad=5, color="#173B63")
            if row == 3:
                ax.set_xlabel(var_cfg["axis"], fontsize=8.3 if manuscript_layout else 10.7, labelpad=1)
            else:
                ax.tick_params(labelbottom=False)

    title = var_cfg["title"]
    fig.text(
        0.105 if manuscript_layout else 0.055,
        0.975 if manuscript_layout else 0.955,
        title,
        fontsize=13.2 if manuscript_layout else 23.8,
        fontweight="bold",
        ha="left",
        va="top",
        color="#111827",
    )
    if data_cache is not None:
        meta = data_cache.get("metadata", {})
        data_line = (
            f"Interim data subset: {meta.get('complete_runs', '?')} complete runs, "
            f"{meta.get('complete_root_files_used', '?')} completed ROOT chunks; final merge will replace this when the tail drains."
        )
    else:
        data_line = var_cfg["lead"]
    if manuscript_layout:
        fig.text(
            0.105,
            0.952,
            "Historical diagnostic: common 15--35 GeV table-QA contract; unit-area shapes",
            fontsize=8.0,
            ha="left",
            va="top",
            color="#4B5563",
        )
    else:
        draw_arrow_text(fig, 0.058, 0.895, data_line)
    if include_npb and not manuscript_layout:
        draw_arrow_text(
            fig,
            0.058,
            0.855,
            "Green is the same raw cut4 NPB-tagged data sideband in pp and AuAu; entries are printed in every before-preselection panel.",
        )
        draw_arrow_text(
            fig,
            0.058,
            0.818,
            "The green height is compressed into a diagnostic strip, so sparse AuAu sidebands do not look like high-stat unit-area shapes.",
        )
    elif not manuscript_layout:
        draw_arrow_text(fig, 0.058, 0.855, var_cfg["takeaway"])
        if VAR == "bdt":
            draw_arrow_text(
                fig,
                0.058,
                0.818,
                "AuAu uses the new 14-feature default BDT/WP80 production; the pp row is a PPG12/baseV3E reference-score shape, not a shared score calibration.",
            )

    handles = [
        plt.Line2D([0], [0], color=SAMPLE_COLORS["Signal MC"], lw=2.4, label="Signal MC"),
        plt.Line2D([0], [0], color=SAMPLE_COLORS["Inclusive MC"], lw=2.4, label="Inclusive MC"),
        plt.Line2D([0], [0], color=SAMPLE_COLORS["Data"], marker="o", markersize=6.5, lw=0, markerfacecolor=SAMPLE_COLORS["Data"], markeredgecolor="white", label="Data"),
    ]
    if VAR == "bdt":
        handles.append(plt.Line2D([0], [0], color="#4B5563", lw=1.2, ls=(0, (3.2, 2.6)), label="AuAu WP80"))
    if include_npb:
        handles.append(
            plt.Line2D(
            [0],
            [0],
            color=SAMPLE_COLORS["NPB-tagged data"],
            lw=2.9,
            marker="s",
            markersize=4.4,
            markerfacecolor="white",
            markeredgecolor=SAMPLE_COLORS["NPB-tagged data"],
            label="NPB-tagged sideband strip",
            )
        )
    if manuscript_layout:
        fig.legend(
            handles=handles,
            loc="upper right",
            bbox_to_anchor=(0.985, 0.925),
            frameon=False,
            ncol=4,
            fontsize=6.8,
            handlelength=1.4,
            columnspacing=0.8,
        )
    else:
        fig.legend(handles=handles, loc="upper right", bbox_to_anchor=(0.985, 0.800), frameon=False, ncol=4, fontsize=9.2, handlelength=1.6, columnspacing=1.0)

    if manuscript_layout:
        out_stem = f"auau_pp_{var_cfg['slug']}_stage_matrix_historical"
    else:
        out_stem = "current_auau_tableqa_e11_to_e33_data_mc_overlay_slide" if VAR == "e11_to_e33" else f"current_auau_tableqa_{var_cfg['slug']}_data_mc_overlay_slide"
    out_png = outdir / f"{out_stem}.png"
    out_manifest = outdir / f"{out_stem}_manifest.json"
    out_script = outdir / f"{out_stem}_speaker_script.md"
    fig.savefig(out_png, dpi=SLIDE_DPI)
    plt.close(fig)

    manifest["png"] = str(out_png)
    manifest["speaker_script"] = str(out_script)
    out_manifest.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    out_script.write_text(
        "\n".join(
            [
                "# Speaker Script",
                "",
                f"This slide is the {var_cfg['axis']} check from the current AuAu table-QA chain.",
                "Each AuAu row is one centrality bin. In every panel, black markers are data, the red line is matched embedded signal MC, and the blue line is matched embedded inclusive MC.",
                "The columns show the selection flow: before preselection, after preselection, and after the tight BDT ID.",
                "The bottom row is the validated pp reference using the repaired table-QA plotting path, so the audience can compare the AuAu behavior against the known pp photon-ID pattern.",
                "The green curve is omitted for the BDT-score slide." if not include_npb else "The green curve is the raw cut4 NPB-tagged data sideband in both systems. It is compressed into a diagnostic strip and does not share the unit-area y-scale of the data/MC shape overlays.",
                "This rendered version uses the interim completed-run data subset if a data cache is recorded in the manifest. The final merged data ROOT can be substituted without changing the slide layout.",
                "",
                var_cfg["takeaway"],
                f"All curves are normalized within the plotted {var_cfg['axis']} range, so this is a shape comparison rather than a yield comparison.",
                "",
            ]
        )
    )
    return {"png": str(out_png), "manifest": str(out_manifest), "speaker_script": str(out_script)}


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", type=Path, default=DEFAULT_DATA_ROOT)
    ap.add_argument("--var", choices=sorted(VARIABLE_CONFIG), default="e11_to_e33")
    ap.add_argument("--all-basev3e-shapes", action="store_true")
    ap.add_argument("--layout", choices=("slide", "manuscript"), default="slide")
    ap.add_argument("--data-cache", type=Path, default=None)
    ap.add_argument("--use-default-interim-cache", action="store_true")
    ap.add_argument("--signal-root", type=Path, default=DEFAULT_SIGNAL_ROOT)
    ap.add_argument("--inclusive-root", type=Path, default=DEFAULT_INCLUSIVE_ROOT)
    ap.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    ap.add_argument(
        "--data-trigger-regex",
        default=r"photon_12_plus_MBD_NS_geq_2_vtx_lt_150/",
        help="Regex used to choose the AuAu data trigger directory when duplicate table-QA histograms exist.",
    )
    ap.add_argument(
        "--npb-display-coarsen",
        type=int,
        default=2,
        help="Extra display-only bin grouping for the sparse NPB-tagged data sideband.",
    )
    ap.add_argument(
        "--include-npb-sideband",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Draw the cut4 NPB sideband strip when the selected variable supports it.",
    )
    args = ap.parse_args()
    if args.use_default_interim_cache:
        args.data_cache = DEFAULT_INTERIM_DATA_CACHE
    return args


def main() -> int:
    args = parse_args()
    variables = (
        ["weta_cogx", "wphi_cogx", "e11_to_e33", "et1", "et2", "et3", "et4", "e32_to_e35"]
        if args.all_basev3e_shapes
        else [args.var]
    )
    for variable in variables:
        args.var = variable
        outputs = render(args)
        for value in outputs.values():
            print(value)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
