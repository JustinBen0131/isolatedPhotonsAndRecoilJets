#!/usr/bin/env python3
"""Build the current AuAu table-QA E11/E33 data/MC overlay slide."""

from __future__ import annotations

import argparse
import importlib.util
import json
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

PP_CAMPAIGN = REPO / "dataOutput/ppg12TableQA/THE42_ppg12_tableqa_v1_basev3e_20260611"
PP_ROOT_DIR = PP_CAMPAIGN / "merged_roots"
PP_INCLUSIVE_CACHE = PP_CAMPAIGN / "inclusive_sample_hist_cache/the42_ppg12_tableqa_v1_inclusive_sample_projectx_hists.json"
PP_PLOTTER_PATH = REPO / "scripts/plotting/pp_currentian/make_the42_ppg12_tableqa_v1_tables.py"
PP_DATA_ROOT = PP_ROOT_DIR / "RecoilJets_pp_ALL_preselectionNewPPG12_tightNewPPG12_nonTightNewPPG12.root"
PP_SIGNAL_ROOT = PP_ROOT_DIR / "RecoilJets_photonjet5plus10plus20_MERGED.root"
PP_TOPDIR_DATA = "Photon_4_GeV_plus_MBD_NS_geq_1"

DEFAULT_OUTDIR = REPO / "dataOutput/ppg12TableQA/THE42_current_auau_tableqa_20260614/e11_overlay_slide"
VAR = "e11_to_e33"
PT_TOKEN = "1535"
PT_LABEL = r"$15<E_T<35$ GeV"

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


def load_tableqa_arrays(root: Path, index: dict[str, list[str]], cent: str, stage: str) -> Arrays:
    name = tableqa_name(cent, stage)
    paths = index.get(name, [])
    if not paths:
        raise KeyError(f"Missing {name} in {root}")
    f = open_root(root)
    try:
        obj = f.Get(paths[0])
        if not obj or not obj.InheritsFrom("TH1"):
            raise TypeError(f"{paths[0]} is not a TH1 in {root}")
        h = obj.Clone(f"{name}_clone")
        h.SetDirectory(0)
    finally:
        f.Close()
    return hist_to_arrays(h, source=paths[0], rebin=4)


def load_pp_plotter():
    spec = importlib.util.spec_from_file_location("pp_tableqa_plotter", PP_PLOTTER_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not import pp plotter from {PP_PLOTTER_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_pp_arrays(stage: str) -> dict[str, Arrays]:
    plotter = load_pp_plotter()
    files = {
        "data": plotter.open_root(PP_DATA_ROOT),
        "signal_mc": plotter.open_root(PP_SIGNAL_ROOT),
    }
    inclusive_cache = plotter.load_inclusive_cache(PP_INCLUSIVE_CACHE, use_stitched_inclusive=False)
    xlim, rebin = plotter.ppg12_axis_settings(VAR)

    out: dict[str, Arrays] = {}
    data = plotter.norm_arrays(plotter.get_hist(files["data"], PP_TOPDIR_DATA, VAR, PT_TOKEN, stage), rebin, xlim)
    sig = plotter.norm_arrays(plotter.get_hist(files["signal_mc"], "SIM", VAR, PT_TOKEN, stage), rebin, xlim)
    inc_payload = plotter.cache_hist(inclusive_cache, "current_ian_jet8to40", VAR, PT_TOKEN, stage)
    inc = plotter.norm_payload(inc_payload, rebin, xlim)
    out["Data"] = Arrays(*data, integral=float(np.sum(data[1])), source=f"{PP_DATA_ROOT}:{stage}")
    out["Signal MC"] = Arrays(*sig, integral=float(np.sum(sig[1])), source=f"{PP_SIGNAL_ROOT}:{stage}")
    out["Inclusive MC"] = Arrays(*inc, integral=float(np.sum(inc[1])), source=f"{PP_INCLUSIVE_CACHE}:current_ian_jet8to40:{stage}")
    if stage == "cut0":
        npb_scale = plotter.npb_tail_scale(
            files,
            {
                "pt_token": PT_TOKEN,
                "pt_label": r"$15<E_T<35$ GeV",
                "cut": stage,
                "cut_label": "all candidates",
                "include_npb_template": True,
            },
        )
        npb = plotter.norm_arrays(plotter.get_hist(files["data"], PP_TOPDIR_DATA, VAR, PT_TOKEN, "cut4"), rebin, xlim)
        npb = plotter.scale_arrays(npb, npb_scale)
        out["NPB-tagged data"] = Arrays(*npb, integral=float(np.sum(npb[1])), source=f"{PP_DATA_ROOT}:cut4 scaled")
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
    else:
        ax.step(arrays.x, arrays.y, where="mid", color=color, lw=1.85, alpha=0.97, label=label)
        if label != "NPB-tagged data":
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


def render(args: argparse.Namespace) -> dict:
    setup_style()
    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    data_root = args.data_root
    signal_root = args.signal_root
    inclusive_root = args.inclusive_root
    for path in (data_root, signal_root, inclusive_root):
        if not path.exists():
            raise FileNotFoundError(path)

    data_index = build_index(data_root)
    signal_index = build_index(signal_root)
    inclusive_index = build_index(inclusive_root)

    fig, axes = plt.subplots(4, 3, figsize=slide_figsize(), constrained_layout=False)
    fig.patch.set_facecolor("white")
    fig.subplots_adjust(left=0.070, right=0.990, top=0.755, bottom=0.070, wspace=0.120, hspace=0.235)

    manifest: dict = {
        "schema": "CURRENT_AUAU_TABLEQA_E11_OVERLAY_SLIDE_V1",
        "variable": VAR,
        "pt_token": PT_TOKEN,
        "pt_label": "15 < E_T < 35 GeV",
        "data_root": str(data_root),
        "signal_root": str(signal_root),
        "inclusive_root": str(inclusive_root),
        "pp_reference": {
            "data": str(PP_DATA_ROOT),
            "signal": str(PP_SIGNAL_ROOT),
            "inclusive_sample_cache": str(PP_INCLUSIVE_CACHE),
        },
        "curves": [],
    }

    rows = CENTRALITIES + [("pp", "pp reference")]
    for row, (cent_key, row_label) in enumerate(rows):
        for col, (stage, stage_label) in enumerate(STAGES):
            ax = axes[row, col]
            plotted: list[Arrays] = []
            if cent_key == "pp":
                curves = load_pp_arrays(stage)
            else:
                curves = {
                    "Data": load_tableqa_arrays(data_root, data_index, cent_key, stage),
                    "Signal MC": load_tableqa_arrays(signal_root, signal_index, cent_key, stage),
                    "Inclusive MC": load_tableqa_arrays(inclusive_root, inclusive_index, cent_key, stage),
                }
            for sample in ("Data", "Signal MC", "Inclusive MC", "NPB-tagged data"):
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
            ax.set_xlim(0.0, 1.0)
            ax.set_ylim(0.0, max(0.025, ymax * 1.16))
            ax.grid(True, axis="y", color="#E5E7EB", lw=0.52, alpha=0.78)
            ax.tick_params(labelsize=7.8, pad=1, direction="in", top=True, right=True)
            if col == 0:
                ax.set_ylabel(row_label.replace(" ", "\n", 1), fontsize=10.4, labelpad=8)
                add_sphenix_label(ax, system="pp" if cent_key == "pp" else "AuAu")
            if row == 0:
                ax.set_title(stage_label, fontsize=14.7, fontweight="bold", pad=6, color="#173B63")
            if row == 3:
                ax.set_xlabel(r"$E_{11}/E_{33}$", fontsize=10.7, labelpad=1)
            else:
                ax.tick_params(labelbottom=False)

    title = r"$E_{11}/E_{33}$ shower-shape overlay: current AuAu table-QA data with matched MC"
    fig.text(0.055, 0.955, title, fontsize=23.8, fontweight="bold", ha="left", va="top", color="#111827")
    draw_arrow_text(fig, 0.058, 0.895, "Rows show AuAu centrality bins and the validated pp reference; columns follow the photon-ID selection flow.")
    draw_arrow_text(fig, 0.058, 0.855, "Each panel overlays data markers with matched signal and inclusive MC shapes, normalized within the visible range.")

    handles = [
        plt.Line2D([0], [0], color=SAMPLE_COLORS["Signal MC"], lw=2.4, label="Signal MC"),
        plt.Line2D([0], [0], color=SAMPLE_COLORS["Inclusive MC"], lw=2.4, label="Inclusive MC"),
        plt.Line2D([0], [0], color=SAMPLE_COLORS["Data"], marker="o", markersize=6.5, lw=0, markerfacecolor=SAMPLE_COLORS["Data"], markeredgecolor="white", label="Data"),
        plt.Line2D([0], [0], color=SAMPLE_COLORS["NPB-tagged data"], lw=2.0, label="NPB-tagged data (pp cut0 only)"),
    ]
    fig.legend(handles=handles, loc="upper right", bbox_to_anchor=(0.985, 0.810), frameon=False, ncol=4, fontsize=9.4, handlelength=1.6, columnspacing=1.1)

    out_png = outdir / "current_auau_tableqa_e11_to_e33_data_mc_overlay_slide.png"
    out_manifest = outdir / "current_auau_tableqa_e11_to_e33_data_mc_overlay_slide_manifest.json"
    out_script = outdir / "current_auau_tableqa_e11_to_e33_data_mc_overlay_slide_speaker_script.md"
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
                "This slide is the E11 over E33 shower-shape check from the current AuAu table-QA chain.",
                "Each AuAu row is one centrality bin. In every panel, black markers are data, the red line is matched embedded signal MC, and the blue line is matched embedded inclusive MC.",
                "The columns show the selection flow: before preselection, after preselection, and after the tight BDT ID.",
                "The bottom row is the validated pp reference using the repaired table-QA plotting path, so the audience can compare the AuAu behavior against the known pp photon-ID pattern.",
                "",
                "The point to emphasize is whether the tight-ID column moves the data toward the signal-like shower-shape region while retaining a coherent inclusive-MC comparison.",
                "All curves are normalized within the plotted E11 over E33 range, so this is a shape comparison rather than a yield comparison.",
                "",
            ]
        )
    )
    return {"png": str(out_png), "manifest": str(out_manifest), "speaker_script": str(out_script)}


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-root", type=Path, default=DEFAULT_DATA_ROOT)
    ap.add_argument("--signal-root", type=Path, default=DEFAULT_SIGNAL_ROOT)
    ap.add_argument("--inclusive-root", type=Path, default=DEFAULT_INCLUSIVE_ROOT)
    ap.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    return ap.parse_args()


def main() -> int:
    outputs = render(parse_args())
    for value in outputs.values():
        print(value)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
