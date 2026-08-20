#!/usr/bin/env python3
"""Render the THE-100 complement isolation-template JSTG candidate.

The Au+Au panels use only the stored sliding-R=0.4 ABCD isolation views:
``A+B`` is the tight-BDT population and ``C+D`` is the broad non-tight BDT
complement.  This avoids treating THE-100's unsuffixed fixed-isolation
histograms as the nominal Au+Au selection.

The Au+Au horizontal range is symmetric about zero (``--auau-xmin`` /
``--auau-xmax``) so the full UE-subtraction tail is shown on both sides and no
isolation-energy population is hidden at the panel edges.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import uproot
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Rectangle


REPO = Path(__file__).resolve().parents[3]
OUTDIR = REPO / "dataOutput/slides/the45_jstg_20260720/the100_complement_isolation_templates"
PNG = OUTDIR / "the100_complement_isolation_template_comparison_pp_current_auau_sliding.png"
MANIFEST = OUTDIR / "the100_complement_isolation_template_comparison_manifest.json"
SPEAKER = OUTDIR / "the100_complement_isolation_template_comparison_speaker_script.md"
LAYOUT = OUTDIR / "the100_complement_isolation_template_comparison_layout_nodes.json"

THE100_DATA = REPO / (
    "InputFiles/the100_auau_dualview_20260714/complement/data/"
    "RecoilJets_auau_ALL_preselectionNewPPG12_tightAuAuCentInputBase3x3BDT_"
    "nonTightAuAuBDTComplement_baseVariant.root"
)
THE100_SIGNAL = REPO / (
    "InputFiles/the100_auau_dualview_20260714/complement/signal/"
    "RecoilJets_embeddedPhoton12plus20_MERGED.root"
)
THE100_INCLUSIVE_REMOTE = (
    "/sphenix/u/patsfan753/scratch/thesisAnalysis/runs/recoiljets/current/"
    "the100_auau_dualview_20260714/inclusive/simembeddedinclusive/"
    "preselectionNewPPG12_tightAuAuCentInputBase3x3BDT_nonTightAuAuBDTComplement_baseVariant/"
    "embeddedJet12and20and30and40merged_SIM/"
    "RecoilJets_embeddedJet12plus20plus30plus40_MERGED.root"
)

DATA_TOP = "MBD_NS_geq_2_vtx_lt_150"
SIM_TOP = "SIM"
PP_BINS = ("16_18", "18_20", "20_22", "22_24", "24_26", "26_35")
AUAU_BINS = ("15_17", "17_19", "19_21", "21_23", "23_26", "26_35")
CENTRALITIES = (("0_20", "0–20%"), ("20_50", "20–50%"), ("50_80", "50–80%"))

INK = "#142235"
MUTED = "#52657A"
GRID = "#D7E0EA"
RED_FILL = "#F3B0AD"
RED_EDGE = "#B33D38"
BLUE_FILL = "#AAADEE"
BLUE_EDGE = "#2F51A7"
BLACK = "#111827"
BLUE_BULLET = "#2468A8"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def current_root(name: str) -> Path:
    pointer = REPO / "dataOutput/current_recoiljets_artifacts/current" / name / "current.json"
    payload = json.loads(pointer.read_text())
    for field in ("root", "root_path", "path"):
        if field in payload:
            return Path(payload[field])
    if payload.get("root_paths"):
        return Path(payload["root_paths"][0])
    raise KeyError(f"{pointer} has no ROOT path field")


def read_hist(handle: uproot.ReadOnlyDirectory, key: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    hist = handle[key]
    values, edges = hist.to_numpy(flow=False)
    variances = hist.variances(flow=False)
    if variances is None:
        variances = np.clip(values, 0.0, None)
    return np.asarray(values, dtype=float), np.asarray(variances, dtype=float), np.asarray(edges, dtype=float)


def sum_histograms(handle: uproot.ReadOnlyDirectory, keys: list[str]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    values_sum = variances_sum = edges_ref = None
    for key in keys:
        values, variances, edges = read_hist(handle, key)
        if values_sum is None:
            values_sum, variances_sum, edges_ref = values.copy(), variances.copy(), edges
        else:
            if not np.array_equal(edges, edges_ref):
                raise ValueError(f"incompatible isolation binning: {key}")
            values_sum += values
            variances_sum += variances
    assert values_sum is not None and variances_sum is not None and edges_ref is not None
    return values_sum, variances_sum, edges_ref


def rebin_variable(values: np.ndarray, variances: np.ndarray, edges: np.ndarray, low_group: int, high_group: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    out_values, out_variances, out_edges = [], [], [float(edges[0])]
    index = 0
    while index < len(values):
        group = low_group if edges[index] < 2.5 else high_group
        stop = min(index + group, len(values))
        out_values.append(float(np.sum(values[index:stop])))
        out_variances.append(float(np.sum(variances[index:stop])))
        out_edges.append(float(edges[stop]))
        index = stop
    return np.asarray(out_values), np.asarray(out_variances), np.asarray(out_edges)


def density(values: np.ndarray, variances: np.ndarray, edges: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    widths = np.diff(edges)
    return values / widths, variances / widths**2


def tail(values: np.ndarray, edges: np.ndarray, low: float = 6.0) -> float:
    first = max(0, min(len(values) - 1, int(np.searchsorted(edges, low, side="right") - 1)))
    return float(np.sum(values[first:]))


def build_stack(tight_raw: np.ndarray, tight_var_raw: np.ndarray, non_tight_raw: np.ndarray, non_tight_var_raw: np.ndarray, signal_raw: np.ndarray, signal_var_raw: np.ndarray, raw_edges: np.ndarray, low_group: int, high_group: int) -> dict[str, np.ndarray | float]:
    tight, tight_var, edges = rebin_variable(tight_raw, tight_var_raw, raw_edges, low_group, high_group)
    non_tight, non_tight_var, _ = rebin_variable(non_tight_raw, non_tight_var_raw, raw_edges, low_group, high_group)
    signal, signal_var, _ = rebin_variable(signal_raw, signal_var_raw, raw_edges, low_group, high_group)
    tight, tight_var = density(tight, tight_var, edges)
    non_tight, non_tight_var = density(non_tight, non_tight_var, edges)
    signal, signal_var = density(signal, signal_var, edges)
    scale_background = tail(tight, edges) / tail(non_tight, edges)
    non_tight *= scale_background
    non_tight_var *= scale_background**2
    signal_target = float(np.sum(tight) - np.sum(non_tight))
    scale_signal = signal_target / float(np.sum(signal))
    signal *= scale_signal
    signal_var *= scale_signal**2
    return {
        "edges": edges,
        "tight": tight,
        "tight_err": np.sqrt(np.clip(tight_var, 0.0, None)),
        "background": non_tight,
        "signal": signal,
        "background_tail_scale": scale_background,
        "signal_scale": scale_signal,
    }


def pp_panel(pp_data: Path, pp_signal: Path) -> tuple[dict[str, np.ndarray | float], dict[str, object]]:
    with uproot.open(pp_data) as data, uproot.open(pp_signal) as signal:
        tight_keys = [f"PPG12_scaledtrigger30/h_Eiso_tight_pT_{pt}" for pt in PP_BINS]
        non_tight_keys = [f"PPG12_scaledtrigger30/h_Eiso_nonTight_pT_{pt}" for pt in PP_BINS]
        signal_keys = [f"SIM/h_Eiso_tight_pT_{pt}" for pt in PP_BINS]
        tight, tight_var, edges = sum_histograms(data, tight_keys)
        non_tight, non_tight_var, _ = sum_histograms(data, non_tight_keys)
        sig, sig_var, _ = sum_histograms(signal, signal_keys)
    return build_stack(tight, tight_var, non_tight, non_tight_var, sig, sig_var, edges, 5, 10), {"tight": tight_keys, "non_tight": non_tight_keys, "signal": signal_keys}


def auau_panel(cent: str) -> tuple[dict[str, np.ndarray | float], dict[str, object]]:
    data_tight, data_nontight, signal_tight = [], [], []
    with uproot.open(THE100_DATA) as data, uproot.open(THE100_SIGNAL) as signal:
        for pt in AUAU_BINS:
            prefix_data = f"{DATA_TOP}/h_Eiso_ABCD_{{region}}_isoR40_isSliding_pT_{pt}_cent_{cent}"
            prefix_sim = f"{SIM_TOP}/h_Eiso_ABCD_{{region}}_isoR40_isSliding_pT_{pt}_cent_{cent}"
            data_tight.extend([prefix_data.format(region="A"), prefix_data.format(region="B")])
            data_nontight.extend([prefix_data.format(region="C"), prefix_data.format(region="D")])
            signal_tight.extend([prefix_sim.format(region="A"), prefix_sim.format(region="B")])
        tight, tight_var, edges = sum_histograms(data, data_tight)
        non_tight, non_tight_var, _ = sum_histograms(data, data_nontight)
        sig, sig_var, _ = sum_histograms(signal, signal_tight)
    return build_stack(tight, tight_var, non_tight, non_tight_var, sig, sig_var, edges, 5, 10), {"tight_A_plus_B": data_tight, "non_tight_C_plus_D": data_nontight, "signal_A_plus_B": signal_tight}


def draw_band(ax: plt.Axes, edges: np.ndarray, lower: np.ndarray, upper: np.ndarray, face: str, edge: str, label: str) -> None:
    x = np.ravel(np.column_stack([edges[:-1], edges[1:]]))
    ax.fill_between(x, np.repeat(lower, 2), np.repeat(upper, 2), color=face, alpha=0.70, linewidth=0, label=label)
    ax.stairs(upper, edges, baseline=lower, color=edge, linewidth=1.05)


def draw_panel(ax: plt.Axes, result: dict[str, np.ndarray | float], title: str, collision: str, pt: str, *, pp: bool = False, xmin: float = -2.0, xmax: float = 15.0) -> None:
    edges = np.asarray(result["edges"])
    tight = np.asarray(result["tight"])
    errors = np.asarray(result["tight_err"])
    bkg = np.asarray(result["background"])
    signal = np.asarray(result["signal"])
    centres = 0.5 * (edges[:-1] + edges[1:])
    draw_band(ax, edges, np.zeros_like(bkg), bkg, RED_FILL, RED_EDGE, "non-tight data control")
    draw_band(ax, edges, bkg, bkg + signal, BLUE_FILL, BLUE_EDGE, "Signal MC")
    ax.errorbar(centres, tight, yerr=errors, fmt="o", color=BLACK, ecolor=BLACK, markersize=3.8, markeredgecolor="white", markeredgewidth=0.45, elinewidth=0.85, zorder=5, label="tight-photon data")
    ymax = max(float(np.max(tight + errors)), float(np.max(bkg + signal)))
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(0, ymax * 1.24)
    ax.set_title(title, fontsize=17.0, fontweight="bold", pad=7, color=INK)
    ax.set_xlabel(r"$E_T^{\mathrm{iso,reco}}\;[\mathrm{GeV}]$", fontsize=12.5)
    ax.tick_params(which="both", direction="in", top=True, right=True, labelsize=10.5, length=5)
    ax.tick_params(which="minor", length=2.6)
    ax.minorticks_on()
    if pp:
        ax.text(0.94, 0.92, r"$\it{\bf{sPHENIX}}$ Internal", transform=ax.transAxes, ha="right", va="top", fontsize=11.0)
        ax.text(0.94, 0.82, collision, transform=ax.transAxes, ha="right", va="top", fontsize=10.0)
        ax.text(0.94, 0.73, pt, transform=ax.transAxes, ha="right", va="top", fontsize=10.0)
        ax.text(0.94, 0.64, r"$|\eta^\gamma|<0.7$", transform=ax.transAxes, ha="right", va="top", fontsize=10.0)
    else:
        ax.text(0.035, 0.945, r"$\it{\bf{sPHENIX}}$ Internal", transform=ax.transAxes, ha="left", va="top", fontsize=10.5)
        ax.text(0.035, 0.877, collision, transform=ax.transAxes, ha="left", va="top", fontsize=9.4)
        ax.text(0.035, 0.815, pt, transform=ax.transAxes, ha="left", va="top", fontsize=9.4)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--auau-xmin", type=float, default=-20.0, help="Left edge of the AuAu panels (symmetric with --auau-xmax).")
    parser.add_argument("--auau-xmax", type=float, default=20.0, help="Right edge of the AuAu panels; symmetric with --auau-xmin so the full isolation tail is shown and no data is hidden.")
    parser.add_argument("--variant-suffix", default="", help="Suffix appended to the output filenames, e.g. '_fulltail'.")
    args = parser.parse_args()

    OUTDIR.mkdir(parents=True, exist_ok=True)
    png = PNG if not args.variant_suffix else PNG.with_name(PNG.stem + args.variant_suffix + PNG.suffix)
    pp_data = current_root("pp_data_merged")
    pp_signal = current_root("pp_sim_photonjet_merged")
    pp, pp_keys = pp_panel(pp_data, pp_signal)
    auau = []
    for cent, label in CENTRALITIES:
        result, keys = auau_panel(cent)
        auau.append((label, result, keys))

    plt.rcParams.update({"font.family": "serif", "font.serif": ["Times New Roman", "Times", "DejaVu Serif"], "mathtext.fontset": "dejavuserif", "axes.linewidth": 1.15})
    fig = plt.figure(figsize=(16, 9), dpi=160)
    fig.patch.set_facecolor("white")
    fig.text(0.05, 0.953, "Isolation comparisons in p+p and Au+Au", ha="left", va="top", fontsize=34.5, fontweight="bold", color=INK)
    fig.text(0.057, 0.871, "▶", ha="center", va="top", fontsize=15.5, color=BLUE_BULLET, fontfamily="DejaVu Sans")
    fig.text(0.073, 0.871, "Tail-normalized data control and prompt-photon simulation across p+p and Au+Au centrality.", ha="left", va="top", fontsize=20.0, color=INK)
    fig.legend(handles=[Line2D([0],[0], marker="o", color=BLACK, markerfacecolor=BLACK, linestyle="", label="Tight-photon data"), Rectangle((0,0),1,1, facecolor=RED_FILL, edgecolor=RED_EDGE, label="Non-tight data control"), Rectangle((0,0),1,1, facecolor=BLUE_FILL, edgecolor=BLUE_EDGE, label="Signal MC")], loc="upper center", bbox_to_anchor=(0.5,0.836), ncol=3, frameon=False, fontsize=20.0, handlelength=1.7, columnspacing=3.0, handletextpad=0.7)
    axes = fig.subplots(1, 4, gridspec_kw={"left":0.06,"right":0.982,"bottom":0.245,"top":0.735,"wspace":0.28})
    draw_panel(axes[0], pp, r"$p{+}p$", r"$p{+}p\;\sqrt{s}=200\;\mathrm{GeV}$", r"$16<E_T^\gamma<35\;\mathrm{GeV}$", pp=True)
    for ax, (label, result, _) in zip(axes[1:], auau):
        draw_panel(ax, result, f"Au+Au {label}", r"$\mathrm{Au+Au}\;\sqrt{s_{NN}}=200\;\mathrm{GeV}$", r"$15<E_T^\gamma<35\;\mathrm{GeV}$", xmin=args.auau_xmin, xmax=args.auau_xmax)
    axes[0].set_ylabel("Counts / Bin Width", fontsize=14.0)
    fig.text(0.5, 0.181, "Prompt-photon simulation: Photon5+10+20 in p+p; embedded Photon12+20 in Au+Au.", ha="center", va="top", fontsize=18.5, color=INK)
    fig.text(0.5, 0.117, "Au+Au sliding isolation cut:", ha="center", va="top", fontsize=17.0, color=INK, fontweight="bold")
    fig.text(0.5, 0.069, r"$E_T^{\mathrm{iso}}(R=0.4)<\left(7.57-0.0658\,c\right)\ \mathrm{GeV},\qquad c=\mathrm{centrality\ percentile}$", ha="center", va="top", fontsize=16.7, color=INK)
    fig.savefig(png, dpi=160)
    plt.close(fig)

    manifest = {
        "schema": "THE100_COMPLEMENT_ISOLATION_TEMPLATE_SLIDE_V1",
        "png": str(png),
        "layout_nodes": str(LAYOUT),
        "auau_xrange": [args.auau_xmin, args.auau_xmax],
        "pp": {"data_root": str(pp_data), "data_sha256": sha256(pp_data), "signal_root": str(pp_signal), "signal_sha256": sha256(pp_signal), "pt_bins": list(PP_BINS), "keys": pp_keys},
        "auau": {"campaign": "the100_auau_dualview_20260714", "view": "unrestricted BDT complement", "isolation_view": "stored isoR40_isSliding ABCD histograms", "data_root": str(THE100_DATA), "data_sha256": sha256(THE100_DATA), "signal_root": str(THE100_SIGNAL), "signal_sha256": sha256(THE100_SIGNAL), "inclusive_triplet_member_remote_root": THE100_INCLUSIVE_REMOTE, "pt_bins": list(AUAU_BINS), "centralities": {label: keys for label, _, keys in auau}, "component_definition": "tight data A+B; non-tight data C+D; photon embedding A+B"},
        "normalization": "PPG12-style variable rebinning, bin-width scaling, tail-match non-tight above 6 GeV, and photon-embedding residual normalization.",
        "caveat": "THE-100 is a completed historical baseline using its pre-THE-111 model/reconstruction contract. It is appropriate for the requested clean JSTG baseline view but is not a claim about the pending corrected canonical production.",
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n")
    layout = {
        "canvas": [0, 0, 2560, 1440],
        "title_axis_x": 128,
        "title_axis_tolerance_px": 14,
        "minimum_audience_font_px": 26,
        "minimum_title_font_px": 60,
        "nodes": [
            {"kind": "text", "name": "slide title", "role": "title", "bbox": [128, 50, 2440, 136], "font_px": 69, "title_axis_align": "left"},
            {"kind": "text", "name": "subtitle", "role": "subtitle", "bbox": [130, 178, 2320, 227], "font_px": 38, "title_axis_align": "left"},
            {"kind": "legend", "name": "shared legend", "role": "legend", "bbox": [460, 235, 2100, 300], "font_px": 33},
            {"kind": "plot", "name": "pp panel frame", "role": "plot", "bbox": [128, 382, 822, 1080], "title_axis_align": "left"},
            {"kind": "plot", "name": "AuAu 0-20 panel", "role": "plot", "bbox": [991, 353, 1285, 1190]},
            {"kind": "plot", "name": "AuAu 20-50 panel", "role": "plot", "bbox": [1395, 353, 1689, 1190]},
            {"kind": "plot", "name": "AuAu 50-80 panel", "role": "plot", "bbox": [1799, 353, 2514, 1190]},
            {"kind": "text", "name": "simulation definition", "role": "caption", "bbox": [128, 1170, 1760, 1230], "font_px": 36, "title_axis_align": "left"},
            {"kind": "text", "name": "isolation definition", "role": "caption", "bbox": [128, 1270, 760, 1335], "font_px": 36, "title_axis_align": "left"},
            {"kind": "equation", "name": "canonical sliding isolation window", "role": "equation", "bbox": [852, 1260, 2439, 1367], "font_px": 37},
        ],
    }
    LAYOUT.write_text(json.dumps(layout, indent=2) + "\n")
    SPEAKER.write_text("# JSTG Slide Script - Isolation-energy templates\n\nHere I show the reconstructed isolation-energy distribution for the selected tight-photon sample. The black points are data. The red component comes from the broad non-tight photon control region, normalized in the high-isolation tail, and the blue component is prompt-photon embedding normalized to the remaining tight yield.\n\nThe point is that this decomposition remains well behaved from central to peripheral Au+Au events, using the same centrality-dependent sliding isolation definition throughout. This is the completed baseline used to motivate the next photon-selection update.\n")
    print(png)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
