#!/usr/bin/env python3
"""Build current THE-42 AuAu shower-shape and isolation debug plots."""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


REPO = Path(__file__).resolve().parents[3]
DEFAULT_DATA_JSON = (
    REPO
    / "dataOutput/ppg12TableQA/THE42_current_auau_tableqa_20260614/debug/current_auau_completed_run_shape_iso_hists.json"
)
DEFAULT_SIGNAL_ROOT = REPO / "InputFiles/the42_current_auau_tableqa/RecoilJets_embeddedPhoton12plus20_MERGED.root"
DEFAULT_INCLUSIVE_ROOT = REPO / "InputFiles/the42_current_auau_tableqa/RecoilJets_embeddedJet12plus20plus30plus40_MERGED.root"
DEFAULT_OUTDIR = REPO / "dataOutput/ppg12TableQA/THE42_current_auau_tableqa_20260614/debug"

CENTRALITIES = [
    ("cent0_20", "0-20%"),
    ("cent20_50", "20-50%"),
    ("cent50_80", "50-80%"),
]
STAGES = [("cut0", "Before preselection"), ("cut1", "After preselection"), ("cut2", "After tight ID")]
SHAPE_VARS = [
    ("e11_to_e33", r"$E_{11}/E_{33}$", (0.0, 1.0), 4),
    ("weta_cogx", r"$w_{\eta}^{cogX}$", (0.0, 0.6), 2),
    ("e32_to_e35", r"$E_{32}/E_{35}$", (0.0, 1.0), 4),
    ("bdt", "BDT score", (0.0, 1.0), 4),
]
ISO_MODES = [
    ("isoR40_all", "All candidates"),
    ("isoR40_tight", "Tight ID candidates"),
]
SAMPLE_STYLE = {
    "Data": {"color": "#111827", "marker": "o"},
    "Signal MC": {"color": "#C22F2F"},
    "Inclusive MC": {"color": "#2E63D4"},
}


@dataclass
class HistArrays:
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
            "axes.linewidth": 0.9,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
        }
    )


def rebin(x: np.ndarray, y: np.ndarray, e: np.ndarray, factor: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if factor <= 1 or len(x) < factor:
        return x, y, e
    n = len(x) // factor
    trim = n * factor
    return (
        x[:trim].reshape(n, factor).mean(axis=1),
        y[:trim].reshape(n, factor).sum(axis=1),
        np.sqrt((e[:trim].reshape(n, factor) ** 2).sum(axis=1)),
    )


def norm_payload(payload: dict, *, xlim: tuple[float, float], rebin_factor: int, source: str) -> HistArrays:
    x = np.asarray(payload["x"], dtype=float)
    y = np.asarray(payload["y"], dtype=float)
    e = np.asarray(payload.get("e", np.sqrt(np.maximum(y, 0.0))), dtype=float)
    x, y, e = rebin(x, y, e, rebin_factor)
    mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(e) & (x >= xlim[0]) & (x <= xlim[1])
    x, y, e = x[mask], y[mask], e[mask]
    integral = float(y.sum())
    if integral > 0:
        y = y / integral
        e = e / integral
    return HistArrays(x=x, y=y, e=e, integral=integral, source=source)


def open_root(path: Path):
    import ROOT

    ROOT.gROOT.SetBatch(True)
    ROOT.TH1.AddDirectory(False)
    f = ROOT.TFile.Open(str(path), "READ")
    if not f or f.IsZombie():
        raise OSError(path)
    return f


def root_hist_payload(root: Path, hist_path: str) -> dict:
    f = open_root(root)
    try:
        h = f.Get(hist_path)
        if not h:
            return {"x": [], "y": [], "e": []}
        h = h.Clone("tmp_clone")
        h.SetDirectory(0)
    finally:
        f.Close()
    n = h.GetNbinsX()
    return {
        "x": [float(h.GetXaxis().GetBinCenter(i)) for i in range(1, n + 1)],
        "y": [float(h.GetBinContent(i)) for i in range(1, n + 1)],
        "e": [float(h.GetBinError(i)) for i in range(1, n + 1)],
    }


def add_payload(a: dict, b: dict) -> dict:
    if not a.get("x"):
        return {"x": list(b.get("x", [])), "y": list(b.get("y", [])), "e": list(b.get("e", []))}
    if not b.get("x"):
        return a
    ay = np.asarray(a["y"], dtype=float)
    ae = np.asarray(a["e"], dtype=float)
    by = np.asarray(b["y"], dtype=float)
    be = np.asarray(b["e"], dtype=float)
    return {"x": a["x"], "y": (ay + by).tolist(), "e": np.sqrt(ae * ae + be * be).tolist()}


def mc_shape_payload(root: Path, variable: str, cent: str, stage: str) -> dict:
    return root_hist_payload(root, f"SIM/h1d_{variable}_eta0_pt1535_{cent}_{stage}")


PT_BINS = ["15_18", "18_20", "20_22", "22_24", "24_26", "26_28", "28_30", "30_35"]
CENT_FINE = {
    "cent0_20": ["0_10", "10_20"],
    "cent20_50": ["20_30", "30_40", "40_50"],
    "cent50_80": ["50_60", "60_80"],
}


def mc_iso_payload(root: Path, mode: str, cent: str) -> dict:
    if mode == "isoR40_all":
        pattern = "SIM/h_Eiso_isoR40_pT_{pt}_cent_{cent}"
    elif mode == "isoR40_tight":
        pattern = "SIM/h_Eiso_tight_isoR40_pT_{pt}_cent_{cent}"
    else:
        raise ValueError(mode)
    out = {"x": [], "y": [], "e": []}
    for pt in PT_BINS:
        for fine in CENT_FINE[cent]:
            out = add_payload(out, root_hist_payload(root, pattern.format(pt=pt, cent=fine)))
    return out


def draw_hist(ax, arr: HistArrays, label: str) -> None:
    style = SAMPLE_STYLE[label]
    color = style["color"]
    if label == "Data":
        ax.errorbar(
            arr.x,
            arr.y,
            yerr=arr.e,
            fmt="o",
            color=color,
            mfc=color,
            mec="white",
            mew=0.55,
            ms=4.2,
            elinewidth=0.7,
            capsize=1.3,
            label=label,
        )
    else:
        ax.step(arr.x, arr.y, where="mid", color=color, lw=1.9, label=label)


def fraction_below(payload: dict, threshold: float) -> float:
    x = np.asarray(payload.get("x", []), dtype=float)
    y = np.asarray(payload.get("y", []), dtype=float)
    total = float(y.sum())
    if total <= 0:
        return float("nan")
    return float(y[x < threshold].sum() / total)


def make_shape_table(data: dict, signal_root: Path, inclusive_root: Path, outdir: Path) -> list[dict]:
    rows = []
    fig, axes = plt.subplots(len(SHAPE_VARS), len(CENTRALITIES), figsize=(17.0, 12.0), constrained_layout=False)
    fig.subplots_adjust(left=0.065, right=0.990, top=0.900, bottom=0.075, wspace=0.18, hspace=0.28)
    for r, (var, var_label, xlim, rb) in enumerate(SHAPE_VARS):
        for c, (cent, cent_label) in enumerate(CENTRALITIES):
            ax = axes[r, c]
            curves = {
                "Data": norm_payload(data["shape_hists"][var][cent]["cut1"], xlim=xlim, rebin_factor=rb, source="data aggregate"),
                "Signal MC": norm_payload(mc_shape_payload(signal_root, var, cent, "cut1"), xlim=xlim, rebin_factor=rb, source="signal"),
                "Inclusive MC": norm_payload(mc_shape_payload(inclusive_root, var, cent, "cut1"), xlim=xlim, rebin_factor=rb, source="inclusive"),
            }
            ymax = 0.0
            for label, arr in curves.items():
                draw_hist(ax, arr, label)
                if len(arr.y):
                    ymax = max(ymax, float(np.nanmax(arr.y + arr.e)))
            ax.set_xlim(*xlim)
            ax.set_ylim(0, max(0.025, ymax * 1.18))
            ax.grid(True, axis="y", color="#E5E7EB", lw=0.6)
            if r == 0:
                ax.set_title(cent_label, fontsize=16, fontweight="bold")
            if c == 0:
                ax.set_ylabel(f"{var_label}\nnormalized", fontsize=13)
            if r == len(SHAPE_VARS) - 1:
                ax.set_xlabel(var_label, fontsize=13)
            ax.tick_params(labelsize=10)
            for label, payload in (
                ("Data", data["shape_hists"][var][cent]["cut1"]),
                ("Signal MC", mc_shape_payload(signal_root, var, cent, "cut1")),
                ("Inclusive MC", mc_shape_payload(inclusive_root, var, cent, "cut1")),
            ):
                rows.append(
                    {
                        "variable": var,
                        "centrality": cent,
                        "stage": "cut1",
                        "sample": label,
                        "integral": float(np.sum(payload.get("y", []))),
                        "frac_x_lt_0p05": fraction_below(payload, 0.05),
                        "frac_x_lt_0p10": fraction_below(payload, 0.10),
                    }
                )
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper right", bbox_to_anchor=(0.986, 0.965), frameon=False, ncol=3, fontsize=13)
    fig.text(0.060, 0.965, "Current AuAu table-QA shower-shape sanity: after preselection", fontsize=21, fontweight="bold", ha="left", va="top")
    fig.text(0.060, 0.930, "Data uses completed-run aggregate from the current campaign; MC uses matched table-QA merged signal/inclusive outputs.", fontsize=13.5, ha="left", va="top")
    out = outdir / "current_auau_shape_table_after_preselection.png"
    fig.savefig(out, dpi=180)
    plt.close(fig)
    return rows


def make_e11_flow(data: dict, signal_root: Path, inclusive_root: Path, outdir: Path) -> None:
    fig, axes = plt.subplots(len(CENTRALITIES), len(STAGES), figsize=(17.0, 10.5), constrained_layout=False)
    fig.subplots_adjust(left=0.070, right=0.990, top=0.875, bottom=0.075, wspace=0.15, hspace=0.24)
    for r, (cent, cent_label) in enumerate(CENTRALITIES):
        for c, (stage, stage_label) in enumerate(STAGES):
            ax = axes[r, c]
            curves = {
                "Data": norm_payload(data["shape_hists"]["e11_to_e33"][cent][stage], xlim=(0, 1), rebin_factor=2, source="data aggregate"),
                "Signal MC": norm_payload(mc_shape_payload(signal_root, "e11_to_e33", cent, stage), xlim=(0, 1), rebin_factor=2, source="signal"),
                "Inclusive MC": norm_payload(mc_shape_payload(inclusive_root, "e11_to_e33", cent, stage), xlim=(0, 1), rebin_factor=2, source="inclusive"),
            }
            ymax = 0
            for label, arr in curves.items():
                draw_hist(ax, arr, label)
                if len(arr.y):
                    ymax = max(ymax, float(np.nanmax(arr.y + arr.e)))
            ax.axvspan(0, 0.05, color="#FDE68A", alpha=0.20, lw=0)
            ax.set_xlim(0, 1)
            ax.set_ylim(0, max(0.03, ymax * 1.18))
            ax.grid(True, axis="y", color="#E5E7EB", lw=0.6)
            if r == 0:
                ax.set_title(stage_label, fontsize=15, fontweight="bold")
            if c == 0:
                ax.set_ylabel(f"{cent_label}\nnormalized", fontsize=13)
            if r == len(CENTRALITIES) - 1:
                ax.set_xlabel(r"$E_{11}/E_{33}$", fontsize=13)
            ax.tick_params(labelsize=10)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper right", bbox_to_anchor=(0.986, 0.945), frameon=False, ncol=3, fontsize=13)
    fig.text(0.060, 0.960, r"Current AuAu $E_{11}/E_{33}$ flow with low-edge region highlighted", fontsize=21, fontweight="bold", ha="left", va="top")
    fig.text(0.060, 0.920, "Yellow band marks x < 0.05. A real low-edge excess should persist by sample/stage; a hot-tower-like artifact would also distort data sharply.", fontsize=13.5, ha="left", va="top")
    fig.savefig(outdir / "current_auau_e11_low_edge_flow.png", dpi=180)
    plt.close(fig)


def make_iso_table(data: dict, signal_root: Path, inclusive_root: Path, outdir: Path) -> None:
    fig, axes = plt.subplots(len(ISO_MODES), len(CENTRALITIES), figsize=(17.0, 7.6), constrained_layout=False)
    fig.subplots_adjust(left=0.070, right=0.990, top=0.845, bottom=0.100, wspace=0.17, hspace=0.28)
    for r, (mode, mode_label) in enumerate(ISO_MODES):
        for c, (cent, cent_label) in enumerate(CENTRALITIES):
            ax = axes[r, c]
            curves = {
                "Data": norm_payload(data["iso_hists"][mode][cent], xlim=(-20, 40), rebin_factor=3, source="data aggregate"),
                "Signal MC": norm_payload(mc_iso_payload(signal_root, mode, cent), xlim=(-20, 40), rebin_factor=3, source="signal"),
                "Inclusive MC": norm_payload(mc_iso_payload(inclusive_root, mode, cent), xlim=(-20, 40), rebin_factor=3, source="inclusive"),
            }
            ymax = 0
            for label, arr in curves.items():
                draw_hist(ax, arr, label)
                if len(arr.y):
                    ymax = max(ymax, float(np.nanmax(arr.y + arr.e)))
            ax.axvline(0, color="#6B7280", lw=0.85, ls="--")
            ax.set_xlim(-20, 40)
            ax.set_ylim(0, max(0.025, ymax * 1.18))
            ax.grid(True, axis="y", color="#E5E7EB", lw=0.6)
            if r == 0:
                ax.set_title(cent_label, fontsize=16, fontweight="bold")
            if c == 0:
                ax.set_ylabel(f"{mode_label}\nnormalized", fontsize=13)
            if r == len(ISO_MODES) - 1:
                ax.set_xlabel(r"$E_{\mathrm{iso}}$ (GeV), R=0.4 sliding", fontsize=13)
            ax.tick_params(labelsize=10)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper right", bbox_to_anchor=(0.986, 0.945), frameon=False, ncol=3, fontsize=13)
    fig.text(0.060, 0.960, "Current AuAu isolation sanity: R=0.4 sliding isolation", fontsize=21, fontweight="bold", ha="left", va="top")
    fig.text(0.060, 0.913, "Data aggregate is from completed current-campaign chunks; curves are normalized shapes over -20 to 40 GeV.", fontsize=13.5, ha="left", va="top")
    fig.savefig(outdir / "current_auau_isolation_overlay.png", dpi=180)
    plt.close(fig)


def write_metrics(rows: list[dict], outdir: Path) -> None:
    csv = outdir / "current_auau_shape_low_edge_metrics.csv"
    keys = ["variable", "centrality", "stage", "sample", "integral", "frac_x_lt_0p05", "frac_x_lt_0p10"]
    with csv.open("w") as f:
        f.write(",".join(keys) + "\n")
        for row in rows:
            f.write(",".join(str(row.get(k, "")) for k in keys) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data-json", type=Path, default=DEFAULT_DATA_JSON)
    ap.add_argument("--signal-root", type=Path, default=DEFAULT_SIGNAL_ROOT)
    ap.add_argument("--inclusive-root", type=Path, default=DEFAULT_INCLUSIVE_ROOT)
    ap.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    args = ap.parse_args()
    setup_style()
    args.outdir.mkdir(parents=True, exist_ok=True)
    data = json.loads(args.data_json.read_text())
    rows = make_shape_table(data, args.signal_root, args.inclusive_root, args.outdir)
    make_e11_flow(data, args.signal_root, args.inclusive_root, args.outdir)
    make_iso_table(data, args.signal_root, args.inclusive_root, args.outdir)
    write_metrics(rows, args.outdir)
    manifest = {
        "schema": "THE42_CURRENT_AUAU_SHAPE_ISO_DEBUG_V1",
        "data_json": str(args.data_json),
        "signal_root": str(args.signal_root),
        "inclusive_root": str(args.inclusive_root),
        "outputs": [
            str(args.outdir / "current_auau_shape_table_after_preselection.png"),
            str(args.outdir / "current_auau_e11_low_edge_flow.png"),
            str(args.outdir / "current_auau_isolation_overlay.png"),
            str(args.outdir / "current_auau_shape_low_edge_metrics.csv"),
        ],
        "data_metadata": data.get("metadata", {}),
    }
    (args.outdir / "current_auau_shape_iso_debug_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    for item in manifest["outputs"]:
        print(item)
    print(str(args.outdir / "current_auau_shape_iso_debug_manifest.json"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
