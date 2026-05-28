#!/usr/bin/env python3
"""Make centrality-sliced scaled-trigger summary panels from merged ROOT output."""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import ROOT  # type: ignore  # noqa: E402


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
DEFAULT_INPUT = (
    REPO
    / "InputFiles/auau25"
    / "RecoilJets_auau_ALL_jetMinPt5_7pi_8_vz60_isoR40_isSliding_baseVariant_preselectionReference_tightReference_nonTightReference_scaledTriggerCentStudy_cent0_20_50_80.root"
)
DEFAULT_OUTDIR = REPO / "dataOutput/auau/scaledTriggerCentStudy"
FIT_X_MAX = 16.0
FIT_ERROR_FLOOR = 0.04

TRIGGERS = {
    "mbd": (
        "MBD_NS_geq_2_vtx_lt_150",
        "h_maxEnergyClus_NewTriggerFilling_perRunCorrected_MBD_NS_geq_2_vtx_lt_150",
        "MBD N&S >= 2, |vz| < 150 cm",
        "black",
    ),
    "p10": (
        "Photon_10",
        "h_maxEnergyClus_NewTriggerFilling_perRunCorrected_Photon_10",
        "Photon 10, scaled",
        "#1f4aff",
    ),
    "p12": (
        "Photon_12",
        "h_maxEnergyClus_NewTriggerFilling_perRunCorrected_Photon_12",
        "Photon 12, scaled",
        "#e31a1c",
    ),
}


@dataclass
class Hist:
    values: np.ndarray
    errors: np.ndarray
    edges: np.ndarray

    @property
    def centers(self) -> np.ndarray:
        return 0.5 * (self.edges[:-1] + self.edges[1:])

    @property
    def widths(self) -> np.ndarray:
        return self.edges[1:] - self.edges[:-1]


@dataclass
class TurnOnFit:
    floor: float
    plateau: float
    x50: float
    width: float
    x90: float
    used_points: int
    loss: float

    def value(self, xval: float) -> float:
        arg = max(min(-(xval - self.x50) / self.width, 60.0), -60.0)
        return self.floor + (self.plateau - self.floor) / (1.0 + math.exp(arg))


def parse_cent_bins(text: str) -> list[tuple[str, str]]:
    vals = [int(v.strip()) for v in text.replace(":", ",").split(",") if v.strip()]
    if len(vals) < 2:
        raise ValueError("--cent-edges needs at least two comma-separated edges")
    if vals != sorted(vals) or len(set(vals)) != len(vals):
        raise ValueError("--cent-edges must be strictly increasing")
    return [(f"cent{lo}_{hi}", f"{lo}-{hi}%") for lo, hi in zip(vals[:-1], vals[1:])]


def read_hist(root_file: "ROOT.TFile", trigger: str, suffix: str) -> Hist:
    directory, prefix, _, _ = TRIGGERS[trigger]
    hist_name = f"{prefix}_{suffix}" if suffix else prefix
    path = f"{directory}/{hist_name}"
    obj = root_file.Get(path)
    if not obj:
        raise KeyError(f"missing histogram: {path}")
    nbins = obj.GetNbinsX()
    axis = obj.GetXaxis()
    values = np.array([obj.GetBinContent(i) for i in range(1, nbins + 1)], dtype=float)
    errors = np.array([obj.GetBinError(i) for i in range(1, nbins + 1)], dtype=float)
    edges = np.array([axis.GetBinLowEdge(i) for i in range(1, nbins + 1)] + [axis.GetBinUpEdge(nbins)], dtype=float)
    return Hist(values, errors, edges)


def ratio_hist(num: Hist, den: Hist) -> Hist:
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.divide(num.values, den.values, out=np.full_like(num.values, np.nan), where=den.values > 0)
        rel_num = np.divide(num.errors, num.values, out=np.zeros_like(num.values), where=num.values > 0)
        rel_den = np.divide(den.errors, den.values, out=np.zeros_like(den.values), where=den.values > 0)
        err = ratio * np.sqrt(rel_num * rel_num + rel_den * rel_den)
    err = np.where(np.isfinite(err), err, np.nan)
    return Hist(ratio, err, den.edges)


def integral(hist: Hist, lo: float, hi: float) -> float:
    mask = (hist.centers >= lo) & (hist.centers < hi)
    return float(np.nansum(hist.values[mask]))


def ratio_err(num: float, den: float) -> tuple[float, float]:
    if den <= 0:
        return math.nan, math.nan
    value = num / den
    if num <= 0:
        return value, math.sqrt(max(num, 1.0)) / den
    return value, value * math.sqrt(1.0 / num + 1.0 / den)


def fit_turnon(num: Hist, den: Hist) -> TurnOnFit:
    data: list[tuple[float, float, float]] = []
    for xval, nval, dval in zip(den.centers, num.values, den.values):
        if not (1.0 <= xval <= FIT_X_MAX) or dval <= 0:
            continue
        value, err = ratio_err(float(nval), float(dval))
        if math.isfinite(value):
            data.append((float(xval), value, max(err, FIT_ERROR_FLOOR)))

    def model(params: list[float], xval: float) -> float:
        floor, plateau, x50, width = params
        arg = max(min(-(xval - x50) / width, 60.0), -60.0)
        return floor + (plateau - floor) / (1.0 + math.exp(arg))

    def loss(params: list[float]) -> float:
        floor, plateau, x50, width = params
        if not (0.0 <= floor <= 0.10 and 0.50 <= plateau <= 1.15 and 3.0 <= x50 <= 10.5 and 0.35 <= width <= 5.0):
            return 1.0e99
        if plateau <= floor + 0.20:
            return 1.0e99
        total = 0.0
        for xval, yval, err in data:
            residual = (model(params, xval) - yval) / err
            total += 2.0 * (math.sqrt(1.0 + residual * residual) - 1.0)
        return total

    starts = [
        [floor, plateau, x50, width]
        for floor in [0.0, 0.004, 0.01, 0.02, 0.05]
        for plateau in [0.70, 0.82, 0.90, 0.96, 1.02]
        for x50 in [5.0, 5.8, 6.6, 7.4, 8.2]
        for width in [0.7, 1.1, 1.6, 2.2, 3.0]
    ]
    best_params = starts[0]
    best_loss = loss(best_params)
    for params in starts[1:]:
        current = loss(params)
        if current < best_loss:
            best_loss = current
            best_params = params

    steps = [0.010, 0.015, 0.12, 0.12]
    for _ in range(100):
        improved = False
        for idx, step in enumerate(steps):
            for sign in (-1.0, 1.0):
                trial = best_params[:]
                trial[idx] += sign * step
                current = loss(trial)
                if current < best_loss:
                    best_loss = current
                    best_params = trial
                    improved = True
        if not improved:
            steps = [step * 0.70 for step in steps]

    floor, plateau, x50, width = best_params
    x90 = x50 + width * math.log(0.9 / 0.1)
    return TurnOnFit(floor, plateau, x50, width, x90, len(data), best_loss)


def plot_panel(out_path: Path, label: str, hists: dict[str, Hist]) -> dict[str, float]:
    ratios = {"p10": ratio_hist(hists["p10"], hists["mbd"]), "p12": ratio_hist(hists["p12"], hists["mbd"])}
    fits = {"p10": fit_turnon(hists["p10"], hists["mbd"]), "p12": fit_turnon(hists["p12"], hists["mbd"])}

    fig, (ax_l, ax_r) = plt.subplots(1, 2, figsize=(15.6, 7.0), dpi=180)
    fig.subplots_adjust(left=0.075, right=0.985, bottom=0.12, top=0.90, wspace=0.25)

    for key in ("mbd", "p10", "p12"):
        _, _, legend, color = TRIGGERS[key]
        h = hists[key]
        ax_l.step(h.edges[:-1], h.values, where="post", lw=2.0, color=color, label=legend)

    positives = np.concatenate([h.values[h.values > 0] for h in hists.values()])
    ymin = max(1.0, float(np.nanmin(positives)) * 0.5) if positives.size else 1.0
    ymax = max(10.0, float(max(np.nanmax(h.values) for h in hists.values())) * 1.8)
    ax_l.set_yscale("log")
    ax_l.set_xlim(1.0, 20.0)
    ax_l.set_ylim(ymin, ymax)
    ax_l.set_xlabel(r"Max cluster energy [GeV], $E_{\mathrm{clus}}>1$ GeV")
    ax_l.set_ylabel("Live/scaled-corrected counts")
    ax_l.legend(frameon=False, fontsize=10, loc="upper right")
    ax_l.text(0.04, 0.95, rf"$\it{{sPHENIX}}$ Internal  Au+Au, $\sqrt{{s_{{NN}}}}=200$ GeV",
              transform=ax_l.transAxes, fontsize=10, va="top")
    ax_l.text(0.04, 0.89, f"Centrality {label}", transform=ax_l.transAxes, fontsize=10, va="top")

    xfit = np.linspace(1.0, FIT_X_MAX, 220)
    for key, marker in (("p10", "o"), ("p12", "s")):
        _, _, _, color = TRIGGERS[key]
        r = ratios[key]
        mask = np.isfinite(r.values)
        ax_r.errorbar(
            r.centers[mask],
            r.values[mask],
            yerr=r.errors[mask],
            xerr=0.5 * r.widths[mask],
            fmt=marker,
            ms=3.0,
            lw=0.9,
            capsize=1.6,
            color=color,
            label=f"{'Photon 10' if key == 'p10' else 'Photon 12'} / MBD",
        )
        ax_r.plot(xfit, [fits[key].value(float(x)) for x in xfit], color=color, lw=2.0)

    ax_r.axhline(1.0, color="0.75", lw=1.0, zorder=0)
    ax_r.set_xlim(1.0, 20.0)
    max_ratio = max(float(np.nanmax(r.values[np.isfinite(r.values)])) if np.any(np.isfinite(r.values)) else 1.0 for r in ratios.values())
    ax_r.set_ylim(0.0, min(2.2, max(1.25, max_ratio * 1.15)))
    ax_r.set_xlabel(r"Max cluster energy [GeV], $E_{\mathrm{clus}}>1$ GeV")
    ax_r.set_ylabel("Trigger / MBD")
    ax_r.legend(frameon=False, fontsize=10, loc="lower right")
    ax_r.text(0.04, 0.95, "Free-plateau robust sigmoid fit over 1-16 GeV",
              transform=ax_r.transAxes, fontsize=10, va="top")
    ax_r.text(
        0.04,
        0.88,
        f"P10 plateau={fits['p10'].plateau:.3f}, x90={fits['p10'].x90:.2f} GeV\n"
        f"P12 plateau={fits['p12'].plateau:.3f}, x90={fits['p12'].x90:.2f} GeV",
        transform=ax_r.transAxes,
        fontsize=10,
        va="top",
    )

    for ax in (ax_l, ax_r):
        ax.tick_params(direction="in", top=True, right=True)

    fig.savefig(out_path)
    plt.close(fig)

    mbd_tail = integral(hists["mbd"], 15.0, 20.0)
    p10_tail = integral(hists["p10"], 15.0, 20.0)
    p12_tail = integral(hists["p12"], 15.0, 20.0)
    return {
        "mbd_total1_20": integral(hists["mbd"], 1.0, 20.0),
        "mbd_tail15_20": mbd_tail,
        "p10_tail15_20": p10_tail,
        "p12_tail15_20": p12_tail,
        "r10_tail15_20": p10_tail / mbd_tail if mbd_tail > 0 else math.nan,
        "r12_tail15_20": p12_tail / mbd_tail if mbd_tail > 0 else math.nan,
        "p10_plateau_fit": fits["p10"].plateau,
        "p12_plateau_fit": fits["p12"].plateau,
        "p10_x90_fit": fits["p10"].x90,
        "p12_x90_fit": fits["p12"].x90,
        "p10_fit_points": fits["p10"].used_points,
        "p12_fit_points": fits["p12"].used_points,
    }


def write_summary(path: Path, rows: Iterable[dict[str, object]]) -> None:
    rows = list(rows)
    if not rows:
        return
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--cent-edges", default="0,20,50,80")
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    cent_bins = parse_cent_bins(args.cent_edges)

    rows: list[dict[str, object]] = []
    ROOT.gROOT.SetBatch(True)
    root_file = ROOT.TFile.Open(str(args.input), "READ")
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"could not open ROOT input: {args.input}")
    try:
        for suffix, label in cent_bins:
            hists = {key: read_hist(root_file, key, suffix) for key in TRIGGERS}
            out_png = args.outdir / f"scaledTriggerCentStudy_{suffix}_1x2.png"
            metrics = plot_panel(out_png, label, hists)
            rows.append({"centrality": label, "suffix": suffix, "png": str(out_png), **metrics})
    finally:
        root_file.Close()

    summary_path = args.outdir / "scaledTriggerCentStudy_centrality_summary.csv"
    write_summary(summary_path, rows)
    print(f"Wrote {len(rows)} panels to {args.outdir}")
    print(f"Wrote summary CSV: {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
