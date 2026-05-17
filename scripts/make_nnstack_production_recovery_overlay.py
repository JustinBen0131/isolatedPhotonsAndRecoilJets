#!/usr/bin/env python3
"""Make slide-ready recovery overlays for the NN-stack WP80 production output."""

from __future__ import annotations

import csv
import math
import re
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import ROOT

import plot_target80_deep_dive as deep


ROOT.gROOT.SetBatch(True)

REPO = Path(__file__).resolve().parents[1]
OUT_DIR = REPO / "dataOutput/auauBDTMLPStackProductionPlots/nnstack_wp080_prod_20260514_200206"

NNSTACK_ROOT = (
    REPO
    / "dataOutput/combinedSimOnlyEMBEDDED/bdt_mlp_nnstack_wp080_nnstack_wp080_prod_20260514_200206"
)
NNSTACK_CFG = "preselectionReference_tightAuAuBDTMLPStack_nonTightAuAuBDTMLPStackComplement_baseVariant"
NNSTACK_SIGNAL = NNSTACK_ROOT / NNSTACK_CFG / deep.SIGNAL_REL
NNSTACK_BKG = NNSTACK_ROOT / NNSTACK_CFG / deep.BKG_REL


@dataclass(frozen=True)
class Curve:
    key: str
    label: str
    signal_path: Path
    bkg_path: Path
    color: str
    marker: str
    shift: float


def curves() -> list[Curve]:
    bdt = deep.width_variants("pt15to35", 15.0, 35.0)[1]
    return [
        Curve(
            "box",
            "Box-cuts",
            deep.REFERENCE.signal_path(),
            deep.REFERENCE.bkg_path(),
            "#111111",
            "o",
            -0.14,
        ),
        Curve(
            "bdt_wp80",
            "15-35 BDT, target 80%",
            bdt.signal_path(),
            bdt.bkg_path(),
            "#D55E00",
            "D",
            0.0,
        ),
        Curve(
            "nnstack_wp80",
            "BDT+MLP NN stack, target 80%",
            NNSTACK_SIGNAL,
            NNSTACK_BKG,
            "#6F3FB5",
            "^",
            0.14,
        ),
    ]


def setup_axes(ax, ylabel: str | None, xlabel: str | None, ylim: tuple[float, float]):
    ax.tick_params(direction="in", which="both", top=True, right=True, length=6)
    ax.minorticks_on()
    ax.grid(True, color="#D0D0D0", alpha=0.45, linewidth=0.7)
    ax.set_xlim(14.7, 35.4)
    ax.set_ylim(*ylim)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=15)
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=15)
    ax.tick_params(labelsize=12)


def draw_header(ax, y=0.905):
    ax.text(
        0.055,
        y,
        "sPHENIX",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=13,
        fontweight="bold",
        fontstyle="italic",
    )
    ax.text(0.265, y, "Internal", transform=ax.transAxes, ha="left", va="top", fontsize=13)
    ax.text(
        0.055,
        y - 0.095,
        r"Pythia overlay, $\sqrt{s_{NN}} = 200$ GeV",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=11,
    )
    ax.text(
        0.055,
        y - 0.160,
        r"$R = 0.3$, sliding isolation",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=11,
    )


def points_for(kind: str, c: Curve, cent: str) -> list[dict]:
    if kind == "recovery":
        return deep.recovery_points(c.signal_path, cent, 15.0, 35.0)
    if kind == "inclusive_tight":
        return tight_fraction_points_from_file(c.bkg_path, cent)
    raise ValueError(kind)


def sum_hist(h) -> tuple[float, float]:
    if not h:
        return 0.0, 0.0
    total = 0.0
    err2 = 0.0
    for ib in range(0, h.GetNbinsX() + 2):
        total += h.GetBinContent(ib)
        err2 += h.GetBinError(ib) ** 2
    return total, math.sqrt(err2)


def tight_fraction_points_from_file(path: Path, cent: str) -> list[dict]:
    f, d = deep.open_sim(path)
    if not d:
        return []
    pattern = re.compile(rf"^h_Eiso_tight_isoR30_pT_([0-9]+)_([0-9]+)_cent_{cent}$")
    bins: list[tuple[float, float, str]] = []
    keys = d.GetListOfKeys()
    for ik in range(keys.GetEntries()):
        name = keys.At(ik).GetName()
        match = pattern.match(name)
        if not match:
            continue
        lo = float(match.group(1))
        hi = float(match.group(2))
        if hi <= 15.0 or lo >= 35.0:
            continue
        bins.append((lo, hi, name))
    points: list[dict] = []
    for lo, hi, tight_name in sorted(bins):
        non_name = tight_name.replace("h_Eiso_tight_", "h_Eiso_nonTight_", 1)
        h_tight = d.Get(tight_name)
        h_nontight = d.Get(non_name)
        tight, e_tight = sum_hist(h_tight)
        nontight, e_nontight = sum_hist(h_nontight)
        value, err, den = deep.ratio_point(tight, nontight, e_tight, e_nontight)
        if math.isfinite(value):
            points.append({"x": 0.5 * (lo + hi), "y": value, "ey": err, "num": tight, "den": den})
    f.Close()
    return points


def integrated_tight_fraction_from_file(path: Path, cent: str) -> tuple[float, float, float]:
    pts = tight_fraction_points_from_file(path, cent)
    num = sum(p["num"] for p in pts)
    den = sum(p["den"] for p in pts)
    if den <= 0:
        return math.nan, math.nan, 0.0
    value = num / den
    err = math.sqrt(max(value * (1.0 - value), 0.0) / den)
    return value, err, den


def draw_plot(kind: str, out_name: str, ylabel: str, xlabel: str, ylim: tuple[float, float]) -> Path:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(17.0, 5.2), sharey=True)
    rows: list[dict] = []
    for ic, (cent, title) in enumerate(deep.CENTS):
        ax = axes[ic]
        for c in curves():
            pts = points_for(kind, c, cent)
            if not pts:
                continue
            x = np.array([p["x"] + c.shift for p in pts], dtype=float)
            y = np.array([p["y"] for p in pts], dtype=float)
            ey = np.array([p["ey"] for p in pts], dtype=float)
            ax.errorbar(
                x,
                y,
                yerr=ey,
                fmt=c.marker,
                linestyle="none",
                color=c.color,
                markerfacecolor=c.color,
                markeredgecolor=c.color,
                markersize=6.8,
                elinewidth=1.3,
                capsize=2.3,
                label=c.label,
                zorder=3,
            )
            for p in pts:
                rows.append(
                    {
                        "plot": kind,
                        "model_key": c.key,
                        "model_label": c.label,
                        "centrality": cent,
                        "pt_mid": p["x"],
                        "value": p["y"],
                        "error": p["ey"],
                        "numerator": p["num"],
                        "denominator": p["den"],
                        "source_signal": str(c.signal_path.relative_to(REPO)),
                        "source_background": str(c.bkg_path.relative_to(REPO)),
                    }
                )
        ax.set_title(title, fontsize=17, pad=12)
        setup_axes(ax, ylabel if ic == 0 else None, xlabel if ic == 1 else None, ylim)
        if ic == 0:
            draw_header(ax)
        if ic == 2:
            leg = ax.legend(
                loc="upper right",
                frameon=True,
                framealpha=0.94,
                fontsize=10.5,
                borderpad=0.55,
                handletextpad=0.45,
            )
            leg.get_frame().set_linewidth(0.0)
    fig.tight_layout(w_pad=1.5)
    path = OUT_DIR / out_name
    fig.savefig(path, dpi=190)
    plt.close(fig)

    csv_path = OUT_DIR / f"{Path(out_name).stem}.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()) if rows else [])
        if rows:
            writer.writeheader()
            writer.writerows(rows)
    return path


def write_integrated_summary() -> Path:
    rows = []
    for c in curves():
        for cent, _ in deep.CENTS:
            rec, rec_err, rec_den = deep.integrated_recovery(c.signal_path, cent, 15.0, 35.0)
            leak, leak_err, leak_den = integrated_tight_fraction_from_file(c.bkg_path, cent)
            rows.append(
                {
                    "model_key": c.key,
                    "model_label": c.label,
                    "centrality": cent,
                    "signal_recovery_15_35": rec,
                    "signal_recovery_err": rec_err,
                    "signal_recovery_den": rec_den,
                    "inclusive_jet_tight_fraction_15_35": leak,
                    "inclusive_jet_tight_fraction_err": leak_err,
                    "inclusive_jet_tight_fraction_den": leak_den,
                    "source_signal": str(c.signal_path.relative_to(REPO)),
                    "source_background": str(c.bkg_path.relative_to(REPO)),
                }
            )
    out = OUT_DIR / "nnstack_wp80_integrated_summary_15to35.csv"
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    return out


def main() -> int:
    missing = [str(path) for c in curves() for path in (c.signal_path, c.bkg_path) if not path.exists()]
    if missing:
        raise SystemExit("Missing required ROOT outputs:\n" + "\n".join(missing))
    recovery = draw_plot(
        "recovery",
        "nnstack_wp80_truth_photon_recovery_box_bdt_stack_1x3.png",
        "Truth-photon recovery",
        r"Truth photon $E_T$ [GeV]",
        (0.0, 0.64),
    )
    leakage = draw_plot(
        "inclusive_tight",
        "nnstack_wp80_inclusive_jet_tight_fraction_box_bdt_stack_1x3.png",
        "Inclusive-jet tight fraction",
        r"Photon candidate $E_T$ [GeV]",
        (0.0, 0.58),
    )
    summary = write_integrated_summary()
    print(f"[nnstackOverlay] recovery={recovery}")
    print(f"[nnstackOverlay] inclusive_tight={leakage}")
    print(f"[nnstackOverlay] summary={summary}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
