#!/usr/bin/env python3
"""Make the THE-85 AuAu 0-20% photon efficiency-stage overlay."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import ROOT


ROOT.gROOT.SetBatch(True)

REPO = Path(__file__).resolve().parents[3]
DEFAULT_ROOT = (
    REPO
    / "InputFiles/the85_signal_eff_stage_20260628/"
    / "RecoilJets_embeddedPhoton12plus20_MERGED.root"
)
DEFAULT_OUTDIR = (
    REPO
    / "dataOutput/the85_auau_xjgamma_unfolding_push/efficiency_stage"
)


def require_hist(directory, name: str):
    hist = directory.Get(name)
    if not hist:
        raise RuntimeError(f"Missing histogram: SIM/{name}")
    hist.SetDirectory(0)
    return hist


def ratio_error(num: float, den: float, enum: float, eden: float) -> float:
    if den <= 0.0 or num < 0.0:
        return math.nan
    # Conservative weighted-hist ratio propagation. The categories are not
    # strictly independent for all curves, so this is used as a visual QA error.
    return math.sqrt((enum / den) ** 2 + ((num * eden) / (den * den)) ** 2)


def read_points(root_path: Path) -> tuple[list[dict], dict]:
    handle = ROOT.TFile.Open(str(root_path), "READ")
    if not handle or handle.IsZombie():
        raise RuntimeError(f"Could not open ROOT file: {root_path}")
    sim = handle.Get("SIM")
    if not sim:
        raise RuntimeError(f"Missing SIM directory in {root_path}")

    h_den = require_hist(sim, "h_photonEffTruthDen_pTgamma_cent_0_20")
    h_reco = require_hist(sim, "h_photonEffReco_pTgamma_cent_0_20")
    h_reco_iso = require_hist(
        sim, "h_photonEffRecoIso_pTgamma_isoR40_isSliding_cent_0_20"
    )
    h_tight_iso = require_hist(
        sim, "h_photonEffRecoTightIso_pTgamma_isoR40_isSliding_cent_0_20"
    )

    fixed_iso_present = bool(sim.Get("h_photonEffRecoIso_pTgamma_isoR30_fixedIso4GeV_cent_0_20"))
    sliding_eff_stage_present = bool(h_reco_iso and h_tight_iso)

    rows: list[dict] = []
    for ib in range(1, h_den.GetNbinsX() + 1):
        lo = float(h_den.GetXaxis().GetBinLowEdge(ib))
        hi = float(h_den.GetXaxis().GetBinUpEdge(ib))
        mid = 0.5 * (lo + hi)
        if lo < 15.0 or hi > 35.0:
            continue

        den = float(h_den.GetBinContent(ib))
        eden = float(h_den.GetBinError(ib))
        if den <= 0.0:
            continue

        def content(hist):
            jb = hist.GetXaxis().FindBin(mid)
            return float(hist.GetBinContent(jb)), float(hist.GetBinError(jb))

        reco, ereco = content(h_reco)
        reco_iso, ereco_iso = content(h_reco_iso)
        tight_iso, etight_iso = content(h_tight_iso)

        reco_eff = reco / den
        ereco_eff = ratio_error(reco, den, ereco, eden)
        id_eff = tight_iso / reco_iso if reco_iso > 0.0 else math.nan
        eid_eff = (
            ratio_error(tight_iso, reco_iso, etight_iso, ereco_iso)
            if reco_iso > 0.0
            else math.nan
        )
        reco_id_eff = reco_eff * id_eff if math.isfinite(id_eff) else math.nan
        ereco_id_eff = (
            math.sqrt((id_eff * ereco_eff) ** 2 + (reco_eff * eid_eff) ** 2)
            if math.isfinite(id_eff) and math.isfinite(eid_eff)
            else math.nan
        )
        reco_id_iso_eff = tight_iso / den
        ereco_id_iso_eff = ratio_error(tight_iso, den, etight_iso, eden)

        for stage, label, value, err, num, den_for_stage, source in [
            (
                "reco",
                "Reco",
                reco_eff,
                ereco_eff,
                reco,
                den,
                "h_photonEffReco_pTgamma_cent_0_20 / h_photonEffTruthDen_pTgamma_cent_0_20",
            ),
            (
                "reco_id",
                "Reco x ID",
                reco_id_eff,
                ereco_id_eff,
                tight_iso,
                reco_iso,
                "(h_photonEffReco_pTgamma / h_photonEffTruthDen_pTgamma) * "
                "(h_photonEffRecoTightIso_pTgamma_isoR40_isSliding / "
                "h_photonEffRecoIso_pTgamma_isoR40_isSliding)",
            ),
            (
                "reco_id_iso",
                "Reco x ID x iso",
                reco_id_iso_eff,
                ereco_id_iso_eff,
                tight_iso,
                den,
                "h_photonEffRecoTightIso_pTgamma_isoR40_isSliding_cent_0_20 / "
                "h_photonEffTruthDen_pTgamma_cent_0_20",
            ),
        ]:
            rows.append(
                {
                    "centrality": "0-20",
                    "stage": stage,
                    "label": label,
                    "pt_low": lo,
                    "pt_high": hi,
                    "pt_mid": mid,
                    "numerator": num,
                    "denominator": den_for_stage,
                    "efficiency": value,
                    "efficiency_error": err,
                    "source": source,
                }
            )

    if not rows:
        raise RuntimeError(
            "No 15-35 GeV points were produced. Check the ROOT photon-efficiency binning."
        )

    meta = {
        "root_path": str(root_path),
        "denominator": "SIM/h_photonEffTruthDen_pTgamma_cent_0_20",
        "reco_source": "SIM/h_photonEffReco_pTgamma_cent_0_20",
        "reco_iso_source": "SIM/h_photonEffRecoIso_pTgamma_isoR40_isSliding_cent_0_20",
        "reco_tight_iso_source": "SIM/h_photonEffRecoTightIso_pTgamma_isoR40_isSliding_cent_0_20",
        "reco_id_iso_source": (
            "SIM/h_photonEffRecoTightIso_pTgamma_isoR40_isSliding_cent_0_20"
        ),
        "reco_id_semantics": (
            "PPG12-style product: reco efficiency times tight-ID efficiency "
            "conditional on matched reco photon passing the same sliding-isolation view."
        ),
        "isolation_contract": "isoR40_isSliding; Eiso < 7.57 - 0.0658*c; sideGap=0",
        "fixed_iso_eff_stage_present": fixed_iso_present,
        "sliding_eff_stage_present": sliding_eff_stage_present,
        "note": (
            "No unfolding-miss fallback is used. The final iso curve is the direct "
            "sliding-isolation tight numerator divided by the truth-isolated denominator."
        ),
    }
    handle.Close()
    return rows, meta


def write_csv(rows: list[dict], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def stage_arrays(rows: list[dict], stage: str):
    sub = [r for r in rows if r["stage"] == stage]
    return (
        np.array([r["pt_mid"] for r in sub], dtype=float),
        np.array([(r["pt_high"] - r["pt_low"]) * 0.5 for r in sub], dtype=float),
        np.array([r["efficiency"] for r in sub], dtype=float),
        np.array([r["efficiency_error"] for r in sub], dtype=float),
    )


def render(rows: list[dict], meta: dict, out_png: Path) -> None:
    colors = {
        "reco": "#000000",
        "reco_id": "#CC33A0",
        "reco_id_iso": "#269B2F",
    }
    markers = {"reco": "s", "reco_id": "o", "reco_id_iso": "P"}
    labels = {
        "reco": r"$\epsilon_{\mathrm{reco}}$",
        "reco_id": r"$\epsilon_{\mathrm{reco}}\times\epsilon_{\mathrm{ID}}$",
        "reco_id_iso": r"$\epsilon_{\mathrm{reco}}\times\epsilon_{\mathrm{ID}}\times\epsilon_{\mathrm{iso}}$",
    }

    fig, ax = plt.subplots(figsize=(8.4, 7.4), dpi=220)
    for stage in ["reco", "reco_id", "reco_id_iso"]:
        x, ex, y, ey = stage_arrays(rows, stage)
        ax.errorbar(
            x,
            y,
            xerr=ex,
            yerr=ey,
            fmt=markers[stage],
            color=colors[stage],
            markerfacecolor=colors[stage],
            markeredgecolor=colors[stage],
            markeredgewidth=1.1,
            markersize=8.5,
            elinewidth=1.6,
            capsize=0.0,
            linestyle="none",
            label=labels[stage],
            zorder=4,
        )

    ax.set_xlim(15.0, 35.0)
    ax.set_ylim(0.0, 1.20)
    ax.set_xlabel(r"$E_T^{\gamma,\mathrm{truth}}\ \mathrm{[GeV]}$", fontsize=22, ha="right", x=1.0)
    ax.set_ylabel("Efficiency", fontsize=22)
    ax.tick_params(direction="in", which="both", top=True, right=True, labelsize=18, length=8, width=1.2)
    ax.tick_params(which="minor", length=4, width=1.0)
    ax.minorticks_on()
    for spine in ax.spines.values():
        spine.set_linewidth(1.25)
    ax.legend(
        loc="lower left",
        bbox_to_anchor=(0.075, 0.105),
        frameon=False,
        facecolor="white",
        fontsize=20,
        handlelength=1.4,
        borderaxespad=0.0,
    )
    ax.text(
        0.05,
        0.965,
        r"$\it{\bf{sPHENIX}}$ Internal",
        transform=ax.transAxes,
        fontsize=22,
        va="top",
        ha="left",
    )
    ax.text(
        0.05,
        0.900,
        r"Au+Au embedded $\gamma$ MC, 0-20%",
        transform=ax.transAxes,
        fontsize=21,
        va="top",
        ha="left",
        color="#000000",
    )
    ax.text(
        0.95,
        0.965,
        "PYTHIA8",
        transform=ax.transAxes,
        fontsize=22,
        va="top",
        ha="right",
        color="#000000",
    )
    ax.text(
        0.95,
        0.900,
        r"$|\eta^\gamma| < 0.7$",
        transform=ax.transAxes,
        fontsize=21,
        va="top",
        ha="right",
        color="#000000",
    )

    fig.tight_layout(pad=1.0)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--stem", default="auau020_photon_efficiency_stage_overlay")
    args = parser.parse_args()

    rows, meta = read_points(args.root)
    args.outdir.mkdir(parents=True, exist_ok=True)
    png = args.outdir / f"{args.stem}.png"
    csv_path = args.outdir / f"{args.stem}.csv"
    manifest = args.outdir / f"{args.stem}_manifest.json"

    write_csv(rows, csv_path)
    render(rows, meta, png)
    meta.update(
        {
            "png": str(png),
            "csv": str(csv_path),
            "stages": ["reco", "reco_id", "reco_id_iso"],
            "pt_window": "15-35 GeV; bins are accepted only if low edge >=15 and high edge <=35",
        }
    )
    manifest.write_text(json.dumps(meta, indent=2, sort_keys=True) + "\n")
    print(png)
    print(csv_path)
    print(manifest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
