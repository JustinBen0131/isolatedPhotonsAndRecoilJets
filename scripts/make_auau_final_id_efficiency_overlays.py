#!/usr/bin/env python3
"""Final AuAu photon-ID efficiency overlays from completed SIM productions.

The plot inputs are the final merged local isSimEmbedded/isSimEmbeddedInclusive
outputs for box cuts, best BDT, clean MLP, and the BDT+MLP stack.  The script
writes the PNGs, CSV point tables, and a JSON manifest with the exact ROOT
files used so the comparison is reproducible.
"""

from __future__ import annotations

import csv
import json
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import ROOT


ROOT.gROOT.SetBatch(True)

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

import plot_target80_deep_dive as deep  # noqa: E402


OUT_DIR = (
    REPO
    / "dataOutput/auauFinalIDOverlays/"
    / "box_bdt_mlp_stack_wp80_20260515"
)

SIGNAL_REL = deep.SIGNAL_REL
BKG_REL = deep.BKG_REL
PT_MIN = 15.0
PT_MAX = 35.0


@dataclass(frozen=True)
class ModelOutput:
    key: str
    label: str
    root_dir: Path
    cfg: str
    color: str
    marker: str
    shift: float
    note: str

    @property
    def signal_path(self) -> Path:
        return self.root_dir / self.cfg / SIGNAL_REL

    @property
    def bkg_path(self) -> Path:
        return self.root_dir / self.cfg / BKG_REL


def model_outputs() -> list[ModelOutput]:
    target80 = REPO / "dataOutput/target80_all_available/bdt_target80_gated_20260512_001012"
    etfine = target80 / "analysis_config_etfine_15to35_target80"
    combined = REPO / "dataOutput/combinedSimOnlyEMBEDDED"
    mlp_root = combined / "mlp_wp80_release_20260512_152538"
    stack_root = combined / "bdt_mlp_nnstack_wp080_nnstack_wp080_prod_20260514_200206"
    return [
        ModelOutput(
            "box",
            "Box cuts",
            etfine / "simembedded",
            "preselectionNewPPG12_tightReference_nonTightReference_baseVariant",
            "#111111",
            "o",
            -0.24,
            "Reference box-cut photon ID from the BDT target80 production namespace.",
        ),
        ModelOutput(
            "bdt",
            r"Best BDT, WP80",
            etfine / "simembedded",
            "preselectionNewPPG12_tightAuAuEtFineCent7BDT_nonTightAuAuBDTComplement_baseVariant",
            "#0072B2",
            "s",
            -0.08,
            "Fine E_T x 7-centrality BDT, target-80 configuration.",
        ),
        ModelOutput(
            "mlp",
            r"MLP, WP80",
            mlp_root,
            "preselectionNewPPG12_tightAuAuCentInputBase3x3MLP_nonTightAuAuMLPComplement_baseVariant",
            "#009E73",
            "D",
            0.08,
            "Clean centrality-input base+3x3 MLP WP80 release.",
        ),
        ModelOutput(
            "stack",
            r"BDT+MLP stack, WP80",
            stack_root,
            "preselectionReference_tightAuAuBDTMLPStack_nonTightAuAuBDTMLPStackComplement_baseVariant",
            "#6F3FB5",
            "^",
            0.24,
            "NN-stack WP80 production; BDT and MLP scores are stack inputs.",
        ),
    ]


def background_model_outputs(models: list[ModelOutput]) -> list[ModelOutput]:
    out = []
    for model in models:
        root_dir = model.root_dir
        if root_dir.name == "simembedded":
            root_dir = root_dir.parent / "simembeddedinclusive"
        out.append(
            ModelOutput(
                model.key,
                model.label,
                root_dir,
                model.cfg,
                model.color,
                model.marker,
                model.shift,
                model.note,
            )
        )
    return out


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
    points = []
    pattern = re.compile(rf"^h_Eiso_tight_isoR30_pT_([0-9]+)_([0-9]+)_cent_{cent}$")
    keys = d.GetListOfKeys()
    bins = []
    for ik in range(keys.GetEntries()):
        name = keys.At(ik).GetName()
        match = pattern.match(name)
        if not match:
            continue
        lo = float(match.group(1))
        hi = float(match.group(2))
        mid = 0.5 * (lo + hi)
        if not (PT_MIN <= mid <= PT_MAX):
            continue
        bins.append((lo, hi, name))
    for lo, hi, tight_name in sorted(bins):
        h_tight = deep.get_hist(d, tight_name)
        h_nontight = deep.get_hist(d, tight_name.replace("h_Eiso_tight_", "h_Eiso_nonTight_", 1))
        tight, e_tight = sum_hist(h_tight)
        nontight, e_nontight = sum_hist(h_nontight)
        value, err, den = deep.ratio_point(tight, nontight, e_tight, e_nontight)
        if math.isfinite(value):
            points.append(
                {
                    "x": 0.5 * (lo + hi),
                    "y": value,
                    "ey": err,
                    "num": tight,
                    "den": den,
                }
            )
    f.Close()
    return points


def signal_truth_matched_tight_fraction_points_from_file(path: Path, cent: str) -> list[dict]:
    f, d = deep.open_sim(path)
    if not d:
        return []
    points = []
    pattern = re.compile(
        rf"^h_EisoReco_truthSigMatched_tight_isoR30_pT_([0-9]+)_([0-9]+)_cent_{cent}$"
    )
    keys = d.GetListOfKeys()
    bins = []
    for ik in range(keys.GetEntries()):
        name = keys.At(ik).GetName()
        match = pattern.match(name)
        if not match:
            continue
        lo = float(match.group(1))
        hi = float(match.group(2))
        mid = 0.5 * (lo + hi)
        if not (PT_MIN <= mid <= PT_MAX):
            continue
        bins.append((lo, hi, name))
    for lo, hi, tight_name in sorted(bins):
        h_tight = deep.get_hist(d, tight_name)
        h_all = deep.get_hist(d, tight_name.replace("h_EisoReco_truthSigMatched_tight_", "h_EisoReco_truthSigMatched_", 1))
        tight, e_tight = sum_hist(h_tight)
        all_reco, e_all = sum_hist(h_all)
        non_tight = max(all_reco - tight, 0.0)
        e_non_tight = math.hypot(e_all, e_tight)
        value, err, den = deep.ratio_point(tight, non_tight, e_tight, e_non_tight)
        if math.isfinite(value):
            points.append(
                {
                    "x": 0.5 * (lo + hi),
                    "y": value,
                    "ey": err,
                    "num": tight,
                    "den": den,
                }
            )
    f.Close()
    return points


def integrated_from_points(points: list[dict]) -> tuple[float, float, float]:
    num = sum(p["num"] for p in points)
    den = sum(p["den"] for p in points)
    if den <= 0:
        return math.nan, math.nan, 0.0
    value = num / den
    err = math.sqrt(max(value * (1.0 - value), 0.0) / den)
    return value, err, den


def style_axes(ax, ylabel: str | None, xlabel: str | None, ylim: tuple[float, float]):
    ax.tick_params(direction="in", which="both", top=True, right=True, length=6)
    ax.minorticks_on()
    ax.grid(True, color="#D0D0D0", alpha=0.45, linewidth=0.7)
    ax.set_xlim(14.5, 35.6)
    ax.set_ylim(*ylim)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=15)
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=15)
    ax.tick_params(labelsize=12)


def draw_figure_header(fig):
    fig.text(
        0.022,
        0.982,
        "sPHENIX",
        ha="left",
        va="top",
        fontsize=13,
        fontweight="bold",
        fontstyle="italic",
    )
    fig.text(0.084, 0.982, "Internal", ha="left", va="top", fontsize=13)
    fig.text(
        0.022,
        0.928,
        r"Au+Au embedded validation, $\sqrt{s_{NN}}=200$ GeV",
        ha="left",
        va="top",
        fontsize=10.7,
    )
    fig.text(
        0.022,
        0.885,
        r"$15 < E_T < 35$ GeV",
        ha="left",
        va="top",
        fontsize=10.7,
    )


def draw_1x3(
    *,
    models: list[ModelOutput],
    out_name: str,
    csv_name: str,
    ylabel: str,
    xlabel: str,
    ylim: tuple[float, float],
    point_reader: Callable[[ModelOutput, str], list[dict]],
):
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(17.2, 5.2), sharey=True)
    rows = []
    for ic, (cent, title) in enumerate(deep.CENTS):
        ax = axes[ic]
        for model in models:
            pts = point_reader(model, cent)
            if not pts:
                continue
            x = np.array([p["x"] + model.shift for p in pts], dtype=float)
            y = np.array([p["y"] for p in pts], dtype=float)
            ey = np.array([p["ey"] for p in pts], dtype=float)
            ax.errorbar(
                x,
                y,
                yerr=ey,
                fmt=model.marker,
                linestyle="none",
                color=model.color,
                markerfacecolor=model.color,
                markeredgecolor=model.color,
                markersize=6.8,
                elinewidth=1.3,
                capsize=2.4,
                label=model.label,
                zorder=3,
            )
            for p in pts:
                rows.append(
                    {
                        "plot": Path(out_name).stem,
                        "model_key": model.key,
                        "model_label": model.label,
                        "centrality": cent,
                        "pt_mid": p["x"],
                        "value": p["y"],
                        "error": p["ey"],
                        "numerator": p["num"],
                        "denominator": p["den"],
                        "signal_root": str(model.signal_path.relative_to(REPO)) if model.signal_path.exists() else "",
                        "background_root": str(model.bkg_path.relative_to(REPO)) if model.bkg_path.exists() else "",
                        "note": model.note,
                    }
                )
        ax.set_title(title, fontsize=17, pad=12)
        style_axes(ax, ylabel if ic == 0 else None, xlabel if ic == 1 else None, ylim)
        if ic == 2:
            leg = ax.legend(
                loc="upper right",
                frameon=True,
                framealpha=0.94,
                fontsize=10.4,
                borderpad=0.55,
                handletextpad=0.45,
            )
            leg.get_frame().set_linewidth(0.0)
    draw_figure_header(fig)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.88), w_pad=1.5)
    png = OUT_DIR / out_name
    fig.savefig(png, dpi=190)
    plt.close(fig)

    csv_path = OUT_DIR / csv_name
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()) if rows else [])
        if rows:
            writer.writeheader()
            writer.writerows(rows)
    return png, csv_path


def write_manifest(models: list[ModelOutput], bkg_models: list[ModelOutput]) -> Path:
    manifest = {
        "schema": "AUAU_FINAL_ID_OVERLAY_V1",
        "pt_range": [PT_MIN, PT_MAX],
        "centrality_bins": [{"suffix": c, "label": label} for c, label in deep.CENTS],
        "definitions": {
            "truth_photon_recovery": "1 - h_unfoldTruthPhoMisses_pTgamma_isoR30_isSliding / h_unfoldTruthPho_pTgamma from isSimEmbedded",
            "signal_candidate_tight_fraction": "h_EisoReco_truthSigMatched_tight / h_EisoReco_truthSigMatched from isSimEmbedded",
            "inclusive_jet_tight_fraction": "h_Eiso_tight / (h_Eiso_tight + h_Eiso_nonTight) from isSimEmbeddedInclusive",
        },
        "models": [
            {
                "key": model.key,
                "label": model.label,
                "signal_root": str(model.signal_path.relative_to(REPO)),
                "background_root": str(bkg.bkg_path.relative_to(REPO)),
                "cfg": model.cfg,
                "note": model.note,
            }
            for model, bkg in zip(models, bkg_models)
        ],
    }
    out = OUT_DIR / "final_id_overlay_manifest.json"
    out.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return out


def write_integrated_summary(
    models: list[ModelOutput],
    bkg_models: list[ModelOutput],
) -> Path:
    rows = []
    for signal_model, bkg_model in zip(models, bkg_models):
        for cent, _ in deep.CENTS:
            rec_pts = deep.recovery_points(signal_model.signal_path, cent, PT_MIN, PT_MAX)
            sig_tight_pts = signal_truth_matched_tight_fraction_points_from_file(signal_model.signal_path, cent)
            bkg_tight_pts = tight_fraction_points_from_file(bkg_model.bkg_path, cent)
            rec, rec_err, rec_den = integrated_from_points(rec_pts)
            sig, sig_err, sig_den = integrated_from_points(sig_tight_pts)
            bkg, bkg_err, bkg_den = integrated_from_points(bkg_tight_pts)
            rows.append(
                {
                    "model_key": signal_model.key,
                    "model_label": signal_model.label,
                    "centrality": cent,
                    "truth_photon_recovery_15_35": rec,
                    "truth_photon_recovery_err": rec_err,
                    "truth_photon_recovery_den": rec_den,
                    "signal_candidate_tight_fraction_15_35": sig,
                    "signal_candidate_tight_fraction_err": sig_err,
                    "signal_candidate_tight_fraction_den": sig_den,
                    "inclusive_jet_tight_fraction_15_35": bkg,
                    "inclusive_jet_tight_fraction_err": bkg_err,
                    "inclusive_jet_tight_fraction_den": bkg_den,
                    "signal_root": str(signal_model.signal_path.relative_to(REPO)),
                    "background_root": str(bkg_model.bkg_path.relative_to(REPO)),
                }
            )
    out = OUT_DIR / "final_id_overlay_integrated_summary_15to35.csv"
    with out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    return out


def main() -> int:
    models = model_outputs()
    bkg_models = background_model_outputs(models)
    missing = []
    for model, bkg in zip(models, bkg_models):
        for path in (model.signal_path, bkg.bkg_path):
            if not path.exists():
                missing.append(str(path))
    if missing:
        raise SystemExit("Missing required ROOT files:\n" + "\n".join(missing))

    recovery_png, recovery_csv = draw_1x3(
        models=models,
        out_name="truth_photon_recovery_box_bdt_mlp_stack_1x3.png",
        csv_name="truth_photon_recovery_box_bdt_mlp_stack_1x3.csv",
        ylabel="Truth-photon recovery",
        xlabel=r"Truth photon $E_T$ [GeV]",
        ylim=(0.0, 0.70),
        point_reader=lambda model, cent: deep.recovery_points(model.signal_path, cent, PT_MIN, PT_MAX),
    )
    signal_tight_png, signal_tight_csv = draw_1x3(
        models=models,
        out_name="signal_candidate_tight_fraction_box_bdt_mlp_stack_1x3.png",
        csv_name="signal_candidate_tight_fraction_box_bdt_mlp_stack_1x3.csv",
        ylabel="Signal candidate tight fraction",
        xlabel=r"Photon candidate $E_T$ [GeV]",
        ylim=(0.0, 1.05),
        point_reader=lambda model, cent: signal_truth_matched_tight_fraction_points_from_file(model.signal_path, cent),
    )
    bkg_tight_png, bkg_tight_csv = draw_1x3(
        models=bkg_models,
        out_name="inclusive_jet_tight_fraction_box_bdt_mlp_stack_1x3.png",
        csv_name="inclusive_jet_tight_fraction_box_bdt_mlp_stack_1x3.csv",
        ylabel="Inclusive-jet tight fraction",
        xlabel=r"Photon candidate $E_T$ [GeV]",
        ylim=(0.0, 0.62),
        point_reader=lambda model, cent: tight_fraction_points_from_file(model.bkg_path, cent),
    )
    summary = write_integrated_summary(models, bkg_models)
    manifest = write_manifest(models, bkg_models)
    for path in [
        recovery_png,
        signal_tight_png,
        bkg_tight_png,
        recovery_csv,
        signal_tight_csv,
        bkg_tight_csv,
        summary,
        manifest,
    ]:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
