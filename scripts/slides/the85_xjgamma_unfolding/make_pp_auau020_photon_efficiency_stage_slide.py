#!/usr/bin/env python3
"""Make a pp vs AuAu 0-20 photon-efficiency stage slide for THE-85."""

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
from PIL import Image, ImageChops


REPO = Path(__file__).resolve().parents[3]
DEFAULT_OUTDIR = REPO / "dataOutput/the85_auau_xjgamma_unfolding_push/efficiency_stage"
DEFAULT_PP_ROOT = (
    REPO
    / "InputFiles/pp24/ppg12_globalmbd_mbddigi_componentmix_20260629/"
    / "jetMinPtScan_dphiScan_vz60_isoR40_isSliding_preselectionNewPPG12_tightNewPPG12_nonTightNewPPG12/"
    / "photonJet5and10and20componentmix_SIM/RecoilJets_photonjet5plus10plus20_componentmix_MERGED.root"
)
DEFAULT_AUAU_PNG = DEFAULT_OUTDIR / "auau020_photon_efficiency_stage_overlay.png"

PP_PANEL = DEFAULT_OUTDIR / "pp_current_recoiljets_photon_efficiency_stage_overlay.png"
PP_CSV = DEFAULT_OUTDIR / "pp_current_recoiljets_efficiency_stage_overlay.csv"
SLIDE_PNG = DEFAULT_OUTDIR / "pp_vs_auau020_photon_efficiency_stage_1x2_slide.png"
MANIFEST = DEFAULT_OUTDIR / "pp_vs_auau020_photon_efficiency_stage_1x2_manifest.json"

TITLE_COLOR = "#111827"
LABEL_COLOR = "#111827"
BODY_COLOR = "#334155"


ROOT.gROOT.SetBatch(True)


def require_hist(directory, name: str):
    hist = directory.Get(name) if directory else None
    if not hist:
        raise RuntimeError(f"Missing histogram: SIM/{name}")
    hist.SetDirectory(0)
    return hist


def optional_hist(directory, name: str):
    hist = directory.Get(name) if directory else None
    if not hist:
        return None
    hist.SetDirectory(0)
    return hist


def ratio_error(num: float, den: float, enum: float, eden: float) -> float:
    if den <= 0.0 or num < 0.0:
        return math.nan
    return math.sqrt((enum / den) ** 2 + ((num * eden) / (den * den)) ** 2)


def read_pp_recoiljets_points(root_path: Path) -> tuple[list[dict], dict]:
    handle = ROOT.TFile.Open(str(root_path), "READ")
    if not handle or handle.IsZombie():
        raise RuntimeError(f"Could not open ROOT file: {root_path}")
    sim = handle.Get("SIM")
    if not sim:
        raise RuntimeError(f"Missing SIM directory in {root_path}")

    h_den = require_hist(sim, "h_photonEffTruthDen_pTgamma_0")
    h_reco = require_hist(sim, "h_photonEffReco_pTgamma_0")
    h_tight = require_hist(sim, "h_photonEffRecoTight_pTgamma_0")
    h_iso = require_hist(sim, "h_photonEffRecoIso_pTgamma_0")
    h_tight_iso = require_hist(sim, "h_photonEffRecoTightIso_pTgamma_0")
    h_fig6_reco = optional_hist(sim, "h_photonEffPpg12Fig6Reco_pTgamma_0")
    h_fig6_iso = optional_hist(sim, "h_photonEffPpg12Fig6RecoIso_pTgamma_0")
    h_fig6_tight_iso = optional_hist(sim, "h_photonEffPpg12Fig6RecoTightIso_pTgamma_0")
    use_fig6_semantics = bool(h_fig6_reco and h_fig6_iso and h_fig6_tight_iso)

    rows: list[dict] = []
    for ib in range(1, h_den.GetNbinsX() + 1):
        lo = float(h_den.GetXaxis().GetBinLowEdge(ib))
        hi = float(h_den.GetXaxis().GetBinUpEdge(ib))
        mid = 0.5 * (lo + hi)
        if hi <= 15.0 or lo >= 35.0:
            continue
        den = float(h_den.GetBinContent(ib))
        eden = float(h_den.GetBinError(ib))
        if den <= 0.0:
            continue

        if use_fig6_semantics:
            jb_reco = h_fig6_reco.GetXaxis().FindBin(mid)
            jb_iso = h_fig6_iso.GetXaxis().FindBin(mid)
            jb_all = h_fig6_tight_iso.GetXaxis().FindBin(mid)
            reco_num = float(h_fig6_reco.GetBinContent(jb_reco))
            reco_err = float(h_fig6_reco.GetBinError(jb_reco))
            iso_num = float(h_fig6_iso.GetBinContent(jb_iso))
            iso_err = float(h_fig6_iso.GetBinError(jb_iso))
            all_num = float(h_fig6_tight_iso.GetBinContent(jb_all))
            all_err = float(h_fig6_tight_iso.GetBinError(jb_all))

            reco_eff = reco_num / den
            reco_eff_err = ratio_error(reco_num, den, reco_err, eden)
            id_cond_eff = all_num / iso_num if iso_num > 0.0 else math.nan
            id_cond_err = ratio_error(all_num, iso_num, all_err, iso_err) if iso_num > 0.0 else math.nan
            reco_id_eff = reco_eff * id_cond_eff if math.isfinite(id_cond_eff) else math.nan
            reco_id_err = (
                math.sqrt((id_cond_eff * reco_eff_err) ** 2 + (reco_eff * id_cond_err) ** 2)
                if math.isfinite(id_cond_eff) and math.isfinite(id_cond_err)
                else math.nan
            )

            stage_values = [
                (
                    "reco",
                    "Reco",
                    reco_num,
                    den,
                    reco_eff,
                    reco_eff_err,
                    "h_photonEffPpg12Fig6Reco_pTgamma_0 / h_photonEffTruthDen_pTgamma_0",
                ),
                (
                    "reco_id",
                    "Reco x ID",
                    all_num,
                    iso_num,
                    reco_id_eff,
                    reco_id_err,
                    "(Fig6 reco / truth) * (Fig6 tight+iso / Fig6 iso)",
                ),
                (
                    "reco_iso",
                    "Reco x iso",
                    iso_num,
                    den,
                    iso_num / den,
                    ratio_error(iso_num, den, iso_err, eden),
                    "h_photonEffPpg12Fig6RecoIso_pTgamma_0 / h_photonEffTruthDen_pTgamma_0",
                ),
                (
                    "reco_id_iso",
                    "Reco x ID x iso",
                    all_num,
                    den,
                    all_num / den,
                    ratio_error(all_num, den, all_err, eden),
                    "h_photonEffPpg12Fig6RecoTightIso_pTgamma_0 / h_photonEffTruthDen_pTgamma_0",
                ),
            ]
        else:
            stage_values = []
            for stage, label, hist, source in [
                ("reco", "Reco", h_reco, "h_photonEffReco_pTgamma_0 / h_photonEffTruthDen_pTgamma_0"),
                (
                    "reco_id",
                    "Reco x ID",
                    h_tight,
                    "h_photonEffRecoTight_pTgamma_0 / h_photonEffTruthDen_pTgamma_0",
                ),
                (
                    "reco_iso",
                    "Reco x iso",
                    h_iso,
                    "h_photonEffRecoIso_pTgamma_0 / h_photonEffTruthDen_pTgamma_0",
                ),
                (
                    "reco_id_iso",
                    "Reco x ID x iso",
                    h_tight_iso,
                    "h_photonEffRecoTightIso_pTgamma_0 / h_photonEffTruthDen_pTgamma_0",
                ),
            ]:
                jb = hist.GetXaxis().FindBin(mid)
                num = float(hist.GetBinContent(jb))
                enum = float(hist.GetBinError(jb))
                stage_values.append((stage, label, num, den, num / den, ratio_error(num, den, enum, eden), source))

        for stage, label, num, stage_den, eff, eff_err, source in stage_values:
            rows.append(
                {
                    "system": "pp",
                    "stage": stage,
                    "label": label,
                    "pt_low": lo,
                    "pt_high": hi,
                    "pt_mid": mid,
                    "numerator": num,
                    "denominator": stage_den,
                    "efficiency": eff,
                    "efficiency_error": eff_err,
                    "source": source,
                }
            )

    if not rows:
        raise RuntimeError(f"No pp RecoilJets efficiency points found in {root_path}")

    meta = {
        "root_path": str(root_path),
        "source": "our current local RecoilJets pp signal-MC output",
        "row": "jetMinPtScan_dphiScan_vz60_isoR40_isSliding_preselectionNewPPG12_tightNewPPG12_nonTightNewPPG12",
        "denominator": "SIM/h_photonEffTruthDen_pTgamma_0",
        "uses_ppg12_fig6_semantics": use_fig6_semantics,
        "reco_source": "SIM/h_photonEffPpg12Fig6Reco_pTgamma_0" if use_fig6_semantics else "SIM/h_photonEffReco_pTgamma_0",
        "reco_iso_source": "SIM/h_photonEffPpg12Fig6RecoIso_pTgamma_0" if use_fig6_semantics else "SIM/h_photonEffRecoIso_pTgamma_0",
        "reco_id_iso_source": "SIM/h_photonEffPpg12Fig6RecoTightIso_pTgamma_0" if use_fig6_semantics else "SIM/h_photonEffRecoTightIso_pTgamma_0",
        "reco_id_source": (
            "(SIM/h_photonEffPpg12Fig6Reco_pTgamma_0 / SIM/h_photonEffTruthDen_pTgamma_0) * "
            "(SIM/h_photonEffPpg12Fig6RecoTightIso_pTgamma_0 / SIM/h_photonEffPpg12Fig6RecoIso_pTgamma_0)"
            if use_fig6_semantics else "SIM/h_photonEffRecoTight_pTgamma_0"
        ),
        "isolation_contract": "pp baseline isoR40_isSliding in current RecoilJets output",
        "note": "Truth-binned same-definition efficiency-stage histograms from our RecoilJets pp signal-MC output; not digitized from the PPG12 paper and not read from Shuhang PPG12 efficiency ROOT files.",
    }
    handle.Close()
    return rows, meta


def write_rows_csv(rows: list[dict], out_csv: Path) -> None:
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with out_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def stage_xy(rows: list[dict], stage: str) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    ordered = sorted([row for row in rows if row["stage"] == stage], key=lambda row: float(row["pt_mid"]))
    return (
        np.array([float(row["pt_mid"]) for row in ordered], dtype=float),
        np.array([0.5 * (float(row["pt_high"]) - float(row["pt_low"])) for row in ordered], dtype=float),
        np.array([float(row["efficiency"]) for row in ordered], dtype=float),
        np.array([float(row["efficiency_error"]) for row in ordered], dtype=float),
    )


def render_pp_panel(rows: list[dict], meta: dict, out_png: Path) -> dict:
    colors = {
        "reco": "#000000",
        "reco_id": "#CC33A0",
        "reco_id_iso": "#269B2F",
    }
    labels = {
        "reco": r"$\epsilon_{\mathrm{reco}}$",
        "reco_id": r"$\epsilon_{\mathrm{reco}}\times\epsilon_{\mathrm{ID}}$",
        "reco_id_iso": r"$\epsilon_{\mathrm{reco}}\times\epsilon_{\mathrm{ID}}\times\epsilon_{\mathrm{iso}}$",
    }

    plt.rcParams.update(
        {
            "font.family": "Times New Roman",
            "font.serif": ["Times New Roman", "Times", "Nimbus Roman", "DejaVu Serif"],
            "mathtext.fontset": "stix",
            "axes.linewidth": 1.25,
        }
    )
    fig, ax = plt.subplots(figsize=(8.4, 7.4), dpi=220)
    fig.patch.set_facecolor("white")

    markers = {"reco": "s", "reco_id": "o", "reco_id_iso": "P"}
    for stage in ["reco", "reco_id", "reco_id_iso"]:
        px, pex, py, pey = stage_xy(rows, stage)
        ax.errorbar(
            px,
            py,
            xerr=pex,
            yerr=pey,
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
        r"$p{+}p\ \sqrt{s}=200\ \mathrm{GeV}$",
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

    meta = dict(meta)
    if meta.get("uses_ppg12_fig6_semantics"):
        curve_semantics = {
            "black": "SIM/h_photonEffPpg12Fig6Reco_pTgamma_0 / SIM/h_photonEffTruthDen_pTgamma_0",
            "magenta": "(Fig6 reco / truth) * (Fig6 tight+iso / Fig6 iso)",
            "green": "SIM/h_photonEffPpg12Fig6RecoTightIso_pTgamma_0 / SIM/h_photonEffTruthDen_pTgamma_0",
        }
    else:
        curve_semantics = {
            "black": "SIM/h_photonEffReco_pTgamma_0 / SIM/h_photonEffTruthDen_pTgamma_0",
            "magenta": "SIM/h_photonEffRecoTight_pTgamma_0 / SIM/h_photonEffTruthDen_pTgamma_0",
            "green": "SIM/h_photonEffRecoTightIso_pTgamma_0 / SIM/h_photonEffTruthDen_pTgamma_0",
        }
    meta.update(
        {
            "png": str(out_png),
            "pp_curve_semantics": curve_semantics,
            "pt_window": "15-35 GeV; plotted bins overlap the window",
        }
    )
    return meta


def crop_white_border(path: Path, pad: int = 16) -> Image.Image:
    image = Image.open(path).convert("RGBA")
    background = Image.new("RGBA", image.size, (255, 255, 255, 255))
    diff = ImageChops.difference(image, background).convert("L")
    bbox = diff.point(lambda px: 255 if px > 8 else 0).getbbox()
    if not bbox:
        return image
    left = max(0, bbox[0] - pad)
    upper = max(0, bbox[1] - pad)
    right = min(image.width, bbox[2] + pad)
    lower = min(image.height, bbox[3] + pad)
    return image.crop((left, upper, right, lower))


def paste_fit_bottom(canvas: Image.Image, image: Image.Image, box: tuple[int, int, int, int]) -> dict:
    x0, y0, x1, y1 = box
    max_w = x1 - x0
    max_h = y1 - y0
    scale = min(max_w / image.width, max_h / image.height)
    new_size = (int(round(image.width * scale)), int(round(image.height * scale)))
    resized = image.resize(new_size, Image.Resampling.LANCZOS)
    paste_x = x0 + (max_w - new_size[0]) // 2
    paste_y = y1 - new_size[1]
    canvas.alpha_composite(resized, (paste_x, paste_y))
    return {
        "box": list(box),
        "source_size": [image.width, image.height],
        "placed_size": list(new_size),
        "placed_xy": [paste_x, paste_y],
    }


def render_slide(
    pp_png: Path,
    auau_png: Path,
    out_png: Path,
    *,
    pp_panel_title: str,
    subtitle: str,
) -> dict:
    width, height = 2560, 1440
    canvas = Image.new("RGBA", (width, height), (255, 255, 255, 255))

    title_fig = plt.figure(figsize=(16, 9), dpi=160)
    title_fig.patch.set_facecolor("none")
    title_fig.text(
        0.055,
        0.955,
        r"Photon efficiency stages, $15<E_T^\gamma<35$ GeV",
        ha="left",
        va="top",
        fontsize=37,
        fontweight="bold",
        color=TITLE_COLOR,
    )
    title_fig.text(
        0.055,
        0.900,
        subtitle,
        ha="left",
        va="top",
        fontsize=20,
        color=BODY_COLOR,
    )
    title_fig.text(
        0.255,
        0.820,
        pp_panel_title,
        ha="center",
        va="center",
        fontsize=25,
        fontweight="bold",
        color=LABEL_COLOR,
    )
    title_fig.text(
        0.745,
        0.820,
        r"Au+Au embedded $\gamma$ MC, 0-20%",
        ha="center",
        va="center",
        fontsize=25,
        fontweight="bold",
        color=LABEL_COLOR,
    )
    title_fig.canvas.draw()
    title_buf = np.asarray(title_fig.canvas.buffer_rgba())
    title_img = Image.fromarray(title_buf)
    plt.close(title_fig)
    canvas.alpha_composite(title_img, (0, 0))

    pp_img = crop_white_border(pp_png)
    auau_img = crop_white_border(auau_png)
    left_box = (92, 255, 1244, 1360)
    right_box = (1316, 255, 2468, 1360)
    placement = {
        "pp": paste_fit_bottom(canvas, pp_img, left_box),
        "auau_0_20": paste_fit_bottom(canvas, auau_img, right_box),
    }

    out_png.parent.mkdir(parents=True, exist_ok=True)
    canvas.convert("RGB").save(out_png, quality=95)
    return {
        "png": str(out_png),
        "canvas": [width, height],
        "placement": placement,
        "layout": "1x2; plots bottom-aligned inside equal-width boxes",
        "pp_panel_title": pp_panel_title,
        "subtitle": subtitle,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pp-root", type=Path, default=DEFAULT_PP_ROOT)
    parser.add_argument("--auau-png", type=Path, default=DEFAULT_AUAU_PNG)
    parser.add_argument("--pp-png", type=Path, default=PP_PANEL)
    parser.add_argument("--pp-csv", type=Path, default=PP_CSV)
    parser.add_argument("--slide-png", type=Path, default=SLIDE_PNG)
    parser.add_argument("--manifest", type=Path, default=MANIFEST)
    parser.add_argument("--pp-panel-title", default=r"$p{+}p$ current RecoilJets output")
    parser.add_argument(
        "--subtitle",
        default=(
            "Truth-denominator signal-MC efficiencies from current RecoilJets outputs; "
            "both panels use sliding-isolation baseline rows."
        ),
    )
    args = parser.parse_args()

    if not args.pp_root.exists():
        raise FileNotFoundError(args.pp_root)
    if not args.auau_png.exists():
        raise FileNotFoundError(args.auau_png)

    rows, read_meta = read_pp_recoiljets_points(args.pp_root)
    write_rows_csv(rows, args.pp_csv)
    pp_meta = render_pp_panel(rows, read_meta, args.pp_png)
    pp_meta["csv"] = str(args.pp_csv)
    slide_meta = render_slide(
        args.pp_png,
        args.auau_png,
        args.slide_png,
        pp_panel_title=args.pp_panel_title,
        subtitle=args.subtitle,
    )
    manifest = {
        "pp": pp_meta,
        "auau_0_20": {
            "png": str(args.auau_png),
            "source_manifest": str(args.auau_png.with_name("auau020_photon_efficiency_stage_overlay_manifest.json")),
        },
        "slide": slide_meta,
    }
    args.manifest.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(args.pp_png)
    print(args.slide_png)
    print(args.manifest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
