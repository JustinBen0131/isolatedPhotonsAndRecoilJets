#!/usr/bin/env python3
"""PPG12 Fig. 6 efficiency overlay for the current photon+jet artifact.

Resolves the promoted photon+jet ROOT through the current-artifact pointer and
compares its dedicated PPG12 Fig. 6 efficiency-stage histograms with the direct
PPG12 TEfficiency readback.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import ROOT


REPO = Path("/Users/patsfan753/Desktop/ThesisAnalysis")
DEFAULT_CAMPAIGN = "the76_ppg12_fig24_photonjet_fix_20260702_014217"
DEFAULT_ROOT = (
    REPO
    / "dataOutput/ppg12Parity"
    / DEFAULT_CAMPAIGN
    / "final_roots/photonjet/RecoilJets_photonjet5plus10plus20_MERGED.root"
)
DEFAULT_CURRENT_POINTER = (
    REPO
    / "dataOutput/current_recoiljets_artifacts/current/pp_sim_photonjet_merged/current.json"
)
DEFAULT_PPG12_CSV = (
    REPO
    / "dataOutput/the85_auau_xjgamma_unfolding_push/efficiency_stage"
    / "ppg12_bdt_nom_efficiency_tefficiency_readback_full_remote_clean.csv"
)

HIST_DEN = "h_photonEffPpg12Fig6TruthDen_pTgamma_0"
HIST_RECO = "h_photonEffPpg12Fig6Reco_pTgamma_0"
HIST_RECO_ID = "h_photonEffPpg12Fig6RecoTight_pTgamma_0"
HIST_RECO_ISO = "h_photonEffPpg12Fig6RecoIso_pTgamma_0"
HIST_ALL = "h_photonEffPpg12Fig6RecoTightIso_pTgamma_0"

PLOT_XMIN = 10.0
PLOT_XMAX = 35.0

CURVES = [
    ("reco", r"$\varepsilon_{\mathrm{reco}}$", "#111111", "o"),
    ("reco_id", r"$\varepsilon_{\mathrm{reco}}\times\varepsilon_{\mathrm{ID}}$", "#d62d91", "o"),
    ("reco_id_iso", r"$\varepsilon_{\mathrm{reco}}\times\varepsilon_{\mathrm{ID}}\times\varepsilon_{\mathrm{iso}}$", "#2f9638", "o"),
]


def weighted_subset_efficiency(num: float, den: float, den_err: float) -> tuple[float, float]:
    """Weighted-binomial efficiency/error for a nested numerator and denominator."""
    if den <= 0.0 or num < 0.0 or num > den or den_err <= 0.0:
        return math.nan, math.nan
    eff = num / den
    n_eff = (den / den_err) ** 2
    err = math.sqrt(max(0.0, eff * (1.0 - eff) / n_eff)) if n_eff > 0.0 else math.nan
    return eff, err


def resolve_default_root(path_arg: Path | None) -> tuple[Path, str, str]:
    if path_arg is not None:
        return path_arg, "explicit --root", "explicit_root"
    if DEFAULT_CURRENT_POINTER.exists():
        payload = json.loads(DEFAULT_CURRENT_POINTER.read_text())
        roots = payload.get("root_paths") or []
        if roots:
            return Path(roots[0]), str(DEFAULT_CURRENT_POINTER), payload.get("campaign_tag", "current")
    return DEFAULT_ROOT, "fallback hardcoded fixed-rerun ROOT", DEFAULT_CAMPAIGN


def read_ppg12(csv_path: Path) -> dict[str, list[dict[str, float]]]:
    raw: dict[str, dict[float, dict[str, float]]] = {}
    with csv_path.open() as handle:
        for row in csv.DictReader(handle):
            lo = float(row["pt_low"])
            hi = float(row["pt_high"])
            if hi <= PLOT_XMIN or lo >= PLOT_XMAX:
                continue
            mid = float(row["pt_mid"])
            raw.setdefault(row["stage"], {})[mid] = {
                "pt_low": lo,
                "pt_high": hi,
                "pt_mid": mid,
                "eff": float(row["efficiency"]),
                "err_low": float(row["err_low"]),
                "err_high": float(row["err_high"]),
                "source_object": row.get("source_object", ""),
            }

    out = {"reco": [], "reco_id": [], "reco_id_iso": []}
    for mid in sorted(raw["reco"]):
        reco = raw["reco"][mid]
        ident = raw["id"][mid]
        product = reco["eff"] * ident["eff"]
        rel_low = math.hypot(reco["err_low"] / reco["eff"], ident["err_low"] / ident["eff"])
        rel_high = math.hypot(reco["err_high"] / reco["eff"], ident["err_high"] / ident["eff"])
        out["reco"].append(reco)
        out["reco_id"].append(
            {
                "pt_low": reco["pt_low"],
                "pt_high": reco["pt_high"],
                "pt_mid": mid,
                "eff": product,
                "err_low": product * rel_low,
                "err_high": product * rel_high,
                "source_object": "eff_reco_eta_0 * eff_id_eta_0",
            }
        )
        out["reco_id_iso"].append(raw["all"][mid])
    return out


def require_hist(sim_dir: ROOT.TDirectory, name: str) -> ROOT.TH1:
    hist = sim_dir.Get(name)
    if not hist:
        raise RuntimeError(f"Missing required histogram SIM/{name}")
    return hist


def hist_ratio_row(num: ROOT.TH1, den: ROOT.TH1, ppg12_row: dict[str, float]) -> dict[str, float]:
    mid = ppg12_row["pt_mid"]
    ib_num = num.GetXaxis().FindBin(mid)
    ib_den = den.GetXaxis().FindBin(mid)
    n = float(num.GetBinContent(ib_num))
    d = float(den.GetBinContent(ib_den))
    de = float(den.GetBinError(ib_den))
    eff, err = weighted_subset_efficiency(n, d, de)
    return {
        "pt_low": ppg12_row["pt_low"],
        "pt_high": ppg12_row["pt_high"],
        "pt_mid": mid,
        "eff": eff,
        "err_low": err,
        "err_high": err,
        "num": n,
        "den": d,
    }


def product_ratio_row(
    reco: ROOT.TH1,
    den: ROOT.TH1,
    tight_iso: ROOT.TH1,
    reco_iso: ROOT.TH1,
    ppg12_row: dict[str, float],
) -> dict[str, float]:
    """Mirror PPG12 Fig. 6 magenta: reco efficiency times conditional ID."""
    mid = ppg12_row["pt_mid"]
    ir = reco.GetXaxis().FindBin(mid)
    iden = den.GetXaxis().FindBin(mid)
    ia = tight_iso.GetXaxis().FindBin(mid)
    ii = reco_iso.GetXaxis().FindBin(mid)
    reco_num = float(reco.GetBinContent(ir))
    truth_den = float(den.GetBinContent(iden))
    all_num = float(tight_iso.GetBinContent(ia))
    iso_den = float(reco_iso.GetBinContent(ii))
    reco_eff, reco_err = weighted_subset_efficiency(reco_num, truth_den, float(den.GetBinError(iden)))
    id_eff, id_err = weighted_subset_efficiency(all_num, iso_den, float(reco_iso.GetBinError(ii)))
    eff = reco_eff * id_eff
    err = math.sqrt((id_eff * reco_err) ** 2 + (reco_eff * id_err) ** 2)
    return {
        "pt_low": ppg12_row["pt_low"],
        "pt_high": ppg12_row["pt_high"],
        "pt_mid": mid,
        "eff": eff,
        "err_low": err,
        "err_high": err,
        "num": all_num,
        "den": iso_den,
    }


def read_current(root_path: Path, ppg12: dict[str, list[dict[str, float]]]) -> dict[str, list[dict[str, float]]]:
    handle = ROOT.TFile.Open(str(root_path), "READ")
    if not handle or handle.IsZombie():
        raise RuntimeError(f"Could not open current ROOT: {root_path}")
    sim = handle.Get("SIM")
    if not sim:
        raise RuntimeError(f"Missing SIM directory in {root_path}")

    h_den = require_hist(sim, HIST_DEN)
    h_reco = require_hist(sim, HIST_RECO)
    require_hist(sim, HIST_RECO_ID)
    h_reco_iso = require_hist(sim, HIST_RECO_ISO)
    h_all = require_hist(sim, HIST_ALL)

    out = {
        "reco": [hist_ratio_row(h_reco, h_den, r) for r in ppg12["reco"]],
        "reco_id": [product_ratio_row(h_reco, h_den, h_all, h_reco_iso, r) for r in ppg12["reco_id"]],
        "reco_id_iso": [hist_ratio_row(h_all, h_den, r) for r in ppg12["reco_id_iso"]],
    }
    handle.Close()
    return out


def read_current_csv(csv_path: Path, ppg12: dict[str, list[dict[str, float]]]) -> dict[str, list[dict[str, float]]]:
    """Read a compact, provenance-pinned Fig. 6 efficiency extraction."""
    out: dict[str, list[dict[str, float]]] = {stage: [] for stage, *_rest in CURVES}
    with csv_path.open() as handle:
        for row in csv.DictReader(handle):
            stage = row["stage"]
            if stage not in out:
                raise RuntimeError(f"Unexpected current-efficiency stage {stage!r} in {csv_path}")
            out[stage].append(
                {
                    "pt_low": float(row["pt_low"]),
                    "pt_high": float(row["pt_high"]),
                    "pt_mid": float(row["pt_mid"]),
                    "eff": float(row["current_eff"]),
                    "err_low": float(row["current_err"]),
                    "err_high": float(row["current_err"]),
                    "num": float(row["current_num"]),
                    "den": float(row["current_den"]),
                }
            )

    for stage in out:
        out[stage].sort(key=lambda row: row["pt_mid"])
        expected = [row["pt_mid"] for row in ppg12[stage]]
        observed = [row["pt_mid"] for row in out[stage]]
        if observed != expected:
            raise RuntimeError(
                f"Current compact extraction does not match the PPG12 {stage} bin centers: "
                f"observed={observed}, expected={expected}"
            )
    return out


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_csv(path: Path, ppg12: dict[str, list[dict[str, float]]], current: dict[str, list[dict[str, float]]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "stage",
                "pt_low",
                "pt_high",
                "pt_mid",
                "ppg12_eff",
                "current_eff",
                "current_over_ppg12",
                "current_num",
                "current_den",
            ],
        )
        writer.writeheader()
        for stage, _label, _color, _marker in CURVES:
            by_mid = {r["pt_mid"]: r for r in current[stage]}
            for ref in ppg12[stage]:
                row = by_mid[ref["pt_mid"]]
                writer.writerow(
                    {
                        "stage": stage,
                        "pt_low": ref["pt_low"],
                        "pt_high": ref["pt_high"],
                        "pt_mid": ref["pt_mid"],
                        "ppg12_eff": ref["eff"],
                        "current_eff": row["eff"],
                        "current_over_ppg12": row["eff"] / ref["eff"] if ref["eff"] else "",
                        "current_num": row["num"],
                        "current_den": row["den"],
                    }
                )


def draw(
    path: Path,
    ppg12: dict[str, list[dict[str, float]]],
    current: dict[str, list[dict[str, float]]],
    *,
    pixel_width: int = 1544,
    pixel_height: int = 1996,
    dpi: int = 200,
) -> dict[str, dict[str, float]]:
    font_scale = max(0.55, min(1.0, pixel_width / 1544.0, pixel_height / 1996.0))

    def fs(value: float) -> float:
        return value * font_scale

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "mathtext.fontset": "dejavuserif",
            "axes.linewidth": 1.35,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.major.size": 7,
            "ytick.major.size": 7,
            "xtick.minor.size": 3.5,
            "ytick.minor.size": 3.5,
        }
    )
    fig, (ax, rax) = plt.subplots(
        2,
        1,
        figsize=(pixel_width / dpi, pixel_height / dpi),
        dpi=dpi,
        sharex=True,
        gridspec_kw={"height_ratios": [3.25, 1.0], "hspace": 0.04},
    )
    summary: dict[str, dict[str, float]] = {}

    for stage, label, color, marker in CURVES:
        ref = ppg12[stage]
        cur = current[stage]
        x = np.asarray([r["pt_mid"] for r in ref])
        xlo = np.asarray([r["pt_mid"] - r["pt_low"] for r in ref])
        xhi = np.asarray([r["pt_high"] - r["pt_mid"] for r in ref])
        y_ref = np.asarray([r["eff"] for r in ref])
        y_ref_err = np.asarray([[r["err_low"] for r in ref], [r["err_high"] for r in ref]])
        y_cur = np.asarray([r["eff"] for r in cur])
        y_cur_err = np.asarray([[r["err_low"] for r in cur], [r["err_high"] for r in cur]])
        ratio = y_cur / y_ref
        ratio_err = ratio * np.hypot(y_cur_err[1] / y_cur, y_ref_err[1] / y_ref)

        ax.errorbar(
            x,
            y_ref,
            xerr=[xlo, xhi],
            yerr=y_ref_err,
            fmt=marker,
            ms=6.0,
            mfc="white",
            mec=color,
            mew=1.45,
            ecolor=color,
            elinewidth=0.95,
            capsize=0,
            linestyle="none",
            zorder=5,
        )
        ax.errorbar(
            x,
            y_cur,
            xerr=[xlo, xhi],
            yerr=y_cur_err,
            fmt=marker,
            ms=4.6,
            mfc=color,
            mec=color,
            mew=0.9,
            ecolor=color,
            elinewidth=0.95,
            capsize=0,
            linestyle="none",
            zorder=6,
        )
        rax.errorbar(
            x,
            ratio,
            xerr=[xlo, xhi],
            yerr=ratio_err,
            fmt=marker,
            ms=4.4,
            mfc=color,
            mec=color,
            ecolor=color,
            elinewidth=0.9,
            capsize=0,
            linestyle="none",
        )
        summary[stage] = {
            "min_ratio": float(np.nanmin(ratio)),
            "max_ratio": float(np.nanmax(ratio)),
            "mean_ratio": float(np.nanmean(ratio)),
        }

    ax.set_xlim(PLOT_XMIN, PLOT_XMAX)
    ax.set_ylim(0.0, 1.15)
    ax.set_ylabel("Efficiency", fontsize=fs(19))
    ax.minorticks_on()
    ax.tick_params(which="both", top=True, right=True, labelsize=fs(13))
    ax.text(0.055, 0.935, "sPHENIX", transform=ax.transAxes, fontsize=fs(16.5), fontstyle="italic", fontweight="bold")
    ax.text(0.300, 0.935, "Internal", transform=ax.transAxes, fontsize=fs(16.5))
    ax.text(0.055, 0.865, r"$p$+$p$ $\sqrt{s}=200$ GeV", transform=ax.transAxes, fontsize=fs(14.0))
    ax.text(0.735, 0.935, "PYTHIA8", transform=ax.transAxes, fontsize=fs(16.5))
    ax.text(0.735, 0.865, r"$|\eta^\gamma|<0.7$", transform=ax.transAxes, fontsize=fs(14.0))
    shape_handles = [
        Line2D([0], [0], marker="o", color="0.15", mfc="white", mec="0.15", mew=1.45, lw=0, ms=6.0, label="PPG12 SDCC source"),
        Line2D([0], [0], marker="o", color="0.15", mfc="0.15", mec="0.15", lw=0, ms=5.0, label="Current output"),
    ]
    curve_handles = [
        Line2D([0], [0], marker="o", color=color, mfc=color, mec=color, lw=0, ms=5.4, label=label)
        for _stage, label, color, _marker in CURVES
    ]
    leg1 = ax.legend(handles=shape_handles, loc="lower left", bbox_to_anchor=(0.035, 0.235),
                     frameon=False, fontsize=fs(12.0), handlelength=1.0, labelspacing=0.5, borderpad=0.2)
    ax.add_artist(leg1)
    ax.legend(handles=curve_handles, loc="lower left", bbox_to_anchor=(0.035, 0.045),
              frameon=False, fontsize=fs(12.0), handlelength=1.0, labelspacing=0.55, borderpad=0.2)

    rax.axhline(1.0, color="0.45", lw=1.0, ls=(0, (4, 4)))
    rax.set_ylim(0.55, 1.45)
    rax.set_ylabel("Current / PPG12", fontsize=fs(13.5))
    rax.set_xlabel(r"$E_{\mathrm{T}}^{\gamma,\mathrm{truth}}$ [GeV]", fontsize=fs(17))
    rax.minorticks_on()
    rax.tick_params(which="both", top=True, right=True, labelsize=fs(12))
    fig.subplots_adjust(left=0.145, right=0.97, top=0.98, bottom=0.09)
    fig.savefig(path)
    plt.close(fig)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    source = parser.add_mutually_exclusive_group()
    source.add_argument("--root", type=Path, default=None)
    source.add_argument("--current-csv", type=Path, default=None)
    parser.add_argument("--ppg12-csv", type=Path, default=DEFAULT_PPG12_CSV)
    parser.add_argument("--out-dir", type=Path, default=None)
    parser.add_argument("--campaign-tag", default=None)
    parser.add_argument("--current-source-root", default=None)
    parser.add_argument("--current-source-sha256", default=None)
    parser.add_argument("--status-label", default="current production comparison")
    parser.add_argument("--production-artifact-exception", default=None)
    parser.add_argument("--pixel-width", type=int, default=1544)
    parser.add_argument("--pixel-height", type=int, default=1996)
    parser.add_argument("--dpi", type=int, default=200)
    args = parser.parse_args()

    current_csv_sha256 = None
    if args.current_csv is not None:
        if not args.current_source_root or not args.current_source_sha256:
            parser.error("--current-csv requires --current-source-root and --current-source-sha256")
        current_root = args.current_source_root
        root_resolution = "compact read-only histogram extraction from explicit remote ROOT"
        campaign_tag = args.campaign_tag or args.current_csv.parent.name
        current_csv_sha256 = sha256_file(args.current_csv)
    else:
        resolved_root, root_resolution, resolved_campaign = resolve_default_root(args.root)
        current_root = str(resolved_root)
        campaign_tag = args.campaign_tag or resolved_campaign
    if args.out_dir is None:
        args.out_dir = REPO / "dataOutput/ppg12Parity" / campaign_tag / "fig6_efficiency_current_direct"
    args.out_dir.mkdir(parents=True, exist_ok=True)
    out_png = args.out_dir / "ppg12_fig6_sdcc_vs_current_dedicated_overlay.png"
    out_csv = args.out_dir / "ppg12_fig6_sdcc_vs_current_dedicated_overlay_points.csv"
    out_manifest = args.out_dir / "ppg12_fig6_sdcc_vs_current_dedicated_overlay_manifest.json"

    ppg12 = read_ppg12(args.ppg12_csv)
    current = read_current_csv(args.current_csv, ppg12) if args.current_csv else read_current(Path(current_root), ppg12)
    write_csv(out_csv, ppg12, current)
    summary = draw(
        out_png,
        ppg12,
        current,
        pixel_width=args.pixel_width,
        pixel_height=args.pixel_height,
        dpi=args.dpi,
    )
    out_manifest.write_text(
        json.dumps(
            {
                "artifact": str(out_png),
                "comparison_csv": str(out_csv),
                "campaign_tag": campaign_tag,
                "status": args.status_label,
                "current_root": str(current_root),
                "current_root_resolution": root_resolution,
                "current_source_sha256": args.current_source_sha256,
                "current_compact_csv": str(args.current_csv) if args.current_csv else None,
                "current_compact_csv_sha256": current_csv_sha256,
                "production_artifact_exception": args.production_artifact_exception,
                "ppg12_csv": str(args.ppg12_csv),
                "ppg12_csv_sha256": sha256_file(args.ppg12_csv),
                "ppg12_source_root": "/sphenix/user/shuhangli/ppg12/efficiencytool/results/MC_efficiency_bdt_nom.root",
                "current_histograms": {
                    "denominator": f"SIM/{HIST_DEN}",
                    "reco": f"SIM/{HIST_RECO}",
                    "reco_id_direct_crosscheck": f"SIM/{HIST_RECO_ID}",
                    "reco_iso": f"SIM/{HIST_RECO_ISO}",
                    "reco_id": f"({HIST_RECO}/{HIST_DEN})*({HIST_ALL}/{HIST_RECO_ISO})",
                    "reco_id_iso": f"SIM/{HIST_ALL}",
                },
                "definition_note": "Exact PPG12 Fig. 6 construction: reco = Fig6Reco/Fig6TruthDen; reco*ID = (Fig6Reco/Fig6TruthDen)*(Fig6RecoTightIso/Fig6RecoIso); reco*ID*iso = Fig6RecoTightIso/Fig6TruthDen.",
                "current_statistical_errors": "Weighted-binomial subset errors using denominator effective entries; magenta product propagated as in the PPG12 macro.",
                "ratio_panel": "Current output / PPG12 SDCC source",
                "ratio_summary": summary,
                "pixel_dimensions": {
                    "width": args.pixel_width,
                    "height": args.pixel_height,
                    "dpi": args.dpi,
                },
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    print(out_png)
    print(out_csv)
    print(out_manifest)
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
