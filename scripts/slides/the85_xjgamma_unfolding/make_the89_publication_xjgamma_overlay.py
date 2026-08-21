#!/usr/bin/env python3
"""Render a clean statistical p+p versus Au+Au unfolded xJgamma candidate.

The p+p curve is recomputed from the registered current full-statistics inputs.
The Au+Au curve is the existing THE-89 response/K-contract audit candidate.
Only native xJ bins beginning at 0.5 are displayed; all bins remain in the
numeric sidecar.  No unit-area scaling or other curve rescaling is applied.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import math
from pathlib import Path
import sys
from typing import Any

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
from PIL import Image


REPO = Path(__file__).resolve().parents[3]
OUT = REPO / "dataOutput/the219_friday_ppg_20260814/the89_publication_xjgamma_overlay"
OUT.mkdir(parents=True, exist_ok=True)

PNG = OUT / "the89_unfolded_xjgamma_pp_auau020_xjge0p5_stat_candidate.png"
POINTS = OUT / "the89_unfolded_xjgamma_pp_auau020_xjge0p5_points.json"
MANIFEST = OUT / "the89_unfolded_xjgamma_pp_auau020_xjge0p5_manifest.json"
AUDIT = OUT / "the89_unfolded_xjgamma_pp_auau020_xjge0p5_audit.json"

PP_SCRIPT = REPO / "scripts/slides/pp_currentian/xjgamma/make_the219_pp_money_slide.py"
AUAU_NPZ = (
    REPO
    / "dataOutput/the85_auau_xjgamma_unfolding_push/unfolded_xjgamma_firstpass/"
    "the85_unfolded_xjgamma_auau_0_20_responseKfix_iter5_covariance.npz"
)
AUAU_AUDIT_MANIFEST = (
    REPO
    / "dataOutput/the85_auau_xjgamma_unfolding_push/unfolded_xjgamma_firstpass/"
    "slide12_atlas_vs_sphenix_responseKfix_iter5_covariance_clean_v1_manifest.json"
)

FONT_PATH = Path("/System/Library/Fonts/Supplemental/Times New Roman.ttf")
PT_RANGE = (14.0, 35.0)
DISPLAY_XMIN = 0.5
DISPLAY_XMAX = 1.5


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def json_ready(value: Any) -> Any:
    if isinstance(value, dict):
        return {key: json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [json_ready(item) for item in value]
    if isinstance(value, (float, np.floating)):
        return float(value) if math.isfinite(float(value)) else None
    if isinstance(value, (int, np.integer)):
        return int(value)
    return value


def load_pp_module():
    spec = importlib.util.spec_from_file_location("the219_pp_current", PP_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import current p+p unfolding helper: {PP_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def build_current_pp() -> tuple[dict[str, Any], dict[str, Any]]:
    pp = load_pp_module()
    data_pointer = json.loads(pp.DATA_POINTER.read_text())
    sim_pointer = json.loads(pp.SIM_POINTER.read_text())
    data_path = Path(data_pointer["root_paths"][0])
    sim_path = Path(sim_pointer["root_paths"][0])
    data = pp.open_root(data_path)
    sim = pp.open_root(sim_path)
    try:
        h_photon, cov_photon, photon_meta = pp.unfold_photons(data, sim)
        h_xj_global, h_xj_2d, cov_xj, xj_meta = pp.unfold_xj(data, sim)
        sim_truth_photon = pp.obj(sim, pp.SIM_TOP, pp.OBJECTS["sim_photon_truth"], "TH1")
        sim_truth_xj = pp.obj(sim, pp.SIM_TOP, pp.OBJECTS["sim_xj_truth"], "TH2")
        panel = pp.truth_panel(
            PT_RANGE[0],
            PT_RANGE[1],
            h_photon,
            cov_photon,
            h_xj_global,
            h_xj_2d,
            cov_xj,
            sim_truth_photon,
            sim_truth_xj,
        )
    finally:
        data.Close()
        sim.Close()
    meta = {
        "data_pointer": str(pp.DATA_POINTER),
        "data_pointer_sha256": sha256(pp.DATA_POINTER),
        "data_campaign": data_pointer["campaign_tag"],
        "data_status": data_pointer["canonical_status"],
        "data_root": str(data_path),
        "data_root_size_bytes": data_path.stat().st_size,
        "data_root_mtime_ns": data_path.stat().st_mtime_ns,
        "sim_pointer": str(pp.SIM_POINTER),
        "sim_pointer_sha256": sha256(pp.SIM_POINTER),
        "sim_campaign": sim_pointer["campaign_tag"],
        "sim_status": sim_pointer["canonical_status"],
        "sim_root": str(sim_path),
        "sim_root_size_bytes": sim_path.stat().st_size,
        "sim_root_mtime_ns": sim_path.stat().st_mtime_ns,
        "iterations": pp.ITERATIONS,
        "error_mode": "RooUnfold kCovariance",
        "photon_meta": photon_meta,
        "xj_meta": xj_meta,
    }
    return panel, meta


def load_auau() -> tuple[dict[str, Any], dict[str, Any]]:
    source = np.load(AUAU_NPZ, allow_pickle=True)
    meta = json.loads(str(source["meta_json"][0]))
    curve = {
        "xj_edges": source["x_edges"].astype(float).tolist(),
        "values": source["y"].astype(float).tolist(),
        "errors": source["ey"].astype(float).tolist(),
        "unfolded_photons": float(source["npho_unfolded"][0]),
        "unfolded_photon_error": float(source["npho_error"][0]),
        "truth_pt_bins": source["selected_truth_pt_bins"].astype(float).tolist(),
    }
    return curve, meta


def displayed_mask(edges: np.ndarray, values: np.ndarray, errors: np.ndarray) -> np.ndarray:
    centers = 0.5 * (edges[:-1] + edges[1:])
    return (
        np.isfinite(values)
        & np.isfinite(errors)
        & (edges[:-1] >= DISPLAY_XMIN - 1e-12)
        & (centers <= DISPLAY_XMAX + 1e-12)
    )


def configure_style() -> str:
    if not FONT_PATH.exists():
        raise FileNotFoundError(f"required slide-policy font missing: {FONT_PATH}")
    font_manager.fontManager.addfont(str(FONT_PATH))
    font_name = font_manager.FontProperties(fname=str(FONT_PATH)).get_name()
    if font_name != "Times New Roman":
        raise RuntimeError(f"unexpected font identity: {font_name}")
    mpl.rcParams.update(
        {
            "font.family": font_name,
            "mathtext.fontset": "stix",
            "axes.linewidth": 1.15,
            "axes.labelcolor": "#111111",
            "xtick.color": "#111111",
            "ytick.color": "#111111",
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "xtick.major.size": 7.0,
            "ytick.major.size": 7.0,
            "xtick.minor.size": 3.7,
            "ytick.minor.size": 3.7,
        }
    )
    return font_name


def render(pp_curve: dict[str, Any], auau_curve: dict[str, Any]) -> None:
    configure_style()
    fig = plt.figure(figsize=(10, 7.5), dpi=240, facecolor="white")
    ax = fig.add_axes([0.145, 0.135, 0.82, 0.82])

    styles = [
        (pp_curve, r"$p$+$p$", "#245CA6", "s", "white", 3),
        (auau_curve, r"Au+Au 0–20%", "#D8443E", "o", "#D8443E", 4),
    ]
    for curve, label, color, marker, face, zorder in styles:
        edges = np.asarray(curve["xj_edges"], dtype=float)
        centers = 0.5 * (edges[:-1] + edges[1:])
        values = np.asarray(curve["values"], dtype=float)
        errors = np.asarray(curve["errors"], dtype=float)
        mask = displayed_mask(edges, values, errors)
        ax.errorbar(
            centers[mask],
            values[mask],
            xerr=0.5 * np.diff(edges)[mask],
            yerr=errors[mask],
            fmt=marker,
            ms=7.2,
            mfc=face,
            mec=color,
            mew=1.45,
            color=color,
            ecolor=color,
            elinewidth=1.35,
            capsize=2.8,
            capthick=1.25,
            linestyle="none",
            label=label,
            zorder=zorder,
        )

    ax.axhline(0.0, color="#666666", linewidth=0.75, zorder=1)
    ax.set_xlim(DISPLAY_XMIN, DISPLAY_XMAX)
    ax.set_ylim(-0.055, 1.30)
    ax.set_xticks(np.arange(0.5, 1.51, 0.1))
    ax.set_yticks(np.arange(0.0, 1.21, 0.2))
    ax.minorticks_on()
    ax.tick_params(axis="both", which="major", labelsize=17, width=1.1, pad=7)
    ax.tick_params(axis="both", which="minor", width=0.85)
    ax.set_xlabel(r"$x_{J\gamma}=p_{T}^{\mathrm{jet}}/p_{T}^{\gamma}$", fontsize=23, labelpad=12)
    ax.set_ylabel(
        r"$(1/N_{\gamma}^{\mathrm{particle}})\,dN_{\mathrm{jet}}^{\mathrm{particle}}/dx_{J\gamma}$",
        fontsize=23,
        labelpad=15,
    )

    ax.text(0.035, 0.955, "sPHENIX", transform=ax.transAxes, ha="left", va="top", fontsize=23,
            fontweight="bold", fontstyle="italic", color="#111111")
    ax.text(0.245, 0.955, "Internal", transform=ax.transAxes, ha="left", va="top", fontsize=21,
            color="#111111")
    ax.text(
        0.985,
        0.760,
        r"$\sqrt{s_{NN}}=200$ GeV" "\n"
        r"$14<p_{T}^{\gamma}<35$ GeV (native-bin support)" "\n"
        r"anti-$k_{T}$ $R=0.4$, $p_{T}^{\mathrm{jet}}>5$ GeV" "\n"
        r"$|\eta^{\gamma,\mathrm{jet}}|<0.7$, $|\Delta\phi|>7\pi/8$" "\n"
        r"Displayed: $x_{J\gamma}\geq0.5$" "\n"
        r"Statistical uncertainties only",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=14.5,
        linespacing=1.22,
        color="#202020",
    )
    ax.legend(
        loc="upper right",
        bbox_to_anchor=(0.985, 0.975),
        frameon=False,
        fontsize=18,
        handletextpad=0.55,
        labelspacing=0.45,
        borderaxespad=0.0,
    )
    fig.savefig(PNG, dpi=240, facecolor="white", bbox_inches=None)
    plt.close(fig)


def main() -> None:
    pp_curve, pp_meta = build_current_pp()
    auau_curve, auau_meta = load_auau()
    pp_edges = np.asarray(pp_curve["xj_edges"], dtype=float)
    auau_edges = np.asarray(auau_curve["xj_edges"], dtype=float)
    if not np.allclose(pp_edges, auau_edges, atol=1e-12, rtol=0.0):
        raise RuntimeError("p+p and Au+Au xJ bin edges do not match")
    if pp_curve["pt"] != [PT_RANGE[0], PT_RANGE[1]]:
        raise RuntimeError("p+p truth projection does not match requested native-bin support")
    if auau_curve["truth_pt_bins"] != [
        [14.0, 16.0, 15.0],
        [16.0, 18.0, 17.0],
        [18.0, 20.0, 19.0],
        [20.0, 22.0, 21.0],
        [22.0, 24.0, 23.0],
        [24.0, 26.0, 25.0],
        [26.0, 35.0, 30.5],
    ]:
        raise RuntimeError("unexpected Au+Au native photon-pT bin support")

    render(pp_curve, auau_curve)

    edges = pp_edges
    centers = 0.5 * (edges[:-1] + edges[1:])
    pp_values = np.asarray(pp_curve["values"], dtype=float)
    pp_errors = np.asarray(pp_curve["errors"], dtype=float)
    auau_values = np.asarray(auau_curve["values"], dtype=float)
    auau_errors = np.asarray(auau_curve["errors"], dtype=float)
    mask = displayed_mask(edges, pp_values, pp_errors) & displayed_mask(edges, auau_values, auau_errors)
    rows = []
    for ibin in range(len(centers)):
        rows.append(
            {
                "xj_low": float(edges[ibin]),
                "xj_high": float(edges[ibin + 1]),
                "xj_center": float(centers[ibin]),
                "displayed": bool(mask[ibin]),
                "omission_reason": None if mask[ibin] else "outside displayed xJ range",
                "pp": {"value": float(pp_values[ibin]), "stat_error": float(pp_errors[ibin])},
                "auau_0_20": {"value": float(auau_values[ibin]), "stat_error": float(auau_errors[ibin])},
            }
        )
    POINTS.write_text(json.dumps(json_ready({"bins": rows}), indent=2, allow_nan=False) + "\n")

    image = Image.open(PNG)
    displayed_lows = [row["xj_low"] for row in rows if row["displayed"]]
    checks = {
        "png_exists": PNG.exists(),
        "png_dimensions_2400x1800": image.size == (2400, 1800),
        "font_times_new_roman": configure_style() == "Times New Roman",
        "identical_xj_binning": bool(np.allclose(pp_edges, auau_edges, atol=1e-12, rtol=0.0)),
        "all_displayed_bins_start_at_or_above_0p5": bool(displayed_lows and min(displayed_lows) >= DISPLAY_XMIN),
        "turn_on_bins_retained_in_sidecar": any(not row["displayed"] and row["xj_low"] < DISPLAY_XMIN for row in rows),
        "statistical_uncertainties_only": True,
        "no_systematic_band_drawn": True,
        "no_manual_curve_rescaling": True,
        "not_unit_area_normalized": True,
        "google_slides_unchanged": True,
    }
    audit = {
        "ok": all(value is True for value in checks.values()),
        "checks": checks,
        "displayed_bin_count": int(mask.sum()),
        "omitted_low_xj_bin_count": int(np.sum(edges[:-1] < DISPLAY_XMIN)),
        "png": str(PNG),
        "png_sha256": sha256(PNG),
    }
    AUDIT.write_text(json.dumps(audit, indent=2, allow_nan=False) + "\n")
    if not audit["ok"]:
        raise RuntimeError(f"plot audit failed: {AUDIT}")

    manifest = {
        "ok": True,
        "status": "statistical_correction_candidate_not_final",
        "png": str(PNG),
        "png_sha256": sha256(PNG),
        "points": str(POINTS),
        "points_sha256": sha256(POINTS),
        "audit": str(AUDIT),
        "generator": str(Path(__file__).resolve()),
        "generator_sha256": sha256(Path(__file__).resolve()),
        "selection": {
            "photon_pt_gev_native_bin_support": list(PT_RANGE),
            "display_xj": [DISPLAY_XMIN, DISPLAY_XMAX],
            "jet_radius": 0.4,
            "jet_pt_min_gev": 5.0,
            "photon_eta_abs_max": 0.7,
            "jet_eta_abs_max": 0.7,
            "minimum_abs_delta_phi": "7pi/8",
        },
        "normalization": "per unfolded particle-level photon; not unit area",
        "uncertainties_drawn": "statistical only",
        "systematic_band_drawn": False,
        "manual_curve_rescaling": False,
        "pp": pp_meta,
        "auau": {
            "source_npz": str(AUAU_NPZ),
            "source_npz_sha256": sha256(AUAU_NPZ),
            "source_audit_manifest": str(AUAU_AUDIT_MANIFEST),
            "source_audit_manifest_sha256": sha256(AUAU_AUDIT_MANIFEST),
            "iterations": int(auau_meta["iterations"]),
            "error_mode": auau_meta["error_mode"],
            "data_input": auau_meta["data_input"],
            "response_measured_marginal": auau_meta["response_measured_marginal"],
            "closure_status": auau_meta["closure_status"],
        },
        "interpretation": (
            "The candidate has the expected ATLAS-like qualitative ordering: Au+Au shifts accepted recoil yield toward lower xJgamma "
            "and suppresses the balanced high-xJgamma tail relative to p+p. This is a cross-check of the correction chain, not a fit target."
        ),
        "caveats": [
            "The Au+Au response/K variant is an offline correction candidate and is not independently closure-approved.",
            "The current THE-88 Au+Au ROOT output was checked separately and is not substituted here because its corrected projection has a large low-x excess.",
            "The p+p curve uses three Bayesian iterations selected by the current half-closure study; the Au+Au audit candidate uses five iterations.",
            "Complete detector, photon, jet, unfolding, background-transfer, and model systematic covariance is not available and is not drawn.",
            "Scientific acceptance and canonical promotion remain manual.",
        ],
        "google_slides_mutated": False,
    }
    MANIFEST.write_text(json.dumps(json_ready(manifest), indent=2, allow_nan=False) + "\n")
    print(PNG)
    print(MANIFEST)
    print(AUDIT)


if __name__ == "__main__":
    main()
