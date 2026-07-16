#!/usr/bin/env python3
"""Local forensic comparison for PPG12 Fig. 13 E11/E33 vs current full pp.

This is intentionally local/read-only for inputs.  It compares the validated
PPG12 SDCC Fig. 13 JSON extraction against the completed RecoilJets pp ROOT
shower-shape histograms, including nearby ET-window variants, and writes a
compact evidence JSON/CSV for deciding whether a new pp pass needs exact
trigger/run-coverage instrumentation.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Iterable

import numpy as np
import ROOT


BASE = Path("dataOutput/ppg12PhotonYield/ppg12_photon_yield_v1_data_20260620")
DEFAULT_PPG12_JSON = (
    BASE
    / "shower_shape_reference_validation/fig13_e11_e33/"
    / "ppg12_sdcc_fig13_e11_to_e33_histograms.json"
)
DEFAULT_CURRENT_ROOT = Path(
    "InputFiles/pp24/ppg12_photon_yield_v1_data_20260620/pp/"
    "RecoilJets_pp_ALL_jetMinPtScan_dphiScan_vz60_isoR40_isSliding_"
    "preselectionNewPPG12_tightNewPPG12_nonTightNewPPG12.root"
)
DEFAULT_OUT_JSON = (
    BASE
    / "shower_shape_reference_validation/fig13_e11_e33/"
    / "ppg12_current_fullpp_e11e33_local_forensic_summary.json"
)
DEFAULT_OUT_CSV = (
    BASE
    / "shower_shape_reference_validation/fig13_e11_e33/"
    / "ppg12_current_fullpp_e11e33_local_forensic_windows.csv"
)

HIST_KEY = "e11e33"
TAG_KEY = "inclusive"
DIR_HINT = "PPG12_scaledtrigger30"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ppg12-json", type=Path, default=DEFAULT_PPG12_JSON)
    parser.add_argument("--current-root", type=Path, default=DEFAULT_CURRENT_ROOT)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--out-csv", type=Path, default=DEFAULT_OUT_CSV)
    return parser.parse_args()


def load_ppg12(path: Path) -> dict[str, np.ndarray | str | float]:
    with path.open() as f:
        payload = json.load(f)
    data = payload["data"]
    y = np.asarray(data["values"], dtype=float)
    err = np.asarray(data["errors"], dtype=float)
    nonzero_err = err > 0
    n_eff_values = np.divide(
        y[nonzero_err],
        np.square(err[nonzero_err]),
        out=np.full(np.count_nonzero(nonzero_err), np.nan),
        where=err[nonzero_err] > 0,
    )
    n_eff_values = n_eff_values[np.isfinite(n_eff_values)]
    return {
        "x": np.asarray(data["centers"], dtype=float),
        "y": y,
        "err": err,
        "n_eff_visible": float(np.median(n_eff_values)) if len(n_eff_values) else float("nan"),
        "source": payload.get("source_paths", {}).get("data", ""),
        "histname": payload.get("histname", ""),
        "chi2_ndf": payload.get("chi2_ndf", float("nan")),
    }


def get_hist(root_file: ROOT.TFile, directory: str, low: int, high: int):
    name = f"{directory}/h_ss_{HIST_KEY}_{TAG_KEY}_pT_{low}_{high}"
    hist = root_file.Get(name)
    if not hist:
        return None, name
    return hist, name


def find_candidate_dirs(root_file: ROOT.TFile) -> list[str]:
    dirs: list[str] = []
    for key in root_file.GetListOfKeys():
        obj = key.ReadObj()
        if not obj.InheritsFrom("TDirectory"):
            continue
        directory = key.GetName()
        d = root_file.Get(directory)
        if not d:
            continue
        found = False
        for subkey in d.GetListOfKeys():
            name = subkey.GetName()
            if name.startswith(f"h_ss_{HIST_KEY}_{TAG_KEY}_pT_"):
                found = True
                break
        if found:
            dirs.append(directory)
    return sorted(dirs)


def find_matching_keys(root_file: ROOT.TFile, needles: Iterable[str], limit: int = 40) -> list[str]:
    needle_list = [needle.lower() for needle in needles]
    matches: list[str] = []

    def walk(directory, prefix: str = "") -> None:
        if len(matches) >= limit:
            return
        for key in directory.GetListOfKeys():
            name = key.GetName()
            full_name = f"{prefix}/{name}" if prefix else name
            low = full_name.lower()
            if any(needle in low for needle in needle_list):
                matches.append(full_name)
                if len(matches) >= limit:
                    return
            obj = key.ReadObj()
            if obj and obj.InheritsFrom("TDirectory"):
                walk(obj, full_name)
                if len(matches) >= limit:
                    return

    walk(root_file)
    return matches


def read_counter(root_file: ROOT.TFile, name: str) -> dict[str, float] | None:
    obj = root_file.Get(name)
    if not obj:
        return None
    return {
        "integral": float(obj.Integral()) if hasattr(obj, "Integral") else float("nan"),
        "entries": float(obj.GetEntries()) if hasattr(obj, "GetEntries") else float("nan"),
    }


def combine_window(root_file: ROOT.TFile, directory: str, edges: Iterable[int]):
    edge_list = list(edges)
    if len(edge_list) < 2:
        raise ValueError("Need at least two edges")
    combined = None
    raw_integral = 0.0
    missing: list[str] = []
    used: list[str] = []
    for low, high in zip(edge_list[:-1], edge_list[1:]):
        hist, name = get_hist(root_file, directory, low, high)
        if not hist:
            missing.append(name)
            continue
        raw_integral += float(hist.Integral())
        used.append(name)
        if combined is None:
            combined = hist.Clone(f"h_{directory}_{low}_{high}_tmp".replace("/", "_"))
            combined.SetDirectory(0)
        else:
            combined.Add(hist)
    if combined is None:
        return None, raw_integral, used, missing
    return combined, raw_integral, used, missing


def normalized_004(hist) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    rebinned = hist.Rebin(4, f"{hist.GetName()}_rebin004")
    rebinned.SetDirectory(0)
    axis = rebinned.GetXaxis()
    first = axis.FindBin(0.000001)
    last = axis.FindBin(0.999999)
    norm = float(rebinned.Integral(first, last))
    if norm <= 0:
        raise RuntimeError(f"Histogram {hist.GetName()} has no visible 0-1 integral")
    xs: list[float] = []
    ys: list[float] = []
    errs: list[float] = []
    for ibin in range(first, last + 1):
        xs.append(0.5 * (axis.GetBinLowEdge(ibin) + axis.GetBinUpEdge(ibin)))
        ys.append(float(rebinned.GetBinContent(ibin)) / norm)
        errs.append(float(rebinned.GetBinError(ibin)) / norm)
    return np.asarray(xs), np.asarray(ys), np.asarray(errs), norm


def metrics(ref: dict[str, np.ndarray | str | float], x, y, err) -> dict[str, float]:
    ref_x = ref["x"]
    ref_y = ref["y"]
    ref_err = ref["err"]
    if len(x) != len(ref_x) or not np.allclose(x, ref_x, atol=1e-6):
        raise RuntimeError("Bin centers do not match PPG12 reference")

    ratio = np.divide(y, ref_y, out=np.full_like(y, np.nan), where=ref_y > 0)
    stable = np.isfinite(ratio) & (ref_y > 0.01)
    core = stable & (ref_x >= 0.18) & (ref_x <= 0.90)
    high_tail = stable & (ref_x >= 0.58) & (ref_x <= 0.90)
    peak = stable & (ref_x >= 0.26) & (ref_x <= 0.42)

    variance = np.square(err) + np.square(ref_err)
    chi_mask = stable & (variance > 0)
    chi2 = float(np.sum(np.square(y[chi_mask] - ref_y[chi_mask]) / variance[chi_mask]))
    ndf = int(np.count_nonzero(chi_mask))

    def mean_abs(mask):
        return float(np.mean(np.abs(ratio[mask] - 1.0))) if np.any(mask) else float("nan")

    def mean_ratio(mask):
        return float(np.mean(ratio[mask])) if np.any(mask) else float("nan")

    def max_abs(mask):
        return float(np.max(np.abs(ratio[mask] - 1.0))) if np.any(mask) else float("nan")

    return {
        "chi2_stable": chi2,
        "ndf_stable": float(ndf),
        "chi2_ndf_stable": chi2 / ndf if ndf else float("nan"),
        "mean_abs_ratio_minus1_stable": mean_abs(stable),
        "max_abs_ratio_minus1_stable": max_abs(stable),
        "mean_ratio_stable": mean_ratio(stable),
        "mean_abs_ratio_minus1_core_0p18_0p90": mean_abs(core),
        "mean_ratio_peak_0p26_0p42": mean_ratio(peak),
        "mean_ratio_high_tail_0p58_0p90": mean_ratio(high_tail),
        "ratio_at_0p34": float(ratio[np.argmin(np.abs(ref_x - 0.34))]),
        "ratio_at_0p82": float(ratio[np.argmin(np.abs(ref_x - 0.82))]),
        "ratio_at_0p90": float(ratio[np.argmin(np.abs(ref_x - 0.90))]),
    }


def main() -> None:
    args = parse_args()
    ref = load_ppg12(args.ppg12_json)
    root_file = ROOT.TFile.Open(str(args.current_root))
    if not root_file or root_file.IsZombie():
        raise RuntimeError(f"Could not open current ROOT: {args.current_root}")

    dirs = find_candidate_dirs(root_file)
    current_root_audit = {
        "candidate_dirs_with_h_ss": dirs,
        "event_counter": read_counter(
            root_file,
            f"{DIR_HINT}/cnt_{DIR_HINT}",
        ),
        "vertex_counter": read_counter(root_file, f"{DIR_HINT}/h_vertexZ"),
        "do_not_scale_probe_counter": read_counter(
            root_file,
            f"{DIR_HINT}/h_maxEnergyClus_NewTriggerFilling_doNotScale_{DIR_HINT}",
        ),
        "matching_keys_run_segment_trigger_audit": find_matching_keys(
            root_file,
            ["scaledtrigger", "trigger30", "bit30", "runnumber", "run_number", "segment"],
        ),
    }
    windows = {
        "ppg12_nominal_22_28_from_three_2gev_bins": [22, 24, 26, 28],
        "neighbor_low_20_28": [20, 22, 24, 26, 28],
        "neighbor_high_22_32": [22, 24, 26, 28, 32],
        "wide_20_32": [20, 22, 24, 26, 28, 32],
        "lower_half_22_26": [22, 24, 26],
        "upper_half_24_28": [24, 26, 28],
        "central_24_26": [24, 26],
    }

    rows = []
    for directory in dirs:
        for label, edges in windows.items():
            hist, raw_integral, used, missing = combine_window(root_file, directory, edges)
            row = {
                "directory": directory,
                "window_label": label,
                "edges": edges,
                "raw_integral": raw_integral,
                "used": used,
                "missing": missing,
            }
            if hist is not None and not missing:
                x, y, err, visible_norm = normalized_004(hist)
                row["visible_norm_0_to_1"] = visible_norm
                row.update(metrics(ref, x, y, err))
            rows.append(row)

    root_file.Close()

    sortable = [
        r for r in rows
        if "mean_abs_ratio_minus1_stable" in r
        and math.isfinite(float(r["mean_abs_ratio_minus1_stable"]))
    ]
    sortable.sort(key=lambda r: float(r["mean_abs_ratio_minus1_stable"]))
    best = sortable[0] if sortable else {}
    ppg12_n_eff = float(ref["n_eff_visible"])
    current_visible = float(best.get("visible_norm_0_to_1", float("nan")))
    entry_ratio = current_visible / ppg12_n_eff if ppg12_n_eff > 0 else float("nan")

    payload = {
        "inputs": {
            "ppg12_json": str(args.ppg12_json.resolve()),
            "current_root": str(args.current_root.resolve()),
            "ppg12_data_source": ref["source"],
            "ppg12_histname": ref["histname"],
            "ppg12_json_chi2_ndf": ref["chi2_ndf"],
            "ppg12_effective_visible_entries_from_errors": ppg12_n_eff,
        },
        "summary": {
            "best_local_window": best,
            "current_visible_entries_for_best_window": current_visible,
            "current_visible_entries_over_ppg12_effective_visible_entries": entry_ratio,
            "interpretation": (
                "The nearest current full-pp diagnostic is the exact 22-28 GeV "
                "sum of h_ss_e11e33_inclusive pT bins. Its normalized shape is "
                "close to PPG12, and nearby ET windows are worse, so the residual "
                "is not explained by a simple ET-window choice. The current ROOT "
                "does not contain direct scaledtrigger[30] or run/segment audit "
                "objects, so the remaining non-1:1 difference cannot be certified "
                "as shower-shape arithmetic rather than event/trigger/run coverage."
            ),
        },
        "ppg12_contract_from_code": {
            "source_macro": "ppg12codeGit/efficiencytool/ShowerShapeCheck.C",
            "config": "ppg12codeGit/efficiencytool/config_showershape_0rad.yaml",
            "data_glob": "/sphenix/user/shuhangli/ppg12/anatreemaker/macro_maketree/data/ana521/condorout/part_*_with_bdt_split.root",
            "tree": "slimtree",
            "trigger_used": 30,
            "eta": "abs(cluster_Eta) < 0.7",
            "vertex_cut": "abs(vertexz) < 60",
            "pt_window": "22 < cluster_Et < 28",
            "fill": "h2d_e11_to_e33 cut0 before common/tight/non-tight cuts",
        },
        "recoiljets_current_contract_from_code": {
            "source": "src/RecoilJets.cc fillSSSpectra",
            "current_saved_directory_hint": DIR_HINT,
            "trigger_gate": "TriggerAnalyzer::didTriggerFire(\"Photon 4 GeV + MBD NS >= 1\")",
            "fill": "h_ss_e11e33_inclusive_pT_* before NPB/preselection/tight classification",
            "note": "Current full-pp ROOT stores only the friendly trigger directory listed in candidate_dirs; it does not store a scaledtrigger[30] audit directory.",
        },
        "current_root_audit": current_root_audit,
        "ranked_windows": sortable,
        "all_rows": rows,
    }

    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(payload, indent=2, sort_keys=True))

    fieldnames = [
        "directory",
        "window_label",
        "edges",
        "raw_integral",
        "visible_norm_0_to_1",
        "mean_abs_ratio_minus1_stable",
        "max_abs_ratio_minus1_stable",
        "chi2_ndf_stable",
        "mean_ratio_stable",
        "mean_ratio_peak_0p26_0p42",
        "mean_ratio_high_tail_0p58_0p90",
        "ratio_at_0p34",
        "ratio_at_0p82",
        "ratio_at_0p90",
        "missing",
    ]
    with args.out_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fieldnames})

    print(f"Wrote {args.out_json}")
    print(f"Wrote {args.out_csv}")
    print("Top local windows by mean |ratio-1| in stable bins:")
    for row in sortable[:8]:
        print(
            f"  {row['directory']} {row['window_label']}: "
            f"mean_abs={row['mean_abs_ratio_minus1_stable']:.4f}, "
            f"tail_mean={row['mean_ratio_high_tail_0p58_0p90']:.4f}, "
            f"raw={row['raw_integral']:.0f}"
        )


if __name__ == "__main__":
    main()
