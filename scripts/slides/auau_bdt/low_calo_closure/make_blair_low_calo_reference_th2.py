#!/usr/bin/env python3
"""Reproduce Blair-facing low-calo centrality-vs-energy TH2 reference plots.

The input is the full-count event histogram JSON produced from the THE-32
score-cache audit. This script does not read candidate rows or recompute the
low-calo decision; it makes the final compact TH2 views directly from the saved
event-level retained+removed histogram bins.
"""

from __future__ import annotations

import argparse
from array import array
import json
from pathlib import Path

import ROOT


DEFAULT_HIST_JSON = Path(
    "dataOutput/auauTightBDTValidation/THE32_lowCaloDiagnosticClosure_20260603/"
    "the32_low_calo_full_count_histograms_v1.json"
)
DEFAULT_OUTDIR = Path(
    "dataOutput/auauTightBDTValidation/THE32_lowCaloDiagnosticClosure_20260603/"
    "blair_reference_20260702"
)


def load_payload(path: Path) -> dict:
    payload = json.loads(path.read_text())
    required = [
        "schema",
        "source_cache_list",
        "cache_count",
        "event_rows_after_per_cache_dedup",
        "raw_candidate_rows",
        "energy_edges",
        "panels",
        "cut_variable",
    ]
    missing = [key for key in required if key not in payload]
    if missing:
        raise SystemExit(f"{path} is missing required keys: {missing}")
    return payload


def root_array(values: list[float]) -> array:
    return array("d", [float(v) for v in values])


def make_hist(payload: dict, *, use_gev_axis: bool, hist_name: str) -> tuple[ROOT.TH2D, list[float], int, int]:
    panels = payload["panels"]
    x_edges = [float(panels[0]["cent_lo"])] + [float(panel["cent_hi"]) for panel in panels]
    log_edges = [float(v) for v in payload["energy_edges"]]
    y_edges = [(10.0**v - 1.0) for v in log_edges] if use_gev_axis else log_edges

    selection = str(payload.get("tower_selection", ""))
    title = "Centrality vs total calorimeter energy"
    if "get_isGood" in selection:
        title = "Centrality vs total calorimeter energy (good towers)"

    h = ROOT.TH2D(
        hist_name,
        title,
        len(x_edges) - 1,
        root_array(x_edges),
        len(y_edges) - 1,
        root_array(y_edges),
    )
    h.SetDirectory(0)
    h.GetXaxis().SetTitle("Centrality percentile [%]")
    if use_gev_axis:
        h.GetYaxis().SetTitle("E^{total}_{calo} = E_{CEMC}+E_{IHCal}+E_{OHCal} [GeV]")
    else:
        h.GetYaxis().SetTitle("log_{10}(E^{total}_{calo}+1)")
    h.GetZaxis().SetTitle("Event count / bin")

    histogrammed_total = 0
    nonzero_bins = 0
    for ix, panel in enumerate(panels, start=1):
        retained = panel["retained_hist"]
        removed = panel["removed_hist"]
        if len(retained) != len(removed) or len(retained) != len(y_edges) - 1:
            raise SystemExit(
                f"Panel {panel['cent_lo']}-{panel['cent_hi']} has inconsistent histogram length"
            )
        for iy, (keep, cut) in enumerate(zip(retained, removed), start=1):
            value = float(keep) + float(cut)
            if value <= 0.0:
                continue
            h.SetBinContent(ix, iy, value)
            histogrammed_total += int(round(value))
            nonzero_bins += 1
    h.SetEntries(nonzero_bins)
    return h, y_edges, histogrammed_total, nonzero_bins


def draw_hist(h: ROOT.TH2D, out_png: Path, *, log_z: bool = True) -> None:
    ROOT.gROOT.SetBatch(True)
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPalette(ROOT.kBird)
    canvas = ROOT.TCanvas("c_blair_low_calo", "", 1350, 900)
    canvas.SetRightMargin(0.15)
    canvas.SetLeftMargin(0.14)
    canvas.SetBottomMargin(0.13)
    canvas.SetTopMargin(0.08)
    canvas.SetLogz(bool(log_z))

    if log_z:
        h.SetMinimum(1.0)
    else:
        h.SetMinimum(0.0)
    h.GetXaxis().SetTitleSize(0.045)
    h.GetYaxis().SetTitleSize(0.045)
    h.GetZaxis().SetTitleSize(0.042)
    h.GetXaxis().SetLabelSize(0.040)
    h.GetYaxis().SetLabelSize(0.040)
    h.GetZaxis().SetLabelSize(0.036)
    h.GetYaxis().SetTitleOffset(1.35)
    h.GetZaxis().SetTitleOffset(1.25)
    h.Draw("COLZ")
    canvas.SaveAs(str(out_png))


def write_outputs(
    payload: dict,
    hist_json: Path,
    outdir: Path,
    *,
    use_gev_axis: bool,
    log_z: bool = True,
    output_tag: str = "",
) -> Path:
    outdir.mkdir(parents=True, exist_ok=True)
    suffix = "total_caloE_GeV_ROOT_kBird_clean_noline" if use_gev_axis else "log10_total_caloE_ROOT_kBird_clean"
    if not log_z:
        suffix = f"{suffix}_linearZ"
    if output_tag:
        suffix = f"{suffix}_{output_tag}"
    stem = f"blair_cent_vs_{suffix}_20260702"
    hist_name = "h_cent_vs_total_caloE_GeV_noline" if use_gev_axis else "h_cent_vs_log10_total_caloE"
    out_png = outdir / f"{stem}.png"
    out_root = outdir / f"{stem}.root"
    out_manifest = outdir / f"{stem}.json"

    h, y_edges, histogrammed_total, nonzero_bins = make_hist(
        payload, use_gev_axis=use_gev_axis, hist_name=hist_name
    )
    with ROOT.TFile(str(out_root), "RECREATE") as f:
        h.Write()
    draw_hist(h, out_png, log_z=log_z)

    manifest = {
        "png": str(out_png),
        "root": str(out_root),
        "hist_name": hist_name,
        "source_hist_json": str(hist_json),
        "source_cache_list": payload["source_cache_list"],
        "cache_count": payload["cache_count"],
        "raw_candidate_rows": payload["raw_candidate_rows"],
        "event_rows_after_per_cache_dedup": payload["event_rows_after_per_cache_dedup"],
        "histogrammed_event_count": histogrammed_total,
        "nonzero_bins": nonzero_bins,
        "centrality_domain": payload["centrality_domain"],
        "cut_variable": payload["cut_variable"],
        "tower_selection": payload.get("tower_selection"),
        "display_energy_axis": "GeV" if use_gev_axis else "log10(E_calo+1)",
        "y_axis_conversion": "E_calo_GeV = 10**log10(E_calo+1) - 1"
        if use_gev_axis
        else "none",
        "y_axis_edges_first_last": [float(y_edges[0]), float(y_edges[-1])],
        "palette": "ROOT.kBird / TColor kBird = 57",
        "z_axis_scale": "log" if log_z else "linear",
        "caveat": (
            "Embedded overlay/score-cache binned event sample used for the low-calo audit, "
            "not raw AuAu data and not a weighted physics mixture."
        ),
    }
    out_manifest.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return out_png


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--hist-json", type=Path, default=DEFAULT_HIST_JSON)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--linear-z", action="store_true", help="Use a linear color scale instead of log-z.")
    parser.add_argument("--gev-only", action="store_true", help="Only write the GeV-y-axis plot.")
    parser.add_argument("--output-tag", default="", help="Optional suffix tag added to output filenames.")
    args = parser.parse_args()

    payload = load_payload(args.hist_json)
    axes = (True,) if args.gev_only else (False, True)
    for use_gev_axis in axes:
        png = write_outputs(
            payload,
            args.hist_json,
            args.outdir,
            use_gev_axis=use_gev_axis,
            log_z=not args.linear_z,
            output_tag=args.output_tag,
        )
        print(png)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
