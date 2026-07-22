#!/usr/bin/env python3
"""Render the final THE97 pp-data purity comparison to the PPG12 SDCC result.

The current pp-data and photon+jet-simulation inputs are resolved through the
artifact registry.  The plot uses the PPG12 ABCD quadratic closure with the
contract-locked TRandom3(42), 20,000-toy uncertainty treatment.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path

import ROOT

from make_ppg12_fig3_purity_sim_sdcc_vs_current_overlay import (
    effective_count,
    leak_ratio,
    ppg12_toy_estimate,
)


ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetEndErrorSize(3)

REPO = Path(__file__).resolve().parents[3]
CONTRACT_ID = "the97_ppg12_data_purity_leakage_v1"
EXPECTED_EDGES = (10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 32.0, 36.0)

DATA_SAMPLE = "pp_data_merged"
DATA_CAMPAIGN = "the97_ppg12_final_accepted_triple_full_20260714_1550"
DATA_SHA256 = "5854abfa21629482ca3aef87c421b790945c093ca97dcee6093a7d6af53f9a6f"
SIM_SAMPLE = "pp_sim_photonjet_merged"
SIM_CAMPAIGN = "the97_ppg12_si_contract_restore_full_20260715_1420"
SIM_SHA256 = "c51e9d140dfa309d4bf22b9ebce300c80b9b121c27d631d1343f93da26abfd3e"

CURRENT_ROOT = REPO / "dataOutput/current_recoiljets_artifacts/current"
REFERENCE_CSV = REPO / "dataOutput/ppg12Parity/the76_ppg12_parity_full_20260701_003024/fig29_purity_datathief_audit/ppg12_fig29_purity_sdcc_points.csv"
OUT_DIR = REPO / "dataOutput/ppg12Parity/the97_ppg12_final_accepted_triple_full_20260714_1550/final_pp_data_canonical_20260717/purity_overlay_current"

DATA_DIR = "PPG12_scaledtrigger30"
SIM_DIR = "SIM"
DATA_NAMES = {
    "A": "h_tight_iso_cluster_0",
    "B": "h_tight_noniso_cluster_0",
    "C": "h_nontight_iso_cluster_0",
    "D": "h_nontight_noniso_cluster_0",
}
SIM_NAMES = {
    "A": "h_tight_iso_cluster_signal_0",
    "B": "h_tight_noniso_cluster_signal_0",
    "C": "h_nontight_iso_cluster_signal_0",
    "D": "h_nontight_noniso_cluster_signal_0",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def resolve_current(sample: str, campaign: str, expected_sha256: str) -> tuple[Path, dict[str, object]]:
    pointer_path = CURRENT_ROOT / sample / "current.json"
    pointer = json.loads(pointer_path.read_text())
    if pointer.get("sample_key") != sample:
        raise RuntimeError(f"current pointer sample mismatch: {pointer_path}")
    if pointer.get("campaign_tag") != campaign:
        raise RuntimeError(
            f"current pointer campaign mismatch for {sample}: "
            f"{pointer.get('campaign_tag')} != {campaign}"
        )
    roots = pointer.get("root_paths")
    if not isinstance(roots, list) or len(roots) != 1:
        raise RuntimeError(f"current pointer must resolve exactly one ROOT: {pointer_path}")
    root_path = Path(str(roots[0]))
    if not root_path.is_file():
        raise RuntimeError(f"registered current ROOT is missing: {root_path}")
    observed_sha256 = sha256(root_path)
    if observed_sha256 != expected_sha256:
        raise RuntimeError(
            f"registered current ROOT SHA mismatch for {sample}: "
            f"{observed_sha256} != {expected_sha256}"
        )
    return root_path, pointer


def open_root(path: Path) -> ROOT.TFile:
    root_file = ROOT.TFile.Open(str(path), "READ")
    if not root_file or root_file.IsZombie() or root_file.TestBit(ROOT.TFile.kRecovered):
        raise RuntimeError(f"invalid ROOT input: {path}")
    return root_file


def require_hist(root_file: ROOT.TFile, path: str) -> ROOT.TH1:
    hist = root_file.Get(path)
    if not hist or not hist.InheritsFrom("TH1"):
        raise RuntimeError(f"missing TH1 input: {path}")
    clone = hist.Clone(path.replace("/", "_") + "_final_plot")
    clone.SetDirectory(0)
    return clone


def histogram_edges(hist: ROOT.TH1) -> tuple[float, ...]:
    axis = hist.GetXaxis()
    return tuple(float(axis.GetBinLowEdge(index)) for index in range(1, hist.GetNbinsX() + 1)) + (
        float(axis.GetBinUpEdge(hist.GetNbinsX())),
    )


def assert_expected_binning(histograms: list[ROOT.TH1]) -> None:
    for hist in histograms:
        observed = histogram_edges(hist)
        if len(observed) != len(EXPECTED_EDGES) or any(
            abs(left - right) > 1.0e-9 for left, right in zip(observed, EXPECTED_EDGES)
        ):
            raise RuntimeError(f"unexpected ET binning in {hist.GetName()}: {observed}")


def load_reference(path: Path) -> dict[str, list[dict[str, float]]]:
    result: dict[str, list[dict[str, float]]] = {"raw": [], "leakage_corrected": []}
    with path.open(newline="") as handle:
        for record in csv.DictReader(handle):
            series = record["series"]
            if series not in result:
                continue
            result[series].append(
                {
                    key: float(record[key])
                    for key in ("x", "purity", "ex_low", "ex_high", "ey_low", "ey_high")
                }
            )
    if any(len(points) != 11 for points in result.values()):
        raise RuntimeError(f"PPG12 reference must contain 11 raw and 11 corrected points: {path}")
    for points in result.values():
        points.sort(key=lambda row: row["x"])
    return result


def calculate_current(data_root: Path, sim_root: Path) -> list[dict[str, float]]:
    data_file = open_root(data_root)
    sim_file = open_root(sim_root)
    try:
        data = {key: require_hist(data_file, f"{DATA_DIR}/{name}") for key, name in DATA_NAMES.items()}
        signal = {key: require_hist(sim_file, f"{SIM_DIR}/{name}") for key, name in SIM_NAMES.items()}
    finally:
        data_file.Close()
        sim_file.Close()

    assert_expected_binning(list(data.values()) + list(signal.values()))
    leakage = {
        "B": leak_ratio(signal["B"], signal["A"], "final_cB"),
        "C": leak_ratio(signal["C"], signal["A"], "final_cC"),
        "D": leak_ratio(signal["D"], signal["A"], "final_cD"),
    }
    rng = ROOT.TRandom3(42)
    rows: list[dict[str, float]] = []
    for index in range(1, 12):
        values = tuple(float(data[key].GetBinContent(index)) for key in "ABCD")
        counts = tuple(effective_count(data[key], index) for key in "ABCD")
        leak = tuple(float(leakage[key].GetBinContent(index)) for key in "BCD")
        leak_errors = tuple(float(leakage[key].GetBinError(index)) for key in "BCD")
        raw, raw_error, corrected, corrected_error, diagnostics = ppg12_toy_estimate(
            rng,
            values,
            counts,
            leak,
            leak_errors,
            f"the97_final_pp_bin_{index}",
        )
        rows.append(
            {
                "pt_lo": EXPECTED_EDGES[index - 1],
                "pt_hi": EXPECTED_EDGES[index],
                "center": 0.5 * (EXPECTED_EDGES[index - 1] + EXPECTED_EDGES[index]),
                "half_width": 0.5 * (EXPECTED_EDGES[index] - EXPECTED_EDGES[index - 1]),
                "A": values[0],
                "B": values[1],
                "C": values[2],
                "D": values[3],
                "n_eff_A": counts[0],
                "n_eff_B": counts[1],
                "n_eff_C": counts[2],
                "n_eff_D": counts[3],
                "cB": leak[0],
                "cC": leak[1],
                "cD": leak[2],
                "cB_error": leak_errors[0],
                "cC_error": leak_errors[1],
                "cD_error": leak_errors[2],
                "current_raw": raw,
                "current_raw_error": raw_error,
                "current_corrected": corrected,
                "current_corrected_error": corrected_error,
                **diagnostics,
            }
        )
    return rows


def merge_reference(
    current: list[dict[str, float]], reference: dict[str, list[dict[str, float]]]
) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    for index, base in enumerate(current):
        row = dict(base)
        for series, current_key, current_error_key in (
            ("raw", "current_raw", "current_raw_error"),
            ("leakage_corrected", "current_corrected", "current_corrected_error"),
        ):
            point = reference[series][index]
            if abs(point["x"] - base["center"]) > 1.0e-9:
                raise RuntimeError(f"PPG12/current x mismatch in {series} bin {index + 1}")
            prefix = "raw" if series == "raw" else "corrected"
            ppg12 = point["purity"]
            current_value = row[current_key]
            current_error = row[current_error_key]
            ratio = current_value / ppg12
            ratio_error_low = math.hypot(
                current_error / ppg12,
                current_value * point["ey_high"] / (ppg12 * ppg12),
            )
            ratio_error_high = math.hypot(
                current_error / ppg12,
                current_value * point["ey_low"] / (ppg12 * ppg12),
            )
            row.update(
                {
                    f"ppg12_{prefix}": ppg12,
                    f"ppg12_{prefix}_error_low": point["ey_low"],
                    f"ppg12_{prefix}_error_high": point["ey_high"],
                    f"current_over_ppg12_{prefix}": ratio,
                    f"current_over_ppg12_{prefix}_error_low": ratio_error_low,
                    f"current_over_ppg12_{prefix}_error_high": ratio_error_high,
                }
            )
        rows.append(row)
    return rows


def symmetric_graph(
    name: str,
    rows: list[dict[str, float]],
    value_key: str,
    error_key: str,
    color: int,
    marker: int,
) -> ROOT.TGraphErrors:
    graph = ROOT.TGraphErrors(len(rows))
    graph.SetName(name)
    graph.SetMarkerColor(color)
    graph.SetLineColor(color)
    graph.SetMarkerStyle(marker)
    graph.SetMarkerSize(1.15)
    graph.SetLineWidth(2)
    for index, row in enumerate(rows):
        graph.SetPoint(index, row["center"], row[value_key])
        graph.SetPointError(index, row["half_width"], row[error_key])
    return graph


def asymmetric_graph(
    name: str,
    rows: list[dict[str, float]],
    value_key: str,
    error_low_key: str,
    error_high_key: str,
    color: int,
    marker: int,
) -> ROOT.TGraphAsymmErrors:
    graph = ROOT.TGraphAsymmErrors(len(rows))
    graph.SetName(name)
    graph.SetMarkerColor(color)
    graph.SetLineColor(color)
    graph.SetMarkerStyle(marker)
    graph.SetMarkerSize(1.15)
    graph.SetLineWidth(2)
    for index, row in enumerate(rows):
        graph.SetPoint(index, row["center"], row[value_key])
        graph.SetPointError(
            index,
            row["half_width"],
            row["half_width"],
            row[error_low_key],
            row[error_high_key],
        )
    return graph


def render(rows: list[dict[str, float]], output: Path) -> tuple[float, float]:
    canvas = ROOT.TCanvas("the97_current_pp_purity", "", 860, 900)
    canvas.SetFillColor(ROOT.kWhite)
    top = ROOT.TPad("purity_top", "", 0.0, 0.34, 1.0, 1.0)
    bottom = ROOT.TPad("purity_ratio", "", 0.0, 0.0, 1.0, 0.34)
    for pad in (top, bottom):
        pad.SetFillColor(ROOT.kWhite)
        pad.SetLeftMargin(0.14)
        pad.SetRightMargin(0.045)
        pad.SetTicks(1, 1)
    top.SetTopMargin(0.045)
    top.SetBottomMargin(0.018)
    bottom.SetTopMargin(0.025)
    bottom.SetBottomMargin(0.27)
    top.Draw()
    bottom.Draw()

    ppg12_corrected = asymmetric_graph(
        "g_ppg12_leakage_corrected",
        rows,
        "ppg12_corrected",
        "ppg12_corrected_error_low",
        "ppg12_corrected_error_high",
        ROOT.kAzure + 6,
        25,
    )
    current_corrected = symmetric_graph(
        "g_current_output_corrected",
        rows,
        "current_corrected",
        "current_corrected_error",
        ROOT.kBlue + 1,
        21,
    )
    ppg12_raw = asymmetric_graph(
        "g_ppg12_raw",
        rows,
        "ppg12_raw",
        "ppg12_raw_error_low",
        "ppg12_raw_error_high",
        ROOT.kGray + 2,
        24,
    )
    current_raw = symmetric_graph(
        "g_current_output_raw",
        rows,
        "current_raw",
        "current_raw_error",
        ROOT.kBlack,
        20,
    )
    corrected_ratio = asymmetric_graph(
        "g_current_over_ppg12_corrected",
        rows,
        "current_over_ppg12_corrected",
        "current_over_ppg12_corrected_error_low",
        "current_over_ppg12_corrected_error_high",
        ROOT.kBlue + 1,
        21,
    )
    raw_ratio = asymmetric_graph(
        "g_current_over_ppg12_raw",
        rows,
        "current_over_ppg12_raw",
        "current_over_ppg12_raw_error_low",
        "current_over_ppg12_raw_error_high",
        ROOT.kBlack,
        20,
    )

    top.cd()
    frame = ROOT.TH1F("purity_frame", "", 26, 10.0, 36.0)
    frame.SetDirectory(0)
    frame.SetStats(False)
    # Headroom raised from 1.20 so the legend can sit directly under the header
    # block.  The highest drawn point (E_T=34 GeV, 1.146 incl. error) then lands
    # at NDC 0.782, clear of the legend band below 0.80 for x < 0.835.
    frame.GetYaxis().SetRangeUser(0.25, 1.35)
    frame.GetYaxis().SetTitle("Purity")
    frame.GetYaxis().SetTitleSize(0.064)
    frame.GetYaxis().SetTitleOffset(0.88)
    frame.GetYaxis().SetLabelSize(0.048)
    frame.GetYaxis().SetNdivisions(508)
    frame.GetXaxis().SetLabelSize(0.0)
    frame.Draw("AXIS")
    ppg12_corrected.Draw("PZ SAME")
    current_corrected.Draw("PZ SAME")
    ppg12_raw.Draw("PZ SAME")
    current_raw.Draw("PZ SAME")

    label = ROOT.TLatex()
    label.SetNDC(True)
    label.SetTextFont(42)
    label.SetTextColor(ROOT.kBlack)
    # Held off the top frame edge and in from the left, at larger type.
    label.SetTextSize(0.052)
    label.DrawLatex(0.195, 0.880, "#it{#bf{sPHENIX}} Internal")
    label.SetTextSize(0.040)
    label.DrawLatex(0.195, 0.815, "p+p, #sqrt{s} = 200 GeV, |#eta^{#gamma}| < 0.7")

    # Sits directly beneath the header block, two columns, right edge held at
    # x=0.835 so it clears the tall E_T=34 GeV error bar on the far right.
    legend = ROOT.TLegend(0.185, 0.645, 0.835, 0.795)
    legend.SetBorderSize(0)
    legend.SetFillColor(ROOT.kWhite)
    legend.SetFillStyle(1001)
    legend.SetTextFont(42)
    legend.SetTextSize(0.034)
    legend.SetNColumns(2)
    legend.SetColumnSeparation(0.01)
    legend.SetMargin(0.16)
    legend.AddEntry(ppg12_corrected, "PPG12 leakage corrected", "p")
    legend.AddEntry(current_corrected, "This analysis, corrected", "p")
    legend.AddEntry(ppg12_raw, "PPG12 raw", "p")
    legend.AddEntry(current_raw, "This analysis, raw", "p")
    legend.Draw()
    top.RedrawAxis()

    ratio_extents: list[float] = []
    for row in rows:
        for prefix in ("corrected", "raw"):
            value = row[f"current_over_ppg12_{prefix}"]
            ratio_extents.extend(
                (
                    value - row[f"current_over_ppg12_{prefix}_error_low"],
                    value + row[f"current_over_ppg12_{prefix}_error_high"],
                )
            )
    ratio_low = max(0.0, min(0.90, min(ratio_extents) - 0.08))
    ratio_high = min(1.80, max(1.10, max(ratio_extents) + 0.08))

    bottom.cd()
    ratio_frame = ROOT.TH1F("purity_ratio_frame", "", 26, 10.0, 36.0)
    ratio_frame.SetDirectory(0)
    ratio_frame.SetStats(False)
    ratio_frame.GetYaxis().SetRangeUser(ratio_low, ratio_high)
    ratio_frame.GetXaxis().SetTitle("E_{T}^{#gamma,rec} [GeV]")
    ratio_frame.GetYaxis().SetTitle("This analysis / PPG12")
    ratio_frame.GetXaxis().SetTitleSize(0.100)
    ratio_frame.GetXaxis().SetTitleOffset(1.10)
    ratio_frame.GetXaxis().SetLabelSize(0.076)
    ratio_frame.GetYaxis().SetTitleSize(0.078)
    ratio_frame.GetYaxis().SetTitleOffset(0.72)
    ratio_frame.GetYaxis().SetLabelSize(0.061)
    ratio_frame.GetYaxis().SetNdivisions(505)
    ratio_frame.Draw("AXIS")
    unity = ROOT.TLine(10.0, 1.0, 36.0, 1.0)
    unity.SetLineColor(ROOT.kGray + 1)
    unity.SetLineStyle(7)
    unity.SetLineWidth(2)
    unity.Draw()
    corrected_ratio.Draw("PZ SAME")
    raw_ratio.Draw("PZ SAME")
    bottom.RedrawAxis()

    canvas.SaveAs(str(output))
    return ratio_low, ratio_high


def write_csv(path: Path, rows: list[dict[str, float]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    data_root, data_pointer = resolve_current(DATA_SAMPLE, DATA_CAMPAIGN, DATA_SHA256)
    sim_root, sim_pointer = resolve_current(SIM_SAMPLE, SIM_CAMPAIGN, SIM_SHA256)
    reference_sha256 = sha256(REFERENCE_CSV)
    current = calculate_current(data_root, sim_root)
    rows = merge_reference(current, load_reference(REFERENCE_CSV))

    png = OUT_DIR / "the97_pp_data_purity_ppg12_vs_current_fullstat.png"
    points_csv = OUT_DIR / "the97_pp_data_purity_ppg12_vs_current_fullstat_points.csv"
    manifest_path = OUT_DIR / "the97_pp_data_purity_ppg12_vs_current_fullstat_manifest.json"
    write_csv(points_csv, rows)
    ratio_ylim = render(rows, png)

    script_path = Path(__file__).resolve()
    solver_path = script_path.with_name("make_ppg12_fig3_purity_sim_sdcc_vs_current_overlay.py")
    manifest = {
        "schema_version": 1,
        "figure_status": "candidate_current_inputs",
        "contract_id": CONTRACT_ID,
        "full_stat_pp_data": True,
        "partial_coverage_label": False,
        "data": {
            "sample_key": DATA_SAMPLE,
            "current_entry_id": data_pointer["current_entry_id"],
            "campaign_tag": DATA_CAMPAIGN,
            "root": str(data_root),
            "sha256": DATA_SHA256,
            "namespace": DATA_DIR,
            "histograms": DATA_NAMES,
        },
        "signal_leakage": {
            "sample_key": SIM_SAMPLE,
            "current_entry_id": sim_pointer["current_entry_id"],
            "campaign_tag": SIM_CAMPAIGN,
            "root": str(sim_root),
            "sha256": SIM_SHA256,
            "namespace": SIM_DIR,
            "histograms": SIM_NAMES,
            "definition": "fully combined photon+jet _signal leakage; corrected SI plus unchanged valid DI",
        },
        "ppg12_reference": {
            "csv": str(REFERENCE_CSV),
            "sha256": reference_sha256,
            "source": "/sphenix/user/shuhangli/ppg12/efficiencytool/results/Photon_final_bdt_nom.root -> gpurity/gpurity_leak",
        },
        "estimator": {
            "solver": "PPG12 CalculatePhotonYield quadratic physical branch",
            "rng": "ROOT.TRandom3(42), one sequential stream across all 11 bins",
            "toys_per_bin": 20000,
            "data_toys": "effective-Poisson",
            "leakage_toys": "Gaussian cB/cC/cD throws",
            "bin_edges_gev": EXPECTED_EDGES,
        },
        "legend": [
            "PPG12 leakage corrected",
            "Current output corrected",
            "PPG12 raw",
            "Current output raw",
        ],
        "ratio_panel": "Current / PPG12 for leakage-corrected and raw purity",
        "ratio_ylim": ratio_ylim,
        "outputs": {
            "png": str(png),
            "png_sha256": sha256(png),
            "points_csv": str(points_csv),
            "points_csv_sha256": sha256(points_csv),
        },
        "code": {
            "renderer": str(script_path),
            "renderer_sha256": sha256(script_path),
            "ppg12_toy_solver": str(solver_path),
            "ppg12_toy_solver_sha256": sha256(solver_path),
        },
        "visual_qa": {
            "status": "pending_manual_inspection",
            "note": "Updated only after the rendered PNG is visually inspected.",
        },
        "slides_mutated": False,
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"png": str(png), "csv": str(points_csv), "manifest": str(manifest_path)}, indent=2))


if __name__ == "__main__":
    main()
