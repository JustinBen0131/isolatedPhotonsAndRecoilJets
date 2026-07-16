#!/usr/bin/env python3
"""Overlay PPG12 IAN Fig. 3 SIM purities with current RecoilJets SIM.

Top panel: truth, raw ABCD, and signal-leakage-corrected purity for the
authoritative PPG12 SDCC graph ROOT and the current RecoilJets inclusive-SIM
ROOT. Bottom panel: PPG12 SDCC / current for each definition.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path

import ROOT


ROOT.gROOT.SetBatch(True)

REPO = Path(__file__).resolve().parents[3]
DEFAULT_REFERENCE_ROOT = (
    REPO
    / "dataOutput/ppg12PhotonYield/ppg12_photon_yield_v1_data_20260620"
    / "reference_roots/fig3_purity_sim/Photon_final_bdt_nom_mc.root"
)
DEFAULT_CURRENT_JSON = (
    REPO
    / "dataOutput/current_recoiljets_artifacts/current"
    / "pp_sim_inclusivejet_merged/current.json"
)
DEFAULT_PHOTON_JSON = (
    REPO
    / "dataOutput/current_recoiljets_artifacts/current"
    / "pp_sim_photonjet_merged/current.json"
)
DEFAULT_OUTDIR = (
    REPO
    / "dataOutput/ppg12Parity/the97_ppg12_final_parity_full_20260709_2230"
    / "fig3_purity_sim"
)

SERIES = (
    {"key": "truth", "label": "truth", "reference": "g_purity_truth", "color": ROOT.kRed},
    {"key": "raw", "label": "ABCD w/o leakage correction", "reference": "gpurity", "color": ROOT.kBlack},
    {"key": "corrected", "label": "ABCD w/ leakage correction", "reference": "gpurity_leak", "color": ROOT.kBlue},
)


@dataclass(frozen=True)
class Point:
    x: float
    ex_low: float
    ex_high: float
    y: float
    ey_low: float
    ey_high: float


def open_root(path: Path) -> ROOT.TFile:
    f = ROOT.TFile.Open(str(path), "READ")
    if not f or f.IsZombie() or f.TestBit(ROOT.TFile.kRecovered):
        raise RuntimeError(f"invalid ROOT file: {path}")
    return f


def resolve_current_root(pointer: Path) -> Path:
    payload = json.loads(pointer.read_text())
    roots = payload.get("root_paths") or []
    if not roots:
        raise RuntimeError(f"no root_paths in {pointer}")
    path = Path(roots[0])
    if not path.is_absolute():
        path = REPO / path
    if not path.exists():
        raise FileNotFoundError(path)
    return path


def require_hist(directory: ROOT.TDirectory, name: str) -> ROOT.TH1:
    obj = directory.Get(name)
    if not obj or not obj.InheritsFrom("TH1"):
        raise RuntimeError(f"missing histogram SIM/{name}")
    return obj


def graph_points(graph: ROOT.TGraph) -> list[Point]:
    out: list[Point] = []
    for i in range(graph.GetN()):
        x = float(graph.GetPointX(i))
        y = float(graph.GetPointY(i))
        ex_low = float(graph.GetErrorXlow(i)) if graph.InheritsFrom("TGraphAsymmErrors") else float(graph.GetErrorX(i))
        ex_high = float(graph.GetErrorXhigh(i)) if graph.InheritsFrom("TGraphAsymmErrors") else float(graph.GetErrorX(i))
        ey_low = float(graph.GetErrorYlow(i)) if graph.InheritsFrom("TGraphAsymmErrors") else float(graph.GetErrorY(i))
        ey_high = float(graph.GetErrorYhigh(i)) if graph.InheritsFrom("TGraphAsymmErrors") else float(graph.GetErrorY(i))
        out.append(Point(x, ex_low, ex_high, y, ey_low, ey_high))
    return out


def assert_same_binning(histograms: dict[str, ROOT.TH1]) -> int:
    _, first = next(iter(histograms.items()))
    nbins = first.GetNbinsX()
    reference_edges = [float(first.GetXaxis().GetBinLowEdge(i)) for i in range(1, nbins + 1)]
    reference_edges.append(float(first.GetXaxis().GetBinUpEdge(nbins)))
    for key, obj in histograms.items():
        if obj.GetNbinsX() != nbins:
            raise RuntimeError(f"current purity histogram {key} has mismatched bin count")
        edges = [float(obj.GetXaxis().GetBinLowEdge(i)) for i in range(1, nbins + 1)]
        edges.append(float(obj.GetXaxis().GetBinUpEdge(nbins)))
        if any(abs(left - right) > 1.0e-9 for left, right in zip(edges, reference_edges)):
            raise RuntimeError(f"current purity histogram {key} has mismatched bin edges")
    return nbins


def leak_ratio(numerator: ROOT.TH1, denominator: ROOT.TH1, name: str) -> ROOT.TH1:
    """Match PPG12's TH1::Divide leakage fractions, including its errors."""
    ratio = numerator.Clone(name)
    ratio.SetDirectory(0)
    ratio.Divide(denominator)
    return ratio


def abcd_root(a: float, b: float, c: float, d: float,
              c_b: float, c_c: float, c_d: float,
              lower: float, upper: float) -> float:
    """Physical root of the PPG12 CalculatePhotonYield closure equation.

    This is the analytic quadratic form of the same equation evaluated by
    TF1::GetX in CalculatePhotonYield.C.  Selecting the root closest to the
    no-leakage ABCD solution reproduces the physical branch and is checked
    against TF1::GetX before use.
    """
    if a <= 0.0 or d == 0.0:
        return float("nan")
    no_leak = a - b * c / d
    quadratic = c_b * c_c - c_d
    linear = d + a * c_d - b * c_c - c_b * c
    constant = b * c - a * d
    if abs(quadratic) < 1.0e-20:
        roots = [] if linear == 0.0 else [-constant / linear]
    else:
        discriminant = linear * linear - 4.0 * quadratic * constant
        roots = [] if discriminant < 0.0 else [
            (-linear + math.sqrt(discriminant)) / (2.0 * quadratic),
            (-linear - math.sqrt(discriminant)) / (2.0 * quadratic),
        ]
    physical = [root for root in roots if math.isfinite(root) and lower <= root <= upper]
    return min(physical, key=lambda root: abs(root - no_leak)) if physical else float("nan")


def effective_count(histogram: ROOT.TH1, bin_index: int) -> float:
    value = float(histogram.GetBinContent(bin_index))
    sumw2 = float(histogram.GetSumw2().At(bin_index))
    if value <= 0.0 or sumw2 <= 0.0:
        raise RuntimeError(f"nonpositive effective count in {histogram.GetName()} bin {bin_index}")
    return value * value / sumw2


def ppg12_toy_estimate(rng: ROOT.TRandom3, values: tuple[float, float, float, float],
                        effective_counts: tuple[float, float, float, float],
                        leakage: tuple[float, float, float],
                        leakage_errors: tuple[float, float, float],
                        label: str) -> tuple[float, float, float, float, dict[str, float]]:
    """Reproduce PPG12's 20k effective-Poisson + leakage-Gaussian toys.

    PPG12 fills a [-1, 2] toy-purity histogram and takes a Gaussian fit over
    mean-RMS to mean+1.5RMS.  Using the same estimator avoids the unstable
    finite-difference treatment of a nonlinear leakage root.
    """
    raw_toys = ROOT.TH1D(f"h_{label}_raw", "", 1000, -1.0, 2.0)
    corr_toys = ROOT.TH1D(f"h_{label}_corr", "", 1000, -1.0, 2.0)
    raw_toys.SetDirectory(0)
    corr_toys.SetDirectory(0)
    for _ in range(20000):
        toy = tuple(value * rng.PoissonD(n_eff) / n_eff for value, n_eff in zip(values, effective_counts))
        raw_signal = abcd_root(*toy, 0.0, 0.0, 0.0, -0.5 * toy[0], 2.0 * toy[0])
        if math.isfinite(raw_signal):
            raw_toys.Fill(raw_signal / toy[0])
        c_b, c_c, c_d = (rng.Gaus(value, error) for value, error in zip(leakage, leakage_errors))
        corrected_signal = abcd_root(*toy, c_b, c_c, c_d, -0.5 * toy[0], 2.0 * toy[0])
        if math.isfinite(corrected_signal):
            corr_toys.Fill(corrected_signal / toy[0])

    def fit(histogram: ROOT.TH1, name: str) -> tuple[float, float]:
        if histogram.GetEntries() < 100.0:
            raise RuntimeError(f"too few accepted PPG12 toy throws for {name}")
        low = histogram.GetMean() - histogram.GetRMS()
        high = histogram.GetMean() + 1.5 * histogram.GetRMS()
        fit_function = ROOT.TF1(f"f_{name}", "gaus", low, high)
        histogram.Fit(fit_function, "QRMN", "", low, high)
        return float(fit_function.GetParameter(1)), float(fit_function.GetParameter(2))

    raw, raw_error = fit(raw_toys, f"{label}_raw")
    corrected, corrected_error = fit(corr_toys, f"{label}_corr")
    diagnostics = {
        "raw_toy_entries": float(raw_toys.GetEntries()),
        "raw_toy_underflow": float(raw_toys.GetBinContent(0)),
        "raw_toy_overflow": float(raw_toys.GetBinContent(raw_toys.GetNbinsX() + 1)),
        "corrected_toy_entries": float(corr_toys.GetEntries()),
        "corrected_toy_underflow": float(corr_toys.GetBinContent(0)),
        "corrected_toy_overflow": float(corr_toys.GetBinContent(corr_toys.GetNbinsX() + 1)),
    }
    return raw, raw_error, corrected, corrected_error, diagnostics


def current_points(
    inclusive_path: Path,
    photon_path: Path,
    abcd_population: str,
) -> tuple[dict[str, list[Point]], list[dict[str, object]]]:
    inclusive_file = open_root(inclusive_path)
    photon_file = open_root(photon_path)
    inclusive = inclusive_file.Get("SIM")
    photon = photon_file.Get("SIM")
    if not inclusive or not photon:
        raise RuntimeError("missing SIM directory in current inclusive or photon ROOT")

    # PPG12's isMC=true closure uses an inclusive-jet MC-as-data side.  In
    # RecoilJets the mutually exclusive signal and notmatch histograms are the
    # canonical one-for-one truth-class partition of that population.  The
    # unsuffixed histograms are retained only as a QA/display family because
    # they are not the union of the two truth classes in current outputs.
    unsuffixed_names = {
        "a": "h_tight_iso_cluster_0", "b": "h_tight_noniso_cluster_0",
        "c": "h_nontight_iso_cluster_0", "d": "h_nontight_noniso_cluster_0",
    }
    classed_names = {
        "a": ("h_tight_iso_cluster_signal_0", "h_tight_iso_cluster_notmatch_0"),
        "b": ("h_tight_noniso_cluster_signal_0", "h_tight_noniso_cluster_notmatch_0"),
        "c": ("h_nontight_iso_cluster_signal_0", "h_nontight_iso_cluster_notmatch_0"),
        "d": ("h_nontight_noniso_cluster_signal_0", "h_nontight_noniso_cluster_notmatch_0"),
    }
    leak_names = {
        "sa": "h_tight_iso_cluster_signal_0", "sb": "h_tight_noniso_cluster_signal_0",
        "sc": "h_nontight_iso_cluster_signal_0", "sd": "h_nontight_noniso_cluster_signal_0",
    }
    unsuffixed = {key: require_hist(inclusive, name) for key, name in unsuffixed_names.items()}
    classed_signal = {
        key: require_hist(inclusive, names[0]) for key, names in classed_names.items()
    }
    classed_notmatch = {
        key: require_hist(inclusive, names[1]) for key, names in classed_names.items()
    }
    classed = {}
    for key in classed_names:
        histogram = classed_signal[key].Clone(f"h_current_classed_{key}")
        histogram.SetDirectory(0)
        histogram.Add(classed_notmatch[key])
        classed[key] = histogram
    if abcd_population == "classed":
        abcd = classed
    elif abcd_population == "unsuffixed":
        abcd = unsuffixed
    else:
        raise RuntimeError(f"unknown current ABCD population: {abcd_population}")
    leak_hists = {key: require_hist(photon, name) for key, name in leak_names.items()}
    nbins = assert_same_binning({
        **{f"unsuffixed_{key}": value for key, value in unsuffixed.items()},
        **{f"classed_signal_{key}": value for key, value in classed_signal.items()},
        **{f"classed_notmatch_{key}": value for key, value in classed_notmatch.items()},
        **{f"leak_{key}": value for key, value in leak_hists.items()},
    })
    leak_b = leak_ratio(leak_hists["sb"], leak_hists["sa"], "h_current_leak_b")
    leak_c = leak_ratio(leak_hists["sc"], leak_hists["sa"], "h_current_leak_c")
    leak_d = leak_ratio(leak_hists["sd"], leak_hists["sa"], "h_current_leak_d")

    points = {spec["key"]: [] for spec in SERIES}
    table: list[dict[str, object]] = []
    rng = ROOT.TRandom3(42)
    for i in range(1, nbins + 1):
        values = tuple(float(abcd[key].GetBinContent(i)) for key in ("a", "b", "c", "d"))
        counts = tuple(effective_count(abcd[key], i) for key in ("a", "b", "c", "d"))
        leakage = tuple(float(histogram.GetBinContent(i)) for histogram in (leak_b, leak_c, leak_d))
        leakage_errors = tuple(float(histogram.GetBinError(i)) for histogram in (leak_b, leak_c, leak_d))
        raw, raw_error, corrected, corrected_error, toy_diagnostics = ppg12_toy_estimate(
            rng, values, counts, leakage, leakage_errors, f"current_purity_bin{i}")

        signal_a = float(classed_signal["a"].GetBinContent(i))
        notmatch_a = float(classed_notmatch["a"].GetBinContent(i))
        signal_error = float(classed_signal["a"].GetBinError(i))
        notmatch_error = float(classed_notmatch["a"].GetBinError(i))
        truth_total = signal_a + notmatch_a
        truth = signal_a / truth_total if truth_total > 0.0 else float("nan")
        truth_error = math.hypot(notmatch_a * signal_error, signal_a * notmatch_error) / (truth_total * truth_total) if truth_total > 0.0 else 0.0
        axis = abcd["a"].GetXaxis()
        low = float(axis.GetBinLowEdge(i))
        high = float(axis.GetBinUpEdge(i))
        x = 0.5 * (low + high)
        ex = 0.5 * (high - low)
        point_values = {"truth": (truth, truth_error), "raw": (raw, raw_error), "corrected": (corrected, corrected_error)}
        for key, (value, error) in point_values.items():
            points[key].append(Point(x, ex, ex, value, error, error))
        table.append(
            {
                "bin": i,
                "x_low_gev": low,
                "x_high_gev": high,
                "abcd_population": abcd_population,
                "abcd_a": values[0],
                "abcd_b": values[1],
                "abcd_c": values[2],
                "abcd_d": values[3],
                "abcd_a_effective_count": counts[0],
                "abcd_b_effective_count": counts[1],
                "abcd_c_effective_count": counts[2],
                "abcd_d_effective_count": counts[3],
                "unsuffixed_a": float(unsuffixed["a"].GetBinContent(i)),
                "unsuffixed_b": float(unsuffixed["b"].GetBinContent(i)),
                "unsuffixed_c": float(unsuffixed["c"].GetBinContent(i)),
                "unsuffixed_d": float(unsuffixed["d"].GetBinContent(i)),
                "classed_a": float(classed["a"].GetBinContent(i)),
                "classed_b": float(classed["b"].GetBinContent(i)),
                "classed_c": float(classed["c"].GetBinContent(i)),
                "classed_d": float(classed["d"].GetBinContent(i)),
                "inclusive_signal_a": signal_a,
                "inclusive_notmatch_a": notmatch_a,
                "inclusive_truth_classed_a": truth_total,
                "unsuffixed_minus_classed_a_qa": float(unsuffixed["a"].GetBinContent(i)) - truth_total,
                "photon_template_b_over_a": leakage[0],
                "photon_template_c_over_a": leakage[1],
                "photon_template_d_over_a": leakage[2],
                "photon_template_b_over_a_error": leakage_errors[0],
                "photon_template_c_over_a_error": leakage_errors[1],
                "photon_template_d_over_a_error": leakage_errors[2],
                "current_truth": truth,
                "current_truth_error": truth_error,
                "current_raw": raw,
                "current_raw_error": raw_error,
                "current_corrected": corrected,
                "current_corrected_error": corrected_error,
                **toy_diagnostics,
            }
        )
    inclusive_file.Close()
    photon_file.Close()
    return points, table


def reference_points(path: Path) -> dict[str, list[Point]]:
    f = open_root(path)
    out: dict[str, list[Point]] = {}
    for spec in SERIES:
        graph = f.Get(spec["reference"])
        if not graph or not graph.InheritsFrom("TGraph"):
            raise RuntimeError(f"missing graph {spec['reference']} in {path}")
        out[spec["key"]] = graph_points(graph)
    f.Close()
    return out


def match_reference(reference: list[Point], current: list[Point]) -> list[tuple[Point, Point]]:
    pairs = []
    for ref in reference:
        candidates = sorted(current, key=lambda cur: abs(cur.x - ref.x))
        if candidates and abs(candidates[0].x - ref.x) < 0.6:
            pairs.append((ref, candidates[0]))
    return pairs


def ratio_points(reference: list[Point], current: list[Point]) -> list[Point]:
    out = []
    for ref, cur in match_reference(reference, current):
        if cur.y == 0.0 or not math.isfinite(cur.y):
            continue
        value = ref.y / cur.y
        ref_low = ref.ey_low / abs(ref.y) if ref.y else 0.0
        ref_high = ref.ey_high / abs(ref.y) if ref.y else 0.0
        cur_rel = 0.5 * (cur.ey_low + cur.ey_high) / abs(cur.y)
        out.append(
            Point(
                ref.x,
                0.0,
                0.0,
                value,
                abs(value) * math.sqrt(ref_low * ref_low + cur_rel * cur_rel),
                abs(value) * math.sqrt(ref_high * ref_high + cur_rel * cur_rel),
            )
        )
    return out


def make_graph(points: list[Point], name: str, color: int, marker: int, x_shift: float = 0.0) -> ROOT.TGraphAsymmErrors:
    graph = ROOT.TGraphAsymmErrors(len(points))
    graph.SetName(name)
    graph.SetMarkerColor(color)
    graph.SetLineColor(color)
    graph.SetMarkerStyle(marker)
    graph.SetMarkerSize(0.9)
    graph.SetLineWidth(1)
    for i, point in enumerate(points):
        graph.SetPoint(i, point.x + x_shift, point.y)
        graph.SetPointError(i, point.ex_low, point.ex_high, point.ey_low, point.ey_high)
    return graph


def render(reference: dict[str, list[Point]], current: dict[str, list[Point]], output: Path) -> dict[str, list[Point]]:
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.SetTextFont(42)
    ROOT.gStyle.SetLabelFont(42, "XYZ")
    ROOT.gStyle.SetTitleFont(42, "XYZ")

    canvas = ROOT.TCanvas("c_fig3_purity_overlay", "", 720, 760)
    canvas.SetFillColor(ROOT.kWhite)
    canvas.SetFillStyle(1001)
    top = ROOT.TPad("top_fig3", "", 0.0, 0.31, 1.0, 1.0)
    bot = ROOT.TPad("bot_fig3", "", 0.0, 0.0, 1.0, 0.31)
    top.SetLeftMargin(0.14); top.SetRightMargin(0.04); top.SetTopMargin(0.06); top.SetBottomMargin(0.02)
    bot.SetLeftMargin(0.14); bot.SetRightMargin(0.04); bot.SetTopMargin(0.03); bot.SetBottomMargin(0.32)
    top.SetFillColor(ROOT.kWhite); top.SetFillStyle(1001)
    bot.SetFillColor(ROOT.kWhite); bot.SetFillStyle(1001)
    top.Draw(); bot.Draw()

    ref_graphs = {}
    cur_graphs = {}
    ratios = {}
    ratio_graphs = {}
    for spec in SERIES:
        key = spec["key"]
        ref_graphs[key] = make_graph(reference[key], f"g_ref_{key}", spec["color"], 24, -0.12)
        cur_graphs[key] = make_graph(current[key], f"g_cur_{key}", spec["color"], 20, 0.12)
        ratios[key] = ratio_points(reference[key], current[key])
        ratio_graphs[key] = make_graph(ratios[key], f"g_ratio_{key}", spec["color"], 20)

    top.cd()
    frame = ROOT.TH1F("frame_fig3_top", "", 26, 10.0, 36.0)
    frame.SetStats(False)
    frame.GetYaxis().SetRangeUser(0.0, 1.2)
    frame.GetYaxis().SetTitle("Purity")
    frame.GetYaxis().SetTitleSize(0.060); frame.GetYaxis().SetLabelSize(0.047); frame.GetYaxis().SetTitleOffset(0.93)
    frame.GetXaxis().SetLabelSize(0.0); frame.Draw("axis")
    for spec in SERIES:
        ref_graphs[spec["key"]].Draw("same p")
        cur_graphs[spec["key"]].Draw("same p")

    text = ROOT.TLatex(); text.SetNDC(True); text.SetTextFont(42); text.SetTextSize(0.035)
    text.DrawLatex(0.16, 0.87, "#bf{#it{sPHENIX}} Internal")
    text.DrawLatex(0.16, 0.82, "#it{p}+#it{p}  #sqrt{#it{s}} = 200 GeV")
    text.DrawLatex(0.16, 0.77, "PYTHIA8 inclusive MC")
    text.DrawLatex(0.16, 0.72, "bdt_nom")

    series_legend = ROOT.TLegend(0.15, 0.10, 0.60, 0.29)
    series_legend.SetBorderSize(0); series_legend.SetFillStyle(0); series_legend.SetTextFont(42); series_legend.SetTextSize(0.032)
    series_legend.SetMargin(0.13); series_legend.SetEntrySeparation(0.0)
    for spec in SERIES:
        series_legend.AddEntry(cur_graphs[spec["key"]], spec["label"], "p")
    series_legend.Draw()
    source_legend = ROOT.TLegend(0.55, 0.10, 0.94, 0.23)
    source_legend.SetBorderSize(0); source_legend.SetFillStyle(0); source_legend.SetTextFont(42); source_legend.SetTextSize(0.032)
    source_legend.AddEntry(ref_graphs["raw"], "PPG12 SDCC (open)", "p")
    source_legend.AddEntry(cur_graphs["raw"], "Current output (filled)", "p")
    source_legend.Draw(); top.RedrawAxis()

    ratio_values = [p.y for rows in ratios.values() for p in rows if math.isfinite(p.y) and p.y > 0.0]
    low = min(ratio_values + [1.0]); high = max(ratio_values + [1.0]); span = max(high - low, 0.2)
    ratio_min = max(0.0, math.floor((low - 0.10 * span) * 10.0) / 10.0)
    ratio_max = math.ceil((high + 0.10 * span) * 10.0) / 10.0
    bot.cd()
    ratio_frame = ROOT.TH1F("frame_fig3_ratio", "", 26, 10.0, 36.0)
    ratio_frame.SetStats(False); ratio_frame.GetYaxis().SetRangeUser(ratio_min, ratio_max)
    ratio_frame.GetXaxis().SetTitle("#it{E}_{T}^{#gamma,rec} [GeV]"); ratio_frame.GetYaxis().SetTitle("SDCC / Current")
    ratio_frame.GetXaxis().SetTitleSize(0.095); ratio_frame.GetYaxis().SetTitleSize(0.080)
    ratio_frame.GetXaxis().SetLabelSize(0.078); ratio_frame.GetYaxis().SetLabelSize(0.065)
    ratio_frame.GetXaxis().SetTitleOffset(1.0); ratio_frame.GetYaxis().SetTitleOffset(0.78); ratio_frame.Draw("axis")
    unity = ROOT.TLine(10.0, 1.0, 36.0, 1.0); unity.SetLineStyle(7); unity.SetLineColor(ROOT.kGray + 2); unity.Draw()
    for spec in SERIES:
        ratio_graphs[spec["key"]].Draw("same p")
    bot.RedrawAxis()
    output.parent.mkdir(parents=True, exist_ok=True)
    canvas.SaveAs(str(output))
    return ratios


def write_table(path: Path, reference: dict[str, list[Point]], current: dict[str, list[Point]], ratios: dict[str, list[Point]]) -> None:
    fields = ["series", "x_gev", "ppg12_sdcc", "ppg12_error_low", "ppg12_error_high", "current", "current_error", "sdcc_over_current", "ratio_error_low", "ratio_error_high"]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields); writer.writeheader()
        for spec in SERIES:
            key = spec["key"]
            ratio_by_x = {round(p.x, 6): p for p in ratios[key]}
            for ref, cur in match_reference(reference[key], current[key]):
                ratio = ratio_by_x.get(round(ref.x, 6))
                writer.writerow({"series": key, "x_gev": ref.x, "ppg12_sdcc": ref.y, "ppg12_error_low": ref.ey_low, "ppg12_error_high": ref.ey_high, "current": cur.y, "current_error": 0.5 * (cur.ey_low + cur.ey_high), "sdcc_over_current": ratio.y if ratio else "", "ratio_error_low": ratio.ey_low if ratio else "", "ratio_error_high": ratio.ey_high if ratio else ""})


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-root", type=Path, default=DEFAULT_REFERENCE_ROOT)
    parser.add_argument("--current-json", type=Path, default=DEFAULT_CURRENT_JSON, help="inclusive-jet current pointer")
    parser.add_argument("--photon-json", type=Path, default=DEFAULT_PHOTON_JSON, help="photon+jet current pointer for leakage templates")
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument(
        "--current-abcd-population",
        choices=("classed", "unsuffixed"),
        default="classed",
        help="classed=signal+notmatch in every ABCD region (canonical); unsuffixed=legacy QA/display family",
    )
    args = parser.parse_args()

    current_root = resolve_current_root(args.current_json)
    photon_root = resolve_current_root(args.photon_json)
    reference = reference_points(args.reference_root)
    current, current_counts = current_points(
        current_root, photon_root, args.current_abcd_population
    )
    stem = "ppg12_ian_fig3_purity_sim_sdcc_vs_current_overlay_ratio"
    if args.current_abcd_population == "unsuffixed":
        stem += "_legacy_unsuffixed_abcd"
    png = args.outdir / f"{stem}.png"
    csv_path = args.outdir / f"{stem}_points.csv"
    counts_path = args.outdir / f"{stem}_current_counts.json"
    manifest_path = args.outdir / f"{stem}_manifest.json"
    args.outdir.mkdir(parents=True, exist_ok=True)
    ratios = render(reference, current, png)
    write_table(csv_path, reference, current, ratios)
    counts_path.write_text(json.dumps(current_counts, indent=2) + "\n")
    manifest = {
        "png": str(png),
        "points_csv": str(csv_path),
        "current_counts_json": str(counts_path),
        "ppg12_reference_root": str(args.reference_root),
        "ppg12_reference_objects": {spec["key"]: spec["reference"] for spec in SERIES},
        "inclusive_artifact_pointer": str(args.current_json),
        "inclusive_root": str(current_root),
        "photon_artifact_pointer": str(args.photon_json),
        "photon_root": str(photon_root),
        "top_panel": "PPG12 SDCC open markers and current inclusive-SIM filled markers for truth, raw ABCD, and signal-leakage-corrected purity",
        "bottom_panel": "PPG12 SDCC / current output for all three purity definitions",
        "current_abcd_population": args.current_abcd_population,
        "raw_definition": (
            "PPG12 effective-Poisson toy result using mutually exclusive current "
            "inclusive (signal + notmatch) counts in every ABCD region"
            if args.current_abcd_population == "classed"
            else "legacy diagnostic using current unsuffixed inclusive A/B/C/D; not canonical because this family is not signal+notmatch"
        ),
        "corrected_definition": "PPG12 effective-Poisson A/B/C/D toys plus Gaussian cB/cC/cD throws; photon+jet signal templates define cB/cC/cD",
        "truth_definition": "A_signal / (A_signal + A_notmatch)",
        "current_mc_population_mapping": (
            "Canonical current MC closure: inclusive-jet signal+notmatch supplies "
            "the mutually exclusive MC-as-data population in A/B/C/D and photon+jet "
            "signal templates supply cB/cC/cD. Truth uses the same classed A population."
            if args.current_abcd_population == "classed"
            else "Legacy diagnostic: unsuffixed inclusive A/B/C/D supplies the MC-as-data side while truth uses signal+notmatch; populations are not one-for-one."
        ),
        "uncertainty_note": "Raw/corrected errors exactly follow the PPG12 20,000 effective-Poisson toy procedure, Gaussian leakage-fraction throws, the [-1,2] toy histogram, and the PPG12 Gaussian-fit estimator. Truth error is the independent weighted signal/notmatch ratio propagation.",
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    print(png); print(csv_path); print(counts_path); print(manifest_path)


if __name__ == "__main__":
    main()
