#!/usr/bin/env python3
"""Produce hash-bound evidence for the PPG12 stitched-purity closure gate.

Two commands are intentionally kept together because they are the two
numerical evidence products consumed by the closure workflow:

``purity-repeat``
    Load the exact 32 validated lane extracts, sum the unsuffixed inclusive
    A/B/C/D cells and photon+jet signal-leakage cells, and execute the PPG12
    truth/raw/leakage-corrected estimator twice.  Each execution starts one
    ``TRandom3(42)`` stream and uses exactly 20,000 toys per bin.  The output
    is emitted only when both canonical payload hashes are identical.

``historical-comparison``
    Read the historical final-purity graph, the upstream historical ABCD
    histograms, and the candidate inclusive/photon ROOT objects directly,
    then compute the six production-verification metrics.  ABCD and leakage
    coherent-trend significances are kept separate with per-series evidence.

The program is local and read-only with respect to physics inputs.  It never
submits jobs, merges ROOT files, changes current pointers, or promotes an
artifact.  ``purity-repeat`` requires the analysis Python environment because
the executable PPG12 uncertainty estimator uses PyROOT.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import sys
from array import array
from pathlib import Path
from typing import Any, Callable, Iterable


MODULE_DIR = Path(__file__).resolve().parent
if str(MODULE_DIR) not in sys.path:
    sys.path.insert(0, str(MODULE_DIR))

import assemble_ppg12_stitched_purity_manifest as assembler


REPO = Path(__file__).resolve().parents[3]
DEFAULT_CONTRACT = (
    REPO / "agent_context/analysis_contracts/ppg12_stitched_purity_closure.yaml"
)
PPG12_ESTIMATOR_SOURCE = REPO / "ppg12codeGit/efficiencytool/CalculatePhotonYield.C"
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
REQUIRED_HISTORICAL_SOURCE_ROLES = {
    "historical_purity",
    "historical_abcd",
    "historical_leakage",
    "candidate_purity",
    "candidate_inclusive",
    "candidate_photon",
}
_ROOT_SOLVER_DECLARED = False


class EvidenceError(RuntimeError):
    """The requested evidence cannot be produced without violating a gate."""


def _canonical_bytes(payload: Any) -> bytes:
    return json.dumps(
        payload, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode()


def _payload_sha256(payload: Any) -> str:
    return hashlib.sha256(_canonical_bytes(payload)).hexdigest()


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
    except FileNotFoundError as exc:
        raise EvidenceError(f"missing input file: {path}") from exc
    return digest.hexdigest()


def _read_json(path: Path, label: str) -> dict[str, Any]:
    try:
        payload = json.loads(path.read_text())
    except FileNotFoundError as exc:
        raise EvidenceError(f"missing {label}: {path}") from exc
    except json.JSONDecodeError as exc:
        raise EvidenceError(f"invalid JSON in {label} {path}: {exc}") from exc
    if not isinstance(payload, dict):
        raise EvidenceError(f"{label} must contain one JSON object: {path}")
    return payload


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    temporary.replace(path)


def _prepare_output(path: Path) -> None:
    if path.exists():
        path.unlink()


def _finite_number(value: Any) -> bool:
    return (
        isinstance(value, (int, float))
        and not isinstance(value, bool)
        and math.isfinite(float(value))
    )


def _finite_vector(
    value: Any,
    label: str,
    *,
    length: int | None = None,
    nonnegative: bool = False,
) -> list[float]:
    if not isinstance(value, list) or not value:
        raise EvidenceError(f"{label} must be a non-empty list")
    if not all(_finite_number(item) for item in value):
        raise EvidenceError(f"{label} contains a non-finite or non-numeric value")
    output = [float(item) for item in value]
    if length is not None and len(output) != length:
        raise EvidenceError(f"{label} has length {len(output)}; expected {length}")
    if nonnegative and any(item < 0.0 for item in output):
        raise EvidenceError(f"{label} contains a negative value")
    return output


def _load_root() -> Any:
    try:
        import ROOT
    except ModuleNotFoundError as exc:
        raise EvidenceError(
            "PyROOT is required for the canonical PPG12 toy/Gaussian-fit estimator; "
            "run with /Users/patsfan753/Desktop/analysis/env/bin/python3"
        ) from exc
    ROOT.gROOT.SetBatch(True)
    ROOT.TH1.AddDirectory(False)
    return ROOT


def _sum_observable(
    lanes: Iterable[dict[str, Any]], family: str, observable: str
) -> tuple[list[float], list[float]]:
    selected = [lane for lane in lanes if lane["family"] == family]
    if not selected:
        raise EvidenceError(f"no {family} lanes were supplied")
    first = selected[0]["observables"][observable]
    size = len(first["sumw"])
    sumw = [0.0] * size
    sumw2 = [0.0] * size
    for lane in selected:
        histogram = lane["observables"].get(observable)
        if not isinstance(histogram, dict):
            raise EvidenceError(f"{lane['lane_id']} lacks {observable}")
        values = _finite_vector(
            histogram.get("sumw"), f"{lane['lane_id']}:{observable}.sumw", length=size
        )
        variances = _finite_vector(
            histogram.get("sumw2"),
            f"{lane['lane_id']}:{observable}.sumw2",
            length=size,
            nonnegative=True,
        )
        sumw = [left + right for left, right in zip(sumw, values)]
        sumw2 = [left + right for left, right in zip(sumw2, variances)]
    return sumw, sumw2


def _aggregate_estimator_inputs(
    lanes: list[dict[str, Any]], bin_edges: list[float]
) -> dict[str, Any]:
    inclusive = {
        name: _sum_observable(lanes, "inclusive", name)
        for name in (
            "A",
            "B",
            "C",
            "D",
            "A_signal",
            "A_notmatch",
        )
    }
    photon = {
        name: _sum_observable(lanes, "photon", name)
        for name in ("A_signal", "B_signal", "C_signal", "D_signal")
    }
    return {
        "bin_edges": [float(item) for item in bin_edges],
        "inclusive": inclusive,
        "photon": photon,
    }


def _make_histogram(
    root: Any,
    name: str,
    edges: list[float],
    values: list[float],
    variances: list[float],
) -> Any:
    histogram = root.TH1D(name, "", len(edges) - 1, array("d", edges))
    histogram.SetDirectory(0)
    histogram.Sumw2()
    for index, (value, variance) in enumerate(zip(values, variances), start=1):
        histogram.SetBinContent(index, value)
        histogram.SetBinError(index, math.sqrt(variance))
    return histogram


def _declare_root_solver(root: Any) -> None:
    """Declare the executable PPG12 root function once in the ROOT process."""
    global _ROOT_SOLVER_DECLARED
    if _ROOT_SOLVER_DECLARED:
        return
    declared = root.gInterpreter.Declare(
        r"""
        double ppg12_stitched_purity_myfunc(double *x, double *params)
        {
            double NsigA = x[0];
            double NA = params[0];
            double NB = params[1];
            double NC = params[2];
            double ND = params[3];
            double cB = params[4];
            double cC = params[5];
            double cD = params[6];
            double R = params[7];
            double numerator = NB - cB * NsigA;
            double denominator = ND - cD * NsigA;
            if (denominator == 0) denominator = 1e-12;
            double fraction = (NC - cC * NsigA) / denominator;
            return NsigA - (NA - R * numerator * fraction);
        }
        """
    )
    if not declared:
        raise EvidenceError("ROOT failed to declare the executable PPG12 solver")
    _ROOT_SOLVER_DECLARED = True


def _new_root_solver(root: Any, name: str, nominal_a: float) -> Any:
    """Create the same eight-parameter TF1 used by CalculatePhotonYield.C."""
    _declare_root_solver(root)
    function = root.TF1(
        name,
        root.ppg12_stitched_purity_myfunc,
        0.0,
        nominal_a,
        8,
    )
    if not function or not function.IsValid():
        raise EvidenceError("ROOT failed to construct the executable PPG12 TF1")
    return function


def _tf1_root(
    function: Any,
    values: tuple[float, float, float, float],
    leakage: tuple[float, float, float],
) -> float:
    """Evaluate the exact PPG12 ``TF1::GetX`` root-search semantics."""
    for index, value in enumerate((*values, *leakage, 1.0)):
        function.SetParameter(index, value)
    return float(function.GetX(0.0, -0.5 * values[0], 2.0 * values[0]))


def _effective_count(value: float, variance: float, label: str) -> float:
    if value <= 0.0 or variance <= 0.0:
        raise EvidenceError(f"{label} has nonpositive sumw or Sumw2")
    output = value * value / variance
    if not math.isfinite(output) or output <= 0.0:
        raise EvidenceError(f"{label} has invalid effective count")
    return output


def _fit_toy_histogram(
    root: Any, histogram: Any, name: str
) -> tuple[float, float, int]:
    if histogram.GetEntries() < 100.0:
        raise EvidenceError(f"too few accepted PPG12 toy throws for {name}")
    low = float(histogram.GetMean() - histogram.GetRMS())
    high = float(histogram.GetMean() + 1.5 * histogram.GetRMS())
    if not math.isfinite(low) or not math.isfinite(high) or high <= low:
        raise EvidenceError(f"invalid PPG12 Gaussian-fit range for {name}")
    function = root.TF1(f"f_{name}", "gaus", low, high)
    status = int(histogram.Fit(function, "REMQN", "", low, high))
    mean = float(function.GetParameter(1))
    sigma = abs(float(function.GetParameter(2)))
    if not math.isfinite(mean) or not math.isfinite(sigma):
        raise EvidenceError(f"PPG12 Gaussian fit is non-finite for {name}")
    # CalculatePhotonYield.C consumes the fitted parameters without rejecting a
    # nonzero ROOT fit status.  Retain that status as deterministic diagnostic
    # evidence while matching the executable behavior exactly.
    return mean, sigma, status


def _draw_toy_state(
    rng: Any,
    values: tuple[float, float, float, float],
    effective_counts: tuple[float, float, float, float],
    leakage: tuple[float, float, float],
    leakage_errors: tuple[float, float, float],
) -> tuple[
    tuple[float, float, float, float],
    tuple[float, float, float],
]:
    """Consume one toy throw in the preserved PPG12 RNG order.

    ``CalculatePhotonYield.C`` always draws four Poisson deviates followed by
    the three Gaussian leakage deviates.  In particular, the Gaussian draws
    are still consumed when the Poisson-drawn region-A count is zero.  Keeping
    this in one helper makes that stream contract testable without PyROOT.
    """
    toy = tuple(
        value * rng.PoissonD(effective) / effective
        for value, effective in zip(values, effective_counts)
    )
    leakage_toy = tuple(
        rng.Gaus(value, error)
        for value, error in zip(leakage, leakage_errors)
    )
    return toy, leakage_toy


def _toy_estimate(
    root: Any,
    rng: Any,
    values: tuple[float, float, float, float],
    effective_counts: tuple[float, float, float, float],
    leakage: tuple[float, float, float],
    leakage_errors: tuple[float, float, float],
    *,
    toy_count: int,
    label: str,
) -> tuple[float, float, float, float, dict[str, float]]:
    raw_hist = root.TH1D(f"h_{label}_raw", "", 1000, -1.0, 2.0)
    corrected_hist = root.TH1D(f"h_{label}_corrected", "", 1000, -1.0, 2.0)
    raw_hist.SetDirectory(0)
    corrected_hist.SetDirectory(0)
    solver = _new_root_solver(root, f"solver_{label}", values[0])
    for _ in range(toy_count):
        toy, leakage_toy = _draw_toy_state(
            rng,
            values,
            effective_counts,
            leakage,
            leakage_errors,
        )
        if toy[0] <= 0.0:
            continue
        raw_signal = _tf1_root(solver, toy, (0.0, 0.0, 0.0))
        if math.isfinite(raw_signal):
            raw_hist.Fill(raw_signal / toy[0])
        c_b, c_c, c_d = leakage_toy
        corrected_signal = _tf1_root(solver, toy, (c_b, c_c, c_d))
        if math.isfinite(corrected_signal):
            corrected_hist.Fill(corrected_signal / toy[0])
    raw, raw_error, raw_fit_status = _fit_toy_histogram(
        root, raw_hist, f"{label}_raw"
    )
    corrected, corrected_error, corrected_fit_status = _fit_toy_histogram(
        root, corrected_hist, f"{label}_corrected"
    )
    diagnostics = {
        "raw_entries": float(raw_hist.GetEntries()),
        "raw_underflow": float(raw_hist.GetBinContent(0)),
        "raw_overflow": float(raw_hist.GetBinContent(raw_hist.GetNbinsX() + 1)),
        "corrected_entries": float(corrected_hist.GetEntries()),
        "corrected_underflow": float(corrected_hist.GetBinContent(0)),
        "corrected_overflow": float(
            corrected_hist.GetBinContent(corrected_hist.GetNbinsX() + 1)
        ),
        "raw_fit_status": raw_fit_status,
        "corrected_fit_status": corrected_fit_status,
    }
    return raw, raw_error, corrected, corrected_error, diagnostics


def _truth_series(root: Any, inputs: dict[str, Any], run_label: str) -> dict[str, list[float]]:
    edges = inputs["bin_edges"]
    signal, signal_variance = inputs["inclusive"]["A_signal"]
    denominator, denominator_variance = inputs["inclusive"]["A"]
    numerator_hist = _make_histogram(
        root,
        f"h_{run_label}_truth_signal",
        edges,
        signal,
        signal_variance,
    )
    denominator_hist = _make_histogram(
        root,
        f"h_{run_label}_truth_denominator",
        edges,
        denominator,
        denominator_variance,
    )
    graph = root.TGraphAsymmErrors(numerator_hist, denominator_hist)
    if graph.GetN() != len(edges) - 1:
        raise EvidenceError(
            "PPG12 truth-purity graph does not cover every reported bin"
        )
    values: list[float] = []
    errors: list[float] = []
    for index in range(graph.GetN()):
        expected_x = 0.5 * (edges[index] + edges[index + 1])
        observed_x = float(graph.GetPointX(index))
        if abs(observed_x - expected_x) > 1.0e-9:
            raise EvidenceError("PPG12 truth-purity graph bin order changed")
        values.append(float(graph.GetPointY(index)))
        errors.append(
            0.5
            * (
                float(graph.GetErrorYlow(index))
                + float(graph.GetErrorYhigh(index))
            )
        )
    return {"value": values, "error": errors}


def _leakage_histograms(root: Any, inputs: dict[str, Any], run_label: str) -> dict[str, Any]:
    edges = inputs["bin_edges"]
    a_values, a_variances = inputs["photon"]["A_signal"]
    denominator = _make_histogram(
        root,
        f"h_{run_label}_photon_A_signal",
        edges,
        a_values,
        a_variances,
    )
    output: dict[str, Any] = {}
    for region in ("B", "C", "D"):
        values, variances = inputs["photon"][f"{region}_signal"]
        numerator = _make_histogram(
            root,
            f"h_{run_label}_photon_{region}_signal",
            edges,
            values,
            variances,
        )
        ratio = numerator.Clone(f"h_{run_label}_c{region}")
        ratio.SetDirectory(0)
        ratio.Divide(denominator)
        output[f"c{region}"] = ratio
    return output


def _run_estimator_once(
    inputs: dict[str, Any], contract: dict[str, Any], run_label: str
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    root = _load_root()
    seed = int(contract["random_seed"])
    toy_count = int(contract["toy_count"])
    if seed != 42 or toy_count != 20000:
        raise EvidenceError("canonical estimator requires seed 42 and exactly 20000 toys")
    rng = root.TRandom3(seed)
    truth = _truth_series(root, inputs, run_label)
    leakage_histograms = _leakage_histograms(root, inputs, run_label)
    edges = inputs["bin_edges"]
    raw_values: list[float] = []
    raw_errors: list[float] = []
    corrected_values: list[float] = []
    corrected_errors: list[float] = []
    diagnostics: list[dict[str, Any]] = []
    for index in range(len(edges) - 1):
        values = tuple(
            float(inputs["inclusive"][region][0][index])
            for region in ("A", "B", "C", "D")
        )
        variances = tuple(
            float(inputs["inclusive"][region][1][index])
            for region in ("A", "B", "C", "D")
        )
        effective_counts = tuple(
            _effective_count(value, variance, f"bin {index + 1} region {region}")
            for value, variance, region in zip(
                values, variances, ("A", "B", "C", "D")
            )
        )
        leakage = tuple(
            float(leakage_histograms[f"c{region}"].GetBinContent(index + 1))
            for region in ("B", "C", "D")
        )
        leakage_errors = tuple(
            float(leakage_histograms[f"c{region}"].GetBinError(index + 1))
            for region in ("B", "C", "D")
        )
        raw, raw_error, corrected, corrected_error, toy_diagnostics = _toy_estimate(
            root,
            rng,
            values,
            effective_counts,
            leakage,
            leakage_errors,
            toy_count=toy_count,
            label=f"{run_label}_bin{index + 1}",
        )
        raw_values.append(raw)
        raw_errors.append(raw_error)
        corrected_values.append(corrected)
        corrected_errors.append(corrected_error)
        diagnostics.append(
            {
                "bin": index + 1,
                "x_low": float(edges[index]),
                "x_high": float(edges[index + 1]),
                "A": values[0],
                "B": values[1],
                "C": values[2],
                "D": values[3],
                "A_effective": effective_counts[0],
                "B_effective": effective_counts[1],
                "C_effective": effective_counts[2],
                "D_effective": effective_counts[3],
                "cB": leakage[0],
                "cC": leakage[1],
                "cD": leakage[2],
                "cB_error": leakage_errors[0],
                "cC_error": leakage_errors[1],
                "cD_error": leakage_errors[2],
                **toy_diagnostics,
            }
        )
    purity = {
        "bin_edges": [float(item) for item in edges],
        "truth": truth,
        "raw": {"value": raw_values, "error": raw_errors},
        "corrected": {"value": corrected_values, "error": corrected_errors},
    }
    for name in ("truth", "raw", "corrected"):
        if not all(
            math.isfinite(item)
            for field in ("value", "error")
            for item in purity[name][field]
        ):
            raise EvidenceError(f"{name} purity contains a non-finite result")
    return purity, diagnostics


def _build_repeated_purity_payload(
    inputs: dict[str, Any],
    contract: dict[str, Any],
    *,
    evaluator: Callable[
        [dict[str, Any], dict[str, Any], str],
        tuple[dict[str, Any], list[dict[str, Any]]],
    ] = _run_estimator_once,
) -> dict[str, Any]:
    first, first_diagnostics = evaluator(inputs, contract, "repeat1")
    repeated, repeated_diagnostics = evaluator(inputs, contract, "repeat2")
    first_sha = _payload_sha256(first)
    repeated_sha = _payload_sha256(repeated)
    if first_sha != repeated_sha or _canonical_bytes(first) != _canonical_bytes(repeated):
        raise EvidenceError(
            "fixed-seed estimator repetition was not byte-deterministic"
        )
    first_diagnostics_sha = _payload_sha256(first_diagnostics)
    repeated_diagnostics_sha = _payload_sha256(repeated_diagnostics)
    if (
        first_diagnostics_sha != repeated_diagnostics_sha
        or _canonical_bytes(first_diagnostics)
        != _canonical_bytes(repeated_diagnostics)
    ):
        raise EvidenceError(
            "fixed-seed estimator diagnostics were not byte-deterministic"
        )
    return {
        "schema": "ppg12-stitched-purity-purity/v1",
        "random_seed": int(contract["random_seed"]),
        "toy_count": int(contract["toy_count"]),
        "purity": first,
        "fixed_seed_repetition": {
            "first_output_sha256": first_sha,
            "repeated_output_sha256": repeated_sha,
        },
        "run_diagnostics": {
            "first": first_diagnostics,
            "repeated": repeated_diagnostics,
            "diagnostics_sha256": first_diagnostics_sha,
            "repeated_diagnostics_sha256": repeated_diagnostics_sha,
        },
    }


def _command_purity_repeat(args: argparse.Namespace) -> None:
    contract_path = Path(args.contract).resolve()
    contract = assembler._read_contract(contract_path)
    output = Path(args.output).resolve()
    _prepare_output(output)
    index_path = Path(args.lane_index).resolve() if args.lane_index else None
    links = assembler._load_lane_links(
        index_path=index_path,
        lane_paths=[Path(item).resolve() for item in (args.lane_json or [])],
        contract=contract,
    )
    lanes, lane_links, bin_edges = assembler._load_and_validate_lanes(
        links, contract, require_merge_input=False
    )
    assembler._validate_stitched_coverage(lanes, bin_edges)
    inputs = _aggregate_estimator_inputs(lanes, bin_edges)
    payload = _build_repeated_purity_payload(inputs, contract)
    payload["algorithm"] = {
        "name": "PPG12 CalculatePhotonYield truth/raw/leakage-corrected estimator",
        "random_stream": "one TRandom3(42) stream per full repeated evaluation",
        "toys_per_bin": 20000,
        "abcd_toys": "effective-Poisson throws using sumw^2/sumw2",
        "leakage_toys": "Gaussian cB/cC/cD throws from photon signal TH1::Divide",
        "root_solver": "executable PPG12 myfunc evaluated with ROOT TF1::GetX over [-0.5*A,2*A]",
        "toy_histogram": "1000 bins over [-1,2]",
        "fit_window": "mean-RMS to mean+1.5*RMS",
        "fit_model": "ROOT Gaussian with REMQN options",
        "truth": "inclusive A_signal divided by unsuffixed inclusive A via TGraphAsymmErrors",
        "source": {
            "path": str(PPG12_ESTIMATOR_SOURCE.resolve()),
            "sha256": _file_sha256(PPG12_ESTIMATOR_SOURCE.resolve()),
        },
        "producer": {
            "path": str(Path(__file__).resolve()),
            "sha256": _file_sha256(Path(__file__).resolve()),
        },
    }
    payload["inputs"] = {
        "contract": {
            "path": str(contract_path),
            "sha256": assembler._payload_sha256(contract),
        },
        "lane_index": (
            {"path": str(index_path), "sha256": _file_sha256(index_path)}
            if index_path is not None
            else None
        ),
        "lanes": lane_links,
        "lane_set_sha256": _payload_sha256(lane_links),
    }
    _write_json(output, payload)


def _open_root_file(root: Any, path: Path, label: str) -> Any:
    handle = root.TFile.Open(str(path), "READ")
    if (
        not handle
        or handle.IsZombie()
        or handle.TestBit(root.TFile.kRecovered)
    ):
        raise EvidenceError(f"invalid or recovered {label} ROOT file: {path}")
    return handle


def _require_root_histogram(handle: Any, object_name: str, label: str) -> Any:
    histogram = handle.Get(object_name)
    if not histogram or not histogram.InheritsFrom("TH1"):
        raise EvidenceError(f"missing {label} histogram: {object_name}")
    clone = histogram.Clone(
        f"evidence_{hashlib.sha256((label + object_name).encode()).hexdigest()[:16]}"
    )
    clone.SetDirectory(0)
    return clone


def _histogram_edges(histogram: Any) -> list[float]:
    axis = histogram.GetXaxis()
    count = histogram.GetNbinsX()
    return [float(axis.GetBinLowEdge(index)) for index in range(1, count + 1)] + [
        float(axis.GetBinUpEdge(count))
    ]


def _require_edges(observed: list[float], expected: list[float], label: str) -> None:
    if len(observed) != len(expected) or any(
        abs(left - right) > 1.0e-9 for left, right in zip(observed, expected)
    ):
        raise EvidenceError(f"historical/candidate binning mismatch: {label}")


def _histogram_series(histogram: Any, edges: list[float], label: str) -> dict[str, Any]:
    _require_edges(_histogram_edges(histogram), edges, label)
    return {
        "value": [
            float(histogram.GetBinContent(index))
            for index in range(1, histogram.GetNbinsX() + 1)
        ],
        "error": [
            abs(float(histogram.GetBinError(index)))
            for index in range(1, histogram.GetNbinsX() + 1)
        ],
    }


def _graph_series(handle: Any, name: str, edges: list[float]) -> dict[str, Any]:
    graph = handle.Get(name)
    if not graph or not graph.InheritsFrom("TGraph"):
        raise EvidenceError(f"missing historical purity graph: {name}")
    points: list[tuple[float, float, float]] = []
    for index in range(graph.GetN()):
        x = float(graph.GetPointX(index))
        value = float(graph.GetPointY(index))
        if graph.InheritsFrom("TGraphAsymmErrors"):
            error = 0.5 * (
                abs(float(graph.GetErrorYlow(index)))
                + abs(float(graph.GetErrorYhigh(index)))
            )
        else:
            error = abs(float(graph.GetErrorY(index)))
        points.append((x, value, error))
    values: list[float] = []
    errors: list[float] = []
    used: set[int] = set()
    for bin_index, (low, high) in enumerate(zip(edges, edges[1:]), start=1):
        center = 0.5 * (low + high)
        matches = [
            (index, row)
            for index, row in enumerate(points)
            if abs(row[0] - center) <= max(1.0e-9, 1.0e-6 * (high - low))
        ]
        if len(matches) != 1:
            raise EvidenceError(
                f"historical purity graph has {len(matches)} points for bin {bin_index}"
            )
        point_index, (_, value, error) = matches[0]
        if point_index in used:
            raise EvidenceError("historical purity graph point was matched twice")
        used.add(point_index)
        values.append(value)
        errors.append(error)
    if len(used) != graph.GetN():
        raise EvidenceError(
            "historical purity graph contains points outside candidate binning"
        )
    return {"value": values, "error": errors}


def _load_candidate_corrected_purity(path: Path) -> tuple[list[float], dict[str, Any]]:
    payload = _read_json(path, "candidate purity evidence")
    if payload.get("schema") != "ppg12-stitched-purity-purity/v1":
        raise EvidenceError("candidate purity evidence uses an unsupported schema")
    purity = payload.get("purity")
    if not isinstance(purity, dict):
        raise EvidenceError("candidate purity evidence lacks purity")
    edges = _finite_vector(purity.get("bin_edges"), "candidate purity bin_edges")
    if len(edges) < 2 or any(right <= left for left, right in zip(edges, edges[1:])):
        raise EvidenceError("candidate purity bin_edges must be strictly increasing")
    corrected = purity.get("corrected")
    if not isinstance(corrected, dict):
        raise EvidenceError("candidate purity evidence lacks corrected series")
    series = {
        "value": _finite_vector(
            corrected.get("value"),
            "candidate corrected purity value",
            length=len(edges) - 1,
        ),
        "error": _finite_vector(
            corrected.get("error"),
            "candidate corrected purity error",
            length=len(edges) - 1,
            nonnegative=True,
        ),
    }
    canonical_purity = {
        "bin_edges": edges,
        **{
            name: purity.get(name)
            for name in ("truth", "raw", "corrected")
        },
    }
    repetition = payload.get("fixed_seed_repetition")
    canonical_sha = _payload_sha256(canonical_purity)
    if (
        payload.get("random_seed") != 42
        or payload.get("toy_count") != 20000
        or not isinstance(repetition, dict)
        or repetition.get("first_output_sha256") != canonical_sha
        or repetition.get("repeated_output_sha256") != canonical_sha
    ):
        raise EvidenceError(
            "candidate purity is not bound to a deterministic seed-42/20000-toy repetition"
        )
    return edges, series


def _historical_input_from_root(args: argparse.Namespace) -> tuple[dict[str, Any], Path]:
    root = _load_root()
    candidate_purity_path = Path(args.candidate_purity).resolve()
    historical_final_path = Path(args.historical_final_root).resolve()
    historical_abcd_path = Path(args.historical_abcd_root).resolve()
    candidate_inclusive_path = Path(args.candidate_inclusive_root).resolve()
    candidate_photon_path = Path(args.candidate_photon_root).resolve()
    edges, candidate_corrected = _load_candidate_corrected_purity(
        candidate_purity_path
    )
    historical_final = _open_root_file(
        root, historical_final_path, "historical final-purity"
    )
    historical_abcd = _open_root_file(
        root, historical_abcd_path, "historical ABCD"
    )
    candidate_inclusive = _open_root_file(
        root, candidate_inclusive_path, "candidate inclusive"
    )
    candidate_photon = _open_root_file(
        root, candidate_photon_path, "candidate photon"
    )
    try:
        final_purity = {
            "reference": _graph_series(historical_final, "gpurity_leak", edges),
            "candidate": candidate_corrected,
        }
        abcd_names = {
            "A": "h_tight_iso_cluster_0",
            "B": "h_tight_noniso_cluster_0",
            "C": "h_nontight_iso_cluster_0",
            "D": "h_nontight_noniso_cluster_0",
        }
        abcd: dict[str, Any] = {}
        for region, name in abcd_names.items():
            reference = _require_root_histogram(
                historical_abcd, name, f"historical ABCD {region}"
            )
            candidate = _require_root_histogram(
                candidate_inclusive, f"SIM/{name}", f"candidate ABCD {region}"
            )
            abcd[region] = {
                "reference": _histogram_series(
                    reference, edges, f"historical ABCD {region}"
                ),
                "candidate": _histogram_series(
                    candidate, edges, f"candidate ABCD {region}"
                ),
            }
        candidate_signal_a = _require_root_histogram(
            candidate_photon,
            "SIM/h_tight_iso_cluster_signal_0",
            "candidate leakage denominator",
        )
        leakage: dict[str, Any] = {}
        leakage_names = {
            "cB": "h_tight_noniso_cluster_signal_0",
            "cC": "h_nontight_iso_cluster_signal_0",
            "cD": "h_nontight_noniso_cluster_signal_0",
        }
        for ratio_name, candidate_name in leakage_names.items():
            reference = _require_root_histogram(
                historical_final,
                f"h_leak_{ratio_name[-1]}",
                f"historical leakage {ratio_name}",
            )
            numerator = _require_root_histogram(
                candidate_photon,
                f"SIM/{candidate_name}",
                f"candidate leakage {ratio_name}",
            )
            ratio = numerator.Clone(f"candidate_{ratio_name}")
            ratio.SetDirectory(0)
            ratio.Divide(candidate_signal_a)
            leakage[ratio_name] = {
                "reference": _histogram_series(
                    reference, edges, f"historical leakage {ratio_name}"
                ),
                "candidate": _histogram_series(
                    ratio, edges, f"candidate leakage {ratio_name}"
                ),
            }
    finally:
        historical_final.Close()
        historical_abcd.Close()
        candidate_inclusive.Close()
        candidate_photon.Close()

    source_rows = (
        ("historical_purity", historical_final_path),
        ("historical_abcd", historical_abcd_path),
        ("historical_leakage", historical_final_path),
        ("candidate_purity", candidate_purity_path),
        ("candidate_inclusive", candidate_inclusive_path),
        ("candidate_photon", candidate_photon_path),
    )
    return (
        {
            "schema": "ppg12-stitched-purity-historical-input/v1",
            "final_purity_series": "corrected",
            "bin_edges": edges,
            "final_purity": final_purity,
            "abcd": abcd,
            "leakage": leakage,
            "source_links": [
                {
                    "role": role,
                    "path": str(path),
                    "sha256": _file_sha256(path),
                }
                for role, path in source_rows
            ],
        },
        candidate_purity_path,
    )


def _normalize_source_links(
    payload: dict[str, Any], parent: Path
) -> list[dict[str, str]]:
    raw_links = payload.get("source_links")
    if not isinstance(raw_links, list):
        raise EvidenceError("historical input source_links must be a list")
    observed: dict[str, dict[str, str]] = {}
    for index, link in enumerate(raw_links):
        if not isinstance(link, dict):
            raise EvidenceError(f"source_links[{index}] must be an object")
        role = link.get("role")
        raw_path = link.get("path")
        expected_sha = link.get("sha256")
        if not isinstance(role, str) or not role:
            raise EvidenceError(f"source_links[{index}].role is required")
        if role in observed:
            raise EvidenceError(f"duplicate historical source role: {role}")
        if not isinstance(raw_path, str) or not raw_path:
            raise EvidenceError(f"source_links[{index}].path is required")
        if not isinstance(expected_sha, str) or SHA256_RE.fullmatch(expected_sha) is None:
            raise EvidenceError(f"source_links[{index}].sha256 is invalid")
        path = Path(raw_path).expanduser()
        if not path.is_absolute():
            path = (parent / path).resolve()
        observed_sha = _file_sha256(path)
        if observed_sha != expected_sha:
            raise EvidenceError(
                f"historical source hash mismatch for {role}: "
                f"expected {expected_sha}, observed {observed_sha}"
            )
        observed[role] = {
            "role": role,
            "path": str(path),
            "sha256": observed_sha,
        }
    missing = sorted(REQUIRED_HISTORICAL_SOURCE_ROLES - set(observed))
    if missing:
        raise EvidenceError(f"historical source_links are missing roles: {missing}")
    return [observed[role] for role in sorted(observed)]


def _comparison_vectors(
    payload: Any, label: str, bin_count: int
) -> tuple[list[float], list[float], list[float], list[float]]:
    if not isinstance(payload, dict):
        raise EvidenceError(f"{label} must be an object")
    reference = payload.get("reference")
    candidate = payload.get("candidate")
    if not isinstance(reference, dict) or not isinstance(candidate, dict):
        raise EvidenceError(f"{label} requires reference and candidate objects")
    return (
        _finite_vector(reference.get("value"), f"{label}.reference.value", length=bin_count),
        _finite_vector(
            reference.get("error"),
            f"{label}.reference.error",
            length=bin_count,
            nonnegative=True,
        ),
        _finite_vector(candidate.get("value"), f"{label}.candidate.value", length=bin_count),
        _finite_vector(
            candidate.get("error"),
            f"{label}.candidate.error",
            length=bin_count,
            nonnegative=True,
        ),
    )


def _series_metrics(
    payload: Any, label: str, bin_edges: list[float]
) -> dict[str, Any]:
    bin_count = len(bin_edges) - 1
    reference, reference_error, candidate, candidate_error = _comparison_vectors(
        payload, label, bin_count
    )
    points: list[dict[str, float | int]] = []
    chi2 = 0.0
    maximum_pull = 0.0
    weighted_ratio_sum = 0.0
    ratio_weight_sum = 0.0
    for index, (ref, ref_error, cand, cand_error) in enumerate(
        zip(reference, reference_error, candidate, candidate_error), start=1
    ):
        if ref == 0.0:
            raise EvidenceError(f"{label} reference is zero in bin {index}")
        variance = ref_error * ref_error + cand_error * cand_error
        if variance <= 0.0 or not math.isfinite(variance):
            raise EvidenceError(f"{label} has no finite uncertainty in bin {index}")
        pull = (cand - ref) / math.sqrt(variance)
        ratio = cand / ref
        ratio_variance = (cand_error / ref) ** 2 + (
            cand * ref_error / (ref * ref)
        ) ** 2
        if ratio_variance <= 0.0 or not math.isfinite(ratio_variance):
            raise EvidenceError(f"{label} has invalid ratio uncertainty in bin {index}")
        ratio_error = math.sqrt(ratio_variance)
        ratio_weight = 1.0 / ratio_variance
        weighted_ratio_sum += ratio_weight * ratio
        ratio_weight_sum += ratio_weight
        chi2 += pull * pull
        maximum_pull = max(maximum_pull, abs(pull))
        points.append(
            {
                "bin": index,
                "x_low": float(bin_edges[index - 1]),
                "x_high": float(bin_edges[index]),
                "reference": ref,
                "reference_error": ref_error,
                "candidate": cand,
                "candidate_error": cand_error,
                "candidate_over_reference": ratio,
                "ratio_error": ratio_error,
                "pull": pull,
            }
        )
    weighted_mean = weighted_ratio_sum / ratio_weight_sum
    weighted_error = 1.0 / math.sqrt(ratio_weight_sum)
    return {
        "point_count": bin_count,
        "chi2": chi2,
        "ndf": bin_count,
        "chi2_ndf": chi2 / bin_count,
        "max_abs_pull": maximum_pull,
        "weighted_mean_ratio": weighted_mean,
        "weighted_mean_ratio_error": weighted_error,
        "coherent_trend_sigma": abs(weighted_mean - 1.0) / weighted_error,
        "points": points,
    }


def _historical_metrics_payload(
    source: dict[str, Any], source_path: Path, contract: dict[str, Any]
) -> dict[str, Any]:
    if source.get("schema") != "ppg12-stitched-purity-historical-input/v1":
        raise EvidenceError("unsupported historical comparison input schema")
    if source.get("final_purity_series") != "corrected":
        raise EvidenceError(
            "historical comparison must explicitly identify final_purity_series=corrected"
        )
    edges = _finite_vector(source.get("bin_edges"), "historical bin_edges")
    if len(edges) < 2 or any(right <= left for left, right in zip(edges, edges[1:])):
        raise EvidenceError("historical bin_edges must be strictly increasing")
    final_metrics = _series_metrics(source.get("final_purity"), "final_purity", edges)
    abcd = source.get("abcd")
    leakage = source.get("leakage")
    if not isinstance(abcd, dict) or set(abcd) != {"A", "B", "C", "D"}:
        raise EvidenceError("historical ABCD coverage must be exactly A, B, C, D")
    if not isinstance(leakage, dict) or set(leakage) != {"cB", "cC", "cD"}:
        raise EvidenceError("historical leakage coverage must be exactly cB, cC, cD")
    abcd_metrics = {
        name: _series_metrics(abcd[name], f"abcd.{name}", edges)
        for name in ("A", "B", "C", "D")
    }
    leakage_metrics = {
        name: _series_metrics(leakage[name], f"leakage.{name}", edges)
        for name in ("cB", "cC", "cD")
    }
    source_links = _normalize_source_links(source, source_path.parent)
    abcd_trend = max(row["coherent_trend_sigma"] for row in abcd_metrics.values())
    leakage_trend = max(
        row["coherent_trend_sigma"] for row in leakage_metrics.values()
    )
    payload = {
        "schema": "ppg12-stitched-purity-historical-comparison/v1",
        "final_purity_series": "corrected",
        "chi2_ndf": final_metrics["chi2_ndf"],
        "max_abs_pull": final_metrics["max_abs_pull"],
        "weighted_mean_ratio": final_metrics["weighted_mean_ratio"],
        "weighted_mean_ratio_error": final_metrics["weighted_mean_ratio_error"],
        "abcd_coherent_trend_max_sigma": abcd_trend,
        "leakage_coherent_trend_max_sigma": leakage_trend,
        "metric_definition": {
            "ratio": "candidate/reference with independent propagated errors",
            "pull": "(candidate-reference)/sqrt(candidate_error^2+reference_error^2)",
            "chi2_ndf": "sum(pull^2)/N for the final corrected-purity points; no fitted parameters",
            "weighted_mean_ratio": "inverse-ratio-variance weighted mean over final corrected-purity bins",
            "coherent_trend_sigma": "absolute weighted-mean ratio displacement from unity divided by its error; maximum reported separately over A-D and cB-cD",
        },
        "final_purity": final_metrics,
        "abcd": abcd_metrics,
        "leakage": leakage_metrics,
        "source_links": source_links,
        "source_set_sha256": _payload_sha256(source_links),
        "input": {
            "path": str(source_path),
            "sha256": _file_sha256(source_path),
        },
    }
    tolerances = contract["tolerances"]
    unity_sigma = float(
        contract["production_historical_requirements"]
        ["weighted_mean_ratio_compatible_with_unity_sigma"]
    )
    payload["threshold_evaluation"] = {
        "status": (
            "PASS"
            if payload["chi2_ndf"] <= float(tolerances["historical_chi2_ndf_max"])
            and payload["max_abs_pull"] <= float(tolerances["historical_pull_max_abs"])
            and abs(payload["weighted_mean_ratio"] - 1.0)
            <= unity_sigma * payload["weighted_mean_ratio_error"]
            and abcd_trend <= float(tolerances["historical_coherent_trend_max_sigma"])
            and leakage_trend
            <= float(tolerances["historical_coherent_trend_max_sigma"])
            else "FAIL"
        ),
        "contract_sha256": assembler._payload_sha256(contract),
    }
    return payload


def _command_historical_comparison(args: argparse.Namespace) -> None:
    contract = assembler._read_contract(Path(args.contract).resolve())
    output = Path(args.output).resolve()
    _prepare_output(output)
    source, source_path = _historical_input_from_root(args)
    payload = _historical_metrics_payload(source, source_path, contract)
    payload["input"]["mode"] = "direct_root_objects"
    payload["root_objects"] = {
        "historical_final_purity": "gpurity_leak",
        "historical_abcd": {
            "A": "h_tight_iso_cluster_0",
            "B": "h_tight_noniso_cluster_0",
            "C": "h_nontight_iso_cluster_0",
            "D": "h_nontight_noniso_cluster_0",
        },
        "historical_leakage": ["h_leak_B", "h_leak_C", "h_leak_D"],
        "candidate_abcd_namespace": "SIM",
        "candidate_leakage_namespace": "SIM",
    }
    payload["producer"] = {
        "path": str(Path(__file__).resolve()),
        "sha256": _file_sha256(Path(__file__).resolve()),
    }
    _write_json(output, payload)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--contract", default=str(DEFAULT_CONTRACT), help="JSON-compatible YAML contract"
    )
    commands = parser.add_subparsers(dest="command", required=True)

    repeat = commands.add_parser("purity-repeat")
    repeat_inputs = repeat.add_mutually_exclusive_group(required=True)
    repeat_inputs.add_argument("--lane-index")
    repeat_inputs.add_argument("--lane-json", action="append")
    repeat.add_argument("--output", required=True)

    historical = commands.add_parser("historical-comparison")
    historical.add_argument("--candidate-purity", required=True)
    historical.add_argument("--historical-final-root", required=True)
    historical.add_argument("--historical-abcd-root", required=True)
    historical.add_argument("--candidate-inclusive-root", required=True)
    historical.add_argument("--candidate-photon-root", required=True)
    historical.add_argument("--output", required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        if args.command == "purity-repeat":
            _command_purity_repeat(args)
        else:
            _command_historical_comparison(args)
    except (EvidenceError, assembler.AssemblyError, OSError, ValueError) as exc:
        output = Path(args.output).resolve()
        _prepare_output(output)
        print(json.dumps({"status": "FAIL", "error": str(exc)}, sort_keys=True))
        return 2
    output = Path(args.output).resolve()
    print(
        json.dumps(
            {"status": "PASS", "output": str(output), "sha256": _file_sha256(output)},
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
