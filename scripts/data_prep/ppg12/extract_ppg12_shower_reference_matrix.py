#!/usr/bin/env python3
"""Extract the complete PPG12 shower-shape reference matrix from ROOT.

The historical PPG12 plotting macro rebins the shower-shape axis, restricts
the visible range, projects over the full isolation axis (including its ROOT
underflow and overflow bins), and normalizes the visible projection to unit
area.  This helper reproduces that contract in a compact, renderer-friendly
JSON without modifying the source ROOT files.

The three input ROOT files are read-only recovery inputs copied from the
PPG12 SDCC results area.  Their expected SHA-256 digests are intentionally
locked below so a same-named but different file cannot silently become the
reference.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import uproot


DEFAULT_INPUT_DIR = Path("/private/tmp/ppg12_shower_source_20260716")
DEFAULT_OUTPUT_DIR = Path(
    "dataOutput/ppg12PhotonYield/ppg12_shower_reference_matrix_20260716"
)

CUT_REPAIR_COMMIT = {
    "repository": "ppg12codeGit",
    "commit": "29f8223bd9b36dffab07961b597afa94185bbdf1",
    "committed_at": "2026-05-14T13:02:28-04:00",
    "subject": "BDT cut fix (ET=42 crossover) + silent-merge recovery + audit/oneforall patches",
    "interpretation": (
        "The nominal canonical BDT-cut definitions entered Git at this commit. "
        "A ROOT-file modification time after this boundary is supporting timing "
        "evidence, but is not by itself an embedded configuration manifest."
    ),
}

SHOWER_SHAPE_DIAGNOSTIC_BOUNDARY = {
    "source_config": "ppg12codeGit/efficiencytool/config_showershape.yaml",
    "archived_tight_threshold": "score > 0.86 - 0.02 E_T",
    "nominal_threshold_fixed_2026_05_14": (
        "score > 0.815625 - 0.0015625 E_T"
    ),
    "interpretation": (
        "The archived shower-shape cut2 family uses the diagnostic threshold "
        "stored in config_showershape.yaml. Git history shows that threshold "
        "was not replaced by the May-14 nominal-cut repair."
    ),
}

SOURCES: dict[str, dict[str, str]] = {
    "signal": {
        "filename": "MC_efficiencyshower_shape_signal_combined_showershape.root",
        "remote_path": (
            "/sphenix/user/shuhangli/ppg12/efficiencytool/results/"
            "MC_efficiencyshower_shape_signal_combined_showershape.root"
        ),
        "remote_mtime": "2026-06-03T17:33:15-04:00",
        "sha256": "cd5b04a20b5588519ae9b65a3f1cc80c2fd74739b747ab806793320cb9fdbce8",
        "cut_boundary_relation": "postdates_2026_05_14_cut_repair_by_mtime",
    },
    "inclusive": {
        "filename": "MC_efficiencyshower_shape_jet_inclusive_combined_showershape.root",
        "remote_path": (
            "/sphenix/user/shuhangli/ppg12/efficiencytool/results/"
            "MC_efficiencyshower_shape_jet_inclusive_combined_showershape.root"
        ),
        "remote_mtime": "2026-06-03T17:33:25-04:00",
        "sha256": "b5115e793538c717e4efed31252e2db8def19520e84435c57241a5d6690bba84",
        "cut_boundary_relation": "postdates_2026_05_14_cut_repair_by_mtime",
    },
    "data": {
        "filename": "data_histoshower_shape_showershape.root",
        "remote_path": (
            "/sphenix/user/shuhangli/ppg12/efficiencytool/results/"
            "data_histoshower_shape_showershape.root"
        ),
        "remote_mtime": "2026-04-21T12:23:39.406100393-04:00",
        "sha256": "02dcd31009094915e7265b4902ee309dd741a5ec72e747ae5565fdc6881a6e7e",
        "cut_boundary_relation": "predates_2026_05_14_cut_repair",
    },
}

VARIABLES: dict[str, dict[str, Any]] = {
    "weta_cogx": {
        "root_base": "h2d_weta_cogx",
        "axis_label": "w_eta^COGX",
        "xmin": 0.0,
        "xmax": 2.0,
        "rebin": 4,
    },
    "wphi_cogx": {
        "root_base": "h2d_wphi_cogx",
        "axis_label": "w_phi^COGX",
        "xmin": 0.0,
        "xmax": 2.0,
        "rebin": 4,
    },
    "e11_to_e33": {
        "root_base": "h2d_e11_to_e33",
        "axis_label": "E_1x1 / E_3x3",
        "xmin": 0.0,
        "xmax": 1.0,
        "rebin": 4,
    },
    "e32_to_e35": {
        "root_base": "h2d_e32_to_e35",
        "axis_label": "E_3x2 / E_3x5",
        "xmin": 0.4,
        "xmax": 1.0,
        "rebin": 1,
    },
    "et1": {
        "root_base": "h2d_et1",
        "axis_label": "et1",
        "xmin": 0.3,
        "xmax": 1.0,
        "rebin": 1,
    },
    "et2": {
        "root_base": "h2d_et2",
        "axis_label": "et2",
        "xmin": 0.0,
        "xmax": 1.0,
        "rebin": 4,
    },
    "et3": {
        "root_base": "h2d_et3",
        "axis_label": "et3",
        "xmin": 0.0,
        "xmax": 1.0,
        "rebin": 4,
    },
    "et4": {
        "root_base": "h2d_et4",
        "axis_label": "et4",
        "xmin": 0.0,
        "xmax": 0.3,
        "rebin": 1,
    },
}

STAGES: dict[str, dict[str, Any]] = {
    "before_ncb": {
        "pt_index": 3,
        "cut_index": 0,
        "pt_gev": [22.0, 28.0],
        "label": "Before NCB preselection",
        "cut_lineage": (
            "cut0 precedes the common NCB and shower-quality preselection; it "
            "does not evaluate the tight BDT threshold."
        ),
    },
    "after_ncb": {
        "pt_index": 2,
        "cut_index": 1,
        "pt_gev": [18.0, 22.0],
        "label": "Complete PPG12 preselection",
        "cut_lineage": (
            "cut1 follows the complete common preselection, including the NCB "
            "requirement and loose shower-quality requirements; it is not an "
            "NCB-only selection and does not evaluate the tight BDT threshold."
        ),
    },
    "tight_id": {
        "pt_index": 0,
        "cut_index": 2,
        "pt_gev": [10.0, 14.0],
        "label": "Archived tight-stage diagnostic",
        "cut_lineage": (
            "cut2 uses the archived shower-shape diagnostic threshold score > "
            "0.86 - 0.02 E_T, not the nominal score > 0.815625 - 0.0015625 E_T "
            "definition fixed on 2026-05-14. The panel is historical selection-"
            "shape context and is not a canonical tight-ID closure test."
        ),
    },
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def regroup(values: np.ndarray, factor: int) -> np.ndarray:
    """Apply ROOT-compatible fixed-factor rebinning to normal-axis bins."""

    if factor == 1:
        return values.astype(np.float64, copy=True)
    if values.size % factor:
        raise ValueError(
            f"cannot rebin {values.size} bins by factor {factor} without remainder"
        )
    return values.reshape(values.size // factor, factor).sum(axis=1, dtype=np.float64)


def curve_from_histogram(
    root_file: uproot.ReadOnlyFile,
    histogram: str,
    *,
    xmin: float,
    xmax: float,
    rebin: int,
) -> dict[str, Any]:
    """Project one TH2 using the historical PPG12 display contract."""

    if histogram not in root_file:
        raise KeyError(f"missing ROOT object: {histogram}")
    obj = root_file[histogram]
    if not obj.classname.startswith("TH2"):
        raise TypeError(f"{histogram} is {obj.classname}, expected TH2")

    # ROOT's default ProjectionX includes the Y underflow and overflow bins.
    values_flow = np.asarray(obj.values(flow=True), dtype=np.float64)
    variances_flow = obj.variances(flow=True)
    if variances_flow is None:
        raise ValueError(f"{histogram} has no Sumw2/variance payload")
    variances_flow = np.asarray(variances_flow, dtype=np.float64)

    # Remove X under/overflow from the returned TH1 normal bins, but retain all
    # Y flow bins in the projection, matching ProjectionX(name) defaults.
    projected = values_flow[1:-1, :].sum(axis=1, dtype=np.float64)
    projected_var = variances_flow[1:-1, :].sum(axis=1, dtype=np.float64)
    edges = np.asarray(obj.axis(0).edges(flow=False), dtype=np.float64)

    projected = regroup(projected, rebin)
    projected_var = regroup(projected_var, rebin)
    rebinned_edges = edges[::rebin]
    if rebinned_edges.size != projected.size + 1:
        raise ValueError(f"edge/count mismatch after rebinning {histogram}")

    centers = 0.5 * (rebinned_edges[:-1] + rebinned_edges[1:])
    tolerance = 1.0e-10
    visible = (centers >= xmin - tolerance) & (centers <= xmax + tolerance)
    if not np.any(visible):
        raise ValueError(f"display range [{xmin}, {xmax}] selects no bins in {histogram}")

    visible_values = projected[visible]
    visible_variances = projected_var[visible]
    normalization = float(np.sum(visible_values, dtype=np.float64))
    if not math.isfinite(normalization) or normalization <= 0.0:
        raise ValueError(f"non-positive visible integral for {histogram}: {normalization}")

    normalized = visible_values / normalization
    errors = np.sqrt(np.clip(visible_variances, 0.0, None)) / normalization
    visible_edges = np.concatenate(
        [rebinned_edges[:-1][visible], [rebinned_edges[1:][visible][-1]]]
    )

    return {
        "edges": visible_edges.tolist(),
        "centers": centers[visible].tolist(),
        "values": normalized.tolist(),
        "errors": errors.tolist(),
        "raw_integral_all_bins": float(np.sum(projected, dtype=np.float64)),
        "raw_integral_display": normalization,
        "display_integral_after_unit_normalization": float(
            np.sum(normalized, dtype=np.float64)
        ),
    }


def ncb_tail_scale(data_file: uproot.ReadOnlyFile, pt_index: int) -> dict[str, Any]:
    """Recover the PPG12 cut4 normalization from the w_eta > 1.4 data tail."""

    curves: dict[int, tuple[np.ndarray, np.ndarray]] = {}
    for cut_index in (0, 4):
        name = f"h2d_weta_cogx_eta0_pt{pt_index}_cut{cut_index}"
        obj = data_file[name]
        values = np.asarray(obj.values(flow=True), dtype=np.float64)
        projected = regroup(values[1:-1, :].sum(axis=1), 4)
        edges = np.asarray(obj.axis(0).edges(flow=False), dtype=np.float64)[::4]
        centers = 0.5 * (edges[:-1] + edges[1:])
        curves[cut_index] = (centers, projected)

    centers_all, values_all = curves[0]
    centers_ncb, values_ncb = curves[4]
    if not np.array_equal(centers_all, centers_ncb):
        raise ValueError("cut0/cut4 weta axes differ")
    tail = centers_all >= 1.400001
    all_total = float(np.sum(values_all))
    ncb_total = float(np.sum(values_ncb))
    all_tail = float(np.sum(values_all[tail]))
    ncb_tail = float(np.sum(values_ncb[tail]))
    if min(all_total, ncb_total, ncb_tail) <= 0.0:
        raise ValueError("invalid NCB tail normalization integrals")
    fraction_all = all_tail / all_total
    fraction_ncb = ncb_tail / ncb_total
    scale = fraction_all / fraction_ncb
    return {
        "definition": "(cut0 data fraction with w_eta > 1.4) / (cut4 NCB-tagged data fraction with w_eta > 1.4)",
        "raw_histograms": [
            f"h2d_weta_cogx_eta0_pt{pt_index}_cut0",
            f"h2d_weta_cogx_eta0_pt{pt_index}_cut4",
        ],
        "all_total": all_total,
        "all_tail": all_tail,
        "ncb_total": ncb_total,
        "ncb_tail": ncb_tail,
        "fraction_all": fraction_all,
        "fraction_ncb": fraction_ncb,
        "scale": scale,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    input_dir = args.input_dir.resolve()
    output_dir = args.output_dir.resolve()

    source_records: dict[str, dict[str, Any]] = {}
    source_paths: dict[str, Path] = {}
    for lane, metadata in SOURCES.items():
        path = input_dir / metadata["filename"]
        if not path.is_file():
            raise FileNotFoundError(path)
        actual_hash = sha256(path)
        if actual_hash != metadata["sha256"]:
            raise ValueError(
                f"SHA-256 mismatch for {path}: expected {metadata['sha256']}, got {actual_hash}"
            )
        source_paths[lane] = path
        source_records[lane] = {
            **metadata,
            "local_recovery_input": str(path),
            "size_bytes": path.stat().st_size,
            "sha256_verified": True,
        }

    roots = {lane: uproot.open(path) for lane, path in source_paths.items()}
    variables_payload: dict[str, Any] = {}
    simulation_cells = 0
    data_cells = 0
    ncb_cells = 0
    max_norm_deviation = 0.0
    try:
        ncb_scale = ncb_tail_scale(roots["data"], STAGES["before_ncb"]["pt_index"])
        for variable, settings in VARIABLES.items():
            stage_payload: dict[str, Any] = {}
            for stage_name, stage in STAGES.items():
                histogram = (
                    f"{settings['root_base']}_eta0_pt{stage['pt_index']}_cut{stage['cut_index']}"
                )
                samples: dict[str, Any] = {}
                for lane in ("signal", "inclusive", "data"):
                    curve = curve_from_histogram(
                        roots[lane],
                        histogram,
                        xmin=settings["xmin"],
                        xmax=settings["xmax"],
                        rebin=settings["rebin"],
                    )
                    samples[lane] = curve
                    max_norm_deviation = max(
                        max_norm_deviation,
                        abs(curve["display_integral_after_unit_normalization"] - 1.0),
                    )
                    if lane in {"signal", "inclusive"}:
                        simulation_cells += 1
                    else:
                        data_cells += 1

                if stage_name == "before_ncb":
                    ncb_histogram = (
                        f"{settings['root_base']}_eta0_pt{stage['pt_index']}_cut4"
                    )
                    curve = curve_from_histogram(
                        roots["data"],
                        ncb_histogram,
                        xmin=settings["xmin"],
                        xmax=settings["xmax"],
                        rebin=settings["rebin"],
                    )
                    curve["tail_match_scale"] = ncb_scale["scale"]
                    curve["tail_scaled_values"] = (
                        np.asarray(curve["values"]) * ncb_scale["scale"]
                    ).tolist()
                    curve["tail_scaled_errors"] = (
                        np.asarray(curve["errors"]) * ncb_scale["scale"]
                    ).tolist()
                    samples["ncb_tagged_data"] = curve
                    ncb_cells += 1

                # All source lanes must have identical visible binning.
                reference_edges = np.asarray(samples["signal"]["edges"])
                for lane, curve in samples.items():
                    if not np.array_equal(reference_edges, np.asarray(curve["edges"])):
                        raise ValueError(
                            f"visible binning mismatch for {variable}/{stage_name}/{lane}"
                        )
                    for field in ("values", "errors", "edges", "centers"):
                        if not np.all(np.isfinite(np.asarray(curve[field], dtype=float))):
                            raise ValueError(
                                f"non-finite {field} in {variable}/{stage_name}/{lane}"
                            )

                stage_payload[stage_name] = {
                    **stage,
                    "histogram": histogram,
                    "eta_index": 0,
                    "eta_range": [-0.7, 0.7],
                    "samples": samples,
                }

            variables_payload[variable] = {
                "root_base": settings["root_base"],
                "axis_label": settings["axis_label"],
                "display": {
                    "xmin": settings["xmin"],
                    "xmax": settings["xmax"],
                    "rebin": settings["rebin"],
                },
                "stages": stage_payload,
            }
    finally:
        for root_file in roots.values():
            root_file.close()

    expected_simulation_cells = len(VARIABLES) * len(STAGES) * 2
    if simulation_cells != expected_simulation_cells:
        raise RuntimeError(
            f"extracted {simulation_cells} simulation cells, expected {expected_simulation_cells}"
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    matrix_path = output_dir / "ppg12_sdcc_shower_shape_reference_matrix.json"
    provenance_path = output_dir / "provenance.json"
    generated_at = datetime.now(timezone.utc).isoformat()
    validation = {
        "status": "passed",
        "variables": len(VARIABLES),
        "stages_per_variable": len(STAGES),
        "expected_signal_plus_inclusive_cells": expected_simulation_cells,
        "extracted_signal_plus_inclusive_cells": simulation_cells,
        "extracted_data_cells": data_cells,
        "extracted_ncb_tagged_data_cells": ncb_cells,
        "identical_visible_binning_within_each_cell": True,
        "all_values_and_errors_finite": True,
        "maximum_unit_integral_deviation": max_norm_deviation,
    }
    matrix_payload = {
        "schema_version": "PPG12_SDCC_SHOWER_SHAPE_REFERENCE_MATRIX_V1",
        "created_utc": generated_at,
        "contract": {
            "observable_scope": "eight PPG12 baseV3E shower-shape inputs",
            "samples": {
                "signal": "PPG12 combined photon+jet simulation",
                "inclusive": "PPG12 combined inclusive-jet simulation",
                "data": "PPG12 data source used only for historical source-recovery panels",
                "ncb_tagged_data": "historical cut4 non-collision-background-tagged data diagnostic",
            },
            "projection": "TH2 ProjectionX over all Y bins, including Y underflow and overflow",
            "normalization": "unit integral over the displayed X range after the historical PPG12 rebin",
            "x_underflow_overflow": "excluded from the normal-bin projection and unit integral",
            "arrays": "only visible X bins are serialized",
            "historical_macro": "ppg12codeGit/plotting/plot_showershapes_selections.C",
            "stage_order": list(STAGES),
            "variable_order": list(VARIABLES),
            "ncb_tail_normalization": ncb_scale,
        },
        "sources": source_records,
        "canonical_cut_boundary": CUT_REPAIR_COMMIT,
        "archived_shower_shape_tight_boundary": SHOWER_SHAPE_DIAGNOSTIC_BOUNDARY,
        "variables": variables_payload,
        "validation": validation,
    }
    matrix_path.write_text(
        json.dumps(matrix_payload, indent=2, sort_keys=False, allow_nan=False) + "\n"
    )
    matrix_hash = sha256(matrix_path)

    provenance_payload = {
        "schema_version": "PPG12_SDCC_SHOWER_SHAPE_REFERENCE_PROVENANCE_V1",
        "created_utc": generated_at,
        "matrix": {
            "path": str(matrix_path),
            "sha256": matrix_hash,
            "validation": validation,
        },
        "source_files": source_records,
        "canonical_cut_boundary": CUT_REPAIR_COMMIT,
        "archived_shower_shape_tight_boundary": SHOWER_SHAPE_DIAGNOSTIC_BOUNDARY,
        "timestamp_interpretation": {
            "before_ncb_and_after_ncb": (
                "cut0 and cut1 do not apply the tight BDT threshold, so the May-14 "
                "tight-cut change is not directly active in those two stages. The "
                "historical data and simulation inputs nevertheless have different "
                "timestamps and are retained as source-recovery evidence."
            ),
            "tight_id": (
                "The archived shower-shape cut2 family uses score > 0.86 - 0.02 E_T, "
                "which differs from the nominal score > 0.815625 - 0.0015625 E_T "
                "definition fixed on May 14. Tight-stage panels are historical "
                "diagnostic evidence, not canonical tight-ID closure."
            ),
            "current_simulation_comparison": (
                "The recovered PPG12 signal and inclusive arrays are the numerical "
                "historical reference. Current-analysis overlays must separately record "
                "the July-16 registered RecoilJets artifact hashes."
            ),
        },
        "audience_language": (
            "Use NCB (non-collision background) in manuscript prose. Raw historical "
            "object and code names containing 'npb' may appear only as provenance."
        ),
    }
    provenance_path.write_text(
        json.dumps(provenance_payload, indent=2, sort_keys=False, allow_nan=False) + "\n"
    )

    print(f"matrix={matrix_path}")
    print(f"matrix_sha256={matrix_hash}")
    print(f"provenance={provenance_path}")
    print(
        "validation="
        f"{simulation_cells}/{expected_simulation_cells} simulation cells; "
        f"{data_cells} data cells; {ncb_cells} NCB cells; "
        f"max_norm_deviation={max_norm_deviation:.3e}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
