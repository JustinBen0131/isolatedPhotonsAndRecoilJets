#!/usr/bin/env python3
"""Compare Justin/RecoilJets pp stitching histograms against the PPG12 contract.

This intentionally treats PPG12/Shuhang files as references only.  The
"analysis" side must come from RecoilJets in-situ ROOT outputs and must contain
the ppg12TruthSpectrum metadata histograms written by src/RecoilJets.cc.
"""

from __future__ import annotations
# Keep purpose-folder helpers runnable when invoked directly.
import sys as _codex_sys
from pathlib import Path as _CodexPath
_CODEX_THIS_FILE = _CodexPath(__file__).resolve()
_CODEX_SCRIPTS_DIR = next((p for p in _CODEX_THIS_FILE.parents if p.name == "scripts"), _CODEX_THIS_FILE.parent)
_CODEX_SCRIPTS_DIR_STR = str(_CODEX_SCRIPTS_DIR)
if _CODEX_SCRIPTS_DIR_STR not in _codex_sys.path:
    _codex_sys.path.append(_CODEX_SCRIPTS_DIR_STR)
del _CODEX_THIS_FILE, _CODEX_SCRIPTS_DIR, _CODEX_SCRIPTS_DIR_STR

import argparse
import csv
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

try:
    import uproot
except ModuleNotFoundError:  # keep --help usable on lightweight local Pythons
    uproot = None


@dataclass(frozen=True)
class Sample:
    group: str
    sample: str
    xsec_pb: float
    window_low: float
    window_high: float
    sample_bin: int
    all_hist: str
    kept_hist: str
    metadata_hist: str


SAMPLES: tuple[Sample, ...] = (
    Sample("photon", "run28_photonjet5", 146359.3, 0.0, 14.0, 1, "SIM/h_ppPhotonStitch_ppg12Fig5_maxPhotonPt_all", "SIM/h_ppPhotonStitch_ppg12Fig5_maxPhotonPt_kept", "SIM/h_ppPhotonStitch_ppg12TruthSpectrum_metadata"),
    Sample("photon", "run28_photonjet10", 6944.675, 14.0, 22.0, 2, "SIM/h_ppPhotonStitch_ppg12Fig5_maxPhotonPt_all", "SIM/h_ppPhotonStitch_ppg12Fig5_maxPhotonPt_kept", "SIM/h_ppPhotonStitch_ppg12TruthSpectrum_metadata"),
    Sample("photon", "run28_photonjet20", 130.4461, 22.0, 200.0, 3, "SIM/h_ppPhotonStitch_ppg12Fig5_maxPhotonPt_all", "SIM/h_ppPhotonStitch_ppg12Fig5_maxPhotonPt_kept", "SIM/h_ppPhotonStitch_ppg12TruthSpectrum_metadata"),
    Sample("jet", "run28_jet8", 1.3013e7, 9.0, 14.0, 2, "SIM/h_ppInclusiveJetStitch_ppg12TruthSpectrum_r04_maxTruthJetPt_all", "SIM/h_ppInclusiveJetStitch_ppg12TruthSpectrum_r04_maxTruthJetPt_kept", "SIM/h_ppInclusiveJetStitch_ppg12TruthSpectrum_metadata"),
    Sample("jet", "run28_jet12", 1.4903e6, 14.0, 21.0, 3, "SIM/h_ppInclusiveJetStitch_ppg12TruthSpectrum_r04_maxTruthJetPt_all", "SIM/h_ppInclusiveJetStitch_ppg12TruthSpectrum_r04_maxTruthJetPt_kept", "SIM/h_ppInclusiveJetStitch_ppg12TruthSpectrum_metadata"),
    Sample("jet", "run28_jet20", 6.2623e4, 21.0, 32.0, 4, "SIM/h_ppInclusiveJetStitch_ppg12TruthSpectrum_r04_maxTruthJetPt_all", "SIM/h_ppInclusiveJetStitch_ppg12TruthSpectrum_r04_maxTruthJetPt_kept", "SIM/h_ppInclusiveJetStitch_ppg12TruthSpectrum_metadata"),
    Sample("jet", "run28_jet30", 2.5298e3, 32.0, 42.0, 5, "SIM/h_ppInclusiveJetStitch_ppg12TruthSpectrum_r04_maxTruthJetPt_all", "SIM/h_ppInclusiveJetStitch_ppg12TruthSpectrum_r04_maxTruthJetPt_kept", "SIM/h_ppInclusiveJetStitch_ppg12TruthSpectrum_metadata"),
    Sample("jet", "run28_jet40", 1.3553e2, 42.0, 100.0, 6, "SIM/h_ppInclusiveJetStitch_ppg12TruthSpectrum_r04_maxTruthJetPt_all", "SIM/h_ppInclusiveJetStitch_ppg12TruthSpectrum_r04_maxTruthJetPt_kept", "SIM/h_ppInclusiveJetStitch_ppg12TruthSpectrum_metadata"),
)


def detect_sample(path: Path) -> Sample | None:
    text = str(path).lower()
    name = path.name.lower()
    # The comparator is meant to read per-sample RecoilJets outputs.  Final
    # stitched products contain names like photonjet5plus10plus20 and would be
    # misclassified as photonjet5 if they are scanned recursively.
    if "plus" in name or "merged.root" in name or any(part.lower().endswith("merged_sim") for part in path.parts):
        return None
    for sample in SAMPLES:
        sample_name = sample.sample.lower()
        sample_short = sample_name.replace("run28_", "")
        if sample_name in text or sample_short in text:
            return sample
    return None


def read_metadata(root_file: uproot.ReadOnlyDirectory, spec: Sample) -> dict[str, float]:
    if spec.metadata_hist not in root_file:
        raise KeyError(f"missing metadata histogram {spec.metadata_hist}")
    values, _ = root_file[spec.metadata_hist].to_numpy(flow=False)
    if len(values) < 8:
        raise ValueError(f"metadata histogram {spec.metadata_hist} has {len(values)} bins")
    return {
        "events_processed": float(values[0]),
        "window_low_GeV": float(values[1]),
        "window_high_GeV": float(values[2]),
        "upper_edge_inclusive": float(values[3]),
        "bin_width_GeV": float(values[4]),
        "xsec_pb": float(values[5]),
        "truth_def_code": float(values[6]),
        "sample_bin": float(values[7]),
    }


def add_hist(acc: dict[str, np.ndarray], values: np.ndarray, variances: np.ndarray | None, edges: np.ndarray) -> None:
    if "values" not in acc:
        acc["values"] = np.zeros_like(values, dtype=float)
        acc["variances"] = np.zeros_like(values, dtype=float)
        acc["edges"] = edges.astype(float)
    if not np.allclose(acc["edges"], edges):
        raise ValueError("inconsistent histogram bin edges within one sample")
    acc["values"] += values
    acc["variances"] += np.clip(variances if variances is not None else values, 0, None)


def collect(root_dir: Path, allow_legacy: bool) -> tuple[dict[str, dict[str, object]], list[str]]:
    if uproot is None:
        raise RuntimeError("uproot is required to read ROOT files; use the SDCC ML Python or install uproot locally")
    out: dict[str, dict[str, object]] = {}
    problems: list[str] = []
    for path in sorted(root_dir.rglob("*.root")):
        spec = detect_sample(path)
        if spec is None:
            continue
        rec = out.setdefault(
            spec.sample,
            {
                "spec": spec,
                "files": 0,
                "metadata_files": 0,
                "metadata": None,
                "hist_all": {},
                "hist_kept": {},
            },
        )
        try:
            with uproot.open(path) as f:
                if spec.all_hist not in f:
                    problems.append(f"{path}: missing {spec.all_hist}")
                    continue
                if spec.kept_hist not in f:
                    problems.append(f"{path}: missing {spec.kept_hist}")
                    continue
                meta = read_metadata(f, spec)
                if rec["metadata"] is None:
                    rec["metadata"] = dict(meta)
                else:
                    old = rec["metadata"]
                    assert isinstance(old, dict)
                    old["events_processed"] = float(old["events_processed"]) + float(meta["events_processed"])
                rec["metadata_files"] = int(rec["metadata_files"]) + 1
                values_all, edges = f[spec.all_hist].to_numpy(flow=False)
                add_hist(rec["hist_all"], values_all.astype(float), f[spec.all_hist].variances(flow=False), edges.astype(float))
                values_kept, kept_edges = f[spec.kept_hist].to_numpy(flow=False)
                add_hist(rec["hist_kept"], values_kept.astype(float), f[spec.kept_hist].variances(flow=False), kept_edges.astype(float))
                rec["files"] = int(rec["files"]) + 1
        except KeyError as exc:
            if allow_legacy:
                problems.append(f"{path}: legacy/no metadata: {exc}")
                continue
            raise
    return out, problems


def ppg12_bin_center_window(centers: np.ndarray, spec: Sample) -> np.ndarray:
    # Match ppg12codeGit/plotting/plot_combine_uncut.C: lower edge is
    # included and upper edge is excluded, evaluated on the bin center.
    return (centers >= spec.window_low) & (centers < spec.window_high)


def write_outputs(records: dict[str, dict[str, object]], out_csv: Path, out_json: Path, selected_groups: set[str]) -> int:
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    sample_summaries: list[dict[str, object]] = []
    for sample_name, rec in sorted(records.items()):
        spec = rec["spec"]
        assert isinstance(spec, Sample)
        if spec.group not in selected_groups:
            continue
        meta = rec["metadata"]
        hist_all = rec["hist_all"]
        hist_kept = rec["hist_kept"]
        if (
            not isinstance(meta, dict)
            or not isinstance(hist_all, dict)
            or not isinstance(hist_kept, dict)
            or "values" not in hist_all
            or "values" not in hist_kept
        ):
            continue
        all_values = hist_all["values"]
        all_variances = hist_all["variances"]
        kept_values = hist_kept["values"]
        kept_variances = hist_kept["variances"]
        edges = hist_all["edges"]
        kept_edges = hist_kept["edges"]
        assert isinstance(all_values, np.ndarray)
        assert isinstance(all_variances, np.ndarray)
        assert isinstance(kept_values, np.ndarray)
        assert isinstance(kept_variances, np.ndarray)
        assert isinstance(edges, np.ndarray)
        assert isinstance(kept_edges, np.ndarray)
        if not np.allclose(edges, kept_edges):
            raise ValueError(f"{sample_name}: all/kept histogram edges differ")
        centers = 0.5 * (edges[:-1] + edges[1:])
        widths = np.diff(edges)
        window_mask = ppg12_bin_center_window(centers, spec)
        values = np.where(window_mask, all_values, 0.0)
        variances = np.where(window_mask, all_variances, 0.0)
        events_processed = float(meta["events_processed"])
        if events_processed <= 0:
            raise ValueError(f"{sample_name}: non-positive events_processed metadata")
        density = values * spec.xsec_pb / events_processed / widths
        density_err = np.sqrt(np.clip(variances, 0, None)) * spec.xsec_pb / events_processed / widths
        for lo, hi, center, count, all_count, kept_count, windowed, y, ey in zip(
            edges[:-1],
            edges[1:],
            centers,
            values,
            all_values,
            kept_values,
            window_mask,
            density,
            density_err,
        ):
            rows.append(
                {
                    "group": spec.group,
                    "sample": sample_name,
                    "bin_low": f"{lo:.8g}",
                    "bin_high": f"{hi:.8g}",
                    "bin_center": f"{center:.8g}",
                    "raw_kept_events": f"{count:.12g}",
                    "raw_all_events": f"{all_count:.12g}",
                    "raw_event_kept_events": f"{kept_count:.12g}",
                    "ppg12_bin_center_window": int(windowed),
                    "events_processed_metadata": f"{events_processed:.12g}",
                    "xsec_pb": f"{spec.xsec_pb:.12g}",
                    "density_pb_per_gev": f"{y:.12g}",
                    "density_err_pb_per_gev": f"{ey:.12g}",
                    "stitch_window_low": f"{spec.window_low:.8g}",
                    "stitch_window_high": f"{spec.window_high:.8g}",
                    "metadata_bin_width_GeV": f"{float(meta['bin_width_GeV']):.8g}",
                    "metadata_truth_def_code": f"{float(meta['truth_def_code']):.8g}",
                    "source_files": int(rec["files"]),
                }
            )
        sample_summaries.append(
            {
                "group": spec.group,
                "sample": sample_name,
                "source_files": int(rec["files"]),
                "metadata_files": int(rec["metadata_files"]),
                "events_processed_metadata": events_processed,
                "raw_display_integral_ppg12_bin_center_window": float(np.sum(values)),
                "raw_all_integral": float(np.sum(all_values)),
                "raw_event_kept_integral": float(np.sum(kept_values)),
                "xsec_pb": spec.xsec_pb,
                "stitch_window": [spec.window_low, spec.window_high],
                "window_application": "PPG12 plot contract: full _all histogram, zero bins by bin center with lower-inclusive upper-exclusive window",
                "hist_bin_width": float(widths[0]) if len(widths) else math.nan,
                "metadata": meta,
            }
        )
    fields = [
        "group",
        "sample",
        "bin_low",
        "bin_high",
        "bin_center",
        "raw_kept_events",
        "raw_all_events",
        "raw_event_kept_events",
        "ppg12_bin_center_window",
        "events_processed_metadata",
        "xsec_pb",
        "density_pb_per_gev",
        "density_err_pb_per_gev",
        "stitch_window_low",
        "stitch_window_high",
        "metadata_bin_width_GeV",
        "metadata_truth_def_code",
        "source_files",
    ]
    with out_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    out_json.write_text(
        json.dumps(
            {
                "schema": "PP_CURRENTIAN_INSITU_PPG12_STITCH_CONTRACT_V1",
                "normalization": "density = ppg12_bin_center_windowed_all_counts * xsec_pb / events_processed_metadata / bin_width",
                "important_guard": "Do not normalize by histogram integral; metadata event counter is required.",
                "display_contract": "PPG12 plotting macros build the stitched display from the full weighted per-sample truth histogram and zero bins outside the window by bin center. Event-level kept/rejected histograms remain audit metadata, not the displayed Fig.5/6 parity spectrum.",
                "samples": sample_summaries,
            },
            indent=2,
        )
    )
    return len(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root-dir", required=True, type=Path, help="Directory containing Justin/RecoilJets ROOT outputs")
    parser.add_argument("--out-csv", required=True, type=Path)
    parser.add_argument("--out-json", required=True, type=Path)
    parser.add_argument("--groups", default="photon,jet", help="Comma-separated sample groups to require/write, e.g. photon or photon,jet")
    parser.add_argument("--allow-legacy", action="store_true", help="Do not fail on files missing the new metadata")
    args = parser.parse_args()

    selected_groups = {x.strip() for x in args.groups.split(",") if x.strip()}
    known_groups = {s.group for s in SAMPLES}
    unknown_groups = sorted(selected_groups - known_groups)
    if unknown_groups:
        print(f"ERROR unknown groups: {', '.join(unknown_groups)}", file=sys.stderr)
        return 1

    records, problems = collect(args.root_dir, args.allow_legacy)
    missing = [s.sample for s in SAMPLES if s.group in selected_groups and s.sample not in records]
    if missing:
        print(f"ERROR missing samples: {', '.join(missing)}", file=sys.stderr)
        return 1
    if problems and not args.allow_legacy:
        print("\n".join(problems[:50]), file=sys.stderr)
        return 1
    nrows = write_outputs(records, args.out_csv, args.out_json, selected_groups)
    print(args.out_csv)
    print(args.out_json)
    print(f"rows={nrows}")
    if problems:
        print(f"legacy_or_skipped_files={len(problems)}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
