#!/usr/bin/env python3
"""Build pp current-IAN manifests without opening every ROOT file.

The full pipeline manifest builder performs a per-file ROOT QA pass, which is
useful but too slow when we already have accepted extraction/stitching roots.
This helper keeps the manifests full-statistics and limits ROOT readback to a
small representative probe per sample.
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
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import uproot


SIGNAL_SAMPLES = ("run28_photonjet5", "run28_photonjet10", "run28_photonjet20")
TRAIN_JETS = ("run28_jet8", "run28_jet12", "run28_jet20", "run28_jet30")
FULL_INCLUSIVE_JETS = ("run28_jet8", "run28_jet12", "run28_jet20", "run28_jet30", "run28_jet40")
ALL_SAMPLES = SIGNAL_SAMPLES + FULL_INCLUSIVE_JETS


def infer_sample(path: str) -> str:
    for sample in sorted(set(ALL_SAMPLES), key=len, reverse=True):
        aliases = (sample, sample.replace("run28_", ""))
        if any(alias in path for alias in aliases):
            return sample
    return "unknown"


def list_roots(root: Path) -> list[str]:
    return sorted(str(path) for path in root.rglob("*.root") if path.is_file())


def probe_file(path: str, tree_name: str) -> dict:
    with uproot.open(path) as handle:
        if tree_name not in handle:
            return {"path": path, "ok": False, "error": f"missing tree {tree_name}"}
        tree = handle[tree_name]
        keys = set(tree.keys())
        required = {"is_signal", "run", "evt", "npb_score"}
        missing = sorted(required - keys)
        if missing:
            return {"path": path, "ok": False, "error": "missing branches " + ",".join(missing)}
        arrays = tree.arrays(sorted(required), library="np")
    labels = arrays["is_signal"].astype("int32")
    npb = arrays["npb_score"].astype("float64")
    run = arrays["run"].astype("int64")
    evt = arrays["evt"].astype("int64")
    real_npb = np.isfinite(npb) & (npb >= 0.0) & (npb <= 1.0)
    unique_events = len(set(zip(run.tolist(), evt.tolist())))
    return {
        "path": path,
        "ok": True,
        "entries": int(len(labels)),
        "signal_rows": int(np.sum(labels == 1)),
        "background_rows": int(np.sum(labels == 0)),
        "real_npb_rows": int(np.sum(real_npb)),
        "unique_events": int(unique_events),
    }


def write_manifest(path: Path, rows: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(rows) + ("\n" if rows else ""))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", required=True, type=Path)
    parser.add_argument("--photon-root", required=True, type=Path)
    parser.add_argument("--jet-root", required=True, type=Path)
    parser.add_argument("--tree", default="AuAuPhotonIDTrainingTree")
    parser.add_argument("--probe-per-sample", default=3, type=int)
    parser.add_argument("--max-probe-search", default=50, type=int)
    args = parser.parse_args()

    if not args.photon_root.is_dir():
        raise SystemExit(f"missing photon root: {args.photon_root}")
    if not args.jet_root.is_dir():
        raise SystemExit(f"missing jet root: {args.jet_root}")

    photon_roots = list_roots(args.photon_root)
    jet_roots = list_roots(args.jet_root)
    by_sample: dict[str, list[str]] = defaultdict(list)
    for path in photon_roots:
        sample = infer_sample(path)
        if sample in SIGNAL_SAMPLES:
            by_sample[sample].append(path)
    for path in jet_roots:
        sample = infer_sample(path)
        if sample in FULL_INCLUSIVE_JETS:
            by_sample[sample].append(path)

    train_paths: list[str] = []
    signal_paths: list[str] = []
    inclusive_paths: list[str] = []
    for sample in SIGNAL_SAMPLES:
        paths = by_sample.get(sample, [])
        signal_paths.extend(paths)
        train_paths.extend(paths)
    for sample in TRAIN_JETS:
        paths = by_sample.get(sample, [])
        train_paths.extend(paths)
        inclusive_paths.extend(paths)
    for sample in FULL_INCLUSIVE_JETS:
        if sample not in TRAIN_JETS:
            inclusive_paths.extend(by_sample.get(sample, []))

    train_manifest = args.run_root / "training_roots_currentIAN.list"
    signal_manifest = args.run_root / "signal_roots_currentIAN.list"
    inclusive_manifest = args.run_root / "inclusive_roots_currentIAN.list"
    write_manifest(train_manifest, sorted(train_paths))
    write_manifest(signal_manifest, sorted(signal_paths))
    write_manifest(inclusive_manifest, sorted(inclusive_paths))

    probes: dict[str, list[dict]] = {}
    failures: list[dict] = []
    for sample in ALL_SAMPLES:
        probes[sample] = []
        searched = 0
        for path in by_sample.get(sample, [])[: args.max_probe_search]:
            searched += 1
            try:
                result = probe_file(path, args.tree)
            except Exception as exc:  # noqa: BLE001 - persisted as QA evidence
                result = {"path": path, "ok": False, "error": str(exc)}
            if result.get("ok"):
                probes[sample].append(result)
                if len(probes[sample]) >= args.probe_per_sample:
                    break
            else:
                failures.append({"sample": sample, **result})
        if len(probes[sample]) < min(args.probe_per_sample, len(by_sample.get(sample, []))):
            failures.append(
                {
                    "sample": sample,
                    "ok": False,
                    "error": f"only {len(probes[sample])} good probes after {searched} files",
                }
            )

    sample_file_counts = {sample: len(by_sample.get(sample, [])) for sample in ALL_SAMPLES}
    summary = {
        "schema": "PP_CURRENTIAN_FAST_MANIFEST_QA_V1",
        "run_root": str(args.run_root),
        "photon_root": str(args.photon_root),
        "jet_root": str(args.jet_root),
        "tree": args.tree,
        "sample_file_counts": sample_file_counts,
        "training_files": len(train_paths),
        "signal_files": len(signal_paths),
        "inclusive_files": len(inclusive_paths),
        "training_manifest": str(train_manifest),
        "signal_manifest": str(signal_manifest),
        "inclusive_manifest": str(inclusive_manifest),
        "probe_per_sample": args.probe_per_sample,
        "probes": probes,
        "failures": failures[:50],
    }
    missing = [sample for sample in ALL_SAMPLES if sample_file_counts.get(sample, 0) <= 0]
    bad_probe_samples = [sample for sample in ALL_SAMPLES if not probes.get(sample)]
    summary["missing_samples"] = missing
    summary["bad_probe_samples"] = bad_probe_samples

    qa_path = args.run_root / "manifest_tree_qa_currentIAN_fast.json"
    qa_path.parent.mkdir(parents=True, exist_ok=True)
    qa_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, sort_keys=True))

    fatal = []
    if missing:
        fatal.append("missing samples " + ",".join(missing))
    if bad_probe_samples:
        fatal.append("bad probe samples " + ",".join(bad_probe_samples))
    if not train_paths or not signal_paths or not inclusive_paths:
        fatal.append("empty manifest")
    if fatal:
        raise SystemExit("fast manifest QA failed: " + "; ".join(fatal))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
