#!/usr/bin/env python3
"""Add only selected TH1 objects across an exact THE-110 raw ROOT population.

This is a diagnostic accelerator, not the canonical campaign merge ladder.  It
preserves ROOT TH1 additive/Sumw2 semantics while avoiding unrelated objects.
The caller supplies the exact input root, expected file count, namespace, and
object names; exclusions are explicit basenames recorded in the manifest.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path

import ROOT


ROOT.gROOT.SetBatch(True)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def input_digest(paths: list[Path]) -> str:
    digest = hashlib.sha256()
    for path in paths:
        stat = path.stat()
        digest.update(f"{path}\t{stat.st_size}\n".encode())
    return digest.hexdigest()


def edges(hist: ROOT.TH1) -> list[float]:
    axis = hist.GetXaxis()
    values = [float(axis.GetBinLowEdge(i)) for i in range(1, hist.GetNbinsX() + 1)]
    values.append(float(axis.GetBinUpEdge(hist.GetNbinsX())))
    return values


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-root", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--namespace", required=True)
    parser.add_argument("--object", action="append", dest="objects", required=True)
    parser.add_argument("--exclude-path", type=Path, action="append", default=[])
    parser.add_argument("--expected-inputs", type=int, required=True)
    parser.add_argument("--status", required=True)
    args = parser.parse_args()

    excluded = {str(path) for path in args.exclude_path}
    all_roots = sorted(
        path
        for input_root in args.input_root
        for path in input_root.rglob("*.root")
    )
    selected = [path for path in all_roots if str(path) not in excluded]
    seen_excluded = sorted(str(path) for path in all_roots if str(path) in excluded)
    if seen_excluded != sorted(excluded):
        raise RuntimeError(f"explicit exclusions not found exactly once: {seen_excluded} vs {sorted(excluded)}")
    if len(selected) != args.expected_inputs:
        raise RuntimeError(f"selected {len(selected)} inputs, expected {args.expected_inputs}")

    sums: dict[str, ROOT.TH1] = {}
    missing = {name: 0 for name in args.objects}
    input_bytes = 0
    for index, path in enumerate(selected, start=1):
        input_bytes += path.stat().st_size
        root_file = ROOT.TFile.Open(str(path), "READ")
        if not root_file or root_file.IsZombie() or root_file.TestBit(ROOT.TFile.kRecovered):
            raise RuntimeError(f"invalid ROOT input: {path}")
        try:
            for name in args.objects:
                obj = root_file.Get(f"{args.namespace}/{name}")
                if not obj:
                    missing[name] += 1
                    continue
                if not obj.InheritsFrom("TH1"):
                    raise RuntimeError(f"non-TH1 target {args.namespace}/{name} in {path}")
                if name not in sums:
                    clone = obj.Clone(name)
                    clone.SetDirectory(0)
                    sums[name] = clone
                elif not sums[name].Add(obj):
                    raise RuntimeError(f"TH1::Add failed for {args.namespace}/{name} from {path}")
        finally:
            root_file.Close()
        if index % 1000 == 0 or index == len(selected):
            print(f"processed {index}/{len(selected)}", flush=True)

    absent = sorted(set(args.objects) - set(sums))
    if absent:
        raise RuntimeError(f"targets absent from every input: {absent}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    tmp = args.output.with_suffix(args.output.suffix + f".tmp.{os.getpid()}")
    out = ROOT.TFile.Open(str(tmp), "RECREATE")
    if not out or out.IsZombie():
        raise RuntimeError(f"cannot create output: {tmp}")
    directory = out.mkdir(args.namespace)
    directory.cd()
    for name in args.objects:
        sums[name].Write(name)
    out.cd()
    ROOT.TNamed("the110_diagnostic_status", args.status).Write()
    ROOT.TNamed("the110_input_list_sha256", input_digest(selected)).Write()
    out.Close()
    tmp.replace(args.output)

    check = ROOT.TFile.Open(str(args.output), "READ")
    if not check or check.IsZombie() or check.TestBit(ROOT.TFile.kRecovered):
        raise RuntimeError(f"invalid written output: {args.output}")
    check.Close()

    manifest = {
        "schema": "THE110_TARGET_HISTOGRAM_ADDITIVE_DIAGNOSTIC_V1",
        "status": args.status,
        "canonical_merge": False,
        "input_roots": [str(path) for path in args.input_root],
        "input_count": len(selected),
        "input_bytes": input_bytes,
        "input_list_sha256": input_digest(selected),
        "explicit_excluded_paths": sorted(excluded),
        "namespace": args.namespace,
        "objects": {
            name: {
                "missing_input_count_treated_as_zero_like_hadd": missing[name],
                "bin_edges": edges(sums[name]),
                "integral_with_flow": float(sums[name].Integral(0, sums[name].GetNbinsX() + 1)),
                "sumw2_size": int(sums[name].GetSumw2N()),
            }
            for name in args.objects
        },
        "output": str(args.output),
        "output_bytes": args.output.stat().st_size,
        "output_sha256": sha256(args.output),
        "addition_contract": "Each selected raw ROOT once; TH1::Add only; no external weight, normalization, fit, or retuning.",
    }
    args.manifest.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"output": str(args.output), "manifest": str(args.manifest)}, indent=2))


if __name__ == "__main__":
    main()
