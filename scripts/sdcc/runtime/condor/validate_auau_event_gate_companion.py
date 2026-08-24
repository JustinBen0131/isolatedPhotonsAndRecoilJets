#!/usr/bin/env python3
"""Validate one inline AuAuEventGateV1 companion against its schema-10 base."""

from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import json
import os
from pathlib import Path
from typing import Any


GATE_SCHEMA = "THE243AuAuScaledPhoton10MinimumBiasGateV1"
RECEIPT_SCHEMA = "AuAuInlineEventGateValidationReceiptV1"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def tree_rows(tree: Any, fields: tuple[str, ...]) -> list[dict[str, int]]:
    return [
        {field: int(getattr(entry, field)) for field in fields}
        for entry in tree
    ]


def named_title(root_file: Any, key: str) -> str:
    value = root_file.Get(key)
    if value is None or not value.InheritsFrom("TNamed"):
        return ""
    return str(value.GetTitle())


def validate(args: argparse.Namespace) -> dict[str, Any]:
    try:
        import ROOT  # type: ignore[import-not-found]
    except ImportError as error:  # pragma: no cover - SDCC integration
        raise RuntimeError("PyROOT is required for event-gate validation") from error

    gate_file = ROOT.TFile.Open(str(args.gate), "READ")
    if not gate_file or gate_file.IsZombie() or gate_file.TestBit(ROOT.TFile.kRecovered):
        raise ValueError(f"invalid event-gate ROOT: {args.gate}")
    try:
        if named_title(gate_file, "schema") != GATE_SCHEMA:
            raise ValueError("event-gate schema differs")
        if named_title(gate_file, "contract_family") != "AuAuEventGateV1":
            raise ValueError("event-gate contract family differs")
        if named_title(gate_file, "producer") != "INLINE_SCHEMA10_AUAU_BASE_EVENT_LOOP_V1":
            raise ValueError("event-gate producer differs")
        if named_title(gate_file, "row_id") != args.row_id:
            raise ValueError("event-gate row identity differs")

        summary_tree = gate_file.Get("AuAuEventGateSummaryV1")
        gate_tree = gate_file.Get("AuAuEventGateV1")
        if summary_tree is None or gate_tree is None or int(summary_tree.GetEntries()) != 1:
            raise ValueError("event-gate trees differ")
        summary = tree_rows(
            summary_tree,
            (
                "run",
                "photon10_bit",
                "source_pairs",
                "processed_events",
                "gl1_missing_events",
                "minimum_bias_missing_events",
                "photon10_pass_events",
                "minimum_bias_pass_events",
                "gate_pass_events",
            ),
        )[0]
        gate_rows = tree_rows(
            gate_tree,
            (
                "run",
                "event_count",
                "event_header_sequence",
                "raw_trigger_bits",
                "live_trigger_bits",
                "scaled_trigger_bits",
                "photon10_bit",
                "photon10_pass",
                "minimum_bias_pass",
            ),
        )
    finally:
        gate_file.Close()

    if summary["run"] != args.run or summary["photon10_bit"] != args.photon10_bit:
        raise ValueError("run or Photon10 bit differs")
    if summary["source_pairs"] != args.expected_source_pairs:
        raise ValueError("source-pair count differs")
    if summary["processed_events"] != args.expected_processed:
        raise ValueError("processed-event count differs")
    if summary["gl1_missing_events"] != 0 or summary["minimum_bias_missing_events"] != 0:
        raise ValueError("required GL1 or MinimumBiasInfo node was missing")
    if summary["gate_pass_events"] != len(gate_rows):
        raise ValueError("gate-pass count differs from AuAuEventGateV1 entries")

    event_counts = [row["event_count"] for row in gate_rows]
    if len(event_counts) != len(set(event_counts)):
        raise ValueError("duplicate event_count in event-gate pass list")
    bit_mask = 1 << args.photon10_bit
    for row in gate_rows:
        if (
            row["run"] != args.run
            or row["photon10_bit"] != args.photon10_bit
            or row["photon10_pass"] != 1
            or row["minimum_bias_pass"] != 1
            or not (row["scaled_trigger_bits"] & bit_mask)
        ):
            raise ValueError("pass-list row violates scaled Photon10 plus minimum-bias contract")

    base_file = ROOT.TFile.Open(str(args.base), "READ")
    if not base_file or base_file.IsZombie() or base_file.TestBit(ROOT.TFile.kRecovered):
        raise ValueError(f"invalid schema-10 base ROOT: {args.base}")
    try:
        foundation = base_file.Get("ReplayFoundationV1")
        event_tree = foundation.Get("RJEventV1") if foundation else None
        if event_tree is None:
            raise ValueError("schema-10 base has no ReplayFoundationV1/RJEventV1")
        base_rows = tree_rows(
            event_tree,
            (
                "event_sequence",
                "physical_event_sequence",
                "run",
                "trigger_bits",
                "live_trigger_bits",
            ),
        )
    finally:
        base_file.Close()

    base_by_ordinal: dict[int, list[dict[str, int]]] = defaultdict(list)
    for row in base_rows:
        base_by_ordinal[row["event_sequence"]].append(row)
    if any(len(rows) != 1 for rows in base_by_ordinal.values()):
        raise ValueError("schema-10 base event_sequence is not unique")

    run_mismatches: list[int] = []
    physical_mismatches: list[int] = []
    witness_mismatches: list[int] = []
    intersection = 0
    for gate_row in gate_rows:
        matches = base_by_ordinal.get(gate_row["event_count"], [])
        if not matches:
            continue
        base_row = matches[0]
        if (
            base_row["trigger_bits"] != gate_row["raw_trigger_bits"]
            or base_row["live_trigger_bits"] != gate_row["live_trigger_bits"]
        ):
            witness_mismatches.append(gate_row["event_count"])
            continue
        intersection += 1
        if base_row["run"] not in (0, args.run):
            run_mismatches.append(gate_row["event_count"])
        if base_row["physical_event_sequence"] not in (
            0,
            gate_row["event_header_sequence"],
        ):
            physical_mismatches.append(gate_row["event_count"])
    if run_mismatches or physical_mismatches or witness_mismatches:
        raise ValueError("event-gate to schema-10 witness join differs")
    if intersection != len(base_rows):
        raise ValueError(
            "event-gate companion does not cover every retained schema-10 base event"
        )

    payload = {
        "schema": RECEIPT_SCHEMA,
        "status": "PASS",
        "row_id": args.row_id,
        "run": args.run,
        "photon10_bit": args.photon10_bit,
        "base_path": str(args.base),
        "base_sha256": sha256(args.base),
        "gate_path": str(args.gate),
        "gate_sha256": sha256(args.gate),
        "summary": summary,
        "join": {
            "base_events": len(base_rows),
            "gate_events": len(gate_rows),
            "intersection_events": intersection,
            "run_mismatches": 0,
            "physical_event_sequence_mismatches": 0,
            "trigger_witness_mismatches": 0,
            "coordinate": "descriptor_row_plus_global_event_ordinal_plus_raw_live_witness",
        },
        "publication_contract": "BASE_ROOT_PUBLISHED_LAST_AFTER_COMPANION_V1",
    }
    descriptor = os.open(args.receipt, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o440)
    with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
        json.dump(payload, stream, sort_keys=True, separators=(",", ":"))
        stream.write("\n")
    return payload


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gate", type=Path, required=True)
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--row-id", required=True)
    parser.add_argument("--run", type=int, required=True)
    parser.add_argument("--photon10-bit", type=int, default=22)
    parser.add_argument("--expected-processed", type=int, required=True)
    parser.add_argument("--expected-source-pairs", type=int, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    payload = validate(parse_args())
    print(json.dumps(payload, sort_keys=True, separators=(",", ":")))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
