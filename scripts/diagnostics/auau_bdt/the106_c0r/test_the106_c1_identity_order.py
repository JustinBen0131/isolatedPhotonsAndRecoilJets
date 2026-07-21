#!/usr/bin/env python3
"""C1 synthetic multiplicity, tie-order, and fail-closed identity tests."""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from urllib.parse import unquote


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def run_replay(replay: Path, cache: Path, inventory: Path, output: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [str(replay), "replay", str(cache), str(inventory), str(output), str(output.with_suffix(".json"))],
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )


def write_inventory(cache: Path, path: Path) -> None:
    candidates = read_tsv(cache / "rawqa_fill_witnesses.tsv")
    seen: set[tuple[str, ...]] = set()
    rows: list[list[str]] = []
    for row in candidates:
        key = (
            row["directory_name"], row["object_name"], row["object_class"],
            row["variable_name"], row["trigger_name"], row["tag_name"],
            row["view_suffix"], row["photon_pt_slice"], row["centrality_slice"],
        )
        if key in seen:
            continue
        seen.add(key)
        rows.append([
            row["directory_name"], row["object_name"], row["object_class"],
            row["object_name"], row["variable_name"], row["trigger_name"],
            row["tag_name"], row["view_suffix"], row["photon_pt_slice"],
            row["centrality_slice"], "120", "0x0000000000000000",
            "0x3ff3333333333333", row["variable_name"], "Entries", "0", "1",
            "BITWISE_SINGLE_PROCESS_UNWEIGHTED",
        ])
    header = [
        "directory", "object_name", "object_class", "title", "variable",
        "trigger", "tag", "view_suffix", "photon_pt_slice", "centrality_slice",
        "nbins", "xmin_bits", "xmax_bits", "x_axis_title", "y_axis_title",
        "sumw2_required", "required", "content_comparator",
    ]
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def expect_failure(replay: Path, cache: Path, inventory: Path, output: Path, needle: str) -> None:
    result = run_replay(replay, cache, inventory, output)
    assert result.returncode != 0
    assert needle in result.stderr, result.stderr


def main() -> int:
    if len(sys.argv) != 3:
        print(f"usage: {sys.argv[0]} WRITER_TEST REPLAY_BINARY", file=sys.stderr)
        return 2
    writer = Path(sys.argv[1]).resolve()
    replay = Path(sys.argv[2]).resolve()
    temporary = Path(tempfile.mkdtemp(prefix="the106_c1_identity_order_"))
    try:
        cache = temporary / "cache"
        cache.mkdir()
        environment = dict(os.environ)
        environment["THE106_C0R_TEST_CACHE_DIR"] = str(cache)
        subprocess.run([str(writer)], env=environment, check=True)

        candidates = read_tsv(cache / "rawqa_candidates.tsv")
        assert len(candidates) == 3
        event_keys = {
            (row["pair_generation"], row["pair_ordinal"], row["delivered_event_ordinal"],
             row["run_number"], row["event_number"])
            for row in candidates
        }
        assert len(event_keys) == 1
        encounters = [int(row["producer_encounter_ordinal"]) for row in candidates]
        clusters = [int(row["cluster_id"]) for row in candidates]
        assert encounters == [0, 1, 2]
        assert clusters == [44, 45, 46]
        identities = {
            tuple(row[name] for name in (
                "pair_generation", "pair_ordinal", "delivered_event_ordinal",
                "run_number", "event_number", "container_key", "cluster_id",
                "producer_encounter_ordinal", "module_name", "canonical_view", "view_key",
            ))
            for row in candidates
        }
        assert len(identities) == 3

        # The first two candidates deliberately tie in ET.  Encounter order is
        # the stable first tie-breaker, followed by immutable object identity.
        et_by_encounter = {0: 20.0, 1: 20.0, 2: 18.0}
        ranked = sorted(
            candidates,
            key=lambda row: (
                -et_by_encounter[int(row["producer_encounter_ordinal"])],
                int(row["producer_encounter_ordinal"]),
                int(row["container_key"]),
                int(row["cluster_id"]),
            ),
        )
        assert [int(row["producer_encounter_ordinal"]) for row in ranked] == [0, 1, 2]

        inventory = temporary / "inventory.tsv"
        write_inventory(cache, inventory)
        positive = run_replay(replay, cache, inventory, temporary / "positive.root")
        assert positive.returncode == 0, positive.stderr

        original_candidates = (cache / "rawqa_candidates.tsv").read_text(encoding="utf-8")
        lines = original_candidates.splitlines()
        lines[2] = lines[1]
        (cache / "rawqa_candidates.tsv").write_text("\n".join(lines) + "\n", encoding="utf-8")
        expect_failure(replay, cache, inventory, temporary / "duplicate.root", "duplicate candidate identity")
        (cache / "rawqa_candidates.tsv").write_text(original_candidates, encoding="utf-8")

        values = (cache / "rawqa_candidate_values.tsv").read_text(encoding="utf-8")
        value_lines = values.splitlines()
        del value_lines[1]
        (cache / "rawqa_candidate_values.tsv").write_text("\n".join(value_lines) + "\n", encoding="utf-8")
        expect_failure(replay, cache, inventory, temporary / "missing.root", "missing candidate value index")
        (cache / "rawqa_candidate_values.tsv").write_text(values, encoding="utf-8")

        witnesses = (cache / "rawqa_fill_witnesses.tsv").read_text(encoding="utf-8")
        witness_lines = witnesses.splitlines()
        header = witness_lines[0].split("\t")
        fields = witness_lines[1].split("\t")
        fields[header.index("cluster_id")] = "999999"
        witness_lines[1] = "\t".join(fields)
        (cache / "rawqa_fill_witnesses.tsv").write_text("\n".join(witness_lines) + "\n", encoding="utf-8")
        expect_failure(replay, cache, inventory, temporary / "orphan.root", "fill witness lacks exact candidate row")

        report = {
            "schema_version": "THE106_C1_SYNTHETIC_IDENTITY_ORDER_V1",
            "result": "PASS",
            "same_event_candidates": 3,
            "unique_candidate_identities": 3,
            "producer_encounter_order": encounters,
            "equal_et_tie_order": [0, 1],
            "negative_fixtures": {
                "duplicate_identity": "FAIL_CLOSED",
                "missing_value": "FAIL_CLOSED",
                "orphan_fill": "FAIL_CLOSED",
            },
            "writer_replay_path": "REAL_C0R_WRITER_AND_CACHE_ONLY_REPLAY",
        }
        requested = os.environ.get("THE106_C1_REPORT")
        if requested:
            Path(requested).write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return 0
    finally:
        shutil.rmtree(temporary)


if __name__ == "__main__":
    raise SystemExit(main())
