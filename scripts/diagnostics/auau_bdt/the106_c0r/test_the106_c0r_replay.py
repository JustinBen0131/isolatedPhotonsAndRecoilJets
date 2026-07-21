#!/usr/bin/env python3
"""Synthetic cache-only C0-R replay/closure test.

This test processes no DST and invokes no Fun4All code.  It asks the C++
writer fixture to produce a three-candidate same-event cache, freezes a small explicit
inventory, replays that inventory from raw candidate values, and exercises the
bitwise ROOT comparator and a fail-closed witness mutation.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
import shutil
import struct
import subprocess
import sys
import tempfile
from urllib.parse import quote


def bits64(value: float) -> str:
    return "0x" + struct.pack(">d", value).hex()


def encoded(value: str) -> str:
    return quote(value, safe="._-/")


def main() -> int:
    if len(sys.argv) != 3:
        print(f"usage: {sys.argv[0]} WRITER_TEST REPLAY_BINARY", file=sys.stderr)
        return 2
    writer_test = Path(sys.argv[1]).resolve()
    replay = Path(sys.argv[2]).resolve()
    temporary = Path(tempfile.mkdtemp(prefix="the106_c0r_replay_test_"))
    try:
        cache = temporary / "cache"
        cache.mkdir()
        environment = dict(os.environ)
        environment["THE106_C0R_TEST_CACHE_DIR"] = str(cache)
        subprocess.run([str(writer_test)], env=environment, check=True)

        inventory = temporary / "inventory.tsv"
        header = [
            "directory", "object_name", "object_class", "title", "variable",
            "trigger", "tag", "view_suffix", "photon_pt_slice",
            "centrality_slice", "nbins", "xmin_bits", "xmax_bits",
            "x_axis_title", "y_axis_title", "sumw2_required", "required",
            "content_comparator",
        ]
        variables = [
            "weta", "wphi", "weta33", "wphi33", "weta35", "wphi53",
            "et1", "e11e33", "e32e35",
        ]
        rows: list[list[str]] = []
        for tag in ("pre", "tight"):
            for variable in variables:
                object_name = f"h_ss_{variable}_{tag}_pT_5_10_cent_30_50"
                # ROOT parses the semicolon booking string into the persisted
                # object title plus independent axis titles.
                title = f"h_ss_{variable}_{tag}"
                rows.append(
                    [
                        encoded("MBDNS2"), encoded(object_name), encoded("TH1F"),
                        encoded(title), encoded(variable), encoded("MBDNS2"),
                        encoded(tag), encoded("canonical"), "0", "1", "120",
                        bits64(0.0), bits64(1.2), encoded(variable),
                        encoded("Entries"), "0", "1",
                        "BITWISE_SINGLE_PROCESS_UNWEIGHTED",
                    ]
                )
        inventory.write_text(
            "\t".join(header) + "\n" +
            "".join("\t".join(row) + "\n" for row in rows),
            encoding="utf-8",
        )

        output = temporary / "replay.root"
        replay_report = temporary / "replay_report.json"
        subprocess.run(
            [str(replay), "replay", str(cache), str(inventory), str(output),
             str(replay_report)],
            check=True,
        )
        report = json.loads(replay_report.read_text(encoding="utf-8"))
        assert report["result"] == "PASS"
        assert report["candidate_count"] == 3
        assert report["predicted_fill_count"] == 54
        assert report["object_count"] == 18

        comparison = temporary / "comparison.json"
        subprocess.run(
            [str(replay), "compare", str(inventory), str(output), str(output),
             str(comparison)],
            check=True,
        )
        assert json.loads(comparison.read_text(encoding="utf-8"))["result"] == "PASS"

        witness_path = cache / "rawqa_fill_witnesses.tsv"
        witness_lines = witness_path.read_text(encoding="utf-8").splitlines()
        fields = witness_lines[1].split("\t")
        weight_column = witness_lines[0].split("\t").index("weight_bits")
        fields[weight_column] = bits64(2.0)
        witness_lines[1] = "\t".join(fields)
        witness_path.write_text("\n".join(witness_lines) + "\n", encoding="utf-8")
        failed = subprocess.run(
            [str(replay), "replay", str(cache), str(inventory),
             str(temporary / "must_not_close.root"),
             str(temporary / "must_not_close.json")],
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        assert failed.returncode != 0
        assert "lacks exact fill witness" in failed.stderr
        return 0
    finally:
        shutil.rmtree(temporary)


if __name__ == "__main__":
    raise SystemExit(main())
