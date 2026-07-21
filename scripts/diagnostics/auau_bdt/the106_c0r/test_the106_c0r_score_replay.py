#!/usr/bin/env python3
"""Fail-closed parser/provenance tests for the bounded score validator.

The positive model evaluation is the real C0-R gate.  These synthetic tests
process no DST and deliberately stop before constructing an RBDT model.
"""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from urllib.parse import quote


def encoded(value: str) -> str:
    return quote(value, safe="._-/")


def replace_metadata(cache: Path, replacements: dict[str, str]) -> None:
    path = cache / "metadata.tsv"
    rows = [line.split("\t") for line in path.read_text(encoding="utf-8").splitlines()]
    for row in rows[1:]:
        if row[0] in replacements:
            row[1] = encoded(replacements[row[0]])
    path.write_text("".join("\t".join(row) + "\n" for row in rows), encoding="utf-8")


def replace_score_fields(cache: Path, replacements: dict[str, str]) -> None:
    path = cache / "score_observations.tsv"
    lines = path.read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t")
    fields = lines[1].split("\t")
    for name, value in replacements.items():
        fields[header.index(name)] = encoded(value)
    lines[1] = "\t".join(fields)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def expect_failure(helper: Path, cache: Path, report: Path, needle: str) -> bytes:
    result = subprocess.run(
        [str(helper), str(cache), str(report)],
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    assert result.returncode == 1, result
    payload = report.read_bytes()
    parsed = json.loads(payload)
    assert parsed["result"] == "FAIL"
    assert needle in parsed["error"]
    assert needle in result.stderr
    return payload


def main() -> int:
    if len(sys.argv) != 4:
        print(f"usage: {sys.argv[0]} WRITER_TEST SCORE_HELPER ROOT_VERSION", file=sys.stderr)
        return 2
    writer_test = Path(sys.argv[1]).resolve()
    helper = Path(sys.argv[2]).resolve()
    root_version = sys.argv[3]
    temporary = Path(tempfile.mkdtemp(prefix="the106_c0r_score_test_"))
    try:
        cache = temporary / "cache"
        cache.mkdir()
        environment = dict(os.environ)
        environment["THE106_C0R_TEST_CACHE_DIR"] = str(cache)
        subprocess.run([str(writer_test)], env=environment, check=True)

        # Make every enclosing authority internally real enough that the
        # deliberately changed score-row model SHA is the first bad boundary.
        model = temporary / "model.root"
        model.write_bytes(b"not-an-rbdt: parser-provenance-test\n")
        model_sha = hashlib.sha256(model.read_bytes()).hexdigest()
        feature_sha = hashlib.sha256(b"minus_zero\none\n").hexdigest()
        runtime = f"ROOT/{root_version}|TMVA::Experimental::RBDT|model_key=myBDT"
        replace_metadata(
            cache,
            {
                "model_reference": str(model),
                "model_sha256": model_sha,
                "feature_order_sha256": feature_sha,
                "runtime_provider_identity": runtime,
            },
        )
        replace_score_fields(
            cache,
            {
                "model_reference": str(model),
                "model_sha256": "0" * 64,
                "feature_order_sha256": feature_sha,
                "runtime_provider_identity": runtime,
            },
        )

        first = expect_failure(
            helper,
            cache,
            temporary / "provenance_1.json",
            "score-row provenance differs",
        )
        second = expect_failure(
            helper,
            cache,
            temporary / "provenance_2.json",
            "score-row provenance differs",
        )
        assert first == second, "identical evidence did not produce byte-identical JSON"

        # Once model provenance is repaired, a changed ordered-name hash must
        # still fail before the intentionally non-RBDT fixture is opened.
        replace_score_fields(cache, {"model_sha256": model_sha})
        replace_metadata(cache, {"feature_order_sha256": "1" * 64})
        replace_score_fields(cache, {"feature_order_sha256": "1" * 64})
        expect_failure(
            helper,
            cache,
            temporary / "feature_order.json",
            "observed ordered feature names do not match",
        )
        return 0
    finally:
        shutil.rmtree(temporary)


if __name__ == "__main__":
    raise SystemExit(main())
