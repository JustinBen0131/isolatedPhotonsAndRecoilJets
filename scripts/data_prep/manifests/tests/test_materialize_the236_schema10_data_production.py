#!/usr/bin/env python3
"""Regression tests for the THE-236 data materializer.

Every case here corresponds to a defect that reached production in campaign
the236_schema10_data_prod_20260818_813a2538_user04 and cost real farm time or
real science:

  pp run-boundary crossing   1,605 of 10,007 pp rows, up to 8 h burned each
  calo half never checked    1,956 Au+Au rows, ~20 s each
  records written pre-filter latent processed-vs-expected mismatch
  silent row disappearance   11 whole GRL runs left the campaign unrecorded
"""

from __future__ import annotations

import importlib.util
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

REPO = next(p for p in Path(__file__).resolve().parents if (p / "AGENTS.md").is_file())
TARGET = REPO / "scripts/sdcc/workflows/submit/materialize_the236_schema10_data_production.py"

spec = importlib.util.spec_from_file_location("mat", TARGET)
mat = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(mat)


def jet(run: int, seg: int) -> str:
    return f"DST_Jet_run3auau_pro001_2025p012_v001-{run:08d}-{seg:05d}.root"


def calo(run: int, seg: int) -> str:
    return f"DST_JETCALO_run3auau_pro001_2025p012_v001-{run:08d}-{seg:05d}.root"


class MaterializerContract(unittest.TestCase):
    def build(self, *, system, runs, events, resolvable, prefix, segments=25, restrict=None):
        """Materialize `runs` and return (rows, records-content-by-row)."""

        self.tmp = TemporaryDirectory()
        root = Path(self.tmp.name)
        list_dir = root / "lists"
        list_dir.mkdir()
        for run in runs:
            lines = [f"{jet(run, s)}\t{calo(run, s)}" for s in range(segments)]
            (list_dir / f"{prefix}{run:08d}.list").write_text("\n".join(lines) + "\n")

        rows = mat.materialize_system(
            system=system,
            list_dir=list_dir,
            packet_root=root / "packet",
            output_root=root / "out",
            evidence_root=root / "ev",
            events=events,
            list_prefix=prefix,
            resolvable=resolvable,
            restrict=restrict,
        )
        contents = {
            r["row_id"]: Path(r["records_path"]).read_text().splitlines() for r in rows
        }
        return rows, contents

    def setUp(self):
        for name in ("_SKIPPED", "_REASONS", "_RECOVERED", "_DROPPED"):
            setattr(mat, name, {})

    def tearDown(self):
        if getattr(self, "tmp", None):
            self.tmp.cleanup()

    def test_pp_groups_never_span_runs(self):
        """Fun4AllSyncManager refuses mixed run numbers, so a group that spans a
        run boundary is guaranteed to die partway through."""

        runs = [50070, 50073]
        allnames = {n(r, s) for r in runs for s in range(25) for n in (jet, calo)}
        rows, contents = self.build(
            system="pp", runs=runs, segments=25,
            events={jet(r, s): 100 for r in runs for s in range(25)},
            resolvable=allnames, prefix="dst_ppg12_pair-",
        )
        self.assertTrue(rows)
        for row_id, lines in contents.items():
            seen = {line.split("\t")[0].split("-")[1] for line in lines}
            self.assertEqual(len(seen), 1, f"{row_id} spans runs {seen}")

    def test_missing_calo_half_discards_the_source(self):
        """Checking only the jet is why 1,956 Au+Au rows died at runtime."""

        run = 76153
        names = {jet(run, s) for s in range(25)}
        names |= {calo(run, s) for s in range(25) if s != 7}  # calo seg 7 removed
        rows, contents = self.build(
            system="auau", runs=[run], segments=25,
            events={jet(run, s): 100 for s in range(25)},
            resolvable=names, prefix="dst_auau_jet_pair-",
        )
        shipped = [line.split("\t")[0] for lines in contents.values() for line in lines]
        self.assertNotIn(jet(run, 7), shipped)
        self.assertEqual(len(shipped), 24)

    def test_records_match_expected_processed(self):
        """The records file and expected_processed must describe the SAME set."""

        run = 68491
        names = {jet(run, s) for s in range(25)} | {calo(run, s) for s in range(25)}
        counts = {jet(run, s): 100 for s in range(25)}
        del counts[jet(run, 3)]  # no catalog count, and unrecoverable in this env
        rows, contents = self.build(
            system="auau", runs=[run], segments=25,
            events=counts, resolvable=names, prefix="dst_auau_jet_pair-",
        )
        for row in rows:
            lines = contents[row["row_id"]]
            self.assertEqual(
                row["expected_processed"], 100 * len(lines),
                f"{row['row_id']} declares {row['expected_processed']} for {len(lines)} sources",
            )
            self.assertEqual(row["source_count"], len(lines))

    def test_restrict_sources_prevents_double_counting(self):
        """A recovery rerun must cover ONLY sources no successful row consumed.
        Reprocessing an already-published source would double count its events
        in the final dataset, which is a silent physics error, not a waste."""

        run = 47289
        todo = {jet(run, s) for s in range(5, 25)}          # 20 still to do
        names = {jet(run, s) for s in range(25)} | {calo(run, s) for s in range(25)}
        rows, contents = self.build(
            system="pp", runs=[run], segments=25,
            events={jet(run, s): 100 for s in range(25)},
            resolvable=names, prefix="dst_ppg12_pair-", restrict=todo,
        )
        shipped = [line.split("\t")[0] for lines in contents.values() for line in lines]
        self.assertEqual(sorted(shipped), sorted(todo))
        self.assertEqual(len(rows), 1, "20 survivors should regroup into one compact row")
        for done in (jet(run, 0), jet(run, 4)):
            self.assertNotIn(done, shipped, "already-published source must not be reprocessed")

    def test_fully_unusable_row_is_recorded_not_silently_dropped(self):
        """11 whole GRL runs vanished leaving only a count in a manifest."""

        run = 76098
        mat._DROPPED = {}
        rows, _ = self.build(
            system="auau", runs=[run], segments=10,
            events={}, resolvable=set(), prefix="dst_auau_jet_pair-",
        )
        self.assertEqual(rows, [])
        self.assertTrue(mat._DROPPED.get("auau"), "a row that vanishes must be recorded")
        self.assertEqual(mat._DROPPED["auau"][0]["run"], f"{run:08d}")


if __name__ == "__main__":
    unittest.main()
