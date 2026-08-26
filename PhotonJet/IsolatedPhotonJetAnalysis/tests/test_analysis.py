from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
import uproot

from photonjet.analysis.purity import purity_counts
from photonjet.analysis.reduce import RecoilSelection, recoil_histogram, write_recoil_skim


ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "tests" / "fixtures"


class AnalysisTest(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.output = FIXTURES / "photonjet_trees_pp.root"

    def test_histogram_is_deterministic_and_complete(self) -> None:
        selection = RecoilSelection(region="inclusive")
        first = recoil_histogram([self.output], selection)
        second = recoil_histogram([self.output], selection)
        self.assertEqual(first, second)
        self.assertEqual(len(first["sumw"]), len(first["axis"]["edges"]) - 1)
        self.assertEqual(first["selected_pairs"], int(sum(first["sumw"])))

    def test_event_leading_skim_preserves_selection(self) -> None:
        selection = RecoilSelection(region="A")
        skim = Path(self.temporary.name) / "skim.root"
        receipt = write_recoil_skim([self.output], skim, selection)
        with uproot.open(skim) as root:
            self.assertEqual(root["recoilPairs"].num_entries, receipt["selected_pairs"])
            values = root["recoilPairs"].arrays(["photon_et", "jet_pt", "delta_phi"], library="np")
        self.assertTrue(np.all(values["photon_et"] >= selection.photon_et_min))
        self.assertTrue(np.all(values["photon_et"] < selection.photon_et_max))
        self.assertTrue(np.all(values["jet_pt"] > selection.jet_pt_min))
        self.assertTrue(np.all(values["delta_phi"] > selection.delta_phi_min))
        self.assertEqual(receipt["output"]["name"], "skim.root")
        self.assertEqual(receipt["output"]["role"], "recoil_skim")
        serialized = json.dumps(receipt, sort_keys=True)
        self.assertNotIn(self.temporary.name, serialized)
        root_bytes = skim.read_bytes()
        for marker in (
            b"/" + b"Users/",
            b"/" + b"home/",
            b"/private" + b"/" + b"tmp/",
            b"/" + b"tmp/",
            b"/" + b"sphenix/",
            b"/" + b"gpfs",
        ):
            self.assertNotIn(marker, root_bytes)

    def test_purity_boundaries_are_explicit(self) -> None:
        bounded = purity_counts(
            [self.output],
            non_tight_definition="bounded",
            isolation_radius=0.4,
        )
        complement = purity_counts(
            [self.output],
            non_tight_definition="complement",
            isolation_radius=0.4,
        )
        self.assertEqual(set(bounded["weighted_counts"]), set("ABCD"))
        self.assertGreaterEqual(
            complement["weighted_counts"]["C"],
            bounded["weighted_counts"]["C"],
        )


if __name__ == "__main__":
    unittest.main()
