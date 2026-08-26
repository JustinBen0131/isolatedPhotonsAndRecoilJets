from __future__ import annotations

import unittest
from pathlib import Path

import numpy as np
import uproot

from photonjet.io.tree_validation import validate


ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "tests" / "fixtures"


class TreeContractTest(unittest.TestCase):
    def test_pp_tree_contract(self) -> None:
        with uproot.open(FIXTURES / "photonjet_trees_pp.root") as root:
            centrality = root["events"]["centrality"].array(library="np")
        self.assertTrue(np.all(centrality == -1.0))
        lines = validate(
            [FIXTURES / "photonjet_trees_pp.root"],
            model_input_count=11,
        )
        self.assertIn("tree_schema_exact=PASS", lines)
        self.assertIn("eventTree_flat_equivalence=PASS", lines)
        self.assertIn("truth_link_type_event_index=PASS", lines)
        self.assertIn(
            "counts=events:2,eventTree:2,photons:4,jets:4,photonJets:8,"
            "truthPhotons:2,truthJets:2,recoTruthLinks:4",
            lines,
        )

    def test_auau_tree_contract(self) -> None:
        lines = validate(
            [FIXTURES / "photonjet_trees_auau.root"],
            model_input_count=14,
        )
        self.assertIn("eventTree_executable_leaders=PASS", lines)
        self.assertIn("numeric_contract=PASS", lines)


if __name__ == "__main__":
    unittest.main()
