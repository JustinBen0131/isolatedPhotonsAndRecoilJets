from __future__ import annotations

import unittest

from scripts.data_prep.recoiljets.reduce_the243_direct_histograms import (
    _direct_inputs,
    _missing_root_object,
    extract_histogram,
    histogram_contract,
    merge_histogram,
)


class FakeAxis:
    def __init__(self, edges):
        self.edges = edges

    def GetNbins(self):
        return len(self.edges) - 1

    def GetBinLowEdge(self, index):
        return self.edges[index - 1]

    def GetBinUpEdge(self, index):
        return self.edges[index]


class FakeHist2:
    def __init__(self, scale=1.0, xedges=(10.0, 20.0, 35.0), yedges=(0.0, 0.5, 1.0)):
        self.scale = scale
        self.xaxis = FakeAxis(xedges)
        self.yaxis = FakeAxis(yedges)

    def GetDimension(self):
        return 2

    def GetXaxis(self):
        return self.xaxis

    def GetYaxis(self):
        return self.yaxis

    def GetBinContent(self, xbin, ybin):
        return self.scale * (10 * xbin + ybin)

    def GetBinError(self, xbin, ybin):
        return self.scale

    def GetEntries(self):
        return 12.0 * self.scale

    def ClassName(self):
        return "TH2F"


class DirectHistogramReducerTest(unittest.TestCase):
    def test_missing_root_object_accepts_null_proxy(self):
        class NullRootProxy:
            def __bool__(self):
                return False

        self.assertTrue(_missing_root_object(None))
        self.assertTrue(_missing_root_object(NullRootProxy()))
        self.assertFalse(_missing_root_object(FakeHist2()))

    def test_system_contracts_are_exact(self):
        pp_top, pp = histogram_contract("pp")
        au_top, au = histogram_contract("auau")
        self.assertEqual(pp_top, "Photon_4_GeV_plus_MBD_NS_geq_1")
        self.assertEqual(au_top, "photon_10_plus_MBD_NS_geq_2_vtx_lt_150")
        self.assertEqual(pp["recoil_xj_region_a"], "h2_unfoldReco_pTgamma_xJ_incl_r04")
        self.assertEqual(
            au["recoil_xj_region_c"],
            "h2_unfoldReco_pTgamma_xJ_incl_sidebandC_r04_cent_0_20",
        )
        self.assertEqual(
            pp["abcd_leading_a_pt15_20"],
            "h_xJpurityLead_isIsolated_isTight_pT_15_20",
        )
        self.assertEqual(
            au["abcd_leading_d_pt25_35"],
            "h_xJpurityLead_notIsolated_notTight_pT_25_35_cent_0_20",
        )

    def test_extract_includes_root_underflow_and_overflow(self):
        row = extract_histogram(FakeHist2())
        self.assertEqual(row["axis_edges"], [[10.0, 20.0, 35.0], [0.0, 0.5, 1.0]])
        self.assertEqual(row["shape_with_flow"], [4, 4])
        self.assertEqual(row["sumw_with_flow"][3][3], 33.0)
        self.assertEqual(row["sumw2_with_flow"][3][3], 1.0)

    def test_merge_adds_sumw_sumw2_and_entries(self):
        left = merge_histogram(None, extract_histogram(FakeHist2(1.0)))
        merged = merge_histogram(left, extract_histogram(FakeHist2(2.0)))
        self.assertEqual(merged["contributing_files"], 2)
        self.assertEqual(merged["entries"], 36.0)
        self.assertEqual(merged["sumw_with_flow"][2][1], 63.0)
        self.assertEqual(merged["sumw2_with_flow"][2][1], 5.0)

    def test_merge_rejects_axis_drift(self):
        left = merge_histogram(None, extract_histogram(FakeHist2()))
        with self.assertRaisesRegex(ValueError, "axis edges differ"):
            merge_histogram(left, extract_histogram(FakeHist2(yedges=(0.0, 0.4, 1.0))))

    def test_direct_input_identity_is_preserved(self):
        from unittest.mock import patch

        path = (
            "/sphenix/tg/tg01/bulk/jbennett/thesisAna/recoiljets/outputs/"
            "the236_schema10_data_prod_20260818_813a2538_user04/pp/shard_0000/"
            "pp_data_000001.root"
        )
        digest = "a" * 64
        with patch("pathlib.Path.stat") as stat:
            stat.return_value.st_size = 123
            rows = _direct_inputs([path], "pp", [digest], [123])
        self.assertEqual(
            rows,
            [{"path": path, "sha256": digest, "size_bytes": 123}],
        )

    def test_direct_input_identity_cardinality_is_exact(self):
        path = (
            "/sphenix/tg/tg01/bulk/jbennett/thesisAna/recoiljets/outputs/"
            "the236_schema10_data_prod_20260818_813a2538_user04/pp/shard_0000/"
            "pp_data_000001.root"
        )
        with self.assertRaisesRegex(ValueError, "one direct-input SHA-256"):
            _direct_inputs([path], "pp", ["a" * 64, "b" * 64], [1])


if __name__ == "__main__":
    unittest.main()
