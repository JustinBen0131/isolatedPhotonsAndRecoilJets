#include "../../../ppg12codeGit/plotting/plotcommon.h"

void MakePPG12Fig3PuritySimSelectionLabelLeft()
{
  init_plot();
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  const TString input =
      "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/ppg12PhotonYield/"
      "ppg12_photon_yield_v1_data_20260620/reference_roots/fig3_purity_sim/"
      "Photon_final_bdt_nom_mc.root";
  const TString output =
      "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/ppg12PhotonYield/"
      "ppg12_photon_yield_v1_data_20260620/shower_shape_reference_validation/"
      "fig3_purity_sim/"
      "ppg12_ian_fig3_purity_sim_bdt_nom_sdcc_reproduction_label_left.png";

  TFile fdata(input, "READ");
  if (fdata.IsZombie()) {
    Error("MakePPG12Fig3PuritySimSelectionLabelLeft", "could not open %s", input.Data());
    return;
  }

  TGraphErrors *gpurity = (TGraphErrors *) fdata.Get("gpurity");
  TGraphErrors *gpurity_leak = (TGraphErrors *) fdata.Get("gpurity_leak");
  TGraphAsymmErrors *g_purity_truth = (TGraphAsymmErrors *) fdata.Get("g_purity_truth");
  if (!gpurity || !gpurity_leak || !g_purity_truth) {
    Error("MakePPG12Fig3PuritySimSelectionLabelLeft",
          "missing required graph(s): gpurity=%p gpurity_leak=%p g_purity_truth=%p",
          (void *) gpurity, (void *) gpurity_leak, (void *) g_purity_truth);
    return;
  }

  TCanvas *c1 = new TCanvas("c1", "c1", 600, 600);
  frame_et_rec->SetYTitle("Purity");
  frame_et_rec->GetYaxis()->SetRangeUser(0.0, 1.2);
  frame_et_rec->GetXaxis()->SetRangeUser(10, 35);
  frame_et_rec->Draw("axis");

  gpurity->SetMarkerColor(kBlack);
  gpurity->SetMarkerStyle(20);
  gpurity->SetMarkerSize(1.5);
  gpurity->SetLineColor(kBlack);
  gpurity->Draw("P same");

  gpurity_leak->SetMarkerColor(kBlue);
  gpurity_leak->SetMarkerStyle(20);
  gpurity_leak->SetMarkerSize(1.5);
  gpurity_leak->SetLineColor(kBlue);
  gpurity_leak->Draw("P same");

  g_purity_truth->SetMarkerColor(kRed);
  g_purity_truth->SetMarkerStyle(20);
  g_purity_truth->SetMarkerSize(1.5);
  g_purity_truth->SetLineColor(kRed);
  g_purity_truth->Draw("P same");

  myText(0.19, 0.88, 1, strleg1.c_str(), 0.04);
  myText(0.19, 0.83, 1, strleg2.c_str(), 0.04);
  myText(0.19, 0.78, 1, strMC.c_str(), 0.04);
  myText(0.19, 0.73, 1, "bdt_nom", 0.04);

  myMarkerLineText(0.30, 0.25, 1, kBlack, 20, kBlack, 1,
                   "w/o signal leakage correction", 0.05, true);
  myMarkerLineText(0.30, 0.20, 1, kBlue, 20, kBlue, 1,
                   "w/ signal leakage correction", 0.05, true);
  myMarkerLineText(0.30, 0.30, 1, kRed, 20, kRed, 1, "truth", 0.05, true);

  c1->SaveAs(output);
  printf("wrote %s\n", output.Data());
}
