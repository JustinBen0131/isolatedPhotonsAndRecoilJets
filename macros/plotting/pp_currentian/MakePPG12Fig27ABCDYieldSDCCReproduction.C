R__ADD_INCLUDE_PATH(/Users/patsfan753/Desktop/ThesisAnalysis/ppg12codeGit/plotting)
#include "../../../ppg12codeGit/plotting/plotcommon.h"

namespace
{
TH1D *requireHist(TFile &file, const char *name)
{
  auto *hist = dynamic_cast<TH1D *>(file.Get(name));
  if (!hist)
  {
    Error("MakePPG12Fig27ABCDYieldSDCCReproduction", "Missing histogram %s", name);
    return nullptr;
  }
  hist = dynamic_cast<TH1D *>(hist->Clone(Form("%s_clone", name)));
  hist->SetDirectory(nullptr);
  hist->SetStats(false);
  return hist;
}

void styleABCD(TH1D *hist, int color)
{
  hist->SetLineColor(color);
  hist->SetMarkerColor(color);
}
}  // namespace

void MakePPG12Fig27ABCDYieldSDCCReproduction(
    const char *input_root =
        "dataOutput/ppg12PhotonYield/ppg12_photon_yield_v1_data_20260620/reference_roots/fig27_abcd_yield/data_histo_bdt_nom.root",
    const char *output_dir =
        "dataOutput/ppg12PhotonYield/ppg12_photon_yield_v1_data_20260620/shower_shape_reference_validation/fig27_abcd_yield")
{
  gSystem->mkdir(output_dir, kTRUE);
  init_plot();
  gStyle->SetOptStat(0);

  TFile file(input_root, "READ");
  if (file.IsZombie() || file.TestBit(TFile::kRecovered))
  {
    Error("MakePPG12Fig27ABCDYieldSDCCReproduction",
          "Input ROOT is zombie or recovered: %s", input_root);
    return;
  }

  TH1D *h_tight_iso_cluster = requireHist(file, "h_tight_iso_cluster_0");
  TH1D *h_tight_noniso_cluster = requireHist(file, "h_tight_noniso_cluster_0");
  TH1D *h_nontight_iso_cluster = requireHist(file, "h_nontight_iso_cluster_0");
  TH1D *h_nontight_noniso_cluster = requireHist(file, "h_nontight_noniso_cluster_0");
  if (!h_tight_iso_cluster || !h_tight_noniso_cluster || !h_nontight_iso_cluster ||
      !h_nontight_noniso_cluster)
  {
    return;
  }

  styleABCD(h_tight_iso_cluster, kBlack);
  styleABCD(h_tight_noniso_cluster, kRed);
  styleABCD(h_nontight_iso_cluster, kBlue);
  styleABCD(h_nontight_noniso_cluster, kMagenta);

  TCanvas *canvas = new TCanvas("c_ppg12_fig27_abcd_yield", "c_ppg12_fig27_abcd_yield", 600, 600);
  frame_et_rec->Draw("axis");
  frame_et_rec->GetXaxis()->SetRangeUser(10, 35);
  frame_et_rec->GetYaxis()->SetRangeUser(1, 5e5);

  h_tight_iso_cluster->Draw("same e");
  h_tight_noniso_cluster->Draw("same e");
  h_nontight_iso_cluster->Draw("same e");
  h_nontight_noniso_cluster->Draw("same e");

  myText(0.5, 0.9, 1, strleg1.c_str(), 0.04);
  myText(0.5, 0.85, 1, strleg2.c_str(), 0.04);
  myText(0.5, 0.80, 1, Form("Data   %s", strleg3.c_str()), 0.04);
  myText(0.18, 0.75, 1, "bdt_nom", 0.04);
  myMarkerLineText(0.55, 0.75, 0, kBlack, 0, kBlack, 1, "A: tight iso", 0.05, true);
  myMarkerLineText(0.55, 0.70, 0, kRed, 0, kRed, 1, "B: tight noniso", 0.05, true);
  myMarkerLineText(0.55, 0.65, 0, kBlue, 0, kBlue, 1, "C: nontight iso", 0.05, true);
  myMarkerLineText(0.55, 0.60, 0, kMagenta, 0, kMagenta, 1, "D: nontight noniso", 0.05, true);

  gPad->SetLogy();
  canvas->RedrawAxis();

  const TString output_png = Form("%s/ppg12_fig27_abcd_yield_bdt_nom_sdcc_reproduction.png", output_dir);
  canvas->SaveAs(output_png);

  printf("wrote=%s\n", output_png.Data());
  printf("source_root=%s\n", input_root);
  printf("source_code=ppg12codeGit/plotting/plot_sideband_selection.C(\"bdt_nom\")\n");
}
