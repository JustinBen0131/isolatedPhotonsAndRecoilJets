#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
#include <iostream>
#include <memory>
#include <string>

namespace
{
const std::string kTag = "jetMinPt5_7pi_8_vz60_isoR40_isSliding_baseVariant_preselectionReference_tightReference_nonTightReference";
const std::string kTrigger = "photon_10_plus_MBD_NS_geq_2_vtx_lt_150";
const std::string kMbdTrigger = "MBD_NS_geq_2_vtx_lt_150";
const std::string kDataFile =
    "InputFiles/auau25/RecoilJets_auau_ALL_" + kTag + ".root";
const std::string kSimFile =
    "dataOutput/combinedSimOnlyEMBEDDED/" + kTag +
    "/photonJet12and20merged_SIM/RecoilJets_embeddedPhoton12plus20_MERGED.root";
const std::string kSimPhoton12File =
    "InputFiles/simEmbedded/RecoilJets_embeddedPhoton12_ALL_" + kTag + ".root";
const std::string kSimPhoton20File =
    "InputFiles/simEmbedded/RecoilJets_embeddedPhoton20_ALL_" + kTag + ".root";
const std::string kOutDir =
    "dataOutput/auau/" + kTag + "/" + kTrigger;
const std::string kQaOutDir =
    "dataOutput/auau/" + kTag + "/vertexCentralityQA_rebin2_vz30";

TH1* GetHist(TFile* f, const std::string& dirName, const std::string& histName)
{
  if (!f) return nullptr;
  TDirectory* d = f->GetDirectory(dirName.c_str());
  if (!d)
  {
    std::cerr << "[ERROR] Missing directory " << dirName << " in " << f->GetName() << "\n";
    return nullptr;
  }
  TH1* h = dynamic_cast<TH1*>(d->Get(histName.c_str()));
  if (!h)
  {
    std::cerr << "[ERROR] Missing histogram " << histName << " in " << dirName
              << " from " << f->GetName() << "\n";
    return nullptr;
  }
  TH1* out = dynamic_cast<TH1*>(h->Clone((histName + "_" + dirName + "_clone").c_str()));
  if (out) out->SetDirectory(nullptr);
  return out;
}

void ShapeNormalize(TH1* h)
{
  if (!h) return;
  const double integral = h->Integral(1, h->GetNbinsX());
  if (integral > 0.0) h->Scale(1.0 / integral);
}

void StyleData(TH1* h)
{
  h->SetMarkerStyle(20);
  h->SetMarkerSize(0.75);
  h->SetMarkerColor(kBlack);
  h->SetLineColor(kBlack);
  h->SetLineWidth(1);
  h->SetStats(0);
}

void StyleSim(TH1* h)
{
  h->SetMarkerStyle(1);
  h->SetMarkerSize(0.0);
  h->SetMarkerColor(kRed + 1);
  h->SetLineColor(kRed + 1);
  h->SetLineWidth(2);
  h->SetStats(0);
}

double MaxBin(const TH1* h)
{
  if (!h) return 0.0;
  double m = 0.0;
  for (int i = 1; i <= h->GetNbinsX(); ++i) m = std::max(m, h->GetBinContent(i));
  return m;
}

void DrawPanel(TH1* hData, TH1* hSim,
               const std::string& title,
               const std::string& xTitle,
               const std::string& yTitle,
               double xMin, double xMax,
               bool drawInfoBlock,
               const std::string& selectionText = "Photon 10 GeV + MBD N&S #geq 2",
               const std::string& scopeText = "Centrality inclusive",
               const std::string& simLegend = "Photon+Jet (12+20) Embedded")
{
  hData->GetXaxis()->SetRangeUser(xMin, xMax);
  hSim->GetXaxis()->SetRangeUser(xMin, xMax);
  hData->SetTitle(title.c_str());
  hData->GetXaxis()->SetTitle(xTitle.c_str());
  hData->GetYaxis()->SetTitle(yTitle.c_str());
  hData->GetXaxis()->SetTitleSize(0.050);
  hData->GetYaxis()->SetTitleSize(0.047);
  hData->GetXaxis()->SetLabelSize(0.044);
  hData->GetYaxis()->SetLabelSize(0.044);
  hData->GetYaxis()->SetTitleOffset(1.35);
  hData->SetMaximum(1.25 * std::max(MaxBin(hData), MaxBin(hSim)));
  hData->SetMinimum(0.0);
  hData->Draw("E1");
  hSim->Draw("HIST SAME");
  hData->Draw("E1 SAME");

  TLegend* leg = new TLegend(0.18, 0.78, 0.72, 0.90);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.030);
  leg->AddEntry(hData, "Au+Au data", "pe");
  leg->AddEntry(hSim, simLegend.c_str(), "l");
  leg->Draw();

  if (!drawInfoBlock) return;
  TLatex tx;
  tx.SetNDC();
  tx.SetTextFont(42);
  tx.SetTextAlign(31);
  tx.SetTextSize(0.030);
  tx.DrawLatex(0.90, 0.66, selectionText.c_str());
  tx.DrawLatex(0.90, 0.60, "|v_{z}| < 60 cm");
  if (!scopeText.empty()) tx.DrawLatex(0.90, 0.54, scopeText.c_str());
  tx.SetTextSize(0.035);
  tx.DrawLatex(0.90, 0.44, "#bf{sPHENIX} #it{Internal}");
  tx.SetTextSize(0.030);
  tx.DrawLatex(0.90, 0.38, "Au+Au  #sqrt{s_{NN}} = 200 GeV");
}

void DrawTriggerTable(TFile* fData,
                      TFile* fSim,
                      const std::string& dataDir,
                      const std::string& outPng,
                      const std::string& selectionText,
                      const std::string& scopeText,
                      const std::string& simLegend,
                      int vertexRebin = 1,
                      double vertexXMin = -60.0,
                      double vertexXMax = 60.0)
{
  std::unique_ptr<TH1> hVzData(GetHist(fData, dataDir, "h_vertexZ"));
  std::unique_ptr<TH1> hCentData(GetHist(fData, dataDir, "h_centrality"));
  std::unique_ptr<TH1> hVzSim(GetHist(fSim, "SIM", "h_vertexZ"));
  std::unique_ptr<TH1> hCentSim(GetHist(fSim, "SIM", "h_centrality"));
  if (!hVzData || !hCentData || !hVzSim || !hCentSim) return;

  if (vertexRebin > 1)
  {
    hVzData->Rebin(vertexRebin);
    hVzSim->Rebin(vertexRebin);
  }

  for (TH1* h : {hVzData.get(), hCentData.get(), hVzSim.get(), hCentSim.get()})
  {
    ShapeNormalize(h);
  }
  StyleData(hVzData.get());
  StyleData(hCentData.get());
  StyleSim(hVzSim.get());
  StyleSim(hCentSim.get());

  TCanvas c(("c_" + dataDir + "_vertexZ_centrality_DATA_embedded12plus20_1x2").c_str(),
            ("c_" + dataDir + "_vertexZ_centrality_DATA_embedded12plus20_1x2").c_str(), 1600, 720);
  c.Divide(2, 1, 0.01, 0.01);

  c.cd(1);
  gPad->SetLeftMargin(0.16);
  gPad->SetRightMargin(0.035);
  gPad->SetTopMargin(0.10);
  gPad->SetBottomMargin(0.14);
  gPad->SetTicks(1, 1);
  DrawPanel(hVzData.get(), hVzSim.get(),
            "z_{vtx}: Au+Au data and embedded SIM",
            "v_{z} [cm]",
            "Shape normalized",
            vertexXMin, vertexXMax,
            false,
            selectionText,
            scopeText,
            simLegend);

  c.cd(2);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.035);
  gPad->SetTopMargin(0.10);
  gPad->SetBottomMargin(0.14);
  gPad->SetTicks(1, 1);
  DrawPanel(hCentData.get(), hCentSim.get(),
            "Centrality: Au+Au data and embedded SIM",
            "Centrality [%]",
            "Shape normalized",
            0.0, 100.0,
            true,
            selectionText,
            scopeText,
            simLegend);

  c.SaveAs(outPng.c_str());
  std::cout << "[DONE] Wrote " << outPng << "\n";
}
}

void MakeAuAuReferenceDataEmbeddedVertexCentralityTable()
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetTitleSize(0.040, "t");

  std::unique_ptr<TFile> fData(TFile::Open(kDataFile.c_str(), "READ"));
  std::unique_ptr<TFile> fSim(TFile::Open(kSimFile.c_str(), "READ"));
  if (!fData || fData->IsZombie() || !fSim || fSim->IsZombie())
  {
    std::cerr << "[ERROR] Could not open input files:\n"
              << "  data: " << kDataFile << "\n"
              << "  sim : " << kSimFile << "\n";
    return;
  }

  gSystem->mkdir(kOutDir.c_str(), true);

  DrawTriggerTable(fData.get(), fSim.get(), kMbdTrigger,
                   kOutDir + "/vertexZ_centrality_DATA_embedded12plus20_1x2_reference_MBD_NS_geq_2_vtx_lt_150.png",
                   "MBD N&S #geq 2, vtx < 150 cm",
                   "",
                   "Photon+Jet (12+20) Embedded");

  DrawTriggerTable(fData.get(), fSim.get(), kTrigger,
                   kOutDir + "/vertexZ_centrality_DATA_embedded12plus20_1x2_reference_photon10_MBD_NS_geq_2_vtx_lt_150.png",
                   "Photon 10 GeV + MBD N&S #geq 2",
                   "Photon-trigger selection",
                   "Photon+Jet (12+20) Embedded");

  DrawTriggerTable(fData.get(), fSim.get(), kMbdTrigger,
                   kOutDir + "/vertexZ_rebin2_vz30_centrality_DATA_embedded12plus20_1x2_reference_MBD_NS_geq_2_vtx_lt_150.png",
                   "MBD N&S #geq 2, vtx < 150 cm",
                   "",
                   "Photon+Jet (12+20) Embedded",
                   2, -30.0, 30.0);

  DrawTriggerTable(fData.get(), fSim.get(), kTrigger,
                   kOutDir + "/vertexZ_rebin2_vz30_centrality_DATA_embedded12plus20_1x2_reference_photon10_MBD_NS_geq_2_vtx_lt_150.png",
                   "Photon 10 GeV + MBD N&S #geq 2",
                   "Photon-trigger selection",
                   "Photon+Jet (12+20) Embedded",
                   2, -30.0, 30.0);

  std::unique_ptr<TFile> fSim12(TFile::Open(kSimPhoton12File.c_str(), "READ"));
  std::unique_ptr<TFile> fSim20(TFile::Open(kSimPhoton20File.c_str(), "READ"));
  if (!fSim12 || fSim12->IsZombie() || !fSim20 || fSim20->IsZombie())
  {
    std::cerr << "[ERROR] Could not open individual embedded inputs:\n"
              << "  photon12: " << kSimPhoton12File << "\n"
              << "  photon20: " << kSimPhoton20File << "\n";
    return;
  }

  const std::string combinedDir = kQaOutDir + "/combined12plus20";
  const std::string photon12Dir = kQaOutDir + "/photon12";
  const std::string photon20Dir = kQaOutDir + "/photon20";
  for (const std::string& dir : {combinedDir, photon12Dir, photon20Dir})
  {
    gSystem->mkdir((dir + "/MBD_NS_geq_2_vtx_lt_150").c_str(), true);
    gSystem->mkdir((dir + "/photon10_MBD_NS_geq_2_vtx_lt_150").c_str(), true);
  }

  DrawTriggerTable(fData.get(), fSim.get(), kMbdTrigger,
                   combinedDir + "/MBD_NS_geq_2_vtx_lt_150/vertexZ_rebin2_vz30_centrality_DATA_embedded12plus20.png",
                   "MBD N&S #geq 2, vtx < 150 cm",
                   "",
                   "Photon+Jet (12+20) Embedded",
                   2, -30.0, 30.0);

  DrawTriggerTable(fData.get(), fSim.get(), kTrigger,
                   combinedDir + "/photon10_MBD_NS_geq_2_vtx_lt_150/vertexZ_rebin2_vz30_centrality_DATA_embedded12plus20.png",
                   "Photon 10 GeV + MBD N&S #geq 2",
                   "Photon-trigger selection",
                   "Photon+Jet (12+20) Embedded",
                   2, -30.0, 30.0);

  DrawTriggerTable(fData.get(), fSim12.get(), kMbdTrigger,
                   photon12Dir + "/MBD_NS_geq_2_vtx_lt_150/vertexZ_rebin2_vz30_centrality_DATA_embeddedPhoton12.png",
                   "MBD N&S #geq 2, vtx < 150 cm",
                   "",
                   "Photon+Jet 12 Embedded",
                   2, -30.0, 30.0);

  DrawTriggerTable(fData.get(), fSim12.get(), kTrigger,
                   photon12Dir + "/photon10_MBD_NS_geq_2_vtx_lt_150/vertexZ_rebin2_vz30_centrality_DATA_embeddedPhoton12.png",
                   "Photon 10 GeV + MBD N&S #geq 2",
                   "Photon-trigger selection",
                   "Photon+Jet 12 Embedded",
                   2, -30.0, 30.0);

  DrawTriggerTable(fData.get(), fSim20.get(), kMbdTrigger,
                   photon20Dir + "/MBD_NS_geq_2_vtx_lt_150/vertexZ_rebin2_vz30_centrality_DATA_embeddedPhoton20.png",
                   "MBD N&S #geq 2, vtx < 150 cm",
                   "",
                   "Photon+Jet 20 Embedded",
                   2, -30.0, 30.0);

  DrawTriggerTable(fData.get(), fSim20.get(), kTrigger,
                   photon20Dir + "/photon10_MBD_NS_geq_2_vtx_lt_150/vertexZ_rebin2_vz30_centrality_DATA_embeddedPhoton20.png",
                   "Photon 10 GeV + MBD N&S #geq 2",
                   "Photon-trigger selection",
                   "Photon+Jet 20 Embedded",
                   2, -30.0, 30.0);
}
