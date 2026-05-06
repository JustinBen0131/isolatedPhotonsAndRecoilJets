#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>

#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TStyle.h>

namespace
{
std::string envOrDefault(const char* name, const std::string& fallback)
{
  if (const char* value = std::getenv(name))
  {
    if (*value) return value;
  }
  return fallback;
}

std::string joinPath(const std::string& a, const std::string& b)
{
  if (a.empty()) return b;
  if (a.back() == '/') return a + b;
  return a + "/" + b;
}

TH1* loadHist(TFile& f, const std::string& key)
{
  const std::string name = "h_maxEnergyClus_NewTriggerFilling_perRunCorrected_" + key;
  TDirectory* dir = dynamic_cast<TDirectory*>(f.Get(key.c_str()));
  if (!dir) return nullptr;
  TH1* src = dynamic_cast<TH1*>(dir->Get(name.c_str()));
  if (!src) return nullptr;
  TH1* h = dynamic_cast<TH1*>(src->Clone((name + "_plotClone").c_str()));
  if (h) h->SetDirectory(nullptr);
  return h;
}

TH1* normalizedClone(TH1* src, const char* name)
{
  if (!src) return nullptr;
  TH1* h = dynamic_cast<TH1*>(src->Clone(name));
  if (!h) return nullptr;
  h->SetDirectory(nullptr);
  const double integral = h->Integral(0, h->GetNbinsX() + 1);
  if (integral > 0.0) h->Scale(1.0 / integral);
  return h;
}

void styleHist(TH1* h, int color, int markerStyle)
{
  if (!h) return;
  h->SetLineColor(color);
  h->SetMarkerColor(color);
  h->SetLineWidth(3);
  h->SetMarkerStyle(markerStyle);
  h->SetMarkerSize(1.0);
  h->SetStats(false);
}

void saveCanvas(TCanvas& c, const std::string& path)
{
  c.SaveAs(path.c_str());
  std::cout << "[WROTE] " << path << std::endl;
}
}

int PlotScaledTriggerStudy(const char* inputPathArg = "", const char* outDirArg = "")
{
  const std::string cfg =
    "jetMinPt5_7pi_8_vz60_isoR40_isSliding_baseVariant_"
    "preselectionReference_tightReference_nonTightReference_scaledTriggerStudy";

  const std::string inputPath = inputPathArg && *inputPathArg
    ? inputPathArg
    : envOrDefault("RJ_SCALED_TRIGGER_INPUT",
                   "InputFiles/auau25/RecoilJets_auau_ALL_" + cfg + ".root");

  const std::string outDir = outDirArg && *outDirArg
    ? outDirArg
    : envOrDefault("RJ_SCALED_TRIGGER_PLOT_DIR",
                   "dataOutput/auau/" + cfg + "/scaledTrigQA");

  gSystem->mkdir(outDir.c_str(), true);
  gStyle->SetOptStat(0);

  TFile f(inputPath.c_str(), "READ");
  if (!f.IsOpen() || f.IsZombie())
  {
    std::cerr << "[ERROR] Could not open " << inputPath << std::endl;
    return 1;
  }

  std::unique_ptr<TH1> hBase(loadHist(f, "MBD_NS_geq_2_vtx_lt_150"));
  std::unique_ptr<TH1> hPho10(loadHist(f, "Photon_10"));
  std::unique_ptr<TH1> hPho12(loadHist(f, "Photon_12"));
  if (!hBase || !hPho10 || !hPho12)
  {
    std::cerr << "[ERROR] Missing one or more perRunCorrected scaled-trigger histograms in "
              << inputPath << std::endl;
    return 2;
  }

  std::unique_ptr<TH1> hBaseNorm(normalizedClone(hBase.get(), "hBaseNorm_scaledTrigQA"));
  std::unique_ptr<TH1> hPho10Norm(normalizedClone(hPho10.get(), "hPho10Norm_scaledTrigQA"));
  std::unique_ptr<TH1> hPho12Norm(normalizedClone(hPho12.get(), "hPho12Norm_scaledTrigQA"));

  styleHist(hBaseNorm.get(), kBlack, 20);
  styleHist(hPho10Norm.get(), kRed + 1, 21);
  styleHist(hPho12Norm.get(), kBlue + 1, 22);

  {
    TCanvas c("c_scaledTrigQA_groupOverlay", "c_scaledTrigQA_groupOverlay", 1100, 800);
    c.cd();
    gPad->SetLogy();
    gPad->SetTicks(1, 1);

    const double maxY = std::max(hBaseNorm->GetMaximum(),
                                 std::max(hPho10Norm->GetMaximum(), hPho12Norm->GetMaximum()));
    hBaseNorm->GetXaxis()->SetTitle("Max cluster energy [GeV]");
    hBaseNorm->GetYaxis()->SetTitle("Normalized counts");
    hBaseNorm->SetMinimum(1e-6);
    hBaseNorm->SetMaximum(std::max(1e-4, 1.5 * maxY));
    hBaseNorm->Draw("hist");
    hPho10Norm->Draw("hist same");
    hPho12Norm->Draw("hist same");

    TLegend leg(0.56, 0.68, 0.90, 0.88);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.034);
    leg.AddEntry(hBaseNorm.get(), "MBD N&S #geq 2, |v_{z}| < 150", "l");
    leg.AddEntry(hPho10Norm.get(), "Photon 10 (per-run corrected)", "l");
    leg.AddEntry(hPho12Norm.get(), "Photon 12 (per-run corrected)", "l");
    leg.Draw();

    TLatex t;
    t.SetNDC(true);
    t.SetTextSize(0.034);
    t.DrawLatex(0.14, 0.92, "Scaled-trigger QA");
    t.DrawLatex(0.14, 0.875, "Runs from scaledEffRuns_MBD_NS_geq_2_vtx_lt_150__Pho10_12.list");
    t.DrawLatex(0.14, 0.83, "Per-run histograms scaled by live/scaled before merge");

    saveCanvas(c, joinPath(outDir, "hMaxClusterEnergy_groupOverlay.png"));
  }

  std::unique_ptr<TH1> hPho10Eff(dynamic_cast<TH1*>(hPho10->Clone("hPho10Eff_scaledTrigQA")));
  std::unique_ptr<TH1> hPho12Eff(dynamic_cast<TH1*>(hPho12->Clone("hPho12Eff_scaledTrigQA")));
  if (hPho10Eff) hPho10Eff->SetDirectory(nullptr);
  if (hPho12Eff) hPho12Eff->SetDirectory(nullptr);
  hPho10Eff->Divide(hBase.get());
  hPho12Eff->Divide(hBase.get());
  styleHist(hPho10Eff.get(), kRed + 1, 21);
  styleHist(hPho12Eff.get(), kBlue + 1, 22);

  {
    TCanvas c("c_scaledTrigQA_groupTurnOnOverlay", "c_scaledTrigQA_groupTurnOnOverlay", 1100, 800);
    c.cd();
    gPad->SetTicks(1, 1);

    const double maxY = std::max(hPho10Eff->GetMaximum(), hPho12Eff->GetMaximum());
    hPho10Eff->GetXaxis()->SetTitle("Max cluster energy [GeV]");
    hPho10Eff->GetYaxis()->SetTitle("Efficiency");
    hPho10Eff->SetMinimum(0.0);
    hPho10Eff->SetMaximum(std::max(1.2, 1.15 * maxY));
    hPho10Eff->Draw("hist");
    hPho12Eff->Draw("hist same");

    TLegend leg(0.56, 0.72, 0.90, 0.88);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetTextSize(0.034);
    leg.AddEntry(hPho10Eff.get(), "Photon 10 / MBD", "l");
    leg.AddEntry(hPho12Eff.get(), "Photon 12 / MBD", "l");
    leg.Draw();

    TLatex t;
    t.SetNDC(true);
    t.SetTextSize(0.034);
    t.DrawLatex(0.14, 0.92, "Scaled-trigger QA turn-on");
    t.DrawLatex(0.14, 0.875, "Common selected-run subset only");

    saveCanvas(c, joinPath(outDir, "hMaxClusterEnergy_groupTurnOnOverlay.png"));
  }

  return 0;
}
