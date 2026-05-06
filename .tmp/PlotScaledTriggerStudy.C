#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TGaxis.h"
#include "TH1.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"

namespace
{
std::string JoinPath(const std::string& a, const std::string& b)
{
  if (a.empty()) return b;
  if (a.back() == '/') return a + b;
  return a + "/" + b;
}

TH1* LoadHist(TFile& f, const std::string& key)
{
  const std::string prefix = "h_maxEnergyClus_NewTriggerFilling_perRunCorrected_";
  TDirectory* dir = dynamic_cast<TDirectory*>(f.Get(key.c_str()));
  if (!dir)
  {
    std::cerr << "[ERROR] Missing directory: " << key << "\n";
    return nullptr;
  }
  TH1* h = dynamic_cast<TH1*>(dir->Get((prefix + key).c_str()));
  if (!h)
  {
    std::cerr << "[ERROR] Missing histogram: " << key << "/" << prefix << key << "\n";
    return nullptr;
  }
  TH1* out = dynamic_cast<TH1*>(h->Clone((key + "_scaledTriggerStudy_clone").c_str()));
  if (out) out->SetDirectory(nullptr);
  return out;
}

TH1* NormalizedClone(TH1* h, const char* name)
{
  TH1* out = dynamic_cast<TH1*>(h->Clone(name));
  if (!out) return nullptr;
  out->SetDirectory(nullptr);
  const double integral = out->Integral(0, out->GetNbinsX() + 1);
  if (integral > 0.0) out->Scale(1.0 / integral);
  return out;
}

void StyleHist(TH1* h, int color, int marker)
{
  if (!h) return;
  h->SetStats(false);
  h->SetTitle("");
  h->SetLineColor(color);
  h->SetMarkerColor(color);
  h->SetLineWidth(3);
  h->SetMarkerStyle(marker);
  h->SetMarkerSize(0.9);
}

void StyleAxes(TH1* h, const char* xTitle, const char* yTitle)
{
  h->GetXaxis()->SetTitle(xTitle);
  h->GetYaxis()->SetTitle(yTitle);
  h->GetXaxis()->SetTitleSize(0.046);
  h->GetYaxis()->SetTitleSize(0.046);
  h->GetXaxis()->SetLabelSize(0.039);
  h->GetYaxis()->SetLabelSize(0.039);
  h->GetYaxis()->SetTitleOffset(1.12);
  h->GetXaxis()->SetTitleOffset(1.02);
}

void DrawLabelBlock(const char* title, const char* subtitle)
{
  TLatex t;
  t.SetNDC(true);
  t.SetTextFont(42);
  t.SetTextSize(0.036);
  t.DrawLatex(0.135, 0.925, "sPHENIX Internal");
  t.SetTextSize(0.033);
  t.DrawLatex(0.135, 0.878, title);
  t.SetTextSize(0.030);
  t.DrawLatex(0.135, 0.835, subtitle);
}
}  // namespace

void PlotScaledTriggerStudy()
{
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);
  TGaxis::SetMaxDigits(3);

  const std::string input =
      "InputFiles/auau25/RecoilJets_auau_ALL_jetMinPt5_7pi_8_vz60_isoR40_isSliding_baseVariant_preselectionReference_tightReference_nonTightReference_scaledTriggerStudy.root";
  const std::string outDir =
      "dataOutput/auau/jetMinPt5_7pi_8_vz60_isoR40_isSliding_baseVariant_preselectionReference_tightReference_nonTightReference_scaledTriggerStudy/scaledTrigQA";
  gSystem->mkdir(outDir.c_str(), true);

  TFile f(input.c_str(), "READ");
  if (f.IsZombie())
  {
    std::cerr << "[FATAL] Could not open " << input << "\n";
    return;
  }

  TH1* hBaseRaw = LoadHist(f, "MBD_NS_geq_2_vtx_lt_150");
  TH1* hPho10Raw = LoadHist(f, "Photon_10");
  TH1* hPho12Raw = LoadHist(f, "Photon_12");
  if (!hBaseRaw || !hPho10Raw || !hPho12Raw) return;

  std::cout << "RECOILJETS_SCALED_TRIGGER_PLOT_INPUTS\n";
  std::cout << "  base_integral=" << hBaseRaw->Integral(0, hBaseRaw->GetNbinsX() + 1) << "\n";
  std::cout << "  photon10_integral=" << hPho10Raw->Integral(0, hPho10Raw->GetNbinsX() + 1) << "\n";
  std::cout << "  photon12_integral=" << hPho12Raw->Integral(0, hPho12Raw->GetNbinsX() + 1) << "\n";

  TH1* hBaseCounts = dynamic_cast<TH1*>(hBaseRaw->Clone("hBaseCounts_scaledTriggerStudy"));
  TH1* hPho10Counts = dynamic_cast<TH1*>(hPho10Raw->Clone("hPho10Counts_scaledTriggerStudy"));
  TH1* hPho12Counts = dynamic_cast<TH1*>(hPho12Raw->Clone("hPho12Counts_scaledTriggerStudy"));
  hBaseCounts->SetDirectory(nullptr);
  hPho10Counts->SetDirectory(nullptr);
  hPho12Counts->SetDirectory(nullptr);

  StyleHist(hBaseCounts, kBlack, 20);
  StyleHist(hPho10Counts, kRed + 1, 21);
  StyleHist(hPho12Counts, kBlue + 1, 22);

  {
    TCanvas c("c_scaledTriggerStudy_maxCluster", "", 1200, 900);
    c.SetLeftMargin(0.12);
    c.SetRightMargin(0.04);
    c.SetTopMargin(0.08);
    c.SetBottomMargin(0.12);
    c.SetTicks(1, 1);
    c.SetLogy();

    const double maxY = std::max({hBaseCounts->GetMaximum(), hPho10Counts->GetMaximum(), hPho12Counts->GetMaximum()});
    StyleAxes(hBaseCounts, "Max cluster energy [GeV]", "Scaled Counts");
    hBaseCounts->SetMinimum(0.8);
    hBaseCounts->SetMaximum(std::max(10.0, 1.8 * maxY));
    hBaseCounts->GetXaxis()->SetRangeUser(0.0, 32.0);
    hBaseCounts->Draw("hist");
    hPho10Counts->Draw("hist same");
    hPho12Counts->Draw("hist same");

    TLegend leg(0.58, 0.68, 0.90, 0.88);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetTextFont(42);
    leg.SetTextSize(0.031);
    leg.AddEntry(hBaseCounts, "MBD N&S #geq 2, |v_{z}| < 150 cm", "l");
    leg.AddEntry(hPho10Counts, "Photon 10, corrected", "l");
    leg.AddEntry(hPho12Counts, "Photon 12, corrected", "l");
    leg.Draw();

    DrawLabelBlock("Run3 Au+Au scaled-trigger study", "Max cluster energy, common selected-run subset");
    const std::string out = JoinPath(outDir, "scaledTriggerStudy_maxClusterEnergyDistribution.png");
    c.SaveAs(out.c_str());
    std::cout << "[WROTE] " << out << "\n";
  }

  TH1* hPho10Eff = dynamic_cast<TH1*>(hPho10Raw->Clone("hPho10Eff_scaledTriggerStudy"));
  TH1* hPho12Eff = dynamic_cast<TH1*>(hPho12Raw->Clone("hPho12Eff_scaledTriggerStudy"));
  hPho10Eff->SetDirectory(nullptr);
  hPho12Eff->SetDirectory(nullptr);
  hPho10Eff->Divide(hBaseRaw);
  hPho12Eff->Divide(hBaseRaw);
  StyleHist(hPho10Eff, kRed + 1, 21);
  StyleHist(hPho12Eff, kBlue + 1, 22);

  {
    TCanvas c("c_scaledTriggerStudy_efficiency", "", 1200, 900);
    c.SetLeftMargin(0.12);
    c.SetRightMargin(0.04);
    c.SetTopMargin(0.08);
    c.SetBottomMargin(0.12);
    c.SetTicks(1, 1);

    StyleAxes(hPho10Eff, "Max cluster energy [GeV]", "Trigger / MBD efficiency");
    hPho10Eff->SetMinimum(0.0);
    hPho10Eff->SetMaximum(1.25);
    hPho10Eff->GetXaxis()->SetRangeUser(0.0, 32.0);
    hPho10Eff->Draw("hist");
    hPho12Eff->Draw("hist same");

    TLine line90(0.0, 0.9, 32.0, 0.9);
    line90.SetLineStyle(2);
    line90.SetLineWidth(2);
    line90.SetLineColor(kGray + 2);
    line90.Draw("same");

    TLegend leg(0.58, 0.72, 0.90, 0.88);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetTextFont(42);
    leg.SetTextSize(0.031);
    leg.AddEntry(hPho10Eff, "Photon 10 / MBD", "l");
    leg.AddEntry(hPho12Eff, "Photon 12 / MBD", "l");
    leg.AddEntry(&line90, "90% reference", "l");
    leg.Draw();

    DrawLabelBlock("Run3 Au+Au scaled-trigger study", "Per-run live/scaled corrected turn-on");
    const std::string out = JoinPath(outDir, "scaledTriggerStudy_triggerEfficiency.png");
    c.SaveAs(out.c_str());
    std::cout << "[WROTE] " << out << "\n";
  }

  delete hBaseRaw;
  delete hPho10Raw;
  delete hPho12Raw;
  delete hBaseCounts;
  delete hPho10Counts;
  delete hPho12Counts;
  delete hPho10Eff;
  delete hPho12Eff;
}
