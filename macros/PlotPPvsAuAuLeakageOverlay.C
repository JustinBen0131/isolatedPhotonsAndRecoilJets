#include "AnalyzeRecoilJets.h"

#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH1F.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TStyle.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace
{
constexpr double kSigmaEmbeddedPhoton12To20_pb = 2598.12425;
constexpr double kSigmaEmbeddedPhoton20Plus_pb = 133.317866;

const std::string kPPTag =
    "jetMinPt5_7pi_8_vz60_isoR40_fixedIso2GeV_preselectionReference_tightReference_nonTightReference";
const std::string kAATag =
    "jetMinPt5_7pi_8_vz60_isoR40_isSliding_baseVariant_preselectionReference_tightReference_nonTightReference";
const std::string kComboPP = "photonJet5and10and20merged_SIM";
const std::string kComboAA = "photonJet12and20merged_SIM";
const std::string kMergedPPName = "RecoilJets_photonjet5plus10plus20_MERGED.root";
const std::string kMergedAAName = "RecoilJets_embeddedPhoton12plus20_MERGED.root";

struct Series
{
  std::string label;
  std::string centSuffix;
  int color = kBlack;
  int marker = 20;
  bool isPP = false;
};

struct LeakSpec
{
  int regionBin = 2;
  const char* yTitle = "";
  const char* title = "";
  const char* outName = "";
};

std::string PPInput(const std::string& sample)
{
  return ARJ::kInputBase + "/simPhotonJet/RecoilJets_" + sample + "_ALL_" + kPPTag + ".root";
}

std::string AAInput(const std::string& sample)
{
  return ARJ::kInputBase + "/simEmbedded/RecoilJets_" + sample + "_ALL_" + kAATag + ".root";
}

std::string PPMergedPath()
{
  return ARJ::kOutputBase + "/combinedSimOnly/" + kPPTag + "/" + kComboPP + "/" + kMergedPPName;
}

std::string AAMergedPath()
{
  return ARJ::kInputBase + "/simEmbedded/merged/" + kAATag + "/" + kComboAA + "/" + kMergedAAName;
}

std::string OutputDir()
{
  return ARJ::kOutputBase + "/combinedSimOnlyEMBEDDED/" + kAATag + "/" + kComboAA;
}

double RatioErr(const double num, const double den)
{
  if (num <= 0.0 || den <= 0.0) return 0.0;
  const double r = num / den;
  return r * std::sqrt(1.0 / num + 1.0 / den);
}

bool RebuildInputs()
{
  const bool okPP = ARJ::BuildMergedSIMFile_PhotonSlices(
      {PPInput("photonjet5"), PPInput("photonjet10"), PPInput("photonjet20")},
      {ARJ::kSigmaPhoton5_pb, ARJ::kSigmaPhoton10_pb, ARJ::kSigmaPhoton20_pb},
      PPMergedPath(),
      ARJ::kDirSIM,
      {"photonJet5", "photonJet10", "photonJet20"});

  const bool okAA = ARJ::BuildMergedSIMFile_PhotonSlices(
      {AAInput("embeddedPhoton12"), AAInput("embeddedPhoton20")},
      {kSigmaEmbeddedPhoton12To20_pb, kSigmaEmbeddedPhoton20Plus_pb},
      AAMergedPath(),
      ARJ::kDirSIM,
      {"embeddedPhoton12", "embeddedPhoton20"});

  return okPP && okAA;
}

std::unique_ptr<TGraphErrors> BuildLeakageGraph(TDirectory* dir,
                                                const Series& series,
                                                const LeakSpec& spec)
{
  if (!dir) return nullptr;

  std::vector<double> x;
  std::vector<double> ex;
  std::vector<double> y;
  std::vector<double> ey;
  x.reserve(ARJ::kNPtBins);
  ex.reserve(ARJ::kNPtBins);
  y.reserve(ARJ::kNPtBins);
  ey.reserve(ARJ::kNPtBins);

  for (int i = 0; i < ARJ::kNPtBins; ++i)
  {
    const ARJ::PtBin& b = ARJ::PtBins()[i];
    if (b.lo < 10) continue;

    const std::string hname = "h_sigABCD_MC" + b.suffix + series.centSuffix;
    TH1* h = dynamic_cast<TH1*>(dir->Get(hname.c_str()));
    if (!h)
    {
      std::cerr << "[WARN] Missing " << hname << " for " << series.label << "\n";
      continue;
    }

    const double a = h->GetBinContent(1);
    const double n = h->GetBinContent(spec.regionBin);
    if (a <= 0.0) continue;

    x.push_back(0.5 * (ARJ::kPtEdges[(std::size_t)i] + ARJ::kPtEdges[(std::size_t)i + 1]));
    ex.push_back(0.0);
    y.push_back(n / a);
    ey.push_back(RatioErr(n, a));
  }

  if (x.empty()) return nullptr;

  std::unique_ptr<TGraphErrors> g(new TGraphErrors((int)x.size(), x.data(), y.data(), ex.data(), ey.data()));
  g->SetName(("g_" + series.label + "_" + std::to_string(spec.regionBin)).c_str());
  g->SetLineColor(series.color);
  g->SetMarkerColor(series.color);
  g->SetMarkerStyle(series.isPP ? 24 : 20);
  g->SetMarkerSize(1.05);
  g->SetLineWidth(2);
  return g;
}

double GraphYMax(const std::vector<std::unique_ptr<TGraphErrors>>& graphs)
{
  double ymax = 0.0;
  for (const auto& g : graphs)
  {
    if (!g) continue;
    for (int i = 0; i < g->GetN(); ++i)
    {
      double gx = 0.0;
      double gy = 0.0;
      g->GetPoint(i, gx, gy);
      ymax = std::max(ymax, gy + g->GetErrorY(i));
    }
  }
  return ymax;
}

void DrawOne(TDirectory* ppDir, TDirectory* aaDir, const LeakSpec& spec)
{
  const std::vector<Series> series = {
      {"pp", "", kRed + 1, 24, true},
      {"0-20%", "_cent_0_20", kBlack, 20, false},
      {"20-50%", "_cent_20_50", kBlue + 1, 20, false},
      {"50-80%", "_cent_50_80", kOrange + 7, 20, false},
  };

  std::vector<std::unique_ptr<TGraphErrors>> graphs;
  graphs.reserve(series.size());
  for (const Series& s : series)
  {
    graphs.emplace_back(BuildLeakageGraph(s.isPP ? ppDir : aaDir, s, spec));
  }

  const double ymax = std::max(0.08, GraphYMax(graphs));
  TCanvas c(("c_" + std::string(spec.outName)).c_str(), spec.title, 1100, 850);
  c.SetLeftMargin(0.14);
  c.SetRightMargin(0.04);
  c.SetTopMargin(0.12);
  c.SetBottomMargin(0.14);
  c.SetTicks(1, 1);

  TH1F frame(("hframe_" + std::string(spec.outName)).c_str(), "", 100, 10.0, 35.0);
  frame.SetDirectory(nullptr);
  frame.SetStats(0);
  frame.SetMinimum(0.0);
  const double yScale = (spec.regionBin == 2) ? 1.22 : 1.35;
  frame.SetMaximum(yScale * ymax);
  frame.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
  frame.GetYaxis()->SetTitle(spec.yTitle);
  frame.GetXaxis()->SetTitleSize(0.048);
  frame.GetYaxis()->SetTitleSize(0.048);
  frame.GetXaxis()->SetLabelSize(0.040);
  frame.GetYaxis()->SetLabelSize(0.040);
  frame.Draw();

  for (const auto& g : graphs)
  {
    if (g) g->Draw("P SAME");
  }

  const bool isBLeakage = (spec.regionBin == 2);
  TLegend leg(isBLeakage ? 0.56 : 0.56,
              isBLeakage ? 0.20 : 0.57,
              isBLeakage ? 0.92 : 0.92,
              isBLeakage ? 0.38 : 0.75);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextFont(42);
  leg.SetTextSize(0.032);
  leg.SetNColumns(2);
  for (std::size_t i = 0; i < graphs.size(); ++i)
  {
    if (graphs[i]) leg.AddEntry(graphs[i].get(), series[i].label.c_str(), "pe");
  }
  leg.Draw();

  TLatex tx;
  tx.SetNDC();
  tx.SetTextFont(42);
  tx.SetTextAlign(22);
  tx.SetTextSize(0.030);
  tx.DrawLatex(0.50, 0.895, spec.title);

  tx.SetTextAlign(13);
  tx.SetTextSize(0.026);
  tx.DrawLatex(0.17, 0.84, "Reference preselection, reference tight/non-tight");
  tx.DrawLatex(0.17, 0.80, "Au+Au: embedded #gamma+jet 12+20, sliding iso");
  tx.DrawLatex(0.17, 0.76, "pp: #gamma+jet 5+10+20, fixed E_{T}^{iso} < 2 GeV");

  tx.SetTextAlign(31);
  tx.SetTextSize(0.032);
  tx.DrawLatex(0.92, 0.84, "#bf{sPHENIX} #it{Internal}");
  tx.DrawLatex(0.92, 0.80, "Pythia #sqrt{s} = 200 GeV");

  const std::string outPath = OutputDir() + "/" + spec.outName;
  c.SaveAs(outPath.c_str());
  std::cout << "[WROTE] " << outPath << "\n";
}
}

void PlotPPvsAuAuLeakageOverlay()
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gSystem->mkdir(OutputDir().c_str(), true);

  if (!RebuildInputs())
  {
    std::cerr << "[ERROR] Failed to build one or both merged inputs.\n";
    return;
  }

  std::unique_ptr<TFile> fPP(TFile::Open(PPMergedPath().c_str(), "READ"));
  std::unique_ptr<TFile> fAA(TFile::Open(AAMergedPath().c_str(), "READ"));
  if (!fPP || fPP->IsZombie() || !fAA || fAA->IsZombie())
  {
    std::cerr << "[ERROR] Cannot open merged pp or AuAu input.\n";
    return;
  }

  TDirectory* ppDir = fPP->GetDirectory(ARJ::kDirSIM.c_str());
  TDirectory* aaDir = fAA->GetDirectory(ARJ::kDirSIM.c_str());
  if (!ppDir || !aaDir)
  {
    std::cerr << "[ERROR] Missing SIM top directory in merged inputs.\n";
    return;
  }

  const std::vector<LeakSpec> specs = {
      {2, "f_{B} = B_{sig}/A_{sig}", "B leakage, pp vs Au+Au centrality, #gamma+jet simulation", "leakageFactor_fB_pp_AuAuCentrality_overlay.png"},
      {3, "f_{C} = C_{sig}/A_{sig}", "C leakage, pp vs Au+Au centrality, #gamma+jet simulation", "leakageFactor_fC_pp_AuAuCentrality_overlay.png"},
      {4, "f_{D} = D_{sig}/A_{sig}", "D leakage, pp vs Au+Au centrality, #gamma+jet simulation", "leakageFactor_fD_pp_AuAuCentrality_overlay.png"},
  };

  for (const LeakSpec& spec : specs)
  {
    DrawOne(ppDir, aaDir, spec);
  }
}
