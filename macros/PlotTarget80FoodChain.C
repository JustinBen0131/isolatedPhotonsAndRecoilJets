#include "sPhenixStyle.C"

#include <TCanvas.h>
#include <TColor.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2.h>
#include <TH2F.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace
{
const std::string kTargetBase =
    "dataOutput/target80_all_available/bdt_target80_gated_20260512_001012/"
    "analysis_config_expanded_5to40_target80";
const std::string kOutDir = kTargetBase + "/food_chain_plots";
const std::string kSignalRel =
    "photonJet12and20merged_SIM/RecoilJets_embeddedPhoton12plus20_MERGED.root";
const std::string kPPFile =
    "dataOutput/combinedSimOnly/"
    "jetMinPt5_7pi_8_vz60_isoR40_fixedIso2GeV_preselectionReference_tightReference_nonTightReference/"
    "photonJet5and10and20merged_SIM/RecoilJets_photonjet5plus10plus20_MERGED.root";

struct CentSpec
{
  std::string suffix;
  std::string title;
};

struct PtBin
{
  int lo = 0;
  int hi = 0;
};

struct ModelSpec
{
  std::string label;
  std::string cfg;
};

const std::vector<CentSpec> kCents = {
    {"0_20", "0-20% central"},
    {"20_50", "20-50% mid-central"},
    {"50_80", "50-80% peripheral"}};

const std::vector<PtBin> kPtBins = {
    {14, 16}, {16, 18}, {18, 20}, {20, 22}, {22, 24}, {24, 26}, {26, 35}};

const std::vector<ModelSpec> kModels = {
    {"Box-cuts", "preselectionNewPPG12_tightReference_nonTightReference_baseVariant"},
    {"Centrality-input BDT target-80", "preselectionNewPPG12_tightAuAuCentInputBDT_nonTightAuAuBDTComplement_baseVariant"}};

TDirectory* simDir(TFile* f)
{
  return f ? dynamic_cast<TDirectory*>(f->Get("SIM")) : nullptr;
}

TH1* getHist(TDirectory* d, const std::string& name)
{
  TH1* h = d ? dynamic_cast<TH1*>(d->Get(name.c_str())) : nullptr;
  if (h) h->SetStats(false);
  return h;
}

TH2* getHist2(TDirectory* d, const std::string& name)
{
  TH2* h = d ? dynamic_cast<TH2*>(d->Get(name.c_str())) : nullptr;
  if (h) h->SetStats(false);
  return h;
}

double integralAndError(TH1* h, double& err)
{
  err = 0.0;
  if (!h) return 0.0;
  h->SetStats(false);
  return h->IntegralAndError(0, h->GetNbinsX() + 1, err);
}

std::string signalPath(const ModelSpec& model)
{
  return kTargetBase + "/simembedded/" + model.cfg + "/" + kSignalRel;
}

std::string ptSuffix(const PtBin& pt)
{
  return "_pT_" + std::to_string(pt.lo) + "_" + std::to_string(pt.hi);
}

std::string centSuffix(const CentSpec& cent)
{
  return "_cent_" + cent.suffix;
}

void drawHeader(double x, double y, const std::string& line2)
{
  TLatex tx;
  tx.SetNDC();
  tx.SetTextFont(42);
  tx.SetTextAlign(13);
  tx.SetTextSize(0.036);
  tx.DrawLatex(x, y, "#it{#bf{sPHENIX}} Internal");
  tx.SetTextSize(0.028);
  tx.DrawLatex(x, y - 0.055, line2.c_str());
}

double truthSigAbcdCount(TDirectory* d, const PtBin& pt, const CentSpec& cent, int regionBin)
{
  TH1* h = getHist(d, "h_sigABCD_MC_isoR30_isSliding" + ptSuffix(pt) + centSuffix(cent));
  if (!h) return 0.0;
  return h->GetBinContent(regionBin);
}

std::array<double, 4> truthSigAbcdFractions(TDirectory* d, const CentSpec& cent)
{
  std::array<double, 4> counts = {0.0, 0.0, 0.0, 0.0};
  for (const auto& pt : kPtBins)
  {
    for (int ir = 0; ir < 4; ++ir)
    {
      counts[ir] += truthSigAbcdCount(d, pt, cent, ir + 1);
    }
  }
  const double sum = counts[0] + counts[1] + counts[2] + counts[3];
  if (sum > 0.0)
  {
    for (double& v : counts) v /= sum;
  }
  return counts;
}

double abcdCandidateYield(TDirectory* d, const PtBin& pt, const CentSpec& cent, char region)
{
  TH1* h = getHist(d, std::string("h_Eiso_ABCD_") + region + "_isoR30_isSliding" + ptSuffix(pt) + centSuffix(cent));
  double err = 0.0;
  return integralAndError(h, err);
}

std::array<double, 4> candidateAbcdFractions(TDirectory* d, const CentSpec& cent)
{
  std::array<double, 4> counts = {0.0, 0.0, 0.0, 0.0};
  const std::array<char, 4> reg = {'A', 'B', 'C', 'D'};
  for (const auto& pt : kPtBins)
  {
    for (int ir = 0; ir < 4; ++ir)
    {
      counts[ir] += abcdCandidateYield(d, pt, cent, reg[ir]);
    }
  }
  const double sum = counts[0] + counts[1] + counts[2] + counts[3];
  if (sum > 0.0)
  {
    for (double& v : counts) v /= sum;
  }
  return counts;
}

std::unique_ptr<TH1D> projectXJ(TH2* h2, const std::string& name, double ptMin, double ptMax)
{
  if (!h2) return nullptr;
  const int ixLo = h2->GetXaxis()->FindBin(ptMin + 1e-6);
  const int ixHi = h2->GetXaxis()->FindBin(ptMax - 1e-6);
  std::unique_ptr<TH1D> h(h2->ProjectionY(name.c_str(), ixLo, ixHi, "e"));
  if (!h) return nullptr;
  h->SetDirectory(nullptr);
  h->SetStats(false);
  return h;
}

void normalize(TH1* h)
{
  if (!h) return;
  double err = 0.0;
  const double sum = h->IntegralAndError(0, h->GetNbinsX() + 1, err);
  if (sum > 0.0) h->Scale(1.0 / sum);
}

void setHistStyle(TH1* h, int color, int marker)
{
  h->SetLineColor(color);
  h->SetMarkerColor(color);
  h->SetMarkerStyle(marker);
  h->SetMarkerSize(1.1);
  h->SetLineWidth(2);
}

void drawAbcdHeatmap(bool truthSignal)
{
  std::vector<std::unique_ptr<TFile>> files;
  std::vector<TDirectory*> dirs;
  for (const auto& model : kModels)
  {
    std::unique_ptr<TFile> f(TFile::Open(signalPath(model).c_str(), "READ"));
    if (!f || f->IsZombie())
    {
      std::cerr << "[WARN] missing " << signalPath(model) << "\n";
      continue;
    }
    TDirectory* d = simDir(f.get());
    dirs.push_back(d);
    files.push_back(std::move(f));
  }

  TCanvas c(truthSignal ? "c_abcd_truth" : "c_abcd_candidates",
            truthSignal ? "c_abcd_truth" : "c_abcd_candidates", 2250, 760);
  c.Divide(3, 1, 0.012, 0.0);
  std::vector<std::unique_ptr<TH2F>> heatmaps;
  std::ofstream csv(kOutDir + (truthSignal ? "/truth_signal_abcd_fractions.csv" : "/candidate_abcd_fractions.csv"));
  csv << "source,centrality,region,fraction\n";

  const int blue = TColor::GetColor("#2166AC");
  const int red = TColor::GetColor("#B2182B");
  const std::array<std::string, 4> regions = {"A", "B", "C", "D"};

  for (size_t ic = 0; ic < kCents.size(); ++ic)
  {
    c.cd(ic + 1);
    gPad->SetTicks(1, 1);
    gPad->SetLeftMargin(0.18);
    gPad->SetRightMargin(ic == 2 ? 0.05 : 0.02);
    gPad->SetTopMargin(0.26);
    gPad->SetBottomMargin(0.13);

    auto h = std::make_unique<TH2F>(("h_abcd_" + kCents[ic].suffix + (truthSignal ? "_truth" : "_cand")).c_str(),
                                    "", 4, 0.5, 4.5, 2, 0.5, 2.5);
    h->SetDirectory(nullptr);
    h->SetStats(false);
    h->GetXaxis()->SetTitle("ABCD region");
    h->GetYaxis()->SetTitle("");
    h->GetZaxis()->SetRangeUser(0.0, 1.0);
    for (int ir = 0; ir < 4; ++ir) h->GetXaxis()->SetBinLabel(ir + 1, regions[ir].c_str());
    h->GetYaxis()->SetBinLabel(1, ic == 0 ? "Box-cuts" : "");
    h->GetYaxis()->SetBinLabel(2, ic == 0 ? "BDT target-80" : "");
    h->GetXaxis()->SetLabelSize(0.060);
    h->GetYaxis()->SetLabelSize(0.048);
    h->GetXaxis()->SetTitleSize(0.050);

    for (size_t im = 0; im < kModels.size() && im < dirs.size(); ++im)
    {
      const auto vals = truthSignal ? truthSigAbcdFractions(dirs[im], kCents[ic])
                                    : candidateAbcdFractions(dirs[im], kCents[ic]);
      for (int ir = 0; ir < 4; ++ir)
      {
        h->SetBinContent(ir + 1, im + 1, vals[ir]);
        csv << kModels[im].label << "," << kCents[ic].suffix << "," << regions[ir] << ","
            << std::setprecision(10) << vals[ir] << "\n";
      }
    }

    h->Draw("COL");
    TLatex text;
    text.SetTextFont(42);
    text.SetTextAlign(22);
    for (int iy = 1; iy <= 2; ++iy)
    {
      for (int ix = 1; ix <= 4; ++ix)
      {
        const double val = h->GetBinContent(ix, iy);
        text.SetTextSize(0.054);
        text.SetTextColor(val > 0.55 ? kWhite : kBlack);
        text.DrawLatex(ix, iy, Form("%.0f%%", 100.0 * val));
      }
    }

    TLatex title;
    title.SetNDC();
    title.SetTextFont(42);
    title.SetTextAlign(22);
    title.SetTextSize(0.052);
    title.DrawLatex(0.52, 0.955, kCents[ic].title.c_str());
    if (ic == 0)
    {
      drawHeader(0.22, 0.85, truthSignal
                               ? "Truth-matched photon candidates, 15 #leq E_{T}^{#gamma} < 35 GeV"
                               : "All reconstructed candidates, 15 #leq E_{T}^{#gamma} < 35 GeV");
    }
    heatmaps.push_back(std::move(h));
  }

  const std::string out = kOutDir + (truthSignal ? "/abcd_truth_signal_region_fractions_1x3.png"
                                                 : "/abcd_candidate_region_fractions_1x3.png");
  c.SaveAs(out.c_str());
}

void drawRecoXJOverlay()
{
  const ModelSpec model = kModels[1];
  std::unique_ptr<TFile> fAA(TFile::Open(signalPath(model).c_str(), "READ"));
  std::unique_ptr<TFile> fPP(TFile::Open(kPPFile.c_str(), "READ"));
  TDirectory* dAA = simDir(fAA.get());
  TDirectory* dPP = simDir(fPP.get());
  if (!dAA || !dPP)
  {
    std::cerr << "[WARN] missing directories for xJ overlay\n";
    return;
  }

  const double ptMin = 15.0;
  const double ptMax = 35.0;
  auto hCent = projectXJ(getHist2(dAA, "h2_unfoldReco_pTgamma_xJ_incl_r03_jetPt7_dphiPi2_isoR30_isSliding_cent_0_20"),
                         "h_xj_cent", ptMin, ptMax);
  auto hPeri = projectXJ(getHist2(dAA, "h2_unfoldReco_pTgamma_xJ_incl_r03_jetPt7_dphiPi2_isoR30_isSliding_cent_50_80"),
                         "h_xj_peri", ptMin, ptMax);
  auto hPP = projectXJ(getHist2(dPP, "h2_unfoldReco_pTgamma_xJ_incl_r03"),
                       "h_xj_pp", ptMin, ptMax);
  if (!hCent || !hPeri || !hPP)
  {
    std::cerr << "[WARN] missing xJ histograms for overlay\n";
    return;
  }
  normalize(hCent.get());
  normalize(hPeri.get());
  normalize(hPP.get());
  setHistStyle(hCent.get(), TColor::GetColor("#D55E00"), 20);
  setHistStyle(hPeri.get(), TColor::GetColor("#0072B2"), 21);
  setHistStyle(hPP.get(), kBlack, 24);

  TCanvas c("c_reco_xj_shape", "c_reco_xj_shape", 950, 760);
  gPad->SetTicks(1, 1);
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.04);
  gPad->SetTopMargin(0.12);
  gPad->SetBottomMargin(0.13);

  auto frame = std::make_unique<TH1F>("frame_xj", "", 100, 0.0, 2.0);
  frame->SetDirectory(nullptr);
  frame->SetStats(false);
  frame->SetMinimum(0.0);
  frame->SetMaximum(0.34);
  frame->GetXaxis()->SetTitle("x_{J#gamma}");
  frame->GetYaxis()->SetTitle("Normalized reco-level yield");
  frame->GetXaxis()->SetTitleSize(0.046);
  frame->GetYaxis()->SetTitleSize(0.046);
  frame->GetXaxis()->SetLabelSize(0.038);
  frame->GetYaxis()->SetLabelSize(0.038);
  frame->GetYaxis()->SetTitleOffset(1.25);
  frame->Draw();

  hPP->Draw("E SAME");
  hPeri->Draw("E SAME");
  hCent->Draw("E SAME");

  drawHeader(0.18, 0.87, "15 #leq E_{T}^{#gamma} < 35 GeV, R = 0.3, p_{T}^{jet} > 7 GeV");

  TLegend leg(0.52, 0.66, 0.89, 0.85);
  leg.SetBorderSize(0);
  leg.SetFillStyle(1001);
  leg.SetFillColorAlpha(kWhite, 0.88);
  leg.SetTextFont(42);
  leg.SetTextSize(0.035);
  leg.AddEntry(hPP.get(), "pp reference", "pe");
  leg.AddEntry(hPeri.get(), "Au+Au 50-80%", "pe");
  leg.AddEntry(hCent.get(), "Au+Au 0-20%", "pe");
  leg.Draw();

  c.SaveAs((kOutDir + "/reco_xj_shape_pp_peripheral_central_target80.png").c_str());
}

void drawFakeFraction()
{
  const ModelSpec model = kModels[1];
  std::unique_ptr<TFile> fAA(TFile::Open(signalPath(model).c_str(), "READ"));
  TDirectory* dAA = simDir(fAA.get());
  if (!dAA) return;

  TCanvas c("c_fake_fraction", "c_fake_fraction", 2250, 760);
  c.Divide(3, 1, 0.010, 0.0);
  std::vector<std::unique_ptr<TH1D>> keep;
  std::ofstream csv(kOutDir + "/reco_xj_fake_fraction.csv");
  csv << "centrality,xJ,fake,matched,fake_fraction\n";

  for (size_t ic = 0; ic < kCents.size(); ++ic)
  {
    c.cd(ic + 1);
    gPad->SetTicks(1, 1);
    gPad->SetLeftMargin(ic == 0 ? 0.13 : 0.065);
    gPad->SetRightMargin(ic == 2 ? 0.035 : 0.015);
    gPad->SetTopMargin(0.26);
    gPad->SetBottomMargin(0.13);

    const std::string suffix = centSuffix(kCents[ic]);
    auto hFake = projectXJ(getHist2(dAA, "h2_unfoldRecoFakes_pTgamma_xJ_incl_r03_jetPt7_dphiPi2_isoR30_isSliding" + suffix),
                           "h_fake_" + kCents[ic].suffix, 15.0, 35.0);
    auto hMatch = projectXJ(getHist2(dAA, "h2_unfoldRecoMatched_pTgamma_xJ_incl_r03_jetPt7_dphiPi2_isoR30_isSliding" + suffix),
                            "h_match_" + kCents[ic].suffix, 15.0, 35.0);
    auto hFrac = std::make_unique<TH1D>(("h_fakefrac_" + kCents[ic].suffix).c_str(), "", 20, 0.0, 2.0);
    hFrac->SetDirectory(nullptr);
    hFrac->SetStats(false);
    hFrac->SetMinimum(0.0);
    hFrac->SetMaximum(0.90);
    hFrac->GetXaxis()->SetTitle("x_{J#gamma}");
    hFrac->GetYaxis()->SetTitle(ic == 0 ? "Reco fake fraction" : "");
    hFrac->GetXaxis()->SetTitleSize(0.046);
    hFrac->GetYaxis()->SetTitleSize(0.046);
    hFrac->GetXaxis()->SetLabelSize(0.038);
    hFrac->GetYaxis()->SetLabelSize(0.038);
    hFrac->GetYaxis()->SetTitleOffset(ic == 0 ? 1.20 : 1.0);

    if (hFake && hMatch)
    {
      for (int ib = 1; ib <= hFrac->GetNbinsX(); ++ib)
      {
        const double x = hFrac->GetXaxis()->GetBinCenter(ib);
        const int jf = hFake->GetXaxis()->FindBin(x);
        const int jm = hMatch->GetXaxis()->FindBin(x);
        const double fake = hFake->GetBinContent(jf);
        const double match = hMatch->GetBinContent(jm);
        const double den = fake + match;
        const double frac = den > 0.0 ? fake / den : 0.0;
        hFrac->SetBinContent(ib, frac);
        hFrac->SetBinError(ib, den > 0.0 ? std::sqrt(frac * (1.0 - frac) / den) : 0.0);
        csv << kCents[ic].suffix << "," << x << "," << fake << "," << match << "," << frac << "\n";
      }
    }

    hFrac->SetMarkerStyle(20);
    hFrac->SetMarkerColor(TColor::GetColor("#CC79A7"));
    hFrac->SetLineColor(TColor::GetColor("#CC79A7"));
    hFrac->SetMarkerSize(1.1);
    hFrac->SetLineWidth(2);
    hFrac->Draw("E");

    TLatex title;
    title.SetNDC();
    title.SetTextFont(42);
    title.SetTextAlign(22);
    title.SetTextSize(0.052);
    title.DrawLatex(0.52, 0.955, kCents[ic].title.c_str());
    if (ic == 0) drawHeader(0.17, 0.86, "Centrality-input BDT target-80, 15 #leq E_{T}^{#gamma} < 35 GeV");
    keep.push_back(std::move(hFrac));
  }

  c.SaveAs((kOutDir + "/reco_xj_fake_fraction_1x3_target80.png").c_str());
}
}  // namespace

void PlotTarget80FoodChain()
{
  SetsPhenixStyle();
  gStyle->SetOptStat(0);
  gSystem->mkdir(kOutDir.c_str(), true);

  drawAbcdHeatmap(true);
  drawAbcdHeatmap(false);
  drawRecoXJOverlay();
  drawFakeFraction();

  std::cout << "[DONE] wrote target80 food-chain plots to " << kOutDir << "\n";
}
