#include "sPhenixStyle.C"

#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH1F.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TSystem.h>

#include <algorithm>
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
const std::string kFlatBase = "dataOutput/combinedSimOnlyEMBEDDED";
const std::string kOutDir = kTargetBase + "/first_look_plots";
const std::string kSignalRel =
    "photonJet12and20merged_SIM/RecoilJets_embeddedPhoton12plus20_MERGED.root";
const std::string kBkgRel =
    "embeddedJet12and20merged_SIM/RecoilJets_embeddedJet12plus20_MERGED.root";

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

struct Point
{
  double x = 0.0;
  double ex = 0.0;
  double y = 0.0;
  double ey = 0.0;
  double num = 0.0;
  double den = 0.0;
};

struct CurveSpec
{
  std::string key;
  std::string label;
  std::string cfg;
  int color = kBlack;
  int marker = 20;
  double shift = 0.0;
  bool target80 = true;
};

const std::vector<CentSpec> kCents = {
    {"0_20", "0-20% central"},
    {"20_50", "20-50% mid-central"},
    {"50_80", "50-80% peripheral"}};

const std::vector<PtBin> kPtBins = {
    {5, 8}, {8, 10}, {10, 12}, {12, 14}, {14, 16}, {16, 18},
    {18, 20}, {20, 22}, {22, 24}, {24, 26}, {26, 35}};

const std::vector<CurveSpec> kFamilyCurves = {
    {"box", "Box-cuts", "preselectionNewPPG12_tightReference_nonTightReference_baseVariant", kBlack, 20, -0.30, true},
    {"noCent", "BDT: no centrality input", "preselectionNewPPG12_tightAuAuNoCentBDT_nonTightAuAuBDTComplement_baseVariant", kGray + 2, 20, -0.20, true},
    {"centInput", "BDT: centrality as input", "preselectionNewPPG12_tightAuAuCentInputBDT_nonTightAuAuBDTComplement_baseVariant", kAzure + 2, 20, -0.10, true},
    {"minOpt", "BDT: minority-balanced", "preselectionNewPPG12_tightAuAuCentInputMinOptBDT_nonTightAuAuBDTComplement_baseVariant", kGreen + 2, 20, 0.00, true},
    {"cent3", "BDT: 3 centrality bins", "preselectionNewPPG12_tightAuAuCent3BDT_nonTightAuAuBDTComplement_baseVariant", kOrange + 7, 20, 0.10, true},
    {"cent7", "BDT: 7 centrality bins", "preselectionNewPPG12_tightAuAuCent7BDT_nonTightAuAuBDTComplement_baseVariant", kMagenta + 1, 20, 0.20, true},
    {"ptCent3", "BDT: E_{T} #times 3 cent.", "preselectionNewPPG12_tightAuAuPtCent3BDT_nonTightAuAuBDTComplement_baseVariant", kRed + 1, 20, 0.30, true}};

const std::vector<CurveSpec> kFlatVsTargetCurves = {
    {"box", "Box-cuts", "preselectionNewPPG12_tightReference_nonTightReference_baseVariant", kBlack, 20, -0.16, true},
    {"flat050", "Centrality-input BDT, score > 0.50", "preselectionNewPPG12_tightAuAuCentInputBDT_nonTightAuAuBDTComplement_baseVariant", kAzure + 2, 20, 0.00, false},
    {"target80", "Centrality-input BDT, target 80% signal efficiency", "preselectionNewPPG12_tightAuAuCentInputBDT_nonTightAuAuBDTComplement_baseVariant", kRed + 1, 20, 0.16, true}};

std::string signalPath(const CurveSpec& c)
{
  if (c.target80) return kTargetBase + "/simembedded/" + c.cfg + "/" + kSignalRel;
  return kFlatBase + "/" + c.cfg + "/" + kSignalRel;
}

std::string bkgPath(const CurveSpec& c)
{
  if (c.target80) return kTargetBase + "/simembeddedinclusive/" + c.cfg + "/" + kBkgRel;
  return kFlatBase + "/" + c.cfg + "/" + kBkgRel;
}

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

double integralAndError(TH1* h, double& err)
{
  err = 0.0;
  if (!h) return 0.0;
  h->SetStats(false);
  return h->IntegralAndError(0, h->GetNbinsX() + 1, err);
}

Point ratioPoint(double num, double denPart, double eNum, double eDenPart, double x)
{
  Point p;
  p.x = x;
  p.num = num;
  p.den = num + denPart;
  if (p.den <= 0.0) return p;
  p.y = num / p.den;
  p.ey = std::sqrt(denPart * denPart * eNum * eNum + num * num * eDenPart * eDenPart) / (p.den * p.den);
  if (!std::isfinite(p.ey)) p.ey = 0.0;
  return p;
}

std::vector<Point> readRecovery(TDirectory* d, const CentSpec& cent, double ptMin, double ptMax)
{
  std::vector<Point> pts;
  TH1* hTruth = getHist(d, "h_unfoldTruthPho_pTgamma_cent_" + cent.suffix);
  TH1* hMiss = getHist(d, "h_unfoldTruthPhoMisses_pTgamma_isoR30_isSliding_cent_" + cent.suffix);
  if (!hTruth || !hMiss) return pts;

  for (int ib = 1; ib <= hTruth->GetNbinsX(); ++ib)
  {
    const double lo = hTruth->GetXaxis()->GetBinLowEdge(ib);
    const double hi = hTruth->GetXaxis()->GetBinUpEdge(ib);
    const double mid = 0.5 * (lo + hi);
    if (mid < ptMin || mid > ptMax) continue;
    const int im = hMiss->GetXaxis()->FindBin(mid);
    if (im < 1 || im > hMiss->GetNbinsX()) continue;
    const double truth = hTruth->GetBinContent(ib);
    const double miss = hMiss->GetBinContent(im);
    if (truth <= 0.0) continue;

    Point p;
    p.x = mid;
    p.y = 1.0 - miss / truth;
    p.num = truth - miss;
    p.den = truth;
    const double eTruth = hTruth->GetBinError(ib);
    const double eMiss = hMiss->GetBinError(im);
    if (miss > 0.0)
    {
      const double missFrac = miss / truth;
      p.ey = missFrac * std::hypot(eMiss / miss, eTruth / truth);
    }
    if (!std::isfinite(p.ey)) p.ey = 0.0;
    pts.push_back(p);
  }
  return pts;
}

Point readIntegratedRecovery(TDirectory* d, const CentSpec& cent, double ptMin, double ptMax)
{
  Point out;
  TH1* hTruth = getHist(d, "h_unfoldTruthPho_pTgamma_cent_" + cent.suffix);
  TH1* hMiss = getHist(d, "h_unfoldTruthPhoMisses_pTgamma_isoR30_isSliding_cent_" + cent.suffix);
  if (!hTruth || !hMiss) return out;

  double truth = 0.0;
  double miss = 0.0;
  double eTruth2 = 0.0;
  double eMiss2 = 0.0;
  for (int ib = 1; ib <= hTruth->GetNbinsX(); ++ib)
  {
    const double lo = hTruth->GetXaxis()->GetBinLowEdge(ib);
    const double hi = hTruth->GetXaxis()->GetBinUpEdge(ib);
    const double mid = 0.5 * (lo + hi);
    if (mid < ptMin || mid > ptMax) continue;
    const int im = hMiss->GetXaxis()->FindBin(mid);
    if (im < 1 || im > hMiss->GetNbinsX()) continue;
    truth += hTruth->GetBinContent(ib);
    miss += hMiss->GetBinContent(im);
    eTruth2 += hTruth->GetBinError(ib) * hTruth->GetBinError(ib);
    eMiss2 += hMiss->GetBinError(im) * hMiss->GetBinError(im);
  }
  out.num = truth - miss;
  out.den = truth;
  if (truth <= 0.0) return out;
  out.y = 1.0 - miss / truth;
  const double eTruth = std::sqrt(eTruth2);
  const double eMiss = std::sqrt(eMiss2);
  if (miss > 0.0)
  {
    const double missFrac = miss / truth;
    out.ey = missFrac * std::hypot(eMiss / miss, eTruth / truth);
  }
  if (!std::isfinite(out.ey)) out.ey = 0.0;
  return out;
}

Point readCandidateTightFraction(TDirectory* d, const PtBin& pt, const CentSpec& cent)
{
  const std::string tag = "isoR30_pT_" + std::to_string(pt.lo) + "_" +
                          std::to_string(pt.hi) + "_cent_" + cent.suffix;
  double eTight = 0.0;
  double eNonTight = 0.0;
  const double tight = integralAndError(getHist(d, "h_Eiso_tight_" + tag), eTight);
  const double nonTight = integralAndError(getHist(d, "h_Eiso_nonTight_" + tag), eNonTight);
  return ratioPoint(tight, nonTight, eTight, eNonTight, 0.5 * (pt.lo + pt.hi));
}

std::vector<Point> readCandidateTightFractions(TDirectory* d, const CentSpec& cent, double ptMin, double ptMax)
{
  std::vector<Point> pts;
  for (const auto& pt : kPtBins)
  {
    const double mid = 0.5 * (pt.lo + pt.hi);
    if (mid < ptMin || mid > ptMax) continue;
    Point p = readCandidateTightFraction(d, pt, cent);
    if (p.den > 0.0) pts.push_back(p);
  }
  return pts;
}

Point readIntegratedCandidateTightFraction(TDirectory* d, const CentSpec& cent, double ptMin, double ptMax)
{
  double tight = 0.0;
  double nonTight = 0.0;
  double eTight2 = 0.0;
  double eNonTight2 = 0.0;
  for (const auto& pt : kPtBins)
  {
    const double mid = 0.5 * (pt.lo + pt.hi);
    if (mid < ptMin || mid > ptMax) continue;
    const std::string tag = "isoR30_pT_" + std::to_string(pt.lo) + "_" +
                            std::to_string(pt.hi) + "_cent_" + cent.suffix;
    double eT = 0.0;
    double eN = 0.0;
    tight += integralAndError(getHist(d, "h_Eiso_tight_" + tag), eT);
    nonTight += integralAndError(getHist(d, "h_Eiso_nonTight_" + tag), eN);
    eTight2 += eT * eT;
    eNonTight2 += eN * eN;
  }
  return ratioPoint(tight, nonTight, std::sqrt(eTight2), std::sqrt(eNonTight2), 0.0);
}

std::unique_ptr<TGraphErrors> graphFromPoints(const CurveSpec& curve,
                                              const std::vector<Point>& pts,
                                              std::ofstream& csv,
                                              const std::string& quantity,
                                              const std::string& cent)
{
  auto g = std::make_unique<TGraphErrors>();
  g->SetMarkerStyle(curve.marker);
  g->SetMarkerColor(curve.color);
  g->SetLineColor(curve.color);
  g->SetMarkerSize(1.14);
  g->SetLineWidth(2);
  for (const auto& p : pts)
  {
    const int i = g->GetN();
    g->SetPoint(i, p.x + curve.shift, p.y);
    g->SetPointError(i, p.ex, p.ey);
    csv << quantity << "," << curve.key << "," << curve.label << "," << cent << ","
        << std::setprecision(10) << p.x << "," << p.num << "," << p.den << ","
        << p.y << "," << p.ey << "\n";
  }
  return g;
}

void drawHeader(double x = 0.125, double y = 0.86)
{
  TLatex tx;
  tx.SetNDC();
  tx.SetTextFont(42);
  tx.SetTextAlign(13);
  tx.SetTextSize(0.036);
  tx.DrawLatex(x, y, "#it{#bf{sPHENIX}} Internal");
  tx.SetTextSize(0.028);
  tx.DrawLatex(x, y - 0.055, "Pythia overlay, #sqrt{s_{NN}} = 200 GeV");
  tx.DrawLatex(x, y - 0.100, "R = 0.3, sliding isolation");
}

void drawPanelTitle(const std::string& title)
{
  TLatex tx;
  tx.SetNDC();
  tx.SetTextFont(42);
  tx.SetTextAlign(22);
  tx.SetTextSize(0.052);
  tx.DrawLatex(0.50, 0.955, title.c_str());
}

std::unique_ptr<TFile> openFile(const std::string& path)
{
  std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
  if (!f || f->IsZombie())
  {
    std::cerr << "[WARN] could not open " << path << "\n";
    return nullptr;
  }
  return f;
}

void drawLegend(const std::vector<CurveSpec>& curves, double x1, double y1, double x2, double y2, int nColumns,
                std::vector<std::unique_ptr<TGraphErrors>>& keep)
{
  TLegend* leg = new TLegend(x1, y1, x2, y2);
  leg->SetBorderSize(0);
  leg->SetFillStyle(1001);
  leg->SetFillColorAlpha(kWhite, 0.88);
  leg->SetTextFont(42);
  leg->SetTextSize(0.031);
  leg->SetNColumns(nColumns);
  for (const auto& c : curves)
  {
    auto g = std::make_unique<TGraphErrors>(1);
    g->SetMarkerStyle(c.marker);
    g->SetMarkerSize(1.14);
    g->SetMarkerColor(c.color);
    g->SetLineColor(c.color);
    leg->AddEntry(g.get(), c.label.c_str(), "pe");
    keep.push_back(std::move(g));
  }
  leg->Draw();
}

void drawOneByThree(const std::vector<CurveSpec>& curves,
                    const std::string& outName,
                    const std::string& yTitle,
                    double yMin,
                    double yMax,
                    bool useRecovery,
                    double ptMin,
                    double ptMax)
{
  TCanvas c(("c_" + outName).c_str(), ("c_" + outName).c_str(), 2300, 780);
  c.Divide(3, 1, 0.010, 0.0);
  std::vector<std::unique_ptr<TFile>> files;
  std::vector<std::unique_ptr<TGraphErrors>> keep;
  std::vector<std::unique_ptr<TH1F>> frames;
  std::vector<std::vector<std::unique_ptr<TGraphErrors>>> graphs(kCents.size());
  std::ofstream csv(kOutDir + "/" + outName + ".csv");
  csv << "quantity,model_key,model_label,centrality,pt_mid,numerator,denominator,value,error\n";

  for (size_t ic = 0; ic < kCents.size(); ++ic)
  {
    c.cd(ic + 1);
    gPad->SetTicks(1, 1);
    gPad->SetLeftMargin(ic == 0 ? 0.115 : 0.055);
    gPad->SetRightMargin(ic == kCents.size() - 1 ? 0.035 : 0.015);
    gPad->SetTopMargin(0.30);
    gPad->SetBottomMargin(0.14);

    auto frame = std::make_unique<TH1F>(("frame_" + outName + "_" + kCents[ic].suffix).c_str(), "", 100, 4.5, 35.5);
    frame->SetDirectory(nullptr);
    frame->SetStats(false);
    frame->SetMinimum(yMin);
    frame->SetMaximum(yMax);
    frame->GetXaxis()->SetTitle(useRecovery ? "Truth photon E_{T} [GeV]" : "Photon candidate E_{T} [GeV]");
    frame->GetYaxis()->SetTitle(ic == 0 ? yTitle.c_str() : "");
    frame->GetXaxis()->SetTitleSize(0.044);
    frame->GetYaxis()->SetTitleSize(0.044);
    frame->GetXaxis()->SetLabelSize(0.037);
    frame->GetYaxis()->SetLabelSize(0.037);
    frame->GetYaxis()->SetTitleOffset(ic == 0 ? 1.15 : 1.0);
    frame->Draw();
    frames.push_back(std::move(frame));
    drawPanelTitle(kCents[ic].title);
    if (ic == 0) drawHeader();

    for (const auto& curve : curves)
    {
      auto f = openFile(useRecovery ? signalPath(curve) : bkgPath(curve));
      if (!f) continue;
      TDirectory* d = simDir(f.get());
      if (!d)
      {
        std::cerr << "[WARN] missing SIM directory for " << curve.label << "\n";
        continue;
      }
      const auto pts = useRecovery ? readRecovery(d, kCents[ic], ptMin, ptMax)
                                   : readCandidateTightFractions(d, kCents[ic], ptMin, ptMax);
      auto g = graphFromPoints(curve, pts, csv, useRecovery ? "truth_photon_recovery" : "inclusive_jet_tight_fraction", kCents[ic].suffix);
      if (g->GetN() > 0)
      {
        g->Draw("P SAME");
        graphs[ic].push_back(std::move(g));
      }
      files.push_back(std::move(f));
    }

    if (ic == 2)
    {
      drawLegend(curves, 0.05, 0.735, 0.96, 0.94, curves.size() > 4 ? 2 : 1, keep);
    }
  }

  c.SaveAs((kOutDir + "/" + outName + ".png").c_str());
}

void drawTradeoff(const std::vector<CurveSpec>& curves, double ptMin, double ptMax)
{
  TCanvas c("c_target80_tradeoff", "c_target80_tradeoff", 2250, 760);
  c.Divide(3, 1, 0.010, 0.0);
  std::ofstream csv(kOutDir + "/target80_signal_recovery_vs_background_tight_fraction.csv");
  csv << "model_key,model_label,centrality,signal_recovery,signal_recovery_err,background_tight_fraction,background_tight_fraction_err\n";

  std::vector<std::unique_ptr<TFile>> files;
  std::vector<std::unique_ptr<TGraphErrors>> legendKeep;
  std::vector<std::unique_ptr<TGraphErrors>> graphs;
  std::vector<std::unique_ptr<TH1F>> frames;

  for (size_t ic = 0; ic < kCents.size(); ++ic)
  {
    c.cd(ic + 1);
    gPad->SetTicks(1, 1);
    gPad->SetLeftMargin(ic == 0 ? 0.13 : 0.070);
    gPad->SetRightMargin(ic == kCents.size() - 1 ? 0.035 : 0.020);
    gPad->SetTopMargin(0.28);
    gPad->SetBottomMargin(0.14);

    auto frame = std::make_unique<TH1F>(("frame_tradeoff_" + kCents[ic].suffix).c_str(), "", 100, 0.0, 0.46);
    frame->SetDirectory(nullptr);
    frame->SetStats(false);
    frame->SetMinimum(0.00);
    frame->SetMaximum(0.64);
    frame->GetXaxis()->SetTitle("Inclusive-jet tight fraction");
    frame->GetYaxis()->SetTitle(ic == 0 ? "Truth photon recovery" : "");
    frame->GetXaxis()->SetTitleSize(0.043);
    frame->GetYaxis()->SetTitleSize(0.043);
    frame->GetXaxis()->SetLabelSize(0.036);
    frame->GetYaxis()->SetLabelSize(0.036);
    frame->GetYaxis()->SetTitleOffset(ic == 0 ? 1.35 : 1.00);
    frame->Draw();
    frames.push_back(std::move(frame));
    drawPanelTitle(kCents[ic].title);
    if (ic == 0) drawHeader(0.15, 0.88);

    TLatex arrow;
    arrow.SetNDC();
    arrow.SetTextFont(42);
    arrow.SetTextSize(0.030);
    arrow.SetTextColor(kGray + 2);
    arrow.DrawLatex(0.53, 0.19, "better: higher recovery, lower jet leakage");

    for (const auto& curve : curves)
    {
      auto fs = openFile(signalPath(curve));
      auto fb = openFile(bkgPath(curve));
      if (!fs || !fb) continue;
      TDirectory* ds = simDir(fs.get());
      TDirectory* db = simDir(fb.get());
      if (!ds || !db) continue;

      const Point rec = readIntegratedRecovery(ds, kCents[ic], ptMin, ptMax);
      const Point leak = readIntegratedCandidateTightFraction(db, kCents[ic], ptMin, ptMax);
      csv << curve.key << "," << curve.label << "," << kCents[ic].suffix << ","
          << std::setprecision(10) << rec.y << "," << rec.ey << ","
          << leak.y << "," << leak.ey << "\n";

      auto g = std::make_unique<TGraphErrors>(1);
      g->SetPoint(0, leak.y, rec.y);
      g->SetPointError(0, leak.ey, rec.ey);
      g->SetMarkerStyle(curve.marker);
      g->SetMarkerSize(1.45);
      g->SetMarkerColor(curve.color);
      g->SetLineColor(curve.color);
      g->SetLineWidth(2);
      g->Draw("P SAME");
      graphs.push_back(std::move(g));
      files.push_back(std::move(fs));
      files.push_back(std::move(fb));
    }

    if (ic == 2)
    {
      drawLegend(curves, 0.05, 0.735, 0.96, 0.94, 2, legendKeep);
    }
  }

  c.SaveAs((kOutDir + "/target80_signal_recovery_vs_background_tight_fraction.png").c_str());
}

void writeMetricSummary(const std::vector<CurveSpec>& curves, double ptMin, double ptMax)
{
  std::ofstream csv(kOutDir + "/target80_integrated_metric_summary_15to35.csv");
  csv << "model_key,model_label,centrality,signal_recovery,signal_recovery_err,signal_candidate_tight_fraction,signal_candidate_tight_fraction_err,background_tight_fraction,background_tight_fraction_err\n";

  for (const auto& curve : curves)
  {
    auto fs = openFile(signalPath(curve));
    auto fb = openFile(bkgPath(curve));
    if (!fs || !fb) continue;
    TDirectory* ds = simDir(fs.get());
    TDirectory* db = simDir(fb.get());
    if (!ds || !db) continue;
    for (const auto& cent : kCents)
    {
      const Point rec = readIntegratedRecovery(ds, cent, ptMin, ptMax);
      const Point sigTight = readIntegratedCandidateTightFraction(ds, cent, ptMin, ptMax);
      const Point bkgTight = readIntegratedCandidateTightFraction(db, cent, ptMin, ptMax);
      csv << curve.key << "," << curve.label << "," << cent.suffix << ","
          << std::setprecision(10) << rec.y << "," << rec.ey << ","
          << sigTight.y << "," << sigTight.ey << ","
          << bkgTight.y << "," << bkgTight.ey << "\n";
    }
  }
}
}  // namespace

void PlotTarget80Expanded5to40Campaign()
{
  SetsPhenixStyle();
  gStyle->SetOptStat(0);
  gSystem->mkdir(kOutDir.c_str(), true);

  const double ptMin = 5.0;
  const double ptMax = 35.0;
  const double tradePtMin = 15.0;
  const double tradePtMax = 35.0;

  drawOneByThree(kFamilyCurves,
                 "target80_model_family_truth_photon_recovery_1x3",
                 "Truth photon recovery",
                 0.0,
                 0.62,
                 true,
                 ptMin,
                 ptMax);

  drawOneByThree(kFamilyCurves,
                 "target80_model_family_inclusive_jet_tight_fraction_1x3",
                 "Inclusive-jet tight fraction",
                 0.0,
                 0.52,
                 false,
                 ptMin,
                 ptMax);

  drawOneByThree(kFlatVsTargetCurves,
                 "centinput_reference_flat050_target80_truth_photon_recovery_1x3",
                 "Truth photon recovery",
                 0.0,
                 0.62,
                 true,
                 ptMin,
                 ptMax);

  drawOneByThree(kFlatVsTargetCurves,
                 "centinput_reference_flat050_target80_inclusive_jet_tight_fraction_1x3",
                 "Inclusive-jet tight fraction",
                 0.0,
                 0.52,
                 false,
                 ptMin,
                 ptMax);

  drawTradeoff(kFamilyCurves, tradePtMin, tradePtMax);
  writeMetricSummary(kFamilyCurves, tradePtMin, tradePtMax);

  std::cout << "[DONE] wrote target80 expanded 5-40 first-look plots to " << kOutDir << "\n";
}
