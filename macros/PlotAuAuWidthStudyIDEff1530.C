#include "sPhenixStyle.C"

#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH1F.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TPad.h>
#include <TStyle.h>
#include <TSystem.h>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace
{
const std::string kBaseDir =
    "dataOutput/auau_widthstudy_pt1530_wp050/combinedSimOnlyEMBEDDED";
const std::string kSignalRel =
    "photonJet12and20merged_SIM/RecoilJets_embeddedPhoton12plus20_MERGED.root";
const std::string kOutDir =
    "dataOutput/jstg_slide_candidates/widthstudy_id_efficiency_15to30";

struct CurveSpec
{
  std::string baseDir;
  std::string cfg;
  std::string label;
  int color = kBlack;
  int marker = 20;
  double xShift = 0.0;
  double ptMin = 15.0;
  double ptMax = 30.0;
};

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

struct EffPoint
{
  double x = 0.0;
  double y = 0.0;
  double ey = 0.0;
  double truth = 0.0;
  double miss = 0.0;
  double lo = 0.0;
  double hi = 0.0;
};

enum class EfficiencyMode
{
  TruthPhotonEt,
  CandidateEt
};

const std::vector<CentSpec> kCents = {
    {"0_20", "0-20% central"},
    {"20_50", "20-50% mid-central"},
    {"50_80", "50-80% peripheral"}};

const std::vector<PtBin> kCandidatePtBins = {
    {5, 8}, {8, 10}, {10, 12}, {12, 14}, {14, 16}, {16, 18},
    {18, 20}, {20, 22}, {22, 24}, {24, 26}, {26, 35}};

const std::string kOlderBaseDir = "dataOutput/combinedSimOnlyEMBEDDED";

const CurveSpec kBoxCuts = {
    kBaseDir,
    "preselectionNewPPG12_tightReference_nonTightReference_baseVariant",
    "box-cuts", kBlack, 20, -0.12, 5.0, 35.0};

const CurveSpec kBoxCutsWide = {
    kBaseDir,
    "preselectionNewPPG12_tightReference_nonTightReference_baseVariant",
    "box-cuts", kBlack, 20, -0.12, 5.0, 35.0};

const CurveSpec kCurrentBase = {
    kBaseDir,
    "preselectionNewPPG12_tightAuAuCentInputBDT_nonTightAuAuBDTComplement_baseVariant",
    "base widths only, 15-30 training", kAzure + 2, 20, 0.00, 15.0, 30.0};

const CurveSpec kCurrentBase3x3 = {
    kBaseDir,
    "preselectionNewPPG12_tightAuauCentInputBase3x3BDT_nonTightAuAuBDTComplement_baseVariant",
    "base + 3#times3 widths, 15-30 training", kRed + 1, 20, 0.12, 15.0, 30.0};

const CurveSpec kOlderBase = {
    kOlderBaseDir,
    "preselectionNewPPG12_tightAuAuCentInputBDT_nonTightAuAuBDTComplement_baseVariant",
    "base widths only, 5-40 training", kAzure + 2, 20, 0.08, 5.0, 35.0};

const CurveSpec kCurrentBaseSlide33 = {
    kBaseDir,
    "preselectionNewPPG12_tightAuAuCentInputBDT_nonTightAuAuBDTComplement_baseVariant",
    "base widths only, 15-30 training", kRed + 1, 20, 0.16, 15.0, 30.0};

const CurveSpec kSlide33BoxCuts = {
    kOlderBaseDir,
    "preselectionNewPPG12_tightReference_nonTightReference_baseVariant",
    "box-cuts", kBlack, 20, -0.16, 5.0, 35.0};

const std::vector<CurveSpec> kCurves1530 = {kBoxCuts, kCurrentBase, kCurrentBase3x3};
const std::vector<CurveSpec> kCurvesWide = {kBoxCutsWide, kCurrentBase, kCurrentBase3x3, kOlderBase};
const std::vector<CurveSpec> kSlide33Requested = {kSlide33BoxCuts, kCurrentBaseSlide33, kOlderBase};

std::string signalPath(const CurveSpec& c)
{
  return c.baseDir + "/" + c.cfg + "/" + kSignalRel;
}

TH1* getHist(TDirectory* d, const std::string& name)
{
  if (!d) return nullptr;
  TH1* h = dynamic_cast<TH1*>(d->Get(name.c_str()));
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

std::vector<EffPoint> readPoints(TDirectory* d, const CentSpec& cent, const CurveSpec& curve)
{
  std::vector<EffPoint> pts;
  TH1* hTruth = getHist(d, "h_unfoldTruthPho_pTgamma_cent_" + cent.suffix);
  TH1* hMiss = getHist(d, "h_unfoldTruthPhoMisses_pTgamma_isoR30_isSliding_cent_" + cent.suffix);
  if (!hTruth || !hMiss)
  {
    std::cerr << "[WARN] missing efficiency histograms for " << cent.suffix << "\n";
    return pts;
  }

  for (int ib = 1; ib <= hTruth->GetNbinsX(); ++ib)
  {
    const double lo = hTruth->GetXaxis()->GetBinLowEdge(ib);
    const double hi = hTruth->GetXaxis()->GetBinUpEdge(ib);
    if (lo < curve.ptMin || hi > curve.ptMax) continue;

    const double mid = 0.5 * (lo + hi);
    const int im = hMiss->GetXaxis()->FindBin(mid);
    if (im < 1 || im > hMiss->GetNbinsX()) continue;

    const double truth = hTruth->GetBinContent(ib);
    const double miss = hMiss->GetBinContent(im);
    if (truth <= 0.0) continue;

    const double eTruth = hTruth->GetBinError(ib);
    const double eMiss = hMiss->GetBinError(im);
    const double missFrac = miss / truth;

    EffPoint p;
    p.x = mid;
    p.y = 1.0 - missFrac;
    p.ey = miss > 0.0 ? missFrac * std::hypot(eMiss / miss, eTruth / truth) : 0.0;
    if (!std::isfinite(p.ey)) p.ey = 0.0;
    p.truth = truth;
    p.miss = miss;
    p.lo = lo;
    p.hi = hi;
    pts.push_back(p);
  }
  return pts;
}

std::vector<EffPoint> readCandidateEtPoints(TDirectory* d,
                                            const CentSpec& cent,
                                            const CurveSpec& curve)
{
  std::vector<EffPoint> pts;
  for (const auto& bin : kCandidatePtBins)
  {
    if (bin.lo < curve.ptMin || bin.hi > curve.ptMax) continue;

    const std::string tag =
        "isoR30_pT_" + std::to_string(bin.lo) + "_" + std::to_string(bin.hi) +
        "_cent_" + cent.suffix;
    double eTight = 0.0;
    double eNonTight = 0.0;
    const double tight = integralAndError(getHist(d, "h_Eiso_tight_" + tag), eTight);
    const double nonTight = integralAndError(getHist(d, "h_Eiso_nonTight_" + tag), eNonTight);
    const double den = tight + nonTight;
    if (den <= 0.0) continue;

    const double eDen = std::hypot(eTight, eNonTight);
    EffPoint p;
    p.x = 0.5 * (bin.lo + bin.hi);
    p.y = tight / den;
    p.ey = tight > 0.0 ? p.y * std::hypot(eTight / tight, eDen / den) : 0.0;
    if (!std::isfinite(p.ey)) p.ey = 0.0;
    p.truth = tight;
    p.miss = nonTight;
    p.lo = bin.lo;
    p.hi = bin.hi;
    pts.push_back(p);
  }
  return pts;
}

std::unique_ptr<TGraphErrors> makeGraph(const CurveSpec& curve,
                                        TDirectory* d,
                                        const CentSpec& cent,
                                        std::ofstream& csv,
                                        EfficiencyMode mode,
                                        double markerScale = 1.0)
{
  auto pts = (mode == EfficiencyMode::TruthPhotonEt)
                 ? readPoints(d, cent, curve)
                 : readCandidateEtPoints(d, cent, curve);
  if (pts.empty()) return nullptr;

  auto g = std::make_unique<TGraphErrors>();
  g->SetName(("g_id_eff_" + cent.suffix + "_" + curve.label).c_str());
  g->SetMarkerStyle(curve.marker);
  g->SetMarkerSize(1.48 * markerScale);
  g->SetMarkerColor(curve.color);
  g->SetLineColor(curve.color);
  g->SetLineWidth(2);

  for (const auto& p : pts)
  {
    const int ip = g->GetN();
    g->SetPoint(ip, p.x + curve.xShift, p.y);
    g->SetPointError(ip, 0.0, p.ey);
    csv << curve.cfg << "," << curve.label << "," << cent.suffix << ","
        << std::setprecision(8) << p.lo << "," << p.hi << "," << p.x << ","
        << p.truth << "," << p.miss << "," << p.y << "," << p.ey << "\n";
  }
  return g;
}

void drawPanelLabel(const std::string& text, double y = 0.955)
{
  TLatex tx;
  tx.SetNDC();
  tx.SetTextFont(42);
  tx.SetTextAlign(22);
  tx.SetTextSize(0.050);
  tx.DrawLatex(0.50, y, text.c_str());
}

void drawHeader(const std::string& scope, double x = 0.135)
{
  TLatex tx;
  tx.SetNDC();
  tx.SetTextFont(42);
  tx.SetTextAlign(13);
  tx.SetTextSize(0.040);
  tx.DrawLatex(x, 0.835, "#it{#bf{sPHENIX}} Internal");
  tx.SetTextSize(0.031);
  tx.DrawLatex(x, 0.765, "Pythia overlay, #sqrt{s_{NN}} = 200 GeV");
  tx.DrawLatex(x, 0.700, scope.c_str());
}

void drawLegend(const std::vector<CurveSpec>& curves, bool compact)
{
  auto leg = compact ? new TLegend(0.39, 0.700, 0.94, 0.890)
                     : new TLegend(0.38, 0.680, 0.94, 0.895);
  leg->SetBorderSize(0);
  leg->SetFillStyle(1001);
  leg->SetFillColorAlpha(kWhite, 0.88);
  leg->SetTextFont(42);
  leg->SetTextSize(compact ? 0.030 : 0.027);
  leg->SetNColumns(1);

  static std::vector<std::unique_ptr<TGraphErrors>> keepAlive;
  keepAlive.clear();
  for (const auto& curve : curves)
  {
    auto g = std::make_unique<TGraphErrors>(1);
    g->SetMarkerStyle(curve.marker);
    g->SetMarkerSize(1.48);
    g->SetMarkerColor(curve.color);
    g->SetLineColor(curve.color);
    leg->AddEntry(g.get(), curve.label.c_str(), "pe");
    keepAlive.push_back(std::move(g));
  }
  leg->Draw();
}

void drawGlobalLegend(const std::vector<CurveSpec>& curves, double markerScale = 1.0)
{
  auto leg = new TLegend(0.300, 0.905, 0.970, 0.992);
  leg->SetBorderSize(0);
  leg->SetFillStyle(1001);
  leg->SetFillColorAlpha(kWhite, 0.94);
  leg->SetTextFont(42);
  leg->SetTextSize(0.034);
  leg->SetNColumns(3);

  static std::vector<std::unique_ptr<TGraphErrors>> keepAlive;
  keepAlive.clear();
  for (const auto& curve : curves)
  {
    auto g = std::make_unique<TGraphErrors>(1);
    g->SetMarkerStyle(curve.marker);
    g->SetMarkerSize(1.95 * markerScale);
    g->SetMarkerColor(curve.color);
    g->SetLineColor(curve.color);
    leg->AddEntry(g.get(), curve.label.c_str(), "pe");
    keepAlive.push_back(std::move(g));
  }
  leg->Draw();
}

void drawGlobalHeader()
{
  TLatex tx;
  tx.SetNDC();
  tx.SetTextFont(42);
  tx.SetTextAlign(13);
  tx.SetTextSize(0.037);
  tx.DrawLatex(0.040, 0.980, "#it{#bf{sPHENIX}} Internal");
  tx.SetTextSize(0.029);
  tx.DrawLatex(0.040, 0.935, "Pythia overlay, #sqrt{s_{NN}} = 200 GeV");
}

void drawOverlay(const std::string& outStem,
                 const std::vector<CurveSpec>& curves,
                 double xMin,
                 double xMax,
                 double yMin,
                 double yMax,
                 const std::string& headerScope,
                 const std::string& csvName,
                 double headerX = 0.135,
                 bool globalLegend = false,
                 double markerScale = 1.0,
                 EfficiencyMode mode = EfficiencyMode::TruthPhotonEt,
                 const std::string& xAxisTitle = "Truth photon E_{T} [GeV]")
{
  std::ofstream csv(kOutDir + "/" + csvName);
  if (mode == EfficiencyMode::TruthPhotonEt)
    csv << "config_key,label,cent,pt_low,pt_high,pt_mid,truth_sum,miss_sum,id_eff,id_eff_err\n";
  else
    csv << "config_key,label,cent,pt_low,pt_high,pt_mid,tight_sum,nontight_sum,id_eff,id_eff_err\n";

  TCanvas c(("c_" + outStem).c_str(), ("c_" + outStem).c_str(), 3590, 1159);
  if (!globalLegend) c.Divide(3, 1, 0.010, 0.0);

  std::vector<std::unique_ptr<TFile>> openFiles;
  std::vector<std::unique_ptr<TGraphErrors>> graphs;
  std::vector<std::unique_ptr<TH1F>> frames;
  std::vector<std::unique_ptr<TPad>> pads;

  for (size_t ic = 0; ic < kCents.size(); ++ic)
  {
    if (globalLegend)
    {
      const double x0 = static_cast<double>(ic) / 3.0 + (ic == 0 ? 0.030 : 0.012);
      const double x1 = static_cast<double>(ic + 1) / 3.0 - (ic == 2 ? 0.018 : 0.010);
      auto pad = std::make_unique<TPad>(("pad_" + outStem + "_" + kCents[ic].suffix).c_str(), "", x0, 0.030, x1, 0.890);
      pad->SetFillStyle(4000);
      pad->SetFrameFillStyle(4000);
      c.cd();
      pad->Draw();
      pad->cd();
      pads.push_back(std::move(pad));
    }
    else
    {
      c.cd(ic + 1);
    }
    gPad->SetTicks(1, 1);
    gPad->SetGrid(0, 0);
    gPad->SetLeftMargin(ic == 0 ? 0.150 : 0.060);
    gPad->SetRightMargin(ic == kCents.size() - 1 ? 0.035 : 0.020);
    gPad->SetTopMargin(globalLegend ? 0.170 : 0.145);
    gPad->SetBottomMargin(0.145);

    auto frame = std::make_unique<TH1F>(("hframe_" + outStem + "_" + kCents[ic].suffix).c_str(), "", 100, xMin, xMax);
    frame->SetDirectory(nullptr);
    frame->SetStats(false);
    frame->SetMinimum(yMin);
    frame->SetMaximum(yMax);
    frame->GetXaxis()->SetTitle(xAxisTitle.c_str());
    frame->GetYaxis()->SetTitle(ic == 0 ? "Tight-ID efficiency" : "");
    frame->GetXaxis()->SetTitleSize(0.046);
    frame->GetYaxis()->SetTitleSize(0.046);
    frame->GetXaxis()->SetLabelSize(0.038);
    frame->GetYaxis()->SetLabelSize(ic == 0 ? 0.038 : 0.0);
    frame->GetYaxis()->SetTitleOffset(ic == 0 ? 1.18 : 0.95);
    frame->GetXaxis()->SetNdivisions(506);
    frame->GetYaxis()->SetNdivisions(506);
    frame->Draw("AXIS");
    frames.push_back(std::move(frame));

    drawPanelLabel(kCents[ic].title, globalLegend ? 0.925 : 0.955);
    if (ic == 0 && !globalLegend) drawHeader(headerScope, headerX);
    if (ic == 2 && !globalLegend) drawLegend(curves, curves.size() <= 3);

    for (const auto& curve : curves)
    {
      const std::string path = signalPath(curve);
      auto f = std::unique_ptr<TFile>(TFile::Open(path.c_str(), "READ"));
      if (!f || f->IsZombie())
      {
        std::cerr << "[WARN] missing file: " << path << "\n";
        continue;
      }
      TDirectory* d = dynamic_cast<TDirectory*>(f->Get("SIM"));
      if (!d)
      {
        std::cerr << "[WARN] missing SIM directory: " << path << "\n";
        continue;
      }
      auto g = makeGraph(curve, d, kCents[ic], csv, mode, markerScale);
      if (g && g->GetN() > 0)
      {
        g->Draw("P SAME");
        graphs.push_back(std::move(g));
      }
      openFiles.push_back(std::move(f));
    }

    gPad->RedrawAxis();
  }

  if (globalLegend)
  {
    c.cd();
    drawGlobalHeader();
    drawGlobalLegend(curves, markerScale);
  }

  const std::string out = kOutDir + "/" + outStem + ".png";
  c.SaveAs(out.c_str());
  std::cout << "[DONE] " << out << "\n";
}
}  // namespace

void PlotAuAuWidthStudyIDEff1530()
{
  SetsPhenixStyle();
  gStyle->SetOptStat(0);
  gStyle->SetPadGridX(false);
  gStyle->SetPadGridY(false);
  gSystem->mkdir(kOutDir.c_str(), true);

  drawOverlay("id_efficiency_boxcuts_base_vs_base3x3_wp050_15to30",
              kCurves1530,
              4.7,
              35.3,
              0.00,
              0.62,
              "box-cuts shown over full range; BDT points use 15-30 training domain",
              "id_efficiency_boxcuts_base_vs_base3x3_wp050_15to30_points.csv",
              0.180);

  drawOverlay("id_efficiency_context_older5to35_vs_widthstudy1530_wp050",
              kCurvesWide,
              4.7,
              35.3,
              0.00,
              0.62,
              "current 15-30 BDT points with earlier 5-35 output, BDT > 0.50",
              "id_efficiency_context_older5to35_vs_widthstudy1530_wp050_points.csv",
              0.180);

  drawOverlay("id_efficiency_slide33_boxcuts_base1530_base540_wp050",
              kSlide33Requested,
              4.7,
              35.3,
              0.00,
              0.62,
              "full range: box-cuts and 5-40 BDT; red: 15-30 BDT",
              "id_efficiency_slide33_boxcuts_base1530_base540_wp050_points.csv",
              0.180,
              true,
              1.15);

  drawOverlay("candidate_et_id_efficiency_slide33_boxcuts_base1530_base540_wp050",
              kSlide33Requested,
              4.7,
              35.3,
              0.00,
              1.04,
              "full range: box-cuts and 5-40 BDT; red: 15-30 BDT",
              "candidate_et_id_efficiency_slide33_boxcuts_base1530_base540_wp050_points.csv",
              0.180,
              true,
              1.15,
              EfficiencyMode::CandidateEt,
              "Photon candidate E_{T} [GeV]");
}
