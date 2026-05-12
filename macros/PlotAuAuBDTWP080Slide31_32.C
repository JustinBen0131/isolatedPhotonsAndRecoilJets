#include "sPhenixStyle.C"

#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH1F.h>
#include <TLatex.h>
#include <TLegend.h>
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
    "dataOutput/auau_bdt_mc_validation/wp080_no3x3_fanout1_clean/combinedSimOnlyEMBEDDED";
const std::string kSignalRel =
    "photonJet12and20merged_SIM/RecoilJets_embeddedPhoton12plus20_MERGED.root";
const std::string kOutDir =
    "dataOutput/auau_bdt_mc_validation/wp080_no3x3_fanout1_clean/photon_efficiency_qa";

struct CurveSpec
{
  std::string cfg;
  std::string label;
  int color = kBlack;
  int marker = 20;
  double xShift = 0.0;
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
  double a = 0.0;
  double b = 0.0;
};

const std::vector<CentSpec> kCents = {
    {"0_20", "0-20% central"},
    {"20_50", "20-50% mid-central"},
    {"50_80", "50-80% peripheral"}};

const std::vector<PtBin> kPtBins = {
    {5, 8}, {8, 10}, {10, 12}, {12, 14}, {14, 16}, {16, 18},
    {18, 20}, {20, 22}, {22, 24}, {24, 26}, {26, 35}};

const std::vector<CurveSpec> kCurves = {
    {"preselectionNewPPG12_tightReference_nonTightReference_baseVariant",
     "Box-cuts", kBlack, 20, -0.24},
    {"preselectionNewPPG12_tightAuAuCentInputBDT_nonTightAuAuBDTComplement_baseVariant",
     "Centrality as input", kAzure + 2, 20, -0.12},
    {"preselectionNewPPG12_tightAuAuCent3BDT_nonTightAuAuBDTComplement_baseVariant",
     "3 centrality-bin BDTs", kOrange + 7, 20, 0.00},
    {"preselectionNewPPG12_tightAuAuCent7BDT_nonTightAuAuBDTComplement_baseVariant",
     "7 centrality-bin BDTs", kPink + 1, 20, 0.12},
    {"preselectionNewPPG12_tightAuAuPtCent7BDT_nonTightAuAuBDTComplement_baseVariant",
     "E_{T} #times 7 centrality-bin BDTs", kViolet + 2, 20, 0.24}};

const std::vector<CurveSpec> kSimpleCurves = {
    {"preselectionNewPPG12_tightReference_nonTightReference_baseVariant",
     "Box-cuts", kBlack, 20, -0.08},
    {"preselectionNewPPG12_tightAuAuCentInputBDT_nonTightAuAuBDTComplement_baseVariant",
     "Centrality as input BDT", kAzure + 2, 20, 0.08}};

std::string signalPath(const CurveSpec& c)
{
  return kBaseDir + "/" + c.cfg + "/" + kSignalRel;
}

TH1* getHist(TDirectory* d, const std::string& name)
{
  if (!d) return nullptr;
  return dynamic_cast<TH1*>(d->Get(name.c_str()));
}

double integralAndError(TH1* h, double& err)
{
  err = 0.0;
  if (!h) return 0.0;
  return h->IntegralAndError(0, h->GetNbinsX() + 1, err);
}

std::vector<EffPoint> readRecoPoints(TDirectory* d, const CentSpec& cent)
{
  std::vector<EffPoint> pts;
  TH1* hTruth = getHist(d, "h_unfoldTruthPho_pTgamma_cent_" + cent.suffix);
  TH1* hMiss = getHist(d, "h_unfoldTruthPhoMisses_pTgamma_isoR30_isSliding_cent_" + cent.suffix);
  if (!hTruth || !hMiss) return pts;

  for (int ib = 1; ib <= hTruth->GetNbinsX(); ++ib)
  {
    const double lo = hTruth->GetXaxis()->GetBinLowEdge(ib);
    const double hi = hTruth->GetXaxis()->GetBinUpEdge(ib);
    if (lo < 5.0 || lo >= 35.0) continue;
    const double mid = 0.5 * (lo + hi);
    const int im = hMiss->GetXaxis()->FindBin(mid);
    if (im < 1 || im > hMiss->GetNbinsX()) continue;

    const double truth = hTruth->GetBinContent(ib);
    const double miss = hMiss->GetBinContent(im);
    if (truth <= 0.0) continue;
    const double eTruth = hTruth->GetBinError(ib);
    const double eMiss = hMiss->GetBinError(im);
    const double r = miss / truth;

    EffPoint p;
    p.x = mid;
    p.y = 1.0 - r;
    p.ey = miss > 0.0 ? r * std::hypot(eMiss / miss, eTruth / truth) : 0.0;
    if (!std::isfinite(p.ey)) p.ey = 0.0;
    p.a = truth;
    p.b = miss;
    pts.push_back(p);
  }
  return pts;
}

EffPoint readIDPoint(TDirectory* d, const PtBin& pt, const CentSpec& cent)
{
  const std::string tag =
      "isoR30_pT_" + std::to_string(pt.lo) + "_" + std::to_string(pt.hi) +
      "_cent_" + cent.suffix;
  double eTight = 0.0;
  double eNonTight = 0.0;
  const double tight = integralAndError(getHist(d, "h_Eiso_tight_" + tag), eTight);
  const double nonTight = integralAndError(getHist(d, "h_Eiso_nonTight_" + tag), eNonTight);
  const double den = tight + nonTight;

  EffPoint p;
  p.x = 0.5 * (pt.lo + pt.hi);
  p.a = tight;
  p.b = nonTight;
  if (den > 0.0)
  {
    p.y = tight / den;
    const double eDen = std::hypot(eTight, eNonTight);
    const double relT = tight > 0.0 ? eTight / tight : 0.0;
    const double relD = eDen / den;
    p.ey = p.y * std::hypot(relT, relD);
    if (!std::isfinite(p.ey)) p.ey = 0.0;
  }
  return p;
}

std::unique_ptr<TGraphErrors> makeGraph(const CurveSpec& curve,
                                        const std::vector<EffPoint>& points,
                                        std::ofstream& csv,
                                        const std::string& quantity,
                                        const std::string& cent)
{
  auto g = std::make_unique<TGraphErrors>();
  g->SetName(("g_" + quantity + "_" + curve.label + "_" + cent).c_str());
  g->SetMarkerStyle(curve.marker);
  g->SetMarkerSize(1.22);
  g->SetMarkerColor(curve.color);
  g->SetLineColor(curve.color);
  g->SetLineWidth(0);
  for (const auto& p : points)
  {
    const int ip = g->GetN();
    g->SetPoint(ip, p.x + curve.xShift, p.y);
    g->SetPointError(ip, 0.0, p.ey);
    csv << quantity << "," << curve.cfg << "," << curve.label << "," << cent
        << "," << std::setprecision(8) << p.x << "," << p.a << "," << p.b
        << "," << p.y << "," << p.ey << "\n";
  }
  return g;
}

void drawPanelLabel(const std::string& text)
{
  TLatex tx;
  tx.SetNDC();
  tx.SetTextFont(42);
  tx.SetTextAlign(22);
  tx.SetTextSize(0.050);
  tx.DrawLatex(0.50, 0.955, text.c_str());
}

void drawHeader(const std::string& extra)
{
  TLatex tx;
  tx.SetNDC();
  tx.SetTextFont(42);
  tx.SetTextAlign(13);
  tx.SetTextSize(0.038);
  tx.DrawLatex(0.145, 0.875, "#it{#bf{sPHENIX}} Internal");
  tx.SetTextSize(0.030);
  tx.DrawLatex(0.145, 0.805, "Pythia overlay, #sqrt{s_{NN}} = 200 GeV");
  tx.DrawLatex(0.145, 0.744, extra.c_str());
}

void drawLegend(const std::vector<CurveSpec>& curves, bool compact)
{
  std::unique_ptr<TLegend> leg(
      compact ? new TLegend(0.58, 0.770, 0.93, 0.900)
              : new TLegend(0.43, 0.745, 0.93, 0.905));
  leg->SetBorderSize(0);
  leg->SetFillStyle(1001);
  leg->SetFillColorAlpha(kWhite, 0.90);
  leg->SetTextFont(42);
  leg->SetTextSize(compact ? 0.031 : 0.028);
  leg->SetNColumns(compact ? 1 : 2);
  static std::vector<std::unique_ptr<TGraphErrors>> keepAlive;
  for (const auto& curve : curves)
  {
    auto g = std::make_unique<TGraphErrors>(1);
    g->SetMarkerStyle(curve.marker);
    g->SetMarkerSize(1.22);
    g->SetMarkerColor(curve.color);
    leg->AddEntry(g.get(), curve.label.c_str(), "p");
    keepAlive.push_back(std::move(g));
  }
  leg->Draw();
  leg.release();
}

void drawOverlay(const std::string& quantity,
                 const std::string& yTitle,
                 double yMin,
                 double yMax,
                 const std::string& xTitle,
                 const std::string& headerExtra,
                 const std::string& outName,
                 const std::vector<CurveSpec>& curves = kCurves,
                 bool compactLegend = false)
{
  TCanvas c(("c_" + quantity + "_wp080").c_str(), ("c_" + quantity + "_wp080").c_str(), 3590, 1159);
  c.Divide(3, 1, 0.006, 0.0);

  std::ofstream csv(kOutDir + "/" + quantity + "_wp080_points.csv");
  csv << "quantity,config_key,label,cent,pt_mid_or_bin_center,numerator_or_truth,denominator_piece_or_miss,value,error\n";

  std::vector<std::unique_ptr<TFile>> openFiles;
  std::vector<std::unique_ptr<TGraphErrors>> graphs;
  std::vector<std::unique_ptr<TH1F>> frames;

  for (size_t ic = 0; ic < kCents.size(); ++ic)
  {
    c.cd(ic + 1);
    gPad->SetTicks(1, 1);
    gPad->SetGrid(0, 0);
    gPad->SetLeftMargin(ic == 0 ? 0.155 : 0.055);
    gPad->SetRightMargin(ic == kCents.size() - 1 ? 0.030 : 0.018);
    gPad->SetTopMargin(0.240);
    gPad->SetBottomMargin(0.14);

    auto frame = std::make_unique<TH1F>(("frame_" + quantity + "_" + kCents[ic].suffix).c_str(), "", 100, 4.7, 35.3);
    frame->SetDirectory(nullptr);
    frame->SetStats(false);
    frame->SetMinimum(yMin);
    frame->SetMaximum(yMax);
    frame->GetXaxis()->SetTitle(xTitle.c_str());
    frame->GetYaxis()->SetTitle(ic == 0 ? yTitle.c_str() : "");
    frame->GetXaxis()->SetTitleSize(0.045);
    frame->GetYaxis()->SetTitleSize(0.043);
    frame->GetXaxis()->SetLabelSize(0.036);
    frame->GetYaxis()->SetLabelSize(0.036);
    frame->GetYaxis()->SetTitleOffset(ic == 0 ? 1.28 : 1.02);
    frame->Draw("AXIS");
    frames.push_back(std::move(frame));

    drawPanelLabel(kCents[ic].title);
    if (ic == 0) drawHeader(headerExtra);

    for (const auto& curve : curves)
    {
      auto f = std::unique_ptr<TFile>(TFile::Open(signalPath(curve).c_str(), "READ"));
      if (!f || f->IsZombie())
      {
        std::cerr << "[WARN] missing file: " << signalPath(curve) << "\n";
        continue;
      }
      TDirectory* d = dynamic_cast<TDirectory*>(f->Get("SIM"));
      if (!d)
      {
        std::cerr << "[WARN] missing SIM directory: " << signalPath(curve) << "\n";
        continue;
      }

      std::vector<EffPoint> pts;
      if (quantity == "reco_efficiency")
        pts = readRecoPoints(d, kCents[ic]);
      else
      {
        for (const auto& pt : kPtBins)
        {
          auto p = readIDPoint(d, pt, kCents[ic]);
          if (p.a + p.b > 0.0) pts.push_back(p);
        }
      }
      if (!pts.empty())
      {
        auto g = makeGraph(curve, pts, csv, quantity, kCents[ic].suffix);
        g->Draw("P SAME");
        graphs.push_back(std::move(g));
      }
      openFiles.push_back(std::move(f));
    }

    if (ic == 2) drawLegend(curves, compactLegend);
  }

  const std::string out = kOutDir + "/" + outName;
  c.SaveAs(out.c_str());
  std::cout << "[DONE] wrote " << out << "\n";
}
}  // namespace

void PlotAuAuBDTWP080Slide31_32()
{
  SetsPhenixStyle();
  gStyle->SetOptStat(0);
  gStyle->SetPadGridX(false);
  gStyle->SetPadGridY(false);
  gSystem->mkdir(kOutDir.c_str(), true);

  drawOverlay("reco_efficiency",
              "Reconstruction efficiency",
              0.0,
              0.64,
              "Truth photon E_{T} [GeV]",
              "R = 0.3, sliding isolation; BDT rows use score > 0.80",
              "reco_efficiency_slide31_style_bdt_score_gt080.png");

  drawOverlay("id_efficiency",
              "Tight-ID efficiency",
              0.0,
              1.04,
              "Photon candidate E_{T} [GeV]",
              "R = 0.3, sliding isolation; BDT rows use score > 0.80",
              "id_efficiency_slide32_style_bdt_score_gt080.png");

  drawOverlay("id_efficiency_box_vs_centinput",
              "Tight-ID efficiency",
              0.0,
              1.04,
              "Photon candidate E_{T} [GeV]",
              "R = 0.3, sliding isolation; centrality-input BDT uses score > 0.80",
              "id_efficiency_slide32_box_vs_centinput_bdt_score_gt080.png",
              kSimpleCurves,
              true);
}
