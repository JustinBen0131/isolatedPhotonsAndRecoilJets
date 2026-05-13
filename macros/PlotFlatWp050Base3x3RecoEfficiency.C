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
const std::string kInputBase =
    "dataOutput/auau_widthstudy_pt1530_wp050/combinedSimOnlyEMBEDDED";
const std::string kOutputDir =
    "dataOutput/target80_first_offline/bdt_target80_gated_20260512_001012/"
    "analysis_config_etfine_15to35_target80/efficiency_qa";
const std::string kSignalRel =
    "photonJet12and20merged_SIM/RecoilJets_embeddedPhoton12plus20_MERGED.root";

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

struct Point
{
  double x = 0.0;
  double y = 0.0;
  double ey = 0.0;
  double truth = 0.0;
  double miss = 0.0;
};

const std::vector<CentSpec> kCents = {
    {"0_20", "0-20% central"},
    {"20_50", "20-50% mid-central"},
    {"50_80", "50-80% peripheral"}};

const std::vector<CurveSpec> kCurves = {
    {"preselectionNewPPG12_tightReference_nonTightReference_baseVariant",
     "Box-cuts", kBlack, 20, -0.10},
    {"preselectionNewPPG12_tightAuauCentInputBase3x3BDT_nonTightAuAuBDTComplement_baseVariant",
     "BDT > 0.5, base + 3x3 widths", kAzure + 2, 20, 0.10}};

TH1* getHist(TDirectory* d, const std::string& name)
{
  if (!d) return nullptr;
  TH1* h = dynamic_cast<TH1*>(d->Get(name.c_str()));
  if (h) h->SetStats(false);
  return h;
}

std::string sigPath(const CurveSpec& curve)
{
  return kInputBase + "/" + curve.cfg + "/" + kSignalRel;
}

std::vector<Point> readRecoPoints(TDirectory* d, const CentSpec& cent)
{
  std::vector<Point> pts;
  TH1* hTruth = getHist(d, "h_unfoldTruthPho_pTgamma_cent_" + cent.suffix);
  TH1* hMiss = getHist(d, "h_unfoldTruthPhoMisses_pTgamma_isoR30_isSliding_cent_" + cent.suffix);
  if (!hTruth || !hMiss)
  {
    std::cerr << "[WARN] missing truth/miss histogram for centrality " << cent.suffix << "\n";
    return pts;
  }

  for (int ib = 1; ib <= hTruth->GetNbinsX(); ++ib)
  {
    const double lo = hTruth->GetXaxis()->GetBinLowEdge(ib);
    const double hi = hTruth->GetXaxis()->GetBinUpEdge(ib);
    if (lo < 15.0 || hi > 30.0) continue;
    const double mid = 0.5 * (lo + hi);
    const int im = hMiss->GetXaxis()->FindBin(mid);
    if (im < 1 || im > hMiss->GetNbinsX()) continue;

    const double truth = hTruth->GetBinContent(ib);
    const double miss = hMiss->GetBinContent(im);
    if (truth <= 0.0) continue;
    const double eTruth = hTruth->GetBinError(ib);
    const double eMiss = hMiss->GetBinError(im);
    const double missFrac = miss / truth;

    Point p;
    p.x = mid;
    p.y = 1.0 - missFrac;
    p.ey = miss > 0.0 ? missFrac * std::hypot(eMiss / miss, eTruth / truth) : 0.0;
    if (!std::isfinite(p.ey)) p.ey = 0.0;
    p.truth = truth;
    p.miss = miss;
    pts.push_back(p);
  }
  return pts;
}

std::unique_ptr<TGraphErrors> makeGraph(const CurveSpec& curve,
                                        const std::vector<Point>& pts,
                                        std::ofstream& csv,
                                        const std::string& cent)
{
  auto g = std::make_unique<TGraphErrors>();
  g->SetMarkerStyle(curve.marker);
  g->SetMarkerSize(1.25);
  g->SetMarkerColor(curve.color);
  g->SetLineColor(curve.color);
  g->SetLineWidth(2);
  for (const auto& p : pts)
  {
    const int ip = g->GetN();
    g->SetPoint(ip, p.x + curve.xShift, p.y);
    g->SetPointError(ip, 0.0, p.ey);
    csv << curve.cfg << "," << curve.label << "," << cent << ","
        << std::setprecision(9) << p.x << "," << p.truth << "," << p.miss
        << "," << p.y << "," << p.ey << "\n";
  }
  return g;
}

void drawHeader()
{
  TLatex tx;
  tx.SetNDC();
  tx.SetTextFont(42);
  tx.SetTextAlign(13);
  tx.SetTextSize(0.040);
  tx.DrawLatex(0.145, 0.875, "#it{#bf{sPHENIX}} Internal");
  tx.SetTextSize(0.031);
  tx.DrawLatex(0.145, 0.815, "Pythia overlay, #sqrt{s_{NN}} = 200 GeV");
  tx.DrawLatex(0.145, 0.760, "15 #leq E_{T}^{cluster} < 30 GeV; R = 0.3, sliding isolation");
}

void drawPanelLabel(const std::string& text)
{
  TLatex tx;
  tx.SetNDC();
  tx.SetTextFont(42);
  tx.SetTextAlign(22);
  tx.SetTextSize(0.052);
  tx.DrawLatex(0.50, 0.950, text.c_str());
}

void drawLegend()
{
  auto leg = new TLegend(0.42, 0.775, 0.94, 0.925);
  leg->SetBorderSize(0);
  leg->SetFillStyle(1001);
  leg->SetFillColorAlpha(kWhite, 0.92);
  leg->SetTextFont(42);
  leg->SetTextSize(0.031);
  static std::vector<std::unique_ptr<TGraphErrors>> keepAlive;
  for (const auto& curve : kCurves)
  {
    auto g = std::make_unique<TGraphErrors>(1);
    g->SetMarkerStyle(curve.marker);
    g->SetMarkerSize(1.25);
    g->SetMarkerColor(curve.color);
    leg->AddEntry(g.get(), curve.label.c_str(), "p");
    keepAlive.push_back(std::move(g));
  }
  leg->Draw();
}
}  // namespace

void PlotFlatWp050Base3x3RecoEfficiency()
{
  SetsPhenixStyle();
  gStyle->SetOptStat(0);
  gStyle->SetPadGridX(false);
  gStyle->SetPadGridY(false);
  gSystem->mkdir(kOutputDir.c_str(), true);

  std::ofstream csv(kOutputDir + "/wp050_base3x3_pt1530_reco_efficiency_points.csv");
  csv << "config_key,label,cent,pt_mid,n_truth,n_miss,recovery,error\n";

  TCanvas c("c_wp050_base3x3_reco", "c_wp050_base3x3_reco", 3590, 1159);
  c.Divide(3, 1, 0.008, 0.0);

  std::vector<std::unique_ptr<TFile>> files;
  std::vector<std::unique_ptr<TGraphErrors>> graphs;
  std::vector<std::unique_ptr<TH1F>> frames;

  for (size_t ic = 0; ic < kCents.size(); ++ic)
  {
    c.cd(ic + 1);
    gPad->SetTicks(1, 1);
    gPad->SetGrid(0, 0);
    gPad->SetLeftMargin(ic == 0 ? 0.155 : 0.055);
    gPad->SetRightMargin(ic == kCents.size() - 1 ? 0.030 : 0.018);
    gPad->SetTopMargin(0.320);
    gPad->SetBottomMargin(0.14);

    auto frame = std::make_unique<TH1F>(("frame_wp050_" + kCents[ic].suffix).c_str(), "", 100, 14.4, 30.6);
    frame->SetDirectory(nullptr);
    frame->SetStats(false);
    frame->SetMinimum(0.0);
    frame->SetMaximum(0.70);
    frame->GetXaxis()->SetTitle("Truth photon E_{T} [GeV]");
    frame->GetYaxis()->SetTitle(ic == 0 ? "Truth-photon recovery" : "");
    frame->GetXaxis()->SetTitleSize(0.045);
    frame->GetYaxis()->SetTitleSize(0.043);
    frame->GetXaxis()->SetLabelSize(0.036);
    frame->GetYaxis()->SetLabelSize(0.036);
    frame->GetYaxis()->SetTitleOffset(ic == 0 ? 1.28 : 1.02);
    frame->Draw("AXIS");
    frames.push_back(std::move(frame));

    drawPanelLabel(kCents[ic].title);
    if (ic == 0) drawHeader();

    for (const auto& curve : kCurves)
    {
      const std::string path = sigPath(curve);
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
      auto pts = readRecoPoints(d, kCents[ic]);
      auto g = makeGraph(curve, pts, csv, kCents[ic].suffix);
      g->Draw("P SAME");
      graphs.push_back(std::move(g));
      files.push_back(std::move(f));
    }
    if (ic == 2) drawLegend();
  }

  const std::string out = kOutputDir + "/wp050_base3x3_pt1530_reco_efficiency_1x3.png";
  c.SaveAs(out.c_str());
  std::cout << "[DONE] wrote " << out << "\n";
}
