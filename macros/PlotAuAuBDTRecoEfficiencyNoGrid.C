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
const std::string kBaseDir = "dataOutput/combinedSimOnlyEMBEDDED";
const std::string kMergedRel = "photonJet12and20merged_SIM/RecoilJets_embeddedPhoton12plus20_MERGED.root";
const std::string kOutDir = "dataOutput/auau_bdt_mc_validation/photon_efficiency_qa";

struct CurveSpec
{
  std::string cfg;
  std::string label;
  int color = kBlack;
  double xShift = 0.0;
};

struct CentSpec
{
  std::string suffix;
  std::string title;
};

struct EffPoint
{
  double x = 0.0;
  double y = 0.0;
  double ey = 0.0;
  double truth = 0.0;
  double miss = 0.0;
};

std::string mergedPath(const CurveSpec& c)
{
  return kBaseDir + "/" + c.cfg + "/" + kMergedRel;
}

TH1* getHist(TDirectory* d, const std::string& name)
{
  if (!d) return nullptr;
  return dynamic_cast<TH1*>(d->Get(name.c_str()));
}

std::vector<EffPoint> readPoints(TDirectory* d, const CentSpec& cent)
{
  std::vector<EffPoint> pts;
  TH1* hTruth = getHist(d, "h_unfoldTruthPho_pTgamma_cent_" + cent.suffix);
  TH1* hMiss = getHist(d, "h_unfoldTruthPhoMisses_pTgamma_isoR30_isSliding_cent_" + cent.suffix);
  if (!hTruth || !hMiss) return pts;

  const int nb = hTruth->GetNbinsX();
  for (int ib = 1; ib <= nb; ++ib)
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
    p.truth = truth;
    p.miss = miss;
    p.y = 1.0 - r;
    p.ey = (miss > 0.0) ? r * std::hypot(eMiss / miss, eTruth / truth) : 0.0;
    if (!std::isfinite(p.ey)) p.ey = 0.0;
    pts.push_back(p);
  }
  return pts;
}

std::unique_ptr<TGraphErrors> buildGraph(TDirectory* d,
                                         const CurveSpec& curve,
                                         const CentSpec& cent,
                                         std::ofstream& csv)
{
  const auto pts = readPoints(d, cent);
  if (pts.empty()) return nullptr;

  std::unique_ptr<TGraphErrors> g(new TGraphErrors());
  g->SetName(("g_reco_eff_" + curve.label + "_" + cent.suffix).c_str());
  g->SetMarkerStyle(24);
  g->SetMarkerSize(1.30);
  g->SetMarkerColor(curve.color);
  g->SetLineColor(curve.color);
  g->SetLineWidth(2);

  for (const auto& p : pts)
  {
    const int ip = g->GetN();
    g->SetPoint(ip, p.x + curve.xShift, p.y);
    g->SetPointError(ip, 0.0, p.ey);
    csv << curve.cfg << "," << curve.label << "," << cent.suffix << ","
        << std::setprecision(8) << p.x << "," << p.truth << "," << p.miss
        << "," << p.y << "," << p.ey << "\n";
  }
  return g;
}

void drawPanelLabel(const char* text)
{
  TLatex tx;
  tx.SetNDC();
  tx.SetTextFont(42);
  tx.SetTextAlign(22);
  tx.SetTextSize(0.050);
  tx.DrawLatex(0.50, 0.955, text);
}

void drawHeader()
{
  TLatex tx;
  tx.SetNDC();
  tx.SetTextFont(42);
  tx.SetTextAlign(13);
  tx.SetTextSize(0.041);
  tx.DrawLatex(0.135, 0.815, "#it{#bf{sPHENIX}} Internal");
  tx.SetTextSize(0.032);
  tx.DrawLatex(0.135, 0.745, "Pythia overlay, #sqrt{s_{NN}} = 200 GeV");
  tx.DrawLatex(0.135, 0.680, "R = 0.3, sliding isolation");
}
}  // namespace

void PlotAuAuBDTRecoEfficiencyNoGrid()
{
  SetsPhenixStyle();
  gStyle->SetOptStat(0);
  gStyle->SetPadGridX(false);
  gStyle->SetPadGridY(false);
  gSystem->mkdir(kOutDir.c_str(), true);

  const std::vector<CentSpec> cents = {
      {"0_20", "0-20% central"},
      {"20_50", "20-50% mid-central"},
      {"50_80", "50-80% peripheral"}};

  const std::vector<CurveSpec> curves = {
      {"preselectionNewPPG12_tightReference_nonTightReference_baseVariant",
       "Box-cuts", kBlack, -0.24},
      {"preselectionNewPPG12_tightAuAuCentInputBDT_nonTightAuAuBDTComplement_baseVariant",
       "Cent. as input", kAzure + 2, -0.12},
      {"preselectionNewPPG12_tightAuAuCent3BDT_nonTightAuAuBDTComplement_baseVariant",
       "3 cent.-bin BDTs", kOrange + 7, 0.00},
      {"preselectionNewPPG12_tightAuAuCent7BDT_nonTightAuAuBDTComplement_baseVariant",
       "7 cent.-bin BDTs", kPink + 1, 0.12},
      {"preselectionNewPPG12_tightAuAuPtCent7BDT_nonTightAuAuBDTComplement_baseVariant",
       "E_{T} #times 7 cent. BDTs", kViolet + 2, 0.24}};

  std::ofstream csv(kOutDir + "/reco_efficiency_1x3_nogrid_points.csv");
  csv << "config_key,label,cent,pt_mid,truth_sum,miss_sum,reco_eff,reco_eff_err\n";

  TCanvas c("c_auau_bdt_reco_efficiency_1x3_nogrid",
            "c_auau_bdt_reco_efficiency_1x3_nogrid", 3590, 1159);
  c.Divide(3, 1, 0.014, 0.0);

  std::vector<std::unique_ptr<TFile>> openFiles;
  std::vector<std::vector<std::unique_ptr<TGraphErrors>>> allGraphs(cents.size());
  std::vector<std::unique_ptr<TGraphErrors>> legendKeepAlive;
  std::vector<std::unique_ptr<TLegend>> legends;
  std::vector<std::unique_ptr<TH1F>> frames;

  for (size_t ic = 0; ic < cents.size(); ++ic)
  {
    c.cd(ic + 1);
    gPad->SetTicks(1, 1);
    gPad->SetGrid(0, 0);
    gPad->SetLeftMargin(ic == 0 ? 0.155 : 0.055);
    gPad->SetRightMargin(ic == cents.size() - 1 ? 0.030 : 0.018);
    gPad->SetTopMargin(0.145);
    gPad->SetBottomMargin(0.14);

    std::unique_ptr<TH1F> frame(new TH1F(("hframe_reco_eff_nogrid_" + cents[ic].suffix).c_str(), "", 100, 4.7, 35.3));
    frame->SetDirectory(nullptr);
    frame->SetStats(false);
    frame->SetMinimum(0.0);
    frame->SetMaximum(0.64);
    frame->GetXaxis()->SetTitle("Truth photon E_{T} [GeV]");
    frame->GetYaxis()->SetTitle(ic == 0 ? "Reconstruction efficiency" : "");
    frame->GetXaxis()->SetTitleSize(0.045);
    frame->GetYaxis()->SetTitleSize(0.043);
    frame->GetXaxis()->SetLabelSize(0.036);
    frame->GetYaxis()->SetLabelSize(0.036);
    frame->GetYaxis()->SetTitleOffset(ic == 0 ? 1.28 : 1.02);
    frame->Draw("AXIS");
    frames.push_back(std::move(frame));

    drawPanelLabel(cents[ic].title.c_str());
    if (ic == 0) drawHeader();

    for (const auto& curve : curves)
    {
      const std::string path = mergedPath(curve);
      std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
      if (!f || f->IsZombie())
      {
        std::cerr << "[WARN] missing file: " << path << "\n";
        continue;
      }
      TDirectory* d = dynamic_cast<TDirectory*>(f->Get("SIM"));
      if (!d)
      {
        std::cerr << "[WARN] missing SIM directory in " << path << "\n";
        continue;
      }
      auto g = buildGraph(d, curve, cents[ic], csv);
      if (g && g->GetN() > 0)
      {
        g->Draw("P SAME");
        allGraphs[ic].push_back(std::move(g));
      }
      openFiles.push_back(std::move(f));
    }

    if (ic == 2)
    {
      std::unique_ptr<TLegend> leg(new TLegend(0.16, 0.13, 0.72, 0.30));
      leg->SetBorderSize(0);
      leg->SetFillStyle(1001);
      leg->SetFillColorAlpha(kWhite, 0.88);
      leg->SetTextFont(42);
      leg->SetTextSize(0.031);
      leg->SetNColumns(2);
      for (const auto& curve : curves)
      {
        std::unique_ptr<TGraphErrors> g(new TGraphErrors(1));
        g->SetMarkerStyle(24);
        g->SetMarkerSize(1.30);
        g->SetMarkerColor(curve.color);
        g->SetLineColor(curve.color);
        leg->AddEntry(g.get(), curve.label.c_str(), "pe");
        legendKeepAlive.push_back(std::move(g));
      }
      leg->Draw();
      legends.push_back(std::move(leg));
    }
  }

  c.cd();
  const std::string out = kOutDir + "/reco_efficiency_1x3_centrality_table_box_vs_bdt_enhanced_biglegend_isoR30_isSliding_5to35_nogrid.png";
  c.SaveAs(out.c_str());
  std::cout << "[DONE] wrote " << out << "\n";
  std::cout << "[DONE] wrote " << kOutDir << "/reco_efficiency_1x3_nogrid_points.csv\n";
}
