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
const std::string kBaseDir = "dataOutput/combinedSimOnlyEMBEDDED";
const std::string kMergedRel = "photonJet12and20merged_SIM/RecoilJets_embeddedPhoton12plus20_MERGED.root";
const std::string kOutDir = "dataOutput/auau_bdt_mc_validation/photon_efficiency_qa";

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

struct Point
{
  double x = 0.0;
  double ex = 0.0;
  double y = 0.0;
  double ey = 0.0;
  double tight = 0.0;
  double nonTight = 0.0;
  double tightErr = 0.0;
  double nonTightErr = 0.0;
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

double integralAndError(TH1* h, double& err)
{
  err = 0.0;
  if (!h) return 0.0;
  return h->IntegralAndError(0, h->GetNbinsX() + 1, err);
}

Point buildPoint(TDirectory* d, const PtBin& pt, const CentSpec& cent)
{
  const std::string tag = "isoR30_pT_" + std::to_string(pt.lo) + "_" + std::to_string(pt.hi) + "_cent_" + cent.suffix;
  TH1* hTight = getHist(d, "h_Eiso_tight_" + tag);
  TH1* hNonTight = getHist(d, "h_Eiso_nonTight_" + tag);

  double eTight = 0.0;
  double eNonTight = 0.0;
  const double nTight = integralAndError(hTight, eTight);
  const double nNonTight = integralAndError(hNonTight, eNonTight);
  const double den = nTight + nNonTight;

  Point p;
  p.x = 0.5 * (pt.lo + pt.hi);
  p.ex = 0.0;
  p.tight = nTight;
  p.nonTight = nNonTight;
  p.tightErr = eTight;
  p.nonTightErr = eNonTight;
  if (den > 0.0)
  {
    p.y = nTight / den;
    // Ratio uncertainty using histogram sumw2 errors; conservative for weighted MC.
    const double eDen = std::hypot(eTight, eNonTight);
    const double relT = (nTight > 0.0) ? eTight / nTight : 0.0;
    const double relD = eDen / den;
    p.ey = p.y * std::hypot(relT, relD);
    if (!std::isfinite(p.ey)) p.ey = 0.0;
  }
  return p;
}

std::unique_ptr<TGraphErrors> buildGraph(TDirectory* d,
                                         const CurveSpec& curve,
                                         const CentSpec& cent,
                                         const std::vector<PtBin>& bins,
                                         std::ofstream& csv)
{
  std::unique_ptr<TGraphErrors> g(new TGraphErrors());
  g->SetName(("g_" + curve.label + "_" + cent.suffix).c_str());
  g->SetMarkerStyle(curve.marker);
  g->SetMarkerSize(1.15);
  g->SetMarkerColor(curve.color);
  g->SetLineColor(curve.color);
  g->SetLineWidth(2);

  for (const auto& pt : bins)
  {
    Point p = buildPoint(d, pt, cent);
    if (p.tight + p.nonTight <= 0.0) continue;

    const int ip = g->GetN();
    g->SetPoint(ip, p.x + curve.xShift, p.y);
    g->SetPointError(ip, p.ex, p.ey);

    csv << curve.cfg << "," << curve.label << "," << cent.suffix << ","
        << pt.lo << "," << pt.hi << ","
        << std::setprecision(10) << p.tight << "," << p.nonTight << ","
        << p.y << "," << p.ey << "\n";
  }
  return g;
}

void drawPanelLabel(const char* text)
{
  TLatex tx;
  tx.SetNDC();
  tx.SetTextFont(42);
  tx.SetTextAlign(22);
  tx.SetTextSize(0.052);
  tx.DrawLatex(0.50, 0.955, text);
}

void drawHeader()
{
  TLatex tx;
  tx.SetNDC();
  tx.SetTextFont(42);
  tx.SetTextAlign(13);
  tx.SetTextSize(0.037);
  tx.DrawLatex(0.135, 0.855, "#it{#bf{sPHENIX}} Internal");
  tx.SetTextSize(0.029);
  tx.DrawLatex(0.135, 0.785, "Pythia overlay, #sqrt{s_{NN}} = 200 GeV");
}
}  // namespace

void PlotAuAuBDTIDEfficiency()
{
  SetsPhenixStyle();
  gStyle->SetOptStat(0);
  gSystem->mkdir(kOutDir.c_str(), true);

  const std::vector<PtBin> bins = {
      {5, 8}, {8, 10}, {10, 12}, {12, 14}, {14, 16}, {16, 18},
      {18, 20}, {20, 22}, {22, 24}, {24, 26}, {26, 35}};

  const std::vector<CentSpec> cents = {
      {"0_20", "0-20% central"},
      {"20_50", "20-50% mid-central"},
      {"50_80", "50-80% peripheral"}};

  const std::vector<CurveSpec> curves = {
      {"preselectionNewPPG12_tightReference_nonTightReference_baseVariant",
       "Box-cuts", kBlack, 20, -0.20},
      {"preselectionNewPPG12_tightAuAuCentInputBDT_nonTightAuAuBDTComplement_baseVariant",
       "Centrality as input", kAzure + 2, 20, -0.10},
      {"preselectionNewPPG12_tightAuAuCent3BDT_nonTightAuAuBDTComplement_baseVariant",
       "3 centrality-bin BDTs", kOrange + 7, 20, 0.00},
      {"preselectionNewPPG12_tightAuAuCent7BDT_nonTightAuAuBDTComplement_baseVariant",
       "7 centrality-bin BDTs", kPink + 1, 20, 0.10},
      {"preselectionNewPPG12_tightAuAuPtCent7BDT_nonTightAuAuBDTComplement_baseVariant",
       "E_{T} #times 7 centrality-bin BDTs", kViolet + 2, 20, 0.20}};

  std::ofstream csv(kOutDir + "/id_efficiency_summary_5to35.csv");
  csv << "config_key,label,cent,pt_lo,pt_hi,tight_sum,nontight_sum,id_eff,id_eff_err\n";

  TCanvas c("c_auau_bdt_id_efficiency_1x3",
            "c_auau_bdt_id_efficiency_1x3", 2200, 760);
  c.Divide(3, 1, 0.012, 0.0);

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
    gPad->SetLeftMargin(ic == 0 ? 0.115 : 0.055);
    gPad->SetRightMargin(ic == cents.size() - 1 ? 0.035 : 0.015);
    gPad->SetTopMargin(0.305);
    gPad->SetBottomMargin(0.14);

    std::unique_ptr<TH1F> frame(new TH1F(("hframe_id_eff_" + cents[ic].suffix).c_str(), "", 100, 4.7, 35.3));
    frame->SetDirectory(nullptr);
    frame->SetStats(false);
    frame->SetMinimum(0.0);
    frame->SetMaximum(1.04);
    frame->GetXaxis()->SetTitle("Photon candidate E_{T} [GeV]");
    frame->GetYaxis()->SetTitle(ic == 0 ? "Tight-ID efficiency" : "");
    frame->GetXaxis()->SetTitleSize(0.044);
    frame->GetYaxis()->SetTitleSize(0.044);
    frame->GetXaxis()->SetLabelSize(0.037);
    frame->GetYaxis()->SetLabelSize(0.037);
    frame->GetYaxis()->SetTitleOffset(ic == 0 ? 1.16 : 1.02);
    frame->Draw();
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
      auto g = buildGraph(d, curve, cents[ic], bins, csv);
      if (g && g->GetN() > 0)
      {
        g->Draw("P SAME");
        allGraphs[ic].push_back(std::move(g));
      }
      openFiles.push_back(std::move(f));
    }

    if (ic == 2)
    {
      std::unique_ptr<TLegend> leg(new TLegend(0.15, 0.13, 0.90, 0.43));
      leg->SetBorderSize(0);
      leg->SetFillStyle(1001);
      leg->SetFillColorAlpha(kWhite, 0.86);
      leg->SetTextFont(42);
      leg->SetTextSize(0.034);
      leg->SetNColumns(2);
      for (const auto& curve : curves)
      {
        std::unique_ptr<TGraphErrors> g(new TGraphErrors(1));
        g->SetMarkerStyle(curve.marker);
        g->SetMarkerSize(1.15);
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
  const std::string out = kOutDir + "/id_efficiency_1x3_centrality_table_box_vs_bdt_isoR30_isSliding_5to35.png";
  c.SaveAs(out.c_str());
  std::cout << "[DONE] wrote " << out << "\n";
  std::cout << "[DONE] wrote " << kOutDir << "/id_efficiency_summary_5to35.csv\n";
}
