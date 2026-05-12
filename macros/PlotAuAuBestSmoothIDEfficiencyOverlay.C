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
const std::string kBaseDir = "dataOutput/combinedSimOnlyEMBEDDED";
const std::string kMergedRel =
    "photonJet12and20merged_SIM/RecoilJets_embeddedPhoton12plus20_MERGED.root";
const std::string kOutDir =
    "dataOutput/jstg_slide_candidates/best_smooth_id_efficiency";

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

const std::vector<CentSpec> kCents = {
    {"0_20", "0-20% central"},
    {"20_50", "20-50% mid-central"},
    {"50_80", "50-80% peripheral"}};

const std::vector<CurveSpec> kCurves = {
    {"preselectionNewPPG12_tightReference_nonTightReference_baseVariant",
     "box-cuts", kBlack, -0.14},
    {"preselectionNewPPG12_tightAuAuCentInputBDT_nonTightAuAuBDTComplement_baseVariant",
     "centrality as input BDT", kAzure + 2, 0.00},
    {"preselectionNewPPG12_tightAuAuPtCent7BDT_nonTightAuAuBDTComplement_baseVariant",
     "E_{T} #times 7 centrality-bin BDT", kViolet + 2, 0.14}};

std::string mergedPath(const CurveSpec& curve)
{
  return kBaseDir + "/" + curve.cfg + "/" + kMergedRel;
}

TH1* getHist(TDirectory* d, const std::string& name)
{
  if (!d) return nullptr;
  TH1* h = dynamic_cast<TH1*>(d->Get(name.c_str()));
  if (h) h->SetStats(false);
  return h;
}

std::vector<EffPoint> readPoints(TDirectory* d, const CentSpec& cent)
{
  std::vector<EffPoint> pts;
  TH1* hTruth = getHist(d, "h_unfoldTruthPho_pTgamma_cent_" + cent.suffix);
  TH1* hMiss = getHist(d, "h_unfoldTruthPhoMisses_pTgamma_isoR30_isSliding_cent_" + cent.suffix);
  if (!hTruth || !hMiss)
  {
    std::cerr << "[WARN] missing truth/reco efficiency histograms for " << cent.suffix << "\n";
    return pts;
  }

  for (int ib = 1; ib <= hTruth->GetNbinsX(); ++ib)
  {
    const double lo = hTruth->GetXaxis()->GetBinLowEdge(ib);
    const double hi = hTruth->GetXaxis()->GetBinUpEdge(ib);
    if (lo < 5.0 || hi > 35.0) continue;

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
    pts.push_back(p);
  }
  return pts;
}

std::unique_ptr<TGraphErrors> makeGraph(TDirectory* d,
                                        const CurveSpec& curve,
                                        const CentSpec& cent,
                                        std::ofstream& csv)
{
  const auto pts = readPoints(d, cent);
  if (pts.empty()) return nullptr;

  auto g = std::make_unique<TGraphErrors>();
  g->SetName(("g_best_eff_" + cent.suffix + "_" + curve.label).c_str());
  g->SetMarkerStyle(20);
  g->SetMarkerSize(1.55);
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

void drawHeader()
{
  TLatex tx;
  tx.SetNDC();
  tx.SetTextFont(42);
  tx.SetTextAlign(13);
  tx.SetTextSize(0.038);
  tx.DrawLatex(0.040, 0.980, "#it{#bf{sPHENIX}} Internal");
  tx.SetTextSize(0.029);
  tx.DrawLatex(0.040, 0.935, "Pythia overlay, #sqrt{s_{NN}} = 200 GeV");
}

void drawPanelLabel(const std::string& title)
{
  TLatex tx;
  tx.SetNDC();
  tx.SetTextFont(42);
  tx.SetTextAlign(22);
  tx.SetTextSize(0.052);
  tx.DrawLatex(0.50, 0.925, title.c_str());
}

void drawGlobalLegend()
{
  auto leg = new TLegend(0.295, 0.905, 0.965, 0.992);
  leg->SetBorderSize(0);
  leg->SetFillStyle(1001);
  leg->SetFillColorAlpha(kWhite, 0.94);
  leg->SetTextFont(42);
  leg->SetTextSize(0.034);
  leg->SetNColumns(3);

  static std::vector<std::unique_ptr<TGraphErrors>> keepAlive;
  keepAlive.clear();
  for (const auto& curve : kCurves)
  {
    auto g = std::make_unique<TGraphErrors>(1);
    g->SetMarkerStyle(20);
    g->SetMarkerSize(1.75);
    g->SetMarkerColor(curve.color);
    g->SetLineColor(curve.color);
    leg->AddEntry(g.get(), curve.label.c_str(), "pe");
    keepAlive.push_back(std::move(g));
  }
  leg->Draw();
}
}  // namespace

void PlotAuAuBestSmoothIDEfficiencyOverlay()
{
  SetsPhenixStyle();
  gStyle->SetOptStat(0);
  gStyle->SetPadGridX(false);
  gStyle->SetPadGridY(false);
  gSystem->mkdir(kOutDir.c_str(), true);

  std::ofstream csv(kOutDir + "/best_smooth_id_efficiency_overlay_points.csv");
  csv << "config_key,label,cent,pt_mid,truth_sum,miss_sum,reco_eff,reco_eff_err\n";

  TCanvas c("c_best_smooth_id_efficiency_overlay",
            "c_best_smooth_id_efficiency_overlay", 3590, 1159);

  std::vector<std::unique_ptr<TPad>> pads;
  std::vector<std::unique_ptr<TFile>> files;
  std::vector<std::unique_ptr<TGraphErrors>> graphs;
  std::vector<std::unique_ptr<TH1F>> frames;

  for (size_t ic = 0; ic < kCents.size(); ++ic)
  {
    const double x0 = static_cast<double>(ic) / 3.0 + (ic == 0 ? 0.030 : 0.012);
    const double x1 = static_cast<double>(ic + 1) / 3.0 - (ic == 2 ? 0.018 : 0.010);
    auto pad = std::make_unique<TPad>(("pad_best_" + kCents[ic].suffix).c_str(), "", x0, 0.030, x1, 0.890);
    pad->SetFillStyle(4000);
    pad->SetFrameFillStyle(4000);
    c.cd();
    pad->Draw();
    pad->cd();
    pad->SetTicks(1, 1);
    pad->SetGrid(0, 0);
    pad->SetLeftMargin(ic == 0 ? 0.155 : 0.060);
    pad->SetRightMargin(ic == kCents.size() - 1 ? 0.035 : 0.020);
    pad->SetTopMargin(0.170);
    pad->SetBottomMargin(0.145);

    auto frame = std::make_unique<TH1F>(("hframe_best_" + kCents[ic].suffix).c_str(), "", 100, 4.7, 35.3);
    frame->SetDirectory(nullptr);
    frame->SetStats(false);
    frame->SetMinimum(0.0);
    frame->SetMaximum(0.64);
    frame->GetXaxis()->SetTitle("Truth photon E_{T} [GeV]");
    frame->GetYaxis()->SetTitle(ic == 0 ? "Truth-photon recovery efficiency" : "");
    frame->GetXaxis()->SetTitleSize(0.046);
    frame->GetYaxis()->SetTitleSize(0.044);
    frame->GetXaxis()->SetLabelSize(0.038);
    frame->GetYaxis()->SetLabelSize(ic == 0 ? 0.038 : 0.0);
    frame->GetYaxis()->SetTitleOffset(ic == 0 ? 1.28 : 1.0);
    frame->GetXaxis()->SetNdivisions(506);
    frame->GetYaxis()->SetNdivisions(506);
    frame->Draw("AXIS");
    frames.push_back(std::move(frame));

    drawPanelLabel(kCents[ic].title);

    for (const auto& curve : kCurves)
    {
      auto f = std::unique_ptr<TFile>(TFile::Open(mergedPath(curve).c_str(), "READ"));
      if (!f || f->IsZombie())
      {
        std::cerr << "[WARN] missing file: " << mergedPath(curve) << "\n";
        continue;
      }
      TDirectory* d = dynamic_cast<TDirectory*>(f->Get("SIM"));
      if (!d)
      {
        std::cerr << "[WARN] missing SIM directory: " << mergedPath(curve) << "\n";
        continue;
      }
      auto g = makeGraph(d, curve, kCents[ic], csv);
      if (g && g->GetN() > 0)
      {
        g->Draw("P SAME");
        graphs.push_back(std::move(g));
      }
      files.push_back(std::move(f));
    }
    pad->RedrawAxis();
    pads.push_back(std::move(pad));
  }

  c.cd();
  drawHeader();
  drawGlobalLegend();

  const std::string out = kOutDir + "/best_smooth_id_efficiency_box_centinput_ptcent7_overlay.png";
  c.SaveAs(out.c_str());
  std::cout << "[DONE] wrote " << out << "\n";
  std::cout << "[DONE] wrote " << kOutDir << "/best_smooth_id_efficiency_overlay_points.csv\n";
}
