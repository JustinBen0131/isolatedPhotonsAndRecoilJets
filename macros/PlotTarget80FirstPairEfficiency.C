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
#include <cctype>
#include <string>
#include <vector>

namespace
{
const std::string kBaseDir =
    "dataOutput/target80_first_offline/bdt_target80_gated_20260512_001012/"
    "analysis_config_etfine_15to35_target80";
const std::string kOutDir =
    "dataOutput/target80_first_offline/bdt_target80_gated_20260512_001012/"
    "analysis_config_etfine_15to35_target80/efficiency_qa";
const std::string kSignalRel =
    "photonJet12and20merged_SIM/RecoilJets_embeddedPhoton12plus20_MERGED.root";
const std::string kBkgRel =
    "embeddedJet12and20merged_SIM/RecoilJets_embeddedJet12plus20_MERGED.root";

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
    {15, 17}, {17, 19}, {19, 21}, {21, 23}, {23, 25}, {25, 27}, {27, 30}, {30, 35}};

const std::vector<CurveSpec> kCurves = {
    {"preselectionNewPPG12_tightReference_nonTightReference_baseVariant",
     "Box-cuts", kBlack, 20, -0.12},
    {"preselectionNewPPG12_tightAuAuCentInputBase3x3BDT_nonTightAuAuBDTComplement_baseVariant",
     "Au+Au BDT target-80", kAzure + 2, 20, 0.12}};

const std::vector<CurveSpec> kReferenceVsBdtCurves = {
    {"preselectionNewPPG12_tightReference_nonTightReference_baseVariant",
     "Reference box-cuts", kBlack, 20, -0.08},
    {"preselectionNewPPG12_tightAuAuCentInputBase3x3BDT_nonTightAuAuBDTComplement_baseVariant",
     "Au+Au BDT target-80", kAzure + 2, 20, 0.08}};

std::string sigPath(const CurveSpec& c) { return kBaseDir + "/simembedded/" + c.cfg + "/" + kSignalRel; }
std::string bkgPath(const CurveSpec& c) { return kBaseDir + "/simembeddedinclusive/" + c.cfg + "/" + kBkgRel; }

std::string slugLabel(const std::string& label)
{
  std::string out;
  for (char ch : label)
  {
    if (std::isalnum(static_cast<unsigned char>(ch))) out.push_back(std::tolower(static_cast<unsigned char>(ch)));
    else if (!out.empty() && out.back() != '_') out.push_back('_');
  }
  if (!out.empty() && out.back() == '_') out.pop_back();
  return out;
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

Point ratioPoint(double num, double denPart, double eNum, double eDenPart, double x)
{
  Point p;
  p.x = x;
  p.a = num;
  p.b = denPart;
  const double den = num + denPart;
  if (den <= 0.0) return p;
  p.y = num / den;
  p.ey = std::sqrt(denPart * denPart * eNum * eNum + num * num * eDenPart * eDenPart) / (den * den);
  if (!std::isfinite(p.ey)) p.ey = 0.0;
  return p;
}

std::vector<Point> readRecoPoints(TDirectory* d, const CentSpec& cent)
{
  std::vector<Point> pts;
  TH1* hTruth = getHist(d, "h_unfoldTruthPho_pTgamma_cent_" + cent.suffix);
  TH1* hMiss = getHist(d, "h_unfoldTruthPhoMisses_pTgamma_isoR30_isSliding_cent_" + cent.suffix);
  if (!hTruth || !hMiss) return pts;

  for (int ib = 1; ib <= hTruth->GetNbinsX(); ++ib)
  {
    const double lo = hTruth->GetXaxis()->GetBinLowEdge(ib);
    const double hi = hTruth->GetXaxis()->GetBinUpEdge(ib);
    if (lo < 15.0 || hi > 35.0) continue;
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
    p.a = truth;
    p.b = miss;
    pts.push_back(p);
  }
  return pts;
}

Point readTightFraction(TDirectory* d, const PtBin& pt, const CentSpec& cent)
{
  const std::string tag =
      "isoR30_pT_" + std::to_string(pt.lo) + "_" + std::to_string(pt.hi) +
      "_cent_" + cent.suffix;
  double eT = 0.0;
  double eN = 0.0;
  const double tight = integralAndError(getHist(d, "h_Eiso_tight_" + tag), eT);
  const double nonTight = integralAndError(getHist(d, "h_Eiso_nonTight_" + tag), eN);
  return ratioPoint(tight, nonTight, eT, eN, 0.5 * (pt.lo + pt.hi));
}

std::vector<Point> readTightFractionPoints(TDirectory* d, const CentSpec& cent)
{
  std::vector<Point> pts;
  for (const auto& pt : kPtBins)
  {
    auto p = readTightFraction(d, pt, cent);
    if (p.a + p.b > 0.0) pts.push_back(p);
  }
  return pts;
}

std::unique_ptr<TGraphErrors> makeGraph(const CurveSpec& curve,
                                        const std::vector<Point>& points,
                                        std::ofstream& csv,
                                        const std::string& quantity,
                                        const std::string& cent)
{
  auto g = std::make_unique<TGraphErrors>();
  g->SetMarkerStyle(curve.marker);
  g->SetMarkerSize(1.30);
  g->SetMarkerColor(curve.color);
  g->SetLineColor(curve.color);
  g->SetLineWidth(2);
  for (const auto& p : points)
  {
    const int ip = g->GetN();
    g->SetPoint(ip, p.x + curve.xShift, p.y);
    g->SetPointError(ip, 0.0, p.ey);
    csv << quantity << "," << curve.cfg << "," << curve.label << "," << cent
        << "," << std::setprecision(9) << p.x << "," << p.a << "," << p.b
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
  tx.SetTextSize(0.052);
  tx.DrawLatex(0.50, 0.950, text.c_str());
}

void drawHeader(const std::string& extra)
{
  TLatex tx;
  tx.SetNDC();
  tx.SetTextFont(42);
  tx.SetTextAlign(13);
  tx.SetTextSize(0.040);
  tx.DrawLatex(0.145, 0.865, "#it{#bf{sPHENIX}} Internal");
  tx.SetTextSize(0.031);
  tx.DrawLatex(0.145, 0.795, "Pythia overlay, #sqrt{s_{NN}} = 200 GeV");
  tx.DrawLatex(0.145, 0.735, extra.c_str());
}

void drawLegend()
{
  auto leg = new TLegend(0.48, 0.775, 0.93, 0.925);
  leg->SetBorderSize(0);
  leg->SetFillStyle(1001);
  leg->SetFillColorAlpha(kWhite, 0.92);
  leg->SetTextFont(42);
  leg->SetTextSize(0.031);
  leg->SetNColumns(1);
  static std::vector<std::unique_ptr<TGraphErrors>> keepAlive;
  for (const auto& curve : kCurves)
  {
    auto g = std::make_unique<TGraphErrors>(1);
    g->SetMarkerStyle(curve.marker);
    g->SetMarkerSize(1.30);
    g->SetMarkerColor(curve.color);
    leg->AddEntry(g.get(), curve.label.c_str(), "p");
    keepAlive.push_back(std::move(g));
  }
  leg->Draw();
}

void drawThreePanel(const std::string& quantity,
                    const std::string& yTitle,
                    double yMin,
                    double yMax,
                    const std::string& xTitle,
                    const std::string& headerExtra,
                    const std::string& outName,
                    bool useBackground,
                    bool useReco)
{
  TCanvas c(("c_" + quantity).c_str(), ("c_" + quantity).c_str(), 3590, 1159);
  c.Divide(3, 1, 0.008, 0.0);

  std::ofstream csv(kOutDir + "/" + quantity + "_points.csv");
  csv << "quantity,config_key,label,cent,pt_mid,numerator_or_truth,denominator_piece_or_miss,value,error\n";

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
    gPad->SetTopMargin(0.235);
    gPad->SetBottomMargin(0.14);

    auto frame = std::make_unique<TH1F>(("frame_" + quantity + "_" + kCents[ic].suffix).c_str(), "", 100, 14.4, 35.6);
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

    for (const auto& curve : kCurves)
    {
      const std::string path = useBackground ? bkgPath(curve) : sigPath(curve);
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
      auto pts = useReco ? readRecoPoints(d, kCents[ic]) : readTightFractionPoints(d, kCents[ic]);
      if (!pts.empty())
      {
        auto g = makeGraph(curve, pts, csv, quantity, kCents[ic].suffix);
        g->Draw("P SAME");
        graphs.push_back(std::move(g));
      }
      files.push_back(std::move(f));
    }
    if (ic == 2) drawLegend();
  }

  const std::string out = kOutDir + "/" + outName;
  c.SaveAs(out.c_str());
  std::cout << "[DONE] wrote " << out << "\n";
}

double meanValue(const std::vector<Point>& pts)
{
  double sum = 0.0;
  int n = 0;
  for (const auto& p : pts)
  {
    if (std::isfinite(p.y))
    {
      sum += p.y;
      ++n;
    }
  }
  return n > 0 ? sum / n : 0.0;
}

void drawSignalVsBackground()
{
  TCanvas c("c_signal_vs_background_target80", "c_signal_vs_background_target80", 2100, 700);
  c.Divide(3, 1, 0.010, 0.0);

  std::ofstream csv(kOutDir + "/signal_vs_background_tight_fraction_points.csv");
  csv << "config_key,label,cent,pt_mid,signal_tight_eff,background_tight_fraction\n";

  std::vector<std::unique_ptr<TFile>> files;
  std::vector<std::unique_ptr<TGraphErrors>> graphs;
  std::vector<std::unique_ptr<TH1F>> frames;

  for (size_t ic = 0; ic < kCents.size(); ++ic)
  {
    c.cd(ic + 1);
    gPad->SetTicks(1, 1);
    gPad->SetGrid(0, 0);
    gPad->SetLeftMargin(ic == 0 ? 0.155 : 0.060);
    gPad->SetRightMargin(ic == kCents.size() - 1 ? 0.035 : 0.020);
    gPad->SetTopMargin(0.245);
    gPad->SetBottomMargin(0.14);

    auto frame = std::make_unique<TH1F>(("frame_sigbkg_" + kCents[ic].suffix).c_str(), "", 100, 0.0, 1.02);
    frame->SetDirectory(nullptr);
    frame->SetStats(false);
    frame->SetMinimum(0.0);
    frame->SetMaximum(1.02);
    frame->GetXaxis()->SetTitle("Inclusive-jet tight fraction");
    frame->GetYaxis()->SetTitle(ic == 0 ? "Signal tight-ID efficiency" : "");
    frame->GetXaxis()->SetTitleSize(0.043);
    frame->GetYaxis()->SetTitleSize(0.043);
    frame->GetXaxis()->SetLabelSize(0.035);
    frame->GetYaxis()->SetLabelSize(0.035);
    frame->GetYaxis()->SetTitleOffset(ic == 0 ? 1.28 : 1.02);
    frame->Draw("AXIS");
    frames.push_back(std::move(frame));

    drawPanelLabel(kCents[ic].title);
    if (ic == 0) drawHeader("15 #leq E_{T}^{cluster} < 35 GeV; R = 0.3, sliding isolation");

    for (const auto& curve : kCurves)
    {
      auto fs = std::unique_ptr<TFile>(TFile::Open(sigPath(curve).c_str(), "READ"));
      auto fb = std::unique_ptr<TFile>(TFile::Open(bkgPath(curve).c_str(), "READ"));
      if (!fs || fs->IsZombie() || !fb || fb->IsZombie()) continue;
      TDirectory* ds = dynamic_cast<TDirectory*>(fs->Get("SIM"));
      TDirectory* db = dynamic_cast<TDirectory*>(fb->Get("SIM"));
      auto sig = readTightFractionPoints(ds, kCents[ic]);
      auto bkg = readTightFractionPoints(db, kCents[ic]);
      auto g = std::make_unique<TGraphErrors>();
      g->SetMarkerStyle(curve.marker);
      g->SetMarkerSize(1.20);
      g->SetMarkerColor(curve.color);
      g->SetLineColor(curve.color);
      g->SetLineWidth(2);
      const size_t n = std::min(sig.size(), bkg.size());
      for (size_t i = 0; i < n; ++i)
      {
        const int ip = g->GetN();
        g->SetPoint(ip, bkg[i].y, sig[i].y);
        g->SetPointError(ip, bkg[i].ey, sig[i].ey);
        csv << curve.cfg << "," << curve.label << "," << kCents[ic].suffix
            << "," << sig[i].x << "," << sig[i].y << "," << bkg[i].y << "\n";
      }
      g->Draw("P SAME");
      graphs.push_back(std::move(g));
      files.push_back(std::move(fs));
      files.push_back(std::move(fb));
    }
    if (ic == 2) drawLegend();
  }
  const std::string out = kOutDir + "/signal_efficiency_vs_background_tight_fraction_target80.png";
  c.SaveAs(out.c_str());
  std::cout << "[DONE] wrote " << out << "\n";
}

void drawReferenceVsBdtRecoAndId()
{
  TCanvas c("c_reference_vs_bdt_reco_id", "c_reference_vs_bdt_reco_id", 2600, 1450);
  c.Divide(3, 2, 0.006, 0.0);

  std::ofstream csv(kOutDir + "/reference_vs_bdt_reco_and_id_points.csv");
  csv << "quantity,config_key,label,cent,pt_mid,numerator_or_truth,denominator_piece_or_miss,value,error\n";

  std::vector<std::unique_ptr<TFile>> files;
  std::vector<std::unique_ptr<TGraphErrors>> graphs;
  std::vector<std::unique_ptr<TH1F>> frames;

  for (int row = 0; row < 2; ++row)
  {
    const bool useReco = row == 0;
    const std::string quantity = useReco ? "reco_efficiency_box_vs_bdt" : "signal_id_efficiency_box_vs_bdt";
    const std::string yTitle = useReco ? "Truth-photon recovery" : "Candidate tight-ID fraction";
    const double yMax = useReco ? 0.68 : 1.04;

    for (size_t ic = 0; ic < kCents.size(); ++ic)
    {
      c.cd(row * 3 + ic + 1);
      gPad->SetTicks(1, 1);
      gPad->SetGrid(0, 0);
      gPad->SetLeftMargin(ic == 0 ? 0.155 : 0.055);
      gPad->SetRightMargin(ic == kCents.size() - 1 ? 0.030 : 0.018);
      gPad->SetTopMargin(row == 0 ? 0.245 : 0.105);
      gPad->SetBottomMargin(row == 1 ? 0.160 : 0.045);

      auto frame = std::make_unique<TH1F>(
          ("frame_ref_bdt_" + std::to_string(row) + "_" + kCents[ic].suffix).c_str(), "", 100, 14.4, 35.6);
      frame->SetDirectory(nullptr);
      frame->SetStats(false);
      frame->SetMinimum(0.0);
      frame->SetMaximum(yMax);
      frame->GetXaxis()->SetTitle(row == 1 ? "Photon E_{T} [GeV]" : "");
      frame->GetYaxis()->SetTitle(ic == 0 ? yTitle.c_str() : "");
      frame->GetXaxis()->SetTitleSize(0.050);
      frame->GetYaxis()->SetTitleSize(0.048);
      frame->GetXaxis()->SetLabelSize(row == 1 ? 0.041 : 0.0);
      frame->GetYaxis()->SetLabelSize(0.039);
      frame->GetYaxis()->SetTitleOffset(ic == 0 ? 1.22 : 1.02);
      frame->Draw("AXIS");
      frames.push_back(std::move(frame));

      if (row == 0) drawPanelLabel(kCents[ic].title);
      if (row == 0 && ic == 0) drawHeader("15 #leq E_{T}^{cluster} < 35 GeV; R = 0.3, sliding isolation");

      for (const auto& curve : kReferenceVsBdtCurves)
      {
        auto f = std::unique_ptr<TFile>(TFile::Open(sigPath(curve).c_str(), "READ"));
        if (!f || f->IsZombie())
        {
          std::cerr << "[WARN] missing file: " << sigPath(curve) << "\n";
          continue;
        }
        TDirectory* d = dynamic_cast<TDirectory*>(f->Get("SIM"));
        if (!d)
        {
          std::cerr << "[WARN] missing SIM directory: " << sigPath(curve) << "\n";
          continue;
        }
        auto pts = useReco ? readRecoPoints(d, kCents[ic]) : readTightFractionPoints(d, kCents[ic]);
        if (!pts.empty())
        {
          auto g = makeGraph(curve, pts, csv, quantity, kCents[ic].suffix);
          g->Draw("P SAME");
          graphs.push_back(std::move(g));
        }
        files.push_back(std::move(f));
      }

      if (row == 0 && ic == 2)
      {
        auto leg = new TLegend(0.49, 0.815, 0.93, 0.925);
        leg->SetBorderSize(0);
        leg->SetFillStyle(1001);
        leg->SetFillColorAlpha(kWhite, 0.94);
        leg->SetTextFont(42);
        leg->SetTextSize(0.036);
        static std::vector<std::unique_ptr<TGraphErrors>> keepAlive;
        for (const auto& curve : kReferenceVsBdtCurves)
        {
          auto g = std::make_unique<TGraphErrors>(1);
          g->SetMarkerStyle(curve.marker);
          g->SetMarkerSize(1.30);
          g->SetMarkerColor(curve.color);
          leg->AddEntry(g.get(), curve.label.c_str(), "p");
          keepAlive.push_back(std::move(g));
        }
        leg->Draw();
      }
    }
  }

  const std::string out = kOutDir + "/target80_reference_vs_bdt_reco_and_id_efficiency_2x3.png";
  c.SaveAs(out.c_str());
  std::cout << "[DONE] wrote " << out << "\n";
}

void drawReferenceVsBdtIdOnly()
{
  TCanvas c("c_reference_vs_bdt_id_only", "c_reference_vs_bdt_id_only", 3000, 1020);
  c.Divide(3, 1, 0.006, 0.0);

  std::ofstream csv(kOutDir + "/reference_vs_bdt_signal_id_efficiency_points.csv");
  csv << "quantity,config_key,label,cent,pt_mid,tight_count,non_tight_count,value,error\n";

  std::vector<std::unique_ptr<TFile>> files;
  std::vector<std::unique_ptr<TGraphErrors>> graphs;
  std::vector<std::unique_ptr<TH1F>> frames;

  for (size_t ic = 0; ic < kCents.size(); ++ic)
  {
    c.cd(ic + 1);
    gPad->SetTicks(1, 1);
    gPad->SetGrid(0, 0);
    gPad->SetLeftMargin(ic == 0 ? 0.150 : 0.055);
    gPad->SetRightMargin(ic == kCents.size() - 1 ? 0.030 : 0.018);
    gPad->SetTopMargin(0.280);
    gPad->SetBottomMargin(0.150);

    auto frame = std::make_unique<TH1F>(("frame_id_only_" + kCents[ic].suffix).c_str(), "", 100, 14.4, 35.6);
    frame->SetDirectory(nullptr);
    frame->SetStats(false);
    frame->SetMinimum(0.0);
    frame->SetMaximum(1.05);
    frame->GetXaxis()->SetTitle("Photon candidate E_{T} [GeV]");
    frame->GetYaxis()->SetTitle(ic == 0 ? "Tight-ID efficiency" : "");
    frame->GetXaxis()->SetTitleSize(0.048);
    frame->GetYaxis()->SetTitleSize(0.046);
    frame->GetXaxis()->SetLabelSize(0.039);
    frame->GetYaxis()->SetLabelSize(0.039);
    frame->GetYaxis()->SetTitleOffset(ic == 0 ? 1.25 : 1.00);
    frame->Draw("AXIS");
    frames.push_back(std::move(frame));

    drawPanelLabel(kCents[ic].title);
    if (ic == 0) drawHeader("newPPG12 preselection; 15 #leq E_{T}^{cluster} < 35 GeV");

    for (const auto& curve : kReferenceVsBdtCurves)
    {
      auto f = std::unique_ptr<TFile>(TFile::Open(sigPath(curve).c_str(), "READ"));
      if (!f || f->IsZombie())
      {
        std::cerr << "[WARN] missing file: " << sigPath(curve) << "\n";
        continue;
      }
      TDirectory* d = dynamic_cast<TDirectory*>(f->Get("SIM"));
      if (!d)
      {
        std::cerr << "[WARN] missing SIM directory: " << sigPath(curve) << "\n";
        continue;
      }
      auto pts = readTightFractionPoints(d, kCents[ic]);
      if (!pts.empty())
      {
        auto g = makeGraph(curve, pts, csv, "signal_tight_id_efficiency", kCents[ic].suffix);
        g->Draw("P SAME");
        graphs.push_back(std::move(g));
      }
      files.push_back(std::move(f));
    }

    if (ic == 2)
    {
      auto leg = new TLegend(0.48, 0.805, 0.93, 0.925);
      leg->SetBorderSize(0);
      leg->SetFillStyle(1001);
      leg->SetFillColorAlpha(kWhite, 0.94);
      leg->SetTextFont(42);
      leg->SetTextSize(0.035);
      static std::vector<std::unique_ptr<TGraphErrors>> keepAlive;
      for (const auto& curve : kReferenceVsBdtCurves)
      {
        auto g = std::make_unique<TGraphErrors>(1);
        g->SetMarkerStyle(curve.marker);
        g->SetMarkerSize(1.30);
        g->SetMarkerColor(curve.color);
        leg->AddEntry(g.get(), curve.label.c_str(), "p");
        keepAlive.push_back(std::move(g));
      }
      leg->Draw();
    }
  }

  const std::string out = kOutDir + "/target80_reference_vs_best_bdt_tight_id_efficiency_1x3.png";
  c.SaveAs(out.c_str());
  std::cout << "[DONE] wrote " << out << "\n";
}

void drawSingleCurveThreePanel(const CurveSpec& curve,
                               const std::string& quantity,
                               const std::string& yTitle,
                               double yMin,
                               double yMax,
                               const std::string& xTitle,
                               const std::string& headerExtra,
                               const std::string& outName,
                               bool useBackground,
                               bool useReco)
{
  TCanvas c(("c_single_" + quantity + "_" + slugLabel(curve.label)).c_str(), "", 2800, 980);
  c.Divide(3, 1, 0.006, 0.0);

  std::ofstream csv(kOutDir + "/" + quantity + "_" + slugLabel(curve.label) + "_points.csv");
  csv << "quantity,config_key,label,cent,pt_mid,numerator_or_truth,denominator_piece_or_miss,value,error\n";

  std::vector<std::unique_ptr<TFile>> files;
  std::vector<std::unique_ptr<TGraphErrors>> graphs;
  std::vector<std::unique_ptr<TH1F>> frames;

  for (size_t ic = 0; ic < kCents.size(); ++ic)
  {
    c.cd(ic + 1);
    gPad->SetTicks(1, 1);
    gPad->SetGrid(0, 0);
    gPad->SetLeftMargin(ic == 0 ? 0.150 : 0.055);
    gPad->SetRightMargin(ic == kCents.size() - 1 ? 0.030 : 0.018);
    gPad->SetTopMargin(0.235);
    gPad->SetBottomMargin(0.150);

    auto frame = std::make_unique<TH1F>(("frame_single_" + quantity + "_" + kCents[ic].suffix).c_str(), "", 100, 14.4, 35.6);
    frame->SetDirectory(nullptr);
    frame->SetStats(false);
    frame->SetMinimum(yMin);
    frame->SetMaximum(yMax);
    frame->GetXaxis()->SetTitle(xTitle.c_str());
    frame->GetYaxis()->SetTitle(ic == 0 ? yTitle.c_str() : "");
    frame->GetXaxis()->SetTitleSize(0.048);
    frame->GetYaxis()->SetTitleSize(0.046);
    frame->GetXaxis()->SetLabelSize(0.039);
    frame->GetYaxis()->SetLabelSize(0.039);
    frame->GetYaxis()->SetTitleOffset(ic == 0 ? 1.25 : 1.00);
    frame->Draw("AXIS");
    frames.push_back(std::move(frame));

    drawPanelLabel(kCents[ic].title);
    if (ic == 0) drawHeader(headerExtra);

    const std::string path = useBackground ? bkgPath(curve) : sigPath(curve);
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
    auto pts = useReco ? readRecoPoints(d, kCents[ic]) : readTightFractionPoints(d, kCents[ic]);
    if (!pts.empty())
    {
      CurveSpec drawCurve = curve;
      drawCurve.xShift = 0.0;
      auto g = makeGraph(drawCurve, pts, csv, quantity, kCents[ic].suffix);
      g->Draw("P SAME");
      graphs.push_back(std::move(g));
    }
    files.push_back(std::move(f));

    if (ic == 2)
    {
      auto leg = new TLegend(0.48, 0.825, 0.93, 0.925);
      leg->SetBorderSize(0);
      leg->SetFillStyle(1001);
      leg->SetFillColorAlpha(kWhite, 0.94);
      leg->SetTextFont(42);
      leg->SetTextSize(0.035);
      static std::vector<std::unique_ptr<TGraphErrors>> keepAlive;
      auto g = std::make_unique<TGraphErrors>(1);
      g->SetMarkerStyle(curve.marker);
      g->SetMarkerSize(1.30);
      g->SetMarkerColor(curve.color);
      leg->AddEntry(g.get(), curve.label.c_str(), "p");
      keepAlive.push_back(std::move(g));
      leg->Draw();
    }
  }

  const std::string out = kOutDir + "/" + outName;
  c.SaveAs(out.c_str());
  std::cout << "[DONE] wrote " << out << "\n";
}

void drawEveryVariationSingles()
{
  for (const auto& curve : kCurves)
  {
    const std::string slug = slugLabel(curve.label);
    drawSingleCurveThreePanel(curve,
                              "truth_photon_recovery",
                              "Truth-photon recovery",
                              0.0,
                              0.70,
                              "Truth photon E_{T} [GeV]",
                              "15 #leq E_{T}^{cluster} < 35 GeV; R = 0.3, sliding isolation",
                              "single_" + slug + "_truth_photon_recovery_1x3.png",
                              false,
                              true);

    drawSingleCurveThreePanel(curve,
                              "candidate_tight_fraction",
                              "Candidate tight-ID fraction",
                              0.0,
                              1.05,
                              "Photon candidate E_{T} [GeV]",
                              "15 #leq E_{T}^{cluster} < 35 GeV; R = 0.3, sliding isolation",
                              "single_" + slug + "_candidate_tight_fraction_1x3.png",
                              false,
                              false);

    drawSingleCurveThreePanel(curve,
                              "inclusive_jet_tight_fraction",
                              "Inclusive-jet tight fraction",
                              0.0,
                              1.05,
                              "Photon candidate E_{T} [GeV]",
                              "15 #leq E_{T}^{cluster} < 35 GeV; R = 0.3, sliding isolation",
                              "single_" + slug + "_inclusive_jet_tight_fraction_1x3.png",
                              true,
                              false);
  }
}
}  // namespace

void PlotTarget80FirstPairEfficiency()
{
  SetsPhenixStyle();
  gStyle->SetOptStat(0);
  gStyle->SetPadGridX(false);
  gStyle->SetPadGridY(false);
  gSystem->mkdir(kOutDir.c_str(), true);

  drawThreePanel("reco_efficiency",
                 "Reconstruction efficiency",
                 0.0,
                 0.70,
                 "Truth photon E_{T} [GeV]",
                 "15 #leq E_{T}^{cluster} < 35 GeV; R = 0.3, sliding isolation",
                 "target80_reco_efficiency_1x3.png",
                 false,
                 true);

  drawThreePanel("signal_tight_id_efficiency",
                 "Signal tight-ID efficiency",
                 0.0,
                 1.04,
                 "Photon candidate E_{T} [GeV]",
                 "15 #leq E_{T}^{cluster} < 35 GeV; R = 0.3, sliding isolation",
                 "target80_signal_id_efficiency_1x3.png",
                 false,
                 false);

  drawThreePanel("inclusive_jet_tight_fraction",
                 "Inclusive-jet tight fraction",
                 0.0,
                 1.04,
                 "Photon candidate E_{T} [GeV]",
                 "15 #leq E_{T}^{cluster} < 35 GeV; R = 0.3, sliding isolation",
                 "target80_inclusive_jet_tight_fraction_1x3.png",
                 true,
                 false);

  drawSignalVsBackground();

  drawReferenceVsBdtRecoAndId();

  drawReferenceVsBdtIdOnly();

  drawEveryVariationSingles();
}
