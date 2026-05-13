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
const std::string kTargetBase =
    "dataOutput/target80_first_offline/bdt_target80_gated_20260512_001012/"
    "analysis_config_etfine_15to35_target80";
const std::string kOutDir =
    "dataOutput/target80_first_offline/bdt_target80_gated_20260512_001012/"
    "analysis_config_etfine_15to35_target80/efficiency_qa";
const std::string kSignalRel =
    "photonJet12and20merged_SIM/RecoilJets_embeddedPhoton12plus20_MERGED.root";
const std::string kBkgRel =
    "embeddedJet12and20merged_SIM/RecoilJets_embeddedJet12plus20_MERGED.root";
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
  double lo = 0.0;
  double hi = 0.0;
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
  std::string label;
  std::string cfg;
  int color = kBlack;
  int marker = 20;
  double shift = 0.0;
};

const std::vector<CentSpec> kCents = {
    {"0_20", "0-20% central"},
    {"20_50", "20-50% mid-central"},
    {"50_80", "50-80% peripheral"}};

const std::vector<PtBin> kPtBins = {
    {15, 17}, {17, 19}, {19, 21}, {21, 23}, {23, 25}, {25, 27}, {27, 30}, {30, 35}};

const CurveSpec kAuAuBox = {
    "Au+Au box-cuts",
    "preselectionNewPPG12_tightReference_nonTightReference_baseVariant",
    kBlack,
    20,
    -0.13};

const CurveSpec kAuAuTarget80 = {
    "Au+Au BDT target-80",
    "preselectionNewPPG12_tightAuAuCentInputBase3x3BDT_nonTightAuAuBDTComplement_baseVariant",
    kAzure + 2,
    20,
    0.13};

const CurveSpec kPPReference = {
    "pp reference ID",
    "",
    kGreen + 2,
    21,
    0.0};

std::string signalPath(const CurveSpec& c)
{
  return kTargetBase + "/simembedded/" + c.cfg + "/" + kSignalRel;
}

std::string bkgPath(const CurveSpec& c)
{
  return kTargetBase + "/simembeddedinclusive/" + c.cfg + "/" + kBkgRel;
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

std::vector<Point> readAuAuRecovery(TDirectory* d, const CentSpec& cent)
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

    Point p;
    p.x = mid;
    p.ex = 0.0;
    p.y = 1.0 - miss / truth;
    p.num = truth - miss;
    p.den = truth;
    const double eTruth = hTruth->GetBinError(ib);
    const double eMiss = hMiss->GetBinError(im);
    const double missFrac = miss / truth;
    p.ey = miss > 0.0 ? missFrac * std::hypot(eMiss / miss, eTruth / truth) : 0.0;
    if (!std::isfinite(p.ey)) p.ey = 0.0;
    pts.push_back(p);
  }
  return pts;
}

std::vector<Point> readPPRecovery(TDirectory* d)
{
  std::vector<Point> pts;
  TH1* hTruth = getHist(d, "h_unfoldTruthPho_pTgamma_ppg12obj");
  TH1* hMiss = getHist(d, "h_unfoldTruthPhoMisses_pTgamma_ppg12obj");
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

    Point p;
    p.x = mid;
    p.ex = 0.0;
    p.y = 1.0 - miss / truth;
    p.num = truth - miss;
    p.den = truth;
    const double eTruth = hTruth->GetBinError(ib);
    const double eMiss = hMiss->GetBinError(im);
    const double missFrac = miss / truth;
    p.ey = miss > 0.0 ? missFrac * std::hypot(eMiss / miss, eTruth / truth) : 0.0;
    if (!std::isfinite(p.ey)) p.ey = 0.0;
    pts.push_back(p);
  }
  return pts;
}

Point readCandidateTightFraction(TDirectory* d, const PtBin& pt, const CentSpec& cent)
{
  const std::string tag = "isoR30_pT_" + std::to_string((int)pt.lo) + "_" +
                          std::to_string((int)pt.hi) + "_cent_" + cent.suffix;
  double eTight = 0.0;
  double eNonTight = 0.0;
  const double tight = integralAndError(getHist(d, "h_Eiso_tight_" + tag), eTight);
  const double nonTight = integralAndError(getHist(d, "h_Eiso_nonTight_" + tag), eNonTight);
  return ratioPoint(tight, nonTight, eTight, eNonTight, 0.5 * (pt.lo + pt.hi));
}

std::vector<Point> readCandidateTightFractions(TDirectory* d, const CentSpec& cent)
{
  std::vector<Point> pts;
  for (const auto& pt : kPtBins)
  {
    Point p = readCandidateTightFraction(d, pt, cent);
    if (p.den > 0.0) pts.push_back(p);
  }
  return pts;
}

std::unique_ptr<TGraphErrors> graphFromPoints(const CurveSpec& curve,
                                              const std::vector<Point>& points,
                                              std::ofstream* csv,
                                              const std::string& quantity,
                                              const std::string& cent)
{
  auto g = std::make_unique<TGraphErrors>();
  g->SetMarkerStyle(curve.marker);
  g->SetMarkerColor(curve.color);
  g->SetLineColor(curve.color);
  g->SetMarkerSize(1.25);
  g->SetLineWidth(2);
  for (const auto& p : points)
  {
    const int i = g->GetN();
    g->SetPoint(i, p.x + curve.shift, p.y);
    g->SetPointError(i, p.ex, p.ey);
    if (csv)
    {
      (*csv) << quantity << "," << curve.label << "," << cent << ","
             << std::setprecision(9) << p.x << "," << p.num << "," << p.den << ","
             << p.y << "," << p.ey << "\n";
    }
  }
  return g;
}

void drawExperimentLabel(double x = 0.15, double y = 0.875, double size = 0.040)
{
  TLatex tx;
  tx.SetNDC();
  tx.SetTextFont(42);
  tx.SetTextAlign(13);
  tx.SetTextSize(size);
  tx.DrawLatex(x, y, "#it{#bf{sPHENIX}} Internal");
}

void drawCollisionLabel(const char* extra, double x = 0.15, double y = 0.800, double size = 0.030)
{
  TLatex tx;
  tx.SetNDC();
  tx.SetTextFont(42);
  tx.SetTextAlign(13);
  tx.SetTextSize(size);
  tx.DrawLatex(x, y, "Pythia overlay, #sqrt{s_{NN}} = 200 GeV");
  tx.DrawLatex(x, y - 0.055, extra);
}

void drawPanelTitle(const std::string& title)
{
  TLatex tx;
  tx.SetNDC();
  tx.SetTextFont(42);
  tx.SetTextSize(0.052);
  tx.SetTextAlign(22);
  tx.DrawLatex(0.50, 0.955, title.c_str());
}

double meanY(const std::vector<Point>& pts)
{
  double s = 0.0;
  int n = 0;
  for (const auto& p : pts)
  {
    if (std::isfinite(p.y))
    {
      s += p.y;
      ++n;
    }
  }
  return n ? s / n : 0.0;
}

void drawPPOverlay()
{
  auto fBox = std::unique_ptr<TFile>(TFile::Open(signalPath(kAuAuBox).c_str(), "READ"));
  auto fBdt = std::unique_ptr<TFile>(TFile::Open(signalPath(kAuAuTarget80).c_str(), "READ"));
  auto fPP = std::unique_ptr<TFile>(TFile::Open(kPPFile.c_str(), "READ"));
  if (!fBox || fBox->IsZombie() || !fBdt || fBdt->IsZombie() || !fPP || fPP->IsZombie())
  {
    std::cerr << "[ERROR] Failed to open one or more input files.\n";
    return;
  }

  TDirectory* dBox = simDir(fBox.get());
  TDirectory* dBdt = simDir(fBdt.get());
  TDirectory* dPP = simDir(fPP.get());
  if (!dBox || !dBdt || !dPP)
  {
    std::cerr << "[ERROR] Missing SIM directory.\n";
    return;
  }

  gSystem->mkdir(kOutDir.c_str(), true);
  std::ofstream csv(kOutDir + "/target80_vs_pp_reco_efficiency_points.csv");
  csv << "quantity,label,cent,pt_mid,n_pass_like,n_truth,value,error\n";

  const std::vector<Point> ppPts = readPPRecovery(dPP);

  TCanvas c("c_target80_vs_pp_reco", "c_target80_vs_pp_reco", 3370, 1120);
  c.Divide(3, 1, 0.006, 0.0);

  std::vector<std::unique_ptr<TGraphErrors>> graphs;
  std::vector<std::unique_ptr<TH1F>> frames;
  for (size_t ic = 0; ic < kCents.size(); ++ic)
  {
    c.cd(ic + 1);
    gPad->SetTicks(1, 1);
    gPad->SetGrid(0, 0);
    gPad->SetLeftMargin(ic == 0 ? 0.150 : 0.058);
    gPad->SetRightMargin(ic == kCents.size() - 1 ? 0.030 : 0.018);
    gPad->SetTopMargin(0.300);
    gPad->SetBottomMargin(0.145);

    auto frame = std::make_unique<TH1F>(("frame_pp_overlay_" + kCents[ic].suffix).c_str(), "", 100, 14.3, 35.7);
    frame->SetDirectory(nullptr);
    frame->SetStats(false);
    frame->SetMinimum(0.0);
    frame->SetMaximum(1.02);
    frame->GetXaxis()->SetTitle("Truth photon E_{T} [GeV]");
    frame->GetYaxis()->SetTitle(ic == 0 ? "Truth-photon recovery" : "");
    frame->GetXaxis()->SetTitleSize(0.047);
    frame->GetYaxis()->SetTitleSize(0.045);
    frame->GetXaxis()->SetLabelSize(0.039);
    frame->GetYaxis()->SetLabelSize(0.039);
    frame->GetYaxis()->SetTitleOffset(ic == 0 ? 1.25 : 1.0);
    frame->Draw("AXIS");
    frames.push_back(std::move(frame));

    drawPanelTitle(kCents[ic].title);
    if (ic == 0)
    {
      drawExperimentLabel(0.145, 0.870);
      drawCollisionLabel("Au+Au: R = 0.3 sliding iso; pp: reference ID, R = 0.4 fixed iso", 0.145, 0.810);
    }

    auto boxPts = readAuAuRecovery(dBox, kCents[ic]);
    auto bdtPts = readAuAuRecovery(dBdt, kCents[ic]);

    auto gBox = graphFromPoints(kAuAuBox, boxPts, &csv, "truth_photon_recovery", kCents[ic].suffix);
    auto gBdt = graphFromPoints(kAuAuTarget80, bdtPts, &csv, "truth_photon_recovery", kCents[ic].suffix);
    auto gPP = graphFromPoints(kPPReference, ppPts, &csv, "truth_photon_recovery", "pp_reference");
    gPP->SetMarkerSize(1.35);
    gPP->SetLineStyle(1);

    gBox->Draw("P SAME");
    gBdt->Draw("P SAME");
    gPP->Draw("P SAME");

    graphs.push_back(std::move(gBox));
    graphs.push_back(std::move(gBdt));
    graphs.push_back(std::move(gPP));

    if (ic == 2)
    {
      TLegend* leg = new TLegend(0.42, 0.740, 0.94, 0.925);
      leg->SetBorderSize(0);
      leg->SetFillStyle(1001);
      leg->SetFillColorAlpha(kWhite, 0.94);
      leg->SetTextFont(42);
      leg->SetTextSize(0.033);
      leg->SetNColumns(1);
      leg->AddEntry(graphs[graphs.size() - 3].get(), kAuAuBox.label.c_str(), "p");
      leg->AddEntry(graphs[graphs.size() - 2].get(), kAuAuTarget80.label.c_str(), "p");
      leg->AddEntry(graphs[graphs.size() - 1].get(), kPPReference.label.c_str(), "p");
      leg->Draw();
    }
  }

  const std::string out = kOutDir + "/target80_vs_pp_truth_photon_recovery_1x3.png";
  c.SaveAs(out.c_str());
  std::cout << "[DONE] wrote " << out << "\n";
}

void drawLossDecomposition()
{
  auto fBdt = std::unique_ptr<TFile>(TFile::Open(signalPath(kAuAuTarget80).c_str(), "READ"));
  auto fBkg = std::unique_ptr<TFile>(TFile::Open(bkgPath(kAuAuTarget80).c_str(), "READ"));
  if (!fBdt || fBdt->IsZombie() || !fBkg || fBkg->IsZombie())
  {
    std::cerr << "[ERROR] Failed to open target80 signal/background files.\n";
    return;
  }
  TDirectory* dSig = simDir(fBdt.get());
  TDirectory* dBkg = simDir(fBkg.get());
  if (!dSig || !dBkg) return;

  std::ofstream csv(kOutDir + "/target80_efficiency_loss_decomposition.csv");
  csv << "cent,pt_mid,truth_photon_recovery,candidate_tight_fraction,inferred_pre_tight_survival,inclusive_jet_tight_fraction\n";

  TCanvas c("c_target80_loss_decomposition", "c_target80_loss_decomposition", 3370, 1120);
  c.Divide(3, 1, 0.006, 0.0);

  std::vector<std::unique_ptr<TGraphErrors>> keep;
  std::vector<std::unique_ptr<TH1F>> frames;
  for (size_t ic = 0; ic < kCents.size(); ++ic)
  {
    c.cd(ic + 1);
    gPad->SetTicks(1, 1);
    gPad->SetGrid(0, 0);
    gPad->SetLeftMargin(ic == 0 ? 0.150 : 0.058);
    gPad->SetRightMargin(ic == kCents.size() - 1 ? 0.030 : 0.018);
    gPad->SetTopMargin(0.300);
    gPad->SetBottomMargin(0.145);

    auto frame = std::make_unique<TH1F>(("frame_loss_" + kCents[ic].suffix).c_str(), "", 100, 14.3, 35.7);
    frame->SetDirectory(nullptr);
    frame->SetStats(false);
    frame->SetMinimum(0.0);
    frame->SetMaximum(1.02);
    frame->GetXaxis()->SetTitle("Photon E_{T} [GeV]");
    frame->GetYaxis()->SetTitle(ic == 0 ? "Efficiency or tight fraction" : "");
    frame->GetXaxis()->SetTitleSize(0.047);
    frame->GetYaxis()->SetTitleSize(0.045);
    frame->GetXaxis()->SetLabelSize(0.039);
    frame->GetYaxis()->SetLabelSize(0.039);
    frame->GetYaxis()->SetTitleOffset(ic == 0 ? 1.25 : 1.0);
    frame->Draw("AXIS");
    frames.push_back(std::move(frame));

    drawPanelTitle(kCents[ic].title);
    if (ic == 0)
    {
      drawExperimentLabel(0.145, 0.870);
      drawCollisionLabel("Au+Au BDT target-80; R = 0.3 sliding isolation", 0.145, 0.810);
    }

    std::vector<Point> recovery = readAuAuRecovery(dSig, kCents[ic]);
    std::vector<Point> id = readCandidateTightFractions(dSig, kCents[ic]);
    std::vector<Point> bkg = readCandidateTightFractions(dBkg, kCents[ic]);
    std::vector<Point> preLike;
    const size_t n = std::min(recovery.size(), id.size());
    preLike.reserve(n);
    for (size_t i = 0; i < n; ++i)
    {
      Point p;
      p.x = recovery[i].x;
      p.y = id[i].y > 0.0 ? recovery[i].y / id[i].y : 0.0;
      p.ey = p.y * std::hypot(recovery[i].y > 0.0 ? recovery[i].ey / recovery[i].y : 0.0,
                              id[i].y > 0.0 ? id[i].ey / id[i].y : 0.0);
      if (!std::isfinite(p.ey)) p.ey = 0.0;
      preLike.push_back(p);
      if (i < bkg.size())
      {
        csv << kCents[ic].suffix << "," << p.x << ","
            << recovery[i].y << "," << id[i].y << "," << p.y << "," << bkg[i].y << "\n";
      }
    }

    CurveSpec rec{"Truth-photon recovery", "", kAzure + 2, 20, -0.18};
    CurveSpec idc{"Candidate tight fraction", "", kViolet + 1, 21, -0.06};
    CurveSpec pre{"Recovery / tight fraction", "", kOrange + 7, 22, 0.06};
    CurveSpec bgc{"Inclusive-jet tight fraction", "", kGray + 2, 33, 0.18};
    auto gRec = graphFromPoints(rec, recovery, nullptr, "", "");
    auto gId = graphFromPoints(idc, id, nullptr, "", "");
    auto gPre = graphFromPoints(pre, preLike, nullptr, "", "");
    auto gBkg = graphFromPoints(bgc, bkg, nullptr, "", "");

    gRec->Draw("P SAME");
    gId->Draw("P SAME");
    gPre->Draw("P SAME");
    gBkg->Draw("P SAME");

    if (ic == 2)
    {
      TLegend* leg = new TLegend(0.37, 0.710, 0.94, 0.925);
      leg->SetBorderSize(0);
      leg->SetFillStyle(1001);
      leg->SetFillColorAlpha(kWhite, 0.94);
      leg->SetTextFont(42);
      leg->SetTextSize(0.030);
      leg->AddEntry(gRec.get(), "Truth-photon recovery", "p");
      leg->AddEntry(gId.get(), "Candidate tight fraction", "p");
      leg->AddEntry(gPre.get(), "Recovery / tight fraction", "p");
      leg->AddEntry(gBkg.get(), "Inclusive-jet tight fraction", "p");
      leg->Draw();
    }

    keep.push_back(std::move(gRec));
    keep.push_back(std::move(gId));
    keep.push_back(std::move(gPre));
    keep.push_back(std::move(gBkg));
  }

  const std::string out = kOutDir + "/target80_efficiency_loss_decomposition_1x3.png";
  c.SaveAs(out.c_str());
  std::cout << "[DONE] wrote " << out << "\n";
}
}  // namespace

void PlotTarget80VsPPAndLoss()
{
  SetsPhenixStyle();
  gStyle->SetOptStat(0);
  drawPPOverlay();
  drawLossDecomposition();
}
