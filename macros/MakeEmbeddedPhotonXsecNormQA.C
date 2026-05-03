#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TH1.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace
{
// IMPORTANT:
//   After RecoilJets_AuAu stitching, PhotonJet12 must use the EXCLUSIVE
//   generator cross section for 12 <= pT_filter^gamma < 20 GeV.
constexpr double kSigmaEmbeddedPhoton12To20_pb = 2598.12425;
constexpr double kSigmaEmbeddedPhoton20Plus_pb = 133.317866;

const std::string kConfigTag =
    "jetMinPt5_7pi_8_vz60_isoR40_fixedIso2GeV_baseVariant_"
    "preselectionReference_tightReference_nonTightReference";

const std::string kFile12 =
    "InputFiles/simEmbedded/RecoilJets_embeddedPhoton12_ALL_" + kConfigTag + ".root";
const std::string kFile20 =
    "InputFiles/simEmbedded/RecoilJets_embeddedPhoton20_ALL_" + kConfigTag + ".root";

const std::string kOutDir =
    "dataOutput/combinedSimOnlyEMBEDDED/" + kConfigTag + "/photonJet12and20merged_SIM";

int VzCutFromConfigTag()
{
  if (kConfigTag.find("vz60") != std::string::npos) return 60;
  if (kConfigTag.find("vz30") != std::string::npos) return 30;
  std::cerr << "[WARN] Could not infer vz cut from config tag: " << kConfigTag
            << ". Falling back to 30 cm." << std::endl;
  return 30;
}

double IsoConeRFromConfigTag()
{
  if (kConfigTag.find("isoR40") != std::string::npos) return 0.4;
  if (kConfigTag.find("isoR30") != std::string::npos) return 0.3;
  std::cerr << "[WARN] Could not infer isolation cone radius from config tag: " << kConfigTag
            << ". Falling back to 0.3." << std::endl;
  return 0.3;
}

std::string VzEtaLabel()
{
  return TString::Format("|v_{z}| < %d cm,  |#eta^{#gamma}| < 0.7", VzCutFromConfigTag()).Data();
}

std::string IsoConeLabel()
{
  return TString::Format("#DeltaR_{cone} < %.1f", IsoConeRFromConfigTag()).Data();
}

double ReadEventCount(TFile* f)
{
  if (!f) return 0.0;
  TDirectory* d = f->GetDirectory("SIM");
  if (!d) return 0.0;

  TH1* cnt = dynamic_cast<TH1*>(d->Get("cnt_SIM"));
  if (!cnt) return 0.0;

  const double bin1 = cnt->GetBinContent(1);
  if (bin1 > 0.0) return bin1;

  const double integral = cnt->Integral(0, cnt->GetNbinsX() + 1);
  if (integral > 0.0) return integral;

  return cnt->GetEntries();
}

TH1* CloneFromDir(TDirectory* d, const std::string& name, const std::string& cloneName)
{
  if (!d) return nullptr;
  TH1* h = dynamic_cast<TH1*>(d->Get(name.c_str()));
  if (!h) return nullptr;

  TH1* c = dynamic_cast<TH1*>(h->Clone(cloneName.c_str()));
  if (!c) return nullptr;
  c->SetDirectory(nullptr);
  if (c->GetSumw2N() == 0) c->Sumw2();
  return c;
}

std::unique_ptr<TH1> BuildRegionPtSpectrum(TFile* f,
                                           const std::vector<std::string>& regions,
                                           const std::string& cloneName)
{
  if (!f) return nullptr;
  TDirectory* d = f->GetDirectory("SIM");
  if (!d) return nullptr;

  auto addSet = [&](const std::vector<std::string>& suffixes) -> std::unique_ptr<TH1>
  {
    std::unique_ptr<TH1> sum;
    for (const std::string& suffix : suffixes)
    {
      for (const std::string& region : regions)
      {
        const std::string name = "h_pTgamma_ABCD_" + region + suffix;
        std::unique_ptr<TH1> h(CloneFromDir(d, name, cloneName + "_" + region + suffix));
        if (!h) continue;

        if (!sum)
        {
          sum.reset(dynamic_cast<TH1*>(h->Clone(cloneName.c_str())));
          if (!sum) return nullptr;
          sum->SetDirectory(nullptr);
          sum->Reset("ICES");
          if (sum->GetSumw2N() == 0) sum->Sumw2();
        }
        sum->Add(h.get());
      }
    }
    return sum;
  };

  std::unique_ptr<TH1> inclusive = addSet({""});
  if (inclusive && inclusive->Integral(0, inclusive->GetNbinsX() + 1) > 0.0)
  {
    return inclusive;
  }

  return addSet({"_cent_0_20", "_cent_20_50", "_cent_50_80"});
}

std::unique_ptr<TH1> BuildKeptFilterPtSpectrum(TFile* f, const std::string& cloneName)
{
  if (!f) return nullptr;
  TDirectory* d = f->GetDirectory("SIM");
  return std::unique_ptr<TH1>(CloneFromDir(d, "h_embedStitch_filterPhotonPt_kept", cloneName));
}

void StyleHist(TH1* h, int color, int marker)
{
  if (!h) return;
  h->SetLineColor(color);
  h->SetMarkerColor(color);
  h->SetMarkerStyle(marker);
  h->SetMarkerSize(0.85);
  h->SetLineWidth(2);
}

std::string Sci(double v, int p = 3)
{
  std::ostringstream os;
  os << std::scientific << std::setprecision(p) << v;
  return os.str();
}

std::unique_ptr<TH1> MakeSmoothReference(const TH1* h, const std::string& name)
{
  if (!h) return nullptr;
  std::unique_ptr<TH1> ref(dynamic_cast<TH1*>(h->Clone(name.c_str())));
  if (!ref) return nullptr;
  ref->SetDirectory(nullptr);
  ref->Smooth(2);

  for (int ib = 1; ib <= ref->GetNbinsX(); ++ib)
  {
    if (ref->GetBinContent(ib) <= 0.0 && h->GetBinContent(ib) > 0.0)
    {
      ref->SetBinContent(ib, h->GetBinContent(ib));
    }
  }
  return ref;
}

std::unique_ptr<TH1> MakeRatioToSmooth(const TH1* h, const TH1* smooth, const std::string& name)
{
  if (!h || !smooth) return nullptr;
  std::unique_ptr<TH1> ratio(dynamic_cast<TH1*>(h->Clone(name.c_str())));
  if (!ratio) return nullptr;
  ratio->SetDirectory(nullptr);
  ratio->Reset("ICES");
  if (ratio->GetSumw2N() == 0) ratio->Sumw2();

  for (int ib = 1; ib <= h->GetNbinsX(); ++ib)
  {
    const double den = smooth->GetBinContent(ib);
    const double num = h->GetBinContent(ib);
    if (den <= 0.0 || num <= 0.0) continue;
    ratio->SetBinContent(ib, num / den);
    ratio->SetBinError(ib, h->GetBinError(ib) / den);
  }
  return ratio;
}

std::unique_ptr<TH1> WeightedClone(std::unique_ptr<TH1> h, double weight)
{
  if (!h) return nullptr;
  h->Scale(weight);
  return h;
}

std::unique_ptr<TH1> BuildWeightedPairHistogram(TFile* f12,
                                                TFile* f20,
                                                const std::string& histName,
                                                const std::string& cloneName,
                                                double w12,
                                                double w20)
{
  TDirectory* d12 = f12 ? f12->GetDirectory("SIM") : nullptr;
  TDirectory* d20 = f20 ? f20->GetDirectory("SIM") : nullptr;

  std::unique_ptr<TH1> h12(CloneFromDir(d12, histName, cloneName + "_12"));
  std::unique_ptr<TH1> h20(CloneFromDir(d20, histName, cloneName + "_20"));
  if (!h12 && !h20) return nullptr;

  std::unique_ptr<TH1> hSum;
  if (h12)
  {
    hSum.reset(dynamic_cast<TH1*>(h12->Clone(cloneName.c_str())));
    if (!hSum) return nullptr;
    hSum->SetDirectory(nullptr);
    hSum->Scale(w12);
  }
  else if (h20)
  {
    hSum.reset(dynamic_cast<TH1*>(h20->Clone(cloneName.c_str())));
    if (!hSum) return nullptr;
    hSum->SetDirectory(nullptr);
    hSum->Reset("ICES");
    if (hSum->GetSumw2N() == 0) hSum->Sumw2();
  }

  if (h20)
  {
    h20->Scale(w20);
    hSum->Add(h20.get());
  }

  return hSum;
}

bool FindEfficiencyCut(TH1* hIn, double eff, double& cut, double& cutErr)
{
  cut = 0.0;
  cutErr = 0.0;
  if (!hIn) return false;

  const int nb = hIn->GetNbinsX();
  const double total = hIn->Integral(1, nb);
  if (!(total > 0.0)) return false;

  double running = 0.0;
  int ibCut = nb;
  for (int ib = 1; ib <= nb; ++ib)
  {
    running += hIn->GetBinContent(ib);
    if ((running / total) >= eff)
    {
      ibCut = ib;
      break;
    }
  }

  const double binLo = hIn->GetXaxis()->GetBinLowEdge(ibCut);
  const double binHi = hIn->GetXaxis()->GetBinUpEdge(ibCut);
  const double prev = running - hIn->GetBinContent(ibCut);
  const double binC = hIn->GetBinContent(ibCut);
  const double target = eff * total;

  if (binC > 0.0)
  {
    const double frac = std::min(1.0, std::max(0.0, (target - prev) / binC));
    cut = binLo + frac * (binHi - binLo);
  }
  else
  {
    cut = 0.5 * (binLo + binHi);
  }

  cutErr = 0.5 * (binHi - binLo);
  return std::isfinite(cut);
}

void StyleEffGraph(TGraphErrors& g, int color, int marker)
{
  g.SetLineWidth(2);
  g.SetLineColor(color);
  g.SetMarkerColor(color);
  g.SetMarkerStyle(marker);
  g.SetMarkerSize(marker == 22 ? 1.45 : 1.35);
}

struct IsoFlatCentResult
{
  int lo = 0;
  int hi = 0;
  double center = 0.0;
  double halfWidth = 0.0;
  bool have70 = false;
  bool have80 = false;
  bool have90 = false;
  double flat70 = 0.0;
  double err70 = 0.0;
  double flat80 = 0.0;
  double err80 = 0.0;
  double flat90 = 0.0;
  double err90 = 0.0;
};

bool ComputeFlatGraphAverage(TGraphErrors& g, double xLo, double xHi, double& value, double& error)
{
  value = 0.0;
  error = 0.0;

  double sumW = 0.0;
  double sumWY = 0.0;
  int nUsed = 0;

  for (int ip = 0; ip < g.GetN(); ++ip)
  {
    double x = 0.0;
    double y = 0.0;
    g.GetPoint(ip, x, y);
    if (x < xLo || x > xHi) continue;

    const double ey = g.GetErrorY(ip);
    const double w = (ey > 0.0) ? (1.0 / (ey * ey)) : 1.0;
    sumW += w;
    sumWY += w * y;
    ++nUsed;
  }

  if (nUsed <= 0 || !(sumW > 0.0)) return false;
  value = sumWY / sumW;
  error = std::sqrt(1.0 / sumW);
  return std::isfinite(value);
}

std::vector<IsoFlatCentResult> DrawIsoEfficiencyCutoffQA(TFile* f12, TFile* f20, double w12, double w20)
{
  const std::vector<double> ptEdges = {8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 35};
  const double fitXLo = 10.0;
  const double fitXHi = 35.0;
  struct CentBin
  {
    int lo;
    int hi;
    std::string suffix;
    std::string tag;
  };
  const std::vector<CentBin> centBins = {
      {0, 20, "_cent_0_20", "cent_0_20"},
      {20, 50, "_cent_20_50", "cent_20_50"},
      {50, 80, "_cent_50_80", "cent_50_80"},
  };

  std::vector<IsoFlatCentResult> flatResults;

  for (const auto& cb : centBins)
  {
    std::vector<double> x;
    std::vector<double> ex;
    std::vector<double> y70, ey70;
    std::vector<double> y80, ey80;
    std::vector<double> y90, ey90;
    double yMin = std::numeric_limits<double>::max();
    double yMax = -std::numeric_limits<double>::max();

    for (std::size_t ipt = 1; ipt + 1 < ptEdges.size(); ++ipt)
    {
      const int lo = static_cast<int>(std::llround(ptEdges[ipt]));
      const int hi = static_cast<int>(std::llround(ptEdges[ipt + 1]));
      const std::string histName =
          "h_EisoReco_truthSigMatched_pT_" + std::to_string(lo) + "_" + std::to_string(hi) + cb.suffix;

      std::unique_ptr<TH1> h = BuildWeightedPairHistogram(
          f12, f20, histName, "h_weightedIso_" + cb.tag + "_" + std::to_string(lo) + "_" + std::to_string(hi), w12, w20);
      if (!h || h->Integral(1, h->GetNbinsX()) <= 0.0) continue;

      double c70 = 0.0, e70 = 0.0;
      double c80 = 0.0, e80 = 0.0;
      double c90 = 0.0, e90 = 0.0;
      if (!FindEfficiencyCut(h.get(), 0.70, c70, e70)) continue;
      if (!FindEfficiencyCut(h.get(), 0.80, c80, e80)) continue;
      if (!FindEfficiencyCut(h.get(), 0.90, c90, e90)) continue;

      x.push_back(0.5 * (ptEdges[ipt] + ptEdges[ipt + 1]));
      ex.push_back(0.0);
      y70.push_back(c70);
      ey70.push_back(e70);
      y80.push_back(c80);
      ey80.push_back(e80);
      y90.push_back(c90);
      ey90.push_back(e90);

      yMin = std::min(yMin, std::min(c70 - e70, std::min(c80 - e80, c90 - e90)));
      yMax = std::max(yMax, std::max(c70 + e70, std::max(c80 + e80, c90 + e90)));
    }

    if (x.empty())
    {
      std::cerr << "[WARN] No iso-efficiency points for " << cb.tag << std::endl;
      continue;
    }

    const double pad = (yMax > yMin) ? 0.70 * (yMax - yMin) : 0.5;
    const double yLo = std::max(0.0, yMin - pad);
    const double yHi = yMax + pad;

    TCanvas c(("c_embeddedPhoton_stitchedIsoCutEfficiency_" + cb.tag).c_str(),
              "stitched embedded photon iso-cut efficiency", 1100, 800);
    c.SetLeftMargin(0.14);
    c.SetRightMargin(0.04);
    c.SetBottomMargin(0.13);
    c.SetTopMargin(0.10);

    TH1F frame(("hFrame_embeddedPhoton_stitchedIsoCutEfficiency_" + cb.tag).c_str(),
               "", 100, 10.0, 35.0);
    frame.SetDirectory(nullptr);
    frame.SetStats(0);
    frame.SetMinimum(yLo);
    frame.SetMaximum(yHi);
    frame.GetXaxis()->SetTitle("Cluster p_{T} [GeV]");
    frame.GetYaxis()->SetTitle("E_{T}^{iso} Cutoff [GeV]");
    frame.GetXaxis()->SetTitleSize(0.060);
    frame.GetYaxis()->SetTitleSize(0.060);
    frame.GetXaxis()->SetLabelSize(0.050);
    frame.GetYaxis()->SetLabelSize(0.050);
    frame.GetYaxis()->SetTitleOffset(1.05);
    frame.Draw();

    TGraphErrors g90(static_cast<int>(x.size()), x.data(), y90.data(), ex.data(), ey90.data());
    TGraphErrors g80(static_cast<int>(x.size()), x.data(), y80.data(), ex.data(), ey80.data());
    TGraphErrors g70(static_cast<int>(x.size()), x.data(), y70.data(), ex.data(), ey70.data());
    StyleEffGraph(g90, kMagenta + 1, 20);
    StyleEffGraph(g80, kGreen + 2, 21);
    StyleEffGraph(g70, kBlue + 1, 22);

    IsoFlatCentResult flat;
    flat.lo = cb.lo;
    flat.hi = cb.hi;
    flat.center = 0.5 * (cb.lo + cb.hi);
    flat.halfWidth = 0.5 * (cb.hi - cb.lo);
    flat.have70 = ComputeFlatGraphAverage(g70, fitXLo, fitXHi, flat.flat70, flat.err70);
    flat.have80 = ComputeFlatGraphAverage(g80, fitXLo, fitXHi, flat.flat80, flat.err80);
    flat.have90 = ComputeFlatGraphAverage(g90, fitXLo, fitXHi, flat.flat90, flat.err90);
    flatResults.push_back(flat);

    TF1 f90(("f_embeddedIsoFlat90_" + cb.tag).c_str(), "[0]", fitXLo, fitXHi);
    TF1 f80(("f_embeddedIsoFlat80_" + cb.tag).c_str(), "[0]", fitXLo, fitXHi);
    TF1 f70(("f_embeddedIsoFlat70_" + cb.tag).c_str(), "[0]", fitXLo, fitXHi);
    if (flat.have90) f90.SetParameter(0, flat.flat90);
    if (flat.have80) f80.SetParameter(0, flat.flat80);
    if (flat.have70) f70.SetParameter(0, flat.flat70);
    f90.SetLineColor(kMagenta + 1);
    f80.SetLineColor(kGreen + 2);
    f70.SetLineColor(kBlue + 1);
    f90.SetLineWidth(3);
    f80.SetLineWidth(3);
    f70.SetLineWidth(3);
    f90.SetLineStyle(2);
    f80.SetLineStyle(2);
    f70.SetLineStyle(2);

    g90.Draw("PE1 SAME");
    g80.Draw("PE1 SAME");
    g70.Draw("PE1 SAME");
    if (flat.have90) f90.Draw("SAME");
    if (flat.have80) f80.Draw("SAME");
    if (flat.have70) f70.Draw("SAME");

    TLegend leg(0.20, 0.16, 0.78, 0.24);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetTextFont(42);
    leg.SetTextSize(0.030);
    leg.AddEntry(&g70, flat.have70 ? TString::Format("70%% Efficiency, E_{T}^{iso} #sim %.2f", flat.flat70).Data() : "70% Efficiency", "ep");
    leg.AddEntry(&g80, flat.have80 ? TString::Format("80%% Efficiency, E_{T}^{iso} #sim %.2f", flat.flat80).Data() : "80% Efficiency", "ep");
    leg.AddEntry(&g90, flat.have90 ? TString::Format("90%% Efficiency, E_{T}^{iso} #sim %.2f", flat.flat90).Data() : "90% Efficiency", "ep");
    leg.SetNColumns(3);
    leg.Draw();

    TLatex info;
    info.SetNDC(true);
    info.SetTextFont(42);
    info.SetTextAlign(13);
    info.SetTextSize(0.030);
    info.DrawLatex(0.18, 0.88, "Photon+Jet 12+20 Embedded SIM");
    info.DrawLatex(0.18, 0.83, TString::Format("%d-%d%% centrality", cb.lo, cb.hi).Data());
    info.DrawLatex(0.18, 0.78, VzEtaLabel().c_str());
    info.DrawLatex(0.18, 0.73, IsoConeLabel().c_str());

    TLatex sph;
    sph.SetNDC(true);
    sph.SetTextFont(42);
    sph.SetTextAlign(33);
    sph.SetTextSize(0.042);
    sph.DrawLatex(0.92, 0.88, "#bf{sPHENIX} #it{Internal}");
    sph.SetTextSize(0.034);
    sph.DrawLatex(0.92, 0.83, "Pythia Overlay  #sqrt{s_{NN}} = 200 GeV");

    const std::string outPng = kOutDir + "/embeddedPhoton_stitchedIsoCutEfficiency_" + cb.tag + ".png";
    c.SaveAs(outPng.c_str());
    std::cout << "[DONE] Wrote " << outPng << std::endl;
  }

  return flatResults;
}

void DrawIsoFlatCutoffVsCentralityQA(const std::vector<IsoFlatCentResult>& flatResults)
{
  if (flatResults.empty())
  {
    std::cerr << "[WARN] No flat iso cutoff values available for centrality summary" << std::endl;
    return;
  }

  std::vector<double> x70, y70, ex70, ey70;
  std::vector<double> x80, y80, ex80, ey80;
  std::vector<double> x90, y90, ex90, ey90;

  double yMin = std::numeric_limits<double>::max();
  double yMax = -std::numeric_limits<double>::max();
  auto addPoint = [&](std::vector<double>& x,
                      std::vector<double>& y,
                      std::vector<double>& ex,
                      std::vector<double>& ey,
                      const IsoFlatCentResult& r,
                      double value,
                      double err)
  {
    x.push_back(r.center);
    y.push_back(value);
    ex.push_back(0.0);
    ey.push_back(err);
    yMin = std::min(yMin, value - err);
    yMax = std::max(yMax, value + err);
  };

  for (const auto& r : flatResults)
  {
    if (r.have70) addPoint(x70, y70, ex70, ey70, r, r.flat70, r.err70);
    if (r.have80) addPoint(x80, y80, ex80, ey80, r, r.flat80, r.err80);
    if (r.have90) addPoint(x90, y90, ex90, ey90, r, r.flat90, r.err90);
  }

  if (x70.empty() && x80.empty() && x90.empty()) return;

  const double pad = (yMax > yMin) ? 0.65 * (yMax - yMin) : 0.75;
  TCanvas c("c_embeddedPhoton_stitchedIsoCutEfficiencyFits_vsCentrality",
            "stitched embedded photon flat iso-cut fits vs centrality", 900, 700);
  c.SetLeftMargin(0.14);
  c.SetRightMargin(0.04);
  c.SetBottomMargin(0.13);
  c.SetTopMargin(0.10);

  TH1F frame("hFrame_embeddedPhoton_stitchedIsoCutEfficiencyFits_vsCentrality",
             "", 100, 0.0, 80.0);
  frame.SetDirectory(nullptr);
  frame.SetStats(0);
  frame.SetMinimum(std::max(0.0, yMin - pad));
  frame.SetMaximum(yMax + pad);
  frame.GetXaxis()->SetTitle("Centrality [%]");
  frame.GetYaxis()->SetTitle("E_{T}^{iso} Cutoff [GeV]");
  frame.GetXaxis()->SetTitleSize(0.055);
  frame.GetYaxis()->SetTitleSize(0.055);
  frame.GetXaxis()->SetLabelSize(0.045);
  frame.GetYaxis()->SetLabelSize(0.045);
  frame.GetYaxis()->SetTitleOffset(1.15);
  frame.Draw();

  TLegend leg(0.48, 0.62, 0.92, 0.78);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextFont(42);
  leg.SetTextSize(0.028);

  auto drawFitGraph = [&](std::vector<double>& x,
                          std::vector<double>& y,
                          std::vector<double>& ex,
                          std::vector<double>& ey,
                          const char* name,
                          const char* label,
                          int color,
                          int marker)
  {
    if (x.empty()) return;
    TGraphErrors* g = new TGraphErrors(static_cast<int>(x.size()), x.data(), y.data(), ex.data(), ey.data());
    g->SetLineWidth(2);
    g->SetLineColor(color);
    g->SetMarkerColor(color);
    g->SetMarkerStyle(marker);
    g->SetMarkerSize(marker == 22 ? 1.5 : 1.2);
    g->Draw("P SAME");

    TF1* fit = new TF1((std::string("fit_") + name).c_str(), "pol1", 0.0, 80.0);
    fit->SetLineColor(color);
    fit->SetLineWidth(2);
    fit->SetLineStyle(2);
    g->Fit(fit, "QNR");
    fit->Draw("SAME");
    leg.AddEntry(g, TString::Format("%s: y = %.4fx %+.2f", label, fit->GetParameter(1), fit->GetParameter(0)).Data(), "lp");
  };

  drawFitGraph(x90, y90, ex90, ey90, "90", "90% Eff", kMagenta + 1, 20);
  drawFitGraph(x80, y80, ex80, ey80, "80", "80% Eff", kGreen + 2, 21);
  drawFitGraph(x70, y70, ex70, ey70, "70", "70% Eff", kBlue + 1, 22);
  leg.Draw();

  TLatex title;
  title.SetNDC(true);
  title.SetTextFont(42);
  title.SetTextAlign(23);
  title.SetTextSize(0.042);
  title.DrawLatex(0.50, 0.98, "Flat E_{T}^{iso} cutoff vs centrality");

  TLatex info;
  info.SetNDC(true);
  info.SetTextFont(42);
  info.SetTextAlign(13);
  info.SetTextSize(0.030);
  info.DrawLatex(0.18, 0.88, "Photon+Jet 12+20 Embedded SIM");
  info.DrawLatex(0.18, 0.83, "constant fit over full plotted p_{T}^{#gamma} range");
  info.DrawLatex(0.18, 0.78, VzEtaLabel().c_str());
  info.DrawLatex(0.18, 0.73, IsoConeLabel().c_str());

  TLatex sph;
  sph.SetNDC(true);
  sph.SetTextFont(42);
  sph.SetTextAlign(33);
  sph.SetTextSize(0.042);
  sph.DrawLatex(0.92, 0.88, "#bf{sPHENIX} #it{Internal}");
  sph.SetTextSize(0.034);
  sph.DrawLatex(0.92, 0.83, "Pythia Overlay  #sqrt{s_{NN}} = 200 GeV");

  const std::string outPng = kOutDir + "/embeddedPhoton_stitchedIsoCutEfficiencyFits_vsCentrality.png";
  c.SaveAs(outPng.c_str());
  std::cout << "[DONE] Wrote " << outPng << std::endl;
}

void DrawIsoCentralityOverlay10To12QA(TFile* f12, TFile* f20, double w12, double w20)
{
  struct CentBin
  {
    int lo;
    int hi;
    std::string suffix;
    std::string tag;
    int color;
  };

  const std::vector<CentBin> centBins = {
      {0, 20, "_cent_0_20", "cent_0_20", kBlack},
      {20, 50, "_cent_20_50", "cent_20_50", kBlue + 1},
      {50, 80, "_cent_50_80", "cent_50_80", kOrange + 1},
  };

  std::vector<std::unique_ptr<TH1>> hOwned;
  std::vector<TH1*> hCents;
  std::vector<std::string> labels;
  double yMax = 0.0;

  for (const auto& cb : centBins)
  {
    const std::string histName = "h_Eiso_pT_10_12" + cb.suffix;
    std::unique_ptr<TH1> h = BuildWeightedPairHistogram(
        f12, f20, histName, "h_combinedIsoCentOverlay_10_12_" + cb.tag, w12, w20);
    if (!h || h->Integral(1, h->GetNbinsX()) <= 0.0)
    {
      std::cerr << "[WARN] Missing or empty " << histName << " for iso centrality overlay" << std::endl;
      continue;
    }

    h->Rebin(10);
    if (h->GetSumw2N() == 0) h->Sumw2();
    const double integral = h->Integral(1, h->GetNbinsX());
    if (integral <= 0.0) continue;
    h->Scale(1.0 / integral);

    h->SetTitle("");
    h->SetStats(0);
    h->SetLineColor(cb.color);
    h->SetMarkerColor(cb.color);
    h->SetMarkerStyle(20);
    h->SetMarkerSize(0.9);
    h->SetLineWidth(2);
    h->SetFillStyle(0);

    yMax = std::max(yMax, h->GetMaximum());
    labels.push_back(std::to_string(cb.lo) + "-" + std::to_string(cb.hi) + "%");
    hCents.push_back(h.get());
    hOwned.push_back(std::move(h));
  }

  if (hCents.empty())
  {
    std::cerr << "[WARN] No centrality histograms available for embeddedPhoton_stitchedIsoCentralityOverlay_pT_10_12.png" << std::endl;
    return;
  }

  TCanvas c("c_embeddedPhoton_stitchedIsoCentralityOverlay_pT_10_12",
            "stitched embedded photon isolation centrality overlay", 900, 700);
  c.SetLeftMargin(0.14);
  c.SetRightMargin(0.04);
  c.SetBottomMargin(0.13);
  c.SetTopMargin(0.10);

  hCents[0]->GetXaxis()->SetTitle("E_{T}^{iso} [GeV]");
  hCents[0]->GetYaxis()->SetTitle("Normalized to Unit Area");
  hCents[0]->GetXaxis()->SetTitleSize(0.048);
  hCents[0]->GetYaxis()->SetTitleSize(0.050);
  hCents[0]->GetXaxis()->SetLabelSize(0.040);
  hCents[0]->GetYaxis()->SetLabelSize(0.050);
  hCents[0]->GetYaxis()->SetTitleOffset(1.15);
  hCents[0]->GetXaxis()->SetRangeUser(-10.0, 50.0);
  hCents[0]->SetMinimum(0.0);
  hCents[0]->SetMaximum((yMax > 0.0) ? 1.25 * yMax : 1.0);
  hCents[0]->Draw("E1");
  for (std::size_t ih = 1; ih < hCents.size(); ++ih)
  {
    hCents[ih]->Draw("E1 SAME");
  }

  TLegend leg(0.56, 0.39, 0.89, 0.61);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextFont(42);
  leg.SetTextSize(0.038);
  leg.SetNColumns(2);
  for (std::size_t ih = 0; ih < hCents.size(); ++ih)
  {
    leg.AddEntry(hCents[ih], labels[ih].c_str(), "ep");
  }
  leg.Draw();

  TLatex title;
  title.SetNDC(true);
  title.SetTextFont(42);
  title.SetTextAlign(23);
  title.SetTextSize(0.052);
  title.DrawLatex(0.50, 0.99, "Embedded photon centrality overlays, p_{T}^{#gamma} = 10-12 GeV");

  TLatex info;
  info.SetNDC(true);
  info.SetTextFont(42);
  info.SetTextAlign(33);
  info.SetTextSize(0.040);
  info.DrawLatex(0.90, 0.88, "PhotonJet12+20 stitched SIM");
  info.DrawLatex(0.90, 0.84, "with UE Sub");
  info.DrawLatex(0.90, 0.80, IsoConeLabel().c_str());

  TLatex sph;
  sph.SetNDC(true);
  sph.SetTextFont(42);
  sph.SetTextAlign(33);
  sph.SetTextSize(0.052);
  sph.DrawLatex(0.90, 0.32, "#bf{sPHENIX} #it{Internal}");
  sph.SetTextSize(0.042);
  sph.DrawLatex(0.90, 0.26, "Au+Au  #sqrt{s_{NN}} = 200 GeV");

  const std::string outPng = kOutDir + "/embeddedPhoton_stitchedIsoCentralityOverlay_pT_10_12.png";
  c.SaveAs(outPng.c_str());
  std::cout << "[DONE] Wrote " << outPng << std::endl;
}

void DrawSpectrumSmoothQA(std::unique_ptr<TH1> h12,
                          std::unique_ptr<TH1> h20,
                          const std::string& outputName,
                          const std::string& title,
                          const std::string& xTitle,
                          const std::string& note,
                          double w12,
                          double w20)
{
  if (!h12 || !h20)
  {
    std::cerr << "[WARN] Missing histogram(s) for " << outputName << std::endl;
    return;
  }

  const bool isFilterPtQA = (outputName == "embeddedPhoton_stitchedTruthFilterPtSpectrum");
  const double xMin = isFilterPtQA ? 10.0 : h12->GetXaxis()->GetXmin();
  const double xMax = isFilterPtQA ? 45.0 : h12->GetXaxis()->GetXmax();
  const double ratioMin = isFilterPtQA ? 0.85 : 0.45;
  const double ratioMax = isFilterPtQA ? 1.15 : 1.55;

  h12->SetTitle("");
  h20->SetTitle("");
  if (isFilterPtQA)
  {
    h12->Rebin(2);
    h20->Rebin(2);
  }

  StyleHist(h12.get(), kBlue + 1, 20);
  StyleHist(h20.get(), kRed + 1, 21);

  std::unique_ptr<TH1> hSum(dynamic_cast<TH1*>(h12->Clone((outputName + "_sum").c_str())));
  hSum->SetDirectory(nullptr);
  hSum->Add(h20.get());
  hSum->SetTitle("");
  StyleHist(hSum.get(), kBlack, 24);

  std::unique_ptr<TH1> hSmooth = MakeSmoothReference(hSum.get(), outputName + "_smooth");
  std::unique_ptr<TH1> hRatio = MakeRatioToSmooth(hSum.get(), hSmooth.get(), outputName + "_ratio");
  hRatio->SetTitle("");
  StyleHist(hRatio.get(), kBlack, 20);

  hSum->GetXaxis()->SetTitle(xTitle.c_str());
  hSum->GetXaxis()->SetRangeUser(xMin, xMax);
  hSum->GetYaxis()->SetTitle("#sigma_{eff}/N scaled entries [pb / bin]");
  hSum->GetYaxis()->SetTitleOffset(1.12);

  hRatio->GetXaxis()->SetTitle(xTitle.c_str());
  hRatio->GetXaxis()->SetRangeUser(xMin, xMax);
  hRatio->GetYaxis()->SetTitle("sum / smooth");
  hRatio->GetYaxis()->SetRangeUser(ratioMin, ratioMax);
  hRatio->GetYaxis()->SetNdivisions(505);
  hRatio->GetYaxis()->SetTitleSize(0.085);
  hRatio->GetYaxis()->SetLabelSize(0.075);
  hRatio->GetYaxis()->SetTitleOffset(0.50);
  hRatio->GetXaxis()->SetTitleSize(0.090);
  hRatio->GetXaxis()->SetLabelSize(0.080);

  TCanvas c(("c_" + outputName).c_str(), outputName.c_str(), 1050, 850);
  TPad top("top", "top", 0.0, 0.30, 1.0, 1.0);
  TPad bot("bot", "bot", 0.0, 0.0, 1.0, 0.31);
  top.SetBottomMargin(0.025);
  top.SetLeftMargin(0.12);
  top.SetRightMargin(0.04);
  top.SetLogy();
  bot.SetTopMargin(0.03);
  bot.SetBottomMargin(0.28);
  bot.SetLeftMargin(0.12);
  bot.SetRightMargin(0.04);
  top.Draw();
  bot.Draw();

  top.cd();
  const double ymax = std::max({hSum->GetMaximum(), h12->GetMaximum(), h20->GetMaximum()});
  h12->GetXaxis()->SetRangeUser(xMin, xMax);
  h20->GetXaxis()->SetRangeUser(xMin, xMax);
  hSum->SetMinimum(std::max(1.0e-8, ymax * 2.0e-5));
  hSum->SetMaximum(std::max(1.0e-6, ymax * (isFilterPtQA ? 85.0 : 18.0)));
  hSum->GetXaxis()->SetLabelSize(0.0);
  hSum->GetXaxis()->SetTitleSize(0.0);
  hSum->Draw("E1");
  h12->Draw("E1 SAME");
  h20->Draw("E1 SAME");
  if (hSmooth)
  {
    hSmooth->SetTitle("");
    hSmooth->GetXaxis()->SetRangeUser(xMin, xMax);
    hSmooth->SetLineColor(kGray + 2);
    hSmooth->SetLineStyle(2);
    hSmooth->SetLineWidth(2);
    hSmooth->SetMarkerSize(0);
    hSmooth->Draw("HIST SAME");
  }

  TLegend leg(isFilterPtQA ? 0.15 : 0.55,
              isFilterPtQA ? 0.16 : 0.60,
              isFilterPtQA ? 0.50 : 0.91,
              isFilterPtQA ? 0.41 : 0.87);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextSize(0.038);
  leg.AddEntry(h12.get(), isFilterPtQA ? "PhotonJet12 stitched" : "weighted PhotonJet12", "lep");
  leg.AddEntry(h20.get(), isFilterPtQA ? "PhotonJet20 stitched" : "weighted PhotonJet20", "lep");
  leg.AddEntry(hSum.get(), isFilterPtQA ? "Combined" : "weighted sum", "lep");
  leg.AddEntry(hSmooth.get(), isFilterPtQA ? "Smoothed reference" : "smoothed sum reference", "l");
  leg.Draw();

  TLatex lat;
  lat.SetNDC();
  lat.SetTextSize(0.040);
  lat.DrawLatex(0.15, 0.86, title.c_str());
  lat.SetTextSize(0.031);
  lat.DrawLatex(0.15, 0.80, note.c_str());
  if (isFilterPtQA)
  {
    lat.DrawLatex(0.15, 0.75, "PhotonJet12: 12 #leq p_{T,filter}^{#gamma} < 20 GeV; PhotonJet20: p_{T,filter}^{#gamma} #geq 20 GeV");
    lat.DrawLatex(0.15, 0.70, ("w_{12#rightarrow20}=" + Sci(w12) + " pb/event, w_{20+}=" + Sci(w20) + " pb/event").c_str());
  }
  else
  {
    lat.DrawLatex(0.15, 0.75, ("w_{12#rightarrow20}=" + Sci(w12) + " pb/event, w_{20+}=" + Sci(w20) + " pb/event").c_str());
  }

  if (isFilterPtQA)
  {
    TLatex tSphM;
    tSphM.SetNDC(true);
    tSphM.SetTextFont(42);
    tSphM.SetTextAlign(33);
    tSphM.SetTextSize(0.042);
    tSphM.DrawLatex(0.92, 0.58, "#bf{sPHENIX} #it{Internal}");
    tSphM.SetTextSize(0.034);
    tSphM.DrawLatex(0.92, 0.53, "Pythia Overlay #sqrt{s_{NN}} = 200 GeV");
  }

  bot.cd();
  hRatio->Draw("E1");
  TLine one(xMin, 1.0, xMax, 1.0);
  one.SetLineColor(kGray + 1);
  one.SetLineStyle(2);
  one.Draw("SAME");

  const std::string outPng = kOutDir + "/" + outputName + ".png";
  c.SaveAs(outPng.c_str());
  std::cout << "[DONE] Wrote " << outPng << std::endl;
}

void DrawCompositionQA(std::unique_ptr<TH1> h12,
                       std::unique_ptr<TH1> h20,
                       double w12,
                       double w20)
{
  if (!h12 || !h20)
  {
    std::cerr << "[WARN] Missing ABCD histograms for composition QA" << std::endl;
    return;
  }

  StyleHist(h12.get(), kBlue + 1, 20);
  StyleHist(h20.get(), kRed + 1, 21);

  std::unique_ptr<TH1> hSum(dynamic_cast<TH1*>(h12->Clone("h_composition_sum")));
  hSum->SetDirectory(nullptr);
  hSum->Add(h20.get());
  StyleHist(hSum.get(), kBlack, 24);

  std::unique_ptr<TH1> frac12(dynamic_cast<TH1*>(h12->Clone("h_fractionPhoton12")));
  std::unique_ptr<TH1> frac20(dynamic_cast<TH1*>(h20->Clone("h_fractionPhoton20")));
  frac12->SetDirectory(nullptr);
  frac20->SetDirectory(nullptr);
  frac12->Divide(hSum.get());
  frac20->Divide(hSum.get());
  StyleHist(frac12.get(), kBlue + 1, 20);
  StyleHist(frac20.get(), kRed + 1, 21);

  hSum->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
  hSum->GetYaxis()->SetTitle("#sigma_{eff}/N scaled entries [pb / bin]");
  hSum->GetYaxis()->SetTitleOffset(1.12);

  frac12->GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
  frac12->GetYaxis()->SetTitle("sample fraction");
  frac12->GetYaxis()->SetRangeUser(0.0, 1.08);
  frac12->GetYaxis()->SetNdivisions(505);
  frac12->GetYaxis()->SetTitleSize(0.085);
  frac12->GetYaxis()->SetLabelSize(0.075);
  frac12->GetYaxis()->SetTitleOffset(0.50);
  frac12->GetXaxis()->SetTitleSize(0.090);
  frac12->GetXaxis()->SetLabelSize(0.080);

  TCanvas c("c_embeddedPhoton_stitchedSampleComposition_ABCDsum",
            "embedded photon stitched sample composition ABCD sum", 1050, 850);
  TPad top("top", "top", 0.0, 0.30, 1.0, 1.0);
  TPad bot("bot", "bot", 0.0, 0.0, 1.0, 0.31);
  top.SetBottomMargin(0.025);
  top.SetLeftMargin(0.12);
  top.SetRightMargin(0.04);
  top.SetLogy();
  bot.SetTopMargin(0.03);
  bot.SetBottomMargin(0.28);
  bot.SetLeftMargin(0.12);
  bot.SetRightMargin(0.04);
  top.Draw();
  bot.Draw();

  top.cd();
  const double ymax = std::max({hSum->GetMaximum(), h12->GetMaximum(), h20->GetMaximum()});
  hSum->SetMinimum(std::max(1.0e-8, ymax * 2.0e-5));
  hSum->SetMaximum(std::max(1.0e-6, ymax * 18.0));
  hSum->Draw("E1");
  h12->Draw("E1 SAME");
  h20->Draw("E1 SAME");

  TLegend leg(0.55, 0.63, 0.91, 0.87);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextSize(0.040);
  leg.AddEntry(hSum.get(), "weighted ABCD sum", "lep");
  leg.AddEntry(h12.get(), "weighted PhotonJet12 ABCD", "lep");
  leg.AddEntry(h20.get(), "weighted PhotonJet20 ABCD", "lep");
  leg.Draw();

  TLatex lat;
  lat.SetNDC();
  lat.SetTextSize(0.040);
  lat.DrawLatex(0.15, 0.86, "Embedded photon stitched sample composition");
  lat.SetTextSize(0.032);
  lat.DrawLatex(0.15, 0.80, "Reference ID; ABCD-summed reco photon spectrum");
  lat.DrawLatex(0.15, 0.75, ("w_{12#rightarrow20}=" + Sci(w12) + " pb/event, w_{20+}=" + Sci(w20) + " pb/event").c_str());

  bot.cd();
  frac12->Draw("E1");
  frac20->Draw("E1 SAME");
  TLine half(frac12->GetXaxis()->GetXmin(), 0.5, frac12->GetXaxis()->GetXmax(), 0.5);
  half.SetLineColor(kGray + 1);
  half.SetLineStyle(2);
  half.Draw("SAME");

  TLegend leg2(0.58, 0.70, 0.91, 0.93);
  leg2.SetBorderSize(0);
  leg2.SetFillStyle(0);
  leg2.SetTextSize(0.075);
  leg2.AddEntry(frac12.get(), "PhotonJet12 / sum", "lep");
  leg2.AddEntry(frac20.get(), "PhotonJet20 / sum", "lep");
  leg2.Draw();

  const std::string outPng = kOutDir + "/embeddedPhoton_stitchedSampleComposition_ABCDsum.png";
  c.SaveAs(outPng.c_str());
  std::cout << "[DONE] Wrote " << outPng << std::endl;

  const std::string legacyNormPng = kOutDir + "/embeddedPhoton12to20_plusPhoton20_xsecNormalizationQA.png";
  c.SaveAs(legacyNormPng.c_str());
  std::cout << "[DONE] Wrote " << legacyNormPng << std::endl;
}
}

void MakeEmbeddedPhotonXsecNormQA()
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetErrorX(0.5);

  std::unique_ptr<TFile> f12(TFile::Open(kFile12.c_str(), "READ"));
  std::unique_ptr<TFile> f20(TFile::Open(kFile20.c_str(), "READ"));
  if (!f12 || f12->IsZombie())
  {
    std::cerr << "[ERROR] Cannot open " << kFile12 << std::endl;
    return;
  }
  if (!f20 || f20->IsZombie())
  {
    std::cerr << "[ERROR] Cannot open " << kFile20 << std::endl;
    return;
  }

  if (kSigmaEmbeddedPhoton12To20_pb <= 0.0 || kSigmaEmbeddedPhoton20Plus_pb <= 0.0)
  {
    std::cerr << "[ERROR] Embedded photon cross sections are not set." << std::endl;
    return;
  }

  const double n12 = ReadEventCount(f12.get());
  const double n20 = ReadEventCount(f20.get());
  if (n12 <= 0.0 || n20 <= 0.0)
  {
    std::cerr << "[ERROR] Bad event counts: N12=" << n12 << " N20=" << n20 << std::endl;
    return;
  }

  const double w12 = kSigmaEmbeddedPhoton12To20_pb / n12;
  const double w20 = kSigmaEmbeddedPhoton20Plus_pb / n20;

  gSystem->mkdir(kOutDir.c_str(), true);

  DrawSpectrumSmoothQA(
      WeightedClone(BuildKeptFilterPtSpectrum(f12.get(), "h_filterPt12_weighted"), w12),
      WeightedClone(BuildKeptFilterPtSpectrum(f20.get(), "h_filterPt20_weighted"), w20),
      "embeddedPhoton_stitchedTruthFilterPtSpectrum",
      "Embedded photon generator stitching spectrum",
      "p_{T,filter}^{#gamma} [GeV]",
      "Uses SIM/h_embedStitch_filterPhotonPt_kept",
      w12,
      w20);

  const std::vector<IsoFlatCentResult> flatIsoResults = DrawIsoEfficiencyCutoffQA(f12.get(), f20.get(), w12, w20);
  DrawIsoFlatCutoffVsCentralityQA(flatIsoResults);
  DrawIsoCentralityOverlay10To12QA(f12.get(), f20.get(), w12, w20);

  std::cout << "[INFO] N12=" << n12 << " N20=" << n20
            << " w12=" << w12 << " pb/event"
            << " w20=" << w20 << " pb/event" << std::endl;
}
