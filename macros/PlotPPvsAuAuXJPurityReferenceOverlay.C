#include "AnalyzeRecoilJets.h"

#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2.h>
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
#include <string>
#include <vector>

namespace
{
constexpr double kSigmaEmbeddedPhoton12To20_pb = 2598.12425;
constexpr double kSigmaEmbeddedPhoton20Plus_pb = 133.317866;

const std::string kPPTag =
    "jetMinPt5_7pi_8_vz60_isoR40_fixedIso2GeV_preselectionReference_tightReference_nonTightReference";
const std::string kAATag =
    "jetMinPt5_7pi_8_vz60_isoR40_isSliding_baseVariant_preselectionReference_tightReference_nonTightReference";
const std::string kComboPP = "photonJet5and10and20merged_SIM";
const std::string kComboAA = "photonJet12and20merged_SIM";
const std::string kMergedPPName = "RecoilJets_photonjet5plus10plus20_MERGED.root";
const std::string kMergedAAName = "RecoilJets_embeddedPhoton12plus20_MERGED.root";

struct Series
{
  std::string label;
  std::string suffix;
  int color = kBlack;
  int marker = 20;
  bool isPP = false;
};

std::string PPInput(const std::string& sample)
{
  return ARJ::kInputBase + "/simPhotonJet/RecoilJets_" + sample + "_ALL_" + kPPTag + ".root";
}

std::string AAInput(const std::string& sample)
{
  return ARJ::kInputBase + "/simEmbedded/RecoilJets_" + sample + "_ALL_" + kAATag + ".root";
}

std::string PPMergedPath()
{
  return ARJ::kOutputBase + "/combinedSimOnly/" + kPPTag + "/" + kComboPP + "/" + kMergedPPName;
}

std::string AAMergedPath()
{
  return ARJ::kInputBase + "/simEmbedded/merged/" + kAATag + "/" + kComboAA + "/" + kMergedAAName;
}

std::string PPDataPath()
{
  return ARJ::kInputBase + "/pp24/RecoilJets_pp_ALL_" + kPPTag + ".root";
}

std::string AADataPath()
{
  return ARJ::kInputBase + "/auau25/RecoilJets_auau_ALL_" + kAATag + ".root";
}

std::string OutputDir()
{
  return ARJ::kOutputBase + "/combinedSimOnlyEMBEDDED/" + kAATag + "/" + kComboAA;
}

std::string Sanitize(const std::string& s)
{
  std::string out;
  out.reserve(s.size());
  for (const char ch : s)
  {
    if ((ch >= 'a' && ch <= 'z') ||
        (ch >= 'A' && ch <= 'Z') ||
        (ch >= '0' && ch <= '9'))
    {
      out.push_back(ch);
    }
    else
    {
      out.push_back('_');
    }
  }
  return out;
}

bool RebuildInputs()
{
  const bool okPP = ARJ::BuildMergedSIMFile_PhotonSlices(
      {PPInput("photonjet5"), PPInput("photonjet10"), PPInput("photonjet20")},
      {ARJ::kSigmaPhoton5_pb, ARJ::kSigmaPhoton10_pb, ARJ::kSigmaPhoton20_pb},
      PPMergedPath(),
      ARJ::kDirSIM,
      {"photonJet5", "photonJet10", "photonJet20"});

  const bool okAA = ARJ::BuildMergedSIMFile_PhotonSlices(
      {AAInput("embeddedPhoton12"), AAInput("embeddedPhoton20")},
      {kSigmaEmbeddedPhoton12To20_pb, kSigmaEmbeddedPhoton20Plus_pb},
      AAMergedPath(),
      ARJ::kDirSIM,
      {"embeddedPhoton12", "embeddedPhoton20"});

  return okPP && okAA;
}

void StyleHist(TH1* h, const Series& s)
{
  if (!h) return;
  h->SetDirectory(nullptr);
  h->SetStats(0);
  h->SetLineColor(s.color);
  h->SetMarkerColor(s.color);
  h->SetMarkerStyle(s.marker);
  h->SetMarkerSize(1.05);
  h->SetLineWidth(2);
}

std::unique_ptr<TH1D> ProjectXJ(TDirectory* dir,
                                const Series& s,
                                const double ptMin,
                                const double ptMax,
                                const bool maskFirstXJBin)
{
  if (!dir) return nullptr;

  const std::string hname = "h2_unfoldReco_pTgamma_xJ_incl_r03" + s.suffix;
  TH2* h2 = dynamic_cast<TH2*>(dir->Get(hname.c_str()));
  if (!h2)
  {
    std::cerr << "[WARN] Missing " << hname << "\n";
    return nullptr;
  }

  const int ixLo = h2->GetXaxis()->FindBin(ptMin + 1e-6);
  const int ixHi = h2->GetXaxis()->FindBin(ptMax - 1e-6);
  std::unique_ptr<TH1D> h(h2->ProjectionY(
      TString::Format("h_xJ_%s_%.0f_%.0f", s.label.c_str(), ptMin, ptMax).Data(),
      ixLo, ixHi, "e"));
  if (!h) return nullptr;

  StyleHist(h.get(), s);
  h->GetXaxis()->SetRangeUser(0.0, 2.0);
  if (maskFirstXJBin && h->GetNbinsX() >= 1)
  {
    h->SetBinContent(1, 0.0);
    h->SetBinError(1, 0.0);
  }
  const double integral = h->Integral(0, h->GetNbinsX() + 1);
  if (integral > 0.0) h->Scale(1.0 / integral);
  return h;
}

double RawPurityError(double A,
                      double B,
                      double C,
                      double D,
                      double eA,
                      double eB,
                      double eC,
                      double eD)
{
  if (A <= 0.0 || D <= 0.0) return 0.0;
  const double dPdA =  (B * C) / (A * A * D);
  const double dPdB = -(C) / (A * D);
  const double dPdC = -(B) / (A * D);
  const double dPdD =  (B * C) / (A * D * D);
  double var = 0.0;
  var += dPdA * dPdA * eA * eA;
  var += dPdB * dPdB * eB * eB;
  var += dPdC * dPdC * eC * eC;
  var += dPdD * dPdD * eD * eD;
  return (var > 0.0) ? std::sqrt(var) : 0.0;
}

double LeakageCorrectedPurityValue(const double A,
                                   const double B,
                                   const double C,
                                   const double D,
                                   const double fB,
                                   const double fC,
                                   const double fD)
{
  if (A <= 0.0) return 0.0;
  double SA = 0.0;
  const bool ok = ARJ::SolveLeakageCorrectedSA(A, B, C, D, fB, fC, fD, SA);
  if (ok) return SA / A;

  if (D <= 0.0) return 0.0;
  return std::max(0.0, A - B * (C / D)) / A;
}

double LeakageCorrectedPurityError(const double A,
                                   const double B,
                                   const double C,
                                   const double D,
                                   const double eA,
                                   const double eB,
                                   const double eC,
                                   const double eD,
                                   const double fB,
                                   const double fC,
                                   const double fD)
{
  if (A <= 0.0) return 0.0;

  auto Value = [&](const double a, const double b, const double c, const double d) {
    return LeakageCorrectedPurityValue(a, b, c, d, fB, fC, fD);
  };

  const double dA = std::max(eA, 0.0);
  const double dB = std::max(eB, 0.0);
  const double dC = std::max(eC, 0.0);
  const double dD = std::max(eD, 0.0);

  auto Deriv = [&](const double up, const double dn) {
    return (up > dn) ? 1.0 / (up - dn) : 0.0;
  };

  const double Aup = A + dA;
  const double Adn = std::max(0.0, A - dA);
  const double Bup = B + dB;
  const double Bdn = std::max(0.0, B - dB);
  const double Cup = C + dC;
  const double Cdn = std::max(0.0, C - dC);
  const double Dup = D + dD;
  const double Ddn = std::max(0.0, D - dD);

  const double pA = (Value(Aup, B, C, D) - Value(Adn, B, C, D)) * Deriv(Aup, Adn);
  const double pB = (Value(A, Bup, C, D) - Value(A, Bdn, C, D)) * Deriv(Bup, Bdn);
  const double pC = (Value(A, B, Cup, D) - Value(A, B, Cdn, D)) * Deriv(Cup, Cdn);
  const double pD = (Value(A, B, C, Dup) - Value(A, B, C, Ddn)) * Deriv(Dup, Ddn);

  const double var = pA * pA * eA * eA + pB * pB * eB * eB + pC * pC * eC * eC + pD * pD * eD * eD;
  return (var > 0.0) ? std::sqrt(var) : 0.0;
}

std::unique_ptr<TGraphErrors> BuildPurityGraph(TDirectory* obsDir,
                                               TDirectory* leakageDir,
                                               const Series& s,
                                               const bool leakageCorrected)
{
  if (!obsDir) return nullptr;

  std::vector<double> x;
  std::vector<double> ex;
  std::vector<double> y;
  std::vector<double> ey;

  for (int i = 0; i < ARJ::kNPtBins; ++i)
  {
    const ARJ::PtBin& b = ARJ::PtBins()[i];

    auto GetHist = [&](const std::string& base) -> TH1* {
      TH1* h = dynamic_cast<TH1*>(obsDir->Get((base + b.suffix + s.suffix).c_str()));
      return h;
    };

    TH1* hA = GetHist("h_isIsolated_isTight");
    TH1* hB = GetHist("h_notIsolated_isTight");
    TH1* hC = GetHist("h_isIsolated_notTight");
    TH1* hD = GetHist("h_notIsolated_notTight");

    const double A = hA ? hA->GetBinContent(1) : 0.0;
    const double B = hB ? hB->GetBinContent(1) : 0.0;
    const double C = hC ? hC->GetBinContent(1) : 0.0;
    const double D = hD ? hD->GetBinContent(1) : 0.0;
    const double eA = hA ? hA->GetBinError(1) : 0.0;
    const double eB = hB ? hB->GetBinError(1) : 0.0;
    const double eC = hC ? hC->GetBinError(1) : 0.0;
    const double eD = hD ? hD->GetBinError(1) : 0.0;
    if (A <= 0.0) continue;

    double fB = 0.0;
    double fC = 0.0;
    double fD = 0.0;
    if (leakageCorrected)
    {
      TH1* hLeak = leakageDir
          ? dynamic_cast<TH1*>(leakageDir->Get(("h_sigABCD_MC" + b.suffix + s.suffix).c_str()))
          : nullptr;
      if (hLeak)
      {
        const double AsigMC = hLeak->GetBinContent(1);
        if (AsigMC > 0.0)
        {
          fB = hLeak->GetBinContent(2) / AsigMC;
          fC = hLeak->GetBinContent(3) / AsigMC;
          fD = hLeak->GetBinContent(4) / AsigMC;
        }
      }
    }

    double purity = 0.0;
    double purityErr = 0.0;
    if (leakageCorrected)
    {
      purity = LeakageCorrectedPurityValue(A, B, C, D, fB, fC, fD);
      purityErr = LeakageCorrectedPurityError(A, B, C, D, eA, eB, eC, eD, fB, fC, fD);
    }
    else if (D > 0.0)
    {
      const double Asig = std::max(0.0, A - B * (C / D));
      purity = Asig / A;
      purityErr = RawPurityError(A, B, C, D, eA, eB, eC, eD);
    }

    const double ptLo = ARJ::kPtEdges[(std::size_t)i];
    const double ptHi = ARJ::kPtEdges[(std::size_t)i + 1];
    x.push_back(0.5 * (ptLo + ptHi));
    ex.push_back(0.5 * (ptHi - ptLo));
    y.push_back(purity);
    ey.push_back(purityErr);
  }

  if (x.empty()) return nullptr;

  std::unique_ptr<TGraphErrors> g(new TGraphErrors((int)x.size(), x.data(), y.data(), ex.data(), ey.data()));
  g->SetName(("g_purity_" + s.label).c_str());
  g->SetLineColor(s.color);
  g->SetMarkerColor(s.color);
  g->SetMarkerStyle(s.marker);
  g->SetMarkerSize(1.0);
  g->SetLineWidth(2);
  return g;
}

void AuditPuritySeries(TDirectory* obsDir,
                       TDirectory* leakageDir,
                       const Series& s,
                       const std::string& obsPath,
                       const std::string& obsDirName,
                       const std::string& leakagePath,
                       const std::string& leakageDirName)
{
  if (!obsDir) return;

  const std::string outPath = OutputDir() + "/purityLeakageAudit_" + Sanitize(s.label) + ".txt";
  std::ofstream out(outPath.c_str());
  if (!out)
  {
    std::cerr << "[WARN] Could not open audit file for writing: " << outPath << "\n";
    return;
  }

  out << "Purity/leakage correction audit for " << s.label << "\n";
  out << "Observed ABCD DATA file: " << obsPath << "\n";
  out << "Observed ABCD DATA directory: " << obsDirName << "\n";
  out << "Leakage SIM file: " << leakagePath << "\n";
  out << "Leakage SIM directory: " << leakageDirName << "\n";
  out << "Formula raw: Sraw = max(0, A - B*C/D), Praw = Sraw/A\n";
  out << "Formula corrected: solve S_A = A - (B - fB*S_A)*(C - fC*S_A)/(D - fD*S_A), Pcorr = S_A/A\n";
  out << "Notes: A/B/C/D are DATA counts. fB/fC/fD are from SIM h_sigABCD_MC bin ratios to A_sig_MC.\n\n";

  out << std::scientific << std::setprecision(8);
  out
      << "idx ptLo ptHi suffix "
      << "A eA B eB C eC D eD "
      << "rawSub Sraw Praw PrawErr "
      << "AsigMC BsigMC CsigMC DsigMC fB fC fD "
      << "solverOk Scorr Pcorr PcorrErr "
      << "leakHistFound\n";

  std::cout << "[LEAKAGE AUDIT] " << s.label << " -> " << outPath << "\n";

  for (int i = 0; i < ARJ::kNPtBins; ++i)
  {
    const ARJ::PtBin& b = ARJ::PtBins()[i];

    auto GetHist = [&](const std::string& base) -> TH1* {
      return dynamic_cast<TH1*>(obsDir->Get((base + b.suffix + s.suffix).c_str()));
    };

    TH1* hA = GetHist("h_isIsolated_isTight");
    TH1* hB = GetHist("h_notIsolated_isTight");
    TH1* hC = GetHist("h_isIsolated_notTight");
    TH1* hD = GetHist("h_notIsolated_notTight");

    const double A = hA ? hA->GetBinContent(1) : 0.0;
    const double B = hB ? hB->GetBinContent(1) : 0.0;
    const double C = hC ? hC->GetBinContent(1) : 0.0;
    const double D = hD ? hD->GetBinContent(1) : 0.0;
    const double eA = hA ? hA->GetBinError(1) : 0.0;
    const double eB = hB ? hB->GetBinError(1) : 0.0;
    const double eC = hC ? hC->GetBinError(1) : 0.0;
    const double eD = hD ? hD->GetBinError(1) : 0.0;

    const double rawSub = (D > 0.0) ? B * (C / D) : 0.0;
    const double Sraw = (A > 0.0 && D > 0.0) ? std::max(0.0, A - rawSub) : 0.0;
    const double Praw = (A > 0.0) ? Sraw / A : 0.0;
    const double PrawErr = RawPurityError(A, B, C, D, eA, eB, eC, eD);

    TH1* hLeak = leakageDir
        ? dynamic_cast<TH1*>(leakageDir->Get(("h_sigABCD_MC" + b.suffix + s.suffix).c_str()))
        : nullptr;
    const double AsigMC = hLeak ? hLeak->GetBinContent(1) : 0.0;
    const double BsigMC = hLeak ? hLeak->GetBinContent(2) : 0.0;
    const double CsigMC = hLeak ? hLeak->GetBinContent(3) : 0.0;
    const double DsigMC = hLeak ? hLeak->GetBinContent(4) : 0.0;
    const double fB = (AsigMC > 0.0) ? BsigMC / AsigMC : 0.0;
    const double fC = (AsigMC > 0.0) ? CsigMC / AsigMC : 0.0;
    const double fD = (AsigMC > 0.0) ? DsigMC / AsigMC : 0.0;

    double Scorr = 0.0;
    const bool solverOk = ARJ::SolveLeakageCorrectedSA(A, B, C, D, fB, fC, fD, Scorr);
    const double Pcorr = (solverOk && A > 0.0) ? Scorr / A : Praw;
    const double PcorrErr = LeakageCorrectedPurityError(A, B, C, D, eA, eB, eC, eD, fB, fC, fD);

    const double ptLo = ARJ::kPtEdges[(std::size_t)i];
    const double ptHi = ARJ::kPtEdges[(std::size_t)i + 1];
    out << i << " " << ptLo << " " << ptHi << " " << b.suffix << " "
        << A << " " << eA << " " << B << " " << eB << " "
        << C << " " << eC << " " << D << " " << eD << " "
        << rawSub << " " << Sraw << " " << Praw << " " << PrawErr << " "
        << AsigMC << " " << BsigMC << " " << CsigMC << " " << DsigMC << " "
        << fB << " " << fC << " " << fD << " "
        << (solverOk ? 1 : 0) << " " << Scorr << " " << Pcorr << " " << PcorrErr << " "
        << (hLeak ? 1 : 0) << "\n";

    if (i < 4 || i == ARJ::kNPtBins - 1)
    {
      std::cout << "  " << s.label << " " << b.suffix
                << " DATA(A,B,C,D)=(" << A << ", " << B << ", " << C << ", " << D << ")"
                << " Praw=" << Praw
                << " f=(" << fB << ", " << fC << ", " << fD << ")"
                << " ok=" << solverOk
                << " Pcorr=" << Pcorr << "\n";
    }
  }
}

void DrawTopLabels(const char* title)
{
  TLatex tx;
  tx.SetNDC(true);
  tx.SetTextFont(42);
  tx.SetTextAlign(22);
  tx.SetTextSize(0.032);
  tx.DrawLatex(0.50, 0.94, title);

  tx.SetTextAlign(31);
  tx.SetTextSize(0.034);
  tx.DrawLatex(0.92, 0.86, "#bf{sPHENIX} #it{Internal}");
  tx.DrawLatex(0.92, 0.82, "Pythia overlay, #sqrt{s_{NN}} = 200 GeV");
}

void DrawXJOverlay(TDirectory* ppDir,
                   TDirectory* aaDir,
                   const double ptMin,
                   const double ptMax,
                   const char* outSuffix)
{
  const std::vector<Series> series = {
      {"0-20% Au+Au", "_cent_0_20", kBlack, 20, false},
      {"50-80% Au+Au", "_cent_50_80", kBlue + 1, 20, false},
      {"pp", "", kRed + 1, 24, true},
  };

  std::vector<std::unique_ptr<TH1D>> hists;
  hists.push_back(ProjectXJ(aaDir, series[0], ptMin, ptMax, true));
  hists.push_back(ProjectXJ(aaDir, series[1], ptMin, ptMax, true));
  hists.push_back(ProjectXJ(ppDir, series[2], ptMin, ptMax, true));

  double ymax = 0.0;
  for (const auto& h : hists)
  {
    if (h) ymax = std::max(ymax, h->GetMaximum());
  }
  if (ymax <= 0.0) ymax = 1.0;

  TCanvas c("c_xJ_reference_pp_AuAu_overlay", "c_xJ_reference_pp_AuAu_overlay", 1000, 760);
  c.SetLeftMargin(0.13);
  c.SetRightMargin(0.04);
  c.SetTopMargin(0.11);
  c.SetBottomMargin(0.13);
  c.SetTicks(1, 1);

  TH1F frame("hFrameXJRef", "", 100, 0.0, 2.0);
  frame.SetDirectory(nullptr);
  frame.SetStats(0);
  frame.SetMinimum(0.0);
  frame.SetMaximum(1.45 * ymax);
  frame.GetXaxis()->SetTitle("x_{J#gamma}");
  frame.GetYaxis()->SetTitle("Normalized to unit area");
  frame.GetXaxis()->SetTitleSize(0.048);
  frame.GetYaxis()->SetTitleSize(0.048);
  frame.GetXaxis()->SetLabelSize(0.040);
  frame.GetYaxis()->SetLabelSize(0.040);
  frame.Draw();

  for (auto& h : hists)
  {
    if (h) h->Draw("E1 SAME");
  }

  TLegend leg(0.58, 0.54, 0.88, 0.72);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextFont(42);
  leg.SetTextSize(0.033);
  for (std::size_t i = 0; i < hists.size(); ++i)
  {
    if (hists[i]) leg.AddEntry(hists[i].get(), series[i].label.c_str(), "pe");
  }
  leg.Draw();

  TLatex tx;
  tx.SetNDC(true);
  tx.SetTextFont(42);
  tx.SetTextAlign(13);
  tx.SetTextSize(0.028);
  tx.DrawLatex(0.16, 0.87,
               TString::Format("Inclusive x_{J#gamma}, R = 0.3, p_{T}^{#gamma} = %.0f-%.0f GeV", ptMin, ptMax).Data());
  tx.DrawLatex(0.16, 0.83, "Reference preselection, reference tight/non-tight");
  tx.DrawLatex(0.16, 0.79, "Au+Au: Run25 DATA, sliding iso");
  tx.DrawLatex(0.16, 0.75, "pp: Run24 DATA, fixed E_{T}^{iso} < 2 GeV");

  tx.SetTextAlign(31);
  tx.SetTextSize(0.034);
  tx.DrawLatex(0.92, 0.82, "#bf{sPHENIX} #it{Internal}");
  tx.DrawLatex(0.92, 0.78, "#sqrt{s_{NN}} = 200 GeV");

  const std::string out = OutputDir() + "/xJ_recoInclusive_reference_pp_AuAu_" + std::string(outSuffix) + "_overlay.png";
  c.SaveAs(out.c_str());
  std::cout << "[WROTE] " << out << "\n";
}

void DrawPurityOverlay(TDirectory* ppDataDir,
                       TDirectory* aaDataDir,
                       TDirectory* ppLeakageDir,
                       TDirectory* aaLeakageDir,
                       const bool leakageCorrected)
{
  const std::vector<Series> series = {
      {"pp", "", kRed + 1, 24, true},
      {"Au+Au 0-20%", "_cent_0_20", kBlack, 20, false},
      {"Au+Au 20-50%", "_cent_20_50", kBlue + 1, 20, false},
      {"Au+Au 50-80%", "_cent_50_80", kOrange + 7, 20, false},
  };

  std::vector<std::unique_ptr<TGraphErrors>> graphs;
  graphs.push_back(BuildPurityGraph(ppDataDir, ppLeakageDir, series[0], leakageCorrected));
  graphs.push_back(BuildPurityGraph(aaDataDir, aaLeakageDir, series[1], leakageCorrected));
  graphs.push_back(BuildPurityGraph(aaDataDir, aaLeakageDir, series[2], leakageCorrected));
  graphs.push_back(BuildPurityGraph(aaDataDir, aaLeakageDir, series[3], leakageCorrected));

  TCanvas c("c_purity_reference_pp_AuAu_overlay", "c_purity_reference_pp_AuAu_overlay", 1000, 760);
  c.SetLeftMargin(0.13);
  c.SetRightMargin(0.04);
  c.SetTopMargin(0.11);
  c.SetBottomMargin(0.13);
  c.SetTicks(1, 1);

  TH1F frame(leakageCorrected ? "hFramePurityCorrRef" : "hFramePurityRef", "", 100, 5.0, 35.0);
  frame.SetDirectory(nullptr);
  frame.SetStats(0);
  frame.SetMinimum(0.0);
  frame.SetMaximum(1.50);
  frame.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
  frame.GetYaxis()->SetTitle(leakageCorrected ? "Purity (leakage corrected)" : "Purity (raw ABCD)");
  frame.GetXaxis()->SetTitleSize(0.048);
  frame.GetYaxis()->SetTitleSize(0.048);
  frame.GetXaxis()->SetLabelSize(0.040);
  frame.GetYaxis()->SetLabelSize(0.040);
  frame.Draw();

  for (auto& g : graphs)
  {
    if (g) g->Draw("P SAME");
  }

  TLine one(5.0, 1.0, 35.0, 1.0);
  one.SetLineStyle(2);
  one.SetLineColor(kGray + 2);
  one.Draw("SAME");

  TLegend leg(0.58, 0.66, 0.91, 0.79);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextFont(42);
  leg.SetTextSize(0.032);
  leg.SetNColumns(2);
  for (std::size_t i = 0; i < graphs.size(); ++i)
  {
    if (graphs[i]) leg.AddEntry(graphs[i].get(), series[i].label.c_str(), "pe");
  }
  leg.Draw();

  TLatex tx;
  tx.SetNDC(true);
  tx.SetTextFont(42);
  tx.SetTextAlign(22);
  tx.SetTextSize(0.032);
  tx.DrawLatex(0.50, 0.94,
               leakageCorrected
                   ? "ABCD leakage-corrected purity, DATA with SIM leakage factors"
                   : "ABCD raw purity, DATA pp vs Au+Au centrality");

  tx.SetTextAlign(13);
  tx.SetTextSize(0.027);
  tx.DrawLatex(0.16, 0.86, "Reference preselection, reference tight/non-tight");
  tx.DrawLatex(0.16, 0.82, "#DeltaR_{cone} < 0.4");
  tx.DrawLatex(0.16, 0.78, "Observed ABCD: Run25 Au+Au DATA and Run24 pp DATA");
  tx.DrawLatex(0.16, 0.74, "Leakage factors: stitched photon+jet SIM");

  tx.SetTextAlign(31);
  tx.SetTextSize(0.034);
  tx.DrawLatex(0.92, 0.86, "#bf{sPHENIX} #it{Internal}");
  tx.DrawLatex(0.92, 0.82, "#sqrt{s_{NN}} = 200 GeV");

  const std::string out = OutputDir() + (leakageCorrected
      ? "/purity_leakageCorrected_reference_pp_AuAuCentrality_overlay.png"
      : "/purity_rawABCD_reference_pp_AuAuCentrality_overlay.png");
  c.SaveAs(out.c_str());
  std::cout << "[WROTE] " << out << "\n";
}
}

void PlotPPvsAuAuXJPurityReferenceOverlay()
{
  gStyle->SetOptStat(0);
  gSystem->mkdir(OutputDir().c_str(), true);

  if (!RebuildInputs())
  {
    std::cerr << "[ERROR] Failed to rebuild one or both merged SIM inputs.\n";
    return;
  }

  std::unique_ptr<TFile> fPP(TFile::Open(PPMergedPath().c_str(), "READ"));
  std::unique_ptr<TFile> fAA(TFile::Open(AAMergedPath().c_str(), "READ"));
  std::unique_ptr<TFile> fPPData(TFile::Open(PPDataPath().c_str(), "READ"));
  std::unique_ptr<TFile> fAAData(TFile::Open(AADataPath().c_str(), "READ"));
  if (!fPP || fPP->IsZombie() || !fAA || fAA->IsZombie())
  {
    std::cerr << "[ERROR] Cannot open merged pp/AuAu inputs.\n";
    return;
  }
  if (!fPPData || fPPData->IsZombie() || !fAAData || fAAData->IsZombie())
  {
    std::cerr << "[ERROR] Cannot open pp/AuAu DATA inputs:\n"
              << "  pp   = " << PPDataPath() << "\n"
              << "  AuAu = " << AADataPath() << "\n";
    return;
  }

  TDirectory* ppDir = fPP->GetDirectory(ARJ::kDirSIM.c_str());
  TDirectory* aaDir = fAA->GetDirectory(ARJ::kDirSIM.c_str());
  TDirectory* ppDataDir = fPPData->GetDirectory(ARJ::kTriggerPP.c_str());
  TDirectory* aaDataDir = fAAData->GetDirectory("photon_10_plus_MBD_NS_geq_2_vtx_lt_150");
  if (!ppDir || !aaDir)
  {
    std::cerr << "[ERROR] Missing SIM directory in merged pp/AuAu inputs.\n";
    return;
  }
  if (!ppDataDir || !aaDataDir)
  {
    std::cerr << "[ERROR] Missing DATA trigger directories for purity overlay.\n";
    return;
  }

  const std::vector<Series> auditSeries = {
      {"pp", "", kRed + 1, 24, true},
      {"Au+Au 0-20%", "_cent_0_20", kBlack, 20, false},
      {"Au+Au 20-50%", "_cent_20_50", kBlue + 1, 20, false},
      {"Au+Au 50-80%", "_cent_50_80", kOrange + 7, 20, false},
  };
  AuditPuritySeries(ppDataDir, ppDir, auditSeries[0],
                    PPDataPath(), ARJ::kTriggerPP,
                    PPMergedPath(), ARJ::kDirSIM);
  for (std::size_t i = 1; i < auditSeries.size(); ++i)
  {
    AuditPuritySeries(aaDataDir, aaDir, auditSeries[i],
                      AADataPath(), "photon_10_plus_MBD_NS_geq_2_vtx_lt_150",
                      AAMergedPath(), ARJ::kDirSIM);
  }

  DrawXJOverlay(ppDataDir, aaDataDir, 16.0, 40.0, "16_40");
  DrawXJOverlay(ppDataDir, aaDataDir, 20.0, 35.0, "20_35");
  DrawPurityOverlay(ppDataDir, aaDataDir, ppDir, aaDir, false);
  DrawPurityOverlay(ppDataDir, aaDataDir, ppDir, aaDir, true);
}
