#include <TBox.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace
{
const std::string kBase = "/Users/patsfan753/Desktop/ThesisAnalysis/";
const std::string kOutDir = kBase + "dataOutput/stitchDiagnostics/focus21_photon_variants_20260521";
const std::string kCsv = kOutDir + "/photon_stitch_variants_bins.csv";

struct Point
{
  double x = 0.0;
  double ex = 0.0;
  double y = 0.0;
  double ey = 0.0;
};

std::vector<std::string> Split(const std::string& s, char delim)
{
  std::vector<std::string> out;
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, delim)) out.push_back(item);
  return out;
}

std::vector<Point> ReadVariant(const std::string& variant)
{
  std::ifstream in(kCsv);
  std::vector<Point> pts;
  std::string line;
  while (std::getline(in, line))
  {
    if (line.empty() || line[0] == '#') continue;
    if (line.rfind("variant,", 0) == 0) continue;
    const auto tok = Split(line, ',');
    if (tok.size() != 7 || tok[0] != variant) continue;
    const double lo = std::stod(tok[2]);
    const double hi = std::stod(tok[3]);
    Point p;
    p.x = std::stod(tok[4]);
    p.ex = 0.5 * (hi - lo);
    p.y = std::stod(tok[5]);
    p.ey = std::stod(tok[6]);
    if (p.x >= 10.0 && p.x <= 35.0 && p.y > 0.0) pts.push_back(p);
  }
  return pts;
}

std::vector<Point> SelectRange(const std::vector<Point>& pts, double lo, double hi)
{
  std::vector<Point> out;
  for (const auto& p : pts)
  {
    if (p.x >= lo && p.x < hi) out.push_back(p);
  }
  return out;
}

TGraphErrors* MakeGraph(const std::vector<Point>& pts, const char* name, int color, int marker)
{
  auto* g = new TGraphErrors();
  g->SetName(name);
  for (size_t i = 0; i < pts.size(); ++i)
  {
    g->SetPoint(i, pts[i].x, pts[i].y);
    g->SetPointError(i, pts[i].ex, pts[i].ey);
  }
  g->SetLineColor(color);
  g->SetMarkerColor(color);
  g->SetMarkerStyle(marker);
  g->SetMarkerSize(0.86);
  g->SetLineWidth(1);
  return g;
}

TGraphErrors* MakeLogGraph(const std::vector<Point>& pts)
{
  auto* g = new TGraphErrors();
  int ip = 0;
  for (const auto& p : pts)
  {
    if (p.y <= 0.0) continue;
    g->SetPoint(ip, p.x, std::log(p.y));
    g->SetPointError(ip, p.ex, p.ey > 0.0 ? p.ey / p.y : 0.0);
    ++ip;
  }
  return g;
}

double RefEval(const TF1* logFit, double x)
{
  const double logy = logFit->Eval(x);
  return std::isfinite(logy) ? std::exp(logy) : 0.0;
}

TGraph* MakeReference(const TF1* logFit, int color = kBlack)
{
  auto* g = new TGraph();
  int ip = 0;
  for (double x = 12.0; x <= 35.0001; x += 0.05)
  {
    g->SetPoint(ip, x, RefEval(logFit, x));
    ++ip;
  }
  g->SetLineColor(color);
  g->SetLineStyle(2);
  g->SetLineWidth(3);
  return g;
}

TGraphErrors* MakeRatio(const std::vector<Point>& pts, const TF1* logFit, const char* name, int color, int marker)
{
  auto* g = new TGraphErrors();
  g->SetName(name);
  int ip = 0;
  for (const auto& p : pts)
  {
    const double fit = RefEval(logFit, p.x);
    if (fit <= 0.0) continue;
    g->SetPoint(ip, p.x, p.y / fit);
    g->SetPointError(ip, p.ex, p.ey / fit);
    ++ip;
  }
  g->SetLineColor(color);
  g->SetMarkerColor(color);
  g->SetMarkerStyle(marker);
  g->SetMarkerSize(0.86);
  g->SetLineWidth(1);
  return g;
}

double MeanRatio(const std::vector<Point>& pts, const TF1* logFit, double lo, double hi)
{
  double sum = 0.0;
  int n = 0;
  for (const auto& p : pts)
  {
    if (p.x < lo || p.x >= hi) continue;
    const double fit = RefEval(logFit, p.x);
    if (fit <= 0.0) continue;
    sum += p.y / fit;
    ++n;
  }
  return n > 0 ? sum / n : std::numeric_limits<double>::quiet_NaN();
}

void SetupPad(TPad& top, TPad& bot)
{
  top.SetLeftMargin(0.095);
  top.SetRightMargin(0.035);
  top.SetTopMargin(0.20);
  top.SetBottomMargin(0.02);
  top.SetLogy(true);
  bot.SetLeftMargin(0.095);
  bot.SetRightMargin(0.035);
  bot.SetTopMargin(0.035);
  bot.SetBottomMargin(0.28);
}

void DrawBoundary(double ylo, double yhi)
{
  TLine* line = new TLine(21.0, ylo, 21.0, yhi);
  line->SetLineColor(kGray + 2);
  line->SetLineStyle(3);
  line->SetLineWidth(3);
  line->Draw("SAME");
}

void DrawTopLabel(const char* title, const char* subtitle)
{
  TLatex lat;
  lat.SetNDC(true);
  lat.SetTextFont(42);
  lat.SetTextSize(0.038);
  lat.DrawLatex(0.115, 0.955, title);
  lat.SetTextSize(0.026);
  lat.DrawLatex(0.115, 0.905, subtitle);
  lat.DrawLatex(0.115, 0.855, "x-axis capped at 35 GeV; reference is one common modified-power-law fit to the bounded candidate");
  lat.SetTextAlign(33);
  lat.SetTextSize(0.035);
  lat.DrawLatex(0.955, 0.955, "#it{#bf{sPHENIX}} Internal");
  lat.SetTextSize(0.027);
  lat.DrawLatex(0.955, 0.905, "PYTHIA8 embedded photon+jet");
  lat.SetTextAlign(11);
}

void DrawTakeaway(double x1, double y1, double x2, double y2, const char* line1,
                  const char* line2, int fillColor, int textColor)
{
  TLatex lat;
  lat.SetNDC(true);
  lat.SetTextFont(42);
  lat.SetTextColor(textColor);
  lat.SetTextSize(0.064);
  lat.DrawLatex(x1 + 0.02, y2 - 0.09, line1);
  lat.SetTextSize(0.046);
  lat.DrawLatex(x1 + 0.02, y2 - 0.18, line2);
}

void DrawLowerOnlyTable()
{
  TLatex lat;
  lat.SetNDC(true);
  lat.SetTextFont(42);
  lat.SetTextColor(kBlack);
  lat.SetTextSize(0.030);
  lat.DrawLatex(0.135, 0.402, "No-upper-bound diagnostic");
  lat.SetTextSize(0.024);
  lat.DrawLatex(0.135, 0.360, "sample");
  lat.DrawLatex(0.250, 0.360, "gate used");
  lat.DrawLatex(0.370, 0.360, "#sigma_{eff} pb / scale");

  lat.SetTextSize(0.025);
  lat.SetTextColor(kBlue + 1);
  lat.DrawLatex(0.135, 0.322, "PhotonJet12");
  lat.DrawLatex(0.250, 0.322, "p_{T}^{#gamma,filter} #geq 12");
  lat.DrawLatex(0.370, 0.322, "2741.7 / 21.66");

  lat.SetTextColor(kOrange + 7);
  lat.DrawLatex(0.135, 0.286, "PhotonJet20");
  lat.DrawLatex(0.250, 0.286, "p_{T}^{#gamma,filter} #geq 21");
  lat.DrawLatex(0.370, 0.286, "89.96 / 1.00");

  lat.SetTextColor(kBlack);
  lat.SetTextSize(0.030);
  lat.DrawLatex(0.135, 0.238, "Bounded reference used for gray points and fit");
  lat.SetTextSize(0.024);
  lat.DrawLatex(0.135, 0.200, "sample");
  lat.DrawLatex(0.250, 0.200, "gate used");
  lat.DrawLatex(0.370, 0.200, "#sigma_{eff} pb / scale");
  lat.SetTextColor(kBlue + 1);
  lat.DrawLatex(0.135, 0.164, "PhotonJet12");
  lat.DrawLatex(0.250, 0.164, "12 #leq p_{T}^{#gamma,filter} < 21");
  lat.DrawLatex(0.370, 0.164, "2663.5 / 21.27");
  lat.SetTextColor(kRed + 1);
  lat.DrawLatex(0.135, 0.128, "PhotonJet20");
  lat.DrawLatex(0.250, 0.128, "p_{T}^{#gamma,filter} #geq 21");
  lat.DrawLatex(0.370, 0.128, "92.05 / 1.00");
}

void DrawRatioFitNote()
{
  auto* box = new TPaveText(0.125, 0.735, 0.555, 0.925, "NDC");
  box->SetFillColorAlpha(kWhite, 0.88);
  box->SetLineColor(kGray + 1);
  box->SetLineWidth(2);
  box->Draw("SAME");

  TLatex lat;
  lat.SetNDC(true);
  lat.SetTextFont(42);
  lat.SetTextSize(0.064);
  lat.SetTextColor(kBlack);
  lat.DrawLatex(0.145, 0.855, "Ratio denominator: bounded fit");
  lat.SetTextSize(0.046);
  lat.DrawLatex(0.145, 0.785, "Fit form: modified power law A*(1/x)^(b+c ln x+d x)");
  lat.DrawLatex(0.145, 0.728, "orange near 2 means both samples are counted");
}

void DrawCandidate(const std::vector<Point>& bounded, TF1* logFit)
{
  const auto bounded12 = SelectRange(bounded, 10.0, 21.0);
  const auto bounded20 = SelectRange(bounded, 21.0, 35.1);
  auto* gBounded12 = MakeGraph(bounded12, "gBoundedCandidatePhoton12", kBlue + 1, 20);
  auto* gBounded20 = MakeGraph(bounded20, "gBoundedCandidatePhoton20", kRed + 1, 21);
  auto* gReference = MakeReference(logFit);
  auto* rBounded12 = MakeRatio(bounded12, logFit, "rBoundedCandidatePhoton12", kBlue + 1, 20);
  auto* rBounded20 = MakeRatio(bounded20, logFit, "rBoundedCandidatePhoton20", kRed + 1, 21);

  const double below = MeanRatio(bounded, logFit, 19.0, 21.0);
  const double above = MeanRatio(bounded, logFit, 21.0, 23.0);
  const double jump = above / below;

  TCanvas c("c_blair_bounded21", "bounded 21 GeV PhotonJet stitch", 1600, 900);
  TPad top("top_bounded", "top_bounded", 0.0, 0.32, 1.0, 1.0);
  TPad bot("bot_bounded", "bot_bounded", 0.0, 0.0, 1.0, 0.32);
  SetupPad(top, bot);
  top.Draw();
  bot.Draw();

  top.cd();
  TH1F frame("frame_bounded", "", 100, 10.0, 35.0);
  frame.SetMinimum(2.0e3);
  frame.SetMaximum(1.1e8);
  frame.GetXaxis()->SetLabelSize(0.0);
  frame.GetYaxis()->SetTitle("Weighted stitched entries / bin");
  frame.GetYaxis()->SetTitleSize(0.047);
  frame.GetYaxis()->SetLabelSize(0.038);
  frame.GetYaxis()->SetTitleOffset(0.92);
  frame.Draw("AXIS");
  gReference->Draw("L SAME");
  gBounded12->Draw("P SAME");
  gBounded20->Draw("P SAME");
  DrawBoundary(frame.GetMinimum(), frame.GetMaximum());
  DrawTopLabel("PhotonJet12+20 bounded 21 GeV ownership is smooth",
               "PhotonJet12: 12 #leq p_{T}^{#gamma,filter} < 21   |   PhotonJet20: p_{T}^{#gamma,filter} #geq 21");

  TLegend leg(0.62, 0.63, 0.94, 0.78);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextFont(42);
  leg.SetTextSize(0.031);
  leg.AddEntry(gBounded12, "PhotonJet12-owned: 12-21 GeV", "lep");
  leg.AddEntry(gBounded20, "PhotonJet20-owned: #geq21 GeV", "lep");
  leg.AddEntry(gReference, "Common modified-power-law fit", "l");
  leg.Draw();

  bot.cd();
  TH1F rframe("rframe_bounded", "", 100, 10.0, 35.0);
  rframe.SetMinimum(0.94);
  rframe.SetMaximum(1.06);
  rframe.GetXaxis()->SetTitle("Truth photon filter p_{T} [GeV]");
  rframe.GetYaxis()->SetTitle("stitched / fit");
  rframe.GetXaxis()->SetTitleSize(0.084);
  rframe.GetXaxis()->SetLabelSize(0.070);
  rframe.GetYaxis()->SetTitleSize(0.075);
  rframe.GetYaxis()->SetLabelSize(0.060);
  rframe.GetYaxis()->SetTitleOffset(0.55);
  rframe.GetYaxis()->SetNdivisions(505);
  rframe.Draw("AXIS");
  TLine one(10.0, 1.0, 35.0, 1.0);
  one.SetLineColor(kGray + 2);
  one.SetLineStyle(2);
  one.Draw("SAME");
  DrawBoundary(0.94, 1.06);
  rBounded12->Draw("P SAME");
  rBounded20->Draw("P SAME");
  DrawTakeaway(0.50, 0.62, 0.94, 0.92, Form("Smooth handoff: jump = %.3f", jump),
               Form("<R>19-21 = %.3f,  <R>21-23 = %.3f", below, above), kBlue + 1, kBlue + 2);

  c.SaveAs((kOutDir + "/photon12plus20_bounded21_blair_takeaway.png").c_str());
}

void DrawLowerOnlyDiagnostic(const std::vector<Point>& bounded, const std::vector<Point>& lower, TF1* logFit)
{
  auto* gBounded = MakeGraph(bounded, "gBoundedReference", kGray + 2, 24);
  const auto lowerBelow = SelectRange(lower, 10.0, 21.0);
  const auto lowerOverlap = SelectRange(lower, 21.0, 35.1);
  auto* gLowerBelow = MakeGraph(lowerBelow, "gLowerOnlyBelowBoundary", kBlue + 1, 20);
  auto* gLowerOverlap = MakeGraph(lowerOverlap, "gLowerOnlyOverlap", kOrange + 7, 21);
  auto* gReference = MakeReference(logFit);
  auto* rBounded = MakeRatio(bounded, logFit, "rBoundedReference", kGray + 2, 24);
  auto* rLowerBelow = MakeRatio(lowerBelow, logFit, "rLowerOnlyBelowBoundary", kBlue + 1, 20);
  auto* rLowerOverlap = MakeRatio(lowerOverlap, logFit, "rLowerOnlyOverlap", kOrange + 7, 21);

  const double below = MeanRatio(lower, logFit, 19.0, 21.0);
  const double above = MeanRatio(lower, logFit, 21.0, 23.0);
  const double jump = above / below;

  TCanvas c("c_blair_loweronly", "lower-only PhotonJet stitch diagnostic", 1600, 900);
  TPad top("top_lower", "top_lower", 0.0, 0.32, 1.0, 1.0);
  TPad bot("bot_lower", "bot_lower", 0.0, 0.0, 1.0, 0.32);
  SetupPad(top, bot);
  top.Draw();
  bot.Draw();

  top.cd();
  TH1F frame("frame_lower", "", 100, 10.0, 35.0);
  frame.SetMinimum(2.0e3);
  frame.SetMaximum(1.1e8);
  frame.GetXaxis()->SetLabelSize(0.0);
  frame.GetYaxis()->SetTitle("Weighted stitched entries / bin");
  frame.GetYaxis()->SetTitleSize(0.047);
  frame.GetYaxis()->SetLabelSize(0.038);
  frame.GetYaxis()->SetTitleOffset(0.92);
  frame.Draw("AXIS");
  gReference->Draw("L SAME");
  gBounded->Draw("P SAME");
  gLowerBelow->Draw("P SAME");
  gLowerOverlap->Draw("P SAME");
  DrawBoundary(frame.GetMinimum(), frame.GetMaximum());
  DrawTopLabel("No-upper-bound test fails: samples overlap above 21 GeV",
               "Lower-only gates let PhotonJet12 and PhotonJet20 both contribute for p_{T}^{#gamma,filter} #geq 21");

  TLegend leg(0.58, 0.60, 0.94, 0.79);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextFont(42);
  leg.SetTextSize(0.030);
  leg.AddEntry(gLowerBelow, "Below 21: PhotonJet12 region", "lep");
  leg.AddEntry(gLowerOverlap, "Above 21: overlapping lower-only sum", "lep");
  leg.AddEntry(gBounded, "Bounded candidate reference", "lep");
  leg.AddEntry(gReference, "Bounded fit: modified power law A(1/x)^{b+c ln x+d x}", "l");
  leg.Draw();
  DrawLowerOnlyTable();

  bot.cd();
  TH1F rframe("rframe_lower", "", 100, 10.0, 35.0);
  rframe.SetMinimum(0.82);
  rframe.SetMaximum(2.20);
  rframe.GetXaxis()->SetTitle("Truth photon filter p_{T} [GeV]");
  rframe.GetYaxis()->SetTitle("stitched / fit");
  rframe.GetXaxis()->SetTitleSize(0.084);
  rframe.GetXaxis()->SetLabelSize(0.070);
  rframe.GetYaxis()->SetTitleSize(0.075);
  rframe.GetYaxis()->SetLabelSize(0.060);
  rframe.GetYaxis()->SetTitleOffset(0.55);
  rframe.GetYaxis()->SetNdivisions(505);
  rframe.Draw("AXIS");
  TLine one(10.0, 1.0, 35.0, 1.0);
  one.SetLineColor(kGray + 2);
  one.SetLineStyle(2);
  one.Draw("SAME");
  DrawBoundary(0.82, 2.20);
  rBounded->Draw("P SAME");
  rLowerBelow->Draw("P SAME");
  rLowerOverlap->Draw("P SAME");
  c.SaveAs((kOutDir + "/photon12plus20_loweronly_overlap_blair_takeaway.png").c_str());
}
}

void MakeFocus21PhotonStitchBlairPNGs()
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetErrorX(0.5);

  const auto bounded = ReadVariant("bounded_12_21_ge21");
  const auto lower = ReadVariant("loweronly_ge12_ge21");
  if (bounded.empty() || lower.empty())
  {
    std::cerr << "[ERROR] Missing input points from " << kCsv << std::endl;
    return;
  }

  auto* gLogBounded = MakeLogGraph(bounded);
  TF1 logFit("log_fit_ppg12_modified_power_blair",
             "[0] + ([1] + [2]*TMath::Log(x) + [3]*x)*TMath::Log(1.0/x)",
             12.0, 35.0);
  logFit.SetParameters(25.0, 2.0, 0.0, 0.0);
  gLogBounded->Fit(&logFit, "QNR");

  gSystem->mkdir(kOutDir.c_str(), true);
  DrawCandidate(bounded, &logFit);
  DrawLowerOnlyDiagnostic(bounded, lower, &logFit);

  std::cout << "[DONE] " << kOutDir << "/photon12plus20_bounded21_blair_takeaway.png" << std::endl;
  std::cout << "[DONE] " << kOutDir << "/photon12plus20_loweronly_overlap_blair_takeaway.png" << std::endl;
}
