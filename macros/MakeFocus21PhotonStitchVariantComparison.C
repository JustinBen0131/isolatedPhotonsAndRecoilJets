#include <TCanvas.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
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
    Point p;
    const double lo = std::stod(tok[2]);
    const double hi = std::stod(tok[3]);
    p.x = std::stod(tok[4]);
    p.ex = 0.5 * (hi - lo);
    p.y = std::stod(tok[5]);
    p.ey = std::stod(tok[6]);
    if (p.x >= 10.0 && p.x <= 35.0 && p.y > 0.0) pts.push_back(p);
  }
  return pts;
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
  g->SetMarkerSize(0.72);
  g->SetLineWidth(1);
  return g;
}

TGraphErrors* MakeLogGraph(const std::vector<Point>& pts, const char* name)
{
  auto* g = new TGraphErrors();
  g->SetName(name);
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

TGraph* MakeReferenceGraph(const TF1* logFit)
{
  auto* g = new TGraph();
  g->SetName("gReferenceFit");
  int ip = 0;
  for (double x = 12.0; x <= 35.0001; x += 0.05)
  {
    g->SetPoint(ip, x, RefEval(logFit, x));
    ++ip;
  }
  g->SetLineColor(kBlack);
  g->SetLineStyle(2);
  g->SetLineWidth(2);
  return g;
}

TGraphErrors* MakeRatioGraph(const std::vector<Point>& pts, TF1* logFit, const char* name, int color, int marker)
{
  auto* g = new TGraphErrors();
  g->SetName(name);
  int ip = 0;
  for (const auto& p : pts)
  {
    const double fy = RefEval(logFit, p.x);
    if (fy <= 0.0) continue;
    g->SetPoint(ip, p.x, p.y / fy);
    g->SetPointError(ip, p.ex, p.ey / fy);
    ++ip;
  }
  g->SetLineColor(color);
  g->SetMarkerColor(color);
  g->SetMarkerStyle(marker);
  g->SetMarkerSize(0.72);
  g->SetLineWidth(1);
  return g;
}

double MeanRatio(const std::vector<Point>& pts, TF1* logFit, double lo, double hi)
{
  double sum = 0.0;
  int n = 0;
  for (const auto& p : pts)
  {
    if (p.x < lo || p.x >= hi) continue;
    const double fy = RefEval(logFit, p.x);
    if (fy <= 0.0) continue;
    sum += p.y / fy;
    ++n;
  }
  return n > 0 ? sum / n : std::numeric_limits<double>::quiet_NaN();
}
}

void MakeFocus21PhotonStitchVariantComparison()
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetErrorX(0.5);

  auto bounded = ReadVariant("bounded_12_21_ge21");
  auto lower = ReadVariant("loweronly_ge12_ge21");
  if (bounded.empty() || lower.empty())
  {
    std::cerr << "[ERROR] Missing input points from " << kCsv << std::endl;
    return;
  }

  auto* gBounded = MakeGraph(bounded, "gBounded", kBlue + 1, 20);
  auto* gLower = MakeGraph(lower, "gLower", kOrange + 7, 21);

  auto* gLogBounded = MakeLogGraph(bounded, "gLogBounded");
  TF1 logFit("log_fit_ppg12_modified_power",
             "[0] + ([1] + [2]*TMath::Log(x) + [3]*x)*TMath::Log(1.0/x)",
             12.0, 35.0);
  logFit.SetParameters(25.0, 2.0, 0.0, 0.0);
  gLogBounded->Fit(&logFit, "QNR");
  auto* gReference = MakeReferenceGraph(&logFit);

  auto* rBounded = MakeRatioGraph(bounded, &logFit, "rBounded", kBlue + 1, 20);
  auto* rLower = MakeRatioGraph(lower, &logFit, "rLower", kOrange + 7, 21);

  const double belowBounded = MeanRatio(bounded, &logFit, 19.0, 21.0);
  const double aboveBounded = MeanRatio(bounded, &logFit, 21.0, 23.0);
  const double belowLower = MeanRatio(lower, &logFit, 19.0, 21.0);
  const double aboveLower = MeanRatio(lower, &logFit, 21.0, 23.0);

  gSystem->mkdir(kOutDir.c_str(), true);

  TCanvas c("c_focus21_photon_stitch", "focus21 photon stitch comparison", 1500, 930);
  TPad top("top", "top", 0.0, 0.32, 1.0, 1.0);
  TPad bot("bot", "bot", 0.0, 0.0, 1.0, 0.32);
  top.SetLeftMargin(0.105);
  top.SetRightMargin(0.035);
  top.SetTopMargin(0.075);
  top.SetBottomMargin(0.02);
  top.SetLogy(true);
  bot.SetLeftMargin(0.105);
  bot.SetRightMargin(0.035);
  bot.SetTopMargin(0.03);
  bot.SetBottomMargin(0.28);
  top.Draw();
  bot.Draw();

  top.cd();
  TH1F frame("frame", "", 100, 10.0, 35.0);
  frame.SetMinimum(2.0e3);
  frame.SetMaximum(1.2e8);
  frame.GetXaxis()->SetLabelSize(0.0);
  frame.GetYaxis()->SetTitle("Weighted stitched entries / bin");
  frame.GetYaxis()->SetTitleSize(0.047);
  frame.GetYaxis()->SetLabelSize(0.038);
  frame.GetYaxis()->SetTitleOffset(1.04);
  frame.Draw("AXIS");
  gBounded->Draw("P SAME");
  gLower->Draw("P SAME");
  gReference->Draw("L SAME");

  TLine boundary(21.0, frame.GetMinimum(), 21.0, frame.GetMaximum());
  boundary.SetLineColor(kGray + 2);
  boundary.SetLineStyle(3);
  boundary.SetLineWidth(2);
  boundary.Draw("SAME");

  TLatex lat;
  lat.SetNDC(true);
  lat.SetTextFont(42);
  lat.SetTextSize(0.039);
  lat.DrawLatex(0.125, 0.90, "Embedded PhotonJet12+20 stitching variants");
  lat.SetTextSize(0.028);
  lat.DrawLatex(0.125, 0.845, "Final stitched RecoilJets outputs; common modified-power-law reference fit to bounded 21 GeV candidate");
  lat.DrawLatex(0.125, 0.795, "Boundary test: move PhotonJet12/20 handoff from 20 to 21 GeV");
  lat.SetTextAlign(33);
  lat.SetTextSize(0.038);
  lat.DrawLatex(0.955, 0.90, "#it{#bf{sPHENIX}} Internal");
  lat.SetTextSize(0.029);
  lat.DrawLatex(0.955, 0.85, "PYTHIA8 embedded photon+jet");
  lat.SetTextAlign(11);

  TLegend leg(0.50, 0.18, 0.94, 0.40);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextFont(42);
  leg.SetTextSize(0.029);
  leg.AddEntry(gBounded, "Bounded: PhotonJet12 12-21, PhotonJet20 #geq21", "lep");
  leg.AddEntry(gLower, "Diagnostic lower-only: PhotonJet12 #geq12, PhotonJet20 #geq21", "lep");
  leg.AddEntry(gReference, "Fit: A(1/x)^{b + c ln(x) + d x}", "l");
  leg.Draw();

  bot.cd();
  TH1F rframe("rframe", "", 100, 10.0, 35.0);
  rframe.SetMinimum(0.86);
  rframe.SetMaximum(2.16);
  rframe.GetXaxis()->SetTitle("Truth photon filter p_{T} [GeV]");
  rframe.GetYaxis()->SetTitle("stitched / fit");
  rframe.GetXaxis()->SetTitleSize(0.090);
  rframe.GetXaxis()->SetLabelSize(0.075);
  rframe.GetYaxis()->SetTitleSize(0.082);
  rframe.GetYaxis()->SetLabelSize(0.068);
  rframe.GetYaxis()->SetTitleOffset(0.55);
  rframe.GetYaxis()->SetNdivisions(505);
  rframe.Draw("AXIS");
  TLine one(10.0, 1.0, 35.0, 1.0);
  one.SetLineColor(kGray + 2);
  one.SetLineStyle(2);
  one.Draw("SAME");
  TLine b2(21.0, 0.86, 21.0, 2.16);
  b2.SetLineColor(kGray + 2);
  b2.SetLineStyle(3);
  b2.SetLineWidth(2);
  b2.Draw("SAME");
  rBounded->Draw("P SAME");
  rLower->Draw("P SAME");

  TLatex note;
  note.SetNDC(true);
  note.SetTextFont(42);
  note.SetTextSize(0.043);
  note.SetTextColor(kOrange + 7);
  note.DrawLatex(0.62, 0.77, "lower-only overlaps samples above 21 GeV");

  TLatex metric;
  metric.SetNDC(true);
  metric.SetTextFont(42);
  metric.SetTextSize(0.040);
  metric.DrawLatex(0.135, 0.18,
                   Form("bounded: <R>_{19-21}=%.3f, <R>_{21-23}=%.3f, jump=%.3f",
                        belowBounded, aboveBounded, aboveBounded / belowBounded));
  metric.SetTextColor(kOrange + 7);
  metric.DrawLatex(0.135, 0.075,
                   Form("lower-only: <R>_{19-21}=%.3f, <R>_{21-23}=%.3f, jump=%.3f",
                        belowLower, aboveLower, aboveLower / belowLower));

  const std::string png = kOutDir + "/photon12plus20_focus21_variant_comparison_x35.png";
  c.SaveAs(png.c_str());

  std::ofstream json(kOutDir + "/photon12plus20_focus21_variant_comparison_metrics.json");
  json << std::setprecision(12)
       << "{\n"
       << "  \"fit_reference\": \"bounded_12_21_ge21\",\n"
       << "  \"fit_function\": \"A*(1/x)^(b + c ln(x) + d x)\",\n"
       << "  \"bounded_mean_ratio_19_21\": " << belowBounded << ",\n"
       << "  \"bounded_mean_ratio_21_23\": " << aboveBounded << ",\n"
       << "  \"bounded_jump_21_over_19_21\": " << aboveBounded / belowBounded << ",\n"
       << "  \"loweronly_mean_ratio_19_21\": " << belowLower << ",\n"
       << "  \"loweronly_mean_ratio_21_23\": " << aboveLower << ",\n"
       << "  \"loweronly_jump_21_over_19_21\": " << aboveLower / belowLower << "\n"
       << "}\n";

  std::cout << "[DONE] " << png << std::endl;
}
