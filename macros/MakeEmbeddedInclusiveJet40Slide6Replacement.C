#include <TCanvas.h>
#include <TColor.h>
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
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace
{
const std::string kBase = "/Users/patsfan753/Desktop/ThesisAnalysis/";
const std::string kOutDir = kBase + "dataOutput/stitchDiagnostics/jet40_slide6_spectrum_20260526";
const std::string kCsv = kOutDir + "/inclusive4_stitch_bins.csv";
const std::string kVariant = "inclusive4";
const double kXMax = 50.0;
const int kJet12Color = kBlue + 1;
const int kJet20Color = TColor::GetColor("#ff7f0e");
const int kJet30Color = TColor::GetColor("#cc3399");
const int kJet40Color = TColor::GetColor("#2ca02c");

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

std::vector<Point> ReadVariant()
{
  std::ifstream in(kCsv);
  std::vector<Point> pts;
  std::string line;
  while (std::getline(in, line))
  {
    if (line.empty() || line[0] == '#') continue;
    if (line.rfind("variant,", 0) == 0) continue;
    const auto tok = Split(line, ',');
    if (tok.size() != 7 || tok[0] != kVariant) continue;
    const double lo = std::stod(tok[2]);
    const double hi = std::stod(tok[3]);
    Point p;
    p.x = std::stod(tok[4]);
    p.ex = 0.5 * (hi - lo);
    p.y = std::stod(tok[5]);
    p.ey = std::stod(tok[6]);
    if (p.x >= 12.0 && p.x <= kXMax && p.y > 0.0) pts.push_back(p);
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
  g->SetMarkerSize(1.05);
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

TGraph* MakeReference(const TF1* logFit)
{
  auto* g = new TGraph();
  int ip = 0;
  for (double x = 12.0; x <= kXMax + 0.0001; x += 0.05)
  {
    g->SetPoint(ip, x, RefEval(logFit, x));
    ++ip;
  }
  g->SetLineColor(kBlack);
  g->SetLineStyle(2);
  g->SetLineWidth(2);
  return g;
}

TGraphErrors* MakeRatio(const std::vector<Point>& pts, const TF1* logFit, const char* name, int color, int marker)
{
  auto* g = new TGraphErrors();
  g->SetName(name);
  int ip = 0;
  for (const auto& p : pts)
  {
    const double ref = RefEval(logFit, p.x);
    if (ref <= 0.0) continue;
    g->SetPoint(ip, p.x, p.y / ref);
    g->SetPointError(ip, p.ex, p.ey / ref);
    ++ip;
  }
  g->SetLineColor(color);
  g->SetMarkerColor(color);
  g->SetMarkerStyle(marker);
  g->SetMarkerSize(1.05);
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
    const double ref = RefEval(logFit, p.x);
    if (ref <= 0.0) continue;
    sum += p.y / ref;
    ++n;
  }
  return n > 0 ? sum / n : std::numeric_limits<double>::quiet_NaN();
}

void AddBox(double x1, double y1, double x2, double y2)
{
  auto* p = new TPaveText(x1, y1, x2, y2, "NDC");
  p->SetFillColor(TColor::GetColor("#f7f8fa"));
  p->SetLineColor(TColor::GetColor("#d9dde4"));
  p->SetShadowColor(0);
  p->SetBorderSize(1);
  p->Draw();
}

void Text(double x, double y, const char* s, double size, int color = kBlack,
          int align = 11, int font = 132)
{
  TLatex lat;
  lat.SetNDC(true);
  lat.SetTextFont(font);
  lat.SetTextSize(size);
  lat.SetTextColor(color);
  lat.SetTextAlign(align);
  lat.DrawLatex(x, y, s);
}
}

void MakeEmbeddedInclusiveJet40Slide6Replacement()
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetErrorX(0.5);
  gStyle->SetEndErrorSize(2);

  const auto pts = ReadVariant();
  if (pts.empty())
  {
    std::cerr << "[ERROR] No points read for " << kVariant << " from " << kCsv << std::endl;
    return;
  }

  std::unique_ptr<TGraphErrors> logGraph(MakeLogGraph(pts));
  TF1 logFit("log_modified_power_law",
             "[0] + ([1] + [2]*TMath::Log(x) + [3]*x)*TMath::Log(1.0/x)",
             12.0, kXMax);
  logFit.SetParameters(22.0, 5.0, 0.0, 0.0);
  logGraph->Fit(&logFit, "QNR");

  const auto p12 = SelectRange(pts, 12.0, 21.0);
  const auto p20 = SelectRange(pts, 21.0, 31.0);
  const auto p30 = SelectRange(pts, 31.0, 41.0);
  const auto p40 = SelectRange(pts, 41.0, kXMax + 0.1);
  std::unique_ptr<TGraphErrors> g12(MakeGraph(p12, "g_jet12", kJet12Color, 20));
  std::unique_ptr<TGraphErrors> g20(MakeGraph(p20, "g_jet20", kJet20Color, 21));
  std::unique_ptr<TGraphErrors> g30(MakeGraph(p30, "g_jet30", kJet30Color, 22));
  std::unique_ptr<TGraphErrors> g40(MakeGraph(p40, "g_jet40", kJet40Color, 23));
  std::unique_ptr<TGraph> ref(MakeReference(&logFit));
  std::unique_ptr<TGraphErrors> r12(MakeRatio(p12, &logFit, "r_jet12", kJet12Color, 20));
  std::unique_ptr<TGraphErrors> r20(MakeRatio(p20, &logFit, "r_jet20", kJet20Color, 21));
  std::unique_ptr<TGraphErrors> r30(MakeRatio(p30, &logFit, "r_jet30", kJet30Color, 22));
  std::unique_ptr<TGraphErrors> r40(MakeRatio(p40, &logFit, "r_jet40", kJet40Color, 23));

  const double r19_21 = MeanRatio(pts, &logFit, 19.0, 21.0);
  const double r21_23 = MeanRatio(pts, &logFit, 21.0, 23.0);
  const double jump21 = r21_23 / r19_21;
  const double r29_31 = MeanRatio(pts, &logFit, 29.0, 31.0);
  const double r31_33 = MeanRatio(pts, &logFit, 31.0, 33.0);
  const double jump31 = r31_33 / r29_31;
  const double r39_41 = MeanRatio(pts, &logFit, 39.0, 41.0);
  const double r41_43 = MeanRatio(pts, &logFit, 41.0, 43.0);
  const double jump41 = r41_43 / r39_41;

  gSystem->mkdir(kOutDir.c_str(), true);
  const std::string png = kOutDir + "/embeddedInclusiveJet12plus20plus30plus40_slide6_replacement.png";
  const std::string json = kOutDir + "/embeddedInclusiveJet12plus20plus30plus40_slide6_replacement_summary.json";

  TCanvas c("c_slide6_inclusive4", "Embedded Inclusive Jet 12+20+30+40 Stitching", 1600, 900);
  c.SetWindowSize(1600 + (1600 - c.GetWw()), 900 + (900 - c.GetWh()));
  c.SetCanvasSize(1600, 900);
  c.SetFillColor(kWhite);
  c.SetMargin(0, 0, 0, 0);

  Text(0.025, 0.935, "Embedded Inclusive Jet 12+20+30+40 Stitching", 0.055, kBlack, 11, 22);

  AddBox(0.030, 0.425, 0.470, 0.850);
  Text(0.042, 0.810, "Stitching rule", 0.037, kBlack, 11, 22);
  Text(0.042, 0.755, "Build one inclusive-jet generator spectrum from four", 0.024, kBlack, 11, 132);
  Text(0.042, 0.718, "exclusive samples. Each event is assigned by the", 0.024, kBlack, 11, 132);
  Text(0.042, 0.681, "highest-p_{T} generator jet passing the production filter.", 0.024, kBlack, 11, 132);
  Text(0.042, 0.615, "Ownership windows", 0.029, kBlack, 11, 22);
  Text(0.075, 0.565, "12-21", 0.030, kJet12Color, 21, 22);
  Text(0.160, 0.565, "21-31", 0.030, kJet20Color, 21, 22);
  Text(0.245, 0.565, "31-41", 0.030, kJet30Color, 21, 22);
  Text(0.330, 0.565, "#geq41", 0.030, kJet40Color, 21, 22);
  Text(0.075, 0.526, "Jet12", 0.022, kJet12Color, 21, 132);
  Text(0.160, 0.526, "Jet20", 0.022, kJet20Color, 21, 132);
  Text(0.245, 0.526, "Jet30", 0.022, kJet30Color, 21, 132);
  Text(0.330, 0.526, "Jet40", 0.022, kJet40Color, 21, 132);
  Text(0.075, 0.475, "p_{T}^{jet,filter} [GeV]", 0.020, kGray + 2, 11, 132);

  AddBox(0.030, 0.145, 0.470, 0.405);
  Text(0.042, 0.372, "Relative weights", 0.037, kBlack, 11, 22);
  Text(0.042, 0.332, "Each slice is scaled by #sigma_{eff}/N_{merged}; Jet40 is the anchor.", 0.020, kGray + 2, 11, 132);

  AddBox(0.052, 0.236, 0.230, 0.304);
  AddBox(0.250, 0.236, 0.428, 0.304);
  AddBox(0.052, 0.158, 0.230, 0.226);
  AddBox(0.250, 0.158, 0.428, 0.226);
  Text(0.070, 0.281, "Jet12", 0.024, kJet12Color, 11, 22);
  Text(0.155, 0.281, "9881#times", 0.024, kJet12Color, 11, 22);
  Text(0.070, 0.253, "12 #leq p_{T} < 21", 0.018, kBlack, 11, 132);
  Text(0.268, 0.281, "Jet20", 0.024, kJet20Color, 11, 22);
  Text(0.353, 0.281, "425.5#times", 0.024, kJet20Color, 11, 22);
  Text(0.268, 0.253, "21 #leq p_{T} < 31", 0.018, kBlack, 11, 132);
  Text(0.070, 0.203, "Jet30", 0.024, kJet30Color, 11, 22);
  Text(0.155, 0.203, "18.43#times", 0.024, kJet30Color, 11, 22);
  Text(0.070, 0.175, "31 #leq p_{T} < 41", 0.018, kBlack, 11, 132);
  Text(0.268, 0.203, "Jet40", 0.024, kJet40Color, 11, 22);
  Text(0.353, 0.203, "1.0#times", 0.024, kJet40Color, 11, 22);
  Text(0.268, 0.175, "p_{T} #geq 41", 0.018, kBlack, 11, 132);

  auto* top = new TPad("plot_top", "", 0.480, 0.390, 0.965, 0.850);
  auto* bot = new TPad("plot_bot", "", 0.480, 0.245, 0.965, 0.380);
  top->SetLeftMargin(0.12);
  top->SetRightMargin(0.035);
  top->SetTopMargin(0.10);
  top->SetBottomMargin(0.015);
  top->SetLogy(true);
  bot->SetLeftMargin(0.12);
  bot->SetRightMargin(0.035);
  bot->SetTopMargin(0.02);
  bot->SetBottomMargin(0.31);
  top->Draw();
  bot->Draw();

  top->cd();
  TH1F frame("frame", "", 100, 12.0, kXMax);
  frame.SetMinimum(1.0);
  frame.SetMaximum(1.2e6);
  frame.GetYaxis()->SetTitle("#sigma_{eff} #times N scaled entries [pb / bin]");
  frame.GetYaxis()->SetTitleFont(132);
  frame.GetYaxis()->SetLabelFont(132);
  frame.GetYaxis()->SetTitleSize(0.060);
  frame.GetYaxis()->SetLabelSize(0.050);
  frame.GetYaxis()->SetTitleOffset(1.00);
  frame.GetXaxis()->SetLabelSize(0.0);
  frame.GetXaxis()->SetTitleSize(0.0);
  frame.SetStats(false);
  frame.Draw("AXIS");
  ref->Draw("L SAME");
  g12->Draw("PZ SAME");
  g20->Draw("PZ SAME");
  g30->Draw("PZ SAME");
  g40->Draw("PZ SAME");

  TLine b21top(21.0, 1.0, 21.0, 1.2e6);
  b21top.SetX1(21.0);
  b21top.SetX2(21.0);
  b21top.SetY1(1.0);
  b21top.SetY2(1.2e6);
  b21top.SetLineColor(kGray + 1);
  b21top.SetLineStyle(3);
  b21top.SetLineWidth(2);
  b21top.Draw("SAME");
  TLine b31top(31.0, 1.0, 31.0, 1.2e6);
  b31top.SetLineColor(kGray + 1);
  b31top.SetLineStyle(3);
  b31top.SetLineWidth(2);
  b31top.Draw("SAME");
  TLine b41top(41.0, 1.0, 41.0, 1.2e6);
  b41top.SetLineColor(kGray + 1);
  b41top.SetLineStyle(3);
  b41top.SetLineWidth(2);
  b41top.Draw("SAME");

  TLatex plat;
  plat.SetNDC(true);
  plat.SetTextFont(132);
  plat.SetTextSize(0.044);
  plat.DrawLatex(0.150, 0.265, "Embedded inclusive-jet stitch");
  plat.SetTextSize(0.031);
  plat.DrawLatex(0.150, 0.210, "Jet12: 12 #leq p_{T}^{jet,filter} < 21 GeV");
  plat.DrawLatex(0.150, 0.165, "Jet20: 21 #leq p_{T}^{jet,filter} < 31 GeV");
  plat.DrawLatex(0.150, 0.120, "Jet30: 31 #leq p_{T}^{jet,filter} < 41 GeV");
  plat.DrawLatex(0.150, 0.075, "Jet40: p_{T}^{jet,filter} #geq 41 GeV");
  plat.SetTextAlign(11);
  plat.SetTextSize(0.048);
  plat.DrawLatex(0.790, 0.925, "#it{#bf{sPHENIX}}");
  plat.SetTextSize(0.046);
  plat.DrawLatex(0.890, 0.925, "Internal");
  plat.SetTextAlign(33);
  plat.SetTextSize(0.040);
  plat.DrawLatex(0.965, 0.865, "PYTHIA8 embedded inclusive jets, #sqrt{s_{NN}} = 200 GeV");
  plat.SetTextAlign(11);

  TLegend leg(0.61, 0.36, 0.94, 0.72);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextFont(132);
  leg.SetTextSize(0.035);
  leg.AddEntry(g12.get(), "Jet12 stitched", "lep");
  leg.AddEntry(g20.get(), "Jet20 stitched", "lep");
  leg.AddEntry(g30.get(), "Jet30 stitched", "lep");
  leg.AddEntry(g40.get(), "Jet40 stitched", "lep");
  leg.AddEntry(ref.get(), "Fit: modified power law", "l");
  leg.Draw();

  bot->cd();
  TH1F rframe("rframe", "", 100, 12.0, kXMax);
  rframe.SetMinimum(0.88);
  rframe.SetMaximum(1.12);
  rframe.GetYaxis()->SetTitle("stitched / fit");
  rframe.GetXaxis()->SetTitle("p_{T}^{jet,filter} [GeV]");
  rframe.GetYaxis()->SetTitleFont(132);
  rframe.GetYaxis()->SetLabelFont(132);
  rframe.GetXaxis()->SetTitleFont(132);
  rframe.GetXaxis()->SetLabelFont(132);
  rframe.GetYaxis()->SetTitleSize(0.115);
  rframe.GetYaxis()->SetLabelSize(0.100);
  rframe.GetYaxis()->SetTitleOffset(0.46);
  rframe.GetYaxis()->SetNdivisions(505);
  rframe.GetXaxis()->SetTitleSize(0.125);
  rframe.GetXaxis()->SetLabelSize(0.100);
  rframe.GetXaxis()->SetTitleOffset(1.00);
  rframe.SetStats(false);
  rframe.Draw("AXIS");
  TLine one(12.0, 1.0, kXMax, 1.0);
  one.SetLineColor(kGray + 2);
  one.SetLineStyle(2);
  one.Draw("SAME");
  r12->Draw("PZ SAME");
  r20->Draw("PZ SAME");
  r30->Draw("PZ SAME");
  r40->Draw("PZ SAME");
  TLine b21bot(21.0, 0.88, 21.0, 1.12);
  b21bot.SetLineColor(kGray + 1);
  b21bot.SetLineStyle(3);
  b21bot.SetLineWidth(2);
  b21bot.Draw("SAME");
  TLine b31bot(31.0, 0.88, 31.0, 1.12);
  b31bot.SetLineColor(kGray + 1);
  b31bot.SetLineStyle(3);
  b31bot.SetLineWidth(2);
  b31bot.Draw("SAME");
  TLine b41bot(41.0, 0.88, 41.0, 1.12);
  b41bot.SetLineColor(kGray + 1);
  b41bot.SetLineStyle(3);
  b41bot.SetLineWidth(2);
  b41bot.Draw("SAME");

  c.cd();
  AddBox(0.030, 0.030, 0.955, 0.125);
  auto* sep = new TLine(0.03, 0.145, 0.955, 0.145);
  sep->SetNDC(true);
  sep->SetLineColor(TColor::GetColor("#e3e5e8"));
  sep->Draw();
  Text(0.500, 0.094,
       "Bounded 21/31/41 GeV ownership removes double counting while keeping the stitched spectrum continuous.",
       0.024, kBlack, 21, 132);
  Text(0.500, 0.055,
       Form("Boundary checks: 21 GeV jump = %.3f; 31 GeV jump = %.3f; 41 GeV jump = %.3f.",
            jump21, jump31, jump41),
       0.019, kGray + 2, 21, 132);

  c.SaveAs(png.c_str());

  std::ofstream out(json);
  out << std::setprecision(12)
      << "{\n"
      << "  \"variant\": \"inclusive4\",\n"
      << "  \"source_base\": \"/sphenix/u/patsfan753/scratch/thesisAnalysis/output/simembeddedinclusive/preselectionReference_tightReference_nonTightReference_baseVariant\",\n"
      << "  \"source_root\": \"/sphenix/u/patsfan753/scratch/thesisAnalysis/output/simembeddedinclusive/preselectionReference_tightReference_nonTightReference_baseVariant/embeddedJet12and20and30and40merged_SIM/RecoilJets_embeddedJet12plus20plus30plus40_MERGED.root\",\n"
      << "  \"csv\": \"" << kCsv << "\",\n"
      << "  \"jet12_gate\": \"12 <= pT_jet_filter < 21\",\n"
      << "  \"jet20_gate\": \"21 <= pT_jet_filter < 31\",\n"
      << "  \"jet30_gate\": \"31 <= pT_jet_filter < 41\",\n"
      << "  \"jet40_gate\": \"pT_jet_filter >= 41\",\n"
      << "  \"jet12_merged_entries\": 7213068,\n"
      << "  \"jet20_merged_entries\": 5295096,\n"
      << "  \"jet30_merged_entries\": 5471684,\n"
      << "  \"jet40_merged_entries\": 5842568,\n"
      << "  \"jet12_sigma_eff_pb\": 1227724.77,\n"
      << "  \"jet20_sigma_eff_pb\": 38811.785,\n"
      << "  \"jet30_sigma_eff_pb\": 1736.65908,\n"
      << "  \"jet40_sigma_eff_pb\": 100.642312,\n"
      << "  \"jet12_relative_weight_to_jet40\": 9881.074239833762,\n"
      << "  \"jet20_relative_weight_to_jet40\": 425.5131140895668,\n"
      << "  \"jet30_relative_weight_to_jet40\": 18.425391902309496,\n"
      << "  \"jet40_relative_weight_to_jet40\": 1.0,\n"
      << "  \"fit_function\": \"A*(1/x)^(b + c ln(x) + d x)\",\n"
      << "  \"mean_ratio_19_21\": " << r19_21 << ",\n"
      << "  \"mean_ratio_21_23\": " << r21_23 << ",\n"
      << "  \"jump_21_over_19_21\": " << jump21 << ",\n"
      << "  \"mean_ratio_29_31\": " << r29_31 << ",\n"
      << "  \"mean_ratio_31_33\": " << r31_33 << ",\n"
      << "  \"jump_31_over_29_31\": " << jump31 << ",\n"
      << "  \"mean_ratio_39_41\": " << r39_41 << ",\n"
      << "  \"mean_ratio_41_43\": " << r41_43 << ",\n"
      << "  \"jump_41_over_39_41\": " << jump41 << ",\n"
      << "  \"png\": \"" << png << "\"\n"
      << "}\n";

  std::cout << "[DONE] " << png << std::endl;
  std::cout << "[SUMMARY] 21 jump=" << jump21 << " 31 jump=" << jump31 << " 41 jump=" << jump41 << std::endl;
}
