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
const std::string kOutDir = kBase + "dataOutput/stitchDiagnostics/focus21_inclusive_variants_20260521";
const std::string kCsv = kOutDir + "/inclusive_exclusive31_stitch_bins.csv";
const std::string kVariant = "exclusive31";
const double kXMax = 40.0;

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

void MakeEmbeddedInclusiveJetSlide5Replacement()
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
  const auto p30 = SelectRange(pts, 31.0, kXMax + 0.1);
  std::unique_ptr<TGraphErrors> g12(MakeGraph(p12, "g_jet12", kBlue + 1, 20));
  std::unique_ptr<TGraphErrors> g20(MakeGraph(p20, "g_jet20", TColor::GetColor("#ff7f0e"), 21));
  std::unique_ptr<TGraphErrors> g30(MakeGraph(p30, "g_jet30", TColor::GetColor("#cc3399"), 22));
  std::unique_ptr<TGraph> ref(MakeReference(&logFit));
  std::unique_ptr<TGraphErrors> r12(MakeRatio(p12, &logFit, "r_jet12", kBlue + 1, 20));
  std::unique_ptr<TGraphErrors> r20(MakeRatio(p20, &logFit, "r_jet20", TColor::GetColor("#ff7f0e"), 21));
  std::unique_ptr<TGraphErrors> r30(MakeRatio(p30, &logFit, "r_jet30", TColor::GetColor("#cc3399"), 22));

  const double r19_21 = MeanRatio(pts, &logFit, 19.0, 21.0);
  const double r21_23 = MeanRatio(pts, &logFit, 21.0, 23.0);
  const double jump21 = r21_23 / r19_21;
  const double r29_31 = MeanRatio(pts, &logFit, 29.0, 31.0);
  const double r31_33 = MeanRatio(pts, &logFit, 31.0, 33.0);
  const double jump31 = r31_33 / r29_31;

  gSystem->mkdir(kOutDir.c_str(), true);
  const std::string png = kOutDir + "/embeddedInclusiveJet12plus20plus30_exclusive31_slide5_replacement.png";
  const std::string json = kOutDir + "/embeddedInclusiveJet12plus20plus30_exclusive31_slide5_replacement_summary.json";

  TCanvas c("c_slide5_inclusive31", "Embedded Inclusive Jet 12+20+30 Stitching", 1600, 900);
  c.SetWindowSize(1600 + (1600 - c.GetWw()), 900 + (900 - c.GetWh()));
  c.SetCanvasSize(1600, 900);
  c.SetFillColor(kWhite);
  c.SetMargin(0, 0, 0, 0);

  Text(0.025, 0.935, "Embedded Inclusive Jet 12+20+30 Stitching", 0.058, kBlack, 11, 22);

  AddBox(0.030, 0.448, 0.470, 0.850);
  Text(0.042, 0.818, "Generator reproduction + event ownership", 0.031, kBlack, 11, 22);
  Text(0.042, 0.775, "Base path:", 0.020, kBlack, 11, 132);
  Text(0.115, 0.775, "/sphenix/tg/tg01/commissioning/CaloCalibWG/bseidlitz/embed_2025", 0.018, kGray + 2, 11, 132);
  Text(0.042, 0.732, "Configs:", 0.024, kBlack, 11, 132);
  Text(0.063, 0.699, "#bullet  Jet12 #rightarrow phpythia8_10GeV_JS_MDC2.cfg", 0.020, kBlack, 11, 132);
  Text(0.063, 0.667, "#bullet  Jet20 #rightarrow phpythia8_20GeV_JS_MDC2.cfg", 0.020, kBlack, 11, 132);
  Text(0.063, 0.635, "#bullet  Jet30 #rightarrow phpythia8_30GeV_JS_MDC2.cfg", 0.020, kBlack, 11, 132);
  Text(0.042, 0.592, "Jet filter: generator-level inclusive-jet sample", 0.019, kBlack, 11, 132);
  Text(0.042, 0.562, "Define: p_{T}(jet, filter) = highest-p_{T} generator jet passing the filter", 0.019, kBlack, 11, 132);
  Text(0.042, 0.524, "Ownership:", 0.024, kBlack, 11, 132);
  Text(0.063, 0.498, "#bullet  Jet12: 12 #leq p_{T}(jet, filter) < 21", 0.018, kBlue + 1, 11, 132);
  Text(0.063, 0.475, "#bullet  Jet20: 21 #leq p_{T}(jet, filter) < 31", 0.018, TColor::GetColor("#cc5f00"), 11, 132);
  Text(0.063, 0.452, "#bullet  Jet30: p_{T}(jet, filter) #geq 31", 0.018, TColor::GetColor("#99006f"), 11, 132);

  AddBox(0.030, 0.195, 0.470, 0.435);
  Text(0.042, 0.402, "Derived stitched weights", 0.031, kBlack, 11, 22);
  Text(0.042, 0.362, "Sample", 0.020, kBlack, 11, 132);
  Text(0.126, 0.362, "Ownership region", 0.020, kBlack, 11, 132);
  Text(0.245, 0.362, "Npass", 0.020, kBlack, 11, 132);
  Text(0.318, 0.362, "#sigma_{eff}", 0.020, kBlack, 11, 132);
  auto* hline = new TLine(0.042, 0.346, 0.445, 0.346);
  hline->SetNDC(true);
  hline->SetLineColor(kGray + 1);
  hline->Draw();
  Text(0.042, 0.316, "Jet12", 0.019, kBlue + 1, 11, 132);
  Text(0.126, 0.316, "12 #leq p_{T} < 21", 0.019, kBlack, 11, 132);
  Text(0.245, 0.316, "487,343/50M", 0.019, kBlack, 11, 132);
  Text(0.318, 0.316, "1.229e6 #pm 1.76e3 pb", 0.019, kBlack, 11, 132);
  Text(0.042, 0.281, "Jet20", 0.019, TColor::GetColor("#cc5f00"), 11, 132);
  Text(0.126, 0.281, "21 #leq p_{T} < 31", 0.019, kBlack, 11, 132);
  Text(0.245, 0.281, "371,220/50M", 0.019, kBlack, 11, 132);
  Text(0.318, 0.281, "3.878e4 #pm 63.6 pb", 0.019, kBlack, 11, 132);
  Text(0.042, 0.246, "Jet30", 0.019, TColor::GetColor("#99006f"), 11, 132);
  Text(0.126, 0.246, "p_{T} #geq 31", 0.019, kBlack, 11, 132);
  Text(0.245, 0.246, "203,163/50M", 0.019, kBlack, 11, 132);
  Text(0.318, 0.246, "1793.0 #pm 4.0 pb", 0.019, kBlack, 11, 132);
  Text(0.060, 0.210, "#sigma_{eff} = #sigma_{gen} #times N_{pass} / N_{thrown}", 0.019, kBlack, 11, 132);

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
  frame.SetMinimum(3.0e1);
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

  TLine b21top(21.0, 3.0e1, 21.0, 1.2e6);
  b21top.SetX1(21.0);
  b21top.SetX2(21.0);
  b21top.SetY1(3.0e1);
  b21top.SetY2(1.2e6);
  b21top.SetLineColor(kGray + 1);
  b21top.SetLineStyle(3);
  b21top.SetLineWidth(2);
  b21top.Draw("SAME");
  TLine b31top(31.0, 3.0e1, 31.0, 1.2e6);
  b31top.SetLineColor(kGray + 1);
  b31top.SetLineStyle(3);
  b31top.SetLineWidth(2);
  b31top.Draw("SAME");

  TLatex plat;
  plat.SetNDC(true);
  plat.SetTextFont(132);
  plat.SetTextSize(0.044);
  plat.DrawLatex(0.150, 0.265, "Embedded inclusive-jet stitch");
  plat.SetTextSize(0.031);
  plat.DrawLatex(0.150, 0.210, "Jet12: 12 #leq p_{T}^{jet,filter} < 21 GeV");
  plat.DrawLatex(0.150, 0.165, "Jet20: 21 #leq p_{T}^{jet,filter} < 31 GeV");
  plat.DrawLatex(0.150, 0.120, "Jet30: p_{T}^{jet,filter} #geq 31 GeV");
  plat.SetTextAlign(11);
  plat.SetTextSize(0.048);
  plat.DrawLatex(0.790, 0.925, "#it{#bf{sPHENIX}}");
  plat.SetTextSize(0.046);
  plat.DrawLatex(0.890, 0.925, "Internal");
  plat.SetTextAlign(33);
  plat.SetTextSize(0.040);
  plat.DrawLatex(0.965, 0.865, "PYTHIA8 embedded inclusive jets, #sqrt{s_{NN}} = 200 GeV");
  plat.SetTextAlign(11);

  TLegend leg(0.58, 0.39, 0.94, 0.72);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextFont(132);
  leg.SetTextSize(0.035);
  leg.AddEntry(g12.get(), "Jet12 stitched", "lep");
  leg.AddEntry(g20.get(), "Jet20 stitched", "lep");
  leg.AddEntry(g30.get(), "Jet30 stitched", "lep");
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

  c.cd();
  AddBox(0.030, 0.045, 0.955, 0.165);
  auto* sep = new TLine(0.03, 0.185, 0.955, 0.185);
  sep->SetNDC(true);
  sep->SetLineColor(TColor::GetColor("#e3e5e8"));
  sep->Draw();
  Text(0.500, 0.112,
       "The bounded 21/31 GeV ownership rule removes double counting across the Jet12 #rightarrow Jet20 #rightarrow Jet30 stitched spectrum.",
       0.028, kBlack, 21, 132);
  Text(0.500, 0.070,
       Form("Boundary checks: 21 GeV jump = %.3f (+%.1f%%); 31 GeV jump = %.3f (+%.1f%%).",
            jump21, 100.0 * (jump21 - 1.0), jump31, 100.0 * (jump31 - 1.0)),
       0.022, kGray + 2, 21, 132);

  c.SaveAs(png.c_str());

  std::ofstream out(json);
  out << std::setprecision(12)
      << "{\n"
      << "  \"variant\": \"exclusive31\",\n"
      << "  \"source_root\": \"/sphenix/u/patsfan753/scratch/thesisAnalysis/output_focus21_stitch_20260521_123833_exclusive31/simembeddedinclusive/preselectionReference_tightReference_nonTightReference_baseVariant/embeddedJet12and20and30merged_SIM/RecoilJets_embeddedJet12plus20plus30_MERGED.root\",\n"
      << "  \"jet12_gate\": \"12 <= pT_jet_filter < 21\",\n"
      << "  \"jet20_gate\": \"21 <= pT_jet_filter < 31\",\n"
      << "  \"jet30_gate\": \"pT_jet_filter >= 31\",\n"
      << "  \"jet12_npass\": 487343,\n"
      << "  \"jet20_npass\": 371220,\n"
      << "  \"jet30_npass\": 203163,\n"
      << "  \"n_thrown_each\": 50000000,\n"
      << "  \"jet12_sigma_eff_pb\": 1229410.38,\n"
      << "  \"jet12_sigma_eff_stat_unc_pb\": 1760.7,\n"
      << "  \"jet20_sigma_eff_pb\": 38784.6122,\n"
      << "  \"jet20_sigma_eff_stat_unc_pb\": 63.6,\n"
      << "  \"jet30_sigma_eff_pb\": 1793.04310,\n"
      << "  \"jet30_sigma_eff_stat_unc_pb\": 4.0,\n"
      << "  \"jet12_merge_scale\": 549.955492198160,\n"
      << "  \"jet20_merge_scale\": 23.633949843035,\n"
      << "  \"jet30_merge_scale\": 1.0,\n"
      << "  \"fit_function\": \"A*(1/x)^(b + c ln(x) + d x)\",\n"
      << "  \"mean_ratio_19_21\": " << r19_21 << ",\n"
      << "  \"mean_ratio_21_23\": " << r21_23 << ",\n"
      << "  \"jump_21_over_19_21\": " << jump21 << ",\n"
      << "  \"mean_ratio_29_31\": " << r29_31 << ",\n"
      << "  \"mean_ratio_31_33\": " << r31_33 << ",\n"
      << "  \"jump_31_over_29_31\": " << jump31 << ",\n"
      << "  \"png\": \"" << png << "\"\n"
      << "}\n";

  std::cout << "[DONE] " << png << std::endl;
  std::cout << "[SUMMARY] 21 jump=" << jump21 << " 31 jump=" << jump31 << std::endl;
}
