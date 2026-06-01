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
#include <sstream>
#include <string>
#include <vector>

namespace
{
const std::string kBase = "/Users/patsfan753/Desktop/ThesisAnalysis/";
const std::string kOutDir = kBase + "dataOutput/stitchDiagnostics/focus21_photon_variants_20260521";
const std::string kCsv = kOutDir + "/photon_bounded21_extended_bins.csv";
const std::string kVariant = "bounded_12_21_ge21";
const double kXMax = 40.0;
const double kPlotBinWidth = 1.0;

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
  std::vector<Point> rebinned;
  for (double lo = 12.0; lo < kXMax - 1.0e-6; lo += kPlotBinWidth)
  {
    const double hi = lo + kPlotBinWidth;
    Point q;
    q.x = 0.5 * (lo + hi);
    q.ex = 0.5 * (hi - lo);
    double err2 = 0.0;
    for (const auto& p : pts)
    {
      const double plo = p.x - p.ex;
      const double phi = p.x + p.ex;
      if (plo >= lo - 1.0e-6 && phi <= hi + 1.0e-6)
      {
        q.y += p.y;
        err2 += p.ey * p.ey;
      }
    }
    q.ey = std::sqrt(err2);
    if (q.y > 0.0) rebinned.push_back(q);
  }
  return rebinned;
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
  g->SetMarkerSize(0.90);
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
  g->SetMarkerSize(0.90);
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

void MakePhotonJet12Plus20Slide4Replacement()
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
  logFit.SetParameters(24.0, 5.0, 0.0, 0.0);
  logGraph->Fit(&logFit, "QNR");

  const auto p12 = SelectRange(pts, 12.0, 21.0);
  const auto p20 = SelectRange(pts, 21.0, kXMax + 0.1);
  std::unique_ptr<TGraphErrors> g12(MakeGraph(p12, "g_photon12", kBlue + 1, 20));
  std::unique_ptr<TGraphErrors> g20(MakeGraph(p20, "g_photon20", kRed + 1, 21));
  std::unique_ptr<TGraph> ref(MakeReference(&logFit));
  std::unique_ptr<TGraphErrors> r12(MakeRatio(p12, &logFit, "r_photon12", kBlue + 1, 20));
  std::unique_ptr<TGraphErrors> r20(MakeRatio(p20, &logFit, "r_photon20", kRed + 1, 21));

  const double rBelow = MeanRatio(pts, &logFit, 19.0, 21.0);
  const double rAbove = MeanRatio(pts, &logFit, 21.0, 23.0);
  const double jump = rAbove / rBelow;

  gSystem->mkdir(kOutDir.c_str(), true);
  const std::string png = kOutDir + "/photon12plus20_bounded21_slide4_replacement.png";
  const std::string json = kOutDir + "/photon12plus20_bounded21_slide4_replacement_summary.json";

  TCanvas c("c_slide4_photon21", "Embedded Photon Jet 12+20 Stitching", 1600, 900);
  c.SetWindowSize(1600 + (1600 - c.GetWw()), 900 + (900 - c.GetWh()));
  c.SetCanvasSize(1600, 900);
  c.SetFillColor(kWhite);
  c.SetMargin(0, 0, 0, 0);

  Text(0.025, 0.935, "Embedded Photon Jet 12+20 Stitching", 0.058, kBlack, 11, 22);

  AddBox(0.030, 0.472, 0.470, 0.850);
  Text(0.042, 0.818, "Generator reproduction + event ownership", 0.031, kBlack, 11, 22);
  Text(0.042, 0.775, "Base path:", 0.020, kBlack, 11, 132);
  Text(0.115, 0.775, "/sphenix/tg/tg01/commissioning/CaloCalibWG/bseidlitz/embed_2025", 0.018, kGray + 2, 11, 132);
  Text(0.042, 0.732, "Configs:", 0.024, kBlack, 11, 132);
  Text(0.063, 0.699, "#bullet  PhotonJet12 #rightarrow phpythia8_10GeV_JS_MDC2.cfg", 0.020, kBlack, 11, 132);
  Text(0.063, 0.667, "#bullet  PhotonJet20 #rightarrow phpythia8_20GeV_JS_MDC2.cfg", 0.020, kBlack, 11, 132);
  Text(0.042, 0.624, "Photon filter: PDG = 22, |#eta| < 1.5, 1 #leq |mother PDG| #leq 22", 0.019, kBlack, 11, 132);
  Text(0.042, 0.594, "Define: p_{T}(#gamma, filter) = highest-p_{T} generator photon passing the filter", 0.019, kBlack, 11, 132);
  Text(0.042, 0.550, "Ownership:", 0.024, kBlack, 11, 132);
  Text(0.063, 0.519, "#bullet  PhotonJet12: 12 #leq p_{T}(#gamma, filter) < 21", 0.020, kBlue + 1, 11, 132);
  Text(0.063, 0.489, "#bullet  PhotonJet20: p_{T}(#gamma, filter) #geq 21", 0.020, kRed + 1, 11, 132);

  AddBox(0.030, 0.245, 0.470, 0.465);
  Text(0.042, 0.430, "Derived stitched weights", 0.031, kBlack, 11, 22);
  Text(0.042, 0.390, "Sample", 0.020, kBlack, 11, 132);
  Text(0.150, 0.390, "Ownership region", 0.020, kBlack, 11, 132);
  Text(0.280, 0.390, "Npass", 0.020, kBlack, 11, 132);
  Text(0.340, 0.390, "#sigma_{eff}", 0.020, kBlack, 11, 132);
  auto* hline = new TLine(0.042, 0.374, 0.445, 0.374);
  hline->SetNDC(true);
  hline->SetLineColor(kGray + 1);
  hline->Draw();
  Text(0.042, 0.344, "PhotonJet12", 0.020, kBlue + 1, 11, 132);
  Text(0.150, 0.344, "12 #leq p_{T}(#gamma, filter) < 21", 0.020, kBlack, 11, 132);
  Text(0.280, 0.344, "1,056/50M", 0.020, kBlack, 11, 132);
  Text(0.340, 0.344, "2663.5 #pm 82.0 pb", 0.020, kBlack, 11, 132);
  Text(0.042, 0.307, "PhotonJet20", 0.020, kRed + 1, 11, 132);
  Text(0.150, 0.307, "p_{T}(#gamma, filter) #geq 21", 0.020, kBlack, 11, 132);
  Text(0.280, 0.307, "881/50M", 0.020, kBlack, 11, 132);
  Text(0.340, 0.307, "92.05 #pm 3.10 pb", 0.020, kBlack, 11, 132);
  Text(0.060, 0.260, "#sigma_{eff} = #sigma_{gen} #times N_{pass} / N_{thrown}", 0.026, kBlack, 11, 132);

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
  frame.SetMinimum(2.0e3);
  frame.SetMaximum(1.6e8);
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

  TLine btop(21.0, 2.0e3, 21.0, 1.6e8);
  btop.SetLineColor(kGray + 1);
  btop.SetLineStyle(3);
  btop.SetLineWidth(2);
  btop.Draw("SAME");

  TLatex plat;
  plat.SetNDC(true);
  plat.SetTextFont(132);
  plat.SetTextSize(0.047);
  plat.DrawLatex(0.155, 0.245, "Embedded photon generator stitching spectrum");
  plat.SetTextSize(0.034);
  plat.DrawLatex(0.155, 0.185, "PhotonJet12: 12 #leq p_{T}^{#gamma,filter} < 21 GeV");
  plat.DrawLatex(0.155, 0.135, "PhotonJet20: p_{T}^{#gamma,filter} #geq 21 GeV");
  plat.SetTextAlign(11);
  plat.SetTextSize(0.048);
  plat.DrawLatex(0.790, 0.925, "#it{#bf{sPHENIX}}");
  plat.SetTextSize(0.046);
  plat.DrawLatex(0.890, 0.925, "Internal");
  plat.SetTextAlign(33);
  plat.SetTextSize(0.040);
  plat.DrawLatex(0.965, 0.865, "PYTHIA8 embedded photon+jet, #sqrt{s_{NN}} = 200 GeV");
  plat.SetTextAlign(11);

  TLegend leg(0.60, 0.42, 0.94, 0.71);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextFont(132);
  leg.SetTextSize(0.036);
  leg.AddEntry(g12.get(), "PhotonJet12 stitched", "lep");
  leg.AddEntry(g20.get(), "PhotonJet20 stitched", "lep");
  leg.AddEntry(ref.get(), "Fit: modified power law", "l");
  leg.Draw();

  bot->cd();
  TH1F rframe("rframe", "", 100, 12.0, kXMax);
  rframe.SetMinimum(0.88);
  rframe.SetMaximum(1.12);
  rframe.GetYaxis()->SetTitle("stitched / fit");
  rframe.GetXaxis()->SetTitle("p_{T}^{#gamma,filter} [GeV]");
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
  TLine bbot(21.0, 0.88, 21.0, 1.12);
  bbot.SetLineColor(kGray + 1);
  bbot.SetLineStyle(3);
  bbot.SetLineWidth(2);
  bbot.Draw("SAME");

  c.cd();
  AddBox(0.030, 0.045, 0.955, 0.165);
  auto* sep = new TLine(0.03, 0.185, 0.955, 0.185);
  sep->SetNDC(true);
  sep->SetLineColor(TColor::GetColor("#e3e5e8"));
  sep->Draw();
  Text(0.500, 0.112,
       "The 21 GeV ownership rule removes double counting and gives a smooth PhotonJet12 #rightarrow PhotonJet20 transition.",
       0.030, kBlack, 21, 132);
  Text(0.500, 0.070,
       Form("Boundary check: 19-21 GeV average = %.3f, 21-23 GeV average = %.3f, jump = %.3f (+%.1f%%).",
            rBelow, rAbove, jump, 100.0 * (jump - 1.0)),
       0.022, kGray + 2, 21, 132);

  c.SaveAs(png.c_str());

  std::ofstream out(json);
  out << std::setprecision(12)
      << "{\n"
      << "  \"variant\": \"bounded_12_21_ge21\",\n"
      << "  \"photonjet12_gate\": \"12 <= pT_gamma_filter < 21\",\n"
      << "  \"photonjet20_gate\": \"pT_gamma_filter >= 21\",\n"
      << "  \"photonjet12_npass\": 1056,\n"
      << "  \"photonjet20_npass\": 881,\n"
      << "  \"n_thrown_each\": 50000000,\n"
      << "  \"photonjet12_sigma_eff_pb\": 2663.51030,\n"
      << "  \"photonjet12_sigma_eff_stat_unc_pb\": 82.0,\n"
      << "  \"photonjet20_sigma_eff_pb\": 92.0516019,\n"
      << "  \"photonjet20_sigma_eff_stat_unc_pb\": 3.10,\n"
      << "  \"photonjet12_merge_scale\": 21.265641283645,\n"
      << "  \"photonjet20_merge_scale\": 1.0,\n"
      << "  \"fit_function\": \"A*(1/x)^(b + c ln(x) + d x)\",\n"
      << "  \"mean_ratio_19_21\": " << rBelow << ",\n"
      << "  \"mean_ratio_21_23\": " << rAbove << ",\n"
      << "  \"jump_21_over_19_21\": " << jump << ",\n"
      << "  \"png\": \"" << png << "\"\n"
      << "}\n";

  std::cout << "[DONE] " << png << std::endl;
  std::cout << "[SUMMARY] 19-21=" << rBelow << " 21-23=" << rAbove
            << " jump=" << jump << std::endl;
}
