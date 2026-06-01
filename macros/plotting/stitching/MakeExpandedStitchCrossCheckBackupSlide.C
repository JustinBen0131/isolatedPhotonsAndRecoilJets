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
#include <TObject.h>

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
const std::string kOutDir = kBase + "dataOutput/stitchDiagnostics/focus21_backup_crosscheck_20260521";
const std::string kPhotonCsv = kBase + "dataOutput/stitchDiagnostics/focus21_photon_variants_20260521/photon_bounded21_extended_bins.csv";
const std::string kInclusiveCsv = kBase + "dataOutput/stitchDiagnostics/focus21_inclusive_variants_20260521/inclusive_exclusive31_stitch_bins.csv";

struct Point
{
  double x = 0.0;
  double ex = 0.0;
  double y = 0.0;
  double ey = 0.0;
};

struct PanelConfig
{
  std::string csv;
  std::string variant;
  std::string title;
  std::string subtitle;
  std::string xtitle;
  double xmax = 40.0;
  double ymin = 1.0;
  double ymax = 1.0e8;
  double boundary1 = 21.0;
  double boundary2 = -1.0;
  double rebinWidth = 0.0;
  std::string leg1;
  std::string leg2;
  std::string leg3;
  std::string gate1;
  std::string gate2;
  std::string gate3;
  int color1 = kBlue + 1;
  int color2 = kRed + 1;
  int color3 = kMagenta + 1;
};

std::vector<std::string> Split(const std::string& s, char delim)
{
  std::vector<std::string> out;
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, delim)) out.push_back(item);
  return out;
}

std::vector<Point> ReadVariant(const PanelConfig& cfg)
{
  std::ifstream in(cfg.csv);
  std::vector<Point> pts;
  std::string line;
  while (std::getline(in, line))
  {
    if (line.empty() || line[0] == '#') continue;
    if (line.rfind("variant,", 0) == 0) continue;
    const auto tok = Split(line, ',');
    if (tok.size() != 7 || tok[0] != cfg.variant) continue;
    const double lo = std::stod(tok[2]);
    const double hi = std::stod(tok[3]);
    Point p;
    p.x = std::stod(tok[4]);
    p.ex = 0.5 * (hi - lo);
    p.y = std::stod(tok[5]);
    p.ey = std::stod(tok[6]);
    if (p.x >= 12.0 && p.x <= cfg.xmax && p.y > 0.0) pts.push_back(p);
  }
  if (cfg.rebinWidth > 0.0)
  {
    std::vector<Point> rebinned;
    for (double lo = 12.0; lo < cfg.xmax - 1.0e-6; lo += cfg.rebinWidth)
    {
      const double hi = lo + cfg.rebinWidth;
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
  g->SetMarkerSize(0.85);
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

TGraph* MakeReference(const TF1* logFit, double xmax)
{
  auto* g = new TGraph();
  int ip = 0;
  for (double x = 12.0; x <= xmax + 1.0e-4; x += 0.05)
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
  g->SetMarkerSize(0.85);
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

void DrawPanel(const PanelConfig& cfg, double x1, double y1, double x2, double y2,
               double& jump1, double& jump2)
{
  const auto pts = ReadVariant(cfg);
  if (pts.empty())
  {
    std::cerr << "[ERROR] No points for " << cfg.variant << " from " << cfg.csv << std::endl;
    return;
  }

  TGraphErrors* logGraph = MakeLogGraph(pts);
  static std::vector<TObject*> keep;
  keep.push_back(logGraph);
  TF1 logFit(("fit_" + cfg.variant).c_str(),
             "[0] + ([1] + [2]*TMath::Log(x) + [3]*x)*TMath::Log(1.0/x)",
             12.0, cfg.xmax);
  logFit.SetParameters(22.0, 5.0, 0.0, 0.0);
  logGraph->Fit(&logFit, "QNR");

  const auto p1 = SelectRange(pts, 12.0, cfg.boundary1);
  const auto p2 = SelectRange(pts, cfg.boundary1, cfg.boundary2 > 0.0 ? cfg.boundary2 : cfg.xmax + 0.1);
  const auto p3 = cfg.boundary2 > 0.0 ? SelectRange(pts, cfg.boundary2, cfg.xmax + 0.1) : std::vector<Point>();

  TGraphErrors* g1 = MakeGraph(p1, "g1", cfg.color1, 20);
  keep.push_back(g1);
  TGraphErrors* g2 = MakeGraph(p2, "g2", cfg.color2, 21);
  keep.push_back(g2);
  TGraphErrors* g3 = cfg.boundary2 > 0.0 ? MakeGraph(p3, "g3", cfg.color3, 22) : nullptr;
  if (g3) keep.push_back(g3);
  TGraph* ref = MakeReference(&logFit, cfg.xmax);
  keep.push_back(ref);
  TGraphErrors* r1 = MakeRatio(p1, &logFit, "r1", cfg.color1, 20);
  keep.push_back(r1);
  TGraphErrors* r2 = MakeRatio(p2, &logFit, "r2", cfg.color2, 21);
  keep.push_back(r2);
  TGraphErrors* r3 = cfg.boundary2 > 0.0 ? MakeRatio(p3, &logFit, "r3", cfg.color3, 22) : nullptr;
  if (r3) keep.push_back(r3);

  jump1 = MeanRatio(pts, &logFit, cfg.boundary1 - 2.0, cfg.boundary1) > 0.0
            ? MeanRatio(pts, &logFit, cfg.boundary1, cfg.boundary1 + 2.0) /
                MeanRatio(pts, &logFit, cfg.boundary1 - 2.0, cfg.boundary1)
            : std::numeric_limits<double>::quiet_NaN();
  jump2 = std::numeric_limits<double>::quiet_NaN();
  if (cfg.boundary2 > 0.0)
  {
    jump2 = MeanRatio(pts, &logFit, cfg.boundary2 - 2.0, cfg.boundary2) > 0.0
              ? MeanRatio(pts, &logFit, cfg.boundary2, cfg.boundary2 + 2.0) /
                  MeanRatio(pts, &logFit, cfg.boundary2 - 2.0, cfg.boundary2)
              : std::numeric_limits<double>::quiet_NaN();
  }

  auto* top = new TPad((cfg.variant + "_top").c_str(), "", x1, y1 + 0.17 * (y2 - y1), x2, y2);
  auto* bot = new TPad((cfg.variant + "_bot").c_str(), "", x1, y1, x2, y1 + 0.16 * (y2 - y1));
  keep.push_back(top);
  keep.push_back(bot);
  top->SetLeftMargin(0.14);
  top->SetRightMargin(0.035);
  top->SetTopMargin(0.13);
  top->SetBottomMargin(0.015);
  top->SetLogy(true);
  bot->SetLeftMargin(0.14);
  bot->SetRightMargin(0.035);
  bot->SetTopMargin(0.02);
  bot->SetBottomMargin(0.34);
  top->Draw();
  bot->Draw();

  top->cd();
  auto* frame = new TH1F((cfg.variant + "_frame").c_str(), "", 100, 12.0, cfg.xmax);
  keep.push_back(frame);
  frame->SetMinimum(cfg.ymin);
  frame->SetMaximum(cfg.ymax);
  frame->GetYaxis()->SetTitle("#sigma_{eff} #times N scaled entries [pb / bin]");
  frame->GetYaxis()->SetTitleFont(132);
  frame->GetYaxis()->SetLabelFont(132);
  frame->GetYaxis()->SetTitleSize(0.058);
  frame->GetYaxis()->SetLabelSize(0.048);
  frame->GetYaxis()->SetTitleOffset(1.10);
  frame->GetXaxis()->SetLabelSize(0.0);
  frame->GetXaxis()->SetTitleSize(0.0);
  frame->SetStats(false);
  frame->Draw("AXIS");
  ref->Draw("L SAME");
  g1->Draw("PZ SAME");
  g2->Draw("PZ SAME");
  if (g3) g3->Draw("PZ SAME");

  auto* b1 = new TLine(cfg.boundary1, cfg.ymin, cfg.boundary1, cfg.ymax);
  keep.push_back(b1);
  b1->SetLineColor(kGray + 1);
  b1->SetLineStyle(3);
  b1->SetLineWidth(2);
  b1->Draw("SAME");
  if (cfg.boundary2 > 0.0)
  {
    auto* b2 = new TLine(cfg.boundary2, cfg.ymin, cfg.boundary2, cfg.ymax);
    keep.push_back(b2);
    b2->SetLineColor(kGray + 1);
    b2->SetLineStyle(3);
    b2->SetLineWidth(2);
    b2->Draw("SAME");
  }

  TLatex lat;
  lat.SetNDC(true);
  lat.SetTextFont(132);
  lat.SetTextSize(0.054);
  lat.DrawLatex(0.165, 0.870, cfg.title.c_str());
  lat.SetTextSize(0.038);
  lat.DrawLatex(0.165, 0.800, cfg.subtitle.c_str());
  lat.SetTextSize(0.030);
  lat.DrawLatex(0.165, 0.245, cfg.gate1.c_str());
  lat.DrawLatex(0.165, 0.198, cfg.gate2.c_str());
  if (!cfg.gate3.empty()) lat.DrawLatex(0.165, 0.151, cfg.gate3.c_str());

  auto* leg = new TLegend(0.59, 0.63, 0.94, 0.85);
  keep.push_back(leg);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(132);
  leg->SetTextSize(0.032);
  leg->AddEntry(g1, cfg.leg1.c_str(), "lep");
  leg->AddEntry(g2, cfg.leg2.c_str(), "lep");
  if (g3) leg->AddEntry(g3, cfg.leg3.c_str(), "lep");
  leg->AddEntry(ref, "Fit: modified power law", "l");
  leg->Draw();

  bot->cd();
  auto* rframe = new TH1F((cfg.variant + "_rframe").c_str(), "", 100, 12.0, cfg.xmax);
  keep.push_back(rframe);
  rframe->SetMinimum(0.84);
  rframe->SetMaximum(1.16);
  rframe->GetYaxis()->SetTitle("stitched / fit");
  rframe->GetXaxis()->SetTitle(cfg.xtitle.c_str());
  rframe->GetYaxis()->SetTitleFont(132);
  rframe->GetYaxis()->SetLabelFont(132);
  rframe->GetXaxis()->SetTitleFont(132);
  rframe->GetXaxis()->SetLabelFont(132);
  rframe->GetYaxis()->SetTitleSize(0.120);
  rframe->GetYaxis()->SetLabelSize(0.092);
  rframe->GetYaxis()->SetTitleOffset(0.48);
  rframe->GetYaxis()->SetNdivisions(505);
  rframe->GetXaxis()->SetTitleSize(0.120);
  rframe->GetXaxis()->SetLabelSize(0.092);
  rframe->GetXaxis()->SetTitleOffset(1.05);
  rframe->SetStats(false);
  rframe->Draw("AXIS");
  auto* one = new TLine(12.0, 1.0, cfg.xmax, 1.0);
  keep.push_back(one);
  one->SetLineColor(kGray + 2);
  one->SetLineStyle(2);
  one->Draw("SAME");
  r1->Draw("PZ SAME");
  r2->Draw("PZ SAME");
  if (r3) r3->Draw("PZ SAME");
  auto* rb1 = new TLine(cfg.boundary1, 0.84, cfg.boundary1, 1.16);
  keep.push_back(rb1);
  rb1->SetLineColor(kGray + 1);
  rb1->SetLineStyle(3);
  rb1->SetLineWidth(2);
  rb1->Draw("SAME");
  if (cfg.boundary2 > 0.0)
  {
    auto* rb2 = new TLine(cfg.boundary2, 0.84, cfg.boundary2, 1.16);
    keep.push_back(rb2);
    rb2->SetLineColor(kGray + 1);
    rb2->SetLineStyle(3);
    rb2->SetLineWidth(2);
    rb2->Draw("SAME");
  }
}
}

void MakeExpandedStitchCrossCheckBackupSlide()
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetErrorX(0.5);
  gStyle->SetEndErrorSize(2);

  gSystem->mkdir(kOutDir.c_str(), true);
  const std::string png = kOutDir + "/expanded_stitch_crosscheck_signal_bkg_1x2.png";
  const std::string json = kOutDir + "/expanded_stitch_crosscheck_signal_bkg_1x2_summary.json";

  TCanvas c("c_expanded_stitch_crosscheck", "Expanded Stitching Cross-Check", 1600, 900);
  c.SetWindowSize(1600 + (1600 - c.GetWw()), 900 + (900 - c.GetWh()));
  c.SetCanvasSize(1600, 900);
  c.SetFillColor(kWhite);
  c.SetMargin(0, 0, 0, 0);

  Text(0.035, 0.944, "Expanded Stitching Cross-Check", 0.054, kBlack, 11, 22);
  Text(0.035, 0.902,
       "Same corrected bounded ownership logic, shown over a wider generator-filter p_{T} range for backup QA.",
       0.026, kGray + 2, 11, 132);
  Text(0.965, 0.944, "#it{#bf{sPHENIX}} Internal", 0.030, kBlack, 31, 132);

  PanelConfig photon;
  photon.csv = kPhotonCsv;
  photon.variant = "bounded_12_21_ge21";
  photon.title = "Signal: PhotonJet12+20";
  photon.subtitle = "bounded 12-21 / #geq21";
  photon.xtitle = "p_{T}^{#gamma,filter} [GeV]";
  photon.xmax = 40.0;
  photon.rebinWidth = 1.0;
  photon.ymin = 1.5e3;
  photon.ymax = 8.0e7;
  photon.boundary1 = 21.0;
  photon.boundary2 = -1.0;
  photon.leg1 = "PhotonJet12 stitched";
  photon.leg2 = "PhotonJet20 stitched";
  photon.gate1 = "PhotonJet12: 12 #leq p_{T}^{#gamma,filter} < 21 GeV";
  photon.gate2 = "PhotonJet20: p_{T}^{#gamma,filter} #geq 21 GeV";
  photon.color1 = kBlue + 1;
  photon.color2 = kRed + 1;

  PanelConfig inclusive;
  inclusive.csv = kInclusiveCsv;
  inclusive.variant = "exclusive31";
  inclusive.title = "Background: inclusive Jet12+20+30";
  inclusive.subtitle = "bounded 12-21 / 21-31 / #geq31";
  inclusive.xtitle = "p_{T}^{jet,filter} [GeV]";
  inclusive.xmax = 40.0;
  inclusive.ymin = 1.0;
  inclusive.ymax = 1.2e6;
  inclusive.boundary1 = 21.0;
  inclusive.boundary2 = 31.0;
  inclusive.leg1 = "Jet12 stitched";
  inclusive.leg2 = "Jet20 stitched";
  inclusive.leg3 = "Jet30 stitched";
  inclusive.gate1 = "Jet12: 12 #leq p_{T}^{jet,filter} < 21 GeV";
  inclusive.gate2 = "Jet20: 21 #leq p_{T}^{jet,filter} < 31 GeV";
  inclusive.gate3 = "Jet30: p_{T}^{jet,filter} #geq 31 GeV";
  inclusive.color1 = kBlue + 1;
  inclusive.color2 = TColor::GetColor("#ff7f0e");
  inclusive.color3 = TColor::GetColor("#cc3399");

  double photonJump = 0.0;
  double photonUnused = 0.0;
  double inclJump21 = 0.0;
  double inclJump31 = 0.0;
  DrawPanel(photon, 0.055, 0.210, 0.485, 0.865, photonJump, photonUnused);
  c.cd();
  DrawPanel(inclusive, 0.535, 0.210, 0.965, 0.865, inclJump21, inclJump31);
  c.cd();

  c.cd();
  auto* band = new TPaveText(0.055, 0.060, 0.965, 0.160, "NDC");
  band->SetFillColor(TColor::GetColor("#f7f8fa"));
  band->SetLineColor(TColor::GetColor("#d9dde4"));
  band->SetBorderSize(1);
  band->SetShadowColor(0);
  band->Draw();
  Text(0.500, 0.118,
       "Backup readout: bounded ownership stays smooth at the sample handoffs; no lower-only overlap is used for these final-candidate spectra.",
       0.026, kBlack, 21, 132);
  Text(0.500, 0.080,
       Form("Boundary jumps: signal 21 GeV = %.3f; background 21 GeV = %.3f, 31 GeV = %.3f.",
            photonJump, inclJump21, inclJump31),
       0.022, kGray + 2, 21, 132);

  c.SaveAs(png.c_str());

  std::ofstream out(json);
  out << std::setprecision(12)
      << "{\n"
      << "  \"png\": \"" << png << "\",\n"
      << "  \"signal_csv\": \"" << kPhotonCsv << "\",\n"
      << "  \"signal_variant\": \"bounded_12_21_ge21\",\n"
      << "  \"signal_xmax_gev\": 40,\n"
      << "  \"signal_jump_21\": " << photonJump << ",\n"
      << "  \"background_csv\": \"" << kInclusiveCsv << "\",\n"
      << "  \"background_variant\": \"exclusive31\",\n"
      << "  \"background_xmax_gev\": 40,\n"
      << "  \"background_jump_21\": " << inclJump21 << ",\n"
      << "  \"background_jump_31\": " << inclJump31 << ",\n"
      << "  \"fit_function\": \"A*(1/x)^(b + c ln(x) + d x)\"\n"
      << "}\n";

  std::cout << "[DONE] " << png << std::endl;
  std::cout << "[SUMMARY] signal21=" << photonJump
            << " background21=" << inclJump21
            << " background31=" << inclJump31 << std::endl;
}
