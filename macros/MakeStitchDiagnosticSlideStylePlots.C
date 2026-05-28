#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TLatex.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace
{
struct Table
{
  std::map<std::string, std::vector<double>> cols;
};

std::vector<std::string> Split(const std::string& line)
{
  std::vector<std::string> out;
  std::stringstream ss(line);
  std::string item;
  while (std::getline(ss, item, ',')) out.push_back(item);
  return out;
}

Table ReadCsv(const std::string& path)
{
  std::ifstream in(path);
  if (!in)
  {
    std::cerr << "[ERROR] Cannot open " << path << std::endl;
    return {};
  }
  std::string headerLine;
  std::getline(in, headerLine);
  const auto headers = Split(headerLine);
  Table table;
  for (const auto& h : headers) table.cols[h] = {};

  std::string line;
  while (std::getline(in, line))
  {
    if (line.empty()) continue;
    const auto cells = Split(line);
    for (size_t i = 0; i < headers.size(); ++i)
    {
      double value = 0.0;
      if (i < cells.size() && !cells[i].empty()) value = std::stod(cells[i]);
      table.cols[headers[i]].push_back(value);
    }
  }
  return table;
}

bool HasCol(const Table& t, const std::string& key)
{
  return t.cols.find(key) != t.cols.end();
}

double Col(const Table& t, const std::string& key, size_t i)
{
  const auto it = t.cols.find(key);
  if (it == t.cols.end() || i >= it->second.size()) return 0.0;
  return it->second[i];
}

int NRows(const Table& t)
{
  const auto it = t.cols.find("pt_center");
  return it == t.cols.end() ? 0 : static_cast<int>(it->second.size());
}

struct Series
{
  std::string label;
  std::string yKey;
  std::string eyKey;
  int color = kBlack;
  int marker = 20;
};

std::unique_ptr<TGraphErrors> MakeGraph(
    const Table& table,
    const Series& series,
    double xMin,
    double xMax,
    bool ratio)
{
  auto graph = std::make_unique<TGraphErrors>();
  int point = 0;
  for (int i = 0; i < NRows(table); ++i)
  {
    const double x = Col(table, "pt_center", i);
    const double y = Col(table, series.yKey, i);
    const double ey = Col(table, series.eyKey, i);
    const double fit = Col(table, "fit_pb_per_GeV", i);
    if (x < xMin || x > xMax || y <= 0.0 || !std::isfinite(y)) continue;
    double yy = y;
    double eyy = ey;
    if (ratio)
    {
      if (fit <= 0.0 || !std::isfinite(fit)) continue;
      yy = y / fit;
      eyy = ey / fit;
    }
    graph->SetPoint(point, x, yy);
    graph->SetPointError(point, 0.0, eyy);
    ++point;
  }
  graph->SetLineColor(series.color);
  graph->SetMarkerColor(series.color);
  graph->SetMarkerStyle(series.marker);
  graph->SetMarkerSize(0.78);
  graph->SetLineWidth(1);
  return graph;
}

std::unique_ptr<TGraph> MakeFitGraph(const Table& table, double xMin, double xMax)
{
  auto graph = std::make_unique<TGraph>();
  int point = 0;
  for (int i = 0; i < NRows(table); ++i)
  {
    const double x = Col(table, "pt_center", i);
    const double y = Col(table, "fit_pb_per_GeV", i);
    if (x < xMin || x > xMax || y <= 0.0 || !std::isfinite(y)) continue;
    graph->SetPoint(point++, x, y);
  }
  graph->SetLineColor(kGray + 2);
  graph->SetLineStyle(2);
  graph->SetLineWidth(2);
  return graph;
}

double PositiveMin(const Table& table, const std::vector<Series>& series, double xMin, double xMax)
{
  double out = 1.0e300;
  for (int i = 0; i < NRows(table); ++i)
  {
    const double x = Col(table, "pt_center", i);
    if (x < xMin || x > xMax) continue;
    for (const auto& s : series)
    {
      const double y = Col(table, s.yKey, i);
      if (y > 0.0 && std::isfinite(y)) out = std::min(out, y);
    }
  }
  return out < 1.0e299 ? out : 1.0;
}

double PositiveMax(const Table& table, const std::vector<Series>& series, double xMin, double xMax)
{
  double out = 1.0;
  for (int i = 0; i < NRows(table); ++i)
  {
    const double x = Col(table, "pt_center", i);
    if (x < xMin || x > xMax) continue;
    for (const auto& s : series)
    {
      const double y = Col(table, s.yKey, i);
      if (y > 0.0 && std::isfinite(y)) out = std::max(out, y);
    }
  }
  return out;
}

void DrawOne(
    const std::string& csvPath,
    const std::string& outPath,
    const std::vector<Series>& series,
    const std::vector<double>& boundaries,
    const std::string& title,
    const std::string& ownership,
    const std::string& yTitle,
    const std::string& xTitle,
    const std::string& footnote,
    double xMax,
    double ratioMin,
    double ratioMax)
{
  const Table table = ReadCsv(csvPath);
  if (NRows(table) <= 0) return;

  const double xMin = 10.0;
  const double yMin = std::max(PositiveMin(table, series, xMin, xMax) * 0.35, 1.0e-3);
  const double yMax = PositiveMax(table, series, xMin, xMax) * 3.5;

  std::vector<std::unique_ptr<TGraphErrors>> topGraphs;
  std::vector<std::unique_ptr<TGraphErrors>> ratioGraphs;
  for (const auto& s : series)
  {
    topGraphs.push_back(MakeGraph(table, s, xMin, xMax, false));
    ratioGraphs.push_back(MakeGraph(table, s, xMin, xMax, true));
  }
  auto fit = MakeFitGraph(table, xMin, xMax);

  gSystem->mkdir(gSystem->DirName(outPath.c_str()), true);

  TCanvas canvas("c", "stitch", 1350, 900);
  canvas.SetFillColor(kWhite);
  TPad top("top", "", 0.0, 0.30, 1.0, 1.0);
  TPad bottom("bottom", "", 0.0, 0.0, 1.0, 0.30);
  top.SetLeftMargin(0.12);
  top.SetRightMargin(0.04);
  top.SetTopMargin(0.08);
  top.SetBottomMargin(0.02);
  top.SetLogy(true);
  top.SetTicks(1, 1);
  bottom.SetLeftMargin(0.12);
  bottom.SetRightMargin(0.04);
  bottom.SetTopMargin(0.02);
  bottom.SetBottomMargin(0.30);
  bottom.SetTicks(1, 1);
  top.Draw();
  bottom.Draw();

  top.cd();
  auto* frame = top.DrawFrame(xMin, yMin, xMax, yMax);
  frame->SetTitle("");
  frame->GetXaxis()->SetLabelSize(0.0);
  frame->GetXaxis()->SetTitleSize(0.0);
  frame->GetYaxis()->SetTitle(yTitle.c_str());
  frame->GetYaxis()->SetTitleSize(0.052);
  frame->GetYaxis()->SetLabelSize(0.045);
  frame->GetYaxis()->SetTitleOffset(1.08);

  for (const double boundary : boundaries)
  {
    TLine line(boundary, yMin, boundary, yMax);
    line.SetLineColor(kGray + 1);
    line.SetLineStyle(3);
    line.DrawClone("SAME");
  }
  fit->Draw("L SAME");
  for (auto& graph : topGraphs) graph->Draw("P SAME");

  TLatex text;
  text.SetNDC(true);
  text.SetTextFont(42);
  text.SetTextSize(0.032);
  text.DrawLatex(0.66, 0.52, "#it{#bf{sPHENIX}} Internal");

  TLegend leg(0.58, 0.58, 0.94, 0.88);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextFont(42);
  leg.SetTextSize(0.034);
  leg.AddEntry(fit.get(), "Fit: A(1+p_{T}/p_{0})^{-n}", "l");
  for (size_t i = 0; i < series.size(); ++i)
  {
    leg.AddEntry(topGraphs[i].get(), series[i].label.c_str(), "ep");
  }
  leg.Draw();

  bottom.cd();
  auto* rframe = bottom.DrawFrame(xMin, ratioMin, xMax, ratioMax);
  rframe->SetTitle("");
  rframe->GetXaxis()->SetTitle(xTitle.c_str());
  rframe->GetXaxis()->SetTitleSize(0.110);
  rframe->GetXaxis()->SetLabelSize(0.088);
  rframe->GetYaxis()->SetTitle("stitched / fit");
  rframe->GetYaxis()->SetTitleSize(0.090);
  rframe->GetYaxis()->SetLabelSize(0.078);
  rframe->GetYaxis()->SetTitleOffset(0.60);
  rframe->GetYaxis()->SetNdivisions(505);

  TLine one(xMin, 1.0, xMax, 1.0);
  one.SetLineColor(kGray + 2);
  one.SetLineStyle(2);
  one.Draw("SAME");
  for (const double boundary : boundaries)
  {
    TLine line(boundary, ratioMin, boundary, ratioMax);
    line.SetLineColor(kGray + 1);
    line.SetLineStyle(3);
    line.DrawClone("SAME");
  }
  for (auto& graph : ratioGraphs) graph->Draw("P SAME");

  canvas.cd();
  TLatex foot;
  foot.SetNDC(true);
  foot.SetTextFont(42);
  foot.SetTextSize(0.025);
  foot.SetTextColor(kGray + 2);
  foot.DrawLatex(0.12, 0.012, footnote.c_str());

  canvas.SaveAs(outPath.c_str());
  std::cout << "[WROTE] " << outPath << std::endl;
}
}  // namespace

void MakeStitchDiagnosticSlideStylePlots()
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(2);

  DrawOne(
      "dataOutput/stitchDiagnostics/blair21_stitch_20260521_0938/blair21_combined_spectrum.csv",
      "dataOutput/stitchDiagnostics/blair21_stitch_20260521_0938/blair21_corrected_embedded_inclusive_stitch_slide_style.png",
      {
          {"Jet12: 12 #leq p_{T} < 21", "run28_embeddedJet12_pb_per_GeV", "run28_embeddedJet12_err_pb_per_GeV", kBlue + 1, 20},
          {"Jet20: 21 #leq p_{T} < 30", "run28_embeddedJet20_pb_per_GeV", "run28_embeddedJet20_err_pb_per_GeV", kOrange + 7, 21},
          {"Jet30: p_{T} #geq 30", "run28_embeddedJet30_pb_per_GeV", "run28_embeddedJet30_err_pb_per_GeV", kMagenta + 1, 22},
      },
      {21.0, 30.0},
      "Embedded inclusive-jet generator stitching spectrum",
      "Corrected ownership: 12-21, 21-30, #geq30 GeV",
      "d#sigma/dp_{T}^{truth jet} [pb / GeV]",
      "Leading truth jet p_{T} [GeV]",
      "50 #times 1M generator shards; worker QA: no kept entries outside assigned ownership windows.",
      45.0,
      0.82,
      1.14);

  DrawOne(
      "dataOutput/stitchDiagnostics/photon12_20_stitch_20260521/embeddedPhoton12plus20_stitch_x35.csv",
      "dataOutput/stitchDiagnostics/photon12_20_stitch_20260521/embeddedPhoton12plus20_stitch_x35_slide_style.png",
      {
          {"PhotonJet12: 12 #leq p_{T} < 20", "photon12_pb_per_GeV", "photon12_err_pb_per_GeV", kBlue + 1, 20},
          {"PhotonJet20: p_{T} #geq 20", "photon20_pb_per_GeV", "photon20_err_pb_per_GeV", kRed + 1, 21},
      },
      {20.0},
      "Embedded photon generator stitching spectrum",
      "Current ownership: PhotonJet12 12-20, PhotonJet20 #geq20 GeV",
      "d#sigma/dp_{T}^{#gamma, filter} [pb / GeV]",
      "Leading truth photon p_{T} [GeV]",
      "Current 20 GeV handoff check; mean ratio jump from 18-20 to 20-22 GeV = 1.055.",
      35.0,
      0.88,
      1.10);
}
