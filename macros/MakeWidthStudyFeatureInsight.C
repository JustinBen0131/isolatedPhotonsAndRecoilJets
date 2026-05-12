#include "sPhenixStyle.C"

#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace
{
struct FeatureRow
{
  std::string feature;
  std::string scope;
  std::string klass;
  int entries = 0;
  double mean = NAN;
  double std = NAN;
};

std::vector<std::string> SplitCsvLine(const std::string& line)
{
  std::vector<std::string> out;
  std::string cur;
  bool quoted = false;
  for (const char c : line)
  {
    if (c == '"')
    {
      quoted = !quoted;
      continue;
    }
    if (c == ',' && !quoted)
    {
      out.push_back(cur);
      cur.clear();
      continue;
    }
    cur.push_back(c);
  }
  out.push_back(cur);
  return out;
}

std::map<std::string, FeatureRow> ReadFeatureSummary(const std::string& path)
{
  std::ifstream in(path);
  std::map<std::string, FeatureRow> rows;
  std::string line;
  std::getline(in, line);
  while (std::getline(in, line))
  {
    const auto cols = SplitCsvLine(line);
    if (cols.size() < 6) continue;
    FeatureRow r;
    r.feature = cols[0];
    r.scope = cols[1];
    r.klass = cols[2];
    r.entries = std::stoi(cols[3]);
    r.mean = std::stod(cols[4]);
    r.std = std::stod(cols[5]);
    rows[r.feature + "|" + r.scope + "|" + r.klass] = r;
  }
  return rows;
}

const FeatureRow& Row(const std::map<std::string, FeatureRow>& rows,
                      const std::string& feature,
                      const std::string& scope,
                      const std::string& klass)
{
  return rows.at(feature + "|" + scope + "|" + klass);
}

double CohenD(const FeatureRow& signal, const FeatureRow& background)
{
  const double pooled = std::sqrt(0.5 * (signal.std * signal.std + background.std * background.std));
  if (pooled <= 0) return NAN;
  return std::fabs(background.mean - signal.mean) / pooled;
}

double CentralityVariationPercent(const std::map<std::string, FeatureRow>& rows,
                                  const std::string& feature,
                                  const std::string& klass)
{
  const std::vector<std::string> scopes = {"cent_0_20", "cent_20_50", "cent_50_80"};
  std::vector<double> means;
  for (const auto& scope : scopes)
  {
    means.push_back(Row(rows, feature, scope, klass).mean);
  }
  const double inclusive = Row(rows, feature, "inclusive", klass).mean;
  if (inclusive == 0) return NAN;
  const auto [mn, mx] = std::minmax_element(means.begin(), means.end());
  return (*mx - *mn) / inclusive * 100.0;
}

void DrawSphenixLabel(const double x, const double y)
{
  TLatex text;
  text.SetNDC();
  text.SetTextAlign(11);
  text.SetTextFont(42);
  text.SetTextSize(0.033);
  text.DrawLatex(x, y, "#it{#bf{sPHENIX}} Internal");
  text.SetTextSize(0.025);
  text.DrawLatex(x, y - 0.043,
                 "Pythia overlay, #sqrt{s_{NN}} = 200 GeV; 15 #leq E_{T}^{cluster} < 30 GeV");
}

void DrawComplementarityPlot(const std::string& outDir,
                             const std::map<std::string, FeatureRow>& rows)
{
  const std::vector<std::string> features = {
    "cluster_weta_cogx",
    "cluster_weta33_cogx",
    "cluster_wphi_cogx",
    "cluster_wphi33_cogx",
  };
  const std::vector<std::string> labels = {
    "Full #eta width",
    "Local 3#times3 #eta width",
    "Full #phi width",
    "Local 3#times3 #phi width",
  };
  const int colors[] = {kBlue + 1, kAzure + 7, kOrange + 7, kRed + 1};

  TCanvas c("c_feature_complementarity", "c_feature_complementarity", 1500, 760);
  c.SetTopMargin(0.0);
  c.SetBottomMargin(0.0);
  c.SetLeftMargin(0.0);
  c.SetRightMargin(0.0);

  TLatex header;
  header.SetNDC();
  header.SetTextFont(42);
  header.SetTextAlign(22);
  header.SetTextSize(0.038);
  header.DrawLatex(0.53, 0.955, "Why do full-cluster and local 3#times3 widths work best together?");
  header.SetTextAlign(11);
  header.SetTextSize(0.028);
  header.DrawLatex(0.075, 0.925, "#it{#bf{sPHENIX}} Internal");
  header.SetTextSize(0.021);
  header.DrawLatex(0.075, 0.892,
                   "Pythia overlay, #sqrt{s_{NN}} = 200 GeV; 15 #leq E_{T}^{cluster} < 30 GeV");

  TPad leftPad("leftPad", "leftPad", 0.04, 0.08, 0.51, 0.82);
  TPad rightPad("rightPad", "rightPad", 0.53, 0.08, 0.98, 0.82);
  leftPad.Draw();
  rightPad.Draw();

  leftPad.cd();
  gPad->SetLeftMargin(0.20);
  gPad->SetRightMargin(0.04);
  gPad->SetTopMargin(0.10);
  gPad->SetBottomMargin(0.14);
  TH2D frame1("frame1", ";Signal/background separation  |Cohen d|;Shower-width variable",
              10, 0, 1.15, 4, 0, 4);
  frame1.SetStats(false);
  frame1.GetXaxis()->SetTitleSize(0.045);
  frame1.GetXaxis()->SetLabelSize(0.038);
  frame1.GetYaxis()->SetTitleSize(0.045);
  frame1.GetYaxis()->SetLabelSize(0.034);
  frame1.GetYaxis()->SetTitleOffset(2.10);
  for (int i = 0; i < 4; ++i)
  {
    frame1.GetYaxis()->SetBinLabel(4 - i, labels[i].c_str());
  }
  frame1.Draw("AXIS");

  TLatex text;
  text.SetTextFont(42);
  text.SetTextAlign(12);
  for (int i = 0; i < 4; ++i)
  {
    const double d = CohenD(Row(rows, features[i], "inclusive", "signal"),
                            Row(rows, features[i], "inclusive", "background"));
    const double y = 3.5 - i;
    auto* bar = new TLine(0, y, d, y);
    bar->SetLineWidth(16);
    bar->SetLineColor(colors[i]);
    bar->Draw();
    auto* marker = new TGraphErrors(1);
    marker->SetPoint(0, d, y);
    marker->SetMarkerStyle(20);
    marker->SetMarkerSize(1.7);
    marker->SetMarkerColor(colors[i]);
    marker->SetLineColor(colors[i]);
    marker->Draw("P");
    text.SetTextSize(0.033);
    text.DrawLatex(d + 0.035, y, Form("%.2f", d));
  }
  text.SetNDC();
  text.SetTextAlign(22);
  text.SetTextSize(0.043);
  text.DrawLatex(0.58, 0.955, "Signal/background separation");

  rightPad.cd();
  gPad->SetLeftMargin(0.20);
  gPad->SetRightMargin(0.04);
  gPad->SetTopMargin(0.10);
  gPad->SetBottomMargin(0.14);
  TH2D frame2("frame2", ";Centrality variation of mean [%];Shower-width variable",
              10, 0, 46, 4, 0, 4);
  frame2.SetStats(false);
  frame2.GetXaxis()->SetTitleSize(0.045);
  frame2.GetXaxis()->SetLabelSize(0.038);
  frame2.GetYaxis()->SetTitleSize(0.045);
  frame2.GetYaxis()->SetLabelSize(0.034);
  frame2.GetYaxis()->SetTitleOffset(2.10);
  for (int i = 0; i < 4; ++i)
  {
    frame2.GetYaxis()->SetBinLabel(4 - i, labels[i].c_str());
  }
  frame2.Draw("AXIS");

  for (int i = 0; i < 4; ++i)
  {
    const double sigVar = CentralityVariationPercent(rows, features[i], "signal");
    const double bkgVar = CentralityVariationPercent(rows, features[i], "background");
    const double y = 3.5 - i;
    auto* sig = new TGraphErrors(1);
    sig->SetPoint(0, sigVar, y + 0.10);
    sig->SetMarkerStyle(20);
    sig->SetMarkerSize(1.6);
    sig->SetMarkerColor(colors[i]);
    sig->SetLineColor(colors[i]);
    sig->Draw("P");
    auto* bkg = new TGraphErrors(1);
    bkg->SetPoint(0, bkgVar, y - 0.10);
    bkg->SetMarkerStyle(24);
    bkg->SetMarkerSize(1.6);
    bkg->SetMarkerColor(colors[i]);
    bkg->SetLineColor(colors[i]);
    bkg->Draw("P");
    text.SetNDC(false);
    text.SetTextAlign(12);
    text.SetTextSize(0.030);
    text.DrawLatex(sigVar + 1.0, y + 0.10, Form("%.1f", sigVar));
    text.DrawLatex(bkgVar + 1.0, y - 0.10, Form("%.1f", bkgVar));
  }
  text.SetNDC();
  text.SetTextAlign(22);
  text.SetTextSize(0.040);
  text.DrawLatex(0.58, 0.955, "Centrality stability  (filled: signal, open: background)");

  c.SaveAs((outDir + "/widthstudy_feature_complementarity_2panel.png").c_str());
}

void DrawCentralityTrendPlot(const std::string& outDir,
                             const std::map<std::string, FeatureRow>& rows)
{
  const std::vector<std::string> features = {
    "cluster_weta_cogx",
    "cluster_weta33_cogx",
    "cluster_wphi_cogx",
    "cluster_wphi33_cogx",
  };
  const std::vector<std::string> labels = {
    "Full #eta width",
    "Local 3#times3 #eta width",
    "Full #phi width",
    "Local 3#times3 #phi width",
  };
  const int colors[] = {kBlue + 1, kAzure + 7, kOrange + 7, kRed + 1};
  const std::vector<std::string> scopes = {"cent_0_20", "cent_20_50", "cent_50_80"};
  const double x[3] = {10, 35, 65};
  const double ex[3] = {0, 0, 0};

  TCanvas c("c_width_centrality_trends", "c_width_centrality_trends", 1420, 760);
  c.SetLeftMargin(0.12);
  c.SetRightMargin(0.05);
  c.SetTopMargin(0.18);
  c.SetBottomMargin(0.14);

  TH2D frame("frame", ";Centrality percentile;Mean shower-width value",
             10, 0, 80, 10, 0.09, 0.39);
  frame.SetStats(false);
  frame.GetXaxis()->SetTitleSize(0.045);
  frame.GetXaxis()->SetLabelSize(0.038);
  frame.GetYaxis()->SetTitleSize(0.045);
  frame.GetYaxis()->SetLabelSize(0.038);
  frame.Draw("AXIS");

  TLatex text;
  text.SetTextFont(42);
  text.SetTextSize(0.030);
  text.SetTextAlign(22);
  text.DrawLatex(10, 0.382, "0-20%");
  text.DrawLatex(35, 0.382, "20-50%");
  text.DrawLatex(65, 0.382, "50-80%");
  for (double v : {20.0, 50.0})
  {
    TLine line(v, 0.09, v, 0.375);
    line.SetLineStyle(2);
    line.SetLineColor(kGray + 1);
    line.Draw();
  }

  auto* leg = new TLegend(0.56, 0.66, 0.92, 0.86);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.031);

  for (int i = 0; i < 4; ++i)
  {
    double y[3];
    double ey[3];
    for (int j = 0; j < 3; ++j)
    {
      const auto& r = Row(rows, features[i], scopes[j], "signal");
      y[j] = r.mean;
      ey[j] = r.std / std::sqrt(std::max(1, r.entries));
    }
    auto* g = new TGraphErrors(3, x, y, ex, ey);
    g->SetLineColor(colors[i]);
    g->SetMarkerColor(colors[i]);
    g->SetLineWidth(3);
    g->SetMarkerStyle(i < 2 ? 20 : 21);
    g->SetMarkerSize(1.5);
    g->Draw("PL");
    leg->AddEntry(g, labels[i].c_str(), "pl");
  }
  leg->Draw();

  text.SetNDC();
  text.SetTextAlign(22);
  text.SetTextSize(0.043);
  text.DrawLatex(0.53, 0.945, "Local 3#times3 widths are more stable across centrality");
  DrawSphenixLabel(0.15, 0.875);

  c.SaveAs((outDir + "/widthstudy_width_means_vs_centrality.png").c_str());
}
}

void MakeWidthStudyFeatureInsight()
{
  SetsPhenixStyle();
  gStyle->SetOptStat(0);
  const std::string input = "dataOutput/auauTightBDTValidation/model_validation_condor_20260511_110255/validation_feature_summary.csv";
  const std::string outDir = "dataOutput/jstg_slide_candidates/slide16_widthstudy_feature_insight";
  gSystem->mkdir(outDir.c_str(), true);

  const auto rows = ReadFeatureSummary(input);
  DrawComplementarityPlot(outDir, rows);
  DrawCentralityTrendPlot(outDir, rows);
  std::cout << "[MakeWidthStudyFeatureInsight] wrote " << outDir << "\n";
}
