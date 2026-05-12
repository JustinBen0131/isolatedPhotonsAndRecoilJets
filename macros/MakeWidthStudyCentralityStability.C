#include "sPhenixStyle.C"

#include <TBox.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <fstream>
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
  std::string cls;
  double mean = NAN;
};

struct AucRow
{
  std::string cent;
  std::string product;
  double auc = NAN;
  double gainPct = NAN;
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

std::map<std::string, FeatureRow> ReadFeatureRows(const std::string& path)
{
  std::ifstream in(path);
  std::map<std::string, FeatureRow> rows;
  std::string line;
  std::getline(in, line);
  while (std::getline(in, line))
  {
    const auto cols = SplitCsvLine(line);
    if (cols.size() < 5) continue;
    FeatureRow row;
    row.feature = cols[0];
    row.scope = cols[1];
    row.cls = cols[2];
    row.mean = std::stod(cols[4]);
    rows[row.feature + "|" + row.scope + "|" + row.cls] = row;
  }
  return rows;
}

std::map<std::string, AucRow> ReadAucRows(const std::string& path)
{
  std::ifstream in(path);
  std::map<std::string, AucRow> rows;
  std::string line;
  std::getline(in, line);
  while (std::getline(in, line))
  {
    const auto cols = SplitCsvLine(line);
    if (cols.size() < 6) continue;
    AucRow row;
    row.cent = cols[0];
    row.product = cols[2];
    row.auc = std::stod(cols[4]);
    row.gainPct = std::stod(cols[5]);
    rows[row.cent + "|" + row.product] = row;
  }
  return rows;
}

const FeatureRow& F(const std::map<std::string, FeatureRow>& rows,
                    const std::string& feature,
                    const std::string& scope,
                    const std::string& cls)
{
  return rows.at(feature + "|" + scope + "|" + cls);
}

const AucRow& A(const std::map<std::string, AucRow>& rows,
                const std::string& cent,
                const std::string& product)
{
  return rows.at(cent + "|" + product);
}

double DriftPct(const std::map<std::string, FeatureRow>& rows,
                const std::string& feature,
                const std::string& cls)
{
  const double central = F(rows, feature, "cent_0_20", cls).mean;
  const double peripheral = F(rows, feature, "cent_50_80", cls).mean;
  return 100.0 * (peripheral - central) / central;
}

void DrawSphenixLabel(const double x, const double y)
{
  TLatex text;
  text.SetNDC();
  text.SetTextAlign(11);
  text.SetTextFont(42);
  text.SetTextSize(0.032);
  text.DrawLatex(x, y, "#it{#bf{sPHENIX}} Internal");
  text.SetTextSize(0.024);
  text.DrawLatex(x, y - 0.042,
                 "Pythia overlay, #sqrt{s_{NN}} = 200 GeV; 15 #leq E_{T}^{cluster} < 30 GeV");
}

void DrawBar(TBox*& box,
             const double x,
             const double width,
             const double value,
             const int color)
{
  const double y0 = std::min(0.0, value);
  const double y1 = std::max(0.0, value);
  box = new TBox(x - width / 2.0, y0, x + width / 2.0, y1);
  box->SetFillColor(color);
  box->SetLineColor(color);
  box->Draw();
}

void DrawStabilitySummary(const std::string& outDir,
                          const std::map<std::string, FeatureRow>& features)
{
  struct Item
  {
    std::string feature;
    std::string label;
  };
  const std::vector<Item> items = {
      {"cluster_weta_cogx", "base #eta"},
      {"cluster_weta33_cogx", "3#times3 #eta"},
      {"cluster_wphi_cogx", "base #phi"},
      {"cluster_wphi33_cogx", "3#times3 #phi"}};

  TCanvas c("c_width_drift", "c_width_drift", 1240, 760);
  c.SetLeftMargin(0.13);
  c.SetRightMargin(0.04);
  c.SetTopMargin(0.18);
  c.SetBottomMargin(0.17);

  TH2D frame("frame_width_drift", ";Width variable;Mean shift from 0-20% to 50-80% [%]",
             4, 0, 4, 10, -36.0, 8.0);
  frame.SetStats(false);
  frame.GetXaxis()->SetTitleSize(0.044);
  frame.GetXaxis()->SetLabelSize(0.037);
  frame.GetYaxis()->SetTitleSize(0.044);
  frame.GetYaxis()->SetLabelSize(0.038);
  frame.GetYaxis()->SetTitleOffset(1.15);
  for (int i = 0; i < 4; ++i) frame.GetXaxis()->SetBinLabel(i + 1, items[i].label.c_str());
  frame.Draw("AXIS");

  TLine zero(0, 0, 4, 0);
  zero.SetLineColor(kGray + 2);
  zero.SetLineWidth(2);
  zero.Draw();

  const int signalColor = kAzure + 2;
  const int backgroundColor = kOrange + 7;
  const double barW = 0.24;
  TLatex value;
  value.SetTextFont(42);
  value.SetTextAlign(22);
  value.SetTextSize(0.028);
  for (int i = 0; i < 4; ++i)
  {
    const double x = i + 0.5;
    const double sig = DriftPct(features, items[i].feature, "signal");
    const double bkg = DriftPct(features, items[i].feature, "background");
    TBox* sigBox = nullptr;
    TBox* bkgBox = nullptr;
    DrawBar(sigBox, x - 0.15, barW, sig, signalColor);
    DrawBar(bkgBox, x + 0.15, barW, bkg, backgroundColor);
    value.SetTextColor(kBlack);
    value.DrawLatex(x - 0.15, sig - 2.6, Form("%.0f%%", sig));
    value.DrawLatex(x + 0.15, bkg - 2.6, Form("%.0f%%", bkg));
  }

  auto* leg = new TLegend(0.63, 0.75, 0.93, 0.84);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetFillStyle(1001);
  leg->SetTextFont(42);
  leg->SetTextSize(0.031);
  auto* sig = new TGraph();
  sig->SetFillColor(signalColor);
  auto* bkg = new TGraph();
  bkg->SetFillColor(backgroundColor);
  leg->AddEntry(sig, "truth-matched photons", "f");
  leg->AddEntry(bkg, "inclusive-jet background", "f");
  leg->Draw();

  TLatex title;
  title.SetNDC();
  title.SetTextFont(42);
  title.SetTextAlign(22);
  title.SetTextSize(0.039);
  title.DrawLatex(0.54, 0.955, "Centrality dependence of shower-width means");
  DrawSphenixLabel(0.16, 0.885);

  c.SaveAs((outDir + "/widthstudy_width_mean_centrality_drift.png").c_str());
}

void DrawAbsoluteStabilitySummary(const std::string& outDir,
                                  const std::map<std::string, FeatureRow>& features)
{
  struct Item
  {
    std::string feature;
    std::string label;
  };
  const std::vector<Item> items = {
      {"cluster_weta_cogx", "base #eta"},
      {"cluster_weta33_cogx", "3#times3 #eta"},
      {"cluster_wphi_cogx", "base #phi"},
      {"cluster_wphi33_cogx", "3#times3 #phi"}};

  TCanvas c("c_width_abs_drift", "c_width_abs_drift", 1240, 760);
  c.SetLeftMargin(0.13);
  c.SetRightMargin(0.04);
  c.SetTopMargin(0.18);
  c.SetBottomMargin(0.17);

  TH2D frame("frame_width_abs_drift", ";Width variable;Central-to-peripheral mean shift [%]",
             4, 0, 4, 10, 0.0, 42.0);
  frame.SetStats(false);
  frame.GetXaxis()->SetTitleSize(0.044);
  frame.GetXaxis()->SetLabelSize(0.037);
  frame.GetYaxis()->SetTitleSize(0.044);
  frame.GetYaxis()->SetLabelSize(0.038);
  frame.GetYaxis()->SetTitleOffset(1.08);
  for (int i = 0; i < 4; ++i) frame.GetXaxis()->SetBinLabel(i + 1, items[i].label.c_str());
  frame.Draw("AXIS");

  const int signalColor = kAzure + 2;
  const int backgroundColor = kOrange + 7;
  const double barW = 0.24;
  TLatex value;
  value.SetTextFont(42);
  value.SetTextAlign(22);
  value.SetTextSize(0.029);
  for (int i = 0; i < 4; ++i)
  {
    const double x = i + 0.5;
    const double sig = std::abs(DriftPct(features, items[i].feature, "signal"));
    const double bkg = std::abs(DriftPct(features, items[i].feature, "background"));
    TBox* sigBox = nullptr;
    TBox* bkgBox = nullptr;
    DrawBar(sigBox, x - 0.15, barW, sig, signalColor);
    DrawBar(bkgBox, x + 0.15, barW, bkg, backgroundColor);
    value.DrawLatex(x - 0.15, sig + 1.4, Form("%.0f%%", sig));
    value.DrawLatex(x + 0.15, bkg + 1.4, Form("%.0f%%", bkg));
  }

  auto* leg = new TLegend(0.62, 0.845, 0.93, 0.925);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.029);
  auto* sig = new TGraph();
  sig->SetFillColor(signalColor);
  auto* bkg = new TGraph();
  bkg->SetFillColor(backgroundColor);
  leg->AddEntry(sig, "truth-matched photons", "f");
  leg->AddEntry(bkg, "jet-background candidates", "f");
  leg->Draw();

  TLatex title;
  title.SetNDC();
  title.SetTextFont(42);
  title.SetTextAlign(22);
  title.SetTextSize(0.039);
  title.DrawLatex(0.54, 0.955, "Central-to-peripheral drift of shower-width means");
  DrawSphenixLabel(0.16, 0.885);

  c.SaveAs((outDir + "/widthstudy_width_mean_centrality_stability_abs.png").c_str());
}

void DrawStabilityAndAuc(const std::string& outDir,
                         const std::map<std::string, FeatureRow>& features,
                         const std::map<std::string, AucRow>& aucs)
{
  TCanvas c("c_stability_auc", "c_stability_auc", 1280, 760);
  c.Divide(2, 1, 0.02, 0.0);

  c.cd(1);
  gPad->SetLeftMargin(0.17);
  gPad->SetRightMargin(0.04);
  gPad->SetTopMargin(0.20);
  gPad->SetBottomMargin(0.22);
  TH2D left("left", ";Width type;Absolute mean shift [%]",
            2, 0, 2, 10, 0.0, 34.0);
  left.SetStats(false);
  left.GetXaxis()->SetTitleSize(0.049);
  left.GetXaxis()->SetLabelSize(0.043);
  left.GetYaxis()->SetTitleSize(0.049);
  left.GetYaxis()->SetLabelSize(0.041);
  left.GetYaxis()->SetTitleOffset(1.35);
  left.GetXaxis()->SetBinLabel(1, "base widths");
  left.GetXaxis()->SetBinLabel(2, "3#times3 widths");
  left.Draw("AXIS");
  auto meanAbsDrift = [&](const std::string& etaFeature, const std::string& phiFeature, const std::string& cls) {
    return 0.5 * (std::abs(DriftPct(features, etaFeature, cls)) + std::abs(DriftPct(features, phiFeature, cls)));
  };
  const double sigBase = meanAbsDrift("cluster_weta_cogx", "cluster_wphi_cogx", "signal");
  const double bkgBase = meanAbsDrift("cluster_weta_cogx", "cluster_wphi_cogx", "background");
  const double sig33 = meanAbsDrift("cluster_weta33_cogx", "cluster_wphi33_cogx", "signal");
  const double bkg33 = meanAbsDrift("cluster_weta33_cogx", "cluster_wphi33_cogx", "background");
  const int signalColor = kAzure + 2;
  const int backgroundColor = kOrange + 7;
  TBox* box = nullptr;
  DrawBar(box, 0.38, 0.22, sigBase, signalColor);
  DrawBar(box, 0.62, 0.22, bkgBase, backgroundColor);
  DrawBar(box, 1.38, 0.22, sig33, signalColor);
  DrawBar(box, 1.62, 0.22, bkg33, backgroundColor);
  TLatex lab;
  lab.SetTextFont(42);
  lab.SetTextAlign(22);
  lab.SetTextSize(0.030);
  for (const auto& point : std::vector<std::pair<double, double>>{{0.38, sigBase}, {0.62, bkgBase}, {1.38, sig33}, {1.62, bkg33}})
  {
    lab.DrawLatex(point.first, point.second + 1.3, Form("%.0f%%", point.second));
  }
  auto* leg = new TLegend(0.26, 0.76, 0.88, 0.86);
  leg->SetBorderSize(0);
  leg->SetFillColor(kWhite);
  leg->SetFillStyle(1001);
  leg->SetTextFont(42);
  leg->SetTextSize(0.032);
  auto* sigGraph = new TGraph();
  sigGraph->SetFillColor(signalColor);
  auto* bkgGraph = new TGraph();
  bkgGraph->SetFillColor(backgroundColor);
  leg->AddEntry(sigGraph, "truth-matched photons", "f");
  leg->AddEntry(bkgGraph, "inclusive-jet background", "f");
  leg->Draw();
  TLatex panel;
  panel.SetNDC();
  panel.SetTextFont(42);
  panel.SetTextAlign(22);
  panel.SetTextSize(0.043);
  panel.DrawLatex(0.55, 0.94, "Centrality stability of width means");

  c.cd(2);
  gPad->SetLeftMargin(0.16);
  gPad->SetRightMargin(0.05);
  gPad->SetTopMargin(0.20);
  gPad->SetBottomMargin(0.22);
  TH2D right("right", ";Centrality;ROC-AUC gain [%]",
             3, 0, 3, 10, 0.0, 1.25);
  right.SetStats(false);
  right.GetXaxis()->SetTitleSize(0.049);
  right.GetXaxis()->SetLabelSize(0.043);
  right.GetYaxis()->SetTitleSize(0.049);
  right.GetYaxis()->SetLabelSize(0.041);
  right.GetYaxis()->SetTitleOffset(1.18);
  right.GetXaxis()->SetBinLabel(1, "0-20%");
  right.GetXaxis()->SetBinLabel(2, "20-50%");
  right.GetXaxis()->SetBinLabel(3, "50-80%");
  right.Draw("AXIS");

  const std::vector<std::string> cents = {"0_20", "20_50", "50_80"};
  const int bothColor = kGreen + 2;
  lab.SetTextSize(0.031);
  for (int i = 0; i < 3; ++i)
  {
    const double gain = A(aucs, cents[i], "centAsFeatBase3x3_pt15to30").gainPct;
    TBox* gainBox = nullptr;
    DrawBar(gainBox, i + 0.5, 0.42, gain, bothColor);
    lab.DrawLatex(i + 0.5, gain + 0.08, Form("+%.2f%%", gain));
  }
  panel.DrawLatex(0.54, 0.94, "ROC-AUC gain from adding 3#times3 widths");

  c.cd();
  TLatex head;
  head.SetNDC();
  head.SetTextFont(42);
  head.SetTextAlign(11);
  head.SetTextSize(0.026);
  head.DrawLatex(0.06, 0.965, "#it{#bf{sPHENIX}} Internal");
  head.SetTextSize(0.020);
  head.DrawLatex(0.06, 0.935,
                 "Pythia overlay, #sqrt{s_{NN}} = 200 GeV; 15 #leq E_{T}^{cluster} < 30 GeV");

  c.SaveAs((outDir + "/widthstudy_centrality_stability_and_auc_gain.png").c_str());
}
}

void MakeWidthStudyCentralityStability()
{
  SetsPhenixStyle();
  gStyle->SetOptStat(0);
  const std::string outDir = "dataOutput/jstg_slide_candidates/slide17_widthstudy_centrality_gain";
  gSystem->mkdir(outDir.c_str(), true);
  const auto features = ReadFeatureRows("dataOutput/auauTightBDTValidation/model_validation_condor_20260511_110255/validation_feature_summary.csv");
  const auto aucs = ReadAucRows(outDir + "/widthstudy_centrality_auc_summary.csv");
  DrawStabilitySummary(outDir, features);
  DrawAbsoluteStabilitySummary(outDir, features);
  DrawStabilityAndAuc(outDir, features, aucs);
  std::cout << "[MakeWidthStudyCentralityStability] wrote " << outDir << "\n";
}
