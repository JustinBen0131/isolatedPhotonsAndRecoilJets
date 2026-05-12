#include "sPhenixStyle.C"

#include <TBox.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TGraph.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>

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
struct Importance
{
  std::string product;
  std::string modelLabel;
  std::string feature;
  std::string featureLabel;
  double gainFraction = 0.0;
  double splitFraction = 0.0;
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

std::map<std::string, Importance> ReadImportance(const std::string& path)
{
  std::ifstream in(path);
  std::map<std::string, Importance> rows;
  std::string line;
  std::getline(in, line);
  while (std::getline(in, line))
  {
    const auto cols = SplitCsvLine(line);
    if (cols.size() < 9) continue;
    Importance r;
    r.product = cols[0];
    r.modelLabel = cols[1];
    r.feature = cols[2];
    r.featureLabel = cols[3];
    r.gainFraction = std::stod(cols[5]);
    r.splitFraction = std::stod(cols[7]);
    rows[r.product + "|" + r.feature] = r;
  }
  return rows;
}

double GainPercent(const std::map<std::string, Importance>& rows,
                   const std::string& product,
                   const std::string& feature)
{
  const auto it = rows.find(product + "|" + feature);
  if (it == rows.end()) return NAN;
  return 100.0 * it->second.gainFraction;
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

void DrawFeatureGainHeatmap(const std::string& outDir,
                            const std::map<std::string, Importance>& rows)
{
  const std::vector<std::string> products = {
    "centAsFeat_pt15to30",
    "centAsFeat3x3_pt15to30",
    "centAsFeatBase3x3_pt15to30",
  };
  const std::vector<std::string> productLabels = {
    "Base widths",
    "3#times3 widths",
    "Base + 3#times3",
  };
  const std::vector<std::string> features = {
    "e32_over_e35",
    "cluster_et1",
    "cluster_weta_cogx",
    "cluster_wphi_cogx",
    "cluster_weta33_cogx",
    "cluster_wphi33_cogx",
    "centrality",
    "cluster_Eta",
    "e11_over_e33",
    "cluster_Et",
  };
  const std::vector<std::string> featureLabels = {
    "E_{3#times2}/E_{3#times5}",
    "Leading tower E_{T}",
    "Full #eta width",
    "Full #phi width",
    "Local 3#times3 #eta width",
    "Local 3#times3 #phi width",
    "Centrality",
    "Cluster #eta",
    "E_{1#times1}/E_{3#times3}",
    "Cluster E_{T}",
  };

  const Int_t nStops = 5;
  Double_t stops[nStops] = {0.00, 0.20, 0.45, 0.72, 1.00};
  Double_t red[nStops] = {0.98, 0.89, 0.42, 0.10, 0.02};
  Double_t green[nStops] = {0.98, 0.92, 0.76, 0.45, 0.19};
  Double_t blue[nStops] = {0.98, 0.82, 0.62, 0.70, 0.45};
  TColor::CreateGradientColorTable(nStops, stops, red, green, blue, 255);
  gStyle->SetNumberContours(255);

  TCanvas c("c_feature_gain_heatmap", "c_feature_gain_heatmap", 1180, 900);
  c.SetLeftMargin(0.28);
  c.SetRightMargin(0.17);
  c.SetTopMargin(0.18);
  c.SetBottomMargin(0.17);

  TH2D h("h_gain", ";BDT trained with;Feature used by the BDT",
         3, 0, 3, features.size(), 0, features.size());
  h.SetStats(false);
  h.SetMinimum(0.0);
  h.SetMaximum(34.0);
  h.GetXaxis()->SetTitleSize(0.042);
  h.GetXaxis()->SetTitleOffset(1.28);
  h.GetXaxis()->SetLabelSize(0.036);
  h.GetYaxis()->SetTitleSize(0.042);
  h.GetYaxis()->SetTitleOffset(3.25);
  h.GetYaxis()->SetLabelSize(0.030);
  h.GetZaxis()->SetTitle("Fraction of split gain [%]");
  h.GetZaxis()->SetTitleSize(0.036);
  h.GetZaxis()->SetTitleOffset(1.20);
  h.GetZaxis()->SetLabelSize(0.032);

  for (int ix = 0; ix < 3; ++ix)
  {
    h.GetXaxis()->SetBinLabel(ix + 1, productLabels[ix].c_str());
  }
  for (size_t iy = 0; iy < features.size(); ++iy)
  {
    h.GetYaxis()->SetBinLabel(features.size() - iy, featureLabels[iy].c_str());
    for (int ix = 0; ix < 3; ++ix)
    {
      const double v = GainPercent(rows, products[ix], features[iy]);
      h.SetBinContent(ix + 1, features.size() - iy, std::isnan(v) ? 0.0 : v);
    }
  }
  h.Draw("COLZ");

  TLatex text;
  text.SetTextFont(42);
  text.SetTextAlign(22);
  for (size_t iy = 0; iy < features.size(); ++iy)
  {
    for (int ix = 0; ix < 3; ++ix)
    {
      const double v = GainPercent(rows, products[ix], features[iy]);
      const double x = ix + 0.5;
      const double y = features.size() - iy - 0.5;
      if (std::isnan(v))
      {
        text.SetTextSize(0.028);
        text.SetTextColor(kGray + 2);
        text.DrawLatex(x, y, "not used");
      }
      else
      {
        text.SetTextSize(0.031);
        text.SetTextColor(v > 19 ? kWhite : kBlack);
        text.DrawLatex(x, y, Form("%.1f", v));
      }
    }
  }
  text.SetTextColor(kBlack);
  for (int ix = 1; ix < 3; ++ix)
  {
    TLine line(ix, 0, ix, features.size());
    line.SetLineColor(kWhite);
    line.SetLineWidth(4);
    line.Draw();
  }

  TLatex ndc;
  ndc.SetNDC();
  ndc.SetTextFont(42);
  ndc.SetTextAlign(22);
  ndc.SetTextSize(0.039);
  ndc.DrawLatex(0.56, 0.955, "Which inputs drive the BDT score in each width study?");
  DrawSphenixLabel(0.31, 0.900);

  c.SaveAs((outDir + "/widthstudy_feature_gain_heatmap.png").c_str());
}

void DrawWidthOnlyBars(const std::string& outDir,
                       const std::map<std::string, Importance>& rows)
{
  const std::vector<std::string> products = {
    "centAsFeat_pt15to30",
    "centAsFeat3x3_pt15to30",
    "centAsFeatBase3x3_pt15to30",
  };
  const std::vector<std::string> productLabels = {
    "Base widths",
    "3#times3 widths",
    "Base + 3#times3",
  };
  const std::vector<int> colors = {kBlue + 1, kOrange + 7, kGreen + 2};
  const std::vector<std::string> features = {
    "cluster_weta_cogx",
    "cluster_wphi_cogx",
    "cluster_weta33_cogx",
    "cluster_wphi33_cogx",
  };
  const std::vector<std::string> labels = {
    "Full #eta width",
    "Full #phi width",
    "Local 3#times3 #eta width",
    "Local 3#times3 #phi width",
  };

  TCanvas c("c_width_gain_bars", "c_width_gain_bars", 1250, 760);
  c.SetLeftMargin(0.14);
  c.SetRightMargin(0.04);
  c.SetTopMargin(0.18);
  c.SetBottomMargin(0.20);

  TH2D frame("frame", ";Width feature;Fraction of split gain [%]",
             4, 0, 4, 10, 0, 15.0);
  frame.SetStats(false);
  frame.GetXaxis()->SetTitleSize(0.044);
  frame.GetXaxis()->SetTitleOffset(1.24);
  frame.GetXaxis()->SetLabelSize(0.034);
  frame.GetYaxis()->SetTitleSize(0.044);
  frame.GetYaxis()->SetLabelSize(0.038);
  for (int i = 0; i < 4; ++i)
  {
    frame.GetXaxis()->SetBinLabel(i + 1, labels[i].c_str());
  }
  frame.Draw("AXIS");

  auto* leg = new TLegend(0.58, 0.20, 0.93, 0.36);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.032);

  TLatex text;
  text.SetTextFont(42);
  text.SetTextAlign(22);
  for (int ip = 0; ip < 3; ++ip)
  {
    auto* g = new TGraph();
    g->SetMarkerStyle(20 + ip);
    g->SetMarkerSize(1.8);
    g->SetMarkerColor(colors[ip]);
    g->SetLineColor(colors[ip]);
    for (int ifeature = 0; ifeature < 4; ++ifeature)
    {
      const double x = ifeature + 0.5 + 0.16 * (ip - 1);
      const double v = GainPercent(rows, products[ip], features[ifeature]);
      g->SetPoint(ifeature, x, std::isnan(v) ? 0.0 : v);
      if (!std::isnan(v) && v > 0.05)
      {
        text.SetTextSize(0.027);
        text.DrawLatex(x, v + 0.55, Form("%.1f", v));
      }
    }
    g->Draw("P");
    leg->AddEntry(g, productLabels[ip].c_str(), "p");
  }
  leg->Draw();

  TLatex ndc;
  ndc.SetNDC();
  ndc.SetTextFont(42);
  ndc.SetTextAlign(22);
  ndc.SetTextSize(0.040);
  ndc.DrawLatex(0.54, 0.955, "Width-feature importance from XGBoost split gain");
  DrawSphenixLabel(0.17, 0.885);

  c.SaveAs((outDir + "/widthstudy_width_feature_gain_bars.png").c_str());
}

void DrawWidthContributionStack(const std::string& outDir,
                                const std::map<std::string, Importance>& rows)
{
  const std::vector<std::string> products = {
    "centAsFeat_pt15to30",
    "centAsFeat3x3_pt15to30",
    "centAsFeatBase3x3_pt15to30",
  };
  const std::vector<std::string> productLabels = {
    "Base widths",
    "3#times3 widths",
    "Base + 3#times3",
  };
  const std::vector<std::string> features = {
    "cluster_weta_cogx",
    "cluster_wphi_cogx",
    "cluster_weta33_cogx",
    "cluster_wphi33_cogx",
  };
  const std::vector<std::string> featureLabels = {
    "Full #eta width",
    "Full #phi width",
    "Local 3#times3 #eta width",
    "Local 3#times3 #phi width",
  };
  const std::vector<int> colors = {
    kBlue + 1,
    kAzure + 7,
    kOrange + 7,
    kRed + 1,
  };

  TCanvas c("c_width_gain_stack", "c_width_gain_stack", 1250, 760);
  c.SetLeftMargin(0.12);
  c.SetRightMargin(0.30);
  c.SetTopMargin(0.20);
  c.SetBottomMargin(0.18);

  TH2D frame("frame_stack", ";BDT trained with;Fraction of split gain from width inputs [%]",
             3, 0, 3, 10, 0, 32.0);
  frame.SetStats(false);
  frame.GetXaxis()->SetTitleSize(0.043);
  frame.GetXaxis()->SetTitleOffset(1.22);
  frame.GetXaxis()->SetLabelSize(0.037);
  frame.GetYaxis()->SetTitleSize(0.043);
  frame.GetYaxis()->SetLabelSize(0.037);
  for (int i = 0; i < 3; ++i)
  {
    frame.GetXaxis()->SetBinLabel(i + 1, productLabels[i].c_str());
  }
  frame.Draw("AXIS");

  auto* leg = new TLegend(0.72, 0.56, 0.98, 0.78);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.030);

  const double width = 0.56;
  for (int ip = 0; ip < 3; ++ip)
  {
    const double x0 = ip + 0.5 - 0.5 * width;
    const double x1 = ip + 0.5 + 0.5 * width;
    double y0 = 0.0;
    for (int ifeature = 0; ifeature < 4; ++ifeature)
    {
      const double v = GainPercent(rows, products[ip], features[ifeature]);
      if (std::isnan(v) || v < 0.05) continue;
      auto* box = new TBox(x0, y0, x1, y0 + v);
      box->SetFillColor(colors[ifeature]);
      box->SetLineColor(kWhite);
      box->SetLineWidth(2);
      box->Draw();
      if (ip == 0 || (ip == 1 && ifeature >= 2))
      {
        leg->AddEntry(box, featureLabels[ifeature].c_str(), "f");
      }
      y0 += v;
    }
  }
  leg->Draw();

  TLatex ndc;
  ndc.SetNDC();
  ndc.SetTextFont(42);
  ndc.SetTextAlign(22);
  ndc.SetTextSize(0.040);
  ndc.DrawLatex(0.54, 0.955, "XGBoost split-gain fraction from shower-width inputs");
  DrawSphenixLabel(0.15, 0.885);

  c.SaveAs((outDir + "/widthstudy_width_feature_gain_stacked.png").c_str());
}

void DrawWidthSeparationBars(const std::string& outDir)
{
  const std::string featureSummary =
    "dataOutput/auauTightBDTValidation/model_validation_condor_20260511_110255/validation_feature_summary.csv";
  std::ifstream in(featureSummary);
  std::map<std::string, std::map<std::string, std::pair<double, double>>> moments;
  std::string line;
  std::getline(in, line);
  while (std::getline(in, line))
  {
    const auto cols = SplitCsvLine(line);
    if (cols.size() < 6) continue;
    if (cols[1] != "inclusive") continue;
    moments[cols[0]][cols[2]] = {std::stod(cols[4]), std::stod(cols[5])};
  }

  const std::vector<std::string> features = {
    "cluster_weta_cogx",
    "cluster_wphi_cogx",
    "cluster_weta33_cogx",
    "cluster_wphi33_cogx",
  };
  const std::vector<std::string> labels = {
    "Full #eta width",
    "Full #phi width",
    "Local 3#times3 #eta width",
    "Local 3#times3 #phi width",
  };
  const std::vector<int> colors = {
    kBlue + 1,
    kAzure + 7,
    kOrange + 7,
    kRed + 1,
  };

  std::vector<double> sep;
  for (const auto& feature : features)
  {
    const auto sig = moments[feature]["signal"];
    const auto bkg = moments[feature]["background"];
    const double pooled = std::sqrt(0.5 * (sig.second * sig.second + bkg.second * bkg.second));
    sep.push_back(pooled > 0 ? std::fabs(bkg.first - sig.first) / pooled : 0.0);
  }

  TCanvas c("c_width_separation", "c_width_separation", 1180, 760);
  c.SetLeftMargin(0.25);
  c.SetRightMargin(0.05);
  c.SetTopMargin(0.20);
  c.SetBottomMargin(0.14);

  TH2D frame("frame_sep", ";Single-variable signal/background separation;",
             10, 0, 1.05, 4, 0, 4);
  frame.SetStats(false);
  frame.GetXaxis()->SetTitleSize(0.043);
  frame.GetXaxis()->SetLabelSize(0.038);
  frame.GetYaxis()->SetLabelSize(0.0);
  frame.GetYaxis()->SetTickLength(0.0);
  frame.GetYaxis()->SetNdivisions(0);
  frame.Draw("AXIS");

  TLatex text;
  text.SetTextFont(42);
  text.SetTextAlign(32);
  text.SetTextSize(0.034);
  for (int i = 0; i < 4; ++i)
  {
    const int row = 3 - i;
    const double y = row + 0.5;
    auto* box = new TBox(0.0, y - 0.28, sep[i], y + 0.28);
    box->SetFillColor(colors[i]);
    box->SetLineColor(colors[i]);
    box->Draw();
    text.DrawLatex(-0.025, y, labels[i].c_str());
  }

  TLine grid;
  grid.SetLineColorAlpha(kGray + 1, 0.25);
  for (double x = 0.2; x < 1.05; x += 0.2) grid.DrawLine(x, 0.0, x, 4.0);

  TLatex ndc;
  ndc.SetNDC();
  ndc.SetTextFont(42);
  ndc.SetTextAlign(22);
  ndc.SetTextSize(0.040);
  ndc.DrawLatex(0.56, 0.955, "Single-variable separation for combined-width BDT inputs");
  DrawSphenixLabel(0.29, 0.885);

  c.SaveAs((outDir + "/widthstudy_combined_width_single_variable_separation.png").c_str());
}

void DrawSlideTablePng(const std::string& outDir,
                       const std::map<std::string, Importance>& rows)
{
  const std::vector<std::string> products = {
    "centAsFeat_pt15to30",
    "centAsFeat3x3_pt15to30",
    "centAsFeatBase3x3_pt15to30",
  };
  const std::vector<std::string> labels = {
    "Base widths",
    "3#times3 widths",
    "Base + 3#times3",
  };
  const std::vector<double> auc = {0.7694, 0.7683, 0.7738};

  TCanvas c("c_width_table", "c_width_table", 980, 420);
  c.SetMargin(0, 0, 0, 0);

  auto drawCell = [](double x0, double y0, double x1, double y1, int color) {
    auto* box = new TBox(x0, y0, x1, y1);
    box->SetFillColor(color);
    box->SetLineColor(kGray + 1);
    box->SetLineWidth(1);
    box->Draw();
  };

  const double x[] = {0.02, 0.35, 0.55, 0.75, 0.94};
  const double yTop = 0.78;
  const double rowH = 0.20;
  const int headerColor = TColor::GetColor("#e8e8e8");
  const int rowColor = TColor::GetColor("#f8f8f8");

  drawCell(x[0], yTop, x[1], yTop + rowH, headerColor);
  drawCell(x[1], yTop, x[2], yTop + rowH, headerColor);
  drawCell(x[2], yTop, x[3], yTop + rowH, headerColor);
  drawCell(x[3], yTop, x[4], yTop + rowH, headerColor);

  TLatex t;
  t.SetNDC();
  t.SetTextFont(42);
  t.SetTextAlign(22);
  t.SetTextFont(62);
  t.SetTextSize(0.034);
  t.DrawLatex(0.185, yTop + 0.095, "BDT input set");
  t.DrawLatex(0.450, yTop + 0.095, "Full widths");
  t.DrawLatex(0.650, yTop + 0.095, "3#times3 widths");
  t.DrawLatex(0.845, yTop + 0.095, "ROC-AUC");
  t.SetTextFont(42);

  for (int i = 0; i < 3; ++i)
  {
    const double y1 = yTop - (i + 1) * rowH;
    const double y0 = y1 + rowH;
    for (int col = 0; col < 4; ++col) drawCell(x[col], y1, x[col + 1], y0, rowColor);
    const double fullGain = GainPercent(rows, products[i], "cluster_weta_cogx") + GainPercent(rows, products[i], "cluster_wphi_cogx");
    const double localGain = GainPercent(rows, products[i], "cluster_weta33_cogx") + GainPercent(rows, products[i], "cluster_wphi33_cogx");
    t.SetTextAlign(22);
    t.SetTextSize(0.034);
    t.DrawLatex(0.185, y1 + 0.095, labels[i].c_str());
    t.DrawLatex(0.450, y1 + 0.095, fullGain > 0.05 ? Form("%.1f%%", fullGain) : "not used");
    t.DrawLatex(0.650, y1 + 0.095, localGain > 0.05 ? Form("%.1f%%", localGain) : "not used");
    t.DrawLatex(0.845, y1 + 0.095, Form("%.4f", auc[i]));
  }

  c.SaveAs((outDir + "/widthstudy_slide16_table.png").c_str());
}
}

void MakeWidthStudyFeatureImportance()
{
  SetsPhenixStyle();
  gStyle->SetOptStat(0);

  const std::string outDir = "dataOutput/jstg_slide_candidates/slide16_widthstudy_feature_importance";
  const std::string csvPath = outDir + "/widthstudy_xgb_feature_importance.csv";
  gSystem->mkdir(outDir.c_str(), true);

  const auto rows = ReadImportance(csvPath);
  DrawFeatureGainHeatmap(outDir, rows);
  DrawWidthOnlyBars(outDir, rows);
  DrawWidthContributionStack(outDir, rows);
  DrawWidthSeparationBars(outDir);
  DrawSlideTablePng(outDir, rows);
  std::cout << "[MakeWidthStudyFeatureImportance] wrote " << outDir << "\n";
}
