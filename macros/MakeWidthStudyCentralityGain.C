#include "sPhenixStyle.C"

#include <TArrow.h>
#include <TBox.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
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
struct Row
{
  std::string cent;
  std::string centLabel;
  std::string product;
  std::string modelLabel;
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

std::map<std::string, Row> ReadRows(const std::string& path)
{
  std::ifstream in(path);
  std::map<std::string, Row> rows;
  std::string line;
  std::getline(in, line);
  while (std::getline(in, line))
  {
    const auto cols = SplitCsvLine(line);
    if (cols.size() < 7) continue;
    Row r;
    r.cent = cols[0];
    r.centLabel = cols[1];
    r.product = cols[2];
    r.modelLabel = cols[3];
    r.auc = std::stod(cols[4]);
    r.gainPct = std::stod(cols[5]);
    rows[r.cent + "|" + r.product] = r;
  }
  return rows;
}

const Row& Get(const std::map<std::string, Row>& rows,
               const std::string& cent,
               const std::string& product)
{
  return rows.at(cent + "|" + product);
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

void DrawGainBars(const std::string& outDir,
                  const std::map<std::string, Row>& rows)
{
  const std::vector<std::string> cents = {"0_20", "20_50", "50_80"};
  const std::vector<std::string> labels = {"0-20%", "20-50%", "50-80%"};
  TCanvas c("c_gain_bars", "c_gain_bars", 1180, 760);
  c.SetLeftMargin(0.13);
  c.SetRightMargin(0.05);
  c.SetTopMargin(0.18);
  c.SetBottomMargin(0.15);

  TH2D frame("frame", ";Centrality;ROC-AUC gain relative to base-width BDT [%]",
             3, 0, 3, 10, -0.95, 1.35);
  frame.SetStats(false);
  frame.GetXaxis()->SetTitleSize(0.044);
  frame.GetXaxis()->SetLabelSize(0.038);
  frame.GetYaxis()->SetTitleSize(0.044);
  frame.GetYaxis()->SetLabelSize(0.038);
  for (int i = 0; i < 3; ++i) frame.GetXaxis()->SetBinLabel(i + 1, labels[i].c_str());
  frame.Draw("AXIS");

  TLine zero(0, 0, 3, 0);
  zero.SetLineColor(kGray + 2);
  zero.SetLineWidth(2);
  zero.Draw();

  const double barW = 0.26;
  const int colorOnly = kOrange + 7;
  const int colorBoth = kGreen + 2;
  for (int i = 0; i < 3; ++i)
  {
    const double x = i + 0.5;
    const double gOnly = Get(rows, cents[i], "centAsFeat3x3_pt15to30").gainPct;
    const double gBoth = Get(rows, cents[i], "centAsFeatBase3x3_pt15to30").gainPct;
    for (const auto& item : {std::pair<double, int>{gOnly, colorOnly}, std::pair<double, int>{gBoth, colorBoth}})
    {
      const bool both = item.second == colorBoth;
      const double xc = x + (both ? 0.16 : -0.16);
      const double y0 = std::min(0.0, item.first);
      const double y1 = std::max(0.0, item.first);
      auto* box = new TBox(xc - barW / 2.0, y0, xc + barW / 2.0, y1);
      box->SetFillColor(item.second);
      box->SetLineColor(item.second);
      box->Draw();
    }
  }

  auto* leg = new TLegend(0.60, 0.72, 0.93, 0.83);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.032);
  auto* g1 = new TGraph();
  g1->SetMarkerColor(colorOnly);
  g1->SetFillColor(colorOnly);
  auto* g2 = new TGraph();
  g2->SetMarkerColor(colorBoth);
  g2->SetFillColor(colorBoth);
  leg->AddEntry(g1, "3#times3 widths only", "f");
  leg->AddEntry(g2, "Base + 3#times3 widths", "f");
  leg->Draw();

  TLatex text;
  text.SetNDC();
  text.SetTextFont(42);
  text.SetTextAlign(22);
  text.SetTextSize(0.040);
  text.DrawLatex(0.54, 0.955, "Centrality-binned ROC-AUC gain from local 3#times3 widths");
  DrawSphenixLabel(0.17, 0.885);

  c.SaveAs((outDir + "/widthstudy_auc_gain_vs_centrality_bars.png").c_str());
}

void DrawArrowPlot(const std::string& outDir,
                   const std::map<std::string, Row>& rows)
{
  const std::vector<std::string> cents = {"0_20", "20_50", "50_80"};
  const std::vector<std::string> labels = {"0-20%", "20-50%", "50-80%"};
  TCanvas c("c_arrow", "c_arrow", 1180, 760);
  c.SetLeftMargin(0.12);
  c.SetRightMargin(0.05);
  c.SetTopMargin(0.18);
  c.SetBottomMargin(0.15);

  TH2D frame("frame_arrow", ";Centrality;ROC-AUC",
             3, 0, 3, 10, 0.745, 0.785);
  frame.SetStats(false);
  frame.GetXaxis()->SetTitleSize(0.044);
  frame.GetXaxis()->SetLabelSize(0.038);
  frame.GetYaxis()->SetTitleSize(0.044);
  frame.GetYaxis()->SetLabelSize(0.038);
  for (int i = 0; i < 3; ++i) frame.GetXaxis()->SetBinLabel(i + 1, labels[i].c_str());
  frame.Draw("AXIS");

  const int baseColor = kGray + 2;
  const int bothColor = kGreen + 2;
  auto* baseGraph = new TGraph();
  auto* bothGraph = new TGraph();
  baseGraph->SetMarkerStyle(20);
  baseGraph->SetMarkerSize(1.8);
  baseGraph->SetMarkerColor(baseColor);
  baseGraph->SetLineColor(baseColor);
  bothGraph->SetMarkerStyle(20);
  bothGraph->SetMarkerSize(1.8);
  bothGraph->SetMarkerColor(bothColor);
  bothGraph->SetLineColor(bothColor);

  TLatex lab;
  lab.SetTextFont(42);
  lab.SetTextAlign(22);
  lab.SetTextSize(0.030);
  for (int i = 0; i < 3; ++i)
  {
    const double xBase = i + 0.38;
    const double xBoth = i + 0.62;
    const double aucBase = Get(rows, cents[i], "centAsFeat_pt15to30").auc;
    const double aucBoth = Get(rows, cents[i], "centAsFeatBase3x3_pt15to30").auc;
    const double gain = Get(rows, cents[i], "centAsFeatBase3x3_pt15to30").gainPct;
    auto* arrow = new TArrow(xBase + 0.03, aucBase, xBoth - 0.03, aucBoth, 0.018, "|>");
    arrow->SetLineColor(kGreen + 3);
    arrow->SetFillColor(kGreen + 3);
    arrow->SetLineWidth(3);
    arrow->Draw();
    baseGraph->SetPoint(i, xBase, aucBase);
    bothGraph->SetPoint(i, xBoth, aucBoth);
    lab.SetTextColor(kGreen + 3);
    lab.DrawLatex(i + 0.5, std::max(aucBase, aucBoth) + 0.0025, Form("+%.2f%%", gain));
  }
  baseGraph->Draw("P");
  bothGraph->Draw("P");
  lab.SetTextColor(kBlack);

  auto* leg = new TLegend(0.58, 0.72, 0.92, 0.83);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.032);
  leg->AddEntry(baseGraph, "Base widths", "p");
  leg->AddEntry(bothGraph, "Base + 3#times3 widths", "p");
  leg->Draw();

  TLatex text;
  text.SetNDC();
  text.SetTextFont(42);
  text.SetTextAlign(22);
  text.SetTextSize(0.040);
  text.DrawLatex(0.54, 0.955, "Centrality-binned ROC-AUC before and after adding 3#times3 widths");
  DrawSphenixLabel(0.16, 0.885);

  c.SaveAs((outDir + "/widthstudy_auc_before_after_arrows.png").c_str());
}

void DrawTwoPanel(const std::string& outDir,
                  const std::map<std::string, Row>& rows)
{
  const std::vector<std::string> cents = {"0_20", "20_50", "50_80"};
  const std::vector<std::string> labels = {"0-20%", "20-50%", "50-80%"};
  TCanvas c("c_twopanel", "c_twopanel", 1180, 900);
  c.Divide(1, 2, 0.0, 0.02);

  c.cd(1);
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.05);
  gPad->SetTopMargin(0.20);
  gPad->SetBottomMargin(0.06);
  TH2D top("top", "; ;ROC-AUC", 3, 0, 3, 10, 0.745, 0.785);
  top.SetStats(false);
  top.GetXaxis()->SetLabelSize(0.0);
  top.GetXaxis()->SetTickLength(0.0);
  top.GetYaxis()->SetTitleSize(0.050);
  top.GetYaxis()->SetLabelSize(0.043);
  top.Draw("AXIS");

  auto* gBase = new TGraph();
  auto* gOnly = new TGraph();
  auto* gBoth = new TGraph();
  for (int i = 0; i < 3; ++i)
  {
    gBase->SetPoint(i, i + 0.5, Get(rows, cents[i], "centAsFeat_pt15to30").auc);
    gOnly->SetPoint(i, i + 0.5, Get(rows, cents[i], "centAsFeat3x3_pt15to30").auc);
    gBoth->SetPoint(i, i + 0.5, Get(rows, cents[i], "centAsFeatBase3x3_pt15to30").auc);
  }
  gBase->SetMarkerStyle(20); gBase->SetMarkerSize(1.7); gBase->SetMarkerColor(kGray + 2); gBase->SetLineColor(kGray + 2); gBase->SetLineWidth(3);
  gOnly->SetMarkerStyle(21); gOnly->SetMarkerSize(1.7); gOnly->SetMarkerColor(kOrange + 7); gOnly->SetLineColor(kOrange + 7); gOnly->SetLineWidth(3);
  gBoth->SetMarkerStyle(20); gBoth->SetMarkerSize(1.9); gBoth->SetMarkerColor(kGreen + 2); gBoth->SetLineColor(kGreen + 2); gBoth->SetLineWidth(4);
  gBase->Draw("PL"); gOnly->Draw("PL"); gBoth->Draw("PL");
  auto* leg = new TLegend(0.56, 0.66, 0.92, 0.82);
  leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextFont(42); leg->SetTextSize(0.035);
  leg->AddEntry(gBase, "Base widths", "pl");
  leg->AddEntry(gOnly, "3#times3 widths only", "pl");
  leg->AddEntry(gBoth, "Base + 3#times3 widths", "pl");
  leg->Draw();
  TLatex text; text.SetNDC(); text.SetTextFont(42); text.SetTextAlign(22); text.SetTextSize(0.045);
  text.DrawLatex(0.54, 0.94, "ROC-AUC vs centrality for width-study BDTs");
  DrawSphenixLabel(0.16, 0.86);

  c.cd(2);
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.05);
  gPad->SetTopMargin(0.05);
  gPad->SetBottomMargin(0.24);
  TH2D bot("bot", ";Centrality;ROC-AUC gain relative to base-width BDT [%]", 3, 0, 3, 10, -0.95, 1.35);
  bot.SetStats(false);
  bot.GetXaxis()->SetTitleSize(0.050); bot.GetXaxis()->SetLabelSize(0.043);
  bot.GetYaxis()->SetTitleSize(0.050); bot.GetYaxis()->SetLabelSize(0.043);
  for (int i = 0; i < 3; ++i) bot.GetXaxis()->SetBinLabel(i + 1, labels[i].c_str());
  bot.Draw("AXIS");
  TLine zero(0, 0, 3, 0); zero.SetLineColor(kGray + 2); zero.SetLineWidth(2); zero.Draw();
  auto* dgOnly = new TGraph();
  auto* dgBoth = new TGraph();
  for (int i = 0; i < 3; ++i)
  {
    dgOnly->SetPoint(i, i + 0.5, Get(rows, cents[i], "centAsFeat3x3_pt15to30").gainPct);
    dgBoth->SetPoint(i, i + 0.5, Get(rows, cents[i], "centAsFeatBase3x3_pt15to30").gainPct);
  }
  dgOnly->SetMarkerStyle(21); dgOnly->SetMarkerSize(1.7); dgOnly->SetMarkerColor(kOrange + 7); dgOnly->SetLineColor(kOrange + 7); dgOnly->SetLineWidth(3);
  dgBoth->SetMarkerStyle(20); dgBoth->SetMarkerSize(1.9); dgBoth->SetMarkerColor(kGreen + 2); dgBoth->SetLineColor(kGreen + 2); dgBoth->SetLineWidth(4);
  dgOnly->Draw("PL"); dgBoth->Draw("PL");

  c.SaveAs((outDir + "/widthstudy_auc_centrality_twopanel.png").c_str());
}
}

void MakeWidthStudyCentralityGain()
{
  SetsPhenixStyle();
  gStyle->SetOptStat(0);
  const std::string outDir = "dataOutput/jstg_slide_candidates/slide17_widthstudy_centrality_gain";
  gSystem->mkdir(outDir.c_str(), true);
  const auto rows = ReadRows(outDir + "/widthstudy_centrality_auc_summary.csv");
  DrawGainBars(outDir, rows);
  DrawArrowPlot(outDir, rows);
  DrawTwoPanel(outDir, rows);
  std::cout << "[MakeWidthStudyCentralityGain] wrote " << outDir << "\n";
}
