#include "sPhenixStyle.C"

#include <TBox.h>
#include <TCanvas.h>
#include <TColor.h>
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
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace
{
struct Ranking
{
  std::string label;
  double aucInclusive = NAN;
  double aucPtMean = NAN;
  double aucCentMean = NAN;
  double aucPtCentMean = NAN;
};

struct AucCell
{
  double auc = NAN;
  int entries = 0;
};

std::vector<std::string> SplitCsvLine(const std::string& line)
{
  std::vector<std::string> out;
  std::string cur;
  bool quoted = false;
  for (char c : line)
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

double ToDouble(const std::string& s)
{
  if (s == "nan" || s == "NaN" || s.empty()) return NAN;
  return std::stod(s);
}

std::map<std::string, Ranking> ReadRankings(const std::string& path)
{
  std::ifstream in(path);
  std::map<std::string, Ranking> rows;
  std::string line;
  std::getline(in, line);
  while (std::getline(in, line))
  {
    const auto cols = SplitCsvLine(line);
    if (cols.size() < 8) continue;
    Ranking r;
    r.label = cols[1];
    r.aucInclusive = ToDouble(cols[2]);
    r.aucPtMean = ToDouble(cols[3]);
    r.aucCentMean = ToDouble(cols[4]);
    r.aucPtCentMean = ToDouble(cols[5]);
    rows[cols[0]] = r;
  }
  return rows;
}

std::map<std::string, AucCell> ReadAucCells(const std::string& path)
{
  std::ifstream in(path);
  std::map<std::string, AucCell> rows;
  std::string line;
  std::getline(in, line);
  while (std::getline(in, line))
  {
    const auto cols = SplitCsvLine(line);
    if (cols.size() < 5) continue;
    AucCell c;
    c.auc = ToDouble(cols[3]);
    c.entries = std::stoi(cols[4]);
    rows[cols[0] + "|" + cols[1] + "|" + cols[2]] = c;
  }
  return rows;
}

void DrawSphenixLabel(double x, double y)
{
  TLatex text;
  text.SetNDC();
  text.SetTextAlign(11);
  text.SetTextFont(42);
  text.SetTextSize(0.035);
  text.DrawLatex(x, y, "#it{#bf{sPHENIX}} Internal");
  text.SetTextSize(0.028);
  text.DrawLatex(x, y - 0.045,
                 "Pythia overlay, #sqrt{s_{NN}} = 200 GeV; 15 #leq E_{T}^{cluster} < 30 GeV");
}

void DrawSummaryPlot(const std::string& outDir,
                     const std::map<std::string, Ranking>& rank)
{
  const std::vector<std::string> products = {
    "centAsFeat_pt15to30",
    "centAsFeat3x3_pt15to30",
    "centAsFeatBase3x3_pt15to30",
  };
  const std::vector<std::string> labels = {
    "Base widths",
    "3#times3 widths only",
    "Base + 3#times3 widths",
  };
  const int colors[] = {kBlue + 1, kOrange + 7, kGreen + 2};
  TCanvas c("c", "c", 1320, 760);
  c.SetLeftMargin(0.13);
  c.SetRightMargin(0.05);
  c.SetTopMargin(0.20);
  c.SetBottomMargin(0.20);

  TH2D frame("frame", ";Width-feature option;Inclusive ROC AUC", 3, 0, 3, 10, 0.760, 0.778);
  frame.SetStats(false);
  frame.GetXaxis()->SetTitleSize(0.046);
  frame.GetXaxis()->SetLabelSize(0.038);
  frame.GetYaxis()->SetTitleSize(0.046);
  frame.GetYaxis()->SetLabelSize(0.038);
  frame.GetXaxis()->SetBinLabel(1, "Base widths");
  frame.GetXaxis()->SetBinLabel(2, "3#times3 widths only");
  frame.GetXaxis()->SetBinLabel(3, "Base + 3#times3 widths");
  frame.Draw("AXIS");

  TLatex text;
  text.SetNDC();
  text.SetTextFont(42);
  const double base = rank.at(products[0]).aucInclusive;
  int bestIdx = 0;
  double best = base;
  for (size_t j = 0; j < products.size(); ++j)
  {
    const double val = rank.at(products[j]).aucInclusive;
    if (val > best) { best = val; bestIdx = j; }
    auto* stem = new TLine(j + 0.5, 0.760, j + 0.5, val);
    stem->SetLineColor(colors[j]);
    stem->SetLineWidth(5);
    stem->Draw();
    auto* marker = new TGraph(1);
    marker->SetPoint(0, j + 0.5, val);
    marker->SetMarkerStyle(20);
    marker->SetMarkerSize(3.0);
    marker->SetMarkerColor(colors[j]);
    marker->SetLineColor(colors[j]);
    marker->Draw("P");

    const double pct = (val - base) / base * 100.0;
    std::ostringstream line1;
    std::ostringstream line2;
    line1 << "AUC " << std::fixed << std::setprecision(4) << val;
    if (j > 0)
    {
      line2 << std::showpos << std::setprecision(2) << pct << "% vs base";
    }
    else
    {
      line2 << "reference";
    }
    text.SetTextAlign(22);
    text.SetTextSize(0.029);
    text.SetTextColor(kBlack);
    const double xText = 0.25 + j * 0.275;
    const double yText = j == 2 ? 0.56 : (j == 1 ? 0.45 : 0.50);
    text.DrawLatex(xText, yText, line1.str().c_str());
    text.SetTextSize(0.026);
    text.DrawLatex(xText, yText - 0.043, line2.str().c_str());
  }

  text.SetTextAlign(22);
  text.SetTextSize(0.044);
  text.DrawLatex(0.50, 0.945, "15-30 GeV centrality-input BDT: width-feature comparison");
  DrawSphenixLabel(0.16, 0.875);
  text.SetTextAlign(11);
  text.SetTextSize(0.031);
  text.SetTextColor(colors[bestIdx]);
  text.DrawLatex(0.54, 0.875, (std::string("Winner: ") + labels[bestIdx]).c_str());
  text.SetTextColor(kGray + 2);
  text.SetTextSize(0.028);
  text.DrawLatex(0.54, 0.830,
                 "Best model combines full-cluster widths with local 3#times3 widths.");
  text.SetTextColor(kBlack);

  c.SaveAs((outDir + "/widthstudy_auc_summary_oneplot.png").c_str());
}

void DrawDeltaHeatmap(const std::string& outDir,
                      const std::map<std::string, AucCell>& cells,
                      const std::vector<std::string>& pt = {"15_20", "20_25", "25_35"},
                      const std::vector<std::string>& ptLabel = {"15-20", "20-25", "25-35"},
                      const std::string& outName = "widthstudy_auc_delta_heatmap_oneplot.png")
{
  const std::vector<std::string> cent = {"0_20", "20_50", "50_80"};
  const std::vector<std::string> centLabel = {"0-20%", "20-50%", "50-80%"};
  const std::vector<std::string> products = {
    "centAsFeat_pt15to30",
    "centAsFeat3x3_pt15to30",
    "centAsFeatBase3x3_pt15to30",
  };
  const std::vector<std::string> labels = {
    "3#times3 widths only",
    "Base + 3#times3 widths",
  };

  const Int_t nStops = 5;
  Double_t stops[nStops] = {0.00, 0.2857, 0.48, 0.72, 1.00};
  Double_t red[nStops] = {0.129, 0.965, 0.996, 0.878, 0.698};
  Double_t green[nStops] = {0.400, 0.965, 0.878, 0.486, 0.094};
  Double_t blue[nStops] = {0.675, 0.965, 0.565, 0.286, 0.168};
  TColor::CreateGradientColorTable(nStops, stops, red, green, blue, 255);
  gStyle->SetNumberContours(255);

  TCanvas c("c2", "c2", 1550, 760);
  c.SetLeftMargin(0.18);
  c.SetRightMargin(0.14);
  c.SetTopMargin(0.24);
  c.SetBottomMargin(0.17);

  std::vector<std::vector<std::string>> cellText(2, std::vector<std::string>(9));
  std::vector<std::vector<double>> deltaText(2, std::vector<double>(9, 0.0));
  std::vector<int> winnerRow(9, 0);
  for (int ic = 0; ic < 3; ++ic)
  {
    for (int ip = 0; ip < 3; ++ip)
    {
      const int xbin = ic * 3 + ip + 1;
      const std::string keyBase = products[0] + "|" + cent[ic] + "|" + pt[ip];
      const double base = cells.at(keyBase).auc;
      double bestVal = base;
      int bestRow = -1;
      for (int row = 0; row < 2; ++row)
      {
        const std::string p = products[row + 1];
        const double val = cells.at(p + "|" + cent[ic] + "|" + pt[ip]).auc;
        const double delta = (val - base) / base * 100.0;
        deltaText[row][xbin - 1] = delta;
        std::ostringstream ss;
        ss << std::showpos << std::fixed << std::setprecision(1) << delta << "%";
        cellText[row][xbin - 1] = ss.str();
        if (val > bestVal)
        {
          bestVal = val;
          bestRow = row;
        }
      }
      winnerRow[xbin - 1] = bestRow;
    }
  }

  auto* h = new TH2D("h_delta", ";Photon candidate E_{T} bin within each centrality group [GeV];",
                     9, 0, 9, 2, 0, 2);
  h->SetDirectory(nullptr);
  h->SetStats(false);
  h->SetMinimum(-1.0);
  h->SetMaximum(2.5);
  h->GetXaxis()->SetTitleSize(0.045);
  h->GetXaxis()->SetLabelSize(0.035);
  h->GetXaxis()->SetTitleOffset(1.10);
  h->GetYaxis()->SetTitleSize(0.045);
  h->GetYaxis()->SetLabelSize(0.038);
  h->GetYaxis()->SetTitleOffset(1.95);
  h->GetZaxis()->SetTitle("Percent gain wrt base ROC-AUC [%]");
  h->GetZaxis()->SetTitleSize(0.037);
  h->GetZaxis()->SetTitleOffset(1.12);
  h->GetZaxis()->SetLabelSize(0.033);
  h->GetYaxis()->SetBinLabel(2, labels[0].c_str());
  h->GetYaxis()->SetBinLabel(1, labels[1].c_str());
  for (int ic = 0; ic < 3; ++ic)
  {
    for (int ip = 0; ip < 3; ++ip)
    {
      const int globalCol = ic * 3 + ip;
      h->GetXaxis()->SetBinLabel(globalCol + 1, ptLabel[ip].c_str());
      h->SetBinContent(globalCol + 1, 2, deltaText[0][globalCol]);
      h->SetBinContent(globalCol + 1, 1, deltaText[1][globalCol]);
    }
  }
  h->SetContour(255);
  h->Draw("COLZ");

  TLatex text;
  text.SetTextFont(42);
  text.SetTextAlign(22);
  for (int row = 0; row < 2; ++row)
  {
    for (int col = 0; col < 9; ++col)
    {
      const double v = deltaText[row][col];
      text.SetTextSize(0.039);
      text.SetTextColor(std::fabs(v) > 1.8 ? kWhite : kBlack);
      text.DrawLatex(col + 0.5, 1.5 - row, cellText[row][col].c_str());
    }
  }
  text.SetTextColor(kBlack);

  TLine rowLine(0, 1, 9, 1);
  rowLine.SetLineWidth(2);
  rowLine.SetLineColor(kGray + 2);
  rowLine.Draw();

  TLine outerTop(0, 2, 9, 2);
  TLine outerBottom(0, 0, 9, 0);
  TLine outerLeft(0, 0, 0, 2);
  TLine outerRight(9, 0, 9, 2);
  for (auto* l : {&outerTop, &outerBottom, &outerLeft, &outerRight})
  {
    l->SetLineWidth(2);
    l->SetLineColor(kGray + 2);
    l->Draw();
  }

  for (double x : {3.0, 6.0})
  {
    auto* gutter = new TBox(x - 0.045, 0.0, x + 0.045, 2.0);
    gutter->SetFillColor(kWhite);
    gutter->SetLineColor(kWhite);
    gutter->Draw();
    TLine divider(x, 0.0, x, 2.0);
    divider.SetLineWidth(2);
    divider.SetLineColor(kGray + 1);
    divider.SetLineStyle(2);
    divider.Draw();
  }

  TLatex dataText;
  dataText.SetTextFont(42);
  dataText.SetTextAlign(22);
  dataText.SetTextSize(0.036);
  for (int ic = 0; ic < 3; ++ic)
  {
    const double x0 = ic * 3.0;
    const double x1 = x0 + 3.0;
    const double xc = 0.5 * (x0 + x1);
    TLine bracketTop(x0 + 0.08, 2.06, x1 - 0.08, 2.06);
    TLine bracketL(x0 + 0.08, 2.02, x0 + 0.08, 2.06);
    TLine bracketR(x1 - 0.08, 2.02, x1 - 0.08, 2.06);
    for (auto* l : {&bracketTop, &bracketL, &bracketR})
    {
      l->SetLineWidth(2);
      l->SetLineColor(kGray + 2);
      l->Draw();
    }
    dataText.DrawLatex(xc, 2.16, centLabel[ic].c_str());
  }

  TLatex ndc;
  ndc.SetNDC();
  ndc.SetTextFont(42);
  ndc.SetTextAlign(22);
  ndc.SetTextSize(0.043);
  ndc.DrawLatex(0.50, 0.965, "Percent gain in ROC-AUC from local 3#times3 shower-width inputs");
  DrawSphenixLabel(0.20, 0.905);

  c.SaveAs((outDir + "/" + outName).c_str());
}
}

void MakeWidthStudyBDTComparison()
{
  SetsPhenixStyle();
  gStyle->SetOptStat(0);
  gSystem->mkdir("dataOutput/jstg_slide_candidates/slide15_widthstudy_bdt_comparison", true);

  const std::string baseDir = "dataOutput/auauTightBDTValidation/model_validation_condor_20260511_110255";
  const std::string outDir = "dataOutput/jstg_slide_candidates/slide15_widthstudy_bdt_comparison";
  const auto rank = ReadRankings(baseDir + "/validation_model_rankings.csv");
  const auto cells = ReadAucCells(baseDir + "/validation_auc_table.csv");
  const auto cells2530 = ReadAucCells(baseDir + "/validation_auc_table_pt15_20_20_25_25_30.csv");

  DrawSummaryPlot(outDir, rank);
  DrawDeltaHeatmap(outDir, cells);
  DrawDeltaHeatmap(outDir, cells2530,
                   {"15_20", "20_25", "25_30"},
                   {"15-20", "20-25", "25-30"},
                   "widthstudy_auc_delta_heatmap_oneplot_pt25to30_rootstyle.png");
  std::cout << "[MakeWidthStudyBDTComparison] wrote " << outDir << "\n";
}
