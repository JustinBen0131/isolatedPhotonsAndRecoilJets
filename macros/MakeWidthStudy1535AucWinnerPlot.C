#include "sPhenixStyle.C"

#include <TBox.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TGaxis.h>
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
#include <tuple>
#include <vector>

namespace
{
struct AucRow
{
  double auc = NAN;
  int entries = 0;
  int signal = 0;
  int background = 0;
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

using AucKey = std::tuple<std::string, std::string, std::string>;

std::map<AucKey, AucRow> ReadAucRows(const std::string& path)
{
  std::ifstream in(path);
  std::map<AucKey, AucRow> rows;
  if (!in)
  {
    std::cerr << "Could not open " << path << "\n";
    return rows;
  }

  std::string line;
  std::getline(in, line);
  while (std::getline(in, line))
  {
    const auto cols = SplitCsvLine(line);
    if (cols.size() < 7) continue;
    AucRow r;
    r.auc = cols[3] == "nan" ? NAN : std::stod(cols[3]);
    r.entries = std::stoi(cols[4]);
    r.signal = std::stoi(cols[5]);
    r.background = std::stoi(cols[6]);
    rows[{cols[0], cols[1], cols[2]}] = r;
  }
  return rows;
}

const AucRow& GetAuc(const std::map<AucKey, AucRow>& rows,
                     const std::string& product,
                     const std::string& cent,
                     const std::string& pt)
{
  return rows.at({product, cent, pt});
}

double GainPercent(const double val, const double base)
{
  if (!std::isfinite(val) || !std::isfinite(base) || base == 0) return NAN;
  return 100.0 * (val - base) / base;
}

void DrawSphenixLabel(const double x, const double y)
{
  TLatex text;
  text.SetNDC();
  text.SetTextAlign(11);
  text.SetTextFont(42);
  text.SetTextSize(0.030);
  text.DrawLatex(x, y, "#it{#bf{sPHENIX}} Internal");
  text.SetTextSize(0.022);
  text.DrawLatex(x, y - 0.040, "Pythia overlay, #sqrt{s_{NN}} = 200 GeV");
}

std::string FormatAuc(const double v)
{
  std::ostringstream ss;
  ss << std::fixed << std::setprecision(3) << v;
  return ss.str();
}

std::string FormatGain(const double v)
{
  std::ostringstream ss;
  ss << std::showpos << std::fixed << std::setprecision(1) << v << "%";
  return ss.str();
}

std::string FormatDeltaAuc(const double v)
{
  std::ostringstream ss;
  ss << std::showpos << std::fixed << std::setprecision(3) << v;
  return ss.str();
}
}

void MakeWidthStudy1535AucWinnerPlot()
{
  SetsPhenixStyle();
  gStyle->SetOptStat(0);
  TGaxis::SetMaxDigits(3);

  const std::string report = "dataOutput/auauTightBDTValidation/model_validation_condor_20260511_171943";
  const std::string outDir = "dataOutput/jstg_slide_candidates/slide10_widthstudy_1535";
  gSystem->mkdir(outDir.c_str(), true);

  const auto rows = ReadAucRows(report + "/validation_auc_table.csv");

  const std::vector<std::string> products = {
      "centAsFeat_pt15to35",
      "centAsFeat3x3_pt15to35",
      "centAsFeatBase3x3_pt15to35"};
  const std::vector<std::string> modelLabels = {
      "Base widths",
      "3#times3 widths",
      "Base + 3#times3"};
  const std::vector<std::string> modelCsvLabels = {
      "Base widths",
      "3x3 widths",
      "Base + 3x3"};
  const std::vector<int> modelColors = {
      kGray + 1,
      kOrange + 7,
      kTeal + 3};

  const std::vector<std::string> cents = {"0_20", "20_50", "50_80"};
  const std::vector<std::string> centLabels = {"0-20% central", "20-50% mid-cent.", "50-80% peripheral"};
  const std::vector<std::string> pts = {"15_20", "20_25", "25_35"};
  const std::vector<std::string> ptLabels = {"15-20", "20-25", "25-35"};

  std::vector<std::vector<int>> winner(3, std::vector<int>(3, 0));
  std::vector<std::vector<double>> bestAuc(3, std::vector<double>(3, NAN));
  std::vector<std::vector<double>> bestDelta(3, std::vector<double>(3, NAN));
  std::vector<double> sumAuc(3, 0.0);
  std::vector<int> winCount(3, 0);
  int nCells = 0;

  std::ofstream csv(outDir + "/width_inputs_1535_auc_winner_summary.csv");
  csv << "centrality,pt_bin,base_width_auc,width3x3_auc,base_plus_3x3_auc,winner,winner_auc,winner_delta_auc_vs_base,entries,signal_entries,background_entries\n";

  for (int ic = 0; ic < 3; ++ic)
  {
    for (int ip = 0; ip < 3; ++ip)
    {
      std::vector<double> auc(3, NAN);
      for (int im = 0; im < 3; ++im)
      {
        auc[im] = GetAuc(rows, products[im], cents[ic], pts[ip]).auc;
        sumAuc[im] += auc[im];
      }
      ++nCells;

      int best = 0;
      for (int im = 1; im < 3; ++im)
      {
        if (auc[im] > auc[best]) best = im;
      }

      const auto& ref = GetAuc(rows, products[2], cents[ic], pts[ip]);
      winner[ic][ip] = best;
      bestAuc[ic][ip] = auc[best];
      bestDelta[ic][ip] = auc[best] - auc[0];
      ++winCount[best];

      csv << cents[ic] << "," << pts[ip] << ","
          << std::setprecision(10) << auc[0] << "," << auc[1] << "," << auc[2] << ","
          << modelCsvLabels[best] << "," << auc[best] << "," << bestDelta[ic][ip] << ","
          << ref.entries << "," << ref.signal << "," << ref.background << "\n";
    }
  }

  std::vector<double> meanAuc(3, NAN);
  for (int im = 0; im < 3; ++im)
  {
    meanAuc[im] = nCells > 0 ? sumAuc[im] / nCells : NAN;
  }

  TCanvas c("c_width_inputs_1535_auc_winner", "c_width_inputs_1535_auc_winner", 1600, 900);
  c.SetLeftMargin(0.08);
  c.SetRightMargin(0.05);
  c.SetTopMargin(0.17);
  c.SetBottomMargin(0.11);

  TLatex ndc;
  ndc.SetNDC();
  ndc.SetTextFont(42);
  ndc.SetTextColor(kBlack);
  ndc.SetTextAlign(22);
  ndc.SetTextSize(0.036);
  ndc.DrawLatex(0.53, 0.962, "15-35 GeV width-input BDT comparison");
  ndc.SetTextSize(0.023);
  ndc.DrawLatex(0.53, 0.922, "Same training window and labels; only the shower-width input set changes. Higher ROC AUC is better.");
  DrawSphenixLabel(0.085, 0.922);

  TPad left("left", "left", 0.06, 0.13, 0.67, 0.82);
  left.SetLeftMargin(0.16);
  left.SetRightMargin(0.04);
  left.SetTopMargin(0.08);
  left.SetBottomMargin(0.15);
  left.Draw();
  left.cd();

  TH2D frame("frame", ";Truth photon E_{T} bin [GeV];Centrality bin", 3, 0, 3, 3, 0, 3);
  frame.SetStats(false);
  frame.GetXaxis()->SetTitleSize(0.045);
  frame.GetXaxis()->SetLabelSize(0.040);
  frame.GetXaxis()->SetTitleOffset(1.12);
  frame.GetYaxis()->SetTitleSize(0.045);
  frame.GetYaxis()->SetLabelSize(0.038);
  frame.GetYaxis()->SetTitleOffset(1.55);
  for (int ip = 0; ip < 3; ++ip)
  {
    frame.GetXaxis()->SetBinLabel(ip + 1, ptLabels[ip].c_str());
  }
  for (int ic = 0; ic < 3; ++ic)
  {
    frame.GetYaxis()->SetBinLabel(3 - ic, centLabels[ic].c_str());
  }
  frame.Draw("AXIS");

  TBox box;
  box.SetLineColor(kWhite);
  box.SetLineWidth(2);
  TLatex txt;
  txt.SetTextFont(42);
  txt.SetTextAlign(22);

  for (int ic = 0; ic < 3; ++ic)
  {
    for (int ip = 0; ip < 3; ++ip)
    {
      const int best = winner[ic][ip];
      const double x0 = ip;
      const double x1 = ip + 1;
      const double y0 = 2 - ic;
      const double y1 = 3 - ic;
      box.SetFillColor(modelColors[best]);
      box.DrawBox(x0, y0, x1, y1);

      txt.SetTextColor(kWhite);
      txt.SetTextSize(0.050);
      txt.DrawLatex(0.5 * (x0 + x1), y0 + 0.68, modelLabels[best].c_str());
      txt.SetTextSize(0.038);
      txt.DrawLatex(0.5 * (x0 + x1), y0 + 0.42, ("AUC " + FormatAuc(bestAuc[ic][ip])).c_str());
      txt.SetTextSize(0.034);
      txt.DrawLatex(0.5 * (x0 + x1), y0 + 0.21, ("#DeltaAUC vs base " + FormatDeltaAuc(bestDelta[ic][ip])).c_str());
    }
  }

  TLine line;
  line.SetLineColor(kWhite);
  line.SetLineWidth(3);
  for (int i = 0; i <= 3; ++i)
  {
    line.DrawLine(i, 0, i, 3);
    line.DrawLine(0, i, 3, i);
  }
  frame.Draw("AXIS SAME");

  c.cd();
  TPad right("right", "right", 0.69, 0.17, 0.96, 0.78);
  right.SetLeftMargin(0.34);
  right.SetRightMargin(0.04);
  right.SetTopMargin(0.18);
  right.SetBottomMargin(0.17);
  right.Draw();
  right.cd();

  std::vector<int> summaryOrder = {2, 1, 0};
  TH2D bars("bars", ";Mean ROC AUC over 9 validation cells;", 1, 0.700, 0.755, 3, 0, 3);
  bars.SetStats(false);
  bars.GetXaxis()->SetTitleSize(0.050);
  bars.GetXaxis()->SetLabelSize(0.040);
  bars.GetXaxis()->SetNdivisions(505);
  bars.GetXaxis()->SetTitleOffset(1.10);
  bars.GetYaxis()->SetLabelSize(0.049);
  bars.GetYaxis()->SetTickLength(0);
  for (int i = 0; i < 3; ++i)
  {
    bars.GetYaxis()->SetBinLabel(3 - i, modelLabels[summaryOrder[i]].c_str());
  }
  bars.Draw("AXIS");

  TLatex title;
  title.SetNDC();
  title.SetTextFont(42);
  title.SetTextAlign(22);
  title.SetTextSize(0.070);
  title.DrawLatex(0.58, 0.95, "Mean AUC summary");
  title.SetTextSize(0.036);
  title.DrawLatex(0.58, 0.87, "Dashed line: base-width mean");

  TLine xline;
  xline.SetLineColor(kGray + 1);
  xline.SetLineStyle(2);
  xline.SetLineWidth(2);
  xline.DrawLine(meanAuc[0], 0.12, meanAuc[0], 2.88);

  TLatex rtxt;
  rtxt.SetTextFont(42);
  rtxt.SetTextAlign(12);
  for (int i = 0; i < 3; ++i)
  {
    const int im = summaryOrder[i];
    const double yCenter = 2.5 - i;
    const double y0 = yCenter - 0.22;
    const double y1 = yCenter + 0.22;
    box.SetFillColor(modelColors[im]);
    box.SetLineColor(modelColors[im]);
    box.DrawBox(0.700, y0, meanAuc[im], y1);
    rtxt.SetTextColor(kBlack);
    rtxt.SetTextSize(0.052);
    rtxt.DrawLatex(meanAuc[im] + 0.0010, yCenter + 0.06, Form("%.3f", meanAuc[im]));
    rtxt.SetTextSize(0.042);
    rtxt.DrawLatex(meanAuc[im] + 0.0010, yCenter - 0.12, Form("%d/9 best", winCount[im]));
  }

  c.cd();
  TLegend leg(0.705, 0.790, 0.950, 0.885);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextFont(42);
  leg.SetTextSize(0.025);
  TGraph dummyBase;
  TGraph dummy33;
  TGraph dummyBoth;
  dummyBase.SetMarkerStyle(21);
  dummyBase.SetMarkerColor(modelColors[0]);
  dummy33.SetMarkerStyle(21);
  dummy33.SetMarkerColor(modelColors[1]);
  dummyBoth.SetMarkerStyle(21);
  dummyBoth.SetMarkerColor(modelColors[2]);
  leg.AddEntry(&dummyBase, "Winner: base widths", "p");
  leg.AddEntry(&dummy33, "Winner: 3#times3 widths", "p");
  leg.AddEntry(&dummyBoth, "Winner: base + 3#times3", "p");
  leg.Draw();

  ndc.SetTextAlign(22);
  ndc.SetTextSize(0.026);
  ndc.DrawLatex(0.365, 0.080, "Winner map: base + 3#times3 is best in 8 of 9 validation cells");

  const std::string png = outDir + "/width_inputs_1535_auc_winner_map.png";
  c.SaveAs(png.c_str());
  c.SaveAs((outDir + "/width_inputs_1535_auc_winner_map_bar_panel.png").c_str());
  std::cout << "[MakeWidthStudy1535AucWinnerPlot] wrote " << png << "\n";
}
