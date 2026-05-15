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

struct RankRow
{
  double aucInclusive = NAN;
  double aucPtCentMean = NAN;
  int eligibleEntries = 0;
  double signalMean = NAN;
  double backgroundMean = NAN;
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

std::map<std::string, RankRow> ReadRankRows(const std::string& path)
{
  std::ifstream in(path);
  std::map<std::string, RankRow> rows;
  std::string line;
  std::getline(in, line);
  while (std::getline(in, line))
  {
    const auto cols = SplitCsvLine(line);
    if (cols.size() < 10) continue;
    RankRow r;
    r.aucInclusive = std::stod(cols[2]);
    r.aucPtCentMean = std::stod(cols[5]);
    r.eligibleEntries = std::stoi(cols[7]);
    r.signalMean = std::stod(cols[8]);
    r.backgroundMean = std::stod(cols[9]);
    rows[cols[0]] = r;
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
  text.SetTextSize(0.034);
  text.DrawLatex(x, y, "#it{#bf{sPHENIX}} Internal");
  text.SetTextSize(0.026);
  text.DrawLatex(x, y - 0.043, "Pythia overlay, #sqrt{s_{NN}} = 200 GeV");
}

void WriteSummaryCsv(const std::map<std::string, RankRow>& ranks, const std::string& outPath)
{
  const std::vector<std::string> windows = {"pt5to35", "pt10to35", "pt15to35"};
  const std::vector<std::string> labels = {"Base widths", "3x3 widths", "Base + 3x3 widths"};
  const std::vector<std::string> prefixes = {"centAsFeat", "centAsFeat3x3", "centAsFeatBase3x3"};

  std::ofstream out(outPath);
  out << "pt_window,width_set,product,inclusive_auc,pt_cent_mean_auc,inclusive_gain_vs_base_percent,pt_cent_mean_gain_vs_base_percent,eligible_entries,signal_background_mean_score_gap\n";
  for (const auto& window : windows)
  {
    const auto& base = ranks.at("centAsFeat_" + window);
    for (int i = 0; i < static_cast<int>(prefixes.size()); ++i)
    {
      const std::string product = prefixes[i] + "_" + window;
      const auto& row = ranks.at(product);
      out << window << "," << labels[i] << "," << product << ","
          << std::setprecision(10) << row.aucInclusive << "," << row.aucPtCentMean << ","
          << GainPercent(row.aucInclusive, base.aucInclusive) << ","
          << GainPercent(row.aucPtCentMean, base.aucPtCentMean) << ","
          << row.eligibleEntries << "," << row.signalMean - row.backgroundMean << "\n";
    }
  }
}

void WriteCellGainCsv(const std::map<AucKey, AucRow>& rows, const std::string& outPath)
{
  const std::vector<std::string> windows = {"pt5to35", "pt10to35", "pt15to35"};
  const std::vector<std::string> cents = {"0_20", "20_50", "50_80"};
  const std::vector<std::string> pts = {"6_10", "10_15", "15_20", "20_25", "25_35"};

  std::ofstream out(outPath);
  out << "pt_window,centrality,pt_bin,base_width_auc,width3x3_auc,base_plus_3x3_auc,width3x3_gain_percent,base_plus_3x3_gain_percent,entries,signal_entries,background_entries\n";
  for (const auto& window : windows)
  {
    for (const auto& cent : cents)
    {
      for (const auto& pt : pts)
      {
        const auto& base = GetAuc(rows, "centAsFeat_" + window, cent, pt);
        const auto& w33 = GetAuc(rows, "centAsFeat3x3_" + window, cent, pt);
        const auto& both = GetAuc(rows, "centAsFeatBase3x3_" + window, cent, pt);
        if (!std::isfinite(base.auc) || !std::isfinite(w33.auc) || !std::isfinite(both.auc)) continue;
        out << window << "," << cent << "," << pt << ","
            << std::setprecision(10) << base.auc << "," << w33.auc << "," << both.auc << ","
            << GainPercent(w33.auc, base.auc) << "," << GainPercent(both.auc, base.auc) << ","
            << both.entries << "," << both.signal << "," << both.background << "\n";
      }
    }
  }
}

int TextColorForGain(const double v)
{
  return std::abs(v) > 1.9 ? kWhite : kBlack;
}
}

void MakeWidthStudyWindowAucComparisons()
{
  SetsPhenixStyle();
  gStyle->SetOptStat(0);

  const std::string report = "dataOutput/auauTightBDTValidation/model_validation_condor_20260511_171943";
  const std::string outDir = "dataOutput/jstg_slide_candidates/slide25_widthstudy_windows";
  gSystem->mkdir(outDir.c_str(), true);

  const auto aucRows = ReadAucRows(report + "/validation_auc_table.csv");
  const auto rankRows = ReadRankRows(report + "/validation_model_rankings.csv");

  WriteSummaryCsv(rankRows, outDir + "/widthstudy_window_auc_summary.csv");
  WriteCellGainCsv(aucRows, outDir + "/widthstudy_window_cell_gains.csv");

  const std::vector<std::string> windows = {"pt5to35", "pt10to35", "pt15to35"};
  const std::vector<std::string> windowLabels = {"5-35", "10-35", "15-35"};
  const std::vector<std::string> widthLabels = {"3x3 widths only", "Base + 3x3 widths"};
  const std::vector<std::string> widthProducts = {"centAsFeat3x3", "centAsFeatBase3x3"};

  {
    TCanvas c("c_window_summary", "c_window_summary", 1280, 760);
    c.SetLeftMargin(0.17);
    c.SetRightMargin(0.05);
    c.SetTopMargin(0.25);
    c.SetBottomMargin(0.16);

    TH2D frame("frame", ";Training window [GeV];ROC-AUC gain vs base widths [%]",
               3, 0.5, 3.5, 1, -1.5, 1.6);
    frame.SetStats(false);
    frame.GetXaxis()->SetBinLabel(1, "5-35");
    frame.GetXaxis()->SetBinLabel(2, "10-35");
    frame.GetXaxis()->SetBinLabel(3, "15-35");
    frame.GetXaxis()->SetTitleSize(0.043);
    frame.GetXaxis()->SetLabelSize(0.040);
    frame.GetYaxis()->SetTitleSize(0.041);
    frame.GetYaxis()->SetLabelSize(0.035);
    frame.GetYaxis()->SetTitleOffset(1.55);
    frame.Draw();

    TLine zero(0.5, 0, 3.5, 0);
    zero.SetLineColor(kGray + 2);
    zero.SetLineStyle(2);
    zero.SetLineWidth(2);
    zero.Draw();

    TGraph g33;
    TGraph gboth;
    TGraph gmean;
    for (int iw = 0; iw < 3; ++iw)
    {
      const std::string baseProduct = "centAsFeat_" + windows[iw];
      const auto& base = rankRows.at(baseProduct);
      const auto& w33 = rankRows.at("centAsFeat3x3_" + windows[iw]);
      const auto& both = rankRows.at("centAsFeatBase3x3_" + windows[iw]);
      g33.SetPoint(iw, iw + 1 - 0.10, GainPercent(w33.aucInclusive, base.aucInclusive));
      gboth.SetPoint(iw, iw + 1 + 0.10, GainPercent(both.aucInclusive, base.aucInclusive));
      gmean.SetPoint(iw, iw + 1 + 0.10, GainPercent(both.aucPtCentMean, base.aucPtCentMean));
    }
    g33.SetMarkerStyle(24);
    g33.SetMarkerSize(1.8);
    g33.SetMarkerColor(kAzure + 2);
    g33.SetLineColor(kAzure + 2);
    gboth.SetMarkerStyle(20);
    gboth.SetMarkerSize(1.8);
    gboth.SetMarkerColor(kRed + 1);
    gboth.SetLineColor(kRed + 1);
    gmean.SetMarkerStyle(25);
    gmean.SetMarkerSize(1.6);
    gmean.SetMarkerColor(kRed + 1);
    gmean.SetLineColor(kRed + 1);
    gmean.SetLineStyle(2);
    g33.Draw("P");
    gboth.Draw("P");
    gmean.Draw("P");

    TLatex text;
    text.SetTextFont(42);
    text.SetTextAlign(22);
    text.SetTextSize(0.033);
    for (int iw = 0; iw < 3; ++iw)
    {
      double x = 0, y = 0;
      gboth.GetPoint(iw, x, y);
      text.DrawLatex(x, y + 0.18, Form("%+.2f%%", y));
      g33.GetPoint(iw, x, y);
      text.DrawLatex(x, y - 0.18, Form("%+.2f%%", y));
    }

    TLegend leg(0.59, 0.70, 0.92, 0.82);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetTextFont(42);
    leg.SetTextSize(0.028);
    leg.AddEntry(&g33, "3#times3 only, inclusive", "p");
    leg.AddEntry(&gboth, "Base + 3#times3, inclusive", "p");
    leg.AddEntry(&gmean, "Base + 3#times3, mean E_{T}#timescent", "p");
    leg.Draw();

    TLatex ndc;
    ndc.SetNDC();
    ndc.SetTextFont(42);
    ndc.SetTextAlign(22);
    ndc.SetTextSize(0.043);
    ndc.DrawLatex(0.50, 0.960, "AUC gain from 3#times3 shower-width information");
    ndc.SetTextSize(0.027);
    ndc.DrawLatex(0.50, 0.910, "Each point compares to the base-width model trained in the same p_{T} window");
    DrawSphenixLabel(0.15, 0.870);

    const std::string png = outDir + "/widthstudy_window_auc_gain_summary.png";
    c.SaveAs(png.c_str());
  }

  {
    const std::vector<std::string> cents = {"0_20", "20_50", "50_80"};
    const std::vector<std::string> centLabels = {"0-20%", "20-50%", "50-80%"};
    const std::vector<std::string> pts = {"6_10", "10_15", "15_20", "20_25", "25_35"};
    const std::vector<std::string> ptLabels = {"6-10", "10-15", "15-20", "20-25", "25-35"};

    TCanvas c("c_base3x3_gain", "c_base3x3_gain", 1550, 760);
    c.SetLeftMargin(0.14);
    c.SetRightMargin(0.14);
    c.SetTopMargin(0.24);
    c.SetBottomMargin(0.17);

    const double centGap = 0.45;
    const double groupWidth = 5.0;
    const Double_t xEdges[] = {
      0.0, 1.0, 2.0, 3.0, 4.0, 5.0,
      5.0 + centGap,
      6.0 + centGap, 7.0 + centGap, 8.0 + centGap, 9.0 + centGap, 10.0 + centGap,
      10.0 + 2.0 * centGap,
      11.0 + 2.0 * centGap, 12.0 + 2.0 * centGap, 13.0 + 2.0 * centGap,
      14.0 + 2.0 * centGap, 15.0 + 2.0 * centGap};
    auto groupStart = [centGap, groupWidth](const int ic) { return ic * (groupWidth + centGap); };
    auto xBinFor = [](const int ic, const int ip) { return 1 + ic * 6 + ip; };
    auto xCenterFor = [groupStart](const int ic, const int ip) { return groupStart(ic) + ip + 0.5; };

    auto* h = new TH2D("h_base3x3_gain",
                       ";Photon candidate E_{T} validation bin within each centrality group [GeV];Training p_{T} window",
                       17, xEdges, 3, 0, 3);
    h->SetDirectory(nullptr);
    h->SetStats(false);
    h->SetMinimum(-0.5);
    h->SetMaximum(2.0);
    h->GetXaxis()->SetTitleSize(0.040);
    h->GetXaxis()->SetLabelSize(0.026);
    h->GetXaxis()->SetTitleOffset(1.18);
    h->GetYaxis()->SetTitleSize(0.041);
    h->GetYaxis()->SetLabelSize(0.034);
    h->GetYaxis()->SetTitleOffset(2.50);
    h->GetZaxis()->SetTitle("ROC-AUC gain vs base widths [%]");
    h->GetZaxis()->SetTitleSize(0.032);
    h->GetZaxis()->SetTitleOffset(1.25);
    h->GetZaxis()->SetLabelSize(0.029);

    const Int_t nStops = 5;
    Double_t stops[nStops] = {0.00, 0.20, 0.38, 0.68, 1.00};
    Double_t red[nStops] = {0.129, 0.965, 0.996, 0.878, 0.698};
    Double_t green[nStops] = {0.400, 0.965, 0.878, 0.486, 0.094};
    Double_t blue[nStops] = {0.675, 0.965, 0.565, 0.286, 0.168};
    TColor::CreateGradientColorTable(nStops, stops, red, green, blue, 255);
    gStyle->SetNumberContours(255);

    std::vector<std::vector<std::string>> cellText(3, std::vector<std::string>(15, ""));
    std::vector<std::vector<double>> cellGain(3, std::vector<double>(15, NAN));
    for (int iw = 0; iw < 3; ++iw)
    {
      h->GetYaxis()->SetBinLabel(3 - iw, (windowLabels[iw] + " GeV").c_str());
      for (int ic = 0; ic < 3; ++ic)
      {
        for (int ip = 0; ip < 5; ++ip)
        {
          const int col = ic * 5 + ip;
          const int xbin = xBinFor(ic, ip);
          h->GetXaxis()->SetBinLabel(xbin, ptLabels[ip].c_str());
          const auto& base = GetAuc(aucRows, "centAsFeat_" + windows[iw], cents[ic], pts[ip]);
          const auto& both = GetAuc(aucRows, "centAsFeatBase3x3_" + windows[iw], cents[ic], pts[ip]);
          const double gain = GainPercent(both.auc, base.auc);
          if (std::isfinite(gain))
          {
            h->SetBinContent(xbin, 3 - iw, gain);
            cellGain[iw][col] = gain;
            std::ostringstream ss;
            ss << std::showpos << std::fixed << std::setprecision(1) << gain << "%";
            cellText[iw][col] = ss.str();
          }
          else
          {
            h->SetBinContent(xbin, 3 - iw, -999);
            cellText[iw][col] = "";
          }
        }
      }
    }

    h->SetContour(255);
    h->Draw("COLZ");

    TBox gapBox;
    gapBox.SetFillColor(kWhite);
    gapBox.SetLineColor(kWhite);
    gapBox.DrawBox(groupStart(1) - centGap, 0.0, groupStart(1), 3.0);
    gapBox.DrawBox(groupStart(2) - centGap, 0.0, groupStart(2), 3.0);

    for (const double xmid : {groupStart(1) - 0.5 * centGap, groupStart(2) - 0.5 * centGap})
    {
      TLine guide(xmid, 0.0, xmid, 3.0);
      guide.SetLineWidth(2);
      guide.SetLineColor(kGray + 1);
      guide.SetLineStyle(2);
      guide.Draw();
    }

    TLatex text;
    text.SetTextFont(42);
    text.SetTextAlign(22);
    text.SetTextSize(0.026);
    for (int iw = 0; iw < 3; ++iw)
    {
      for (int col = 0; col < 15; ++col)
      {
        if (!std::isfinite(cellGain[iw][col])) continue;
        text.SetTextColor(TextColorForGain(cellGain[iw][col]));
        const int ic = col / 5;
        const int ip = col % 5;
        text.DrawLatex(xCenterFor(ic, ip), 2.5 - iw, cellText[iw][col].c_str());
      }
    }
    text.SetTextColor(kBlack);

    for (int y = 0; y <= 3; ++y)
    {
      for (int ic = 0; ic < 3; ++ic)
      {
        const double x0 = groupStart(ic);
        TLine line(x0, y, x0 + groupWidth, y);
        line.SetLineWidth(y == 0 || y == 3 ? 2 : 1);
        line.SetLineColor(kGray + 2);
        line.Draw();
      }
    }
    for (int ic = 0; ic < 3; ++ic)
    {
      const double x0 = groupStart(ic);
      for (int ip = 0; ip <= 5; ++ip)
      {
        TLine line(x0 + ip, 0, x0 + ip, 3);
        line.SetLineWidth(ip == 0 || ip == 5 ? 2 : 1);
        line.SetLineColor(kGray + 2);
        line.Draw();
      }
    }

    TLatex header;
    header.SetTextFont(42);
    header.SetTextAlign(22);
    header.SetTextSize(0.034);
    for (int ic = 0; ic < 3; ++ic)
    {
      const double x0 = groupStart(ic);
      const double x1 = x0 + groupWidth;
      header.DrawLatex(0.5 * (x0 + x1), 3.24, centLabels[ic].c_str());
    }

    TLatex ndc;
    ndc.SetNDC();
    ndc.SetTextFont(42);
    ndc.SetTextAlign(22);
    ndc.SetTextSize(0.043);
    ndc.DrawLatex(0.50, 0.965, "AUC gain from combining base and 3#times3 width inputs");
    ndc.SetTextSize(0.027);
    ndc.DrawLatex(0.50, 0.912, "Percent change relative to the base-width model in the same training window and validation bin");
    DrawSphenixLabel(0.18, 0.875);

    const std::string png = outDir + "/widthstudy_base_plus_3x3_gain_by_window.png";
    c.SaveAs(png.c_str());
  }

  std::cout << "[MakeWidthStudyWindowAucComparisons] wrote outputs under " << outDir << "\n";
}
