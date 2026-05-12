#include "sPhenixStyle.C"

#include <TBox.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TH2D.h>
#include <TLatex.h>
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
struct Row
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

using Key = std::tuple<std::string, std::string, std::string>;

std::map<Key, Row> ReadRows(const std::string& path)
{
  std::ifstream in(path);
  std::map<Key, Row> rows;
  std::string line;
  std::getline(in, line);
  while (std::getline(in, line))
  {
    const auto cols = SplitCsvLine(line);
    if (cols.size() < 7) continue;
    Row r;
    r.auc = std::stod(cols[3]);
    r.entries = std::stoi(cols[4]);
    r.signal = std::stoi(cols[5]);
    r.background = std::stoi(cols[6]);
    rows[{cols[0], cols[1], cols[2]}] = r;
  }
  return rows;
}

const Row& Get(const std::map<Key, Row>& rows,
               const std::string& product,
               const std::string& cent,
               const std::string& pt)
{
  return rows.at({product, cent, pt});
}

void DrawSphenixLabel(const double x, const double y)
{
  TLatex text;
  text.SetNDC();
  text.SetTextAlign(11);
  text.SetTextFont(42);
  text.SetTextSize(0.035);
  text.DrawLatex(x, y, "#it{#bf{sPHENIX}} Internal");
  text.SetTextSize(0.028);
  text.DrawLatex(x, y - 0.045, "Pythia overlay, #sqrt{s_{NN}} = 200 GeV");
}

void WriteGainCsv(const std::map<Key, Row>& rows, const std::string& outPath)
{
  const std::vector<std::string> products = {"centDepBDTs_pt5to40", "centDepFineBDTs_pt5to40"};
  const std::vector<std::string> cents = {"0_20", "20_50", "50_80"};
  const std::vector<std::string> pts = {"6_10", "10_15", "15_20", "20_25", "25_35"};
  std::ofstream out(outPath);
  out << "model,centrality,pt_bin,baseline_auc,dedicated_auc,auc_gain_percent,entries,signal_entries,background_entries\n";
  for (const auto& product : products)
  {
    for (const auto& cent : cents)
    {
      for (const auto& pt : pts)
      {
        const auto& base = Get(rows, "centAsFeat_pt5to40", cent, pt);
        const auto& val = Get(rows, product, cent, pt);
        const double gain = 100.0 * (val.auc - base.auc) / base.auc;
        out << product << "," << cent << "," << pt << ","
            << std::setprecision(10) << base.auc << "," << val.auc << ","
            << gain << "," << val.entries << "," << val.signal << "," << val.background << "\n";
      }
    }
  }
}
}

void MakeCentBinnedAucGainHeatmap()
{
  SetsPhenixStyle();
  gStyle->SetOptStat(0);
  const std::string report = "dataOutput/auauTightBDTValidation/model_validation_condor_20260509_192942/validation_auc_table.csv";
  const std::string outDir = "dataOutput/jstg_slide_candidates/slide17_ptcent_auc_gain";
  gSystem->mkdir(outDir.c_str(), true);
  const auto rows = ReadRows(report);

  const std::vector<std::string> cent = {"0_20", "20_50", "50_80"};
  const std::vector<std::string> centLabel = {"0-20%", "20-50%", "50-80%"};
  const std::vector<std::string> pt = {"6_10", "10_15", "15_20", "20_25", "25_35"};
  const std::vector<std::string> ptLabel = {"6-10", "10-15", "15-20", "20-25", "25-35"};
  const std::vector<std::string> products = {"centDepBDTs_pt5to40", "centDepFineBDTs_pt5to40"};
  const std::vector<std::string> labels = {
      "3 centrality-bin BDTs",
      "7 centrality-bin BDTs",
  };

  const Int_t nStops = 5;
  Double_t stops[nStops] = {0.00, 1.0 / 23.0, 0.36, 0.68, 1.00};
  Double_t red[nStops] = {0.129, 0.965, 0.996, 0.878, 0.698};
  Double_t green[nStops] = {0.400, 0.965, 0.878, 0.486, 0.094};
  Double_t blue[nStops] = {0.675, 0.965, 0.565, 0.286, 0.168};
  TColor::CreateGradientColorTable(nStops, stops, red, green, blue, 255);
  gStyle->SetNumberContours(255);

  std::vector<std::vector<std::string>> cellText(2, std::vector<std::string>(15));
  std::vector<std::vector<double>> deltaText(2, std::vector<double>(15, 0.0));
  for (int ic = 0; ic < 3; ++ic)
  {
    for (int ip = 0; ip < 5; ++ip)
    {
      const int xbin = ic * 5 + ip + 1;
      const double base = Get(rows, "centAsFeat_pt5to40", cent[ic], pt[ip]).auc;
      for (int row = 0; row < 2; ++row)
      {
        const double val = Get(rows, products[row], cent[ic], pt[ip]).auc;
        const double delta = (val - base) / base * 100.0;
        deltaText[row][xbin - 1] = delta;
        std::ostringstream ss;
        ss << std::showpos << std::fixed << std::setprecision(1) << delta << "%";
        cellText[row][xbin - 1] = ss.str();
      }
    }
  }

  TCanvas c("c_centbinned_auc_gain", "c_centbinned_auc_gain", 1550, 760);
  c.SetLeftMargin(0.18);
  c.SetRightMargin(0.14);
  c.SetTopMargin(0.24);
  c.SetBottomMargin(0.17);

  auto* h = new TH2D("h_centbinned_delta", ";Photon candidate E_{T} bin within each centrality group [GeV];",
                     15, 0, 15, 2, 0, 2);
  h->SetDirectory(nullptr);
  h->SetStats(false);
  h->SetMinimum(-1.0);
  h->SetMaximum(22.0);
  h->GetXaxis()->SetTitleSize(0.045);
  h->GetXaxis()->SetLabelSize(0.027);
  h->GetXaxis()->SetTitleOffset(1.10);
  h->GetYaxis()->SetTitleSize(0.045);
  h->GetYaxis()->SetLabelSize(0.036);
  h->GetYaxis()->SetTitleOffset(2.85);
  h->GetZaxis()->SetTitle("Percent gain wrt centrality-as-input ROC-AUC [%]");
  h->GetZaxis()->SetTitleSize(0.034);
  h->GetZaxis()->SetTitleOffset(1.25);
  h->GetZaxis()->SetLabelSize(0.030);
  h->GetYaxis()->SetBinLabel(2, labels[0].c_str());
  h->GetYaxis()->SetBinLabel(1, labels[1].c_str());
  for (int ic = 0; ic < 3; ++ic)
  {
    for (int ip = 0; ip < 5; ++ip)
    {
      const int globalCol = ic * 5 + ip;
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
    for (int col = 0; col < 15; ++col)
    {
      const double v = deltaText[row][col];
      text.SetTextSize(0.027);
      text.SetTextColor(v > 8.0 ? kWhite : kBlack);
      text.DrawLatex(col + 0.5, 1.5 - row, cellText[row][col].c_str());
    }
  }
  text.SetTextColor(kBlack);

  TLine rowLine(0, 1, 15, 1);
  rowLine.SetLineWidth(2);
  rowLine.SetLineColor(kGray + 2);
  rowLine.Draw();

  TLine outerTop(0, 2, 15, 2);
  TLine outerBottom(0, 0, 15, 0);
  TLine outerLeft(0, 0, 0, 2);
  TLine outerRight(15, 0, 15, 2);
  for (auto* l : {&outerTop, &outerBottom, &outerLeft, &outerRight})
  {
    l->SetLineWidth(2);
    l->SetLineColor(kGray + 2);
    l->Draw();
  }

  for (double x : {5.0, 10.0})
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
    const double x0 = ic * 5.0;
    const double x1 = x0 + 5.0;
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
  ndc.DrawLatex(0.50, 0.965, "Percent gain in ROC-AUC from centrality-binned BDTs");
  DrawSphenixLabel(0.20, 0.905);

  const std::string png = outDir + "/centbinned_auc_gain_vs_centinput_heatmap.png";
  c.SaveAs(png.c_str());
  WriteGainCsv(rows, outDir + "/centbinned_auc_gain_vs_centinput_table.csv");
  std::cout << "[MakeCentBinnedAucGainHeatmap] wrote " << png << "\n";
}
