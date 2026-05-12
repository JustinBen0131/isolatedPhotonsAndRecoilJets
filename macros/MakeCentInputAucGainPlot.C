#include "sPhenixStyle.C"

#include <TCanvas.h>
#include <TColor.h>
#include <TH2D.h>
#include <TLatex.h>
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
  text.SetTextSize(0.034);
  text.DrawLatex(x, y, "#it{#bf{sPHENIX}} Internal");
  text.SetTextSize(0.026);
  text.DrawLatex(x, y - 0.043, "Pythia overlay, #sqrt{s_{NN}} = 200 GeV");
}

void WriteGainCsv(const std::map<Key, Row>& rows, const std::string& outPath)
{
  const std::vector<std::string> cents = {"0_20", "20_50", "50_80"};
  const std::vector<std::string> pts = {"6_10", "10_15", "15_20", "20_25", "25_35"};

  std::ofstream out(outPath);
  out << "centrality,pt_bin,no_centrality_input_auc,centrality_as_input_auc,auc_gain_percent,entries,signal_entries,background_entries\n";
  for (const auto& cent : cents)
  {
    for (const auto& pt : pts)
    {
      const auto& base = Get(rows, "centINDcontrol_allRange", cent, pt);
      const auto& val = Get(rows, "centAsFeat_allRange", cent, pt);
      const double gain = 100.0 * (val.auc - base.auc) / base.auc;
      out << cent << "," << pt << ","
          << std::setprecision(10) << base.auc << "," << val.auc << ","
          << gain << "," << val.entries << "," << val.signal << "," << val.background << "\n";
    }
  }
}
}

void MakeCentInputAucGainPlot()
{
  SetsPhenixStyle();
  gStyle->SetOptStat(0);

  const std::string report = "dataOutput/auauTightBDTValidation/model_validation_condor_20260509_192942/validation_auc_table.csv";
  const std::string outDir = "dataOutput/jstg_slide_candidates/slide24_centinput_gain";
  gSystem->mkdir(outDir.c_str(), true);

  const auto rows = ReadRows(report);
  const std::vector<std::string> cent = {"0_20", "20_50", "50_80"};
  const std::vector<std::string> centLabel = {"0-20%", "20-50%", "50-80%"};
  const std::vector<std::string> pt = {"6_10", "10_15", "15_20", "20_25", "25_35"};
  const std::vector<std::string> ptLabel = {"6-10", "10-15", "15-20", "20-25", "25-35"};

  const Int_t nStops = 5;
  Double_t stops[nStops] = {0.00, 0.30, 0.50, 0.70, 1.00};
  Double_t red[nStops] = {0.129, 0.698, 0.965, 0.878, 0.698};
  Double_t green[nStops] = {0.400, 0.820, 0.965, 0.486, 0.094};
  Double_t blue[nStops] = {0.675, 0.933, 0.965, 0.286, 0.168};
  TColor::CreateGradientColorTable(nStops, stops, red, green, blue, 255);
  gStyle->SetNumberContours(255);

  TCanvas c("c_centinput_auc_gain", "c_centinput_auc_gain", 1280, 760);
  c.SetLeftMargin(0.13);
  c.SetRightMargin(0.20);
  c.SetTopMargin(0.25);
  c.SetBottomMargin(0.17);

  auto* h = new TH2D("h_centinput_auc_gain",
                     ";Photon candidate E_{T} validation bin [GeV];Centrality",
                     pt.size(), 0, pt.size(), cent.size(), 0, cent.size());
  h->SetDirectory(nullptr);
  h->SetStats(false);
  h->SetMinimum(-2.5);
  h->SetMaximum(2.5);
  h->GetXaxis()->SetTitleSize(0.042);
  h->GetXaxis()->SetLabelSize(0.034);
  h->GetXaxis()->SetTitleOffset(1.05);
  h->GetYaxis()->SetTitleSize(0.042);
  h->GetYaxis()->SetLabelSize(0.036);
  h->GetYaxis()->SetTitleOffset(1.25);
  h->GetZaxis()->SetTitle("ROC-AUC gain [%]");
  h->GetZaxis()->SetTitleSize(0.033);
  h->GetZaxis()->SetTitleOffset(1.15);
  h->GetZaxis()->SetLabelSize(0.030);

  std::vector<std::vector<double>> gains(cent.size(), std::vector<double>(pt.size(), 0.0));
  std::vector<std::vector<std::string>> labels(cent.size(), std::vector<std::string>(pt.size()));
  for (int ic = 0; ic < static_cast<int>(cent.size()); ++ic)
  {
    h->GetYaxis()->SetBinLabel(cent.size() - ic, centLabel[ic].c_str());
    for (int ip = 0; ip < static_cast<int>(pt.size()); ++ip)
    {
      h->GetXaxis()->SetBinLabel(ip + 1, ptLabel[ip].c_str());
      const double base = Get(rows, "centINDcontrol_allRange", cent[ic], pt[ip]).auc;
      const double val = Get(rows, "centAsFeat_allRange", cent[ic], pt[ip]).auc;
      const double gain = 100.0 * (val - base) / base;
      gains[ic][ip] = gain;
      std::ostringstream ss;
      ss << std::showpos << std::fixed << std::setprecision(1) << gain << "%";
      labels[ic][ip] = ss.str();
      h->SetBinContent(ip + 1, cent.size() - ic, gain);
    }
  }

  h->SetContour(255);
  h->Draw("COLZ");

  TLatex text;
  text.SetTextFont(42);
  text.SetTextAlign(22);
  for (int ic = 0; ic < static_cast<int>(cent.size()); ++ic)
  {
    const double y = cent.size() - ic - 0.5;
    for (int ip = 0; ip < static_cast<int>(pt.size()); ++ip)
    {
      const double gain = gains[ic][ip];
      text.SetTextSize(0.040);
      text.SetTextColor(std::abs(gain) > 1.4 ? kWhite : kBlack);
      text.DrawLatex(ip + 0.5, y, labels[ic][ip].c_str());
    }
  }
  text.SetTextColor(kBlack);

  for (int ix = 0; ix <= static_cast<int>(pt.size()); ++ix)
  {
    TLine line(ix, 0, ix, cent.size());
    line.SetLineColor(kGray + 2);
    line.SetLineWidth(ix == 0 || ix == static_cast<int>(pt.size()) ? 2 : 1);
    line.Draw();
  }
  for (int iy = 0; iy <= static_cast<int>(cent.size()); ++iy)
  {
    TLine line(0, iy, pt.size(), iy);
    line.SetLineColor(kGray + 2);
    line.SetLineWidth(iy == 0 || iy == static_cast<int>(cent.size()) ? 2 : 1);
    line.Draw();
  }

  TLatex ndc;
  ndc.SetNDC();
  ndc.SetTextFont(42);
  ndc.SetTextAlign(22);
  ndc.SetTextSize(0.043);
  ndc.DrawLatex(0.50, 0.960, "ROC-AUC gain from adding centrality as a BDT input");
  ndc.SetTextSize(0.027);
  ndc.DrawLatex(0.50, 0.910, "Single centrality-input model compared to a single no-centrality-input model");
  DrawSphenixLabel(0.15, 0.875);

  const std::string png = outDir + "/centinput_vs_nocent_auc_gain_heatmap.png";
  c.SaveAs(png.c_str());
  WriteGainCsv(rows, outDir + "/centinput_vs_nocent_auc_gain_table.csv");
  std::cout << "[MakeCentInputAucGainPlot] wrote " << png << "\n";
}
