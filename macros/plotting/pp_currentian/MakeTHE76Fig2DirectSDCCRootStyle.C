#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "TAxis.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TROOT.h"
#include "TStyle.h"

namespace {

std::vector<std::string> split_csv_line(const std::string& line)
{
  std::vector<std::string> out;
  std::string cur;
  bool in_quote = false;
  for (char ch : line)
  {
    if (ch == '"')
    {
      in_quote = !in_quote;
      continue;
    }
    if (ch == ',' && !in_quote)
    {
      out.push_back(cur);
      cur.clear();
    }
    else
    {
      cur.push_back(ch);
    }
  }
  out.push_back(cur);
  return out;
}

void set_axis_style(TH1F* h)
{
  h->GetXaxis()->SetTitle("E_{T}^{iso} Cutoff [GeV]");
  h->GetYaxis()->SetTitle("Fraction of Events");
  h->GetXaxis()->SetTitleOffset(1.12);
  h->GetYaxis()->SetTitleOffset(1.40);
  h->GetXaxis()->SetTitleSize(0.045);
  h->GetYaxis()->SetTitleSize(0.045);
  h->GetXaxis()->SetLabelSize(0.042);
  h->GetYaxis()->SetLabelSize(0.042);
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
  h->GetYaxis()->SetRangeUser(0.9, 1.05);
}

void draw_text(double x, double y, const char* text, double size = 0.04)
{
  TLatex latex;
  latex.SetNDC();
  latex.SetTextFont(42);
  latex.SetTextSize(size);
  latex.DrawLatex(x, y, text);
}

}  // namespace

void MakeTHE76Fig2DirectSDCCRootStyle(
    const char* csv_path =
        "dataOutput/ppg12Parity/the76_ppg12_fig24_photonjet_fix_20260702_014217/"
        "fig2_truth_iso_fraction_sdcc_direct/ppg12_fig2_direct_points_from_sdcc_root.csv",
    const char* out_png =
        "dataOutput/ppg12Parity/the76_ppg12_fig24_photonjet_fix_20260702_014217/"
        "fig2_truth_iso_fraction_sdcc_direct/ppg12_fig2_direct_sdcc_root_reproduction_root_ppg12_style.png")
{
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  std::map<std::string, TH1F*> hists;
  hists["pt10_15"] = new TH1F("h_eff_direct_0", "", 20, 0, 20);
  hists["pt15_20"] = new TH1F("h_eff_direct_1", "", 20, 0, 20);
  hists["pt25_30"] = new TH1F("h_eff_direct_2", "", 20, 0, 20);

  std::ifstream input(csv_path);
  if (!input)
  {
    std::cerr << "Could not open " << csv_path << std::endl;
    return;
  }

  std::string line;
  std::getline(input, line);  // header
  while (std::getline(input, line))
  {
    if (line.empty())
    {
      continue;
    }
    const auto cols = split_csv_line(line);
    if (cols.size() < 9)
    {
      continue;
    }
    const std::string series = cols[2];
    const int bin = std::stoi(cols[5]);
    const double fraction = std::stod(cols[8]);
    auto iter = hists.find(series);
    if (iter == hists.end())
    {
      continue;
    }
    iter->second->SetBinContent(bin, fraction);
    iter->second->SetBinError(bin, 0.0);
  }

  const std::vector<std::string> order = {"pt10_15", "pt15_20", "pt25_30"};
  const std::vector<int> colors = {kPink + 8, kSpring - 7, kAzure - 3};
  const std::vector<int> markers = {20, 21, 22};
  const std::vector<const char*> labels = {
      "Direct #gamma p_{T} [10, 15]",
      "Direct #gamma p_{T} [15, 20]",
      "Direct #gamma p_{T} [25, 30]"};

  TCanvas* can = new TCanvas("can", "", 800, 600);
  can->SetTopMargin(0.12);
  can->SetLeftMargin(0.15);
  can->SetBottomMargin(0.14);
  can->SetRightMargin(0.08);
  can->SetTicks(1, 1);

  for (std::size_t i = 0; i < order.size(); ++i)
  {
    TH1F* h = hists[order[i]];
    h->SetMarkerStyle(markers[i]);
    h->SetMarkerColor(colors[i]);
    h->SetLineColor(colors[i]);
    h->SetMarkerSize(0.9);
    set_axis_style(h);
    if (i == 0)
    {
      h->Draw("p");
    }
    else
    {
      h->Draw("p same");
    }
  }

  draw_text(0.065, 0.970, "#bf{#it{sPHENIX}} Internal", 0.04);
  draw_text(0.065, 0.920, "Photon+Jet Samples", 0.04);
  draw_text(0.150, 0.850, "Pythia, #sqrt{s}=200 GeV", 0.035);
  draw_text(0.150, 0.800, "vtx |z| < 30 cm, |#eta^{#gamma}| < 0.7", 0.035);

  TLegend* leg = new TLegend(0.60, 0.65, 0.90, 0.82);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.035);
  for (std::size_t i = 0; i < order.size(); ++i)
  {
    leg->AddEntry(hists[order[i]], labels[i], "p");
  }
  leg->Draw();

  can->SaveAs(out_png);
  delete can;
}
