R__ADD_INCLUDE_PATH(/Users/patsfan753/Desktop/ThesisAnalysis/ppg12codeGit/plotting)
#include "/Users/patsfan753/Desktop/ThesisAnalysis/ppg12codeGit/plotting/plotcommon.h"

#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TLine.h>
#include <TMath.h>
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
struct Row
{
  std::string sample;
  double x = 0;
  double sdcc = 0;
  double sdcc_err = 0;
  double fit = 0;
  double sdcc_over_fit = 0;
  double sdcc_over_fit_err = 0;
  double current = 0;
  double current_err = 0;
  double current_over_fit = 0;
  double current_over_fit_err = 0;
  double current_over_sdcc = 0;
};

std::vector<std::string> split_csv_line(const std::string& line)
{
  std::vector<std::string> out;
  std::string field;
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
      out.push_back(field);
      field.clear();
      continue;
    }
    field.push_back(ch);
  }
  out.push_back(field);
  return out;
}

double to_double(const std::vector<std::string>& vals,
                 const std::map<std::string, int>& col,
                 const std::string& key)
{
  auto it = col.find(key);
  if (it == col.end() || it->second < 0 || it->second >= static_cast<int>(vals.size())) return 0.0;
  return std::atof(vals[it->second].c_str());
}

std::vector<Row> read_rows(const std::string& csv_path)
{
  std::ifstream in(csv_path);
  if (!in)
  {
    std::cerr << "ERROR: cannot open " << csv_path << std::endl;
    return {};
  }

  std::string line;
  if (!std::getline(in, line)) return {};
  auto headers = split_csv_line(line);
  std::map<std::string, int> col;
  for (int i = 0; i < static_cast<int>(headers.size()); ++i) col[headers[i]] = i;

  std::vector<Row> rows;
  while (std::getline(in, line))
  {
    if (line.empty()) continue;
    auto vals = split_csv_line(line);
    Row r;
    r.sample = vals[col["sample"]];
    r.x = to_double(vals, col, "bin_center");
    r.sdcc = to_double(vals, col, "ppg12_sdcc_value");
    r.sdcc_err = to_double(vals, col, "ppg12_sdcc_error");
    r.fit = to_double(vals, col, "ppg12_fit_value");
    r.sdcc_over_fit = to_double(vals, col, "ppg12_sdcc_over_fit");
    r.sdcc_over_fit_err = to_double(vals, col, "ppg12_sdcc_over_fit_error");
    r.current = to_double(vals, col, "current_value");
    r.current_err = to_double(vals, col, "current_error");
    r.current_over_fit = to_double(vals, col, "current_over_fit");
    r.current_over_fit_err = to_double(vals, col, "current_over_fit_error");
    r.current_over_sdcc = to_double(vals, col, "current_over_sdcc");
    rows.push_back(r);
  }
  return rows;
}

std::vector<Row> by_sample(const std::vector<Row>& rows, const std::string& sample)
{
  std::vector<Row> out;
  for (const auto& r : rows)
  {
    if (r.sample == sample) out.push_back(r);
  }
  std::sort(out.begin(), out.end(), [](const Row& a, const Row& b) { return a.x < b.x; });
  return out;
}

TGraphErrors* make_graph(const std::vector<Row>& rows,
                         const std::string& name,
                         double Row::*y,
                         double Row::*ey)
{
  auto* g = new TGraphErrors();
  g->SetName(name.c_str());
  for (int i = 0; i < static_cast<int>(rows.size()); ++i)
  {
    g->SetPoint(i, rows[i].x, rows[i].*y);
    g->SetPointError(i, 0.0, rows[i].*ey);
  }
  return g;
}

TGraphErrors* make_fit_graph(const std::vector<Row>& rows)
{
  std::vector<Row> sorted = rows;
  std::sort(sorted.begin(), sorted.end(), [](const Row& a, const Row& b) { return a.x < b.x; });
  auto* g = new TGraphErrors();
  g->SetName("g_ppg12_fit");
  for (int i = 0; i < static_cast<int>(sorted.size()); ++i)
  {
    g->SetPoint(i, sorted[i].x, sorted[i].fit);
    g->SetPointError(i, 0.0, 0.0);
  }
  return g;
}

void style_sample_graph(TGraphErrors* g, int color, int marker)
{
  g->SetMarkerStyle(marker);
  g->SetMarkerColor(color);
  g->SetLineColor(color);
  g->SetMarkerSize(0.78);
  g->SetLineWidth(1);
}

void style_ratio_graph(TGraphErrors* g, int marker, int fill_color)
{
  g->SetMarkerStyle(marker);
  g->SetMarkerColor(kBlack);
  g->SetLineColor(kBlack);
  g->SetMarkerSize(0.72);
  g->SetLineWidth(1);
  if (fill_color == 0) g->SetFillStyle(0);
}

}  // namespace

void MakeTHE76PhotonFig5OverlayRootStyle(const char* csv_override = "",
                                         const char* png_override = "",
                                         const char* manifest_override = "",
                                         const char* comparison_scope_override = "")
{
  const std::string base =
      "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/ppg12Parity/"
      "the76_ppg12_parity_full_20260701_003024/strict_stitched_photonjet";
  const std::string csv_path = std::string(csv_override).empty()
                                   ? base + "/photon_data_over_fit_sdcc_vs_current_overlay_corrected_source_scope_points.csv"
                                   : std::string(csv_override);
  const std::string png_path = std::string(png_override).empty()
                                   ? base + "/photon_data_over_fit_sdcc_vs_current_overlay_root_ppg12_style_corrected_source_scope.png"
                                   : std::string(png_override);
  const std::string manifest_path = std::string(manifest_override).empty()
                                        ? base + "/photon_data_over_fit_sdcc_vs_current_overlay_root_ppg12_style_corrected_source_scope_manifest.json"
                                        : std::string(manifest_override);
  const std::string comparison_scope = std::string(comparison_scope_override).empty()
                                           ? "Corrected photon Fig.5 source-scope comparison: PPG12 SDCC open circles vs July1 current strict no-scale photon stitch points"
                                           : std::string(comparison_scope_override);

  auto rows = read_rows(csv_path);
  if (rows.empty())
  {
    std::cerr << "ERROR: no rows read from " << csv_path << std::endl;
    return;
  }

  init_plot();
  gStyle->SetOptStat(0);

  const std::vector<std::string> samples = {"photon5", "photon10", "photon20"};
  const std::map<std::string, int> colors = {
      {"photon5", kPink + 5},
      {"photon10", kGreen - 2},
      {"photon20", kAzure + 7},
  };

  std::map<std::string, TGraphErrors*> g_sdcc_top;
  std::map<std::string, TGraphErrors*> g_current_top;
  std::map<std::string, TGraphErrors*> g_sdcc_ratio;
  std::map<std::string, TGraphErrors*> g_current_ratio;

  double current_over_sdcc_min = 1e9;
  double current_over_sdcc_max = -1e9;
  std::vector<double> current_over_sdcc;
  std::vector<double> fit_diffs;

  for (const auto& s : samples)
  {
    auto rs = by_sample(rows, s);
    g_sdcc_top[s] = make_graph(rs, ("g_sdcc_top_" + s).c_str(), &Row::sdcc, &Row::sdcc_err);
    g_current_top[s] = make_graph(rs, ("g_current_top_" + s).c_str(), &Row::current, &Row::current_err);
    g_sdcc_ratio[s] = make_graph(rs, ("g_sdcc_ratio_" + s).c_str(), &Row::sdcc_over_fit, &Row::sdcc_over_fit_err);
    g_current_ratio[s] = make_graph(rs, ("g_current_ratio_" + s).c_str(), &Row::current_over_fit, &Row::current_over_fit_err);
    style_sample_graph(g_sdcc_top[s], colors.at(s), 24);
    style_sample_graph(g_current_top[s], colors.at(s), 20);
    style_ratio_graph(g_sdcc_ratio[s], 24, 0);
    style_ratio_graph(g_current_ratio[s], 20, kBlack);

    for (const auto& r : rs)
    {
      current_over_sdcc_min = std::min(current_over_sdcc_min, r.current_over_sdcc);
      current_over_sdcc_max = std::max(current_over_sdcc_max, r.current_over_sdcc);
      current_over_sdcc.push_back(r.current_over_sdcc);
      fit_diffs.push_back(r.current_over_fit - r.sdcc_over_fit);
    }
  }

  auto* g_fit = make_fit_graph(rows);
  g_fit->SetLineColor(kRed - 4);
  g_fit->SetLineWidth(2);
  g_fit->SetMarkerSize(0);

  TCanvas* c4 = new TCanvas("can_the76_ppg12_root_overlay", "", 800, 889);
  c4->SetCanvasSize(800, 889);
  c4->SetWindowSize(800 + (800 - c4->GetWw()), 889 + (889 - c4->GetWh()));
  c4->Divide(1, 2);

  TPad* p1 = static_cast<TPad*>(c4->cd(1));
  p1->SetPad(0, 0.4, 1, 1);
  p1->SetTopMargin(0.05);
  p1->SetLeftMargin(0.13);
  p1->SetBottomMargin(0.002);
  p1->SetRightMargin(0.08);
  p1->SetLogy();
  p1->SetTickx(1);
  p1->SetTicky(1);

  frame_et_rec->SetYTitle("d#sigma / d#it{E}_{T}^{#gamma} [pb / GeV]");
  frame_et_rec->GetYaxis()->SetRangeUser(0.025, 5.0e4);
  frame_et_rec->GetXaxis()->SetRangeUser(10, 40);
  frame_et_rec->GetXaxis()->SetTitleOffset(1.05);
  frame_et_rec->GetYaxis()->SetTitleOffset(1.05);
  frame_et_rec->GetYaxis()->SetTitleSize(0.053);
  frame_et_rec->GetXaxis()->SetLabelSize(0.050);
  frame_et_rec->GetYaxis()->SetLabelSize(0.050);
  frame_et_rec->GetXaxis()->SetLabelOffset(2);
  frame_et_rec->GetXaxis()->SetNdivisions(505);
  frame_et_rec->SetStats(false);
  frame_et_rec->Draw("axis");

  for (const auto& s : samples) g_sdcc_top[s]->Draw("P same");
  for (const auto& s : samples) g_current_top[s]->Draw("P same");
  g_fit->Draw("L same");

  myText(0.5, 0.9, 1, strleg1.c_str(), 0.05);
  myText(0.5, 0.85, 1, strleg2.c_str(), 0.05);
  myText(0.5, 0.80, 1, strMC.c_str(), 0.05);

  myText(0.20, 0.25, 1, "sample", 0.038);
  myMarkerLineText(0.20, 0.20, 1, kPink + 5, 20, kPink + 5, 1, "photon5", 0.040, true);
  myMarkerLineText(0.20, 0.15, 1, kGreen - 2, 20, kGreen - 2, 1, "photon10", 0.040, true);
  myMarkerLineText(0.20, 0.10, 1, kAzure + 7, 20, kAzure + 7, 1, "photon20", 0.040, true);
  myText(0.42, 0.25, 1, "source", 0.038);
  myMarkerLineText(0.42, 0.20, 1, kBlack, 24, kBlack, 1, "PPG12 SDCC", 0.040, true);
  myMarkerLineText(0.42, 0.15, 1, kBlack, 20, kBlack, 1, "Current output", 0.040, true);
  myMarkerLineText(0.42, 0.10, 0, kRed - 4, 0, kRed - 4, 1, "PPG12 fit", 0.040, true);

  TPad* p2 = static_cast<TPad*>(c4->cd(2));
  p2->SetPad(0, 0, 1, 0.4);
  p2->SetTopMargin(0.023);
  p2->SetLeftMargin(0.13);
  p2->SetBottomMargin(0.25);
  p2->SetRightMargin(0.08);
  p2->SetTickx(1);
  p2->SetTicky(1);

  frame_et_truth->SetYTitle("Data / Fit");
  frame_et_truth->SetXTitle("Leading #it{E}_{T}^{#gamma} [GeV]");
  frame_et_truth->GetYaxis()->SetNdivisions(506);
  frame_et_truth->GetYaxis()->SetRangeUser(0.85, 1.15);
  frame_et_truth->GetXaxis()->SetRangeUser(10, 40);
  frame_et_truth->GetXaxis()->SetTitleOffset(frame_et_rec->GetXaxis()->GetTitleOffset() * 4 / 6. * 1.4);
  frame_et_truth->GetYaxis()->SetTitleOffset(frame_et_rec->GetYaxis()->GetTitleOffset() * 4 / 6.);
  frame_et_truth->GetYaxis()->SetLabelOffset(frame_et_rec->GetYaxis()->GetLabelOffset() * 4 / 6.);
  frame_et_truth->GetXaxis()->SetLabelSize(frame_et_rec->GetXaxis()->GetLabelSize() * 6 / 4.);
  frame_et_truth->GetYaxis()->SetLabelSize(frame_et_rec->GetYaxis()->GetLabelSize() * 6 / 4.);
  frame_et_truth->GetXaxis()->SetTitleSize(frame_et_rec->GetXaxis()->GetTitleSize() * 6 / 4. * 1.2);
  frame_et_truth->GetYaxis()->SetTitleSize(frame_et_rec->GetYaxis()->GetTitleSize() * 6 / 4.);
  frame_et_truth->GetXaxis()->SetNdivisions(505);
  frame_et_truth->SetStats(false);
  frame_et_truth->Draw("axis");

  for (const auto& s : samples) g_current_ratio[s]->Draw("PE same");
  for (const auto& s : samples) g_sdcc_ratio[s]->Draw("PE same");
  lineone->Draw("L");

  c4->SaveAs(png_path.c_str());

  double median = 0;
  if (!current_over_sdcc.empty())
  {
    auto tmp = current_over_sdcc;
    std::sort(tmp.begin(), tmp.end());
    const size_t mid = tmp.size() / 2;
    median = (tmp.size() % 2) ? tmp[mid] : 0.5 * (tmp[mid - 1] + tmp[mid]);
  }
  double rms = 0;
  for (double d : fit_diffs) rms += d * d;
  if (!fit_diffs.empty()) rms = std::sqrt(rms / fit_diffs.size());

  std::ofstream out(manifest_path);
  out << std::setprecision(12);
  out << "{\n";
  out << "  \"status\": \"ok_root_ppg12_style_overlay_corrected_source_scope\",\n";
  out << "  \"plot_png\": \"" << png_path << "\",\n";
  out << "  \"points_csv\": \"" << csv_path << "\",\n";
  out << "  \"renderer\": \"ROOT with ppg12codeGit/plotting/plotcommon.h and sPhenixStyle.C\",\n";
  out << "  \"pixel_size\": [800, 889],\n";
  out << "  \"ppg12_checked_code\": \"ppg12codeGit/plotting/plot_combine_uncut.C lines 129-202 plus plotcommon.h/init_plot style helpers\",\n";
  out << "  \"style_contract\": \"Direct ROOT render using PPG12 plot_combine_uncut.C canvas, pad, margin, axis size, SetNdivisions, TLatex helper calls, and exact 10-40 GeV x range; no Matplotlib axis emulation.\",\n";
  out << "  \"normalization\": \"no global scale applied; PPG12 photon source convention raw*xsec/N\",\n";
  out << "  \"comparison_scope\": \"" << comparison_scope << "\",\n";
  out << "  \"n_matched_bins\": " << rows.size() << ",\n";
  out << "  \"current_over_sdcc_median\": " << median << ",\n";
  out << "  \"current_over_sdcc_min\": " << current_over_sdcc_min << ",\n";
  out << "  \"current_over_sdcc_max\": " << current_over_sdcc_max << ",\n";
  out << "  \"rms_current_minus_sdcc_over_fit\": " << rms << "\n";
  out << "}\n";
  out.close();
}
