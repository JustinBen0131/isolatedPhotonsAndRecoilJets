R__ADD_INCLUDE_PATH(/Users/patsfan753/Desktop/ThesisAnalysis/ppg12codeGit/plotting)
#include "/Users/patsfan753/Desktop/ThesisAnalysis/ppg12codeGit/plotting/plotcommon.h"

#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TLine.h>

#include <algorithm>
#include <cstdlib>
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
  double y = 0;
  double ey = 0;
  double fit = 0;
  double y_over_fit = 0;
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
    r.y = to_double(vals, col, "ppg12_sdcc_value");
    r.ey = to_double(vals, col, "ppg12_sdcc_error");
    r.fit = to_double(vals, col, "ppg12_fit_value");
    r.y_over_fit = to_double(vals, col, "ppg12_sdcc_over_fit");
    if (r.x <= 50.0) rows.push_back(r);
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

TGraphErrors* make_ratio_graph(const std::vector<Row>& rows, const std::string& name)
{
  auto* g = new TGraphErrors();
  g->SetName(name.c_str());
  for (int i = 0; i < static_cast<int>(rows.size()); ++i)
  {
    const double ey = rows[i].fit > 0.0 ? rows[i].ey / rows[i].fit : 0.0;
    g->SetPoint(i, rows[i].x, rows[i].y_over_fit);
    g->SetPointError(i, 0.0, ey);
  }
  return g;
}

TGraphErrors* make_fit_graph(const std::vector<Row>& rows)
{
  std::vector<Row> sorted = rows;
  std::sort(sorted.begin(), sorted.end(), [](const Row& a, const Row& b) { return a.x < b.x; });
  auto* g = new TGraphErrors();
  g->SetName("g_ppg12_jet_fit_source_scope");
  for (int i = 0; i < static_cast<int>(sorted.size()); ++i)
  {
    g->SetPoint(i, sorted[i].x, sorted[i].fit);
    g->SetPointError(i, 0.0, 0.0);
  }
  return g;
}

void style_top(TGraphErrors* g, int color)
{
  g->SetMarkerStyle(20);
  g->SetMarkerColor(color);
  g->SetLineColor(color);
  g->SetMarkerSize(0.78);
  g->SetLineWidth(1);
}

void style_ratio(TGraphErrors* g)
{
  g->SetMarkerStyle(20);
  g->SetMarkerColor(kBlack);
  g->SetLineColor(kBlack);
  g->SetMarkerSize(0.72);
  g->SetLineWidth(1);
}
}  // namespace

void MakeTHE94JetFig6PPG12SourceOnlyRootStyle(const char* csv_override = "",
                                              const char* png_override = "",
                                              const char* manifest_override = "")
{
  const std::string base =
      "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/ppg12Parity/"
      "the94_ppg12_inclusivejet_fig6_fixed_20260702_204032/ian_source_reproduction";
  const std::string input_csv =
      "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/ppg12Parity/"
      "the76_ppg12_parity_full_20260701_003024/strict_stitched_inclusivejet/"
      "jet_sdcc_over_current_shape_overlay_ppg12_jet8_xsecfix_points.csv";
  const std::string csv_path = std::string(csv_override).empty() ? input_csv : std::string(csv_override);
  const std::string png_path = std::string(png_override).empty()
                                   ? base + "/ppg12_fig6_inclusivejet_ian_source_reproduction_rootstyle.png"
                                   : std::string(png_override);
  const std::string manifest_path = std::string(manifest_override).empty()
                                        ? base + "/ppg12_fig6_inclusivejet_ian_source_reproduction_rootstyle_manifest.json"
                                        : std::string(manifest_override);

  gSystem->mkdir(base.c_str(), true);

  auto rows = read_rows(csv_path);
  if (rows.empty())
  {
    std::cerr << "ERROR: no rows read from " << csv_path << std::endl;
    return;
  }

  init_plot();
  gStyle->SetOptStat(0);

  const std::vector<std::string> samples = {"jet8", "jet12", "jet20", "jet30", "jet40"};
  const std::map<std::string, int> colors = {
      {"jet8", kPink + 5},
      {"jet12", kGreen - 2},
      {"jet20", kAzure + 7},
      {"jet30", kOrange + 7},
      {"jet40", kPink + 5},
  };

  std::map<std::string, TGraphErrors*> g_top;
  std::map<std::string, TGraphErrors*> g_ratio;
  for (const auto& s : samples)
  {
    auto rs = by_sample(rows, s);
    g_top[s] = make_graph(rs, ("g_ppg12_ian_top_" + s).c_str(), &Row::y, &Row::ey);
    g_ratio[s] = make_ratio_graph(rs, ("g_ppg12_ian_ratio_" + s).c_str());
    style_top(g_top[s], colors.at(s));
    style_ratio(g_ratio[s]);
  }

  auto* g_fit = make_fit_graph(rows);
  g_fit->SetLineColor(kRed - 4);
  g_fit->SetLineWidth(2);
  g_fit->SetMarkerSize(0);

  TCanvas* c4 = new TCanvas("can_the94_ppg12_jet_ian_source_only", "", 800, 889);
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

  frame_et_rec->SetYTitle("counts");
  frame_et_rec->GetYaxis()->SetRangeUser(5e3, 2.0e12);
  frame_et_rec->GetXaxis()->SetRangeUser(9, 50);
  frame_et_rec->GetXaxis()->SetTitleOffset(1.05);
  frame_et_rec->GetYaxis()->SetTitleOffset(1.05);
  frame_et_rec->GetYaxis()->SetTitleSize(0.053);
  frame_et_rec->GetXaxis()->SetLabelSize(0.050);
  frame_et_rec->GetYaxis()->SetLabelSize(0.050);
  frame_et_rec->GetXaxis()->SetLabelOffset(2);
  frame_et_rec->GetXaxis()->SetNdivisions(505);
  frame_et_rec->SetStats(false);
  frame_et_rec->Draw("axis");

  for (const auto& s : samples) g_top[s]->Draw("PE same");
  g_fit->Draw("L same");

  myText(0.47, 0.93, 1, strleg1.c_str(), 0.050);
  myText(0.47, 0.875, 1, strleg2.c_str(), 0.050);
  myText(0.47, 0.82, 1, strMC.c_str(), 0.050);

  myMarkerLineText(0.57, 0.75, 1, colors.at("jet8"), 20, colors.at("jet8"), 1, "jet8", 0.045, true);
  myMarkerLineText(0.57, 0.70, 1, colors.at("jet12"), 20, colors.at("jet12"), 1, "jet12", 0.045, true);
  myMarkerLineText(0.57, 0.65, 1, colors.at("jet20"), 20, colors.at("jet20"), 1, "jet20", 0.045, true);
  myMarkerLineText(0.57, 0.60, 1, colors.at("jet30"), 20, colors.at("jet30"), 1, "jet30", 0.045, true);
  myMarkerLineText(0.57, 0.55, 1, colors.at("jet40"), 20, colors.at("jet40"), 1, "jet40", 0.045, true);

  TPad* p2 = static_cast<TPad*>(c4->cd(2));
  p2->SetPad(0, 0, 1, 0.4);
  p2->SetTopMargin(0.023);
  p2->SetLeftMargin(0.13);
  p2->SetBottomMargin(0.25);
  p2->SetRightMargin(0.08);
  p2->SetTickx(1);
  p2->SetTicky(1);

  frame_et_truth->SetYTitle("MC / Fit");
  frame_et_truth->SetXTitle("Leading #it{p}_{T}^{jet} [GeV]");
  frame_et_truth->GetYaxis()->SetNdivisions(506);
  frame_et_truth->GetYaxis()->SetRangeUser(0.85, 1.15);
  frame_et_truth->GetXaxis()->SetRangeUser(9, 50);
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

  for (const auto& s : samples) g_ratio[s]->Draw("PE same");
  lineone->Draw("L");

  c4->SaveAs(png_path.c_str());

  std::ofstream out(manifest_path);
  out << std::setprecision(12);
  out << "{\n";
  out << "  \"status\": \"ok_ppg12_ian_source_reproduction\",\n";
  out << "  \"plot_png\": \"" << png_path << "\",\n";
  out << "  \"input_points_csv\": \"" << csv_path << "\",\n";
  out << "  \"renderer\": \"ROOT with ppg12codeGit/plotting/plotcommon.h\",\n";
  out << "  \"pixel_size\": [800, 889],\n";
  out << "  \"reference_scope\": \"PPG12 IAN/no-suffix source-scope inclusive-jet Fig.6 reproduction, not true period-combined THE-94 parity target\",\n";
  out << "  \"fit_source\": \"ppg12_fit_value column from the IAN/source-scope point table previously extracted from the PPG12 source reproduction\",\n";
  out << "  \"displayed_bins\": " << rows.size() << "\n";
  out << "}\n";
  out.close();
}
