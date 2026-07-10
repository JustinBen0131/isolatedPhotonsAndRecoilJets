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
  double ppg12 = 0;
  double ppg12_err = 0;
  double fit = 0;
  double ppg12_over_fit = 0;
  double ppg12_over_fit_err = 0;
  double current = 0;
  double current_err = 0;
  double current_over_fit = 0;
  double current_over_fit_err = 0;
  double current_over_ppg12 = 0;
  double display_factor = 1.0;
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
    r.ppg12 = to_double(vals, col, "ppg12_source_value");
    r.ppg12_err = to_double(vals, col, "ppg12_source_error");
    r.fit = to_double(vals, col, "ppg12_fit_value");
    r.ppg12_over_fit = to_double(vals, col, "ppg12_source_over_fit");
    r.ppg12_over_fit_err = to_double(vals, col, "ppg12_source_over_fit_error");
    r.current = to_double(vals, col, "current_display_value");
    r.current_err = to_double(vals, col, "current_display_error");
    r.current_over_fit = to_double(vals, col, "current_display_over_fit");
    r.current_over_fit_err = to_double(vals, col, "current_display_over_fit_error");
    r.current_over_ppg12 = to_double(vals, col, "current_display_over_ppg12_source");
    r.display_factor = to_double(vals, col, "source_scope_display_factor");
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

TGraphErrors* make_fit_graph(const std::vector<Row>& rows)
{
  std::vector<Row> sorted = rows;
  std::sort(sorted.begin(), sorted.end(), [](const Row& a, const Row& b) { return a.x < b.x; });
  auto* g = new TGraphErrors();
  g->SetName("g_ppg12_jet_fit_ian_overlay");
  for (int i = 0; i < static_cast<int>(sorted.size()); ++i)
  {
    g->SetPoint(i, sorted[i].x, sorted[i].fit);
    g->SetPointError(i, 0.0, 0.0);
  }
  return g;
}

void style_top(TGraphErrors* g, int color, int marker)
{
  g->SetMarkerStyle(marker);
  g->SetMarkerColor(color);
  g->SetLineColor(color);
  g->SetMarkerSize(marker == 24 ? 0.86 : 0.74);
  g->SetLineWidth(1);
  if (marker == 24) g->SetFillStyle(0);
}

void style_ratio(TGraphErrors* g, int marker)
{
  g->SetMarkerStyle(marker);
  g->SetMarkerColor(kBlack);
  g->SetLineColor(kBlack);
  g->SetMarkerSize(marker == 24 ? 0.80 : 0.66);
  g->SetLineWidth(1);
  if (marker == 24) g->SetFillStyle(0);
}
}  // namespace

void MakeTHE94JetFig6IANSourceOverlayRootStyle(const char* csv_override = "",
                                               const char* png_override = "",
                                               const char* manifest_override = "")
{
  const std::string base =
      "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/ppg12Parity/"
      "the94_ppg12_inclusivejet_fig6_fixed_20260702_204032/ian_source_overlay_the94_current";
  const std::string input_csv = base + "/inclusivejet_fig6_ian_source_overlay_the94_current_source_aligned_points.csv";
  const std::string csv_path = std::string(csv_override).empty() ? input_csv : std::string(csv_override);
  const std::string png_path = std::string(png_override).empty()
                                   ? base + "/inclusivejet_fig6_ian_source_overlay_the94_current_rootstyle.png"
                                   : std::string(png_override);
  const std::string manifest_path = std::string(manifest_override).empty()
                                        ? base + "/inclusivejet_fig6_ian_source_overlay_the94_current_rootstyle_manifest.json"
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

  std::map<std::string, TGraphErrors*> g_ppg12_top;
  std::map<std::string, TGraphErrors*> g_current_top;
  std::map<std::string, TGraphErrors*> g_ppg12_ratio;
  std::map<std::string, TGraphErrors*> g_current_ratio;

  for (const auto& s : samples)
  {
    auto rs = by_sample(rows, s);
    g_current_top[s] = make_graph(rs, ("g_the94_ian_top_" + s).c_str(), &Row::current, &Row::current_err);
    g_ppg12_top[s] = make_graph(rs, ("g_ppg12_ian_top_" + s).c_str(), &Row::ppg12, &Row::ppg12_err);
    g_current_ratio[s] = make_graph(rs, ("g_the94_ian_ratio_" + s).c_str(), &Row::current_over_fit, &Row::current_over_fit_err);
    g_ppg12_ratio[s] = make_graph(rs, ("g_ppg12_ian_ratio_" + s).c_str(), &Row::ppg12_over_fit, &Row::ppg12_over_fit_err);
    style_top(g_current_top[s], colors.at(s), 20);
    style_top(g_ppg12_top[s], colors.at(s), 24);
    style_ratio(g_current_ratio[s], 20);
    style_ratio(g_ppg12_ratio[s], 24);
  }

  auto* g_fit = make_fit_graph(rows);
  g_fit->SetLineColor(kRed - 4);
  g_fit->SetLineWidth(2);
  g_fit->SetMarkerSize(0);

  TCanvas* c4 = new TCanvas("can_the94_ppg12_jet_ian_overlay", "", 800, 889);
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

  for (const auto& s : samples) g_current_top[s]->Draw("PE same");
  for (const auto& s : samples) g_ppg12_top[s]->Draw("PE same");
  g_fit->Draw("L same");

  myText(0.47, 0.93, 1, strleg1.c_str(), 0.050);
  myText(0.47, 0.875, 1, strleg2.c_str(), 0.050);
  myText(0.47, 0.82, 1, strMC.c_str(), 0.050);

  myMarkerLineText(0.57, 0.75, 1, colors.at("jet8"), 20, colors.at("jet8"), 1, "jet8", 0.045, true);
  myMarkerLineText(0.57, 0.70, 1, colors.at("jet12"), 20, colors.at("jet12"), 1, "jet12", 0.045, true);
  myMarkerLineText(0.57, 0.65, 1, colors.at("jet20"), 20, colors.at("jet20"), 1, "jet20", 0.045, true);
  myMarkerLineText(0.57, 0.60, 1, colors.at("jet30"), 20, colors.at("jet30"), 1, "jet30", 0.045, true);
  myMarkerLineText(0.57, 0.55, 1, colors.at("jet40"), 20, colors.at("jet40"), 1, "jet40", 0.045, true);

  myText(0.19, 0.31, 1, "source", 0.036);
  myMarkerLineText(0.19, 0.26, 1, kBlack, 24, kBlack, 1, "PPG12 IAN", 0.038, true);
  myMarkerLineText(0.19, 0.21, 1, kBlack, 20, kBlack, 1, "THE-94 current", 0.038, true);
  myMarkerLineText(0.19, 0.16, 0, kRed - 4, 0, kRed - 4, 1, "PPG12 fit", 0.038, true);

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

  for (const auto& s : samples) g_current_ratio[s]->Draw("PE same");
  for (const auto& s : samples) g_ppg12_ratio[s]->Draw("PE same");
  lineone->Draw("L");

  c4->SaveAs(png_path.c_str());

  std::ofstream out(manifest_path);
  out << std::setprecision(12);
  out << "{\n";
  out << "  \"status\": \"ok_ian_source_overlay_the94_current\",\n";
  out << "  \"plot_png\": \"" << png_path << "\",\n";
  out << "  \"input_points_csv\": \"" << csv_path << "\",\n";
  out << "  \"renderer\": \"ROOT with ppg12codeGit/plotting/plotcommon.h\",\n";
  out << "  \"pixel_size\": [800, 889],\n";
  out << "  \"reference_scope\": \"PPG12 IAN/no-suffix source-scope overlay using exact IAN-style axes and MC/Fit panel\",\n";
  out << "  \"display_treatment\": \"THE-94 current points are displayed after per-sample median source-scope alignment; nominal current values remain in the companion CSV.\",\n";
  out << "  \"interpretation\": \"presentation/source-scope comparison, not nominal no-scale parity\",\n";
  out << "  \"displayed_bins\": " << rows.size() << "\n";
  out << "}\n";
  out.close();
}
