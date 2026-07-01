#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "TBox.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TGaxis.h"
#include "TGraph.h"
#include "TH1.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPad.h"
#include "TPaveText.h"
#include "TStyle.h"

namespace {

struct Row {
  double lo = 0;
  double hi = 0;
  double center = 0;
  double curB = 0;
  double ppgB = 0;
  double ratB = 0;
  double curC = 0;
  double ppgC = 0;
  double ratC = 0;
  double curD = 0;
  double ppgD = 0;
  double ratD = 0;
};

std::vector<std::string> split_csv_line(const std::string& line) {
  std::vector<std::string> out;
  std::stringstream ss(line);
  std::string item;
  while (std::getline(ss, item, ',')) out.push_back(item);
  return out;
}

std::vector<Row> read_rows(const std::string& path) {
  std::ifstream in(path);
  std::vector<Row> rows;
  if (!in) {
    std::cerr << "Cannot open " << path << std::endl;
    return rows;
  }

  std::string line;
  std::getline(in, line);  // header
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    auto vals = split_csv_line(line);
    if (vals.size() < 12) continue;
    Row r;
    r.lo = std::stod(vals[0]);
    r.hi = std::stod(vals[1]);
    r.center = std::stod(vals[2]);
    r.curB = std::stod(vals[3]);
    r.ppgB = std::stod(vals[4]);
    r.ratB = std::stod(vals[5]);
    r.curC = std::stod(vals[6]);
    r.ppgC = std::stod(vals[7]);
    r.ratC = std::stod(vals[8]);
    r.curD = std::stod(vals[9]);
    r.ppgD = std::stod(vals[10]);
    r.ratD = std::stod(vals[11]);
    rows.push_back(r);
  }
  return rows;
}

TGraph* make_graph(const std::vector<Row>& rows, double Row::*member) {
  auto* g = new TGraph(static_cast<int>(rows.size()));
  for (int i = 0; i < static_cast<int>(rows.size()); ++i) {
    g->SetPoint(i, rows[i].center, rows[i].*member);
  }
  return g;
}

void style_graph(TGraph* g, int color, int marker, bool open) {
  g->SetMarkerColor(color);
  g->SetLineColor(color);
  g->SetMarkerStyle(marker);
  g->SetMarkerSize(1.25);
  g->SetLineWidth(2);
  if (open) g->SetFillStyle(0);
}

void draw_text(double x, double y, const char* text, double size, int font = 42,
               int color = kBlack, int align = 11) {
  TLatex t;
  t.SetNDC();
  t.SetTextFont(font);
  t.SetTextColor(color);
  t.SetTextSize(size);
  t.SetTextAlign(align);
  t.DrawLatex(x, y, text);
}

void draw_note_box(double x1, double y1, double x2, double y2, int fill_color,
                   int line_color) {
  auto* b = new TPaveText(x1, y1, x2, y2, "NDC");
  b->SetFillColor(fill_color);
  b->SetLineColor(line_color);
  b->SetLineWidth(2);
  b->SetBorderSize(1);
  b->Draw();
}

void configure_frame(TH1* frame, const char* ytitle, bool show_x_title,
                     bool left_axis) {
  frame->GetXaxis()->SetLimits(9.5, 36.5);
  frame->GetXaxis()->SetTitle(show_x_title ? "reco cluster E_{T} [GeV]" : "");
  frame->GetXaxis()->SetTitleSize(0.070);
  frame->GetXaxis()->SetLabelSize(0.060);
  frame->GetXaxis()->SetTitleOffset(0.95);
  frame->GetXaxis()->SetNdivisions(506);
  frame->GetYaxis()->SetTitle(ytitle);
  frame->GetYaxis()->SetTitleSize(left_axis ? 0.068 : 0.0);
  frame->GetYaxis()->SetLabelSize(0.058);
  frame->GetYaxis()->SetTitleOffset(1.0);
  frame->GetYaxis()->SetNdivisions(505);
  frame->GetYaxis()->SetNoExponent(kTRUE);
  frame->GetYaxis()->SetDecimals(kTRUE);
}

void draw_component_panel(TPad* pad, const std::vector<Row>& rows,
                          const char* title, double Row::*cur,
                          double Row::*ppg, double ymax, int color,
                          int marker, bool show_y_title) {
  pad->cd();
  pad->SetFillColor(kWhite);
  pad->SetFrameLineWidth(2);
  pad->SetLeftMargin(show_y_title ? 0.18 : 0.10);
  pad->SetRightMargin(0.045);
  pad->SetTopMargin(0.15);
  pad->SetBottomMargin(0.04);
  pad->SetTicks(1, 1);
  auto* frame = pad->DrawFrame(9.5, 0.0, 36.5, ymax);
  configure_frame(frame, show_y_title ? "leakage / A" : "", false,
                  show_y_title);
  frame->SetTitle("");

  auto* g_ppg = make_graph(rows, ppg);
  auto* g_cur = make_graph(rows, cur);
  style_graph(g_ppg, kBlack, 20, false);
  style_graph(g_cur, color, marker, true);
  g_ppg->Draw("P SAME");
  g_cur->Draw("P SAME");

  draw_text(0.08, 0.90, title, 0.075, 62, kBlack, 11);
  if (std::string(title).find("B") != std::string::npos) {
    auto* leg = new TLegend(0.52, 0.70, 0.95, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.052);
    leg->AddEntry(g_ppg, "PPG12 reference", "p");
    leg->AddEntry(g_cur, "Current rerun", "p");
    leg->Draw();
  }
}

void draw_ratio_panel(TPad* pad, const std::vector<Row>& rows,
                      double Row::*ratio, int color, int marker,
                      bool show_y_title) {
  pad->cd();
  pad->SetFillColor(kWhite);
  pad->SetFrameLineWidth(2);
  pad->SetLeftMargin(show_y_title ? 0.18 : 0.10);
  pad->SetRightMargin(0.045);
  pad->SetTopMargin(0.04);
  pad->SetBottomMargin(0.25);
  pad->SetTicks(1, 1);
  auto* frame = pad->DrawFrame(9.5, 0.30, 36.5, 1.10);
  configure_frame(frame, show_y_title ? "current / PPG12" : "",
                  true, show_y_title);
  frame->SetTitle("");

  TLine one(9.5, 1.0, 36.5, 1.0);
  one.SetLineColor(kGray + 2);
  one.SetLineStyle(2);
  one.SetLineWidth(2);
  one.Draw();

  auto* g = make_graph(rows, ratio);
  style_graph(g, color, marker, true);
  g->Draw("P SAME");
}

}  // namespace

void build_where_stuck_slide() {
  const std::string base =
      "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/ppg12PhotonYield/"
      "ppg12_photon_yield_v1_data_20260620/purity_fig29_comparison/"
      "ppg12_ratio_diagnostic/leakage_components";
  const std::string csv_path =
      base + "/current_vs_ppg12_fig29_leakage_components_bcd.csv";
  const std::string out_dir = base + "/slide_candidates";
  const std::string out_png =
      out_dir + "/the76_where_stuck_ppg12_leakage_bcd_slide_20260629.png";

  auto rows = read_rows(csv_path);
  if (rows.empty()) {
    std::cerr << "No rows read from " << csv_path << std::endl;
    return;
  }

  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(0);
  gStyle->SetLineScalePS(1);
  TGaxis::SetMaxDigits(5);

  auto* c = new TCanvas("c", "THE-76 where stuck", 2560, 1440);
  c->SetCanvasSize(2560, 1440);
  c->SetFillColor(kWhite);
  c->cd();

  const int blue = TColor::GetColor("#1f77b4");
  const int green = TColor::GetColor("#2ca02c");
  const int orange = TColor::GetColor("#d95f02");
  const int lightBlue = TColor::GetColor("#EEF5FF");
  const int lightGold = TColor::GetColor("#FFF7E8");
  const int outlineBlue = TColor::GetColor("#6E8CB8");
  const int outlineGold = TColor::GetColor("#D49A42");

  draw_text(0.045, 0.948, "Where we are stuck: PPG12 leakage parity", 0.041,
            62);
  draw_text(0.955, 0.952, "THE-76 pp-SIM signal-only diagnostic", 0.020, 42,
            kGray + 2, 31);

  draw_note_box(0.045, 0.805, 0.505, 0.915, lightBlue, outlineBlue);
  draw_text(0.062, 0.885,
            "#bf{Remaining issue:} B leakage tracks PPG12, but C and D are "
            "still low.",
            0.027, 42, kBlack, 11);
  draw_text(0.062, 0.846,
            "Core bins: B/A = 0.88-0.95, C/A = 0.66-0.73, D/A = 0.71-0.85; "
            "high-ET C/D fall further.",
            0.023, 42, kBlack, 11);

  draw_note_box(0.525, 0.805, 0.955, 0.915, lightGold, outlineGold);
  draw_text(0.542, 0.885,
            "#bf{Likely discrepancy:} non-tight BDT sideband / candidate "
            "classification",
            0.025, 42, kBlack, 11);
  draw_text(0.542, 0.846,
            "before isolation. Next: same-cluster feature/score parity vs. "
            "PPG12.",
            0.023, 42, kBlack, 11);

  draw_text(0.045, 0.770,
            "Current overlay after global-MBD component-mix rerun: the deficit "
            "is the total non-tight signal population, not the B isolation "
            "region.",
            0.024, 42, kGray + 2, 11);

  const double left = 0.045;
  const double right = 0.955;
  const double gap = 0.024;
  const double colw = (right - left - 2 * gap) / 3.0;
  const double top_y1 = 0.405;
  const double top_y2 = 0.745;
  const double bot_y1 = 0.090;
  const double bot_y2 = 0.390;

  std::array<TPad*, 6> pads;
  for (int i = 0; i < 3; ++i) {
    double x1 = left + i * (colw + gap);
    double x2 = x1 + colw;
    pads[i] = new TPad(Form("top_%d", i), "", x1, top_y1, x2, top_y2);
    pads[i + 3] = new TPad(Form("bot_%d", i), "", x1, bot_y1, x2, bot_y2);
    pads[i]->Draw();
    pads[i + 3]->Draw();
  }

  draw_component_panel(pads[0], rows, "Region B", &Row::curB, &Row::ppgB,
                       0.165, blue, 24, true);
  draw_component_panel(pads[1], rows, "Region C", &Row::curC, &Row::ppgC,
                       0.90, green, 25, false);
  draw_component_panel(pads[2], rows, "Region D", &Row::curD, &Row::ppgD,
                       0.095, orange, 26, false);

  draw_ratio_panel(pads[3], rows, &Row::ratB, blue, 24, true);
  draw_ratio_panel(pads[4], rows, &Row::ratC, green, 25, false);
  draw_ratio_panel(pads[5], rows, &Row::ratD, orange, 26, false);

  c->cd();
  draw_text(0.045, 0.035,
            "Readout: if the same clusters receive too-tight BDT scores, A and "
            "B can look aligned while C+D remains suppressed.",
            0.023, 42, kGray + 2, 11);

  gSystem->mkdir(out_dir.c_str(), true);
  c->SaveAs(out_png.c_str());
  std::cout << out_png << std::endl;
}
