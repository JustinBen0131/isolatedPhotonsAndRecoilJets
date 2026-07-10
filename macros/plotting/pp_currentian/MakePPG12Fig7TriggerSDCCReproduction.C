R__ADD_INCLUDE_PATH(/Users/patsfan753/Desktop/ThesisAnalysis/ppg12codeGit/plotting)
#include "/Users/patsfan753/Desktop/ThesisAnalysis/ppg12codeGit/plotting/plotcommon.h"

#include <TCanvas.h>
#include <TEfficiency.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TStyle.h>
#include <TSystem.h>

#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace
{
struct Row
{
  int bin = 0;
  double lo = 0.0;
  double hi = 0.0;
  double x = 0.0;
  double width = 1.0;
  double den = 0.0;
  double n29 = 0.0;
  double n30 = 0.0;
  double n31 = 0.0;
};

struct PlateauRow
{
  int bin = 0;
  double lo = 0.0;
  double hi = 0.0;
  double x = 0.0;
  double width = 1.0;
  double eff = 0.0;
  double err = 0.0;
  double err_low = 0.0;
  double err_high = 0.0;
};

std::vector<std::string> split_csv(const std::string& line)
{
  std::vector<std::string> out;
  std::string field;
  bool quote = false;
  for (const char ch : line)
  {
    if (ch == '"')
    {
      quote = !quote;
      continue;
    }
    if (ch == ',' && !quote)
    {
      out.push_back(field);
      field.clear();
      continue;
    }
    field.push_back(ch);
  }
  out.push_back(field);
  for (auto& value : out)
  {
    while (!value.empty() && (value.back() == '\r' || value.back() == '\n' || value.back() == ' ' || value.back() == '\t'))
    {
      value.pop_back();
    }
    size_t first = 0;
    while (first < value.size() && (value[first] == ' ' || value[first] == '\t')) ++first;
    if (first > 0) value.erase(0, first);
  }
  return out;
}

double atod(const std::vector<std::string>& v, const std::map<std::string, int>& col, const std::string& key)
{
  const auto it = col.find(key);
  if (it == col.end() || it->second < 0 || it->second >= static_cast<int>(v.size())) return 0.0;
  return std::atof(v[it->second].c_str());
}

std::vector<Row> read_rows(const std::string& path)
{
  std::ifstream in(path);
  if (!in)
  {
    std::cerr << "ERROR: cannot open " << path << std::endl;
    return {};
  }

  std::string line;
  if (!std::getline(in, line)) return {};
  auto header = split_csv(line);
  std::map<std::string, int> col;
  for (int i = 0; i < static_cast<int>(header.size()); ++i) col[header[i]] = i;

  std::vector<Row> rows;
  while (std::getline(in, line))
  {
    if (line.empty()) continue;
    auto v = split_csv(line);
    Row r;
    r.bin = static_cast<int>(atod(v, col, "bin"));
    r.lo = atod(v, col, "xlow");
    r.hi = atod(v, col, "xhigh");
    r.x = atod(v, col, "xcenter");
    r.width = atod(v, col, "width");
    r.den = atod(v, col, "den_scaled10");
    r.n29 = atod(v, col, "num_live29");
    r.n30 = atod(v, col, "num_live30");
    r.n31 = atod(v, col, "num_live31");
    rows.push_back(r);
  }
  return rows;
}

std::vector<PlateauRow> read_plateau_rows(const std::string& path)
{
  std::ifstream in(path);
  if (!in)
  {
    std::cerr << "ERROR: cannot open " << path << std::endl;
    return {};
  }

  std::string line;
  if (!std::getline(in, line)) return {};
  auto header = split_csv(line);
  std::map<std::string, int> col;
  for (int i = 0; i < static_cast<int>(header.size()); ++i) col[header[i]] = i;

  std::vector<PlateauRow> rows;
  while (std::getline(in, line))
  {
    if (line.empty()) continue;
    auto v = split_csv(line);
    PlateauRow r;
    r.bin = static_cast<int>(atod(v, col, "bin"));
    r.lo = atod(v, col, "xlow");
    r.hi = atod(v, col, "xhigh");
    r.x = atod(v, col, "xcenter");
    r.width = atod(v, col, "width");
    r.eff = atod(v, col, "eff");
    const bool has_asymm = col.find("eylow") != col.end() && col.find("eyhigh") != col.end();
    r.err_low = has_asymm ? atod(v, col, "eylow") : atod(v, col, "eff_error");
    r.err_high = has_asymm ? atod(v, col, "eyhigh") : atod(v, col, "eff_error");
    r.err = 0.5 * (r.err_low + r.err_high);
    rows.push_back(r);
  }
  return rows;
}

std::vector<double> edges_from_rows(const std::vector<Row>& rows)
{
  std::vector<double> edges;
  if (rows.empty()) return edges;
  edges.push_back(rows.front().lo);
  for (const auto& r : rows) edges.push_back(r.hi);
  return edges;
}

TH1D* make_hist(const std::vector<Row>& rows, const std::vector<double>& edges, const char* name, double Row::*value)
{
  auto* h = new TH1D(name, "", static_cast<int>(edges.size()) - 1, edges.data());
  h->Sumw2();
  for (int i = 0; i < static_cast<int>(rows.size()); ++i)
  {
    const double y = rows[i].*value;
    h->SetBinContent(i + 1, y);
    h->SetBinError(i + 1, y > 0.0 ? std::sqrt(y) : 0.0);
  }
  return h;
}

void style_hist(TH1D* h, int color, int marker, int line_style, int line_width = 2)
{
  h->SetStats(false);
  h->SetLineColor(color);
  h->SetMarkerColor(color);
  h->SetLineStyle(line_style);
  h->SetLineWidth(line_width);
  h->SetMarkerStyle(marker);
  h->SetMarkerSize(0.72);
}

TGraphAsymmErrors* make_eff_graph(const std::vector<Row>& rows,
                                  const char* name,
                                  double Row::*pass_value,
                                  int color,
                                  int marker)
{
  auto* g = new TGraphAsymmErrors();
  g->SetName(name);
  for (int i = 0; i < static_cast<int>(rows.size()); ++i)
  {
    const double den = rows[i].den;
    const double pass = rows[i].*pass_value;
    if (den <= 0.0) continue;
    const double eff = pass / den;
    const double lo = TEfficiency::ClopperPearson(static_cast<int>(std::llround(den)),
                                                  static_cast<int>(std::llround(pass)),
                                                  0.682689492137,
                                                  false);
    const double hi = TEfficiency::ClopperPearson(static_cast<int>(std::llround(den)),
                                                  static_cast<int>(std::llround(pass)),
                                                  0.682689492137,
                                                  true);
    const int p = g->GetN();
    g->SetPoint(p, rows[i].x, eff);
    g->SetPointError(p, 0.0, 0.0, std::max(0.0, eff - lo), std::max(0.0, hi - eff));
  }
  g->SetLineColor(color);
  g->SetMarkerColor(color);
  g->SetMarkerStyle(marker);
  g->SetMarkerSize(0.86);
  g->SetLineWidth(1);
  return g;
}

TH1D* make_eff_fit_hist(const std::vector<Row>& rows,
                        const std::vector<double>& edges,
                        const char* name,
                        double Row::*pass_value)
{
  auto* h = new TH1D(name, "", static_cast<int>(edges.size()) - 1, edges.data());
  for (int i = 0; i < static_cast<int>(rows.size()); ++i)
  {
    const double den = rows[i].den;
    const double pass = rows[i].*pass_value;
    if (den <= 0.0) continue;
    const double eff = pass / den;
    const double lo = TEfficiency::ClopperPearson(static_cast<int>(std::llround(den)),
                                                  static_cast<int>(std::llround(pass)),
                                                  0.682689492137,
                                                  false);
    const double hi = TEfficiency::ClopperPearson(static_cast<int>(std::llround(den)),
                                                  static_cast<int>(std::llround(pass)),
                                                  0.682689492137,
                                                  true);
    const double err = std::max(1e-12, 0.5 * ((eff - lo) + (hi - eff)));
    h->SetBinContent(i + 1, eff);
    h->SetBinError(i + 1, err);
  }
  return h;
}

TGraphAsymmErrors* make_plateau_graph(const std::vector<PlateauRow>& rows)
{
  auto* g = new TGraphAsymmErrors();
  g->SetName("g_sdcc_fig7_bit30_plateau");
  for (const auto& r : rows)
  {
    const int p = g->GetN();
    g->SetPoint(p, r.x, r.eff);
    g->SetPointError(p, 0.0, 0.0, r.err_low, r.err_high);
  }
  g->SetLineColor(kRed + 1);
  g->SetMarkerColor(kRed + 1);
  g->SetMarkerStyle(20);
  g->SetMarkerSize(0.86);
  g->SetLineWidth(1);
  return g;
}

TH1D* make_plateau_fit_hist(const std::vector<PlateauRow>& rows)
{
  std::vector<double> edges;
  if (rows.empty()) return nullptr;
  edges.push_back(rows.front().lo);
  for (const auto& r : rows) edges.push_back(r.hi);
  auto* h = new TH1D("h_sdcc_fig7_bit30_plateau_fit", "", static_cast<int>(edges.size()) - 1, edges.data());
  for (int i = 0; i < static_cast<int>(rows.size()); ++i)
  {
    h->SetBinContent(i + 1, rows[i].eff);
    h->SetBinError(i + 1, rows[i].err > 0.0 ? rows[i].err : 1e-12);
  }
  return h;
}

double plateau_weighted_mean(const std::vector<PlateauRow>& rows, const double xmin, const double xmax)
{
  double num = 0.0;
  double den = 0.0;
  for (const auto& r : rows)
  {
    if (r.x < xmin || r.x > xmax || r.err <= 0.0) continue;
    const double w = 1.0 / (r.err * r.err);
    num += w * r.eff;
    den += w;
  }
  return den > 0.0 ? num / den : 0.996;
}

void draw_labels(double x, double y, double dy, double size)
{
  myText(x, y, 1, strleg1.c_str(), size, 0);
  myText(x, y - dy, 1, strleg2_1.c_str(), size, 0);
  myText(x, y - 2.0 * dy, 1, strleg3.c_str(), size, 0);
}
}  // namespace

void MakePPG12Fig7TriggerSDCCReproduction(const char* csv_override = "", const char* outdir_override = "")
{
  const std::string base =
      "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/ppg12PhotonYield/"
      "ppg12_photon_yield_v1_data_20260620/shower_shape_reference_validation/fig7_trigger";
  const std::string csv_path = std::string(csv_override).empty()
                                   ? base + "/ppg12_fig7_sdcc_trigger_arrays.csv"
                                   : std::string(csv_override);
  const std::string plateau_asymm_csv = base + "/ppg12_fig7_sdcc_bit30_plateau_trigger_root_asymm_from_numden.csv";
  const std::string plateau_csv = !gSystem->AccessPathName(plateau_asymm_csv.c_str())
                                      ? plateau_asymm_csv
                                      : base + "/ppg12_fig7_sdcc_bit30_plateau_trigger_root.csv";
  const std::string outdir = std::string(outdir_override).empty() ? base : std::string(outdir_override);
  gSystem->mkdir(outdir.c_str(), true);

  auto rows = read_rows(csv_path);
  if (rows.empty())
  {
    std::cerr << "ERROR: no rows from " << csv_path << std::endl;
    return;
  }
  const auto edges = edges_from_rows(rows);
  auto plateau_rows = read_plateau_rows(plateau_csv);

  init_plot();
  gStyle->SetOptStat(0);
  gStyle->SetErrorX(0.0);
  gStyle->SetEndErrorSize(3);

  auto* h_den = make_hist(rows, edges, "h_sdcc_fig7_scaled10", &Row::den);
  auto* h_29 = make_hist(rows, edges, "h_sdcc_fig7_live29", &Row::n29);
  auto* h_30 = make_hist(rows, edges, "h_sdcc_fig7_live30", &Row::n30);
  auto* h_31 = make_hist(rows, edges, "h_sdcc_fig7_live31", &Row::n31);

  auto* e29 = make_eff_graph(rows, "g_eff_bit29", &Row::n29, kAzure + 2, 25);
  auto* e30 = make_eff_graph(rows, "g_eff_bit30", &Row::n30, kRed + 1, 20);
  auto* e31 = make_eff_graph(rows, "g_eff_bit31", &Row::n31, kSpring - 6, 26);

  // Panel 1: trigger-count overlay, matching the log-y Fig. 7 left panel.
  {
    auto* c = new TCanvas("c_fig7_overlay", "", 600, 600);
    c->SetLogy();
    c->SetLeftMargin(0.17);
    c->SetRightMargin(0.04);
    c->SetTopMargin(0.06);
    c->SetBottomMargin(0.15);

    auto* d_den = static_cast<TH1D*>(h_den->Clone("h_sdcc_fig7_scaled10_density"));
    auto* d_29 = static_cast<TH1D*>(h_29->Clone("h_sdcc_fig7_live29_density"));
    auto* d_30 = static_cast<TH1D*>(h_30->Clone("h_sdcc_fig7_live30_density"));
    auto* d_31 = static_cast<TH1D*>(h_31->Clone("h_sdcc_fig7_live31_density"));
    d_den->Scale(1.0, "width");
    d_29->Scale(1.0, "width");
    d_30->Scale(1.0, "width");
    d_31->Scale(1.0, "width");
    style_hist(d_den, kBlack, 20, 1, 2);
    style_hist(d_29, kAzure + 2, 25, 2, 2);
    style_hist(d_30, kRed + 1, 24, 2, 3);
    style_hist(d_31, kSpring - 6, 26, 2, 2);

    d_den->GetXaxis()->SetRangeUser(5, 15);
    d_den->GetYaxis()->SetRangeUser(0.5, 5.0 * d_den->GetMaximum());
    d_den->GetXaxis()->SetTitle("#it{E}_{T}^{cluster} [GeV]");
    d_den->GetYaxis()->SetTitle("Clusters / GeV");
    d_den->GetYaxis()->SetTitleOffset(1.45);
    d_den->GetXaxis()->SetTitleOffset(1.12);
    d_den->Draw("E1");
    d_29->Draw("E1 same");
    d_30->Draw("E1 same");
    d_31->Draw("E1 same");

    auto* leg = new TLegend(0.42, 0.62, 0.95, 0.86);
    legStyle(leg, 0.20, 0.034);
    leg->AddEntry(d_den, "scaled[10]  (MBD N&S)", "lep");
    leg->AddEntry(d_29, "scaled[10] & live[29]", "lep");
    leg->AddEntry(d_30, "scaled[10] & live[30]", "lep");
    leg->AddEntry(d_31, "scaled[10] & live[31]", "lep");
    leg->Draw("same");

    const float xpos = 0.22, xpos2 = 0.93, ypos = 0.88;
    const float dy = 0.054, fontsize = 0.040;
    myText(xpos, ypos - 0 * dy, 1, strleg1.c_str(), fontsize, 0);
    myText(xpos, ypos - 1 * dy, 1, strleg2_1.c_str(), fontsize, 0);
    myText(xpos2, ypos - 0 * dy, 1, strleg3.c_str(), fontsize, 1);
    c->SaveAs((outdir + "/ppg12_fig7_cluster_et_overlay_sdcc_reproduction.png").c_str());
  }

  // Panel 2: turn-on efficiency curves for bits 29, 30, and 31.
  {
    auto* c = new TCanvas("c_fig7_turnon", "", 640, 620);
    c->SetLeftMargin(0.16);
    c->SetRightMargin(0.04);
    c->SetTopMargin(0.06);
    c->SetBottomMargin(0.15);

    auto* frame = new TH2F("frame_fig7_turnon", "", 100, 5, 15, 100, 0.0, 1.05);
    frame->SetStats(false);
    frame->GetXaxis()->SetTitle("#it{E}_{T}^{cluster} [GeV]");
    frame->GetYaxis()->SetTitle("#varepsilon_{L1}(Photon N GeV | MBD N&S)");
    frame->GetYaxis()->SetTitleOffset(1.30);
    frame->Draw("axis");

    struct Spec
    {
      TGraphAsymmErrors* e;
      const char* name;
      const char* label;
      int color;
      double Row::*pass_value;
    };
    std::vector<Spec> specs = {
        {e29, "bit29", "Photon 3 GeV (bit 29)", kAzure + 2, &Row::n29},
        {e30, "bit30", "Photon 4 GeV (bit 30, nom.)", kRed + 1, &Row::n30},
        {e31, "bit31", "Photon 5 GeV (bit 31)", kSpring - 6, &Row::n31},
    };

    auto* leg = new TLegend(0.48, 0.35, 0.93, 0.54);
    legStyle(leg, 0.20, 0.030);
    std::vector<TF1*> fits;
    for (const auto& s : specs)
    {
      s.e->Draw("p same");
      auto* hf = make_eff_fit_hist(rows, edges, Form("h_fit_%s", s.name), s.pass_value);
      auto* f = new TF1(Form("fgumbel_%s", s.name), "[0]*TMath::Exp(-TMath::Exp(-(x-[1])/[2]))", 5.0, 15.0);
      f->SetParameters(1.0, -2.0, 3.33);
      f->SetParLimits(0, 0.5, 1.02);
      f->SetParLimits(1, -50.0, 8.0);
      f->SetParLimits(2, 0.1, 50.0);
      f->SetLineColor(s.color);
      f->SetLineWidth(2);
      hf->Fit(f, "R Q N");
      f->Draw("same");
      fits.push_back(f);
      leg->AddEntry(s.e, s.label, "pl");
    }
    leg->Draw("same");
    draw_labels(0.21, 0.31, 0.050, 0.038);
    c->SaveAs((outdir + "/ppg12_fig7_trigger_turnon_sdcc_reproduction.png").c_str());
  }

  // Panel 3: bit-30 plateau zoom with the same Gumbel turn-on form used by
  // PPG12 for Fig. 7.  This panel must use the SDCC plateau source table:
  // the lower-stat trigger-array CSV is a separate diagnostic source and does
  // not reproduce the IAN plateau points.
  {
    auto* c = new TCanvas("c_fig7_plateau", "", 600, 600);
    c->SetLeftMargin(0.20);
    c->SetRightMargin(0.04);
    c->SetTopMargin(0.06);
    c->SetBottomMargin(0.15);

    auto* frame = new TH2F("frame_fig7_plateau", "", 100, 7, 30, 100, 0.985, 1.005);
    frame->SetStats(false);
    frame->GetXaxis()->SetTitle("#it{E}_{T}^{cluster} [GeV]");
    frame->GetYaxis()->SetTitle("#varepsilon_{4}(Photon 4 GeV | MBD N&S)");
    frame->GetYaxis()->SetTitleOffset(1.58);
    frame->GetYaxis()->SetTitleSize(0.047);
    frame->GetYaxis()->SetLabelSize(0.044);
    frame->GetXaxis()->SetTitleSize(0.047);
    frame->GetXaxis()->SetLabelSize(0.044);
    frame->GetYaxis()->SetNdivisions(505);
    frame->Draw("axis");

    auto* g_plateau = make_plateau_graph(plateau_rows.empty() ? std::vector<PlateauRow>{} : plateau_rows);
    g_plateau->SetMarkerSize(0.78);
    auto* h_plateau = make_plateau_fit_hist(plateau_rows);

    // Match the PPG12 turn-on fit form:
    //   eps(ET) = p0 * exp(-exp(-(ET - mu) / beta)).
    auto* h_fit = static_cast<TH1D*>(h_plateau->Clone("h_sdcc_fig7_bit30_plateau_gumbel_fit_input"));
    for (int i = 1; i <= h_fit->GetNbinsX(); ++i)
    {
      const double x = h_fit->GetBinCenter(i);
      const double err = h_fit->GetBinError(i);
      if (x < 7.0 || x > 20.0 || err <= 1e-10)
      {
        h_fit->SetBinContent(i, 0.0);
        h_fit->SetBinError(i, 0.0);
      }
    }

    auto* fgumbel = new TF1("fgumbel_sdcc_bit30_plateau",
                            "[0]*TMath::Exp(-TMath::Exp(-(x-[1])/[2]))",
                            7.0, 30.0);
    fgumbel->SetParameters(0.9966, 1.0, 1.7);
    fgumbel->SetParLimits(0, 0.990, 1.002);
    fgumbel->SetParLimits(1, -20.0, 8.0);
    fgumbel->SetParLimits(2, 0.1, 20.0);
    fgumbel->SetLineColor(kRed + 1);
    fgumbel->SetLineWidth(2);
    TFitResultPtr fr = h_fit->Fit(fgumbel, "R Q N S");

    TGraphErrors* gci = nullptr;
    if (fr.Get() != nullptr && fr->IsValid())
    {
      const int n_ci = 200;
      std::vector<double> x_arr(n_ci), ci_arr(n_ci);
      for (int i = 0; i < n_ci; ++i) x_arr[i] = 7.0 + i * 23.0 / (n_ci - 1);
      fr->GetConfidenceIntervals(n_ci, 1, 1, x_arr.data(), ci_arr.data(), 0.683, false);
      gci = new TGraphErrors(n_ci);
      for (int i = 0; i < n_ci; ++i)
      {
        gci->SetPoint(i, x_arr[i], fgumbel->Eval(x_arr[i]));
        gci->SetPointError(i, 0.0, ci_arr[i]);
      }
      gci->SetFillColorAlpha(kRed + 1, 0.25);
      gci->SetFillStyle(1001);
      gci->SetLineWidth(0);
      gci->Draw("3 same");
    }
    fgumbel->Draw("same");
    g_plateau->Draw("p same");

    auto* leg = new TLegend(0.42, 0.78, 0.98, 0.91);
    legStyle(leg, 0.20, 0.031);
    leg->AddEntry(g_plateau, "Photon 4 GeV (bit 30, nom.)", "pl");
    leg->AddEntry(gci ? static_cast<TObject*>(gci) : static_cast<TObject*>(fgumbel),
                  Form("Gumbel fit, plateau = %.4f", fgumbel->GetParameter(0)),
                  gci ? "f" : "l");
    leg->Draw("same");
    draw_labels(0.23, 0.29, 0.048, 0.037);
    c->SaveAs((outdir + "/ppg12_fig7_bit30_plateau_sdcc_reproduction.png").c_str());
    c->SaveAs((outdir + "/ppg12_fig7_bit30_plateau_gumbel_sdcc_reproduction.png").c_str());

    std::cout << "SDCC bit30 plateau Gumbel fit:"
              << " plateau=" << fgumbel->GetParameter(0)
              << " +/- " << fgumbel->GetParError(0)
              << " mu=" << fgumbel->GetParameter(1)
              << " beta=" << fgumbel->GetParameter(2)
              << " chi2/ndf=" << fgumbel->GetChisquare()
              << "/" << fgumbel->GetNDF() << std::endl;
  }
}
