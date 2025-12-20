// gaugeFeasibilityGrid_grouped.C
// UE false-fail gauge using PPG04 (Randomized ηφ, Iterative) Gaussian tails.
// (1) Compute P(E_T^iso > E_cut) per 9 centrality slices from (mu, sigma).
// (2) Compress to ATLAS-style bins: 0–10, 10–30, 30–50, 50–80 (PPG04 provides up to 70–80).
// (3) Print a tidy table and draw a 2×2 summary. Optionally also draw the original 3×3.
// (4) NEW: Write a publication-grade single-panel PNG for the 0–10% bin only.
//
// Inputs (from PPG04 Table A.5, "Randomized ηφ, Iterative"):
//   mu  = {0.31, 0.19, 0.09, 0.05, 0.016, -0.0022, -0.000055, 0.0019, 0.0067} GeV
//   sig = {4.80, 4.40, 3.80, 3.10, 2.60,  2.10,     1.70,      1.40,   1.20 } GeV
//
// NOTE: 50–90% is not available in PPG04; we use 50–80% as the widest peripheral bin.

#include "TCanvas.h"
#include "TF1.h"
#include "TMath.h"
#include "TLine.h"
#include "TBox.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TColor.h"
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <algorithm> // std::max

// ----------------------------- tunables -----------------------------
static constexpr double kCut   = 4.0;    // GeV, E_T^iso threshold (vertical line in plots)
static constexpr bool   kDraw3x3 = false; // also draw the 3×3 grid (set true if you want it)
static constexpr double kAlpha =  0.28;  // shading alpha
static constexpr double kXmin  = -10.0;  // GeV
static constexpr double kXmax  =  30.0;  // GeV

// ----------------------------- utils -------------------------------
inline double P_right(double cut, double mu, double sig) {
  // P(X > cut) for X~N(mu, sig): 0.5 * erfc( (cut - mu)/(sqrt(2)*sig) )
  return (sig>0.0) ? 0.5 * TMath::Erfc( (cut - mu)/(TMath::Sqrt2()*sig) ) : 0.0;
}

// Normalized Gaussian pdf
double pdf(double *x, double *p) {
  const double mu=p[0], s=p[1];
  if (s<=0) return 0.0;
  const double z=(x[0]-mu)/s;
  const double norm = 1.0/( TMath::Sqrt(2.0*TMath::Pi()) * s );
  return norm * TMath::Exp(-0.5*z*z);
}

// Shade area under f from xa to xb
void Shade(const TF1* f, double xa, double xb, int col, float a=0.30, int N=300){
  const double dx=(xb-xa)/N;
  for(int i=0;i<N;++i){
    const double x1=xa+i*dx, x2=x1+dx;
    const double y = std::max(f->Eval(x1), f->Eval(x2));
    auto b = new TBox(x1,0,x2,y);
    b->SetFillColorAlpha(col,a);
    b->SetLineColorAlpha(col,0);
    b->Draw("f");
  }
}

// Weighted average helper on probabilities
double weighted_mean(const std::vector<double>& val, const std::vector<double>& w) {
  double num = 0.0, den = 0.0;
  for (size_t i=0;i<val.size();++i) { num += w[i]*val[i]; den += w[i]; }
  return (den>0.0) ? num/den : 0.0;
}

// ----------------------------- publication-grade single panel -----------------------------
void DrawPublicationPanel_0to10(double muW, double sigW, double Pgt, double Plt, double ratioRL) {
  // Global style for publication-quality
  gStyle->SetOptStat(0);
  gStyle->SetTitleFont(42, "");  // Helvetica/Arial
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetTitleSize(0.05, "XYZ");
  gStyle->SetLabelSize(0.044, "XYZ");
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  TCanvas cpub("cpub","UE isolation (0–10%)", 900, 700);
  cpub.SetMargin(0.13, 0.04, 0.13, 0.08); // L, R, B, T

  // Function and axes
  auto f = new TF1("f_pub", pdf, kXmin, kXmax, 2);
  f->SetParameters(muW, sigW);
  f->SetLineWidth(4);
  f->SetLineColor(kBlack);
  f->SetTitle("0-10% Centrality;UE cone E_{T}^{iso}  [GeV];Probability density");
  f->Draw();

  // Cut line
  auto L = new TLine(kCut, 0, kCut, 1.25*f->Eval(kCut));
  L->SetLineColor(kBlack);
  L->SetLineWidth(4);
  L->SetLineStyle(7);
  L->Draw();

  // Shading
  Shade(f, kXmin, kCut, kAzure+1, 0.30);    // pass region (E<cut)
  Shade(f, kCut,  kXmax, kOrange+7, 0.30);  // fail region (E>cut)

  // Legend-like annotation box
  TPaveText pt(0.5, 0.58, 0.92, 0.88, "NDC");
  pt.SetFillColor(0);
  pt.SetFillStyle(0);
  pt.SetLineColor(0);
  pt.SetTextAlign(12);
  pt.SetTextFont(42);
  pt.SetTextSize(0.037);
  pt.AddText("#bf{PPG04 Randomized #eta#varphi, Iterative}");
  pt.AddText(Form("#mu = %.3g GeV,  #sigma = %.3g GeV", muW, sigW));
  pt.AddText(Form("E_{cut} = %.1f GeV", kCut));
  pt.AddText(Form("P(E < cut) = %.4f", Plt));
  pt.AddText(Form("P(E > cut) = %.4f", Pgt));
  pt.AddText(Form("Right/Left  = %.3f", ratioRL));
  pt.Draw();

  // Axis cosmetics
  f->GetXaxis()->SetNdivisions(510);
  f->GetYaxis()->SetNdivisions(506);
  f->GetXaxis()->SetTitleOffset(1.00);
  f->GetYaxis()->SetTitleOffset(1.10);

  cpub.Update();
  cpub.SaveAs("rc_0to10_pub.png"); // NEW single-panel, publication-grade PNG
}

// ----------------------------- main -------------------------------
void gaugeFeasibility() {
  gStyle->SetOptStat(0);

  // 9 centrality slices (widths in percentage points)
  const std::vector<std::string> lab9 = {
    "0-5%","5-10%","10-20%","20-30%","30-40%","40-50%","50-60%","60-70%","70-80%"
  };
  const std::vector<double> w9 = {5,5,10,10,10,10,10,10,10}; // weights for grouping

  // μ, σ from PPG04 (Randomized ηφ, Iterative)
  const std::vector<double> mu  = { 0.31,  0.19,  0.09,  0.05,  0.016, -0.0022, -5.5e-5, 0.0019, 0.0067 };
  const std::vector<double> sig = { 4.80,  4.40,  3.80,  3.10,  2.60,   2.10,     1.70,   1.40,   1.20  };

  // Precompute P_right, P_left, ratio for each of the 9 slices
  std::vector<double> PR(9), PL(9), RR(9);
  for (size_t i=0;i<9;++i) {
    PR[i] = P_right(kCut, mu[i], sig[i]);
    PL[i] = 1.0 - PR[i];
    RR[i] = (PL[i] > 0.0) ? PR[i]/PL[i] : 0.0;
  }

  // Grouping: 0–10 (0–5, 5–10), 10–30 (10–20, 20–30), 30–50 (30–40, 40–50), 50–80 (50–60, 60–70, 70–80)
  struct Group { std::string name; std::vector<int> idx; };
  const std::vector<Group> groups = {
    {"0-10%", {0,1}},
    {"10-30%",{2,3}},
    {"30-50%",{4,5}},
    {"50-80%",{6,7,8}}
  };

  // Aggregate probabilities (weighting by centrality width)
  struct OutRow { std::string name; double PR, PL, R; double muW, sigW; };
  std::vector<OutRow> rows; rows.reserve(groups.size());

  for (const auto& g : groups) {
    std::vector<double> pr, pl, ww, mu_g, sig_g;
    pr.reserve(g.idx.size()); pl.reserve(g.idx.size()); ww.reserve(g.idx.size());
    mu_g.reserve(g.idx.size()); sig_g.reserve(g.idx.size());
    for (int k : g.idx) {
      pr.push_back(PR[k]); pl.push_back(PL[k]); ww.push_back(w9[k]);
      mu_g.push_back(mu[k]); sig_g.push_back(sig[k]);
    }
    const double PRg = weighted_mean(pr, ww);
    const double PLg = weighted_mean(pl, ww);
    const double Rg  = (PLg>0.0) ? PRg/PLg : 0.0;

    // Weighted means of mu and sigma for DRAWING ONLY (probabilities already aggregated above)
    auto wavg = [&](const std::vector<double>& v){
      double num=0, den=0; for (size_t i=0;i<v.size();++i){ num += ww[i]*v[i]; den += ww[i]; }
      return (den>0)? num/den : 0.0;
    };
    const double muW  = wavg(mu_g);
    const double sigW = wavg(sig_g);

    rows.push_back({g.name, PRg, PLg, Rg, muW, sigW});
  }

  // --------------------- print a tidy table ---------------------
  std::cout.setf(std::ios::fixed);
  std::cout << "\nUE-only isolation false-fail summary (Gaussian tails, PPG04 Randomized ηφ, Iterative)\n";
  std::cout << "E_cut = " << kCut << " GeV\n";
  std::cout << "--------------------------------------------------------------------------\n";
  std::cout << "Centrality   P(E>cut)     P(E<cut)     Right/Left (=P> / P<)\n";
  std::cout << "--------------------------------------------------------------------------\n";
  std::cout << std::setprecision(4);
  for (const auto& r : rows) {
    std::cout << std::setw(9) << r.name << "   "
              << std::setw(9) << r.PR   << "   "
              << std::setw(9) << r.PL   << "   "
              << std::setw(9) << r.R    << "\n";
  }
  std::cout << "--------------------------------------------------------------------------\n";

  // --------------------- draw 2×2 summary panels (grouped) ---------------------
  TCanvas c("c","UE false-fail (grouped centralities)", 1100, 900);
  c.Divide(2,2,0.004,0.004);

  for (int ig=0; ig<4; ++ig) {
    c.cd(ig+1);
    const auto& entry = rows[ig];

    // DRAWING function uses weighted mu/sigma (for visualization only)
    auto f = new TF1(Form("f_grp_%d",ig), pdf, kXmin, kXmax, 2);
    f->SetParameters(entry.muW, entry.sigW);
    f->SetLineWidth(3);
    f->SetTitle(Form("%s;UE cone E_{T}^{iso}  [GeV];Probability density", entry.name.c_str()));
    f->Draw();

    auto L = new TLine(kCut, 0, kCut, 1.2*f->Eval(kCut));
    L->SetLineColor(kBlack); L->SetLineWidth(3); L->SetLineStyle(7); L->Draw();

    Shade(f, kXmin, kCut, kAzure+1, kAlpha);    // pass region (E<cut)
    Shade(f, kCut,  kXmax, kGreen+2, kAlpha);   // fail region (E>cut)

    TLatex tx; tx.SetNDC(true); tx.SetTextAlign(13); tx.SetTextFont(42); tx.SetTextSize(0.050);
    tx.DrawLatex(0.14, 0.86, Form("P(E>cut)  = %.4f", entry.PR));
    tx.DrawLatex(0.14, 0.79, Form("P(E<cut)  = %.4f", entry.PL));
    tx.DrawLatex(0.14, 0.72, Form("Right/Left = %.3f", entry.R));
    tx.DrawLatex(0.14, 0.64, Form("#mu_{draw}=%.3g, #sigma_{draw}=%.3g", entry.muW, entry.sigW));
  }

  c.Update();
  c.SaveAs("rc_grouped_2x2.png");

  // --------------------- NEW: publication-grade panel for 0–10% only ---------------------
  // Grab the aggregated 0–10% entry (rows[0])
  {
    const auto& r010 = rows[0];
    DrawPublicationPanel_0to10(r010.muW, r010.sigW, r010.PR, r010.PL, r010.R);
  }

  // --------------------- optionally also draw the original 3×3 ---------------------
  if (kDraw3x3) {
    TCanvas c9("c9","UE tail gauge: 9 centrality slices (PPG04 Rand. #eta#varphi, Iterative)", 1200, 1000);
    c9.Divide(3,3,0.002,0.002);

    for (int k=0;k<9;++k){
      c9.cd(k+1);
      auto f = new TF1(Form("f_%d",k), pdf, kXmin, kXmax, 2);
      f->SetParameters(mu[k], sig[k]);
      f->SetLineWidth(3);
      f->SetTitle(Form("%s;UE cone E_{T}^{iso}  [GeV];Probability density", lab9[k].c_str()));
      f->Draw();

      auto L = new TLine(kCut, 0, kCut, 1.2*f->Eval(kCut));
      L->SetLineColor(kBlack); L->SetLineWidth(3); L->SetLineStyle(7); L->Draw();

      Shade(f, kXmin, kCut, kAzure+1, kAlpha);
      Shade(f, kCut,  kXmax, kGreen+2, kAlpha);

      TLatex tx; tx.SetNDC(true); tx.SetTextAlign(33); tx.SetTextFont(42); tx.SetTextSize(0.048);
      tx.DrawLatex(0.86, 0.88, Form("#mu=%.3g, #sigma=%.3g", mu[k], sig[k]));
      tx.DrawLatex(0.86, 0.80, Form("P(E>cut)=%.4g", PR[k]));
      tx.DrawLatex(0.86, 0.73, Form("P(E<cut)=%.4g", PL[k]));
      tx.DrawLatex(0.86, 0.66, Form("Right/Left=%.3g", RR[k]));
    }
    c9.Update();
    c9.SaveAs("rc_grid_3x3.png");
  }

  // Done
  std::cout << "Saved rc_grouped_2x2.png and rc_0to10_pub.png";
  if (kDraw3x3) std::cout << " and rc_grid_3x3.png";
  std::cout << std::endl;
}
