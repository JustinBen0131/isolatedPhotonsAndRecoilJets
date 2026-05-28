#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>

namespace
{
constexpr double kSigmaPhoton12_pb = 2598.12425;   // 12 <= pT_filter^gamma < 20
constexpr double kSigmaPhoton20_pb = 133.317866;   // pT_filter^gamma >= 20

const char* kFile12 =
    "InputFiles/simEmbedded/"
    "RecoilJets_embeddedPhoton12_ALL_preselectionReference_tightReference_nonTightReference_baseVariant.root";
const char* kFile20 =
    "InputFiles/simEmbedded/"
    "RecoilJets_embeddedPhoton20_ALL_preselectionReference_tightReference_nonTightReference_baseVariant.root";
const char* kHistName = "SIM/h_embedStitch_filterPhotonPt_kept";
const char* kCountName = "SIM/cnt_SIM";
const char* kOutDir = "dataOutput/stitchDiagnostics/photon12_20_stitch_20260521";

std::unique_ptr<TH1> GetHist(TFile* f, const char* name, const char* cloneName)
{
  TH1* h = dynamic_cast<TH1*>(f ? f->Get(name) : nullptr);
  if (!h)
  {
    std::cerr << "[ERROR] Missing histogram " << name << std::endl;
    return nullptr;
  }
  std::unique_ptr<TH1> out(dynamic_cast<TH1*>(h->Clone(cloneName)));
  out->SetDirectory(nullptr);
  out->Sumw2(false);
  out->Sumw2();
  return out;
}

double GetCount(TFile* f)
{
  TH1* h = dynamic_cast<TH1*>(f ? f->Get(kCountName) : nullptr);
  if (!h) return 0.0;
  return h->GetBinContent(1);
}

double MeanRatio(const TH1* ratio, double lo, double hi)
{
  double sum = 0.0;
  int n = 0;
  for (int ib = 1; ib <= ratio->GetNbinsX(); ++ib)
  {
    const double x = ratio->GetBinCenter(ib);
    const double y = ratio->GetBinContent(ib);
    if (x < lo || x >= hi) continue;
    if (!std::isfinite(y) || y <= 0.0) continue;
    sum += y;
    ++n;
  }
  return n > 0 ? sum / n : std::numeric_limits<double>::quiet_NaN();
}

void StyleHist(TH1* h, int color, int marker, double size = 0.85)
{
  h->SetLineColor(color);
  h->SetMarkerColor(color);
  h->SetMarkerStyle(marker);
  h->SetMarkerSize(size);
  h->SetLineWidth(1);
  h->SetStats(false);
}
}

void MakeEmbeddedPhoton12Plus20Stitch35()
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetErrorX(0.5);

  std::unique_ptr<TFile> f12(TFile::Open(kFile12, "READ"));
  std::unique_ptr<TFile> f20(TFile::Open(kFile20, "READ"));
  if (!f12 || f12->IsZombie() || !f20 || f20->IsZombie())
  {
    std::cerr << "[ERROR] Cannot open input ROOT files." << std::endl;
    return;
  }

  const double n12 = GetCount(f12.get());
  const double n20 = GetCount(f20.get());
  if (n12 <= 0.0 || n20 <= 0.0)
  {
    std::cerr << "[ERROR] Bad event counts: n12=" << n12 << " n20=" << n20 << std::endl;
    return;
  }

  std::unique_ptr<TH1> h12 = GetHist(f12.get(), kHistName, "h_photon12_weighted");
  std::unique_ptr<TH1> h20 = GetHist(f20.get(), kHistName, "h_photon20_weighted");
  if (!h12 || !h20) return;

  h12->Rebin(2);
  h20->Rebin(2);
  h12->Scale(kSigmaPhoton12_pb / n12, "width");
  h20->Scale(kSigmaPhoton20_pb / n20, "width");

  std::unique_ptr<TH1> hSum(dynamic_cast<TH1*>(h12->Clone("h_photon12_20_sum")));
  hSum->SetDirectory(nullptr);
  hSum->Add(h20.get());
  h12->SetTitle("");
  h20->SetTitle("");
  hSum->SetTitle("");

  StyleHist(h12.get(), kBlue + 1, 20);
  StyleHist(h20.get(), kOrange + 7, 21);
  StyleHist(hSum.get(), kBlack, 24, 0.95);

  const double xMin = 10.0;
  const double xMax = 35.0;
  TF1 fit("fit_modified_power_law", "[0]*TMath::Power(1.0 + x/[1], -[2])", 12.0, xMax);
  fit.SetParameters(std::max(1.0, hSum->GetMaximum()) * 100.0, 3.0, 6.0);
  fit.SetParLimits(0, 0.0, 1.0e12);
  fit.SetParLimits(1, 0.1, 200.0);
  fit.SetParLimits(2, 0.1, 50.0);
  fit.SetLineColor(kBlack);
  fit.SetLineStyle(2);
  fit.SetLineWidth(2);
  TFitResultPtr fitResult = hSum->Fit(&fit, "QRS0");

  std::unique_ptr<TH1> hRatio(dynamic_cast<TH1*>(hSum->Clone("h_photon12_20_ratio")));
  hRatio->SetDirectory(nullptr);
  hRatio->Reset("ICES");
  for (int ib = 1; ib <= hSum->GetNbinsX(); ++ib)
  {
    const double x = hSum->GetBinCenter(ib);
    const double y = hSum->GetBinContent(ib);
    const double ey = hSum->GetBinError(ib);
    const double fy = fit.Eval(x);
    if (x < xMin || x > xMax || y <= 0.0 || fy <= 0.0) continue;
    hRatio->SetBinContent(ib, y / fy);
    hRatio->SetBinError(ib, ey / fy);
  }

  gSystem->mkdir(kOutDir, true);

  TCanvas c("c_photon12_20_stitch35", "Embedded photon 12+20 stitching", 1450, 960);
  TPad p1("p1", "", 0.0, 0.30, 1.0, 1.0);
  TPad p2("p2", "", 0.0, 0.0, 1.0, 0.30);
  p1.SetLeftMargin(0.12);
  p1.SetRightMargin(0.04);
  p1.SetTopMargin(0.08);
  p1.SetBottomMargin(0.02);
  p1.SetLogy(true);
  p2.SetLeftMargin(0.12);
  p2.SetRightMargin(0.04);
  p2.SetTopMargin(0.02);
  p2.SetBottomMargin(0.30);
  p1.Draw();
  p2.Draw();

  p1.cd();
  hSum->GetXaxis()->SetRangeUser(xMin, xMax);
  hSum->GetXaxis()->SetTitle("");
  hSum->GetXaxis()->SetTitleSize(0.0);
  hSum->GetXaxis()->SetLabelSize(0.0);
  h12->GetXaxis()->SetRangeUser(xMin, xMax);
  h20->GetXaxis()->SetRangeUser(xMin, xMax);
  hSum->GetYaxis()->SetTitle("d#sigma/dp_{T}^{#gamma, filter} [pb / GeV]");
  hSum->GetYaxis()->SetTitleSize(0.052);
  hSum->GetYaxis()->SetLabelSize(0.043);
  hSum->GetYaxis()->SetTitleOffset(1.10);
  hSum->SetMinimum(std::max(1.0e-3, hSum->GetMinimum(1.0e-30) * 0.4));
  hSum->SetMaximum(hSum->GetMaximum() * 3.8);
  hSum->Draw("E1");
  h12->Draw("E1 SAME");
  h20->Draw("E1 SAME");
  fit.Draw("SAME");
  hSum->Draw("E1 SAME");

  TLine boundary(20.0, hSum->GetMinimum(), 20.0, hSum->GetMaximum());
  boundary.SetLineColor(kGray + 1);
  boundary.SetLineStyle(3);
  boundary.Draw("SAME");

  TLatex text;
  text.SetNDC(true);
  text.SetTextFont(42);
  text.SetTextSize(0.038);
  text.DrawLatex(0.15, 0.30, "#it{#bf{sPHENIX}} Internal");
  text.DrawLatex(0.15, 0.24, "PYTHIA8 embedded photon+jet samples");
  text.DrawLatex(0.15, 0.18, "Ownership: Photon12 12-20, Photon20 #geq20 GeV");

  TLegend leg(0.62, 0.67, 0.93, 0.90);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextFont(42);
  leg.SetTextSize(0.036);
  leg.AddEntry(&fit, "Fit: A(1+p_{T}/p_{0})^{-n}", "l");
  leg.AddEntry(h12.get(), "Photon12: 12 #leq p_{T} < 20", "ep");
  leg.AddEntry(h20.get(), "Photon20: p_{T} #geq 20", "ep");
  leg.AddEntry(hSum.get(), "Weighted stitched sum", "ep");
  leg.Draw();

  p2.cd();
  hRatio->SetTitle("");
  hRatio->GetXaxis()->SetTitle("Leading truth photon p_{T} [GeV]");
  hRatio->GetYaxis()->SetTitle("stitched / fit");
  hRatio->GetXaxis()->SetRangeUser(xMin, xMax);
  hRatio->GetYaxis()->SetRangeUser(0.82, 1.18);
  hRatio->GetXaxis()->SetTitleSize(0.115);
  hRatio->GetXaxis()->SetLabelSize(0.090);
  hRatio->GetYaxis()->SetTitleSize(0.095);
  hRatio->GetYaxis()->SetLabelSize(0.082);
  hRatio->GetYaxis()->SetTitleOffset(0.55);
  hRatio->GetYaxis()->SetNdivisions(505);
  StyleHist(hRatio.get(), kBlack, 20, 0.85);
  hRatio->Draw("E1");
  TLine one(xMin, 1.0, xMax, 1.0);
  one.SetLineColor(kGray + 2);
  one.SetLineStyle(2);
  one.Draw("SAME");
  TLine b2(20.0, 0.82, 20.0, 1.18);
  b2.SetLineColor(kGray + 1);
  b2.SetLineStyle(3);
  b2.Draw("SAME");

  const std::string png = std::string(kOutDir) + "/embeddedPhoton12plus20_stitch_x35.png";
  c.SaveAs(png.c_str());

  const double below = MeanRatio(hRatio.get(), 18.0, 20.0);
  const double above = MeanRatio(hRatio.get(), 20.0, 22.0);
  const double jump = above / below;

  std::ofstream csv(std::string(kOutDir) + "/embeddedPhoton12plus20_stitch_x35.csv");
  csv << "pt_center,pt_low,pt_high,photon12_pb_per_GeV,photon12_err_pb_per_GeV,"
         "photon20_pb_per_GeV,photon20_err_pb_per_GeV,stitched_pb_per_GeV,"
         "stitched_err_pb_per_GeV,fit_pb_per_GeV,ratio,ratio_err\n";
  for (int ib = 1; ib <= hSum->GetNbinsX(); ++ib)
  {
    const double x = hSum->GetBinCenter(ib);
    if (x < xMin || x > xMax) continue;
    const double fy = fit.Eval(x);
    csv << std::setprecision(12)
        << x << "," << hSum->GetBinLowEdge(ib) << "," << hSum->GetBinLowEdge(ib + 1) << ","
        << h12->GetBinContent(ib) << "," << h12->GetBinError(ib) << ","
        << h20->GetBinContent(ib) << "," << h20->GetBinError(ib) << ","
        << hSum->GetBinContent(ib) << "," << hSum->GetBinError(ib) << ","
        << fy << "," << hRatio->GetBinContent(ib) << "," << hRatio->GetBinError(ib) << "\n";
  }

  std::ofstream json(std::string(kOutDir) + "/embeddedPhoton12plus20_stitch_x35_metrics.json");
  json << std::setprecision(12)
       << "{\n"
       << "  \"input_file12\": \"" << kFile12 << "\",\n"
       << "  \"input_file20\": \"" << kFile20 << "\",\n"
       << "  \"hist_key\": \"" << kHistName << "\",\n"
       << "  \"n12\": " << n12 << ",\n"
       << "  \"n20\": " << n20 << ",\n"
       << "  \"sigma_photon12_pb\": " << kSigmaPhoton12_pb << ",\n"
       << "  \"sigma_photon20_pb\": " << kSigmaPhoton20_pb << ",\n"
       << "  \"mean_ratio_18_20\": " << below << ",\n"
       << "  \"mean_ratio_20_22\": " << above << ",\n"
       << "  \"jump_20_over_18_20\": " << jump << ",\n"
       << "  \"fit_status\": " << fitResult->Status() << ",\n"
       << "  \"fit_chi2_ndf\": " << (fit.GetNDF() > 0 ? fit.GetChisquare() / fit.GetNDF() : -1.0) << "\n"
       << "}\n";

  std::cout << "[DONE] " << png << std::endl;
  std::cout << "[METRIC] mean ratio 18-20=" << below
            << " 20-22=" << above
            << " jump=" << jump << std::endl;
}
