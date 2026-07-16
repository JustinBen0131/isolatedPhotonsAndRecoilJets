// Render the current merged inclusive-MC counterpart of PPG12 IAN Fig. 25.
//
// This intentionally mirrors ppg12codeGit/plotting/plot_showershapes.C:
//   h2d_e11_to_e33_eta0_pt1_cut1 -> RebinX(4) -> ProfileX().
// The caller supplies the ROOT input resolved from the current-artifact pointer
// so the output cannot silently use an older campaign artifact.

#include <TCanvas.h>
#include <TError.h>
#include <TFile.h>
#include <TH2.h>
#include <TLatex.h>
#include <TProfile.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

#include <fstream>
#include <iomanip>

namespace
{
void writeManifest(const char *manifestPath,
                   const char *inputRoot,
                   const char *outputPng,
                   const char *currentPointer,
                   const char *rootSha256,
                   const char *histPath,
                   const int rebinX,
                   const double correlation)
{
  std::ofstream out(manifestPath);
  out << "{\n"
      << "  \"figure\": \"PPG12 IAN Fig. 25 current counterpart\",\n"
      << "  \"source_reference\": \"ppg12codeGit/plotting/plot_showershapes.C\",\n"
      << "  \"artifact_pointer\": \"" << currentPointer << "\",\n"
      << "  \"resolved_input_root\": \"" << inputRoot << "\",\n"
      << "  \"resolved_input_sha256\": \"" << rootSha256 << "\",\n"
      << "  \"histogram\": \"" << histPath << "\",\n"
      << "  \"selection\": \"inclusive background MC; |eta|<0.7; 14<pT<18 GeV; cut1 (NPB/preselection)\",\n"
      << "  \"mixture\": \"current canonical inclusive SIM: jet8-40, SI+DI, 0mrad+1p5mrad, additive period/component merge\",\n"
      << "  \"transformation\": \"TH2 RebinX(" << rebinX
      << "), x in [0,1], then ProfileX over all isolation bins\",\n"
      << "  \"correlation_factor\": " << std::fixed << std::setprecision(6) << correlation << ",\n"
      << "  \"output_png\": \"" << outputPng << "\"\n"
      << "}\n";
}
}  // namespace

void MakePPG12Fig25Current(const char *inputRoot,
                           const char *outputPng,
                           const char *manifestPath,
                           const char *currentPointer,
                           const char *rootSha256,
                           const char *histPath = "SIM/h2d_e11_to_e33_eta0_pt1_cut1",
                           const char *xAxisTitle = "e11_to_e33",
                           const int rebinX = 4)
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetPadColor(kWhite);
  gStyle->SetFrameFillColor(kWhite);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetEndErrorSize(0);

  TFile input(inputRoot, "READ");
  if (input.IsZombie())
  {
    Error("MakePPG12Fig25Current", "Cannot open %s", inputRoot);
    return;
  }
  auto *source = dynamic_cast<TH2 *>(input.Get(histPath));
  if (!source)
  {
    Error("MakePPG12Fig25Current", "Missing %s", histPath);
    return;
  }

  auto *hist = dynamic_cast<TH2 *>(source->Clone("h2d_current_ppg12_fig25_profile"));
  hist->SetDirectory(nullptr);
  hist->RebinX(rebinX);  // Exact per-variable display/profile rebin from plot_showershapes.C.
  hist->GetXaxis()->SetRangeUser(0.0, 1.0);
  auto *profile = hist->ProfileX("pfx_current_ppg12_fig25_profile", 1, -1, "");
  profile->SetDirectory(nullptr);
  profile->SetLineColor(kBlue);
  profile->SetMarkerColor(kBlue);
  profile->SetMarkerStyle(20);
  profile->SetMarkerSize(0.45);
  profile->SetLineWidth(1);
  profile->SetTitle("");
  profile->GetXaxis()->SetRangeUser(0.0, 1.0);
  profile->GetXaxis()->SetTitle(xAxisTitle);
  profile->GetYaxis()->SetTitle("<#it{E}_{T}^{iso}> [GeV]");
  profile->GetYaxis()->SetTitleOffset(1.25);
  profile->GetXaxis()->SetNdivisions(505);

  const double ymax = 1.30 * profile->GetMaximum();
  profile->GetYaxis()->SetRangeUser(0.0, ymax > 0.0 ? ymax : 1.0);
  const double corr = source->GetCorrelationFactor();

  TCanvas canvas("c_ppg12_fig25_current", "PPG12 Fig. 25 current", 600, 600);
  canvas.SetLeftMargin(0.18);
  canvas.SetRightMargin(0.06);
  canvas.SetTopMargin(0.06);
  canvas.SetBottomMargin(0.15);
  profile->Draw("HIST");
  profile->Draw("EX0 SAME");

  TLatex text;
  text.SetNDC();
  text.SetTextFont(42);
  text.SetTextSize(0.042);
  text.DrawLatex(0.20, 0.94, "#it{#bf{sPHENIX}} Internal");
  text.DrawLatex(0.20, 0.89, "p+p #sqrt{s} = 200 GeV");
  text.DrawLatex(0.20, 0.84, "|#eta| < 0.7");
  // This lower-left region is empty for this profile; keep the full sample
  // annotation inside the panel without obscuring the curve or its errors.
  TPaveText legend(0.22, 0.21, 0.55, 0.45, "NDC");
  legend.SetFillColor(kWhite);
  legend.SetFillStyle(1001);
  legend.SetBorderSize(1);
  legend.SetTextAlign(12);
  legend.SetTextFont(42);
  legend.SetTextSize(0.030);
  legend.SetMargin(0.12);
  legend.AddText("14 < p_{T} < 18 GeV");
  legend.AddText("w/ nbkg cut");
  legend.AddText("Current combined");
  legend.AddText("background MC");
  legend.AddText(Form("Correlation: %.3f", corr));
  legend.Draw();

  gSystem->mkdir(gSystem->DirName(outputPng), kTRUE);
  canvas.SaveAs(outputPng);
  writeManifest(manifestPath, inputRoot, outputPng, currentPointer, rootSha256,
                histPath, rebinX, corr);

  Info("MakePPG12Fig25Current", "wrote %s (correlation %.6f)", outputPng, corr);
}
