#if defined(__CLING__)
R__ADD_INCLUDE_PATH(/Users/patsfan753/Desktop/ThesisAnalysis/macros)
#endif

#include "../../../ppg12codeGit/plotting/plotcommon.h"

#include <TCanvas.h>
#include <TFile.h>
#include <TH2.h>
#include <TSystem.h>

#include <cmath>
#include <iostream>
#include <string>

namespace
{
constexpr double kPriorParams[5] = {4.05592, -0.984728, -0.478818, 0.0723232, 0.0522681};

double ppg12TruthPriorWeight(double x)
{
  return (kPriorParams[0] + kPriorParams[1] * x + kPriorParams[3] * x * x) /
         (1.0 + kPriorParams[2] * x + kPriorParams[4] * x * x);
}

TH2 *cloneWithTruthPrior(TH2 *input, const char *name)
{
  TH2 *out = static_cast<TH2 *>(input->Clone(name));
  out->SetDirectory(nullptr);

  for (int ix = 1; ix <= out->GetNbinsX(); ++ix)
  {
    for (int iy = 1; iy <= out->GetNbinsY(); ++iy)
    {
      const double yCenter = out->GetYaxis()->GetBinCenter(iy);
      const double w = ppg12TruthPriorWeight(yCenter);
      out->SetBinContent(ix, iy, out->GetBinContent(ix, iy) * w);
      out->SetBinError(ix, iy, out->GetBinError(ix, iy) * std::abs(w));
    }
  }

  return out;
}

void drawOne(TH2 *h, const char *outPath, bool reweighted, bool ppg12ExactText)
{
  init_plot();

  h->SetDirectory(nullptr);
  h->SetXTitle(frame_response->GetXaxis()->GetTitle());
  h->SetYTitle(frame_response->GetYaxis()->GetTitle());
  h->SetMinimum(1.0);
  h->SetMaximum(1.0e7);

  TCanvas *c = new TCanvas("c_fig36_response", "c_fig36_response", 600, 600);
  c->SetLogz();
  gPad->SetRightMargin(0.15);
  h->Draw("colz");

  if (ppg12ExactText)
  {
    myText(0.5, 0.9, 1, strleg1.c_str(), 0.04);
    myText(0.5, 0.85, 1, strleg2.c_str(), 0.04);
    myText(0.3, 0.75, 1, "response matrix", 0.04);
  }
  else
  {
    myText(0.20, 0.90, 1, strleg1.c_str(), 0.04);
    myText(0.20, 0.85, 1, strleg2.c_str(), 0.04);
    myText(0.20, 0.80, 1, "response matrix", 0.04);
  }

  c->SaveAs(outPath);
  delete c;
}
}  // namespace

void MakeTHE76Fig36ResponseMatrixRootStyle(
    const char *inputRoot = "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/ppg12Parity/the76_ppg12_fig24_photonjet_fix_20260702_014217/final_roots/photonjet/RecoilJets_photonjet5plus10plus20_MERGED.root",
    const char *outDir = "/Users/patsfan753/Desktop/ThesisAnalysis/dataOutput/ppg12Parity/the76_ppg12_fig24_photonjet_fix_20260702_014217/fig36_response_matrix",
    const char *objectName = "SIM/h_response_full_0")
{
  gSystem->mkdir(outDir, true);

  TFile *fin = TFile::Open(inputRoot, "READ");
  if (!fin || fin->IsZombie())
  {
    std::cerr << "Could not open input ROOT: " << inputRoot << std::endl;
    return;
  }

  TH2 *raw = dynamic_cast<TH2 *>(fin->Get(objectName));
  if (!raw)
  {
    std::cerr << "Could not find object: " << objectName << " in " << inputRoot << std::endl;
    fin->Close();
    return;
  }

  TH2 *rawClone = static_cast<TH2 *>(raw->Clone("h_response_full_0_current_raw_rootstyle"));
  rawClone->SetDirectory(nullptr);
  TH2 *priorClone = cloneWithTruthPrior(raw, "h_response_full_0_current_offline_ppg12_prior_rootstyle");
  fin->Close();

  const std::string outBase(outDir);
  drawOne(rawClone,
          (outBase + "/fig36_current_h_response_full_raw_root_ppg12_exact_style.png").c_str(),
          false,
          true);
  drawOne(priorClone,
          (outBase + "/fig36_current_h_response_full_offline_truth_prior_reweighted_root_ppg12_exact_style.png").c_str(),
          true,
          true);
  drawOne(priorClone,
          (outBase + "/fig36_current_h_response_full_offline_truth_prior_reweighted_root_ppg12_lhs_text.png").c_str(),
          true,
          false);

  std::cout << "Wrote ROOT-style Fig36 response plots under " << outDir << std::endl;
}
