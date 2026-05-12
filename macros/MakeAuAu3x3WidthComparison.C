#include "sPhenixStyle.C"

#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TPad.h>
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace
{
const std::vector<std::pair<std::string, std::string>> kCentBins = {
  {"0_20", "0-20%"},
  {"20_50", "20-50%"},
  {"50_80", "50-80%"},
};

const int kCentColors[] = {kBlue + 1, kOrange + 7, kGreen + 2};

bool StartsWith(const std::string& s, const std::string& prefix)
{
  return s.size() >= prefix.size() && s.compare(0, prefix.size(), prefix) == 0;
}

bool EndsWith(const std::string& s, const std::string& suffix)
{
  return s.size() >= suffix.size() &&
         s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
}

std::unique_ptr<TH1> SumCentralityTemplates(TDirectory* dir,
                                            const std::string& var,
                                            const std::string& stage,
                                            const std::string& truthTag,
                                            const std::string& centKey,
                                            const std::string& name)
{
  if (!dir) return nullptr;

  const std::string prefix = "h_ss_" + var + "_" + stage + "_" + truthTag + "_pT_";
  const std::string suffix = "_cent_" + centKey;
  std::unique_ptr<TH1> sum;

  for (auto* obj : *dir->GetListOfKeys())
  {
    const std::string key = obj->GetName();
    if (!StartsWith(key, prefix) || !EndsWith(key, suffix)) continue;

    TH1* h = dynamic_cast<TH1*>(dir->Get(key.c_str()));
    if (!h) continue;

    if (!sum)
    {
      sum.reset(dynamic_cast<TH1*>(h->Clone(name.c_str())));
      sum->SetDirectory(nullptr);
      sum->Reset("ICES");
    }
    sum->Add(h);
  }

  if (!sum || sum->Integral() <= 0.0)
  {
    std::cerr << "[MakeAuAu3x3WidthComparison] missing/empty template for "
              << var << " " << stage << " " << truthTag << " cent=" << centKey << "\n";
    return sum;
  }

  sum->Scale(1.0 / sum->Integral());
  sum->SetStats(false);
  sum->SetLineWidth(4);
  sum->SetFillStyle(0);
  return sum;
}

void StyleHist(TH1* h, int iCent)
{
  if (!h) return;
  h->SetLineColor(kCentColors[iCent]);
  h->SetMarkerColor(kCentColors[iCent]);
  h->SetLineWidth(4);
  h->SetStats(false);
  h->GetXaxis()->SetLabelSize(0.060);
  h->GetYaxis()->SetLabelSize(0.060);
  h->GetXaxis()->SetTitleSize(0.070);
  h->GetYaxis()->SetTitleSize(0.070);
  h->GetXaxis()->SetTitleOffset(0.95);
  h->GetYaxis()->SetTitleOffset(1.05);
  h->GetXaxis()->SetNdivisions(505);
  h->GetYaxis()->SetNdivisions(505);
}

double MaxOf(const std::vector<std::unique_ptr<TH1>>& hs)
{
  double ymax = 0.0;
  for (const auto& h : hs)
  {
    if (h) ymax = std::max(ymax, h->GetMaximum());
  }
  return ymax;
}

void DrawPanel(TPad* pad,
               TDirectory* sigDir,
               TDirectory* bkgDir,
               const std::string& var,
               const std::string& stage,
               const std::string& truthTag,
               const std::string& panelTitle,
               const std::string& xTitle,
               const std::string& yTitle,
               bool drawSphenixBlock)
{
  pad->cd();
  pad->SetGrid(0, 0);
  pad->SetTicks(1, 1);
  pad->SetLeftMargin(0.145);
  pad->SetRightMargin(0.035);
  pad->SetBottomMargin(0.155);
  pad->SetTopMargin(0.145);

  TDirectory* dir = (truthTag == "sig") ? sigDir : bkgDir;
  std::vector<std::unique_ptr<TH1>> hs;
  hs.reserve(kCentBins.size());

  for (size_t i = 0; i < kCentBins.size(); ++i)
  {
    auto h = SumCentralityTemplates(dir, var, stage, truthTag, kCentBins[i].first,
                                    Form("h_%s_%s_%s_%s", var.c_str(), stage.c_str(),
                                         truthTag.c_str(), kCentBins[i].first.c_str()));
    if (h)
    {
      StyleHist(h.get(), static_cast<int>(i));
      h->GetXaxis()->SetTitle(xTitle.c_str());
      h->GetYaxis()->SetTitle(yTitle.c_str());
      h->GetXaxis()->SetRangeUser(0.0, 0.70);
    }
    hs.push_back(std::move(h));
  }

  const double ymax = MaxOf(hs);
  bool first = true;
  for (const auto& h : hs)
  {
    if (!h) continue;
    h->GetYaxis()->SetRangeUser(0.0, std::max(0.001, ymax * 1.36));
    h->Draw(first ? "hist" : "hist same");
    first = false;
  }

  TLatex tx;
  tx.SetNDC();
  tx.SetTextFont(42);
  tx.SetTextAlign(23);
  tx.SetTextSize(0.072);
  tx.DrawLatex(0.50, 0.982, panelTitle.c_str());

  auto* leg = new TLegend(drawSphenixBlock ? 0.58 : 0.57,
                          drawSphenixBlock ? 0.52 : 0.58,
                          0.94,
                          drawSphenixBlock ? 0.70 : 0.77);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(drawSphenixBlock ? 0.046 : 0.054);
  for (size_t i = 0; i < kCentBins.size(); ++i)
  {
    if (hs[i]) leg->AddEntry(hs[i].get(), kCentBins[i].second.c_str(), "l");
  }
  leg->Draw();

  if (drawSphenixBlock)
  {
    tx.SetTextAlign(13);
    tx.SetTextSize(0.044);
    tx.DrawLatex(0.58, 0.835, "#it{#bf{sPHENIX}} Internal");
    tx.SetTextSize(0.032);
    tx.DrawLatex(0.58, 0.780, "Pythia overlay, #sqrt{s_{NN}}=200 GeV");
  }

  // ROOT canvases keep object pointers, so keep these cloned templates alive
  // after this helper returns.
  for (auto& h : hs) h.release();

  pad->Modified();
}

void DrawComparison(const std::string& reportDir,
                    const std::string& outDir,
                    const std::string& stage,
                    const std::string& stageLabel,
                    const std::vector<std::string>& subtitleLines,
                    const std::string& outStem)
{
  const std::string sigPath =
    "dataOutput/combinedSimOnlyEMBEDDED/preselectionNewPPG12_tightAuAuCentInputBDT_nonTightAuAuBDTComplement_baseVariant/photonJet12and20merged_SIM/RecoilJets_embeddedPhoton12plus20_MERGED.root";
  const std::string bkgPath =
    "dataOutput/combinedSimOnlyEMBEDDED/preselectionNewPPG12_tightAuAuCentInputBDT_nonTightAuAuBDTComplement_baseVariant/embeddedJet12and20merged_SIM/RecoilJets_embeddedJet12plus20_MERGED.root";

  std::unique_ptr<TFile> sigFile(TFile::Open(sigPath.c_str(), "READ"));
  std::unique_ptr<TFile> bkgFile(TFile::Open(bkgPath.c_str(), "READ"));
  if (!sigFile || sigFile->IsZombie() || !bkgFile || bkgFile->IsZombie())
  {
    std::cerr << "[MakeAuAu3x3WidthComparison] cannot open input files\n";
    return;
  }

  TDirectory* sigDir = dynamic_cast<TDirectory*>(sigFile->Get("SIM"));
  TDirectory* bkgDir = dynamic_cast<TDirectory*>(bkgFile->Get("SIM"));
  if (!sigDir || !bkgDir)
  {
    std::cerr << "[MakeAuAu3x3WidthComparison] missing SIM directory\n";
    return;
  }

  TCanvas c(Form("c_%s", outStem.c_str()), "", 2200, 1350);
  c.SetFillColor(0);
  c.SetBorderMode(0);

  TPad p1("p1", "", 0.055, 0.500, 0.495, 0.825);
  TPad p2("p2", "", 0.555, 0.500, 0.995, 0.825);
  TPad p3("p3", "", 0.055, 0.080, 0.495, 0.405);
  TPad p4("p4", "", 0.555, 0.080, 0.995, 0.405);
  p1.Draw();
  p2.Draw();
  p3.Draw();
  p4.Draw();

  DrawPanel(&p1, sigDir, bkgDir, "weta", stage, "sig",
            "Signal: full-cluster width",
            "w_{#eta}^{COG}", "Unit-normalized candidates",
            true);
  DrawPanel(&p2, sigDir, bkgDir, "weta33", stage, "sig",
            "Signal: 3#times3 tower-window width",
            "w_{#eta}^{3#times3,COG}", "Unit-normalized candidates",
            false);
  DrawPanel(&p3, sigDir, bkgDir, "weta", stage, "bkg",
            "Background: full-cluster width",
            "w_{#eta}^{COG}", "Unit-normalized candidates",
            false);
  DrawPanel(&p4, sigDir, bkgDir, "weta33", stage, "bkg",
            "Background: 3#times3 tower-window width",
            "w_{#eta}^{3#times3,COG}", "Unit-normalized candidates",
            false);

  c.cd();
  TLatex title;
  title.SetNDC();
  title.SetTextFont(62);
  title.SetTextAlign(23);
  title.SetTextSize(0.033);
  title.DrawLatex(0.50, 0.982, stageLabel.c_str());
  title.SetTextFont(42);
  title.SetTextSize(0.023);
  double ySubtitle = 0.940;
  for (const auto& line : subtitleLines)
  {
    title.DrawLatex(0.50, ySubtitle, line.c_str());
    ySubtitle -= 0.032;
  }

  gSystem->mkdir(outDir.c_str(), true);
  const std::string png = outDir + "/" + outStem + ".png";
  const std::string root = outDir + "/" + outStem + ".root";
  c.SaveAs(png.c_str());
  c.SaveAs(root.c_str());
  std::cout << "[MakeAuAu3x3WidthComparison] wrote " << png << "\n";
  std::cout << "[MakeAuAu3x3WidthComparison] wrote " << root << "\n";
}
}

void MakeAuAu3x3WidthComparison(
  const char* reportDir = "dataOutput/auauTightBDTValidation/model_validation_condor_cent3x3_20260510_221419",
  const char* outDir = "dataOutput/auauTightBDTValidation/model_validation_condor_cent3x3_20260510_221419/money_plots_3x3")
{
  SetsPhenixStyle();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  const std::string report(reportDir);
  const std::string out(outDir);

  DrawComparison(report, out, "pre",
                 "Full-cluster vs 3#times3 shower-width templates after common photon preselection",
                 {
                   "Preselection: NPB > 0.5; 0.6 < e_{T1} < 1.0; E_{11}/E_{33} < 0.98",
                   "0.8 < E_{32}/E_{35} < 1.0; no w_{#eta} cut; no tight-BDT requirement",
                 },
                 "shower_width_weta_vs_3x3_centrality_preselected_2x2_root_v4");

  DrawComparison(report, out, "inclusive",
                 "Full-cluster vs 3#times3 shower-width templates before common photon preselection",
                 {
                   "Before common preselection; truth split from embedded MC; no tight-BDT requirement",
                 },
                 "shower_width_weta_vs_3x3_centrality_inclusive_2x2_root_v4");
}
