#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TKey.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TNamed.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace
{
const std::string kMergedBaseDir = "dataOutput/combinedSimOnlyEMBEDDED";
const std::string kPPMergedBaseDir = "dataOutput/combinedSimOnly";
const std::string kOutDir = "dataOutput/combinedSimOnlyEMBEDDED/photonEfficiencyCompare";
const std::string kPPOutDir = "dataOutput/pp";
const std::string kTopDir = "SIM";
const std::string kMergedFilename = "RecoilJets_embeddedPhoton12plus20_MERGED.root";
const std::string kPPMergedFilename = "RecoilJets_photonjet5plus10plus20_MERGED.root";

const std::string kReferencePreselectionLabel = "reference preselection";
const std::string kVariantAPreselectionLabel = "New PPG12 Presel";
const std::string kVariantCPreselectionLabel = "NPB Only Presel";
const std::string kVariantDPreselectionLabel = "reference + NPB presel";
const std::string kNoPreselectionLabel = "no preselection";

struct Variant
{
  std::string tag;
  std::string label;
};

struct CentBin
{
  std::string label;
  std::string suffix;
  int color = kBlack;
};

struct EffPoint
{
  double lo = 0.0;
  double hi = 0.0;
  double truth = 0.0;
  double miss = 0.0;
  double eTruth = 0.0;
  double eMiss = 0.0;
  double eff = 0.0;
  double eEff = 0.0;
};

std::ofstream gSummary;

const std::vector<CentBin>& CentBins()
{
  static const std::vector<CentBin> bins = {
      {"0-20%", "_cent_0_20", kBlack},
      {"20-50%", "_cent_20_50", kBlue + 1},
      {"50-80%", "_cent_50_80", kOrange + 7},
  };
  return bins;
}

void Log(const std::string& line)
{
  std::cout << line << "\n";
  if (gSummary.is_open()) gSummary << line << "\n";
}

std::string MergedOutDir(const std::string& tag)
{
  return kMergedBaseDir + "/" + tag + "/photonJet12and20merged_SIM";
}

std::string MergedFileName(const std::string& tag)
{
  return MergedOutDir(tag) + "/" + kMergedFilename;
}

std::string PPMergedOutDir(const std::string& tag)
{
  return kPPMergedBaseDir + "/" + tag + "/photonJet5and10and20merged_SIM";
}

std::string PPMergedFileName(const std::string& tag)
{
  return PPMergedOutDir(tag) + "/" + kPPMergedFilename;
}

bool RequireMergedFile(const std::string& tag)
{
  const std::string path = MergedFileName(tag);
  if (gSystem && !gSystem->AccessPathName(path.c_str())) return true;
  Log("[ERROR] Missing canonical merged embedded photon SIM file for " + tag);
  Log("        expected: " + path);
  Log("        build it with: ./scripts/sftp_get_recoiljets_outputs.sh mergeLocalSim isSimEmbedded");
  return false;
}

bool RequirePPMergedFile(const std::string& tag)
{
  const std::string path = PPMergedFileName(tag);
  if (gSystem && !gSystem->AccessPathName(path.c_str())) return true;
  Log("[ERROR] Missing canonical merged pp photon SIM file for " + tag);
  Log("        expected: " + path);
  Log("        build it with: ./scripts/sftp_get_recoiljets_outputs.sh mergeLocalSim isSim");
  return false;
}

std::vector<EffPoint> ReadEfficiencyPoints(TFile* f, const Variant& variant, const CentBin& cent)
{
  std::vector<EffPoint> pts;
  if (!f) return pts;

  TDirectory* d = f->GetDirectory(kTopDir.c_str());
  if (!d) return pts;

  TH1* hTruth = dynamic_cast<TH1*>(d->Get(("h_unfoldTruthPho_pTgamma" + cent.suffix).c_str()));
  TH1* hMiss = dynamic_cast<TH1*>(d->Get(("h_unfoldTruthPhoMisses_pTgamma" + cent.suffix).c_str()));
  if (!hTruth || !hMiss)
  {
    Log("[WARN] Missing photon efficiency histograms for " + variant.label + " " + cent.label);
    return pts;
  }

  Log("");
  Log("[PHOTON EFFICIENCY] " + variant.label + "  centrality " + cent.label);
  Log("  file: " + MergedFileName(variant.tag));
  Log("  pTbin        N_truth       N_miss       epsilon");
  Log("  ------------------------------------------------");

  const int nb = hTruth->GetNbinsX();
  for (int ib = 1; ib <= nb; ++ib)
  {
    const double lo = hTruth->GetXaxis()->GetBinLowEdge(ib);
    const double hi = hTruth->GetXaxis()->GetBinUpEdge(ib);
    if (lo < 10.0 || lo >= 35.0) continue;

    const int im = hMiss->GetXaxis()->FindBin(0.5 * (lo + hi));
    if (im < 1 || im > hMiss->GetNbinsX()) continue;

    EffPoint p;
    p.lo = lo;
    p.hi = hi;
    p.truth = hTruth->GetBinContent(ib);
    p.eTruth = hTruth->GetBinError(ib);
    p.miss = hMiss->GetBinContent(im);
    p.eMiss = hMiss->GetBinError(im);
    if (p.truth <= 0.0) continue;

    const double r = p.miss / p.truth;
    p.eff = 1.0 - r;
    double er = 0.0;
    if (p.miss > 0.0)
    {
      er = r * std::sqrt(std::pow(p.eMiss / p.miss, 2) + std::pow(p.eTruth / p.truth, 2));
    }
    p.eEff = er;

    std::ostringstream line;
    line << "  " << std::setw(5) << std::right << std::fixed << std::setprecision(0) << lo
         << "-" << std::setw(2) << hi
         << std::setw(14) << std::setprecision(3) << p.truth
         << std::setw(13) << p.miss
         << std::setw(14) << std::setprecision(6) << p.eff;
    Log(line.str());

    pts.push_back(p);
  }

  return pts;
}

std::vector<EffPoint> ReadPPEfficiencyPoints(TFile* f, const Variant& variant)
{
  std::vector<EffPoint> pts;
  if (!f) return pts;

  TDirectory* d = f->GetDirectory(kTopDir.c_str());
  if (!d) return pts;

  TH1* hTruth = dynamic_cast<TH1*>(d->Get("h_unfoldTruthPho_pTgamma_ppg12obj"));
  TH1* hMiss = dynamic_cast<TH1*>(d->Get("h_unfoldTruthPhoMisses_pTgamma_ppg12obj"));
  if (!hTruth || !hMiss)
  {
    Log("[WARN] Missing pp photon efficiency histograms for " + variant.label);
    return pts;
  }

  Log("");
  Log("[PP PHOTON EFFICIENCY] " + variant.label);
  Log("  file: " + PPMergedFileName(variant.tag));
  Log("  pTbin        N_truth       N_miss       epsilon");
  Log("  ------------------------------------------------");

  const int nb = hTruth->GetNbinsX();
  for (int ib = 1; ib <= nb; ++ib)
  {
    const double lo = hTruth->GetXaxis()->GetBinLowEdge(ib);
    const double hi = hTruth->GetXaxis()->GetBinUpEdge(ib);
    if (lo < 10.0 || lo >= 35.0) continue;

    const int im = hMiss->GetXaxis()->FindBin(0.5 * (lo + hi));
    if (im < 1 || im > hMiss->GetNbinsX()) continue;

    EffPoint p;
    p.lo = lo;
    p.hi = hi;
    p.truth = hTruth->GetBinContent(ib);
    p.eTruth = hTruth->GetBinError(ib);
    p.miss = hMiss->GetBinContent(im);
    p.eMiss = hMiss->GetBinError(im);
    if (p.truth <= 0.0) continue;

    const double r = p.miss / p.truth;
    p.eff = 1.0 - r;
    p.eEff = (p.miss > 0.0)
                 ? r * std::sqrt(std::pow(p.eMiss / p.miss, 2) + std::pow(p.eTruth / p.truth, 2))
                 : 0.0;

    std::ostringstream line;
    line << "  " << std::setw(5) << std::right << std::fixed << std::setprecision(0) << lo
         << "-" << std::setw(2) << hi
         << std::setw(14) << std::setprecision(3) << p.truth
         << std::setw(13) << p.miss
         << std::setw(14) << std::setprecision(6) << p.eff;
    Log(line.str());

    pts.push_back(p);
  }

  return pts;
}

TGraphErrors* BuildEfficiencyGraph(TFile* f, const Variant& variant, const CentBin& cent,
                                   bool openMarker, double xShift = 0.0)
{
  const std::vector<EffPoint> pts = ReadEfficiencyPoints(f, variant, cent);
  if (pts.empty()) return nullptr;

  std::vector<double> x, ex, y, ey;
  x.reserve(pts.size());
  ex.reserve(pts.size());
  y.reserve(pts.size());
  ey.reserve(pts.size());

  for (const EffPoint& p : pts)
  {
    x.push_back(0.5 * (p.lo + p.hi) + xShift);
    ex.push_back(0.5 * (p.hi - p.lo));
    y.push_back(p.eff);
    ey.push_back(p.eEff);
  }

  TGraphErrors* g = new TGraphErrors((int)x.size(), x.data(), y.data(), ex.data(), ey.data());
  g->SetName(("gPhoEff_" + variant.label + "_" + cent.label).c_str());
  g->SetLineColor(cent.color);
  g->SetMarkerColor(cent.color);
  g->SetMarkerStyle(openMarker ? 24 : 20);
  g->SetMarkerSize(1.1);
  g->SetLineWidth(2);
  g->SetLineStyle(openMarker ? 2 : 1);
  return g;
}

TGraphErrors* BuildPPEfficiencyGraph(TFile* f, const Variant& variant,
                                     bool openMarker, double xShift = 0.0,
                                     int color = kRed + 1)
{
  const std::vector<EffPoint> pts = ReadPPEfficiencyPoints(f, variant);
  if (pts.empty()) return nullptr;

  std::vector<double> x, ex, y, ey;
  x.reserve(pts.size());
  ex.reserve(pts.size());
  y.reserve(pts.size());
  ey.reserve(pts.size());

  for (const EffPoint& p : pts)
  {
    x.push_back(0.5 * (p.lo + p.hi) + xShift);
    ex.push_back(0.5 * (p.hi - p.lo));
    y.push_back(p.eff);
    ey.push_back(p.eEff);
  }

  TGraphErrors* g = new TGraphErrors((int)x.size(), x.data(), y.data(), ex.data(), ey.data());
  g->SetName(("gPhoEffPP_" + variant.label).c_str());
  g->SetLineColor(color);
  g->SetMarkerColor(color);
  g->SetMarkerStyle(openMarker ? 24 : 20);
  g->SetMarkerSize(1.05);
  g->SetLineWidth(2);
  g->SetLineStyle(openMarker ? 2 : 1);
  return g;
}

std::map<std::string, EffPoint> PointMap(TFile* f, const Variant& variant, const CentBin& cent)
{
  std::map<std::string, EffPoint> m;
  for (const EffPoint& p : ReadEfficiencyPoints(f, variant, cent))
  {
    std::ostringstream key;
    key << std::fixed << std::setprecision(0) << p.lo << "_" << p.hi;
    m[key.str()] = p;
  }
  return m;
}

double RatioError(double num, double enum_, double den, double eden)
{
  if (num <= 0.0 || den <= 0.0) return 0.0;
  const double r = num / den;
  return r * std::sqrt(std::pow(enum_ / num, 2) + std::pow(eden / den, 2));
}

TGraphErrors* BuildRatioGraph(TFile* fNum, const Variant& vNum,
                              TFile* fDen, const Variant& vDen,
                              const CentBin& cent,
                              int color, int marker, bool openMarker,
                              double xShift = 0.0)
{
  const auto nmap = PointMap(fNum, vNum, cent);
  const auto dmap = PointMap(fDen, vDen, cent);

  std::vector<double> x, ex, y, ey;
  for (const auto& kv : nmap)
  {
    auto it = dmap.find(kv.first);
    if (it == dmap.end()) continue;
    const EffPoint& n = kv.second;
    const EffPoint& d = it->second;
    if (d.eff <= 0.0) continue;
    x.push_back(0.5 * (n.lo + n.hi) + xShift);
    ex.push_back(0.0);
    y.push_back(n.eff / d.eff);
    ey.push_back(RatioError(n.eff, n.eEff, d.eff, d.eEff));
  }
  if (x.empty()) return nullptr;

  TGraphErrors* g = new TGraphErrors((int)x.size(), x.data(), y.data(), ex.data(), ey.data());
  g->SetLineColor(color);
  g->SetMarkerColor(color);
  g->SetMarkerStyle(openMarker ? 24 : marker);
  g->SetMarkerSize(1.05);
  g->SetLineWidth(2);
  return g;
}

void DrawHeader(const std::string& title,
                double xRight = 0.92,
                double yTop = 0.89,
                double yExperiment = 0.675,
                double yCollision = 0.635,
                int experimentAlign = 31,
                double titleSize = 0.030)
{
  TLatex tx;
  tx.SetNDC();
  tx.SetTextFont(42);
  tx.SetTextAlign(22);
  tx.SetTextSize(titleSize);
  tx.DrawLatex(0.50, yTop, title.c_str());

  tx.SetTextAlign(experimentAlign);
  tx.SetTextSize(0.032);
  tx.DrawLatex(xRight, yExperiment, "#bf{sPHENIX} #it{Internal}");
  tx.DrawLatex(xRight, yCollision, "Pythia Overlay  #sqrt{s_{NN}} = 200 GeV");
}

void DrawEfficiencyPair(const std::string& title,
                        const std::string& outStem,
                        const Variant& first,
                        const Variant& second,
                        const std::string& firstLegend,
                        const std::string& secondLegend,
                        const std::string& outputSubdir,
                        bool forceNoPPReferenceOverlays = false)
{
  if (!RequireMergedFile(first.tag)) return;
  if (!RequireMergedFile(second.tag)) return;

  std::unique_ptr<TFile> fFirst(TFile::Open(MergedFileName(first.tag).c_str(), "READ"));
  std::unique_ptr<TFile> fSecond(TFile::Open(MergedFileName(second.tag).c_str(), "READ"));
  if (!fFirst || fFirst->IsZombie() || !fSecond || fSecond->IsZombie()) return;

  const bool addPPReferenceOverlays = !forceNoPPReferenceOverlays &&
                                      (outputSubdir == "04_reference_plus_npb_id_overlays" ||
                                       outputSubdir == "04_reference_preselection_id_overlays" ||
                                       outputSubdir == "04_variantA_preselection_id_overlays");
  const Variant ppReference{
      "jetMinPt5_7pi_8_vz60_isoR40_fixedIso2GeV_preselectionReference_tightReference_nonTightReference",
      "pp Reference/Reference/Reference"};
  const Variant ppVariantA{
      "jetMinPt5_7pi_8_vz60_isoR40_fixedIso2GeV_preselectionVariantA_tightVariantA_nonTightVariantA",
      "pp VariantA/VariantA/VariantA"};
  std::unique_ptr<TFile> fPPReference;
  std::unique_ptr<TFile> fPPVariantA;
  if (addPPReferenceOverlays)
  {
    if (!RequirePPMergedFile(ppReference.tag)) return;
    if (!RequirePPMergedFile(ppVariantA.tag)) return;
    fPPReference.reset(TFile::Open(PPMergedFileName(ppReference.tag).c_str(), "READ"));
    fPPVariantA.reset(TFile::Open(PPMergedFileName(ppVariantA.tag).c_str(), "READ"));
    if (!fPPReference || fPPReference->IsZombie() || !fPPVariantA || fPPVariantA->IsZombie()) return;
  }

  std::vector<std::unique_ptr<TGraphErrors>> graphs;
  for (const CentBin& cent : CentBins())
  {
    graphs.emplace_back(BuildEfficiencyGraph(fFirst.get(), first, cent, false));
    graphs.emplace_back(BuildEfficiencyGraph(fSecond.get(), second, cent, true));
  }
  if (addPPReferenceOverlays)
  {
    graphs.emplace_back(BuildPPEfficiencyGraph(fPPReference.get(), ppReference, false, -0.08));
    graphs.emplace_back(BuildPPEfficiencyGraph(fPPVariantA.get(), ppVariantA, true, 0.08));
  }

  double ymin = 1.0;
  double ymax = 0.0;
  for (const auto& owned : graphs)
  {
    TGraphErrors* g = owned.get();
    if (!g) continue;
    for (int i = 0; i < g->GetN(); ++i)
    {
      double gx = 0.0, gy = 0.0;
      g->GetPoint(i, gx, gy);
      ymin = std::min(ymin, gy - g->GetErrorY(i));
      ymax = std::max(ymax, gy + g->GetErrorY(i));
    }
  }
  if (ymax <= ymin)
  {
    ymin = 0.0;
    ymax = 1.0;
  }
  ymin = std::max(0.0, ymin - 0.08);
  ymax = std::min(1.60, std::max(addPPReferenceOverlays ? 1.42 : 1.08,
                                  ymax + (addPPReferenceOverlays ? 0.48 : 0.22)));

  TCanvas c(("c_" + outStem).c_str(), ("c_" + outStem).c_str(), 1100, 850);
  c.SetLeftMargin(0.14);
  c.SetRightMargin(0.04);
  c.SetTopMargin(0.13);
  c.SetBottomMargin(0.14);
  c.SetTicks(1, 1);

  TH1F frame(("hframe_" + outStem).c_str(), "", 100, 10.0, 35.0);
  frame.SetDirectory(nullptr);
  frame.SetStats(0);
  frame.SetMinimum(ymin);
  frame.SetMaximum(ymax);
  frame.GetXaxis()->SetTitle("p_{T}^{#gamma,truth} [GeV]");
  frame.GetYaxis()->SetTitle("#epsilon_{#gamma}^{MC} = 1 - N_{miss}^{truth}/N_{#gamma}^{truth}");
  frame.GetXaxis()->SetTitleSize(0.048);
  frame.GetYaxis()->SetTitleSize(0.048);
  frame.GetXaxis()->SetLabelSize(0.040);
  frame.GetYaxis()->SetLabelSize(0.040);
  frame.Draw();

  for (const auto& owned : graphs)
  {
    if (owned) owned->Draw("P SAME");
  }

  TGraphErrors legFirst(1);
  legFirst.SetMarkerStyle(20);
  legFirst.SetMarkerColor(kBlack);
  legFirst.SetLineColor(kBlack);
  legFirst.SetLineWidth(2);
  TGraphErrors legSecond(1);
  legSecond.SetMarkerStyle(24);
  legSecond.SetMarkerColor(kBlack);
  legSecond.SetLineColor(kBlack);
  legSecond.SetLineStyle(2);
  legSecond.SetLineWidth(2);

  const bool longEmbeddedOnlyLegend =
      forceNoPPReferenceOverlays &&
      (outStem.find("originalBoxCutsPlusNPBPresel") != std::string::npos ||
       outStem.find("refPreselBDTTightLooseNT") != std::string::npos);
  const double varLegY1 = addPPReferenceOverlays ? 0.61 : 0.72;
  const double varLegX1 = addPPReferenceOverlays ? 0.49 : (longEmbeddedOnlyLegend ? 0.42 : 0.58);
  TLegend varLeg(varLegX1, varLegY1, 0.985, 0.87);
  varLeg.SetBorderSize(0);
  varLeg.SetFillStyle(0);
  varLeg.SetTextFont(42);
  varLeg.SetTextSize(addPPReferenceOverlays ? 0.0215 : (longEmbeddedOnlyLegend ? 0.024 : 0.030));
  varLeg.AddEntry(&legFirst, firstLegend.c_str(), "pe");
  varLeg.AddEntry(&legSecond, secondLegend.c_str(), "pe");

  TGraphErrors legPPReference(1);
  legPPReference.SetMarkerStyle(20);
  legPPReference.SetMarkerColor(kRed + 1);
  legPPReference.SetLineColor(kRed + 1);
  legPPReference.SetLineWidth(2);
  TGraphErrors legPPVariantA(1);
  legPPVariantA.SetMarkerStyle(24);
  legPPVariantA.SetMarkerColor(kRed + 1);
  legPPVariantA.SetLineColor(kRed + 1);
  legPPVariantA.SetLineStyle(2);
  legPPVariantA.SetLineWidth(2);
  if (addPPReferenceOverlays)
  {
    varLeg.AddEntry(&legPPReference, "MC original PPG12 box-cuts", "pe");
    varLeg.AddEntry(&legPPVariantA, "MC updated PPG12 BDT tight/nontight", "pe");
  }
  varLeg.Draw();

  std::vector<std::unique_ptr<TGraphErrors>> centLegendGraphs;
  TLegend centLeg(0.18, 0.72, 0.39, 0.86);
  centLeg.SetBorderSize(0);
  centLeg.SetFillStyle(0);
  centLeg.SetTextFont(42);
  centLeg.SetTextSize(0.030);
  for (const CentBin& cent : CentBins())
  {
    std::unique_ptr<TGraphErrors> g(new TGraphErrors(1));
    g->SetMarkerStyle(20);
    g->SetMarkerColor(cent.color);
    g->SetLineColor(cent.color);
    g->SetLineWidth(2);
    centLeg.AddEntry(g.get(), cent.label.c_str(), "pe");
    centLegendGraphs.push_back(std::move(g));
  }
  centLeg.Draw();

  DrawHeader(title, 0.22, 0.89, 0.685, 0.640, 11, addPPReferenceOverlays ? 0.023 : 0.030);

  const std::string outDir = kOutDir + "/" + outputSubdir;
  gSystem->mkdir(outDir.c_str(), true);
  const std::string outPng = outDir + "/" + outStem + ".png";
  c.SaveAs(outPng.c_str());
  Log("[DONE] Wrote " + outPng);
}

void DrawPPFixedIsoEfficiencyPair()
{
  const Variant boxCuts{
      "jetMinPt5_7pi_8_vz60_isoR40_fixedIso2GeV_preselectionReference_tightReference_nonTightReference",
      "Box-cuts"};
  const Variant bdtNpb{
      "jetMinPt5_7pi_8_vz60_isoR40_fixedIso2GeV_preselectionVariantA_tightVariantA_nonTightVariantA",
      "BDT + NPB"};

  if (!RequirePPMergedFile(boxCuts.tag)) return;
  if (!RequirePPMergedFile(bdtNpb.tag)) return;

  std::unique_ptr<TFile> fBox(TFile::Open(PPMergedFileName(boxCuts.tag).c_str(), "READ"));
  std::unique_ptr<TFile> fBdt(TFile::Open(PPMergedFileName(bdtNpb.tag).c_str(), "READ"));
  if (!fBox || fBox->IsZombie() || !fBdt || fBdt->IsZombie()) return;

  std::unique_ptr<TGraphErrors> gBox(BuildPPEfficiencyGraph(fBox.get(), boxCuts, false, -0.06, kBlack));
  std::unique_ptr<TGraphErrors> gBdt(BuildPPEfficiencyGraph(fBdt.get(), bdtNpb, false, 0.06, kBlue + 1));
  if (!gBox && !gBdt) return;
  if (gBox) gBox->SetMarkerStyle(20);
  if (gBdt) gBdt->SetMarkerStyle(21);

  double ymin = 1.0;
  double ymax = 0.0;
  for (TGraphErrors* g : {gBox.get(), gBdt.get()})
  {
    if (!g) continue;
    for (int i = 0; i < g->GetN(); ++i)
    {
      double x = 0.0, y = 0.0;
      g->GetPoint(i, x, y);
      ymin = std::min(ymin, y - g->GetErrorY(i));
      ymax = std::max(ymax, y + g->GetErrorY(i));
    }
  }
  if (ymax <= ymin)
  {
    ymin = 0.0;
    ymax = 1.0;
  }
  ymin = std::max(0.0, ymin - 0.08);
  ymax = std::min(1.30, std::max(1.05, ymax + 0.16));

  TCanvas c("c_ppFixedIsoPhotonEfficiencyReferenceVsVariantA",
            "c_ppFixedIsoPhotonEfficiencyReferenceVsVariantA", 1100, 850);
  c.SetLeftMargin(0.14);
  c.SetRightMargin(0.04);
  c.SetTopMargin(0.13);
  c.SetBottomMargin(0.14);
  c.SetTicks(1, 1);

  TH1F frame("hframe_ppFixedIsoPhotonEfficiencyReferenceVsVariantA", "", 100, 10.0, 35.0);
  frame.SetDirectory(nullptr);
  frame.SetStats(0);
  frame.SetMinimum(ymin);
  frame.SetMaximum(ymax);
  frame.GetXaxis()->SetTitle("p_{T}^{#gamma,truth} [GeV]");
  frame.GetYaxis()->SetTitle("#epsilon_{#gamma}^{MC} = 1 - N_{miss}^{truth}/N_{#gamma}^{truth}");
  frame.GetXaxis()->SetTitleSize(0.048);
  frame.GetYaxis()->SetTitleSize(0.044);
  frame.GetXaxis()->SetLabelSize(0.040);
  frame.GetYaxis()->SetLabelSize(0.040);
  frame.Draw();

  if (gBox) gBox->Draw("P SAME");
  if (gBdt) gBdt->Draw("P SAME");

  TLegend leg(0.66, 0.72, 0.91, 0.86);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextFont(42);
  leg.SetTextSize(0.032);
  if (gBox) leg.AddEntry(gBox.get(), "Box-cuts", "pe");
  if (gBdt) leg.AddEntry(gBdt.get(), "BDT + NPB", "pe");
  leg.Draw();

  TLatex tx;
  tx.SetNDC();
  tx.SetTextFont(42);
  tx.SetTextAlign(22);
  tx.SetTextSize(0.034);
  tx.DrawLatex(0.50, 0.92, "MC #gamma-efficiency - Run24pp SIM");

  tx.SetTextAlign(13);
  tx.SetTextSize(0.030);
  tx.DrawLatex(0.18, 0.845, "Photon 4 GeV + MBD N&S #geq 1");
  tx.DrawLatex(0.18, 0.800, "|v_{z}| < 60 cm, #Delta R_{cone} < 0.4, E_{T}^{iso} < 2 GeV");

  tx.SetTextAlign(31);
  tx.SetTextSize(0.030);
  tx.DrawLatex(0.91, 0.655, "#bf{sPHENIX} #it{Internal}");
  tx.DrawLatex(0.91, 0.615, "Pythia  #sqrt{s} = 200 GeV");

  gSystem->mkdir(kPPOutDir.c_str(), true);
  const std::string outPng = kPPOutDir + "/phoEff_overlay_fixedIso2GeV_reference_vs_variantA.png";
  c.SaveAs(outPng.c_str());
  Log("[DONE] Wrote " + outPng);
}

void DrawRatioSingle(const std::string& title,
                     const std::string& outStem,
                     const Variant& numerator,
                     const Variant& denominator,
                     const std::string& legendText)
{
  if (!RequireMergedFile(numerator.tag)) return;
  if (!RequireMergedFile(denominator.tag)) return;

  std::unique_ptr<TFile> fNum(TFile::Open(MergedFileName(numerator.tag).c_str(), "READ"));
  std::unique_ptr<TFile> fDen(TFile::Open(MergedFileName(denominator.tag).c_str(), "READ"));
  if (!fNum || fNum->IsZombie() || !fDen || fDen->IsZombie()) return;

  std::vector<std::unique_ptr<TGraphErrors>> graphs;
  for (const CentBin& cent : CentBins())
  {
    graphs.emplace_back(BuildRatioGraph(fNum.get(), numerator, fDen.get(), denominator,
                                        cent, cent.color, 20, false));
  }

  TCanvas c(("c_" + outStem).c_str(), ("c_" + outStem).c_str(), 1100, 850);
  c.SetLeftMargin(0.14);
  c.SetRightMargin(0.04);
  c.SetTopMargin(0.13);
  c.SetBottomMargin(0.14);
  c.SetTicks(1, 1);

  TH1F frame(("hframe_" + outStem).c_str(), "", 100, 10.0, 35.0);
  frame.SetDirectory(nullptr);
  frame.SetStats(0);
  frame.SetMinimum(0.80);
  frame.SetMaximum(1.20);
  frame.GetXaxis()->SetTitle("p_{T}^{#gamma,truth} [GeV]");
  frame.GetYaxis()->SetTitle("#epsilon_{#gamma}^{numerator} / #epsilon_{#gamma}^{denominator}");
  frame.GetXaxis()->SetTitleSize(0.048);
  frame.GetYaxis()->SetTitleSize(0.048);
  frame.GetXaxis()->SetLabelSize(0.040);
  frame.GetYaxis()->SetLabelSize(0.040);
  frame.Draw();

  TLine one(10.0, 1.0, 35.0, 1.0);
  one.SetLineStyle(2);
  one.SetLineColor(kGray + 2);
  one.Draw("same");

  for (const auto& owned : graphs)
  {
    if (owned) owned->Draw("P SAME");
  }

  std::vector<std::unique_ptr<TGraphErrors>> centLegendGraphs;
  TLegend centLeg(0.18, 0.72, 0.39, 0.86);
  centLeg.SetBorderSize(0);
  centLeg.SetFillStyle(0);
  centLeg.SetTextFont(42);
  centLeg.SetTextSize(0.030);
  for (const CentBin& cent : CentBins())
  {
    std::unique_ptr<TGraphErrors> g(new TGraphErrors(1));
    g->SetMarkerStyle(20);
    g->SetMarkerColor(cent.color);
    g->SetLineColor(cent.color);
    g->SetLineWidth(2);
    centLeg.AddEntry(g.get(), cent.label.c_str(), "pe");
    centLegendGraphs.push_back(std::move(g));
  }
  centLeg.Draw();

  TLegend varLeg(0.42, 0.78, 0.90, 0.86);
  varLeg.SetBorderSize(0);
  varLeg.SetFillStyle(0);
  varLeg.SetTextFont(42);
  varLeg.SetTextSize(0.027);
  TGraphErrors dummy(1);
  dummy.SetMarkerStyle(20);
  dummy.SetMarkerColor(kBlack);
  dummy.SetLineColor(kBlack);
  varLeg.AddEntry(&dummy, legendText.c_str(), "pe");
  varLeg.Draw();

  DrawHeader(title);

  const std::string outDir = kOutDir + "/05_ratio_plots";
  gSystem->mkdir(outDir.c_str(), true);
  const std::string outPng = outDir + "/" + outStem + ".png";
  c.SaveAs(outPng.c_str());
  Log("[DONE] Wrote " + outPng);
}

void DrawEmbeddedEfficiencyRatio(const Variant& denominator, const Variant& numerator,
                                 const std::string& ratioTitle,
                                 const std::string& numeratorAxisLabel,
                                 const std::string& outSubdir,
                                 const std::string& outStem,
                                 double yMin = 0.965,
                                 double yMax = 1.030)
{
  std::unique_ptr<TFile> fDen(TFile::Open(MergedFileName(denominator.tag).c_str(), "READ"));
  std::unique_ptr<TFile> fNum(TFile::Open(MergedFileName(numerator.tag).c_str(), "READ"));
  if (!fDen || fDen->IsZombie() || !fNum || fNum->IsZombie())
  {
    Log("[ERROR] Cannot open embedded merged files for efficiency ratio " + outStem);
    return;
  }

  std::vector<std::unique_ptr<TGraphErrors>> graphs;
  int iCent = 0;
  for (const CentBin& cent : CentBins())
  {
    const double xShift = (iCent - 1) * 0.18;
    std::unique_ptr<TGraphErrors> g(BuildRatioGraph(fNum.get(), numerator,
                                                    fDen.get(), denominator,
                                                    cent, cent.color, 21, false, xShift));
    if (g)
    {
      g->SetMarkerStyle(21);
      g->SetMarkerSize(1.45);
      g->SetLineWidth(2);
    }
    graphs.push_back(std::move(g));
    ++iCent;
  }

  TCanvas c(("c_" + outStem).c_str(), ("c_" + outStem).c_str(), 1100, 850);
  c.SetLeftMargin(0.14);
  c.SetRightMargin(0.04);
  c.SetTopMargin(0.15);
  c.SetBottomMargin(0.14);
  c.SetTicks(1, 1);

  TH1F frame(("hframe_" + outStem).c_str(), "", 100, 10.0, 35.0);
  frame.SetDirectory(nullptr);
  frame.SetStats(0);
  frame.SetMinimum(yMin);
  frame.SetMaximum(yMax);
  frame.GetXaxis()->SetTitle("p_{T}^{#gamma,truth} [GeV]");
  const std::string yTitle = "#epsilon_{#gamma}^{MC}(" + numeratorAxisLabel + ") / #epsilon_{#gamma}^{MC}(box-cuts)";
  frame.GetYaxis()->SetTitle(yTitle.c_str());
  frame.GetXaxis()->SetTitleSize(0.048);
  frame.GetYaxis()->SetTitleSize(0.038);
  frame.GetXaxis()->SetLabelSize(0.040);
  frame.GetYaxis()->SetLabelSize(0.040);
  frame.Draw();

  TLine one(10.0, 1.0, 35.0, 1.0);
  one.SetLineStyle(2);
  one.SetLineColor(kGray + 2);
  if (yMin <= 1.0 && 1.0 <= yMax) one.Draw("same");

  for (const auto& owned : graphs)
  {
    if (owned) owned->Draw("P SAME");
  }

  std::vector<std::unique_ptr<TGraphErrors>> centLegendGraphs;
  const bool ratioBelowUnity = yMax < 1.0;
  TLegend centLeg(0.18, ratioBelowUnity ? 0.18 : 0.70,
                  0.39, ratioBelowUnity ? 0.32 : 0.84);
  centLeg.SetBorderSize(0);
  centLeg.SetFillStyle(0);
  centLeg.SetTextFont(42);
  centLeg.SetTextSize(0.028);
  for (const CentBin& cent : CentBins())
  {
    std::unique_ptr<TGraphErrors> g(new TGraphErrors(1));
    g->SetMarkerStyle(21);
    g->SetMarkerSize(1.45);
    g->SetMarkerColor(cent.color);
    g->SetLineColor(cent.color);
    g->SetLineWidth(2);
    centLeg.AddEntry(g.get(), cent.label.c_str(), "pe");
    centLegendGraphs.push_back(std::move(g));
  }
  centLeg.Draw();

  TLatex tx;
  tx.SetNDC();
  tx.SetTextFont(42);
  tx.SetTextSize(0.027);
  tx.SetTextAlign(22);
  tx.DrawLatex(0.50, 0.925, ratioTitle.c_str());
  tx.SetTextAlign(31);
  tx.SetTextSize(0.030);
  tx.DrawLatex(0.92, 0.83, "#bf{sPHENIX} #it{Internal}");
  tx.DrawLatex(0.92, 0.79, "Pythia Overlay  #sqrt{s_{NN}} = 200 GeV");

  const std::string outDir = kOutDir + "/" + outSubdir;
  gSystem->mkdir(outDir.c_str(), true);
  const std::string outPng = outDir + "/" + outStem + ".png";
  c.SaveAs(outPng.c_str());
  Log("[DONE] Wrote " + outPng);
}

void DrawRatioToNoPre020(const Variant& variantA, const Variant& variantC, const Variant& variantD, const Variant& variantB)
{
  for (const Variant& v : {variantA, variantC, variantD, variantB})
  {
    if (!RequireMergedFile(v.tag)) return;
  }

  std::unique_ptr<TFile> fA(TFile::Open(MergedFileName(variantA.tag).c_str(), "READ"));
  std::unique_ptr<TFile> fC(TFile::Open(MergedFileName(variantC.tag).c_str(), "READ"));
  std::unique_ptr<TFile> fD(TFile::Open(MergedFileName(variantD.tag).c_str(), "READ"));
  std::unique_ptr<TFile> fB(TFile::Open(MergedFileName(variantB.tag).c_str(), "READ"));
  if (!fA || !fC || !fD || !fB) return;

  const CentBin cent = CentBins().front();
  std::unique_ptr<TGraphErrors> gA(BuildRatioGraph(fA.get(), variantA, fB.get(), variantB, cent, kBlack, 20, false, -0.20));
  std::unique_ptr<TGraphErrors> gC(BuildRatioGraph(fC.get(), variantC, fB.get(), variantB, cent, kBlue + 1, 20, true, 0.0));
  std::unique_ptr<TGraphErrors> gD(BuildRatioGraph(fD.get(), variantD, fB.get(), variantB, cent, kRed + 1, 21, false, 0.20));

  TCanvas c("c_effRatioVariantOverNoPre020", "c_effRatioVariantOverNoPre020", 1100, 850);
  c.SetLeftMargin(0.14);
  c.SetRightMargin(0.04);
  c.SetTopMargin(0.13);
  c.SetBottomMargin(0.14);
  c.SetTicks(1, 1);

  TH1F frame("hframe_effRatioVariantOverNoPre020", "", 100, 10.0, 35.0);
  frame.SetDirectory(nullptr);
  frame.SetStats(0);
  frame.SetMinimum(0.80);
  frame.SetMaximum(1.20);
  frame.GetXaxis()->SetTitle("p_{T}^{#gamma,truth} [GeV]");
  frame.GetYaxis()->SetTitle("#epsilon_{#gamma}^{variant} / #epsilon_{#gamma}^{no preselection}");
  frame.GetXaxis()->SetTitleSize(0.048);
  frame.GetYaxis()->SetTitleSize(0.048);
  frame.GetXaxis()->SetLabelSize(0.040);
  frame.GetYaxis()->SetLabelSize(0.040);
  frame.Draw();

  TLine one(10.0, 1.0, 35.0, 1.0);
  one.SetLineStyle(2);
  one.SetLineColor(kGray + 2);
  one.Draw("same");

  if (gA) gA->Draw("P SAME");
  if (gC) gC->Draw("P SAME");
  if (gD) gD->Draw("P SAME");

  TLegend leg(0.20, 0.70, 0.76, 0.86);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetTextFont(42);
  leg.SetTextSize(0.030);
  if (gA) leg.AddEntry(gA.get(), "New PPG12 Presel / no preselection", "pe");
  if (gC) leg.AddEntry(gC.get(), "NPB Only Presel / no preselection", "pe");
  if (gD) leg.AddEntry(gD.get(), "reference + NPB presel / no preselection", "pe");
  leg.Draw();

  DrawHeader("Photon efficiency ratios to no preselection, 0-20%, Photon+Jet Embedded Pythia (12 + 20) GeV");

  const std::string outDir = kOutDir + "/05_ratio_plots";
  gSystem->mkdir(outDir.c_str(), true);
  const std::string outPng = outDir + "/phoEff_ratio_variantA_variantC_variantD_overNoPre_cent0_20_isSliding.png";
  c.SaveAs(outPng.c_str());
  Log("[DONE] Wrote " + outPng);

  TCanvas cTwo("c_effRatioVariantAOverBVariantCOverB020", "c_effRatioVariantAOverBVariantCOverB020", 1100, 850);
  cTwo.SetLeftMargin(0.14);
  cTwo.SetRightMargin(0.04);
  cTwo.SetTopMargin(0.13);
  cTwo.SetBottomMargin(0.14);
  cTwo.SetTicks(1, 1);

  TH1F frameTwo("hframe_effRatioVariantAOverBVariantCOverB020", "", 100, 10.0, 35.0);
  frameTwo.SetDirectory(nullptr);
  frameTwo.SetStats(0);
  frameTwo.SetMinimum(0.80);
  frameTwo.SetMaximum(1.20);
  frameTwo.GetXaxis()->SetTitle("p_{T}^{#gamma,truth} [GeV]");
  frameTwo.GetYaxis()->SetTitle("#epsilon_{#gamma}^{variant} / #epsilon_{#gamma}^{no preselection}");
  frameTwo.GetXaxis()->SetTitleSize(0.048);
  frameTwo.GetYaxis()->SetTitleSize(0.048);
  frameTwo.GetXaxis()->SetLabelSize(0.040);
  frameTwo.GetYaxis()->SetLabelSize(0.040);
  frameTwo.Draw();

  TLine oneTwo(10.0, 1.0, 35.0, 1.0);
  oneTwo.SetLineStyle(2);
  oneTwo.SetLineColor(kGray + 2);
  oneTwo.Draw("same");

  if (gA) gA->Draw("P SAME");
  if (gC) gC->Draw("P SAME");

  TLegend legTwo(0.20, 0.74, 0.76, 0.86);
  legTwo.SetBorderSize(0);
  legTwo.SetFillStyle(0);
  legTwo.SetTextFont(42);
  legTwo.SetTextSize(0.030);
  if (gA) legTwo.AddEntry(gA.get(), "New PPG12 Presel / no preselection", "pe");
  if (gC) legTwo.AddEntry(gC.get(), "NPB Only Presel / no preselection", "pe");
  legTwo.Draw();

  DrawHeader("Photon efficiency ratios to no preselection, 0-20%, Photon+Jet Embedded Pythia (12 + 20) GeV");

  const std::string outPngTwo = outDir + "/phoEff_ratio_variantAOverB_variantCOverB_cent0_20_isSliding.png";
  cTwo.SaveAs(outPngTwo.c_str());
  Log("[DONE] Wrote " + outPngTwo);
}
}

void MakeEmbeddedPhotonEfficiencyCompare()
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetErrorX(0.5);

  gSystem->mkdir(kOutDir.c_str(), true);
  gSystem->mkdir((kOutDir + "/01_reference_vs_preselection_variants").c_str(), true);
  gSystem->mkdir((kOutDir + "/02_vs_no_preselection").c_str(), true);
  gSystem->mkdir((kOutDir + "/03_preselection_variant_overlays").c_str(), true);
  gSystem->mkdir((kOutDir + "/04_reference_plus_npb_id_overlays").c_str(), true);
  gSystem->mkdir((kOutDir + "/04_reference_preselection_id_overlays").c_str(), true);
  gSystem->mkdir((kOutDir + "/04_variantA_preselection_id_overlays").c_str(), true);
  gSystem->mkdir((kOutDir + "/05_ratio_plots").c_str(), true);
  gSystem->mkdir((kOutDir + "/06_text_summaries").c_str(), true);
  gSystem->mkdir(kPPOutDir.c_str(), true);

  if (gSummary.is_open()) gSummary.close();
  gSummary.open((kOutDir + "/06_text_summaries/summary_photonEfficiencyCompare_fromMergedFiles.txt").c_str());
  Log("[PHOTON EFFICIENCY DEBUG] epsilon_gamma^MC = 1 - N_truth_miss / N_truth_signal");
  Log("  merged base: " + kMergedBaseDir);
  Log("  plot base  : " + kOutDir);
  Log("  pp plot dir: " + kPPOutDir);

  const std::string base = "jetMinPt5_7pi_8_vz60_isoR40_isSliding_baseVariant_";
  const std::string tail = "_tightReference_nonTightReference";
  const Variant reference{base + "preselectionReference" + tail, "sliding_reference"};
  const Variant variantA{base + "preselectionVariantA" + tail, "sliding_variantA"};
  const Variant variantB{base + "preselectionVariantB" + tail, "sliding_variantB"};
  const Variant variantC{base + "preselectionVariantC" + tail, "sliding_variantC"};
  const Variant variantD{base + "preselectionVariantD" + tail, "sliding_variantD"};

  DrawPPFixedIsoEfficiencyPair();

  DrawEfficiencyPair(
      "Photon efficiency, reference vs New PPG12 Presel, Photon+Jet Embedded Pythia (12 + 20) GeV",
      "phoEff_referenceVsNPBPreselection_isSliding",
      reference,
      variantA,
      kReferencePreselectionLabel,
      kVariantAPreselectionLabel,
      "01_reference_vs_preselection_variants");

  DrawEfficiencyPair(
      "Photon efficiency, reference vs no preselection, Photon+Jet Embedded Pythia (12 + 20) GeV",
      "phoEff_referenceVsNoPreselection_isSliding",
      reference,
      variantB,
      kReferencePreselectionLabel,
      kNoPreselectionLabel,
      "01_reference_vs_preselection_variants");

  DrawEfficiencyPair(
      "Photon efficiency, reference vs NPB Only Presel, Photon+Jet Embedded Pythia (12 + 20) GeV",
      "phoEff_referenceVsNPBOnlyPresel_isSliding",
      reference,
      variantC,
      kReferencePreselectionLabel,
      kVariantCPreselectionLabel,
      "01_reference_vs_preselection_variants");

  DrawEfficiencyPair(
      "Photon efficiency, reference vs reference + NPB presel, Photon+Jet Embedded Pythia (12 + 20) GeV",
      "phoEff_referenceVsReferencePlusNPBPresel_isSliding",
      reference,
      variantD,
      kReferencePreselectionLabel,
      kVariantDPreselectionLabel,
      "01_reference_vs_preselection_variants");

  DrawEfficiencyPair(
      "MC #gamma-efficiency, original PPG12 box-cuts with and without NPB presel, Emb Photon+Jet (12 + 20)",
      "phoEff_embOriginalBoxCuts_vs_originalBoxCutsPlusNPBPresel_isSliding",
      reference,
      variantD,
      "Emb MC original PPG12 box-cuts",
      "Emb MC original PPG12 box-cuts + NPB presel",
      "04_reference_plus_npb_id_overlays",
      true);

  DrawEmbeddedEfficiencyRatio(
      reference,
      variantD,
      "MC #gamma-efficiency ratio, original PPG12 box-cuts + NPB presel over original box-cuts",
      "box-cuts + NPB presel",
      "04_reference_plus_npb_id_overlays",
      "phoEff_ratio_embOriginalBoxCutsPlusNPBPresel_overOriginalBoxCuts_isSliding");

  DrawEfficiencyPair(
      "Photon efficiency, New PPG12 Presel vs no preselection, Photon+Jet Embedded Pythia (12 + 20) GeV",
      "phoEff_NPBVsNoPreselection_isSliding",
      variantA,
      variantB,
      kVariantAPreselectionLabel,
      kNoPreselectionLabel,
      "02_vs_no_preselection");

  DrawEfficiencyPair(
      "Photon efficiency, NPB Only Presel vs no preselection, Photon+Jet Embedded Pythia (12 + 20) GeV",
      "phoEff_NPBOnlyPreselVsNoPreselection_isSliding",
      variantC,
      variantB,
      kVariantCPreselectionLabel,
      kNoPreselectionLabel,
      "02_vs_no_preselection");

  DrawEfficiencyPair(
      "Photon efficiency, New PPG12 Presel vs NPB Only Presel, Photon+Jet Embedded Pythia (12 + 20) GeV",
      "phoEff_NewPPG12PreselVsNPBOnlyPresel_isSliding",
      variantA,
      variantC,
      kVariantAPreselectionLabel,
      kVariantCPreselectionLabel,
      "03_preselection_variant_overlays");

  const std::string baseD = "jetMinPt5_7pi_8_vz60_isoR40_isSliding_baseVariant_preselectionVariantD_";
  const Variant refId{baseD + "tightReference_nonTightReference", "reference + NPB presel, Reference tight/nontight"};
  const Variant bdtLooseNt{baseD + "tightVariantA_nonTightReference", "reference + NPB presel, BDT tight/loose nontight"};
  const Variant bdtSidebandNt{baseD + "tightVariantA_nonTightVariantA", "reference + NPB presel, BDT tight/BDT nontight"};
  const std::string idOverlayTitle =
      "MC #gamma-efficiency, box-cuts versus BDT, Photon+Jet (5 + 10 + 20) and Emb Photon+Jet (12 + 20)";

  DrawEfficiencyPair(
      idOverlayTitle,
      "phoEff_refPlusNPB_referenceID_vs_BDTTightLooseNT_isSliding",
      refId,
      bdtLooseNt,
      "Emb MC original PPG12 box-cuts",
      "Emb MC updated PPG12 BDT tight/nontight",
      "04_reference_plus_npb_id_overlays");

  DrawEfficiencyPair(
      idOverlayTitle,
      "phoEff_refPlusNPB_referenceID_vs_BDTTightBDTNT_isSliding",
      refId,
      bdtSidebandNt,
      "Emb MC original PPG12 box-cuts",
      "Emb MC updated PPG12 BDT tight/nontight",
      "04_reference_plus_npb_id_overlays");

  const std::string baseRef = "jetMinPt5_7pi_8_vz60_isoR40_isSliding_baseVariant_preselectionReference_";
  const Variant refPreRefId{baseRef + "tightReference_nonTightReference", "reference presel, Reference tight/nontight"};
  const Variant refPreBdtLooseNt{baseRef + "tightVariantA_nonTightReference", "reference presel, BDT tight/loose nontight"};
  const Variant refPreBdtSidebandNt{baseRef + "tightVariantA_nonTightVariantA", "reference presel, BDT tight/BDT nontight"};

  DrawEfficiencyPair(
      "MC #gamma-efficiency, original PPG12 box-cuts versus ref tight ID + BDT, Emb Photon+Jet (12 + 20)",
      "phoEff_embOriginalBoxCuts_vs_refPreselBDTTightLooseNT_isSliding",
      refPreRefId,
      refPreBdtLooseNt,
      "Emb MC original PPG12 box-cuts",
      "Emb MC ref tight ID + BDT",
      "04_reference_preselection_id_overlays",
      true);

  DrawEmbeddedEfficiencyRatio(
      refPreRefId,
      refPreBdtLooseNt,
      "MC #gamma-efficiency ratio, ref tight ID + BDT over original PPG12 box-cuts",
      "ref tight ID + BDT",
      "04_reference_preselection_id_overlays",
      "phoEff_ratio_refPreselBDTTightLooseNT_overOriginalBoxCuts_isSliding",
      0.0,
      0.85);

  DrawEfficiencyPair(
      idOverlayTitle,
      "phoEff_referencePresel_referenceID_vs_BDTTightLooseNT_isSliding",
      refPreRefId,
      refPreBdtLooseNt,
      "Emb MC original PPG12 box-cuts",
      "Emb MC updated PPG12 BDT tight/nontight",
      "04_reference_preselection_id_overlays");

  DrawEfficiencyPair(
      idOverlayTitle,
      "phoEff_referencePresel_referenceID_vs_BDTTightBDTNT_isSliding",
      refPreRefId,
      refPreBdtSidebandNt,
      "Emb MC original PPG12 box-cuts",
      "Emb MC updated PPG12 BDT tight/nontight",
      "04_reference_preselection_id_overlays");

  const std::string baseVarA = "jetMinPt5_7pi_8_vz60_isoR40_isSliding_baseVariant_preselectionVariantA_";
  const Variant varAPreRefId{baseVarA + "tightReference_nonTightReference", "variantA presel, Reference tight/nontight"};
  const Variant varAPreBdtLooseNt{baseVarA + "tightVariantA_nonTightReference", "variantA presel, BDT tight/loose nontight"};
  const Variant varAPreBdtSidebandNt{baseVarA + "tightVariantA_nonTightVariantA", "variantA presel, BDT tight/BDT nontight"};

  DrawEfficiencyPair(
      idOverlayTitle,
      "phoEff_variantAPresel_referenceID_vs_BDTTightLooseNT_isSliding",
      varAPreRefId,
      varAPreBdtLooseNt,
      "Emb MC original PPG12 box-cuts",
      "Emb MC updated PPG12 BDT tight/nontight",
      "04_variantA_preselection_id_overlays");

  DrawEfficiencyPair(
      idOverlayTitle,
      "phoEff_variantAPresel_referenceID_vs_BDTTightBDTNT_isSliding",
      varAPreRefId,
      varAPreBdtSidebandNt,
      "Emb MC original PPG12 box-cuts",
      "Emb MC updated PPG12 BDT tight/nontight",
      "04_variantA_preselection_id_overlays");

  DrawRatioSingle(
      "Photon efficiency ratio, reference / NPB Only Presel, Photon+Jet Embedded Pythia (12 + 20) GeV",
      "phoEff_ratio_referenceOverNPBOnlyPresel_isSliding",
      reference,
      variantC,
      "reference preselection / NPB Only Presel");

  DrawRatioSingle(
      "Photon efficiency ratio, reference / New PPG12 Presel, Photon+Jet Embedded Pythia (12 + 20) GeV",
      "phoEff_ratio_referenceOverNPBPreselection_isSliding",
      reference,
      variantA,
      "reference preselection / New PPG12 Presel");

  DrawRatioToNoPre020(variantA, variantC, variantD, variantB);

  if (gSummary.is_open()) gSummary.close();
}
