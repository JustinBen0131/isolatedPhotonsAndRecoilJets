#include <TCanvas.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TKey.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TNamed.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TH1.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace
{
constexpr double kSigmaEmbeddedPhoton12To20_pb = 2598.12425;
constexpr double kSigmaEmbeddedPhoton20Plus_pb = 133.317866;

const std::string kInputDir = "InputFiles/simEmbedded";
const std::string kMergedDir = "InputFiles/simEmbedded/merged";
const std::string kOutDir = "dataOutput/combinedSimOnlyEMBEDDED/b-leakageCompare";
const std::string kCOutDir = "dataOutput/combinedSimOnlyEMBEDDED/c-leakageCompare";
const std::string kTopDir = "SIM";
const std::string kMergedFilename = "RecoilJets_embeddedPhoton12plus20_MERGED.root";

struct PtBin
{
  double lo = 0.0;
  double hi = 0.0;
  std::string suffix;
};

struct Variant
{
  std::string tag;
  std::string label;
  int color = kBlack;
  int marker = 20;
};

struct CentBin
{
  std::string label;
  std::string suffix;
  int color = kBlack;
};

struct PointRecord
{
  PtBin bin;
  std::string centLabel;
  double A = 0.0;
  double B = 0.0;
  double C = 0.0;
  double D = 0.0;
  double eA = 0.0;
  double eB = 0.0;
  double eC = 0.0;
  double fB = 0.0;
  double fC = 0.0;
  double fD = 0.0;
};

struct VariantRecord
{
  std::string label;
  std::string tag;
  std::string mergedPath;
  std::vector<PointRecord> points;
  double sumA = 0.0;
  double sumB = 0.0;
  double sumC = 0.0;
  double sumD = 0.0;
};

std::vector<VariantRecord> gRecords;
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

std::string InputFileName(const std::string& sample, const std::string& tag)
{
  return kInputDir + "/RecoilJets_" + sample + "_ALL_" + tag + ".root";
}

std::string MergedOutDir(const std::string& tag)
{
  return kMergedDir + "/" + tag + "/photonJet12and20merged_SIM";
}

std::string MergedFileName(const std::string& tag)
{
  return MergedOutDir(tag) + "/" + kMergedFilename;
}

double ReadEventCount(TFile* f)
{
  if (!f) return 0.0;
  TDirectory* d = f->GetDirectory(kTopDir.c_str());
  if (!d) return 0.0;
  TH1* cnt = dynamic_cast<TH1*>(d->Get("cnt_SIM"));
  if (!cnt) return 0.0;
  return cnt->GetBinContent(1);
}

bool IsDirClass(const std::string& cls)
{
  return cls == "TDirectoryFile" || cls == "TDirectory";
}

void AddScaledRecursive(TDirectory* outDir, TDirectory* inDir, double weight)
{
  if (!outDir || !inDir) return;

  TIter next(inDir->GetListOfKeys());
  while (TKey* key = dynamic_cast<TKey*>(next()))
  {
    const std::string name = key->GetName();
    const std::string cls = key->GetClassName();

    if (IsDirClass(cls))
    {
      TDirectory* subIn = dynamic_cast<TDirectory*>(inDir->Get(name.c_str()));
      if (!subIn) continue;
      outDir->cd();
      TDirectory* subOut = outDir->GetDirectory(name.c_str());
      if (!subOut) subOut = outDir->mkdir(name.c_str());
      AddScaledRecursive(subOut, subIn, weight);
      continue;
    }

    TObject* objIn = inDir->Get(name.c_str());
    if (!objIn) continue;

    if (TH1* hIn = dynamic_cast<TH1*>(objIn))
    {
      outDir->cd();
      TH1* hOut = dynamic_cast<TH1*>(outDir->Get(name.c_str()));
      if (!hOut)
      {
        TH1* hNew = dynamic_cast<TH1*>(hIn->Clone(name.c_str()));
        if (!hNew) continue;
        hNew->SetDirectory(outDir);
        if (hNew->GetSumw2N() == 0) hNew->Sumw2();
        hNew->Scale(weight);
        hNew->Write(name.c_str(), TObject::kOverwrite);
      }
      else
      {
        std::unique_ptr<TH1> tmp(dynamic_cast<TH1*>(hIn->Clone((name + "_tmpAdd").c_str())));
        if (!tmp) continue;
        tmp->SetDirectory(nullptr);
        if (tmp->GetSumw2N() == 0) tmp->Sumw2();
        tmp->Scale(weight);
        hOut->Add(tmp.get());
        hOut->Write(name.c_str(), TObject::kOverwrite);
      }
      continue;
    }

    outDir->cd();
    if (!outDir->Get(name.c_str()))
    {
      TObject* clone = objIn->Clone(name.c_str());
      if (clone) clone->Write(name.c_str(), TObject::kOverwrite);
    }
  }
}

bool BuildMergedFile(const std::string& tag)
{
  const std::string file12 = InputFileName("embeddedPhoton12", tag);
  const std::string file20 = InputFileName("embeddedPhoton20", tag);
  const std::string outPath = MergedFileName(tag);

  std::unique_ptr<TFile> f12(TFile::Open(file12.c_str(), "READ"));
  std::unique_ptr<TFile> f20(TFile::Open(file20.c_str(), "READ"));
  if (!f12 || f12->IsZombie() || !f20 || f20->IsZombie())
  {
    Log("[ERROR] Cannot open both embedded inputs for " + tag);
    Log("        " + file12);
    Log("        " + file20);
    return false;
  }

  TDirectory* d12 = f12->GetDirectory(kTopDir.c_str());
  TDirectory* d20 = f20->GetDirectory(kTopDir.c_str());
  if (!d12 || !d20)
  {
    Log("[ERROR] Missing SIM directory while merging " + tag);
    return false;
  }

  const double n12 = ReadEventCount(f12.get());
  const double n20 = ReadEventCount(f20.get());
  if (n12 <= 0.0 || n20 <= 0.0)
  {
    std::ostringstream os;
    os << "[ERROR] Bad raw event counters for " << tag << ": N12=" << n12 << " N20=" << n20;
    Log(os.str());
    return false;
  }

  const double perEvent12 = kSigmaEmbeddedPhoton12To20_pb / n12;
  const double perEvent20 = kSigmaEmbeddedPhoton20Plus_pb / n20;
  const double w12 = perEvent12 / perEvent20;
  const double w20 = 1.0;

  gSystem->mkdir(MergedOutDir(tag).c_str(), true);
  std::unique_ptr<TFile> fout(TFile::Open(outPath.c_str(), "RECREATE"));
  if (!fout || fout->IsZombie())
  {
    Log("[ERROR] Cannot create merged file: " + outPath);
    return false;
  }

  fout->cd();
  TDirectory* outTop = fout->mkdir(kTopDir.c_str());
  if (!outTop)
  {
    Log("[ERROR] Cannot create SIM directory in merged file: " + outPath);
    return false;
  }

  AddScaledRecursive(outTop, d12, w12);
  AddScaledRecursive(outTop, d20, w20);

  TTree* tOutED = nullptr;
  for (TFile* fin : {f12.get(), f20.get()})
  {
    TTree* tInED = dynamic_cast<TTree*>(fin->Get("EventDisplayTree"));
    if (!tInED)
    {
      TDirectory* din = fin->GetDirectory(kTopDir.c_str());
      if (din) tInED = dynamic_cast<TTree*>(din->Get("EventDisplayTree"));
    }
    if (!tInED) continue;
    if (!tOutED)
    {
      outTop->cd();
      tOutED = tInED->CloneTree(0);
      if (tOutED) tOutED->SetDirectory(outTop);
    }
    if (tOutED) tOutED->CopyEntries(tInED);
  }
  if (tOutED)
  {
    outTop->cd();
    tOutED->Write("EventDisplayTree", TObject::kOverwrite);
  }

  std::ostringstream meta;
  meta << "Merged embedded PhotonJet12+20. N12=" << std::fixed << std::setprecision(0) << n12
       << " N20=" << n20
       << " sigma12_pb=" << std::setprecision(12) << kSigmaEmbeddedPhoton12To20_pb
       << " sigma20_pb=" << kSigmaEmbeddedPhoton20Plus_pb
       << " w12_relative=" << w12
       << " w20_relative=" << w20
       << " convention=(sigma/Nraw)/(sigma20/Nraw20)";
  outTop->cd();
  TNamed mergeInfo("MERGE_INFO", meta.str().c_str());
  mergeInfo.Write("MERGE_INFO", TObject::kOverwrite);

  fout->Write("", TObject::kOverwrite);
  fout->Close();

  std::ostringstream os;
  os << "[MERGE] " << tag << "\n"
     << "        out = " << outPath << "\n"
     << "        N12=" << std::fixed << std::setprecision(0) << n12
     << " N20=" << n20
     << " w12=" << std::setprecision(8) << w12
     << " w20=" << w20;
  Log(os.str());
  return true;
}

bool ParsePtSuffix(const std::string& name, PtBin& out)
{
  const std::string prefix = "h_sigABCD_MC_pT_";
  if (name.rfind(prefix, 0) != 0) return false;
  if (name.find("_cent_") != std::string::npos) return false;

  const std::string rest = name.substr(prefix.size());
  const std::size_t sep = rest.find('_');
  if (sep == std::string::npos) return false;

  char* end = nullptr;
  const double lo = std::strtod(rest.substr(0, sep).c_str(), &end);
  if (!end || *end != '\0') return false;
  const double hi = std::strtod(rest.substr(sep + 1).c_str(), &end);
  if (!end || *end != '\0') return false;
  if (!(hi > lo)) return false;

  out.lo = lo;
  out.hi = hi;
  out.suffix = "_pT_" + rest;
  return true;
}

std::vector<PtBin> DiscoverPtBins(TFile* f)
{
  std::vector<PtBin> bins;
  if (!f) return bins;
  TDirectory* d = f->GetDirectory(kTopDir.c_str());
  if (!d) return bins;

  std::set<std::string> seen;
  TIter it(d->GetListOfKeys());
  while (TKey* key = dynamic_cast<TKey*>(it()))
  {
    PtBin b;
    if (!ParsePtSuffix(key->GetName(), b)) continue;
    if (seen.insert(b.suffix).second) bins.push_back(b);
  }

  std::sort(bins.begin(), bins.end(), [](const PtBin& a, const PtBin& b) {
    if (a.lo != b.lo) return a.lo < b.lo;
    return a.hi < b.hi;
  });
  return bins;
}

bool ReadABCD(TFile* f, const PtBin& b, const CentBin& cent, PointRecord& p)
{
  if (!f) return false;
  TDirectory* d = f->GetDirectory(kTopDir.c_str());
  if (!d) return false;

  const std::string name = "h_sigABCD_MC" + b.suffix + cent.suffix;
  TH1* h = dynamic_cast<TH1*>(d->Get(name.c_str()));
  if (!h || h->GetNbinsX() < 4) return false;

  p.bin = b;
  p.centLabel = cent.label;
  p.A = h->GetBinContent(1);
  p.B = h->GetBinContent(2);
  p.C = h->GetBinContent(3);
  p.D = h->GetBinContent(4);
  p.eA = h->GetBinError(1);
  p.eB = h->GetBinError(2);
  p.eC = h->GetBinError(3);
  p.fB = (p.A > 0.0) ? p.B / p.A : 0.0;
  p.fC = (p.A > 0.0) ? p.C / p.A : 0.0;
  p.fD = (p.A > 0.0) ? p.D / p.A : 0.0;
  return p.A > 0.0;
}

TGraphErrors* BuildBLeakageGraphFromMerged(TFile* f, const Variant& variant, const CentBin& cent, bool openMarker)
{
  std::vector<PtBin> bins = DiscoverPtBins(f);
  if (bins.empty())
  {
    Log("[WARN] No merged h_sigABCD_MC pT bins found for " + variant.tag);
    return nullptr;
  }

  VariantRecord rec;
  rec.label = variant.label + " " + cent.label;
  rec.tag = variant.tag;
  rec.mergedPath = MergedFileName(variant.tag);

  std::vector<double> x, ex, y, ey;
  Log("");
  Log("[LEAKAGE FROM MERGED FILE] " + variant.label + "  centrality " + cent.label);
  Log("  file: " + rec.mergedPath);
  Log("  pTbin        A_sig        B_sig        C_sig        D_sig          fB          fC          fD");
  Log("  ----------------------------------------------------------------------------------------------");

  for (const PtBin& pb : bins)
  {
    if (pb.lo < 10.0) continue;

    PointRecord p;
    if (!ReadABCD(f, pb, cent, p)) continue;

    rec.points.push_back(p);
    rec.sumA += p.A;
    rec.sumB += p.B;
    rec.sumC += p.C;
    rec.sumD += p.D;

    std::ostringstream line;
    line << "  " << std::setw(5) << std::right << std::fixed << std::setprecision(0) << pb.lo
         << "-" << std::setw(2) << pb.hi
         << std::setw(13) << std::setprecision(3) << p.A
         << std::setw(13) << p.B
         << std::setw(13) << p.C
         << std::setw(13) << p.D
         << std::setw(12) << std::setprecision(6) << p.fB
         << std::setw(12) << p.fC
         << std::setw(12) << p.fD;
    Log(line.str());

    x.push_back(0.5 * (pb.lo + pb.hi));
    ex.push_back(0.5 * (pb.hi - pb.lo));
    y.push_back(p.fB);
    double err = 0.0;
    if (p.A > 0.0 && p.B > 0.0)
    {
      err = p.fB * std::sqrt((p.eB / p.B) * (p.eB / p.B) + (p.eA / p.A) * (p.eA / p.A));
    }
    ey.push_back(err);
  }

  if (rec.sumA > 0.0)
  {
    std::ostringstream line;
    line << "  overall: fB=" << std::fixed << std::setprecision(6) << rec.sumB / rec.sumA
         << " fC=" << rec.sumC / rec.sumA
         << " fD=" << rec.sumD / rec.sumA
         << " binsWithA=" << rec.points.size();
    Log(line.str());
  }
  gRecords.push_back(rec);

  if (x.empty()) return nullptr;
  TGraphErrors* g = new TGraphErrors(static_cast<int>(x.size()), x.data(), y.data(), ex.data(), ey.data());
  g->SetName(("gBLeak_" + variant.label + "_" + cent.label).c_str());
  g->SetLineColor(cent.color);
  g->SetMarkerColor(cent.color);
  g->SetMarkerStyle(openMarker ? 24 : 20);
  g->SetMarkerSize(1.1);
  g->SetLineWidth(2);
  g->SetLineStyle(openMarker ? 2 : 1);
  return g;
}

TGraphErrors* BuildCLeakageGraphFromMerged(TFile* f, const Variant& variant, const CentBin& cent, bool openMarker)
{
  std::vector<PtBin> bins = DiscoverPtBins(f);
  if (bins.empty())
  {
    Log("[WARN] No merged h_sigABCD_MC pT bins found for " + variant.tag);
    return nullptr;
  }

  VariantRecord rec;
  rec.label = variant.label + " " + cent.label;
  rec.tag = variant.tag;
  rec.mergedPath = MergedFileName(variant.tag);

  std::vector<double> x, ex, y, ey;
  Log("");
  Log("[C LEAKAGE FROM MERGED FILE] " + variant.label + "  centrality " + cent.label);
  Log("  file: " + rec.mergedPath);
  Log("  pTbin        A_sig        B_sig        C_sig        D_sig          fB          fC          fD");
  Log("  ----------------------------------------------------------------------------------------------");

  for (const PtBin& pb : bins)
  {
    if (pb.lo < 10.0) continue;

    PointRecord p;
    if (!ReadABCD(f, pb, cent, p)) continue;

    rec.points.push_back(p);
    rec.sumA += p.A;
    rec.sumB += p.B;
    rec.sumC += p.C;
    rec.sumD += p.D;

    std::ostringstream line;
    line << "  " << std::setw(5) << std::right << std::fixed << std::setprecision(0) << pb.lo
         << "-" << std::setw(2) << pb.hi
         << std::setw(13) << std::setprecision(3) << p.A
         << std::setw(13) << p.B
         << std::setw(13) << p.C
         << std::setw(13) << p.D
         << std::setw(12) << std::setprecision(6) << p.fB
         << std::setw(12) << p.fC
         << std::setw(12) << p.fD;
    Log(line.str());

    x.push_back(0.5 * (pb.lo + pb.hi));
    ex.push_back(0.5 * (pb.hi - pb.lo));
    y.push_back(p.fC);
    double err = 0.0;
    if (p.A > 0.0 && p.C > 0.0)
    {
      err = p.fC * std::sqrt((p.eC / p.C) * (p.eC / p.C) + (p.eA / p.A) * (p.eA / p.A));
    }
    ey.push_back(err);
  }

  if (rec.sumA > 0.0)
  {
    std::ostringstream line;
    line << "  overall: fB=" << std::fixed << std::setprecision(6) << rec.sumB / rec.sumA
         << " fC=" << rec.sumC / rec.sumA
         << " fD=" << rec.sumD / rec.sumA
         << " binsWithA=" << rec.points.size();
    Log(line.str());
  }
  gRecords.push_back(rec);

  if (x.empty()) return nullptr;
  TGraphErrors* g = new TGraphErrors(static_cast<int>(x.size()), x.data(), y.data(), ex.data(), ey.data());
  g->SetName(("gCLeak_" + variant.label + "_" + cent.label).c_str());
  g->SetLineColor(cent.color);
  g->SetMarkerColor(cent.color);
  g->SetMarkerStyle(openMarker ? 24 : 20);
  g->SetMarkerSize(1.1);
  g->SetLineWidth(2);
  g->SetLineStyle(openMarker ? 2 : 1);
  return g;
}

void DrawOne(const std::string& title, const std::string& outStem,
             const Variant& fixed, const Variant& sliding)
{
  if (!BuildMergedFile(fixed.tag)) return;
  if (!BuildMergedFile(sliding.tag)) return;

  std::unique_ptr<TFile> fFixed(TFile::Open(MergedFileName(fixed.tag).c_str(), "READ"));
  std::unique_ptr<TFile> fSliding(TFile::Open(MergedFileName(sliding.tag).c_str(), "READ"));
  if (!fFixed || fFixed->IsZombie() || !fSliding || fSliding->IsZombie())
  {
    Log("[ERROR] Cannot open one or both merged files for " + outStem);
    return;
  }

  std::vector<std::unique_ptr<TGraphErrors>> graphs;
  graphs.reserve(2 * CentBins().size());
  for (const CentBin& cent : CentBins())
  {
    graphs.emplace_back(BuildBLeakageGraphFromMerged(fFixed.get(), fixed, cent, false));
    graphs.emplace_back(BuildBLeakageGraphFromMerged(fSliding.get(), sliding, cent, true));
  }

  double ymax = 0.0;
  for (const auto& owned : graphs)
  {
    TGraphErrors* g = owned.get();
    if (!g) continue;
    for (int i = 0; i < g->GetN(); ++i)
    {
      double gx = 0.0, gy = 0.0;
      g->GetPoint(i, gx, gy);
      ymax = std::max(ymax, gy + g->GetErrorY(i));
    }
  }
  if (ymax <= 0.0) ymax = 1.0;

  TCanvas c(("c_" + outStem).c_str(), ("c_" + outStem).c_str(), 1100, 850);
  c.SetLeftMargin(0.14);
  c.SetRightMargin(0.04);
  c.SetTopMargin(0.13);
  c.SetBottomMargin(0.14);
  c.SetTicks(1, 1);

  TH1F frame(("hframe_" + outStem).c_str(), "", 100, 10.0, 35.0);
  frame.SetDirectory(nullptr);
  frame.SetStats(0);
  frame.SetMinimum(0.0);
  frame.SetMaximum(std::max(1.00, 1.55 * ymax));
  frame.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
  frame.GetYaxis()->SetTitle("Region B signal leakage, f_{B} = B_{sig}/A_{sig}");
  frame.GetXaxis()->SetTitleSize(0.048);
  frame.GetYaxis()->SetTitleSize(0.048);
  frame.GetXaxis()->SetLabelSize(0.040);
  frame.GetYaxis()->SetLabelSize(0.040);
  frame.Draw();

  for (const auto& owned : graphs)
  {
    if (owned) owned->Draw("P SAME");
  }

  TGraphErrors legFixed(1);
  legFixed.SetMarkerStyle(20);
  legFixed.SetMarkerColor(kBlack);
  legFixed.SetLineColor(kBlack);
  legFixed.SetLineStyle(1);
  legFixed.SetLineWidth(2);
  TGraphErrors legSliding(1);
  legSliding.SetMarkerStyle(24);
  legSliding.SetMarkerColor(kBlack);
  legSliding.SetLineColor(kBlack);
  legSliding.SetLineStyle(2);
  legSliding.SetLineWidth(2);

  TLegend isoLeg(0.44, 0.72, 0.92, 0.86);
  isoLeg.SetBorderSize(0);
  isoLeg.SetFillStyle(0);
  isoLeg.SetTextFont(42);
  isoLeg.SetTextSize(0.030);
  isoLeg.AddEntry(&legFixed, "Fixed E_{T}^{iso} < 2 GeV", "pe");
  isoLeg.AddEntry(&legSliding, "E_{T}^{iso}(cent) < 7.57 - 0.0673[cent %]", "pe");
  isoLeg.Draw();

  std::vector<std::unique_ptr<TGraphErrors>> centLegendGraphs;
  centLegendGraphs.reserve(CentBins().size());
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

  TLatex tx;
  tx.SetNDC();
  tx.SetTextFont(42);
  tx.SetTextAlign(22);
  tx.SetTextSize(0.030);
  tx.DrawLatex(0.50, 0.890, title.c_str());

  tx.SetTextAlign(31);
  tx.SetTextSize(0.032);
  tx.DrawLatex(0.92, 0.675, "#bf{sPHENIX} #it{Internal}");
  tx.DrawLatex(0.92, 0.635, "Pythia Overlay  #sqrt{s_{NN}} = 200 GeV");

  const std::string outPng = kOutDir + "/" + outStem + ".png";
  c.SaveAs(outPng.c_str());
  Log("[DONE] Wrote " + outPng);
}

void DrawCPreselectionPair(const std::string& title,
                           const std::string& outStem,
                           const Variant& first,
                           const Variant& second,
                           const std::string& firstLegend,
                           const std::string& secondLegend)
{
  if (!BuildMergedFile(first.tag)) return;
  if (!BuildMergedFile(second.tag)) return;

  std::unique_ptr<TFile> fFirst(TFile::Open(MergedFileName(first.tag).c_str(), "READ"));
  std::unique_ptr<TFile> fSecond(TFile::Open(MergedFileName(second.tag).c_str(), "READ"));
  if (!fFirst || fFirst->IsZombie() || !fSecond || fSecond->IsZombie())
  {
    Log("[ERROR] Cannot open one or both merged files for " + outStem);
    return;
  }

  std::vector<std::unique_ptr<TGraphErrors>> graphs;
  graphs.reserve(2 * CentBins().size());
  for (const CentBin& cent : CentBins())
  {
    graphs.emplace_back(BuildCLeakageGraphFromMerged(fFirst.get(), first, cent, false));
    graphs.emplace_back(BuildCLeakageGraphFromMerged(fSecond.get(), second, cent, true));
  }

  double ymax = 0.0;
  for (const auto& owned : graphs)
  {
    TGraphErrors* g = owned.get();
    if (!g) continue;
    for (int i = 0; i < g->GetN(); ++i)
    {
      double gx = 0.0, gy = 0.0;
      g->GetPoint(i, gx, gy);
      ymax = std::max(ymax, gy + g->GetErrorY(i));
    }
  }
  if (ymax <= 0.0) ymax = 1.0;

  TCanvas c(("c_" + outStem).c_str(), ("c_" + outStem).c_str(), 1100, 850);
  c.SetLeftMargin(0.14);
  c.SetRightMargin(0.04);
  c.SetTopMargin(0.13);
  c.SetBottomMargin(0.14);
  c.SetTicks(1, 1);

  TH1F frame(("hframe_" + outStem).c_str(), "", 100, 10.0, 35.0);
  frame.SetDirectory(nullptr);
  frame.SetStats(0);
  frame.SetMinimum(0.0);
  frame.SetMaximum(std::max(4.0, 1.65 * ymax));
  frame.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
  frame.GetYaxis()->SetTitle("Region C signal leakage, f_{C} = C_{sig}/A_{sig}");
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
  legFirst.SetLineStyle(1);
  legFirst.SetLineWidth(2);
  TGraphErrors legSecond(1);
  legSecond.SetMarkerStyle(24);
  legSecond.SetMarkerColor(kBlack);
  legSecond.SetLineColor(kBlack);
  legSecond.SetLineStyle(2);
  legSecond.SetLineWidth(2);

  TLegend varLeg(0.44, 0.72, 0.92, 0.86);
  varLeg.SetBorderSize(0);
  varLeg.SetFillStyle(0);
  varLeg.SetTextFont(42);
  varLeg.SetTextSize(0.030);
  varLeg.AddEntry(&legFirst, firstLegend.c_str(), "pe");
  varLeg.AddEntry(&legSecond, secondLegend.c_str(), "pe");
  varLeg.Draw();

  std::vector<std::unique_ptr<TGraphErrors>> centLegendGraphs;
  centLegendGraphs.reserve(CentBins().size());
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

  TLatex tx;
  tx.SetNDC();
  tx.SetTextFont(42);
  tx.SetTextAlign(22);
  tx.SetTextSize(0.030);
  tx.DrawLatex(0.50, 0.890, title.c_str());

  tx.SetTextAlign(31);
  tx.SetTextSize(0.032);
  tx.DrawLatex(0.92, 0.675, "#bf{sPHENIX} #it{Internal}");
  tx.DrawLatex(0.92, 0.635, "Pythia Overlay  #sqrt{s_{NN}} = 200 GeV");

  const std::string outPng = kCOutDir + "/" + outStem + ".png";
  c.SaveAs(outPng.c_str());
  Log("[DONE] Wrote " + outPng);
}

const VariantRecord* FindRecord(const std::string& label)
{
  for (const VariantRecord& r : gRecords)
  {
    if (r.label == label) return &r;
  }
  return nullptr;
}

void PrintComparisonSummary()
{
  Log("");
  Log("[OVERALL C-LEAKAGE RANKING FROM MERGED FILES]");
  std::vector<const VariantRecord*> ranked;
  for (const VariantRecord& r : gRecords)
  {
    if (r.sumA > 0.0) ranked.push_back(&r);
  }
  std::sort(ranked.begin(), ranked.end(), [](const VariantRecord* a, const VariantRecord* b) {
    return (a->sumC / a->sumA) < (b->sumC / b->sumA);
  });
  for (std::size_t i = 0; i < ranked.size(); ++i)
  {
    const VariantRecord* r = ranked[i];
    std::ostringstream os;
    os << "  " << (i + 1) << ". " << std::setw(17) << std::left << r->label
       << " fC=" << std::fixed << std::setprecision(6) << r->sumC / r->sumA
       << " fB=" << r->sumB / r->sumA
       << " bins=" << r->points.size();
    Log(os.str());
  }

  auto compare = [&](const std::string& a, const std::string& b)
  {
    const VariantRecord* ra = FindRecord(a);
    const VariantRecord* rb = FindRecord(b);
    if (!ra || !rb) return;
    double maxDB = 0.0;
    double maxDC = 0.0;
    for (const PointRecord& pa : ra->points)
    {
      for (const PointRecord& pb : rb->points)
      {
        if (pa.bin.suffix != pb.bin.suffix) continue;
        maxDB = std::max(maxDB, std::abs(pa.fB - pb.fB));
        maxDC = std::max(maxDC, std::abs(pa.fC - pb.fC));
      }
    }
    std::ostringstream os;
    os << "  " << a << " vs " << b
       << ": max |delta fB|=" << std::fixed << std::setprecision(8) << maxDB
       << " max |delta fC|=" << maxDC;
    Log(os.str());
  };

  Log("");
  Log("[PAIRWISE DIFFERENCE CHECK]");
  for (const CentBin& cent : CentBins())
  {
    Log("  centrality " + cent.label);
    compare("fixed_reference " + cent.label, "fixed_variantA " + cent.label);
    compare("fixed_reference " + cent.label, "fixed_variantB " + cent.label);
    compare("fixed_variantA " + cent.label, "fixed_variantB " + cent.label);
    compare("sliding_reference " + cent.label, "sliding_variantA " + cent.label);
    compare("sliding_reference " + cent.label, "sliding_variantB " + cent.label);
    compare("sliding_variantA " + cent.label, "sliding_variantB " + cent.label);
  }
}

void PrintCPreselectionSummary()
{
  Log("");
  Log("[OVERALL C-LEAKAGE RANKING FROM MERGED SLIDING FILES]");
  std::vector<const VariantRecord*> ranked;
  for (const VariantRecord& r : gRecords)
  {
    if (r.sumA > 0.0) ranked.push_back(&r);
  }
  std::sort(ranked.begin(), ranked.end(), [](const VariantRecord* a, const VariantRecord* b) {
    return (a->sumC / a->sumA) < (b->sumC / b->sumA);
  });
  for (std::size_t i = 0; i < ranked.size(); ++i)
  {
    const VariantRecord* r = ranked[i];
    std::ostringstream os;
    os << "  " << (i + 1) << ". " << std::setw(24) << std::left << r->label
       << " fC=" << std::fixed << std::setprecision(6) << r->sumC / r->sumA
       << " fB=" << r->sumB / r->sumA
       << " bins=" << r->points.size();
    Log(os.str());
  }
}

void MakeEmbeddedCLeakageCompare_SlidingOnly()
{
  gRecords.clear();
  gSystem->mkdir(kCOutDir.c_str(), true);
  if (gSummary.is_open()) gSummary.close();
  gSummary.open((kCOutDir + "/summary_cLeakageSlidingPreselection_fromMergedFiles.txt").c_str());
  Log("[C-LEAKAGE DEBUG] Building canonical merged ROOT files, then reading centrality-binned C leakage from sliding-iso merged files.");
  Log("  merged base: " + kMergedDir);
  Log("  plot base  : " + kCOutDir);

  const std::string base = "jetMinPt5_7pi_8_vz60_isoR40_isSliding_baseVariant_";
  const std::string tail = "_tightReference_nonTightReference";
  const Variant reference{base + "preselectionReference" + tail, "sliding_reference", kBlack, 20};
  const Variant variantA{base + "preselectionVariantA" + tail, "sliding_variantA", kBlack, 20};
  const Variant variantB{base + "preselectionVariantB" + tail, "sliding_variantB", kBlack, 20};

  DrawCPreselectionPair(
      "C-leakage, reference vs NPB preselection, Photon+Jet Embedded Pythia (12 + 20) GeV",
      "cLeakage_referenceVsNPBPreselection_isSliding",
      reference,
      variantA,
      "reference preselection",
      "NPB preselection");

  DrawCPreselectionPair(
      "C-leakage, reference vs no preselection, Photon+Jet Embedded Pythia (12 + 20) GeV",
      "cLeakage_referenceVsNoPreselection_isSliding",
      reference,
      variantB,
      "reference preselection",
      "no preselection");

  DrawCPreselectionPair(
      "C-leakage, NPB vs no preselection, Photon+Jet Embedded Pythia (12 + 20) GeV",
      "cLeakage_NPBVsNoPreselection_isSliding",
      variantA,
      variantB,
      "NPB preselection",
      "no preselection");

  PrintCPreselectionSummary();
  if (gSummary.is_open()) gSummary.close();
}
}

void MakeEmbeddedBLeakageCompare()
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);
  gStyle->SetErrorX(0.5);
  gSystem->mkdir(kOutDir.c_str(), true);
  gSystem->mkdir(kMergedDir.c_str(), true);

  gRecords.clear();
  gSummary.open((kOutDir + "/summary_bLeakage_fromMergedFiles.txt").c_str());
  Log("[B-LEAKAGE DEBUG] Building canonical merged ROOT files, then reading leakage from those merged files.");
  Log("  merged base: " + kMergedDir);
  Log("  plot base  : " + kOutDir);

  const std::string base = "jetMinPt5_7pi_8_vz60_isoR40_";
  const std::string mid = "_baseVariant_";
  const std::string tail = "_tightReference_nonTightReference";

  DrawOne(
      "B-leakage, reference preselection, Photon+Jet Embedded Pythia (12 + 20) GeV",
      "bLeakage_referencePreselection_fixedVsSliding",
      {base + "fixedIso2GeV" + mid + "preselectionReference" + tail, "fixed_reference", kBlack, 20},
      {base + "isSliding" + mid + "preselectionReference" + tail, "sliding_reference", kBlue + 1, 24});

  DrawOne(
      "B-leakage, NPB preselection, Photon+Jet Embedded Pythia (12 + 20) GeV",
      "bLeakage_NPBPreselection_fixedVsSliding",
      {base + "fixedIso2GeV" + mid + "preselectionVariantA" + tail, "fixed_variantA", kBlack, 20},
      {base + "isSliding" + mid + "preselectionVariantA" + tail, "sliding_variantA", kBlue + 1, 24});

  DrawOne(
      "B-leakage, no preselection, Photon+Jet Embedded Pythia (12 + 20) GeV",
      "bLeakage_noPreselection_fixedVsSliding",
      {base + "fixedIso2GeV" + mid + "preselectionVariantB" + tail, "fixed_variantB", kBlack, 20},
      {base + "isSliding" + mid + "preselectionVariantB" + tail, "sliding_variantB", kBlue + 1, 24});

  PrintComparisonSummary();
  if (gSummary.is_open()) gSummary.close();

  MakeEmbeddedCLeakageCompare_SlidingOnly();
}
