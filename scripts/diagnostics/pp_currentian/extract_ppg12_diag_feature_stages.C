#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>
#include <TKey.h>
#include <TROOT.h>
#include <TString.h>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace
{
TH1* findHistRecursive(TDirectory* dir, const TString& name)
{
  if (!dir) return nullptr;

  if (auto* obj = dir->Get(name))
  {
    if (obj->InheritsFrom(TH1::Class())) return static_cast<TH1*>(obj);
  }

  TIter next(dir->GetListOfKeys());
  while (auto* keyObj = next())
  {
    auto* key = static_cast<TKey*>(keyObj);
    if (!key) continue;
    TObject* obj = key->ReadObj();
    if (!obj) continue;
    if (obj->InheritsFrom(TDirectory::Class()))
    {
      if (auto* found = findHistRecursive(static_cast<TDirectory*>(obj), name))
      {
        return found;
      }
    }
  }
  return nullptr;
}

double binValue(TH1* h, int ibin)
{
  return h ? h->GetBinContent(ibin) : 0.0;
}

double safeDiv(double num, double den)
{
  return den != 0.0 ? num / den : std::numeric_limits<double>::quiet_NaN();
}

void printFeatureSummary(TFile* file, const char* histName)
{
  TH1* h = findHistRecursive(file, histName);
  if (!h)
  {
    std::cout << "feature," << histName << ",missing,0,nan,nan\n";
    return;
  }
  std::cout << "feature," << histName
            << ",present," << std::setprecision(12) << h->Integral()
            << "," << h->GetMean()
            << "," << h->GetRMS()
            << "\n";
}
}

void extract_ppg12_diag_feature_stages(const char* rootPath)
{
  TFile* file = TFile::Open(rootPath, "READ");
  if (!file || file->IsZombie())
  {
    std::cerr << "ERROR cannot open " << (rootPath ? rootPath : "(null)") << "\n";
    return;
  }

  const std::vector<std::pair<std::string, std::string>> histNames = {
      {"raw_cluster", "h_ppg12_diag_raw_signal_cluster_0"},
      {"raw_edge_accept", "h_ppg12_diag_raw_signal_edge_accept_0"},
      {"raw_unmasked", "h_ppg12_diag_raw_signal_unmasked_0"},
      {"raw_response_window", "h_ppg12_diag_raw_signal_response_window_0"},
      {"selected_signal", "h_ppg12_diag_signal_cluster_0"},
      {"preselection_fail", "h_ppg12_diag_preselection_fail_signal_0"},
      {"common", "h_ppg12_diag_common_signal_0"},
      {"bdt_tight", "h_ppg12_diag_bdt_tight_window_signal_0"},
      {"bdt_nontight", "h_ppg12_diag_bdt_nontight_window_signal_0"},
      {"nontight_shape", "h_ppg12_diag_nontight_shape_signal_0"},
      {"tag_tight", "h_ppg12_diag_tag_tight_signal_0"},
      {"tag_nontight", "h_ppg12_diag_tag_nontight_signal_0"},
      {"tag_neither", "h_ppg12_diag_tag_neither_signal_0"},
      {"all", "h_all_cluster_signal_0"},
      {"tight", "h_tight_cluster_signal_0"},
      {"A", "h_tight_iso_cluster_signal_0"},
      {"B", "h_tight_noniso_cluster_signal_0"},
      {"C", "h_nontight_iso_cluster_signal_0"},
      {"D", "h_nontight_noniso_cluster_signal_0"},
  };

  std::map<std::string, TH1*> hists;
  for (const auto& item : histNames)
  {
    hists[item.first] = findHistRecursive(file, item.second.c_str());
  }

  TH1* axisHist = hists["common"];
  if (!axisHist) axisHist = hists["A"];
  if (!axisHist)
  {
    std::cerr << "ERROR missing both common and A histograms in " << rootPath << "\n";
    file->Close();
    return;
  }

  std::cout << std::setprecision(12);
  std::cout << "source," << rootPath << "\n";
  std::cout << "hist_status";
  for (const auto& item : histNames)
  {
    std::cout << "," << item.first << ":" << (hists[item.first] ? "present" : "missing");
  }
  std::cout << "\n";

  std::cout
      << "row,pt_lo,pt_hi"
      << ",raw_cluster,raw_edge_accept,raw_unmasked,raw_response_window"
      << ",selected_signal,preselection_fail,common"
      << ",bdt_tight,bdt_nontight,nontight_shape"
      << ",tag_tight,tag_nontight,tag_neither"
      << ",all,tight,A,B,C,D"
      << ",selected_over_raw_unmasked,common_over_selected,preselection_fail_over_selected"
      << ",bdt_tight_over_common,bdt_nontight_over_common,nontight_shape_over_common"
      << ",tag_tight_over_common,tag_nontight_over_common,tag_neither_over_common"
      << ",A_over_common,B_over_common,C_over_common,D_over_common,nt_over_common,D_share_nt"
      << "\n";

  double sumCommon = 0.0;
  double sumSelected = 0.0;
  double sumTagNontight = 0.0;
  double sumBDTNontight = 0.0;
  double sumC = 0.0;
  double sumD = 0.0;
  double sumA = 0.0;
  double sumB = 0.0;

  for (int ibin = 1; ibin <= axisHist->GetNbinsX(); ++ibin)
  {
    const double lo = axisHist->GetXaxis()->GetBinLowEdge(ibin);
    const double hi = axisHist->GetXaxis()->GetBinUpEdge(ibin);
    const double rawCluster = binValue(hists["raw_cluster"], ibin);
    const double rawEdgeAccept = binValue(hists["raw_edge_accept"], ibin);
    const double rawUnmasked = binValue(hists["raw_unmasked"], ibin);
    const double rawResponseWindow = binValue(hists["raw_response_window"], ibin);
    const double selectedSignal = binValue(hists["selected_signal"], ibin);
    const double preselectionFail = binValue(hists["preselection_fail"], ibin);
    const double common = binValue(hists["common"], ibin);
    const double bdtTight = binValue(hists["bdt_tight"], ibin);
    const double bdtNontight = binValue(hists["bdt_nontight"], ibin);
    const double nontightShape = binValue(hists["nontight_shape"], ibin);
    const double tagTight = binValue(hists["tag_tight"], ibin);
    const double tagNontight = binValue(hists["tag_nontight"], ibin);
    const double tagNeither = binValue(hists["tag_neither"], ibin);
    const double all = binValue(hists["all"], ibin);
    const double tight = binValue(hists["tight"], ibin);
    const double A = binValue(hists["A"], ibin);
    const double B = binValue(hists["B"], ibin);
    const double C = binValue(hists["C"], ibin);
    const double D = binValue(hists["D"], ibin);
    const double nt = C + D;

    if (hi <= 26.0)
    {
      sumCommon += common;
      sumSelected += selectedSignal;
      sumTagNontight += tagNontight;
      sumBDTNontight += bdtNontight;
      sumA += A;
      sumB += B;
      sumC += C;
      sumD += D;
    }

    std::cout << "bin," << lo << "," << hi
              << "," << rawCluster
              << "," << rawEdgeAccept
              << "," << rawUnmasked
              << "," << rawResponseWindow
              << "," << selectedSignal
              << "," << preselectionFail
              << "," << common
              << "," << bdtTight
              << "," << bdtNontight
              << "," << nontightShape
              << "," << tagTight
              << "," << tagNontight
              << "," << tagNeither
              << "," << all
              << "," << tight
              << "," << A
              << "," << B
              << "," << C
              << "," << D
              << "," << safeDiv(selectedSignal, rawUnmasked)
              << "," << safeDiv(common, selectedSignal)
              << "," << safeDiv(preselectionFail, selectedSignal)
              << "," << safeDiv(bdtTight, common)
              << "," << safeDiv(bdtNontight, common)
              << "," << safeDiv(nontightShape, common)
              << "," << safeDiv(tagTight, common)
              << "," << safeDiv(tagNontight, common)
              << "," << safeDiv(tagNeither, common)
              << "," << safeDiv(A, common)
              << "," << safeDiv(B, common)
              << "," << safeDiv(C, common)
              << "," << safeDiv(D, common)
              << "," << safeDiv(nt, common)
              << "," << safeDiv(D, nt)
              << "\n";
  }

  std::cout << "stable_summary,pt_hi_le_26"
            << ",selected_over_common=" << safeDiv(sumSelected, sumCommon)
            << ",tag_nontight_over_common=" << safeDiv(sumTagNontight, sumCommon)
            << ",bdt_nontight_over_common=" << safeDiv(sumBDTNontight, sumCommon)
            << ",B_over_A=" << safeDiv(sumB, sumA)
            << ",C_over_A=" << safeDiv(sumC, sumA)
            << ",D_over_A=" << safeDiv(sumD, sumA)
            << ",nt_over_A=" << safeDiv(sumC + sumD, sumA)
            << ",D_share_nt=" << safeDiv(sumD, sumC + sumD)
            << "\n";

  printFeatureSummary(file, "h_ppg12_diag_common_feature_bdt_score_0");
  printFeatureSummary(file, "h_ppg12_diag_common_tag_tight_feature_bdt_score_0");
  printFeatureSummary(file, "h_ppg12_diag_common_tag_nontight_feature_bdt_score_0");
  printFeatureSummary(file, "h_ppg12_diag_common_tag_neither_feature_bdt_score_0");
  printFeatureSummary(file, "h_ppg12_diag_common_feature_eiso_0");
  printFeatureSummary(file, "h_ppg12_diag_common_tag_nontight_feature_eiso_0");

  file->Close();
}
