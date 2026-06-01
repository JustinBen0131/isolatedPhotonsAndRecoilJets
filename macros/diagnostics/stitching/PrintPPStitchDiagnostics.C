#include <TDirectory.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace
{
struct Sample
{
  std::string label;
  std::string path;
  double sigmaPb = 1.0;
  TFile* file = nullptr;
  TDirectory* dir = nullptr;
  double nRaw = 0.0;
  double relWeight = 1.0;
};

struct PtBin
{
  std::string label;
  std::string suffix;
};

const std::vector<PtBin>& PtBins()
{
  static const std::vector<PtBin> bins = {
      {"8-10",  "_pT_8_10"},
      {"10-12", "_pT_10_12"},
      {"12-14", "_pT_12_14"},
      {"14-16", "_pT_14_16"},
      {"16-18", "_pT_16_18"},
      {"18-20", "_pT_18_20"},
      {"20-22", "_pT_20_22"},
      {"22-24", "_pT_22_24"},
      {"24-26", "_pT_24_26"},
      {"26-35", "_pT_26_35"},
  };
  return bins;
}

const std::vector<std::string>& Variants()
{
  static const std::vector<std::string> variants = {
      "g4Stored",
      "hepmcParent",
      "hepmcParentStable",
      "hepmcAnyPhoton",
      "g4NoEmbed",
  };
  return variants;
}

TH1* GetH1(TDirectory* dir, const std::string& name)
{
  return dir ? dynamic_cast<TH1*>(dir->Get(name.c_str())) : nullptr;
}

TH2* GetH2(TDirectory* dir, const std::string& name)
{
  return dir ? dynamic_cast<TH2*>(dir->Get(name.c_str())) : nullptr;
}

double Integral(TH1* h)
{
  return h ? h->Integral(0, h->GetNbinsX() + 1) : 0.0;
}

double BinContentByPtBinIndex(TH1* h, const std::size_t ptBinIndex)
{
  if (!h) return 0.0;
  const int bin = static_cast<int>(ptBinIndex) + 1;
  if (bin < 1 || bin > h->GetNbinsX()) return 0.0;
  return h->GetBinContent(bin);
}

double SumPtBinRange(TH1* h, const std::size_t firstBin, const std::size_t lastBinInclusive)
{
  if (!h) return 0.0;
  double sum = 0.0;
  for (std::size_t i = firstBin; i <= lastBinInclusive && i < PtBins().size(); ++i)
  {
    sum += BinContentByPtBinIndex(h, i);
  }
  return sum;
}

double CountFromCntSim(TDirectory* dir)
{
  TH1* h = GetH1(dir, "cnt_SIM");
  const double n = Integral(h);
  return (std::isfinite(n) && n > 0.0) ? n : 0.0;
}

bool Quantile(TH1* h, const double prob, double& q)
{
  q = std::numeric_limits<double>::quiet_NaN();
  if (!h || Integral(h) <= 0.0) return false;
  double p[1] = {prob};
  double x[1] = {0.0};
  h->GetQuantiles(1, x, p);
  if (!std::isfinite(x[0])) return false;
  q = x[0];
  return true;
}

std::string Fmt(const double x, const int prec = 3)
{
  if (!std::isfinite(x)) return "nan";
  std::ostringstream os;
  os << std::fixed << std::setprecision(prec) << x;
  return os.str();
}

std::unique_ptr<TH1> CloneScaled(TH1* h, const std::string& name, const double w)
{
  if (!h) return nullptr;
  std::unique_ptr<TH1> out(static_cast<TH1*>(h->Clone(name.c_str())));
  if (!out) return nullptr;
  out->SetDirectory(nullptr);
  out->Scale(w);
  return out;
}

std::unique_ptr<TH1> PseudoMerge(const std::vector<Sample>& samples,
                                 const std::string& histName,
                                 std::vector<double>& contributions)
{
  contributions.assign(samples.size(), 0.0);
  std::unique_ptr<TH1> merged;

  for (std::size_t i = 0; i < samples.size(); ++i)
  {
    TH1* h = GetH1(samples[i].dir, histName);
    if (!h) continue;

    contributions[i] = Integral(h) * samples[i].relWeight;
    auto scaled = CloneScaled(h, histName + "_scaled_" + samples[i].label, samples[i].relWeight);
    if (!scaled) continue;

    if (!merged)
    {
      merged = std::move(scaled);
      merged->SetName((histName + "_pseudoMerged").c_str());
    }
    else
    {
      merged->Add(scaled.get());
    }
  }

  return merged;
}

void PrintHeader(const std::string& s)
{
  std::cout << "\n============================================================\n"
            << s << "\n"
            << "============================================================\n";
}

void PrintIndependentTable(const std::vector<Sample>& samples)
{
  PrintHeader("Independent sample q80 from h_EisoReco_truthSigMatched");
  std::cout << std::left << std::setw(9) << "pTbin";
  for (const auto& s : samples)
  {
    std::cout << std::right << std::setw(15) << (s.label + " q80")
              << std::setw(15) << (s.label + " int");
  }
  std::cout << "\n";

  for (const auto& b : PtBins())
  {
    std::cout << std::left << std::setw(9) << b.label;
    for (const auto& s : samples)
    {
      TH1* h = GetH1(s.dir, "h_EisoReco_truthSigMatched" + b.suffix);
      double q80 = std::numeric_limits<double>::quiet_NaN();
      Quantile(h, 0.80, q80);
      std::cout << std::right << std::setw(15) << Fmt(q80)
                << std::setw(15) << Fmt(Integral(h), 1);
    }
    std::cout << "\n";
  }
}

void PrintFlowTable(const std::vector<Sample>& samples, const double targetRaw)
{
  PrintHeader("Truth-signal match flow counters, raw no-vtx-weight entries");

  const std::vector<std::pair<std::string, std::string>> flowHists = {
      {"truthIso truth pT", "h_ppStitchDiag_flow_truthIso_truthPt_noVtxW"},
      {"truthIso+G4 truth pT", "h_ppStitchDiag_flow_truthIsoG4_truthPt_noVtxW"},
      {"raw cluster dR<0.05 truth pT", "h_ppStitchDiag_flow_rawClusterDR005_truthPt_noVtxW"},
      {"raw cluster dR<0.10 truth pT", "h_ppStitchDiag_flow_rawClusterDR010_truthPt_noVtxW"},
      {"raw PPG12 track truth pT", "h_ppStitchDiag_flow_rawPPG12Track_truthPt_noVtxW"},
      {"raw PPG12 track fid truth pT", "h_ppStitchDiag_flow_rawPPG12TrackFid_truthPt_noVtxW"},
      {"selected photon dR<0.05 truth pT", "h_ppStitchDiag_flow_selectedPhotonDR005_truthPt_noVtxW"},
      {"selected photon dR<0.10 truth pT", "h_ppStitchDiag_flow_selectedPhotonDR010_truthPt_noVtxW"},
      {"selected PPG12 track truth pT", "h_ppStitchDiag_flow_selectedPPG12Track_truthPt_noVtxW"},
      {"selected PPG12 track fid truth pT", "h_ppStitchDiag_flow_selectedPPG12TrackFid_truthPt_noVtxW"},
      {"legacy barcode+dR<0.05 truth pT", "h_ppStitchDiag_flow_legacyBarcodeDR005_truthPt_noVtxW"},
      {"current PPG12 reco pT", "h_ppStitchDiag_flow_recoMatched_recoPt_noVtxW"},
      {"validEiso reco pT", "h_ppStitchDiag_flow_validEiso_recoPt_noVtxW"},
      {"hepmcParent kept validEiso reco pT", "h_ppStitchDiag_flow_stitch_hepmcParent_validEiso_recoPt_noVtxW"},
  };

  std::cout << std::left << std::setw(39) << "counter"
            << std::setw(10) << "sample"
            << std::right << std::setw(12) << "total"
            << std::setw(12) << "14-16"
            << std::setw(12) << "16-18"
            << std::setw(12) << "18-20"
            << std::setw(12) << "20-22"
            << std::setw(12) << "14-20"
            << "\n";

  for (const auto& fh : flowHists)
  {
    for (const auto& s : samples)
    {
      TH1* h = GetH1(s.dir, fh.second);
      const double b1416 = BinContentByPtBinIndex(h, 3);
      const double b1618 = BinContentByPtBinIndex(h, 4);
      const double b1820 = BinContentByPtBinIndex(h, 5);
      const double b2022 = BinContentByPtBinIndex(h, 6);
      const double sum1420 = b1416 + b1618 + b1820;

      std::cout << std::left << std::setw(39) << fh.first
                << std::setw(10) << s.label
                << std::right << std::setw(12) << Fmt(Integral(h), 1)
                << std::setw(12) << Fmt(b1416, 1)
                << std::setw(12) << Fmt(b1618, 1)
                << std::setw(12) << Fmt(b1820, 1)
                << std::setw(12) << Fmt(b2022, 1)
                << std::setw(12) << Fmt(sum1420, 1)
                << "\n";
    }
  }

  PrintHeader("Adaptive localSTITCHtest decision inputs");

  const auto& pj10 = samples.at(1);
  const auto& pj20 = samples.at(2);
  TH1* pj10StitchKept = GetH1(pj10.dir, "h_ppPhotonStitch_hepmcParent_maxPhotonPt_kept");
  TH1* pj10TruthG4 = GetH1(pj10.dir, "h_ppStitchDiag_flow_truthIsoG4_truthPt_noVtxW");
  TH1* pj10RawDR005 = GetH1(pj10.dir, "h_ppStitchDiag_flow_rawClusterDR005_truthPt_noVtxW");
  TH1* pj10RawDR010 = GetH1(pj10.dir, "h_ppStitchDiag_flow_rawClusterDR010_truthPt_noVtxW");
  TH1* pj10RawPPG12 = GetH1(pj10.dir, "h_ppStitchDiag_flow_rawPPG12Track_truthPt_noVtxW");
  TH1* pj10RawPPG12Fid = GetH1(pj10.dir, "h_ppStitchDiag_flow_rawPPG12TrackFid_truthPt_noVtxW");
  TH1* pj10SelDR005 = GetH1(pj10.dir, "h_ppStitchDiag_flow_selectedPhotonDR005_truthPt_noVtxW");
  TH1* pj10SelDR010 = GetH1(pj10.dir, "h_ppStitchDiag_flow_selectedPhotonDR010_truthPt_noVtxW");
  TH1* pj10SelPPG12 = GetH1(pj10.dir, "h_ppStitchDiag_flow_selectedPPG12Track_truthPt_noVtxW");
  TH1* pj10SelPPG12Fid = GetH1(pj10.dir, "h_ppStitchDiag_flow_selectedPPG12TrackFid_truthPt_noVtxW");
  TH1* pj10Legacy = GetH1(pj10.dir, "h_ppStitchDiag_flow_legacyBarcodeDR005_truthPt_noVtxW");
  TH1* pj10ValidReco = GetH1(pj10.dir, "h_ppStitchDiag_flow_validEiso_recoPt_noVtxW");
  TH1* pj10ParentKeptReco = GetH1(pj10.dir, "h_ppStitchDiag_flow_stitch_hepmcParent_validEiso_recoPt_noVtxW");
  TH1* pj20ParentKeptReco = GetH1(pj20.dir, "h_ppStitchDiag_flow_stitch_hepmcParent_validEiso_recoPt_noVtxW");

  const double pj10KeptEvents = Integral(pj10StitchKept);
  const double pj10TruthG4_14_22 = SumPtBinRange(pj10TruthG4, 3, 6);
  const double pj10RawDR005_14_22 = SumPtBinRange(pj10RawDR005, 3, 6);
  const double pj10RawDR010_14_22 = SumPtBinRange(pj10RawDR010, 3, 6);
  const double pj10RawPPG12_14_22 = SumPtBinRange(pj10RawPPG12, 3, 6);
  const double pj10RawPPG12Fid_14_22 = SumPtBinRange(pj10RawPPG12Fid, 3, 6);
  const double pj10SelDR005_14_22 = SumPtBinRange(pj10SelDR005, 3, 6);
  const double pj10SelDR010_14_22 = SumPtBinRange(pj10SelDR010, 3, 6);
  const double pj10SelPPG12_14_22 = SumPtBinRange(pj10SelPPG12, 3, 6);
  const double pj10SelPPG12Fid_14_22 = SumPtBinRange(pj10SelPPG12Fid, 3, 6);
  const double pj10Legacy_14_22 = SumPtBinRange(pj10Legacy, 3, 6);
  const double pj10ValidReco_14_20 = SumPtBinRange(pj10ValidReco, 3, 5);
  const double pj10ParentReco_14_20 = SumPtBinRange(pj10ParentKeptReco, 3, 5);
  const double pj20ParentReco_14_20 = SumPtBinRange(pj20ParentKeptReco, 3, 5);
  const double totalParentReco_14_20 = pj10ParentReco_14_20 + pj20ParentReco_14_20;

  std::cout << "LOCALSTITCHTEST_TARGET_RAW=" << Fmt(targetRaw, 1) << "\n";
  std::cout << "LOCALSTITCHTEST_METRIC pj10_hepmcParent_kept_events=" << Fmt(pj10KeptEvents, 1) << "\n";
  std::cout << "LOCALSTITCHTEST_METRIC pj10_truthIsoG4_truthPt_14_22=" << Fmt(pj10TruthG4_14_22, 1) << "\n";
  std::cout << "LOCALSTITCHTEST_METRIC pj10_rawClusterDR005_truthPt_14_22=" << Fmt(pj10RawDR005_14_22, 1) << "\n";
  std::cout << "LOCALSTITCHTEST_METRIC pj10_rawClusterDR010_truthPt_14_22=" << Fmt(pj10RawDR010_14_22, 1) << "\n";
  std::cout << "LOCALSTITCHTEST_METRIC pj10_rawPPG12Track_truthPt_14_22=" << Fmt(pj10RawPPG12_14_22, 1) << "\n";
  std::cout << "LOCALSTITCHTEST_METRIC pj10_rawPPG12TrackFid_truthPt_14_22=" << Fmt(pj10RawPPG12Fid_14_22, 1) << "\n";
  std::cout << "LOCALSTITCHTEST_METRIC pj10_selectedPhotonDR005_truthPt_14_22=" << Fmt(pj10SelDR005_14_22, 1) << "\n";
  std::cout << "LOCALSTITCHTEST_METRIC pj10_selectedPhotonDR010_truthPt_14_22=" << Fmt(pj10SelDR010_14_22, 1) << "\n";
  std::cout << "LOCALSTITCHTEST_METRIC pj10_selectedPPG12Track_truthPt_14_22=" << Fmt(pj10SelPPG12_14_22, 1) << "\n";
  std::cout << "LOCALSTITCHTEST_METRIC pj10_selectedPPG12TrackFid_truthPt_14_22=" << Fmt(pj10SelPPG12Fid_14_22, 1) << "\n";
  std::cout << "LOCALSTITCHTEST_METRIC pj10_legacyBarcodeDR005_truthPt_14_22=" << Fmt(pj10Legacy_14_22, 1) << "\n";
  std::cout << "LOCALSTITCHTEST_METRIC pj10_validEiso_recoPt_14_20=" << Fmt(pj10ValidReco_14_20, 1) << "\n";
  std::cout << "LOCALSTITCHTEST_METRIC pj10_hepmcParent_kept_validEiso_recoPt_14_20=" << Fmt(pj10ParentReco_14_20, 1) << "\n";
  std::cout << "LOCALSTITCHTEST_METRIC pj20_hepmcParent_kept_validEiso_recoPt_14_20=" << Fmt(pj20ParentReco_14_20, 1) << "\n";
  std::cout << "LOCALSTITCHTEST_METRIC total_hepmcParent_kept_validEiso_recoPt_14_20=" << Fmt(totalParentReco_14_20, 1) << "\n";

  const double minSignal = std::max(3.0, targetRaw / 5.0);
  const bool enoughStitchOwnership = (pj10KeptEvents >= targetRaw);
  const bool enoughPj10Truth = (pj10TruthG4_14_22 >= targetRaw);
  const bool enoughTransitionReco = (totalParentReco_14_20 >= targetRaw);
  const bool decisivePj10Branch =
      enoughPj10Truth &&
      (pj10RawDR010_14_22 < minSignal ||
       pj10RawPPG12Fid_14_22 < minSignal ||
       pj10SelDR010_14_22 < minSignal ||
       pj10SelPPG12Fid_14_22 < minSignal ||
       pj10ValidReco_14_20 < minSignal ||
       pj10ParentReco_14_20 < minSignal ||
       enoughTransitionReco);

  if (decisivePj10Branch || (enoughTransitionReco && enoughStitchOwnership))
  {
    std::cout << "LOCALSTITCHTEST_STATUS=PASS\n";
  }
  else
  {
    std::cout << "LOCALSTITCHTEST_STATUS=NEED_MORE_STATS\n";
  }

  if (enoughStitchOwnership && pj10TruthG4_14_22 < minSignal)
  {
    std::cout << "LOCALSTITCHTEST_DIAGNOSIS=pj10 has stitch-owned events, but very few PPG12 truthIso+G4 signal photons in 14-22 GeV; investigate truth-signal definition vs stitch variable.\n";
  }
  else if (enoughPj10Truth && pj10RawDR010_14_22 < minSignal)
  {
    std::cout << "LOCALSTITCHTEST_DIAGNOSIS=pj10 has PPG12 truthIso+G4 signal photons, but no nearby raw CEMC cluster; investigate sample content/input pairing or truth photon acceptance vs reco cluster formation.\n";
  }
  else if (enoughPj10Truth && pj10RawDR010_14_22 >= minSignal && pj10RawPPG12Fid_14_22 < minSignal)
  {
    std::cout << "LOCALSTITCHTEST_DIAGNOSIS=pj10 has nearby raw clusters, but few fiducial raw clusters whose max truth primary is the PPG12 signal photon; investigate cluster truth ownership, not stitching weights.\n";
  }
  else if (enoughPj10Truth && pj10RawPPG12Fid_14_22 >= minSignal && pj10SelDR010_14_22 < minSignal)
  {
    std::cout << "LOCALSTITCHTEST_DIAGNOSIS=pj10 raw PPG12 signal clusters exist, but PhotonClusterBuilder does not produce selected photons near them; investigate photon-builder cuts/preselection.\n";
  }
  else if (enoughPj10Truth && pj10SelDR010_14_22 >= minSignal && pj10SelPPG12Fid_14_22 < minSignal)
  {
    std::cout << "LOCALSTITCHTEST_DIAGNOSIS=pj10 selected photons exist near signal truth photons, but their PPG12 max-truth track ownership differs; compare PPG12 track match vs legacy barcode+dR.\n";
  }
  else if (enoughPj10Truth && pj10SelPPG12Fid_14_22 >= minSignal && pj10ValidReco_14_20 < minSignal)
  {
    std::cout << "LOCALSTITCHTEST_DIAGNOSIS=pj10 has selected PPG12 truth-track photons, but very few valid-Eiso reco photons in 14-20; investigate isolation calculation/bin migration.\n";
  }
  else if (pj10ValidReco_14_20 >= minSignal && pj10ParentReco_14_20 < minSignal)
  {
    std::cout << "LOCALSTITCHTEST_DIAGNOSIS=pj10 has valid reco truth-signal photons, but stitch ownership removes them; investigate event ownership window/variable.\n";
  }
  else
  {
    std::cout << "LOCALSTITCHTEST_DIAGNOSIS=transition-region statistics are sufficient or still accumulating; inspect flow table above.\n";
  }
}

void PrintStitchSummary(const std::vector<Sample>& samples)
{
  PrintHeader("Stitch variable QA by sample");
  std::cout << std::left << std::setw(19) << "variant"
            << std::setw(10) << "sample"
            << std::right << std::setw(12) << "all"
            << std::setw(12) << "kept"
            << std::setw(12) << "rejected"
            << std::setw(12) << "keepFrac"
            << std::setw(12) << "meanKept"
            << "\n";

  for (const auto& v : Variants())
  {
    for (const auto& s : samples)
    {
      TH1* hAll = GetH1(s.dir, "h_ppPhotonStitch_" + v + "_maxPhotonPt_all");
      TH1* hKept = GetH1(s.dir, "h_ppPhotonStitch_" + v + "_maxPhotonPt_kept");
      TH1* hRejected = GetH1(s.dir, "h_ppPhotonStitch_" + v + "_maxPhotonPt_rejected");
      const double nAll = Integral(hAll);
      const double nKept = Integral(hKept);
      const double nRej = Integral(hRejected);
      const double keepFrac = (nAll > 0.0) ? nKept / nAll : std::numeric_limits<double>::quiet_NaN();
      const double meanKept = hKept ? hKept->GetMean() : std::numeric_limits<double>::quiet_NaN();

      std::cout << std::left << std::setw(19) << v
                << std::setw(10) << s.label
                << std::right << std::setw(12) << Fmt(nAll, 1)
                << std::setw(12) << Fmt(nKept, 1)
                << std::setw(12) << Fmt(nRej, 1)
                << std::setw(12) << Fmt(keepFrac, 3)
                << std::setw(12) << Fmt(meanKept, 3)
                << "\n";
    }
  }
}

void PrintVariantMergeTable(const std::vector<Sample>& samples,
                            const std::string& variant,
                            const std::string& weightMode)
{
  const std::string base = "h_EisoReco_truthSigMatched_stitch_" + variant + "_" + weightMode;

  PrintHeader("Pseudo-merged q70/q80/q90 and sample fractions: " + variant + " " + weightMode);
  std::cout << std::left << std::setw(9) << "pTbin"
            << std::right << std::setw(10) << "q70"
            << std::setw(10) << "q80"
            << std::setw(10) << "q90"
            << std::setw(12) << "int"
            << std::setw(10) << "pj5"
            << std::setw(10) << "pj10"
            << std::setw(10) << "pj20"
            << "\n";

  double roughness = 0.0;
  double prevQ80 = std::numeric_limits<double>::quiet_NaN();
  int nSteps = 0;

  for (const auto& b : PtBins())
  {
    std::vector<double> contrib;
    std::unique_ptr<TH1> h = PseudoMerge(samples, base + b.suffix, contrib);

    double q70 = std::numeric_limits<double>::quiet_NaN();
    double q80 = std::numeric_limits<double>::quiet_NaN();
    double q90 = std::numeric_limits<double>::quiet_NaN();
    Quantile(h.get(), 0.70, q70);
    Quantile(h.get(), 0.80, q80);
    Quantile(h.get(), 0.90, q90);

    const double total = h ? Integral(h.get()) : 0.0;
    std::cout << std::left << std::setw(9) << b.label
              << std::right << std::setw(10) << Fmt(q70)
              << std::setw(10) << Fmt(q80)
              << std::setw(10) << Fmt(q90)
              << std::setw(12) << Fmt(total, 1);

    for (std::size_t i = 0; i < samples.size(); ++i)
    {
      const double frac = (total > 0.0 && i < contrib.size()) ? contrib[i] / total : std::numeric_limits<double>::quiet_NaN();
      std::cout << std::setw(10) << Fmt(frac, 3);
    }
    std::cout << "\n";

    if (std::isfinite(prevQ80) && std::isfinite(q80))
    {
      roughness += std::fabs(q80 - prevQ80);
      ++nSteps;
    }
    if (std::isfinite(q80)) prevQ80 = q80;
  }

  std::cout << "roughness_sum_abs_delta_q80=" << Fmt(roughness, 3)
            << " over " << nSteps << " adjacent finite steps\n";
}

void PrintCompare2D(const std::vector<Sample>& samples)
{
  PrintHeader("G4 stored vs alternative stitch variable 2D summaries");
  std::cout << std::left << std::setw(21) << "hist"
            << std::setw(10) << "sample"
            << std::right << std::setw(12) << "entries"
            << std::setw(12) << "meanX"
            << std::setw(12) << "meanY"
            << std::setw(12) << "meanY-X"
            << "\n";

  const std::vector<std::string> keys = {"hepmcParent", "hepmcParentStable", "hepmcAnyPhoton", "g4NoEmbed"};
  for (const auto& key : keys)
  {
    const std::string hname = "h2_ppPhotonStitch_" + key + "_vs_g4Stored";
    for (const auto& s : samples)
    {
      TH2* h2 = GetH2(s.dir, hname);
      if (!h2 || h2->GetEntries() <= 0.0)
      {
        std::cout << std::left << std::setw(21) << key
                  << std::setw(10) << s.label
                  << std::right << std::setw(12) << "missing"
                  << "\n";
        continue;
      }
      const double meanX = h2->GetMean(1);
      const double meanY = h2->GetMean(2);
      std::cout << std::left << std::setw(21) << key
                << std::setw(10) << s.label
                << std::right << std::setw(12) << Fmt(h2->GetEntries(), 0)
                << std::setw(12) << Fmt(meanX, 3)
                << std::setw(12) << Fmt(meanY, 3)
                << std::setw(12) << Fmt(meanY - meanX, 3)
                << "\n";
    }
  }
}

bool OpenSamples(std::vector<Sample>& samples)
{
  for (auto& s : samples)
  {
    s.file = TFile::Open(s.path.c_str(), "READ");
    if (!s.file || s.file->IsZombie())
    {
      std::cerr << "[ERROR] Cannot open " << s.label << ": " << s.path << "\n";
      return false;
    }
    s.dir = s.file->GetDirectory("SIM");
    if (!s.dir) s.dir = s.file;
    if (!s.dir)
    {
      std::cerr << "[ERROR] Missing SIM directory in " << s.path << "\n";
      return false;
    }
    s.nRaw = CountFromCntSim(s.dir);
    if (s.nRaw <= 0.0)
    {
      std::cerr << "[ERROR] Missing/zero SIM/cnt_SIM in " << s.path << "\n";
      return false;
    }
  }

  const double ref = samples.back().sigmaPb / samples.back().nRaw;
  if (ref <= 0.0) return false;
  for (auto& s : samples)
  {
    s.relWeight = (s.sigmaPb / s.nRaw) / ref;
  }
  return true;
}
}  // namespace

void PrintPPStitchDiagnostics(const char* photon5Path,
                              const char* photon10Path,
                              const char* photon20Path,
                              double targetRaw = 30.0)
{
  std::vector<Sample> samples = {
      {"pj5", photon5Path ? photon5Path : "", 146359.3, nullptr, nullptr, 0.0, 1.0},
      {"pj10", photon10Path ? photon10Path : "", 6944.675, nullptr, nullptr, 0.0, 1.0},
      {"pj20", photon20Path ? photon20Path : "", 130.4461, nullptr, nullptr, 0.0, 1.0},
  };

  std::cout << std::fixed << std::setprecision(3);
  PrintHeader("PP stitch local diagnostic inputs");
  for (const auto& s : samples)
  {
    std::cout << s.label << " : " << s.path << "\n";
  }

  if (!OpenSamples(samples)) return;

  PrintHeader("Pseudo-merge weights");
  std::cout << std::left << std::setw(10) << "sample"
            << std::right << std::setw(14) << "Nraw"
            << std::setw(16) << "sigma_pb"
            << std::setw(16) << "relWeight"
            << "\n";
  for (const auto& s : samples)
  {
    std::cout << std::left << std::setw(10) << s.label
              << std::right << std::setw(14) << Fmt(s.nRaw, 1)
              << std::setw(16) << Fmt(s.sigmaPb, 6)
              << std::setw(16) << Fmt(s.relWeight, 6)
              << "\n";
  }

  PrintIndependentTable(samples);
  PrintFlowTable(samples, targetRaw);
  PrintStitchSummary(samples);
  PrintCompare2D(samples);

  for (const auto& variant : Variants())
  {
    PrintVariantMergeTable(samples, variant, "vtxW");
    PrintVariantMergeTable(samples, variant, "noVtxW");
  }

  PrintHeader("Interpretation checklist");
  std::cout << "1. If g4Stored has a large q80 jump but hepmcParent does not, the pp stitch variable is the culprit.\n"
            << "2. If vtxW and noVtxW are nearly identical, vertex reweighting is not driving the artifact.\n"
            << "3. If independent pj20 q80 is high in 16-20 GeV but hepmcParent keeps pj20 out of that reco bin, the artifact is sample ownership.\n"
            << "4. Use sample fractions to identify which sample dominates each pseudo-merged reco-pT bin.\n";
}
