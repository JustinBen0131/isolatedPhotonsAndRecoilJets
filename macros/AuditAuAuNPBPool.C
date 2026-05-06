#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "TDirectory.h"
#include "TFile.h"
#include "TH1.h"
#include "TString.h"

namespace
{
const std::vector<std::pair<int,int>> kPtBins = {
  {10,12}, {12,14}, {14,16}, {16,18}, {18,20}, {20,22}, {22,24}, {24,26}, {26,35}
};
const std::vector<std::string> kCents = {"0_20", "20_50", "50_80"};

double Count(TDirectory* d, const std::string& name)
{
  TH1* h = d ? dynamic_cast<TH1*>(d->Get(name.c_str())) : nullptr;
  return h ? h->GetBinContent(1) : 0.0;
}

TH1* Hist(TDirectory* d, const std::string& name)
{
  return d ? dynamic_cast<TH1*>(d->Get(name.c_str())) : nullptr;
}

std::string Suffix(int lo, int hi, const std::string& cent)
{
  return TString::Format("_pT_%d_%d_cent_%s", lo, hi, cent.c_str()).Data();
}

double RawPurity(double A, double B, double C, double D)
{
  if (A <= 0.0 || D <= 0.0) return 0.0;
  double s = A - B * C / D;
  if (s < 0.0) s = 0.0;
  return s / A;
}

struct ABCD
{
  double A = 0.0, B = 0.0, C = 0.0, D = 0.0;
  double purity() const { return RawPurity(A, B, C, D); }
};

ABCD ReadABCD(TDirectory* d, const std::string& suffix)
{
  ABCD x;
  x.A = Count(d, "h_isIsolated_isTight" + suffix);
  x.B = Count(d, "h_notIsolated_isTight" + suffix);
  x.C = Count(d, "h_isIsolated_notTight" + suffix);
  x.D = Count(d, "h_notIsolated_notTight" + suffix);
  return x;
}

struct Matrix
{
  double pp = 0.0; // npb pass, reference pass
  double pf = 0.0; // npb pass, reference fail
  double fp = 0.0; // npb fail, reference pass
  double ff = 0.0; // npb fail, reference fail
  double total() const { return pp + pf + fp + ff; }
};

Matrix ReadMatrix(TDirectory* d, const std::string& suffix)
{
  Matrix m;
  m.pp = Count(d, "h_preRefMatrix_npbPass_refPass" + suffix);
  m.pf = Count(d, "h_preRefMatrix_npbPass_refFail" + suffix);
  m.fp = Count(d, "h_preRefMatrix_npbFail_refPass" + suffix);
  m.ff = Count(d, "h_preRefMatrix_npbFail_refFail" + suffix);
  return m;
}

struct MaskEntry
{
  int mask = 0;
  double n = 0.0;
};

std::vector<MaskEntry> ReadMaskEntries(TDirectory* d, const std::string& name)
{
  std::vector<MaskEntry> out;
  TH1* h = Hist(d, name);
  if (!h) return out;
  for (int ib = 1; ib <= h->GetNbinsX(); ++ib)
  {
    const double n = h->GetBinContent(ib);
    if (n <= 0.0) continue;
    int mask = static_cast<int>(std::llround(h->GetBinCenter(ib)));
    out.push_back({mask, n});
  }
  return out;
}

std::string MaskLabel(int mask)
{
  if (mask == 0) return "ref-pass";
  std::vector<std::string> bits;
  if (mask & 1) bits.push_back("e11/e33");
  if (mask & 2) bits.push_back("et1");
  if (mask & 4) bits.push_back("e32/e35");
  if (mask & 8) bits.push_back("weta");
  std::ostringstream os;
  for (std::size_t i = 0; i < bits.size(); ++i)
  {
    if (i) os << "+";
    os << bits[i];
  }
  return os.str();
}

void AddInto(std::unique_ptr<TH1>& acc, TH1* h)
{
  if (!h) return;
  if (!acc)
  {
    acc.reset(dynamic_cast<TH1*>(h->Clone()));
    acc->SetDirectory(nullptr);
  }
  else acc->Add(h);
}

std::string CentLabel(const std::string& cent)
{
  std::string out = cent;
  std::replace(out.begin(), out.end(), '_', '-');
  return out + "%";
}
}

void AuditAuAuNPBPool()
{
  const std::string refFile =
      "InputFiles/auau25/RecoilJets_auau_ALL_jetMinPt5_7pi_8_vz60_isoR40_isSliding_baseVariant_preselectionReference_tightReference_nonTightReference.root";
  const std::string npbFile =
      "InputFiles/auau25/RecoilJets_auau_ALL_jetMinPt5_7pi_8_vz60_isoR40_isSliding_baseVariant_preselectionVariantA_tightReference_nonTightReference.root";
  const std::string trig = "photon_10_plus_MBD_NS_geq_2_vtx_lt_150";

  TFile fRef(refFile.c_str(), "READ");
  TFile fNPB(npbFile.c_str(), "READ");
  TDirectory* dRef = fRef.GetDirectory(trig.c_str());
  TDirectory* dNPB = fNPB.GetDirectory(trig.c_str());
  if (!dRef || !dNPB)
  {
    std::cerr << "Missing trigger directory\n";
    return;
  }

  std::cout << "\n=== NPB vs old reference preselection matrix, integrated pT 10-35 ===\n";
  std::cout << std::setw(8) << "cent"
            << std::setw(12) << "refPass"
            << std::setw(12) << "npbPass"
            << std::setw(12) << "agree%"
            << std::setw(14) << "recovered"
            << std::setw(14) << "rejectRef"
            << std::setw(14) << "netPass%"
            << "\n";

  for (const auto& cent : kCents)
  {
    Matrix sum;
    for (const auto& p : kPtBins)
    {
      const auto m = ReadMatrix(dNPB, Suffix(p.first, p.second, cent));
      sum.pp += m.pp; sum.pf += m.pf; sum.fp += m.fp; sum.ff += m.ff;
    }
    const double refPass = sum.pp + sum.fp;
    const double npbPass = sum.pp + sum.pf;
    const double total = sum.total();
    const double agree = total > 0.0 ? (sum.pp + sum.ff) / total : 0.0;
    const double recovered = refPass > 0.0 ? sum.pf / refPass : 0.0;
    const double rejectRef = refPass > 0.0 ? sum.fp / refPass : 0.0;
    const double netPass = refPass > 0.0 ? (npbPass / refPass - 1.0) : 0.0;

    std::cout << std::setw(8) << CentLabel(cent)
              << std::setw(12) << std::fixed << std::setprecision(0) << refPass
              << std::setw(12) << npbPass
              << std::setw(11) << std::setprecision(1) << 100.0 * agree
              << "%"
              << std::setw(13) << 100.0 * recovered << "%"
              << std::setw(13) << 100.0 * rejectRef << "%"
              << std::setw(13) << 100.0 * netPass << "%"
              << "\n";
  }

  std::cout << "\n=== Matrix by pT: recovered old-fail and rejected old-pass fractions ===\n";
  for (const auto& cent : kCents)
  {
    std::cout << "\ncentrality " << CentLabel(cent) << "\n";
    std::cout << std::setw(8) << "pT"
              << std::setw(12) << "refPass"
              << std::setw(13) << "recover%"
              << std::setw(13) << "reject%"
              << std::setw(13) << "netPass%"
              << "\n";
    for (const auto& p : kPtBins)
    {
      const auto m = ReadMatrix(dNPB, Suffix(p.first, p.second, cent));
      const double refPass = m.pp + m.fp;
      const double npbPass = m.pp + m.pf;
      if (m.total() <= 0.0) continue;
      std::cout << std::setw(8) << TString::Format("%d-%d", p.first, p.second).Data()
                << std::setw(12) << std::fixed << std::setprecision(0) << refPass
                << std::setw(12) << std::setprecision(1) << (refPass > 0 ? 100.0*m.pf/refPass : 0.0) << "%"
                << std::setw(12) << (refPass > 0 ? 100.0*m.fp/refPass : 0.0) << "%"
                << std::setw(12) << (refPass > 0 ? 100.0*(npbPass/refPass - 1.0) : 0.0) << "%"
                << "\n";
    }
  }

  std::cout << "\n=== Downstream ABCD effect: NPB presel / box-cuts, integrated pT 10-35 ===\n";
  std::cout << std::setw(8) << "cent"
            << std::setw(10) << "A ratio"
            << std::setw(10) << "B ratio"
            << std::setw(10) << "C ratio"
            << std::setw(10) << "D ratio"
            << std::setw(12) << "pur box"
            << std::setw(12) << "pur NPB"
            << std::setw(12) << "dPur"
            << "\n";
  for (const auto& cent : kCents)
  {
    ABCD r, n;
    for (const auto& p : kPtBins)
    {
      const auto sr = ReadABCD(dRef, Suffix(p.first, p.second, cent));
      const auto sn = ReadABCD(dNPB, Suffix(p.first, p.second, cent));
      r.A += sr.A; r.B += sr.B; r.C += sr.C; r.D += sr.D;
      n.A += sn.A; n.B += sn.B; n.C += sn.C; n.D += sn.D;
    }
    std::cout << std::setw(8) << CentLabel(cent)
              << std::setw(10) << std::fixed << std::setprecision(3) << (r.A > 0 ? n.A/r.A : 0)
              << std::setw(10) << (r.B > 0 ? n.B/r.B : 0)
              << std::setw(10) << (r.C > 0 ? n.C/r.C : 0)
              << std::setw(10) << (r.D > 0 ? n.D/r.D : 0)
              << std::setw(12) << r.purity()
              << std::setw(12) << n.purity()
              << std::setw(12) << n.purity() - r.purity()
              << "\n";
  }

  std::cout << "\n=== Top old-reference fail patterns among NPB-pass recovered candidates ===\n";
  for (const auto& cent : kCents)
  {
    std::map<int,double> maskCounts;
    for (const auto& p : kPtBins)
    {
      for (const auto& e : ReadMaskEntries(dNPB, "h_preRefFailMask_npbPass" + Suffix(p.first, p.second, cent)))
        if (e.mask != 0) maskCounts[e.mask] += e.n;
    }
    std::vector<MaskEntry> v;
    for (const auto& kv : maskCounts) v.push_back({kv.first, kv.second});
    std::sort(v.begin(), v.end(), [](const MaskEntry& a, const MaskEntry& b){ return a.n > b.n; });
    double total = 0.0;
    for (const auto& e : v) total += e.n;
    std::cout << "\ncentrality " << CentLabel(cent) << ", recovered old-fail total=" << total << "\n";
    for (std::size_t i = 0; i < std::min<std::size_t>(6, v.size()); ++i)
    {
      std::cout << "  " << std::setw(22) << MaskLabel(v[i].mask)
                << "  " << std::setw(10) << std::fixed << std::setprecision(0) << v[i].n
                << "  " << std::setw(6) << std::setprecision(1) << (total > 0 ? 100.0*v[i].n/total : 0.0) << "%\n";
    }
  }

  std::cout << "\n=== NPB pass/fail shower-shape means, integrated pT 10-35 ===\n";
  const std::vector<std::string> vars = {"npbScore", "weta", "wphi", "et1", "e11e33", "e32e35", "weta35", "wphi53"};
  for (const auto& cent : kCents)
  {
    std::cout << "\ncentrality " << CentLabel(cent) << "\n";
    std::cout << std::setw(10) << "var"
              << std::setw(14) << "passMean"
              << std::setw(14) << "failMean"
              << std::setw(14) << "delta"
              << std::setw(14) << "passN"
              << std::setw(14) << "failN"
              << "\n";
    for (const auto& var : vars)
    {
      std::unique_ptr<TH1> hp, hf;
      for (const auto& p : kPtBins)
      {
        const std::string s = Suffix(p.first, p.second, cent);
        if (var == "npbScore")
        {
          AddInto(hp, Hist(dNPB, "h_ss_npbScore_npbPass" + s));
          AddInto(hf, Hist(dNPB, "h_ss_npbScore_npbFail" + s));
        }
        else
        {
          AddInto(hp, Hist(dNPB, "h_ss_" + var + "_npbPass" + s));
          AddInto(hf, Hist(dNPB, "h_ss_" + var + "_npbFail" + s));
        }
      }
      if (!hp || !hf || hp->GetEntries() <= 0 || hf->GetEntries() <= 0) continue;
      std::cout << std::setw(10) << var
                << std::setw(14) << std::fixed << std::setprecision(4) << hp->GetMean()
                << std::setw(14) << hf->GetMean()
                << std::setw(14) << hp->GetMean() - hf->GetMean()
                << std::setw(14) << std::setprecision(0) << hp->GetEntries()
                << std::setw(14) << hf->GetEntries()
                << "\n";
    }
  }
}
