#include "AnalyzeRecoilJets.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TSystem.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace
{
struct VariantSpec
{
    std::string label;
    std::string tag;
    int color;
    int marker;
};

double ReadBin1(TDirectory* dir, const std::string& name)
{
    if (!dir) return 0.0;
    TH1* h = dynamic_cast<TH1*>(dir->Get(name.c_str()));
    return h ? h->GetBinContent(1) : 0.0;
}

double RawPurityError(double A, double B, double C, double D)
{
    if (A <= 0.0 || D <= 0.0) return 0.0;

    const double dPdA =  (B * C) / (A * A * D);
    const double dPdB = -(C) / (A * D);
    const double dPdC = -(B) / (A * D);
    const double dPdD =  (B * C) / (A * D * D);

    double var = 0.0;
    if (A > 0.0) var += dPdA * dPdA * A;
    if (B > 0.0) var += dPdB * dPdB * B;
    if (C > 0.0) var += dPdC * dPdC * C;
    if (D > 0.0) var += dPdD * dPdD * D;

    return (var > 0.0) ? std::sqrt(var) : 0.0;
}

std::string PPPath(const std::string& tag)
{
    return ARJ::kInputBase + "/pp24/RecoilJets_pp_ALL_" + tag + ".root";
}

std::string SimPath(const std::string& sample, const std::string& tag)
{
    return ARJ::kInputBase + "/simPhotonJet/RecoilJets_" + sample + "_ALL_" + tag + ".root";
}

std::string MergedSimPath(const std::string& tag)
{
    return ARJ::kOutputBase + "/combinedSimOnly/" + tag +
           "/photonJet5and10and20merged_SIM/RecoilJets_photonjet5plus10plus20_MERGED.root";
}

std::string ExtractTokenAfter(const std::string& text, const std::string& key)
{
    const std::size_t pos = text.find(key);
    if (pos == std::string::npos) return "";
    const std::size_t begin = pos + key.size();
    std::size_t end = begin;
    while (end < text.size() && text[end] != '_') ++end;
    return text.substr(begin, end - begin);
}

std::string FormatFixed(double value)
{
    std::ostringstream os;
    os << std::fixed << std::setprecision(1) << value;
    std::string out = os.str();
    while (out.size() > 1 && out.back() == '0') out.pop_back();
    if (!out.empty() && out.back() == '.') out.pop_back();
    return out;
}

std::string PPTriggerLabel()
{
    const std::string& trig = ARJ::kTriggerPP;
    std::size_t photonPos = trig.find("Photon_");
    std::size_t gevPos = trig.find("_GeV", photonPos == std::string::npos ? 0 : photonPos);
    std::string photon = "4";
    if (photonPos != std::string::npos && gevPos != std::string::npos)
        photon = trig.substr(photonPos + 7, gevPos - (photonPos + 7));

    std::string mbd = "1";
    const std::string mbdKey = "MBD_NS_geq_";
    std::size_t mbdPos = trig.find(mbdKey);
    if (mbdPos != std::string::npos)
        mbd = ExtractTokenAfter(trig.substr(mbdPos), mbdKey);

    return "Photon " + photon + " GeV + MBD N&S #geq " + mbd;
}

std::vector<std::string> CutLinesFromInputTag(const std::string& tag)
{
    int vzCut = ARJ::kVzCut;
    const std::string vzText = ExtractTokenAfter(tag, "vz");
    if (!vzText.empty()) vzCut = std::atoi(vzText.c_str());

    double coneR = (ARJ::kIsoConeR == "isoR40") ? 0.4 : 0.3;
    const std::string isoRText = ExtractTokenAfter(tag, "isoR");
    if (!isoRText.empty()) coneR = std::atoi(isoRText.c_str()) / 100.0;

    std::string isoLabel = "sliding E_{T}^{iso} cut";
    const std::string fixedKey = "fixedIso";
    const std::size_t fixedPos = tag.find(fixedKey);
    if (fixedPos != std::string::npos)
    {
        const std::size_t valBegin = fixedPos + fixedKey.size();
        const std::size_t valEnd = tag.find("GeV", valBegin);
        if (valEnd != std::string::npos && valEnd > valBegin)
        {
            const double fixedIso = std::atof(tag.substr(valBegin, valEnd - valBegin).c_str());
            isoLabel = "E_{T}^{iso} < " + FormatFixed(fixedIso) + " GeV";
        }
    }

    return {
        PPTriggerLabel(),
        "|v_{z}| < " + std::to_string(vzCut) + " cm, #Delta R_{cone} < " +
            FormatFixed(coneR) + ", " + isoLabel
    };
}

bool EnsureMergedSim(const VariantSpec& spec)
{
    const std::string out = MergedSimPath(spec.tag);
    if (!gSystem->AccessPathName(out.c_str()))
    {
        std::cout << "[OK] merged SIM exists for " << spec.label << ": " << out << "\n";
        return true;
    }

    const std::vector<std::string> inputs = {
        SimPath("photonjet5", spec.tag),
        SimPath("photonjet10", spec.tag),
        SimPath("photonjet20", spec.tag)
    };

    for (const auto& path : inputs)
    {
        if (gSystem->AccessPathName(path.c_str()))
        {
            std::cerr << "[ERROR] Missing SIM input for " << spec.label << ": " << path << "\n";
            return false;
        }
    }

    std::cout << "[MERGE] Building 5+10+20 merged SIM for " << spec.label << "\n";
    return ARJ::BuildMergedSIMFile_PhotonSlices(
        inputs,
        {ARJ::kSigmaPhoton5_pb, ARJ::kSigmaPhoton10_pb, ARJ::kSigmaPhoton20_pb},
        out,
        ARJ::kDirSIM,
        {"photonJet5", "photonJet10", "photonJet20"});
}

bool BuildPurityGraph(const VariantSpec& spec,
                      std::vector<double>& x,
                      std::vector<double>& ex,
                      std::vector<double>& y,
                      std::vector<double>& ey)
{
    TFile f(PPPath(spec.tag).c_str(), "READ");
    if (f.IsZombie())
    {
        std::cerr << "[ERROR] Cannot open pp input for " << spec.label << ": " << PPPath(spec.tag) << "\n";
        return false;
    }

    TDirectory* top = f.GetDirectory(ARJ::kTriggerPP.c_str());
    if (!top) top = &f;

    const int n = ARJ::kNPtBins - 1;
    x.assign(n, 0.0);
    ex.assign(n, 0.0);
    y.assign(n, 0.0);
    ey.assign(n, 0.0);

    for (int i = 1; i < ARJ::kNPtBins; ++i)
    {
        const int j = i - 1;
        const auto& b = ARJ::PtBins()[i];
        const std::string& suf = b.suffix;

        const double A = ReadBin1(top, "h_isIsolated_isTight" + suf);
        const double B = ReadBin1(top, "h_notIsolated_isTight" + suf);
        const double C = ReadBin1(top, "h_isIsolated_notTight" + suf);
        const double D = ReadBin1(top, "h_notIsolated_notTight" + suf);

        const double ptLo = ARJ::kPtEdges[(std::size_t)i];
        const double ptHi = ARJ::kPtEdges[(std::size_t)i + 1];
        x[j] = 0.5 * (ptLo + ptHi);
        ex[j] = 0.5 * (ptHi - ptLo);

        if (A > 0.0 && D > 0.0)
        {
            double sig = A - B * (C / D);
            if (sig < 0.0) sig = 0.0;
            y[j] = sig / A;
            ey[j] = RawPurityError(A, B, C, D);
        }
    }

    return true;
}

bool BuildLeakageCGraph(const VariantSpec& spec,
                        std::vector<double>& x,
                        std::vector<double>& ex,
                        std::vector<double>& y,
                        std::vector<double>& ey)
{
    if (!EnsureMergedSim(spec)) return false;

    TFile f(MergedSimPath(spec.tag).c_str(), "READ");
    if (f.IsZombie())
    {
        std::cerr << "[ERROR] Cannot open merged SIM for " << spec.label << ": " << MergedSimPath(spec.tag) << "\n";
        return false;
    }

    TDirectory* top = f.GetDirectory(ARJ::kDirSIM.c_str());
    if (!top) top = &f;

    const int n = ARJ::kNPtBins - 1;
    x.assign(n, 0.0);
    ex.assign(n, 0.0);
    y.assign(n, 0.0);
    ey.assign(n, 0.0);

    for (int i = 1; i < ARJ::kNPtBins; ++i)
    {
        const int j = i - 1;
        const auto& b = ARJ::PtBins()[i];
        const std::string hname = "h_sigABCD_MC" + b.suffix;
        TH1* h = dynamic_cast<TH1*>(top->Get(hname.c_str()));

        const double ptLo = ARJ::kPtEdges[(std::size_t)i];
        const double ptHi = ARJ::kPtEdges[(std::size_t)i + 1];
        x[j] = 0.5 * (ptLo + ptHi);
        ex[j] = 0.5 * (ptHi - ptLo);

        if (!h)
        {
            std::cerr << "[WARN] Missing " << hname << " in merged SIM for " << spec.label << "\n";
            continue;
        }

        const double As = h->GetBinContent(1);
        const double Cs = h->GetBinContent(3);
        y[j] = (As > 0.0) ? (Cs / As) : 0.0;
    }

    return true;
}

double CorrectedPurityValue(double A, double B, double C, double D,
                            double fB, double fC, double fD)
{
    double raw = 0.0;
    if (A > 0.0 && D > 0.0)
    {
        double sig = A - B * (C / D);
        if (sig < 0.0) sig = 0.0;
        raw = sig / A;
    }

    double SA = 0.0;
    const bool ok = ARJ::SolveLeakageCorrectedSA(A, B, C, D, fB, fC, fD, SA);
    if (ok && A > 0.0) return SA / A;
    return raw;
}

double CorrectedPurityError(double A, double B, double C, double D,
                            double fB, double fC, double fD,
                            bool hasCorrection)
{
    if (!hasCorrection) return RawPurityError(A, B, C, D);
    if (A <= 0.0) return 0.0;

    const double dA = std::sqrt(std::max(A, 1.0));
    const double dB = std::sqrt(std::max(B, 1.0));
    const double dC = std::sqrt(std::max(C, 1.0));
    const double dD = std::sqrt(std::max(D, 1.0));

    const double Aup = A + dA;
    const double Adn = std::max(0.0, A - dA);
    const double Bup = B + dB;
    const double Bdn = std::max(0.0, B - dB);
    const double Cup = C + dC;
    const double Cdn = std::max(0.0, C - dC);
    const double Dup = D + dD;
    const double Ddn = std::max(0.0, D - dD);

    const double dPdA = (Aup > Adn)
        ? (CorrectedPurityValue(Aup, B,   C,   D,   fB, fC, fD) -
           CorrectedPurityValue(Adn, B,   C,   D,   fB, fC, fD)) / (Aup - Adn)
        : 0.0;
    const double dPdB = (Bup > Bdn)
        ? (CorrectedPurityValue(A,   Bup, C,   D,   fB, fC, fD) -
           CorrectedPurityValue(A,   Bdn, C,   D,   fB, fC, fD)) / (Bup - Bdn)
        : 0.0;
    const double dPdC = (Cup > Cdn)
        ? (CorrectedPurityValue(A,   B,   Cup, D,   fB, fC, fD) -
           CorrectedPurityValue(A,   B,   Cdn, D,   fB, fC, fD)) / (Cup - Cdn)
        : 0.0;
    const double dPdD = (Dup > Ddn)
        ? (CorrectedPurityValue(A,   B,   C,   Dup, fB, fC, fD) -
           CorrectedPurityValue(A,   B,   C,   Ddn, fB, fC, fD)) / (Dup - Ddn)
        : 0.0;

    double var = 0.0;
    if (A > 0.0) var += dPdA * dPdA * A;
    if (B > 0.0) var += dPdB * dPdB * B;
    if (C > 0.0) var += dPdC * dPdC * C;
    if (D > 0.0) var += dPdD * dPdD * D;

    return (var > 0.0) ? std::sqrt(var) : 0.0;
}

bool BuildLeakageCorrectedPurityGraph(const VariantSpec& spec,
                                      std::vector<double>& x,
                                      std::vector<double>& ex,
                                      std::vector<double>& y,
                                      std::vector<double>& ey)
{
    if (!EnsureMergedSim(spec)) return false;

    TFile fPP(PPPath(spec.tag).c_str(), "READ");
    if (fPP.IsZombie())
    {
        std::cerr << "[ERROR] Cannot open pp input for " << spec.label << ": " << PPPath(spec.tag) << "\n";
        return false;
    }

    TDirectory* ppTop = fPP.GetDirectory(ARJ::kTriggerPP.c_str());
    if (!ppTop) ppTop = &fPP;

    TFile fSim(MergedSimPath(spec.tag).c_str(), "READ");
    if (fSim.IsZombie())
    {
        std::cerr << "[ERROR] Cannot open merged SIM for " << spec.label << ": " << MergedSimPath(spec.tag) << "\n";
        return false;
    }

    TDirectory* simTop = fSim.GetDirectory(ARJ::kDirSIM.c_str());
    if (!simTop) simTop = &fSim;

    const int n = ARJ::kNPtBins - 1;
    x.assign(n, 0.0);
    ex.assign(n, 0.0);
    y.assign(n, 0.0);
    ey.assign(n, 0.0);

    for (int i = 1; i < ARJ::kNPtBins; ++i)
    {
        const int j = i - 1;
        const auto& b = ARJ::PtBins()[i];
        const std::string& suf = b.suffix;

        const double A = ReadBin1(ppTop, "h_isIsolated_isTight" + suf);
        const double B = ReadBin1(ppTop, "h_notIsolated_isTight" + suf);
        const double C = ReadBin1(ppTop, "h_isIsolated_notTight" + suf);
        const double D = ReadBin1(ppTop, "h_notIsolated_notTight" + suf);

        const double ptLo = ARJ::kPtEdges[(std::size_t)i];
        const double ptHi = ARJ::kPtEdges[(std::size_t)i + 1];
        x[j] = 0.5 * (ptLo + ptHi);
        ex[j] = 0.5 * (ptHi - ptLo);

        double fB = 0.0;
        double fC = 0.0;
        double fD = 0.0;
        bool hasCorrection = false;
        TH1* hSig = dynamic_cast<TH1*>(simTop->Get(("h_sigABCD_MC" + suf).c_str()));
        if (hSig)
        {
            const double As = hSig->GetBinContent(1);
            const double Bs = hSig->GetBinContent(2);
            const double Cs = hSig->GetBinContent(3);
            const double Ds = hSig->GetBinContent(4);
            if (As > 0.0)
            {
                fB = Bs / As;
                fC = Cs / As;
                fD = Ds / As;

                double SA = 0.0;
                hasCorrection = ARJ::SolveLeakageCorrectedSA(A, B, C, D, fB, fC, fD, SA) && A > 0.0;
                y[j] = hasCorrection ? (SA / A) : CorrectedPurityValue(A, B, C, D, fB, fC, fD);
            }
            else
            {
                y[j] = CorrectedPurityValue(A, B, C, D, fB, fC, fD);
            }
        }
        else
        {
            std::cerr << "[WARN] Missing h_sigABCD_MC" << suf << " in merged SIM for " << spec.label
                      << "; using raw purity fallback\n";
            y[j] = CorrectedPurityValue(A, B, C, D, fB, fC, fD);
        }

        ey[j] = CorrectedPurityError(A, B, C, D, fB, fC, fD, hasCorrection);
    }

    return true;
}

void DrawOverlay(const std::vector<VariantSpec>& specs,
                 const std::vector<std::vector<double>>& xs,
                 const std::vector<std::vector<double>>& exs,
                 const std::vector<std::vector<double>>& ys,
                 const std::vector<std::vector<double>>& eys,
                 const std::string& yTitle,
                 const std::string& title,
                 const std::string& outPath,
                 double yMax,
                 const std::string& selectionTag = "")
{
    TCanvas c("c_overlay", title.c_str(), 1000, 760);
    c.SetLeftMargin(0.13);
    c.SetRightMargin(0.04);
    c.SetBottomMargin(0.13);
    c.SetTopMargin(0.11);
    c.SetTicks(1, 1);

    TH1F frame("frame", "", 100, ARJ::kPtEdges[1], ARJ::kPtEdges.back());
    frame.SetDirectory(nullptr);
    frame.SetStats(0);
    frame.SetMinimum(0.0);
    frame.SetMaximum(yMax);
    frame.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
    frame.GetYaxis()->SetTitle(yTitle.c_str());
    frame.GetXaxis()->SetTitleSize(0.045);
    frame.GetYaxis()->SetTitleSize(0.045);
    frame.GetXaxis()->SetLabelSize(0.040);
    frame.GetYaxis()->SetLabelSize(0.040);
    frame.Draw();

    std::vector<TGraphErrors*> graphs;
    for (std::size_t i = 0; i < specs.size(); ++i)
    {
        if (xs[i].empty()) continue;
        auto* g = new TGraphErrors((int)xs[i].size(),
                                   const_cast<double*>(xs[i].data()),
                                   const_cast<double*>(ys[i].data()),
                                   const_cast<double*>(exs[i].data()),
                                   const_cast<double*>(eys[i].data()));
        g->SetLineColor(specs[i].color);
        g->SetMarkerColor(specs[i].color);
        g->SetLineWidth(2);
        g->SetMarkerStyle(specs[i].marker);
        g->SetMarkerSize(1.15);
        g->Draw("PE1 SAME");
        graphs.push_back(g);
    }

    TLegend leg(0.70, 0.70, 0.96, 0.89);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetTextFont(42);
    leg.SetTextSize(0.034);
    for (std::size_t i = 0; i < graphs.size() && i < specs.size(); ++i)
        leg.AddEntry(graphs[i], specs[i].label.c_str(), "pe");
    leg.Draw();

    TLatex text;
    text.SetNDC(true);
    text.SetTextFont(42);
    text.SetTextAlign(22);
    text.SetTextSize(0.040);
    text.DrawLatex(0.50, 0.935, title.c_str());

    const std::string tagForLabels =
        selectionTag.empty() ? (specs.empty() ? std::string() : specs.front().tag) : selectionTag;
    const std::vector<std::string> cutLines = CutLinesFromInputTag(tagForLabels);

    text.SetTextAlign(13);
    text.SetTextSize(0.033);
    text.DrawLatex(0.18, 0.855, cutLines[0].c_str());
    text.DrawLatex(0.18, 0.805, cutLines[1].c_str());

    c.SaveAs(outPath.c_str());
    std::cout << "[WROTE] " << outPath << "\n";

    for (auto* g : graphs) delete g;
}
}

void OverlayPPPurityAndLeakageVariants()
{
    const std::vector<VariantSpec> specs = {
        {"box-cuts",
         "jetMinPt5_7pi_8_vz60_isoR40_isSliding_preselectionReference_tightReference_nonTightReference",
         kBlack, 20},
        {"BDT with NPB pre",
         "jetMinPt5_7pi_8_vz60_isoR40_isSliding_preselectionVariantA_tightVariantA_nonTightReference",
         kBlue + 1, 21},
        {"BDT only",
         "jetMinPt5_7pi_8_vz60_isoR40_isSliding_preselectionVariantB_tightVariantA_nonTightReference",
         kRed + 1, 22}
    };

    gSystem->mkdir((ARJ::kOutputBase + "/pp").c_str(), true);

    std::vector<std::vector<double>> xs(specs.size()), exs(specs.size());
    std::vector<std::vector<double>> ys(specs.size()), eys(specs.size());
    for (std::size_t i = 0; i < specs.size(); ++i)
    {
        if (!BuildPurityGraph(specs[i], xs[i], exs[i], ys[i], eys[i]))
            std::cerr << "[WARN] Skipping purity graph for " << specs[i].label << "\n";
    }
    DrawOverlay(specs, xs, exs, ys, eys,
                "Purity (raw ABCD)",
                "Raw ABCD purity overlay",
                ARJ::kOutputBase + "/pp/rawPurity_overlay_isSliding_photonIDVariants.png",
                1.05);

    for (auto& v : xs) v.clear();
    for (auto& v : exs) v.clear();
    for (auto& v : ys) v.clear();
    for (auto& v : eys) v.clear();

    for (std::size_t i = 0; i < specs.size(); ++i)
    {
        if (!BuildLeakageCGraph(specs[i], xs[i], exs[i], ys[i], eys[i]))
            std::cerr << "[WARN] Skipping region-C leakage graph for " << specs[i].label << "\n";
    }
    DrawOverlay(specs, xs, exs, ys, eys,
                "Region C signal leakage, f_{C} = C_{sig}/A_{sig}",
                "Truth-signal leakage into ABCD region C",
                ARJ::kOutputBase + "/pp/signalLeakage_regionC_overlay_isSliding_photonIDVariants.png",
                1.05);
}

void OverlayPPFixedIsoLeakageReferenceVsVariantA()
{
    const std::vector<VariantSpec> specs = {
        {"Box-cuts",
         "jetMinPt5_7pi_8_vz60_isoR40_fixedIso2GeV_preselectionReference_tightReference_nonTightReference",
         kBlack, 20},
        {"BDT + NPB",
         "jetMinPt5_7pi_8_vz60_isoR40_fixedIso2GeV_preselectionVariantA_tightVariantA_nonTightVariantA",
         kBlue + 1, 21}
    };

    gSystem->mkdir((ARJ::kOutputBase + "/pp").c_str(), true);

    std::vector<std::vector<double>> xs(specs.size()), exs(specs.size());
    std::vector<std::vector<double>> ys(specs.size()), eys(specs.size());

    for (std::size_t i = 0; i < specs.size(); ++i)
    {
        if (!BuildLeakageCGraph(specs[i], xs[i], exs[i], ys[i], eys[i]))
            std::cerr << "[WARN] Skipping region-C leakage graph for " << specs[i].label << "\n";
    }

    DrawOverlay(specs, xs, exs, ys, eys,
                "Region C signal leakage, f_{C} = C_{sig}/A_{sig}",
                "Truth-signal Leakage into ABCD Region C - Run24pp",
                ARJ::kOutputBase + "/pp/signalLeakage_regionC_overlay_fixedIso2GeV_reference_vs_variantA.png",
                1.05,
                specs.front().tag);
}

void OverlayPPFixedIsoRawPurityReferenceVsVariantA()
{
    const std::vector<VariantSpec> specs = {
        {"Box-cuts",
         "jetMinPt5_7pi_8_vz60_isoR40_fixedIso2GeV_preselectionReference_tightReference_nonTightReference",
         kBlack, 20},
        {"BDT + NPB",
         "jetMinPt5_7pi_8_vz60_isoR40_fixedIso2GeV_preselectionVariantA_tightVariantA_nonTightVariantA",
         kBlue + 1, 21}
    };

    gSystem->mkdir((ARJ::kOutputBase + "/pp").c_str(), true);

    std::vector<std::vector<double>> xs(specs.size()), exs(specs.size());
    std::vector<std::vector<double>> ys(specs.size()), eys(specs.size());

    for (std::size_t i = 0; i < specs.size(); ++i)
    {
        if (!BuildPurityGraph(specs[i], xs[i], exs[i], ys[i], eys[i]))
            std::cerr << "[WARN] Skipping raw purity graph for " << specs[i].label << "\n";
    }

    DrawOverlay(specs, xs, exs, ys, eys,
                "Purity (raw ABCD)",
                "Raw ABCD Purity - Run24pp",
                ARJ::kOutputBase + "/pp/rawPurity_overlay_fixedIso2GeV_reference_vs_variantA.png",
                1.05,
                specs.front().tag);
}

void OverlayPPFixedIsoCorrectedPurityReferenceVsVariantA()
{
    const std::vector<VariantSpec> specs = {
        {"Box-cuts",
         "jetMinPt5_7pi_8_vz60_isoR40_fixedIso2GeV_preselectionReference_tightReference_nonTightReference",
         kBlack, 20},
        {"BDT + NPB",
         "jetMinPt5_7pi_8_vz60_isoR40_fixedIso2GeV_preselectionVariantA_tightVariantA_nonTightVariantA",
         kBlue + 1, 21}
    };

    gSystem->mkdir((ARJ::kOutputBase + "/pp").c_str(), true);

    std::vector<std::vector<double>> xs(specs.size()), exs(specs.size());
    std::vector<std::vector<double>> ys(specs.size()), eys(specs.size());

    for (std::size_t i = 0; i < specs.size(); ++i)
    {
        if (!BuildLeakageCorrectedPurityGraph(specs[i], xs[i], exs[i], ys[i], eys[i]))
            std::cerr << "[WARN] Skipping leakage-corrected purity graph for " << specs[i].label << "\n";
    }

    DrawOverlay(specs, xs, exs, ys, eys,
                "Purity (leakage-corrected)",
                "Leakage Corrected ABCD Purity - Run24pp",
                ARJ::kOutputBase + "/pp/purity_leakageCorrected_overlay_fixedIso2GeV_reference_vs_variantA.png",
                1.05,
                specs.front().tag);
}
