// ===================================================================
// Au+Au DATA isolation QA: UE-subtraction variant overlay comparison
//   Output:
//     auau/<CfgTagAA>/<trigger>/<cent>/isoQA/UEcomparisons/<variant>/<pT>/
// ===================================================================
// embeddedMode: 0 = AuAu data, 1 = photon+jet embedded, 2 = inclusive jet embedded
void RunIsoQA_UEComparisons_AuAu(int embeddedMode = 0)
{
    const bool forEmbeddedSim = (embeddedMode != 0);
    const string modeLabel = (embeddedMode == 2) ? "Inclusive Jet Embedded SIM"
    : (embeddedMode == 1) ? "Photon+Jet Embedded SIM"
    : "AuAu";
    
    cout << ANSI_BOLD_CYN << "\n==============================\n"
    << "[ISO QA] " << modeLabel
    << " UE-subtraction variant overlay comparisons\n"
    << "==============================" << ANSI_RESET << "\n";
    
    const SimSample activeEmbeddedSample = CurrentSimSample();
    
    // Determine folder name for embedded output
    string activeEmbeddedSimFolder;
    if (embeddedMode == 2)
    {
        if (bothInclusiveJet10and20simEmbedded)
            activeEmbeddedSimFolder = "embeddedJet10and20merged_SIM";
        else if (isInclusiveJet20Embedded)
            activeEmbeddedSimFolder = "embeddedJet20_SIM";
        else if (isInclusiveJet10Embedded)
            activeEmbeddedSimFolder = "embeddedJet10_SIM";
        else
            activeEmbeddedSimFolder = "embeddedJet20_SIM";
    }
    else
    {
        activeEmbeddedSimFolder =
        (activeEmbeddedSample == SimSample::kEmbeddedPhoton10) ? "embeddedPhoton10_SIM" :
        (activeEmbeddedSample == SimSample::kEmbeddedPhoton10And20Merged) ? "photonJet10and20merged_SIM" :
        "embeddedPhoton20_SIM";
    }
    
    auto EmbeddedVariantInput = [&](const string& ueVariant) -> string
    {
        if (embeddedMode == 2)
        {
            if (isInclusiveJet20Embedded)
                return InputInclusiveJetEmbeddedSample("embeddedJet20", ueVariant);
            if (isInclusiveJet10Embedded)
                return InputInclusiveJetEmbeddedSample("embeddedJet10", ueVariant);
            return "";
        }
        
        if (activeEmbeddedSample == SimSample::kEmbeddedPhoton10)
        {
            return InputSimEmbeddedSample("embeddedPhoton10", ueVariant);
        }
        if (activeEmbeddedSample == SimSample::kEmbeddedPhoton10And20Merged)
        {
            return MergedSimEmbeddedPath(
                                         CfgTagWithUEFor(kAA_JetPtMin, kAA_B2BCut, kAA_VzCut, kAA_IsoConeR, kAA_IsoMode, ueVariant),
                                         "photonJet10and20merged_SIM",
                                         "RecoilJets_embeddedPhoton10plus20_MERGED.root"
                                         );
        }
        return InputSimEmbedded(ueVariant);
    };
    
    const string outRoot = kOutputBase + "/auau/" + CfgTagAA();
    
    const vector<string> ueVariants = {"noSub", "baseVariant", "variantA", "variantB"};
    const vector<string> ueLabels   = {"No UE Sub", "with UE Sub", "Variant A", "Variant B"};
    
    const auto& centBins = CentBins();
    if (centBins.empty())
    {
        cout << ANSI_BOLD_YEL << "[WARN] No centrality bins defined — skipping UE comparisons\n"
        << ANSI_RESET;
        return;
    }
    
    TFile* fPP = TFile::Open(InputPP(isRun25pp).c_str(), "READ");
    if (!fPP || fPP->IsZombie())
    {
        if (fPP) { fPP->Close(); delete fPP; fPP = nullptr; }
        cout << ANSI_BOLD_YEL
        << "[WARN] Missing PP reference input for UE comparisons: "
        << InputPP(isRun25pp)
        << ANSI_RESET << "\n";
        return;
    }
    
    TDirectory* ppTop = fPP->GetDirectory(kTriggerPP.c_str());
    if (!ppTop) ppTop = fPP;
    
    // pp 5+10+20 merged SIM file (always open, independent of sim toggles)
    TFile* fPPSim510_20 = nullptr;
    TDirectory* ppSimTop = nullptr;
    {
        const string ppSimPath = SimInputPathForSample(SimSample::kPhotonJet5And10And20Merged);
        if (!ppSimPath.empty())
        {
            fPPSim510_20 = TFile::Open(ppSimPath.c_str(), "READ");
            if (fPPSim510_20 && !fPPSim510_20->IsZombie())
                ppSimTop = fPPSim510_20->GetDirectory(kDirSIM.c_str());
            else { if (fPPSim510_20) { fPPSim510_20->Close(); delete fPPSim510_20; fPPSim510_20 = nullptr; } }
        }
    }
    
    auto CountEmbeddedSelection = [](bool a, bool b, bool c) -> int
    {
        return (a ? 1 : 0) + (b ? 1 : 0) + (c ? 1 : 0);
    };
    
    const int nInclusiveEmbeddedSel =
    CountEmbeddedSelection(isInclusiveJet10Embedded, isInclusiveJet20Embedded, bothInclusiveJet10and20simEmbedded);
    const int nPhotonEmbeddedSel =
    CountEmbeddedSelection(isPhotonJet10Embedded, isPhotonJet20Embedded, bothPhoton10and20simEmbedded);
    
    auto InclusiveEmbeddedShortLabel = [&]() -> string
    {
        if (bothInclusiveJet10and20simEmbedded) return "(10+20)";
        if (isInclusiveJet10Embedded) return "10";
        if (isInclusiveJet20Embedded) return "20";
        return "";
    };
    
    auto PhotonEmbeddedShortLabel = [&]() -> string
    {
        if (bothPhoton10and20simEmbedded) return "(10+20)";
        if (isPhotonJet10Embedded) return "10";
        if (isPhotonJet20Embedded) return "20";
        return "";
    };
    
    const string inclusiveEmbeddedTitleTag =
    InclusiveEmbeddedShortLabel().empty() ? "" : ("inclusive jet " + InclusiveEmbeddedShortLabel());
    const string photonEmbeddedTitleTag =
    PhotonEmbeddedShortLabel().empty() ? "" : ("photon+jet " + PhotonEmbeddedShortLabel());
    
    auto ResolveInclusiveJetVariantInput = [&](const string& ueVariant) -> string
    {
        if (nInclusiveEmbeddedSel != 1)
        {
            return "";
        }
        
        if (bothInclusiveJet10and20simEmbedded)
        {
            const string cfgTag = CfgTagWithUEFor(
                                                  kAA_JetPtMin, kAA_B2BCut, kAA_VzCut, kAA_IsoConeR, kAA_IsoMode, ueVariant
                                                  );
            const string mergedPath = MergedSimEmbeddedPath(
                                                            cfgTag,
                                                            "embeddedJet10and20merged_SIM",
                                                            "RecoilJets_embeddedJet10plus20_MERGED.root"
                                                            );
            
            if (doPhotonJetMerge)
            {
                const string in10 = InputInclusiveJetEmbeddedSample("embeddedJet10", ueVariant);
                const string in20 = InputInclusiveJetEmbeddedSample("embeddedJet20", ueVariant);
                
                if (!gSystem->AccessPathName(in10.c_str()) &&
                    !gSystem->AccessPathName(in20.c_str()))
                {
                    const bool okMerge = BuildMergedSIMFile_PhotonSlices(
                                                                         {in10, in20},
                                                                         {kSigmaPhoton10_pb, kSigmaPhoton20_pb},
                                                                         mergedPath,
                                                                         kDirSIM,
                                                                         {"embeddedJet10", "embeddedJet20"}
                                                                         );
                    if (!okMerge) return "";
                }
                else
                {
                    return "";
                }
            }
            
            if (gSystem->AccessPathName(mergedPath.c_str())) return "";
            return mergedPath;
        }
        
        if (isInclusiveJet10Embedded)
        {
            const string in10 = InputInclusiveJetEmbeddedSample("embeddedJet10", ueVariant);
            if (gSystem->AccessPathName(in10.c_str())) return "";
            return in10;
        }
        
        if (isInclusiveJet20Embedded)
        {
            const string in20 = InputInclusiveJetEmbeddedSample("embeddedJet20", ueVariant);
            if (gSystem->AccessPathName(in20.c_str())) return "";
            return in20;
        }
        
        return "";
    };
    
    auto ResolvePhotonJetVariantInput = [&](const string& ueVariant) -> string
    {
        if (nPhotonEmbeddedSel != 1)
        {
            return "";
        }
        
        if (bothPhoton10and20simEmbedded)
        {
            const string cfgTag = CfgTagWithUEFor(
                                                  kAA_JetPtMin, kAA_B2BCut, kAA_VzCut, kAA_IsoConeR, kAA_IsoMode, ueVariant
                                                  );
            const string mergedPath = MergedSimEmbeddedPath(
                                                            cfgTag,
                                                            "photonJet10and20merged_SIM",
                                                            "RecoilJets_embeddedPhoton10plus20_MERGED.root"
                                                            );
            
            if (doPhotonJetMerge)
            {
                const string in10 = InputSimEmbeddedSample("embeddedPhoton10", ueVariant);
                const string in20 = InputSimEmbeddedSample("embeddedPhoton20", ueVariant);
                
                if (!gSystem->AccessPathName(in10.c_str()) &&
                    !gSystem->AccessPathName(in20.c_str()))
                {
                    const bool okMerge = BuildMergedSIMFile_PhotonSlices(
                                                                         {in10, in20},
                                                                         {kSigmaPhoton10_pb, kSigmaPhoton20_pb},
                                                                         mergedPath,
                                                                         kDirSIM,
                                                                         {"embeddedPhoton10", "embeddedPhoton20"}
                                                                         );
                    if (!okMerge) return "";
                }
                else
                {
                    return "";
                }
            }
            
            if (gSystem->AccessPathName(mergedPath.c_str())) return "";
            return mergedPath;
        }
        
        if (isPhotonJet10Embedded)
        {
            const string in10 = InputSimEmbeddedSample("embeddedPhoton10", ueVariant);
            if (gSystem->AccessPathName(in10.c_str())) return "";
            return in10;
        }
        
        if (isPhotonJet20Embedded)
        {
            const string in20 = InputSimEmbeddedSample("embeddedPhoton20", ueVariant);
            if (gSystem->AccessPathName(in20.c_str())) return "";
            return in20;
        }
        
        return "";
    };
    
    const auto& trigLoop = kTriggersAuAu;
    for (const auto& trigAA : trigLoop)
    {
        const string trigOutBase = JoinPath(outRoot, trigAA);
        const string sourceTopDirName = forEmbeddedSim ? kDirSIM : trigAA;
        const string ueCompModeBase = forEmbeddedSim
        ? JoinPath(trigOutBase, "isoQA/UEcomparisons_" + activeEmbeddedSimFolder)
        : JoinPath(trigOutBase, "isoQA/UEcomparisons");
        
        auto TriggerDisplayLabel = [&]() -> string
        {
            if (trigAA.find("photon_10") != string::npos) return "Photon 10 GeV + MBD NS #geq 2, vtx < 150 cm";
            if (trigAA.find("photon_12") != string::npos) return "Photon 12 GeV + MBD NS #geq 2, vtx < 150 cm";
            return "MBD NS #geq 2, vtx < 150 cm";
        };
        const string trigDisplayLabel = TriggerDisplayLabel();
        
        struct VariantHandle
        {
            string variant;
            string label;
            TFile* file = nullptr;
        };
        
        vector<VariantHandle> handles;
        handles.reserve(ueVariants.size());
        
        for (std::size_t iv = 0; iv < ueVariants.size(); ++iv)
        {
            VariantHandle H;
            H.variant = ueVariants[iv];
            H.label   = ueLabels[iv];
            const string varInput = forEmbeddedSim
            ? EmbeddedVariantInput(H.variant) : InputAuAu(H.variant);
            H.file    = TFile::Open(varInput.c_str(), "READ");
            
            if (!H.file || H.file->IsZombie())
            {
                if (H.file) { H.file->Close(); delete H.file; H.file = nullptr; }
                cout << ANSI_BOLD_YEL
                << "[WARN] Missing " << (forEmbeddedSim ? "embedded SIM" : "AuAu")
                << " UE variant input: "
                << varInput
                << ANSI_RESET << "\n";
            }
            
            handles.push_back(std::move(H));
        }
        
        
        const string ueCompBase = ueCompModeBase;
        const string centralitySummaryBase = JoinPath(ueCompBase, "centralitySummaryPerPt");
        const string meanIsoSummaryDir = JoinPath(centralitySummaryBase, "summaryOutput");
        const string ptSummaryBase = JoinPath(ueCompBase, "pTsummaryPerCentrality");
        const string ptMeanIsoSummaryDir = JoinPath(ptSummaryBase, "summaryOutput");
        const string perVariantOverlayBase = JoinPath(ueCompBase, "perVariantOverlays");
        EnsureDir(ueCompBase);
        EnsureDir(centralitySummaryBase);
        EnsureDir(meanIsoSummaryDir);
        EnsureDir(ptSummaryBase);
        EnsureDir(ptMeanIsoSummaryDir);
        EnsureDir(perVariantOverlayBase);
        
        if ((skipToCentralityAndPtOverlaysWithSSQA || SSoverlayPerVAR_processONLY) && !generateUEcomparisonSSQA)
        {
            cerr << ANSI_BOLD_RED
            << "[FATAL] skipToCentralityAndPtOverlaysWithSSQA / SSoverlayPerVAR_processONLY requires generateUEcomparisonSSQA = true"
            << ANSI_RESET << "\n";
            std::exit(1);
        }
        
        if (generateUEcomparisonSSQA && (skipToCentralityAndPtOverlaysWithSSQA || SSoverlayPerVAR_processONLY)) {
#include "AnalyzeRecoilJets_SSQA.cpp"
        }
        
        // ── Gaussian fit accumulators [ivH][ipt][ic] ──
        const std::size_t nVarSlots = handles.size();
        vector<vector<vector<double>>> gaussMean(nVarSlots, vector<vector<double>>(kNPtBins, vector<double>(centBins.size(), 0.0)));
        vector<vector<vector<double>>> gaussSigma(nVarSlots, vector<vector<double>>(kNPtBins, vector<double>(centBins.size(), 0.0)));
        vector<vector<vector<double>>> gaussMeanErr(nVarSlots, vector<vector<double>>(kNPtBins, vector<double>(centBins.size(), 0.0)));
        vector<vector<vector<double>>> gaussSigmaErr(nVarSlots, vector<vector<double>>(kNPtBins, vector<double>(centBins.size(), 0.0)));
        vector<vector<vector<bool>>>   gaussFilled(nVarSlots, vector<vector<bool>>(kNPtBins, vector<bool>(centBins.size(), false)));
        
        // Restricted-range Gaussian around the mode (ALICE/ATLAS isolation approach):
        // seed from peak bin, fit asymmetrically weighting the left (signal) side,
        // truncating the right-side combinatoric tail.
        auto FitGaussianIterative = [](TH1* h, double& outMean, double& outSigma,
                                       double& outMeanErr, double& outSigmaErr) -> bool
        {
            if (!h || h->GetEntries() < 30) return false;
            // Seed from mode (peak bin), not arithmetic mean
            const int maxBin = h->GetMaximumBin();
            double mu  = h->GetBinCenter(maxBin);
            // Initial σ from left-side half-width at half-max
            const double halfMax = 0.5 * h->GetBinContent(maxBin);
            double sig = 0.0;
            for (int ib = maxBin - 1; ib >= 1; --ib)
            {
                if (h->GetBinContent(ib) <= halfMax)
                {
                    sig = std::abs(mu - h->GetBinCenter(ib));
                    break;
                }
            }
            if (sig <= 0.5) sig = h->GetRMS() * 0.5;
            // 4 narrowing passes: asymmetric range (wider left, tighter right)
            for (int pass = 0; pass < 4; ++pass)
            {
                const double nL = (pass < 2) ? 2.0 : 2.5;
                const double nR = (pass < 2) ? 1.2 : 1.5;
                const double lo = mu - nL * sig;
                const double hi = mu + nR * sig;
                TF1 fg("fg_iter", "gaus", lo, hi);
                fg.SetParameters(h->GetBinContent(maxBin), mu, sig);
                fg.SetParLimits(2, 0.01, 10.0 * std::abs(sig));
                h->Fit(&fg, "QNR0", "", lo, hi);
                mu  = fg.GetParameter(1);
                sig = std::abs(fg.GetParameter(2));
                if (sig < 0.01) sig = 0.5;
            }
            // Final fit: asymmetric range centered on converged peak
            const double loF = mu - 2.5 * sig;
            const double hiF = mu + 1.5 * sig;
            TF1 gf("gf_final", "gaus", loF, hiF);
            gf.SetParameters(h->GetBinContent(h->FindBin(mu)), mu, sig);
            gf.SetParLimits(2, 0.01, 10.0 * std::abs(sig));
            const int st = h->Fit(&gf, "QNR0", "", loF, hiF);
            if (st != 0 && st != 4000) return false;
            outMean     = gf.GetParameter(1);
            outSigma    = std::abs(gf.GetParameter(2));
            outMeanErr  = gf.GetParError(1);
            outSigmaErr = gf.GetParError(2);
            return true;
        };
        
        // Draw a Gaussian fit curve on the current pad matching a histogram's color
        auto DrawGaussFitCurve = [&](TH1* h, int color) -> TF1*
        {
            if (!h || h->GetEntries() < 30) return nullptr;
            double mu, sig, muE, sigE;
            TH1* hTmp = (TH1*)h->Clone("hTmp_gaussFit");
            hTmp->SetDirectory(nullptr);
            bool ok = FitGaussianIterative(hTmp, mu, sig, muE, sigE);
            delete hTmp;
            if (!ok) return nullptr;
            TF1* fDraw = new TF1(TString::Format("fGauss_%p", (void*)h).Data(),
                                 "gaus", mu - 3.0 * sig, mu + 3.0 * sig);
            fDraw->SetParameters(h->GetMaximum() * 0.95, mu, sig);
            // scale amplitude to match histogram bin width
            const double binW = h->GetBinWidth(1);
            const double area = h->Integral(h->FindBin(mu - 3.0*sig), h->FindBin(mu + 3.0*sig), "width");
            if (area > 0.0) fDraw->SetParameter(0, area / (sig * std::sqrt(2.0 * TMath::Pi())));
            else fDraw->SetParameter(0, h->GetMaximum());
            fDraw->SetParameter(1, mu);
            fDraw->SetParameter(2, sig);
            fDraw->SetLineColor(color);
            fDraw->SetLineStyle(2);
            fDraw->SetLineWidth(2);
            fDraw->SetNpx(200);
            fDraw->Draw("SAME");
            return fDraw;
        };
        
        for (std::size_t ivH = 0; ivH < handles.size(); ++ivH)
        {
            auto& H = handles[ivH];
            if (!H.file) continue;
            
            const string variantDir = JoinPath(ueCompBase, H.variant);
            const string variantCentralitySummaryDir = JoinPath(centralitySummaryBase, H.variant);
            const string variantPtSummaryDir = JoinPath(ptSummaryBase, H.variant);
            const bool doVariantDetailPlots = !forEmbeddedSim;
            if (doVariantDetailPlots) EnsureDir(variantDir);
            EnsureDir(variantCentralitySummaryDir);
            EnsureDir(variantPtSummaryDir);
            
            TDirectory* aaTop = H.file->GetDirectory(sourceTopDirName.c_str());
            if (!aaTop) continue;
            
            // -- open MC files for inclusive Eiso overlay (DATA mode only) --
            TFile* fIncMCvar = nullptr;
            TFile* fPhoMCvar = nullptr;
            TDirectory* incMCvarTop = nullptr;
            TDirectory* phoMCvarTop = nullptr;
            if (!forEmbeddedSim)
            {
                const string incIn = ResolveInclusiveJetVariantInput(H.variant);
                if (!incIn.empty())
                {
                    fIncMCvar = TFile::Open(incIn.c_str(), "READ");
                    if (fIncMCvar && !fIncMCvar->IsZombie())
                        incMCvarTop = fIncMCvar->GetDirectory(kDirSIM.c_str());
                    else { if (fIncMCvar) { fIncMCvar->Close(); delete fIncMCvar; fIncMCvar = nullptr; } }
                }
                const string phoIn = ResolvePhotonJetVariantInput(H.variant);
                if (!phoIn.empty())
                {
                    fPhoMCvar = TFile::Open(phoIn.c_str(), "READ");
                    if (fPhoMCvar && !fPhoMCvar->IsZombie())
                        phoMCvarTop = fPhoMCvar->GetDirectory(kDirSIM.c_str());
                    else { if (fPhoMCvar) { fPhoMCvar->Close(); delete fPhoMCvar; fPhoMCvar = nullptr; } }
                }
            }
            
            // -- accumulators for <E_T^iso> vs centrality (per pT bin) --
            vector<vector<double>> vsCent_yPP(kNPtBins);
            vector<vector<double>> vsCent_eyPP(kNPtBins);
            vector<vector<double>> vsCent_yAA(kNPtBins);
            vector<vector<double>> vsCent_eyAA(kNPtBins);
            vector<vector<bool>>   vsCent_filled(kNPtBins, vector<bool>(centBins.size(), false));
            for (int ip = 0; ip < kNPtBins; ++ip)
            {
                vsCent_yPP[ip].resize(centBins.size(), 0.0);
                vsCent_eyPP[ip].resize(centBins.size(), 0.0);
                vsCent_yAA[ip].resize(centBins.size(), 0.0);
                vsCent_eyAA[ip].resize(centBins.size(), 0.0);
            }
            
            // MC mean iso accumulators [ipt][ic] for vs-centrality and vs-pT summary plots
            vector<vector<double>> vsCent_yIncMC(kNPtBins, vector<double>(centBins.size(), 0.0));
            vector<vector<double>> vsCent_eyIncMC(kNPtBins, vector<double>(centBins.size(), 0.0));
            vector<vector<double>> vsCent_yPhoMC(kNPtBins, vector<double>(centBins.size(), 0.0));
            vector<vector<double>> vsCent_eyPhoMC(kNPtBins, vector<double>(centBins.size(), 0.0));
            vector<vector<bool>>   vsCent_filledIncMC(kNPtBins, vector<bool>(centBins.size(), false));
            vector<vector<bool>>   vsCent_filledPhoMC(kNPtBins, vector<bool>(centBins.size(), false));
            
            // Signal (tight data + truth-sig-matched MC) mean iso accumulators
            vector<vector<double>> vsCent_yTightData(kNPtBins, vector<double>(centBins.size(), 0.0));
            vector<vector<double>> vsCent_eyTightData(kNPtBins, vector<double>(centBins.size(), 0.0));
            vector<vector<bool>>   vsCent_filledTightData(kNPtBins, vector<bool>(centBins.size(), false));
            vector<vector<double>> vsCent_ySigIncMC(kNPtBins, vector<double>(centBins.size(), 0.0));
            vector<vector<double>> vsCent_eySigIncMC(kNPtBins, vector<double>(centBins.size(), 0.0));
            vector<vector<bool>>   vsCent_filledSigIncMC(kNPtBins, vector<bool>(centBins.size(), false));
            vector<vector<double>> vsCent_ySigPhoMC(kNPtBins, vector<double>(centBins.size(), 0.0));
            vector<vector<double>> vsCent_eySigPhoMC(kNPtBins, vector<double>(centBins.size(), 0.0));
            vector<vector<bool>>   vsCent_filledSigPhoMC(kNPtBins, vector<bool>(centBins.size(), false));
            
            for (std::size_t ic = 0; ic < centBins.size(); ++ic)
            {
                const auto& cb = centBins[ic];
                const string centDir = JoinPath(variantDir, cb.folder);
                if (doVariantDetailPlots) EnsureDir(centDir);
                
                vector<double> xPt;
                vector<double> exPt;
                vector<double> yPP;
                vector<double> eyPP;
                vector<double> yAA;
                vector<double> eyAA;
                
                xPt.reserve(kNPtBins);
                exPt.reserve(kNPtBins);
                yPP.reserve(kNPtBins);
                eyPP.reserve(kNPtBins);
                yAA.reserve(kNPtBins);
                eyAA.reserve(kNPtBins);
                
                for (int ipt = 0; ipt < kNPtBins; ++ipt)
                {
                    const PtBin& b = PtBins()[ipt];
                    const string ptDir = JoinPath(centDir, b.folder);
                    if (doVariantDetailPlots) EnsureDir(ptDir);
                    
                    const string hPPName = "h_Eiso" + b.suffix;
                    const string hAAName = "h_Eiso" + b.suffix + cb.suffix;
                    
                    TH1* hPPsrc = dynamic_cast<TH1*>(ppTop->Get(hPPName.c_str()));
                    TH1* hAAsrc = dynamic_cast<TH1*>(aaTop->Get(hAAName.c_str()));
                    if (!hPPsrc || !hAAsrc) continue;
                    
                    TH1* hPP = CloneTH1(
                                        hPPsrc,
                                        TString::Format("hPP_isoOverlay_%s_%s_%s",
                                                        trigAA.c_str(), cb.folder.c_str(), b.folder.c_str()).Data()
                                        );
                    TH1* hAA = CloneTH1(
                                        hAAsrc,
                                        TString::Format("hAA_isoOverlay_%s_%s_%s_%s",
                                                        trigAA.c_str(), H.variant.c_str(), cb.folder.c_str(), b.folder.c_str()).Data()
                                        );
                    
                    if (!hPP || !hAA)
                    {
                        if (hPP) delete hPP;
                        if (hAA) delete hAA;
                        continue;
                    }
                    
                    EnsureSumw2(hPP);
                    EnsureSumw2(hAA);
                    
                    xPt.push_back(0.5 * (kPtEdges[(std::size_t)ipt] + kPtEdges[(std::size_t)ipt + 1]));
                    exPt.push_back(0.5 * (kPtEdges[(std::size_t)ipt + 1] - kPtEdges[(std::size_t)ipt]));
                    yPP.push_back(hPP->GetMean());
                    yAA.push_back(hAA->GetMean());
                    eyPP.push_back((hPP->GetEntries() > 0.0) ? (hPP->GetRMS() / std::sqrt(hPP->GetEntries())) : 0.0);
                    eyAA.push_back((hAA->GetEntries() > 0.0) ? (hAA->GetRMS() / std::sqrt(hAA->GetEntries())) : 0.0);
                    
                    // accumulate for vs-centrality plot
                    vsCent_yPP[ipt][ic]  = hPP->GetMean();
                    vsCent_eyPP[ipt][ic] = (hPP->GetEntries() > 0.0) ? (hPP->GetRMS() / std::sqrt(hPP->GetEntries())) : 0.0;
                    vsCent_yAA[ipt][ic]  = hAA->GetMean();
                    vsCent_eyAA[ipt][ic] = (hAA->GetEntries() > 0.0) ? (hAA->GetRMS() / std::sqrt(hAA->GetEntries())) : 0.0;
                    vsCent_filled[ipt][ic] = true;
                    
                    // Gaussian fit on raw AuAu histogram (rebin for stability)
                    {
                        TH1* hFitTmp = CloneTH1(hAAsrc,
                            TString::Format("hGaussFit_%s_%s_%s_%zu",
                                H.variant.c_str(), cb.folder.c_str(), b.folder.c_str(), ivH).Data());
                        if (hFitTmp)
                        {
                            hFitTmp->Rebin(10);
                            EnsureSumw2(hFitTmp);
                            double gM, gS, gME, gSE;
                            if (FitGaussianIterative(hFitTmp, gM, gS, gME, gSE))
                            {
                                gaussMean[ivH][ipt][ic]     = gM;
                                gaussSigma[ivH][ipt][ic]    = gS;
                                gaussMeanErr[ivH][ipt][ic]  = gME;
                                gaussSigmaErr[ivH][ipt][ic] = gSE;
                                gaussFilled[ivH][ipt][ic]   = true;
                            }
                            delete hFitTmp;
                        }
                    }
                    
                    // accumulate MC embedded means for summary plots
                    if (incMCvarTop)
                    {
                        TH1* hIncMC = dynamic_cast<TH1*>(incMCvarTop->Get(hAAName.c_str()));
                        if (hIncMC)
                        {
                            vsCent_yIncMC[ipt][ic]  = hIncMC->GetMean();
                            vsCent_eyIncMC[ipt][ic] = (hIncMC->GetEntries() > 0.0)
                            ? (hIncMC->GetRMS() / std::sqrt(hIncMC->GetEntries())) : 0.0;
                            vsCent_filledIncMC[ipt][ic] = true;
                        }
                    }
                    if (phoMCvarTop)
                    {
                        TH1* hPhoMC = dynamic_cast<TH1*>(phoMCvarTop->Get(hAAName.c_str()));
                        if (hPhoMC)
                        {
                            vsCent_yPhoMC[ipt][ic]  = hPhoMC->GetMean();
                            vsCent_eyPhoMC[ipt][ic] = (hPhoMC->GetEntries() > 0.0)
                            ? (hPhoMC->GetRMS() / std::sqrt(hPhoMC->GetEntries())) : 0.0;
                            vsCent_filledPhoMC[ipt][ic] = true;
                        }
                    }
                    
                    // accumulate signal (tight data + truth-sig-matched MC) means
                    {
                        const string hTightName = "h_Eiso_tight" + b.suffix + cb.suffix;
                        TH1* hTightSrc = dynamic_cast<TH1*>(aaTop->Get(hTightName.c_str()));
                        if (hTightSrc && hTightSrc->GetEntries() > 0.0)
                        {
                            vsCent_yTightData[ipt][ic]  = hTightSrc->GetMean();
                            vsCent_eyTightData[ipt][ic] = hTightSrc->GetMeanError();
                            vsCent_filledTightData[ipt][ic] = true;
                        }
                        const string hSigMCName = "h_EisoReco_truthSigMatched" + b.suffix + cb.suffix;
                        if (incMCvarTop)
                        {
                            TH1* hSigInc = dynamic_cast<TH1*>(incMCvarTop->Get(hSigMCName.c_str()));
                            if (hSigInc && hSigInc->GetEntries() > 0.0)
                            {
                                vsCent_ySigIncMC[ipt][ic]  = hSigInc->GetMean();
                                vsCent_eySigIncMC[ipt][ic] = hSigInc->GetMeanError();
                                vsCent_filledSigIncMC[ipt][ic] = true;
                            }
                        }
                        if (phoMCvarTop)
                        {
                            TH1* hSigPho = dynamic_cast<TH1*>(phoMCvarTop->Get(hSigMCName.c_str()));
                            if (hSigPho && hSigPho->GetEntries() > 0.0)
                            {
                                vsCent_ySigPhoMC[ipt][ic]  = hSigPho->GetMean();
                                vsCent_eySigPhoMC[ipt][ic] = hSigPho->GetMeanError();
                                vsCent_filledSigPhoMC[ipt][ic] = true;
                            }
                        }
                    }
                    
                    if (generateISOpTcentOverlaysONLY)
                    {
                        delete hPP;
                        delete hAA;
                        continue;  // skip detail plots, keep accumulation
                    }
                    
                    const double intPP = hPP->Integral(0, hPP->GetNbinsX() + 1);
                    const double intAA = hAA->Integral(0, hAA->GetNbinsX() + 1);
                    if (intPP > 0.0) hPP->Scale(1.0 / intPP);
                    if (intAA > 0.0) hAA->Scale(1.0 / intAA);
                    
                    hPP->SetLineColor(kRed + 1);
                    hPP->SetMarkerColor(kRed + 1);
                    hPP->SetMarkerStyle(24);
                    hPP->SetMarkerSize(1.1);
                    hPP->SetLineWidth(2);
                    hPP->SetFillStyle(0);
                    
                    hAA->SetLineColor(kBlack);
                    hAA->SetMarkerColor(kBlack);
                    hAA->SetMarkerStyle(20);
                    hAA->SetMarkerSize(1.1);
                    hAA->SetLineWidth(2);
                    hAA->SetFillStyle(0);
                    
                    const double ymax = std::max(hPP->GetMaximum(), hAA->GetMaximum());
                    
                    TCanvas c(
                              TString::Format("c_isoOverlayPPAuAu_%s_%s_%s_%s",
                                              trigAA.c_str(), H.variant.c_str(), cb.folder.c_str(), b.folder.c_str()).Data(),
                              "c_isoOverlayPPAuAu", 900, 700
                              );
                    ApplyCanvasMargins1D(c);
                    c.cd();
                    
                    hAA->SetTitle("");
                    hAA->GetXaxis()->SetTitle("E_{T}^{iso} [GeV]");
                    hAA->GetYaxis()->SetTitle("Normalized to unit area");
                    hAA->GetXaxis()->SetTitleSize(0.055);
                    hAA->GetYaxis()->SetTitleSize(0.055);
                    hAA->GetXaxis()->SetLabelSize(0.045);
                    hAA->GetYaxis()->SetLabelSize(0.045);
                    hAA->GetYaxis()->SetTitleOffset(1.15);
                    hAA->SetMinimum(0.0);
                    hAA->SetMaximum((ymax > 0.0) ? (1.25 * ymax) : 1.0);
                    
                    hAA->Draw("E1");
                    hPP->Draw("E1 SAME");
                    
                    TLegend leg(0.56, 0.68, 0.92, 0.88);
                    leg.SetBorderSize(0);
                    leg.SetFillStyle(0);
                    leg.SetTextFont(42);
                    leg.SetTextSize(0.032);
                    leg.AddEntry(hPP, "pp data", "ep");
                    leg.AddEntry(hAA, TString::Format("AuAu data (%s)", H.label.c_str()).Data(), "ep");
                    leg.Draw();
                    
                    TLatex tTitle;
                    tTitle.SetNDC(true);
                    tTitle.SetTextFont(42);
                    tTitle.SetTextAlign(23);
                    tTitle.SetTextSize(0.045);
                    tTitle.DrawLatex(0.50, 0.955,
                                     TString::Format("isolation overlay pp and %d-%d%% Cent AuAu", cb.lo, cb.hi).Data());
                    
                    TLatex t;
                    t.SetNDC(true);
                    t.SetTextFont(42);
                    t.SetTextAlign(13);
                    t.SetTextSize(0.034);
                    t.DrawLatex(0.16, 0.88, TString::Format("Trigger: %s", trigAA.c_str()).Data());
                    t.DrawLatex(0.16, 0.84, TString::Format("|v_{z}| < %d cm", kAA_VzCut).Data());
                    t.DrawLatex(0.16, 0.80, TString::Format("UE subtraction: %s", H.label.c_str()).Data());
                    t.DrawLatex(0.16, 0.76, TString::Format("p_{T}^{#gamma}: %d-%d GeV", b.lo, b.hi).Data());
                    
                    {
                        const string ppAuAuDir = JoinPath(ptDir, "ppAuAuOverlay");
                        EnsureDir(ppAuAuDir);
                        const string furRealDir = JoinPath(ppAuAuDir, "ppAuAuOverlayFURREAL");
                        EnsureDir(furRealDir);
                        if (doVariantDetailPlots) SaveCanvas(c, JoinPath(furRealDir, "isolationOverlay_pp_vs_auau.png"));
                        
                        // ── pp data vs pp SIM (5+10+20) inclusive overlay ──
                        if (ppSimTop)
                        {
                            const string ppIncDir = JoinPath(ppAuAuDir, "inclusiveIsoOverlays");
                            EnsureDir(ppIncDir);
                            const string hPPincName = "h_Eiso" + b.suffix;
                            TH1* hPPdSrc = dynamic_cast<TH1*>(ppTop->Get(hPPincName.c_str()));
                            TH1* hPPmSrc = dynamic_cast<TH1*>(ppSimTop->Get(hPPincName.c_str()));
                            if (hPPdSrc && hPPmSrc)
                            {
                                TH1* hPPd = CloneTH1(hPPdSrc, TString::Format("hPPd_inc_%s_%s", cb.folder.c_str(), b.folder.c_str()).Data());
                                TH1* hPPm = CloneTH1(hPPmSrc, TString::Format("hPPm_inc_%s_%s", cb.folder.c_str(), b.folder.c_str()).Data());
                                if (hPPd && hPPm)
                                {
                                    EnsureSumw2(hPPd); EnsureSumw2(hPPm);
                                    hPPd->Rebin(10); hPPm->Rebin(10);
                                    const double iD = hPPd->Integral(0, hPPd->GetNbinsX()+1);
                                    const double iM = hPPm->Integral(0, hPPm->GetNbinsX()+1);
                                    if (iD > 0.0 && iM > 0.0)
                                    {
                                        hPPd->Scale(1.0/iD); hPPm->Scale(1.0/iM);
                                        for (int ib=0; ib<=hPPm->GetNbinsX()+1; ++ib) hPPm->SetBinError(ib,0.0);
                                        hPPm->SetTitle(""); hPPm->SetLineColor(kBlack); hPPm->SetLineWidth(2);
                                        hPPm->SetFillStyle(0); hPPm->SetMarkerSize(0.0);
                                        hPPd->SetLineColor(kBlue+1); hPPd->SetMarkerColor(kBlue+1);
                                        hPPd->SetMarkerStyle(20); hPPd->SetMarkerSize(1.0); hPPd->SetLineWidth(2); hPPd->SetFillStyle(0);
                                        const double ym = std::max(hPPd->GetMaximum(), hPPm->GetMaximum());
                                        TCanvas cPI(TString::Format("c_ppInc_%s_%s", cb.folder.c_str(), b.folder.c_str()).Data(), "c_ppInc", 900, 700);
                                        ApplyCanvasMargins1D(cPI); cPI.cd();
                                        hPPm->GetXaxis()->SetTitle("E_{T}^{iso} [GeV]");
                                        hPPm->GetYaxis()->SetTitle("Normalized to unit area");
                                        hPPm->GetXaxis()->SetTitleSize(0.055); hPPm->GetYaxis()->SetTitleSize(0.055);
                                        hPPm->GetXaxis()->SetLabelSize(0.045); hPPm->GetYaxis()->SetLabelSize(0.045);
                                        hPPm->GetYaxis()->SetTitleOffset(1.15);
                                        hPPm->SetMinimum(0.0); hPPm->SetMaximum((ym>0)?(1.25*ym):1.0);
                                        hPPm->Draw("hist"); hPPd->Draw("E1 SAME");
                                        TLegend lg(0.56,0.72,0.92,0.88); lg.SetBorderSize(0); lg.SetFillStyle(0); lg.SetTextFont(42); lg.SetTextSize(0.032);
                                        lg.AddEntry(hPPd, "pp data", "ep"); lg.AddEntry(hPPm, "pp MC (5+10+20)", "l"); lg.Draw();
                                        TLatex ts; ts.SetNDC(true); ts.SetTextFont(42); ts.SetTextAlign(33); ts.SetTextSize(0.048);
                                        ts.DrawLatex(0.92,0.6,"#bf{sPHENIX} #it{Internal}"); ts.SetTextSize(0.038);
                                        ts.DrawLatex(0.92,0.54,"p+p  #sqrt{s} = 200 GeV");
                                        TLatex ti; ti.SetNDC(true); ti.SetTextFont(42); ti.SetTextAlign(13); ti.SetTextSize(0.045);
                                        ti.DrawLatex(0.18,0.89,TString::Format("p_{T}^{#gamma} = %d-%d GeV", b.lo, b.hi).Data());
                                        SaveCanvas(cPI, JoinPath(ppIncDir, "Eiso_dataMC_overlay.png"));
                                    }
                                    delete hPPd; delete hPPm;
                                }
                                else { if (hPPd) delete hPPd; if (hPPm) delete hPPm; }
                            }
                            
                            // ── pp signal/bkg decomposition ──
                            const string ppSigDir = JoinPath(ppAuAuDir, "signalIsoOverlays");
                            EnsureDir(ppSigDir);
                            const string hTName = "h_Eiso_tight" + b.suffix;
                            const string hNTName = "h_Eiso_nonTight" + b.suffix;
                            const string hSigName = "h_EisoReco_truthSigMatched" + b.suffix;
                            TH1* hTdSrc  = dynamic_cast<TH1*>(ppTop->Get(hTName.c_str()));
                            TH1* hNTdSrc = dynamic_cast<TH1*>(ppTop->Get(hNTName.c_str()));
                            TH1* hSigSrc = dynamic_cast<TH1*>(ppSimTop->Get(hSigName.c_str()));
                            if (hTdSrc && hNTdSrc && hSigSrc)
                            {
                                TH1* hTd  = CloneTH1(hTdSrc,  TString::Format("hTd_ppSig_%s_%s", cb.folder.c_str(), b.folder.c_str()).Data());
                                TH1* hNTd = CloneTH1(hNTdSrc, TString::Format("hNTd_ppSig_%s_%s", cb.folder.c_str(), b.folder.c_str()).Data());
                                TH1* hSig = CloneTH1(hSigSrc, TString::Format("hSig_ppSig_%s_%s", cb.folder.c_str(), b.folder.c_str()).Data());
                                if (hTd && hNTd && hSig)
                                {
                                    hTd->Rebin(10); hNTd->Rebin(10); hSig->Rebin(10);
                                    {
                                        const double iTdN=hTd->Integral(0,hTd->GetNbinsX()+1); const double iNTdN=hNTd->Integral(0,hNTd->GetNbinsX()+1); const double iSigN=hSig->Integral(0,hSig->GetNbinsX()+1);
                                        if (iTdN>0) hTd->Scale(1.0/iTdN); if (iNTdN>0) hNTd->Scale(1.0/iNTdN); if (iSigN>0) hSig->Scale(1.0/iSigN);
                                        hSig->SetLineColor(kBlue-9); hSig->SetFillColorAlpha(kBlue-9,0.5); hSig->SetFillStyle(1001); hSig->SetLineWidth(1); hSig->SetMarkerSize(0.0);
                                        hNTd->SetLineColor(kPink-4); hNTd->SetFillColorAlpha(kPink-4,0.5); hNTd->SetFillStyle(1001); hNTd->SetLineWidth(1); hNTd->SetMarkerSize(0.0);
                                        for (int ib=0; ib<=hSig->GetNbinsX()+1; ++ib) hSig->SetBinError(ib,0.0);
                                        for (int ib=0; ib<=hNTd->GetNbinsX()+1; ++ib) hNTd->SetBinError(ib,0.0);
                                        hTd->SetLineColor(kBlack); hTd->SetMarkerColor(kBlack); hTd->SetMarkerStyle(20); hTd->SetMarkerSize(1.0); hTd->SetLineWidth(2); hTd->SetFillStyle(0);
                                        TH1* hStackPP = (TH1*)hSig->Clone(TString::Format("hStack_ppSig_%s_%s", cb.folder.c_str(), b.folder.c_str()).Data());
                                        hStackPP->SetDirectory(nullptr); hStackPP->Add(hNTd);
                                        const double ym = std::max(hTd->GetMaximum(), hStackPP->GetMaximum());
                                        TCanvas cSB(TString::Format("c_ppSigBkg_%s_%s", cb.folder.c_str(), b.folder.c_str()).Data(), "c_ppSigBkg", 900, 700);
                                        ApplyCanvasMargins1D(cSB); cSB.cd();
                                        hStackPP->SetTitle(""); hStackPP->GetXaxis()->SetTitle("E_{T}^{iso} [GeV]");
                                        hStackPP->GetYaxis()->SetTitle("Normalized to unit area");
                                        hTd->GetXaxis()->SetTitleSize(0.055); hTd->GetYaxis()->SetTitleSize(0.055);
                                        hTd->GetXaxis()->SetLabelSize(0.045); hTd->GetYaxis()->SetLabelSize(0.045);
                                        hTd->GetYaxis()->SetTitleOffset(1.15);
                                        hStackPP->GetXaxis()->SetTitleSize(0.055); hStackPP->GetYaxis()->SetTitleSize(0.055);
                                        hStackPP->GetXaxis()->SetLabelSize(0.045); hStackPP->GetYaxis()->SetLabelSize(0.045);
                                        hStackPP->GetYaxis()->SetTitleOffset(1.15);
                                        hStackPP->SetMinimum(0.0); hStackPP->SetMaximum((ym>0)?(1.3*ym):1.0);
                                        hStackPP->SetLineColor(kPink-4); hStackPP->SetFillColorAlpha(kPink-4,0.5); hStackPP->SetFillStyle(1001);
                                        hStackPP->Draw("hist");
                                        hSig->Draw("hist SAME");
                                        hTd->Draw("E1 SAME");
                                        TLatex tSBtitle; tSBtitle.SetNDC(true); tSBtitle.SetTextFont(42); tSBtitle.SetTextAlign(23); tSBtitle.SetTextSize(0.040);
                                        tSBtitle.DrawLatex(0.50, 0.97,
                                                           "E_{T}^{iso}, p+p data versus Photon+Jet (5+10+20) Pythia");
                                        TLegend lg(0.56,0.65,0.92,0.88); lg.SetBorderSize(0); lg.SetFillStyle(0); lg.SetTextFont(42); lg.SetTextSize(0.032);
                                        lg.AddEntry(hTd, "Data (Signal)", "ep"); lg.AddEntry(hNTd, "Data (Background)", "f"); lg.AddEntry(hSig, "Signal Embedded MC", "f"); lg.Draw();
                                        TLatex ts; ts.SetNDC(true); ts.SetTextFont(42); ts.SetTextAlign(33); ts.SetTextSize(0.048);
                                        ts.DrawLatex(0.92,0.6,"#bf{sPHENIX} #it{Internal}"); ts.SetTextSize(0.038);
                                        ts.DrawLatex(0.92,0.54,"p+p  #sqrt{s} = 200 GeV");
                                        TLatex tsPPinfo; tsPPinfo.SetNDC(true); tsPPinfo.SetTextFont(42); tsPPinfo.SetTextAlign(33); tsPPinfo.SetTextSize(0.034);
                                        tsPPinfo.DrawLatex(0.92,0.46,"Photon 4 GeV + MBD NS #geq 1");
                                        tsPPinfo.DrawLatex(0.92,0.42,TString::Format("|v_{z}| < %d cm", kAA_VzCut).Data());
                                        tsPPinfo.DrawLatex(0.92,0.38,TString::Format("#DeltaR_{cone} < %.1f", (kAA_IsoConeR == "isoR40") ? 0.4 : 0.3).Data());
                                        TLatex ti; ti.SetNDC(true); ti.SetTextFont(42); ti.SetTextAlign(13); ti.SetTextSize(0.045);
                                        ti.DrawLatex(0.18,0.89,TString::Format("p_{T}^{#gamma} = %d-%d GeV", b.lo, b.hi).Data());
                                        SaveCanvas(cSB, JoinPath(ppSigDir, "Eiso_sigBkg_overlay.png"));
                                        delete hStackPP;
                                    }
                                }
                                if (hTd) delete hTd; if (hNTd) delete hNTd; if (hSig) delete hSig;
                            }
                        }
                    }
                    
                    delete hPP;
                    delete hAA;
                    // ── inclusive Eiso: DATA vs MC overlays ──
                    if (!forEmbeddedSim)
                    {
                        struct MCOverlayCfg {
                            TDirectory*  mcTop;
                            string       folder;
                            string       titleTag;
                            string       mcLegend;
                        };
                        const vector<MCOverlayCfg> mcOvCfgs = {
                            { incMCvarTop, "inclusiveMCoverlays", inclusiveEmbeddedTitleTag, "inclusive MC" },
                            { phoMCvarTop, "photonJetOverlays",   photonEmbeddedTitleTag,    "photon+jet MC" }
                        };
                        
                        for (const auto& mcCfg : mcOvCfgs)
                        {
                            if (!mcCfg.mcTop || mcCfg.titleTag.empty()) continue;
                            
                            const string mcOvDir = JoinPath(ptDir, mcCfg.folder);
                            EnsureDir(mcOvDir);
                            
                            const string hAAName = "h_Eiso" + b.suffix + cb.suffix;
                            TH1* hDataSrc = dynamic_cast<TH1*>(aaTop->Get(hAAName.c_str()));
                            TH1* hMCSrc   = dynamic_cast<TH1*>(mcCfg.mcTop->Get(hAAName.c_str()));
                            if (!hDataSrc || !hMCSrc) continue;
                            
                            TH1* hData = CloneTH1(hDataSrc,
                                                  TString::Format("hData_mcOv_%s_%s_%s_%s",
                                                                  mcCfg.folder.c_str(), H.variant.c_str(), cb.folder.c_str(), b.folder.c_str()).Data());
                            TH1* hMC = CloneTH1(hMCSrc,
                                                TString::Format("hMC_mcOv_%s_%s_%s_%s",
                                                                mcCfg.folder.c_str(), H.variant.c_str(), cb.folder.c_str(), b.folder.c_str()).Data());
                            if (!hData || !hMC) { delete hData; delete hMC; continue; }
                            
                            EnsureSumw2(hData);
                            EnsureSumw2(hMC);
                            
                            hData->Rebin(10);
                            hMC->Rebin(10);
                            
                            const double intData = hData->Integral(0, hData->GetNbinsX() + 1);
                            const double intMC   = hMC->Integral(0, hMC->GetNbinsX() + 1);
                            if (!(intData > 0.0) || !(intMC > 0.0)) { delete hData; delete hMC; continue; }
                            
                            hData->Scale(1.0 / intData);
                            hMC->Scale(1.0 / intMC);
                            
                            hMC->SetTitle("");
                            hMC->SetLineColor(kBlack);
                            hMC->SetLineWidth(2);
                            hMC->SetFillStyle(0);
                            hMC->SetMarkerSize(0.0);
                            for (int ib = 0; ib <= hMC->GetNbinsX() + 1; ++ib) hMC->SetBinError(ib, 0.0);
                            
                            hData->SetLineColor(kBlue + 1);
                            hData->SetMarkerColor(kBlue + 1);
                            hData->SetMarkerStyle(20);
                            hData->SetMarkerSize(1.0);
                            hData->SetLineWidth(2);
                            hData->SetFillStyle(0);
                            
                            const double ymx = std::max(hData->GetMaximum(), hMC->GetMaximum());
                            
                            TCanvas cMcOv(
                                          TString::Format("c_mcOv_%s_%s_%s_%s",
                                                          mcCfg.folder.c_str(), H.variant.c_str(), cb.folder.c_str(), b.folder.c_str()).Data(),
                                          "c_mcOv", 900, 700);
                            ApplyCanvasMargins1D(cMcOv);
                            cMcOv.cd();
                            
                            hMC->GetXaxis()->SetTitle("E_{T}^{iso} [GeV]");
                            hMC->GetYaxis()->SetTitle("Normalized to unit area");
                            hMC->GetXaxis()->SetTitleSize(0.055);
                            hMC->GetYaxis()->SetTitleSize(0.055);
                            hMC->GetXaxis()->SetLabelSize(0.045);
                            hMC->GetYaxis()->SetLabelSize(0.045);
                            hMC->GetYaxis()->SetTitleOffset(1.15);
                            hMC->SetMinimum(0.0);
                            hMC->SetMaximum((ymx > 0.0) ? (1.25 * ymx) : 1.0);
                            
                            hMC->Draw("hist");
                            hData->Draw("E1 SAME");
                            
                            TLegend legMcOv(0.56, 0.72, 0.92, 0.88);
                            legMcOv.SetBorderSize(0);
                            legMcOv.SetFillStyle(0);
                            legMcOv.SetTextFont(42);
                            legMcOv.SetTextSize(0.032);
                            legMcOv.AddEntry(hData, TString::Format("AuAu data (%s)", H.label.c_str()).Data(), "ep");
                            legMcOv.AddEntry(hMC, mcCfg.mcLegend.c_str(), "l");
                            legMcOv.Draw();
                            
                            TLatex tMcOvTitle;
                            tMcOvTitle.SetNDC(true);
                            tMcOvTitle.SetTextFont(42);
                            tMcOvTitle.SetTextAlign(23);
                            tMcOvTitle.SetTextSize(0.040);
                            tMcOvTitle.DrawLatex(0.50, 0.97,
                                                 TString::Format("E_{T}^{iso} overlay: AuAu data vs %s Embedded MC", mcCfg.titleTag.c_str()).Data());
                            
                            TLatex tMcOvInfo;
                            tMcOvInfo.SetNDC(true);
                            tMcOvInfo.SetTextFont(42);
                            tMcOvInfo.SetTextAlign(13);
                            tMcOvInfo.SetTextSize(0.045);
                            tMcOvInfo.DrawLatex(0.18, 0.89, TString::Format("%d-%d%%", cb.lo, cb.hi).Data());
                            tMcOvInfo.DrawLatex(0.18, 0.85, TString::Format("p_{T}^{#gamma} = %d-%d GeV", b.lo, b.hi).Data());
                            
                            TLatex tSph;
                            tSph.SetNDC(true);
                            tSph.SetTextFont(42);
                            tSph.SetTextAlign(33);
                            tSph.SetTextSize(0.048);
                            tSph.DrawLatex(0.92, 0.6, "#bf{sPHENIX} #it{Internal}");
                            tSph.SetTextSize(0.038);
                            tSph.DrawLatex(0.92, 0.54, "Au+Au  #sqrt{s_{NN}} = 200 GeV");
                            
                            TLatex tUE;
                            tUE.SetNDC(true);
                            tUE.SetTextFont(42);
                            tUE.SetTextAlign(33);
                            tUE.SetTextSize(0.034);
                            tUE.DrawLatex(0.92, 0.46, trigDisplayLabel.c_str());
                            tUE.DrawLatex(0.92, 0.42, TString::Format("|v_{z}| < %d cm", kAA_VzCut).Data());
                            tUE.DrawLatex(0.92, 0.38, TString::Format("UE: %s", H.label.c_str()).Data());
                            tUE.DrawLatex(0.92, 0.34, TString::Format("#DeltaR_{cone} < %.1f", (kAA_IsoConeR == "isoR40") ? 0.4 : 0.3).Data());
                            
                            SaveCanvas(cMcOv, JoinPath(mcOvDir, "Eiso_dataMC_overlay.png"));
                            
                            // individual unnormalized distributions
                            for (const auto& rv : std::vector<std::pair<TH1*,std::string>>{
                                {hDataSrc, "Eiso_data_raw"}, {hMCSrc, "Eiso_mc_raw"}})
                            {
                                TH1* hRaw = dynamic_cast<TH1*>(rv.first->Clone(
                                                                               TString::Format("hRaw_%s_%s_%s_%s_%s", rv.second.c_str(),
                                                                                               mcCfg.folder.c_str(), H.variant.c_str(), cb.folder.c_str(), b.folder.c_str()).Data()));
                                if (!hRaw) continue;
                                hRaw->SetDirectory(nullptr);
                                hRaw->SetTitle("");
                                hRaw->GetXaxis()->SetTitle("E_{T}^{iso} [GeV]");
                                hRaw->GetYaxis()->SetTitle("Counts");
                                hRaw->GetXaxis()->SetTitleSize(0.055);
                                hRaw->GetYaxis()->SetTitleSize(0.055);
                                hRaw->GetXaxis()->SetLabelSize(0.045);
                                hRaw->GetYaxis()->SetLabelSize(0.045);
                                hRaw->GetYaxis()->SetTitleOffset(1.15);
                                const bool isMC = (rv.second.find("mc") != std::string::npos);
                                hRaw->SetLineColor(isMC ? kBlack : (kBlue + 1));
                                hRaw->SetLineWidth(2);
                                hRaw->SetFillStyle(0);
                                TCanvas cRaw(
                                             TString::Format("c_raw_%s_%s_%s_%s_%s", rv.second.c_str(),
                                                             mcCfg.folder.c_str(), H.variant.c_str(), cb.folder.c_str(), b.folder.c_str()).Data(),
                                             "c_raw", 900, 700);
                                ApplyCanvasMargins1D(cRaw);
                                cRaw.cd();
                                if (isMC)
                                { hRaw->SetMarkerSize(0.0); hRaw->Draw("hist"); }
                                else
                                { hRaw->SetMarkerColor(kBlue + 1); hRaw->SetMarkerStyle(20); hRaw->SetMarkerSize(1.0); hRaw->Draw("E1"); }
                                TLatex tRawInfo;
                                tRawInfo.SetNDC(true);
                                tRawInfo.SetTextFont(42);
                                tRawInfo.SetTextAlign(13);
                                tRawInfo.SetTextSize(0.045);
                                tRawInfo.DrawLatex(0.16, 0.88, TString::Format("%d-%d%%", cb.lo, cb.hi).Data());
                                tRawInfo.DrawLatex(0.16, 0.82, TString::Format("p_{T}^{#gamma} = %d-%d GeV", b.lo, b.hi).Data());
                                SaveCanvas(cRaw, JoinPath(mcOvDir, rv.second + ".png"));
                                delete hRaw;
                            }
                            
                            delete hData;
                            delete hMC;
                            
                            // ── signal/bkg decomposition: AuAu data + embedded MC ──
                            {
                                const string sigOvDir = JoinPath(mcOvDir, "signalIsoOverlays");
                                EnsureDir(sigOvDir);
                                const string hTName  = "h_Eiso_tight" + b.suffix + cb.suffix;
                                const string hNTName = "h_Eiso_nonTight" + b.suffix + cb.suffix;
                                const string hSigName = "h_EisoReco_truthSigMatched" + b.suffix + cb.suffix;
                                TH1* hTdSrc  = dynamic_cast<TH1*>(aaTop->Get(hTName.c_str()));
                                TH1* hNTdSrc = dynamic_cast<TH1*>(aaTop->Get(hNTName.c_str()));
                                TH1* hSigSrc = dynamic_cast<TH1*>(mcCfg.mcTop->Get(hSigName.c_str()));
                                if (hTdSrc && hNTdSrc && hSigSrc)
                                {
                                    TH1* hTd  = CloneTH1(hTdSrc,  TString::Format("hTd_aaSig_%s_%s_%s_%s", mcCfg.folder.c_str(), H.variant.c_str(), cb.folder.c_str(), b.folder.c_str()).Data());
                                    TH1* hNTd = CloneTH1(hNTdSrc, TString::Format("hNTd_aaSig_%s_%s_%s_%s", mcCfg.folder.c_str(), H.variant.c_str(), cb.folder.c_str(), b.folder.c_str()).Data());
                                    TH1* hSig = CloneTH1(hSigSrc, TString::Format("hSig_aaSig_%s_%s_%s_%s", mcCfg.folder.c_str(), H.variant.c_str(), cb.folder.c_str(), b.folder.c_str()).Data());
                                    if (hTd && hNTd && hSig)
                                    {
                                        hTd->Rebin(10); hNTd->Rebin(10); hSig->Rebin(10);
                                        { const double iTdAA=hTd->Integral(0,hTd->GetNbinsX()+1); const double iNTdAA=hNTd->Integral(0,hNTd->GetNbinsX()+1); const double iSigAA=hSig->Integral(0,hSig->GetNbinsX()+1);
                                            if (iTdAA>0) hTd->Scale(1.0/iTdAA); if (iNTdAA>0) hNTd->Scale(1.0/iNTdAA); if (iSigAA>0) hSig->Scale(1.0/iSigAA); }
                                        hSig->SetLineColor(kBlue-9); hSig->SetFillColorAlpha(kBlue-9,0.5); hSig->SetFillStyle(1001); hSig->SetLineWidth(1); hSig->SetMarkerSize(0.0);
                                        hNTd->SetLineColor(kPink-4); hNTd->SetFillColorAlpha(kPink-4,0.5); hNTd->SetFillStyle(1001); hNTd->SetLineWidth(1); hNTd->SetMarkerSize(0.0);
                                        for (int ib=0; ib<=hSig->GetNbinsX()+1; ++ib) hSig->SetBinError(ib,0.0);
                                        for (int ib=0; ib<=hNTd->GetNbinsX()+1; ++ib) hNTd->SetBinError(ib,0.0);
                                        hTd->SetLineColor(kBlack); hTd->SetMarkerColor(kBlack); hTd->SetMarkerStyle(20); hTd->SetMarkerSize(1.0); hTd->SetLineWidth(2); hTd->SetFillStyle(0);
                                        TH1* hStack = (TH1*)hSig->Clone(TString::Format("hStack_aaSig_%s_%s_%s_%s", mcCfg.folder.c_str(), H.variant.c_str(), cb.folder.c_str(), b.folder.c_str()).Data());
                                        hStack->SetDirectory(nullptr);
                                        hStack->Add(hNTd);
                                        const double ym = std::max(hTd->GetMaximum(), hStack->GetMaximum());
                                        TCanvas cSB(TString::Format("c_aaSigBkg_%s_%s_%s_%s", mcCfg.folder.c_str(), H.variant.c_str(), cb.folder.c_str(), b.folder.c_str()).Data(), "c_aaSigBkg", 900, 700);
                                        ApplyCanvasMargins1D(cSB); cSB.cd();
                                        hStack->SetTitle(""); hStack->GetXaxis()->SetTitle("E_{T}^{iso} [GeV]");
                                        hStack->GetYaxis()->SetTitle("Normalized to unit area");
                                        hStack->GetXaxis()->SetTitleSize(0.055); hStack->GetYaxis()->SetTitleSize(0.055);
                                        hStack->GetXaxis()->SetLabelSize(0.045); hStack->GetYaxis()->SetLabelSize(0.045);
                                        hStack->GetYaxis()->SetTitleOffset(1.15);
                                        hStack->SetMinimum(0.0); hStack->SetMaximum((ym>0)?(1.3*ym):1.0);
                                        hStack->SetLineColor(kPink-4); hStack->SetFillColorAlpha(kPink-4,0.5); hStack->SetFillStyle(1001);
                                        hStack->Draw("hist");
                                        hSig->Draw("hist SAME");
                                        hTd->Draw("E1 SAME");
                                        TLatex tSBtitle; tSBtitle.SetNDC(true); tSBtitle.SetTextFont(42); tSBtitle.SetTextAlign(23); tSBtitle.SetTextSize(0.040);
                                        tSBtitle.DrawLatex(0.50, 0.97,
                                                           TString::Format("E_{T}^{iso} overlay: AuAu data vs %s Embedded MC", mcCfg.titleTag.c_str()).Data());
                                        TLegend lg(0.56,0.65,0.92,0.88); lg.SetBorderSize(0); lg.SetFillStyle(0); lg.SetTextFont(42); lg.SetTextSize(0.032);
                                        lg.AddEntry(hTd, "Data (Signal)", "ep"); lg.AddEntry(hNTd, "Data (Background)", "f"); lg.AddEntry(hSig, "Signal Embedded MC", "f"); lg.Draw();
                                        TLatex ts; ts.SetNDC(true); ts.SetTextFont(42); ts.SetTextAlign(33); ts.SetTextSize(0.048);
                                        ts.DrawLatex(0.92,0.6,"#bf{sPHENIX} #it{Internal}"); ts.SetTextSize(0.038);
                                        ts.DrawLatex(0.92,0.54,"Au+Au  #sqrt{s_{NN}} = 200 GeV");
                                        TLatex ti; ti.SetNDC(true); ti.SetTextFont(42); ti.SetTextAlign(13); ti.SetTextSize(0.034);
                                        ti.DrawLatex(0.18,0.89,TString::Format("%d-%d%%", cb.lo, cb.hi).Data());
                                        ti.DrawLatex(0.18,0.85,TString::Format("p_{T}^{#gamma} = %d-%d GeV", b.lo, b.hi).Data());
                                        TLatex tInfo2; tInfo2.SetNDC(true); tInfo2.SetTextFont(42); tInfo2.SetTextAlign(33); tInfo2.SetTextSize(0.034);
                                        tInfo2.DrawLatex(0.92,0.46, trigDisplayLabel.c_str());
                                        tInfo2.DrawLatex(0.92,0.42, TString::Format("|v_{z}| < %d cm", kAA_VzCut).Data());
                                        tInfo2.DrawLatex(0.92,0.38, TString::Format("UE: %s", H.label.c_str()).Data());
                                        tInfo2.DrawLatex(0.92,0.34, TString::Format("#DeltaR_{cone} < %.1f", (kAA_IsoConeR == "isoR40") ? 0.4 : 0.3).Data());
                                        SaveCanvas(cSB, JoinPath(sigOvDir, "Eiso_sigBkg_overlay.png"));
                                    }
                                    if (hTd) delete hTd; if (hNTd) delete hNTd; if (hSig) delete hSig;
                                }
                            }
                        }
                    }
                }
                
                // ── Merged pT 20–35 bin: combine last N bins where lo >= 20 ──
                {
                    const int mergedLo = 20;
                    const int mergedHi = 35;
                    const string mergedFolder = TString::Format("pT_%d_%d", mergedLo, mergedHi).Data();
                    const string mergedSuffix = TString::Format("_pT_%d_%d", mergedLo, mergedHi).Data();
                    const string mergedPtDir = JoinPath(centDir, mergedFolder);
                    EnsureDir(mergedPtDir);
                    
                    // identify which pT bins to merge
                    vector<int> mergeBins;
                    for (int ip = 0; ip < kNPtBins; ++ip)
                    {
                        const PtBin& bp = PtBins()[ip];
                        if (bp.lo >= mergedLo && bp.hi <= mergedHi) mergeBins.push_back(ip);
                    }
                    
                    if (!mergeBins.empty())
                    {
                        // helper: merge histograms from a TDirectory by summing bins
                        auto MergeHists = [&](TDirectory* dir, const string& baseName, const string& centSuffix, const string& cloneName) -> TH1*
                        {
                            TH1* hSum = nullptr;
                            for (int ip : mergeBins)
                            {
                                const PtBin& bp = PtBins()[ip];
                                const string hName = baseName + bp.suffix + centSuffix;
                                TH1* hSrc = dynamic_cast<TH1*>(dir->Get(hName.c_str()));
                                if (!hSrc) continue;
                                if (!hSum)
                                {
                                    hSum = CloneTH1(hSrc, cloneName.c_str());
                                    if (hSum) EnsureSumw2(hSum);
                                }
                                else
                                {
                                    hSum->Add(hSrc);
                                }
                            }
                            return hSum;
                        };
                        
                        // --- pp vs AuAu overlay ---
                        {
                            TH1* hPPm = MergeHists(ppTop, "h_Eiso", "", "hPP_merged_" + mergedFolder + "_" + cb.folder);
                            TH1* hAAm = MergeHists(aaTop, "h_Eiso", cb.suffix, "hAA_merged_" + mergedFolder + "_" + cb.folder);
                            if (hPPm && hAAm)
                            {
                                const double iPP = hPPm->Integral(0, hPPm->GetNbinsX()+1);
                                const double iAA = hAAm->Integral(0, hAAm->GetNbinsX()+1);
                                if (iPP > 0.0) hPPm->Scale(1.0/iPP);
                                if (iAA > 0.0) hAAm->Scale(1.0/iAA);
                                hPPm->SetLineColor(kRed+1); hPPm->SetMarkerColor(kRed+1); hPPm->SetMarkerStyle(24); hPPm->SetMarkerSize(1.1); hPPm->SetLineWidth(2); hPPm->SetFillStyle(0);
                                hAAm->SetLineColor(kBlack); hAAm->SetMarkerColor(kBlack); hAAm->SetMarkerStyle(20); hAAm->SetMarkerSize(1.1); hAAm->SetLineWidth(2); hAAm->SetFillStyle(0);
                                const double ym = std::max(hPPm->GetMaximum(), hAAm->GetMaximum());
                                TCanvas cM("c_merged_ppAuAu","c_merged_ppAuAu",900,700); ApplyCanvasMargins1D(cM); cM.cd();
                                hAAm->SetTitle(""); hAAm->GetXaxis()->SetTitle("E_{T}^{iso} [GeV]"); hAAm->GetYaxis()->SetTitle("Normalized to unit area");
                                hAAm->GetXaxis()->SetTitleSize(0.055); hAAm->GetYaxis()->SetTitleSize(0.055); hAAm->GetXaxis()->SetLabelSize(0.045); hAAm->GetYaxis()->SetLabelSize(0.045); hAAm->GetYaxis()->SetTitleOffset(1.15);
                                hAAm->SetMinimum(0.0); hAAm->SetMaximum((ym>0)?(1.25*ym):1.0);
                                hAAm->Draw("E1"); hPPm->Draw("E1 SAME");
                                TLegend lgM(0.56,0.68,0.92,0.88); lgM.SetBorderSize(0); lgM.SetFillStyle(0); lgM.SetTextFont(42); lgM.SetTextSize(0.032);
                                lgM.AddEntry(hPPm,"pp data","ep"); lgM.AddEntry(hAAm,TString::Format("AuAu data (%s)",H.label.c_str()).Data(),"ep"); lgM.Draw();
                                TLatex tMt; tMt.SetNDC(true); tMt.SetTextFont(42); tMt.SetTextAlign(23); tMt.SetTextSize(0.045);
                                tMt.DrawLatex(0.50,0.955,TString::Format("isolation overlay pp and %d-%d%% Cent AuAu",cb.lo,cb.hi).Data());
                                TLatex tMi; tMi.SetNDC(true); tMi.SetTextFont(42); tMi.SetTextAlign(13); tMi.SetTextSize(0.034);
                                tMi.DrawLatex(0.16,0.88,trigDisplayLabel.c_str());
                                tMi.DrawLatex(0.16,0.84,TString::Format("|v_{z}| < %d cm",kAA_VzCut).Data());
                                tMi.DrawLatex(0.16,0.80,TString::Format("UE subtraction: %s",H.label.c_str()).Data());
                                tMi.DrawLatex(0.16,0.76,TString::Format("p_{T}^{#gamma}: %d-%d GeV",mergedLo,mergedHi).Data());
                                const string ppAuAuMDir = JoinPath(mergedPtDir,"ppAuAuOverlay");
                                EnsureDir(ppAuAuMDir);
                                const string furRealMDir = JoinPath(ppAuAuMDir,"ppAuAuOverlayFURREAL");
                                EnsureDir(furRealMDir);
                                SaveCanvas(cM,JoinPath(furRealMDir,"isolationOverlay_pp_vs_auau.png"));
                            }
                            if (hPPm) delete hPPm;
                            if (hAAm) delete hAAm;
                        }
                        
                        // --- MC overlays: inclusive + photon+jet ---
                        if (!forEmbeddedSim)
                        {
                            struct MCOverlayCfgM {
                                TDirectory* mcTop;
                                string folder;
                                string titleTag;
                                string mcLegend;
                            };
                            const vector<MCOverlayCfgM> mcOvCfgsM = {
                                { incMCvarTop, "inclusiveMCoverlays", inclusiveEmbeddedTitleTag, "inclusive MC" },
                                { phoMCvarTop, "photonJetOverlays",   photonEmbeddedTitleTag,    "photon+jet MC" }
                            };
                            
                            for (const auto& mcCfg : mcOvCfgsM)
                            {
                                if (!mcCfg.mcTop || mcCfg.titleTag.empty()) continue;
                                const string mcOvDirM = JoinPath(mergedPtDir, mcCfg.folder);
                                EnsureDir(mcOvDirM);
                                
                                TH1* hDataM = MergeHists(aaTop, "h_Eiso", cb.suffix, "hData_merged_" + mcCfg.folder + "_" + cb.folder);
                                TH1* hMCM   = MergeHists(mcCfg.mcTop, "h_Eiso", cb.suffix, "hMC_merged_" + mcCfg.folder + "_" + cb.folder);
                                if (hDataM && hMCM)
                                {
                                    hDataM->Rebin(10); hMCM->Rebin(10);
                                    const double iD = hDataM->Integral(0,hDataM->GetNbinsX()+1);
                                    const double iM = hMCM->Integral(0,hMCM->GetNbinsX()+1);
                                    if (iD > 0.0 && iM > 0.0)
                                    {
                                        hDataM->Scale(1.0/iD); hMCM->Scale(1.0/iM);
                                        for (int ib=0;ib<=hMCM->GetNbinsX()+1;++ib) hMCM->SetBinError(ib,0.0);
                                        hMCM->SetTitle(""); hMCM->SetLineColor(kBlack); hMCM->SetLineWidth(2); hMCM->SetFillStyle(0); hMCM->SetMarkerSize(0.0);
                                        hDataM->SetLineColor(kBlue+1); hDataM->SetMarkerColor(kBlue+1); hDataM->SetMarkerStyle(20); hDataM->SetMarkerSize(1.0); hDataM->SetLineWidth(2); hDataM->SetFillStyle(0);
                                        const double ymx = std::max(hDataM->GetMaximum(),hMCM->GetMaximum());
                                        TCanvas cMov("c_merged_mcOv","c_merged_mcOv",900,700); ApplyCanvasMargins1D(cMov); cMov.cd();
                                        hMCM->GetXaxis()->SetTitle("E_{T}^{iso} [GeV]"); hMCM->GetYaxis()->SetTitle("Normalized to unit area");
                                        hMCM->GetXaxis()->SetTitleSize(0.055); hMCM->GetYaxis()->SetTitleSize(0.055); hMCM->GetXaxis()->SetLabelSize(0.045); hMCM->GetYaxis()->SetLabelSize(0.045); hMCM->GetYaxis()->SetTitleOffset(1.15);
                                        hMCM->SetMinimum(0.0); hMCM->SetMaximum((ymx>0)?(1.25*ymx):1.0);
                                        hMCM->Draw("hist"); hDataM->Draw("E1 SAME");
                                        TLegend lgOv(0.56,0.72,0.92,0.88); lgOv.SetBorderSize(0); lgOv.SetFillStyle(0); lgOv.SetTextFont(42); lgOv.SetTextSize(0.032);
                                        lgOv.AddEntry(hDataM,TString::Format("AuAu data (%s)",H.label.c_str()).Data(),"ep");
                                        lgOv.AddEntry(hMCM,mcCfg.mcLegend.c_str(),"l"); lgOv.Draw();
                                        TLatex tOvT; tOvT.SetNDC(true); tOvT.SetTextFont(42); tOvT.SetTextAlign(23); tOvT.SetTextSize(0.040);
                                        tOvT.DrawLatex(0.50,0.97,TString::Format("E_{T}^{iso} overlay: AuAu data vs %s Embedded MC",mcCfg.titleTag.c_str()).Data());
                                        TLatex tOvI; tOvI.SetNDC(true); tOvI.SetTextFont(42); tOvI.SetTextAlign(13); tOvI.SetTextSize(0.045);
                                        tOvI.DrawLatex(0.18,0.89,TString::Format("%d-%d%%",cb.lo,cb.hi).Data());
                                        tOvI.DrawLatex(0.18,0.85,TString::Format("p_{T}^{#gamma} = %d-%d GeV",mergedLo,mergedHi).Data());
                                        TLatex tSphM; tSphM.SetNDC(true); tSphM.SetTextFont(42); tSphM.SetTextAlign(33); tSphM.SetTextSize(0.048);
                                        tSphM.DrawLatex(0.92,0.6,"#bf{sPHENIX} #it{Internal}"); tSphM.SetTextSize(0.038);
                                        tSphM.DrawLatex(0.92,0.54,"Au+Au  #sqrt{s_{NN}} = 200 GeV");
                                        TLatex tUEM; tUEM.SetNDC(true); tUEM.SetTextFont(42); tUEM.SetTextAlign(33); tUEM.SetTextSize(0.034);
                                        tUEM.DrawLatex(0.92,0.46,trigDisplayLabel.c_str());
                                        tUEM.DrawLatex(0.92,0.42,TString::Format("|v_{z}| < %d cm",kAA_VzCut).Data());
                                        tUEM.DrawLatex(0.92,0.38,TString::Format("UE: %s",H.label.c_str()).Data());
                                        tUEM.DrawLatex(0.92,0.34,TString::Format("#DeltaR_{cone} < %.1f",(kAA_IsoConeR=="isoR40")?0.4:0.3).Data());
                                        SaveCanvas(cMov,JoinPath(mcOvDirM,"Eiso_dataMC_overlay.png"));
                                    }
                                }
                                if (hDataM) delete hDataM;
                                if (hMCM) delete hMCM;
                                
                                // --- signal/bkg decomposition for merged bin ---
                                {
                                    const string sigOvDirM = JoinPath(mcOvDirM,"signalIsoOverlays");
                                    EnsureDir(sigOvDirM);
                                    TH1* hTdM  = MergeHists(aaTop,      "h_Eiso_tight",    cb.suffix, "hTd_merged_"+mcCfg.folder+"_"+cb.folder);
                                    TH1* hNTdM = MergeHists(aaTop,      "h_Eiso_nonTight", cb.suffix, "hNTd_merged_"+mcCfg.folder+"_"+cb.folder);
                                    TH1* hSigM = MergeHists(mcCfg.mcTop,"h_EisoReco_truthSigMatched", cb.suffix, "hSig_merged_"+mcCfg.folder+"_"+cb.folder);
                                    if (hTdM && hNTdM && hSigM)
                                    {
                                        hTdM->Rebin(10); hNTdM->Rebin(10); hSigM->Rebin(10);
                                        { const double iTdMN=hTdM->Integral(0,hTdM->GetNbinsX()+1); const double iNTdMN=hNTdM->Integral(0,hNTdM->GetNbinsX()+1); const double iSigMN=hSigM->Integral(0,hSigM->GetNbinsX()+1);
                                            if (iTdMN>0) hTdM->Scale(1.0/iTdMN); if (iNTdMN>0) hNTdM->Scale(1.0/iNTdMN); if (iSigMN>0) hSigM->Scale(1.0/iSigMN); }
                                        hSigM->SetLineColor(kBlue-9); hSigM->SetFillColorAlpha(kBlue-9,0.5); hSigM->SetFillStyle(1001); hSigM->SetLineWidth(1); hSigM->SetMarkerSize(0.0);
                                        hNTdM->SetLineColor(kPink-4); hNTdM->SetFillColorAlpha(kPink-4,0.5); hNTdM->SetFillStyle(1001); hNTdM->SetLineWidth(1); hNTdM->SetMarkerSize(0.0);
                                        for (int ib=0;ib<=hSigM->GetNbinsX()+1;++ib) hSigM->SetBinError(ib,0.0);
                                        for (int ib=0;ib<=hNTdM->GetNbinsX()+1;++ib) hNTdM->SetBinError(ib,0.0);
                                        hTdM->SetLineColor(kBlack); hTdM->SetMarkerColor(kBlack); hTdM->SetMarkerStyle(20); hTdM->SetMarkerSize(1.0); hTdM->SetLineWidth(2); hTdM->SetFillStyle(0);
                                        TH1* hStackM = (TH1*)hSigM->Clone("hStack_merged"); hStackM->SetDirectory(nullptr); hStackM->Add(hNTdM);
                                        const double ymM = std::max(hTdM->GetMaximum(),hStackM->GetMaximum());
                                        TCanvas cSBm("c_merged_sigBkg","c_merged_sigBkg",900,700); ApplyCanvasMargins1D(cSBm); cSBm.cd();
                                        hStackM->SetTitle(""); hStackM->GetXaxis()->SetTitle("E_{T}^{iso} [GeV]"); hStackM->GetYaxis()->SetTitle("Normalized to unit area");
                                        hStackM->GetXaxis()->SetTitleSize(0.055); hStackM->GetYaxis()->SetTitleSize(0.055); hStackM->GetXaxis()->SetLabelSize(0.045); hStackM->GetYaxis()->SetLabelSize(0.045); hStackM->GetYaxis()->SetTitleOffset(1.15);
                                        hStackM->SetMinimum(0.0); hStackM->SetMaximum((ymM>0)?(1.3*ymM):1.0);
                                        hStackM->SetLineColor(kPink-4); hStackM->SetFillColorAlpha(kPink-4,0.5); hStackM->SetFillStyle(1001);
                                        hStackM->Draw("hist"); hSigM->Draw("hist SAME"); hTdM->Draw("E1 SAME");
                                        TLatex tSBtM; tSBtM.SetNDC(true); tSBtM.SetTextFont(42); tSBtM.SetTextAlign(23); tSBtM.SetTextSize(0.040);
                                        tSBtM.DrawLatex(0.50,0.97,TString::Format("E_{T}^{iso} overlay: AuAu data vs %s Embedded MC",mcCfg.titleTag.c_str()).Data());
                                        TLegend lgSBm(0.56,0.65,0.92,0.88); lgSBm.SetBorderSize(0); lgSBm.SetFillStyle(0); lgSBm.SetTextFont(42); lgSBm.SetTextSize(0.032);
                                        lgSBm.AddEntry(hTdM,"Data (Signal)","ep"); lgSBm.AddEntry(hNTdM,"Data (Background)","f"); lgSBm.AddEntry(hSigM,"Signal Embedded MC","f"); lgSBm.Draw();
                                        TLatex tsM; tsM.SetNDC(true); tsM.SetTextFont(42); tsM.SetTextAlign(33); tsM.SetTextSize(0.048);
                                        tsM.DrawLatex(0.92,0.6,"#bf{sPHENIX} #it{Internal}"); tsM.SetTextSize(0.038);
                                        tsM.DrawLatex(0.92,0.54,"Au+Au  #sqrt{s_{NN}} = 200 GeV");
                                        TLatex tiM; tiM.SetNDC(true); tiM.SetTextFont(42); tiM.SetTextAlign(13); tiM.SetTextSize(0.034);
                                        tiM.DrawLatex(0.18,0.89,TString::Format("%d-%d%%",cb.lo,cb.hi).Data());
                                        tiM.DrawLatex(0.18,0.85,TString::Format("p_{T}^{#gamma} = %d-%d GeV",mergedLo,mergedHi).Data());
                                        TLatex tI2M; tI2M.SetNDC(true); tI2M.SetTextFont(42); tI2M.SetTextAlign(33); tI2M.SetTextSize(0.034);
                                        tI2M.DrawLatex(0.92,0.46,trigDisplayLabel.c_str());
                                        tI2M.DrawLatex(0.92,0.42,TString::Format("|v_{z}| < %d cm",kAA_VzCut).Data());
                                        tI2M.DrawLatex(0.92,0.38,TString::Format("UE: %s",H.label.c_str()).Data());
                                        tI2M.DrawLatex(0.92,0.34,TString::Format("#DeltaR_{cone} < %.1f",(kAA_IsoConeR=="isoR40")?0.4:0.3).Data());
                                        SaveCanvas(cSBm,JoinPath(sigOvDirM,"Eiso_sigBkg_overlay.png"));
                                        delete hStackM;
                                    }
                                    if (hTdM) delete hTdM;
                                    if (hNTdM) delete hNTdM;
                                    if (hSigM) delete hSigM;
                                }
                            }
                        }
                    }
                }
                if (!xPt.empty())
                {
                    TCanvas cMean(
                                  TString::Format("c_meanIsoEtPPAuAu_%s_%s_%s",
                                                  trigAA.c_str(), H.variant.c_str(), cb.folder.c_str()).Data(),
                                  "c_meanIsoEtPPAuAu", 900, 700
                                  );
                    ApplyCanvasMargins1D(cMean);
                    cMean.cd();
                    
                    double yMin = std::numeric_limits<double>::max();
                    double yMax = -std::numeric_limits<double>::max();
                    for (std::size_t i = 0; i < xPt.size(); ++i)
                    {
                        yMin = std::min(yMin, std::min(yPP[i] - eyPP[i], yAA[i] - eyAA[i]));
                        yMax = std::max(yMax, std::max(yPP[i] + eyPP[i], yAA[i] + eyAA[i]));
                    }
                    if (!std::isfinite(yMin) || !std::isfinite(yMax))
                    {
                        yMin = 0.0;
                        yMax = 1.0;
                    }
                    const double pad = (yMax > yMin) ? (0.7 * (yMax - yMin)) : 0.25;
                    
                    TH1F hFrame(
                                TString::Format("hFrame_meanIsoEtPPAuAu_%s_%s_%s",
                                                trigAA.c_str(), H.variant.c_str(), cb.folder.c_str()).Data(),
                                "", 100, kPtEdges.front(), kPtEdges.back()
                                );
                    hFrame.SetDirectory(nullptr);
                    hFrame.SetStats(0);
                    hFrame.SetMinimum(std::max(0.0, yMin - pad));
                    hFrame.SetMaximum(yMax + pad);
                    hFrame.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                    hFrame.GetYaxis()->SetTitle("<E_{T}^{iso}> [GeV]");
                    hFrame.GetXaxis()->SetTitleSize(0.055);
                    hFrame.GetYaxis()->SetTitleSize(0.055);
                    hFrame.GetXaxis()->SetLabelSize(0.045);
                    hFrame.GetYaxis()->SetLabelSize(0.045);
                    hFrame.GetYaxis()->SetTitleOffset(1.15);
                    hFrame.Draw();
                    
                    TGraphErrors gPP((int)xPt.size(), &xPt[0], &yPP[0], &exPt[0], &eyPP[0]);
                    gPP.SetLineWidth(2);
                    gPP.SetLineColor(kRed + 1);
                    gPP.SetMarkerColor(kRed + 1);
                    gPP.SetMarkerStyle(24);
                    gPP.SetMarkerSize(1.2);
                    gPP.Draw("PE1 SAME");
                    
                    TGraphErrors gAA((int)xPt.size(), &xPt[0], &yAA[0], &exPt[0], &eyAA[0]);
                    gAA.SetLineWidth(2);
                    gAA.SetLineColor(kBlack);
                    gAA.SetMarkerColor(kBlack);
                    gAA.SetMarkerStyle(20);
                    gAA.SetMarkerSize(1.2);
                    gAA.Draw("PE1 SAME");
                    
                    TLegend leg(0.56, 0.68, 0.92, 0.88);
                    leg.SetBorderSize(0);
                    leg.SetFillStyle(0);
                    leg.SetTextFont(42);
                    leg.SetTextSize(0.032);
                    leg.AddEntry(&gPP, "pp data", "ep");
                    leg.AddEntry(&gAA, TString::Format("AuAu data (%s)", H.label.c_str()).Data(), "ep");
                    leg.Draw();
                    
                    TLatex tTitle;
                    tTitle.SetNDC(true);
                    tTitle.SetTextFont(42);
                    tTitle.SetTextAlign(23);
                    tTitle.SetTextSize(0.045);
                    tTitle.DrawLatex(0.50, 0.955,
                                     TString::Format("<E_{T}^{iso}> overlay pp and %d-%d%% Cent AuAu", cb.lo, cb.hi).Data());
                    
                    TLatex t;
                    t.SetNDC(true);
                    t.SetTextFont(42);
                    t.SetTextAlign(13);
                    t.SetTextSize(0.034);
                    t.DrawLatex(0.16, 0.88, TString::Format("Trigger: %s", trigAA.c_str()).Data());
                    t.DrawLatex(0.16, 0.84, TString::Format("|v_{z}| < %d cm", kAA_VzCut).Data());
                    t.DrawLatex(0.16, 0.80, TString::Format("UE subtraction: %s", H.label.c_str()).Data());
                    
                    {
                        const string ppAuAuVsPtDir = JoinPath(centDir, "ppAuAu_meanIsoEt_versusPT");
                        EnsureDir(ppAuAuVsPtDir);
                        SaveCanvas(cMean, JoinPath(ppAuAuVsPtDir,
                                                   TString::Format("meanIsoEt_pp_vs_auau_vs_pT_%s.png", H.variant.c_str()).Data()));
                    }
                }
                
                // ── MC embedded <E_T^iso> vs pT per centrality folder ──
                {
                    struct MCvsPtCfg {
                        const char* folderName;
                        string mcLabel;
                        const vector<vector<double>>& yMC;
                        const vector<vector<double>>& eyMC;
                        const vector<vector<bool>>&   filledMC;
                        const vector<vector<double>>& ySigMC;
                        const vector<vector<double>>& eySigMC;
                        const vector<vector<bool>>&   filledSigMC;
                    };
                    const string incJetSampleName = InclusiveEmbeddedShortLabel().empty()
                    ? "" : "Inclusive Jet " + InclusiveEmbeddedShortLabel();
                    const string phoJetSampleName = PhotonEmbeddedShortLabel().empty()
                    ? "" : "Photon+Jet " + PhotonEmbeddedShortLabel();
                    const vector<MCvsPtCfg> mcPtCfgs = {
                        {"AuAu_InclusiveJetEmbedded_meanIsoEt_versusPT",
                            incJetSampleName,
                            vsCent_yIncMC, vsCent_eyIncMC, vsCent_filledIncMC,
                            vsCent_ySigIncMC, vsCent_eySigIncMC, vsCent_filledSigIncMC},
                        {"AuAu_PhotonJetEmbedded_meanIsoEt_versusPT",
                            phoJetSampleName,
                            vsCent_yPhoMC, vsCent_eyPhoMC, vsCent_filledPhoMC,
                            vsCent_ySigPhoMC, vsCent_eySigPhoMC, vsCent_filledSigPhoMC}
                    };
                    for (const auto& mcPt : mcPtCfgs)
                    {
                        vector<double> xPtMC, exPtMC, yPPmc, eyPPmc, yMCmc, eyMCmc;
                        for (int ipt = 0; ipt < kNPtBins; ++ipt)
                        {
                            if (!mcPt.filledMC[ipt][ic]) continue;
                            const PtBin& bp = PtBins()[ipt];
                            xPtMC.push_back(0.5 * (kPtEdges[(std::size_t)ipt] + kPtEdges[(std::size_t)ipt + 1]));
                            exPtMC.push_back(0.5 * (kPtEdges[(std::size_t)ipt + 1] - kPtEdges[(std::size_t)ipt]));
                            yPPmc.push_back(vsCent_yAA[ipt][ic]);
                            eyPPmc.push_back(vsCent_eyAA[ipt][ic]);
                            yMCmc.push_back(mcPt.yMC[ipt][ic]);
                            eyMCmc.push_back(mcPt.eyMC[ipt][ic]);
                        }
                        if (xPtMC.empty()) continue;
                        
                        const string mcPtDir = JoinPath(centDir, mcPt.folderName);
                        EnsureDir(mcPtDir);
                        
                        double yMinMC = std::numeric_limits<double>::max();
                        double yMaxMC = -std::numeric_limits<double>::max();
                        for (std::size_t i = 0; i < xPtMC.size(); ++i)
                        {
                            yMinMC = std::min(yMinMC, std::min(yPPmc[i] - eyPPmc[i], yMCmc[i] - eyMCmc[i]));
                            yMaxMC = std::max(yMaxMC, std::max(yPPmc[i] + eyPPmc[i], yMCmc[i] + eyMCmc[i]));
                        }
                        if (!std::isfinite(yMinMC) || !std::isfinite(yMaxMC)) { yMinMC = 0.0; yMaxMC = 1.0; }
                        const double padMC = (yMaxMC > yMinMC) ? (0.7 * (yMaxMC - yMinMC)) : 0.25;
                        
                        TCanvas cMCpt(
                                      TString::Format("c_meanIsoMCvsPt_%s_%s_%s", mcPt.folderName, H.variant.c_str(), cb.folder.c_str()).Data(),
                                      "c_meanIsoMCvsPt", 900, 700);
                        ApplyCanvasMargins1D(cMCpt); cMCpt.cd();
                        
                        TH1F hFrMC(
                                   TString::Format("hFr_meanIsoMCvsPt_%s_%s_%s", mcPt.folderName, H.variant.c_str(), cb.folder.c_str()).Data(),
                                   "", 100, kPtEdges.front(), kPtEdges.back());
                        hFrMC.SetDirectory(nullptr); hFrMC.SetStats(0);
                        hFrMC.SetMinimum(std::max(0.0, yMinMC - padMC));
                        hFrMC.SetMaximum(yMaxMC + padMC);
                        hFrMC.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                        hFrMC.GetYaxis()->SetTitle("<E_{T}^{iso}> [GeV]");
                        hFrMC.GetXaxis()->SetTitleSize(0.055); hFrMC.GetYaxis()->SetTitleSize(0.055);
                        hFrMC.GetXaxis()->SetLabelSize(0.045); hFrMC.GetYaxis()->SetLabelSize(0.045);
                        hFrMC.GetYaxis()->SetTitleOffset(1.15);
                        hFrMC.Draw();
                        
                        TGraphErrors gPPmc((int)xPtMC.size(), &xPtMC[0], &yPPmc[0], &exPtMC[0], &eyPPmc[0]);
                        gPPmc.SetLineWidth(2); gPPmc.SetLineColor(kRed+1); gPPmc.SetMarkerColor(kRed+1);
                        gPPmc.SetMarkerStyle(24); gPPmc.SetMarkerSize(1.2);
                        gPPmc.Draw("PE1 SAME");
                        
                        TGraphErrors gMCmc((int)xPtMC.size(), &xPtMC[0], &yMCmc[0], &exPtMC[0], &eyMCmc[0]);
                        gMCmc.SetLineWidth(2); gMCmc.SetLineColor(kBlue+1); gMCmc.SetMarkerColor(kBlue+1);
                        gMCmc.SetMarkerStyle(20); gMCmc.SetMarkerSize(1.2);
                        gMCmc.Draw("PE1 SAME");
                        
                        TLegend lgMC(0.58, 0.58, 0.92, 0.72);
                        lgMC.SetBorderSize(0); lgMC.SetFillStyle(0); lgMC.SetTextFont(42); lgMC.SetTextSize(0.032);
                        lgMC.AddEntry(&gPPmc, TString::Format("AuAu data (%s)", H.label.c_str()).Data(), "ep");
                        lgMC.AddEntry(&gMCmc, TString::Format("%s Embedded MC", mcPt.mcLabel.c_str()).Data(), "ep");
                        lgMC.Draw();
                        
                        TLatex tTmc; tTmc.SetNDC(true); tTmc.SetTextFont(42); tTmc.SetTextAlign(23); tTmc.SetTextSize(0.038);
                        tTmc.DrawLatex(0.50, 0.955,
                                       TString::Format("<E_{T}^{iso}> Overlay, Au+Au & %s Embedded MC", mcPt.mcLabel.c_str()).Data());
                        
                        TLatex tMc; tMc.SetNDC(true); tMc.SetTextFont(42); tMc.SetTextAlign(13); tMc.SetTextSize(0.034);
                        tMc.DrawLatex(0.16, 0.88, trigDisplayLabel.c_str());
                        tMc.DrawLatex(0.16, 0.84, TString::Format("|v_{z}| < %d cm", kAA_VzCut).Data());
                        tMc.DrawLatex(0.16, 0.80, TString::Format("UE subtraction: %s", H.label.c_str()).Data());
                        tMc.DrawLatex(0.16, 0.76, TString::Format("#DeltaR_{cone} < %.1f", (kAA_IsoConeR == "isoR40") ? 0.4 : 0.3).Data());
                        tMc.DrawLatex(0.16, 0.72, TString::Format("Cent: %d-%d%%", cb.lo, cb.hi).Data());
                        
                        SaveCanvas(cMCpt, JoinPath(mcPtDir,
                                                   TString::Format("meanIsoEt_pp_vs_%s_vs_pT_%s.png",
                                                                   mcPt.folderName, H.variant.c_str()).Data()));
                        
                        // ── Signal-only <E_T^iso> vs pT ──
                        {
                            vector<double> xSig, exSig, yTDpt, eyTDpt, ySMpt, eySMpt;
                            for (int ipt = 0; ipt < kNPtBins; ++ipt)
                            {
                                if (!mcPt.filledSigMC[ipt][ic] && !vsCent_filledTightData[ipt][ic]) continue;
                                xSig.push_back(0.5 * (kPtEdges[(std::size_t)ipt] + kPtEdges[(std::size_t)ipt + 1]));
                                exSig.push_back(0.5 * (kPtEdges[(std::size_t)ipt + 1] - kPtEdges[(std::size_t)ipt]));
                                yTDpt.push_back(vsCent_filledTightData[ipt][ic] ? vsCent_yTightData[ipt][ic] : 0.0);
                                eyTDpt.push_back(vsCent_filledTightData[ipt][ic] ? vsCent_eyTightData[ipt][ic] : 0.0);
                                ySMpt.push_back(mcPt.filledSigMC[ipt][ic] ? mcPt.ySigMC[ipt][ic] : 0.0);
                                eySMpt.push_back(mcPt.filledSigMC[ipt][ic] ? mcPt.eySigMC[ipt][ic] : 0.0);
                            }
                            if (!xSig.empty())
                            {
                                const string sigPtDir = JoinPath(centDir, string(mcPt.folderName) + "_SIGNALonly");
                                EnsureDir(sigPtDir);
                                
                                double yMinS = std::numeric_limits<double>::max();
                                double yMaxS = -std::numeric_limits<double>::max();
                                for (std::size_t i = 0; i < xSig.size(); ++i)
                                {
                                    yMinS = std::min(yMinS, std::min(yTDpt[i] - eyTDpt[i], ySMpt[i] - eySMpt[i]));
                                    yMaxS = std::max(yMaxS, std::max(yTDpt[i] + eyTDpt[i], ySMpt[i] + eySMpt[i]));
                                }
                                if (!std::isfinite(yMinS) || !std::isfinite(yMaxS)) { yMinS = 0.0; yMaxS = 1.0; }
                                const double padS = (yMaxS > yMinS) ? (0.7 * (yMaxS - yMinS)) : 0.25;
                                
                                TCanvas cSigPt(
                                               TString::Format("c_meanIsoSigVsPt_%s_%s_%s", mcPt.folderName, H.variant.c_str(), cb.folder.c_str()).Data(),
                                               "c_meanIsoSigVsPt", 900, 700);
                                ApplyCanvasMargins1D(cSigPt); cSigPt.cd();
                                
                                TH1F hFrS(
                                          TString::Format("hFr_meanIsoSigVsPt_%s_%s_%s", mcPt.folderName, H.variant.c_str(), cb.folder.c_str()).Data(),
                                          "", 100, kPtEdges.front(), kPtEdges.back());
                                hFrS.SetDirectory(nullptr); hFrS.SetStats(0);
                                hFrS.SetMinimum(std::max(0.0, yMinS - padS));
                                hFrS.SetMaximum(yMaxS + padS);
                                hFrS.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                                hFrS.GetYaxis()->SetTitle("<E_{T}^{iso}> [GeV]");
                                hFrS.GetXaxis()->SetTitleSize(0.055); hFrS.GetYaxis()->SetTitleSize(0.055);
                                hFrS.GetXaxis()->SetLabelSize(0.045); hFrS.GetYaxis()->SetLabelSize(0.045);
                                hFrS.GetYaxis()->SetTitleOffset(1.15);
                                hFrS.Draw();
                                
                                TGraphErrors gTDs((int)xSig.size(), &xSig[0], &yTDpt[0], &exSig[0], &eyTDpt[0]);
                                gTDs.SetLineWidth(2); gTDs.SetLineColor(kBlack); gTDs.SetMarkerColor(kBlack);
                                gTDs.SetMarkerStyle(20); gTDs.SetMarkerSize(1.2);
                                gTDs.Draw("PE1 SAME");
                                
                                TGraphErrors gSMs((int)xSig.size(), &xSig[0], &ySMpt[0], &exSig[0], &eySMpt[0]);
                                gSMs.SetLineWidth(2); gSMs.SetLineColor(kBlue+1); gSMs.SetMarkerColor(kBlue+1);
                                gSMs.SetMarkerStyle(25); gSMs.SetMarkerSize(1.2);
                                gSMs.Draw("PE1 SAME");
                                
                                TLegend lgS(0.58, 0.58, 0.92, 0.72);
                                lgS.SetBorderSize(0); lgS.SetFillStyle(0); lgS.SetTextFont(42); lgS.SetTextSize(0.032);
                                lgS.AddEntry(&gTDs, "Data (Signal)", "ep");
                                lgS.AddEntry(&gSMs, TString::Format("%s Signal MC", mcPt.mcLabel.c_str()).Data(), "ep");
                                lgS.Draw();
                                
                                TLatex tTs; tTs.SetNDC(true); tTs.SetTextFont(42); tTs.SetTextAlign(23); tTs.SetTextSize(0.038);
                                tTs.DrawLatex(0.54, 0.975,
                                              TString::Format("Signal Overlay, <E_{T}^{iso}>, Au+Au & %s Embedded MC", mcPt.mcLabel.c_str()).Data());
                                
                                TLatex tIs; tIs.SetNDC(true); tIs.SetTextFont(42); tIs.SetTextAlign(13); tIs.SetTextSize(0.034);
                                tIs.DrawLatex(0.2, 0.88, trigDisplayLabel.c_str());
                                tIs.DrawLatex(0.2, 0.84, TString::Format("|v_{z}| < %d cm", kAA_VzCut).Data());
                                tIs.DrawLatex(0.2, 0.80, TString::Format("UE subtraction: %s", H.label.c_str()).Data());
                                tIs.DrawLatex(0.2, 0.76, TString::Format("#DeltaR_{cone} < %.1f", (kAA_IsoConeR == "isoR40") ? 0.4 : 0.3).Data());
                                tIs.DrawLatex(0.2, 0.72, TString::Format("Cent: %d-%d%%", cb.lo, cb.hi).Data());
                                
                                SaveCanvas(cSigPt, JoinPath(sigPtDir,
                                                            TString::Format("meanIsoEt_signal_vs_pT_%s.png", H.variant.c_str()).Data()));
                            }
                        }
                    }
                }
            }
            
            if (!generateISOpTcentOverlaysONLY)
            {
            // ── <E_T^iso> vs centrality (one per pT bin, this variant) ──
            for (int ipt = 0; ipt < kNPtBins; ++ipt)
            {
                const PtBin& b = PtBins()[ipt];
                
                vector<double> xCent, exCent, yCentPP, eyCentPP, yCentAA, eyCentAA;
                for (std::size_t ic = 0; ic < centBins.size(); ++ic)
                {
                    if (!vsCent_filled[ipt][ic]) continue;
                    const auto& cb = centBins[ic];
                    xCent.push_back(0.5 * (cb.lo + cb.hi));
                    exCent.push_back(0.5 * (cb.hi - cb.lo));
                    yCentPP.push_back(vsCent_yPP[ipt][ic]);
                    eyCentPP.push_back(vsCent_eyPP[ipt][ic]);
                    yCentAA.push_back(vsCent_yAA[ipt][ic]);
                    eyCentAA.push_back(vsCent_eyAA[ipt][ic]);
                }
                if (xCent.empty()) continue;
                
                TCanvas cVsCent(
                                TString::Format("c_meanIsoEtVsCent_%s_%s_%s",
                                                trigAA.c_str(), H.variant.c_str(), b.folder.c_str()).Data(),
                                "c_meanIsoEtVsCent", 900, 700
                                );
                ApplyCanvasMargins1D(cVsCent);
                cVsCent.cd();
                
                double yMin = std::numeric_limits<double>::max();
                double yMax = -std::numeric_limits<double>::max();
                for (std::size_t i = 0; i < xCent.size(); ++i)
                {
                    yMin = std::min(yMin, std::min(yCentPP[i] - eyCentPP[i], yCentAA[i] - eyCentAA[i]));
                    yMax = std::max(yMax, std::max(yCentPP[i] + eyCentPP[i], yCentAA[i] + eyCentAA[i]));
                }
                if (!std::isfinite(yMin) || !std::isfinite(yMax))
                {
                    yMin = 0.0;
                    yMax = 1.0;
                }
                const double pad = (yMax > yMin) ? (0.15 * (yMax - yMin)) : 0.25;
                
                const double centLo = centBins.front().lo;
                const double centHi = centBins.back().hi;
                
                TH1F hFrame(
                            TString::Format("hFrame_meanIsoEtVsCent_%s_%s_%s",
                                            trigAA.c_str(), H.variant.c_str(), b.folder.c_str()).Data(),
                            "", 100, centLo, centHi
                            );
                hFrame.SetDirectory(nullptr);
                hFrame.SetStats(0);
                hFrame.SetMinimum(std::max(0.0, yMin - pad));
                hFrame.SetMaximum(yMax + pad);
                hFrame.GetXaxis()->SetTitle("Centrality [%]");
                hFrame.GetYaxis()->SetTitle("<E_{T}^{iso}> [GeV]");
                hFrame.GetXaxis()->SetTitleSize(0.055);
                hFrame.GetYaxis()->SetTitleSize(0.055);
                hFrame.GetXaxis()->SetLabelSize(0.045);
                hFrame.GetYaxis()->SetLabelSize(0.045);
                hFrame.GetYaxis()->SetTitleOffset(1.15);
                hFrame.Draw();
                
                TGraphErrors gPP((int)xCent.size(), &xCent[0], &yCentPP[0], &exCent[0], &eyCentPP[0]);
                gPP.SetLineWidth(2);
                gPP.SetLineColor(kRed + 1);
                gPP.SetMarkerColor(kRed + 1);
                gPP.SetMarkerStyle(24);
                gPP.SetMarkerSize(1.2);
                gPP.Draw("PE1 SAME");
                
                TGraphErrors gAA((int)xCent.size(), &xCent[0], &yCentAA[0], &exCent[0], &eyCentAA[0]);
                gAA.SetLineWidth(2);
                gAA.SetLineColor(kBlack);
                gAA.SetMarkerColor(kBlack);
                gAA.SetMarkerStyle(20);
                gAA.SetMarkerSize(1.2);
                gAA.Draw("PE1 SAME");
                
                TLegend leg(0.56, 0.68, 0.92, 0.88);
                leg.SetBorderSize(0);
                leg.SetFillStyle(0);
                leg.SetTextFont(42);
                leg.SetTextSize(0.032);
                leg.AddEntry(&gPP, "pp data", "ep");
                leg.AddEntry(&gAA, TString::Format("AuAu data (%s)", H.label.c_str()).Data(), "ep");
                leg.Draw();
                
                TLatex tTitle;
                tTitle.SetNDC(true);
                tTitle.SetTextFont(42);
                tTitle.SetTextAlign(23);
                tTitle.SetTextSize(0.045);
                tTitle.DrawLatex(0.50, 0.955,
                                 TString::Format("<E_{T}^{iso}> vs Centrality  p_{T}^{#gamma} %d-%d GeV", b.lo, b.hi).Data());
                
                TLatex t;
                t.SetNDC(true);
                t.SetTextFont(42);
                t.SetTextAlign(13);
                t.SetTextSize(0.034);
                t.DrawLatex(0.16, 0.88, TString::Format("Trigger: %s", trigAA.c_str()).Data());
                t.DrawLatex(0.16, 0.84, TString::Format("|v_{z}| < %d cm", kAA_VzCut).Data());
                t.DrawLatex(0.16, 0.80, TString::Format("UE subtraction: %s", H.label.c_str()).Data());
                
                {
                    const string ppAuAuVsCentDir = JoinPath(variantDir, "ppAuAu_meanIsoEt_versusCentralityPerPtBin");
                    EnsureDir(ppAuAuVsCentDir);
                    SaveCanvas(cVsCent, JoinPath(ppAuAuVsCentDir,
                                                 TString::Format("meanIsoEt_pp_vs_auau_vs_cent_%s.png", b.folder.c_str()).Data()));
                }
            }
            
            // ── MC embedded <E_T^iso> vs centrality (one per pT bin, at variant level) ──
            {
                struct MCvsCentCfg {
                    const char* folderName;
                    string mcLabel;
                    const vector<vector<double>>& yMC;
                    const vector<vector<double>>& eyMC;
                    const vector<vector<bool>>&   filledMC;
                    const vector<vector<double>>& ySigMC;
                    const vector<vector<double>>& eySigMC;
                    const vector<vector<bool>>&   filledSigMC;
                };
                const string incJetCentName = InclusiveEmbeddedShortLabel().empty()
                ? "" : "Inclusive Jet " + InclusiveEmbeddedShortLabel();
                const string phoJetCentName = PhotonEmbeddedShortLabel().empty()
                ? "" : "Photon+Jet " + PhotonEmbeddedShortLabel();
                const vector<MCvsCentCfg> mcCentCfgs = {
                    {"AuAu_InclusiveJetEmbedded_meanIsoEt_versusCentralityPerPtBin",
                        incJetCentName,
                        vsCent_yIncMC, vsCent_eyIncMC, vsCent_filledIncMC,
                        vsCent_ySigIncMC, vsCent_eySigIncMC, vsCent_filledSigIncMC},
                    {"AuAu_PhotonJetEmbedded_meanIsoEt_versusCentralityPerPtBin",
                        phoJetCentName,
                        vsCent_yPhoMC, vsCent_eyPhoMC, vsCent_filledPhoMC,
                        vsCent_ySigPhoMC, vsCent_eySigPhoMC, vsCent_filledSigPhoMC}
                };
                for (const auto& mcC : mcCentCfgs)
                {
                    const string mcCentDir = JoinPath(variantDir, mcC.folderName);
                    EnsureDir(mcCentDir);
                    
                    for (int ipt = 0; ipt < kNPtBins; ++ipt)
                    {
                        const PtBin& bmc = PtBins()[ipt];
                        
                        vector<double> xC, exC, yAAc, eyAAc, yMCc, eyMCc;
                        for (std::size_t ic = 0; ic < centBins.size(); ++ic)
                        {
                            if (!mcC.filledMC[ipt][ic] || !vsCent_filled[ipt][ic]) continue;
                            const auto& cbmc = centBins[ic];
                            xC.push_back(0.5 * (cbmc.lo + cbmc.hi));
                            exC.push_back(0.5 * (cbmc.hi - cbmc.lo));
                            yAAc.push_back(vsCent_yAA[ipt][ic]);
                            eyAAc.push_back(vsCent_eyAA[ipt][ic]);
                            yMCc.push_back(mcC.yMC[ipt][ic]);
                            eyMCc.push_back(mcC.eyMC[ipt][ic]);
                        }
                        if (xC.empty()) continue;
                        
                        double yMinC = std::numeric_limits<double>::max();
                        double yMaxC = -std::numeric_limits<double>::max();
                        for (std::size_t i = 0; i < xC.size(); ++i)
                        {
                            yMinC = std::min(yMinC, std::min(yAAc[i] - eyAAc[i], yMCc[i] - eyMCc[i]));
                            yMaxC = std::max(yMaxC, std::max(yAAc[i] + eyAAc[i], yMCc[i] + eyMCc[i]));
                        }
                        if (!std::isfinite(yMinC) || !std::isfinite(yMaxC)) { yMinC = 0.0; yMaxC = 1.0; }
                        const double padC = (yMaxC > yMinC) ? (0.15 * (yMaxC - yMinC)) : 0.25;
                        
                        const double centLoC = centBins.front().lo;
                        const double centHiC = centBins.back().hi;
                        
                        TCanvas cMCcent(
                                        TString::Format("c_meanIsoMCvsCent_%s_%s_%s", mcC.folderName, H.variant.c_str(), bmc.folder.c_str()).Data(),
                                        "c_meanIsoMCvsCent", 900, 700);
                        ApplyCanvasMargins1D(cMCcent); cMCcent.cd();
                        
                        TH1F hFrC(
                                  TString::Format("hFr_meanIsoMCvsCent_%s_%s_%s", mcC.folderName, H.variant.c_str(), bmc.folder.c_str()).Data(),
                                  "", 100, centLoC, centHiC);
                        hFrC.SetDirectory(nullptr); hFrC.SetStats(0);
                        hFrC.SetMinimum(std::max(0.0, yMinC - padC));
                        hFrC.SetMaximum(yMaxC + padC);
                        hFrC.GetXaxis()->SetTitle("Centrality [%]");
                        hFrC.GetYaxis()->SetTitle("<E_{T}^{iso}> [GeV]");
                        hFrC.GetXaxis()->SetTitleSize(0.055); hFrC.GetYaxis()->SetTitleSize(0.055);
                        hFrC.GetXaxis()->SetLabelSize(0.045); hFrC.GetYaxis()->SetLabelSize(0.045);
                        hFrC.GetYaxis()->SetTitleOffset(1.15);
                        hFrC.Draw();
                        
                        TGraphErrors gAAc((int)xC.size(), &xC[0], &yAAc[0], &exC[0], &eyAAc[0]);
                        gAAc.SetLineWidth(2); gAAc.SetLineColor(kBlack); gAAc.SetMarkerColor(kBlack);
                        gAAc.SetMarkerStyle(20); gAAc.SetMarkerSize(1.2);
                        gAAc.Draw("PE1 SAME");
                        
                        TGraphErrors gMCc((int)xC.size(), &xC[0], &yMCc[0], &exC[0], &eyMCc[0]);
                        gMCc.SetLineWidth(2); gMCc.SetLineColor(kBlue+1); gMCc.SetMarkerColor(kBlue+1);
                        gMCc.SetMarkerStyle(25); gMCc.SetMarkerSize(1.2);
                        gMCc.Draw("PE1 SAME");
                        
                        TLegend lgC(0.58, 0.58, 0.92, 0.72);
                        lgC.SetBorderSize(0); lgC.SetFillStyle(0); lgC.SetTextFont(42); lgC.SetTextSize(0.032);
                        lgC.AddEntry(&gAAc, TString::Format("AuAu data (%s)", H.label.c_str()).Data(), "ep");
                        lgC.AddEntry(&gMCc, TString::Format("%s Embedded MC", mcC.mcLabel.c_str()).Data(), "ep");
                        lgC.Draw();
                        
                        TLatex tTc; tTc.SetNDC(true); tTc.SetTextFont(42); tTc.SetTextAlign(23); tTc.SetTextSize(0.038);
                        tTc.DrawLatex(0.50, 0.955,
                                      TString::Format("<E_{T}^{iso}> Overlay, Au+Au & %s Embedded MC", mcC.mcLabel.c_str()).Data());
                        
                        TLatex tCi; tCi.SetNDC(true); tCi.SetTextFont(42); tCi.SetTextAlign(13); tCi.SetTextSize(0.034);
                        tCi.DrawLatex(0.16, 0.88, trigDisplayLabel.c_str());
                        tCi.DrawLatex(0.16, 0.84, TString::Format("|v_{z}| < %d cm", kAA_VzCut).Data());
                        tCi.DrawLatex(0.16, 0.80, TString::Format("UE subtraction: %s", H.label.c_str()).Data());
                        tCi.DrawLatex(0.16, 0.76, TString::Format("#DeltaR_{cone} < %.1f", (kAA_IsoConeR == "isoR40") ? 0.4 : 0.3).Data());
                        tCi.DrawLatex(0.16, 0.72, TString::Format("p_{T}^{#gamma}: %d-%d GeV", bmc.lo, bmc.hi).Data());
                        
                        SaveCanvas(cMCcent, JoinPath(mcCentDir,
                                                     TString::Format("meanIsoEt_auau_vs_%s_vs_cent_%s.png",
                                                                     mcC.folderName, bmc.folder.c_str()).Data()));
                        
                        // ── Signal-only <E_T^iso> vs centrality ──
                        {
                            vector<double> xTDc, exTDc, yTDc, eyTDc, ySMCc2, eySMCc2;
                            for (std::size_t ic = 0; ic < centBins.size(); ++ic)
                            {
                                if (!mcC.filledSigMC[ipt][ic] && !vsCent_filledTightData[ipt][ic]) continue;
                                const auto& cbmc2 = centBins[ic];
                                xTDc.push_back(0.5 * (cbmc2.lo + cbmc2.hi));
                                exTDc.push_back(0.5 * (cbmc2.hi - cbmc2.lo));
                                yTDc.push_back(vsCent_filledTightData[ipt][ic] ? vsCent_yTightData[ipt][ic] : 0.0);
                                eyTDc.push_back(vsCent_filledTightData[ipt][ic] ? vsCent_eyTightData[ipt][ic] : 0.0);
                                ySMCc2.push_back(mcC.filledSigMC[ipt][ic] ? mcC.ySigMC[ipt][ic] : 0.0);
                                eySMCc2.push_back(mcC.filledSigMC[ipt][ic] ? mcC.eySigMC[ipt][ic] : 0.0);
                            }
                            if (!xTDc.empty())
                            {
                                const string sigCentDir = JoinPath(variantDir, string(mcC.folderName) + "_SIGNALonly");
                                EnsureDir(sigCentDir);
                                
                                double yMinS = std::numeric_limits<double>::max();
                                double yMaxS = -std::numeric_limits<double>::max();
                                for (std::size_t i = 0; i < xTDc.size(); ++i)
                                {
                                    yMinS = std::min(yMinS, std::min(yTDc[i] - eyTDc[i], ySMCc2[i] - eySMCc2[i]));
                                    yMaxS = std::max(yMaxS, std::max(yTDc[i] + eyTDc[i], ySMCc2[i] + eySMCc2[i]));
                                }
                                if (!std::isfinite(yMinS) || !std::isfinite(yMaxS)) { yMinS = 0.0; yMaxS = 1.0; }
                                const double padS = (yMaxS > yMinS) ? (0.15 * (yMaxS - yMinS)) : 0.25;
                                
                                const double centLoS = centBins.front().lo;
                                const double centHiS = centBins.back().hi;
                                
                                TCanvas cSigC(
                                              TString::Format("c_meanIsoSigVsCent_%s_%s_%s", mcC.folderName, H.variant.c_str(), bmc.folder.c_str()).Data(),
                                              "c_meanIsoSigVsCent", 900, 700);
                                ApplyCanvasMargins1D(cSigC); cSigC.cd();
                                
                                TH1F hFrSC(
                                           TString::Format("hFr_meanIsoSigVsCent_%s_%s_%s", mcC.folderName, H.variant.c_str(), bmc.folder.c_str()).Data(),
                                           "", 100, centLoS, centHiS);
                                hFrSC.SetDirectory(nullptr); hFrSC.SetStats(0);
                                hFrSC.SetMinimum(std::max(0.0, yMinS - padS));
                                hFrSC.SetMaximum(yMaxS + padS);
                                hFrSC.GetXaxis()->SetTitle("Centrality [%]");
                                hFrSC.GetYaxis()->SetTitle("<E_{T}^{iso}> [GeV]");
                                hFrSC.GetXaxis()->SetTitleSize(0.055); hFrSC.GetYaxis()->SetTitleSize(0.055);
                                hFrSC.GetXaxis()->SetLabelSize(0.045); hFrSC.GetYaxis()->SetLabelSize(0.045);
                                hFrSC.GetYaxis()->SetTitleOffset(1.15);
                                hFrSC.Draw();
                                
                                TGraphErrors gTDs((int)xTDc.size(), &xTDc[0], &yTDc[0], &exTDc[0], &eyTDc[0]);
                                gTDs.SetLineWidth(2); gTDs.SetLineColor(kBlack); gTDs.SetMarkerColor(kBlack);
                                gTDs.SetMarkerStyle(20); gTDs.SetMarkerSize(1.2);
                                gTDs.Draw("PE1 SAME");
                                
                                TGraphErrors gSMs((int)xTDc.size(), &xTDc[0], &ySMCc2[0], &exTDc[0], &eySMCc2[0]);
                                gSMs.SetLineWidth(2); gSMs.SetLineColor(kBlue+1); gSMs.SetMarkerColor(kBlue+1);
                                gSMs.SetMarkerStyle(25); gSMs.SetMarkerSize(1.2);
                                gSMs.Draw("PE1 SAME");
                                
                                TLegend lgS(0.58, 0.58, 0.92, 0.72);
                                lgS.SetBorderSize(0); lgS.SetFillStyle(0); lgS.SetTextFont(42); lgS.SetTextSize(0.032);
                                lgS.AddEntry(&gTDs, "Data (Signal)", "ep");
                                lgS.AddEntry(&gSMs, TString::Format("%s Signal MC", mcC.mcLabel.c_str()).Data(), "ep");
                                lgS.Draw();
                                
                                TLatex tTs; tTs.SetNDC(true); tTs.SetTextFont(42); tTs.SetTextAlign(23); tTs.SetTextSize(0.038);
                                tTs.DrawLatex(0.50, 0.955,
                                              TString::Format("Signal Overlay, <E_{T}^{iso}>, Au+Au & %s Embedded MC", mcC.mcLabel.c_str()).Data());
                                
                                TLatex tIs; tIs.SetNDC(true); tIs.SetTextFont(42); tIs.SetTextAlign(13); tIs.SetTextSize(0.034);
                                tIs.DrawLatex(0.16, 0.88, trigDisplayLabel.c_str());
                                tIs.DrawLatex(0.16, 0.84, TString::Format("|v_{z}| < %d cm", kAA_VzCut).Data());
                                tIs.DrawLatex(0.16, 0.80, TString::Format("UE subtraction: %s", H.label.c_str()).Data());
                                tIs.DrawLatex(0.16, 0.76, TString::Format("#DeltaR_{cone} < %.1f", (kAA_IsoConeR == "isoR40") ? 0.4 : 0.3).Data());
                                tIs.DrawLatex(0.16, 0.72, TString::Format("p_{T}^{#gamma}: %d-%d GeV", bmc.lo, bmc.hi).Data());
                                
                                SaveCanvas(cSigC, JoinPath(sigCentDir,
                                                           TString::Format("meanIsoEt_signal_vs_cent_%s.png", bmc.folder.c_str()).Data()));
                            }
                        }
                    }
                }
            }
            
            } // end if (!generateISOpTcentOverlaysONLY) — vs-centrality & MC sections
            
            // ── Centrality overlay of E_T^iso (unnormalized counts) per pT bin, this variant ──
            {
                const int centOvColors[] = {kRed+1, kBlue+1, kGreen+2, kMagenta+1, kOrange+1,
                    kCyan+1, kYellow+2, kViolet+1};
                for (int ipt = 0; ipt < kNPtBins; ++ipt)
                {
                    const PtBin& b = PtBins()[ipt];
                    
                    vector<TH1*> hCents;
                    vector<string> cLabels;
                    double yMaxCO = 0.0;
                    
                    for (std::size_t ic = 0; ic < centBins.size(); ++ic)
                    {
                        const auto& cb = centBins[ic];
                        const string hAAName = "h_Eiso" + b.suffix + cb.suffix;
                        TH1* hAAsrc = dynamic_cast<TH1*>(aaTop->Get(hAAName.c_str()));
                        if (!hAAsrc) continue;
                        
                        TH1* hAA = CloneTH1(hAAsrc,
                                            TString::Format("hAA_centOvlay_%s_%s_%s_%s",
                                                            trigAA.c_str(), H.variant.c_str(), b.folder.c_str(), cb.folder.c_str()).Data());
                        if (!hAA) continue;
                        
                        hAA->Rebin(10);
                        EnsureSumw2(hAA);
                        const int ci = (ic < 8) ? centOvColors[ic] : kBlack;
                        hAA->SetLineColor(ci);
                        hAA->SetMarkerColor(ci);
                        hAA->SetMarkerStyle(20);
                        hAA->SetMarkerSize(0.9);
                        hAA->SetLineWidth(2);
                        hAA->SetFillStyle(0);
                        
                        yMaxCO = std::max(yMaxCO, hAA->GetMaximum());
                        hCents.push_back(hAA);
                        cLabels.push_back(TString::Format("%d-%d%%", cb.lo, cb.hi).Data());
                    }
                    
                    if (hCents.empty()) continue;
                    
                    // --- Collect Gaussian means per centrality for subpanel ---
                    vector<double> subX, subY, subEY;
                    for (std::size_t ih = 0; ih < hCents.size(); ++ih)
                    {
                        double gM, gS, gME, gSE;
                        if (FitGaussianIterative(hCents[ih], gM, gS, gME, gSE))
                        {
                            subX.push_back((double)ih);
                            subY.push_back(gM);
                            subEY.push_back(gME);
                        }
                    }
                    
                    TCanvas cCO(
                                TString::Format("c_isoCentOverlay_%s_%s_%s",
                                                trigAA.c_str(), H.variant.c_str(), b.folder.c_str()).Data(),
                                "c_isoCentOverlay", 900, 850
                                );
                    cCO.cd();
                    
                    // Upper pad: histogram overlay
                    TPad* padUp = new TPad(TString::Format("padUp_%s_%s", H.variant.c_str(), b.folder.c_str()).Data(),
                                           "padUp", 0.0, 0.25, 1.0, 1.0);
                    padUp->SetBottomMargin(0.02);
                    padUp->SetLeftMargin(0.14);
                    padUp->SetRightMargin(0.04);
                    padUp->SetTopMargin(0.06);
                    padUp->Draw();
                    padUp->cd();
                    
                    hCents[0]->SetTitle("");
                    hCents[0]->GetXaxis()->SetTitle("");
                    hCents[0]->GetXaxis()->SetLabelSize(0.0);
                    hCents[0]->GetYaxis()->SetTitle("Counts");
                    hCents[0]->GetXaxis()->SetTitleSize(0.055);
                    hCents[0]->GetYaxis()->SetTitleSize(0.060);
                    hCents[0]->GetYaxis()->SetLabelSize(0.050);
                    hCents[0]->GetYaxis()->SetTitleOffset(1.05);
                    hCents[0]->GetXaxis()->SetRangeUser(-10.0, 50.0);
                    hCents[0]->SetMinimum(0.0);
                    {
                        const double yScale = (H.variant == "noSub") ? 1.35
                                            : (H.variant == "baseVariant") ? 1.15
                                            : 1.25;
                        hCents[0]->SetMaximum((yMaxCO > 0.0) ? (yScale * yMaxCO) : 1.0);
                    }
                    hCents[0]->Draw("E1");
                    for (std::size_t ih = 1; ih < hCents.size(); ++ih) hCents[ih]->Draw("E1 SAME");
                    
                    // Gaussian fit curves per centrality bin
                    vector<TF1*> gaussCurvesC;
                    for (std::size_t ih = 0; ih < hCents.size(); ++ih)
                    {
                        TF1* fc = DrawGaussFitCurve(hCents[ih], hCents[ih]->GetLineColor());
                        if (fc) gaussCurvesC.push_back(fc);
                    }
                    
                    const bool legTopLeft = (H.variant == "noSub");
                    TLegend legCO(legTopLeft ? 0.19 : 0.60,
                                  legTopLeft ? 0.68 : 0.40,
                                  legTopLeft ? 0.45 : 0.90,
                                  legTopLeft ? 0.90 : 0.62);
                    legCO.SetColumnSeparation(-0.05);
                    legCO.SetBorderSize(0);
                    legCO.SetFillStyle(0);
                    legCO.SetTextFont(42);
                    legCO.SetTextSize(0.035);
                    legCO.SetNColumns(2);
                    for (std::size_t ih = 0; ih < hCents.size(); ++ih)
                        legCO.AddEntry(hCents[ih], cLabels[ih].c_str(), "ep");
                    legCO.Draw();
                    
                    TLatex tTitleCO;
                    tTitleCO.SetNDC(true);
                    tTitleCO.SetTextFont(42);
                    tTitleCO.SetTextAlign(23);
                    tTitleCO.SetTextSize(0.045);
                    tTitleCO.DrawLatex(0.50, 0.97,
                                       TString::Format("Au+Au (counts), centrality overlays, p_{T}^{#gamma} = %d-%d GeV", b.lo, b.hi).Data());
                    
                    TLatex tInfoCO;
                    tInfoCO.SetNDC(true);
                    tInfoCO.SetTextFont(42);
                    tInfoCO.SetTextAlign(33);
                    tInfoCO.SetTextSize(0.035);
                    tInfoCO.DrawLatex(0.92, 0.90,
                                      forEmbeddedSim
                                      ? (embeddedMode == 2 ? activeEmbeddedSimFolder.substr(0, activeEmbeddedSimFolder.find("_SIM")).c_str()
                                         : SimSampleLabel(CurrentSimSample()).c_str())
                                      : trigDisplayLabel.c_str());
                    tInfoCO.DrawLatex(0.92, 0.86, H.label.c_str());
                    {
                        const double coneRValCO = (kAA_IsoConeR == "isoR40") ? 0.4 : 0.3;
                        tInfoCO.DrawLatex(0.92, 0.82, TString::Format("#DeltaR_{cone} < %.1f", coneRValCO).Data());
                    }
                    {
                        TLatex tSph;
                        tSph.SetNDC(true);
                        tSph.SetTextFont(42);
                        tSph.SetTextAlign(33);
                        tSph.SetTextSize(0.052);
                        const double sphY1 = (H.variant == "noSub") ? 0.58 : 0.34;
                        tSph.DrawLatex(0.92, sphY1,        "#bf{sPHENIX} #it{Internal}");
                        tSph.SetTextSize(0.042);
                        tSph.DrawLatex(0.92, sphY1 - 0.06, "Au+Au  #sqrt{s_{NN}} = 200 GeV");
                    }
                    
                    // Lower pad: Gaussian mean vs centrality bin index
                    cCO.cd();
                    TPad* padLo = new TPad(TString::Format("padLo_%s_%s", H.variant.c_str(), b.folder.c_str()).Data(),
                                           "padLo", 0.0, 0.0, 1.0, 0.25);
                    padLo->SetTopMargin(0.02);
                    padLo->SetBottomMargin(0.35);
                    padLo->SetLeftMargin(0.14);
                    padLo->SetRightMargin(0.04);
                    padLo->Draw();
                    padLo->cd();
                    
                    if (!subX.empty())
                    {
                        double yLoSub = 1e30, yHiSub = -1e30;
                        for (std::size_t is = 0; is < subY.size(); ++is)
                        {
                            yLoSub = std::min(yLoSub, subY[is] - subEY[is]);
                            yHiSub = std::max(yHiSub, subY[is] + subEY[is]);
                        }
                        const double subPad = (yHiSub > yLoSub) ? 0.25 * (yHiSub - yLoSub) : 1.0;
                        
                        TH1F hFrSub(TString::Format("hFrSub_centOv_%s_%s", H.variant.c_str(), b.folder.c_str()).Data(),
                                    "", (int)centBins.size(), -0.5, (double)centBins.size() - 0.5);
                        hFrSub.SetDirectory(nullptr);
                        hFrSub.SetStats(0);
                        hFrSub.SetMinimum(yLoSub - subPad);
                        hFrSub.SetMaximum(yHiSub + subPad);
                        hFrSub.GetYaxis()->SetTitle("Gauss #mu");
                        hFrSub.GetYaxis()->SetTitleSize(0.16);
                        hFrSub.GetYaxis()->SetTitleOffset(0.35);
                        hFrSub.GetYaxis()->SetLabelSize(0.13);
                        hFrSub.GetYaxis()->SetNdivisions(505);
                        hFrSub.GetXaxis()->SetTitle("Centrality Bin");
                        hFrSub.GetXaxis()->SetTitleSize(0.16);
                        hFrSub.GetXaxis()->SetTitleOffset(0.95);
                        hFrSub.GetXaxis()->SetLabelSize(0.13);
                        // Label bins with centrality %
                        for (std::size_t ic = 0; ic < centBins.size(); ++ic)
                            hFrSub.GetXaxis()->SetBinLabel((int)ic + 1,
                                TString::Format("%d-%d%%", centBins[ic].lo, centBins[ic].hi).Data());
                        hFrSub.Draw();
                        
                        vector<double> exSub(subX.size(), 0.0);
                        TGraphErrors gSub((int)subX.size(), &subX[0], &subY[0], &exSub[0], &subEY[0]);
                        gSub.SetMarkerStyle(20);
                        gSub.SetMarkerSize(1.0);
                        gSub.SetMarkerColor(kBlack);
                        gSub.SetLineColor(kBlack);
                        gSub.SetLineWidth(2);
                        gSub.Draw("PE1 SAME");
                    }
                    
                    SaveCanvas(cCO, JoinPath(variantCentralitySummaryDir,
                                             TString::Format("isoCentOverlay_counts_%s.png", b.folder.c_str()).Data()));
                    
                    for (auto* f : gaussCurvesC) delete f;
                    for (auto* h : hCents) delete h;
                }
            }
            
            // ── pT overlay of E_T^iso (unnormalized counts) per centrality bin, this variant ──
            {
                const int ptOvColors[] = {kRed+1, kBlue+1, kGreen+2, kMagenta+1, kOrange+1,
                    kCyan+1, kViolet+1, kYellow+2, kSpring+2, kAzure+1,
                    kTeal+1, kPink+1};
                for (std::size_t ic = 0; ic < centBins.size(); ++ic)
                {
                    const auto& cb = centBins[ic];
                    const string centSubDir = JoinPath(variantDir, cb.folder);
                    EnsureDir(centSubDir);
                    
                    vector<TH1*> hPts;
                    vector<string> ptLabels;
                    double yMaxPO = 0.0;
                    
                    for (int ipt = 0; ipt < kNPtBins; ++ipt)
                    {
                        const PtBin& b = PtBins()[ipt];
                        const string hAAName = "h_Eiso" + b.suffix + cb.suffix;
                        TH1* hAAsrc = dynamic_cast<TH1*>(aaTop->Get(hAAName.c_str()));
                        if (!hAAsrc) continue;
                        
                        TH1* hAA = CloneTH1(hAAsrc,
                                            TString::Format("hAA_ptOvlay_%s_%s_%s_%s",
                                                            trigAA.c_str(), H.variant.c_str(), cb.folder.c_str(), b.folder.c_str()).Data());
                        if (!hAA) continue;
                        
                        hAA->Rebin(10);
                        EnsureSumw2(hAA);
                        const int ci = (ipt < 12) ? ptOvColors[ipt] : kBlack;
                        hAA->SetLineColor(ci);
                        hAA->SetMarkerColor(ci);
                        hAA->SetMarkerStyle(20);
                        hAA->SetMarkerSize(0.9);
                        hAA->SetLineWidth(2);
                        hAA->SetFillStyle(0);
                        
                        yMaxPO = std::max(yMaxPO, hAA->GetMaximum());
                        hPts.push_back(hAA);
                        ptLabels.push_back(TString::Format("%d-%d GeV", b.lo, b.hi).Data());
                    }
                    
                    if (hPts.empty()) continue;
                    
                    // --- Collect Gaussian means per pT bin for subpanel ---
                    vector<double> subPtX, subPtY, subPtEY;
                    {
                        int ptIdx = 0;
                        for (int ipt2 = 0; ipt2 < kNPtBins; ++ipt2)
                        {
                            const PtBin& bp = PtBins()[ipt2];
                            const string hN2 = "h_Eiso" + bp.suffix + cb.suffix;
                            TH1* hSrc2 = dynamic_cast<TH1*>(aaTop->Get(hN2.c_str()));
                            if (!hSrc2) continue;
                            TH1* hTmp2 = CloneTH1(hSrc2, TString::Format("hTmp_ptSub_%s_%s_%d", H.variant.c_str(), cb.folder.c_str(), ipt2).Data());
                            if (!hTmp2) { ++ptIdx; continue; }
                            hTmp2->Rebin(10); EnsureSumw2(hTmp2);
                            double gM2, gS2, gME2, gSE2;
                            if (FitGaussianIterative(hTmp2, gM2, gS2, gME2, gSE2))
                            {
                                subPtX.push_back(0.5 * (kPtEdges[(std::size_t)ipt2] + kPtEdges[(std::size_t)ipt2 + 1]));
                                subPtY.push_back(gM2);
                                subPtEY.push_back(gME2);
                            }
                            delete hTmp2;
                            ++ptIdx;
                        }
                    }
                    
                    TCanvas cPO(
                                TString::Format("c_isoPtOverlay_%s_%s_%s",
                                                trigAA.c_str(), H.variant.c_str(), cb.folder.c_str()).Data(),
                                "c_isoPtOverlay", 900, 850
                                );
                    cPO.cd();
                    
                    // Upper pad: histogram overlay
                    TPad* padUpPO = new TPad(TString::Format("padUpPO_%s_%s", H.variant.c_str(), cb.folder.c_str()).Data(),
                                             "padUpPO", 0.0, 0.25, 1.0, 1.0);
                    padUpPO->SetBottomMargin(0.02);
                    padUpPO->SetLeftMargin(0.14);
                    padUpPO->SetRightMargin(0.04);
                    padUpPO->SetTopMargin(0.06);
                    padUpPO->Draw();
                    padUpPO->cd();
                    
                    hPts[0]->SetTitle("");
                    hPts[0]->GetXaxis()->SetTitle("");
                    hPts[0]->GetXaxis()->SetLabelSize(0.0);
                    hPts[0]->GetYaxis()->SetTitle("Counts");
                    hPts[0]->GetXaxis()->SetTitleSize(0.055);
                    hPts[0]->GetYaxis()->SetTitleSize(0.060);
                    hPts[0]->GetYaxis()->SetLabelSize(0.050);
                    hPts[0]->GetYaxis()->SetTitleOffset(1.05);
                    hPts[0]->SetMinimum(0.0);
                    hPts[0]->SetMaximum((yMaxPO > 0.0) ? (1.25 * yMaxPO) : 1.0);
                    hPts[0]->Draw("E1");
                    for (std::size_t ih = 1; ih < hPts.size(); ++ih) hPts[ih]->Draw("E1 SAME");
                    
                    // Gaussian fit curves per pT bin
                    vector<TF1*> gaussCurvesP;
                    for (std::size_t ih = 0; ih < hPts.size(); ++ih)
                    {
                        TF1* fp = DrawGaussFitCurve(hPts[ih], hPts[ih]->GetLineColor());
                        if (fp) gaussCurvesP.push_back(fp);
                    }
                    
                    TLegend legPO(0.72, 0.48, 0.92, 0.90);
                    legPO.SetBorderSize(0);
                    legPO.SetFillStyle(0);
                    legPO.SetTextFont(42);
                    legPO.SetTextSize(0.032);
                    for (std::size_t ih = 0; ih < hPts.size(); ++ih)
                        legPO.AddEntry(hPts[ih], ptLabels[ih].c_str(), "ep");
                    legPO.Draw();
                    
                    TLatex tTitlePO;
                    tTitlePO.SetNDC(true);
                    tTitlePO.SetTextFont(42);
                    tTitlePO.SetTextAlign(23);
                    tTitlePO.SetTextSize(0.045);
                    tTitlePO.DrawLatex(0.50, 0.97,
                                       TString::Format("Au+Au (counts), p_{T}^{#gamma} overlays, cent = %d-%d%%", cb.lo, cb.hi).Data());
                    
                    TLatex tInfoPO;
                    tInfoPO.SetNDC(true);
                    tInfoPO.SetTextFont(42);
                    tInfoPO.SetTextAlign(33);
                    tInfoPO.SetTextSize(0.035);
                    tInfoPO.DrawLatex(0.55, 0.90,
                                      forEmbeddedSim
                                      ? (embeddedMode == 2 ? activeEmbeddedSimFolder.substr(0, activeEmbeddedSimFolder.find("_SIM")).c_str()
                                         : SimSampleLabel(CurrentSimSample()).c_str())
                                      : trigDisplayLabel.c_str());
                    tInfoPO.DrawLatex(0.55, 0.86, H.label.c_str());
                    {
                        const double coneRValPO = (kAA_IsoConeR == "isoR40") ? 0.4 : 0.3;
                        tInfoPO.DrawLatex(0.55, 0.82, TString::Format("#DeltaR_{cone} < %.1f", coneRValPO).Data());
                    }
                    {
                        TLatex tSphPO;
                        tSphPO.SetNDC(true);
                        tSphPO.SetTextFont(42);
                        tSphPO.SetTextAlign(33);
                        tSphPO.SetTextSize(0.052);
                        tSphPO.DrawLatex(0.55, 0.45, "#bf{sPHENIX} #it{Internal}");
                        tSphPO.SetTextSize(0.042);
                        tSphPO.DrawLatex(0.55, 0.39, "Au+Au  #sqrt{s_{NN}} = 200 GeV");
                    }
                    
                    // Lower pad: Gaussian mean vs pT
                    cPO.cd();
                    TPad* padLoPO = new TPad(TString::Format("padLoPO_%s_%s", H.variant.c_str(), cb.folder.c_str()).Data(),
                                             "padLoPO", 0.0, 0.0, 1.0, 0.25);
                    padLoPO->SetTopMargin(0.02);
                    padLoPO->SetBottomMargin(0.35);
                    padLoPO->SetLeftMargin(0.14);
                    padLoPO->SetRightMargin(0.04);
                    padLoPO->Draw();
                    padLoPO->cd();
                    
                    if (!subPtX.empty())
                    {
                        double yLoPt = 1e30, yHiPt = -1e30;
                        for (std::size_t is = 0; is < subPtY.size(); ++is)
                        {
                            yLoPt = std::min(yLoPt, subPtY[is] - subPtEY[is]);
                            yHiPt = std::max(yHiPt, subPtY[is] + subPtEY[is]);
                        }
                        const double subPadPt = (yHiPt > yLoPt) ? 0.25 * (yHiPt - yLoPt) : 1.0;
                        
                        TH1F hFrSubPt(TString::Format("hFrSub_ptOv_%s_%s", H.variant.c_str(), cb.folder.c_str()).Data(),
                                      "", 100, kPtEdges.front(), kPtEdges.back());
                        hFrSubPt.SetDirectory(nullptr);
                        hFrSubPt.SetStats(0);
                        hFrSubPt.SetMinimum(yLoPt - subPadPt);
                        hFrSubPt.SetMaximum(yHiPt + subPadPt);
                        hFrSubPt.GetYaxis()->SetTitle("Gauss #mu");
                        hFrSubPt.GetYaxis()->SetTitleSize(0.16);
                        hFrSubPt.GetYaxis()->SetTitleOffset(0.35);
                        hFrSubPt.GetYaxis()->SetLabelSize(0.13);
                        hFrSubPt.GetYaxis()->SetNdivisions(505);
                        hFrSubPt.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                        hFrSubPt.GetXaxis()->SetTitleSize(0.16);
                        hFrSubPt.GetXaxis()->SetTitleOffset(0.95);
                        hFrSubPt.GetXaxis()->SetLabelSize(0.13);
                        hFrSubPt.Draw();
                        
                        vector<double> exSubPt(subPtX.size(), 0.0);
                        TGraphErrors gSubPt((int)subPtX.size(), &subPtX[0], &subPtY[0], &exSubPt[0], &subPtEY[0]);
                        gSubPt.SetMarkerStyle(20);
                        gSubPt.SetMarkerSize(1.0);
                        gSubPt.SetMarkerColor(kBlack);
                        gSubPt.SetLineColor(kBlack);
                        gSubPt.SetLineWidth(2);
                        gSubPt.Draw("PE1 SAME");
                    }
                    
                    SaveCanvas(cPO, JoinPath(centSubDir, "isoPtOverlay_counts.png"));
                    SaveCanvas(cPO, JoinPath(variantPtSummaryDir,
                                             TString::Format("isoPtOverlay_counts_%s.png", cb.folder.c_str()).Data()));
                    
                    for (auto* f : gaussCurvesP) delete f;
                    for (auto* h : hPts) delete h;
                }
            }
            
            if (fIncMCvar) { fIncMCvar->Close(); delete fIncMCvar; fIncMCvar = nullptr; }
            if (fPhoMCvar) { fPhoMCvar->Close(); delete fPhoMCvar; fPhoMCvar = nullptr; }
            
            cout << ANSI_DIM
            << "  -> UE comparison overlays written under:\n"
            << "     " << variantDir << "/\n"
            << ANSI_RESET;
        }
        
        
        const int variantColors[4]  = {kBlack, kBlue + 1, kOrange + 7, kGreen + 2};
        const int variantMarkers[4] = {20, 20, 20, 20};
        
        if (!generateISOpTcentOverlaysONLY)
        {
            for (std::size_t ic = 0; ic < centBins.size(); ++ic)
            {
                const auto& cb = centBins[ic];
                const string centOverlayDir = JoinPath(perVariantOverlayBase, cb.folder);
                EnsureDir(centOverlayDir);
                
                for (int ipt = 0; ipt < kNPtBins; ++ipt)
                {
                    const PtBin& b = PtBins()[ipt];
                    const string ptOverlayDir = JoinPath(centOverlayDir, b.folder);
                    EnsureDir(ptOverlayDir);
                    
                    vector<TH1*> hVars;
                    vector<std::size_t> hVarIndices;
                    double yMax = 0.0;
                    
                    for (std::size_t iv = 0; iv < handles.size(); ++iv)
                    {
                        auto& H = handles[iv];
                        if (!H.file) continue;
                        
                        TDirectory* aaTop = H.file->GetDirectory(trigAA.c_str());
                        if (!aaTop) continue;
                        
                        const string hAAName = "h_Eiso" + b.suffix + cb.suffix;
                        TH1* hAAsrc = dynamic_cast<TH1*>(aaTop->Get(hAAName.c_str()));
                        if (!hAAsrc) continue;
                        
                        TH1* hAA = CloneTH1(
                                            hAAsrc,
                                            TString::Format("hAA_allVariantOverlay_%s_%s_%s_%s",
                                                            trigAA.c_str(), H.variant.c_str(), cb.folder.c_str(), b.folder.c_str()).Data()
                                            );
                        if (!hAA) continue;
                        
                        hAA->Rebin(10);
                        EnsureSumw2(hAA);
                        const double integral = hAA->Integral(0, hAA->GetNbinsX() + 1);
                        if (integral > 0.0) hAA->Scale(1.0 / integral);
                        
                        hAA->SetLineColor((iv < 4) ? variantColors[iv] : kBlack);
                        hAA->SetMarkerColor((iv < 4) ? variantColors[iv] : kBlack);
                        hAA->SetMarkerStyle((iv < 4) ? variantMarkers[iv] : 20);
                        hAA->SetMarkerSize(1.1);
                        hAA->SetLineWidth(2);
                        hAA->SetFillStyle(0);
                        
                        yMax = std::max(yMax, hAA->GetMaximum());
                        hVars.push_back(hAA);
                        hVarIndices.push_back(iv);
                    }
                    
                    if (hVars.empty()) continue;
                    
                    TCanvas cOverlay(
                                     TString::Format("c_isoOverlayAllVariants_%s_%s_%s",
                                                     trigAA.c_str(), cb.folder.c_str(), b.folder.c_str()).Data(),
                                     "c_isoOverlayAllVariants", 900, 700
                                     );
                    ApplyCanvasMargins1D(cOverlay);
                    cOverlay.cd();
                    
                    hVars[0]->SetTitle("");
                    hVars[0]->GetXaxis()->SetTitle("E_{T}^{iso} [GeV]");
                    hVars[0]->GetYaxis()->SetTitle("Normalized to unit area");
                    hVars[0]->GetXaxis()->SetTitleSize(0.055);
                    hVars[0]->GetYaxis()->SetTitleSize(0.055);
                    hVars[0]->GetXaxis()->SetLabelSize(0.045);
                    hVars[0]->GetYaxis()->SetLabelSize(0.045);
                    hVars[0]->GetYaxis()->SetTitleOffset(1.15);
                    hVars[0]->SetMinimum(0.0);
                    hVars[0]->SetMaximum((yMax > 0.0) ? (1.25 * yMax) : 1.0);
                    hVars[0]->Draw("E1");
                    for (std::size_t ih = 1; ih < hVars.size(); ++ih) hVars[ih]->Draw("E1 SAME");
                    
                    TLegend leg(0.56, 0.62, 0.92, 0.88);
                    leg.SetBorderSize(0);
                    leg.SetFillStyle(0);
                    leg.SetTextFont(42);
                    leg.SetTextSize(0.032);
                    for (std::size_t ih = 0; ih < hVars.size(); ++ih)
                    {
                        leg.AddEntry(hVars[ih], handles[hVarIndices[ih]].label.c_str(), "ep");
                    }
                    leg.Draw();
                    
                    TLatex tTitle;
                    tTitle.SetNDC(true);
                    tTitle.SetTextFont(42);
                    tTitle.SetTextAlign(23);
                    tTitle.SetTextSize(0.045);
                    tTitle.DrawLatex(0.50, 0.98,
                                     TString::Format("Isolation overlay, all UE variants, %d-%d%% Cent AuAu", cb.lo, cb.hi).Data());
                    
                    const double coneRValOv = (kAA_IsoConeR == "isoR40") ? 0.4 : 0.3;
                    const int    isoEtMaxOv = (kAA_IsoMode == "fixedIso5GeV") ? 5 : -1;
                    
                    TLatex t;
                    t.SetNDC(true);
                    t.SetTextFont(42);
                    t.SetTextAlign(13);
                    t.SetTextSize(0.028);
                    t.DrawLatex(0.18, 0.89, "Trigger = Photon 10 GeV + MBD NS #geq 2, vtx < 150 cm");
                    t.DrawLatex(0.18, 0.84, TString::Format("#DeltaR_{cone} < %.1f", coneRValOv).Data());
                    if (isoEtMaxOv > 0)
                        t.DrawLatex(0.18, 0.8, TString::Format("E_{T}^{iso} < %d GeV", isoEtMaxOv).Data());
                    t.DrawLatex(0.18, 0.76, TString::Format("p_{T}^{#gamma}: %d-%d GeV", b.lo, b.hi).Data());
                    
                    for (std::size_t ih = 0; ih < hVars.size(); ++ih)
                    {
                        if (hVarIndices[ih] == 1) continue;  // skip baseVariant mean line
                        const double meanX = hVars[ih]->GetMean();
                        TLine lMean(meanX, 0.0, meanX, hVars[0]->GetMaximum());
                        lMean.SetLineColor(hVars[ih]->GetLineColor());
                        lMean.SetLineStyle(2);
                        lMean.SetLineWidth(2);
                        lMean.DrawClone();
                    }
                    
                    SaveCanvas(cOverlay, JoinPath(ptOverlayDir, "isolationOverlay_allVariants.png"));
                    
                    for (auto* h : hVars) delete h;
                    
                    // ------ noSub + baseVariant + variantA only overlay ------
                    {
                        vector<TH1*> h3Vars;
                        vector<std::size_t> h3Indices;
                        double yMax3 = 0.0;
                        
                        for (std::size_t iv : {std::size_t(0), std::size_t(1), std::size_t(2)})
                        {
                            if (iv >= handles.size()) continue;
                            auto& H = handles[iv];
                            if (!H.file) continue;
                            
                            TDirectory* aaTop3 = H.file->GetDirectory(trigAA.c_str());
                            if (!aaTop3) continue;
                            
                            const string hAAName3 = "h_Eiso" + b.suffix + cb.suffix;
                            TH1* hAAsrc3 = dynamic_cast<TH1*>(aaTop3->Get(hAAName3.c_str()));
                            if (!hAAsrc3) continue;
                            
                            TH1* hAA3 = CloneTH1(
                                                 hAAsrc3,
                                                 TString::Format("hAA_3varOverlay_%s_%s_%s_%s",
                                                                 trigAA.c_str(), H.variant.c_str(), cb.folder.c_str(), b.folder.c_str()).Data()
                                                 );
                            if (!hAA3) continue;
                            
                            hAA3->Rebin(10);
                            EnsureSumw2(hAA3);
                            const double integral3 = hAA3->Integral(0, hAA3->GetNbinsX() + 1);
                            if (integral3 > 0.0) hAA3->Scale(1.0 / integral3);
                            
                            hAA3->SetLineColor(variantColors[iv]);
                            hAA3->SetMarkerColor(variantColors[iv]);
                            hAA3->SetMarkerStyle(variantMarkers[iv]);
                            hAA3->SetMarkerSize(1.1);
                            hAA3->SetLineWidth(2);
                            hAA3->SetFillStyle(0);
                            
                            yMax3 = std::max(yMax3, hAA3->GetMaximum());
                            h3Vars.push_back(hAA3);
                            h3Indices.push_back(iv);
                        }
                        
                        if (!h3Vars.empty())
                        {
                            TCanvas c3V("c_isoOverlay_noSub_baseVariant_variantA",
                                        "c_isoOverlay_noSub_baseVariant_variantA", 900, 700);
                            ApplyCanvasMargins1D(c3V);
                            c3V.cd();
                            
                            h3Vars[0]->SetTitle("");
                            h3Vars[0]->GetXaxis()->SetTitle("E_{T}^{iso} [GeV]");
                            h3Vars[0]->GetYaxis()->SetTitle("Normalized to unit area");
                            h3Vars[0]->GetXaxis()->SetTitleSize(0.055);
                            h3Vars[0]->GetYaxis()->SetTitleSize(0.055);
                            h3Vars[0]->GetXaxis()->SetLabelSize(0.045);
                            h3Vars[0]->GetYaxis()->SetLabelSize(0.045);
                            h3Vars[0]->GetYaxis()->SetTitleOffset(1.15);
                            h3Vars[0]->SetMinimum(0.0);
                            h3Vars[0]->SetMaximum((yMax3 > 0.0) ? (1.25 * yMax3) : 1.0);
                            h3Vars[0]->Draw("E1");
                            for (std::size_t ih = 1; ih < h3Vars.size(); ++ih) h3Vars[ih]->Draw("E1 SAME");
                            
                            TLegend leg3V(0.56, 0.68, 0.92, 0.88);
                            leg3V.SetBorderSize(0);
                            leg3V.SetFillStyle(0);
                            leg3V.SetTextFont(42);
                            leg3V.SetTextSize(0.032);
                            for (std::size_t ih = 0; ih < h3Vars.size(); ++ih)
                                leg3V.AddEntry(h3Vars[ih], handles[h3Indices[ih]].label.c_str(), "ep");
                            leg3V.Draw();
                            
                            TLatex tTitle3V;
                            tTitle3V.SetNDC(true);
                            tTitle3V.SetTextFont(42);
                            tTitle3V.SetTextAlign(23);
                            tTitle3V.SetTextSize(0.045);
                            tTitle3V.DrawLatex(0.50, 0.98,
                                               TString::Format("Isolation overlay, noSub + baseVariant + variantA, %d-%d%% Cent AuAu", cb.lo, cb.hi).Data());
                            
                            const double coneRVal3V = (kAA_IsoConeR == "isoR40") ? 0.4 : 0.3;
                            const int    isoEtMax3V = (kAA_IsoMode == "fixedIso5GeV") ? 5 : -1;
                            
                            TLatex t3V;
                            t3V.SetNDC(true);
                            t3V.SetTextFont(42);
                            t3V.SetTextAlign(13);
                            t3V.SetTextSize(0.028);
                            t3V.DrawLatex(0.18, 0.89, "Trigger = Photon 10 GeV + MBD NS #geq 2, vtx < 150 cm");
                            t3V.DrawLatex(0.18, 0.84, TString::Format("#DeltaR_{cone} < %.1f", coneRVal3V).Data());
                            if (isoEtMax3V > 0)
                                t3V.DrawLatex(0.18, 0.8, TString::Format("E_{T}^{iso} < %d GeV", isoEtMax3V).Data());
                            t3V.DrawLatex(0.18, 0.76, TString::Format("p_{T}^{#gamma}: %d-%d GeV", b.lo, b.hi).Data());
                            
                            SaveCanvas(c3V, JoinPath(ptOverlayDir, "isolationOverlay_noSub_baseVariant_variantA.png"));
                        }
                        
                        for (auto* h : h3Vars) delete h;
                    }
                }
            }
            
            // ------ noSub vs variantA: <E_T^iso> vs centrality per pT bin ------
            {
                const string noSubVsVarADir = JoinPath(perVariantOverlayBase, "noSub_vsVariantA");
                EnsureDir(noSubVsVarADir);
                
                // Use UE-consistent colors (variantA = kOrange+7, not kRed+1)
                const int nvColors[4] = {kBlack, kBlue + 1, kOrange + 7, kGreen + 2};
                
                const double centLo2 = centBins.front().lo;
                const double centHi2 = centBins.back().hi;
                
                for (int ipt = 0; ipt < kNPtBins; ++ipt)
                {
                    const PtBin& b = PtBins()[ipt];
                    
                    struct VarEntry { std::size_t idx; vector<double> x, y, ey; };
                    vector<VarEntry> entries;
                    
                    double yMinNV = std::numeric_limits<double>::max();
                    double yMaxNV = -std::numeric_limits<double>::max();
                    
                    for (std::size_t iv : {std::size_t(0), std::size_t(2)})
                    {
                        if (iv >= handles.size()) continue;
                        auto& H = handles[iv];
                        if (!H.file) continue;
                        
                        TDirectory* aaTopNV = H.file->GetDirectory(trigAA.c_str());
                        if (!aaTopNV) continue;
                        
                        VarEntry E;
                        E.idx = iv;
                        
                        for (std::size_t ic = 0; ic < centBins.size(); ++ic)
                        {
                            const auto& cb = centBins[ic];
                            const string hName = "h_Eiso" + b.suffix + cb.suffix;
                            TH1* hSrc = dynamic_cast<TH1*>(aaTopNV->Get(hName.c_str()));
                            if (!hSrc) continue;
                            
                            const double mean = hSrc->GetMean();
                            const double err  = (hSrc->GetEntries() > 0.0) ? hSrc->GetMeanError() : 0.0;
                            
                            E.x.push_back(0.5 * (cb.lo + cb.hi));
                            E.y.push_back(mean);
                            E.ey.push_back(err);
                            
                            yMinNV = std::min(yMinNV, mean - err);
                            yMaxNV = std::max(yMaxNV, mean + err);
                        }
                        
                        if (!E.x.empty()) entries.push_back(std::move(E));
                    }
                    
                    if (entries.empty()) continue;
                    
                    double ppMean = 0.0, ppErr = 0.0;
                    bool havePP = false;
                    if (ppTop)
                    {
                        const string hPPName = "h_Eiso" + b.suffix;
                        TH1* hPPsrc = dynamic_cast<TH1*>(ppTop->Get(hPPName.c_str()));
                        if (hPPsrc && hPPsrc->GetEntries() > 0.0)
                        {
                            ppMean = hPPsrc->GetMean();
                            ppErr  = hPPsrc->GetMeanError();
                            yMinNV = std::min(yMinNV, ppMean - ppErr);
                            yMaxNV = std::max(yMaxNV, ppMean + ppErr);
                            havePP = true;
                        }
                    }
                    
                    if (!std::isfinite(yMinNV) || !std::isfinite(yMaxNV))
                    {
                        yMinNV = 0.0;
                        yMaxNV = 1.0;
                    }
                    const double padNV = (yMaxNV > yMinNV) ? (0.40 * (yMaxNV - yMinNV)) : 0.25;
                    
                    TCanvas cNV(
                                TString::Format("c_meanIsoEt_noSubVsVarA_%s_%s",
                                                trigAA.c_str(), b.folder.c_str()).Data(),
                                "c_meanIsoEt_noSubVsVarA", 900, 700
                                );
                    ApplyCanvasMargins1D(cNV);
                    cNV.cd();
                    
                    TH1F hFrameNV(
                                  TString::Format("hFrame_meanIsoEt_noSubVsVarA_%s_%s",
                                                  trigAA.c_str(), b.folder.c_str()).Data(),
                                  "", 100, centLo2, centHi2
                                  );
                    hFrameNV.SetDirectory(nullptr);
                    hFrameNV.SetStats(0);
                    hFrameNV.SetMinimum(yMinNV - padNV);
                    hFrameNV.SetMaximum(yMaxNV + padNV);
                    hFrameNV.GetXaxis()->SetTitle("Centrality [%]");
                    hFrameNV.GetYaxis()->SetTitle("<E_{T}^{iso}> [GeV]");
                    hFrameNV.GetXaxis()->SetTitleSize(0.055);
                    hFrameNV.GetYaxis()->SetTitleSize(0.055);
                    hFrameNV.GetXaxis()->SetLabelSize(0.045);
                    hFrameNV.GetYaxis()->SetLabelSize(0.045);
                    hFrameNV.GetYaxis()->SetTitleOffset(1.15);
                    hFrameNV.Draw();
                    
                    vector<TGraphErrors*> keepNV;
                    
                    for (const auto& E : entries)
                    {
                        vector<double> exZero(E.x.size(), 0.0);
                        TGraphErrors* g = new TGraphErrors(
                                                           (int)E.x.size(),
                                                           &E.x[0], &E.y[0],
                                                           &exZero[0], &E.ey[0]
                                                           );
                        g->SetLineWidth(2);
                        g->SetLineColor(nvColors[E.idx]);
                        g->SetMarkerColor(nvColors[E.idx]);
                        g->SetMarkerStyle(20);
                        g->SetMarkerSize(1.2);
                        g->Draw("PE1 SAME");
                        keepNV.push_back(g);
                    }
                    
                    TGraphErrors* gPP = nullptr;
                    if (havePP && !entries.empty())
                    {
                        const auto& xRef = entries.front().x;
                        vector<double> ppY(xRef.size(), ppMean);
                        vector<double> ppEY(xRef.size(), ppErr);
                        vector<double> exZeroPP(xRef.size(), 0.0);
                        gPP = new TGraphErrors(
                                               (int)xRef.size(),
                                               &xRef[0], &ppY[0],
                                               &exZeroPP[0], &ppEY[0]
                                               );
                        gPP->SetLineWidth(2);
                        gPP->SetLineColor(kRed+1);
                        gPP->SetMarkerColor(kRed+1);
                        gPP->SetMarkerStyle(24);
                        gPP->SetMarkerSize(1.2);
                        gPP->Draw("PE1 SAME");
                    }
                    
                    TLegend legNV(0.56, 0.72, 0.92, 0.88);
                    legNV.SetBorderSize(0);
                    legNV.SetFillStyle(0);
                    legNV.SetTextFont(42);
                    legNV.SetTextSize(0.032);
                    for (std::size_t ig = 0; ig < keepNV.size() && ig < entries.size(); ++ig)
                        legNV.AddEntry(keepNV[ig], handles[entries[ig].idx].label.c_str(), "ep");
                    if (gPP) legNV.AddEntry(gPP, "pp reference", "ep");
                    legNV.Draw();
                    
                    TLatex tTitleNV;
                    tTitleNV.SetNDC(true);
                    tTitleNV.SetTextFont(42);
                    tTitleNV.SetTextAlign(23);
                    tTitleNV.SetTextSize(0.045);
                    tTitleNV.DrawLatex(0.50, 0.98,
                                       TString::Format("<E_{T}^{iso}> vs Centrality, noSub vs variantA, p_{T}^{#gamma} %d-%d GeV", b.lo, b.hi).Data());
                    
                    const double coneRValNV = (kAA_IsoConeR == "isoR40") ? 0.4 : 0.3;
                    const int    isoEtMaxNV = (kAA_IsoMode == "fixedIso5GeV") ? 5 : -1;
                    
                    TLatex tNV;
                    tNV.SetNDC(true);
                    tNV.SetTextFont(42);
                    tNV.SetTextAlign(13);
                    tNV.SetTextSize(0.028);
                    tNV.DrawLatex(0.20, 0.89, "Trigger = Photon 10 GeV + MBD NS #geq 2, vtx < 150 cm");
                    tNV.DrawLatex(0.20, 0.85, TString::Format("#DeltaR_{cone} < %.1f", coneRValNV).Data());
                    if (isoEtMaxNV > 0)
                        tNV.DrawLatex(0.20, 0.81, TString::Format("E_{T}^{iso} < %d GeV", isoEtMaxNV).Data());
                    
                    SaveCanvas(cNV, JoinPath(noSubVsVarADir,
                                             TString::Format("meanIsoEt_vs_cent_noSub_vsVariantA_%s.png", b.folder.c_str()).Data()));
                    
                    for (auto* g : keepNV) delete g;
                    if (gPP) { delete gPP; gPP = nullptr; }
                }
            }
            
            // ------ all 4 UE variants + PP: <E_T^iso> vs centrality per pT bin ------
            {
                const string fourWayDir = JoinPath(perVariantOverlayBase, "isoEtMEAN_4wayCompare");
                EnsureDir(fourWayDir);
                
                const int fwColors[4] = {kBlack, kBlue + 1, kOrange + 7, kGreen + 2};
                
                const double centLo4 = centBins.front().lo;
                const double centHi4 = centBins.back().hi;
                
                for (int ipt = 0; ipt < kNPtBins; ++ipt)
                {
                    const PtBin& b = PtBins()[ipt];
                    
                    struct VarEntry4 { std::size_t idx; vector<double> x, y, ey; };
                    vector<VarEntry4> entries4;
                    
                    double yMin4 = std::numeric_limits<double>::max();
                    double yMax4 = -std::numeric_limits<double>::max();
                    
                    for (std::size_t iv = 0; iv < handles.size(); ++iv)
                    {
                        auto& H = handles[iv];
                        if (!H.file) continue;
                        
                        TDirectory* aaTop4 = H.file->GetDirectory(trigAA.c_str());
                        if (!aaTop4) continue;
                        
                        VarEntry4 E;
                        E.idx = iv;
                        
                        for (std::size_t ic = 0; ic < centBins.size(); ++ic)
                        {
                            const auto& cb = centBins[ic];
                            const string hName = "h_Eiso" + b.suffix + cb.suffix;
                            TH1* hSrc = dynamic_cast<TH1*>(aaTop4->Get(hName.c_str()));
                            if (!hSrc) continue;
                            
                            const double mean = hSrc->GetMean();
                            const double err  = (hSrc->GetEntries() > 0.0) ? hSrc->GetMeanError() : 0.0;
                            
                            E.x.push_back(0.5 * (cb.lo + cb.hi));
                            E.y.push_back(mean);
                            E.ey.push_back(err);
                            
                            yMin4 = std::min(yMin4, mean - err);
                            yMax4 = std::max(yMax4, mean + err);
                        }
                        
                        if (!E.x.empty()) entries4.push_back(std::move(E));
                    }
                    
                    if (entries4.empty()) continue;
                    
                    double ppMean4 = 0.0, ppErr4 = 0.0;
                    bool havePP4 = false;
                    if (ppTop)
                    {
                        const string hPPName = "h_Eiso" + b.suffix;
                        TH1* hPPsrc = dynamic_cast<TH1*>(ppTop->Get(hPPName.c_str()));
                        if (hPPsrc && hPPsrc->GetEntries() > 0.0)
                        {
                            ppMean4 = hPPsrc->GetMean();
                            ppErr4  = hPPsrc->GetMeanError();
                            yMin4 = std::min(yMin4, ppMean4 - ppErr4);
                            yMax4 = std::max(yMax4, ppMean4 + ppErr4);
                            havePP4 = true;
                        }
                    }
                    
                    if (!std::isfinite(yMin4) || !std::isfinite(yMax4))
                    {
                        yMin4 = 0.0;
                        yMax4 = 1.0;
                    }
                    const double pad4 = (yMax4 > yMin4) ? (0.15 * (yMax4 - yMin4)) : 0.25;
                    
                    TCanvas c4W(
                                TString::Format("c_meanIsoEt_4way_%s_%s",
                                                trigAA.c_str(), b.folder.c_str()).Data(),
                                "c_meanIsoEt_4way", 900, 700
                                );
                    ApplyCanvasMargins1D(c4W);
                    c4W.cd();
                    
                    TH1F hFrame4W(
                                  TString::Format("hFrame_meanIsoEt_4way_%s_%s",
                                                  trigAA.c_str(), b.folder.c_str()).Data(),
                                  "", 100, centLo4, centHi4
                                  );
                    hFrame4W.SetDirectory(nullptr);
                    hFrame4W.SetStats(0);
                    hFrame4W.SetMinimum(yMin4 - pad4);
                    hFrame4W.SetMaximum(yMax4 + pad4);
                    hFrame4W.GetXaxis()->SetTitle("Centrality [%]");
                    hFrame4W.GetYaxis()->SetTitle("<E_{T}^{iso}> [GeV]");
                    hFrame4W.GetXaxis()->SetTitleSize(0.055);
                    hFrame4W.GetYaxis()->SetTitleSize(0.055);
                    hFrame4W.GetXaxis()->SetLabelSize(0.045);
                    hFrame4W.GetYaxis()->SetLabelSize(0.045);
                    hFrame4W.GetYaxis()->SetTitleOffset(1.15);
                    hFrame4W.Draw();
                    
                    vector<TGraphErrors*> keep4W;
                    
                    for (const auto& E : entries4)
                    {
                        vector<double> exZero(E.x.size(), 0.0);
                        TGraphErrors* g = new TGraphErrors(
                                                           (int)E.x.size(),
                                                           &E.x[0], &E.y[0],
                                                           &exZero[0], &E.ey[0]
                                                           );
                        g->SetLineWidth(2);
                        g->SetLineColor((E.idx < 4) ? fwColors[E.idx] : kBlack);
                        g->SetMarkerColor((E.idx < 4) ? fwColors[E.idx] : kBlack);
                        g->SetMarkerStyle(20);
                        g->SetMarkerSize(1.2);
                        g->Draw("PE1 SAME");
                        keep4W.push_back(g);
                    }
                    
                    TGraphErrors* gPP4 = nullptr;
                    if (havePP4 && !entries4.empty())
                    {
                        const auto& xRef = entries4.front().x;
                        vector<double> ppY(xRef.size(), ppMean4);
                        vector<double> ppEY(xRef.size(), ppErr4);
                        vector<double> exZeroPP(xRef.size(), 0.0);
                        gPP4 = new TGraphErrors(
                                                (int)xRef.size(),
                                                &xRef[0], &ppY[0],
                                                &exZeroPP[0], &ppEY[0]
                                                );
                        gPP4->SetLineWidth(2);
                        gPP4->SetLineColor(kRed+1);
                        gPP4->SetMarkerColor(kRed+1);
                        gPP4->SetMarkerStyle(24);
                        gPP4->SetMarkerSize(1.2);
                        gPP4->Draw("PE1 SAME");
                    }
                    
                    TLegend leg4W(0.56, 0.68, 0.92, 0.88);
                    leg4W.SetBorderSize(0);
                    leg4W.SetFillStyle(0);
                    leg4W.SetTextFont(42);
                    leg4W.SetTextSize(0.032);
                    for (std::size_t ig = 0; ig < keep4W.size() && ig < entries4.size(); ++ig)
                        leg4W.AddEntry(keep4W[ig], handles[entries4[ig].idx].label.c_str(), "ep");
                    if (gPP4) leg4W.AddEntry(gPP4, "pp reference", "ep");
                    leg4W.Draw();
                    
                    TLatex tTitle4W;
                    tTitle4W.SetNDC(true);
                    tTitle4W.SetTextFont(42);
                    tTitle4W.SetTextAlign(23);
                    tTitle4W.SetTextSize(0.045);
                    tTitle4W.DrawLatex(0.50, 0.98,
                                       TString::Format("<E_{T}^{iso}> vs Centrality, all UE variants + pp, p_{T}^{#gamma} %d-%d GeV", b.lo, b.hi).Data());
                    
                    const double coneRVal4W = (kAA_IsoConeR == "isoR40") ? 0.4 : 0.3;
                    const int    isoEtMax4W = (kAA_IsoMode == "fixedIso5GeV") ? 5 : -1;
                    
                    TLatex t4W;
                    t4W.SetNDC(true);
                    t4W.SetTextFont(42);
                    t4W.SetTextAlign(13);
                    t4W.SetTextSize(0.028);
                    t4W.DrawLatex(0.20, 0.58, "Trigger = Photon 10 GeV + MBD NS #geq 2, vtx < 150 cm");
                    t4W.DrawLatex(0.20, 0.54, TString::Format("#DeltaR_{cone} < %.1f", coneRVal4W).Data());
                    if (isoEtMax4W > 0)
                        t4W.DrawLatex(0.20, 0.50, TString::Format("E_{T}^{iso} < %d GeV", isoEtMax4W).Data());
                    
                    SaveCanvas(c4W, JoinPath(fourWayDir,
                                             TString::Format("meanIsoEt_vs_cent_allVariants_pp_%s.png", b.folder.c_str()).Data()));
                    
                    for (auto* g : keep4W) delete g;
                    if (gPP4) { delete gPP4; gPP4 = nullptr; }
                }
            }
            
            // ------ noSub + variantA + variantB + PP: <E_T^iso> vs centrality per pT bin ------
            {
                const string threeWayDir = JoinPath(perVariantOverlayBase, "isoEtMEAN_3wayCompare");
                EnsureDir(threeWayDir);
                
                const int fwColors[4] = {kBlack, kBlue + 1, kOrange + 7, kGreen + 2};
                
                const double centLo3 = centBins.front().lo;
                const double centHi3 = centBins.back().hi;
                
                for (int ipt = 0; ipt < kNPtBins; ++ipt)
                {
                    const PtBin& b = PtBins()[ipt];
                    
                    struct VarEntry3 { std::size_t idx; vector<double> x, y, ey; };
                    vector<VarEntry3> entries3;
                    
                    double yMin3 = std::numeric_limits<double>::max();
                    double yMax3 = -std::numeric_limits<double>::max();
                    
                    // indices 0=noSub, 2=variantA, 3=variantB (skip 1=baseVariant)
                    for (std::size_t iv : {std::size_t(0), std::size_t(2), std::size_t(3)})
                    {
                        if (iv >= handles.size()) continue;
                        auto& H = handles[iv];
                        if (!H.file) continue;
                        
                        TDirectory* aaTop3 = H.file->GetDirectory(trigAA.c_str());
                        if (!aaTop3) continue;
                        
                        VarEntry3 E;
                        E.idx = iv;
                        
                        for (std::size_t ic = 0; ic < centBins.size(); ++ic)
                        {
                            const auto& cb = centBins[ic];
                            const string hName = "h_Eiso" + b.suffix + cb.suffix;
                            TH1* hSrc = dynamic_cast<TH1*>(aaTop3->Get(hName.c_str()));
                            if (!hSrc) continue;
                            
                            const double mean = hSrc->GetMean();
                            const double err  = (hSrc->GetEntries() > 0.0) ? hSrc->GetMeanError() : 0.0;
                            
                            E.x.push_back(0.5 * (cb.lo + cb.hi));
                            E.y.push_back(mean);
                            E.ey.push_back(err);
                            
                            yMin3 = std::min(yMin3, mean - err);
                            yMax3 = std::max(yMax3, mean + err);
                        }
                        
                        if (!E.x.empty()) entries3.push_back(std::move(E));
                    }
                    
                    if (entries3.empty()) continue;
                    
                    double ppMean3 = 0.0, ppErr3 = 0.0;
                    bool havePP3 = false;
                    if (ppTop)
                    {
                        const string hPPName = "h_Eiso" + b.suffix;
                        TH1* hPPsrc = dynamic_cast<TH1*>(ppTop->Get(hPPName.c_str()));
                        if (hPPsrc && hPPsrc->GetEntries() > 0.0)
                        {
                            ppMean3 = hPPsrc->GetMean();
                            ppErr3  = hPPsrc->GetMeanError();
                            yMin3 = std::min(yMin3, ppMean3 - ppErr3);
                            yMax3 = std::max(yMax3, ppMean3 + ppErr3);
                            havePP3 = true;
                        }
                    }
                    
                    if (!std::isfinite(yMin3) || !std::isfinite(yMax3))
                    {
                        yMin3 = 0.0;
                        yMax3 = 1.0;
                    }
                    const double pad3 = (yMax3 > yMin3) ? (0.45 * (yMax3 - yMin3)) : 0.25;
                    
                    TCanvas c3W(
                                TString::Format("c_meanIsoEt_3way_%s_%s",
                                                trigAA.c_str(), b.folder.c_str()).Data(),
                                "c_meanIsoEt_3way", 900, 700
                                );
                    ApplyCanvasMargins1D(c3W);
                    c3W.cd();
                    
                    TH1F hFrame3W(
                                  TString::Format("hFrame_meanIsoEt_3way_%s_%s",
                                                  trigAA.c_str(), b.folder.c_str()).Data(),
                                  "", 100, centLo3, centHi3
                                  );
                    hFrame3W.SetDirectory(nullptr);
                    hFrame3W.SetStats(0);
                    hFrame3W.SetMinimum(yMin3 - pad3);
                    hFrame3W.SetMaximum(yMax3 + pad3);
                    hFrame3W.GetXaxis()->SetTitle("Centrality [%]");
                    hFrame3W.GetYaxis()->SetTitle("<E_{T}^{iso}> [GeV]");
                    hFrame3W.GetXaxis()->SetTitleSize(0.055);
                    hFrame3W.GetYaxis()->SetTitleSize(0.055);
                    hFrame3W.GetXaxis()->SetLabelSize(0.045);
                    hFrame3W.GetYaxis()->SetLabelSize(0.045);
                    hFrame3W.GetYaxis()->SetTitleOffset(1.15);
                    hFrame3W.Draw();
                    
                    vector<TGraphErrors*> keep3W;
                    
                    for (const auto& E : entries3)
                    {
                        vector<double> exZero(E.x.size(), 0.0);
                        TGraphErrors* g = new TGraphErrors(
                                                           (int)E.x.size(),
                                                           &E.x[0], &E.y[0],
                                                           &exZero[0], &E.ey[0]
                                                           );
                        g->SetLineWidth(2);
                        g->SetLineColor((E.idx < 4) ? fwColors[E.idx] : kBlack);
                        g->SetMarkerColor((E.idx < 4) ? fwColors[E.idx] : kBlack);
                        g->SetMarkerStyle(20);
                        g->SetMarkerSize(1.2);
                        g->Draw("PE1 SAME");
                        keep3W.push_back(g);
                    }
                    
                    TGraphErrors* gPP3 = nullptr;
                    if (havePP3 && !entries3.empty())
                    {
                        const auto& xRef = entries3.front().x;
                        vector<double> ppY(xRef.size(), ppMean3);
                        vector<double> ppEY(xRef.size(), ppErr3);
                        vector<double> exZeroPP(xRef.size(), 0.0);
                        gPP3 = new TGraphErrors(
                                                (int)xRef.size(),
                                                &xRef[0], &ppY[0],
                                                &exZeroPP[0], &ppEY[0]
                                                );
                        gPP3->SetLineWidth(2);
                        gPP3->SetLineColor(kRed+1);
                        gPP3->SetMarkerColor(kRed+1);
                        gPP3->SetMarkerStyle(24);
                        gPP3->SetMarkerSize(1.2);
                        gPP3->Draw("PE1 SAME");
                    }
                    
                    TLegend leg3W(0.56, 0.68, 0.92, 0.88);
                    leg3W.SetBorderSize(0);
                    leg3W.SetFillStyle(0);
                    leg3W.SetTextFont(42);
                    leg3W.SetTextSize(0.032);
                    for (std::size_t ig = 0; ig < keep3W.size() && ig < entries3.size(); ++ig)
                        leg3W.AddEntry(keep3W[ig], handles[entries3[ig].idx].label.c_str(), "ep");
                    if (gPP3) leg3W.AddEntry(gPP3, "pp reference", "ep");
                    leg3W.Draw();
                    
                    TLatex tTitle3W;
                    tTitle3W.SetNDC(true);
                    tTitle3W.SetTextFont(42);
                    tTitle3W.SetTextAlign(23);
                    tTitle3W.SetTextSize(0.045);
                    tTitle3W.DrawLatex(0.50, 0.98,
                                       TString::Format("<E_{T}^{iso}> vs Centrality, noSub + varA + varB + pp, p_{T}^{#gamma} %d-%d GeV", b.lo, b.hi).Data());
                    
                    const double coneRVal3W = (kAA_IsoConeR == "isoR40") ? 0.4 : 0.3;
                    const int    isoEtMax3W = (kAA_IsoMode == "fixedIso5GeV") ? 5 : -1;
                    
                    TLatex t3W;
                    t3W.SetNDC(true);
                    t3W.SetTextFont(42);
                    t3W.SetTextAlign(13);
                    t3W.SetTextSize(0.028);
                    t3W.DrawLatex(0.20, 0.89, "Trigger = Photon 10 GeV + MBD NS #geq 2, vtx < 150 cm");
                    t3W.DrawLatex(0.20, 0.85, TString::Format("#DeltaR_{cone} < %.1f", coneRVal3W).Data());
                    if (isoEtMax3W > 0)
                        t3W.DrawLatex(0.20, 0.81, TString::Format("E_{T}^{iso} < %d GeV", isoEtMax3W).Data());
                    
                    SaveCanvas(c3W, JoinPath(threeWayDir,
                                             TString::Format("meanIsoEt_vs_cent_3variants_pp_%s.png", b.folder.c_str()).Data()));
                    
                    for (auto* g : keep3W) delete g;
                    if (gPP3) { delete gPP3; gPP3 = nullptr; }
                }
            }
            
        } // end if (!generateISOpTcentOverlaysONLY) — perVariantOverlay + multi-variant sections
        
        for (int ipt = 0; ipt < kNPtBins; ++ipt)
        {
            const PtBin& b = PtBins()[ipt];
            
            vector<TGraphErrors*> graphs;
            vector<std::size_t> graphIndices;
            double yMin = std::numeric_limits<double>::max();
            double yMax = -std::numeric_limits<double>::max();
            
            const double centLo = centBins.front().lo;
            const double centHi = centBins.back().hi;
            
            for (std::size_t iv = 0; iv < handles.size(); ++iv)
            {
                auto& H = handles[iv];
                if (!H.file) continue;
                
                TDirectory* aaTop = H.file->GetDirectory(trigAA.c_str());
                if (!aaTop) continue;
                
                vector<double> xCent;
                vector<double> exCent;
                vector<double> yCent;
                vector<double> eyCent;
                
                for (std::size_t ic = 0; ic < centBins.size(); ++ic)
                {
                    const auto& cb = centBins[ic];
                    const string hAAName = "h_Eiso" + b.suffix + cb.suffix;
                    TH1* hAAsrc = dynamic_cast<TH1*>(aaTop->Get(hAAName.c_str()));
                    if (!hAAsrc) continue;
                    
                    xCent.push_back(0.5 * (cb.lo + cb.hi));
                    exCent.push_back(0.5 * (cb.hi - cb.lo));
                    yCent.push_back(hAAsrc->GetMean());
                    eyCent.push_back((hAAsrc->GetEntries() > 0.0) ? hAAsrc->GetMeanError() : 0.0);
                    yMin = std::min(yMin, yCent.back() - eyCent.back());
                    yMax = std::max(yMax, yCent.back() + eyCent.back());
                }
                
                if (xCent.empty()) continue;
                
                TGraphErrors* g = new TGraphErrors((int)xCent.size(), &xCent[0], &yCent[0], &exCent[0], &eyCent[0]);
                g->SetLineWidth(2);
                g->SetLineColor((iv < 4) ? variantColors[iv] : kBlack);
                g->SetMarkerColor((iv < 4) ? variantColors[iv] : kBlack);
                g->SetMarkerStyle((iv < 4) ? variantMarkers[iv] : 20);
                g->SetMarkerSize(1.2);
                graphs.push_back(g);
                graphIndices.push_back(iv);
            }
            
            if (graphs.empty()) continue;
            
            if (!std::isfinite(yMin) || !std::isfinite(yMax))
            {
                yMin = 0.0;
                yMax = 1.0;
            }
            const double pad = (yMax > yMin) ? (0.15 * (yMax - yMin)) : 0.25;
            
            TCanvas cSummary(
                             TString::Format("c_meanIsoEtAllVariantsVsCent_%s_%s",
                                             trigAA.c_str(), b.folder.c_str()).Data(),
                             "c_meanIsoEtAllVariantsVsCent", 900, 700
                             );
            ApplyCanvasMargins1D(cSummary);
            cSummary.cd();
            
            TH1F hFrame(
                        TString::Format("hFrame_meanIsoEtAllVariantsVsCent_%s_%s",
                                        trigAA.c_str(), b.folder.c_str()).Data(),
                        "", 100, centLo, centHi
                        );
            hFrame.SetDirectory(nullptr);
            hFrame.SetStats(0);
            hFrame.SetMinimum(yMin - pad);
            hFrame.SetMaximum(yMax + pad);
            hFrame.GetXaxis()->SetTitle("Centrality [%]");
            hFrame.GetYaxis()->SetTitle("<E_{T}^{iso}> [GeV]");
            hFrame.GetXaxis()->SetTitleSize(0.055);
            hFrame.GetYaxis()->SetTitleSize(0.055);
            hFrame.GetXaxis()->SetLabelSize(0.045);
            hFrame.GetYaxis()->SetLabelSize(0.045);
            hFrame.GetYaxis()->SetTitleOffset(1.15);
            hFrame.Draw();
            
            for (auto* g : graphs) g->Draw("PE1 SAME");
            
            TLegend leg(0.56, 0.62, 0.92, 0.88);
            leg.SetBorderSize(0);
            leg.SetFillStyle(0);
            leg.SetTextFont(42);
            leg.SetTextSize(0.032);
            for (std::size_t ig = 0; ig < graphs.size(); ++ig)
            {
                leg.AddEntry(graphs[ig], handles[graphIndices[ig]].label.c_str(), "ep");
            }
            leg.Draw();
            
            TLatex tTitle;
            tTitle.SetNDC(true);
            tTitle.SetTextFont(42);
            tTitle.SetTextAlign(23);
            tTitle.SetTextSize(0.045);
            tTitle.DrawLatex(0.50, 0.98,
                             TString::Format("<E_{T}^{iso}> vs Centrality, all UE variants, p_{T}^{#gamma} %d-%d GeV", b.lo, b.hi).Data());
            
            const double coneRVal = (kAA_IsoConeR == "isoR40") ? 0.4 : 0.3;
            const int    isoEtMax = (kAA_IsoMode == "fixedIso5GeV") ? 5 : -1;
            
            TLatex t;
            t.SetNDC(true);
            t.SetTextFont(42);
            t.SetTextAlign(13);
            t.SetTextSize(0.028);
            t.DrawLatex(0.20, 0.58, TString::Format("Trigger = %s", trigDisplayLabel.c_str()).Data());
            t.DrawLatex(0.20, 0.54, TString::Format("#DeltaR_{cone} < %.1f", coneRVal).Data());
            if (isoEtMax > 0)
                t.DrawLatex(0.20, 0.50, TString::Format("E_{T}^{iso} < %d GeV", isoEtMax).Data());
            
            SaveCanvas(cSummary, JoinPath(meanIsoSummaryDir,
                                          TString::Format("meanIsoEt_allVariants_vs_cent_%s.png", b.folder.c_str()).Data()));
            
            for (auto* g : graphs) delete g;
        }
        
        // ====== centralitySummaryPerPt/summaryOutput: Gaussian mean & sigma vs centrality ======
        {
            const string gaussMeanDir  = JoinPath(meanIsoSummaryDir, "meanSummary");
            const string gaussSigmaDir = JoinPath(meanIsoSummaryDir, "sigmaSummary");
            EnsureDir(gaussMeanDir);
            EnsureDir(gaussSigmaDir);
            
            for (int ipt = 0; ipt < kNPtBins; ++ipt)
            {
                const PtBin& b = PtBins()[ipt];
                
                // --- Gaussian Mean vs Centrality ---
                {
                    vector<TGraphErrors*> gGraphs;
                    vector<std::size_t>   gIndices;
                    double yMinG = 1e30, yMaxG = -1e30;
                    
                    for (std::size_t iv = 0; iv < handles.size(); ++iv)
                    {
                        if (!handles[iv].file) continue;
                        vector<double> xC, exC, yC, eyC;
                        for (std::size_t ic = 0; ic < centBins.size(); ++ic)
                        {
                            if (!gaussFilled[iv][ipt][ic]) continue;
                            xC.push_back(0.5 * (centBins[ic].lo + centBins[ic].hi));
                            exC.push_back(0.5 * (centBins[ic].hi - centBins[ic].lo));
                            yC.push_back(gaussMean[iv][ipt][ic]);
                            eyC.push_back(gaussMeanErr[iv][ipt][ic]);
                            yMinG = std::min(yMinG, yC.back() - eyC.back());
                            yMaxG = std::max(yMaxG, yC.back() + eyC.back());
                        }
                        if (xC.empty()) continue;
                        TGraphErrors* g = new TGraphErrors((int)xC.size(), &xC[0], &yC[0], &exC[0], &eyC[0]);
                        g->SetLineWidth(2);
                        g->SetLineColor((iv < 4) ? variantColors[iv] : kBlack);
                        g->SetMarkerColor((iv < 4) ? variantColors[iv] : kBlack);
                        g->SetMarkerStyle((iv < 4) ? variantMarkers[iv] : 20);
                        g->SetMarkerSize(1.2);
                        gGraphs.push_back(g);
                        gIndices.push_back(iv);
                    }
                    if (!gGraphs.empty())
                    {
                        if (!std::isfinite(yMinG) || !std::isfinite(yMaxG)) { yMinG = 0; yMaxG = 1; }
                        const double padG = (yMaxG > yMinG) ? 0.20 * (yMaxG - yMinG) : 0.5;
                        TCanvas cGM(TString::Format("c_gaussMeanVsCent_%s", b.folder.c_str()).Data(), "", 900, 700);
                        ApplyCanvasMargins1D(cGM); cGM.cd();
                        TH1F hFrGM(TString::Format("hFr_gaussMeanVsCent_%s", b.folder.c_str()).Data(),
                                   "", 100, centBins.front().lo, centBins.back().hi);
                        hFrGM.SetDirectory(nullptr); hFrGM.SetStats(0);
                        hFrGM.SetMinimum(yMinG - padG); hFrGM.SetMaximum(yMaxG + padG);
                        hFrGM.GetXaxis()->SetTitle("Centrality [%]");
                        hFrGM.GetYaxis()->SetTitle("Gaussian #mu [GeV]");
                        hFrGM.GetXaxis()->SetTitleSize(0.055); hFrGM.GetYaxis()->SetTitleSize(0.055);
                        hFrGM.GetXaxis()->SetLabelSize(0.045); hFrGM.GetYaxis()->SetLabelSize(0.045);
                        hFrGM.GetYaxis()->SetTitleOffset(1.15);
                        hFrGM.Draw();
                        for (auto* g : gGraphs) g->Draw("PE1 SAME");
                        TLegend lgGM(0.56, 0.62, 0.92, 0.88);
                        lgGM.SetBorderSize(0); lgGM.SetFillStyle(0); lgGM.SetTextFont(42); lgGM.SetTextSize(0.032);
                        for (std::size_t ig = 0; ig < gGraphs.size(); ++ig)
                            lgGM.AddEntry(gGraphs[ig], handles[gIndices[ig]].label.c_str(), "ep");
                        lgGM.Draw();
                        TLatex tGM; tGM.SetNDC(true); tGM.SetTextFont(42); tGM.SetTextAlign(23); tGM.SetTextSize(0.042);
                        tGM.DrawLatex(0.50, 0.98,
                            TString::Format("Gaussian #mu vs Centrality, p_{T}^{#gamma} %d-%d GeV", b.lo, b.hi).Data());
                        TLatex tGMi; tGMi.SetNDC(true); tGMi.SetTextFont(42); tGMi.SetTextAlign(13); tGMi.SetTextSize(0.028);
                        tGMi.DrawLatex(0.20, 0.58, TString::Format("Trigger = %s", trigDisplayLabel.c_str()).Data());
                        tGMi.DrawLatex(0.20, 0.54, TString::Format("#DeltaR_{cone} < %.1f", (kAA_IsoConeR == "isoR40") ? 0.4 : 0.3).Data());
                        SaveCanvas(cGM, JoinPath(gaussMeanDir,
                            TString::Format("gaussMean_allVariants_vs_cent_%s.png", b.folder.c_str()).Data()));
                        for (auto* g : gGraphs) delete g;
                    }
                }
                
                // --- Gaussian Sigma vs Centrality ---
                {
                    vector<TGraphErrors*> sGraphs;
                    vector<std::size_t>   sIndices;
                    double yMinS = 1e30, yMaxS = -1e30;
                    
                    for (std::size_t iv = 0; iv < handles.size(); ++iv)
                    {
                        if (!handles[iv].file) continue;
                        vector<double> xC, exC, yC, eyC;
                        for (std::size_t ic = 0; ic < centBins.size(); ++ic)
                        {
                            if (!gaussFilled[iv][ipt][ic]) continue;
                            xC.push_back(0.5 * (centBins[ic].lo + centBins[ic].hi));
                            exC.push_back(0.5 * (centBins[ic].hi - centBins[ic].lo));
                            yC.push_back(gaussSigma[iv][ipt][ic]);
                            eyC.push_back(gaussSigmaErr[iv][ipt][ic]);
                            yMinS = std::min(yMinS, yC.back() - eyC.back());
                            yMaxS = std::max(yMaxS, yC.back() + eyC.back());
                        }
                        if (xC.empty()) continue;
                        TGraphErrors* g = new TGraphErrors((int)xC.size(), &xC[0], &yC[0], &exC[0], &eyC[0]);
                        g->SetLineWidth(2);
                        g->SetLineColor((iv < 4) ? variantColors[iv] : kBlack);
                        g->SetMarkerColor((iv < 4) ? variantColors[iv] : kBlack);
                        g->SetMarkerStyle((iv < 4) ? variantMarkers[iv] : 20);
                        g->SetMarkerSize(1.2);
                        sGraphs.push_back(g);
                        sIndices.push_back(iv);
                    }
                    if (!sGraphs.empty())
                    {
                        if (!std::isfinite(yMinS) || !std::isfinite(yMaxS)) { yMinS = 0; yMaxS = 5; }
                        const double padS = (yMaxS > yMinS) ? 0.20 * (yMaxS - yMinS) : 0.5;
                        TCanvas cGS(TString::Format("c_gaussSigmaVsCent_%s", b.folder.c_str()).Data(), "", 900, 700);
                        ApplyCanvasMargins1D(cGS); cGS.cd();
                        TH1F hFrGS(TString::Format("hFr_gaussSigmaVsCent_%s", b.folder.c_str()).Data(),
                                   "", 100, centBins.front().lo, centBins.back().hi);
                        hFrGS.SetDirectory(nullptr); hFrGS.SetStats(0);
                        hFrGS.SetMinimum(std::max(0.0, yMinS - padS)); hFrGS.SetMaximum(yMaxS + padS);
                        hFrGS.GetXaxis()->SetTitle("Centrality [%]");
                        hFrGS.GetYaxis()->SetTitle("Gaussian #sigma [GeV]");
                        hFrGS.GetXaxis()->SetTitleSize(0.055); hFrGS.GetYaxis()->SetTitleSize(0.055);
                        hFrGS.GetXaxis()->SetLabelSize(0.045); hFrGS.GetYaxis()->SetLabelSize(0.045);
                        hFrGS.GetYaxis()->SetTitleOffset(1.15);
                        hFrGS.Draw();
                        for (auto* g : sGraphs) g->Draw("PE1 SAME");
                        TLegend lgGS(0.56, 0.62, 0.92, 0.88);
                        lgGS.SetBorderSize(0); lgGS.SetFillStyle(0); lgGS.SetTextFont(42); lgGS.SetTextSize(0.032);
                        for (std::size_t ig = 0; ig < sGraphs.size(); ++ig)
                            lgGS.AddEntry(sGraphs[ig], handles[sIndices[ig]].label.c_str(), "ep");
                        lgGS.Draw();
                        TLatex tGS; tGS.SetNDC(true); tGS.SetTextFont(42); tGS.SetTextAlign(23); tGS.SetTextSize(0.042);
                        tGS.DrawLatex(0.50, 0.98,
                            TString::Format("Gaussian #sigma vs Centrality, p_{T}^{#gamma} %d-%d GeV", b.lo, b.hi).Data());
                        TLatex tGSi; tGSi.SetNDC(true); tGSi.SetTextFont(42); tGSi.SetTextAlign(13); tGSi.SetTextSize(0.028);
                        tGSi.DrawLatex(0.20, 0.58, TString::Format("Trigger = %s", trigDisplayLabel.c_str()).Data());
                        tGSi.DrawLatex(0.20, 0.54, TString::Format("#DeltaR_{cone} < %.1f", (kAA_IsoConeR == "isoR40") ? 0.4 : 0.3).Data());
                        SaveCanvas(cGS, JoinPath(gaussSigmaDir,
                            TString::Format("gaussSigma_allVariants_vs_cent_%s.png", b.folder.c_str()).Data()));
                        for (auto* g : sGraphs) delete g;
                    }
                }
            }
        }
        
        // ====== pTsummaryPerCentrality: <E_T^iso> vs pT for each centrality bin, all UE variants ======
        {
            
            for (std::size_t ic = 0; ic < centBins.size(); ++ic)
            {
                const auto& cb = centBins[ic];
                
                vector<TGraphErrors*> graphs;
                vector<std::size_t> graphIndices;
                double yMin = std::numeric_limits<double>::max();
                double yMax = -std::numeric_limits<double>::max();
                
                for (std::size_t iv = 0; iv < handles.size(); ++iv)
                {
                    if (iv == 1) continue;  // skip baseVariant
                    auto& H = handles[iv];
                    if (!H.file) continue;
                    
                    TDirectory* aaTop = H.file->GetDirectory(trigAA.c_str());
                    if (!aaTop) continue;
                    
                    vector<double> xPt;
                    vector<double> exPt;
                    vector<double> yPt;
                    vector<double> eyPt;
                    
                    for (int ipt = 1; ipt < kNPtBins; ++ipt)
                    {
                        const PtBin& b = PtBins()[ipt];
                        const string hAAName = "h_Eiso" + b.suffix + cb.suffix;
                        TH1* hAAsrc = dynamic_cast<TH1*>(aaTop->Get(hAAName.c_str()));
                        if (!hAAsrc) continue;
                        
                        xPt.push_back(0.5 * (kPtEdges[(std::size_t)ipt] + kPtEdges[(std::size_t)ipt + 1]));
                        exPt.push_back(0.0);
                        yPt.push_back(hAAsrc->GetMean());
                        eyPt.push_back((hAAsrc->GetEntries() > 0.0) ? hAAsrc->GetMeanError() : 0.0);
                        
                        yMin = std::min(yMin, yPt.back() - eyPt.back());
                        yMax = std::max(yMax, yPt.back() + eyPt.back());
                    }
                    
                    if (xPt.empty()) continue;
                    
                    TGraphErrors* g = new TGraphErrors((int)xPt.size(), &xPt[0], &yPt[0], &exPt[0], &eyPt[0]);
                    g->SetLineWidth(2);
                    g->SetLineColor((iv < 4) ? variantColors[iv] : kBlack);
                    g->SetMarkerColor((iv < 4) ? variantColors[iv] : kBlack);
                    g->SetMarkerStyle((iv < 4) ? variantMarkers[iv] : 20);
                    g->SetMarkerSize(1.2);
                    graphs.push_back(g);
                    graphIndices.push_back(iv);
                }
                
                if (graphs.empty()) continue;
                
                // -- collect PP reference vs pT --
                vector<double> xPtPP, exPtPP, yPtPP, eyPtPP;
                if (ppTop)
                {
                    for (int ipt = 0; ipt < kNPtBins; ++ipt)
                    {
                        const PtBin& bp = PtBins()[ipt];
                        const string hPPName = "h_Eiso" + bp.suffix;
                        TH1* hPPsrc = dynamic_cast<TH1*>(ppTop->Get(hPPName.c_str()));
                        if (!hPPsrc || hPPsrc->GetEntries() <= 0.0) continue;
                        
                        xPtPP.push_back(0.5 * (kPtEdges[(std::size_t)ipt] + kPtEdges[(std::size_t)ipt + 1]));
                        exPtPP.push_back(0.0);
                        yPtPP.push_back(hPPsrc->GetMean());
                        eyPtPP.push_back(hPPsrc->GetMeanError());
                        
                        yMin = std::min(yMin, yPtPP.back() - eyPtPP.back());
                        yMax = std::max(yMax, yPtPP.back() + eyPtPP.back());
                    }
                }
                
                if (!std::isfinite(yMin) || !std::isfinite(yMax))
                {
                    yMin = 0.0;
                    yMax = 1.0;
                }
                const double pad = (yMax > yMin) ? (0.55 * (yMax - yMin)) : 0.25;
                
                TCanvas cPtSummary(
                                   TString::Format("c_meanIsoEtAllVariantsVsPt_%s_%s",
                                                   trigAA.c_str(), cb.folder.c_str()).Data(),
                                   "c_meanIsoEtAllVariantsVsPt", 900, 700
                                   );
                ApplyCanvasMargins1D(cPtSummary);
                cPtSummary.cd();
                
                TH1F hFrame(
                            TString::Format("hFrame_meanIsoEtAllVariantsVsPt_%s_%s",
                                            trigAA.c_str(), cb.folder.c_str()).Data(),
                            "", 100, kPtEdges[1], kPtEdges.back()
                            );
                hFrame.SetDirectory(nullptr);
                hFrame.SetStats(0);
                hFrame.SetMinimum(yMin - pad);
                hFrame.SetMaximum(yMax + pad);
                hFrame.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                hFrame.GetYaxis()->SetTitle("<E_{T}^{iso}> [GeV]");
                hFrame.GetXaxis()->SetTitleSize(0.055);
                hFrame.GetYaxis()->SetTitleSize(0.055);
                hFrame.GetXaxis()->SetLabelSize(0.045);
                hFrame.GetYaxis()->SetLabelSize(0.045);
                hFrame.GetYaxis()->SetTitleOffset(1.15);
                hFrame.Draw();
                
                for (auto* g : graphs) g->Draw("PE1 SAME");
                
                TGraphErrors* gPPpt = nullptr;
                if (!xPtPP.empty())
                {
                    gPPpt = new TGraphErrors((int)xPtPP.size(), &xPtPP[0], &yPtPP[0], &exPtPP[0], &eyPtPP[0]);
                    gPPpt->SetLineWidth(2);
                    gPPpt->SetLineColor(kRed + 1);
                    gPPpt->SetMarkerColor(kRed + 1);
                    gPPpt->SetMarkerStyle(24);
                    gPPpt->SetMarkerSize(1.2);
                    gPPpt->Draw("PE1 SAME");
                }
                
                TLegend leg(0.20, 0.75, 0.60, 0.88);
                leg.SetBorderSize(0);
                leg.SetFillStyle(0);
                leg.SetTextFont(42);
                leg.SetTextSize(0.026);
                leg.SetNColumns(2);
                for (std::size_t ig = 0; ig < graphs.size(); ++ig)
                {
                    leg.AddEntry(graphs[ig], handles[graphIndices[ig]].label.c_str(), "ep");
                }
                if (gPPpt) leg.AddEntry(gPPpt, "pp reference", "ep");
                leg.Draw();
                
                TLatex tTitle;
                tTitle.SetNDC(true);
                tTitle.SetTextFont(42);
                tTitle.SetTextAlign(23);
                tTitle.SetTextSize(0.045);
                tTitle.DrawLatex(0.50, 0.98,
                                 TString::Format("<E_{T}^{iso}> vs p_{T}^{#gamma}, all UE variants, %d-%d%% Cent AuAu", cb.lo, cb.hi).Data());
                
                const double coneRVal = (kAA_IsoConeR == "isoR40") ? 0.4 : 0.3;
                const int    isoEtMax = (kAA_IsoMode == "fixedIso5GeV") ? 5 : -1;
                
                TLatex t;
                t.SetNDC(true);
                t.SetTextFont(42);
                t.SetTextAlign(33);
                t.SetTextSize(0.028);
                t.DrawLatex(0.92, 0.89, "Trigger = Photon 10 GeV + MBD NS #geq 2, vtx < 150 cm");
                t.DrawLatex(0.92, 0.85, TString::Format("#DeltaR_{cone} < %.1f", coneRVal).Data());
                if (isoEtMax > 0)
                    t.DrawLatex(0.92, 0.81, TString::Format("E_{T}^{iso} < %d GeV", isoEtMax).Data());
                
                SaveCanvas(cPtSummary, JoinPath(ptMeanIsoSummaryDir,
                                                TString::Format("meanIsoEt_allVariants_vs_pT_%s.png", cb.folder.c_str()).Data()));
                
                for (auto* g : graphs) delete g;
                if (gPPpt) delete gPPpt;
            }
        }
        
        // ====== pTsummaryPerCentrality/summaryOutput: Gaussian mean & sigma vs pT ======
        {
            const string gaussMeanPtDir  = JoinPath(ptMeanIsoSummaryDir, "meanSummary");
            const string gaussSigmaPtDir = JoinPath(ptMeanIsoSummaryDir, "sigmaSummary");
            EnsureDir(gaussMeanPtDir);
            EnsureDir(gaussSigmaPtDir);
            
            for (std::size_t ic = 0; ic < centBins.size(); ++ic)
            {
                const auto& cb = centBins[ic];
                
                // --- Gaussian Mean vs pT ---
                {
                    vector<TGraphErrors*> gGraphs;
                    vector<std::size_t>   gIdx;
                    double yMinG = 1e30, yMaxG = -1e30;
                    for (std::size_t iv = 0; iv < handles.size(); ++iv)
                    {
                        if (!handles[iv].file) continue;
                        vector<double> xP, exP, yP, eyP;
                        for (int ipt = 0; ipt < kNPtBins; ++ipt)
                        {
                            if (!gaussFilled[iv][ipt][ic]) continue;
                            xP.push_back(0.5 * (kPtEdges[(std::size_t)ipt] + kPtEdges[(std::size_t)ipt + 1]));
                            exP.push_back(0.5 * (kPtEdges[(std::size_t)ipt + 1] - kPtEdges[(std::size_t)ipt]));
                            yP.push_back(gaussMean[iv][ipt][ic]);
                            eyP.push_back(gaussMeanErr[iv][ipt][ic]);
                            yMinG = std::min(yMinG, yP.back() - eyP.back());
                            yMaxG = std::max(yMaxG, yP.back() + eyP.back());
                        }
                        if (xP.empty()) continue;
                        TGraphErrors* g = new TGraphErrors((int)xP.size(), &xP[0], &yP[0], &exP[0], &eyP[0]);
                        g->SetLineWidth(2);
                        g->SetLineColor((iv < 4) ? variantColors[iv] : kBlack);
                        g->SetMarkerColor((iv < 4) ? variantColors[iv] : kBlack);
                        g->SetMarkerStyle((iv < 4) ? variantMarkers[iv] : 20);
                        g->SetMarkerSize(1.2);
                        gGraphs.push_back(g);
                        gIdx.push_back(iv);
                    }
                    if (!gGraphs.empty())
                    {
                        if (!std::isfinite(yMinG) || !std::isfinite(yMaxG)) { yMinG = 0; yMaxG = 1; }
                        const double padG = (yMaxG > yMinG) ? 0.20 * (yMaxG - yMinG) : 0.5;
                        TCanvas cGM(TString::Format("c_gaussMeanVsPt_%s", cb.folder.c_str()).Data(), "", 900, 700);
                        ApplyCanvasMargins1D(cGM); cGM.cd();
                        TH1F hFr(TString::Format("hFr_gaussMeanVsPt_%s", cb.folder.c_str()).Data(),
                                 "", 100, kPtEdges.front(), kPtEdges.back());
                        hFr.SetDirectory(nullptr); hFr.SetStats(0);
                        hFr.SetMinimum(yMinG - padG); hFr.SetMaximum(yMaxG + padG);
                        hFr.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                        hFr.GetYaxis()->SetTitle("Gaussian #mu [GeV]");
                        hFr.GetXaxis()->SetTitleSize(0.055); hFr.GetYaxis()->SetTitleSize(0.055);
                        hFr.GetXaxis()->SetLabelSize(0.045); hFr.GetYaxis()->SetLabelSize(0.045);
                        hFr.GetYaxis()->SetTitleOffset(1.15);
                        hFr.Draw();
                        for (auto* g : gGraphs) g->Draw("PE1 SAME");
                        TLegend lgGM(0.56, 0.62, 0.92, 0.88);
                        lgGM.SetBorderSize(0); lgGM.SetFillStyle(0); lgGM.SetTextFont(42); lgGM.SetTextSize(0.032);
                        for (std::size_t ig = 0; ig < gGraphs.size(); ++ig)
                            lgGM.AddEntry(gGraphs[ig], handles[gIdx[ig]].label.c_str(), "ep");
                        lgGM.Draw();
                        TLatex tGM; tGM.SetNDC(true); tGM.SetTextFont(42); tGM.SetTextAlign(23); tGM.SetTextSize(0.042);
                        tGM.DrawLatex(0.50, 0.98,
                            TString::Format("Gaussian #mu vs p_{T}^{#gamma}, %d-%d%% Cent", cb.lo, cb.hi).Data());
                        TLatex tGMi; tGMi.SetNDC(true); tGMi.SetTextFont(42); tGMi.SetTextAlign(13); tGMi.SetTextSize(0.028);
                        tGMi.DrawLatex(0.20, 0.58, TString::Format("Trigger = %s", trigDisplayLabel.c_str()).Data());
                        tGMi.DrawLatex(0.20, 0.54, TString::Format("#DeltaR_{cone} < %.1f", (kAA_IsoConeR == "isoR40") ? 0.4 : 0.3).Data());
                        SaveCanvas(cGM, JoinPath(gaussMeanPtDir,
                            TString::Format("gaussMean_allVariants_vs_pT_%s.png", cb.folder.c_str()).Data()));
                        for (auto* g : gGraphs) delete g;
                    }
                }
                
                // --- Gaussian Sigma vs pT ---
                {
                    vector<TGraphErrors*> sGraphs;
                    vector<std::size_t>   sIdx;
                    double yMinS = 1e30, yMaxS = -1e30;
                    for (std::size_t iv = 0; iv < handles.size(); ++iv)
                    {
                        if (!handles[iv].file) continue;
                        vector<double> xP, exP, yP, eyP;
                        for (int ipt = 0; ipt < kNPtBins; ++ipt)
                        {
                            if (!gaussFilled[iv][ipt][ic]) continue;
                            xP.push_back(0.5 * (kPtEdges[(std::size_t)ipt] + kPtEdges[(std::size_t)ipt + 1]));
                            exP.push_back(0.5 * (kPtEdges[(std::size_t)ipt + 1] - kPtEdges[(std::size_t)ipt]));
                            yP.push_back(gaussSigma[iv][ipt][ic]);
                            eyP.push_back(gaussSigmaErr[iv][ipt][ic]);
                            yMinS = std::min(yMinS, yP.back() - eyP.back());
                            yMaxS = std::max(yMaxS, yP.back() + eyP.back());
                        }
                        if (xP.empty()) continue;
                        TGraphErrors* g = new TGraphErrors((int)xP.size(), &xP[0], &yP[0], &exP[0], &eyP[0]);
                        g->SetLineWidth(2);
                        g->SetLineColor((iv < 4) ? variantColors[iv] : kBlack);
                        g->SetMarkerColor((iv < 4) ? variantColors[iv] : kBlack);
                        g->SetMarkerStyle((iv < 4) ? variantMarkers[iv] : 20);
                        g->SetMarkerSize(1.2);
                        sGraphs.push_back(g);
                        sIdx.push_back(iv);
                    }
                    if (!sGraphs.empty())
                    {
                        if (!std::isfinite(yMinS) || !std::isfinite(yMaxS)) { yMinS = 0; yMaxS = 5; }
                        const double padS = (yMaxS > yMinS) ? 0.20 * (yMaxS - yMinS) : 0.5;
                        TCanvas cGS(TString::Format("c_gaussSigmaVsPt_%s", cb.folder.c_str()).Data(), "", 900, 700);
                        ApplyCanvasMargins1D(cGS); cGS.cd();
                        TH1F hFr(TString::Format("hFr_gaussSigmaVsPt_%s", cb.folder.c_str()).Data(),
                                 "", 100, kPtEdges.front(), kPtEdges.back());
                        hFr.SetDirectory(nullptr); hFr.SetStats(0);
                        hFr.SetMinimum(std::max(0.0, yMinS - padS)); hFr.SetMaximum(yMaxS + padS);
                        hFr.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                        hFr.GetYaxis()->SetTitle("Gaussian #sigma [GeV]");
                        hFr.GetXaxis()->SetTitleSize(0.055); hFr.GetYaxis()->SetTitleSize(0.055);
                        hFr.GetXaxis()->SetLabelSize(0.045); hFr.GetYaxis()->SetLabelSize(0.045);
                        hFr.GetYaxis()->SetTitleOffset(1.15);
                        hFr.Draw();
                        for (auto* g : sGraphs) g->Draw("PE1 SAME");
                        TLegend lgGS(0.56, 0.62, 0.92, 0.88);
                        lgGS.SetBorderSize(0); lgGS.SetFillStyle(0); lgGS.SetTextFont(42); lgGS.SetTextSize(0.032);
                        for (std::size_t ig = 0; ig < sGraphs.size(); ++ig)
                            lgGS.AddEntry(sGraphs[ig], handles[sIdx[ig]].label.c_str(), "ep");
                        lgGS.Draw();
                        TLatex tGS; tGS.SetNDC(true); tGS.SetTextFont(42); tGS.SetTextAlign(23); tGS.SetTextSize(0.042);
                        tGS.DrawLatex(0.50, 0.98,
                            TString::Format("Gaussian #sigma vs p_{T}^{#gamma}, %d-%d%% Cent", cb.lo, cb.hi).Data());
                        TLatex tGSi; tGSi.SetNDC(true); tGSi.SetTextFont(42); tGSi.SetTextAlign(13); tGSi.SetTextSize(0.028);
                        tGSi.DrawLatex(0.20, 0.58, TString::Format("Trigger = %s", trigDisplayLabel.c_str()).Data());
                        tGSi.DrawLatex(0.20, 0.54, TString::Format("#DeltaR_{cone} < %.1f", (kAA_IsoConeR == "isoR40") ? 0.4 : 0.3).Data());
                        SaveCanvas(cGS, JoinPath(gaussSigmaPtDir,
                            TString::Format("gaussSigma_allVariants_vs_pT_%s.png", cb.folder.c_str()).Data()));
                        for (auto* g : sGraphs) delete g;
                    }
                }
            }
        }
        
        // ====== ueCompBase: Gaussian mean & sigma vs pT, centrality overlay, per variant ======
        {
            const int centOvColors[] = {kRed+1, kBlue+1, kGreen+2, kMagenta+1, kOrange+1,
                kCyan+1, kYellow+2, kViolet+1};
            
            const string gaussCentOvBase = JoinPath(ueCompBase, "gaussianSummary_centOverlay");
            const string gaussMeanCentOvDir  = JoinPath(gaussCentOvBase, "meanVsPt");
            const string gaussSigmaCentOvDir = JoinPath(gaussCentOvBase, "sigmaVsPt");
            EnsureDir(gaussCentOvBase);
            EnsureDir(gaussMeanCentOvDir);
            EnsureDir(gaussSigmaCentOvDir);
            
            for (std::size_t iv = 0; iv < handles.size(); ++iv)
            {
                if (!handles[iv].file) continue;
                const auto& H = handles[iv];
                
                // --- Gaussian Mean vs pT, centrality overlay ---
                {
                    vector<TGraphErrors*> gCent;
                    vector<string> cLbls;
                    double yMinG = 1e30, yMaxG = -1e30;
                    
                    for (std::size_t ic = 0; ic < centBins.size(); ++ic)
                    {
                        const auto& cb = centBins[ic];
                        vector<double> xP, exP, yP, eyP;
                        for (int ipt = 0; ipt < kNPtBins; ++ipt)
                        {
                            if (!gaussFilled[iv][ipt][ic]) continue;
                            xP.push_back(0.5 * (kPtEdges[(std::size_t)ipt] + kPtEdges[(std::size_t)ipt + 1]));
                            exP.push_back(0.5 * (kPtEdges[(std::size_t)ipt + 1] - kPtEdges[(std::size_t)ipt]));
                            yP.push_back(gaussMean[iv][ipt][ic]);
                            eyP.push_back(gaussMeanErr[iv][ipt][ic]);
                            yMinG = std::min(yMinG, yP.back() - eyP.back());
                            yMaxG = std::max(yMaxG, yP.back() + eyP.back());
                        }
                        if (xP.empty()) continue;
                        TGraphErrors* g = new TGraphErrors((int)xP.size(), &xP[0], &yP[0], &exP[0], &eyP[0]);
                        const int ci = (ic < 8) ? centOvColors[ic] : kBlack;
                        g->SetLineWidth(2);
                        g->SetLineColor(ci);
                        g->SetMarkerColor(ci);
                        g->SetMarkerStyle(20);
                        g->SetMarkerSize(1.2);
                        gCent.push_back(g);
                        cLbls.push_back(TString::Format("%d-%d%%", cb.lo, cb.hi).Data());
                    }
                    
                    if (!gCent.empty())
                    {
                        if (!std::isfinite(yMinG) || !std::isfinite(yMaxG)) { yMinG = 0; yMaxG = 10; }
                        const double padG = (yMaxG > yMinG) ? 0.20 * (yMaxG - yMinG) : 0.5;
                        TCanvas cGM(TString::Format("c_gaussMeanVsPt_centOv_%s", H.variant.c_str()).Data(), "", 900, 700);
                        ApplyCanvasMargins1D(cGM); cGM.cd();
                        TH1F hFr(TString::Format("hFr_gaussMeanVsPt_centOv_%s", H.variant.c_str()).Data(),
                                 "", 100, kPtEdges.front(), kPtEdges.back());
                        hFr.SetDirectory(nullptr); hFr.SetStats(0);
                        hFr.SetMinimum(yMinG - padG); hFr.SetMaximum(yMaxG + padG);
                        hFr.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                        hFr.GetYaxis()->SetTitle("Gaussian #mu [GeV]");
                        hFr.GetXaxis()->SetTitleSize(0.055); hFr.GetYaxis()->SetTitleSize(0.055);
                        hFr.GetXaxis()->SetLabelSize(0.045); hFr.GetYaxis()->SetLabelSize(0.045);
                        hFr.GetYaxis()->SetTitleOffset(1.15);
                        hFr.Draw();
                        for (auto* g : gCent) g->Draw("PE1 SAME");
                        TLegend lgGM(0.60, 0.62, 0.92, 0.88);
                        lgGM.SetBorderSize(0); lgGM.SetFillStyle(0); lgGM.SetTextFont(42); lgGM.SetTextSize(0.030);
                        for (std::size_t ig = 0; ig < gCent.size(); ++ig)
                            lgGM.AddEntry(gCent[ig], cLbls[ig].c_str(), "ep");
                        lgGM.Draw();
                        TLatex tGM; tGM.SetNDC(true); tGM.SetTextFont(42); tGM.SetTextAlign(23); tGM.SetTextSize(0.042);
                        tGM.DrawLatex(0.50, 0.98,
                            TString::Format("Gaussian #mu vs p_{T}^{#gamma}, centrality overlay (%s)", H.label.c_str()).Data());
                        TLatex tGMi; tGMi.SetNDC(true); tGMi.SetTextFont(42); tGMi.SetTextAlign(13); tGMi.SetTextSize(0.028);
                        tGMi.DrawLatex(0.18, 0.89, TString::Format("Trigger = %s", trigDisplayLabel.c_str()).Data());
                        tGMi.DrawLatex(0.18, 0.85, TString::Format("#DeltaR_{cone} < %.1f", (kAA_IsoConeR == "isoR40") ? 0.4 : 0.3).Data());
                        tGMi.DrawLatex(0.18, 0.81, TString::Format("UE: %s", H.label.c_str()).Data());
                        {
                            TLatex tSph; tSph.SetNDC(true); tSph.SetTextFont(42); tSph.SetTextAlign(33);
                            tSph.SetTextSize(0.048);
                            tSph.DrawLatex(0.92, 0.55, "#bf{sPHENIX} #it{Internal}");
                            tSph.SetTextSize(0.038);
                            tSph.DrawLatex(0.92, 0.49, "Au+Au  #sqrt{s_{NN}} = 200 GeV");
                        }
                        SaveCanvas(cGM, JoinPath(gaussMeanCentOvDir,
                            TString::Format("gaussMean_vs_pT_centOverlay_%s.png", H.variant.c_str()).Data()));
                        for (auto* g : gCent) delete g;
                    }
                }
                
                // --- Gaussian Sigma vs pT, centrality overlay ---
                {
                    vector<TGraphErrors*> sCent;
                    vector<string> sLbls;
                    double yMinS = 1e30, yMaxS = -1e30;
                    
                    for (std::size_t ic = 0; ic < centBins.size(); ++ic)
                    {
                        const auto& cb = centBins[ic];
                        vector<double> xP, exP, yP, eyP;
                        for (int ipt = 0; ipt < kNPtBins; ++ipt)
                        {
                            if (!gaussFilled[iv][ipt][ic]) continue;
                            xP.push_back(0.5 * (kPtEdges[(std::size_t)ipt] + kPtEdges[(std::size_t)ipt + 1]));
                            exP.push_back(0.5 * (kPtEdges[(std::size_t)ipt + 1] - kPtEdges[(std::size_t)ipt]));
                            yP.push_back(gaussSigma[iv][ipt][ic]);
                            eyP.push_back(gaussSigmaErr[iv][ipt][ic]);
                            yMinS = std::min(yMinS, yP.back() - eyP.back());
                            yMaxS = std::max(yMaxS, yP.back() + eyP.back());
                        }
                        if (xP.empty()) continue;
                        TGraphErrors* g = new TGraphErrors((int)xP.size(), &xP[0], &yP[0], &exP[0], &eyP[0]);
                        const int ci = (ic < 8) ? centOvColors[ic] : kBlack;
                        g->SetLineWidth(2);
                        g->SetLineColor(ci);
                        g->SetMarkerColor(ci);
                        g->SetMarkerStyle(20);
                        g->SetMarkerSize(1.2);
                        sCent.push_back(g);
                        sLbls.push_back(TString::Format("%d-%d%%", cb.lo, cb.hi).Data());
                    }
                    
                    if (!sCent.empty())
                    {
                        if (!std::isfinite(yMinS) || !std::isfinite(yMaxS)) { yMinS = 0; yMaxS = 5; }
                        const double padS = (yMaxS > yMinS) ? 0.20 * (yMaxS - yMinS) : 0.5;
                        TCanvas cGS(TString::Format("c_gaussSigmaVsPt_centOv_%s", H.variant.c_str()).Data(), "", 900, 700);
                        ApplyCanvasMargins1D(cGS); cGS.cd();
                        TH1F hFr(TString::Format("hFr_gaussSigmaVsPt_centOv_%s", H.variant.c_str()).Data(),
                                 "", 100, kPtEdges.front(), kPtEdges.back());
                        hFr.SetDirectory(nullptr); hFr.SetStats(0);
                        hFr.SetMinimum(std::max(0.0, yMinS - padS)); hFr.SetMaximum(yMaxS + padS);
                        hFr.GetXaxis()->SetTitle("p_{T}^{#gamma} [GeV]");
                        hFr.GetYaxis()->SetTitle("Gaussian #sigma [GeV]");
                        hFr.GetXaxis()->SetTitleSize(0.055); hFr.GetYaxis()->SetTitleSize(0.055);
                        hFr.GetXaxis()->SetLabelSize(0.045); hFr.GetYaxis()->SetLabelSize(0.045);
                        hFr.GetYaxis()->SetTitleOffset(1.15);
                        hFr.Draw();
                        for (auto* g : sCent) g->Draw("PE1 SAME");
                        TLegend lgGS(0.60, 0.62, 0.92, 0.88);
                        lgGS.SetBorderSize(0); lgGS.SetFillStyle(0); lgGS.SetTextFont(42); lgGS.SetTextSize(0.030);
                        for (std::size_t ig = 0; ig < sCent.size(); ++ig)
                            lgGS.AddEntry(sCent[ig], sLbls[ig].c_str(), "ep");
                        lgGS.Draw();
                        TLatex tGS; tGS.SetNDC(true); tGS.SetTextFont(42); tGS.SetTextAlign(23); tGS.SetTextSize(0.042);
                        tGS.DrawLatex(0.50, 0.98,
                            TString::Format("Gaussian #sigma vs p_{T}^{#gamma}, centrality overlay (%s)", H.label.c_str()).Data());
                        TLatex tGSi; tGSi.SetNDC(true); tGSi.SetTextFont(42); tGSi.SetTextAlign(13); tGSi.SetTextSize(0.028);
                        tGSi.DrawLatex(0.18, 0.89, TString::Format("Trigger = %s", trigDisplayLabel.c_str()).Data());
                        tGSi.DrawLatex(0.18, 0.85, TString::Format("#DeltaR_{cone} < %.1f", (kAA_IsoConeR == "isoR40") ? 0.4 : 0.3).Data());
                        tGSi.DrawLatex(0.18, 0.81, TString::Format("UE: %s", H.label.c_str()).Data());
                        {
                            TLatex tSph; tSph.SetNDC(true); tSph.SetTextFont(42); tSph.SetTextAlign(33);
                            tSph.SetTextSize(0.048);
                            tSph.DrawLatex(0.92, 0.55, "#bf{sPHENIX} #it{Internal}");
                            tSph.SetTextSize(0.038);
                            tSph.DrawLatex(0.92, 0.49, "Au+Au  #sqrt{s_{NN}} = 200 GeV");
                        }
                        SaveCanvas(cGS, JoinPath(gaussSigmaCentOvDir,
                            TString::Format("gaussSigma_vs_pT_centOverlay_%s.png", H.variant.c_str()).Data()));
                        for (auto* g : sCent) delete g;
                    }
                }
            }
        }
        
        if (!forEmbeddedSim && !generateISOpTcentOverlaysONLY)
        {
            const string tightNonTightBase = JoinPath(ueCompBase, "tightNonTightOverlays");
            EnsureDir(tightNonTightBase);
            
            auto MakeTightNonTightOverlays =
            [&](const string& mcOverlayFolder,
                const string& mcTitleTag,
                const string& tightMcLegend,
                const string& nonTightMcLegend,
                auto&& ResolveMCInput)
            {
                if (mcTitleTag.empty()) return;
                
                vector<VariantHandle> mcHandles;
                mcHandles.reserve(ueVariants.size());
                
                for (std::size_t iv = 0; iv < ueVariants.size(); ++iv)
                {
                    VariantHandle M;
                    M.variant = ueVariants[iv];
                    M.label   = ueLabels[iv];
                    
                    const string mcInput = ResolveMCInput(M.variant);
                    if (!mcInput.empty())
                    {
                        M.file = TFile::Open(mcInput.c_str(), "READ");
                        if (!M.file || M.file->IsZombie())
                        {
                            if (M.file) { M.file->Close(); delete M.file; M.file = nullptr; }
                        }
                    }
                    
                    mcHandles.push_back(std::move(M));
                }
                
                for (std::size_t iv = 0; iv < handles.size(); ++iv)
                {
                    if (iv >= mcHandles.size()) continue;
                    
                    auto& dataH = handles[iv];
                    auto& mcH   = mcHandles[iv];
                    
                    if (!dataH.file || !mcH.file) continue;
                    
                    TDirectory* dataTop = dataH.file->GetDirectory(trigAA.c_str());
                    TDirectory* mcTop   = mcH.file->GetDirectory(kDirSIM.c_str());
                    if (!dataTop || !mcTop) continue;
                    
                    const string variantBaseDir = JoinPath(JoinPath(tightNonTightBase, dataH.variant), mcOverlayFolder);
                    EnsureDir(JoinPath(tightNonTightBase, dataH.variant));
                    EnsureDir(variantBaseDir);
                    
                    for (const auto& cb : centBins)
                    {
                        const string centDir = JoinPath(variantBaseDir, cb.folder);
                        EnsureDir(centDir);
                        
                        for (const auto& b : PtBins())
                        {
                            const string ptDir = JoinPath(centDir, b.folder);
                            EnsureDir(ptDir);
                            
                            const string hTightName    = "h_Eiso_tight" + b.suffix + cb.suffix;
                            const string hNonTightName = "h_Eiso_nonTight" + b.suffix + cb.suffix;
                            
                            TH1* hTDataSrc  = dynamic_cast<TH1*>(dataTop->Get(hTightName.c_str()));
                            TH1* hNTDataSrc = dynamic_cast<TH1*>(dataTop->Get(hNonTightName.c_str()));
                            TH1* hTMcSrc    = dynamic_cast<TH1*>(mcTop->Get(hTightName.c_str()));
                            TH1* hNTMcSrc   = dynamic_cast<TH1*>(mcTop->Get(hNonTightName.c_str()));
                            if (!hTDataSrc || !hNTDataSrc || !hTMcSrc || !hNTMcSrc) continue;
                            
                            TH1* hTData = CloneTH1(
                                                   hTDataSrc,
                                                   TString::Format("hTightData_%s_%s_%s_%s",
                                                                   mcOverlayFolder.c_str(), dataH.variant.c_str(), cb.folder.c_str(), b.folder.c_str()).Data()
                                                   );
                            TH1* hNTData = CloneTH1(
                                                    hNTDataSrc,
                                                    TString::Format("hNonTightData_%s_%s_%s_%s",
                                                                    mcOverlayFolder.c_str(), dataH.variant.c_str(), cb.folder.c_str(), b.folder.c_str()).Data()
                                                    );
                            TH1* hTMc = CloneTH1(
                                                 hTMcSrc,
                                                 TString::Format("hTightMC_%s_%s_%s_%s",
                                                                 mcOverlayFolder.c_str(), dataH.variant.c_str(), cb.folder.c_str(), b.folder.c_str()).Data()
                                                 );
                            TH1* hNTMc = CloneTH1(
                                                  hNTMcSrc,
                                                  TString::Format("hNonTightMC_%s_%s_%s_%s",
                                                                  mcOverlayFolder.c_str(), dataH.variant.c_str(), cb.folder.c_str(), b.folder.c_str()).Data()
                                                  );
                            
                            if (!hTData || !hNTData || !hTMc || !hNTMc)
                            {
                                if (hTData) delete hTData;
                                if (hNTData) delete hNTData;
                                if (hTMc) delete hTMc;
                                if (hNTMc) delete hNTMc;
                                continue;
                            }
                            
                            EnsureSumw2(hTData);
                            EnsureSumw2(hNTData);
                            EnsureSumw2(hTMc);
                            EnsureSumw2(hNTMc);
                            
                            hTData->Rebin(10);
                            hNTData->Rebin(10);
                            hTMc->Rebin(10);
                            hNTMc->Rebin(10);
                            
                            const double intTData  = hTData->Integral(0, hTData->GetNbinsX() + 1);
                            const double intNTData = hNTData->Integral(0, hNTData->GetNbinsX() + 1);
                            const double intTMc    = hTMc->Integral(0, hTMc->GetNbinsX() + 1);
                            const double intNTMc   = hNTMc->Integral(0, hNTMc->GetNbinsX() + 1);
                            
                            if (!(intTData > 0.0) || !(intNTData > 0.0) || !(intTMc > 0.0) || !(intNTMc > 0.0))
                            {
                                delete hTData;
                                delete hNTData;
                                delete hTMc;
                                delete hNTMc;
                                continue;
                            }
                            
                            hTData->Scale(1.0 / intTData);
                            hNTData->Scale(1.0 / intNTData);
                            hTMc->Scale(1.0 / intTMc);
                            hNTMc->Scale(1.0 / intNTMc);
                            
                            hTMc->SetTitle("");
                            hTMc->SetLineColor(kBlack);
                            hTMc->SetLineWidth(2);
                            hTMc->SetFillStyle(0);
                            hTMc->SetMarkerSize(0.0);
                            
                            hTData->SetLineColor(kBlack);
                            hTData->SetMarkerColor(kBlack);
                            hTData->SetMarkerStyle(20);
                            hTData->SetMarkerSize(1.0);
                            hTData->SetLineWidth(2);
                            hTData->SetFillStyle(0);
                            
                            hNTMc->SetLineColor(kRed + 1);
                            hNTMc->SetLineWidth(2);
                            hNTMc->SetFillStyle(0);
                            hNTMc->SetMarkerSize(0.0);
                            
                            hNTData->SetLineColor(kRed + 1);
                            hNTData->SetMarkerColor(kRed + 1);
                            hNTData->SetMarkerStyle(24);
                            hNTData->SetMarkerSize(1.0);
                            hNTData->SetLineWidth(2);
                            hNTData->SetFillStyle(0);
                            
                            const double yMaxTNT = std::max(
                                                            std::max(hTData->GetMaximum(), hNTData->GetMaximum()),
                                                            std::max(hTMc->GetMaximum(), hNTMc->GetMaximum())
                                                            );
                            
                            TCanvas cTNT(
                                         TString::Format("c_tightNonTight_%s_%s_%s_%s",
                                                         mcOverlayFolder.c_str(), dataH.variant.c_str(), cb.folder.c_str(), b.folder.c_str()).Data(),
                                         "c_tightNonTight", 900, 700
                                         );
                            ApplyCanvasMargins1D(cTNT);
                            cTNT.cd();
                            
                            hTMc->GetXaxis()->SetTitle("E_{T}^{iso} [GeV]");
                            hTMc->GetYaxis()->SetTitle("Normalized to unit area");
                            hTMc->GetXaxis()->SetTitleSize(0.055);
                            hTMc->GetYaxis()->SetTitleSize(0.055);
                            hTMc->GetXaxis()->SetLabelSize(0.045);
                            hTMc->GetYaxis()->SetLabelSize(0.045);
                            hTMc->GetYaxis()->SetTitleOffset(1.15);
                            hTMc->SetMinimum(0.0);
                            hTMc->SetMaximum((yMaxTNT > 0.0) ? (1.25 * yMaxTNT) : 1.0);
                            
                            for (int ib = 0; ib <= hTMc->GetNbinsX() + 1; ++ib) hTMc->SetBinError(ib, 0.0);
                            for (int ib = 0; ib <= hNTMc->GetNbinsX() + 1; ++ib) hNTMc->SetBinError(ib, 0.0);
                            hTMc->Draw("hist");
                            hNTMc->Draw("hist SAME");
                            hTData->Draw("E1 SAME");
                            hNTData->Draw("E1 SAME");
                            
                            TLegend legTNT(0.58, 0.68, 0.92, 0.88);
                            legTNT.SetBorderSize(0);
                            legTNT.SetFillStyle(0);
                            legTNT.SetTextFont(42);
                            legTNT.SetTextSize(0.032);
                            legTNT.AddEntry(hTData, "tight data", "ep");
                            legTNT.AddEntry(hTMc, tightMcLegend.c_str(), "l");
                            legTNT.AddEntry(hNTData, "nontight data", "ep");
                            legTNT.AddEntry(hNTMc, nonTightMcLegend.c_str(), "l");
                            legTNT.Draw();
                            
                            TLatex tTitleTNT;
                            tTitleTNT.SetNDC(true);
                            tTitleTNT.SetTextFont(42);
                            tTitleTNT.SetTextAlign(23);
                            tTitleTNT.SetTextSize(0.040);
                            tTitleTNT.DrawLatex(0.50, 0.97,
                                                TString::Format("E_{T}^{iso} tight/nontight, auau data and %s embedded MC",
                                                                mcTitleTag.c_str()).Data());
                            
                            TLatex tInfoTNT;
                            tInfoTNT.SetNDC(true);
                            tInfoTNT.SetTextFont(42);
                            tInfoTNT.SetTextAlign(13);
                            tInfoTNT.SetTextSize(0.050);
                            tInfoTNT.DrawLatex(0.22, 0.88, TString::Format("%d-%d%%", cb.lo, cb.hi).Data());
                            tInfoTNT.DrawLatex(0.22, 0.82, TString::Format("p_{T}^{#gamma} = %d-%d GeV", b.lo, b.hi).Data());
                            
                            // ---- sPHENIX Internal + Au+Au between legend and trigger info ----
                            TLatex tSph;
                            tSph.SetNDC(true);
                            tSph.SetTextFont(42);
                            tSph.SetTextAlign(33);
                            tSph.SetTextSize(0.042);
                            tSph.DrawLatex(0.92, 0.64, "#bf{sPHENIX} #it{Internal}");
                            tSph.SetTextSize(0.034);
                            tSph.DrawLatex(0.92, 0.59, "Au+Au  #sqrt{s_{NN}} = 200 GeV");
                            
                            TLatex tUE;
                            tUE.SetNDC(true);
                            tUE.SetTextFont(42);
                            tUE.SetTextAlign(33);
                            tUE.SetTextSize(0.030);
                            tUE.DrawLatex(0.92, 0.53, trigDisplayLabel.c_str());
                            tUE.DrawLatex(0.92, 0.49, TString::Format("|v_{z}| < %d cm", kAA_VzCut).Data());
                            tUE.DrawLatex(0.92, 0.45, TString::Format("UE: %s", dataH.label.c_str()).Data());
                            tUE.DrawLatex(0.92, 0.41, TString::Format("#DeltaR_{cone} < %.1f", (kAA_IsoConeR == "isoR40") ? 0.4 : 0.3).Data());
                            
                            SaveCanvas(cTNT, JoinPath(ptDir, "Eiso_tightNonTight_overlay.png"));
                            
                            // --- save 4 individual unnormalized distributions ---
                            {
                                const std::vector<std::pair<TH1*,std::string>> rawVec = {
                                    {hTDataSrc,  "Eiso_tight_data"},
                                    {hNTDataSrc, "Eiso_nonTight_data"},
                                    {hTMcSrc,    "Eiso_tight_mc"},
                                    {hNTMcSrc,   "Eiso_nonTight_mc"}
                                };
                                for (const auto& rv : rawVec)
                                {
                                    TH1* hRaw = dynamic_cast<TH1*>(rv.first->Clone(
                                                                                   TString::Format("hRaw_%s_%s_%s", rv.second.c_str(), cb.folder.c_str(), b.folder.c_str()).Data()));
                                    if (!hRaw) continue;
                                    hRaw->SetDirectory(nullptr);
                                    hRaw->SetTitle("");
                                    hRaw->GetXaxis()->SetTitle("E_{T}^{iso} [GeV]");
                                    hRaw->GetYaxis()->SetTitle("Counts");
                                    hRaw->GetXaxis()->SetTitleSize(0.055);
                                    hRaw->GetYaxis()->SetTitleSize(0.055);
                                    hRaw->GetXaxis()->SetLabelSize(0.045);
                                    hRaw->GetYaxis()->SetLabelSize(0.045);
                                    hRaw->GetYaxis()->SetTitleOffset(1.15);
                                    const bool isMC = (rv.second.find("_mc") != std::string::npos);
                                    const bool isNT = (rv.second.find("nonTight") != std::string::npos);
                                    const Color_t col = isNT ? (kRed + 1) : kBlack;
                                    hRaw->SetLineColor(col);
                                    hRaw->SetLineWidth(2);
                                    hRaw->SetFillStyle(0);
                                    TCanvas cRaw(
                                                 TString::Format("c_raw_%s_%s_%s", rv.second.c_str(), cb.folder.c_str(), b.folder.c_str()).Data(),
                                                 "c_raw", 900, 700);
                                    ApplyCanvasMargins1D(cRaw);
                                    cRaw.cd();
                                    if (isMC)
                                    {
                                        hRaw->SetMarkerSize(0.0);
                                        hRaw->Draw("hist");
                                    }
                                    else
                                    {
                                        hRaw->SetMarkerColor(col);
                                        hRaw->SetMarkerStyle(isNT ? 24 : 20);
                                        hRaw->SetMarkerSize(1.0);
                                        hRaw->Draw("E1");
                                    }
                                    TLatex tRawInfo;
                                    tRawInfo.SetNDC(true);
                                    tRawInfo.SetTextFont(42);
                                    tRawInfo.SetTextAlign(13);
                                    tRawInfo.SetTextSize(0.045);
                                    tRawInfo.DrawLatex(0.16, 0.88, TString::Format("%d-%d%%", cb.lo, cb.hi).Data());
                                    tRawInfo.DrawLatex(0.16, 0.82, TString::Format("p_{T}^{#gamma} = %d-%d GeV", b.lo, b.hi).Data());
                                    SaveCanvas(cRaw, JoinPath(ptDir, rv.second + ".png"));
                                    delete hRaw;
                                }
                            }
                            
                            delete hTData;
                            delete hNTData;
                            delete hTMc;
                            delete hNTMc;
                        }
                    }
                }
                
                for (auto& M : mcHandles)
                {
                    if (M.file)
                    {
                        M.file->Close();
                        delete M.file;
                        M.file = nullptr;
                    }
                }
            };
            
            MakeTightNonTightOverlays(
                                      "inclusiveMCoverlays",
                                      inclusiveEmbeddedTitleTag,
                                      "tight inclusive MC",
                                      "nontight inclusive MC",
                                      ResolveInclusiveJetVariantInput
                                      );
            
            MakeTightNonTightOverlays(
                                      "photonJetOverlays",
                                      photonEmbeddedTitleTag,
                                      "tight photon+jet MC",
                                      "nontight photon+jet MC",
                                      ResolvePhotonJetVariantInput
                                      );
        }
        if (generateUEcomparisonSSQA && !skipToCentralityAndPtOverlaysWithSSQA && !SSoverlayPerVAR_processONLY) {
#include "AnalyzeRecoilJets_SSQA.cpp"
        }
        for (auto& H : handles)
        {
            if (H.file)
            {
                H.file->Close();
                delete H.file;
                H.file = nullptr;
            }
        }
    }
    
    if (fPPSim510_20)
    {
        fPPSim510_20->Close();
        delete fPPSim510_20;
        fPPSim510_20 = nullptr;
    }
    if (fPP)
    {
        fPP->Close();
        delete fPP;
        fPP = nullptr;
    }
}
